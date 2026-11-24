"""Card C-63 (DEC-71, DEC-72): stage ART-RT-enriched GlobDB r226 bacterial phyla as GART1 and GART2.

    python c63_globdb_art_phyla.py stage NON_GTDB_STAGES_TSV GENOME_FASTA_TAR_GZ OUT_ROOT
        [--fasta-dir DIR] [--final-dir STAGE=DIR ...] [--root ROOT] [--root-stages S1,P1,...]
        [--owner-inputs NAME=DIR ...] [--retrieved-at 2026-10-23T19:51:00Z]
        [--workers 2] [--gart1-screen-max-bp 500000] [--bgzip]
        [--work WORK_DIR] [--threads 2] [--memory-limit 4GB] [--report-out JSON]

Stage split (stage_driver.sh requires ^[A-Z][A-Za-z0-9]{0,15}$, no underscores):
- GART1: Myxococcota + Myxococcota_A + Cyanobacteriota + Omnitrophota + Planctomycetota
  (9,234 non-GTDB species representatives, 39.413 Gbp, 5,218,293 contigs).
- GART2: Patescibacteriota (15,320 non-GTDB CPR species representatives, 11.853 Gbp, 1,213,095 contigs).

Streams the local `globdb_r226_genome_fasta.tar.gz` archive ONCE in bounded memory, verifies each
selected genome's record and base counts against `non_gtdb_stages.tsv`, screens genetic codes
(pyrodigal `translation_table=11` vs `25` vs `4`, assigning code 25/4 only when coding density
gain >= 0.10 and mean ORF length ratio >= 1.25), consolidates each stage into one deterministic
`.fna.gz` file (>> 64 MB, AGENTS.md §7.3), deduplicates `contig_uid` ownership against earlier
stages (DEC-26 / DEC-27), and writes `genomes.json`, `objects.json`, `records.tsv`,
`repeated_records.tsv`, `shared_records.tsv`, and `genetic_codes.tsv`.
"""

from __future__ import annotations

import argparse
import array
import csv
import gzip
import io
import json
import math
import shutil
import subprocess
import sys
import tarfile
import time
from collections import Counter, defaultdict
from collections.abc import Iterable
from concurrent.futures import Future, ProcessPoolExecutor
from pathlib import Path

import duckdb
import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parent))

import make_inputs  # noqa: E402
from biobox import fasta, ids, stage as bb_stage  # noqa: E402
from biobox.l0seq import store  # noqa: E402
from biobox.l1 import genecall  # noqa: E402

Refused = make_inputs.Refused

GART1_PHYLA = ("Myxococcota", "Myxococcota_A", "Cyanobacteriota", "Omnitrophota", "Planctomycetota")
GART2_PHYLA = ("Patescibacteriota",)
STAGE_PHYLA: dict[str, tuple[str, ...]] = {
    "GART1": GART1_PHYLA,
    "GART2": GART2_PHYLA,
}
STAGE_FASTA: dict[str, str] = {
    "GART1": "gart1_globdb.fna.gz",
    "GART2": "gart2_globdb_cpr.fna.gz",
}

GLOBDB_LICENSE = (
    "CC BY-SA 4.0 (GlobDB r226, Speth et al. 2025, doi:10.1093/bioadv/vbaf280; "
    "propagates licenses of underlying 13 MAG catalogues)"
)
GLOBDB_SOURCE_BASE = {
    "source_name": "globdb",
    "source_release": "r226",
    "license": GLOBDB_LICENSE,
    "redistributable": True,
}

# Known CPR UGA->Gly (code 25) taxonomic markers in GTDB R226 taxonomy
_CODE25_TAX_MARKERS = ("c__Gracilibacteria", "c__Absconditabacteria", "o__Absconditabacterales", "o__BD1-5")

RECORDS_HEAD = [
    "ordinal", "genome", "record_ordinal", "contig_uid", "length",
    "phylum", "resource", "assembly_quality", "genetic_code", "c41_stage",
]
GC_HEAD = [
    "genome", "stage", "phylum", "tax_class", "tax_order", "resource",
    "genome_size", "contig_number", "assembly_quality",
    "genetic_code", "reason", "tested_alt",
    "cd_11", "mean_orf_11", "n_genes_11",
    "cd_25", "mean_orf_25", "n_genes_25",
    "cd_4", "mean_orf_4", "n_genes_4",
    "tga_per_kb",
]


def _sql_lit(path: Path | str) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def _sql_list(files: Iterable[Path]) -> str:
    return "[" + ", ".join(_sql_lit(f) for f in files) + "]"


def mag_quality(comp: float, cont: float) -> str:
    """MIMAG completeness/contamination tier mapped to `genome_input.assembly_quality`."""
    if comp > 90.0 and cont < 5.0:
        return "mag_hq"
    if comp >= 50.0 and cont < 10.0:
        return "mag_mq"
    return "contig"


def parse_tax_ranks(taxonomy: str) -> dict[str, str]:
    """Extract GTDB/GlobDB taxonomic ranks (d, p, c, o, f, g, s) from a semicolon-delimited string."""
    ranks: dict[str, str] = {}
    for tok in taxonomy.split(";"):
        tok = tok.strip()
        if "__" in tok:
            k, v = tok.split("__", 1)
            ranks[k] = v
        elif "_" in tok:
            # Handle GlobDB's `p_Unclassified` quirk (single underscore)
            k, v = tok.split("_", 1)
            ranks[k] = v
    return ranks


def parse_phylum(taxonomy: str) -> str:
    return parse_tax_ranks(taxonomy).get("p", "")


def select_stage_genomes(
    stages_tsv: Path,
    stages: Iterable[str] = ("GART1", "GART2"),
) -> dict[str, list[dict]]:
    """Select non-GTDB species representatives for each requested GART stage, sorted by (taxonomy, ID)."""
    stage_list = list(stages)
    phylum_to_stage: dict[str, str] = {}
    for st in stage_list:
        if st not in STAGE_PHYLA:
            raise Refused(f"unknown stage {st!r}; expected subset of {sorted(STAGE_PHYLA)}")
        for ph in STAGE_PHYLA[st]:
            phylum_to_stage[ph] = st

    selected: dict[str, list[dict]] = {st: [] for st in stage_list}
    seen_ids: set[str] = set()
    with open(stages_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        need = {
            "ID", "stage", "resource", "genome_size", "contig_number",
            "checkm2_completeness", "checkm2_contamination", "taxonomy",
        }
        missing = sorted(need - set(reader.fieldnames or ()))
        if missing:
            raise Refused(f"{stages_tsv}: missing column {missing[0]!r}")
        for r in reader:
            gid = r["ID"]
            if gid in seen_ids:
                raise Refused(f"{stages_tsv}: duplicate ID {gid!r}")
            seen_ids.add(gid)
            tax = r["taxonomy"]
            ranks = parse_tax_ranks(tax)
            phylum = ranks.get("p", "")
            st = phylum_to_stage.get(phylum)
            if st is None:
                continue
            comp = float(r["checkm2_completeness"])
            cont = float(r["checkm2_contamination"])
            selected[st].append({
                "ID": gid,
                "c41_stage": r["stage"],
                "resource": r["resource"],
                "genome_size": int(r["genome_size"]),
                "contig_number": int(r["contig_number"]),
                "checkm2_completeness": comp,
                "checkm2_contamination": cont,
                "assembly_quality": mag_quality(comp, cont),
                "taxonomy": tax,
                "phylum": phylum,
                "tax_class": ranks.get("c", ""),
                "tax_order": ranks.get("o", ""),
            })

    for st in stage_list:
        selected[st].sort(key=lambda r: (r["taxonomy"], r["ID"]))
    return selected


def _call_stats(seq_strs: list[str], total_bp: int, table: int, *, count_tga: bool = False) -> tuple[float, float, int, float]:
    """Train pyrodigal single-mode on `seq_strs` with `table`, return (coding_density, mean_orf_bp, n_genes, tga_per_kb)."""
    finder = pyrodigal.GeneFinder(meta=False)
    finder.train(*seq_strs, translation_table=table)
    coding_bp = 0
    n_genes = 0
    internal_tga = 0
    for s in seq_strs:
        for g in finder.find_genes(s):
            gl = g.end - g.begin + 1
            coding_bp += gl
            n_genes += 1
            if count_tga:
                nuc = g.sequence()
                # Count in-frame internal TGA codons (excluding terminal stop codon)
                for pos in range(0, len(nuc) - 3, 3):
                    if nuc[pos : pos + 3] == "TGA":
                        internal_tga += 1
    cd = coding_bp / total_bp if total_bp > 0 else 0.0
    mean_orf = coding_bp / n_genes if n_genes > 0 else 0.0
    tga_per_kb = (1000.0 * internal_tga / coding_bp) if (count_tga and coding_bp > 0) else 0.0
    return cd, mean_orf, n_genes, tga_per_kb


def _prefix_seqs_for_screen(seqs: list[bytes], max_bp: int) -> tuple[list[str], int, bool]:
    """Select longest contigs up to at least max_bp (and >= 100 kb) for fast code-11 pre-screen."""
    total_bp = sum(len(s) for s in seqs)
    if max_bp <= 0 or total_bp <= max_bp:
        return [s.upper().decode("ascii") for s in seqs], total_bp, True
    ordered = sorted(seqs, key=len, reverse=True)
    kept: list[str] = []
    acc = 0
    for s in ordered:
        kept.append(s.upper().decode("ascii"))
        acc += len(s)
        if acc >= max_bp and acc >= genecall.SINGLE_MODE_MIN_BP:
            break
    return kept, acc, False


def screen_genetic_code(
    seqs: list[bytes],
    taxonomy: str = "",
    *,
    cd_gain_min: float = 0.10,
    orf_ratio_min: float = 1.25,
    cd11_trigger: float = 0.82,
    orf11_trigger: float = 600.0,
    screen_max_bp: int = 0,
) -> dict:
    """Screen a prokaryotic genome (>= 100 kb) for UGA recoding (code 11 vs 25 vs 4) using pyrodigal."""
    pre_strs, pre_bp, is_full = _prefix_seqs_for_screen(seqs, screen_max_bp)
    cd_11, mean_orf_11, n_11, _ = _call_stats(pre_strs, pre_bp, 11)

    force_tax = any(m in taxonomy for m in _CODE25_TAX_MARKERS)
    triggered = force_tax or (cd_11 < cd11_trigger) or (mean_orf_11 < orf11_trigger)

    if not triggered:
        return {
            "genetic_code": 11,
            "reason": "code11_pass",
            "tested_alt": False,
            "cd_11": round(cd_11, 5),
            "mean_orf_11": round(mean_orf_11, 2),
            "n_genes_11": n_11,
            "cd_25": None,
            "mean_orf_25": None,
            "n_genes_25": None,
            "cd_4": None,
            "mean_orf_4": None,
            "n_genes_4": None,
            "tga_per_kb": None,
        }

    # Re-run code 11 on full genome if the pre-screen only used a prefix
    if not is_full:
        full_strs = [s.upper().decode("ascii") for s in seqs]
        full_bp = sum(len(s) for s in seqs)
        cd_11, mean_orf_11, n_11, _ = _call_stats(full_strs, full_bp, 11)
    else:
        full_strs, full_bp = pre_strs, pre_bp

    cd_25, mean_orf_25, n_25, _ = _call_stats(full_strs, full_bp, 25, count_tga=False)
    cd_gain_25 = cd_25 - cd_11
    orf_ratio_25 = (mean_orf_25 / mean_orf_11) if mean_orf_11 > 0 else float("inf")

    # Only run code 4 and TGA-codon density when code 25 actually improves over code 11
    cd_4: float | None = None
    mean_orf_4: float | None = None
    n_4: int | None = None
    tga_per_kb: float | None = None

    if cd_gain_25 >= cd_gain_min and orf_ratio_25 >= orf_ratio_min:
        _, _, _, tga_per_kb = _call_stats(full_strs, full_bp, 25, count_tga=True)
        cd_4, mean_orf_4, n_4, _ = _call_stats(full_strs, full_bp, 4)
        is_cpr = "p__Patescibacteriota" in taxonomy
        if is_cpr:
            code = 25
        elif cd_4 > cd_25 + 1e-4:
            code = 4
        else:
            code = 25
        reason = f"uga_recoded_code{code}"
    else:
        code = 11
        reason = "alt_tested_kept_11"

    return {
        "genetic_code": code,
        "reason": reason,
        "tested_alt": True,
        "cd_11": round(cd_11, 5),
        "mean_orf_11": round(mean_orf_11, 2),
        "n_genes_11": n_11,
        "cd_25": round(cd_25, 5),
        "mean_orf_25": round(mean_orf_25, 2),
        "n_genes_25": n_25,
        "cd_4": None if cd_4 is None else round(cd_4, 5),
        "mean_orf_4": None if mean_orf_4 is None else round(mean_orf_4, 2),
        "n_genes_4": n_4,
        "tga_per_kb": None if tga_per_kb is None else round(tga_per_kb, 3),
    }


def _screen_worker(args: tuple[str, list[bytes], str, int]) -> tuple[str, dict]:
    gid, seqs, taxonomy, max_bp = args
    return gid, screen_genetic_code(seqs, taxonomy, screen_max_bp=max_bp)


def _iter_fasta_bytes(raw_fa: bytes) -> list[tuple[str, bytes, int]]:
    """Parse FASTA bytes into (header, seq, line_width) tuples."""
    recs: list[tuple[str, bytes, int]] = []
    header: str | None = None
    seq_lines: list[bytes] = []
    for raw in raw_fa.splitlines():
        if not raw.strip():
            continue
        if raw.startswith(b">"):
            if header is not None:
                recs.append((header, b"".join(seq_lines), fasta._width(seq_lines)))
                seq_lines = []
            header = raw[1:].decode("utf-8")
        else:
            if header is None:
                raise Refused("FASTA sequence before header")
            seq_lines.append(raw)
    if header is not None:
        recs.append((header, b"".join(seq_lines), fasta._width(seq_lines)))
    return recs


def _open_tar_stream(tar_gz: Path):
    """Open local .tar.gz for streaming, preferring pigz/gzip subprocess for speed."""
    cmd = shutil.which("pigz") or shutil.which("gzip")
    if cmd and tar_gz.stat().st_size > 10_000_000:
        args = [cmd, "-dc", "-p", "2", str(tar_gz)] if Path(cmd).name == "pigz" else [cmd, "-dc", str(tar_gz)]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, bufsize=1 << 22)
        assert proc.stdout is not None
        tf = tarfile.open(fileobj=proc.stdout, mode="r|")

        def _close():
            tf.close()
            proc.stdout.close()
            rc = proc.wait()
            if rc != 0:
                raise Refused(f"{args[0]} -dc {tar_gz} exited {rc}")

        return tf, _close
    tf = tarfile.open(tar_gz, mode="r|gz")
    return tf, tf.close


class _GzipWriter:
    """Deterministic .fna.gz writer using bgzip -c (if requested) or gzip.GzipFile(mtime=0)."""

    def __init__(self, tmp_path: Path, use_bgzip: bool):
        self.tmp_path = tmp_path
        self._proc = None
        self._raw_fh = open(tmp_path, "wb", buffering=1 << 22)
        bgz = shutil.which("bgzip") if use_bgzip else None
        if bgz:
            self._proc = subprocess.Popen(
                [bgz, "-c"],
                stdin=subprocess.PIPE,
                stdout=self._raw_fh,
                bufsize=1 << 22,
            )
            assert self._proc.stdin is not None
            self._fh = self._proc.stdin
        else:
            self._fh = gzip.GzipFile(filename="", mode="wb", fileobj=self._raw_fh, mtime=0, compresslevel=6)

    def write(self, data: bytes) -> None:
        self._fh.write(data)

    def close(self) -> None:
        if self._proc is not None:
            self._fh.close()
            rc = self._proc.wait()
            self._raw_fh.close()
            if rc != 0:
                raise Refused(f"bgzip -c exited {rc}")
        else:
            self._fh.close()
            self._raw_fh.close()


def _owners_from_records(
    records_tsv: Path,
    root: Path | None,
    root_stages: list[str],
    owner_inputs: dict[str, Path],
    *,
    work: Path,
    threads: int = 2,
    memory_limit: str = "4GB",
) -> dict[bytes, str]:
    """Bounded-memory DuckDB lookup of earlier-stage contig_uid owners in build order (DEC-26)."""
    if root is None and not owner_inputs:
        return {}
    duck_dir = work / f"duck_own_{records_tsv.parent.name}"
    duck_dir.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(
        config={"threads": int(threads), "memory_limit": memory_limit, "temp_directory": str(duck_dir)}
    )
    try:
        con.execute(
            f"CREATE TEMP TABLE v AS SELECT DISTINCT from_hex(contig_uid) AS contig_uid "
            f"FROM read_csv({_sql_lit(records_tsv)}, delim='\t', header=true, quote='', all_varchar=true)"
        )
        out: dict[bytes, str] = {}
        stages = list(root_stages)
        if root is not None and not stages:
            stages = list(store.stages_on_disk(root))
        for s in stages:
            if root is None:
                continue
            part_dir = Path(root) / store.INDEX_DIR / f"stage={s}"
            if not part_dir.is_dir():
                continue
            files = sorted(part_dir.glob("part-*.parquet"))
            if not files:
                continue
            for (cu,) in con.execute(
                f"SELECT DISTINCT contig_uid FROM read_parquet({_sql_list(files)}, hive_partitioning = false) "
                "WHERE is_default_variant AND contig_uid IN (SELECT contig_uid FROM v) ORDER BY 1"
            ).fetchall():
                out.setdefault(bytes(cu), s)
        for name, d in owner_inputs.items():
            gj = Path(d) / "genomes.json"
            if gj.is_file():
                p = _sql_lit(gj)
                for (cu,) in con.execute(
                    f"SELECT DISTINCT from_hex(u) FROM (SELECT unnest(contig_uids) AS u FROM read_json({p}, "
                    "format='array', columns={'contig_uids': 'VARCHAR[]'})) "
                    "WHERE from_hex(u) IN (SELECT contig_uid FROM v) ORDER BY 1"
                ).fetchall():
                    out.setdefault(bytes(cu), name)
                continue
            rec = Path(d) / "records.tsv"
            if rec.is_file():
                for (cu,) in con.execute(
                    f"SELECT DISTINCT from_hex(contig_uid) FROM read_csv({_sql_lit(rec)}, delim='\t', "
                    "header=true, quote='', all_varchar=true) "
                    "WHERE from_hex(contig_uid) IN (SELECT contig_uid FROM v) ORDER BY 1"
                ).fetchall():
                    out.setdefault(bytes(cu), name)
        return out
    finally:
        con.close()


def _estimate_called_crispr(
    fasta_uncompressed_bytes: int,
    new_contigs: int,
    new_bases: int,
    held_bases: int,
) -> dict:
    """Compute --depth called --crispr CPU estimates using biobox.stage and biobox.l1.genecall constants."""
    ingest_min = (
        bb_stage.INGEST_FIXED_CPU_MIN
        + (max(0, fasta_uncompressed_bytes - held_bases) / bb_stage.INGEST_BP) * bb_stage.INGEST_CPU_MIN
        + (fasta_uncompressed_bytes / bb_stage.INGEST_BP) * bb_stage.INGEST_GZ_CPU_MIN
    )
    genecall_min = (
        (new_bases / 1_000_000.0) * genecall.CALL_S_PER_MBP[("single", False)] * genecall.STEP_OVERHEAD / 60.0
    )
    genecall_batches = max(1, math.ceil(genecall_min / 240.0)) if genecall_min > 240.0 else 1
    crispr_min = (
        (new_bases / 1_000_000_000.0) * bb_stage.CRISPR_CPU_S_PER_GBP / 60.0
        + (new_contigs / 1_000.0) * bb_stage.CRISPR_CPU_S_PER_K_CONTIGS / 60.0
    )
    return {
        "ingest_cpu_min": round(ingest_min, 2),
        "genecall_single_cpu_min": round(genecall_min, 2),
        "genecall_batches": genecall_batches,
        "crispr_default_cpu_min": round(crispr_min, 2),
        "total_called_crispr_cpu_min": round(ingest_min + genecall_min + crispr_min, 2),
        "root_gb_lo": round((new_bases / 1e9) * 0.72, 2),
        "root_gb_hi": round((new_bases / 1e9) * 0.85, 2),
    }


def _fmt_opt(val: float | int | None) -> str:
    return "" if val is None else str(val)


def stage_gart(
    stages_tsv: Path,
    tar_gz: Path,
    out_root: Path,
    *,
    stages: Iterable[str] = ("GART1", "GART2"),
    fasta_dir: Path | None = None,
    final_dirs: dict[str, Path] | None = None,
    root: Path | None = None,
    root_stages: Iterable[str] = (),
    owner_inputs: dict[str, Path] | None = None,
    retrieved_at: str = "2026-10-23T19:51:00Z",
    workers: int = 1,
    gart1_screen_max_bp: int = 500_000,
    use_bgzip: bool = False,
    work: Path | None = None,
    threads: int = 2,
    memory_limit: str = "4GB",
) -> dict:
    """Stream GlobDB r226 tar.gz once, screen genetic codes, write consolidated .fna.gz + stage metadata."""
    stage_list = list(stages)
    out_root = Path(out_root)
    for st in stage_list:
        if (out_root / st).exists():
            raise Refused(f"{out_root / st} already exists")

    selected = select_stage_genomes(stages_tsv, stage_list)
    by_id: dict[str, tuple[str, dict]] = {}
    for st in stage_list:
        for row in selected[st]:
            by_id[row["ID"]] = (st, row)

    out_root.mkdir(parents=True, exist_ok=True)
    work_dir = Path(work) if work is not None else out_root / "_work"
    work_dir.mkdir(parents=True, exist_ok=True)
    fa_dir = Path(fasta_dir) if fasta_dir is not None else None
    if fa_dir is not None:
        fa_dir.mkdir(parents=True, exist_ok=True)

    stage_tmp: dict[str, Path] = {}
    writers: dict[str, _GzipWriter] = {}
    fasta_tmp_paths: dict[str, Path] = {}
    fasta_final_paths: dict[str, Path] = {}
    fasta_uncomp_bytes: dict[str, int] = {st: 0 for st in stage_list}

    for st in stage_list:
        tmp_d = out_root / f".{st}.tmp"
        shutil.rmtree(tmp_d, ignore_errors=True)
        tmp_d.mkdir(parents=True)
        stage_tmp[st] = tmp_d
        dst_dir = fa_dir if fa_dir is not None else tmp_d
        fname = STAGE_FASTA[st]
        tmp_fa = dst_dir / f".tmp.{fname}"
        tmp_fa.unlink(missing_ok=True)
        fasta_tmp_paths[st] = tmp_fa
        fasta_final_paths[st] = (fa_dir / fname) if fa_dir is not None else (out_root / st / fname)
        writers[st] = _GzipWriter(tmp_fa, use_bgzip=use_bgzip)

    # Compact per-genome storage: gid -> (concatenated 16-byte contig_uids, array('I', lengths))
    genome_contigs: dict[str, tuple[bytes, array.array]] = {}
    gc_results: dict[str, dict] = {}
    members_seen = 0
    t0 = time.monotonic()
    last_log = t0

    def _drain_futures(pending: list[Future], *, drain_all: bool = False) -> None:
        keep: list[Future] = []
        for fut in pending:
            if drain_all or fut.done():
                gid_done, res_done = fut.result()
                gc_results[gid_done] = res_done
            else:
                keep.append(fut)
        pending[:] = keep
        if not drain_all and pending:
            gid_done, res_done = pending.pop(0).result()
            gc_results[gid_done] = res_done

    pool = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None
    pending: list[Future] = []
    max_in_flight = max(4, workers * 4)

    try:
        tf, close_tf = _open_tar_stream(Path(tar_gz))
        try:
            for m in tf:
                if not m.isfile():
                    continue
                members_seen += 1
                fname = m.name.rsplit("/", 1)[-1]
                if not fname.endswith(".fa.gz"):
                    continue
                gid = fname[:-6]
                hit = by_id.get(gid)
                if hit is None or gid in genome_contigs:
                    continue
                st, meta = hit
                ef = tf.extractfile(m)
                assert ef is not None
                raw_fa = gzip.decompress(ef.read())
                recs = _iter_fasta_bytes(raw_fa)
                rec_bp = sum(len(s) for _, s, _ in recs)
                if len(recs) != meta["contig_number"] or rec_bp != meta["genome_size"]:
                    raise Refused(
                        f"{gid}: FASTA has {len(recs)} records / {rec_bp} bp, "
                        f"table says {meta['contig_number']} / {meta['genome_size']}"
                    )

                cu_buf = bytearray()
                len_arr = array.array("I")
                seqs_for_screen: list[bytes] = []
                w = writers[st]
                for hdr, seq, lw in recs:
                    cu_buf.extend(ids.contig_uid(seq))
                    len_arr.append(len(seq))
                    seqs_for_screen.append(seq)
                    out_bytes = b"\n".join(fasta.format_record(hdr, seq, lw)) + b"\n"
                    w.write(out_bytes)
                    fasta_uncomp_bytes[st] += len(out_bytes)

                genome_contigs[gid] = (bytes(cu_buf), len_arr)
                max_bp = gart1_screen_max_bp if st == "GART1" else 0
                if pool is not None:
                    pending.append(pool.submit(_screen_worker, (gid, seqs_for_screen, meta["taxonomy"], max_bp)))
                    if len(pending) >= max_in_flight:
                        _drain_futures(pending, drain_all=False)
                else:
                    gc_results[gid] = screen_genetic_code(seqs_for_screen, meta["taxonomy"], screen_max_bp=max_bp)

                now = time.monotonic()
                if now - last_log >= 30.0:
                    print(
                        f"[c63] members={members_seen} found={len(genome_contigs)}/{len(by_id)} "
                        f"screened={len(gc_results)} elapsed={now - t0:.1f}s",
                        file=sys.stderr,
                        flush=True,
                    )
                    last_log = now
        finally:
            close_tf()

        if pool is not None:
            _drain_futures(pending, drain_all=True)
            pool.shutdown(wait=True)
            pool = None

        missing = sorted(set(by_id) - set(genome_contigs))
        if missing:
            raise Refused(f"missing {len(missing)} selected genomes in {tar_gz}; first {missing[0]}")

        for st in stage_list:
            writers[st].close()
            if fa_dir is not None:
                fasta_tmp_paths[st].replace(fasta_final_paths[st])
            else:
                fasta_tmp_paths[st].replace(stage_tmp[st] / STAGE_FASTA[st])

        # Build metadata, check DEC-26/DEC-27 ownership, and finalize each stage in order
        reports: dict[str, dict] = {}
        eff_owners = dict(owner_inputs or {})
        r_stages = list(root_stages)

        for st in stage_list:
            tmp_d = stage_tmp[st]
            rows = selected[st]

            # Check for within-stage shared-contig genetic_code conflicts (DEC-27 requires agreement)
            contig_code_map: dict[bytes, tuple[int, str]] = {}
            code_conflicts: list[dict] = []
            resolved_code: dict[str, int] = {r["ID"]: int(gc_results[r["ID"]]["genetic_code"]) for r in rows}

            # First pass: compute genome_uid for each genome so least genome_uid wins if a conflict ever occurs
            genome_uids: dict[str, bytes] = {}
            for r in rows:
                gid = r["ID"]
                cu_raw, _ = genome_contigs[gid]
                uids = [cu_raw[i : i + 16] for i in range(0, len(cu_raw), 16)]
                genome_uids[gid] = ids.genome_uid(list(dict.fromkeys(uids)))

            for r in sorted(rows, key=lambda x: genome_uids[x["ID"]]):
                gid = r["ID"]
                cu_raw, _ = genome_contigs[gid]
                c = resolved_code[gid]
                conflict_target: int | None = None
                for i in range(0, len(cu_raw), 16):
                    cu = cu_raw[i : i + 16]
                    prev = contig_code_map.get(cu)
                    if prev is not None and prev[0] != c:
                        conflict_target = prev[0]
                        code_conflicts.append({
                            "contig_uid": cu.hex(),
                            "first_genome": prev[1],
                            "first_code": prev[0],
                            "second_genome": gid,
                            "second_code": c,
                        })
                if conflict_target is not None:
                    resolved_code[gid] = conflict_target
                    c = conflict_target
                for i in range(0, len(cu_raw), 16):
                    contig_code_map.setdefault(cu_raw[i : i + 16], (c, gid))

            # Write records.tsv and genetic_codes.tsv
            records_path = tmp_d / "records.tsv"
            gc_path = tmp_d / "genetic_codes.tsv"
            with open(records_path, "w", encoding="utf-8", buffering=1 << 20) as rf, open(
                gc_path, "w", encoding="utf-8", buffering=1 << 20
            ) as gf:
                rf.write("\t".join(RECORDS_HEAD) + "\n")
                gf.write("\t".join(GC_HEAD) + "\n")
                ordinal = 0
                for r in rows:
                    gid = r["ID"]
                    cu_raw, len_arr = genome_contigs[gid]
                    code = resolved_code[gid]
                    gc = gc_results[gid]
                    gf.write(
                        "\t".join([
                            gid, st, r["phylum"], r["tax_class"], r["tax_order"], r["resource"],
                            str(r["genome_size"]), str(r["contig_number"]), r["assembly_quality"],
                            str(code), gc["reason"], str(gc["tested_alt"]).lower(),
                            _fmt_opt(gc["cd_11"]), _fmt_opt(gc["mean_orf_11"]), _fmt_opt(gc["n_genes_11"]),
                            _fmt_opt(gc["cd_25"]), _fmt_opt(gc["mean_orf_25"]), _fmt_opt(gc["n_genes_25"]),
                            _fmt_opt(gc["cd_4"]), _fmt_opt(gc["mean_orf_4"]), _fmt_opt(gc["n_genes_4"]),
                            _fmt_opt(gc["tga_per_kb"]),
                        ]) + "\n"
                    )
                    for rec_ord, length in enumerate(len_arr):
                        cu = cu_raw[rec_ord * 16 : (rec_ord + 1) * 16]
                        rf.write(
                            f"{ordinal}\t{gid}\t{rec_ord}\t{cu.hex()}\t{length}\t"
                            f"{r['phylum']}\t{r['resource']}\t{r['assembly_quality']}\t{code}\t{r['c41_stage']}\n"
                        )
                        ordinal += 1

            owners = _owners_from_records(
                records_path, root, r_stages, eff_owners,
                work=work_dir, threads=threads, memory_limit=memory_limit,
            )

            listing = make_inputs.Listing()
            genomes_list: list[dict] = []
            repeated: list[tuple] = []
            shared: list[tuple] = []
            first_of: dict[bytes, str] = {}
            contig_len: dict[bytes, int] = {}
            ordinal = 0
            for r in rows:
                gid = r["ID"]
                cu_raw, len_arr = genome_contigs[gid]
                uids = [cu_raw[i : i + 16] for i in range(0, len(cu_raw), 16)]
                code = resolved_code[gid]
                for rec_ord, (cu, length) in enumerate(zip(uids, len_arr)):
                    contig_len.setdefault(cu, int(length))
                    if cu in owners:
                        shared.append((ordinal, gid, rec_ord, cu.hex(), int(length), owners[cu]))
                    ordinal += 1
                g = listing.add(uids, code, "prokaryote", r["assembly_quality"])
                guid = genome_uids[gid]
                if g is None:
                    repeated.append((ordinal - len(uids), gid, guid.hex(), first_of[guid]))
                    continue
                genomes_list.append(g)
                first_of[guid] = gid

            make_inputs._write_tsv(
                tmp_d / "repeated_records.tsv",
                ["ordinal", "genome", "genome_uid", "first_genome"],
                repeated,
            )
            make_inputs._write_tsv(
                tmp_d / "shared_records.tsv",
                ["ordinal", "genome", "record_ordinal", "contig_uid", "length", "owner_stage"],
                shared,
            )

            final_obj_dir = (final_dirs or {}).get(st) or (fa_dir if fa_dir is not None else (out_root / st))
            objects = [{
                "path": str(final_obj_dir / STAGE_FASTA[st]),
                **GLOBDB_SOURCE_BASE,
                "retrieved_at": retrieved_at,
            }]
            (tmp_d / "objects.json").write_text(json.dumps(objects, indent=1))
            if final_dirs and st in final_dirs:
                staging_objects = [{
                    "path": str(fasta_final_paths[st]),
                    **GLOBDB_SOURCE_BASE,
                    "retrieved_at": retrieved_at,
                }]
                (tmp_d / "objects_staging.json").write_text(json.dumps(staging_objects, indent=1))
            (tmp_d / "genomes.json").write_text(json.dumps(genomes_list))

            stage_out = out_root / st
            tmp_d.rename(stage_out)

            # Register this stage as an earlier owner for subsequent stages (e.g. GART1 before GART2)
            eff_owners[st] = stage_out

            owner_contigs: Counter = Counter()
            owner_bases: Counter = Counter()
            for cu, ow in owners.items():
                if cu in contig_len:
                    owner_contigs[ow] += 1
                    owner_bases[ow] += contig_len[cu]

            new_contigs = sum(1 for cu in contig_len if cu not in owners)
            new_bases = sum(n for cu, n in contig_len.items() if cu not in owners)
            held_bases = sum(owner_bases.values())

            # Per-phylum and per-class/order genetic-code summary
            phylum_counts: dict[str, dict] = defaultdict(lambda: {"genomes": 0, "bp": 0, "contigs": 0})
            tax_code_counts: dict[str, Counter] = defaultdict(Counter)
            for r in rows:
                gid = r["ID"]
                ph = r["phylum"]
                phylum_counts[ph]["genomes"] += 1
                phylum_counts[ph]["bp"] += r["genome_size"]
                phylum_counts[ph]["contigs"] += r["contig_number"]
                key = f"{ph};{r['tax_class']};{r['tax_order']}"
                tax_code_counts[key][str(resolved_code[gid])] += 1

            meta_names = [
                "genomes.json", "objects.json", "records.tsv",
                "repeated_records.tsv", "shared_records.tsv", "genetic_codes.tsv",
            ]
            sha_map = {n: ids.sha256_file(stage_out / n).hex() for n in meta_names}
            sha_map[STAGE_FASTA[st]] = ids.sha256_file(fasta_final_paths[st]).hex()
            bytes_map = {n: (stage_out / n).stat().st_size for n in meta_names}
            bytes_map[STAGE_FASTA[st]] = fasta_final_paths[st].stat().st_size

            reports[st] = {
                "stage": st,
                "phyla": list(STAGE_PHYLA[st]),
                "selected_genomes": len(rows),
                "selected_bp": sum(r["genome_size"] for r in rows),
                "selected_contigs": sum(r["contig_number"] for r in rows),
                "genomes": len(genomes_list),
                "contigs": len(contig_len),
                "contig_bases": sum(contig_len.values()),
                "repeated_genomes": len(repeated),
                "shared_records": len(shared),
                "owner_contigs": dict(sorted(owner_contigs.items())),
                "owner_bases": dict(sorted(owner_bases.items())),
                "new_contigs": new_contigs,
                "new_bases": new_bases,
                "fasta_uncompressed_bytes": fasta_uncomp_bytes[st],
                "codes": {
                    str(c): sum(1 for g in genomes_list if g["genetic_code"] == c)
                    for c in sorted({g["genetic_code"] for g in genomes_list})
                },
                "tested_alt_genomes": sum(1 for r in rows if gc_results[r["ID"]]["tested_alt"]),
                "code_conflicts": code_conflicts,
                "quality": {
                    q: sum(1 for g in genomes_list if g["assembly_quality"] == q)
                    for q in sorted({g["assembly_quality"] for g in genomes_list})
                },
                "by_phylum": dict(sorted(phylum_counts.items())),
                "by_tax_code": {k: dict(sorted(v.items())) for k, v in sorted(tax_code_counts.items())},
                "estimates": _estimate_called_crispr(fasta_uncomp_bytes[st], new_contigs, new_bases, held_bases),
                "sha256": sha_map,
                "bytes": bytes_map,
                **listing.report(),
            }

        if work is None:
            shutil.rmtree(work_dir, ignore_errors=True)
        return reports
    except BaseException:
        if pool is not None:
            pool.shutdown(wait=False, cancel_futures=True)
        for st in stage_list:
            try:
                writers[st].close()
            except Exception:
                pass
            fasta_tmp_paths[st].unlink(missing_ok=True)
            shutil.rmtree(stage_tmp[st], ignore_errors=True)
        raise


def _parse_kv_paths(items: list[str], flag: str) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for item in items:
        if "=" not in item:
            raise Refused(f"{flag} {item!r}: expected NAME=DIR")
        name, path = item.split("=", 1)
        out[name] = Path(path)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Stage GlobDB r226 ART-RT phyla as GART1 and GART2 (card C-63)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("stage")
    p.add_argument("stages_tsv", type=Path)
    p.add_argument("tar_gz", type=Path)
    p.add_argument("out_root", type=Path)
    p.add_argument("--stages", default="GART1,GART2")
    p.add_argument("--fasta-dir", type=Path, default=None)
    p.add_argument("--final-dir", action="append", default=[])
    p.add_argument("--root", type=Path, default=None)
    p.add_argument("--root-stages", default="")
    p.add_argument("--owner-inputs", action="append", default=[])
    p.add_argument("--retrieved-at", default="2026-10-23T19:51:00Z")
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--gart1-screen-max-bp", type=int, default=500_000)
    p.add_argument("--bgzip", action="store_true")
    p.add_argument("--work", type=Path, default=None)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--memory-limit", default="4GB")
    p.add_argument("--report-out", type=Path, default=None)
    a = ap.parse_args(argv)
    try:
        stages = [s for s in a.stages.split(",") if s]
        r_stages = [s for s in a.root_stages.split(",") if s]
        rep = stage_gart(
            a.stages_tsv, a.tar_gz, a.out_root,
            stages=stages,
            fasta_dir=a.fasta_dir,
            final_dirs=_parse_kv_paths(a.final_dir, "--final-dir"),
            root=a.root,
            root_stages=r_stages,
            owner_inputs=_parse_kv_paths(a.owner_inputs, "--owner-inputs"),
            retrieved_at=a.retrieved_at,
            workers=a.workers,
            gart1_screen_max_bp=a.gart1_screen_max_bp,
            use_bgzip=a.bgzip,
            work=a.work,
            threads=a.threads,
            memory_limit=a.memory_limit,
        )
    except Refused as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2
    txt = json.dumps(rep, indent=2, sort_keys=True)
    if a.report_out is not None:
        a.report_out.parent.mkdir(parents=True, exist_ok=True)
        a.report_out.write_text(txt + "\n")
    print(txt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
