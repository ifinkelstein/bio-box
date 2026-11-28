#!/usr/bin/env python3
"""Card C-87 (DEC-80, DEC-88): stream-census and stage pod-resident JGI IUIG and IMG-UMAG prokaryote archives.

Subcommands:
    python c87_jgi_ug_inputs.py census
        --iuig-zip PATH --umag-zip PATH --meta-zip PATH --r232-index DIR --out-dir DIR
        [--gtdb-meta-dir DIR] [--gems-meta-tar PATH] [--img-m-tar PATH]
        [--non-gtdb-stages PATH] [--grare-selection PATH]
        [--max-stage-bp 35000000000] [--workers 2] [--cpu-guard-min 280] [--report-out JSON]

    python c87_jgi_ug_inputs.py stage
        --census-dir DIR --iuig-zip PATH --umag-zip PATH --meta-out DIR --fasta-out DIR
        [--root DIR] [--root-stages S1,...,SPX1] [--chain-home /tmp/if2453/s2/inputs]
        [--later STAGE=DIR ...] [--inputs-home /tmp/if2453/s2/inputs]
        [--workers 4] [--max-tasks-per-child 250] [--screen-max-bp 500000]
        [--cpu-guard-min 280] [--move-large-records-to-stor] [--report-out JSON]

    python c87_jgi_ug_inputs.py write-conf
        --meta-dir DIR --fasta-dir DIR --out-conf PATH --note REL_PATH --after STAGE

Supports both JGI Data Portal nested wrapper `.zip` archives (`zip -> outer .tar -> per-genome <id>.tar.gz`,
skipping inner `.tar.gz` decompression for unselected entries and caching per-metagenome `mbin_datafile` tables
so only one bin bundle per metagenome is opened during census) and direct per-genome `.zip` archives (seeking
in O(1) directly to selected members via `ZipFile.open(info)` without reading unselected entries).
"""

from __future__ import annotations

import argparse
import array
import csv
import gzip
import io
import json
import multiprocessing
import os
import re
import resource
import shutil
import subprocess
import sys
import tarfile
import time
import zipfile
import zlib
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import a26_r26_11_4_chain as a26  # noqa: E402
import a29_r26_11_5_chain as a29  # noqa: E402
import c58_sags_symbionts_niche as c58  # noqa: E402
import c77_sag2_inputs as c77  # noqa: E402
import make_inputs  # noqa: E402
import r232_index  # noqa: E402
from biobox import fasta, ids  # noqa: E402

import pyarrow as pa  # noqa: E402
import pyarrow.csv as pacsv  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

Refused = make_inputs.Refused

STAGE = "JGIUG1"
FASTA_NAME = "jgiug1_genomes.fna.gz"
STOR_FLOOR_BYTES = 64 * 1024 * 1024  # 67,108,864 B (AGENTS.md §7.3)
DEFAULT_MAX_STAGE_BP = 35_000_000_000  # <= 35 Gbp cap for JGIUG1

ROOT_STAGES = a29.ROOT_STAGES
JGIUG_RETRIEVED = "2026-11-19T10:29:56Z"
JGIUG_SOURCE = {
    "source_name": "jgi_iuig_umag",
    "source_release": (
        "JGI Custom_MPI-IUIG (2023-11-14_1.0) and Custom_MPI-IMG-UMAG (2024-02-06), "
        "jgi_iuig_genomes_20261118.zip and jgi_img_umag_prok_bins_20261118.zip (card C-87)"
    ),
    "license": "JGI Data Utilization Policy (dataset marked Unrestricted on portal; accepted by the PI on 2026-11-18)",
    "redistributable": False,
}

IUPAC = b"ACGTURYSWKMBDHVNacgturyswkmbdhvn"
GAPS = b"-."
GCA_RE = re.compile(r"(GC[AF]_\d{9})(?:\.\d+)?")
BIOSAMPLE_RE = re.compile(r"SAM[NED][A-Z]?\d+")
WGS_RE = re.compile(r"([A-Z]{4,6})\d{8,10}(?:\.\d+)?")

NON_SPECIES_SECONDS = frozenset({
    "sp", "sp.", "spp", "spp.", "bacterium", "archaeon", "uncultured", "unclassified",
    "endosymbiont", "symbiont", "metagenome", "isolate", "candidatus", "enrichment",
    "clone", "strain", "group", "cluster", "none", "na",
})

DPANN_ASGARD_MARKERS = (
    "dpann", "nanoarchae", "nanobdell", "micrarchae", "aenigmarchae", "aenigmatarchae",
    "undinarchae", "iainarchae", "diapherotrites", "altiarchae", "huberarchae", "parvarchae",
    "woesearchae", "pacearchae", "asgard", "heimdallarchae", "lokiarchae", "thorarchae", "odinarchae",
)
CPR_MARKERS = (
    "patescibacteri", "cpr", "paceibacteri", "pacebacteria", "parcubacteri", "microgenomat",
    "saccharibacteri", "gracilibacteri", "absconditabacteri", "peregrinibacteri", "dojkabacteri",
    "kazan", "berkelbacteri", "cand. div.", "candidate division", "bd1-5", "sr1", "tm7", "od1", "op11",
)
ART_ENRICHED_PHYLA = frozenset({
    "myxococcota", "myxococcota_a", "cyanobacteriota", "cyanobacteria",
    "omnitrophota", "candidatus omnitrophota", "planctomycetota", "planctomycetes",
})
MAJOR_COMMON_PHYLA = frozenset({
    "pseudomonadota", "proteobacteria",
    "bacillota", "bacillota_a", "bacillota_b", "bacillota_c", "bacillota_d", "bacillota_e",
    "bacillota_f", "bacillota_g", "bacillota_h", "bacillota_i",
    "firmicutes", "firmicutes_a", "firmicutes_b", "firmicutes_c", "firmicutes_d", "firmicutes_e",
    "firmicutes_f", "firmicutes_g", "firmicutes_h", "firmicutes_i",
    "bacteroidota", "bacteroidota_a", "bacteroidetes",
    "actinomycetota", "actinobacteriota", "actinobacteria",
    "chloroflexota", "chloroflexi",
    "acidobacteriota", "acidobacteria",
    "verrucomicrobiota", "verrucomicrobia",
    "desulfobacterota", "desulfobacterota_a", "desulfobacterota_b", "desulfobacterota_g", " desulfobacterota_i",
    "gemmatimonadota", "gemmatimonadetes",
    "campylobacterota", "epsilonbacteraeota",
    "spirochaetota", "spirochaetes",
})

SELECTION_HEAD = [
    "rank", "genome_id", "source", "metagenome_id", "domain", "phylum", "lineage",
    "quality_tier", "assembly_quality", "completeness", "contamination", "quality_score",
    "bp", "contigs", "tgz_bytes", "lineage_group", "novelty_tier",
    "r232_status", "gtdb_r232_acc", "globdb_r232_id", "held_globdb_stage", "ecosystem",
]
LEFT_HEAD = [*SELECTION_HEAD, "defer_reason"]
RECORDS_HEAD = [
    "ordinal", "genome", "source", "record_ordinal", "contig_uid", "length",
    "domain", "phylum", "quality_tier", "assembly_quality", "genetic_code", "r232_status",
]
GC_HEAD = [
    "genome", "source", "domain", "phylum", "lineage", "bp", "contigs",
    "assembly_quality", "genetic_code", "reason", "tested_alt",
    "cd_11", "mean_orf_11", "n_genes_11",
    "cd_25", "mean_orf_25", "n_genes_25",
    "cd_4", "mean_orf_4", "n_genes_4",
    "tga_per_kb",
]
STATUS_HEAD = [
    "genome", "source", "domain", "phylum", "quality_tier", "bp", "contigs",
    "r232_status", "gtdb_r232_acc", "globdb_r232_id", "held_globdb_stage", "lineage",
]
LEFT_OUT_HEAD = ["genome_id", "source", "records", "bp", "reason"]


def _cpu_min_self() -> float:
    ru = resource.getrusage(resource.RUSAGE_SELF)
    return (ru.ru_utime + ru.ru_stime) / 60.0


def _peak_rss_bytes() -> int:
    try:
        for ln in Path("/proc/self/status").read_text().splitlines():
            if ln.startswith("VmHWM:"):
                return int(ln.split()[1]) * 1024
    except OSError:
        pass
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss if sys.platform == "darwin" else rss * 1024


def _stem_gca(val: str) -> str:
    m = GCA_RE.search(val or "")
    return m.group(1) if m else ""


def _wgs_proj(val: str) -> str:
    m = WGS_RE.fullmatch((val or "").strip())
    return m.group(1) if m else ""


def _clean_binomial(name: str) -> str:
    s = (name or "").strip()
    if s.startswith("s__"):
        s = s[3:].strip()
    if not s or s.lower() in ("none", "unclassified", "na"):
        return ""
    toks = s.split()
    if tok0_cand := (toks[0].lower() == "candidatus"):
        if len(toks) < 3:
            return ""
        g, sp = f"{toks[0]} {toks[1]}", toks[2]
    else:
        if len(toks) < 2:
            return ""
        g, sp = toks[0], toks[1]
    if sp.lower().rstrip(".") in NON_SPECIES_SECONDS or any(ch.isdigit() for ch in sp):
        return ""
    return f"{g} {sp}"


def _parse_lineage_domain_phylum(lineage: str, fallback_domain: str = "", fallback_phylum: str = "") -> tuple[str, str, str]:
    """Return (domain, phylum, species_token) from GTDB-Tk or IMG semicolon lineage."""
    parts = [p.strip() for p in (lineage or "").split(";") if p.strip()]
    dom = fallback_domain
    phy = fallback_phylum
    sp = ""
    if any("__" in p for p in parts):
        for p in parts:
            if p.startswith("d__") and len(p) > 3:
                dom = p[3:].strip()
            elif p.startswith("p__") and len(p) > 3:
                phy = p[3:].strip()
            elif p.startswith("s__") and len(p) > 3:
                sp = p[3:].strip()
    elif parts:
        if parts[0] in ("Bacteria", "Archaea", "Eukaryota", "Viruses"):
            dom = parts[0]
        if len(parts) >= 2 and parts[1] not in ("None", "unclassified", "Unclassified"):
            phy = parts[1]
        if len(parts) >= 7 and parts[6] not in ("None", "unclassified", "Unclassified"):
            sp = parts[6]
    return dom or fallback_domain or "Unclassified", phy or fallback_phylum or "Unclassified", sp


def _classify_lineage_group(domain: str, phylum: str, lineage: str, ecosystem: str) -> tuple[int, str]:
    low = f"{domain};{phylum};{lineage}".lower()
    if domain == "Archaea" or "d__archaea" in low:
        if any(m in low for m in DPANN_ASGARD_MARKERS):
            return 0, "archaea_dpann_asgard"
        return 1, "archaea_other"
    if any(m in low for m in CPR_MARKERS):
        return 2, "cpr_patescibacteria"
    p_low = (phylum or "").strip().lower()
    if p_low.startswith("p__"):
        p_low = p_low[3:]
    if not p_low or p_low in ("unclassified", "none", "unknown"):
        return 3, "bacteria_unclassified_phylum"
    if p_low not in MAJOR_COMMON_PHYLA:
        if p_low in ART_ENRICHED_PHYLA:
            return 5, "bacteria_art_enriched_phylum"
        return 4, "bacteria_rare_or_candidate_phylum"
    eco_low = (ecosystem or "").lower()
    if any(k in eco_low for k in ("human", "clinical", "built environment")):
        return 7, "bacteria_major_phylum_host_human"
    return 6, "bacteria_major_phylum_environmental"


def _member_genome_id(member_name: str) -> str:
    base = member_name.rsplit("/", 1)[-1]
    for suf in (".tar.gz", ".tgz", ".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"):
        if base.endswith(suf):
            return base[: -len(suf)]
    return base


def inspect_zip_entries(zip_path: Path) -> dict:
    """Read ONLY the ZIP central directory (`ZipFile.infolist()`) without decompressing any member."""
    zip_path = Path(zip_path)
    with zipfile.ZipFile(zip_path, "r") as zf:
        infos = zf.infolist()

    files = [i for i in infos if not i.is_dir() and not i.filename.startswith("__MACOSX/")]
    outer_tars = [i for i in files if i.filename.endswith(".tar")]
    manifests = [i for i in files if i.filename.endswith((".csv", ".txt", ".tsv", ".md"))]
    genome_members = [
        i for i in files
        if i.filename.endswith((".tar.gz", ".tgz", ".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"))
    ]

    if len(outer_tars) == 1 and not genome_members:
        tinfo = outer_tars[0]
        return {
            "zip_path": str(zip_path),
            "zip_bytes": zip_path.stat().st_size,
            "layout": "nested_tar",
            "zip_entries": len(infos),
            "inner_tar": tinfo.filename,
            "inner_tar_uncompressed_bytes": tinfo.file_size,
            "inner_tar_compressed_bytes": tinfo.compress_size,
            "inner_tar_compress_type": tinfo.compress_type,
            "manifest_files": [m.filename for m in manifests],
        }

    return {
        "zip_path": str(zip_path),
        "zip_bytes": zip_path.stat().st_size,
        "layout": "direct_members",
        "zip_entries": len(infos),
        "genome_member_count": len(genome_members),
        "genome_uncompressed_bytes": sum(i.file_size for i in genome_members),
        "genome_compressed_bytes": sum(i.compress_size for i in genome_members),
        "manifest_files": [m.filename for m in manifests],
    }


def _read_fasta_stream(fh, label: str) -> list[tuple[str, bytes, int]]:
    """Parse FASTA records `(header, seq_bytes, line_width)` from binary stream `fh`."""
    out: list[tuple[str, bytes, int]] = []
    header: str | None = None
    lines: list[bytes] = []
    for raw in fh:
        ln = raw.rstrip(b"\r\n")
        if not ln:
            continue
        if ln.startswith(b">"):
            if header is not None:
                out.append((header, b"".join(lines), fasta._width(lines) if lines else 80))
            try:
                header = ln[1:].decode("utf-8")
            except UnicodeDecodeError as exc:
                raise ValueError(f"{label}: header_not_utf8: {exc}") from exc
            lines = []
        elif header is None:
            raise ValueError(f"{label}: sequence_before_first_header")
        else:
            lines.append(ln)
    if header is not None:
        out.append((header, b"".join(lines), fasta._width(lines) if lines else 80))
    return out


def _extract_fna_from_inner_tgz(fileobj, gid: str) -> list[tuple[str, bytes, int]]:
    with tarfile.open(fileobj=fileobj, mode="r|gz") as itf:
        for sub in itf:
            if not sub.isfile():
                continue
            sname = sub.name.rsplit("/", 1)[-1]
            if sname.endswith((".fna", ".fa", ".fasta")) and not sname.endswith((".genes.fna", ".intergenic.fna")):
                sf = itf.extractfile(sub)
                if sf is not None:
                    return _read_fasta_stream(sf, f"{gid}:{sub.name}")
    return []


def iter_selected_from_zip(
    zip_path: Path,
    selected_ids: set[str],
    *,
    cpu_guard_min: float = 280.0,
):
    """Yield `(genome_id, [(header, seq, width), ...])` for `selected_ids` ONLY.

    - In `direct_members` layout: seeks in O(1) directly to selected `ZipInfo` entries via `zf.open(info)`,
      never touching unselected entries.
    - In `nested_tar` layout: streams the outer `.tar` sequentially (`mode="r|"`) and skips unselected
      inner `<id>.tar.gz` bundles without decompressing them.
    """
    if not selected_ids:
        return
    zip_path = Path(zip_path)
    info = inspect_zip_entries(zip_path)
    remaining = set(selected_ids)

    with zipfile.ZipFile(zip_path, "r") as zf:
        if info["layout"] == "direct_members":
            by_gid: dict[str, zipfile.ZipInfo] = {}
            for zi in zf.infolist():
                if zi.is_dir() or zi.filename.startswith("__MACOSX/"):
                    continue
                if zi.filename.endswith((".tar.gz", ".tgz", ".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta")):
                    gid = _member_genome_id(zi.filename)
                    if gid in remaining:
                        by_gid[gid] = zi
            # Read in archive central-directory order for sequential disk locality and determinism
            for zi in zf.infolist():
                gid = _member_genome_id(zi.filename)
                if by_gid.get(gid) is not zi:
                    continue
                if _cpu_min_self() > cpu_guard_min:
                    raise Refused(f"CPU guard exceeded ({_cpu_min_self():.1f} > {cpu_guard_min} CPU-min)")
                with zf.open(zi, "r") as raw:
                    if zi.filename.endswith((".tar.gz", ".tgz")):
                        recs = _extract_fna_from_inner_tgz(raw, gid)
                    elif zi.filename.endswith(".gz"):
                        with gzip.GzipFile(fileobj=raw, mode="rb") as gz:
                            recs = _read_fasta_stream(gz, gid)
                    else:
                        recs = _read_fasta_stream(raw, gid)
                yield gid, recs
                remaining.discard(gid)
                if not remaining:
                    break
        else:
            inner_tar = info["inner_tar"]
            with zf.open(inner_tar, "r") as outer_fh:
                with tarfile.open(fileobj=outer_fh, mode="r|") as tf:
                    for m in tf:
                        if not m.isfile():
                            continue
                        if _cpu_min_self() > cpu_guard_min:
                            raise Refused(f"CPU guard exceeded ({_cpu_min_self():.1f} > {cpu_guard_min} CPU-min)")
                        gid = _member_genome_id(m.name)
                        if gid not in remaining:
                            continue
                        ef = tf.extractfile(m)
                        if ef is None:
                            continue
                        if m.name.endswith((".tar.gz", ".tgz")):
                            recs = _extract_fna_from_inner_tgz(ef, gid)
                        elif m.name.endswith(".gz"):
                            with gzip.GzipFile(fileobj=ef, mode="rb") as gz:
                                recs = _read_fasta_stream(gz, gid)
                        else:
                            recs = _read_fasta_stream(ef, gid)
                        yield gid, recs
                        remaining.discard(gid)
                        if not remaining:
                            break


def _parse_all_mbin_rows(raw: bytes) -> dict[str, dict]:
    """Parse every bin row in `mbin_datafile_<taxon_oid>.txt` so one read serves the entire metagenome."""
    out: dict[str, dict] = {}
    header: list[str] | None = None
    for line in raw.decode("utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if header is None:
            header = [c.strip() for c in cols]
            continue
        row = dict(zip(header, cols))
        bid = (row.get("IMG Bin ID") or cols[0]).strip()
        if not bid:
            continue
        try:
            comp = float(row.get("Bin Completeness") or 0.0)
        except ValueError:
            comp = 0.0
        cont_s = row.get("Bin Contamintation") or row.get("Bin Contamination") or "0.0"
        try:
            cont = float(cont_s)
        except ValueError:
            cont = 0.0
        try:
            bp = int(float(row.get("Total Number of Bases") or 0))
        except ValueError:
            bp = 0
        scaf_str = (row.get("IMG scaffold IDs of members") or "").strip()
        n_scaf = (scaf_str.count(",") + 1) if scaf_str else 0
        out[bid] = {
            "bin_quality": (row.get("Bin Quality") or "").strip(),
            "bin_lineage": (row.get("Bin Lineage") or "").strip(),
            "gtdb_lineage": (row.get("GTDB-TK lineage") or "").strip(),
            "completeness": comp,
            "contamination": cont,
            "bp": bp,
            "contigs": n_scaf,
        }
    return out


def _load_r232_crosswalk(
    r232_index_dir: Path,
    gtdb_meta_dir: Path | None = None,
    gems_meta_tar: Path | None = None,
    img_m_tar: Path | None = None,
    non_gtdb_stages_tsv: Path | None = None,
    grare_selection_tsv: Path | None = None,
) -> dict:
    """Load C-74 r232 index + GTDB R232 metadata + GEMs OTU/CheckM + IMG/M 70col + held GlobDB stages."""
    r232_index_dir = Path(r232_index_dir)
    gtdb_tbl = pq.read_table(
        r232_index_dir / "gtdb_r232.parquet",
        columns=["accession", "genbank_accession", "species", "is_representative", "representative"],
    )
    gtdb_by_acc: dict[str, str] = {}
    gtdb_species_set: set[str] = set()
    for acc, gb, sp in zip(
        gtdb_tbl.column("accession").to_pylist(),
        gtdb_tbl.column("genbank_accession").to_pylist(),
        gtdb_tbl.column("species").to_pylist(),
    ):
        if acc:
            gtdb_by_acc[acc] = acc
        if gb:
            gtdb_by_acc.setdefault(gb, acc)
        if sp:
            gtdb_species_set.add(sp)
            b = _clean_binomial(sp)
            if b:
                gtdb_species_set.add(b)

    glob_tbl = pq.read_table(
        r232_index_dir / "globdb_r232.parquet",
        columns=["globdb_id", "original_id", "cluster_id", "insdc_accession"],
    )
    glob_ids: set[str] = set()
    glob_by_acc: dict[str, str] = {}
    glob_by_orig: dict[str, str] = {}
    glob_by_cluster: dict[str, str] = {}
    for gid, orig, cl, ins in zip(
        glob_tbl.column("globdb_id").to_pylist(),
        glob_tbl.column("original_id").to_pylist(),
        glob_tbl.column("cluster_id").to_pylist(),
        glob_tbl.column("insdc_accession").to_pylist(),
    ):
        glob_ids.add(gid)
        if ins:
            glob_by_acc.setdefault(ins, gid)
        if orig:
            glob_by_orig.setdefault(orig, gid)
        if cl:
            glob_by_cluster.setdefault(cl, gid)

    dict_tbl = pq.read_table(
        r232_index_dir / "globdb_r232_dictionary.parquet",
        columns=["globdb_id", "original_id", "cluster_id"],
    )
    dict_orig: dict[str, str] = {}
    dict_cluster: dict[str, str] = {}
    for gid, orig, cl in zip(
        dict_tbl.column("globdb_id").to_pylist(),
        dict_tbl.column("original_id").to_pylist(),
        dict_tbl.column("cluster_id").to_pylist(),
    ):
        if orig:
            dict_orig.setdefault(orig, gid)
        if cl:
            dict_cluster.setdefault(cl, gid)

    # Held GlobDB r226 stages (G1, GART1, GART2, GRARE1..3) mapped through globdb_r226_keys.parquet
    r226_stage_by_id: dict[str, str] = {}
    if non_gtdb_stages_tsv is not None and Path(non_gtdb_stages_tsv).is_file():
        with open(non_gtdb_stages_tsv, encoding="utf-8", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rid = r["ID"]
                st = r["stage"]
                phy = r232_index._lineage_field(r["taxonomy"], "p__") or ""
                if st == "G1":
                    r226_stage_by_id[rid] = "G1"
                elif phy == "Patescibacteriota":
                    r226_stage_by_id[rid] = "GART2"
                elif phy in r232_index.GART_PHYLA:
                    r226_stage_by_id[rid] = "GART1"
    if grare_selection_tsv is not None and Path(grare_selection_tsv).is_file():
        with open(grare_selection_tsv, encoding="utf-8", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                r226_stage_by_id[r["genome"]] = r["stage"]

    held_globdb_by_orig_or_cluster: dict[str, str] = {}
    keys_pq = r232_index_dir / "globdb_r226_keys.parquet"
    if r226_stage_by_id and keys_pq.is_file():
        kt = pq.read_table(keys_pq, columns=["r226_id", "original_id", "cluster_id"])
        for rid, orig, cl in zip(
            kt.column("r226_id").to_pylist(),
            kt.column("original_id").to_pylist(),
            kt.column("cluster_id").to_pylist(),
        ):
            hst = r226_stage_by_id.get(rid)
            if hst:
                if orig:
                    held_globdb_by_orig_or_cluster.setdefault(orig, hst)
                if cl:
                    held_globdb_by_orig_or_cluster.setdefault(cl, hst)

    gtdb_by_bs: dict[str, str] = {}
    gtdb_by_wgs: dict[str, str] = {}
    gtdb_by_taxid: dict[str, str] = {}
    if gtdb_meta_dir is not None and Path(gtdb_meta_dir).is_dir():
        want_cols = ["accession", "ncbi_biosample", "ncbi_wgs_master", "ncbi_taxid", "ncbi_organism_name"]
        for key in ("ar53", "bac120"):
            p = Path(gtdb_meta_dir) / r232_index.FILES[key]
            if not p.is_file():
                continue
            with gzip.open(p, "rt", encoding="utf-8", errors="replace") as gz:
                avail = set(gz.readline().rstrip("\r\n").split("\t"))
            use = [c for c in want_cols if c in avail]
            if "accession" not in use:
                continue
            rd = pacsv.open_csv(
                str(p),
                read_options=pacsv.ReadOptions(block_size=1 << 24),
                parse_options=pacsv.ParseOptions(delimiter="\t"),
                convert_options=pacsv.ConvertOptions(
                    include_columns=use,
                    column_types={c: pa.string() for c in use},
                ),
            )
            for b in rd:
                cols_py = {c: b.column(c).to_pylist() for c in use}
                n_rows = len(cols_py["accession"])
                for idx_r in range(n_rows):
                    acc = r232_index.nover(cols_py["accession"][idx_r] or "")
                    bs = (cols_py["ncbi_biosample"][idx_r] or "").strip() if "ncbi_biosample" in cols_py else ""
                    if bs and bs != "none":
                        gtdb_by_bs.setdefault(bs, acc)
                    wgs = _wgs_proj(cols_py["ncbi_wgs_master"][idx_r] or "") if "ncbi_wgs_master" in cols_py else ""
                    if wgs:
                        gtdb_by_wgs.setdefault(wgs, acc)
                    tid = (cols_py["ncbi_taxid"][idx_r] or "").strip() if "ncbi_taxid" in cols_py else ""
                    if tid and tid != "none":
                        gtdb_by_taxid.setdefault(tid, acc)
                    org = (cols_py["ncbi_organism_name"][idx_r] or "").strip() if "ncbi_organism_name" in cols_py else ""
                    bname = _clean_binomial(org)
                    if bname:
                        gtdb_species_set.add(bname)

    gem_otu_of: dict[str, str] = {}
    gem_checkm: dict[str, tuple[float, float]] = {}
    if gems_meta_tar is not None and Path(gems_meta_tar).is_file():
        with tarfile.open(gems_meta_tar, "r") as tf:
            names = {m.name: m for m in tf.getmembers() if m.isfile()}
            if "gems/otus/genome_to_otu.tsv" in names:
                ef = tf.extractfile(names["gems/otus/genome_to_otu.tsv"])
                assert ef is not None
                for r in csv.DictReader(io.TextIOWrapper(ef, encoding="utf-8", errors="replace"), delimiter="\t"):
                    gid = (r.get("genome_id") or "").strip()
                    oid = (r.get("otu_id") or "").strip()
                    if gid and oid:
                        gem_otu_of[gid] = oid
            if "gems/otus/genome_to_checkm.tsv" in names:
                ef = tf.extractfile(names["gems/otus/genome_to_checkm.tsv"])
                assert ef is not None
                for r in csv.DictReader(io.TextIOWrapper(ef, encoding="utf-8", errors="replace"), delimiter="\t"):
                    gid = (r.get("genome_id") or "").strip()
                    if gid:
                        gem_checkm[gid] = (float(r["completeness"]), float(r["contamination"]))

    img_meta: dict[str, dict] = {}
    if img_m_tar is not None and Path(img_m_tar).is_file():
        with tarfile.open(img_m_tar, "r") as tf:
            for m in tf.getmembers():
                if not m.isfile() or not m.name.endswith(".tsv.gz"):
                    continue
                ef = tf.extractfile(m)
                if ef is None:
                    continue
                with gzip.open(ef, "rt", encoding="utf-8", errors="replace") as gz:
                    rd = csv.DictReader(gz, delimiter="\t")
                    fnames = rd.fieldnames or []
                    if "IMG Genome ID" not in fnames:
                        continue
                    eco_cols = [c for c in ("Ecosystem", "Ecosystem Category", "Ecosystem Type") if c in fnames]
                    for r in rd:
                        oid = (r.get("IMG Genome ID") or "").strip()
                        if not oid:
                            continue
                        gca_s = _stem_gca(r.get("NCBI Assembly Accession") or "")
                        bs_s = (r.get("NCBI Biosample Accession") or "").strip()
                        gsize_s = (r.get("Genome Size * assembled") or "").strip()
                        scaf_s = (r.get("Scaffold Count * assembled") or "").strip()
                        unc_s = (r.get("Uncultured Type") or "").strip()
                        eco_s = " / ".join((r.get(c) or "").strip() for c in eco_cols if (r.get(c) or "").strip())
                        img_meta[oid] = {
                            "gca": gca_s,
                            "bs": bs_s if BIOSAMPLE_RE.fullmatch(bs_s) else "",
                            "bp": int(gsize_s) if gsize_s.isdigit() else 0,
                            "contigs": int(scaf_s) if scaf_s.isdigit() else 0,
                            "uncultured_type": unc_s,
                            "ecosystem": eco_s,
                        }

    return {
        "gtdb_by_acc": gtdb_by_acc,
        "gtdb_by_bs": gtdb_by_bs,
        "gtdb_by_wgs": gtdb_by_wgs,
        "gtdb_by_taxid": gtdb_by_taxid,
        "gtdb_species_set": gtdb_species_set,
        "glob_ids": glob_ids,
        "glob_by_acc": glob_by_acc,
        "glob_by_orig": glob_by_orig,
        "glob_by_cluster": glob_by_cluster,
        "dict_orig": dict_orig,
        "dict_cluster": dict_cluster,
        "held_globdb_by_orig_or_cluster": held_globdb_by_orig_or_cluster,
        "gem_otu_of": gem_otu_of,
        "gem_checkm": gem_checkm,
        "img_meta": img_meta,
    }


def _classify_r232(
    genome_id: str,
    *,
    gca: str,
    biosample: str,
    taxid: str,
    species_str: str,
    refs: dict,
) -> tuple[str, str, str, str, int]:
    """Return `(r232_status, gtdb_acc, globdb_id, held_globdb_stage, novelty_tier)`.

    `novelty_tier`:
      0 = completely absent from GTDB R232, GlobDB r232, GEM OTUs, and GTDB species binomials (unclassified/novel species)
      1 = absent from GTDB R232 & GlobDB r232 accessions/OTUs, carries a non-GTDB species string
      2 = absent from GTDB R232 & GlobDB r232 accessions/OTUs, but species string/taxid matches a known GTDB R232 species
      9 = directly in GTDB R232, GlobDB r232, or collapsed GEM OTU
    """
    gt_acc = ""
    if gca:
        st = _stem_gca(gca)
        gt_acc = refs["gtdb_by_acc"].get(st, "")
    if not gt_acc and biosample:
        gt_acc = refs["gtdb_by_bs"].get(biosample, "")
    if not gt_acc and species_str:
        st_sp = _stem_gca(species_str)
        if st_sp:
            gt_acc = refs["gtdb_by_acc"].get(st_sp, "")

    otu = refs["gem_otu_of"].get(genome_id, "")
    held_st = (
        refs["held_globdb_by_orig_or_cluster"].get(genome_id, "")
        or (refs["held_globdb_by_orig_or_cluster"].get(otu, "") if otu else "")
    )

    gl_id = (
        refs["glob_by_orig"].get(genome_id, "")
        or (refs["glob_by_cluster"].get(otu, "") if otu else "")
        or (refs["glob_by_acc"].get(_stem_gca(gca), "") if gca else "")
        or (refs["glob_by_acc"].get(gt_acc, "") if gt_acc else "")
    )
    if gt_acc and gt_acc in refs["glob_ids"] and not gl_id:
        gl_id = gt_acc

    in_gt = bool(gt_acc)
    in_gl = bool(gl_id)
    if in_gt and in_gl:
        return "in_gtdb_and_globdb_r232", gt_acc, gl_id, held_st, 9
    if in_gt:
        return "in_gtdb_r232", gt_acc, "", held_st, 9
    if in_gl:
        return "in_globdb_r232", "", gl_id, held_st, 9

    # Check if evaluated and collapsed by GlobDB r232's GEMOTU / other dictionaries
    collapsed_gid = (
        refs["dict_orig"].get(genome_id, "")
        or (refs["dict_cluster"].get(otu, "") if otu else "")
    )
    if collapsed_gid or otu:
        return "globdb_r232_collapsed_gem_otu", "", collapsed_gid or otu, held_st, 9

    # Taxid or species binomial match to GTDB R232 (strain of a represented species)
    if taxid and taxid in refs["gtdb_by_taxid"]:
        return "gtdb_r232_taxid_match", refs["gtdb_by_taxid"][taxid], "", held_st, 2
    binom = _clean_binomial(species_str)
    if binom and binom in refs["gtdb_species_set"]:
        return "gtdb_r232_species_match", "", "", held_st, 2

    if gca or biosample:
        return "not_in_r232", "", "", held_st, (1 if binom else 0)
    return "not_in_r232", "", "", held_st, (1 if binom else 0)


def _scan_iuig_archive_worker(args: tuple) -> dict:
    """Scan IUIG zip archive: list all member IDs & tgz sizes, and measure exact bp/contigs for novel prok candidates."""
    iuig_zip, measure_fna_ids, cpu_guard_min = Path(args[0]), set(args[1]), float(args[2])
    info = inspect_zip_entries(iuig_zip)
    members: dict[str, tuple[int, int, int]] = {}  # gid -> (tgz_bytes, fna_bp, fna_contigs)
    other_members = 0

    with zipfile.ZipFile(iuig_zip, "r") as zf:
        if info["layout"] == "direct_members":
            for zi in zf.infolist():
                if zi.is_dir() or zi.filename.startswith("__MACOSX/"):
                    continue
                if not zi.filename.endswith((".tar.gz", ".tgz", ".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta")):
                    other_members += 1
                    continue
                gid = _member_genome_id(zi.filename)
                bp, n_c = 0, 0
                if gid in measure_fna_ids:
                    with zf.open(zi, "r") as raw:
                        if zi.filename.endswith((".tar.gz", ".tgz")):
                            recs = _extract_fna_from_inner_tgz(raw, gid)
                        elif zi.filename.endswith(".gz"):
                            with gzip.GzipFile(fileobj=raw, mode="rb") as gz:
                                recs = _read_fasta_stream(gz, gid)
                        else:
                            recs = _read_fasta_stream(raw, gid)
                    bp = sum(len(s) for _, s, _ in recs)
                    n_c = len(recs)
                members[gid] = (zi.file_size, bp, n_c)
        else:
            t_start = time.monotonic()
            t_last = t_start
            with zf.open(info["inner_tar"], "r") as outer_fh:
                with tarfile.open(fileobj=outer_fh, mode="r|") as tf:
                    for m in tf:
                        if not m.isfile():
                            continue
                        if _cpu_min_self() > cpu_guard_min:
                            raise Refused(f"IUIG census CPU guard exceeded ({_cpu_min_self():.1f} > {cpu_guard_min})")
                        if not m.name.endswith((".tar.gz", ".tgz", ".fna.gz", ".fna")):
                            other_members += 1
                            continue
                        gid = _member_genome_id(m.name)
                        bp, n_c = 0, 0
                        if gid in measure_fna_ids:
                            ef = tf.extractfile(m)
                            if ef is not None:
                                if m.name.endswith((".tar.gz", ".tgz")):
                                    recs = _extract_fna_from_inner_tgz(ef, gid)
                                elif m.name.endswith(".gz"):
                                    with gzip.GzipFile(fileobj=ef, mode="rb") as gz:
                                        recs = _read_fasta_stream(gz, gid)
                                else:
                                    recs = _read_fasta_stream(ef, gid)
                                bp = sum(len(s) for _, s, _ in recs)
                                n_c = len(recs)
                        members[gid] = (m.size, bp, n_c)
                        now = time.monotonic()
                        if now - t_last >= 60.0:
                            print(
                                f"[c87-iuig-census] seen={len(members)} measured_fna={sum(1 for _, b, _ in members.values() if b > 0)} "
                                f"cpu_min={_cpu_min_self():.2f} elapsed={now - t_start:.0f}s",
                                file=sys.stderr,
                                flush=True,
                            )
                            t_last = now

    return {
        "zip_info": info,
        "members": members,
        "other_members": other_members,
        "cpu_min": round(_cpu_min_self(), 3),
        "peak_rss_bytes": _peak_rss_bytes(),
    }


def _scan_umag_archive_worker(args: tuple) -> dict:
    """Scan uMAG zip archive: cache `mbin_datafile_<mg_id>.txt` once per metagenome, skipping >90% of inner .tar.gz."""
    umag_zip, cpu_guard_min = Path(args[0]), float(args[1])
    info = inspect_zip_entries(umag_zip)
    bins: dict[str, dict] = {}
    mbin_cache: dict[str, dict[str, dict]] = {}
    mbin_reads = 0
    other_members = 0

    def _read_mbin_from_tgz(fileobj) -> dict[str, dict]:
        with tarfile.open(fileobj=fileobj, mode="r|gz") as itf:
            for sub in itf:
                if not sub.isfile():
                    continue
                sname = sub.name.rsplit("/", 1)[-1]
                if sname.startswith("mbin_datafile_") and sname.endswith(".txt"):
                    sf = itf.extractfile(sub)
                    if sf is not None:
                        return _parse_all_mbin_rows(sf.read())
        return {}

    with zipfile.ZipFile(umag_zip, "r") as zf:
        if info["layout"] == "direct_members":
            for zi in zf.infolist():
                if zi.is_dir() or zi.filename.startswith("__MACOSX/"):
                    continue
                if not zi.filename.endswith((".tar.gz", ".tgz", ".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta")):
                    other_members += 1
                    continue
                bid = _member_genome_id(zi.filename)
                mg_id = bid.split("_", 1)[0] if "_" in bid else bid
                meta = mbin_cache.get(mg_id, {}).get(bid)
                if meta is None and zi.filename.endswith((".tar.gz", ".tgz")):
                    with zf.open(zi, "r") as raw:
                        parsed = _read_mbin_from_tgz(raw)
                    mbin_reads += 1
                    mbin_cache.setdefault(mg_id, {}).update(parsed)
                    meta = mbin_cache[mg_id].get(bid)
                bins[bid] = {"tgz_bytes": zi.file_size, "meta": meta}
        else:
            t_start = time.monotonic()
            t_last = t_start
            with zf.open(info["inner_tar"], "r") as outer_fh:
                with tarfile.open(fileobj=outer_fh, mode="r|") as tf:
                    for m in tf:
                        if not m.isfile():
                            continue
                        if _cpu_min_self() > cpu_guard_min:
                            raise Refused(f"uMAG census CPU guard exceeded ({_cpu_min_self():.1f} > {cpu_guard_min})")
                        if not m.name.endswith((".tar.gz", ".tgz")):
                            other_members += 1
                            continue
                        bid = _member_genome_id(m.name)
                        mg_id = bid.split("_", 1)[0] if "_" in bid else bid
                        meta = mbin_cache.get(mg_id, {}).get(bid)
                        if meta is None:
                            ef = tf.extractfile(m)
                            if ef is not None:
                                parsed = _read_mbin_from_tgz(ef)
                                mbin_reads += 1
                                mbin_cache.setdefault(mg_id, {}).update(parsed)
                                meta = mbin_cache[mg_id].get(bid)
                        bins[bid] = {"tgz_bytes": m.size, "meta": meta}
                        now = time.monotonic()
                        if now - t_last >= 60.0:
                            print(
                                f"[c87-umag-census] bins={len(bins)} mgs={len(mbin_cache)} mbin_reads={mbin_reads} "
                                f"cpu_min={_cpu_min_self():.2f} elapsed={now - t_start:.0f}s",
                                file=sys.stderr,
                                flush=True,
                            )
                            t_last = now

    return {
        "zip_info": info,
        "bins": bins,
        "mbin_tarballs_opened": mbin_reads,
        "metagenomes_seen": len(mbin_cache),
        "other_members": other_members,
        "cpu_min": round(_cpu_min_self(), 3),
        "peak_rss_bytes": _peak_rss_bytes(),
    }


def _rank_and_partition_candidates(
    candidates: list[dict],
    *,
    max_stage_bp: int = DEFAULT_MAX_STAGE_BP,
) -> tuple[list[dict], list[tuple[dict, str]], int]:
    """Deterministic novelty + lineage + quality ranking and cap partitioning for JGIUG1."""

    def _sort_key(c: dict) -> tuple:
        is_redundant = 1 if (c["r232_status"] != "not_in_r232" or bool(c["held_globdb_stage"])) else 0
        if c["quality_tier"] == "HQ":
            q_ord = 0
        elif c.get("_quality_scored", True):
            q_ord = 1
        else:
            q_ord = 2
        return (
            is_redundant,
            c["novelty_tier"],
            c["lineage_group_num"],
            q_ord,
            -c["quality_score"],
            -c["bp"],
            c["genome_id"],
        )

    ordered = sorted(candidates, key=_sort_key)
    selected: list[dict] = []
    left_for_later: list[tuple[dict, str]] = []
    acc_bp = 0

    for rank_idx, c in enumerate(ordered, 1):
        c["rank"] = rank_idx
        is_redundant = c["r232_status"] != "not_in_r232" or bool(c["held_globdb_stage"])
        if is_redundant:
            reason = f"held_{c['held_globdb_stage']}" if c["held_globdb_stage"] else c["r232_status"]
            left_for_later.append((c, reason))
        elif acc_bp + c["bp"] <= max_stage_bp:
            selected.append(c)
            acc_bp += c["bp"]
        else:
            left_for_later.append((c, "spillover_stage_cap"))

    return selected, left_for_later, acc_bp


def census_and_select(
    *,
    iuig_zip: Path,
    umag_zip: Path,
    meta_zip: Path,
    r232_index_dir: Path,
    out_dir: Path,
    gtdb_meta_dir: Path | None = None,
    gems_meta_tar: Path | None = None,
    img_m_tar: Path | None = None,
    non_gtdb_stages_tsv: Path | None = None,
    grare_selection_tsv: Path | None = None,
    max_stage_bp: int = DEFAULT_MAX_STAGE_BP,
    min_unscored_bp: int = 0,
    workers: int = 1,
    cpu_guard_min: float = 280.0,
) -> dict:
    """Stream-census IUIG + uMAG, cross-join against C-74 r232 index, select <= max_stage_bp for JGIUG1."""
    out_dir = Path(out_dir)
    if out_dir.exists():
        raise Refused(f"{out_dir} already exists")

    t0 = time.monotonic()
    refs = _load_r232_crosswalk(
        r232_index_dir=Path(r232_index_dir),
        gtdb_meta_dir=gtdb_meta_dir,
        gems_meta_tar=gems_meta_tar,
        img_m_tar=img_m_tar,
        non_gtdb_stages_tsv=non_gtdb_stages_tsv,
        grare_selection_tsv=grare_selection_tsv,
    )

    # 1. Parse metadata zip (IUIG genome rows + uMAG metagenome ecosystem rows)
    iuig_meta_rows: dict[str, dict] = {}
    non_prok_counts: Counter[str] = Counter()
    non_prok_rows: list[tuple[str, str, str, str]] = []
    umag_eco_by_mg: dict[str, str] = {}

    with zipfile.ZipFile(meta_zip, "r") as zf:
        iuig_name = next((n for n in zf.namelist() if "img_genomes_tarballs_metadata" in n), None)
        umag_name = next((n for n in zf.namelist() if "unrestricted_metagenomes_for_bin_tarballs" in n), None)
        if iuig_name is None or umag_name is None:
            raise Refused(f"{meta_zip}: missing IUIG or uMAG metadata table")

        with zf.open(iuig_name) as f:
            for r in csv.DictReader(io.TextIOWrapper(f, encoding="utf-8", errors="replace"), delimiter="\t"):
                oid = (r.get("taxon_oid") or "").strip()
                dom = (r.get("domain") or "").strip() or "Unclassified"
                if dom not in ("Bacteria", "Archaea"):
                    non_prok_counts[dom] += 1
                    non_prok_rows.append((oid, dom, (r.get("phylum") or "").strip(), (r.get("seq_status") or "").strip()))
                    continue
                iuig_meta_rows[oid] = {
                    "domain": dom,
                    "phylum": (r.get("phylum") or "").strip() or "Unclassified",
                    "ir_class": (r.get("ir_class") or "").strip(),
                    "ir_order": (r.get("ir_order") or "").strip(),
                    "family": (r.get("family") or "").strip(),
                    "genus": (r.get("genus") or "").strip(),
                    "species": (r.get("species") or "").strip(),
                    "seq_status": (r.get("seq_status") or "").strip(),
                    "taxid": (r.get("ncbi_taxon_id") or "").strip(),
                }

        with zf.open(umag_name) as f:
            for r in csv.DictReader(io.TextIOWrapper(f, encoding="utf-8", errors="replace"), delimiter="\t"):
                oid = (r.get("taxon_oid") or "").strip()
                eco_parts = [(r.get(k) or "").strip() for k in ("phylum", "ir_class", "ir_order") if (r.get(k) or "").strip()]
                umag_eco_by_mg[oid] = " / ".join(eco_parts) if eco_parts else "Unclassified"

    # Pre-classify IUIG prokaryote rows so we only measure .fna lengths for novel non-r232 candidates missing bp!
    iuig_pre: dict[str, tuple[str, str, str, str, int]] = {}
    measure_iuig_fna_ids: set[str] = set()
    for oid, mr in iuig_meta_rows.items():
        im = refs["img_meta"].get(oid, {})
        st_info = _classify_r232(
            oid,
            gca=im.get("gca", ""),
            biosample=im.get("bs", ""),
            taxid=mr["taxid"],
            species_str=mr["species"],
            refs=refs,
        )
        iuig_pre[oid] = st_info
        if st_info[4] <= 2:
            measure_iuig_fna_ids.add(oid)

    # 2. Stream-scan IUIG and uMAG archives (parallelized across 2 child processes when workers >= 2)
    if workers >= 2:
        ctx = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=2, mp_context=ctx) as pool:
            fut_iuig = pool.submit(_scan_iuig_archive_worker, (str(iuig_zip), sorted(measure_iuig_fna_ids), cpu_guard_min))
            fut_umag = pool.submit(_scan_umag_archive_worker, (str(umag_zip), cpu_guard_min))
            iuig_scan = fut_iuig.result()
            umag_scan = fut_umag.result()
    else:
        iuig_scan = _scan_iuig_archive_worker((str(iuig_zip), sorted(measure_iuig_fna_ids), cpu_guard_min))
        umag_scan = _scan_umag_archive_worker((str(umag_zip), cpu_guard_min))

    # 3. Build unified census & candidate list across IUIG and uMAG
    all_candidates: list[dict] = []
    stats_by_src = {
        "iuig": {"archive_members": len(iuig_scan["members"]), "by_domain": Counter(), "by_quality": Counter(), "by_r232_status": Counter(), "by_held_globdb": Counter(), "total_est_bp": 0},
        "umag": {"archive_members": len(umag_scan["bins"]), "by_domain": Counter(), "by_quality": Counter(), "by_r232_status": Counter(), "by_held_globdb": Counter(), "total_est_bp": 0},
    }

    # Process IUIG archive members
    iuig_prok_missing_from_archive = 0
    for oid, mr in sorted(iuig_meta_rows.items()):
        if oid not in iuig_scan["members"]:
            iuig_prok_missing_from_archive += 1
            continue
        tgz_bytes, fna_bp, fna_contigs = iuig_scan["members"][oid]
        im = refs["img_meta"].get(oid, {})
        r232_st, gt_acc, gl_id, held_st, nov_tier = iuig_pre[oid]
        bp = fna_bp or im.get("bp", 0) or tgz_bytes
        n_c = fna_contigs or im.get("contigs", 0) or 1

        if oid in refs["gem_checkm"]:
            comp, cont = refs["gem_checkm"][oid]
            if comp < 50.0 or cont > 10.0:
                qtier = "LQ"
                aq = "contig"
            elif comp >= 90.0 and cont <= 5.0:
                qtier = "HQ"
                aq = "mag_hq"
            else:
                qtier = "MQ"
                aq = "mag_mq"
            qscore = round(comp - 5.0 * cont, 4)
            q_scored = True
        else:
            if mr["seq_status"] == "Finished":
                qtier = "HQ"
                aq = "complete" if n_c == 1 else "scaffold"
                comp, cont, qscore = 99.0, 0.0, 99.0
                q_scored = True
            elif min_unscored_bp > 0 and bp < min_unscored_bp:
                qtier = "LQ"
                aq = "contig"
                comp, cont, qscore = 0.0, 0.0, 0.0
                q_scored = False
            else:
                qtier = "MQ"
                aq = "mag_mq" if im.get("uncultured_type") else "contig"
                comp, cont, qscore = 75.0, 1.0, 70.0
                q_scored = False

        lineage = ";".join([
            f"d__{mr['domain']}",
            f"p__{mr['phylum']}",
            f"c__{mr['ir_class']}",
            f"o__{mr['ir_order']}",
            f"f__{mr['family']}",
            f"g__{mr['genus']}",
            f"s__{mr['species']}",
        ])
        eco = im.get("ecosystem") or "Unclassified"
        lgrp_num, lgrp_name = _classify_lineage_group(mr["domain"], mr["phylum"], lineage, eco)

        s_iuig = stats_by_src["iuig"]
        s_iuig["by_domain"][mr["domain"]] += 1
        s_iuig["by_quality"][qtier] += 1
        s_iuig["by_r232_status"][r232_st] += 1
        if held_st:
            s_iuig["by_held_globdb"][held_st] += 1
        s_iuig["total_est_bp"] += bp

        if qtier in ("HQ", "MQ"):
            all_candidates.append({
                "genome_id": oid,
                "source": "iuig",
                "metagenome_id": oid,
                "domain": mr["domain"],
                "phylum": mr["phylum"],
                "lineage": lineage,
                "quality_tier": qtier,
                "assembly_quality": aq,
                "completeness": round(comp, 2),
                "contamination": round(cont, 2),
                "quality_score": qscore,
                "bp": bp,
                "contigs": n_c,
                "tgz_bytes": tgz_bytes,
                "lineage_group_num": lgrp_num,
                "lineage_group": lgrp_name,
                "novelty_tier": nov_tier,
                "r232_status": r232_st,
                "gtdb_r232_acc": gt_acc,
                "globdb_r232_id": gl_id,
                "held_globdb_stage": held_st,
                "ecosystem": eco,
                "_quality_scored": q_scored,
            })

    # Process uMAG archive bins
    for bid, binfo in sorted(umag_scan["bins"].items()):
        tgz_bytes = binfo["tgz_bytes"]
        meta = binfo["meta"] or {}
        mg_id = bid.split("_", 1)[0] if "_" in bid else bid
        eco = umag_eco_by_mg.get(mg_id) or refs["img_meta"].get(mg_id, {}).get("ecosystem") or "Unclassified"

        gt_lin = meta.get("gtdb_lineage") or ""
        bin_lin = meta.get("bin_lineage") or ""
        use_lin = gt_lin if gt_lin and gt_lin != "None" else bin_lin
        dom, phy, sp_tok = _parse_lineage_domain_phylum(use_lin, fallback_domain="Bacteria", fallback_phylum="Unclassified")
        if dom not in ("Bacteria", "Archaea"):
            dom = "Bacteria"

        comp = float(meta.get("completeness") or 0.0)
        cont = float(meta.get("contamination") or 0.0)
        bq = (meta.get("bin_quality") or "").upper()
        if not meta:
            qtier, aq, qscore = "MQ", "mag_mq", 50.0
        elif bq == "LQ" or comp < 50.0 or cont > 10.0:
            qtier, aq, qscore = "LQ", "contig", round(comp - 5.0 * cont, 4)
        elif bq == "HQ" or (comp >= 90.0 and cont <= 5.0):
            qtier, aq, qscore = "HQ", "mag_hq", round(comp - 5.0 * cont, 4)
        else:
            qtier, aq, qscore = "MQ", "mag_mq", round(comp - 5.0 * cont, 4)

        bp = int(meta.get("bp") or 0) or tgz_bytes
        n_c = int(meta.get("contigs") or 0) or 1

        r232_st, gt_acc, gl_id, held_st, nov_tier = _classify_r232(
            bid,
            gca=_stem_gca(gt_lin),
            biosample="",
            taxid="",
            species_str=sp_tok,
            refs=refs,
        )
        lgrp_num, lgrp_name = _classify_lineage_group(dom, phy, use_lin, eco)

        s_umag = stats_by_src["umag"]
        s_umag["by_domain"][dom] += 1
        s_umag["by_quality"][qtier] += 1
        s_umag["by_r232_status"][r232_st] += 1
        if held_st:
            s_umag["by_held_globdb"][held_st] += 1
        s_umag["total_est_bp"] += bp

        if qtier in ("HQ", "MQ"):
            all_candidates.append({
                "genome_id": bid,
                "source": "umag",
                "metagenome_id": mg_id,
                "domain": dom,
                "phylum": phy,
                "lineage": use_lin,
                "quality_tier": qtier,
                "assembly_quality": aq,
                "completeness": round(comp, 2),
                "contamination": round(cont, 2),
                "quality_score": qscore,
                "bp": bp,
                "contigs": n_c,
                "tgz_bytes": tgz_bytes,
                "lineage_group_num": lgrp_num,
                "lineage_group": lgrp_name,
                "novelty_tier": nov_tier,
                "r232_status": r232_st,
                "gtdb_r232_acc": gt_acc,
                "globdb_r232_id": gl_id,
                "held_globdb_stage": held_st,
                "ecosystem": eco,
                "_quality_scored": bool(meta),
            })

    selected, left_for_later, acc_bp = _rank_and_partition_candidates(all_candidates, max_stage_bp=max_stage_bp)

    if not selected:
        raise Refused("census_and_select: 0 genomes selected for JGIUG1")

    tmp = out_dir.with_name(f".{out_dir.name}.tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)

    with open(tmp / "selection.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SELECTION_HEAD, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for c in selected:
            w.writerow(c)

    left_gz_path = tmp / "JGIUG_left_for_later.tsv.gz"
    with open(left_gz_path, "wb") as raw_fh:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_fh, mtime=0) as gz_fh:
            txt_fh = io.TextIOWrapper(gz_fh, encoding="utf-8", newline="")
            w = csv.DictWriter(txt_fh, fieldnames=LEFT_HEAD, delimiter="\t", extrasaction="ignore", lineterminator="\n")
            w.writeheader()
            for c, reason in left_for_later:
                row = dict(c)
                row["defer_reason"] = reason
                w.writerow(row)
            txt_fh.flush()

    make_inputs._write_tsv(
        tmp / "iuig_non_prok_skipped.tsv",
        ["taxon_oid", "domain", "phylum", "seq_status"],
        sorted(non_prok_rows),
    )

    rep = {
        "max_stage_bp": max_stage_bp,
        "iuig": {
            "zip": iuig_scan["zip_info"],
            "metadata_prok_rows": len(iuig_meta_rows),
            "prok_missing_from_archive": iuig_prok_missing_from_archive,
            "archive_members": stats_by_src["iuig"]["archive_members"],
            "non_prokaryote_rows": dict(sorted(non_prok_counts.items())),
            "by_domain": dict(sorted(stats_by_src["iuig"]["by_domain"].items())),
            "by_quality": dict(sorted(stats_by_src["iuig"]["by_quality"].items())),
            "by_r232_status": dict(sorted(stats_by_src["iuig"]["by_r232_status"].items())),
            "by_held_globdb": dict(sorted(stats_by_src["iuig"]["by_held_globdb"].items())),
            "total_est_bp": stats_by_src["iuig"]["total_est_bp"],
            "scan_cpu_min": iuig_scan["cpu_min"],
            "scan_peak_rss_bytes": iuig_scan["peak_rss_bytes"],
        },
        "umag": {
            "zip": umag_scan["zip_info"],
            "metadata_metagenome_rows": len(umag_eco_by_mg),
            "metagenomes_seen": umag_scan["metagenomes_seen"],
            "mbin_tarballs_opened": umag_scan["mbin_tarballs_opened"],
            "archive_members": stats_by_src["umag"]["archive_members"],
            "by_domain": dict(sorted(stats_by_src["umag"]["by_domain"].items())),
            "by_quality": dict(sorted(stats_by_src["umag"]["by_quality"].items())),
            "by_r232_status": dict(sorted(stats_by_src["umag"]["by_r232_status"].items())),
            "by_held_globdb": dict(sorted(stats_by_src["umag"]["by_held_globdb"].items())),
            "total_est_bp": stats_by_src["umag"]["total_est_bp"],
            "scan_cpu_min": umag_scan["cpu_min"],
            "scan_peak_rss_bytes": umag_scan["peak_rss_bytes"],
        },
        "selected": {
            "genomes": len(selected),
            "bp": acc_bp,
            "by_source": dict(sorted(Counter(c["source"] for c in selected).items())),
            "by_domain": dict(sorted(Counter(c["domain"] for c in selected).items())),
            "by_quality": dict(sorted(Counter(c["quality_tier"] for c in selected).items())),
            "by_lineage_group": dict(sorted(Counter(c["lineage_group"] for c in selected).items())),
            "by_novelty_tier": dict(sorted(Counter(str(c["novelty_tier"]) for c in selected).items())),
            "top_phyla": dict(Counter(c["phylum"] for c in selected).most_common(25)),
        },
        "left_for_later": {
            "genomes": len(left_for_later),
            "bp": sum(c["bp"] for c, _ in left_for_later),
            "by_reason": dict(sorted(Counter(reason for _, reason in left_for_later).items())),
            "spillover_novel_genomes": sum(1 for _, r in left_for_later if r == "spillover_stage_cap"),
            "spillover_novel_bp": sum(c["bp"] for c, r in left_for_later if r == "spillover_stage_cap"),
        },
        "sha256": {
            "selection.tsv": ids.sha256_file(tmp / "selection.tsv").hex(),
            "JGIUG_left_for_later.tsv.gz": ids.sha256_file(left_gz_path).hex(),
            "iuig_non_prok_skipped.tsv": ids.sha256_file(tmp / "iuig_non_prok_skipped.tsv").hex(),
        },
        "wall_s": round(time.monotonic() - t0, 2),
        "peak_rss_bytes": _peak_rss_bytes(),
    }
    (tmp / "census_summary.json").write_text(json.dumps(rep, indent=1, sort_keys=True) + "\n")
    tmp.rename(out_dir)
    return rep


def reselect_census(
    census_dir: Path,
    *,
    max_stage_bp: int = DEFAULT_MAX_STAGE_BP,
    min_unscored_bp: int = 100_000,
) -> dict:
    """Re-apply `min_unscored_bp` + strict `not_in_r232` selection in place on an existing census directory."""
    census_dir = Path(census_dir)
    summary_path = census_dir / "census_summary.json"
    sel_path = census_dir / "selection.tsv"
    left_path = census_dir / "JGIUG_left_for_later.tsv.gz"
    skip_path = census_dir / "iuig_non_prok_skipped.tsv"
    for req in (summary_path, sel_path, left_path, skip_path):
        if not req.is_file():
            raise Refused(f"reselect_census: missing {req}")

    rep = json.loads(summary_path.read_text(encoding="utf-8"))
    raw_rows: list[dict[str, str]] = []
    with open(sel_path, encoding="utf-8", newline="") as fh:
        raw_rows.extend(csv.DictReader(fh, delimiter="\t"))
    with gzip.open(left_path, "rt", encoding="utf-8", newline="") as fh:
        raw_rows.extend(csv.DictReader(fh, delimiter="\t"))

    candidates: list[dict] = []
    n_iuig_unscored_lq = 0
    for r in raw_rows:
        bp = int(r["bp"])
        comp = float(r["completeness"])
        cont = float(r["contamination"])
        qscore = float(r["quality_score"])
        is_unscored_iuig = (
            r["source"] == "iuig"
            and r["quality_tier"] == "MQ"
            and comp == 75.0
            and cont == 1.0
            and qscore == 70.0
        )
        if is_unscored_iuig and min_unscored_bp > 0 and bp < min_unscored_bp:
            n_iuig_unscored_lq += 1
            continue
        dom = r["domain"]
        phy = r["phylum"]
        lin = r["lineage"]
        eco = r["ecosystem"]
        lgrp_num, lgrp_name = _classify_lineage_group(dom, phy, lin, eco)
        candidates.append({
            "genome_id": r["genome_id"],
            "source": r["source"],
            "metagenome_id": r["metagenome_id"],
            "domain": dom,
            "phylum": phy,
            "lineage": lin,
            "quality_tier": r["quality_tier"],
            "assembly_quality": r["assembly_quality"],
            "completeness": comp,
            "contamination": cont,
            "quality_score": qscore,
            "bp": bp,
            "contigs": int(r["contigs"]),
            "tgz_bytes": int(r["tgz_bytes"]),
            "lineage_group_num": lgrp_num,
            "lineage_group": lgrp_name,
            "novelty_tier": int(r["novelty_tier"]),
            "r232_status": r["r232_status"],
            "gtdb_r232_acc": r["gtdb_r232_acc"],
            "globdb_r232_id": r["globdb_r232_id"],
            "held_globdb_stage": r["held_globdb_stage"],
            "ecosystem": eco,
            "_quality_scored": not is_unscored_iuig,
        })

    selected, left_for_later, acc_bp = _rank_and_partition_candidates(candidates, max_stage_bp=max_stage_bp)
    if not selected:
        raise Refused("reselect_census: 0 genomes selected")

    if n_iuig_unscored_lq > 0:
        by_q = dict(rep["iuig"]["by_quality"])
        by_q["MQ"] = max(0, int(by_q.get("MQ", 0)) - n_iuig_unscored_lq)
        by_q["LQ"] = int(by_q.get("LQ", 0)) + n_iuig_unscored_lq
        rep["iuig"]["by_quality"] = dict(sorted(by_q.items()))
    rep["min_unscored_bp"] = min_unscored_bp
    rep["max_stage_bp"] = max_stage_bp

    tmp = census_dir.with_name(f".{census_dir.name}.reselect.tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    shutil.copyfile(skip_path, tmp / "iuig_non_prok_skipped.tsv")

    with open(tmp / "selection.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SELECTION_HEAD, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for c in selected:
            w.writerow(c)

    left_gz_path = tmp / "JGIUG_left_for_later.tsv.gz"
    with open(left_gz_path, "wb") as raw_fh:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_fh, mtime=0) as gz_fh:
            txt_fh = io.TextIOWrapper(gz_fh, encoding="utf-8", newline="")
            w = csv.DictWriter(txt_fh, fieldnames=LEFT_HEAD, delimiter="\t", extrasaction="ignore", lineterminator="\n")
            w.writeheader()
            for c, reason in left_for_later:
                row = dict(c)
                row["defer_reason"] = reason
                w.writerow(row)
            txt_fh.flush()

    rep["selected"] = {
        "genomes": len(selected),
        "bp": acc_bp,
        "by_source": dict(sorted(Counter(c["source"] for c in selected).items())),
        "by_domain": dict(sorted(Counter(c["domain"] for c in selected).items())),
        "by_quality": dict(sorted(Counter(c["quality_tier"] for c in selected).items())),
        "by_lineage_group": dict(sorted(Counter(c["lineage_group"] for c in selected).items())),
        "by_novelty_tier": dict(sorted(Counter(str(c["novelty_tier"]) for c in selected).items())),
        "top_phyla": dict(Counter(c["phylum"] for c in selected).most_common(25)),
    }
    rep["left_for_later"] = {
        "genomes": len(left_for_later),
        "bp": sum(c["bp"] for c, _ in left_for_later),
        "by_reason": dict(sorted(Counter(reason for _, reason in left_for_later).items())),
        "spillover_novel_genomes": sum(1 for _, r in left_for_later if r == "spillover_stage_cap"),
        "spillover_novel_bp": sum(c["bp"] for c, r in left_for_later if r == "spillover_stage_cap"),
    }
    rep["sha256"] = {
        "selection.tsv": ids.sha256_file(tmp / "selection.tsv").hex(),
        "JGIUG_left_for_later.tsv.gz": ids.sha256_file(left_gz_path).hex(),
        "iuig_non_prok_skipped.tsv": ids.sha256_file(tmp / "iuig_non_prok_skipped.tsv").hex(),
    }
    (tmp / "census_summary.json").write_text(json.dumps(rep, indent=1, sort_keys=True) + "\n")
    for fname in ("selection.tsv", "JGIUG_left_for_later.tsv.gz", "iuig_non_prok_skipped.tsv", "census_summary.json"):
        os.replace(tmp / fname, census_dir / fname)
    shutil.rmtree(tmp, ignore_errors=True)
    return rep


def _validate_genome_records(gid: str, recs: list[tuple[str, bytes, int]]) -> str | None:
    """Return `None` if every record of `gid` is non-empty clean IUPAC DNA, else a refusal reason for `left_out.tsv`."""
    if not recs:
        return "empty_genome"
    for hdr, seq, _ in recs:
        if "\t" in hdr:
            hdr_clean = hdr.replace("\t", " ")
        else:
            hdr_clean = hdr
        if not hdr_clean.strip():
            return "empty_header"
        if not seq:
            return "empty_record"
        bad = seq.translate(None, IUPAC)
        if bad:
            if any(b in GAPS for b in bad):
                return "gap_byte"
            return f"non_iupac_byte_0x{bad[0]:02x}"
    return None


def _resolve_component_codes(
    streamed_order: list[str],
    unpack_fn,
    guid_map: dict[str, bytes],
    initial_codes: dict[str, int],
) -> tuple[dict[str, int], list[dict]]:
    """Resolve shared-contig `genetic_code` disagreements across full connected components (DEC-27).

    Any two genomes sharing at least one `contig_uid` are unioned into one connected component;
    every genome in a component takes the initial `genetic_code` of the component's canonical leader
    `min(members, key=lambda g: (guid_map[g], g))`.
    """
    parent: dict[str, str] = {g: g for g in streamed_order}

    def _find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def _union(a: str, b: str) -> None:
        ra, rb = _find(a), _find(b)
        if ra != rb:
            parent[rb] = ra

    first_owner: dict[bytes, str] = {}
    sample_shared_contig: dict[str, tuple[bytes, str]] = {}
    for gid in streamed_order:
        for u, _ in unpack_fn(gid):
            prev_g = first_owner.setdefault(u, gid)
            if prev_g != gid:
                _union(prev_g, gid)
                sample_shared_contig.setdefault(gid, (u, prev_g))
                sample_shared_contig.setdefault(prev_g, (u, gid))

    comps: dict[str, list[str]] = {}
    for gid in streamed_order:
        comps.setdefault(_find(gid), []).append(gid)

    code_map: dict[str, int] = {}
    conflicts: list[dict] = []
    for members in sorted(comps.values(), key=lambda ms: min((guid_map[g], g) for g in ms)):
        ordered = sorted(members, key=lambda g: (guid_map[g], g))
        leader = ordered[0]
        target = int(initial_codes[leader])
        for g in ordered:
            orig = int(initial_codes[g])
            if orig != target:
                sh = sample_shared_contig.get(g)
                conflicts.append({
                    "contig_uid": sh[0].hex() if sh else "",
                    "first_genome": leader,
                    "first_code": target,
                    "second_genome": g,
                    "second_code": orig,
                })
            code_map[g] = target
    return code_map, conflicts


def _run_s2_check(stage: str, meta_dir: Path, objects_json: Path) -> dict:
    p = subprocess.run(
        [
            sys.executable,
            str(_HERE / "s2_check_inputs.py"),
            stage,
            str(meta_dir / "genomes.json"),
            str(objects_json),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    doc = json.loads(p.stdout) if p.stdout.strip().startswith("{") else {"stdout": p.stdout}
    return {
        "exit": p.returncode,
        "problems": doc.get("problems"),
        "genomes": doc.get("genomes"),
        "genome_input_rows": doc.get("genome_input_rows"),
        "codes": doc.get("codes"),
        "byte_pass": doc.get("byte_pass"),
        "stderr": p.stderr[-2000:],
    }


def stage_jgiug1(
    *,
    census_dir: Path,
    iuig_zip: Path,
    umag_zip: Path,
    meta_out: Path,
    fasta_out: Path,
    root: Path | None = None,
    root_stages: list[str] | tuple[str, ...] = ROOT_STAGES,
    owner_inputs: dict[str, Path] | None = None,
    later_inputs: dict[str, Path] | None = None,
    inputs_home: Path = Path("/tmp/if2453/s2/inputs"),
    retrieved_at: str = JGIUG_RETRIEVED,
    work: Path | None = None,
    threads: int = 2,
    memory_limit: str = "4GB",
    workers: int = 1,
    max_tasks_per_child: int = 250,
    screen_max_bp: int = 500_000,
    cpu_guard_min: float = 280.0,
    move_large_records_to_stor: bool = False,
) -> dict:
    """Stream selected JGIUG1 genomes, screen genetic codes (11, 25, 4), deduplicate contig_uid, write BGZF & meta."""
    census_dir = Path(census_dir)
    meta_dir = Path(meta_out) / STAGE
    fa_dir = Path(fasta_out) / STAGE
    final_fa = fa_dir / FASTA_NAME
    for p in (meta_dir, final_fa):
        if p.exists():
            raise Refused(f"{p} already exists")

    owner_inputs = dict(owner_inputs or {})
    later_inputs = dict(later_inputs or {})
    for st_name, d in [*owner_inputs.items(), *later_inputs.items()]:
        if not (Path(d) / "genomes.json").is_file():
            raise Refused(f"stage {st_name}: missing genomes.json under {d}")

    with open(census_dir / "selection.tsv", encoding="utf-8", newline="") as fh:
        sel_rows = list(csv.DictReader(fh, delimiter="\t"))
    if not sel_rows:
        raise Refused(f"{census_dir / 'selection.tsv'}: empty selection")

    sel_by_id = {r["genome_id"]: r for r in sel_rows}
    want_iuig = {r["genome_id"] for r in sel_rows if r["source"] == "iuig"}
    want_umag = {r["genome_id"] for r in sel_rows if r["source"] == "umag"}

    tmp = Path(meta_out) / f".{STAGE}.tmp"
    partial = fa_dir / f".{FASTA_NAME}.partial"
    work_dir = Path(work) if work is not None else Path(meta_out) / f".{STAGE}.work"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    fa_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.monotonic()
    pool = (
        ProcessPoolExecutor(
            max_workers=workers,
            mp_context=multiprocessing.get_context("spawn"),
            max_tasks_per_child=max_tasks_per_child,
        )
        if workers > 1
        else None
    )
    pending: set = set()
    gc_results: dict[str, dict] = {}
    streamed_order: list[str] = []
    genome_contigs: dict[str, tuple[bytes, array.array]] = {}
    left_out: list[tuple[str, str, int, int, str]] = []
    seen_headers: set[bytes] = set()
    n_renamed_headers = 0

    def _unpack_contigs(gid: str) -> list[tuple[bytes, int]]:
        raw_b, lens = genome_contigs[gid]
        return [(raw_b[i * 16 : (i + 1) * 16], int(lens[i])) for i in range(len(lens))]

    def _submit_gc(gid: str, seqs: list[bytes], hint: str) -> None:
        nonlocal pending
        if pool is None:
            gc_results[gid] = c58.screen_genetic_code(seqs, hint, screen_max_bp=screen_max_bp)
            return
        pending.add(pool.submit(c58._screen_worker, (gid, seqs, hint, screen_max_bp)))
        if len(pending) >= 4 * workers:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for fut in done:
                g_done, res = fut.result()
                gc_results[g_done] = res

    writer = c77.BgzfWriter(partial)
    try:
        for src_label, zpath, want_ids in (
            ("iuig", Path(iuig_zip), want_iuig),
            ("umag", Path(umag_zip), want_umag),
        ):
            if not want_ids:
                continue
            for gid, recs in iter_selected_from_zip(zpath, want_ids, cpu_guard_min=cpu_guard_min):
                bad_reason = _validate_genome_records(gid, recs)
                if bad_reason is not None:
                    left_out.append((gid, src_label, len(recs), sum(len(s) for _, s, _ in recs), bad_reason))
                    continue

                cu_chunks: list[bytes] = []
                len_arr = array.array("I")
                seqs_for_gc: list[bytes] = []
                for hdr, seq, width in recs:
                    hdr_clean = hdr.replace("\t", " ").strip()
                    first_tok = hdr_clean.split(None, 1)[0]
                    if not first_tok.startswith(gid + "_") and not first_tok.startswith(gid + "|"):
                        hdr_clean = f"{gid}|{hdr_clean}"
                        n_renamed_headers += 1
                    first_tok = hdr_clean.split(None, 1)[0]
                    tok_key = ids.contig_uid(first_tok.encode("utf-8"))
                    if tok_key in seen_headers:
                        hdr_clean = f"{gid}_{len(len_arr)}|{hdr_clean}"
                        first_tok = hdr_clean.split(None, 1)[0]
                        tok_key = ids.contig_uid(first_tok.encode("utf-8"))
                    seen_headers.add(tok_key)
                    writer.write(b"\n".join(fasta.format_record(hdr_clean, seq, width)) + b"\n")
                    cu_chunks.append(ids.contig_uid(seq))
                    len_arr.append(len(seq))
                    seqs_for_gc.append(seq)

                streamed_order.append(gid)
                genome_contigs[gid] = (b"".join(cu_chunks), len_arr)
                r_meta = sel_by_id[gid]
                _submit_gc(gid, seqs_for_gc, f"{r_meta['domain']};{r_meta['phylum']};{r_meta['lineage']}")
                if len(streamed_order) % 500 == 0:
                    print(
                        f"[c87-stage] src={src_label} streamed={len(streamed_order)} gc_done={len(gc_results)} "
                        f"cpu_min={_cpu_min_self():.2f} elapsed={time.monotonic() - t0:.0f}s",
                        file=sys.stderr,
                        flush=True,
                    )

        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for fut in done:
                g_done, res = fut.result()
                gc_results[g_done] = res
    except BaseException:
        if pool is not None:
            pool.shutdown(wait=False, cancel_futures=True)
            pool = None
        writer.close()
        partial.unlink(missing_ok=True)
        shutil.rmtree(tmp, ignore_errors=True)
        if work is None:
            shutil.rmtree(work_dir, ignore_errors=True)
        raise
    finally:
        writer.close()
        if pool is not None:
            pool.shutdown(wait=True)

    if not streamed_order:
        partial.unlink(missing_ok=True)
        shutil.rmtree(tmp, ignore_errors=True)
        raise Refused("stage_jgiug1: no valid genome survived IUPAC validation")

    # Copy JGIUG_left_for_later.tsv.gz into stage meta dir so stage mirror carries it
    left_later_src = census_dir / "JGIUG_left_for_later.tsv.gz"
    if left_later_src.is_file():
        shutil.copyfile(left_later_src, tmp / "JGIUG_left_for_later.tsv.gz")

    # Resolve shared-contig genetic_code conflicts across full connected components by least genome_uid (DEC-27)
    guid_map = {
        gid: ids.genome_uid(list(dict.fromkeys(u for u, _ in _unpack_contigs(gid))))
        for gid in streamed_order
    }
    initial_code_map = {gid: int(gc_results[gid]["genetic_code"]) for gid in streamed_order}
    code_map, code_conflicts = _resolve_component_codes(
        streamed_order, _unpack_contigs, guid_map, initial_code_map
    )

    # Write records.tsv, genetic_codes.tsv, r232_status.tsv, left_out.tsv
    records_path = tmp / "records.tsv"
    gc_path = tmp / "genetic_codes.tsv"
    status_path = tmp / "r232_status.tsv"
    ordinal = 0

    def _fmt(v) -> str:
        return "" if v is None else str(v)

    with (
        open(records_path, "w", encoding="utf-8", buffering=1 << 20) as rf,
        open(gc_path, "w", encoding="utf-8", buffering=1 << 20) as gf,
        open(status_path, "w", encoding="utf-8", buffering=1 << 20) as sf,
    ):
        rf.write("\t".join(RECORDS_HEAD) + "\n")
        gf.write("\t".join(GC_HEAD) + "\n")
        sf.write("\t".join(STATUS_HEAD) + "\n")
        for gid in streamed_order:
            rm = sel_by_id[gid]
            s = gc_results[gid]
            c_list = _unpack_contigs(gid)
            bp_actual = sum(n for _, n in c_list)
            n_c_actual = len(c_list)
            code_val = code_map[gid]

            gf.write("\t".join([
                gid, rm["source"], rm["domain"], rm["phylum"], rm["lineage"],
                str(bp_actual), str(n_c_actual), rm["assembly_quality"],
                str(code_val), s["reason"], "true" if s["tested_alt"] else "false",
                _fmt(s["cd_11"]), _fmt(s["mean_orf_11"]), _fmt(s["n_genes_11"]),
                _fmt(s["cd_25"]), _fmt(s["mean_orf_25"]), _fmt(s["n_genes_25"]),
                _fmt(s["cd_4"]), _fmt(s["mean_orf_4"]), _fmt(s["n_genes_4"]),
                _fmt(s["tga_per_kb"]),
            ]) + "\n")

            sf.write("\t".join([
                gid, rm["source"], rm["domain"], rm["phylum"], rm["quality_tier"],
                str(bp_actual), str(n_c_actual),
                rm["r232_status"], rm["gtdb_r232_acc"], rm["globdb_r232_id"],
                rm["held_globdb_stage"], rm["lineage"],
            ]) + "\n")

            for rec_ord, (u, length) in enumerate(c_list):
                rf.write("\t".join([
                    str(ordinal), gid, rm["source"], str(rec_ord), u.hex(), str(length),
                    rm["domain"], rm["phylum"], rm["quality_tier"], rm["assembly_quality"],
                    str(code_val), rm["r232_status"],
                ]) + "\n")
                ordinal += 1

    make_inputs._write_tsv(tmp / "left_out.tsv", LEFT_OUT_HEAD, left_out)

    # Exact contig_uid ownership against root + r26.11.4 (owner_inputs) and r26.11.5 (later_inputs)
    owners = c58._owners_from_records(
        records_path,
        None if root is None else Path(root),
        list(root_stages),
        owner_inputs,
        work=work_dir,
        threads=threads,
        memory_limit=memory_limit,
    )
    later = (
        c58._owners_from_records(
            records_path,
            None,
            [],
            later_inputs,
            work=work_dir,
            threads=threads,
            memory_limit=memory_limit,
        )
        if later_inputs
        else {}
    )

    listing = make_inputs.Listing()
    genomes: list[dict] = []
    repeated: list[tuple] = []
    shared: list[tuple] = []
    held_later: list[tuple] = []
    first_of_guid: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}
    ord_cur = 0

    for gid in streamed_order:
        rm = sel_by_id[gid]
        c_list = _unpack_contigs(gid)
        raw_uids = [u for u, _ in c_list]
        for rec_ord, (u, length) in enumerate(c_list):
            contig_len.setdefault(u, length)
            if u in owners:
                shared.append((ord_cur, gid, rec_ord, u.hex(), length, owners[u]))
            elif u in later:
                held_later.append((ord_cur, gid, rec_ord, u.hex(), length, later[u]))
            ord_cur += 1

        entry = listing.add(raw_uids, code_map[gid], "prokaryote", rm["assembly_quality"])
        g_uid = guid_map[gid]
        if entry is None:
            repeated.append((ord_cur - len(c_list), gid, g_uid.hex(), first_of_guid[g_uid]))
            continue
        genomes.append(entry)
        first_of_guid[g_uid] = gid

    make_inputs._write_tsv(tmp / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
    make_inputs._write_tsv(
        tmp / "shared_records.tsv",
        ["ordinal", "genome", "record_ordinal", "contig_uid", "length", "owner_stage"],
        shared,
    )
    make_inputs._write_tsv(
        tmp / "later_records.tsv",
        ["ordinal", "genome", "record_ordinal", "contig_uid", "length", "holder"],
        held_later,
    )

    src_obj = {**JGIUG_SOURCE, "retrieved_at": retrieved_at}
    (tmp / "objects.json").write_text(
        json.dumps([{"path": str(Path(inputs_home) / STAGE / FASTA_NAME), **src_obj}], indent=1) + "\n"
    )
    (tmp / "objects_staging.json").write_text(
        json.dumps([{"path": str(final_fa), **src_obj}], indent=1) + "\n"
    )
    (tmp / "genomes.json").write_text(json.dumps(genomes))

    if str(fa_dir).startswith("/stor/") and partial.stat().st_size < STOR_FLOOR_BYTES:
        partial.unlink(missing_ok=True)
        shutil.rmtree(tmp, ignore_errors=True)
        raise Refused(
            f"{partial}: {partial.stat().st_size} B is under the 64 MiB /stor floor ({STOR_FLOOR_BYTES} B; AGENTS.md §7.3)"
        )

    partial.rename(final_fa)

    records_final_dir = meta_dir
    if move_large_records_to_stor and records_path.stat().st_size >= STOR_FLOOR_BYTES:
        dst_rec = fa_dir / "records.tsv"
        tmp_rec = fa_dir / ".records.tsv.partial"
        shutil.copyfile(records_path, tmp_rec)
        tmp_rec.rename(dst_rec)
        records_path.unlink()
        records_final_dir = fa_dir

    tmp.rename(meta_dir)
    if work is None:
        shutil.rmtree(work_dir, ignore_errors=True)

    s2_res = _run_s2_check(STAGE, meta_dir, meta_dir / "objects_staging.json")

    own_c: Counter[str] = Counter()
    own_b: Counter[str] = Counter()
    lat_c: Counter[str] = Counter()
    lat_b: Counter[str] = Counter()
    for cu, length in contig_len.items():
        if cu in owners:
            own_c[owners[cu]] += 1
            own_b[owners[cu]] += length
        elif cu in later:
            lat_c[later[cu]] += 1
            lat_b[later[cu]] += length

    tracked_files: dict[str, Path] = {
        "genomes.json": meta_dir / "genomes.json",
        "objects.json": meta_dir / "objects.json",
        "objects_staging.json": meta_dir / "objects_staging.json",
        "records.tsv": records_final_dir / "records.tsv",
        "repeated_records.tsv": meta_dir / "repeated_records.tsv",
        "shared_records.tsv": meta_dir / "shared_records.tsv",
        "later_records.tsv": meta_dir / "later_records.tsv",
        "genetic_codes.tsv": meta_dir / "genetic_codes.tsv",
        "r232_status.tsv": meta_dir / "r232_status.tsv",
        "left_out.tsv": meta_dir / "left_out.tsv",
        FASTA_NAME: final_fa,
    }
    if (meta_dir / "JGIUG_left_for_later.tsv.gz").is_file():
        tracked_files["JGIUG_left_for_later.tsv.gz"] = meta_dir / "JGIUG_left_for_later.tsv.gz"

    missing_selected = sorted(set(sel_by_id) - set(streamed_order) - {x[0] for x in left_out})

    rep = {
        "stage": STAGE,
        "selected_genomes": len(sel_rows),
        "missing_from_archives": missing_selected,
        "non_iupac_dropped_genomes": len(left_out),
        "streamed_valid_genomes": len(streamed_order),
        "genomes": len(genomes),
        "genomes_repeated": len(repeated),
        "records_total": ord_cur,
        "renamed_headers": n_renamed_headers,
        "contigs_distinct": len(contig_len),
        "bases_distinct": sum(contig_len.values()),
        "held_by": dict(sorted(own_c.items())),
        "held_bases_by": dict(sorted(own_b.items())),
        "later_held_by": dict(sorted(lat_c.items())),
        "later_held_bases_by": dict(sorted(lat_b.items())),
        "contigs_new": sum(1 for cu in contig_len if cu not in owners),
        "bases_new": sum(n for cu, n in contig_len.items() if cu not in owners),
        "contigs_new_after_later": sum(1 for cu in contig_len if cu not in owners and cu not in later),
        "bases_new_after_later": sum(n for cu, n in contig_len.items() if cu not in owners and cu not in later),
        "codes": dict(sorted(Counter(str(g["genetic_code"]) for g in genomes).items())),
        "tested_alt_genomes": sum(1 for g in streamed_order if gc_results[g]["tested_alt"]),
        "code_conflicts": code_conflicts,
        "assembly_quality": dict(sorted(Counter(g["assembly_quality"] for g in genomes).items())),
        "by_source": dict(sorted(Counter(sel_by_id[g]["source"] for g in streamed_order).items())),
        "by_domain": dict(sorted(Counter(sel_by_id[g]["domain"] for g in streamed_order).items())),
        "r232_status": dict(sorted(Counter(sel_by_id[g]["r232_status"] for g in streamed_order).items())),
        "records_location": str(records_final_dir / "records.tsv"),
        "sha256": {k: ids.sha256_file(p).hex() for k, p in tracked_files.items()},
        "bytes": {k: p.stat().st_size for k, p in tracked_files.items()},
        "s2_check": s2_res,
        "cpu_min_self": round(_cpu_min_self(), 3),
        "wall_s": round(time.monotonic() - t0, 2),
        "peak_rss_bytes": _peak_rss_bytes(),
        **{"listing_" + k: v for k, v in listing.report().items()},
    }
    (meta_dir / "stage_report.json").write_text(json.dumps(rep, indent=1, sort_keys=True) + "\n")
    return rep


def reconcile_stage_codes(
    meta_dir: Path,
    records_in_path: Path,
) -> dict:
    """Reconcile shared-contig `genetic_code` components in `meta_dir` without touching `/stor` FASTA or records."""
    meta_dir = Path(meta_dir)
    records_in_path = Path(records_in_path)
    gc_path = meta_dir / "genetic_codes.tsv"
    rep_path = meta_dir / "stage_report.json"
    objs_stg = meta_dir / "objects_staging.json"
    for req in (gc_path, rep_path, objs_stg, records_in_path):
        if not req.is_file():
            raise Refused(f"reconcile_stage_codes: missing {req}")

    rep = json.loads(rep_path.read_text(encoding="utf-8"))
    gc_rows: list[dict[str, str]] = []
    initial_codes: dict[str, int] = {}
    streamed_order: list[str] = []
    with open(gc_path, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gid = r["genome"]
            streamed_order.append(gid)
            gc_rows.append(r)
            rsn = r["reason"]
            if rsn == "uga_recoded_code25":
                initial_codes[gid] = 25
            elif rsn == "uga_recoded_code4":
                initial_codes[gid] = 4
            else:
                initial_codes[gid] = 11

    genome_contigs: dict[str, tuple[bytes, array.array]] = {}
    rec_meta: dict[str, tuple[str, str, str, str, str, str]] = {}
    cur_gid: str | None = None
    cu_chunks: list[bytes] = []
    len_arr = array.array("I")
    with open(records_in_path, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gid = r["genome"]
            if gid != cur_gid:
                if cur_gid is not None:
                    genome_contigs[cur_gid] = (b"".join(cu_chunks), len_arr)
                cur_gid = gid
                cu_chunks = []
                len_arr = array.array("I")
                rec_meta[gid] = (
                    r["source"], r["domain"], r["phylum"],
                    r["quality_tier"], r["assembly_quality"], r["r232_status"],
                )
            cu_chunks.append(bytes.fromhex(r["contig_uid"]))
            len_arr.append(int(r["length"]))
        if cur_gid is not None:
            genome_contigs[cur_gid] = (b"".join(cu_chunks), len_arr)

    def _unpack(gid: str) -> list[tuple[bytes, int]]:
        raw_b, lens = genome_contigs[gid]
        return [(raw_b[i * 16 : (i + 1) * 16], int(lens[i])) for i in range(len(lens))]

    guid_map = {
        gid: ids.genome_uid(list(dict.fromkeys(u for u, _ in _unpack(gid))))
        for gid in streamed_order
    }
    code_map, code_conflicts = _resolve_component_codes(streamed_order, _unpack, guid_map, initial_codes)

    # Atomically write updated genetic_codes.tsv, records.tsv (inside meta_dir), and genomes.json
    gc_tmp = meta_dir / ".genetic_codes.tsv.tmp"
    with open(gc_tmp, "w", encoding="utf-8", newline="", buffering=1 << 20) as gf:
        w = csv.DictWriter(gf, fieldnames=GC_HEAD, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in gc_rows:
            row = dict(r)
            row["genetic_code"] = str(code_map[r["genome"]])
            w.writerow(row)

    rec_out = meta_dir / "records.tsv"
    rec_tmp = meta_dir / ".records.tsv.tmp"
    ordinal = 0
    with open(rec_tmp, "w", encoding="utf-8", buffering=1 << 20) as rf:
        rf.write("\t".join(RECORDS_HEAD) + "\n")
        for gid in streamed_order:
            src_s, dom_s, phy_s, qt_s, aq_s, r232_s = rec_meta[gid]
            c_str = str(code_map[gid])
            for rec_ord, (u, length) in enumerate(_unpack(gid)):
                rf.write("\t".join([
                    str(ordinal), gid, src_s, str(rec_ord), u.hex(), str(length),
                    dom_s, phy_s, qt_s, aq_s, c_str, r232_s,
                ]) + "\n")
                ordinal += 1

    listing = make_inputs.Listing()
    genomes: list[dict] = []
    for gid in streamed_order:
        _, _, _, _, aq_s, _ = rec_meta[gid]
        raw_uids = [u for u, _ in _unpack(gid)]
        entry = listing.add(raw_uids, code_map[gid], "prokaryote", aq_s)
        if entry is not None:
            genomes.append(entry)

    gen_tmp = meta_dir / ".genomes.json.tmp"
    gen_tmp.write_text(json.dumps(genomes))

    os.replace(gc_tmp, gc_path)
    os.replace(rec_tmp, rec_out)
    os.replace(gen_tmp, meta_dir / "genomes.json")

    s2_res = _run_s2_check(STAGE, meta_dir, objs_stg)
    final_fa = Path(json.loads(objs_stg.read_text(encoding="utf-8"))[0]["path"])

    for k, p in (
        ("genomes.json", meta_dir / "genomes.json"),
        ("genetic_codes.tsv", gc_path),
        ("records.tsv", rec_out),
    ):
        rep["sha256"][k] = ids.sha256_file(p).hex()
        rep["bytes"][k] = p.stat().st_size
    if FASTA_NAME not in rep["sha256"] and final_fa.is_file():
        rep["sha256"][FASTA_NAME] = ids.sha256_file(final_fa).hex()
        rep["bytes"][FASTA_NAME] = final_fa.stat().st_size

    rep["genomes"] = len(genomes)
    rep["codes"] = dict(sorted(Counter(str(g["genetic_code"]) for g in genomes).items()))
    rep["code_conflicts"] = code_conflicts
    rep["records_location"] = str(rec_out)
    rep["s2_check"] = s2_res
    (meta_dir / "stage_report.json").write_text(json.dumps(rep, indent=1, sort_keys=True) + "\n")
    return rep


def write_conf_template(
    out_conf: Path,
    *,
    stage: str = STAGE,
    note: str,
    after: str = "GBRMGD",
    sums: dict[str, str],
    expect_held_by: dict[str, int],
    summary_comment: str = "",
) -> Path:
    """Write `stages/templates/<stage>.conf` and return its path."""
    out_conf = Path(out_conf)
    out_conf.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"# Stage {stage}: pod-resident JGI uncultivated prokaryote archives (IUIG + IMG-UMAG prokaryote bins; card C-87, DEC-80, DEC-88).",
    ]
    if summary_comment:
        for ln in summary_comment.splitlines():
            lines.append(f"# {ln}")
    lines += [
        "# PRE-BUILD TEMPLATE (promoted to stages/ when the r26.11.6 chain order is finalized).",
        "# Read by stage_driver.sh through `stage_est.py conf`; every value is from the note below unless said.",
        f"STAGE={stage}",
        f"NOTE={note}",
        f"INPUTS=/tmp/if2453/s2/inputs/{stage}",
    ]
    for name in ("genomes.json", "objects.json", FASTA_NAME):
        if name not in sums:
            raise Refused(f"write_conf_template: missing sha256 for {name}")
        lines.append(f"SUM={sums[name]} {name}")
    lines += [
        "DEPTH=called",
        "CRISPR=yes",
        f"AFTER={after}",
    ]
    for owner in sorted(expect_held_by):
        cnt = int(expect_held_by[owner])
        if cnt > 0:
            lines.append(f"EXPECT=held_by {owner} {cnt} {note}")
    lines.append("# No least contig length (DEC-60 (1): MAGs and uncultivated genomes stay whole). No JUDGE or ACCEPTED line.")
    tmp = out_conf.with_name(f".{out_conf.name}.tmp.{os.getpid()}")
    tmp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    os.replace(tmp, out_conf)
    return out_conf


def _parse_kv_pairs(items: list[str], flag: str) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for item in items:
        k, eq, v = item.partition("=")
        if not eq or not k or not v:
            raise Refused(f"{flag} {item!r}: expected STAGE=DIR")
        out[k] = Path(v)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Card C-87 JGI IUIG + uMAG census and staging")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_cen = sub.add_parser("census")
    p_cen.add_argument("--iuig-zip", type=Path, required=True)
    p_cen.add_argument("--umag-zip", type=Path, required=True)
    p_cen.add_argument("--meta-zip", type=Path, required=True)
    p_cen.add_argument("--r232-index", type=Path, required=True)
    p_cen.add_argument("--out-dir", type=Path, required=True)
    p_cen.add_argument("--gtdb-meta-dir", type=Path, default=None)
    p_cen.add_argument("--gems-meta-tar", type=Path, default=None)
    p_cen.add_argument("--img-m-tar", type=Path, default=None)
    p_cen.add_argument("--non-gtdb-stages", type=Path, default=None)
    p_cen.add_argument("--grare-selection", type=Path, default=None)
    p_cen.add_argument("--max-stage-bp", type=int, default=DEFAULT_MAX_STAGE_BP)
    p_cen.add_argument("--min-unscored-bp", type=int, default=100_000)
    p_cen.add_argument("--workers", type=int, default=2)
    p_cen.add_argument("--cpu-guard-min", type=float, default=280.0)
    p_cen.add_argument("--report-out", type=Path, default=None)

    p_res = sub.add_parser("reselect")
    p_res.add_argument("--census-dir", type=Path, required=True)
    p_res.add_argument("--max-stage-bp", type=int, default=DEFAULT_MAX_STAGE_BP)
    p_res.add_argument("--min-unscored-bp", type=int, default=100_000)
    p_res.add_argument("--report-out", type=Path, default=None)

    p_stg = sub.add_parser("stage")
    p_stg.add_argument("--census-dir", type=Path, required=True)
    p_stg.add_argument("--iuig-zip", type=Path, required=True)
    p_stg.add_argument("--umag-zip", type=Path, required=True)
    p_stg.add_argument("--meta-out", type=Path, required=True)
    p_stg.add_argument("--fasta-out", type=Path, required=True)
    p_stg.add_argument("--root", type=Path, default=None)
    p_stg.add_argument("--root-stages", default=",".join(ROOT_STAGES))
    p_stg.add_argument("--chain-home", type=Path, default=None)
    p_stg.add_argument("--use-a29-later", action="store_true", help="include all 22 A-29 r26.11.5 stages from STAGING_DIRS")
    p_stg.add_argument("--later", action="append", default=[])
    p_stg.add_argument("--inputs-home", type=Path, default=Path("/tmp/if2453/s2/inputs"))
    p_stg.add_argument("--retrieved-at", default=JGIUG_RETRIEVED)
    p_stg.add_argument("--work", type=Path, default=None)
    p_stg.add_argument("--threads", type=int, default=2)
    p_stg.add_argument("--memory-limit", default="4GB")
    p_stg.add_argument("--workers", type=int, default=4)
    p_stg.add_argument("--max-tasks-per-child", type=int, default=250)
    p_stg.add_argument("--screen-max-bp", type=int, default=500_000)
    p_stg.add_argument("--cpu-guard-min", type=float, default=280.0)
    p_stg.add_argument("--move-large-records-to-stor", action="store_true")
    p_stg.add_argument("--report-out", type=Path, default=None)

    p_cnf = sub.add_parser("write-conf")
    p_cnf.add_argument("--stage-report", type=Path, required=True)
    p_cnf.add_argument("--out-conf", type=Path, required=True)
    p_cnf.add_argument("--note", default="planning/measurements/2026-11-27_c87_jgi_ug_inputs.md")
    p_cnf.add_argument("--after", default="GBRMGD")
    p_cnf.add_argument("--include-later-expect", action="store_true")

    p_rec = sub.add_parser("reconcile-codes")
    p_rec.add_argument("--meta-dir", type=Path, required=True)
    p_rec.add_argument("--records-in", type=Path, required=True)
    p_rec.add_argument("--report-out", type=Path, default=None)

    args = ap.parse_args(argv)
    try:
        if args.cmd == "census":
            rep = census_and_select(
                iuig_zip=args.iuig_zip,
                umag_zip=args.umag_zip,
                meta_zip=args.meta_zip,
                r232_index_dir=args.r232_index,
                out_dir=args.out_dir,
                gtdb_meta_dir=args.gtdb_meta_dir,
                gems_meta_tar=args.gems_meta_tar,
                img_m_tar=args.img_m_tar,
                non_gtdb_stages_tsv=args.non_gtdb_stages,
                grare_selection_tsv=args.grare_selection,
                max_stage_bp=args.max_stage_bp,
                min_unscored_bp=args.min_unscored_bp,
                workers=args.workers,
                cpu_guard_min=args.cpu_guard_min,
            )
            txt = json.dumps(rep, indent=1, sort_keys=True)
            if args.report_out is not None:
                args.report_out.parent.mkdir(parents=True, exist_ok=True)
                args.report_out.write_text(txt + "\n")
            print(txt)
            return 0

        if args.cmd == "reselect":
            rep = reselect_census(
                args.census_dir,
                max_stage_bp=args.max_stage_bp,
                min_unscored_bp=args.min_unscored_bp,
            )
            txt = json.dumps(rep, indent=1, sort_keys=True)
            if args.report_out is not None:
                args.report_out.parent.mkdir(parents=True, exist_ok=True)
                args.report_out.write_text(txt + "\n")
            print(txt)
            return 0

        if args.cmd == "stage":
            owners = {st: args.chain_home / st for st, _ in a26.R26_11_4_CHAIN} if args.chain_home is not None else {}
            later: dict[str, Path] = {}
            if args.use_a29_later:
                for st_name in a29.R26_11_5_CHAIN:
                    later[st_name] = Path(a29.STAGING_DIRS[st_name])
            later.update(_parse_kv_pairs(args.later, "--later"))
            roots = [s for s in args.root_stages.split(",") if s]
            rep = stage_jgiug1(
                census_dir=args.census_dir,
                iuig_zip=args.iuig_zip,
                umag_zip=args.umag_zip,
                meta_out=args.meta_out,
                fasta_out=args.fasta_out,
                root=args.root,
                root_stages=roots,
                owner_inputs=owners,
                later_inputs=later,
                inputs_home=args.inputs_home,
                retrieved_at=args.retrieved_at,
                work=args.work,
                threads=args.threads,
                memory_limit=args.memory_limit,
                workers=args.workers,
                max_tasks_per_child=args.max_tasks_per_child,
                screen_max_bp=args.screen_max_bp,
                cpu_guard_min=args.cpu_guard_min,
                move_large_records_to_stor=args.move_large_records_to_stor,
            )
            txt = json.dumps(rep, indent=1, sort_keys=True)
            if args.report_out is not None:
                args.report_out.parent.mkdir(parents=True, exist_ok=True)
                args.report_out.write_text(txt + "\n")
            print(txt)
            return 0 if rep["s2_check"]["exit"] == 0 and rep["s2_check"]["problems"] == [] else 1

        if args.cmd == "write-conf":
            rep = json.loads(args.stage_report.read_text(encoding="utf-8"))
            expect = dict(rep.get("held_by", {}))
            if args.include_later_expect:
                for k, v in rep.get("later_held_by", {}).items():
                    expect[k] = expect.get(k, 0) + v
            summary = (
                f"{rep['genomes']} distinct genomes ({rep['contigs_distinct']} distinct contigs, "
                f"{rep['bases_distinct']} bp; {rep['contigs_new_after_later']} contigs / "
                f"{rep['bases_new_after_later']} bp net new after S1..GBRMGD); redistributable: false."
            )
            write_conf_template(
                args.out_conf,
                stage=STAGE,
                note=args.note,
                after=args.after,
                sums=rep["sha256"],
                expect_held_by=expect,
                summary_comment=summary,
            )
            print(str(args.out_conf))
            return 0

        if args.cmd == "reconcile-codes":
            rep = reconcile_stage_codes(args.meta_dir, args.records_in)
            txt = json.dumps(rep, indent=1, sort_keys=True)
            if args.report_out is not None:
                args.report_out.parent.mkdir(parents=True, exist_ok=True)
                args.report_out.write_text(txt + "\n")
            print(txt)
            return 0 if rep["s2_check"]["exit"] == 0 and rep["s2_check"]["problems"] == [] else 1

        return 2
    except Refused as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
