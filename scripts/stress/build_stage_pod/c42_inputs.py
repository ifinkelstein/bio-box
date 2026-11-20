"""Card C-42 step (4): direct stage input builder (`genomes.json`, `objects.json`) and depth-1 cost analyzer.

    python -I c42_inputs.py plan-stages --selection S --sample-stats T --fasta-dir D --fetch-log L --out-dir O
        [--min-contig-len 0] [--max-stage-bp 25000000000] [--max-stage-contigs 12000000]
    python -I c42_inputs.py summarize-depth1 --scratch-root R --stages-manifest M --out-dir O

Writes `genomes.json` and `objects.json` directly (without modifying `make_inputs.py`):
- One `GenomeInput` per distinct `contig_uid` in the stage (`kind="prokaryote"`, `genetic_code=11`,
  `circular=[]`, `assembly_quality="contig"`), sorted deterministically by `contig_uid` hex so output
  bytes are independent of input sample order (AGENTS §4.3).
- Supports `--min-contig-len L` (default 0): when `L > 0`, writes deterministic gzip (`mtime=0`)
  filtered per-sample FASTAs and drops `< L` contigs consistently from both `objects.json` and
  `genomes.json`, recording dropped records and bases.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import heapq
import json
import shutil
import sys
from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path

import duckdb
import pyarrow.parquet as pq

from biobox import fasta, ids, io

SPIRE_SOURCE = {
    "source_name": "spire",
    "source_release": "v1",
    "license": "CC BY-SA 4.0",
    "redistributable": True,
}
DEFAULT_RETRIEVED = "2026-11-18T00:00:00Z"


class Refused(Exception):
    """Invalid stage inputs (refused before writing)."""


def stage_root(base: Path, stage: str) -> Path:
    """Resolve either a per-stage root `<base>/<stage>` or a shared root `<base>`."""
    base = Path(base)
    cand = base / stage
    return cand if (cand / "l0seq").exists() or (cand / "runs").exists() else base


def _iter_stage_owned_uids(root: Path, stage: str) -> Iterator[bytes]:
    sr = stage_root(root, stage)
    files = io.part_files(sr / "l0seq" / "contig_index" / f"stage={stage}")
    for f in files:
        pf = pq.ParquetFile(f)
        for b in pf.iter_batches(batch_size=65_536, columns=["contig_uid", "is_default_variant"]):
            for cu, is_def in zip(b.column(0).to_pylist(), b.column(1).to_pylist()):
                if is_def:
                    yield cu


def write_owned_uids_bin(roots_base: Path, stages: list[str], dst_bin: Path) -> int:
    """Stream sorted 16-byte owned `contig_uid`s across `stages` into a compact binary file (`16 * N` bytes)."""
    dst_bin = Path(dst_bin)
    dst_bin.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst_bin.with_name(f".{dst_bin.name}.tmp")
    n = 0
    last = None
    buf = bytearray()
    with open(tmp, "wb") as fh:
        for cu in heapq.merge(*(_iter_stage_owned_uids(roots_base, s) for s in stages)):
            if cu == last:
                continue
            last = cu
            buf.extend(cu)
            n += 1
            if len(buf) >= (1 << 20):
                fh.write(buf)
                buf.clear()
        if buf:
            fh.write(buf)
            buf.clear()
    tmp.replace(dst_bin)
    return n


class _Bytes16View:
    """Zero-copy sequence view over sorted 16-byte binary UIDs (`bisect` probes O(log N) slices)."""

    __slots__ = ("_buf", "_n")

    def __init__(self, buf: bytes):
        self._buf = buf
        self._n = len(buf) // 16

    def __len__(self) -> int:
        return self._n

    def __getitem__(self, i: int) -> bytes:
        o = i * 16
        return self._buf[o : o + 16]


def _in_sorted_uids(view: _Bytes16View | None, cu: bytes) -> bool:
    if view is None or len(view) == 0:
        return False
    idx = bisect.bisect_left(view, cu)
    return idx < len(view) and view[idx] == cu


def _write_deterministic_fa(records: list[tuple[str, bytes, int]], dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(f".{dst.name}.tmp")
    with open(tmp, "wb") as out_fh:
        for hdr, seq, width in records:
            for line in fasta.format_record(hdr, seq, width):
                out_fh.write(line + b"\n")
    tmp.replace(dst)


def build_stage_inputs(
    samples: list[dict],
    fasta_by_sample: dict[str, Path],
    out_dir: Path,
    *,
    min_contig_len: int = 0,
    exclude_uids_bin: Path | None = None,
) -> dict:
    """Write `genomes.json` and `objects.json` for a stage of SPIRE samples into `out_dir`."""
    out_dir = Path(out_dir)
    if out_dir.exists():
        raise Refused(f"{out_dir} already exists")
    if min_contig_len < 0:
        raise Refused(f"min_contig_len must be >= 0, got {min_contig_len}")

    excl_view: _Bytes16View | None = None
    if exclude_uids_bin is not None and Path(exclude_uids_bin).is_file() and Path(exclude_uids_bin).stat().st_size > 0:
        excl_view = _Bytes16View(Path(exclude_uids_bin).read_bytes())

    seen_sids: set[str] = set()
    ordered = sorted(samples, key=lambda r: r["sample_id"])
    for s in ordered:
        sid = s["sample_id"]
        if sid in seen_sids:
            raise Refused(f"duplicate sample_id {sid!r}")
        seen_sids.add(sid)
        fa = fasta_by_sample.get(sid)
        if fa is None or not Path(fa).is_file():
            raise Refused(f"missing FASTA for sample {sid!r}: {fa}")

    work = out_dir.with_name(f".{out_dir.name}.tmp")
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    contig_lens: dict[bytes, int] = {}
    objects: list[dict] = []
    cross_stage_rows: list[tuple[str, str, str, int]] = []
    total_records = kept_records = dropped_records = held_records = 0
    total_bp = kept_bp = dropped_bp = held_bp = 0
    empty_samples = 0

    for s in ordered:
        sid = s["sample_id"]
        src_fa = Path(fasta_by_sample[sid]).resolve()
        retrieved = s.get("retrieved_at") or DEFAULT_RETRIEVED

        if min_contig_len == 0 and excl_view is None:
            s_recs = 0
            for rec in fasta.iter_records(src_fa):
                ln = len(rec.seq)
                total_records += 1
                kept_records += 1
                total_bp += ln
                kept_bp += ln
                s_recs += 1
                cu = ids.contig_uid(rec.seq)
                contig_lens.setdefault(cu, ln)
            if s_recs == 0:
                empty_samples += 1
                continue
            objects.append({"path": str(src_fa), **SPIRE_SOURCE, "retrieved_at": retrieved})
        else:
            kept_for_sample: list[tuple[str, bytes, int]] = []
            for rec in fasta.iter_records(src_fa):
                ln = len(rec.seq)
                total_records += 1
                total_bp += ln
                if ln < min_contig_len:
                    dropped_records += 1
                    dropped_bp += ln
                    continue
                cu = ids.contig_uid(rec.seq)
                if _in_sorted_uids(excl_view, cu):
                    held_records += 1
                    held_bp += ln
                    acc = rec.header.split(None, 1)[0] if rec.header else ""
                    cross_stage_rows.append((sid, cu.hex(), acc, ln))
                    continue
                kept_records += 1
                kept_bp += ln
                contig_lens.setdefault(cu, ln)
                kept_for_sample.append((rec.header, rec.seq, rec.line_width))
            if not kept_for_sample:
                empty_samples += 1
                continue
            rel_dst = Path("filtered") / f"{sid}.fa"
            _write_deterministic_fa(kept_for_sample, work / rel_dst)
            objects.append({"path": str((out_dir / rel_dst).resolve()), **SPIRE_SOURCE, "retrieved_at": retrieved})

    sorted_hex = sorted(cu.hex() for cu in contig_lens)
    genomes = [
        {
            "contig_uids": [h],
            "genetic_code": 11,
            "kind": "prokaryote",
            "circular": [],
            "assembly_quality": "contig",
        }
        for h in sorted_hex
    ]

    (work / "objects.json").write_text(json.dumps(objects, indent=1) + "\n")
    (work / "genomes.json").write_text(json.dumps(genomes) + "\n")
    cross_stage_rows.sort()
    (work / "cross_stage_occurrences.tsv").write_text(
        "sample_id\tcontig_uid\taccession\tcontig_len\n"
        + "".join(f"{sid}\t{ch}\t{acc}\t{ln}\n" for sid, ch, acc, ln in cross_stage_rows)
    )

    report = {
        "samples_in": len(ordered),
        "samples_with_kept_contigs": len(objects),
        "empty_or_filtered_out_samples": empty_samples,
        "min_contig_len": min_contig_len,
        "total_records": total_records,
        "kept_records": kept_records,
        "dropped_records": dropped_records,
        "held_by_earlier_stage_records": held_records,
        "total_bp": total_bp,
        "kept_bp": kept_bp,
        "dropped_bp": dropped_bp,
        "held_by_earlier_stage_bp": held_bp,
        "distinct_contigs": len(genomes),
        "distinct_contig_bp": sum(contig_lens.values()),
    }
    (work / "input_report.json").write_text(json.dumps(report, indent=2) + "\n")
    work.rename(out_dir)
    return report


HAB_SHORT = {
    "host-associated": "H",
    "aquatic": "Q",
    "terrestrial": "T",
    "anthropogenic": "N",
    "air": "A",
    "unlabelled": "U",
}


def partition_samples_into_stages(
    samples: list[dict],
    stats_by_sample: dict[str, dict],
    *,
    min_contig_len: int = 0,
    max_stage_bp: int = 25_000_000_000,
    max_stage_contigs: int = 12_000_000,
) -> list[tuple[str, str, list[dict]]]:
    """Group samples by habitat into stages `SP<hab><k>` within `max_stage_bp` and `max_stage_contigs`."""
    by_hab: dict[str, list[dict]] = defaultdict(list)
    for r in sorted(samples, key=lambda x: (x["habitat"], int(x.get("rank_in_stratum", 0)), x["sample_id"])):
        st = stats_by_sample.get(r["sample_id"], {})
        if st.get("records", 0) == 0:
            continue
        by_hab[r["habitat"]].append(r)

    stages: list[tuple[str, str, list[dict]]] = []
    for hab in sorted(by_hab):
        tag = HAB_SHORT.get(hab, "X")
        cur: list[dict] = []
        cur_bp = cur_n = 0
        idx = 1
        for r in by_hab[hab]:
            st = stats_by_sample[r["sample_id"]]
            if min_contig_len >= 1500:
                s_bp, s_n = int(st["bp_ge_1500"]), int(st["n_ge_1500"])
            elif min_contig_len >= 1000:
                s_bp, s_n = int(st["bp_ge_1000"]), int(st["n_ge_1000"])
            else:
                s_bp, s_n = int(st["total_bp"]), int(st["records"])
            if s_n == 0:
                continue
            if cur and (cur_bp + s_bp > max_stage_bp or cur_n + s_n > max_stage_contigs):
                stages.append((f"SP{tag}{idx}", hab, cur))
                idx += 1
                cur, cur_bp, cur_n = [], 0, 0
            cur.append(r)
            cur_bp += s_bp
            cur_n += s_n
        if cur:
            stages.append((f"SP{tag}{idx}", hab, cur))
    return stages


def cmd_plan_stages(args: argparse.Namespace) -> int:
    sel = []
    with open(args.selection, newline="") as fh:
        sel = list(csv.DictReader(fh, delimiter="\t"))
    stats_by_s = {}
    with open(args.sample_stats, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            stats_by_s[r["sample_id"]] = r

    retrieved_by_s = {}
    if args.fetch_log and Path(args.fetch_log).is_file():
        with open(args.fetch_log) as fh:
            for line in fh:
                if line.strip():
                    rec = json.loads(line)
                    if rec.get("sample_id") and rec.get("fetched_utc"):
                        retrieved_by_s[rec["sample_id"]] = rec["fetched_utc"]

    for r in sel:
        r["retrieved_at"] = retrieved_by_s.get(r["sample_id"], DEFAULT_RETRIEVED)

    fdir = Path(args.fasta_dir)
    fastas = {r["sample_id"]: fdir / f"{r['sample_id']}.fa.gz" for r in sel}
    stages = partition_samples_into_stages(
        sel,
        stats_by_s,
        min_contig_len=args.min_contig_len,
        max_stage_bp=args.max_stage_bp,
        max_stage_contigs=args.max_stage_contigs,
    )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for st_name, hab, st_samples in stages:
        st_out = out_dir / st_name
        rep = build_stage_inputs(st_samples, fastas, st_out, min_contig_len=args.min_contig_len)
        manifest_rows.append({"stage": st_name, "habitat": hab, **rep})
        (st_out / "samples.tsv").write_text(
            "sample_id\thabitat\thabitat2\tstudy\n"
            + "".join(f"{s['sample_id']}\t{s['habitat']}\t{s['habitat2']}\t{s['study']}\n" for s in st_samples)
        )

    (out_dir / "stages_manifest.json").write_text(json.dumps(manifest_rows, indent=2) + "\n")
    print(json.dumps({"stages": len(manifest_rows), "details": manifest_rows}))
    return 0


def _dir_bytes(base: Path) -> int:
    if not base.exists():
        return 0
    if base.is_file():
        return base.stat().st_size
    return sum(p.stat().st_size for p in base.rglob("*") if p.is_file())


def summarize_depth1_root(scratch_root: Path, stages_manifest: list[dict], work_tmp: Path) -> dict:
    """Extract step-by-step CPU-minutes, bytes/Gbp, genes/Gbp, and fragment share for depth-1 stages."""
    scratch_root = Path(scratch_root)
    work_tmp.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(config={"threads": 1, "memory_limit": "4GB", "temp_directory": str(work_tmp)})

    per_stage: list[dict] = []
    try:
        all_ci = []
        all_occ = []
        for st_info in stages_manifest:
            st = st_info["stage"]
            sr = stage_root(scratch_root, st)
            all_ci.extend(io.part_files(sr / "l0seq" / "contig_index" / f"stage={st}"))
            all_occ.extend(sorted((sr / "l1" / "protein_occurrence").glob(f"gene_call_version=*/stage={st}/part-*.parquet")))

        band_rows = []
        if all_ci and all_occ:
            ci_all_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in all_ci)
            occ_all_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in all_occ)
            raw_bands = con.execute(
                f"""
                WITH ci AS (
                    SELECT
                        contig_uid,
                        max(contig_len)::BIGINT AS clen,
                        CASE
                            WHEN max(contig_len) < 500 THEN '1_200_to_499bp'
                            WHEN max(contig_len) < 1000 THEN '2_500_to_999bp'
                            WHEN max(contig_len) < 1500 THEN '3_1000_to_1499bp'
                            WHEN max(contig_len) < 5000 THEN '4_1500_to_4999bp'
                            WHEN max(contig_len) < 10000 THEN '5_5000_to_9999bp'
                            ELSE '6_ge_10000bp'
                        END AS band
                    FROM read_parquet([{ci_all_lits}], hive_partitioning=false)
                    WHERE is_default_variant
                    GROUP BY contig_uid
                ),
                c_agg AS (
                    SELECT band, count(*)::BIGINT AS contigs, coalesce(sum(clen), 0)::BIGINT AS bp
                    FROM ci
                    GROUP BY band
                ),
                g_agg AS (
                    SELECT
                        ci.band,
                        count(*)::BIGINT AS genes,
                        coalesce(sum(CASE WHEN o.dist_to_contig_end_bp = 0 THEN 1 ELSE 0 END), 0)::BIGINT AS frag_genes
                    FROM read_parquet([{occ_all_lits}], hive_partitioning=false) o
                    JOIN ci USING (contig_uid)
                    GROUP BY ci.band
                )
                SELECT
                    c.band,
                    c.contigs,
                    c.bp,
                    coalesce(g.genes, 0)::BIGINT,
                    coalesce(g.frag_genes, 0)::BIGINT
                FROM c_agg c
                LEFT JOIN g_agg g USING (band)
                ORDER BY c.band
                """
            ).fetchall()
            for bname, cnt, bp, g, fg in raw_bands:
                gbp_b = bp / 1e9 if bp else 0.0
                band_rows.append({
                    "band": bname,
                    "contigs": int(cnt),
                    "bp": int(bp),
                    "genes": int(g),
                    "genes_per_gbp": round(g / gbp_b, 1) if gbp_b else 0.0,
                    "fragment_genes": int(fg),
                    "fragment_share": round(fg / g, 6) if g else 0.0,
                })

        for st_info in stages_manifest:
            st = st_info["stage"]
            hab = st_info["habitat"]
            sr = stage_root(scratch_root, st)
            runs_dir = sr / "runs" / f"stage={st}"
            cpu_by_cat = {"inputs": 0.0, "ingest": 0.0, "genecall": 0.0, "crispr": 0.0, "igr": 0.0, "tapes": 0.0}
            wall_by_cat = {"inputs": 0.0, "ingest": 0.0, "genecall": 0.0, "crispr": 0.0, "igr": 0.0, "tapes": 0.0}
            peak_rss_mb = 0.0
            has_crispr_run = False
            for rf in sorted(runs_dir.glob("*.json")):
                if rf.name.endswith(".publishing.json"):
                    continue
                rec = json.loads(rf.read_text())
                step_name = rec.get("step", rf.stem)
                cpu_m = float(rec["cpu_min"]) if "cpu_min" in rec else float(rec.get("cpu_seconds", 0.0)) / 60.0
                wall_s = float(rec["wall_s"]) if "wall_s" in rec else float(rec.get("wall_seconds", 0.0))
                rss_mb = float(rec.get("peak_rss_bytes") or rec.get("max_rss_bytes") or 0) / (1024 * 1024)
                peak_rss_mb = max(peak_rss_mb, rss_mb)
                if step_name.startswith("crispr"):
                    has_crispr_run = True
                for k in cpu_by_cat:
                    if step_name == k or step_name.startswith(k + ".") or step_name.startswith(k + "-"):
                        cpu_by_cat[k] += cpu_m
                        wall_by_cat[k] += wall_s
                        break

            # Stored bases & new contigs in this stage's contig_index (default variants)
            ci_files = io.part_files(sr / "l0seq" / "contig_index" / f"stage={st}")
            ci_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in ci_files)
            owned_contigs, stored_bp = con.execute(
                f"SELECT count(*)::BIGINT, coalesce(sum(contig_len), 0)::BIGINT "
                f"FROM read_parquet([{ci_lits}], hive_partitioning=false) WHERE is_default_variant"
            ).fetchone()

            # Genes, fragments (dist_to_contig_end_bp == 0), distinct proteins, mean protein length
            occ_files = sorted((sr / "l1" / "protein_occurrence").glob(f"gene_call_version=*/stage={st}/part-*.parquet"))
            pseq_files = sorted((sr / "l0seq" / "protein_seq").glob(f"gene_call_version=*/stage={st}/protein_prefix=*/part-*.parquet"))
            crispr_files = sorted((sr / "l1" / "crispr").glob(f"producer=*/stage={st}/part-*.parquet"))

            genes = frag_genes = 0
            if occ_files:
                occ_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in occ_files)
                genes, frag_genes = con.execute(
                    f"SELECT count(*)::BIGINT, coalesce(sum(CASE WHEN dist_to_contig_end_bp = 0 THEN 1 ELSE 0 END), 0)::BIGINT "
                    f"FROM read_parquet([{occ_lits}], hive_partitioning=false)"
                ).fetchone()

            distinct_prots = total_aa = 0
            if pseq_files:
                ps_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for p in pseq_files for f in [p])
                distinct_prots, total_aa = con.execute(
                    f"SELECT count(*)::BIGINT, coalesce(sum(length), 0)::BIGINT "
                    f"FROM read_parquet([{ps_lits}], hive_partitioning=false)"
                ).fetchone()

            crispr_arrays = 0
            if crispr_files:
                cr_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in crispr_files)
                (crispr_arrays,) = con.execute(
                    f"SELECT count(DISTINCT feature_uid)::BIGINT FROM read_parquet([{cr_lits}], hive_partitioning=false)"
                ).fetchone()

            # File bytes by layer
            b_shards = _dir_bytes(sr / "l0seq" / "seq_shards" / f"stage={st}")
            b_d0_tables = sum(
                _dir_bytes(sr / rel / f"stage={st}")
                for rel in (
                    "l0/l0_manifest",
                    "l0seq/contig_index",
                    "l0seq/contig_label",
                    "l0seq/contig_occurrence",
                    "l0seq/genome_input",
                    "l0seq/seq_locator",
                )
            ) + _dir_bytes(sr / "ledger" / f"stage={st}" / "step=ingest")
            b_igr = sum(_dir_bytes(p) for p in (sr / "l0seq" / "igr").glob(f"igr_version=*/stage={st}"))
            b_d1_other = (
                sum(_dir_bytes(p) for p in (sr / "l0seq" / "protein_seq").glob(f"gene_call_version=*/stage={st}"))
                + sum(_dir_bytes(p) for p in (sr / "l1" / "feature").glob(f"producer=*/stage={st}"))
                + sum(_dir_bytes(p) for p in (sr / "l1" / "feature_protein").glob(f"gene_call_version=*/stage={st}"))
                + sum(_dir_bytes(p) for p in (sr / "l1" / "protein_occurrence").glob(f"gene_call_version=*/stage={st}"))
                + sum(_dir_bytes(p) for p in (sr / "l1" / "crispr").glob(f"producer=*/stage={st}"))
                + sum(_dir_bytes(p) for p in (sr / "l3" / "tapes").glob(f"tape_version=*/stage={st}"))
                + sum(_dir_bytes(p) for p in (sr / "l3" / "contig_tapes").glob(f"tape_version=*/stage={st}"))
                + sum(
                    _dir_bytes(p)
                    for p in (sr / "ledger" / f"stage={st}").glob("step=*")
                    if p.name != "step=ingest"
                )
            )
            b_stored = b_shards + b_d0_tables
            b_called_added = b_igr + b_d1_other
            b_total = b_stored + b_called_added
            b_without_nt_igr = b_total - b_shards - b_igr

            gbp = stored_bp / 1e9 if stored_bp else 0.0
            fasta_gbp = st_info["kept_bp"] / 1e9 if st_info.get("kept_bp") else gbp
            per_stage.append({
                "stage": st,
                "habitat": hab,
                "samples": st_info["samples_with_kept_contigs"],
                "min_contig_len": st_info["min_contig_len"],
                "input_bp": st_info["kept_bp"],
                "owned_contigs": int(owned_contigs),
                "stored_bp": int(stored_bp),
                "genes": int(genes),
                "fragment_genes": int(frag_genes),
                "fragment_share": round(frag_genes / genes, 6) if genes else 0.0,
                "genes_per_gbp": round(genes / gbp, 1) if gbp else 0.0,
                "distinct_proteins": int(distinct_prots),
                "distinct_proteins_per_gbp": round(distinct_prots / gbp, 1) if gbp else 0.0,
                "mean_protein_aa": round(total_aa / distinct_prots, 2) if distinct_prots else 0.0,
                "crispr_ran": has_crispr_run,
                "crispr_arrays": int(crispr_arrays) if has_crispr_run else None,
                "crispr_per_gbp": round(crispr_arrays / gbp, 2) if (has_crispr_run and gbp) else None,
                "cpu_min_inputs": round(cpu_by_cat["inputs"], 3),
                "cpu_min_ingest": round(cpu_by_cat["ingest"], 3),
                "cpu_min_genecall": round(cpu_by_cat["genecall"], 3),
                "cpu_min_crispr": round(cpu_by_cat["crispr"], 3),
                "cpu_min_igr": round(cpu_by_cat["igr"], 3),
                "cpu_min_tapes": round(cpu_by_cat["tapes"], 3),
                "cpu_min_total": round(sum(cpu_by_cat.values()), 3),
                "ingest_cpu_min_per_input_gbp": round(cpu_by_cat["ingest"] / fasta_gbp, 3) if fasta_gbp else 0.0,
                "ingest_cpu_min_per_stored_gbp": round(cpu_by_cat["ingest"] / gbp, 3) if gbp else 0.0,
                "genecall_cpu_min_per_gbp": round(cpu_by_cat["genecall"] / gbp, 3) if gbp else 0.0,
                "crispr_cpu_min_per_gbp": round(cpu_by_cat["crispr"] / gbp, 3) if gbp else 0.0,
                "igr_cpu_min_per_gbp": round(cpu_by_cat["igr"] / gbp, 3) if gbp else 0.0,
                "tapes_cpu_min_per_gbp": round(cpu_by_cat["tapes"] / gbp, 3) if gbp else 0.0,
                "bytes_seq_shards": b_shards,
                "bytes_depth0_stored": b_stored,
                "bytes_igr": b_igr,
                "bytes_depth1_added": b_called_added,
                "bytes_depth1_total": b_total,
                "bytes_without_nt_and_igr": b_without_nt_igr,
                "gb_per_gbp_stored": round((b_stored / 1e9) / gbp, 4) if gbp else 0.0,
                "gb_per_gbp_called_total": round((b_total / 1e9) / gbp, 4) if gbp else 0.0,
                "gb_per_gbp_without_nt_igr": round((b_without_nt_igr / 1e9) / gbp, 4) if gbp else 0.0,
                "peak_rss_mb": round(peak_rss_mb, 1),
            })
    finally:
        con.close()

    tot_stored_bp = sum(r["stored_bp"] for r in per_stage)
    tot_input_bp = sum(r["input_bp"] for r in per_stage)
    tot_gbp = tot_stored_bp / 1e9 if tot_stored_bp else 0.0
    tot_in_gbp = tot_input_bp / 1e9 if tot_input_bp else 0.0
    tot_genes = sum(r["genes"] for r in per_stage)
    tot_frags = sum(r["fragment_genes"] for r in per_stage)
    tot_b_stored = sum(r["bytes_depth0_stored"] for r in per_stage)
    tot_b_total = sum(r["bytes_depth1_total"] for r in per_stage)
    tot_b_no_nt = sum(r["bytes_without_nt_and_igr"] for r in per_stage)
    tot_ing = sum(r["cpu_min_ingest"] for r in per_stage)
    tot_gc = sum(r["cpu_min_genecall"] for r in per_stage)
    tot_cr = sum(r["cpu_min_crispr"] for r in per_stage)
    tot_igr = sum(r["cpu_min_igr"] for r in per_stage)
    tot_tp = sum(r["cpu_min_tapes"] for r in per_stage)
    tot_inp = sum(r["cpu_min_inputs"] for r in per_stage)

    crispr_stages = [r for r in per_stage if r["crispr_ran"]]
    cr_gbp = sum(r["stored_bp"] for r in crispr_stages) / 1e9 if crispr_stages else 0.0
    cr_arrays_tot = sum(r["crispr_arrays"] or 0 for r in crispr_stages)

    overall = {
        "stages": len(per_stage),
        "samples": sum(r["samples"] for r in per_stage),
        "input_bp": tot_input_bp,
        "stored_bp": tot_stored_bp,
        "owned_contigs": sum(r["owned_contigs"] for r in per_stage),
        "genes": tot_genes,
        "fragment_genes": tot_frags,
        "fragment_share": round(tot_frags / tot_genes, 6) if tot_genes else 0.0,
        "genes_per_gbp": round(tot_genes / tot_gbp, 1) if tot_gbp else 0.0,
        "stage_sum_distinct_proteins": sum(r["distinct_proteins"] for r in per_stage),
        "crispr_stages_count": len(crispr_stages),
        "crispr_stored_bp": sum(r["stored_bp"] for r in crispr_stages),
        "crispr_arrays": cr_arrays_tot,
        "crispr_per_gbp": round(cr_arrays_tot / cr_gbp, 2) if cr_gbp else None,
        "cpu_min_inputs": round(tot_inp, 3),
        "cpu_min_ingest": round(tot_ing, 3),
        "cpu_min_genecall": round(tot_gc, 3),
        "cpu_min_crispr": round(tot_cr, 3),
        "cpu_min_igr": round(tot_igr, 3),
        "cpu_min_tapes": round(tot_tp, 3),
        "cpu_min_total": round(tot_inp + tot_ing + tot_gc + tot_cr + tot_igr + tot_tp, 3),
        "ingest_cpu_min_per_input_gbp": round(tot_ing / tot_in_gbp, 3) if tot_in_gbp else 0.0,
        "ingest_cpu_min_per_stored_gbp": round(tot_ing / tot_gbp, 3) if tot_gbp else 0.0,
        "genecall_cpu_min_per_gbp": round(tot_gc / tot_gbp, 3) if tot_gbp else 0.0,
        "crispr_cpu_min_per_gbp": round(tot_cr / cr_gbp, 3) if cr_gbp else None,
        "igr_cpu_min_per_gbp": round(tot_igr / tot_gbp, 3) if tot_gbp else 0.0,
        "tapes_cpu_min_per_gbp": round(tot_tp / tot_gbp, 3) if tot_gbp else 0.0,
        "gb_per_gbp_stored": round((tot_b_stored / 1e9) / tot_gbp, 4) if tot_gbp else 0.0,
        "gb_per_gbp_called_total": round((tot_b_total / 1e9) / tot_gbp, 4) if tot_gbp else 0.0,
        "gb_per_gbp_without_nt_igr": round((tot_b_no_nt / 1e9) / tot_gbp, 4) if tot_gbp else 0.0,
        "survey_s4_1_comparison": {
            "ingest_survey_cpu_min_per_gbp": 4.5,
            "ingest_b52_cpu_min_per_gbp": 4.11,
            "ingest_measured_vs_survey_ratio": round((tot_ing / tot_gbp) / 4.5, 4) if tot_gbp else None,
            "genecall_survey_range_cpu_min_per_gbp": [60.0, 91.0],
            "genecall_c38_linear_meta_cpu_min_per_gbp": 60.4,
            "genecall_measured_vs_c38_ratio": round((tot_gc / tot_gbp) / 60.4, 4) if tot_gbp else None,
            "crispr_survey_cpu_min_per_gbp": 3.0,
            "crispr_b53_cpu_min_per_gbp": round(168.1 / 60.0, 4),
            "crispr_measured_vs_survey_ratio": round((tot_cr / cr_gbp) / 3.0, 4) if cr_gbp else None,
            "storage_survey_gb_per_gbp": 0.71,
            "storage_measured_vs_survey_ratio": round(((tot_b_total / 1e9) / tot_gbp) / 0.71, 4) if tot_gbp else None,
            "storage_no_nt_survey_gb_per_gbp": 0.37,
            "storage_no_nt_measured_vs_survey_ratio": round(((tot_b_no_nt / 1e9) / tot_gbp) / 0.37, 4) if tot_gbp else None,
        },
    }
    return {"overall": overall, "by_length_band": band_rows, "per_stage": per_stage}


def cmd_summarize_depth1(args: argparse.Namespace) -> int:
    manifest = json.loads(Path(args.stages_manifest).read_text())
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    res = summarize_depth1_root(Path(args.scratch_root), manifest, out / "duck_tmp")
    shutil.rmtree(out / "duck_tmp", ignore_errors=True)
    (out / "depth1_summary.json").write_text(json.dumps(res, indent=2) + "\n")
    if res["per_stage"]:
        cols = list(res["per_stage"][0].keys())
        with open(out / "depth1_by_stage.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(res["per_stage"])
    print(json.dumps(res["overall"]))
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_plan = sub.add_parser("plan-stages")
    p_plan.add_argument("--selection", required=True)
    p_plan.add_argument("--sample-stats", required=True)
    p_plan.add_argument("--fasta-dir", required=True)
    p_plan.add_argument("--fetch-log", default=None)
    p_plan.add_argument("--out-dir", required=True)
    p_plan.add_argument("--min-contig-len", type=int, default=0)
    p_plan.add_argument("--max-stage-bp", type=int, default=25_000_000_000)
    p_plan.add_argument("--max-stage-contigs", type=int, default=12_000_000)

    p_sum = sub.add_parser("summarize-depth1")
    p_sum.add_argument("--scratch-root", required=True)
    p_sum.add_argument("--stages-manifest", required=True)
    p_sum.add_argument("--out-dir", required=True)

    args = ap.parse_args(argv)
    if args.cmd == "plan-stages":
        return cmd_plan_stages(args)
    if args.cmd == "summarize-depth1":
        return cmd_summarize_depth1(args)
    return 2


if __name__ == "__main__":
    sys.exit(main())
