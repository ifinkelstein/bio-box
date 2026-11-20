"""Card C-42 step (3): streaming contig length distributions, bin shares, and exact `contig_uid` repeats.

    python -I c42_contigs.py run --selection /tmp/if2453/s3/spire/select/pilot_selection.tsv
        --fasta-dir /tmp/if2453/s3/spire/dl/assemblies --work-dir /tmp/if2453/s3/spire/contigs
        [--root /tmp/if2453/r10/root] [--duck-mem 32GB] [--threads 10] [--reuse-chunks]

Bounded memory & zero spill (AGENTS §4.1):
1. Streams one sample FASTA at a time per worker, flushing compact zstd Parquet row groups of
   `(contig_uid, canon_uid, sample_idx, hab_idx, length)` to `work_dir/chunks/`.
2. Aggregates `contig_uid` and `canon_uid` repeats across 16 disjoint hex-nibble prefix buckets
   (`0x00..0x0f`, ..., `0xf0..0xff`), keeping each hash table ~1/16th of corpus size (< 4 GB even at
   450M contigs) so DuckDB never spills tens of GB to `/tmp`.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from biobox import fasta, ids
from biobox.l0seq import store

CHUNK_SCHEMA = pa.schema([
    ("contig_uid", pa.binary(16)),
    ("canon_uid", pa.binary(16)),
    ("sample_idx", pa.uint32()),
    ("hab_idx", pa.uint8()),
    ("length", pa.uint64()),
])
FLUSH_ROWS = 250_000
_RC_TRANS = bytes.maketrans(b"ACGTRYKMSWBDHVNacgtrykmswbdhvn", b"TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def _canon_uid(seq: bytes, cu: bytes) -> bytes:
    rc_cu = ids.contig_uid(seq.translate(_RC_TRANS)[::-1])
    return cu if cu <= rc_cu else rc_cu


def _n50_hist(len_counts: list[tuple[int, int]]) -> int:
    if not len_counts:
        return 0
    s = sorted(len_counts, key=lambda x: x[0], reverse=True)
    half = sum(ln * c for ln, c in s) / 2.0
    acc = 0
    for ln, c in s:
        acc += ln * c
        if acc >= half:
            return ln
    return s[-1][0]


def _empty_summary() -> dict:
    return {
        "records": 0, "distinct_contigs": 0, "within_sample_dup_records": 0,
        "total_bp": 0, "min_len": 0, "max_len": 0, "n50": 0,
        "n_ge_500": 0, "bp_ge_500": 0, "n_ge_1000": 0, "bp_ge_1000": 0,
        "n_ge_1500": 0, "bp_ge_1500": 0, "n_ge_5000": 0, "bp_ge_5000": 0,
        "n_ge_10000": 0, "bp_ge_10000": 0, "n_ge_100000": 0, "bp_ge_100000": 0,
    }


def _scan_worker_shard(tasks: list[tuple[int, int, str, str | None]], out_parquet: Path) -> list[tuple[int, dict]]:
    """Stream assigned sample FASTAs, flushing `(contig_uid, canon_uid, sample_idx, hab_idx, length)` to Parquet."""
    buf_cu: list[bytes] = []
    buf_can: list[bytes] = []
    buf_sidx: list[int] = []
    buf_hidx: list[int] = []
    buf_len: list[int] = []
    writer: pq.ParquetWriter | None = None

    def flush() -> None:
        nonlocal writer
        if not buf_cu:
            return
        if writer is None:
            writer = pq.ParquetWriter(out_parquet, CHUNK_SCHEMA, compression="zstd")
        tbl = pa.Table.from_arrays(
            [
                pa.array(buf_cu, type=pa.binary(16)),
                pa.array(buf_can, type=pa.binary(16)),
                pa.array(buf_sidx, type=pa.uint32()),
                pa.array(buf_hidx, type=pa.uint8()),
                pa.array(buf_len, type=pa.uint64()),
            ],
            schema=CHUNK_SCHEMA,
        )
        writer.write_table(tbl)
        buf_cu.clear()
        buf_can.clear()
        buf_sidx.clear()
        buf_hidx.clear()
        buf_len.clear()

    summaries: list[tuple[int, dict]] = []
    try:
        for s_idx, h_idx, _sid, fa_str in tasks:
            if not fa_str or not Path(fa_str).is_file():
                summaries.append((s_idx, _empty_summary()))
                continue
            seen: set[bytes] = set()
            len_hist: dict[int, int] = defaultdict(int)
            n_rec = bp_total = bp_500 = bp_1000 = bp_1500 = bp_5000 = bp_10000 = bp_100000 = 0
            n_500 = n_1000 = n_1500 = n_5000 = n_10000 = n_100000 = 0

            for rec in fasta.iter_records(Path(fa_str)):
                ln = len(rec.seq)
                cu = ids.contig_uid(rec.seq)
                can = _canon_uid(rec.seq, cu)
                seen.add(cu)
                len_hist[ln] += 1
                n_rec += 1
                bp_total += ln
                if ln >= 500:
                    n_500 += 1
                    bp_500 += ln
                if ln >= 1000:
                    n_1000 += 1
                    bp_1000 += ln
                if ln >= 1500:
                    n_1500 += 1
                    bp_1500 += ln
                if ln >= 5000:
                    n_5000 += 1
                    bp_5000 += ln
                if ln >= 10000:
                    n_10000 += 1
                    bp_10000 += ln
                if ln >= 100000:
                    n_100000 += 1
                    bp_100000 += ln

                buf_cu.append(cu)
                buf_can.append(can)
                buf_sidx.append(s_idx)
                buf_hidx.append(h_idx)
                buf_len.append(ln)
                if len(buf_cu) >= FLUSH_ROWS:
                    flush()

            summaries.append((
                s_idx,
                {
                    "records": n_rec,
                    "distinct_contigs": len(seen),
                    "within_sample_dup_records": n_rec - len(seen),
                    "total_bp": bp_total,
                    "min_len": min(len_hist) if len_hist else 0,
                    "max_len": max(len_hist) if len_hist else 0,
                    "n50": _n50_hist(list(len_hist.items())),
                    "n_ge_500": n_500, "bp_ge_500": bp_500,
                    "n_ge_1000": n_1000, "bp_ge_1000": bp_1000,
                    "n_ge_1500": n_1500, "bp_ge_1500": bp_1500,
                    "n_ge_5000": n_5000, "bp_ge_5000": bp_5000,
                    "n_ge_10000": n_10000, "bp_ge_10000": bp_10000,
                    "n_ge_100000": n_100000, "bp_ge_100000": bp_100000,
                },
            ))
        flush()
    finally:
        if writer is not None:
            writer.close()
    return summaries


def _bucket_where(col: str, b: int) -> str:
    lo = f"from_hex('{b * 16:02x}')"
    if b == 15:
        return f"{col} >= {lo}"
    hi = f"from_hex('{(b + 1) * 16:02x}')"
    return f"({col} >= {lo} AND {col} < {hi})"


def _reconstruct_sample_summaries(
    con: duckdb.DuckDBPyConnection,
    lits: str,
    n_samples: int,
    distinct_by_sidx: dict[int, int],
) -> dict[int, dict]:
    """Reconstruct exact per-sample length stats & N50 from `chunks/part-*.parquet` without re-reading FASTAs."""
    out = {i: _empty_summary() for i in range(n_samples)}
    hist_by_s: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for s_idx, ln, cnt in con.execute(
        f"SELECT sample_idx::INTEGER, length::BIGINT, count(*)::BIGINT FROM read_parquet([{lits}], hive_partitioning=false) GROUP BY sample_idx, length"
    ).fetchall():
        hist_by_s[int(s_idx)].append((int(ln), int(cnt)))

    for s_idx, items in hist_by_s.items():
        n_rec = sum(c for _, c in items)
        dist_c = distinct_by_sidx.get(s_idx, 0)
        tot_bp = sum(ln * c for ln, c in items)
        out[s_idx] = {
            "records": n_rec,
            "distinct_contigs": dist_c,
            "within_sample_dup_records": n_rec - dist_c,
            "total_bp": tot_bp,
            "min_len": min(ln for ln, _ in items),
            "max_len": max(ln for ln, _ in items),
            "n50": _n50_hist(items),
            "n_ge_500": sum(c for ln, c in items if ln >= 500),
            "bp_ge_500": sum(ln * c for ln, c in items if ln >= 500),
            "n_ge_1000": sum(c for ln, c in items if ln >= 1000),
            "bp_ge_1000": sum(ln * c for ln, c in items if ln >= 1000),
            "n_ge_1500": sum(c for ln, c in items if ln >= 1500),
            "bp_ge_1500": sum(ln * c for ln, c in items if ln >= 1500),
            "n_ge_5000": sum(c for ln, c in items if ln >= 5000),
            "bp_ge_5000": sum(ln * c for ln, c in items if ln >= 5000),
            "n_ge_10000": sum(c for ln, c in items if ln >= 10000),
            "bp_ge_10000": sum(ln * c for ln, c in items if ln >= 10000),
            "n_ge_100000": sum(c for ln, c in items if ln >= 100000),
            "bp_ge_100000": sum(ln * c for ln, c in items if ln >= 100000),
        }
    return out


def _agg_habitat(per_sample: list[dict], key: str) -> list[dict]:
    by_h: dict[str, list[dict]] = defaultdict(list)
    for r in per_sample:
        by_h[r[key]].append(r)
    out: list[dict] = []
    for h in sorted(by_h):
        rs = by_h[h]
        tot = sum(r["total_bp"] for r in rs)
        mag_bp = sum(int(r.get("mag_bp") or 0) for r in rs)
        bf_num = sum(
            r["total_bp"] * float(r["binned_fraction"])
            for r in rs
            if r.get("binned_fraction") not in (None, "") and r["total_bp"] > 0
        )
        bf_den = sum(
            r["total_bp"]
            for r in rs
            if r.get("binned_fraction") not in (None, "") and r["total_bp"] > 0
        )
        b1000 = sum(r["bp_ge_1000"] for r in rs)
        b1500 = sum(r["bp_ge_1500"] for r in rs)
        b5000 = sum(r["bp_ge_5000"] for r in rs)
        b10000 = sum(r["bp_ge_10000"] for r in rs)
        b100000 = sum(r["bp_ge_100000"] for r in rs)
        out.append({
            key: h,
            "samples": len(rs),
            "nonempty_samples": sum(1 for r in rs if r["records"] > 0),
            "records": sum(r["records"] for r in rs),
            "total_bp": tot,
            "n_ge_1000": sum(r["n_ge_1000"] for r in rs),
            "bp_ge_1000": b1000,
            "share_ge_1000": round(b1000 / tot, 6) if tot else 0.0,
            "n_ge_1500": sum(r["n_ge_1500"] for r in rs),
            "bp_ge_1500": b1500,
            "share_ge_1500": round(b1500 / tot, 6) if tot else 0.0,
            "n_ge_5000": sum(r["n_ge_5000"] for r in rs),
            "bp_ge_5000": b5000,
            "share_ge_5000": round(b5000 / tot, 6) if tot else 0.0,
            "n_ge_10000": sum(r["n_ge_10000"] for r in rs),
            "bp_ge_10000": b10000,
            "share_ge_10000": round(b10000 / tot, 6) if tot else 0.0,
            "n_ge_100000": sum(r["n_ge_100000"] for r in rs),
            "bp_ge_100000": b100000,
            "share_ge_100000": round(b100000 / tot, 6) if tot else 0.0,
            "mag_bp": mag_bp,
            "share_in_spire_mags": round(mag_bp / tot, 6) if tot else 0.0,
            "weighted_metabat2_binned_fraction": round(bf_num / bf_den, 6) if bf_den else None,
        })
    return out


def analyze_pilot_contigs(
    sample_meta: list[dict],
    fasta_by_sample: dict[str, Path],
    work_dir: Path,
    *,
    root: Path | None = None,
    duck_mem: str = "4GB",
    threads: int = 1,
    reuse_chunks: bool = False,
    keep_chunks: bool = False,
) -> dict:
    """Stream all pilot sample FASTAs, compute per-habitat length/bin shares and exact `contig_uid` repeats."""
    work_dir = Path(work_dir)
    chunks_dir = work_dir / "chunks"
    spill_dir = work_dir / "duckdb_spill"
    spill_dir.mkdir(parents=True, exist_ok=True)

    ordered = sorted(sample_meta, key=lambda r: r["sample_id"])
    hab_names = sorted({r["habitat"] for r in ordered})
    hab_to_idx = {h: i for i, h in enumerate(hab_names)}

    summary_by_sidx: dict[int, dict] = {}
    existing_parts = sorted(chunks_dir.glob("part-*.parquet")) if chunks_dir.is_dir() else []

    if not (reuse_chunks and existing_parts):
        if chunks_dir.exists():
            shutil.rmtree(chunks_dir)
        chunks_dir.mkdir(parents=True, exist_ok=True)
        all_tasks: list[tuple[int, int, str, str | None]] = []
        for s_idx, meta in enumerate(ordered):
            sid = meta["sample_id"]
            fa = fasta_by_sample.get(sid)
            all_tasks.append((s_idx, hab_to_idx[meta["habitat"]], sid, str(fa) if (fa and fa.is_file()) else None))

        n_workers = max(1, min(threads, len(all_tasks) or 1))
        shards = [all_tasks[w::n_workers] for w in range(n_workers)]
        if n_workers == 1:
            for s_idx, s_sum in _scan_worker_shard(shards[0], chunks_dir / "part-00000.parquet"):
                summary_by_sidx[s_idx] = s_sum
        else:
            with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as pool:
                futs = [
                    pool.submit(_scan_worker_shard, shards[w], chunks_dir / f"part-{w:05d}.parquet")
                    for w in range(n_workers)
                ]
                for fut in futs:
                    for s_idx, s_sum in fut.result():
                        summary_by_sidx[s_idx] = s_sum

    part_files = sorted(chunks_dir.glob("part-*.parquet"))
    cu_all = [0] * 8
    cu_1500 = [0] * 8
    can_all = [0] * 5
    can_1500 = [0] * 5
    hab_acc: dict[int, list[int]] = {i: [0, 0, 0, 0] for i in range(len(hab_names))}
    distinct_by_sidx: dict[int, int] = defaultdict(int)
    root_overlap: dict[str, dict] = {}

    if part_files:
        lits = ", ".join("'" + str(p).replace("'", "''") + "'" for p in part_files)
        con = duckdb.connect(config={"threads": max(1, threads), "memory_limit": duck_mem, "temp_directory": str(spill_dir)})
        try:
            # 16 disjoint hex-nibble buckets on `contig_uid`: ~1/16th RAM per pass, zero spill!
            for b in range(16):
                bw = _bucket_where("contig_uid", b)
                r = con.execute(
                    f"""
                    WITH per_cu AS (
                        SELECT
                            contig_uid,
                            max(length)::BIGINT AS len,
                            count(*)::BIGINT AS n_rec,
                            count(DISTINCT sample_idx)::BIGINT AS n_samp,
                            count(DISTINCT hab_idx)::BIGINT AS n_hab
                        FROM read_parquet([{lits}], hive_partitioning = false)
                        WHERE {bw}
                        GROUP BY contig_uid
                    )
                    SELECT
                        coalesce(sum(n_rec), 0)::BIGINT,
                        count(*)::BIGINT,
                        coalesce(sum(len), 0)::BIGINT,
                        coalesce(sum(CASE WHEN n_samp > 1 THEN 1 ELSE 0 END), 0)::BIGINT,
                        coalesce(sum(CASE WHEN n_samp > 1 AND n_hab = 1 THEN 1 ELSE 0 END), 0)::BIGINT,
                        coalesce(sum(CASE WHEN n_hab > 1 THEN 1 ELSE 0 END), 0)::BIGINT,
                        coalesce(sum(n_samp - 1), 0)::BIGINT,
                        coalesce(sum((n_samp - 1) * len), 0)::BIGINT,
                        coalesce(sum(CASE WHEN len >= 1500 THEN n_rec ELSE 0 END), 0)::BIGINT,
                        coalesce(sum(CASE WHEN len >= 1500 THEN 1 ELSE 0 END), 0)::BIGINT,
                        coalesce(sum(CASE WHEN len >= 1500 THEN len ELSE 0 END), 0)::BIGINT,
                        coalesce(sum(CASE WHEN len >= 1500 AND n_samp > 1 THEN 1 ELSE 0 END), 0)::BIGINT,
                        coalesce(sum(CASE WHEN len >= 1500 AND n_samp > 1 AND n_hab = 1 THEN 1 ELSE 0 END), 0)::BIGINT,
                        coalesce(sum(CASE WHEN len >= 1500 AND n_hab > 1 THEN 1 ELSE 0 END), 0)::BIGINT,
                        coalesce(sum(CASE WHEN len >= 1500 THEN n_samp - 1 ELSE 0 END), 0)::BIGINT,
                        coalesce(sum(CASE WHEN len >= 1500 THEN (n_samp - 1) * len ELSE 0 END), 0)::BIGINT
                    FROM per_cu
                    """
                ).fetchone()
                for i in range(8):
                    cu_all[i] += int(r[i])
                    cu_1500[i] += int(r[8 + i])

                # Per-habitat repeats in this bucket
                for h_idx, dc, dbp, mc, rbp in con.execute(
                    f"""
                    WITH per_hab_cu AS (
                        SELECT hab_idx, contig_uid, max(length)::BIGINT AS len, count(DISTINCT sample_idx)::BIGINT AS n_samp
                        FROM read_parquet([{lits}], hive_partitioning = false)
                        WHERE {bw}
                        GROUP BY hab_idx, contig_uid
                    )
                    SELECT
                        hab_idx::INTEGER,
                        count(*)::BIGINT,
                        coalesce(sum(len), 0)::BIGINT,
                        coalesce(sum(CASE WHEN n_samp > 1 THEN 1 ELSE 0 END), 0)::BIGINT,
                        coalesce(sum((n_samp - 1) * len), 0)::BIGINT
                    FROM per_hab_cu
                    GROUP BY hab_idx
                    """
                ).fetchall():
                    acc = hab_acc[int(h_idx)]
                    acc[0] += int(dc)
                    acc[1] += int(dbp)
                    acc[2] += int(mc)
                    acc[3] += int(rbp)

                if not summary_by_sidx:
                    for s_idx, dc in con.execute(
                        f"""
                        SELECT sample_idx::INTEGER, count(DISTINCT contig_uid)::BIGINT
                        FROM read_parquet([{lits}], hive_partitioning = false)
                        WHERE {bw}
                        GROUP BY sample_idx
                        """
                    ).fetchall():
                        distinct_by_sidx[int(s_idx)] += int(dc)

            # 16 disjoint hex-nibble buckets on `canon_uid`
            for b in range(16):
                bw = _bucket_where("canon_uid", b)
                cr = con.execute(
                    f"""
                    WITH per_can AS (
                        SELECT
                            canon_uid,
                            max(length)::BIGINT AS len,
                            count(DISTINCT sample_idx)::BIGINT AS n_samp,
                            count(DISTINCT hab_idx)::BIGINT AS n_hab
                        FROM read_parquet([{lits}], hive_partitioning = false)
                        WHERE {bw}
                        GROUP BY canon_uid
                    )
                    SELECT
                        count(*)::BIGINT,
                        coalesce(sum(len), 0)::BIGINT,
                        coalesce(sum(CASE WHEN n_samp > 1 THEN 1 ELSE 0 END), 0)::BIGINT,
                        coalesce(sum(CASE WHEN n_hab > 1 THEN 1 ELSE 0 END), 0)::BIGINT,
                        coalesce(sum((n_samp - 1) * len), 0)::BIGINT,
                        coalesce(sum(CASE WHEN len >= 1500 THEN 1 ELSE 0 END), 0)::BIGINT,
                        coalesce(sum(CASE WHEN len >= 1500 THEN len ELSE 0 END), 0)::BIGINT,
                        coalesce(sum(CASE WHEN len >= 1500 AND n_samp > 1 THEN 1 ELSE 0 END), 0)::BIGINT,
                        coalesce(sum(CASE WHEN len >= 1500 AND n_hab > 1 THEN 1 ELSE 0 END), 0)::BIGINT,
                        coalesce(sum(CASE WHEN len >= 1500 THEN (n_samp - 1) * len ELSE 0 END), 0)::BIGINT
                    FROM per_can
                    """
                ).fetchone()
                for i in range(5):
                    can_all[i] += int(cr[i])
                    can_1500[i] += int(cr[5 + i])

            if not summary_by_sidx:
                summary_by_sidx = _reconstruct_sample_summaries(con, lits, len(ordered), distinct_by_sidx)

            if root is not None and Path(root).is_dir():
                idx_files = store.index_files(Path(root))
                if idx_files:
                    idx_lits = ", ".join("'" + str(p).replace("'", "''") + "'" for p in idx_files)
                    # Root has < 1M contigs (`S1..PL1`), so build hash table on `root_cu` (tiny!) and stream `chunks`!
                    ov_rows = con.execute(
                        f"""
                        WITH root_cu AS (
                            SELECT contig_uid, regexp_extract(filename, 'stage=([^/]+)', 1) AS stage
                            FROM read_parquet([{idx_lits}], hive_partitioning = false, filename = true)
                            WHERE is_default_variant
                        ),
                        matched AS (
                            SELECT DISTINCT p.contig_uid, p.length::BIGINT AS len, r.stage
                            FROM read_parquet([{lits}], hive_partitioning = false) p
                            JOIN root_cu r USING (contig_uid)
                        )
                        SELECT stage, count(*)::BIGINT, coalesce(sum(len), 0)::BIGINT
                        FROM matched
                        GROUP BY stage
                        ORDER BY stage
                        """
                    ).fetchall()
                    root_overlap = {
                        st: {"shared_contigs": int(cnt), "shared_bp": int(bp)}
                        for st, cnt, bp in ov_rows
                    }
        finally:
            con.close()

    if not summary_by_sidx:
        summary_by_sidx = {i: _empty_summary() for i in range(len(ordered))}

    per_sample: list[dict] = []
    for s_idx, meta in enumerate(ordered):
        per_sample.append({
            "sample_id": meta["sample_id"],
            "habitat": meta["habitat"],
            "habitat2": meta.get("habitat2", meta["habitat"]),
            "mag_bp": int(meta.get("mag_bp") or 0),
            "binned_fraction": meta.get("binned_fraction"),
            **summary_by_sidx[s_idx],
        })

    repeats = {
        "total_records": cu_all[0],
        "distinct_contig_uids": cu_all[1],
        "distinct_bp": cu_all[2],
        "contigs_in_multiple_samples": cu_all[3],
        "contigs_within_single_habitat_multi_sample": cu_all[4],
        "contigs_across_habitats": cu_all[5],
        "repeated_sample_occurrences": cu_all[6],
        "repeated_bp_between_samples": cu_all[7],
        "distinct_strand_canonical_uids": can_all[0],
        "distinct_strand_canonical_bp": can_all[1],
        "strand_canonical_contigs_in_multiple_samples": can_all[2],
        "strand_canonical_contigs_across_habitats": can_all[3],
        "strand_canonical_repeated_bp_between_samples": can_all[4],
        "within_sample_dup_records": sum(r["within_sample_dup_records"] for r in per_sample),
    }
    repeats_ge_1500 = {
        "total_records": cu_1500[0],
        "distinct_contig_uids": cu_1500[1],
        "distinct_bp": cu_1500[2],
        "contigs_in_multiple_samples": cu_1500[3],
        "contigs_within_single_habitat_multi_sample": cu_1500[4],
        "contigs_across_habitats": cu_1500[5],
        "repeated_sample_occurrences": cu_1500[6],
        "repeated_bp_between_samples": cu_1500[7],
        "distinct_strand_canonical_uids": can_1500[0],
        "distinct_strand_canonical_bp": can_1500[1],
        "strand_canonical_contigs_in_multiple_samples": can_1500[2],
        "strand_canonical_contigs_across_habitats": can_1500[3],
        "strand_canonical_repeated_bp_between_samples": can_1500[4],
    }
    by_hab_repeats = [
        {
            "habitat": hab_names[h_idx],
            "distinct_contigs": hab_acc[h_idx][0],
            "distinct_bp": hab_acc[h_idx][1],
            "multi_sample_contigs_in_habitat": hab_acc[h_idx][2],
            "repeated_bp_within_habitat": hab_acc[h_idx][3],
        }
        for h_idx in range(len(hab_names))
        if hab_acc[h_idx][0] > 0
    ]

    if not keep_chunks:
        shutil.rmtree(chunks_dir, ignore_errors=True)
    shutil.rmtree(spill_dir, ignore_errors=True)

    return {
        "per_sample": per_sample,
        "by_habitat": _agg_habitat(per_sample, "habitat"),
        "by_habitat2": _agg_habitat(per_sample, "habitat2"),
        "Overall": _agg_habitat([{**r, "_all": "all"} for r in per_sample], "_all")[0] if per_sample else {},
        "repeats": repeats,
        "repeats_ge_1500": repeats_ge_1500,
        "repeats_by_habitat": by_hab_repeats,
        "root_overlap": root_overlap,
    }


def _write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    cols = list(rows[0].keys())
    tmp = path.with_name(f".{path.name}.tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in cols})
    tmp.replace(path)


def cmd_run(args: argparse.Namespace) -> int:
    sel_rows = []
    with open(args.selection, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["mag_bp"] = int(r["mag_bp"] or 0)
            r["binned_fraction"] = float(r["binned_fraction"]) if r.get("binned_fraction") not in (None, "") else None
            sel_rows.append(r)

    fdir = Path(args.fasta_dir)
    fastas = {r["sample_id"]: fdir / f"{r['sample_id']}.fa.gz" for r in sel_rows}
    wdir = Path(args.work_dir)
    wdir.mkdir(parents=True, exist_ok=True)
    root = Path(args.root) if args.root else None

    res = analyze_pilot_contigs(
        sel_rows,
        fastas,
        wdir,
        root=root,
        duck_mem=args.duck_mem,
        threads=args.threads,
        reuse_chunks=args.reuse_chunks,
    )
    _write_tsv(wdir / "sample_contig_stats.tsv", res["per_sample"])
    _write_tsv(wdir / "contigs_by_habitat1.tsv", res["by_habitat"])
    _write_tsv(wdir / "contigs_by_habitat2.tsv", res["by_habitat2"])
    if res["repeats_by_habitat"]:
        _write_tsv(wdir / "repeats_by_habitat1.tsv", res["repeats_by_habitat"])

    summary = {
        "overall": res["Overall"],
        "by_habitat": res["by_habitat"],
        "repeats": res["repeats"],
        "repeats_ge_1500": res["repeats_ge_1500"],
        "repeats_by_habitat": res["repeats_by_habitat"],
        "root_overlap": res["root_overlap"],
    }
    (wdir / "contigs_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({
        "overall": res["Overall"],
        "repeats": res["repeats"],
        "repeats_ge_1500": res["repeats_ge_1500"],
        "root_overlap": res["root_overlap"],
    }))
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    run_p = sub.add_parser("run")
    run_p.add_argument("--selection", required=True)
    run_p.add_argument("--fasta-dir", required=True)
    run_p.add_argument("--work-dir", required=True)
    run_p.add_argument("--root", default=None)
    run_p.add_argument("--duck-mem", default="4GB")
    run_p.add_argument("--threads", type=int, default=1)
    run_p.add_argument("--reuse-chunks", action="store_true")
    args = ap.parse_args(argv)
    if args.cmd == "run":
        return cmd_run(args)
    return 2


if __name__ == "__main__":
    sys.exit(main())
