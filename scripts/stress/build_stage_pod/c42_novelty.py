"""Card C-42 step (5): MMseqs2 50%-identity linear clustering novelty against `r26.11.1` / `rM.2`.

    python -I c42_novelty.py export-faa --ref-root /tmp/if2453/r10/root --ref-release r26.11.1
        --pilot-root /tmp/if2453/s3/spire/root --pilot-stages SPH1,... --work-dir /tmp/if2453/s3/spire/novelty
    python -I c42_novelty.py linclust-50 --work-dir /tmp/if2453/s3/spire/novelty --threads 32 [--align-parts 4]
        [--split-mem 80G]
    python -I c42_novelty.py analyze --work-dir /tmp/if2453/s3/spire/novelty --pilot-root /tmp/if2453/s3/spire/root
        --inputs-dir /tmp/if2453/s3/spire/inputs --selection /tmp/if2453/s3/spire/select/pilot_selection.tsv

Bounded memory (AGENTS §4.1):
- `export-faa` k-way merges already-sorted `protein_seq` streams from the reference release and the
  pilot stages without holding the protein set in RAM, writing `combined.faa` and `ref_uids.parquet`.
- `linclust-50` runs the exact MMseqs2 50%-identity linear clustering pipeline (`kmermatcher` with
  `--split-memory-limit`, `rescorediagonal`, `clust`, `createsubdb`, `filterdb`, `rescorediagonal`,
  sharded `align` so every process stays < 280 CPU-min, `clust`, `mergeclusters`, `createtsv`).
- `analyze` uses spilled DuckDB joins + a compact bitset (`1 bit` per novel cluster) to compute exact
  per-habitat novel NR50 families per Gbp and sample-by-sample accumulation curves.
"""

from __future__ import annotations

import argparse
from array import array
import csv
import heapq
import json
import os
import random
import shlex
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

from biobox import catalog, io
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.l2 import cluster as l2cluster


def compute_novelty(
    clu_tsv: Path,
    ref_uids: Iterable[bytes],
    samples: list[dict],
    work_dir: Path,
    *,
    duck_mem: str = "4GB",
    n_perm: int = 10,
    seed: int = 20261118,
) -> dict:
    """Compute exact NR50 novelty and accumulation curves (unit-testable core).

    `samples`: list of `{"sample_id": str, "habitat": str, "bp": int, "proteins": Iterable[bytes]}` in
    deterministic pilot order.
    """
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    spill = work_dir / "spill"
    spill.mkdir(parents=True, exist_ok=True)

    # Convert clu_tsv (hex rep \t hex member) and ref_uids to Parquet for spilled DuckDB processing
    clu_pq = work_dir / "clu.parquet"
    ref_pq = work_dir / "ref.parquet"

    ref_list = sorted(set(ref_uids))
    pq.write_table(pa.table({"member_uid": pa.array(ref_list, type=pa.binary(16))}), ref_pq, compression="zstd")

    reps_b: list[bytes] = []
    mems_b: list[bytes] = []
    with open(clu_tsv) as fh:
        for line in fh:
            if not line.strip():
                continue
            r_hex, m_hex = line.rstrip("\n").split("\t")
            reps_b.append(bytes.fromhex(r_hex))
            mems_b.append(bytes.fromhex(m_hex))
    pq.write_table(
        pa.table({
            "rep_uid": pa.array(reps_b, type=pa.binary(16)),
            "member_uid": pa.array(mems_b, type=pa.binary(16)),
        }),
        clu_pq,
        compression="zstd",
    )

    # Write pilot distinct proteins
    pilot_set: set[bytes] = set()
    for s in samples:
        pilot_set.update(s["proteins"])
    pilot_pq = work_dir / "pilot.parquet"
    pq.write_table(
        pa.table({"member_uid": pa.array(sorted(pilot_set), type=pa.binary(16))}),
        pilot_pq,
        compression="zstd",
    )

    con = duckdb.connect(config={"threads": 1, "memory_limit": duck_mem, "temp_directory": str(spill)})
    try:
        summary_row = con.execute(
            f"""
            WITH m AS (
                SELECT
                    c.rep_uid,
                    c.member_uid,
                    r.member_uid IS NOT NULL AS in_ref,
                    p.member_uid IS NOT NULL AS in_pilot
                FROM read_parquet('{clu_pq}') c
                LEFT JOIN read_parquet('{ref_pq}') r USING (member_uid)
                LEFT JOIN read_parquet('{pilot_pq}') p USING (member_uid)
            ),
            per_clu AS (
                SELECT
                    rep_uid,
                    count(*)::BIGINT AS sz,
                    bool_or(in_ref) AS has_ref,
                    bool_or(in_pilot) AS has_pilot
                FROM m
                GROUP BY rep_uid
            )
            SELECT
                count(*)::BIGINT AS total_clusters,
                coalesce(sum(CASE WHEN has_ref THEN 1 ELSE 0 END), 0)::BIGINT AS ref_clusters,
                coalesce(sum(CASE WHEN has_ref AND has_pilot THEN 1 ELSE 0 END), 0)::BIGINT AS shared_clusters,
                coalesce(sum(CASE WHEN NOT has_ref AND has_pilot THEN 1 ELSE 0 END), 0)::BIGINT AS novel_clusters,
                coalesce(sum(CASE WHEN NOT has_ref AND has_pilot AND sz >= 2 THEN 1 ELSE 0 END), 0)::BIGINT AS novel_nonsingleton
            FROM per_clu
            """
        ).fetchone()

        # Assign compact 1..K IDs to novel clusters and map each pilot member_uid -> (cid, is_nonsingleton)
        novel_rows = con.execute(
            f"""
            WITH m AS (
                SELECT
                    c.rep_uid,
                    c.member_uid,
                    r.member_uid IS NOT NULL AS in_ref,
                    p.member_uid IS NOT NULL AS in_pilot
                FROM read_parquet('{clu_pq}') c
                LEFT JOIN read_parquet('{ref_pq}') r USING (member_uid)
                LEFT JOIN read_parquet('{pilot_pq}') p USING (member_uid)
            ),
            novel_clu AS (
                SELECT
                    rep_uid,
                    row_number() OVER (ORDER BY rep_uid)::INTEGER AS cid,
                    (count(*) >= 2) AS nonsingle
                FROM m
                GROUP BY rep_uid
                HAVING NOT bool_or(in_ref) AND bool_or(in_pilot)
            )
            SELECT m.member_uid, n.cid, n.nonsingle
            FROM m
            JOIN novel_clu n USING (rep_uid)
            WHERE m.in_pilot
            """
        ).fetchall()
    finally:
        con.close()
        shutil.rmtree(spill, ignore_errors=True)

    k_novel = int(summary_row[3])
    prot_to_cid: dict[bytes, tuple[int, bool]] = {u: (int(cid), bool(ns)) for u, cid, ns in novel_rows}

    def _run_curve(ordered_samples: list[dict]) -> list[dict]:
        seen = bytearray((k_novel + 8) // 8)
        cum_bp = cum_new = cum_ns = 0
        curve = []
        for idx, s in enumerate(ordered_samples, 1):
            new_c = new_ns = 0
            for pu in s["proteins"]:
                info = prot_to_cid.get(pu)
                if info is None:
                    continue
                cid, is_ns = info
                byte_i, mask = cid >> 3, 1 << (cid & 7)
                if not (seen[byte_i] & mask):
                    seen[byte_i] |= mask
                    new_c += 1
                    if is_ns:
                        new_ns += 1
            bp = int(s["bp"])
            cum_bp += bp
            cum_new += new_c
            cum_ns += new_ns
            curve.append({
                "step": idx,
                "sample_id": s["sample_id"],
                "habitat": s["habitat"],
                "sample_bp": bp,
                "cum_bp": cum_bp,
                "new_clusters": new_c,
                "cum_new_clusters": cum_new,
                "new_nonsingleton_clusters": new_ns,
                "cum_new_nonsingleton_clusters": cum_ns,
                "marginal_new_per_gbp": round(new_c / (bp / 1e9), 2) if bp > 0 else 0.0,
                "cum_new_per_gbp": round(cum_new / (cum_bp / 1e9), 2) if cum_bp > 0 else 0.0,
            })
        return curve

    global_curve = _run_curve(samples)

    by_hab_samples: dict[str, list[dict]] = defaultdict(list)
    for s in samples:
        by_hab_samples[s["habitat"]].append(s)

    by_habitat: list[dict] = []
    habitat_curves: dict[str, list[dict]] = {}
    rng = random.Random(seed)
    for hab in sorted(by_hab_samples):
        hs = by_hab_samples[hab]
        h_curve = _run_curve(hs)
        habitat_curves[hab] = h_curve
        last = h_curve[-1]
        gbp = last["cum_bp"] / 1e9 if last["cum_bp"] > 0 else 0.0
        # First-half vs second-half rate in permutation average to measure saturation slope cleanly
        first_rates, last_rates = [], []
        if len(hs) >= 2 and gbp > 0:
            half = max(1, len(hs) // 2)
            for _ in range(n_perm):
                perm = list(hs)
                rng.shuffle(perm)
                pc = _run_curve(perm)
                bp1 = pc[half - 1]["cum_bp"] / 1e9
                bp2 = (pc[-1]["cum_bp"] - pc[half - 1]["cum_bp"]) / 1e9
                if bp1 > 0:
                    first_rates.append(pc[half - 1]["cum_new_clusters"] / bp1)
                if bp2 > 0:
                    last_rates.append((pc[-1]["cum_new_clusters"] - pc[half - 1]["cum_new_clusters"]) / bp2)
        by_habitat.append({
            "habitat": hab,
            "samples": len(hs),
            "total_bp": last["cum_bp"],
            "novel_clusters_in_habitat": last["cum_new_clusters"],
            "novel_nonsingleton_in_habitat": last["cum_new_nonsingleton_clusters"],
            "novel_per_gbp": round(last["cum_new_clusters"] / gbp, 2) if gbp > 0 else 0.0,
            "novel_nonsingleton_per_gbp": round(last["cum_new_nonsingleton_clusters"] / gbp, 2) if gbp > 0 else 0.0,
            "first_half_mean_novel_per_gbp": round(sum(first_rates) / len(first_rates), 2) if first_rates else None,
            "second_half_mean_novel_per_gbp": round(sum(last_rates) / len(last_rates), 2) if last_rates else None,
        })

    return {
        "summary": {
            "total_clusters": int(summary_row[0]),
            "ref_clusters": int(summary_row[1]),
            "shared_clusters": int(summary_row[2]),
            "novel_clusters": k_novel,
            "novel_nonsingleton_clusters": int(summary_row[4]),
        },
        "by_habitat": by_habitat,
        "global_curve": global_curve,
        "habitat_curves": habitat_curves,
    }


# --- Large-scale pod subcommands --------------------------------------------------------------


def _stage_root(base: Path, stage: str) -> Path:
    cand = base / stage
    return cand if (cand / "l0seq").exists() or (cand / "runs").exists() else base


def cmd_export_faa(args: argparse.Namespace) -> int:
    ref_root = Path(args.ref_root)
    pilot_root = Path(args.pilot_root)
    wdir = Path(args.work_dir)
    wdir.mkdir(parents=True, exist_ok=True)

    # Pick reference stages from release manifest if present, else explicit --ref-stages
    if args.ref_stages:
        ref_stages = [s.strip() for s in args.ref_stages.split(",") if s.strip()]
        ref_label = args.ref_release or ",".join(ref_stages)
    else:
        man_path = ref_root / "manifests" / f"release={args.ref_release}.json"
        if man_path.is_file():
            ref_stages = catalog.release_stages(ref_root, args.ref_release)
            ref_label = args.ref_release
        else:
            fallback = "rM.2"
            ref_stages = catalog.release_stages(ref_root, fallback)
            ref_label = fallback

    pilot_stages = [s.strip() for s in args.pilot_stages.split(",") if s.strip()]
    ref_files = [l2cluster._protein_files(_stage_root(ref_root, s), GENE_CALL_VERSION, [s])[0] for s in ref_stages]
    pilot_files = [l2cluster._protein_files(_stage_root(pilot_root, s), GENE_CALL_VERSION, [s])[0] for s in pilot_stages]

    # 1. Stream distinct reference protein_uids to ref_uids.parquet
    ref_pq = wdir / "ref_uids.parquet"
    ref_writer = pq.ParquetWriter(ref_pq, pa.schema([("member_uid", pa.binary(16))]), compression="zstd")
    buf: list[bytes] = []
    n_ref = 0
    last = None
    for u, _ in heapq.merge(*(l2cluster._stage_proteins(f, seq=False) for f in ref_files), key=lambda r: r[0]):
        if u == last:
            continue
        last = u
        n_ref += 1
        buf.append(u)
        if len(buf) >= 250_000:
            ref_writer.write_table(pa.table({"member_uid": pa.array(buf, type=pa.binary(16))}))
            buf.clear()
    if buf:
        ref_writer.write_table(pa.table({"member_uid": pa.array(buf, type=pa.binary(16))}))
        buf.clear()
    ref_writer.close()

    # 2. Stream distinct pilot protein_uids to pilot_uids.parquet
    pilot_pq = wdir / "pilot_uids.parquet"
    pilot_writer = pq.ParquetWriter(pilot_pq, pa.schema([("member_uid", pa.binary(16))]), compression="zstd")
    n_pilot = 0
    last = None
    for u, _ in heapq.merge(*(l2cluster._stage_proteins(f, seq=False) for f in pilot_files), key=lambda r: r[0]):
        if u == last:
            continue
        last = u
        n_pilot += 1
        buf.append(u)
        if len(buf) >= 250_000:
            pilot_writer.write_table(pa.table({"member_uid": pa.array(buf, type=pa.binary(16))}))
            buf.clear()
    if buf:
        pilot_writer.write_table(pa.table({"member_uid": pa.array(buf, type=pa.binary(16))}))
        buf.clear()
    pilot_writer.close()

    # 3. K-way merge all reference + pilot proteins into combined.faa (or directly into mmseqs createdb stdin)
    faa_path = wdir / "combined.faa"
    n_combined = 0
    faa_bytes = 0
    last = None
    streams = [*(l2cluster._stage_proteins(f, seq=True) for f in ref_files),
               *(l2cluster._stage_proteins(f, seq=True) for f in pilot_files)]
    if args.direct_createdb:
        for rel in ("db", "db_h"):
            l2cluster._clear(wdir, rel)
        done_file = wdir / "steps" / "00_createdb.json"
        done_file.parent.mkdir(parents=True, exist_ok=True)
        done_file.unlink(missing_ok=True)
        argv = ["mmseqs", "createdb", "stdin", "db", "--dbtype", "1", "--shuffle", "0", "--write-lookup", "0", "-v", "2"]
        t0 = time.perf_counter()
        with open(wdir / "mmseqs.log", "a") as log:
            log.write(f"=== 00_createdb (stdin): {shlex.join(argv)} ===\n")
            log.flush()
            proc = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=log, stderr=subprocess.STDOUT, cwd=wdir)
            assert proc.stdin is not None
            chunk: list[str] = []
            try:
                for u, seq in heapq.merge(*streams, key=lambda r: r[0]):
                    if u == last:
                        continue
                    last = u
                    rec_s = f">{u.hex()}\n{seq}\n"
                    faa_bytes += len(rec_s)
                    chunk.append(rec_s)
                    n_combined += 1
                    if len(chunk) >= 16_384:
                        proc.stdin.write("".join(chunk).encode())
                        chunk.clear()
                if chunk:
                    proc.stdin.write("".join(chunk).encode())
                    chunk.clear()
            finally:
                proc.stdin.close()
            _, status, ru = os.wait4(proc.pid, 0)
        wall = time.perf_counter() - t0
        code = os.waitstatus_to_exitcode(status)
        if code != 0:
            raise RuntimeError(f"createdb stdin failed with exit code {code}")
        done_file.write_text(json.dumps({
            "label": "00_createdb",
            "argv": argv,
            "wall_s": round(wall, 3),
            "cpu_min": round((ru.ru_utime + ru.ru_stime) / 60.0, 3),
            "max_rss_mb": round(ru.ru_maxrss / 1024.0, 1),
        }, indent=2) + "\n")
    else:
        tmp_faa = faa_path.with_name(".combined.faa.tmp")
        with open(tmp_faa, "w", buffering=1 << 22) as fh:
            for u, seq in heapq.merge(*streams, key=lambda r: r[0]):
                if u == last:
                    continue
                last = u
                fh.write(f">{u.hex()}\n{seq}\n")
                n_combined += 1
        tmp_faa.replace(faa_path)
        faa_bytes = faa_path.stat().st_size

    meta = {
        "ref_release": ref_label,
        "ref_stages": ref_stages,
        "pilot_stages": pilot_stages,
        "ref_distinct_proteins": n_ref,
        "pilot_distinct_proteins": n_pilot,
        "exact_shared_proteins": n_ref + n_pilot - n_combined,
        "combined_distinct_proteins": n_combined,
        "combined_faa_bytes": faa_bytes,
        "direct_createdb": bool(args.direct_createdb),
    }
    (wdir / "export_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta))
    return 0


def _run_mmseqs(wdir: Path, label: str, argv: list[str], stop_files: list[Path]) -> dict:
    for sf in stop_files:
        if sf.exists():
            raise SystemExit(f"STOP file exists: {sf}")
    done_file = wdir / "steps" / f"{label}.json"
    if done_file.is_file():
        return json.loads(done_file.read_text())
    done_file.parent.mkdir(parents=True, exist_ok=True)
    log_path = wdir / "mmseqs.log"
    with open(log_path, "a") as log:
        log.write(f"=== {label}: {shlex.join(argv)} ===\n")
        log.flush()
        t0 = time.perf_counter()
        p = subprocess.Popen(argv, stdout=log, stderr=subprocess.STDOUT, cwd=wdir)
        _, status, ru = os.wait4(p.pid, 0)
        wall = time.perf_counter() - t0
    code = os.waitstatus_to_exitcode(status)
    if code != 0:
        raise RuntimeError(f"{label} failed with exit code {code}: {shlex.join(argv)}")
    rec = {
        "label": label,
        "argv": argv,
        "wall_s": round(wall, 3),
        "cpu_min": round((ru.ru_utime + ru.ru_stime) / 60.0, 3),
        "max_rss_mb": round(ru.ru_maxrss / 1024.0, 1),
    }
    done_file.write_text(json.dumps(rec, indent=2) + "\n")
    return rec


def cmd_linclust_50(args: argparse.Namespace) -> int:
    wdir = Path(args.work_dir).resolve()
    lin = wdir / "lin50"
    lin.mkdir(parents=True, exist_ok=True)
    t = str(args.threads)
    stops = [Path("/tmp/if2453/r10/root/STOP"), Path("/tmp/if2453/s3/spire/STOP")]
    records = []

    records.append(_run_mmseqs(
        wdir, "00_createdb",
        ["mmseqs", "createdb", "combined.faa", "db", "--dbtype", "1", "--shuffle", "0", "--write-lookup", "0", "-v", "2"],
        stops,
    ))
    if args.cleanup_db and (wdir / "db.dbtype").is_file() and (wdir / "combined.faa").is_file():
        (wdir / "combined.faa").unlink()

    km_args = list(l2cluster._kmermatcher("db", "lin50/pref", "0.5"))
    idx_sm = km_args.index("--split-memory-limit")
    km_args[idx_sm + 1] = args.split_mem
    records.append(_run_mmseqs(wdir, "01_kmermatcher", ["mmseqs", *km_args, "--threads", t, "-v", "2"], stops))
    records.append(_run_mmseqs(
        wdir, "02_rescore1",
        ["mmseqs", *l2cluster._rescore("db", "lin50/pref", "lin50/pref_rescore1", "0", "0", "0.5"), "--threads", t, "-v", "2"],
        stops,
    ))
    records.append(_run_mmseqs(
        wdir, "03_pre_clust",
        ["mmseqs", *l2cluster._clust("db", "lin50/pref_rescore1", "lin50/pre_clust"), "--threads", t, "-v", "2"],
        stops,
    ))
    if args.cleanup_db and (lin / "pre_clust.dbtype").is_file():
        l2cluster._clear(wdir, "lin50/pref_rescore1")

    # Write order_redundancy = first column of lin50/pre_clust.index (as MMseqs2's linclust.sh does)
    order_file = lin / "order_redundancy"
    if not order_file.is_file():
        with open(lin / "pre_clust.index", "rb") as r_fh, open(order_file, "wb") as w_fh:
            for ln in r_fh:
                w_fh.write(ln.split(b"\t", 1)[0] + b"\n")

    records.append(_run_mmseqs(
        wdir, "04_subdb_seq",
        ["mmseqs", *l2cluster._subdb("lin50/order_redundancy", "db", "lin50/input_step_redundancy"), "-v", "2"],
        stops,
    ))
    records.append(_run_mmseqs(
        wdir, "05_subdb_pref",
        ["mmseqs", *l2cluster._subdb("lin50/order_redundancy", "lin50/pref", "lin50/pref_filter1"), "-v", "2"],
        stops,
    ))

    records.append(_run_mmseqs(
        wdir, "06_filterdb",
        ["mmseqs", "filterdb", "lin50/pref_filter1", "lin50/pref_filter2", "--filter-file", "lin50/order_redundancy",
         "--compressed", "0", "--threads", t, "-v", "2"],
        stops,
    ))
    if args.cleanup_db and (lin / "pref_filter2.dbtype").is_file():
        # Note: `pref_filter1` uses `--subdb-mode 1` (soft-links `pref`), so clear `pref` only AFTER `06_filterdb`!
        l2cluster._clear(wdir, "lin50/pref_filter1")
        l2cluster._clear(wdir, "lin50/pref")

    records.append(_run_mmseqs(
        wdir, "07_rescore2",
        ["mmseqs", *l2cluster._rescore("lin50/input_step_redundancy", "lin50/pref_filter2", "lin50/pref_rescore2",
                                       "1", "1", "0.5"), "--threads", t, "-v", "2"],
        stops,
    ))
    if args.cleanup_db and (lin / "pref_rescore2.dbtype").is_file():
        l2cluster._clear(wdir, "lin50/pref_filter2")

    join_done = wdir / "steps" / "08_align_join.json"
    if args.align_parts <= 1:
        records.append(_run_mmseqs(
            wdir, "08_align",
            ["mmseqs", *l2cluster._align("lin50/input_step_redundancy", "lin50/pref_rescore2", "lin50/aln", "0.5", "1"),
             "--threads", t, "-v", "2"],
            stops,
        ))
    elif join_done.is_file():
        for pf in sorted((wdir / "steps").glob("08_align_part_*.json")):
            records.append(json.loads(pf.read_text()))
        records.append(json.loads(join_done.read_text()))
    else:
        cut_parts = [
            cp for cp in l2cluster.cut(
                lambda: ((k, n) for k, _, n in l2cluster.index_rows(lin / "pref_rescore2")),
                1.0,
                parts=args.align_parts,
            )
            if cp["keys"] > 0
        ]
        if len(cut_parts) <= 1:
            records.append(_run_mmseqs(
                wdir, "08_align",
                ["mmseqs", *l2cluster._align("lin50/input_step_redundancy", "lin50/pref_rescore2", "lin50/aln", "0.5", "1"),
                 "--threads", t, "-v", "2"],
                stops,
            ))
        else:
            parts_dir = lin / "aln_parts"
            parts_dir.mkdir(parents=True, exist_ok=True)
            out_paths: list[Path] = []
            for k, cp in enumerate(cut_parts):
                part_done = wdir / "steps" / f"08_align_part_{k:04d}.json"
                rel_ids = f"lin50/aln_parts/ids_{k:04d}"
                rel_in = f"lin50/aln_parts/in_{k:04d}"
                rel_out = f"lin50/aln_parts/out_{k:04d}"
                if not part_done.is_file():
                    with open(wdir / rel_ids, "w") as id_fh:
                        for key, _, _ in l2cluster.index_rows(lin / "pref_rescore2"):
                            if cp["first"] <= key <= cp["last"]:
                                id_fh.write(f"{key}\n")
                    l2cluster._clear(wdir, rel_in)
                    l2cluster._clear(wdir, rel_out)
                    _run_mmseqs(
                        wdir, f"08_subdb_part_{k:04d}",
                        ["mmseqs", "createsubdb", rel_ids, "lin50/pref_rescore2", rel_in, "--subdb-mode", "0", "-v", "2"],
                        stops,
                    )
                records.append(_run_mmseqs(
                    wdir, f"08_align_part_{k:04d}",
                    ["mmseqs", *l2cluster._align("lin50/input_step_redundancy", rel_in, rel_out, "0.5", "1"),
                     "--threads", t, "-v", "2"],
                    stops,
                ))
                if args.cleanup_db:
                    l2cluster._clear(wdir, rel_in)
                    (wdir / rel_ids).unlink(missing_ok=True)
                out_paths.append(wdir / rel_out)

            t_j0 = time.perf_counter()
            l2cluster.join_dbs(out_paths, lin / "aln")
            join_rec = {
                "label": "08_align_join",
                "argv": ["biobox.l2.cluster.join_dbs", f"parts={len(out_paths)}"],
                "wall_s": round(time.perf_counter() - t_j0, 3),
                "cpu_min": 0.0,
                "max_rss_mb": 0.0,
            }
            join_done.write_text(json.dumps(join_rec, indent=2) + "\n")
            records.append(join_rec)

    if args.cleanup_db and (lin / "aln.dbtype").is_file():
        shutil.rmtree(lin / "aln_parts", ignore_errors=True)
        l2cluster._clear(wdir, "lin50/pref_rescore2")

    records.append(_run_mmseqs(
        wdir, "09_clust",
        ["mmseqs", *l2cluster._clust("lin50/input_step_redundancy", "lin50/aln", "lin50/clust"), "--threads", t, "-v", "2"],
        stops,
    ))
    if args.cleanup_db and (lin / "clust.dbtype").is_file():
        l2cluster._clear(wdir, "lin50/aln")
        l2cluster._clear(wdir, "lin50/input_step_redundancy")

    records.append(_run_mmseqs(
        wdir, "10_mergeclusters",
        ["mmseqs", "mergeclusters", "db", "clu50", "lin50/pre_clust", "lin50/clust", "--compressed", "0",
         "--threads", t, "-v", "2"],
        stops,
    ))
    records.append(_run_mmseqs(
        wdir, "11_createtsv",
        ["mmseqs", "createtsv", "db", "db", "clu50", "clu50.tsv", "--threads", t, "-v", "2"],
        stops,
    ))

    # Clean intermediate MMseqs2 DBs + combined.faa once clu50.tsv is written to save scratch space!
    if (wdir / "clu50.tsv").is_file() and args.cleanup_db:
        shutil.rmtree(lin, ignore_errors=True)
        for rel in ("db", "db_h", "clu50"):
            l2cluster._clear(wdir, rel)
        if (wdir / "combined.faa").is_file():
            (wdir / "combined.faa").unlink()

    (wdir / "linclust_steps.json").write_text(json.dumps(records, indent=2) + "\n")
    print(json.dumps({"steps": len(records), "total_cpu_min": round(sum(r["cpu_min"] for r in records), 2),
                      "peak_rss_mb": max(r["max_rss_mb"] for r in records)}))
    return 0


def cmd_analyze(args: argparse.Namespace) -> int:
    wdir = Path(args.work_dir).resolve()
    pilot_root = Path(args.pilot_root).resolve()
    inputs_dir = Path(args.inputs_dir).resolve()
    spill = wdir / "analyze_spill"
    spill.mkdir(parents=True, exist_ok=True)

    # 1. Convert clu50.tsv to Parquet in blocks (bounded memory)
    clu_pq = wdir / "clu50.parquet"
    if not clu_pq.is_file():
        tmp_pq = clu_pq.with_name(".clu50.parquet.tmp")
        writer = pq.ParquetWriter(
            tmp_pq,
            pa.schema([("rep_uid", pa.binary(16)), ("member_uid", pa.binary(16))]),
            compression="zstd",
        )
        reader = pacsv.open_csv(
            wdir / "clu50.tsv",
            read_options=pacsv.ReadOptions(column_names=["rep", "mem"], block_size=1 << 24),
            parse_options=pacsv.ParseOptions(delimiter="\t"),
            convert_options=pacsv.ConvertOptions(column_types={"rep": pa.string(), "mem": pa.string()}),
        )
        for batch in reader:
            reps = [bytes.fromhex(x) for x in batch.column(0).to_pylist()]
            mems = [bytes.fromhex(x) for x in batch.column(1).to_pylist()]
            writer.write_table(pa.table({
                "rep_uid": pa.array(reps, type=pa.binary(16)),
                "member_uid": pa.array(mems, type=pa.binary(16)),
            }))
        writer.close()
        tmp_pq.replace(clu_pq)
    if args.cleanup_tsv and clu_pq.is_file() and (wdir / "clu50.tsv").is_file():
        (wdir / "clu50.tsv").unlink()

    ref_pq = wdir / "ref_uids.parquet"
    pilot_pq = wdir / "pilot_uids.parquet"

    con = duckdb.connect(config={"threads": args.threads, "memory_limit": args.duck_mem, "temp_directory": str(spill)})
    try:
        summary_row = con.execute(
            f"""
            WITH m AS (
                SELECT
                    c.rep_uid,
                    c.member_uid,
                    r.member_uid IS NOT NULL AS in_ref,
                    p.member_uid IS NOT NULL AS in_pilot
                FROM read_parquet('{clu_pq}') c
                LEFT JOIN read_parquet('{ref_pq}') r USING (member_uid)
                LEFT JOIN read_parquet('{pilot_pq}') p USING (member_uid)
            ),
            per_clu AS (
                SELECT
                    rep_uid,
                    count(*)::BIGINT AS sz,
                    bool_or(in_ref) AS has_ref,
                    bool_or(in_pilot) AS has_pilot
                FROM m
                GROUP BY rep_uid
            )
            SELECT
                count(*)::BIGINT AS total_clusters,
                coalesce(sum(CASE WHEN has_ref THEN 1 ELSE 0 END), 0)::BIGINT AS ref_clusters,
                coalesce(sum(CASE WHEN has_ref AND has_pilot THEN 1 ELSE 0 END), 0)::BIGINT AS shared_clusters,
                coalesce(sum(CASE WHEN NOT has_ref AND has_pilot THEN 1 ELSE 0 END), 0)::BIGINT AS novel_clusters,
                coalesce(sum(CASE WHEN NOT has_ref AND has_pilot AND sz >= 2 THEN 1 ELSE 0 END), 0)::BIGINT AS novel_nonsingleton
            FROM per_clu
            """
        ).fetchone()

        # Write mapping of pilot member_uid -> (cid, nonsingleton) to Parquet on disk
        map_pq = wdir / "pilot_novel_map.parquet"
        con.execute(
            f"""
            COPY (
                WITH m AS (
                    SELECT
                        c.rep_uid,
                        c.member_uid,
                        r.member_uid IS NOT NULL AS in_ref,
                        p.member_uid IS NOT NULL AS in_pilot
                    FROM read_parquet('{clu_pq}') c
                    LEFT JOIN read_parquet('{ref_pq}') r USING (member_uid)
                    LEFT JOIN read_parquet('{pilot_pq}') p USING (member_uid)
                ),
                novel_clu AS (
                    SELECT
                        rep_uid,
                        row_number() OVER (ORDER BY rep_uid)::INTEGER AS cid,
                        (count(*) >= 2) AS nonsingle
                    FROM m
                    GROUP BY rep_uid
                    HAVING NOT bool_or(in_ref) AND bool_or(in_pilot)
                )
                SELECT m.member_uid AS protein_uid, n.cid, n.nonsingle
                FROM m
                JOIN novel_clu n USING (rep_uid)
                WHERE m.in_pilot
            ) TO '{map_pq}' (FORMAT PARQUET, COMPRESSION ZSTD)
            """
        )

        # Map every sample -> distinct novel (cid, nonsingle) via contig_occurrence + cross_stage_occurrences + protein_occurrence
        stages_man = json.loads((inputs_dir / "stages_manifest.json").read_text())
        pilot_stages = [s["stage"] for s in stages_man]
        c_occ_files = []
        p_occ_files = []
        l0_man_files = []
        cross_rows: list[tuple[str, bytes]] = []
        for st in pilot_stages:
            sr = _stage_root(pilot_root, st)
            c_occ_files.extend(io.part_files(sr / "l0seq" / "contig_occurrence" / f"stage={st}", "source=*/part-*.parquet"))
            p_occ_files.extend(sorted((sr / "l1" / "protein_occurrence").glob(f"gene_call_version=*/stage={st}/part-*.parquet")))
            l0_man_files.extend(io.part_files(sr / "l0" / "l0_manifest" / f"stage={st}"))
            ctsv = inputs_dir / st / "cross_stage_occurrences.tsv"
            if ctsv.is_file():
                with open(ctsv, newline="") as fh:
                    for r in csv.DictReader(fh, delimiter="\t"):
                        cross_rows.append((r["sample_id"], bytes.fromhex(r["contig_uid"])))

        cross_pq = wdir / "cross_stage_occ.parquet"
        pq.write_table(
            pa.table({
                "sample_id": pa.array([x[0] for x in cross_rows], type=pa.string()),
                "contig_uid": pa.array([x[1] for x in cross_rows], type=pa.binary(16)),
            }),
            cross_pq,
            compression="zstd",
        )

        c_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in c_occ_files)
        p_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in p_occ_files)
        m_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in l0_man_files)

        ns_cid_rows = con.execute(
            f"SELECT DISTINCT cid FROM read_parquet('{map_pq}') WHERE nonsingle"
        ).fetchall()

        sample_cids_pq = wdir / "sample_novel_cids.parquet"
        con.execute(
            f"""
            COPY (
                WITH obj_sample AS (
                    SELECT DISTINCT
                        object_sha256 AS l0_object_sha256,
                        regexp_replace(regexp_extract(object_path, '([^/]+)\\.fa(\\.gz)?$', 1), '\\.fa(\\.gz)?$', '') AS sample_id
                    FROM read_parquet([{m_lits}], hive_partitioning=false)
                ),
                sample_contigs AS (
                    SELECT o.sample_id, c.contig_uid
                    FROM read_parquet([{c_lits}], hive_partitioning=false) c
                    JOIN obj_sample o USING (l0_object_sha256)
                    UNION
                    SELECT sample_id, contig_uid FROM read_parquet('{cross_pq}')
                ),
                contig_prots AS (
                    SELECT DISTINCT contig_uid, protein_uid
                    FROM read_parquet([{p_lits}], hive_partitioning=false)
                )
                SELECT DISTINCT sc.sample_id, nm.cid
                FROM sample_contigs sc
                JOIN contig_prots cp USING (contig_uid)
                JOIN read_parquet('{map_pq}') nm USING (protein_uid)
                ORDER BY sc.sample_id, nm.cid
            ) TO '{sample_cids_pq}' (FORMAT PARQUET, COMPRESSION ZSTD)
            """
        )

        # Also count sample bp actually kept in the scratch root (sum of contig lengths in each sample)
        ci_files = []
        for st in pilot_stages:
            sr = _stage_root(pilot_root, st)
            ci_files.extend(io.part_files(sr / "l0seq" / "contig_index" / f"stage={st}"))
        ci_lits = ", ".join("'" + str(f).replace("'", "''") + "'" for f in ci_files)
        s_bp_rows = con.execute(
            f"""
            WITH obj_sample AS (
                SELECT DISTINCT
                    object_sha256 AS l0_object_sha256,
                    regexp_replace(regexp_extract(object_path, '([^/]+)\\.fa(\\.gz)?$', 1), '\\.fa(\\.gz)?$', '') AS sample_id
                FROM read_parquet([{m_lits}], hive_partitioning=false)
            ),
            sample_contigs AS (
                SELECT o.sample_id, c.contig_uid
                FROM read_parquet([{c_lits}], hive_partitioning=false) c
                JOIN obj_sample o USING (l0_object_sha256)
                UNION
                SELECT sample_id, contig_uid FROM read_parquet('{cross_pq}')
            ),
            ci AS (
                SELECT contig_uid, max(contig_len)::BIGINT AS len
                FROM read_parquet([{ci_lits}], hive_partitioning=false)
                WHERE is_default_variant
                GROUP BY contig_uid
            )
            SELECT sc.sample_id, coalesce(sum(ci.len), 0)::BIGINT
            FROM sample_contigs sc
            JOIN ci USING (contig_uid)
            GROUP BY sc.sample_id
            """
        ).fetchall()
    finally:
        con.close()
        shutil.rmtree(spill, ignore_errors=True)

    k_novel = int(summary_row[3])
    ns_bits = bytearray((k_novel + 8) // 8)
    for (cid_val,) in ns_cid_rows:
        c_i = int(cid_val)
        ns_bits[c_i >> 3] |= 1 << (c_i & 7)
    del ns_cid_rows

    bp_by_sample = {sid: int(bp) for sid, bp in s_bp_rows}
    # Load per-sample packed uint32 arrays (`4 bytes` per cid instead of Python object pointers)
    cids_by_sample: dict[str, array] = defaultdict(lambda: array("I"))
    pf = pq.ParquetFile(sample_cids_pq)
    for rg in range(pf.num_row_groups):
        t = pf.read_row_group(rg)
        sids_col = t.column(0).to_pylist()
        cids_col = t.column(1).to_pylist()
        n_r = len(sids_col)
        i = 0
        while i < n_r:
            sid = sids_col[i]
            j = i + 1
            while j < n_r and sids_col[j] == sid:
                j += 1
            cids_by_sample[sid].extend(cids_col[i:j])
            i = j

    sel_rows = []
    with open(args.selection, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["sample_id"] in bp_by_sample:
                sel_rows.append(r)
    sel_rows.sort(key=lambda r: (int(r["rank_in_stratum"]), r["habitat"], r["habitat2"], r["sample_id"]))

    empty_arr = array("I")

    def run_curve(ordered_rows: list[dict]) -> list[dict]:
        seen = bytearray((k_novel + 8) // 8)
        cum_bp = cum_new = cum_ns = 0
        curve = []
        for step_i, r in enumerate(ordered_rows, 1):
            sid = r["sample_id"]
            cids = cids_by_sample.get(sid, empty_arr)
            new_c = new_ns = 0
            for cid in cids:
                bi, mask = cid >> 3, 1 << (cid & 7)
                if not (seen[bi] & mask):
                    seen[bi] |= mask
                    new_c += 1
                    if ns_bits[bi] & mask:
                        new_ns += 1
            bp = bp_by_sample.get(sid, 0)
            cum_bp += bp
            cum_new += new_c
            cum_ns += new_ns
            curve.append({
                "step": step_i,
                "sample_id": sid,
                "habitat": r["habitat"],
                "habitat2": r["habitat2"],
                "sample_bp": bp,
                "cum_bp": cum_bp,
                "new_clusters": new_c,
                "cum_new_clusters": cum_new,
                "new_nonsingleton_clusters": new_ns,
                "cum_new_nonsingleton_clusters": cum_ns,
                "marginal_new_per_gbp": round(new_c / (bp / 1e9), 2) if bp > 0 else 0.0,
                "cum_new_per_gbp": round(cum_new / (cum_bp / 1e9), 2) if cum_bp > 0 else 0.0,
            })
        return curve

    global_curve = run_curve(sel_rows)

    def summarize_strata(key: str) -> tuple[list[dict], list[dict]]:
        groups: dict[str, list[dict]] = defaultdict(list)
        for r in sel_rows:
            groups[r[key]].append(r)
        rng = random.Random(20261118)
        rows_out = []
        all_curves = []
        for g in sorted(groups):
            gs = sorted(groups[g], key=lambda x: (int(x["rank_in_stratum"]), x["sample_id"]))
            c = run_curve(gs)
            all_curves.extend([{key: g, **pt} for pt in c])
            last = c[-1]
            gbp = last["cum_bp"] / 1e9 if last["cum_bp"] > 0 else 0.0
            first_rates, last_rates = [], []
            if len(gs) >= 2 and gbp > 0:
                half = max(1, len(gs) // 2)
                for _ in range(20):
                    perm = list(gs)
                    rng.shuffle(perm)
                    pc = run_curve(perm)
                    bp1 = pc[half - 1]["cum_bp"] / 1e9
                    bp2 = (pc[-1]["cum_bp"] - pc[half - 1]["cum_bp"]) / 1e9
                    if bp1 > 0:
                        first_rates.append(pc[half - 1]["cum_new_clusters"] / bp1)
                    if bp2 > 0:
                        last_rates.append((pc[-1]["cum_new_clusters"] - pc[half - 1]["cum_new_clusters"]) / bp2)
            rows_out.append({
                key: g,
                "samples": len(gs),
                "total_bp": last["cum_bp"],
                "novel_clusters": last["cum_new_clusters"],
                "novel_nonsingleton_clusters": last["cum_new_nonsingleton_clusters"],
                "novel_per_gbp": round(last["cum_new_clusters"] / gbp, 2) if gbp > 0 else 0.0,
                "novel_nonsingleton_per_gbp": round(last["cum_new_nonsingleton_clusters"] / gbp, 2) if gbp > 0 else 0.0,
                "first_half_mean_novel_per_gbp": round(sum(first_rates) / len(first_rates), 2) if first_rates else None,
                "second_half_mean_novel_per_gbp": round(sum(last_rates) / len(last_rates), 2) if last_rates else None,
            })
        return rows_out, all_curves

    hab1_nov, hab1_curves = summarize_strata("habitat")
    hab2_nov, _ = summarize_strata("habitat2")

    def write_tsv(p: Path, rs: list[dict]) -> None:
        if not rs:
            return
        cols = list(rs[0].keys())
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            for r in rs:
                w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in cols})

    write_tsv(wdir / "novelty_by_habitat1.tsv", hab1_nov)
    write_tsv(wdir / "novelty_by_habitat2.tsv", hab2_nov)
    write_tsv(wdir / "accumulation_global.tsv", global_curve)
    write_tsv(wdir / "accumulation_by_habitat1.tsv", hab1_curves)

    tot_bp = global_curve[-1]["cum_bp"] if global_curve else 0
    tot_gbp = tot_bp / 1e9 if tot_bp else 0.0
    summary = {
        "total_clusters": int(summary_row[0]),
        "ref_clusters": int(summary_row[1]),
        "shared_clusters": int(summary_row[2]),
        "novel_clusters": k_novel,
        "novel_nonsingleton_clusters": int(summary_row[4]),
        "pilot_samples_in_root": len(sel_rows),
        "pilot_bp_in_root": tot_bp,
        "overall_novel_per_gbp": round(k_novel / tot_gbp, 2) if tot_gbp > 0 else 0.0,
        "overall_novel_nonsingleton_per_gbp": round(int(summary_row[4]) / tot_gbp, 2) if tot_gbp > 0 else 0.0,
        "by_habitat1": hab1_nov,
    }
    (wdir / "novelty_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "by_habitat1"}))
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_exp = sub.add_parser("export-faa")
    p_exp.add_argument("--ref-root", required=True)
    p_exp.add_argument("--ref-release", default="r26.11.1")
    p_exp.add_argument("--ref-stages", default=None)
    p_exp.add_argument("--pilot-root", required=True)
    p_exp.add_argument("--pilot-stages", required=True)
    p_exp.add_argument("--work-dir", required=True)
    p_exp.add_argument("--direct-createdb", action="store_true")

    p_lin = sub.add_parser("linclust-50")
    p_lin.add_argument("--work-dir", required=True)
    p_lin.add_argument("--threads", type=int, default=16)
    p_lin.add_argument("--align-parts", type=int, default=4)
    p_lin.add_argument("--split-mem", default="80G")
    p_lin.add_argument("--cleanup-db", action="store_true")

    p_ana = sub.add_parser("analyze")
    p_ana.add_argument("--work-dir", required=True)
    p_ana.add_argument("--pilot-root", required=True)
    p_ana.add_argument("--inputs-dir", required=True)
    p_ana.add_argument("--selection", required=True)
    p_ana.add_argument("--threads", type=int, default=8)
    p_ana.add_argument("--duck-mem", default="8GB")
    p_ana.add_argument("--cleanup-tsv", action="store_true")

    args = ap.parse_args(argv)
    if args.cmd == "export-faa":
        return cmd_export_faa(args)
    if args.cmd == "linclust-50":
        return cmd_linclust_50(args)
    if args.cmd == "analyze":
        return cmd_analyze(args)
    return 2


if __name__ == "__main__":
    sys.exit(main())
