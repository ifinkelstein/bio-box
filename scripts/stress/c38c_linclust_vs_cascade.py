#!/usr/bin/env python3
"""Card C-38c: Linclust alone against our cascade on r26.11.1's proteins.

Implements:
1. Extraction of the exact redundancy-phase linclust module calls from
   `biobox.l2.cluster.workflow` at NR90 and NR50 (no prefilter calls; `--shuffle 0`,
   release identity/coverage/cluster-mode parameters).
2. Bounded DuckDB comparison of linclust vs published cascade (`cluster_member`,
   `cluster_member_aln`): clusters/singletons per level, cascade NR50 -> linclust
   split distribution, missed-join identity-band histogram, and `N19`/`N20` 20%-rule
   checks (`AGENTS.md` §6).
3. Scratch release view builder (symlinking read-only stage/dim files, writing linclust
   `cluster_member` and `member_cluster` + manifest) so `bio-box nbr summarize` runs
   the exact production code path on the TniQ and RNA_pol worked examples.
4. Neighbour-summary comparison helper and extrapolations to 350M and 1B proteins.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from collections.abc import Sequence
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
from blake3 import blake3

from biobox import io, layout
from biobox.l2 import cluster
from biobox.registry.schemas import table

# Design note constants (`planning/design/C-38_import_depths.md` §1 & §4.2)
N19_NR90_SHARE_EST = 0.074  # 7.4% of NR90 cascade CPU (rN.3: 73.0 / 992 CPU-min)
N19_BOTH_SHARE_EST = 0.06  # ~6% of both levels (100k sample)
N19_CPU_PER_M_EST = 14.0  # ~640 CPU-min at 45.8M proteins (twice 6.9 CPU-min/M NR90 rate)
N20_BYTES_PER_SEQ_EST = 336.0  # MMseqs2 wiki linclust memory formula

IDENTITY_BANDS = ("<50%", "50-60%", "60-70%", "70-80%", "80-90%", ">=90%")
SPLIT_BUCKETS = ("1", "2", "3-5", "6-10", "11-50", ">50")


def createdb_argv(faa: str = "in.faa", db: str = "db") -> tuple[str, ...]:
    """Exact `createdb` call from `biobox.l2.cluster._make_db` (`--shuffle 0`)."""
    return ("createdb", faa, db, "--dbtype", "1", "--shuffle", "0", "--write-lookup", "1")


def linclust_redundancy_calls(level: str) -> list[cluster.Call]:
    """Extract the exact 10 redundancy-phase linclust calls from `cluster.workflow(level)`.

    Stops at `mergeclusters` into `{w}/clu_redundancy` (excludes the cascade's
    `input_step_redundancy` `createsubdb` and all `prefilter`/`align`/`clust` cascade steps).
    """
    wf = cluster.workflow(level)
    calls: list[cluster.Call] = []
    for i, argv in enumerate(wf):
        c = cluster.Call(f"linclust.{level}.{i:02d}.{argv[0]}", argv)
        if c.module == "prefilter":
            raise AssertionError(f"unexpected prefilter call in redundancy phase: {c}")
        calls.append(c)
        if c.module == "mergeclusters" and c.output.endswith("/clu_redundancy"):
            break
    if len(calls) != 10:
        raise RuntimeError(f"{level}: expected 10 redundancy-phase calls, got {len(calls)}")
    return calls


def linclust_pipeline_calls() -> list[cluster.Call]:
    """Full linclust-alone two-level pipeline (NR90 -> rep90 -> NR50 -> mergeclusters)."""
    out: list[cluster.Call] = []
    for level, final_out in (("nr90", "clu90"), ("nr50", "clu50r")):
        red = linclust_redundancy_calls(level)
        for idx, c in enumerate(red):
            if idx == len(red) - 1:
                argv = list(c.argv)
                argv[1 + c.out_index] = final_out
                out.append(cluster.Call(c.name, tuple(argv)))
            else:
                out.append(c)
        if level == "nr90":
            out.append(cluster.Call("cluster.rep90", ("createsubdb", "clu90", "db", "rep90")))
    out.append(cluster.Call("cluster.merge", ("mergeclusters", "db", "clu50", "clu90", "clu50r")))
    return out


def identity_band(pident: float) -> str:
    """Map fractional sequence identity (`cluster_member_aln.pident` in [0, 1]) to a band."""
    if pident < 0.5:
        return "<50%"
    if pident < 0.6:
        return "50-60%"
    if pident < 0.7:
        return "60-70%"
    if pident < 0.8:
        return "70-80%"
    if pident < 0.9:
        return "80-90%"
    return ">=90%"


def check_twenty_percent_rule(measured: float, estimate: float) -> dict:
    """Compare a measured quantity against a design-note estimate (`AGENTS.md` §6)."""
    if estimate <= 0:
        raise ValueError(f"estimate must be positive, got {estimate}")
    ratio = measured / estimate
    pct_diff = (measured - estimate) / estimate * 100.0
    return {
        "measured": measured,
        "estimate": estimate,
        "ratio": ratio,
        "pct_diff": pct_diff,
        "exceeds_20pct": abs(pct_diff) > 20.0,
    }


def extrapolate_linclust(
    measured_n: int,
    measured_cpu_min: float,
    measured_peak_rss_bytes: int,
    target_ns: Sequence[int] = (350_000_000, 1_000_000_000),
    *,
    cascade_n0: int = 2_053_075,
    cascade_t0: float = 284.10,
    cascade_exp: float = 1.244,
    cascade_rss_bytes_per_prot: float | None = None,
) -> list[dict]:
    """Extrapolate linclust CPU-minutes (linear and n log n) and peak RSS (linear) to target_ns.

    Compares against the cascade N17 power law (`284.10 * (n / 2,053,075)^1.244`) and N18 linear
    cascade memory (`178 GB / 45.8M` bytes/protein by default).
    """
    if measured_n <= 0 or measured_cpu_min <= 0 or measured_peak_rss_bytes <= 0:
        raise ValueError("measured_n, measured_cpu_min and measured_peak_rss_bytes must be positive")
    cas_bpp = (178.0 * 1e9 / 45_815_843) if cascade_rss_bytes_per_prot is None else cascade_rss_bytes_per_prot
    lin_bpp = measured_peak_rss_bytes / measured_n
    rows = []
    for n in target_ns:
        scale = n / measured_n
        nlogn_scale = (n * math.log2(n)) / (measured_n * math.log2(measured_n))
        cpu_lin = measured_cpu_min * scale
        cpu_nlogn = measured_cpu_min * nlogn_scale
        rss_gib = (measured_peak_rss_bytes * scale) / (1024**3)
        cas_cpu = cascade_t0 * (n / cascade_n0) ** cascade_exp
        cas_rss_gib = (cas_bpp * n) / (1024**3)
        rows.append({
            "target_proteins": int(n),
            "scale_vs_measured": scale,
            "linclust_cpu_min_linear": cpu_lin,
            "linclust_cpu_min_nlogn": cpu_nlogn,
            "linclust_cpu_hours_linear": cpu_lin / 60.0,
            "linclust_cpu_hours_nlogn": cpu_nlogn / 60.0,
            "linclust_peak_rss_gib_linear": rss_gib,
            "linclust_bytes_per_protein": lin_bpp,
            "cascade_cpu_min_n17": cas_cpu,
            "cascade_cpu_hours_n17": cas_cpu / 60.0,
            "cascade_peak_rss_gib_n18_linear": cas_rss_gib,
            "speedup_vs_cascade_linear": cas_cpu / cpu_lin,
            "speedup_vs_cascade_nlogn": cas_cpu / cpu_nlogn,
            "memory_reduction_factor": cas_rss_gib / rss_gib,
            "min_parts_240_cpu_min_linear": math.ceil(cpu_lin / 240.0),
            "min_parts_240_cpu_min_nlogn": math.ceil(cpu_nlogn / 240.0),
            "label": "derived (extrapolation)",
        })
    return rows


# ---------------------------------------------------------------------------
# Bounded DuckDB comparison of linclust vs published cascade
# ---------------------------------------------------------------------------


def _sql_list(paths: Sequence[Path]) -> str:
    return "[" + ", ".join("'" + str(Path(p)).replace("'", "''") + "'" for p in paths) + "]"


def _sql_path(p: Path) -> str:
    return "'" + str(Path(p)).replace("'", "''") + "'"


def _register_linclust_view(con: duckdb.DuckDBPyConnection, name: str, path: Path) -> None:
    p = Path(path)
    if p.suffix == ".parquet":
        con.execute(f"CREATE OR REPLACE VIEW {name} AS SELECT rep_uid, member_uid FROM read_parquet({_sql_path(p)})")
    else:
        con.execute(
            f"CREATE OR REPLACE VIEW {name} AS "
            f"SELECT unhex(rep_hex) AS rep_uid, unhex(member_hex) AS member_uid "
            f"FROM read_csv({_sql_path(p)}, delim='\t', header=false, "
            f"columns={{'rep_hex': 'VARCHAR', 'member_hex': 'VARCHAR'}})"
        )


def compare_clusterings_duckdb(
    cascade_member_paths: Sequence[Path],
    cascade_aln_paths: Sequence[Path],
    linclust_by_level: dict[str, Path],
    *,
    memory_limit: str = "8GB",
    threads: int = 2,
    tmp_dir: Path | None = None,
) -> dict:
    """Compare linclust against cascade in bounded DuckDB queries.

    Returns level summaries, cascade NR50 -> linclust split distributions (overall and by size
    strata), and missed-join histograms stratified by `cluster_member_aln.pident` identity bands.
    """
    with tempfile.TemporaryDirectory(prefix="c38c-duck-", dir=tmp_dir) as d:
        con = duckdb.connect(config={"threads": threads, "memory_limit": memory_limit, "temp_directory": d})
        try:
            con.execute(
                f"CREATE VIEW cas_mem AS SELECT member_uid, level, cluster_uid, rep_uid "
                f"FROM read_parquet({_sql_list(cascade_member_paths)}, hive_partitioning=false)"
            )
            con.execute(
                f"CREATE VIEW cas_aln AS SELECT member_uid, level, rep_uid, pident, alnlen "
                f"FROM read_parquet({_sql_list(cascade_aln_paths)}, hive_partitioning=false)"
            )
            for lvl in ("nr90", "nr50"):
                _register_linclust_view(con, f"lin_{lvl}", linclust_by_level[lvl])

            levels_out: dict[str, dict] = {}
            missed_out: dict[str, list[dict]] = {}
            reverse_out: dict[str, dict] = {}

            for lvl in ("nr90", "nr50"):
                c_prot, c_clu, c_sing = con.execute(
                    f"""
                    WITH g AS (
                        SELECT rep_uid, count(*) AS sz
                        FROM cas_mem WHERE level = '{lvl}'
                        GROUP BY rep_uid
                    )
                    SELECT coalesce(sum(sz), 0), count(*), count(*) FILTER (WHERE sz = 1) FROM g
                    """
                ).fetchone()
                l_prot, l_clu, l_sing = con.execute(
                    f"""
                    WITH g AS (
                        SELECT rep_uid, count(*) AS sz
                        FROM lin_{lvl}
                        GROUP BY rep_uid
                    )
                    SELECT coalesce(sum(sz), 0), count(*), count(*) FILTER (WHERE sz = 1) FROM g
                    """
                ).fetchone()
                if c_prot != l_prot:
                    raise ValueError(f"{lvl}: protein count mismatch (cascade={c_prot}, linclust={l_prot})")
                levels_out[lvl] = {
                    "n_proteins": int(c_prot),
                    "cascade_clusters": int(c_clu),
                    "cascade_singletons": int(c_sing),
                    "cascade_singleton_share": c_sing / c_clu if c_clu else 0.0,
                    "cascade_singleton_protein_share": c_sing / c_prot if c_prot else 0.0,
                    "linclust_clusters": int(l_clu),
                    "linclust_singletons": int(l_sing),
                    "linclust_singleton_share": l_sing / l_clu if l_clu else 0.0,
                    "linclust_singleton_protein_share": l_sing / l_prot if l_prot else 0.0,
                    "delta_clusters": int(l_clu - c_clu),
                    "ratio_clusters": l_clu / c_clu if c_clu else 0.0,
                    "delta_singletons": int(l_sing - c_sing),
                    "ratio_singletons": l_sing / c_sing if c_sing else 0.0,
                }

                # Missed joins stratified by cluster_member_aln.pident identity bands
                band_rows = con.execute(
                    f"""
                    WITH lin_sz AS (
                        SELECT rep_uid, count(*) AS lin_size FROM lin_{lvl} GROUP BY rep_uid
                    ),
                    nonrep AS (
                        SELECT cm.member_uid, cm.rep_uid AS cas_rep, al.pident,
                               lm.rep_uid AS m_lin_rep, ls.lin_size AS m_lin_size,
                               lr.rep_uid AS r_lin_rep
                        FROM cas_mem cm
                        JOIN cas_aln al ON al.member_uid = cm.member_uid AND al.level = cm.level
                        JOIN lin_{lvl} lm ON lm.member_uid = cm.member_uid
                        JOIN lin_sz ls ON ls.rep_uid = lm.rep_uid
                        JOIN lin_{lvl} lr ON lr.member_uid = cm.rep_uid
                        WHERE cm.level = '{lvl}' AND cm.member_uid <> cm.rep_uid
                    ),
                    binned AS (
                        SELECT
                            CASE
                                WHEN pident < 0.5 THEN '<50%'
                                WHEN pident < 0.6 THEN '50-60%'
                                WHEN pident < 0.7 THEN '60-70%'
                                WHEN pident < 0.8 THEN '70-80%'
                                WHEN pident < 0.9 THEN '80-90%'
                                ELSE '>=90%'
                            END AS band,
                            (m_lin_rep <> r_lin_rep) AS missed,
                            (m_lin_rep <> r_lin_rep AND m_lin_size = 1) AS missed_singleton
                        FROM nonrep
                    )
                    SELECT band,
                           count(*) AS cascade_nonrep,
                           count(*) FILTER (WHERE NOT missed) AS joined_by_linclust,
                           count(*) FILTER (WHERE missed) AS missed_by_linclust,
                           count(*) FILTER (WHERE missed_singleton) AS missed_to_singleton
                    FROM binned
                    GROUP BY band
                    """
                ).fetchall()
                by_band = {r[0]: r for r in band_rows}
                missed_list = []
                for b in IDENTITY_BANDS:
                    _, tot, joined, missed, sing = by_band.get(b, (b, 0, 0, 0, 0))
                    missed_list.append({
                        "band": b,
                        "cascade_nonrep": int(tot),
                        "joined_by_linclust": int(joined),
                        "missed_by_linclust": int(missed),
                        "missed_to_singleton": int(sing),
                        "missed_share": (missed / tot) if tot else 0.0,
                        "singleton_share_of_missed": (sing / missed) if missed else 0.0,
                    })
                missed_out[lvl] = missed_list

                # Reverse check: linclust non-rep members whose cascade cluster differs from linclust rep's
                rev_tot, rev_diff = con.execute(
                    f"""
                    SELECT count(*) AS lin_nonrep,
                           count(*) FILTER (WHERE cm_m.cluster_uid <> cm_r.cluster_uid) AS not_in_same_cascade
                    FROM lin_{lvl} l
                    JOIN cas_mem cm_m ON cm_m.member_uid = l.member_uid AND cm_m.level = '{lvl}'
                    JOIN cas_mem cm_r ON cm_r.member_uid = l.rep_uid AND cm_r.level = '{lvl}'
                    WHERE l.member_uid <> l.rep_uid
                    """
                ).fetchone()
                reverse_out[lvl] = {
                    "linclust_nonrep": int(rev_tot),
                    "not_in_same_cascade_cluster": int(rev_diff),
                    "share": (rev_diff / rev_tot) if rev_tot else 0.0,
                }

            # Cascade NR50 -> linclust NR50 split distribution
            con.execute(
                """
                CREATE TEMP TABLE nr50_cluster_splits AS
                SELECT cm.cluster_uid,
                       count(*) AS cas_size,
                       count(DISTINCT l.rep_uid) AS k_lin
                FROM cas_mem cm
                JOIN lin_nr50 l ON l.member_uid = cm.member_uid
                WHERE cm.level = 'nr50'
                GROUP BY cm.cluster_uid
                """
            )
            splits_out: dict[str, dict] = {}
            for stratum, where in (
                ("all", "cas_size >= 1"),
                ("size>=2", "cas_size >= 2"),
                ("size>=10", "cas_size >= 10"),
                ("size>=100", "cas_size >= 100"),
            ):
                row = con.execute(
                    f"""
                    SELECT count(*) AS n_clu,
                           coalesce(sum(cas_size), 0) AS n_prot,
                           count(*) FILTER (WHERE k_lin = 1) AS unsplit,
                           count(*) FILTER (WHERE k_lin >= 2) AS split,
                           coalesce(avg(k_lin), 0.0) AS mean_k,
                           coalesce(quantile_cont(k_lin, 0.5), 0.0) AS med_k,
                           coalesce(quantile_cont(k_lin, 0.95), 0.0) AS p95_k,
                           coalesce(max(k_lin), 0) AS max_k
                    FROM nr50_cluster_splits WHERE {where}
                    """
                ).fetchone()
                hist_rows = con.execute(
                    f"""
                    SELECT CASE
                             WHEN k_lin = 1 THEN '1'
                             WHEN k_lin = 2 THEN '2'
                             WHEN k_lin BETWEEN 3 AND 5 THEN '3-5'
                             WHEN k_lin BETWEEN 6 AND 10 THEN '6-10'
                             WHEN k_lin BETWEEN 11 AND 50 THEN '11-50'
                             ELSE '>50'
                           END AS bucket,
                           count(*) AS n_clu,
                           coalesce(sum(cas_size), 0) AS n_prot
                    FROM nr50_cluster_splits WHERE {where}
                    GROUP BY bucket
                    """
                ).fetchall()
                hmap = {r[0]: {"clusters": int(r[1]), "proteins": int(r[2])} for r in hist_rows}
                n_clu, n_prot, unsplit, split, mean_k, med_k, p95_k, max_k = row
                splits_out[stratum] = {
                    "n_cascade_clusters": int(n_clu),
                    "n_proteins": int(n_prot),
                    "unsplit_clusters": int(unsplit),
                    "split_clusters": int(split),
                    "split_share": (split / n_clu) if n_clu else 0.0,
                    "mean_linclust_clusters": float(mean_k),
                    "median_linclust_clusters": float(med_k),
                    "p95_linclust_clusters": float(p95_k),
                    "max_linclust_clusters": int(max_k),
                    "histogram": {b: hmap.get(b, {"clusters": 0, "proteins": 0}) for b in SPLIT_BUCKETS},
                }

            return {
                "levels": levels_out,
                "nr50_splits": splits_out,
                "missed_joins_by_band": missed_out,
                "reverse_joins": reverse_out,
            }
        finally:
            con.close()


# ---------------------------------------------------------------------------
# Pod execution helpers: prepare, run calls, build scratch view, compare nbr
# ---------------------------------------------------------------------------


def _run_mmseqs_timed(w: Path, argv: Sequence[str], threads: int | None, log_path: Path) -> dict:
    """Run one MMseqs2 module call and record wall time, user+sys CPU-minutes, and child ru_maxrss."""
    exe = shutil.which("mmseqs")
    if exe is None:
        raise FileNotFoundError("mmseqs not on PATH")
    cmd = [exe, *argv, "-v", "3", *(["--threads", str(threads)] if threads is not None else [])]
    r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    t0 = time.monotonic()
    with open(log_path, "w") as fh:
        fh.write("### " + " ".join(cmd) + "\n")
        fh.flush()
        proc = subprocess.run(cmd, cwd=w, stdout=fh, stderr=subprocess.STDOUT)
    wall_s = time.monotonic() - t0
    r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    if proc.returncode != 0:
        raise RuntimeError(f"mmseqs {argv[0]} exited {proc.returncode}; see {log_path}")
    cpu_s = (r1.ru_utime - r0.ru_utime) + (r1.ru_stime - r0.ru_stime)
    # Note: ru_maxrss on Linux is in KiB, cumulative high-water mark over all children in this process.
    return {
        "wall_s": round(wall_s, 3),
        "cpu_s": round(cpu_s, 3),
        "cpu_min": round(cpu_s / 60.0, 4),
        "children_maxrss_kib": int(r1.ru_maxrss),
    }


def cmd_prepare(args: argparse.Namespace) -> None:
    """Read r26.11.1 protein_seq files read-only, write sorted FASTA + createdb, remove FASTA."""
    root = Path(args.root)
    w = Path(args.work)
    w.mkdir(parents=True, exist_ok=True)
    (w / "logs").mkdir(exist_ok=True)
    plan_rec = json.loads((root / f"l2/cluster_plan/release={args.release}/cluster.plan.json").read_text())
    groups = [[root / f["path"] for f in g] for g in plan_rec["files"]]
    t0 = time.monotonic()
    n_prot = cluster._write_fasta(groups, w / "in.faa", w / "prot.parquet")
    fasta_s = time.monotonic() - t0
    cdb = _run_mmseqs_timed(w, createdb_argv("in.faa", "db"), None, w / "logs" / "createdb.log")
    (w / "in.faa").unlink()
    rec = {
        "step": "prepare",
        "release": args.release,
        "proteins": n_prot,
        "fasta_wall_s": round(fasta_s, 3),
        "createdb": cdb,
    }
    (w / "logs" / "prepare.json").write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec))


def _ensure_parent_dirs(w: Path, c: cluster.Call) -> None:
    (w / c.output).parent.mkdir(parents=True, exist_ok=True)


def _write_order_redundancy(w: Path, c: cluster.Call) -> None:
    if not c.order_file:
        return
    order = w / c.order_file
    tmp = order.with_name(f".{order.name}.tmp-{os.getpid()}")
    with open(w / f"{c.output}.index") as src, open(tmp, "w") as dst:
        for line in src:
            dst.write(line.split(None, 1)[0] + "\n" if line.strip() else "\n")
    os.replace(tmp, order)


def cmd_run_call(args: argparse.Namespace) -> None:
    """Execute one step of `linclust_pipeline_calls()` (by index or name) in `work`, writing a JSON record."""
    w = Path(args.work)
    (w / "logs").mkdir(parents=True, exist_ok=True)
    pipe = linclust_pipeline_calls()
    by_name = {c.name: c for c in pipe}
    if args.step.isdigit():
        c = pipe[int(args.step)]
    elif args.step in by_name:
        c = by_name[args.step]
    else:
        raise ValueError(f"unknown step {args.step!r}")
    _ensure_parent_dirs(w, c)
    cluster._clear(w, c.output)
    threads = None if c.module in cluster.UNTHREADED else args.threads
    log_path = w / "logs" / f"{c.name}.log"
    timing = _run_mmseqs_timed(w, c.argv, threads, log_path)
    _write_order_redundancy(w, c)
    # Count output entries
    idx_path = Path(f"{w / c.output}.index")
    n_keys = sum(1 for _ in open(idx_path)) if idx_path.is_file() else 0
    rec = {
        "step": c.name,
        "module": c.module,
        "argv": list(c.argv),
        "threads": threads,
        "output": c.output,
        "output_keys": n_keys,
        **timing,
    }
    (w / "logs" / f"{c.name}.json").write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec))


def cmd_run_split_call(args: argparse.Namespace) -> None:
    """Split a per-query call (align/rescorediagonal) into N parts using `cluster.cut`, or run one part / join."""
    w = Path(args.work)
    (w / "logs").mkdir(parents=True, exist_ok=True)
    pipe = {c.name: c for c in linclust_pipeline_calls()}
    c = pipe[args.step]
    if not c.split:
        raise ValueError(f"{c.name} ({c.module}) is not a per-query split call")
    cut_path = w / "logs" / f"{c.name}.cut.json"
    rel = f"parts/{c.name}"
    if args.action == "cut":
        pos = cluster.PER_QUERY[c.module]
        src = c.positional[pos]

        def rows():
            for k, _, n in cluster.index_rows(w / src):
                yield k, n

        parts = cluster.cut(rows, 1.0, parts=args.parts)
        doc = {"step": c.name, "src": src, "parts": parts}
        cut_path.write_text(json.dumps(doc, indent=2) + "\n")
        print(json.dumps({"step": f"{c.name}.cut", "parts": len(parts)}))
        return
    cut_doc = json.loads(cut_path.read_text())
    if args.action == "part":
        k = args.part
        p = cut_doc["parts"][k]
        pos = cluster.PER_QUERY[c.module]
        src = c.positional[pos]
        (w / rel).mkdir(parents=True, exist_ok=True)
        ids_path = w / rel / f"ids-{k:05d}"
        with open(ids_path, "w") as fh:
            for key, _, _ in cluster.index_rows(w / src):
                if p["first"] <= key <= p["last"]:
                    fh.write(f"{key}\n")
        pin, pout = f"{rel}/in-{k:05d}", f"{rel}/out-{k:05d}"
        cluster._clear(w, pin)
        cluster._clear(w, pout)
        log_path = w / "logs" / f"{c.name}.part-{k:05d}.log"
        t_sub = _run_mmseqs_timed(
            w, ("createsubdb", f"{rel}/ids-{k:05d}", src, pin, "--subdb-mode", "0"), None, log_path
        )
        argv = list(c.argv)
        argv[1 + pos] = pin
        argv[1 + c.out_index] = pout
        t_call = _run_mmseqs_timed(w, argv, args.threads, w / "logs" / f"{c.name}.part-{k:05d}.call.log")
        rec = {
            "step": f"{c.name}.part-{k:05d}",
            "module": c.module,
            "part": k,
            "keys": p["keys"],
            "subdb": t_sub,
            **t_call,
        }
        (w / "logs" / f"{c.name}.part-{k:05d}.json").write_text(json.dumps(rec, indent=2) + "\n")
        print(json.dumps(rec))
        return
    if args.action == "join":
        _ensure_parent_dirs(w, c)
        t0 = time.monotonic()
        part_paths = [w / f"{rel}/out-{k:05d}" for k in range(len(cut_doc["parts"]))]
        cluster.join_dbs(part_paths, w / c.output)
        shutil.rmtree(w / rel, ignore_errors=True)
        rec = {"step": f"{c.name}.join", "parts": len(part_paths), "wall_s": round(time.monotonic() - t0, 3)}
        (w / "logs" / f"{c.name}.join.json").write_text(json.dumps(rec, indent=2) + "\n")
        print(json.dumps(rec))


def cmd_createtsv(args: argparse.Namespace) -> None:
    """Convert `clu90` and `clu50` MMseqs2 databases to TSVs and flat Parquet."""
    w = Path(args.work)
    (w / "logs").mkdir(parents=True, exist_ok=True)
    out = {}
    for level, clu in (("nr90", "clu90"), ("nr50", "clu50")):
        t = _run_mmseqs_timed(w, ("createtsv", "db", "db", clu, f"{clu}.tsv"), args.threads, w / "logs" / f"createtsv_{level}.log")
        out[level] = t
    (w / "logs" / "createtsv.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out))


def cmd_compare(args: argparse.Namespace) -> None:
    """Run bounded DuckDB comparison against r26.11.1's published cascade."""
    root = Path(args.root)
    manifest = json.loads((root / f"manifests/release={args.release}.json").read_text())
    mem_paths = sorted(root / f["path"] for f in manifest["files"] if f["table"] == "cluster_member")
    aln_paths = sorted(root / f["path"] for f in manifest["files"] if f["table"] == "cluster_member_aln")
    res = compare_clusterings_duckdb(
        mem_paths,
        aln_paths,
        {"nr90": Path(args.clu90_tsv), "nr50": Path(args.clu50_tsv)},
        memory_limit=args.memory_limit,
        threads=args.threads,
        tmp_dir=Path(args.tmp) if args.tmp else None,
    )
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res["levels"], indent=2))


def cmd_build_view(args: argparse.Namespace) -> None:
    """Build a scratch release view over r26.11.1 with linclust cluster_member & member_cluster.

    Symlinks `r26.11.1`'s read-only `l0seq`, `l1`, `l3`, `dim` directories so `bio-box nbr summarize`
    runs the exact production code path without copying stage files.
    """
    root = Path(args.root)
    view = Path(args.view)
    release = args.release
    view.mkdir(parents=True, exist_ok=True)
    for sub in ("l0seq", "l1", "l3", "dim"):
        dst = view / sub
        if not dst.exists() and not dst.is_symlink():
            os.symlink(root / sub, dst)

    # Convert clu90.tsv / clu50.tsv into cluster_member + member_cluster under view/l2/
    mem_dir = view / f"l2/cluster_member/release={release}"
    mc_dir = view / f"l2/member_cluster/release={release}"
    if mem_dir.exists():
        shutil.rmtree(mem_dir)
    if mc_dir.exists():
        shutil.rmtree(mc_dir)

    with tempfile.TemporaryDirectory(prefix="c38c-view-", dir=args.tmp) as td:
        f = Path(td)
        for lvl, tsv in (("nr90", Path(args.clu90_tsv)), ("nr50", Path(args.clu50_tsv))):
            cluster._tsv_to_parquet(tsv, f / f"clu_{lvl}.parquet", {"rep": pa.string(), "member": pa.string()})
        con = duckdb.connect(
            config={
                "threads": 1,
                "memory_limit": args.memory_limit,
                "temp_directory": str(f / "duckdb"),
                "preserve_insertion_order": True,
            }
        )
        try:
            ph = {lvl: hashlib.sha256(f"linclust-only-{lvl}".encode()).hexdigest()[:16] for lvl in ("nr90", "nr50")}
            counts = {lvl: cluster._cluster_members(con, f, lvl, ph[lvl]) for lvl in ("nr90", "nr50")}
            mem_spec = table("cluster_member")

            def members():
                for lvl in ("nr50", "nr90"):
                    yield from io.duck_batches(
                        con,
                        f"SELECT unhex(member) AS member_uid, '{lvl}' AS level, cluster_uid, unhex(rep) AS rep_uid, "
                        f"'{ph[lvl]}' AS cluster_params_hash FROM read_parquet({_sql_path(f / f'mem_{lvl}.parquet')}) "
                        "ORDER BY cluster_uid, member_uid",
                        mem_spec,
                    )

            mem_paths = io.write_parts(members(), mem_spec, mem_dir, cluster.EST_MEMBER_ROW, layout.MAX_FILE_BYTES)
        finally:
            con.close()

        cluster.write_member_cluster(
            mem_paths, mc_dir, tmp=f, max_file_bytes=layout.MAX_FILE_BYTES, memory_limit=args.memory_limit
        )

    # Write scratch manifest replacing cluster_member and member_cluster entries
    src_man = json.loads((root / f"manifests/release={release}.json").read_text())
    new_files = [
        entry for entry in src_man["files"] if entry["table"] not in ("cluster_member", "member_cluster")
    ]
    for tbl, d in (("cluster_member", mem_dir), ("member_cluster", mc_dir)):
        for p in sorted(d.glob("part-*.parquet")):
            rel_p = p.relative_to(view).as_posix()
            new_files.append({
                "table": tbl,
                "path": rel_p,
                "partitions": {"release": release},
                "rows": pq.read_metadata(p).num_rows,
                "bytes": p.stat().st_size,
            })
    src_man["files"] = new_files
    man_dir = view / "manifests"
    man_dir.mkdir(parents=True, exist_ok=True)
    (man_dir / f"release={release}.json").write_text(json.dumps(src_man, indent=1) + "\n")
    print(json.dumps({"view": str(view), "counts": counts, "member_parts": len(mem_paths)}))


def parse_nbr_summary_tsv(path: Path) -> tuple[dict, list[dict]]:
    """Parse `bio-box nbr summarize --format tsv` output into (header_meta, rows)."""
    meta: dict[str, str] = {}
    lines = Path(path).read_text().splitlines()
    data_lines = []
    for ln in lines:
        if ln.startswith("#"):
            parts = ln[1:].strip().split("\t", 1)
            if len(parts) == 2:
                meta[parts[0].strip()] = parts[1].strip()
        elif ln.strip():
            data_lines.append(ln)
    rows = list(csv.DictReader(data_lines, delimiter="\t"))
    return meta, rows


def summarize_nbr_comparison(cascade_tsv: Path, linclust_tsv: Path, top_k: int = 15) -> dict:
    """Compare cascade vs linclust `nbr summarize --by nr50` TSV outputs."""
    c_meta, c_rows = parse_nbr_summary_tsv(cascade_tsv)
    l_meta, l_rows = parse_nbr_summary_tsv(linclust_tsv)

    def nr50_keyed(rows: list[dict]) -> list[dict]:
        return [r for r in rows if r.get("level") == "nr50" and r.get("unannotated", "").lower() != "true" and r.get("key")]

    ck = nr50_keyed(c_rows)
    lk = nr50_keyed(l_rows)
    c_loci = int(c_meta.get("n_anchor_loci", "0").strip('"'))
    l_loci = int(l_meta.get("n_anchor_loci", "0").strip('"'))

    def top_list(rows: list[dict]) -> list[dict]:
        return [
            {
                "rank": i + 1,
                "cluster_uid": r["key"],
                "model_name": r.get("model_name", ""),
                "n_loci": int(r["n_loci"]),
                "share_of_loci": float(r["share_of_loci"]),
                "n_neighbors": int(r["n_neighbors"]),
                "n_genomes": int(r["n_genomes"]),
                "modal_offset": r.get("modal_offset", ""),
            }
            for i, r in enumerate(rows[:top_k])
        ]

    # Group by dominant model_name to quantify partner fragmentation (how many NR50 clusters carry each Pfam)
    def by_model(rows: list[dict]) -> dict[str, dict]:
        out: dict[str, dict] = {}
        for r in rows:
            m = r.get("model_name", "") or "(unannotated)"
            e = out.setdefault(m, {"n_clusters": 0, "total_neighbors": 0, "top_cluster_loci": 0})
            e["n_clusters"] += 1
            e["total_neighbors"] += int(r["n_neighbors"])
            e["top_cluster_loci"] = max(e["top_cluster_loci"], int(r["n_loci"]))
        return out

    return {
        "n_anchor_loci_cascade": c_loci,
        "n_anchor_loci_linclust": l_loci,
        "n_nr50_keys_cascade": len(ck),
        "n_nr50_keys_linclust": len(lk),
        "nr50_key_ratio": (len(lk) / len(ck)) if ck else 0.0,
        "top_cascade": top_list(ck),
        "top_linclust": top_list(lk),
        "cascade_by_model": by_model(ck),
        "linclust_by_model": by_model(lk),
    }


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Card C-38c linclust vs cascade benchmark & analysis")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_prep = sub.add_parser("prepare")
    p_prep.add_argument("--root", required=True)
    p_prep.add_argument("--release", default="r26.11.1")
    p_prep.add_argument("--work", required=True)
    p_prep.set_defaults(func=cmd_prepare)

    p_call = sub.add_parser("run-call")
    p_call.add_argument("--work", required=True)
    p_call.add_argument("--step", required=True)
    p_call.add_argument("--threads", type=int, default=12)
    p_call.set_defaults(func=cmd_run_call)

    p_split = sub.add_parser("run-split-call")
    p_split.add_argument("--work", required=True)
    p_split.add_argument("--step", required=True)
    p_split.add_argument("--action", choices=("cut", "part", "join"), required=True)
    p_split.add_argument("--parts", type=int, default=4)
    p_split.add_argument("--part", type=int, default=0)
    p_split.add_argument("--threads", type=int, default=12)
    p_split.set_defaults(func=cmd_run_split_call)

    p_tsv = sub.add_parser("createtsv")
    p_tsv.add_argument("--work", required=True)
    p_tsv.add_argument("--threads", type=int, default=12)
    p_tsv.set_defaults(func=cmd_createtsv)

    p_cmp = sub.add_parser("compare")
    p_cmp.add_argument("--root", required=True)
    p_cmp.add_argument("--release", default="r26.11.1")
    p_cmp.add_argument("--clu90-tsv", required=True)
    p_cmp.add_argument("--clu50-tsv", required=True)
    p_cmp.add_argument("--memory-limit", default="8GB")
    p_cmp.add_argument("--threads", type=int, default=4)
    p_cmp.add_argument("--tmp", default=None)
    p_cmp.add_argument("--out", required=True)
    p_cmp.set_defaults(func=cmd_compare)

    p_view = sub.add_parser("build-view")
    p_view.add_argument("--root", required=True)
    p_view.add_argument("--release", default="r26.11.1")
    p_view.add_argument("--clu90-tsv", required=True)
    p_view.add_argument("--clu50-tsv", required=True)
    p_view.add_argument("--view", required=True)
    p_view.add_argument("--memory-limit", default="8GB")
    p_view.add_argument("--tmp", default=None)
    p_view.set_defaults(func=cmd_build_view)

    args = ap.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
