"""Short-contig measurements over the pilot store, in a fresh process (measurements file section 6c). Throwaway.

    python probe_cami.py WORKDIR RELEASE [--threads 8] [--n 1000] [--flank 1000]

Reads WORKDIR/result.json (from pilot_cami.py) and writes WORKDIR/probe.json after every part.
Card B-18 moved it to the phase-B layout: tapes and tape tokens under `tape_version=`, every table
read through its part files, `seq_locator` a directory per stage, and the `SeqStore` that no longer
loads its index (card B-13), so the open cost is measured with the contig counts read from the index
files.
"""

import argparse
import json
import platform
import random
import resource
import sys
import threading
import time
import traceback
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("workdir", type=Path)
ap.add_argument("release")
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--n", type=int, default=1000)
ap.add_argument("--flank", type=int, default=1000)
ap.add_argument("--query-timeout", type=int, default=1800)
a = ap.parse_args()

KIB = 1024 if sys.platform == "linux" else 1024 * 1024


def rss_now_mb() -> float:
    with open("/proc/self/statm") as fh:
        return round(int(fh.read().split()[1]) * resource.getpagesize() / 2**20, 1)


def peak_mb() -> float:
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / KIB, 1)


W = a.workdir
ROOT = W / "root"
OUT = W / "probe.json"
P: dict = {"host": platform.node(), "started_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
           "parts": {}, "failed": {}}
P["rss_at_start_mb"] = rss_now_mb()

# A. SeqStore open, before anything else is loaded
import pyarrow as pa  # noqa: E402

pa.set_cpu_count(a.threads)
pa.set_io_thread_count(a.threads)
import pyarrow.compute as pc  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from biobox import io  # noqa: E402
from biobox.l0seq.store import SeqStore  # noqa: E402

rss0 = rss_now_mb()
t = time.perf_counter()
store = SeqStore(ROOT)
open_s = time.perf_counter() - t
rss1 = rss_now_mb()
n_loc = sum(pq.ParquetFile(f).metadata.num_rows for f in store.locator_files)
n_default = sum(pc.sum(pq.read_table(f, columns=["is_default_variant"]).column(0)).as_py() or 0 for f in store.index_files)
P["seqstore_open"] = {"open_s": round(open_s, 3), "rss_before_mb": rss0, "rss_after_mb": rss1,
                      "rss_growth_mb": round(rss1 - rss0, 1), "peak_after_mb": peak_mb(),
                      "locator_entries": n_loc, "default_contigs": n_default,
                      "index_footers_read": store.locator.footers_read + store.index.footers_read,
                      "index_row_groups_read_on_open": store.locator.row_groups_read + store.index.row_groups_read}
n_m = n_default / 1e6
P["seqstore_open"]["open_s_per_million_contigs"] = round(open_s / n_m, 3)
P["seqstore_open"]["rss_growth_mb_per_million_contigs"] = round((rss1 - rss0) / n_m, 1)


def save():
    tmp = OUT.with_suffix(".tmp")
    tmp.write_text(json.dumps(P, indent=1, default=str))
    tmp.replace(OUT)


save()

from biobox import catalog  # noqa: E402
from biobox.l1.genecall import GENE_CALL_VERSION  # noqa: E402
from biobox.l3.tapes import tape_version  # noqa: E402

V = GENE_CALL_VERSION
TV = tape_version([V])  # the gene-call producer alone (card B-07)
R = json.loads((W / "result.json").read_text())
STAGE = R["stage"]
plan = R.get("genecall_plan", {})


def part(name, fn):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {name}: start", flush=True)
    t0 = time.perf_counter()
    try:
        out, err = fn(), None
    except BaseException as e:  # noqa: BLE001
        out, err = None, f"{type(e).__name__}: {e}"
        P["failed"][name] = {"error": err, "traceback": traceback.format_exc()}
    P["parts"][name] = {"wall_s": round(time.perf_counter() - t0, 2), "peak_rss_mb": peak_mb(),
                        "rss_now_mb": rss_now_mb(), "result": out, "error": err}
    save()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {name}: {P['parts'][name]['wall_s']} s err={err}", flush=True)
    return out


# the gene-called contigs (all, or the driver's seeded subset)
idx = [r for f in io.part_files(ROOT / f"l0seq/contig_index/stage={STAGE}")
       for r in pq.read_table(f, columns=["contig_uid", "contig_len", "is_default_variant"]).to_pylist()]
default = sorted((r["contig_uid"], r["contig_len"]) for r in idx if r["is_default_variant"])
del idx
frac = plan.get("fraction", 1.0)
pick = default if frac >= 1 else sorted(random.Random(1109).sample(default, round(frac * len(default))))
called_len = dict(pick)
P["called"] = {"fraction": frac, "contigs": len(pick), "bp": sum(called_len.values()),
               "all_contigs": len(default)}


# A2. catalog open
def cat_open():
    out = {}
    cat = ROOT / f"catalog/release={a.release}.duckdb"
    r0 = rss_now_mb()
    t0 = time.perf_counter()
    con = catalog.open_catalog(cat, threads=a.threads, memory_limit="32GB")
    out["open_catalog_s"] = round(time.perf_counter() - t0, 3)
    out["rss_growth_mb"] = round(rss_now_mb() - r0, 1)
    t0 = time.perf_counter()
    out["contig_index_rows"] = con.execute("SELECT count(*) FROM contig_index").fetchall()[0][0]
    out["first_count_s"] = round(time.perf_counter() - t0, 3)
    out["catalog_bytes"] = cat.stat().st_size
    con.close()
    return out


part("catalog_open", cat_open)


# B. index and small-table sizes per contig
def sizes():
    n = len(default)
    out = {"contigs": n}
    sh = sorted((ROOT / f"l0seq/seq_shards/stage={STAGE}").glob("shard-*.fa.gz"))
    out["shards"] = len(sh)
    groups = {
        "shards": sh,
        "fai": [Path(str(p) + ".fai") for p in sh],
        "gzi": [Path(str(p) + ".gzi") for p in sh],
        "seq_locator": io.part_files(ROOT / "l0seq/seq_locator", "stage=*/part-*.parquet"),
        "contig_index": sorted((ROOT / "l0seq/contig_index").glob("stage=*/part-*.parquet")),
        "contig_occurrence": sorted((ROOT / "l0seq/contig_occurrence").glob("stage=*/source=*/part-*.parquet")),
    }
    for k, fs in groups.items():
        b = sum(f.stat().st_size for f in fs)
        out[k] = {"files": len(fs), "bytes": b, "bytes_per_contig": round(b / n, 2)}
    out["gzi_entries"] = sum((f.stat().st_size - 8) // 16 for f in groups["gzi"])
    return out


part("sizes", sizes)


# C. genes per contig, partial genes, genes per kb (over the called contigs)
def genes():
    files = io.part_files(ROOT / f"l3/tapes/tape_version={TV}", "stage=*/part-*.parquet")
    t = pa.concat_tables(pq.read_table(f, columns=["contig_uid", "contig_len", "n_slots"]) for f in files)
    cu, ln, ns = (t.column(c).to_pylist() for c in ("contig_uid", "contig_len", "n_slots"))
    bins = Counter()
    by_len = {}
    tot_g = tot_bp = 0
    edges = [(0, 1000, "<1kb"), (1000, 10000, "1-10kb"), (10000, 100000, "10-100kb"), (100000, 10**12, ">=100kb")]
    for c, L, n in zip(cu, ln, ns):
        if c not in called_len:
            continue
        b = "0" if n == 0 else "1" if n == 1 else "2-5" if n <= 5 else "6-20" if n <= 20 else ">20"
        bins[b] += 1
        tot_g += n
        tot_bp += L
        for lo, hi, lab in edges:
            if lo <= L < hi:
                d = by_len.setdefault(lab, {"contigs": 0, "genes": 0, "bp": 0, "no_gene": 0})
                d["contigs"] += 1
                d["genes"] += n
                d["bp"] += L
                d["no_gene"] += n == 0
    for d in by_len.values():
        d["genes_per_contig"] = round(d["genes"] / d["contigs"], 3)
        d["genes_per_kb"] = round(d["genes"] / (d["bp"] / 1e3), 3)
        d["fraction_no_gene"] = round(d["no_gene"] / d["contigs"], 4)
    nc = sum(bins.values())
    occ = pa.concat_tables(pq.read_table(f, columns=["dist_to_contig_end_bp"]) for f in
                           io.part_files(ROOT / f"l1/protein_occurrence/gene_call_version={V}/stage={STAGE}"))
    partial = pc.sum(pc.equal(occ.column("dist_to_contig_end_bp"), 0)).as_py() or 0
    return {"contigs": nc, "genes": tot_g, "bp": tot_bp,
            "genes_per_contig_bins": {k: bins.get(k, 0) for k in ("0", "1", "2-5", "6-20", ">20")},
            "fraction_contigs_no_gene": round(bins.get("0", 0) / nc, 4),
            "genes_flagged_partial": partial, "fraction_genes_partial": round(partial / max(1, occ.num_rows), 4),
            "genes_per_kb": round(tot_g / (tot_bp / 1e3), 4), "mean_genes_per_contig": round(tot_g / nc, 3),
            "by_contig_length": by_len, "tape_rows_all_contigs": t.num_rows}


part("genes", genes)


# E. tapes layout
def tape_layout():
    out = {}
    for name, glob in (("contig_tapes", f"l3/tapes/tape_version={TV}"),
                       ("tape_tokens", f"l3/tape_tokens/tape_version={TV}")):
        files = sorted((ROOT / glob).rglob("part-*.parquet"))
        rg_rows, b, rows = [], 0, 0
        for f in files:
            md = pq.ParquetFile(f).metadata
            rows += md.num_rows
            b += f.stat().st_size
            rg_rows += [md.row_group(i).num_rows for i in range(md.num_row_groups)]
        out[name] = {"files": len(files), "rows": rows, "bytes": b, "row_groups": len(rg_rows),
                     "contigs_per_row_group_mean": round(rows / max(1, len(rg_rows)), 1),
                     "contigs_per_row_group_max": max(rg_rows) if rg_rows else 0,
                     "bytes_per_contig": round(b / max(1, rows), 2)}
    slots = (R["steps"].get("tapes", {}).get("result") or {}).get("slots") or 0
    for k in out:
        out[k]["bytes_per_slot"] = round(out[k]["bytes"] / max(1, slots), 3)
    out["slots"] = slots
    return out


part("tape_layout", tape_layout)

rng = random.Random(4242)
with_genes = None


def summarize(prof):
    by = {}
    for s in prof["scans"]:
        names = " ".join(s["files"])
        key = next((k for k in ("contig_index", "seq_locator", "contig_occurrence", "l3/tapes", "tape_tokens",
                                "protein_occurrence", "token_postings", "catalog_info") if k in names), names[:80] or "other")
        d = by.setdefault(key, {"scans": 0, "row_groups_scanned": 0, "row_groups_total": 0, "rows_scanned": 0})
        d["scans"] += 1
        d["row_groups_scanned"] += s["row_groups_scanned"] or 0
        d["row_groups_total"] += s["row_groups_total"] or 0
        d["rows_scanned"] += s["rows_scanned"] or 0
    return {"by_table": by, "bytes_read": prof["bytes_read"], "bgzf": prof["bgzf"], "n_scans": len(prof["scans"])}


def with_timeout(con, fn):
    timer = threading.Timer(a.query_timeout, con.interrupt)
    timer.start()
    try:
        return fn()
    finally:
        timer.cancel()


# tapes_of: 1 contig and 1,000 random contigs (from the called contigs that have genes)
def tapes_of():
    files = io.part_files(ROOT / f"l3/tapes/tape_version={TV}", "stage=*/part-*.parquet")
    t = pa.concat_tables(pq.read_table(f, columns=["contig_uid", "n_slots"]) for f in files)
    pool = sorted(c for c, n in zip(t.column("contig_uid").to_pylist(), t.column("n_slots").to_pylist())
                  if n and c in called_len)
    out = {"pool_contigs_with_genes": len(pool)}
    for label, k in (("one_contig", 1), ("contigs_1000", a.n)):
        cu = rng.sample(pool, k)
        con = catalog.open_catalog(ROOT / f"catalog/release={a.release}.duckdb", threads=a.threads, memory_limit="32GB")
        t0 = time.perf_counter()
        prof = with_timeout(con, lambda: catalog.profile_query(con, "SELECT sum(n_slots) FROM tapes_of(?)", [cu]))
        wall = time.perf_counter() - t0
        out[label] = {"contigs": k, "wall_s": round(wall, 3), "sum_n_slots": prof["result"][0][0], **summarize(prof)}
        con.close()
    return out


part("tapes_of", tapes_of)

# F/G. 1,000 random gene-plus-flank fetches
occ = pa.concat_tables(pq.read_table(f, columns=["contig_uid", "start", "end", "strand"]) for f in
                       io.part_files(ROOT / f"l1/protein_occurrence/gene_call_version={V}/stage={STAGE}"))
gi = sorted(rng.sample(range(occ.num_rows), min(a.n, occ.num_rows)))
sub = occ.take(pa.array(gi))
del occ
picks = []
for c, s, e, st in zip(*(sub.column(k).to_pylist() for k in ("contig_uid", "start", "end", "strand"))):
    L = called_len[c]
    picks.append((c, max(0, s - a.flank), min(L, e + a.flank), st))
rng.shuffle(picks)
P["fetch_picks"] = {"n": len(picks), "flank_bp": a.flank,
                    "mean_bp": round(sum(e - s for _, s, e, _ in picks) / max(1, len(picks)), 1)}
ref: list = []


def fetch_store():
    st = SeqStore(ROOT)  # fresh: empty block cache
    c0 = st.counters()
    t0 = time.perf_counter()
    for c, s, e, strand in picks:
        ref.append(st.get(c, s, e, strand))
    wall = time.perf_counter() - t0
    c1 = st.counters()
    d = {k: c1[k] - c0[k] for k in c1}
    return {"wall_s": round(wall, 3), "ms_per_fetch": round(1000 * wall / len(picks), 3),
            "bases": sum(len(x) for x in ref), **d, "blocks_per_fetch": round(d["blocks_read"] / len(picks), 3)}


part("fetch_seqstore", fetch_store)


def fetch_catalog():
    tbl = pa.table({"i": pa.array(range(len(picks)), pa.int64()),
                    "cuid": pa.array([p[0] for p in picks], pa.binary()),
                    "s": pa.array([p[1] for p in picks], pa.int64()),
                    "e": pa.array([p[2] for p in picks], pa.int64()),
                    "strand": pa.array([p[3] for p in picks], pa.int64())})
    sql = "SELECT i, seq_get(cuid, s, e, strand) AS seq FROM picks ORDER BY i"
    out = {}
    cat = ROOT / f"catalog/release={a.release}.duckdb"
    con = catalog.open_catalog(cat, threads=a.threads, memory_limit="32GB")
    con.register("picks", tbl)
    t0 = time.perf_counter()
    rows = with_timeout(con, lambda: con.execute(sql).fetchall())
    out["plain_wall_s"] = round(time.perf_counter() - t0, 3)
    out["plain_bgzf"] = catalog.seq_reader(con).counters()
    out["rows"] = len(rows)
    if ref:
        out["mismatches_vs_seqstore"] = sum(1 for (i, s), r in zip(rows, ref) if s != r.decode("ascii"))
    con.close()
    con = catalog.open_catalog(cat, threads=a.threads, memory_limit="32GB")
    con.register("picks", tbl)
    t0 = time.perf_counter()
    prof = with_timeout(con, lambda: catalog.profile_query(con, sql))
    out["profiled_wall_s"] = round(time.perf_counter() - t0, 3)
    out["profile"] = summarize(prof)
    con.close()
    return out


part("fetch_catalog", fetch_catalog)


# H. clusters
def clusters():
    t = pa.concat_tables(pq.read_table(f, columns=["level", "cluster_uid"]) for f in
                         io.part_files(ROOT / f"l2/cluster_member/release={a.release}"))
    g = t.group_by(["level", "cluster_uid"]).aggregate([([], "count_all")])
    out = {}
    for lvl in ("nr90", "nr50"):
        sz = g.filter(pc.equal(g.column("level"), lvl)).column("count_all")
        n = len(sz)
        single = pc.sum(pc.equal(sz, 1)).as_py() or 0
        out[lvl] = {"clusters": n, "members": pc.sum(sz).as_py(), "singletons": single,
                    "fraction_singleton": round(single / max(1, n), 4), "max_size": pc.max(sz).as_py()}
    return out


part("clusters", clusters)
P["finished_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
save()
print("done", flush=True)
