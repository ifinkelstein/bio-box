"""SeqStore open cost and 1,000 random +-10-gene neighborhood fetches, in a fresh process.

    python store_probe.py ROOT OUT_JSON [--n 1000] [--k 10] [--release R | --stages S1,P1]

Neighbourhoods are drawn from the tapes of the given stages (card B-25: readers take an explicit stage
list): `--release` reads them from that release's `release_def`, `--stages` lists them; with neither,
every stage directory under the gene-call tape version is used and the list is printed in the output.
"""

import argparse
import json
import platform
import random
import resource
import sys
import time
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("root", type=Path)
ap.add_argument("out", type=Path)
ap.add_argument("--n", type=int, default=1000)
ap.add_argument("--k", type=int, default=10)
g = ap.add_mutually_exclusive_group()
g.add_argument("--release", help="take the stages from this release's release_def")
g.add_argument("--stages", help="comma-separated stage list")
a = ap.parse_args()

KIB = 1024 if sys.platform == "linux" else 1024 * 1024


def rss_now_mb() -> float:
    with open("/proc/self/statm") as fh:
        return round(int(fh.read().split()[1]) * resource.getpagesize() / 2**20, 1)


def peak_mb() -> float:
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / KIB, 1)


from biobox.l0seq.store import SeqStore  # noqa: E402
from biobox.l1.genecall import GENE_CALL_VERSION  # noqa: E402
from biobox.l3.tokens import read_tapes  # noqa: E402

from biobox import catalog  # noqa: E402
from biobox.l3.tapes import tape_version  # noqa: E402

if a.release:
    stages = catalog.release_stages(a.root, a.release)
elif a.stages:
    stages = a.stages.split(",")
else:
    stages = sorted(d.name.removeprefix("stage=") for d in
                    (a.root / f"l3/tapes/tape_version={tape_version([GENE_CALL_VERSION])}").glob("stage=*"))

R = {"host": platform.node(), "root": str(a.root)}
R["stages"] = stages
R["rss_before_open_mb"] = rss_now_mb()
R["peak_before_open_mb"] = peak_mb()
t = time.perf_counter()
store = SeqStore(a.root, stages)  # the stages of the release or of --stages (card B-28)
R["open_s"] = round(time.perf_counter() - t, 3)
R["rss_after_open_mb"] = rss_now_mb()
R["peak_after_open_mb"] = peak_mb()
R["index_footers_read"] = store.locator.footers_read + store.index.footers_read
R["index_row_groups_read_on_open"] = store.locator.row_groups_read + store.index.row_groups_read
# counted from the index files after the timed open (the store itself holds no dict of contigs)
import pyarrow.compute as pc  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

loc = pq.read_table(store.locator_files, columns=["shard"])
R["locator_entries"] = loc.num_rows
R["default_contigs"] = sum(pc.sum(pq.read_table(f, columns=["is_default_variant"]).column(0)).as_py() or 0
                           for f in store.index_files)
R["shards"] = sorted(set(loc.column("shard").to_pylist()))
del loc

tp = read_tapes(a.root, GENE_CALL_VERSION, ["contig_uid", "contig_len", "n_slots", "start_bps", "end_bps"], stages=stages)
tp = [t for t in tp if t["n_slots"]]
cum, tot = [], 0
for t in tp:
    tot += t["n_slots"]
    cum.append(tot)
rng = random.Random(1109)
import bisect  # noqa: E402

picks = []
for _ in range(a.n):
    g = rng.randrange(tot)
    ci = bisect.bisect_right(cum, g)
    t = tp[ci]
    i = g - (cum[ci - 1] if ci else 0)
    lo, hi = max(0, i - a.k), min(t["n_slots"] - 1, i + a.k)
    s = min(t["start_bps"][lo : hi + 1])
    e = min(t["contig_len"], max(t["end_bps"][lo : hi + 1]))
    picks.append((t["contig_uid"], s, e))
R["neighborhood_bp"] = {"mean": round(sum(e - s for _, s, e in picks) / len(picks), 1),
                        "max": max(e - s for _, s, e in picks)}


def run(st: SeqStore, label: str) -> None:
    c0 = st.counters()
    t0 = time.perf_counter()
    nbytes = 0
    for cu, s, e in picks:
        nbytes += len(st.get(cu, s, e))
    wall = time.perf_counter() - t0
    c1 = st.counters()
    d = {k: c1[k] - c0[k] for k in c1}
    R[label] = {"wall_s": round(wall, 3), "ms_per_fetch": round(1000 * wall / len(picks), 3),
                "bases_returned": nbytes, **d,
                "blocks_per_fetch": round(d["blocks_read"] / len(picks), 3),
                "compressed_bytes_per_fetch": round(d["compressed_bytes_read"] / len(picks), 1)}


run(store, "fetch_cache64_first_pass")
run(SeqStore(a.root, stages, cache_blocks=0), "fetch_cache0")
R["peak_end_mb"] = peak_mb()
a.out.write_text(json.dumps(R, indent=1, default=str))
print(json.dumps(R, indent=1, default=str))
