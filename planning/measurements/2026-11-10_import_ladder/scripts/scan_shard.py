"""Scan one protein shard with one HMM library through biobox.l1.domains.scan.

    python scan_shard.py SHARD_INPUT.parquet HMM_PATH CUTOFF OUT.parquet [--cpus 4]

Writes OUT.parquet (l1_domains schema, sorted as domains.scan sorts) and OUT with .json:
wall time, CPU-minutes (all threads of this process), peak RSS, counts. On failure the JSON
holds the traceback and the exit code is 1.
"""

import argparse
import json
import platform
import resource
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from biobox.l1 import domains
from biobox.registry.schemas import table

ap = argparse.ArgumentParser()
ap.add_argument("inp", type=Path)
ap.add_argument("hmm", type=Path)
ap.add_argument("cutoff")
ap.add_argument("out", type=Path)
ap.add_argument("--cpus", type=int, default=4)
a = ap.parse_args()

meta_path = a.out.with_suffix(".json")
R = {"host": platform.node(), "input": str(a.inp), "hmm": str(a.hmm), "cutoff": a.cutoff, "cpus": a.cpus,
     "started_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), "ok": False}


def cpu_s() -> float:
    u = resource.getrusage(resource.RUSAGE_SELF)
    return u.ru_utime + u.ru_stime


def save() -> None:
    tmp = meta_path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(R, indent=1))
    tmp.replace(meta_path)


save()
t0 = time.perf_counter()
try:
    t = pq.read_table(a.inp, columns=["protein_uid", "seq"])
    proteins = dict(zip(t.column("protein_uid").to_pylist(), t.column("seq").to_pylist()))
    R["proteins"] = len(proteins)
    R["residues"] = sum(len(s.replace("*", "")) for s in proteins.values())
    rows = domains.scan(proteins, a.hmm, a.cutoff, cpus=a.cpus)
    t_scan = time.perf_counter() - t0
    out_t = pa.Table.from_pylist(rows, schema=table("l1_domains").schema)
    tmp = a.out.with_suffix(".parquet.tmp")
    pq.write_table(out_t, tmp, compression="zstd", row_group_size=122_880)
    tmp.replace(a.out)
    R.update({"ok": True, "domains": len(rows), "proteins_with_hit": len({r["protein_uid"] for r in rows}),
              "models_with_hit": len({r["model_acc"] for r in rows}), "scan_wall_s": round(t_scan, 1)})
except BaseException as e:  # noqa: BLE001 - recorded as a finding
    R["error"] = f"{type(e).__name__}: {e}"
    R["traceback"] = traceback.format_exc()
R["wall_s"] = round(time.perf_counter() - t0, 1)
R["cpu_min"] = round(cpu_s() / 60, 2)
R["peak_rss_mb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
R["finished_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
if R["ok"] and R["residues"]:
    R["cpu_min_per_mres"] = round(R["cpu_min"] / (R["residues"] / 1e6), 3)
save()
print(json.dumps({k: v for k, v in R.items() if k != "traceback"}))
raise SystemExit(0 if R["ok"] else 1)
