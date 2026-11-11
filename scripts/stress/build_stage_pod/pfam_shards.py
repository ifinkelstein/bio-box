#!/usr/bin/env python3
"""Run the scan shards of a pfam_plan.py plan, at most --parallel at a time (card C-10). Host python3, stdlib.

    python3 pfam_shards.py WORK HMM CUTOFF --image SIF --pythonpath P --scan SCAN_SHARD_PY --stop ROOT/STOP
        [--parallel 8] [--cpus 4] [--bind SPEC ...]

Run it inside a memory-capped unit (systemd-run --user --unit=... -p MemoryMax=... -p MemorySwapMax=0
-p CPUQuota=...). Each shard is its own process, `nice -n 19 singularity exec ... python SCAN_SHARD_PY
WORK/input/shard-<k>.parquet HMM CUTOFF WORK/scan/shard-<k>.parquet --cpus N` (the ladder's
scan_shard.py, which calls biobox.l1.domains.scan); a shard whose WORK/scan/shard-<k>.json says ok is
not run again. No shard is started while the STOP file exists. Wall, CPU-minutes and peak RSS of each
shard process (os.wait4) go to WORK/shards.json, rewritten after each shard. Exit 0 when every shard
is ok, else 1.
"""

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("work", type=Path)
ap.add_argument("hmm")
ap.add_argument("cutoff")
ap.add_argument("--image", required=True)
ap.add_argument("--pythonpath", required=True)
ap.add_argument("--scan", required=True)
ap.add_argument("--stop", type=Path, required=True)
ap.add_argument("--parallel", type=int, default=8)
ap.add_argument("--cpus", type=int, default=4)
ap.add_argument("--bind", action="append", default=[])
a = ap.parse_args()

plan = json.loads((a.work / "input/plan.json").read_text())
(a.work / "scan").mkdir(exist_ok=True)
(a.work / "logs").mkdir(exist_ok=True)
rec_path = a.work / "shards.json"
R = json.loads(rec_path.read_text()) if rec_path.exists() else {"runs": [], "shards": {}}
R["runs"].append({"started": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "parallel": a.parallel, "cpus": a.cpus})


def ok(k: int) -> bool:
    j = a.work / f"scan/shard-{k}.json"
    return j.exists() and json.loads(j.read_text()).get("ok", False)


def save() -> None:
    tmp = rec_path.with_suffix(".tmp")
    tmp.write_text(json.dumps(R, indent=1))
    tmp.replace(rec_path)


def argv(k: int) -> list[str]:
    binds = [x for b in a.bind for x in ("--bind", b)]
    return ["nice", "-n", "19", "singularity", "exec", "--cleanenv", "--env", f"PYTHONPATH={a.pythonpath}",
            "--env", f"OMP_NUM_THREADS={a.cpus}", *binds, a.image, "python", a.scan,
            str(a.work / f"input/shard-{k}.parquet"), a.hmm, a.cutoff, str(a.work / f"scan/shard-{k}.parquet"),
            "--cpus", str(a.cpus)]


todo = [s["shard"] for s in plan["shards"] if not ok(s["shard"])]
running: dict[int, tuple[int, float]] = {}
t0 = time.time()
while todo or running:
    while todo and len(running) < a.parallel and not a.stop.exists():
        k = todo.pop(0)
        with open(a.work / f"logs/shard-{k}.log", "a") as lf:
            p = subprocess.Popen(argv(k), stdout=lf, stderr=subprocess.STDOUT)
        running[p.pid] = (k, time.monotonic())
        print(f"shard {k} started pid {p.pid}", flush=True)
    if not running:
        print("STOP present; shards not started: " + ",".join(map(str, todo)), flush=True)
        break
    pid, status, ru = os.wait4(-1, 0)
    if pid not in running:
        continue
    k, ts = running.pop(pid)
    R["shards"][str(k)] = {"wall_s": round(time.monotonic() - ts, 1),
                           "cpu_min": round((ru.ru_utime + ru.ru_stime) / 60, 2),
                           "peak_rss_mb": round(ru.ru_maxrss / 1024, 1), "rc": os.waitstatus_to_exitcode(status)}
    print(f"shard {k} done {json.dumps(R['shards'][str(k)])}", flush=True)
    save()
bad = [s["shard"] for s in plan["shards"] if not ok(s["shard"])]
R["runs"][-1].update({"wall_s": round(time.time() - t0, 1), "not_ok": bad})
save()
print(json.dumps(R["runs"][-1]), flush=True)
sys.exit(1 if bad else 0)
