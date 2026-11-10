"""Print step tables from the ladder result files copied to the laptop.

    python3 summarize.py RESULTS_DIR   (holds A/result.json, B/result.json, ...)
"""

import json
import sys
from pathlib import Path

d = Path(sys.argv[1])
for key in ("A", "B", "Bp", "C"):
    p = d / key / "result.json"
    if not p.exists():
        continue
    R = json.loads(p.read_text())
    print(f"\n## {key} stage={R['stage']} release={R['release']} finished={R.get('finished')} "
          f"failed_hard={R.get('failed_hard_step')} soft={R.get('soft_failures')}")
    print(f"{'step':14s} {'wall_s':>9s} {'thr':>4s} {'cpu_min':>8s} {'peakMB':>9s}  error")
    for name, e in R["steps"].items():
        if name.endswith("_progress"):
            continue
        print(f"{name:14s} {e.get('wall_s', 0):9.1f} {e.get('threads', 1):4d} {e.get('cpu_min', e.get('driver_cpu_s', 0) / 60):8.1f} "
              f"{e.get('peak_rss_mb', 0):9.1f}  {e.get('error') or ''}")
    for name, e in R["steps"].items():
        if e.get("error"):
            print(f"--- {name} traceback:\n{(e.get('traceback') or '')[-2500:]}")
