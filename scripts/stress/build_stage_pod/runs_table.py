#!/usr/bin/env python3
"""Print the run records of a data root as a Markdown table (card C-10). Host python3, stdlib, read only.

    python3 runs_table.py ROOT

One row per ROOT/runs/stage=<S>/<step>.json, in start order: stage, step, started (UTC), wall s,
CPU-minutes, peak RSS MiB, threads, over budget, mem_cap (how, bytes, cgroup), outputs, image digest.
"""

import json
import sys
from pathlib import Path

if len(sys.argv) != 2 or sys.argv[1] in ("-h", "--help"):
    print(__doc__)
    sys.exit(0 if sys.argv[1:] in (["-h"], ["--help"]) else 2)
recs = [json.loads(p.read_text()) for p in sorted(Path(sys.argv[1]).glob("runs/stage=*/*.json"))
        if not p.name.endswith(".publishing.json")]
print("| stage | step | started (UTC) | wall s | CPU-min | peak RSS MiB | threads | over budget | mem_cap | outputs |")
print("|---|---|---|---|---|---|---|---|---|---|")
for r in sorted(recs, key=lambda r: r["started_at"]):
    c = r.get("mem_cap", {})
    cap = c.get("how", "?") + (f" {c['memory_max_bytes']} B at {c.get('limited_at')}" if "memory_max_bytes" in c else "")
    print(f"| {r['stage']} | {r['step']} | {r['started_at']} | {r['wall_s']:.1f} | {r['cpu_min']:.2f} | "
          f"{r['peak_rss_bytes'] / 2**20:,.0f} | {r['threads']} | {r['over_cpu_budget']} | {cap} | {len(r['outputs'])} |")
digests = sorted({r.get("image_digest") for r in recs})
print(f"\n{len(recs)} records; image digests: {', '.join(map(str, digests))}; hosts: "
      f"{', '.join(sorted({r['host'] for r in recs}))}")
