"""Print the cost table of card 0-09a from the run records (host python3, reads only).

    python3 costs.py <m09a_dir> <root> <stage>

<m09a_dir>/runs/m09a-*.json and <m09a_dir>/data/*/build-*.json are timed.py records; the stage's hmmsearch
cost is the sum of `cpu_min` over the producer's shard run records
`<root>/runs/stage=<S>/domains-pfam-38.2.shard-*.json` (plus the plan and merge records, listed apart).
"""

from __future__ import annotations

import glob
import json
import sys
from pathlib import Path


def main() -> int:
    m, root, stage = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
    shards = sorted(glob.glob(str(root / f"runs/stage={stage}/domains-pfam-38.2.shard-*.json")))
    cpu = [json.load(open(f))["cpu_min"] for f in shards]
    walls = [json.load(open(f))["wall_s"] for f in shards]
    other = {}
    for name in ("plan", ""):
        f = root / f"runs/stage={stage}/domains-pfam-38.2{'.' + name if name else ''}.json"
        if f.exists():
            other[name or "merge"] = json.load(open(f))["cpu_min"]
    print(f"hmmsearch {stage}: {len(cpu)} shards, sum cpu_min {sum(cpu):.2f}, max shard {max(cpu):.2f}, "
          f"shard wall {min(walls):.0f}-{max(walls):.0f} s; plan/merge cpu_min {other}")
    print("| record | wall s | CPU-min | peak RSS MiB | exit |")
    print("|---|---|---|---|---|")
    recs = sorted(glob.glob(str(m / "runs/m09a-*.json"))) + sorted(glob.glob(str(m / "data/*/build-*.json")))
    for f in recs:
        r = json.load(open(f))
        print(f"| {Path(f).relative_to(m)} | {r['wall_s']:.1f} | {r['cpu_min']:.2f} | {r['peak_rss_mib']:.0f} | "
              f"{r['exit_code']} |")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
