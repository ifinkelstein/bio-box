"""Card B-49: compare two releases' clustering records, read-only.

    python b49_compare_plans.py ROOT R1 R2

Prints, as one JSON document: for every `l2/cluster_plan/release=<R>/*.json` record of the two releases, whether the
two are equal once `release` and `plan_sha256` (the release's name, and the hash of the plan record that covers it)
are set aside, and the keys that differ otherwise; the names present in one release only; and both releases'
clustering ledger rows (`ledger/release=<R>/step=cluster`). Writes nothing.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pyarrow
import pyarrow.parquet as pq

NAMED = ("release", "plan_sha256")  # the fields that name the release


def records(root: Path, release: str) -> dict[str, dict]:
    d = root / "l2/cluster_plan" / f"release={release}"
    return {p.name: json.loads(p.read_bytes()) for p in sorted(d.glob("*.json"))}


def main(root: str, r1: str, r2: str) -> None:
    root = Path(root)
    a, b = records(root, r1), records(root, r2)
    differ, named = {}, {}
    for name in sorted(set(a) & set(b)):
        x, y = a[name], b[name]
        keys = sorted(k for k in set(x) | set(y) if x.get(k) != y.get(k))
        named[name] = [k for k in keys if k in NAMED]
        if [k for k in keys if k not in NAMED]:
            differ[name] = [k for k in keys if k not in NAMED]
    ledgers = {}
    for r in (r1, r2):
        p = root / f"ledger/release={r}/step=cluster/part-0.parquet"
        ledgers[r] = ([{k: (v.decode() if isinstance(v, bytes) else v) for k, v in row.items()}
                       for row in pq.read_table(p).to_pylist()] if p.exists() else None)
    named_sets = sorted({tuple(v) for v in named.values()})
    print(json.dumps({"records": {r1: len(a), r2: len(b)}, "common": len(set(a) & set(b)),
                      "only": {r1: sorted(set(a) - set(b)), r2: sorted(set(b) - set(a))},
                      "differ_beyond_the_name": differ, "named_fields_that_differ": named_sets,
                      "ledgers": ledgers, "pyarrow_here": pyarrow.__version__}, indent=1, sort_keys=True))


if __name__ == "__main__":
    main(*sys.argv[1:])
