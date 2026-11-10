"""Bytes on disk per table and per column, from the Parquet footers of a release's manifest.

    python column_bytes.py MANIFEST_JSON [--root ROOT] [--where KEY=VALUE ...] > bytes.json

Reads only file footers (no data pages). For every table in the manifest: files, rows, file bytes and,
per column, the compressed bytes of its column chunks summed over row groups, with bytes per row.
`--where stage=P1` keeps the files whose partition `stage` is `P1` (repeatable; all must match), so
one stage of a multi-stage store can be measured alone. Column chunks compress independently, so the
bytes a table would have without some of its columns are its file bytes minus those columns' chunks
(plus a small footer change): card B-18 uses this to correct per-unit sizes measured before card B-17
removed the legacy columns. `--root` overrides the manifest's root (another host's path).
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("manifest", type=Path)
ap.add_argument("--root", type=Path)
ap.add_argument("--where", action="append", default=[], help="KEY=VALUE partition filter")
a = ap.parse_args()

import pyarrow.parquet as pq  # noqa: E402

man = json.loads(a.manifest.read_text())
root = a.root or Path(man["root"])
want = dict(w.split("=", 1) for w in a.where)
out: dict = defaultdict(lambda: {"files": 0, "rows": 0, "bytes": 0, "columns": defaultdict(int)})
for f in man["files"]:
    if any(f["partitions"].get(k) != v for k, v in want.items()):
        continue
    t = out[f["table"]]
    t["files"] += 1
    t["rows"] += f["rows"]
    t["bytes"] += f["bytes"]
    md = pq.ParquetFile(root / f["path"]).metadata
    for g in range(md.num_row_groups):
        rg = md.row_group(g)
        for c in range(rg.num_columns):
            col = rg.column(c)
            t["columns"][col.path_in_schema] += col.total_compressed_size
res = {}
for name, t in sorted(out.items()):
    n = t["rows"] or 1
    res[name] = {"files": t["files"], "rows": t["rows"], "bytes": t["bytes"], "bytes_per_row": round(t["bytes"] / n, 3),
                 "columns": {c: {"bytes": b, "bytes_per_row": round(b / n, 3)} for c, b in sorted(t["columns"].items())}}
json.dump({"manifest": str(a.manifest), "release": man.get("release_id"), "where": want, "tables": res}, sys.stdout, indent=1)
print()
