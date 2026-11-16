"""Card B-49: check a release built by reuse in a scratch view, read-only.

    python b49_check.py VIEW NEW SOURCE

Prints one JSON document: for every table file the source's run record of its step `cluster` lists, the sha256 of the
new release's file of the same name (computed here, a block at a time) against the record's and whether it is a real
file of the view (not a link); the new release's run records (step, cpu_min, wall_s, peak_rss_bytes, argv's source);
its clustering ledger; its reuse record's names and hashes; and the bytes of the view's real files (links not
counted). Writes nothing.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

import pyarrow.parquet as pq


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while b := fh.read(1 << 23):
            h.update(b)
    return h.hexdigest()


def main(view: str, new: str, source: str) -> None:
    view = Path(view)
    rec = json.loads((view / f"runs/stage={source}/cluster.json").read_text())
    files = []
    for o in rec["outputs"]:
        if not o["path"].startswith("l2/"):
            continue
        p = view / o["path"].replace(f"/release={source}/", f"/release={new}/", 1)
        got = sha256(p) if p.exists() else None
        files.append({"path": str(p.relative_to(view)), "bytes": o["bytes"], "sha256_equal": got == o["sha256"],
                      "real_file": p.is_file() and not p.is_symlink(), "nlink": p.stat().st_nlink if p.exists() else None})
    runs = {}
    for r in sorted((view / f"runs/stage={new}").glob("*.json")):
        d = json.loads(r.read_text())
        argv = d.get("argv", [])
        runs[d["step"]] = {"cpu_min": d["cpu_min"], "wall_s": d["wall_s"], "peak_rss_bytes": d["peak_rss_bytes"],
                           "from": argv[argv.index("--from") + 1] if "--from" in argv else None,
                           "outputs": len(d["outputs"])}
    led = [{k: (v.decode() if isinstance(v, bytes) else v) for k, v in row.items()} for row in
           pq.read_table(view / f"ledger/release={new}/step=cluster/part-0.parquet").to_pylist()]
    reuse = json.loads((view / f"l2/cluster_plan/release={new}/cluster.reuse.json").read_text())
    later = {}  # the later steps' files against the source release's same steps, by their run records' sha256
    for step in ("tokens", "postings", "dim"):
        mine, theirs = (view / f"runs/stage={new}/{step}.json"), (view / f"runs/stage={source}/{step}.json")
        if not mine.exists() or not theirs.exists():
            continue
        want = {o["path"]: o["sha256"] for o in json.loads(theirs.read_text())["outputs"]}
        got = {o["path"].replace(f"release={new}", f"release={source}"): o["sha256"]
               for o in json.loads(mine.read_text())["outputs"]}
        later[step] = {"files": len(got), "same_paths": sorted(got) == sorted(want),
                       "equal": sum(got[p] == want.get(p) for p in got)}
    source_cost = {"records": 0, "cpu_min": 0.0}  # the source's clustering: every cluster* run record of its release
    for r in sorted((view / f"runs/stage={source}").glob("cluster*.json")):
        d = json.loads(r.read_text())
        source_cost["records"] += 1
        source_cost["cpu_min"] += d["cpu_min"]
    source_cost["cpu_min"] = round(source_cost["cpu_min"], 3)
    real = 0
    for d, _, fs in os.walk(view):
        for f in fs:
            p = os.path.join(d, f)
            if not os.path.islink(p):
                real += os.path.getsize(p)
    print(json.dumps({"files": files, "all_equal": all(f["sha256_equal"] and f["real_file"] for f in files),
                      "copied_bytes": sum(f["bytes"] for f in files), "runs": runs, "ledger": led,
                      "reuse": {k: reuse[k] for k in ("from", "candidates", "plan_sha256", "from_plan_sha256",
                                                       "inputs_sha256", "from_inputs_sha256")},
                      "records": len(reuse["records"]), "view_real_bytes": real, "source_cost": source_cost,
                      "later_steps_against_source": later},
                     indent=1, sort_keys=True))


if __name__ == "__main__":
    main(*sys.argv[1:])
