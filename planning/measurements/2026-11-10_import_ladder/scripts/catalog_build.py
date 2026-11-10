"""Manifest, catalog and verification for one store root.

    python catalog_build.py ROOT RELEASE OUT_JSON [--build-dir LOCAL_DIR]

build_manifest -> write_manifest(<root>/manifests/release=<R>.json) -> generate_catalog
(<root>/catalog/release=<R>.duckdb) -> verify_manifest, plus a reopen of the catalog with
open_catalog and row counts per view. With --build-dir (for a root on NFS), the DuckDB file is
generated on local disk under LOCAL_DIR and then copied into <root>/catalog/ (tmp name + rename),
so DuckDB never writes a database file over NFS; the views still point at ROOT.
"""

import argparse
import json
import os
import shutil
import time
from collections import defaultdict
from pathlib import Path

from biobox import catalog
from biobox.registry.schemas import TABLES

ap = argparse.ArgumentParser()
ap.add_argument("root", type=Path)
ap.add_argument("release")
ap.add_argument("out", type=Path)
ap.add_argument("--build-dir", type=Path)
ap.add_argument("--cache-dir", type=Path, default=Path("/tmp/if2453/ladder/catalog_cache/scratch"))
a = ap.parse_args()

R = {"root": str(a.root), "release": a.release, "times_s": {}}


def timed(name, fn):
    t = time.perf_counter()
    out = fn()
    R["times_s"][name] = round(time.perf_counter() - t, 2)
    return out


man = timed("build_manifest", lambda: catalog.build_manifest(a.root, a.release))
mpath = a.root / f"manifests/release={a.release}.json"
timed("write_manifest", lambda: catalog.write_manifest(man, mpath))
R["manifest"] = str(mpath)
R["manifest_hash"] = catalog.manifest_hash(man)
R["manifest_root"] = man["root"]
per = defaultdict(lambda: {"files": 0, "rows": 0, "bytes": 0})
for f in man["files"]:
    per[f["table"]]["files"] += 1
    per[f["table"]]["rows"] += f["rows"]
    per[f["table"]]["bytes"] += f["bytes"]
R["manifest_tables"] = dict(sorted(per.items()))
R["manifest_artifacts"] = len(man["artifacts"])
cpath = a.root / f"catalog/release={a.release}.duckdb"
if a.build_dir:
    a.build_dir.mkdir(parents=True, exist_ok=True)
    local = a.build_dir / cpath.name
    timed("generate_catalog", lambda: catalog.generate_catalog(man, local))
    cpath.parent.mkdir(parents=True, exist_ok=True)
    tmp = cpath.with_name(f".{cpath.name}.tmp-{os.getpid()}")
    timed("copy_catalog_to_root", lambda: shutil.copyfile(local, tmp))
    os.replace(tmp, cpath)
else:
    timed("generate_catalog", lambda: catalog.generate_catalog(man, cpath))
R["catalog"] = str(cpath)
R["catalog_bytes"] = cpath.stat().st_size
problems = timed("verify_manifest", lambda: catalog.verify_manifest(catalog.read_manifest(mpath), a.root))
R["verify_problems"] = problems
con = timed("open_catalog", lambda: catalog.open_catalog(cpath, local_dir=a.cache_dir, threads=4))
R["catalog_info"] = dict(con.execute("SELECT key, value FROM catalog_info").fetchall())
counts = {}
t = time.perf_counter()
for name in sorted(TABLES):
    try:
        counts[name] = con.execute(f"SELECT count(*) FROM {name}").fetchall()[0][0]
    except Exception as e:  # noqa: BLE001
        counts[name] = f"{type(e).__name__}: {e}"
R["times_s"]["count_all_views"] = round(time.perf_counter() - t, 2)
R["view_counts"] = counts
mismatch = {k: (counts.get(k), v["rows"]) for k, v in per.items() if counts.get(k) != v["rows"]}
R["view_vs_manifest_row_mismatch"] = mismatch
a.out.write_text(json.dumps(R, indent=1))
print(json.dumps(R))
raise SystemExit(1 if problems or mismatch else 0)
