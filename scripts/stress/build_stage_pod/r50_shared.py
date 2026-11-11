"""Contigs of later stages whose sequence an earlier stage stores: r50_shared.py ROOT STAGE=GENOMES_JSON... Read only.

50 Gbp build (planning/measurements/2026-11-11_build_r50.md). Stages are given in build order. For each, the
contig uids its `--genomes` list names (a contig uid is the hash of the sequence, card B-28) that are already in
the `contig_index` of a stage built before it in the root, or named by an earlier stage of the list. Such a
contig is not stored again by ingest, so gene calling must read the earlier stage (`genecall --seq-stages`).
Prints JSON: per stage, contigs named, distinct, and the shared ones with the stage that stores each first.
"""
import json
import sys
from pathlib import Path

import duckdb

root = Path(sys.argv[1])
con = duckdb.connect(config={"threads": 2})
first: dict[str, str] = {}
for f in sorted((root / "l0seq" / "contig_index").glob("stage=*/*.parquet")):
    st = f.parent.name.removeprefix("stage=")
    for (u,) in con.execute(f"SELECT lower(hex(contig_uid)) FROM read_parquet('{f}')").fetchall():
        first.setdefault(u, st)
in_root = dict(first)
out = {"root_contigs": len(in_root), "stages": {}}
for arg in sys.argv[2:]:
    stage, path = arg.split("=", 1)
    named = [u.lower() for g in json.loads(Path(path).read_text()) for u in g["contig_uids"]]
    distinct = sorted(set(named))
    shared = {u: first[u] for u in distinct if u in first and first[u] != stage}
    by: dict[str, int] = {}
    for s in shared.values():
        by[s] = by.get(s, 0) + 1
    out["stages"][stage] = {"named": len(named), "distinct": len(distinct), "shared": len(shared), "by_stage": by,
                            "examples": dict(list(shared.items())[:5])}
    for u in distinct:
        first.setdefault(u, stage)
print(json.dumps(out, indent=1))
