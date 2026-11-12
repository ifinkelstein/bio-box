"""Compare the clusters of two releases row by row, not by bytes (50 Gbp build, checks A1 and A3). Read only.

    python r50_compare_clusters.py NEW_ROOT NEW_RELEASE OLD_ROOT OLD_RELEASE

r10_compare_clusters.py (branch build/r10) with the members that move. Files are the `cluster_member` files of
each release's own partition (`release=<R>`) in its manifest (explicit lists, AGENTS §3.9). Prints JSON: per
level, rows and distinct clusters on each side; the (member_uid, level, cluster_uid, rep_uid) rows only in one
side (`EXCEPT ALL`); and per level the members whose cluster or representative differs, or that are on one side
only, with up to five examples.
"""

import argparse
import json
from pathlib import Path

import duckdb

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for n in ("new_root", "new_release", "old_root", "old_release"):
    ap.add_argument(n)
a = ap.parse_args()


def files(root: str, release: str) -> list[str]:
    man = json.loads((Path(root) / f"manifests/release={release}.json").read_text())
    return [str(Path(root) / f["path"]) for f in man["files"]
            if f["table"] == "cluster_member" and f["partitions"].get("release") == release]


def src(fs: list[str]) -> str:
    return "read_parquet([" + ", ".join("'" + f.replace("'", "''") + "'" for f in fs) + "])"


new, old = src(files(a.new_root, a.new_release)), src(files(a.old_root, a.old_release))
con = duckdb.connect(config={"threads": 4, "memory_limit": "16GB"})
cols = "member_uid, level, cluster_uid, rep_uid"
res = {"new": {"release": a.new_release, "files": len(files(a.new_root, a.new_release))},
       "old": {"release": a.old_release, "files": len(files(a.old_root, a.old_release))}}
for side, s in (("new", new), ("old", old)):
    res[side]["levels"] = {lv: {"rows": n, "clusters": c} for lv, n, c in con.execute(
        f"SELECT level, count(*), count(DISTINCT cluster_uid) FROM {s} GROUP BY level ORDER BY level").fetchall()}
res["only_new"] = con.execute(f"SELECT count(*) FROM (SELECT {cols} FROM {new} EXCEPT ALL SELECT {cols} FROM {old})").fetchone()[0]
res["only_old"] = con.execute(f"SELECT count(*) FROM (SELECT {cols} FROM {old} EXCEPT ALL SELECT {cols} FROM {new})").fetchone()[0]
if res["only_new"] or res["only_old"]:
    moved = {}
    for (lv,) in con.execute(f"SELECT DISTINCT level FROM {new} UNION SELECT DISTINCT level FROM {old} ORDER BY 1").fetchall():
        q = (f"SELECT lower(hex(coalesce(n.member_uid, o.member_uid))) AS member, lower(hex(o.cluster_uid)) AS old_cluster, "
             f"lower(hex(n.cluster_uid)) AS new_cluster, lower(hex(o.rep_uid)) AS old_rep, lower(hex(n.rep_uid)) AS new_rep "
             f"FROM (SELECT * FROM {new} WHERE level = ?) n FULL OUTER JOIN (SELECT * FROM {old} WHERE level = ?) o "
             f"ON n.member_uid = o.member_uid WHERE n.cluster_uid IS DISTINCT FROM o.cluster_uid OR n.rep_uid IS DISTINCT FROM o.rep_uid")
        rows = con.execute(f"SELECT count(*) FROM ({q})", [lv, lv]).fetchone()[0]
        ex = con.execute(f"{q} ORDER BY member LIMIT 5", [lv, lv]).fetchall()
        moved[lv] = {"members": rows, "examples": [dict(zip(("member", "old_cluster", "new_cluster", "old_rep", "new_rep"), e)) for e in ex]}
    res["moved"] = moved
print(json.dumps(res, indent=1, default=str))
