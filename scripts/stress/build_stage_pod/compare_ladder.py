"""Compare one stage of two releases by rows and coordinates, not bytes (card C-10). Read only.

    python compare_ladder.py NEW_ROOT NEW_RELEASE OLD_ROOT OLD_RELEASE STAGE

Files come from each release's manifest (explicit lists, AGENTS §3.9). Prints JSON: rows per table of
the stage in each manifest; for `feature` the (contig_uid, start, end, strand, kind) rows only in one
side and the matched rows whose feature_uid differs; for `protein_seq` the protein_uids only in one
side and matched uids whose sequence differs; for `igr` the (contig_uid, start, end, strand) rows only
in one side; and the first five genes by sort key from both sides.
"""

import argparse
import json
from collections import Counter
from pathlib import Path

import duckdb

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for n in ("new_root", "new_release", "old_root", "old_release", "stage"):
    ap.add_argument(n)
a = ap.parse_args()


def files(root: str, release: str) -> tuple[dict[str, list[str]], dict[str, int]]:
    man = json.loads((Path(root) / f"manifests/release={release}.json").read_text())
    out: dict[str, list[str]] = {}
    rows: Counter = Counter()
    for f in man["files"]:
        if f["partitions"].get("stage") == a.stage:
            out.setdefault(f["table"], []).append(str(Path(root) / f["path"]))
            rows[f["table"]] += f["rows"]
    return out, dict(sorted(rows.items()))


new, new_rows = files(a.new_root, a.new_release)
old, old_rows = files(a.old_root, a.old_release)
con = duckdb.connect(config={"threads": 4, "memory_limit": "8GB"})


def src(fs: list[str]) -> str:
    return "read_parquet([" + ", ".join("'" + f.replace("'", "''") + "'" for f in fs) + "])"


def one(sql: str):
    return con.execute(sql).fetchone()[0]


def only(t: str, cols: str) -> dict:
    n, o = src(new[t]), src(old[t])
    return {"only_new": one(f"SELECT count(*) FROM (SELECT {cols} FROM {n} EXCEPT ALL SELECT {cols} FROM {o})"),
            "only_old": one(f"SELECT count(*) FROM (SELECT {cols} FROM {o} EXCEPT ALL SELECT {cols} FROM {n})")}


key = "contig_uid, start, \"end\", strand, kind"
res = {"stage": a.stage, "rows_new": new_rows, "rows_old": old_rows,
       "feature": only("feature", key), "igr": only("igr", 'contig_uid, start, "end", strand'),
       "protein_seq": only("protein_seq", "protein_uid")}
res["feature"]["matched_uid_differs"] = one(
    f"SELECT count(*) FROM {src(new['feature'])} n JOIN {src(old['feature'])} o USING ({key}) "
    "WHERE n.feature_uid <> o.feature_uid")
res["protein_seq"]["matched_seq_differs"] = one(
    f"SELECT count(*) FROM {src(new['protein_seq'])} n JOIN {src(old['protein_seq'])} o USING (protein_uid) "
    "WHERE n.seq <> o.seq")
res["first_genes"] = {side: [[r[0].hex(), *r[1:]] for r in con.execute(
    f"SELECT contig_uid, start, \"end\", strand, kind FROM {src(fs['feature'])} ORDER BY {key} LIMIT 5").fetchall()]
    for side, fs in (("new", new), ("old", old))}
print(json.dumps(res, indent=1))
