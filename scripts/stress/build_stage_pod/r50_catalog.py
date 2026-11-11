"""Open a release's catalog with the current code and report what readers find: r50_catalog.py ROOT RELEASE.

50 Gbp build (planning/measurements/2026-11-11_build_r50.md), check A1. Prints one JSON object: the time to open,
`catalog_info`, the views and tables the catalog holds, the row count of the lookup views current readers use
(`contig_index`, `contig_label`, `member_cluster`, `cluster_member`, `token_postings`, `l1_domains`), the columns
of `tape_tokens` and whether `domain_family_tokens` is among them, and the family-level postings. Writes nothing
in the root (DuckDB's spill goes to the system temporary directory, card W-18).
"""
import json
import sys
import time
from pathlib import Path

from biobox import catalog

root, release = Path(sys.argv[1]), sys.argv[2]
t0 = time.monotonic()
con = catalog.open_catalog(root / "catalog" / f"release={release}.duckdb", threads=4)
out = {"release": release, "open_s": round(time.monotonic() - t0, 3)}
out["catalog_info"] = dict(con.execute("SELECT key, value FROM catalog_info ORDER BY key").fetchall())
names = [r[0] for r in con.execute(
    "SELECT table_name FROM information_schema.tables ORDER BY table_name").fetchall()]
out["relations"] = names
counts = {}
for v in ("contig_index", "contig_label", "member_cluster", "cluster_member", "token_postings", "l1_domains",
          "tape_tokens", "token_dict"):
    if v not in names:
        counts[v] = "absent"
        continue
    t = time.monotonic()
    try:
        n = con.execute(f"SELECT count(*) FROM {v}").fetchone()[0]
        counts[v] = {"rows": n, "s": round(time.monotonic() - t, 3)}
    except Exception as e:  # report, do not stop: this is a compatibility probe
        counts[v] = {"error": f"{type(e).__name__}: {e}"[:400]}
out["counts"] = counts
if "tape_tokens" in names:
    cols = [r[0] for r in con.execute(
        "SELECT column_name FROM information_schema.columns WHERE table_name = 'tape_tokens' "
        "ORDER BY ordinal_position").fetchall()]
    out["tape_tokens_columns"] = cols
    out["has_domain_family_tokens"] = "domain_family_tokens" in cols
if "token_postings" in names:
    try:
        cols = [r[0] for r in con.execute(
            "SELECT column_name FROM information_schema.columns WHERE table_name = 'token_postings' "
            "ORDER BY ordinal_position").fetchall()]
        out["token_postings_columns"] = cols
        if "level" in cols:
            out["postings_by_level"] = dict(con.execute(
                "SELECT level, count(*) FROM token_postings GROUP BY level ORDER BY level").fetchall())
    except Exception as e:
        out["postings_by_level"] = f"{type(e).__name__}: {e}"[:400]
print(json.dumps(out, default=str))
