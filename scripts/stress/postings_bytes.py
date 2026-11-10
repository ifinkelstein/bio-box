"""Bytes per column of `token_postings` for one release, from Parquet footers (card B-26).

    python scripts/stress/postings_bytes.py ROOT RELEASE [--json]

ROOT is the data root, opened read only. File lists come from `ROOT/manifests/release=<RELEASE>.json`
(AGENTS.md section 3.9), never from a glob. For each column: compressed and uncompressed bytes summed
over every row group of every file, the encodings seen, and bytes per row. Rows per token level
(level = token_id >> 61, PLAN.md section 4.1) are counted from the row groups' token_id min/max
statistics; only a row group whose min and max fall in different levels has its token_id column read.
No other data page is read. One thread.
"""
import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import pyarrow.compute as pc
import pyarrow.parquet as pq

LEVEL_NAMES = {1: "family", 2: "nr50", 3: "nr90"}


def measure(root: Path, release: str) -> dict:
    manifest = json.loads((root / f"manifests/release={release}.json").read_text())
    paths = sorted(root / f["path"] for f in manifest["files"] if f["table"] == "token_postings")
    if not paths:
        raise SystemExit(f"no token_postings files in manifests/release={release}.json")
    cols: dict[str, dict] = {}
    rows = row_groups = file_bytes = 0
    levels: Counter = Counter()
    rg_read = 0
    for p in paths:
        pf = pq.ParquetFile(p)
        file_bytes += p.stat().st_size
        md = pf.metadata
        names = [md.row_group(0).column(j).path_in_schema for j in range(md.num_columns)] if md.num_row_groups else []
        for i in range(md.num_row_groups):
            rg = md.row_group(i)
            rows += rg.num_rows
            row_groups += 1
            for j, name in enumerate(names):
                c = rg.column(j)
                d = cols.setdefault(name, {"compressed": 0, "uncompressed": 0, "encodings": set(),
                                           "codec": set(), "dict_row_groups": 0})
                d["compressed"] += c.total_compressed_size
                d["uncompressed"] += c.total_uncompressed_size
                d["encodings"].update(c.encodings)
                d["codec"].add(c.compression)
                d["dict_row_groups"] += int(c.has_dictionary_page)
                if name == "token_id":
                    st = c.statistics
                    if st is not None and st.has_min_max and (st.min >> 61) == (st.max >> 61):
                        levels[st.min >> 61] += rg.num_rows
                    else:
                        t = pf.read_row_group(i, columns=["token_id"]).column(0)
                        rg_read += 1
                        for lv, n in Counter(pc.shift_right(t, 61).to_pylist()).items():
                            levels[lv] += n
    return {"release": release, "files": len(paths), "row_groups": row_groups, "rows": rows,
            "file_bytes": file_bytes, "token_id_row_groups_read": rg_read,
            "rows_by_level": {LEVEL_NAMES.get(k, str(k)): v for k, v in sorted(levels.items())},
            "columns": {k: {"compressed": v["compressed"], "uncompressed": v["uncompressed"],
                            "compressed_per_row": v["compressed"] / rows,
                            "uncompressed_per_row": v["uncompressed"] / rows,
                            "encodings": sorted(v["encodings"]), "codec": sorted(v["codec"]),
                            "row_groups_with_dictionary": v["dict_row_groups"]} for k, v in cols.items()}}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path, help="data root (read only)")
    ap.add_argument("release", help="release name, for example rL.2")
    ap.add_argument("--json", action="store_true", help="print JSON instead of a table")
    a = ap.parse_args(argv)
    r = measure(a.root, a.release)
    if a.json:
        print(json.dumps(r, indent=1, sort_keys=True))
        return 0
    print(f"release {r['release']}: {r['files']} files, {r['row_groups']} row groups, {r['rows']} rows, "
          f"{r['file_bytes']} file bytes ({r['file_bytes'] / r['rows']:.2f} per row)")
    print(f"rows by level: {r['rows_by_level']} (token_id row groups read: {r['token_id_row_groups_read']})")
    print("column\tcompressed\tuncompressed\tcompressed_B_per_row\tuncompressed_B_per_row\tencodings\tcodec\tdict_rgs")
    tot_c = tot_u = 0
    for k, v in r["columns"].items():
        tot_c += v["compressed"]
        tot_u += v["uncompressed"]
        print(f"{k}\t{v['compressed']}\t{v['uncompressed']}\t{v['compressed_per_row']:.3f}\t"
              f"{v['uncompressed_per_row']:.3f}\t{','.join(v['encodings'])}\t{','.join(v['codec'])}\t"
              f"{v['row_groups_with_dictionary']}")
    print(f"TOTAL\t{tot_c}\t{tot_u}\t{tot_c / r['rows']:.3f}\t{tot_u / r['rows']:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
