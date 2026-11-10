"""Do the catalog's sequence-lookup macros prune row groups? 1e6-contig index (9 row groups).

    python macro_scan.py LOCATOR_ROOT
"""
import json
import sys
import time
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from biobox import catalog, io

if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
    print(__doc__)
    sys.exit(0 if sys.argv[1:2] in (["-h"], ["--help"]) else 2)

root = Path(sys.argv[1])
man = catalog.build_manifest(root, "rX")
cat = catalog.generate_catalog(man, root.parent / "locator.duckdb")
t = time.perf_counter()
con = catalog.open_catalog(cat)
open_s = round(time.perf_counter() - t, 2)
idx = pa.concat_tables(pq.read_table(f, columns=["contig_uid"]) for f in io.part_files(root / "l0seq/contig_index/stage=S0"))
uids = idx.column("contig_uid").to_pylist()
out = {"open_catalog_s": open_s, "contig_index_row_groups": man["files"][0]["row_groups"] if man["files"][0]["table"] == "contig_index" else None}


def run(name, sql, params=None):
    t = time.perf_counter()
    r = catalog.profile_query(con, sql, params)
    out[name] = {"seconds": round(time.perf_counter() - t, 3),
                 "scans": [f"{Path(s['files'][0]).parts[-3] if 'contig_index' in s['files'][0] else 'seq_locator'}: row groups "
                           f"{s['row_groups_scanned']}/{s['row_groups_total']}, rows {s['rows_scanned']}" for s in r["scans"] if s["files"]]}


one = uids[500_000]
run("literal: default_digest(?) + digest_len", "SELECT digest_len(default_digest(?))", [one])
con.register("g1", pa.table({"contig_uid": pa.array([one], pa.binary())}))
run("join, 1 contig: digest_len(default_digest(contig_uid)) FROM g", "SELECT sum(digest_len(default_digest(contig_uid))) FROM g1")
con.register("g1000", pa.table({"contig_uid": pa.array(uids[:: len(uids) // 1000][:1000], pa.binary())}))
run("join, 1,000 contigs", "SELECT sum(digest_len(default_digest(contig_uid))) FROM g1000")
print(json.dumps(out, indent=1))
