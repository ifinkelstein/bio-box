"""Where is a contig's sequence stored? r50_dup.py ROOT GENOMES_JSON CONTIG_UID_HEX... (50 Gbp build). Read only.

For each contig uid: the genome entries of GENOMES_JSON (a build-stage `--genomes` list) that name it, and the
rows of every stage's `contig_index` and `contig_occurrence` files under ROOT that hold it (stage, accession,
length). Prints JSON.
"""
import json
import sys
from pathlib import Path

import duckdb

root, genomes, uids = Path(sys.argv[1]), Path(sys.argv[2]), [u.lower() for u in sys.argv[3:]]
con = duckdb.connect(config={"threads": 2})
gs = json.loads(genomes.read_text())
out = {}
for u in uids:
    named = [{k: g[k] for k in g if k != "contigs"} | {"n_contigs": len(g.get("contigs", []))}
             for g in gs if u in json.dumps(g).lower()]
    res = {"genomes_naming_it": named}
    for table in ("contig_index", "contig_occurrence"):
        files = sorted(str(p) for p in (root / "l0seq" / table).glob("stage=*/*.parquet"))
        if not files:
            res[table] = "no files"
            continue
        lst = "[" + ", ".join("'" + f + "'" for f in files) + "]"
        cur = con.execute(f"SELECT * FROM read_parquet({lst}, filename = true) WHERE lower(hex(contig_uid)) = ?", [u])
        cols = [d[0] for d in cur.description]
        res[table] = [{c: (v.hex() if isinstance(v, bytes) else v) for c, v in zip(cols, r)} for r in cur.fetchall()]
    out[u] = res
print(json.dumps(out, indent=1, default=str))
