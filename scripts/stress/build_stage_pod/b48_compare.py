"""Card B-48: a stage's gene-call files in two roots, by sha256, and by rows where the bytes differ (read-only on both).

    python b48_compare.py ROOT_A ROOT_B STAGE OUT_JSON [RECORD_A RECORD_B]

Every file under the stage's gene-call directories (feature of the gene caller, feature_protein, protein_occurrence,
protein_seq, and the ledger of step genecall) in either root: its sha256 in each, and whether they are equal. For a file
in both whose bytes differ: whether its rows are equal (the same table, row for row, through pyarrow), and its rows
counted. With two run records (their gene-call `genecall.json`), whether their outputs (path, bytes, sha256) are the
same. The summary (files, equal, differing, missing on either side) goes to stdout, the rest to OUT_JSON.
"""
import hashlib
import json
import sys
from pathlib import Path

import pyarrow.parquet as pq

from biobox.l1.genecall import GENE_CALL_VERSION as V

a, b, stage, out = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], Path(sys.argv[4])
DIRS = [f"l1/feature/producer={V}/stage={stage}", f"l1/feature_protein/gene_call_version={V}/stage={stage}",
        f"l1/protein_occurrence/gene_call_version={V}/stage={stage}",
        f"l0seq/protein_seq/gene_call_version={V}/stage={stage}", f"ledger/stage={stage}/step=genecall"]


def files(root: Path) -> dict[str, Path]:
    found = {}
    for d in DIRS:
        base = root / d
        if base.is_dir():
            for p in sorted(base.rglob("*")):
                if p.is_file():
                    found[p.relative_to(root).as_posix()] = p
    return found


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while chunk := fh.read(1 << 22):
            h.update(chunk)
    return h.hexdigest()


fa, fb = files(a), files(b)
rows = []
for rel in sorted(set(fa) | set(fb)):
    r = {"path": rel, "a": sha(fa[rel]) if rel in fa else None, "b": sha(fb[rel]) if rel in fb else None}
    r["equal"] = r["a"] is not None and r["a"] == r["b"]
    if r["a"] and r["b"] and not r["equal"]:
        ta, tb = pq.read_table(fa[rel]), pq.read_table(fb[rel])
        r["rows"] = [ta.num_rows, tb.num_rows]
        r["rows_equal"] = ta.equals(tb)
        r["metadata"] = [str(pq.read_metadata(fa[rel]).created_by), str(pq.read_metadata(fb[rel]).created_by)]
    rows.append(r)
res = {"a": str(a), "b": str(b), "stage": stage, "files": len(rows), "equal": sum(r["equal"] for r in rows),
       "differ": [r["path"] for r in rows if r["a"] and r["b"] and not r["equal"]],
       "only_a": [r["path"] for r in rows if r["b"] is None], "only_b": [r["path"] for r in rows if r["a"] is None]}
if len(sys.argv) > 6:
    ra, rb = (json.loads(Path(p).read_text())["outputs"] for p in sys.argv[5:7])
    res["record_outputs_equal"] = ra == rb
out.write_text(json.dumps({**res, "rows": rows}, indent=1))
print(json.dumps(res))
