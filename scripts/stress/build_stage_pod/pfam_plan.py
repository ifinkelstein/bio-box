"""Cut a stage's distinct proteins into scan shards (card C-10; the ladder's pfam_plan, step.py).

    python pfam_plan.py ROOT STAGE GENE_CALL_VERSION WORK [--residues 16000000]

Reads the `protein_seq` files that the stage's gene-call run record lists (runs/stage=<S>/genecall.json,
an explicit list), which are the proteins `bio-box stage domains` scans (`domains.load_proteins`). Sorts
them by protein_uid and cuts them into shards of about --residues residues (without '*', as
`domains.scan` counts). Writes WORK/input/shard-<k>.parquet (protein_uid, seq) and WORK/input/plan.json
and prints the plan without its shard list. Refuses when WORK/input exists. Writes nothing in ROOT.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from biobox import ids

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("root", type=Path)
ap.add_argument("stage")
ap.add_argument("gene_call_version")
ap.add_argument("work", type=Path)
ap.add_argument("--residues", type=int, default=16_000_000)
a = ap.parse_args()

rec = json.loads((a.root / f"runs/stage={a.stage}/genecall.json").read_text())
prefix = f"l0seq/protein_seq/gene_call_version={a.gene_call_version}/stage={a.stage}/"
files = sorted(o["path"] for o in rec["outputs"] if o["path"].startswith(prefix))
inp = a.work / "input"
inp.mkdir(parents=True)  # refuses an existing plan
t = pa.concat_tables([pq.read_table(a.root / f, columns=["protein_uid", "seq"]) for f in files])
t = t.cast(pa.schema([("protein_uid", pa.binary(16)), ("seq", pa.string())]))
t = t.take(pc.sort_indices(t, sort_keys=[("protein_uid", "ascending")]))
if pc.count_distinct(t.column("protein_uid")).as_py() != t.num_rows:
    raise SystemExit("protein_uid repeats across protein_seq files")
lens = pc.subtract(pc.utf8_length(t.column("seq")), pc.count_substring(t.column("seq"), "*")).to_numpy()
cum = np.cumsum(lens)
shards, k, start = [], 0, 0
while start < t.num_rows:
    base = int(cum[start - 1]) if start else 0
    end = min(max(int(np.searchsorted(cum, base + a.residues, side="right")), start + 1), t.num_rows)
    out = inp / f"shard-{k}.parquet"
    pq.write_table(t.slice(start, end - start), out, compression="zstd")
    shards.append({"shard": k, "proteins": end - start, "residues": int(cum[end - 1]) - base,
                   "first_uid": t.column("protein_uid")[start].as_py().hex(),
                   "sha256": ids.sha256_file(out).hex()})
    k, start = k + 1, end
plan = {"root": str(a.root), "stage": a.stage, "gene_call_version": a.gene_call_version,
        "protein_seq_files": files, "proteins": t.num_rows, "residues": int(cum[-1]) if len(cum) else 0,
        "residues_per_shard": a.residues, "shards": shards}
(inp / "plan.json").write_text(json.dumps(plan, indent=1))
print(json.dumps({k: v for k, v in plan.items() if k not in ("shards", "protein_seq_files")} | {"n_shards": len(shards)}))
