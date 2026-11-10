"""Merge the scan shards of a pfam_plan.py plan into the files `bio-box stage domains` writes (card C-10).

    python pfam_merge.py WORK LIBRARY VERSION [--out DIR]

Run as the step command of `bio-box run-step --stage <S> --step domains-<LIBRARY>-<VERSION>`: it writes
under --out, else $BIOBOX_STAGING, never into the plan's root. Checks that the plan lists exactly the
`protein_seq` files of the stage's gene-call run record and that every shard's scan is ok with the
plan's protein count and input sha256. The shards are disjoint runs of proteins sorted by uid and each
shard's rows are sorted as `domains.scan` sorts them, so their concatenation in shard order is the row
order one scan of all proteins gives (checked). Then it writes the rows and the ledger exactly as
`domains.run` does: `io.write_parts` of `ROW_GROUP_ROWS` batches under
l1/domains/library=<LIBRARY>-<VERSION>/stage=<S>/, and the two ledger rows with n_in = the stage's
proteins. Prints a JSON report.
"""

import argparse
import json
import os
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from biobox import ids, layout, ledger, provenance
from biobox import io as bio
from biobox.l1 import domains
from biobox.registry.schemas import table

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("work", type=Path)
ap.add_argument("library")
ap.add_argument("version")
ap.add_argument("--out", type=Path)
a = ap.parse_args()
out = a.out or (Path(os.environ["BIOBOX_STAGING"]) if os.environ.get("BIOBOX_STAGING") else None)
if out is None:
    sys.exit("pfam_merge: no --out and no $BIOBOX_STAGING; refusing to write into the root")

plan = json.loads((a.work / "input/plan.json").read_text())
root, stage, gcv = Path(plan["root"]), plan["stage"], plan["gene_call_version"]
rec = json.loads((root / f"runs/stage={stage}/genecall.json").read_text())
prefix = f"l0seq/protein_seq/gene_call_version={gcv}/stage={stage}/"
if sorted(o["path"] for o in rec["outputs"] if o["path"].startswith(prefix)) != plan["protein_seq_files"]:
    sys.exit("pfam_merge: the plan's protein_seq files are not the gene-call record's")
lib_dir = f"l1/domains/library={a.library}-{a.version}/stage={stage}"
step = f"domains-{a.library}-{a.version}"
ledger.refuse_written(root, out, [lib_dir, f"ledger/stage={stage}/step={step}"], f"{step} of stage {stage}")

spec = table("l1_domains")
parts, metas = [], []
for s in plan["shards"]:
    k = s["shard"]
    m = json.loads((a.work / f"scan/shard-{k}.json").read_text())
    if not m.get("ok") or m["proteins"] != s["proteins"]:
        sys.exit(f"pfam_merge: shard {k} not ok or its protein count differs: {m.get('error')}")
    if ids.sha256_file(a.work / f"input/shard-{k}.parquet").hex() != s["sha256"]:
        sys.exit(f"pfam_merge: shard {k} input changed since the plan")
    metas.append(m)
    parts.append(pq.read_table(a.work / f"scan/shard-{k}.parquet").cast(spec.schema))
t = pa.concat_tables(parts).combine_chunks()
order = pc.sort_indices(t, sort_keys=[("protein_uid", "ascending"), ("env_start", "ascending"),
                                      ("model_acc", "ascending")])
if order.to_pylist() != list(range(t.num_rows)):
    sys.exit("pfam_merge: the concatenated shards are not in domains.scan order")
n = bio.ROW_GROUP_ROWS
bio.write_parts([t.slice(i, n).to_batches()[0] for i in range(0, t.num_rows, n)],  # domains._batches' cuts
                spec, out / lib_dir, domains.EST_ROW_BYTES, layout.MAX_FILE_BYTES)
hit = pc.count_distinct(t.column("protein_uid")).as_py() if t.num_rows else 0
lib = f"{a.library}-{a.version}".encode()
tools = provenance.tool_versions("pyhmmer", "pyarrow")
proteins = plan["proteins"]
ledger.write(out, stage, step,
             sorted([ledger.row(lib, "library", proteins, hit, tools=tools),
                     ledger.row(lib, "library", proteins, hit, n_rejected=proteins - hit, reason="no_hit", tools=tools)],
                    key=lambda r: (r["object_key"], r["reason_code"])))
print(json.dumps({"stage": stage, "step": step, "proteins": proteins, "domains": t.num_rows, "proteins_with_hit": hit,
                  "models_with_hit": pc.count_distinct(t.column("model_acc")).as_py() if t.num_rows else 0,
                  "shards": len(metas), "shard_cpu_min": round(sum(m["cpu_min"] for m in metas), 2),
                  "shard_peak_rss_mb_max": max(m["peak_rss_mb"] for m in metas)}))
