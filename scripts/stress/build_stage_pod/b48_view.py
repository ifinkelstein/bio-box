"""Card B-48: a scratch root in which a stage of the data root looks ingested and not yet gene-called (read-only on it;
the pattern of cards C-33 and C-34, `c33_view.py`, `c34_view.py`).

    python3 -I b48_view.py SRC DST STAGE

DST gets real directories holding one symlink per directory of what the stage's `inputs` and `ingest` steps wrote in
SRC, the data root (its genome_input, contig_index, contig_occurrence, seq_locator, seq_shards, l0_manifest, the
ingest ledger) and store_params.json, and copies of those two steps' run records under DST/runs/stage=<STAGE>/, so that
`build-stage` passes over them (run-step's resume: each listed output present with its size) and reaches the gene-call
step. Only that stage's sequences are linked: `--seq-stages` is then the stage alone, as when it was first called (a
stage that owns every contig it lists). Nothing of its gene calls is linked: those directories, the gene-call ledger
and every later step's land in DST's real directories, so nothing is written through a link. Refuses when DST exists,
when the stage has no such record, or when one of its outputs is not under a linked directory.
"""
import json
import os
import shutil
import sys

src, dst, stage = sys.argv[1], sys.argv[2], sys.argv[3]
recs = {s: os.path.join(src, "runs", f"stage={stage}", f"{s}.json") for s in ("inputs", "ingest")}
for p in recs.values():
    if not os.path.isfile(p):
        sys.exit(f"{p}: no such run record")
DIRS = [f"l0seq/{t}/stage={stage}" for t in ("genome_input", "contig_index", "contig_occurrence", "seq_locator",
                                             "seq_shards")] + [f"l0/l0_manifest/stage={stage}",
                                                               f"ledger/stage={stage}/step=ingest"]
os.makedirs(dst)  # refuses an existing DST
for d in DIRS:
    if not os.path.isdir(os.path.join(src, d)):
        sys.exit(f"{os.path.join(src, d)}: no such directory")
    os.makedirs(os.path.dirname(os.path.join(dst, d)), exist_ok=True)
    os.symlink(os.path.join(src, d), os.path.join(dst, d))
os.symlink(os.path.join(src, "l0seq", "store_params.json"), os.path.join(dst, "l0seq", "store_params.json"))
os.makedirs(os.path.join(dst, "runs", f"stage={stage}"))
for s, p in recs.items():
    for o in json.loads(open(p).read())["outputs"]:
        if not (o["path"] == "l0seq/store_params.json" or any(o["path"].startswith(d + "/") for d in DIRS)):
            sys.exit(f"{s}: output {o['path']} is not under a linked directory")
    shutil.copy2(p, os.path.join(dst, "runs", f"stage={stage}", f"{s}.json"))
print(f"{len(DIRS)} directories and store_params.json of stage {stage} linked from {src} into {dst}; the records of "
      f"its inputs and ingest steps copied")
