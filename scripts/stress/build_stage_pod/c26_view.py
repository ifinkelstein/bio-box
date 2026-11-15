"""Card C-26 part 1: the rehearsal root for a sample of the PLSDB stage, with links to the stages its records touch.

    python3 -I c26_view.py SRC DST P1,P2

DST gets real directories holding one symlink per stage directory of SRC, the data root, that a new stage's steps
and its cluster gate read, for every stage named: the l0seq tables contig_index, seq_locator, seq_shards and
genome_input (ingest's anti-join, gene calling's and the CRISPR producer's --seq-stages, dim's ownership check), the
gene calls l1/feature_protein/gene_call_version=<v>/stage=<S> and l0seq/protein_seq/gene_call_version=<v>/stage=<S>
and the run directory runs/stage=<S> (the cluster gate counts each stage's proteins from its gene-call run record,
`stage.release_proteins`), and store_params.json (the pattern of c33_view.py). The new stage's files, ledgers and run
records land in DST's real directories, so nothing is written through a link; the launches bind SRC read-only too.
Refuses when DST exists.
"""
import os
import sys

src, dst, stages = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
os.makedirs(dst)  # refuses an existing DST
n = 0


def link(rel_dir: str, s: str) -> None:
    global n
    p = os.path.join(src, rel_dir, f"stage={s}")
    if not os.path.isdir(p):
        sys.exit(f"{p}: no such directory")
    os.makedirs(os.path.join(dst, rel_dir), exist_ok=True)
    os.symlink(p, os.path.join(dst, rel_dir, f"stage={s}"))
    n += 1


gcv = sorted(os.listdir(os.path.join(src, "l1", "feature_protein")))
if len(gcv) != 1:
    sys.exit(f"gene call versions {gcv}: expected one")
for s in stages:
    for t in ("contig_index", "seq_locator", "seq_shards", "genome_input"):
        link(os.path.join("l0seq", t), s)
    link(os.path.join("l1", "feature_protein", gcv[0]), s)
    link(os.path.join("l0seq", "protein_seq", gcv[0]), s)
    link("runs", s)
os.symlink(os.path.join(src, "l0seq", "store_params.json"), os.path.join(dst, "l0seq", "store_params.json"))
for d in ("ledger", os.path.join("l1", "feature"), "dim", "release_def"):
    os.makedirs(os.path.join(dst, d), exist_ok=True)
print(f"{n} stage directories and store_params.json linked into {dst} (stages {stages}, {gcv[0]})")
