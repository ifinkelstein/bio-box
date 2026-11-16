"""Card C-26: the rehearsal root for a sample of the PLSDB stage, with links to the stages its records touch.

    python3 -I c26_view.py SRC DST P1,P2 [LIBRARY]

DST gets real directories holding one symlink per stage directory of SRC, the data root, that a new stage's steps
and its cluster gate read, for every stage named: the l0seq tables contig_index, seq_locator, seq_shards and
genome_input (ingest's anti-join, gene calling's and the CRISPR producer's --seq-stages, dim's ownership check), the
gene calls l1/feature_protein/gene_call_version=<v>/stage=<S> and l0seq/protein_seq/gene_call_version=<v>/stage=<S>
and the run directory runs/stage=<S> (the cluster gate counts each stage's proteins from its gene-call run record,
`stage.release_proteins`), and store_params.json (the pattern of c33_view.py). The new stage's files, ledgers and run
records land in DST's real directories, so nothing is written through a link; the launches bind SRC read-only too.
Refuses when DST exists.

Part 2 (card B-47b in the rehearsal): with LIBRARY (for example pfam-38.2), also each named stage's completed domain
scan of that library, so that the new stage's scan plan can copy the rows those stages measured (`domains.sources`
reads the scan's run record, its shard plan, its domain files and its ledger): l1/domain_scan/library=<L>/stage=<S>,
l1/domains/library=<L>/stage=<S> and ledger/stage=<S> (the source's ledger rows name its tool versions); and the
root's l3/tape_version_def (one link), so that the new stage's tapes step finds its tape version defined, as in the
data root, and writes no file outside its stage (the first rehearsal of part 2 wrote one without it).
"""
import os
import sys

src, dst, stages = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
library = sys.argv[4] if len(sys.argv) > 4 else None
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
    if library is not None:
        link(os.path.join("l1", "domain_scan", f"library={library}"), s)
        link(os.path.join("l1", "domains", f"library={library}"), s)
        link("ledger", s)
os.symlink(os.path.join(src, "l0seq", "store_params.json"), os.path.join(dst, "l0seq", "store_params.json"))
if library is not None and os.path.isdir(os.path.join(src, "l3", "tape_version_def")):
    # the root's tape version definitions (the stages' tape version is defined there once, 2026-11-10): without them
    # the new stage's tapes step writes l3/tape_version_def/... in the view, a file the data root already holds
    os.makedirs(os.path.join(dst, "l3"), exist_ok=True)
    os.symlink(os.path.join(src, "l3", "tape_version_def"), os.path.join(dst, "l3", "tape_version_def"))
    n += 1
for d in ("ledger", os.path.join("l1", "feature"), "dim", "release_def"):
    os.makedirs(os.path.join(dst, d), exist_ok=True)
print(f"{n} stage directories and store_params.json linked into {dst} (stages {stages}, {gcv[0]}, "
      f"library {library})")
