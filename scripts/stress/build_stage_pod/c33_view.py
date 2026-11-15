"""Card C-33: the rehearsal root for the two CRISPR producers on one stage (the pattern of cards E-12, E-15).

    python3 -I c33_view.py SRC DST STAGE S1,P1,P2,P3,P4,P5,P6

DST gets real directories holding one symlink per stage directory of what the producers read from SRC, the data root,
read-only: for every stage named (so that `--seq-stages` can be the real run's, every stage of the root) the l0seq
tables contig_index, seq_locator, seq_shards and genome_input, and store_params.json; for STAGE alone the gene calls
l1/feature_protein/gene_call_version=<v>/stage=<STAGE> and the Pfam scan l1/domains/library=pfam-38.2/stage=<STAGE>.
The producers' own files, ledgers and run records land in DST's real directories (l1/feature, l1/crispr, ledger, runs),
so nothing is written through a link. Refuses when DST exists.
"""
import os
import sys

src, dst, stage, stages = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split(",")
if stage not in stages:
    sys.exit(f"{stage} is not one of {stages}")
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
link(os.path.join("l1", "feature_protein", gcv[0]), stage)
link(os.path.join("l1", "domains", "library=pfam-38.2"), stage)
os.symlink(os.path.join(src, "l0seq", "store_params.json"), os.path.join(dst, "l0seq", "store_params.json"))
for d in ("runs", "ledger", os.path.join("l1", "feature"), os.path.join("l1", "crispr")):
    os.makedirs(os.path.join(dst, d), exist_ok=True)
print(f"{n} stage directories and store_params.json linked into {dst} (sequences of {stages}; gene calls and "
      f"Pfam of {stage}, {gcv[0]})")
