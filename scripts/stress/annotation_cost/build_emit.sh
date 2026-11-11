#!/bin/sh
# Build the MMseqs2 profile database of card 0-09a from Pfam-A.hmm alone (inside the image):
#
#   sh build_emit.sh <Pfam-A.hmm.gz> <out_dir> <n_seqs> <seed>
#
# MMseqs2 18.8cc5c cannot read HMMER3 models: its profile importers are `convertmsa` + `msa2profile`
# (Stockholm MSAs, e.g. Pfam-A.seed or Pfam-A.full, which are not on the pod; fetching them is a download)
# and `convertprofiledb` (HH-suite HHM). So each model's alignment is sampled from the model itself:
#   1. hmmemit -a -N <n_seqs> --seed <seed>: an alignment of n_seqs sequences emitted by each HMM's core
#      model, whose match columns are the model's match states (`pfam_models.py sto` adds the ID and AC
#      lines and writes inserts as '-'-gapped upper case, since MMseqs2 reads the rows as FASTA);
#   2. mmseqs convertmsa (identifier from AC);
#   3. mmseqs msa2profile --match-mode 1 --match-ratio 0.5 (the settings MMseqs2's own Pfam-A.full
#      download uses), default filtering and pseudocounts.
# Each step is timed by timed.py into <out_dir>/build-<step>.json.
set -eu
HMM=$1 OUT=$2 N=$3 SEED=$4
T=${M09A_THREADS:-4}
BIN=$(dirname "$0")
mkdir -p "$OUT"
python "$BIN/timed.py" "$OUT/build-emit.json" -- sh -c \
  "gzip -dc '$HMM' | hmmemit -a -N $N --seed $SEED - | python '$BIN/pfam_models.py' sto '$HMM' > '$OUT/emit.sto'"
python "$BIN/timed.py" "$OUT/build-convertmsa.json" -- mmseqs convertmsa "$OUT/emit.sto" "$OUT/msa" --identifier-field 1
python "$BIN/timed.py" "$OUT/build-msa2profile.json" -- mmseqs msa2profile "$OUT/msa" "$OUT/prof" \
  --match-mode 1 --match-ratio 0.5 --threads "$T"
rm -f "$OUT/emit.sto"
wc -l "$OUT/prof.index"
