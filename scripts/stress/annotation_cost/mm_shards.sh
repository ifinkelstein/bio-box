#!/bin/sh
# card 0-09a: one MMseqs2 search of the proteins (queries) against the profiles, cut into <n> query shards
# run one after another, so that no process passes 300 CPU-minutes (AGENTS §7.2). Inside the image.
#
#   sh mm_shards.sh <seqDB> <profileDB> <out_dir> <sensitivity> <evalue> <n>
#
# The shards are consecutive ranges of the query DB's keys (`mmseqs createsubdb`, sequences and headers);
# each runs `mm_search.sh ... seq` under timed.py (<out_dir>/shard-<k>.json); <out_dir>/hits.tsv is their
# hits in shard order. A protein's hits do not depend on the other queries, except through the E-value's
# database size, which is the profile DB's and the same in every shard.
set -eu
Q=$1 P=$2 OUT=$3 S=$4 E=$5 N=$6
BIN=$(dirname "$0")
mkdir -p "$OUT"
total=$(wc -l < "$Q.index")
per=$(( (total + N - 1) / N ))
: > "$OUT/hits.tsv"
k=0
while [ "$k" -lt "$N" ]; do
  D="$OUT/shard-$k"
  mkdir -p "$D"
  cut -f1 "$Q.index" | sort -n | awk -v lo=$((k * per)) -v hi=$(((k + 1) * per)) 'NR > lo && NR <= hi' > "$D/keys"
  mmseqs createsubdb "$D/keys" "$Q" "$D/q" -v 1
  mmseqs createsubdb "$D/keys" "${Q}_h" "$D/q_h" -v 1
  python "$BIN/timed.py" "$OUT/shard-$k.json" -- sh "$BIN/mm_search.sh" "$D/q" "$P" "$D" "$S" "$E" seq 300
  cat "$D/hits.tsv" >> "$OUT/hits.tsv"
  k=$((k + 1))
done
wc -l "$OUT/hits.tsv"
