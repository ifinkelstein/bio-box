#!/bin/sh
# One MMseqs2 profile search of card 0-09a (inside the image, under timed.py):
#
#   sh mm_search.sh <seqDB> <profileDB> <out_dir> <sensitivity> <evalue> [seq|prof] [max_seqs]
#
# seq (default): proteins are the queries and the Pfam profiles the target (`mmseqs search seqDB profDB`).
# prof: the profiles are the queries and the proteins the target (`mmseqs search profDB seqDB`); then
#   max_seqs (default 300, MMseqs2's) caps the prefilter hits per profile, and E-values are relative to the
#   protein DB's size.
# Either way `convertalis` writes <out_dir>/hits.tsv with the protein first: protein uid, Pfam accession,
# protein start, end (1-based inclusive), profile start, end, evalue, bits. Thread count from
# M09A_THREADS (launch.sh sets it). The search's temporary directory is removed at the end.
set -eu
Q=$1 P=$2 OUT=$3 S=$4 E=$5 DIR=${6:-seq} MAXS=${7:-300}
T=${M09A_THREADS:-4}
mkdir -p "$OUT"
if [ "$DIR" = seq ]; then
  # M09A_EXTRA: extra search options (the diagnostic run uses --exhaustive-search)
  # shellcheck disable=SC2086
  mmseqs search "$Q" "$P" "$OUT/aln" "$OUT/tmp" -s "$S" -e "$E" --max-seqs "$MAXS" --threads "$T" ${M09A_EXTRA:-}
  mmseqs convertalis "$Q" "$P" "$OUT/aln" "$OUT/hits.tsv" \
    --format-output query,target,qstart,qend,tstart,tend,evalue,bits --threads "$T"
else
  mmseqs search "$P" "$Q" "$OUT/aln" "$OUT/tmp" -s "$S" -e "$E" --max-seqs "$MAXS" --threads "$T"
  mmseqs convertalis "$P" "$Q" "$OUT/aln" "$OUT/hits.tsv" \
    --format-output target,query,tstart,tend,qstart,qend,evalue,bits --threads "$T"
fi
rm -rf "$OUT/tmp"
wc -l "$OUT/hits.tsv"
