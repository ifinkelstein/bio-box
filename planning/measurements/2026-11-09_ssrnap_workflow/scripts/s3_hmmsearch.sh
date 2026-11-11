#!/bin/sh
# Step 2 (direct path): HMMER's own hmmsearch at the project's fixed search space.
# usage: sh s3_hmmsearch.sh <query.hmm> <tag>
Q="$1"
TAG="$2"
D=/tmp/if2453/ssrnap/out
hmmsearch --cpu 8 -Z 1e8 --domZ 1e8 -T 10 --domT 10 --noali \
  --tblout "$D/$TAG.tbl" --domtblout "$D/$TAG.domtbl" -o /dev/null "$Q" "$D/proteins.faa"
