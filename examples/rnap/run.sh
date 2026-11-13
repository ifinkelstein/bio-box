#!/bin/bash
# The RNA polymerase example as a release check (examples/rnap/README.md; GUIDE.md sections 3 to 7 teach the same
# commands). Run ON THE POD HOST (not in the image), from any directory:
#
#   TREE=... PROFILE=... IMAGE=... DIGEST=sha256:... ROOT=... RELEASE=rN.3 OUT=/tmp/if2453/<you>/rnap-rN.3 \
#     bash "$TREE/examples/rnap/run.sh"
#
# Required: TREE (a checkout of the code, also the image's PYTHONPATH), PROFILE (host profile), IMAGE (the .sif,
# with matplotlib), DIGEST (sha256: of the image), ROOT (data root, bound read-only), RELEASE, OUT (a new or empty
# scratch directory outside ROOT). Optional: PFAM (the Pfam-A.hmm.gz to fetch the model from), THREADS (2), CAP
# (the account's thread cap; default 40 Mon-Fri 08:00-20:00, else 56), BIND (the directory bound read-write into
# the image, holding TREE, OUT, the pack directory and the image; default /tmp/if2453), UNIT (unit name prefix;
# default rnap), MEM (pod_run.sh --mem; default 16G).
#
# What it does, each bio-box step through scripts/pod_run.sh (a capped unit, nice 19), each other tool in a capped
# scope (the launcher, examples/_shared/launch.sh): fetch Pfam's RNA_pol (PF00940) by name; find the control, T7
# gene 1; search hmm at 30 bits with 5 genes a side; nbr flanks with figures; the SQL checks (sql/checks.sql with
# the run, sql/checks_coord.sql without it); the hmmsearch oracle (examples/_shared/oracle.py) on the hits and the
# Pfam-annotated proteins; the values read from the run's files (observe.py); then the shared check.py against
# expected.tsv and recorded/$RELEASE.tsv. Writes nothing under ROOT. Exit: check.py's (1 when an invariant fails).
set -euo pipefail

: "${TREE:?set TREE}" "${PROFILE:?set PROFILE}" "${IMAGE:?set IMAGE}" "${DIGEST:?set DIGEST}"
: "${ROOT:?set ROOT}" "${RELEASE:?set RELEASE}" "${OUT:?set OUT}"
PFAM=${PFAM:-/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2/Pfam-A.hmm.gz}
THREADS=${THREADS:-2}
BIND=${BIND:-/tmp/if2453}
UNIT=${UNIT:-rnap}
MEM=${MEM:-16G}
EX="$TREE/examples/rnap"
SHARED="$TREE/examples/_shared"
PFAM_DIR=$(dirname "$PFAM")
. "$SHARED/launch.sh"

# The positive control: T7 gene 1, the RNA polymerase (GenBank NC_001604 3171..5822), 0-based half-open.
POS_T7=NC_001604:3170-5822
T7_FEATURE=39f9ef263e234b012d2671634468546e
BITS=30  # GUIDE.md section 4's cutoff, above the model's gathering cutoff (22 in Pfam 38.2)

[ ! -e "$ROOT/STOP" ] || die "a STOP file exists in $ROOT; not starting (AGENTS 7.4)"
case "$OUT/" in "$ROOT"/*) die "OUT is under ROOT" ;; esac
[ -f "$PFAM" ] || die "no Pfam file $PFAM"
[ -f "$EX/expected.tsv" ] || die "TREE has no $EX/expected.tsv"
[ -f "$SHARED/check.py" ] || die "TREE has no $SHARED/check.py"
if [ -e "$OUT" ] && [ -n "$(ls -A "$OUT")" ]; then
  die "OUT $OUT is not empty"
fi
mkdir -p "$OUT/logs" "$OUT/hmm" "$OUT/tmp"

run_id() {  # the run id a search's log names
  sed -n '/search run: [0-9a-f]\{32\}/{s/.*search run: \([0-9a-f]\{32\}\).*/\1/p;q;}' "$1"
}

note "release $RELEASE, root $ROOT, out $OUT, $THREADS threads"

# 1. The model, by name from Pfam (read-only; nothing is downloaded).
tool "$OUT/logs/hmmfetch.time" hmmfetch -o "$OUT/hmm/RNA_pol.hmm" "$PFAM" RNA_pol > "$OUT/logs/hmmfetch.log"
sha256sum "$OUT/hmm/RNA_pol.hmm" > "$OUT/hmm/SHA256SUMS"

# 2. The positive control must name one gene in this release.
bb find-t7 find --root "$ROOT" --release "$RELEASE" --accession "$POS_T7"

# 3. The search, as GUIDE.md section 4 runs it: 30 bits, 5 genes on each side.
S="--root $ROOT --release $RELEASE --profile $PROFILE --runs $OUT/runs --work $OUT/work --neighbors 5"
S="$S --cpus $THREADS --mem-gb 14 --threads $THREADS"
bb search search hmm $S --hmm "$OUT/hmm/RNA_pol.hmm" --bits "$BITS" --positive "$POS_T7" --out "$OUT/rnap"
RUN=$(run_id "$OUT/logs/search.log")
[ -n "$RUN" ] || die "no run id in the search log"
note "run $RUN"

# 4. Loci grouped by their flanks, as GUIDE.md section 6 does it (NR50 tokens, 5 genes a side, one group per
#    anchor NR90 cluster), with one figure per group.
F="--root $ROOT --profile $PROFILE --runs $OUT/runs --threads $THREADS"
bb flanks nbr flanks "$RUN" $F --window 5g --level nr50 --figure --out "$OUT/flanks"

# 5. The SQL checks, on the release's catalog with the run attached, and the control's neighbours again by
#    coordinates from the release alone (no run, no nbr).
Q="--root $ROOT --release $RELEASE --runs $OUT/runs --threads $THREADS"
bb sql-checks sql $Q "$EX/sql/checks.sql" --out "$OUT/observed_sql.tsv"
bb sql-coord sql --root "$ROOT" --release "$RELEASE" --threads "$THREADS" "$EX/sql/checks_coord.sql" \
  --out "$OUT/observed_coord.tsv"

# 6. The oracle: hmmsearch (HMMER, in the image) on every hit and every Pfam-annotated RNA_pol protein.
bb sql-uids sql $Q "$EX/sql/oracle_uids.sql" --out "$OUT/oracle_uids.txt"
bb seq-oracle seq --root "$ROOT" --release "$RELEASE" --proteins-file "$OUT/oracle_uids.txt" --out "$OUT/oracle.faa"
tool "$OUT/logs/hmmsearch.time" hmmsearch --cpu "$THREADS" -Z 1e8 --domZ 1e8 -T 0 --noali \
  --tblout "$OUT/oracle.tbl" -o "$OUT/logs/hmmsearch.out" "$OUT/hmm/RNA_pol.hmm" "$OUT/oracle.faa"
tool "$OUT/logs/oracle.time" python -I "$SHARED/oracle.py" --hmm "$OUT/hmm/RNA_pol.hmm" --tblout "$OUT/oracle.tbl" \
  --prefix rnap --model RNA_pol --bits "$BITS" "$OUT/rnap/hits.tsv" > "$OUT/observed_oracle.tsv"

# 7. The values read from the run's files.
tool "$OUT/logs/observe.time" python -I "$EX/observe.py" --hmm "$OUT/hmm/RNA_pol.hmm" --hits "$OUT/rnap/hits.tsv" \
  --neighbors "$OUT/rnap/neighbors.tsv" --flanks "$OUT/flanks" --control "$T7_FEATURE" > "$OUT/observed_files.tsv"

# 8. One observed.tsv (check_id, value), then the comparison. Its lines are the check's result; recorded_new.tsv
#    holds this run's release-specific values, to add as recorded/$RELEASE.tsv for a release that has none.
cp "$OUT/observed_sql.tsv" "$OUT/observed.tsv"
tail -n +2 "$OUT/observed_coord.tsv" >> "$OUT/observed.tsv"
tail -n +2 "$OUT/observed_oracle.tsv" >> "$OUT/observed.tsv"
tail -n +2 "$OUT/observed_files.tsv" >> "$OUT/observed.tsv"
set +e
tool "$OUT/logs/check.time" python -I "$SHARED/check.py" --expected "$EX/expected.tsv" \
  --recorded "$EX/recorded/$RELEASE.tsv" --write-recorded "$OUT/recorded_new.tsv" \
  "$OUT/observed.tsv" > "$OUT/check.txt"
rc=$?
set -e
cat "$OUT/check.txt"
note "done: check.py exited $rc; outputs in $OUT"
exit "$rc"
