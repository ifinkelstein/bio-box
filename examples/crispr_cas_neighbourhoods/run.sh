#!/bin/bash
# The CRISPR-Cas neighbourhood example as a release check (examples/crispr_cas_neighbourhoods/README.md).
# Run ON THE POD HOST (not in the image), from any directory:
#
#   TREE=... PROFILE=... IMAGE=... DIGEST=sha256:... ROOT=... RELEASE=rN.3 OUT=/tmp/if2453/<you>/crispr-rN.3 \
#     bash "$TREE/examples/crispr_cas_neighbourhoods/run.sh"
#
# Required: TREE (a checkout of the code, also the image's PYTHONPATH), PROFILE (host profile), IMAGE (the .sif,
# with matplotlib), DIGEST (sha256: of the image), ROOT (data root, bound read-only), RELEASE, OUT (a new or empty
# scratch directory outside ROOT). Optional: PFAM (the Pfam-A.hmm.gz to fetch the anchors from), THREADS (2),
# CAP (the account's thread cap; default 40 Mon-Fri 08:00-20:00, else 56), BIND (the directory bound read-write
# into the image, holding TREE, OUT, the pack directory and the image; default /tmp/if2453), UNIT (unit name
# prefix; default crispr), MEM (pod_run.sh --mem for the searches; default 16G).
#
# What it does, each bio-box step through scripts/pod_run.sh (a capped unit, nice 19), each other tool in a
# capped scope: fetch the anchor HMMs TniQ (PF06527) and Cas_Csy3 (Cas7f, PF09615) from Pfam by name; search
# both with their gathering cutoffs; nbr summarize at the paper's windows (20 kb for TniQ, 10 kb for Cas7f); nbr
# flanks with figures; the SQL checks (sql/checks.sql); the hmmsearch oracle (oracle.py); then the shared
# check.py against expected.tsv and recorded/$RELEASE.tsv. Writes nothing under ROOT. Exit: check.py's (1 when an
# invariant fails). The launcher (room, bb, tool) and check.py are shared by every example: examples/_shared/.
set -euo pipefail

: "${TREE:?set TREE}" "${PROFILE:?set PROFILE}" "${IMAGE:?set IMAGE}" "${DIGEST:?set DIGEST}"
: "${ROOT:?set ROOT}" "${RELEASE:?set RELEASE}" "${OUT:?set OUT}"
PFAM=${PFAM:-/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2/Pfam-A.hmm.gz}
THREADS=${THREADS:-2}
BIND=${BIND:-/tmp/if2453}
UNIT=${UNIT:-crispr}
MEM=${MEM:-16G}
EX="$TREE/examples/crispr_cas_neighbourhoods"
SHARED="$TREE/examples/_shared"
PFAM_DIR=$(dirname "$PFAM")
. "$SHARED/launch.sh"

# The positive controls: the TniQ and the Cas7f of one transposon-associated type I-F element of Photobacterium
# damselae subsp. damselae 04Ya311, chromosome 2 (NZ_AP026781.1), 0-based half-open (README section 3).
POS_TNIQ=NZ_AP026781.1:426340-427525
POS_CAS7F=NZ_AP026781.1:423377-424433

[ ! -e "$ROOT/STOP" ] || die "a STOP file exists in $ROOT; not starting (AGENTS 7.4)"
case "$OUT/" in "$ROOT"/*) die "OUT is under ROOT" ;; esac
[ -f "$PFAM" ] || die "no Pfam file $PFAM"
[ -f "$EX/expected.tsv" ] || die "TREE has no $EX/expected.tsv"
[ -f "$SHARED/check.py" ] || die "TREE has no $SHARED/check.py"
if [ -e "$OUT" ] && [ -n "$(ls -A "$OUT")" ]; then
  die "OUT $OUT is not empty"
fi
mkdir -p "$OUT/logs" "$OUT/hmm" "$OUT/tmp" "$OUT/summary" "$OUT/flanks"

ga() {  # the sequence gathering cutoff of an HMMER3 file of one model
  sed -n 's/^GA  *\([0-9.]*\) .*/\1/p' "$1"
}

run_id() {  # the run id a search's log names
  sed -n '/search run: [0-9a-f]\{32\}/{s/.*search run: \([0-9a-f]\{32\}\).*/\1/p;q;}' "$1"
}

note "release $RELEASE, root $ROOT, out $OUT, $THREADS threads"

# 1. The anchors, by name from Pfam (one pass over the file), then one file per model.
printf 'TniQ\nCas_Csy3\n' > "$OUT/hmm/anchors.keys"
tool "$OUT/logs/hmmfetch.time" hmmfetch -o "$OUT/hmm/anchors.hmm" -f "$PFAM" "$OUT/hmm/anchors.keys" > "$OUT/logs/hmmfetch.log"
tool "$OUT/logs/hmmfetch.time" hmmfetch -o "$OUT/hmm/TniQ.hmm" "$OUT/hmm/anchors.hmm" TniQ >> "$OUT/logs/hmmfetch.log"
tool "$OUT/logs/hmmfetch.time" hmmfetch -o "$OUT/hmm/Cas_Csy3.hmm" "$OUT/hmm/anchors.hmm" Cas_Csy3 >> "$OUT/logs/hmmfetch.log"
GA_TNIQ=$(ga "$OUT/hmm/TniQ.hmm")
GA_CAS7F=$(ga "$OUT/hmm/Cas_Csy3.hmm")
[ -n "$GA_TNIQ" ] && [ -n "$GA_CAS7F" ] || die "no GA line in the fetched models"
sha256sum "$OUT"/hmm/*.hmm > "$OUT/hmm/SHA256SUMS"
note "anchors: TniQ GA $GA_TNIQ, Cas_Csy3 GA $GA_CAS7F"

# 2. The positive controls must name one gene each in this release.
bb find-tniq find --root "$ROOT" --release "$RELEASE" --accession "$POS_TNIQ"
bb find-cas7f find --root "$ROOT" --release "$RELEASE" --accession "$POS_CAS7F"

# 3. The searches, at the models' gathering cutoffs, with 10 genes on each side in the run's neighbourhoods.
S="--root $ROOT --release $RELEASE --profile $PROFILE --runs $OUT/runs --work $OUT/work --neighbors 10"
S="$S --cpus $THREADS --mem-gb 14 --threads $THREADS"
bb search-tniq search hmm $S --hmm "$OUT/hmm/TniQ.hmm" --bits "$GA_TNIQ" --positive "$POS_TNIQ" --out "$OUT/tniq"
bb search-cas7f search hmm $S --hmm "$OUT/hmm/Cas_Csy3.hmm" --bits "$GA_CAS7F" --positive "$POS_CAS7F" \
  --out "$OUT/cas7f"
RUN_TNIQ=$(run_id "$OUT/logs/search-tniq.log")
RUN_CAS7F=$(run_id "$OUT/logs/search-cas7f.log")
[ -n "$RUN_TNIQ" ] && [ -n "$RUN_CAS7F" ] || die "no run id in the search logs"
note "runs: TniQ $RUN_TNIQ, Cas7f $RUN_CAS7F"

# 4. Neighbours at the paper's windows (Peters et al. 2017, Methods), every family and cluster.
N="--root $ROOT --profile $PROFILE --runs $OUT/runs --by pfam,nr50 --top 0 --format tsv"
bb sum-tniq nbr summarize "$RUN_TNIQ" $N --window 20kb --out "$OUT/summary/tniq_20kb.tsv"
bb sum-cas7f nbr summarize "$RUN_CAS7F" $N --window 10kb --out "$OUT/summary/cas7f_10kb.tsv"

# 5. Loci grouped by their flanks. (a) Every TniQ locus as the guide does it (NR50 tokens, one group per anchor
#    NR90 cluster), no figures. (b) Every Cas7f locus in one group, NR50 tokens, 6 genes a side, one figure.
#    (c) Only the TniQ loci with a Cas domain within 10 genes (their uids from sql/uids_tniq_cas.sql), Pfam
#    tokens, 10 genes a side, one figure.
F="--root $ROOT --profile $PROFILE --runs $OUT/runs --threads $THREADS"
Q="--root $ROOT --release $RELEASE --runs $OUT/runs --threads $THREADS"
TV=$(sed -n 's/^# tape_version\t"\(.*\)"$/\1/p' "$OUT/tniq/neighbors.tsv")
[ -n "$TV" ] || die "no tape_version in $OUT/tniq/neighbors.tsv"
bb flanks-tniq nbr flanks "$RUN_TNIQ" $F --window 10g --level nr50 --out "$OUT/flanks/tniq_nr50"
bb flanks-cas7f nbr flanks "$RUN_CAS7F" $F --window 6g --level nr50 --group-by none --figure \
  --out "$OUT/flanks/cas7f_all_nr50"
bb sql-tniq-cas sql $Q "$EX/sql/uids_tniq_cas.sql" --out "$OUT/uids_tniq_cas.txt"
bb flanks-tniq-cas nbr flanks --anchors-file "$OUT/uids_tniq_cas.txt" --release "$RELEASE" --tape-version "$TV" $F \
  --window 10g --level pfam --group-by none --figure --out "$OUT/flanks/tniq_cas_pfam"
bb sql-tniq-cas-loci sql $Q "$EX/sql/q03_tniq_cas_loci.sql" --out "$OUT/tniq_cas_loci.tsv"

# 6. The SQL checks, on the release's catalog with both runs attached, and again by base-pair windows from the
#    release's Pfam scan alone (no run, no nbr).
bb sql-checks sql $Q "$EX/sql/checks.sql" --out "$OUT/observed_sql.tsv"
bb sql-checks-bp sql --root "$ROOT" --release "$RELEASE" --threads "$THREADS" "$EX/sql/checks_bp.sql" \
  --out "$OUT/observed_bp.tsv"

# 7. The oracle: hmmsearch (HMMER, in the image) on every hit and every Pfam-annotated protein of the two models.
bb sql-uids sql $Q "$EX/sql/oracle_uids.sql" --out "$OUT/oracle_uids.txt"
bb seq-oracle seq --root "$ROOT" --release "$RELEASE" --proteins-file "$OUT/oracle_uids.txt" --out "$OUT/oracle.faa"
tool "$OUT/logs/hmmsearch.time" hmmsearch --cpu "$THREADS" -Z 1e8 --domZ 1e8 -T 0 --noali \
  --tblout "$OUT/oracle.tbl" -o "$OUT/logs/hmmsearch.out" "$OUT/hmm/anchors.hmm" "$OUT/oracle.faa"
tool "$OUT/logs/oracle.time" python -I "$EX/oracle.py" "$OUT/hmm/anchors.hmm" "$OUT/oracle.tbl" \
  tniq="$OUT/tniq/hits.tsv" cas7f="$OUT/cas7f/hits.tsv" > "$OUT/observed_oracle.tsv"

# 8. One observed.tsv (check_id, value), then the comparison. Its lines are the check's result; recorded_new.tsv
#    holds this run's release-specific values, to add as recorded/$RELEASE.tsv for a release that has none.
cp "$OUT/observed_sql.tsv" "$OUT/observed.tsv"
tail -n +2 "$OUT/observed_bp.tsv" >> "$OUT/observed.tsv"
tail -n +2 "$OUT/observed_oracle.tsv" >> "$OUT/observed.tsv"
set +e
tool "$OUT/logs/check.time" python -I "$SHARED/check.py" --expected "$EX/expected.tsv" \
  --recorded "$EX/recorded/$RELEASE.tsv" --write-recorded "$OUT/recorded_new.tsv" \
  "$OUT/observed.tsv" > "$OUT/check.txt"
rc=$?
set -e
cat "$OUT/check.txt"
note "done: check.py exited $rc; outputs in $OUT"
exit "$rc"
