#!/bin/bash
# Card C-52 (DEC-63): stages PL2 and PL3 from card C-40's staging, on the pod. Nothing is written in the data root, in
# card C-40's staging or in its fetched files (all three bound read-only into the image), or under /stor. The tree is a
# copy synced to $W/tree. Each mode in a small detached unit at 1 CPU, e.g.
#
#   systemd-run --user --unit=c52-survey --collect -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=100% \
#     nice -n 19 ionice -c 2 -n 7 bash /tmp/if2453/c52/tree/scripts/stress/build_stage_pod/c52_run.sh survey
#
# Modes:
#   survey    c52_pl23.py survey over card C-40's staging (and PL1's circular genomes) into $W/survey.
#   deferred  c52_pl23.py deferred into $W/deferred/gut_plasmids_yu2024_circular_non_representatives.tsv, and its
#             sha256 and rows into $W/deferred/SUMS.
#   stage     c52_pl23.py stage into $OUT2 and $OUT3 (objects.json naming $FINAL2 and $FINAL3, default the same),
#             with $STAGE_ARGS (--ambiguous ..., --fragments ...) as given; then inputs_SUMS in each. C-40's staging is
#             reachable through the /tmp/if2453 bind here (a hard link cannot cross into a read-only bind); the script
#             only reads it, and `verify_old` checks its files against card C-40's SUMS after the run.
#   check     s2_check_inputs.py on $OUT2 and $OUT3 (the inputs step's parser and checks; writes nothing).
#   overlap   c52_pl23.py overlap of $OUT2 and $OUT3 against the root and the staged inputs of V2, GV1, G1, SPX1.
#   verify_old  sha256 of card C-40's staging against its SUMS (planning/measurements/2026-11-16_pl2_inputs/SUMS).
# A STOP file in the data root stops every mode before it starts.
set -u

W=/tmp/if2453/c52
ROOT=/tmp/if2453/r10/root
OLD=/tmp/if2453/s2/inputs/PL2
PL1=/tmp/if2453/s2/inputs/PL1
PLD=/tmp/if2453/s3/plasmids
TREE=${TREE:-$W/tree}
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
BSP=$TREE/scripts/stress/build_stage_pod
LOG=$W/logs/driver.log
mkdir -p "$W/logs"

log() { echo "$(date '+%F %T %Z') $*" >> "$LOG"; }
fail() { log "FAILED: $*"; exit 1; }

py() {  # python in the runtime image, the tree's src on the path; the root, C-40's staging and its files read-only
  singularity exec --cleanenv --env PYTHONPATH="$TREE/src" --env OMP_NUM_THREADS=1 --env PYTHONUNBUFFERED=1 \
    --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$OLD:$OLD:ro" --bind "$PL1:$PL1:ro" --bind "$PLD:$PLD:ro" \
    "$IMAGE" python "$@"
}

survey() {
  local rc
  py "$BSP/c52_pl23.py" survey --inputs "$OLD" --out "$W/survey" --pl1 "$PL1" > "$W/logs/survey.out" 2> "$W/logs/survey.err"
  rc=$?
  log "survey: exit $rc: $(cut -c1-600 "$W/logs/survey.out") $(tail -n 2 "$W/logs/survey.err")"
  [ "$rc" = 0 ] || fail "survey exit $rc"
}

deferred() {
  local rc f=gut_plasmids_yu2024_circular_non_representatives.tsv
  mkdir -p "$W/deferred"
  py "$BSP/c52_pl23.py" deferred --gut-table "$PLD/fasta/gut/predicted_plasmids_summary.txt" \
    --skipped "$OLD/skipped_records.tsv" --out "$W/deferred/$f" > "$W/logs/deferred.out" 2> "$W/logs/deferred.err"
  rc=$?
  log "deferred: exit $rc: $(cat "$W/logs/deferred.out") $(tail -n 2 "$W/logs/deferred.err")"
  [ "$rc" = 0 ] || fail "deferred exit $rc"
  (cd "$W/deferred" && sha256sum "$f" && wc -lc "$f") > "$W/deferred/SUMS"
}

OUT2=${OUT2:-$W/dec65/PL2}
OUT3=${OUT3:-$W/dec65/PL3}
FINAL2=${FINAL2:-/tmp/if2453/s2/inputs/PL2}
FINAL3=${FINAL3:-/tmp/if2453/s2/inputs/PL3}
STAGE_ARGS=${STAGE_ARGS:---rule dec65}
ROOT_STAGES=${ROOT_STAGES:-S1,P1,P2,P3,P4,P5,P6,H1,PL1,C5,V1}

py_rw() {  # as py, but C-40's staging only through the /tmp/if2453 bind, so that its objects can be hard-linked
  singularity exec --cleanenv --env PYTHONPATH="$TREE/src" --env OMP_NUM_THREADS=2 --env PYTHONUNBUFFERED=1 \
    --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$PLD:$PLD:ro" --bind /stor:/stor:ro "$IMAGE" python "$@"
}

stage() {
  local rc d
  # shellcheck disable=SC2086  # STAGE_ARGS is a list of options
  py_rw "$BSP/c52_pl23.py" stage --inputs "$OLD" --survey "$W/survey" --root "$ROOT" --root-stages "$ROOT_STAGES" \
    --out-pl2 "$OUT2" --out-pl3 "$OUT3" --final-pl2 "$FINAL2" --final-pl3 "$FINAL3" --work "$W/work/stage" $STAGE_ARGS \
    > "$W/logs/stage.out" 2> "$W/logs/stage.err"
  rc=$?
  log "stage: exit $rc: $(cut -c1-600 "$W/logs/stage.out") $(tail -n 2 "$W/logs/stage.err")"
  [ "$rc" = 0 ] || fail "stage exit $rc"
  for d in "$OUT2" "$OUT3"; do
    (cd "$d" && sha256sum genomes.json objects.json records.tsv trims.tsv skipped_records.tsv stage_report.json \
      ./*.fasta) > "$d/.inputs_SUMS.tmp" && mv "$d/.inputs_SUMS.tmp" "$d/inputs_SUMS" || fail "sums of $d"
    log "stage: $d/inputs_SUMS written"
  done
  verify_old
}

check() {
  local rc d s
  for s in PL2 PL3; do
    if [ "$s" = PL2 ]; then d=$OUT2; else d=$OUT3; fi
    py "$TREE/scripts/stress/build_stage_pod/s2_check_inputs.py" "$s" "$d/genomes.json" "$d/objects.json" \
      > "$W/logs/check_$s.json" 2> "$W/logs/check_$s.err"
    rc=$?
    log "check $s: exit $rc: $(tr -d '\n' < "$W/logs/check_$s.json" | cut -c1-600)"
    [ "$rc" = 0 ] || fail "s2_check_inputs.py $s exit $rc"
  done
}

overlap() {
  local rc
  py "$BSP/c52_pl23.py" overlap --staging "PL2=$OUT2" --staging "PL3=$OUT3" --root "$ROOT" \
    --root-stages "$ROOT_STAGES" \
    --other G1=/tmp/if2453/s2/inputs/G1 --other SPX1=/tmp/if2453/s2/inputs/SPX1 \
    --other V2=/tmp/if2453/s2/inputs/V2 --other GV1=/tmp/if2453/s2/inputs/GV1 \
    --work "$W/work/overlap" --out "$W/overlap.json" \
    > "$W/logs/overlap.out" 2> "$W/logs/overlap.err"
  rc=$?
  log "overlap: exit $rc: $(tail -n 2 "$W/logs/overlap.err")"
  [ "$rc" = 0 ] || fail "overlap exit $rc"
}

verify_old() {  # C-40's SUMS: genomes.json objects.json records.tsv skipped_records.tsv and the ten FASTAs
  (cd "$OLD" && sha256sum -c "$TREE/planning/measurements/2026-11-16_pl2_inputs/SUMS") > "$W/logs/verify_old.txt" 2>&1 \
    || fail "card C-40's staging differs from its SUMS: $(grep -v ': OK$' "$W/logs/verify_old.txt" | head -c 400)"
  log "verify_old: $(grep -c ': OK$' "$W/logs/verify_old.txt") files of $OLD equal card C-40's SUMS"
}

MODE=${1:-}
[ -e "$ROOT/STOP" ] && fail "$ROOT/STOP exists; mode '$MODE' not started"
log "c52: mode '$MODE' (pid $$), tree $TREE, image $IMAGE"
case "$MODE" in
  survey) survey ;;
  deferred) deferred ;;
  stage) stage ;;
  check) check ;;
  overlap) overlap ;;
  verify_old) verify_old ;;
  *) echo "usage: c52_run.sh survey|deferred|stage|check|overlap|verify_old" >&2; exit 2 ;;
esac
log "c52: mode '$MODE' done"
