#!/bin/bash
# Card C-36: the inputs of the six S2C phase-2 stages C1 to C6 and what they share, on the pod. Nothing is written in
# the data root (bound read-only into the image) or under /stor (the two archives are read in place and extracted
# under /tmp/if2453/s2/phase2). Copies of this script, make_inputs.py, s2_check_inputs.py and the c36_*.py scripts
# run from /tmp/if2453/c36/bin. Each mode in a small unit, for example
#
#   systemd-run --user --unit=c36-assign --collect -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=100% \
#     -p StandardOutput=append:/tmp/if2453/c36/driver.out -p StandardError=append:/tmp/if2453/c36/driver.out \
#     nice -n 19 ionice -c 2 -n 7 bash /tmp/if2453/c36/bin/c36_inputs.sh assign
#
# Modes (each refuses to start while ROOT/STOP exists, and logs to C/driver.log):
#   extract  The sha256 of each archive on /stor against sources.yaml, then `tar -xf` of both into X.tmp (host tar,
#            no owners kept), renamed to X when both are out. Refuses when X exists, or when /tmp would keep less than
#            FLOOR_GB (400) free after the archives' bytes. Logs `df` and `du` before and after.
#   assign   c36_stages.py assign: each genome's stage by the proposal's rule, the per-stage table against
#            stage_proposal.tsv, each stage's selection.tsv and lock.json under C/stages (report C/stages/assign.json).
#   verify   c36_verify.py for the six stages at once (one CPU each, scopes c36-verify-C<k>): every genome's .fna.gz
#            and FASTA against its lock, and C/stages/C<k>/contigs.tsv.
#   inputs   make_inputs.py prok for the six stages at once (scopes c36-inputs-C<k>), into IN/C<k>.tmp renamed to
#            IN/C<k>; refuses a stage whose IN/C<k> exists. Reports in C/inputs/C<k>.json.
#   check    s2_check_inputs.py for each stage (the inputs step's parser and checks; writes nothing) into
#            C/inputs/C<k>.check.json, and the sha256 of each stage's two files into C/inputs/SUMS.
#   analyze  c36_analyze.py into C/analyze (read only on ROOT and PL1's inputs); with REHEARSAL=<check json>, also the
#            protein and Pfam estimates; with PL1_STORED=1, the comparison with the stored contig_index of PL1.
set -u

C=/tmp/if2453/c36
ROOT=/tmp/if2453/r10/root
TREE=${TREE:-/tmp/if2453/ci/C-36}
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
STOR=/stor/scratch/Wilke/if2453/all-the-data/provisional/s2c_phase2
TAR1=$STOR/s2c_phase2_genomes.tar
TAR1_SHA=af3ae7c348328e78e7091e7790ba20e3103a1f0bef21eef1be18d36c2e36fcb9
TAR2=$STOR/s2c_phase2_rest_genomes.tar
TAR2_SHA=94cc747e86c619f80404941508ca5b8b3cffcf5d7302426f32830bd3db7c19df
X=/tmp/if2453/s2/phase2
IN=/tmp/if2453/s2/inputs
PL1=$IN/PL1
REMAINDER=/tmp/if2453/c27/out_ani/remainder.tsv
REMAINDER_SHA=48e6072846f6b40e396eea2639db61a2c02528be18e704feafaacc7462b93098
MEAS=$TREE/planning/measurements/2026-11-13_s2c_phase2_fetch
RELEASE_LABEL=datasets-v2-2026-11-13
FLOOR_GB=${FLOOR_GB:-400}
BIN=$C/bin
LOG=$C/driver.log
ST=$C/stages
mkdir -p "$C/inputs" "$C/analyze" "$ST"

log() { echo "$(date '+%F %T %Z') $*" >> "$LOG"; }
fail() { log "FAILED: $*"; exit 1; }
disk() { log "disk: /tmp $(df -Pk /tmp | awk 'NR==2 {print $4}') KiB free"; }

py() {  # $1 scope name, then python args: the runtime image, the tree's src on the path, the root bound read-only
  local unit=$1
  shift
  systemd-run --user --scope --quiet --unit="$unit" -p MemoryMax=12G -p MemorySwapMax=0 -p CPUQuota=100% \
    nice -n 19 ionice -c 2 -n 7 \
    singularity exec --cleanenv --env PYTHONPATH="$TREE/src" --env OMP_NUM_THREADS=1 --env PYTHONUNBUFFERED=1 \
    --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" "$IMAGE" python "$@"
}

extract() {
  local free need sha
  [ -e "$X" ] && fail "$X exists"
  disk
  free=$(df -Pk /tmp | awk 'NR==2 {print $4}')
  need=$(( ($(stat -c %s "$TAR1") + $(stat -c %s "$TAR2")) / 1024 + FLOOR_GB * 1000 * 1000 * 1000 / 1024 ))
  [ "$free" -gt "$need" ] || fail "/tmp has $free KiB free; the archives and the floor of $FLOOR_GB GB need $need"
  for t in "$TAR1:$TAR1_SHA" "$TAR2:$TAR2_SHA"; do
    sha=$(sha256sum "${t%%:*}" | cut -d' ' -f1)
    log "extract: sha256 of ${t%%:*}: $sha"
    [ "$sha" = "${t##*:}" ] || fail "sha256 of ${t%%:*} is not the one sources.yaml records"
  done
  rm -rf "$X.tmp"
  mkdir -p "$X.tmp"
  for t in "$TAR1" "$TAR2"; do
    [ -e "$ROOT/STOP" ] && fail "$ROOT/STOP appeared; stopping before $t"
    tar -C "$X.tmp" --no-same-owner --no-same-permissions -xf "$t" || fail "tar -xf $t"
    log "extract: $t extracted; $(find "$X.tmp" -type f | wc -l) files so far"
    disk
  done
  mv "$X.tmp" "$X"
  log "extract: $X: $(du -sb "$X" | cut -f1) bytes, $(find "$X" -name '*.fna.gz' | wc -l) .fna.gz files"
  disk
}

assign() {
  local sha rc
  sha=$(sha256sum "$REMAINDER" | cut -d' ' -f1)
  log "assign: sha256 of $REMAINDER: $sha"
  [ "$sha" = "$REMAINDER_SHA" ] || fail "sha256 of $REMAINDER is not the one the C-27 note records"
  py c36-assign-py "$BIN/c36_stages.py" assign "$REMAINDER" "$MEAS/not_served.txt" "$MEAS/stage_proposal.tsv" "$X" \
    "$ST" > "$ST/assign.json" 2> "$ST/assign.err"
  rc=$?
  log "assign: exit $rc: equal_to_proposal $(grep -o '"equal_to_proposal": [a-z]*' "$ST/assign.json") $(tail -n 2 "$ST/assign.err")"
  [ "$rc" = 0 ] || fail "c36_stages.py assign exit $rc"
}

each() {  # $1 mode: the six stages at once, one scope each; fails when any fails
  local mode=$1 k pids="" rc=0 p
  for k in 1 2 3 4 5 6; do
    "${mode}_one" "C$k" &
    pids="$pids $!"
  done
  for p in $pids; do
    wait "$p" || rc=1
  done
  [ "$rc" = 0 ] || fail "$mode: a stage failed (above)"
}

verify_one() {
  local s=$1 rc
  py "c36-verify-$s" "$BIN/c36_verify.py" "$ST/$s" "$ST/$s/contigs.tsv" > "$ST/$s/verify.json" 2> "$ST/$s/verify.err"
  rc=$?
  log "verify: $s exit $rc: $(tr -d '\n ' < "$ST/$s/verify.json" | cut -c1-400) $(tail -n 2 "$ST/$s/verify.err")"
  return $rc
}

inputs_one() {
  local s=$1 rc
  [ -e "$IN/$s" ] && { log "inputs: $IN/$s exists; refused"; return 1; }
  rm -rf "$IN/$s.tmp"
  py "c36-inputs-$s" "$BIN/make_inputs.py" prok "$ST/$s/lock.json" "$ST/$s/selection.tsv" "$IN/$s.tmp" \
    --source-release "$RELEASE_LABEL" > "$C/inputs/$s.json" 2> "$C/inputs/$s.err"
  rc=$?
  log "inputs: $s make_inputs.py prok exit $rc: $(cut -c1-400 "$C/inputs/$s.json") $(tail -n 2 "$C/inputs/$s.err")"
  [ "$rc" = 0 ] || return $rc
  mv "$IN/$s.tmp" "$IN/$s"
}

check() {
  local k s rc all=0
  : > "$C/inputs/SUMS"
  for k in 1 2 3 4 5 6; do
    s=C$k
    py "c36-check-$s" "$BIN/s2_check_inputs.py" "$s" "$IN/$s/genomes.json" "$IN/$s/objects.json" \
      > "$C/inputs/$s.check.json" 2> "$C/inputs/$s.check.err"
    rc=$?
    log "check: $s exit $rc: $(tr -d '\n' < "$C/inputs/$s.check.json" | sed 's/"objects": \[.*\]/"objects": [...]/' | cut -c1-500)"
    [ "$rc" = 0 ] || all=1
    (cd "$IN" && sha256sum "$s/genomes.json" "$s/objects.json") >> "$C/inputs/SUMS"
  done
  log "check: SUMS: $(tr '\n' ' ' < "$C/inputs/SUMS")"
  [ "$all" = 0 ] || fail "s2_check_inputs.py failed for a stage"
}

analyze() {
  local rc extra=()
  [ -n "${REHEARSAL:-}" ] && extra=(--rehearsal "$REHEARSAL")
  [ -n "${PL1_STORED:-}" ] && extra+=(--pl1-stored)
  py c36-analyze-py "$BIN/c36_analyze.py" --root "$ROOT" --pl1 "$PL1" --stages "$ST" --inputs "$IN" --out "$C/analyze" \
    "${extra[@]}" > "$C/analyze/out.json" 2> "$C/analyze/err.txt"
  rc=$?
  log "analyze: exit $rc: $(cut -c1-400 "$C/analyze/out.json") $(tail -n 2 "$C/analyze/err.txt")"
  [ "$rc" = 0 ] || fail "c36_analyze.py exit $rc"
}

MODE=${1:-}
log "c36 inputs: mode '$MODE' (pid $$), tree $TREE, image $IMAGE"
[ -e "$ROOT/STOP" ] && { log "$ROOT/STOP exists; not starting mode '$MODE'"; exit 75; }
case "$MODE" in
  extract) extract ;;
  assign) assign ;;
  verify) each verify ;;
  inputs) each inputs ;;
  check) check ;;
  analyze) analyze ;;
  *) echo "usage: c36_inputs.sh extract|assign|verify|inputs|check|analyze" >&2; exit 2 ;;
esac
log "c36 inputs: mode '$MODE' done"
