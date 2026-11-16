#!/bin/bash
# Card C-36 (6): a rehearsal of the first N (300) genomes of stage STAGE, in the stage's order, as stage ${STAGE}R in a
# scratch view of the data root: the view links the stored stages S1, P1 to P6 and H1 (c26_view.py), every launch
# binds the data root read-only, and the sample stage's files land in the view's own directories. Only the stage's
# first three steps run, each a `bio-box run-step` launched by pod_run.sh in the runtime image (by digest), in the
# argv `build-stage` gives them (`stage.steps`, `with_seq_stages`): inputs, ingest, gene calling. Copies of this script
# and c26_view.py, c36_stages.py, c36_rh_check.py, make_inputs.py run from /tmp/if2453/c36/bin. Each mode in its own
# small unit, for example
#
#   systemd-run --user --unit=c36-rh-driver-run --collect -p MemoryMax=4G -p MemorySwapMax=0 \
#     -p StandardOutput=append:/tmp/if2453/c36/driver.out -p StandardError=append:/tmp/if2453/c36/driver.out \
#     --setenv=STAGE=C5 nice -n 19 bash /tmp/if2453/c36/bin/c36_rehearse.sh run
#
# Modes (each refuses to start while the data root's STOP exists):
#   prep    The view (c26_view.py ROOT RH/root S1,P1,...,H1), the sample (c36_stages.py sample C/stages/STAGE N), its
#           inputs (make_inputs.py prok, the stage's rule) and a scratch profile (T threads, the real budget).
#   run     The three run-steps one after another (units c36-rh-inputs, c36-rh-ingest, c36-rh-genecall); a step with
#           a run record is passed over; each must exit 0.
#   check   c36_rh_check.py: the rehearsal's files against the sample's own inputs and the stored stages' contig_index.
set -u

C=/tmp/if2453/c36
STAGE=${STAGE:?STAGE=C<k>}
N=${N:-300}
S=${STAGE}R
RH=$C/rh
ROOT=/tmp/if2453/r10/root
TREE=${TREE:-/tmp/if2453/ci/C-36}
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
DIGEST=sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a
STORED=S1,P1,P2,P3,P4,P5,P6,H1
T=${T:-8}
BIN=$C/bin
LOGS=$C/logs
LOG=$C/driver.log
mkdir -p "$LOGS"

log() { echo "$(date '+%F %T %Z') $*" >> "$LOG"; }
fail() { log "FAILED: $*"; exit 1; }

py() {  # python in the runtime image, the tree's src on the path, the data root bound read-only
  systemd-run --user --scope --quiet --unit="c36-rh-py-$1" -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% \
    nice -n 19 singularity exec --cleanenv --env PYTHONPATH="$TREE/src" --env OMP_NUM_THREADS=2 \
    --env PYTHONUNBUFFERED=1 --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" "$IMAGE" python "${@:2}"
}

exit_of() { sed -n "s/^pod_run: $1 exited \([0-9]*\)$/\1/p" "$LOGS/$1.log" | tail -n 1; }

prep() {
  [ -e "$RH/inputs/genomes.json" ] && { log "rh prep: inputs present; nothing to do"; return; }
  python3 -I "$BIN/c26_view.py" "$ROOT" "$RH/root" "$STORED" >> "$LOG" 2>&1 || fail "rh prep: view"
  mkdir -p "$RH/tmp" "$RH/durable" "$RH/bulk" "$RH/inputs"
  py sample "$BIN/c36_stages.py" sample "$C/stages/$STAGE" "$N" "$RH/sample" > "$RH/sample.json" 2>> "$LOG" \
    || fail "rh prep: sample"
  log "rh prep: sample: $(tr -d '\n' < "$RH/sample.json")"
  py inputs "$BIN/make_inputs.py" prok "$RH/sample/lock.json" "$RH/sample/selection.tsv" "$RH/inputs" \
    --source-release datasets-v2-2026-11-13 > "$RH/inputs/report.json" 2>> "$LOG" || fail "rh prep: make_inputs.py"
  log "rh prep: sample inputs: $(cat "$RH/inputs/report.json")"
  cat > "$RH/profile.yaml" << EOF
max_threads_day: $T
max_threads_night: $T
mem_scope: 32G
tmp: $RH/tmp
durable: $RH/durable
bulk: $RH/bulk
nice_monolithic: 19
nice_fanout: 15
day_start_hour: 8
day_end_hour: 20
cpu_min_per_process: 300
cluster_cpu_min_max: 20000
small_file_floor_bytes: 1000
require_mem_cap: true
EOF
  log "rh prep: done"
}

step() {  # $1 step name, $2 unit, then the stage subcommand's arguments: one run-step through pod_run.sh, waited on
  local name=$1 unit=$2 rc
  shift 2
  if [ -e "$RH/root/runs/stage=$S/$name.json" ]; then log "rh run: $S/$name recorded; passed over"; return; fi
  [ -e "$ROOT/STOP" ] && fail "$ROOT/STOP exists; not starting $S.$name"
  log "rh run: launch $unit: run-step $S.$name: stage $*"
  bash "$TREE/scripts/pod_run.sh" --profile "$RH/profile.yaml" --image "$IMAGE" --digest "$DIGEST" --mem 32G \
    --threads "$T" --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --env "PYTHONPATH=$TREE/src" --log "$LOGS/$unit.log" \
    "$unit" run-step --profile "$RH/profile.yaml" --root "$RH/root" --stage "$S" --step "$name" --threads "$T" -- \
    python -m biobox.cli stage "$@" >> "$LOG" 2>&1 || fail "pod_run $unit"
  sleep 2
  while systemctl --user is-active --quiet "$unit"; do sleep 5; done
  rc=$(exit_of "$unit")
  log "rh run: $unit exit '$rc': $(grep -E 'published|failed|refused|usage' "$LOGS/$unit.log" | tail -n 2 | tr '\n' ' ')"
  [ "$rc" = 0 ] || fail "$unit exit '$rc'"
}

run() {
  local seq
  [ -e "$RH/profile.yaml" ] || fail "rh run: no prep"
  mkdir -p "$RH/tmp/scratch"
  step inputs c36-rh-inputs inputs --root "$RH/root" --stage "$S" --genomes "$RH/inputs/genomes.json"
  rm -rf "$RH/tmp/scratch/$S.ingest"
  mkdir -p "$RH/tmp/scratch/$S.ingest"
  step ingest c36-rh-ingest ingest --root "$RH/root" --stage "$S" --objects "$RH/inputs/objects.json" \
    --tmp "$RH/tmp/scratch/$S.ingest"
  # the stage, then every stage of the view whose seq_locator exists (`stage.with_seq_stages`)
  seq=$S,$STORED
  rm -rf "$RH/tmp/scratch/$S.genecall"
  mkdir -p "$RH/tmp/scratch/$S.genecall"
  step genecall c36-rh-genecall genecall --root "$RH/root" --stage "$S" --tmp "$RH/tmp/scratch/$S.genecall" \
    --workers "$((T - 1))" --seq-stages "$seq"
  log "rh run: done"
}

check() {
  local rc
  py check "$BIN/c36_rh_check.py" --root "$RH/root" --data-root "$ROOT" --stage "$S" --inputs "$RH/inputs" \
    --sample "$RH/sample" --contigs "$C/stages/$STAGE/contigs.tsv" --stored "$STORED" --pl1 /tmp/if2453/s2/inputs/PL1 \
    --logs "$LOGS" --out "$RH/check" \
    > "$RH/check.out" 2>&1
  rc=$?
  log "rh check: exit $rc: $(tail -n 1 "$RH/check.out" | cut -c1-600)"
  return $rc
}

MODE=${1:-}
log "c36 rehearsal: mode '$MODE' (pid $$), stage $STAGE as $S, tree $TREE"
[ -e "$ROOT/STOP" ] && { log "$ROOT/STOP exists; not starting rehearsal mode '$MODE'"; exit 75; }
case "$MODE" in
  prep) prep ;;
  run) run ;;
  check) check ;;
  *) echo "usage: STAGE=C<k> c36_rehearse.sh prep|run|check" >&2; exit 2 ;;
esac
log "c36 rehearsal: mode '$MODE' done"
