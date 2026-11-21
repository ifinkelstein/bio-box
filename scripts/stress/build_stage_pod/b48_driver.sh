#!/bin/bash
# Card B-48: gene calling in batch steps, on the pod; read-only on the data root (bound read-only into the image).
# Everything is written under B=/tmp/if2453/b48. Each mode runs in a small unit, for example
#
#   systemd-run --user --unit=b48-stats --collect -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% \
#     -p StandardOutput=append:/tmp/if2453/b48/driver.out -p StandardError=append:/tmp/if2453/b48/driver.out \
#     nice -n 19 ionice -c 2 -n 7 bash /tmp/if2453/b48/tree/scripts/stress/build_stage_pod/b48_driver.sh stats
#
# TREE is a copy of the branch's tree made for the session (`rsync` to /tmp/if2453/b48/tree); python runs in the
# runtime image with TREE/src alone on the path.
# Modes:
#   stats    b48_stage_inputs.py over every gene-called stage of ROOT, into B/stats/inputs.json.
#   pl1      b48_pl1_lengths.py over the PL1 FASTA and inputs (read only), into B/pl1/ (lengths.parquet, sample.json).
#   rate S   b48_call_rate.py: every EVERY-th (default 30) genome of stage S called in one thread, the calls timed
#            alone, into B/rate/S.json (the calls-only rate by mode and topology, as card C-26 measured PL1's).
#   pl1plan  b48_pl1_plan.py: PL1's estimate and batch plan on its real inputs at the real budget, into B/pl1/plan.json.
#   batches  A view of S1 (b48_view.py: ingested, not gene-called) in B/s1/view, and build-stage of S1 in it with a
#            scratch profile whose cpu_min_per_process (60) puts S1 over the margin, through pod_run.sh (unit
#            b48-s1-batches, T threads, J batches at once of CPUS workers); a watcher writes STOP in the view once the
#            genecall record or its intent file appears, so the launch ends before igr (exit 75).
#   one      S1 gene-called in one process in its own view B/s1/view-one (run-step, unit b48-s1-one, T workers).
#   compare A B NAME   b48_compare.py on S1's gene-call files of roots A and B, into B/s1/compare-NAME.json.
set -u

B=/tmp/if2453/b48
ROOT=/tmp/if2453/r10/root
TREE=${TREE:-$B/tree}
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
IN=/tmp/if2453/s2/inputs/PL1
FA=$IN/plsdb_2024_05_31_v2_sequences.fasta
BIN=$TREE/scripts/stress/build_stage_pod
LOG=$B/driver.log
mkdir -p "$B/stats" "$B/pl1"

log() { echo "$(date '+%F %T %Z') $*" >> "$LOG"; }
fail() { log "FAILED: $*"; exit 1; }

py() {  # python in the runtime image, the tree's src on the path, the root and the PL1 inputs bound read-only
  singularity exec --cleanenv --env PYTHONPATH="$TREE/src" --env OMP_NUM_THREADS=1 --env PYTHONUNBUFFERED=1 \
    --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$IN:$IN:ro" "$IMAGE" python "$@"
}

stats() {
  local rc
  py "$BIN/b48_stage_inputs.py" "$ROOT" "$B/stats/inputs.json" > "$B/stats/out.json" 2> "$B/stats/err.txt"
  rc=$?
  log "stats: exit $rc: $(cut -c1-600 "$B/stats/out.json") $(tail -n 3 "$B/stats/err.txt")"
  [ "$rc" = 0 ] || fail "b48_stage_inputs.py exit $rc"
}

pl1() {
  local rc
  py "$BIN/b48_pl1_lengths.py" "$FA" "$IN" "$B/pl1" > "$B/pl1/out.json" 2> "$B/pl1/err.txt"
  rc=$?
  log "pl1: exit $rc: $(cut -c1-600 "$B/pl1/out.json") $(tail -n 3 "$B/pl1/err.txt")"
  [ "$rc" = 0 ] || fail "b48_pl1_lengths.py exit $rc"
}

rate() {  # rate STAGE: b48_call_rate.py on every EVERY-th genome of STAGE, one thread, into B/rate/STAGE.json
  local rc s=$1
  mkdir -p "$B/rate"
  py "$BIN/b48_call_rate.py" "$ROOT" "$s" "$B/rate/$s.json" "${EVERY:-30}" > "$B/rate/$s.out" 2> "$B/rate/$s.err"
  rc=$?
  log "rate $s: exit $rc: $(cut -c1-800 "$B/rate/$s.out") $(tail -n 3 "$B/rate/$s.err")"
  [ "$rc" = 0 ] || fail "b48_call_rate.py $s exit $rc"
}

pl1plan() {  # the estimate and batch plan of PL1 on its real inputs (b48_pl1_plan.py), into B/pl1/plan.json
  local rc
  py "$BIN/b48_pl1_plan.py" "$ROOT" "$IN" "$B/pl1/lengths.parquet" "$B/pl1/plan.json" 300 > "$B/pl1/plan.out" \
    2> "$B/pl1/plan.err"
  rc=$?
  log "pl1plan: exit $rc: $(tr '\n' ' ' < "$B/pl1/plan.out" | cut -c1-2000) $(tail -n 3 "$B/pl1/plan.err")"
  [ "$rc" = 0 ] || fail "b48_pl1_plan.py exit $rc"
}

profile() {  # $1 file, $2 cpu_min_per_process, $3 tmp: a scratch profile of T threads, the pod's caps otherwise
  mkdir -p "$3" "$B/durable" "$B/bulk"
  cat > "$1" << EOF
max_threads_day: $T
max_threads_night: $T
mem_scope: 32G
tmp: $3
durable: $B/durable
bulk: $B/bulk
nice_monolithic: 19
nice_fanout: 15
day_start_hour: 8
day_end_hour: 20
cpu_min_per_process: $2
cluster_cpu_min_max: 20000
small_file_floor_bytes: 1000
require_mem_cap: true
EOF
}

exit_of() { sed -n "s/^pod_run: $1 exited \([0-9]*\)$/\1/p" "$LOGS/$1.log" | tail -n 1; }

batches() {  # S1 gene-called by build-stage in batches in the view B/s1/view, stopped by a STOP once its record is in
  local unit=b48-s1-batches v=$B/s1/view seen="" rc
  [ -e "$v/runs/stage=S1/genecall.json" ] && { log "batches: the view's genecall record exists; nothing to do"; return; }
  [ -e "$v" ] || python3 -I "$BIN/b48_view.py" "$ROOT" "$v" S1 >> "$LOG" 2>&1 || fail "batches: view"
  # cpu_min_per_process 60: S1's one-process estimate (about 144 CPU-minutes) is over 80% of it, so build-stage plans
  # batches of at most 48 (card B-48); every other key as the pod's profile
  profile "$B/s1/profile.yaml" 60 "$B/s1/ptmp"
  echo '[]' > "$B/none.json"
  log "batches: launch $unit: build-stage --root $v --stage S1 --threads $T --scan-jobs $J --scan-cpus $CPUS"
  bash "$TREE/scripts/pod_run.sh" --profile "$B/s1/profile.yaml" --image "$IMAGE" --digest "$DIGEST" --mem 32G \
    --threads "$T" --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --env "PYTHONPATH=$TREE/src" --log "$LOGS/$unit.log" \
    "$unit" build-stage --profile "$B/s1/profile.yaml" --root "$v" --stage S1 --release rB48.s1 --stages S1 \
    --genomes "$B/none.json" --objects "$B/none.json" --threads "$T" --scan-jobs "$J" --scan-cpus "$CPUS" \
    >> "$LOG" 2>&1 || fail "batches: pod_run"
  sleep 2
  while systemctl --user is-active --quiet "$unit"; do
    if [ -z "$seen" ] && { [ -e "$v/runs/stage=S1/genecall.json" ] || [ -e "$v/runs/stage=S1/genecall.publishing.json" ]; }; then
      echo "b48 watcher $(date '+%F %T.%N')" > "$v/STOP"
      seen=1
      log "batches: the genecall record or intent file appeared; wrote $v/STOP"
    fi
    sleep 0.5
  done
  rc=$(exit_of "$unit")
  log "batches: $unit exit '$rc'; $(grep -E 'STOP|refused|genecall' "$LOGS/$unit.log" | tail -n 4 | tr '\n' ' ')"
  log "batches: usage $(grep -E 'pod_run: .* usage' "$LOGS/$unit.log" | tail -n 1)"
  [ -e "$v/runs/stage=S1/genecall.json" ] || fail "batches: no genecall record in the view"
  grep -q "^b48 watcher" "$v/STOP" && rm "$v/STOP"
}

one() {  # S1 gene-called in one process (run-step, the current image) in its own view B/s1/view-one
  local unit=b48-s1-one v=$B/s1/view-one
  [ -e "$v/runs/stage=S1/genecall.json" ] && { log "one: the view's genecall record exists; nothing to do"; return; }
  [ -e "$v" ] || python3 -I "$BIN/b48_view.py" "$ROOT" "$v" S1 >> "$LOG" 2>&1 || fail "one: view"
  profile "$B/s1/profile-one.yaml" 300 "$B/s1/ptmp-one"
  mkdir -p "$B/s1/ptmp-one/scratch/S1.genecall"
  log "one: launch $unit: run-step genecall, one process, $T workers"
  bash "$TREE/scripts/pod_run.sh" --profile "$B/s1/profile-one.yaml" --image "$IMAGE" --digest "$DIGEST" --mem 32G \
    --threads "$T" --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --env "PYTHONPATH=$TREE/src" --log "$LOGS/$unit.log" \
    "$unit" run-step --profile "$B/s1/profile-one.yaml" --root "$v" --stage S1 --step genecall --threads "$T" -- \
    python -m biobox.cli stage genecall --root "$v" --gene-call-version pyrodigal-3.7.1-gv-0.3.2-p2 --stage S1 \
    --tmp "$B/s1/ptmp-one/scratch/S1.genecall" --workers "$T" --seq-stages S1 >> "$LOG" 2>&1 || fail "one: pod_run"
  sleep 5
  while systemctl --user is-active --quiet "$unit"; do sleep 5; done
  log "one: $unit exit '$(exit_of "$unit")'; usage $(grep -E 'pod_run: .* usage' "$LOGS/$unit.log" | tail -n 1)"
  [ "$(exit_of "$unit")" = 0 ] || fail "one: exit '$(exit_of "$unit")'"
}

compare() {  # compare A B NAME: b48_compare.py over S1's gene-call files of roots A and B, into B/s1/compare-NAME.json
  local rc
  py "$BIN/b48_compare.py" "$1" "$2" S1 "$B/s1/compare-$3.json" "$1/runs/stage=S1/genecall.json" \
    "$2/runs/stage=S1/genecall.json" > "$B/s1/compare-$3.out" 2> "$B/s1/compare-$3.err"
  rc=$?
  log "compare $3: exit $rc: $(cut -c1-1500 "$B/s1/compare-$3.out") $(tail -n 2 "$B/s1/compare-$3.err")"
  [ "$rc" = 0 ] || fail "b48_compare.py exit $rc"
}

MODE=${1:-}
T=${T:-8}
J=${J:-2}
CPUS=${CPUS:-4}
DIGEST=sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a
LOGS=$B/logs
mkdir -p "$LOGS" "$B/s1"
log "b48: mode '$MODE' (pid $$), tree $TREE, image $IMAGE"
case "$MODE" in
  stats) stats ;;
  pl1) pl1 ;;
  rate) rate "$2" ;;
  pl1plan) pl1plan ;;
  batches) batches ;;
  one) one ;;
  compare) compare "$2" "$3" "$4" ;;
  *) echo "usage: b48_driver.sh stats|pl1|rate STAGE|pl1plan|batches|one|compare A B NAME" >&2; exit 2 ;;
esac
# The mode's own status is the driver's exit status (card A-11; on 2026-11-20 c47_driver.sh `check` returned 1 and the
# driver exited 0, the status of the log line that followed).
rc=$?
if [ "$rc" = 0 ]; then
  log "b48: mode '$MODE' done"
else
  log "b48: mode '$MODE' FAILED with status $rc"
fi
exit "$rc"
