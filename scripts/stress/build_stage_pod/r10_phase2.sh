#!/bin/bash
# Phase 2 of the 10 Gbp rebuild (planning/measurements/2026-11-10_rebuild_10gbp.md): the Pfam scan of
# stages S1, P1 and P2 through `build-stage`'s sharded path, then release rN.2 over S1,P1,P2. Host side,
# bash. Launch it in its own small unit on the pod host:
#
#   systemd-run --user --unit=r10-p2-driver --collect -p MemoryMax=32G -p MemorySwapMax=0 \
#     nice -n 19 bash /tmp/if2453/r10/bin/r10_phase2.sh
#
# (/tmp/if2453/r10/bin holds copies of this script and r10_queries.py, so the mirrored tree that the
# build runs from is not touched while it runs.)
#
# For each stage in order it launches `build-stage` with scripts/pod_run.sh (units r10-p2-<S>-<n>, a new
# n per launch), waits for the unit to end (polled once a minute) and goes on when it exits 0.
#
# The shard count follows the clock, because build-stage fixes --scan-jobs when it starts (card C-12 not
# merged): DAY_JOBS shards (16 threads) when the day window (the profile's day_start_hour to
# day_end_hour, Monday to Friday) holds now or LEAD_MIN minutes from now, else NIGHT_JOBS (32 threads).
# When the wanted count changes while shards are still waiting to be launched, the driver switches:
# it writes ROOT/STOP, waits (log polled every 5 s) for build-stage's "launching no more shards" line,
# removes STOP at once so the shards still running publish, waits for the unit to end (exit 75) and
# relaunches with the new count. The shard that ends first after STOP sees STOP and publishes nothing
# (run-step's rule); its next launch reruns it. Non-scan steps run at --threads 16, legal at any hour,
# and are never switched (a STOP during the 2.6 h clustering would discard it).
#
# Files under STATE=/tmp/if2453/r10/phase2: driver.log (every launch with load and top, every switch),
# units (unit, stage, release, jobs), DONE or FAILED at the end, verify.json and queries-rN.2.json.
# A file STATE/HALT makes the driver exit before its next launch (it does not stop a running unit).
# Rerunning the driver after it exited resumes: build-stage skips every step with a run record.
set -u

W=/tmp/if2453/r10
ROOT=$W/root
TREE=/tmp/if2453/ci/r10
PROFILE=$W/profile.yaml
IMAGE=/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif
DIGEST=sha256:c30be58557d756095346308c58cc3098af4bbe66404f67ab87be3f64f147c898
PYPATH=$TREE/src:/tmp/if2453/pilot_s1/pydeps
PFAMDIR=/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2
PFAM=$PFAMDIR/Pfam-A.hmm.gz
UNITLOGS=$W/tmp/logs
STATE=$W/phase2
DAY_JOBS=${DAY_JOBS:-4}
NIGHT_JOBS=${NIGHT_JOBS:-8}
CPUS=${CPUS:-4}
LEAD_MIN=${LEAD_MIN:-90}
STEP_THREADS=16
STEP=domains-pfam-38.2
mkdir -p "$STATE"

log() { echo "$(date '+%F %T %Z') $*" >> "$STATE/driver.log"; }
prof() { sed -n "s/^$1:[[:space:]]*\([^#[:space:]]*\).*/\1/p" "$PROFILE"; }
DAY_START=$(prof day_start_hour)
DAY_END=$(prof day_end_hour)

is_day() {  # $1: seconds since the epoch; local time
  local dow h
  dow=$(date -d "@$1" +%u)
  h=$(date -d "@$1" +%-H)
  [ "$dow" -le 5 ] && [ "$h" -ge "$DAY_START" ] && [ "$h" -lt "$DAY_END" ]
}
want_jobs() {  # $1 (tests): seconds since the epoch, default now
  local t=${1:-$(date +%s)}
  if is_day "$t" || is_day $((t + LEAD_MIN * 60)); then echo "$DAY_JOBS"; else echo "$NIGHT_JOBS"; fi
}
pending_shards() {  # shards of stage $1 neither recorded nor running under $2 jobs; -1 when no plan yet
  local plan
  plan=$(find "$ROOT/l1" -path "*library=pfam-38.2/stage=$1/plan.json" 2> /dev/null | head -n 1)
  [ -n "$plan" ] || { echo -1; return; }
  [ -e "$ROOT/runs/stage=$1/$STEP.json" ] && { echo 0; return; }
  python3 -c 'import json, sys, glob
n = len(json.load(open(sys.argv[1]))["shards"])
done = len(glob.glob(sys.argv[2] + "/" + sys.argv[3] + ".shard-*.json")) - len(glob.glob(sys.argv[2] + "/" + sys.argv[3] + ".shard-*.publishing.json"))
print(max(0, n - done - int(sys.argv[4])))' "$plan" "$ROOT/runs/stage=$1" "$STEP" "$2"
}
snapshot() {
  log "uptime: $(uptime)"
  log "account CPU now (top, %CPU of if2453 processes): $(top -b -n 1 -u if2453 | awk 'NR > 7 {s += $9} END {print s}')"
  top -b -n 1 -o %CPU -w 150 | sed -n '1,5p;7,20p' >> "$STATE/driver.log"
  log "df /tmp: $(df -h /tmp | tail -n 1)"
}

n=$(wc -l < "$STATE/units" 2> /dev/null || echo 0)
for spec in "S1 rN.1a S1" "P1 rN.1b S1,P1" "P2 rN.2 S1,P1,P2"; do
  set -- $spec
  S=$1 REL=$2 STAGES=$3
  while :; do
    if [ -e "$STATE/HALT" ]; then log "HALT present; exiting before launching $S"; exit 0; fi
    if [ -e "$ROOT/STOP" ]; then log "$ROOT/STOP present and not ours; exiting"; touch "$STATE/FAILED"; exit 1; fi
    J=$(want_jobs)
    Q=$((J * CPUS > STEP_THREADS ? J * CPUS : STEP_THREADS))
    n=$((n + 1))
    UNIT=r10-p2-$S-$n
    snapshot
    log "launch $UNIT: stage $S release $REL stages $STAGES, --scan-jobs $J --scan-cpus $CPUS, --threads $STEP_THREADS, pod_run --threads $Q"
    echo "$UNIT $S $REL $J" >> "$STATE/units"
    (cd "$W" && bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" \
      --threads "$Q" --bind /tmp/if2453 --bind "$PFAMDIR:$PFAMDIR:ro" --env "PYTHONPATH=$PYPATH" \
      "$UNIT" build-stage --profile "$PROFILE" --root "$ROOT" --stage "$S" --release "$REL" --stages "$STAGES" \
      --genomes "$W/inputs/$S/genomes.json" --objects "$W/inputs/$S/objects.json" \
      --library pfam 38.2 "$PFAM" gathering --threads "$STEP_THREADS" --scan-jobs "$J" --scan-cpus "$CPUS") \
      >> "$STATE/driver.log" 2>&1
    sleep 5
    switching=0
    while systemctl --user is-active --quiet "$UNIT"; do  # the unit: once a minute
      for _ in 1 2 3 4 5 6 7 8 9 10 11 12; do  # the log, during a switch: every 5 s
        if [ "$switching" = 1 ] && grep -q "launching no more shards" "$UNITLOGS/$UNIT.log"; then
          rm -f "$ROOT/STOP"
          log "switch: build-stage launches no more shards; $ROOT/STOP removed; waiting for the running shards"
          switching=2
        fi
        sleep 5
      done
      if [ "$switching" = 0 ] && [ "$(want_jobs)" != "$J" ]; then
        p=$(pending_shards "$S" "$J")
        if [ "$p" -gt 0 ]; then
          log "switch: the hour wants $(want_jobs) shards at once, $UNIT runs $J, $p shards wait; writing $ROOT/STOP"
          echo "r10 phase-2 driver: switching from $J to $(want_jobs) scan shards at once" > "$ROOT/STOP"
          switching=1
        fi
      fi
    done
    if [ "$switching" = 1 ]; then rm -f "$ROOT/STOP"; log "switch: unit ended before the line; $ROOT/STOP removed"; fi
    rc=$(sed -n "s/^pod_run: $UNIT exited \([0-9]*\)$/\1/p" "$UNITLOGS/$UNIT.log" | tail -n 1)
    log "$UNIT ended, exit '${rc}'"
    if [ "$rc" = 0 ]; then break; fi
    if [ "$rc" = 75 ] && [ "$switching" != 0 ]; then continue; fi
    log "$UNIT: exit '${rc}' is not 0 and not a switch; stopping (read $UNITLOGS/$UNIT.log)"
    touch "$STATE/FAILED"
    exit 1
  done
done

log "all three stages done; verifying rN.1 and rN.2, then the query set on rN.2"
singularity exec --cleanenv --bind /tmp/if2453 --env "PYTHONPATH=$PYPATH" "$IMAGE" python \
  "$(dirname "$0")/r10_verify.py" "$ROOT" rN.1a rN.1b rN.1 rN.2 > "$STATE/verify.json" 2>> "$STATE/driver.log"
log "verify exit $? (one JSON line per release in $STATE/verify.json)"
singularity exec --cleanenv --bind /tmp/if2453 --env "PYTHONPATH=$PYPATH" "$IMAGE" python \
  "$(dirname "$0")/r10_queries.py" "$ROOT" rN.2 "$STATE/queries-rN.2.json" > /dev/null 2>> "$STATE/driver.log"
log "queries exit $?"
touch "$STATE/DONE"
log "DONE"
