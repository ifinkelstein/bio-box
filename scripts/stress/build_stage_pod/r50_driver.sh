#!/bin/bash
# The 50 Gbp build, parts B1 to B3 (planning/measurements/2026-11-11_build_r50.md), unattended. Host side, bash.
# Launch it in its own small unit on the pod host (copies of this script and r50_watch.py, r50_bb.sh in
# /tmp/if2453/r50/bin, so the mirrored tree the build runs from is not touched while it runs):
#
#   systemd-run --user --unit=r50-driver --collect -p MemoryMax=8G -p MemorySwapMax=0 \
#     -p StandardOutput=append:/tmp/if2453/r50/logs/r50-driver.out \
#     -p StandardError=append:/tmp/if2453/r50/logs/r50-driver.out \
#     nice -n 19 bash /tmp/if2453/r50/bin/r50_driver.sh
#
# Root /tmp/if2453/r10/root (stages S1, P1, P2 and their releases are there already). Release rM.1 is over
# S1,P1,P2,P3,P4,P5,P6 without a new library scan; rM.2 the same stages after Pfam 38.2 is scanned on P3 to P6.
# Every build is `build-stage` through scripts/pod_run.sh, units r50-<stage>-<n> (a new n per launch), and the
# driver waits for each unit to end (polled once a minute) and goes on when it ended as expected:
#
#   B1  P3, P4, P5: `build-stage --stage S --release rM.1 --stages <all seven>`, no --library. build-stage has
#       no stage-only mode: it runs the stage's steps (inputs, ingest, gene calling, IGRs, tapes) and then
#       refuses at rM.1's cluster step, exit 2, because the next stage has no gene-call record yet (the
#       documented gate of card B-36, `_cluster_gate` in src/biobox/stage.py; RUNBOOK §2). That refusal, with
#       exactly that message and the stage's tapes record present, is the expected end; anything else stops
#       the driver. --threads 16 (gene calling), --scan-jobs 4 --scan-cpus 4 (ingest batches, if any).
#   B2  P6 the same, but every stage is then gene-called, so it runs on into rM.1's release steps: clustering of
#       about 47 M proteins in parts (card B-41), tokens, postings, dim, manifest, catalog. --threads 32,
#       --scan-jobs 8 --scan-cpus 4 (32 threads: the day budget of this work), fixed for the whole launch.
#       Then `release verify` of rM.1 must pass.
#   B3  P3, P4, P5, P6: the Pfam scan, `build-stage --stage S --release rM.1 ... --library pfam 38.2 <file>
#       gathering` (rM.1's release steps all have run records, so only the scan runs, in shards), then
#       `build-stage --stage P6 --release rM.2 ... --library ...` for rM.2's release steps at 32 threads.
#       Then `release verify` of rM.2 must pass (rM.1 and rN.3 are verified too and only logged: card C-19).
#
# Scan shards follow the clock (as r10_phase2.sh did). build-stage fixes --scan-jobs when it starts; since card
# C-12 it re-reads the hour's allowance before each launch, but that allowance is the profile's (40 by day, 56 or
# 40 at night), not this work's (32 by day, 40 at night), so a launch at night with 10 shards would go on with 40
# threads into the day. The driver therefore runs DAY_JOBS shards (8, 32 threads) when the day window (the
# profile's day_start_hour to day_end_hour, Monday to Friday) holds now or LEAD_MIN minutes from now, else
# NIGHT_JOBS (10, 40 threads), and switches between them during a scan: it writes ROOT/STOP, waits for
# build-stage's "launching no more shards" line, removes STOP at once so the shards still running publish, waits
# for the unit to end (exit 75) and relaunches with the new count. The shard that ends first after STOP publishes
# nothing and is rerun. Only scan launches switch; a STOP during clustering would discard a step.
#
# Guard (r50_watch.py, once a minute while a unit runs): a run record of this driver's time or a step command still
# running with more than 300 CPU-minutes, scratch over 300 GB, or /tmp under 200 GB free: the driver writes
# ROOT/STOP with the reason (build-stage stops at its next step boundary), waits for the unit and stops, FAILED.
# The STOP is left for a person to read and remove.
#
# Before a launch that clusters or scans, the driver waits while a unit r50-rN3-* or r50-compat-* (part A) is
# active, so the two parts never pass 32 threads together; B1's stage launches (16 threads) do not wait.
#
# Files: LOG=/tmp/if2453/r50/driver.log (every launch with load, top and df, every switch, guard lines every 10
# minutes); STATE=/tmp/if2453/r50/driver: units (unit, stage, release, jobs), guard.jsonl (every guard line),
# scratch_peak, verify-<release>.log, and at the end DONE or FAILED (with the reason inside).
# Pause: touch STATE/HALT: the driver exits before its next launch (a running unit goes on to its end; to stop it
# at once as well, `touch ROOT/STOP` or `systemctl --user stop <unit>`, RUNBOOK §3). Resume after any stop: remove
# HALT (and STOP, FAILED), launch the driver again: it skips each phase whose result exists and build-stage skips
# every step with a run record.
set -u

W=/tmp/if2453/r50
ROOT=/tmp/if2453/r10/root
TREE=/tmp/if2453/ci/r50
PROFILE=$W/profile.yaml
IMAGE=/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif
DIGEST=sha256:c30be58557d756095346308c58cc3098af4bbe66404f67ab87be3f64f147c898
PYPATH=$TREE/src:/tmp/if2453/pilot_s1/pydeps
PFAMDIR=/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2
PFAM=$PFAMDIR/Pfam-A.hmm.gz
UNITLOGS=$W/tmp/logs
SCRATCH=$W/tmp/scratch
STATE=$W/driver
LOG=$W/driver.log
BIN=$(cd "$(dirname "$0")" && pwd)
ALL=S1,P1,P2,P3,P4,P5,P6
STAGE_THREADS=${STAGE_THREADS:-16}
BIG_THREADS=${BIG_THREADS:-32}
DAY_JOBS=${DAY_JOBS:-8}
NIGHT_JOBS=${NIGHT_JOBS:-10}
CPUS=4
LEAD_MIN=${LEAD_MIN:-90}
INPUT_WAIT_S=${INPUT_WAIT_S:-3600}
STEP=domains-pfam-38.2
mkdir -p "$STATE"
[ -e "$STATE/since" ] || touch "$STATE/since"
GUARD_STOPPED=0

log() { echo "$(date '+%F %T %Z') $*" >> "$LOG"; }
fail() { log "FAILED: $*"; echo "$*" > "$STATE/FAILED"; exit 1; }
prof() { sed -n "s/^$1:[[:space:]]*\([^#[:space:]]*\).*/\1/p" "$PROFILE"; }
DAY_START=$(prof day_start_hour)
DAY_END=$(prof day_end_hour)

is_day() {  # $1: seconds since the epoch; local time
  local dow h
  dow=$(date -d "@$1" +%u)
  h=$(date -d "@$1" +%-H)
  [ "$dow" -le 5 ] && [ "$h" -ge "$DAY_START" ] && [ "$h" -lt "$DAY_END" ]
}
want_jobs() {
  local t
  t=$(date +%s)
  if is_day "$t" || is_day $((t + LEAD_MIN * 60)); then echo "$DAY_JOBS"; else echo "$NIGHT_JOBS"; fi
}
pending_shards() {  # shards of stage $1 neither recorded nor running under $2 jobs; -1 when no plan yet
  local plan
  plan=$(find "$ROOT/l1" -path "*library=pfam-38.2/stage=$1/plan.json" 2> /dev/null | head -n 1)
  [ -n "$plan" ] || { echo -1; return; }
  [ -e "$ROOT/runs/stage=$1/$STEP.json" ] && { echo 0; return; }
  python3 -I -c 'import json, sys, glob
n = len(json.load(open(sys.argv[1]))["shards"])
done = len(glob.glob(sys.argv[2] + "/" + sys.argv[3] + ".shard-*.json")) - len(glob.glob(sys.argv[2] + "/" + sys.argv[3] + ".shard-*.publishing.json"))
print(max(0, n - done - int(sys.argv[4])))' "$plan" "$ROOT/runs/stage=$1" "$STEP" "$2"
}
snapshot() {
  log "uptime: $(uptime)"
  log "account CPU now (top, %CPU of if2453 processes): $(top -b -n 1 -u if2453 | awk 'NR > 7 {s += $9} END {print s}')"
  top -b -n 1 -o %CPU -w 150 | sed -n '1,5p;7,20p' >> "$LOG"
  log "df /tmp: $(df -h /tmp | tail -n 1)"
  log "our units: $(systemctl --user list-units --type=service --state=active --no-legend --plain 'r50-*' | awk '{print $1}' | paste -sd' ' -)"
}
others_active() {  # part A's units (not their watchers)
  systemctl --user list-units --type=service --state=active --no-legend --plain 'r50-rN3-*' 'r50-compat-*' \
    | awk '{print $1}' | grep -v -- '-watch\.service$' | paste -sd' ' -
}
halt_check() {
  if [ -e "$STATE/HALT" ]; then log "HALT present; exiting before $1"; exit 0; fi
  if [ -e "$ROOT/STOP" ]; then fail "$ROOT/STOP present before $1 (not written by this run of the driver): $(head -c 300 "$ROOT/STOP")"; fi
}
wait_quiet() {  # $1: what waits
  local o said=""
  while o=$(others_active); [ -n "$o" ]; do
    halt_check "$1"
    [ "$said" = "$o" ] || log "waiting for part A's units ($o) to end before $1 (thread budget)"
    said=$o
    sleep 300
  done
}
wait_inputs() {  # $1: stage
  local t0
  t0=$(date +%s)
  until [ -s "$W/inputs/$1/genomes.json" ] && [ -s "$W/inputs/$1/objects.json" ]; do
    [ $(($(date +%s) - t0)) -lt "$INPUT_WAIT_S" ] || fail "inputs of $1 not in $W/inputs/$1 after $INPUT_WAIT_S s"
    halt_check "inputs of $1"
    sleep 60
  done
}
guard() {  # run the watcher; on a reason, write STOP once
  local line rc
  line=$(python3 -I "$BIN/r50_watch.py" --root "$ROOT" --scratch "$SCRATCH" --since "$STATE/since" --state "$STATE")
  rc=$?
  echo "$line" >> "$STATE/guard.jsonl"
  if [ "$rc" = 3 ] && [ "$GUARD_STOPPED" = 0 ]; then
    echo "r50 driver guard: $line" > "$ROOT/STOP"
    GUARD_STOPPED=1
    log "guard: writing $ROOT/STOP: $line"
  elif [ "$rc" != 0 ] && [ "$rc" != 3 ]; then
    log "guard: r50_watch.py exited $rc: $line"
  fi
}
next_unit() {  # $1: tag
  local n
  n=$(($(wc -l < "$STATE/units" 2> /dev/null || echo 0) + 1))
  echo "r50-$1-$n"
}
launch() {  # $1 unit, $2 pod_run --threads, rest: build-stage arguments
  local unit=$1 q=$2
  shift 2
  snapshot
  log "launch $unit (pod_run --threads $q): build-stage $*"
  (cd "$W" && bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" \
    --threads "$q" --bind /tmp/if2453 --bind "$PFAMDIR:$PFAMDIR:ro" --env "PYTHONPATH=$PYPATH" \
    "$unit" build-stage --profile "$PROFILE" --root "$ROOT" "$@") >> "$LOG" 2>&1
}
exit_of() { sed -n "s/^pod_run: $1 exited \([0-9]*\)$/\1/p" "$UNITLOGS/$1.log" | tail -n 1; }
wait_unit() {  # $1 unit, $2 stage, $3 jobs (switch shard counts during a scan when $3 is not empty)
  local unit=$1 s=$2 j=${3:-} k=0 p
  SWITCHING=0
  sleep 5
  while systemctl --user is-active --quiet "$unit"; do
    for _ in 1 2 3 4 5 6 7 8 9 10 11 12; do  # the log, during a switch: every 5 s
      if [ "$SWITCHING" = 1 ] && grep -q "launching no more shards" "$UNITLOGS/$unit.log"; then
        [ "$GUARD_STOPPED" = 1 ] || rm -f "$ROOT/STOP"
        log "switch: build-stage launches no more shards; $ROOT/STOP removed; waiting for the running shards"
        SWITCHING=2
      fi
      sleep 5
    done
    guard
    k=$((k + 1))
    [ $((k % 10)) != 1 ] || log "guard: $(tail -n 1 "$STATE/guard.jsonl" | cut -c1-400); unit log: $(grep -E 'run-step|build-stage' "$UNITLOGS/$unit.log" | tail -n 1 | cut -c1-200)"
    if [ -n "$j" ] && [ "$SWITCHING" = 0 ] && [ "$GUARD_STOPPED" = 0 ] && [ "$(want_jobs)" != "$j" ]; then
      p=$(pending_shards "$s" "$j")
      if [ "$p" -gt 0 ]; then
        log "switch: the hour wants $(want_jobs) shards at once, $unit runs $j, $p shards wait; writing $ROOT/STOP"
        echo "r50 driver: switching from $j to $(want_jobs) scan shards at once" > "$ROOT/STOP"
        SWITCHING=1
      fi
    fi
  done
  if [ "$SWITCHING" = 1 ] && [ "$GUARD_STOPPED" = 0 ]; then rm -f "$ROOT/STOP"; log "switch: unit ended before the line; $ROOT/STOP removed"; fi
  RC=$(exit_of "$unit")
  log "$unit ended, exit '${RC}'; $(grep -E 'pod_run: .* usage' "$UNITLOGS/$unit.log" | tail -n 1)"
  [ "$GUARD_STOPPED" = 0 ] || fail "the guard stopped $unit: $(head -c 400 "$ROOT/STOP")"
}
verify() {  # $1 release; writes STATE/verify-<release>.log; returns verify's exit code
  RO=none bash "$BIN/r50_bb.sh" cli release verify --root "$ROOT" --release "$1" > "$STATE/verify-$1.log" 2>&1
  local rc=$?
  log "verify $1: exit $rc: $(grep -E '^release ' "$STATE/verify-$1.log" | tail -n 1)"
  return $rc
}

# B1: stage steps only (the refusal at the cluster gate ends each launch)
stage_only() {  # $1 stage
  local s=$1 unit
  if [ -e "$ROOT/runs/stage=$s/tapes.json" ]; then log "B1 $s: tapes record present; nothing to do"; return; fi
  wait_inputs "$s"
  halt_check "B1 $s"
  unit=$(next_unit "$s")
  echo "$unit $s rM.1 stage-only" >> "$STATE/units"
  launch "$unit" "$STAGE_THREADS" --stage "$s" --release rM.1 --stages "$ALL" \
    --genomes "$W/inputs/$s/genomes.json" --objects "$W/inputs/$s/objects.json" \
    --threads "$STAGE_THREADS" --scan-jobs 4 --scan-cpus "$CPUS"
  wait_unit "$unit" "$s"
  if grep -q "$s.ingest: .* it runs as pass-1 shards" "$UNITLOGS/$unit.log"; then
    log "B1 $s: ingest ran in batch processes (card B-30)"
  else
    log "B1 $s: ingest ran as one process (no batch line in the log)"
  fi
  if [ "$RC" = 2 ] && [ -e "$ROOT/runs/stage=$s/tapes.json" ] \
     && grep -Eq "build-stage: refused before rM\.1\.cluster: release rM\.1: stage P[3-6] has no gene-call run record" "$UNITLOGS/$unit.log"; then
    log "B1 $s: stage steps done; build-stage stopped at rM.1's cluster gate as expected: $(grep 'refused before rM.1.cluster' "$UNITLOGS/$unit.log" | tail -n 1 | cut -c1-200)"
    return
  fi
  [ "$RC" = 0 ] && [ -e "$ROOT/runs/stage=$s/tapes.json" ] && { log "B1 $s: exit 0"; return; }
  fail "B1 $s: $unit exit '$RC' is not the expected cluster-gate refusal (read $UNITLOGS/$unit.log)"
}

# B2 and the end of B3: a launch that runs a release's steps (rM.1 with P6's stage steps first; rM.2)
release_phase() {  # $1 release, $2 stage, $3 tag, $4 "lib" for --library
  local rel=$1 s=$2 tag=$3 unit lib=()
  [ "${4:-}" = lib ] && lib=(--library pfam 38.2 "$PFAM" gathering)
  if [ -e "$ROOT/runs/stage=$rel/catalog.json" ]; then log "$rel: catalog record present; nothing to do"; return; fi
  wait_inputs "$s"
  wait_quiet "$rel"
  halt_check "$rel"
  unit=$(next_unit "$tag")
  echo "$unit $s $rel release ${4:-}" >> "$STATE/units"
  launch "$unit" "$BIG_THREADS" --stage "$s" --release "$rel" --stages "$ALL" \
    --genomes "$W/inputs/$s/genomes.json" --objects "$W/inputs/$s/objects.json" "${lib[@]}" \
    --threads "$BIG_THREADS" --scan-jobs $((BIG_THREADS / CPUS)) --scan-cpus "$CPUS"
  wait_unit "$unit" "$s"
  [ "$RC" = 0 ] || fail "$rel: $unit exit '$RC' (read $UNITLOGS/$unit.log)"
  log "$rel: $(grep 'steps done' "$UNITLOGS/$unit.log" | tail -n 1 | cut -c1-200)"
}

# B3: the Pfam scan of one stage, the shard count following the clock
scan_phase() {  # $1 stage
  local s=$1 unit j q
  while :; do
    if [ -e "$ROOT/runs/stage=$s/$STEP.json" ]; then log "B3 $s: scan record present; nothing to do"; return; fi
    wait_quiet "B3 $s"
    halt_check "B3 $s"
    j=$(want_jobs)
    q=$((j * CPUS > BIG_THREADS ? j * CPUS : BIG_THREADS))
    unit=$(next_unit "$s")
    echo "$unit $s rM.1 scan $j" >> "$STATE/units"
    launch "$unit" "$q" --stage "$s" --release rM.1 --stages "$ALL" \
      --genomes "$W/inputs/$s/genomes.json" --objects "$W/inputs/$s/objects.json" \
      --library pfam 38.2 "$PFAM" gathering --threads "$BIG_THREADS" --scan-jobs "$j" --scan-cpus "$CPUS"
    wait_unit "$unit" "$s" "$j"
    if [ "$RC" = 0 ]; then log "B3 $s: $(grep "$s.$STEP " "$UNITLOGS/$unit.log" | tail -n 1 | cut -c1-200)"; return; fi
    if [ "$RC" = 75 ] && [ "$SWITCHING" != 0 ]; then continue; fi
    fail "B3 $s: $unit exit '$RC' is not 0 and not a switch (read $UNITLOGS/$unit.log)"
  done
}

log "driver start (pid $$), state $STATE, root $ROOT"
rm -f "$STATE/DONE"
for s in P3 P4 P5; do stage_only "$s"; done
release_phase rM.1 P6 P6
verify rM.1 || fail "rM.1 does not verify (read $STATE/verify-rM.1.log)"
for s in P3 P4 P5 P6; do scan_phase "$s"; done
release_phase rM.2 P6 rM2 lib
verify rM.2 || fail "rM.2 does not verify (read $STATE/verify-rM.2.log)"
verify rM.1 || log "rM.1 reports problems after the scans (card C-19 if they are only the new domain files)"
verify rN.3 || log "rN.3 reports problems after the scans (read $STATE/verify-rN.3.log)"
touch "$STATE/DONE"
log "DONE"
