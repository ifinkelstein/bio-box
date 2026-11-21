#!/bin/bash
# Card A-12: build one stage unattended, its driver's modes in order, each on the marker the mode before it wrote.
#
#   systemd-run --user --unit=<stage>-queue --collect -p MemoryMax=4G -p MemorySwapMax=0 \
#     -p StandardOutput=append:/tmp/if2453/stages/<stage>/queue.out \
#     -p StandardError=append:/tmp/if2453/stages/<stage>/queue.out \
#     nice -n 19 bash "$TREE/scripts/stress/build_stage_pod/stage_queue.sh" STAGE
#
# 1. Every POLL seconds (60), until MAX_WAIT_H hours (48) have passed: start nothing while the data root has a STOP
#    file (exit 3); wait for each stage the conf names in AFTER until its driver's backup/DONE is there, and stop
#    (exit 3) when that stage's queue is not running and one of its modes left FAILED, REFUSED or STOPPED; wait while a
#    build unit of ours is active in the root (`stage_driver.sh STAGE info`, the driver's own lock).
# 2. Then `stage_driver.sh STAGE <mode>` for estimate, rehearse, before, run, after, check and backup, stopping at the
#    first mode that exits non-zero (with its status) or that exits 0 without writing <mode>/DONE (exit 6): on
#    2026-11-20 a mode failed while its driver exited 0, and a queue that judged by the exit status alone went on
#    (STATUS "6m"). A waiter (scripts/orchestrator/rel_wait.sh) reads the markers, never this script's last line.
# The driver keeps every guard (STOP, free space, 280 CPU-minutes per process, DEC-43 and DEC-53); a guard's stop
# writes run/FAILED, which only a person moves aside. Nothing here lifts a stop.
set -u

STAGE=${1:-}
[[ "$STAGE" =~ ^[A-Z][A-Za-z0-9]{0,15}$ ]] || { echo "usage: stage_queue.sh STAGE" >&2; exit 2; }
LC=${STAGE,,}
BIN=${STAGE_BIN:-$(cd "$(dirname "$0")" && pwd)}
STAGES_HOME=${STAGES_HOME:-/tmp/if2453/stages}
C=$STAGES_HOME/$LC
ROOT=/tmp/if2453/r10/root
POLL=${POLL:-60}
MAX_WAIT_H=${MAX_WAIT_H:-48}
MODES="estimate rehearse before run after check backup"

say() { echo "$(date '+%F %T %Z') stage_queue $STAGE: $*"; }

mkdir -p "$C" || exit 1
say "start (pid $$); driver $BIN/stage_driver.sh; poll ${POLL}s, at most ${MAX_WAIT_H}h of waiting"
t0=$(date +%s)
k=0
while :; do
  if [ -e "$ROOT/STOP" ]; then
    say "STOPPED: $ROOT/STOP exists; starting nothing"
    exit 3
  fi
  if ! bash "$BIN/stage_driver.sh" "$STAGE" info > "$C/queue.info" 2> "$C/queue.info.err"; then
    say "STOPPED: stage_driver.sh $STAGE info failed: $(tr '\n' ' ' < "$C/queue.info.err" | cut -c1-400)"
    exit 4
  fi
  after=$(sed -n 's/^AFTER=//p' "$C/queue.info")
  busy=$(sed -n 's/^BUSY=//p' "$C/queue.info")
  waiting=""
  for a in $after; do
    ad=$STAGES_HOME/${a,,}
    [ -e "$ad/backup/DONE" ] && continue
    if ! systemctl --user is-active --quiet "${a,,}-queue"; then
      bad=$(find "$ad" -mindepth 2 -maxdepth 2 \( -name FAILED -o -name REFUSED -o -name STOPPED \) 2> /dev/null \
        | LC_ALL=C sort | head -n 1)
      if [ -n "$bad" ]; then
        say "STOPPED: stage $a, which $STAGE comes after, did not pass ($bad: $(head -c 300 "$bad")); starting nothing"
        exit 3
      fi
    fi
    waiting="$waiting $a"
  done
  if [ -z "$waiting" ] && [ -z "$busy" ]; then
    say "ready: every stage it comes after is backed up (${after:-none}) and no build unit is active"
    break
  fi
  if [ $(( $(date +%s) - t0 )) -ge $(( MAX_WAIT_H * 3600 )) ]; then
    say "STOPPED: waited ${MAX_WAIT_H}h (for:${waiting:- nothing}; active: ${busy:-none}); starting nothing"
    exit 5
  fi
  [ $((k % 30)) = 0 ] && say "waiting for:${waiting:- nothing}; active build units: ${busy:-none}"
  k=$((k + 1))
  sleep "$POLL"
done

for m in $MODES; do
  say "mode $m: start"
  bash "$BIN/stage_driver.sh" "$STAGE" "$m"
  rc=$?
  say "mode $m: exit $rc"
  if [ "$rc" != 0 ]; then
    say "STOPPED at mode $m (exit $rc); read $C/driver.log and $C/$m/"
    exit "$rc"
  fi
  if [ ! -e "$C/$m/DONE" ]; then
    say "STOPPED at mode $m: it exited 0 but left no $C/$m/DONE; read $C/driver.log"
    exit 6
  fi
done
say "done: estimate, rehearse, before, run, after, check and backup all passed, each with its marker"
