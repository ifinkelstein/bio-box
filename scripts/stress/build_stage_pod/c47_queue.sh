#!/bin/bash
# Card C-47: the resumed build of stage C5 (DEC-56), queued behind stage V1 (card C-48c part 1) as one detached unit,
# so that it needs neither an agent nor the laptop. One build runs in the root at a time.
#
#   systemd-run --user --unit=c47-queue --collect -p MemoryMax=4G -p MemorySwapMax=0 \
#     -p StandardOutput=append:/tmp/if2453/c47/queue.out -p StandardError=append:/tmp/if2453/c47/queue.out \
#     nice -n 19 bash /tmp/if2453/c47/bin/c47_queue.sh
#
# 1. Every POLL seconds, wait until V1's `after` snapshot is complete (/tmp/if2453/c48c/after/DONE) and no `c48c-*`
#    service is active; also wait while a release-train unit (`c46-*`) is active (c47_driver.sh's own build_units
#    refuses then). If V1 failed instead (/tmp/if2453/c48c/run/FAILED or /tmp/if2453/c48c/after/FAILED), or the root
#    has a STOP file, exit non-zero with one clear line and start nothing: a stopped V1 is the orchestrator's to read.
# 2. Then the driver's modes rebefore, run, after, check, backup, in order, stopping at the first that fails.
#    c47_driver.sh run itself refuses while any build unit (`c48c-*` among them) is active: the second lock.
# The driver's own guards are unchanged: 280 CPU-minutes per process, free space, DEC-43, DEC-53; DEC-56 accepts
# C5's one ingest record. A guard stop stays a stop.
set -u

C=${C47_DIR:-/tmp/if2453/c47}  # the overrides exist for tests/r1/test_c47_driver.py; the defaults are the pod's
ROOT=${C47_ROOT:-/tmp/if2453/r10/root}
V=${C47_V1_DIR:-/tmp/if2453/c48c}
BIN=${C47_BIN:-$C/bin}
POLL=${POLL:-60}
MAX_WAIT_H=${MAX_WAIT_H:-48}
# pod_run's CPUQuota for the build unit, as the first launch c47-c5-1 had it; build-stage's own allowance holds the
# account's cap of the hour (40 by day, 56 at night and at weekends) inside it.
export Q=${Q:-56}

say() { echo "$(date '+%F %T %Z') c47_queue: $*"; }

active() {
  systemctl --user list-units --type=service --state=active,activating,deactivating,reloading --no-legend --plain \
    "$@" | awk '{print $1}' | paste -sd' ' -
}

say "start (pid $$); waiting for $V/after/DONE and no active c48c-* or c46-* service; poll ${POLL}s, at most ${MAX_WAIT_H}h"
t0=$(date +%s)
k=0
while :; do
  if [ -e "$ROOT/STOP" ]; then
    say "STOPPED: $ROOT/STOP exists; starting nothing"
    exit 3
  fi
  for f in "$V/run/FAILED" "$V/after/FAILED"; do
    if [ -e "$f" ]; then
      say "STOPPED: V1 failed ($f: $(head -c 400 "$f")); starting nothing: the orchestrator reads it"
      exit 3
    fi
  done
  v1=$(active 'c48c-*')
  tr=$(active 'c46-*')
  if [ -e "$V/after/DONE" ] && [ -z "$v1" ] && [ -z "$tr" ]; then
    say "V1's after snapshot is done and no c48c-* or c46-* service is active"
    break
  fi
  if [ $(( $(date +%s) - t0 )) -gt $(( MAX_WAIT_H * 3600 )) ]; then
    say "STOPPED: waited over ${MAX_WAIT_H}h; starting nothing"
    exit 5
  fi
  [ $((k % 30)) = 0 ] && say "waiting: after/DONE $([ -e "$V/after/DONE" ] && echo yes || echo no); c48c active: '${v1}'; c46 active: '${tr}'"
  k=$((k + 1))
  sleep "$POLL"
done

if [ -e "$C/run/FAILED" ] && [ ! -e "$C/run/DONE" ]; then
  # The first launch's stop record (c47-c5-1, 2026-11-19, the ingest record DEC-56 accepts). It never blocked `run`;
  # it is kept under the unit's name so that a FAILED in run/ means this launch.
  mv -n "$C/run/FAILED" "$C/run/FAILED.c47-c5-1"
  say "kept the first launch's stop record as $C/run/FAILED.c47-c5-1"
fi

# Each mode must exit 0 AND leave its marker: on 2026-11-20 `check` failed while the driver exited 0, and the queue,
# which then judged by the exit status alone, ran `backup` and reported "all passed".
marker() {
  case "$1" in
    rebefore) echo "$C/before/REBEFORE" ;;
    run) echo "$C/run/DONE" ;;
    after) echo "$C/after/DONE" ;;
    check) echo "$C/check/DONE" ;;
    backup) echo "$C/backup/DONE" ;;
  esac
}

for m in rebefore run after check backup; do
  say "driver mode $m: start"
  bash "$BIN/c47_driver.sh" "$m"
  rc=$?
  say "driver mode $m: exit $rc"
  if [ "$rc" != 0 ]; then
    say "STOPPED at mode $m (exit $rc); read $C/driver.log"
    exit "$rc"
  fi
  if [ ! -e "$(marker "$m")" ]; then
    say "STOPPED at mode $m: it exited 0 but left no $(marker "$m"); read $C/driver.log"
    exit 6
  fi
done
say "done: rebefore, run, after, check and backup all passed"
