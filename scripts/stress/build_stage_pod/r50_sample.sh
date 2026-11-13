#!/bin/bash
# Watch one build unit of the 50 Gbp build (planning/measurements/2026-11-11_build_r50.md) once a minute with
# r50_watch.py, until the unit is no longer active. Host side; launch it in its own small unit:
#
#   systemd-run --user --unit=<unit>-watch --collect -p MemoryMax=2G -p MemorySwapMax=0 -p CPUQuota=50% \
#     nice -n 19 bash r50_sample.sh <unit> <root> <state dir>
#
# Appends one JSON line per minute to <state dir>/watch.jsonl. When the watcher names a reason to stop (a step
# over 300 CPU-minutes, scratch over 300 GB, /tmp under 200 GB free) it writes <root>/STOP with the reasons,
# once, and goes on watching: build-stage then stops at its next step boundary (RUNBOOK §3).
set -u
UNIT=$1
ROOT=$2
STATE=$3
BIN=$(dirname "$0")
mkdir -p "$STATE"
[ -e "$STATE/since" ] || touch "$STATE/since"
sleep 5
while systemctl --user is-active --quiet "$UNIT"; do
  line=$(python3 -I "$BIN/r50_watch.py" --root "$ROOT" --scratch /tmp/if2453/r50/tmp/scratch \
    --since "$STATE/since" --state "$STATE")
  rc=$?
  echo "$line" >> "$STATE/watch.jsonl"
  if [ "$rc" = 3 ] && [ ! -e "$ROOT/STOP" ]; then
    echo "r50_sample.sh: stopping $UNIT: $line" > "$ROOT/STOP"
    echo "$(date '+%F %T %Z') wrote $ROOT/STOP" >> "$STATE/watch.jsonl"
  fi
  sleep 60
done
echo "$(date '+%F %T %Z') $UNIT no longer active" >> "$STATE/watch.jsonl"
