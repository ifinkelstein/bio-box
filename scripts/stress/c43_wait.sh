#!/bin/sh
# Wait on the pod until a systemd user unit is no longer active, then print how it ended (card C-43).
#   sh c43_wait.sh UNIT OUT_DIR
# Exits on every terminal state of the unit (finished, failed, stopped, never started).
UNIT=$1
OUT=$2
while systemctl --user is-active -q "$UNIT"; do
  sleep 30
done
echo "unit $UNIT: $(systemctl --user show -p Result -p ExecMainStatus "$UNIT" | tr '\n' ' ')"
ls "$OUT"/DONE 2>/dev/null || echo "no DONE in $OUT"
wc -l < "$OUT"/fetch_log.jsonl
