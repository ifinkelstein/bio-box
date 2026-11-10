#!/bin/sh
# Card C-10: trigger one exercise mid-step. Waits (polling every 10 s, at most 6 h) until PATH exists,
# then sleeps DELAY seconds and runs the action, logging each event with its time:
#   sh watch_then.sh PATH DELAY stop-file ROOT        touch ROOT/STOP
#   sh watch_then.sh PATH DELAY stop-unit UNIT        systemctl --user stop UNIT
# Launch it in its own unit: systemd-run --user --unit=c10-watch-<x> --collect -p MemoryMax=1G sh watch_then.sh ...
set -eu
[ $# -eq 4 ] || { echo "usage: watch_then.sh PATH DELAY stop-file ROOT | stop-unit UNIT" >&2; exit 2; }
P=$1 DELAY=$2 ACTION=$3 ARG=$4
i=0
until [ -e "$P" ]; do
  i=$((i + 1))
  [ "$i" -le 2160 ] || { echo "$(date -Is) gave up waiting for $P"; exit 1; }
  sleep 10
done
echo "$(date -Is) $P exists; acting in $DELAY s"
sleep "$DELAY"
case "$ACTION" in
  stop-file) touch "$ARG/STOP"; echo "$(date -Is) created $ARG/STOP" ;;
  stop-unit) echo "$(date -Is) systemctl --user stop $ARG"; systemctl --user stop "$ARG"; echo "$(date -Is) stopped $ARG" ;;
  *) echo "unknown action $ACTION" >&2; exit 2 ;;
esac
