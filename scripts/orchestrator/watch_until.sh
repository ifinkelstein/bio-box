#!/bin/sh
# Poll the pod with one short ssh call per interval until a remote test succeeds, then print a remote report.
# A long-held ssh session drops when the laptop sleeps or the network blips (exit 255); short calls survive that.
# Usage: watch_until.sh '<remote test: exit 0 when done>' '<remote report>' [interval_s]
# Exit 0 when the test passed; 2 when the pod was unreachable for 6 intervals in a row.
HOST=if2453@wilkcomp02.ccbb.utexas.edu
TEST=$1
REPORT=$2
INTERVAL=${3:-300}
fails=0
while true; do
  ssh -o BatchMode=yes -o ConnectTimeout=20 "$HOST" "$TEST"
  rc=$?
  if [ "$rc" -eq 0 ]; then
    ssh -o BatchMode=yes -o ConnectTimeout=20 "$HOST" "$REPORT"
    exit 0
  fi
  if [ "$rc" -eq 255 ]; then
    fails=$((fails + 1))
    if [ "$fails" -ge 6 ]; then
      echo "POD UNREACHABLE for 6 intervals"
      exit 2
    fi
  else
    fails=0
  fi
  sleep "$INTERVAL"
done
