#!/bin/sh
# Exit with one report when any ladder unit on the pod changes state, or when all logs go quiet.
# Usage: watch_ladder.sh [stall_minutes]
HOST=if2453@wilkcomp02.ccbb.utexas.edu
STALL=${1:-60}
UNITS="ladder-B ladder-Bp ladder-C"
snap() {
  ssh -o BatchMode=yes -o ConnectTimeout=20 "$HOST" "for u in $UNITS; do printf '%s=%s/%s ' \$u \$(systemctl --user show -p ActiveState --value \$u) \$(systemctl --user show -p Result --value \$u); done; echo; find /tmp/if2453/ladder/logs -type f -mmin -$STALL | wc -l"
}
prev=""
fails=0
while true; do
  out=$(snap 2>/dev/null) || out=""
  if [ -z "$out" ]; then
    fails=$((fails + 1))
    if [ "$fails" -ge 6 ]; then echo "POD UNREACHABLE for 30 min"; exit 2; fi
    sleep 300
    continue
  fi
  fails=0
  states=$(printf '%s\n' "$out" | sed -n 1p)
  fresh=$(printf '%s\n' "$out" | sed -n 2p | tr -d ' ')
  if [ -n "$prev" ] && [ "$states" != "$prev" ]; then
    echo "STATE CHANGE: was [$prev] now [$states]"
    ssh -o BatchMode=yes "$HOST" "ls -lt /tmp/if2453/ladder/logs | head -8; for f in \$(ls -t /tmp/if2453/ladder/logs/* | head -3); do echo \"--- \$f\"; tail -n 6 \$f; done"
    exit 0
  fi
  case "$states" in
    *=activ*)
      if [ "$fresh" = "0" ]; then
        echo "STALL: units [$states] but no log written for $STALL min"
        exit 3
      fi ;;
    *)
      if [ -z "$prev" ]; then echo "NO ACTIVE UNITS at start: [$states]"; exit 0; fi ;;
  esac
  prev=$states
  sleep 300
done
