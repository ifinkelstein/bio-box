#!/bin/sh
# Laptop-side: wait until any of the given files exists on the pod, checking every INTERVAL seconds.
#   sh waitfor.sh INTERVAL FILE [FILE ...]
INTERVAL=$1
shift
while true; do
  for F in "$@"; do
    if ssh -o BatchMode=yes -o ConnectTimeout=20 if2453@wilkcomp02.ccbb.utexas.edu "test -e $F"; then
      echo "present: $F"
      exit 0
    fi
  done
  sleep "$INTERVAL"
done
