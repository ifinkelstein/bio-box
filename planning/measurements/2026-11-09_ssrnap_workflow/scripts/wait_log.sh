#!/bin/sh
# Server-side: wait up to $2 seconds for a line starting "exit" in log $1, then print the log.
n=0
until grep -q '^exit' "$1"; do
  sleep 10
  n=$((n + 10))
  if [ "$n" -ge "$2" ]; then
    echo "(still running after $2 s)"
    break
  fi
done
cat "$1"
