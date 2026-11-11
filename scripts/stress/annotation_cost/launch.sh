#!/bin/bash
# Launch one card 0-09a job on the pod, detached, inside the image (run ON THE POD HOST).
#
#   bash launch.sh <name> <command...>
#
# Unit m09a-<name>, MemoryMax=32G, MemorySwapMax=0, CPUQuota=<threads*100>%, nice 19 (AGENTS §7.2 and the
# card's limits). The command runs inside the image under timed.py, which writes runs/m09a-<name>.json
# (wall, CPU-minutes, peak RSS); stdout and stderr go to logs/m09a-<name>.log, whose last line is
# "launch: m09a-<name> exited <code>".
#
# Before starting, the account's CPU use (sum of `ps -u if2453 -o pcpu=`, the ps process excluded) plus this job's threads*100 must
# be at most 5600% between 20:00 and 05:00 CST, or 4000% otherwise; else it polls every 4 minutes and
# gives up (exit 3) after 2 hours. M09A_THREADS sets the thread count (default 4, the card's maximum).
set -euo pipefail
M=/tmp/if2453/m09a
IMAGE=/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif
PFAM=/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2
THREADS=${M09A_THREADS:-4}
[ "$THREADS" -le 4 ] || { echo "launch: at most 4 threads" >&2; exit 2; }
[ $# -ge 2 ] || { echo "usage: launch.sh <name> <command...>" >&2; exit 2; }
UNIT=m09a-$1
shift
mkdir -p "$M/logs" "$M/runs"
deadline=$(( $(date +%s) + 7200 ))
while :; do
  # ps's own line is left out: pcpu is CPU time over lifetime, and the ps process reads several hundred %
  used=$(ps -u if2453 -o pcpu=,comm= | awk '$2 != "ps" {s+=$1} END {printf "%d", s}')
  h=$(TZ=America/Chicago date +%H)
  h=$((10#$h))
  if [ "$h" -ge 20 ] || [ "$h" -lt 5 ]; then cap=5600; else cap=4000; fi
  if [ $((used + THREADS * 100)) -le "$cap" ]; then break; fi
  echo "launch: $UNIT waiting at $(date -Is): account ${used}% + ${THREADS}00% > ${cap}%" >> "$M/logs/launch.log"
  if [ "$(date +%s)" -ge "$deadline" ]; then
    echo "launch: $UNIT gave up after 2 hours" | tee -a "$M/logs/launch.log" >&2
    exit 3
  fi
  sleep 240
done
echo "launch: $UNIT started at $(date -Is): account ${used}% + ${THREADS}00% <= ${cap}%" | tee -a "$M/logs/launch.log"
systemd-run --user --unit="$UNIT" --collect --quiet --same-dir \
  -p MemoryMax=32G -p MemorySwapMax=0 -p CPUQuota="$((THREADS * 100))%" \
  -p StandardOutput="append:$M/logs/$UNIT.log" -p StandardError="append:$M/logs/$UNIT.log" \
  /bin/sh -c 'unit=$1; shift; "$@"; rc=$?; echo "launch: $unit exited $rc"; exit $rc' sh "$UNIT" \
  nice -n 19 singularity exec --cleanenv --bind /tmp/if2453 --bind "$PFAM:$PFAM:ro" \
  --env PYTHONPATH=/tmp/if2453/ci/r10/src:/tmp/if2453/pilot_s1/pydeps --env OMP_NUM_THREADS="$THREADS" \
  --env M09A_THREADS="$THREADS" --env M09A_EXTRA="${M09A_EXTRA:-}" --pwd "$M/work" "$IMAGE" \
  python "$M/bin/timed.py" "$M/runs/$UNIT.json" -- "$@"
