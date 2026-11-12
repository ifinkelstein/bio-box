#!/bin/bash
# Launch one card C-21f job on the pod, detached, inside the image (run ON THE POD HOST).
#
#   bash launch.sh <name> <command...>
#
# Unit c21f-<name>: MemoryMax=16G, MemorySwapMax=0, CPUQuota=200%, nice 19. The data root is bound read-only; the
# work directory /tmp/if2453/c21f and the tree /tmp/if2453/ci/c21f (rsync of the branch) are the only other binds
# that matter. The command runs under annotation_cost/timed.py, which writes runs/c21f-<name>.json (wall, CPU-minutes,
# peak RSS); stdout and stderr go to logs/c21f-<name>.log, whose last line is "launch: c21f-<name> exited <code>".
#
# Refuses to start when a STOP file is in the data root or the work directory, when two c21f units are already
# active (the brief's limit), or when the account's CPU use (sum of `ps -u if2453 -o pcpu=`, ps excluded) plus 200%
# would pass the cap: 4000% from 08:00 to 20:00 CST on weekdays, 5600% otherwise (AGENTS §7.2; the cap is a total
# over every job of the account). It does not wait: the caller retries later.
set -euo pipefail
W=/tmp/if2453/c21f
ROOT=/tmp/if2453/r10/root
TREE=/tmp/if2453/ci/c21f
IMAGE=/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif
PYDEPS=/tmp/if2453/pilot_s1/pydeps
[ $# -ge 2 ] || { echo "usage: launch.sh <name> <command...>" >&2; exit 2; }
UNIT=c21f-$1
shift
mkdir -p "$W/logs" "$W/runs" "$W/work" "$W/tmp"
for s in "$ROOT/STOP" "$W/STOP"; do
  [ ! -e "$s" ] || { echo "launch: $s exists; not starting $UNIT" >&2; exit 3; }
done
active=$(systemctl --user list-units --state=active,activating --no-legend 'c21f-*' | wc -l)
[ "$active" -lt 2 ] || { echo "launch: $active c21f units active; not starting $UNIT" >&2; exit 3; }
used=$(ps -u if2453 -o pcpu=,comm= | awk '$2 != "ps" {s+=$1} END {printf "%d", s}')
h=$((10#$(TZ=America/Chicago date +%H)))
dow=$(TZ=America/Chicago date +%u)
if [ "$dow" -le 5 ] && [ "$h" -ge 8 ] && [ "$h" -lt 20 ]; then cap=4000; else cap=5600; fi
[ $((used + 200)) -le "$cap" ] || { echo "launch: account ${used}% + 200% > ${cap}%; not starting $UNIT" >&2; exit 3; }
echo "launch: $UNIT started at $(date -Is): account ${used}% + 200% <= ${cap}%" | tee -a "$W/logs/launch.log"
systemd-run --user --unit="$UNIT" --collect --quiet --same-dir \
  -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% \
  -p StandardOutput="append:$W/logs/$UNIT.log" -p StandardError="append:$W/logs/$UNIT.log" \
  /bin/sh -c 'unit=$1; shift; "$@"; rc=$?; echo "launch: $unit exited $rc"; exit $rc' sh "$UNIT" \
  nice -n 19 singularity exec --cleanenv \
  --bind "$ROOT:$ROOT:ro" --bind "$W:$W" --bind "$TREE:$TREE:ro" --bind "$PYDEPS:$PYDEPS:ro" \
  --env PYTHONPATH="$TREE/src:$PYDEPS" --env OMP_NUM_THREADS=1 --pwd "$W/work" "$IMAGE" \
  python "$TREE/scripts/stress/annotation_cost/timed.py" "$W/runs/$UNIT.json" -- "$@"
