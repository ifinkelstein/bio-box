#!/bin/bash
# Short commands of the 50 Gbp build (planning/measurements/2026-11-11_build_r50.md), in the image, in their own
# capped scope at nice 19 (RUNBOOK §0's bio-box and biopy), with the 10 Gbp root bound read-only on top of
# /tmp/if2453 so nothing run here can write into it.
#
#   r50_bb.sh cli <bio-box args...>     python -m biobox.cli ...
#   r50_bb.sh py <script> <args...>     python <script> ...
#
# MEM (default 32G) and THREADS (default 4, OMP_NUM_THREADS and the scope's CPUQuota) may be set in the
# environment. Long jobs go through scripts/pod_run.sh instead.
set -u
TREE=/tmp/if2453/ci/r50
IMAGE=/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif
PYPATH=$TREE/src:/tmp/if2453/pilot_s1/pydeps
PFAMDIR=/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2
MEM=${MEM:-32G}
THREADS=${THREADS:-4}
RO=${RO:-/tmp/if2453/r10/root}
what=$1
shift
case "$what" in
  cli) set -- python -m biobox.cli "$@" ;;
  py) set -- python "$@" ;;
  *) echo "usage: r50_bb.sh cli|py ..." >&2; exit 2 ;;
esac
ro=()
[ "$RO" = none ] || ro=(--bind "$RO:$RO:ro")
exec systemd-run --user --scope --quiet -p MemoryMax="$MEM" -p MemorySwapMax=0 -p CPUQuota="${THREADS}00%" \
  nice -n 19 singularity exec --cleanenv --env PYTHONPATH="$PYPATH" --env OMP_NUM_THREADS="$THREADS" \
  --env PYTHONUNBUFFERED=1 --bind /tmp/if2453 "${ro[@]}" --bind "$PFAMDIR:$PFAMDIR:ro" \
  --bind /sys/fs/cgroup:/sys/fs/cgroup:ro "$IMAGE" "$@"
