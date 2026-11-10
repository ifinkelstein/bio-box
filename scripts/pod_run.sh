#!/bin/bash
# The one host-side launcher for bio-box on the pod (card C-07). Run it ON THE POD HOST, not in
# the image: systemd-run is not visible inside the .sif, so the memory cap is opened here, around
# `singularity exec`, and run-step and the search runner inside find it in their cgroup.
#
#   bash scripts/pod_run.sh [options] <unit> <bio-box args...>
#
# It starts, detached,
#   systemd-run --user --unit=<unit> --collect -p MemoryMax=<mem> -p MemorySwapMax=0
#     -p CPUQuota=<threads*100>% -p OOMPolicy=continue
#     nice -n <nice> singularity exec --cleanenv --bind /sys/fs/cgroup:/sys/fs/cgroup:ro ...
#     <sif> python -m biobox.cli <bio-box args...>
# and returns at once. MemorySwapMax=0: the pod has swap, and a cap that can spill into swap is
# not a cap. OOMPolicy=continue: when the kernel kills the step command or a search worker at the
# cap, the unit is not torn down, so run-step or the runner lives to record the failure.
# /sys/fs/cgroup is bound read-only because inside the image it is otherwise an empty directory.
#
# Options (defaults from the profile, except the image):
#   --profile FILE   host profile (default: profiles/wilkcomp02.yaml beside this script)
#   --image SIF      the image (default: $BIOBOX_POD_IMAGE; one of the two is required)
#   --digest sha256:HEX   refuse unless the image's sha256 is this (it is computed either way
#                    and passed in as BIOBOX_IMAGE_DIGEST, which run records carry)
#   --mem SIZE       MemoryMax (default: the profile's mem_scope)
#   --threads N      CPUQuota of N CPUs (default: the profile's max_threads_night, the most any
#                    hour allows). The quota is fixed when the unit starts, but a build runs
#                    across hours, so it is only a ceiling: the limit of the hour is held inside
#                    (card C-12). build-stage reads the allowance again before each step and
#                    before each shard launch and starts nothing that would pass it, and each
#                    run-step chooses its own thread count by the hour when it starts.
#   --nice N         (default: the profile's nice_monolithic)
#   --bind SPEC      an extra `singularity --bind` (repeatable)
#   --env K=V        an extra `singularity --env` (repeatable), e.g. PYTHONPATH=<tree>/src
#   --log FILE       stdout and stderr of the unit (default: <profile tmp>/logs/<unit>.log);
#                    its last line is "pod_run: <unit> exited <code>"
# Watch it with `systemctl --user status <unit>` and the log; stop it with
# `systemctl --user stop <unit>` (run-step forwards the SIGTERM to its step).
set -euo pipefail

die() { echo "pod_run: $*" >&2; exit 2; }

HERE=$(cd "$(dirname "$0")" && pwd)
PROFILE="$HERE/../profiles/wilkcomp02.yaml"
IMAGE=${BIOBOX_POD_IMAGE:-}
DIGEST="" MEM="" THREADS="" NICE="" LOG=""
BINDS=() ENVS=()
while [ $# -gt 0 ]; do
  case "$1" in
    --profile) PROFILE=$2; shift 2 ;;
    --image) IMAGE=$2; shift 2 ;;
    --digest) DIGEST=$2; shift 2 ;;
    --mem) MEM=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    --nice) NICE=$2; shift 2 ;;
    --bind) BINDS+=(--bind "$2"); shift 2 ;;
    --env) ENVS+=(--env "$2"); shift 2 ;;
    --log) LOG=$2; shift 2 ;;
    --) shift; break ;;
    -*) die "unknown option $1" ;;
    *) break ;;
  esac
done
[ $# -ge 2 ] || die "usage: pod_run.sh [options] <unit> <bio-box args...>"
UNIT=$1
shift
[[ "$UNIT" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || die "bad unit name '$UNIT'"

prof() {  # the value of a top-level scalar key of the profile
  [ -f "$PROFILE" ] || die "profile $PROFILE not found"
  sed -n "s/^$1:[[:space:]]*\([^#[:space:]]*\).*/\1/p" "$PROFILE"
}
MEM=${MEM:-$(prof mem_scope)}
THREADS=${THREADS:-$(prof max_threads_night)}
NICE=${NICE:-$(prof nice_monolithic)}
TMPDIR_PROFILE=$(prof tmp)
LOG=${LOG:-$TMPDIR_PROFILE/logs/$UNIT.log}
[ -n "$MEM" ] || die "no memory cap: give --mem or a profile with mem_scope"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "bad thread count '$THREADS'"
[[ "$NICE" =~ ^[0-9]+$ ]] || die "bad nice value '$NICE'"
[ -n "$IMAGE" ] || die "no image: give --image or set BIOBOX_POD_IMAGE"
[ -f "$IMAGE" ] || die "image $IMAGE not found"

SUM=$(sha256sum "$IMAGE")
SUM=sha256:${SUM%% *}
if [ -n "$DIGEST" ] && [ "$DIGEST" != "$SUM" ]; then
  die "image $IMAGE is $SUM, not $DIGEST"
fi

SING=$(command -v singularity) || die "singularity not found on this host"
NICE_BIN=$(command -v nice) || die "nice not found on this host"
command -v systemd-run > /dev/null || die "systemd-run not found: run pod_run.sh on the host, not in the image"
mkdir -p "$(dirname "$LOG")"

# A small sh is the unit's main process, so the log ends with the exit code even though
# --collect unloads the unit as soon as it ends.
systemd-run --user --unit="$UNIT" --collect --quiet --same-dir \
  -p MemoryMax="$MEM" -p MemorySwapMax=0 -p CPUQuota="$((THREADS * 100))%" -p OOMPolicy=continue \
  -p StandardOutput="append:$LOG" -p StandardError="append:$LOG" \
  /bin/sh -c 'unit=$1; shift; "$@"; rc=$?; echo "pod_run: $unit exited $rc"; exit $rc' sh "$UNIT" \
  "$NICE_BIN" -n "$NICE" "$SING" exec --cleanenv --bind /sys/fs/cgroup:/sys/fs/cgroup:ro "${BINDS[@]}" \
  --env "BIOBOX_IMAGE_DIGEST=$SUM" "${ENVS[@]}" "$IMAGE" python -m biobox.cli "$@"

echo "pod_run: started $UNIT (MemoryMax=$MEM, MemorySwapMax=0, CPUQuota=$((THREADS * 100))%, nice $NICE)"
echo "pod_run: image $IMAGE $SUM"
echo "pod_run: log $LOG; watch: systemctl --user status $UNIT"
