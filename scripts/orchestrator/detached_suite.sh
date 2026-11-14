#!/bin/bash
# Run the test suite of an already-mirrored tree ON THE POD as detached units, so a dropped ssh
# connection (a sleeping laptop) cannot hang it. Run this script on the pod:
#
#   bash detached_suite.sh <name> <shards> [pytest args...]
#   bash detached_suite.sh --own <name> <shards> [pytest args...]
#
# <name> is the mirror under /tmp/if2453/ci/<name> that scripts/pod_pytest.sh made (run
# `pod_pytest.sh <name> tests/r0/test_runbook.py` first to sync the tree); its src/ alone is put on
# PYTHONPATH. --own runs the image's own suite instead (the image check, RUNBOOK section 8): the
# source and tests built into the image under /opt/bio-box, with no PYTHONPATH; <name> only names
# the units and logs. Each shard is a unit ci-<name>-<k>; its log is
# /tmp/if2453/ci-tmp/<name>-detached/log.<k>, whose last line is pytest's summary (the skip reasons
# just before it), and junit.<k>.xml beside it has each test's outcome. Poll with:
# systemctl --user is-active 'ci-<name>-*'; tail -n 1 .../log.*
#
# The image is the test image built from `main` at faec1f4 (card C-32; it holds MinCED), whose
# environment is the locked one. As in pod_pytest.sh: BIOBOX_POD_IMAGE names another image; BIOBOX_POD_PYDEPS adds
# directories after the tree's src/ (none by default: one put on the path would shadow the image's
# packages); BIOBOX_PEER_MANIFEST (a file under /tmp/if2453) is passed to the cross-host test (D-04).
set -eu
OWN=0
if [ "${1:-}" = --own ]; then
  OWN=1
  shift
fi
NAME=$1
N=$2
shift 2
OUT=/tmp/if2453/ci-tmp/$NAME-detached
IMAGE=${BIOBOX_POD_IMAGE:-/stor/home/if2453/bio-box/images/bio-box-test_faec1f4.sif}
PYDEPS=${BIOBOX_POD_PYDEPS:-}
ENVS=(--env OMP_NUM_THREADS=2)
[ -z "${BIOBOX_PEER_MANIFEST:-}" ] || ENVS+=(--env "BIOBOX_PEER_MANIFEST=$BIOBOX_PEER_MANIFEST")
case "$OWN" in
  1) DST=/opt/bio-box
     FILES=$(singularity exec --cleanenv --pwd "$DST" "$IMAGE" find tests -name 'test_*.py' | sort) ;;
  *) DST=/tmp/if2453/ci/$NAME
     ENVS+=(--env "PYTHONPATH=$DST/src${PYDEPS:+:$PYDEPS}")
     FILES=$(cd "$DST" && find tests -name 'test_*.py' | sort) ;;
esac
mkdir -p "$OUT"
rm -f "$OUT"/files.* "$OUT"/log.* "$OUT"/junit.*
i=0
for f in $FILES; do
  printf '%s ' "$f" >> "$OUT/files.$((i % N))"
  i=$((i + 1))
done
for k in $(seq 0 $((N - 1))); do
  # shellcheck disable=SC2046
  systemd-run --user --unit="ci-$NAME-$k" --collect --quiet -p MemoryMax=32G -p MemorySwapMax=0 \
    -p StandardOutput="file:$OUT/log.$k" -p StandardError="file:$OUT/log.$k" \
    nice -n 19 singularity exec --cleanenv "${ENVS[@]}" \
    --bind /tmp/if2453 --bind /sys/fs/cgroup:/sys/fs/cgroup:ro --pwd "$DST" "$IMAGE" \
    python -m pytest -q -p no:cacheprovider -rs "--junitxml=$OUT/junit.$k.xml" "--basetemp=/tmp/if2453/ci-tmp/$NAME-$k" "$@" \
    $(cat "$OUT/files.$k")
done
echo "started $N units ci-$NAME-0..$((N - 1)) in $IMAGE; logs in $OUT"
