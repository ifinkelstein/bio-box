#!/bin/bash
# Run the test suite of an already-mirrored tree ON THE POD as detached units, so a dropped ssh
# connection (a sleeping laptop) cannot hang it. Run this script on the pod:
#
#   bash detached_suite.sh <name> <shards>
#
# <name> is the mirror under /tmp/if2453/ci/<name> that scripts/pod_pytest.sh made (run
# `pod_pytest.sh <name> tests/r0/test_runbook.py` first to sync the tree). Each shard is a unit
# ci-<name>-<k>; its log is /tmp/if2453/ci-tmp/<name>-detached/log.<k>, whose last line is pytest's
# summary. Poll with: systemctl --user is-active 'ci-<name>-*'; tail -n 1 .../log.*
set -eu
NAME=$1
N=$2
DST=/tmp/if2453/ci/$NAME
OUT=/tmp/if2453/ci-tmp/$NAME-detached
IMAGE=/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif
PYDEPS=/tmp/if2453/pilot_s1/pydeps
mkdir -p "$OUT"
rm -f "$OUT"/files.* "$OUT"/log.*
cd "$DST"
i=0
for f in $(find tests -name 'test_*.py' | sort); do
  printf '%s ' "$f" >> "$OUT/files.$((i % N))"
  i=$((i + 1))
done
for k in $(seq 0 $((N - 1))); do
  # shellcheck disable=SC2046
  systemd-run --user --unit="ci-$NAME-$k" --collect --quiet -p MemoryMax=32G -p MemorySwapMax=0 \
    -p StandardOutput="file:$OUT/log.$k" -p StandardError="file:$OUT/log.$k" \
    nice -n 19 singularity exec --cleanenv --env "PYTHONPATH=$DST/src:$PYDEPS" --env OMP_NUM_THREADS=2 \
    --bind /tmp/if2453 --bind /sys/fs/cgroup:/sys/fs/cgroup:ro --pwd "$DST" "$IMAGE" \
    python -m pytest -q -p no:cacheprovider "--basetemp=/tmp/if2453/ci-tmp/$NAME-$k" $(cat "$OUT/files.$k")
done
echo "started $N units ci-$NAME-0..$((N - 1)); logs in $OUT"
