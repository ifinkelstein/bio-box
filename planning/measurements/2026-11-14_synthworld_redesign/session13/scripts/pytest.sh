#!/bin/sh
# pytest.sh TREE BASETEMP_NAME [pytest args...]: pytest in the test image on a tree mirrored to the pod (card F-02
# session 13), started by run.sh as a detached unit whose CPUQuota caps the threads; the null module's seed comes
# from BIOBOX_SYNTHWORLD_NULL_SEED, passed as the first pytest argument "seed=<n>" when given.
T=$1
B=$2
shift 2
cd "$T" || exit 1
case "$1" in
  seed=*) export BIOBOX_SYNTHWORLD_NULL_SEED="${1#seed=}"; shift ;;
esac
exec python -m pytest -q -s -p no:cacheprovider --basetemp="/tmp/if2453/ci-tmp/F-02-13-$B" "$@"
