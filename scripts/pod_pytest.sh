#!/bin/sh
# Run pytest for a working tree ON THE POD, inside the pinned image (AGENTS.md §5.5, §7).
# No tests or measurements run on the laptop.
#
#   scripts/pod_pytest.sh <name> [--golden] [--shards N] [pytest args...]
#
# <name> is a short label (the task ID). The working tree, committed or not, is mirrored to
# $BIOBOX_POD_CI/<name> on the pod and pytest runs there with that tree's src/ first on the path.
# --golden also sets BIOBOX_GOLDEN_DIR, so ring R3 runs.
# --shards N (only without explicit test paths) splits the test files over N pytest processes
#   that run at once; each prints its own summary line and the script fails if any shard fails.
#   Timing (measured, wilkcomp02, 2026-11-13, card C-29): 2 shards of 2 threads took 3,047 s and 4,761 s
#   (51 and 80 minutes) while a 32-thread build ran beside them. The shard count must fit the account's
#   thread cap (AGENTS.md §7.2, a total over everything the account runs); 8 shards do not fit by day.
# The image is the test image built from `main` at 64f1bf0 (card C-29): the locked Python environment (with
#   matplotlib and zombi2), MMseqs2, HMMER, Infernal and skani, so tests that skip on a laptop run here. Only the
#   tree's src/ goes on PYTHONPATH: a directory put before it would shadow the image's locked packages.
# BIOBOX_POD_IMAGE: another image (a .sif path on the pod). BIOBOX_POD_PYDEPS: directories (colon-separated) put
#   after the tree's src/ on PYTHONPATH, for an image that lacks a package; none by default.
# BIOBOX_PEER_MANIFEST (a path on the pod, under /tmp/if2453) is passed into the image when set: the
#   cross-host test (card D-04) compares the fixture's manifest hash with the one another host recorded.
set -eu
NAME=$1
shift
GOLDEN_ENV=""
SHARDS=1
while [ $# -gt 0 ]; do
  case "$1" in
    --golden) GOLDEN_ENV="--env BIOBOX_GOLDEN_DIR=${BIOBOX_POD_GOLDEN:-/tmp/if2453/golden}"; shift ;;
    --shards) SHARDS=$2; shift 2 ;;
    *) break ;;
  esac
done
HOST=${BIOBOX_POD:-if2453@wilkcomp02.ccbb.utexas.edu}
CI=${BIOBOX_POD_CI:-/tmp/if2453/ci}
IMAGE=${BIOBOX_POD_IMAGE:-/stor/home/if2453/bio-box/images/bio-box-test_64f1bf0.sif}
PYDEPS=${BIOBOX_POD_PYDEPS:-}
PEER_ENV=""
[ -z "${BIOBOX_PEER_MANIFEST:-}" ] || PEER_ENV="--env BIOBOX_PEER_MANIFEST=$BIOBOX_PEER_MANIFEST"
TOP=${BIOBOX_TREE:-$(git rev-parse --show-toplevel)}  # BIOBOX_TREE: test a tree other than the current directory's
DST="$CI/$NAME"
ARGS=""
for a in "$@"; do
  ARGS="$ARGS '$(printf %s "$a" | sed "s/'/'\\\\''/g")'"
done
ssh -o BatchMode=yes "$HOST" "mkdir -p '$DST' /tmp/if2453/ci-tmp"
rsync -a --delete -e "ssh -o BatchMode=yes" \
  --exclude .git --exclude .venv --exclude .hypothesis --exclude .pytest_cache --exclude mutants \
  --exclude .claude --exclude __pycache__ "$TOP/" "$HOST:$DST/"

# /sys/fs/cgroup is bound read-only so that run-step and the search runner inside the image see
# the scope's memory cap (card C-07: inside the image the directory is otherwise empty).
run() {  # run <basetemp suffix> <pytest args string>
  ssh -o BatchMode=yes "$HOST" "systemd-run --user --scope --quiet -p MemoryMax=32G -p MemorySwapMax=0 nice -n 19 \
    singularity exec --cleanenv --env PYTHONPATH=$DST/src${PYDEPS:+:$PYDEPS} $GOLDEN_ENV $PEER_ENV \
    --env OMP_NUM_THREADS=2 --bind /tmp/if2453 --bind /sys/fs/cgroup:/sys/fs/cgroup:ro --pwd '$DST' '$IMAGE' \
    python -m pytest -q -p no:cacheprovider --basetemp=/tmp/if2453/ci-tmp/$NAME-$1 $2"
}

if [ "$SHARDS" -le 1 ] || [ -n "$ARGS" ]; then
  run 0 "$ARGS"
  exit $?
fi

OUT=$(mktemp -d)
FILES=$(cd "$TOP" && find tests -name 'test_*.py' | sort)
i=0
for f in $FILES; do
  k=$((i % SHARDS))
  printf "%s " "$f" >> "$OUT/files.$k"
  i=$((i + 1))
done
PIDS=""
for k in $(seq 0 $((SHARDS - 1))); do
  [ -f "$OUT/files.$k" ] || continue
  ( set +e; run "$k" "$(cat "$OUT/files.$k")" > "$OUT/log.$k" 2>&1; echo $? > "$OUT/rc.$k" ) &
  PIDS="$PIDS $!"
done
for p in $PIDS; do
  wait "$p" || true
done
RC=0
for k in $(seq 0 $((SHARDS - 1))); do
  [ -f "$OUT/log.$k" ] || continue
  if [ "$(cat "$OUT/rc.$k")" != "0" ]; then
    RC=1
    echo "=== shard $k FAILED ==="
    tail -n 60 "$OUT/log.$k"
  else
    echo "shard $k: $(tail -n 1 "$OUT/log.$k")"
  fi
done
rm -rf "$OUT"
exit $RC
