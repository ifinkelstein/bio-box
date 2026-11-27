#!/bin/bash
# Gentle post-retry20 download chain for `tower` (card C-85, DEC-88 (1), (2)).
#
# Waits for an in-flight `tower` job PID (default usage: the PID of `after_chain_retry20.sh`) to exit cleanly, then
# runs the prepared `r26.11.6+` job lists strictly one at a time under `ionice -c 3` and `nice -n 19`:
#   - multi-file lists (`23`, `24`, `25`, `27`): `--jobs 2 --pause 1` (at most 2 connections, 1 s pause per worker)
#   - single large tarball (`26_r26_11_6_globdb_r232.tsv`): `--jobs 1 --pause 1`
# Checks `$DEST/STOP` before waiting and before every list so `touch $DEST/STOP` halts the chain cleanly.
#
# Usage (on `tower`, launched by the orchestrating session after review — NOT launched by card C-85):
#   bash scripts/orchestrator/after_retry20_gentle_chain.sh [--dest DIR] [--sources FILE] [--lists-dir DIR] [--dry-run] WAIT_PID
set -u
export LC_ALL=C

die() { echo "after_retry20_gentle_chain: $*" >&2; exit 2; }

DEST=/mnt/user/research-bulk/bio-box/sources_incoming
HERE=${0%/*}
[ "$HERE" != "$0" ] || HERE=.
SOURCES=$HERE/../../sources.yaml
LISTS_DIR=$HERE/tower_jobs
DRY_RUN=0

while [ $# -gt 0 ]; do
  case "$1" in
    --dest) [ $# -ge 2 ] || die "--dest needs a value"; DEST=$2; shift 2 ;;
    --sources) [ $# -ge 2 ] || die "--sources needs a value"; SOURCES=$2; shift 2 ;;
    --lists-dir) [ $# -ge 2 ] || die "--lists-dir needs a value"; LISTS_DIR=$2; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    --) shift; break ;;
    -*) die "unknown option $1" ;;
    *) break ;;
  esac
done

[ $# -eq 1 ] || die "usage: after_retry20_gentle_chain.sh [--dest DIR] [--sources FILE] [--lists-dir DIR] [--dry-run] WAIT_PID"
WAIT_PID=$1
case "$WAIT_PID" in ''|*[!0-9]*|0) die "WAIT_PID must be a positive integer pid: $WAIT_PID" ;; esac
DEST=${DEST%/}
case "$DEST" in /*) ;; *) die "--dest must be an absolute path: $DEST" ;; esac

if [ -e "$DEST/STOP" ]; then
  echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) STOP file found in $DEST before waiting: gentle chain ended cleanly"
  exit 0
fi

echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) waiting for pid $WAIT_PID to finish"
if [ "$DRY_RUN" = 0 ]; then
  if tail --pid="$WAIT_PID" -f /dev/null 2>/dev/null; then
    :
  else
    while kill -0 "$WAIT_PID" 2>/dev/null; do
      sleep 5
    done
  fi
fi

# Ordered specs: `<jobs> <pause> <list_filename>` (DEC-88 (1): <= 2 jobs with 1 s pause, 1 job for GlobDB r232 tarball)
SPECS=(
  "2 1 23_r26_11_6_asg1_fna.tsv"
  "2 1 24_r26_11_6_dsv70_prjna1244896.tsv"
  "2 1 25_r26_11_6_mgnify_nonhuman_reps.tsv"
  "1 1 26_r26_11_6_globdb_r232.tsv"
  "2 1 27_r26_11_6_gtdb_r232_delta.tsv"
)

OVERALL_RC=0
for spec in "${SPECS[@]}"; do
  set -- $spec
  jobs=$1
  pause=$2
  list_name=$3
  list_path=$LISTS_DIR/$list_name

  if [ -e "$DEST/STOP" ]; then
    echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) STOP file found in $DEST: gentle chain ended cleanly before $list_name"
    exit 0
  fi

  if [ ! -f "$list_path" ]; then
    echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) skip absent list $list_name ($list_path)"
    continue
  fi

  echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) start $list_name (ionice -c 3 nice -n 19, --jobs $jobs --pause $pause)"
  if [ "$DRY_RUN" = 1 ]; then
    echo "DRY-RUN: ionice -c 3 nice -n 19 /bin/bash $HERE/fetch_to_tower.sh --dest $DEST --sources $SOURCES --jobs $jobs --pause $pause $list_path"
    rc=0
  else
    ionice -c 3 nice -n 19 /bin/bash "$HERE/fetch_to_tower.sh" \
      --dest "$DEST" \
      --sources "$SOURCES" \
      --jobs "$jobs" \
      --pause "$pause" \
      "$list_path"
    rc=$?
  fi
  echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) end $list_name exit $rc"
  [ "$rc" = 0 ] || OVERALL_RC=1
done

echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) gentle chain done overall_rc=$OVERALL_RC"
exit "$OVERALL_RC"
