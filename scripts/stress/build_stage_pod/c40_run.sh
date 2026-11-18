#!/bin/bash
# Card C-40: stage PL2's inputs and what they share, on the pod; nothing is written in the data root (bound read-only
# into the image) or under /stor. Copies of this script and the c40_*.py scripts run from /tmp/if2453/s3/plasmids/bin;
# make_inputs.py and s2_check_inputs.py from the tree synced by `pod_pytest.sh C-40`. Each mode in a small unit, e.g.
#
#   systemd-run --user --unit=c40-inputs --collect -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=100% \
#     nice -n 19 ionice -c 2 -n 7 bash /tmp/if2453/s3/plasmids/bin/c40_run.sh inputs
#
# Modes:
#   concat SRC  The archive's sha256 checked against dl/SUMS, then c40_concat.py into fasta/SRC/ (new, empty), with
#               python -I (no environment, no script directory on the path): the archive is read as a stream, never
#               unpacked. SRC "gut" instead checks the gut table's sha256 and decompresses it into fasta/gut/.
#   spec        c40_spec.py into work/spec.json.
#   inputs      make_inputs.py plasmid-sources work/spec.json IN --root ROOT (refuses when IN exists).
#   check       s2_check_inputs.py PL2 (the inputs step's parser and checks, writes nothing) and the sha256 of IN's files.
#   analyze     c40_analyze.py into work/analyze.
#   est         c40_est.py into work/est (one thread; EVERY, default 50).
# A STOP file in the data root stops every mode before it starts (it is the release build's).
set -u

W=/tmp/if2453/s3/plasmids
ROOT=/tmp/if2453/r10/root
TREE=${TREE:-/tmp/if2453/ci/C-40}
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
IN=/tmp/if2453/s2/inputs/PL2
BIN=$W/bin
DL=$W/dl
FA=$W/fasta
LOG=$W/logs/driver.log
mkdir -p "$W/work" "$FA" "$W/logs"

log() { echo "$(date '+%F %T %Z') $*" >> "$LOG"; }
fail() { log "FAILED: $*"; exit 1; }

py() {  # python in the runtime image, the tree's src on the path, the root and the downloads bound read-only
  singularity exec --cleanenv --env PYTHONPATH="$TREE/src" --env OMP_NUM_THREADS=1 --env PYTHONUNBUFFERED=1 \
    --pwd "$BIN" --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$DL:$DL:ro" "$IMAGE" python "$@"
}

sha_ok() {  # file relative to DL: its sha256 equals the one c40_fetch.sh recorded
  local want got
  want=$(awk -v f="$1" '$1 == f {print $7}' "$DL/SUMS")
  got=$(sha256sum "$DL/$1" | cut -d' ' -f1)
  log "sha256 $1: $got (fetched $want)"
  [ -n "$want" ] && [ "$want" = "$got" ]
}

concat() {
  local src=$1 rc
  if [ "$src" = gut ]; then
    sha_ok gut/predicted_plasmids_summary.txt.gz || fail "gut table sha256"
    sha_ok gut/predicted_plasmids.fa.gz || fail "gut FASTA sha256"
    mkdir -p "$FA/gut"
    [ -e "$FA/gut/predicted_plasmids_summary.txt" ] && fail "$FA/gut/predicted_plasmids_summary.txt exists"
    zcat "$DL/gut/predicted_plasmids_summary.txt.gz" > "$FA/gut/.summary.tmp" || fail "zcat"
    mv "$FA/gut/.summary.tmp" "$FA/gut/predicted_plasmids_summary.txt"
    log "concat gut: $(wc -l < "$FA/gut/predicted_plasmids_summary.txt") lines"
    return
  fi
  sha_ok "plasmidscope/$src.fasta.tar.gz" || fail "$src archive sha256"
  sha_ok "plasmidscope/$src.plasmid_list.download.tsv" || fail "$src table sha256"
  [ -e "$FA/$src" ] && fail "$FA/$src exists"
  mkdir "$FA/$src"
  singularity exec --cleanenv --pwd "$BIN" --bind /tmp/if2453 --bind "$DL:$DL:ro" "$IMAGE" \
    python -I "$BIN/c40_concat.py" "$src" "$DL/plasmidscope/$src.fasta.tar.gz" "$FA/$src" \
    > "$W/logs/concat_$src.json" 2> "$W/logs/concat_$src.err"
  rc=$?
  log "concat $src: exit $rc: $(cut -c1-400 "$W/logs/concat_$src.json") $(tail -n 2 "$W/logs/concat_$src.err")"
  [ "$rc" = 0 ] || fail "c40_concat.py $src exit $rc"
}

spec() {
  py "$BIN/c40_spec.py" "$DL" "$FA" "$W/work/spec.json" > /dev/null 2> "$W/logs/spec.err" || fail "c40_spec.py"
  log "spec: $(sha256sum "$W/work/spec.json")"
}

inputs() {
  local rc
  [ -e "$IN" ] && fail "$IN exists"
  py "$TREE/scripts/stress/build_stage_pod/make_inputs.py" plasmid-sources "$W/work/spec.json" "$IN" --root "$ROOT" \
    > "$W/work/report.json" 2> "$W/work/report.err"
  rc=$?
  log "inputs: exit $rc: $(cut -c1-600 "$W/work/report.json") $(tail -n 2 "$W/work/report.err")"
  [ "$rc" = 0 ] || fail "make_inputs.py plasmid-sources exit $rc"
}

check() {
  local rc
  py "$TREE/scripts/stress/build_stage_pod/s2_check_inputs.py" PL2 "$IN/genomes.json" "$IN/objects.json" \
    > "$W/work/check.json" 2> "$W/work/check.err"
  rc=$?
  log "check: exit $rc: $(tr -d '\n' < "$W/work/check.json" | cut -c1-600)"
  (cd "$IN" && sha256sum genomes.json objects.json records.tsv skipped_records.tsv ./*.fasta) > "$W/work/SUMS"
  log "check: SUMS written"
  [ "$rc" = 0 ] || fail "s2_check_inputs.py exit $rc"
}

analyze() {
  local rc
  py "$BIN/c40_analyze.py" --root "$ROOT" --inputs "$IN" --scope-plsdb "$FA/PLSDB/PLSDB.fasta" \
    --scope-dir "$DL/plasmidscope" --gut-table "$FA/gut/predicted_plasmids_summary.txt" --out "$W/work/analyze" \
    > "$W/work/analyze.out" 2> "$W/work/analyze.err"
  rc=$?
  log "analyze: exit $rc: $(tail -n 2 "$W/work/analyze.err")"
  [ "$rc" = 0 ] || fail "c40_analyze.py exit $rc"
}

est() {
  local rc
  py "$BIN/c40_est.py" --root "$ROOT" --inputs "$IN" --profile "$TREE/profiles/wilkcomp02.yaml" \
    --every "${EVERY:-50}" --out "$W/work/est" > "$W/work/est.out" 2> "$W/work/est.err"
  rc=$?
  log "est: exit $rc: $(tail -n 2 "$W/work/est.err")"
  [ "$rc" = 0 ] || fail "c40_est.py exit $rc"
}

MODE=${1:-}
[ -e "$ROOT/STOP" ] && fail "$ROOT/STOP exists; mode '$MODE' not started"
log "c40: mode '$MODE' ${2:-} (pid $$), tree $TREE, image $IMAGE"
case "$MODE" in
  concat) concat "${2:?concat needs a source}" ;;
  spec) spec ;;
  inputs) inputs ;;
  check) check ;;
  analyze) analyze ;;
  est) est ;;
  *) echo "usage: c40_run.sh concat SRC|spec|inputs|check|analyze|est" >&2; exit 2 ;;
esac
log "c40: mode '$MODE' ${2:-} done"
