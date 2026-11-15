#!/bin/bash
# Card C-26 part 1: the inputs of the PLSDB stage PL1 and what they add, on the pod; nothing is written in the data
# root (bound read-only into the image). Copies of this script, make_inputs.py, c26_analyze.py and s2_check_inputs.py
# run from /tmp/if2453/c26/bin. Each mode in a small unit, for example
#
#   systemd-run --user --unit=c26-inputs --collect -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=100% \
#     -p StandardOutput=append:/tmp/if2453/c26/driver.out -p StandardError=append:/tmp/if2453/c26/driver.out \
#     nice -n 19 ionice -c 2 -n 7 bash /tmp/if2453/c26/bin/c26_inputs.sh inputs
#
# Modes:
#   inputs   Checks the PLSDB FASTA's sha256 against sources.yaml, links it (a hard link, same file system) into
#            IN as plsdb_2024_05_31_v2_sequences.fasta, the name of its copy on /stor, and runs
#            `make_inputs.py plasmid FA nuccore.csv IN --root ROOT`; its report goes to C/inputs/report.json.
#            Refuses when IN already holds genomes.json.
#   check    s2_check_inputs.py PL1 (the inputs step's parsing and checks, writes nothing) into C/inputs/check.json,
#            and the sha256 of every file IN holds into C/inputs/SUMS.
#   analyze  c26_analyze.py into C/analyze (read-only on ROOT and the metadata).
#   rate     c26_rate.py: one thread, every EVERY-th (default 30) record the stage stores called with the code's
#            call path, into C/rate (the gene-calling rate in meta mode, for the estimate of the gene-call step).
set -u

C=/tmp/if2453/c26
ROOT=/tmp/if2453/r10/root
TREE=${TREE:-/tmp/if2453/ci/C-26}
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
S2=/tmp/if2453/s2
IN=$S2/inputs/PL1
SRC=$S2/dl/plsdb_download_fasta
SRC_SHA=e654c106d8a15a2d5d55a8b680087377cb84e654e63d41bca181fc69a23e587a
FA=$IN/plsdb_2024_05_31_v2_sequences.fasta
META=$S2/plsdb/meta
BIN=$C/bin
LOG=$C/driver.log
mkdir -p "$C/inputs" "$C/analyze"

log() { echo "$(date '+%F %T %Z') $*" >> "$LOG"; }
fail() { log "FAILED: $*"; exit 1; }

py() {  # python in the runtime image, the tree's src on the path, the root and the metadata bound read-only
  singularity exec --cleanenv --env PYTHONPATH="$TREE/src" --env OMP_NUM_THREADS=1 --env PYTHONUNBUFFERED=1 \
    --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$META:$META:ro" --bind "$S2/dl:$S2/dl:ro" "$IMAGE" python "$@"
}

inputs() {
  local sha rc
  [ -e "$IN/genomes.json" ] && fail "$IN/genomes.json exists"
  [ -e "$ROOT/STOP" ] && log "note: $ROOT/STOP exists (nothing here writes the root)"
  sha=$(sha256sum "$SRC" | cut -d' ' -f1)
  log "inputs: sha256 of $SRC: $sha"
  [ "$sha" = "$SRC_SHA" ] || fail "sha256 of $SRC is not the one sources.yaml records"
  mkdir -p "$IN"
  [ -e "$FA" ] || ln "$SRC" "$FA" || fail "ln $SRC $FA"
  log "inputs: $(stat -c '%n %s bytes, %h links, inode %i' "$FA")"
  py "$BIN/make_inputs.py" plasmid "$FA" "$META/nuccore.csv" "$IN" --root "$ROOT" \
    > "$C/inputs/report.json" 2> "$C/inputs/report.err"
  rc=$?
  log "inputs: make_inputs.py plasmid exit $rc: $(cut -c1-600 "$C/inputs/report.json") $(tail -n 2 "$C/inputs/report.err")"
  [ "$rc" = 0 ] || fail "make_inputs.py plasmid exit $rc"
}

check() {
  local rc
  py "$BIN/s2_check_inputs.py" PL1 "$IN/genomes.json" "$IN/objects.json" > "$C/inputs/check.json" 2> "$C/inputs/check.err"
  rc=$?
  log "check: s2_check_inputs.py exit $rc: $(tr -d '\n' < "$C/inputs/check.json" | cut -c1-600)"
  (cd "$IN" && sha256sum genomes.json objects.json repeated_records.tsv shared_records.tsv) > "$C/inputs/SUMS"
  log "check: SUMS: $(tr '\n' ' ' < "$C/inputs/SUMS")"
  [ "$rc" = 0 ] || fail "s2_check_inputs.py exit $rc"
}

analyze() {
  local rc
  py "$BIN/c26_analyze.py" --root "$ROOT" --inputs "$IN" --meta "$META" --out "$C/analyze" \
    > "$C/analyze/out.json" 2> "$C/analyze/err.txt"
  rc=$?
  log "analyze: exit $rc: $(cut -c1-400 "$C/analyze/out.json") $(tail -n 2 "$C/analyze/err.txt")"
  [ "$rc" = 0 ] || fail "c26_analyze.py exit $rc"
}

rate() {
  local rc
  mkdir -p "$C/rate"
  py "$BIN/c26_rate.py" "$FA" "$IN" "$C/rate/rate.json" "${EVERY:-30}" > "$C/rate/out.json" 2> "$C/rate/err.txt"
  rc=$?
  log "rate: exit $rc: $(cut -c1-600 "$C/rate/out.json") $(tail -n 2 "$C/rate/err.txt")"
  [ "$rc" = 0 ] || fail "c26_rate.py exit $rc"
}

MODE=${1:-}
log "c26 inputs: mode '$MODE' (pid $$), tree $TREE, image $IMAGE"
case "$MODE" in
  inputs) inputs ;;
  check) check ;;
  analyze) analyze ;;
  rate) rate ;;
  *) echo "usage: c26_inputs.sh inputs|check|analyze|rate" >&2; exit 2 ;;
esac
log "c26 inputs: mode '$MODE' done"
