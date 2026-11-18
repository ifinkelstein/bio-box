#!/bin/bash
# Card C-39: the inputs of stage V1 (UHGV vOTU representatives) and what they share, on the pod. Nothing is written in
# the data root (bound read-only into the image) or under /stor. Copies of this script and the c39_*.py scripts,
# make_inputs.py and s2_check_inputs.py run from /tmp/if2453/s3/bin; downloads are in /tmp/if2453/s3/dl (bound
# read-only). Each mode in a small unit named c39-<mode>, for example
#
#   systemd-run --user --unit=c39-inputs --collect -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=100% \
#     nice -n 19 ionice -c 2 -n 7 bash /tmp/if2453/s3/bin/c39_run.sh inputs
#
# Modes:
#   inputs   c39_inputs.py uhgv (the FASTA decompressed into IN, genomes.json, objects.json, record tables) with
#            --root ROOT, report to C/inputs/report.json. Refuses when IN exists.
#   rerun    the inputs again into C/rerun/V1, each file compared with IN's (cmp), then removed: the determinism check.
#   check    s2_check_inputs.py V1 (the inputs step's parsing and checks; writes nothing) into C/inputs/check.json,
#            and the sha256 of every file IN holds into C/inputs/SUMS.
#   overlap  c39_overlap.py: by contig_uid against every stage's contig_index and genome_input (read only), into
#            C/overlap.
#   ani      c39_ani.py: skani at the vOTU threshold on a sample of V1 against S1 and H1 (THREADS, default 8), into
#            C/ani.
#   est      c39_est.py: the code's own estimates and a gene-calling sample (one thread), into C/est.
#   vire     c39_vire.py: VIRE's table summarised, into C/vire.
set -u

C=/tmp/if2453/s3/work
ROOT=/tmp/if2453/r10/root
TREE=${TREE:-/tmp/if2453/ci/C-39}
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
DL=/tmp/if2453/s3/dl
IN=/tmp/if2453/s2/inputs/V1
BIN=/tmp/if2453/s3/bin
LOG=/tmp/if2453/s3/logs/run.log
S1FA=/tmp/if2453/r10/inputs/S1/phage_first34076.fa
H1FA=/tmp/if2453/s2/genomes/H1/HPGC.7230.genomes.final.fasta
mkdir -p "$C"

log() { echo "$(date '+%F %T %Z') $*" >> "$LOG"; }
fail() { log "FAILED: $*"; exit 1; }
for f in "$ROOT/STOP" /tmp/if2453/s3/STOP; do [ -e "$f" ] && { log "STOP present: $f; not starting"; exit 3; }; done

py() {  # python in the runtime image, the tree's src on the path, the root and the downloads bound read-only
  singularity exec --cleanenv --env PYTHONPATH="$TREE/src" --env OMP_NUM_THREADS=1 --env PYTHONUNBUFFERED=1 \
    --pwd "$BIN" --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$DL:$DL:ro" "$IMAGE" python "$@"
}

step() {  # name out_dir command...: run, log the exit and the head of the output
  local name=$1 out=$2 rc
  shift 2
  mkdir -p "$out"
  log "$name: start"
  "$@" > "$out/out.json" 2> "$out/err.txt"
  rc=$?
  log "$name: exit $rc: $(cut -c1-600 "$out/out.json") $(tail -n 2 "$out/err.txt")"
  [ "$rc" = 0 ] || fail "$name exit $rc"
}

MODE=${1:-}
log "c39: mode '$MODE' (pid $$), tree $TREE, image $IMAGE"
case "$MODE" in
  inputs)
    [ -e "$IN" ] && fail "$IN exists"
    step inputs "$C/inputs" py "$BIN/c39_inputs.py" uhgv "$DL/uhgv/votus_full.fna.zst" \
      "$DL/uhgv/votus_metadata_extended.tsv.zst" "$IN" 2026-11-16T19:53:06Z --root "$ROOT"
    cp "$C/inputs/out.json" "$C/inputs/report.json" ;;
  rerun)  # the same inputs again into scratch, compared file by file with IN (determinism)
    rm -rf "$C/rerun"
    mkdir -p "$C/rerun"
    step rerun "$C/rerun/log" py "$BIN/c39_inputs.py" uhgv "$DL/uhgv/votus_full.fna.zst" \
      "$DL/uhgv/votus_metadata_extended.tsv.zst" "$C/rerun/V1" 2026-11-16T19:53:06Z --root "$ROOT"
    for f in genomes.json objects.json records.tsv repeated_records.tsv shared_records.tsv uhgv_v1.0_votus_full.fna; do
      if [ "$f" = objects.json ]; then
        log "rerun: objects.json names its own path; compared after the path: $(sed "s#$C/rerun/V1#$IN#" \
          "$C/rerun/V1/$f" | cmp - "$IN/$f" && echo equal)"
      else
        log "rerun: $f $(cmp "$C/rerun/V1/$f" "$IN/$f" && echo equal)"
      fi
    done
    rm -rf "$C/rerun/V1" ;;
  check)
    step check "$C/check" py "$BIN/s2_check_inputs.py" V1 "$IN/genomes.json" "$IN/objects.json"
    (cd "$IN" && sha256sum genomes.json objects.json records.tsv repeated_records.tsv shared_records.tsv \
      uhgv_v1.0_votus_full.fna) > "$C/check/SUMS"
    log "check: SUMS: $(tr '\n' ' ' < "$C/check/SUMS")" ;;
  overlap)
    step overlap "$C/overlap" py "$BIN/c39_overlap.py" --root "$ROOT" --inputs "$IN" \
      --uhgv-meta "$DL/uhgv/uhgv_metadata.tsv.zst" --out "$C/overlap" ;;
  ani)
    step ani "$C/ani" py "$BIN/c39_ani.py" --inputs "$IN" --s1 "$S1FA" --h1 "$H1FA" --out "$C/ani" \
      --threads "${THREADS:-8}" --every "${EVERY:-20}" ${REUSE:+--reuse} ;;
  est)
    step est "$C/est" py "$BIN/c39_est.py" --root "$ROOT" --inputs "$IN" --faa "$DL/uhgv/votus_full.faa.zst" \
      --out "$C/est" --every "${EVERY:-30}" ;;
  vire)
    step vire "$C/vire" py "$BIN/c39_vire.py" "$DL/vire/vire_all_samples_metadata.txt" "$C/vire" ;;
  *) echo "usage: c39_run.sh inputs|rerun|check|overlap|ani|est|vire" >&2; exit 2 ;;
esac
log "c39: mode '$MODE' done"
