#!/bin/bash
# Card C-48c part 2: the inputs of stages V2 (IMG/VR v4.1 HC high-quality vOTU representatives + Inovirus) and GV1
# (GVMAGs V2), on the pod. Reads the archives on /stor (bound read-only), the data root (bound read-only) and V1's
# inputs; writes only under $C (/tmp/if2453/c48c2) and the two input directories. Each mode runs in one detached unit,
# for example
#
#   systemd-run --user --unit=c48c2-look-gv --collect -p MemoryMax=32G -p MemorySwapMax=0 -p CPUQuota=100% \
#     nice -n 19 ionice -c 2 -n 7 bash /tmp/if2453/c48c2/tree/scripts/stress/build_stage_pod/c48c_v2_gv1_run.sh look-gv
#
# and leaves $C/logs/<mode>.DONE or <mode>.FAILED (with the exit code) beside the log $C/logs/run.log.
#
# Card C-48c part 3: the gv1-* modes stage GV1 (every GVMAGs V2 genome) into $GV1 after V2; their scratch is
# /tmp/if2453/c48c3 and V2's inputs are bound read-only. They read part 2's GVMAGs table from /tmp/if2453/c48c2.
set -u

MODE=${1:-}
case "$MODE" in
  gv1-*) C=/tmp/if2453/c48c3; RO_V2=(--bind /tmp/if2453/s2/inputs/V2:/tmp/if2453/s2/inputs/V2:ro) ;;
  *) C=/tmp/if2453/c48c2; RO_V2=() ;;
esac
ROOT=/tmp/if2453/r10/root
TREE=${TREE:-$C/tree}
BIN=$TREE/scripts/stress/build_stage_pod
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
JGI=/stor/work/Wilke/if2453/all-the-data/jgi_20261118
HCZIP=$JGI/jgi_img_vr_v4.1_hc_20261118.zip
METAZIP=$JGI/jgi_inovirus_umag_iuig_meta_20261118.zip
GVTAR=$JGI/gvmags_v2_20261118.tar
LOG=$C/logs/run.log
V2=/tmp/if2453/s2/inputs/V2
# the stages whose contig_index owns a contig before V2 (DEC-26), read only: every stage stored in the root but V1,
# which is being built (its staged inputs are read instead)
OWNER_STAGES=S1,P1,P2,P3,P4,P5,P6,H1,PL1,C5
V2_FILES="genomes.json objects.json records.tsv repeated_records.tsv shared_records.tsv imgvr_v4.1_hc_hq_votu_reps.fna inovirus_roux2019.fna"
GV1=/tmp/if2453/s2/inputs/GV1  # not G1 (GlobDB archaea, another stage)
GV1_FILES="genomes.json objects.json records.tsv repeated_records.tsv shared_records.tsv skipped_members.tsv gvmags_v2_all_18kgenomes.fna"
mkdir -p "$C/logs" "$C/work"

log() { echo "$(date '+%F %T %Z') $*" >> "$LOG"; }
done_() { log "$MODE: done"; echo 0 > "$C/logs/$MODE.DONE"; exit 0; }
fail() { log "FAILED: $*"; echo "$*" > "$C/logs/$MODE.FAILED"; exit 1; }
rm -f "$C/logs/$MODE.DONE" "$C/logs/$MODE.FAILED"
for f in "$ROOT/STOP" "$C/STOP"; do [ -e "$f" ] && { log "STOP present: $f; not starting $MODE"; exit 3; }; done

py() {  # python in the runtime image, the tree's src on the path, the root, the archives and V1's inputs read-only
  singularity exec --cleanenv --env PYTHONPATH="$TREE/src" --env OMP_NUM_THREADS=1 --env PYTHONUNBUFFERED=1 \
    --pwd "$BIN" --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$JGI:$JGI:ro" \
    --bind /tmp/if2453/s2/inputs/V1:/tmp/if2453/s2/inputs/V1:ro "${RO_V2[@]}" "$IMAGE" python "$@"
}

step() {  # name out_dir command...: run, log the exit and the head of the output
  local name=$1 out=$2 rc
  shift 2
  mkdir -p "$out"
  log "$name: start"
  "$@" > "$out/out.json" 2> "$out/err.txt"
  rc=$?
  log "$name: exit $rc: $(cut -c1-400 "$out/out.json") $(tail -n 2 "$out/err.txt")"
  [ "$rc" = 0 ] || fail "$name exit $rc"
}

log "c48c2: mode '$MODE' (pid $$), tree $TREE, image $IMAGE"
case "$MODE" in
  look-gv)  step look-gv "$C/work/look-gv" py "$BIN/c48c_v2_gv1_look.py" gvmags "$GVTAR" "$C/work/look-gv" ;;
  look-ino) step look-ino "$C/work/look-ino" py "$BIN/c48c_v2_gv1_look.py" inovirus "$METAZIP" "$C/work/look-ino" ;;
  look-vr)  step look-vr "$C/work/look-vr" py "$BIN/c48c_v2_gv1_look.py" imgvr "$HCZIP" "$C/work/imgvr" ;;
  inputs)   # stage V2's inputs into $V2 (written as $V2.tmp, renamed when complete); refused when $V2 exists
    [ -e "$V2" ] && fail "$V2 exists"
    step inputs "$C/work/inputs" py "$BIN/c48c_virus_inputs.py" v2 "$HCZIP" "$METAZIP" "$V2" --root "$ROOT" \
      --root-stages "$OWNER_STAGES" --owner-inputs V1=/tmp/if2453/s2/inputs/V1
    cp "$C/work/inputs/out.json" "$C/work/inputs/report.json" ;;
  rerun)    # the same command into scratch; every file compared with $V2's (cmp), then the copy removed
    rm -rf "$C/work/rerun"
    mkdir -p "$C/work/rerun"
    step rerun "$C/work/rerun/log" py "$BIN/c48c_virus_inputs.py" v2 "$HCZIP" "$METAZIP" "$C/work/rerun/V2" \
      --root "$ROOT" --root-stages "$OWNER_STAGES" --owner-inputs V1=/tmp/if2453/s2/inputs/V1
    for f in $V2_FILES; do
      if [ "$f" = objects.json ]; then
        log "rerun: objects.json names its own paths; compared after them: $(sed "s#$C/work/rerun/V2#$V2#" \
          "$C/work/rerun/V2/$f" | cmp - "$V2/$f" && echo equal)"
      else
        log "rerun: $f $(cmp "$C/work/rerun/V2/$f" "$V2/$f" && echo equal)"
      fi
    done
    rm -rf "$C/work/rerun/V2" ;;
  check)    # the inputs step's own parsing and checks (writes nothing), and the sha256 of every file of $V2
    step check "$C/work/check" py "$BIN/s2_check_inputs.py" V2 "$V2/genomes.json" "$V2/objects.json"
    (cd "$V2" && sha256sum $V2_FILES) > "$C/work/check/inputs_SUMS"
    log "check: SUMS: $(tr '\n' ' ' < "$C/work/check/inputs_SUMS")" ;;
  overlap)  # by contig_uid: V2 and GV1's candidates against every stored stage, V1's inputs and each other
    step overlap "$C/work/overlap" py "$BIN/c48c_v2_gv1_overlap.py" --root "$ROOT" --stages "$OWNER_STAGES" \
      --v1 /tmp/if2453/s2/inputs/V1 --v2 "$V2" --gvmags "$GVTAR" --gv-genomes "$C/work/look-gv/gvmags_genomes.tsv" \
      --out "$C/work/overlap" ;;
  ani)      # skani at the vOTU threshold, V2's sample against S1 and H1 (their stored shards) and V1's inputs
    step ani "$C/work/ani" py "$BIN/c48c_v2_gv1_ani.py" --v2 "$V2" \
      --s1 "$ROOT/l0seq/seq_shards/stage=S1/shard-0000.fa.gz" --h1 "$ROOT/l0seq/seq_shards/stage=H1/shard-0000.fa.gz" \
      --v1 /tmp/if2453/s2/inputs/V1/uhgv_v1.0_votus_full.fna --out "$C/work/ani" --threads "${THREADS:-3}" \
      --every "${EVERY:-20}" ${REUSE:+--reuse} ;;
  est-v2)
    step est-v2 "$C/work/est" py "$BIN/c48c_v2_gv1_est.py" v2 --root "$ROOT" --inputs "$V2" --out "$C/work/est" \
      --every "${EVERY:-40}" ;;
  est-gv1)
    step est-gv1 "$C/work/est-gv1" py "$BIN/c48c_v2_gv1_est.py" gv1 --root "$ROOT" --gvmags "$GVTAR" \
      --contigs "$C/work/overlap/gv1_candidate_contigs.tsv" --genomes "$C/work/look-gv/gvmags_genomes.tsv" \
      --held "$C/work/overlap/gv1_held_contigs.txt" --out "$C/work/est-gv1" --every "${EVERY:-50}" ;;
  est-full) # IMG/VR v4.0's whole table, metadata only (its 48 GB FASTA is not read)
    step est-full "$C/work/est-full" py "$BIN/c48c_v2_gv1_est.py" imgvr-full --full-zip "$JGI/jgi_img_vr_v4_full_20261118.zip" \
      --hc-reps "$V2" --out "$C/work/est-full" ;;
  left-v2)  # DEC-59: the HC rows V2 does not take, beside V2's inputs (not one of the inputs build-stage reads)
    step left-v2 "$C/work/left-v2" py "$BIN/c48c_v2_gv1_left.py" v2 \
      "$C/work/imgvr/IMGVR_all_Sequence_information-high_confidence.tsv" "$V2" "$V2/V2_left_for_later_hc.tsv.gz" ;;
  left-rerun) # the HC list again into scratch, compared byte for byte with the one beside V2's inputs
    mkdir -p "$C/work/left-rerun"
    step left-rerun "$C/work/left-rerun/log" py "$BIN/c48c_v2_gv1_left.py" v2 \
      "$C/work/imgvr/IMGVR_all_Sequence_information-high_confidence.tsv" "$V2" "$C/work/left-rerun/hc.tsv.gz"
    log "left-rerun: $(cmp "$C/work/left-rerun/hc.tsv.gz" "$V2/V2_left_for_later_hc.tsv.gz" && echo equal)"
    rm -f "$C/work/left-rerun/hc.tsv.gz" ;;
  left-put) # DEC-59: v4.0's putative rows (in no v4.1 high-confidence file); moved beside V2's inputs once they exist
    mkdir -p "$C/work/left"
    step left-put "$C/work/left-put" py "$BIN/c48c_v2_gv1_left.py" v4-putative \
      "$C/work/est-full/IMGVR_all_Sequence_information.tsv" "$C/work/left/V2_left_for_later_v4.0_putative.tsv.gz" ;;
  left-gv1) # DEC-59: GV1's candidates, none taken yet
    mkdir -p "$C/work/left"
    step left-gv1 "$C/work/left-gv1" py "$BIN/c48c_v2_gv1_left.py" gv1 "$C/work/look-gv/gvmags_genomes.tsv" \
      "$C/work/left/GV1_left_for_later.tsv.gz" ;;
  gv1-inputs) # stage GV1 into $GV1 (written as $GV1.tmp, renamed when complete); refused when $GV1 exists
    [ -e "$GV1" ] && fail "$GV1 exists"
    step gv1-inputs "$C/work/inputs" py "$BIN/c48c_virus_inputs.py" gv1 "$GVTAR" "$GV1" --root "$ROOT" \
      --root-stages "$OWNER_STAGES" --owner-inputs V1=/tmp/if2453/s2/inputs/V1 --owner-inputs V2="$V2"
    cp "$C/work/inputs/out.json" "$C/work/inputs/report.json" ;;
  gv1-rerun)  # the same command into scratch; every file compared with $GV1's (cmp), then the copy removed
    rm -rf "$C/work/rerun"
    mkdir -p "$C/work/rerun"
    step gv1-rerun "$C/work/rerun/log" py "$BIN/c48c_virus_inputs.py" gv1 "$GVTAR" "$C/work/rerun/GV1" \
      --root "$ROOT" --root-stages "$OWNER_STAGES" --owner-inputs V1=/tmp/if2453/s2/inputs/V1 --owner-inputs V2="$V2"
    for f in $GV1_FILES; do
      if [ "$f" = objects.json ]; then
        log "gv1-rerun: objects.json names its own path; compared after it: $(sed "s#$C/work/rerun/GV1#$GV1#" \
          "$C/work/rerun/GV1/$f" | cmp - "$GV1/$f" && echo equal)"
      else
        log "gv1-rerun: $f $(cmp "$C/work/rerun/GV1/$f" "$GV1/$f" && echo equal)"
      fi
    done
    rm -rf "$C/work/rerun/GV1" ;;
  gv1-check)  # the inputs step's own parsing and checks (writes nothing), and the sha256 of every file of $GV1
    step gv1-check "$C/work/check" py "$BIN/s2_check_inputs.py" GV1 "$GV1/genomes.json" "$GV1/objects.json"
    (cd "$GV1" && sha256sum $GV1_FILES) > "$C/work/check/inputs_SUMS"
    log "gv1-check: SUMS: $(tr '\n' ' ' < "$C/work/check/inputs_SUMS")" ;;
  gv1-overlap) # by contig_uid: GV1's staged contigs against the stored stages, C5, V1's and V2's inputs
    step gv1-overlap "$C/work/overlap" py "$BIN/c48c_v2_gv1_overlap.py" --root "$ROOT" --stages "$OWNER_STAGES" \
      --v1 /tmp/if2453/s2/inputs/V1 --v2 "$V2" --gvmags "$GVTAR" \
      --gv-genomes /tmp/if2453/c48c2/work/look-gv/gvmags_genomes.tsv --gv1-inputs "$GV1" --out "$C/work/overlap" ;;
  gv1-subset) # what fna_files.tar.gz is relative to GV1's genomes
    step gv1-subset "$C/work/subset" py "$BIN/c48c_v2_gv1_look.py" gvsubset "$GVTAR" "$GV1" "$C/work/subset" ;;
  gv1-est)    # the code's estimates on GV1's inputs, with a gene-calling sample (one thread)
    step gv1-est "$C/work/est" py "$BIN/c48c_v2_gv1_est.py" stage --root "$ROOT" --inputs "$GV1" --out "$C/work/est" \
      --every "${EVERY:-25}" ;;
  *) echo "usage: c48c_v2_gv1_run.sh look-gv|look-ino|look-vr|inputs|rerun|check|overlap|ani|est-v2|est-gv1|est-full|left-v2|left-put|left-gv1|gv1-inputs|gv1-rerun|gv1-check|gv1-overlap|gv1-subset|gv1-est" >&2
     exit 2 ;;
esac
done_
