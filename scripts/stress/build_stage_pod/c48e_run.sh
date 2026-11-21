#!/bin/sh
# Card C-48e part 1: the tranche's pod jobs, one mode per detached unit, each leaving <W>/<mode>.DONE or .FAILED.
# Run inside the pinned image under a capped unit at nice 19 (the data root bound read-only), for example:
#
#   systemd-run --user --unit=c48e1-fetch --collect --quiet -p MemoryMax=32G -p MemorySwapMax=0 -p CPUQuota=100% \
#     -p StandardOutput=append:/tmp/if2453/c48e1/logs/fetch.log -p StandardError=append:/tmp/if2453/c48e1/logs/fetch.log \
#     nice -n 19 singularity exec --cleanenv --env PYTHONPATH=/tmp/if2453/c48e1/tree/src --env OMP_NUM_THREADS=1 \
#     --bind /tmp/if2453 --bind /tmp/if2453/r10/root:/tmp/if2453/r10/root:ro \
#     /stor/home/if2453/bio-box/images/bio-box_faec1f4.sif sh /tmp/if2453/c48e1/tree/scripts/stress/build_stage_pod/c48e_run.sh fetch
#
# Modes: fetch (the assemblies, resumable), inputs (stage SPX1's inputs, its dry check and sums), rerun (the inputs
# again into scratch, compared byte for byte, then removed), est (overlap, measured sizes, costs, remainder).
set -u
W=${C48E_WORK:-/tmp/if2453/c48e1}
TREE=${C48E_TREE:-$W/tree}
S=$TREE/scripts/stress/build_stage_pod
ROOT=/tmp/if2453/r10/root
STAGE=${C48E_STAGE:-SPX1}
IN=/tmp/if2453/s2/inputs/$STAGE
now() { date -u +%FT%TZ; }

inner() {
  set -eu
  for f in "$ROOT/STOP" "$W/STOP"; do
    if [ -e "$f" ]; then echo "[$(now)] STOP present: $f" >&2; exit 75; fi
  done
  case "$1" in
    fetch)
      python "$S/c48e_fetch.py" --selection "$W/select/selection.tsv" --out-dir "$W/dl" --interval 1.0 --tries 3 \
        --max-dl-gb 150 --min-free-gb 400 --stop-file "$ROOT/STOP" --stop-file "$W/STOP"
      sha256sum "$W/dl/fetch_log.jsonl" > "$W/dl/fetch_log.sha256"
      cat "$W/dl/fetch_log.sha256"
      ;;
    inputs)
      python "$S/c48e_inputs.py" build --selection "$W/select/selection.tsv" --fetch-log "$W/dl/fetch_log.jsonl" \
        --dl-dir "$W/dl" --stage "$STAGE" --final-dir "$IN" --out-dir "$IN" --min-contig-len 5000 \
        --uids-out "$W/work/contigs_ge5000.parquet" --tmp "$W/tmp" > "$W/out/inputs_report.json"
      cat "$W/out/inputs_report.json"
      python "$S/s2_check_inputs.py" "$STAGE" "$IN/genomes.json" "$IN/objects.json" > "$W/out/check_full.json"
      python "$S/c48e_inputs.py" check-summary "$W/out/check_full.json" "$W/out/check.json"
      (cd "$IN" && sha256sum genomes.json objects.json samples.tsv input_report.json assemblies.sha256) > "$W/out/inputs_SUMS"
      cat "$W/out/inputs_SUMS"
      ;;
    rerun)
      rm -rf "$W/rerun"
      mkdir -p "$W/rerun"
      python "$S/c48e_inputs.py" build --selection "$W/select/selection.tsv" --fetch-log "$W/dl/fetch_log.jsonl" \
        --dl-dir "$W/dl" --stage "$STAGE" --final-dir "$IN" --out-dir "$W/rerun/$STAGE" --min-contig-len 5000 \
        --uids-out "$W/rerun/contigs_ge5000.parquet" --tmp "$W/tmp" > "$W/rerun/inputs_report.json"
      for f in genomes.json objects.json samples.tsv input_report.json assemblies.sha256; do
        cmp "$IN/$f" "$W/rerun/$STAGE/$f"
        echo "identical: $f"
      done
      cmp "$W/work/contigs_ge5000.parquet" "$W/rerun/contigs_ge5000.parquet"
      echo "identical: contigs_ge5000.parquet"
      cmp "$W/out/inputs_report.json" "$W/rerun/inputs_report.json"
      echo "identical: inputs_report.json (stdout)"
      n=0
      for f in "$IN"/assemblies/*.fa.gz; do
        b=$(basename "$f")
        cmp "$f" "$W/rerun/$STAGE/assemblies/$b"
        n=$((n + 1))
      done
      echo "identical: $n assemblies"
      rm -rf "$W/rerun"
      ;;
    est)
      python "$S/c48e_est.py" --root "$ROOT" --inputs "$IN" --uids "$W/work/contigs_ge5000.parquet" \
        --selection "$W/select/selection.tsv" --select-dir "$W/select" --fetch-log "$W/dl/fetch_log.jsonl" \
        --pilot-depth1 /tmp/if2453/s3/spire/depth1/depth1_by_stage.tsv \
        --pilot-depth1-unfiltered /tmp/if2453/s3/spire/depth1_unfiltered/depth1_by_stage.tsv \
        --pilot-pack /tmp/if2453/s3/spire/depth1/pack_probe_SPA1.json \
        --profile "$TREE/profiles/wilkcomp02.yaml" --other-inputs /tmp/if2453/s2/inputs --tmp "$W/tmp" --out "$W/out"
      ;;
    *) echo "unknown mode $1" >&2; exit 2 ;;
  esac
}

if [ "${1:-}" = "_inner" ]; then
  inner "$2"
  exit $?
fi
MODE=${1:?mode}
mkdir -p "$W/logs" "$W/out" "$W/work" "$W/tmp"
rm -f "$W/$MODE.DONE" "$W/$MODE.FAILED"
echo "[$(now)] start $MODE"
sh "$0" _inner "$MODE"
rc=$?
if [ "$rc" -eq 0 ]; then
  now > "$W/$MODE.DONE"
else
  echo "rc=$rc $(now)" > "$W/$MODE.FAILED"
fi
echo "[$(now)] end $MODE rc=$rc"
exit "$rc"
