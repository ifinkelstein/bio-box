#!/bin/sh
# Card C-41: build stage G1 inputs in /tmp/if2453/s2/inputs/G1/, dry-check them, verify determinism, and run c41_est.py.
# Run inside the pinned image under a memory- and CPU-capped systemd unit (at most 2 threads, root bound read-only):
#
#   systemd-run --user --unit=c41-run-G1 --collect --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% \
#     -p StandardOutput=append:/tmp/if2453/s3/globdb/G1_out/run.log \
#     -p StandardError=append:/tmp/if2453/s3/globdb/G1_out/run.log \
#     nice -n 19 singularity exec --cleanenv --env PYTHONUNBUFFERED=1 --env PYTHONPATH=/tmp/if2453/ci/C-41/src \
#     --env OMP_NUM_THREADS=2 --bind /tmp/if2453/ci/C-41:/tmp/if2453/ci/C-41:ro --bind /tmp/if2453/s3/globdb \
#     --bind /tmp/if2453/s2/inputs --bind /tmp/if2453/r10/root:/tmp/if2453/r10/root:ro \
#     --bind /sys/fs/cgroup:/sys/fs/cgroup:ro /stor/home/if2453/bio-box/images/bio-box_faec1f4.sif \
#     sh /tmp/if2453/ci/C-41/scripts/stress/build_stage_pod/c41_run.sh all
set -eu

TREE=/tmp/if2453/ci/C-41
SCRIPTS=$TREE/scripts/stress/build_stage_pod
W=/tmp/if2453/s3/globdb
GEN=$W/genomes/G1
OUT=$W/G1_out
IN=/tmp/if2453/s2/inputs/G1
ROOT=/tmp/if2453/r10/root

now() { date -u +%FT%TZ; }
check_stop() {
  for f in "$ROOT/STOP" "$W/STOP"; do
    if [ -e "$f" ]; then
      echo "STOP present: $f $(now)" >&2
      exit 75
    fi
  done
}

mkdir -p "$OUT"

do_inputs() {
  check_stop
  echo "[$(now)] start make_inputs.py prok -> $IN"
  rm -rf "$IN"
  python "$SCRIPTS/make_inputs.py" prok \
    "$GEN/lock.json" "$GEN/selection.tsv" "$IN" \
    --source-name globdb --source-release r226 \
    --license "CC BY-SA 4.0 (GlobDB; underlying MAG catalogues)" \
    --quality-column assembly_quality \
    --root "$ROOT" > "$OUT/report.json"
  cat "$OUT/report.json"

  check_stop
  echo "[$(now)] start s2_check_inputs.py G1"
  python "$SCRIPTS/s2_check_inputs.py" G1 "$IN/genomes.json" "$IN/objects.json" > "$OUT/check_full.json"
  python -c '
import json, sys
from pathlib import Path
d = json.loads(Path(sys.argv[1]).read_text())
objs = d.pop("objects")
d["objects_summary"] = {
    "count": len(objs),
    "total_bytes": sum(o["bytes"] for o in objs),
    "source_names": sorted({o["source_name"] for o in objs}),
    "source_releases": sorted({o["source_release"] for o in objs}),
    "licenses": sorted({o["license"] for o in objs}),
    "redistributable": sorted({o["redistributable"] for o in objs}),
    "first": objs[0] if objs else None,
    "last": objs[-1] if objs else None,
}
Path(sys.argv[2]).write_text(json.dumps(d, indent=1) + "\n")
print("problems:", d["problems"])
' "$OUT/check_full.json" "$OUT/check.json"

  check_stop
  echo "[$(now)] determinism check (second run into scratch)"
  RERUN=$W/G1_rerun
  rm -rf "$RERUN"
  python "$SCRIPTS/make_inputs.py" prok \
    "$GEN/lock.json" "$GEN/selection.tsv" "$RERUN" \
    --source-name globdb --source-release r226 \
    --license "CC BY-SA 4.0 (GlobDB; underlying MAG catalogues)" \
    --quality-column assembly_quality \
    --root "$ROOT" > "$OUT/report_rerun.json"
  cmp -s "$IN/genomes.json" "$RERUN/genomes.json"
  cmp -s "$IN/objects.json" "$RERUN/objects.json"
  cmp -s "$IN/shared_records.tsv" "$RERUN/shared_records.tsv"
  cmp -s "$OUT/report.json" "$OUT/report_rerun.json"
  rm -rf "$RERUN" "$OUT/report_rerun.json"
  echo "[$(now)] determinism check: byte-identical"

  sha256sum \
    "$IN/genomes.json" \
    "$IN/objects.json" \
    "$IN/shared_records.tsv" \
    "$GEN/lock.json" \
    "$GEN/selection.tsv" \
    "$GEN/stream_summary.json" > "$OUT/inputs_SUMS"
  cat "$OUT/inputs_SUMS"
}

do_est() {
  check_stop
  echo "[$(now)] start c41_est.py"
  python "$SCRIPTS/c41_est.py" \
    --root "$ROOT" \
    --inputs "$IN" \
    --lock "$GEN/lock.json" \
    --selection "$GEN/selection.tsv" \
    --table-summary "$W/tables/table_summary.json" \
    --profile "$TREE/profiles/wilkcomp02.yaml" \
    --pl2-inputs /tmp/if2453/s2/inputs/PL2 \
    --v1-inputs /tmp/if2453/s2/inputs/V1 \
    --every 20 \
    --crispr-every 100 \
    --out "$OUT"
  echo "[$(now)] done c41_est.py"
}

case "${1:-all}" in
  inputs) do_inputs ;;
  est) do_est ;;
  all) do_inputs; do_est ;;
  *) echo "unknown mode $1" >&2; exit 2 ;;
esac
echo "[$(now)] c41_run.sh finished ($1)"
