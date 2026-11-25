#!/bin/bash
# Card C-38c pod driver: runs linclust alone on r26.11.1's 49.1M proteins, compares against
# r26.11.1's published cascade in bounded DuckDB, builds a scratch release view, and compares
# TniQ & RNA_pol `nbr summarize` outputs under cascade vs linclust.
#
# Run ON wilkcomp02.ccbb.utexas.edu (detached via systemd-run --user --unit=c38c-driver --collect):
#   TREE=/tmp/if2453/ci/C-38c bash "$TREE/scripts/stress/c38c_pod_driver.sh"
set -euo pipefail

TREE=${TREE:-/tmp/if2453/ci/C-38c}
ROOT=${ROOT:-/tmp/if2453/r10/root}
PACKS=${PACKS:-/tmp/if2453/bio-box/packs}
PFAM_DIR=${PFAM_DIR:-/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2}
IMAGE=${IMAGE:-/stor/home/if2453/bio-box/images/bio-box-test_faec1f4.sif}
RELEASE=${RELEASE:-r26.11.1}
WORK=${WORK:-/tmp/if2453/c38c}
THREADS=${THREADS:-12}
MEM=${MEM:-64G}
ALIGN_PARTS=${ALIGN_PARTS:-4}
PROFILE="$TREE/profiles/wilkcomp02.yaml"
SCRIPT="$TREE/scripts/stress/c38c_linclust_vs_cascade.py"

[ ! -e "$ROOT/STOP" ] || { echo "STOP file present in $ROOT" >&2; exit 2; }
mkdir -p "$WORK/logs" "$WORK/work" "$WORK/summary" "$WORK/hmm" "$WORK/runs" "$WORK/search_work" "$WORK/tmp"
METRICS="$WORK/metrics.tsv"
if [ ! -f "$METRICS" ]; then
  printf 'step\twall_s\tcpu_nsec\tcpu_min\tmemory_peak_bytes\texit_code\n' > "$METRICS"
fi

IMAGE_DIGEST="sha256:$(sha256sum "$IMAGE" | awk '{print $1}')"

# Run one command in a capped systemd user unit, capturing cgroup CPUUsageNSec and MemoryPeak.
run_unit() {
  local name=$1
  shift
  local unit="c38c-${name}"
  local ulog="$WORK/logs/unit-${name}.log"
  local t0 t1 wall cpu_ns mem_peak rc
  t0=$(date +%s)
  set +e
  systemd-run --user --unit="$unit" --wait --collect --quiet \
    -p MemoryMax="$MEM" -p MemorySwapMax=0 -p CPUQuota="$((THREADS * 100))%" -p OOMPolicy=continue \
    -p StandardOutput="truncate:$ulog" -p StandardError="append:$ulog" \
    /bin/sh -c 'unit=$1; shift; "$@"; rc=$?; use=$(systemctl --user show "$unit" -p CPUUsageNSec -p MemoryPeak 2>/dev/null | paste -sd" " -); echo "C38C_USAGE $use"; echo "C38C_EXIT $rc"; exit $rc' \
    sh "$unit" \
    nice -n 19 singularity exec --cleanenv \
      --env PYTHONUNBUFFERED=1 \
      --env "PYTHONPATH=$TREE/src" \
      --env "BIOBOX_IMAGE_DIGEST=$IMAGE_DIGEST" \
      --env "BIOBOX_POD_UNIT=$unit" \
      --bind /sys/fs/cgroup:/sys/fs/cgroup:ro \
      --bind "$TREE:$TREE:ro" \
      --bind "$ROOT:$ROOT:ro" \
      --bind "$PACKS:$PACKS:ro" \
      --bind "$PFAM_DIR:$PFAM_DIR:ro" \
      --bind "$WORK:$WORK" \
      "$IMAGE" "$@"
  rc=$?
  set -e
  t1=$(date +%s)
  wall=$((t1 - t0))
  cpu_ns=$(sed -n 's/^C38C_USAGE .*CPUUsageNSec=\([0-9]*\).*/\1/p' "$ulog" | tail -n 1)
  mem_peak=$(sed -n 's/^C38C_USAGE .*MemoryPeak=\([0-9]*\).*/\1/p' "$ulog" | tail -n 1)
  cpu_ns=${cpu_ns:-0}
  mem_peak=${mem_peak:-0}
  local cpu_min
  cpu_min=$(awk -v ns="$cpu_ns" 'BEGIN { printf "%.4f", ns / 60000000000.0 }')
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$wall" "$cpu_ns" "$cpu_min" "$mem_peak" "$rc" >> "$METRICS"
  echo "[$(date -Is)] step=$name rc=$rc wall=${wall}s cpu=${cpu_min}m mem_gb=$(awk -v b="$mem_peak" 'BEGIN{printf "%.2f", b/1e9}')"
  if [ "$rc" -ne 0 ]; then
    echo "FAILED step=$name rc=$rc; tail of $ulog:" >&2
    tail -n 40 "$ulog" >&2
    exit "$rc"
  fi
}

log_du() {
  local label=$1
  local bytes
  bytes=$(du -sb "$WORK" | awk '{print $1}')
  echo "[$(date -Is)] scratch_du label=$label bytes=$bytes gib=$(awk -v b="$bytes" 'BEGIN{printf "%.2f", b/1073741824}')"
}

echo "[$(date -Is)] C-38c starting: release=$RELEASE threads=$THREADS mem=$MEM align_parts=$ALIGN_PARTS"

# 1. Prepare FASTA + createdb (FASTA removed immediately inside prepare)
run_unit "prepare" python "$SCRIPT" prepare --root "$ROOT" --release "$RELEASE" --work "$WORK/work"
log_du "after_prepare"

# 2. NR90 linclust (steps 0..9; step 7 align split into ALIGN_PARTS parts)
for step in \
  linclust.nr90.00.kmermatcher \
  linclust.nr90.01.rescorediagonal \
  linclust.nr90.02.clust \
  linclust.nr90.03.createsubdb \
  linclust.nr90.04.createsubdb \
  linclust.nr90.05.filterdb \
  linclust.nr90.06.rescorediagonal; do
  run_unit "${step//./-}" python "$SCRIPT" run-call --work "$WORK/work" --step "$step" --threads "$THREADS"
done

run_unit "nr90-07-align-cut" python "$SCRIPT" run-split-call --work "$WORK/work" --step linclust.nr90.07.align --action cut --parts "$ALIGN_PARTS"
for p in $(seq 0 $((ALIGN_PARTS - 1))); do
  run_unit "nr90-07-align-p${p}" python "$SCRIPT" run-split-call --work "$WORK/work" --step linclust.nr90.07.align --action part --part "$p" --threads "$THREADS"
done
run_unit "nr90-07-align-join" python "$SCRIPT" run-split-call --work "$WORK/work" --step linclust.nr90.07.align --action join

for step in linclust.nr90.08.clust linclust.nr90.09.mergeclusters; do
  run_unit "${step//./-}" python "$SCRIPT" run-call --work "$WORK/work" --step "$step" --threads "$THREADS"
done
log_du "after_nr90_before_clean"
rm -rf "$WORK/work/tmp90"
log_du "after_nr90_clean"

# 3. rep90 + NR50 linclust (steps 0..9; step 7 align split into ALIGN_PARTS parts) + merge
run_unit "cluster-rep90" python "$SCRIPT" run-call --work "$WORK/work" --step cluster.rep90 --threads "$THREADS"

for step in \
  linclust.nr50.00.kmermatcher \
  linclust.nr50.01.rescorediagonal \
  linclust.nr50.02.clust \
  linclust.nr50.03.createsubdb \
  linclust.nr50.04.createsubdb \
  linclust.nr50.05.filterdb \
  linclust.nr50.06.rescorediagonal; do
  run_unit "${step//./-}" python "$SCRIPT" run-call --work "$WORK/work" --step "$step" --threads "$THREADS"
done

run_unit "nr50-07-align-cut" python "$SCRIPT" run-split-call --work "$WORK/work" --step linclust.nr50.07.align --action cut --parts "$ALIGN_PARTS"
for p in $(seq 0 $((ALIGN_PARTS - 1))); do
  run_unit "nr50-07-align-p${p}" python "$SCRIPT" run-split-call --work "$WORK/work" --step linclust.nr50.07.align --action part --part "$p" --threads "$THREADS"
done
run_unit "nr50-07-align-join" python "$SCRIPT" run-split-call --work "$WORK/work" --step linclust.nr50.07.align --action join

for step in linclust.nr50.08.clust linclust.nr50.09.mergeclusters cluster.merge; do
  run_unit "${step//./-}" python "$SCRIPT" run-call --work "$WORK/work" --step "$step" --threads "$THREADS"
done
log_du "after_nr50_before_clean"
rm -rf "$WORK/work/tmp50" "$WORK/work/rep90"* "$WORK/work/clu50r"*

# 4. createtsv for clu90 and clu50, then remove raw MMseqs DBs to keep scratch lean
run_unit "createtsv" python "$SCRIPT" createtsv --work "$WORK/work" --threads "$THREADS"
mv "$WORK/work/clu90.tsv" "$WORK/clu90.tsv"
mv "$WORK/work/clu50.tsv" "$WORK/clu50.tsv"
rm -rf "$WORK/work/db"* "$WORK/work/clu90"* "$WORK/work/clu50"*
log_du "after_createtsv_clean"

# 5. Bounded DuckDB comparison against r26.11.1 cascade
run_unit "compare" python "$SCRIPT" compare \
  --root "$ROOT" --release "$RELEASE" \
  --clu90-tsv "$WORK/clu90.tsv" --clu50-tsv "$WORK/clu50.tsv" \
  --memory-limit 8GB --threads 4 --tmp "$WORK/tmp" \
  --out "$WORK/compare.json"

# 6. Build scratch release view with linclust cluster_member & member_cluster
run_unit "build-view" python "$SCRIPT" build-view \
  --root "$ROOT" --release "$RELEASE" \
  --clu90-tsv "$WORK/clu90.tsv" --clu50-tsv "$WORK/clu50.tsv" \
  --view "$WORK/view" --memory-limit 8GB --tmp "$WORK/tmp"
log_du "after_build_view"

# 7. Worked examples: TniQ and RNA_pol searches on r26.11.1 (read-only), then nbr summarize on cascade vs linclust
PFAM="$PFAM_DIR/Pfam-A.hmm.gz"
run_unit "hmmfetch" /bin/sh -c "hmmfetch -o '$WORK/hmm/TniQ.hmm' '$PFAM' TniQ && hmmfetch -o '$WORK/hmm/RNA_pol.hmm' '$PFAM' RNA_pol"
GA_TNIQ=$(sed -n 's/^GA  *\([0-9.]*\) .*/\1/p' "$WORK/hmm/TniQ.hmm")
echo "[$(date -Is)] GA_TNIQ=$GA_TNIQ"

run_unit "search-tniq" python -m biobox.cli search hmm \
  --root "$ROOT" --release "$RELEASE" --profile "$PROFILE" \
  --runs "$WORK/runs" --work "$WORK/search_work" --neighbors 10 \
  --cpus "$THREADS" --mem-gb 14 --threads "$THREADS" \
  --hmm "$WORK/hmm/TniQ.hmm" --bits "$GA_TNIQ" \
  --positive NZ_AP026781.1:426340-427525 --out "$WORK/tniq_out"

run_unit "search-rnap" python -m biobox.cli search hmm \
  --root "$ROOT" --release "$RELEASE" --profile "$PROFILE" \
  --runs "$WORK/runs" --work "$WORK/search_work" --neighbors 5 \
  --cpus "$THREADS" --mem-gb 14 --threads "$THREADS" \
  --hmm "$WORK/hmm/RNA_pol.hmm" --bits 30 \
  --positive NC_001604:3170-5822 --out "$WORK/rnap_out"

RUN_TNIQ=$(sed -n 's/^\([0-9a-f]\{32\}\)\tok$/\1/p' "$WORK/logs/unit-search-tniq.log" | head -n 1)
RUN_RNAP=$(sed -n 's/^\([0-9a-f]\{32\}\)\tok$/\1/p' "$WORK/logs/unit-search-rnap.log" | head -n 1)
echo "[$(date -Is)] RUN_TNIQ=$RUN_TNIQ RUN_RNAP=$RUN_RNAP"

# Cascade vs Linclust nbr summarize (TniQ: 10g & 20kb; RNA_pol: 10g & 5g)
for mode in cascade linclust; do
  if [ "$mode" = "cascade" ]; then
    R="$ROOT"
  else
    R="$WORK/view"
  fi
  run_unit "nbr-tniq-10g-${mode}" python -m biobox.cli nbr summarize "$RUN_TNIQ" \
    --root "$R" --profile "$PROFILE" --runs "$WORK/runs" \
    --by pfam,nr50 --top 0 --format tsv --window 10g --out "$WORK/summary/tniq_10g_${mode}.tsv"
  run_unit "nbr-tniq-20kb-${mode}" python -m biobox.cli nbr summarize "$RUN_TNIQ" \
    --root "$R" --profile "$PROFILE" --runs "$WORK/runs" \
    --by pfam,nr50 --top 0 --format tsv --window 20kb --out "$WORK/summary/tniq_20kb_${mode}.tsv"
  run_unit "nbr-rnap-10g-${mode}" python -m biobox.cli nbr summarize "$RUN_RNAP" \
    --root "$R" --profile "$PROFILE" --runs "$WORK/runs" \
    --by pfam,nr50 --top 0 --format tsv --window 10g --out "$WORK/summary/rnap_10g_${mode}.tsv"
  run_unit "nbr-rnap-5g-${mode}" python -m biobox.cli nbr summarize "$RUN_RNAP" \
    --root "$R" --profile "$PROFILE" --runs "$WORK/runs" \
    --by pfam,nr50 --top 0 --format tsv --window 5g --out "$WORK/summary/rnap_5g_${mode}.tsv"
done

log_du "final_before_cleanup"
echo "[$(date -Is)] C-38c ALL DONE"
