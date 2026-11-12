#!/bin/bash
# Card 0-08b: one job lane (run ON THE POD HOST, detached, in its own capped unit):
#
#   systemd-run --user --unit=m08b-lane<k> --collect --quiet -p MemoryMax=16G -p MemorySwapMax=0 \
#     -p CPUQuota=200% -p StandardOutput=append:/tmp/if2453/m08b/logs/lane<k>.log \
#     -p StandardError=append:/tmp/if2453/m08b/logs/lane<k>.log nice -n 19 bash ncrna_set_lane.sh <k> [args]
#
# Each lane runs its jobs one after the other; at most two lanes run at once, so at most two processes of two
# threads (the unit's CPUQuota holds a lane to 2 CPUs). Before every job the account's CPU use (the sum of
# `ps -u if2453 -o pcpu=`, the ps line left out) plus 400% must be at most 4000% (the card's rule, at any hour);
# else the lane polls every 4 minutes and gives up after 3 hours of waiting for one job (exit 3). A STOP file in
# /tmp/if2453/m08b ends the lane before its next job. A job whose record runs/<name>.json has exit_code 0 is
# skipped, so a lane can be relaunched.
#
# Every job runs in the image with the data it reads bound read-only: the release root /tmp/if2453/r10/root
# (a 50 Gbp build is writing there; this card never writes there) and measurement 0-09c's /tmp/if2453/m09c.
set -uo pipefail
M=/tmp/if2453/m08b
IMG=/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif
TREE=/tmp/if2453/ci/m08b
BIN=$TREE/scripts/stress/annotation_cost
M9=/tmp/if2453/m09c
R10=/tmp/if2453/r10/root
CMGZ=/tmp/if2453/golden/rfam/Rfam.cm.gz
LANE=${1:?lane}
mkdir -p "$M/logs" "$M/runs" "$M/work" "$M/m2/fa" "$M/m2/out" "$M/m5" "$M/tmp" "$M/wrap"
cp "$BIN/cmsearch_wrap.sh" "$M/wrap/cmsearch"
chmod +x "$M/wrap/cmsearch"

log() { echo "lane$LANE: $(date -Is) $*"; }

img() {  # a command inside the image; M5 jobs see the cmsearch wrapper first on PATH
  singularity exec --cleanenv --bind /tmp/if2453 --bind "$R10:$R10:ro" --bind "$M9:$M9:ro" \
    --env PYTHONPATH="$TREE/src:/tmp/if2453/pilot_s1/pydeps" --env OMP_NUM_THREADS=2 \
    --env PATH="$M/wrap:/opt/bio-box/.venv/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
    --env M08B_REC="${M08B_REC:-$M/m5/calls}" --env M08B_LABEL="${M08B_LABEL:-x}" --env M08B_TIMED="$BIN/timed.py" \
    --pwd "$M/work" "$IMG" "$@"
}

gate() {
  local t0
  t0=$(date +%s)
  while :; do
    if [ -e "$M/STOP" ]; then log "STOP file: exiting"; exit 0; fi
    used=$(ps -u if2453 -o pcpu=,comm= | awk '$2 != "ps" {s+=$1} END {printf "%d", s}')
    if [ $((used + 400)) -le 4000 ]; then
      log "gate open for $1: account ${used}% + 400% <= 4000%"
      return 0
    fi
    if [ $(( $(date +%s) - t0 )) -ge 10800 ]; then log "gave up on $1 after 3 hours of waiting"; exit 3; fi
    log "gate waiting for $1: account ${used}% + 400% > 4000%"
    sleep 240
  done
}

done_ok() { [ -f "$M/runs/$1.json" ] && grep -q '"exit_code": 0' "$M/runs/$1.json"; }

job() {  # job <name> <command...>: gated, timed, inside the image
  local name=$1
  shift
  if done_ok "$name"; then log "skip $name (done)"; return 0; fi
  gate "$name"
  log "start $name"
  img python "$BIN/timed.py" "$M/runs/$name.json" -- "$@" >> "$M/logs/$name.log" 2>&1
  local rc=$?
  log "end $name exit $rc"
  return $rc
}

wait_for() {  # wait_for <file>: poll every minute, at most 8 hours
  local i=0
  until [ -e "$1" ]; do
    if [ -e "$M/STOP" ]; then log "STOP file: exiting"; exit 0; fi
    i=$((i + 1))
    if [ "$i" -gt 480 ]; then log "gave up waiting for $1"; exit 3; fi
    sleep 60
  done
}

trna_scan() {  # trna_scan <name> <records stem> [extra cmsearch flags]: the producer's argv, the CM $TRNA_CM
  local name=$1 stem=$2
  shift 2
  job "$name" /opt/conda/bin/cmsearch --cpu 2 -Z 100000 --cut_ga --noali --tformat fasta --tblout "$M/m2/out/$name.tbl" \
    -o "$M/m2/out/$name.out" "$@" "${TRNA_CM:-$M/m2/trna.cm}" "$stem.fna" || return $?
  img python "$BIN/ncrna_set.py" collect "$M/m2/out/$name.tbl" "$stem.map.json" "$M9/fa/wg_all.map.json" \
    "$M/m2/out/$name.hits.jsonl" > "$M/m2/out/$name.collect.json" 2>> "$M/logs/$name.log"
  img python "$BIN/infernal_igr.py" permodel "$M/m2/out/$name.out" "$M/m2/out/$name.permodel.json" \
    > "$M/m2/out/$name.permodel.sum.json" 2>> "$M/logs/$name.log"
  log "collected $name: $(cat "$M/m2/out/$name.collect.json")"
}

m5part() {  # m5part <k>: the producer on part k of the sample, into its own scratch root
  local k=$1
  wait_for "$M/m5/sample.json"
  export M08B_LABEL="part$k" M08B_REC="$M/m5/calls"
  job "m5_part$k" python "$BIN/ncrna_set.py" m5run "$R10" S1 "$M/m5/sample.json" "$k" "$CMGZ" 2 \
    "$M/m5/out$k" "$M/tmp/m5_part$k" "$M/m5/result$k.json"
}

case "$LANE" in
1)
  # M2, phase 1: the tRNA CM, the default (filtered) scan of the 8 genomes whole, and a --max calibration on the
  # two phages lambda and T4; then M5 part 1
  if [ ! -f "$M/m2/trna.cm" ]; then
    printf 'RF00005\nRF01852\n' > "$M/m2/trna.keys"
    img cmfetch -f -o "$M/m2/trna.cm.tmp" "$M9/lib/Rfam.cm" "$M/m2/trna.keys" >> "$M/logs/m2_cm.log" 2>&1 || exit 1
    mv "$M/m2/trna.cm.tmp" "$M/m2/trna.cm"
  fi
  [ -f "$M/m2/fa/wg_whole.map.json" ] || img python "$BIN/ncrna_set.py" records "$M9/fa/wg_all.fna" \
    "$M9/fa/wg_all.map.json" all 1000000000 0 "$M/m2/fa/wg_whole.fna" "$M/m2/fa/wg_whole.map.json" \
    >> "$M/logs/m2_records.log" 2>&1 || exit 1
  [ -f "$M/m2/fa/cal.map.json" ] || img python "$BIN/ncrna_set.py" records "$M9/fa/wg_all.fna" \
    "$M9/fa/wg_all.map.json" lambda,t4 1000000000 0 "$M/m2/fa/cal.fna" "$M/m2/fa/cal.map.json" \
    >> "$M/logs/m2_records.log" 2>&1 || exit 1
  trna_scan trna_default "$M/m2/fa/wg_whole"
  trna_scan trna_max_cal "$M/m2/fa/cal" --max
  m5part 1
  ;;
2)
  # M5: the sample, its regions, part 0
  job m5_sample python "$BIN/ncrna_set.py" sample "$R10" S1 50000000 "$M/m5/sample.json" || exit 1
  job m5_regions python "$BIN/ncrna_set.py" m5regions "$R10" S1 "$M/m5/sample.json" "$M9/fa/igr_m200.map.json" \
    "$M/m5/regions.json"
  m5part 0
  ;;
3)
  # M2, phase 2: --max with the tRNA model RF00005 alone (the calibration put both models over the 8 genomes at
  # 2.3 times the design note's estimate; RF00005 alone fits it), on the genomes $4 (a comma list, or all) in
  # chunks of $2 bases overlapping by $3 (a chunk larger than a record leaves it whole); name suffix $5
  CH=${2:?chunk} OV=${3:?overlap} GS=${4:?genomes} SUF=${5:?suffix}
  if [ ! -f "$M/m2/trna1.cm" ]; then
    printf 'RF00005\n' > "$M/m2/trna1.keys"
    img cmfetch -f -o "$M/m2/trna1.cm.tmp.$SUF" "$M9/lib/Rfam.cm" "$M/m2/trna1.keys" >> "$M/logs/m2_cm.log" 2>&1 || exit 1
    mv "$M/m2/trna1.cm.tmp.$SUF" "$M/m2/trna1.cm"
  fi
  [ -f "$M/m2/fa/max_$SUF.map.json" ] || img python "$BIN/ncrna_set.py" records "$M9/fa/wg_all.fna" \
    "$M9/fa/wg_all.map.json" "$GS" "$CH" "$OV" "$M/m2/fa/max_$SUF.fna" "$M/m2/fa/max_$SUF.map.json" \
    >> "$M/logs/m2_records.log" 2>&1 || exit 1
  export TRNA_CM="$M/m2/trna1.cm"
  trna_scan "trna_max_$SUF" "$M/m2/fa/max_$SUF" --max
  ;;
*)
  echo "lane 1 to 3" >&2
  exit 2
  ;;
esac
log "lane done"
touch "$M/lane$LANE.done"
