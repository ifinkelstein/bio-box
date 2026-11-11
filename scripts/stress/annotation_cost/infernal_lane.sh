#!/bin/bash
# Card 0-09c: one of the two job lanes (run ON THE POD HOST, detached, in its own capped unit):
#
#   systemd-run --user --unit=m09c-lane<k> --collect --quiet -p MemoryMax=16G -p MemorySwapMax=0 \
#     -p CPUQuota=200% -p StandardOutput=append:/tmp/if2453/m09c/logs/lane<k>.log \
#     -p StandardError=append:/tmp/if2453/m09c/logs/lane<k>.log nice -n 19 bash infernal_lane.sh <k>
#
# Each lane runs its jobs one after the other, so at most two cmsearch processes run at once, each with
# --cpu 2, and the unit's CPUQuota holds each lane to 2 CPUs (4 in all). Before every job the account's CPU
# use (the sum of `ps -u if2453 -o pcpu=`, the ps line left out) plus 200% must be at most 4000% between
# 08:00 and 20:00 CST on a weekday, else 5600%; otherwise it polls every 4 minutes (AGENTS §7.2, the
# card's limits). A STOP file in /tmp/if2453/m09c ends the lane before its next job. A job whose
# record runs/<name>.json has exit_code 0 is skipped, so a lane can be relaunched.
#
# Every cmsearch is the producer's argv (biobox.l1.ncrna._cmsearch_argv) with -o <file> for the per-model
# CPU lines; a variant appends its flags. Each runs under timed.py (wall, CPU-minutes, peak RSS).
set -uo pipefail
M=/tmp/if2453/m09c
IMG=/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif
TREE=$M/tree
BIN=$TREE/scripts/stress/annotation_cost
CM=$M/lib/Rfam.cm
G04=$M/work/g04
LANE=${1:?lane 1 or 2}
mkdir -p "$M/logs" "$M/runs" "$M/out" "$M/work" "$M/fa" "$M/scratch" "$M/lib"

log() { echo "lane$LANE: $(date -Is) $*"; }

img() {  # a command inside the image
  singularity exec --cleanenv --bind /tmp/if2453 \
    --env PYTHONPATH="$TREE/src:/tmp/if2453/pilot_s1/pydeps" --env OMP_NUM_THREADS=2 --pwd "$M/work" "$IMG" "$@"
}

gate() {
  while :; do
    if [ -e "$M/STOP" ]; then log "STOP file: exiting"; exit 0; fi
    used=$(ps -u if2453 -o pcpu=,comm= | awk '$2 != "ps" {s+=$1} END {printf "%d", s}')
    h=$((10#$(TZ=America/Chicago date +%H)))
    dow=$(TZ=America/Chicago date +%u)
    if [ "$dow" -le 5 ] && [ "$h" -ge 8 ] && [ "$h" -lt 20 ]; then cap=4000; else cap=5600; fi
    if [ $((used + 200)) -le "$cap" ]; then
      log "gate open for $1: account ${used}% + 200% <= ${cap}%"
      return 0
    fi
    log "gate waiting for $1: account ${used}% + 200% > ${cap}%"
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

scan() {  # scan <name> <cm> <fasta stem> [extra cmsearch flags...]
  local name=$1 cm=$2 stem=$3
  shift 3
  job "$name" cmsearch --cpu 2 -Z 100000 --cut_ga --noali --tformat fasta --tblout "$M/out/$name.tbl" \
    -o "$M/out/$name.out" "$@" "$cm" "$stem.fna" || return $?
  img python "$BIN/infernal_igr.py" collect "$M/out/$name.tbl" "$stem.map.json" "$M/out/$name.hits.jsonl" \
    > "$M/out/$name.collect.json" 2>> "$M/logs/$name.log"
  img python "$BIN/infernal_igr.py" permodel "$M/out/$name.out" "$M/out/$name.permodel.json" \
    > "$M/out/$name.permodel.sum.json" 2>> "$M/logs/$name.log"
  log "collected $name: $(cat "$M/out/$name.collect.json")"
}

case "$LANE" in
1)
  job prepare python "$BIN/infernal_igr.py" prepare "$TREE" /tmp/if2453/golden "$M/root" "$M/scratch" || exit 1
  wait_for "$M/lib/done"
  for b in A B; do
    [ -f "$M/fa/wg_$b.map.json" ] || img python "$BIN/infernal_igr.py" fasta "$M/root" "$CM" wg "$b" \
      "$M/fa/wg_$b.fna" "$M/fa/wg_$b.map.json" >> "$M/logs/fasta.log" 2>&1 || exit 1
  done
  for m in 0 50 200; do
    [ -f "$M/fa/igr_m$m.map.json" ] || img python "$BIN/infernal_igr.py" fasta "$M/root" "$CM" igr all \
      "$M/fa/igr_m$m.fna" "$M/fa/igr_m$m.map.json" "$m" >> "$M/logs/fasta.log" 2>&1 || exit 1
  done
  [ -f "$M/fa/wg_all.map.json" ] || img python "$BIN/infernal_igr.py" fasta "$M/root" "$CM" wg all \
    "$M/fa/wg_all.fna" "$M/fa/wg_all.map.json" >> "$M/logs/fasta.log" 2>&1 || exit 1
  touch "$M/fa/done"
  scan wgA_default "$CM" "$M/fa/wg_A"
  scan igr_m0_default "$CM" "$M/fa/igr_m0"
  scan igr_m200_default "$CM" "$M/fa/igr_m200"
  scan wgA_rfam "$CM" "$M/fa/wg_A" --rfam
  wait_for "$M/out/wgB_default.hits.jsonl"
  if [ ! -f "$M/lib/subset.cm" ]; then
    img python "$BIN/infernal_igr.py" subset "$TREE" "$M/out/wgA_default.hits.jsonl" "$M/out/wgB_default.hits.jsonl" \
      -- "$M/lib/subset.keys" > "$M/lib/subset.json" 2>> "$M/logs/subset.log" || exit 1
    img cmfetch -f -o "$M/lib/subset.cm.tmp" "$CM" "$M/lib/subset.keys" >> "$M/logs/subset.log" 2>&1 || exit 1
    mv "$M/lib/subset.cm.tmp" "$M/lib/subset.cm"
  fi
  scan wgA_subset "$M/lib/subset.cm" "$M/fa/wg_A"
  wait_for "$G04/region.map.json"
  scan wgAll_panel20 "$G04/panel.cm" "$M/fa/wg_all"
  scan igr_m200_rfam "$CM" "$M/fa/igr_m200" --rfam
  ;;
2)
  if [ ! -f "$M/lib/done" ]; then
    gunzip -c /tmp/if2453/golden/rfam/Rfam.cm.gz > "$M/lib/Rfam.cm.tmp" || exit 1
    mv "$M/lib/Rfam.cm.tmp" "$CM"
    sha256sum /tmp/if2453/golden/rfam/Rfam.cm.gz > "$M/lib/Rfam.cm.gz.sha256"
    touch "$M/lib/done"
  fi
  if [ ! -f "$M/fa/tiny.map.json" ]; then
    img python "$BIN/infernal_igr.py" tiny "$M/fa/tiny.fna" || exit 1
    echo '{"mode": "igr", "map": {"tiny": {"contig": "tiny", "genome": "tiny", "len": 1000, "circular": false, "start": 0, "end": 1000}}}' > "$M/fa/tiny.map.json"
  fi
  scan tiny_default "$CM" "$M/fa/tiny"
  wait_for "$M/fa/done"
  scan wgB_default "$CM" "$M/fa/wg_B"
  scan igr_m50_default "$CM" "$M/fa/igr_m50"
  scan wgB_rfam "$CM" "$M/fa/wg_B" --rfam
  job g04_emit python "$BIN/infernal_spikein.py" emit "$TREE" "$CM" "$G04"
  scan g04_default "$G04/panel.cm" "$G04/region"
  scan g04_rfam "$G04/panel.cm" "$G04/region" --rfam
  for v in default rfam; do
    img python "$BIN/infernal_spikein.py" score "$TREE" "$G04" "$M/out/g04_$v.hits.jsonl" > "$M/out/g04_$v.score.json" \
      2>> "$M/logs/g04_score.log"
  done
  wait_for "$M/lib/subset.cm"
  scan wgB_subset "$M/lib/subset.cm" "$M/fa/wg_B"
  scan igr_m50_rfam "$CM" "$M/fa/igr_m50" --rfam
  scan igr_m0_rfam "$CM" "$M/fa/igr_m0" --rfam
  ;;
3)
  # second plan, after wgB_rfam showed --rfam identical to the default at Z = 1e5 (every filter-stage count of
  # every model and every hit): lane 2 was stopped while it waited, and its --rfam jobs were dropped
  scan wgAll_panel20 "$G04/panel.cm" "$M/fa/wg_all"
  scan wgB_hmmonly "$CM" "$M/fa/wg_B" --hmmonly
  scan g04_hmmonly "$G04/panel.cm" "$G04/region" --hmmonly
  img python "$BIN/infernal_spikein.py" score "$TREE" "$G04" "$M/out/g04_hmmonly.hits.jsonl" \
    > "$M/out/g04_hmmonly.score.json" 2>> "$M/logs/g04_score.log"
  wait_for "$M/lib/subset.cm"
  scan wgB_subset "$M/lib/subset.cm" "$M/fa/wg_B"
  ;;
4)
  # replaces lane 1's igr_m200_rfam (stopped at its start): intergenic without the truncated-hit re-search
  scan igr_m0_notrunc "$CM" "$M/fa/igr_m0" --notrunc
  scan igr_m200_notrunc "$CM" "$M/fa/igr_m200" --notrunc
  ;;
*)
  echo "lane 1 to 4" >&2
  exit 2
  ;;
esac
log "lane done"
touch "$M/lane$LANE.done"
