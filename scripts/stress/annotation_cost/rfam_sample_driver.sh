#!/bin/bash
# Card 0-08a: the unattended driver of the sample scan (run ON THE POD HOST, detached, in its own capped unit):
#
#   systemd-run --user --unit=m08a-driver --collect --quiet -p MemoryMax=4G -p MemorySwapMax=0 -p CPUQuota=100% \
#     -p StandardOutput=append:/tmp/if2453/m08a/driver.log -p StandardError=append:/tmp/if2453/m08a/driver.log \
#     nice -n 19 bash /tmp/if2453/m08a/tree/scripts/stress/annotation_cost/rfam_sample_driver.sh
#
# Works through the processes of $M/sample/shards.txt (written by `rfam_sample.py sample`). Each process k is a
# unit m08a-s<k> (MemoryMax=16G, CPUQuota=200%, nice 19, cmsearch --cpu 2 in the image: rfam_sample_shard.sh).
# At most $M/MAXP processes run at once (default 2; the file may be edited while the driver runs, at most 4).
# Before EVERY launch: the sum of `ps -u if2453 -o pcpu=` (the ps line left out) plus 400% must be at most 4000%;
# otherwise it polls every 4 minutes. A HALT file in $M, or a STOP file in the data root (AGENTS 7.4, only
# read), ends the driver before its next launch (running processes finish; it writes no DONE). A process is
# done when out/s<k>.ok exists; one that failed (out/s<k>.fail, or its unit gone without a marker after its
# second launch) is not launched again in this run. When every process is done the driver builds the report
# inputs (`rfam_sample_report.py`) and writes DONE; when nothing is left to run and a process failed, it
# writes FAILED (listing them). Run again, it resumes: done processes are skipped, failed ones retried.
set -uo pipefail
M=/tmp/if2453/m08a
IMG=/tmp/if2453/images/bio-box_80d4c55.sif
TREE=$M/tree
BIN=$TREE/scripts/stress/annotation_cost
DATA_STOP=/tmp/if2453/r10/root/STOP
mkdir -p "$M/logs" "$M/runs" "$M/out" "$M/work"
rm -f "$M/DONE" "$M/FAILED"
for f in "$M"/out/s*.fail; do [ -e "$f" ] && mv "$f" "$f.prev"; done
for f in "$M"/out/s*.launched; do [ -e "$f" ] && rm -f "$f"; done

log() { echo "driver: $(date -Is) $*"; }

halted() {
  if [ -e "$M/HALT" ]; then log "HALT file: exiting"; return 0; fi
  if [ -e "$DATA_STOP" ]; then log "STOP file in the data root: exiting"; return 0; fi
  return 1
}

account() { ps -u if2453 -o pcpu=,comm= | awk '$2 != "ps" {s+=$1} END {printf "%d", s}'; }

maxp() {
  local p=2
  [ -f "$M/MAXP" ] && p=$(tr -dc 0-9 < "$M/MAXP")
  [ -n "$p" ] || p=2
  [ "$p" -gt 4 ] && p=4
  [ "$p" -lt 1 ] && p=1
  echo "$p"
}

active() { systemctl --user is-active --quiet "m08a-s$1"; }

log "start: $(wc -l < "$M/sample/shards.txt") processes planned"
while :; do
  running=0 pending=() failed=() done_n=0
  while read -r k bases; do
    if [ -e "$M/out/s$k.ok" ]; then done_n=$((done_n + 1)); continue; fi
    if active "$k"; then running=$((running + 1)); continue; fi
    if [ -e "$M/out/s$k.fail" ]; then failed+=("$k"); continue; fi
    n=0
    [ -f "$M/out/s$k.launched" ] && n=$(cat "$M/out/s$k.launched")
    if [ "$n" -ge 2 ]; then failed+=("$k"); continue; fi
    pending+=("$k")
  done < "$M/sample/shards.txt"
  total=$(wc -l < "$M/sample/shards.txt")
  if [ "$done_n" -eq "$total" ]; then
    log "all $total processes done; building the report"
    singularity exec --cleanenv --bind /tmp/if2453 --env PYTHONPATH="$TREE/src" --pwd "$M/work" "$IMG" \
      python "$BIN/rfam_sample_report.py" "$M" "$TREE" "$M/report" >> "$M/logs/report.log" 2>&1
    rc=$?
    if [ "$rc" -eq 0 ]; then log "report written to $M/report"; touch "$M/DONE"; exit 0; fi
    log "report exit $rc (see logs/report.log)"
    echo "report exit $rc" > "$M/FAILED"
    exit 1
  fi
  if [ "$running" -eq 0 ] && [ "${#pending[@]}" -eq 0 ]; then
    log "nothing left to run; failed: ${failed[*]}"
    echo "failed processes: ${failed[*]}" > "$M/FAILED"
    exit 1
  fi
  if halted; then exit 0; fi
  if [ "${#pending[@]}" -gt 0 ] && [ "$running" -lt "$(maxp)" ]; then
    used=$(account)
    if [ $((used + 400)) -le 4000 ]; then
      k=${pending[0]}
      # a process may have ended between the listing and now
      if [ -e "$M/out/s$k.ok" ] || [ -e "$M/out/s$k.fail" ] || active "$k"; then sleep 5; continue; fi
      n=0
      [ -f "$M/out/s$k.launched" ] && n=$(cat "$M/out/s$k.launched")
      echo $((n + 1)) > "$M/out/s$k.launched"
      log "launch s$k (done $done_n, running $running, pending ${#pending[@]}): account ${used}% + 400% <= 4000%"
      systemd-run --user --unit="m08a-s$k" --collect --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% \
        -p StandardOutput="append:$M/logs/s$k.log" -p StandardError="append:$M/logs/s$k.log" \
        nice -n 19 bash "$BIN/rfam_sample_shard.sh" "$k"
      sleep 30
      continue
    fi
    log "waiting (done $done_n, running $running, pending ${#pending[@]}): account ${used}% + 400% > 4000%"
  fi
  sleep 240
done
