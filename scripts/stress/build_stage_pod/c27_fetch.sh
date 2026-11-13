#!/bin/sh
# Card C-27 (3): fetch S2C phase 2's ANI genomes on wilkcomp02, one detached unit per batch.
#
#     sh /tmp/if2453/c27/bin/c27_fetch.sh launch B00 chain    # start (or resume) at batch B00, then each next one
#     sh /tmp/if2453/c27/bin/c27_fetch.sh launch B05          # one batch only
#
# launch refuses (exit 1, reason in $W/logs/launch.log and $W/genomes/HALT) when STOP is present, when another
# c27-fetch-* unit is active, or when /tmp would keep less than 250 GB free after the batch's expected size
# (0.30 bytes per bp, c27_fetch.py need). A batch already DONE is skipped. The unit is c27-fetch-<kk> (kk the
# batch number), at nice 19, MemoryMax=8G, CPUQuota=200% (network-bound; at most 2 CPUs in all). It runs
# `run`: c27_fetch.py fetch for the batch, logs df and du before and after, and with `chain`, when the batch
# ended DONE, launches the next batch through the same checks. A batch that ends FAILED (refusals, STOP, or more
# than 5% of it not served; B00, the 7 phase-1 retries, may fail whole) stops the chain.
# Progress: $W/logs/fetch.log (one line per genome); per batch $W/genomes/<B>/DONE or FAILED; $W/logs/disk.log.
# Writes only under /tmp/if2453/c27.
set -u
W=/tmp/if2453/c27
BIN=$W/bin
STOPF=/stor/scratch/Wilke/if2453/all-the-data/STOP
MIN_FREE_KB=244140625  # 250 GB = 250e9 bytes, in KiB
CMD=$1
B=$2
CHAIN=${3:-}
K=$(echo "$B" | sed 's/^B//')
UNIT=c27-fetch-$K
mkdir -p "$W/logs" "$W/genomes"

log() { echo "$(date -Is) $*" >> "$W/logs/launch.log"; echo "$*"; }
disk() { echo "$(date -Is) $1 $B df_avail_kb=$(df -Pk /tmp | awk 'NR==2 {print $4}') du_genomes_b=$(du -sb "$W/genomes" | cut -f1)" >> "$W/logs/disk.log"; }
halt() { log "HALT $B: $*"; echo "$(date -Is) $B: $*" > "$W/genomes/HALT"; exit 1; }

case "$CMD" in
launch)
  TSV=$W/batches/batch_$B.tsv
  [ -f "$TSV" ] || { log "no batch $B: end of chain"; exit 0; }
  if [ -f "$W/genomes/$B/DONE" ]; then
    log "$B already DONE"
    if [ "$CHAIN" = chain ]; then
      NEXT=$(printf 'B%02d' "$(expr "$K" + 1)")
      exec sh "$BIN/c27_fetch.sh" launch "$NEXT" chain
    fi
    exit 0
  fi
  [ ! -e "$STOPF" ] || halt "STOP present"
  ACTIVE=$(systemctl --user list-units --type=service --state=active --plain --no-legend 'c27-fetch-*' | awk '{print $1}' | grep -v "^${C27_SELF:-none}.service$")
  [ -z "$ACTIVE" ] || halt "already active: $ACTIVE"
  AVAIL_KB=$(df -Pk /tmp | awk 'NR==2 {print $4}')
  NEED_KB=$(python3 -I "$BIN/c27_fetch.py" need "$TSV")
  [ $((AVAIL_KB - NEED_KB)) -ge $MIN_FREE_KB ] || halt "/tmp would keep $((AVAIL_KB - NEED_KB)) KiB free, under 250 GB"
  rm -f "$W/genomes/HALT"
  systemctl --user reset-failed "$UNIT" 2>/dev/null
  log "launch $UNIT for $B: avail_kb=$AVAIL_KB need_kb=$NEED_KB chain=$CHAIN"
  exec systemd-run --user --unit="$UNIT" --collect --quiet \
    -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=200% \
    -p StandardOutput=append:$W/logs/$UNIT.out -p StandardError=append:$W/logs/$UNIT.out \
    --setenv=C27_SELF="$UNIT" --setenv=XDG_RUNTIME_DIR="${XDG_RUNTIME_DIR:-/run/user/$(id -u)}" \
    --setenv=DBUS_SESSION_BUS_ADDRESS="${DBUS_SESSION_BUS_ADDRESS:-unix:path=/run/user/$(id -u)/bus}" \
    nice -n 19 sh "$BIN/c27_fetch.sh" run "$B" "$CHAIN"
  ;;
run)
  FRAC=0.05
  [ "$B" != B00 ] || FRAC=1.0
  disk start
  python3 -I "$BIN/c27_fetch.py" fetch "$W/batches/batch_$B.tsv" "$W/genomes/$B" "$W/logs/fetch.log" --max-fail-frac $FRAC
  RC=$?
  disk end
  echo "$(date -Is) $B exit $RC" >> "$W/logs/launch.log"
  if [ $RC -eq 0 ] && [ "$CHAIN" = chain ]; then
    NEXT=$(printf 'B%02d' "$(expr "$K" + 1)")
    exec sh "$BIN/c27_fetch.sh" launch "$NEXT" chain
  fi
  exit $RC
  ;;
*)
  echo "usage: c27_fetch.sh launch|run <B> [chain]"; exit 2 ;;
esac
