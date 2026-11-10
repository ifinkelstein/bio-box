#!/bin/sh
# Download the stage B (P1) then stage C (P2) genomes; markers logs/fetch-<S>.done or .failed.
W=/tmp/if2453/ladder
exec >> $W/logs/fetch.log 2>&1 < /dev/null
echo "$$" > $W/run/fetch.pid
echo "1 $$" > $W/run/fetch.threads
for S in P1 P2; do
  rm -f $W/logs/fetch-$S.failed
  if [ -e $W/logs/fetch-$S.done ]; then echo "fetch $S already done"; continue; fi
  if [ -e /stor/scratch/Wilke/if2453/all-the-data/STOP ]; then echo "STOP present"; touch $W/logs/fetch-$S.failed; exit 75; fi
  echo "fetch $S started $(date)"
  nice -n 19 python3 $W/scripts/fetch.py $W/selection/stage_$S.tsv $W/genomes/$S
  RC=$?
  echo "fetch $S finished $(date) rc=$RC"
  if [ $RC = 0 ]; then touch $W/logs/fetch-$S.done; else touch $W/logs/fetch-$S.failed; fi
done
rm -f $W/run/fetch.threads
