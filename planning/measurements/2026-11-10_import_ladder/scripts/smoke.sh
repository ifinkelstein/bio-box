#!/bin/sh
# Smoke test of step.py's stage-A path on the first 300 INPHARED records, in its own work dir.
W=/tmp/if2453/ladder
SM=$W/smoke
mkdir -p $SM
exec > $SM/smoke.log 2>&1 < /dev/null
if [ -e /stor/scratch/Wilke/if2453/all-the-data/STOP ]; then echo "STOP present"; exit 75; fi
zcat /tmp/if2453/pilot_s1/src/14Apr2025_genomes.fa.gz | awk '/^>/{n++} n>300{exit} {print}' > $SM/src.fa
for ST in ingest genome_input genecall igr tapes cluster tokens postings dim catalog equiv queries; do
  echo "=== $ST $(date)"
  sh $W/scripts/bbn.sh $W/scripts/step.py A $ST --threads 4 --work $SM --src $SM/src.fa
  echo "=== $ST rc=$?"
done
echo "smoke finished $(date)"
