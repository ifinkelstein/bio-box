#!/bin/sh
# Smoke test of the stage-B path (and Bp's pfam_plan, cluster_copy) on the smoke root plus 3 P1 genomes.
W=/tmp/if2453/ladder
SM=$W/smoke2
mkdir -p $SM/genomes/P1 $SM/selection
exec > $SM/smoke.log 2>&1 < /dev/null
if [ -e /stor/scratch/Wilke/if2453/all-the-data/STOP ]; then echo "STOP present"; exit 75; fi
cp -a $W/smoke/root $SM/root
cp $W/selection/stage_P1.tsv $SM/selection/
python3 -c "
import json
d = json.load(open('$W/genomes/P1/lock.json'))
keep = dict(sorted(d['genomes'].items())[:3])
json.dump({'genomes': keep}, open('$SM/genomes/P1/lock.json', 'w'), indent=1)
"
for ST in ingest genome_input genecall igr tapes cluster tokens postings dim catalog appendonly queries; do
  echo "=== B $ST $(date)"
  sh $W/scripts/bbn.sh $W/scripts/step.py B $ST --threads 4 --work $SM
  echo "=== B $ST rc=$?"
done
for ST in pfam_plan cluster_copy; do
  echo "=== Bp $ST $(date)"
  sh $W/scripts/bbn.sh $W/scripts/step.py Bp $ST --threads 4 --work $SM
  echo "=== Bp $ST rc=$?"
done
echo "smoke finished $(date)"
