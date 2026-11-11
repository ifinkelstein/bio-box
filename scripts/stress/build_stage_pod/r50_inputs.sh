#!/bin/sh
# Write the build-stage inputs (genomes.json, objects.json) of stages P3..P6 of the 50 Gbp build,
# one stage after another, in the image, inside the r50-inputs unit (run this ON THE POD):
#
#     systemd-run --user --unit=r50-inputs --collect -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=100% \
#       nice -n 19 sh /tmp/if2453/r50/bin/r50_inputs.sh
#
# Each stage goes to /tmp/if2453/r50/inputs/<S>/ (written under <S>.tmp, renamed when complete);
# objects carry source_release datasets-v2-2026-11-11, the date of this fetch. Log:
# /tmp/if2453/r50/logs/inputs.log (one JSON report per stage); at the end inputs/DONE or inputs/FAILED.
W=/tmp/if2453/r50
SIF=/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif
exec >> $W/logs/inputs.log 2>&1 < /dev/null
mkdir -p $W/inputs
rm -f $W/inputs/DONE $W/inputs/FAILED
for S in P3 P4 P5 P6; do
  if [ -e $W/inputs/$S/genomes.json ]; then echo "$S already written"; continue; fi
  rm -rf $W/inputs/$S.tmp
  echo "$S start $(date)"
  singularity exec --cleanenv --bind /tmp/if2453 --env PYTHONPATH=$W/tree/src $SIF \
    python $W/bin/make_inputs.py prok $W/genomes/$S/lock.json $W/selection/stage_$S.tsv $W/inputs/$S.tmp \
    --source-release datasets-v2-2026-11-11
  RC=$?
  if [ $RC != 0 ]; then echo "$S failed rc=$RC $(date)"; echo "$S rc=$RC" > $W/inputs/FAILED; exit $RC; fi
  mv $W/inputs/$S.tmp $W/inputs/$S
  echo "$S done $(date)"
done
date > $W/inputs/DONE
