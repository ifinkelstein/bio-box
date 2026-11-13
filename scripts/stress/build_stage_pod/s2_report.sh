#!/bin/sh
# The read-only reports after s2_seqstats.py (run ON THE POD):
#
#     systemd-run --user --unit=s2-report --collect -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=100% \
#       nice -n 19 ionice -c 2 -n 7 sh /tmp/if2453/s2/bin/s2_report.sh
#
# s2_overlap_labels.py for HPGC and PLSDB (work/<src>/overlap_labels.tsv), then s2_plsdb_meta.py
# (work/plsdb/meta.json). The r10 root is bound read-only. Log: logs/report.log.
W=/tmp/if2453/s2
SIF=/tmp/if2453/images/bio-box_80d4c55.sif
R=/tmp/if2453/r10/root
exec >> $W/logs/report.log 2>&1 < /dev/null
run() { singularity exec --cleanenv --bind /tmp/if2453 --bind $R:$R:ro --env PYTHONPATH=$W/tree/src $SIF python "$@"; }
echo "report start $(date)"
for s in hpgc plsdb; do
  run $W/bin/s2_overlap_labels.py $R $W/work/stage_uids.tsv $W/work/$s/overlap.tsv $W/work/$s/overlap_labels.tsv
done
run $W/bin/s2_plsdb_meta.py $W/plsdb/meta $W/work/plsdb $W/work/stage_uids.tsv $R > $W/work/plsdb/.meta.json.tmp \
  && mv $W/work/plsdb/.meta.json.tmp $W/work/plsdb/meta.json
echo "report done rc=$? $(date)"
