#!/bin/sh
# Stage H1's build-stage inputs from the HPGC download, then two checks; run ON THE POD:
#
#     systemd-run --user --unit=s2-h1-inputs --collect -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=100% \
#       nice -n 19 ionice -c 2 -n 7 sh /tmp/if2453/s2/bin/s2_h1.sh
#
# 1. s2_h1_inputs.py: genomes/H1/HPGC.7230.genomes.final.fasta (the download decompressed byte for
#    byte) and inputs/H1/{genomes.json, objects.json, dropped_records.tsv}.
# 2. make_inputs.py phage on the same download, into work/make_inputs_phage (a check that H1's
#    genomes.json is what make_inputs.py itself gives).
# 3. s2_check_inputs.py: the inputs step's own parsing and checks, writing nothing.
# Log: logs/h1.log. Image bio-box_80d4c55.sif, code PYTHONPATH=tree/src (plan/task-cards at 1038dd8).
W=/tmp/if2453/s2
SIF=/tmp/if2453/images/bio-box_80d4c55.sif
SRC=$W/dl/HPGC.7230.genomes.final.fasta.gz
exec >> $W/logs/h1.log 2>&1 < /dev/null
for f in /tmp/if2453/r10/root/STOP $W/STOP; do [ -e "$f" ] && { echo "STOP present: $f"; exit 3; }; done
run() { singularity exec --cleanenv --bind /tmp/if2453 --env PYTHONPATH=$W/tree/src $SIF python "$@"; }
mkdir -p $W/inputs
echo "h1 inputs start $(date)"
run $W/bin/s2_h1_inputs.py $SRC $W/genomes/H1/HPGC.7230.genomes.final.fasta $W/inputs/H1 2026-11-11T23:09:28Z || exit 1
echo "make_inputs.py phage check start $(date)"
rm -rf $W/work/make_inputs_phage
run $W/bin/make_inputs.py phage $SRC 7230 $W/work/make_inputs_phage || exit 1
cmp $W/work/make_inputs_phage/genomes.json $W/inputs/H1/genomes.json && echo "genomes.json identical to make_inputs.py phage"
cmp $W/work/make_inputs_phage/phage_first7230.fa $W/genomes/H1/HPGC.7230.genomes.final.fasta && echo "FASTA copy identical to make_inputs.py phage"
rm -f $W/work/make_inputs_phage/phage_first7230.fa
echo "check start $(date)"
run $W/bin/s2_check_inputs.py H1 $W/inputs/H1/genomes.json $W/inputs/H1/objects.json
echo "check rc=$? $(date)"
