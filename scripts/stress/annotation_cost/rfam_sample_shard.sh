#!/bin/bash
# Card 0-08a: one scan process of the plan (run ON THE POD HOST by rfam_sample_driver.sh, in its own unit):
#
#   systemd-run --user --unit=m08a-s<k> --collect --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% \
#     ... nice -n 19 bash rfam_sample_shard.sh <k>
#
# The producer's whole-contig cmsearch (biobox.l1.ncrna._cmsearch_argv: --cpu, -Z 1e5, --cut_ga, --noali,
# --tformat fasta, --tblout; not --notrunc) over the full Rfam 15.1 library, with -o <file> instead of /dev/null
# so that each model's own "# CPU time" line is kept, under timed.py (wall, CPU-minutes, peak RSS) inside the
# image. Then the hits in contig coordinates (infernal_igr.py collect: the producer's leftmost-base rule) and the
# per-model CPU (infernal_igr.py permodel). Ends with out/s<k>.ok, or out/s<k>.fail holding the exit code.
set -uo pipefail
M=/tmp/if2453/m08a
IMG=/tmp/if2453/images/bio-box_80d4c55.sif
TREE=$M/tree
BIN=$TREE/scripts/stress/annotation_cost
CM=$M/lib/Rfam.cm
K=${1:?shard number, two digits}
S=s$K
if [ -e "$M/out/$S.ok" ]; then echo "shard $S: already done"; exit 0; fi
rm -f "$M/out/$S.fail"

img() {
  singularity exec --cleanenv --bind /tmp/if2453 --env PYTHONPATH="$TREE/src" --env OMP_NUM_THREADS=2 \
    --pwd "$M/work" "$IMG" "$@"
}

echo "shard $S: start $(date -Is)"
img python "$BIN/timed.py" "$M/runs/$S.json" -- cmsearch --cpu 2 -Z 100000 --cut_ga --noali --tformat fasta \
  --tblout "$M/out/$S.tbl" -o "$M/out/$S.out" "$CM" "$M/fa/$S.fna"
rc=$?
if [ "$rc" -ne 0 ]; then echo "$rc" > "$M/out/$S.fail"; echo "shard $S: cmsearch exit $rc"; exit "$rc"; fi
img python "$BIN/infernal_igr.py" collect "$M/out/$S.tbl" "$M/fa/$S.map.json" "$M/out/$S.hits.jsonl" \
  > "$M/out/$S.collect.json"
rc=$?
if [ "$rc" -ne 0 ]; then echo "collect $rc" > "$M/out/$S.fail"; exit "$rc"; fi
img python "$BIN/infernal_igr.py" permodel "$M/out/$S.out" "$M/out/$S.permodel.json" > "$M/out/$S.permodel.sum.json"
rc=$?
if [ "$rc" -ne 0 ]; then echo "permodel $rc" > "$M/out/$S.fail"; exit "$rc"; fi
gzip -f "$M/out/$S.out"
touch "$M/out/$S.ok"
echo "shard $S: done $(date -Is) $(cat "$M/out/$S.collect.json")"
