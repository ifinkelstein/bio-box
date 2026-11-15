#!/bin/sh
# igrs1 (card F-02 session 12, DEC-44 (5)): default-seed S1-S5 worlds for the cross-phylum IGR count and the motif
# scan, candidate B with plasmids 0.1 / 0.3 and decoy transfer 8.0, background transfer T (loss (-0.01 + 2 x 0.012
# + 3 T) / 2, hub multiplier 2.0 / T, origination 0.035); one unit per candidate.
R=/tmp/if2453/f0212/scripts/run.sh
T0=/tmp/if2453/f0212/tree-a
for C in "0.045,0.076,44.444444" "0.05,0.082,40" "0.055,0.0895,36.363636"; do
  TR=$(echo "$C" | cut -d, -f1)
  L=$(echo "$C" | cut -d, -f2)
  H=$(echo "$C" | cut -d, -f3)
  sh "$R" "f0212-igr-t$TR" "$T0" "/tmp/if2453/f0212/logs/igr_t$TR.log" python /tmp/if2453/f0212/scripts/igr.py \
    "/tmp/if2453/f0212/igr/t$TR" 20261112 "genomes.transfer=$TR" genomes.duplication=0.012 "genomes.loss=$L" \
    genomes.origination_per_root_family=0.035 "hubs.transfer_multiplier=$H" plasmids.transfer=0.1 plasmids.loss=0.3 \
    decoys.transfer=8.0
done
