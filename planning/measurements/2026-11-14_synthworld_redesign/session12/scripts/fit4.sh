#!/bin/sh
# batch 4 (card F-02 session 12, DEC-44 (2), (5)): candidate B with plasmids 0.1 / 0.3 and decoy transfer 8.0 over the
# 22 seeds at background transfer 0.05 and 0.055 (loss (-0.01 + 2 x 0.012 + 3 T) / 2, hub multiplier 2.0 / T,
# origination 0.035); 5 processes (3 threads are the igrs1 units).
T0=/tmp/if2453/f0212/tree-a
O=/tmp/if2453/f0212/scan4
cd "$T0" || exit 1
for C in "0.05,0.082,40" "0.055,0.0895,36.363636"; do
  [ -f /tmp/if2453/f0212/STOP ] && exit 0
  TR=$(echo "$C" | cut -d, -f1)
  L=$(echo "$C" | cut -d, -f2)
  H=$(echo "$C" | cut -d, -f3)
  python scripts/stress/f0210_seed_scan.py scan "$O" "e_t${TR}_l${L}_dec8.0" --jobs 5 --seeds default,20261113,n2 \
    --set "genomes.transfer=$TR" --set genomes.duplication=0.012 --set "genomes.loss=$L" \
    --set genomes.origination_per_root_family=0.035 --set "hubs.transfer_multiplier=$H" \
    --set plasmids.transfer=0.1 --set plasmids.loss=0.3 --set decoys.transfer=8.0
done
