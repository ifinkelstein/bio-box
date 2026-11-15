#!/bin/sh
# batch 3 (card F-02 session 12, DEC-44 (2)): candidate B with plasmids 0.1 / 0.3 over the 22 seeds, decoy transfer
# 4, 8, 16 per copy (decoys.py on scan2: a decoy fails when its two crown copies die before transferring, or when its
# partner module is held by one lineage from about t 1, so the rate must give 3 transfers before then).
T0=/tmp/if2453/f0212/tree-a
O=/tmp/if2453/f0212/scan3
cd "$T0" || exit 1
for D in 4.0 8.0 16.0; do
  [ -f /tmp/if2453/f0212/STOP ] && exit 0
  python scripts/stress/f0210_seed_scan.py scan "$O" "e_t0.045_l0.076_dec$D" --jobs 7 --seeds default,20261113,n2 \
    --set genomes.transfer=0.045 --set genomes.duplication=0.012 --set genomes.loss=0.076 \
    --set genomes.origination_per_root_family=0.035 --set hubs.transfer_multiplier=44.444444 \
    --set plasmids.transfer=0.1 --set plasmids.loss=0.3 --set "decoys.transfer=$D"
done
