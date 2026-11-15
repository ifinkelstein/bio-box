#!/bin/sh
# batch 1 (card F-02 session 12, DEC-44): session 11's candidate B (duplication 0.012, origination 0.035) with the
# background transfer T, loss derived as (-0.01 + 2 dup + 3 T) / 2 (session 11's rule), hub multiplier 2.0 / T,
# plasmids 0.1 / 0.3 (DEC-44 (3)) and the decoy transfer rate D (DEC-44 (2)). 22 seeds, 7 processes.
T0=/tmp/if2453/f0212/tree-a
O=/tmp/if2453/f0212/scan1
cd "$T0" || exit 1
for C in "0.045,0.076,44.444444,0.4" "0.045,0.076,44.444444,0.8" "0.045,0.076,44.444444,1.2" \
         "0.05,0.082,40,0.8" "0.05,0.082,40,1.2"; do
  [ -f /tmp/if2453/f0212/STOP ] && exit 0
  TR=$(echo "$C" | cut -d, -f1)
  L=$(echo "$C" | cut -d, -f2)
  H=$(echo "$C" | cut -d, -f3)
  D=$(echo "$C" | cut -d, -f4)
  python scripts/stress/f0210_seed_scan.py scan "$O" "e_t${TR}_l${L}_dec$D" --jobs 7 --seeds default,20261113,n2 \
    --set "genomes.transfer=$TR" --set genomes.duplication=0.012 --set "genomes.loss=$L" \
    --set genomes.origination_per_root_family=0.035 --set "hubs.transfer_multiplier=$H" \
    --set plasmids.transfer=0.1 --set plasmids.loss=0.3 --set "decoys.transfer=$D"
done
