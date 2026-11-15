#!/bin/sh
# batch 2 (card F-02 session 12, DEC-44 (2)): candidate B with plasmids 0.1 / 0.3 at the seven seeds where the
# fewest decoy transfers was below 3 at decoy transfer 1.2 (scan1), kept for decoys.py; decoy transfer 1.2, 2.0, 4.0.
T0=/tmp/if2453/f0212/tree-a
O=/tmp/if2453/f0212/scan2
S=9254619518473077140,4701969716356831982,14843108908847119356,1995372076782357328,13514522956726874792,4734920099712499969,3830344756309592292
cd "$T0" || exit 1
for D in 1.2 2.0 4.0; do
  [ -f /tmp/if2453/f0212/STOP ] && exit 0
  python scripts/stress/f0210_seed_scan.py scan "$O" "z_dec$D" --jobs 7 --seeds "$S" --keep \
    --set genomes.transfer=0.045 --set genomes.duplication=0.012 --set genomes.loss=0.076 \
    --set genomes.origination_per_root_family=0.035 --set hubs.transfer_multiplier=44.444444 \
    --set plasmids.transfer=0.1 --set plasmids.loss=0.3 --set "decoys.transfer=$D"
done
