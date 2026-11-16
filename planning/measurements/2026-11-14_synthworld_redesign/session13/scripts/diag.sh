#!/bin/sh
# diag (card F-02 session 13): why a decoy has fewer than 3 transfers with its partner at the top level. Keeps the
# backbones of p2_dec1.2 at four of its six failing seeds and replays the decoys with session 12's decoys.py
# (transfers, losses, when the last copy died, when the partner was last held by at most one lineage).
T0=/tmp/if2453/f0213/tree-a
O=/tmp/if2453/f0213/diag
cd "$T0" || exit 1
python scripts/stress/f0210_seed_scan.py scan "$O" p2_dec1.2_keep --jobs 4 --keep \
  --seeds 9254619518473077140,7689157899245528466,15390672246448918260,13514522956726874792 \
  --set "decoys.partner_modules=[6,7,15,16,24]" --set decoys.transfer=1.2
python planning/measurements/2026-11-14_synthworld_redesign/session12/scripts/decoys.py "$O"/worlds/p2_dec1.2_keep/*
