#!/bin/sh
# fit1.sh OUTDIR TAG PARTNERS JOBS RATE... (card F-02 session 13, DEC-48 (1)): the session-12 preset with the decoy
# partners PARTNERS (a JSON list of module indices) and decoy transfer RATE, over the 22 seeds, one scan per rate,
# label TAG_dec<RATE>. Member transfer levels are (i // 3) % 3 of module i, so 6, 7, 8, 15, 16, 17, 24 are the top
# level (0.4); session 12's partners 1, 9, 11, 19, 22 sat at the two lowest. JOBS processes of one thread each.
T0=/tmp/if2453/f0213/tree-a
O=$1
TAG=$2
P=$3
JOBS=$4
shift 4
cd "$T0" || exit 1
for D in "$@"; do
  [ -f /tmp/if2453/f0213/STOP ] && exit 0
  python scripts/stress/f0210_seed_scan.py scan "$O" "${TAG}_dec$D" --jobs "$JOBS" --seeds default,20261113,n2 \
    --set "decoys.partner_modules=$P" --set "decoys.transfer=$D"
done
