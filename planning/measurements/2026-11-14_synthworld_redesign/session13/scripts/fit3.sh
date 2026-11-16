#!/bin/sh
# fit3 (card F-02 session 13, DEC-48 (1) fallback search): after fit1 and fit2 no placement holds "every decoy >= 3
# transfers at every seed" below rate 4.0 (l1b), and at every rate that holds decoys are 34-47 genes per genome. The
# smallest holding rate is looked for between 2.0 and 6.0 at the two middle-level placements, with --n0 (N0-full and
# N0-strict: the background size ratio and the root-presence KS at every seed, for the 22-seed bound).
T0=/tmp/if2453/f0213/tree-a
O=/tmp/if2453/f0213/scan3
cd "$T0" || exit 1
for C in "l1b [5,12,14,22,23] 3.0" "l1a [4,12,13,21,22] 3.0" "l1b [5,12,14,22,23] 4.0" "l1a [4,12,13,21,22] 5.0" \
         "l1b [5,12,14,22,23] 5.0" "l1b [5,12,14,22,23] 6.0"; do
  [ -f /tmp/if2453/f0213/STOP ] && exit 0
  set -- $C
  python scripts/stress/f0210_seed_scan.py scan "$O" "n_$1_dec$3" --jobs 4 --seeds default,20261113,n2 --n0 \
    --set "decoys.partner_modules=$2" --set "decoys.transfer=$3"
done
