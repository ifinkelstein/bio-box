#!/bin/sh
# fit4 (card F-02 session 13, DEC-48 (1) fallback, second attempt): l1b at 6.0 holds every scan check, 202 cross-phylum
# IGR pairs and the null module at the default seed, but at 20261113 the null module's positive control on flank
# contexts fails (planted pairs above tau 2 of 75). The next rates up, with --n0.
T0=/tmp/if2453/f0213/tree-a
O=/tmp/if2453/f0213/scan4
cd "$T0" || exit 1
for C in "l1b [5,12,14,22,23] 7.0" "l1a [4,12,13,21,22] 7.0" "l1b [5,12,14,22,23] 8.0"; do
  [ -f /tmp/if2453/f0213/STOP ] && exit 0
  set -- $C
  python scripts/stress/f0210_seed_scan.py scan "$O" "n_$1_dec$3" --jobs 4 --seeds default,20261113,n2 --n0 \
    --set "decoys.partner_modules=$2" --set "decoys.transfer=$3"
done
