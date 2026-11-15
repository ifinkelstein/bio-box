#!/bin/sh
# final (card F-02 session 12): the preset as committed (tree-b), no overrides, over the 22 seeds with N0-full.
T0=/tmp/if2453/f0212/tree-b
O=/tmp/if2453/f0212/scanFinal
cd "$T0" || exit 1
[ -f /tmp/if2453/f0212/STOP ] && exit 0
python scripts/stress/f0210_seed_scan.py scan "$O" final12 --jobs "${1:-6}" --seeds default,20261113,n2 --n0
