#!/bin/sh
# the default world S1-S11 of the preset as committed (tree-b) and its summary (card F-02 session 13). 1 process.
T0=/tmp/if2453/f0213/tree-b
O=/tmp/if2453/f0213/scanFinal
cd "$T0" || exit 1
[ -f /tmp/if2453/f0213/STOP ] && exit 0
mkdir -p "$O"
rm -rf /tmp/if2453/f0213/world-final
python scripts/stress/f0210_world_summary.py build /tmp/if2453/f0213/world-final 20261112
python scripts/stress/f0210_world_summary.py summary /tmp/if2453/f0213/world-final "$O/default_world_summary_after13.json"
