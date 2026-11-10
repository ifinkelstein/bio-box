#!/bin/sh
# detached: systemd-run --user --unit=ssrnap-pack --collect -p MemoryMax=20G sh /tmp/if2453/ssrnap/scripts/pack.sh [--code15-as-11]
sh /tmp/if2453/ssrnap/scripts/run.sh python /tmp/if2453/ssrnap/scripts/s7_pack.py "$@" > /tmp/if2453/ssrnap/logs/pack2.log 2>&1
echo "exit $?" >> /tmp/if2453/ssrnap/logs/pack2.log
