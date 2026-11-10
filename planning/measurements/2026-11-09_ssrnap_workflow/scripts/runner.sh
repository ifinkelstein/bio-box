#!/bin/sh
# detached: systemd-run --user --unit=ssrnap-runner --collect -p MemoryMax=20G sh /tmp/if2453/ssrnap/scripts/runner.sh
cd /tmp/if2453/ssrnap/scripts
/usr/bin/time -v sh /tmp/if2453/ssrnap/scripts/run.sh python /tmp/if2453/ssrnap/scripts/s8_runner.py "$@" > /tmp/if2453/ssrnap/logs/runner.log 2>&1
echo "exit $?" >> /tmp/if2453/ssrnap/logs/runner.log
