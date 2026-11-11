#!/bin/sh
# Launch (or resume) the 50 Gbp build's genome download on wilkcomp02 as a detached unit.
#
#     sh /tmp/if2453/r50/bin/r50_fetch.sh
#
# Refuses when a unit r50-fetch is active, when STOP is present, or when / would keep less than
# 200 GB free after the remaining download (estimate: 41 GB for all four stages). Progress: one line
# per genome in /tmp/if2453/r50/fetch.log; at the end /tmp/if2453/r50/DONE or /tmp/if2453/r50/FAILED.
# Rerunning resumes: genomes already in a stage's lock.json are not fetched again.
W=/tmp/if2453/r50
if systemctl --user is-active --quiet r50-fetch; then echo "r50-fetch is already active"; exit 1; fi
if [ -e /stor/scratch/Wilke/if2453/all-the-data/STOP ]; then echo "STOP present"; exit 75; fi
AVAIL_KB=$(df -Pk /tmp | awk 'NR==2 {print $4}')
if [ "$AVAIL_KB" -lt 251658240 ]; then echo "less than 240 GiB free on /tmp: $AVAIL_KB KiB"; exit 1; fi
mkdir -p $W/genomes
systemctl --user reset-failed r50-fetch 2>/dev/null
exec systemd-run --user --unit=r50-fetch --collect --quiet \
  -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=200% \
  -p StandardOutput=append:$W/logs/r50-fetch.out -p StandardError=append:$W/logs/r50-fetch.out \
  nice -n 19 python3 -I $W/bin/r50_fetch.py $W/selection $W/genomes $W/fetch.log --stages P3,P4,P5,P6
