#!/bin/sh
# One ladder stage as a detached unit:  systemd-run --user --unit=ladder-<key> --collect -p MemoryMax=200G sh stage.sh <KEY>
# The driver waits for its prerequisites itself; marker logs/<KEY>.done or logs/<KEY>.failed at the end.
W=/tmp/if2453/ladder
KEY=$1
exec >> $W/logs/$KEY.log 2>&1 < /dev/null
echo "stage $KEY wrapper started $(date) on $(hostname)"
rm -f $W/logs/$KEY.failed
if [ -e /stor/scratch/Wilke/if2453/all-the-data/STOP ]; then echo "STOP present"; touch $W/logs/$KEY.failed; exit 75; fi
nice -n 19 python3 $W/scripts/driver.py $KEY
RC=$?
echo "stage $KEY wrapper finished $(date) rc=$RC"
if [ $RC = 0 ]; then touch $W/logs/$KEY.done; else touch $W/logs/$KEY.failed; fi
exit $RC
