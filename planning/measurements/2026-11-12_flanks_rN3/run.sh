#!/bin/bash
# Card W-25 real run: one bio-box command in the image, inside the systemd unit named by $1 (started by the caller:
#   systemd-run --user --unit w25-<name> -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=400% -p RemainAfterExit=yes
#     -p StandardOutput=file:... nice -n 19 bash /tmp/if2453/w25/run.sh w25-<name> <bio-box args...>)
# Reads the root read-only; writes only under /tmp/if2453/w25. Prints wall time and the unit's usage at the end.
UNIT=$1
shift
if [ -e /tmp/if2453/r10/root/STOP ]; then echo "w25: STOP file present, not running"; exit 0; fi
T0=$(date +%s.%N)
/usr/bin/time -v singularity exec --cleanenv \
  --env PYTHONPATH=/tmp/if2453/w25/tree/src:/tmp/if2453/pilot_s1/pydeps \
  --env TMPDIR=/tmp/if2453/w25/tmp --env OMP_NUM_THREADS=4 --env PYTHONUNBUFFERED=1 \
  --bind /tmp/if2453 --bind /tmp/if2453/r10:/tmp/if2453/r10:ro --bind /sys/fs/cgroup:/sys/fs/cgroup:ro \
  /stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif python -m biobox.cli "$@"
RC=$?
T1=$(date +%s.%N)
echo "w25: $UNIT wall_s $(echo "$T1 - $T0" | bc)"
echo "w25: $UNIT usage $(systemctl --user show "$UNIT" -p CPUUsageNSec -p MemoryPeak | paste -sd' ')"
echo "w25: $UNIT exited $RC"
exit $RC
