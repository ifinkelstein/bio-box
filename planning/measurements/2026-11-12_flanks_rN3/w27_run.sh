#!/bin/bash
# Card W-27 real run: one bio-box command in the image, inside the systemd unit named by $1 (started by the caller:
#   systemd-run --user --unit w27-<name> -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% -p RemainAfterExit=yes
#     -p StandardOutput=file:... nice -n 19 bash /tmp/if2453/w27/run.sh w27-<name> <bio-box args...>)
# The tree is task/W-27 mirrored to /tmp/if2453/w27/tree; matplotlib and its new dependencies are the uv.lock wheels
# unpacked into /tmp/if2453/pydeps_w27 (sha256 checked, fetch_wheels.sh). Reads the root read-only; writes only under
# /tmp/if2453/w27. Prints wall time and the unit's usage at the end.
UNIT=$1
shift
if [ -e /tmp/if2453/r10/root/STOP ]; then echo "w27: STOP file present, not running"; exit 0; fi
mkdir -p /tmp/if2453/w27/tmp
T0=$(date +%s.%N)
/usr/bin/time -v singularity exec --cleanenv \
  --env PYTHONPATH=/tmp/if2453/w27/tree/src:/tmp/if2453/pilot_s1/pydeps:/tmp/if2453/pydeps_w27 \
  --env TMPDIR=/tmp/if2453/w27/tmp --env OMP_NUM_THREADS=2 --env PYTHONUNBUFFERED=1 \
  --bind /tmp/if2453 --bind /tmp/if2453/r10:/tmp/if2453/r10:ro --bind /sys/fs/cgroup:/sys/fs/cgroup:ro \
  /stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif python -m biobox.cli "$@"
RC=$?
T1=$(date +%s.%N)
echo "w27: $UNIT wall_s $(echo "$T1 - $T0" | bc)"
echo "w27: $UNIT usage $(systemctl --user show "$UNIT" -p CPUUsageNSec -p MemoryPeak | paste -sd' ')"
echo "w27: $UNIT exited $RC"
exit $RC
