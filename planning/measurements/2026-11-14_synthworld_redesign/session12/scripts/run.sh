#!/bin/sh
# run.sh UNIT TREE LOG command...   (card F-02 session 12: a detached unit at nice 19 with a memory cap, in the image)
UNIT=$1
TREE=$2
LOG=$3
shift 3
IMAGE=/stor/home/if2453/bio-box/images/bio-box-test_faec1f4.sif
mkdir -p /tmp/if2453/f0212/logs
if [ -f /tmp/if2453/f0212/STOP ]; then
  echo "STOP file present; not starting $UNIT"
  exit 0
fi
exec systemd-run --user --unit="$UNIT" -p MemoryMax=32G -p MemorySwapMax=0 \
  -p StandardOutput=file:"$LOG" -p StandardError=file:"$LOG.err" \
  nice -n 19 singularity exec --cleanenv --env PYTHONPATH="$TREE/src" --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 --bind /tmp/if2453 --pwd "$TREE" "$IMAGE" "$@"
