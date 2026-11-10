#!/bin/sh
# Continue the S1 pilot unattended (pod only). Waits for an earlier pilot process if a PID is given,
# refuses to start if the STOP file exists, then runs pilot_resume.py in the image inside a memory-capped
# systemd scope at nice 19 (AGENTS.md section 7.2; card B-18 added the scope). Steps that already
# finished are skipped. Launch it detached:
#   systemd-run --user --unit=pilot-s1-resume --collect sh scripts/stress/pilot_s1/resume.sh [PID]
# Settings come from the environment; the defaults are the 2026-11-09 run's paths:
#   PILOT_WORK   work directory (root/, result.json)          /tmp/if2453/pilot_s1
#   PILOT_TREE   checkout whose src/ and scripts are used       /tmp/if2453/ci/B-18
#   PILOT_SRC    INPHARED FASTA                                 $PILOT_WORK/src/14Apr2025_genomes.fa.gz
#   PILOT_REL    release name                                   rS1.pilot1
#   PILOT_IMAGE  test image                                     /stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif
#   PILOT_PYDEPS extra Python packages                          /tmp/if2453/pilot_s1/pydeps
#   PILOT_MEM    MemoryMax of the scope                         64G
#   PILOT_THREADS MMseqs2 threads and gene-calling workers      16
case "${1:-}" in
  -h|--help)
    echo "usage: resume.sh [PID_OF_EARLIER_PROCESS]  (pod only; continues the S1 pilot, see scripts/stress/README.md)"
    exit 0 ;;
esac
W=${PILOT_WORK:-/tmp/if2453/pilot_s1}
TREE=${PILOT_TREE:-/tmp/if2453/ci/B-18}
SRC=${PILOT_SRC:-$W/src/14Apr2025_genomes.fa.gz}
REL=${PILOT_REL:-rS1.pilot1}
IMAGE=${PILOT_IMAGE:-/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif}
PYDEPS=${PILOT_PYDEPS:-/tmp/if2453/pilot_s1/pydeps}
MEM=${PILOT_MEM:-64G}
THREADS=${PILOT_THREADS:-16}
exec >> "$W/resume.log" 2>&1 < /dev/null
echo "resume wrapper started $(date)"
if [ -n "${1:-}" ]; then
  tail --pid="$1" -f /dev/null
  echo "earlier process $1 gone $(date)"
fi
if [ -e /stor/scratch/Wilke/if2453/all-the-data/STOP ]; then
  echo "STOP file present; not starting"
  exit 75
fi
exec systemd-run --user --scope --quiet -p MemoryMax="$MEM" -p MemorySwapMax=0 nice -n 19 \
  singularity exec --cleanenv --env PYTHONPATH="$TREE/src:$PYDEPS" --bind /tmp/if2453 \
  --bind /sys/fs/cgroup:/sys/fs/cgroup:ro --pwd "$TREE" "$IMAGE" \
  python "$TREE/scripts/stress/pilot_s1/pilot_resume.py" "$SRC" "$W" "$REL" --threads "$THREADS" \
  --workers "$THREADS" --step-timeout 12600
