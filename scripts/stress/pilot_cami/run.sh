#!/bin/sh
# Short-contig pilot (pod only): pilot_cami.py, then probe_cami.py, in the pinned image, each inside a
# memory-capped systemd scope at nice 19 (AGENTS.md section 7.2). Output goes to WORKDIR/run.log.
#   usage: run.sh WORKDIR SRC_FASTA RELEASE
# Launch it detached:  systemd-run --user --unit=pilot-cami --collect sh scripts/stress/pilot_cami/run.sh ...
# Settings from the environment (defaults in brackets): CAMI_TREE checkout whose src/ and scripts run
# [/tmp/if2453/ci/B-18]; CAMI_IMAGE [/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif];
# CAMI_PYDEPS [/tmp/if2453/pilot_s1/pydeps]; CAMI_MEM scope MemoryMax [64G]; CAMI_THREADS [8].
case "${1:-}" in
  -h|--help|"")
    echo "usage: run.sh WORKDIR SRC_FASTA RELEASE  (pod only; see scripts/stress/README.md)"
    exit 0 ;;
esac
W=$1
SRC=$2
REL=$3
TREE=${CAMI_TREE:-/tmp/if2453/ci/B-18}
IMAGE=${CAMI_IMAGE:-/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif}
PYDEPS=${CAMI_PYDEPS:-/tmp/if2453/pilot_s1/pydeps}
MEM=${CAMI_MEM:-64G}
THREADS=${CAMI_THREADS:-8}
mkdir -p "$W"
exec >> "$W/run.log" 2>&1 < /dev/null
echo "run started $(date -u +%FT%TZ) W=$W SRC=$SRC REL=$REL"
uptime
in_image() {
  systemd-run --user --scope --quiet -p MemoryMax="$MEM" -p MemorySwapMax=0 nice -n 19 \
    singularity exec --cleanenv --env PYTHONPATH="$TREE/src:$PYDEPS" --env OMP_NUM_THREADS="$THREADS" \
    --bind /tmp/if2453 --bind /sys/fs/cgroup:/sys/fs/cgroup:ro --pwd "$TREE" "$IMAGE" python "$@"
}
if [ -e /stor/scratch/Wilke/if2453/all-the-data/STOP ]; then
  echo "STOP file present; not starting"
  exit 75
fi
in_image "$TREE/scripts/stress/pilot_cami/pilot_cami.py" "$SRC" "$W" "$REL" --threads "$THREADS" --workers "$THREADS"
echo "driver exit $? $(date -u +%FT%TZ)"
if [ -e /stor/scratch/Wilke/if2453/all-the-data/STOP ]; then
  echo "STOP file present; not probing"
  exit 75
fi
in_image "$TREE/scripts/stress/pilot_cami/probe_cami.py" "$W" "$REL" --threads "$THREADS"
echo "probe exit $? $(date -u +%FT%TZ)"
# MMseqs2 temp: cluster.run removes its own TemporaryDirectory; drop the (ours) parent
if [ -d "$W/mmseqs_tmp" ]; then
  rm -rf "$W/mmseqs_tmp"
fi
uptime
echo "run finished $(date -u +%FT%TZ)"
