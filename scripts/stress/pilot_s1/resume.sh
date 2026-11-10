#!/bin/sh
# Continue the S1 pilot unattended. Waits for the gene-calling process started earlier (it loses
# its stdout when its ssh session ends and stops after saving the gene-calling step), then runs the
# remaining steps. If gene calling did not finish, pilot_resume.py runs it again.
case "${1:-}" in
  -h|--help)
    echo "usage: resume.sh [PID_OF_EARLIER_PROCESS]  (pod only; continues the S1 pilot, see scripts/stress/README.md)"
    exit 0 ;;
esac
W=/tmp/if2453/pilot_s1
exec >> "$W/resume.log" 2>&1 < /dev/null
echo "resume wrapper started $(date)"
if [ -n "$1" ]; then
  tail --pid="$1" -f /dev/null
  echo "earlier process $1 gone $(date)"
fi
if [ -e /stor/scratch/Wilke/if2453/all-the-data/STOP ]; then
  echo "STOP file present; not starting"
  exit 75
fi
exec nice -n 19 singularity exec --cleanenv \
  --bind /stor/home/if2453/bio-box/src:/opt/bio-box/src --bind /tmp/if2453 \
  /stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif \
  python "$W/scripts/pilot_resume.py" "$W/src/14Apr2025_genomes.fa.gz" "$W" rS1.pilot1 --threads 16 --step-timeout 12600
