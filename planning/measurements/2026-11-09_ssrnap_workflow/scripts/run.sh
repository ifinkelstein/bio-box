#!/bin/sh
# Run a command inside the pinned bio-box image (nice 19), server checkout source bound in.
# usage: sh run.sh <program> [args...]   e.g. sh run.sh python s1.py
exec nice -n 19 singularity exec --cleanenv \
  --env PYTHONPATH=/tmp/if2453/pilot_s1/pydeps \
  --bind /stor/home/if2453/bio-box/src:/opt/bio-box/src \
  --bind /tmp/if2453 \
  --bind /stor/scratch/Wilke/if2453/all-the-data \
  --bind /stor/scratch/Wilke/if2453/all-the-data/l0/refdb:/refdb:ro \
  /stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif "$@"
