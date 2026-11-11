#!/bin/sh
# Run a Python script inside the pinned bio-box image with the NEW source (bio-box-next, f9af5b6) bound in.
#   sh bbn.sh <script.py> [args...]
# OMP_NUM_THREADS caps Arrow's CPU pool (Arrow reads it); builders set their own DuckDB/MMseqs2/pyrodigal threads.
exec nice -n 19 singularity exec --cleanenv --env PYTHONPATH=/tmp/if2453/pilot_s1/pydeps --env OMP_NUM_THREADS=4 \
  --bind /stor/home/if2453/bio-box-next/src:/opt/bio-box/src --bind /tmp/if2453 \
  --bind /stor/scratch/Wilke/if2453/all-the-data \
  /stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif python "$@"
