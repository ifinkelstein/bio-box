"""`pilot.py --resume`: the S1 pilot, skipping the steps an earlier process finished without error.

    python pilot_resume.py SRC_FASTA WORKDIR RELEASE [--limit N] [--threads 16] [--workers 1] [--step-timeout 10800]

Same arguments as `pilot.py`. Until card B-18 this file was a copy of `pilot.py` with the resume logic
added; it is now that script run with `--resume`, so the two cannot drift apart.
"""

import runpy
import sys
from pathlib import Path

if "--resume" not in sys.argv:
    sys.argv.append("--resume")
sys.argv[0] = str(Path(__file__).with_name("pilot.py"))
runpy.run_path(sys.argv[0], run_name="__main__")
