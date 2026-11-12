"""Print the caller versions in the image and whether the data root is writable from inside it (card C-21f).

    python envcheck.py ROOT WORK

The root must be read-only (exit 1 otherwise); WORK must be writable. `os.access` asks the kernel and writes
nothing.
"""

import os
import sys

import pyrodigal
import pyrodigal_gv

import biobox
from biobox.l1 import genecall

root, work = sys.argv[1], sys.argv[2]
ro = not os.access(root, os.W_OK) and not os.access(os.path.join(root, "l1"), os.W_OK)
print(f"pyrodigal {pyrodigal.__version__} pyrodigal-gv {pyrodigal_gv.__version__} biobox {biobox.__file__}")
print(f"GENE_CALL_VERSION {genecall.GENE_CALL_VERSION}")
print(f"root read-only: {ro}; work writable: {os.access(work, os.W_OK)}")
sys.exit(0 if ro and os.access(work, os.W_OK) else 1)
