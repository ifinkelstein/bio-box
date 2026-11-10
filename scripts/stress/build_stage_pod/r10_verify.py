"""Verify releases of a data root by their manifests (10 Gbp rebuild). Read only.

    python r10_verify.py ROOT RELEASE [RELEASE ...]

For each release: `catalog.read_manifest(catalog.manifest_path(ROOT, R))`, then
`catalog.verify_manifest(manifest)` (every listed file re-read: size and sha256; unlisted files in the
release's scope reported), as RUNBOOK §6 does. Prints one JSON line per release (stages, files, manifest
hash, problems, seconds) and exits 1 if any release has a problem.
"""

import json
import sys
import time
from pathlib import Path

from biobox import catalog

if len(sys.argv) < 3 or sys.argv[1] in ("-h", "--help"):
    print(__doc__)
    sys.exit(0 if sys.argv[1:2] in (["-h"], ["--help"]) else 2)
root, bad = Path(sys.argv[1]), 0
for rel in sys.argv[2:]:
    t = time.perf_counter()
    m = catalog.read_manifest(catalog.manifest_path(root, rel))
    problems = catalog.verify_manifest(m)
    bad += bool(problems)
    print(json.dumps({"release": m["release_id"], "stages": m["stages"], "files": len(m["files"]),
                      "hash": catalog.manifest_hash(m), "problems": problems,
                      "verify_s": round(time.perf_counter() - t, 1)}), flush=True)
sys.exit(1 if bad else 0)
