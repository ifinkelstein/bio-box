"""Re-verify the scratch copy against each of its manifests (sha256 and size of every listed file).

    python verify_scratch.py OUT_JSON RELEASE [RELEASE ...]
"""

import json
import sys
import time
from pathlib import Path

from biobox import catalog

SCR = Path("/stor/scratch/Wilke/if2453/all-the-data/provisional/rL")
LOCAL = Path("/tmp/if2453/ladder/root")
out = {}
for rel in sys.argv[2:]:
    t = time.perf_counter()
    m = catalog.read_manifest(SCR / f"manifests/release={rel}.json")
    probs = catalog.verify_manifest(m, SCR)
    lm = catalog.read_manifest(LOCAL / f"manifests/release={rel}.json")
    out[rel] = {"files": len(m["files"]), "artifacts": len(m["artifacts"]),
                "bytes": sum(f["bytes"] for f in m["files"] + m["artifacts"]),
                "missing_or_changed": [p for p in probs if not p.startswith("not in manifest")],
                "later_files_not_in_this_manifest": sum(p.startswith("not in manifest") for p in probs),
                "scratch_manifest_hash": catalog.manifest_hash(m), "local_manifest_hash": catalog.manifest_hash(lm),
                "same_file_sha256_as_local": sorted((f["path"], f["sha256"]) for f in m["files"] + m["artifacts"]) ==
                                             sorted((f["path"], f["sha256"]) for f in lm["files"] + lm["artifacts"]),
                "verify_s": round(time.perf_counter() - t, 1)}
Path(sys.argv[1]).write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
