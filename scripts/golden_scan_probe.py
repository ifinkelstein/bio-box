"""Time library domain scans on an existing golden store (one-off S0 measurement).

    python scripts/golden_scan_probe.py STORE_ROOT REFDB_ROOT CPUS
"""

import json
import sys
import time
from pathlib import Path

import yaml

from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION

root, refdb_root, cpus = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
refdb = yaml.safe_load((Path(__file__).resolve().parents[1] / "refdb.yaml").read_text())
out = {"cpus": cpus}
for db in refdb["databases"]:
    if db["name"] not in ("pfam", "ncbifam"):
        continue
    path = refdb_root / db["name"] / db["version"] / db["file"]
    t0 = time.perf_counter()
    rep = domains.run(root, "S0", GENE_CALL_VERSION, db["name"], db["version"], path, db["cutoff"], cpus=cpus)
    rep["seconds"] = round(time.perf_counter() - t0, 1)
    out[db["name"]] = rep
    print(json.dumps(out), flush=True)
