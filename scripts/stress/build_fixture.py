"""Build the micro-S0 fixture store (stage A) and report its manifest hash: the cross-host check.

    python build_fixture.py REPO ROOT RELEASE > report.json

Card B-18: the laptop no longer computes, so cross-host equality compares a GitHub Actions runner
(`.github/workflows/cross_host.yml`, set up as `ci.yml` is) with the pod's image. Both run this script
on `REPO/tests/fixtures/micro_s0` (three genomes, about 0.8 Mbp, no network, no MMseqs2): ingest,
genome inputs, gene calling, IGRs, tapes and `contig_dim`, then `catalog.build_manifest`. ROOT must not
exist yet. The report has the shape of `build_s0.py`'s, so `diff_reports.py A.json B.json` compares two
runs file by file; it also names the versions of the packages whose bytes reach the files.
"""

import argparse
import json
import platform
import resource
import sys
import time
from collections import defaultdict
from importlib.metadata import version
from pathlib import Path

import yaml

ap = argparse.ArgumentParser()
ap.add_argument("repo", type=Path)
ap.add_argument("root", type=Path)
ap.add_argument("release")
a = ap.parse_args()
if a.root.exists():
    sys.exit(f"{a.root}: exists; the fixture store is built into a fresh directory")

from biobox import catalog, dim, fasta, ids  # noqa: E402
from biobox.l0seq import SourceObject, ingest  # noqa: E402
from biobox.l0seq.inputs import write_genome_input  # noqa: E402
from biobox.l1 import igr  # noqa: E402
from biobox.l1.genecall import GenomeInput, run as call_genes  # noqa: E402
from biobox.l3 import tapes  # noqa: E402

fx = a.repo / "tests/fixtures/micro_s0"
src = yaml.safe_load((fx / "sources.yaml").read_text())
timings, steps = {}, {}


def timed(name, fn):
    t = time.perf_counter()
    steps[name] = fn()
    timings[name] = round(time.perf_counter() - t, 2)


objs = [SourceObject(fx / o["file"], src["source_name"], src["source_release"], src["license"], src["redistributable"])
        for o in src["objects"]]
timed("ingest", lambda: ingest(objs, a.root, "S0"))
genomes, total_bp = [], 0
for o in src["objects"]:
    seq = fasta.read(fx / o["file"])[0].seq
    cu = ids.contig_uid(seq)
    total_bp += len(seq)
    genomes.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"],
                               circular=(cu,) if o["topology"] == "circular" else ()))
timed("genome_input", lambda: write_genome_input(a.root, "S0", genomes))
timed("genecall", lambda: call_genes(a.root, "S0"))
timed("igr", lambda: igr.run(a.root, "S0"))
timed("tapes", lambda: tapes.run(a.root, "S0"))
timed("contig_dim", lambda: dim.run(a.root, a.release, ["S0"]))

man = catalog.build_manifest(a.root, a.release)
per = defaultdict(lambda: {"files": 0, "rows": 0, "bytes": 0, "row_groups": 0})
for f in man["files"]:
    p = per[f["table"]]
    for k in ("rows", "bytes", "row_groups"):
        p[k] += f[k]
    p["files"] += 1
for f in man["artifacts"]:
    p = per["artifact:" + f["artifact"]]
    p["files"] += 1
    p["bytes"] += f["bytes"]
print(json.dumps({
    "host": platform.node(), "machine": platform.machine(), "python": sys.version.split()[0],
    "packages": {p: version(p) for p in ("pyarrow", "duckdb", "pyrodigal", "pyrodigal-gv", "blake3", "numpy")},
    "release": a.release, "manifest_hash": catalog.manifest_hash(man),
    "total_bp": total_bp, "contigs": len(genomes), "genomes": len(genomes),
    "peak_rss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 if sys.platform == "linux" else 1024 * 1024)),
    "timings_s": timings, "steps": steps, "tables": dict(sorted(per.items())),
    "file_sha256": {f["path"]: f["sha256"] for f in man["files"] + man["artifacts"]},
}, indent=1, default=str))
