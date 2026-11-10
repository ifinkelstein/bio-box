"""Build S0 from the golden set with the repo's own builders and report what was written.

    python build_s0.py REPO GOLDEN_DIR ROOT RELEASE [--cluster THREADS] [--pfam PATH] [--ncbifam PATH] [--cpus N]
                       [--workers N]

Stage A (always): ingest, genome inputs, gene calling, IGRs, tapes.
Stage B (--cluster, needs mmseqs): clustering, tokens, token postings, contig_dim (card B-18: the
pairs step and the `slot_pairs` table are gone, AGENTS §3.5). --pfam/--ncbifam add domain scans
before tokens. Writes ROOT/manifests/release=<R>.json and prints a JSON report.
"""

import argparse
import json
import platform
import resource
import sys
import time
from collections import defaultdict
from pathlib import Path

import yaml

ap = argparse.ArgumentParser()
ap.add_argument("repo", type=Path)
ap.add_argument("golden_dir", type=Path)
ap.add_argument("root", type=Path)
ap.add_argument("release")
ap.add_argument("--cluster", type=int, default=0, help="MMseqs2 threads; 0 skips stage B")
ap.add_argument("--pfam", type=Path)
ap.add_argument("--ncbifam", type=Path)
ap.add_argument("--cpus", type=int, default=1, help="threads for the domain scans")
ap.add_argument("--workers", type=int, default=1, help="genomes called at once (output bytes do not depend on it)")
a = ap.parse_args()

from biobox import catalog, fasta, ids  # noqa: E402
from biobox.l0seq import SourceObject, ingest  # noqa: E402
from biobox.l0seq.inputs import write_genome_input  # noqa: E402
from biobox.l1 import igr  # noqa: E402
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, run as call_genes  # noqa: E402
from biobox.l3 import tapes  # noqa: E402

spec = yaml.safe_load((a.repo / "golden.yaml").read_text())
lock = yaml.safe_load((a.repo / "golden.lock.yaml").read_text())
timings, steps = {}, {}


def timed(name, fn):
    t = time.perf_counter()
    steps[name] = fn()
    timings[name] = round(time.perf_counter() - t, 2)


objs = [SourceObject(a.golden_dir / lock["files"][e["id"]]["path"], "golden", str(spec["version"])) for e in spec["genomes"]]
timed("ingest", lambda: ingest(objs, a.root, "S0"))
genomes, total_bp, n_contigs = [], 0, 0
for e in spec["genomes"]:
    recs = fasta.read(a.golden_dir / lock["files"][e["id"]]["path"])
    uids = tuple(ids.contig_uid(r.seq) for r in recs)
    total_bp += sum(len(r.seq) for r in recs)
    n_contigs += len(recs)
    circ = tuple(u for u, r in zip(uids, e["records"]) if r["topology"] == "circular")
    genomes.append(GenomeInput(uids, e["genetic_code"], kind=e["kind"], circular=circ))
circular = {u for g in genomes for u in g.circular}
timed("genome_input", lambda: write_genome_input(a.root, "S0", genomes))
timed("genecall", lambda: call_genes(a.root, "S0", genomes, workers=a.workers))
timed("igr", lambda: igr.run(a.root, "S0", circular))
timed("tapes", lambda: tapes.run(a.root, "S0", circular))
if a.pfam or a.ncbifam:
    from biobox.l1 import domains

    if a.pfam:
        timed("pfam", lambda: domains.run(a.root, "S0", GENE_CALL_VERSION, "pfam", "38.2", a.pfam, "gathering", cpus=a.cpus))
    if a.ncbifam:
        timed("ncbifam", lambda: domains.run(a.root, "S0", GENE_CALL_VERSION, "ncbifam", "20.0", a.ncbifam, "trusted", cpus=a.cpus))
if a.cluster:
    from biobox import dim
    from biobox.l2 import cluster
    from biobox.l3 import postings, tokens

    tmp = a.root.parent / "mmseqs_tmp"
    tmp.mkdir(exist_ok=True)
    timed("cluster", lambda: cluster.run(a.root, a.release, GENE_CALL_VERSION, threads=a.cluster, tmp=tmp))
    timed("tokens", lambda: tokens.run(a.root, a.release))
    timed("postings", lambda: postings.run(a.root, a.release))
    timed("contig_dim", lambda: dim.run(a.root, a.release, ["S0"]))

man = catalog.build_manifest(a.root, a.release)
catalog.write_manifest(man, a.root / f"manifests/release={a.release}.json")
per = defaultdict(lambda: {"files": 0, "rows": 0, "bytes": 0, "row_groups": 0})
for f in man["files"]:
    p = per[f["table"]]
    p["files"] += 1
    p["rows"] += f["rows"]
    p["bytes"] += f["bytes"]
    p["row_groups"] += f["row_groups"]
for f in man["artifacts"]:
    p = per["artifact:" + f["artifact"]]
    p["files"] += 1
    p["bytes"] += f["bytes"]
print(json.dumps({
    "host": platform.node(), "machine": platform.machine(), "python": sys.version.split()[0],
    "release": a.release, "manifest_hash": catalog.manifest_hash(man),
    "total_bp": total_bp, "contigs": n_contigs, "genomes": len(genomes),
    "peak_rss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 if sys.platform == "linux" else 1024 * 1024)),
    "timings_s": timings, "steps": steps, "tables": dict(sorted(per.items())),
    "file_sha256": {f["path"]: f["sha256"] for f in man["files"] + man["artifacts"]},
}, indent=1, default=str))
