"""Time ingest, gene calling and clustering of the golden set (one-off S0 measurement).

    python scripts/golden_cluster_probe.py GOLDEN_DIR WORK_DIR THREADS
"""

import json
import sys
import time
from pathlib import Path

import yaml

from biobox import fasta, ids
from biobox.l0seq import SourceObject, ingest
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, run as call_genes
from biobox.l2.cluster import run as cluster

golden_dir, work, threads = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
repo = Path(__file__).resolve().parents[1]
spec = yaml.safe_load((repo / "golden.yaml").read_text())
lock = yaml.safe_load((repo / "golden.lock.yaml").read_text())
t0 = time.perf_counter()
ingest([SourceObject(golden_dir / lock["files"][e["id"]]["path"], "golden", "1") for e in spec["genomes"]], work, "S0")
t1 = time.perf_counter()
genomes = []
for e in spec["genomes"]:
    recs = fasta.read(golden_dir / lock["files"][e["id"]]["path"])
    uids = tuple(ids.contig_uid(r.seq) for r in recs)
    circ = tuple(u for u, r in zip(uids, e["records"]) if r["topology"] == "circular")
    genomes.append(GenomeInput(uids, e["genetic_code"], kind=e["kind"], circular=circ))
calls = call_genes(work, "S0", genomes)
t2 = time.perf_counter()
clus = cluster(work, "rS0.0-dev", GENE_CALL_VERSION, threads=threads, tmp=work)
t3 = time.perf_counter()
print(json.dumps({"ingest_s": round(t1 - t0, 1), "genecall_s": round(t2 - t1, 1), "cluster_s": round(t3 - t2, 1),
                  "threads": threads, "genes": calls, "cluster": clus}, indent=2))
