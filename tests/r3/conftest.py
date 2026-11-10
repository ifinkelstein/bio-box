"""Ring R3 shared fixtures: the golden store is built once per session and used by every R3 module."""

from __future__ import annotations

import os
import shutil
import time
from pathlib import Path

import pytest
import yaml

from biobox import fasta, ids
from biobox.l0seq import SourceObject, ingest
from biobox.l1 import igr
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, run
from biobox.l3 import postings, tapes, tokens

ROOT = Path(__file__).resolve().parents[2]


def golden_input(golden_dir: Path, lock: dict, entry: dict, code=None, kind=None) -> GenomeInput:
    recs = fasta.read(golden_dir / lock["files"][entry["id"]]["path"])
    uids = tuple(ids.contig_uid(r.seq) for r in recs)
    circ = tuple(u for u, spec in zip(uids, entry["records"]) if spec["topology"] == "circular")
    return GenomeInput(uids, code or entry["genetic_code"], kind=kind or entry["kind"], circular=circ)


@pytest.fixture(scope="session")
def store_root(tmp_path_factory):
    golden = yaml.safe_load((ROOT / "golden.yaml").read_text())
    lock_path = ROOT / "golden.lock.yaml"
    directory = os.environ.get("BIOBOX_GOLDEN_DIR")
    if not directory or not lock_path.exists():
        pytest.skip("golden set not fetched (set BIOBOX_GOLDEN_DIR)")
    directory = Path(directory)
    lock = yaml.safe_load(lock_path.read_text())
    root = tmp_path_factory.mktemp("golden")
    objs = [SourceObject(directory / lock["files"][e["id"]]["path"], "golden", str(golden["version"]))
            for e in golden["genomes"]]
    t0 = time.perf_counter()
    ingest(objs, root, "S0")
    t1 = time.perf_counter()
    inputs = [golden_input(directory, lock, e) for e in golden["genomes"]]
    report = run(root, "S0", inputs)
    t2 = time.perf_counter()
    igrs = igr.run(root, "S0", {u for gi in inputs for u in gi.circular})
    tp = tapes.run(root, "S0", {u for gi in inputs for u in gi.circular})
    print(f"\ngolden: ingest {t1 - t0:.1f}s, gene calling {t2 - t1:.1f}s, {report}, igr {igrs}, tapes {tp}")
    if shutil.which("mmseqs"):
        from biobox.l2 import cluster
        cluster.run(root, "rS0.golden", GENE_CALL_VERSION, threads=4, tmp=root.parent)
        tokens.run(root, "rS0.golden")
        print("postings", postings.run(root, "rS0.golden"))
    return root
