"""Ring R2: K-13 as a regression test on micro-S0 (data_layer_and_tests.md §2A.1, DEC-17).

The decision ran on the golden set (`scripts/k13_bakeoff.py`). This reruns the same two
layouts, Parquet 64 KB chunks x 16 per row group (the best Parquet layout then) and the
store's bgzip shard, on the micro set, so a change to the BGZF writer or reader that gave up
the margin K-13 was decided on fails here, with a `test_run` row for K-13.
"""

from __future__ import annotations

import importlib.util
import random
from pathlib import Path

import pytest
import yaml

from biobox import fasta, ids
from biobox.l0seq.bgzf import write_bgzf_fasta

pytestmark = pytest.mark.r2

REPO = Path(__file__).resolve().parents[2]
FIX = REPO / "tests" / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())


def _bakeoff():
    spec = importlib.util.spec_from_file_location("k13_bakeoff", REPO / "scripts" / "k13_bakeoff.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.covers("K-13", "LY-06")
def test_k13_bgzip_reads_fewer_bytes_than_parquet_per_slice(tmp_path):
    k13 = _bakeoff()
    contigs = sorted((ids.raw_digest(r.seq), r.seq) for o in CONTRACT["objects"] for r in fasta.read(FIX / o["file"]))
    qs = k13.queries(contigs, random.Random(20261030), 60, 20)
    del qs["contig"]  # no micro contig reaches the 1 Mbp whole-contig class
    pq_path = tmp_path / "p.parquet"
    k13.write_parquet(contigs, pq_path, 65_536, 16)
    parquet = k13.run_parquet(pq_path, contigs, qs, 65_536)
    bg_path = tmp_path / "s.fa.gz"
    write_bgzf_fasta(((rd.hex(), s) for rd, s in contigs), bg_path)
    bgzf = k13.run_bgzf(bg_path, contigs, qs)
    for label in qs:
        p, b = parquet["classes"][label], bgzf["classes"][label]
        # golden: 13-17x fewer bytes read; micro files are smaller, so the margin is smaller
        assert b["bytes_read"]["mean"] * 5 < p["bytes_read"]["mean"], (label, p["bytes_read"], b["bytes_read"])
        assert b["bytes_decoded_upper"]["mean"] < p["bytes_decoded_upper"]["mean"], label
        assert b["units"]["max"] <= 2  # a 1 kb or 20 kb slice spans at most two BGZF blocks
    assert bgzf["bytes_per_bp"] < 0.35 and parquet["bytes_per_bp"] < 0.35  # storage is a wash
