"""Ring R2, card C-38h (AGENTS §4.1, §4.2): bounded memory for domain planning with shared owned contigs.

`scripts/stress/synth_scale.py` is run unchanged at 10x (`SMALL`) and 100x (`BIG`) copies of the micro-S0 fixture.
Stage `S0` acts as an earlier `called` stage (owning all contigs in `feature_protein` and `protein_seq`, with NO
domain scan). Stage `T` has `genome_input` listing every contig owned by `S0` (and 0 own `protein_seq` files).
`domains.plan` for `T` with `seq_stages=["T", "S0"]` runs in a child process measured via `runstep.measured` and
`runstep.step_peak`: it must stream `T`'s shared-contig proteins from `S0` without holding the stage in Python
collections, staying under `GROWTH_BOUND` (64 MiB growth from 10x to 100x) and `PEAK_BOUND` (512 MiB absolute).
"""

from __future__ import annotations

import contextlib
import hashlib
import io as stdio
import json
import runpy
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import io as bio
from biobox import runstep
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.registry.schemas import table

pytestmark = [pytest.mark.r2, pytest.mark.covers("C-38h", "MR-08")]

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts/stress/synth_scale.py"
SMALL, BIG = 10, 100
GROWTH_BOUND = 64 * 2**20
PEAK_BOUND = 512 * 2**20
LIB, VER = "synfam", "1"

CHILD = r"""
import json, sys
from pathlib import Path
from biobox.l1 import domains
root, hmm, out = (Path(x) for x in sys.argv[1:4])
p = domains.plan(root, "T", sys.argv[4], "synfam", "1", hmm, "gathering", 10**12, seq_stages=["T", "S0"], out=out)
print(json.dumps(p, sort_keys=True))
"""


def _store(out: Path, copies: int) -> Path:
    argv = sys.argv
    sys.argv = [str(SCRIPT), str(REPO), str(out), str(copies), "--shard-bases", str(2**26), "--seed", "1"]
    try:
        with contextlib.redirect_stdout(stdio.StringIO()):
            g = runpy.run_path(str(SCRIPT), run_name="__main__")
    finally:
        sys.argv = argv
    g["con"].close()
    return Path(g["root"])


def _prepare_called_s0_and_shared_t(root: Path) -> int:
    """Write `protein_seq` for `S0` and `genome_input` for `T` referencing all contigs of `S0`.
    Returns the number of distinct proteins in `S0`'s `feature_protein`."""
    fp_files = bio.part_files(root / f"l1/feature_protein/gene_call_version={V}/stage=S0")
    uids_set: set[bytes] = set()
    contigs_set: set[bytes] = set()
    for f in fp_files:
        t = pq.read_table(f, columns=["contig_uid", "protein_uid"], partitioning=None)
        contigs_set.update(t.column("contig_uid").to_pylist())
        uids_set.update(t.column("protein_uid").to_pylist())
    uids = sorted(uids_set)
    seqs = ["M" * (60 + u[0] % 200) for u in uids]
    ps_table = pa.table({
        "protein_uid": uids,
        "protein_digest": [hashlib.sha256(s.encode()).digest() for s in seqs],
        "length": [len(s) for s in seqs],
        "seq": seqs,
    })
    bio.write_table(
        ps_table,
        table("protein_seq"),
        root / f"l0seq/protein_seq/gene_call_version={V}/stage=S0/protein_prefix=all/part-0.parquet",
    )
    contigs = sorted(contigs_set)
    gid = b"\x01" * 16
    gi_table = pa.table({
        "contig_uid": contigs,
        "genome_uid": [gid] * len(contigs),
        "is_circular": [False] * len(contigs),
        "kind": ["virus"] * len(contigs),
        "genetic_code": pa.array([11] * len(contigs), type=pa.uint8()),
        "assembly_quality": ["high"] * len(contigs),
    })
    bio.write_table(gi_table, table("genome_input"), root / "l0seq/genome_input/stage=T/part-0.parquet")
    return len(uids)


def _measured(root: Path, hmm: Path, tmp: Path) -> tuple[int, str, dict]:
    name = root.parent.name
    peak, out = tmp / f"peak-{name}.json", tmp / f"out-{name}"
    cmd = runstep.measured([sys.executable, "-c", CHILD, str(root), str(hmm), str(out), V], peak, False)
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    assert p.returncode == 0, p.stderr[-4000:]
    rss, source = runstep.step_peak(peak)
    assert rss is not None, source
    res = json.loads(p.stdout.strip().splitlines()[-1])
    return rss, source, res


@pytest.fixture(scope="module")
def stores(tmp_path_factory):
    out = {}
    hmm = tmp_path_factory.mktemp("c38h-hmm") / "stand-in.hmm"
    hmm.write_bytes(b"a stand-in for an HMM file: only its sha256 is read by a plan\n")
    for n in (SMALL, BIG):
        root = _store(tmp_path_factory.mktemp(f"c38hx{n}"), n)
        n_prots = _prepare_called_s0_and_shared_t(root)
        out[n] = (root, n_prots)
    return out, hmm


def test_c38h_shared_contig_domain_plan_peak_memory_is_bounded(stores, tmp_path, record_property):
    roots, hmm = stores
    got = {n: _measured(root, hmm, tmp_path) for n, (root, _) in roots.items()}
    record_property("c38h_peak_bytes", {n: g[0] for n, g in got.items()})
    record_property("c38h_runs", {n: g[2] for n, g in got.items()})
    print(f"C38h peaks {{copies: (bytes, source, result)}}: {got}")
    for n, (_, n_prots) in roots.items():
        res = got[n][2]
        assert res["proteins"] == n_prots and res["scanned_proteins"] == n_prots
        assert res["copied_proteins"] == 0 and res["shards"] >= 1
    assert got[BIG][2]["proteins"] >= 9 * got[SMALL][2]["proteins"]
    growth = got[BIG][0] - got[SMALL][0]
    assert growth <= GROWTH_BOUND, (
        f"peak RSS {got[SMALL][0] / 2**20:.1f} MiB at {SMALL}x, {got[BIG][0] / 2**20:.1f} MiB at {BIG}x: "
        f"grew {growth / 2**20:.1f} MiB, bound {GROWTH_BOUND / 2**20:.0f}"
    )
    assert max(g[0] for g in got.values()) <= PEAK_BOUND, got
