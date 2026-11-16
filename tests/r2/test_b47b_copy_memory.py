"""Ring R2, card B-47b: a scan that copies measured rows plans and merges in bounded memory (AGENTS §4.1, §4.2), in
the form of `tests/r2/test_e15_near_cas_memory.py`.

`scripts/stress/synth_scale.py` is run unchanged at 10 and at 100 copies of the micro-S0 fixture. Its stage S0 has gene
calls (`feature_protein`, synthetic protein uids) and the rows of a synthetic library `synfam-1`, but no `protein_seq`;
the test writes one row per distinct protein of `feature_protein` (a sequence of its length in `M`), and makes S0 a
source as a completed sharded scan leaves one: a plan (`domains.plan`, with a stand-in HMM file), a ledger with this
process's tool versions, and a run record listing the domain files and the ledger with their sha256. Stage T holds
the same proteins. T's plan (naming T and S0) and merge then run in a process of their own, launched through
`runstep.measured`: every protein of T is copied from S0, no shard is scanned, and T's domain file must equal S0's
byte for byte (the same rows through the same writer). The peak is `runstep.step_peak`.

The bound. The plan and the merge read the stage's and the source's `protein_seq` in batches of `_READ_ROWS` rows,
and the merge holds one row group (122,880 rows) of the source's domain file and of the shards' runs at a time, with
its uids as a Python list, plus the writer's row group: never a set of the stage's proteins. The 100x stores are about
one row group of each table, so the absolute bound matters as much as the growth: the 100x run may exceed the 10x run
by at most `GROWTH_BOUND` (64 MiB) and neither may pass `PEAK_BOUND` (512 MiB). Measured peaks are printed and
recorded as test properties.
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

from biobox import ids, ledger, runstep
from biobox import io as bio
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.registry.schemas import table

pytestmark = [pytest.mark.r2, pytest.mark.covers("MR-08")]

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts/stress/synth_scale.py"
SMALL, BIG = 10, 100
GROWTH_BOUND = 64 * 2**20
PEAK_BOUND = 512 * 2**20
LIB, VER = "synfam", "1"
STEP = f"domains-{LIB}-{VER}"

CHILD = r"""
import json, sys
from pathlib import Path
from biobox.l1 import domains
root, hmm, out, tmp = (Path(x) for x in sys.argv[1:5])
p = domains.plan(root, "T", sys.argv[5], "synfam", "1", hmm, "gathering", 10**12, seq_stages=["T", "S0"], out=out)
m = domains.merge(root, "T", "synfam", "1", tmp=tmp, out=out)
print(json.dumps({"plan": p, "merge": m}, sort_keys=True))
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


def _protein_seq(root: Path, stage: str) -> None:
    """One protein_seq row per distinct protein of S0's feature_protein, in uid order, a sequence of `M`s."""
    fp = bio.part_files(root / f"l1/feature_protein/gene_call_version={V}/stage=S0")
    uids = sorted({u for f in fp for u in pq.read_table(f, columns=["protein_uid"], partitioning=None)
                   .column(0).to_pylist()})
    seqs = ["M" * (60 + u[0] % 200) for u in uids]
    t = pa.table({"protein_uid": uids, "protein_digest": [hashlib.sha256(s.encode()).digest() for s in seqs],
                  "length": [len(s) for s in seqs], "seq": seqs})
    bio.write_table(t, table("protein_seq"),
                    root / f"l0seq/protein_seq/gene_call_version={V}/stage={stage}/protein_prefix=all/part-0.parquet")


def _as_completed_source(root: Path, hmm: Path) -> None:
    """S0 as a completed sharded scan of `synfam-1` leaves it: plan, ledger, run record (module docstring)."""
    from biobox.l1 import domains

    # the plan is made in a root that holds only S0's proteins (a plan refuses a stage whose domain files exist),
    # then put where the scan's plan step would have published it: its paths are relative to the root
    side = root.parent / "plan-side"
    seq = f"l0seq/protein_seq/gene_call_version={V}/stage=S0"
    (side / seq).parent.mkdir(parents=True)
    (side / seq).symlink_to(root / seq)
    domains.plan(side, "S0", V, LIB, VER, hmm, "gathering", 10**12)
    rel = domains.scan_dir("S0", LIB, VER) + "/plan.json"
    (root / rel).parent.mkdir(parents=True)
    (root / rel).write_bytes((side / rel).read_bytes())
    files =bio.part_files(root / f"l1/domains/library={LIB}-{VER}/stage=S0")
    n_in = pq.read_metadata(bio.part_files(root / f"l0seq/protein_seq/gene_call_version={V}/stage=S0/protein_prefix=all")[0]).num_rows
    hit = len({u for f in files for u in pq.read_table(f, columns=["protein_uid"], partitioning=None).column(0).to_pylist()})
    tools = domains._tools()
    lib = f"{LIB}-{VER}".encode()
    ledger.write(root, "S0", STEP, [ledger.row(lib, "library", n_in, hit, n_rejected=n_in - hit, reason="no_hit",
                                               tools=tools), ledger.row(lib, "library", n_in, hit, tools=tools)])
    led = bio.part_files(root / f"ledger/stage=S0/step={STEP}")
    outputs = [{"path": p.relative_to(root).as_posix(), "bytes": p.stat().st_size, "sha256": ids.sha256_file(p).hex()}
               for p in files + led]
    rec = root / f"runs/stage=S0/{STEP}.json"
    rec.parent.mkdir(parents=True, exist_ok=True)
    rec.write_text(json.dumps({"stage": "S0", "step": STEP, "outputs": sorted(outputs, key=lambda o: o["path"])},
                              sort_keys=True, indent=1) + "\n")


def _measured(root: Path, hmm: Path, tmp: Path) -> tuple[int, str, dict]:
    name = root.parent.name
    peak, out, scratch = tmp / f"peak-{name}.json", tmp / f"out-{name}", tmp / f"scratch-{name}"
    scratch.mkdir()
    cmd = runstep.measured([sys.executable, "-c", CHILD, str(root), str(hmm), str(out), str(scratch), V], peak, False)
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    assert p.returncode == 0, p.stderr[-4000:]
    rss, source = runstep.step_peak(peak)
    assert rss is not None, source
    res = json.loads(p.stdout.strip().splitlines()[-1])
    src = bio.part_files(root / f"l1/domains/library={LIB}-{VER}/stage=S0")
    got = bio.part_files(out / f"l1/domains/library={LIB}-{VER}/stage=T")
    assert [f.read_bytes() for f in got] == [f.read_bytes() for f in src]  # every row copied, the same bytes
    return rss, source, res


@pytest.fixture(scope="module")
def stores(tmp_path_factory):
    out = {}
    hmm = tmp_path_factory.mktemp("b47b-hmm") / "stand-in.hmm"
    hmm.write_bytes(b"a stand-in for an HMM file: only its sha256 is read by a plan\n")
    for n in (SMALL, BIG):
        root = _store(tmp_path_factory.mktemp(f"b47bx{n}"), n)
        _protein_seq(root, "S0")
        _protein_seq(root, "T")
        _as_completed_source(root, hmm)
        out[n] = root
    return out, hmm


def test_b47b_copy_plan_and_merge_peak_memory_does_not_grow_with_the_stage(stores, tmp_path, record_property):
    roots, hmm = stores
    got = {n: _measured(root, hmm, tmp_path) for n, root in roots.items()}
    record_property("b47b_peak_bytes", {n: g[0] for n, g in got.items()})
    record_property("b47b_runs", {n: g[2] for n, g in got.items()})
    print(f"B47b peaks {{copies: (bytes, source, result)}}: { {n: g for n, g in got.items()} }")
    for n, (_, _, res) in got.items():
        assert res["plan"]["shards"] == 0 and res["plan"]["copied_proteins"] == res["plan"]["proteins"] > 0, n
    assert got[BIG][2]["plan"]["proteins"] >= 9 * got[SMALL][2]["plan"]["proteins"]
    growth = got[BIG][0] - got[SMALL][0]
    assert growth <= GROWTH_BOUND, (f"peak RSS {got[SMALL][0] / 2**20:.1f} MiB at {SMALL}x, {got[BIG][0] / 2**20:.1f} "
                                    f"MiB at {BIG}x: grew {growth / 2**20:.1f} MiB, bound {GROWTH_BOUND / 2**20:.0f}")
    assert max(g[0] for g in got.values()) <= PEAK_BOUND, got
