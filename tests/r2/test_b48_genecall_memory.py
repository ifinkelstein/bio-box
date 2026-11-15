"""Ring R2: gene calling in batches in bounded memory (card B-48; AGENTS §4.1, §4.2), in the form of card W-16's test.

`scripts/stress/synth_scale.py` is run unchanged (by path, as `__main__`) at 10 and at 100 copies of the micro-S0
fixture (3% substitutions per copy: 30 and 300 one-contig genomes, about 8 and 80 Mbp). Its stage S0 is then
gene-called in batches under a version of its own (`b48-mem`: the script wrote gene calls under the real version), the
10x stage as one batch and the 100x stage as ten, so that every batch holds about the 10x stage's genomes. Each batch
and the merge run as `python -m biobox.cli stage genecall batch|merge` in a process of its own, launched through
`runstep.measured`, and their peak memory is read with `runstep.step_peak` (card C-15: the largest `ru_maxrss` in the
process tree).

The bound. A batch holds its genomes' sequences one genome (or one task of short genomes) at a time, the rows of one
task, and the scratch writers' buffers (`SCRATCH_FLUSH_ROWS` rows); it streams `genome_input` and holds only the
stage's shared contigs and its own skips. So a batch's peak follows its own size, not the stage's: a batch of the 100x
stage may exceed the 10x stage's one batch by at most `GROWTH_BOUND` (32 MiB). The merge sorts in DuckDB under
`DUCK_MEMORY_LIMIT` (2 GB, spilling past it) and holds the skips: at 100x it may exceed the 10x merge by at most
`GROWTH_BOUND` too, and no process may pass `PEAK_BOUND`, the card's stated peak: DuckDB's limit and 512 MiB for the
interpreter, pyrodigal and the writers (derived: 2,048 + 512 = 2,560 MiB). The merge's run records the step's output
(the files) for the 100x stage, which the byte-identity tests check elsewhere; here only its counts.
"""

from __future__ import annotations

import contextlib
import io as stdio
import json
import runpy
import subprocess
import sys
from pathlib import Path

import pytest

from biobox import runstep
from biobox.l1 import genecall

pytestmark = [pytest.mark.r2, pytest.mark.covers("B-48")]

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts/stress/synth_scale.py"
SMALL, BIG = 10, 100  # copies of the fixture
GROWTH_BOUND = 32 * 2**20
PEAK_BOUND = 2560 * 2**20
VERSION = "b48-mem"


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


@pytest.fixture(scope="module")
def stores(tmp_path_factory):
    return {n: _store(tmp_path_factory.mktemp(f"b48x{n}"), n) for n in (SMALL, BIG)}


def _measured(args: list[str], peak: Path) -> tuple[int, dict]:
    cmd = runstep.measured([sys.executable, "-m", "biobox.cli", *args], peak, False)
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    assert p.returncode == 0, p.stderr[-4000:]
    rss, source = runstep.step_peak(peak)
    assert rss is not None, source
    return rss, json.loads(p.stdout.strip().splitlines()[-1])["result"]


def test_b48_a_batchs_and_the_merges_peak_memory_do_not_grow_with_the_stage(stores, tmp_path, record_property):
    got: dict[int, dict] = {}
    for copies, n in ((SMALL, 1), (BIG, 10)):
        root = stores[copies]
        plan = genecall.plan(root, "S0", version=VERSION, batch_cpu_min=240.0, batches=n, tmp=tmp_path)
        assert plan["genomes"] == 3 * copies
        path = genecall.write_plan(plan, genecall.work_dir(tmp_path / f"w{copies}", root, "S0", plan))
        peaks = []
        for j in range(n):
            rss, rec = _measured(["stage", "genecall", "batch", "--root", str(root), "--stage", "S0", "--gene-call-version",
                                  VERSION, "--plan", str(path), "--batch", str(j), "--workers", "2"],
                                 tmp_path / f"peak-{copies}-{j}.json")
            assert rec["genomes"] == plan["batches"][j]["genomes"]
            peaks.append(rss)
        rss, res = _measured(["stage", "genecall", "merge", "--root", str(root), "--stage", "S0", "--gene-call-version",
                              VERSION, "--plan", str(path), "--out", str(tmp_path / f"out{copies}"), "--tmp",
                              str(tmp_path)], tmp_path / f"peak-{copies}-merge.json")
        got[copies] = {"batches": peaks, "merge": rss, "genes": res["genes"]}
    record_property("b48_peaks", got)
    print(f"B48 peaks (bytes): {got}")
    assert got[BIG]["genes"] >= 9 * got[SMALL]["genes"] > 0
    small_batch = got[SMALL]["batches"][0]
    for j, rss in enumerate(got[BIG]["batches"]):
        assert rss - small_batch <= GROWTH_BOUND, (j, rss / 2**20, small_batch / 2**20)
    assert got[BIG]["merge"] - got[SMALL]["merge"] <= GROWTH_BOUND, (got[BIG]["merge"] / 2**20,
                                                                     got[SMALL]["merge"] / 2**20)
    assert max(max(g["batches"] + [g["merge"]]) for g in got.values()) <= PEAK_BOUND, got
