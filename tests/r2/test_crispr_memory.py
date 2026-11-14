"""Ring R2: the CRISPR producer in bounded memory (card E-12; AGENTS §4.1, §4.2), in the form of card W-16's test.

`scripts/stress/synth_scale.py` is run unchanged (by path, as `__main__`) at 10 and at 100 copies of the micro-S0
fixture (3% substitutions per copy, so 100x is about 80 Mbp in about 300 contigs). The producer then runs on each
store's stage in a process of its own, launched through `runstep.measured`, with one MinCED process at a time, and
its peak memory is read with `runstep.step_peak` (card C-15): the largest `ru_maxrss` in the process tree, so the
Python producer or one JVM, whichever is larger.

The bound. The producer holds one contig's sequence, one batch's parsed rows and its writers' buffers; a JVM holds
MinCED's fixed heap (`crispr.JAVA_HEAP`) and one record. So the peak does not grow with the stage: the 100x run may
exceed the 10x run by at most `GROWTH_BOUND` (32 MiB), and neither may pass `PEAK_BOUND`, the card's stated peak:
the JVM's heap plus 256 MiB for its own code and the Python side (derived: 512 MiB + 256 MiB = 768 MiB).
"""

from __future__ import annotations

import contextlib
import io as stdio
import json
import runpy
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from biobox import runstep

pytestmark = [pytest.mark.r2, pytest.mark.covers("E-12"),
              pytest.mark.skipif(shutil.which("minced") is None, reason="minced (MinCED) not on PATH; it is in the "
                                                                       "bio-box image, not in GitHub Actions")]

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts/stress/synth_scale.py"
SMALL, BIG = 10, 100  # copies of the fixture
GROWTH_BOUND = 32 * 2**20
PEAK_BOUND = 768 * 2**20

CHILD = r"""
import json, sys
from pathlib import Path
from biobox.l1 import crispr
res = crispr.run(Path(sys.argv[1]), "S0", 1, tmp=Path(sys.argv[2]), out=Path(sys.argv[3]))
print(json.dumps({k: res[k] for k in ("contigs", "bases", "batches", "arrays", "repeats", "spacers")}))
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


def _measured(root: Path, tmp: Path) -> tuple[int, str, dict]:
    peak = tmp / f"peak-{root.parent.name}.json"
    scratch, out = tmp / f"scratch-{root.parent.name}", tmp / f"out-{root.parent.name}"
    scratch.mkdir()
    cmd = runstep.measured([sys.executable, "-c", CHILD, str(root), str(scratch), str(out)], peak, False)
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    assert p.returncode == 0, p.stderr[-4000:]
    rss, source = runstep.step_peak(peak)
    assert rss is not None, source
    return rss, source, json.loads(p.stdout.strip().splitlines()[-1])


@pytest.fixture(scope="module")
def stores(tmp_path_factory):
    return {n: _store(tmp_path_factory.mktemp(f"e12x{n}"), n) for n in (SMALL, BIG)}


def test_crispr_peak_memory_does_not_grow_with_the_stage(stores, tmp_path, record_property):
    from biobox.l1 import crispr

    assert crispr.JAVA_HEAP == "512m"  # PEAK_BOUND is derived from it
    got = {n: _measured(root, tmp_path) for n, root in stores.items()}
    record_property("e12_peak_bytes", {n: g[0] for n, g in got.items()})
    record_property("e12_runs", {n: g[2] for n, g in got.items()})
    print(f"E12 peaks {{copies: (bytes, source, result)}}: { {n: g for n, g in got.items()} }")
    assert got[BIG][2]["bases"] >= 9 * got[SMALL][2]["bases"] > 0
    growth = got[BIG][0] - got[SMALL][0]
    assert growth <= GROWTH_BOUND, (f"peak RSS {got[SMALL][0] / 2**20:.1f} MiB at {SMALL}x, {got[BIG][0] / 2**20:.1f} "
                                    f"MiB at {BIG}x: grew {growth / 2**20:.1f} MiB, bound {GROWTH_BOUND / 2**20:.0f}")
    assert max(g[0] for g in got.values()) <= PEAK_BOUND, got
