"""Ring R2: the search pack build in bounded memory (card W-16; AGENTS §4.1, §4.2).

`scripts/stress/synth_scale.py` is run unchanged (by path, as `__main__`) at 10 and at 200 copies of the
micro-S0 fixture, with the real row-group caps; the IGRs of each store are added with the repo's own
`igr.run`. The pack of each store is then built and verified in a process of its own, launched through
`runstep.measured`, and its peak memory read with `runstep.step_peak` (card C-15), so the test process's
own image is not counted.

The bound. The pack build holds, by design, at most (`biobox.search.pack`, module docstring): one external
sort's buffer (`pack.RUN_ROWS` rows), `pack.FAN_IN` runs' read batches (`pack.MERGE_BATCH` rows each), the meta
file's pending row group (at most 122,880 rows, about 12 MB), one contig's sequence and ORFs, and a batch of
each input file. A store of test size fills the sort's buffer only at the larger scale, so the measured child
scales the sort to the store, as card B-29's test scaled the postings limit: runs of 4,096 rows, merged 4 at
a time, read back 256 rows at a time (the 10x store's regions fill a run; the 200x store's take several merge
passes). With those caps, the 200x store has twenty times the regions (IGRs and dark ORFs), contigs, genes
and sequence of the 10x store, and a build whose memory does not grow with them has the same peak at both
scales, up to the caps that the larger store reaches and the smaller does not (the meta row group, about
11 MB more; the input batches) and allocator noise. So the bound is on the growth: the 200x build's peak may
exceed the 10x build's by at most `GROWTH_BOUND` (32 MiB), on any host, and it does not grow with the number
of regions. Before card W-16 the build held every region, gene and strand of the release in Python
objects: 329.6 MiB at 10x (11,631 regions) and 516.7 MiB at 200x (231,152 regions), 187.2 MiB more (measured,
wilkcomp02, 2026-11-10, this test on the tree before W-16); on the 10 Gbp release `rN.1` it filled a 32 GiB
cap (measured, wilkcomp02, 2026-11-10, planning/measurements/2026-11-10_trial_search_rN1.md). With the
streaming build: 330.2 MiB at 10x and 330.4 MiB at 200x, 0.2 MiB more (measured, wilkcomp02, 2026-11-10).
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
from biobox.l1 import igr

pytestmark = [pytest.mark.r2, pytest.mark.covers("W-16")]

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts/stress/synth_scale.py"
REL = "rS0.syn"  # synth_scale's release
SMALL, BIG = 10, 200  # copies of the fixture
GROWTH_BOUND = 32 * 2**20  # bytes the BIG build's peak may exceed the SMALL build's by
RUN_ROWS, FAN_IN, MERGE_BATCH = 4_096, 4, 256  # the sort's caps in the measured child, scaled to the store

# the measured process: build the release's pack under the scaled caps, verify it, print its target counts
CHILD = r"""
import json, sys
from pathlib import Path
from biobox.search import pack
root, packs = Path(sys.argv[1]), Path(sys.argv[2])
pack.RUN_ROWS, pack.FAN_IN, pack.MERGE_BATCH = (int(x) for x in sys.argv[3:6])
man = pack.build(root, "rS0.syn", packs=packs)
problems = pack.verify(root, "rS0.syn", packs=packs, min_shard_bytes=0)
print(json.dumps({"targets": {t: v["n_seq"] for t, v in man["targets"].items()}, "problems": problems}))
"""


def _store(out: Path, copies: int) -> Path:
    """`synth_scale.py` at `copies` into `out`, then the IGRs of its stage; the data root."""
    argv = sys.argv
    sys.argv = [str(SCRIPT), str(REPO), str(out), str(copies), "--shard-bases", str(2**26), "--seed", "1"]
    try:
        with contextlib.redirect_stdout(stdio.StringIO()):
            g = runpy.run_path(str(SCRIPT), run_name="__main__")
    finally:
        sys.argv = argv
    g["con"].close()
    root = Path(g["root"])
    igr.run(root, "S0")  # circular contigs from the stage's genome_input
    return root


def _measured_build(root: Path, packs: Path, tmp: Path) -> tuple[int, str, dict]:
    peak = tmp / f"peak-{root.parent.name}.json"
    cmd = runstep.measured([sys.executable, "-c", CHILD, str(root), str(packs), str(RUN_ROWS), str(FAN_IN),
                            str(MERGE_BATCH)], peak, False)
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    assert p.returncode == 0, p.stderr[-4000:]
    rss, source = runstep.step_peak(peak)
    assert rss is not None, source
    return rss, source, json.loads(p.stdout.strip().splitlines()[-1])


@pytest.fixture(scope="module")
def stores(tmp_path_factory):
    return {n: _store(tmp_path_factory.mktemp(f"w16x{n}"), n) for n in (SMALL, BIG)}


def test_w16_pack_build_peak_memory_does_not_grow_with_the_regions_of_the_release(stores, tmp_path, record_property):
    got = {n: _measured_build(root, tmp_path / f"packs{n}", tmp_path) for n, root in stores.items()}
    for g in got.values():
        assert g[2]["problems"] == [], g[2]
    regions = {n: g[2]["targets"]["igr"] + g[2]["targets"]["dark_orf"] for n, g in got.items()}
    record_property("w16_peak_bytes", {n: g[0] for n, g in got.items()})
    record_property("w16_regions", regions)
    print(f"W16 peaks {{copies: (bytes, source, targets)}}: { {n: (g[0], g[1], g[2]['targets']) for n, g in got.items()} }")
    # the stores differ as intended: about twenty times the regions, more than a run holds at the larger scale
    assert regions[BIG] >= 18 * regions[SMALL] > 0 and regions[BIG] > 4 * RUN_ROWS, regions
    growth = got[BIG][0] - got[SMALL][0]
    assert growth <= GROWTH_BOUND, (
        f"peak RSS of the pack build ({got[BIG][1]}): {got[SMALL][0] / 2**20:.1f} MiB at {SMALL}x "
        f"({regions[SMALL]} regions), {got[BIG][0] / 2**20:.1f} MiB at {BIG}x ({regions[BIG]} regions): grew "
        f"{growth / 2**20:.1f} MiB, bound {GROWTH_BOUND / 2**20:.0f} MiB")
