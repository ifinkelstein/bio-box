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
scales the sort's runs to the store, as card B-29's test scaled the postings limit: runs of 4,096 rows, read
back 256 rows at a time, merged 64 at a time as in production (the 10x store's dark ORFs make 3 runs, the 200x
store's 43, all open in one merge). With those caps, the 200x store has twenty times the regions (IGRs and
dark ORFs), contigs, genes and sequence of the 10x store, and a build whose memory does not grow with them
has the same peak at both
scales, up to the caps that the larger store reaches and the smaller does not (the meta row group, about
11 MB more; the input batches) and allocator noise. So the bound is on the growth: the 200x build's peak may
exceed the 10x build's by at most `GROWTH_BOUND` (32 MiB), on any host, and it does not grow with the number
of regions. Before card W-16 the build held every region, gene and strand of the release in Python
objects: 329.6 MiB at 10x (11,631 regions) and 516.7 MiB at 200x (231,152 regions), 187.2 MiB more (measured,
wilkcomp02, 2026-11-10, this test on the tree before W-16); on the 10 Gbp release `rN.1` it filled a 32 GiB
cap (measured, wilkcomp02, 2026-11-10, planning/measurements/2026-11-10_trial_search_rN1.md). With the
streaming build: 330.4 MiB at 10x and 330.6 MiB at 200x, 0.1 MiB more (measured, wilkcomp02, 2026-11-10).

These runs are too small for one defect to show: a run reader that reads ahead. pyarrow's
`ParquetFile.iter_batches` reads a whole file ahead, so a 64-way merge of real runs held one run per reader
(1.8 GB VmHWM in the unique part on `rN.1`, measured, wilkcomp02, 2026-11-10). The second test checks that
directly: the Arrow memory of a 64-way merge of runs larger than its read batches.
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
RUN_ROWS, FAN_IN, MERGE_BATCH = 4_096, 64, 256  # the sort's caps in the measured child (module docstring)

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


# Arrow bytes a 64-way merge may reach: 64 readers x one 256-row group of about 420-byte rows is 6.9 MB; the runs
# themselves are 64 x 4,096 rows, 110 MB
MERGE_ARROW_BOUND = 16 * 2**20

MERGE_CHILD = r"""
import json, random, sys
from operator import itemgetter
from pathlib import Path
import pyarrow as pa
from biobox.search import pack
d = Path(sys.argv[1])
pack.RUN_ROWS, pack.MERGE_BATCH, pack.FAN_IN = 4_096, 256, 64
rng = random.Random(1)
rows = [(rng.randbytes(16) if i % 7 else b"\0" * 16, f"{i:08d}" + "A" * 392) for i in range(64 * 4_096)]
s = pack._Sorter(d / "runs", pack._PROTEIN_SORT, itemgetter(0), None)
for r in rows:
    s.add(r)
runs = len(list((d / "runs").iterdir()))
before = pa.default_memory_pool().max_memory()
got = list(s.rows())
s.close()
print(json.dumps({"runs": runs, "before": before, "peak": pa.default_memory_pool().max_memory(),
                  "stable": got == sorted(rows, key=itemgetter(0)), "left": (d / "runs").exists()}))
"""


def test_w16_a_64_way_merge_holds_a_row_group_per_run_and_sorts_stably(tmp_path):
    """64 runs of 4,096 rows (protein uid, 400 characters: 110 MB in all) merged at once, each read back 256 rows at
    a time, in a fresh process: the Arrow memory pool's high-water mark over the whole merge stays under
    `MERGE_ARROW_BOUND` (a reader holds its row group, not its run), and the merge equals a stable sort of the
    rows (the oracle is Python's `sorted`, which is stable; every 7th row has the same uid)."""
    p = subprocess.run([sys.executable, "-c", MERGE_CHILD, str(tmp_path)], capture_output=True, text=True,
                       timeout=600)
    assert p.returncode == 0, p.stderr[-4000:]
    got = json.loads(p.stdout.strip().splitlines()[-1])
    assert got["runs"] == 64 and got["stable"] and not got["left"], got
    assert got["peak"] <= MERGE_ARROW_BOUND, f"the merge reached {got['peak'] / 2**20:.1f} MiB of Arrow memory"
