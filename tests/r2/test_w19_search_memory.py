"""Ring R2: a search in bounded memory (card W-19; AGENTS §4.1, §4.2).

`scripts/stress/synth_scale.py` is run unchanged (by path, as `__main__`) at 10 and at 100 copies of the micro-S0
fixture, the IGRs of each store are added with the repo's own `igr.run`, and each store's pack is built (card W-16).
The synthetic stores hold no `protein_seq`, so their `unique` target is empty; their `dark_orf` target scales with
the copies, and it is what this test searches: a `phmmer` of one dark ORF of the pack against the whole target,
with that ORF as its positive control and the shuffled seed as its negative control (two passes over the target).

The search runs in a process of its own, launched through `runstep.measured` (card C-15), so the test process's
image is not counted: that process is the runner (`runner.run` with a pending record, as `search hmm` runs it),
and the runner measures its worker in turn (card W-14: the run record's `peak_rss_bytes`). The child prints its
own `ru_maxrss` (the runner alone), the peak of its Python heap over the run (`tracemalloc`) and the run record;
the helper's file gives the largest process of the tree. The runner's resident peak is its imports' (about 330
MiB at either scale, measured, wilkcomp02, 2026-11-11, this test on the tree before W-19), so a set of every
target name, built later at a lower level, does not show in it: the runner is bounded by its heap instead.

The bound. A search holds, by design (`biobox.search.engine`, `biobox.search.runner`, module docstrings): one block
of the target (`engine.BLOCK_RESIDUES` residues, digitized), its hits, and the hits of the run; the runner holds
the pack manifest and one read chunk of one shard while it looks up the positive control. A store of test size
is far smaller than a production block, so the measured child scales the block to the store, as the W-16 test
scaled its sort: blocks of 2**16 residues, so that both stores are searched in many blocks. With that, the 100x
store has ten times the dark ORFs of the 10x store, and a search whose memory does not grow with the target has
the same peak at both scales, up to allocator noise and the hits (one per copy of the seed's ORF: about ten times
as many, a few hundred bytes each). So the bound is on the growth: the 100x search's peak may exceed the 10x
search's by at most `GROWTH_BOUND` (8 MiB), in the worker's resident memory and in the runner's Python heap.
Before card W-19 the worker held the whole target three times over (a dict of `str`, a `TextSequenceBlock` and
its digitized copy) and the runner a set of every target name and a decoded copy of every shard: on the 10 Gbp
release `rN.1` 26.27 GiB and 5.65 GiB (measured, wilkcomp02, 2026-11-11,
planning/measurements/2026-11-11_trial_search_rN1_second.md); here, the worker 323.0 MiB at 10x (8,701 dark ORFs)
and 352.6 MiB at 100x (86,427), 29.6 MiB more (measured, wilkcomp02, 2026-11-11, this test on the tree before
W-19).
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
from biobox.search import fastz, pack

pytestmark = [pytest.mark.r2, pytest.mark.covers("W-19")]

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts/stress/synth_scale.py"
REL = "rS0.syn"  # synth_scale's release
SMALL, BIG = 10, 100  # copies of the fixture
GROWTH_BOUND = 8 * 2**20  # bytes the BIG search's peak may exceed the SMALL search's by, per process
BLOCK_RESIDUES = 2**16  # the target's block in the measured child (module docstring)

# the measured process: the runner, with a pending record, as `search hmm` runs it
CHILD = r"""
import json, resource, sys, tracemalloc
from pathlib import Path
from biobox.search import engine, runner
root, packs, base, recipe, block = sys.argv[1:6]
engine.BLOCK_RESIDUES = int(block)
base = Path(base)
tracemalloc.start()
row = runner.run(Path(root), recipe, "rS0.syn", limits=runner.Limits(cpus=2, wall_s=1800), work_base=base / "work",
                 packs=Path(packs), runs=base / "runs", pending=True)
heap = tracemalloc.get_traced_memory()[1]
tracemalloc.stop()
self_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print(json.dumps({"status": row["status"], "error": row["error"], "worker_peak": row["peak_rss_bytes"],
                  "runner_peak": self_rss if sys.platform == "darwin" else self_rss * 1024, "runner_heap": heap,
                  "target_size": row["target_size"], "controls": json.loads(row["controls_json"] or "null")}))
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


def _recipe(packs: Path) -> str:
    """phmmer of the first dark ORF of the pack, with itself as the positive control."""
    pdir = pack.find(packs, REL)
    first = next(s for s in pack.read_manifest(pdir)[0]["shards"] if s["target"] == "dark_orf")
    name, aa = next(fastz.read_records(pdir / first["path"]))
    return (f"recipe: w19_dark\nwritten_against: {REL}\nseed: [{{sequence: {aa}}}]\ntarget: dark_orf\n"
            f"tool: phmmer\niterations: {{max: 1, inclusion_bits: 20}}\ncontrols: {{positive: [{{region_uid: {name}}}]}}\n")


def _measured_search(root: Path, packs: Path, base: Path) -> tuple[int, str, dict]:
    peak = base / "peak.json"
    cmd = runstep.measured([sys.executable, "-c", CHILD, str(root), str(packs), str(base), _recipe(packs),
                            str(BLOCK_RESIDUES)], peak, False)
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    assert p.returncode == 0, p.stderr[-4000:]
    rss, source = runstep.step_peak(peak)
    assert rss is not None, source
    return rss, source, json.loads(p.stdout.strip().splitlines()[-1])


@pytest.fixture(scope="module")
def stores(tmp_path_factory):
    out = {}
    for n in (SMALL, BIG):
        root = _store(tmp_path_factory.mktemp(f"w19x{n}"), n)
        packs = root.parent / "packs"
        pack.build(root, REL, packs=packs)
        out[n] = (root, packs)
    return out


def test_w19_search_peak_memory_does_not_grow_with_the_target(stores, tmp_path, record_property):
    got = {}
    for n, (root, packs) in stores.items():
        base = tmp_path / f"x{n}"
        base.mkdir()
        got[n] = _measured_search(root, packs, base)
    for n, g in got.items():
        assert g[2]["status"] == "ok", (n, g[2])
        assert g[2]["controls"]["positive"][0]["recovered"], (n, g[2])
    size = {n: g[2]["target_size"] for n, g in got.items()}
    worker = {n: g[2]["worker_peak"] for n, g in got.items()}
    runner_rss = {n: g[2]["runner_peak"] for n, g in got.items()}
    runner_heap = {n: g[2]["runner_heap"] for n, g in got.items()}
    tree = {n: g[0] for n, g in got.items()}
    record_property("w19_target_size", size)
    record_property("w19_worker_peak_bytes", worker)
    record_property("w19_runner_peak_bytes", runner_rss)
    record_property("w19_runner_heap_peak_bytes", runner_heap)
    record_property("w19_tree_peak_bytes", tree)
    print(f"W19 {{copies: target size, worker RSS, runner RSS, runner heap, tree RSS (source)}}: "
          f"{ {n: (size[n], worker[n], runner_rss[n], runner_heap[n], tree[n], got[n][1]) for n in got} }")
    # the targets differ as intended: about ten times the dark ORFs, many blocks at either scale
    assert size[BIG] >= 8 * size[SMALL] > 0, size
    for what, peaks in (("worker's peak RSS", worker), ("runner's peak Python heap", runner_heap)):
        growth = peaks[BIG] - peaks[SMALL]
        assert growth <= GROWTH_BOUND, (
            f"the search {what}: {peaks[SMALL] / 2**20:.1f} MiB at {SMALL}x ({size[SMALL]} targets), "
            f"{peaks[BIG] / 2**20:.1f} MiB at {BIG}x ({size[BIG]} targets): grew {growth / 2**20:.1f} MiB, "
            f"bound {GROWTH_BOUND / 2**20:.0f} MiB")
