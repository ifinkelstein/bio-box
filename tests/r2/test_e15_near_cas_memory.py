"""Ring R2, card E-15: the filtered CRISPR producer in bounded memory (AGENTS §4.1, §4.2), in the form of
`tests/r2/test_crispr_memory.py`.

`scripts/stress/synth_scale.py` is run unchanged at 10 and at 100 copies of the micro-S0 fixture (about 8 and 80 Mbp).
Its stage S0 has gene calls (`feature_protein`) and a synthetic domain library; the test copies that library's rows
into a `pfam-38.2` partition with every fifth family renamed to Cas6f (PF09618), so the filter has Cas genes to find,
and runs the default producer first (in this process, not measured). The filtered producer then runs on each store in a
process of its own, launched through `runstep.measured`, with one MinCED process at a time; its peak memory is
`runstep.step_peak`: the largest `ru_maxrss` in the process tree, so the Python side (with DuckDB) or one JVM.

The bound. The producer holds one contig's sequence, one batch's parsed rows and its writers' buffers; a JVM holds
MinCED's fixed heap; the filter's tables live in one DuckDB connection with a memory cap that spills to scratch. The
100x run may exceed the 10x run by at most `GROWTH_BOUND` (32 MiB) and neither may pass `PEAK_BOUND` (768 MiB: the JVM's
heap plus 256 MiB, derived as for the default producer).
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

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import io, runstep
from biobox.registry.schemas import table

pytestmark = [pytest.mark.r2, pytest.mark.covers("E-15"),
              pytest.mark.skipif(shutil.which("minced") is None, reason="minced (MinCED) not on PATH; it is in the "
                                                                       "bio-box image, not in GitHub Actions")]

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts/stress/synth_scale.py"
SMALL, BIG = 10, 100
GROWTH_BOUND = 32 * 2**20
PEAK_BOUND = 768 * 2**20

CHILD = r"""
import json, sys
from pathlib import Path
from biobox.l1 import crispr
res = crispr.run_near_cas(Path(sys.argv[1]), "S0", 1, library_version="38.2", tmp=Path(sys.argv[2]),
                          out=Path(sys.argv[3]))
print(json.dumps({k: res[k] for k in ("contigs", "bases", "batches", "arrays", "kept")}))
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


def _as_pfam(root: Path) -> None:
    """The synthetic library's rows as `pfam-38.2`, every fifth family (by its number) renamed PF09618.1 (Cas6f)."""
    spec = table("l1_domains")
    src = io.part_files(root / "l1/domains/library=synfam-1/stage=S0")

    def rename(acc: str) -> str:
        return "PF09618.1" if acc.startswith("F") and int(acc[1:]) % 5 == 0 else acc

    def batches():
        for f in src:
            for b in pq.ParquetFile(f).iter_batches(batch_size=65_536):
                t = pa.Table.from_batches([b])
                k = t.schema.get_field_index("model_acc")
                t = t.set_column(k, "model_acc", pa.array([rename(a) for a in t.column(k).to_pylist()], pa.string()))
                yield t.cast(spec.schema)

    io.write_table(batches(), spec, root / "l1/domains/library=pfam-38.2/stage=S0/part-0.parquet")


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
    from biobox.l1 import crispr

    out = {}
    for n in (SMALL, BIG):
        root = _store(tmp_path_factory.mktemp(f"e15x{n}"), n)
        _as_pfam(root)
        crispr.run(root, "S0", 1, tmp=tmp_path_factory.mktemp(f"e15x{n}-scratch"))
        out[n] = root
    return out


def test_near_cas_peak_memory_does_not_grow_with_the_stage(stores, tmp_path, record_property):
    from biobox.l1 import crispr

    assert crispr.JAVA_HEAP == "512m"  # PEAK_BOUND is derived from it
    got = {n: _measured(root, tmp_path) for n, root in stores.items()}
    record_property("e15_peak_bytes", {n: g[0] for n, g in got.items()})
    record_property("e15_runs", {n: g[2] for n, g in got.items()})
    print(f"E15 peaks {{copies: (bytes, source, result)}}: { {n: g for n, g in got.items()} }")
    assert got[BIG][2]["bases"] >= 9 * got[SMALL][2]["bases"] > 0
    growth = got[BIG][0] - got[SMALL][0]
    assert growth <= GROWTH_BOUND, (f"peak RSS {got[SMALL][0] / 2**20:.1f} MiB at {SMALL}x, {got[BIG][0] / 2**20:.1f} "
                                    f"MiB at {BIG}x: grew {growth / 2**20:.1f} MiB, bound {GROWTH_BOUND / 2**20:.0f}")
    assert max(g[0] for g in got.values()) <= PEAK_BOUND, got
