"""B-30 acceptance probe: ingest of the 1,000x synthetic stage by batch processes (D-01; ring R2).

Not collected by the suite (its name is not `test_*.py`): run it by path on the pod, about 15 minutes on one
thread,

    sh scripts/pod_pytest.sh B-30 tests/r2/probe_b30_ingest_1000x.py -s

The source objects are those of `scripts/stress/synth_scale.py`, run unchanged (by path, as `__main__`) at
COPIES = 1000 until it calls `ingest` on them: the fixture's contigs copied 1,000 times with 3% substitutions,
50 copies per object. That call is stopped there (the script's later builders are not needed). The oracle is the
single-process ingest of those objects, itself a `run-step`; the batch build is every `stage ingest` action, each
its own `run-step`, at the batch size of `test_scale_layout` (2**26 bases). Each run record's CPU-minutes must be
under the profile's budget, and the files must equal the oracle's byte for byte. The script prints one JSON line
with each step's CPU-minutes, wall seconds and peak RSS (from the run records).
"""

from __future__ import annotations

import json
import runpy
import sys
from pathlib import Path

import pytest

import biobox.l0seq
from biobox.l0seq import store

sys.path.insert(0, str(Path(__file__).parent))
import test_builders_out as c08  # noqa: E402
import test_ingest_batches as t  # noqa: E402

pytestmark = [pytest.mark.r2, pytest.mark.covers("D-01")]

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts/stress/synth_scale.py"
COPIES = 1000
SHARD_BASES = 2**26  # as test_scale_layout builds its synthetic store
PASS1_BYTES = 200 * 2**20  # bytes of FASTA per pass-1 shard: four shards of the 20 objects


class Stop(Exception):
    """Raised in place of the script's ingest of the synthetic objects."""


def synthetic_objects(out: Path, monkeypatch) -> list:
    """The 1,000x objects of synth_scale.py (its first `ingest`, of the fixture itself, runs as written)."""
    real, got = biobox.l0seq.ingest, []

    def ingest(objs, root, stage, **kw):
        if Path(root).parent.name == f"x{COPIES}":
            got.extend(objs)
            raise Stop
        return real(objs, root, stage, **kw)

    monkeypatch.setattr(biobox.l0seq, "ingest", ingest)
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), str(REPO), str(out), str(COPIES), "--shard-bases", str(SHARD_BASES)])
    with pytest.raises(Stop):
        runpy.run_path(str(SCRIPT), run_name="__main__")
    return got


def test_b30_ingest_of_the_1000x_stage_by_batch_processes(tmp_path, monkeypatch):
    objs = synthetic_objects(tmp_path / "synth", monkeypatch)
    assert len(objs) == COPIES // 50
    objs_file = t.objects_file(tmp_path / "objects", [o.path for o in objs])
    # synth_scale's objects: source "synth", release "1", the SourceObject defaults
    objs_file.write_text(json.dumps([{"path": str(o.path), "source_name": o.source_name,
                                      "source_release": o.source_release} for o in objs]))
    prof = c08.profile_at(tmp_path / "prof")
    scratch = tmp_path / "scratch"
    scratch.mkdir()

    oracle = tmp_path / "oracle"
    one = [sys.executable, "-m", "biobox.cli", "stage", "ingest", "--root", str(oracle), "--stage", "S0",
           "--objects", str(objs_file), "--tmp", str(scratch), "--shard-bases", str(SHARD_BASES)]
    assert t.step(prof, oracle, "S0", "ingest", one) == 0

    root = tmp_path / "root"
    sb = ("--shard-bases", str(SHARD_BASES))
    names = []
    for k in range(len(store.pass1_shards(objs, PASS1_BYTES))):
        names.append(f"ingest.pass1-{k:05d}")
        assert t.step(prof, root, "S0", names[-1], t.cmd("pass1", root, "S0", objs_file, scratch, "--shard", str(k),
                                                          "--pass1-bytes", str(PASS1_BYTES), *sb)) == 0
    names.append("ingest.plan")
    assert t.step(prof, root, "S0", "ingest.plan", t.cmd("plan", root, "S0", objs_file, scratch, *sb)) == 0
    for j in range(len(store.read_ingest_plan(root, "S0")["batches"])):
        names.append(f"ingest.batch-{j:04d}")
        assert t.step(prof, root, "S0", names[-1], t.cmd("batch", root, "S0", objs_file, scratch, "--batch", str(j))) == 0
    names.append("ingest")
    assert t.step(prof, root, "S0", "ingest", t.cmd("merge", root, "S0", objs_file, scratch)) == 0

    def summary(r: Path, name: str) -> dict:
        rec = t.record(r, "S0", name)
        return {"cpu_min": rec["cpu_min"], "wall_s": rec["wall_s"], "peak_rss_mb": rec["peak_rss_bytes"] / 2**20,
                "host": rec["host"]}

    report = {"copies": COPIES, "shard_bases": SHARD_BASES, "pass1_bytes": PASS1_BYTES, "objects": len(objs),
              "bases": sum(d["bases"] for d in store.read_ingest_plan(root, "S0")["objects"]),
              "single_process": summary(oracle, "ingest"), "batch_steps": {n: summary(root, n) for n in names}}
    print(json.dumps(report, sort_keys=True))
    assert sum(n.startswith("ingest.batch-") for n in names) > 1 and sum(n.startswith("ingest.pass1-") for n in names) > 1
    assert all(s["cpu_min"] < prof.cpu_min_per_process for s in report["batch_steps"].values())
    assert t.data(t.tree(root)) == t.tree(oracle)  # file for file
