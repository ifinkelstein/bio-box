"""B-23: an ingest that dies mid-stage, then is run again (D-03; ring R2).

Under `bio-box run-step` the step is a real process: `ingest` reads the earlier stages from
`BIOBOX_ROOT`, writes this stage under `BIOBOX_STAGING`, and is killed by SIGKILL from a hook at
a named point (no sleeps, as in the C-03 tests). run-step publishes nothing from a killed step and
empties staging before the next attempt, so the rerun must give the same root as an uninterrupted
ingest. The oracle is a direct, uninterrupted `ingest` into its own root, compared byte for byte.

`l0seq/store_params.json` is not a path any registry table or artifact claims, so run-step cannot
publish it (LY-10). The run-step root is given the oracle's copy beforehand, as a store whose
parameters were fixed earlier holds it; ingest then checks it and writes none (open question of
card B-23).

A direct call that dies leaves its files in the stage directories; the rerun must refuse before it
reads a source, change nothing, and name the cure (a clean stage directory).
"""

from __future__ import annotations

import hashlib
import shutil
import sys
from datetime import datetime
from pathlib import Path

import pytest
import yaml

from biobox import fasta, profile, runstep
from biobox.l0seq import SourceObject, ingest

pytestmark = [pytest.mark.r2, pytest.mark.covers("D-03")]

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
PATHS = [FIX / o["file"] for o in CONTRACT["objects"]]
SHARD_BASES = 200_000  # the three fixture genomes land in three shards
NOON_MONDAY = datetime(2026, 11, 9, 12, 0)
# every point the hook is called at, in order; a kill at the last one is a finished ingest
POINTS = ["pass1", "shard:0", "shard:1", "shard:2", "contig_occurrence", "contig_index", "seq_locator",
          "ledger", "l0_manifest"]
STAGE_DIRS = ["l0seq/seq_shards/stage=S0", "l0seq/contig_occurrence/stage=S0", "l0seq/contig_index/stage=S0",
              "l0seq/seq_locator/stage=S0", "ledger/stage=S0/step=ingest"]

STEP = r"""
import os, signal, sys
from pathlib import Path
from biobox.l0seq import SourceObject, ingest
kill_at, scratch, *paths = sys.argv[1:]

def hook(point):
    if point == kill_at:
        os.kill(os.getpid(), signal.SIGKILL)

objs = [SourceObject(Path(p), "micro", "1", "INSDC (public)", True) for p in paths]
ingest(objs, os.environ["BIOBOX_ROOT"], "S0", shard_bases=200_000, tmp=scratch,
       out=os.environ["BIOBOX_STAGING"], _hook=hook)
"""


def objects() -> list[SourceObject]:
    return [SourceObject(p, "micro", "1", "INSDC (public)", True) for p in PATHS]


def tree(root: Path) -> dict[str, str]:
    """sha256 of every file under the root except run-step's `runs/`."""
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file() and p.relative_to(root).parts[0] != "runs"}


class Crash(Exception):
    """Raised by the hook to stop a direct ingest dead at a point."""


def crash_at(point: str):
    def hook(name: str) -> None:
        if name == point:
            raise Crash(name)
    return hook


@pytest.fixture(scope="module")
def oracle(tmp_path_factory) -> Path:
    root = tmp_path_factory.mktemp("oracle") / "root"
    ingest(objects(), root, "S0", shard_bases=SHARD_BASES)
    return root


@pytest.fixture
def prof(tmp_path) -> profile.Profile:
    data = {
        "max_threads_day": 4, "max_threads_night": 6, "mem_scope": "4G",
        "tmp": str(tmp_path / "tmp"), "durable": str(tmp_path / "durable"), "bulk": str(tmp_path / "bulk"),
        "nice_monolithic": 19, "nice_fanout": 15, "day_start_hour": 8, "day_end_hour": 20,
        "cpu_min_per_process": 300, "small_file_floor_bytes": 1000,
    }
    (tmp_path / "profile.yaml").write_text(yaml.safe_dump(data))
    return profile.load(tmp_path / "profile.yaml")


def run_ingest(prof, root: Path, scratch: Path, kill_at: str) -> int:
    scratch.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-c", STEP, kill_at, str(scratch), *map(str, PATHS)]
    return runstep.run_step(prof, root, "S0", "ingest", cmd, clock=lambda: NOON_MONDAY, load=lambda: 0.0)


def test_hook_sees_every_point_in_order(tmp_path, oracle):
    seen: list[str] = []
    ingest(objects(), tmp_path / "root", "S0", shard_bases=SHARD_BASES, _hook=seen.append)
    assert seen == POINTS
    assert tree(tmp_path / "root") == tree(oracle)


@pytest.mark.parametrize("point", POINTS)
def test_killed_under_run_step_then_rerun_is_byte_identical(prof, tmp_path, oracle, point):
    root = tmp_path / "root"
    params = root / "l0seq/store_params.json"
    params.parent.mkdir(parents=True)
    shutil.copy(oracle / "l0seq/store_params.json", params)
    staging = Path(prof.tmp) / "staging" / "S0.ingest"

    assert run_ingest(prof, root, tmp_path / "scratch", point) == 1  # killed: run-step's exit code
    assert tree(root) == {"l0seq/store_params.json": tree(oracle)["l0seq/store_params.json"]}  # nothing published
    left = sorted(p.relative_to(staging).as_posix() for p in staging.rglob("*") if p.is_file())
    if point.startswith("shard:"):  # died mid-stage: some of the stage is written, not all
        assert f"l0seq/seq_shards/stage=S0/shard-000{point[-1]}.fa.gz" in left
        assert not any(p.startswith("l0seq/contig_occurrence") for p in left)

    assert run_ingest(prof, root, tmp_path / "scratch", "none") == 0
    assert tree(root) == tree(oracle)
    assert not staging.exists()


def test_uninterrupted_run_step_equals_direct_ingest(prof, tmp_path, oracle):
    root = tmp_path / "root"
    (root / "l0seq").mkdir(parents=True)
    shutil.copy(oracle / "l0seq/store_params.json", root / "l0seq/store_params.json")
    assert run_ingest(prof, root, tmp_path / "scratch", "none") == 0
    assert tree(root) == tree(oracle)


@pytest.mark.parametrize("point", POINTS[:-1])
def test_direct_rerun_over_a_crashed_stage_refuses_before_reading(tmp_path, oracle, monkeypatch, point):
    root = tmp_path / "root"
    with pytest.raises(Crash):
        ingest(objects(), root, "S0", shard_bases=SHARD_BASES, _hook=crash_at(point))
    crashed = tree(root)
    reads = []
    orig = fasta.iter_records
    monkeypatch.setattr(fasta, "iter_records", lambda p: (reads.append(p), orig(p))[1])
    if point != "pass1":  # pass 1 writes nothing to the stage; a crash there leaves it clean
        with pytest.raises(FileExistsError, match="clean stage directory"):
            ingest(objects(), root, "S0", shard_bases=SHARD_BASES)
        assert reads == []
        assert tree(root) == crashed
        for d in STAGE_DIRS:  # the cure the error names
            shutil.rmtree(root / d, ignore_errors=True)
    ingest(objects(), root, "S0", shard_bases=SHARD_BASES)
    assert tree(root) == tree(oracle)


def test_a_finished_stage_is_not_ingested_again(tmp_path, oracle, monkeypatch):
    root = tmp_path / "root"
    shutil.copytree(oracle, root)
    reads = []
    orig = fasta.iter_records
    monkeypatch.setattr(fasta, "iter_records", lambda p: (reads.append(p), orig(p))[1])
    with pytest.raises(FileExistsError, match="clean stage directory"):
        ingest(objects(), root, "S0", shard_bases=SHARD_BASES)
    with pytest.raises(FileExistsError, match="clean stage directory"):  # a staging run checks the root too
        ingest(objects(), root, "S0", shard_bases=SHARD_BASES, out=tmp_path / "staging")
    assert reads == [] and tree(root) == tree(oracle)
    assert not (tmp_path / "staging").exists()


def test_leftover_temporary_file_counts(tmp_path):
    root = tmp_path / "root"
    shard_dir = root / "l0seq/seq_shards/stage=S0"
    shard_dir.mkdir(parents=True)
    (shard_dir / ".shard-0000.fa.gz.tmp-123").write_bytes(b"half a shard from a killed writer")
    with pytest.raises(FileExistsError, match=r"\.shard-0000\.fa\.gz\.tmp-123"):
        ingest(objects(), root, "S0", shard_bases=SHARD_BASES)


def test_scratch_may_not_lie_under_the_output(tmp_path):
    out = tmp_path / "staging"
    with pytest.raises(ValueError, match="scratch"):
        ingest(objects(), tmp_path / "root", "S0", tmp=out / "scratch", out=out)
