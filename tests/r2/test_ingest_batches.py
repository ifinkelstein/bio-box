"""B-30: a stage's ingest as pass-1 shards, a plan, one process per batch and a merge (D-01, D-03; ring R2).

The oracle is the single-process `ingest` (cards B-15, B-22, B-23), unchanged: the batch steps must give
its files byte for byte. The batch size (`shard_bases`) is forced small, so the fixture's stages hold more
than one batch (one shard each, DEC-21), and pass 1 is cut into more than one shard by source object.

Each action is the command build-stage runs, `python -m biobox.cli stage ingest pass1|plan|batch|merge`,
started by `runstep.run_step` as a real subprocess: its run record holds its CPU-minutes, and it publishes
only its own files (a pass-1 shard and the plan their `ingest_plan` records, a batch its shard, the merge
every table). The scratch runs stay in scratch; a run whose scratch was lost is made again
(`python -m biobox.l0seq.store rescan`) and must reproduce the sha256 its published record holds.

Then build-stage (card C-04): an ingest estimated over `SHARD_MARGIN` of the budget runs as those steps, a
killed build resumes at the next batch, lost scratch is rescanned, and each step's CPU-minutes are in its
run record. The build is stopped after ingest (an injected hook), so no MMseqs2 is needed.
"""

from __future__ import annotations

import dataclasses
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pytest
import yaml

from biobox import profile, runstep
from biobox import stage as bs
from biobox.l0seq import SourceObject, ingest, store
from biobox.l0seq.inputs import write_genome_input

sys.path.insert(0, str(Path(__file__).parent))
import test_builders_out as c08  # noqa: E402

pytestmark = [pytest.mark.r2, pytest.mark.covers("D-01", "D-03")]

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
PATHS = [FIX / o["file"] for o in CONTRACT["objects"]]  # lambda, T4, M. genitalium
SHARD_BASES = 200_000  # forced small: S0's two genomes are two batches
STAGES = {"S0": PATHS[:2], "S1": PATHS[1:]}  # S1 repeats T4 (stored already) and adds M. genitalium
NOON_MONDAY = datetime(2026, 11, 9, 12, 0)
PLAN_RECORDS = "l0seq/ingest_plan/"  # the pass-1 and plan records: a build's own, never in a manifest
BIG = 10**12  # pass-1 bytes that hold every object in one shard


def tree(root: Path) -> dict[str, str]:
    """sha256 of every file under the root except run-step's `runs/`."""
    root = Path(root)
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file() and p.relative_to(root).parts[0] != "runs"}


def data(t: dict[str, str]) -> dict[str, str]:
    return {k: v for k, v in t.items() if not k.startswith(PLAN_RECORDS)}


def objects_of(paths) -> list[SourceObject]:
    return [SourceObject(Path(p), "micro", "1", "INSDC (public)", True) for p in paths]


def objects_file(base: Path, paths, name: str = "objects.json") -> Path:
    base.mkdir(parents=True, exist_ok=True)
    f = base / name
    f.write_text(json.dumps([{"path": str(p), "source_name": "micro", "source_release": "1",
                              "license": "INSDC (public)", "redistributable": True} for p in paths]))
    return f


@pytest.fixture
def prof(tmp_path) -> profile.Profile:
    return c08.profile_at(tmp_path / "prof")


@pytest.fixture(scope="module")
def oracle(tmp_path_factory) -> Path:
    """Two stages ingested by the single-process `ingest`, one after the other, into one root."""
    root = tmp_path_factory.mktemp("b30-oracle") / "root"
    for stage, paths in STAGES.items():
        ingest(objects_of(paths), root, stage, shard_bases=SHARD_BASES)
    return root


def cmd(action: str, root: Path, stage: str, objs: Path, scratch: Path, *extra: str) -> list[str]:
    return [sys.executable, "-m", "biobox.cli", "stage", "ingest", action, "--root", str(root), "--stage", stage,
            "--objects", str(objs), "--tmp", str(scratch), *extra]


def rescan_cmd(root: Path, stage: str, objs: Path, scratch: Path, *what: str) -> list[str]:
    return [sys.executable, "-m", "biobox.l0seq.store", "rescan", "--root", str(root), "--stage", stage,
            "--objects", str(objs), "--tmp", str(scratch), *what]


def step(prof, root: Path, stage: str, name: str, command: list[str]) -> int:
    return runstep.run_step(prof, root, stage, name, command, clock=lambda: NOON_MONDAY, load=lambda: 0.0)


def pass1_and_plan(prof, root: Path, stage: str, objs: Path, scratch: Path, pass1_bytes: int) -> list[str]:
    """Every pass-1 shard of the stage, then the plan, each a run-step; the step names, in order."""
    n = len(store.pass1_shards(objects_of(o["path"] for o in json.loads(objs.read_text())), pass1_bytes))
    names = []
    for k in range(n):
        names.append(f"ingest.pass1-{k:05d}")
        assert step(prof, root, stage, names[-1], cmd("pass1", root, stage, objs, scratch, "--shard", str(k),
                                                     "--pass1-bytes", str(pass1_bytes),
                                                     "--shard-bases", str(SHARD_BASES))) == 0, names[-1]
    names.append("ingest.plan")
    assert step(prof, root, stage, "ingest.plan", cmd("plan", root, stage, objs, scratch,
                                                       "--shard-bases", str(SHARD_BASES))) == 0
    return names


def batches_and_merge(prof, root: Path, stage: str, objs: Path, scratch: Path) -> list[str]:
    names = []
    for j in range(len(store.read_ingest_plan(root, stage)["batches"])):
        names.append(f"ingest.batch-{j:04d}")
        assert step(prof, root, stage, names[-1], cmd("batch", root, stage, objs, scratch, "--batch", str(j))) == 0
    names.append("ingest")
    assert step(prof, root, stage, "ingest", cmd("merge", root, stage, objs, scratch)) == 0
    return names


def batch_ingest(prof, root: Path, stage: str, objs: Path, scratch: Path, pass1_bytes: int) -> list[str]:
    return pass1_and_plan(prof, root, stage, objs, scratch, pass1_bytes) + batches_and_merge(prof, root, stage, objs,
                                                                                              scratch)


def record(root: Path, stage: str, name: str) -> dict:
    return json.loads((Path(root) / f"runs/stage={stage}/{name}.json").read_text())


# --- the batch steps equal the single-process ingest -----------------------------------------------------


@pytest.mark.covers("R-09", "LY-10")
@pytest.mark.parametrize("pass1_bytes", [1, BIG], ids=["pass1-per-object", "pass1-one-shard"])
def test_batch_steps_equal_the_single_process_ingest_file_for_file(prof, tmp_path, oracle, pass1_bytes):
    root, scratch = tmp_path / "root", tmp_path / "scratch"
    for stage, paths in STAGES.items():
        objs = objects_file(tmp_path / stage, paths)
        names = batch_ingest(prof, root, stage, objs, scratch, pass1_bytes)
        n_pass1 = sum(n.startswith("ingest.pass1-") for n in names)
        n_batch = sum(n.startswith("ingest.batch-") for n in names)
        assert n_pass1 == (len(paths) if pass1_bytes == 1 else 1)
        assert n_batch == len(store.read_ingest_plan(root, stage)["batches"])
        if stage == "S0":
            assert n_batch == 2  # the batch size is forced small: one shard per batch (DEC-21)
        for name in names:  # each step its own process, its CPU-minutes in its run record, under the budget
            rec = record(root, stage, name)
            assert isinstance(rec["cpu_min"], float) and 0 <= rec["cpu_min"] <= prof.cpu_min_per_process
            assert rec["over_cpu_budget"] is False
            outs = [o["path"] for o in rec["outputs"]]
            if name.startswith("ingest.pass1-"):
                assert outs == [f"{PLAN_RECORDS}stage={stage}/pass1-{name[-5:]}.json"]
            elif name == "ingest.plan":
                assert outs == [f"{PLAN_RECORDS}stage={stage}/plan.json"]
            elif name.startswith("ingest.batch-"):
                shard = f"l0seq/seq_shards/stage={stage}/shard-{name[-4:]}.fa.gz"
                assert outs == [shard, f"{shard}.fai", f"{shard}.gzi"]
            else:
                assert not any(o.startswith((PLAN_RECORDS, "l0seq/seq_shards/")) for o in outs)
    got = tree(root)
    assert data(got) == tree(oracle)  # file for file, by sha256
    assert set(got) - set(tree(oracle)) and all(k.startswith(PLAN_RECORDS) for k in set(got) - set(tree(oracle)))


def test_input_order_and_pass1_cut_change_no_output(prof, tmp_path, oracle):
    """The batches are a function of the sorted inputs and the batch size alone: the objects listed in another
    order give the same files, the records included; another pass-1 cut gives the same batches and outputs."""
    stage, paths = "S0", STAGES["S0"]
    trees, plans = [], []
    for name, order, pass1_bytes in (("a", paths, 1), ("b", paths[::-1], 1), ("c", paths, BIG)):
        root = tmp_path / name / "root"
        batch_ingest(prof, root, stage, objects_file(tmp_path / name, order), tmp_path / name / "scratch", pass1_bytes)
        trees.append(tree(root))
        plans.append(store.read_ingest_plan(root, stage))
    assert trees[0] == trees[1]  # records and all
    assert data(trees[0]) == data(trees[2]) == {k: v for k, v in tree(oracle).items() if "stage=S1" not in k}
    assert plans[0]["batches"] == plans[2]["batches"]
    text = json.dumps(plans[0]) + "".join((tmp_path / "a/root" / k).read_text() for k in trees[0] if k.startswith(PLAN_RECORDS))
    assert str(tmp_path) not in text and str(FIX) not in text  # no absolute path in a published record


# --- lost scratch, changed sources, missing batches ------------------------------------------------------


def test_scratch_lost_after_the_plan_is_rescanned_to_the_same_bytes(prof, tmp_path, oracle, capfd):
    root, scratch = tmp_path / "root", tmp_path / "scratch"
    objs = objects_file(tmp_path / "S0", STAGES["S0"])
    pass1_and_plan(prof, root, "S0", objs, scratch, 1)
    runs = store.ingest_runs_dir(scratch, root, "S0")
    want = {p.name: p.read_bytes() for p in sorted(runs.iterdir())}
    assert sorted(want) == ["new_seq.parquet", "pass1-00000.parquet", "pass1-00001.parquet"]
    assert store.ingest_run_problems(root, "S0", scratch) == []
    shutil.rmtree(runs)  # a reboot cleared /tmp
    assert [(w, k) for w, k, _ in store.ingest_run_problems(root, "S0", scratch)] == [
        ("pass1", 0), ("pass1", 1), ("plan", None)]
    capfd.readouterr()
    assert step(prof, root, "S0", "ingest.batch-0000", cmd("batch", root, "S0", objs, scratch, "--batch", "0")) != 0
    assert "rescan" in capfd.readouterr().err
    before = tree(root)
    for what in (["--pass1", "0"], ["--pass1", "1"], ["--plan"]):
        res = subprocess.run(rescan_cmd(root, "S0", objs, scratch, *what), capture_output=True, text=True)
        assert res.returncode == 0, res.stderr
    assert tree(root) == before  # a rescan writes nothing into the root
    assert {p.name: p.read_bytes() for p in sorted(runs.iterdir())} == want
    assert store.ingest_run_problems(root, "S0", scratch) == []
    batches_and_merge(prof, root, "S0", objs, scratch)
    assert data(tree(root)) == {k: v for k, v in tree(oracle).items() if "stage=S1" not in k}


def test_a_rescan_that_differs_from_its_record_is_refused_naming_the_shard(prof, tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    paths = [shutil.copy(p, src / p.name) for p in STAGES["S0"]]
    root, scratch = tmp_path / "root", tmp_path / "scratch"
    objs = objects_file(tmp_path / "S0", paths)
    pass1_and_plan(prof, root, "S0", objs, scratch, 1)
    runs = store.ingest_runs_dir(scratch, root, "S0")
    shutil.rmtree(runs)
    k = next(i for i, sh in enumerate(store.pass1_shards(objects_of(paths), 1))
             if store.pass1_order(objects_of(paths))[sh[0]].path.name == PATHS[0].name)
    # the source is another one now: its pass-1 run cannot be the published one
    Path(paths[0]).write_bytes(gzip.compress(b">x changed\nACGTACGTACGT\n", mtime=0))
    res = subprocess.run(rescan_cmd(root, "S0", objs, scratch, "--pass1", str(k)), capture_output=True, text=True)
    assert res.returncode == 2
    assert f"pass-1 shard {k}" in res.stderr
    assert not (runs / f"pass1-{k:05d}.parquet").exists()


def test_a_source_changed_after_pass1_is_refused_by_its_batch(prof, tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    paths = [shutil.copy(p, src / p.name) for p in STAGES["S0"]]
    root, scratch = tmp_path / "root", tmp_path / "scratch"
    objs = objects_file(tmp_path / "S0", paths)
    pass1_and_plan(prof, root, "S0", objs, scratch, BIG)
    for p in paths:  # the same size, other bytes
        b = bytearray(Path(p).read_bytes())
        b[-9] ^= 0xFF
        Path(p).write_bytes(bytes(b))
    with pytest.raises(ValueError, match="changed since pass 1"):
        store.ingest_batch(objects_of(paths), root, "S0", 0, tmp=scratch, out=tmp_path / "out")
    assert not (tmp_path / "out" / "l0seq/seq_shards").exists()


def test_merge_refuses_until_every_batch_shard_is_published(prof, tmp_path):
    root, scratch = tmp_path / "root", tmp_path / "scratch"
    objs = objects_file(tmp_path / "S0", STAGES["S0"])
    pass1_and_plan(prof, root, "S0", objs, scratch, 1)
    assert step(prof, root, "S0", "ingest.batch-0000", cmd("batch", root, "S0", objs, scratch, "--batch", "0")) == 0
    before = tree(root)
    with pytest.raises(FileNotFoundError, match="shard-0001"):
        store.ingest_merge(objects_of(STAGES["S0"]), root, "S0", tmp=scratch, out=tmp_path / "out")
    assert tree(root) == before and not (tmp_path / "out").exists()


def test_every_action_refuses_a_stage_already_ingested(tmp_path, oracle):
    root = tmp_path / "root"
    shutil.copytree(oracle, root)
    objs, scratch, out = objects_of(STAGES["S0"]), tmp_path / "scratch", tmp_path / "out"
    with pytest.raises(FileExistsError, match="clean stage directory|already"):
        store.ingest_pass1(objs, root, "S0", 0, 1, shard_bases=SHARD_BASES, tmp=scratch, out=out)
    with pytest.raises(FileExistsError, match="clean stage directory|already"):
        store.ingest_plan(objs, root, "S0", shard_bases=SHARD_BASES, tmp=scratch, out=out)
    assert tree(root) == tree(oracle) and not out.exists()


def test_a_later_stage_with_other_store_parameters_is_refused_by_pass1(prof, tmp_path, oracle):
    root = tmp_path / "root"
    shutil.copytree(oracle, root)
    with pytest.raises(ValueError, match="store_params"):
        store.ingest_pass1(objects_of(PATHS[2:]), root, "S2", 0, 1, shard_bases=SHARD_BASES + 1,
                           tmp=tmp_path / "scratch", out=tmp_path / "out")
    assert not (tmp_path / "out").exists()


# --- build-stage runs the batch steps (card C-04) ----------------------------------------------------------

REL = "rS0.b30"
BUILD_SHARD_BASES = 200_000  # the S0 of test_builders_out (lambda, M. genitalium): two batches


class Killed(Exception):
    """Raised by the injected hook: the build process dies between two steps."""


@pytest.fixture(scope="module")
def work(tmp_path_factory) -> c08.Work:
    w = c08.Work(tmp_path_factory.mktemp("b30-work"))
    w.prepare()
    return w


@pytest.fixture(scope="module")
def build_oracle(work, tmp_path_factory) -> Path:
    """inputs and ingest of S0, as the existing code writes them in one process each."""
    root = tmp_path_factory.mktemp("b30-build-oracle") / "root"
    write_genome_input(root, "S0", c08.genome_inputs("S0"))
    ingest(c08.source_objects("S0"), root, "S0", shard_bases=BUILD_SHARD_BASES)
    return root


def build_of(work, root: Path) -> bs.Build:
    return bs.Build(root=root, stage="S0", release=REL, stages=("S0",), genomes=work.genomes("S0"),
                    objects=work.objects("S0"), threads=1, shard_bases=BUILD_SHARD_BASES, tmp=work.base / "scratch-b")


def runner(prof, root, stage, step_, command, **kw) -> int:
    return runstep.run_step(prof, root, stage, step_, command, clock=lambda: NOON_MONDAY, load=lambda: 0.0, **kw)


def stop_after_ingest(seen: list[str]):
    def hook(point: str) -> None:
        seen.append(point)
        if point == "after:ingest":
            raise Killed(point)
    return hook


def build(prof, b: bs.Build, hook) -> None:
    with pytest.raises(Killed):
        bs.build_stage(prof, b, run=runner, clock=lambda: NOON_MONDAY, load=lambda: 0.0, _hook=hook)


@pytest.fixture
def over_budget(monkeypatch, prof):
    """The measured ingest rate raised so that M. genitalium alone is estimated at 90% of the budget: S0's
    ingest is over `SHARD_MARGIN` of it in one process, and each object still fits in one."""
    mgen = next(p for p in PATHS if p.name.startswith("NC_000908"))
    rate = 0.9 * prof.cpu_min_per_process / len(gzip.decompress(mgen.read_bytes()))
    monkeypatch.setattr(bs, "INGEST_CPU_MIN", rate * bs.INGEST_BP)
    return prof


def ingest_records(root: Path) -> list[str]:
    return sorted(p.name.removesuffix(".json") for p in (Path(root) / "runs/stage=S0").glob("ingest*.json"))


def test_build_stage_runs_an_ingest_over_the_budget_as_pass1_plan_batch_and_merge_steps(
        work, build_oracle, over_budget, tmp_path):
    prof, root, seen = over_budget, tmp_path / "root", []
    build(prof, build_of(work, root), stop_after_ingest(seen))
    want = ["ingest", "ingest.batch-0000", "ingest.batch-0001", "ingest.pass1-00000", "ingest.pass1-00001",
            "ingest.plan"]
    assert ingest_records(root) == want
    # the pass-1 shards, then the plan, then the batches (several may run at once), then the merge
    assert seen[0] == "after:inputs" and seen[3] == "after:ingest.plan" and seen[-1] == "after:ingest"
    assert set(seen[1:3]) == {"after:ingest.pass1-00000", "after:ingest.pass1-00001"}
    assert set(seen[4:6]) == {"after:ingest.batch-0000", "after:ingest.batch-0001"} and len(seen) == 7
    for name in want:  # each step is its own process; its CPU-minutes are in its run record, under the budget
        rec = record(root, "S0", name)
        assert 0 <= rec["cpu_min"] <= prof.cpu_min_per_process and rec["over_cpu_budget"] is False
        assert rec["threads"] == 1
    assert data(tree(root)) == tree(build_oracle)
    assert not store.ingest_runs_dir(build_of(work, root).tmp, root, "S0").exists()  # scratch removed after the merge
    # each step had its own scratch directory (card C-13), removed once it exited 0
    for name in want:
        argv = record(root, "S0", name)["argv"]
        assert argv[argv.index("--tmp") + 1] == str(Path(build_of(work, root).tmp).resolve() / f"S0.{name}")
    assert not list(Path(build_of(work, root).tmp).glob("S0.ingest*"))


def test_build_stage_within_the_margin_runs_ingest_in_one_process(work, build_oracle, prof, tmp_path):
    root = tmp_path / "root"
    build(prof, build_of(work, root), stop_after_ingest([]))
    assert ingest_records(root) == ["ingest"]
    assert tree(root) == tree(build_oracle)


def test_build_stage_killed_between_batches_resumes_at_the_next_batch(work, build_oracle, over_budget, tmp_path):
    prof, root = over_budget, tmp_path / "root"
    one_at_a_time = dataclasses.replace(build_of(work, root), scan_jobs=1, scan_cpus=1)  # 1 ingest process at once

    def kill(point: str) -> None:
        if point == "after:ingest.batch-0000":
            raise Killed(point)

    build(prof, one_at_a_time, kill)
    done = {n: record(root, "S0", n) for n in ingest_records(root)}
    assert sorted(done) == ["ingest.batch-0000", "ingest.pass1-00000", "ingest.pass1-00001", "ingest.plan"]
    started: list[str] = []

    def recording(prof_, root_, stage, step_, command, **kw):
        if not (Path(root_) / "runs" / f"stage={stage}" / f"{step_}.json").exists():
            started.append(step_)
        return runner(prof_, root_, stage, step_, command, **kw)

    with pytest.raises(Killed):
        bs.build_stage(prof, build_of(work, root), run=recording, clock=lambda: NOON_MONDAY, load=lambda: 0.0,
                       _hook=stop_after_ingest([]))
    assert started == ["ingest.batch-0001", "ingest"]
    assert {n: record(root, "S0", n) for n in done} == done
    assert data(tree(root)) == tree(build_oracle)


def test_build_stage_rescans_lost_ingest_scratch_and_ends_byte_identical(work, build_oracle, over_budget, tmp_path):
    prof, root = over_budget, tmp_path / "root"
    b = build_of(work, root)

    def kill(point: str) -> None:
        if point == "after:ingest.plan":
            raise Killed(point)

    build(prof, b, kill)
    shutil.rmtree(store.ingest_runs_dir(b.tmp, root, "S0"))  # a reboot cleared /tmp
    seen: list[str] = []
    build(prof, b, stop_after_ingest(seen))
    # the published steps are passed over; the lost runs are made again (pass 1 first, at once), then the batches
    assert seen[:4] == ["after:inputs", "after:ingest.pass1-00000", "after:ingest.pass1-00001", "after:ingest.plan"]
    assert set(seen[4:6]) == {"after:ingest.pass1-00000.rescan", "after:ingest.pass1-00001.rescan"}
    assert seen[6] == "after:ingest.plan.rescan"
    assert set(seen[7:9]) == {"after:ingest.batch-0000", "after:ingest.batch-0001"} and seen[9:] == ["after:ingest"]
    assert data(tree(root)) == tree(build_oracle)


def test_an_object_over_the_budget_is_refused_before_any_step(work, monkeypatch, prof, tmp_path, capsys):
    mgen = next(p for p in PATHS if p.name.startswith("NC_000908"))
    rate = 1.5 * prof.cpu_min_per_process / len(gzip.decompress(mgen.read_bytes()))
    monkeypatch.setattr(bs, "INGEST_CPU_MIN", rate * bs.INGEST_BP)
    root = tmp_path / "root"
    assert bs.build_stage(prof, build_of(work, root), run=runner, clock=lambda: NOON_MONDAY,
                          load=lambda: 0.0) == runstep.EXIT_REFUSED
    err = capsys.readouterr().err
    assert "B-30" in err and mgen.name in err
    assert not root.exists() or tree(root) == {}


def test_the_step_list_passes_the_batch_size_to_ingest(work, prof, tmp_path):
    b = build_of(work, tmp_path / "root")
    ing = next(s for s in bs.steps(prof, b) if s.name == "ingest")
    assert ing.argv[ing.argv.index("--shard-bases") + 1] == str(BUILD_SHARD_BASES)
    plain = dataclasses.replace(b, shard_bases=None)
    assert "--shard-bases" not in next(s for s in bs.steps(prof, plain) if s.name == "ingest").argv
