"""Ring R2: a run record for each gene-calling batch (card B-50).

The world is card B-48's (`test_b48_genecall_batches.py`): S0 stored first, S1 the stage gene-called in batches. Here
build-stage calls S1 in 3 batches (`b50_batched_tree.build_in_batches`: the rates raised and the plan cut into 3 in this
process, where build-stage makes the plan). Oracles, each independent of the code under test:

- the files: the one-process call (`genecall.run`) for the gene-call files, and the whole root as card B-48's code built
  it in 3 batches before this card: `PARENT_DIGEST`, the sha256 over every file's path and sha256 but `runs/` and
  `l1/genecall_plan/`, printed by `tests/r2/b50_batched_tree.py` run with the `src` of commit c775930 (the tip of
  `plan/task-cards` before this card, whose code is card B-48's merge 466874a) on the path, in the test image on
  wilkcomp02, 2026-11-15 15:07 CST (*measured*: `PYTHONPATH=/tmp/if2453/b50/parent-c775930/src python
  tests/r2/b50_batched_tree.py /tmp/if2453/b50/parent-run`, through `singularity exec` of `bio-box-test_faec1f4.sif`);
- the CPU: this test process's own accounting of its children (`getrusage(RUSAGE_CHILDREN)` around the call that runs
  only gene calling), and each batch's own measurement of itself (`done.json`, `getrusage` inside the batch);
- the record's fields: the run record of S1's ingest from the same build (another step, the same `run-step`).
"""

from __future__ import annotations

import hashlib
import json
import resource
import shutil
import signal
import subprocess
import sys
from pathlib import Path

import pytest

from biobox import catalog
from biobox import stage as bs
from biobox.l0seq import ingest
from biobox.l0seq.inputs import write_genome_input
from biobox.l1 import genecall
from biobox.l1.genecall import GENE_CALL_VERSION as V

sys.path.insert(0, str(Path(__file__).parent))
import b50_batched_tree as b50  # noqa: E402
import test_b48_genecall_batches as b48  # noqa: E402

pytestmark = [pytest.mark.r2, pytest.mark.covers("B-50", "B-48", "D-01")]

# Moved once in card C-38d from 530c05e069f2aa51dd3d188c85fdb1164eadfdb50919b01f253dfdeaef6b523b when the two
# appended nullable columns (`n_bases_in`, `n_bases_rejected`, LY-09) entered the 3 `ledger/` Parquet files; all 23
# non-ledger files in `b50_batched_tree.tree(root)` kept their exact pre-C-38d sha256 (*measured*, wilkcomp02,
# 2026-11-18).
PARENT_DIGEST = "c139304e7e18f8d2c7416103ea0f8cc1ba93aedb7ab12338d4ebc2761bbf451c"
RUNS = "runs/stage=S1"
RECORDS = "l1/genecall_plan/stage=S1"
BATCHES = [f"genecall.batch-{j:05d}" for j in range(3)]
LISTING = sorted([*(f"{s}.json" for s in BATCHES), "genecall.json", "genecall.plan.json", "ingest.json", "inputs.json"])
# the fields of a published batch record that the batch's done.json holds too (the rest: what names its plan)
DONE_FIELDS = ("batch", "first", "genomes", "digest", "genes", "wrapping_genes", "modes", "skipped", "estimate_cpu_min")


class Stop(b48.Killed):
    pass


def stop_at(point: str, seen: list[str] | None = None):
    def hook(p: str) -> None:
        if seen is not None:
            seen.append(p)
        if p == point:
            raise Stop(p)
    return hook


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """The inputs of S0 and S1, and the one-process call of S1 (card B-48's oracle) into `one`."""
    d = tmp_path_factory.mktemp("b50w")
    w = b50.world_at(d)
    root = d / "root"
    for s in ("S0", "S1"):
        write_genome_input(root, s, b48.genome_inputs(s))
        ingest(b48.source_objects(d / "in", s), root, s)
    (d / "scratch").mkdir()
    w["root"] = root
    w["result"] = genecall.run(root, "S1", seq_stages=b48.SEQ_STAGES, out=d / "one", tmp=d / "scratch")
    w["one"] = b48.tree(d / "one")
    return w


def work_files(scratch: Path, name: str) -> list[Path]:
    return sorted((scratch / "genecall").glob(f"*/*/{name}"))


def build_measured(world, base: Path, jobs: int) -> dict:
    """S1 ingested by build-stage (stopped after ingest), then built again to the end of gene calling in 3 batches,
    this process's children's CPU taken around that second call; the work directory's plan and each batch's done.json
    read as each step ends (the merge removes them)."""
    root, code = b50.build_in_batches(world, base, jobs=jobs, hook=stop_at("after:ingest"))
    assert code is None
    plans: list[bytes] = []
    done: dict[int, dict] = {}
    names: list[str] = []

    def hook(p: str) -> None:
        if p == "after:genecall.plan":
            plans.extend(f.read_bytes() for f in work_files(base / "scratch", "plan.json"))
        if p.startswith("after:genecall.batch-"):
            j = int(p.rsplit("-", 1)[1])
            (f,) = work_files(base / "scratch", f"{genecall.batch_name(j)}/done.json")
            done[j] = json.loads(f.read_text())
        if p == "after:genecall":
            raise Stop(p)

    def batch_run(prof_, root_, stage, step, cmd, *, threads):
        names.append(step)
        return bs._batch_process(prof_, root_, stage, step, cmd, threads=threads)

    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    _, code = b50.build_in_batches(world, base, jobs=jobs, hook=hook, batch_run=batch_run, root=root)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    assert code is None
    cpu = (after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime) / 60
    return {"root": root, "children_cpu_min": cpu, "plans": plans, "done": done, "names": names, "base": base}


@pytest.fixture(scope="module")
def batched(world, tmp_path_factory):
    return build_measured(world, tmp_path_factory.mktemp("b50b"), jobs=1)


def rec(root: Path, step: str) -> dict:
    return json.loads((root / RUNS / f"{step}.json").read_text())


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


# --- the records -------------------------------------------------------------------------------------------------------


def test_b50_a_stage_in_3_batches_records_the_plan_each_batch_and_the_merge_and_their_cpu_sums_to_the_step(
        world, batched):
    """runs/stage=S1 holds genecall.plan, the three genecall.batch-<j> and genecall; each batch's record is a run-step
    record like ingest's (same fields, the same code and image) and lists its published record; the CPU-minutes of the
    four gene-call records and the plan's are the gene calling this process's children did, each batch's at least what
    the batch measured of itself and at most a second more (the measuring helper)."""
    root = batched["root"]
    assert sorted(p.name for p in (root / RUNS).iterdir()) == LISTING
    assert sorted(p.name for p in (root / RECORDS).iterdir()) == ["batch-00000.json", "batch-00001.json",
                                                                  "batch-00002.json", "plan.json"]
    assert batched["names"] == BATCHES  # each launched once, as its own run-step process
    ing = rec(root, "ingest")
    plan_bytes = (root / RECORDS / "plan.json").read_bytes()
    assert batched["plans"] == [plan_bytes]  # the plan published is the plan the batches read
    plan = json.loads(plan_bytes)
    assert len(plan["batches"]) == 3
    assert rec(root, "genecall.plan")["outputs"] == [{"path": f"{RECORDS}/plan.json", "bytes": len(plan_bytes),
                                                     "sha256": hashlib.sha256(plan_bytes).hexdigest()}]
    for j, step in enumerate(BATCHES):
        r = rec(root, step)
        assert set(r) == set(ing), set(r) ^ set(ing)
        for k in ("code_id", "code_version", "code_commit", "image_digest", "host", "nice"):
            assert r[k] == ing[k], k
        assert r["stage"] == "S1" and r["step"] == step and r["exit_code"] == 0 and not r["over_cpu_budget"]
        assert r["peak_rss_bytes"] > 0 and r["wall_s"] > 0
        argv = r["argv"][r["argv"].index("biobox.cli") + 1:]
        assert argv[:3] == ["stage", "genecall", "batch"] and argv[argv.index("--batch") + 1] == str(j)
        path = root / RECORDS / f"batch-{j:05d}.json"
        assert r["outputs"] == [{"path": f"{RECORDS}/batch-{j:05d}.json", "bytes": path.stat().st_size,
                                 "sha256": sha(path)}]
        published, done = json.loads(path.read_text()), batched["done"][j]
        assert {k: published[k] for k in DONE_FIELDS} == {k: done[k] for k in DONE_FIELDS}
        assert published["plan_sha256"] == hashlib.sha256(plan_bytes).hexdigest()
        assert published["estimate_cpu_min"] == plan["batches"][j]["estimate_cpu_min"]
        assert done["cpu_min"] - 0.0006 <= r["cpu_min"] <= done["cpu_min"] + 1 / 60, (j, r["cpu_min"], done["cpu_min"])
    assert sum(batched["done"][j]["genes"] for j in range(3)) == world["result"]["genes"]
    steps = ["genecall.plan", *BATCHES, "genecall"]
    total = sum(rec(root, s)["cpu_min"] for s in steps)
    children = batched["children_cpu_min"]
    print(f"B50 cpu_min: records {[rec(root, s)['cpu_min'] for s in steps]} sum {total:.6f}; children {children:.6f}; "
          f"the run-step processes' own {children - total:.6f}")
    # the children are the five steps' commands and the three batches' run-step processes around them
    assert total <= children + 1e-4
    assert children - total <= 3 * 0.5


def test_b50_a_stage_called_in_batches_has_the_files_it_had_before_this_card(world, batched):
    """Every file of the root but runs/ and the build records is the one card B-48's code wrote in 3 batches
    (`PARENT_DIGEST`), and the gene-call files are the one-process call's."""
    t = b50.tree(batched["root"])
    assert b48.genecall_files(t) == world["one"]
    assert b50.tree_digest(t) == PARENT_DIGEST


def test_b50_the_records_and_files_do_not_depend_on_how_many_batches_run_at_once(world, batched, tmp_path):
    other = build_measured(world, tmp_path, jobs=3)
    a, b = (b48.tree(x["root"]) for x in (batched, other))
    assert a == b  # every file but runs/, the build records l1/genecall_plan/ among them
    assert sorted(p.name for p in (other["root"] / RUNS).iterdir()) == LISTING


# --- kill and resume -----------------------------------------------------------------------------------------------------


def test_b50_a_batch_killed_mid_way_leaves_no_record_and_the_next_call_runs_only_the_batches_without_one(
        world, tmp_path):
    """Batch 1's run-step is sent SIGTERM while the batch calls genes: run-step forwards it and publishes nothing, so
    batch 1 has neither a run record nor a published record, and batch 2 is not launched. The next call runs batches 1
    and 2 only, then the merge, and the root is the one before this card."""
    flag = tmp_path / "held"

    def batch_run(prof_, root_, stage, step, cmd, *, threads):
        if step != "genecall.batch-00001":
            return bs._batch_process(prof_, root_, stage, step, cmd, threads=threads)
        held = [sys.executable, "-c", b48._HELD, *cmd[3:], str(flag), "calling"]
        argv = [sys.executable, "-m", "biobox.cli", "run-step", "--profile", str(bs._profile_file(prof_, tmp_path)),
                "--root", str(root_), "--stage", stage, "--step", step, "--threads", str(threads), "--", *held]
        p = subprocess.Popen(argv)
        try:
            assert b48._poll(flag.exists), "batch 1 never reached its calls"
            p.send_signal(signal.SIGTERM)
            return p.wait(timeout=120)
        finally:
            if p.poll() is None:
                p.kill()

    root, code = b50.build_in_batches(world, tmp_path, batch_run=batch_run)
    assert code == 143, code
    assert (root / RUNS / "genecall.batch-00000.json").is_file() and (root / RECORDS / "batch-00000.json").is_file()
    for j in (1, 2):
        assert not (root / RUNS / f"genecall.batch-{j:05d}.json").exists()
        assert not (root / RECORDS / f"batch-{j:05d}.json").exists()
    assert not (root / RUNS / "genecall.json").exists() and not b48.genecall_files(b48.tree(root))
    again: list[str] = []

    def count(prof_, root_, stage, step, cmd, *, threads):
        again.append(step)
        return bs._batch_process(prof_, root_, stage, step, cmd, threads=threads)

    _, code = b50.build_in_batches(world, tmp_path, batch_run=count, root=root)
    assert code is None
    assert again == ["genecall.batch-00001", "genecall.batch-00002"]
    assert sorted(p.name for p in (root / RUNS).iterdir()) == LISTING
    assert b50.tree_digest(b50.tree(root)) == PARENT_DIGEST


def test_b50_with_every_batch_recorded_and_no_merge_record_the_next_call_runs_only_the_merge(world, tmp_path):
    root, code = b50.build_in_batches(world, tmp_path, hook=stop_at("after:genecall.batch-00002"))
    assert code is None
    assert all((root / RUNS / f"{s}.json").is_file() for s in BATCHES)
    assert not (root / RUNS / "genecall.json").exists()
    runs = {p.name: p.read_bytes() for p in (root / RUNS).iterdir()}
    again: list[str] = []
    seen: list[str] = []

    def count(prof_, root_, stage, step, cmd, *, threads):
        again.append(step)
        return bs._batch_process(prof_, root_, stage, step, cmd, threads=threads)

    _, code = b50.build_in_batches(world, tmp_path, batch_run=count, root=root, hook=stop_at("after:genecall", seen))
    assert code is None
    assert again == [] and not [p for p in seen if p.endswith(".redo")]
    assert {k: (root / RUNS / k).read_bytes() for k in runs} == runs  # the batches' records stay as they were
    assert sorted(p.name for p in (root / RUNS).iterdir()) == LISTING
    assert b50.tree_digest(b50.tree(root)) == PARENT_DIGEST


@pytest.mark.parametrize("tampered", [False, True])
def test_b50_a_batch_whose_scratch_rows_are_lost_is_redone_and_checked_against_its_record(world, tmp_path, capfd,
                                                                                          tampered):
    """Every batch recorded, then the scratch directory lost (a reboot clears /tmp): the next call does every batch
    again into scratch, outside run-step (its records stay as published), each checked against its published record's
    counts and genome digest, then merges. A record that a batch does not reproduce (here: its gene count changed in
    the root) stops the build before the merge, naming the batch."""
    root, code = b50.build_in_batches(world, tmp_path, hook=stop_at("after:genecall.batch-00002"))
    assert code is None
    shutil.rmtree(tmp_path / "scratch" / "genecall")
    if tampered:
        p = root / RECORDS / "batch-00001.json"
        doc = json.loads(p.read_text())
        n = doc["genes"]
        doc["genes"] = n + 1 if len(str(n + 1)) == len(str(n)) else n - 1
        p.write_text(json.dumps(doc, sort_keys=True, indent=1) + "\n")
    published = {p.name: p.read_bytes() for p in (root / RECORDS).iterdir()}
    runs = {p.name: p.read_bytes() for p in (root / RUNS).iterdir()}
    capfd.readouterr()
    seen: list[str] = []
    again: list[str] = []

    def count(prof_, root_, stage, step, cmd, *, threads):
        again.append(step)
        return bs._batch_process(prof_, root_, stage, step, cmd, threads=threads)

    _, code = b50.build_in_batches(world, tmp_path, batch_run=count, root=root, hook=stop_at("after:genecall", seen))
    err = capfd.readouterr().err
    assert again == []  # no batch is run again as a run-step
    assert {p.name: p.read_bytes() for p in (root / RECORDS).iterdir()} == published
    if tampered:
        assert code == 2, code
        assert "batch 1" in err and "record" in err, err[-3000:]
        assert not (root / RUNS / "genecall.json").exists() and not b48.genecall_files(b48.tree(root))
        assert {p.name: p.read_bytes() for p in (root / RUNS).iterdir()} == runs
        return
    assert code is None
    assert [p for p in seen if p.endswith(".redo")] == [f"after:{s}.redo" for s in BATCHES]
    assert {k: (root / RUNS / k).read_bytes() for k in runs} == runs
    assert sorted(p.name for p in (root / RUNS).iterdir()) == LISTING
    assert b50.tree_digest(b50.tree(root)) == PARENT_DIGEST


# --- the commands --------------------------------------------------------------------------------------------------------


def test_b50_plan_publishes_the_plan_it_is_given_and_a_batch_its_record_only_under_out(world, tmp_path):
    """`stage genecall plan --plan P` writes P's bytes as the plan record under --out after checking its batches'
    digests against the stage's genome_input (a plan of other inputs is refused, nothing written); `stage genecall
    batch` writes its record under --out (again, from its output, when that is in scratch), and with neither --out nor
    run-step's staging, nothing; the root is never written."""
    before = b48.tree(world["root"])
    plan = b48.make_plan(world, tmp_path, 3)
    out = tmp_path / "out"
    b48.cli("stage", "genecall", "plan", "--root", str(world["root"]), "--stage", "S1", "--plan", str(plan),
            "--out", str(out), "--tmp", str(tmp_path))
    assert (out / RECORDS / "plan.json").read_bytes() == plan.read_bytes()
    bad = json.loads(plan.read_text())
    bad["batches"][1]["digest"] = "0" * 32
    other = tmp_path / "other" / "plan.json"
    other.parent.mkdir()
    other.write_text(json.dumps(bad, sort_keys=True, indent=1))
    p = b48.cli("stage", "genecall", "plan", "--root", str(world["root"]), "--stage", "S1", "--plan", str(other),
                "--out", str(tmp_path / "out2"), "--tmp", str(tmp_path), check=False)
    assert p.returncode == 2 and "digest" in p.stderr and "batch 1" in p.stderr, p.stderr
    assert not (tmp_path / "out2").exists()
    b48.cli(*b48.batch_args(world["root"], plan, 0), "--out", str(out))
    first = (out / RECORDS / "batch-00000.json").read_bytes()
    assert json.loads(first)["batch"] == 0
    b48.cli(*b48.batch_args(world["root"], plan, 0), "--out", str(tmp_path / "out3"))  # its output is in scratch
    assert (tmp_path / "out3" / RECORDS / "batch-00000.json").read_bytes() == first
    b48.cli(*b48.batch_args(world["root"], plan, 1))
    assert b48.tree(world["root"]) == before


# --- the release ---------------------------------------------------------------------------------------------------------


@pytest.mark.skipif(b48.c08.MMSEQS is None, reason="build-stage clusters with MMseqs2 (in the image)")
def test_b50_release_verify_and_the_manifest_treat_the_batch_records_as_ingests(world, tmp_path):
    """S0 and its release, then S1 and its release with S1 gene-called in 3 batches: the release verifies, its manifest
    lists no build record (`l1/genecall_plan/`, like `l0seq/ingest_plan/`, is never a release file), and `built_by`
    counts the plan's and the batches' run records with every other step's."""
    g, o = world["inputs"]["S0"]
    prof = b48.c08.profile_at(tmp_path / "prof")
    root = tmp_path / "root"
    b = bs.Build(root=root, stage="S0", release="rB48.0", stages=("S0",), genomes=g, objects=o, threads=1,
                 tmp=tmp_path / "scratch", scan_jobs=1, scan_cpus=1)
    assert bs.build_stage(prof, b, run=b48.runner, clock=lambda: b48.NOON_MONDAY, load=lambda: 0.0,
                          cpu=lambda: b48.IDLE) == 0
    seen: list[str] = []
    _, code = b50.build_in_batches(world, tmp_path, root=root, hook=seen.append)
    assert code == 0
    assert sorted(p.name for p in (root / RUNS).iterdir() if p.name.startswith("genecall")) == sorted(
        f"{s}.json" for s in [*BATCHES, "genecall", "genecall.plan"])
    m = catalog.read_manifest(catalog.manifest_path(root, b48.REL))
    assert catalog.verify_manifest(m, root) == []
    assert not [f["path"] for f in m["files"] + m["artifacts"] if "genecall_plan" in f["path"]]
    n = sum(1 for s in ("S0", "S1", b48.REL) for p in (root / "runs" / f"stage={s}").glob("*.json")
            if not p.name.endswith(".publishing.json"))
    got = json.loads(catalog.built_by(root, b48.REL, ["S0", "S1"]))
    assert sum(c["steps"] for c in got["code"]) == n and got["unknown_steps"] == 0
    assert b48.genecall_files(b48.tree(root)) == world["one"]
    assert bs.release_proteins(root, ["S0", "S1"], V)[0] is not None
