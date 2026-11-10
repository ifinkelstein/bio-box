"""C-04: `bio-box build-stage` runs a stage and its release as a fixed list of `run-step` steps.

The oracle is the direct build of card C-08's test (`test_builders_out.direct`): every builder
called as the existing code calls it, one after the other, in one root. One `build-stage` call on
an empty root must give the same files (sha256, run records excluded) and the same manifest hash.
Its steps are the same commands that test runs one by one: `python -m biobox.cli stage ...` and
`python -m biobox.cli release ...`, each started by `runstep.run_step` as a real subprocess.

Then: a second call starts no step command; a build stopped after a middle step (an injected
hook, or a `STOP` file) resumes on the next call and ends byte-identical; the CPU-minute check
refuses ingest before any step (card B-30); a scan estimated over the budget runs as a plan step,
one step per shard and a merge step, whose files equal the unsharded build's byte for byte, and a
build killed between shards resumes at the next shard (card B-24); a second stage and a release
over both stages, built by a second call, change no file of the first.

The test profile's small-file floor is 1,000 bytes, as in `test_builders_out`: at the pod's real
64 MiB floor a second release's manifest is refused until card C-09 merges.
"""

from __future__ import annotations

import dataclasses
import gzip
import json
import shutil
import sys
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from biobox import catalog, cli, runstep
from biobox import stage as bs
from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION as V

sys.path.insert(0, str(Path(__file__).parent))
import test_builders_out as c08  # noqa: E402

direct = c08.direct  # the module-scoped direct build of S0 and release c08.REL (a fixture)

pytestmark = [
    pytest.mark.r2,
    pytest.mark.covers("D-01", "D-03"),
    pytest.mark.skipif(c08.MMSEQS is None, reason="build-stage clusters with MMseqs2 (in the image)"),
]

REL, REL2 = c08.REL, c08.REL2  # same release names as the oracle, so the manifests can be equal
NAMES = ["inputs", "ingest", "genecall", "features", "igr", "tapes", "domains-test-1",
         "cluster", "tokens", "postings", "dim", "manifest", "catalog"]


class Killed(Exception):
    """Raised by the injected hook: the build process dies between two steps."""


def profile_path(direct) -> Path:
    return Path(direct["prof"].tmp).parent / "profile.yaml"  # where c08.profile_at wrote it


def build_of(direct, root: Path, stage: str = "S0", release: str = REL, stages=("S0",)) -> bs.Build:
    w = direct["work"]
    return bs.Build(root=Path(root), stage=stage, release=release, stages=None if stages is None else tuple(stages),
                    genomes=w.genomes(stage), objects=w.objects(stage),
                    libraries=(bs.Library("test", "1", w.lib, "gathering"),), threads=1)


def argv_of(direct, root: Path, stage: str = "S0", release: str = REL, stages: str = "S0") -> list[str]:
    w = direct["work"]
    return ["build-stage", "--profile", str(profile_path(direct)), "--root", str(root), "--stage", stage,
            "--release", release, "--stages", stages, "--genomes", str(w.genomes(stage)),
            "--objects", str(w.objects(stage)), "--library", "test", "1", str(w.lib), "gathering",
            "--threads", "1"]


def runner(prof, root, stage, step, cmd, **kw) -> int:
    """run-step at a fixed hour and load, as `test_builders_out.run_cli` runs it."""
    return runstep.run_step(prof, root, stage, step, cmd, clock=lambda: c08.NOON_MONDAY, load=lambda: 0.0, **kw)


def run(direct, prof, root, hook=None, **kw) -> int:
    return bs.build_stage(prof, build_of(direct, root, **kw), run=runner, _hook=hook)


def records(root: Path) -> dict[str, bytes]:
    runs = Path(root) / "runs"
    return {p.relative_to(runs).as_posix(): p.read_bytes() for p in sorted(runs.rglob("*")) if p.is_file()} \
        if runs.exists() else {}


def record_names(root: Path) -> set[str]:
    return {k.removesuffix(".json") for k in records(root)}


def check_equals_direct(direct, root: Path) -> None:
    got, want = c08.tree(root), c08.tree(direct["root"])
    assert sorted(got) == sorted(want)
    assert c08.data(got) == c08.data(want)
    man = c08.check_manifest(root, REL, direct["manifest"])
    assert catalog.verify_manifest(man) == []
    c08.check_catalog(root, REL, man, Path(direct["prof"].tmp) / "staging")


# --- one call builds the stage and its release -----------------------------------------------------


@pytest.fixture(scope="module")
def built(direct, tmp_path_factory):
    root = tmp_path_factory.mktemp("c04-built") / "root"
    assert cli.main(argv_of(direct, root)) == 0
    return root


def test_the_step_list_is_fixed_and_a_pure_function_of_the_command_line(direct, tmp_path):
    prof, w = direct["prof"], direct["work"]
    b = build_of(direct, tmp_path / "root")
    got = bs.steps(prof, b)
    assert [s.name for s in got] == NAMES
    assert got == bs.steps(prof, build_of(direct, tmp_path / "root"))  # nothing but its arguments decides it
    assert [s.run_stage for s in got] == ["S0"] * 7 + [REL] * 6
    assert next(s for s in got if s.name == "features").argv is None  # the empty slot: no command
    # one scan step per library, in the order given
    two = dataclasses.replace(b, libraries=(bs.Library("b", "2", w.lib, "gathering"), bs.Library("a", "1", w.lib, "none")))
    assert [s.name for s in bs.steps(prof, two) if s.name.startswith("domains")] == ["domains-b-2", "domains-a-1"]
    # the release steps take the stages from release_def: --stages is passed only when given
    none = dataclasses.replace(b, stages=None)
    assert all("--stages" not in s.argv for s in bs.steps(prof, none) if s.argv)
    assert [s.name for s in got if s.argv and "--stages" in s.argv] == ["cluster", "tokens", "postings", "dim"]


def test_one_call_builds_the_fixture_as_the_direct_calls_do(direct, built):
    check_equals_direct(direct, built)
    # every step published its own files and nothing else: the run records partition the root
    steps = [s for s in bs.steps(direct["prof"], build_of(direct, built)) if s.argv is not None]
    outs = [p for s in steps for p in c08.published(built, s)]
    assert sorted(outs) == sorted(c08.tree(built)) and len(outs) == len(set(outs))


def test_a_second_call_starts_no_step_command(direct, built, monkeypatch):
    before, files = records(built), c08.tree(built)

    def no_command(*a, **k):
        raise AssertionError("a step command was started on the second call")

    monkeypatch.setattr(runstep, "_execute", no_command)
    assert cli.main(argv_of(direct, built)) == 0
    assert records(built) == before  # byte-identical run records, and no new one
    assert c08.tree(built) == files


# --- stopped and resumed -----------------------------------------------------------------------------


@pytest.mark.covers("D-03")
def test_a_build_killed_after_a_middle_step_resumes_and_ends_byte_identical(direct, tmp_path):
    root = tmp_path / "root"

    def kill(point: str) -> None:
        if point == "after:tapes":
            raise Killed(point)

    with pytest.raises(Killed):
        run(direct, direct["prof"], root, hook=kill)
    done = records(root)
    assert {k.removesuffix(".json") for k in done} == {f"stage=S0/{n}" for n in ("inputs", "ingest", "genecall",
                                                                                  "igr", "tapes")}
    seen: list[str] = []
    assert run(direct, direct["prof"], root, hook=seen.append) == 0
    assert seen == [f"after:{n}" for n in NAMES]
    assert {k: v for k, v in records(root).items() if k in done} == done  # the earlier steps did not run again
    check_equals_direct(direct, root)


@pytest.mark.covers("D-03")
def test_a_stop_file_stops_before_the_next_step_and_a_later_call_continues(direct, tmp_path, capsys):
    root = tmp_path / "root"

    def stop(point: str) -> None:
        if point == "after:genecall":
            (root / "STOP").write_text("")

    assert run(direct, direct["prof"], root, hook=stop) == runstep.EXIT_STOP
    assert "STOP" in capsys.readouterr().err
    assert record_names(root) == {"stage=S0/inputs", "stage=S0/ingest", "stage=S0/genecall"}
    before = c08.tree(root)
    assert run(direct, direct["prof"], root) == runstep.EXIT_STOP  # still there: nothing runs
    assert c08.tree(root) == before
    (root / "STOP").unlink()
    assert run(direct, direct["prof"], root) == 0
    check_equals_direct(direct, root)


# --- the CPU-minute checks -----------------------------------------------------------------------------


def fixture_fasta_bytes(direct, stage: str) -> int:
    """Bytes of FASTA the stage's objects hold, decompressed here independently of `stage.py`."""
    objs = json.loads(direct["work"].objects(stage).read_text())
    return sum(len(gzip.decompress(Path(o["path"]).read_bytes())) for o in objs)


def fixture_residues(root: Path, stage: str) -> int:
    """Residues of the stage's distinct proteins, from the sequences themselves."""
    files = sorted((Path(root) / f"l0seq/protein_seq/gene_call_version={V}/stage={stage}").rglob("*.parquet"))
    return sum(len(s) for f in files for s in pq.read_table(f, columns=["seq"]).column("seq").to_pylist())


def test_the_ingest_estimate_refuses_before_any_step(direct, tmp_path, capsys):
    root = tmp_path / "root"
    n = fixture_fasta_bytes(direct, "S0")
    assert n > 0
    prof = dataclasses.replace(direct["prof"], cpu_min_per_process=0)
    assert run(direct, prof, root) == runstep.EXIT_REFUSED
    err = capsys.readouterr().err
    assert "B-30" in err and f"{n:,}" in err
    assert c08.tree(root) == {} and records(root) == {}


# --- a scan over the budget runs in shard steps (card B-24) ---------------------------------------------

SCAN = "domains-test-1"
SHARD_RECORDS = "l1/domain_scan/"  # the plan and shard records: a build's own, never in a manifest


def small_budget(direct):
    """A profile whose budget lies between the two estimates on the fixture (measured rates in stage.py)."""
    residues = fixture_residues(direct["root"], "S0")
    assert fixture_fasta_bytes(direct, "S0") * bs.INGEST_CPU_MIN / bs.INGEST_BP < 1 < residues * bs.SCAN_CPU_MIN_PER_M / 1e6
    return dataclasses.replace(direct["prof"], cpu_min_per_process=1)


def shard_names(root: Path) -> list[str]:
    return sorted(n.removeprefix("stage=S0/") for n in record_names(root) if n.startswith(f"stage=S0/{SCAN}.shard-"))


def check_sharded_equals_direct(direct, root: Path) -> None:
    """Every file but the shard records is the direct (unsharded) build's, byte for byte; so are the manifest
    and the catalog, which never list the shard records."""
    got = {k: v for k, v in c08.tree(root).items() if not k.startswith(SHARD_RECORDS)}
    want = c08.tree(direct["root"])
    assert sorted(got) == sorted(want)
    assert c08.data(got) == c08.data(want)
    man = c08.check_manifest(root, REL, direct["manifest"])
    assert catalog.verify_manifest(man) == []
    c08.check_catalog(root, REL, man, Path(direct["prof"].tmp) / "staging")


@pytest.fixture(scope="module")
def sharded(direct, tmp_path_factory):
    """One build-stage call with a budget below the fixture's scan estimate, and the hook points it passed."""
    root = tmp_path_factory.mktemp("c04-sharded") / "root"
    seen: list[str] = []
    assert run(direct, small_budget(direct), root, hook=seen.append) == 0
    return root, seen


@pytest.mark.covers("D-01")
def test_a_scan_over_the_budget_runs_in_shard_steps_whose_merge_equals_the_unsharded_build(direct, sharded):
    root, seen = sharded
    budget = small_budget(direct).cpu_min_per_process
    shards = shard_names(root)
    assert len(shards) >= 2 and shards == [f"{SCAN}.shard-{k:05d}" for k in range(len(shards))]
    stage_steps = ["inputs", "ingest", "genecall", "igr", "tapes"]
    assert record_names(root) == ({f"stage=S0/{n}" for n in stage_steps + [f"{SCAN}.plan", *shards, SCAN]}
                                  | {f"stage={REL}/{n}" for n in NAMES[7:]})
    # each shard is its own step, in order, between the plan and the merge (the scan step's own record)
    assert seen == ([f"after:{n}" for n in NAMES[:6]] + [f"after:{SCAN}.plan"] + [f"after:{s}" for s in shards]
                    + [f"after:{n}" for n in NAMES[6:]])
    # every shard fits the budget at the measured rate, and together they hold every protein once
    plan = domains.read_plan(root, "S0", "test", "1")
    assert len(plan["shards"]) == len(shards)
    assert all(s["residues"] * bs.SCAN_CPU_MIN_PER_M / 1e6 <= budget for s in plan["shards"])
    assert sum(s["residues"] for s in plan["shards"]) == fixture_residues(root, "S0")
    check_sharded_equals_direct(direct, root)


def test_a_second_call_after_a_sharded_scan_starts_nothing(direct, sharded, tmp_path, monkeypatch):
    root = tmp_path / "root"
    shutil.copytree(sharded[0], root, symlinks=True)
    before, files = records(root), c08.tree(root)

    def no_command(*a, **k):
        raise AssertionError("a step command was started on the second call")

    monkeypatch.setattr(runstep, "_execute", no_command)
    assert run(direct, small_budget(direct), root) == 0
    assert records(root) == before and c08.tree(root) == files


@pytest.mark.covers("D-03")
def test_a_build_killed_between_shards_resumes_at_the_next_shard(direct, sharded, tmp_path):
    root = tmp_path / "root"
    first = f"after:{SCAN}.shard-00000"

    def kill(point: str) -> None:
        if point == first:
            raise Killed(point)

    with pytest.raises(Killed):
        run(direct, small_budget(direct), root, hook=kill)
    done = records(root)
    assert {k.removesuffix(".json") for k in done if k.startswith(f"stage=S0/{SCAN}")} == {
        f"stage=S0/{SCAN}.plan", f"stage=S0/{SCAN}.shard-00000"}
    started: list[str] = []  # the steps that had no run record when build-stage reached them

    def recording(prof, root_, stage, step, cmd, **kw):
        if not (Path(root_) / "runs" / f"stage={stage}" / f"{step}.json").exists():
            started.append(step)
        return runner(prof, root_, stage, step, cmd, **kw)

    assert bs.build_stage(small_budget(direct), build_of(direct, root), run=recording) == 0
    assert started[0] == f"{SCAN}.shard-00001"  # the next shard, then the rest of the build
    assert started == [*shard_names(root)[1:], SCAN, *NAMES[7:]]
    assert {k: v for k, v in records(root).items() if k in done} == done  # the plan and shard 0 did not run again
    check_sharded_equals_direct(direct, root)


# --- refusals before any step ---------------------------------------------------------------------------


def test_a_release_without_the_stage_or_without_a_definition_is_refused_before_any_step(direct, tmp_path, capsys):
    root = tmp_path / "root"
    assert run(direct, direct["prof"], root, stages=("S1",)) == runstep.EXIT_REFUSED
    assert "S0" in capsys.readouterr().err
    assert run(direct, direct["prof"], root, stages=None) == runstep.EXIT_REFUSED  # no --stages, no release_def
    assert "release_def" in capsys.readouterr().err
    assert c08.tree(root) == {} and records(root) == {}


# --- a second stage and release ---------------------------------------------------------------------------


@pytest.mark.covers("R-09")
def test_a_second_stage_and_release_change_no_file_of_the_first(direct, built, tmp_path):
    w = direct["work"]
    root = tmp_path / "root"
    shutil.copytree(built, root, symlinks=True)
    first, first_runs = c08.tree(root), records(root)
    assert cli.main(argv_of(direct, root, stage="S1", release=REL2, stages="S0,S1")) == 0
    after = c08.tree(root)
    assert {k: after[k] for k in first} == first  # every first-stage and first-release file unchanged
    assert {k: v for k, v in records(root).items() if k in first_runs} == first_runs
    later = [s for s in bs.steps(direct["prof"], build_of(direct, root, "S1", REL2, ("S0", "S1"))) if s.argv]
    assert all("stage=S0" not in p for s in later for p in c08.published(root, s))
    # the same steps called directly on the direct store give the same files and the same manifest
    oracle = tmp_path / "direct"
    shutil.copytree(direct["root"], oracle, symlinks=True)
    for step in c08.stage_steps("S1", w) + c08.release_steps(REL2, ["S0", "S1"], w):
        step.direct(oracle)
    assert c08.data(after) == c08.data(c08.tree(oracle))
    assert sorted(after) == sorted(c08.tree(oracle))
    man = c08.check_manifest(root, REL2, catalog.read_manifest(catalog.manifest_path(oracle, REL2)))
    assert catalog.verify_manifest(man) == []
    c08.check_catalog(root, REL2, man, Path(direct["prof"].tmp) / "staging")
