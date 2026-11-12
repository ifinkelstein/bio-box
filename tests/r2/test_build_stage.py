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
over both stages, built by a second call, change no file of the first. Clustering in parts (card
B-41): forced into parts, the build equals the one-process build; it is chosen above 80% of the
budget; a scratch database lost after its record is made again and the build ends byte-identical.

The test profile's small-file floor is 1,000 bytes, as in `test_builders_out`: at the pod's real
64 MiB floor a second release's manifest is refused until card C-09 merges.
"""

from __future__ import annotations

import dataclasses
import gzip
import json
import shutil
import subprocess
import sys
import threading
from datetime import datetime
from pathlib import Path

import pyarrow.parquet as pq
import pytest
import yaml

from biobox import catalog, cli, runstep
from biobox import stage as bs
from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l2 import cluster as l2cluster

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


# --- shards several at a time, planned with a margin (card C-11) -----------------------------------------

NOON = c08.NOON_MONDAY  # the day cap of the test profile (4) holds: the allowance of every test below
SUNDAY_NIGHT = datetime(2026, 11, 8, 23, 0)  # the night cap (6) holds when the host is idle
TIMEOUT = 120  # seconds a shard waits for the others it must run beside before the test fails


def idle() -> float:
    return 0.0


def shard(k: int) -> str:
    return f"{SCAN}.shard-{k:05d}"


def patched_residues(direct):
    """Residues per shard that cut the fixture into six or seven shards, so four jobs have work at once (the
    plan's residues come from `shard_residues`; the measured rate gives the fixture at most two shards)."""
    n = fixture_residues(direct["root"], "S0") // 6 + 1
    return lambda prof: n


class Counting:
    """A runner that counts the shard steps running at once (`most`) and records each shard call with its
    thread count. The first `meet` shard calls wait at a barrier until `meet` of them are running, so a build
    that runs fewer at once breaks the barrier (after TIMEOUT) instead of passing by timing. A shard named in
    `fail` returns that code (after the barrier) without running."""

    def __init__(self, meet: int = 1, fail: dict[str, int] | None = None):
        self.lock = threading.Lock()
        self.alive = self.most = 0
        self.calls: list[tuple[str, int]] = []
        self.meet, self.fail = meet, dict(fail or {})
        self.barrier = threading.Barrier(meet, timeout=TIMEOUT) if meet > 1 else None

    def __call__(self, prof, root, stage, step, cmd, **kw) -> int:
        if ".shard-" not in step:
            return runner(prof, root, stage, step, cmd, **kw)
        with self.lock:
            self.alive += 1
            self.most = max(self.most, self.alive)
            self.calls.append((step, kw["threads"]))
            first = len(self.calls) <= self.meet
        try:
            if first and self.barrier is not None:
                self.barrier.wait()
            if step in self.fail:
                return self.fail[step]
            return runner(prof, root, stage, step, cmd, **kw)
        finally:
            with self.lock:
                self.alive -= 1


def jobs_of(direct, root: Path, jobs: int, cpus: int) -> bs.Build:
    return dataclasses.replace(build_of(direct, root), scan_jobs=jobs, scan_cpus=cpus)


def profile_file(direct, prof, path: Path) -> Path:
    path.write_text(yaml.safe_dump(dataclasses.asdict(prof)))
    return path


def cli_argv(direct, root: Path, prof_file: Path, *extra: str) -> list[str]:
    argv = argv_of(direct, root)
    argv[argv.index("--profile") + 1] = str(prof_file)
    return argv + list(extra)


@pytest.fixture(scope="module")
def by_jobs(direct, tmp_path_factory):
    """The fixture built with 1, 2 and 4 scan jobs over six or seven shards. 1 job of 4 threads and 4 jobs of
    1 thread run through a counting runner (the shards' run-steps in this process); 2 jobs of 2 threads run
    through the command line, where each shard is its own `bio-box run-step` process."""
    out = {}
    prof = small_budget(direct)
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(bs, "shard_residues", patched_residues(direct))
        for jobs, cpus in ((1, 4), (4, 1)):
            root = tmp_path_factory.mktemp(f"c11-j{jobs}") / "root"
            count = Counting(meet=jobs)
            assert bs.build_stage(prof, jobs_of(direct, root, jobs, cpus), run=count, clock=lambda: NOON,
                                  load=idle) == 0
            out[jobs] = {"root": root, "count": count, "cpus": cpus}
        base = tmp_path_factory.mktemp("c11-j2")
        pf = profile_file(direct, prof, base / "profile.yaml")
        assert cli.main(cli_argv(direct, base / "root", pf, "--scan-jobs", "2", "--scan-cpus", "2")) == 0
        out[2] = {"root": base / "root", "count": None, "cpus": 2}
    return out


@pytest.mark.covers("D-01")
def test_with_1_2_and_4_scan_jobs_the_root_equals_the_unsharded_build(direct, by_jobs):
    trees = {}
    for jobs, b in sorted(by_jobs.items()):
        root = b["root"]
        assert len(shard_names(root)) >= 5, "the fixture must give four jobs work at once"
        check_sharded_equals_direct(direct, root)
        assert not (root / ".incoming").exists()  # nothing of the publishing area is left behind
        trees[jobs] = c08.data(c08.tree(root))  # every file, the plan and the shard records included
    assert trees[1] == trees[2] == trees[4]


def test_the_plan_is_the_same_file_whatever_the_jobs_and_threads(direct, by_jobs):
    plan = domains.scan_dir("S0", "test", "1") + "/plan.json"
    assert (by_jobs[1]["root"] / plan).read_bytes() == (by_jobs[4]["root"] / plan).read_bytes()
    assert (by_jobs[2]["root"] / plan).read_bytes() == (by_jobs[4]["root"] / plan).read_bytes()


def test_shard_steps_alive_at_once_never_exceed_the_jobs_and_each_takes_its_threads(direct, by_jobs):
    for jobs in (1, 4):
        count, root = by_jobs[jobs]["count"], by_jobs[jobs]["root"]
        assert count.most == jobs  # the barrier made `jobs` run at once; the count shows never more
        assert sorted(s for s, _ in count.calls) == shard_names(root)  # each shard started once
        assert {t for _, t in count.calls} == {by_jobs[jobs]["cpus"]}
    # through the command line each shard is a run-step of its own, with --scan-cpus threads
    root = by_jobs[2]["root"]
    recs = [json.loads((root / "runs/stage=S0" / f"{n}.json").read_text()) for n in shard_names(root)]
    assert {r["threads"] for r in recs} == {2}
    assert all(r["argv"][r["argv"].index("--cpus") + 1] == "2" for r in recs)


def test_shards_are_planned_at_80_percent_of_the_budget(direct, sharded):
    assert bs.SHARD_MARGIN == 0.8
    # derived: 300 CPU-minutes × 0.8 / 7.96 CPU-minutes per million residues = 30,150,753.8 residues (the rate of
    # C-10's P1 scan; card C-13 migrated this from 7.30 and 32,876,712)
    assert bs.shard_residues(dataclasses.replace(direct["prof"], cpu_min_per_process=300)) == 30_150_753
    prof = small_budget(direct)
    plan = domains.read_plan(sharded[0], "S0", "test", "1")
    assert plan["shard_residues"] == bs.shard_residues(prof)
    assert all(s["residues"] * bs.SCAN_CPU_MIN_PER_M / 1e6 <= 0.8 * prof.cpu_min_per_process for s in plan["shards"])


def test_a_scan_estimated_above_80_percent_of_the_budget_is_sharded(direct, sharded, tmp_path, monkeypatch):
    """The same margin decides whether a scan is sharded: an estimate of 90% of the budget is."""
    root = tmp_path / "root"
    shutil.copytree(sharded[0], root, symlinks=True, ignore=shutil.ignore_patterns("domain_scan"))
    prof, residues = small_budget(direct), fixture_residues(root, "S0")
    b = build_of(direct, root)
    scan = next(s for s in bs.steps(prof, b) if s.name == SCAN)
    for share, want in ((0.9, True), (0.7, False)):
        monkeypatch.setattr(bs, "SCAN_CPU_MIN_PER_M", share * prof.cpu_min_per_process * 1e6 / residues)
        assert bs._sharded(prof, b, root, scan) is want


@pytest.mark.covers("D-03")
def test_a_failed_shard_stops_launches_lets_running_shards_publish_and_the_next_call_runs_the_rest(
        direct, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(bs, "shard_residues", patched_residues(direct))
    prof, root = small_budget(direct), tmp_path / "root"
    count = Counting(meet=4, fail={shard(3): 3})
    b = jobs_of(direct, root, 4, 1)
    assert bs.build_stage(prof, b, run=count, clock=lambda: NOON, load=idle) == 3  # the failed shard's code
    assert shard(3) in capsys.readouterr().err
    n = len(domains.read_plan(root, "S0", "test", "1")["shards"])
    assert n >= 5
    assert sorted(s for s, _ in count.calls) == [shard(k) for k in range(4)]  # nothing launched after it
    assert shard_names(root) == [shard(k) for k in range(3)]  # the shards running beside it published
    assert not any(k.startswith((f"stage=S0/{SCAN}.json", f"stage={REL}/")) for k in records(root))
    done = records(root)
    started: list[str] = []  # the steps that had no run record when build-stage reached them

    def recording(prof_, root_, stage, step, cmd, **kw):
        if not (Path(root_) / "runs" / f"stage={stage}" / f"{step}.json").exists():
            started.append(step)
        return runner(prof_, root_, stage, step, cmd, **kw)

    assert bs.build_stage(prof, b, run=recording, clock=lambda: NOON, load=idle) == 0
    assert sorted(started[:n - 3]) == [shard(k) for k in range(3, n)]  # only the shards without a record
    assert started[n - 3:] == [SCAN, *NAMES[7:]]  # then the merge and the release
    assert {k: v for k, v in records(root).items() if k in done} == done
    check_sharded_equals_direct(direct, root)


@pytest.mark.covers("D-03")
def test_stop_is_checked_before_each_shard_launch(direct, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(bs, "shard_residues", patched_residues(direct))
    prof, root = small_budget(direct), tmp_path / "root"

    def stop(point: str) -> None:
        if point.startswith(f"after:{SCAN}.shard-"):
            (root / "STOP").write_text("")

    count = Counting(meet=2)
    assert bs.build_stage(prof, jobs_of(direct, root, 2, 2), run=count, clock=lambda: NOON, load=idle,
                          _hook=stop) == runstep.EXIT_STOP
    assert "STOP" in capsys.readouterr().err
    assert sorted(s for s, _ in count.calls) == [shard(0), shard(1)]  # the two running; none launched after
    assert 1 <= len(shard_names(root)) <= 2 and set(shard_names(root)) <= {shard(0), shard(1)}
    (root / "STOP").unlink()
    assert bs.build_stage(prof, jobs_of(direct, root, 2, 2), run=runner, clock=lambda: NOON, load=idle) == 0
    check_sharded_equals_direct(direct, root)


@pytest.mark.covers("D-03")
def test_shards_publishing_at_once_do_not_remove_each_others_incoming_area(direct, tmp_path, monkeypatch):
    """Before card C-12, run-step removed `<root>/.incoming` after its renames when it listed it empty
    (`_remove_leftovers`); a shard that copied its outputs into `.incoming` between that listing and the `rmdir`
    made the `rmdir` fail, after the first shard's files were in the root and before its run record. Ordered
    here without timing: shard 1, before its copy (rsync, or copytree where rsync is absent, as in the image),
    waits until shard 0 has looked at `.incoming` (then shard 0 waits for that copy before it goes on) or, as
    run-step since card C-12 never lists it, until shard 0 has ended. Both must publish. (Migrated by card
    C-12: the wait for shard 0's listing became a wait for its listing or its end.)"""
    monkeypatch.setattr(bs, "shard_residues", patched_residues(direct))
    prof, root = small_budget(direct), tmp_path / "root"
    incoming = root.resolve() / ".incoming"
    who = threading.local()
    listed, copied = threading.Event(), threading.Event()  # listed: shard 0 looked at `.incoming`, or ended
    real_iterdir, real_run, real_copytree = Path.iterdir, subprocess.run, shutil.copytree

    def iterdir(self):
        if self == incoming and getattr(who, "step", None) == shard(0) and not listed.is_set():
            names = list(real_iterdir(self))
            listed.set()
            assert copied.wait(TIMEOUT), "shard 1 never published beside shard 0"
            return iter(names)
        return real_iterdir(self)

    def before_copy(dst) -> None:
        if getattr(who, "step", None) == shard(1) and not copied.is_set() and shard(1) in str(dst):
            assert listed.wait(TIMEOUT), "shard 0 never cleaned up or ended beside shard 1"

    def run_cmd(argv, *a, **k):
        before_copy(argv[-1] if argv and str(argv[0]).endswith("rsync") else "")
        try:
            return real_run(argv, *a, **k)
        finally:
            if getattr(who, "step", None) == shard(1) and listed.is_set():
                copied.set()

    def copytree(src, dst, *a, **k):
        before_copy(dst)
        try:
            return real_copytree(src, dst, *a, **k)
        finally:
            if getattr(who, "step", None) == shard(1) and listed.is_set():
                copied.set()

    def run_(prof_, root_, stage, step, cmd, **kw):
        who.step = step
        try:
            return runner(prof_, root_, stage, step, cmd, **kw)
        finally:
            who.step = None
            if step == shard(1):
                copied.set()  # whatever happened, shard 0 waits no longer
            if step == shard(0):
                listed.set()  # shard 0 has ended: shard 1 waits no longer

    def after_scan(point: str) -> None:
        if point == f"after:{SCAN}":
            raise Killed(point)

    monkeypatch.setattr(Path, "iterdir", iterdir)
    monkeypatch.setattr(runstep.subprocess, "run", run_cmd)
    monkeypatch.setattr(runstep.shutil, "copytree", copytree)
    with pytest.raises(Killed):
        bs.build_stage(prof, jobs_of(direct, root, 2, 2), run=run_, clock=lambda: NOON, load=idle, _hook=after_scan)
    assert listed.is_set() and copied.is_set()
    n = len(domains.read_plan(root, "S0", "test", "1")["shards"])
    assert shard_names(root) == [shard(k) for k in range(n)] and f"stage=S0/{SCAN}" in record_names(root)
    assert not list(root.rglob("*.publishing.json"))  # no publish was left half done
    assert not (root / ".incoming").exists()


def test_jobs_times_threads_above_the_hours_allowance_is_refused_before_any_step(direct, tmp_path, capsys):
    root = tmp_path / "root"
    prof = direct["prof"]  # day cap 4, night cap 6
    for when, load, jobs, cpus, allowance in ((NOON, 0.0, 3, 2, 4), (SUNDAY_NIGHT, 20.0, 5, 1, 4),
                                              (SUNDAY_NIGHT, 0.0, 7, 1, 6), (NOON, 0.0, 1, 5, 4)):
        b = jobs_of(direct, root, jobs, cpus)
        assert bs.build_stage(prof, b, run=runner, clock=lambda: when, load=lambda: load) == runstep.EXIT_REFUSED
        err = capsys.readouterr().err
        assert f"--scan-jobs {jobs}" in err and f"--scan-cpus {cpus}" in err and f"{jobs * cpus}" in err
        assert f"allows {allowance}" in err
    for jobs, cpus in ((0, 1), (1, 0)):
        assert bs.build_stage(prof, jobs_of(direct, root, jobs, cpus), run=runner, clock=lambda: NOON,
                              load=idle) == runstep.EXIT_REFUSED
        assert "at least 1" in capsys.readouterr().err
    # at night on an idle host 3 × 2 fits the night cap of 6: refused only for the missing release definition
    b = dataclasses.replace(jobs_of(direct, root, 3, 2), stages=None)
    assert bs.build_stage(prof, b, run=runner, clock=lambda: SUNDAY_NIGHT, load=idle) == runstep.EXIT_REFUSED
    assert "--scan-jobs" not in capsys.readouterr().err
    assert c08.tree(root) == {} and records(root) == {}
    # the command line refuses the same way (at any hour the test profile allows 4 or 6)
    argv = cli_argv(direct, root, profile_path(direct), "--scan-jobs", "4", "--scan-cpus", "2")
    assert cli.main(argv) == runstep.EXIT_REFUSED
    assert "--scan-jobs 4" in capsys.readouterr().err
    assert c08.tree(root) == {} and records(root) == {}


def test_the_defaults_use_the_hours_allowance_and_no_more(direct, tmp_path):
    b = build_of(direct, tmp_path / "root")
    assert bs.scan_slots(b, 40) == (10, 4)  # the day cap
    assert bs.scan_slots(b, 56) == (14, 4)  # the night cap
    assert bs.scan_slots(b, 6) == (1, 4)
    assert bs.scan_slots(b, 2) == (1, 2)  # fewer threads than the default: one job of all of them
    assert bs.scan_slots(dataclasses.replace(b, scan_cpus=8), 40) == (5, 8)
    assert bs.scan_slots(dataclasses.replace(b, scan_jobs=3), 40) == (3, 4)
    assert bs.scan_allowance(direct["prof"], NOON, 0.0) == 4
    assert bs.scan_allowance(direct["prof"], SUNDAY_NIGHT, 0.0) == 6
    assert bs.scan_allowance(direct["prof"], SUNDAY_NIGHT, 20.0) == 4  # busy host: the day cap


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


# --- long builds: the hour and lost scratch (card C-12) ----------------------------------------------------
#
# Each test continues a build from a shared fixture instead of building from nothing: `planned` is a build stopped
# after its scan's plan step (six or seven shards), `scanned` is `planned` continued through every shard at night and
# stopped before the merge. The oracle for the files is the direct (unsharded) build; for the plan and the shard
# records it is `scanned`, whose shards all ran at the night cap with the allowance never changing.

DAY = datetime(2026, 11, 9, 9, 0)  # Monday 09:00: the test profile's day cap (4) holds


class Clock:
    """An injected clock: `before` until it is armed, then one time of `script` per call (its last time for every
    call past its end). The runners below arm it while the scan's plan step runs, so the calls it counts are the ones
    build-stage makes after the plan: one before each shard launch, and one at each check while it holds a launch."""

    def __init__(self, before: datetime, script: list[datetime]):
        self.before, self.script, self.armed, self.calls = before, list(script), False, 0

    def __call__(self) -> datetime:
        if not self.armed:
            return self.before
        self.calls += 1
        return self.script[min(self.calls, len(self.script)) - 1]


def index(step: str) -> int:
    return int(step.rsplit("-", 1)[1])


class InFlight:
    """A runner that records, for each scan shard, the threads of the shards running when it starts (its own
    included) and the clock calls made by then, and runs every step through `run-step` at the time `when`. The shards
    numbered below `hold` stay running until build-stage holds a launch for the allowance (its `hold:` hook point,
    `hook`) or TIMEOUT passes: a build that launches past the allowance is then seen doing it, without timing."""

    def __init__(self, clock: Clock, when: datetime, hold: int = 0):
        self.clock, self.when, self.hold = clock, when, hold
        self.lock = threading.Lock()
        self.alive: dict[str, int] = {}
        self.starts: list[tuple[str, int, int]] = []  # (shard step, threads in flight at its start, clock calls)
        self.held = threading.Event()

    def hook(self, point: str) -> None:
        if point.startswith("hold:"):
            self.held.set()

    def __call__(self, prof, root, stage, step, cmd, **kw) -> int:
        def go() -> int:
            return runstep.run_step(prof, root, stage, step, cmd, clock=lambda: self.when, load=idle, **kw)

        if step.endswith(".plan"):
            self.clock.armed = True
        if ".shard-" not in step:
            return go()
        with self.lock:
            self.alive[step] = kw["threads"]
            self.starts.append((step, sum(self.alive.values()), self.clock.calls))
        try:
            if index(step) < self.hold:
                self.held.wait(TIMEOUT)
            return go()
        finally:
            with self.lock:
                del self.alive[step]


def copy_root(src: Path, base: Path) -> Path:
    shutil.copytree(src, base / "root", symlinks=True)
    return base / "root"


def scan_records(root: Path) -> dict[str, str]:
    return {k: v for k, v in c08.tree(root).items() if k.startswith(SHARD_RECORDS)}


def check_stage_files(direct, root: Path) -> None:
    """For a build stopped after its scan: every file but the shard records is the direct build's file at that path,
    byte for byte, and every domain file of the direct build is there."""
    got = {k: v for k, v in c08.tree(root).items() if not k.startswith(SHARD_RECORDS)}
    want = c08.tree(direct["root"])
    assert got == {k: want.get(k) for k in got}
    assert {k for k in want if k.startswith("l1/domains/")} <= set(got)


@pytest.fixture(scope="module")
def planned(direct, tmp_path_factory):
    """A build stopped after its scan's plan step, the scan cut into six or seven shards."""
    root = tmp_path_factory.mktemp("c12-planned") / "root"

    def stop(point: str) -> None:
        if point == f"after:{SCAN}.plan":
            raise Killed(point)

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(bs, "shard_residues", patched_residues(direct))
        with pytest.raises(Killed):
            bs.build_stage(small_budget(direct), build_of(direct, root), run=runner, clock=lambda: NOON, load=idle,
                           _hook=stop)
    return root


@pytest.fixture(scope="module")
def scanned(direct, planned, tmp_path_factory):
    """`planned` continued through every shard at night (four jobs of one thread) and stopped before the merge: the
    root with every shard's record, and the scratch directory with every shard's run."""
    base = tmp_path_factory.mktemp("c12-scanned")
    root = copy_root(planned, base)

    def until_merge(prof, root_, stage, step, cmd, **kw):
        if step == SCAN:
            raise Killed(step)
        return runner(prof, root_, stage, step, cmd, **kw)

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(bs, "shard_residues", patched_residues(direct))
        b = dataclasses.replace(jobs_of(direct, root, 4, 1), tmp=base / "scratch")
        with pytest.raises(Killed):
            bs.build_stage(small_budget(direct), b, run=until_merge, clock=lambda: SUNDAY_NIGHT, load=idle)
    n = len(domains.read_plan(root, "S0", "test", "1")["shards"])
    assert n >= 5 and shard_names(root) == [shard(k) for k in range(n)]
    return root, base / "scratch"


@pytest.mark.covers("D-01")
def test_c12_after_the_clock_moves_from_night_to_day_no_shard_launch_passes_the_day_cap(
        direct, planned, scanned, tmp_path, monkeypatch):
    """Six jobs of one thread fit the night cap (6) when the build starts. The clock reads night for the first three
    shard launches, then day (cap 4). After the change no shard may start while it would put more than 4 threads in
    flight; the shards started at night run on and finish. Counted through the runner. The build ends byte-identical
    to the direct build, with the shard records of `scanned`, whose shards all ran at night."""
    monkeypatch.setattr(bs, "shard_residues", patched_residues(direct))
    monkeypatch.setattr(bs, "RECHECK_S", 1.0, raising=False)
    prof, root = small_budget(direct), copy_root(planned, tmp_path)
    clock = Clock(SUNDAY_NIGHT, [SUNDAY_NIGHT] * 3 + [DAY])
    count = InFlight(clock, SUNDAY_NIGHT, hold=3)
    b = dataclasses.replace(jobs_of(direct, root, 6, 1), tmp=tmp_path / "scratch")
    assert bs.build_stage(prof, b, run=count, clock=clock, load=idle, _hook=count.hook) == 0
    n = len(domains.read_plan(root, "S0", "test", "1")["shards"])
    assert sorted(s for s, _, _ in count.starts) == [shard(k) for k in range(n)]  # each shard started once
    after = {s: t for s, t, _ in count.starts if index(s) >= 3}  # the shards launched after the change
    assert after and all(t <= prof.max_threads_day for t in after.values()), count.starts
    assert count.held.is_set() and clock.calls > 3  # a launch was held, and the clock moved
    check_sharded_equals_direct(direct, root)
    assert scan_records(root) == scan_records(scanned[0])


def test_c12_with_the_allowance_below_one_shards_threads_the_build_waits_then_proceeds(
        direct, planned, scanned, tmp_path, monkeypatch, capsys):
    """One job of 5 threads fits the night cap (6) when the build starts. The clock then reads day (cap 4) for two
    checks and night after. No shard starts while it reads day: build-stage holds the launch, waits RECHECK_S (at
    least a second) on the shards it runs, of which there are none, checks again, and launches when the clock
    allows."""
    monkeypatch.setattr(bs, "shard_residues", patched_residues(direct))
    monkeypatch.setattr(bs, "RECHECK_S", 1.0, raising=False)
    prof, root = small_budget(direct), copy_root(planned, tmp_path)
    clock = Clock(SUNDAY_NIGHT, [DAY, DAY, SUNDAY_NIGHT])
    count = InFlight(clock, SUNDAY_NIGHT)
    seen: list[str] = []

    def hook(point: str) -> None:
        seen.append(point)
        if point == f"after:{SCAN}":
            raise Killed(point)

    b = dataclasses.replace(jobs_of(direct, root, 1, 5), tmp=tmp_path / "scratch")
    with pytest.raises(Killed):
        bs.build_stage(prof, b, run=count, clock=clock, load=idle, _hook=hook)
    assert [p for p in seen if p.startswith("hold:")] == [f"hold:{shard(0)}"] * 2
    assert f"allows {prof.max_threads_day}" in capsys.readouterr().err
    assert [c for s, _, c in count.starts if s == shard(0)] == [3]  # launched at the third check, the first at night
    assert {t for _, t, _ in count.starts} == {5}
    assert bs.RECHECK_S >= 1
    check_stage_files(direct, root)
    assert scan_records(root) == scan_records(scanned[0])


def lost_scratch(scanned, tmp_path: Path) -> tuple[Path, Path]:
    """Copies of `scanned`'s root and scratch directory, for a test to damage."""
    root = copy_root(scanned[0], tmp_path)
    shutil.copytree(scanned[1], tmp_path / "scratch")
    return root, tmp_path / "scratch"


def shard_run(root: Path, scratch: Path, k: int) -> Path:
    return domains.run_path(scratch, root, "S0", "test", "1", k)


def shard_record_path(root: Path, k: int) -> Path:
    return root / domains.scan_dir("S0", "test", "1") / f"shard-{k:05d}.json"


@pytest.mark.covers("D-01", "D-03")
@pytest.mark.parametrize("damage", ["deleted", "corrupted"])
def test_c12_a_shard_whose_scratch_run_was_lost_is_scanned_again_and_the_build_ends_byte_identical(
        direct, scanned, tmp_path, monkeypatch, capsys, damage):
    """Scratch is damaged between the shards and their merge (a reboot clears /tmp): shard 2's run deleted, or one
    of its bytes changed (same size: only the sha256 its record holds can tell). build-stage scans that shard again
    into scratch, and no other (counted by the `after:<shard>.rescan` hook points), without publishing its record
    again, then merges. The deleted case runs to the end and equals the direct build; the corrupted case stops
    after the scan step and its files equal the direct build's."""
    monkeypatch.setattr(bs, "shard_residues", patched_residues(direct))
    prof = small_budget(direct)
    root, scratch = lost_scratch(scanned, tmp_path)
    run2 = shard_run(root, scratch, 2)
    if damage == "deleted":
        run2.unlink()
    else:
        data = bytearray(run2.read_bytes())
        data[len(data) // 2] ^= 0xFF
        run2.write_bytes(bytes(data))
    rec2 = shard_record_path(root, 2).read_bytes()
    before = records(root)
    seen: list[str] = []

    def hook(point: str) -> None:
        seen.append(point)
        if point == f"after:{SCAN}" and damage == "corrupted":
            raise Killed(point)

    b = dataclasses.replace(jobs_of(direct, root, 4, 1), tmp=scratch)
    if damage == "deleted":
        assert bs.build_stage(prof, b, run=runner, clock=lambda: NOON, load=idle, _hook=hook) == 0
        check_sharded_equals_direct(direct, root)
    else:
        with pytest.raises(Killed):
            bs.build_stage(prof, b, run=runner, clock=lambda: NOON, load=idle, _hook=hook)
        check_stage_files(direct, root)
    assert [p for p in seen if p.endswith(".rescan")] == [f"after:{shard(2)}.rescan"]  # that shard, and only it
    assert shard(2) in capsys.readouterr().err
    assert shard_record_path(root, 2).read_bytes() == rec2  # the record stayed
    assert scan_records(root) == scan_records(scanned[0])
    assert {k: v for k, v in records(root).items() if k in before} == before  # no shard step ran or published again


def test_c12_a_rescan_that_gives_another_checksum_stops_with_an_error_naming_the_shard(
        direct, scanned, tmp_path, monkeypatch, capfd):
    """Shard 1's record is made to hold a run sha256 that its inputs do not give. Its run in scratch fails the
    record, so build-stage scans it again; the new run does not give the record's sha256 either, so the build stops
    before the merge, with an error that names the shard, and leaves the run in scratch as it was."""
    monkeypatch.setattr(bs, "shard_residues", patched_residues(direct))
    prof = small_budget(direct)
    root, scratch = lost_scratch(scanned, tmp_path)
    path = shard_record_path(root, 1)
    rec = json.loads(path.read_text())
    rec["run_sha256"] = "0" * 64
    path.write_text(json.dumps(rec, sort_keys=True, indent=1) + "\n")
    run1 = shard_run(root, scratch, 1).read_bytes()
    b = dataclasses.replace(jobs_of(direct, root, 4, 1), tmp=scratch)
    code = bs.build_stage(prof, b, run=runner, clock=lambda: NOON, load=idle)
    assert code not in (0, runstep.EXIT_STOP)
    err = capfd.readouterr().err
    assert f"{shard(1)}.rescan exited {code}" in err  # build-stage names the shard step
    assert "shard 1 of domains-test-1" in err and "0" * 64 in err  # the rescan names the shard and the checksums
    assert f"stage=S0/{SCAN}.json" not in records(root)  # no merge
    assert not (root / "l1/domains").exists()
    assert shard_run(root, scratch, 1).read_bytes() == run1


def test_c12_a_record_without_a_run_checksum_is_scanned_again_and_accepted_with_a_note(
        direct, scanned, tmp_path, monkeypatch, capsys):
    """A shard record that holds no run sha256 (`domains.scan_shard` has written one into every record since card
    B-24; the case is handled for a record written otherwise) cannot be checked against its run: build-stage scans
    the shard again, accepts the new run when its domain and hit counts equal the record's, says so in the log, and
    merges."""
    monkeypatch.setattr(bs, "shard_residues", patched_residues(direct))
    prof = small_budget(direct)
    root, scratch = lost_scratch(scanned, tmp_path)
    path = shard_record_path(root, 3)
    rec = json.loads(path.read_text())
    del rec["run_sha256"]
    path.write_text(json.dumps(rec, sort_keys=True, indent=1) + "\n")
    seen: list[str] = []

    def hook(point: str) -> None:
        seen.append(point)
        if point == f"after:{SCAN}":
            raise Killed(point)

    b = dataclasses.replace(jobs_of(direct, root, 4, 1), tmp=scratch)
    with pytest.raises(Killed):
        bs.build_stage(prof, b, run=runner, clock=lambda: NOON, load=idle, _hook=hook)
    assert [p for p in seen if p.endswith(".rescan")] == [f"after:{shard(3)}.rescan"]
    err = capsys.readouterr().err
    assert shard(3) in err and "no run sha256" in err
    check_stage_files(direct, root)


# --- defects of the first real build (card C-13, evidence planning/measurements/2026-11-10_build_stage_pod.md) -----


def snap_after(direct, name: str) -> dict[str, str]:
    """The direct build's root as it stood after step `name`."""
    return c08.tree(direct["snaps"][c08.NAMES.index(name)])


@pytest.mark.covers("C-13", "D-03")
def test_c13_each_step_has_its_own_scratch_and_its_next_attempt_removes_what_a_killed_attempt_left(direct, tmp_path):
    """C-10: a clustering step stopped by `systemctl --user stop` left 1.2 GB in the shared scratch directory, which
    nothing removed (*measured*). Each step that takes `--tmp` now has its own directory, `<scratch>/<stage>.<step>`;
    the next attempt of that step removes it before the step starts, and a step that succeeds leaves none. Here the
    first attempt of ingest is killed (it leaves a builder's scratch directory and exits 143, as run-step does after
    a SIGTERM); the second finds its directory empty and publishes."""
    root, prof = tmp_path / "root", direct["prof"]
    b = build_of(direct, root)
    scratch = bs.scratch_of(prof, b)
    with_tmp = [s for s in bs.steps(prof, b) if s.argv and "--tmp" in s.argv]
    assert {s.name for s in with_tmp} >= {"ingest", "genecall", "igr", "tapes", "cluster", "tokens", "postings", "dim"}
    for s in with_tmp:
        assert s.argv[s.argv.index("--tmp") + 1] == str(scratch / f"{s.run_stage}.{s.name}"), s.name
    mine = scratch / "S0.ingest"
    left = mine / "biobox-ingest-killed"
    attempts: list[str] = []

    def killed_once(prof_, root_, stage, step, cmd, **kw):
        if step == "ingest":
            attempts.append(step)
            if len(attempts) == 1:
                left.mkdir(parents=True)
                (left / "pass1.parquet").write_bytes(b"x" * 4096)  # what a builder killed mid-step leaves
                return 143
            assert mine.is_dir() and not any(mine.iterdir()), "the earlier attempt's scratch is still there"
        return runner(prof_, root_, stage, step, cmd, **kw)

    def stop(point: str) -> None:
        if point == "after:ingest":
            raise Killed(point)

    assert bs.build_stage(prof, b, run=killed_once) == 143
    assert left.is_dir()  # left until the step's next attempt (or by hand)
    with pytest.raises(Killed):
        bs.build_stage(prof, b, run=killed_once, _hook=stop)
    assert attempts == ["ingest", "ingest"]
    assert not mine.exists()  # the step succeeded: its scratch is gone
    assert record_names(root) == {"stage=S0/inputs", "stage=S0/ingest"}


@pytest.mark.covers("C-13", "D-03")
def test_c13_a_stop_inside_a_step_discards_its_finished_work_and_the_next_call_runs_it_again(direct, tmp_path, capsys):
    """C-10 found no fixture test with `STOP` created *during* a step (the one above creates it between steps). The
    documented rule (card C-02, step 6): the step's command finishes, run-step publishes nothing of it, leaves its
    staging and exits 75; build-stage stops there; the next call, `STOP` removed, runs that step again from the start
    (on the pod: S1a's gene calling, 57.7 s and 4.98 CPU-minutes, ran twice). Here genecall's command creates `STOP`
    as it starts and then runs to the end."""
    root, prof = tmp_path / "root", direct["prof"]
    stop_then = ("import pathlib, subprocess, sys; pathlib.Path(sys.argv[1]).write_text(''); "
                 "sys.exit(subprocess.call(sys.argv[2:]))")

    def stop_inside_genecall(prof_, root_, stage, step, cmd, **kw):
        if step == "genecall":
            cmd = [sys.executable, "-c", stop_then, str(Path(root_) / "STOP"), *cmd]
        return runner(prof_, root_, stage, step, cmd, **kw)

    assert bs.build_stage(prof, build_of(direct, root), run=stop_inside_genecall) == runstep.EXIT_STOP
    assert "STOP appeared; not publishing S0.genecall" in capsys.readouterr().err
    assert record_names(root) == {"stage=S0/inputs", "stage=S0/ingest"}
    assert {k: v for k, v in c08.tree(root).items() if k != "STOP"} == snap_after(direct, "ingest")
    staged = Path(prof.tmp) / "staging" / "S0.genecall"
    assert any(p.is_file() for p in staged.rglob("*")), "the finished work is left in staging, not published"
    (root / "STOP").unlink()
    started: list[str] = []

    def recording(prof_, root_, stage, step, cmd, **kw):
        if not (Path(root_) / "runs" / f"stage={stage}" / f"{step}.json").exists():
            started.append(step)
        return runner(prof_, root_, stage, step, cmd, **kw)

    def after_genecall(point: str) -> None:
        if point == "after:genecall":
            raise Killed(point)

    with pytest.raises(Killed):
        bs.build_stage(prof, build_of(direct, root), run=recording, _hook=after_genecall)
    assert started == ["genecall"]  # inputs and ingest skipped; genecall run again from the start
    assert c08.tree(root) == snap_after(direct, "genecall")


def run_cli(*argv: str) -> int:
    return cli.main([str(a) for a in argv])


def test_c23_a_built_release_names_the_code_that_built_it(built, tmp_path, capsys):
    """Card C-23: `release verify` prints one `built by:` line (information, exit code unchanged); a copy whose run
    records predate the card says so, and still verifies; the catalog's `built_by` is the same JSON."""
    import json

    import duckdb

    from biobox import codeid

    me = codeid.code_identity()["code_id"]
    assert run_cli("release", "verify", "--root", built, "--release", REL) == 0
    lines = [x for x in capsys.readouterr().out.splitlines() if x.startswith("built by:")]
    assert len(lines) == 1 and me[:12] in lines[0]
    root = tmp_path / "old"
    shutil.copytree(built, root, symlinks=True)
    for p in (root / "runs").rglob("*.json"):  # records as written before the card
        rec = json.loads(p.read_text())
        for k in ("code_id", "code_version", "code_commit"):
            rec.pop(k, None)
        p.write_text(json.dumps(rec))
    assert run_cli("release", "verify", "--root", root, "--release", REL) == 0
    lines = [x for x in capsys.readouterr().out.splitlines() if x.startswith("built by:")]
    assert lines == ["built by: unknown: built before card C-23"]
    con = duckdb.connect(str(catalog.catalog_path(built, REL)), read_only=True)
    try:
        info = dict(con.execute("SELECT key, value FROM catalog_info").fetchall())
    finally:
        con.close()
    got = json.loads(info["built_by"])
    assert [c["code_id"] for c in got["code"]] == [me] and got["code"][0]["steps"] > 3
    assert got["unknown_steps"] == 0


@pytest.mark.covers("C-13", "LY-10")
def test_c13_release_verify_passes_a_built_release_and_fails_a_damaged_copy(built, tmp_path, capsys):
    assert run_cli("release", "verify", "--root", built, "--release", REL) == 0
    out = capsys.readouterr().out
    assert f"release {REL}: verified" in out
    root = tmp_path / "root"
    shutil.copytree(built, root, symlinks=True)
    victim = catalog.read_manifest(catalog.manifest_path(root, REL))["files"][0]["path"]
    data = bytearray((root / victim).read_bytes())
    data[len(data) // 2] ^= 0xFF
    (root / victim).write_bytes(bytes(data))
    assert run_cli("release", "verify", "--root", root, "--release", REL) == 1
    assert f"sha256 differs: {victim}" in capsys.readouterr().out
    assert run_cli("release", "verify", "--root", root, "--release", "rNo.1") == 1  # no manifest: a failure
    assert "manifests/release=rNo.1.json" in capsys.readouterr().out


@pytest.mark.covers("C-13", "R-09")
def test_c13_release_verify_checks_a_later_release_against_a_named_previous_one(built, tmp_path, capsys):
    """`--previous P`: every file of P's manifest outside P's own `release=` partition is in the release's manifest
    with the same sha256 (`contracts.check_append_only`). The later release here is a second release over the same
    stage, defined and given its manifest directly; then P's manifest is made to list one stage file with another
    sha256, as a rewritten file would show."""
    root = tmp_path / "root"
    shutil.copytree(built, root, symlinks=True)
    later = "rS0.c13"
    catalog.define_release(root, later, ["S0"])
    catalog.write_manifest(catalog.build_manifest(root, later, ["S0"]))
    argv = ("release", "verify", "--root", root, "--release", later, "--previous", REL)
    assert run_cli(*argv) == 0
    out = capsys.readouterr().out
    assert f"release {later}: verified" in out and f"append-only against {REL}: ok" in out
    path = catalog.manifest_path(root, REL)
    man = catalog.read_manifest(path)
    f = next(f for f in man["files"] if f["partitions"].get("stage") == "S0")
    f["sha256"] = "0" * 64
    path.write_text(json.dumps(man, sort_keys=True, indent=1) + "\n")
    assert run_cli(*argv) == 1
    assert f"{f['path']}: rewritten" in capsys.readouterr().out
    assert run_cli("release", "verify", "--root", root, "--release", later, "--previous", "rNo.1") == 1


# --- clustering's own cap (card B-36, DEC-25) ---------------------------------------------------------------------

STAGE_STEPS = ["inputs", "ingest", "genecall", "igr", "tapes", SCAN]


def fixture_proteins(root: Path, stage: str) -> int:
    """Rows of the stage's `protein_seq` files, counted here from the files themselves."""
    files = sorted((Path(root) / f"l0seq/protein_seq/gene_call_version={V}/stage={stage}").rglob("*.parquet"))
    return sum(pq.read_metadata(f).num_rows for f in files)


class Calls:
    """A runner that records each step it is asked to run (and its keywords), then runs it as `runner` does."""

    def __init__(self):
        self.seen: list[tuple[str, dict]] = []

    def __call__(self, prof, root, stage, step, cmd, **kw) -> int:
        self.seen.append((step, dict(kw)))
        return runner(prof, root, stage, step, cmd, **kw)


@pytest.fixture(scope="module")
def capped(direct, tmp_path_factory):
    """A build on an empty root whose profile caps clustering at 0 CPU-minutes: the stage's proteins are known only
    after its gene calling, so the stage steps run and the cluster step is refused before its command."""
    root = tmp_path_factory.mktemp("b36-capped") / "root"
    calls = Calls()
    prof = dataclasses.replace(direct["prof"], cluster_cpu_min_max=0)
    code = bs.build_stage(prof, build_of(direct, root), run=calls)
    return root, code, calls


@pytest.mark.covers("B-36")
def test_b36_a_release_over_the_cap_is_refused_before_the_cluster_steps_command(direct, capped):
    root, code, calls = capped
    assert code == runstep.EXIT_REFUSED
    assert [s for s, _ in calls.seen] == STAGE_STEPS  # the cluster step was never handed to run-step
    assert record_names(root) == {f"stage=S0/{n}" for n in STAGE_STEPS}
    assert not any(k.startswith("l2/") or f"release={REL}" in k for k in c08.tree(root))


@pytest.mark.covers("B-36")
def test_b36_once_every_stage_is_gene_called_a_release_over_the_cap_is_refused_before_any_step(
        direct, capped, tmp_path, capsys):
    root = tmp_path / "root"
    shutil.copytree(capped[0], root, symlinks=True)
    before, files = records(root), c08.tree(root)
    n = fixture_proteins(root, "S0")
    assert n > 0
    calls, seen = Calls(), []
    prof = dataclasses.replace(direct["prof"], cluster_cpu_min_max=0)
    assert bs.build_stage(prof, build_of(direct, root), run=calls, _hook=seen.append) == runstep.EXIT_REFUSED
    err = capsys.readouterr().err
    assert "refused before any step" in err and "cluster_cpu_min_max" in err and f"{n:,} proteins" in err
    assert calls.seen == [] and seen == []
    assert records(root) == before and c08.tree(root) == files


@pytest.mark.covers("B-36")
def test_b36_without_the_cap_clustering_over_cpu_min_per_process_is_refused(direct, capped, tmp_path, capsys):
    """No `cluster_cpu_min_max`: no exemption. A budget of 0 puts the fixture's clustering over it; ingest and the
    scan have their run records, so neither is estimated again and the refusal is clustering's."""
    root = tmp_path / "root"
    shutil.copytree(capped[0], root, symlinks=True)
    before = records(root)
    calls = Calls()
    prof = dataclasses.replace(direct["prof"], cpu_min_per_process=0, cluster_cpu_min_max=None)
    assert bs.build_stage(prof, build_of(direct, root), run=calls) == runstep.EXIT_REFUSED
    err = capsys.readouterr().err
    assert "cpu_min_per_process" in err and "no cluster_cpu_min_max" in err and "DEC-25" in err
    assert calls.seen == [] and records(root) == before


@pytest.mark.covers("B-36")
def test_b36_a_release_under_the_cap_runs_and_its_cluster_record_names_dec25(direct, capped, tmp_path):
    root = tmp_path / "root"
    shutil.copytree(capped[0], root, symlinks=True)
    calls = Calls()
    prof = dataclasses.replace(direct["prof"], cluster_cpu_min_max=12000)
    assert bs.build_stage(prof, build_of(direct, root), run=calls) == 0
    assert [s for s, _ in calls.seen] == STAGE_STEPS + NAMES[7:]  # the stage steps are resumed (skipped) by run-step
    assert [s for s, kw in calls.seen if "budget_reason" in kw] == ["cluster"]
    assert dict(calls.seen)["cluster"]["budget_reason"] == "exempt: DEC-25"
    rec = json.loads((root / "runs" / f"stage={REL}" / "cluster.json").read_bytes())
    assert rec["cpu_budget_reason"] == "exempt: DEC-25"
    assert rec["over_cpu_budget"] is False  # kept: the fixture's clustering is far under 300 CPU-minutes
    others = {k: json.loads(v) for k, v in records(root).items() if k != f"stage={REL}/cluster.json"}
    assert others and not any("cpu_budget_reason" in r for r in others.values())  # only clustering is exempt
    check_equals_direct(direct, root)


@pytest.mark.covers("B-36")
def test_b36_a_profile_without_the_cap_writes_no_reason(built):
    """The default build of this file (a profile without `cluster_cpu_min_max`): no record carries a reason."""
    assert not any("cpu_budget_reason" in json.loads(v) for v in records(built).values())


# --- clustering in parts (card B-41) --------------------------------------------------------------------------------

PLAN_RECORDS = "l2/cluster_plan/"  # the split clustering's records: a build's own, never in a manifest


def check_parts_equal_direct(direct, root: Path) -> None:
    """Every file but the split clustering's records is the direct (one-process) build's, byte for byte; so are the
    manifest and the catalog, which never list those records."""
    got = {k: v for k, v in c08.tree(root).items() if not k.startswith(PLAN_RECORDS)}
    want = c08.tree(direct["root"])
    assert sorted(got) == sorted(want)
    assert c08.data(got) == c08.data(want)
    man = c08.check_manifest(root, REL, direct["manifest"])
    assert catalog.verify_manifest(man) == []
    c08.check_catalog(root, REL, man, Path(direct["prof"].tmp) / "staging")


def release_steps_done(root: Path) -> list[str]:
    return sorted(n.removeprefix(f"stage={REL}/") for n in record_names(root) if n.startswith(f"stage={REL}/"))


def parts_of(direct, root: Path) -> bs.Build:
    return dataclasses.replace(build_of(direct, root), cluster_parts=2)


@pytest.fixture(scope="module")
def in_parts(direct, tmp_path_factory):
    """One build-stage call with every per-query call of clustering forced into 2 parts, and its hook points."""
    root = tmp_path_factory.mktemp("b41-parts") / "root"
    seen: list[str] = []
    assert bs.build_stage(direct["prof"], parts_of(direct, root), run=runner, _hook=seen.append) == 0
    return root, seen


@pytest.mark.covers("B-37", "MR-16")
def test_b41_a_release_clustered_in_parts_equals_the_build_that_clusters_in_one_process(direct, in_parts):
    root, seen = in_parts
    check_parts_equal_direct(direct, root)
    done = release_steps_done(root)
    calls = l2cluster.calls()
    want = {l2cluster.PLAN_STEP, l2cluster.FINAL_STEP, *(c.name for c in calls),
            *(f"{c.name}.cut" for c in calls if c.split)}
    parts = [n for n in done if ".part-" in n]
    assert set(done) == want | set(parts) | set(NAMES[8:])
    # every per-query call with two keys or more ran as two part steps; each part is a run-step of its own
    by_call = {c.name: [p for p in parts if p.startswith(f"{c.name}.part-")] for c in calls if c.split}
    assert all(len(v) in (0, 2) for v in by_call.values()) and sum(len(v) == 2 for v in by_call.values()) >= 12
    # each step but the tables left its record in the cluster_plan artifact; none is in the manifest (checked above)
    recs = {p.name.removesuffix(".json") for p in (root / PLAN_RECORDS / f"release={REL}").glob("*.json")}
    assert recs == set(done) - {l2cluster.FINAL_STEP} - set(NAMES[8:])
    # in order: the plan, each call (its parts before its join), the tables under the cluster step's name
    order = [p.removeprefix("after:") for p in seen]
    i, j = order.index(l2cluster.PLAN_STEP), order.index(l2cluster.FINAL_STEP)
    assert order[j - 1] == calls[-1].name and set(order[i:j]) == recs
    for c in calls:
        k = order.index(c.name)
        assert all(order.index(p) < k for p in by_call.get(c.name, []))
    # the work directory and the parts are gone
    assert not l2cluster.work_dir(bs.scratch_of(direct["prof"], parts_of(direct, root)), root.resolve(), REL).exists()


def test_b41_a_second_call_after_clustering_in_parts_starts_nothing(direct, in_parts, tmp_path, monkeypatch):
    root = tmp_path / "root"
    shutil.copytree(in_parts[0], root, symlinks=True)
    before, files = records(root), c08.tree(root)

    def no_command(*a, **k):
        raise AssertionError("a step command was started on the second call")

    monkeypatch.setattr(runstep, "_execute", no_command)
    assert bs.build_stage(direct["prof"], parts_of(direct, root), run=runner) == 0
    assert records(root) == before and c08.tree(root) == files


@pytest.mark.covers("B-37")
def test_b41_clustering_runs_in_parts_when_its_estimate_passes_80_percent_of_the_budget(direct, capped, tmp_path):
    """The routing (`_cluster_split`): in one process while the estimate fits 80% of the budget (the default, as an
    unsharded scan), in parts above that, when forced, or when a plan is already published."""
    root = capped[0]
    b = build_of(direct, root)
    est = bs.cluster_cpu_min(fixture_proteins(root, "S0"))
    assert est > 0
    over = dataclasses.replace(direct["prof"], cpu_min_per_process=est / 0.8 * 0.999)
    under = dataclasses.replace(direct["prof"], cpu_min_per_process=est / 0.8 * 1.001)
    assert bs._cluster_split(over, b, root) is True and bs._cluster_split(under, b, root) is False
    assert bs._cluster_split(direct["prof"], b, root) is False  # the fixture at 300 CPU-minutes: one process
    assert bs._cluster_split(direct["prof"], dataclasses.replace(b, cluster_parts=1), root) is True
    copy = tmp_path / "root"
    shutil.copytree(root, copy, symlinks=True)
    plan = copy / PLAN_RECORDS / f"release={REL}" / f"{l2cluster.PLAN_STEP}.json"
    plan.parent.mkdir(parents=True)
    plan.write_text("{}")
    assert bs._cluster_split(direct["prof"], build_of(direct, copy), copy) is True
    # the plan step's cut: the forced count, else 80% of the budget per part
    argv = bs.cluster_steps(direct["prof"], b)("plan").argv
    assert argv[argv.index("--part-cpu-min") + 1] == repr(0.8 * 300)
    argv = bs.cluster_steps(direct["prof"], dataclasses.replace(b, cluster_parts=4))("plan").argv
    assert argv[argv.index("--parts") + 1] == "4" and "--part-cpu-min" not in argv


@pytest.mark.covers("B-37", "D-03")
def test_b41_a_scratch_database_lost_after_its_record_is_made_again_and_the_build_ends_byte_identical(direct, tmp_path):
    """A build killed after the NR50 cascade's last prefilter; then the joined NR90 cascade prefilter is deleted from
    scratch (its parts were removed when the join was published). The next call makes the parts and the join again
    (each to its record), then goes on from the step after the kill, and ends byte-identical."""
    root = tmp_path / "root"
    b = parts_of(direct, root)
    pref90 = next(c for c in l2cluster.calls() if c.name.startswith("cluster.nr90.") and c.module == "prefilter")
    last50 = [c for c in l2cluster.calls() if c.name.startswith("cluster.nr50.") and c.module == "prefilter"][-1]

    def kill(point: str) -> None:
        if point == f"after:{last50.name}":
            raise Killed(point)

    with pytest.raises(Killed):
        bs.build_stage(direct["prof"], b, run=runner, _hook=kill)
    w = l2cluster.work_dir(bs.scratch_of(direct["prof"], b), root.resolve(), REL)
    gone = list((w / pref90.output).parent.glob(Path(pref90.output).name + "*"))
    assert gone
    for p in gone:
        p.unlink()
    parts = sorted(n for n in release_steps_done(root) if n.startswith(f"{pref90.name}.part-"))
    assert len(parts) == 2
    started: list[str] = []

    def recording(prof, root_, stage, step, cmd, **kw):
        if not (Path(root_) / "runs" / f"stage={stage}" / f"{step}.json").exists():
            started.append(step)
        return runner(prof, root_, stage, step, cmd, **kw)

    seen: list[str] = []
    assert bs.build_stage(direct["prof"], b, run=recording, _hook=seen.append) == 0
    redone = [p.removeprefix("after:").removesuffix(".redo") for p in seen if p.endswith(".redo")]
    assert redone == [*parts, pref90.name]
    after = l2cluster.calls()[[c.name for c in l2cluster.calls()].index(last50.name) + 1]
    assert started[0] == (f"{after.name}.cut" if after.split else after.name)  # nothing published ran again
    check_parts_equal_direct(direct, root)


# --- one process for a clustering call that cannot be split (card B-45, DEC-29) -------------------------------------
#
# The forced per-part table (`rebuild_cpu_min`) is a plan parameter that build-stage never passes, so these builds
# start from a plan step run through run-step by hand (`plan_forced`, `cluster.plan` with the forced numbers, as the
# plan step would publish it); build-stage then finds the plan published, clusters in parts and makes every cut from
# that plan. The numbers are relative to the fixture's own one-process estimate T(n), so the cap lies under T(n) (the
# release in one process is refused) and over the one-process calls' estimates.

FAST = ["cluster.nr90.11.prefilter", "cluster.nr50.11.prefilter"]  # the prefilters whose parts rebuild their table

_PLAN_FORCED = (
    "import os, sys\nfrom pathlib import Path\nfrom biobox.l2 import cluster\na = sys.argv[1:]\n"
    "cluster.plan(Path(a[0]), a[1], a[2], [a[3]], Path(a[4]), part_cpu_min=float(a[5]), "
    "one_process_cpu_min_max=float(a[6]), rebuild_cpu_min=float(a[7]), out=Path(os.environ['BIOBOX_STAGING']))\n")


def forced_numbers(root: Path) -> tuple[int, float, float, float]:
    """(proteins, T(n), part budget, forced table): a part budget twice the shared index's estimate (so the sensitive
    prefilters' index fits it) and a forced table half as large again (so it fills the budget)."""
    n = fixture_proteins(root, "S0")
    part = 2 * l2cluster.index_cpu_min(n)
    return n, bs.cluster_cpu_min(n), part, 1.5 * part


def plan_forced(direct, root: Path, plan_cap: float) -> None:
    b = build_of(direct, root)
    scratch = bs.scratch_of(direct["prof"], b)
    scratch.mkdir(parents=True, exist_ok=True)
    _, _, part, rebuild = forced_numbers(root)
    cmd = [sys.executable, "-c", _PLAN_FORCED, str(Path(root).resolve()), REL, V, "S0", str(scratch), repr(part),
           repr(plan_cap), repr(rebuild)]
    assert runner(direct["prof"], root, REL, l2cluster.PLAN_STEP, cmd) == 0


@pytest.mark.covers("B-45", "B-36")
def test_b45_the_plan_step_carries_the_profiles_one_process_cap(direct, tmp_path):
    b = build_of(direct, tmp_path / "root")
    assert direct["prof"].cluster_cpu_min_max is None
    assert "--one-process-cpu-min-max" not in bs.cluster_steps(direct["prof"], b)("plan").argv
    for parts in (None, 4):
        prof = dataclasses.replace(direct["prof"], cluster_cpu_min_max=20000)
        argv = bs.cluster_steps(prof, dataclasses.replace(b, cluster_parts=parts))("plan").argv
        assert argv[argv.index("--one-process-cpu-min-max") + 1] == repr(20000.0)
        # the other steps carry nothing of it: the cut reads it from the plan
        assert "--one-process-cpu-min-max" not in bs.cluster_steps(prof, b)("step", f"{FAST[0]}.cut").argv


@pytest.mark.covers("B-45")
def test_b45_release_cluster_plan_records_the_cap_and_shows_the_one_process_calls(direct, capped, tmp_path, capsys):
    root, out, tmp = capped[0], tmp_path / "out", tmp_path / "tmp"
    assert cli.main(["release", "cluster", "plan", "--root", str(root), "--gene-call-version", V, "--release", REL,
                     "--stages", "S0", "--tmp", str(tmp), "--part-cpu-min", "240", "--one-process-cpu-min-max",
                     "20000", "--out", str(out)]) == 0
    res = json.loads(capsys.readouterr().out)["result"]
    assert res["cut"] == {"part_cpu_min": 240.0, "one_process_cpu_min_max": 20000.0}
    assert sorted(res["one_process"]["calls"]) == sorted(FAST)
    assert not any(v["possible"] for v in res["one_process"]["calls"].values())  # the fixture is far too small
    assert "DEC-29" in res["one_process"]["why"] and res["one_process"]["cpu_min_max"] == 20000.0
    assert l2cluster.read_record(root, REL, l2cluster.PLAN_STEP, out)["cut"] == res["cut"]


@pytest.mark.covers("B-45", "B-36")
def test_b45_in_parts_a_release_over_the_whole_release_cap_runs_its_unsplittable_calls_as_one_process(
        direct, capped, tmp_path):
    root = tmp_path / "root"
    shutil.copytree(capped[0], root, symlinks=True)
    n, whole, part, rebuild = forced_numbers(root)
    cap = (rebuild + whole) / 2
    assert rebuild * 1.05 < cap < whole
    prof = dataclasses.replace(direct["prof"], cluster_cpu_min_max=cap)
    assert bs.cluster_refusal(prof, REL, n) is not None  # in one process this release is refused: T(n) > cap
    plan_forced(direct, root, cap)
    calls = Calls()
    assert bs.build_stage(prof, build_of(direct, root), run=calls) == 0
    exempt = {s: kw["budget_reason"] for s, kw in calls.seen if "budget_reason" in kw}
    assert exempt == {name: "exempt: DEC-29" for name in FAST}
    for name in FAST:
        cut = l2cluster.read_record(root, REL, f"{name}.cut")
        assert len(cut["parts"]) == 1 and cut["one_process"]["cpu_min"] <= cap
        rec = json.loads((root / "runs" / f"stage={REL}" / f"{name}.json").read_bytes())
        assert rec["cpu_budget_reason"] == "exempt: DEC-29"
    others = {k: json.loads(v) for k, v in records(root).items()
              if k not in {f"stage={REL}/{name}.json" for name in FAST}}
    assert not any("cpu_budget_reason" in r for r in others.values())
    check_parts_equal_direct(direct, root)


@pytest.mark.covers("B-45", "B-36")
def test_b45_a_one_process_call_over_the_profiles_cap_is_refused_before_its_launch(direct, capped, tmp_path, capsys):
    """The plan allowed the call (its own cap was high); the profile's cap, read before the call's launch, is under
    the call's estimate: refused (exit 2) after its cut, before its step. A second call is refused before any step,
    since the published cut is over the cap; and a profile without `cluster_cpu_min_max` is refused too."""
    root = tmp_path / "root"
    shutil.copytree(capped[0], root, symlinks=True)
    _, _, _, rebuild = forced_numbers(root)
    plan_forced(direct, root, 1e9)
    first = FAST[0]
    prof = dataclasses.replace(direct["prof"], cluster_cpu_min_max=rebuild / 2)
    calls = Calls()
    assert bs.build_stage(prof, build_of(direct, root), run=calls) == runstep.EXIT_REFUSED
    err = capsys.readouterr().err
    assert f"refused before {REL}.{first}" in err and "DEC-29" in err and "cluster_cpu_min_max" in err
    assert l2cluster.read_record(root, REL, f"{first}.cut")["one_process"]["decision"] == "DEC-29"
    assert f"stage={REL}/{first}.json" not in records(root)
    seen = [s for s, _ in calls.seen]
    assert seen[-1] == f"{first}.cut" and first not in seen
    before = records(root)
    calls = Calls()
    assert bs.build_stage(prof, build_of(direct, root), run=calls) == runstep.EXIT_REFUSED
    err = capsys.readouterr().err
    assert "refused before any step" in err and first in err and "DEC-29" in err
    assert calls.seen == [] and records(root) == before
    calls = Calls()
    none = dataclasses.replace(direct["prof"], cluster_cpu_min_max=None)
    assert bs.build_stage(none, build_of(direct, root), run=calls) == runstep.EXIT_REFUSED
    err = capsys.readouterr().err
    assert f"refused before {REL}.{first}" in err and "no cluster_cpu_min_max" in err and "DEC-29" in err
    assert f"stage={REL}/{first}.json" not in records(root)
