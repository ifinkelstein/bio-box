"""W-10: `bio-box search run` takes the host profile; W-03: and the directories it writes, outside the data root
(ring R1 by directory)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from biobox import cli, runstep
from biobox.search import runner

RECIPE = "no-such-recipe.yaml"  # never read: every case here is refused or intercepted first


def _profile(tmp_path, **over):
    data = {"max_threads_day": 4, "max_threads_night": 4, "mem_scope": "200G", "tmp": "/tmp/x",
            "durable": "/d", "bulk": "/b", "nice_monolithic": 19, "nice_fanout": 15,
            "day_start_hour": 8, "day_end_hour": 20, "cpu_min_per_process": 300,
            "small_file_floor_bytes": 67108864, "require_mem_cap": True}
    data.update(over)
    p = tmp_path / "profile.yaml"
    p.write_text(yaml.safe_dump(data))
    return p


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)
    monkeypatch.delenv("BIOBOX_POD_UNIT", raising=False)


@pytest.fixture
def no_cap(tmp_path, monkeypatch):
    """No systemd-run on PATH and an empty cgroup tree: neither a scope to open nor a cap above."""
    monkeypatch.setattr(runstep.shutil, "which", lambda name: None)
    monkeypatch.setattr(runstep, "PROC_ROOT", tmp_path / "proc")
    monkeypatch.setattr(runstep, "CGROUP_ROOT", tmp_path / "cgroup")


@pytest.fixture
def captured(monkeypatch):
    seen = {}

    def fake_submit(root, recipe, release, limits, work, fault, **kw):
        seen["limits"] = limits
        seen["work"], seen["kw"] = work, kw  # card W-03: the work, pack and runs directories
        return "runid"

    monkeypatch.setattr(runner, "submit", fake_submit)
    return seen


def _argv(tmp_path, *extra):
    return ["search", "run", RECIPE, "--root", str(tmp_path), "--release", "rS0.1", *extra]


@pytest.mark.covers("W-10")
def test_required_cap_without_scope_or_cgroup_is_refused(tmp_path, no_cap, capsys):
    code = cli.main(_argv(tmp_path, "--profile", str(_profile(tmp_path))))
    assert code == 2
    assert "scripts/pod_run.sh" in capsys.readouterr().err


@pytest.mark.covers("W-10")
def test_profile_from_environment(tmp_path, no_cap, monkeypatch, capsys):
    monkeypatch.setenv("BIOBOX_PROFILE", str(_profile(tmp_path)))
    assert cli.main(_argv(tmp_path)) == 2
    assert "scripts/pod_run.sh" in capsys.readouterr().err


@pytest.mark.covers("W-10")
def test_profile_that_does_not_require_a_cap_runs(tmp_path, no_cap, captured):
    prof = _profile(tmp_path, require_mem_cap=False)
    assert cli.main(_argv(tmp_path, "--profile", str(prof))) == 0
    assert captured["limits"].require_cap is False


@pytest.mark.covers("W-10")
def test_cpus_above_profile_cap_is_refused(tmp_path, captured, capsys):
    code = cli.main(_argv(tmp_path, "--profile", str(_profile(tmp_path)), "--cpus", "5"))
    assert code == 2
    assert "4" in capsys.readouterr().err
    assert "limits" not in captured


@pytest.mark.covers("W-10")
def test_mem_gb_above_profile_cap_is_refused(tmp_path, captured, capsys):
    code = cli.main(_argv(tmp_path, "--profile", str(_profile(tmp_path, mem_scope="8G")), "--mem-gb", "9"))
    assert code == 2
    assert "mem" in capsys.readouterr().err
    assert "limits" not in captured


@pytest.mark.covers("W-10")
def test_profile_sets_limits_and_flags_may_lower(tmp_path, captured):
    prof = _profile(tmp_path, require_mem_cap=False, max_threads_day=4, max_threads_night=4)
    assert cli.main(_argv(tmp_path, "--profile", str(prof))) == 0
    lim = captured["limits"]
    assert (lim.cpus, lim.nice, lim.require_cap) == (4, 19, False)
    assert cli.main(_argv(tmp_path, "--profile", str(prof), "--cpus", "2", "--mem-gb", "3")) == 0
    lim = captured["limits"]
    assert (lim.cpus, lim.mem_gb) == (2, 3.0)


@pytest.mark.covers("W-10")
def test_profile_cap_is_the_default_memory(tmp_path, captured):
    prof = _profile(tmp_path, require_mem_cap=False, mem_scope="8G")
    assert cli.main(_argv(tmp_path, "--profile", str(prof))) == 0
    assert captured["limits"].mem_gb == 8 * 1024 ** 3 / runner.GB


@pytest.mark.covers("W-10")
def test_without_profile_limits_are_unchanged(tmp_path, captured):
    assert cli.main(_argv(tmp_path)) == 0
    assert captured["limits"] == runner.Limits(cpus=8, cpu_min=300.0, wall_s=6 * 3600.0)
    assert cli.main(_argv(tmp_path, "--cpus", "64")) == 0
    assert captured["limits"].cpus == 64


# --- W-03: where a search writes -------------------------------------------------------------------------------


@pytest.mark.covers("W-03")
def test_w03_default_directories_with_and_without_a_profile(tmp_path, monkeypatch):
    t = tmp_path / "ptmp"
    assert runner.default_runs_dir(t, user="alice") == t / "search" / "alice" / "runs"
    assert runner.default_work_base(t, user="alice") == t / "search" / "alice" / "work"
    assert runner.default_runs_dir(t, user="bob") != runner.default_runs_dir(t, user="alice")
    assert runner.default_pack_base(t) == t / "packs"  # one per host: a pack is a cache of a release
    monkeypatch.setenv("LOGNAME", "carol")
    assert runner.default_runs_dir(t) == t / "search" / "carol" / "runs"
    monkeypatch.delenv("BIOBOX_SEARCH_WORK", raising=False)
    assert runner.default_work_base() == Path("/tmp/atd/search")
    assert runner.default_runs_dir() == Path("/tmp/atd/search/runs")
    assert runner.default_pack_base() == Path("/tmp/atd/packs")
    monkeypatch.setenv("BIOBOX_SEARCH_WORK", str(tmp_path / "w"))
    assert (runner.default_work_base(), runner.default_runs_dir(), runner.default_pack_base()) == (
        tmp_path / "w", tmp_path / "w" / "runs", tmp_path / "packs")


@pytest.mark.covers("W-03")
def test_w03_search_run_takes_its_directories_from_the_profile_and_the_flags(tmp_path, captured, monkeypatch):
    monkeypatch.setenv("LOGNAME", "alice")
    t = tmp_path / "ptmp"
    prof = _profile(tmp_path, require_mem_cap=False, tmp=str(t))
    argv = ["search", "run", RECIPE, "--root", str(tmp_path / "root"), "--release", "rS0.1", "--profile", str(prof)]
    assert cli.main(argv) == 0
    assert captured["work"] == t / "search" / "alice" / "work"
    assert captured["kw"] == {"packs": t / "packs", "runs": t / "search" / "alice" / "runs", "out": None}
    flags = ["--work", str(tmp_path / "w"), "--pack-dir", str(tmp_path / "p"), "--runs", str(tmp_path / "r")]
    assert cli.main(argv + flags) == 0
    assert captured["work"] == tmp_path / "w"
    assert captured["kw"] == {"packs": tmp_path / "p", "runs": tmp_path / "r", "out": None}
    assert cli.main(argv + flags + ["--out", str(tmp_path / "one-run")]) == 0  # --out is the run's own directory
    assert captured["kw"] == {"packs": tmp_path / "p", "runs": None, "out": tmp_path / "one-run"}


# --- W-13: a search in a pod_run.sh unit -----------------------------------------------------------------------


@pytest.fixture
def waited(monkeypatch):
    """`runner.run` intercepted: what `search run --wait` asked of it."""
    seen = {}

    def fake_run(root, recipe, release, **kw):
        seen.update(kw)
        return {"run_id": b"\x13" * 16, "status": "ok"}

    monkeypatch.setattr(runner, "run", fake_run)
    return seen


@pytest.mark.covers("W-13")
def test_w13_search_run_without_wait_in_a_pod_run_unit_is_refused(tmp_path, captured, monkeypatch, capsys):
    """Card W-13 (1): in a `pod_run.sh` unit (a transient service started with `systemd-run --collect`) the unit
    ends with its main process and systemd stops every process left in it, the detached runner included (measured
    on the pod: the run ended `failed: runner exited without closing the run`). So it is refused, before anything
    is checked or recorded, with a message that says to use --wait."""
    monkeypatch.setenv("BIOBOX_POD_UNIT", "w13-x")
    assert cli.main(_argv(tmp_path)) == 2
    err = capsys.readouterr().err
    assert "--wait" in err and "w13-x" in err and "pod_run.sh" in err, err
    assert captured == {}  # submit was never called


@pytest.mark.covers("W-13")
def test_w13_search_run_with_wait_in_a_pod_run_unit_runs_and_records_its_pending_record(tmp_path, waited,
                                                                                       monkeypatch, capsys):
    monkeypatch.setenv("BIOBOX_POD_UNIT", "w13-x")
    assert cli.main(_argv(tmp_path, "--wait")) == 0
    assert waited["pending"] is True
    assert capsys.readouterr().out.strip() == "13" * 16 + "\tok"


@pytest.mark.covers("W-13")
def test_w13_the_runner_that_submit_starts_leaves_the_pending_record_to_submit(tmp_path, waited):
    """`submit` writes the pending record itself, before it returns; the runner it starts is told so (`--submitted`,
    checked on a real submit in tests/r2/test_micro_s0_search.py) and writes none of its own: two writers of one
    record would race."""
    assert cli.main(_argv(tmp_path, "--wait", "--submitted")) == 0
    assert waited["pending"] is False


# --- W-07: `bio-box search hmm`, one command; its refusals come before any run --------------------------------


def _hmm_argv(tmp_path, *extra):
    return ["search", "hmm", "--root", str(tmp_path / "root"), "--release", "rS0.1", "--hmm", str(tmp_path / "x.hmm"),
            "--bits", "30", "--out", str(tmp_path / "out"), "--pack-dir", str(tmp_path / "packs"),
            "--runs", str(tmp_path / "runs"), "--work", str(tmp_path / "work"), *extra]


@pytest.fixture
def no_run(monkeypatch):
    """`runner.run` replaced by a tripwire: these refusals must come before a run starts."""
    def tripwire(*a, **k):
        raise AssertionError("a run was started")

    monkeypatch.setattr(runner, "run", tripwire)
    monkeypatch.setattr(runner, "submit", tripwire)


def _nothing_written(tmp_path):
    return not any((tmp_path / n).exists() for n in ("out", "runs", "work", "packs"))


@pytest.mark.covers("W-07", "S-04")
def test_w07_search_hmm_without_a_positive_control_is_refused(tmp_path, no_run, capsys):
    (tmp_path / "root").mkdir()
    assert cli.main(_hmm_argv(tmp_path)) == 2
    err = capsys.readouterr().err
    assert "positive control" in err and "--positive" in err, err
    assert _nothing_written(tmp_path)


@pytest.mark.covers("W-07", "W-03")
def test_w07_search_hmm_refuses_an_output_directory_under_the_root_or_one_that_holds_files(tmp_path, no_run, capsys):
    root = tmp_path / "root"
    root.mkdir()
    argv = _hmm_argv(tmp_path, "--positive", "NC_1.1:0-90")
    assert cli.main(argv + ["--out", str(root / "out")]) == 2
    assert "data root" in capsys.readouterr().err
    assert not (root / "out").exists()
    (tmp_path / "out").mkdir()
    (tmp_path / "out" / "old.tsv").write_text("x\n")
    assert cli.main(argv) == 2
    assert "not empty" in capsys.readouterr().err
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == ["old.tsv"]


@pytest.mark.covers("W-07", "W-03")
def test_w07_search_hmm_without_the_releases_pack_names_the_pack_build_command(tmp_path, no_run, capsys):
    """A pack build is a long job: the command says which one to run and exits 2; it never builds the pack."""
    (tmp_path / "root").mkdir()
    assert cli.main(_hmm_argv(tmp_path, "--positive", "NC_1.1:0-90")) == 2
    err = capsys.readouterr().err
    assert f"bio-box pack build --root {tmp_path / 'root'} --release rS0.1 --pack-dir {tmp_path / 'packs'}" in err, err
    assert _nothing_written(tmp_path)


@pytest.mark.covers("W-07", "W-13")
def test_w07_search_hmm_always_waits_for_its_run_with_a_pending_record(tmp_path, monkeypatch, capsys):
    """Inside a `pod_run.sh` unit a search must run in the foreground (card W-13): `search hmm` runs its search in
    this process, with the pending record that lets `search status` and `search cancel` find it."""
    from biobox.search import pack

    (tmp_path / "root").mkdir()
    (tmp_path / "x.hmm").write_text("")
    monkeypatch.setenv("BIOBOX_POD_UNIT", "w07-x")
    monkeypatch.setattr(pack, "find", lambda packs, release, **k: tmp_path / "packs" / "p")
    monkeypatch.setattr(pack, "read_manifest", lambda pdir: ({"tape_version": "t-x"}, b""))
    monkeypatch.setattr(pack, "manifest_tape_version", lambda m: "t-x")
    monkeypatch.setattr(cli, "_release_files_or_refusal", lambda a, root: (object(), []))
    monkeypatch.setattr(cli, "_positive_uid", lambda *a, **k: ("ab" * 16, []))
    seen = {}

    def fake_run(root, recipe, release, **kw):
        seen.update(kw, recipe=recipe)
        return {"run_id": b"\x07" * 16, "status": "failed", "error": "positive control not recovered"}

    monkeypatch.setattr(runner, "run", fake_run)
    monkeypatch.setattr(runner, "submit", lambda *a, **k: (_ for _ in ()).throw(AssertionError("submitted")))
    assert cli.main(_hmm_argv(tmp_path, "--positive", "NC_1.1:0-90")) == 1
    assert seen["pending"] is True and seen["on_pending"] is not None
    assert seen["runs"] == tmp_path / "runs" and seen["packs"] == tmp_path / "packs"
    out = capsys.readouterr()
    assert out.out.splitlines()[0] == "07" * 16 + "\tfailed"
    assert "positive control not recovered" in out.err
    # a run that is not ok leaves the recipe it ran in DIR, and nothing else
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == ["recipe.yaml"]
    assert (tmp_path / "out" / "recipe.yaml").read_text() == seen["recipe"]
