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
