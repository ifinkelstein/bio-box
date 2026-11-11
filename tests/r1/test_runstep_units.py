"""C-07: the memory cap a step or a search worker runs under (ring R1 by directory).

`systemd-run` is not visible inside the image, so the cap is opened on the host around
`singularity exec` (`scripts/pod_run.sh`) and run-step and the search runner, inside, find it
in their cgroup. Three cases on a fake /proc and cgroup v2 tree: systemd-run available and no
cap (a scope is opened, with MemorySwapMax=0); already capped by the unit around the process
(no systemd-run, the inherited cap is recorded); neither (refused when the profile requires a
cap, naming scripts/pod_run.sh). One test reads the real tree and runs only where a cap is
visible, as under scripts/pod_pytest.sh on the pod.
"""

from __future__ import annotations

import json
import os
import resource
import shutil
import sys
from datetime import datetime
from pathlib import Path

import pytest
import yaml

from biobox import profile, runstep
from biobox.search import runner
from biobox.search.recipe import RecipeError

OCC = "l1/protein_occurrence/gene_call_version=x/stage=S0"
NOON_MONDAY = datetime(2026, 11, 9, 12, 0)
MIB = 2**20
UNIT = "user.slice/user-1000.slice/user@1000.service/app.slice/c07-test.service"
LEAF = "/" + UNIT
UNCAPPED = {f"{UNIT}/memory.max": "max\n", "user.slice/user-1000.slice/memory.max": "max\n",
            "user.slice/memory.max": "max\n"}
CAPPED = {**UNCAPPED, f"{UNIT}/memory.max": f"{512 * MIB}\n"}
CAPPED_ABOVE = {**UNCAPPED, "user.slice/user-1000.slice/memory.max": f"{768 * MIB}\n"}  # on an ancestor

# writes one output under $BIOBOX_STAGING and a marker: argv[1] the output path, argv[2] the marker
STEP = ("import os, pathlib, sys; s = pathlib.Path(os.environ['BIOBOX_STAGING']) / sys.argv[1]; "
        "s.parent.mkdir(parents=True); s.write_bytes(b'x' * 10); pathlib.Path(sys.argv[2]).write_text('ran')")


def _tree(tmp_path: Path, files: dict[str, str], line: str = f"0::{LEAF}\n") -> tuple[Path, Path]:
    """Fake `/proc` and `/sys/fs/cgroup` roots; `files` maps cgroup-relative paths to their text."""
    proc, cg = tmp_path / "proc", tmp_path / "cg"
    (proc / "self").mkdir(parents=True)
    (proc / "self" / "cgroup").write_text(line)
    cg.mkdir()
    for rel, text in files.items():
        (cg / rel).parent.mkdir(parents=True, exist_ok=True)
        (cg / rel).write_text(text)
    return proc, cg


def _fake_systemd_run(tmp_path: Path) -> tuple[Path, Path]:
    """A stand-in `systemd-run` on PATH: logs its arguments, drops its own options, execs the rest."""
    bindir, log = tmp_path / "bin", tmp_path / "systemd-run.log"
    bindir.mkdir()
    exe = bindir / "systemd-run"
    exe.write_text(f'#!/bin/sh\nprintf "%s\\n" "$@" > "{log}"\nwhile [ "$1" != nice ]; do shift; done\nexec "$@"\n')
    exe.chmod(0o755)
    return bindir, log


def _prof(tmp_path: Path, require: bool | None, mem_scope: str = "1G") -> profile.Profile:
    data = {
        "max_threads_day": 2, "max_threads_night": 2, "mem_scope": mem_scope,
        "tmp": str(tmp_path / "tmp"), "durable": str(tmp_path / "durable"), "bulk": str(tmp_path / "bulk"),
        "nice_monolithic": 19, "nice_fanout": 15, "day_start_hour": 8, "day_end_hour": 20,
        "cpu_min_per_process": 300, "small_file_floor_bytes": 1000,
    }
    if require is not None:
        data["require_mem_cap"] = require
    p = tmp_path / "profile.yaml"
    p.write_text(yaml.safe_dump(data))
    return profile.load(p)


def _step(prof, root, proc, cg, marker, **kw) -> int:
    return runstep.run_step(prof, root, "S0", "occ",
                            [sys.executable, "-c", STEP, f"{OCC}/part-0.parquet", str(marker)],
                            clock=lambda: NOON_MONDAY, load=lambda: 0.0, proc_root=proc, cgroup_root=cg, **kw)


def _record(root: Path) -> dict:
    return json.loads((root / "runs" / "stage=S0" / "occ.json").read_text())


def _no_systemd(name: str) -> str | None:
    return None if name == "systemd-run" else shutil.which(name)


# --- detecting a cap ----------------------------------------------------------------------------


@pytest.mark.covers("C-07")
def test_c07_detect_cap_from_the_leaf_or_an_ancestor(tmp_path):
    proc, cg = _tree(tmp_path / "a", CAPPED)
    assert runstep.detect_cap(proc, cg) == {"memory_max_bytes": 512 * MIB, "limited_at": LEAF, "cgroup": LEAF}
    proc, cg = _tree(tmp_path / "b", CAPPED_ABOVE)
    got = runstep.detect_cap(proc, cg)
    assert got["memory_max_bytes"] == 768 * MIB and got["limited_at"] == "/user.slice/user-1000.slice"
    assert got["cgroup"] == LEAF


@pytest.mark.covers("C-07")
def test_c07_no_cap_detected(tmp_path):
    assert runstep.detect_cap(*_tree(tmp_path / "max", UNCAPPED)) is None  # every level says max
    assert runstep.detect_cap(*_tree(tmp_path / "v1", {}, "12:memory:/x\n")) is None
    # inside the image without the bind: /sys/fs/cgroup is an empty directory
    assert runstep.detect_cap(*_tree(tmp_path / "empty", {})) is None
    assert runstep.detect_cap(tmp_path / "no-proc", tmp_path / "no-cg") is None


# --- deciding how the command is capped ---------------------------------------------------------


def _found(name: str) -> str:
    return f"/usr/bin/{name}"


@pytest.mark.covers("C-07")
@pytest.mark.parametrize("files,which,mem,how", [
    (UNCAPPED, _found, "1G", "opened"),             # scope available, no cap: open one
    (CAPPED, _found, "1G", "inherited"),            # capped tighter than asked: keep it, open nothing
    (CAPPED, _found, "256M", "opened"),             # asked for less than the unit gives: a tighter scope
    (CAPPED, lambda n: None, "1G", "inherited"),    # inside the image: no systemd-run, the unit's cap
    (CAPPED, lambda n: None, "256M", "inherited"),
    (CAPPED_ABOVE, lambda n: None, None, "inherited"),
])
def test_c07_memory_cap_how(tmp_path, files, which, mem, how):
    proc, cg = _tree(tmp_path, files)
    got = runstep.memory_cap(mem, require=True, system="Linux", which=which, proc_root=proc, cgroup_root=cg)
    assert got["how"] == how
    if how == "opened":
        assert got["memory_max"] == mem
    else:
        assert got["memory_max_bytes"] == runstep.detect_cap(proc, cg)["memory_max_bytes"]


@pytest.mark.covers("C-07")
def test_c07_memory_cap_neither(tmp_path):
    proc, cg = _tree(tmp_path, UNCAPPED)
    kw = {"system": "Linux", "which": lambda n: None, "proc_root": proc, "cgroup_root": cg}
    assert runstep.memory_cap("1G", require=False, **kw) == {"how": "none"}
    with pytest.raises(runstep.NoCap, match="scripts/pod_run.sh"):
        runstep.memory_cap("1G", require=True, **kw)
    with pytest.raises(runstep.NoCap, match="scripts/pod_run.sh"):  # no cgroups at all (macOS)
        runstep.memory_cap("1G", require=True, system="Darwin", which=_found, proc_root=proc, cgroup_root=cg)


@pytest.mark.covers("C-07")
def test_c07_an_inherited_cap_opens_no_scope():
    cmd = ["mmseqs", "easy-search"]
    cap = {"how": "inherited", "memory_max_bytes": 512 * MIB, "limited_at": LEAF, "cgroup": LEAF}
    assert runstep.scope_argv(cmd, 19, mem_max="1G", inherited=cap, system="Linux", which=_found) == (
        ["nice", "-n", "19", *cmd], "inherited")
    assert runstep.scope_argv(cmd, 19, mem_max="1G", inherited=cap, system="Linux", which=lambda n: None) == (
        cmd, "inherited")


# --- run-step in the three cases ----------------------------------------------------------------


@pytest.mark.covers("C-07")
def test_c07_run_step_opens_a_scope_when_nothing_caps_it(tmp_path, monkeypatch):
    bindir, log = _fake_systemd_run(tmp_path)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    proc, cg = _tree(tmp_path, UNCAPPED)
    root, marker = tmp_path / "root", tmp_path / "marker"
    root.mkdir()
    which = lambda n: str(bindir / n) if n == "systemd-run" else shutil.which(n)  # noqa: E731
    assert _step(_prof(tmp_path, True), root, proc, cg, marker, system="Linux", which=which) == 0
    rec = _record(root)
    assert rec["mem_cap"] == {"how": "opened", "memory_max": "1G"}
    assert rec["argv"][:8] == ["systemd-run", "--user", "--scope", "--quiet", "-p", "MemoryMax=1G",
                               "-p", "MemorySwapMax=0"]
    assert log.read_text().splitlines()[:6] == ["--user", "--scope", "--quiet", "-p", "MemoryMax=1G", "-p"]
    assert marker.read_text() == "ran"


@pytest.mark.covers("C-07")
@pytest.mark.parametrize("files,limited_at,cap", [(CAPPED, LEAF, 512 * MIB),
                                                  (CAPPED_ABOVE, "/user.slice/user-1000.slice", 768 * MIB)])
def test_c07_run_step_in_a_capped_unit_records_the_cap_and_opens_nothing(tmp_path, monkeypatch, files,
                                                                         limited_at, cap):
    bindir, log = _fake_systemd_run(tmp_path)  # on PATH, but not found by `which`: as inside the image
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    proc, cg = _tree(tmp_path, files)
    root, marker = tmp_path / "root", tmp_path / "marker"
    root.mkdir()
    assert _step(_prof(tmp_path, True), root, proc, cg, marker, system="Linux", which=_no_systemd) == 0
    rec = _record(root)
    assert rec["mem_cap"] == {"how": "inherited", "memory_max_bytes": cap, "limited_at": limited_at, "cgroup": LEAF}
    assert "systemd-run" not in rec["argv"] and rec["argv"][:3] == ["nice", "-n", "19"]
    assert not log.exists()
    assert marker.read_text() == "ran"


@pytest.mark.covers("C-07")
def test_c07_run_step_refuses_without_a_scope_or_a_cap(tmp_path, capsys):
    proc, cg = _tree(tmp_path, UNCAPPED)
    root, marker = tmp_path / "root", tmp_path / "marker"
    root.mkdir()
    prof = _prof(tmp_path, True)
    assert _step(prof, root, proc, cg, marker, system="Linux", which=_no_systemd) == runstep.EXIT_REFUSED
    err = capsys.readouterr().err
    assert "scripts/pod_run.sh" in err and "require_mem_cap" in err
    assert not marker.exists()
    assert not (Path(prof.tmp) / "staging" / "S0.occ").exists()
    assert sorted(root.rglob("*")) == []


@pytest.mark.covers("C-07")
@pytest.mark.parametrize("require", [False, None])
def test_c07_run_step_without_a_cap_runs_when_the_profile_does_not_require_one(tmp_path, require):
    proc, cg = _tree(tmp_path, UNCAPPED)
    root, marker = tmp_path / "root", tmp_path / "marker"
    root.mkdir()
    assert _step(_prof(tmp_path, require), root, proc, cg, marker, system="Linux", which=_no_systemd) == 0
    assert _record(root)["mem_cap"] == {"how": "none"}
    assert marker.read_text() == "ran"


# --- the search runner --------------------------------------------------------------------------


CAP = {"how": "inherited", "memory_max_bytes": 8 * 10**9, "limited_at": LEAF, "cgroup": LEAF}


@pytest.mark.covers("C-07")
def test_c07_worker_in_a_capped_unit_runs_under_nice_with_no_scope():
    job = Path("/w/run/job.json")
    lim = runner.Limits(cpus=2, mem_gb=18, nice=15)
    assert runner.worker_argv(job, lim, inherited=CAP, system="Linux", which=_found) == (
        ["nice", "-n", "15", *runner.worker_command(job)], "inherited")


@pytest.mark.covers("C-07")
def test_c07_worker_limits_under_an_inherited_cap():
    cpu = {resource.RLIMIT_CPU: (120, 125)}
    cap = CAP["memory_max_bytes"]
    # the unit's cap is the memory limit: no RLIMIT_AS, and the oom fault allocates past the cap
    for lim in (runner.Limits(cpu_min=2), runner.Limits(cpu_min=2, mem_gb=8), runner.Limits(cpu_min=2, mem_gb=20)):
        assert runner.rlimits(lim, "inherited", "Linux", cap_bytes=cap) == cpu
        env = runner.worker_env(lim, "inherited", "Linux", {}, cap_bytes=cap)
        assert env["BIOBOX_MEMORY_MAX"] == str(cap) and "BIOBOX_RLIMIT_AS" not in env
    # asked for less than the unit gives, with no scope to open: RLIMIT_AS at the request as well
    lim = runner.Limits(cpu_min=2, mem_gb=4)
    assert runner.rlimits(lim, "inherited", "Linux", cap_bytes=cap) == {**cpu, resource.RLIMIT_AS: (4 * 10**9,) * 2}
    env = runner.worker_env(lim, "inherited", "Linux", {}, cap_bytes=cap)
    assert env["BIOBOX_RLIMIT_AS"] == str(4 * 10**9) and "BIOBOX_MEMORY_MAX" not in env


@pytest.mark.covers("C-07")
def test_c07_runner_memory_cap_three_cases(tmp_path):
    proc, cg = _tree(tmp_path / "u", UNCAPPED)
    lim = runner.Limits(mem_gb=4, require_cap=True)
    assert runner.memory_cap(lim, system="Linux", which=_found, proc_root=proc, cgroup_root=cg) == {
        "how": "opened", "memory_max": str(4 * 10**9)}
    with pytest.raises(RecipeError, match="scripts/pod_run.sh"):
        runner.memory_cap(lim, system="Linux", which=lambda n: None, proc_root=proc, cgroup_root=cg)
    assert runner.memory_cap(runner.Limits(mem_gb=4), system="Linux", which=lambda n: None, proc_root=proc,
                             cgroup_root=cg) == {"how": "none"}
    proc, cg = _tree(tmp_path / "c", CAPPED)
    assert runner.memory_cap(lim, system="Linux", which=lambda n: None, proc_root=proc, cgroup_root=cg)["how"] == \
        "inherited"


@pytest.mark.covers("C-07")
def test_c07_runner_refuses_before_it_records_anything(tmp_path, monkeypatch):
    proc, cg = _tree(tmp_path, UNCAPPED)
    monkeypatch.setattr(runstep, "PROC_ROOT", proc)
    monkeypatch.setattr(runstep, "CGROUP_ROOT", cg)
    empty = tmp_path / "empty-bin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))  # no systemd-run
    root, work = tmp_path / "root", tmp_path / "work"
    with pytest.raises(RecipeError, match="scripts/pod_run.sh"):
        runner.run(root, "not even a recipe", "rX", limits=runner.Limits(require_cap=True), work_base=work)
    with pytest.raises(RecipeError, match="scripts/pod_run.sh"):
        runner.submit(root, tmp_path / "recipe.yaml", "rX", runner.Limits(require_cap=True), work)
    assert not root.exists() and not work.exists()


@pytest.mark.covers("C-07")
def test_c07_oom_evidence_in_a_shared_unit_counts_only_new_kills(tmp_path):
    work, cg = tmp_path / "work", tmp_path / "cg"
    work.mkdir()
    (work / runner.CGROUP_FILE).write_text(LEAF + "\n")
    (cg / UNIT).mkdir(parents=True)
    (cg / UNIT / "memory.events").write_text("low 0\nhigh 0\nmax 3\noom 1\noom_kill 1\n")
    assert runner.scope_oom_killed(work, cgroup_root=cg, systemctl=None) is True
    assert runner.scope_oom_killed(work, cgroup_root=cg, systemctl=None, baseline=1) is False
    assert runner.scope_oom_killed(work, cgroup_root=cg, systemctl=None, baseline=0) is True


# --- the real cgroup tree -----------------------------------------------------------------------


def _real_cap_by_hand() -> int | None:
    """The smallest numeric memory.max from this process's cgroup up to the root, walked here
    independently of `doctor.cgroup_memory`."""
    lines = [ln for ln in Path("/proc/self/cgroup").read_text().splitlines() if ln.startswith("0::")]
    top = Path("/sys/fs/cgroup")
    d = top / lines[0][3:].strip().lstrip("/")
    values = []
    while True:
        f = d / "memory.max"
        if f.is_file() and f.read_text().strip() != "max":
            values.append(int(f.read_text().strip()))
        if d == top:
            break
        d = d.parent
    return min(values) if values else None


@pytest.mark.covers("C-07")
def test_c07_detects_the_real_cap_it_runs_under(tmp_path):
    cap = runstep.detect_cap()
    if cap is None:
        pytest.skip("no cgroup v2 memory cap visible to this process: not inside a capped unit, or "
                    "/sys/fs/cgroup is not bound into the image (scripts/pod_pytest.sh binds it)")
    assert cap["memory_max_bytes"] == _real_cap_by_hand()
    root, marker = tmp_path / "root", tmp_path / "marker"
    root.mkdir()
    prof = _prof(tmp_path, True, mem_scope=str(cap["memory_max_bytes"]))
    code = runstep.run_step(prof, root, "S0", "occ",
                            [sys.executable, "-c", STEP, f"{OCC}/part-0.parquet", str(marker)],
                            clock=lambda: NOON_MONDAY, load=lambda: 0.0)
    assert code == 0
    rec = _record(root)
    assert rec["mem_cap"]["how"] == "inherited"
    assert rec["mem_cap"]["memory_max_bytes"] == cap["memory_max_bytes"]
    assert "systemd-run" not in rec["argv"]
