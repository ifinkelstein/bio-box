"""A-10a: evidence of who killed the measuring helper, kept when a peak is a sampled one (ring R1 by directory).

Card A-10 found that a run's peak is the helper's last sample (`*_sampled`) only when the helper itself was killed
before its final read, and that what killed it (the kernel's OOM killer, systemd stopping the scope, systemd-oomd)
cannot be told from the record. The evidence that would tell exists for a moment after the run and was thrown away:
the worker's stderr, the scope's `memory.events`, the scope unit's `Result` and `OOMPolicy`, whether systemd-oomd is
active. `runstep.peak_and_evidence` keeps it when, and only when, the source is a sampled one.

On a fake cgroup v2 tree whose `memory.events` and `memory.peak` the test writes, with a stand-in `systemctl` that
logs its calls: the helper's file as the helper writes it (sampled or final), and run-step end to end with a
stand-in `systemd-run`, a step killed with its helper (sampled: the failed attempt's record carries the evidence)
or alone (final: it does not). Nothing here asserts what the evidence means; it only has to be kept, as read.
"""

from __future__ import annotations

import json
import os
import shutil
import signal
import sys
import threading
import time
from datetime import datetime
from pathlib import Path

import pytest
import yaml

from biobox import profile, runstep
from biobox.search import runner

SCOPE = "user.slice/user-1000.slice/user@1000.service/app.slice/run-a10a.scope"
EVENTS = "low 0\nhigh 0\nmax 17\noom 1\noom_kill 1\noom_group_kill 0\n"
SHOW = "LoadState=loaded\nActiveState=failed\nResult=oom-kill\nOOMPolicy=stop\nMemoryPeak=[not set]\n"
# a stand-in systemctl: logs its arguments, answers `--user show` with SHOW and `is-active` with "inactive" (exit 3)
SYSTEMCTL = """#!/bin/sh
echo "$@" >> "$(dirname "$0")/systemctl.log"
case "$1" in
  --user) printf '%s' '{show}' ;;
  is-active) echo inactive; exit 3 ;;
esac
"""


def _systemctl(tmp_path: Path) -> Path:
    bindir = tmp_path / "bin"
    bindir.mkdir(exist_ok=True)
    exe = bindir / "systemctl"
    exe.write_text(SYSTEMCTL.replace("{show}", SHOW))
    exe.chmod(0o755)
    return exe


def _calls(exe: Path) -> list[str]:
    log = exe.parent / "systemctl.log"
    return log.read_text().splitlines() if log.exists() else []


def _scope(tmp_path: Path, *, events: str | None = EVENTS, peak: int = 10**9, memory_max: str = "1000000000") -> Path:
    cg = tmp_path / "cg"
    (cg / SCOPE).mkdir(parents=True)
    (cg / SCOPE / "memory.max").write_text(memory_max + "\n")
    (cg / SCOPE / "memory.peak").write_text(f"{peak}\n")
    if events is not None:
        (cg / SCOPE / "memory.events").write_text(events)
    return cg


def _helper_file(path: Path, **rec) -> Path:
    path.write_text(json.dumps(rec))
    return path


@pytest.mark.covers("A-10a", "C-16")
def test_a10a_a_sampled_scope_peak_keeps_the_scopes_evidence(tmp_path):
    cg = _scope(tmp_path)
    exe = _systemctl(tmp_path)
    err = tmp_path / "stderr.txt"
    err.write_text("Exception in thread Thread-1 (sample):\nOSError: [Errno 12] Cannot allocate memory\n")
    out = _helper_file(tmp_path / "peak.json", maxrss_bytes=None, memory_peak_bytes=625729536, final=False,
                       cgroup="/" + SCOPE)
    rss, source, ev = runstep.peak_and_evidence(out, cgroup_root=cg, which=lambda n: str(exe), stderr=err)
    assert (rss, source) == (625729536, "scope_memory_peak_sampled")
    assert not out.exists()  # read once, as step_peak reads it
    assert ev["cgroup"] == "/" + SCOPE and ev["cgroup_present"] is True
    assert ev["memory_events"] == {"low": 0, "high": 0, "max": 17, "oom": 1, "oom_kill": 1, "oom_group_kill": 0}
    assert (ev["memory_peak"], ev["memory_max"]) == (10**9, "1000000000")
    assert ev["unit"] == "run-a10a.scope"
    assert ev["unit_show"] == {"LoadState": "loaded", "ActiveState": "failed", "Result": "oom-kill",
                               "OOMPolicy": "stop", "MemoryPeak": "[not set]"}
    assert ev["systemd_oomd"] == "inactive"
    assert ev["stderr_tail"] == err.read_text()
    assert _calls(exe) == ["--user show -p LoadState -p ActiveState -p Result -p OOMPolicy -p MemoryPeak "
                           "run-a10a.scope", "is-active systemd-oomd"]


@pytest.mark.covers("A-10a", "C-16")
@pytest.mark.parametrize("rec,source", [
    ({"maxrss_bytes": 11 << 20, "memory_peak_bytes": 10**9, "final": True, "cgroup": "/" + SCOPE},
     "scope_memory_peak"),
    ({"maxrss_bytes": 11 << 20, "memory_peak_bytes": None, "final": True, "cgroup": None}, "child_maxrss"),
    # a helper file written before this card has no cgroup; it is read as before
    ({"maxrss_bytes": 11 << 20, "memory_peak_bytes": 10**9, "final": True}, "scope_memory_peak"),
])
def test_a10a_a_final_peak_keeps_no_evidence_and_asks_nothing(tmp_path, rec, source):
    cg = _scope(tmp_path)
    exe = _systemctl(tmp_path)
    err = tmp_path / "stderr.txt"
    err.write_text("anything\n")
    rss, got, ev = runstep.peak_and_evidence(_helper_file(tmp_path / "peak.json", **rec), cgroup_root=cg,
                                             which=lambda n: str(exe), stderr=err)
    assert got == source and rss is not None
    assert ev is None
    assert _calls(exe) == []  # a run that ended as usual costs nothing more


@pytest.mark.covers("A-10a", "C-16")
def test_a10a_an_unmeasured_peak_keeps_no_evidence(tmp_path):
    exe = _systemctl(tmp_path)
    assert runstep.peak_and_evidence(tmp_path / "absent.json", which=lambda n: str(exe)) == (None, "unmeasured",
                                                                                             None)
    assert _calls(exe) == []


@pytest.mark.covers("A-10a", "C-16")
def test_a10a_a_scope_already_removed_is_recorded_as_gone(tmp_path):
    cg = tmp_path / "cg"
    cg.mkdir()
    exe = _systemctl(tmp_path)
    out = _helper_file(tmp_path / "peak.json", maxrss_bytes=None, memory_peak_bytes=5 << 20, final=False,
                       cgroup="/" + SCOPE)
    rss, source, ev = runstep.peak_and_evidence(out, cgroup_root=cg, which=lambda n: str(exe))
    assert source == "scope_memory_peak_sampled"
    assert ev["cgroup_present"] is False
    assert (ev["memory_events"], ev["memory_peak"], ev["memory_max"]) == (None, None, None)
    assert ev["unit_show"]["Result"] == "oom-kill"  # the unit outlives its cgroup until reset-failed
    assert "stderr_tail" not in ev  # run-step's step writes its stderr to run-step's own


@pytest.mark.covers("A-10a", "C-16")
def test_a10a_a_sampled_process_peak_keeps_what_there_is_without_a_scope(tmp_path):
    exe = _systemctl(tmp_path)
    err = tmp_path / "stderr.txt"
    err.write_text("x" * 10_000 + "the end\n")
    out = _helper_file(tmp_path / "peak.json", maxrss_bytes=300 << 20, memory_peak_bytes=None, final=False)
    rss, source, ev = runstep.peak_and_evidence(out, which=lambda n: str(exe) if n == "systemctl" else None,
                                                stderr=err)
    assert (rss, source) == (300 << 20, "child_hwm_sampled")
    assert ev["cgroup"] is None and ev["cgroup_present"] is None and ev["unit"] is None and ev["unit_show"] is None
    assert ev["systemd_oomd"] == "inactive"
    assert ev["stderr_tail"] == err.read_text()[-runstep.EVIDENCE_STDERR_CHARS:]
    assert ev["stderr_tail"].endswith("the end\n") and len(ev["stderr_tail"]) == runstep.EVIDENCE_STDERR_CHARS
    assert _calls(exe) == ["is-active systemd-oomd"]


@pytest.mark.covers("A-10a", "C-16")
def test_a10a_without_systemctl_the_unit_and_oomd_are_unknown(tmp_path):
    cg = _scope(tmp_path)
    out = _helper_file(tmp_path / "peak.json", maxrss_bytes=None, memory_peak_bytes=5 << 20, final=False,
                       cgroup="/" + SCOPE)
    _, _, ev = runstep.peak_and_evidence(out, cgroup_root=cg, which=lambda n: None)
    assert ev["unit"] == "run-a10a.scope" and ev["unit_show"] is None and ev["systemd_oomd"] is None
    assert ev["memory_events"]["oom_kill"] == 1


@pytest.mark.covers("A-10a", "C-16")
def test_a10a_the_search_record_holds_the_evidence_only_when_there_is_some():
    lim = runner.Limits(cpus=1, mem_gb=1)
    cap = {"how": "opened"}
    plain = json.loads(runner._limits_json(lim, cap, "scope_memory_peak", ["bio-box", "search"], None))
    assert "peak_evidence" not in plain
    ev = {"cgroup": "/" + SCOPE, "memory_events": {"oom_kill": 1}}
    kept = json.loads(runner._limits_json(lim, cap, "scope_memory_peak_sampled", ["bio-box", "search"], ev))
    assert kept["peak_evidence"] == ev
    assert {k: v for k, v in kept.items() if k != "peak_evidence"} == {**plain,
                                                                        "peak_rss_source": "scope_memory_peak_sampled"}


# --- run-step end to end, on the fake tree, with stand-ins for systemd-run and systemctl -------------------------------

NOON_MONDAY = datetime(2026, 11, 9, 12, 0)
# plays the kernel's part: writes the scope's memory.peak (argv[1]) and memory.events (argv[2]) from argv[3], argv[4];
# then says it is ready in the marker named by its last argument (pid and process group) and waits to be killed
STEP = r"""
import os, signal, sys
open(sys.argv[1], "w").write(sys.argv[3] + "\n")
open(sys.argv[2], "w").write(sys.argv[4].replace("|", "\n"))
with open(sys.argv[-1] + ".tmp", "w") as f:
    f.write(f"{os.getpid()} {os.getpgid(0)}")
os.replace(sys.argv[-1] + ".tmp", sys.argv[-1])
while True:
    signal.pause()
"""


@pytest.fixture
def prof(tmp_path) -> profile.Profile:
    data = {"max_threads_day": 4, "max_threads_night": 6, "mem_scope": "1G", "tmp": str(tmp_path / "tmp"),
            "durable": str(tmp_path / "durable"), "bulk": str(tmp_path / "bulk"), "nice_monolithic": 19,
            "nice_fanout": 15, "day_start_hour": 8, "day_end_hour": 20, "cpu_min_per_process": 300,
            "small_file_floor_bytes": 1000}
    p = tmp_path / "profile.yaml"
    p.write_text(yaml.safe_dump(data))
    return profile.load(p)


def _wait_for(pred, timeout: float) -> bool:
    """Poll `pred` every 0.05 s until it holds or `timeout` seconds pass (AGENTS §5.3)."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if pred():
            return True
        time.sleep(0.05)
    return False


def _read(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return {}


@pytest.mark.covers("A-10a", "C-16")
@pytest.mark.skipif(sys.platform != "linux", reason="the helper samples only where /proc exists (Linux)")
@pytest.mark.parametrize("kill", ["group", "step"])
def test_a10a_run_step_keeps_the_evidence_of_a_step_killed_with_its_helper(prof, tmp_path, monkeypatch, kill):
    """`group`: the step and its helper are killed together (as in the Actions failure): the record's peak is the
    helper's last sample and the failed attempt's record carries the evidence, read from the fake scope. `step`: the
    step alone is killed, the helper reads the final peak, and the record carries no evidence."""
    root = tmp_path / "root"
    root.mkdir()
    # memory.max "max": no cap above run-step, so it opens a scope (the stand-in) and the helper reads memory.peak
    proc, cg = tmp_path / "proc", _scope(tmp_path, events="oom_kill 0\n", peak=1, memory_max="max")
    (proc / "self").mkdir(parents=True)
    (proc / "self" / "cgroup").write_text(f"0::/{SCOPE}\n")
    bindir = tmp_path / "bin"
    bindir.mkdir(exist_ok=True)
    (bindir / "systemd-run").write_text('#!/bin/sh\nwhile [ "$1" != nice ]; do shift; done\nexec "$@"\n')
    (bindir / "systemd-run").chmod(0o755)
    monkeypatch.setenv("PATH", f"{bindir}:{os.environ.get('PATH', '')}")  # the scope's argv names it by name
    exe = _systemctl(tmp_path)
    which = lambda n: str(bindir / n) if n in ("systemd-run", "systemctl") else shutil.which(n)  # noqa: E731
    peak = 7 * 2**30 + 99
    marker = Path(prof.tmp) / "a10a.marker"
    cmd = [sys.executable, "-c", STEP, str(cg / SCOPE / "memory.peak"), str(cg / SCOPE / "memory.events"),
           str(peak), EVENTS.replace("\n", "|"), str(marker)]
    codes: list[int] = []
    t = threading.Thread(target=lambda: codes.append(runstep.run_step(
        prof, root, "S0", "a10a", cmd, threads=1, clock=lambda: NOON_MONDAY, load=lambda: 0.0, system="Linux",
        which=which, proc_root=proc, cgroup_root=cg)))
    t.start()
    pgid = None
    try:
        assert _wait_for(marker.exists, 60), "the step never said it was ready"
        pid, pgid = map(int, marker.read_text().split())
        if kill == "group":
            helper_file = Path(prof.tmp) / "staging" / "S0.a10a.peak.json"
            assert _wait_for(lambda: _read(helper_file).get("memory_peak_bytes") == peak, 30)
            os.killpg(pgid, signal.SIGKILL)
        else:
            os.kill(pid, signal.SIGKILL)
        t.join(60)
    finally:
        if pgid is not None:
            try:
                os.killpg(pgid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass
        t.join(60)
    assert not t.is_alive() and codes == [1]
    att = json.loads((Path(prof.tmp) / "staging" / "S0.a10a.failed.json").read_text())
    assert att["peak_rss_bytes"] == peak
    if kill == "group":
        assert att["peak_rss_source"] == "scope_memory_peak_sampled"
        ev = att["peak_evidence"]
        assert ev["cgroup"] == "/" + SCOPE and ev["cgroup_present"] is True
        assert ev["memory_events"]["oom_kill"] == 1 and ev["memory_peak"] == peak
        assert ev["unit_show"]["Result"] == "oom-kill" and ev["systemd_oomd"] == "inactive"
    else:
        assert att["peak_rss_source"] == "scope_memory_peak"
        assert "peak_evidence" not in att
        assert _calls(exe) == []
    assert not (root / "runs").exists()  # nothing published either way
