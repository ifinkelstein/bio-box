"""C-28: the night thread allowance asks the rule's question (ring R1 by directory).

AGENTS §7.2 allows 56 threads at night "only when `ps` or `top` shows no other compute". Before this card the code
read the load average instead (above 8: the day cap), and our own scan held it near 42, so a build never had more
than 40 at night: on 2026-11-12 at 23:24 CST the account used 41.1 CPUs and every other user together under 0.5 on a
112-core host (*measured*, `ps`, card C-28).

Every reading here is injected: a fake `/proc` tree whose CPU ticks a fake wait advances, or a `CpuUse` handed in.
No test reads the real host. The cases: other users busy, other users idle, system accounts busy (not compute), the
account's own processes elsewhere (subtracted) and in this unit (not subtracted), a failed reading (the day cap), and
a launch at night that could still run when the day window opens (the day cap).
"""

from __future__ import annotations

import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pytest
import yaml

from biobox import profile, runstep

pytestmark = pytest.mark.covers("C-28")

HZ = 100
ME = 601287  # the account's uid
OTHER = 600795  # another login user
SPLUNK = 1002  # a service account above UID_MIN: still another user's process for the rule (its CPU is tiny)
MYPID = 4000
UNIT = "/user.slice/user-601287.slice/user@601287.service/app.slice/r50-P3-12.service"
TESTS = "/user.slice/user-601287.slice/user@601287.service/app.slice/run-u77.scope"
SYSTEM = "/system.slice/ampdaemon.service"
MONDAY_NIGHT = datetime(2026, 11, 9, 23, 0)  # a Monday
MONDAY_NOON = datetime(2026, 11, 9, 12, 0)
SUNDAY_NIGHT = datetime(2026, 11, 8, 23, 0)
DRIVER = Path(__file__).resolve().parents[2] / "scripts" / "stress" / "build_stage_pod" / "r50_driver.sh"


class FakeProc:
    """A fake `/proc`: one directory per process with `stat` (ppid, utime + stime, starttime), `status` (Uid) and
    `cgroup` files, as Linux writes them."""

    def __init__(self, root: Path):
        self.root = root
        root.mkdir(parents=True, exist_ok=True)
        (root / "loadavg").write_text("41.98 42.32 42.39 43/2603 2885733\n")
        self.procs: dict[int, dict] = {}

    def add(self, pid: int, uid: int, ticks: int = 0, *, ppid: int = 1, cgroup: str = SYSTEM, start: int = 1000,
            comm: str = "x y") -> None:
        self.procs[pid] = {"uid": uid, "ticks": ticks, "ppid": ppid, "cgroup": cgroup, "start": start, "comm": comm}
        self._write(pid)

    def run(self, pid: int, ticks: int) -> None:
        self.procs[pid]["ticks"] += ticks
        self._write(pid)

    def gone(self, pid: int) -> None:
        del self.procs[pid]
        for f in (self.root / str(pid)).iterdir():
            f.unlink()
        (self.root / str(pid)).rmdir()

    def _write(self, pid: int) -> None:
        p = self.procs[pid]
        d = self.root / str(pid)
        d.mkdir(exist_ok=True)
        utime, stime = p["ticks"] - p["ticks"] // 3, p["ticks"] // 3
        fields = ["S", p["ppid"], pid, pid, 0, -1, 4194560, 10, 0, 0, 0, utime, stime, 999, 999, 20, 0, 1, 0,
                  p["start"], 1000, 100]
        (d / "stat").write_text(f"{pid} ({p['comm']}) " + " ".join(str(x) for x in fields) + "\n")
        (d / "status").write_text(f"Name:\tx\nPid:\t{pid}\nPPid:\t{p['ppid']}\nUid:\t{p['uid']}\t{p['uid']}\t"
                                  f"{p['uid']}\t{p['uid']}\nGid:\t1\t1\t1\t1\n")
        (d / "cgroup").write_text(f"0::{p['cgroup']}\n")


def host(tmp_path: Path) -> FakeProc:
    """The pod as `ps` saw it: init and a scanner of root's, a system account's daemon, our build-stage (MYPID) with
    two shards in its unit, another login user's idle editor."""
    fp = FakeProc(tmp_path / "proc")
    fp.add(1, 0, 5000)
    fp.add(900, 0, 100, cgroup=SYSTEM)  # a virus scanner of root's
    fp.add(901, 999, 100, cgroup="/system.slice/ampscansvc.service")  # a system account
    fp.add(902, SPLUNK, 100, cgroup="/system.slice/splunk.service")
    fp.add(MYPID, ME, 100, ppid=3999, cgroup=UNIT)
    fp.add(MYPID + 1, ME, 100, ppid=MYPID, cgroup=UNIT)  # a shard's run-step, our child
    fp.add(MYPID + 2, ME, 100, ppid=MYPID + 1, cgroup=UNIT + "/sub")  # its step, below our cgroup
    fp.add(MYPID + 3, ME, 100, ppid=3999, cgroup=UNIT)  # a sibling in our unit (not our descendant)
    fp.add(7000, OTHER, 100, cgroup="/user.slice/user-600795.slice/session-9.scope")
    return fp


def measure(fp: FakeProc, during, window: float = 5.0) -> runstep.CpuUse:
    times = iter([100.0, 100.0 + window])
    return runstep.measure_cpu(fp.root, window_s=window, uid=ME, pid=MYPID, hz=HZ, wait=lambda s: during(),
                               now=lambda: next(times))


def test_cpu_of_other_users_is_measured_over_the_window(tmp_path):
    fp = host(tmp_path)
    use = measure(fp, lambda: fp.run(7000, 5 * HZ))  # one core for 5 s
    assert use.error is None
    assert use.others == pytest.approx(1.0)
    assert use.window_s == pytest.approx(5.0)
    idle = host(tmp_path / "idle")
    assert measure(idle, lambda: idle.run(7000, 10)).others == pytest.approx(0.02)


def test_root_and_system_accounts_are_never_other_compute(tmp_path):
    """A virus scanner (root) at 0.75 of a core and a system account's daemon (uid below 1000) at 0.2: system, not
    compute; a service account above UID_MIN (splunk, 1002) is another user's process and counts, though its own
    CPU is about 0.007 cores on the pod (ps, 2026-11-13)."""
    fp = host(tmp_path)

    def during():
        fp.run(900, 375)
        fp.run(901, 100)

    use = measure(fp, during)
    assert use.system == pytest.approx(0.95)
    assert use.others == 0.0
    assert runstep.choose_threads(PROF, MONDAY_NIGHT, use) == PROF.max_threads_night
    splunk = host(tmp_path / "splunk")
    assert measure(splunk, lambda: splunk.run(902, 5)).others == pytest.approx(0.01)


def test_the_accounts_processes_are_split_into_this_unit_and_elsewhere(tmp_path):
    """Ours in this process's cgroup or below it, or descended from it: the build's own (`own_here`). Ours elsewhere
    (a test suite in its own scope, another agent's build): `own_other`."""
    fp = host(tmp_path)
    fp.add(8000, ME, 0, ppid=3998, cgroup=TESTS)
    fp.add(8001, ME, 0, ppid=MYPID, cgroup=TESTS)  # a child of ours in another cgroup (a scope it opened)

    def during():
        for pid in (MYPID, MYPID + 1, MYPID + 2, MYPID + 3, 8001):
            fp.run(pid, 4 * 5 * HZ)  # four cores each
        fp.run(8000, 3 * 5 * HZ)  # three cores

    use = measure(fp, during)
    assert use.own_here == pytest.approx(20.0)
    assert use.own_other == pytest.approx(3.0)
    assert use.others == 0.0


def test_a_process_born_in_the_window_counts_its_cpu_and_one_that_ended_counts_none(tmp_path):
    fp = host(tmp_path)
    fp.add(7001, OTHER, 0, start=5000)

    def during():
        fp.gone(7000)  # ended: nothing of it is counted
        fp.add(7002, OTHER, 2 * HZ, start=9000)  # born in the window: all its CPU is in it
        fp.gone(7001)
        fp.add(7001, OTHER, HZ, start=9100)  # the pid taken again by a new process: its own CPU only

    use = measure(fp, during)
    assert use.others == pytest.approx((2 * HZ + HZ) / HZ / 5.0)


def test_a_reading_that_cannot_see_other_users_fails(tmp_path):
    """With `hidepid` (or in a namespace) only our own processes are visible: no other user would ever look busy, so
    the reading fails and the day cap applies. So it does without /proc."""
    fp = FakeProc(tmp_path / "proc")
    fp.add(MYPID, ME, 100, cgroup=UNIT)
    use = measure(fp, lambda: fp.run(MYPID, 100))
    assert use.error and "other user" in use.error
    assert runstep.choose_threads(PROF, MONDAY_NIGHT, use) == PROF.max_threads_day
    gone = runstep.measure_cpu(tmp_path / "no-proc", uid=ME, pid=MYPID, hz=HZ, wait=lambda s: None)
    assert gone.error
    assert runstep.choose_threads(PROF, MONDAY_NIGHT, gone) == PROF.max_threads_day


def test_the_window_waits_once_for_its_length(tmp_path):
    fp = host(tmp_path)
    waits: list[float] = []
    times = iter([0.0, 5.0])
    runstep.measure_cpu(fp.root, window_s=5.0, uid=ME, pid=MYPID, hz=HZ, wait=waits.append, now=lambda: next(times))
    assert waits == [5.0]
    assert runstep.CPU_WINDOW_S >= 1.0


# --- the allowance ------------------------------------------------------------------------------------------------

PROF = profile.Profile(max_threads_day=40, max_threads_night=56, mem_scope="200G", tmp="/tmp/x", durable="/d",
                       bulk="/b", nice_monolithic=19, nice_fanout=15, day_start_hour=8, day_end_hour=20,
                       cpu_min_per_process=300, small_file_floor_bytes=1000)
IDLE = runstep.CpuUse(others=0.01, system=0.3, own_other=0.0, own_here=41.1, window_s=5.0)
BUSY = runstep.CpuUse(others=runstep.OTHERS_BUSY_CPUS + 0.5, system=0.3, own_here=41.1, window_s=5.0)
FAILED = runstep.CpuUse(error="no /proc")


@pytest.mark.parametrize("now,use,cap", [
    (MONDAY_NIGHT, IDLE, 56),  # our own 41.1 CPUs do not make the night busy
    (MONDAY_NIGHT, BUSY, 40),  # another user computes: the day cap at night
    (SUNDAY_NIGHT, BUSY, 40),
    (MONDAY_NOON, IDLE, 40),  # the day cap is unchanged
    (MONDAY_NOON, BUSY, 40),
    (MONDAY_NIGHT, FAILED, 40),  # no reading: the day cap
    (MONDAY_NIGHT, runstep.CpuUse(others=runstep.OTHERS_BUSY_CPUS), 56),  # at the threshold: not above it
])
def test_the_night_cap_holds_only_while_no_other_user_computes(now, use, cap):
    assert runstep.choose_threads(PROF, now, use) == cap
    assert runstep.allowance(PROF, now, use).cap == cap


def test_the_allowance_is_the_cap_less_our_other_processes_floored_at_zero_and_at_most_the_ask():
    """A test suite of ours using 6.4 CPUs elsewhere leaves 56 - 6 at night; the build's own 41.1 CPUs (in its unit)
    are not subtracted; never more than asked; never below zero."""
    tests = runstep.CpuUse(others=0.01, own_other=6.4, own_here=41.1, window_s=5.0)
    assert runstep.allowance(PROF, MONDAY_NIGHT, tests).threads == 50
    assert runstep.allowance(PROF, MONDAY_NIGHT, tests, want=16).threads == 16
    assert runstep.allowance(PROF, MONDAY_NOON, tests).threads == 34
    full = runstep.CpuUse(own_other=47.0)
    assert runstep.allowance(PROF, MONDAY_NOON, full).threads == 0
    assert runstep.allowance(PROF, MONDAY_NIGHT, IDLE).threads == 56
    assert runstep.allowance(PROF, MONDAY_NIGHT, runstep.CpuUse(own_other=0.3)).threads == 56  # an idle shell


def test_a_night_launch_that_could_run_into_the_day_window_gets_the_day_cap():
    """Card C-28 (3), beyond card C-12's re-read before each launch: a launch at night of T threads is counted against
    the day cap when its CPU budget (`cpu_min_per_process` / T minutes at full use) reaches the day window. 4 threads:
    75 minutes, so from 06:45 on a weekday."""
    assert runstep.allowance(PROF, datetime(2026, 11, 9, 6, 44), IDLE, threads=4).cap == 56
    assert runstep.allowance(PROF, datetime(2026, 11, 9, 6, 46), IDLE, threads=4).cap == 40
    assert runstep.allowance(PROF, datetime(2026, 11, 9, 7, 30), IDLE, threads=56).cap == 56  # 5.4 minutes
    assert runstep.allowance(PROF, datetime(2026, 11, 9, 7, 56), IDLE, threads=56).cap == 40
    assert runstep.allowance(PROF, datetime(2026, 11, 13, 19, 0), IDLE, threads=4).cap == 40  # Friday, day
    assert runstep.allowance(PROF, datetime(2026, 11, 14, 6, 46), IDLE, threads=4).cap == 56  # Saturday: no window
    assert runstep.allowance(PROF, SUNDAY_NIGHT, IDLE, threads=1).cap == 56  # 300 minutes: Monday 04:00
    assert runstep.allowance(PROF, datetime(2026, 11, 9, 3, 1), IDLE, threads=1).cap == 40  # past 08:00


def test_the_allowance_line_names_the_readings_and_the_rule():
    a = runstep.allowance(PROF, MONDAY_NIGHT, BUSY, want=56, threads=4, load1=41.98)
    assert a.threads == 40
    line = a.line
    assert "\n" not in line
    for part in ("allowance 40", "cap 40", "other users", f"{BUSY.others:.2f}", "system", "ours elsewhere",
                 "ours in this unit", "41.10", "load average 41.98", "5.0 s"):
        assert part in line, (part, line)
    assert "night" in runstep.allowance(PROF, MONDAY_NIGHT, IDLE).line


# --- run-step reports the allowance ----------------------------------------------------------------------------------

OCC = "l1/protein_occurrence/gene_call_version=x/stage=S0"
STEP = ("import os, pathlib, sys; s = pathlib.Path(os.environ['BIOBOX_STAGING']) / sys.argv[1]; "
        "s.parent.mkdir(parents=True); s.write_bytes(b'x' * 10); "
        "pathlib.Path(sys.argv[2]).write_text(os.environ['BIOBOX_THREADS'])")


@pytest.fixture
def small(tmp_path) -> profile.Profile:
    data = {"max_threads_day": 4, "max_threads_night": 6, "mem_scope": "1G", "tmp": str(tmp_path / "tmp"),
            "durable": str(tmp_path / "durable"), "bulk": str(tmp_path / "bulk"), "nice_monolithic": 19,
            "nice_fanout": 15, "day_start_hour": 8, "day_end_hour": 20, "cpu_min_per_process": 300,
            "small_file_floor_bytes": 1000}
    p = tmp_path / "profile.yaml"
    p.write_text(yaml.safe_dump(data))
    return profile.load(p)


def step(prof, tmp_path, name: str, when: datetime, use: runstep.CpuUse, **kw) -> tuple[int, Path]:
    root, marker = tmp_path / "root", tmp_path / f"{name}.threads"
    root.mkdir(exist_ok=True)
    code = runstep.run_step(prof, root, "S0", name, [sys.executable, "-c", STEP, f"{OCC}/{name}.parquet", str(marker)],
                            clock=lambda: when, load=lambda: 41.98, cpu=lambda: use, **kw)
    return code, marker


STAMP = r"^\S+ run-step: "


def test_run_step_reports_its_allowance_and_the_readings_in_one_line(small, tmp_path, capsys):
    code, marker = step(small, tmp_path, "a", MONDAY_NIGHT, runstep.CpuUse(others=0.01, own_other=2.2, window_s=5.0))
    assert code == 0
    assert marker.read_text() == "4"  # the night cap 6, less our 2 CPUs elsewhere
    lines = [ln for ln in capsys.readouterr().err.splitlines() if "allowance" in ln]
    assert len(lines) == 1, lines
    assert re.match(STAMP + r"S0\.a: allowance 4 threads", lines[0]), lines
    assert "ours elsewhere 2.20" in lines[0] and "load average 41.98" in lines[0]


def test_run_step_with_threads_above_the_allowance_runs_and_says_so(small, tmp_path, capsys):
    """A fixed-size job of ours (a test, a measurement, a shard) is not refused because the account's other work
    fills the cap: the launcher that can wait (build-stage) holds for it instead. run-step says the total is over."""
    code, marker = step(small, tmp_path, "b", MONDAY_NOON, runstep.CpuUse(others=0.0, own_other=3.0), threads=2)
    assert code == 0 and marker.read_text() == "2"
    err = capsys.readouterr().err
    assert "allowance 1 thread" in err and "above" in err


def test_run_step_refuses_threads_that_would_run_on_at_the_night_cap_into_the_day(small, tmp_path, capsys):
    """Card C-28 (3): 6 threads (the test profile's night cap) at 07:30 on a Monday could run 50 minutes on a 300
    CPU-minute budget, past 08:00: above the day cap of 4 that then holds, so refused; at 06:30 it is allowed."""
    code, _ = step(small, tmp_path, "c", datetime(2026, 11, 9, 7, 30), IDLE, threads=6)
    assert code == runstep.EXIT_REFUSED
    assert "--threads 6 refused" in capsys.readouterr().err
    code, marker = step(small, tmp_path, "d", datetime(2026, 11, 9, 6, 30), IDLE, threads=6)
    assert code == 0 and marker.read_text() == "6"


def test_a_skipped_step_takes_no_reading(small, tmp_path):
    assert step(small, tmp_path, "e", MONDAY_NOON, IDLE)[0] == 0

    def never() -> runstep.CpuUse:
        raise AssertionError("a step with a run record took a CPU reading")

    root = tmp_path / "root"
    assert runstep.run_step(small, root, "S0", "e", [sys.executable, "-c", "pass"], clock=lambda: MONDAY_NOON,
                            cpu=never) == 0


def test_the_module_prints_the_allowance_for_the_driver(small, tmp_path, capsys):
    p = tmp_path / "profile.yaml"
    use = runstep.CpuUse(others=0.0, own_other=1.0, window_s=5.0)
    assert runstep.main(["--profile", str(p), "--threads", "1"], clock=lambda: MONDAY_NIGHT, cpu=lambda: use) == 0
    out, err = capsys.readouterr()
    assert out == "5\n"
    assert "allowance 5" in err
    assert runstep.main(["--profile", str(p), "--threads", "1"], clock=lambda: MONDAY_NIGHT, cpu=lambda: FAILED) == 0
    out, err = capsys.readouterr()
    assert out == "4\n" and "no /proc" in err  # a failed reading: the day cap


def test_the_test_suite_never_reads_the_real_host():
    """tests/conftest.py replaces the default readers (CPU by owner, and the load average that is only reported) for
    every test in this process."""
    assert runstep.read_cpu() == runstep.CpuUse()
    assert runstep._load1() == 0.0


# --- the driver's night shard count ------------------------------------------------------------------------------------


def night_jobs(tmp_path: Path, stub: str) -> str:
    """The driver's `night_jobs` function, cut out of the script and run with a stand-in `r50_bb.sh`."""
    text = DRIVER.read_text()
    m = re.search(r"^night_jobs\(\) \{.*?^\}$", text, re.M | re.S)
    assert m, "no night_jobs function in the driver"
    bindir = tmp_path / "bin"
    bindir.mkdir(exist_ok=True)
    (bindir / "r50_bb.sh").write_text(stub)
    script = (f"BIN={bindir}\nLOG={tmp_path / 'log'}\nPROFILE={tmp_path / 'p.yaml'}\nCPUS=4\nNIGHT_JOBS=10\n"
              f"NIGHT_MAX_JOBS=14\n{m.group(0)}\nnight_jobs\n")
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, timeout=60, check=True).stdout.strip()


@pytest.mark.parametrize("stub,jobs", [
    ("echo 56\n", "14"),  # an idle night: the night cap
    ("echo 40\n", "10"),  # another user computes: the day cap
    ("echo 22\n", "5"),  # ours elsewhere use 34 of 56
    ("echo 0\n", "1"),  # nothing left: one job, which build-stage holds until it fits
    ("echo 'Traceback' >&2\nexit 1\n", "10"),  # the reading failed: the constant, the smaller count
    ("echo nonsense\n", "10"),
    ("echo 400\n", "14"),  # never above the night cap's count
])
def test_the_drivers_night_shard_count_follows_the_allowance(tmp_path, stub, jobs):
    assert night_jobs(tmp_path, stub) == jobs


def test_the_driver_switches_on_the_hours_mode_not_on_a_count_it_reads_while_its_own_build_runs():
    """During a scan the reading would count the running build (another unit) as ours elsewhere and flap; the driver
    reads the allowance only at a launch, and switches when the hour's mode (day or night) changes."""
    text = DRIVER.read_text()
    loop = re.search(r"^wait_unit\(\) \{.*?^\}$", text, re.M | re.S).group(0)
    assert "night_jobs" not in loop and "want_mode" in loop
    assert "$(want_jobs)" not in loop
    assert re.search(r"NIGHT_JOBS=\$\{NIGHT_JOBS:-10\}", text)  # the fallback constant is kept
