"""Ring R1: a wave of stage builds is a script in the tree (card R-07, audit of 2026-11-25 section 2.1 (2)).

`scripts/stress/build_stage_pod/stage_wave.sh [--resume-at STAGE --from-mode MODE] LOG STAGE...` runs
`stage_queue.sh` for each stage in order, writes a dated line at each start and end, and ends at once with a non-zero
status and a line naming the stage, the mode and the status on a failure, a missing `DONE` marker after an exit 0, or a
`STOP` file of the data root. It is driven here with stub `stage_queue.sh` and `stage_driver.sh` in a copy of the
script's directory whose `/tmp/if2453` is rewritten to the test's temporary directory, in the `World` of
tests/r1/test_stage_driver.py (imported through its directory on sys.path). Nothing of the pod's `/tmp/if2453` is read
or written.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_stage_driver as sd  # noqa: E402  (the sibling module, never `from tests.r1 import`)

pytestmark = pytest.mark.covers("R-07")  # ring R1 by directory (tests/conftest.py)

SCRIPT = sd.BSP / "stage_wave.sh"
MODES = ("estimate", "rehearse", "before", "run", "after", "check", "backup")

QUEUE = """#!/bin/bash
echo "queue $1" >> "$FAKE_CALLS"
lc=$(echo "$1" | tr 'A-Z' 'a-z')
rc=0
v=FAKE_QUEUE_RC_$1
[ -n "${!v:-}" ] && rc=${!v}
if [ "$rc" = 0 ]; then
  case " ${FAKE_QUEUE_NODONE:-} " in *" $1 "*) ;; *) mkdir -p "$STAGES_HOME/$lc/backup"; echo x > "$STAGES_HOME/$lc/backup/DONE" ;; esac
fi
[ "${FAKE_STOP_AFTER:-}" = "$1" ] && touch "$FAKE_ROOT/STOP"
exit "$rc"
"""
DRIVER = """#!/bin/bash
echo "driver $1 $2" >> "$FAKE_CALLS"
lc=$(echo "$1" | tr 'A-Z' 'a-z')
rc=0
case "${FAKE_DRIVER_RC:-}" in "$2":*) rc=${FAKE_DRIVER_RC#*:} ;; esac
if [ "$rc" = 0 ] && [ "${FAKE_DRIVER_NODONE:-}" != "$2" ]; then
  mkdir -p "$STAGES_HOME/$lc/$2"; echo x > "$STAGES_HOME/$lc/$2/DONE"
fi
exit "$rc"
"""


class Wave:
    def __init__(self, tmp: Path):
        self.w = sd.World(tmp)
        self.calls_file = tmp / "calls.txt"
        self.calls_file.write_text("")
        self.log = tmp / "out" / "wave.log"
        self.stages_home = self.w.pod / "stages"
        script = self.w.script("stage_wave.sh")
        for name, text in (("stage_queue.sh", QUEUE), ("stage_driver.sh", DRIVER)):
            dst = self.w.bin / name
            if dst.is_symlink() or dst.exists():
                dst.unlink()  # World linked the real script beside the copy: never write through it
            dst.write_text(text)
        self.script = script

    def run(self, *args: str, **env: str) -> subprocess.CompletedProcess:
        e = self.w.env(FAKE_CALLS=str(self.calls_file), STAGES_HOME=str(self.stages_home), FAKE_ROOT=str(self.w.root),
                       **env)
        return subprocess.run(["bash", str(self.script), *args], env=e, capture_output=True, text=True, timeout=120)

    def calls(self) -> list[str]:
        return self.calls_file.read_text().splitlines()

    def log_lines(self) -> list[str]:
        return self.log.read_text().splitlines() if self.log.exists() else []


@pytest.fixture
def wave(tmp_path):
    return Wave(tmp_path)


def _wait_for(cond, what: str, timeout: float = 10.0) -> None:
    deadline = time.monotonic() + timeout
    while not cond():
        assert time.monotonic() < deadline, what
        time.sleep(0.05)


def test_the_script_is_in_the_tree_and_passes_both_shells_syntax_checks():
    assert SCRIPT.is_file(), f"{SCRIPT} is not in the tree"
    for shell in ("sh", "bash"):
        p = subprocess.run([shell, "-n", str(SCRIPT)], capture_output=True, text=True)
        assert p.returncode == 0, shell + ": " + p.stderr


def test_each_stage_goes_through_the_queue_in_the_order_given_with_a_start_and_an_end_line(wave):
    p = wave.run(str(wave.log), "B1", "A1", "C1")
    assert p.returncode == 0, p.stdout + p.stderr
    assert wave.calls() == ["queue B1", "queue A1", "queue C1"]
    lines = wave.log_lines()
    for s in ("B1", "A1", "C1"):
        assert sum(1 for ln in lines if f"stage {s}: start" in ln) == 1, lines
        assert sum(1 for ln in lines if f"stage {s}: end, exit 0" in ln) == 1, lines
    order = [i for s in ("B1", "A1", "C1") for i, ln in enumerate(lines) if f"stage {s}: start" in ln]
    assert order == sorted(order)
    assert "done: every stage of the wave passed" in lines[-1]
    assert all(" stage_wave: " in ln and ln[:4].isdigit() for ln in lines), lines  # dated lines


def test_the_first_failure_ends_the_wave_with_that_status_and_a_line_naming_it(wave):
    p = wave.run(str(wave.log), "A1", "B1", "C1", FAKE_QUEUE_RC_B1="7")
    assert p.returncode == 7, p.stdout + p.stderr
    assert wave.calls() == ["queue A1", "queue B1"]
    last = wave.log_lines()[-1]
    assert "STOPPED" in last and "stage B1" in last and "mode queue" in last and "exit 7" in last, last
    assert not any("stage C1" in ln for ln in wave.log_lines())


def test_a_stage_that_exits_0_without_its_backup_marker_is_a_failure_with_status_6(wave):
    p = wave.run(str(wave.log), "A1", "B1", "C1", FAKE_QUEUE_NODONE="B1")
    assert p.returncode == 6, p.stdout + p.stderr
    assert wave.calls() == ["queue A1", "queue B1"]
    last = wave.log_lines()[-1]
    assert "STOPPED" in last and "stage B1" in last and "backup" in last and "no " in last, last
    assert not any("stage B1: end" in ln for ln in wave.log_lines())


def test_a_stop_file_before_the_first_stage_starts_nothing(wave):
    (wave.w.root / "STOP").write_text("pause\n")
    p = wave.run(str(wave.log), "A1", "B1")
    assert p.returncode == 3, p.stdout + p.stderr
    assert wave.calls() == []
    last = wave.log_lines()[-1]
    assert "STOPPED" in last and "stage A1" in last and "STOP" in last and "exit 3" in last, last


def test_a_stop_file_that_appears_while_a_stage_runs_is_found_before_the_next_stage(wave):
    p = wave.run(str(wave.log), "A1", "B1", "C1", FAKE_STOP_AFTER="A1")
    assert p.returncode == 3, p.stdout + p.stderr
    assert wave.calls() == ["queue A1"]
    assert any("stage A1: end, exit 0" in ln for ln in wave.log_lines())
    last = wave.log_lines()[-1]
    assert "STOPPED" in last and "stage B1" in last and "exit 3" in last, last


def test_a_stop_file_is_polled_not_slept_for(wave):
    """A STOP written while the wave is between stages is seen without any wait in the script."""
    (wave.w.root / "STOP").write_text("x\n")
    p = wave.run(str(wave.log), "A1")
    assert p.returncode == 3
    _wait_for(lambda: wave.log.exists(), "the log was not written")


def test_resume_runs_the_modes_left_of_the_named_stage_then_the_rest_through_the_queue(wave):
    p = wave.run("--resume-at", "V1", "--from-mode", "run", str(wave.log), "V1", "PL2", "PL3")
    assert p.returncode == 0, p.stdout + p.stderr
    assert wave.calls() == ["driver V1 run", "driver V1 after", "driver V1 check", "driver V1 backup",
                            "queue PL2", "queue PL3"]
    lines = wave.log_lines()
    for m in ("run", "after", "check", "backup"):
        assert any(f"stage V1 mode {m}: end, exit 0, marker present" in ln for ln in lines), lines
    assert any("stage V1: end (resumed from run)" in ln for ln in lines)


def test_resume_from_the_first_mode_runs_every_mode(wave):
    p = wave.run("--resume-at", "V1", "--from-mode", "estimate", str(wave.log), "V1")
    assert p.returncode == 0, p.stdout + p.stderr
    assert wave.calls() == [f"driver V1 {m}" for m in MODES]


def test_resume_stops_at_a_mode_that_fails_with_its_status(wave):
    p = wave.run("--resume-at", "V1", "--from-mode", "run", str(wave.log), "V1", "PL2", FAKE_DRIVER_RC="check:4")
    assert p.returncode == 4, p.stdout + p.stderr
    assert wave.calls() == ["driver V1 run", "driver V1 after", "driver V1 check"]
    last = wave.log_lines()[-1]
    assert "STOPPED" in last and "stage V1" in last and "mode check" in last and "exit 4" in last, last


def test_resume_requires_each_modes_marker_even_after_an_exit_0(wave):
    p = wave.run("--resume-at", "V1", "--from-mode", "run", str(wave.log), "V1", "PL2", FAKE_DRIVER_NODONE="after")
    assert p.returncode == 6, p.stdout + p.stderr
    assert wave.calls() == ["driver V1 run", "driver V1 after"]
    last = wave.log_lines()[-1]
    assert "STOPPED" in last and "stage V1" in last and "mode after" in last and "no " in last, last


def test_resume_looks_for_a_stop_file_before_each_mode(wave):
    (wave.w.root / "STOP").write_text("x\n")
    p = wave.run("--resume-at", "V1", "--from-mode", "run", str(wave.log), "V1", "PL2")
    assert p.returncode == 3, p.stdout + p.stderr
    assert wave.calls() == []
    last = wave.log_lines()[-1]
    assert "STOPPED" in last and "stage V1" in last and "mode run" in last, last


@pytest.mark.parametrize("args,why", [
    (("{log}", "A1", "B1", "A1"), "twice"),
    (("{log}", "A1", "b1"), "not a stage name"),
    (("{log}", "A1", "B 1"), "not a stage name"),
    (("{log}",), "no stage listed"),
    (("--resume-at", "Z9", "--from-mode", "run", "{log}", "A1", "B1"), "must be the first stage"),
    (("--resume-at", "B1", "--from-mode", "run", "{log}", "A1", "B1"), "must be the first stage"),
    (("--resume-at", "A1", "--from-mode", "bogus", "{log}", "A1"), "not one of"),
    (("--resume-at", "A1", "{log}", "A1"), "needs --from-mode"),
    (("--from-mode", "run", "{log}", "A1"), "needs --resume-at"),
])
def test_a_bad_list_is_refused_before_anything_runs(wave, args, why):
    p = wave.run(*[a.format(log=wave.log) for a in args])
    assert p.returncode == 2, p.stdout + p.stderr
    assert why in p.stderr, p.stderr
    assert wave.calls() == []
    assert wave.log_lines() == [], "a refused list wrote to the log"


def test_no_argument_prints_usage_and_exits_2(wave):
    p = wave.run()
    assert p.returncode == 2
    assert "usage" in p.stderr


def _env_stubs(wave: Wave) -> None:
    """The stubs, each also writing the wave's name and the share switch it was given (card A-30)."""
    for name, text in (("stage_queue.sh", QUEUE), ("stage_driver.sh", DRIVER)):
        dst = wave.w.bin / name
        dst.unlink()
        dst.write_text(text.replace('>> "$FAKE_CALLS"', '>> "$FAKE_CALLS"\necho "env $1 wave=${STAGE_WAVE:-} '
                                    'share=${STAGE_SHARE-unset}" >> "$FAKE_CALLS"', 1))


def _env_lines(wave: Wave) -> list[str]:
    return [ln for ln in wave.calls() if ln.startswith("env ")]


@pytest.mark.covers("A-30")
def test_a30_every_stage_of_one_wave_gets_the_same_new_wave_name_and_the_log_names_it(wave):
    _env_stubs(wave)
    p = wave.run(str(wave.log), "A1", "B1", STAGE_WAVE="inherited")
    assert p.returncode == 0, p.stdout + p.stderr
    names = {ln.split()[2] for ln in _env_lines(wave)}
    assert len(names) == 1 and len(_env_lines(wave)) == 2, wave.calls()
    name = names.pop().removeprefix("wave=")
    assert re.fullmatch(r"w[0-9a-f]{16}", name), name
    assert name != "inherited", "a wave never takes a name from the environment"
    assert f"wave {name}" in wave.log_lines()[0] and "sharing on" in wave.log_lines()[0], wave.log_lines()[0]
    wave.calls_file.write_text("")
    p = wave.run(str(wave.log), "C1")
    assert p.returncode == 0
    again = _env_lines(wave)[0].split()[2].removeprefix("wave=")
    assert re.fullmatch(r"w[0-9a-f]{16}", again) and again != name, "each wave has a name of its own"


@pytest.mark.covers("A-30")
def test_a30_the_modes_of_a_resumed_stage_run_in_the_wave_too(wave):
    _env_stubs(wave)
    p = wave.run("--resume-at", "V1", "--from-mode", "after", str(wave.log), "V1", "PL2")
    assert p.returncode == 0, p.stdout + p.stderr
    names = {ln.split()[2] for ln in _env_lines(wave)}
    assert len(names) == 1 and len(_env_lines(wave)) == 4, wave.calls()  # after, check, backup; then PL2's queue
    assert re.fullmatch(r"wave=w[0-9a-f]{16}", names.pop())


@pytest.mark.covers("A-30")
def test_a30_the_share_switch_reaches_every_stage_and_the_log_says_so(wave):
    _env_stubs(wave)
    p = wave.run(str(wave.log), "A1", "B1", STAGE_SHARE="0")
    assert p.returncode == 0, p.stdout + p.stderr
    assert all(ln.endswith(" share=0") for ln in _env_lines(wave)), wave.calls()
    assert "sharing off (STAGE_SHARE=0)" in wave.log_lines()[0], wave.log_lines()[0]


def test_the_wave_writes_nothing_outside_its_log_and_the_stages_own_directories(wave):
    before = sorted(str(q.relative_to(wave.w.tmp)) for q in wave.w.root.rglob("*"))
    p = wave.run(str(wave.log), "A1", "B1")
    assert p.returncode == 0, p.stdout + p.stderr
    after = sorted(str(q.relative_to(wave.w.tmp)) for q in wave.w.root.rglob("*"))
    assert before == after, "the wave changed the data root"
