"""Card C-13 (3): the wrapper shell of `scripts/pod_run.sh` writes its last log line when the unit is stopped.

`pod_run.sh` starts the unit's main process as `/bin/sh -c '<wrapper>' sh <unit> <command...>`; the wrapper runs the
command and then prints `pod_run: <unit> exited <code>`. `systemctl --user stop` sends SIGTERM to every process of
the unit's cgroup at once, the wrapper included (measured in C-10: the log of a stopped unit had no last line). Here
the wrapper is taken from the script itself and run around a command that exits 143 on SIGTERM, as `build-stage`
does, in a process group of its own; the group is sent SIGTERM as systemd sends it to the cgroup. Nothing is timed:
the test waits for the command to say it is ready, then for the wrapper to end, each with a timeout.
"""

from __future__ import annotations

import os
import re
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-13")

POD_RUN = Path(__file__).resolve().parents[2] / "scripts" / "pod_run.sh"
TIMEOUT = 30.0  # seconds

CHILD = r"""
import pathlib, signal, sys
signal.signal(signal.SIGTERM, lambda *a: sys.exit(143))
print("step output", flush=True)
pathlib.Path(sys.argv[1]).write_text("ready")
while True:
    signal.pause()
"""


def wrapper() -> str:
    """The wrapper script that pod_run.sh gives `/bin/sh -c`."""
    m = re.search(r"/bin/sh -c '([^']*)'", POD_RUN.read_text())
    assert m, "pod_run.sh starts no `/bin/sh -c '...'` wrapper"
    return m.group(1)


def wait_for(cond, what: str) -> None:
    end = time.monotonic() + TIMEOUT
    while not cond():
        assert time.monotonic() < end, f"timed out waiting for {what}"
        time.sleep(0.05)


def test_c13_the_wrapper_writes_its_last_line_when_the_unit_is_stopped(tmp_path):
    ready, log = tmp_path / "ready", tmp_path / "unit.log"
    with log.open("wb") as fh:
        proc = subprocess.Popen(["/bin/sh", "-c", wrapper(), "sh", "c13-unit", sys.executable, "-c", CHILD,
                                 str(ready)], stdout=fh, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            wait_for(ready.exists, "the command to start")
            os.killpg(proc.pid, signal.SIGTERM)  # what `systemctl --user stop` does to the unit's cgroup
            proc.wait(TIMEOUT)
        finally:
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait(TIMEOUT)
    lines = log.read_text().splitlines()
    assert lines[0] == "step output"
    assert lines[-1] == "pod_run: c13-unit exited 143", lines
    assert "pod_run: c13-unit received SIGTERM" in lines
    assert proc.returncode == 143


def test_c13_the_wrapper_still_reports_a_normal_exit(tmp_path):
    log = tmp_path / "unit.log"
    with log.open("wb") as fh:
        code = subprocess.run(["/bin/sh", "-c", wrapper(), "sh", "c13-unit", sys.executable, "-c",
                               "import sys; print('hi'); sys.exit(3)"], stdout=fh, stderr=subprocess.STDOUT,
                              timeout=TIMEOUT).returncode
    assert log.read_text().splitlines() == ["hi", "pod_run: c13-unit exited 3"]
    assert code == 3


def test_c13_python_in_the_unit_runs_unbuffered():
    """Card C-13 (2), for the wrapper: Python inside the image runs with PYTHONUNBUFFERED, so a step command's own
    output reaches the log as it is written."""
    text = POD_RUN.read_text()
    exec_line = next(ln for ln in text.splitlines() if "exec --cleanenv" in ln and not ln.lstrip().startswith("#"))
    assert "--env PYTHONUNBUFFERED=1" in exec_line
