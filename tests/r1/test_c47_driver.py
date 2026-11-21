"""Ring R1: card C-47's driver and queue report every failure (2026-11-20).

On the pod `c47_driver.sh check` exited 0 although `c47_check.py` exited 1: the mode's status was lost behind the
driver's last `log` line, and the queue, which judged each mode by its exit status alone, ran `backup` and wrote
"all passed". Here the driver runs in a scratch directory (C47_DIR and friends) with fake `systemd-run` and
`singularity` first on PATH; nothing of the pod's /tmp/if2453/c47 is touched.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-47")  # ring R1 by directory (tests/conftest.py)

HERE = Path(__file__).resolve().parents[2] / "scripts" / "stress" / "build_stage_pod"
DRIVER = HERE / "c47_driver.sh"
QUEUE = HERE / "c47_queue.sh"

SYSTEMD_RUN = """#!/bin/sh
while [ $# -gt 0 ]; do
  case "$1" in --user|--scope|--quiet|--collect|--unit=*) shift ;; -p) shift 2 ;; *) break ;; esac
done
exec "$@"
"""
SYSTEMCTL = "#!/bin/sh\nexit 0\n"
SINGULARITY = """#!/bin/sh
echo "$*" >> "$FAKE_LOG"
case "$*" in
  *c47_check.py*)
    out=""; prev=""
    for a in "$@"; do [ "$prev" = "--out" ] && out=$a; prev=$a; done
    mkdir -p "$out"
    printf 'path\\tbytes\\tsha256\\n%s\\t1\\tab\\n' "${FAKE_ADDED:-stage=C5/x}" > "$out/root_added.tsv"
    echo "RESULT: fake"
    exit "${FAKE_CHECK_RC:-0}" ;;
  *" release verify "*) echo "release R: verified"; exit "${FAKE_VERIFY_RC:-0}" ;;
esac
exit 0
"""


def _env(tmp_path: Path, **extra) -> dict:
    fake = tmp_path / "fakebin"
    fake.mkdir(exist_ok=True)
    for name, text in (("systemd-run", SYSTEMD_RUN), ("systemctl", SYSTEMCTL), ("singularity", SINGULARITY)):
        (fake / name).write_text(text)
        (fake / name).chmod(0o755)
    c = tmp_path / "c47"
    root = tmp_path / "root"
    (root / "l0seq/x/stage=S1").mkdir(parents=True, exist_ok=True)
    (root / "l0seq/x/stage=S1/part-0.parquet").write_bytes(b"x")
    (c / "after").mkdir(parents=True, exist_ok=True)
    (c / "before").mkdir(parents=True, exist_ok=True)
    (c / "units").write_text("c47-c5-2\n")
    return {**os.environ, "PATH": f"{fake}:{os.environ['PATH']}", "FAKE_LOG": str(tmp_path / "fake.log"),
            "C47_DIR": str(c), "C47_ROOT": str(root), "C47_INPUTS": str(tmp_path / "inputs"),
            "C47_BACKUP_DIR": str(tmp_path / "stor"), "C47_V1_DIR": str(tmp_path / "v1"), **extra}


def _drive(mode: str, env: dict) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", str(DRIVER), mode], env=env, capture_output=True, text=True, timeout=120)


def test_a_failed_check_is_the_drivers_exit_status(tmp_path):
    env = _env(tmp_path, FAKE_CHECK_RC="1")
    (Path(env["C47_DIR"]) / "after" / "DONE").touch()
    got = _drive("check", env)
    assert got.returncode != 0
    assert not (Path(env["C47_DIR"]) / "check" / "DONE").exists()


def test_a_passed_check_writes_its_marker(tmp_path):
    env = _env(tmp_path)
    c = Path(env["C47_DIR"])
    (c / "after" / "DONE").touch()
    got = _drive("check", env)
    assert got.returncode == 0, got.stderr
    assert (c / "check" / "DONE").exists()
    assert (c / "root_added.tsv").read_text() == (c / "check" / "root_added.tsv").read_text()


def test_check_refuses_without_the_after_marker(tmp_path):
    env = _env(tmp_path)
    got = _drive("check", env)
    assert got.returncode != 0
    assert not (Path(env["C47_DIR"]) / "check" / "DONE").exists()


def test_a_second_check_into_its_own_directory_must_list_the_same_files(tmp_path):
    env = _env(tmp_path)
    c = Path(env["C47_DIR"])
    (c / "after" / "DONE").touch()
    assert _drive("check", env).returncode == 0
    first = (c / "check" / "check.out").read_text()
    same = _drive("check", {**env, "CHECK_DIR": str(c / "check2")})
    assert same.returncode == 0 and (c / "check2" / "DONE").exists()
    other = _drive("check", {**env, "CHECK_DIR": str(c / "check3"), "FAKE_ADDED": "stage=C5/y"})
    assert other.returncode != 0 and not (c / "check3" / "DONE").exists()
    assert (c / "check" / "check.out").read_text() == first  # the first check's record is kept as it is


def test_a_failed_verify_after_the_build_is_the_drivers_exit_status(tmp_path):
    env = _env(tmp_path, FAKE_VERIFY_RC="1")
    got = _drive("after", env)
    assert got.returncode != 0
    assert not (Path(env["C47_DIR"]) / "after" / "DONE").exists()


def test_backup_refuses_without_a_passed_check(tmp_path):
    env = _env(tmp_path)
    c = Path(env["C47_DIR"])
    (c / "root_added.tsv").write_text("path\tbytes\tsha256\n")
    got = _drive("backup", env)
    assert got.returncode != 0
    assert not (c / "backup" / "DONE").exists()
    assert not Path(env["C47_BACKUP_DIR"]).exists() or not any(Path(env["C47_BACKUP_DIR"]).iterdir())


def test_an_unknown_mode_exits_2(tmp_path):
    assert _drive("nosuchmode", _env(tmp_path)).returncode == 2


def _queue(tmp_path: Path, driver_body: str) -> subprocess.CompletedProcess:
    env = _env(tmp_path)
    v1 = Path(env["C47_V1_DIR"])
    (v1 / "after").mkdir(parents=True)
    (v1 / "after" / "DONE").touch()
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "c47_driver.sh").write_text("#!/bin/bash\n" + driver_body)
    return subprocess.run(["bash", str(QUEUE)], env={**env, "C47_BIN": str(bindir), "POLL": "1"},
                          capture_output=True, text=True, timeout=120)


def test_the_queue_stops_when_a_mode_exits_0_without_its_marker(tmp_path):
    got = _queue(tmp_path, 'echo "$1" >> "$C47_DIR/modes"\nexit 0\n')
    assert got.returncode != 0
    assert (tmp_path / "c47" / "modes").read_text().split() == ["rebefore"]
    assert "all passed" not in got.stdout


def test_the_queue_runs_every_mode_when_each_writes_its_marker(tmp_path):
    body = """C=$C47_DIR
echo "$1" >> "$C/modes"
case "$1" in
  rebefore) mkdir -p "$C/before"; touch "$C/before/REBEFORE" "$C/before/DONE" ;;
  run) mkdir -p "$C/run"; touch "$C/run/DONE" ;;
  after) mkdir -p "$C/after"; touch "$C/after/DONE" ;;
  check) mkdir -p "$C/check"; touch "$C/check/DONE" ;;
  backup) mkdir -p "$C/backup"; touch "$C/backup/DONE" ;;
esac
exit 0
"""
    got = _queue(tmp_path, body)
    assert got.returncode == 0, got.stdout
    assert (tmp_path / "c47" / "modes").read_text().split() == ["rebefore", "run", "after", "check", "backup"]
