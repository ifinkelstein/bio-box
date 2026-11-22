"""Ring R1: card A-19 — DEC-64 rate-guard stop messages and resilient systemd waits in rehearsal and waiters.

1. Transient `systemctl` failures (exit 1, empty stdout / empty ActiveState, before the unit's exit line or
   explicit inactive/failed state) must not prematurely end:
   - `release_train.sh` `rh_launch` (rehearsal unit wait),
   - `scripts/orchestrator/rel_wait.sh` (queue wait, suite shard wait, and `build_units` list-units check),
   - `scripts/stress/build_stage_pod/stage_queue.sh` (`AFTER` queue check),
   while `rh_launch` and `rel_wait.sh` still fail loudly when no readable answer arrives within `WATCH_LIMIT_S`.
2. DEC-64 stop messages in `stage_check.py`, `stage_driver.sh`, and `release_train.sh` name `>50%` / `more than 50%`
   and `DEC-64` while keeping all JSON keys (`over_20pct`, `recorded_over_20pct`, `outside_20pct`, `under_20pct`,
   `accepted_over_20pct`) unchanged.
"""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("A-19", "A-16", "A-12", "C-46")  # ring R1 by its directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
ORCH = REPO / "scripts" / "orchestrator"

NOOP = "#!/bin/sh\nexit 0\n"
NICE = '#!/bin/sh\n[ "$1" = -n ] && shift 2\nexec "$@"\n'
SYSTEMD_RUN = """#!/bin/sh
while [ $# -gt 0 ]; do
  case "$1" in --user|--scope|--quiet|--collect|--unit=*) shift ;; -p) shift 2 ;; *) break ;; esac
done
exec "$@"
"""


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _write_exec(path: Path, text: str) -> None:
    path.write_text(text)
    path.chmod(0o755)


def _lines(path: Path, *rows: str) -> str:
    path.write_text("".join(r + "\n" for r in rows))
    return str(path)


# --------------------------------------------------------------------------------------------------------------
# 1. release_train.sh rh_launch: resilient wait on transient systemctl failure + WATCH_LIMIT_S bound

RH_SYSTEMCTL = """#!/bin/sh
echo "systemctl $*" >> "$FAKE_LOG"
d=$(dirname "$FAKE_LOG")
last=""
for a in "$@"; do last=$a; done
finish() {
  f="$d/pod.$1"
  [ -f "$f" ] || return 0
  read -r log rc < "$f"
  echo "pod_run: $1 exited $rc" >> "$log"
  rm -f "$f"
}
state() {
  if [ -n "${FAKE_STATES:-}" ] && [ -s "$FAKE_STATES" ]; then
    head -n 1 "$FAKE_STATES"
    sed -i 1d "$FAKE_STATES"
  else
    echo inactive
  fi
}
case "$*" in
  *is-active*)
    s=$(state)
    case "$s" in
      active) exit 0 ;;
      fail) exit 1 ;;
      empty) exit 0 ;;
    esac
    finish "$last"
    exit 3 ;;
  *"show -p ActiveState"*)
    s=$(state)
    [ "$s" = fail ] && exit 1
    [ "$s" = empty ] && exit 0
    [ "$s" = inactive ] && finish "$last"
    echo "$s"
    exit 0 ;;
  *" stop "*)
    exit 0 ;;
esac
exit 0
"""

RH_POD_RUN = """#!/bin/sh
log=""; unit=""
while [ $# -gt 0 ]; do
  case "$1" in
    --log) log=$2; shift 2 ;;
    --*) shift 2 ;;
    *) unit=$1; shift; break ;;
  esac
done
echo "pod_run $unit $*" >> "$FAKE_LOG"
mkdir -p "$(dirname "$log")"
echo "pod_run: $unit started" >> "$log"
echo "$log 0" > "$(dirname "$FAKE_LOG")/pod.$unit"
exit 0
"""

RH_SINGULARITY = """#!/bin/sh
echo "singularity $*" >> "$FAKE_LOG"
# Stop cleanly right after rh_launch finishes, when rehearse runs dry-run
echo "rehearse-dry-run-reached" >&2
exit 42
"""


def _setup_rh_world(tmp_path: Path):
    pod = tmp_path / "if2453"
    fake = tmp_path / "fakebin"
    tree = tmp_path / "tree"
    work = pod / "c46" / "train"
    rh = work / "rh"
    for d in (fake, tree / "scripts", tree / "profiles", rh / "inputs", rh / "root", rh / "logs",
              pod / "r10" / "root"):
        d.mkdir(parents=True, exist_ok=True)
    for name, text in (("systemctl", RH_SYSTEMCTL), ("sleep", NOOP), ("systemd-run", SYSTEMD_RUN),
                       ("nice", NICE), ("singularity", RH_SINGULARITY)):
        _write_exec(fake / name, text)
    _write_exec(tree / "scripts" / "pod_run.sh", RH_POD_RUN)
    (tree / "profiles" / "wilkcomp02.yaml").write_text("tmp: /x\nmax_threads_day: 40\n")
    # Pre-populate inputs and S1/S2/S3 sentinel records so do_rehearse invokes only rh_launch c46-rh-s0
    (rh / "inputs" / "objects-S2.json").write_text("[]\n")
    for rel_path in ("runs/stage=r26.11.1/catalog.json", "runs/stage=S2/tapes.json", "runs/stage=S3/ingest.json"):
        p = rh / "root" / rel_path
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("{}\n")
    drv = tmp_path / "release_train.sh"
    drv.write_text((BSP / "release_train.sh").read_text().replace("/tmp/if2453", str(pod)))
    log = tmp_path / "fake.log"
    log.write_text("")
    return pod, fake, tree, work, drv, log


def test_release_train_rh_launch_survives_transient_systemctl_failure(tmp_path):
    pod, fake, tree, work, drv, log = _setup_rh_world(tmp_path)
    states = _lines(tmp_path / "states", "fail", "empty", "active", "fail", "inactive")
    env = {**os.environ, "PATH": f"{fake}:{os.environ['PATH']}", "FAKE_LOG": str(log),
           "FAKE_STATES": states, "TREE": str(tree), "WORK": str(work),
           "ROOT": str(pod / "r10" / "root")}
    got = subprocess.run(["bash", str(drv), "rehearse"], env=env, capture_output=True, text=True, timeout=60)
    # rh_launch c46-rh-s0 must NOT fail with "exited ''" on the transient systemctl failure;
    # execution must advance past rh_launch into the dry-run step (which our stub marks).
    assert "unit c46-rh-s0 exited ''" not in got.stderr, got.stderr
    assert "rehearse: dry-run before" in got.stderr, got.stderr
    calls = log.read_text()
    assert calls.count("show -p ActiveState --value c46-rh-s0") == 5, calls
    assert "no answer from systemd about c46-rh-s0" in (work / "driver.log").read_text()


def test_release_train_rh_launch_fails_loudly_at_watch_limit(tmp_path):
    pod, fake, tree, work, drv, log = _setup_rh_world(tmp_path)
    states = _lines(tmp_path / "states", *["fail"] * 20)
    env = {**os.environ, "PATH": f"{fake}:{os.environ['PATH']}", "FAKE_LOG": str(log),
           "FAKE_STATES": states, "WATCH_LIMIT_S": "0", "TREE": str(tree), "WORK": str(work),
           "ROOT": str(pod / "r10" / "root")}
    got = subprocess.run(["bash", str(drv), "rehearse"], env=env, capture_output=True, text=True, timeout=60)
    assert got.returncode != 0
    assert "no readable answer from systemd about unit c46-rh-s0" in got.stderr, got.stderr
    assert "its state is unknown and it may still be running" in got.stderr, got.stderr


# --------------------------------------------------------------------------------------------------------------
# 2. rel_wait.sh: queue wait, suite shard wait, and build_units resilient to transient systemctl failures

REL_WAIT_SYSTEMCTL = """#!/bin/sh
echo "systemctl $*" >> "$FAKE_LOG"
last=""
for a in "$@"; do last=$a; done
next_line() {
  f=$1
  if [ -n "$f" ] && [ -s "$f" ]; then
    head -n 1 "$f"
    sed -i 1d "$f"
  else
    echo "$2"
  fi
}
case "$*" in
  *is-active*)
    s=$(next_line "${FAKE_Q_STATES:-}" inactive)
    case "$s" in
      active) exit 0 ;;
      fail) exit 1 ;;
      empty) exit 0 ;;
      inactive)
        [ -n "${FAKE_MARK_ON_INACTIVE:-}" ] && mkdir -p "$(dirname "$FAKE_MARK_ON_INACTIVE")" && touch "$FAKE_MARK_ON_INACTIVE"
        exit 3 ;;
    esac
    exit 3 ;;
  *"show -p ActiveState"*)
    s=$(next_line "${FAKE_Q_STATES:-}" inactive)
    [ "$s" = fail ] && exit 1
    [ "$s" = empty ] && exit 0
    if [ "$s" = inactive ] && [ -n "${FAKE_MARK_ON_INACTIVE:-}" ]; then
      mkdir -p "$(dirname "$FAKE_MARK_ON_INACTIVE")"
      touch "$FAKE_MARK_ON_INACTIVE"
    fi
    echo "$s"
    exit 0 ;;
  *list-units*"ci-"*)
    s=$(next_line "${FAKE_SUITE_STATES:-}" done)
    [ "$s" = fail ] && exit 1
    if [ "$s" = running ]; then
      echo "ci-rel1120-0.service loaded active running x"
      exit 0
    fi
    if [ -n "${FAKE_SUITE_FINISH_DIR:-}" ]; then
      for k in 0 1 2 3 4 5 6 7; do
        echo "10 passed in 1s" > "$FAKE_SUITE_FINISH_DIR/log.$k"
      done
    fi
    exit 0 ;;
  *list-units*)
    s=$(next_line "${FAKE_BUILD_STATES:-}" empty)
    [ "$s" = fail ] && exit 1
    [ "$s" != empty ] && echo "$s loaded active running x"
    exit 0 ;;
esac
exit 0
"""


def _setup_rel_wait(tmp_path: Path, *, with_backup_done: bool = True, suite_ready: bool = True):
    pod = tmp_path / "if2453"
    fake = tmp_path / "fakebin"
    drv_dir = pod / "stages" / "c1"
    sd = pod / "ci-tmp" / "rel1120-detached"
    tree = pod / "ci" / "rel1120"
    stages_dir = tree / "scripts" / "stress" / "build_stage_pod" / "stages"
    wait_dir = pod / "c46" / "wait"
    for d in (fake, drv_dir, sd, stages_dir, wait_dir, pod / "r10" / "root"):
        d.mkdir(parents=True, exist_ok=True)
    for name, text in (("systemctl", REL_WAIT_SYSTEMCTL), ("sleep", NOOP)):
        _write_exec(fake / name, text)
    (stages_dir / "C1.conf").write_text("STAGE=C1\n")
    _write_exec(
        tree / "scripts" / "stress" / "build_stage_pod" / "release_train.sh",
        '#!/bin/bash\necho "release_train $*" >> "$FAKE_LOG"\n[ "$1" = --dry-run ] && echo "Action: cut"\nexit 0\n',
    )
    for m in ("before", "run", "after", "check"):
        (drv_dir / m).mkdir(parents=True, exist_ok=True)
        (drv_dir / m / "DONE").touch()
    if with_backup_done:
        (drv_dir / "backup").mkdir(parents=True, exist_ok=True)
        (drv_dir / "backup" / "DONE").touch()
    (drv_dir / "backup.tsv").write_text("x\t1\tab\n")
    (drv_dir / "queue.out").write_text("stage_queue C1: running\n")
    for k in range(8):
        (sd / f"log.{k}").write_text("10 passed in 1s\n" if suite_ready else "running...\n")
    sh = tmp_path / "rel_wait.sh"
    sh.write_text((ORCH / "rel_wait.sh").read_text().replace("/tmp/if2453", str(pod)))
    log = tmp_path / "fake.log"
    log.write_text("")
    return pod, fake, drv_dir, sd, sh, log


def test_rel_wait_queue_survives_transient_systemctl_failure(tmp_path):
    pod, fake, drv_dir, sd, sh, log = _setup_rel_wait(tmp_path, with_backup_done=False)
    q_states = _lines(tmp_path / "qstates", "fail", "empty", "active", "fail", "inactive")
    env = {**os.environ, "PATH": f"{fake}:{os.environ['PATH']}", "FAKE_LOG": str(log),
           "FAKE_Q_STATES": q_states, "FAKE_MARK_ON_INACTIVE": str(drv_dir / "backup" / "DONE")}
    got = subprocess.run(["bash", str(sh), str(drv_dir), "rel1120"], env=env, capture_output=True, text=True, timeout=60)
    assert got.returncode == 0, got.stdout + got.stderr
    assert "release_train run" in log.read_text()
    assert log.read_text().count("show -p ActiveState --value c1-queue") == 5


def test_rel_wait_queue_fails_loudly_at_watch_limit(tmp_path):
    pod, fake, drv_dir, sd, sh, log = _setup_rel_wait(tmp_path, with_backup_done=False)
    q_states = _lines(tmp_path / "qstates", *["fail"] * 20)
    env = {**os.environ, "PATH": f"{fake}:{os.environ['PATH']}", "FAKE_LOG": str(log),
           "FAKE_Q_STATES": q_states, "WATCH_LIMIT_S": "0"}
    got = subprocess.run(["bash", str(sh), str(drv_dir), "rel1120"], env=env, capture_output=True, text=True, timeout=60)
    assert got.returncode != 0
    assert "release_train" not in log.read_text()
    assert "no readable answer from systemd" in got.stdout, got.stdout


def test_rel_wait_suite_wait_survives_transient_list_units_failure(tmp_path):
    pod, fake, drv_dir, sd, sh, log = _setup_rel_wait(tmp_path, with_backup_done=True, suite_ready=False)
    s_states = _lines(tmp_path / "sstates", "fail", "running", "fail", "done")
    env = {**os.environ, "PATH": f"{fake}:{os.environ['PATH']}", "FAKE_LOG": str(log),
           "FAKE_SUITE_STATES": s_states, "FAKE_SUITE_FINISH_DIR": str(sd)}
    got = subprocess.run(["bash", str(sh), str(drv_dir), "rel1120"], env=env, capture_output=True, text=True, timeout=60)
    assert got.returncode == 0, got.stdout + got.stderr
    assert "release_train run" in log.read_text()


def test_rel_wait_suite_wait_fails_loudly_at_watch_limit(tmp_path):
    pod, fake, drv_dir, sd, sh, log = _setup_rel_wait(tmp_path, with_backup_done=True, suite_ready=False)
    s_states = _lines(tmp_path / "sstates", *["fail"] * 20)
    env = {**os.environ, "PATH": f"{fake}:{os.environ['PATH']}", "FAKE_LOG": str(log),
           "FAKE_SUITE_STATES": s_states, "WATCH_LIMIT_S": "0"}
    got = subprocess.run(["bash", str(sh), str(drv_dir), "rel1120"], env=env, capture_output=True, text=True, timeout=60)
    assert got.returncode != 0
    assert "release_train" not in log.read_text()
    assert "no readable answer from systemd" in got.stdout, got.stdout


def test_rel_wait_refuses_to_launch_when_build_units_list_fails(tmp_path):
    pod, fake, drv_dir, sd, sh, log = _setup_rel_wait(tmp_path, with_backup_done=True, suite_ready=True)
    b_states = _lines(tmp_path / "bstates", "fail")
    env = {**os.environ, "PATH": f"{fake}:{os.environ['PATH']}", "FAKE_LOG": str(log),
           "FAKE_BUILD_STATES": b_states}
    got = subprocess.run(["bash", str(sh), str(drv_dir), "rel1120"], env=env, capture_output=True, text=True, timeout=60)
    assert got.returncode == 6, got.stdout + got.stderr
    assert "release_train" not in log.read_text()


# --------------------------------------------------------------------------------------------------------------
# 3. stage_queue.sh: AFTER queue check ignores stale FAILED/REFUSED on transient systemctl failure

QUEUE_DRIVER_STUB = """#!/bin/bash
C=$STAGES_HOME/$(echo "$1" | tr 'A-Z' 'a-z')
case "$2" in
  info) echo "AFTER=V2"; echo "BUSY="; exit 0 ;;
esac
mkdir -p "$C/$2"
echo "$2" >> "$C/modes"
touch "$C/$2/DONE"
exit 0
"""

QUEUE_SYSTEMCTL = """#!/bin/sh
echo "systemctl $*" >> "$FAKE_LOG"
state() {
  if [ -n "${FAKE_STATES:-}" ] && [ -s "$FAKE_STATES" ]; then
    head -n 1 "$FAKE_STATES"
    sed -i 1d "$FAKE_STATES"
  else
    echo inactive
  fi
}
case "$*" in
  *is-active*)
    s=$(state)
    case "$s" in
      active) exit 0 ;;
      fail) exit 1 ;;
      empty) exit 0 ;;
      done_after)
        mkdir -p "$STAGES_HOME/v2/backup"
        touch "$STAGES_HOME/v2/backup/DONE"
        exit 3 ;;
    esac
    exit 3 ;;
  *"show -p ActiveState"*)
    s=$(state)
    [ "$s" = fail ] && exit 1
    [ "$s" = empty ] && exit 0
    if [ "$s" = done_after ]; then
      mkdir -p "$STAGES_HOME/v2/backup"
      touch "$STAGES_HOME/v2/backup/DONE"
      echo inactive
      exit 0
    fi
    echo "$s"
    exit 0 ;;
esac
exit 0
"""


def test_stage_queue_after_check_survives_transient_systemctl_failure_with_stale_failed_file(tmp_path):
    home = tmp_path / "stages"
    (home / "v2" / "before").mkdir(parents=True)
    # Stale REFUSED file left from an earlier attempt while v2-queue is still running
    (home / "v2" / "before" / "REFUSED").write_text("stale refusal\n")
    b = tmp_path / "qbin"
    b.mkdir()
    _write_exec(b / "stage_driver.sh", QUEUE_DRIVER_STUB)
    q = b / "stage_queue.sh"
    q.write_text((BSP / "stage_queue.sh").read_text().replace("/tmp/if2453", str(tmp_path / "pod")))
    (tmp_path / "pod" / "r10" / "root").mkdir(parents=True)
    fake = tmp_path / "fakebin"
    fake.mkdir()
    for name, text in (("systemctl", QUEUE_SYSTEMCTL), ("sleep", NOOP)):
        _write_exec(fake / name, text)
    states = _lines(tmp_path / "states", "fail", "empty", "active", "done_after")
    log = tmp_path / "fake.log"
    log.write_text("")
    env = {**os.environ, "PATH": f"{fake}:{os.environ['PATH']}", "FAKE_LOG": str(log),
           "FAKE_STATES": states, "STAGES_HOME": str(home), "POLL": "0", "MAX_WAIT_H": "48"}
    got = subprocess.run(["bash", str(q), "ZZ1"], env=env, capture_output=True, text=True, timeout=60)
    assert got.returncode == 0, got.stdout + got.stderr
    assert (home / "zz1" / "modes").read_text().split() == [
        "estimate", "rehearse", "before", "run", "after", "check", "backup",
    ]
    assert log.read_text().count("show -p ActiveState --value v2-queue") >= 3


# --------------------------------------------------------------------------------------------------------------
# 4. DEC-64 rate-guard stop messages (and JSON keys unchanged)

def test_dec64_stop_messages_in_stage_check_stage_driver_and_release_train():
    sc = (BSP / "stage_check.py").read_text()
    sd = (BSP / "stage_driver.sh").read_text()
    rt = (BSP / "release_train.sh").read_text()

    assert "rates over +50% (DEC-64):" in sc
    assert "rates over +20% (DEC-53):" not in sc
    assert "no step more than 50% over its estimate" in sc

    assert "more than 50% over its estimate, DEC-64" in sd
    assert "more than 20% over its estimate, DEC-53" not in sd

    assert rt.count(">50% over its estimate, DEC-43, DEC-64") == 2
    assert ">20% over its estimate, DEC-43, DEC-53" not in rt


def test_stage_check_emits_dec64_message_while_keeping_json_keys_unchanged(tmp_path, monkeypatch):
    from biobox.l1 import crispr, genecall

    monkeypatch.setattr(crispr, "minced_version", lambda: "0.4.2")
    td = _load("a19_stage_driver_tests", REPO / "tests" / "r1" / "test_stage_driver.py")
    batch = f"genecall.{genecall.batch_name(0)}"

    def records(S):
        S["records"]["ZQ"][batch] = []
        S["cpu"][batch] = 151.0
        S["gcplan"]["batches"][0]["estimate_cpu_min"] = 100.0

    r = td.build_root(tmp_path, records)
    rc, summary = td.run_check(r, monkeypatch)
    assert rc == 1
    assert any(f"rates over +50% (DEC-64): ['{batch}']" in msg for msg in summary["failed"]), summary["failed"]
    rates = summary["steps"]["rates"]
    for key in ("over_20pct", "recorded_over_20pct", "under_20pct", "accepted_over_20pct"):
        assert key in rates, f"JSON key {key} must remain unchanged"
    assert rates["over_20pct"] == [batch]
