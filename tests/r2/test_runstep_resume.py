"""C-03: `bio-box run-step` killed at each phase, then rerun (D-03; ring R2).

A crash in the publish path is simulated by a hook that raises `Crash` at a named point; `run_step`
must not catch it, so the root, `.incoming`, staging and the intent file are left exactly as a
killed process would leave them. A kill while the command runs is a real SIGTERM or SIGINT sent to
a `python -m biobox.cli run-step` process once its step has written a marker. After every recovery
the root (every file except `runs/`) and the run record's `outputs` must equal an uninterrupted run
byte for byte; the oracle is that uninterrupted run, in its own root.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path

import pytest
import yaml

from biobox import profile, runstep

pytestmark = [pytest.mark.r2, pytest.mark.covers("D-03")]

FAKE = r"""
import os, pathlib, sys
args = sys.argv[1:]
staging = pathlib.Path(os.environ["BIOBOX_STAGING"])
for flag, val in zip(args[::2], args[1::2]):
    if flag == "--out":
        path, n = val.rsplit("=", 1)
        (staging / path).parent.mkdir(parents=True, exist_ok=True)
        (staging / path).write_bytes(bytes(range(256)) * (int(n) // 256) + b"x" * (int(n) % 256))
    elif flag == "--marker":
        pathlib.Path(val).write_text("ran")
"""

# Writes one output, starts a grandchild in its own process group, records both pids, then blocks
# until the grandchild exits; the grandchild blocks opening a FIFO nobody writes.
BLOCKING = r"""
import os, pathlib, subprocess, sys
out, fifo, marker = sys.argv[1:4]
staging = pathlib.Path(os.environ["BIOBOX_STAGING"])
(staging / out).parent.mkdir(parents=True, exist_ok=True)
(staging / out).write_bytes(b"x" * 2000)
kid = subprocess.Popen([sys.executable, "-c", "import sys; open(sys.argv[1]).read()", fifo])
pathlib.Path(marker + ".tmp").write_text(f"{os.getpid()} {kid.pid}")
os.replace(marker + ".tmp", marker)
kid.wait()
"""

OCC = "l1/protein_occurrence/gene_call_version=x/stage=S0"
LEDGER = "ledger/stage=S0/step=a/part-0.parquet"
OUTS = [f"{OCC}/part-0.parquet", f"{OCC}/part-1.parquet", LEDGER]  # sorted: the rename order
SIZES = [2000, 1500, 10]
FLOOR = 1000
NOON_MONDAY = datetime(2026, 11, 9, 12, 0)
POINTS = ["after_intent", "after_rsync", "after_rename:0", "after_rename:1", "after_rename:2",
          "before_record", "after_record"]


class Crash(Exception):
    """Raised by the hook to stop `run_step` dead at a publish point."""


def step_cmd(*extra: str) -> list[str]:
    outs = [a for p, n in zip(OUTS, SIZES) for a in ("--out", f"{p}={n}")]
    return [sys.executable, "-c", FAKE, *outs, *extra]


def make_prof(d: Path) -> profile.Profile:
    data = {
        "max_threads_day": 4, "max_threads_night": 6, "mem_scope": "1G",
        "tmp": str(d / "tmp"), "durable": str(d / "durable"), "bulk": str(d / "bulk"),
        "nice_monolithic": 19, "nice_fanout": 15, "day_start_hour": 8, "day_end_hour": 20,
        "cpu_min_per_process": 300, "small_file_floor_bytes": FLOOR,
    }
    d.mkdir(parents=True, exist_ok=True)
    (d / "profile.yaml").write_text(yaml.safe_dump(data))
    return profile.load(d / "profile.yaml")


def run(prof, root, cmd, **kw) -> int:
    kw.setdefault("clock", lambda: NOON_MONDAY)
    kw.setdefault("load", lambda: 0.0)
    return runstep.run_step(prof, root, "S0", "occ", cmd, **kw)


def crash_at(point: str):
    def hook(name: str) -> None:
        if name == point:
            raise Crash(name)
    return hook


def tree(root: Path) -> dict[str, str]:
    """sha256 of every file under the root except `runs/`."""
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file() and p.relative_to(root).parts[0] != "runs"}


def root_files(root: Path) -> list[str]:
    return sorted(p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file())


@pytest.fixture
def prof(tmp_path) -> profile.Profile:
    return make_prof(tmp_path)


@pytest.fixture
def root(tmp_path) -> Path:
    r = tmp_path / "root"
    r.mkdir()
    return r


@pytest.fixture(scope="module")
def baseline(tmp_path_factory) -> tuple[dict[str, str], list[dict]]:
    """The root tree and run-record outputs of an uninterrupted run."""
    d = tmp_path_factory.mktemp("baseline")
    r = d / "root"
    r.mkdir()
    assert run(make_prof(d), r, step_cmd()) == 0
    rec = json.loads((r / "runs" / "stage=S0" / "occ.json").read_text())
    return tree(r), rec["outputs"]


def paths(prof, root):
    runs = root / "runs" / "stage=S0"
    return (runs / "occ.publishing.json", runs / "occ.json", root / ".incoming" / "S0.occ",
            Path(prof.tmp) / "staging" / "S0.occ")


def assert_like_baseline(root, baseline):
    want_tree, want_outputs = baseline
    assert tree(root) == want_tree
    assert json.loads((root / "runs" / "stage=S0" / "occ.json").read_text())["outputs"] == want_outputs


def test_hook_sees_every_point_in_order(prof, root, baseline):
    seen: list[str] = []
    assert run(prof, root, step_cmd(), _hook=seen.append) == 0
    assert seen == POINTS
    intent, record, incoming, stage_dir = paths(prof, root)
    assert not intent.exists() and record.exists()
    assert not incoming.exists() and not stage_dir.exists()
    assert root_files(root) == OUTS + ["runs/stage=S0/occ.json"]
    assert_like_baseline(root, baseline)


def test_intent_lists_every_output_and_the_measurements(prof, root):
    with pytest.raises(Crash):
        run(prof, root, step_cmd(), _hook=crash_at("after_intent"))
    intent, record, incoming, stage_dir = paths(prof, root)
    assert not record.exists() and not incoming.exists()
    got = json.loads(intent.read_text())
    assert got["outputs"] == [
        {"path": p, "bytes": n, "sha256": hashlib.sha256((stage_dir / p).read_bytes()).hexdigest()}
        for p, n in zip(OUTS, SIZES)]
    for key in ("wall_s", "cpu_min", "peak_rss_bytes", "started_at", "argv", "threads", "nice"):
        assert key in got, key
    assert got["stage"] == "S0" and got["step"] == "occ" and got["threads"] == 4 and got["nice"] == 19


@pytest.mark.parametrize("drop_staging", [False, True], ids=["staging_kept", "staging_deleted"])
@pytest.mark.parametrize("point", POINTS)
def test_crash_at_each_publish_point_then_rerun(prof, root, tmp_path, baseline, point, drop_staging):
    with pytest.raises(Crash):
        run(prof, root, step_cmd(), _hook=crash_at(point))
    intent, record, incoming, stage_dir = paths(prof, root)
    assert intent.exists()
    assert record.exists() == (point == "after_record")
    crashed_intent = json.loads(intent.read_text())
    if drop_staging:
        shutil.rmtree(stage_dir, ignore_errors=True)

    marker = tmp_path / "marker"
    assert run(prof, root, step_cmd("--marker", str(marker))) == 0
    # Only a crash before anything reached `.incoming` or the root needs the command again.
    assert marker.exists() == (point == "after_intent")
    assert_like_baseline(root, baseline)
    assert not intent.exists() and not incoming.exists() and not stage_dir.exists()
    assert not (root / ".incoming").exists()
    rec = json.loads(record.read_text())
    if point != "after_intent":  # recovered: the record is the crashed run's measurements
        assert rec == crashed_intent
    assert root_files(root) == OUTS + ["runs/stage=S0/occ.json"]

    # A third attempt is an ordinary resume.
    assert run(prof, root, step_cmd("--marker", str(tmp_path / "again"))) == 0
    assert not (tmp_path / "again").exists()


@pytest.mark.parametrize("damage", ["truncated", "altered", "missing"])
def test_partial_rsync_is_rerun(prof, root, tmp_path, baseline, damage):
    """Killed during the rsync: `.incoming` holds a partial copy, which recovery must not trust."""
    with pytest.raises(Crash):
        run(prof, root, step_cmd(), _hook=crash_at("after_rsync"))
    intent, record, incoming, stage_dir = paths(prof, root)
    bad = incoming / OUTS[1]
    if damage == "truncated":
        bad.write_bytes(bad.read_bytes()[:700])
    elif damage == "altered":  # same size, other bytes: only the sha256 can tell
        data = bad.read_bytes()
        bad.write_bytes(bytes([data[0] ^ 0xFF]) + data[1:])
    else:
        bad.unlink()
    marker = tmp_path / "marker"
    assert run(prof, root, step_cmd("--marker", str(marker))) == 0
    assert marker.exists()
    assert_like_baseline(root, baseline)
    assert not intent.exists() and not incoming.exists() and not stage_dir.exists()


@pytest.mark.parametrize("drop_staging", [False, True], ids=["staging_kept", "staging_deleted"])
def test_unrecoverable_deletes_only_our_files(prof, root, tmp_path, baseline, drop_staging):
    """Output 0 is in the root, output 1 is lost from `.incoming`: output 0 (ours, by sha256) is
    deleted and the step reruns; a file the intent does not list is never touched."""
    with pytest.raises(Crash):
        run(prof, root, step_cmd(), _hook=crash_at("after_rename:0"))
    intent, record, incoming, stage_dir = paths(prof, root)
    assert (root / OUTS[0]).exists()
    (incoming / OUTS[1]).unlink()
    if drop_staging:
        shutil.rmtree(stage_dir)
    other = root / OCC / "part-9.parquet"
    other.write_bytes(b"o" * 3000)
    unrelated = root / "notes.txt"
    unrelated.write_text("not a step output")

    marker = tmp_path / "marker"
    assert run(prof, root, step_cmd("--marker", str(marker))) == 0
    assert marker.exists()
    assert other.read_bytes() == b"o" * 3000
    assert unrelated.read_text() == "not a step output"
    other.unlink()
    unrelated.unlink()
    assert_like_baseline(root, baseline)
    assert not intent.exists() and not incoming.exists()


@pytest.mark.parametrize("point,where", [("after_intent", 1), ("after_rename:0", 2), ("after_rename:1", 0)])
def test_foreign_file_at_a_listed_path_refused(prof, root, tmp_path, capsys, point, where):
    """A root file at a listed path whose sha256 differs from the intent's is not ours: refuse,
    name it, and change nothing (not the foreign file, not our files, not the intent)."""
    with pytest.raises(Crash):
        run(prof, root, step_cmd(), _hook=crash_at(point))
    intent, record, incoming, stage_dir = paths(prof, root)
    foreign = root / OUTS[where]
    if foreign.exists():
        foreign.unlink()
    foreign.parent.mkdir(parents=True, exist_ok=True)
    foreign.write_bytes(b"F" * SIZES[where])  # same size as ours, other bytes
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest()
              for base in (root, incoming) if base.exists() for p in base.rglob("*") if p.is_file()}
    capsys.readouterr()

    marker = tmp_path / "marker"
    assert run(prof, root, step_cmd("--marker", str(marker))) == runstep.EXIT_REFUSED
    assert OUTS[where] in capsys.readouterr().err
    assert not marker.exists()
    after = {p: hashlib.sha256(p.read_bytes()).hexdigest()
             for base in (root, incoming) if base.exists() for p in base.rglob("*") if p.is_file()}
    assert after == before
    assert foreign.read_bytes() == b"F" * SIZES[where]
    assert intent.exists() and not record.exists()


def _wait_for(cond, proc: subprocess.Popen | None, what: str, timeout: float = 30.0) -> None:
    deadline = time.monotonic() + timeout
    while not cond():
        if proc is not None and proc.poll() is not None:
            pytest.fail(f"run-step exited {proc.returncode} before {what}")
        if time.monotonic() > deadline:
            pytest.fail(f"timed out waiting for {what}")
        if proc is not None:
            try:
                proc.wait(timeout=0.02)
            except subprocess.TimeoutExpired:
                pass
        else:
            threading.Event().wait(0.02)


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="needs a FIFO")
@pytest.mark.parametrize("sig", [signal.SIGTERM, signal.SIGINT], ids=["SIGTERM", "SIGINT"])
def test_signal_while_the_command_runs(tmp_path, sig):
    prof = make_prof(tmp_path)
    root = tmp_path / "root"
    root.mkdir()
    fifo, marker = tmp_path / "fifo", tmp_path / "marker"
    os.mkfifo(fifo)
    out = f"{OCC}/part-0.parquet"
    argv = [sys.executable, "-m", "biobox.cli", "run-step", "--profile", str(tmp_path / "profile.yaml"),
            "--root", str(root), "--stage", "S0", "--step", "occ", "--threads", "1",
            "--", sys.executable, "-c", BLOCKING, out, str(fifo), str(marker)]
    log = tmp_path / "run-step.log"  # a file, not a pipe: a surviving step would hold a pipe open
    with log.open("wb") as fh:
        proc = subprocess.Popen(argv, stdout=fh, stderr=subprocess.STDOUT)
    pids: list[int] = []
    try:
        _wait_for(marker.exists, proc, "the step's marker")
        pids = [int(x) for x in marker.read_text().split()]
        proc.send_signal(sig)
        proc.wait(timeout=30)
        assert proc.returncode == 128 + sig, log.read_text()
        for pid in pids:  # the step and its grandchild both got the signal
            _wait_for(lambda pid=pid: not _alive(pid), None, f"pid {pid} to exit")
        intent, record, incoming, stage_dir = paths(prof, root)
        assert root_files(root) == []
        assert not intent.exists() and not incoming.exists()
        assert (stage_dir / out).read_bytes() == b"x" * 2000  # staging is left for inspection
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()
        for pid in pids:
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass


def test_signal_handlers_restored(prof, root):
    def mine(signum, frame):
        pass
    old_term, old_int = signal.signal(signal.SIGTERM, mine), signal.signal(signal.SIGINT, mine)
    try:
        assert run(prof, root, step_cmd()) == 0
        assert signal.getsignal(signal.SIGTERM) is mine and signal.getsignal(signal.SIGINT) is mine
    finally:
        signal.signal(signal.SIGTERM, old_term)
        signal.signal(signal.SIGINT, old_int)


def test_runs_outside_the_main_thread(prof, root, baseline):
    got: list = []
    t = threading.Thread(target=lambda: got.append(run(prof, root, step_cmd())))
    t.start()
    t.join(timeout=60)
    assert got == [0]
    assert_like_baseline(root, baseline)
