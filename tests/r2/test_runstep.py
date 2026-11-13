"""C-02: `bio-box run-step` with fake step commands (ring R2 by directory).

The fake step is `python -c FAKE ...`. It writes the files named by `--out PATH=BYTES` under
`$BIOBOX_STAGING` and can also dump its environment, touch a marker, create `<root>/STOP`, exit
with a code, or kill itself. Output files are filler bytes, not valid Parquet: the publish rules
under test (claim by path, append-only, small-file floor) look only at paths and sizes.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import threading
from datetime import datetime
from pathlib import Path

import pytest
import yaml

from biobox import catalog, cli, profile, runstep

FAKE = r"""
import json, os, pathlib, signal, sys
args = sys.argv[1:]
staging = pathlib.Path(os.environ["BIOBOX_STAGING"])
code = 0
i = 0
while i < len(args):
    flag, val = args[i], (args[i + 1] if i + 1 < len(args) else None)
    i += 2
    if flag == "--out":
        path, n = val.rsplit("=", 1)
        (staging / path).parent.mkdir(parents=True, exist_ok=True)
        (staging / path).write_bytes(bytes(range(256)) * (int(n) // 256) + b"x" * (int(n) % 256))
    elif flag == "--env":
        pathlib.Path(val).write_text(json.dumps(dict(os.environ)))
    elif flag == "--marker":
        pathlib.Path(val).write_text("ran")
    elif flag == "--stop":
        (pathlib.Path(os.environ["BIOBOX_ROOT"]) / "STOP").write_text("")
        i -= 1
    elif flag == "--exit":
        code = int(val)
    elif flag == "--kill":
        os.kill(os.getpid(), signal.SIGKILL)
    elif flag == "--ledger":
        from biobox import ledger
        ledger.write(staging, "S0", val, [])
sys.exit(code)
"""

OCC = "l1/protein_occurrence/gene_call_version=x/stage=S0"
FLOOR = 1000
NOON_MONDAY = datetime(2026, 11, 9, 12, 0)  # a Monday
NIGHT_MONDAY = datetime(2026, 11, 9, 23, 0)
NOON_SATURDAY = datetime(2026, 11, 14, 12, 0)


def fake(*args: str) -> list[str]:
    return [sys.executable, "-c", FAKE, *args]


@pytest.fixture
def prof(tmp_path) -> profile.Profile:
    data = {
        "max_threads_day": 4, "max_threads_night": 6, "mem_scope": "1G",
        "tmp": str(tmp_path / "tmp"), "durable": str(tmp_path / "durable"), "bulk": str(tmp_path / "bulk"),
        "nice_monolithic": 19, "nice_fanout": 15, "day_start_hour": 8, "day_end_hour": 20,
        "cpu_min_per_process": 300, "small_file_floor_bytes": FLOOR,
    }
    p = tmp_path / "profile.yaml"
    p.write_text(yaml.safe_dump(data))
    return profile.load(p)


@pytest.fixture
def root(tmp_path) -> Path:
    r = tmp_path / "root"
    r.mkdir()
    return r


def run(prof, root, step, cmd, **kw) -> int:
    kw.setdefault("clock", lambda: NOON_MONDAY)
    kw.setdefault("load", lambda: 0.0)
    return runstep.run_step(prof, root, "S0", step, cmd, **kw)


def root_files(root: Path) -> list[str]:
    return sorted(p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file())


def staging(prof, step) -> Path:
    return Path(prof.tmp) / "staging" / f"S0.{step}"


def record_path(root, step) -> Path:
    return root / "runs" / "stage=S0" / f"{step}.json"


def test_outputs_published_and_record_written(prof, root, monkeypatch):
    monkeypatch.setenv("BIOBOX_IMAGE_DIGEST", "sha256:abc")
    a, b = f"{OCC}/part-0.parquet", f"{OCC}/part-1.parquet"
    assert run(prof, root, "occ", fake("--out", f"{a}=2000", "--out", f"{b}=1500")) == 0
    assert root_files(root) == [a, b, "runs/stage=S0/occ.json"]
    rec = json.loads(record_path(root, "occ").read_text())
    want = [{"path": p, "bytes": (root / p).stat().st_size,
             "sha256": hashlib.sha256((root / p).read_bytes()).hexdigest()} for p in (a, b)]
    assert rec["outputs"] == want
    assert want[0]["bytes"] == 2000
    assert rec["stage"] == "S0" and rec["step"] == "occ"
    assert rec["threads"] == 4 and rec["nice"] == 19 and rec["exit_code"] == 0
    assert rec["image_digest"] == "sha256:abc"
    assert rec["argv"][-6:] == fake("--out", f"{a}=2000", "--out", f"{b}=1500")[-6:]
    assert rec["wall_s"] >= 0 and rec["cpu_min"] > 0 and rec["peak_rss_bytes"] > 1 << 20
    assert rec["over_cpu_budget"] is False
    assert rec["host"] and rec["started_at"].endswith("+00:00")
    assert not staging(prof, "occ").exists()
    assert not (root / ".incoming" / "S0.occ").exists()


def test_image_digest_defaults_to_none(prof, root, monkeypatch):
    monkeypatch.delenv("BIOBOX_IMAGE_DIGEST", raising=False)
    assert run(prof, root, "occ", fake("--out", f"{OCC}/part-0.parquet=10")) == 0
    assert json.loads(record_path(root, "occ").read_text())["image_digest"] == "none"


def test_record_names_the_code_that_ran_it(prof, root):
    """Card C-23: the run record's meta has the code identity (`biobox.codeid`), beside the image digest."""
    from biobox import codeid

    assert run(prof, root, "occ", fake("--out", f"{OCC}/part-0.parquet=10")) == 0
    rec = json.loads(record_path(root, "occ").read_text())
    want = codeid.code_identity()
    assert (rec["code_id"], rec["code_version"], rec["code_commit"]) == (
        want["code_id"], want["code_version"], want["code_commit"])
    assert len(rec["code_id"]) == 64


def test_over_cpu_budget_is_a_warning(prof, root, tmp_path):
    p = tmp_path / "zero_budget.yaml"
    data = {f: getattr(prof, f) for f in prof.__dataclass_fields__}
    data["cpu_min_per_process"] = 0
    p.write_text(yaml.safe_dump(data))
    assert run(profile.load(p), root, "occ", fake("--out", f"{OCC}/part-0.parquet=10")) == 0
    assert json.loads(record_path(root, "occ").read_text())["over_cpu_budget"] is True


def test_env_reaches_the_step(prof, root, tmp_path):
    env_file = tmp_path / "env.json"
    assert run(prof, root, "occ", fake("--env", str(env_file), "--out", f"{OCC}/part-0.parquet=10"),
               clock=lambda: NIGHT_MONDAY) == 0
    env = json.loads(env_file.read_text())
    assert env["BIOBOX_STAGING"] == str(staging(prof, "occ"))
    assert env["BIOBOX_ROOT"] == str(root.resolve())
    for var in ("BIOBOX_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "POLARS_MAX_THREADS"):
        assert env[var] == "6", var


def test_staging_is_empty_when_the_step_starts(prof, root, tmp_path):
    probe = tmp_path / "listing.txt"
    code = ("import os, pathlib, sys; s = pathlib.Path(os.environ['BIOBOX_STAGING']); "
            "pathlib.Path(sys.argv[1]).write_text(repr(sorted(os.listdir(s)))); "
            f"(s / 'ledger/stage=S0/step=a').mkdir(parents=True); "
            f"(s / 'ledger/stage=S0/step=a/part-0.parquet').write_bytes(b'x')")
    assert run(prof, root, "a", [sys.executable, "-c", code, str(probe)]) == 0
    assert probe.read_text() == "[]"


def test_stop_before_start_runs_nothing(prof, root, tmp_path):
    (root / "STOP").write_text("")
    marker = tmp_path / "marker"
    assert run(prof, root, "occ", fake("--marker", str(marker), "--out", f"{OCC}/part-0.parquet=10")) == 75
    assert not marker.exists()
    assert root_files(root) == ["STOP"]


def test_stop_during_the_step_blocks_publishing(prof, root):
    out = f"{OCC}/part-0.parquet"
    assert run(prof, root, "occ", fake("--out", f"{out}=10", "--stop")) == 75
    assert root_files(root) == ["STOP"]
    assert (staging(prof, "occ") / out).exists()  # left for the next attempt to discard


def test_failing_command_leaves_the_root_untouched(prof, root):
    out = f"{OCC}/part-0.parquet"
    assert run(prof, root, "occ", fake("--out", f"{out}=10", "--exit", "3")) == 3
    assert root_files(root) == []
    assert (staging(prof, "occ") / out).exists()


def test_killed_command_exits_1_and_leaves_the_root_untouched(prof, root):
    assert run(prof, root, "occ", fake("--out", f"{OCC}/part-0.parquet=10", "--kill")) == 1
    assert root_files(root) == []


def test_empty_staging_refused(prof, root, capsys):
    assert run(prof, root, "occ", fake()) != 0
    assert "empty" in capsys.readouterr().err
    assert root_files(root) == []


def test_existing_final_path_refused(prof, root, capsys):
    out = f"{OCC}/part-0.parquet"
    (root / out).parent.mkdir(parents=True)
    (root / out).write_bytes(b"earlier stage file")
    assert run(prof, root, "occ", fake("--out", f"{out}=2000")) != 0
    assert out in capsys.readouterr().err
    assert (root / out).read_bytes() == b"earlier stage file"
    assert not record_path(root, "occ").exists()


@pytest.mark.covers("LY-04")
def test_small_file_in_a_multi_file_table_refused(prof, root, capsys):
    small, big = f"{OCC}/part-0.parquet", f"{OCC}/part-1.parquet"
    assert run(prof, root, "occ", fake("--out", f"{small}=10", "--out", f"{big}=2000")) != 0
    err = capsys.readouterr().err
    assert small in err and big not in err
    assert root_files(root) == []


@pytest.mark.covers("LY-04")
def test_small_file_refused_against_a_file_already_in_the_root(prof, root, capsys):
    assert run(prof, root, "a", fake("--out", f"{OCC}/part-0.parquet=2000")) == 0
    assert run(prof, root, "b", fake("--out", f"{OCC}/part-1.parquet=10")) != 0
    assert f"{OCC}/part-1.parquet" in capsys.readouterr().err
    assert not (root / f"{OCC}/part-1.parquet").exists()


@pytest.mark.covers("LY-04")
def test_single_small_file_accepted(prof, root):
    assert run(prof, root, "occ", fake("--out", f"{OCC}/part-0.parquet=10")) == 0
    assert (root / f"{OCC}/part-0.parquet").stat().st_size == 10


@pytest.mark.covers("LY-04")
def test_small_ok_table_accepted(prof, root):
    outs = ["ledger/stage=S0/step=a/part-0.parquet", "ledger/stage=S0/step=b/part-0.parquet"]
    assert run(prof, root, "ledgers", fake("--out", f"{outs[0]}=10", "--out", f"{outs[1]}=20")) == 0
    assert root_files(root) == outs + ["runs/stage=S0/ledgers.json"]


@pytest.mark.covers("LY-04")
def test_shard_companions_are_exempt(prof, root):
    base = "l0seq/seq_shards/stage=S0/shard-0000.fa.gz"
    assert run(prof, root, "shards", fake("--out", f"{base}=10", "--out", f"{base}.fai=10",
                                          "--out", f"{base}.gzi=10")) == 0
    assert (root / f"{base}.gzi").exists()


@pytest.mark.covers("LY-10")
@pytest.mark.parametrize("path", ["notes/readme.txt", "l1/protein_occurrence/stage=S0/part-0.parquet"])
def test_unclaimed_path_refused(prof, root, capsys, path):
    assert run(prof, root, "occ", fake("--out", f"{path}=2000")) != 0
    assert path in capsys.readouterr().err
    assert root_files(root) == []


def test_resume_skips_a_finished_step(prof, root, tmp_path, capsys):
    assert run(prof, root, "occ", fake("--out", f"{OCC}/part-0.parquet=10")) == 0
    before = record_path(root, "occ").read_bytes()
    marker = tmp_path / "marker"
    assert run(prof, root, "occ", fake("--marker", str(marker), "--exit", "9")) == 0
    assert "skipped" in capsys.readouterr().err  # migrated by card C-13: run-step's own lines go to stderr
    assert not marker.exists()
    assert record_path(root, "occ").read_bytes() == before


@pytest.mark.parametrize("damage", ["missing", "size"])
def test_resume_with_a_damaged_output_errors(prof, root, tmp_path, capsys, damage):
    out = f"{OCC}/part-0.parquet"
    assert run(prof, root, "occ", fake("--out", f"{out}=10")) == 0
    if damage == "missing":
        (root / out).unlink()
    else:
        (root / out).write_bytes(b"short")
    marker = tmp_path / "marker"
    assert run(prof, root, "occ", fake("--marker", str(marker), "--out", f"{out}=10")) != 0
    assert out in capsys.readouterr().err
    assert not marker.exists()


@pytest.mark.parametrize("now,load,expected", [
    (NOON_MONDAY, 0.0, 4),
    (datetime(2026, 11, 9, 7, 59), 0.0, 6),
    (datetime(2026, 11, 9, 8, 0), 0.0, 4),
    (datetime(2026, 11, 9, 19, 59), 0.0, 4),
    (datetime(2026, 11, 9, 20, 0), 0.0, 6),
    (NOON_SATURDAY, 0.0, 6),
    (NIGHT_MONDAY, runstep.BUSY_LOAD + 0.5, 4),
    (NOON_SATURDAY, runstep.BUSY_LOAD + 0.5, 4),
    (NIGHT_MONDAY, runstep.BUSY_LOAD, 6),
    (NOON_MONDAY, runstep.BUSY_LOAD + 0.5, 4),
])
def test_thread_choice_by_hour_weekday_and_load(prof, now, load, expected):
    assert runstep.choose_threads(prof, now, load) == expected


def test_thread_choice_reaches_the_step_at_a_busy_night(prof, root, tmp_path):
    env_file = tmp_path / "env.json"
    assert run(prof, root, "occ", fake("--env", str(env_file), "--out", f"{OCC}/part-0.parquet=10"),
               clock=lambda: NIGHT_MONDAY, load=lambda: 50.0) == 0
    assert json.loads(env_file.read_text())["OMP_NUM_THREADS"] == "4"
    assert json.loads(record_path(root, "occ").read_text())["threads"] == 4


def test_requested_threads_within_the_cap(prof, root):
    assert run(prof, root, "occ", fake("--out", f"{OCC}/part-0.parquet=10"), threads=2) == 0
    assert json.loads(record_path(root, "occ").read_text())["threads"] == 2


@pytest.mark.parametrize("threads", [5, 0])
def test_requested_threads_outside_the_cap_refused(prof, root, tmp_path, capsys, threads):
    marker = tmp_path / "marker"
    assert run(prof, root, "occ", fake("--marker", str(marker)), threads=threads) != 0
    assert "threads" in capsys.readouterr().err
    assert not marker.exists()


def test_argv_shapes():
    cmd = ["mmseqs", "easy-search", "--threads", "4"]
    found = lambda name: f"/usr/bin/{name}"  # noqa: E731
    assert runstep.build_argv(cmd, 19, "200G", system="Linux", which=found) == [
        "systemd-run", "--user", "--scope", "--quiet", "-p", "MemoryMax=200G", "-p", "MemorySwapMax=0",
        "nice", "-n", "19", *cmd]
    assert runstep.build_argv(cmd, 15, "200G", system="Darwin", which=found) == ["nice", "-n", "15", *cmd]
    no_systemd = lambda name: None if name == "systemd-run" else f"/usr/bin/{name}"  # noqa: E731
    assert runstep.build_argv(cmd, 19, "200G", system="Linux", which=no_systemd) == ["nice", "-n", "19", *cmd]
    assert runstep.build_argv(cmd, 19, "200G", system="Linux", which=lambda name: None) == cmd


def test_fanout_nice(prof, root):
    assert run(prof, root, "occ", fake("--out", f"{OCC}/part-0.parquet=10"), nice="fanout") == 0
    assert json.loads(record_path(root, "occ").read_text())["nice"] == 15


def test_copy_without_rsync(prof, root):
    no_rsync = lambda name: None if name == "rsync" else shutil.which(name)  # noqa: E731
    out = f"{OCC}/part-0.parquet"
    assert run(prof, root, "occ", fake("--out", f"{out}=10"), which=no_rsync) == 0
    assert (root / out).stat().st_size == 10


def test_stale_incoming_and_staging_cleaned(prof, root):
    stale_incoming = root / ".incoming" / "S0.occ"
    stale_incoming.mkdir(parents=True)
    (stale_incoming / "junk").write_text("from a crash")
    staging(prof, "occ").mkdir(parents=True)
    (staging(prof, "occ") / "junk").write_text("from a crash")
    out = f"{OCC}/part-0.parquet"
    assert run(prof, root, "occ", fake("--out", f"{out}=10")) == 0
    assert not stale_incoming.exists()
    assert root_files(root) == [out, "runs/stage=S0/occ.json"]


TIMEOUT = 60  # seconds one step waits for the other before the test fails


@pytest.mark.covers("D-03")
@pytest.mark.parametrize("copier", ["copytree", "rsync"])
@pytest.mark.parametrize("first", ["a", "b"])
def test_two_steps_publishing_at_once_both_succeed(prof, root, copier, first):
    """Card C-12: steps `a` and `b` publish into one root at once, with nothing outside run-step to keep them apart.
    Ordered without timing: the second step creates its incoming area and then waits, before its copy, until the
    first has looked at the shared `<root>/.incoming` to remove it (run-step before card C-12 listed it and removed
    it when it found it empty) or, if the first never looks, until the first has ended; the first, once it has
    looked, waits for the second's copy before it goes on. Both orders of `a` and `b`, on the `copytree` path and,
    where `rsync` exists, on the `rsync` path. Both must publish, and nothing of the publishing area is left."""
    if copier == "rsync" and shutil.which("rsync") is None:
        pytest.skip("rsync is not on PATH")
    which = shutil.which if copier == "rsync" else (lambda name: None if name == "rsync" else shutil.which(name))
    second = "b" if first == "a" else "a"
    incoming = root.resolve() / ".incoming"
    who = threading.local()
    go, copied = threading.Event(), threading.Event()  # go: the first looked at `.incoming`, or ended
    real_iterdir, real_run, real_copytree = Path.iterdir, subprocess.run, shutil.copytree

    def iterdir(self):
        if self == incoming and getattr(who, "step", None) == first and not go.is_set():
            names = list(real_iterdir(self))
            go.set()
            assert copied.wait(TIMEOUT), "the second step never copied"
            return iter(names)
        return real_iterdir(self)

    def copying(copy, *a, **k):
        if getattr(who, "step", None) == second:
            assert go.wait(TIMEOUT), "the first step never ended"
            try:
                return copy(*a, **k)
            finally:
                copied.set()
        return copy(*a, **k)

    def run_cmd(argv, *a, **k):
        return copying(real_run, argv, *a, **k) if argv and str(argv[0]).endswith("rsync") else real_run(argv, *a, **k)

    results: dict[str, object] = {}

    def publish(step: str) -> None:
        who.step = step
        try:
            results[step] = run(prof, root, step, fake("--out", f"ledger/stage=S0/step={step}/part-0.parquet=10"),
                                which=which)
        except BaseException as e:  # noqa: BLE001 - reported below
            results[step] = e
        finally:
            who.step = None
            (go if step == first else copied).set()  # whatever happened, the other waits no longer

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(Path, "iterdir", iterdir)
        mp.setattr(runstep.subprocess, "run", run_cmd)
        mp.setattr(runstep.shutil, "copytree", lambda *a, **k: copying(real_copytree, *a, **k))
        threads = [threading.Thread(target=publish, args=(s,)) for s in ("a", "b")]
        for t in threads:
            t.start()
        for t in threads:
            t.join(2 * TIMEOUT)
    assert results == {"a": 0, "b": 0}
    assert root_files(root) == ["ledger/stage=S0/step=a/part-0.parquet", "ledger/stage=S0/step=b/part-0.parquet",
                                "runs/stage=S0/a.json", "runs/stage=S0/b.json"]
    assert not incoming.exists()
    assert not staging(prof, "a").exists() and not staging(prof, "b").exists()


def test_build_manifest_ignores_runs(prof, root):
    assert "runs" in catalog.EXCLUDED_TREES
    assert run(prof, root, "ledger", fake("--ledger", "a")) == 0
    (root / "runs" / "stage=S0" / "stray.parquet").write_bytes(b"not a table")
    man = catalog.build_manifest(root, "rS0.1", ["S0"])
    assert [f["path"] for f in man["files"]] == ["ledger/stage=S0/step=a/part-0.parquet"]
    assert man["artifacts"] == []


def test_cli_run_step(prof, root, tmp_path):
    p = tmp_path / "profile.yaml"
    out = f"{OCC}/part-0.parquet"
    code = cli.main(["run-step", "--profile", str(p), "--root", str(root), "--stage", "S0", "--step", "occ",
                     "--nice", "fanout", "--threads", "1", "--", *fake("--out", f"{out}=10")])
    assert code == 0
    rec = json.loads(record_path(root, "occ").read_text())
    assert rec["threads"] == 1 and rec["nice"] == 15
    assert rec["outputs"][0]["path"] == out


STAMPED = re.compile(r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d[+-]\d\d:\d\d (run-step|build-stage): ")
UNIT = r"""
import sys
from pathlib import Path
from biobox import profile, runstep, stage
prof, root, cmd = profile.load(sys.argv[1]), Path(sys.argv[2]), sys.argv[3:]
for step in ("a", "b", "a"):  # the third call finds a's run record: skipped
    runstep.run_step(prof, root, "S0", step, [*cmd, step], threads=1)
stage._say("done")
"""
SAYS = r"""
import os, pathlib, sys
print(f"command {sys.argv[1]} ran")
p = pathlib.Path(os.environ["BIOBOX_STAGING"]) / f"ledger/stage=S0/step={sys.argv[1]}/part-0.parquet"
p.parent.mkdir(parents=True)
p.write_bytes(b"x" * 10)
"""


@pytest.mark.covers("C-13")
def test_c13_run_step_and_build_stage_lines_are_timestamped_and_in_time_order_in_one_log(prof, root, tmp_path):
    """Card C-13 (2): a unit log is stdout and stderr of `build-stage` in one file, with `run-step` running in its
    process. `run-step`'s own lines (published, skipped) and `build-stage`'s are written line-buffered to stderr with
    a timestamp, so they fall between the step commands' output in the order things happened and the last line shows
    progress. Measured in C-10: block-buffered on stdout, every `published` line came at the end of the log."""
    log = tmp_path / "unit.log"
    argv = [sys.executable, "-c", UNIT, str(tmp_path / "profile.yaml"), str(root), sys.executable, "-c", SAYS]
    with log.open("wb") as fh:
        assert subprocess.run(argv, stdout=fh, stderr=subprocess.STDOUT, timeout=120).returncode == 0
    lines = log.read_text().splitlines()
    assert [STAMPED.sub(r"\1: ", ln) for ln in lines] == [
        "command a ran",
        "run-step: S0.a published 1 files",
        "command b ran",
        "run-step: S0.b published 1 files",
        "run-step: S0.a skipped: run record and outputs present",
        "build-stage: done",
    ], lines
    assert all(STAMPED.match(ln) for ln in lines if not ln.startswith("command "))


# --- card C-15: a run record's peak memory is the step's, not its launcher's ----------------------------------------
#
# The B-30 probe on the 1,000x stage found every run record at 344.0 MB, the oracle's included (measured, wilkcomp02,
# 2026-11-10). A launcher HOLDER keeps PARENT_MIB of written bytes resident and runs two steps through run-step in its
# own process, as `bio-box run-step` runs inside its CLI process: `big` writes BIG_MIB, `small` writes nothing. Each
# record must show that step's own peak: the two differ by about BIG_MIB, and both are below the launcher's.

MIB = 2**20
PARENT_MIB, BIG_MIB = 512, 192

HOLDER = r"""
import sys
from pathlib import Path
from biobox import profile, runstep
held = b"\x01" * (int(sys.argv[3]) << 20)  # written, so resident when run-step starts a step
prof, root, cmd = profile.load(sys.argv[1]), Path(sys.argv[2]), sys.argv[4:]
codes = [runstep.run_step(prof, root, "S0", step, [*cmd, step, mib], threads=1)
         for step, mib in (("big", "%MIB%"), ("small", "0"))]
assert len(held) == int(sys.argv[3]) << 20
sys.exit(max(codes))
""".replace("%MIB%", str(BIG_MIB))

ALLOCATES = r"""
import os, pathlib, sys
block = b"\x02" * (int(sys.argv[2]) << 20)  # written, so resident
p = pathlib.Path(os.environ["BIOBOX_STAGING"]) / f"ledger/stage=S0/step={sys.argv[1]}/part-0.parquet"
p.parent.mkdir(parents=True)
p.write_bytes(b"x" * 10)
"""


@pytest.mark.covers("C-15")
def test_c15_peak_memory_is_the_steps_own_not_the_launchers(prof, root, tmp_path):
    argv = [sys.executable, "-c", HOLDER, str(tmp_path / "profile.yaml"), str(root), str(PARENT_MIB),
            sys.executable, "-c", ALLOCATES]
    res = subprocess.run(argv, capture_output=True, text=True, timeout=300)
    assert res.returncode == 0, res.stderr
    big, small = (json.loads(record_path(root, s).read_text()) for s in ("big", "small"))
    peaks = {"big": big["peak_rss_bytes"] / MIB, "small": small["peak_rss_bytes"] / MIB}
    print(f"recorded peak MiB: {peaks}")
    assert small["peak_rss_bytes"] < PARENT_MIB * MIB, peaks  # the launcher's image is not the step's
    assert big["peak_rss_bytes"] < PARENT_MIB * MIB, peaks
    assert big["peak_rss_bytes"] >= BIG_MIB * MIB, peaks  # what the step itself wrote is counted
    assert big["peak_rss_bytes"] - small["peak_rss_bytes"] >= (BIG_MIB - 32) * MIB, peaks
    # The source follows the cap the record names: a host with systemd-run and no inherited cap (GitHub Actions)
    # opens a scope for the step and reads its memory.peak; elsewhere (the pod's inherited cap) it is the helper's.
    expected = "scope_memory_peak" if big["mem_cap"]["how"] == "opened" else "child_maxrss"
    assert big["peak_rss_source"] == small["peak_rss_source"] == expected, (big["mem_cap"], big["peak_rss_source"])


# A fake /proc/self/cgroup and cgroup v2 tree, and a stand-in `systemd-run` that drops its own options and execs the
# rest (as the real one does with --scope): run-step opens a scope, and the step's measuring helper, inside it, reads
# the scope's memory.peak.
SCOPE = "user.slice/user-1000.slice/user@1000.service/app.slice/run-c15.scope"


def _cgroup_tree(tmp_path: Path, memory_max: str, memory_peak: int) -> tuple[Path, Path]:
    proc, cg = tmp_path / "proc", tmp_path / "cg"
    (proc / "self").mkdir(parents=True)
    (proc / "self" / "cgroup").write_text(f"0::/{SCOPE}\n")
    (cg / SCOPE).mkdir(parents=True)
    (cg / SCOPE / "memory.max").write_text(memory_max + "\n")
    (cg / SCOPE / "memory.peak").write_text(f"{memory_peak}\n")
    return proc, cg


@pytest.mark.covers("C-15")
def test_c15_an_opened_scope_records_its_memory_peak(prof, root, tmp_path, monkeypatch):
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "systemd-run").write_text('#!/bin/sh\nwhile [ "$1" != nice ]; do shift; done\nexec "$@"\n')
    (bindir / "systemd-run").chmod(0o755)
    monkeypatch.setenv("PATH", f"{bindir}:{os.environ.get('PATH', '')}")
    peak = 3 * 2**30 + 12345
    proc, cg = _cgroup_tree(tmp_path, "max", peak)
    which = lambda n: str(bindir / n) if n == "systemd-run" else shutil.which(n)  # noqa: E731
    assert run(prof, root, "occ", fake("--out", f"{OCC}/part-0.parquet=10"), system="Linux", which=which,
               proc_root=proc, cgroup_root=cg) == 0
    rec = json.loads(record_path(root, "occ").read_text())
    assert rec["mem_cap"]["how"] == "opened"
    assert rec["peak_rss_bytes"] == peak and rec["peak_rss_source"] == "scope_memory_peak"


@pytest.mark.covers("C-15")
def test_c15_an_inherited_cgroup_peak_is_not_the_steps(prof, root, tmp_path):
    """Under `scripts/pod_run.sh` the cgroup with the cap is the whole unit: build-stage, every run-step and every
    shard running at once. Its memory.peak is theirs together over the unit's life, so the record does not take it."""
    unit_peak = 50 * 2**30
    proc, cg = _cgroup_tree(tmp_path, str(32 * 2**30), unit_peak)
    no_systemd = lambda n: None if n == "systemd-run" else shutil.which(n)  # noqa: E731
    assert run(prof, root, "occ", fake("--out", f"{OCC}/part-0.parquet=10"), system="Linux", which=no_systemd,
               proc_root=proc, cgroup_root=cg) == 0
    rec = json.loads(record_path(root, "occ").read_text())
    assert rec["mem_cap"]["how"] == "inherited"
    assert rec["peak_rss_source"] == "child_maxrss"
    assert 1 << 20 < rec["peak_rss_bytes"] < 1 << 30


def test_cli_stop_exit_code(prof, root, tmp_path):
    (root / "STOP").write_text("")
    code = cli.main(["run-step", "--profile", str(tmp_path / "profile.yaml"), "--root", str(root),
                     "--stage", "S0", "--step", "occ", "--", *fake()])
    assert code == 75


# --- C-16: a killed step has a peak ---------------------------------------------------------------------------------

HOLD_MIB = 256
# The end of a step that is to be killed: it says it is ready in the marker named by its last argument (its pid and
# process group, written under a temporary name and renamed), then waits to be killed.
MARK = r"""
with open(sys.argv[-1] + ".tmp", "w") as f:
    f.write(f"{os.getpid()} {os.getpgid(0)}")
os.replace(sys.argv[-1] + ".tmp", sys.argv[-1])
while True:
    signal.pause()
"""
# writes argv[1] MiB, so that they are resident, and holds them
HOLD = "import os, signal, sys\nblock = b'\\x03' * (int(sys.argv[1]) << 20)\n" + MARK
# plays the kernel's part in a stand-in cgroup tree: sets the scope's memory.peak (argv[1]) to argv[2]
SCOPED = "import os, signal, sys\nopen(sys.argv[1], 'w').write(sys.argv[2] + '\\n')\n" + MARK


def _wait_for(pred, timeout: float) -> bool:
    """Poll `pred` every 0.05 s until it holds or `timeout` seconds pass (AGENTS §5.3)."""
    import time

    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if pred():
            return True
        time.sleep(0.05)
    return False


def failed_path(prof, step) -> Path:
    return Path(prof.tmp) / "staging" / f"S0.{step}.failed.json"


def _sampled(prof, step) -> dict:
    """What the measuring helper has written so far for a running step (its file beside staging), or {}."""
    try:
        return json.loads((Path(prof.tmp) / "staging" / f"S0.{step}.peak.json").read_text())
    except (OSError, ValueError):
        return {}


def _hold_and_kill(prof, root, step, cmd, kill, sampled, **kw) -> list[int]:
    """Run `cmd` as step `step` in a thread, with a marker path as its last argument; once the marker exists, kill
    the step's own process (`kill="step"`) or its whole process group, the measuring helper included
    (`kill="group"`), after the helper's file satisfies `sampled`. Whatever the outcome, the step's process group is
    gone when this returns. Returns run-step's exit codes."""
    import signal

    marker = Path(prof.tmp) / f"{step}.marker"
    codes: list[int] = []
    t = threading.Thread(target=lambda: codes.append(run(prof, root, step, [*cmd, str(marker)], threads=1, **kw)))
    t.start()
    pgid = None
    try:
        assert _wait_for(marker.exists, 60), "the step never said it was ready"
        pid, pgid = map(int, marker.read_text().split())
        if kill == "group":
            # Not asserted: before card C-16 the helper writes nothing until the command ends; the test goes on to
            # show what the record then lacks.
            _wait_for(lambda: sampled(_sampled(prof, step)), 30)
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
    assert not t.is_alive()
    return codes


class C16PeakEvidence(UserWarning):
    """Not a fault: a peak a C-16 test saw, put in the warnings summary so that the plain `pytest -q` of GitHub Actions
    (where the scope-only tests run) shows it in the log."""


def _evidence(record_property, name: str, text: str) -> None:
    import warnings

    record_property(name, text)
    print(text)
    warnings.warn(C16PeakEvidence(f"{name}: {text}"), stacklevel=2)


@pytest.mark.covers("C-16")
@pytest.mark.parametrize("kill", ["step", "group"])
def test_c16_a_killed_step_records_its_peak_and_where_it_came_from(prof, root, kill, record_property):
    """A step killed alone (as the kernel's OOM killer takes the largest process of a cgroup: the helper survives and
    measures it) or with its whole process group (a scope stopped after an OOM kill, `systemctl kill`, SIGKILL of the
    group: the helper dies too). Before card C-16 a failed step left no record and its helper's file was deleted,
    and a killed helper wrote nothing. Constructed truth: the step wrote HOLD_MIB MiB before it was killed.

    Without a scope (the pod) the record is the step's own resident high-water mark: the HOLD_MIB MiB and the
    interpreter's pages, within 20% above. In a scope (GitHub Actions) it is the scope's memory.peak, which counts
    what is charged to the scope, not a process's resident set: the HOLD_MIB MiB are anonymous pages the step wrote
    inside the scope (a fresh image from `exec`, MemorySwapMax=0), so they are charged there and bound it from below;
    what else is charged (the helper, page cache the scope read first) is not a share of any process's RSS, so the
    only upper bound is the scope's MemoryMax, which a charge never passes."""
    held = lambda s: max(s.get("maxrss_bytes") or 0, s.get("memory_peak_bytes") or 0) >= HOLD_MIB * MIB  # noqa: E731
    codes = _hold_and_kill(prof, root, "hold", [sys.executable, "-c", HOLD, str(HOLD_MIB)], kill, held)
    assert codes == [1]
    assert root_files(root) == []  # nothing published, as before
    assert failed_path(prof, "hold").exists(), "a killed step left no record of its attempt"
    att = json.loads(failed_path(prof, "hold").read_text())
    _evidence(record_property, f"c16_killed_step_peak_{kill}",
              f"recorded {att['peak_rss_bytes'] and round(att['peak_rss_bytes'] / MIB, 1)} MiB "
              f"({att['peak_rss_source']}, mem_cap {att['mem_cap']['how']}); the step wrote {HOLD_MIB} MiB")
    scope = att["mem_cap"]["how"] == "opened"
    expected = {("step", True): "scope_memory_peak", ("step", False): "child_maxrss",
                ("group", True): "scope_memory_peak_sampled", ("group", False): "child_hwm_sampled"}[kill, scope]
    assert att["peak_rss_source"] == expected, att
    upper = runstep.size_bytes(prof.mem_scope) if scope else 1.2 * HOLD_MIB * MIB
    assert HOLD_MIB * MIB <= att["peak_rss_bytes"] <= upper, att["peak_rss_bytes"] / MIB
    assert (att["stage"], att["step"], att["exit_code"]) == ("S0", "hold", -9)


@pytest.mark.covers("C-16")
def test_c16_a_scope_killed_with_its_helper_keeps_the_last_sampled_memory_peak(prof, root, tmp_path, monkeypatch):
    """An opened scope (the C-15 stand-ins: a `systemd-run` that execs the command, a cgroup tree whose memory.peak the
    step itself sets, as the kernel would) whose every process is killed: the record has the scope's memory.peak as the
    helper last read it, and says so."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "systemd-run").write_text('#!/bin/sh\nwhile [ "$1" != nice ]; do shift; done\nexec "$@"\n')
    (bindir / "systemd-run").chmod(0o755)
    monkeypatch.setenv("PATH", f"{bindir}:{os.environ.get('PATH', '')}")
    proc, cg = _cgroup_tree(tmp_path, "max", 4096)
    peak = 5 * 2**30 + 777
    which = lambda n: str(bindir / n) if n == "systemd-run" else shutil.which(n)  # noqa: E731
    cmd = [sys.executable, "-c", SCOPED, str(cg / SCOPE / "memory.peak"), str(peak)]
    codes = _hold_and_kill(prof, root, "scoped", cmd, "group", lambda s: s.get("memory_peak_bytes") == peak,
                           system="Linux", which=which, proc_root=proc, cgroup_root=cg)
    assert codes == [1]
    assert failed_path(prof, "scoped").exists(), "a killed step left no record of its attempt"
    att = json.loads(failed_path(prof, "scoped").read_text())
    assert att["mem_cap"]["how"] == "opened"
    assert (att["peak_rss_bytes"], att["peak_rss_source"]) == (peak, "scope_memory_peak_sampled")


OOM_CAP = "256M"
GROWS = r"""
blocks = []
for _ in range(64):  # 1 GiB in 16 MiB blocks, each written so that it is resident
    blocks.append(b"\x04" * (16 << 20))
"""


@pytest.mark.covers("C-16")
@pytest.mark.skipif(runstep.memory_cap(OOM_CAP, require=False)["how"] != "opened",
                    reason="needs systemd-run and no cgroup cap above this process, so that run-step opens a scope "
                           "whose MemoryMax the kernel enforces (GitHub Actions; not inside the image on the pod)")
def test_c16_an_oom_killed_step_records_a_peak_within_20_percent_of_its_cap(prof, root, tmp_path, record_property):
    """The card's acceptance with the kernel's own OOM killer: a step that grows past its scope's MemoryMax. The
    record (the scope's memory.peak) is compared with the cap, not with a process's resident set: the kernel OOM-kills
    when a charge cannot be met at memory.max, so the high-water mark is then at the cap and never above it."""
    p = tmp_path / "tight.yaml"
    data = {f: getattr(prof, f) for f in prof.__dataclass_fields__}
    data["mem_scope"] = OOM_CAP
    p.write_text(yaml.safe_dump(data))
    tight = profile.load(p)
    assert run(tight, root, "grows", [sys.executable, "-c", GROWS], threads=1) == 1
    assert failed_path(tight, "grows").exists(), "an OOM-killed step left no record of its attempt"
    att = json.loads(failed_path(tight, "grows").read_text())
    cap = runstep.size_bytes(OOM_CAP)
    _evidence(record_property, "c16_oom_step_peak",
              f"recorded {att['peak_rss_bytes'] and round(att['peak_rss_bytes'] / MIB, 1)} MiB "
              f"({att['peak_rss_source']}); cap {cap / MIB:.1f} MiB")
    assert att["peak_rss_source"] in ("scope_memory_peak", "scope_memory_peak_sampled"), att
    assert 0.8 * cap <= att["peak_rss_bytes"] <= 1.2 * cap, att["peak_rss_bytes"] / MIB


@pytest.mark.covers("C-16")
def test_c16_the_next_attempt_removes_the_record_of_a_failed_one(prof, root):
    assert run(prof, root, "occ", fake("--out", f"{OCC}/part-0.parquet=10", "--exit", "3")) == 3
    assert failed_path(prof, "occ").exists(), "a failed step left no record of its attempt"
    att = json.loads(failed_path(prof, "occ").read_text())
    assert att["exit_code"] == 3 and att["peak_rss_source"] in ("child_maxrss", "scope_memory_peak")
    assert run(prof, root, "occ", fake("--out", f"{OCC}/part-0.parquet=10")) == 0
    assert not failed_path(prof, "occ").exists()
    assert json.loads(record_path(root, "occ").read_text())["exit_code"] == 0
