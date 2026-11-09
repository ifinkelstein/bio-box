"""C-02: `bio-box run-step` with fake step commands (ring R2 by directory).

The fake step is `python -c FAKE ...`. It writes the files named by `--out PATH=BYTES` under
`$BIOBOX_STAGING` and can also dump its environment, touch a marker, create `<root>/STOP`, exit
with a code, or kill itself. Output files are filler bytes, not valid Parquet: the publish rules
under test (claim by path, append-only, small-file floor) look only at paths and sizes.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
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
    assert "skipped" in capsys.readouterr().out
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
        "systemd-run", "--user", "--scope", "--quiet", "-p", "MemoryMax=200G", "nice", "-n", "19", *cmd]
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


def test_build_manifest_ignores_runs(prof, root):
    assert "runs" in catalog.EXCLUDED_TREES
    assert run(prof, root, "ledger", fake("--ledger", "a")) == 0
    (root / "runs" / "stage=S0" / "stray.parquet").write_bytes(b"not a table")
    man = catalog.build_manifest(root, "rS0.1")
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


def test_cli_stop_exit_code(prof, root, tmp_path):
    (root / "STOP").write_text("")
    code = cli.main(["run-step", "--profile", str(tmp_path / "profile.yaml"), "--root", str(root),
                     "--stage", "S0", "--step", "occ", "--", *fake()])
    assert code == 75
