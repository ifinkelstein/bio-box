"""`bio-box run-step`: run one step command under a host profile's caps (card C-02).

    bio-box run-step --profile profiles/wilkcomp02.yaml --root R --stage S0 --step genecall \
        [--nice monolithic|fanout] [--threads N] -- CMD ARGS...

Contract with the step command. It receives `BIOBOX_STAGING` (an empty directory
`<profile.tmp>/staging/<stage>.<step>/`), `BIOBOX_ROOT`, `BIOBOX_THREADS`, and
`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `POLARS_MAX_THREADS` set to the same
count; a step passes `BIOBOX_THREADS` on to DuckDB `threads` and MMseqs2 `--threads` itself. It
writes its outputs under `BIOBOX_STAGING` at their root-relative paths and never into the root.

Order of work:

1. `<root>/STOP` exists: exit 75 (EX_TEMPFAIL) and run nothing.
2. Resume: `<root>/runs/stage=<S>/<step>.json` exists and every output it lists is in the root
   with its recorded size: exit 0. A listed output missing or of another size is an error;
   nothing is rebuilt over existing files.
3. Remove a staging or `<root>/.incoming/<stage>.<step>/` directory left by an earlier attempt.
4. Run the command inside `systemd-run --user --scope -p MemoryMax=<mem_scope>` (Linux, when
   present) at the profile's nice, and measure wall time, CPU-minutes and peak RSS (`os.wait4`).
5. Command failed: exit with its code (1 if killed); staging is left for inspection.
   `STOP` appeared meanwhile: exit 75 without publishing.
6. Publish (all or nothing): refuse empty staging, a final path that exists (append-only,
   AGENTS §3.10), a path no registry table or artifact claims (LY-10), and small files by the
   rule of `contracts.check_file_sizes` (LY-04). Then one `rsync -a` into
   `<root>/.incoming/<stage>.<step>/`, `os.replace` of each file into place, and removal of the
   incoming and staging directories.
7. Write the run record atomically. `runs/` is not part of the release manifest (it holds
   timestamps and host names) and is in `catalog.EXCLUDED_TREES`.
"""

from __future__ import annotations

import json
import os
import platform
import re
import shutil
import socket
import subprocess
import sys
import time
from collections.abc import Callable, Sequence
from datetime import datetime, timezone
from pathlib import Path

from biobox import contracts, ids
from biobox.catalog import LayoutError, classify, classify_artifact
from biobox.profile import Profile
from biobox.registry.schemas import ARTIFACTS, TABLES

EXIT_STOP = 75  # EX_TEMPFAIL (data_layer_and_tests.md §5, STOP file)
EXIT_REFUSED = 2

# Estimate: a 1-minute load average above this means someone else is computing, so the night
# cap (which the footprint agreement grants only when no other compute runs) falls back to the
# day cap. Not measured on wilkcomp02; to be revisited with the pod's idle load.
BUSY_LOAD = 8.0

THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "POLARS_MAX_THREADS")
_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")


class Refused(Exception):
    """A publish or argument check failed; nothing was published."""


def _load1() -> float:
    return os.getloadavg()[0]


def choose_threads(prof: Profile, now: datetime, load1: float) -> int:
    """The thread cap at local time `now`: the profile's cap for the hour and weekday, but never
    more than the day cap while the 1-minute load average is above `BUSY_LOAD`."""
    cap = prof.threads_at(now.hour, now.weekday())
    if load1 > BUSY_LOAD:
        cap = min(cap, prof.max_threads_day)
    return cap


def build_argv(cmd: Sequence[str], nice: int, mem_scope: str, *, system: str | None = None,
               which: Callable[[str], str | None] = shutil.which) -> list[str]:
    """The argv that runs `cmd`: in a memory-capped systemd scope on Linux when `systemd-run`
    exists, else under `nice` when it exists, else bare."""
    system = platform.system() if system is None else system
    if system == "Linux" and which("systemd-run"):
        return ["systemd-run", "--user", "--scope", "--quiet", "-p", f"MemoryMax={mem_scope}",
                "nice", "-n", str(nice), *cmd]
    if which("nice"):
        return ["nice", "-n", str(nice), *cmd]
    return list(cmd)


def _err(msg: str) -> None:
    print(f"run-step: {msg}", file=sys.stderr)


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def _resume(root: Path, record: Path) -> int:
    rec = json.loads(record.read_text())
    bad = []
    for o in rec["outputs"]:
        p = root / o["path"]
        if not p.is_file():
            bad.append(f"{o['path']} missing")
        elif p.stat().st_size != o["bytes"]:
            bad.append(f"{o['path']} is {p.stat().st_size} bytes, recorded {o['bytes']}")
    if bad:
        _err(f"{record} exists but {'; '.join(bad)}; not rebuilding over existing files")
        return EXIT_REFUSED
    print(f"run-step: {rec['stage']}.{rec['step']} skipped: run record and outputs present")
    return 0


def _claim(rel: str) -> tuple[dict, tuple]:
    """The manifest-style entry of an output and its LY-04 group key, as `check_file_sizes` forms them."""
    try:
        hit = classify(rel)
        if hit is not None:
            spec, parts = hit
            key = (spec.name, parts.get("stage") or parts.get("release") or "")
            return {"path": rel, "table": spec.name, "partitions": parts}, key
        art = classify_artifact(rel)
    except LayoutError as e:
        raise Refused(str(e)) from e
    if art is None:
        raise Refused(f"{rel}: no registered table or artifact claims this path (LY-10)")
    spec, parts, _ = art
    return {"path": rel, "artifact": spec.name, "partitions": parts}, (spec.name, parts.get("stage", ""))


def _root_entries(root: Path, entry: dict) -> list[tuple[dict, tuple]]:
    """Files already in the root under the same table or artifact location as `entry`."""
    if "table" in entry:
        loc = TABLES[entry["table"]].location
        paths = [root / loc] if loc.endswith(".parquet") else sorted((root / loc).rglob("*.parquet"))
    else:
        loc = ARTIFACTS[entry["artifact"]].location
        paths = sorted(p for p in (root / loc).rglob("*") if p.is_file())
    out = []
    for p in paths:
        if not p.is_file():
            continue
        try:
            e, key = _claim(p.relative_to(root).as_posix())
        except Refused:
            continue  # a stray is the manifest builder's to report, not this step's
        e["bytes"] = p.stat().st_size
        out.append((e, key))
    return out


def _check_sizes(root: Path, outputs: list[tuple[dict, tuple]], floor: int) -> None:
    """LY-04 with `contracts.check_file_sizes` semantics, over this step's outputs together with
    the files of the same table and stage (or release) already in the root."""
    groups: dict[tuple, list[dict]] = {}
    mine: set[str] = set()
    for e, key in outputs:
        groups.setdefault(key, []).append(e)
        mine.add(e["path"])
    for key in sorted(groups):
        members = groups[key]
        for e, k in _root_entries(root, members[0]):
            if k == key and e["path"] not in mine:
                members.append(e)
        manifest = {"files": [e for e in members if "table" in e],
                    "artifacts": [e for e in members if "artifact" in e]}
        problems = contracts.check_file_sizes(manifest, floor)
        if problems:
            small = [e["path"] for e in members if e["bytes"] < floor and e["path"] in mine]
            small = small or [e["path"] for e in members if e["bytes"] < floor]
            raise Refused(f"{', '.join(sorted(small))}: {problems[0]} (LY-04)")


def _publish(root: Path, stage_dir: Path, incoming: Path, floor: int,
             which: Callable[[str], str | None]) -> list[dict]:
    files = sorted(p for p in stage_dir.rglob("*") if not p.is_dir())
    if not files:
        raise Refused(f"staging {stage_dir} is empty; nothing to publish")
    outputs = []
    for p in files:
        rel = p.relative_to(stage_dir).as_posix()
        if p.is_symlink() or not p.is_file():
            raise Refused(f"{rel}: not a regular file")
        if (root / rel).exists() or (root / rel).is_symlink():
            raise Refused(f"{rel}: already exists in the root; a step never rewrites a file (AGENTS §3.10)")
        e, key = _claim(rel)
        e["bytes"] = p.stat().st_size
        outputs.append((e, key))
    _check_sizes(root, outputs, floor)
    records = [{"path": e["path"], "bytes": e["bytes"], "sha256": ids.sha256_file(stage_dir / e["path"]).hex()}
               for e, _ in outputs]

    incoming.parent.mkdir(parents=True, exist_ok=True)
    rsync = which("rsync")
    if rsync:
        subprocess.run([rsync, "-a", f"{stage_dir}/", f"{incoming}/"], check=True)
    else:
        shutil.copytree(stage_dir, incoming)
    for r in records:
        src, dst = incoming / r["path"], root / r["path"]
        if src.stat().st_size != r["bytes"]:
            raise RuntimeError(f"{r['path']}: copy into {incoming} has the wrong size")
        if dst.exists():  # checked above; only a concurrent writer gets here
            raise RuntimeError(f"{r['path']}: appeared in the root during publishing")
        dst.parent.mkdir(parents=True, exist_ok=True)
        os.replace(src, dst)
    shutil.rmtree(incoming)
    if not any(incoming.parent.iterdir()):
        incoming.parent.rmdir()
    shutil.rmtree(stage_dir)
    return records


def _execute(argv: list[str], env: dict[str, str]) -> tuple[int, float, float, int]:
    """Run argv; return (exit status, wall seconds, CPU minutes, peak RSS bytes) of the child tree."""
    t0 = time.monotonic()
    proc = subprocess.Popen(argv, env=env)
    _, status, ru = os.wait4(proc.pid, 0)
    wall = time.monotonic() - t0
    proc.returncode = os.waitstatus_to_exitcode(status)
    rss = ru.ru_maxrss if sys.platform == "darwin" else ru.ru_maxrss * 1024  # bytes on macOS, KiB on Linux
    return proc.returncode, wall, (ru.ru_utime + ru.ru_stime) / 60.0, rss


def run_step(prof: Profile, root: str | Path, stage: str, step: str, cmd: Sequence[str], *,
             nice: str = "monolithic", threads: int | None = None,
             clock: Callable[[], datetime] = datetime.now, load: Callable[[], float] = _load1,
             which: Callable[[str], str | None] = shutil.which, system: str | None = None) -> int:
    """Run one step; return the exit code for `bio-box run-step` (0, 75 for STOP, else non-zero)."""
    root = Path(root).resolve()
    for what, name in (("stage", stage), ("step", step)):
        if not _NAME.fullmatch(name):
            _err(f"bad {what} name {name!r}")
            return EXIT_REFUSED
    if not cmd:
        _err("no command given after --")
        return EXIT_REFUSED
    if nice not in ("monolithic", "fanout"):
        _err(f"--nice must be monolithic or fanout, got {nice!r}")
        return EXIT_REFUSED
    nice_value = prof.nice_monolithic if nice == "monolithic" else prof.nice_fanout
    cap = choose_threads(prof, clock(), load())
    if threads is not None and not 1 <= threads <= cap:
        _err(f"--threads {threads} refused: the cap now is {cap}")
        return EXIT_REFUSED
    n_threads = cap if threads is None else threads

    if (root / "STOP").exists():
        _err(f"{root / 'STOP'} exists; not starting {stage}.{step}")
        return EXIT_STOP
    record_path = root / "runs" / f"stage={stage}" / f"{step}.json"
    if record_path.exists():
        return _resume(root, record_path)

    stage_dir = Path(prof.tmp) / "staging" / f"{stage}.{step}"
    incoming = root / ".incoming" / f"{stage}.{step}"
    for stale in (incoming, stage_dir):
        if stale.exists():
            shutil.rmtree(stale)
    stage_dir.mkdir(parents=True)

    env = dict(os.environ, BIOBOX_STAGING=str(stage_dir), BIOBOX_ROOT=str(root), BIOBOX_THREADS=str(n_threads))
    env.update({v: str(n_threads) for v in THREAD_VARS})
    argv = build_argv(cmd, nice_value, prof.mem_scope, system=system, which=which)
    started_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    code, wall, cpu_min, rss = _execute(argv, env)

    if code != 0:
        _err(f"{stage}.{step} failed with exit {code}; nothing published, staging left at {stage_dir}")
        return code if code > 0 else 1
    if (root / "STOP").exists():
        _err(f"{root / 'STOP'} appeared; not publishing {stage}.{step}, staging left at {stage_dir}")
        return EXIT_STOP
    try:
        outputs = _publish(root, stage_dir, incoming, prof.small_file_floor_bytes, which)
    except Refused as e:
        _err(f"refused: {e}; nothing published, staging left at {stage_dir}")
        return EXIT_REFUSED

    over = cpu_min > prof.cpu_min_per_process
    if over:
        _err(f"warning: {stage}.{step} used {cpu_min:.1f} CPU-minutes, over {prof.cpu_min_per_process}")
    record = {
        "stage": stage, "step": step, "argv": argv, "threads": n_threads, "nice": nice_value,
        "started_at": started_at, "wall_s": round(wall, 3), "cpu_min": round(cpu_min, 6),
        "peak_rss_bytes": rss, "host": socket.gethostname(), "exit_code": code,
        "image_digest": os.environ.get("BIOBOX_IMAGE_DIGEST", "none"),
        "over_cpu_budget": over, "outputs": sorted(outputs, key=lambda o: o["path"]),
    }
    _atomic_write(record_path, json.dumps(record, sort_keys=True, indent=1).encode() + b"\n")
    print(f"run-step: {stage}.{step} published {len(outputs)} files")
    return 0
