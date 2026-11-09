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
2. Recover (card C-03): an intent file `<root>/runs/stage=<S>/<step>.publishing.json` means an
   earlier attempt died while publishing. Each output it lists is kept if it is in the root with
   the listed sha256, else renamed in from `.incoming` if it is there with the listed sha256.
   If every output is then in place, the run record is written from the intent (no rerun) and
   run-step exits 0. Otherwise the root files that the intent lists and whose sha256 equals the
   intent's are deleted (never any other file), the intent, `.incoming` and staging are removed,
   and the step runs afresh. A root file at a listed path with another sha256 is not this step's:
   exit 2, naming it, and change nothing.
3. Resume: `<root>/runs/stage=<S>/<step>.json` exists and every output it lists is in the root
   with its recorded size: exit 0. A listed output missing or of another size is an error;
   nothing is rebuilt over existing files.
4. Remove a staging or `<root>/.incoming/<stage>.<step>/` directory left by an earlier attempt.
5. Run the command inside `systemd-run --user --scope -p MemoryMax=<mem_scope>` (Linux, when
   present) at the profile's nice, in its own session and process group, and measure wall time,
   CPU-minutes and peak RSS (`os.wait4`). SIGTERM and SIGINT to run-step meanwhile are forwarded
   to that process group; run-step waits for the command, publishes nothing, leaves staging and
   exits 128 + the signal number.
6. Command failed: exit with its code (1 if killed); staging is left for inspection.
   `STOP` appeared meanwhile: exit 75 without publishing.
7. Publish (all or nothing): refuse empty staging, a final path that exists (append-only,
   AGENTS §3.10), a path no registry table or artifact claims (LY-10), and small files by the
   rule of `contracts.check_file_sizes` (LY-04). Then write the intent file atomically (every
   output with bytes and sha256, plus everything the run record holds), one `rsync -a` into
   `<root>/.incoming/<stage>.<step>/`, `os.replace` of each file into place, and removal of the
   incoming and staging directories.
8. Write the run record atomically, then delete the intent file. `runs/` is not part of the
   release manifest (it holds timestamps and host names) and is in `catalog.EXCLUDED_TREES`.
"""

from __future__ import annotations

import json
import os
import platform
import re
import shutil
import signal
import socket
import subprocess
import sys
import threading
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
FORWARDED = (signal.SIGTERM, signal.SIGINT)
_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")

Hook = Callable[[str], None]


def _no_hook(point: str) -> None:
    pass


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


def _record_bytes(record: dict) -> bytes:
    return json.dumps(record, sort_keys=True, indent=1).encode() + b"\n"


def _ours(path: Path, sha256: str) -> bool:
    """`path` is a regular file whose content has this sha256."""
    return path.is_file() and not path.is_symlink() and ids.sha256_file(path).hex() == sha256


def _remove_leftovers(incoming: Path, stage_dir: Path) -> None:
    for d in (incoming, stage_dir):
        if d.exists():
            shutil.rmtree(d)
    if incoming.parent.is_dir() and not any(incoming.parent.iterdir()):
        incoming.parent.rmdir()


def _recover(root: Path, intent_path: Path, record_path: Path, incoming: Path, stage_dir: Path) -> int | None:
    """Finish or undo a publish that an earlier attempt died in (module docstring, step 2).
    Return the exit code when there is nothing left to run, or None to run the step afresh."""
    intent = json.loads(intent_path.read_text())
    keep: list[Path] = []
    move: list[tuple[Path, Path]] = []
    foreign: list[str] = []
    lost: list[str] = []
    for o in intent["outputs"]:
        dst = root / o["path"]
        if (dst.exists() or dst.is_symlink()) and _ours(dst, o["sha256"]):
            keep.append(dst)
        elif dst.exists() or dst.is_symlink():
            foreign.append(o["path"])
        elif _ours(incoming / o["path"], o["sha256"]):
            move.append((incoming / o["path"], dst))
        else:
            lost.append(o["path"])
    if foreign:
        _err(f"{intent_path} lists {', '.join(foreign)}, but the root holds other content there; "
             f"not this step's file, left untouched; nothing done")
        return EXIT_REFUSED
    if not lost:
        for src, dst in move:
            dst.parent.mkdir(parents=True, exist_ok=True)
            os.replace(src, dst)
        _atomic_write(record_path, _record_bytes(intent))
        _remove_leftovers(incoming, stage_dir)
        intent_path.unlink()
        print(f"run-step: {intent['stage']}.{intent['step']} recovered an interrupted publish of "
              f"{len(intent['outputs'])} files")
        return 0
    for p in keep:
        p.unlink()
    intent_path.unlink()
    _remove_leftovers(incoming, stage_dir)
    _err(f"interrupted publish of {intent['stage']}.{intent['step']} cannot be finished "
         f"({', '.join(lost)} lost); removed its {len(keep)} published files; running it again")
    return None


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


def _publish(root: Path, stage_dir: Path, incoming: Path, intent_path: Path, meta: dict, floor: int,
             which: Callable[[str], str | None], hook: Hook) -> dict:
    """Check, then publish the staged outputs; return the run record (`meta` plus `outputs`).
    `hook` is called at each point a crash can leave work half done (tests only)."""
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
    records = sorted(({"path": e["path"], "bytes": e["bytes"],
                       "sha256": ids.sha256_file(stage_dir / e["path"]).hex()} for e, _ in outputs),
                     key=lambda o: o["path"])
    record = {**meta, "outputs": records}
    _atomic_write(intent_path, _record_bytes(record))
    hook("after_intent")

    incoming.parent.mkdir(parents=True, exist_ok=True)
    rsync = which("rsync")
    if rsync:
        subprocess.run([rsync, "-a", f"{stage_dir}/", f"{incoming}/"], check=True)
    else:
        shutil.copytree(stage_dir, incoming)
    hook("after_rsync")
    for i, r in enumerate(records):
        src, dst = incoming / r["path"], root / r["path"]
        if src.stat().st_size != r["bytes"]:
            raise RuntimeError(f"{r['path']}: copy into {incoming} has the wrong size")
        if dst.exists():  # checked above; only a concurrent writer gets here
            raise RuntimeError(f"{r['path']}: appeared in the root during publishing")
        dst.parent.mkdir(parents=True, exist_ok=True)
        os.replace(src, dst)
        hook(f"after_rename:{i}")
    _remove_leftovers(incoming, stage_dir)
    return record


def _signal_group(pgid: int, signum: int) -> None:
    try:
        os.killpg(pgid, signum)
    except (ProcessLookupError, PermissionError):  # the group has already exited
        pass


def _execute(argv: list[str], env: dict[str, str]) -> tuple[int, float, float, int, int | None]:
    """Run argv in its own session (so its process group can be signalled as one; a systemd-run
    scope keeps that group); return (exit status, wall seconds, CPU minutes, peak RSS bytes, the
    signal forwarded or None). SIGTERM and SIGINT received meanwhile are forwarded to the group,
    and the command is still waited for; handlers are installed only in the main thread and the
    previous ones restored after."""
    received: list[int] = []
    group: list[int] = []

    def forward(signum, frame) -> None:
        received.append(signum)
        if group:
            _signal_group(group[0], signum)

    previous = {}
    if threading.current_thread() is threading.main_thread():
        previous = {s: signal.signal(s, forward) for s in FORWARDED}
    try:
        t0 = time.monotonic()
        proc = subprocess.Popen(argv, env=env, start_new_session=True)
        group.append(proc.pid)
        if received:  # arrived before the group existed
            _signal_group(proc.pid, received[-1])
        _, status, ru = os.wait4(proc.pid, 0)
        wall = time.monotonic() - t0
    finally:
        for s, h in previous.items():
            signal.signal(s, h)
    proc.returncode = os.waitstatus_to_exitcode(status)
    rss = ru.ru_maxrss if sys.platform == "darwin" else ru.ru_maxrss * 1024  # bytes on macOS, KiB on Linux
    return proc.returncode, wall, (ru.ru_utime + ru.ru_stime) / 60.0, rss, (received[0] if received else None)


def run_step(prof: Profile, root: str | Path, stage: str, step: str, cmd: Sequence[str], *,
             nice: str = "monolithic", threads: int | None = None,
             clock: Callable[[], datetime] = datetime.now, load: Callable[[], float] = _load1,
             which: Callable[[str], str | None] = shutil.which, system: str | None = None,
             _hook: Hook | None = None) -> int:
    """Run one step; return the exit code for `bio-box run-step` (0, 75 for STOP, else non-zero).
    `_hook` (tests only) is called with each publish point name: after_intent, after_rsync,
    after_rename:<i>, before_record, after_record; an exception it raises is not caught."""
    hook = _hook or _no_hook
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
    intent_path = root / "runs" / f"stage={stage}" / f"{step}.publishing.json"
    stage_dir = Path(prof.tmp) / "staging" / f"{stage}.{step}"
    incoming = root / ".incoming" / f"{stage}.{step}"
    if intent_path.exists():
        code = _recover(root, intent_path, record_path, incoming, stage_dir)
        if code is not None:
            return code
    if record_path.exists():
        return _resume(root, record_path)

    for stale in (incoming, stage_dir):
        if stale.exists():
            shutil.rmtree(stale)
    stage_dir.mkdir(parents=True)

    env = dict(os.environ, BIOBOX_STAGING=str(stage_dir), BIOBOX_ROOT=str(root), BIOBOX_THREADS=str(n_threads))
    env.update({v: str(n_threads) for v in THREAD_VARS})
    argv = build_argv(cmd, nice_value, prof.mem_scope, system=system, which=which)
    started_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    code, wall, cpu_min, rss, sig = _execute(argv, env)

    if sig is not None:
        _err(f"{signal.Signals(sig).name} received; forwarded to {stage}.{step}, which exited {code}; "
             f"nothing published, staging left at {stage_dir}")
        return 128 + sig
    if code != 0:
        _err(f"{stage}.{step} failed with exit {code}; nothing published, staging left at {stage_dir}")
        return code if code > 0 else 1
    if (root / "STOP").exists():
        _err(f"{root / 'STOP'} appeared; not publishing {stage}.{step}, staging left at {stage_dir}")
        return EXIT_STOP
    over = cpu_min > prof.cpu_min_per_process
    meta = {
        "stage": stage, "step": step, "argv": argv, "threads": n_threads, "nice": nice_value,
        "started_at": started_at, "wall_s": round(wall, 3), "cpu_min": round(cpu_min, 6),
        "peak_rss_bytes": rss, "host": socket.gethostname(), "exit_code": code,
        "image_digest": os.environ.get("BIOBOX_IMAGE_DIGEST", "none"), "over_cpu_budget": over,
    }
    try:
        record = _publish(root, stage_dir, incoming, intent_path, meta, prof.small_file_floor_bytes,
                          which, hook)
    except Refused as e:
        _err(f"refused: {e}; nothing published, staging left at {stage_dir}")
        return EXIT_REFUSED

    if over:
        _err(f"warning: {stage}.{step} used {cpu_min:.1f} CPU-minutes, over {prof.cpu_min_per_process}")
    hook("before_record")
    _atomic_write(record_path, _record_bytes(record))
    hook("after_record")
    intent_path.unlink()
    print(f"run-step: {stage}.{step} published {len(record['outputs'])} files")
    return 0
