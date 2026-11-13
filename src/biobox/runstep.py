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
   intent's are deleted (never any other file), the intent, the step's `.incoming/<stage>.<step>/`
   and staging are removed, and the step runs afresh. A root file at a listed path with another sha256 is not this step's:
   exit 2, naming it, and change nothing.
3. Resume: `<root>/runs/stage=<S>/<step>.json` exists and every output it lists is in the root
   with its recorded size: exit 0. A listed output missing or of another size is an error;
   nothing is rebuilt over existing files.
4. Remove a staging or `<root>/.incoming/<stage>.<step>/` directory, and the record of a failed
   attempt (step 6), left by an earlier attempt.
5. Decide the memory cap (card C-07, `memory_cap`). A cgroup v2 `memory.max` already above
   run-step (its own cgroup or an ancestor, read by `doctor.cgroup_memory`) is inherited: no
   scope is opened, unless `systemd-run` is present and `mem_scope` is tighter than that cap.
   Inside the image `systemd-run` is not visible, so the cap is the unit `scripts/pod_run.sh`
   opened on the host around `singularity exec` (which binds /sys/fs/cgroup read-only so it can
   be seen). With neither a cap nor `systemd-run`, a profile with `require_mem_cap` refuses to
   start (exit 2, naming scripts/pod_run.sh); otherwise the command runs under `nice` alone. The
   run record's `mem_cap` says which: `opened`, `inherited` (with the bytes and the cgroup that
   sets them) or `none`.
   Run the command (inside `systemd-run --user --scope -p MemoryMax=<mem_scope> -p
   MemorySwapMax=0` when opened) at the profile's nice, in its own session and process group,
   and measure wall time and CPU-minutes (`os.wait4`). The step's peak memory (card C-15) is not
   run-step's `wait4`: Linux counts the launcher's image, forked into the child before `exec`, in
   the child's `ru_maxrss`. The command runs under a small measuring helper (`_MEASURE`) that
   spawns it and records, as `peak_rss_bytes` with `peak_rss_source`:
   - `scope_memory_peak` when run-step opened a scope for this step: that scope's cgroup
     `memory.peak`, read before the scope goes away. It holds only this step and the helper, and
     counts what the cap counts: every process of the step together, page cache charged to it
     included.
   - `child_maxrss` otherwise, notably under an inherited cap: the `ru_maxrss` of the step's
     process tree as the helper's `wait4` sees it (the largest single process, not the sum; the
     helper's own image as a floor: a step that allocates nothing recorded 10.5 MiB, and one that
     writes 192 MiB recorded 203.0 MiB, from a launcher holding 512 MiB whose steps both recorded
     about 606 MiB before; measured, wilkcomp02, 2026-11-10, `test_c15_peak_memory_is_the_steps_own_not_the_launchers`).
     An inherited cgroup is the whole `scripts/pod_run.sh` unit:
     build-stage, every run-step and every shard running at once, over the unit's life, so its
     `memory.peak` is not one step's and is not used.
   - `scope_memory_peak_sampled` or `child_hwm_sampled` (card C-16) when the helper was killed
     with the command (a stopped scope, SIGKILL of the process group): its last sample, taken
     every `SAMPLE_S` seconds while the command ran, of the scope's `memory.peak` or of the
     largest `VmHWM` (the kernel's per-process high-water mark) in the step's process tree.
   - `unmeasured` (and null) if the helper wrote nothing.
   The command dies when the helper dies (PR_SET_PDEATHSIG, on Linux).
   The record's `argv` is the command under its limits, without the helper. SIGTERM and SIGINT
   to run-step meanwhile are forwarded to that process group; run-step waits for the command,
   publishes nothing, leaves staging and exits 128 + the signal number.
6. Command failed: exit with its code (1 if killed); staging is left for inspection.
   `STOP` appeared meanwhile: exit 75 without publishing.
   A command that failed, was killed or was stopped by a forwarded signal leaves the record of
   its attempt (card C-16), `<profile tmp>/staging/<stage>.<step>.failed.json` (`failed_path`):
   what a run record holds but outputs, `forwarded_signal` beside it, the peak memory and its
   source among them; the log line says the peak too. It is outside the root and is not a run
   record; the next attempt of the step removes it with the earlier staging.
7. Publish (all or nothing): refuse empty staging, a final path that exists (append-only,
   AGENTS §3.10), a path no registry table or artifact claims (LY-10), and small files by the
   rule of `contracts.check_file_sizes` (LY-04). Then write the intent file atomically (every
   output with bytes and sha256, plus everything the run record holds), one `rsync -a` (or
   `copytree` where there is no rsync) into `<root>/.incoming/<stage>.<step>/`, created with its
   parent just before the copy, `os.replace` of each file into place, and removal of the step's own
   incoming and staging directories.
   Directories shared between steps that run at once (card C-12): `<root>/.incoming` is removed only
   when it is empty, by `rmdir` under an exclusive `flock` on the root directory, while a step
   creates its incoming directory under a shared one, so a step publishing beside another never
   loses its parent; it is never listed and then removed. The staging parent `<profile tmp>/staging`,
   `runs/stage=<S>/` and the output directories are created with their parents and never removed.
8. Write the run record atomically, then delete the intent file. `runs/` is not part of the
   release manifest (it holds timestamps and host names) and is in `catalog.EXCLUDED_TREES`.

Every line run-step writes of its own (published, skipped, recovered, refused, failed) goes to stderr
after a timestamp and is flushed at once (`log`, card C-13), so a unit log is a timeline.
"""

from __future__ import annotations

import contextlib
import fcntl
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
from collections.abc import Callable, Iterator, Sequence
from datetime import datetime, timezone
from pathlib import Path

from biobox import codeid, contracts, ids
from biobox.catalog import LayoutError, classify, classify_artifact
from biobox.doctor import cgroup_memory
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

PROC_ROOT = Path("/proc")
CGROUP_ROOT = Path("/sys/fs/cgroup")  # the cgroup v2 mount
LAUNCHER = "scripts/pod_run.sh"  # the host-side launcher that opens the unit around `singularity exec`
_SIZE = re.compile(r"(\d+)([KMGTPE]?)")  # systemd's MemoryMax= byte syntax (suffixes are powers of 1024)

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


class NoCap(Exception):
    """The profile requires a memory cap and there is neither a scope to open nor a cap to inherit."""


def size_bytes(spec: str | None) -> int | None:
    """Bytes of a systemd `MemoryMax=` value (`200G`, `512M`, `4000000000`), or None for anything
    else (None, `infinity`, a percentage)."""
    m = _SIZE.fullmatch(spec.strip()) if spec else None
    if m is None:
        return None
    return int(m.group(1)) * 1024 ** " KMGTPE".index(m.group(2) or " ")


def detect_cap(proc_root: str | Path | None = None, cgroup_root: str | Path | None = None) -> dict | None:
    """The memory cap this process already runs under: the smallest numeric cgroup v2
    `memory.max` on its own cgroup or an ancestor (`doctor.cgroup_memory`, card A-04), as
    {"memory_max_bytes", "limited_at", "cgroup"}; None when there is none (every level `max`,
    cgroup v1, or no tree visible: inside the image /sys/fs/cgroup is empty unless bound in)."""
    return _cap_of(_cgroup(proc_root, cgroup_root))


def _cgroup(proc_root: str | Path | None, cgroup_root: str | Path | None) -> dict:
    return cgroup_memory(PROC_ROOT if proc_root is None else proc_root,
                         CGROUP_ROOT if cgroup_root is None else cgroup_root)


def _cap_of(cg: dict) -> dict | None:
    if not cg["available"] or cg["memory_max_bytes"] is None:
        return None
    return {"memory_max_bytes": cg["memory_max_bytes"], "limited_at": cg["limited_at"], "cgroup": cg["path"]}


def memory_cap(mem_max: str | None, *, require: bool, system: str | None = None,
               which: Callable[[str], str | None] | None = None, proc_root: str | Path | None = None,
               cgroup_root: str | Path | None = None) -> dict:
    """How a command started here is memory-capped (card C-07): {"how": "inherited", ...the cap
    `detect_cap` found}, {"how": "opened", "memory_max": mem_max} (a systemd scope will be opened
    for it), or {"how": "none"}. A cap found is inherited unless `systemd-run` is present and
    `mem_max` is tighter than it. Raises NoCap when `require` and neither is possible."""
    system = platform.system() if system is None else system
    which = shutil.which if which is None else which
    cg = _cgroup(proc_root, cgroup_root) if system == "Linux" else {"available": False, "reason": system}
    found = _cap_of(cg)
    can_open = system == "Linux" and which("systemd-run") is not None
    want = size_bytes(mem_max)
    if found is not None and not (can_open and want is not None and want < found["memory_max_bytes"]):
        return {"how": "inherited", **found}
    if can_open:
        return {"how": "opened", "memory_max": mem_max}
    if require:
        why = (f"no cgroup v2 memory cap is visible ({cg['reason']})" if not cg["available"]
               else f"no memory.max is set on {cg['path']} or its parents")
        raise NoCap(f"the profile requires a memory cap (require_mem_cap), but systemd-run is not on PATH "
                    f"and {why}. Inside the image, start this from the host with {LAUNCHER}, which opens "
                    f"a MemoryMax unit around singularity exec and binds /sys/fs/cgroup so the cap is seen")
    return {"how": "none"}


def scope_argv(cmd: Sequence[str], nice: int, *, mem_max: str | None, cpu_quota: int | None = None,
               inherited: dict | None = None, system: str | None = None,
               which: Callable[[str], str | None] = shutil.which) -> tuple[list[str], str]:
    """The argv that runs `cmd` under a host's limits, and which way it does ("scope",
    "inherited", "nice" or "bare"). With `inherited` (a cap `memory_cap` found above this
    process, card C-07): no scope, under the `nice` binary when it exists. Else on Linux with
    `systemd-run`: a transient user scope (`MemoryMax=<mem_max>` when given, `MemorySwapMax=0`
    so a cap is never stretched by swap, `CPUQuota=<cpu_quota>%` when given) at `nice`;
    systemd-run and nice exec the command, so it keeps the caller's pid, session and rlimits.
    Else under the `nice` binary when it exists, else bare (the caller sets its own niceness
    then). Shared by run-step and the search runner (card C-05)."""
    system = platform.system() if system is None else system
    if inherited is not None:
        return (["nice", "-n", str(nice), *cmd] if which("nice") else list(cmd)), "inherited"
    if system == "Linux" and which("systemd-run"):
        props = [] if mem_max is None else ["-p", f"MemoryMax={mem_max}"]
        props += ["-p", "MemorySwapMax=0"]
        props += [] if cpu_quota is None else ["-p", f"CPUQuota={cpu_quota}%"]
        return ["systemd-run", "--user", "--scope", "--quiet", *props, "nice", "-n", str(nice), *cmd], "scope"
    if which("nice"):
        return ["nice", "-n", str(nice), *cmd], "nice"
    return list(cmd), "bare"


def build_argv(cmd: Sequence[str], nice: int, mem_scope: str, *, inherited: dict | None = None,
               system: str | None = None, which: Callable[[str], str | None] = shutil.which) -> list[str]:
    """The argv that runs a step command: under `nice` inside an inherited cap, else in a
    memory-capped systemd scope on Linux when `systemd-run` exists, else under `nice` when it
    exists, else bare (`scope_argv`)."""
    return scope_argv(cmd, nice, mem_max=mem_scope, inherited=inherited, system=system, which=which)[0]


def log(who: str, msg: str) -> None:
    """One line of run-step's or build-stage's own (card C-13): to stderr, after a timestamp (local time with its UTC
    offset, to the second), flushed at once. A unit log holds stdout and stderr of `build-stage`, in whose process
    `run-step` runs; lines on a block-buffered stdout reached the log only when the process ended (C-10), so the log
    was not a timeline and its last line showed no progress."""
    stamp = datetime.now().astimezone().isoformat(timespec="seconds")
    print(f"{stamp} {who}: {msg}", file=sys.stderr, flush=True)


def _err(msg: str) -> None:
    log("run-step", msg)


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


@contextlib.contextmanager
def _incoming_lock(root: Path, exclusive: bool) -> Iterator[bool]:
    """Hold a `flock` on the data root's directory while `<root>/.incoming` is created into (shared) or removed
    (exclusive) (card C-12); yield whether it is held. Released when the block ends. A filesystem that cannot lock a
    directory (an NFS mount refuses an exclusive lock on a directory opened read-only) gives False."""
    try:
        fd = os.open(root, os.O_RDONLY)
    except OSError:
        yield False
        return
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH)
            held = True
        except OSError:
            held = False
        yield held
    finally:
        os.close(fd)


def _make_incoming(incoming: Path) -> None:
    """Create this step's `<root>/.incoming/<stage>.<step>/`, with its parent, just before copying into it. Under the
    shared lock, so another step's `_remove_leftovers` cannot remove the parent between its creation and the
    step's own directory's; once that exists the parent is not empty and is never removed."""
    with _incoming_lock(incoming.parent.parent, exclusive=False):
        incoming.mkdir(parents=True, exist_ok=True)


def _remove_leftovers(incoming: Path, stage_dir: Path) -> None:
    """Remove this step's incoming and staging directories, then the shared `<root>/.incoming` if it is empty (card
    C-12). Other steps may be publishing into `.incoming` meanwhile, so it is never listed and then removed: `rmdir`
    under the exclusive lock (`_incoming_lock`) removes it only when it is empty and no step is between creating it
    and creating its own directory in it; a non-empty `.incoming` is left to the last step to finish. Where the
    lock cannot be held, `.incoming` is left in place (empty, and outside every manifest)."""
    for d in (incoming, stage_dir):
        if d.exists():
            shutil.rmtree(d)
    if not incoming.parent.exists():
        return
    with _incoming_lock(incoming.parent.parent, exclusive=True) as held:
        if held:
            with contextlib.suppress(OSError):  # not empty (another step is publishing), or removed already
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
        _err(f"{intent['stage']}.{intent['step']} recovered an interrupted publish of {len(intent['outputs'])} files")
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
    _err(f"{rec['stage']}.{rec['step']} skipped: run record and outputs present")
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

    rsync = which("rsync")
    _make_incoming(incoming)
    if rsync:
        subprocess.run([rsync, "-a", f"{stage_dir}/", f"{incoming}/"], check=True)
    else:
        shutil.copytree(stage_dir, incoming, dirs_exist_ok=True)
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


# The measuring helper's sampling interval in seconds (card C-16). Each sample reads one file (a scope's memory.peak)
# or two small files per process of the step (its /proc status and its threads' children lists), and the helper
# writes only when the peak has grown; so a fifth of a second costs little. A step whose helper is killed loses what
# it grew in its last 0.2 s: the `oom` fault, which only touches pages, touched 4 GB in a scope on the pod in 4.5 to
# 5.2 s, its imports included (measured, wilkcomp02, 2026-11-10, card A-07), so at least 0.77 GB/s and at least
# 0.15 GB in 0.2 s (derived); a real step grows more slowly (estimate).
SAMPLE_S = 0.2

# The step's measuring helper (card C-15), run as `python -I -S -c _MEASURE OUT CGROUP_FILE CGROUP_ROOT -- CMD...`
# inside the step's limits (after `nice`, so inside a scope run-step opens). Linux carries a process's memory image
# at `exec` into the `ru_maxrss` of what it execs (fs/exec.c, `exec_mmap`), and a child forked by run-step starts
# with run-step's image, so `wait4` on that child reports at least the launcher's peak: every record of the B-30
# probe read 344.0 MB (measured, wilkcomp02, 2026-11-10). This helper is a fresh, small interpreter (stdlib only,
# no site): it starts the command, waits for it with `wait4` and writes, to OUT, the command tree's `ru_maxrss`
# (the largest single process of it, with the helper's own image as a floor) and, when CGROUP_FILE is given, the
# `memory.peak` of the cgroup it runs in, with `"final": true`. It then ends as the command ended (same exit code,
# or killed by the same signal). SIGTERM, SIGINT and SIGHUP reach the command through its process group, so the
# helper only waits.
# Card C-16: the file is written as the command runs, too, so that a helper killed with its command (a wall timeout
# or a stopped scope kills the whole process group; nothing a process does after SIGKILL) leaves its last sample.
# Every SAMPLE_S seconds a thread reads the scope's `memory.peak` (with CGROUP_FILE; the kernel's high-water mark of
# the cgroup) or else the largest `VmHWM` (the kernel's high-water mark of one process's resident memory) over the
# command's process tree, found through /proc/<pid>/task/<tid>/children, and rewrites OUT with `"final": false`
# when the value has grown. Without /proc (macOS) nothing is sampled. The command is started by fork and exec, not
# posix_spawn, so that on Linux it gets PR_SET_PDEATHSIG (SIGKILL) first: it dies when the helper dies. The search
# runner kills the helper when it dies itself (`search.runner._die_with_runner`), so a worker never outlives its
# runner, even one still importing, before its own parent watch (`search.worker._watch_parent`) is armed.
_MEASURE = r"""
import json, os, signal, sys, threading
out, cgroup_file, cgroup_root, cmd = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[5:]
SAMPLE_S = @SAMPLE_S@
GROUP = (signal.SIGTERM, signal.SIGINT, signal.SIGHUP)
for s in GROUP:
    signal.signal(s, lambda signum, frame: None)
peak_path = None
if cgroup_file:
    try:
        with open(cgroup_file) as f:
            rel = next(ln[3:].strip() for ln in f if ln.startswith("0::"))
        peak_path = os.path.join(cgroup_root, rel.lstrip("/"), "memory.peak")
    except (OSError, StopIteration):
        peak_path = None


def read_int(path):
    try:
        with open(path) as f:
            return int(f.read().strip())
    except (OSError, ValueError):
        return None


def tree_hwm(top):
    best, todo, seen = None, [top], set()
    while todo:
        p = todo.pop()
        if p in seen:
            continue
        seen.add(p)
        try:
            with open(f"/proc/{p}/status") as f:
                for ln in f:
                    if ln.startswith("VmHWM:"):
                        v = int(ln.split()[1]) * 1024
                        best = v if best is None or v > best else best
                        break
            tids = os.listdir(f"/proc/{p}/task")
        except (OSError, ValueError):
            continue
        for t in tids:
            try:
                with open(f"/proc/{p}/task/{t}/children") as f:
                    todo.extend(int(c) for c in f.read().split())
            except (OSError, ValueError):
                pass
    return best


def save(rec):
    with open(out + ".tmp", "w") as f:
        json.dump(rec, f)
    os.replace(out + ".tmp", out)


libc = None
if sys.platform.startswith("linux"):
    try:
        import ctypes
        libc = ctypes.CDLL(None)
    except (ImportError, OSError):
        libc = None
me = os.getpid()
pid = os.fork()
if pid == 0:
    try:
        if libc is not None:
            libc.prctl(1, int(signal.SIGKILL))  # PR_SET_PDEATHSIG: the command dies when this helper dies
            if os.getppid() != me:  # the helper died before the setting took
                os._exit(1)
        for s in (*GROUP, signal.SIGPIPE, signal.SIGXFSZ):
            signal.signal(s, signal.SIG_DFL)
        os.execvp(cmd[0], cmd)
    except BaseException as e:
        try:
            os.write(2, f"run-step: cannot run {cmd[0]}: {e}\n".encode())
        finally:
            os._exit(127)
seen = {"maxrss_bytes": None, "memory_peak_bytes": None}
stop = threading.Event()


def sample():
    key = "memory_peak_bytes" if peak_path else "maxrss_bytes"
    while not stop.wait(SAMPLE_S):
        v = read_int(peak_path) if peak_path else (tree_hwm(pid) if os.path.isdir("/proc/self") else None)
        if v is not None and (seen[key] is None or v > seen[key]):
            seen[key] = v
            save({**seen, "final": False})


sampler = threading.Thread(target=sample, daemon=True)
sampler.start()
_, status, ru = os.wait4(pid, 0)
stop.set()
sampler.join()
peak = read_int(peak_path) if peak_path else None
maxrss = ru.ru_maxrss if sys.platform == "darwin" else ru.ru_maxrss * 1024
save({"maxrss_bytes": maxrss, "memory_peak_bytes": peak, "final": True})
if os.WIFSIGNALED(status):
    sig = os.WTERMSIG(status)
    try:
        import resource
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))  # the helper's own death writes no core file
    except (ImportError, OSError, ValueError):
        pass
    try:
        signal.signal(sig, signal.SIG_DFL)
    except (OSError, ValueError):  # SIGKILL and SIGSTOP cannot be caught and need no reset
        pass
    os.kill(os.getpid(), sig)
    os._exit(128 + sig)
os._exit(os.waitstatus_to_exitcode(status))
""".replace("@SAMPLE_S@", repr(SAMPLE_S))


def measured(cmd: Sequence[str], out: Path, scope: bool, *, proc_root: str | Path | None = None,
             cgroup_root: str | Path | None = None) -> list[str]:
    """`cmd` under the measuring helper (`_MEASURE`), which writes the step's peak memory to `out`. With `scope`
    (run-step opened a scope for this step alone) it also reads that scope's `memory.peak`."""
    cgroup_file = str(Path(PROC_ROOT if proc_root is None else proc_root) / "self" / "cgroup") if scope else ""
    root = str(CGROUP_ROOT if cgroup_root is None else cgroup_root) if scope else ""
    return [sys.executable, "-I", "-S", "-c", _MEASURE, str(out), cgroup_file, root, "--", *cmd]


def step_peak(out: Path) -> tuple[int | None, str]:
    """The step's peak memory in bytes and where it comes from, from the helper's file `out` (removed after).
    Written when the command ended: `scope_memory_peak` (the memory.peak of the scope run-step opened for this step),
    else `child_maxrss` (the largest `ru_maxrss` in the step's process tree, measured from the helper). Card C-16: a
    helper killed before its command ended left only its last sample, taken at most `SAMPLE_S` seconds before:
    `scope_memory_peak_sampled` (the scope's memory.peak) or `child_hwm_sampled` (the largest VmHWM in the step's
    process tree). Else (None, `unmeasured`)."""
    try:
        got = json.loads(out.read_text())
    except (OSError, ValueError):
        return None, "unmeasured"
    finally:
        out.unlink(missing_ok=True)
        out.with_name(out.name + ".tmp").unlink(missing_ok=True)  # a sample a killed helper was writing
    tail ="" if got.get("final", True) else "_sampled"
    if got.get("memory_peak_bytes") is not None:
        return int(got["memory_peak_bytes"]), "scope_memory_peak" + tail
    if got.get("maxrss_bytes") is not None:
        return int(got["maxrss_bytes"]), ("child_maxrss" if not tail else "child_hwm_sampled")
    return None, "unmeasured"


def failed_path(tmp: str | Path, stage: str, step: str) -> Path:
    """The record of a step's last attempt whose command failed or was killed (card C-16):
    `<profile tmp>/staging/<stage>.<step>.failed.json`, beside the staging it left. Outside the root, so it is never a
    run record (`catalog.built_by` reads those) and never marks the step done; the next attempt removes it."""
    return Path(tmp) / "staging" / f"{stage}.{step}.failed.json"


def _execute(argv: list[str], env: dict[str, str]) -> tuple[int, float, float, int, int | None]:
    """Run argv in its own session (so its process group can be signalled as one; a systemd-run
    scope keeps that group); return (exit status, wall seconds, CPU minutes, `ru_maxrss` in bytes
    (at least run-step's own peak: card C-15; `step_peak` gives the step's), the signal forwarded
    or None). SIGTERM and SIGINT received meanwhile are forwarded to the group,
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
             proc_root: str | Path | None = None, cgroup_root: str | Path | None = None,
             budget_reason: str | None = None, _hook: Hook | None = None) -> int:
    """Run one step; return the exit code for `bio-box run-step` (0, 75 for STOP, else non-zero).
    `budget_reason`, when given, is written to the run record as `cpu_budget_reason` (card B-36:
    `build-stage` gives `exempt: DEC-29` for a clustering call that no number of parts fits, card B-45); `over_cpu_budget` is kept as it is.
    `proc_root` and `cgroup_root` (tests only) stand in for /proc and /sys/fs/cgroup.
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
    try:
        cap = memory_cap(prof.mem_scope, require=prof.require_mem_cap, system=system, which=which,
                         proc_root=proc_root, cgroup_root=cgroup_root)
    except NoCap as e:
        _err(f"not starting {stage}.{step}: {e}")
        return EXIT_REFUSED

    failed = failed_path(prof.tmp, stage, step)
    for stale in (incoming, stage_dir):
        if stale.exists():
            shutil.rmtree(stale)
    failed.unlink(missing_ok=True)  # an earlier attempt's (card C-16)
    stage_dir.mkdir(parents=True)

    env = dict(os.environ, BIOBOX_STAGING=str(stage_dir), BIOBOX_ROOT=str(root), BIOBOX_THREADS=str(n_threads))
    env.update({v: str(n_threads) for v in THREAD_VARS})
    inherited = cap if cap["how"] == "inherited" else None
    argv = build_argv(cmd, nice_value, prof.mem_scope, inherited=inherited, system=system, which=which)
    peak_file = stage_dir.parent / f"{stage}.{step}.peak.json"
    peak_file.unlink(missing_ok=True)
    opened = scope_argv(cmd, nice_value, mem_max=prof.mem_scope, inherited=inherited, system=system,
                        which=which)[1] == "scope"
    started_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    code, wall, cpu_min, _, sig = _execute(
        build_argv(measured(cmd, peak_file, opened, proc_root=proc_root, cgroup_root=cgroup_root), nice_value,
                   prof.mem_scope, inherited=inherited, system=system, which=which), env)
    rss, rss_source = step_peak(peak_file)
    over = cpu_min > prof.cpu_min_per_process
    meta = {
        "stage": stage, "step": step, "argv": argv, "threads": n_threads, "nice": nice_value,
        "started_at": started_at, "wall_s": round(wall, 3), "cpu_min": round(cpu_min, 6),
        "peak_rss_bytes": rss, "peak_rss_source": rss_source, "host": socket.gethostname(), "exit_code": code,
        "image_digest": os.environ.get("BIOBOX_IMAGE_DIGEST", "none"), "over_cpu_budget": over,
        "mem_cap": cap, **codeid.code_identity(),  # card C-23: the code that ran, beside the image digest
    }
    if budget_reason is not None:
        meta["cpu_budget_reason"] = budget_reason

    if sig is not None or code != 0:  # card C-16: the attempt's record, its peak memory with it
        _atomic_write(failed, _record_bytes({**meta, "forwarded_signal": signal.Signals(sig).name if sig else None}))
        peak = f"peak memory {rss / 2**20:.1f} MiB ({rss_source})" if rss is not None else "peak memory unmeasured"
    if sig is not None:
        _err(f"{signal.Signals(sig).name} received; forwarded to {stage}.{step}, which exited {code}; {peak}; "
             f"nothing published, staging left at {stage_dir}, attempt recorded in {failed}")
        return 128 + sig
    if code != 0:
        _err(f"{stage}.{step} failed with exit {code}; {peak}; nothing published, staging left at {stage_dir}, "
             f"attempt recorded in {failed}")
        return code if code > 0 else 1
    if (root / "STOP").exists():
        _err(f"{root / 'STOP'} appeared; not publishing {stage}.{step}, staging left at {stage_dir}")
        return EXIT_STOP
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
    _err(f"{stage}.{step} published {len(record['outputs'])} files")
    return 0
