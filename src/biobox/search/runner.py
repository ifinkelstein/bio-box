"""The `bio-box search` runner (operation O-RUN; SEARCH_PLAN SD-01, SD-15, SD-20; data_layer_and_tests.md §2B.3).

A run is a recorded job:

  1. The recipe is checked on paper (`recipe.load`) and against the release (`resolve`):
     the pack must verify (C-12), every seed, model ID and positive control must be in it,
     and nothing in the recipe may be beyond what this build runs. A file seed (`hmm_file`,
     `msa`, `sequences`) is read once here: its sha256 replaces its path in the recipe hash
     and the query hash, and its bytes are staged for the worker. Refusals raise
     RecipeError and record nothing.
  2. The target's shards are copied to `<work>/<run_id>/` and checked against the pack
     manifest's sha256; their total must fit the /tmp quota.
  3. `biobox.search.worker` runs in a child process started the way run-step starts a step
     (`runstep.scope_argv`, card C-05): on Linux with `systemd-run`, inside a transient user
     scope with `MemoryMax` at the memory cap, `MemorySwapMax=0` and `CPUQuota` at the CPU
     count, at nice 15; elsewhere under `nice` (or `os.nice`). A cgroup memory cap already above
     the runner (card C-07: inside the image, the unit `scripts/pod_run.sh` opened on the host)
     is inherited instead (`runstep.memory_cap`): no scope, and the worker's memory limit is
     that cap. With `Limits.require_cap` and neither a scope nor a cap, the run is refused
     before anything is recorded. RLIMIT_CPU at the CPU cap (300 CPU-min; the pod watchdog
     kills near 369) applies in every path, because the watchdog's rule is per process.
     RLIMIT_AS at the memory cap is set only on Linux where no cgroup limit holds the worker to
     `mem_gb` (no scope and no cap, or an inherited cap looser than `mem_gb`): it caps virtual
     memory, not resident memory, so under it a thread-heavy process fails to start threads
     long before it is out of memory (194 threads failed at 5 GB resident under an 18 GB
     RLIMIT_AS on the 112-core host). A wall-clock timeout covers all paths.
  4. The exit decides the status (`classify_exit`): 0 -> the worker's own (`ok`, `failed` for a
     control, `output_cap`); EXIT_OOM (MemoryError, under RLIMIT_AS in the fallback) -> `oom`;
     a SIGKILL the runner did not send -> `oom` only when the scope's cgroup counts an OOM
     kill (`scope_oom_killed`; in an inherited unit, one more than before the worker started),
     else `failed`, killed from outside (card W-04); SIGXCPU or the
     wall timeout -> `timeout`; a cancel (SIGTERM) -> `stopped`; EXIT_FAILED (an uncaught
     exception in any worker thread, traceback in stderr) and anything else -> `failed`.
     Only runs whose worker wrote its result (`ok`, `failed`-by-control, `output_cap`) carry
     hits; the others write `run.parquet` alone, never an empty or partial hits file (S-05).
  5. The run directory is written next to its final place and renamed in.

Where things live (card W-03): nothing in the search path writes under the data root, which a user
who can read a release may not be able to write. The pack is read from a pack directory (`packs`,
`biobox.search.pack`); a run's directory is `out` (`bio-box search run --out DIR`), else
`run=<run_id>/` in a runs directory (`runs`, `--runs DIR`). The pending record of a submitted run and
its runner's stderr are `pending/<run dir name>.json` and `.stderr` beside the run directory
(`pending_path`); the runner's lock is in its work directory. A run started with `bio-box search run
--wait` (`run(pending=True)`, card W-13) takes that lock and writes that record itself, so it is watched
and cancelled the same way. Defaults: with a host profile, a per-user
directory under the profile's `tmp` (`<tmp>/search/<user>/runs` and `.../work`) and `<tmp>/packs`;
without one, `default_work_base()` (`$BIOBOX_SEARCH_WORK`, else /tmp/atd/search), `runs` inside it and
`packs` beside it. A run, work, runs or pack directory under the data root is refused before anything
is read or recorded. A run is found by its directory (`status`, `cancel`, `read_run`, `read_hits`), or
by its id in a runs directory (`run_dir`).
"""

from __future__ import annotations

import copy
import fcntl
import getpass
import hashlib
import json
import os
import platform
import resource
import shutil
import signal
import socket
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from blake3 import blake3

from biobox import ids, provenance, runstep
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.registry.schemas import table
from biobox.search import engine, fastz, pack, recipe as rcp
from biobox.search.recipe import RecipeError
from biobox.search.worker import BOOT, CGROUP_FILE, EXIT_FAILED, EXIT_OOM, THREAD_VARS

GB = 10**9
CGROUP_ROOT = Path("/sys/fs/cgroup")  # the cgroup v2 mount
PROC_ROOT = Path("/proc")
PF_EXITING = 0x4  # task flag (include/linux/sched.h) the kernel sets once a task has begun to exit
_WHICH = object()  # default for scope_oom_killed(systemctl=...): look it up on PATH
LOCK_NAME = "runner.lock"  # flock held by the runner for its lifetime, through an fd inherited from `submit`
RUN_FILES = ("run.parquet", "hits.parquet", "profiles.parquet", "neighborhoods.parquet")


@dataclass(frozen=True)
class Limits:
    cpu_min: float = 300.0       # per search process (SD-15)
    wall_s: float = 6 * 3600.0
    mem_gb: float | None = None  # the scope's MemoryMax; RLIMIT_AS on Linux without systemd-run; else none
    tmp_quota_gb: float = 150.0
    nice: int = 15
    cpus: int = 8
    require_cap: bool = False    # card C-07: refuse to run with neither a scope to open nor a cgroup cap above


def _user_dir(profile_tmp: Path, user: str | None) -> Path:
    return Path(profile_tmp) / "search" / (user or getpass.getuser())


def default_work_base(profile_tmp: Path | None = None, user: str | None = None) -> Path:
    """Where runs' work directories go: with a host profile's `tmp`, `<tmp>/search/<user>/work`; else
    `$BIOBOX_SEARCH_WORK`, else /tmp/atd/search."""
    if profile_tmp is not None:
        return _user_dir(profile_tmp, user) / "work"
    return Path(os.environ.get("BIOBOX_SEARCH_WORK", "/tmp/atd/search"))


def default_runs_dir(profile_tmp: Path | None = None, user: str | None = None) -> Path:
    """Where run directories go by default: with a host profile's `tmp`, `<tmp>/search/<user>/runs` (one per
    user, so two users never see each other's runs); else `runs` in `default_work_base()`."""
    if profile_tmp is not None:
        return _user_dir(profile_tmp, user) / "runs"
    return default_work_base() / "runs"


def default_pack_base(profile_tmp: Path | None = None) -> Path:
    """The default pack directory: `<tmp>/packs` with a host profile (shared: a pack is a cache of a release),
    else `packs` beside `default_work_base()`."""
    if profile_tmp is not None:
        return Path(profile_tmp) / "packs"
    return default_work_base().parent / "packs"


def run_dir(runs: Path, run_id: bytes | str) -> Path:
    """The directory of run `run_id` in the runs directory `runs`."""
    return Path(runs) / f"run={run_id if isinstance(run_id, str) else run_id.hex()}"


def pending_path(d: Path) -> Path:
    """The pending record of the run whose directory is `d` (its runner's stderr has the suffix .stderr)."""
    d = Path(d)
    return d.parent / "pending" / f"{d.name}.json"


def _outside(root: Path, path: Path, what: str) -> Path:
    try:
        return pack.outside_root(root, path, what)
    except ValueError as e:
        raise RecipeError(str(e)) from None


def resolve(root: Path, r: dict, release: str, packs: Path) -> dict:
    """Submit-time checks against the release's pack in the pack directory `packs`; returns the job's seed
    and positives. A missing pack is refused with the command that builds it; it is never built here."""
    gaps = rcp.executable(r)
    if gaps:
        raise RecipeError("not executable in this build: " + "; ".join(gaps))
    try:
        pdir = pack.find(packs, release, root=root)
    except (FileNotFoundError, ValueError) as e:
        raise RecipeError(str(e)) from None
    problems = pack.verify_dir(root, release, pdir)
    if problems:
        raise RecipeError(f"search pack for {release} does not verify (C-12): {problems[:3]}")
    manifest, pid = pack.read_manifest(pdir)
    tname = r["target"]
    shards = [s for s in manifest["shards"] if s["target"] == tname]
    models = [s for s in manifest["shards"] if s["target"] == "models"]
    target_ids = {n for s in shards for n, _ in fastz.read_records(pdir / s["path"])}
    seed = dict(r["seed"][0])
    if seed["kind"] == "protein_uid":
        uniq = {n: q for s in manifest["shards"] if s["target"] == "unique"
                for n, q in fastz.read_records(pdir / s["path"])}
        if seed["value"] not in uniq:
            raise RecipeError(f"seed protein_uid {seed['value']} is not in release {release}")
        seed["sequence"] = uniq[seed["value"]]
    elif seed["kind"] == "sequence":
        seed["sequence"] = seed["value"]
    elif seed["kind"] == "hmm":
        if engine.find_model([pdir / m["path"] for m in models], seed["value"]) is None:
            raise RecipeError(f"model {seed['value']} is not in the release's libraries")
    seed_bytes = None
    if seed["kind"] in rcp.FILE_SEED_KINDS:
        seed, seed_bytes = _read_seed_file(seed)
        r = copy.deepcopy(r)
        r["seed"][0]["sha256"] = seed["sha256"]  # what recipe_hash covers in place of the path
    positives = []
    for p in r["controls"]["positive"]:
        uid = ids.protein_uid(p["value"]).hex() if p["kind"] == "sequence" else p["value"]
        if uid not in target_ids:
            raise RecipeError(f"positive control {p.get('label', p['value'][:20])} is not in target {tname}: "
                              "it could never be recovered")
        positives.append({"kind": p["kind"], "uid": uid, "label": p.get("label")})
    return {"manifest": manifest, "pack_id": pid, "shards": shards, "models": models, "seed": seed,
            "positives": positives, "pack_dir": pdir, "recipe": r, "seed_bytes": seed_bytes}


def _read_seed_file(seed: dict) -> tuple[dict, bytes]:
    """Read a file seed once at submit: its sha256 and what it holds travel with the job, and the
    worker searches with these bytes (staged), not with whatever the path holds later."""
    kind, path = seed["kind"], seed["value"]
    try:
        data = Path(path).read_bytes()
    except FileNotFoundError:
        raise RecipeError(f"seed {kind} file {path} not found") from None
    except OSError as e:
        raise RecipeError(f"seed {kind} file {path} cannot be read: {e}") from None
    out = {k: v for k, v in seed.items() if k != "value"}
    out.update(path=path, sha256=hashlib.sha256(data).hexdigest())
    try:
        if kind == "hmm_file":
            hmm = engine.read_hmm_file(data, seed.get("model"))
            out.update(model_name=engine.text(hmm.name), model_accession=engine.text(hmm.accession),
                       model_length=int(hmm.M))
        elif kind == "sequences":
            recs = engine.read_seed_fasta(data)
            out.update(n_sequences=len(recs), names=[n for n, _ in recs])
    except ValueError as e:
        raise RecipeError(f"seed {kind} file {path}: {e}") from None
    return out, data


def query_hash(seed: dict) -> str:
    """The resolved seed's hash; a file seed counts by its content (sha256 and what it holds), not its path."""
    q = {k: v for k, v in seed.items() if not (seed["kind"] in rcp.FILE_SEED_KINDS and k == "path")}
    return blake3(json.dumps(q, sort_keys=True).encode()).hexdigest()


def worker_command(job: Path) -> list[str]:
    """The worker's own entry (`worker.BOOT`): it sets the thread variables before importing anything."""
    return [sys.executable, "-c", BOOT, str(job)]


def _mem_max(limits: Limits) -> str | None:
    return None if limits.mem_gb is None else str(int(limits.mem_gb * GB))


def memory_cap(limits: Limits, *, system: str | None = None, which=None, proc_root: Path | None = None,
               cgroup_root: Path | None = None) -> dict:
    """How the worker's memory will be capped (`runstep.memory_cap`, card C-07): opened (a scope),
    inherited (a cgroup cap above the runner) or none. With `limits.require_cap` and neither,
    a RecipeError that names scripts/pod_run.sh."""
    try:
        return runstep.memory_cap(_mem_max(limits), require=limits.require_cap, system=system, which=which,
                                  proc_root=proc_root, cgroup_root=cgroup_root)
    except runstep.NoCap as e:
        raise RecipeError(str(e)) from None


def worker_argv(job: Path, limits: Limits, *, inherited: dict | None = None, system: str | None = None,
                which=shutil.which, peak_file: Path | None = None) -> tuple[list[str], str]:
    """The worker's argv and how it is limited ("scope", "inherited", "nice", "bare"), by
    run-step's helper; `inherited` is the cap `memory_cap` found above the runner, if any.
    With `peak_file` (card W-14) the worker runs under run-step's measuring helper
    (`runstep.measured`), inside the limits, which writes the worker's own peak memory there."""
    kw = dict(mem_max=_mem_max(limits), cpu_quota=limits.cpus * 100, inherited=inherited, system=system,
              which=which)
    cmd = worker_command(job)
    argv, kind = runstep.scope_argv(cmd, limits.nice, **kw)
    if peak_file is None:
        return argv, kind
    return runstep.scope_argv(runstep.measured(cmd, peak_file, kind == "scope"), limits.nice, **kw)[0], kind


def _tighter(limits: Limits, cap_bytes: int | None) -> bool:
    """`mem_gb` asks for less than the inherited cap gives (so the cap alone does not hold it)."""
    return limits.mem_gb is not None and (cap_bytes is None or int(limits.mem_gb * GB) < cap_bytes)


def rlimits(limits: Limits, kind: str, system: str, cap_bytes: int | None = None) -> dict[int, tuple[int, int]]:
    """The rlimits the worker starts under. RLIMIT_CPU always (per process); RLIMIT_AS only on
    Linux where no cgroup holds memory to `mem_gb`: no scope and no cap, or an inherited cap
    (`cap_bytes`) looser than `mem_gb` (module docstring, point 3)."""
    cpu = max(1, int(limits.cpu_min * 60))
    out = {resource.RLIMIT_CPU: (cpu, cpu + 5)}
    unheld = kind in ("nice", "bare") or (kind == "inherited" and _tighter(limits, cap_bytes))
    if limits.mem_gb is not None and system == "Linux" and unheld:
        b = int(limits.mem_gb * GB)
        out[resource.RLIMIT_AS] = (b, b)
    return out


def worker_env(limits: Limits, kind: str, system: str, base: dict[str, str],
               cap_bytes: int | None = None) -> dict[str, str]:
    """The worker's environment: `base`, every thread variable at `limits.cpus` (card W-04; the
    worker sets them again itself), and the memory limit it runs under, by name, for the `oom`
    fault to allocate past: BIOBOX_MEMORY_MAX in a scope or an inherited cap (`cap_bytes`),
    BIOBOX_RLIMIT_AS where RLIMIT_AS is the limit (`rlimits`)."""
    env = dict(base)
    env.update({v: str(max(1, limits.cpus)) for v in THREAD_VARS})
    env.pop("BIOBOX_MEMORY_MAX", None)
    env.pop("BIOBOX_RLIMIT_AS", None)
    if kind == "inherited":
        if _tighter(limits, cap_bytes) and system == "Linux":
            env["BIOBOX_RLIMIT_AS"] = str(int(limits.mem_gb * GB))
        elif cap_bytes is not None:
            env["BIOBOX_MEMORY_MAX"] = str(cap_bytes)
    elif limits.mem_gb is not None:
        if kind == "scope":
            env["BIOBOX_MEMORY_MAX"] = str(int(limits.mem_gb * GB))
        elif system == "Linux":
            env["BIOBOX_RLIMIT_AS"] = str(int(limits.mem_gb * GB))
    return env


def _die_with_runner() -> None:
    """Linux: SIGKILL this process when the runner dies. The worker watches its parent (`worker._watch_parent`),
    and since card W-14 its parent is run-step's measuring helper, which would outlive a dead runner; the
    setting survives the `exec` of the helper (through `systemd-run` and `nice`)."""
    if platform.system() != "Linux":
        return
    import ctypes
    parent = os.getppid()
    ctypes.CDLL(None, use_errno=True).prctl(1, signal.SIGKILL)  # PR_SET_PDEATHSIG
    if os.getppid() != parent:  # the runner died before the setting took
        os._exit(1)


def _preexec(limits: Limits, kind: str, system: str, cap_bytes: int | None = None, renice: bool = False):
    """Runs in the child before exec: its own session (so the runner signals the group), its
    niceness when no `nice` binary sets it (`renice`), and its rlimits (inherited through
    systemd-run and nice)."""
    lims = rlimits(limits, kind, system, cap_bytes)

    def fn():
        os.setsid()
        _die_with_runner()
        if renice:
            os.nice(limits.nice)
        for which, value in lims.items():
            resource.setrlimit(which, value)
    return fn


def _oom_kills(memory_events: str) -> int:
    for line in memory_events.splitlines():
        key, _, value = line.partition(" ")
        if key == "oom_kill":
            return int(value)
    return 0


def _oom_kills_now(cgroup: str, cgroup_root: Path = CGROUP_ROOT) -> int:
    """The `oom_kill` count of a cgroup's `memory.events` now; 0 when it cannot be read."""
    try:
        return _oom_kills((Path(cgroup_root) / cgroup.lstrip("/") / "memory.events").read_text())
    except (OSError, ValueError):
        return 0


def scope_oom_killed(work: Path, *, cgroup_root: Path = CGROUP_ROOT, systemctl=_WHICH,
                     baseline: int = 0) -> bool | None:
    """Whether the kernel OOM-killed a process of the worker's scope: True, False, or None when
    there is no evidence either way. The worker recorded its cgroup path (`CGROUP_FILE`); the
    scope's `memory.events` counts OOM kills while the cgroup exists, and a kill is counted when
    it is above `baseline` (in an inherited unit, shared with the runner and with whatever ran
    in it before, the count when the worker started; card C-07). systemd removes the cgroup
    of a scope once its last process has exited, so after that the scope unit's `Result` is
    asked instead (`oom-kill` after an OOM kill; a failed unit stays loaded until reset-failed,
    which is done here once it has been read)."""
    try:
        cg = (Path(work) / CGROUP_FILE).read_text().strip()
    except OSError:
        return None
    if not cg:
        return None
    try:
        return _oom_kills((Path(cgroup_root) / cg.lstrip("/") / "memory.events").read_text()) > baseline
    except (OSError, ValueError):
        pass
    exe = shutil.which("systemctl") if systemctl is _WHICH else systemctl
    unit = cg.rsplit("/", 1)[-1]
    if not exe or not unit.endswith(".scope"):
        return None
    try:
        out = subprocess.run([exe, "--user", "show", "-p", "LoadState", "-p", "Result", unit],
                             capture_output=True, text=True, timeout=10).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    props = dict(line.partition("=")[::2] for line in out.splitlines())
    if props.get("LoadState") != "loaded" or not props.get("Result"):
        return None
    try:
        subprocess.run([exe, "--user", "reset-failed", unit], capture_output=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        pass
    return props["Result"] == "oom-kill"


def _stage(res: dict, work: Path, limits: Limits) -> tuple[list[Path], list[Path], list[bytes]]:
    need = res["shards"] + res["models"]
    total = sum(s["bytes"] for s in need)
    if total > limits.tmp_quota_gb * GB:
        raise _RunFailure("failed", f"staging {total} bytes exceeds the /tmp quota of {limits.tmp_quota_gb} GB")
    staged, sums = {}, []
    for s in need:
        dst = work / "pack" / s["path"]
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(res["pack_dir"] / s["path"], dst)
        got = ids.sha256_file(dst)
        if got.hex() != s["sha256"]:
            raise _RunFailure("failed", f"staged {s['path']} sha256 differs from the pack manifest")
        staged[s["path"]] = dst
        sums.append(got)
    return ([staged[s["path"]] for s in res["shards"]], [staged[s["path"]] for s in res["models"]], sums)


class _RunFailure(Exception):
    def __init__(self, status: str, error: str):
        super().__init__(error)
        self.status, self.error = status, error


def _wait(proc: subprocess.Popen, wall_s: float, cancelled: list) -> tuple[int, resource.struct_rusage | None, bool]:
    """(returncode, the child's rusage, whether the wall timeout fired)."""
    deadline = time.monotonic() + wall_s
    timed_out = False
    while True:
        pid, st, ru = os.wait4(proc.pid, os.WNOHANG)
        if pid:
            proc.returncode = os.waitstatus_to_exitcode(st)
            return proc.returncode, ru, timed_out
        if cancelled == ["requested"] and not timed_out:
            os.killpg(proc.pid, signal.SIGTERM)
            cancelled[0] = "sent"
        if time.monotonic() > deadline and not timed_out:
            timed_out = True
            os.killpg(proc.pid, signal.SIGKILL)
        time.sleep(0.02)


def classify_exit(rc: int, *, timed_out: bool, cancelled: bool, cgroup_oom: bool | None,
                  result: dict | None, stderr: str = "") -> tuple[str, str | None]:
    """(status, error) of a worker exit, from the evidence alone: its exit status `rc`; whether
    the runner itself sent the kill (`timed_out`: SIGKILL at the wall deadline; `cancelled`:
    SIGTERM); whether the scope's cgroup counted an OOM kill (`cgroup_oom`, None without a scope
    or without an answer); the worker's `result.json` (exit 0 only); its stderr. An EXIT_OOM is
    a MemoryError (under RLIMIT_AS in the fallback). A SIGKILL the runner did not send is `oom`
    only when the cgroup says so (card W-04): without that evidence it is `failed`."""
    if timed_out:
        return "timeout", "wall-clock timeout"
    if cancelled:
        return "stopped", "cancelled"
    if rc == 0 and result is not None:
        return result["status"], result["error"]
    if rc == EXIT_OOM:
        return "oom", "out of memory (MemoryError)"
    if rc == -signal.SIGKILL:
        if cgroup_oom:
            return "oom", "killed by the kernel OOM killer at the scope's MemoryMax (memory.events oom_kill)"
        why = "its scope reports no OOM kill" if cgroup_oom is False else "no evidence of an out-of-memory kill"
        return "failed", f"worker was killed by signal 9 from outside ({why})"
    if rc == -signal.SIGXCPU:
        return "timeout", "CPU limit reached (RLIMIT_CPU)"
    tail = stderr.strip().splitlines()[-1:] if stderr.strip() else []
    what = "worker raised an uncaught exception" if rc == EXIT_FAILED else f"worker exited {rc}"
    return "failed", what + (f": {tail[0]}" if tail else "")


def _status(rc: int, timed_out: bool, cancelled: bool, out: Path, stderr: str,
            cgroup_oom: bool | None = None) -> tuple[str, str | None, dict]:
    res = json.loads((out / "result.json").read_text()) if rc == 0 and (out / "result.json").exists() else None
    status, error = classify_exit(rc, timed_out=timed_out, cancelled=cancelled, cgroup_oom=cgroup_oom,
                                  result=res, stderr=stderr)
    return status, error, (res if res is not None and not (timed_out or cancelled) else {})


def new_run_id(recipe_hash: str, release: str, pack_id: bytes, submitted: datetime) -> bytes:
    return blake3(f"{recipe_hash}|{release}|{pack_id.hex()}|{submitted.isoformat()}|".encode()
                  + os.urandom(8)).digest(length=16)


def run(root: Path, recipe_src: Path | str, release: str, *, limits: Limits = Limits(),
        work_base: Path | None = None, argv: list[str] | None = None, run_id: bytes | None = None,
        fault: str | None = None, gene_call_version: str = GENE_CALL_VERSION, keep_work: bool = False,
        tape_version: str | None = None, packs: Path | None = None, runs: Path | None = None,
        out: Path | None = None, pending: bool = False, on_pending=None) -> dict:
    """Run one recipe to completion; returns the search_run row. Refusals raise RecipeError.

    The pack is read from the pack directory `packs` (default `default_pack_base()`); the run is written to
    the directory `out`, else to `run=<run_id>` in `runs` (default `default_runs_dir()`). None of these, nor
    `work_base`, may lie under the root (module docstring). The run reads the tape version of its pack (card
    B-16); `tape_version`, when given, must be it.

    `pending` (card W-13; `bio-box search run --wait`): the run takes its runner lock and writes the pending
    record `submit` writes, with this process as the runner, before its worker starts, so that `status` and
    `cancel` work on it from another process while it runs. A run id whose lock a live runner holds is refused.
    `on_pending(run_id_hex, path)` is called once the record is written. The runner `submit` starts runs
    without it: `submit` writes that record itself."""
    system = platform.system()
    cap = memory_cap(limits, system=system)  # refused (card C-07) before anything is read or recorded
    root = Path(root).resolve()
    packs = _outside(root, Path(packs) if packs is not None else default_pack_base(), "the pack directory")
    work_base = _outside(root, Path(work_base) if work_base is not None else default_work_base(), "the work directory")
    if out is None:
        runs = _outside(root, Path(runs) if runs is not None else default_runs_dir(), "the runs directory")
    else:
        out = _outside(root, Path(out), "the run directory")
        if out.exists():
            raise RecipeError(f"{out} exists: a run directory is written once")
    r, text = rcp.load(Path(recipe_src) if isinstance(recipe_src, Path) else recipe_src)
    res = resolve(root, r, release, packs)
    pack_tv = pack.manifest_tape_version(res["manifest"])
    if tape_version is not None and tape_version != pack_tv:
        raise RecipeError(f"tape version {tape_version} is not the tape version {pack_tv} of the {release} pack: "
                          "a run reads the one tape version its pack's flanks refer to")
    r = res["recipe"]
    submitted = datetime.now(timezone.utc)
    rh = rcp.recipe_hash(r)
    run_id = run_id or new_run_id(rh, release, res["pack_id"], submitted)
    final = out if out is not None else run_dir(runs, run_id)
    if final.exists():
        raise RecipeError(f"{final} exists: a run directory is written once")
    work = work_base / run_id.hex()
    lock_fd = _take_lock(work, run_id.hex()) if pending else None  # before the directory is cleaned
    if work.exists():  # `submit` made this directory and holds runner.lock in it: keep that file
        for child in work.iterdir():
            if child.name != LOCK_NAME:
                shutil.rmtree(child) if child.is_dir() else child.unlink()
    out = work / "out"
    out.mkdir(parents=True)
    cancelled: list = []
    prev = signal.signal(signal.SIGTERM, lambda *_: cancelled.append("requested")) \
        if hasattr(signal, "SIGTERM") and _main_thread() else None
    t0 = time.monotonic()
    sums: list[bytes] = []
    rc, ru, result = None, None, {}
    try:
        if pending:  # after the SIGTERM handler: from here a cancel is recorded as `stopped`
            info = _write_pending(final, run_id.hex(), os.getpid(), work / LOCK_NAME)
            if on_pending is not None:
                on_pending(run_id.hex(), info)
        shards, models, sums = _stage(res, work, limits)
        seed_file = None
        if res["seed_bytes"] is not None:  # the bytes read and hashed at submit, not the path's later content
            seed_file = work / "seed" / Path(res["seed"]["path"]).name
            seed_file.parent.mkdir()
            seed_file.write_bytes(res["seed_bytes"])
        # `release`: the worker reads only this release's files of the root (card B-28)
        job = {"run_id": run_id.hex(), "root": str(root), "release": release, "out": str(out), "recipe": r,
               "seed": res["seed"],
               "seed_file": str(seed_file) if seed_file else None,
               "positives": res["positives"], "target_kind": rcp.RUN_TARGETS[r["target"]],
               "shards": [str(p) for p in shards], "models": [str(p) for p in models],
               "meta": str(res["pack_dir"] / f"{r['target']}/{r['target']}_meta.parquet") if r["target"] != "unique" else None,
               "gene_call_version": gene_call_version, "tape_version": pack_tv, "cpus": limits.cpus,
               "max_output_bytes": int(r["max_output_gb"] * GB), "fault": fault,
               "negative_seed": int.from_bytes(run_id[:4], "big")}
        (work / "job.json").write_text(json.dumps(job, sort_keys=True))
        inherited = cap if cap["how"] == "inherited" else None
        cap_bytes = inherited["memory_max_bytes"] if inherited else None
        argv, kind = worker_argv(work / "job.json", limits, inherited=inherited, system=system)
        peak_file = work / "peak.json"  # the worker's own peak (card W-14), not `wait4`'s on a fork of this runner
        launch, _ = worker_argv(work / "job.json", limits, inherited=inherited, system=system, peak_file=peak_file)
        env = worker_env(limits, kind, system, os.environ, cap_bytes=cap_bytes)
        # an inherited unit is shared with the runner: only kills after this point are the worker's
        baseline = _oom_kills_now(inherited["cgroup"]) if inherited else 0
        renice = kind == "bare" or (kind == "inherited" and argv[0] != "nice")
        with open(work / "stderr.txt", "w") as err:
            proc = subprocess.Popen(launch, stdout=subprocess.DEVNULL, stderr=err,
                                    preexec_fn=_preexec(limits, kind, system, cap_bytes, renice), env=env)
            rc, ru, timed_out = _wait(proc, limits.wall_s, cancelled)
        # only a capped cgroup is evidence (the scope's own, or the unit's the cap was inherited
        # from); outside one, the cgroup is the runner's session's
        cgroup_oom = (scope_oom_killed(work, baseline=baseline)
                      if kind in ("scope", "inherited") and rc == -signal.SIGKILL and not timed_out else None)
        status, error, result = _status(rc, timed_out, bool(cancelled), out, (work / "stderr.txt").read_text(),
                                        cgroup_oom)
    except _RunFailure as f:
        status, error = f.status, f.error
    finally:
        if prev is not None:
            signal.signal(signal.SIGTERM, prev)
    wall = time.monotonic() - t0
    cpu_min = (ru.ru_utime + ru.ru_stime) / 60 if ru is not None else None
    rss = None
    if ru is not None:  # the helper's file; if the worker never ran or was killed without it, `unmeasured`
        rss = runstep.step_peak(work / "peak.json")[0]
    row = {
        "run_id": run_id, "recipe_hash": rh, "recipe_yaml": text, "argv": list(argv or []),
        "tool_versions_json": json.dumps(provenance.tool_versions("bio-box", "pyhmmer", "pyarrow"), sort_keys=True),
        "runtime_image_digest": os.environ.get("BIOBOX_IMAGE_DIGEST", "none"),
        "written_against": r["written_against"], "ran_against": release, "pack_id": res["pack_id"],
        "shard_sha256": sums, "z_prot": pack.Z["z_prot"], "z_nt_mb": pack.Z["z_nt_mb"],
        "query_hash": query_hash(res["seed"]),
        "submitted_at": submitted, "wall_s": wall, "peak_rss_bytes": rss, "cpu_min_max_pid": cpu_min,
        "host": socket.gethostname(), "exit_code": rc, "status": status, "partial_policy": r["partial"],
        "max_output_gb": r["max_output_gb"], "target": r["target"], "target_size": result.get("target_size"),
        "rounds": result.get("rounds"), "final_round": result.get("final_round"),
        "rounds_json": json.dumps(result["rounds_log"], default=str) if result.get("rounds_log") else None,
        "stop_reason": result.get("stop_reason"),
        "controls_json": json.dumps(result["controls"], sort_keys=True, default=str) if result.get("controls") else None,
        "warnings_json": json.dumps(result.get("warnings", [])) if result else None,
        "limits_json": json.dumps({**asdict(limits), "mem_cap": cap}, sort_keys=True), "error": error,
        "seed_json": json.dumps(res["seed"], sort_keys=True),
    }
    try:
        # a worker that died after a round has hits on disk, but no result: they are not the run's
        _publish(final, row, out, with_outputs=bool(result) and status in ("ok", "failed", "output_cap"))
        if not keep_work:
            shutil.rmtree(work, ignore_errors=True)
    finally:
        if lock_fd is not None:  # held until the run is closed: until then `status` says `running`
            os.close(lock_fd)
    return row


def _main_thread() -> bool:
    import threading

    return threading.current_thread() is threading.main_thread()


def _publish(final: Path, row: dict, out: Path, with_outputs: bool) -> None:
    """Write the run directory `final` beside its place under a temporary name and rename it in."""
    final.parent.mkdir(parents=True, exist_ok=True)
    tmp = final.with_name(f".{final.name}.tmp-{os.getpid()}")
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist([row], schema=table("search_run").schema), tmp / "run.parquet")
    if with_outputs:
        for name in RUN_FILES[1:]:
            if (out / name).exists():
                shutil.copyfile(out / name, tmp / name)
    if final.exists():
        raise FileExistsError(f"{final} exists: a run directory is never edited after it closes")
    os.replace(tmp, final)


# --- async submit, status, cancel ----------------------------------------------


def submit(root: Path, recipe_path: Path, release: str, limits: Limits = Limits(),
           work_base: Path | None = None, fault: str | None = None, *, packs: Path | None = None,
           runs: Path | None = None, out: Path | None = None) -> str:
    """Check the recipe now, run it detached, return the run_id at once. `packs`, `runs`, `out` and
    `work_base` as for `run`; the pending record is `pending_path` of the run's directory."""
    memory_cap(limits)  # card C-07: refuse here, before a runner is started or anything is recorded
    root = Path(root).resolve()
    packs = _outside(root, Path(packs) if packs is not None else default_pack_base(), "the pack directory")
    work_base = _outside(root, Path(work_base) if work_base is not None else default_work_base(), "the work directory")
    if out is not None:
        out = _outside(root, Path(out), "the run directory")
    else:
        runs = _outside(root, Path(runs) if runs is not None else default_runs_dir(), "the runs directory")
    r, _ = rcp.load(Path(recipe_path))
    res = resolve(root, r, release, packs)
    run_id = new_run_id(rcp.recipe_hash(res["recipe"]), release, res["pack_id"], datetime.now(timezone.utc))
    final = (out if out is not None else run_dir(runs, run_id)).resolve()
    if final.exists():
        raise RecipeError(f"{final} exists: a run directory is written once")
    info = pending_path(final)
    info.parent.mkdir(parents=True, exist_ok=True)
    # `--submitted`: this function writes the runner's pending record, so the runner does not (card W-13)
    cmd = [sys.executable, "-m", "biobox.cli", "search", "run", str(recipe_path), "--root", str(root),
           "--release", release, "--wait", "--submitted", "--run-id", run_id.hex(), "--cpus", str(limits.cpus),
           "--cpu-min", str(limits.cpu_min), "--wall-s", str(limits.wall_s), "--work", str(work_base.resolve()),
           "--pack-dir", str(packs.resolve()), "--out", str(final)]
    if fault is not None:
        cmd += ["--fault", fault]
    # The runner's liveness is a flock on local disk (never under the data root: flock is unreliable on
    # NFS). It is taken here, before the pending record exists, and handed to the runner as an inherited
    # fd (the lock belongs to the open file description, so it lives exactly as long as the runner's
    # files are open). The pid and its start time are recorded too: a killed multi-threaded runner's files
    # stay open, lock included, until its last thread has torn the process down, while its leader is
    # already a zombie (card A-07); `status` counts that runner as dead (`_runner_live`).
    work = (work_base / run_id.hex()).resolve()
    work.mkdir(parents=True, exist_ok=True)
    lock_path = work / LOCK_NAME
    fd = os.open(lock_path, os.O_RDWR | os.O_CREAT, 0o644)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=open(info.with_suffix(".stderr"), "w"),
                                start_new_session=True, pass_fds=(fd,))
    finally:
        os.close(fd)  # the runner's copy keeps the lock
    _write_pending(final, run_id.hex(), proc.pid, lock_path)  # our unreaped child: its pid cannot be reused yet
    return run_id.hex()


def _take_lock(work: Path, run_id: str) -> int:
    """Take the runner lock of a run started with `pending` (card W-13), as `submit` takes it for its runner: an
    exclusive flock on `<work>/runner.lock`, held through the returned fd for as long as the run is open. The fd
    is not inherited by the worker (Python opens it close-on-exec), so only the runner holds it (card A-07). A
    lock that another runner holds is refused: that run id is running."""
    work.mkdir(parents=True, exist_ok=True)
    fd = os.open(work / LOCK_NAME, os.O_RDWR | os.O_CREAT, 0o644)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        os.close(fd)
        raise RecipeError(f"run {run_id} is running: its runner holds {work / LOCK_NAME}") from None
    return fd


def _write_pending(final: Path, run_id: str, pid: int, lock: Path) -> Path:
    """Write the pending record of the run whose directory is `final` (`pending_path`): its runner's pid and
    start time and the lock that runner holds; written under a temporary name and renamed into place."""
    info = pending_path(final)
    info.parent.mkdir(parents=True, exist_ok=True)
    tmp = info.with_name(f".{info.name}.tmp")
    tmp.write_text(json.dumps({"run_id": run_id, "pid": pid, "start": proc_start_time(pid), "lock": str(lock),
                               "status": "running"}))
    os.replace(tmp, info)
    return info


def status(d: Path) -> dict:
    """The status of the run whose directory is `d`: its record once closed, else from its pending record
    (`pending_path`); KeyError when it has neither."""
    d = Path(d)
    if (d / "run.parquet").exists():
        row = pq.read_table(d / "run.parquet").to_pylist()[0]
        return {"run_id": row["run_id"].hex(), "status": row["status"], "error": row["error"], "rounds": row["rounds"],
                "final_round": row["final_round"], "stop_reason": row["stop_reason"]}
    p = pending_path(d)
    if p.exists():
        info = json.loads(p.read_text())
        run_id = info["run_id"]
        if _runner_live(info):
            return {"run_id": run_id, "status": "running"}
        if (d / "run.parquet").exists():  # closed between the first look and the probes
            return status(d)
        return {"run_id": run_id, "status": "failed", "error": "runner exited without closing the run"}
    raise KeyError(f"no run at {d}")


def _runner_live(info: dict) -> bool:
    """A runner is live when it holds its lock and its process is neither gone, nor exiting, nor another
    process (`runner_process_alive`). The lock alone is not enough: after SIGKILL the leader of a
    multi-threaded runner is a zombie, and gone to `ps`, while its last thread still holds the lock
    (25-47 ms of teardown on the pod; card A-07). The pid alone is not enough either (card A-01)."""
    return _runner_holds(info.get("lock")) and runner_process_alive(info) is not False


def _proc_stat(pid: int, proc_root: Path) -> tuple[str, int, int] | None:
    """(state, flags, start time) from /proc/<pid>/stat (proc(5) fields 3, 9 and 22), or None."""
    try:
        s = (Path(proc_root) / str(int(pid)) / "stat").read_text()
    except (OSError, ValueError):
        return None
    f = s[s.rindex(")") + 2:].split()  # the command name may hold spaces and parentheses
    return f[0], int(f[6]), int(f[19])


def proc_start_time(pid: int, proc_root: Path = PROC_ROOT) -> int | None:
    """The kernel's start time of a pid (clock ticks after boot), which tells a pid's owners apart; None without /proc."""
    st = _proc_stat(pid, proc_root)
    return None if st is None else st[2]


def runner_process_alive(info: dict, proc_root: Path = PROC_ROOT) -> bool | None:
    """Whether the recorded runner process (`pid`, and `start` when recorded) is alive and not exiting:
    False when the pid is gone, a zombie (Z, X), has begun to exit (PF_EXITING: SIGKILL delivered, the
    kernel tearing it down) or now belongs to a process with another start time; None without /proc
    (not Linux), where the lock alone decides."""
    if not (Path(proc_root) / "self" / "stat").exists():
        return None
    st = _proc_stat(info["pid"], proc_root)
    if st is None:
        return False
    state, flags, start = st
    if state in ("Z", "X", "x") or flags & PF_EXITING:
        return False
    return info.get("start") is None or start == info["start"]


def _runner_holds(lock: str | None) -> bool:
    """Whether a live runner holds its lock. If the probe gets the lock the runner is dead; it is released
    at once. A lock file that is gone means the runner finished and cleaned up (or never started)."""
    if not lock:
        return False
    try:
        fd = os.open(lock, os.O_RDWR)
    except FileNotFoundError:
        return False
    try:
        try:  # shared, so two probes at once never block each other into a false "running"
            fcntl.flock(fd, fcntl.LOCK_SH | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(fd, fcntl.LOCK_UN)
        return False
    finally:
        os.close(fd)


def cancel(d: Path) -> None:
    """Ask the live runner of the run whose directory is `d` to stop (its pending record says which). A
    runner that is gone is not signalled: its pid may belong to another process by now."""
    info = json.loads(pending_path(d).read_text())
    if _runner_live(info):
        os.kill(info["pid"], signal.SIGTERM)


def read_run(d: Path) -> dict:
    """The search_run row of the closed run whose directory is `d`."""
    return pq.read_table(Path(d) / "run.parquet").to_pylist()[0]


def read_hits(d: Path, round_: int | None = None, final_only: bool = True) -> list[dict]:
    """The hits of the run whose directory is `d`: by default the kept hits of its final round."""
    d = Path(d)
    if not (d / "hits.parquet").exists():
        return []
    rows = pq.read_table(d / "hits.parquet").to_pylist()
    if round_ is None and final_only:
        round_ = read_run(d)["final_round"]
        rows = [r for r in rows if r["curation_label"] == "keep"]
    if round_ is not None:
        rows = [r for r in rows if r["round"] == round_]
    return rows
