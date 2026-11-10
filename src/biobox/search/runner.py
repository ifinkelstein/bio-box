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
     scope with `MemoryMax` at the memory cap and `CPUQuota` at the CPU count, at nice 15;
     elsewhere under `nice` (or `os.nice`). RLIMIT_CPU at the CPU cap (300 CPU-min; the pod
     watchdog kills near 369) applies in both paths, because the watchdog's rule is per process.
     RLIMIT_AS at the memory cap is set only in the fallback on Linux, where no cgroup limit
     exists: it caps virtual memory, not resident memory, so under it a thread-heavy process
     fails to start threads long before it is out of memory (194 threads failed at 5 GB
     resident under an 18 GB RLIMIT_AS on the 112-core host). A wall-clock timeout covers both.
  4. The exit decides the status (`classify_exit`): 0 -> the worker's own (`ok`, `failed` for a
     control, `output_cap`); EXIT_OOM (MemoryError, under RLIMIT_AS in the fallback) -> `oom`;
     a SIGKILL the runner did not send -> `oom` only when the scope's cgroup counts an OOM
     kill (`scope_oom_killed`), else `failed`, killed from outside (card W-04); SIGXCPU or the
     wall timeout -> `timeout`; a cancel (SIGTERM) -> `stopped`; EXIT_FAILED (an uncaught
     exception in any worker thread, traceback in stderr) and anything else -> `failed`.
     Only runs whose worker wrote its result (`ok`, `failed`-by-control, `output_cap`) carry
     hits; the others write `run.parquet` alone, never an empty or partial hits file (S-05).
  5. The run directory is written next to its final place and renamed in
     (`search/run=<run_id>/`), then `search/run_index.json` is rewritten by atomic rename
     under a lock.
"""

from __future__ import annotations

import copy
import fcntl
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


def default_work_base() -> Path:
    return Path(os.environ.get("BIOBOX_SEARCH_WORK", "/tmp/atd/search"))


def runs_dir(root: Path) -> Path:
    return Path(root) / "search"


def run_dir(root: Path, run_id: bytes | str) -> Path:
    return runs_dir(root) / f"run={run_id if isinstance(run_id, str) else run_id.hex()}"


def resolve(root: Path, r: dict, release: str) -> dict:
    """Submit-time checks against the release's pack; returns the job's seed and positives."""
    gaps = rcp.executable(r)
    if gaps:
        raise RecipeError("not executable in this build: " + "; ".join(gaps))
    problems = pack.verify(root, release)
    if problems:
        raise RecipeError(f"search pack for {release} does not verify (C-12): {problems[:3]}")
    manifest, pid = pack.read_manifest(root, release)
    pdir = pack.pack_dir(root, release)
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


def worker_argv(job: Path, limits: Limits, *, system: str | None = None,
                which=shutil.which) -> tuple[list[str], str]:
    """The worker's argv and how it is limited ("scope", "nice", "bare"), by run-step's helper."""
    mem = None if limits.mem_gb is None else str(int(limits.mem_gb * GB))
    return runstep.scope_argv(worker_command(job), limits.nice, mem_max=mem, cpu_quota=limits.cpus * 100,
                              system=system, which=which)


def rlimits(limits: Limits, kind: str, system: str) -> dict[int, tuple[int, int]]:
    """The rlimits the worker starts under. RLIMIT_CPU always (per process); RLIMIT_AS only on
    Linux where no scope caps memory (module docstring, point 3)."""
    cpu = max(1, int(limits.cpu_min * 60))
    out = {resource.RLIMIT_CPU: (cpu, cpu + 5)}
    if limits.mem_gb is not None and system == "Linux" and kind != "scope":
        b = int(limits.mem_gb * GB)
        out[resource.RLIMIT_AS] = (b, b)
    return out


def worker_env(limits: Limits, kind: str, system: str, base: dict[str, str]) -> dict[str, str]:
    """The worker's environment: `base`, every thread variable at `limits.cpus` (card W-04; the
    worker sets them again itself), and the memory limit it runs under, by name, for the `oom`
    fault to allocate past (BIOBOX_MEMORY_MAX in a scope, BIOBOX_RLIMIT_AS in the fallback)."""
    env = dict(base)
    env.update({v: str(max(1, limits.cpus)) for v in THREAD_VARS})
    env.pop("BIOBOX_MEMORY_MAX", None)
    env.pop("BIOBOX_RLIMIT_AS", None)
    if limits.mem_gb is not None:
        if kind == "scope":
            env["BIOBOX_MEMORY_MAX"] = str(int(limits.mem_gb * GB))
        elif system == "Linux":
            env["BIOBOX_RLIMIT_AS"] = str(int(limits.mem_gb * GB))
    return env


def _preexec(limits: Limits, kind: str, system: str):
    """Runs in the child before exec: its own session (so the runner signals the group), its
    niceness when no `nice` binary sets it, and its rlimits (inherited through systemd-run and nice)."""
    lims = rlimits(limits, kind, system)

    def fn():
        os.setsid()
        if kind == "bare":
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


def scope_oom_killed(work: Path, *, cgroup_root: Path = CGROUP_ROOT, systemctl=_WHICH) -> bool | None:
    """Whether the kernel OOM-killed a process of the worker's scope: True, False, or None when
    there is no evidence either way. The worker recorded its cgroup path (`CGROUP_FILE`); the
    scope's `memory.events` counts OOM kills while the cgroup exists. systemd removes the cgroup
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
        return _oom_kills((Path(cgroup_root) / cg.lstrip("/") / "memory.events").read_text()) > 0
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


def _index_update(root: Path, entry: dict) -> None:
    d = runs_dir(root)
    d.mkdir(parents=True, exist_ok=True)
    with open(d / ".index.lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        path = d / "run_index.json"
        idx = json.loads(path.read_text()) if path.exists() else {"runs": []}
        idx["runs"] = sorted([e for e in idx["runs"] if e["run_id"] != entry["run_id"]] + [entry],
                             key=lambda e: e["run_id"])
        tmp = d / f".run_index.json.tmp-{os.getpid()}"
        tmp.write_text(json.dumps(idx, sort_keys=True, indent=1))
        os.replace(tmp, path)


def new_run_id(recipe_hash: str, release: str, pack_id: bytes, submitted: datetime) -> bytes:
    return blake3(f"{recipe_hash}|{release}|{pack_id.hex()}|{submitted.isoformat()}|".encode()
                  + os.urandom(8)).digest(length=16)


def run(root: Path, recipe_src: Path | str, release: str, *, limits: Limits = Limits(),
        work_base: Path | None = None, argv: list[str] | None = None, run_id: bytes | None = None,
        fault: str | None = None, gene_call_version: str = GENE_CALL_VERSION, keep_work: bool = False,
        tape_version: str | None = None) -> dict:
    """Run one recipe to completion; returns the search_run row. Refusals raise RecipeError.

    The run reads the tape version of its pack (card B-16); `tape_version`, when given, must be it."""
    root = Path(root).resolve()
    r, text = rcp.load(Path(recipe_src) if isinstance(recipe_src, Path) else recipe_src)
    res = resolve(root, r, release)
    pack_tv = pack.manifest_tape_version(res["manifest"])
    if tape_version is not None and tape_version != pack_tv:
        raise RecipeError(f"tape version {tape_version} is not the tape version {pack_tv} of the {release} pack: "
                          "a run reads the one tape version its pack's flanks refer to")
    r = res["recipe"]
    submitted = datetime.now(timezone.utc)
    rh = rcp.recipe_hash(r)
    run_id = run_id or new_run_id(rh, release, res["pack_id"], submitted)
    work = Path(work_base or default_work_base()) / run_id.hex()
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
        shards, models, sums = _stage(res, work, limits)
        seed_file = None
        if res["seed_bytes"] is not None:  # the bytes read and hashed at submit, not the path's later content
            seed_file = work / "seed" / Path(res["seed"]["path"]).name
            seed_file.parent.mkdir()
            seed_file.write_bytes(res["seed_bytes"])
        job = {"run_id": run_id.hex(), "root": str(root), "out": str(out), "recipe": r, "seed": res["seed"],
               "seed_file": str(seed_file) if seed_file else None,
               "positives": res["positives"], "target_kind": rcp.RUN_TARGETS[r["target"]],
               "shards": [str(p) for p in shards], "models": [str(p) for p in models],
               "meta": str(res["pack_dir"] / f"{r['target']}/{r['target']}_meta.parquet") if r["target"] != "unique" else None,
               "gene_call_version": gene_call_version, "tape_version": pack_tv, "cpus": limits.cpus,
               "max_output_bytes": int(r["max_output_gb"] * GB), "fault": fault,
               "negative_seed": int.from_bytes(run_id[:4], "big")}
        (work / "job.json").write_text(json.dumps(job, sort_keys=True))
        system = platform.system()
        argv, kind = worker_argv(work / "job.json", limits, system=system)
        env = worker_env(limits, kind, system, os.environ)
        with open(work / "stderr.txt", "w") as err:
            proc = subprocess.Popen(argv, stdout=subprocess.DEVNULL, stderr=err,
                                    preexec_fn=_preexec(limits, kind, system), env=env)
            rc, ru, timed_out = _wait(proc, limits.wall_s, cancelled)
        # only a scope's own cgroup is evidence: outside one, the cgroup is the runner's session's
        cgroup_oom = scope_oom_killed(work) if kind == "scope" and rc == -signal.SIGKILL and not timed_out else None
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
    if ru is not None:
        rss = ru.ru_maxrss if platform.system() == "Darwin" else ru.ru_maxrss * 1024
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
        "limits_json": json.dumps(asdict(limits), sort_keys=True), "error": error,
        "seed_json": json.dumps(res["seed"], sort_keys=True),
    }
    # a worker that died after a round has hits on disk, but no result: they are not the run's
    _publish(root, run_id, row, out, with_outputs=bool(result) and status in ("ok", "failed", "output_cap"))
    if not keep_work:
        shutil.rmtree(work, ignore_errors=True)
    return row


def _main_thread() -> bool:
    import threading

    return threading.current_thread() is threading.main_thread()


def _publish(root: Path, run_id: bytes, row: dict, out: Path, with_outputs: bool) -> None:
    final = run_dir(root, run_id)
    tmp = final.with_name(f".{final.name}.tmp-{os.getpid()}")
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist([row], schema=table("search_run").schema), tmp / "run.parquet")
    files = ["run.parquet"]
    if with_outputs:
        for name in RUN_FILES[1:]:
            if (out / name).exists():
                shutil.copyfile(out / name, tmp / name)
                files.append(name)
    if final.exists():
        raise FileExistsError(f"{final} exists: a run directory is never edited after it closes")
    os.replace(tmp, final)
    _index_update(root, {"run_id": run_id.hex(), "status": row["status"], "recipe_hash": row["recipe_hash"],
                         "ran_against": row["ran_against"], "files": files})


# --- async submit, status, cancel ----------------------------------------------


def _pending(root: Path) -> Path:
    return runs_dir(root) / "pending"


def submit(root: Path, recipe_path: Path, release: str, limits: Limits = Limits(),
           work_base: Path | None = None, fault: str | None = None) -> str:
    """Check the recipe now, run it detached, return the run_id at once."""
    root = Path(root).resolve()
    r, _ = rcp.load(Path(recipe_path))
    res = resolve(root, r, release)
    run_id = new_run_id(rcp.recipe_hash(res["recipe"]), release, res["pack_id"], datetime.now(timezone.utc))
    p = _pending(root)
    p.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-m", "biobox.cli", "search", "run", str(recipe_path), "--root", str(root),
           "--release", release, "--wait", "--run-id", run_id.hex(), "--cpus", str(limits.cpus),
           "--cpu-min", str(limits.cpu_min), "--wall-s", str(limits.wall_s)]
    if work_base is not None:
        cmd += ["--work", str(work_base)]
    if fault is not None:
        cmd += ["--fault", fault]
    # The runner's liveness is a flock on local disk (never under the data root: flock is unreliable on
    # NFS). It is taken here, before the pending record exists, and handed to the runner as an inherited
    # fd (the lock belongs to the open file description, so it lives exactly as long as the runner).
    work = (Path(work_base or default_work_base()) / run_id.hex()).resolve()
    work.mkdir(parents=True, exist_ok=True)
    lock_path = work / LOCK_NAME
    fd = os.open(lock_path, os.O_RDWR | os.O_CREAT, 0o644)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=open(p / f"{run_id.hex()}.stderr", "w"),
                                start_new_session=True, pass_fds=(fd,))
    finally:
        os.close(fd)  # the runner's copy keeps the lock
    info = p / f"{run_id.hex()}.json"
    tmp = info.with_name(f".{info.name}.tmp")
    tmp.write_text(json.dumps({"run_id": run_id.hex(), "pid": proc.pid, "lock": str(lock_path),
                               "status": "running"}))
    os.replace(tmp, info)
    return run_id.hex()


def status(root: Path, run_id: str) -> dict:
    d = run_dir(root, run_id)
    if (d / "run.parquet").exists():
        row = pq.read_table(d / "run.parquet").to_pylist()[0]
        return {"run_id": run_id, "status": row["status"], "error": row["error"], "rounds": row["rounds"],
                "final_round": row["final_round"], "stop_reason": row["stop_reason"]}
    p = _pending(root) / f"{run_id}.json"
    if p.exists():
        info = json.loads(p.read_text())
        if _runner_holds(info.get("lock")):
            return {"run_id": run_id, "status": "running"}
        if (d / "run.parquet").exists():  # closed between the first look and the lock probe
            return status(root, run_id)
        return {"run_id": run_id, "status": "failed", "error": "runner exited without closing the run"}
    raise KeyError(f"no run {run_id}")


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


def cancel(root: Path, run_id: str) -> None:
    """Ask a live runner to stop. A runner that is gone is not signalled: its pid may belong to
    another process by now."""
    info = json.loads((_pending(root) / f"{run_id}.json").read_text())
    if _runner_holds(info.get("lock")):
        os.kill(info["pid"], signal.SIGTERM)


def read_run(root: Path, run_id: str) -> dict:
    return pq.read_table(run_dir(root, run_id) / "run.parquet").to_pylist()[0]


def read_hits(root: Path, run_id: str, round_: int | None = None, final_only: bool = True) -> list[dict]:
    d = run_dir(root, run_id)
    if not (d / "hits.parquet").exists():
        return []
    rows = pq.read_table(d / "hits.parquet").to_pylist()
    if round_ is None and final_only:
        round_ = read_run(root, run_id)["final_round"]
        rows = [r for r in rows if r["curation_label"] == "keep"]
    if round_ is not None:
        rows = [r for r in rows if r["round"] == round_]
    return rows
