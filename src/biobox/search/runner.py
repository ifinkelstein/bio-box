"""The `bio-box search` runner (operation O-RUN; SEARCH_PLAN SD-01, SD-15, SD-20; data_layer_and_tests.md §2B.3).

A run is a recorded job:

  1. The recipe is checked on paper (`recipe.load`) and against the release (`resolve`):
     the pack must verify (C-12), every seed, model ID and positive control must be in it,
     and nothing in the recipe may be beyond what this build runs. Refusals raise
     RecipeError and record nothing.
  2. The target's shards are copied to `<work>/<run_id>/` and checked against the pack
     manifest's sha256; their total must fit the /tmp quota.
  3. `biobox.search.worker` runs in a child process at nice 15, with RLIMIT_CPU at the CPU
     cap (300 CPU-min; the pod watchdog kills near 369) and, on Linux, RLIMIT_AS at the
     memory cap, under a wall-clock timeout.
  4. The exit decides the status: 0 -> the worker's own (`ok`, `failed` for a control,
     `output_cap`); EXIT_OOM or a SIGKILL the runner did not send -> `oom`; SIGXCPU or the
     wall timeout -> `timeout`; a cancel (SIGTERM) -> `stopped`; anything else -> `failed`.
     Only `ok`, `failed`-by-control and `output_cap` runs carry hits; the others write
     `run.parquet` alone, never an empty hits file (S-05).
  5. The run directory is written next to its final place and renamed in
     (`search/run=<run_id>/`), then `search/run_index.json` is rewritten by atomic rename
     under a lock.
"""

from __future__ import annotations

import fcntl
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

from biobox import ids, provenance
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.registry.schemas import table
from biobox.search import engine, fastz, pack, recipe as rcp
from biobox.search.recipe import RecipeError
from biobox.search.worker import EXIT_OOM

GB = 10**9
# the worker's entry point: a MemoryError while importing under RLIMIT_AS is an OOM too
_BOOT = ("import os, sys\n"
         "try:\n"
         "    from biobox.search import worker\n"
         "    sys.exit(worker.main(sys.argv[1]))\n"
         "except MemoryError:\n"
         f"    os._exit({EXIT_OOM})\n")
RUN_FILES = ("run.parquet", "hits.parquet", "profiles.parquet", "neighborhoods.parquet")


@dataclass(frozen=True)
class Limits:
    cpu_min: float = 300.0       # per search process (SD-15)
    wall_s: float = 6 * 3600.0
    mem_gb: float | None = None  # RLIMIT_AS, Linux only
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
    elif seed["kind"] == "msa":
        if not Path(seed["value"]).exists():
            raise RecipeError(f"seed MSA {seed['value']} not found")
    positives = []
    for p in r["controls"]["positive"]:
        uid = ids.protein_uid(p["value"]).hex() if p["kind"] == "sequence" else p["value"]
        if uid not in target_ids:
            raise RecipeError(f"positive control {p.get('label', p['value'][:20])} is not in target {tname}: "
                              "it could never be recovered")
        positives.append({"kind": p["kind"], "uid": uid, "label": p.get("label")})
    return {"manifest": manifest, "pack_id": pid, "shards": shards, "models": models, "seed": seed,
            "positives": positives, "pack_dir": pdir}


def _set_limits(limits: Limits):
    def fn():
        os.setsid()
        os.nice(limits.nice)
        cpu = max(1, int(limits.cpu_min * 60))
        resource.setrlimit(resource.RLIMIT_CPU, (cpu, cpu + 5))
        if limits.mem_gb is not None and platform.system() == "Linux":
            b = int(limits.mem_gb * GB)
            resource.setrlimit(resource.RLIMIT_AS, (b, b))
    return fn


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


def _status(rc: int, timed_out: bool, cancelled: bool, out: Path, stderr: str) -> tuple[str, str | None, dict]:
    if timed_out:
        return "timeout", "wall-clock timeout", {}
    if cancelled:
        return "stopped", "cancelled", {}
    if rc == 0 and (out / "result.json").exists():
        res = json.loads((out / "result.json").read_text())
        return res["status"], res["error"], res
    if rc == EXIT_OOM or rc == -signal.SIGKILL:
        return "oom", "out of memory" if rc == EXIT_OOM else "killed (SIGKILL; OOM killer)", {}
    if rc == -signal.SIGXCPU:
        return "timeout", "CPU limit reached (RLIMIT_CPU)", {}
    tail = stderr.strip().splitlines()[-1:] if stderr.strip() else []
    return "failed", f"worker exited {rc}" + (f": {tail[0]}" if tail else ""), {}


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
        fault: str | None = None, gene_call_version: str = GENE_CALL_VERSION, keep_work: bool = False) -> dict:
    """Run one recipe to completion; returns the search_run row. Refusals raise RecipeError."""
    root = Path(root).resolve()
    r, text = rcp.load(Path(recipe_src) if isinstance(recipe_src, Path) else recipe_src)
    res = resolve(root, r, release)
    submitted = datetime.now(timezone.utc)
    rh = rcp.recipe_hash(r)
    run_id = run_id or new_run_id(rh, release, res["pack_id"], submitted)
    work = Path(work_base or default_work_base()) / run_id.hex()
    if work.exists():
        shutil.rmtree(work)
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
        job = {"run_id": run_id.hex(), "root": str(root), "out": str(out), "recipe": r, "seed": res["seed"],
               "positives": res["positives"], "target_kind": rcp.RUN_TARGETS[r["target"]],
               "shards": [str(p) for p in shards], "models": [str(p) for p in models],
               "meta": str(res["pack_dir"] / f"{r['target']}/{r['target']}_meta.parquet") if r["target"] != "unique" else None,
               "gene_call_version": gene_call_version, "cpus": limits.cpus,
               "max_output_bytes": int(r["max_output_gb"] * GB), "fault": fault,
               "negative_seed": int.from_bytes(run_id[:4], "big")}
        (work / "job.json").write_text(json.dumps(job, sort_keys=True))
        env = dict(os.environ)
        if limits.mem_gb is not None and platform.system() == "Linux":
            env["BIOBOX_RLIMIT_AS"] = str(int(limits.mem_gb * GB))
        with open(work / "stderr.txt", "w") as err:
            proc = subprocess.Popen([sys.executable, "-c", _BOOT, str(work / "job.json")],
                                    stdout=subprocess.DEVNULL, stderr=err, preexec_fn=_set_limits(limits), env=env)
            rc, ru, timed_out = _wait(proc, limits.wall_s, cancelled)
        status, error, result = _status(rc, timed_out, bool(cancelled), out, (work / "stderr.txt").read_text())
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
        "query_hash": blake3(json.dumps(res["seed"], sort_keys=True).encode()).hexdigest(),
        "submitted_at": submitted, "wall_s": wall, "peak_rss_bytes": rss, "cpu_min_max_pid": cpu_min,
        "host": socket.gethostname(), "exit_code": rc, "status": status, "partial_policy": r["partial"],
        "max_output_gb": r["max_output_gb"], "target": r["target"], "target_size": result.get("target_size"),
        "rounds": result.get("rounds"), "final_round": result.get("final_round"),
        "rounds_json": json.dumps(result["rounds_log"], default=str) if result.get("rounds_log") else None,
        "stop_reason": result.get("stop_reason"),
        "controls_json": json.dumps(result["controls"], sort_keys=True, default=str) if result.get("controls") else None,
        "warnings_json": json.dumps(result.get("warnings", [])) if result else None,
        "limits_json": json.dumps(asdict(limits), sort_keys=True), "error": error,
    }
    _publish(root, run_id, row, out, status)
    if not keep_work:
        shutil.rmtree(work, ignore_errors=True)
    return row


def _main_thread() -> bool:
    import threading

    return threading.current_thread() is threading.main_thread()


def _publish(root: Path, run_id: bytes, row: dict, out: Path, status: str) -> None:
    final = run_dir(root, run_id)
    tmp = final.with_name(f".{final.name}.tmp-{os.getpid()}")
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist([row], schema=table("search_run").schema), tmp / "run.parquet")
    files = ["run.parquet"]
    if status in ("ok", "failed", "output_cap"):
        for name in RUN_FILES[1:]:
            if (out / name).exists():
                shutil.copyfile(out / name, tmp / name)
                files.append(name)
    if final.exists():
        raise FileExistsError(f"{final} exists: a run directory is never edited after it closes")
    os.replace(tmp, final)
    _index_update(root, {"run_id": run_id.hex(), "status": status, "recipe_hash": row["recipe_hash"],
                         "ran_against": row["ran_against"], "files": files})


# --- async submit, status, cancel ----------------------------------------------


def _pending(root: Path) -> Path:
    return runs_dir(root) / "pending"


def submit(root: Path, recipe_path: Path, release: str, limits: Limits = Limits(),
           work_base: Path | None = None) -> str:
    """Check the recipe now, run it detached, return the run_id at once."""
    root = Path(root).resolve()
    r, _ = rcp.load(Path(recipe_path))
    res = resolve(root, r, release)
    run_id = new_run_id(rcp.recipe_hash(r), release, res["pack_id"], datetime.now(timezone.utc))
    p = _pending(root)
    p.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-m", "biobox.cli", "search", "run", str(recipe_path), "--root", str(root),
           "--release", release, "--wait", "--run-id", run_id.hex(), "--cpus", str(limits.cpus),
           "--cpu-min", str(limits.cpu_min), "--wall-s", str(limits.wall_s)]
    if work_base is not None:
        cmd += ["--work", str(work_base)]
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=open(p / f"{run_id.hex()}.stderr", "w"),
                            start_new_session=True)
    (p / f"{run_id.hex()}.json").write_text(json.dumps({"run_id": run_id.hex(), "pid": proc.pid, "status": "running"}))
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
        try:
            os.kill(info["pid"], 0)
            return {"run_id": run_id, "status": "running"}
        except ProcessLookupError:
            return {"run_id": run_id, "status": "failed", "error": "runner exited without closing the run"}
    raise KeyError(f"no run {run_id}")


def cancel(root: Path, run_id: str) -> None:
    info = json.loads((_pending(root) / f"{run_id}.json").read_text())
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
