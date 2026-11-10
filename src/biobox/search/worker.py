"""The search process the runner launches under its limits: `python -c BOOT <job.json>`.

It reads only the staged shards and, of the data root, the files of the job's release (`release_files`,
card B-28: its manifest's, else those of its `release_def` stages), writes `hits.parquet`, `profiles.parquet`,
`neighborhoods.parquet` and `result.json` to the job's `out` directory, and checks the output
cap after every round. A MemoryError exits with EXIT_OOM so the runner records `oom`; any
other failure exits non-zero and the runner records `failed`. Neither writes hits.

Threads (card W-04). `BOOT` sets the thread variables (`THREAD_VARS`) to the job's `cpus`
before it imports this module, because numpy's OpenBLAS starts its threads when pyarrow
imports numpy and Arrow sizes its pools from OMP_NUM_THREADS and ARROW_IO_THREADS (Arrow reads
no ARROW_NUM_THREADS); `main` then sets the Arrow CPU and I/O pools to `cpus` explicitly. The
worker opens no DuckDB connection; one added here must `SET threads` to `cpus`. pyhmmer's
searches get `cpus` through `engine.Settings`, and `engine.families_of` takes `threads=cpus`
(card W-09); no pyhmmer call sizes its pool from the host.

Exit (card W-04). `entry` runs `main` and always ends with `os._exit`: an uncaught exception
in the main thread or in any thread prints its traceback to stderr and exits EXIT_FAILED
(EXIT_OOM for a MemoryError) at once, so no interpreter shutdown waits on a non-daemon pool
thread that lost its peer. That wait is how a worker hung after "can't start new thread".
"""

from __future__ import annotations

import json
import os
import sys
import threading
import time
import traceback
import zlib
from pathlib import Path
from typing import TYPE_CHECKING

import pyarrow as pa
import pyarrow.parquet as pq

from biobox.l3.tapes import tape_version
from biobox.registry.schemas import table
from biobox.search import engine, fastz, nbr

if TYPE_CHECKING:
    from biobox.catalog import ReleaseFiles

# the variables every thread pool the worker loads reads its size from; run-step sets the first four too
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "POLARS_MAX_THREADS",
               "ARROW_IO_THREADS")
EXIT_OOM = 86
EXIT_ORPHAN = 87  # the runner died: nothing would record this run, so stop
EXIT_FAILED = 88  # an uncaught exception in any thread; its traceback is in stderr
ROUND_FAULTS = ("round2", "round2_thread")  # test hooks that fire inside round 2 of an iterative run
# test hook (card A-07): the `oom` fault starts this many seconds late, as on a host slower to reach its limit
OOM_DELAY_VAR = "BIOBOX_OOM_FAULT_DELAY_S"

# The worker's entry point (`python -c BOOT job.json`): thread variables first, then the import.
BOOT = ("import json, os, sys\n"
        "_cpus = str(max(1, int(json.load(open(sys.argv[1])).get('cpus', 1))))\n"
        f"for _v in {THREAD_VARS!r}:\n"
        "    os.environ[_v] = _cpus\n"
        "try:\n"
        "    from biobox.search import worker\n"
        "except MemoryError:\n"
        f"    os._exit({EXIT_OOM})\n"
        "except BaseException:\n"
        "    import traceback\n"
        "    traceback.print_exc()\n"
        "    sys.stderr.flush()\n"
        f"    os._exit({EXIT_FAILED})\n"
        "worker.entry(sys.argv[1])\n")


def _exit(code: int) -> None:
    """End the process now: no atexit, no joining of non-daemon threads (a dead pool's peers wait forever)."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.flush()
        except Exception:
            pass
    os._exit(code)


def _code_of(exc: BaseException) -> int:
    return EXIT_OOM if isinstance(exc, MemoryError) else EXIT_FAILED


def _thread_died(args: threading.ExceptHookArgs) -> None:
    """threading.excepthook: an exception that ends any thread ends the worker, with its traceback."""
    if args.exc_type is SystemExit:
        return
    name = args.thread.name if args.thread is not None else "?"
    print(f"Exception in worker thread {name}:", file=sys.stderr)
    traceback.print_exception(args.exc_type, args.exc_value, args.exc_traceback, file=sys.stderr)
    _exit(_code_of(args.exc_value))


def entry(job_path: str) -> None:
    """Run `main` and exit the process: its return code, EXIT_OOM for a MemoryError, EXIT_FAILED with
    the traceback in stderr for any other uncaught exception, in this thread or any other. Never returns."""
    threading.excepthook = _thread_died
    try:
        code = main(job_path)
    except SystemExit as e:
        code = e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
    except BaseException as e:  # noqa: BLE001 - every exception ends the worker, promptly
        traceback.print_exc()
        code = _code_of(e)
    _exit(code)


def cap_threads(cpus: int) -> None:
    """Size the pools of everything already imported to `cpus` (BOOT set the variables before import)."""
    n = max(1, int(cpus))
    for v in THREAD_VARS:
        os.environ[v] = str(n)
    pa.set_cpu_count(n)
    pa.set_io_thread_count(n)


def os_thread_count() -> int | None:
    """Threads of this process as the kernel counts them (native pools included); Linux only."""
    try:
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith("Threads:"):
                return int(line.split()[1])
    except OSError:
        pass
    return None


class _ThreadPeaks:
    """Thread counts sampled at the end of each round and of the run (a debug field of result.json)."""

    def __init__(self, cpus: int):
        self.cpus, self.python, self.os, self.samples = cpus, 0, None, 0

    def sample(self) -> None:
        self.samples += 1
        self.python = max(self.python, threading.active_count())
        n = os_thread_count()
        if n is not None:
            self.os = max(self.os or 0, n)

    def report(self) -> dict:
        return {"cpus": self.cpus, "python_peak": self.python, "os_peak": self.os, "samples": self.samples,
                "arrow_cpu_count": pa.cpu_count(), "arrow_io_thread_count": pa.io_thread_count(),
                "env": {v: os.environ.get(v) for v in THREAD_VARS}}
CGROUP_FILE = "worker.cgroup"  # beside job.json: the worker's cgroup v2 path, where the runner looks for an OOM kill


def cgroup_v2_path(proc_self_cgroup: str) -> str | None:
    """The unified (v2) hierarchy path in the text of /proc/self/cgroup, or None (no v2 line)."""
    for line in proc_self_cgroup.splitlines():
        if line.startswith("0::"):
            return line[3:]
    return None


def _record_cgroup(job_dir: Path) -> None:
    """Write this process's cgroup v2 path for the runner. Under a systemd scope it is the scope's
    own cgroup, whose memory.events counts the kernel's OOM kills (card C-05)."""
    try:
        path = cgroup_v2_path(Path("/proc/self/cgroup").read_text())
    except OSError:  # not Linux
        return
    if path:
        (job_dir / CGROUP_FILE).write_text(path + "\n")


def _zst(b: bytes | None) -> bytes | None:
    return None if b is None else fastz.zst(b)


def release_files(job: dict) -> ReleaseFiles:
    """The files of the job's release (card B-28): its manifest's, else those of its `release_def` stages.
    Every read of the data root goes through these, so a stage the release does not hold is never read."""
    from biobox.catalog import ReleaseFiles  # here, so importing the worker imports no DuckDB

    if not job.get("release"):
        raise ValueError("job has no release: the worker reads one release's files (card B-28)")
    return ReleaseFiles.of(Path(job["root"]), job["release"])


def _families(files: ReleaseFiles) -> dict[bytes, set[str]]:
    out: dict[bytes, set[str]] = {}
    for f in sorted(files.table("l1_domains"), key=str):
        for r in pq.read_table(f, columns=["protein_uid", "model_acc"]).to_pylist():
            out.setdefault(r["protein_uid"], set()).add(r["model_acc"])
    return out


def _has_domains(files: ReleaseFiles) -> bool:
    """Whether the release holds any l1_domains row at all (not whether these hits have one)."""
    return any(pq.ParquetFile(f).metadata.num_rows > 0 for f in files.table("l1_domains"))


def _partial(files: ReleaseFiles, job: dict, target: engine.Target) -> set[bytes]:
    if target.kind == "protein":
        paths = files.table("protein_occurrence", gene_call_version=job["gene_call_version"])
        return {r["protein_uid"] for f in paths
                for r in pq.read_table(f, columns=["protein_uid", "dist_to_contig_end_bp"]).to_pylist()
                if r["dist_to_contig_end_bp"] == 0}
    meta = Path(job["meta"]) if job.get("meta") else None
    if meta is None or not meta.exists():
        return set()
    return {r["region_uid"] for r in pq.read_table(meta, columns=["region_uid", "partial"]).to_pylist() if r["partial"]}


def _round_fault(kind: str | None, k: int) -> None:
    """Test hooks inside round 2: raise in the main thread, or in a thread the main thread then
    waits on forever, as a caller waits on a pool whose worker thread died."""
    if kind not in ROUND_FAULTS or k != 2:
        return
    if kind == "round2":
        raise RuntimeError("injected in round 2")
    answered = threading.Event()

    def pool_worker():
        raise RuntimeError("injected in a worker thread in round 2")

    threading.Thread(target=pool_worker, name="fault-pool-worker").start()
    answered.wait()  # never set


def _fault(kind: str | None) -> None:
    if kind is None or kind in ROUND_FAULTS:
        return
    if kind == "oom":
        if os.environ.get(OOM_DELAY_VAR):
            time.sleep(float(os.environ[OOM_DELAY_VAR]))
        # a real allocation past the limit the runner set: RLIMIT_AS (Linux without a scope) raises
        # MemoryError; the scope's MemoryMax gets the worker OOM-killed. Elsewhere, the error it would raise.
        limit = os.environ.get("BIOBOX_RLIMIT_AS") or os.environ.get("BIOBOX_MEMORY_MAX")
        if limit:
            _ = bytearray(int(limit) * 2)
        raise MemoryError("injected")
    if kind == "timeout":
        time.sleep(10**6)
    if kind == "cpu":
        x = 0
        while True:
            x = zlib.crc32(x.to_bytes(8, "little"))
    if kind == "crash":
        raise RuntimeError("injected crash")
    raise ValueError(f"unknown fault {kind}")


def _hit_rows(run_id: bytes, rnd: engine.Round, first: dict[bytes, int], target: engine.Target, model_len: int,
              partial: set[bytes]) -> list[dict]:
    rows = []
    for u in sorted(rnd.hits):
        h = rnd.hits[u]
        lab, reason = rnd.labels[u]
        for d in h["domains"]:
            rows.append({"run_id": run_id, "round": rnd.k, "target_kind": target.kind, "target_uid": u,
                         "bitscore": d["bitscore"], "evalue_zref": d["evalue_zref"], "target_size": target.size,
                         "env_start": d["env_start"], "env_end": d["env_end"], "ali_start": d["ali_start"],
                         "ali_end": d["ali_end"], "hmm_from": d["hmm_from"], "hmm_to": d["hmm_to"],
                         "model_len": model_len, "is_partial": u in partial,
                         "first_inclusion_round": first.get(u, rnd.k), "score_kind": "exact",
                         "curation_label": lab, "curation_reason": reason,
                         "seq_bitscore": h["seq_bitscore"], "seq_evalue_zref": h["seq_evalue_zref"]})
    return rows


def _write(rows: list[dict], name: str, path: Path) -> int:
    pq.write_table(pa.Table.from_pylist(rows, schema=table(name).schema), path, compression="zstd")
    return path.stat().st_size


def _degenerate(rows: list[dict]) -> list[str]:
    warn = []
    if len(rows) > 1:
        for col in ("bitscore", "evalue_zref", "env_start"):
            vals = {r[col] for r in rows}
            if len(vals) == 1:
                warn.append(f"degenerate column {col}: every final hit has {next(iter(vals))}")
    return warn


def _watch_parent() -> None:
    """Exit when the runner dies (the worker is reparented), so no orphan runs to the CPU cap."""
    import threading

    parent = os.getppid()

    def loop():
        while True:
            time.sleep(0.5)
            if os.getppid() != parent:
                os._exit(EXIT_ORPHAN)

    threading.Thread(target=loop, daemon=True).start()


def main(job_path: str) -> int:
    job = json.loads(Path(job_path).read_text())
    cpus = max(1, int(job.get("cpus", 1)))
    cap_threads(cpus)
    peaks = _ThreadPeaks(cpus)
    _watch_parent()
    _record_cgroup(Path(job_path).parent)
    (Path(job_path).parent / "worker.pid").write_text(str(os.getpid()))
    _fault(job.get("fault"))
    root, out = Path(job["root"]), Path(job["out"])
    out.mkdir(parents=True, exist_ok=True)
    recipe = job["recipe"]
    run_id = bytes.fromhex(job["run_id"])
    it, cur = recipe["iterations"], recipe["curation"]
    target = engine.load_target(job["target_kind"], (r for p in job["shards"] for r in fastz.read_records(Path(p))))
    seed = job["seed"]
    model_paths = [Path(p) for p in job["models"]]
    seed_seq = seed.get("sequence")
    seed_file = Path(job["seed_file"]) if job.get("seed_file") else None  # staged bytes, hashed at submit
    if seed["kind"] == "hmm":
        hmm = engine.find_model(model_paths, seed["value"])
        msa_text = None
    elif seed["kind"] == "hmm_file":  # used exactly as a library model is
        hmm = engine.read_hmm_file(seed_file.read_bytes(), seed.get("model"))
        msa_text = None
    elif seed["kind"] == "msa":
        msa = engine.read_msa(seed_file)
        msa_text = engine._msa_text(msa)
        hmm = engine.build_from_msa(f"{recipe['recipe']}_r1", msa)
    elif seed["kind"] == "sequences":
        hmm, msa = engine.build_from_sequences(f"{recipe['recipe']}_r1", engine.read_seed_fasta(seed_file.read_bytes()))
        msa_text = engine._msa_text(msa)
    else:
        hmm = engine.build_from_sequence(f"{recipe['recipe']}_r1", seed_seq)
        msa_text = f"# STOCKHOLM 1.0\n\nseed {seed_seq}\n//\n".encode()
    seed_families = engine.families_of(seed_seq or hmm.consensus.upper(), model_paths, threads=cpus)
    if seed["kind"] in ("hmm", "hmm_file"):
        seed_families.add(engine.text(hmm.accession) or engine.text(hmm.name))
    settings = engine.Settings(
        inclusion_bits=it["inclusion_bits"], domain_bits=it["domain_bits"], max_rounds=it["max"],
        stop=tuple(it["stop"]), min_query_coverage=cur["min_query_coverage"], alarm_fraction=it["alarm_fraction"],
        alarm_min=it["alarm_min"], context_retention=it["context_retention"], partial_policy=recipe["partial"],
        cpus=cpus)
    files = release_files(job)
    partial = _partial(files, job, target)
    fams = _families(files) if target.kind == "protein" and cur["foreign_domain_alarm"] else None
    policy = nbr.Policy.parse(**recipe["context"])
    # the one tape version the run reads: its pack's (card B-16), the gene-call producer's for an older job
    tv = job.get("tape_version") or tape_version([job["gene_call_version"]])
    have_tapes = bool(files.table("contig_tapes", tape_version=tv))
    warnings = []
    if fams is None and "drift_alarm" in it["stop"]:
        warnings.append("drift_alarm not evaluated: no stored library hits for this target kind")
    if fams is not None and not fams and it["max"] > 1 and "drift_alarm" in it["stop"] and not _has_domains(files):
        warnings.append("drift_alarm could not be evaluated: the release has no domain annotations (l1_domains is empty)")
    if not (target.kind == "protein" and have_tapes) and "context_lost" in it["stop"]:
        warnings.append("context_lost not evaluated: no tapes for this target kind")

    def families(uids):
        return {u: fams.get(u, set()) for u in uids}

    def context(uids):
        rows, _ = nbr.fetch(root, set(uids), policy, job["gene_call_version"], tape_version=tv, files=files)
        ctx: dict[bytes, set[str]] = {u: set() for u in uids}
        for r in rows:
            if r["offset"] != 0:
                ctx[r["anchor_protein_uid"]] |= fams.get(r["protein_uid"], set()) if fams else set()
        return ctx

    hit_rows, prof_rows, rounds_log, first = [], [], [], {}
    last = None
    cap = job["max_output_bytes"]
    status, error = "ok", None
    for rnd in engine.iterate(hmm, msa_text, seed_seq, target, settings, seed_families,
                              families if fams is not None else None,
                              context if (target.kind == "protein" and have_tapes and fams is not None) else None,
                              partial, recipe["recipe"]):
        _round_fault(job.get("fault"), rnd.k)
        if not rnd.checks.get("revert"):
            for u, (lab, _) in rnd.labels.items():
                if lab == "keep":
                    first.setdefault(u, rnd.k)
        hit_rows += _hit_rows(run_id, rnd, first, target, int(hmm.M), partial)
        prof_rows.append({"run_id": run_id, "round": rnd.k, "kind": "hmm", "profile_zst": _zst(rnd.hmm_text),
                          "msa_zst": _zst(rnd.msa_text), "n_seq": rnd.n_seq, "n_eff": rnd.n_eff})
        rounds_log.append({"round": rnd.k, "hits": len(rnd.hits),
                           "kept": sum(1 for lab, _ in rnd.labels.values() if lab == "keep"),
                           "new": len(rnd.new), **{k: v for k, v in rnd.checks.items() if k != "revert"}})
        last = rnd
        size = _write(hit_rows, "search_hit", out / "hits.parquet") + _write(prof_rows, "search_profile",
                                                                             out / "profiles.parquet")
        peaks.sample()
        if size > cap:
            status, error = "output_cap", f"outputs {size} bytes passed max_output_gb after round {rnd.k}"
            break
    final_round = last.k - 1 if last.checks.get("revert") else last.k
    final = {r["target_uid"] for r in hit_rows if r["round"] == final_round and r["curation_label"] == "keep"}
    controls = {"positive": [], "negative": None}
    for p in job["positives"]:
        found = bytes.fromhex(p["uid"]) in final
        controls["positive"].append({**p, "recovered": found})
    if status == "ok" and not all(p["recovered"] for p in controls["positive"]):
        status, error = "failed", "positive control not recovered"
    if recipe["controls"]["negative"] == "shuffled_seed" and status == "ok":
        neg = engine.negative_control(seed_seq or hmm.consensus.upper(), target, settings, job["negative_seed"])
        controls["negative"] = neg
        if not neg["pass"]:
            status, error = "failed", f"negative control: {neg['n_hits']} hits, expectation {neg['expectation']:.3g}"
    warnings += _degenerate([r for r in hit_rows if r["round"] == final_round and r["curation_label"] == "keep"])
    n_nbr = 0
    if status in ("ok", "failed") and "neighborhoods" in recipe["outputs"] and target.kind == "protein" and have_tapes:
        rows, counts = nbr.fetch(root, final, policy, job["gene_call_version"], tape_version=tv, files=files)
        nbr.write(rows, out / "neighborhoods.parquet")
        n_nbr = counts["rows"]
    peaks.sample()
    result = {"status": status, "error": error, "rounds": last.k, "final_round": final_round,
              "stop_reason": last.checks.get("stop"), "controls": controls, "warnings": warnings,
              "target_size": target.size, "rounds_log": rounds_log, "n_final": len(final), "n_neighborhood_rows": n_nbr,
              "tape_version": tv, "threads": peaks.report()}
    (out / "result.json").write_text(json.dumps(result, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":  # by hand only: pyarrow is imported above before any thread variable is set
    entry(sys.argv[1])
