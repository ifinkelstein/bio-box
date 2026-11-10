"""The search process the runner launches under its limits: `python -m biobox.search.worker <job.json>`.

It reads only the staged shards and the data root, writes `hits.parquet`, `profiles.parquet`,
`neighborhoods.parquet` and `result.json` to the job's `out` directory, and checks the output
cap after every round. A MemoryError exits with EXIT_OOM so the runner records `oom`; any
other failure exits non-zero and the runner records `failed`. Neither writes hits.
"""

from __future__ import annotations

import json
import os
import sys
import time
import zlib
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from biobox.l3.tapes import tape_version
from biobox.registry.schemas import table
from biobox.search import engine, fastz, nbr

EXIT_OOM = 86
EXIT_ORPHAN = 87  # the runner died: nothing would record this run, so stop
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


def _families(root: Path) -> dict[bytes, set[str]]:
    out: dict[bytes, set[str]] = {}
    for f in sorted((root / "l1/domains").glob("library=*/stage=*/part-*.parquet")):
        for r in pq.read_table(f, columns=["protein_uid", "model_acc"]).to_pylist():
            out.setdefault(r["protein_uid"], set()).add(r["model_acc"])
    return out


def _has_domains(root: Path) -> bool:
    """Whether the release holds any l1_domains row at all (not whether these hits have one)."""
    return any(pq.ParquetFile(f).metadata.num_rows > 0
               for f in (root / "l1/domains").glob("library=*/stage=*/part-*.parquet"))


def _partial(root: Path, job: dict, target: engine.Target) -> set[bytes]:
    if target.kind == "protein":
        files = sorted((root / "l1/protein_occurrence").glob(f"gene_call_version={job['gene_call_version']}/stage=*/part-*.parquet"))
        return {r["protein_uid"] for f in files
                for r in pq.read_table(f, columns=["protein_uid", "dist_to_contig_end_bp"]).to_pylist()
                if r["dist_to_contig_end_bp"] == 0}
    meta = Path(job["meta"]) if job.get("meta") else None
    if meta is None or not meta.exists():
        return set()
    return {r["region_uid"] for r in pq.read_table(meta, columns=["region_uid", "partial"]).to_pylist() if r["partial"]}


def _fault(kind: str | None) -> None:
    if kind is None:
        return
    if kind == "oom":
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
    seed_families = engine.families_of(seed_seq or hmm.consensus.upper(), model_paths)
    if seed["kind"] in ("hmm", "hmm_file"):
        seed_families.add(engine.text(hmm.accession) or engine.text(hmm.name))
    settings = engine.Settings(
        inclusion_bits=it["inclusion_bits"], domain_bits=it["domain_bits"], max_rounds=it["max"],
        stop=tuple(it["stop"]), min_query_coverage=cur["min_query_coverage"], alarm_fraction=it["alarm_fraction"],
        alarm_min=it["alarm_min"], context_retention=it["context_retention"], partial_policy=recipe["partial"],
        cpus=job.get("cpus", 1))
    partial = _partial(root, job, target)
    fams = _families(root) if target.kind == "protein" and cur["foreign_domain_alarm"] else None
    policy = nbr.Policy.parse(**recipe["context"])
    # the one tape version the run reads: its pack's (card B-16), the gene-call producer's for an older job
    tv = job.get("tape_version") or tape_version([job["gene_call_version"]])
    have_tapes = any((root / "l3/tapes").glob(f"tape_version={tv}/stage=*/part-*.parquet"))
    warnings = []
    if fams is None and "drift_alarm" in it["stop"]:
        warnings.append("drift_alarm not evaluated: no stored library hits for this target kind")
    if fams is not None and not fams and it["max"] > 1 and "drift_alarm" in it["stop"] and not _has_domains(root):
        warnings.append("drift_alarm could not be evaluated: the release has no domain annotations (l1_domains is empty)")
    if not (target.kind == "protein" and have_tapes) and "context_lost" in it["stop"]:
        warnings.append("context_lost not evaluated: no tapes for this target kind")

    def families(uids):
        return {u: fams.get(u, set()) for u in uids}

    def context(uids):
        rows, _ = nbr.fetch(root, set(uids), policy, job["gene_call_version"], tape_version=tv)
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
        rows, counts = nbr.fetch(root, final, policy, job["gene_call_version"], tape_version=tv)
        nbr.write(rows, out / "neighborhoods.parquet")
        n_nbr = counts["rows"]
    result = {"status": status, "error": error, "rounds": last.k, "final_round": final_round,
              "stop_reason": last.checks.get("stop"), "controls": controls, "warnings": warnings,
              "target_size": target.size, "rounds_log": rounds_log, "n_final": len(final), "n_neighborhood_rows": n_nbr,
              "tape_version": tv}
    (out / "result.json").write_text(json.dumps(result, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1]))
    except MemoryError:
        os._exit(EXIT_OOM)
