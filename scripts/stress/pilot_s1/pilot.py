"""S1 rate pilot on INPHARED: ingest, genome inputs, gene calling (p2, virus, linear), IGRs, tapes,
MMseqs2 clustering, tokens, token postings, contig_dim, manifest. Throwaway.

    python pilot.py SRC_FASTA WORKDIR RELEASE [--limit N] [--threads 16] [--workers 1] [--step-timeout 10800]
                    [--overhead-n 200] [--resume]

WORKDIR gets root/, mmseqs_tmp/, manifest.json and result.json (rewritten after every step so a crash
leaves partial results). With --resume (what `pilot_resume.py` passes), steps that an earlier process
finished without error, as recorded in WORKDIR/result.json, are skipped. Card B-18 replaced the
pre-phase-B steps `pairs_count` and `pairs_sample` (which imported the removed `biobox.l3.pairs` and
wrote the removed `slot_pairs` table) by the postings builder.
"""

import argparse
import gzip
import json
import platform
import resource
import signal
import statistics
import sys
import time
import traceback
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

ap = argparse.ArgumentParser()
ap.add_argument("src", type=Path)
ap.add_argument("workdir", type=Path)
ap.add_argument("release")
ap.add_argument("--limit", type=int, default=0)
ap.add_argument("--threads", type=int, default=16, help="MMseqs2 threads")
ap.add_argument("--workers", type=int, default=1, help="genomes called at once by genecall.run")
ap.add_argument("--step-timeout", type=int, default=10800)
ap.add_argument("--overhead-n", type=int, default=200)
ap.add_argument("--resume", action="store_true", help="skip steps an earlier process finished without error")
a = ap.parse_args()

import pyrodigal_gv  # noqa: E402

from biobox import catalog, io  # noqa: E402
from biobox.l0seq import SourceObject, ingest  # noqa: E402
from biobox.l0seq.inputs import write_genome_input  # noqa: E402
from biobox.l0seq.store import SeqStore  # noqa: E402
from biobox.l1 import igr  # noqa: E402
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, call_genome, run as call_genes  # noqa: E402
from biobox.l3 import tapes  # noqa: E402

W = a.workdir
W.mkdir(parents=True, exist_ok=True)
ROOT = W / "root"
OUT = W / "result.json"
KIB = 1024 if sys.platform == "linux" else 1024 * 1024

R: dict = {
    "host": platform.node(), "machine": platform.machine(), "python": sys.version.split()[0],
    "started_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "gene_call_version": GENE_CALL_VERSION, "release": a.release, "limit": a.limit,
    "threads_mmseqs": a.threads, "steps": {}, "failed": {},
}
# resume: keep the steps an earlier process finished without error and skip them here
if a.resume and OUT.exists():
    prev = json.loads(OUT.read_text())
    prev["steps"] = {k: v for k, v in prev.get("steps", {}).items() if not v.get("error") and "skipped" not in v}
    prev["failed"] = {}
    prev["resumed_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    R = prev


def rss_now_mb() -> float:
    with open("/proc/self/statm") as fh:
        return round(int(fh.read().split()[1]) * resource.getpagesize() / 2**20, 1)


def peak_mb() -> float:
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / KIB, 1)


def child_peak_mb() -> float:
    return round(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / KIB, 1)


def save() -> None:
    tmp = OUT.with_suffix(".tmp")
    tmp.write_text(json.dumps(R, indent=1, default=str))
    tmp.replace(OUT)


class StepTimeout(Exception):
    pass


def _alarm(signum, frame):
    raise StepTimeout(f"step passed {a.step_timeout} s")


signal.signal(signal.SIGALRM, _alarm)


def step(name, fn, threads=1, needs=()):
    if a.resume and name in R["steps"]:
        print(f"{name}: done by an earlier process, skipped", flush=True)
        return R["steps"][name].get("result")
    missing = [n for n in needs if n not in R["steps"] or R["steps"][n].get("error")]
    if missing:
        R["steps"][name] = {"skipped": f"depends on failed/missing {missing}"}
        save()
        return None
    t = time.perf_counter()
    c0 = resource.getrusage(resource.RUSAGE_SELF)
    signal.alarm(a.step_timeout)
    try:
        out = fn()
        err = None
    except BaseException as e:  # noqa: BLE001 - failures at scale are findings
        out, err = None, f"{type(e).__name__}: {e}"
        R["failed"][name] = {"error": err, "traceback": traceback.format_exc()}
    finally:
        signal.alarm(0)
    c1 = resource.getrusage(resource.RUSAGE_SELF)
    R["steps"][name] = {
        "wall_s": round(time.perf_counter() - t, 2),
        "cpu_s": round((c1.ru_utime + c1.ru_stime) - (c0.ru_utime + c0.ru_stime), 2),
        "threads": threads, "peak_rss_mb": peak_mb(), "rss_now_mb": rss_now_mb(),
        "child_peak_rss_mb": child_peak_mb(), "result": out, "error": err,
    }
    save()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {name}: {R['steps'][name]['wall_s']} s err={err}", flush=True)
    return out


# --- source (dry run: first N records into a plain FASTA) ---------------------
src = a.src
if a.limit:
    src = W / f"subset_{a.limit}.fa"
    n = 0
    with gzip.open(a.src, "rb") as fi, open(src, "wb") as fo:
        for ln in fi:
            if ln.startswith(b">"):
                n += 1
                if n > a.limit:
                    break
            fo.write(ln)
R["source_file"] = str(src)
R["source_bytes"] = src.stat().st_size

# --- 1 ingest ------------------------------------------------------------------
obj = SourceObject(src, "inphared", "14Apr2025", license="AGPL-3.0 (repo); data from INSDC")
step("ingest", lambda: ingest([obj], ROOT, "S1"))

# record statistics from the store's own tables
genomes: list = []
if not R["steps"]["ingest"].get("error"):
    def read(base, cols=None):  # every part of a table directory (card B-19 cuts tables into parts)
        return pa.concat_tables(pq.read_table(f, columns=cols, partitioning=None) for f in io.part_files(base))

    occ = read(ROOT / "l0seq/contig_occurrence/stage=S1/source=inphared", ["raw_digest", "contig_uid"])
    loc = read(ROOT / "l0seq/seq_locator/stage=S1", ["raw_digest", "contig_len", "shard"])
    ln_of = dict(zip(loc.column("raw_digest").to_pylist(), loc.column("contig_len").to_pylist()))
    rec_lens = [ln_of[rd] for rd in occ.column("raw_digest").to_pylist()]
    idx = read(ROOT / "l0seq/contig_index/stage=S1").to_pylist()
    default = sorted((r["contig_uid"], r["contig_len"]) for r in idx if r["is_default_variant"])
    R["records"] = {
        "records": len(rec_lens), "total_bp": sum(rec_lens),
        "median_len": statistics.median(rec_lens), "max_len": max(rec_lens), "min_len": min(rec_lens),
        "stored_sequences": loc.num_rows, "stored_bp": sum(ln_of.values()),
        "distinct_contig_uids": len(default), "distinct_contig_bp": sum(n for _, n in default),
        "shards": len(set(loc.column("shard").to_pylist())),
    }
    genomes = [GenomeInput((cu,), 11, kind="virus") for cu, _ in default]
    save()

# --- 2a finder-construction overhead on the first N genomes (sorted uid order ~ random) ----
def overhead():
    store = SeqStore(ROOT)
    sample = genomes[: a.overhead_n]
    seqs = [store.get(g.contig_uids[0]).upper().decode("ascii") for g in sample]
    t_new = t_find = t_tr = 0.0
    per_genome = []
    for s in seqs:
        t0 = time.perf_counter()
        f = pyrodigal_gv.ViralGeneFinder(meta=True)
        t1 = time.perf_counter()
        genes = f.find_genes(s)
        t2 = time.perf_counter()
        for g in genes:
            g.translate(include_stop=False)
        t3 = time.perf_counter()
        t_new += t1 - t0
        t_find += t2 - t1
        t_tr += t3 - t2
        per_genome.append([(g.begin, g.end, g.strand, g.translation_table) for g in genes])
    shared = pyrodigal_gv.ViralGeneFinder(meta=True)
    t0 = time.perf_counter()
    same = 0
    for s, ref in zip(seqs, per_genome):
        genes = shared.find_genes(s)
        same += [(g.begin, g.end, g.strand, g.translation_table) for g in genes] == ref
    t_shared = time.perf_counter() - t0
    t0 = time.perf_counter()
    n_calls = sum(len(call_genome(store, g)) for g in sample)
    t_call = time.perf_counter() - t0
    bp = sum(len(s) for s in seqs)
    return {"genomes": len(seqs), "bp": bp, "genes": n_calls,
            "finder_construct_s": round(t_new, 3), "finder_construct_ms_each": round(1000 * t_new / max(1, len(seqs)), 3),
            "find_genes_s": round(t_find, 3), "translate_s": round(t_tr, 3),
            "construct_share_of_construct_plus_find": round(t_new / max(1e-9, t_new + t_find), 4),
            "shared_finder_find_genes_s": round(t_shared, 3), "shared_finder_identical_genomes": same,
            "call_genome_total_s": round(t_call, 3), "call_genome_mbp_per_s": round(bp / 1e6 / max(1e-9, t_call), 3)}


step("genecall_overhead_probe", overhead, needs=("ingest",))

# --- 2 gene calling --------------------------------------------------------------
step("genome_input", lambda: write_genome_input(ROOT, "S1", genomes), needs=("ingest",))
step("genecall", lambda: call_genes(ROOT, "S1", genomes, workers=a.workers), threads=a.workers, needs=("genome_input",))

# --- 3 IGRs, tapes (all linear: INPHARED topology unknown) ---------------------------
step("igr", lambda: igr.run(ROOT, "S1", set()), needs=("genecall",))
step("tapes", lambda: tapes.run(ROOT, "S1", set()), needs=("genecall",))

# --- 4 clustering ---------------------------------------------------------------
def do_cluster():
    from biobox.l2 import cluster

    tmp = W / "mmseqs_tmp"
    tmp.mkdir(exist_ok=True)
    return cluster.run(ROOT, a.release, GENE_CALL_VERSION, ["S1"], threads=a.threads, tmp=tmp)


step("cluster", do_cluster, threads=a.threads, needs=("genecall",))

# --- 5 tokens ---------------------------------------------------------------------
def do_tokens():
    from biobox.l3 import tokens

    return tokens.run(ROOT, a.release, ["S1"])  # the release's stages (card B-25)


step("tokens", do_tokens, needs=("cluster", "tapes"))

# --- 6 token postings (card B-18: pairs are never materialized, AGENTS §3.5; the pair count and the 2% pair
# sample of the 2026-11-09 run are replaced by the postings builder, which neighbour queries read) ------
def do_postings():
    from biobox.l3 import postings

    out = postings.run(ROOT, a.release, ["S1"])
    files = io.part_files(ROOT / "l3/token_postings", f"release={a.release}/tape_version={out['tape_version']}/token_bucket=*/part-*.parquet")
    nbytes = sum(f.stat().st_size for f in files)
    genes = (R["steps"].get("genecall", {}).get("result") or {}).get("genes") or 0
    return {**out, "bytes": nbytes, "bytes_per_row": round(nbytes / max(1, out["rows"]), 3),
            "rows_per_gene": round(out["rows"] / max(1, genes), 3)}


step("postings", do_postings, needs=("tokens",))


def do_dim():
    from biobox import dim

    return dim.run(ROOT, a.release, ["S1"])


step("contig_dim", do_dim, needs=("ingest",))

# --- 7 manifest -------------------------------------------------------------------
def do_manifest():
    man = catalog.build_manifest(ROOT, a.release, ["S1"])
    catalog.write_manifest(man, W / "manifest.json")
    per = defaultdict(lambda: {"files": 0, "rows": 0, "bytes": 0, "row_groups": 0})
    for f in man["files"]:
        p = per[f["table"]]
        p["files"] += 1
        p["rows"] += f["rows"]
        p["bytes"] += f["bytes"]
        p["row_groups"] += f["row_groups"]
    for f in man["artifacts"]:
        suffix = ".fai" if f["path"].endswith(".fai") else ".gzi" if f["path"].endswith(".gzi") else ""
        p = per[f"artifact:{f['artifact']}{suffix}"]
        p["files"] += 1
        p["bytes"] += f["bytes"]
    return {"manifest_hash": catalog.manifest_hash(man), "tables": dict(sorted(per.items()))}


step("manifest", do_manifest)

# derived per-table ratios
m = R["steps"].get("manifest", {}).get("result")
if m and R.get("records"):
    genes = (R["steps"].get("genecall", {}).get("result") or {}).get("genes") or 0
    bp = R["records"]["stored_bp"]
    for name, p in m["tables"].items():
        p["bytes_per_row"] = round(p["bytes"] / p["rows"], 3) if p["rows"] else None
        p["bytes_per_gene"] = round(p["bytes"] / genes, 3) if genes else None
        p["bytes_per_mbp"] = round(p["bytes"] / (bp / 1e6), 1)
    shard_b = m["tables"].get("artifact:seq_shards", {}).get("bytes")
    idx_b = sum(p["bytes"] for k, p in m["tables"].items() if k.startswith("artifact:"))
    R["derived"] = {"shard_bytes_per_bp": round(shard_b / bp, 4) if shard_b else None,
                    "shard_plus_index_bytes_per_bp": round(idx_b / bp, 4),
                    "genes_per_kb": round(genes / (R["records"]["distinct_contig_bp"] / 1e3), 4) if genes else None}
    gc = R["steps"].get("genecall", {})
    if gc.get("result"):
        R["derived"]["genecall_mbp_per_s"] = round(R["records"]["distinct_contig_bp"] / 1e6 / gc["wall_s"], 3)
        R["derived"]["genecall_proteins_per_s"] = round(genes / gc["wall_s"], 1)
R["finished_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
save()
print("done", flush=True)
