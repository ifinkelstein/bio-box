"""Short-contig pilot on CAMI II marine sample 0 (gold-standard assembly; measurements file section 6c). Throwaway.

Adapted from the S1 phage driver (pilot_s1/pilot.py): ingest, a gene-calling rate probe on a seeded
2% contig sample, genome inputs, gene calling (every contig its own prokaryote genome, linear) on all
contigs or a seeded subset, IGRs, tapes, MMseqs2 clustering, tokens, token postings, contig_dim,
manifest and catalog. The 2026-11-09 run (store at /tmp/if2453/pilot_cami/work/root on the pod, code
at 2556ddd) used the pre-phase-B version of this file, which counted and sampled pair rows; card B-18
added the driver to the repo and moved it to the phase-B model.

    python pilot_cami.py SRC_FASTA WORKDIR RELEASE [--threads 8] [--workers 1] [--step-timeout 10800]
                         [--genecall-timeout 12600] [--full-limit-s 10800] [--subset-target-s 7200]

SRC_FASTA: the extracted gold-standard assembly (`anonymous_gsa.fasta.gz` of
`marmgCAMI2_sample_0_contigs.tar.gz`, sha256 e7639d22...1cb55a; on the pod under
/tmp/if2453/pilot_cami/src/x/). WORKDIR gets root/, mmseqs_tmp/ and result.json (rewritten after
every step; steps an earlier process finished without error are skipped on a rerun). The gene-call
time limits apply to the serial rate (`--workers 1`); with more workers the plan is conservative.
"""

import argparse
import json
import platform
import random
import resource
import signal
import statistics
import sys
import time
import traceback
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("src", type=Path)
ap.add_argument("workdir", type=Path)
ap.add_argument("release")
ap.add_argument("--threads", type=int, default=8, help="MMseqs2 and Arrow threads")
ap.add_argument("--workers", type=int, default=1, help="genomes called at once by genecall.run")
ap.add_argument("--step-timeout", type=int, default=10800)
ap.add_argument("--genecall-timeout", type=int, default=12600)
ap.add_argument("--probe-fraction", type=float, default=0.02)
ap.add_argument("--full-limit-s", type=float, default=10800)
ap.add_argument("--subset-target-s", type=float, default=7200)
ap.add_argument("--stage", default="SM")
ap.add_argument("--source", default="cami2")
a = ap.parse_args()

import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

pa.set_cpu_count(a.threads)
pa.set_io_thread_count(a.threads)

from biobox import catalog, io  # noqa: E402
from biobox.l0seq import SourceObject, ingest  # noqa: E402
from biobox.l0seq.inputs import write_genome_input  # noqa: E402
from biobox.l0seq.store import SeqStore  # noqa: E402
from biobox.l1 import igr  # noqa: E402
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, call_genome, run as call_genes  # noqa: E402
from biobox.l3 import tapes  # noqa: E402


def read(base, cols=None):
    """Every part of a table directory (card B-19 cuts tables into parts)."""
    return pa.concat_tables(pq.read_table(f, columns=cols, partitioning=None) for f in io.part_files(base))

W = a.workdir
W.mkdir(parents=True, exist_ok=True)
ROOT = W / "root"
OUT = W / "result.json"
KIB = 1024 if sys.platform == "linux" else 1024 * 1024
STAGE, SRC = a.stage, a.source

R: dict = {
    "host": platform.node(), "machine": platform.machine(), "python": sys.version.split()[0],
    "started_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "gene_call_version": GENE_CALL_VERSION, "release": a.release, "stage": STAGE, "source": SRC,
    "threads_mmseqs": a.threads, "arrow_cpu_count": pa.cpu_count(), "steps": {}, "failed": {},
}
if OUT.exists():
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
    raise StepTimeout("step passed its time limit")


signal.signal(signal.SIGALRM, _alarm)


def step(name, fn, threads=1, needs=(), timeout=None):
    if name in R["steps"]:
        print(f"{name}: done by an earlier process, skipped", flush=True)
        return R["steps"][name].get("result")
    missing = [n for n in needs if n not in R["steps"] or R["steps"][n].get("error")]
    if missing:
        R["steps"][name] = {"skipped": f"depends on failed/missing {missing}"}
        save()
        return None
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {name}: start", flush=True)
    t = time.perf_counter()
    c0 = resource.getrusage(resource.RUSAGE_SELF)
    signal.alarm(timeout or a.step_timeout)
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


R["source_file"] = str(a.src)
R["source_bytes"] = a.src.stat().st_size
save()

# --- 1 ingest ------------------------------------------------------------------
obj = SourceObject(a.src, SRC, "marmgCAMI2_sample_0_gsa", license="CAMI II (publisso FRL, doi:10.4126/FRL01-006425521)")
step("ingest", lambda: ingest([obj], ROOT, STAGE))


def pct(s, p):
    return s[min(len(s) - 1, int(p / 100 * len(s)))]


default: list = []
if not R["steps"]["ingest"].get("error"):
    occ = read(ROOT / f"l0seq/contig_occurrence/stage={STAGE}/source={SRC}", ["raw_digest", "contig_uid"])
    loc = read(ROOT / f"l0seq/seq_locator/stage={STAGE}", ["raw_digest", "contig_len", "shard"])
    ln_of = dict(zip(loc.column("raw_digest").to_pylist(), loc.column("contig_len").to_pylist()))
    rec_lens = [ln_of[rd] for rd in occ.column("raw_digest").to_pylist()]
    idx = read(ROOT / f"l0seq/contig_index/stage={STAGE}").to_pylist()
    default = sorted((r["contig_uid"], r["contig_len"]) for r in idx if r["is_default_variant"])
    s = sorted(rec_lens)
    R["records"] = {
        "records": len(s), "total_bp": sum(s), "min_len": s[0], "p5_len": pct(s, 5),
        "median_len": statistics.median(s), "p95_len": pct(s, 95), "max_len": s[-1],
        "n_lt_1kb": sum(1 for x in s if x < 1000), "n_1_10kb": sum(1 for x in s if 1000 <= x < 10000),
        "n_10_100kb": sum(1 for x in s if 10000 <= x < 100000), "n_ge_100kb": sum(1 for x in s if x >= 100000),
        "stored_sequences": loc.num_rows, "stored_bp": sum(ln_of.values()),
        "distinct_contig_uids": len(default), "distinct_contig_bp": sum(n for _, n in default),
        "shards": len(set(loc.column("shard").to_pylist())),
    }
    del occ, loc, ln_of, rec_lens, idx
    save()

# --- 2a gene-calling rate on a seeded 2% contig sample, and the plan -----------------
def rate_probe():
    store = SeqStore(ROOT)
    rng = random.Random(20261109)
    sample = rng.sample(default, max(1, round(a.probe_fraction * len(default))))
    t0 = time.perf_counter()
    genes = 0
    modes = defaultdict(lambda: [0, 0, 0.0])  # contigs, bp, s
    for cu, n in sample:
        t1 = time.perf_counter()
        calls = call_genome(store, GenomeInput((cu,), 11, kind="prokaryote"))
        dt = time.perf_counter() - t1
        genes += len(calls)
        m = modes["single" if n >= 100_000 else "meta"]
        m[0] += 1
        m[1] += n
        m[2] += dt
    wall = time.perf_counter() - t0
    bp = sum(n for _, n in sample)
    proj = wall * len(default) / len(sample)
    return {"contigs": len(sample), "bp": bp, "genes": genes, "wall_s": round(wall, 2),
            "mbp_per_s": round(bp / 1e6 / wall, 4), "contigs_per_s": round(len(sample) / wall, 2),
            "projected_full_s_by_count": round(proj, 1),
            "projected_full_s_by_bp": round(wall * R["records"]["distinct_contig_bp"] / bp, 1),
            "by_mode": {k: {"contigs": v[0], "bp": v[1], "s": round(v[2], 2)} for k, v in modes.items()}}


step("genecall_rate_probe", rate_probe, needs=("ingest",))

if "genecall_plan" not in R and R["steps"].get("genecall_rate_probe", {}).get("result"):
    pr = R["steps"]["genecall_rate_probe"]["result"]
    proj = max(pr["projected_full_s_by_count"], pr["projected_full_s_by_bp"])
    frac = 1.0 if proj < a.full_limit_s else round(a.subset_target_s / proj, 4)
    R["genecall_plan"] = {"projected_full_s": proj, "full_limit_s": a.full_limit_s, "fraction": frac,
                          "seed": 1109 if frac < 1 else None}
    save()

genomes: list = []
called: set = set()
if R.get("genecall_plan"):
    frac = R["genecall_plan"]["fraction"]
    pick = default if frac >= 1 else sorted(random.Random(1109).sample(default, round(frac * len(default))))
    genomes = [GenomeInput((cu,), 11, kind="prokaryote") for cu, _ in pick]
    called = {cu for cu, _ in pick}
    R["genecall_plan"]["contigs"] = len(pick)
    R["genecall_plan"]["bp"] = sum(n for _, n in pick)
    save()

# --- 2 gene calling (one serial process; library unchanged) ----------------------
# every default contig gets a genome_input row (contig_dim needs one per contig); genecall is given only
# the planned subset
step("genome_input", lambda: write_genome_input(ROOT, STAGE, [GenomeInput((cu,), 11, kind="prokaryote") for cu, _ in default]),
     needs=("ingest",))
step("genecall", lambda: call_genes(ROOT, STAGE, genomes, workers=a.workers), threads=a.workers,
     needs=("genome_input", "genecall_rate_probe"), timeout=a.genecall_timeout)

# --- 3 IGRs, tapes (all linear) ----------------------------------------------------
step("igr", lambda: igr.run(ROOT, STAGE, set()), needs=("genecall",))
step("tapes", lambda: tapes.run(ROOT, STAGE, set()), needs=("genecall",))


# --- 4 clustering ---------------------------------------------------------------
def do_cluster():
    from biobox.l2 import cluster

    tmp = W / "mmseqs_tmp"
    tmp.mkdir(exist_ok=True)
    return cluster.run(ROOT, a.release, GENE_CALL_VERSION, threads=a.threads, tmp=tmp)


step("cluster", do_cluster, threads=a.threads, needs=("genecall",))


# --- 5 tokens ---------------------------------------------------------------------
def do_tokens():
    from biobox.l3 import tokens

    return tokens.run(ROOT, a.release)


step("tokens", do_tokens, needs=("cluster", "tapes"))


# --- 6 token postings and contig_dim (card B-18: pairs are never materialized, AGENTS §3.5; the pair count and
# 2% pair sample of the 2026-11-09 run, which imported the removed biobox.l3.pairs, are replaced by the postings
# builder that the neighbour queries read) ---------------------------------------------------------------
def do_postings():
    from biobox.l3 import postings

    out = postings.run(ROOT, a.release)
    files = io.part_files(ROOT / "l3/token_postings", f"release={a.release}/tape_version={out['tape_version']}/token_bucket=*/part-*.parquet")
    nbytes = sum(f.stat().st_size for f in files)
    genes = (R["steps"].get("genecall", {}).get("result") or {}).get("genes") or 0
    return {**out, "bytes": nbytes, "bytes_per_row": round(nbytes / max(1, out["rows"]), 3),
            "rows_per_gene": round(out["rows"] / max(1, genes), 3)}


step("postings", do_postings, needs=("tokens",))


def do_dim():
    from biobox import dim

    return dim.run(ROOT, a.release, [STAGE])


step("contig_dim", do_dim, needs=("genome_input",))


# --- 7 manifest and catalog ----------------------------------------------------------
def do_manifest():
    man = catalog.build_manifest(ROOT, a.release)
    catalog.write_manifest(man, ROOT / f"manifests/release={a.release}.json")
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


step("manifest", do_manifest, needs=("ingest",))


def do_catalog():
    man = catalog.read_manifest(ROOT / f"manifests/release={a.release}.json")
    out = catalog.generate_catalog(man, ROOT / f"catalog/release={a.release}.duckdb")
    return {"path": str(out), "bytes": out.stat().st_size}


step("catalog", do_catalog, needs=("manifest",))

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
    gp = R.get("genecall_plan", {})
    R["derived"] = {"shard_bytes_per_bp": round(shard_b / bp, 4) if shard_b else None,
                    "shard_plus_index_bytes_per_bp": round(idx_b / bp, 4),
                    "genes_per_kb_called": round(genes / (gp["bp"] / 1e3), 4) if genes and gp.get("bp") else None}
    gc = R["steps"].get("genecall", {})
    if gc.get("result") and gp.get("bp"):
        R["derived"]["genecall_mbp_per_s"] = round(gp["bp"] / 1e6 / gc["wall_s"], 4)
        R["derived"]["genecall_contigs_per_s"] = round(gp["contigs"] / gc["wall_s"], 2)
        R["derived"]["genecall_proteins_per_s"] = round(genes / gc["wall_s"], 1)
        R["derived"]["distinct_protein_fraction"] = round(gc["result"]["distinct_proteins"] / max(1, genes), 4)
    ing = R["steps"].get("ingest", {})
    if ing.get("result"):
        R["derived"]["ingest_mbp_per_s"] = round(R["records"]["total_bp"] / 1e6 / ing["wall_s"], 3)
R["finished_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
save()
print("done", flush=True)
