"""Card B-58 (amended 2026-11-21): what does the ingest of stage C5 spend its CPU-minutes on?

Arms, each a `build-stage --depth stored` (the `inputs` and `ingest` steps, as card C-47's driver ran them) of a
sample of C5's own inputs in a fresh scratch root under /tmp/if2453/b58/roots (never the data root; nothing under
/stor), as a detached `pod_run.sh` unit of 2 threads at nice 19. The figure of an arm is `cpu_min` of the ingest step's
run record (`runs/stage=<S>/ingest.json`), the quantity the C-47 rate guard judged.

  gz       the sample as `.fna.gz` objects (as C5 read them)            one object per genome
  plain    the same genomes as plain `.fna` copies                      one object per genome
  few      the same bases as one plain object (the copies concatenated)
  view     `gz` in a view of the data root's earlier stages (their contig_index files are read by the anti-join)
  held     `gz` as stage A of a root, then the same objects as stage B of that root: every contig of B is held
  spx      4 of SPX1's files (gzip) with every contig kept
  spx5000  the same 4 files with --min-contig-len 5000 (contigs under 5 kb are read and dropped)
  scale    every 3rd object of C5 as `.fna.gz` (`scale-inputs` writes it; run with --arms scale)

The sample: every 10th object of C5's objects.json (index 0, 10, 20, ...), one genome per object listing the
contigs (ids.contig_uid of each record) not already listed by an earlier sampled genome. SPX1: the 4 files
in the middle of the size order (ranks n/2-2 .. n/2+1).

  python3 b58_ingest_arms.py inputs              (inside the image, PYTHONPATH=<tree>/src: writes the sample)
  python3 b58_ingest_arms.py run [--reps 3]      (on the host: runs arm after arm, interleaved, logs to results.jsonl)
  python3 b58_ingest_arms.py report              (a table from results.jsonl)

Interleaved: rep 1 of every arm, then rep 2 of every arm, and so on. Before each arm /tmp/if2453/r10/root/STOP is
checked (an arm does not start while it exists); the load average and the active units are logged at the start and
end of every arm. A root of an arm is removed once its record is read; the logs stay.

Card B-61 (2026-11-21) adds arms on 12 further files of SPX1, in 6 groups of 2 (`groups-inputs` writes them; run with
B58_DIR=/tmp/if2453/b61 ARMS_TAG=b61 TREE=<the B-61 tree> and --arms spxg1,spxg1m,...,spxg6,spxg6m):

  spxg<k>   group k with every contig kept (as `spx`), then the same objects as stage B of that root: every contig
            held (as `held`), so stage B reads and hashes every record and stores nothing
  spxg<k>m  group k with --min-contig-len 5000 (as `spx5000`)

The 12 files: SPX1's objects in the order of (file size, path), less card B-58's 4 (`spx`), are 300; the files at
ranks floor((i + 0.5) x 300 / 12), i = 0..11, of that order, an even spread over the size range by rank. Group k
is the files i = 2k-2 and 2k-1 (two neighbours in size). The report (`inputs/report.json`) counts each file's FASTA
bytes, records and bases in all, at or over 5,000 bp and under it, and the FASTA bytes of the records under it, by
a loop over the decompressed lines written here (not the code's counter of card B-61).
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

B = Path(os.environ.get("B58_DIR", "/tmp/if2453/b58"))
IN = B / "inputs"
ROOTS = B / "roots"
LOGS = B / "logs"
RESULTS = B / "results.jsonl"
DATA_ROOT = Path("/tmp/if2453/r10/root")
C5 = Path("/tmp/if2453/s2/inputs/C5")
SPX = Path("/tmp/if2453/s2/inputs/SPX1")
TREE = os.environ.get("TREE", "/tmp/if2453/ci/B-58")
IMAGE = "/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif"
DIGEST = "sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a"
EARLIER = "S1,P1,P2,P3,P4,P5,P6,H1,PL1"
ARMS = ["gz", "plain", "few", "view", "held", "spx", "spx5000"]  # `scale` is asked for by name: --arms scale
TAG = os.environ.get("ARMS_TAG", "b58")  # the units' prefix (card B-61: b61)
N_GROUPS, MIN_LEN = 6, 5000  # card B-61's groups of SPX1 files and the least length of its `m` arms
THREADS = 2
PROFILE = """max_threads_day: 40
max_threads_night: 56
mem_scope: 32G
tmp: {tmp}
durable: {b}/durable
bulk: {b}/bulk
nice_monolithic: 19
nice_fanout: 15
day_start_hour: 8
day_end_hour: 20
cpu_min_per_process: 300
cluster_cpu_min_max: 20000
small_file_floor_bytes: 1000
require_mem_cap: true
"""


# --- the sample (inside the image) -----------------------------------------------------------------------------


def _genome(uids: list[bytes], seen: set) -> dict | None:
    kept = []
    for u in uids:
        if u not in seen:
            seen.add(u)
            kept.append(u.hex())
    if not kept:
        return None
    return {"contig_uids": kept, "genetic_code": 11, "kind": "prokaryote", "circular": [],
            "assembly_quality": "complete"}


def make_inputs() -> None:
    import gzip

    from biobox import fasta, ids

    report: dict = {}
    objs = json.load(open(C5 / "objects.json"))
    sample = objs[::10]
    for d in ("c5gz", "c5plain", "c5few", "spx", "spx5000"):
        (IN / d).mkdir(parents=True, exist_ok=True)
    (IN / "c5plain" / "files").mkdir(exist_ok=True)
    seen: set = set()
    genomes, plain_objs = [], []
    n_rec = n_bases = 0
    few = IN / "c5few" / "all.fna"
    with open(few, "wb") as cat:
        for o in sample:
            p = Path(o["path"])
            uids = []
            for r in fasta.iter_records(p):
                uids.append(ids.contig_uid(r.seq))
                n_rec += 1
                n_bases += len(r.seq)
            g = _genome(uids, seen)
            if g:
                genomes.append(g)
            plain = IN / "c5plain" / "files" / p.name.removesuffix(".gz")
            with gzip.open(p, "rb") as src, open(plain, "wb") as dst:
                while chunk := src.read(1 << 20):
                    dst.write(chunk)
                    cat.write(chunk)
            plain_objs.append({**o, "path": str(plain)})
    gz_bytes = sum(Path(o["path"]).stat().st_size for o in sample)
    fa_bytes = sum(Path(o["path"]).stat().st_size for o in plain_objs)
    for d in ("c5gz", "c5plain", "c5few"):
        json.dump(genomes, open(IN / d / "genomes.json", "w"))
    json.dump(sample, open(IN / "c5gz" / "objects.json", "w"))
    json.dump(plain_objs, open(IN / "c5plain" / "objects.json", "w"))
    json.dump([{**sample[0], "path": str(few)}], open(IN / "c5few" / "objects.json", "w"))
    report["c5"] = {"objects": len(sample), "genomes": len(genomes), "gz_bytes": gz_bytes, "fasta_bytes": fa_bytes,
                    "records": n_rec, "bases": n_bases, "rule": "objects.json[::10]"}
    spx_objs = {o["path"]: o for o in json.load(open(SPX / "objects.json"))}
    paths = sorted(spx_objs, key=lambda p: (Path(p).stat().st_size, p))
    mid = len(paths) // 2
    pick = paths[mid - 2: mid + 2]
    for tag, minlen in (("spx", 0), ("spx5000", 5000)):
        seen = set()
        genomes = []
        bases = kept_bases = 0
        for p in pick:
            uids = []
            for r in fasta.iter_records(Path(p)):
                bases += len(r.seq)
                if len(r.seq) >= max(minlen, 1):
                    uids.append(ids.contig_uid(r.seq))
                    kept_bases += len(r.seq)
            g = _genome(uids, seen)
            if g:
                genomes.append(g)
        json.dump(genomes, open(IN / tag / "genomes.json", "w"))
        json.dump([spx_objs[p] for p in pick], open(IN / tag / "objects.json", "w"))
        report[tag] = {"objects": len(pick), "genomes": len(genomes), "gz_bytes": sum(Path(p).stat().st_size for p in pick),
                       "bases_read": bases, "bases_kept": kept_bases, "files": [Path(p).name for p in pick],
                       "min_contig_len": minlen or None}
    json.dump(report, open(IN / "report.json", "w"), indent=1)
    print(json.dumps(report, indent=1))


def make_scale() -> None:
    """The `scale` arm's sample: every 3rd object of C5's objects.json (index 0, 3, 6, ...), `.fna.gz` as C5 read
    them, one genome per object as above: 3.3 times the base sample, to see whether the cost per byte grows with the
    size of the stage."""
    from biobox import fasta, ids

    (IN / "c5x3").mkdir(parents=True, exist_ok=True)
    sample = json.load(open(C5 / "objects.json"))[::3]
    seen: set = set()
    genomes = []
    bases = 0
    for o in sample:
        uids = []
        for r in fasta.iter_records(Path(o["path"])):
            uids.append(ids.contig_uid(r.seq))
            bases += len(r.seq)
        g = _genome(uids, seen)
        if g:
            genomes.append(g)
    json.dump(genomes, open(IN / "c5x3" / "genomes.json", "w"))
    json.dump(sample, open(IN / "c5x3" / "objects.json", "w"))
    rep = json.load(open(IN / "report.json"))
    rep["c5x3"] = {"objects": len(sample), "genomes": len(genomes), "bases": bases,
                   "gz_bytes": sum(Path(o["path"]).stat().st_size for o in sample), "rule": "objects.json[::3]"}
    json.dump(rep, open(IN / "report.json", "w"), indent=1)
    print(json.dumps(rep["c5x3"]))


def line_counts(path: Path, min_len: int) -> dict:
    """FASTA bytes, records and bases of `path` in all and under `min_len` bp, and the FASTA bytes of the records
    under it (a record's bytes: its header line and its sequence lines, newlines included), by a loop over the
    decompressed lines. A base is a byte of a sequence line other than CR and LF; lines that are blank are no bases,
    as `fasta.iter_records` reads them."""
    import gzip

    out = dict.fromkeys(("fasta_bytes", "records", "bases", "records_below", "bases_below", "bytes_below"), 0)
    cur = None  # [bytes, bases] of the record being read

    def close() -> None:
        if cur is not None:
            out["records"] += 1
            out["bases"] += cur[1]
            if cur[1] < min_len:
                out["records_below"] += 1
                out["bases_below"] += cur[1]
                out["bytes_below"] += cur[0]

    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:
        for raw in fh:
            out["fasta_bytes"] += len(raw)
            if raw.startswith(b">"):
                close()
                cur = [len(raw), 0]
            elif cur is not None:
                cur[0] += len(raw)
                if raw.strip():
                    cur[1] += len(raw.rstrip(b"\r\n"))
        close()
    return out


def b61_pick() -> list[list[str]]:
    """Card B-61's 12 files of SPX1 in 6 groups of 2 (the rule in the module's docstring)."""
    objs = json.load(open(SPX / "objects.json"))
    paths = sorted((o["path"] for o in objs), key=lambda p: (Path(p).stat().st_size, p))
    mid = len(paths) // 2
    b58 = set(paths[mid - 2: mid + 2])  # card B-58's `spx` files
    rest = [p for p in paths if p not in b58]
    pick = [rest[int((i + 0.5) * len(rest) / (2 * N_GROUPS))] for i in range(2 * N_GROUPS)]
    return [pick[2 * k: 2 * k + 2] for k in range(N_GROUPS)]


def make_groups() -> None:
    """Card B-61's inputs: for each group k, `spxg<k>` (every contig kept) and `spxg<k>m` (contigs of MIN_LEN or
    more), one genome per object as `make_inputs` makes them, and the counts of each file in report.json."""
    from biobox import fasta, ids

    spx_objs = {o["path"]: o for o in json.load(open(SPX / "objects.json"))}
    rep = json.load(open(IN / "report.json")) if (IN / "report.json").exists() else {}
    IN.mkdir(parents=True, exist_ok=True)
    paths = sorted(spx_objs, key=lambda p: (Path(p).stat().st_size, p))
    mid = len(paths) // 2
    b58 = [line_counts(Path(p), MIN_LEN) | {"name": Path(p).name, "gz_bytes": Path(p).stat().st_size}
           for p in paths[mid - 2: mid + 2]]  # card B-58's `spx` files, counted the same way (no arm is run on them)
    rep["spx_b58"] = {"files": b58, **{key: sum(c[key] for c in b58) for key in b58[0] if key != "name"}}
    for k, files in enumerate(b61_pick(), start=1):
        per_file = []
        for p in files:
            c = line_counts(Path(p), MIN_LEN)
            c.update(name=Path(p).name, gz_bytes=Path(p).stat().st_size)
            per_file.append(c)
        for tag, minlen in ((f"spxg{k}", 0), (f"spxg{k}m", MIN_LEN)):
            (IN / tag).mkdir(parents=True, exist_ok=True)
            seen: set = set()
            genomes = []
            for p in files:
                uids = [ids.contig_uid(r.seq) for r in fasta.iter_records(Path(p)) if len(r.seq) >= max(minlen, 1)]
                g = _genome(uids, seen)
                if g:
                    genomes.append(g)
            json.dump(genomes, open(IN / tag / "genomes.json", "w"))
            json.dump([spx_objs[p] for p in files], open(IN / tag / "objects.json", "w"))
        tot = {key: sum(c[key] for c in per_file) for key in per_file[0] if key != "name"}
        rep[f"spxg{k}"] = {"files": per_file, **tot, "min_contig_len": MIN_LEN}
        json.dump(rep, open(IN / "report.json", "w"), indent=1)
        print(json.dumps({f"spxg{k}": tot}), flush=True)


def spx1_plan() -> None:
    """Card B-61, read-only on SPX1's staged inputs (inside the image, the B-61 tree on PYTHONPATH): (1) the CPU time
    of `stage.fasta_bytes` against `stage.fasta_counts` on card B-61's 12 files (what counting the records' lengths
    adds to the read `object_sizes` makes anyway); (2) SPX1's ingest estimate and first-pass cut by the code, as
    `build-stage --min-contig-len 5000` makes them (`object_sizes`, `below_bytes`, `gz_bytes`, `ingest_cpu_min`,
    `_pass1_bytes` at a budget of 300 CPU-minutes, `store.pass1_shards`), with each shard's estimate. Nothing is
    written but plan.json under B58_DIR. SPX1 holds no contig of the root (card C-48e: 0 shared), so held is 0."""
    import resource

    from biobox import stage
    from biobox.cli import _source_objects
    from biobox.l0seq import store

    def cpu() -> float:
        r = resource.getrusage(resource.RUSAGE_SELF)
        return r.ru_utime + r.ru_stime

    out: dict = {"timing": []}
    for files in b61_pick():
        for p in files:
            t0 = cpu()
            n = stage.fasta_bytes(Path(p))
            t1 = cpu()
            n2, below = stage.fasta_counts(Path(p), MIN_LEN)
            t2 = cpu()
            out["timing"].append({"name": Path(p).name, "fasta_bytes": n, "below_bytes": below, "same": n == n2,
                                  "cpu_s_fasta_bytes": round(t1 - t0, 3), "cpu_s_fasta_counts": round(t2 - t1, 3)})
            print(json.dumps(out["timing"][-1]), flush=True)
    objs_path = SPX / "objects.json"
    t0 = cpu()
    sizes = stage.object_sizes(objs_path, MIN_LEN)
    t_count = cpu() - t0
    n, gz = sum(f for _, _, f in sizes), stage.gz_bytes(sizes)
    below = stage.below_bytes(objs_path, MIN_LEN)
    on_disk = sum(d for _, d, _ in sizes)
    est = stage.ingest_cpu_min(n, 0, gz_bytes=gz, below_bytes=below)
    prof = type("P", (), {"cpu_min_per_process": 300})()
    per_shard = stage._pass1_bytes(prof, on_disk, n, gz, below)
    objects = _source_objects(str(objs_path))
    order = store.pass1_order(objects)
    below_of = stage._BELOW[(stage._object_key(objs_path)[1], MIN_LEN)]
    by_name = {Path(p).name: (d, f, b) for (p, d, f), b in zip(sizes, below_of)}  # the 304 names are distinct
    shards = []
    for members in store.pass1_shards(objects, per_shard):
        d = sum(by_name[Path(order[i].path).name][0] for i in members)
        f = sum(by_name[Path(order[i].path).name][1] for i in members)
        b = sum(by_name[Path(order[i].path).name][2] for i in members)
        shards.append({"objects": len(members), "bytes_on_disk": d, "fasta_bytes": f, "below_bytes": b,
                       "est_cpu_min_one_process": round(stage.ingest_cpu_min(f, 0, gz_bytes=f, below_bytes=b), 1)})
    out.update(objects=len(sizes), bytes_on_disk=on_disk, fasta_bytes=n, gz_bytes=gz, below_bytes=below,
               cpu_s_object_sizes_with_counts=round(t_count, 1), ingest_cpu_min=round(est, 1),
               ingest_cpu_min_before_b61=round(stage.ingest_cpu_min(n, 0), 1), pass1_bytes=per_shard, shards=shards,
               rates={"stored": stage.INGEST_CPU_MIN, "below": stage.INGEST_BELOW_CPU_MIN,
                      "gz": stage.INGEST_GZ_CPU_MIN, "fixed": stage.INGEST_FIXED_CPU_MIN})
    json.dump(out, open(B / "plan.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "timing"}, indent=1))


G1 = Path("/tmp/if2453/s2/inputs/G1")


def make_g1() -> None:
    """Card B-61, at the coordinator's request: the `g1` arm's sample, every 10th object of stage G1's objects.json
    (index 0, 10, 20, ...; GlobDB MAGs as `.fa.gz`, a third genome type beside C5's complete genomes and SPX1's
    unbinned contigs), one genome per object as `make_inputs`, every contig kept; counted by `line_counts`."""
    from biobox import fasta, ids

    (IN / "g1").mkdir(parents=True, exist_ok=True)
    sample = json.load(open(G1 / "objects.json"))[::10]
    seen: set = set()
    genomes, cs = [], []
    for o in sample:
        g = _genome([ids.contig_uid(r.seq) for r in fasta.iter_records(Path(o["path"]))], seen)
        if g:
            genomes.append(g)
        cs.append(line_counts(Path(o["path"]), MIN_LEN))
    json.dump(genomes, open(IN / "g1" / "genomes.json", "w"))
    json.dump(sample, open(IN / "g1" / "objects.json", "w"))
    rep = json.load(open(IN / "report.json"))
    rep["g1"] = {"objects": len(sample), "genomes": len(genomes), **{k: sum(c[k] for c in cs) for k in cs[0]},
                 "gz_bytes": sum(Path(o["path"]).stat().st_size for o in sample), "rule": "G1 objects.json[::10]"}
    json.dump(rep, open(IN / "report.json", "w"), indent=1)
    print(json.dumps(rep["g1"]), flush=True)


def c5_counts() -> None:
    """Card B-61: card B-58's C5 samples (`objects.json[::10]`, `[::3]`) counted by `line_counts` (their bytes of FASTA
    and bases; card B-58's note gives the bytes of the first only), into report.json as `c5` and `c5x3`."""
    rep = json.load(open(IN / "report.json"))
    objs = json.load(open(C5 / "objects.json"))
    for tag, sample in (("c5", objs[::10]), ("c5x3", objs[::3])):
        cs = [line_counts(Path(o["path"]), MIN_LEN) for o in sample]
        rep[tag] = {"objects": len(sample), **{k: sum(c[k] for c in cs) for k in cs[0]},
                    "gz_bytes": sum(Path(o["path"]).stat().st_size for o in sample)}
        print(json.dumps({tag: rep[tag]}), flush=True)
    json.dump(rep, open(IN / "report.json", "w"), indent=1)


F_FIXED, R_STORED = 0.07, 4.11  # card B-52's constants, fixed (every plain stage's estimate stays as it is)


def arm_table() -> list[dict]:
    """Every arm of cards B-58 and B-61: what it read (n bytes of FASTA, z of them out of .gz, h bases held, d bytes in
    records under 5,000 bp) and its records (cpu_min of the ingest step, every run)."""
    rep = json.load(open(IN / "report.json"))
    b58 = [json.loads(ln) for ln in open(Path("/tmp/if2453/b58/results.jsonl"))]
    b61 = [json.loads(ln) for ln in open(RESULTS)]

    def runs(rows, arm, key="ingest", reps=None):
        return [r[key] for r in rows if r["arm"] == arm and r.get(key) is not None and (reps is None or r["rep"] in reps)]

    c5, c5x3, s58 = rep["c5"], rep["c5x3"], rep["spx_b58"]
    n5 = c5["fasta_bytes"]
    t = [dict(arm="gz", n=n5, z=n5, h=0, d=0, runs=runs(b58, "gz")),
         dict(arm="plain", n=n5, z=0, h=0, d=0, runs=runs(b58, "plain")),
         dict(arm="few", n=n5, z=0, h=0, d=0, runs=runs(b58, "few")),
         dict(arm="view", n=n5, z=n5, h=0, d=0, runs=runs(b58, "view")),
         dict(arm="held A", n=n5, z=n5, h=0, d=0, runs=runs(b58, "held")),
         dict(arm="held B", n=n5, z=n5, h=c5["bases"], d=0, runs=runs(b58, "held", "ingest_b")),
         dict(arm="scale", n=c5x3["fasta_bytes"], z=c5x3["fasta_bytes"], h=0, d=0, runs=runs(b58, "scale")),
         dict(arm="spx", n=s58["fasta_bytes"], z=s58["fasta_bytes"], h=0, d=0, runs=runs(b58, "spx")),
         dict(arm="spx5000", n=s58["fasta_bytes"], z=s58["fasta_bytes"], h=0, d=s58["bytes_below"],
              runs=runs(b58, "spx5000"))]
    for k in range(1, N_GROUPS + 1):
        g = rep[f"spxg{k}"]
        n = g["fasta_bytes"]
        t += [dict(arm=f"spxg{k}", n=n, z=n, h=0, d=0, runs=runs(b61, f"spxg{k}")),
              dict(arm=f"spxg{k} B", n=n, z=n, h=g["bases"], d=0, runs=runs(b61, f"spxg{k}", "ingest_b")),
              dict(arm=f"spxg{k}m", n=n, z=n, h=0, d=g["bytes_below"], runs=runs(b61, f"spxg{k}m"))]
    if "g1" in rep:
        n = rep["g1"]["fasta_bytes"]
        t.append(dict(arm="g1", n=n, z=n, h=0, d=0, runs=runs(b61, "g1")))
    for r in t:
        r["mean"] = sum(r["runs"]) / len(r["runs"]) if r["runs"] else None
    return t


def est(r: dict, g: float, d_rate: float) -> float:
    return (F_FIXED + R_STORED * max(0, r["n"] - r["h"] - r["d"]) / 1e9 + d_rate * r["d"] / 1e9
            + g * r["z"] / 1e9)


def fit_below(rows: list[dict], g: float) -> float:
    """The rate per 10^9 bytes under the least length that minimises the squared relative error of `rows`."""
    num = den = 0.0
    for r in rows:
        base = est(r, g, 0.0)
        x = r["d"] / 1e9
        num += x * (r["mean"] - base) / r["mean"] ** 2
        den += x * x / r["mean"] ** 2
    return num / den


def fit_joint(rows: list[dict]) -> tuple[float, float]:
    """(gzip rate, rate under the least length) per 10^9 bytes minimising the squared relative error of `rows`, both
    free (two normal equations); the stored rate and the start-up cost stay card B-52's."""
    a11 = a12 = a22 = b1 = b2 = 0.0
    for r in rows:
        w = 1 / r["mean"] ** 2
        y = r["mean"] - est(r, 0.0, 0.0)
        z, d = r["z"] / 1e9, r["d"] / 1e9
        a11 += w * z * z
        a12 += w * z * d
        a22 += w * d * d
        b1 += w * z * y
        b2 += w * d * y
    det = a11 * a22 - a12 * a12
    return (b1 * a22 - b2 * a12) / det, (a11 * b2 - a12 * b1) / det


# the arms the joint fit uses: every arm whose records are long contigs, or which drops the short ones, and is not
# start-up bound (spxg1, 3 MB). The arms that store or hold SPX's short contigs (spx, spxg<k>, spxg<k> B) carry a cost
# per record the code cannot know before ingest without counting records; they are judged, not fitted.
JOINT = ("gz", "plain", "few", "view", "held A", "held B", "scale", "spx5000", "spxg2m", "spxg3m", "spxg4m", "spxg5m",
         "spxg6m", "g1")


def fit() -> None:
    """Card B-61's fit, printed: the gzip rate from the C5 arms (paired, unpaired), the rate of bytes under the least
    length from the `m` arms at each gzip rate with its leave-one-out errors, and every arm against the estimate."""
    t = arm_table()
    by = {r["arm"]: r for r in t}
    gz, pl = by["gz"]["runs"], by["plain"]["runs"]
    n5 = by["gz"]["n"] / 1e9
    paired = sum(a - b for a, b in zip(gz[:3], pl)) / 3 / n5
    unpaired = (sum(gz) / len(gz) - sum(pl) / len(pl)) / n5
    print(json.dumps({"gz_paired": paired, "gz_unpaired": unpaired,
                      "gz_per_rep": [(a - b) / n5 for a, b in zip(gz[:3], pl)]}))
    t = [r for r in t if r["mean"] is not None]
    m = [r for r in t if r["arm"] in ("spx5000",) or (r["arm"].startswith("spxg") and r["arm"].endswith("m"))]
    for g in (round(paired, 4), round(unpaired, 4)):
        d_all = fit_below(m, g)
        loo = []
        for r in m:
            rest = [x for x in m if x is not r]
            d_r = fit_below(rest, g)
            loo.append((r["arm"], round(d_r, 4), round((r["mean"] - est(r, g, d_r)) / est(r, g, d_r), 4)))
        each = [(r["arm"], round((r["mean"] - est(r, g, 0.0)) / (r["d"] / 1e9), 4)) for r in m]
        print(json.dumps({"g": g, "below_rate": d_all, "each_arm_alone": each, "leave_one_out": loo}))
        print(json.dumps({"g": g, "ratio_record_over_estimate": [(r["arm"], round(r["mean"] / est(r, g, d_all), 3))
                                                                 for r in t]}))
    j = [r for r in t if r["arm"] in JOINT]
    g, d = fit_joint(j)
    loo = []
    for r in j:
        gr, dr = fit_joint([x for x in j if x is not r])
        loo.append((r["arm"], round(gr, 4), round(dr, 4), round((r["mean"] - est(r, gr, dr)) / est(r, gr, dr), 4)))
    print(json.dumps({"joint": {"g": g, "below_rate": d, "arms": [r["arm"] for r in j], "dof": len(j) - 2,
                                "leave_one_out": loo,
                                "ratio": [(r["arm"], round(r["mean"] / est(r, g, d), 3)) for r in t]}}))
    spread = {r["arm"]: (round(min(r["runs"]), 4), round(max(r["runs"]), 4), len(r["runs"])) for r in t}
    print(json.dumps({"runs_min_max_n": spread}))
    for r in t:
        print(json.dumps({k: r[k] for k in ("arm", "n", "z", "h", "d", "runs", "mean")}))


# --- the arms (on the host) ------------------------------------------------------------------------------------


def now() -> str:
    return time.strftime("%F %T")


def sh(*argv: str, check: bool = True) -> str:
    r = subprocess.run(argv, capture_output=True, text=True)
    if check and r.returncode:
        raise SystemExit(f"{argv}: exit {r.returncode}: {r.stderr[-500:]}")
    return r.stdout


def snapshot() -> dict:
    units = sh("systemctl", "--user", "list-units", "--type=service", "--type=scope", "--state=active", "--no-legend",
               "--plain", check=False)
    names = sorted(ln.split()[0] for ln in units.splitlines() if ln.strip())
    return {"time": now(), "load": open("/proc/loadavg").read().split()[:3], "units": names}


def profile(tmp: Path) -> Path:
    p = B / "profile.yaml"
    p.write_text(PROFILE.format(tmp=tmp, b=B))
    return p


def launch(unit: str, root: Path, stage: str, inp: str, extra: list[str], view: bool) -> tuple[int, dict]:
    tmp = B / "tmp"
    for d in (tmp, B / "durable", B / "bulk", LOGS):
        d.mkdir(parents=True, exist_ok=True)
    prof = profile(tmp)
    log = LOGS / f"{unit}.log"
    log.unlink(missing_ok=True)
    binds = ["--bind", "/tmp/if2453"]
    if view:
        binds += ["--bind", f"{DATA_ROOT}:{DATA_ROOT}:ro"]
    sh("bash", f"{TREE}/scripts/pod_run.sh", "--profile", str(prof), "--image", IMAGE, "--digest", DIGEST, "--mem", "32G",
       "--threads", str(THREADS), *binds, "--env", f"PYTHONPATH={TREE}/src", "--log", str(log), unit,
       "build-stage", "--profile", str(prof), "--root", str(root), "--stage", stage, "--depth", "stored",
       "--genomes", str(IN / inp / "genomes.json"), "--objects", str(IN / inp / "objects.json"),
       "--threads", str(THREADS), *extra)
    time.sleep(5)
    t0 = time.time()
    while subprocess.run(["systemctl", "--user", "is-active", "--quiet", unit]).returncode == 0:
        if time.time() - t0 > 4 * 3600:
            raise SystemExit(f"{unit}: still running after 4 hours")
        time.sleep(5)
    rc = -1
    for ln in log.read_text().splitlines():
        if ln.startswith(f"pod_run: {unit} exited "):
            rc = int(ln.rsplit(" ", 1)[1])
    rec = {}
    for step in ("inputs", "ingest"):
        p = root / "runs" / f"stage={stage}" / f"{step}.json"
        if p.exists():
            d = json.loads(p.read_text())
            rec[step] = d.get("cpu_min")
            rec[f"{step}_wall_s"] = d.get("wall_s", d.get("elapsed_s"))
    for ln in log.read_text().splitlines():
        if " usage " in ln and ln.startswith(f"pod_run: {unit} usage"):
            rec["usage"] = ln.split(" usage ", 1)[1]
    return rc, rec


def run_arm(arm: str, rep: int) -> dict:
    if (DATA_ROOT / "STOP").exists():
        raise SystemExit(f"{DATA_ROOT}/STOP exists: not starting {arm}.{rep}")
    root = ROOTS / f"{arm}-{rep}"
    shutil.rmtree(root, ignore_errors=True)
    out = {"arm": arm, "rep": rep, "start": snapshot()}
    view = arm == "view"
    if view:
        subprocess.run(["python3", "-I", f"{TREE}/scripts/stress/build_stage_pod/c26_view.py", str(DATA_ROOT), str(root),
                        EARLIER], check=True)
    else:
        root.mkdir(parents=True)
    unit = f"{TAG}-{arm}-{rep}"
    spec = {"gz": ("c5gz", []), "plain": ("c5plain", []), "few": ("c5few", []), "view": ("c5gz", []), "scale": ("c5x3", []),
            "held": ("c5gz", []), "spx": ("spx", []), "spx5000": ("spx5000", ["--min-contig-len", "5000"])}
    spec["g1"] = ("g1", [])  # card B-61: every 10th object of G1
    for k in range(1, N_GROUPS + 1):  # card B-61
        spec[f"spxg{k}"] = (f"spxg{k}", [])
        spec[f"spxg{k}m"] = (f"spxg{k}m", ["--min-contig-len", str(MIN_LEN)])
    inp, extra = spec[arm]
    rc, rec = launch(unit, root, "BA", inp, extra, view)
    out.update(exit=rc, **rec)
    if (arm == "held" or (arm.startswith("spxg") and not arm.endswith("m"))) and rc == 0:
        rc2, rec2 = launch(f"{unit}-b", root, "BB", inp, extra, False)
        out.update(exit_b=rc2, ingest_b=rec2.get("ingest"), inputs_b=rec2.get("inputs"))
    out["end"] = snapshot()
    shutil.rmtree(root, ignore_errors=True)
    return out


def run(reps: int, arms: list[str], rep0: int = 0) -> None:
    B.mkdir(parents=True, exist_ok=True)
    for rep in range(rep0 + 1, rep0 + reps + 1):
        for arm in arms:
            res = run_arm(arm, rep)
            with open(RESULTS, "a") as fh:
                fh.write(json.dumps(res) + "\n")
            print(now(), arm, rep, res.get("exit"), res.get("ingest"), res.get("ingest_b"), flush=True)


def report() -> None:
    rep = json.load(open(IN / "report.json"))
    rows = [json.loads(ln) for ln in open(RESULTS)]
    print("arm\trep\texit\tingest_cpu_min\tingest_b_cpu_min\tload_start\tload_end\tother_units_start")
    for r in rows:
        mine = f"{TAG}-{r['arm']}-{r['rep']}"
        other = [u for u in r["start"]["units"] if mine not in u and f"{TAG}-" not in u and "session" not in u]
        print(f"{r['arm']}\t{r['rep']}\t{r.get('exit')}\t{r.get('ingest')}\t{r.get('ingest_b')}\t"
              f"{r['start']['load'][0]}\t{r['end']['load'][0]}\t{' '.join(other)}")
    print(json.dumps(rep))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("inputs", "scale-inputs", "groups-inputs", "g1-inputs", "spx1-plan", "c5-counts",
                                     "fit", "run", "report"))
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--rep0", type=int, default=0, help="rep numbers start after this (a second batch of arms)")
    ap.add_argument("--arms", default=",".join(ARMS))
    a = ap.parse_args()
    if a.mode == "inputs":
        make_inputs()
    elif a.mode == "scale-inputs":
        make_scale()
    elif a.mode == "groups-inputs":
        make_groups()
    elif a.mode == "spx1-plan":
        spx1_plan()
    elif a.mode == "c5-counts":
        c5_counts()
    elif a.mode == "g1-inputs":
        make_g1()
    elif a.mode == "fit":
        fit()
    elif a.mode == "run":
        run(a.reps, a.arms.split(","), a.rep0)
    else:
        report()
