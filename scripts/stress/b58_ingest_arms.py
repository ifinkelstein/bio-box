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
    unit = f"b58-{arm}-{rep}"
    spec = {"gz": ("c5gz", []), "plain": ("c5plain", []), "few": ("c5few", []), "view": ("c5gz", []), "scale": ("c5x3", []),
            "held": ("c5gz", []), "spx": ("spx", []), "spx5000": ("spx5000", ["--min-contig-len", "5000"])}
    inp, extra = spec[arm]
    rc, rec = launch(unit, root, "BA", inp, extra, view)
    out.update(exit=rc, **rec)
    if arm == "held" and rc == 0:
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
        mine = f"b58-{r['arm']}-{r['rep']}"
        other = [u for u in r["start"]["units"] if mine not in u and "b58-" not in u and "session" not in u]
        print(f"{r['arm']}\t{r['rep']}\t{r.get('exit')}\t{r.get('ingest')}\t{r.get('ingest_b')}\t"
              f"{r['start']['load'][0]}\t{r['end']['load'][0]}\t{' '.join(other)}")
    print(json.dumps(rep))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("inputs", "scale-inputs", "run", "report"))
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--rep0", type=int, default=0, help="rep numbers start after this (a second batch of arms)")
    ap.add_argument("--arms", default=",".join(ARMS))
    a = ap.parse_args()
    if a.mode == "inputs":
        make_inputs()
    elif a.mode == "scale-inputs":
        make_scale()
    elif a.mode == "run":
        run(a.reps, a.arms.split(","), a.rep0)
    else:
        report()
