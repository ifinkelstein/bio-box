#!/usr/bin/env python3
"""Ladder stage driver (host python3, stdlib only): one chained, resumable run of a stage.

    python3 driver.py KEY        KEY in A, B, Bp, C

Each build step runs as its own process inside the image (sh bbn.sh step.py KEY STEP); the driver
times it and takes its CPU time and peak RSS (whole process tree) from os.wait4. result.json
(/tmp/if2453/ladder/results/<KEY>/result.json) is rewritten after every step; a rerun skips the
steps recorded there without error. Hard steps stop the stage on failure; soft steps (checks,
queries, the durable copy) are recorded and the stage goes on.

Threads: at most 48 in total at night (20:00-08:00 local) and 32 by day, 16 if the machine's
1-minute load excluding our registered threads is above 10. Units register their threads in
run/<unit>.threads ("n pid"); a multi-threaded step waits until its share is free.
"""

import json
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

W = Path("/tmp/if2453/ladder")
S = W / "scripts"
LOGS, RUN = W / "logs", W / "run"
ROOT = W / "root"
STOP = Path("/stor/scratch/Wilke/if2453/all-the-data/STOP")
SCR = Path("/stor/scratch/Wilke/if2453/all-the-data/provisional/rL")
PFAM_HMM = Path("/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2/Pfam-A.hmm.gz")
KEY = sys.argv[1]
RELEASE = {"A": "rL.1", "B": "rL.2", "Bp": "rL.2p", "C": "rL.3"}[KEY]
STAGE = {"A": "S1", "B": "P1", "Bp": "P1", "C": "P2"}[KEY]
STEPS = {
    "A": ["ingest", "genome_input", "genecall", "igr", "tapes", "cluster", "domains", "tokens", "postings", "dim",
          "catalog", "equiv", "queries", "sync"],
    "B": ["wait_inputs", "ingest", "genome_input", "genecall", "igr", "tapes", "cluster", "tokens", "postings", "dim",
          "catalog", "appendonly", "queries", "sync"],
    "Bp": ["wait_inputs", "pfam_plan", "pfam_scan", "pfam_merge", "cluster_copy", "tokens", "postings", "dim",
           "catalog", "appendonly", "queries", "sync"],
    "C": ["wait_inputs", "ingest", "genome_input", "genecall", "igr", "tapes", "cluster", "wait_pfam", "tokens",
          "postings", "dim", "catalog", "appendonly", "queries", "sync"],
}[KEY]
SOFT = {"equiv", "appendonly", "queries", "sync"}
MULTI = {"genecall": 16, "cluster": 16}
OUT = W / f"results/{KEY}/result.json"
OUT.parent.mkdir(parents=True, exist_ok=True)
RUN.mkdir(exist_ok=True)

R = json.loads(OUT.read_text()) if OUT.exists() else {"key": KEY, "stage": STAGE, "release": RELEASE, "steps": {}}
R["steps"] = {k: v for k, v in R["steps"].items() if not v.get("error")}
R["runs"] = R.get("runs", []) + [{"started": datetime.now().isoformat(timespec="seconds"), "pid": os.getpid()}]


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def log(msg: str) -> None:
    print(f"[{now()}] {KEY}: {msg}", flush=True)


def save() -> None:
    t = OUT.with_suffix(".tmp")
    t.write_text(json.dumps(R, indent=1, default=str))
    t.replace(OUT)


def register(n: int, unit: str = KEY) -> None:
    (RUN / f"{unit}.threads").write_text(f"{n} {os.getpid()}\n")


def ours(exclude: str | None = None) -> int:
    total = 0
    for f in RUN.glob("*.threads"):
        if f.stem == exclude:
            continue
        try:
            n, pid = (int(x) for x in f.read_text().split()[:2])
            os.kill(pid, 0)
        except (ValueError, OSError):
            continue
        total += n
    return total


def budget(me: str = KEY) -> tuple[int, dict]:
    """Threads this unit may use now, and the numbers behind it."""
    h = datetime.now().hour
    cap = 48 if (h >= 20 or h < 8) else 32
    load1 = float(Path("/proc/loadavg").read_text().split()[0])
    mine = ours()
    foreign = max(0.0, load1 - mine)
    if foreign > 10:
        cap = 16
    return cap - ours(exclude=me), {"hour": h, "cap": cap, "load1": load1, "ours": mine, "foreign": round(foreign, 1)}


def stopped() -> bool:
    return STOP.exists()


def marker(name: str) -> bool:
    return (LOGS / name).exists()


def wait_for(cond, what: str, timeout_h: float = 24) -> bool:
    t0 = time.time()
    while not cond():
        if stopped():
            raise RuntimeError("STOP present while waiting for " + what)
        if time.time() - t0 > timeout_h * 3600:
            return False
        time.sleep(60)
    return True


def stage_ok(key: str, need: list[str]) -> bool:
    p = W / f"results/{key}/result.json"
    if not p.exists():
        return False
    st = json.loads(p.read_text())["steps"]
    return all(s in st and not st[s].get("error") for s in need)


BUILD_B = ["ingest", "genome_input", "genecall", "igr", "tapes", "cluster", "tokens", "postings", "dim", "catalog"]


def wait_inputs() -> dict:
    """Prerequisites of the stage: earlier stage built (its .done, or .failed after its build steps), the download."""
    t0 = time.time()
    info = {}
    if KEY == "B":
        wait_for(lambda: marker("A.done") or marker("A.failed"), "stage A")
        if not stage_ok("A", ["ingest", "genome_input", "genecall", "igr", "tapes", "cluster", "domains", "tokens",
                              "postings", "dim", "catalog"]):
            raise RuntimeError("stage A did not build; B needs its root")
        wait_for(lambda: marker("fetch-P1.done") or marker("fetch-P1.failed"), "fetch P1")
        if marker("fetch-P1.failed"):
            raise RuntimeError("fetch P1 failed")
    elif KEY == "Bp":
        wait_for(lambda: marker("B.done") or marker("B.failed"), "stage B")
        if not stage_ok("B", BUILD_B):
            raise RuntimeError("stage B did not build")
    elif KEY == "C":
        wait_for(lambda: marker("B.done") or marker("B.failed"), "stage B")
        if not stage_ok("B", BUILD_B):
            raise RuntimeError("stage B did not build; C adds to it")
        wait_for(lambda: marker("fetch-P2.done") or marker("fetch-P2.failed"), "fetch P2")
        if marker("fetch-P2.failed"):
            raise RuntimeError("fetch P2 failed")
    if KEY in ("B", "C"):
        fr = json.loads((W / f"genomes/{STAGE}/fetch_result.json").read_text())
        info = {k: fr.get(k) for k in ("genomes_ok", "genomes_failed", "failed", "bp", "records", "bytes", "wall_s",
                                       "requests", "retries", "fallback_single", "fallback_ftp")}
        info["bp_vs_gtdb_mismatch_n"] = len(fr.get("bp_vs_gtdb_mismatch", {}))
    info["waited_s"] = round(time.time() - t0, 1)
    return info


def wait_pfam() -> dict:
    """C's tokens include stage B's Pfam rows: wait for Bp's merge so rL.3 does not depend on timing."""
    t0 = time.time()
    ok = wait_for(lambda: stage_ok("Bp", ["pfam_merge"]) or marker("Bp.failed") or marker("Bp.done"), "Bp pfam_merge",
                  timeout_h=12)
    return {"waited_s": round(time.time() - t0, 1), "bp_pfam_merged": stage_ok("Bp", ["pfam_merge"]), "timed_out": not ok}


def run_proc(cmd: list[str], logf=None) -> tuple[int, float, dict]:
    t = time.perf_counter()
    p = subprocess.Popen(cmd, stdout=logf, stderr=subprocess.STDOUT if logf else None)
    _, status, ru = os.wait4(p.pid, 0)
    p.returncode = os.waitstatus_to_exitcode(status)
    wall = time.perf_counter() - t
    return p.returncode, wall, {"cpu_s": round(ru.ru_utime + ru.ru_stime, 1), "user_s": round(ru.ru_utime, 1),
                                "sys_s": round(ru.ru_stime, 1), "peak_rss_mb": round(ru.ru_maxrss / 1024, 1)}


def image_step(step: str, threads: int) -> dict:
    out = W / f"results/{KEY}/{step}.json"
    if out.exists():
        out.unlink()  # this driver's own file from an earlier failed attempt
    rc, wall, ru = run_proc(["sh", str(S / "bbn.sh"), str(S / "step.py"), KEY, step, "--threads", str(threads)])
    rec = json.loads(out.read_text()) if out.exists() else {"error": f"no step record (rc={rc})"}
    e = {"wall_s": round(wall, 1), "threads": threads, **ru, "cpu_min": round(ru["cpu_s"] / 60, 1),
         "result": rec.get("result"), "error": rec.get("error") or (None if rc == 0 else f"rc={rc}")}
    if rec.get("traceback"):
        e["traceback"] = rec["traceback"]
    return e


def pfam_scan() -> dict:
    """Shards of 4 threads each, at most 8 at once, inside the thread budget; resumable by shard."""
    pf = W / "pfam_P1"
    plan = json.loads((pf / "input/plan.json").read_text())
    (pf / "scan").mkdir(exist_ok=True)
    (pf / "logs").mkdir(exist_ok=True)

    def ok(k):
        j = pf / f"scan/shard-{k}.json"
        return j.exists() and json.loads(j.read_text()).get("ok")

    todo = [s["shard"] for s in plan["shards"] if not ok(s["shard"])]
    running: dict[int, tuple[int, float]] = {}
    per, t0, peak_conc, budgets = {}, time.time(), 0, []
    while todo or running:
        if not stopped():
            avail, info = budget("Bp")
            slots = min(8, max(0, avail) // 4)
            while todo and len(running) < slots:
                k = todo.pop(0)
                lf = open(pf / f"logs/shard-{k}.log", "a")
                cmd = ["sh", str(S / "bbn.sh"), str(S / "scan_shard.py"), str(pf / f"input/shard-{k}.parquet"),
                       str(PFAM_HMM), "gathering", str(pf / f"scan/shard-{k}.parquet"), "--cpus", "4"]
                p = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT)
                lf.close()
                running[p.pid] = (k, time.perf_counter())
                budgets.append({"t": now(), "started": k, **info, "avail": avail})
            peak_conc = max(peak_conc, len(running))
            register(4 * len(running), "Bp")
        elif not running:
            raise RuntimeError("STOP present; shards not started: " + ",".join(map(str, todo)))
        if not running:
            time.sleep(60)
            continue
        pid, status, ru = os.wait4(-1, 0)
        if pid not in running:
            continue
        k, ts = running.pop(pid)
        per[k] = {"wall_s": round(time.perf_counter() - ts, 1), "cpu_min": round((ru.ru_utime + ru.ru_stime) / 60, 1),
                  "peak_rss_mb": round(ru.ru_maxrss / 1024, 1), "rc": os.waitstatus_to_exitcode(status)}
        log(f"pfam shard {k} done rc={per[k]['rc']} {per[k]['wall_s']} s")
        register(4 * len(running), "Bp")
        R["steps"]["pfam_scan_progress"] = {"done": sorted(per), "todo": todo, "running": [v[0] for v in running.values()]}
        save()
    register(0, "Bp")
    R["steps"].pop("pfam_scan_progress", None)
    bad = [s["shard"] for s in plan["shards"] if not ok(s["shard"])]
    if bad:
        raise RuntimeError(f"shards not ok: {bad}")
    return {"shards": len(plan["shards"]), "wall_s": round(time.time() - t0, 1), "peak_concurrent_shards": peak_conc,
            "cpu_min_total": round(sum(v["cpu_min"] for v in per.values()), 1), "per_shard": per,
            "budget_log": budgets}


README = """PROVISIONAL bio-box store "rL" (ladder build). Will be rebuilt; do not cite or depend on paths.

Built on wilkcomp02 with branch plan/task-cards at f9af5b6 (/stor/home/if2453/bio-box-next) inside
bio-box-test_74b3e05.sif. Mixed key format until phase B finishes: tapes/IGRs keyed by tape version,
legacy slot columns still present (slot_protein, left_slot/right_slot, protein_occurrence.slot_ordinal).

One root that accretes; stages append, releases sit beside each other:
  S1  INPHARED 14Apr2025 phages, every record its own genome (kind virus, code 11, topology unknown -> linear)
  P1  250 GTDB R226 species, one complete genome each, from NCBI (kind prokaryote, all contigs linear)
  P2  2,250 further GTDB R226 species, disjoint from P1
  rL.1   S1                 clusters/tokens/postings/dim; Pfam 38.2 on S1 proteins (reused from the pilot)
  rL.2   S1+P1              Pfam only on S1 proteins
  rL.2p  S1+P1              rL.2's clusters (byte copy) + Pfam on P1's new proteins
  rL.3   S1+P1+P2           Pfam on S1 and P1 proteins only (P2 not scanned)

Open: biobox.catalog.open_catalog("<this dir>/catalog/release=<R>.duckdb", local_dir="/tmp/<you>")
(the catalog's views point at this directory; copy the .duckdb to local disk, never open it over NFS).
Manifests: manifests/release=<R>.json (sha256 of every file). Releases synced so far: {releases}
Last sync: {when}
"""


def sync() -> dict:
    if not SCR.parent.exists():
        raise FileNotFoundError(SCR.parent)
    SCR.mkdir(exist_ok=True)
    res = {}
    lf = LOGS / f"sync-{KEY}.rsync.log"
    with open(lf, "w") as fh:
        rc, wall, ru = run_proc(["nice", "-n", "19", "rsync", "-a", "--itemize-changes", "--stats",
                                 "--exclude=/catalog/", "--exclude=/manifests/", f"{ROOT}/", f"{SCR}/"], fh)
    txt = lf.read_text()
    items = [ln for ln in txt.splitlines() if re.match(r"^[<>ch.][fdL]", ln)]
    res["rsync"] = {"rc": rc, "wall_s": round(wall, 1), **ru,
                    "new_files": sum(1 for ln in items if ln.startswith(">f+++++++")),
                    "updated_files": [ln.split(" ", 1)[1] for ln in items if ln.startswith(">f") and not ln.startswith(">f+++++++")],
                    "transferred": re.findall(r"Total transferred file size: ([\d,]+)", txt)}
    if rc != 0:
        raise RuntimeError(f"rsync rc={rc}")
    synced = sorted({*(json.loads((SCR / ".ladder_releases.json").read_text()) if (SCR / ".ladder_releases.json").exists() else []),
                     RELEASE})
    (SCR / ".ladder_releases.json").write_text(json.dumps(synced))
    (SCR / "README").write_text(README.format(releases=", ".join(synced), when=now()))
    out = W / f"results/{KEY}/scratch_catalog.json"
    rc, wall, ru = run_proc(["sh", str(S / "bbn.sh"), str(S / "catalog_build.py"), str(SCR), RELEASE, str(out),
                             "--build-dir", str(W / "catalog_build_scratch"), "--cache-dir", str(W / "catalog_cache/scratch")])
    cb = json.loads(out.read_text()) if out.exists() else {}
    res["scratch_catalog"] = {"rc": rc, "wall_s": round(wall, 1), **ru, "times_s": cb.get("times_s"),
                              "verify_problems": cb.get("verify_problems"), "manifest_hash": cb.get("manifest_hash"),
                              "mismatch": cb.get("view_vs_manifest_row_mismatch")}
    if rc != 0:
        raise RuntimeError(f"scratch catalog_build rc={rc}")
    return res


def host_step(step: str) -> dict:
    t = time.perf_counter()
    c0 = os.times()
    try:
        out = {"wait_inputs": wait_inputs, "wait_pfam": wait_pfam, "pfam_scan": pfam_scan, "sync": sync}[step]()
        err, tb = None, None
    except Exception as e:  # noqa: BLE001
        out, err, tb = None, f"{type(e).__name__}: {e}", traceback.format_exc()
    c1 = os.times()
    e = {"wall_s": round(time.perf_counter() - t, 1), "threads": 1,
         "driver_cpu_s": round((c1.user + c1.system) - (c0.user + c0.system), 1), "result": out, "error": err}
    if tb:
        e["traceback"] = tb
    return e


def main() -> int:
    log(f"driver start pid {os.getpid()}, steps {STEPS}")
    register(0)
    failed_hard = None
    for step in STEPS:
        if step in R["steps"]:
            log(f"{step}: done earlier, skipped")
            continue
        if stopped():
            R["stopped"] = f"STOP present before {step} at {now()}"
            save()
            log(R["stopped"])
            return 75
        threads = 1
        if step in MULTI:
            want = MULTI[step]
            t0 = time.time()
            while True:
                avail, info = budget()
                if avail >= want or (avail >= 4 and time.time() - t0 > 1800):
                    break
                if stopped():
                    R["stopped"] = f"STOP present while waiting for threads for {step}"
                    save()
                    return 75
                time.sleep(60)
            threads = min(want, avail)
            R.setdefault("thread_decisions", {})[step] = {**info, "avail": avail, "threads": threads,
                                                          "waited_s": round(time.time() - t0)}
        register(threads)
        log(f"{step}: start ({threads} threads)")
        started = now()
        e = host_step(step) if step in ("wait_inputs", "wait_pfam", "pfam_scan", "sync") else image_step(step, threads)
        e["started"], e["finished"] = started, now()
        register(0)
        R["steps"][step] = e
        save()
        log(f"{step}: {e['wall_s']} s err={e['error']}")
        if e["error"] and step not in SOFT:
            failed_hard = step
            break
    R["finished"] = now()
    R["failed_hard_step"] = failed_hard
    R["soft_failures"] = [s for s in STEPS if s in SOFT and s in R["steps"] and R["steps"][s].get("error")]
    save()
    try:
        (RUN / f"{KEY}.threads").unlink()
    except OSError:
        pass
    return 1 if failed_hard else 0


sys.exit(main())
