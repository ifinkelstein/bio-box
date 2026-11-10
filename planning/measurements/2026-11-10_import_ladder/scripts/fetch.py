#!/usr/bin/env python3
"""Download genome FASTAs from NCBI by assembly accession (host python3, stdlib only).

    python3 fetch.py STAGE_TSV OUTDIR [--batch 20]

NCBI Datasets v2 REST API, several accessions per request:
    https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/<A1>,<A2>,.../download?include_annotation_type=GENOME_FASTA
An accession missing from a batch zip is retried alone, then through the genomes FTP over HTTPS
(https://ftp.ncbi.nlm.nih.gov/genomes/all/GCx/nnn/nnn/nnn/<acc>_<asm>/<acc>_<asm>_genomic.fna.gz).
Generic user agent; at most 3 requests per second; retries with exponential backoff; resumable
(OUTDIR/lock.json lists every genome done with its sha256, bytes, records, bp, url, time).
Writes OUTDIR/<acc>.fna (plain FASTA, as NCBI ships it), lock.json and fetch_result.json.
"""

import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
import sys
import time
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

UA = "ladder-fetch/1.0 (python-urllib)"
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}/download?include_annotation_type=GENOME_FASTA"
FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}/{}/{}/"
MIN_INTERVAL = 0.34  # <= 3 requests per second
TRIES = 6

ap = argparse.ArgumentParser()
ap.add_argument("tsv", type=Path)
ap.add_argument("out", type=Path)
ap.add_argument("--batch", type=int, default=20)
a = ap.parse_args()
a.out.mkdir(parents=True, exist_ok=True)
tmpd = a.out / ".tmp"
tmpd.mkdir(exist_ok=True)
lock_path = a.out / "lock.json"
res_path = a.out / "fetch_result.json"

rows = []
with open(a.tsv) as fh:
    head = fh.readline().rstrip("\n").split("\t")
    for ln in fh:
        rows.append(dict(zip(head, ln.rstrip("\n").split("\t"))))
want = {r["accession"]: r for r in rows}
lock = json.loads(lock_path.read_text()) if lock_path.exists() else {"genomes": {}}
# resume: keep entries whose file is still there with the recorded size
lock["genomes"] = {k: v for k, v in lock["genomes"].items()
                   if (a.out / v["file_name"]).exists() and (a.out / v["file_name"]).stat().st_size == v["bytes"]}
R = {"tsv": str(a.tsv), "started_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
     "requests": 0, "retries": 0, "bytes_downloaded": 0, "failed": {}, "fallback_single": [], "fallback_ftp": []}
_last = [0.0]


def save():
    for p, obj in ((lock_path, lock), (res_path, R)):
        t = p.with_suffix(".json.tmp")
        t.write_text(json.dumps(obj, indent=1, sort_keys=True))
        t.replace(p)


def get(url: str, dest: Path | None = None) -> bytes | None:
    """GET with rate limit and backoff; to `dest` when given (returns None), else the body."""
    err = None
    for k in range(TRIES):
        wait = MIN_INTERVAL - (time.time() - _last[0])
        if wait > 0:
            time.sleep(wait)
        _last[0] = time.time()
        R["requests"] += 1
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=600) as resp:
                if dest is None:
                    body = resp.read()
                    R["bytes_downloaded"] += len(body)
                    return body
                part = dest.with_name(dest.name + ".part")
                with open(part, "wb") as fo:
                    shutil.copyfileobj(resp, fo, 1 << 20)
                R["bytes_downloaded"] += part.stat().st_size
                part.replace(dest)
                return None
        except urllib.error.HTTPError as e:
            err = f"HTTP {e.code}"
            if e.code in (400, 404):
                break
        except Exception as e:  # noqa: BLE001 - network errors are retried
            err = f"{type(e).__name__}: {e}"
        R["retries"] += 1
        time.sleep(min(120, 2 ** (k + 1)))
    raise RuntimeError(f"{url}: {err}")


def record(acc: str, path: Path, url: str) -> None:
    h = hashlib.sha256()
    n_rec = bp = 0
    with open(path, "rb") as fh:
        for ln in fh:
            h.update(ln)
            if ln.startswith(b">"):
                n_rec += 1
            else:
                bp += len(ln.strip())
    lock["genomes"][acc] = {"file_name": path.name, "file": str(path), "sha256": h.hexdigest(),
                            "bytes": path.stat().st_size, "records": n_rec, "bp": bp, "url": url,
                            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                            "gtdb_genome_size": int(want[acc]["genome_size"])}


def from_zip(zpath: Path, accs: list[str], url: str) -> set[str]:
    got = set()
    with zipfile.ZipFile(zpath) as z:
        for name in z.namelist():
            m = re.fullmatch(r"ncbi_dataset/data/([^/]+)/\1_.*_genomic\.fna", name)
            if not m or m.group(1) not in accs:
                continue
            acc = m.group(1)
            dest = a.out / f"{acc}.fna"
            part = dest.with_name(dest.name + ".part")
            with z.open(name) as fi, open(part, "wb") as fo:
                shutil.copyfileobj(fi, fo, 1 << 20)
            part.replace(dest)
            record(acc, dest, url)
            got.add(acc)
    return got


def via_ftp(acc: str) -> None:
    pre, num = acc.split("_")[0], acc.split("_")[1].split(".")[0]
    d = FTP.format(pre, num[0:3], num[3:6], num[6:9])
    listing = get(d).decode("utf-8", "replace")
    names = sorted(set(re.findall(r'href="(' + re.escape(acc) + r'_[^"/]+)/?"', listing)))
    if not names:
        raise RuntimeError(f"{d}: no directory for {acc}")
    url = f"{d}{names[0]}/{names[0]}_genomic.fna.gz"
    gz = tmpd / f"{acc}.fna.gz"
    get(url, gz)
    dest = a.out / f"{acc}.fna"
    part = dest.with_name(dest.name + ".part")
    with gzip.open(gz, "rb") as fi, open(part, "wb") as fo:
        shutil.copyfileobj(fi, fo, 1 << 20)
    part.replace(dest)
    gz.unlink()
    record(acc, dest, url)


todo = [acc for acc in sorted(want) if acc not in lock["genomes"]]
print(f"{len(want)} wanted, {len(want) - len(todo)} already done, {len(todo)} to fetch", flush=True)
t0 = time.time()
for i in range(0, len(todo), a.batch):
    accs = todo[i : i + a.batch]
    url = API.format(",".join(accs))
    zpath = tmpd / "batch.zip"
    got = set()
    try:
        get(url, zpath)
        got = from_zip(zpath, accs, url)
    except Exception as e:  # noqa: BLE001
        print(f"batch {i}: {e}", flush=True)
    finally:
        zpath.unlink(missing_ok=True)
    for acc in accs:
        if acc in got:
            continue
        try:
            u1 = API.format(acc)
            get(u1, zpath)
            if from_zip(zpath, [acc], u1):
                R["fallback_single"].append(acc)
                continue
        except Exception as e:  # noqa: BLE001
            print(f"{acc}: single: {e}", flush=True)
        finally:
            zpath.unlink(missing_ok=True)
        try:
            via_ftp(acc)
            R["fallback_ftp"].append(acc)
        except Exception as e:  # noqa: BLE001
            R["failed"][acc] = str(e)
            print(f"{acc}: FAILED: {e}", flush=True)
    save()
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {min(i + a.batch, len(todo))}/{len(todo)} "
          f"{R['bytes_downloaded'] / 1e9:.2f} GB", flush=True)

g = lock["genomes"]
R.update({"finished_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), "wall_s": round(time.time() - t0, 1),
          "genomes_ok": len(g), "genomes_failed": len(R["failed"]), "bp": sum(v["bp"] for v in g.values()),
          "records": sum(v["records"] for v in g.values()), "bytes": sum(v["bytes"] for v in g.values()),
          "bp_vs_gtdb_mismatch": {k: [v["bp"], v["gtdb_genome_size"]] for k, v in g.items()
                                  if v["bp"] != v["gtdb_genome_size"]}})
save()
print(json.dumps({k: v for k, v in R.items() if k != "bp_vs_gtdb_mismatch"}, indent=1), flush=True)
sys.exit(0)
