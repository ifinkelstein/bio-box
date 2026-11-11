#!/usr/bin/env python3
"""Download the 50 Gbp build's new genomes, as the import ladder did (host python3, stdlib only).

    python3 r50_fetch.py SELECTION_DIR GENOMES_DIR LOG [--stages P3,P4,P5,P6] [--batch 20]
    python3 r50_fetch.py selftest ZIP OUT_DIR

The ladder's fetch.py (planning/measurements/2026-11-10_import_ladder/scripts/fetch.py) with the same
source, tool and pace: NCBI Datasets v2 REST API, 20 accessions per request, one connection, at most 3
requests per second, generic user agent, no key, no email; an accession missing from a batch is retried
alone, then (only with --ftp-fallback) through the genomes FTP over HTTPS. The fallback finds the
assembly directory by listing its parent, a URL form sources.yaml does not list, so it is off by
default and such an accession is listed as failed. Added here: every FASTA is checked against the md5 the
source ships beside it (md5sum.txt in the Datasets zip; md5checksums.txt in the FTP directory, on the
.fna.gz); a file that fails is fetched once more alone, then listed as failed. One log line per genome.

Layout (the ladder's): GENOMES_DIR/<stage>/<acc>.fna (plain FASTA as NCBI ships it), lock.json (the
ladder's fields plus md5 and md5_source), fetch_result.json. Resumable: a genome in lock.json whose file
is there with the recorded size is not fetched again. Before each request batch the data root's STOP
file is checked. Five refusals in a row (HTTP 403, 429, 503) stop the run. At the end DONE (every genome
verified or listed as failed) or FAILED (stopped early; the reason inside) in LOG's directory.
"""

import argparse
import gzip
import hashlib
import json
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
REFUSAL = (403, 429, 503)
MAX_REFUSALS = 5
STOP = Path("/stor/scratch/Wilke/if2453/all-the-data/STOP")
FNA = re.compile(r"ncbi_dataset/data/([^/]+)/\1_.*_genomic\.fna")


class Refused(RuntimeError):
    pass


class Stopped(RuntimeError):
    pass


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def write_json(p: Path, obj) -> None:
    t = p.with_name(f".{p.name}.tmp")
    t.write_text(json.dumps(obj, indent=1, sort_keys=True) + "\n")
    t.replace(p)


def scan_fna(path: Path) -> dict:
    """sha256, md5, size, record count and bp of a FASTA file."""
    h, m = hashlib.sha256(), hashlib.md5()
    n_rec = bp = 0
    with open(path, "rb") as fh:
        for ln in fh:
            h.update(ln)
            m.update(ln)
            if ln.startswith(b">"):
                n_rec += 1
            else:
                bp += len(ln.strip())
    return {"sha256": h.hexdigest(), "md5": m.hexdigest(), "bytes": path.stat().st_size, "records": n_rec, "bp": bp}


def md5_table(text: str) -> dict[str, str]:
    out = {}
    for ln in text.splitlines():
        parts = ln.split()
        if len(parts) == 2:
            out[parts[1].lstrip("*").removeprefix("./")] = parts[0].lower()
    return out


def extract_zip(zpath: Path, accs: set[str], dest_dir: Path) -> dict[str, dict]:
    """Extract each wanted accession's FASTA from a Datasets zip to dest_dir/<acc>.fna.

    Returns {acc: {"file": path, "md5_expected": ..., "scan": scan_fna(...)}}; a file whose md5 is not
    in md5sum.txt or differs is extracted, scanned and returned with ok False (the caller decides)."""
    got = {}
    with zipfile.ZipFile(zpath) as z:
        names = z.namelist()
        md5s = md5_table(z.read("md5sum.txt").decode()) if "md5sum.txt" in names else {}
        for name in names:
            mt = FNA.fullmatch(name)
            if not mt or mt.group(1) not in accs:
                continue
            acc = mt.group(1)
            dest = dest_dir / f"{acc}.fna"
            part = dest.with_name(dest.name + ".part")
            with z.open(name) as fi, open(part, "wb") as fo:
                shutil.copyfileobj(fi, fo, 1 << 20)
            sc = scan_fna(part)
            exp = md5s.get(name)
            ok = exp is not None and exp == sc["md5"]
            if ok:
                part.replace(dest)
            else:
                part.unlink()
            got[acc] = {"file": dest, "md5_expected": exp, "scan": sc, "ok": ok, "zip_member": name}
    return got


class Fetcher:
    def __init__(self, log: Path):
        self.log = open(log, "a", buffering=1)
        self._last = 0.0
        self.refusals = 0
        self.R = {"requests": 0, "retries": 0, "bytes_downloaded": 0, "refusal_codes": {}}

    def line(self, *f) -> None:
        self.log.write(" ".join([datetime.now().strftime("%Y-%m-%dT%H:%M:%S")] + [str(x) for x in f]) + "\n")

    def get(self, url: str, dest: Path | None = None) -> bytes | None:
        """GET with rate limit and backoff; to `dest` when given (returns None), else the body."""
        err = None
        for k in range(TRIES):
            wait = MIN_INTERVAL - (time.time() - self._last)
            if wait > 0:
                time.sleep(wait)
            self._last = time.time()
            self.R["requests"] += 1
            try:
                req = urllib.request.Request(url, headers={"User-Agent": UA})
                with urllib.request.urlopen(req, timeout=600) as resp:
                    if dest is None:
                        body = resp.read()
                        self.R["bytes_downloaded"] += len(body)
                        self.refusals = 0
                        return body
                    part = dest.with_name(dest.name + ".part")
                    with open(part, "wb") as fo:
                        shutil.copyfileobj(resp, fo, 1 << 20)
                    self.R["bytes_downloaded"] += part.stat().st_size
                    part.replace(dest)
                    self.refusals = 0
                    return None
            except urllib.error.HTTPError as e:
                err = f"HTTP {e.code}"
                if e.code in REFUSAL:
                    self.refusals += 1
                    self.R["refusal_codes"][str(e.code)] = self.R["refusal_codes"].get(str(e.code), 0) + 1
                    self.line("REFUSED", e.code, url[:160], f"consecutive={self.refusals}")
                    if self.refusals >= MAX_REFUSALS:
                        raise Refused(f"{self.refusals} refusals in a row, last HTTP {e.code} for {url[:160]}")
                if e.code in (400, 404):
                    break
            except Exception as e:  # noqa: BLE001 - network errors are retried
                err = f"{type(e).__name__}: {e}"
            self.R["retries"] += 1
            time.sleep(min(120, 2 ** (k + 1)))
        raise RuntimeError(f"{url[:160]}: {err}")


def run_stage(F: Fetcher, stage: str, tsv: Path, out: Path, batch: int, ftp_fallback: bool) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    tmpd = out / ".tmp"
    tmpd.mkdir(exist_ok=True)
    lock_path, res_path = out / "lock.json", out / "fetch_result.json"
    with open(tsv) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        want = {r["accession"]: r for r in (dict(zip(head, ln.rstrip("\n").split("\t"))) for ln in fh)}
    lock = json.loads(lock_path.read_text()) if lock_path.exists() else {"genomes": {}}
    lock["genomes"] = {k: v for k, v in lock["genomes"].items()
                       if (out / v["file_name"]).exists() and (out / v["file_name"]).stat().st_size == v["bytes"]}
    R = json.loads(res_path.read_text()) if res_path.exists() else {"tsv": str(tsv), "started_utc": now()}
    R.update({"failed": {}, "fallback_single": R.get("fallback_single", []), "fallback_ftp": R.get("fallback_ftp", []),
              "md5_refetched": R.get("md5_refetched", [])})

    def save():
        write_json(lock_path, lock)
        write_json(res_path, {**R, **{"http": F.R}})

    def record(acc: str, path: Path, sc: dict, url: str, md5_source: str) -> None:
        lock["genomes"][acc] = {"file_name": path.name, "file": str(path), "sha256": sc["sha256"],
                                "md5": sc["md5"], "md5_source": md5_source, "bytes": sc["bytes"],
                                "records": sc["records"], "bp": sc["bp"], "url": url, "retrieved_utc": now(),
                                "gtdb_genome_size": int(want[acc]["genome_size"])}
        F.line(stage, acc, "ok", sc["bytes"], sc["bp"], sc["records"], sc["sha256"][:16], md5_source,
               "bp_eq_gtdb" if sc["bp"] == int(want[acc]["genome_size"]) else "bp_NE_gtdb")

    def via_api(accs: list[str]) -> tuple[set[str], set[str]]:
        """(verified, md5-bad) accessions of one Datasets request."""
        url = API.format(",".join(accs))
        zpath = tmpd / "batch.zip"
        try:
            F.get(url, zpath)
            got = extract_zip(zpath, set(accs), out)
        finally:
            zpath.unlink(missing_ok=True)
        ok, bad = set(), set()
        for acc, g in got.items():
            if g["ok"]:
                record(acc, g["file"], g["scan"], url, "datasets-zip md5sum.txt")
                ok.add(acc)
            else:
                bad.add(acc)
                F.line(stage, acc, "md5_mismatch", g["md5_expected"], g["scan"]["md5"])
        return ok, bad

    def via_ftp(acc: str) -> None:
        pre, num = acc.split("_")[0], acc.split("_")[1].split(".")[0]
        d = FTP.format(pre, num[0:3], num[3:6], num[6:9])
        listing = F.get(d).decode("utf-8", "replace")
        names = sorted(set(re.findall(r'href="(' + re.escape(acc) + r'_[^"/]+)/?"', listing)))
        if not names:
            raise RuntimeError(f"{d}: no directory for {acc}")
        base = f"{d}{names[0]}/"
        md5s = md5_table(F.get(base + "md5checksums.txt").decode("utf-8", "replace"))
        fname = f"{names[0]}_genomic.fna.gz"
        for attempt in range(2):
            gz = tmpd / f"{acc}.fna.gz"
            F.get(base + fname, gz)
            m = hashlib.md5()
            with open(gz, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    m.update(chunk)
            if md5s.get(fname) == m.hexdigest():
                break
            F.line(stage, acc, "md5_mismatch_ftp", md5s.get(fname), m.hexdigest(), f"attempt={attempt + 1}")
            gz.unlink()
        else:
            raise RuntimeError(f"{base}{fname}: md5 differs from md5checksums.txt twice")
        dest = out / f"{acc}.fna"
        part = dest.with_name(dest.name + ".part")
        with gzip.open(gz, "rb") as fi, open(part, "wb") as fo:
            shutil.copyfileobj(fi, fo, 1 << 20)
        gz.unlink()
        sc = scan_fna(part)
        part.replace(dest)
        record(acc, dest, sc, base + fname, "ftp md5checksums.txt (of .fna.gz)")

    todo = [acc for acc in sorted(want) if acc not in lock["genomes"]]
    F.line(stage, "start", f"wanted={len(want)}", f"done={len(want) - len(todo)}", f"todo={len(todo)}")
    for i in range(0, len(todo), batch):
        if STOP.exists():
            save()
            raise Stopped(f"STOP present at {STOP}")
        accs = todo[i : i + batch]
        ok, bad = set(), set()
        try:
            ok, bad = via_api(accs)
        except Refused:
            save()
            raise
        except Exception as e:  # noqa: BLE001
            F.line(stage, "batch", i, "error", str(e)[:200])
        for acc in accs:
            if acc in ok:
                continue
            try:
                ok1, bad1 = via_api([acc])
                if acc in ok1:
                    (R["md5_refetched"] if acc in bad else R["fallback_single"]).append(acc)
                    continue
                if acc in bad and acc in bad1:
                    R["failed"][acc] = "md5 differs from md5sum.txt twice"
                    F.line(stage, acc, "FAILED", "md5 twice")
                    continue
            except Refused:
                save()
                raise
            except Exception as e:  # noqa: BLE001
                F.line(stage, acc, "single_error", str(e)[:200])
            if not ftp_fallback:
                R["failed"][acc] = "not returned by Datasets v2 (alone or in its batch); FTP fallback off"
                F.line(stage, acc, "FAILED", "not returned by Datasets; FTP fallback off")
                continue
            try:
                via_ftp(acc)
                R["fallback_ftp"].append(acc)
            except Refused:
                save()
                raise
            except Exception as e:  # noqa: BLE001
                R["failed"][acc] = str(e)[:300]
                F.line(stage, acc, "FAILED", str(e)[:200])
        save()
    g = lock["genomes"]
    R.update({"finished_utc": now(), "genomes_wanted": len(want), "genomes_ok": len(g),
              "genomes_failed": len(R["failed"]), "bp": sum(v["bp"] for v in g.values()),
              "records": sum(v["records"] for v in g.values()), "bytes": sum(v["bytes"] for v in g.values()),
              "bp_vs_gtdb_mismatch": {k: [v["bp"], v["gtdb_genome_size"]] for k, v in g.items()
                                      if v["bp"] != v["gtdb_genome_size"]}})
    save()
    F.line(stage, "end", f"ok={len(g)}", f"failed={len(R['failed'])}", f"bytes={R['bytes']}")
    return {k: R[k] for k in ("genomes_wanted", "genomes_ok", "genomes_failed", "bp", "bytes")}


def selftest(zpath: Path, out: Path) -> int:
    out.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zpath) as z:
        accs = {m.group(1) for m in (FNA.fullmatch(n) for n in z.namelist()) if m}
    got = extract_zip(zpath, accs, out)
    print(json.dumps({a: {"ok": g["ok"], "md5": g["scan"]["md5"], "bp": g["scan"]["bp"]} for a, g in sorted(got.items())}))
    return 0 if got and all(g["ok"] for g in got.values()) else 1


def main(argv: list[str]) -> int:
    if argv[:1] == ["selftest"]:
        return selftest(Path(argv[1]), Path(argv[2]))
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("selection", type=Path)
    ap.add_argument("genomes", type=Path)
    ap.add_argument("log", type=Path)
    ap.add_argument("--stages", default="P3,P4,P5,P6")
    ap.add_argument("--batch", type=int, default=20)
    ap.add_argument("--ftp-fallback", action="store_true",
                    help="use the FTP fallback; it lists the assembly's parent directory, a URL form sources.yaml "
                         "does not list, so it is off by default")
    a = ap.parse_args(argv)
    a.genomes.mkdir(parents=True, exist_ok=True)
    marks = a.log.parent
    for marker in ("DONE", "FAILED"):
        (marks / marker).unlink(missing_ok=True)
    F = Fetcher(a.log)
    F.line("run", "start", f"stages={a.stages}")
    summary = {"started_utc": now(), "stages": {}}
    try:
        for stage in a.stages.split(","):
            summary["stages"][stage] = run_stage(F, stage, a.selection / f"stage_{stage}.tsv", a.genomes / stage, a.batch,
                                                 a.ftp_fallback)
    except (Refused, Stopped, Exception) as e:  # noqa: BLE001
        summary.update({"finished_utc": now(), "stopped": f"{type(e).__name__}: {e}", "http": F.R})
        write_json(marks / "FAILED", summary)
        F.line("run", "FAILED", f"{type(e).__name__}: {e}")
        return 75 if isinstance(e, Stopped) else 1
    summary.update({"finished_utc": now(), "http": F.R})
    write_json(marks / "DONE", summary)
    F.line("run", "DONE", json.dumps({s: [v["genomes_ok"], v["genomes_failed"]] for s, v in summary["stages"].items()}))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
