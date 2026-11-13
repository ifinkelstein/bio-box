#!/usr/bin/env python3
"""Card C-27 (3): fetch the genomes of S2C phase 2's ANI step (host python3, stdlib only).

    python3 -I c27_fetch.py split ANI_FETCH_TSV BATCH_DIR [--batch-gbp 10]
    python3 -I c27_fetch.py fetch BATCH_TSV OUT_DIR LOG [--batch 20] [--max-fail-frac 0.05]
    python3 -I c27_fetch.py need BATCH_TSV          # KiB of disk the batch is expected to take (0.30 bytes per bp)
    python3 -I c27_fetch.py selftest ZIP OUT_DIR

Phase 1's fetch (r50_fetch.py, whose Fetcher, URL and limits this imports unchanged): the NCBI Datasets v2
genome download listed in sources.yaml (`ncbi_datasets_genome`), 20 accessions per request, one connection, at
most 3 requests per second, the same generic user agent, no key, no email; an accession missing from its batch is
asked alone once; every FASTA is checked against the md5sum.txt inside the returned zip, and a file that fails is
asked once more alone, then listed as failed. The FTP fallback is not used (its directory listing is not a
listed URL form). Five refusals in a row (HTTP 403, 429, 503) stop the run; so does the data root's STOP file,
checked before every request batch.

One difference: the FASTA is kept compressed. Each zip member's deflate stream is copied unchanged into a gzip
file, OUT_DIR/<acc>.fna.gz (gzip header with mtime 0, the member's CRC-32 and size as the trailer; a member
stored without compression is gzipped at level 6 with mtime 0), so no CPU is spent recompressing and the bytes
on disk are about the bytes downloaded. The .fna.gz is then read back through gzip (which checks the CRC and
size) and its decompressed bytes are checked against md5sum.txt; lock.json records, as phase 1 did, the sha256,
md5, size, records and bp of the FASTA, plus gz_bytes and gz_sha256 of the file kept.

split: the rows of ani_fetch.tsv (c25_select.py --candidates) whose fetch_reason is phase1_not_served_retry go
to batch B00 (DEC-31: retried once); the rest, sorted by accession, are split as phase 1's stages were: n =
max(1, round(total_bp / batch_bp)); the genome whose predecessors sum to c bp goes to batch B(1 + floor(n c /
total_bp)). Writes BATCH_DIR/batch_B<kk>.tsv (the input's columns) and BATCH_DIR/split.json.

fetch: one batch into OUT_DIR, with lock.json and fetch_result.json (asked, verified, failed, bytes downloaded,
bytes kept, wall time). Resumable: a genome in lock.json whose .fna.gz has the recorded size is not asked again.
At the end OUT_DIR/DONE (exit 0) or OUT_DIR/FAILED (exit 1 refused or error, 75 STOP, 3 when more than
--max-fail-frac of the batch was not served: the chain of batches stops rather than keep asking).
"""

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import struct
import sys
import time
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python3 -I leaves the script's directory off the path
from r50_fetch import API, FNA, STOP, Fetcher, Refused, Stopped, md5_table, now, write_json  # noqa: E402

GZ_HEADER = b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\xff"  # deflate, no flags, mtime 0, no extra flags, OS unknown
BYTES_PER_BP = 0.30  # disk per bp kept compressed; derived: phase 1 downloaded 0.29026 bytes per bp (r50 record)
RETRY = "phase1_not_served_retry"


class TooManyFailed(RuntimeError):
    pass


def _raw_member(zpath: Path, info: zipfile.ZipInfo):
    """Yield the compressed bytes of one zip member, as stored."""
    with open(zpath, "rb") as fh:
        fh.seek(info.header_offset)
        hdr = fh.read(30)
        if hdr[:4] != b"PK\x03\x04":
            raise ValueError(f"{zpath}: no local header for {info.filename}")
        n_name, n_extra = struct.unpack("<HH", hdr[26:30])
        fh.seek(info.header_offset + 30 + n_name + n_extra)
        left = info.compress_size
        while left:
            chunk = fh.read(min(left, 1 << 20))
            if not chunk:
                raise ValueError(f"{zpath}: {info.filename} is truncated")
            left -= len(chunk)
            yield chunk


def write_gz(zpath: Path, z: zipfile.ZipFile, info: zipfile.ZipInfo, dest: Path) -> None:
    with open(dest, "wb") as fo:
        if info.compress_type == zipfile.ZIP_DEFLATED and not info.flag_bits & 1:
            fo.write(GZ_HEADER)
            for chunk in _raw_member(zpath, info):
                fo.write(chunk)
            fo.write(struct.pack("<II", info.CRC & 0xFFFFFFFF, info.file_size & 0xFFFFFFFF))
        else:
            with z.open(info) as fi, gzip.GzipFile(filename="", mode="wb", fileobj=fo, mtime=0, compresslevel=6) as go:
                shutil.copyfileobj(fi, go, 1 << 20)


def scan_gz(path: Path) -> dict:
    """sha256, md5, size, records and bp of the FASTA inside a .fna.gz (gzip checks CRC and size), and the file's own
    size and sha256."""
    h, m = hashlib.sha256(), hashlib.md5()
    n_rec = bp = size = 0
    with gzip.open(path, "rb") as fh:
        for ln in fh:
            h.update(ln)
            m.update(ln)
            size += len(ln)
            if ln.startswith(b">"):
                n_rec += 1
            else:
                bp += len(ln.strip())
    g = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            g.update(chunk)
    return {"sha256": h.hexdigest(), "md5": m.hexdigest(), "bytes": size, "records": n_rec, "bp": bp,
            "gz_bytes": path.stat().st_size, "gz_sha256": g.hexdigest()}


def extract_zip_gz(zpath: Path, accs: set[str], dest_dir: Path) -> dict[str, dict]:
    """Each wanted accession's FASTA in a Datasets zip -> dest_dir/<acc>.fna.gz, kept only if its md5 matches.

    Returns {acc: {"file", "md5_expected", "scan", "ok", "zip_member"}} for every wanted accession in the zip."""
    got = {}
    with zipfile.ZipFile(zpath) as z:
        names = z.namelist()
        md5s = md5_table(z.read("md5sum.txt").decode()) if "md5sum.txt" in names else {}
        for info in z.infolist():
            mt = FNA.fullmatch(info.filename)
            if not mt or mt.group(1) not in accs:
                continue
            acc = mt.group(1)
            dest = dest_dir / f"{acc}.fna.gz"
            part = dest.with_name(dest.name + ".part")
            try:
                write_gz(zpath, z, info, part)
                sc = scan_gz(part)
            except (OSError, ValueError, EOFError) as e:  # gzip.BadGzipFile is an OSError
                part.unlink(missing_ok=True)
                sc = {"error": f"{type(e).__name__}: {e}", "md5": None}
            exp = md5s.get(info.filename)
            ok = exp is not None and exp == sc["md5"]
            if ok:
                part.replace(dest)
            else:
                part.unlink(missing_ok=True)
            got[acc] = {"file": dest, "md5_expected": exp, "scan": sc, "ok": ok, "zip_member": info.filename}
    return got


def read_tsv(path: Path) -> tuple[list[str], list[dict]]:
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        return list(rd.fieldnames or []), list(rd)


def run_batch(F, name: str, tsv: Path, out: Path, batch: int = 20, stop: Path = STOP,
              max_fail_frac: float = 1.0) -> dict:
    """Fetch one batch's genomes into `out` (r50_fetch.run_stage's logic, the FASTA kept as .fna.gz)."""
    t0 = time.time()
    out.mkdir(parents=True, exist_ok=True)
    tmpd = out / ".tmp"
    tmpd.mkdir(exist_ok=True)
    lock_path, res_path = out / "lock.json", out / "fetch_result.json"
    want = {r["accession"]: r for r in read_tsv(tsv)[1]}
    lock = json.loads(lock_path.read_text()) if lock_path.exists() else {"genomes": {}}
    lock["genomes"] = {k: v for k, v in lock["genomes"].items()
                       if (out / v["file_name"]).exists() and (out / v["file_name"]).stat().st_size == v["gz_bytes"]}
    R = json.loads(res_path.read_text()) if res_path.exists() else {"tsv": str(tsv), "started_utc": now()}
    R.update({"failed": {}, "fallback_single": R.get("fallback_single", []), "md5_refetched": R.get("md5_refetched", [])})
    http0 = dict(F.R)

    def save():
        write_json(lock_path, lock)
        write_json(res_path, {**R, "http": F.R})

    def record(acc: str, path: Path, sc: dict, url: str) -> None:
        lock["genomes"][acc] = {"file_name": path.name, "file": str(path), "sha256": sc["sha256"], "md5": sc["md5"],
                                "md5_source": "datasets-zip md5sum.txt", "bytes": sc["bytes"],
                                "records": sc["records"], "bp": sc["bp"], "gz_bytes": sc["gz_bytes"],
                                "gz_sha256": sc["gz_sha256"], "url": url, "retrieved_utc": now(),
                                "gtdb_genome_size": int(want[acc]["genome_size"])}
        F.line(name, acc, "ok", sc["bytes"], sc["gz_bytes"], sc["bp"], sc["records"], sc["sha256"][:16],
               "bp_eq_gtdb" if sc["bp"] == int(want[acc]["genome_size"]) else "bp_NE_gtdb")

    def via_api(accs: list[str]) -> tuple[set[str], set[str]]:
        url = API.format(",".join(accs))
        zpath = tmpd / "batch.zip"
        try:
            F.get(url, zpath)
            got = extract_zip_gz(zpath, set(accs), out)
        finally:
            zpath.unlink(missing_ok=True)
        ok, bad = set(), set()
        for acc, g in got.items():
            if g["ok"]:
                record(acc, g["file"], g["scan"], url)
                ok.add(acc)
            else:
                bad.add(acc)
                F.line(name, acc, "md5_mismatch", g["md5_expected"], g["scan"].get("md5"), g["scan"].get("error", ""))
        return ok, bad

    todo = [acc for acc in sorted(want) if acc not in lock["genomes"]]
    F.line(name, "start", f"wanted={len(want)}", f"done={len(want) - len(todo)}", f"todo={len(todo)}")
    for i in range(0, len(todo), batch):
        if stop.exists():
            save()
            raise Stopped(f"STOP present at {stop}")
        accs = todo[i:i + batch]
        ok, bad = set(), set()
        try:
            ok, bad = via_api(accs)
        except Refused:
            save()
            raise
        except Exception as e:  # noqa: BLE001
            F.line(name, "batch", i, "error", str(e)[:200])
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
                    F.line(name, acc, "FAILED", "md5 twice")
                    continue
            except Refused:
                save()
                raise
            except Exception as e:  # noqa: BLE001
                F.line(name, acc, "single_error", str(e)[:200])
            R["failed"][acc] = "not returned by Datasets v2 (alone or in its batch); FTP fallback off"
            F.line(name, acc, "FAILED", "not returned by Datasets")
        save()
    shutil.rmtree(tmpd, ignore_errors=True)
    g = lock["genomes"]
    R.update({"finished_utc": now(), "genomes_wanted": len(want), "genomes_asked_this_run": len(todo),
              "genomes_ok": len(g), "genomes_failed": len(R["failed"]),
              "bp": sum(v["bp"] for v in g.values()), "records": sum(v["records"] for v in g.values()),
              "bytes": sum(v["bytes"] for v in g.values()), "gz_bytes": sum(v["gz_bytes"] for v in g.values()),
              "requests_this_run": F.R["requests"] - http0.get("requests", 0),
              "bytes_downloaded_this_run": F.R["bytes_downloaded"] - http0.get("bytes_downloaded", 0),
              "wall_s_this_run": round(time.time() - t0, 1),
              "bp_vs_gtdb_mismatch": {k: [v["bp"], v["gtdb_genome_size"]] for k, v in g.items()
                                      if v["bp"] != v["gtdb_genome_size"]}})
    save()
    F.line(name, "end", f"ok={len(g)}", f"failed={len(R['failed'])}", f"gz_bytes={R['gz_bytes']}")
    res = {k: R[k] for k in ("genomes_wanted", "genomes_ok", "genomes_failed", "bp", "bytes", "gz_bytes",
                             "requests_this_run", "bytes_downloaded_this_run", "wall_s_this_run")}
    if want and len(R["failed"]) > max_fail_frac * len(want):
        raise TooManyFailed(f"{len(R['failed'])} of {len(want)} genomes not served (over {max_fail_frac:g})")
    return res


def split(src: Path, out: Path, batch_bp: int = 10**10) -> dict:
    cols, rows = read_tsv(src)
    retry = sorted((r for r in rows if r.get("fetch_reason") == RETRY), key=lambda r: r["accession"])
    rest = sorted((r for r in rows if r.get("fetch_reason") != RETRY), key=lambda r: r["accession"])
    batches: dict[str, list[dict]] = {}
    if retry:
        batches["B00"] = retry
    total = sum(int(r["genome_size"]) for r in rest)
    n = max(1, round(total / batch_bp)) if total else 0
    c = 0
    for r in rest:
        batches.setdefault(f"B{1 + (n * c) // total:02d}", []).append(r)
        c += int(r["genome_size"])
    out.mkdir(parents=True, exist_ok=True)
    plan = {"source": src.name, "source_sha256": hashlib.sha256(src.read_bytes()).hexdigest(), "batch_bp": batch_bp,
            "n": n, "rule": "B00: the phase-1 retries; then B(1 + floor(n * bp_before / total_bp)) over the rest "
                            "sorted by accession, n = max(1, round(total_bp / batch_bp))", "batches": {}}
    for k, rs in sorted(batches.items()):
        text = "\t".join(cols) + "\n" + "".join("\t".join(r[c] for c in cols) + "\n" for r in rs)
        p = out / f"batch_{k}.tsv"
        p.with_name(f".{p.name}.tmp").write_text(text)
        p.with_name(f".{p.name}.tmp").replace(p)
        plan["batches"][k] = {"genomes": len(rs), "bp": sum(int(r["genome_size"]) for r in rs),
                              "first": rs[0]["accession"], "last": rs[-1]["accession"],
                              "sha256": hashlib.sha256(text.encode()).hexdigest()}
    write_json(out / "split.json", plan)
    return plan


def main(argv: list[str]) -> int:
    if argv[:1] == ["selftest"]:
        zpath, out = Path(argv[1]), Path(argv[2])
        out.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zpath) as z:
            accs = {m.group(1) for m in (FNA.fullmatch(n) for n in z.namelist()) if m}
        got = extract_zip_gz(zpath, accs, out)
        print(json.dumps({a: {"ok": g["ok"], "md5": g["scan"].get("md5"), "bp": g["scan"].get("bp"),
                              "gz_bytes": g["scan"].get("gz_bytes")} for a, g in sorted(got.items())}))
        return 0 if got and all(g["ok"] for g in got.values()) else 1
    if argv[:1] == ["need"]:
        print(int(sum(int(r["genome_size"]) for r in read_tsv(Path(argv[1]))[1]) * BYTES_PER_BP / 1024) + 1)
        return 0
    if argv[:1] == ["split"]:
        ap = argparse.ArgumentParser(prog="c27_fetch.py split")
        ap.add_argument("src", type=Path)
        ap.add_argument("out", type=Path)
        ap.add_argument("--batch-gbp", type=float, default=10.0)
        a = ap.parse_args(argv[1:])
        print(json.dumps(split(a.src, a.out, int(a.batch_gbp * 10**9)), indent=1))
        return 0
    ap = argparse.ArgumentParser(prog="c27_fetch.py fetch", description=__doc__.split("\n\n")[0])
    ap.add_argument("cmd", choices=["fetch"])
    ap.add_argument("tsv", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("log", type=Path)
    ap.add_argument("--batch", type=int, default=20)
    ap.add_argument("--max-fail-frac", type=float, default=0.05)
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    for marker in ("DONE", "FAILED"):
        (a.out / marker).unlink(missing_ok=True)
    F = Fetcher(a.log)
    name = a.out.name
    F.line(name, "run", "start", a.tsv)
    try:
        res = run_batch(F, name, a.tsv, a.out, a.batch, STOP, a.max_fail_frac)
    except (Refused, Stopped, TooManyFailed, Exception) as e:  # noqa: BLE001
        write_json(a.out / "FAILED", {"finished_utc": now(), "stopped": f"{type(e).__name__}: {e}", "http": F.R})
        F.line(name, "run", "FAILED", f"{type(e).__name__}: {e}")
        return 75 if isinstance(e, Stopped) else 3 if isinstance(e, TooManyFailed) else 1
    write_json(a.out / "DONE", {"finished_utc": now(), **res, "http": F.R})
    F.line(name, "run", "DONE", json.dumps(res))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
