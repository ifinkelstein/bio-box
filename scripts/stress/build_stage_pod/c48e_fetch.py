"""Card C-48e part 1: fetch a tranche's SPIRE assemblies politely, resumably and verifiably.

    python c48e_fetch.py --selection /tmp/if2453/c48e1/select/selection.tsv --out-dir /tmp/if2453/c48e1/dl \
        [--interval 1.0] [--tries 3] [--max-dl-gb 150] [--min-free-gb 400] \
        [--stop-file /tmp/if2453/r10/root/STOP --stop-file /tmp/if2453/c48e1/STOP]

Rules (AGENTS §8, DEC-50, `sources.yaml` entry `spire_tranche_1`, written before the first request):
- Only `https://spire.embl.de/download_assembly/{sample_id}` (and the provider's own redirect); User-Agent `curl`,
  as the pilot sent (card C-42); one connection; at least `--interval` seconds (never under 0.5) between request
  starts.
- A sample is tried at most `--tries` times on a transient error (timeout, broken or short stream, a gzip stream
  that does not decompress whole, HTTP 500/502/503/504); then it is listed as `failed`, never replaced. HTTP 404 is
  recorded as `not_found` and not retried. HTTP 403 or 429 stops the whole fetch (exit 4) for a report.
- Before each sample: a `STOP` file stops (exit 75); `/tmp` under `--min-free-gb` free or the fetched files over
  `--max-dl-gb` stops (exit 3).
- Each file is written as `.<name>.part`, fsynced and renamed. A plain FASTA answer is stored gzipped (mtime 0).
  Every outcome is a line of `<out-dir>/fetch_log.jsonl` (sizes, md5 and sha256 on the wire and on disk, the
  decompressed bytes, times, the try). A rerun skips a sample whose last line is `ok` and whose file still has the
  logged size and sha256, or whose last line is `not_found`; any other sample is fetched again.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import shutil
import socket
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zlib
from pathlib import Path

URL = "https://spire.embl.de/download_assembly/{sample_id}"
UA = {"User-Agent": "curl"}  # the pilot's client string (c42_fetch.UA): no name, address or project detail
STOP_CODES = {403, 429}
_CHUNK = 1 << 20


class Stop(Exception):
    def __init__(self, code: int, why: str):
        super().__init__(why)
        self.code = code


class _Retry(Exception):
    pass


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while b := fh.read(_CHUNK):
            h.update(b)
    return h.hexdigest()


def read_log(path: Path) -> dict[str, dict]:
    """The last final line (`ok`, `not_found`, `failed`, `refused`) of each sample."""
    last: dict[str, dict] = {}
    if path.is_file():
        for ln in path.read_text().splitlines():
            if ln.strip():
                r = json.loads(ln)
                if r.get("status") in ("ok", "not_found", "failed", "refused"):
                    last[r["sample_id"]] = r
    return last


def verified(rec: dict | None, dst: Path) -> bool:
    return (rec is not None and rec.get("status") == "ok" and dst.is_file()
            and dst.stat().st_size == rec["disk_bytes"] and file_sha256(dst) == rec["disk_sha256"])


class _Gunzip:
    """Decompress a gzip stream of one or more members, counting bytes; raises on a stream that does not end."""

    def __init__(self):
        self.d = zlib.decompressobj(16 + zlib.MAX_WBITS)
        self.started = False  # the current member has been fed bytes
        self.n = 0
        self.head = b""

    def feed(self, data: bytes) -> None:
        while data:
            self.started = True
            out = self.d.decompress(data)
            self.n += len(out)
            if not self.head and out:
                self.head = out[:1]
            if self.d.eof:  # a member ended; what follows starts the next one
                data = self.d.unused_data
                self.d = zlib.decompressobj(16 + zlib.MAX_WBITS)
                self.started = False
            else:
                data = b""

    def complete(self) -> bool:
        return not self.started


class _Pace:
    def __init__(self, interval: float):
        if interval < 0.5:
            raise ValueError(f"interval must be >= 0.5 s (at most 2 requests per second), got {interval}")
        self.interval = interval
        self.last = 0.0

    def wait(self) -> None:
        w = self.last + self.interval - time.monotonic()
        if w > 0:
            time.sleep(w)
        self.last = time.monotonic()


def _try_once(url: str, dst: Path, pace: _Pace, timeout: float) -> dict:
    """One request. Returns an `ok` or `not_found` record; raises `_Retry` or `Stop`."""
    part = dst.with_name(f".{dst.name}.part")
    pace.wait()
    t0 = time.monotonic()
    started = _now()
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as resp:
            final_url, status = resp.geturl(), resp.status
            clen = resp.headers.get("Content-Length")
            first = resp.read(_CHUNK)
            is_gz = first[:2] == b"\x1f\x8b"
            wmd5, wsha = hashlib.md5(), hashlib.sha256()
            wire = 0
            gun = _Gunzip() if is_gz else None
            raw_head = first[:1]
            with open(part, "wb") as out:
                gz = None if is_gz else gzip.GzipFile(filename="", mode="wb", fileobj=out, compresslevel=6, mtime=0)
                chunk = first
                while chunk:
                    wire += len(chunk)
                    wmd5.update(chunk)
                    wsha.update(chunk)
                    if gun is not None:
                        gun.feed(chunk)
                        out.write(chunk)
                    else:
                        gz.write(chunk)
                    chunk = resp.read(_CHUNK)
                if gz is not None:
                    gz.close()
                out.flush()
                os.fsync(out.fileno())
    except urllib.error.HTTPError as e:
        part.unlink(missing_ok=True)
        if e.code == 404:
            return {"status": "not_found", "http": 404, "started_utc": started, "fetched_utc": _now()}
        if e.code in STOP_CODES:
            raise Stop(4, f"HTTP {e.code} on {url}") from None
        raise _Retry(f"HTTP {e.code}") from None
    except (TimeoutError, socket.timeout, urllib.error.URLError, ConnectionError, zlib.error, EOFError, OSError) as e:
        part.unlink(missing_ok=True)
        raise _Retry(f"{type(e).__name__}: {e}") from None
    try:
        if clen is not None and int(clen) != wire:
            raise _Retry(f"short body: {wire} of {clen} bytes")
        if gun is not None:
            if not gun.complete():
                raise _Retry("gzip stream ends inside a member")
            unc, head = gun.n, gun.head
        else:
            unc, head = wire, raw_head
        if unc and head != b">":
            raise _Retry("not FASTA: the first byte is not '>'")
    except _Retry:
        part.unlink(missing_ok=True)
        raise
    dmd5, dsha = hashlib.md5(), hashlib.sha256()
    if is_gz:
        dmd5, dsha = wmd5, wsha
    else:
        with open(part, "rb") as fh:
            while b := fh.read(_CHUNK):
                dmd5.update(b)
                dsha.update(b)
    part.replace(dst)
    return {"status": "ok", "http": status, "final_url": final_url, "wire_gzipped": is_gz, "wire_bytes": wire,
            "disk_bytes": dst.stat().st_size, "uncompressed_bytes": unc, "wire_md5": wmd5.hexdigest(),
            "wire_sha256": wsha.hexdigest(), "disk_md5": dmd5.hexdigest(), "disk_sha256": dsha.hexdigest(),
            "started_utc": started, "fetched_utc": _now(), "seconds": round(time.monotonic() - t0, 3)}


def _dl_bytes(out: Path) -> int:
    return sum(e.stat().st_size for e in os.scandir(out) if e.is_file())


def fetch_all(sample_ids: list[str], out_dir: Path, *, url_template: str = URL, interval: float = 1.0, tries: int = 3,
              retry_wait: float = 15.0, max_dl_bytes: int = 150 * 10**9, min_free_bytes: int = 400 * 10**9,
              stop_files: list[Path] = (), timeout: float = 180.0) -> dict:
    """Fetch `sample_ids` in order into `out_dir` (module docstring). Returns the lists by outcome (sorted)."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    log_path = out / "fetch_log.jsonl"
    last = read_log(log_path)
    pace = _Pace(interval)
    res = {"ok": [], "not_found": [], "failed": [], "skipped_verified": []}

    def log(rec: dict) -> None:
        with open(log_path, "a") as fh:
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
            fh.flush()
            os.fsync(fh.fileno())

    for sid in sample_ids:
        dst = out / f"{sid}.fa.gz"
        prev = last.get(sid)
        if prev is not None and prev["status"] == "not_found":
            res["not_found"].append(sid)
            continue
        if verified(prev, dst):
            res["skipped_verified"].append(sid)
            continue
        for p in stop_files:
            if Path(p).exists():
                raise Stop(75, f"STOP present at {p}")
        free = shutil.disk_usage(out).free
        if free < min_free_bytes:
            raise Stop(3, f"{free} bytes free under {out}, under the floor of {min_free_bytes}")
        used = _dl_bytes(out)
        if used >= max_dl_bytes:
            raise Stop(3, f"{used} bytes fetched under {out}, at or over the cap of {max_dl_bytes}")
        if dst.exists():
            dst.unlink()  # its bytes do not match the log: fetched again
        url = url_template.format(sample_id=urllib.parse.quote(sid, safe=""))
        why = []
        rec = None
        for attempt in range(1, tries + 1):
            try:
                rec = _try_once(url, dst, pace, timeout)
                break
            except _Retry as e:
                why.append(str(e))
                log({"sample_id": sid, "url": url, "status": "retry", "attempt": attempt, "why": str(e), "utc": _now()})
                if attempt < tries:
                    time.sleep(retry_wait * attempt)
            except Stop as e:
                log({"sample_id": sid, "url": url, "status": "refused", "attempt": attempt, "why": str(e),
                     "utc": _now()})
                raise
        if rec is None:
            rec = {"status": "failed", "why": why, "fetched_utc": _now()}
            attempt = tries
        rec.update(sample_id=sid, url=url, attempt=attempt)
        log(rec)
        res[rec["status"]].append(sid)
    return {k: sorted(v) for k, v in res.items()}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--interval", type=float, default=1.0)
    ap.add_argument("--tries", type=int, default=3)
    ap.add_argument("--retry-wait", type=float, default=15.0)
    ap.add_argument("--max-dl-gb", type=float, default=150.0)
    ap.add_argument("--min-free-gb", type=float, default=400.0)
    ap.add_argument("--stop-file", action="append", default=[])
    a = ap.parse_args(argv)
    with open(a.selection, newline="") as fh:
        sids = [r["sample_id"] for r in csv.DictReader(fh, delimiter="\t")]
    try:
        res = fetch_all(sids, Path(a.out_dir), interval=a.interval, tries=a.tries, retry_wait=a.retry_wait,
                        max_dl_bytes=int(a.max_dl_gb * 1e9), min_free_bytes=int(a.min_free_gb * 1e9),
                        stop_files=[Path(p) for p in a.stop_file])
    except Stop as e:
        print(f"{_now()} stopped: {e}", file=sys.stderr, flush=True)
        return e.code
    summary = {k: len(v) for k, v in res.items()} | {"not_found_ids": res["not_found"], "failed_ids": res["failed"],
                                                      "requested": len(sids), "done_utc": _now()}
    (Path(a.out_dir) / "fetch_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
