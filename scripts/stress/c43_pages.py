"""Fetch many small pages of a listed source, politely, and log each one (card C-43).

    python -I c43_pages.py list URLS_TSV OUT_DIR --interval S
    python -I c43_pages.py next START_URL OUT_DIR --interval S [--max-pages N]

`list`: URLS_TSV has two columns, url and a file name; each URL is fetched once into
OUT_DIR/<name> (an existing file is skipped, so a stopped run resumes).
`next`: a JSON:API listing; page k is written to OUT_DIR/page_<k:06d>.json and the next URL
is the page's `links.next`, until it is null. A resumed run starts after the last page on disk.

One connection; at least --interval seconds between request starts; Python's default
User-Agent, no key, no name. A refusal (HTTP 403, 429, 500, 502, 503, 504, or a timeout) is
retried after 10 s times its count in a row; five in a row stop the run with exit code 4.
Every file is written under a temporary name and renamed. OUT_DIR/fetch_log.jsonl gets one
line per file: url, file, bytes, sha256, pod clock (UTC). OUT_DIR/DONE is written at the end.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import socket
import sys
import time
import urllib.error
import urllib.request

REFUSALS = {403, 429, 500, 502, 503, 504}
MAX_IN_A_ROW = 5


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


class Fetcher:
    def __init__(self, interval: float):
        self.interval = interval
        self.last = 0.0

    def get(self, url: str) -> bytes:
        in_a_row = 0
        while True:
            wait = self.last + self.interval - time.monotonic()
            if wait > 0:
                time.sleep(wait)
            self.last = time.monotonic()
            try:
                with urllib.request.urlopen(url, timeout=120) as r:
                    return r.read()
            except urllib.error.HTTPError as e:
                if e.code not in REFUSALS:
                    raise
                why = f"HTTP {e.code}"
            except (TimeoutError, socket.timeout, urllib.error.URLError) as e:
                why = f"{type(e).__name__}: {e}"
            in_a_row += 1
            print(f"{_now()} refusal {in_a_row} in a row ({why}) for {url}", file=sys.stderr, flush=True)
            if in_a_row >= MAX_IN_A_ROW:
                raise SystemExit(4)
            time.sleep(10 * in_a_row)


def _write(out_dir: str, name: str, url: str, body: bytes) -> None:
    path = os.path.join(out_dir, name)
    with open(path + ".part", "wb") as f:
        f.write(body)
    os.replace(path + ".part", path)
    rec = {"url": url, "file": name, "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest(),
           "fetched_utc": _now()}
    with open(os.path.join(out_dir, "fetch_log.jsonl"), "a") as log:
        log.write(json.dumps(rec) + "\n")


def by_list(urls_tsv: str, out_dir: str, f: Fetcher) -> None:
    with open(urls_tsv) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    for url, name in rows:
        if os.path.exists(os.path.join(out_dir, name)):
            continue
        _write(out_dir, name, url, f.get(url))


def by_next(start: str, out_dir: str, f: Fetcher, max_pages: int | None) -> None:
    done = sorted(n for n in os.listdir(out_dir) if n.startswith("page_") and n.endswith(".json"))
    k, url = 1, start
    if done:
        last = done[-1]
        k = int(last[5:11]) + 1
        with open(os.path.join(out_dir, last)) as fh:
            url = json.load(fh)["links"]["next"]
    while url and (max_pages is None or k <= max_pages):
        body = f.get(url)
        _write(out_dir, f"page_{k:06d}.json", url, body)
        url = json.loads(body)["links"]["next"]
        k += 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["list", "next"])
    ap.add_argument("source")
    ap.add_argument("out_dir")
    ap.add_argument("--interval", type=float, required=True)
    ap.add_argument("--max-pages", type=int, default=None)
    a = ap.parse_args(argv)
    if a.interval < 0.34:
        print("refusing: more than 3 requests per second", file=sys.stderr)
        return 2
    os.makedirs(a.out_dir, exist_ok=True)
    f = Fetcher(a.interval)
    if a.mode == "list":
        by_list(a.source, a.out_dir, f)
    else:
        by_next(a.source, a.out_dir, f, a.max_pages)
    if a.max_pages is None:
        with open(os.path.join(a.out_dir, "DONE"), "w") as fh:
            fh.write(_now() + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
