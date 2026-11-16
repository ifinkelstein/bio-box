"""Fetch one table listed in sources.yaml and record what arrived (card C-43).

    python -I c43_fetch.py URL OUT_PATH [--md5 HEX] [--max-bytes N]

Tables only: the card fetches no sequence. One connection, no key, no name and no
identifying User-Agent (Python's own default header is sent). The file is written as
OUT_PATH.part and renamed into place only when complete and, if --md5 is given, only when
the provider's md5 matches; otherwise the part file is removed and the exit code is 1.
A response whose Content-Length exceeds --max-bytes (default 20 GB) is refused before
any byte is written. One JSON line is printed: url, path, bytes, sha256, md5, md5_ok,
the HTTP headers that date and size the file, and the pod's clock at the end (UTC).
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import sys
import urllib.request

LIMIT = 20 * 10**9


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("url")
    ap.add_argument("out")
    ap.add_argument("--md5", default=None)
    ap.add_argument("--max-bytes", type=int, default=LIMIT)
    a = ap.parse_args(argv)
    if os.path.exists(a.out):
        print(f"refusing: {a.out} exists", file=sys.stderr)
        return 2
    part = a.out + ".part"
    sha, md5, n = hashlib.sha256(), hashlib.md5(), 0
    with urllib.request.urlopen(a.url, timeout=120) as r:
        hdr = {k: r.headers.get(k) for k in ("Content-Length", "Last-Modified", "Content-Type", "ETag")}
        clen = hdr["Content-Length"]
        if clen is not None and int(clen) > a.max_bytes:
            print(json.dumps({"url": a.url, "refused": "too large", "content_length": int(clen)}))
            return 3
        with open(part, "wb") as f:
            while True:
                b = r.read(1 << 20)
                if not b:
                    break
                n += len(b)
                if n > a.max_bytes:
                    f.close()
                    os.remove(part)
                    print(json.dumps({"url": a.url, "refused": "too large", "bytes_so_far": n}))
                    return 3
                sha.update(b)
                md5.update(b)
                f.write(b)
            f.flush()
            os.fsync(f.fileno())
        status = r.status
    md5_ok = None if a.md5 is None else (md5.hexdigest() == a.md5.lower())
    rec = {"url": a.url, "path": a.out, "status": status, "bytes": n, "sha256": sha.hexdigest(),
           "md5": md5.hexdigest(), "md5_ok": md5_ok, "headers": hdr,
           "fetched_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
    if md5_ok is False:
        os.remove(part)
        print(json.dumps(rec))
        return 1
    os.replace(part, a.out)
    print(json.dumps(rec))
    return 0


if __name__ == "__main__":
    sys.exit(main())
