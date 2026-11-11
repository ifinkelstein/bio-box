#!/usr/bin/env python3
"""Check the import ladder's fetched genomes before they are reused in place (host python3, stdlib only).

    python3 r50_verify_present.py LADDER_GENOMES_DIR OUT_JSON

For every entry of LADDER_GENOMES_DIR/{P1,P2}/lock.json: the file exists, its size, sha256, record count
and bp equal what the ladder recorded when it fetched it. Writes OUT_JSON (counts and every mismatch) and
exits 1 when any file does not match.
"""

import hashlib
import json
import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    root, out = Path(argv[0]), Path(argv[1])
    res = {"stages": {}, "mismatch": {}}
    for stage in ("P1", "P2"):
        lock = json.loads((root / stage / "lock.json").read_text())["genomes"]
        n = nbytes = bp = 0
        for acc in sorted(lock):
            v = lock[acc]
            p = Path(v["file"])
            if not p.exists():
                res["mismatch"][acc] = "missing"
                continue
            h = hashlib.sha256()
            recs = b = 0
            with open(p, "rb") as fh:
                for ln in fh:
                    h.update(ln)
                    if ln.startswith(b">"):
                        recs += 1
                    else:
                        b += len(ln.strip())
            got = {"sha256": h.hexdigest(), "bytes": p.stat().st_size, "records": recs, "bp": b}
            bad = {k: [got[k], v[k]] for k in got if got[k] != v[k]}
            if bad:
                res["mismatch"][acc] = bad
            n += 1
            nbytes += got["bytes"]
            bp += b
        res["stages"][stage] = {"genomes": len(lock), "checked": n, "bytes": nbytes, "bp": bp}
    res["ok"] = not res["mismatch"]
    tmp = out.with_name(f".{out.name}.tmp")
    tmp.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    tmp.replace(out)
    print(json.dumps(res, sort_keys=True))
    return 0 if res["ok"] else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
