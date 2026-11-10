"""Count `*` in `protein_seq` (card W-11): sequences ending in `*` and sequences with an internal `*`.

Read only. The file list comes from `manifests/release=<RELEASE>.json` under ROOT (AGENTS 3.9). Streams
one record batch at a time with pyarrow, one thread. Prints one JSON object on stdout.

    python scripts/stress/protein_stars.py ROOT RELEASE [--gene-call-version V]
"""
import argparse
import json
from pathlib import Path

import pyarrow.compute as pc
import pyarrow.parquet as pq


def count(root: Path, release: str, gene_call_version: str | None) -> dict:
    manifest = json.loads((root / f"manifests/release={release}.json").read_text())
    files = sorted(f["path"] for f in manifest["files"] if f["table"] == "protein_seq"
                   and (gene_call_version is None or f["partitions"].get("gene_call_version") == gene_call_version))
    out = {"release": release, "files": len(files), "sequences": 0, "ending_in_star": 0, "internal_star": 0,
           "only_star": 0, "any_star": 0}
    for rel in files:
        pf = pq.ParquetFile(root / rel)
        for batch in pf.iter_batches(batch_size=65536, columns=["seq"], use_threads=False):
            seq = batch.column("seq")
            ends = pc.ends_with(seq, "*")
            out["sequences"] += len(seq)
            out["ending_in_star"] += pc.sum(ends).as_py() or 0
            out["any_star"] += pc.sum(pc.match_substring(seq, "*")).as_py() or 0
            stars = pc.count_substring(seq, "*")
            out["internal_star"] += pc.sum(pc.greater(pc.subtract(stars, ends.cast("int64")), 0)).as_py() or 0
            out["only_star"] += pc.sum(pc.equal(seq, "*")).as_py() or 0
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("root", type=Path)
    ap.add_argument("release")
    ap.add_argument("--gene-call-version", default=None)
    a = ap.parse_args()
    print(json.dumps(count(a.root, a.release, a.gene_call_version), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
