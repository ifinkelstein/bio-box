"""Fetch the golden set reproducibly and pin it in a lock file.

    python -m biobox.fetch golden.yaml OUT_DIR [--lock golden.lock.yaml]

Each genome is written to OUT_DIR/<id>/<id>.fna.gz (gzip with mtime 0), one
record per accession in `golden.yaml` order. Every record is checked against
the expected accession and length. The lock records sha256 of the
decompressed FASTA (portable across zlib versions) and of the gzip file. If a
lock already exists, the fetch is verified against it instead of rewriting it.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import sys
import time
import urllib.request
from datetime import date
from pathlib import Path

import yaml

from biobox import fasta

NCBI_DELAY_S = 0.4  # stay under NCBI's 3 requests/s without an API key


def _get(url: str, attempts: int = 4) -> bytes:
    for i in range(attempts):
        try:
            with urllib.request.urlopen(url, timeout=300) as r:
                return r.read()
        except (TimeoutError, OSError) as exc:
            if i == attempts - 1:
                raise
            wait = 5 * 2**i
            print(f"retry {i + 1}/{attempts - 1} in {wait}s after {exc!r}", file=sys.stderr)
            time.sleep(wait)
    raise AssertionError("unreachable")


def _records_from_bytes(raw: bytes, label: str) -> list[fasta.Record]:
    tmp = io.BytesIO(raw)
    records, header, lines = [], None, []
    for ln in tmp:
        ln = ln.rstrip(b"\r\n")
        if not ln.strip():
            continue
        if ln.startswith(b">"):
            if header is not None:
                records.append(fasta.Record(header, b"".join(lines), fasta._width(lines)))
            header, lines = ln[1:].decode(), []
        else:
            if header is None:
                raise ValueError(f"{label}: sequence before header")
            lines.append(ln)
    if header is not None:
        records.append(fasta.Record(header, b"".join(lines), fasta._width(lines)))
    return records


def _ncbi_record(url_tmpl: str, accession: str) -> fasta.Record:
    raw = _get(url_tmpl.format(accession=accession))
    time.sleep(NCBI_DELAY_S)
    recs = _records_from_bytes(raw, accession)
    if len(recs) != 1 or recs[0].accession != accession:
        raise ValueError(f"{accession}: expected one record with that accession, got {[r.accession for r in recs]}")
    return recs[0]


def _genome_records(entry: dict, cfg: dict) -> list[fasta.Record]:
    src = entry.get("source", {"fetch": "ncbi"})
    if src["fetch"] == "ncbi":
        recs = [_ncbi_record(cfg["ncbi"], r["accession"]) for r in entry["records"]]
    elif src["fetch"] == "sdr":
        raw = _get(cfg["sdr"].format(druid=src["druid"], path=src["path"]))
        pool = _records_from_bytes(raw, entry["id"])
        recs = []
        for r in entry["records"]:
            hits = [p for p in pool if p.header.startswith(r["header_prefix"])]
            if len(hits) != 1:
                raise ValueError(f"{entry['id']}: header {r['header_prefix']!r} matched {len(hits)} records")
            recs.append(hits[0])
    else:
        raise ValueError(f"{entry['id']}: unknown fetch kind {src['fetch']}")
    for r, spec in zip(recs, entry["records"]):
        if len(r.seq) != spec["length"]:
            raise ValueError(f"{entry['id']}: {r.header[:40]} length {len(r.seq)} != {spec['length']}")
    return recs


def _write_gz(path: Path, recs: list[fasta.Record]) -> bytes:
    text = b"\n".join(ln for r in recs for ln in fasta.format_record(r.header, r.seq, r.line_width or 70)) + b"\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as fh:
        with gzip.GzipFile(filename="", mode="wb", fileobj=fh, compresslevel=9, mtime=0) as gz:
            gz.write(text)
    return text


def fetch(golden: Path, out_dir: Path, lock_path: Path | None = None) -> dict:
    spec = yaml.safe_load(Path(golden).read_text())
    entries = spec["genomes"] + spec.get("spike_ins", [])
    lock = {"version": spec["version"], "retrieved": date.today().isoformat(), "files": {}}
    for entry in entries:
        entry.setdefault("kind", "spike_in")
        path = Path(out_dir) / entry["id"] / f"{entry['id']}.fna.gz"
        if path.exists():  # resume: re-validate what is already on disk
            recs = fasta.read(path)
            if [len(r.seq) for r in recs] != [r["length"] for r in entry["records"]]:
                raise ValueError(f"{path}: records on disk do not match golden.yaml; delete it and refetch")
            with gzip.open(path, "rb") as fh:
                text = fh.read()
        else:
            recs = _genome_records(entry, spec["fetch"])
            text = _write_gz(path, recs)
        lock["files"][entry["id"]] = {
            "path": f"{entry['id']}/{entry['id']}.fna.gz",
            "sha256_fna": hashlib.sha256(text).hexdigest(),
            "sha256_gz": hashlib.sha256(path.read_bytes()).hexdigest(),
            "records": [{"header": r.header, "length": len(r.seq)} for r in recs],
        }
        print(f"{entry['id']}: {len(recs)} record(s), {sum(len(r.seq) for r in recs):,} bp", file=sys.stderr)
    if lock_path is not None:
        lock_path = Path(lock_path)
        if lock_path.exists():
            old = yaml.safe_load(lock_path.read_text())
            for gid, f in lock["files"].items():
                if old["files"].get(gid, {}).get("sha256_fna") != f["sha256_fna"]:
                    raise SystemExit(f"{gid}: content differs from {lock_path}")
        else:
            lock_path.write_text(yaml.safe_dump(lock, sort_keys=False))
    return lock


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("golden", type=Path)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--lock", type=Path)
    a = ap.parse_args(argv)
    fetch(a.golden, a.out_dir, a.lock)


if __name__ == "__main__":
    main()
