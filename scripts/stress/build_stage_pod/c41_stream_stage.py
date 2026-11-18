"""Card C-41: stream GlobDB r226's genome FASTA archive over one HTTP connection and extract one proposed stage.

    python c41_stream_stage.py --url URL --expected-bytes 251371669660 \
        --expected-md5 a34a10d72cd688ed518eb6c0e2c30d26 \
        --stages-tsv /tmp/if2453/s3/globdb/tables/non_gtdb_stages.tsv \
        --stage G1 --out-dir /tmp/if2453/s3/globdb/genomes/G1

Never stores the 251.4 GB archive on /tmp: streams the response body through a counting MD5+SHA256 reader into
`tarfile.open(..., mode="r|gz")`, writes only the target stage's `<ID>.fa.gz` files (~3 GB for ~10 Gbp), checks every
extracted genome's contig count and base count against `globdb_r226_tax_plus_stats.tsv.gz`, checks `/tmp` free space
(>= `--min-free-gb`, default 450 GiB) and `STOP` files periodically, and writes `lock.json`, `selection.tsv`, and
`stream_summary.json`.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io as pyio
import json
import resource
import shutil
import sys
import tarfile
import time
import traceback
import urllib.error
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from biobox import fasta


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def cpu_sec() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def mag_quality(comp: float, cont: float) -> str:
    """MIMAG completeness/contamination tier mapped to `genome_input.assembly_quality`."""
    if comp > 90.0 and cont < 5.0:
        return "mag_hq"
    if comp >= 50.0 and cont < 10.0:
        return "mag_mq"
    return "contig"


class HashingStream(pyio.RawIOBase):
    """Wrap an HTTP response stream: count bytes, compute MD5 + SHA256, check STOP and /tmp free space."""

    def __init__(self, resp, stop_files: list[Path], tmp_mount: Path, min_free_bytes: int, log_fh):
        self._resp = resp
        self._stop_files = stop_files
        self._tmp_mount = tmp_mount
        self._min_free = min_free_bytes
        self._log = log_fh
        self.bytes_read = 0
        self.md5 = hashlib.md5()
        self.sha256 = hashlib.sha256()
        self._next_check = 0
        self.t0 = time.monotonic()

    def readable(self) -> bool:
        return True

    def readinto(self, b) -> int:
        data = self._resp.read(len(b))
        if not data:
            return 0
        n = len(data)
        b[:n] = data
        self.bytes_read += n
        self.md5.update(data)
        self.sha256.update(data)
        if self.bytes_read >= self._next_check:
            for sf in self._stop_files:
                if sf.exists():
                    raise RuntimeError(f"STOP file present: {sf}")
            free = shutil.disk_usage(self._tmp_mount).free
            if free < self._min_free:
                raise RuntimeError(f"/tmp free {free / (1 << 30):.1f} GiB below floor {self._min_free / (1 << 30):.1f} GiB")
            self._next_check = self.bytes_read + (1 << 30)  # check every 1 GiB streamed
        return n


def load_stage_targets(stages_tsv: Path, stage: str) -> list[dict]:
    rows = []
    with open(stages_tsv, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["stage"] == stage:
                rows.append({
                    "ID": r["ID"],
                    "position": int(r["position_in_stage"]),
                    "resource": r["resource"],
                    "genome_size": int(r["genome_size"]),
                    "contig_number": int(r["contig_number"]),
                    "checkm2_completeness": float(r["checkm2_completeness"]),
                    "checkm2_contamination": float(r["checkm2_contamination"]),
                    "taxonomy": r["taxonomy"],
                })
    rows.sort(key=lambda r: r["position"])
    if not rows:
        raise SystemExit(f"no rows for stage {stage!r} in {stages_tsv}")
    return rows


def verify_and_write_fa_gz(gz_bytes: bytes, dst: Path, want_records: int, want_bp: int, gid: str) -> dict:
    gz_sha256 = hashlib.sha256(gz_bytes).hexdigest()
    gz_md5 = hashlib.md5(gz_bytes).hexdigest()
    raw_fa = gzip.decompress(gz_bytes)
    fa_sha256 = hashlib.sha256(raw_fa).hexdigest()
    fa_md5 = hashlib.md5(raw_fa).hexdigest()
    # Keep `.fa.gz` suffix on temporary file so `biobox.fasta.iter_records` opens it with `gzip.open`.
    tmp = dst.with_name(f".tmp.{dst.name}")
    tmp.write_bytes(gz_bytes)
    rec_count = 0
    bp = 0
    has_n = 0
    try:
        for rec in fasta.iter_records(tmp):
            rec_count += 1
            bp += len(rec.seq)
            if b"N" in rec.seq or b"n" in rec.seq:
                has_n += 1
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    if rec_count != want_records or bp != want_bp:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(
            f"{gid}: FASTA has {rec_count} records / {bp} bp, table says {want_records} / {want_bp}"
        )
    tmp.replace(dst)
    return {
        "gz_bytes": len(gz_bytes),
        "gz_sha256": gz_sha256,
        "gz_md5": gz_md5,
        "bytes": len(raw_fa),
        "sha256": fa_sha256,
        "md5": fa_md5,
        "records": rec_count,
        "bp": bp,
        "records_with_N": has_n,
    }


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", required=True)
    ap.add_argument("--expected-bytes", type=int, required=True)
    ap.add_argument("--expected-md5", required=True)
    ap.add_argument("--stages-tsv", type=Path, required=True)
    ap.add_argument("--stage", required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--min-free-gb", type=float, default=450.0)
    ap.add_argument("--stop-when-found", action="store_true",
                    help="Close HTTP stream once every genome of --stage is extracted and verified.")
    a = ap.parse_args(argv)

    targets = load_stage_targets(a.stages_tsv, a.stage)
    by_id = {r["ID"]: r for r in targets}
    files_dir = a.out_dir / "files"
    files_dir.mkdir(parents=True, exist_ok=True)
    for stale in files_dir.iterdir():
        if stale.name.startswith("."):
            stale.unlink(missing_ok=True)
    log_path = a.out_dir / "stream.log"

    stop_files = [Path("/tmp/if2453/r10/root/STOP"), Path("/tmp/if2453/s3/globdb/STOP")]
    min_free_bytes = int(a.min_free_gb * (1 << 30))

    with open(log_path, "a", buffering=1) as log_fh:
        def log(msg: str) -> None:
            line = f"[{utc_now()}] {msg}"
            print(line, flush=True)
            log_fh.write(line + "\n")

        try:
            return _run_stream(a, targets, by_id, files_dir, stop_files, min_free_bytes, log_fh, log)
        except BaseException:
            tb = traceback.format_exc()
            log(f"EXCEPTION:\n{tb}")
            return 1


def _run_stream(a, targets, by_id, files_dir: Path, stop_files: list[Path], min_free_bytes: int, log_fh, log) -> int:
    if True:

        free0 = shutil.disk_usage("/tmp").free
        log(f"start stage={a.stage} targets={len(by_id)} url={a.url} tmp_free_gib={free0 / (1 << 30):.2f}")
        req = urllib.request.Request(a.url, headers={"User-Agent": "curl"})
        try:
            resp = urllib.request.urlopen(req, timeout=300)
        except urllib.error.HTTPError as e:
            log(f"HTTP error {e.code}: {e.reason}")
            return 1
        cl = resp.headers.get("Content-Length")
        if cl is not None and int(cl) != a.expected_bytes:
            log(f"Content-Length {cl} != expected {a.expected_bytes}")
            return 1

        hstream = HashingStream(resp, stop_files, Path("/tmp"), min_free_bytes, log_fh)
        buf_stream = pyio.BufferedReader(hstream, buffer_size=1 << 20)

        lock_genomes: dict[str, dict] = {}
        chunks_seen: list[str] = []
        chunks_set: set[str] = set()
        members_files = 0
        last_log = time.monotonic()
        stopped_early = False

        with tarfile.open(fileobj=buf_stream, mode="r|gz") as tf:
            for m in tf:
                if not m.isfile():
                    continue
                members_files += 1
                parts = m.name.split("/")
                if len(parts) >= 2 and parts[1] not in chunks_set:
                    chunks_set.add(parts[1])
                    chunks_seen.append(parts[1])
                fname = parts[-1]
                if fname.endswith(".fa.gz"):
                    gid = fname[:-6]
                    target = by_id.get(gid)
                    if target is not None and gid not in lock_genomes:
                        fh = tf.extractfile(m)
                        assert fh is not None
                        gz_bytes = fh.read()
                        dst = files_dir / fname
                        stats = verify_and_write_fa_gz(
                            gz_bytes, dst, target["contig_number"], target["genome_size"], gid
                        )
                        lock_genomes[gid] = {
                            "file": str(dst),
                            "member": m.name,
                            **stats,
                            "retrieved_utc": utc_now(),
                        }
                now_m = time.monotonic()
                if now_m - last_log >= 30.0:
                    dt = now_m - hstream.t0
                    mbps = (hstream.bytes_read / (1 << 20)) / dt if dt > 0 else 0.0
                    free_g = shutil.disk_usage("/tmp").free / (1 << 30)
                    log(
                        f"progress streamed={hstream.bytes_read / (1 << 30):.2f}GiB "
                        f"({100.0 * hstream.bytes_read / a.expected_bytes:.1f}%) "
                        f"rate={mbps:.1f}MiB/s files_seen={members_files} "
                        f"stage_found={len(lock_genomes)}/{len(by_id)} "
                        f"cpu_min={cpu_sec() / 60.0:.2f} tmp_free={free_g:.1f}GiB"
                    )
                    last_log = now_m
                if a.stop_when_found and len(lock_genomes) == len(by_id):
                    stopped_early = True
                    log(f"all {len(by_id)} target genomes extracted and verified; stopping stream early as requested")
                    break

            if not stopped_early:
                # Drain any trailing gzip/tar padding bytes to reach exact EOF for full archive MD5
                while buf_stream.read(1 << 20):
                    pass

        resp.close()
        dt_total = time.monotonic() - hstream.t0
        got_md5 = hstream.md5.hexdigest()
        got_sha256 = hstream.sha256.hexdigest()

        missing = sorted(set(by_id) - set(lock_genomes))
        if missing:
            log(f"FAILED: missing {len(missing)} stage genomes (first: {missing[0]})")
            return 1

        if not stopped_early:
            if hstream.bytes_read != a.expected_bytes:
                log(f"FAILED: stream bytes {hstream.bytes_read} != expected {a.expected_bytes}")
                return 1
            if got_md5 != a.expected_md5:
                log(f"FAILED: stream md5 {got_md5} != expected {a.expected_md5}")
                return 1

        # Write selection.tsv in stage order
        qual_counts = Counter()
        sel_path = a.out_dir / "selection.tsv"
        with open(sel_path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow([
                "accession", "translation_table", "assembly_quality", "resource",
                "genome_size", "contig_number", "checkm2_completeness", "checkm2_contamination", "taxonomy",
            ])
            for r in targets:
                q = mag_quality(r["checkm2_completeness"], r["checkm2_contamination"])
                qual_counts[q] += 1
                w.writerow([
                    r["ID"], "", q, r["resource"], r["genome_size"], r["contig_number"],
                    r["checkm2_completeness"], r["checkm2_contamination"], r["taxonomy"],
                ])

        lock_obj = {
            "source": "globdb_r226",
            "stage": a.stage,
            "archive_url": a.url,
            "archive_complete": not stopped_early,
            "archive_stream_bytes": hstream.bytes_read,
            "archive_stream_md5": got_md5,
            "archive_stream_sha256": got_sha256,
            "archive_expected_bytes": a.expected_bytes,
            "archive_expected_md5": a.expected_md5,
            "genomes": {gid: lock_genomes[gid] for gid in sorted(lock_genomes)},
        }
        (a.out_dir / "lock.json").write_text(json.dumps(lock_obj, indent=1) + "\n")

        summary = {
            "stage": a.stage,
            "genomes_extracted": len(lock_genomes),
            "total_gz_bytes": sum(v["gz_bytes"] for v in lock_genomes.values()),
            "total_fasta_bytes": sum(v["bytes"] for v in lock_genomes.values()),
            "total_records": sum(v["records"] for v in lock_genomes.values()),
            "total_bp": sum(v["bp"] for v in lock_genomes.values()),
            "genomes_with_N_bases": sum(1 for v in lock_genomes.values() if v["records_with_N"] > 0),
            "records_with_N_bases": sum(v["records_with_N"] for v in lock_genomes.values()),
            "assembly_quality_counts": dict(sorted(qual_counts.items())),
            "archive_complete": not stopped_early,
            "archive_stream_bytes": hstream.bytes_read,
            "archive_stream_md5": got_md5,
            "archive_stream_sha256": got_sha256,
            "archive_md5_matches_provider": (not stopped_early and got_md5 == a.expected_md5),
            "archive_files_seen": members_files,
            "archive_chunks_seen": len(chunks_seen),
            "wall_seconds": round(dt_total, 2),
            "cpu_minutes": round(cpu_sec() / 60.0, 3),
            "tmp_free_bytes_after": shutil.disk_usage("/tmp").free,
        }
        (a.out_dir / "stream_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
        log(f"done summary={json.dumps(summary)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
