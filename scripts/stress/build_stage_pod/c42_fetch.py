"""Card C-42 polite single-connection fetcher for SPIRE tables and per-sample assemblies.

    python -I c42_fetch.py tables --out-dir /tmp/if2453/s3/spire/tables [--interval 1.0]
    python -I c42_fetch.py assemblies --selection /tmp/if2453/s3/spire/select/pilot_selection.tsv
        --out-dir /tmp/if2453/s3/spire/dl/assemblies --spire-dir /tmp/if2453/s3/spire
        [--interval 0.5] [--max-scratch-gb 165]
    python -I c42_fetch.py mag-probe --sample SAM --genome-meta GZ --assembly-fa FA_GZ --out-json J

Rules (AGENTS §8, DEC-50):
- Every URL fetched is listed in `sources.yaml` (`spire_tables`, `spire_study_api`, `spire_pilot_assemblies`)
  before any byte is requested.
- Sends `User-Agent: curl` (no name, email address or project detail).
- Single connection; at least `--interval` seconds between request starts (never > 2 req/s).
- Stops immediately on 5 HTTP refusals in a row (403, 429, 500, 502, 503, 504), or if `STOP` appears
  in `/tmp/if2453/r10/root` or `/tmp/if2453/s3/spire`, or if `/tmp/if2453/s3/spire` exceeds `--max-scratch-gb`.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import re
import socket
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zlib
from pathlib import Path

from biobox import fasta, ids

REFUSALS = {403, 429, 500, 502, 503, 504}
MAX_REFUSALS = 5
UA = {"User-Agent": "curl"}

SPIRE_TABLES = [
    ("README", "https://black.embl.de/~fullam/spire/metadata/README", "51152f40061eeab04fabbe18d37c22ac"),
    ("spire_v1_cluster_metadata.tsv.gz", "https://black.embl.de/~fullam/spire/metadata/spire_v1_cluster_metadata.tsv.gz", "3f572296ed0d5918b18916a69b39fd"),
    ("spire_v1_genome_metadata.tsv.gz", "https://black.embl.de/~fullam/spire/metadata/spire_v1_genome_metadata.tsv.gz", "a4ba393df92e9d1a35b8b996a69b39fd"),
    ("ena_drop.v_2_0.microntology.v_1_0.annotations.v_0_2.tsv.gz", "https://black.embl.de/~fullam/spire/microntology/ena_drop.v_2_0.microntology.v_1_0.annotations.v_0_2.tsv.gz", "62056caa199fa3aca9b019f400149527"),
]
# Exact MD5 from sources.yaml for spire_v1_cluster_metadata.tsv.gz:
SPIRE_TABLES[1] = (
    "spire_v1_cluster_metadata.tsv.gz",
    "https://black.embl.de/~fullam/spire/metadata/spire_v1_cluster_metadata.tsv.gz",
    "3f572296ed0d5918b18916b6a9099618",
)


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _check_stop(spire_dir: Path) -> None:
    for p in (Path("/tmp/if2453/r10/root/STOP"), spire_dir / "STOP"):
        if p.exists():
            print(f"{_now()} STOP present at {p}; exiting cleanly", file=sys.stderr, flush=True)
            raise SystemExit(75)


def _dir_bytes(path: Path) -> int:
    total = 0
    if not path.exists():
        return 0
    for root, _, files in os.walk(path):
        for fn in files:
            try:
                total += os.path.getsize(os.path.join(root, fn))
            except OSError:
                pass
    return total


class PoliteClient:
    def __init__(self, interval: float):
        if interval < 0.5:
            raise ValueError(f"interval must be >= 0.5 s (<= 2 req/s), got {interval}")
        self.interval = interval
        self.last = 0.0
        self.refusals_in_a_row = 0

    def fetch_to_file(self, url: str, dst: Path, *, expected_md5: str | None = None, ensure_gzip: bool = False) -> dict:
        part = dst.with_name(f".{dst.name}.part")
        while True:
            wait = self.last + self.interval - time.monotonic()
            if wait > 0:
                time.sleep(wait)
            self.last = time.monotonic()
            req = urllib.request.Request(url, headers=UA)
            try:
                with urllib.request.urlopen(req, timeout=180) as resp:
                    final_url = resp.geturl()
                    status = resp.status
                    first = resp.read(1 << 20)
                    raw_md5 = hashlib.md5()
                    raw_sha = hashlib.sha256()
                    wire_bytes = 0

                    is_gz = first[:2] == b"\x1f\x8b"
                    if ensure_gzip and not is_gz:
                        # Plain FASTA over wire -> store deterministically gzipped on disk
                        disk_sha = hashlib.sha256()
                        disk_md5 = hashlib.md5()
                        uncomp_bytes = 0
                        with open(part, "wb") as raw_out:
                            with gzip.GzipFile(filename="", mode="wb", fileobj=raw_out, compresslevel=6, mtime=0) as gz:
                                chunk = first
                                while chunk:
                                    wire_bytes += len(chunk)
                                    uncomp_bytes += len(chunk)
                                    raw_md5.update(chunk)
                                    raw_sha.update(chunk)
                                    gz.write(chunk)
                                    chunk = resp.read(1 << 20)
                            raw_out.flush()
                            os.fsync(raw_out.fileno())
                        with open(part, "rb") as rf:
                            while b := rf.read(1 << 20):
                                disk_sha.update(b)
                                disk_md5.update(b)
                        disk_bytes = part.stat().st_size
                    else:
                        decomp = zlib.decompressobj(16 + zlib.MAX_WBITS) if is_gz else None
                        uncomp_bytes = 0
                        with open(part, "wb") as out_f:
                            chunk = first
                            while chunk:
                                wire_bytes += len(chunk)
                                raw_md5.update(chunk)
                                raw_sha.update(chunk)
                                if decomp is not None:
                                    uncomp_bytes += len(decomp.decompress(chunk))
                                else:
                                    uncomp_bytes += len(chunk)
                                out_f.write(chunk)
                                chunk = resp.read(1 << 20)
                            if decomp is not None:
                                uncomp_bytes += len(decomp.flush())
                            out_f.flush()
                            os.fsync(out_f.fileno())
                        disk_bytes = wire_bytes
                        disk_sha = raw_sha
                        disk_md5 = raw_md5

                    self.refusals_in_a_row = 0
                    break
            except urllib.error.HTTPError as e:
                part.unlink(missing_ok=True)
                if e.code == 404:
                    self.refusals_in_a_row = 0
                    return {"url": url, "status": 404, "not_found": True, "fetched_utc": _now()}
                if e.code not in REFUSALS:
                    raise
                why = f"HTTP {e.code}"
            except (TimeoutError, socket.timeout, urllib.error.URLError, zlib.error, EOFError) as e:
                part.unlink(missing_ok=True)
                why = f"{type(e).__name__}: {e}"

            self.refusals_in_a_row += 1
            print(f"{_now()} refusal {self.refusals_in_a_row}/{MAX_REFUSALS} ({why}) on {url}", file=sys.stderr, flush=True)
            if self.refusals_in_a_row >= MAX_REFUSALS:
                raise SystemExit(4)
            time.sleep(15 * self.refusals_in_a_row)

        if expected_md5 is not None and raw_md5.hexdigest() != expected_md5.lower():
            part.unlink(missing_ok=True)
            raise RuntimeError(f"MD5 mismatch on {url}: got {raw_md5.hexdigest()}, expected {expected_md5}")

        part.replace(dst)
        return {
            "url": url,
            "final_url": final_url,
            "path": str(dst),
            "status": status,
            "wire_gzipped": is_gz,
            "wire_bytes": wire_bytes,
            "disk_bytes": disk_bytes,
            "uncompressed_bytes": uncomp_bytes,
            "wire_md5": raw_md5.hexdigest(),
            "wire_sha256": raw_sha.hexdigest(),
            "disk_md5": disk_md5.hexdigest(),
            "disk_sha256": disk_sha.hexdigest(),
            "fetched_utc": _now(),
        }


def cmd_tables(args: argparse.Namespace) -> int:
    out = Path(args.out_dir)
    spire_dir = out.parent
    out.mkdir(parents=True, exist_ok=True)
    client = PoliteClient(max(0.5, args.interval))
    log_path = out / "fetch_tables_log.jsonl"

    def record(rec: dict) -> None:
        with open(log_path, "a") as fh:
            fh.write(json.dumps(rec) + "\n")

    # 1. spire_tables
    for fname, url, md5 in SPIRE_TABLES:
        _check_stop(spire_dir)
        dst = out / fname
        if not dst.is_file():
            rec = client.fetch_to_file(url, dst, expected_md5=md5)
            record({"group": "spire_tables", "file": fname, **rec})

    # 2. spire_study_api: studies.html + 715 study TSVs
    _check_stop(spire_dir)
    studies_html = out / "studies.html"
    if not studies_html.is_file():
        rec = client.fetch_to_file("https://spire.embl.de/studies", studies_html)
        record({"group": "spire_study_api", "file": "studies.html", **rec})

    html_txt = studies_html.read_text()
    s0 = html_txt.index("var studiesData")
    s1 = html_txt.index("];", s0)
    studies = sorted({
        s for s, _ in re.findall(r'\{\s*study:\s*"((?:[^"\\]|\\.)*)",\s*samples:\s*(\d+)\s*\}', html_txt[s0:s1])
    })

    study_dir = out / "study_tsv"
    study_dir.mkdir(parents=True, exist_ok=True)
    for study in studies:
        _check_stop(spire_dir)
        enc = urllib.parse.quote(study, safe="")
        dst = study_dir / f"{enc}.tsv"
        if dst.is_file():
            continue
        url = f"https://spire.embl.de/spire/api/study/{enc}?format=tsv"
        rec = client.fetch_to_file(url, dst)
        record({"group": "spire_study_api", "study": study, **rec})

    # 3. spire_pilot_assemblies: api/v1/samples?format=tsv&limit=10000&offset=...
    api_dir = out / "api_samples"
    api_dir.mkdir(parents=True, exist_ok=True)
    offset = 0
    page = 0
    while True:
        _check_stop(spire_dir)
        dst = api_dir / f"page_{page:04d}_offset_{offset}.tsv"
        if not dst.is_file():
            url = f"https://spire.embl.de/api/v1/samples?format=tsv&limit=10000&offset={offset}"
            rec = client.fetch_to_file(url, dst)
            record({"group": "spire_samples_api", "offset": offset, **rec})
        with open(dst, newline="") as fh:
            n_rows = sum(1 for _ in fh) - 1  # minus header
        if n_rows < 10000:
            break
        offset += 10000
        page += 1

    (out / "TABLES_DONE").write_text(_now() + "\n")
    print(json.dumps({"studies": len(studies), "api_sample_pages": page + 1, "done_utc": _now()}))
    return 0


def cmd_assemblies(args: argparse.Namespace) -> int:
    out = Path(args.out_dir)
    spire_dir = Path(args.spire_dir)
    out.mkdir(parents=True, exist_ok=True)
    client = PoliteClient(max(0.5, args.interval))
    log_path = out / "fetch_assemblies.jsonl"
    cap_bytes = int(args.max_scratch_gb * (1024 ** 3))

    with open(args.selection, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    done_sids: set[str] = set()
    if log_path.is_file():
        with open(log_path) as fh:
            for ln in fh:
                if ln.strip():
                    r = json.loads(ln)
                    if r.get("sample_id") and (r.get("not_found") or (out / f"{r['sample_id']}.fa.gz").is_file()):
                        done_sids.add(r["sample_id"])

    fetched = skipped = not_found = 0
    for idx, r in enumerate(rows, 1):
        sid = r["sample_id"]
        dst = out / f"{sid}.fa.gz"
        if sid in done_sids or dst.is_file():
            skipped += 1
            continue
        _check_stop(spire_dir)
        if idx % 20 == 1:
            used = _dir_bytes(spire_dir)
            if used >= cap_bytes:
                print(f"{_now()} scratch {used / 1e9:.2f} GB >= cap {cap_bytes / 1e9:.2f} GB; stopping", file=sys.stderr)
                return 3

        enc = urllib.parse.quote(sid, safe="")
        url = f"https://spire.embl.de/download_assembly/{enc}"
        rec = client.fetch_to_file(url, dst, ensure_gzip=True)
        rec["sample_id"] = sid
        rec["habitat"] = r["habitat"]
        rec["habitat2"] = r["habitat2"]
        if rec.get("not_found"):
            not_found += 1
        else:
            fetched += 1
        with open(log_path, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        if fetched % 50 == 0:
            print(f"{_now()} fetched={fetched} skipped={skipped} not_found={not_found} idx={idx}/{len(rows)}", flush=True)

    (out / "ASSEMBLIES_DONE").write_text(_now() + "\n")
    print(json.dumps({"fetched": fetched, "skipped": skipped, "not_found": not_found, "total": len(rows)}))
    return 0


def cmd_mag_probe(args: argparse.Namespace) -> int:
    """Verify contig_uid identity between one sample's assembly and its SPIRE v1 MAG bins."""
    sid = args.sample
    mags = []
    with gzip.open(args.genome_meta, "rt", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("derived_from_sample") or r.get("sample_id")) == sid:
                gid = r.get("genome_id") or r["genome_name"]
                nc_raw = r["n_contigs"] if "n_contigs" in r else r.get("num_contigs", 0)
                mags.append((gid, int(float(r["genome_size"])), int(float(nc_raw or 0))))

    asm_uids = {ids.contig_uid(rec.seq): len(rec.seq) for rec in fasta.iter_records(Path(args.assembly_fa))}
    client = PoliteClient(0.5)
    out_p = Path(args.out_json)
    tmp_dir = out_p.parent / "mag_probe_files"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    mag_results = []
    all_mag_uids: set[bytes] = set()
    overlap_between_mags = 0
    for mag_id, meta_bp, meta_n in mags:
        dst = tmp_dir / f"{mag_id}.fa.gz"
        if not dst.is_file():
            client.fetch_to_file(f"https://spire.embl.de/download_file/{urllib.parse.quote(mag_id, safe='')}", dst, ensure_gzip=True)
        m_uids = [ids.contig_uid(rec.seq) for rec in fasta.iter_records(dst)]
        m_lens = [len(rec.seq) for rec in fasta.iter_records(dst)]
        in_asm = sum(1 for u in m_uids if u in asm_uids)
        overlap_between_mags += sum(1 for u in m_uids if u in all_mag_uids)
        all_mag_uids.update(m_uids)
        mag_results.append({
            "mag_id": mag_id,
            "metadata_bp": meta_bp,
            "fasta_bp": sum(m_lens),
            "contigs": len(m_uids),
            "min_contig_len": min(m_lens) if m_lens else 0,
            "contigs_in_sample_assembly": in_asm,
        })

    res = {
        "sample_id": sid,
        "sample_assembly_contigs": len(asm_uids),
        "sample_assembly_bp": sum(asm_uids.values()),
        "mags_checked": len(mag_results),
        "total_mag_metadata_bp": sum(r["metadata_bp"] for r in mag_results),
        "total_mag_fasta_bp": sum(r["fasta_bp"] for r in mag_results),
        "all_mag_contigs_in_assembly": all(r["contigs_in_sample_assembly"] == r["contigs"] for r in mag_results),
        "contig_overlap_between_mags": overlap_between_mags,
        "mags": mag_results,
    }
    out_p.write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "mags"}))
    return 0


def cmd_sample_mags(args: argparse.Namespace) -> int:
    """Fetch MAG FASTAs for a small list of samples (e.g. TniQ-hit samples) and record their contig_uid hexes."""
    sids = {ln.strip() for ln in Path(args.samples_file).read_text().splitlines() if ln.strip()}
    mags_by_sample: dict[str, list[str]] = {}
    with gzip.open(args.genome_meta, "rt", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sid = (r.get("derived_from_sample") or r.get("sample_id") or "").strip()
            if sid in sids:
                gid = r.get("genome_id") or r["genome_name"]
                mags_by_sample.setdefault(sid, []).append(gid)

    out_json = Path(args.out_json)
    tmp_dir = out_json.parent / "sample_mags_fa"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    spire_dir = Path(args.spire_dir)
    client = PoliteClient(max(0.5, args.interval))
    log_path = out_json.parent / "fetch_sample_mags.jsonl"

    uid_to_mag: dict[str, str] = {}
    mags_fetched = 0
    for sid in sorted(mags_by_sample):
        for mag_id in sorted(mags_by_sample[sid]):
            _check_stop(spire_dir)
            dst = tmp_dir / f"{mag_id}.fa.gz"
            if not dst.is_file():
                rec = client.fetch_to_file(
                    f"https://spire.embl.de/download_file/{urllib.parse.quote(mag_id, safe='')}",
                    dst,
                    ensure_gzip=True,
                )
                rec["sample_id"] = sid
                rec["mag_id"] = mag_id
                with open(log_path, "a") as lfh:
                    lfh.write(json.dumps(rec) + "\n")
                if rec.get("not_found"):
                    continue
                mags_fetched += 1
            for rec_fa in fasta.iter_records(dst):
                uid_to_mag[ids.contig_uid(rec_fa.seq).hex()] = mag_id
            if args.unlink_after:
                dst.unlink(missing_ok=True)

    payload = {
        "samples_requested": len(sids),
        "samples_with_mags": len(mags_by_sample),
        "mags_total": sum(len(v) for v in mags_by_sample.values()),
        "mags_fetched_now": mags_fetched,
        "distinct_mag_contig_uids": len(uid_to_mag),
        "contig_uid_to_mag": uid_to_mag,
    }
    out_json.write_text(json.dumps(payload) + "\n")
    print(json.dumps({k: v for k, v in payload.items() if k != "contig_uid_to_mag"}))
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_tab = sub.add_parser("tables")
    p_tab.add_argument("--out-dir", required=True)
    p_tab.add_argument("--interval", type=float, default=1.0)

    p_asm = sub.add_parser("assemblies")
    p_asm.add_argument("--selection", required=True)
    p_asm.add_argument("--out-dir", required=True)
    p_asm.add_argument("--spire-dir", required=True)
    p_asm.add_argument("--interval", type=float, default=0.5)
    p_asm.add_argument("--max-scratch-gb", type=float, default=165.0)

    p_mp = sub.add_parser("mag-probe")
    p_mp.add_argument("--sample", required=True)
    p_mp.add_argument("--genome-meta", required=True)
    p_mp.add_argument("--assembly-fa", required=True)
    p_mp.add_argument("--out-json", required=True)

    p_sm = sub.add_parser("sample-mags")
    p_sm.add_argument("--samples-file", required=True)
    p_sm.add_argument("--genome-meta", required=True)
    p_sm.add_argument("--spire-dir", required=True)
    p_sm.add_argument("--out-json", required=True)
    p_sm.add_argument("--interval", type=float, default=0.5)
    p_sm.add_argument("--unlink-after", action="store_true")

    args = ap.parse_args(argv)
    if args.cmd == "tables":
        return cmd_tables(args)
    if args.cmd == "assemblies":
        return cmd_assemblies(args)
    if args.cmd == "mag-probe":
        return cmd_mag_probe(args)
    if args.cmd == "sample-mags":
        return cmd_sample_mags(args)
    return 2


if __name__ == "__main__":
    sys.exit(main())
