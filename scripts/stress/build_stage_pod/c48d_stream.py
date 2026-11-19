#!/usr/bin/env python3
"""Card C-48d stage 2: bounded-memory streaming over JGI GEMs, IMG-UMAG (prok & euk),
and IUIG tar/zip archives without unpacking multi-hundred-GB tars onto `/tmp`.

Computes `contig_uid` digests (biobox.ids.contig_uid), exact overlap against stored
stages (`r10/root`) and staged inputs (`V1`, `PL2`, `G1`), GFF gene/intron diagnostics
for eukaryote bins, and optional bounded extraction of candidate `.fna.gz` files.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import resource
import sys
import tarfile
import time
import zipfile
from collections import Counter
from pathlib import Path

from biobox import ids


def _compute_n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    s = sorted(lengths, reverse=True)
    half = sum(s) / 2.0
    acc = 0
    for x in s:
        acc += x
        if acc >= half:
            return x
    return s[-1]


def _iter_fasta_bytes(fh) -> tuple[list[bytes], list[int], int]:
    """Stream a FASTA byte stream and return (contig_uids, lengths, gc_count)."""
    uids: list[bytes] = []
    lengths: list[int] = []
    gc_total = 0
    chunks: list[bytes] = []

    def _flush() -> None:
        nonlocal gc_total
        if not chunks:
            return
        seq = b"".join(chunks)
        chunks.clear()
        if not seq:
            return
        uids.append(ids.contig_uid(seq))
        lengths.append(len(seq))
        up = seq.upper()
        gc_total += up.count(b"G") + up.count(b"C")

    for line in fh:
        if line.startswith(b">"):
            _flush()
        else:
            s = line.strip()
            if s:
                chunks.append(s)
    _flush()
    return uids, lengths, gc_total


def load_known_contig_uids(root: Path | None = None, staged_inputs: list[Path] | None = None) -> dict[bytes, str]:
    """Load default-variant `contig_uid -> stage` from `r10/root` and staged input dirs."""
    owners: dict[bytes, str] = {}
    if root is not None and root.is_dir():
        import pyarrow.parquet as pq
        from biobox import io as bbio
        from biobox.l0seq import store

        for st in store.stages_on_disk(root):
            for pf in bbio.part_files(root / store.INDEX_DIR / f"stage={st}"):
                t = pq.read_table(pf, columns=["contig_uid", "is_default_variant"])
                for cu, is_def in zip(t.column("contig_uid").to_pylist(), t.column("is_default_variant").to_pylist()):
                    if is_def and cu not in owners:
                        owners[cu] = st
    for inp_dir in staged_inputs or []:
        gjson = inp_dir / "genomes.json"
        if not gjson.is_file():
            continue
        stage_label = f"staged:{inp_dir.name}"
        data = json.loads(gjson.read_text())
        for g in data:
            for hx in g.get("contig_uids", ()):
                cu = bytes.fromhex(hx)
                owners.setdefault(cu, stage_label)
    return owners


def stream_gems_fna_tar(
    tar_path: Path,
    candidate_ids: set[str] | None,
    known_contig_uids: dict[bytes, str],
    extract_dir: Path | None = None,
    max_members: int | None = None,
) -> dict:
    """Stream `jgi_gems_52515_fna_20261118.tar` (`mode='r|'`) in bounded memory."""
    if extract_dir is not None:
        extract_dir.mkdir(parents=True, exist_ok=True)
    members_seen = 0
    members_skipped = 0
    members_processed = 0
    shared_contigs = 0
    new_contigs = 0
    total_bp = 0
    shared_by_stage: Counter[str] = Counter()
    records: list[dict] = []

    with tarfile.open(tar_path, mode="r|") as tf:
        for member in tf:
            if not member.isfile() or not member.name.endswith(".fna.gz"):
                continue
            members_seen += 1
            if max_members is not None and members_seen > max_members:
                break
            gid = Path(member.name).name[:-7]  # strip .fna.gz
            if candidate_ids is not None and gid not in candidate_ids:
                members_skipped += 1
                continue

            fobj = tf.extractfile(member)
            if fobj is None:
                continue
            raw_gz = fobj.read()
            gz_sha256 = hashlib.sha256(raw_gz).hexdigest()
            with gzip.GzipFile(fileobj=io.BytesIO(raw_gz), mode="rb") as gzf:
                uids, lengths, gc = _iter_fasta_bytes(gzf)
            bp = sum(lengths)
            total_bp += bp
            n_shared = 0
            for cu in uids:
                st = known_contig_uids.get(cu)
                if st is not None:
                    n_shared += 1
                    shared_by_stage[st] += 1
            shared_contigs += n_shared
            new_contigs += len(uids) - n_shared
            members_processed += 1

            out_path_s = ""
            if extract_dir is not None:
                dst = extract_dir / f"{gid}.fna.gz"
                tmp = dst.with_name(f".{dst.name}.tmp")
                tmp.write_bytes(raw_gz)
                tmp.replace(dst)
                out_path_s = str(dst)

            records.append({
                "genome_id": gid,
                "bp": bp,
                "contigs": len(uids),
                "n50": _compute_n50(lengths),
                "gc_pct": round(100.0 * gc / bp, 2) if bp else 0.0,
                "shared_contigs": n_shared,
                "gz_bytes": len(raw_gz),
                "gz_sha256": gz_sha256,
                "extracted_path": out_path_s,
            })

    records.sort(key=lambda r: r["genome_id"])
    return {
        "members_seen": members_seen,
        "members_skipped": members_skipped,
        "members_processed": members_processed,
        "total_bp": total_bp,
        "shared_contigs": shared_contigs,
        "new_contigs": new_contigs,
        "shared_by_stage": dict(sorted(shared_by_stage.items())),
        "peak_rss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, 2),
        "records": records,
    }


def _parse_mbin_datafile(raw: bytes, target_bin_id: str) -> dict | None:
    """Parse JGI `mbin_datafile_<taxon_oid>.txt` for the given bin_id."""
    text = raw.decode("utf-8", errors="replace")
    header: list[str] | None = None
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if header is None:
            header = cols
            continue
        if cols and cols[0].strip() == target_bin_id:
            row = dict(zip(header, cols))
            comp = float(row.get("Bin Completeness") or 0.0)
            # JGI header spells "Bin Contamintation" (typo in JGI export); handle both spellings
            cont_s = row.get("Bin Contamintation") or row.get("Bin Contamination") or "0.0"
            cont = float(cont_s)
            return {
                "bin_quality": (row.get("Bin Quality") or "").strip(),
                "bin_lineage": (row.get("Bin Lineage") or "").strip(),
                "gtdb_lineage": (row.get("GTDB-TK lineage") or "").strip(),
                "completeness": comp,
                "contamination": cont,
                "reported_bp": int(float(row.get("Total Number of Bases") or 0)),
                "reported_genes": int(float(row.get("Number of genes") or 0)),
            }
    return None


def _parse_gff_stream(fh) -> dict:
    """Stream an IMG GFF byte stream in bounded memory (never loads full GFF into RAM)."""
    cds_count = 0
    total_cds_bp = 0
    # Sample up to 50,000 CDS lengths for exact/approximate median in bounded memory
    cds_sample: list[int] = []
    tools: Counter[str] = Counter()
    intron_or_exon = 0
    for raw_line in fh:
        if not raw_line or raw_line.startswith(b"#"):
            continue
        parts = raw_line.decode("utf-8", errors="replace").split("\t")
        if len(parts) < 5:
            continue
        tool, ftype = parts[1].strip(), parts[2].strip()
        if ftype == "CDS":
            tools[tool] += 1
            cds_count += 1
            try:
                clen = max(0, int(parts[4]) - int(parts[3]) + 1)
                total_cds_bp += clen
                if len(cds_sample) < 50_000:
                    cds_sample.append(clen)
            except ValueError:
                pass
        elif ftype in ("intron", "exon", "mRNA"):
            intron_or_exon += 1
    cds_sample.sort()
    med = cds_sample[len(cds_sample) // 2] if cds_sample else 0
    return {
        "gff_cds_count": cds_count,
        "gff_median_cds_bp": med,
        "gff_total_cds_bp": total_cds_bp,
        "gff_intron_or_exon_features": intron_or_exon,
        "gff_tools": sorted(tools.keys()),
    }


def _stream_fasta_and_maybe_extract(fh, dst_gz: Path | None) -> tuple[list[bytes], list[int], int, str]:
    """Stream FASTA bytes from `fh`, computing contig_uids/lengths/GC and optionally writing gzipped copy."""
    uids: list[bytes] = []
    lengths: list[int] = []
    gc_total = 0
    chunks: list[bytes] = []
    tmp_path: Path | None = None
    gz_out = None
    hasher = hashlib.sha256() if dst_gz is not None else None

    if dst_gz is not None:
        tmp_path = dst_gz.with_name(f".{dst_gz.name}.tmp")
        gz_out = gzip.GzipFile(filename="", mode="wb", fileobj=open(tmp_path, "wb"), compresslevel=6, mtime=0)

    def _flush() -> None:
        nonlocal gc_total
        if not chunks:
            return
        seq = b"".join(chunks)
        chunks.clear()
        if not seq:
            return
        uids.append(ids.contig_uid(seq))
        lengths.append(len(seq))
        up = seq.upper()
        gc_total += up.count(b"G") + up.count(b"C")

    try:
        for line in fh:
            if gz_out is not None:
                gz_out.write(line)
            if line.startswith(b">"):
                _flush()
            else:
                s = line.strip()
                if s:
                    chunks.append(s)
        _flush()
    finally:
        if gz_out is not None:
            raw_fh = gz_out.fileobj
            gz_out.close()
            raw_fh.close()

    gz_sha = ""
    if dst_gz is not None and tmp_path is not None and hasher is not None:
        with open(tmp_path, "rb") as rf:
            for blk in iter(lambda: rf.read(1 << 20), b""):
                hasher.update(blk)
        gz_sha = hasher.hexdigest()
        tmp_path.replace(dst_gz)
    return uids, lengths, gc_total, gz_sha


def stream_nested_bin_zip(
    zip_path: Path,
    inner_tar_name: str,
    mode: str,
    candidate_ids: set[str] | None,
    known_contig_uids: dict[bytes, str],
    extract_dir: Path | None = None,
    max_members: int | None = None,
    filter_hq_mq: bool = True,
    max_extract: int | None = None,
) -> dict:
    """Stream `.zip -> inner .tar -> per-bin/genome .tar.gz` in O(1) bounded memory.

    mode: "umag_euk", "umag_prok", or "iuig".
    """
    if mode not in ("umag_euk", "umag_prok", "iuig"):
        raise ValueError(f"unknown mode {mode!r}")
    if extract_dir is not None:
        extract_dir.mkdir(parents=True, exist_ok=True)

    members_seen = 0
    members_skipped = 0
    members_processed = 0
    hq_mq_kept = 0
    lq_filtered = 0
    extracted_count = 0
    shared_contigs = 0
    new_contigs = 0
    total_bp = 0
    shared_by_stage: Counter[str] = Counter()
    records: list[dict] = []

    with zipfile.ZipFile(zip_path, "r") as zf:
        with zf.open(inner_tar_name, "r") as inner_fh:
            with tarfile.open(fileobj=inner_fh, mode="r|") as outer_tf:
                for member in outer_tf:
                    if not member.isfile() or not member.name.endswith(".tar.gz"):
                        continue
                    members_seen += 1
                    if max_members is not None and members_seen > max_members:
                        break
                    base = Path(member.name).name
                    bin_id = base[:-7]  # strip .tar.gz
                    if candidate_ids is not None and bin_id not in candidate_ids:
                        members_skipped += 1
                        continue

                    fobj = outer_tf.extractfile(member)
                    if fobj is None:
                        continue
                    tgz_size = member.size
                    uids: list[bytes] = []
                    lengths: list[int] = []
                    gc = 0
                    gff_stats: dict | None = None
                    mbin_info: dict | None = None
                    want_extract = extract_dir is not None and (max_extract is None or extracted_count < max_extract)
                    dst_gz = (extract_dir / f"{bin_id}.fna.gz") if want_extract else None
                    gz_sha256 = ""
                    saw_fna = False

                    with tarfile.open(fileobj=fobj, mode="r|gz") as inner_tf:
                        for sub in inner_tf:
                            if not sub.isfile():
                                continue
                            sname = Path(sub.name).name
                            if sname.endswith(".fna") and not sname.endswith(".genes.fna") and not sname.endswith(".intergenic.fna"):
                                sf = inner_tf.extractfile(sub)
                                if sf is not None:
                                    uids, lengths, gc, gz_sha256 = _stream_fasta_and_maybe_extract(sf, dst_gz)
                                    saw_fna = True
                            elif mode == "umag_euk" and sname.endswith(".gff"):
                                sf = inner_tf.extractfile(sub)
                                if sf is not None:
                                    gff_stats = _parse_gff_stream(sf)
                            elif mode == "umag_prok" and sname.startswith("mbin_datafile_") and sname.endswith(".txt"):
                                sf = inner_tf.extractfile(sub)
                                if sf is not None:
                                    mbin_info = _parse_mbin_datafile(sf.read(), bin_id)

                    if mode == "umag_prok" and filter_hq_mq and mbin_info is not None:
                        if mbin_info["completeness"] < 50.0 or mbin_info["contamination"] > 10.0:
                            lq_filtered += 1
                            if dst_gz is not None and dst_gz.exists():
                                dst_gz.unlink(missing_ok=True)
                            continue
                        hq_mq_kept += 1
                    elif mode == "umag_prok":
                        hq_mq_kept += 1

                    if not saw_fna or not lengths:
                        continue

                    if dst_gz is not None and dst_gz.exists():
                        extracted_count += 1
                    out_path_s = str(dst_gz) if (dst_gz is not None and dst_gz.exists()) else ""

                    bp = sum(lengths)
                    total_bp += bp
                    n_shared = 0
                    for cu in uids:
                        st = known_contig_uids.get(cu)
                        if st is not None:
                            n_shared += 1
                            shared_by_stage[st] += 1
                    shared_contigs += n_shared
                    new_contigs += len(uids) - n_shared
                    members_processed += 1

                    rec: dict = {
                        "bin_id": bin_id,
                        "metagenome_id": bin_id.split("_")[0] if "_" in bin_id else bin_id,
                        "bp": bp,
                        "contigs": len(uids),
                        "n50": _compute_n50(lengths),
                        "gc_pct": round(100.0 * gc / bp, 2) if bp else 0.0,
                        "shared_contigs": n_shared,
                        "tgz_bytes": tgz_size,
                        "gz_sha256": gz_sha256,
                        "extracted_path": out_path_s,
                    }
                    if gff_stats:
                        rec.update(gff_stats)
                        rec["coding_density_pct"] = round(100.0 * gff_stats["gff_total_cds_bp"] / bp, 2) if bp else 0.0
                    if mbin_info:
                        rec.update(mbin_info)
                    records.append(rec)

    records.sort(key=lambda r: r["bin_id"])
    return {
        "mode": mode,
        "members_seen": members_seen,
        "members_skipped": members_skipped,
        "members_processed": members_processed,
        "hq_mq_kept": hq_mq_kept,
        "lq_filtered": lq_filtered,
        "total_bp": total_bp,
        "shared_contigs": shared_contigs,
        "new_contigs": new_contigs,
        "shared_by_stage": dict(sorted(shared_by_stage.items())),
        "peak_rss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, 2),
        "records": records,
    }


def _load_id_list(path: Path | None, col: str | None = None) -> set[str] | None:
    if path is None:
        return None
    out: set[str] = set()
    with open(path, encoding="utf-8") as fh:
        if col:
            for row in csv.DictReader(fh, delimiter="\t"):
                v = (row.get(col) or "").strip()
                if v:
                    out.add(v)
        else:
            for line in fh:
                v = line.strip()
                if v:
                    out.add(v)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Card C-48d bounded-memory archive streamer")
    ap.add_argument("--mode", choices=["gems", "umag_euk", "umag_prok", "iuig"], required=True)
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--inner-tar", default="")
    ap.add_argument("--filter-tsv", type=Path, default=None)
    ap.add_argument("--filter-col", default="genome_id")
    ap.add_argument("--root", type=Path, default=None)
    ap.add_argument("--staged-inputs", type=Path, nargs="*", default=[])
    ap.add_argument("--extract-dir", type=Path, default=None)
    ap.add_argument("--max-members", type=int, default=None)
    ap.add_argument("--max-extract", type=int, default=None)
    ap.add_argument("--out-json", type=Path, required=True)
    ap.add_argument("--out-tsv", type=Path, required=True)
    args = ap.parse_args(argv)

    t0 = time.monotonic()
    known = load_known_contig_uids(args.root, args.staged_inputs)
    cand_ids = _load_id_list(args.filter_tsv, args.filter_col)

    if args.mode == "gems":
        res = stream_gems_fna_tar(
            args.archive, cand_ids, known, extract_dir=args.extract_dir, max_members=args.max_members
        )
    else:
        res = stream_nested_bin_zip(
            args.archive,
            args.inner_tar,
            args.mode,
            cand_ids,
            known,
            extract_dir=args.extract_dir,
            max_members=args.max_members,
            max_extract=args.max_extract,
        )

    res["wall_sec"] = round(time.monotonic() - t0, 2)
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    recs = res.pop("records")
    if recs:
        cols = sorted({k for r in recs for k in r.keys()})
        with open(args.out_tsv, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in recs:
                r_out = {k: (",".join(v) if isinstance(v, list) else v) for k, v in r.items()}
                w.writerow(r_out)
    else:
        args.out_tsv.write_text("")

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(res, indent=2, sort_keys=True) + "\n")
    print(json.dumps(res, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
