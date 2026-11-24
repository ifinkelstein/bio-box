"""Card C-67 (DEC-69, DEC-65, DEC-50, DEC-26, DEC-27):
Fetch, deduplicate, and stage:
  - `LRV1`: compact `redistributable: true` long-read marine & freshwater complete viromes
    (Beaulaurier et al. 2020 Station ALOHA ONT AFVGs PRJNA529454, Kavagutti et al. 2023
    freshwater HiFi complete phages Figshare 20431332, and VirION2 long-read marine viromes
    Zablocki et al. 2021 PeerJ e11088 / PMC8018248) with DEC-65 exact 20..1,000 bp DTR trimming
    into `/tmp/if2453/s2/inputs/LRV1/`.
  - `V4`: `PEONY 1.3` / `PIGEON 3.0` global soil/freshwater viromes (Dryad 10.5061/dryad.prr4xgz35,
    CC0-1.0, `redistributable: true`) streamed into `/stor/scratch/Wilke/if2453/bio-box-staging/V4/`
    with metadata in `/tmp/if2453/c67/V4/` (Part 1 split layout).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import Counter
from collections.abc import Iterable
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # c48c_virus_inputs and make_inputs beside this script

import c48c_virus_inputs as v2  # noqa: E402
import make_inputs  # noqa: E402
from biobox import fasta, ids  # noqa: E402

csv.field_size_limit(1 << 30)
Refused = make_inputs.Refused

VALID_DNA_RE = re.compile(rb"^[ACGTURYKMSWBDHVNacgturykmswbdhvn]+$")

KAVAGUTTI_URL = "https://ndownloader.figshare.com/files/36556458"
KAVAGUTTI_MD5 = "cfbc18185b27bf2d6c2a9c9d1765def7"
VIRION2_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8018248/supplementaryFiles"

DRYAD_FILES = {
    "PEONY_PIGEON3_metadata_Nov2026.csv": (
        5041229, 69471, "fb474db3a2feb9a0f83fcd760aef2e7bcc3396fc9eab8809e46abb0d687eaf05",
    ),
    "PEONY1.3_Nov26.fa.gz": (
        5044872, 8046644960, "0107a228f09fafff66ddedcdaa13f23e42d50058e917515115915a9c269b2928",
    ),
    "PIGEON3.0.fa.gz": (
        5041224, 9125855360, "f12f09a59bf259384c595ae463450d1a886f11bcaf348f45285b3dc17d9cf5b5",
    ),
    "README.md": (
        5041227, 6292, "6e6178f64aa016a0eccf721977bd1ce2f8e4ef0310165249b66aafde648ab53c",
    ),
}

LRV1_TRIM_HEAD = [
    "source", "genome", "accession", "source_contig_uid", "contig_uid",
    "source_length", "repeat_length", "bases_removed", "length",
]
LRV1_RECORD_HEAD = [
    "ordinal", "source", "category", "accession", "contig_uid", "length",
    "dtr_class", "circular", "source_contig_uid", "source_length", "trimmed_bases", "owner_stage",
]
LRV1_META_HEAD = [
    "stage", "source", "category", "genome_id", "contig_uid", "accession",
    "length", "dtr_class", "circular", "trimmed_bases", "owner_stage",
]
V4_RECORD_HEAD = [
    "ordinal", "source", "genome", "contig_uid", "length",
    "dtr_class", "circular", "source_contig_uid", "source_length", "trimmed_bases", "owner_stage",
]
V4_LEFT_HEAD = ["source", "accession", "length", "reason"]


def _open_bin(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return open(path, "rb")


def _parse_owner_inputs(items: list[str]) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for item in items:
        if "=" not in item:
            raise Refused(f"--owner-inputs expects NAME=DIR, got {item!r}")
        name, path = item.split("=", 1)
        if not name or not path:
            raise Refused(f"--owner-inputs expects NAME=DIR, got {item!r}")
        out[name] = Path(path)
    return out


def apply_dec65_metagenomic(seq: bytes) -> tuple[bytes, int, str, bool]:
    """Apply DEC-65 (1) to a metagenome-assembled record:
    exact 20..1,000 bp non-low-complexity DTR -> trim trailing k bases and stage circular;
    all other classes (`none`, `partial`, `long`, `low_complexity`, `unresolved`) -> untrimmed linear.
    """
    e = v2.ends(seq)
    dtr_cls = v2.classify(e)
    if dtr_cls == "exact":
        k = e.exact
        return v2.trim(seq, k), k, "exact", True
    return seq, 0, dtr_cls, False


def _is_fasta_name(name: str) -> bool:
    low = name.lower()
    for suf in (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz"):
        if low.endswith(suf):
            return True
    return False


def _write_fasta_stream_records(raw_stream, out_fh, *, prefix: str = "") -> tuple[int, int]:
    """Stream FASTA records from `raw_stream` into `out_fh`, returning `(records, bases)`."""
    n_rec = 0
    n_bases = 0
    for r in v2.iter_fasta_stream(raw_stream):
        if not r.seq or not VALID_DNA_RE.match(r.seq):
            continue
        hd = f"{prefix}{r.header}" if prefix else r.header
        out_fh.write(b"\n".join(fasta.format_record(hd, r.seq, r.line_width or 70)) + b"\n")
        n_rec += 1
        n_bases += len(r.seq)
    return n_rec, n_bases


def extract_kavagutti_phages(support_tar_gz: Path, out_fa: Path) -> dict:
    """Extract all complete phage FASTAs from `phages-fna.tar.gz` inside Kavagutti's `Support-data.gz`."""
    support_tar_gz = Path(support_tar_gz)
    out_fa = Path(out_fa)
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_fa.with_name(f".{out_fa.name}.tmp")

    n_rec = 0
    n_bases = 0
    members_read = 0
    found_inner = None
    with tarfile.open(support_tar_gz, "r:*") as outer:
        for m in outer.getmembers():
            if m.isfile() and m.name.endswith("phages-fna.tar.gz"):
                found_inner = m
                break
        if found_inner is None:
            raise Refused(f"phages-fna.tar.gz not found inside {support_tar_gz}")
        inner_fobj = outer.extractfile(found_inner)
        if inner_fobj is None:
            raise Refused(f"could not read {found_inner.name} from {support_tar_gz}")
        with tarfile.open(fileobj=inner_fobj, mode="r:*") as inner, open(tmp, "wb") as w:
            inner_members = sorted(
                (im for im in inner.getmembers() if im.isfile() and _is_fasta_name(im.name)),
                key=lambda im: im.name,
            )
            for im in inner_members:
                fobj = inner.extractfile(im)
                if fobj is None:
                    continue
                stream = gzip.GzipFile(fileobj=fobj) if im.name.lower().endswith(".gz") else fobj
                rc, bs = _write_fasta_stream_records(stream, w)
                n_rec += rc
                n_bases += bs
                members_read += 1

    if n_rec == 0:
        tmp.unlink(missing_ok=True)
        raise Refused(f"no FASTA records extracted from {found_inner.name} in {support_tar_gz}")
    tmp.rename(out_fa)
    return {
        "inner_archive": found_inner.name,
        "members": members_read,
        "records": n_rec,
        "bases": n_bases,
        "bytes": out_fa.stat().st_size,
        "sha256": ids.sha256_file(out_fa).hex(),
    }


def _extract_fastas_from_zip(zf: zipfile.ZipFile, out_fh, *, include_filter=None) -> tuple[int, int, list[str]]:
    """Recursively extract FASTA records from a ZipFile (handles nested .zip, .tar.gz, .gz, .fa/.fna/.fasta)."""
    n_rec = 0
    n_bases = 0
    used: list[str] = []
    for name in sorted(zf.namelist()):
        if name.endswith("/") or "__MACOSX" in name:
            continue
        low = name.lower()
        if include_filter is not None and not include_filter(name):
            continue
        if low.endswith(".zip"):
            data = zf.read(name)
            with zipfile.ZipFile(io.BytesIO(data)) as sub:
                rc, bs, sub_used = _extract_fastas_from_zip(sub, out_fh)
                n_rec += rc
                n_bases += bs
                used.extend(f"{name}!{u}" for u in sub_used)
        elif low.endswith((".tar.gz", ".tgz", ".tar")):
            data = zf.read(name)
            with tarfile.open(fileobj=io.BytesIO(data), mode="r:*") as tf:
                for tm in sorted(tf.getmembers(), key=lambda m: m.name):
                    if not tm.isfile() or not _is_fasta_name(tm.name):
                        continue
                    fobj = tf.extractfile(tm)
                    if fobj is None:
                        continue
                    stream = gzip.GzipFile(fileobj=fobj) if tm.name.lower().endswith(".gz") else fobj
                    rc, bs = _write_fasta_stream_records(stream, out_fh)
                    n_rec += rc
                    n_bases += bs
                    used.append(f"{name}!{tm.name}")
        elif _is_fasta_name(name):
            with zf.open(name) as raw:
                stream = gzip.GzipFile(fileobj=raw) if low.endswith(".gz") else raw
                rc, bs = _write_fasta_stream_records(stream, out_fh)
                n_rec += rc
                n_bases += bs
                used.append(name)
    return n_rec, n_bases, used


def extract_virion2_fastas(supp_zip: Path, out_fa: Path) -> dict:
    """Extract viral assembly FASTAs from Europe PMC `PMC8018248_SupplementaryFiles.zip`."""
    supp_zip = Path(supp_zip)
    out_fa = Path(out_fa)
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_fa.with_name(f".{out_fa.name}.tmp")
    with zipfile.ZipFile(supp_zip) as zf, open(tmp, "wb") as w:
        n_rec, n_bases, used = _extract_fastas_from_zip(zf, w)
    if n_rec == 0:
        tmp.unlink(missing_ok=True)
        raise Refused(f"no FASTA records found inside {supp_zip}")
    tmp.rename(out_fa)
    return {
        "members": used,
        "records": n_rec,
        "bases": n_bases,
        "bytes": out_fa.stat().st_size,
        "sha256": ids.sha256_file(out_fa).hex(),
    }


# --------------------------------------------------------------------------------------------------------------
# Stage LRV1 (compact redistributable=True long-read marine & freshwater viromes; DEC-65 DTR trimming)


def stage_lrv1(
    sources_spec: list[dict],
    out: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
    final_dir: Path | None = None,
) -> dict:
    """Stage compact public `LRV1` inputs with DEC-65 exact 20..1,000 bp DTR trimming."""
    out = Path(out)
    if out.exists():
        raise Refused(f"output directory already exists: {out}")
    for spec in sources_spec:
        if not spec.get("redistributable", False):
            raise Refused(f"LRV1 requires redistributable: true, got {spec['name']!r}")

    work = out.with_name(out.name + ".tmp")
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    obj_root = Path(final_dir) if final_dir is not None else out

    trims: list[tuple] = []
    dtr_classes: Counter[str] = Counter()
    staged_recs: list[tuple] = []  # (src, cat, acc, cu, len, dtr_cls, circ, src_cu, src_len, k)
    objects: list[dict] = []
    skipped_invalid = 0
    src_bases = 0
    staged_bases = 0

    try:
        for spec in sources_spec:
            src_name = spec["name"]
            cat = spec["category"]
            in_fa = Path(spec["fasta"])
            out_name = spec["out_name"]
            out_fa = work / out_name
            n_src_rec = 0
            with _open_bin(in_fa) as fh, open(out_fa, "wb") as w:
                for r in v2.iter_fasta_stream(fh):
                    if not r.seq or not VALID_DNA_RE.match(r.seq):
                        skipped_invalid += 1
                        continue
                    src_len = len(r.seq)
                    src_bases += src_len
                    src_cu = ids.contig_uid(r.seq)
                    seq, k, dtr_cls, circ = apply_dec65_metagenomic(r.seq)
                    dtr_classes[dtr_cls] += 1
                    w.write(b"\n".join(fasta.format_record(r.header, seq, r.line_width or 70)) + b"\n")
                    staged_len = len(seq)
                    staged_bases += staged_len
                    cu = ids.contig_uid(seq) if k else src_cu
                    acc = r.accession
                    if k:
                        trims.append((
                            src_name, acc, acc, src_cu.hex(), cu.hex(),
                            src_len, k, k, staged_len,
                        ))
                    staged_recs.append((
                        src_name, cat, acc, cu, staged_len, dtr_cls, circ, src_cu, src_len, k,
                    ))
                    n_src_rec += 1
            if n_src_rec == 0:
                raise Refused(f"source {src_name!r} yielded 0 valid DNA records from {in_fa}")
            objects.append({
                "path": str(obj_root / out_name),
                "source_name": src_name,
                "source_release": spec["source_release"],
                "license": spec["license"],
                "redistributable": True,
                "retrieved_at": retrieved,
            })

        root_stages_list = list(root_stages)
        all_uids = [r[3] for r in staged_recs]
        owners = v2.owners_of(all_uids, root, root_stages_list, owner_inputs)

        listing = make_inputs.Listing()
        genomes: list[dict] = []
        rec_rows: list[tuple] = []
        meta_rows: list[tuple] = []
        repeated: list[tuple] = []
        shared: list[tuple] = []
        first_of: dict[bytes, str] = {}
        contig_len: dict[bytes, int] = {}
        by_source: Counter[str] = Counter()
        by_source_bases: Counter[str] = Counter()

        for ordinal, (src_name, cat, acc, cu, staged_len, dtr_cls, circ, src_cu, src_len, k) in enumerate(staged_recs):
            owner = owners.get(cu, "")
            rec_rows.append((
                ordinal, src_name, cat, acc, cu.hex(), staged_len,
                dtr_cls, str(circ).lower(), src_cu.hex(), src_len, k, owner,
            ))
            if owner:
                shared.append((ordinal, acc, cu.hex(), owner))
            contig_len.setdefault(cu, staged_len)
            g = listing.add([cu], 11, "virus", "complete")
            guid = ids.genome_uid([cu])
            if g is None:
                repeated.append((ordinal, acc, guid.hex(), first_of[guid]))
                continue
            g["circular"] = [cu.hex()] if circ else []
            genomes.append(g)
            first_of[guid] = acc
            by_source[src_name] += 1
            by_source_bases[src_name] += staged_len
            meta_rows.append((
                "LRV1", src_name, cat, guid.hex(), cu.hex(), acc,
                staged_len, dtr_cls, str(circ).lower(), k, owner,
            ))

        trims.sort(key=lambda t: (t[0], t[1], t[2]))
        make_inputs._write_tsv(work / "records.tsv", LRV1_RECORD_HEAD, rec_rows)
        make_inputs._write_tsv(work / "lrv1_metadata.tsv", LRV1_META_HEAD, meta_rows)
        make_inputs._write_tsv(work / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
        make_inputs._write_tsv(work / "shared_records.tsv", ["ordinal", "genome", "contig_uid", "owner_stage"], shared)
        make_inputs._write_tsv(work / "trims.tsv", LRV1_TRIM_HEAD, trims)
        (work / "objects.json").write_text(json.dumps(objects, indent=1))
        (work / "genomes.json").write_text(json.dumps(genomes))

        owner_contigs: Counter[str] = Counter()
        owner_bases: Counter[str] = Counter()
        for cu, n in contig_len.items():
            if cu in owners:
                st = owners[cu]
                owner_contigs[st] += 1
                owner_bases[st] += n

        held_contigs = sum(owner_contigs.values())
        held_bases = sum(owner_bases.values())
        work.rename(out)
    except Exception:
        if work.exists():
            shutil.rmtree(work)
        raise

    return {
        "stage": "LRV1",
        "total_records": len(staged_recs),
        "skipped_invalid": skipped_invalid,
        "source_bases": src_bases,
        "staged_bases": staged_bases,
        "dtr_classes": dict(sorted(dtr_classes.items())),
        "trims": len(trims),
        "trimmed_bases": sum(t[7] for t in trims),
        "genomes": len(genomes),
        "contigs": len(contig_len),
        "contig_bases": sum(contig_len.values()),
        "circular_genomes": sum(1 for g in genomes if g["circular"]),
        "lt_5kb_contigs": sum(1 for n in contig_len.values() if n < 5000),
        "repeated_genomes": len(repeated),
        "shared_records": len(shared),
        "by_source_genomes": dict(sorted(by_source.items())),
        "by_source_bases": dict(sorted(by_source_bases.items())),
        "owner_contigs": dict(sorted(owner_contigs.items())),
        "owner_bases": dict(sorted(owner_bases.items())),
        "held_contigs": held_contigs,
        "held_bases": held_bases,
        "new_contigs": len(contig_len) - held_contigs,
        "new_bases": sum(contig_len.values()) - held_bases,
        "root_stages": root_stages_list,
        "owner_inputs": list(owner_inputs),
        "objects": len(objects),
        **listing.report(),
    }


# --------------------------------------------------------------------------------------------------------------
# Stage V4 (PEONY 1.3 / PIGEON 3.0 streaming staging with Part 1 split layout)


def stage_v4(
    fasta_inputs: list[tuple[str, Path]],
    *,
    fasta_out: Path,
    meta_out: Path,
    final_dir: Path,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
) -> dict:
    """Stream `PEONY 1.3` / `PIGEON 3.0` `.fa.gz` files into a consolidated FASTA on `/stor/scratch`
    (`fasta_out`) and small metadata files in `/tmp` (`meta_out`), applying DEC-65 exact 20..1,000 bp
    DTR trimming and deduplicating `contig_uid` across inputs and earlier stages.
    """
    fasta_out = Path(fasta_out)
    meta_out = Path(meta_out)
    final_dir = Path(final_dir)
    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    meta_out.mkdir(parents=True, exist_ok=True)

    tmp_fa = fasta_out.with_name(f".{fasta_out.name}.tmp")
    trims: list[tuple] = []
    dtr_classes: Counter[str] = Counter()
    left_rows: list[tuple] = []
    staged_recs: list[tuple] = []  # (src, acc, cu, len, dtr_cls, circ, src_cu, src_len, k)
    src_bases = 0
    staged_bases = 0
    seen_written_cu: set[bytes] = set()

    with open(tmp_fa, "wb") as fa_w:
        for src_label, in_path in fasta_inputs:
            in_path = Path(in_path)
            with _open_bin(in_path) as fh:
                for r in v2.iter_fasta_stream(fh):
                    acc = r.accession
                    if not r.seq or not VALID_DNA_RE.match(r.seq):
                        left_rows.append((src_label, acc, len(r.seq), "invalid_or_empty_dna"))
                        continue
                    src_len = len(r.seq)
                    src_bases += src_len
                    src_cu = ids.contig_uid(r.seq)
                    seq, k, dtr_cls, circ = apply_dec65_metagenomic(r.seq)
                    dtr_classes[dtr_cls] += 1
                    staged_len = len(seq)
                    cu = ids.contig_uid(seq) if k else src_cu
                    if k:
                        trims.append((
                            src_label, acc, acc, src_cu.hex(), cu.hex(),
                            src_len, k, k, staged_len,
                        ))
                    if cu not in seen_written_cu:
                        fa_w.write(b"\n".join(fasta.format_record(r.header, seq, r.line_width or 70)) + b"\n")
                        seen_written_cu.add(cu)
                        staged_bases += staged_len
                    staged_recs.append((
                        src_label, acc, cu, staged_len, dtr_cls, circ, src_cu, src_len, k,
                    ))

    if not staged_recs:
        tmp_fa.unlink(missing_ok=True)
        raise Refused("no valid DNA records extracted for V4")
    tmp_fa.rename(fasta_out)

    root_stages_list = list(root_stages)
    all_uids = [r[2] for r in staged_recs]
    owners = v2.owners_of(all_uids, root, root_stages_list, owner_inputs)

    listing = make_inputs.Listing()
    genomes: list[dict] = []
    rec_rows: list[tuple] = []
    repeated: list[tuple] = []
    shared: list[tuple] = []
    first_of: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}

    for ordinal, (src_label, acc, cu, staged_len, dtr_cls, circ, src_cu, src_len, k) in enumerate(staged_recs):
        owner = owners.get(cu, "")
        rec_rows.append((
            ordinal, src_label, acc, cu.hex(), staged_len,
            dtr_cls, str(circ).lower(), src_cu.hex(), src_len, k, owner,
        ))
        if owner:
            shared.append((ordinal, acc, cu.hex(), owner))
        contig_len.setdefault(cu, staged_len)
        g = listing.add([cu], 11, "virus", "vOTU")
        guid = ids.genome_uid([cu])
        if g is None:
            repeated.append((ordinal, acc, guid.hex(), first_of[guid]))
            continue
        g["circular"] = [cu.hex()] if circ else []
        genomes.append(g)
        first_of[guid] = acc

    trims.sort(key=lambda t: (t[0], t[1], t[2]))
    make_inputs._write_tsv(meta_out / "records.tsv", V4_RECORD_HEAD, rec_rows)
    make_inputs._write_tsv(meta_out / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
    make_inputs._write_tsv(meta_out / "shared_records.tsv", ["ordinal", "genome", "contig_uid", "owner_stage"], shared)
    make_inputs._write_tsv(meta_out / "trims.tsv", LRV1_TRIM_HEAD, trims)

    left_path = meta_out / "V4_left_for_later.tsv.gz"
    with gzip.open(left_path, "wt", encoding="utf-8", newline="") as lf:
        w = csv.writer(lf, delimiter="\t", lineterminator="\n")
        w.writerow(V4_LEFT_HEAD)
        for row in left_rows:
            w.writerow(row)

    base_obj = {
        "source_name": "peony_pigeon3_dryad",
        "source_release": "PEONY 1.3 + PIGEON 3.0 (Dryad doi:10.5061/dryad.prr4xgz35, v7 2026-11-06)",
        "license": "CC0-1.0",
        "redistributable": True,
        "retrieved_at": retrieved,
    }
    objects_final = [{"path": str(final_dir / fasta_out.name), **base_obj}]
    objects_staging = [{"path": str(fasta_out), **base_obj}]
    (meta_out / "objects.json").write_text(json.dumps(objects_final, indent=1))
    (meta_out / "objects_staging.json").write_text(json.dumps(objects_staging, indent=1))
    (meta_out / "genomes.json").write_text(json.dumps(genomes))

    owner_contigs: Counter[str] = Counter()
    owner_bases: Counter[str] = Counter()
    for cu, n in contig_len.items():
        if cu in owners:
            st = owners[cu]
            owner_contigs[st] += 1
            owner_bases[st] += n

    held_contigs = sum(owner_contigs.values())
    held_bases = sum(owner_bases.values())

    return {
        "stage": "V4",
        "total_records": len(staged_recs),
        "left_for_later_rows": len(left_rows),
        "source_bases": src_bases,
        "staged_bases": staged_bases,
        "fasta_bytes": fasta_out.stat().st_size,
        "fasta_sha256": ids.sha256_file(fasta_out).hex(),
        "dtr_classes": dict(sorted(dtr_classes.items())),
        "trims": len(trims),
        "trimmed_bases": sum(t[7] for t in trims),
        "genomes": len(genomes),
        "contigs": len(contig_len),
        "contig_bases": sum(contig_len.values()),
        "circular_genomes": sum(1 for g in genomes if g["circular"]),
        "repeated_genomes": len(repeated),
        "shared_records": len(shared),
        "owner_contigs": dict(sorted(owner_contigs.items())),
        "owner_bases": dict(sorted(owner_bases.items())),
        "held_contigs": held_contigs,
        "held_bases": held_bases,
        "new_contigs": len(contig_len) - held_contigs,
        "new_bases": sum(contig_len.values()) - held_bases,
        "root_stages": root_stages_list,
        "owner_inputs": list(owner_inputs),
        "objects": 1,
        **listing.report(),
    }


# --------------------------------------------------------------------------------------------------------------
# Polite public fetchers (User-Agent: curl, <= 2 req/s, stop on 5 consecutive refusals; DEC-50 / AGENTS §8)


def _http_get(
    url: str,
    *,
    min_interval: float = 0.55,
    state: dict | None = None,
    timeout: int = 60,
    headers: dict[str, str] | None = None,
) -> bytes:
    if state is None:
        state = {"last": 0.0, "refusals": 0}
    wait = min_interval - (time.monotonic() - state["last"])
    if wait > 0:
        time.sleep(wait)
    hdrs = {"User-Agent": "curl"}
    if headers:
        hdrs.update(headers)
    req = urllib.request.Request(url, headers=hdrs)
    for attempt in range(3):
        state["last"] = time.monotonic()
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            state["refusals"] = 0
            return data
        except urllib.error.HTTPError as e:
            if e.code in (403, 429, 500, 502, 503, 504):
                state["refusals"] += 1
                if state["refusals"] >= 5:
                    raise Refused(f"stopped after 5 consecutive HTTP refusals (last {e.code} on {url})") from e
                time.sleep(2.0 * (attempt + 1))
                continue
            raise
        except (urllib.error.URLError, TimeoutError) as e:
            state["refusals"] += 1
            if state["refusals"] >= 5:
                raise Refused(f"stopped after 5 consecutive network errors (last {e} on {url})") from e
            time.sleep(2.0 * (attempt + 1))
    raise Refused(f"failed after 3 attempts: {url}")


def _download_file(
    url: str,
    dst: Path,
    *,
    state: dict,
    expected_md5: str | None = None,
    expected_sha256: str | None = None,
    headers: dict[str, str] | None = None,
) -> tuple[int, str, str]:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.is_file() and dst.stat().st_size > 0:
        h_md5 = hashlib.md5()
        h_sha = hashlib.sha256()
        n = 0
        with open(dst, "rb") as fh:
            while True:
                b = fh.read(1 << 20)
                if not b:
                    break
                h_md5.update(b)
                h_sha.update(b)
                n += len(b)
        md5_hex, sha_hex = h_md5.hexdigest(), h_sha.hexdigest()
        if expected_md5 and md5_hex != expected_md5:
            raise Refused(f"{dst}: cached md5 {md5_hex} != expected {expected_md5}")
        if expected_sha256 and sha_hex != expected_sha256:
            raise Refused(f"{dst}: cached sha256 {sha_hex} != expected {expected_sha256}")
        return n, md5_hex, sha_hex

    tmp = dst.with_name(f".{dst.name}.tmp")
    hdrs = {"User-Agent": "curl"}
    if headers:
        hdrs.update(headers)
    req = urllib.request.Request(url, headers=hdrs)
    wait = 0.55 - (time.monotonic() - state["last"])
    if wait > 0:
        time.sleep(wait)
    state["last"] = time.monotonic()
    h_md5 = hashlib.md5()
    h_sha = hashlib.sha256()
    n = 0
    with urllib.request.urlopen(req, timeout=600) as resp, open(tmp, "wb") as fh:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
            h_md5.update(chunk)
            h_sha.update(chunk)
            n += len(chunk)
    md5_hex, sha_hex = h_md5.hexdigest(), h_sha.hexdigest()
    if expected_md5 and md5_hex != expected_md5:
        tmp.unlink(missing_ok=True)
        raise Refused(f"{dst}: downloaded md5 {md5_hex} != expected {expected_md5}")
    if expected_sha256 and sha_hex != expected_sha256:
        tmp.unlink(missing_ok=True)
        raise Refused(f"{dst}: downloaded sha256 {sha_hex} != expected {expected_sha256}")
    tmp.rename(dst)
    state["refusals"] = 0
    return n, md5_hex, sha_hex


def _fetch_beaulaurier_prjna529454(out_dir: Path, *, state: dict, batch_size: int = 200) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    fa_path = out_dir / "beaulaurier_aloha_afvgs.fna"
    meta_path = out_dir / "beaulaurier_aloha_afvgs_meta.tsv"
    if fa_path.is_file() and meta_path.is_file():
        return {
            "fasta": str(fa_path),
            "bytes": fa_path.stat().st_size,
            "sha256": ids.sha256_file(fa_path).hex(),
            "cached": True,
        }

    term = "PRJNA529454[BioProject]"
    uids: list[str] = []
    retstart = 0
    while True:
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore"
            f"&retmax=10000&retstart={retstart}&retmode=json&term={urllib.parse.quote(term)}"
        )
        d = json.loads(_http_get(url, state=state))["esearchresult"]
        got = d.get("idlist") or []
        uids.extend(got)
        total = int(d.get("count") or 0)
        retstart += len(got)
        if not got or retstart >= total:
            break

    meta_rows: list[tuple] = []
    kept_uids: list[str] = []
    for i in range(0, len(uids), batch_size):
        chunk = uids[i:i + batch_size]
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore"
            f"&retmode=json&id={','.join(chunk)}"
        )
        res = json.loads(_http_get(url, state=state)).get("result") or {}
        for uid in chunk:
            info = res.get(uid)
            if not isinstance(info, dict):
                continue
            mol = (info.get("moltype") or "").lower()
            if mol and "dna" not in mol:
                continue
            acc = info.get("accessionversion") or info.get("caption") or ""
            if not acc:
                continue
            topo = (info.get("topology") or "").lower()
            slen = int(info.get("slen") or 0)
            title = (info.get("title") or "").replace("\t", " ").replace("\n", " ")
            kept_uids.append(uid)
            meta_rows.append((uid, acc, slen, topo, title))

    tmp_fa = fa_path.with_name(f".{fa_path.name}.tmp")
    with open(tmp_fa, "wb") as w:
        for i in range(0, len(kept_uids), batch_size):
            chunk = kept_uids[i:i + batch_size]
            url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore"
                f"&rettype=fasta&retmode=text&id={','.join(chunk)}"
            )
            data = _http_get(url, state=state, timeout=180)
            w.write(data)
            if data and not data.endswith(b"\n"):
                w.write(b"\n")
    tmp_fa.rename(fa_path)
    make_inputs._write_tsv(meta_path, ["uid", "accession", "slen", "topology", "title"], meta_rows)
    return {
        "esearch_uids": len(uids),
        "kept_uids": len(kept_uids),
        "fasta": str(fa_path),
        "bytes": fa_path.stat().st_size,
        "sha256": ids.sha256_file(fa_path).hex(),
        "cached": False,
    }


def _cmd_fetch_lrv1(a) -> int:
    dl = Path(a.dl_dir)
    dl.mkdir(parents=True, exist_ok=True)
    state = {"last": 0.0, "refusals": 0}
    rep: dict = {}

    # 1. Beaulaurier et al. 2020 PRJNA529454
    rep["beaulaurier"] = _fetch_beaulaurier_prjna529454(dl / "beaulaurier", state=state)

    # 2. Kavagutti et al. 2023 Figshare 20431332 Support-data.gz -> phages-fna.tar.gz
    kava_dir = dl / "kavagutti"
    kava_tar = kava_dir / "Support-data.gz"
    kb, kmd5, ksha = _download_file(KAVAGUTTI_URL, kava_tar, state=state, expected_md5=KAVAGUTTI_MD5)
    kava_fa = kava_dir / "kavagutti_freshwater_hifi_phages.fna"
    kava_ext = extract_kavagutti_phages(kava_tar, kava_fa)
    rep["kavagutti"] = {"archive_bytes": kb, "archive_md5": kmd5, "archive_sha256": ksha, **kava_ext}

    # 3. VirION2 PMC8018248 SupplementaryFiles.zip (optional; PeerJ behind Cloudflare, PMC lacks FASTA)
    if not getattr(a, "skip_virion2", False):
        v2_dir = dl / "virion2"
        v2_zip = v2_dir / "PMC8018248_SupplementaryFiles.zip"
        vb, vmd5, vsha = _download_file(VIRION2_URL, v2_zip, state=state)
        v2_fa = v2_dir / "virion2_wec_marine_virome.fna"
        v2_ext = extract_virion2_fastas(v2_zip, v2_fa)
        rep["virion2"] = {"archive_bytes": vb, "archive_md5": vmd5, "archive_sha256": vsha, **v2_ext}

    print(json.dumps(rep, indent=2, sort_keys=True))
    return 0


def _cmd_stage_lrv1(a) -> int:
    dl = Path(a.dl_dir)
    sources_spec = [
        {
            "name": "beaulaurier_aloha_afvgs",
            "category": "ont_marine_afvg",
            "fasta": str(dl / "beaulaurier" / "beaulaurier_aloha_afvgs.fna"),
            "out_name": "beaulaurier_aloha_afvgs.fna",
            "source_release": "Beaulaurier et al. 2020 Station ALOHA ONT AFVGs (PRJNA529454, Genome Research)",
            "license": "INSDC / NCBI GenBank (public domain)",
            "redistributable": True,
        },
        {
            "name": "kavagutti_freshwater_hifi_phages",
            "category": "hifi_freshwater_phage",
            "fasta": str(dl / "kavagutti" / "kavagutti_freshwater_hifi_phages.fna"),
            "out_name": "kavagutti_freshwater_hifi_phages.fna",
            "source_release": "Kavagutti et al. 2023 freshwater HiFi phages (Figshare doi:10.6084/m9.figshare.20431332.v1)",
            "license": "CC BY 4.0",
            "redistributable": True,
        },
    ]
    v2_fa = dl / "virion2" / "virion2_wec_marine_virome.fna"
    if v2_fa.is_file() and v2_fa.stat().st_size > 0:
        sources_spec.append({
            "name": "virion2_wec_marine_virome",
            "category": "ont_marine_virion2",
            "fasta": str(v2_fa),
            "out_name": "virion2_wec_marine_virome.fna",
            "source_release": "Zablocki et al. 2021 VirION2 WEC marine virome (PeerJ 9:e11088, PMC8018248)",
            "license": "CC BY 4.0",
            "redistributable": True,
        })
    root = Path(a.root) if a.root else None
    root_stages = [s for s in (a.root_stages.split(",") if a.root_stages else []) if s]
    owner_inputs = _parse_owner_inputs(a.owner_inputs or [])
    rep = stage_lrv1(
        sources_spec,
        Path(a.out),
        root=root,
        root_stages=root_stages,
        owner_inputs=owner_inputs,
        retrieved=a.retrieved,
        final_dir=Path(a.final_dir) if a.final_dir else None,
    )
    print(json.dumps(rep, indent=2, sort_keys=True))
    return 0


def _cmd_fetch_v4(a) -> int:
    token = os.environ.get("DRYAD_TOKEN", "").strip()
    if a.token_file:
        token = Path(a.token_file).read_text().strip()
    if not token:
        raise Refused(
            "Dryad API v2 requires an ORCID Bearer Token (set DRYAD_TOKEN or pass --token-file)"
        )
    dl = Path(a.dl_dir)
    dl.mkdir(parents=True, exist_ok=True)
    state = {"last": 0.0, "refusals": 0}
    headers = {"Authorization": f"Bearer {token}"}
    rep: dict = {}
    for fname, (fid, fbytes, fsha) in DRYAD_FILES.items():
        url = f"https://datadryad.org/api/v2/files/{fid}/download"
        n, md5_hex, sha_hex = _download_file(
            url, dl / fname, state=state, expected_sha256=fsha, headers=headers,
        )
        rep[fname] = {"file_id": fid, "bytes": n, "expected_bytes": fbytes, "md5": md5_hex, "sha256": sha_hex}
    print(json.dumps(rep, indent=2, sort_keys=True))
    return 0


def _cmd_stage_v4(a) -> int:
    dl = Path(a.dl_dir)
    fasta_inputs: list[tuple[str, Path]] = []
    peony = dl / "PEONY1.3_Nov26.fa.gz"
    pigeon = dl / "PIGEON3.0.fa.gz"
    if peony.is_file():
        fasta_inputs.append(("peony1.3", peony))
    if pigeon.is_file():
        fasta_inputs.append(("pigeon3.0", pigeon))
    if not fasta_inputs:
        raise Refused(f"neither PEONY1.3_Nov26.fa.gz nor PIGEON3.0.fa.gz found in {dl}")
    root = Path(a.root) if a.root else None
    root_stages = [s for s in (a.root_stages.split(",") if a.root_stages else []) if s]
    owner_inputs = _parse_owner_inputs(a.owner_inputs or [])
    rep = stage_v4(
        fasta_inputs,
        fasta_out=Path(a.fasta_out),
        meta_out=Path(a.meta_out),
        final_dir=Path(a.final_dir),
        root=root,
        root_stages=root_stages,
        owner_inputs=owner_inputs,
        retrieved=a.retrieved,
    )
    print(json.dumps(rep, indent=2, sort_keys=True))
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Card C-67: LRV1 & V4 (PEONY/PIGEON 3.0) staging")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_fl = sub.add_parser("fetch-lrv1", help="Fetch Beaulaurier, Kavagutti, and VirION2 sources")
    p_fl.add_argument("--dl-dir", required=True)
    p_fl.add_argument("--skip-virion2", action="store_true")
    p_fl.set_defaults(func=_cmd_fetch_lrv1)

    p_sl = sub.add_parser("stage-lrv1", help="Stage LRV1 into /tmp/if2453/s2/inputs/LRV1")
    p_sl.add_argument("--dl-dir", required=True)
    p_sl.add_argument("--out", required=True)
    p_sl.add_argument("--final-dir", default="")
    p_sl.add_argument("--root", default="")
    p_sl.add_argument("--root-stages", default="")
    p_sl.add_argument("--owner-inputs", action="append", default=[])
    p_sl.add_argument("--retrieved", required=True)
    p_sl.set_defaults(func=_cmd_stage_lrv1)

    p_fv = sub.add_parser("fetch-v4", help="Fetch PEONY 1.3 / PIGEON 3.0 from Dryad API v2 (requires bearer token)")
    p_fv.add_argument("--dl-dir", required=True)
    p_fv.add_argument("--token-file", default="")
    p_fv.set_defaults(func=_cmd_fetch_v4)

    p_sv = sub.add_parser("stage-v4", help="Stream and stage V4 (PEONY 1.3 + PIGEON 3.0) in Part-1 split layout")
    p_sv.add_argument("--dl-dir", required=True)
    p_sv.add_argument("--fasta-out", required=True)
    p_sv.add_argument("--meta-out", required=True)
    p_sv.add_argument("--final-dir", default="/tmp/if2453/s2/inputs/V4")
    p_sv.add_argument("--root", default="")
    p_sv.add_argument("--root-stages", default="")
    p_sv.add_argument("--owner-inputs", action="append", default=[])
    p_sv.add_argument("--retrieved", required=True)
    p_sv.set_defaults(func=_cmd_stage_v4)

    args = ap.parse_args(argv)
    try:
        return args.func(args)
    except Refused as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
