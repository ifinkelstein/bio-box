"""Card C-69 (DEC-59, DEC-65, DEC-66, DEC-69, DEC-50, DEC-51, DEC-26, DEC-27):
fetch, deduplicate, and stage *Egovirales* gut giant-virus genomes (`EGV1`) from
Gaïa et al. 2026 (*Nature*, `doi:10.1038/s41586-026-11283-z`), Figshare
`10.6084/m9.figshare.25398895.v14` Item 08 (`ITEM_08_Egovirales_Genomes.zip`,
224 environmental genomes, ~51.3 Mbp, `CC BY 4.0`, `redistributable: true`).

Hard rules enforced:
- Source must be `redistributable: true` (`DEC-51`).
- Every genome is staged byte-identically as deposited (`DEC-65`):
  - `kind="virus"`, `genetic_code=11` (`Prodigal -p meta` / `gv-meta` mode).
  - `topology="linear"` (`circular: []`), NEVER trimmed (`trims.tsv` always 0 rows),
    preserving the 55 bp terminal inverted repeats (ITRs) of complete linear genomes.
  - Exact ITR length (`>= 20` bp) is detected read-only and recorded in `egv1_metadata.tsv`
    and `genome_input.tsv`.
  - No `5 kb` contig length floor: every valid DNA contig is kept as deposited.
- Deterministic outputs regardless of archive member order or filesystem mtimes (`AGENTS.md` §4.3).
- Within-stage `contig_uid` / `genome_uid` deduplication (`DEC-27`) and cross-stage
  build-order ownership (`DEC-26`) against earlier stages (`S1..GART2, GV1`).
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
import shutil
import sys
import tarfile
import time
import urllib.error
import urllib.request
import zipfile
from collections import Counter
from collections.abc import Iterable, Iterator
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # c48c_virus_inputs and make_inputs beside this script

import c48c_virus_inputs as v2  # noqa: E402
import make_inputs  # noqa: E402
from biobox import fasta, ids  # noqa: E402

Refused = make_inputs.Refused

# Verified via Figshare API v2 (https://api.figshare.com/v2/articles/25398895) on 2026-11-24
FIGSHARE_ARTICLE_ID = 25398895
FIGSHARE_DOI = "10.6084/m9.figshare.25398895.v14"
ITEM08_FILE_ID = 66974843
ITEM08_FILENAME = "ITEM_08_Egovirales_Genomes.zip"
ITEM08_URL = "https://ndownloader.figshare.com/files/66974843"
ITEM08_BYTES = 24_753_281
ITEM08_MD5 = "1288ed17163e96afe6332523ff9ac94b"

EGV1_FASTA = "egovirales_gaia_2026.fna"
EGV1_SOURCE = {
    "source_name": "egovirales_gaia_2026",
    "source_release": (
        "figshare:10.6084/m9.figshare.25398895.v14 "
        "(ITEM_08_Egovirales_Genomes.zip; Gaia et al. 2026, Nature doi:10.1038/s41586-026-11283-z)"
    ),
    "license": "CC BY 4.0",
    "redistributable": True,
}

VALID_DNA_RE = re.compile(rb"^[ACGTURYKMSWBDHVNacgturykmswbdhvn]+$")
_COMP_TABLE = bytes.maketrans(
    b"ACGTacgtRYKMSWBDHVNrykmswbdhvn",
    b"TGCAtgcaYRMKWSVHDBNyrmkwsvhdbn",
)

FASTA_SUFFIXES = (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")

RECORD_HEAD = [
    "ordinal", "genome", "record", "accession", "contig_uid", "length",
    "itr_length", "member", "owner_stage",
]
GENOME_INPUT_HEAD = [
    "stage", "genome_uid", "genome_name", "contig_uid", "kind", "genetic_code",
    "topology", "circular", "assembly_quality", "length", "itr_length", "owner_stage",
]
METADATA_HEAD = [
    "stage", "source_name", "genome_name", "genome_uid", "contig_count", "total_bases",
    "min_contig_len", "max_contig_len", "topology", "assembly_quality", "max_itr_length",
    "contigs_with_itr", "held_contigs", "member",
]
TRIM_HEAD = [
    "source", "genome", "accession", "source_contig_uid", "contig_uid",
    "source_length", "repeat_length", "bases_removed", "length",
]
SKIPPED_MEMBER_HEAD = ["member", "bytes", "reason"]
SKIPPED_RECORD_HEAD = ["genome", "accession", "reason", "length"]


# --------------------------------------------------------------------------------------------------------------
# Pure sequence, archive, and DEC-26 prefix-overlap helpers


def detect_itr(seq: bytes, min_len: int = 20, max_len: int = 10_000) -> int:
    """Detect exact terminal inverted repeat (ITR) length at the ends of a linear viral genome.

    Compares `seq[i]` against the Watson-Crick complement of `seq[n - 1 - i]` from `i = 0`
    inward up to `min(max_len, len(seq) // 2)`. Returns the exact match length when
    `>= min_len`, else 0. Read-only: NEVER trims or modifies `seq` (`DEC-65`).
    """
    n = len(seq)
    if n < 2 * min_len:
        return 0
    limit = min(max_len, n // 2)
    head = seq[:limit].upper()
    tail_rc = seq[n - limit :].upper().translate(_COMP_TABLE)[::-1]
    k = 0
    while k < limit and head[k] == tail_rc[k]:
        k += 1
    return k if k >= min_len else 0


def _is_fasta_name(name: str) -> bool:
    base = name.rsplit("/", 1)[-1]
    if not base or base.startswith(".") or "__MACOSX" in name.split("/"):
        return False
    low = base.lower()
    return any(low.endswith(suf) for suf in FASTA_SUFFIXES)


def _genome_name_from_member(member: str) -> str:
    base = member.rsplit("/", 1)[-1]
    low = base.lower()
    for suf in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if low.endswith(suf):
            return base[: -len(suf)]
    return base


def _md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _iter_tar_fasta_members(
    tf: tarfile.TarFile,
    prefix: str,
    skipped_members: list[tuple[str, int, str]],
) -> Iterator[tuple[str, str, bytes]]:
    members = sorted((m for m in tf.getmembers() if m.isfile()), key=lambda m: m.name)
    for m in members:
        full_name = f"{prefix}!{m.name}" if prefix else m.name
        if not _is_fasta_name(m.name):
            skipped_members.append((full_name, m.size, "not_fasta"))
            continue
        fobj = tf.extractfile(m)
        if fobj is None:
            skipped_members.append((full_name, m.size, "unreadable"))
            continue
        raw = fobj.read()
        if m.name.lower().endswith(".gz"):
            raw = gzip.decompress(raw)
        yield _genome_name_from_member(m.name), full_name, raw


def _iter_zip_fasta_members(
    zf: zipfile.ZipFile,
    prefix: str,
    skipped_members: list[tuple[str, int, str]],
) -> Iterator[tuple[str, str, bytes]]:
    infos = sorted((i for i in zf.infolist() if not i.is_dir()), key=lambda i: i.filename)
    for info in infos:
        name = info.filename
        full_name = f"{prefix}!{name}" if prefix else name
        base = name.rsplit("/", 1)[-1]
        if not base or base.startswith(".") or "__MACOSX" in name.split("/"):
            skipped_members.append((full_name, info.file_size, "ignored_system_entry"))
            continue
        low = name.lower()
        if low.endswith(".zip"):
            with zipfile.ZipFile(io.BytesIO(zf.read(info))) as sub:
                yield from _iter_zip_fasta_members(sub, full_name, skipped_members)
        elif low.endswith((".tar.gz", ".tgz", ".tar")):
            with tarfile.open(fileobj=io.BytesIO(zf.read(info)), mode="r:*") as tf:
                yield from _iter_tar_fasta_members(tf, full_name, skipped_members)
        elif _is_fasta_name(name):
            raw = zf.read(info)
            if low.endswith(".gz"):
                raw = gzip.decompress(raw)
            yield _genome_name_from_member(name), full_name, raw
        else:
            skipped_members.append((full_name, info.file_size, "not_fasta"))


def _iter_genome_fastas(
    src: Path,
    skipped_members: list[tuple[str, int, str]],
) -> Iterator[tuple[str, str, bytes]]:
    """Yield `(genome_name, member_label, fasta_bytes)` deterministically from a zip, tar, directory, or FASTA."""
    src = Path(src)
    if src.is_dir():
        files = sorted((p for p in src.rglob("*") if p.is_file()), key=lambda p: p.relative_to(src).as_posix())
        for p in files:
            rel = p.relative_to(src).as_posix()
            if not _is_fasta_name(rel):
                skipped_members.append((rel, p.stat().st_size, "not_fasta"))
                continue
            raw = p.read_bytes()
            if rel.lower().endswith(".gz"):
                raw = gzip.decompress(raw)
            yield _genome_name_from_member(rel), rel, raw
        return

    if not src.is_file():
        raise Refused(f"input source not found: {src}")

    low = src.name.lower()
    if low.endswith(".zip"):
        with zipfile.ZipFile(src) as zf:
            yield from _iter_zip_fasta_members(zf, "", skipped_members)
    elif low.endswith((".tar.gz", ".tgz", ".tar")):
        with tarfile.open(src, "r:*") as tf:
            yield from _iter_tar_fasta_members(tf, "", skipped_members)
    elif _is_fasta_name(src.name):
        raw = src.read_bytes()
        if low.endswith(".gz"):
            raw = gzip.decompress(raw)
        yield _genome_name_from_member(src.name), src.name, raw
    else:
        raise Refused(f"unsupported input format: {src}")


def compute_prefix_overlap(
    contig_uids: Iterable[bytes],
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
) -> dict[bytes, str]:
    """Resolve `DEC-26` first-owner stage for each `contig_uid` across `root_stages` and staged `owner_inputs`."""
    return v2.owners_of(list(contig_uids), root, list(root_stages), owner_inputs)


# --------------------------------------------------------------------------------------------------------------
# Stage EGV1


def stage_egv1(
    src: Path,
    out: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
    final_dir: Path | None = None,
    expected_md5: str | None = None,
    redistributable: bool = True,
    stage_name: str = "EGV1",
) -> dict:
    """Stage `EGV1` (*Egovirales* 224 genomes) byte-identically as linear untrimmed `virus` genomes (`DEC-65`)."""
    if not redistributable:
        raise Refused(f"{stage_name}: refusing non-redistributable source (DEC-51)")
    src = Path(src)
    out = Path(out)
    if out.exists():
        raise Refused(f"{out} exists")
    if not src.exists():
        raise Refused(f"{stage_name}: input path does not exist: {src}")

    got_md5: str | None = None
    if src.is_file():
        got_md5 = _md5_file(src)
        if expected_md5 is not None and got_md5 != expected_md5:
            raise Refused(f"{src.name}: md5 mismatch (got {got_md5}, want {expected_md5})")

    work = out.with_name(f".{out.name}.tmp")
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    obj_dir = Path(final_dir) if final_dir is not None else out

    try:
        rep = _stage_egv1_work(
            src=src,
            work=work,
            obj_dir=obj_dir,
            stage_name=stage_name,
            root=root,
            root_stages=list(root_stages),
            owner_inputs=owner_inputs,
            retrieved=retrieved,
            source_md5=got_md5,
            expected_md5=expected_md5,
        )
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise

    work.rename(out)
    sum_files = sorted(
        p.name for p in out.iterdir()
        if p.is_file() and p.name not in ("inputs_SUMS", "stage_report.json")
    )
    sha_map = {n: ids.sha256_file(out / n).hex() for n in sum_files}
    bytes_map = {n: (out / n).stat().st_size for n in sum_files}
    with open(out / "inputs_SUMS", "w", encoding="utf-8", newline="\n") as w:
        for n in sum_files:
            w.write(f"{sha_map[n]}  {n}\n")
    rep["sha256"] = sha_map
    rep["bytes"] = bytes_map
    (out / "stage_report.json").write_text(json.dumps(rep, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return rep


def _stage_egv1_work(
    *,
    src: Path,
    work: Path,
    obj_dir: Path,
    stage_name: str,
    root: Path | None,
    root_stages: list[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
    source_md5: str | None,
    expected_md5: str | None,
) -> dict:
    skipped_members: list[tuple[str, int, str]] = []
    skipped_records: list[tuple[str, str, str, int]] = []
    trims: list[tuple] = []  # Always empty under DEC-65 (ITRs preserved byte-identically)
    # (genome_name, member, [(contig_uid, length, accession, itr_length)])
    genomes_in: list[tuple[str, str, list[tuple[bytes, int, str, int]]]] = []
    seen_genome_names: set[str] = set()

    with open(work / EGV1_FASTA, "wb") as fa_w:
        for gname, member, raw in _iter_genome_fastas(src, skipped_members):
            if gname in seen_genome_names:
                raise Refused(f"{stage_name}: duplicate genome name {gname!r} ({member})")
            seen_genome_names.add(gname)
            recs: list[tuple[bytes, int, str, int]] = []
            for r in v2.iter_fasta_stream(io.BytesIO(raw)):
                if len(r.seq) == 0:
                    skipped_records.append((gname, r.accession, "empty_sequence", 0))
                    continue
                if not VALID_DNA_RE.match(r.seq):
                    skipped_records.append((gname, r.accession, "non_iupac_dna", len(r.seq)))
                    continue
                seq = r.seq  # Byte-identical, NEVER trimmed (DEC-65)
                itr_len = detect_itr(seq)
                fa_w.write(b"\n".join(fasta.format_record(r.header, seq, r.line_width or 70)) + b"\n")
                cu = ids.contig_uid(seq)
                recs.append((cu, len(seq), r.accession, itr_len))
            if not recs:
                raise Refused(f"{stage_name}: FASTA member {member!r} yielded 0 valid DNA records")
            genomes_in.append((gname, member, recs))

    if not genomes_in:
        raise Refused(f"{stage_name}: no valid genomes found in {src}")

    all_uids = [cu for _, _, recs in genomes_in for cu, _, _, _ in recs]
    owners = compute_prefix_overlap(
        all_uids, root=root, root_stages=root_stages, owner_inputs=owner_inputs
    )

    listing = make_inputs.Listing()
    genomes: list[dict] = []
    rec_rows: list[tuple] = []
    genome_input_rows: list[tuple] = []
    meta_rows: list[tuple] = []
    repeated_rows: list[tuple] = []
    shared_rows: list[tuple] = []
    first_of: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}
    held_whole = 0
    complete_itr_genomes = 0
    ordinal = 0

    for gname, member, recs in genomes_in:
        uids = [cu for cu, _, _, _ in recs]
        for rec_idx, (cu, length, acc, itr_len) in enumerate(recs):
            owner = owners.get(cu, "")
            rec_rows.append((ordinal, gname, rec_idx, acc, cu.hex(), length, itr_len, member, owner))
            if owner:
                shared_rows.append((ordinal, gname, rec_idx, cu.hex(), owner))
            contig_len.setdefault(cu, length)
            ordinal += 1

        distinct_uids = list(dict.fromkeys(uids))
        guid = ids.genome_uid(distinct_uids)
        max_itr = max(itr for _, _, _, itr in recs)
        # Single-contig genomes carrying an exact ITR (>= 20 bp) are complete linear genomes;
        # multi-contig bins or contigs without detected ITRs use schema value "contig".
        quality = "complete" if (len(distinct_uids) == 1 and max_itr > 0) else "contig"
        g_entry = listing.add(uids, 11, "virus", quality)
        if g_entry is None:
            repeated_rows.append((ordinal - len(recs), gname, guid.hex(), first_of[guid]))
            continue
        g_entry["circular"] = []  # Always linear untrimmed (DEC-65)
        genomes.append(g_entry)
        first_of[guid] = gname
        if quality == "complete":
            complete_itr_genomes += 1

        seen_in_g: set[bytes] = set()
        g_held = 0
        g_itr_contigs = 0
        g_lengths: list[int] = []
        for cu, length, _acc, itr_len in recs:
            if cu in seen_in_g:
                continue
            seen_in_g.add(cu)
            owner = owners.get(cu, "")
            if owner:
                g_held += 1
            if itr_len > 0:
                g_itr_contigs += 1
            g_lengths.append(length)
            genome_input_rows.append((
                stage_name, guid.hex(), gname, cu.hex(), "virus", 11,
                "linear", "false", quality, length, itr_len, owner,
            ))

        if all(cu in owners for cu in distinct_uids):
            held_whole += 1

        meta_rows.append((
            stage_name, EGV1_SOURCE["source_name"], gname, guid.hex(),
            len(distinct_uids), sum(g_lengths), min(g_lengths), max(g_lengths),
            "linear", quality, max_itr, g_itr_contigs, g_held, member,
        ))

    meta_rows.sort(key=lambda t: (t[2], t[3]))
    make_inputs._write_tsv(work / "records.tsv", RECORD_HEAD, rec_rows)
    make_inputs._write_tsv(work / "genome_input.tsv", GENOME_INPUT_HEAD, genome_input_rows)
    make_inputs._write_tsv(work / "egv1_metadata.tsv", METADATA_HEAD, meta_rows)
    make_inputs._write_tsv(
        work / "repeated_records.tsv",
        ["ordinal", "genome", "genome_uid", "first_genome"],
        repeated_rows,
    )
    make_inputs._write_tsv(
        work / "shared_records.tsv",
        ["ordinal", "genome", "record", "contig_uid", "owner_stage"],
        shared_rows,
    )
    make_inputs._write_tsv(work / "trims.tsv", TRIM_HEAD, trims)
    make_inputs._write_tsv(work / "skipped_members.tsv", SKIPPED_MEMBER_HEAD, skipped_members)
    make_inputs._write_tsv(work / "skipped_records.tsv", SKIPPED_RECORD_HEAD, skipped_records)

    objects_json = [{
        "path": str(obj_dir / EGV1_FASTA),
        **EGV1_SOURCE,
        "retrieved_at": retrieved,
    }]
    (work / "objects.json").write_text(json.dumps(objects_json, indent=1) + "\n", encoding="utf-8")
    (work / "genomes.json").write_text(json.dumps(genomes) + "\n", encoding="utf-8")

    owner_contigs: Counter[str] = Counter()
    owner_bases: Counter[str] = Counter()
    for cu, length in contig_len.items():
        if cu in owners:
            st = owners[cu]
            owner_contigs[st] += 1
            owner_bases[st] += length

    total_contigs = len(contig_len)
    total_bases = sum(contig_len.values())
    held_contigs = sum(owner_contigs.values())
    held_bases = sum(owner_bases.values())

    return {
        "stage": stage_name,
        "source_name": EGV1_SOURCE["source_name"],
        "license": EGV1_SOURCE["license"],
        "redistributable": True,
        "source_md5": source_md5,
        "expected_md5": expected_md5,
        "fasta_members": len(genomes_in),
        "skipped_members": len(skipped_members),
        "skipped_records": len(skipped_records),
        "records": len(rec_rows),
        "genomes": len(genomes),
        "complete_itr_genomes": complete_itr_genomes,
        "contigs": total_contigs,
        "bases": total_bases,
        "lt_5kb_contigs": sum(1 for n in contig_len.values() if n < 5000),
        "circular_genomes": 0,
        "linear_genomes": len(genomes),
        "trims": 0,
        "repeated_genomes": len(repeated_rows),
        "shared_records": len(shared_rows),
        "genomes_all_contigs_held_earlier": held_whole,
        "owned_by_stage": dict(sorted(owner_contigs.items())),
        "owned_bases_by_stage": dict(sorted(owner_bases.items())),
        "held_contigs": held_contigs,
        "held_bases": held_bases,
        "new_contigs": total_contigs - held_contigs,
        "new_bases": total_bases - held_bases,
        "root_stages": root_stages,
        "owner_inputs": list(owner_inputs),
        **listing.report(),
    }


# --------------------------------------------------------------------------------------------------------------
# Polite public fetcher (DEC-50: User-Agent: curl, <= 2 req/s, stop on 5 consecutive refusals)


class _RateLimiter:
    def __init__(self, min_interval_s: float = 0.55):
        self.min_interval = min_interval_s
        self.last = 0.0
        self.refusals = 0

    def download_file(
        self,
        url: str,
        dst: Path,
        *,
        expected_md5: str | None = None,
        timeout: int = 300,
    ) -> tuple[int, str, str]:
        dst = Path(dst)
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.is_file() and dst.stat().st_size > 0:
            got_md5 = _md5_file(dst)
            if expected_md5 and got_md5 != expected_md5:
                raise Refused(f"{dst.name}: cached md5 {got_md5} != expected {expected_md5}")
            return dst.stat().st_size, ids.sha256_file(dst).hex(), got_md5

        tmp = dst.with_name(f".{dst.name}.part")
        wait = self.min_interval - (time.monotonic() - self.last)
        if wait > 0:
            time.sleep(wait)
        self.last = time.monotonic()
        req = urllib.request.Request(url, headers={"User-Agent": "curl"})
        h_sha = hashlib.sha256()
        h_md5 = hashlib.md5()
        nbytes = 0
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp, open(tmp, "wb") as w:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    w.write(chunk)
                    h_sha.update(chunk)
                    h_md5.update(chunk)
                    nbytes += len(chunk)
            self.refusals = 0
        except urllib.error.HTTPError as e:
            tmp.unlink(missing_ok=True)
            if e.code in (403, 429, 500, 502, 503, 504):
                self.refusals += 1
                if self.refusals >= 5:
                    raise Refused(f"5 consecutive HTTP refusals; last {e.code} on {url}") from e
            raise
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise

        got_md5 = h_md5.hexdigest()
        if expected_md5 and got_md5 != expected_md5:
            tmp.unlink(missing_ok=True)
            raise Refused(f"{dst.name}: downloaded md5 {got_md5} != expected {expected_md5}")
        tmp.rename(dst)
        return nbytes, h_sha.hexdigest(), got_md5


def fetch_egv1(
    dl_dir: Path,
    *,
    url: str = ITEM08_URL,
    filename: str = ITEM08_FILENAME,
    expected_md5: str | None = ITEM08_MD5,
) -> dict:
    """Fetch Figshare Item 08 (`ITEM_08_Egovirales_Genomes.zip`) into `dl_dir` with MD5 + SHA-256 verification."""
    dl_dir = Path(dl_dir)
    dl_dir.mkdir(parents=True, exist_ok=True)
    dst = dl_dir / filename
    rl = _RateLimiter(0.55)
    nbytes, sha_hex, md5_hex = rl.download_file(url, dst, expected_md5=expected_md5)
    rep = {
        "article_id": FIGSHARE_ARTICLE_ID,
        "doi": FIGSHARE_DOI,
        "file_id": ITEM08_FILE_ID,
        "filename": filename,
        "url": url,
        "path": str(dst),
        "bytes": nbytes,
        "expected_bytes": ITEM08_BYTES,
        "md5": md5_hex,
        "expected_md5": expected_md5,
        "sha256": sha_hex,
        "license": EGV1_SOURCE["license"],
        "redistributable": True,
    }
    (dl_dir / "fetch_SUMS").write_text(f"{sha_hex}  {filename}\n", encoding="utf-8")
    (dl_dir / "fetch_report.json").write_text(json.dumps(rep, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return rep


# --------------------------------------------------------------------------------------------------------------
# CLI


def _parse_kv(items: list[str]) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for item in items:
        if "=" not in item:
            raise Refused(f"--owner-inputs expects NAME=DIR, got {item!r}")
        k, v = item.split("=", 1)
        if not k or not v:
            raise Refused(f"--owner-inputs expects NAME=DIR, got {item!r}")
        out[k] = Path(v)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Card C-69: fetch and stage EGV1 (Egovirales gut giant-virus genomes)."
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    fp = sub.add_parser("fetch", help="Download Figshare Item 08 (ITEM_08_Egovirales_Genomes.zip)")
    fp.add_argument("--dl", type=Path, required=True, help="Download directory")
    fp.add_argument("--url", default=ITEM08_URL)
    fp.add_argument("--expected-md5", default=ITEM08_MD5)
    fp.add_argument("--skip-md5", action="store_true")

    sp = sub.add_parser("stage", help="Stage EGV1 into /tmp/if2453/s2/inputs/EGV1")
    sp.add_argument("--src", type=Path, default=None, help="Path to ITEM_08_Egovirales_Genomes.zip or FASTA dir")
    sp.add_argument("--dl", type=Path, default=None, help="Download dir containing ITEM_08_Egovirales_Genomes.zip")
    sp.add_argument("--out", type=Path, required=True, help="Output stage input directory")
    sp.add_argument("--retrieved", required=True, help="ISO-8601 UTC retrieval timestamp")
    sp.add_argument("--root", type=Path, default=None, help="Read-only bio-box data root for DEC-26 overlap")
    sp.add_argument(
        "--root-stages",
        default="S1,P1,P2,P3,P4,P5,P6,H1,PL1,C5,V1,G1,SPX1",
        help="Comma-separated stages already stored in --root",
    )
    sp.add_argument(
        "--owner-inputs",
        action="append",
        default=[],
        help="Staged earlier inputs as STAGE=DIR (repeatable, in build order)",
    )
    sp.add_argument("--final-dir", type=Path, default=None, help="Final object path prefix in objects.json")
    sp.add_argument("--expected-md5", default=ITEM08_MD5)
    sp.add_argument("--skip-md5", action="store_true")
    sp.add_argument("--out-report", type=Path, default=None)

    args = ap.parse_args(argv)
    try:
        if args.cmd == "fetch":
            rep = fetch_egv1(
                args.dl,
                url=args.url,
                expected_md5=None if args.skip_md5 else args.expected_md5,
            )
            print(json.dumps(rep, indent=1, sort_keys=True))
            return 0

        src = args.src
        if src is None:
            if args.dl is None:
                raise Refused("stage requires either --src or --dl")
            src = args.dl / ITEM08_FILENAME
        root_stages = [s.strip() for s in args.root_stages.split(",") if s.strip()]
        owner_inputs = _parse_kv(args.owner_inputs)
        exp_md5 = None if (args.skip_md5 or not src.is_file()) else args.expected_md5
        rep = stage_egv1(
            src,
            args.out,
            root=args.root,
            root_stages=root_stages,
            owner_inputs=owner_inputs,
            retrieved=args.retrieved,
            final_dir=args.final_dir,
            expected_md5=exp_md5,
        )
        if args.out_report is not None:
            args.out_report.write_text(json.dumps(rep, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(rep, indent=1, sort_keys=True))
        return 0
    except Refused as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
