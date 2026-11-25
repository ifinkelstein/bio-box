"""Card C-70 (DEC-59, DEC-65, DEC-66, DEC-69, DEC-50, DEC-51, DEC-26, DEC-27):
fetch, deduplicate, and stage Tier-1 compact high-yield MGE, divergent dsDNA virus,
and post-GlobDB Asgardarchaeota datasets (< 2 Gbp total) for r26.11.5-r26.11.6:

  1. `ECE2` (kind="plasmid", genetic_code=11):
     - `Inocles` giant circular oral *Streptococcus* ECEs (Kiguchi et al. 2025 *Nat. Comms.* 16:7397,
       `PRJDB18431`, `zenodo.16656823` closed circular contigs [CC BY 4.0]; `zenodo.13138926` holds `InoC` marker).
       Preserves circular topology (`circular=[cu]`) untrimmed (`DEC-65`).
     - `Tycheposons` marine tyrosine-recombinase vesicle-forming integrative elements (Hackl & Chisholm 2023 *Cell*,
       GitHub `thackl/pro-tycheposons` [MIT / CC BY-SA 4.0], `PRJNA855972`).
       Always staged linear untrimmed (`circular=[]`, `DEC-65`).

  2. `MIR1` (kind="virus", genetic_code=11):
     - *Mirusviricota* global expansion (`1,257` non-redundant genomes, Medvedeva, Delmont, Krupovic, Koonin,
       Gaia et al. 2025/2026 *Nat. Microbiol.*, Figshare `10.6084/m9.figshare.28955240` [CC BY 4.0]) +
       Gaia et al. 2023 `GOEV` *Mirusviricota* (Figshare `10.6084/m9.figshare.20284713` [CC BY 4.0]).
     - Polinton-like viruses (`PLVs`) & protist EVEs (`543`+ aquatic PLV/virophage genomes, Bellas & Sommaruga 2021
       *Microbiome*, Figshare `9a7a1d16d77ea9d658a1`; Bellas, Hackl et al. 2023 *PNAS*, Figshare
       `10.6084/m9.figshare.21581355`).
       Always staged linear untrimmed (`circular=[]`, `trims=0`, `DEC-65`).

  3. `ASG1` (kind="prokaryote", genetic_code=11):
     - Appler, De Anda & Baker (Feb 2026 *Nature* 652:405-415) `404` post-`GlobDB r226` *Asgardarchaeota* MAGs
       (Figshare `10.6084/m9.figshare.26426707` [CC BY 4.0] + NCBI).
     - Deduplicated at `skani >= 96%` ANI (`AF >= 50%`) against `G1`/GTDB archaea and within-stage, plus
       exact `contig_uid` ownership tracking (`DEC-26`, `DEC-27`).

Subcommands:
    fetch-ece2  --dl DIR
    stage-ece2  --dl DIR --out DIR --retrieved ISO [--root ROOT] [--root-stages S,...] [--owner-inputs K=V ...]
    fetch-mir1  --dl DIR
    stage-mir1  --dl DIR --out DIR --retrieved ISO [--root ROOT] [--root-stages S,...] [--owner-inputs K=V ...]
    fetch-asg1  --dl DIR [--accessions FILE]
    stage-asg1  --dl DIR --out DIR --retrieved ISO [--root ROOT] [--root-stages S,...] [--owner-inputs K=V ...]
                [--skani-ref-tsv TSV] [--skani-ref-list FILE] [--skani-self-dedup] [--skani-ani 96.0]
                [--skani-af 50.0] [--skani-threads 8]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import Counter, defaultdict
from collections.abc import Iterable
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import c48c_virus_inputs as v2  # noqa: E402
import make_inputs  # noqa: E402
from biobox import fasta, ids  # noqa: E402

Refused = make_inputs.Refused

# --------------------------------------------------------------------------------------------------------------
# Constants & verified source endpoints

ECE2_CLASSES = ("inocle", "tycheposon")
MIR1_CLASSES = ("mirusviricota", "mirus_goev", "plv", "plv_eve")
FASTA_SUFFIXES = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")
VALID_DNA_RE = re.compile(rb"^[ACGTURYKMSWBDHVNacgturykmswbdhvn]+$")
_COMP_TABLE = bytes.maketrans(
    b"ACGTacgtRYKMSWBDHVNrykmswbdhvn",
    b"TGCAtgcaYRMKWSVHDBNyrmkwsvhdbn",
)

ZENODO_INOCLES_GENOMES_ID = 16656823  # Closed circular Inocle genomic contigs (CC BY 4.0)
ZENODO_INOCLES_INOC_ID = 13138926     # InoC marker gene deposit
PRJDB18431_ACC = "PRJDB18431"
TYCHEPOSON_REPO = "thackl/pro-tycheposons"
TYCHEPOSON_DIRS = ("elements-pro", "elements-gorg-single-cells", "elements-tara-viral")
PRJNA855972_ACC = "PRJNA855972"

FIGSHARE_MIRUS_2025_ID = 28955240     # 1,257 Mirusviricota genomes (CC BY 4.0)
FIGSHARE_MIRUS_GOEV_ID = 20284713     # Gaia 2023 GOEV / Mirusviricota (CC BY 4.0)
FIGSHARE_PLV_2021_SHARE = "9a7a1d16d77ea9d658a1"  # Bellas & Sommaruga 2021 (543 PLVs)
FIGSHARE_PLV_2023_ID = 21581355       # Bellas, Hackl et al. 2023 PNAS (PLVs & EVEs)
FIGSHARE_ASGARD_2026_ID = 26426707    # Appler, De Anda & Baker 2026 Nature (CC BY 4.0)

DEFAULT_ROOT_STAGES = "S1,P1,P2,P3,P4,P5,P6,H1,PL1,C5,V1,G1,SPX1"

REDISTRIBUTABLE_LICENSE_RE = re.compile(
    r"(?:^|\b)(?:cc[\s-]*by(?:[\s-]*sa)?[\s-]*4\.0|cc[\s-]*by(?:[\s-]*sa)?[\s-]*3\.0|cc0|mit|bsd|apache|insdc|public[\s-]*domain)(?:\b|$)",
    re.IGNORECASE,
)
NON_REDISTRIBUTABLE_LICENSE_RE = re.compile(r"(?:cc[\s-]*by[\s-]*nc|non-?commercial|restricted|all rights reserved)", re.IGNORECASE)

ECE2_RECORD_HEAD = [
    "ordinal", "source", "element_class", "accession", "contig_uid", "length",
    "circular", "tir_length", "owner_stage",
]
ECE2_SIDECAR_HEAD = [
    "stage", "source", "element_class", "genome_uid", "contig_uid", "accession",
    "length", "circular", "tir_length", "owner_stage",
]
MIR1_RECORD_HEAD = [
    "ordinal", "source", "virus_class", "genome_id", "accession", "contig_uid",
    "length", "circular", "tir_length", "owner_stage",
]
MIR1_SIDECAR_HEAD = [
    "stage", "source", "virus_class", "genome_id", "genome_uid", "contigs",
    "bases", "circular", "max_tir_length", "held_contigs", "owner_stages",
]
ASG1_RECORD_HEAD = [
    "ordinal", "source", "genome_id", "record_ordinal", "accession",
    "contig_uid", "length", "assembly_quality", "owner_stage",
]
ASG1_SIDECAR_HEAD = [
    "stage", "source", "genome_id", "genome_uid", "contigs", "bases",
    "assembly_quality", "genetic_code", "circular_contigs", "held_contigs", "owner_stages",
]
TRIM_HEAD = [
    "source", "genome", "accession", "source_contig_uid", "contig_uid",
    "source_length", "repeat_length", "bases_removed", "length",
]
SKIPPED_HEAD = ["source", "accession", "reason", "length"]
SKANI_DROPPED_HEAD = ["genome_id", "source_name", "reason", "ref_match", "ani", "af"]


# --------------------------------------------------------------------------------------------------------------
# Pure helpers: license checking, TIR detection, classification, FASTA I/O


def is_redistributable_license(license_str: str) -> bool:
    """Return True iff `license_str` clearly names a redistributable open license (DEC-51)."""
    s = (license_str or "").strip()
    if not s:
        return False
    if NON_REDISTRIBUTABLE_LICENSE_RE.search(s):
        return False
    return bool(REDISTRIBUTABLE_LICENSE_RE.search(s))


def detect_tir(seq: bytes, min_len: int = 20, max_len: int = 100_000) -> int:
    """Detect exact terminal inverted repeat (TIR) length at the ends of a linear sequence.

    Read-only: NEVER trims or modifies `seq` (`DEC-65`).
    """
    n = len(seq)
    if n < 2 * min_len:
        return 0
    limit = min(max_len, n // 2)
    upper = seq[:limit].upper()
    tail_rc = seq[n - limit :].upper().translate(_COMP_TABLE)[::-1]
    k = 0
    while k < limit and upper[k] == tail_rc[k]:
        k += 1
    return k if k >= min_len else 0


def _is_fasta_name(name: str) -> bool:
    low = name.lower()
    if "/__macosx/" in low or low.startswith("__macosx/") or Path(name).name.startswith("._"):
        return False
    return low.endswith(FASTA_SUFFIXES)


def _strip_fasta_suffix(name: str) -> str:
    low = name.lower()
    for suf in sorted(FASTA_SUFFIXES, key=len, reverse=True):
        if low.endswith(suf):
            return name[: -len(suf)]
    return Path(name).stem


def _safe_float(val: str | float | int | None, default: float = 0.0) -> float:
    try:
        return float(val) if val is not None and val != "" else default
    except (TypeError, ValueError):
        return default


def _open_fasta_in(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return open(path, "rb")


class _DetGzipWriter:
    """Deterministic gzip writer (`mtime=0`, `filename=""`) or plain binary file writer."""

    def __init__(self, path: Path):
        self._path = path
        self._raw = None
        self._gz = None

    def __enter__(self):
        self._raw = open(self._path, "wb")
        if self._path.suffix == ".gz":
            self._gz = gzip.GzipFile(filename="", mode="wb", fileobj=self._raw, compresslevel=6, mtime=0)
            return self._gz
        return self._raw

    def __exit__(self, exc_type, exc, tb):
        if self._gz is not None:
            self._gz.close()
        if self._raw is not None:
            self._raw.close()


def classify_ece2_record(header: str, element_class: str, topology: str | None = None) -> dict:
    """Classify an ECE2 record (`inocle` or `tycheposon`) with strict DEC-65 topology rules.

    - `inocle`: closed circular giant oral *Streptococcus* ECE -> `circular=True` unless
      explicitly marked linear/partial/fragment in `topology` or `header`.
    - `tycheposon`: integrative element -> ALWAYS `circular=False` (linear untrimmed).
    """
    if element_class not in ECE2_CLASSES:
        raise Refused(f"ECE2: unknown element_class {element_class!r} (expected one of {ECE2_CLASSES})")
    tok = header.split(None, 1)[0] if header else ""
    low_h = header.lower()
    topo = (topology or "").strip().lower()

    if element_class == "tycheposon":
        circular = False
    else:
        if topo in ("linear", "fragment", "partial"):
            circular = False
        elif topo == "circular":
            circular = True
        elif any(w in low_h for w in ("linear", "partial", "fragment", "incomplete")):
            circular = False
        else:
            circular = True

    return {
        "accession": tok,
        "element_class": element_class,
        "circular": circular,
        "kind": "plasmid",
    }


def classify_mir1_record(header: str, virus_class: str) -> dict:
    """Classify a MIR1 record (`mirusviricota`, `mirus_goev`, `plv`, `plv_eve`).

    All MIR1 records are staged `kind="virus"`, `circular=False` (linear untrimmed per DEC-65).
    """
    if virus_class not in MIR1_CLASSES:
        raise Refused(f"MIR1: unknown virus_class {virus_class!r} (expected one of {MIR1_CLASSES})")
    tok = header.split(None, 1)[0] if header else ""
    return {
        "accession": tok,
        "virus_class": virus_class,
        "circular": False,
        "kind": "virus",
    }


def _finalize_stage_dir(work: Path, out_dir: Path, rep: dict) -> dict:
    """Atomically rename `work` -> `out_dir`, write `inputs_SUMS` and `stage_report.json`."""
    work.rename(out_dir)
    sum_files = sorted(
        p.name for p in out_dir.iterdir()
        if p.is_file() and p.name not in ("inputs_SUMS", "stage_report.json")
    )
    sha_map = {n: ids.sha256_file(out_dir / n).hex() for n in sum_files}
    bytes_map = {n: (out_dir / n).stat().st_size for n in sum_files}
    with open(out_dir / "inputs_SUMS", "w", encoding="utf-8") as w:
        for n in sum_files:
            w.write(f"{sha_map[n]}  {n}\n")
    rep["sha256"] = sha_map
    rep["bytes"] = bytes_map
    (out_dir / "stage_report.json").write_text(json.dumps(rep, indent=1, sort_keys=True), encoding="utf-8")
    return rep


# --------------------------------------------------------------------------------------------------------------
# Stage 1: ECE2 (Inocles + Tycheposons)


def stage_ece2(
    sources: list[dict],
    out_dir: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
    final_dir: Path | None = None,
    stage_name: str = "ECE2",
) -> dict:
    """Stage `ECE2` (Inocles circular + Tycheposons linear, 0 trims per DEC-65, kind='plasmid')."""
    if not sources:
        raise Refused(f"{stage_name}: no sources to stage")
    out_dir = Path(out_dir)
    if out_dir.exists():
        raise Refused(f"{out_dir} exists")
    for s in sources:
        if not s.get("redistributable", False) or not is_redistributable_license(s.get("license", "")):
            raise Refused(
                f"{stage_name}: refusing non-redistributable source {s.get('source_name')!r} "
                f"(license={s.get('license')!r}, DEC-51)"
            )
        if s.get("element_class") not in ECE2_CLASSES:
            raise Refused(f"{stage_name}: invalid element_class {s.get('element_class')!r}")
        fa = Path(s["fasta"])
        if not fa.is_file():
            raise Refused(f"{stage_name}: FASTA file not found: {fa}")

    work = out_dir.with_name(f".{out_dir.name}.tmp")
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    obj_dir = Path(final_dir) if final_dir is not None else out_dir

    try:
        rep = _stage_ece2_work(
            sources, work, obj_dir, stage_name, root, list(root_stages), owner_inputs, retrieved
        )
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise

    return _finalize_stage_dir(work, out_dir, rep)


def _stage_ece2_work(
    sources: list[dict],
    work: Path,
    obj_dir: Path,
    stage_name: str,
    root: Path | None,
    root_stages: list[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
) -> dict:
    ordered_sources = sorted(sources, key=lambda s: (s["element_class"], s["output_name"], s["source_name"]))
    staged_records: list[tuple[str, str, str, bytes, int, bool, int]] = []
    skipped: list[tuple[str, str, str, int]] = []
    trims: list[tuple] = []  # Always empty under DEC-65
    objects_json: list[dict] = []
    by_source: dict[str, dict] = {}

    for src in ordered_sources:
        sname = src["source_name"]
        out_name = src["output_name"]
        fa_path = Path(src["fasta"])
        elem_cls = src["element_class"]
        topo_map: dict[str, str] = src.get("topology_map") or {}
        dst_fa = work / out_name
        rec_count = 0
        base_count = 0
        skip_count = 0

        with _open_fasta_in(fa_path) as fh, _DetGzipWriter(dst_fa) as w:
            for r in v2.iter_fasta_stream(fh):
                if len(r.seq) == 0:
                    skipped.append((sname, r.accession, "empty_sequence", 0))
                    skip_count += 1
                    continue
                if not VALID_DNA_RE.match(r.seq):
                    skipped.append((sname, r.accession, "non_iupac_dna", len(r.seq)))
                    skip_count += 1
                    continue
                meta = classify_ece2_record(r.header, elem_cls, topo_map.get(r.accession))
                seq = r.seq  # Byte-identical, NEVER trimmed (DEC-65)
                circ = bool(meta["circular"])
                tir_len = 0 if circ else detect_tir(seq)
                w.write(b"\n".join(fasta.format_record(r.header, seq, r.line_width)) + b"\n")
                cu = ids.contig_uid(seq)
                staged_records.append((sname, elem_cls, meta["accession"], cu, len(seq), circ, tir_len))
                rec_count += 1
                base_count += len(seq)

        if rec_count == 0:
            dst_fa.unlink(missing_ok=True)
            continue

        objects_json.append({
            "path": str(obj_dir / out_name),
            "source_name": sname,
            "source_release": src["source_release"],
            "license": src["license"],
            "redistributable": True,
            "retrieved_at": retrieved,
        })
        by_source[out_name] = {
            "source_name": sname,
            "element_class": elem_cls,
            "records": rec_count,
            "bases": base_count,
            "skipped": skip_count,
        }

    if not staged_records:
        raise Refused(f"{stage_name}: no valid records staged")

    all_uids = [x[3] for x in staged_records]
    owners = v2.owners_of(all_uids, root, root_stages, owner_inputs)

    listing = make_inputs.Listing()
    genomes: list[dict] = []
    rec_rows: list[tuple] = []
    sidecar_rows: list[tuple] = []
    repeated_rows: list[tuple] = []
    shared_rows: list[tuple] = []
    first_of: dict[bytes, str] = {}
    class_counts: Counter[str] = Counter()
    genome_class_counts: Counter[str] = Counter()
    circular_genomes = 0
    tir_count = 0

    for ordinal, (sname, elem_cls, acc, cu, length, circ, tir_len) in enumerate(staged_records):
        owner = owners.get(cu, "")
        rec_rows.append((ordinal, sname, elem_cls, acc, cu.hex(), length, str(circ).lower(), tir_len, owner))
        if owner:
            shared_rows.append((ordinal, acc, cu.hex(), owner))
        class_counts[elem_cls] += 1

        qual = "complete" if circ else "contig"
        g = listing.add([cu], 11, "plasmid", qual)
        guid = ids.genome_uid([cu])
        if g is None:
            repeated_rows.append((ordinal, acc, guid.hex(), first_of[guid]))
            continue
        g["circular"] = [cu.hex()] if circ else []
        genomes.append(g)
        first_of[guid] = acc
        genome_class_counts[elem_cls] += 1
        if circ:
            circular_genomes += 1
        if tir_len > 0:
            tir_count += 1
        sidecar_rows.append((
            stage_name, sname, elem_cls, guid.hex(), cu.hex(), acc,
            length, str(circ).lower(), tir_len, owner,
        ))

    sidecar_rows.sort(key=lambda t: (t[2], t[1], t[5], t[4]))
    make_inputs._write_tsv(work / "records.tsv", ECE2_RECORD_HEAD, rec_rows)
    make_inputs._write_tsv(work / "ece2_metadata.tsv", ECE2_SIDECAR_HEAD, sidecar_rows)
    make_inputs._write_tsv(
        work / "repeated_records.tsv",
        ["ordinal", "accession", "genome_uid", "first_accession"],
        repeated_rows,
    )
    make_inputs._write_tsv(
        work / "shared_records.tsv",
        ["ordinal", "accession", "contig_uid", "owner_stage"],
        shared_rows,
    )
    make_inputs._write_tsv(work / "trims.tsv", TRIM_HEAD, trims)
    make_inputs._write_tsv(work / "skipped_records.tsv", SKIPPED_HEAD, skipped)
    (work / "objects.json").write_text(json.dumps(objects_json, indent=1), encoding="utf-8")
    (work / "genomes.json").write_text(json.dumps(genomes), encoding="utf-8")

    contig_len: dict[bytes, int] = {}
    contig_cls: dict[bytes, str] = {}
    for _, elem_cls, _, cu, length, *_ in staged_records:
        contig_len.setdefault(cu, length)
        contig_cls.setdefault(cu, elem_cls)

    owner_contigs: Counter[str] = Counter()
    owner_bases: Counter[str] = Counter()
    bases_by_class: Counter[str] = Counter()
    for cu, length in contig_len.items():
        bases_by_class[contig_cls[cu]] += length
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
        "redistributable": True,
        "records": len(staged_records),
        "skipped_records": len(skipped),
        "genomes": len(genomes),
        "contigs": total_contigs,
        "bases": total_bases,
        "circular_genomes": circular_genomes,
        "linear_genomes": len(genomes) - circular_genomes,
        "genomes_with_tir": tir_count,
        "repeated_records": len(repeated_rows),
        "trims": 0,
        "owned_by_stage": dict(sorted(owner_contigs.items())),
        "owned_bases_by_stage": dict(sorted(owner_bases.items())),
        "held_contigs": held_contigs,
        "held_bases": held_bases,
        "new_contigs": total_contigs - held_contigs,
        "new_bases": total_bases - held_bases,
        "records_by_class": dict(sorted(class_counts.items())),
        "genomes_by_class": dict(sorted(genome_class_counts.items())),
        "bases_by_class": dict(sorted(bases_by_class.items())),
        "by_source": by_source,
        "root_stages": root_stages,
        "owner_inputs": sorted(owner_inputs),
    }


# --------------------------------------------------------------------------------------------------------------
# Stage 2: MIR1 (Mirusviricota + GOEV + PLVs & EVEs)


def stage_mir1(
    sources: list[dict],
    out_dir: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
    final_dir: Path | None = None,
    stage_name: str = "MIR1",
) -> dict:
    """Stage `MIR1` (Mirusviricota + GOEV + PLVs/EVEs; linear untrimmed per DEC-65, kind='virus').

    Each source dict supports two modes:
      - `mode == "single_contig"` (default): `fasta` is one FASTA file; each record is one genome.
      - `mode == "multi_contig_bins"`: `bin_fastas` is a list of `(genome_id, fasta_path)` pairs
        consolidated deterministically into `output_name`, one genome per bin.
    """
    if not sources:
        raise Refused(f"{stage_name}: no sources to stage")
    out_dir = Path(out_dir)
    if out_dir.exists():
        raise Refused(f"{out_dir} exists")
    for s in sources:
        if not s.get("redistributable", False) or not is_redistributable_license(s.get("license", "")):
            raise Refused(
                f"{stage_name}: refusing non-redistributable source {s.get('source_name')!r} "
                f"(license={s.get('license')!r}, DEC-51)"
            )
        if s.get("virus_class") not in MIR1_CLASSES:
            raise Refused(f"{stage_name}: invalid virus_class {s.get('virus_class')!r}")
        mode = s.get("mode", "single_contig")
        if mode == "single_contig":
            if not Path(s["fasta"]).is_file():
                raise Refused(f"{stage_name}: FASTA file not found: {s['fasta']}")
        elif mode == "multi_contig_bins":
            if not s.get("bin_fastas"):
                raise Refused(f"{stage_name}: empty bin_fastas for {s.get('source_name')!r}")
            for _, bp in s["bin_fastas"]:
                if not Path(bp).is_file():
                    raise Refused(f"{stage_name}: bin FASTA not found: {bp}")
        else:
            raise Refused(f"{stage_name}: unknown mode {mode!r}")

    work = out_dir.with_name(f".{out_dir.name}.tmp")
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    obj_dir = Path(final_dir) if final_dir is not None else out_dir

    try:
        rep = _stage_mir1_work(
            sources, work, obj_dir, stage_name, root, list(root_stages), owner_inputs, retrieved
        )
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise

    return _finalize_stage_dir(work, out_dir, rep)


def _stage_mir1_work(
    sources: list[dict],
    work: Path,
    obj_dir: Path,
    stage_name: str,
    root: Path | None,
    root_stages: list[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
) -> dict:
    ordered_sources = sorted(sources, key=lambda s: (s["virus_class"], s["output_name"], s["source_name"]))
    # Each bin entry: (source_name, virus_class, genome_id, quality, list[(acc, cu, length, tir_len)])
    genome_bins: list[tuple[str, str, str, str, list[tuple[str, bytes, int, int]]]] = []
    skipped: list[tuple[str, str, str, int]] = []
    trims: list[tuple] = []  # Always empty under DEC-65
    objects_json: list[dict] = []
    by_source: dict[str, dict] = {}

    for src in ordered_sources:
        sname = src["source_name"]
        out_name = src["output_name"]
        vcls = src["virus_class"]
        mode = src.get("mode", "single_contig")
        qual = src.get("assembly_quality", "mag_mq" if mode == "multi_contig_bins" else "contig")
        dst_fa = work / out_name
        rec_count = 0
        base_count = 0
        skip_count = 0
        bin_count = 0

        with _DetGzipWriter(dst_fa) as w:
            if mode == "single_contig":
                with _open_fasta_in(Path(src["fasta"])) as fh:
                    for r in v2.iter_fasta_stream(fh):
                        if len(r.seq) == 0:
                            skipped.append((sname, r.accession, "empty_sequence", 0))
                            skip_count += 1
                            continue
                        if not VALID_DNA_RE.match(r.seq):
                            skipped.append((sname, r.accession, "non_iupac_dna", len(r.seq)))
                            skip_count += 1
                            continue
                        meta = classify_mir1_record(r.header, vcls)
                        seq = r.seq  # Linear untrimmed (DEC-65)
                        tir_len = detect_tir(seq)
                        w.write(b"\n".join(fasta.format_record(r.header, seq, r.line_width)) + b"\n")
                        cu = ids.contig_uid(seq)
                        genome_bins.append((sname, vcls, meta["accession"], qual, [(meta["accession"], cu, len(seq), tir_len)]))
                        rec_count += 1
                        base_count += len(seq)
                        bin_count += 1
            else:
                ordered_bins = sorted(src["bin_fastas"], key=lambda kv: (kv[0], str(kv[1])))
                for gid, bp in ordered_bins:
                    contigs_in_bin: list[tuple[str, bytes, int, int]] = []
                    with _open_fasta_in(Path(bp)) as fh:
                        for r in v2.iter_fasta_stream(fh):
                            if len(r.seq) == 0:
                                skipped.append((sname, f"{gid}:{r.accession}", "empty_sequence", 0))
                                skip_count += 1
                                continue
                            if not VALID_DNA_RE.match(r.seq):
                                skipped.append((sname, f"{gid}:{r.accession}", "non_iupac_dna", len(r.seq)))
                                skip_count += 1
                                continue
                            meta = classify_mir1_record(r.header, vcls)
                            seq = r.seq  # Linear untrimmed (DEC-65)
                            tir_len = detect_tir(seq)
                            w.write(b"\n".join(fasta.format_record(r.header, seq, r.line_width)) + b"\n")
                            cu = ids.contig_uid(seq)
                            contigs_in_bin.append((meta["accession"], cu, len(seq), tir_len))
                            rec_count += 1
                            base_count += len(seq)
                    if contigs_in_bin:
                        genome_bins.append((sname, vcls, gid, qual, contigs_in_bin))
                        bin_count += 1

        if rec_count == 0:
            dst_fa.unlink(missing_ok=True)
            continue

        objects_json.append({
            "path": str(obj_dir / out_name),
            "source_name": sname,
            "source_release": src["source_release"],
            "license": src["license"],
            "redistributable": True,
            "retrieved_at": retrieved,
        })
        by_source[out_name] = {
            "source_name": sname,
            "virus_class": vcls,
            "mode": mode,
            "genomes_input": bin_count,
            "records": rec_count,
            "bases": base_count,
            "skipped": skip_count,
        }

    if not genome_bins:
        raise Refused(f"{stage_name}: no valid records staged")

    all_uids = [cu for _, _, _, _, clist in genome_bins for _, cu, _, _ in clist]
    owners = v2.owners_of(all_uids, root, root_stages, owner_inputs)

    listing = make_inputs.Listing()
    genomes: list[dict] = []
    rec_rows: list[tuple] = []
    sidecar_rows: list[tuple] = []
    repeated_rows: list[tuple] = []
    shared_rows: list[tuple] = []
    first_of: dict[bytes, str] = {}
    class_records: Counter[str] = Counter()
    class_genomes: Counter[str] = Counter()
    contig_len: dict[bytes, int] = {}
    contig_cls: dict[bytes, str] = {}
    ordinal = 0

    for sname, vcls, gid, qual, clist in genome_bins:
        uids: list[bytes] = []
        held_in_bin = 0
        owner_set: set[str] = set()
        bin_bases = 0
        max_tir = 0
        for acc, cu, length, tir_len in clist:
            owner = owners.get(cu, "")
            rec_rows.append((ordinal, sname, vcls, gid, acc, cu.hex(), length, "false", tir_len, owner))
            if owner:
                shared_rows.append((ordinal, gid, acc, cu.hex(), owner))
                held_in_bin += 1
                owner_set.add(owner)
            uids.append(cu)
            bin_bases += length
            if tir_len > max_tir:
                max_tir = tir_len
            contig_len.setdefault(cu, length)
            contig_cls.setdefault(cu, vcls)
            class_records[vcls] += 1
            ordinal += 1

        g = listing.add(uids, 11, "virus", qual)
        guid = ids.genome_uid(list(dict.fromkeys(uids)))
        if g is None:
            repeated_rows.append((gid, guid.hex(), first_of[guid]))
            continue
        g["circular"] = []  # Always linear untrimmed under DEC-65
        genomes.append(g)
        first_of[guid] = gid
        class_genomes[vcls] += 1
        sidecar_rows.append((
            stage_name, sname, vcls, gid, guid.hex(), len(g["contig_uids"]),
            bin_bases, "false", max_tir, held_in_bin, ",".join(sorted(owner_set)),
        ))

    sidecar_rows.sort(key=lambda t: (t[2], t[1], t[3], t[4]))
    make_inputs._write_tsv(work / "records.tsv", MIR1_RECORD_HEAD, rec_rows)
    make_inputs._write_tsv(work / "mir1_metadata.tsv", MIR1_SIDECAR_HEAD, sidecar_rows)
    make_inputs._write_tsv(
        work / "repeated_records.tsv",
        ["genome_id", "genome_uid", "first_genome_id"],
        repeated_rows,
    )
    make_inputs._write_tsv(
        work / "shared_records.tsv",
        ["ordinal", "genome_id", "accession", "contig_uid", "owner_stage"],
        shared_rows,
    )
    make_inputs._write_tsv(work / "trims.tsv", TRIM_HEAD, trims)
    make_inputs._write_tsv(work / "skipped_records.tsv", SKIPPED_HEAD, skipped)
    (work / "objects.json").write_text(json.dumps(objects_json, indent=1), encoding="utf-8")
    (work / "genomes.json").write_text(json.dumps(genomes), encoding="utf-8")

    owner_contigs: Counter[str] = Counter()
    owner_bases: Counter[str] = Counter()
    bases_by_class: Counter[str] = Counter()
    for cu, length in contig_len.items():
        bases_by_class[contig_cls[cu]] += length
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
        "redistributable": True,
        "records": ordinal,
        "skipped_records": len(skipped),
        "genomes": len(genomes),
        "contigs": total_contigs,
        "bases": total_bases,
        "circular_genomes": 0,
        "linear_genomes": len(genomes),
        "repeated_records": len(repeated_rows),
        "trims": 0,
        "owned_by_stage": dict(sorted(owner_contigs.items())),
        "owned_bases_by_stage": dict(sorted(owner_bases.items())),
        "held_contigs": held_contigs,
        "held_bases": held_bases,
        "new_contigs": total_contigs - held_contigs,
        "new_bases": total_bases - held_bases,
        "records_by_class": dict(sorted(class_records.items())),
        "genomes_by_class": dict(sorted(class_genomes.items())),
        "bases_by_class": dict(sorted(bases_by_class.items())),
        "by_source": by_source,
        "root_stages": root_stages,
        "owner_inputs": sorted(owner_inputs),
    }


# --------------------------------------------------------------------------------------------------------------
# Stage 3: ASG1 (Appler, De Anda & Baker 2026 Asgardarchaeota MAGs + skani >= 96% ANI dedup)


def parse_skani_matches(
    tsv_path: Path,
    *,
    min_ani: float = 96.0,
    min_af: float = 50.0,
) -> dict[str, dict]:
    """Parse `skani search` / `skani dist` TSV and return `query_id -> best_hit` with `ANI >= min_ani` & `AF >= min_af`."""
    best: dict[str, dict] = {}
    with open(tsv_path, encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            qfile = (r.get("Query_file") or "").strip()
            rfile = (r.get("Ref_file") or "").strip()
            qname = _strip_fasta_suffix(Path(qfile).name) if qfile else (r.get("Query_name") or "").strip()
            rname = _strip_fasta_suffix(Path(rfile).name) if rfile else (r.get("Ref_name") or "").strip()
            if not qname:
                continue
            ani = _safe_float(r.get("ANI"))
            af_ref = _safe_float(r.get("Align_fraction_ref"))
            af_q = _safe_float(r.get("Align_fraction_query"))
            af = max(af_ref, af_q)
            if ani >= min_ani and af >= min_af:
                prev = best.get(qname)
                if prev is None or (ani, af, rname) > (prev["ani"], prev["af"], prev["ref"]):
                    best[qname] = {
                        "ref": rname,
                        "ani": round(ani, 4),
                        "af": round(af, 4),
                        "af_ref": round(af_ref, 4),
                        "af_query": round(af_q, 4),
                    }
    return best


def filter_asg1_by_skani_ref(
    mag_entries: list[dict],
    ref_matches: dict[str, dict],
    *,
    min_ani: float = 96.0,
) -> tuple[list[dict], list[dict]]:
    """Drop `ASG1` MAGs matching `G1`/GTDB archaea at `>= min_ani` ANI, returning `(kept, dropped)`."""
    kept: list[dict] = []
    dropped: list[dict] = []
    for e in sorted(mag_entries, key=lambda x: x["genome_id"]):
        gid = e["genome_id"]
        hit = ref_matches.get(gid)
        if hit is not None and hit["ani"] >= min_ani:
            dropped.append({
                "genome_id": gid,
                "source_name": e.get("source_name", "appler_2026_asgard"),
                "reason": f"skani_ref_ge_{min_ani}ani",
                "ref_match": hit["ref"],
                "ani": hit["ani"],
                "af": hit["af"],
            })
        else:
            kept.append(e)
    return kept, dropped


def run_skani_dist_against_ref(
    mag_entries: list[dict],
    ref_list_path: Path,
    work_dir: Path,
    *,
    min_ani: float = 96.0,
    min_af: float = 50.0,
    threads: int = 8,
) -> dict[str, dict]:
    """Run `skani dist` of `ASG1` MAGs against a reference file list (e.g. G1/GTDB archaea)."""
    if not mag_entries:
        return {}
    work_dir.mkdir(parents=True, exist_ok=True)
    qlist = work_dir / "asg1_queries.txt"
    with open(qlist, "w", encoding="utf-8") as fh:
        for e in sorted(mag_entries, key=lambda x: x["genome_id"]):
            fh.write(str(Path(e["fasta"]).resolve()) + "\n")
    out_tsv = work_dir / "asg1_vs_ref_skani.tsv"
    cmd = [
        "skani", "dist",
        "--rl", str(ref_list_path),
        "--ql", str(qlist),
        "-t", str(threads),
        "--min-af", str(min_af),
        "-o", str(out_tsv),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise Refused(f"skani dist failed (rc={proc.returncode}): {proc.stderr[:500]}")
    return parse_skani_matches(out_tsv, min_ani=min_ani, min_af=min_af)


def run_skani_self_dedup(
    mag_entries: list[dict],
    work_dir: Path,
    *,
    min_ani: float = 96.0,
    min_af: float = 50.0,
    threads: int = 8,
) -> tuple[list[dict], list[dict]]:
    """Run `skani triangle` on `ASG1` MAGs to dereplicate within-stage `>= 96%` ANI / `50%` AF."""
    if len(mag_entries) <= 1:
        return list(mag_entries), []
    work_dir.mkdir(parents=True, exist_ok=True)
    ordered = sorted(
        mag_entries,
        key=lambda e: (
            0 if e.get("assembly_quality") == "mag_hq" else 1,
            -int(e.get("size_bytes", Path(e["fasta"]).stat().st_size)),
            e["genome_id"],
        ),
    )
    flist = work_dir / "asg1_self_files.txt"
    path_to_entry: dict[str, dict] = {}
    with open(flist, "w", encoding="utf-8") as fh:
        for e in ordered:
            p = str(Path(e["fasta"]).resolve())
            path_to_entry[p] = e
            fh.write(p + "\n")

    out_tsv = work_dir / "asg1_triangle.tsv"
    cmd = [
        "skani", "triangle",
        "-l", str(flist),
        "-t", str(threads),
        "--sparse",
        "--min-af", str(min_af),
        "-o", str(out_tsv),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise Refused(f"skani triangle failed (rc={proc.returncode}): {proc.stderr[:500]}")

    rank = {str(Path(e["fasta"]).resolve()): idx for idx, e in enumerate(ordered)}
    adj: dict[str, dict[str, tuple[float, float]]] = defaultdict(dict)
    with open(out_tsv, encoding="utf-8", errors="replace", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rf = str(Path(r.get("Ref_file", "")).resolve())
            qf = str(Path(r.get("Query_file", "")).resolve())
            if rf == qf or rf not in rank or qf not in rank:
                continue
            ani = _safe_float(r.get("ANI"))
            af = max(_safe_float(r.get("Align_fraction_ref")), _safe_float(r.get("Align_fraction_query")))
            if ani >= min_ani and af >= min_af:
                adj[rf][qf] = (round(ani, 4), round(af, 4))
                adj[qf][rf] = (round(ani, 4), round(af, 4))

    kept: list[dict] = []
    dropped: list[dict] = []
    eliminated: set[str] = set()
    for e in ordered:
        p = str(Path(e["fasta"]).resolve())
        if p in eliminated:
            continue
        kept.append(e)
        for nb, (ani, af) in sorted(adj.get(p, {}).items(), key=lambda kv: rank[kv[0]]):
            if nb not in eliminated and rank[nb] > rank[p]:
                eliminated.add(nb)
                dropped.append({
                    "genome_id": path_to_entry[nb]["genome_id"],
                    "source_name": path_to_entry[nb].get("source_name", "appler_2026_asgard"),
                    "reason": f"skani_within_stage_ge_{min_ani}ani",
                    "ref_match": e["genome_id"],
                    "ani": ani,
                    "af": af,
                })
    return kept, dropped


def stage_asg1(
    mag_entries: list[dict],
    out_dir: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
    final_dir: Path | None = None,
    stage_name: str = "ASG1",
    skani_ref_matches: dict[str, dict] | None = None,
    skani_min_ani: float = 96.0,
    skani_extra_dropped: list[dict] | None = None,
) -> dict:
    """Stage `ASG1` (Appler 2026 *Asgardarchaeota* MAGs; kind='prokaryote', code=11, skani >= 96% ANI dedup)."""
    if not mag_entries:
        raise Refused(f"{stage_name}: no MAG entries to stage")
    out_dir = Path(out_dir)
    if out_dir.exists():
        raise Refused(f"{out_dir} exists")
    for e in mag_entries:
        if not e.get("redistributable", False) or not is_redistributable_license(e.get("license", "")):
            raise Refused(
                f"{stage_name}: refusing non-redistributable MAG {e.get('genome_id')!r} "
                f"(license={e.get('license')!r}, DEC-51)"
            )
        if not Path(e["fasta"]).is_file():
            raise Refused(f"{stage_name}: MAG FASTA not found: {e['fasta']}")

    skani_dropped: list[dict] = list(skani_extra_dropped or [])
    if skani_ref_matches:
        mag_entries, ref_dropped = filter_asg1_by_skani_ref(mag_entries, skani_ref_matches, min_ani=skani_min_ani)
        skani_dropped.extend(ref_dropped)
    if not mag_entries:
        raise Refused(f"{stage_name}: all MAGs dropped by skani >= {skani_min_ani}% ANI filter")

    work = out_dir.with_name(f".{out_dir.name}.tmp")
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    obj_dir = Path(final_dir) if final_dir is not None else out_dir

    try:
        rep = _stage_asg1_work(
            mag_entries, work, obj_dir, stage_name, root, list(root_stages),
            owner_inputs, retrieved, skani_dropped,
        )
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise

    return _finalize_stage_dir(work, out_dir, rep)


def _stage_asg1_work(
    mag_entries: list[dict],
    work: Path,
    obj_dir: Path,
    stage_name: str,
    root: Path | None,
    root_stages: list[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
    skani_dropped: list[dict],
) -> dict:
    ordered = sorted(mag_entries, key=lambda e: (e["genome_id"], str(e["fasta"])))
    fa_name = "asg1_asgard_mags.fna.gz"
    dst_fa = work / fa_name
    skipped: list[tuple[str, str, str, int]] = []
    trims: list[tuple] = []  # Always empty under DEC-65
    # (genome_id, source_name, quality, list[(acc, cu, length)])
    staged_mags: list[tuple[str, str, str, list[tuple[str, bytes, int]]]] = []

    with _DetGzipWriter(dst_fa) as w:
        for e in ordered:
            gid = e["genome_id"]
            sname = e.get("source_name", "appler_2026_asgard")
            qual = e.get("assembly_quality", "mag_mq")
            if qual not in make_inputs.PROK_QUALITIES:
                raise Refused(f"{stage_name}: invalid assembly_quality {qual!r} for {gid}")
            contigs: list[tuple[str, bytes, int]] = []
            with _open_fasta_in(Path(e["fasta"])) as fh:
                for r in v2.iter_fasta_stream(fh):
                    if len(r.seq) == 0:
                        skipped.append((sname, f"{gid}:{r.accession}", "empty_sequence", 0))
                        continue
                    if not VALID_DNA_RE.match(r.seq):
                        skipped.append((sname, f"{gid}:{r.accession}", "non_iupac_dna", len(r.seq)))
                        continue
                    w.write(b"\n".join(fasta.format_record(r.header, r.seq, r.line_width)) + b"\n")
                    cu = ids.contig_uid(r.seq)
                    contigs.append((r.accession, cu, len(r.seq)))
            if contigs:
                staged_mags.append((gid, sname, qual, contigs))

    if not staged_mags:
        raise Refused(f"{stage_name}: no valid MAG contigs staged")

    all_uids = [cu for _, _, _, clist in staged_mags for _, cu, _ in clist]
    owners = v2.owners_of(all_uids, root, root_stages, owner_inputs)

    listing = make_inputs.Listing()
    genomes: list[dict] = []
    rec_rows: list[tuple] = []
    sidecar_rows: list[tuple] = []
    repeated_rows: list[tuple] = []
    shared_rows: list[tuple] = []
    first_of: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}
    ordinal = 0

    for gid, sname, qual, clist in staged_mags:
        uids: list[bytes] = []
        held_in_mag = 0
        owner_set: set[str] = set()
        mag_bases = 0
        for rec_ord, (acc, cu, length) in enumerate(clist):
            owner = owners.get(cu, "")
            rec_rows.append((ordinal, sname, gid, rec_ord, acc, cu.hex(), length, qual, owner))
            if owner:
                shared_rows.append((ordinal, gid, rec_ord, acc, cu.hex(), owner))
                held_in_mag += 1
                owner_set.add(owner)
            uids.append(cu)
            mag_bases += length
            contig_len.setdefault(cu, length)
            ordinal += 1

        g = listing.add(uids, 11, "prokaryote", qual)
        guid = ids.genome_uid(list(dict.fromkeys(uids)))
        if g is None:
            repeated_rows.append((gid, guid.hex(), first_of[guid]))
            continue
        g["circular"] = []
        genomes.append(g)
        first_of[guid] = gid
        sidecar_rows.append((
            stage_name, sname, gid, guid.hex(), len(g["contig_uids"]),
            mag_bases, qual, 11, 0, held_in_mag, ",".join(sorted(owner_set)),
        ))

    sidecar_rows.sort(key=lambda t: (t[2], t[3]))
    skani_rows = [
        (d["genome_id"], d["source_name"], d["reason"], d["ref_match"], d["ani"], d["af"])
        for d in sorted(skani_dropped, key=lambda x: (x["reason"], x["genome_id"]))
    ]

    first_src = ordered[0]
    objects_json = [{
        "path": str(obj_dir / fa_name),
        "source_name": first_src.get("source_name", "appler_2026_asgard"),
        "source_release": first_src.get("source_release", "figshare.26426707 (Appler et al. 2026 Nature)"),
        "license": first_src.get("license", "CC BY 4.0"),
        "redistributable": True,
        "retrieved_at": retrieved,
    }]

    make_inputs._write_tsv(work / "records.tsv", ASG1_RECORD_HEAD, rec_rows)
    make_inputs._write_tsv(work / "asg1_metadata.tsv", ASG1_SIDECAR_HEAD, sidecar_rows)
    make_inputs._write_tsv(work / "asg1_skani_dropped.tsv", SKANI_DROPPED_HEAD, skani_rows)
    make_inputs._write_tsv(
        work / "repeated_records.tsv",
        ["genome_id", "genome_uid", "first_genome_id"],
        repeated_rows,
    )
    make_inputs._write_tsv(
        work / "shared_records.tsv",
        ["ordinal", "genome_id", "record_ordinal", "accession", "contig_uid", "owner_stage"],
        shared_rows,
    )
    make_inputs._write_tsv(work / "trims.tsv", TRIM_HEAD, trims)
    make_inputs._write_tsv(work / "skipped_records.tsv", SKIPPED_HEAD, skipped)
    (work / "objects.json").write_text(json.dumps(objects_json, indent=1), encoding="utf-8")
    (work / "genomes.json").write_text(json.dumps(genomes), encoding="utf-8")

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
        "redistributable": True,
        "records": ordinal,
        "skipped_records": len(skipped),
        "skani_dropped": len(skani_rows),
        "genomes": len(genomes),
        "contigs": total_contigs,
        "bases": total_bases,
        "circular_genomes": 0,
        "repeated_records": len(repeated_rows),
        "trims": 0,
        "owned_by_stage": dict(sorted(owner_contigs.items())),
        "owned_bases_by_stage": dict(sorted(owner_bases.items())),
        "held_contigs": held_contigs,
        "held_bases": held_bases,
        "new_contigs": total_contigs - held_contigs,
        "new_bases": total_bases - held_bases,
        "root_stages": root_stages,
        "owner_inputs": sorted(owner_inputs),
    }


# --------------------------------------------------------------------------------------------------------------
# Network fetchers (DEC-50: <= 2 req/s, User-Agent: curl, stop after 5 consecutive refusals)


class _RateLimiter:
    def __init__(self, min_interval_s: float = 0.55):
        self.min_interval = min_interval_s
        self.last = 0.0
        self.refusals = 0

    def get_json(self, url: str, timeout: int = 60) -> dict:
        return json.loads(self.get(url, timeout=timeout))

    def get(self, url: str, timeout: int = 120) -> bytes:
        wait = self.min_interval - (time.monotonic() - self.last)
        if wait > 0:
            time.sleep(wait)
        self.last = time.monotonic()
        req = urllib.request.Request(url, headers={"User-Agent": "curl"})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
                self.refusals = 0
                return data
        except urllib.error.HTTPError as e:
            if e.code in (403, 429, 500, 502, 503, 504):
                self.refusals += 1
                if self.refusals >= 5:
                    raise Refused(f"5 consecutive refusals; last HTTP {e.code} on {url}") from e
            raise

    def download_file(
        self,
        url: str,
        dst: Path,
        *,
        expected_md5: str | None = None,
        timeout: int = 600,
    ) -> tuple[int, str, str]:
        dst.parent.mkdir(parents=True, exist_ok=True)
        tmp = dst.with_name(dst.name + ".part")
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
                    raise Refused(f"5 consecutive refusals; last HTTP {e.code} on {url}") from e
            raise
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise
        got_md5 = h_md5.hexdigest()
        if expected_md5 and got_md5 != expected_md5:
            tmp.unlink(missing_ok=True)
            raise Refused(f"{dst.name}: md5 mismatch (got {got_md5}, want {expected_md5})")
        tmp.rename(dst)
        return nbytes, h_sha.hexdigest(), got_md5


def _fetch_zenodo_record(rl: _RateLimiter, record_id: int, out_dir: Path) -> dict:
    """Inspect a Zenodo record via REST API, verify redistributable license, and download its files."""
    meta = rl.get_json(f"https://zenodo.org/api/records/{record_id}")
    lic = (meta.get("metadata", {}).get("license", {}) or {}).get("id", "")
    if not is_redistributable_license(lic):
        raise Refused(f"zenodo.{record_id}: license {lic!r} is not redistributable (DEC-51)")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"zenodo_{record_id}_meta.json").write_text(json.dumps(meta, indent=1, sort_keys=True), encoding="utf-8")
    files_rep: dict[str, dict] = {}
    for f in sorted(meta.get("files", []), key=lambda x: x.get("key", "")):
        fname = f["key"]
        url = f["links"]["self"]
        md5 = (f.get("checksum") or "").removeprefix("md5:") or None
        dst = out_dir / fname
        if not dst.is_file():
            nb, sha, got_md5 = rl.download_file(url, dst, expected_md5=md5)
            files_rep[fname] = {"bytes": nb, "sha256": sha, "md5": got_md5}
        else:
            files_rep[fname] = {"bytes": dst.stat().st_size, "sha256": ids.sha256_file(dst).hex()}
    return {"record_id": record_id, "license": lic, "files": files_rep}


def _fetch_figshare_article(
    rl: _RateLimiter,
    article_id: int,
    out_dir: Path,
    *,
    name_filter=None,
) -> dict:
    """Inspect a Figshare v2 article, verify redistributable license, and download matching files."""
    meta = rl.get_json(f"https://api.figshare.com/v2/articles/{article_id}")
    lic = (meta.get("license") or {}).get("name", "")
    if not is_redistributable_license(lic):
        raise Refused(f"figshare.{article_id}: license {lic!r} is not redistributable (DEC-51)")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"figshare_{article_id}_meta.json").write_text(json.dumps(meta, indent=1, sort_keys=True), encoding="utf-8")
    files_rep: dict[str, dict] = {}
    for f in sorted(meta.get("files", []), key=lambda x: x.get("name", "")):
        fname = f["name"]
        if name_filter is not None and not name_filter(fname):
            continue
        url = f["download_url"]
        md5 = f.get("computed_md5") or f.get("supplied_md5") or None
        dst = out_dir / fname
        if not dst.is_file():
            nb, sha, got_md5 = rl.download_file(url, dst, expected_md5=md5)
            files_rep[fname] = {"bytes": nb, "sha256": sha, "md5": got_md5}
        else:
            files_rep[fname] = {"bytes": dst.stat().st_size, "sha256": ids.sha256_file(dst).hex()}
    return {"article_id": article_id, "license": lic, "files": files_rep}


def fetch_ece2(dl_dir: Path) -> dict:
    """Fetch `ECE2` sources: Inocles (`zenodo.16656823`) and Tycheposons (`thackl/pro-tycheposons`)."""
    dl_dir = Path(dl_dir)
    dl_dir.mkdir(parents=True, exist_ok=True)
    rl = _RateLimiter(0.55)
    rep: dict[str, dict] = {}

    rep["inocles_zenodo"] = _fetch_zenodo_record(rl, ZENODO_INOCLES_GENOMES_ID, dl_dir / "inocles_zenodo")

    tyche_tar = dl_dir / "tycheposons_github" / "pro-tycheposons-main.tar.gz"
    if not tyche_tar.is_file():
        nb, sha, md5 = rl.download_file(
            f"https://github.com/{TYCHEPOSON_REPO}/archive/refs/heads/main.tar.gz",
            tyche_tar,
            timeout=300,
        )
        rep["tycheposons_github"] = {"bytes": nb, "sha256": sha, "md5": md5}
    else:
        rep["tycheposons_github"] = {
            "bytes": tyche_tar.stat().st_size,
            "sha256": ids.sha256_file(tyche_tar).hex(),
        }

    (dl_dir / "fetch_ece2_report.json").write_text(json.dumps(rep, indent=1, sort_keys=True), encoding="utf-8")
    return rep


def fetch_mir1(dl_dir: Path) -> dict:
    """Fetch `MIR1` sources: Mirusviricota (`28955240`, `20284713`) and PLVs (`21581355`)."""
    dl_dir = Path(dl_dir)
    dl_dir.mkdir(parents=True, exist_ok=True)
    rl = _RateLimiter(0.55)
    rep: dict[str, dict] = {}

    def _wanted_nuc(fname: str) -> bool:
        low = fname.lower()
        return any(tok in low for tok in ("genome", "contig", "mirus", "goev", "plv", "eve", "nanopore", "fasta", ".fna", ".fa"))

    rep["mirus_2025"] = _fetch_figshare_article(rl, FIGSHARE_MIRUS_2025_ID, dl_dir / "mirus_2025", name_filter=_wanted_nuc)
    rep["mirus_goev"] = _fetch_figshare_article(rl, FIGSHARE_MIRUS_GOEV_ID, dl_dir / "mirus_goev", name_filter=_wanted_nuc)
    rep["plv_2023"] = _fetch_figshare_article(rl, FIGSHARE_PLV_2023_ID, dl_dir / "plv_2023", name_filter=_wanted_nuc)
    rep["plv_2021_note"] = {
        "share_url": f"https://figshare.com/s/{FIGSHARE_PLV_2021_SHARE}",
        "note": "Private-share Figshare link; resolve numeric article_id & verify license on pod before placing FASTAs in dl/plv_2021/.",
    }
    (dl_dir / "fetch_mir1_report.json").write_text(json.dumps(rep, indent=1, sort_keys=True), encoding="utf-8")
    return rep


def fetch_asg1(dl_dir: Path, *, accessions_file: Path | None = None) -> dict:
    """Fetch `ASG1` sources: Appler 2026 Figshare `26426707` + optional NCBI Datasets v2 accession list."""
    dl_dir = Path(dl_dir)
    dl_dir.mkdir(parents=True, exist_ok=True)
    rl = _RateLimiter(0.55)
    rep: dict[str, dict] = {}

    rep["appler_2026_figshare"] = _fetch_figshare_article(rl, FIGSHARE_ASGARD_2026_ID, dl_dir / "appler_2026_figshare")

    if accessions_file is not None and Path(accessions_file).is_file():
        accs = sorted({
            line.strip().split()[0]
            for line in Path(accessions_file).read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")
        })
        ncbi_dir = dl_dir / "ncbi_batches"
        ncbi_dir.mkdir(parents=True, exist_ok=True)
        batches: dict[str, dict] = {}
        for i in range(0, len(accs), 20):
            chunk = accs[i : i + 20]
            bzip = ncbi_dir / f"asg1_batch_{i // 20:04d}.zip"
            if not bzip.is_file():
                url = (
                    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
                    + ",".join(chunk)
                    + "/download?include_annotation_type=GENOME_FASTA"
                )
                nb, sha, md5 = rl.download_file(url, bzip, timeout=300)
                batches[bzip.name] = {"count": len(chunk), "bytes": nb, "sha256": sha, "md5": md5}
            else:
                batches[bzip.name] = {
                    "count": len(chunk),
                    "bytes": bzip.stat().st_size,
                    "sha256": ids.sha256_file(bzip).hex(),
                }
        rep["ncbi_datasets"] = {"accessions": len(accs), "batches": batches}

    (dl_dir / "fetch_asg1_report.json").write_text(json.dumps(rep, indent=1, sort_keys=True), encoding="utf-8")
    return rep


# --------------------------------------------------------------------------------------------------------------
# Archive extraction helpers for staging from `dl_dir`


def _collect_fastas_recursive(base_dir: Path, scratch: Path, label: str) -> list[Path]:
    """Collect FASTA files from `base_dir`, unpacking any `.zip` or `.tar.gz`/`.tgz` archives into `scratch/label`."""
    found: list[Path] = []
    if not base_dir.is_dir():
        return found
    out_dir = scratch / label
    for p in sorted(base_dir.rglob("*")):
        if not p.is_file():
            continue
        low = p.name.lower()
        if _is_fasta_name(p.name):
            found.append(p)
        elif low.endswith(".zip"):
            out_dir.mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(p) as zf:
                for m in sorted(zf.namelist()):
                    if _is_fasta_name(m):
                        target = out_dir / Path(m).name
                        if not target.exists():
                            with zf.open(m) as src, open(target, "wb") as dst:
                                shutil.copyfileobj(src, dst)
                        found.append(target)
        elif low.endswith((".tar.gz", ".tgz", ".tar")):
            out_dir.mkdir(parents=True, exist_ok=True)
            with tarfile.open(p, "r:*") as tf:
                for m in sorted(tf.getmembers(), key=lambda x: x.name):
                    if m.isfile() and _is_fasta_name(m.name):
                        target = out_dir / Path(m.name).name
                        if not target.exists():
                            ef = tf.extractfile(m)
                            if ef is not None:
                                with open(target, "wb") as dst:
                                    shutil.copyfileobj(ef, dst)
                        found.append(target)
    return sorted(set(found))


def _concat_fastas(paths: list[Path], out_fa: Path) -> Path | None:
    if not paths:
        return None
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with open(out_fa, "wb") as w:
        for p in sorted(paths):
            with _open_fasta_in(p) as fh:
                data = fh.read()
                if data:
                    w.write(data)
                    if not data.endswith(b"\n"):
                        w.write(b"\n")
    return out_fa if out_fa.stat().st_size > 0 else None


def _prepare_ece2_sources(dl_dir: Path, scratch: Path) -> list[dict]:
    sources: list[dict] = []
    inocle_fastas = _collect_fastas_recursive(dl_dir / "inocles_zenodo", scratch, "inocles")
    inocle_cat = _concat_fastas(inocle_fastas, scratch / "inocles_combined.fna")
    if inocle_cat is not None:
        sources.append({
            "source_name": "inocles_kiguchi_2025",
            "source_release": f"10.5281/zenodo.{ZENODO_INOCLES_GENOMES_ID} (Kiguchi et al. 2025 Nat. Comms.)",
            "license": "CC BY 4.0",
            "redistributable": True,
            "fasta": inocle_cat,
            "output_name": "ece2_inocles.fna.gz",
            "element_class": "inocle",
        })

    tyche_fastas = _collect_fastas_recursive(dl_dir / "tycheposons_github", scratch, "tycheposons")
    tyche_cat = _concat_fastas(tyche_fastas, scratch / "tycheposons_combined.fna")
    if tyche_cat is not None:
        sources.append({
            "source_name": "tycheposons_hackl_2023",
            "source_release": f"github:{TYCHEPOSON_REPO} (Hackl & Chisholm 2023 Cell)",
            "license": "CC BY-SA 4.0",
            "redistributable": True,
            "fasta": tyche_cat,
            "output_name": "ece2_tycheposons.fna.gz",
            "element_class": "tycheposon",
        })
    return sources


def _prepare_mir1_sources(dl_dir: Path, scratch: Path) -> list[dict]:
    sources: list[dict] = []
    specs = [
        ("mirus_2025", "mirusviricota_medvedeva_2025", f"10.6084/m9.figshare.{FIGSHARE_MIRUS_2025_ID}", "CC BY 4.0", "mirusviricota", "mir1_mirusviricota_2025.fna.gz", "multi_contig_bins"),
        ("mirus_goev", "mirusviricota_gaia_2023_goev", f"10.6084/m9.figshare.{FIGSHARE_MIRUS_GOEV_ID}", "CC BY 4.0", "mirus_goev", "mir1_mirus_goev_2023.fna.gz", "multi_contig_bins"),
        ("plv_2021", "plv_bellas_2021", f"figshare.s.{FIGSHARE_PLV_2021_SHARE}", "CC BY 4.0", "plv", "mir1_plv_bellas_2021.fna.gz", "single_contig"),
        ("plv_2023", "plv_bellas_hackl_2023", f"10.6084/m9.figshare.{FIGSHARE_PLV_2023_ID}", "CC BY 4.0", "plv_eve", "mir1_plv_eve_2023.fna.gz", "single_contig"),
    ]
    for sub, sname, rel, lic, vcls, out_name, mode in specs:
        fastas = _collect_fastas_recursive(dl_dir / sub, scratch, sub)
        if not fastas:
            continue
        if mode == "multi_contig_bins":
            bins = [(_strip_fasta_suffix(p.name), p) for p in fastas]
            sources.append({
                "source_name": sname,
                "source_release": rel,
                "license": lic,
                "redistributable": True,
                "mode": "multi_contig_bins",
                "bin_fastas": bins,
                "output_name": out_name,
                "virus_class": vcls,
            })
        else:
            cat = _concat_fastas(fastas, scratch / f"{sub}_combined.fna")
            if cat is not None:
                sources.append({
                    "source_name": sname,
                    "source_release": rel,
                    "license": lic,
                    "redistributable": True,
                    "mode": "single_contig",
                    "fasta": cat,
                    "output_name": out_name,
                    "virus_class": vcls,
                })
    return sources


def _prepare_asg1_entries(dl_dir: Path, scratch: Path) -> list[dict]:
    entries: list[dict] = []
    for sub, sname, rel, lic in (
        ("appler_2026_figshare", "appler_2026_asgard", f"10.6084/m9.figshare.{FIGSHARE_ASGARD_2026_ID}", "CC BY 4.0"),
        ("ncbi_batches", "appler_2026_asgard_ncbi", "INSDC (Appler et al. 2026 Nature)", "INSDC (public)"),
    ):
        for p in _collect_fastas_recursive(dl_dir / sub, scratch, sub):
            entries.append({
                "genome_id": _strip_fasta_suffix(p.name),
                "source_name": sname,
                "source_release": rel,
                "license": lic,
                "redistributable": True,
                "assembly_quality": "mag_mq",
                "fasta": p,
                "size_bytes": p.stat().st_size,
            })
    return sorted(entries, key=lambda x: x["genome_id"])


def _parse_kv(items: list[str]) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for item in items:
        if "=" not in item:
            raise Refused(f"--owner-inputs expects NAME=DIR, got {item!r}")
        k, v = item.split("=", 1)
        out[k] = Path(v)
    return out


# --------------------------------------------------------------------------------------------------------------
# CLI


def _add_stage_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--dl", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--retrieved", required=True)
    p.add_argument("--root", type=Path, default=None)
    p.add_argument("--root-stages", default=DEFAULT_ROOT_STAGES)
    p.add_argument("--owner-inputs", action="append", default=[])
    p.add_argument("--final-dir", type=Path, default=None)
    p.add_argument("--out-report", type=Path, default=None)


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description="Card C-70: fetch and stage ECE2, MIR1, and ASG1.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_fe = sub.add_parser("fetch-ece2")
    p_fe.add_argument("--dl", type=Path, required=True)
    p_se = sub.add_parser("stage-ece2")
    _add_stage_args(p_se)

    p_fm = sub.add_parser("fetch-mir1")
    p_fm.add_argument("--dl", type=Path, required=True)
    p_sm = sub.add_parser("stage-mir1")
    _add_stage_args(p_sm)

    p_fa = sub.add_parser("fetch-asg1")
    p_fa.add_argument("--dl", type=Path, required=True)
    p_fa.add_argument("--accessions", type=Path, default=None)
    p_sa = sub.add_parser("stage-asg1")
    _add_stage_args(p_sa)
    p_sa.add_argument("--skani-ref-tsv", type=Path, default=None)
    p_sa.add_argument("--skani-ref-list", type=Path, default=None)
    p_sa.add_argument("--skani-self-dedup", action="store_true")
    p_sa.add_argument("--skani-ani", type=float, default=96.0)
    p_sa.add_argument("--skani-af", type=float, default=50.0)
    p_sa.add_argument("--skani-threads", type=int, default=8)

    args = ap.parse_args(argv)
    if args.cmd == "fetch-ece2":
        print(json.dumps(fetch_ece2(args.dl), indent=1, sort_keys=True))
        return 0
    if args.cmd == "fetch-mir1":
        print(json.dumps(fetch_mir1(args.dl), indent=1, sort_keys=True))
        return 0
    if args.cmd == "fetch-asg1":
        print(json.dumps(fetch_asg1(args.dl, accessions_file=args.accessions), indent=1, sort_keys=True))
        return 0

    root_stages = [s.strip() for s in args.root_stages.split(",") if s.strip()]
    owner_inputs = _parse_kv(args.owner_inputs)
    scratch = args.dl / f".extract_scratch_{args.cmd}"
    shutil.rmtree(scratch, ignore_errors=True)
    try:
        if args.cmd == "stage-ece2":
            sources = _prepare_ece2_sources(args.dl, scratch)
            rep = stage_ece2(
                sources, args.out, root=args.root, root_stages=root_stages,
                owner_inputs=owner_inputs, retrieved=args.retrieved, final_dir=args.final_dir,
            )
        elif args.cmd == "stage-mir1":
            sources = _prepare_mir1_sources(args.dl, scratch)
            rep = stage_mir1(
                sources, args.out, root=args.root, root_stages=root_stages,
                owner_inputs=owner_inputs, retrieved=args.retrieved, final_dir=args.final_dir,
            )
        else:
            entries = _prepare_asg1_entries(args.dl, scratch)
            ref_matches: dict[str, dict] = {}
            if args.skani_ref_tsv is not None:
                ref_matches = parse_skani_matches(args.skani_ref_tsv, min_ani=args.skani_ani, min_af=args.skani_af)
            elif args.skani_ref_list is not None:
                ref_matches = run_skani_dist_against_ref(
                    entries, args.skani_ref_list, scratch / "skani_ref",
                    min_ani=args.skani_ani, min_af=args.skani_af, threads=args.skani_threads,
                )
            extra_dropped: list[dict] = []
            if args.skani_self_dedup and len(entries) > 1:
                entries, extra_dropped = run_skani_self_dedup(
                    entries, scratch / "skani_self",
                    min_ani=args.skani_ani, min_af=args.skani_af, threads=args.skani_threads,
                )
            rep = stage_asg1(
                entries, args.out, root=args.root, root_stages=root_stages,
                owner_inputs=owner_inputs, retrieved=args.retrieved, final_dir=args.final_dir,
                skani_ref_matches=ref_matches, skani_min_ani=args.skani_ani,
                skani_extra_dropped=extra_dropped,
            )
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    if args.out_report is not None:
        args.out_report.write_text(json.dumps(rep, indent=1, sort_keys=True), encoding="utf-8")
    print(json.dumps(rep, indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
