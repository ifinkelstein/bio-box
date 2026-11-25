"""Card C-68 (DEC-50, DEC-59, DEC-65, DEC-69, DEC-71):
Pre-build and pod staging driver for:
  1. `ET1` — Spliced protist & plant mRNA nucleotide transcriptomes (`MMETSP` Johnson et al. 2019
     Zenodo `3247846` / `1212585` + `1KP` Leebens-Mack et al. 2019 GigaDB `100910` / `100627` /
     CyVerse `oneKP_capstone_2019`), >= 1,500 nt cDNA contigs, deterministic exact-sequence and
     optional Trinity isoform deduplication, verified eukaryotic/protist `transl_table` taxon routing
     checked against `biobox.search.orfs.check_code` so `pack build` never crashes, and
     `ET1_left_for_later.tsv.gz` for `< 1,500` nt and deferred-code transcripts per `DEC-59`.
  2. `LGNRNAP1` — Targeted `Logan` v1 SRA contig bait-pull (`logan50` -> `s3://logan-pub/c/`)
     reusing `c64_logan_rt_bait.py`'s streaming architecture with `W-40`'s two-namespace `ssrnap`
     super-panel (`examples/rnap/ssrnap_cutoffs.tsv`, 30 targets + 16 outgroup decoys).

Subcommands:
    fetch-et1         OUT_DIR [--skip-mmetsp] [--skip-onekp] [--onekp-source gigadb|cyverse]
    stage-et1         SAMPLES_TSV OUT_META_DIR OUT_FASTA_DIR [--final-fasta-dir DIR]
                      [--fasta-name et1_transcriptomes.fna.gz] [--min-len 1500]
                      [--dedup-isoforms] [--unsupported-code-policy defer|fallback_table6]
                      [--root ROOT] [--root-stages ...] [--owner-inputs NAME=DIR ...]
    build-ssrnap-bait OUT_DIR [--pfam-hmm PATH] [--hmms-dir PATH]
    bait-lgnrnap1     BAIT_HMM OUT_DIR [--decoy-hmm PATH] [--cutoffs-tsv PATH] [--cpus 6]
                      [--batch-size 50000] [--max-compressed-bytes BYTES] [--max-wall-s SEC]
                      [--start-byte 0] [--index-path PATH] [--max-members-per-cluster 50]
                      [--max-frames 200] [--max-runs 500] [--skip-resolve]
    stage-lgnrnap1    RESOLVED_JSON OUT_META_DIR OUT_FASTA_DIR [--final-fasta-dir DIR]
                      [--fasta-name lgnrnap1_baited.fna.gz] [--min-len 2000]
                      [--root ROOT] [--root-stages ...] [--owner-inputs NAME=DIR ...]
                      [--skip-ena] [--max-run-bytes 100000000] [--workers 4]
"""
from __future__ import annotations

import argparse
import bz2
import csv
import gzip
import hashlib
import io
import json
import re
import subprocess
import sys
import tarfile
import threading
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pyhmmer

sys.path.insert(0, str(Path(__file__).resolve().parent))

import c64_logan_rt_bait as c64  # noqa: E402
import make_inputs  # noqa: E402
from biobox import ids  # noqa: E402
from biobox.search import orfs  # noqa: E402

Refused = make_inputs.Refused

UA = "curl/8.5.0"
MIN_STOR_FILE_BYTES = 64 * 1024 * 1024  # AGENTS.md §7.3

# Primary-source-verified ET1 upstream URLs & checksums (planning/measurements/2026-11-24_c68_et1_prep.md)
MMETSP_ZENODO_RECORD = "3247846"  # Johnson et al. 2019 GigaScience; supersedes v1 record 1212585
MMETSP_FASTA_TAR_NAME = "mmetsp_dib_trinity_zenodo.fasta.tar.gz"
MMETSP_FASTA_TAR_URL = (
    f"https://zenodo.org/api/records/{MMETSP_ZENODO_RECORD}/files/{MMETSP_FASTA_TAR_NAME}/content"
)
MMETSP_FASTA_TAR_MD5 = "78c908575cb2c75c2d0ae71f817175b9"
MMETSP_NAMEMAPS_TAR_NAME = "mmetsp_dib_trinity2.2.0_contigIDnamemaps.tar.gz"
MMETSP_NAMEMAPS_TAR_URL = (
    f"https://zenodo.org/api/records/{MMETSP_ZENODO_RECORD}/files/{MMETSP_NAMEMAPS_TAR_NAME}/content"
)

ONEKP_GIGADB_DOI = "10.5524/100910"  # Scrubbed v2 of 10.5524/100627 (30 contaminated samples removed)
ONEKP_GIGADB_FTP_BASE = "https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100910"
ONEKP_CYVERSE_BASE = (
    "https://data.cyverse.org/dav-anon/iplant/commons/community_released/oneKP_capstone_2019"
)

ET1_MIN_LEN = 1_500
LGNRNAP1_MIN_LEN = 2_000

DEFAULT_PFAM_HMM = Path("/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2/Pfam-A.hmm.gz")
DEFAULT_W40_HMMS_DIR = Path("/tmp/if2453/w40/hmms")
DEFAULT_SSRNAP_CUTOFFS_TSV = (
    Path(__file__).resolve().parents[3] / "examples" / "rnap" / "ssrnap_cutoffs.tsv"
)

# W-40 two-namespace ssRNAP super-panel Pfam accessions (examples/rnap/ssrnap_panel.sh)
SSRNAP_PFAM_PALM_ACCS = (
    "PF00940",  # RNA_pol
    "PF14700",  # RPOL_N
    "PF29426",  # Phage_RNApol
    "PF21769",  # N4 vRNAP quartet
    "PF21867",
    "PF21894",
    "PF21624",
)
SSRNAP_PFAM_DPBB_ACCS = ("PF05183",)  # RdRP / QDE-1-like DPBB
SSRNAP_PFAM_TARGET_ACCS = SSRNAP_PFAM_PALM_ACCS + SSRNAP_PFAM_DPBB_ACCS

SSRNAP_PFAM_DECOY_POL1_ACCS = ("PF00476",)  # DNA_pol_A
SSRNAP_PFAM_DECOY_CELL_RNAP_ACCS = (
    "PF04997", "PF00623", "PF04983", "PF05000", "PF04998", "PF04992", "PF04990",
    "PF04563", "PF04561", "PF04565", "PF04566", "PF10385", "PF04567", "PF00562", "PF04560",
)
SSRNAP_PFAM_DECOY_ACCS = SSRNAP_PFAM_DECOY_POL1_ACCS + SSRNAP_PFAM_DECOY_CELL_RNAP_ACCS

VALID_DNA_RE = re.compile(rb"^[ACGTURYKMSWBDHVNacgturykmswbdhvn]+$")
_TRINITY_LOCUS_RE = re.compile(r"^((?:TRINITY_)?(?:DN\d+_)?c\d+_g\d+)_i\d+\b")

# --------------------------------------------------------------------------------------------------
# Verified eukaryotic/protist nuclear translation-table taxon mapping
# --------------------------------------------------------------------------------------------------

TRANSL_TABLE_DESCRIPTIONS: dict[int, str] = {
    1: "Standard Eukaryotic Nuclear (TAA=*, TAG=*, TGA=*)",
    6: "Ciliate, Dasycladacean and Hexamita Nuclear (TAA=Q, TAG=Q, TGA=*)",
    10: "Euplotid Nuclear (TAA=*, TAG=*, TGA=C)",
    29: "Mesodinium Nuclear (TAA=Y, TAG=Y, TGA=*)",
    30: "Peritrich Nuclear (TAA=E, TAG=E, TGA=*)",
}

# Specificity-ordered regex patterns over lowercased "organism ; taxonomy" text.
# Note: Heterotrichea (Stentor, Blepharisma, Fabrea, Climacostomum, Condylostoma) and Karyorelictea
# are ciliates that do NOT use Table 6 (they use Standard Table 1 or context-dependent stops).
_TABLE_29_RE = re.compile(
    r"\b(?:mesodinium|myrionecta|mesodiniidae)\b",
    re.IGNORECASE,
)
_TABLE_30_RE = re.compile(
    r"\b(?:peritrichia|peritrichida|vorticella|vorticellidae|epistylis|epistylididae|"
    r"carchesium|zoothamnium|opercularia|campanella|pseudovorticella|trichodina)\b",
    re.IGNORECASE,
)
_TABLE_10_RE = re.compile(
    r"\b(?:euplotes|euplotida|euplotidae|euplotia|certesia|diophrys|uronychia|gastrocirrhus|"
    r"aspidisca|moneuplotes)\b",
    re.IGNORECASE,
)
_CILIATE_TABLE_1_EXEMPT_RE = re.compile(
    r"\b(?:heterotrichea|heterotrichida|stentor|blepharisma|fabrea|climacostomum|spirostomum|"
    r"condylostoma|karyorelictea|loxodes|parduczia|nyctotherus)\b",
    re.IGNORECASE,
)
_TABLE_6_RE = re.compile(
    r"\b(?:ciliophora|oligohymenophorea|spirotrichea|litostomatea|phyllopharyngea|"
    r"prostomatea|nassophorea|colpodea|plagiopylea|armophorea|"
    r"tetrahymena|paramecium|oxytricha|stylonychia|sterkiella|pseudocohnilembus|"
    r"ichthyophthirius|uronema|aristerostoma|strombidium|strombidinopsis|tintinnopsis|"
    r"favella|schmidingerella|eutintinnus|tontonia|laboea|balanion|didinium|colpoda|"
    r"dasycladales|dasycladaceae|dasycladus|acetabularia|batophora|cymopolia|"
    r"hexamita|hexamitidae|spironucleus|tretomita)\b",
    re.IGNORECASE,
)

UNSUPPORTED_CODE_POLICIES = ("defer", "fallback_table6")


def check_orfs_code_supported(code: int) -> bool:
    """Return True iff `biobox.search.orfs.check_code(code)` accepts `code` (required by `pack build`)."""
    try:
        orfs.check_code(code)
        return True
    except ValueError:
        return False


def resolve_taxon_genetic_code(
    organism: str,
    taxonomy: str = "",
    *,
    unsupported_code_policy: str = "defer",
) -> dict:
    """Map an organism name and lineage string to its biological NCBI nuclear translation table
    and a `pack build`-safe `effective_transl_table` verified against `orfs.check_code`.

    `pyrodigal >= 3.4` (3.7.1 in `uv.lock`) supports tables 1..6, 9..16, 21..33 (including 29 & 30),
    whereas `biobox.search.orfs.check_code` (pyhmmer Easel + `_EXPLICIT`) currently supports
    1..6, 9..16, 21..25 and raises `ValueError` on 29 & 30 during `pack build`.

    When `orfs.check_code(bio_code)` fails:
      - `unsupported_code_policy="defer"` (default): marks the sample `action="defer"` so it is
        recorded in `ET1_left_for_later.tsv.gz` with zero mistranslated residues.
      - `unsupported_code_policy="fallback_table6"`: maps 29/30 (both TGA-only stop, TAA/TAG sense)
        to Table 6 (also TGA-only stop, TAA/TAG -> Q), preserving unbroken ORF coordinates.
    """
    if unsupported_code_policy not in UNSUPPORTED_CODE_POLICIES:
        raise ValueError(
            f"unsupported_code_policy must be one of {UNSUPPORTED_CODE_POLICIES}, got {unsupported_code_policy!r}"
        )

    text = f"{organism or ''} ; {taxonomy or ''}".strip()

    if _TABLE_29_RE.search(text):
        bio_code, rule = 29, "mesodinium_table29"
    elif _TABLE_30_RE.search(text):
        bio_code, rule = 30, "peritrichia_table30"
    elif _TABLE_10_RE.search(text):
        bio_code, rule = 10, "euplotid_table10"
    elif _CILIATE_TABLE_1_EXEMPT_RE.search(text):
        bio_code, rule = 1, "heterotrich_or_karyorelict_table1"
    elif _TABLE_6_RE.search(text):
        bio_code, rule = 6, "ciliate_dasyclad_hexamita_table6"
    else:
        bio_code, rule = 1, "standard_eukaryote_table1"

    supported = check_orfs_code_supported(bio_code)
    if supported:
        return {
            "bio_transl_table": bio_code,
            "effective_transl_table": bio_code,
            "orfs_supported": True,
            "action": "stage",
            "matched_rule": rule,
            "defer_reason": "",
            "needed_capability": "",
            "code_note": "",
        }

    if unsupported_code_policy == "fallback_table6" and bio_code in (29, 30):
        # Verify Table 6 itself is supported before falling back
        orfs.check_code(6)
        return {
            "bio_transl_table": bio_code,
            "effective_transl_table": 6,
            "orfs_supported": False,
            "action": "stage",
            "matched_rule": rule,
            "defer_reason": "",
            "needed_capability": "orfs_py_explicit_tables_29_30",
            "code_note": f"fallback_from_{bio_code}_to_6_until_orfs_py_supports_{bio_code}",
        }

    return {
        "bio_transl_table": bio_code,
        "effective_transl_table": None,
        "orfs_supported": False,
        "action": "defer",
        "matched_rule": rule,
        "defer_reason": f"transl_table_{bio_code}_unsupported_by_orfs_py",
        "needed_capability": "orfs_py_explicit_tables_29_30",
        "code_note": f"deferred_table_{bio_code}_not_in_orfs_supported_codes",
    }


# --------------------------------------------------------------------------------------------------
# ET1 cDNA filtering, deterministic isoform/exact-sequence dedup, and staging
# --------------------------------------------------------------------------------------------------


def parse_trinity_locus_id(header: str) -> str | None:
    """Extract the Trinity locus prefix (e.g. `TRINITY_DN123_c0_g1` from `TRINITY_DN123_c0_g1_i2`),
    or return None if the header is not a Trinity isoform ID."""
    tok = header.lstrip(">").split()[0] if header.strip() else ""
    m = _TRINITY_LOCUS_RE.match(tok)
    return m.group(1) if m else None


def normalize_cdna_seq(raw_seq: bytes) -> bytes | None:
    """Uppercase, convert RNA `U` -> `T`, strip whitespace, and validate IUPAC DNA alphabet."""
    cleaned = re.sub(rb"\s+", b"", raw_seq).upper().replace(b"U", b"T")
    if not cleaned or not VALID_DNA_RE.match(cleaned):
        return None
    return cleaned


def _iter_fasta_file(path: Path) -> Iterator[tuple[str, bytes]]:
    name = path.name.lower()
    if name.endswith(".gz"):
        fh = gzip.open(path, "rb")
    elif name.endswith(".bz2"):
        fh = bz2.open(path, "rb")
    else:
        fh = open(path, "rb")
    try:
        yield from c64._iter_fasta_bytes_stream(fh)
    finally:
        fh.close()


def filter_and_dedup_cdna_records(
    records: Iterable[tuple[str, bytes]],
    *,
    sample_id: str,
    source_name: str,
    organism: str,
    min_len: int = ET1_MIN_LEN,
    dedup_isoforms: bool = False,
) -> tuple[list[dict], list[dict], dict]:
    """Filter spliced cDNA transcripts >= `min_len` nt and deduplicate deterministically.

    Returns `(kept_records, left_for_later_rows, stats)`:
      - Invalid/empty sequences are logged to `left_for_later_rows` with `reason="invalid_or_empty_seq"`.
      - Sequences `< min_len` nt are logged to `left_for_later_rows` with `reason="below_min_len_1500nt"`
        (parameterized by `min_len`).
      - Exact duplicate sequences within a sample (same `contig_uid`) are collapsed deterministically
        (preferring lexicographically smallest `transcript_id`).
      - When `dedup_isoforms=True`, Trinity isoforms sharing the same locus (`TRINITY_DN..._c..._g...`)
        keep only the longest isoform (tie-break by smallest `contig_uid`, then `transcript_id`),
        logging shorter isoforms to `left_for_later_rows` with `reason="shorter_trinity_isoform"`.
    """
    stats = Counter()
    left_rows: list[dict] = []
    passing: list[dict] = []

    for hdr, raw_seq in records:
        tid = hdr.lstrip(">").split()[0] if hdr.strip() else ""
        seq = normalize_cdna_seq(raw_seq)
        if seq is None:
            stats["invalid_seq"] += 1
            left_rows.append({
                "sample_id": sample_id,
                "source_name": source_name,
                "organism": organism,
                "transcript_id": tid,
                "length": len(raw_seq.strip()),
                "reason": "invalid_or_empty_seq",
                "needed_capability": "manual_review",
            })
            continue
        slen = len(seq)
        stats["total_valid"] += 1
        stats["total_valid_bp"] += slen
        if slen < min_len:
            stats["below_min_len"] += 1
            stats["below_min_len_bp"] += slen
            left_rows.append({
                "sample_id": sample_id,
                "source_name": source_name,
                "organism": organism,
                "transcript_id": tid,
                "length": slen,
                "reason": f"below_min_len_{min_len}nt",
                "needed_capability": "et_short_transcripts_tranche",
            })
            continue
        cu = ids.contig_uid(seq)
        passing.append({
            "transcript_id": tid,
            "header": hdr.lstrip(">").strip(),
            "seq": seq,
            "length": slen,
            "contig_uid": cu,
            "locus_id": parse_trinity_locus_id(hdr),
        })

    # Deterministic exact-sequence deduplication within sample
    passing.sort(key=lambda r: (r["contig_uid"], r["transcript_id"]))
    dedup_exact: list[dict] = []
    seen_cu: set[bytes] = set()
    for rec in passing:
        if rec["contig_uid"] in seen_cu:
            stats["within_sample_exact_dup"] += 1
            continue
        seen_cu.add(rec["contig_uid"])
        dedup_exact.append(rec)

    if not dedup_isoforms:
        dedup_exact.sort(key=lambda r: (r["transcript_id"], r["contig_uid"]))
        stats["kept"] = len(dedup_exact)
        stats["kept_bp"] = sum(r["length"] for r in dedup_exact)
        return dedup_exact, left_rows, dict(stats)

    # Optional Trinity locus isoform dedup: keep longest isoform per locus
    by_locus: dict[str, list[dict]] = defaultdict(list)
    non_trinity: list[dict] = []
    for rec in dedup_exact:
        if rec["locus_id"] is not None:
            by_locus[rec["locus_id"]].append(rec)
        else:
            non_trinity.append(rec)

    kept: list[dict] = list(non_trinity)
    for locus_id in sorted(by_locus):
        isoforms = sorted(
            by_locus[locus_id],
            key=lambda r: (-r["length"], r["contig_uid"], r["transcript_id"]),
        )
        kept.append(isoforms[0])
        for shorter in isoforms[1:]:
            stats["shorter_isoform_dropped"] += 1
            stats["shorter_isoform_bp"] += shorter["length"]
            left_rows.append({
                "sample_id": sample_id,
                "source_name": source_name,
                "organism": organism,
                "transcript_id": shorter["transcript_id"],
                "length": shorter["length"],
                "reason": "shorter_trinity_isoform",
                "needed_capability": "et_all_isoforms_tranche",
            })

    kept.sort(key=lambda r: (r["transcript_id"], r["contig_uid"]))
    stats["kept"] = len(kept)
    stats["kept_bp"] = sum(r["length"] for r in kept)
    return kept, left_rows, dict(stats)


def _write_deterministic_gz_fasta_header(gz_out, contig_uid_hex: str, seq: bytes) -> None:
    gz_out.write(b">" + contig_uid_hex.encode("ascii") + b"\n")
    for i in range(0, len(seq), 80):
        gz_out.write(seq[i : i + 80] + b"\n")


def stage_et1(
    *,
    samples: list[dict],
    out_meta_dir: Path,
    out_fasta_dir: Path,
    final_fasta_dir: Path | None = None,
    fasta_name: str = "et1_transcriptomes.fna.gz",
    min_len: int = ET1_MIN_LEN,
    dedup_isoforms: bool = False,
    unsupported_code_policy: str = "defer",
    root: Path | None = None,
    root_stages: Iterable[str] = (),
    owner_inputs: dict[str, Path] | None = None,
    retrieved_at: str = "2026-11-24T00:00:00Z",
    threads: int = 2,
    memory_limit: str = "4GB",
) -> dict:
    """Stage `ET1` spliced protist & plant cDNA transcriptomes conforming to `s2_check_inputs.py`.

    Each sample dict must provide:
      - `sample_id` (str), `source_name` (str), `source_release` (str), `license` (str)
      - `organism` (str), `taxonomy` (str, optional)
      - either `records` (iterable of `(header, seq_bytes)`) or `fasta_path` (Path/str)

    All kept `>= min_len` nt transcripts of a sample are grouped into one multi-contig `GenomeInput`
    with `kind="prokaryote"`, `genetic_code=effective_transl_table`, `assembly_quality="contig"`,
    and `circular=[]`. Because transcriptome assemblies typically total 1–40 Mbp (`>= 100 kb`),
    `biobox.l1.genecall` trains single-mode Pyrodigal per sample using that sample's `genetic_code`,
    and every emitted `genetic_code` is verified against `biobox.search.orfs.check_code` so
    `pack build` never fails.
    """
    out_meta_dir = Path(out_meta_dir)
    out_fasta_dir = Path(out_fasta_dir)
    out_meta_dir.mkdir(parents=True, exist_ok=True)
    out_fasta_dir.mkdir(parents=True, exist_ok=True)

    earlier_owners = c64._load_earlier_contig_owners(
        root, root_stages, owner_inputs, threads=threads, memory_limit=memory_limit
    )

    ordered_samples = sorted(samples, key=lambda s: (str(s["source_name"]), str(s["sample_id"])))
    fasta_path = out_fasta_dir / fasta_name
    tmp_fasta = fasta_path.with_name(fasta_path.name + ".tmp")

    written_uids: set[bytes] = set()
    uid_code: dict[bytes, int] = {}
    uid_len: dict[bytes, int] = {}
    genomes_out: list[dict] = []
    seen_genome_uids: set[bytes] = set()

    et1_meta_rows: list[list] = []
    record_rows: list[list] = []
    repeated_rows: list[list] = []
    left_rows: list[dict] = []

    counts = Counter()
    code_counts = Counter()
    bio_code_counts = Counter()
    source_counts = Counter()
    total_bp = 0

    with open(tmp_fasta, "wb") as raw_fh:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_fh, compresslevel=6, mtime=0) as gz_out:
            for sample in ordered_samples:
                sid = str(sample["sample_id"])
                src_name = str(sample["source_name"])
                org = str(sample.get("organism", ""))
                tax = str(sample.get("taxonomy", ""))

                routing = resolve_taxon_genetic_code(
                    org, tax, unsupported_code_policy=unsupported_code_policy
                )
                bio_code = routing["bio_transl_table"]
                eff_code = routing["effective_transl_table"]
                bio_code_counts[bio_code] += 1

                if routing["action"] == "defer":
                    counts["deferred_samples_unsupported_code"] += 1
                    left_rows.append({
                        "sample_id": sid,
                        "source_name": src_name,
                        "organism": org,
                        "transcript_id": "*",
                        "length": 0,
                        "reason": routing["defer_reason"],
                        "needed_capability": routing["needed_capability"],
                    })
                    continue

                assert eff_code is not None
                orfs.check_code(eff_code)

                raw_iter = (
                    sample["records"]
                    if "records" in sample and sample["records"] is not None
                    else _iter_fasta_file(Path(sample["fasta_path"]))
                )
                kept, sample_left, st = filter_and_dedup_cdna_records(
                    raw_iter,
                    sample_id=sid,
                    source_name=src_name,
                    organism=org,
                    min_len=min_len,
                    dedup_isoforms=dedup_isoforms,
                )
                left_rows.extend(sample_left)
                counts["below_min_len_transcripts"] += st.get("below_min_len", 0)
                counts["below_min_len_bp"] += st.get("below_min_len_bp", 0)
                counts["within_sample_exact_dups"] += st.get("within_sample_exact_dup", 0)
                counts["shorter_isoforms_dropped"] += st.get("shorter_isoform_dropped", 0)

                sample_cus: list[bytes] = []
                sample_bp = 0
                for rec in kept:
                    cu = rec["contig_uid"]
                    # DEC-27 / s2_check_inputs invariant: a contig shared across two genomes must
                    # agree on genetic_code. If an identical cDNA appears in another sample with a
                    # different genetic_code, drop it from the later sample and log it.
                    if cu in uid_code and uid_code[cu] != eff_code:
                        counts["cross_code_dup_dropped"] += 1
                        repeated_rows.append([
                            sid, rec["transcript_id"], cu.hex(),
                            f"cross_code_conflict_{uid_code[cu]}_vs_{eff_code}",
                        ])
                        continue
                    if cu not in written_uids:
                        written_uids.add(cu)
                        uid_code[cu] = eff_code
                        uid_len[cu] = rec["length"]
                        _write_deterministic_gz_fasta_header(gz_out, cu.hex(), rec["seq"])
                        total_bp += rec["length"]
                    else:
                        counts["cross_sample_exact_dups"] += 1
                        repeated_rows.append([
                            sid, rec["transcript_id"], cu.hex(), "duplicate_contig_uid_across_samples",
                        ])
                    sample_cus.append(cu)
                    sample_bp += rec["length"]
                    record_rows.append([
                        sid, src_name, rec["transcript_id"], cu.hex(), rec["length"],
                        bio_code, eff_code,
                    ])

                if not sample_cus:
                    counts["empty_samples_after_filter"] += 1
                    continue

                # Canonical sorted unique contig_uids for genome_uid
                sample_cus = sorted(set(sample_cus))
                guid = ids.genome_uid(sample_cus)
                if guid in seen_genome_uids:
                    counts["duplicate_sample_genomes"] += 1
                    repeated_rows.append([sid, "*", guid.hex(), "duplicate_genome_uid"])
                    continue
                seen_genome_uids.add(guid)

                genomes_out.append({
                    "contig_uids": [c.hex() for c in sample_cus],
                    "genetic_code": eff_code,
                    "kind": "prokaryote",
                    "circular": [],
                    "assembly_quality": "contig",
                })
                code_counts[eff_code] += 1
                source_counts[src_name] += 1
                et1_meta_rows.append([
                    "ET1", sid, src_name, guid.hex(), org, tax,
                    len(sample_cus), sample_bp, bio_code, eff_code,
                    routing["matched_rule"], routing["code_note"],
                ])

    if not genomes_out:
        tmp_fasta.unlink(missing_ok=True)
        raise Refused("no qualifying >= min_len cDNA transcripts staged for ET1")

    fasta_bytes = tmp_fasta.stat().st_size
    if str(fasta_path.resolve()).startswith("/stor/") and fasta_bytes < MIN_STOR_FILE_BYTES:
        tmp_fasta.unlink(missing_ok=True)
        raise Refused(
            f"Refusing to write {fasta_bytes} B (< 64 MB) to /stor path {fasta_path} (AGENTS.md §7.3); "
            f"stage under /tmp/if2453/s2/inputs/ET1/ instead."
        )
    tmp_fasta.replace(fasta_path)
    fasta_sha256 = c64._sha256_file(fasta_path)

    shared_rows: list[list] = []
    shared_by_stage = Counter()
    shared_bases_by_stage = Counter()
    for uid in sorted(written_uids):
        owner = earlier_owners.get(uid)
        if owner is not None:
            ln = uid_len[uid]
            shared_rows.append([uid.hex(), owner, ln])
            shared_by_stage[owner] += 1
            shared_bases_by_stage[owner] += ln

    genomes_path = out_meta_dir / "genomes.json"
    genomes_path.write_text(json.dumps(genomes_out, separators=(",", ":")) + "\n", encoding="utf-8")

    obj_staging = [{
        "path": str(fasta_path),
        "source_name": "mmetsp_1kp_et1",
        "source_release": "mmetsp-zenodo3247846_1kp-gigadb100910",
        "license": "CC BY 4.0 (MMETSP) / CC0 1.0 (1KP GigaDB)",
        "redistributable": True,
        "retrieved_at": retrieved_at,
    }]
    obj_staging_path = out_meta_dir / "objects_staging.json"
    obj_staging_path.write_text(json.dumps(obj_staging, indent=2) + "\n", encoding="utf-8")

    final_path = (final_fasta_dir / fasta_name) if final_fasta_dir is not None else fasta_path
    obj_final = [{**obj_staging[0], "path": str(final_path)}]
    obj_final_path = out_meta_dir / "objects.json"
    obj_final_path.write_text(json.dumps(obj_final, indent=2) + "\n", encoding="utf-8")

    def _write_tsv(path: Path, header: list[str], rows: list[list]) -> str:
        with open(path, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(header)
            w.writerows(rows)
        return c64._sha256_file(path)

    et1_meta_rows.sort(key=lambda r: (r[2], r[1]))
    record_rows.sort(key=lambda r: (r[1], r[0], r[2], r[3]))
    repeated_rows.sort(key=lambda r: (r[0], r[1], r[2]))
    left_rows.sort(key=lambda r: (r["reason"], r["source_name"], r["sample_id"], r["transcript_id"]))

    meta_sha = _write_tsv(
        out_meta_dir / "et1_metadata.tsv",
        [
            "stage", "sample_id", "source_name", "genome_uid", "organism", "taxonomy",
            "n_contigs", "total_bp", "bio_transl_table", "effective_transl_table",
            "matched_rule", "code_note",
        ],
        et1_meta_rows,
    )
    rec_sha = _write_tsv(
        out_meta_dir / "records.tsv",
        ["sample_id", "source_name", "transcript_id", "contig_uid", "length", "bio_transl_table", "effective_transl_table"],
        record_rows,
    )
    rep_sha = _write_tsv(
        out_meta_dir / "repeated_records.tsv",
        ["sample_id", "transcript_id", "uid", "reason"],
        repeated_rows,
    )
    shr_sha = _write_tsv(
        out_meta_dir / "shared_records.tsv",
        ["contig_uid", "owner_stage", "length"],
        shared_rows,
    )

    left_path = out_meta_dir / "ET1_left_for_later.tsv.gz"
    with open(left_path, "wb") as raw_l:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_l, compresslevel=6, mtime=0) as gz_l:
            text_l = io.TextIOWrapper(gz_l, encoding="utf-8", newline="")
            w = csv.writer(text_l, delimiter="\t")
            w.writerow(["sample_id", "source_name", "organism", "transcript_id", "length", "reason", "needed_capability"])
            for r in left_rows:
                w.writerow([
                    r["sample_id"], r["source_name"], r["organism"],
                    r["transcript_id"], r["length"], r["reason"], r["needed_capability"],
                ])
            text_l.flush()
    left_sha = c64._sha256_file(left_path)

    shared_contigs = len(shared_rows)
    shared_bp = sum(shared_bases_by_stage.values())

    report = {
        "stage": "ET1",
        "min_len": min_len,
        "dedup_isoforms": dedup_isoforms,
        "unsupported_code_policy": unsupported_code_policy,
        "fasta_path": str(fasta_path),
        "fasta_bytes": fasta_bytes,
        "fasta_sha256": fasta_sha256,
        "genomes_json_sha256": c64._sha256_file(genomes_path),
        "objects_json_sha256": c64._sha256_file(obj_final_path),
        "objects_staging_json_sha256": c64._sha256_file(obj_staging_path),
        "et1_metadata_tsv_sha256": meta_sha,
        "records_tsv_sha256": rec_sha,
        "repeated_records_tsv_sha256": rep_sha,
        "shared_records_tsv_sha256": shr_sha,
        "left_for_later_tsv_gz_sha256": left_sha,
        "left_for_later_rows": len(left_rows),
        "genomes": len(genomes_out),
        "distinct_contigs": len(written_uids),
        "total_bp": total_bp,
        "shared_contigs": shared_contigs,
        "shared_bp": shared_bp,
        "net_new_contigs": len(written_uids) - shared_contigs,
        "net_new_bp": total_bp - shared_bp,
        "shared_by_stage": dict(shared_by_stage),
        "shared_bases_by_stage": dict(shared_bases_by_stage),
        "genomes_by_effective_code": {str(k): v for k, v in sorted(code_counts.items())},
        "samples_by_bio_code": {str(k): v for k, v in sorted(bio_code_counts.items())},
        "genomes_by_source": dict(sorted(source_counts.items())),
        "counts": dict(sorted(counts.items())),
    }
    (out_meta_dir / "stage_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


# --------------------------------------------------------------------------------------------------
# LGNRNAP1: W-40 ssRNAP super-panel baiting over Logan logan50 (reusing c64_logan_rt_bait.py)
# --------------------------------------------------------------------------------------------------


def load_ssrnap_cutoffs(cutoffs_tsv: Path = DEFAULT_SSRNAP_CUTOFFS_TSV) -> dict[str, dict]:
    """Load `examples/rnap/ssrnap_cutoffs.tsv` into a lookup keyed by both full `namespace:model`
    and bare `model` name."""
    out: dict[str, dict] = {}
    with open(cutoffs_tsv, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key = (row.get("model_key") or "").strip()
            if not key:
                continue
            entry = {
                "model_key": key,
                "namespace": key.split(":", 1)[0] if ":" in key else "palm",
                "model_name": key.split(":", 1)[1] if ":" in key else key,
                "seq_bits": float(row["seq_bits"]),
                "dom_bits": float(row["dom_bits"]),
                "i_evalue_max": float(row["i_evalue_max"]),
            }
            out[key] = entry
            out[entry["model_name"]] = entry
    return out


def _pfam_base_acc(acc: str) -> str:
    return acc.split(".", 1)[0] if "." in acc else acc


def build_ssrnap_bait_hmm(
    out_dir: Path,
    *,
    pfam_hmm: Path = DEFAULT_PFAM_HMM,
    hmms_dir: Path = DEFAULT_W40_HMMS_DIR,
) -> dict:
    """Assemble W-40's 30 target ssRNAP HMMs (25 palm + 5 dpbb_single) and 16 outgroup decoy HMMs
    (1 DNA_pol_A + 15 cellular rpoB/rpoC) into deterministic files for logan50 baiting."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    target_path = out_dir / "bait_ssrnap30.hmm"
    decoy_path = out_dir / "bait_ssrnap_decoys16.hmm"

    target_models: list[dict] = []
    decoy_models: list[dict] = []
    want_pfam_target = set(SSRNAP_PFAM_TARGET_ACCS)
    want_pfam_decoy = set(SSRNAP_PFAM_DECOY_ACCS)
    found_pfam_target: set[str] = set()
    found_pfam_decoy: set[str] = set()

    extra_target_files = (
        ("vogdb_rnap", "palm", hmms_dir / "vogdb_rnap_palm_clean.hmm"),
        ("phrogs_rnap", "palm", hmms_dir / "phrogs_rnap_palm_clean.hmm"),
        ("ssrnap_custom_seeds", "palm", hmms_dir / "seeds_palm.hmm"),
        ("ssrnap_custom_seeds", "dpbb_single", hmms_dir / "seeds_dpbb.hmm"),
    )
    for _, _, p in extra_target_files:
        if not p.is_file():
            raise Refused(f"Missing W-40 HMM file: {p}")
    if not Path(pfam_hmm).is_file():
        raise Refused(f"Missing Pfam HMM file: {pfam_hmm}")

    with open(target_path, "wb") as tfh, open(decoy_path, "wb") as dfh:
        with pyhmmer.plan7.HMMFile(pfam_hmm) as hf:
            for hmm in hf:
                acc = c64._as_str(hmm.accession)
                base = _pfam_base_acc(acc)
                name = c64._as_str(hmm.name)
                if base in want_pfam_target:
                    ns = "dpbb_single" if base in SSRNAP_PFAM_DPBB_ACCS else "palm"
                    hmm.write(tfh, binary=False)
                    found_pfam_target.add(base)
                    target_models.append({
                        "name": name, "accession": acc, "namespace": ns,
                        "length": hmm.M, "pack": "pfam_38.2",
                    })
                elif base in want_pfam_decoy:
                    ns = "decoy_pol1" if base in SSRNAP_PFAM_DECOY_POL1_ACCS else "decoy_cell_rnap"
                    hmm.write(dfh, binary=False)
                    found_pfam_decoy.add(base)
                    decoy_models.append({
                        "name": name, "accession": acc, "namespace": ns,
                        "length": hmm.M, "pack": "pfam_38.2",
                    })

        if found_pfam_target != want_pfam_target:
            raise Refused(f"Missing Pfam target models: {sorted(want_pfam_target - found_pfam_target)}")
        if found_pfam_decoy != want_pfam_decoy:
            raise Refused(f"Missing Pfam decoy models: {sorted(want_pfam_decoy - found_pfam_decoy)}")

        for pack, ns, p in extra_target_files:
            with pyhmmer.plan7.HMMFile(p) as hf:
                for hmm in hf:
                    name = c64._as_str(hmm.name)
                    hmm.write(tfh, binary=False)
                    target_models.append({
                        "name": name, "accession": c64._as_str(hmm.accession),
                        "namespace": ns, "length": hmm.M, "pack": pack,
                    })

    manifest = {
        "target_hmm": str(target_path),
        "target_hmm_sha256": c64._sha256_file(target_path),
        "target_count": len(target_models),
        "target_models": target_models,
        "decoy_hmm": str(decoy_path),
        "decoy_hmm_sha256": c64._sha256_file(decoy_path),
        "decoy_count": len(decoy_models),
        "decoy_models": decoy_models,
    }
    (out_dir / "ssrnap_bait_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def scan_ssrnap_batch(
    batch: list[tuple[dict, bytes]],
    target_hmms: list[pyhmmer.plan7.HMM],
    decoy_hmms: list[pyhmmer.plan7.HMM],
    alphabet: pyhmmer.easel.Alphabet,
    cutoffs: dict[str, dict],
    *,
    cpus: int = 4,
    default_seq_bits: float = 30.0,
    default_dom_bits: float = 20.0,
    evalue: float = 1e-5,
) -> list[dict]:
    """Scan one batch of `logan50` proteins against the W-40 `ssrnap` panel with inline decoy veto."""
    if not batch:
        return []
    seqs: list[pyhmmer.easel.DigitalSequence] = []
    meta_by_name: dict[str, tuple[dict, bytes]] = {}
    for parsed, aa in batch:
        clean = aa.rstrip(b"*")
        if len(clean) < 80:
            continue
        pid = parsed["protein_id"]
        try:
            tb = pyhmmer.easel.TextSequence(
                name=pid.encode("ascii", "replace"),
                sequence=clean.decode("ascii", "replace"),
            )
            seqs.append(tb.digitize(alphabet))
            meta_by_name[pid] = (parsed, clean)
        except Exception:
            continue
    if not seqs:
        return []

    block = pyhmmer.easel.DigitalSequenceBlock(alphabet, seqs)
    candidate_hits: list[dict] = []
    hit_pids: set[str] = set()

    for th in pyhmmer.hmmer.hmmsearch(target_hmms, block, cpus=cpus, E=evalue, domE=evalue):
        model_name = c64._as_str(th.query.name)
        model_acc = c64._as_str(th.query.accession)
        base_acc = _pfam_base_acc(model_acc) if model_acc else ""
        c_entry = cutoffs.get(model_name) or cutoffs.get(model_acc) or cutoffs.get(base_acc)
        if c_entry is not None:
            min_seq = c_entry["seq_bits"]
            min_dom = c_entry["dom_bits"]
            max_ie = c_entry["i_evalue_max"]
            ns = c_entry["namespace"]
        else:
            ga = th.query.cutoffs.gathering if th.query.cutoffs.gathering_available() else None
            min_seq = float(ga[0]) if ga is not None else default_seq_bits
            min_dom = float(ga[1]) if ga is not None else default_dom_bits
            max_ie = evalue
            ns = "dpbb_single" if base_acc in SSRNAP_PFAM_DPBB_ACCS else "palm"

        for hit in th:
            if not hit.included or hit.score < min_seq:
                continue
            best_dom = max(
                (d for d in hit.domains if d.i_evalue <= max_ie),
                key=lambda d: d.score,
                default=None,
            )
            if best_dom is None or best_dom.score < min_dom:
                continue
            pid = c64._as_str(hit.name)
            parsed, _ = meta_by_name[pid]
            candidate_hits.append({
                "protein_id": parsed["protein_id"],
                "cluster_id": parsed["cluster_id"],
                "run_acc": parsed["run_acc"],
                "contig_id": parsed["contig_id"],
                "gene_idx": parsed["gene_idx"],
                "namespace": ns,
                "model": f"{ns}:{model_name}",
                "score": round(float(hit.score), 2),
                "max_dom_score": round(float(best_dom.score), 2),
                "i_evalue": float(best_dom.i_evalue),
                "evalue": float(hit.evalue),
                "vetoed_by_decoy": False,
            })
            hit_pids.add(pid)

    # Inline outgroup decoy veto (DNA_pol_A + 15 cellular rpoB/rpoC profiles)
    if decoy_hmms and hit_pids:
        hit_seqs = []
        for pid in sorted(hit_pids):
            _, clean = meta_by_name[pid]
            tb = pyhmmer.easel.TextSequence(
                name=pid.encode("ascii", "replace"),
                sequence=clean.decode("ascii", "replace"),
            )
            hit_seqs.append(tb.digitize(alphabet))
        hit_block = pyhmmer.easel.DigitalSequenceBlock(alphabet, hit_seqs)
        decoy_best: dict[str, tuple[float, str]] = {}
        for th in pyhmmer.hmmer.hmmsearch(decoy_hmms, hit_block, cpus=min(cpus, 2), E=1e-3):
            dname = c64._as_str(th.query.name)
            for hit in th:
                pid = c64._as_str(hit.name)
                sc = float(hit.score)
                if sc > decoy_best.get(pid, (0.0, ""))[0]:
                    decoy_best[pid] = (sc, dname)
        for rec in candidate_hits:
            d_sc, d_name = decoy_best.get(rec["protein_id"], (0.0, ""))
            if d_sc >= rec["score"]:
                rec["vetoed_by_decoy"] = True
                rec["decoy_score"] = round(d_sc, 2)
                rec["decoy_model"] = d_name

    return [r for r in candidate_hits if not r["vetoed_by_decoy"]]


def stage_lgnrnap1(
    *,
    baited_contigs: list[dict],
    out_meta_dir: Path,
    out_fasta_dir: Path,
    final_fasta_dir: Path | None = None,
    fasta_name: str = "lgnrnap1_baited.fna.gz",
    root: Path | None = None,
    root_stages: Iterable[str] = (),
    owner_inputs: dict[str, Path] | None = None,
    retrieved_at: str = "2026-11-24T00:00:00Z",
    threads: int = 2,
    memory_limit: str = "4GB",
) -> dict:
    """Deduplicate W-40 ssRNAP-baited Logan contigs by `contig_uid` and write `LGNRNAP1` stage inputs."""
    out_meta_dir = Path(out_meta_dir)
    out_fasta_dir = Path(out_fasta_dir)
    out_meta_dir.mkdir(parents=True, exist_ok=True)
    out_fasta_dir.mkdir(parents=True, exist_ok=True)

    earlier_owners = c64._load_earlier_contig_owners(
        root, root_stages, owner_inputs, threads=threads, memory_limit=memory_limit
    )

    fasta_path = out_fasta_dir / fasta_name
    tmp_fasta = fasta_path.with_name(fasta_path.name + ".tmp")

    written_uids: set[bytes] = set()
    genomes_out: list[dict] = []
    records_rows: list[list] = []
    repeated_rows: list[list] = []
    trims_rows: list[list] = []
    uid_len: dict[bytes, int] = {}

    counts = Counter()
    biomes = Counter()
    models_counter = Counter()
    runs_kept: set[str] = set()
    total_bp = 0

    ordered = sorted(baited_contigs, key=lambda r: (r["run_acc"], r["contig_id"]))
    with open(tmp_fasta, "wb") as raw_fh:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_fh, compresslevel=6, mtime=0) as gz_out:
            for rec in ordered:
                seq: bytes = rec["seq"]
                uid = ids.contig_uid(seq)
                cid = rec["contig_id"]
                run_acc = rec["run_acc"]
                is_circ = bool(rec.get("circular", False))
                trimmed_k = int(rec.get("trimmed_k", 0))
                raw_len = int(rec.get("raw_len", len(seq) + trimmed_k))
                bait_models = sorted(set(rec.get("bait_models", [])))
                cluster_ids = sorted(set(rec.get("cluster_ids", [])))
                biome = str(rec.get("biome", "unknown"))

                if uid in written_uids:
                    counts["repeated_contigs"] += 1
                    repeated_rows.append([cid, run_acc, "duplicate_contig_uid", uid.hex()])
                    continue

                written_uids.add(uid)
                uid_len[uid] = len(seq)
                runs_kept.add(run_acc)
                biomes[biome] += 1
                for bm in bait_models:
                    models_counter[bm] += 1

                _write_deterministic_gz_fasta_header(gz_out, uid.hex(), seq)
                total_bp += len(seq)

                if is_circ:
                    counts["circular_contigs"] += 1
                if trimmed_k > 0:
                    counts["dtr_trimmed"] += 1
                    trims_rows.append([run_acc, cid, raw_len, trimmed_k, len(seq)])

                genomes_out.append({
                    "contig_uids": [uid.hex()],
                    "genetic_code": 11,
                    "kind": "prokaryote",
                    "circular": [uid.hex()] if is_circ else [],
                    "assembly_quality": "contig",
                })
                records_rows.append([
                    cid, run_acc, uid.hex(), len(seq), int(is_circ),
                    ",".join(bait_models), ",".join(cluster_ids), biome,
                ])

    if not genomes_out:
        tmp_fasta.unlink(missing_ok=True)
        raise Refused("no qualifying baited Logan contigs staged for LGNRNAP1")

    fasta_bytes = tmp_fasta.stat().st_size
    if str(fasta_path.resolve()).startswith("/stor/") and fasta_bytes < MIN_STOR_FILE_BYTES:
        tmp_fasta.unlink(missing_ok=True)
        raise Refused(
            f"Refusing to write {fasta_bytes} B (< 64 MB) to /stor path {fasta_path} (AGENTS.md §7.3); "
            f"stage under /tmp/if2453/s2/inputs/LGNRNAP1/ instead."
        )
    tmp_fasta.replace(fasta_path)
    fasta_sha256 = c64._sha256_file(fasta_path)

    shared_rows: list[list] = []
    shared_by_stage = Counter()
    shared_bases_by_stage = Counter()
    for uid in sorted(written_uids):
        owner = earlier_owners.get(uid)
        if owner is not None:
            ln = uid_len[uid]
            shared_rows.append([uid.hex(), owner, ln])
            shared_by_stage[owner] += 1
            shared_bases_by_stage[owner] += ln

    genomes_path = out_meta_dir / "genomes.json"
    genomes_path.write_text(json.dumps(genomes_out, separators=(",", ":")) + "\n", encoding="utf-8")

    obj_staging = [{
        "path": str(fasta_path),
        "source_name": "logan_contigs_baited",
        "source_release": "logan-v1-c1.0-lgnrnap1",
        "license": "INSDC / AWS Open Data",
        "redistributable": True,
        "retrieved_at": retrieved_at,
    }]
    obj_staging_path = out_meta_dir / "objects_staging.json"
    obj_staging_path.write_text(json.dumps(obj_staging, indent=2) + "\n", encoding="utf-8")

    final_path = (final_fasta_dir / fasta_name) if final_fasta_dir is not None else fasta_path
    obj_final = [{**obj_staging[0], "path": str(final_path)}]
    obj_final_path = out_meta_dir / "objects.json"
    obj_final_path.write_text(json.dumps(obj_final, indent=2) + "\n", encoding="utf-8")

    def _write_tsv(path: Path, header: list[str], rows: list[list]) -> str:
        with open(path, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(header)
            w.writerows(rows)
        return c64._sha256_file(path)

    rec_sha = _write_tsv(
        out_meta_dir / "records.tsv",
        ["contig_id", "run_acc", "contig_uid", "length", "circular", "bait_models", "cluster_ids", "biome"],
        records_rows,
    )
    rep_sha = _write_tsv(
        out_meta_dir / "repeated_records.tsv",
        ["contig_id", "run_acc", "reason", "contig_uid"],
        repeated_rows,
    )
    shr_sha = _write_tsv(
        out_meta_dir / "shared_records.tsv",
        ["contig_uid", "owner_stage", "length"],
        shared_rows,
    )
    trm_sha = _write_tsv(
        out_meta_dir / "trims.tsv",
        ["run_acc", "contig_id", "raw_len", "trimmed_k", "final_len"],
        trims_rows,
    )

    shared_contigs = len(shared_rows)
    shared_bp = sum(shared_bases_by_stage.values())

    report = {
        "stage": "LGNRNAP1",
        "fasta_path": str(fasta_path),
        "fasta_bytes": fasta_bytes,
        "fasta_sha256": fasta_sha256,
        "genomes_json_sha256": c64._sha256_file(genomes_path),
        "objects_json_sha256": c64._sha256_file(obj_final_path),
        "objects_staging_json_sha256": c64._sha256_file(obj_staging_path),
        "records_tsv_sha256": rec_sha,
        "repeated_records_tsv_sha256": rep_sha,
        "shared_records_tsv_sha256": shr_sha,
        "trims_tsv_sha256": trm_sha,
        "distinct_contigs": len(written_uids),
        "distinct_runs": len(runs_kept),
        "total_bp": total_bp,
        "repeated_contigs": counts["repeated_contigs"],
        "circular_contigs": counts["circular_contigs"],
        "dtr_trimmed": counts["dtr_trimmed"],
        "shared_contigs": shared_contigs,
        "shared_bp": shared_bp,
        "net_new_contigs": len(written_uids) - shared_contigs,
        "net_new_bp": total_bp - shared_bp,
        "shared_by_stage": dict(shared_by_stage),
        "shared_bases_by_stage": dict(shared_bases_by_stage),
        "contigs_by_model": dict(sorted(models_counter.items())),
        "top_biomes": biomes.most_common(30),
    }
    (out_meta_dir / "stage_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


# --------------------------------------------------------------------------------------------------
# CLI subcommands (pod execution)
# --------------------------------------------------------------------------------------------------


def _download_with_md5(url: str, dest: Path, expected_md5: str | None = None) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    h_md5 = hashlib.md5()
    h_sha = hashlib.sha256()
    n_bytes = 0
    with urllib.request.urlopen(req, timeout=300) as resp, open(tmp, "wb") as out:
        while chunk := resp.read(1 << 20):
            out.write(chunk)
            h_md5.update(chunk)
            h_sha.update(chunk)
            n_bytes += len(chunk)
    md5_hex = h_md5.hexdigest()
    sha_hex = h_sha.hexdigest()
    if expected_md5 and md5_hex != expected_md5:
        tmp.unlink(missing_ok=True)
        raise Refused(f"MD5 mismatch for {dest.name}: expected {expected_md5}, got {md5_hex}")
    tmp.replace(dest)
    return {"path": str(dest), "url": url, "bytes": n_bytes, "md5": md5_hex, "sha256": sha_hex}


def cmd_fetch_et1(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fetched: list[dict] = []

    if not args.skip_mmetsp:
        mmetsp_dir = out_dir / "mmetsp"
        mmetsp_tar = mmetsp_dir / MMETSP_FASTA_TAR_NAME
        if not mmetsp_tar.is_file():
            fetched.append(
                _download_with_md5(MMETSP_FASTA_TAR_URL, mmetsp_tar, expected_md5=MMETSP_FASTA_TAR_MD5)
            )
        else:
            fetched.append({
                "path": str(mmetsp_tar),
                "url": MMETSP_FASTA_TAR_URL,
                "bytes": mmetsp_tar.stat().st_size,
                "sha256": c64._sha256_file(mmetsp_tar),
                "cached": True,
            })

    if not args.skip_onekp:
        onekp_dir = out_dir / "onekp"
        onekp_dir.mkdir(parents=True, exist_ok=True)
        # Write verified upstream manifest for batch retrieval of 1KP GigaDB 100910 / CyVerse assemblies
        manifest = {
            "source": args.onekp_source,
            "gigadb_doi": ONEKP_GIGADB_DOI,
            "gigadb_ftp_base": ONEKP_GIGADB_FTP_BASE,
            "cyverse_base": ONEKP_CYVERSE_BASE,
            "license": "CC0 1.0 (GigaDB 10.5524/100910)",
            "redistributable": True,
        }
        (onekp_dir / "onekp_upstream.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        fetched.append({"path": str(onekp_dir / "onekp_upstream.json"), "manifest": True})

    lock_path = out_dir / "et1_fetch_lock.json"
    lock_path.write_text(json.dumps({"fetched": fetched}, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"lock": str(lock_path), "items": len(fetched)}, indent=2))


def cmd_stage_et1(args: argparse.Namespace) -> None:
    samples: list[dict] = []
    with open(args.samples_tsv, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            samples.append({
                "sample_id": row["sample_id"],
                "source_name": row["source_name"],
                "source_release": row.get("source_release", "et1"),
                "license": row.get("license", "CC BY 4.0"),
                "organism": row.get("organism", ""),
                "taxonomy": row.get("taxonomy", ""),
                "fasta_path": Path(row["fasta_path"]),
            })

    owner_inputs = {}
    for item in args.owner_inputs or []:
        k, v = item.split("=", 1)
        owner_inputs[k] = Path(v)

    rep = stage_et1(
        samples=samples,
        out_meta_dir=Path(args.out_meta_dir),
        out_fasta_dir=Path(args.out_fasta_dir),
        final_fasta_dir=Path(args.final_fasta_dir) if args.final_fasta_dir else None,
        fasta_name=args.fasta_name,
        min_len=args.min_len,
        dedup_isoforms=args.dedup_isoforms,
        unsupported_code_policy=args.unsupported_code_policy,
        root=Path(args.root) if args.root else None,
        root_stages=args.root_stages.split(",") if args.root_stages else (),
        owner_inputs=owner_inputs,
        retrieved_at=args.retrieved_at,
        threads=args.threads,
        memory_limit=args.memory_limit,
    )
    print(json.dumps(rep, indent=2))


def cmd_build_ssrnap_bait(args: argparse.Namespace) -> None:
    manifest = build_ssrnap_bait_hmm(
        Path(args.out_dir),
        pfam_hmm=Path(args.pfam_hmm),
        hmms_dir=Path(args.hmms_dir),
    )
    print(json.dumps(manifest, indent=2))


def cmd_bait_lgnrnap1(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    hits_path = out_dir / "hits.jsonl"
    ckpt_path = out_dir / "checkpoint.json"

    cutoffs = load_ssrnap_cutoffs(Path(args.cutoffs_tsv))
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.plan7.HMMFile(args.bait_hmm) as hf:
        target_hmms = list(hf)
    decoy_hmms: list[pyhmmer.plan7.HMM] = []
    if args.decoy_hmm and Path(args.decoy_hmm).is_file():
        with pyhmmer.plan7.HMMFile(args.decoy_hmm) as hf:
            decoy_hmms = list(hf)

    url = f"{c64.LOGAN_BASE}/{c64.LOGAN50_FASTA_KEY}"
    headers = {"User-Agent": UA}
    if args.start_byte > 0:
        end_byte = (
            args.start_byte + args.max_compressed_bytes - 1
            if args.max_compressed_bytes > 0
            else c64.LOGAN50_FASTA_TOTAL_BYTES - 1
        )
        headers["Range"] = f"bytes={args.start_byte}-{end_byte}"
    elif args.max_compressed_bytes > 0:
        headers["Range"] = f"bytes=0-{args.max_compressed_bytes - 1}"

    t0 = time.monotonic()
    c0 = time.process_time()
    comp_bytes = [0]
    stop_flag = threading.Event()

    proc = subprocess.Popen(
        ["zstd", "-dc", "--long=31"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=1 << 20,
    )

    def _feeder() -> None:
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=120) as resp:
                while not stop_flag.is_set():
                    chunk = resp.read(1 << 18)
                    if not chunk:
                        break
                    comp_bytes[0] += len(chunk)
                    if proc.stdin is not None:
                        try:
                            proc.stdin.write(chunk)
                        except BrokenPipeError:
                            break
                    if args.max_compressed_bytes > 0 and comp_bytes[0] >= args.max_compressed_bytes:
                        break
                    if args.max_wall_s > 0 and (time.monotonic() - t0) >= args.max_wall_s:
                        break
        finally:
            if proc.stdin is not None:
                try:
                    proc.stdin.close()
                except Exception:
                    pass

    feeder = threading.Thread(target=_feeder, daemon=True)
    feeder.start()

    proteins_scanned = 0
    hits_found = 0
    hits_by_model = Counter()
    last_protein_id = ""
    batch: list[tuple[dict, bytes]] = []

    with open(hits_path, "w", encoding="utf-8") as hfh:
        assert proc.stdout is not None
        for hdr, seq in c64._iter_fasta_bytes_stream(proc.stdout):
            try:
                parsed = c64.parse_logan_protein_header(hdr)
            except ValueError:
                continue
            last_protein_id = parsed["protein_id"]
            batch.append((parsed, seq))
            proteins_scanned += 1
            if len(batch) >= args.batch_size:
                hits = scan_ssrnap_batch(
                    batch, target_hmms, decoy_hmms, alphabet, cutoffs,
                    cpus=args.cpus, evalue=args.evalue,
                )
                for h in hits:
                    hfh.write(json.dumps(h) + "\n")
                    hits_found += 1
                    hits_by_model[h["model"]] += 1
                hfh.flush()
                batch.clear()

                wall_s = time.monotonic() - t0
                ckpt = {
                    "url": url,
                    "start_byte": args.start_byte,
                    "compressed_bytes_read": comp_bytes[0],
                    "total_archive_bytes": c64.LOGAN50_FASTA_TOTAL_BYTES,
                    "pct_archive_scanned": round(100.0 * comp_bytes[0] / c64.LOGAN50_FASTA_TOTAL_BYTES, 4),
                    "proteins_scanned": proteins_scanned,
                    "hits_found": hits_found,
                    "hits_by_model": dict(hits_by_model),
                    "last_protein_id": last_protein_id,
                    "wall_s": round(wall_s, 2),
                    "cpus": args.cpus,
                    "status": "running",
                }
                ckpt_path.write_text(json.dumps(ckpt, indent=2) + "\n", encoding="utf-8")
                if args.max_wall_s > 0 and wall_s >= args.max_wall_s:
                    stop_flag.set()
                    break

        if batch:
            hits = scan_ssrnap_batch(
                batch, target_hmms, decoy_hmms, alphabet, cutoffs,
                cpus=args.cpus, evalue=args.evalue,
            )
            for h in hits:
                hfh.write(json.dumps(h) + "\n")
                hits_found += 1
                hits_by_model[h["model"]] += 1
            hfh.flush()

    stop_flag.set()
    try:
        proc.terminate()
    except Exception:
        pass
    proc.wait(timeout=15)
    feeder.join(timeout=15)

    wall_s = time.monotonic() - t0
    cpu_s = time.process_time() - c0
    final_ckpt = {
        "url": url,
        "start_byte": args.start_byte,
        "compressed_bytes_read": comp_bytes[0],
        "total_archive_bytes": c64.LOGAN50_FASTA_TOTAL_BYTES,
        "pct_archive_scanned": round(100.0 * comp_bytes[0] / c64.LOGAN50_FASTA_TOTAL_BYTES, 4),
        "proteins_scanned": proteins_scanned,
        "hits_found": hits_found,
        "hits_by_model": dict(hits_by_model),
        "last_protein_id": last_protein_id,
        "wall_s": round(wall_s, 2),
        "cpu_s": round(cpu_s, 2),
        "cpus": args.cpus,
        "status": "complete",
    }
    ckpt_path.write_text(json.dumps(final_ckpt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(final_ckpt, indent=2))

    if not args.skip_resolve:
        res_args = argparse.Namespace(
            hits_jsonl=str(hits_path),
            out_dir=str(out_dir / "resolved"),
            index_path=args.index_path,
            max_members_per_cluster=args.max_members_per_cluster,
            max_frames=args.max_frames,
            max_runs=args.max_runs,
        )
        c64.cmd_resolve(res_args)


def cmd_stage_lgnrnap1(args: argparse.Namespace) -> None:
    resolved = json.loads(Path(args.resolved_json).read_text(encoding="utf-8"))
    runs_map: dict[str, dict[str, dict]] = resolved["runs"]
    run_list = sorted(runs_map.keys())

    biomes_map: dict[str, dict] = {}
    if not args.skip_ena and run_list:
        biomes_map = c64._fetch_ena_biomes(run_list)

    baited_contigs: list[dict] = []
    fetch_stats = Counter()
    failed_runs: list[str] = []
    deferred_large_runs: list[list] = []
    out_meta = Path(args.out_meta_dir)
    out_meta.mkdir(parents=True, exist_ok=True)
    max_run_bytes = int(args.max_run_bytes)

    def _fetch_one_run(run_acc: str) -> tuple[str, list[dict], dict, str | None, int]:
        wanted_map = runs_map[run_acc]
        wanted_set = set(wanted_map.keys())
        url = f"{c64.LOGAN_BASE}/c/{run_acc}/{run_acc}.contigs.fa.zst"
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        try:
            resp = urllib.request.urlopen(req, timeout=120)
            cl_hdr = resp.headers.get("Content-Length")
            cl_bytes = int(cl_hdr) if cl_hdr and cl_hdr.isdigit() else 0
            if max_run_bytes > 0 and cl_bytes > max_run_bytes:
                resp.close()
                return run_acc, [], {}, "deferred_large_run", cl_bytes

            proc = subprocess.Popen(
                ["zstd", "-dc"],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
            )
            err_box: list[Exception] = []

            def _feed(r=resp, p=proc, eb=err_box) -> None:
                try:
                    with r:
                        while chunk := r.read(1 << 18):
                            if p.stdin is not None:
                                p.stdin.write(chunk)
                except Exception as exc:
                    eb.append(exc)
                finally:
                    if p.stdin is not None:
                        try:
                            p.stdin.close()
                        except Exception:
                            pass

            t = threading.Thread(target=_feed, daemon=True)
            t.start()
            assert proc.stdout is not None
            kept, st = c64.extract_baited_contigs_from_fasta(
                c64._iter_wanted_fasta_bytes_stream(proc.stdout, wanted_set),
                wanted_set,
                min_len=args.min_len,
            )
            proc.wait(timeout=60)
            t.join(timeout=60)
            if err_box:
                raise err_box[0]

            ena_row = biomes_map.get(run_acc, {})
            biome = (ena_row.get("scientific_name") or "unknown").strip() or "unknown"
            lib_src = (ena_row.get("library_source") or "").strip()
            if lib_src and lib_src.upper() == "METATRANSCRIPTOMIC":
                biome = f"{biome} [METATRANSCRIPTOMIC]"

            for rec in kept:
                info = wanted_map.get(rec["contig_id"], {})
                rec["bait_models"] = info.get("models", [])
                rec["cluster_ids"] = info.get("clusters", [])
                rec["biome"] = biome
            return run_acc, kept, st, None, cl_bytes
        except Exception as exc:
            return run_acc, [], {}, str(exc), 0

    workers = max(1, min(int(args.workers), 4))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_fetch_one_run, r) for r in run_list]
        for fut in as_completed(futs):
            run_acc, kept, st, err, cl_bytes = fut.result()
            if err is None:
                fetch_stats["runs_ok"] += 1
                for k, v in st.items():
                    fetch_stats[k] += v
                baited_contigs.extend(kept)
            elif err == "deferred_large_run":
                fetch_stats["runs_deferred_large"] += 1
                deferred_large_runs.append([run_acc, cl_bytes, len(runs_map[run_acc]), "over_max_run_bytes"])
            else:
                failed_runs.append(run_acc)
                fetch_stats["runs_failed"] += 1

    if deferred_large_runs:
        deferred_large_runs.sort(key=lambda r: r[0])
        with open(out_meta / "LGNRNAP1_deferred_large_runs.tsv", "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["run_acc", "compressed_bytes", "wanted_contigs", "reason"])
            w.writerows(deferred_large_runs)

    owner_inputs = {}
    for item in args.owner_inputs or []:
        k, v = item.split("=", 1)
        owner_inputs[k] = Path(v)

    rep = stage_lgnrnap1(
        baited_contigs=baited_contigs,
        out_meta_dir=out_meta,
        out_fasta_dir=Path(args.out_fasta_dir),
        final_fasta_dir=Path(args.final_fasta_dir) if args.final_fasta_dir else None,
        fasta_name=args.fasta_name,
        root=Path(args.root) if args.root else None,
        root_stages=args.root_stages.split(",") if args.root_stages else (),
        owner_inputs=owner_inputs,
        retrieved_at=args.retrieved_at,
        threads=args.threads,
        memory_limit=args.memory_limit,
    )
    rep["fetch_stats"] = dict(fetch_stats)
    rep["failed_runs_count"] = len(failed_runs)
    (out_meta / "stage_report.json").write_text(json.dumps(rep, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(rep, indent=2))


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(
        description="Card C-68: Spliced protist & plant transcriptomes (ET1) and Logan ssRNAP baiting (LGNRNAP1)"
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    p_fe = sub.add_parser("fetch-et1")
    p_fe.add_argument("out_dir")
    p_fe.add_argument("--skip-mmetsp", action="store_true")
    p_fe.add_argument("--skip-onekp", action="store_true")
    p_fe.add_argument("--onekp-source", choices=("gigadb", "cyverse"), default="gigadb")
    p_fe.set_defaults(func=cmd_fetch_et1)

    p_se = sub.add_parser("stage-et1")
    p_se.add_argument("samples_tsv")
    p_se.add_argument("out_meta_dir")
    p_se.add_argument("out_fasta_dir")
    p_se.add_argument("--final-fasta-dir", default=None)
    p_se.add_argument("--fasta-name", default="et1_transcriptomes.fna.gz")
    p_se.add_argument("--min-len", type=int, default=ET1_MIN_LEN)
    p_se.add_argument("--dedup-isoforms", action="store_true")
    p_se.add_argument("--unsupported-code-policy", choices=UNSUPPORTED_CODE_POLICIES, default="defer")
    p_se.add_argument("--root", default=None)
    p_se.add_argument("--root-stages", default="")
    p_se.add_argument("--owner-inputs", action="append", default=[])
    p_se.add_argument("--retrieved-at", default="2026-11-24T00:00:00Z")
    p_se.add_argument("--threads", type=int, default=2)
    p_se.add_argument("--memory-limit", default="4GB")
    p_se.set_defaults(func=cmd_stage_et1)

    p_bb = sub.add_parser("build-ssrnap-bait")
    p_bb.add_argument("out_dir")
    p_bb.add_argument("--pfam-hmm", default=str(DEFAULT_PFAM_HMM))
    p_bb.add_argument("--hmms-dir", default=str(DEFAULT_W40_HMMS_DIR))
    p_bb.set_defaults(func=cmd_build_ssrnap_bait)

    p_bl = sub.add_parser("bait-lgnrnap1")
    p_bl.add_argument("bait_hmm")
    p_bl.add_argument("out_dir")
    p_bl.add_argument("--decoy-hmm", default=None)
    p_bl.add_argument("--cutoffs-tsv", default=str(DEFAULT_SSRNAP_CUTOFFS_TSV))
    p_bl.add_argument("--cpus", type=int, default=6)
    p_bl.add_argument("--batch-size", type=int, default=50_000)
    p_bl.add_argument("--max-compressed-bytes", type=int, default=0)
    p_bl.add_argument("--max-wall-s", type=float, default=0.0)
    p_bl.add_argument("--evalue", type=float, default=1e-5)
    p_bl.add_argument("--start-byte", type=int, default=0)
    p_bl.add_argument("--index-path", default=None)
    p_bl.add_argument("--max-members-per-cluster", type=int, default=50)
    p_bl.add_argument("--max-frames", type=int, default=200)
    p_bl.add_argument("--max-runs", type=int, default=500)
    p_bl.add_argument("--skip-resolve", action="store_true")
    p_bl.set_defaults(func=cmd_bait_lgnrnap1)

    p_sl = sub.add_parser("stage-lgnrnap1")
    p_sl.add_argument("resolved_json")
    p_sl.add_argument("out_meta_dir")
    p_sl.add_argument("out_fasta_dir")
    p_sl.add_argument("--final-fasta-dir", default=None)
    p_sl.add_argument("--fasta-name", default="lgnrnap1_baited.fna.gz")
    p_sl.add_argument("--min-len", type=int, default=LGNRNAP1_MIN_LEN)
    p_sl.add_argument("--root", default=None)
    p_sl.add_argument("--root-stages", default="")
    p_sl.add_argument("--owner-inputs", action="append", default=[])
    p_sl.add_argument("--retrieved-at", default="2026-11-24T00:00:00Z")
    p_sl.add_argument("--skip-ena", action="store_true")
    p_sl.add_argument("--max-run-bytes", type=int, default=100_000_000)
    p_sl.add_argument("--workers", type=int, default=4)
    p_sl.add_argument("--threads", type=int, default=2)
    p_sl.add_argument("--memory-limit", default="4GB")
    p_sl.set_defaults(func=cmd_stage_lgnrnap1)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
