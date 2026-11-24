"""Card C-58 (DEC-69): Survey and stage single-cell SAGs (GORG -> SAG1), genome-reduced
endosymbionts (SYM1), post-R226 prokaryotes, and niche-habitat MAG compendia (NX1..).

Subcommands:
    fetch-tables  OUT_DIR
    diff-ncbi     TABLES_DIR GTDB_TAR GLOBDB_TAX GLOBDB_DICT_TAR OUT_DIR
    fetch-sym1    SYM1_CANDIDATES_TSV OUT_DIR [--max-genomes N]
    stage-sym1    SYM1_MANIFEST_TSV OUT_ROOT [--root ROOT] [--owner-inputs NAME=DIR ...]
    fetch-gorg    OUT_DIR
    stage-sag1    GORG_FASTA_GZ OUT_ROOT [--fasta-dir DIR] [--final-dir DIR]
                  [--gorg-ena-tsv TSV] [--gorg-dark-ena-tsv TSV]
                  [--gtdb-tar TAR] [--globdb-tax GZ] [--iuig-zip ZIP]
                  [--root ROOT] [--owner-inputs NAME=DIR ...] [--workers 2]
"""
from __future__ import annotations

import argparse
import array
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.error
import urllib.request
import zipfile
from collections import Counter, defaultdict
from collections.abc import Iterable
from concurrent.futures import Future, ProcessPoolExecutor
from pathlib import Path

import duckdb
import pyrodigal

sys.path.insert(0, str(Path(__file__).resolve().parent))

import make_inputs  # noqa: E402
from biobox import fasta, ids, stage as bb_stage  # noqa: E402
from biobox.l0seq import store  # noqa: E402
from biobox.l1 import genecall  # noqa: E402

Refused = make_inputs.Refused

NCBI_SUMMARY_URLS = [
    ("refseq", "bacteria", "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"),
    ("refseq", "archaea", "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt"),
    ("genbank", "bacteria", "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"),
    ("genbank", "archaea", "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary.txt"),
]

ENA_GORG_URLS = [
    (
        "PRJEB33281",
        "https://www.ebi.ac.uk/ena/portal/api/search?result=assembly&query=study_accession=%22PRJEB33281%22"
        "&fields=accession,assembly_name,base_count,scientific_name,tax_id,wgs_set,version&limit=0&format=tsv",
    ),
    (
        "PRJNA1259050",
        "https://www.ebi.ac.uk/ena/portal/api/search?result=assembly&query=study_accession=%22PRJNA1259050%22"
        "&fields=accession,assembly_name,base_count,scientific_name,tax_id,wgs_set,version&limit=0&format=tsv",
    ),
]

GORG_OSF_FILES = [
    (
        "GORG_v1.fasta.gz",
        "https://osf.io/download/mqcj4/",
        1441099680,
        "b5c6ac429ed7207b3284b3593cda9a57",
        "ef097ad1dbc4a8da84cc53a52ba828dea230a82445275732085ea2052bd47910",
    ),
    (
        "GORG_v1.tsv.gz",
        "https://osf.io/download/qja57/",
        228096346,
        "d7986fc77bf2b8ed9b6eef3348fd7b56",
        "299e7189b365203cf245d7b22e99e59bb6a4087d5c795177c41f6da77c22a5d5",
    ),
]

SYM1_SOURCE_BASE = {
    "source_name": "ncbi_c58_symbionts",
    "source_release": "2026-11-23",
    "license": "INSDC (public)",
    "redistributable": True,
}

SAG1_SOURCE_BASE = {
    "source_name": "gorg_tropics",
    "source_release": "v1",
    "license": "CC BY-NC 4.0 (GORG-Tropics OSF pcwj9, Pachiadaki et al. 2019, doi:10.1016/j.cell.2019.11.017)",
    "redistributable": False,
}

_ASSEMBLY_LEVEL_RANK = {
    "Complete Genome": 0,
    "Chromosome": 1,
    "Scaffold": 2,
    "Contig": 3,
}

_ASSEMBLY_QUALITY_MAP = {
    "Complete Genome": "complete",
    "Chromosome": "chromosome",
    "Scaffold": "scaffold",
    "Contig": "contig",
}

_SYMBIONT_RE = re.compile(
    r"(endosymbiont|symbiont|carsonella|hodgkinia|nasuia|tremblaya|zinderia|sulcia|buchnera|"
    r"portiera|uzinura|moranella|vidania|karelsulcia|blochmannia|wigglesworthia|baumannia|"
    r"riesia|hamiltonella|regiella|kinetoplastibacterium|endolissoclinum|ruthia|vesicomyosocius|"
    r"nanoarchaeum|mycoplasma|ureaplasma|phytoplasma|spiroplasma|mesoplasma|entomoplasma|"
    r"wolbachia|rickettsia|chlamydia|ehrlichia|anaplasma|bartonella|tropheryma|"
    r"profftella|stammera|ishikawaella|rosenbergiella|nardonella|westella|diaphorina|"
    r"mikella|doolittlea|hoaglandella|gullanella|walczuchella|blattabacterium|sodalis)",
    re.IGNORECASE,
)

_CODE25_TAX_RE = re.compile(
    r"(gracilibacteria|absconditabacteria|absconditabacterales|bd1-5|sr1|jaedam01|patescibacteria)",
    re.IGNORECASE,
)

_CODE4_TAX_RE = re.compile(
    r"(hodgkinia|nasuia|zinderia|mycoplasma|ureaplasma|spiroplasma|mesoplasma|entomoplasma|"
    r"mycoplasmoides|mycoplasmopsis|metamycoplasma|malacoplasma|acholeplasma|phytoplasma)",
    re.IGNORECASE,
)

_GORG_SAG_RE = re.compile(r"^([A-Z]{2,4}-\d+-[A-Z]\d+)(?:[_\s|:/.-]|$)")

GC_HEAD = [
    "genome", "stage", "organism_name", "genome_size", "contig_number", "assembly_quality",
    "genetic_code", "reason", "tested_alt",
    "cd_11", "mean_orf_11", "n_genes_11",
    "cd_25", "mean_orf_25", "n_genes_25",
    "cd_4", "mean_orf_4", "n_genes_4",
    "tga_per_kb",
]


def _sql_lit(path: Path | str) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def _sql_list(files: Iterable[Path]) -> str:
    return "[" + ", ".join(_sql_lit(f) for f in files) + "]"


def _fmt_opt(val: float | int | None) -> str:
    return "" if val is None else str(val)


def strip_acc_version(acc: str) -> str:
    acc = acc.strip()
    if acc.startswith(("RS_", "GB_")):
        acc = acc[3:]
    if "." in acc:
        return acc.split(".", 1)[0]
    return acc


def parse_gorg_sag_id(header: str) -> str:
    """Extract the SAG cell identifier (e.g. 'AG-313-A02') from a GORG FASTA header."""
    tok = header.strip().split()[0]
    m = _GORG_SAG_RE.match(tok)
    if m:
        return m.group(1)
    if "_NODE_" in tok:
        return tok.split("_NODE_", 1)[0]
    if "_contig" in tok:
        return tok.split("_contig", 1)[0]
    return tok.split("_", 1)[0]


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_ncbi_assembly_summary(path: Path):
    """Stream rows from an NCBI assembly_summary.txt(.gz) file as dicts."""
    with _open_text(path) as fh:
        cols: list[str] | None = None
        for line in fh:
            if line.startswith("#"):
                if line.startswith("#assembly_accession\t") or line.startswith("# assembly_accession\t"):
                    cols = line.lstrip("#").strip().split("\t")
                continue
            if cols is None:
                raise Refused(f"{path}: missing #assembly_accession header")
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < len(cols):
                parts.extend([""] * (len(cols) - len(parts)))
            yield dict(zip(cols, parts))


def load_gtdb_base_accessions(gtdb_tar: Path) -> set[str]:
    """Load versionless GCA_/GCF_ accessions from GTDB R226 bac120 + ar53 metadata tar."""
    accs: set[str] = set()
    with tarfile.open(gtdb_tar, "r:*") as tf:
        for m in tf:
            if not m.isfile() or not m.name.endswith(".tsv.gz"):
                continue
            ef = tf.extractfile(m)
            assert ef is not None
            with gzip.open(ef, "rt", encoding="utf-8", errors="replace") as fh:
                hdr = fh.readline().rstrip("\r\n").split("\t")
                idx_acc = hdr.index("accession")
                idx_gb = hdr.index("ncbi_genbank_assembly_accession") if "ncbi_genbank_assembly_accession" in hdr else -1
                for line in fh:
                    parts = line.rstrip("\r\n").split("\t")
                    if idx_acc < len(parts) and parts[idx_acc]:
                        accs.add(strip_acc_version(parts[idx_acc]))
                    if idx_gb >= 0 and idx_gb < len(parts) and parts[idx_gb] and parts[idx_gb] != "none":
                        accs.add(strip_acc_version(parts[idx_gb]))
    return accs


def load_globdb_base_accessions(globdb_tax_gz: Path, globdb_dict_tar: Path | None = None) -> set[str]:
    """Load versionless GCA_/GCF_ accessions from GlobDB r226 tax_plus_stats and dictionaries."""
    accs: set[str] = set()
    with gzip.open(globdb_tax_gz, "rt", encoding="utf-8", errors="replace") as fh:
        hdr = fh.readline().rstrip("\r\n").split("\t")
        idx_id = 0
        for line in fh:
            gid = line.split("\t", 1)[idx_id].strip()
            if gid.startswith(("GCA_", "GCF_", "RS_GCF_", "GB_GCA_")):
                accs.add(strip_acc_version(gid))
    if globdb_dict_tar is not None and globdb_dict_tar.is_file():
        with tarfile.open(globdb_dict_tar, "r:*") as tf:
            for m in tf:
                if not m.isfile():
                    continue
                ef = tf.extractfile(m)
                if ef is None:
                    continue
                for raw in ef:
                    line = raw.decode("utf-8", "replace")
                    for tok in line.split():
                        if tok.startswith(("GCA_", "GCF_")):
                            accs.add(strip_acc_version(tok))
    return accs


def _safe_int(val: str, default: int = 0) -> int:
    try:
        return int(val)
    except (ValueError, TypeError):
        return default


def _asm_rank_key(row: dict) -> tuple:
    lvl = _ASSEMBLY_LEVEL_RANK.get(row["assembly_level"], 9)
    db_rank = 0 if row["db"] == "refseq" else 1
    contigs = row["contig_count"] if row["contig_count"] > 0 else 999_999_999
    # Prefer Complete/Chromosome, RefSeq, fewest contigs, largest genome_size, lexicographically smallest accession
    return (lvl, db_rank, contigs, -row["genome_size"], row["assembly_accession"])


def diff_ncbi_summaries(
    summaries: list[tuple[str, str, Path]],
    *,
    gtdb_base_accs: set[str],
    globdb_base_accs: set[str],
    sym_max_bp: int = 1_000_000,
    sym_min_bp: int = 50_000,
) -> dict:
    """Diff NCBI RefSeq+GenBank prokaryotic summaries against GTDB R226 and GlobDB r226.

    Separates:
      - `sym1_candidates`: < 1 Mbp genome-reduced endosymbionts / obligate symbionts (1 best per species_taxid)
      - `post_r226_species_reps`: remaining absent prokaryotic species reps (1 best per species_taxid)
    """
    held_accs = gtdb_base_accs | globdb_base_accs
    counts: Counter = Counter()
    seen_bases: set[str] = set()
    sym_by_sp: dict[str, dict] = {}
    post_by_sp: dict[str, dict] = {}

    # Process RefSeq before GenBank so RefSeq wins paired GCF/GCA deduplication
    ordered = sorted(summaries, key=lambda x: (0 if x[0] == "refseq" else 1, x[1]))

    for db, domain, path in ordered:
        for r in parse_ncbi_assembly_summary(path):
            counts["total_rows"] += 1
            counts[f"rows_{db}_{domain}"] += 1
            if r.get("version_status") != "latest":
                counts["suppressed_skipped"] += 1
                continue
            if r.get("genome_rep") != "Full":
                counts["partial_skipped"] += 1
                continue
            acc = r["assembly_accession"].strip()
            base = strip_acc_version(acc)
            paired = r.get("gbrs_paired_asm", "").strip()
            paired_base = strip_acc_version(paired) if paired and paired != "na" else ""

            if base in seen_bases or (paired_base and paired_base in seen_bases):
                counts["paired_refseq_genbank_dedup"] += 1
                continue
            seen_bases.add(base)
            if paired_base:
                seen_bases.add(paired_base)

            if base in held_accs or (paired_base and paired_base in held_accs):
                counts["in_gtdb_or_globdb"] += 1
                continue

            gsize = _safe_int(r.get("genome_size", "0"))
            ccnt = _safe_int(r.get("contig_count", "0"))
            sp_taxid = r.get("species_taxid", "").strip() or r.get("taxid", "").strip() or acc
            org = r.get("organism_name", "").strip()
            lvl = r.get("assembly_level", "").strip()
            rel_date = r.get("seq_rel_date", "").strip()
            excl = r.get("excluded_from_refseq", "").strip()
            ftp = r.get("ftp_path", "").strip()

            rec = {
                "assembly_accession": acc,
                "db": db,
                "domain": domain,
                "taxid": r.get("taxid", "").strip(),
                "species_taxid": sp_taxid,
                "organism_name": org,
                "assembly_level": lvl,
                "seq_rel_date": rel_date,
                "genome_size": gsize,
                "contig_count": ccnt,
                "gbrs_paired_asm": paired,
                "excluded_from_refseq": excl,
                "ftp_path": ftp,
            }

            is_closed_small = (lvl in ("Complete Genome", "Chromosome")) and (sym_min_bp <= gsize < sym_max_bp)
            is_named_sym = (
                (sym_min_bp <= gsize < sym_max_bp)
                and bool(_SYMBIONT_RE.search(org))
                and (ccnt == 0 or ccnt <= 50)
            )

            if is_closed_small or is_named_sym:
                counts["sym1_total_assemblies"] += 1
                prev = sym_by_sp.get(sp_taxid)
                if prev is None or _asm_rank_key(rec) < _asm_rank_key(prev):
                    sym_by_sp[sp_taxid] = rec
            else:
                counts["post_r226_total_assemblies"] += 1
                prev = post_by_sp.get(sp_taxid)
                if prev is None or _asm_rank_key(rec) < _asm_rank_key(prev):
                    post_by_sp[sp_taxid] = rec

    sym_list = sorted(sym_by_sp.values(), key=_asm_rank_key)
    post_list = sorted(post_by_sp.values(), key=_asm_rank_key)
    counts["sym1_species_reps"] = len(sym_list)
    counts["post_r226_species_reps"] = len(post_list)
    return {
        "counts": dict(counts),
        "sym1_candidates": sym_list,
        "post_r226_species_reps": post_list,
    }


def _call_stats(
    seq_strs: list[str],
    total_bp: int,
    table: int,
    *,
    count_tga: bool = False,
) -> tuple[float, float, int, float]:
    """Train pyrodigal single-mode on `seq_strs` with `table` (or meta fallback if < 20 kb)."""
    if total_bp < 20_000:
        # pyrodigal single-mode training requires >= 20 kb; pad with a repeat only for table comparison?
        # Better: use meta=True for code 11/4 or train with closed=False if >= 20 kb.
        # For < 20 kb, train() raises ValueError, so fall back to meta=True (evaluates 11 and 4).
        finder = pyrodigal.GeneFinder(meta=True)
    else:
        finder = pyrodigal.GeneFinder(meta=False)
        finder.train(*seq_strs, translation_table=table)
    coding_bp = 0
    n_genes = 0
    internal_tga = 0
    for s in seq_strs:
        for g in finder.find_genes(s):
            gl = g.end - g.begin + 1
            coding_bp += gl
            n_genes += 1
            if count_tga:
                nuc = g.sequence()
                for pos in range(0, len(nuc) - 3, 3):
                    if nuc[pos : pos + 3] == "TGA":
                        internal_tga += 1
    cd = coding_bp / total_bp if total_bp > 0 else 0.0
    mean_orf = coding_bp / n_genes if n_genes > 0 else 0.0
    tga_per_kb = (1000.0 * internal_tga / coding_bp) if (count_tga and coding_bp > 0) else 0.0
    return cd, mean_orf, n_genes, tga_per_kb


def _prefix_seqs_for_screen(seqs: list[bytes], max_bp: int) -> tuple[list[str], int, bool]:
    total_bp = sum(len(s) for s in seqs)
    if max_bp <= 0 or total_bp <= max_bp:
        return [s.upper().decode("ascii") for s in seqs], total_bp, True
    ordered = sorted(seqs, key=len, reverse=True)
    kept: list[str] = []
    acc = 0
    for s in ordered:
        kept.append(s.upper().decode("ascii"))
        acc += len(s)
        if acc >= max_bp and acc >= genecall.SINGLE_MODE_MIN_BP:
            break
    return kept, acc, False


def screen_genetic_code(
    seqs: list[bytes],
    organism_name: str = "",
    *,
    cd_gain_min: float = 0.10,
    orf_ratio_min: float = 1.25,
    cd11_trigger: float = 0.82,
    orf11_trigger: float = 600.0,
    screen_max_bp: int = 0,
) -> dict:
    """Screen a prokaryotic genome for UGA recoding (code 11 vs 4 vs 25) using pyrodigal."""
    pre_strs, pre_bp, is_full = _prefix_seqs_for_screen(seqs, screen_max_bp)
    if pre_bp < 20_000:
        cd_11, mean_orf_11, n_11, _ = _call_stats(pre_strs, pre_bp, 11)
        code = 4 if _CODE4_TAX_RE.search(organism_name) else 11
        return {
            "genetic_code": code,
            "reason": "lt20kb_tax_or_default",
            "tested_alt": False,
            "cd_11": round(cd_11, 5),
            "mean_orf_11": round(mean_orf_11, 2),
            "n_genes_11": n_11,
            "cd_25": None, "mean_orf_25": None, "n_genes_25": None,
            "cd_4": None, "mean_orf_4": None, "n_genes_4": None,
            "tga_per_kb": None,
        }

    cd_11, mean_orf_11, n_11, _ = _call_stats(pre_strs, pre_bp, 11)
    force_tax = bool(_CODE25_TAX_RE.search(organism_name) or _CODE4_TAX_RE.search(organism_name))
    triggered = force_tax or (cd_11 < cd11_trigger) or (mean_orf_11 < orf11_trigger)

    if not triggered:
        return {
            "genetic_code": 11,
            "reason": "code11_pass",
            "tested_alt": False,
            "cd_11": round(cd_11, 5),
            "mean_orf_11": round(mean_orf_11, 2),
            "n_genes_11": n_11,
            "cd_25": None, "mean_orf_25": None, "n_genes_25": None,
            "cd_4": None, "mean_orf_4": None, "n_genes_4": None,
            "tga_per_kb": None,
        }

    if not is_full:
        full_strs = [s.upper().decode("ascii") for s in seqs]
        full_bp = sum(len(s) for s in seqs)
        cd_11, mean_orf_11, n_11, _ = _call_stats(full_strs, full_bp, 11)
    else:
        full_strs, full_bp = pre_strs, pre_bp

    cd_4, mean_orf_4, n_4, _ = _call_stats(full_strs, full_bp, 4)
    cd_25, mean_orf_25, n_25, _ = _call_stats(full_strs, full_bp, 25)
    best_cd = max(cd_4, cd_25)
    best_orf = max(mean_orf_4, mean_orf_25)
    cd_gain = best_cd - cd_11
    orf_ratio = (best_orf / mean_orf_11) if mean_orf_11 > 0 else float("inf")

    tga_per_kb: float | None = None
    if cd_gain >= cd_gain_min and orf_ratio >= orf_ratio_min:
        _, _, _, tga_per_kb = _call_stats(full_strs, full_bp, 4, count_tga=True)
        if _CODE25_TAX_RE.search(organism_name):
            code = 25
        elif _CODE4_TAX_RE.search(organism_name):
            code = 4
        elif cd_25 > cd_4 + 1e-4:
            code = 25
        else:
            code = 4
        reason = f"uga_recoded_code{code}"
    else:
        code = 11
        reason = "alt_tested_kept_11"

    return {
        "genetic_code": code,
        "reason": reason,
        "tested_alt": True,
        "cd_11": round(cd_11, 5),
        "mean_orf_11": round(mean_orf_11, 2),
        "n_genes_11": n_11,
        "cd_25": round(cd_25, 5),
        "mean_orf_25": round(mean_orf_25, 2),
        "n_genes_25": n_25,
        "cd_4": round(cd_4, 5),
        "mean_orf_4": round(mean_orf_4, 2),
        "n_genes_4": n_4,
        "tga_per_kb": None if tga_per_kb is None else round(tga_per_kb, 3),
    }


def _screen_worker(args: tuple[str, list[bytes], str, int]) -> tuple[str, dict]:
    gid, seqs, org, max_bp = args
    return gid, screen_genetic_code(seqs, org, screen_max_bp=max_bp)


def _owners_from_records(
    records_tsv: Path,
    root: Path | None,
    root_stages: list[str],
    owner_inputs: dict[str, Path],
    *,
    work: Path,
    threads: int = 2,
    memory_limit: str = "4GB",
) -> dict[bytes, str]:
    if root is None and not owner_inputs:
        return {}
    duck_dir = work / f"duck_own_{records_tsv.parent.name}"
    duck_dir.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(
        config={"threads": int(threads), "memory_limit": memory_limit, "temp_directory": str(duck_dir)}
    )
    try:
        con.execute(
            f"CREATE TEMP TABLE v AS SELECT DISTINCT from_hex(contig_uid) AS contig_uid "
            f"FROM read_csv({_sql_lit(records_tsv)}, delim='\t', header=true, quote='', all_varchar=true)"
        )
        out: dict[bytes, str] = {}
        stages = list(root_stages)
        if root is not None and not stages:
            stages = list(store.stages_on_disk(root))
        for s in stages:
            if root is None:
                continue
            part_dir = Path(root) / store.INDEX_DIR / f"stage={s}"
            if not part_dir.is_dir():
                continue
            files = sorted(part_dir.glob("part-*.parquet"))
            if not files:
                continue
            for (cu,) in con.execute(
                f"SELECT DISTINCT contig_uid FROM read_parquet({_sql_list(files)}, hive_partitioning = false) "
                "WHERE is_default_variant AND contig_uid IN (SELECT contig_uid FROM v) ORDER BY 1"
            ).fetchall():
                out.setdefault(bytes(cu), s)
        for name, d in owner_inputs.items():
            gj = Path(d) / "genomes.json"
            if gj.is_file():
                p = _sql_lit(gj)
                for (cu,) in con.execute(
                    f"SELECT DISTINCT from_hex(u) FROM (SELECT unnest(contig_uids) AS u FROM read_json({p}, "
                    "format='array', columns={'contig_uids': 'VARCHAR[]'})) "
                    "WHERE from_hex(u) IN (SELECT contig_uid FROM v) ORDER BY 1"
                ).fetchall():
                    out.setdefault(bytes(cu), name)
        return out
    finally:
        con.close()


def _iter_fasta_file(path: Path):
    """Stream (header, seq_bytes, line_width) from a plain or .gz FASTA."""
    fh = gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")
    try:
        header: str | None = None
        seq_lines: list[bytes] = []
        for raw in fh:
            line = raw.rstrip(b"\r\n")
            if not line:
                continue
            if line.startswith(b">"):
                if header is not None:
                    yield header, b"".join(seq_lines), fasta._width(seq_lines)
                    seq_lines = []
                header = line[1:].decode("utf-8", "replace")
            else:
                if header is None:
                    raise Refused(f"{path}: FASTA sequence before header")
                seq_lines.append(line)
        if header is not None:
            yield header, b"".join(seq_lines), fasta._width(seq_lines)
    finally:
        fh.close()


def _resolve_within_stage_code_conflicts(
    ordered_gids: list[str],
    genome_contigs: dict[str, tuple[bytes, array.array]],
    gc_results: dict[str, dict],
) -> tuple[dict[str, int], dict[str, bytes], list[dict]]:
    genome_uids: dict[str, bytes] = {}
    resolved_code: dict[str, int] = {}
    for gid in ordered_gids:
        cu_raw, _ = genome_contigs[gid]
        uids = [cu_raw[i : i + 16] for i in range(0, len(cu_raw), 16)]
        genome_uids[gid] = ids.genome_uid(list(dict.fromkeys(uids)))
        resolved_code[gid] = int(gc_results[gid]["genetic_code"])

    contig_code_map: dict[bytes, tuple[int, str]] = {}
    conflicts: list[dict] = []
    for gid in sorted(ordered_gids, key=lambda x: genome_uids[x]):
        cu_raw, _ = genome_contigs[gid]
        c = resolved_code[gid]
        conflict_target: int | None = None
        for i in range(0, len(cu_raw), 16):
            cu = cu_raw[i : i + 16]
            prev = contig_code_map.get(cu)
            if prev is not None and prev[0] != c:
                conflict_target = prev[0]
                conflicts.append({
                    "contig_uid": cu.hex(),
                    "first_genome": prev[1],
                    "first_code": prev[0],
                    "second_genome": gid,
                    "second_code": c,
                })
        if conflict_target is not None:
            resolved_code[gid] = conflict_target
            c = conflict_target
        for i in range(0, len(cu_raw), 16):
            contig_code_map.setdefault(cu_raw[i : i + 16], (c, gid))
    return resolved_code, genome_uids, conflicts


def stage_sym1(
    selected: list[dict],
    out_root: Path,
    *,
    root: Path | None = None,
    root_stages: Iterable[str] = (),
    owner_inputs: dict[str, Path] | None = None,
    retrieved_at: str = "2026-11-24T00:00:00Z",
    workers: int = 1,
    work: Path | None = None,
    threads: int = 2,
    memory_limit: str = "4GB",
) -> dict:
    """Screen genetic codes and stage SYM1 (< 1 Mbp genome-reduced endosymbionts)."""
    out_root = Path(out_root)
    stage_out = out_root / "SYM1"
    if stage_out.exists():
        raise Refused(f"{stage_out} already exists")
    out_root.mkdir(parents=True, exist_ok=True)
    work_dir = Path(work) if work is not None else out_root / "_work_sym1"
    work_dir.mkdir(parents=True, exist_ok=True)

    tmp_d = out_root / ".SYM1.tmp"
    shutil.rmtree(tmp_d, ignore_errors=True)
    tmp_d.mkdir(parents=True)

    ordered = sorted(selected, key=lambda r: r["assembly_accession"])
    fa_name = "sym1_endosymbionts.fna.gz"
    fa_tmp = tmp_d / f".tmp.{fa_name}"

    genome_contigs: dict[str, tuple[bytes, array.array]] = {}
    genome_circular: dict[str, set[bytes]] = {}
    gc_results: dict[str, dict] = {}
    meta_by_gid: dict[str, dict] = {}

    pool = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None
    pending: list[Future] = []
    max_in_flight = max(4, workers * 4)

    def _drain(all_flag: bool = False) -> None:
        keep: list[Future] = []
        for fut in pending:
            if all_flag or fut.done():
                gid_d, res_d = fut.result()
                gc_results[gid_d] = res_d
            else:
                keep.append(fut)
        pending[:] = keep
        if not all_flag and pending:
            gid_d, res_d = pending.pop(0).result()
            gc_results[gid_d] = res_d

    fasta_uncomp_bytes = 0
    try:
        with open(fa_tmp, "wb") as raw_fh:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw_fh, mtime=0, compresslevel=6) as gz:
                for row in ordered:
                    gid = row["assembly_accession"]
                    fa_path = Path(row["fasta_path"])
                    recs = list(_iter_fasta_file(fa_path))
                    if not recs:
                        continue
                    cu_buf = bytearray()
                    len_arr = array.array("I")
                    seqs_for_screen: list[bytes] = []
                    circ_set: set[bytes] = set()
                    is_complete = row.get("assembly_level") == "Complete Genome"
                    for hdr, seq, lw in recs:
                        if not seq:
                            continue
                        cu = ids.contig_uid(seq)
                        cu_buf.extend(cu)
                        len_arr.append(len(seq))
                        seqs_for_screen.append(seq)
                        if is_complete and ("complete genome" in hdr.lower() or "complete sequence" in hdr.lower()):
                            circ_set.add(cu)
                        out_bytes = b"\n".join(fasta.format_record(hdr, seq, lw)) + b"\n"
                        gz.write(out_bytes)
                        fasta_uncomp_bytes += len(out_bytes)
                    if not seqs_for_screen:
                        continue
                    genome_contigs[gid] = (bytes(cu_buf), len_arr)
                    genome_circular[gid] = circ_set
                    meta_by_gid[gid] = row
                    if pool is not None:
                        pending.append(pool.submit(_screen_worker, (gid, seqs_for_screen, row.get("organism_name", ""), 0)))
                        if len(pending) >= max_in_flight:
                            _drain(False)
                    else:
                        gc_results[gid] = screen_genetic_code(seqs_for_screen, row.get("organism_name", ""))
        if pool is not None:
            _drain(True)
            pool.shutdown(wait=True)
            pool = None
    finally:
        if pool is not None:
            pool.shutdown(wait=False, cancel_futures=True)

    fa_tmp.replace(tmp_d / fa_name)
    kept_gids = [r["assembly_accession"] for r in ordered if r["assembly_accession"] in genome_contigs]
    resolved_code, genome_uids, conflicts = _resolve_within_stage_code_conflicts(kept_gids, genome_contigs, gc_results)

    records_path = tmp_d / "records.tsv"
    gc_path = tmp_d / "genetic_codes.tsv"
    meta_path = tmp_d / "sym1_metadata.tsv"
    with open(records_path, "w", encoding="utf-8") as rf, open(gc_path, "w", encoding="utf-8") as gf, open(
        meta_path, "w", encoding="utf-8"
    ) as mf:
        rf.write("ordinal\tgenome\trecord_ordinal\tcontig_uid\tlength\tassembly_quality\tgenetic_code\n")
        gf.write("\t".join(GC_HEAD) + "\n")
        mf.write(
            "assembly_accession\tspecies_taxid\ttaxid\torganism_name\tassembly_level\tseq_rel_date\t"
            "genome_size\tcontig_count\tgenetic_code\tcode_reason\tcircular_contigs\n"
        )
        ordinal = 0
        for gid in kept_gids:
            r = meta_by_gid[gid]
            cu_raw, len_arr = genome_contigs[gid]
            code = resolved_code[gid]
            gc = gc_results[gid]
            aq = _ASSEMBLY_QUALITY_MAP.get(r.get("assembly_level", ""), "contig")
            tot_bp = sum(len_arr)
            gf.write("\t".join([
                gid, "SYM1", r.get("organism_name", ""), str(tot_bp), str(len(len_arr)), aq,
                str(code), gc["reason"], str(gc["tested_alt"]).lower(),
                _fmt_opt(gc["cd_11"]), _fmt_opt(gc["mean_orf_11"]), _fmt_opt(gc["n_genes_11"]),
                _fmt_opt(gc["cd_25"]), _fmt_opt(gc["mean_orf_25"]), _fmt_opt(gc["n_genes_25"]),
                _fmt_opt(gc["cd_4"]), _fmt_opt(gc["mean_orf_4"]), _fmt_opt(gc["n_genes_4"]),
                _fmt_opt(gc["tga_per_kb"]),
            ]) + "\n")
            mf.write("\t".join([
                gid, str(r.get("species_taxid", "")), str(r.get("taxid", "")), r.get("organism_name", ""),
                r.get("assembly_level", ""), r.get("seq_rel_date", ""),
                str(tot_bp), str(len(len_arr)), str(code), gc["reason"], str(len(genome_circular[gid])),
            ]) + "\n")
            for rec_ord, length in enumerate(len_arr):
                cu = cu_raw[rec_ord * 16 : (rec_ord + 1) * 16]
                rf.write(f"{ordinal}\t{gid}\t{rec_ord}\t{cu.hex()}\t{length}\t{aq}\t{code}\n")
                ordinal += 1

    owners = _owners_from_records(
        records_path, root, list(root_stages), dict(owner_inputs or {}),
        work=work_dir, threads=threads, memory_limit=memory_limit,
    )

    listing = make_inputs.Listing()
    genomes_list: list[dict] = []
    repeated: list[tuple] = []
    shared: list[tuple] = []
    first_of: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}
    ordinal = 0
    for gid in kept_gids:
        r = meta_by_gid[gid]
        cu_raw, len_arr = genome_contigs[gid]
        uids = [cu_raw[i : i + 16] for i in range(0, len(cu_raw), 16)]
        code = resolved_code[gid]
        aq = _ASSEMBLY_QUALITY_MAP.get(r.get("assembly_level", ""), "contig")
        for rec_ord, (cu, length) in enumerate(zip(uids, len_arr)):
            contig_len.setdefault(cu, int(length))
            if cu in owners:
                shared.append((ordinal, gid, rec_ord, cu.hex(), int(length), owners[cu]))
            ordinal += 1
        circ_hex = sorted({cu.hex() for cu in genome_circular[gid]})
        g = listing.add(uids, code, "prokaryote", aq)
        guid = genome_uids[gid]
        if g is None:
            repeated.append((ordinal - len(uids), gid, guid.hex(), first_of[guid]))
            continue
        g["circular"] = circ_hex
        genomes_list.append(g)
        first_of[guid] = gid

    make_inputs._write_tsv(tmp_d / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
    make_inputs._write_tsv(
        tmp_d / "shared_records.tsv",
        ["ordinal", "genome", "record_ordinal", "contig_uid", "length", "owner_stage"],
        shared,
    )

    objects = [{
        "path": str(stage_out / fa_name),
        **SYM1_SOURCE_BASE,
        "retrieved_at": retrieved_at,
    }]
    (tmp_d / "objects.json").write_text(json.dumps(objects, indent=1))
    (tmp_d / "genomes.json").write_text(json.dumps(genomes_list))
    tmp_d.rename(stage_out)

    owner_contigs: Counter = Counter()
    owner_bases: Counter = Counter()
    for cu, ow in owners.items():
        if cu in contig_len:
            owner_contigs[ow] += 1
            owner_bases[ow] += contig_len[cu]

    gc_counts = Counter(str(resolved_code[gid]) for gid in kept_gids)
    return {
        "stage": "SYM1",
        "genomes_selected": len(kept_gids),
        "genomes_distinct": len(genomes_list),
        "genomes_repeated": len(repeated),
        "contigs_distinct": len(contig_len),
        "bases_distinct": sum(contig_len.values()),
        "contigs_held_earlier": sum(owner_contigs.values()),
        "bases_held_earlier": sum(owner_bases.values()),
        "contigs_new": sum(1 for cu in contig_len if cu not in owners),
        "bases_new": sum(n for cu, n in contig_len.items() if cu not in owners),
        "held_by": dict(sorted(owner_contigs.items())),
        "held_bases_by": dict(sorted(owner_bases.items())),
        "genetic_codes": dict(sorted(gc_counts.items())),
        "code_conflicts": len(conflicts),
        "fasta_bytes": (stage_out / fa_name).stat().st_size,
        "fasta_sha256": hashlib.sha256((stage_out / fa_name).read_bytes()).hexdigest(),
        "genomes_json_sha256": hashlib.sha256((stage_out / "genomes.json").read_bytes()).hexdigest(),
        "objects_json_sha256": hashlib.sha256((stage_out / "objects.json").read_bytes()).hexdigest(),
        "fasta_uncomp_bytes": fasta_uncomp_bytes,
    }


def stage_sag1(
    gorg_fasta_gz: Path,
    out_root: Path,
    *,
    fasta_dir: Path | None = None,
    final_dir: Path | None = None,
    sag_taxonomy: dict[str, str] | None = None,
    root: Path | None = None,
    root_stages: Iterable[str] = (),
    owner_inputs: dict[str, Path] | None = None,
    retrieved_at: str = "2026-11-24T00:00:00Z",
    workers: int = 1,
    screen_max_bp: int = 500_000,
    work: Path | None = None,
    threads: int = 2,
    memory_limit: str = "4GB",
) -> dict:
    """Stream GORG_v1.fasta.gz, group contigs by SAG ID, screen genetic codes, and stage SAG1."""
    out_root = Path(out_root)
    stage_out = out_root / "SAG1"
    if stage_out.exists():
        raise Refused(f"{stage_out} already exists")
    out_root.mkdir(parents=True, exist_ok=True)
    work_dir = Path(work) if work is not None else out_root / "_work_sag1"
    work_dir.mkdir(parents=True, exist_ok=True)

    fa_dir = Path(fasta_dir) if fasta_dir is not None else None
    if fa_dir is not None:
        fa_dir.mkdir(parents=True, exist_ok=True)

    tmp_d = out_root / ".SAG1.tmp"
    shutil.rmtree(tmp_d, ignore_errors=True)
    tmp_d.mkdir(parents=True)

    fa_name = "sag1_gorg.fna.gz"
    dst_fa_dir = fa_dir if fa_dir is not None else tmp_d
    fa_tmp = dst_fa_dir / f".tmp.{fa_name}"
    fa_final = (fa_dir / fa_name) if fa_dir is not None else (stage_out / fa_name)

    # First pass: group contigs per SAG ID? In GORG_v1.fasta.gz contigs of the same SAG are usually contiguous,
    # but we handle arbitrary order cleanly by indexing SAG -> contigs in bounded batches if contiguous,
    # or writing consolidated deterministic .fna.gz in sorted SAG order.
    # Actually GORG_v1.fasta.gz is ~1.44 GB compressed (~5.5 GB uncompressed); to stay strictly bounded-memory
    # (AGENTS §4.1), we stream once: write contigs to sag1_gorg.fna.gz as they arrive (deterministic input order)
    # and flush each SAG for genetic-code screening as soon as its contiguous run ends (falling back to a second
    # pass only for any SAG whose contigs were non-contiguous).
    sag_contigs: dict[str, tuple[bytearray, array.array]] = {}
    gc_results: dict[str, dict] = {}
    sag_order: list[str] = []
    non_contiguous_sags: set[str] = set()
    tax_map = sag_taxonomy or {}

    pool = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None
    pending: list[Future] = []
    max_in_flight = max(4, workers * 4)

    def _drain(all_flag: bool = False) -> None:
        keep: list[Future] = []
        for fut in pending:
            if all_flag or fut.done():
                gid_d, res_d = fut.result()
                gc_results[gid_d] = res_d
            else:
                keep.append(fut)
        pending[:] = keep
        if not all_flag and pending:
            gid_d, res_d = pending.pop(0).result()
            gc_results[gid_d] = res_d

    def _submit_screen(gid: str, seqs: list[bytes]) -> None:
        org = tax_map.get(gid, "")
        if pool is not None:
            pending.append(pool.submit(_screen_worker, (gid, seqs, org, screen_max_bp)))
            if len(pending) >= max_in_flight:
                _drain(False)
        else:
            gc_results[gid] = screen_genetic_code(seqs, org, screen_max_bp=screen_max_bp)

    fasta_uncomp_bytes = 0
    cur_sag: str | None = None
    cur_seqs: list[bytes] = []
    finished_sags: set[str] = set()

    try:
        with open(fa_tmp, "wb", buffering=1 << 22) as raw_fh:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw_fh, mtime=0, compresslevel=6) as gz:
                for hdr, seq, lw in _iter_fasta_file(Path(gorg_fasta_gz)):
                    if not seq:
                        continue
                    sag_id = parse_gorg_sag_id(hdr)
                    if sag_id != cur_sag:
                        if cur_sag is not None:
                            if cur_sag not in non_contiguous_sags:
                                _submit_screen(cur_sag, cur_seqs)
                            finished_sags.add(cur_sag)
                        cur_sag = sag_id
                        cur_seqs = []
                        if sag_id in finished_sags:
                            non_contiguous_sags.add(sag_id)
                        elif sag_id not in sag_contigs:
                            sag_order.append(sag_id)
                            sag_contigs[sag_id] = (bytearray(), array.array("I"))
                    cu = ids.contig_uid(seq)
                    buf, larr = sag_contigs[sag_id]
                    buf.extend(cu)
                    larr.append(len(seq))
                    cur_seqs.append(seq)
                    out_bytes = b"\n".join(fasta.format_record(hdr, seq, lw)) + b"\n"
                    gz.write(out_bytes)
                    fasta_uncomp_bytes += len(out_bytes)
                if cur_sag is not None and cur_sag not in non_contiguous_sags:
                    _submit_screen(cur_sag, cur_seqs)

        if non_contiguous_sags:
            # Re-screen any SAGs whose contigs were split across the input FASTA
            split_seqs: dict[str, list[bytes]] = defaultdict(list)
            for hdr, seq, _ in _iter_fasta_file(Path(gorg_fasta_gz)):
                if not seq:
                    continue
                sag_id = parse_gorg_sag_id(hdr)
                if sag_id in non_contiguous_sags:
                    split_seqs[sag_id].append(seq)
            for sag_id in sorted(non_contiguous_sags):
                _submit_screen(sag_id, split_seqs[sag_id])

        if pool is not None:
            _drain(True)
            pool.shutdown(wait=True)
            pool = None
    finally:
        if pool is not None:
            pool.shutdown(wait=False, cancel_futures=True)

    if fa_dir is not None:
        fa_tmp.replace(fa_final)
    else:
        fa_tmp.replace(tmp_d / fa_name)

    ordered_sags = sorted(sag_order)
    frozen_contigs = {gid: (bytes(b), a) for gid, (b, a) in sag_contigs.items()}
    resolved_code, genome_uids, conflicts = _resolve_within_stage_code_conflicts(ordered_sags, frozen_contigs, gc_results)

    records_path = tmp_d / "records.tsv"
    gc_path = tmp_d / "genetic_codes.tsv"
    with open(records_path, "w", encoding="utf-8", buffering=1 << 20) as rf, open(
        gc_path, "w", encoding="utf-8", buffering=1 << 20
    ) as gf:
        rf.write("ordinal\tgenome\trecord_ordinal\tcontig_uid\tlength\tassembly_quality\tgenetic_code\n")
        gf.write("\t".join(GC_HEAD) + "\n")
        ordinal = 0
        for gid in ordered_sags:
            cu_raw, len_arr = frozen_contigs[gid]
            code = resolved_code[gid]
            gc = gc_results[gid]
            tot_bp = sum(len_arr)
            gf.write("\t".join([
                gid, "SAG1", tax_map.get(gid, ""), str(tot_bp), str(len(len_arr)), "sag",
                str(code), gc["reason"], str(gc["tested_alt"]).lower(),
                _fmt_opt(gc["cd_11"]), _fmt_opt(gc["mean_orf_11"]), _fmt_opt(gc["n_genes_11"]),
                _fmt_opt(gc["cd_25"]), _fmt_opt(gc["mean_orf_25"]), _fmt_opt(gc["n_genes_25"]),
                _fmt_opt(gc["cd_4"]), _fmt_opt(gc["mean_orf_4"]), _fmt_opt(gc["n_genes_4"]),
                _fmt_opt(gc["tga_per_kb"]),
            ]) + "\n")
            for rec_ord, length in enumerate(len_arr):
                cu = cu_raw[rec_ord * 16 : (rec_ord + 1) * 16]
                rf.write(f"{ordinal}\t{gid}\t{rec_ord}\t{cu.hex()}\t{length}\tsag\t{code}\n")
                ordinal += 1

    owners = _owners_from_records(
        records_path, root, list(root_stages), dict(owner_inputs or {}),
        work=work_dir, threads=threads, memory_limit=memory_limit,
    )

    listing = make_inputs.Listing()
    genomes_list: list[dict] = []
    repeated: list[tuple] = []
    shared: list[tuple] = []
    first_of: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}
    ordinal = 0
    for gid in ordered_sags:
        cu_raw, len_arr = frozen_contigs[gid]
        uids = [cu_raw[i : i + 16] for i in range(0, len(cu_raw), 16)]
        code = resolved_code[gid]
        for rec_ord, (cu, length) in enumerate(zip(uids, len_arr)):
            contig_len.setdefault(cu, int(length))
            if cu in owners:
                shared.append((ordinal, gid, rec_ord, cu.hex(), int(length), owners[cu]))
            ordinal += 1
        # genome_input.assembly_quality allows 'contig' for SAG drafts (check schemas.py!)
        g = listing.add(uids, code, "prokaryote", "contig")
        guid = genome_uids[gid]
        if g is None:
            repeated.append((ordinal - len(uids), gid, guid.hex(), first_of[guid]))
            continue
        genomes_list.append(g)
        first_of[guid] = gid

    make_inputs._write_tsv(tmp_d / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
    make_inputs._write_tsv(
        tmp_d / "shared_records.tsv",
        ["ordinal", "genome", "record_ordinal", "contig_uid", "length", "owner_stage"],
        shared,
    )

    final_obj_dir = Path(final_dir) if final_dir is not None else (fa_dir if fa_dir is not None else stage_out)
    objects = [{
        "path": str(final_obj_dir / fa_name),
        **SAG1_SOURCE_BASE,
        "retrieved_at": retrieved_at,
    }]
    (tmp_d / "objects.json").write_text(json.dumps(objects, indent=1))
    if final_dir is not None and fa_dir is not None:
        staging_objects = [{
            "path": str(fa_final),
            **SAG1_SOURCE_BASE,
            "retrieved_at": retrieved_at,
        }]
        (tmp_d / "objects_staging.json").write_text(json.dumps(staging_objects, indent=1))
    (tmp_d / "genomes.json").write_text(json.dumps(genomes_list))
    tmp_d.rename(stage_out)

    owner_contigs: Counter = Counter()
    owner_bases: Counter = Counter()
    for cu, ow in owners.items():
        if cu in contig_len:
            owner_contigs[ow] += 1
            owner_bases[ow] += contig_len[cu]

    gc_counts = Counter(str(resolved_code[gid]) for gid in ordered_sags)
    return {
        "stage": "SAG1",
        "genomes_total": len(ordered_sags),
        "genomes_distinct": len(genomes_list),
        "genomes_repeated": len(repeated),
        "contigs_distinct": len(contig_len),
        "bases_distinct": sum(contig_len.values()),
        "contigs_held_earlier": sum(owner_contigs.values()),
        "bases_held_earlier": sum(owner_bases.values()),
        "contigs_new": sum(1 for cu in contig_len if cu not in owners),
        "bases_new": sum(n for cu, n in contig_len.items() if cu not in owners),
        "held_by": dict(sorted(owner_contigs.items())),
        "held_bases_by": dict(sorted(owner_bases.items())),
        "genetic_codes": dict(sorted(gc_counts.items())),
        "code_conflicts": len(conflicts),
        "fasta_path": str(fa_final),
        "fasta_bytes": fa_final.stat().st_size,
        "fasta_sha256": _sha256_file(fa_final),
        "genomes_json_sha256": _sha256_file(stage_out / "genomes.json"),
        "objects_json_sha256": _sha256_file(stage_out / "objects.json"),
        "fasta_uncomp_bytes": fasta_uncomp_bytes,
    }


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def niche_compendia_register() -> list[dict]:
    """Ordered register of niche-habitat & holobiont MAG compendia (NX1..) outside GlobDB r226,
    aligned with card C-66 (which owns NX1).
    """
    return [
        {
            "stage": "NX1",
            "owner_card": "C-66",
            "compendia": (
                "MEER Hadal Trench Microbiome (Cell 2025, OEP004067 / deepseadb, 7,564 piezophile MAGs > 6,000 m), "
                "eHOMD v11.03 (955 oral MAGs + 25,388 oral phages), cyanobacterial mats/biocrusts (figshare.21932371), "
                "serpentinizing subsurface (Lost City zenodo.5798015, Samail zenodo.8108239), marine sponges "
                "(figshare.29695256), lichens (figshare.27054937), and rhizosphere (zenodo.16947911)"
            ),
            "est_gbp": "15-30 Gbp",
            "globdb_overlap": "low (post-GlobDB r226 2024-2025 releases and excluded oral/hadal/holobiont biomes)",
            "dedup_rule": "skani 96% ANI / 50% AF vs GTDB R226 + GlobDB r226; contig_uid vs built/staged",
        },
        {
            "stage": "NX2",
            "owner_card": "C-58 follow-up",
            "compendia": (
                "Deep-sea hydrothermal vent & cold-seep sediment MAG compendia (Global Cold Seep MAG catalogs "
                "figshare.22568107 / figshare.32788992; Guaymas Basin, Mid-Cayman Rise, and Okinawa Trough vent MAGs "
                "via ENA/NCBI BioProjects)"
            ),
            "est_gbp": "10-25 Gbp",
            "globdb_overlap": "low (GlobDB r226 has no dedicated hydrothermal-vent or cold-seep catalogue)",
            "dedup_rule": "skani 96% ANI / 50% AF vs GTDB R226 + GlobDB r226 + NX1; contig_uid dedup",
        },
        {
            "stage": "NX3",
            "owner_card": "C-58 follow-up",
            "compendia": (
                "Hypersaline/soda-lake, acid mine drainage (AMD: PRJNA666025/PRJNA666095 and Iberian Pyrite Belt), "
                "and polar/permafrost/cryosphere MAG compendia beyond GlobDB's GFS glacier set (Southern Ocean, "
                "Antarctic dry valleys, Arctic permafrost outside EMERGE LR2)"
            ),
            "est_gbp": "10-25 Gbp",
            "globdb_overlap": "low (outside GFS glacier-fed streams in GlobDB r226)",
            "dedup_rule": "skani 96% ANI / 50% AF vs GTDB R226 + GlobDB r226 + NX1..NX2; contig_uid dedup",
        },
        {
            "stage": "NX4",
            "owner_card": "C-58 follow-up",
            "compendia": (
                "Insect, termite-gut, coral, and non-human animal holobiont MAGs (termite hindgut Spirochaetota/Fibrobacterota "
                "and protist ecto/endosymbionts, reef coral microbiome MAGs) plus fermented-food MAGs beyond cFMD and "
                "ancient paleomicrobiome MAGs (Pleistocene/Neanderthal dental calculus & coprolite MAGs)"
            ),
            "est_gbp": "10-20 Gbp",
            "globdb_overlap": "low (outside GlobDB r226's human/ruminant/cFMD catalogues)",
            "dedup_rule": "skani 96% ANI / 50% AF vs GTDB R226 + GlobDB r226 + NX1..NX3; contig_uid dedup",
        },
    ]


# ---------------------------------------------------------------------------
# Network fetchers (<= 2 req/s, User-Agent: curl, stop after 5 refusals)
# ---------------------------------------------------------------------------

def _download_stream(
    url: str,
    dest: Path,
    *,
    gzip_output: bool = False,
    min_interval_s: float = 0.6,
    last_req: list[float] | None = None,
) -> dict:
    if last_req is not None:
        wait = min_interval_s - (time.monotonic() - last_req[0])
        if wait > 0:
            time.sleep(wait)
        last_req[0] = time.monotonic()

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.parent / f".tmp.{dest.name}"
    req = urllib.request.Request(url, headers={"User-Agent": "curl"})
    h_raw = hashlib.sha256()
    md5_raw = hashlib.md5()
    raw_bytes = 0
    with urllib.request.urlopen(req, timeout=180) as resp:
        status = resp.status
        last_mod = resp.headers.get("Last-Modified", "")
        if gzip_output:
            with open(tmp, "wb") as out_f, gzip.GzipFile(filename="", mode="wb", fileobj=out_f, mtime=0, compresslevel=4) as gz:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    h_raw.update(chunk)
                    md5_raw.update(chunk)
                    raw_bytes += len(chunk)
                    gz.write(chunk)
        else:
            with open(tmp, "wb") as out_f:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    h_raw.update(chunk)
                    md5_raw.update(chunk)
                    raw_bytes += len(chunk)
                    out_f.write(chunk)
    tmp.replace(dest)
    return {
        "url": url,
        "dest": str(dest),
        "status": status,
        "last_modified": last_mod,
        "raw_bytes": raw_bytes,
        "disk_bytes": dest.stat().st_size,
        "raw_sha256": h_raw.hexdigest(),
        "raw_md5": md5_raw.hexdigest(),
    }


def cmd_fetch_tables(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    last_req = [0.0]
    refusals = 0
    log: list[dict] = []

    for db, domain, url in NCBI_SUMMARY_URLS:
        dest = out_dir / f"ncbi_{db}_{domain}_assembly_summary.txt.gz"
        if dest.is_file() and dest.stat().st_size > 1000:
            print(f"[skip] {dest} already exists ({dest.stat().st_size} B)", file=sys.stderr)
            continue
        try:
            rec = _download_stream(url, dest, gzip_output=True, last_req=last_req)
            refusals = 0
            log.append(rec)
            print(f"[ok] {db}/{domain}: {rec['raw_bytes']} B raw -> {rec['disk_bytes']} B gz sha256={rec['raw_sha256']}", file=sys.stderr)
        except Exception as e:
            refusals += 1
            print(f"[err] {url}: {e}", file=sys.stderr)
            if refusals >= 5:
                raise Refused("5 consecutive refusals fetching tables")

    for bp, url in ENA_GORG_URLS:
        dest = out_dir / f"ena_{bp}_assemblies.tsv"
        if dest.is_file() and dest.stat().st_size > 100:
            print(f"[skip] {dest} already exists ({dest.stat().st_size} B)", file=sys.stderr)
            continue
        try:
            rec = _download_stream(url, dest, gzip_output=False, last_req=last_req)
            refusals = 0
            log.append(rec)
            print(f"[ok] ENA {bp}: {rec['raw_bytes']} B sha256={rec['raw_sha256']}", file=sys.stderr)
        except Exception as e:
            refusals += 1
            print(f"[err] {url}: {e}", file=sys.stderr)
            if refusals >= 5:
                raise Refused("5 consecutive refusals fetching ENA tables")

    log_path = out_dir / "fetch_tables_log.json"
    log_path.write_text(json.dumps(log, indent=2))


def cmd_diff_ncbi(args: argparse.Namespace) -> None:
    tdir = Path(args.tables_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[c58] loading GTDB R226 accessions...", file=sys.stderr, flush=True)
    gtdb_accs = load_gtdb_base_accessions(Path(args.gtdb_tar))
    print(f"[c58] GTDB base accessions: {len(gtdb_accs):,}", file=sys.stderr, flush=True)

    print("[c58] loading GlobDB r226 accessions...", file=sys.stderr, flush=True)
    globdb_accs = load_globdb_base_accessions(Path(args.globdb_tax), Path(args.globdb_dict_tar))
    print(f"[c58] GlobDB base accessions: {len(globdb_accs):,}", file=sys.stderr, flush=True)

    summaries = [
        (db, dom, tdir / f"ncbi_{db}_{dom}_assembly_summary.txt.gz")
        for db, dom, _ in NCBI_SUMMARY_URLS
    ]
    print("[c58] diffing NCBI prokaryotic assembly summaries...", file=sys.stderr, flush=True)
    res = diff_ncbi_summaries(summaries, gtdb_base_accs=gtdb_accs, globdb_base_accs=globdb_accs)

    cols = [
        "assembly_accession", "db", "domain", "taxid", "species_taxid", "organism_name",
        "assembly_level", "seq_rel_date", "genome_size", "contig_count",
        "gbrs_paired_asm", "excluded_from_refseq", "ftp_path",
    ]
    sym_tsv = out_dir / "sym1_candidates.tsv"
    with open(sym_tsv, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in res["sym1_candidates"]:
            w.writerow(r)

    post_gz = out_dir / "SYM1_post_r226_left_for_later.tsv.gz"
    with open(post_gz, "wb") as raw_fh:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_fh, mtime=0, compresslevel=6) as gz:
            tw = io.TextIOWrapper(gz, encoding="utf-8", newline="")
            w = csv.DictWriter(tw, fieldnames=["rank"] + cols, delimiter="\t")
            w.writeheader()
            for idx, r in enumerate(res["post_r226_species_reps"], 1):
                w.writerow({"rank": idx, **r})
            tw.flush()

    summary = {
        "gtdb_base_accessions": len(gtdb_accs),
        "globdb_base_accessions": len(globdb_accs),
        "counts": res["counts"],
        "sym1_candidates_tsv": str(sym_tsv),
        "sym1_candidates_sha256": _sha256_file(sym_tsv),
        "sym1_total_bp": sum(r["genome_size"] for r in res["sym1_candidates"]),
        "post_r226_tsv_gz": str(post_gz),
        "post_r226_bytes": post_gz.stat().st_size,
        "post_r226_sha256": _sha256_file(post_gz),
        "post_r226_total_bp": sum(r["genome_size"] for r in res["post_r226_species_reps"]),
    }
    (out_dir / "ncbi_diff_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


def cmd_fetch_sym1(args: argparse.Namespace) -> None:
    cand_tsv = Path(args.candidates_tsv)
    out_dir = Path(args.out_dir)
    fa_dir = out_dir / "fna"
    fa_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    with open(cand_tsv, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(r)
    if args.max_genomes and args.max_genomes > 0:
        rows = rows[: args.max_genomes]

    by_acc = {r["assembly_accession"]: r for r in rows}
    fetched: dict[str, dict] = {}
    not_served: list[str] = []
    last_req = [0.0]
    refusals = 0

    # Check already downloaded
    todo: list[str] = []
    for r in rows:
        acc = r["assembly_accession"]
        dst = fa_dir / f"{acc}.fna.gz"
        if dst.is_file() and dst.stat().st_size > 100:
            fetched[acc] = {"fasta_path": str(dst), "bytes": dst.stat().st_size}
        else:
            todo.append(acc)

    # Fetch in batches of 20 via NCBI Datasets v2
    batch_size = 20
    for i in range(0, len(todo), batch_size):
        batch = todo[i : i + batch_size]
        wait = 0.6 - (time.monotonic() - last_req[0])
        if wait > 0:
            time.sleep(wait)
        last_req[0] = time.monotonic()

        accs_str = ",".join(batch)
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accs_str}/download?include_annotation_type=GENOME_FASTA"
        req = urllib.request.Request(url, headers={"User-Agent": "curl"})
        try:
            with urllib.request.urlopen(req, timeout=180) as resp:
                zbytes = resp.read()
            refusals = 0
            got_in_batch: set[str] = set()
            with zipfile.ZipFile(io.BytesIO(zbytes)) as zf:
                for name in zf.namelist():
                    if not name.endswith(".fna"):
                        continue
                    # ncbi_dataset/data/<ACC>/<ACC>_*_genomic.fna
                    parts = name.split("/")
                    acc = next((p for p in parts if p in by_acc), None)
                    if acc is None:
                        continue
                    raw = zf.read(name)
                    if not raw.strip():
                        continue
                    dst = fa_dir / f"{acc}.fna.gz"
                    tmp = fa_dir / f".tmp.{acc}.fna.gz"
                    with open(tmp, "wb") as of, gzip.GzipFile(filename="", mode="wb", fileobj=of, mtime=0, compresslevel=6) as gz:
                        gz.write(raw)
                    tmp.replace(dst)
                    fetched[acc] = {"fasta_path": str(dst), "bytes": dst.stat().st_size}
                    got_in_batch.add(acc)
            for acc in batch:
                if acc not in got_in_batch:
                    not_served.append(acc)
        except Exception as e:
            refusals += 1
            print(f"[warn] Datasets batch {i//batch_size} failed: {e}", file=sys.stderr)
            not_served.extend(batch)
            if refusals >= 5:
                print("[warn] 5 consecutive Datasets failures; falling back to FTP for remaining", file=sys.stderr)
                not_served.extend(todo[i + batch_size :])
                break

    # FTP fallback for any accession not served by Datasets
    still_missing: list[str] = []
    refusals = 0
    for acc in not_served:
        if acc in fetched:
            continue
        ftp = by_acc[acc].get("ftp_path", "").strip()
        if not ftp or ftp == "na":
            still_missing.append(acc)
            continue
        if ftp.startswith("ftp://"):
            ftp = "https://" + ftp[6:]
        asm_dir = ftp.rstrip("/").rsplit("/", 1)[-1]
        url = f"{ftp.rstrip('/')}/{asm_dir}_genomic.fna.gz"
        dst = fa_dir / f"{acc}.fna.gz"
        try:
            _download_stream(url, dst, gzip_output=False, last_req=last_req)
            refusals = 0
            fetched[acc] = {"fasta_path": str(dst), "bytes": dst.stat().st_size}
        except Exception as e:
            refusals += 1
            still_missing.append(acc)
            print(f"[warn] FTP fallback {acc}: {e}", file=sys.stderr)
            if refusals >= 5:
                print("[warn] 5 consecutive FTP failures; stopping fallback", file=sys.stderr)
                break

    manifest_tsv = out_dir / "sym1_fetched_manifest.tsv"
    cols = list(rows[0].keys()) + ["fasta_path"] if rows else ["assembly_accession", "fasta_path"]
    with open(manifest_tsv, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            acc = r["assembly_accession"]
            if acc in fetched:
                w.writerow({**r, "fasta_path": fetched[acc]["fasta_path"]})

    (out_dir / "sym1_not_served.txt").write_text("\n".join(sorted(set(still_missing))) + ("\n" if still_missing else ""))
    print(json.dumps({
        "requested": len(rows),
        "fetched": len(fetched),
        "not_served": len(set(still_missing)),
        "manifest_tsv": str(manifest_tsv),
    }, indent=2))


def cmd_stage_sym1(args: argparse.Namespace) -> None:
    rows: list[dict] = []
    with open(args.manifest_tsv, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["genome_size"] = _safe_int(r.get("genome_size", "0"))
            rows.append(r)

    owner_inputs = {}
    for item in args.owner_inputs or []:
        k, v = item.split("=", 1)
        owner_inputs[k] = Path(v)

    rep = stage_sym1(
        rows,
        Path(args.out_root),
        root=Path(args.root) if args.root else None,
        root_stages=args.root_stages.split(",") if args.root_stages else (),
        owner_inputs=owner_inputs,
        retrieved_at=args.retrieved_at,
        workers=args.workers,
        work=Path(args.work) if args.work else None,
        threads=args.threads,
        memory_limit=args.memory_limit,
    )
    if args.report_out:
        Path(args.report_out).write_text(json.dumps(rep, indent=2))
    print(json.dumps(rep, indent=2))


def cmd_fetch_gorg(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    last_req = [0.0]
    log: list[dict] = []
    for fname, url, exp_bytes, exp_md5, exp_sha256 in GORG_OSF_FILES:
        dest = out_dir / fname
        if dest.is_file() and dest.stat().st_size == exp_bytes:
            sha = _sha256_file(dest)
            if sha == exp_sha256:
                print(f"[skip] {dest} verified sha256={sha}", file=sys.stderr)
                log.append({"file": fname, "bytes": exp_bytes, "sha256": sha, "md5": exp_md5, "verified": True})
                continue
        rec = _download_stream(url, dest, gzip_output=False, last_req=last_req)
        if rec["raw_bytes"] != exp_bytes or rec["raw_md5"] != exp_md5 or rec["raw_sha256"] != exp_sha256:
            raise Refused(
                f"{fname}: checksum mismatch got bytes={rec['raw_bytes']} md5={rec['raw_md5']} sha256={rec['raw_sha256']} "
                f"expected bytes={exp_bytes} md5={exp_md5} sha256={exp_sha256}"
            )
        rec["verified"] = True
        log.append(rec)
        print(f"[ok] {fname}: {rec['raw_bytes']} B verified sha256={rec['raw_sha256']}", file=sys.stderr)
    (out_dir / "gorg_fetch_log.json").write_text(json.dumps(log, indent=2))


def compute_gorg_metadata_overlap(
    ena_tsv: Path,
    *,
    gtdb_accs: set[str],
    globdb_accs: set[str],
    iuig_accs: set[str],
) -> dict:
    """Count accession-level overlap of ENA GORG-Tropics / GORG-Dark assemblies vs GTDB, GlobDB, IUIG."""
    total = 0
    in_gtdb = 0
    in_globdb = 0
    in_iuig = 0
    in_any = 0
    total_bp = 0
    with open(ena_tsv, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            total += 1
            base = strip_acc_version(r.get("accession", ""))
            total_bp += _safe_int(r.get("base_count", "0"))
            g1 = base in gtdb_accs
            g2 = base in globdb_accs
            g3 = base in iuig_accs
            if g1:
                in_gtdb += 1
            if g2:
                in_globdb += 1
            if g3:
                in_iuig += 1
            if g1 or g2 or g3:
                in_any += 1
    return {
        "total_assemblies": total,
        "total_bp": total_bp,
        "in_gtdb_r226": in_gtdb,
        "in_globdb_r226": in_globdb,
        "in_iuig": in_iuig,
        "in_any": in_any,
        "novel_vs_all_three": total - in_any,
    }


def load_iuig_base_accessions(iuig_meta_zip: Path) -> set[str]:
    accs: set[str] = set()
    if not iuig_meta_zip.is_file():
        return accs
    with zipfile.ZipFile(iuig_meta_zip) as zf:
        for name in zf.namelist():
            if "img_genomes_tarballs_metadata" in name:
                with zf.open(name) as fh:
                    for raw in fh:
                        line = raw.decode("utf-8", "replace")
                        for m in re.finditer(r"\bGC[AF]_\d+", line):
                            accs.add(m.group(0))
    return accs


def cmd_stage_sag1(args: argparse.Namespace) -> None:
    owner_inputs = {}
    for item in args.owner_inputs or []:
        k, v = item.split("=", 1)
        owner_inputs[k] = Path(v)

    sag_tax: dict[str, str] = {}
    if args.gorg_ena_tsv and Path(args.gorg_ena_tsv).is_file():
        with open(args.gorg_ena_tsv, encoding="utf-8", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                aname = r.get("assembly_name", "").strip()
                if aname:
                    sag_tax[aname] = r.get("scientific_name", "").strip()

    rep = stage_sag1(
        Path(args.gorg_fasta_gz),
        Path(args.out_root),
        fasta_dir=Path(args.fasta_dir) if args.fasta_dir else None,
        final_dir=Path(args.final_dir) if args.final_dir else None,
        sag_taxonomy=sag_tax,
        root=Path(args.root) if args.root else None,
        root_stages=args.root_stages.split(",") if args.root_stages else (),
        owner_inputs=owner_inputs,
        retrieved_at=args.retrieved_at,
        workers=args.workers,
        screen_max_bp=args.screen_max_bp,
        work=Path(args.work) if args.work else None,
        threads=args.threads,
        memory_limit=args.memory_limit,
    )

    if args.gtdb_tar and args.globdb_tax:
        gtdb_accs = load_gtdb_base_accessions(Path(args.gtdb_tar))
        globdb_accs = load_globdb_base_accessions(
            Path(args.globdb_tax), Path(args.globdb_dict_tar) if args.globdb_dict_tar else None
        )
        iuig_accs = load_iuig_base_accessions(Path(args.iuig_zip)) if args.iuig_zip else set()
        overlap = {}
        if args.gorg_ena_tsv and Path(args.gorg_ena_tsv).is_file():
            overlap["gorg_tropics_PRJEB33281"] = compute_gorg_metadata_overlap(
                Path(args.gorg_ena_tsv), gtdb_accs=gtdb_accs, globdb_accs=globdb_accs, iuig_accs=iuig_accs
            )
        if args.gorg_dark_ena_tsv and Path(args.gorg_dark_ena_tsv).is_file():
            overlap["gorg_dark_PRJNA1259050"] = compute_gorg_metadata_overlap(
                Path(args.gorg_dark_ena_tsv), gtdb_accs=gtdb_accs, globdb_accs=globdb_accs, iuig_accs=iuig_accs
            )
        rep["metadata_overlap"] = overlap

    if args.report_out:
        Path(args.report_out).write_text(json.dumps(rep, indent=2))
    print(json.dumps(rep, indent=2))


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Card C-58: GORG SAGs, endosymbionts, post-R226, niche compendia")
    sub = p.add_subparsers(dest="cmd", required=True)

    p_ft = sub.add_parser("fetch-tables")
    p_ft.add_argument("out_dir")
    p_ft.set_defaults(func=cmd_fetch_tables)

    p_dn = sub.add_parser("diff-ncbi")
    p_dn.add_argument("tables_dir")
    p_dn.add_argument("gtdb_tar")
    p_dn.add_argument("globdb_tax")
    p_dn.add_argument("globdb_dict_tar")
    p_dn.add_argument("out_dir")
    p_dn.set_defaults(func=cmd_diff_ncbi)

    p_fs = sub.add_parser("fetch-sym1")
    p_fs.add_argument("candidates_tsv")
    p_fs.add_argument("out_dir")
    p_fs.add_argument("--max-genomes", type=int, default=0)
    p_fs.set_defaults(func=cmd_fetch_sym1)

    p_ss = sub.add_parser("stage-sym1")
    p_ss.add_argument("manifest_tsv")
    p_ss.add_argument("out_root")
    p_ss.add_argument("--root", default=None)
    p_ss.add_argument("--root-stages", default="")
    p_ss.add_argument("--owner-inputs", action="append", default=[])
    p_ss.add_argument("--retrieved-at", default="2026-11-24T00:00:00Z")
    p_ss.add_argument("--workers", type=int, default=2)
    p_ss.add_argument("--work", default=None)
    p_ss.add_argument("--threads", type=int, default=2)
    p_ss.add_argument("--memory-limit", default="4GB")
    p_ss.add_argument("--report-out", default=None)
    p_ss.set_defaults(func=cmd_stage_sym1)

    p_fg = sub.add_parser("fetch-gorg")
    p_fg.add_argument("out_dir")
    p_fg.set_defaults(func=cmd_fetch_gorg)

    p_sg = sub.add_parser("stage-sag1")
    p_sg.add_argument("gorg_fasta_gz")
    p_sg.add_argument("out_root")
    p_sg.add_argument("--fasta-dir", default=None)
    p_sg.add_argument("--final-dir", default=None)
    p_sg.add_argument("--gorg-ena-tsv", default=None)
    p_sg.add_argument("--gorg-dark-ena-tsv", default=None)
    p_sg.add_argument("--gtdb-tar", default=None)
    p_sg.add_argument("--globdb-tax", default=None)
    p_sg.add_argument("--globdb-dict-tar", default=None)
    p_sg.add_argument("--iuig-zip", default=None)
    p_sg.add_argument("--root", default=None)
    p_sg.add_argument("--root-stages", default="")
    p_sg.add_argument("--owner-inputs", action="append", default=[])
    p_sg.add_argument("--retrieved-at", default="2026-11-24T00:00:00Z")
    p_sg.add_argument("--workers", type=int, default=2)
    p_sg.add_argument("--screen-max-bp", type=int, default=500_000)
    p_sg.add_argument("--work", default=None)
    p_sg.add_argument("--threads", type=int, default=2)
    p_sg.add_argument("--memory-limit", default="4GB")
    p_sg.add_argument("--report-out", default=None)
    p_sg.set_defaults(func=cmd_stage_sag1)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
