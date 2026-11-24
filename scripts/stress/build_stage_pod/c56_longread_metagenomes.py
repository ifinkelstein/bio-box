"""Card C-56 (DEC-69, DEC-70, DEC-71): prioritize, fetch, and stage public long-read (PacBio HiFi / ONT)
metagenome assemblies (tranche LR1..).

Subcommands:
    python c56_longread_metagenomes.py crosswalk \
        --read-run-lr PATH --ena-analysis PATH --ena-wgs PATH --ncbi-summary PATH \
        [--img-m PATH] --out-dir DIR [--max-bases 35000000000] [--seed 20261123]

    python c56_longread_metagenomes.py fetch-and-stage \
        --selection LR1_selection.tsv --out-fasta-dir STOR_DIR --out-meta-dir META_DIR \
        --tmp-dl TMP_DL --target-bases 25000000000 --max-contigs 3000000 \
        [--min-contig-len 5000] [--interval 0.6] [--root ROOT] [--root-stages ...] \
        [--owner-inputs NAME=DIR ...]

    python c56_longread_metagenomes.py stage \
        --manifest MANIFEST.json --out-fasta-dir FASTA_DIR --out-meta-dir META_DIR \
        [--final-fasta-dir DIR] [--min-contig-len 5000] [--root ROOT] \
        [--root-stages ...] [--owner-inputs NAME=DIR ...]
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from collections.abc import Iterable
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parent))

import c48c_virus_inputs as v2mod  # noqa: E402
from biobox import fasta, ids  # noqa: E402
from biobox.l0seq import store  # noqa: E402

DEFAULT_MIN_CONTIG_LEN = 5_000
DEFAULT_TARGET_BASES = 25_000_000_000
DEFAULT_MAX_CONTIGS = 3_000_000
DEFAULT_SEED = 20261123
ROW_GROUP = 122_880

LR1_FASTA_NAME = "lr1_longread_ge5kb.fna"
INSDC_LICENSE = (
    "INSDC / ENA / NCBI GenBank public metagenome assemblies "
    "(EMBL-EBI terms of use; DEC-50, DEC-51, DEC-69, DEC-70, DEC-71)"
)

# DEC-71 habitat priority weights (soil/terrestrial, wastewater/sludge, freshwater/cyanobacterial lead)
HABITAT_SPECS: dict[str, tuple[int, float]] = {
    "soil_terrestrial": (1, 3.0),
    "wastewater_sludge": (2, 2.5),
    "freshwater_aquatic": (3, 2.5),
    "marine_sediment": (4, 1.2),
    "host_associated": (5, 0.8),
    "other_unlabelled": (6, 0.4),
}

PLATFORM_RANK: dict[str, int] = {
    "pacbio_ont": 0,
    "pacbio": 1,
    "ont": 2,
    "hybrid_both": 3,
    "hybrid_pacbio": 4,
    "hybrid_ont": 5,
    "short_only": 9,
    "unknown": 10,
}

_CIRCULAR_RE = re.compile(
    r"(?:(?:^|[\s|;,])(?:circular\s*=\s*(?:true|yes|y|1)|circ\s*=\s*(?:true|yes|y|1)|"
    r"is_circular\s*=\s*(?:true|yes|y|1)|topology\s*=\s*circular|circular_contig|dtr\b))",
    re.IGNORECASE,
)

_CONTIG_SCHEMA = pa.schema(
    [
        ("contig_uid", pa.binary(16)),
        ("contig_len", pa.int64()),
        ("circular", pa.bool_()),
        ("asm_idx", pa.int32()),
    ]
)


def _sql_lit(path: Path | str) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def _sha256_str(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def _md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


# --- 1. Platform and habitat classification -------------------------------------------------------------------------

_HIFI_RE = re.compile(r"\b(revio|sequel\s*ii[e]?|vega|hifi)\b", re.IGNORECASE)


def classify_platform(
    platforms: Iterable[str],
    models: Iterable[str] = (),
) -> tuple[str, bool]:
    """Classify ENA/NCBI instrument platforms + models into (platform_class, is_hifi)."""
    has_pb = False
    has_ont = False
    has_short = False
    is_hifi = False

    for p in platforms:
        pu = (p or "").strip().upper()
        if not pu:
            continue
        if "PACBIO" in pu or "SMRT" in pu:
            has_pb = True
        elif "OXFORD" in pu or "NANOPORE" in pu or pu == "ONT":
            has_ont = True
        elif any(k in pu for k in ("ILLUMINA", "BGISEQ", "DNBSEQ", "ION_TORRENT", "LS454", "ELEMENT", "ULTIMA")):
            has_short = True

    for m in models:
        ms = (m or "").strip()
        if not ms:
            continue
        if _HIFI_RE.search(ms):
            has_pb = True
            is_hifi = True
        ml = ms.lower()
        if any(k in ml for k in ("minion", "promethion", "gridion", "flongle", "nanopore")):
            has_ont = True
        elif any(k in ml for k in ("pacbio", "sequel", "rs ii", "rsii")):
            has_pb = True
        elif any(k in ml for k in ("novaseq", "hiseq", "miseq", "nextseq", "illumina", "dnbseq", "bgiseq", "miseq")):
            has_short = True

    if has_pb and has_ont:
        cls = "hybrid_both" if has_short else "pacbio_ont"
    elif has_pb:
        cls = "hybrid_pacbio" if has_short else "pacbio"
    elif has_ont:
        cls = "hybrid_ont" if has_short else "ont"
    elif has_short:
        cls = "short_only"
    else:
        cls = "unknown"
    return cls, is_hifi


def classify_imgm_method(method: str) -> tuple[str, bool]:
    """Classify IMG/M `Sequencing Method` string into (platform_class, is_hifi)."""
    parts = [p.strip() for p in (method or "").split(",") if p.strip()]
    return classify_platform(parts, parts)


def classify_habitat(
    organism: str,
    ecosystem: str = "",
    study_title: str = "",
) -> tuple[str, int]:
    """Map organism/ecosystem/study text to (habitat_stratum, dec71_priority)."""
    text = f"{organism or ''} {ecosystem or ''} {study_title or ''}".lower()

    # Priority 2: wastewater / activated sludge / bioreactor / anaerobic digester
    if any(
        k in text
        for k in (
            "activated sludge",
            "wastewater",
            "sewage",
            "anaerobic digester",
            "bioreactor",
            "sludge",
            "compost",
            "landfill",
            "leachate",
            "biofilm",
            "ferment",
            "mine drainage",
            "tailings",
            "oil sands",
            "drinking water",
        )
    ):
        return "wastewater_sludge", HABITAT_SPECS["wastewater_sludge"][0]

    # Priority 1: soil / terrestrial / rhizosphere / permafrost / peat / subsurface
    if any(
        k in text
        for k in (
            "soil",
            "permafrost",
            "peat",
            "rhizosphere",
            "terrestrial",
            "forest",
            "grassland",
            "tundra",
            "subsurface",
            "desert",
            "bog",
            "fen ",
            "prairie",
            "savanna",
            "myxococc",
            "rock",
            "cave",
            "sand",
            "loam",
            "clay",
        )
    ):
        return "soil_terrestrial", HABITAT_SPECS["soil_terrestrial"][0]

    # Priority 3: freshwater / aquatic / cyanobacterial mats / hot springs / wetland
    if any(
        k in text
        for k in (
            "freshwater",
            "lake",
            "river",
            "stream",
            "pond",
            "wetland",
            "groundwater",
            "aquifer",
            "hot spring",
            "thermal spring",
            "spring",
            "cyanobacter",
            "microbial mat",
            "estuar",
            "brackish",
            "glacier",
            "ice ",
            "snow",
            "lagoon",
            "hypersaline",
            "soda lake",
            "salt lake",
            "aquatic",
        )
    ):
        return "freshwater_aquatic", HABITAT_SPECS["freshwater_aquatic"][0]

    # Priority 4: marine / ocean / sediment / hydrothermal vent
    if any(
        k in text
        for k in (
            "marine",
            "ocean",
            "seawater",
            "pelagic",
            "hydrothermal",
            "cold seep",
            "sediment",
            "benthic",
            "coral",
            "sponge",
            "reef",
            "coastal",
            "deep sea",
            "estuarine",
        )
    ):
        return "marine_sediment", HABITAT_SPECS["marine_sediment"][0]

    # Priority 5: host-associated (gut, rumen, feces, oral, skin, plant, insect, animal)
    if any(
        k in text
        for k in (
            "gut",
            "feces",
            "faeces",
            "fecal",
            "stool",
            "rumen",
            "oral",
            "saliva",
            "skin",
            "vaginal",
            "lung",
            "sputum",
            "respiratory",
            "human",
            "mouse",
            "bovine",
            "pig",
            "swine",
            "chicken",
            "poultry",
            "fish",
            "insect",
            "termite",
            "plant",
            "leaf",
            "phyllosphere",
            "root ",
            "lichen",
            "moss",
            "algae",
            "algal",
            "fungus",
            "animal",
            "host-associated",
            "host associated",
            "milk",
            "blood",
            "urinary",
            "nasal",
            "sheep",
            "goat",
            "canine",
            "feline",
            "equine",
            "bird",
            "invertebrate",
            "mollusc",
            "oyster",
            "crustacean",
            "symbiont",
        )
    ):
        return "host_associated", HABITAT_SPECS["host_associated"][0]

    return "other_unlabelled", HABITAT_SPECS["other_unlabelled"][0]


# --- 2. Crosswalk and deduplication ---------------------------------------------------------------------------------

def _winner_key(row: dict) -> tuple:
    """Deterministic sort key: lower tuple wins."""
    access_rank = 0 if row.get("access") == "open" else 1
    is_emg = 1 if (row.get("center_name") or "").strip().upper() == "EMG" else 0
    p_rank = PLATFORM_RANK.get(row.get("platform_class", "unknown"), 10)
    hifi_rank = 0 if row.get("is_hifi") else 1
    # Prefer submitter long-read assemblies over EMG re-assemblies, pure long-read before hybrid, HiFi first
    pure_rank = 0 if p_rank <= 2 else 1
    est_bp = int(row.get("est_bp") or 0)
    src_rank = {"ncbi_wgs": 0, "ena_analysis": 1, "ena_wgs": 2, "img_m": 3}.get(row.get("source", ""), 4)
    return (access_rank, is_emg, pure_rank, hifi_rank, p_rank, -est_bp, src_rank, row["accession"])


def deduplicate_assemblies(candidates: list[dict]) -> tuple[list[dict], list[dict]]:
    """Deduplicate assemblies sharing a biosample or run accession to one best winner per group.

    Order-invariant: depends only on the set of candidates, never input row order.
    Returns (winners, duplicates).
    """
    if not candidates:
        return [], []

    # Sort canonically first so union-find indexing is deterministic
    items = sorted(candidates, key=lambda r: (r["accession"], r.get("source", "")))
    n = len(items)
    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri != rj:
            if ri < rj:
                parent[rj] = ri
            else:
                parent[ri] = rj

    by_sample: dict[str, int] = {}
    by_run: dict[str, int] = {}
    for idx, row in enumerate(items):
        bs = (row.get("biosample") or "").strip()
        if bs:
            if bs in by_sample:
                union(idx, by_sample[bs])
            else:
                by_sample[bs] = idx
        sec = (row.get("secondary_sample") or "").strip()
        if sec:
            if sec in by_sample:
                union(idx, by_sample[sec])
            else:
                by_sample[sec] = idx
        for run in row.get("runs") or ():
            r = run.strip()
            if not r:
                continue
            if r in by_run:
                union(idx, by_run[r])
            else:
                by_run[r] = idx

    groups: dict[int, list[dict]] = defaultdict(list)
    for idx, row in enumerate(items):
        groups[find(idx)].append(row)

    winners: list[dict] = []
    duplicates: list[dict] = []
    for grp in groups.values():
        grp.sort(key=_winner_key)
        win = dict(grp[0])
        winners.append(win)
        for dup in grp[1:]:
            d = dict(dup)
            d["winner_accession"] = win["accession"]
            d["winner_source"] = win.get("source", "")
            duplicates.append(d)

    winners.sort(key=lambda r: r["accession"])
    duplicates.sort(key=lambda r: r["accession"])
    return winners, duplicates


def select_tranche(
    rows: list[dict],
    *,
    max_bases: int = DEFAULT_TARGET_BASES,
    max_single_bp: int | None = None,
    seed: int = DEFAULT_SEED,
) -> tuple[list[dict], list[dict]]:
    """Order open-access assemblies across habitat strata weighted for DEC-71 ART-RT priority.

    Submitter-deposited long-read assemblies (center_name != 'EMG') are prioritized ahead of
    MGnify/EMG re-assemblies, which often re-assemble companion short-read runs.
    Within each habitat stratum, studies are interleaved round-robin in sha256(seed:study) order,
    preferring HiFi / pure long-read submitter assemblies.
    Across strata, Webster's divisor method with DEC-71 habitat weights picks the next stratum,
    breaking ties by habitat_priority ascending.
    Returns (selected, remainder), each carrying contiguous 1-based `lr_rank`.
    """
    open_rows = [dict(r) for r in rows if r.get("access", "open") == "open"]
    if not open_rows:
        return [], []

    def _order_group(group_rows: list[dict]) -> list[dict]:
        by_hab: dict[str, dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
        for r in group_rows:
            hab = r.get("habitat") or "other_unlabelled"
            st = r.get("study") or r["accession"]
            by_hab[hab][st].append(r)

        queues: dict[str, list[dict]] = {}
        for hab, studies in by_hab.items():
            # Studies with HiFi or NCBI WGS assemblies first, then deterministic hash
            study_order = sorted(
                studies.keys(),
                key=lambda s: (
                    0 if any(r.get("is_hifi") for r in studies[s]) else 1,
                    0 if any(r.get("source") == "ncbi_wgs" for r in studies[s]) else 1,
                    _sha256_str(f"{seed}:study:{hab}:{s}"),
                ),
            )
            study_queues = {}
            for st in study_order:
                items = sorted(
                    studies[st],
                    key=lambda r: (
                        0 if r.get("is_hifi") else 1,
                        0 if r.get("source") == "ncbi_wgs" else 1,
                        PLATFORM_RANK.get(r.get("platform_class", "unknown"), 10),
                        _sha256_str(f"{seed}:acc:{r['accession']}"),
                    ),
                )
                study_queues[st] = items
            q: list[dict] = []
            idx = 0
            while True:
                added = False
                for st in study_order:
                    if idx < len(study_queues[st]):
                        q.append(study_queues[st][idx])
                        added = True
                if not added:
                    break
                idx += 1
            queues[hab] = q

        taken_per_hab: Counter[str] = Counter()
        pos_per_hab: Counter[str] = Counter()
        out: list[dict] = []
        total = len(group_rows)
        while len(out) < total:
            best_hab = None
            best_score = None
            for hab, q in queues.items():
                if pos_per_hab[hab] >= len(q):
                    continue
                prio, weight = HABITAT_SPECS.get(hab, (6, 0.4))
                score = (weight / (2 * taken_per_hab[hab] + 1), -prio, hab)
                if best_score is None or score > best_score:
                    best_score = score
                    best_hab = hab
            assert best_hab is not None
            item = queues[best_hab][pos_per_hab[best_hab]]
            pos_per_hab[best_hab] += 1
            taken_per_hab[best_hab] += 1
            out.append(item)
        return out

    submitter_fit = []
    submitter_oversize = []
    emg_rows = []
    for r in open_rows:
        is_emg = (r.get("center_name") or "").strip().upper() == "EMG"
        bp = int(r.get("est_bp") or 0)
        if is_emg:
            emg_rows.append(r)
        elif max_single_bp is not None and bp > max_single_bp:
            submitter_oversize.append(r)
        else:
            submitter_fit.append(r)

    ordered = _order_group(submitter_fit) + _order_group(submitter_oversize) + _order_group(emg_rows)
    for i, item in enumerate(ordered, start=1):
        item["lr_rank"] = i

    selected: list[dict] = []
    remainder: list[dict] = []
    cum_bp = 0
    stopped = False
    for item in ordered:
        bp = int(item.get("est_bp") or 0)
        if not stopped and (cum_bp + bp <= max_bases or (not selected and bp <= max_bases)):
            selected.append(item)
            cum_bp += bp
        else:
            stopped = True
            remainder.append(item)
    return selected, remainder


def _parse_ena_ftp_url(raw_ftp: str) -> str:
    """Pick the primary FASTA URL from ENA's semicolon-delimited ftp column."""
    if not raw_ftp:
        return ""
    parts = [p.strip() for p in raw_ftp.split(";") if p.strip()]
    # Prefer .fa.gz / .fna.gz / .fasta.gz
    chosen = ""
    for p in parts:
        pl = p.lower()
        if pl.endswith((".fa.gz", ".fna.gz", ".fasta.gz", ".fsa.gz")):
            chosen = p
            break
    if not chosen and parts:
        chosen = parts[0]
    if not chosen:
        return ""
    if chosen.startswith(("http://", "https://")):
        return chosen
    return "https://" + chosen.lstrip("/")


def _parse_ena_md5(raw_ftp: str, raw_md5: str) -> str:
    if not raw_ftp or not raw_md5:
        return ""
    ftps = [p.strip() for p in raw_ftp.split(";") if p.strip()]
    md5s = [m.strip() for m in raw_md5.split(";") if m.strip()]
    if len(ftps) == len(md5s):
        for p, m in zip(ftps, md5s):
            if p.lower().endswith((".fa.gz", ".fna.gz", ".fasta.gz", ".fsa.gz")):
                return m
    return md5s[0] if md5s else ""


def build_crosswalk(
    *,
    read_run_lr: Path,
    ena_analysis: Path,
    ena_wgs: Path,
    ncbi_summary: Path,
    img_m: Path | None,
    out_dir: Path,
    max_bases: int = DEFAULT_TARGET_BASES,
    max_single_bp: int | None = None,
    seed: int = DEFAULT_SEED,
) -> dict:
    """Join ENA long-read runs, ENA primary analyses, ENA WGS sets, NCBI assembly summary, and IMG/M."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Index all long-read metagenomic runs by run_accession and sample_accession / secondary_sample_accession
    run_info: dict[str, tuple[str, str, str, str, int, str]] = {}
    sample_runs: dict[str, list[str]] = defaultdict(list)
    with open(read_run_lr, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            run = (row.get("run_accession") or "").strip()
            if not run:
                continue
            plat = (row.get("instrument_platform") or "").strip()
            model = (row.get("instrument_model") or "").strip()
            sa = (row.get("sample_accession") or "").strip()
            ssa = (row.get("secondary_sample_accession") or "").strip()
            bc = int(row.get("base_count") or 0)
            sci = (row.get("scientific_name") or "").strip()
            run_info[run] = (plat, model, sa, ssa, bc, sci)
            if sa:
                sample_runs[sa].append(run)
            if ssa:
                sample_runs[ssa].append(run)

    def lr_from_runs_or_sample(runs_str: str, sample_acc: str, sec_sample: str = "") -> tuple[list[str], str, bool, int, str]:
        matched_runs: list[str] = []
        seen: set[str] = set()
        for tok in re.split(r"[;,\s]+", runs_str or ""):
            t = tok.strip()
            if t and t in run_info and t not in seen:
                seen.add(t)
                matched_runs.append(t)
        if not matched_runs:
            for s in (sample_acc, sec_sample):
                if s and s in sample_runs:
                    for t in sample_runs[s]:
                        if t not in seen:
                            seen.add(t)
                            matched_runs.append(t)
        if not matched_runs:
            return [], "unknown", False, 0, ""
        plats = [run_info[r][0] for r in matched_runs]
        models = [run_info[r][1] for r in matched_runs]
        p_cls, is_hifi = classify_platform(plats, models)
        read_bp = sum(run_info[r][4] for r in matched_runs)
        sci = next((run_info[r][5] for r in matched_runs if run_info[r][5]), "")
        return sorted(matched_runs), p_cls, is_hifi, read_bp, sci

    candidates: list[dict] = []

    # 2. ENA analysis primary metagenomes (ERZ*)
    ena_analysis_matched = 0
    with open(ena_analysis, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = (row.get("analysis_accession") or "").strip()
            if not acc:
                continue
            atype = (row.get("assembly_type") or "").strip().lower()
            if atype and atype != "primary metagenome":
                continue
            sa = (row.get("sample_accession") or "").strip()
            runs_str = (row.get("run_accession") or "").strip()
            m_runs, p_cls, is_hifi, read_bp, run_sci = lr_from_runs_or_sample(runs_str, sa)
            if not m_runs or p_cls in ("short_only", "unknown"):
                continue
            # Prefer generated_ftp (standardized contig.fa.gz) if present, else submitted_ftp
            gen_ftp = _parse_ena_ftp_url(row.get("generated_ftp") or "")
            sub_ftp = _parse_ena_ftp_url(row.get("submitted_ftp") or "")
            if gen_ftp:
                url = gen_ftp
                md5 = _parse_ena_md5(row.get("generated_ftp") or "", row.get("generated_md5") or "")
                gz_bytes = int((row.get("generated_bytes") or "0").split(";")[0] or 0)
            elif sub_ftp:
                url = sub_ftp
                md5 = _parse_ena_md5(row.get("submitted_ftp") or "", row.get("submitted_md5") or "")
                gz_bytes = int((row.get("submitted_bytes") or "0").split(";")[0] or 0)
            else:
                continue
            if gz_bytes <= 0:
                continue
            # Estimate assembled bases from compressed contig FASTA bytes (~3.4x ratio)
            est_bp = int(gz_bytes * 3.4)
            sci = (row.get("scientific_name") or "").strip() or run_sci
            hab, prio = classify_habitat(sci)
            ena_analysis_matched += 1
            candidates.append(
                {
                    "accession": acc,
                    "source": "ena_analysis",
                    "access": "open",
                    "center_name": (row.get("center_name") or "").strip(),
                    "assembly_type": "primary metagenome",
                    "biosample": sa,
                    "secondary_sample": "",
                    "runs": m_runs,
                    "study": (row.get("study_accession") or "").strip(),
                    "platform_class": p_cls,
                    "is_hifi": is_hifi,
                    "read_bp": read_bp,
                    "gz_bytes": gz_bytes,
                    "est_bp": est_bp,
                    "ftp_url": url,
                    "md5": md5,
                    "scientific_name": sci,
                    "ecosystem": "",
                    "habitat": hab,
                    "habitat_priority": prio,
                }
            )

    # 3. ENA WGS sets + NCBI metagenome assembly summary (excluding MAGs!)
    # First read ENA wgs_set to know which GCA accessions are MAGs vs primary metagenomes
    ena_mag_gcas: set[str] = set()
    wgs_primary_by_gca: dict[str, dict] = {}
    with open(ena_wgs, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            atype = (row.get("assembly_type") or "").strip().lower()
            gca = (row.get("assembly_accession") or "").strip()
            gca_base = gca.split(".")[0] if gca else ""
            if atype not in ("", "primary metagenome"):
                if gca_base:
                    ena_mag_gcas.add(gca_base)
                continue
            if gca_base:
                wgs_primary_by_gca[gca_base] = row

    ncbi_matched = 0
    with open(ncbi_summary, newline="") as fh:
        header: list[str] | None = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line[1:].rstrip("\n").split("\t")
                continue
            if not header:
                continue
            parts = line.rstrip("\n").split("\t")
            row = dict(zip(header, parts))
            acc = (row.get("assembly_accession") or "").strip()
            if not acc or row.get("version_status") != "latest":
                continue
            acc_base = acc.split(".")[0]
            if acc_base in ena_mag_gcas:
                continue
            # Exclude obvious MAGs from NCBI excluded_from_refseq or asm_name
            excl = (row.get("excluded_from_refseq") or "").lower()
            if "derived from metagenome" in excl or "metagenome-assembled" in excl:
                continue
            bs = (row.get("biosample") or "").strip()
            wgs_row = wgs_primary_by_gca.get(acc_base, {})
            runs_str = (wgs_row.get("run_accession") or "").strip()
            sec_sa = (wgs_row.get("secondary_sample_accession") or "").strip()
            m_runs, p_cls, is_hifi, read_bp, run_sci = lr_from_runs_or_sample(runs_str, bs, sec_sa)
            if not m_runs or p_cls in ("short_only", "unknown"):
                continue
            ftp_path = (row.get("ftp_path") or "").strip()
            if not ftp_path or ftp_path == "na":
                continue
            asm_dir = ftp_path.rstrip("/").rsplit("/", 1)[-1]
            url = ftp_path.rstrip("/") + f"/{asm_dir}_genomic.fna.gz"
            if url.startswith("ftp://"):
                url = "https://" + url[6:]
            est_bp = int(row.get("genome_size") or wgs_row.get("base_count") or 0)
            if est_bp < 50_000:
                continue
            sci = (row.get("organism_name") or "").strip() or run_sci
            hab, prio = classify_habitat(sci)
            ncbi_matched += 1
            candidates.append(
                {
                    "accession": acc,
                    "source": "ncbi_wgs",
                    "access": "open",
                    "center_name": (row.get("asm_submitter") or "").strip(),
                    "assembly_type": "primary metagenome",
                    "biosample": bs,
                    "secondary_sample": sec_sa,
                    "runs": m_runs,
                    "study": (row.get("bioproject") or wgs_row.get("study_accession") or "").strip(),
                    "platform_class": p_cls,
                    "is_hifi": is_hifi,
                    "read_bp": read_bp,
                    "gz_bytes": 0,
                    "est_bp": est_bp,
                    "ftp_url": url,
                    "md5": "",
                    "scientific_name": sci,
                    "ecosystem": "",
                    "habitat": hab,
                    "habitat_priority": prio,
                }
            )

    # 4. IMG/M unrestricted long-read metagenomes (JGI login required; inventory + JGI batch only)
    imgm_matched = 0
    if img_m and Path(img_m).is_file():
        opener = gzip.open if str(img_m).endswith(".gz") else open
        with opener(img_m, "rt", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                method = (row.get("Sequencing Method") or "").strip()
                p_cls, is_hifi = classify_imgm_method(method)
                if p_cls in ("short_only", "unknown"):
                    continue
                status = (row.get("JGI Data Utilization Status") or "").strip().lower()
                if status not in ("", "unrestricted"):
                    continue
                oid = (row.get("IMG Genome ID") or "").strip()
                if not oid:
                    continue
                bs = (row.get("NCBI Biosample Accession") or "").strip()
                sra_run = (row.get("SRA Run") or "").strip()
                runs = [r.strip() for r in re.split(r"[;,\s]+", sra_run) if r.strip()]
                eco = ";".join(
                    filter(
                        None,
                        [
                            (row.get("Ecosystem") or "").strip(),
                            (row.get("Ecosystem Category") or "").strip(),
                            (row.get("Ecosystem Type") or "").strip(),
                            (row.get("Ecosystem Subtype") or "").strip(),
                            (row.get("Specific Ecosystem") or "").strip(),
                        ],
                    )
                )
                name = (row.get("Genome Name / Sample Name") or "").strip()
                hab, prio = classify_habitat(name, eco, (row.get("Study Name") or "").strip())
                est_bp = int(row.get("Genome Size * assembled") or 0)
                imgm_matched += 1
                candidates.append(
                    {
                        "accession": oid,
                        "source": "img_m",
                        "access": "jgi_login",
                        "center_name": (row.get("Sequencing Center") or "").strip(),
                        "assembly_type": "primary metagenome",
                        "biosample": bs,
                        "secondary_sample": "",
                        "runs": runs,
                        "study": (row.get("GOLD Study ID") or row.get("NCBI Bioproject Accession") or "").strip(),
                        "platform_class": p_cls,
                        "is_hifi": is_hifi,
                        "read_bp": 0,
                        "gz_bytes": 0,
                        "est_bp": est_bp,
                        "ftp_url": "",
                        "md5": "",
                        "scientific_name": name,
                        "ecosystem": eco,
                        "habitat": hab,
                        "habitat_priority": prio,
                    }
                )

    winners, duplicates = deduplicate_assemblies(candidates)
    open_winners = [w for w in winners if w["access"] == "open"]
    jgi_winners = [w for w in winners if w["access"] == "jgi_login"]

    selected, remainder = select_tranche(open_winners, max_bases=max_bases, max_single_bp=max_single_bp, seed=seed)

    cols = [
        "lr_rank", "accession", "source", "access", "platform_class", "is_hifi",
        "habitat", "habitat_priority", "est_bp", "gz_bytes", "read_bp",
        "biosample", "study", "runs", "center_name", "scientific_name", "ecosystem",
        "ftp_url", "md5",
    ]

    def write_tsv(path: Path, rows_list: list[dict], extra_cols: list[str] = ()) -> None:
        all_cols = list(cols) + list(extra_cols)
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=all_cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            w.writeheader()
            for r in rows_list:
                out = dict(r)
                out["runs"] = ",".join(r.get("runs") or ())
                w.writerow(out)

    write_tsv(out_dir / "inventory_winners.tsv", winners)
    write_tsv(out_dir / "inventory_duplicates.tsv", duplicates, ["winner_accession", "winner_source"])
    write_tsv(out_dir / "LR1_selection.tsv", selected)
    write_tsv(out_dir / "LR1_remainder.tsv", remainder)
    write_tsv(out_dir / "LR_jgi_batch.tsv", jgi_winners)

    def breakdown(rows_list: list[dict], key: str) -> dict[str, dict[str, int]]:
        g: dict[str, dict[str, int]] = {}
        for r in rows_list:
            k = str(r.get(key, ""))
            d = g.setdefault(k, {"assemblies": 0, "est_bp": 0})
            d["assemblies"] += 1
            d["est_bp"] += int(r.get("est_bp") or 0)
        return dict(sorted(g.items()))

    summary = {
        "lr_runs_indexed": len(run_info),
        "raw_candidates": len(candidates),
        "raw_by_source": {"ena_analysis": ena_analysis_matched, "ncbi_wgs": ncbi_matched, "img_m": imgm_matched},
        "dedup_winners": len(winners),
        "dedup_duplicates": len(duplicates),
        "open_winners": len(open_winners),
        "open_winners_est_bp": sum(int(r.get("est_bp") or 0) for r in open_winners),
        "jgi_only_winners": len(jgi_winners),
        "jgi_only_est_bp": sum(int(r.get("est_bp") or 0) for r in jgi_winners),
        "lr1_selected": len(selected),
        "lr1_selected_est_bp": sum(int(r.get("est_bp") or 0) for r in selected),
        "lr1_remainder": len(remainder),
        "lr1_remainder_est_bp": sum(int(r.get("est_bp") or 0) for r in remainder),
        "open_by_habitat": breakdown(open_winners, "habitat"),
        "open_by_platform": breakdown(open_winners, "platform_class"),
        "open_by_source": breakdown(open_winners, "source"),
        "lr1_by_habitat": breakdown(selected, "habitat"),
        "lr1_by_platform": breakdown(selected, "platform_class"),
        "jgi_by_habitat": breakdown(jgi_winners, "habitat"),
        "jgi_by_platform": breakdown(jgi_winners, "platform_class"),
    }
    (out_dir / "crosswalk_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


# --- 3. DEC-65 circular/DTR handling and staging --------------------------------------------------------------------

def is_header_circular(header: str) -> bool:
    """Return True if a long-read assembler FASTA header explicitly marks a circular contig."""
    if not header:
        return False
    if _CIRCULAR_RE.search(header):
        return True
    # hifiasm-meta uses trailing 'c' on ptg/s/u contig IDs (e.g. `s123.ctg000001c`)
    first = header.split()[0]
    if re.match(r"^(?:ptg|s\d+\.ctg|u\d+\.ctg)\d+c$", first):
        return True
    return False


def apply_dec65_longread(header: str, seq: bytes) -> tuple[bytes, bool, str, int]:
    """Apply DEC-65 terminal-repeat rules to a long-read metagenome contig.

    Returns (staged_seq, is_circular, dtr_class, trimmed_bases).
    """
    if not is_header_circular(header):
        return seq, False, "", 0
    e = v2mod.ends(seq)
    dtr_cls = v2mod.classify(e)
    if dtr_cls == "exact":
        return v2mod.trim(seq, e.exact), True, dtr_cls, e.exact
    if dtr_cls == "none":
        return seq, True, dtr_cls, 0
    return seq, False, dtr_cls, 0


def _load_earlier_owners(
    root: Path | None,
    root_stages: list[str],
    owner_inputs: list[tuple[str, Path]],
    scratch: Path,
) -> Path | None:
    """Collect (contig_uid, owner_stage) into a DuckDB-readable Parquet file in scratch."""
    parts: list[tuple[str, list[Path]]] = []
    if root and root_stages:
        for st in root_stages:
            st_dir = Path(root) / "l0seq" / "contig_index" / f"stage={st}"
            pqs = sorted(st_dir.glob("*.parquet"))
            if pqs:
                parts.append((st, pqs))
    if not parts and not owner_inputs:
        return None

    owners_pq = scratch / "earlier_owners.parquet"
    schema = pa.schema([("contig_uid", pa.binary(16)), ("owner", pa.string())])
    writer = pq.ParquetWriter(owners_pq, schema, compression="zstd", write_statistics=False)

    # Track seen uids via DuckDB to keep first owner in build order
    con = duckdb.connect(
        config={"threads": 1, "memory_limit": "2GB", "temp_directory": str(scratch / "duck_own")}
    )
    try:
        for st, pqs in parts:
            lit = "[" + ", ".join(_sql_lit(p) for p in pqs) + "]"
            reader = con.sql(
                f"SELECT DISTINCT contig_uid FROM read_parquet({lit})"
            ).to_arrow_reader(ROW_GROUP)
            for b in reader:
                cu = b.column(0)
                writer.write_batch(
                    pa.RecordBatch.from_arrays(
                        [cu, pa.array([st] * len(cu), type=pa.string())], schema=schema
                    )
                )
        for st, inp_dir in owner_inputs:
            rec_tsv = Path(inp_dir) / "records.tsv"
            gen_json = Path(inp_dir) / "genomes.json"
            if rec_tsv.is_file():
                buf: list[bytes] = []
                with open(rec_tsv, newline="") as fh:
                    for row in csv.DictReader(fh, delimiter="\t"):
                        u = row.get("contig_uid") or ""
                        if len(u) == 32:
                            buf.append(bytes.fromhex(u))
                            if len(buf) >= ROW_GROUP:
                                writer.write_batch(
                                    pa.RecordBatch.from_arrays(
                                        [pa.array(buf, pa.binary(16)), pa.array([st] * len(buf), pa.string())],
                                        schema=schema,
                                    )
                                )
                                buf.clear()
                if buf:
                    writer.write_batch(
                        pa.RecordBatch.from_arrays(
                            [pa.array(buf, pa.binary(16)), pa.array([st] * len(buf), pa.string())],
                            schema=schema,
                        )
                    )
            elif gen_json.is_file():
                data = json.loads(gen_json.read_text())
                buf = []
                for g in data:
                    for u in g.get("contig_uids") or ():
                        buf.append(bytes.fromhex(u))
                        if len(buf) >= ROW_GROUP:
                            writer.write_batch(
                                pa.RecordBatch.from_arrays(
                                    [pa.array(buf, pa.binary(16)), pa.array([st] * len(buf), pa.string())],
                                    schema=schema,
                                )
                            )
                            buf.clear()
                if buf:
                    writer.write_batch(
                        pa.RecordBatch.from_arrays(
                            [pa.array(buf, pa.binary(16)), pa.array([st] * len(buf), pa.string())],
                            schema=schema,
                        )
                    )
    finally:
        writer.close()
        con.close()
    return owners_pq


def _write_fasta_record(fh, header: str, seq: bytes, width: int = 80) -> None:
    fh.write(f">{header}\n".encode("utf-8"))
    for i in range(0, len(seq), width):
        fh.write(seq[i : i + width] + b"\n")


def _finalize_staging(
    *,
    contigs_pq: Path,
    fasta_path: Path,
    final_fasta_path: Path,
    out_meta_dir: Path,
    asm_rows: list[dict],
    trim_rows: list[dict],
    tot: Counter,
    min_contig_len: int,
    root: Path | None,
    root_stages: list[str],
    owner_inputs: list[tuple[str, Path]],
    scratch: Path,
    retrieved_at: str,
) -> dict:
    """Deduplicate staged contigs in DuckDB, check earlier owners, and write genomes.json / objects.json / reports."""
    out_meta_dir.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(
        config={
            "threads": 2,
            "memory_limit": "4GB",
            "temp_directory": str(scratch / "duck"),
            "preserve_insertion_order": False,
        }
    )
    try:
        lit = _sql_lit(contigs_pq)
        con.execute(
            f"""
            CREATE TEMP TABLE lr_distinct AS
            SELECT
                contig_uid,
                max(contig_len)::BIGINT AS contig_len,
                bool_or(circular) AS is_circular,
                count(DISTINCT asm_idx)::INTEGER AS n_assemblies,
                count(*)::INTEGER AS n_records
            FROM read_parquet({lit})
            GROUP BY contig_uid
            ORDER BY contig_uid
            """
        )

        # Write genomes.json sorted by contig_uid (one genome per distinct contig)
        distinct_contigs = 0
        distinct_bp = 0
        circular_contigs = 0
        ge_100k = 0
        ge_100k_bp = 0
        multi_asm = 0

        reader = con.sql(
            "SELECT contig_uid, contig_len, is_circular, n_assemblies FROM lr_distinct ORDER BY contig_uid"
        ).to_arrow_reader(ROW_GROUP)
        with open(out_meta_dir / "genomes.json", "w") as gj:
            gj.write("[\n")
            first = True
            for b in reader:
                cu, ln, circ, na = (b.column(j).to_pylist() for j in range(4))
                lines = []
                for u, n, c, a in zip(cu, ln, circ, na):
                    uhex = bytes(u).hex()
                    distinct_contigs += 1
                    distinct_bp += n
                    if c:
                        circular_contigs += 1
                    if n >= 100_000:
                        ge_100k += 1
                        ge_100k_bp += n
                    if a > 1:
                        multi_asm += 1
                    rec = {
                        "contig_uids": [uhex],
                        "genetic_code": 11,
                        "kind": "prokaryote",
                        "circular": [uhex] if c else [],
                        "assembly_quality": "contig",
                    }
                    lines.append(json.dumps(rec))
                if lines:
                    gj.write(("" if first else ",\n") + ",\n".join(lines))
                    first = False
            gj.write("\n]\n")

        # Owner overlap against earlier stages
        held_by: dict[str, int] = {}
        held_bp: dict[str, int] = {}
        total_held = 0
        total_held_bp = 0
        owners_pq = _load_earlier_owners(root, root_stages, owner_inputs, scratch)
        with open(out_meta_dir / "shared_records.tsv", "w", newline="") as sfh:
            sw = csv.writer(sfh, delimiter="\t", lineterminator="\n")
            sw.writerow(["contig_uid", "owner_stage", "contig_len"])
            if owners_pq is not None and owners_pq.stat().st_size > 0:
                stage_order = list(root_stages) + [name for name, _ in owner_inputs]
                order_case = " ".join(
                    f"WHEN owner = {_sql_lit(st)} THEN {i}" for i, st in enumerate(stage_order)
                )
                own_lit = _sql_lit(owners_pq)
                q = f"""
                    WITH ranked AS (
                        SELECT
                            o.contig_uid,
                            o.owner,
                            d.contig_len,
                            row_number() OVER (
                                PARTITION BY o.contig_uid
                                ORDER BY CASE {order_case} ELSE 9999 END ASC
                            ) AS rn
                        FROM read_parquet({own_lit}) o
                        JOIN lr_distinct d ON o.contig_uid = d.contig_uid
                    )
                    SELECT contig_uid, owner, contig_len
                    FROM ranked
                    WHERE rn = 1
                    ORDER BY contig_uid
                """
                for u, own, ln in con.sql(q).fetchall():
                    sw.writerow([bytes(u).hex(), own, ln])
                    held_by[own] = held_by.get(own, 0) + 1
                    held_bp[own] = held_bp.get(own, 0) + ln
                    total_held += 1
                    total_held_bp += ln
    finally:
        con.close()

    # Write trims.tsv
    with open(out_meta_dir / "trims.tsv", "w", newline="") as tfh:
        tw = csv.DictWriter(
            tfh,
            fieldnames=["accession", "contig_header", "dtr_class", "orig_len", "trimmed_bases", "new_len", "contig_uid"],
            delimiter="\t",
            lineterminator="\n",
        )
        tw.writeheader()
        tw.writerows(trim_rows)

    # Write assemblies.tsv
    asm_cols = [
        "accession", "source", "habitat", "platform_class", "retrieved_at",
        "records_read", "bases_read", "records_below_min_len", "bases_below_min_len",
        "records_kept", "bases_kept", "circular_kept", "dtr_exact_trimmed", "dtr_trimmed_bases",
    ]
    with open(out_meta_dir / "assemblies.tsv", "w", newline="") as afh:
        aw = csv.DictWriter(afh, fieldnames=asm_cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        aw.writeheader()
        aw.writerows(asm_rows)

    # Write objects.json
    objects = [
        {
            "path": str(final_fasta_path),
            "source_name": "lr1_longread_metagenomes",
            "source_release": "ENA/NCBI public long-read metagenome assemblies (tranche LR1; DEC-57, DEC-65, DEC-71)",
            "license": INSDC_LICENSE,
            "redistributable": True,
            "retrieved_at": retrieved_at,
        }
    ]
    (out_meta_dir / "objects.json").write_text(json.dumps(objects, indent=2) + "\n")

    fasta_bytes = fasta_path.stat().st_size if fasta_path.is_file() else 0
    fasta_sha256 = _sha256_file(fasta_path) if fasta_path.is_file() else ""

    rep = {
        "stage": "LR1",
        "min_contig_len": min_contig_len,
        "assemblies_staged": len(asm_rows),
        "fasta_path": str(fasta_path),
        "fasta_bytes": fasta_bytes,
        "fasta_sha256": fasta_sha256,
        "records_read": tot["records_read"],
        "bases_read": tot["bases_read"],
        "records_below_min_len": tot["records_below_min_len"],
        "bases_below_min_len": tot["bases_below_min_len"],
        "records_kept": tot["records_kept"],
        "bases_kept": tot["bases_kept"],
        "distinct_contigs": distinct_contigs,
        "distinct_bp": distinct_bp,
        "circular_contigs": circular_contigs,
        "distinct_contigs_ge_100000": ge_100k,
        "distinct_bp_ge_100000": ge_100k_bp,
        "contigs_in_2plus_assemblies": multi_asm,
        "dtr_exact_trimmed": tot["dtr_exact_trimmed"],
        "dtr_trimmed_bases": tot["dtr_trimmed_bases"],
        "dtr_none_circular": tot["dtr_none_circular"],
        "dtr_demoted_linear": tot["dtr_demoted_linear"],
        "held_by_earlier": held_by,
        "held_bp_by_earlier": held_bp,
        "total_held_contigs": total_held,
        "total_held_bp": total_held_bp,
        "net_new_contigs": distinct_contigs - total_held,
        "net_new_bp": distinct_bp - total_held_bp,
    }
    (out_meta_dir / "input_report.json").write_text(json.dumps(rep, indent=2) + "\n")
    return rep


def _process_one_assembly(
    *,
    asm_idx: int,
    acc: str,
    fa_path: Path,
    meta: dict,
    min_contig_len: int,
    out_fasta_fh,
    pq_writer: pq.ParquetWriter,
    buf: list[list],
    trim_rows: list[dict],
    tot: Counter,
    seen_uids: set[bytes] | None = None,
) -> dict:
    """Stream one assembly FASTA, filter >= min_contig_len, apply DEC-65, write kept records."""
    s = Counter()
    for rec in fasta.iter_records(fa_path):
        raw_seq = rec.seq.upper()
        n = len(raw_seq)
        s["records_read"] += 1
        s["bases_read"] += n
        if n < min_contig_len:
            s["records_below_min_len"] += 1
            s["bases_below_min_len"] += n
            continue
        staged_seq, circ, dtr_cls, trimmed_k = apply_dec65_longread(rec.header, raw_seq)
        if len(staged_seq) < min_contig_len:
            s["records_below_min_len"] += 1
            s["bases_below_min_len"] += n
            continue
        uid = ids.contig_uid(staged_seq)
        if dtr_cls == "exact":
            s["dtr_exact_trimmed"] += 1
            s["dtr_trimmed_bases"] += trimmed_k
            trim_rows.append(
                {
                    "accession": acc,
                    "contig_header": rec.header.split()[0],
                    "dtr_class": dtr_cls,
                    "orig_len": n,
                    "trimmed_bases": trimmed_k,
                    "new_len": len(staged_seq),
                    "contig_uid": uid.hex(),
                }
            )
        elif dtr_cls == "none":
            s["dtr_none_circular"] += 1
        elif dtr_cls:
            s["dtr_demoted_linear"] += 1

        s["records_kept"] += 1
        s["bases_kept"] += len(staged_seq)
        if circ:
            s["circular_kept"] += 1

        # Avoid writing duplicate byte-identical contigs into the concatenated FASTA
        if seen_uids is None or uid not in seen_uids:
            if seen_uids is not None:
                seen_uids.add(uid)
            _write_fasta_record(out_fasta_fh, f"{acc}|{rec.header.split()[0]}", staged_seq)

        buf[0].append(uid)
        buf[1].append(len(staged_seq))
        buf[2].append(circ)
        buf[3].append(asm_idx)
        if len(buf[0]) >= ROW_GROUP:
            pq_writer.write_batch(
                pa.RecordBatch.from_arrays(
                    [
                        pa.array(buf[0], pa.binary(16)),
                        pa.array(buf[1], pa.int64()),
                        pa.array(buf[2], pa.bool_()),
                        pa.array(buf[3], pa.int32()),
                    ],
                    schema=_CONTIG_SCHEMA,
                )
            )
            for b in buf:
                b.clear()

    for k, v in s.items():
        tot[k] += v

    return {
        "accession": acc,
        "source": meta.get("source", ""),
        "habitat": meta.get("habitat", ""),
        "platform_class": meta.get("platform_class", ""),
        "retrieved_at": meta.get("retrieved_at", ""),
        **{k: s[k] for k in s},
    }


def stage_assemblies(
    manifest: list[dict],
    *,
    out_fasta_dir: Path,
    out_meta_dir: Path,
    final_fasta_dir: Path | None = None,
    min_contig_len: int = DEFAULT_MIN_CONTIG_LEN,
    root: Path | None = None,
    root_stages: list[str] = (),
    owner_inputs: list[tuple[str, Path]] = (),
) -> dict:
    """Stage local long-read metagenome assembly FASTAs into concatenated >= min_contig_len FASTA + metadata."""
    out_fasta_dir = Path(out_fasta_dir)
    out_meta_dir = Path(out_meta_dir)
    final_fasta_dir = Path(final_fasta_dir) if final_fasta_dir else out_fasta_dir
    out_fasta_dir.mkdir(parents=True, exist_ok=True)
    out_meta_dir.mkdir(parents=True, exist_ok=True)

    fasta_tmp = out_fasta_dir / f".{LR1_FASTA_NAME}.tmp"
    fasta_final = out_fasta_dir / LR1_FASTA_NAME
    scratch = Path(tempfile.mkdtemp(prefix="c56-stage-", dir=out_meta_dir))
    contigs_pq = scratch / "contigs.parquet"
    pq_writer = pq.ParquetWriter(contigs_pq, _CONTIG_SCHEMA, compression="zstd", write_statistics=False)
    buf: list[list] = [[], [], [], []]
    trim_rows: list[dict] = []
    asm_rows: list[dict] = []
    tot: Counter = Counter()
    seen_uids: set[bytes] = set()

    try:
        items = sorted(manifest, key=lambda r: r["accession"])
        with open(fasta_tmp, "wb") as out_fh:
            for idx, item in enumerate(items):
                row = _process_one_assembly(
                    asm_idx=idx,
                    acc=item["accession"],
                    fa_path=Path(item["path"]),
                    meta=item,
                    min_contig_len=min_contig_len,
                    out_fasta_fh=out_fh,
                    pq_writer=pq_writer,
                    buf=buf,
                    trim_rows=trim_rows,
                    tot=tot,
                    seen_uids=seen_uids,
                )
                asm_rows.append(row)
        if buf[0]:
            pq_writer.write_batch(
                pa.RecordBatch.from_arrays(
                    [
                        pa.array(buf[0], pa.binary(16)),
                        pa.array(buf[1], pa.int64()),
                        pa.array(buf[2], pa.bool_()),
                        pa.array(buf[3], pa.int32()),
                    ],
                    schema=_CONTIG_SCHEMA,
                )
            )
            for b in buf:
                b.clear()
        pq_writer.close()
        fasta_tmp.replace(fasta_final)
        ret_at = items[0].get("retrieved_at") or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        return _finalize_staging(
            contigs_pq=contigs_pq,
            fasta_path=fasta_final,
            final_fasta_path=final_fasta_dir / LR1_FASTA_NAME,
            out_meta_dir=out_meta_dir,
            asm_rows=asm_rows,
            trim_rows=trim_rows,
            tot=tot,
            min_contig_len=min_contig_len,
            root=root,
            root_stages=list(root_stages),
            owner_inputs=list(owner_inputs),
            scratch=scratch,
            retrieved_at=ret_at,
        )
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
        shutil.rmtree(out_fasta_dir / ".tmp_dl", ignore_errors=True)
        if fasta_tmp.exists():
            fasta_tmp.unlink()


# --- 4. Streaming fetch-and-stage (strict /tmp < 1 GB, /stor >= 64 MB) ----------------------------------------------

def _download_one(url: str, dst: Path, timeout: int = 300) -> tuple[int, str, int]:
    """Download URL with User-Agent: curl. Returns (http_status, error_msg, bytes_written)."""
    req = urllib.request.Request(url, headers={"User-Agent": "curl"})
    tmp = dst.with_name(dst.name + ".part")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            status = getattr(resp, "status", 200)
            n = 0
            with open(tmp, "wb") as out:
                while chunk := resp.read(1 << 20):
                    out.write(chunk)
                    n += len(chunk)
        tmp.replace(dst)
        return status, "", n
    except urllib.error.HTTPError as e:
        if tmp.exists():
            tmp.unlink()
        return int(e.code), str(e), 0
    except Exception as e:
        if tmp.exists():
            tmp.unlink()
        return 0, str(e), 0


def fetch_and_stage(
    *,
    selection_tsv: Path,
    out_fasta_dir: Path,
    out_meta_dir: Path,
    tmp_dl: Path,
    final_fasta_dir: Path | None = None,
    target_bases: int = DEFAULT_TARGET_BASES,
    max_contigs: int = DEFAULT_MAX_CONTIGS,
    min_contig_len: int = DEFAULT_MIN_CONTIG_LEN,
    interval: float = 0.6,
    root: Path | None = None,
    root_stages: list[str] = (),
    owner_inputs: list[tuple[str, Path]] = (),
) -> dict:
    """Stream-fetch assemblies one at a time, filter >= min_contig_len + DEC-65, unlink raw .fa.gz immediately."""
    out_fasta_dir = Path(out_fasta_dir)
    out_meta_dir = Path(out_meta_dir)
    tmp_dl = Path(tmp_dl)
    final_fasta_dir = Path(final_fasta_dir) if final_fasta_dir else out_fasta_dir
    out_fasta_dir.mkdir(parents=True, exist_ok=True)
    out_meta_dir.mkdir(parents=True, exist_ok=True)
    tmp_dl.mkdir(parents=True, exist_ok=True)

    with open(selection_tsv, newline="") as fh:
        sel_rows = list(csv.DictReader(fh, delimiter="\t"))

    fasta_tmp = out_fasta_dir / f".{LR1_FASTA_NAME}.tmp"
    fasta_final = out_fasta_dir / LR1_FASTA_NAME
    scratch = Path(tempfile.mkdtemp(prefix="c56-stage-", dir=out_meta_dir))
    contigs_pq = scratch / "contigs.parquet"
    pq_writer = pq.ParquetWriter(contigs_pq, _CONTIG_SCHEMA, compression="zstd", write_statistics=False)
    buf: list[list] = [[], [], [], []]
    trim_rows: list[dict] = []
    asm_rows: list[dict] = []
    tot: Counter = Counter()
    seen_uids: set[bytes] = set()
    log_path = out_meta_dir / "fetch_log.jsonl"

    consec_refusals = 0
    last_req = 0.0
    fetched_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    try:
        with open(fasta_tmp, "wb") as out_fh, open(log_path, "w") as lfh:
            for idx, r in enumerate(sel_rows):
                if tot["bases_kept"] >= target_bases or tot["records_kept"] >= max_contigs:
                    break
                acc = r["accession"]
                url = r["ftp_url"]
                expected_md5 = (r.get("md5") or "").strip()
                if not url:
                    continue

                wait = interval - (time.monotonic() - last_req)
                if wait > 0:
                    time.sleep(wait)
                last_req = time.monotonic()

                suffix = ".fa.gz" if url.endswith(".gz") else ".fa"
                est_gz = int(r.get("gz_bytes") or 0)
                if est_gz >= 500_000_000:
                    stor_tmp = out_fasta_dir / ".tmp_dl"
                    stor_tmp.mkdir(parents=True, exist_ok=True)
                    dl_file = stor_tmp / f"{acc}{suffix}"
                else:
                    dl_file = tmp_dl / f"{acc}{suffix}"
                t0 = time.monotonic()
                status, err, nbytes = _download_one(url, dl_file)
                dt = round(time.monotonic() - t0, 2)
                ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

                if status in (403, 429, 500, 502, 503, 504) or status == 0:
                    consec_refusals += 1
                    lfh.write(
                        json.dumps(
                            {"accession": acc, "url": url, "status": status, "error": err, "utc": ts, "sec": dt}
                        )
                        + "\n"
                    )
                    lfh.flush()
                    if consec_refusals >= 5:
                        raise RuntimeError(f"Stopped after 5 consecutive fetch refusals/errors (last: {acc} HTTP {status}: {err})")
                    continue
                if status != 200:
                    consec_refusals = 0
                    lfh.write(
                        json.dumps(
                            {"accession": acc, "url": url, "status": status, "error": err, "utc": ts, "sec": dt}
                        )
                        + "\n"
                    )
                    lfh.flush()
                    continue

                consec_refusals = 0
                got_md5 = _md5_file(dl_file)
                got_sha = _sha256_file(dl_file)
                if expected_md5 and got_md5 != expected_md5:
                    dl_file.unlink(missing_ok=True)
                    lfh.write(
                        json.dumps(
                            {
                                "accession": acc, "url": url, "status": "md5_mismatch",
                                "expected_md5": expected_md5, "got_md5": got_md5, "utc": ts,
                            }
                        )
                        + "\n"
                    )
                    lfh.flush()
                    continue

                try:
                    row = _process_one_assembly(
                        asm_idx=idx,
                        acc=acc,
                        fa_path=dl_file,
                        meta={**r, "retrieved_at": ts},
                        min_contig_len=min_contig_len,
                        out_fasta_fh=out_fh,
                        pq_writer=pq_writer,
                        buf=buf,
                        trim_rows=trim_rows,
                        tot=tot,
                        seen_uids=seen_uids,
                    )
                    asm_rows.append(row)
                    lfh.write(
                        json.dumps(
                            {
                                "accession": acc, "url": url, "status": "ok", "wire_bytes": nbytes,
                                "md5": got_md5, "sha256": got_sha, "utc": ts, "sec": dt,
                                "records_kept": row["records_kept"], "bases_kept": row["bases_kept"],
                            }
                        )
                        + "\n"
                    )
                    lfh.flush()
                except Exception as exc:
                    lfh.write(
                        json.dumps(
                            {"accession": acc, "url": url, "status": "parse_error", "error": str(exc), "utc": ts}
                        )
                        + "\n"
                    )
                    lfh.flush()
                finally:
                    dl_file.unlink(missing_ok=True)

        if buf[0]:
            pq_writer.write_batch(
                pa.RecordBatch.from_arrays(
                    [
                        pa.array(buf[0], pa.binary(16)),
                        pa.array(buf[1], pa.int64()),
                        pa.array(buf[2], pa.bool_()),
                        pa.array(buf[3], pa.int32()),
                    ],
                    schema=_CONTIG_SCHEMA,
                )
            )
            for b in buf:
                b.clear()
        pq_writer.close()
        fasta_tmp.replace(fasta_final)

        return _finalize_staging(
            contigs_pq=contigs_pq,
            fasta_path=fasta_final,
            final_fasta_path=final_fasta_dir / LR1_FASTA_NAME,
            out_meta_dir=out_meta_dir,
            asm_rows=asm_rows,
            trim_rows=trim_rows,
            tot=tot,
            min_contig_len=min_contig_len,
            root=root,
            root_stages=list(root_stages),
            owner_inputs=list(owner_inputs),
            scratch=scratch,
            retrieved_at=fetched_at,
        )
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
        shutil.rmtree(out_fasta_dir / ".tmp_dl", ignore_errors=True)
        if fasta_tmp.exists():
            fasta_tmp.unlink()


def _parse_owner_inputs(items: list[str] | None) -> list[tuple[str, Path]]:
    out = []
    for item in items or ():
        name, path = item.split("=", 1)
        out.append((name, Path(path)))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    cw = sub.add_parser("crosswalk")
    cw.add_argument("--read-run-lr", required=True)
    cw.add_argument("--ena-analysis", required=True)
    cw.add_argument("--ena-wgs", required=True)
    cw.add_argument("--ncbi-summary", required=True)
    cw.add_argument("--img-m", default="")
    cw.add_argument("--out-dir", required=True)
    cw.add_argument("--max-bases", type=int, default=DEFAULT_TARGET_BASES)
    cw.add_argument("--max-single-bp", type=int, default=0)
    cw.add_argument("--seed", type=int, default=DEFAULT_SEED)

    fs = sub.add_parser("fetch-and-stage")
    fs.add_argument("--selection", required=True)
    fs.add_argument("--out-fasta-dir", required=True)
    fs.add_argument("--out-meta-dir", required=True)
    fs.add_argument("--tmp-dl", required=True)
    fs.add_argument("--final-fasta-dir", default="")
    fs.add_argument("--target-bases", type=int, default=DEFAULT_TARGET_BASES)
    fs.add_argument("--max-contigs", type=int, default=DEFAULT_MAX_CONTIGS)
    fs.add_argument("--min-contig-len", type=int, default=DEFAULT_MIN_CONTIG_LEN)
    fs.add_argument("--interval", type=float, default=0.6)
    fs.add_argument("--root", default="")
    fs.add_argument("--root-stages", default="")
    fs.add_argument("--owner-inputs", action="append", default=[])

    st = sub.add_parser("stage")
    st.add_argument("--manifest", required=True)
    st.add_argument("--out-fasta-dir", required=True)
    st.add_argument("--out-meta-dir", required=True)
    st.add_argument("--final-fasta-dir", default="")
    st.add_argument("--min-contig-len", type=int, default=DEFAULT_MIN_CONTIG_LEN)
    st.add_argument("--root", default="")
    st.add_argument("--root-stages", default="")
    st.add_argument("--owner-inputs", action="append", default=[])

    args = ap.parse_args(argv)
    if args.cmd == "crosswalk":
        summary = build_crosswalk(
            read_run_lr=Path(args.read_run_lr),
            ena_analysis=Path(args.ena_analysis),
            ena_wgs=Path(args.ena_wgs),
            ncbi_summary=Path(args.ncbi_summary),
            img_m=Path(args.img_m) if args.img_m else None,
            out_dir=Path(args.out_dir),
            max_bases=args.max_bases,
            max_single_bp=args.max_single_bp or None,
            seed=args.seed,
        )
        print(json.dumps(summary, indent=2))
        return 0

    root_stages = [s.strip() for s in args.root_stages.split(",") if s.strip()]
    owner_inputs = _parse_owner_inputs(args.owner_inputs)
    if args.cmd == "fetch-and-stage":
        rep = fetch_and_stage(
            selection_tsv=Path(args.selection),
            out_fasta_dir=Path(args.out_fasta_dir),
            out_meta_dir=Path(args.out_meta_dir),
            tmp_dl=Path(args.tmp_dl),
            final_fasta_dir=Path(args.final_fasta_dir) if args.final_fasta_dir else None,
            target_bases=args.target_bases,
            max_contigs=args.max_contigs,
            min_contig_len=args.min_contig_len,
            interval=args.interval,
            root=Path(args.root) if args.root else None,
            root_stages=root_stages,
            owner_inputs=owner_inputs,
        )
        print(json.dumps(rep, indent=2))
        return 0

    manifest = json.loads(Path(args.manifest).read_text())
    rep = stage_assemblies(
        manifest,
        out_fasta_dir=Path(args.out_fasta_dir),
        out_meta_dir=Path(args.out_meta_dir),
        final_fasta_dir=Path(args.final_fasta_dir) if args.final_fasta_dir else None,
        min_contig_len=args.min_contig_len,
        root=Path(args.root) if args.root else None,
        root_stages=root_stages,
        owner_inputs=owner_inputs,
    )
    print(json.dumps(rep, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
