"""Card C-86 (DEC-85 (4), DEC-80, DEC-88): stage light tower-landed r26.11.6 cohorts on wilkcomp02:
- `O4`: protist weird mitochondria, plastids/apicoplasts/chromatophores, cryptophyte nucleomorphs (`O4_ORG`)
        + fungal mitochondrial genomes (`O3_FUNGAL`)
- `GV2`: polar, freshwater, alpine-lake, and WWTP giant viruses (`GV_POLAR_2026`, `GV_POLAR_PITOT`,
         `GV_POLAR_GEVE`, `GV_GKS`, `GV_POLAR_MARINE`, `GV_WWTP`)
- `HROM1`: Human Reference Oral Microbiome species reps incl. novel CPR (`HROM_ORAL`) + alternative-code
           leafhopper symbiont/gut MAGs (`CODE4_GLMC`) + code-4 Actinomycetota (`CODE4_ACTINO1`)

Subcommands:
    python c86_o4_gv2_hrom1_inputs.py census PATH [PATH ...] [--out JSON]
    python c86_o4_gv2_hrom1_inputs.py o4-stage --o4-org-dir DIR --fungal-fa FASTA --fungal-tsv TSV
        --meta-out DIR --fasta-out DIR [common options]
    python c86_o4_gv2_hrom1_inputs.py gv2-stage --dir SUB=DIR ... --meta-out DIR --fasta-out DIR [common options]
    python c86_o4_gv2_hrom1_inputs.py hrom1-stage --hrom-tar TAR --hrom-meta TSV --glmc-zip ZIP --glmc-csv CSV
        --actino-dir DIR --index-dir DIR --meta-out DIR --fasta-out DIR [common options]
    python c86_o4_gv2_hrom1_inputs.py hrom1-ani META_DIR OUT_DIR --fasta FASTA ... [--whole NAME=META:FASTA[,FASTA]]
        [--by-genus NAME=META:FASTA[,FASTA]] [--taxonomy TSV] [--archaea FASTA:LIST] [--threads 8]
"""

from __future__ import annotations

import argparse
import bz2
import csv
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import warnings
import zipfile
from collections import Counter
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import a26_r26_11_4_chain as a26  # noqa: E402
import c39_ani  # noqa: E402
import c48c_virus_inputs as v2  # noqa: E402
import c58_sags_symbionts_niche as c58  # noqa: E402
import c59_refseq_organelles as c59  # noqa: E402
import c62_metavr_inputs as c62  # noqa: E402
import c77_sag2_inputs as c77  # noqa: E402
import c78_compact_tower_sets as c78  # noqa: E402
import c80_kavagutti_inphared as c80  # noqa: E402
import c81_gwv1_sag3_wet1 as c81  # noqa: E402
import make_inputs  # noqa: E402
import r232_index  # noqa: E402
from biobox import fasta, ids  # noqa: E402
from biobox.l1 import genecall  # noqa: E402
from biobox.search import orfs  # noqa: E402

csv.field_size_limit(1 << 30)
Refused = make_inputs.Refused

ROOT_STAGES: tuple[str, ...] = ("S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1")
R26_11_4_STAGES: tuple[str, ...] = tuple(s for s, _ in a26.R26_11_4_CHAIN)
R26_11_5_STAGES: tuple[str, ...] = (
    "GRARE1", "GRARE2", "GRARE3",
    "LR2", "SAG1", "NX1",
    "GV1", "EGV1", "ECE2", "MIR1",
    "SAG2",
    "DPN1", "SYM2", "LRP1", "MIR2",
    "MV2", "MVEUKV",
    "KAVHIFI", "INPH2604",
    "GWV1", "SAG3", "WET1",
    "RMD1", "GBRMGD",
)
DEFAULT_R26_11_5_DIRS: dict[str, str] = {
    "GRARE1": "/tmp/if2453/c75/meta/GRARE1",
    "GRARE2": "/tmp/if2453/c75/meta/GRARE2",
    "GRARE3": "/tmp/if2453/c75/meta/GRARE3",
    "LR2": "/tmp/if2453/r19/LR2",
    "SAG1": "/tmp/if2453/c58/SAG1",
    "NX1": "/tmp/if2453/r20/NX1",
    "GV1": "/tmp/if2453/r23/GV1",
    "EGV1": "/tmp/if2453/c69/EGV1",
    "ECE2": "/tmp/if2453/c70/meta/ECE2",
    "MIR1": "/tmp/if2453/c70/meta/MIR1",
    "SAG2": "/tmp/if2453/c77/meta/SAG2",
    "DPN1": "/tmp/if2453/c78/meta/DPN1",
    "SYM2": "/tmp/if2453/c78/meta/SYM2",
    "LRP1": "/tmp/if2453/c78/meta/LRP1",
    "MIR2": "/tmp/if2453/c78/meta/MIR2",
    "MV2": "/tmp/if2453/c79/MV2",
    "MVEUKV": "/tmp/if2453/c79/MVEUKV",
    "KAVHIFI": "/tmp/if2453/c80/KAVHIFI",
    "INPH2604": "/tmp/if2453/c80/INPH2604",
    "GWV1": "/tmp/if2453/c81/meta/GWV1",
    "SAG3": "/tmp/if2453/c81/meta/SAG3",
    "WET1": "/tmp/if2453/c81/meta/WET1",
    "RMD1": "/tmp/if2453/c83/meta/RMD1",
    "GBRMGD": "/tmp/if2453/c83/meta/GBRMGD",
}

STOR_MIN_BYTES = 64 * 1024 * 1024  # AGENTS.md §7.3: >= 67,108,864 B on /stor
VALID_DNA_RE = re.compile(rb"^[ACGTURYKMSWBDHVNacgturykmswbdhvn]+$")
GAP_RE = re.compile(rb"[-.]")
PROTEIN_ONLY = set(b"EFILPQZJO*efilpqzjo")
CIRCULAR_RE = re.compile(r"circular", re.IGNORECASE)
ACC_GCA_RE = re.compile(r"GC[AF]_\d{9}(?:\.\d+)?")
FASTA_EXTS = (".fa", ".fna", ".fasta", ".fas", ".ffn", ".faa", ".fsa", ".fa.gz", ".fna.gz", ".fasta.gz", ".faa.gz")
DNA_FASTA_EXTS = (".fa", ".fna", ".fasta", ".fas", ".fsa", ".fa.gz", ".fna.gz", ".fasta.gz")

_CODE25_HROM_RE = re.compile(
    r"(patescibacteria|patescibacteriota|saccharibacteria|gracilibacteria|absconditabacteria|"
    r"absconditabacterales|bd1-5|sr1|jaedam01|gn02|tm7|candidatus\s+nanosynbacter)",
    re.IGNORECASE,
)
_CODE4_HROM_RE = re.compile(
    r"(hodgkinia|nasuia|zinderia|mycoplasma|ureaplasma|spiroplasma|mesoplasma|entomoplasma|"
    r"mycoplasmoides|mycoplasmopsis|metamycoplasma|malacoplasma|acholeplasma|phytoplasma|"
    r"equivita|tapirivita|gorillivita|altericodex)",
    re.IGNORECASE,
)


# ==================================================================================================================
# 1. Magic-byte format detection & archive member census
# ==================================================================================================================


def _magic_of_bytes(buf: bytes) -> str:
    if buf.startswith(b"\x1f\x8b"):
        return "gzip"
    if buf.startswith(b"BZh"):
        return "bzip2"
    if buf.startswith(b"\xfd7zXZ\x00"):
        return "xz"
    if buf[:4] in (b"PK\x03\x04", b"PK\x05\x06", b"PK\x07\x08"):
        return "zip"
    if len(buf) >= 262 and buf[257:262] == b"ustar":
        return "tar"
    s = buf.lstrip()
    if s.startswith(b"LOCUS ") or s.startswith(b"LOCUS\t"):
        return "gbff"
    if s.startswith(b">"):
        return "fasta"
    try:
        buf.decode("utf-8")
        return "text"
    except UnicodeDecodeError:
        return "binary"


def magic_format(path: Path | str) -> str:
    with open(path, "rb") as fh:
        return _magic_of_bytes(fh.read(512))


def _looks_fasta(name: str) -> bool:
    base = os.path.basename(name.lower())
    return not base.startswith(".") and base.endswith(FASTA_EXTS)


def _looks_dna_fasta(name: str) -> bool:
    base = os.path.basename(name.lower())
    return not base.startswith(".") and not base.endswith((".proteins.fa", ".proteins.fasta")) and base.endswith(DNA_FASTA_EXTS)


def _seq_alphabet(seq: bytes) -> str:
    if not seq:
        return "empty"
    if GAP_RE.search(seq):
        return "alignment"
    if VALID_DNA_RE.match(seq):
        return "dna"
    if PROTEIN_ONLY & set(seq):
        return "protein"
    return "non_iupac"


def _fasta_member_stats(fobj, name: str, size: int, max_headers: int = 3) -> dict:
    n = bp = circ = non_iupac_recs = 0
    mn, mx = None, 0
    heads: list[str] = []
    kinds: set[str] = set()
    try:
        for r in v2.iter_fasta_stream(fobj):
            n += 1
            bp += len(r.seq)
            mn = len(r.seq) if mn is None else min(mn, len(r.seq))
            mx = max(mx, len(r.seq))
            circ += bool(CIRCULAR_RE.search(r.header))
            if len(heads) < max_headers:
                heads.append(r.header[:160])
            k = _seq_alphabet(r.seq)
            kinds.add(k)
            if k != "dna":
                non_iupac_recs += 1
    except (ValueError, UnicodeDecodeError, OSError, EOFError):
        return {"name": name, "bytes": size, "alphabet": "not_fasta"}
    if n == 0:
        return {"name": name, "bytes": size, "alphabet": "not_fasta"}
    if kinds == {"dna"}:
        alpha = "dna"
    else:
        alpha = next(k for k in ("non_iupac", "alignment", "protein", "empty", "dna") if k in kinds)
    return {
        "name": name,
        "bytes": size,
        "alphabet": alpha,
        "records": n,
        "bases": bp,
        "min_len": mn,
        "max_len": mx,
        "circular_in_header": circ,
        "non_dna_records": non_iupac_recs,
        "first_headers": heads,
    }


def _gbff_stats(path: Path, magic: str) -> dict:
    fh = gzip.open(path, "rt", encoding="utf-8", errors="replace") if magic == "gzip" else open(
        path, "r", encoding="utf-8", errors="replace"
    )
    n_rec = bp = n_circ = 0
    moltypes: Counter[str] = Counter()
    comps: Counter[str] = Counter()
    tables: Counter[str] = Counter()
    spliced_recs = 0
    try:
        for rec in c59.parse_refseq_organelle_gbff(fh, keep_translations=False):
            n_rec += 1
            bp += rec["length"]
            n_circ += bool(rec["circular"])
            moltypes[rec["moltype"] or "unknown"] += 1
            comp = normalize_o4_compartment(rec["organelle"], rec["definition"])
            comps[comp] += 1
            tables[str(rec["transl_table"])] += 1
            if rec["n_spliced_cds"] > 0:
                spliced_recs += 1
    finally:
        fh.close()
    return {
        "magic": magic,
        "bytes": path.stat().st_size,
        "gbff_records": n_rec,
        "gbff_bases": bp,
        "circular": n_circ,
        "moltypes": dict(sorted(moltypes.items())),
        "compartments": dict(sorted(comps.items())),
        "transl_tables": dict(sorted(tables.items())),
        "records_with_spliced_cds": spliced_recs,
    }


def census(paths: Iterable[Path | str]) -> dict:
    out: dict[str, dict] = {}
    for p in map(Path, paths):
        mg = magic_format(p)
        low = p.name.lower()
        if mg == "gbff" or low.endswith((".gbff", ".gbff.gz", ".gbk", ".gbk.gz")):
            out[str(p)] = _gbff_stats(p, mg)
            continue
        members: list[dict] = []
        if mg == "zip":
            with zipfile.ZipFile(p) as zf:
                for info in zf.infolist():
                    if info.is_dir():
                        members.append({"name": info.filename, "bytes": 0, "alphabet": "directory"})
                    elif _looks_fasta(info.filename) and not info.filename.startswith("__MACOSX/"):
                        with zf.open(info) as raw:
                            stream = gzip.open(raw, "rb") if info.filename.lower().endswith(".gz") else raw
                            members.append(_fasta_member_stats(stream, info.filename, info.file_size))
                    else:
                        members.append({"name": info.filename, "bytes": info.file_size, "alphabet": "not_fasta"})
        elif mg == "tar" or (mg in ("gzip", "bzip2", "xz") and (".tar." in low or low.endswith(".tgz"))):
            with tarfile.open(p, "r:*") as tf:
                for m in tf:
                    if m.isdir():
                        members.append({"name": m.name.rstrip("/") + "/", "bytes": 0, "alphabet": "directory"})
                    elif m.isfile() and _looks_fasta(m.name):
                        ef = tf.extractfile(m)
                        assert ef is not None
                        stream = gzip.open(ef, "rb") if m.name.lower().endswith(".gz") else ef
                        members.append(_fasta_member_stats(stream, m.name, m.size))
                    else:
                        members.append({"name": m.name, "bytes": m.size, "alphabet": "not_fasta"})
        else:
            with open(p, "rb") as raw:
                if mg == "gzip":
                    stream = gzip.open(raw, "rb")
                elif mg == "bzip2":
                    stream = bz2.open(raw, "rb")
                else:
                    stream = raw
                members.append(_fasta_member_stats(stream, p.name, p.stat().st_size))
        kinds = Counter(m["alphabet"] for m in members)
        out[str(p)] = {
            "magic": mg,
            "bytes": p.stat().st_size,
            "members": members,
            "by_alphabet": dict(sorted(kinds.items())),
            "dna_members": sum(1 for m in members if m["alphabet"] == "dna"),
            "dna_records": sum(m.get("records", 0) for m in members if m["alphabet"] == "dna"),
            "dna_bases": sum(m.get("bases", 0) for m in members if m["alphabet"] == "dna"),
        }
    return out


# ==================================================================================================================
# 2. Genetic code screens (1, 3, 4, 11, 16, 22, 25)
# ==================================================================================================================


def check_supported_code(code: int) -> None:
    try:
        orfs.check_code(int(code))
    except ValueError as exc:
        raise Refused(str(exc)) from exc


def _train_call_stats(seq_strs: list[str], total_bp: int, table: int) -> dict:
    """Train single-mode pyrodigal on `table` (repeating short genomes >=20 kb solely for training) and score ORFs."""
    import pyrodigal

    if total_bp <= 0:
        return {"cd": 0.0, "mean_orf": 0.0, "n_genes": 0, "n_cod": 0, "n_tga": 0, "n_tgg": 0, "n_tag": 0}
    finder = pyrodigal.GeneFinder(meta=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        if total_bp < 20_000:
            joined = "".join(seq_strs)
            reps = (20_000 // max(1, len(joined))) + 2
            finder.train(joined * reps, translation_table=table)
        else:
            finder.train(*seq_strs, translation_table=table)

    coding_bp = 0
    n_genes = 0
    n_cod = n_tga = n_tgg = n_tag = 0
    for s in seq_strs:
        for g in finder.find_genes(s):
            gl = g.end - g.begin + 1
            coding_bp += gl
            n_genes += 1
            nuc = g.sequence()
            start_off = 0 if g.partial_begin else 3
            end_off = len(nuc) if g.partial_end else max(start_off, len(nuc) - 3)
            for pos in range(start_off, end_off - 2, 3):
                cod = nuc[pos : pos + 3]
                n_cod += 1
                if cod == "TGA":
                    n_tga += 1
                elif cod == "TGG":
                    n_tgg += 1
                elif cod == "TAG":
                    n_tag += 1
    cd = coding_bp / total_bp
    mean_orf = (coding_bp / n_genes) if n_genes > 0 else 0.0
    return {
        "cd": cd,
        "mean_orf": mean_orf,
        "n_genes": n_genes,
        "n_cod": n_cod,
        "n_tga": n_tga,
        "n_tgg": n_tgg,
        "n_tag": n_tag,
    }


def screen_prokaryote_code(
    seqs: list[bytes],
    taxonomy_hint: str = "",
    *,
    force_alt: bool = False,
    screen_max_bp: int = 500_000,
) -> dict:
    """Screen prokaryotic MAGs/isolates across codes 11, 4 (UGA->Trp), and 25 (UGA->Gly)."""
    pre_strs, pre_bp, is_full = c58._prefix_seqs_for_screen(seqs, screen_max_bp)
    if pre_bp < 5_000:
        code = 4 if _CODE4_HROM_RE.search(taxonomy_hint) else (25 if _CODE25_HROM_RE.search(taxonomy_hint) else 11)
        return {
            "genetic_code": code,
            "reason": "lt5kb_tax_or_default",
            "tested_alt": False,
            "cd_11": None, "mean_orf_11": None, "n_genes_11": None,
            "cd_25": None, "mean_orf_25": None, "n_genes_25": None,
            "cd_4": None, "mean_orf_4": None, "n_genes_4": None,
            "tga_per_kb": None,
        }

    s11 = _train_call_stats(pre_strs, pre_bp, 11)
    tax_c25 = bool(_CODE25_HROM_RE.search(taxonomy_hint))
    tax_c4 = bool(_CODE4_HROM_RE.search(taxonomy_hint))
    triggered = force_alt or tax_c25 or tax_c4 or (s11["cd"] < 0.85) or (s11["mean_orf"] < 700.0)
    if not triggered:
        return {
            "genetic_code": 11,
            "reason": "code11_pass",
            "tested_alt": False,
            "cd_11": round(s11["cd"], 5),
            "mean_orf_11": round(s11["mean_orf"], 2),
            "n_genes_11": s11["n_genes"],
            "cd_25": None, "mean_orf_25": None, "n_genes_25": None,
            "cd_4": None, "mean_orf_4": None, "n_genes_4": None,
            "tga_per_kb": None,
        }

    if not is_full:
        full_strs = [s.upper().decode("ascii") for s in seqs]
        full_bp = sum(len(s) for s in seqs)
        s11 = _train_call_stats(full_strs, full_bp, 11)
    else:
        full_strs, full_bp = pre_strs, pre_bp

    s4 = _train_call_stats(full_strs, full_bp, 4)
    s25 = _train_call_stats(full_strs, full_bp, 25)
    best_cd = max(s4["cd"], s25["cd"])
    best_orf = max(s4["mean_orf"], s25["mean_orf"])
    cd_gain = best_cd - s11["cd"]
    orf_ratio = (best_orf / s11["mean_orf"]) if s11["mean_orf"] > 0 else float("inf")

    tga_per_kb: float | None = None
    if cd_gain >= 0.08 and orf_ratio >= 1.20 and s4["n_tga"] >= 3:
        coding_bp_4 = s4["cd"] * full_bp
        tga_per_kb = (1000.0 * s4["n_tga"] / coding_bp_4) if coding_bp_4 > 0 else 0.0
        llr_25 = (
            genecall._code25_llr(s4["n_cod"], s4["n_tga"], s4["n_tgg"])
            if s4["n_cod"] >= genecall.CODE25_MIN_CODONS
            else -100.0
        )
        if tax_c25 and not tax_c4:
            code = 25 if llr_25 >= 0.0 or s25["cd"] >= s4["cd"] - 0.01 else 4
        elif tax_c4 and not tax_c25:
            code = 4
        elif llr_25 >= genecall.CODE25_MIN_LLR:
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
        "cd_11": round(s11["cd"], 5),
        "mean_orf_11": round(s11["mean_orf"], 2),
        "n_genes_11": s11["n_genes"],
        "cd_25": round(s25["cd"], 5),
        "mean_orf_25": round(s25["mean_orf"], 2),
        "n_genes_25": s25["n_genes"],
        "cd_4": round(s4["cd"], 5),
        "mean_orf_4": round(s4["mean_orf"], 2),
        "n_genes_4": s4["n_genes"],
        "tga_per_kb": None if tga_per_kb is None else round(tga_per_kb, 3),
    }


def _prok_screen_worker(args: tuple[str, list[bytes], str, bool]) -> tuple[str, dict]:
    gid, seqs, tax, force_alt = args
    return gid, screen_prokaryote_code(seqs, tax, force_alt=force_alt)


def screen_fungal_mt_code(seq: bytes, taxonomy: str = "") -> dict:
    """Screen an unannotated fungal/protist mitochondrial sequence across tables 1, 3, 4, and 16."""
    s = seq.upper().decode("ascii")
    n = len(s)
    s1 = _train_call_stats([s], n, 1)
    s4 = _train_call_stats([s], n, 4)
    s16 = _train_call_stats([s], n, 16)
    tax_l = taxonomy.lower()

    best_1_4_cd = max(s1["cd"], s4["cd"])
    best_1_4_orf = max(s1["mean_orf"], s4["mean_orf"], 1.0)
    if (s16["cd"] - best_1_4_cd >= 0.05 or (s16["cd"] > best_1_4_cd + 0.02 and s16["mean_orf"] / best_1_4_orf >= 1.25)) and s16["n_tag"] >= 3:
        return {
            "genetic_code": 16,
            "reason": "tag_recoded_code16",
            "cd_1": round(s1["cd"], 5),
            "cd_4": round(s4["cd"], 5),
            "cd_16": round(s16["cd"], 5),
        }

    orf_ratio_4 = s4["mean_orf"] / max(1.0, s1["mean_orf"])
    if s4["n_tga"] >= 2 and ((s4["cd"] - s1["cd"] >= 0.02) or (orf_ratio_4 >= 1.15 and s4["cd"] >= s1["cd"] - 0.01)):
        if "saccharomycetaceae" in tax_l or "saccharomyces" in tax_l or "candida_glabrata" in tax_l or "nakaseomyces" in tax_l:
            s3 = _train_call_stats([s], n, 3)
            if s3["cd"] >= s4["cd"] - 0.02:
                return {
                    "genetic_code": 3,
                    "reason": "yeast_mt_code3",
                    "cd_1": round(s1["cd"], 5),
                    "cd_4": round(s4["cd"], 5),
                    "cd_16": round(s16["cd"], 5),
                }
        return {
            "genetic_code": 4,
            "reason": "tga_recoded_code4",
            "cd_1": round(s1["cd"], 5),
            "cd_4": round(s4["cd"], 5),
            "cd_16": round(s16["cd"], 5),
        }

    return {
        "genetic_code": 1,
        "reason": "standard_code1",
        "cd_1": round(s1["cd"], 5),
        "cd_4": round(s4["cd"], 5),
        "cd_16": round(s16["cd"], 5),
    }


def _fungal_screen_worker(args: tuple[str, bytes, str]) -> tuple[str, dict]:
    gid, seq, tax = args
    return gid, screen_fungal_mt_code(seq, tax)


def check_gv_code25(seqs: list[bytes]) -> dict:
    """Apply genecall.py's exact pyrodigal-gv code-25 LLR test across a viral genome's contigs."""
    import pyrodigal_gv

    finder = pyrodigal_gv.ViralGeneFinder(meta=True)
    total_bp = code4_bp = n_cod = n_tgg = n_tga = 0
    for raw in seqs:
        s = raw.upper().decode("ascii")
        for g in finder.find_genes(s):
            span = g.end - (g.begin - 1)
            total_bp += span
            if g.translation_table == 4:
                code4_bp += span
                aa4 = g.translate(include_stop=False)
                aa25 = g.translate(translation_table=25, include_stop=False)
                body = 0 if g.partial_begin else 1
                n_cod += max(0, len(aa4) - body)
                w4 = aa4.count("W", body)
                w25 = aa25.count("W", body)
                n_tgg += w25
                n_tga += max(0, w4 - w25)
    llr = genecall._code25_llr(n_cod, n_tga, n_tgg) if (n_cod >= genecall.CODE25_MIN_CODONS and n_tga > 0) else None
    use_c25 = bool(
        total_bp > 0
        and code4_bp * 2 >= total_bp
        and n_cod >= genecall.CODE25_MIN_CODONS
        and n_tga >= genecall.CODE25_MIN_TGA
        and llr is not None
        and llr >= genecall.CODE25_MIN_LLR
    )
    return {
        "use_c25": use_c25,
        "total_coding_bp": total_bp,
        "code4_bp": code4_bp,
        "n_cod": n_cod,
        "n_tga": n_tga,
        "n_tgg": n_tgg,
        "llr": None if llr is None else round(llr, 3),
    }


# ==================================================================================================================
# 3. Shared helpers: ownership, 64 MiB StageFiles, and object metadata
# ==================================================================================================================


def _tally_owners(contig_len: dict[bytes, int], holders: dict[bytes, str], skip: dict[bytes, str] | None = None):
    c: Counter[str] = Counter()
    b: Counter[str] = Counter()
    for cu, n in contig_len.items():
        if cu in holders and not (skip and cu in skip):
            c[holders[cu]] += 1
            b[holders[cu]] += n
    return dict(sorted(c.items())), dict(sorted(b.items()))


def _owners_of_uids(
    uids: set[bytes],
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    later_inputs: dict[str, Path],
    work: Path,
    threads: int,
    memory_limit: str,
) -> tuple[dict[bytes, str], dict[bytes, str]]:
    for name, d in later_inputs.items():
        if not (Path(d) / "genomes.json").is_file():
            raise Refused(f"--later {name}: no genomes.json under {d}")
    for name, d in owner_inputs.items():
        if not (Path(d) / "genomes.json").is_file():
            raise Refused(f"--owner {name}: no genomes.json under {d}")
    return c80._owners(uids, root, list(root_stages), owner_inputs, later_inputs, work, threads, memory_limit)


def _write_objects_json(
    sf: c78.StageFiles,
    used_categories: list[str],
    source_meta_by_cat: dict[str, dict],
    inputs_home: Path,
) -> None:
    objs = []
    objs_staging = []
    for k in used_categories:
        fname = sf.fasta_names[k]
        meta = source_meta_by_cat[k]
        objs.append({"path": str(Path(inputs_home) / sf.stage / fname), **meta})
        objs_staging.append({"path": str(sf.final_path(k)), **meta})
    (sf.tmp / "objects.json").write_text(json.dumps(objs, indent=1))
    (sf.tmp / "objects_staging.json").write_text(json.dumps(objs_staging, indent=1))


# ==================================================================================================================
# 4. Stage O4: Protist weird mitochondria, plastids, nucleomorphs (O4_ORG) + Fungal mtDNA (O3_FUNGAL)
# ==================================================================================================================

O4 = "O4"
O4_FASTA = {
    "ncbi": "o4_ncbi_organelles.fna.gz",
    "fungal": "o4_fungal_mitogenomes.fna.gz",
}
O4_SUBCOHORTS = ("O4_CRYP", "O4_PROT_MT", "O4_PROT_PT", "O4_WEIRD_MT")
O4_RETRIEVED_NCBI = "2026-11-26T19:06:55Z"
O4_RETRIEVED_FUNGAL = "2026-11-25T22:05:41Z"

O4_RECORDS_HEAD = [
    "ordinal", "source", "subcohort", "compartment", "genome", "accession", "organism",
    "contig_uid", "length", "circular", "genetic_code", "code_source", "n_cds",
    "n_spliced_cds", "n_introns", "n_homing_rt", "owner_stage", "later_holder",
]
O4_LEFT_HEAD = [
    "accession", "source", "subcohort", "compartment", "organism", "length", "reason",
]


def normalize_o4_compartment(organelle: str, definition: str) -> str:
    org_l = (organelle or "").lower()
    def_l = (definition or "").lower()
    if "nucleomorph" in org_l or "nucleomorph" in def_l:
        return "nucleomorph"
    if "chromatophore" in org_l or "chromatophore" in def_l or "cyanelle" in org_l or "cyanelle" in def_l:
        return "chromatophore"
    if any(w in org_l or w in def_l for w in ("plastid", "chloroplast", "apicoplast", "leucoplast", "rhodoplast", "chromoplast", "elaioplast", "etioplast", "proplastid")):
        return "plastid"
    if any(w in org_l or w in def_l for w in ("mitochondri", "kinetoplast", "maxicircle", "minicircle", "hydrogenosome", "mitosome")):
        return "mitochondrion"
    return "nuclear_or_unmarked"


def classify_o4_gbff_record(rec: dict, fna_seq: bytes | None, *, defer_screen: bool = False) -> dict:
    mol = (rec.get("moltype") or "").upper()
    comp = normalize_o4_compartment(rec.get("organelle", ""), rec.get("definition", ""))
    seq = fna_seq if fna_seq is not None else (rec.get("seq") or b"")
    tax = rec.get("taxonomy") or ""

    if "RNA" in mol:
        return {"keep": False, "compartment": comp, "reason": f"moltype_not_dna:{rec.get('moltype')}"}
    if not seq or not VALID_DNA_RE.match(seq):
        return {"keep": False, "compartment": comp, "reason": "non_iupac_or_empty"}
    if tax and "eukaryota" not in tax.lower():
        return {"keep": False, "compartment": comp, "reason": "not_eukaryota"}
    if comp == "nuclear_or_unmarked":
        return {"keep": False, "compartment": comp, "reason": "nuclear_non_organellar_c60"}
    n_cds = int(rec.get("n_cds") or 0)
    n_spl = int(rec.get("n_spliced_cds") or 0)
    if comp == "nucleomorph" and n_cds >= 5 and (n_spl / n_cds) > 0.5:
        return {"keep": False, "compartment": comp, "reason": "nuclear_spliced_intron_rich_c60"}

    tt = rec.get("transl_table")
    code_source = "gbff_transl_table"
    if tt is None:
        tax_l = tax.lower()
        if comp == "nucleomorph":
            tt = 1
            code_source = "nucleomorph_default_1"
        elif len(seq) < 5_000:
            if comp in ("plastid", "chromatophore"):
                tt = 4 if "apicomplexa" in tax_l else 11
            else:
                tt = 4 if any(w in tax_l for w in (
                    "euglenozoa", "kinetoplastea", "apicomplexa", "alveolata", "fungi",
                    "stramenopiles", "amoebozoa", "haptophyta", "rhizaria",
                )) else 1
            code_source = f"short_lt5kb_lineage_{tt}"
        elif defer_screen:
            return {
                "keep": True,
                "compartment": comp,
                "genetic_code": None,
                "code_source": "needs_screen",
            }
        elif comp in ("plastid", "chromatophore"):
            sc = screen_prokaryote_code([seq], tax, force_alt=True)
            tt = 4 if sc["genetic_code"] == 4 else 11
            code_source = f"screen_{sc['reason']}"
        else:
            sc = screen_fungal_mt_code(seq, tax)
            tt = sc["genetic_code"]
            code_source = f"screen_{sc['reason']}"
    try:
        orfs.check_code(int(tt))
    except ValueError:
        return {"keep": False, "compartment": comp, "reason": f"unsupported_transl_table_{tt}"}

    return {
        "keep": True,
        "compartment": comp,
        "genetic_code": int(tt),
        "code_source": code_source,
    }


def _o4_unann_screen_worker(args: tuple[str, bytes, str, str]) -> tuple[str, int, str]:
    chex, seq, comp, tax = args
    if comp in ("plastid", "chromatophore"):
        sc = screen_prokaryote_code([seq], tax, force_alt=True)
        tt = 4 if sc["genetic_code"] == 4 else 11
        return chex, tt, f"screen_{sc['reason']}"
    sc = screen_fungal_mt_code(seq, tax)
    return chex, int(sc["genetic_code"]), f"screen_{sc['reason']}"


def read_fungal_accessions_tsv(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if "seq_id" not in (rd.fieldnames or []):
            raise Refused(f"{path}: missing seq_id column")
        for r in rd:
            sid = (r.get("seq_id") or "").strip()
            if sid:
                out[sid] = r
    return out


def stage_o4(
    o4_org_dir: Path,
    fungal_fa: Path,
    fungal_tsv: Path,
    meta_out: Path,
    fasta_out: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    later_inputs: dict[str, Path],
    inputs_home: Path,
    retrieved_ncbi: str = O4_RETRIEVED_NCBI,
    retrieved_fungal: str = O4_RETRIEVED_FUNGAL,
    work: Path | None = None,
    threads: int = 2,
    memory_limit: str = "4GB",
    workers: int = 1,
) -> dict:
    o4_org_dir, fungal_fa, fungal_tsv = Path(o4_org_dir), Path(fungal_fa), Path(fungal_tsv)
    for p in (fungal_fa, fungal_tsv):
        if not p.is_file():
            raise Refused(f"missing input file {p}")
    fungal_meta = read_fungal_accessions_tsv(fungal_tsv)

    sf = c78.StageFiles(O4, Path(meta_out), Path(fasta_out), O4_FASTA, work)
    pool = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None
    try:
        staged_recs: list[dict] = []
        left_out: list[tuple] = []
        seen_contig_props: dict[bytes, tuple[int, bool, str]] = {}  # cu -> (code, circular, quality)
        by_subcohort: Counter[str] = Counter()
        by_compartment: Counter[str] = Counter()
        by_code: Counter[str] = Counter()

        # 1. Stream O4_ORG paired .gbff + .fna chunks in deterministic order
        chunk_pairs: list[tuple[str, Path, Path]] = []
        for sub in O4_SUBCOHORTS:
            sdir = o4_org_dir / sub
            if not sdir.is_dir():
                continue
            for gbff in sorted(sdir.glob(f"{sub}_chunk_*.gbff")):
                fna = gbff.with_suffix(".fna")
                if not fna.is_file():
                    raise Refused(f"missing paired .fna for {gbff}")
                chunk_pairs.append((sub, gbff, fna))

        unann_to_screen: dict[bytes, tuple[bytes, str, str]] = {}
        ncbi_pending: list[dict] = []
        for sub, gbff_p, fna_p in chunk_pairs:
            fna_map: dict[str, tuple[str, bytes, int]] = {}
            with open(fna_p, "rb") as ffh:
                for r in v2.iter_fasta_stream(ffh):
                    fna_map[r.accession] = (r.header, r.seq, r.line_width or 70)
                    fna_map.setdefault(c80.base_acc(r.accession), (r.header, r.seq, r.line_width or 70))
            with open(gbff_p, "r", encoding="utf-8", errors="replace") as gfh:
                for rec in c59.parse_refseq_organelle_gbff(gfh, keep_translations=False):
                    acc = rec["accession"]
                    fna_hit = fna_map.get(acc) or fna_map.get(rec["base_accession"])
                    fna_hdr = fna_hit[0] if fna_hit else f"{acc} {rec['definition']}".strip()
                    seq = fna_hit[1] if fna_hit else rec["seq"]
                    width = fna_hit[2] if fna_hit else 70
                    cls = classify_o4_gbff_record(rec, seq, defer_screen=True)
                    if not cls["keep"]:
                        left_out.append((
                            acc, "tower_batch6b_ncbi_o4_org", sub, cls["compartment"],
                            rec["organism"], len(seq), cls["reason"],
                        ))
                        continue
                    cu = ids.contig_uid(seq)
                    circ = bool(rec["circular"])
                    qual = "complete" if circ or "complete" in rec["definition"].lower() else (
                        "chromosome" if cls["compartment"] == "nucleomorph" else "contig"
                    )
                    if cls["genetic_code"] is not None and cu not in seen_contig_props:
                        seen_contig_props[cu] = (cls["genetic_code"], circ, qual)
                    elif cls["genetic_code"] is None and cu not in seen_contig_props and cu not in unann_to_screen:
                        unann_to_screen[cu] = (seq, cls["compartment"], rec["taxonomy"])
                    written = sf.write("ncbi", fna_hdr, seq, width, prefix=sub)
                    ncbi_pending.append({
                        "source": "tower_batch6b_ncbi_o4_org",
                        "category": "ncbi",
                        "subcohort": sub,
                        "compartment": cls["compartment"],
                        "genome": written,
                        "accession": acc,
                        "organism": rec["organism"],
                        "cu": cu,
                        "length": len(seq),
                        "circular": circ,
                        "code": cls["genetic_code"],
                        "code_source": cls["code_source"],
                        "quality": qual,
                        "n_cds": rec["n_cds"],
                        "n_spliced_cds": rec["n_spliced_cds"],
                        "n_introns": rec["n_introns"],
                        "n_homing_rt": rec["n_homing_rt"],
                    })
                    by_subcohort[sub] += 1

        unann_res: dict[bytes, tuple[int, str]] = {}
        if unann_to_screen:
            to_run = [(cu.hex(), sq, cp, tx) for cu, (sq, cp, tx) in sorted(unann_to_screen.items()) if cu not in seen_contig_props]
            if pool is not None and to_run:
                for fut in [pool.submit(_o4_unann_screen_worker, item) for item in to_run]:
                    chex, tt, csrc = fut.result()
                    unann_res[bytes.fromhex(chex)] = (tt, csrc)
            else:
                for item in to_run:
                    chex, tt, csrc = _o4_unann_screen_worker(item)
                    unann_res[bytes.fromhex(chex)] = (tt, csrc)

        for item in ncbi_pending:
            cu = item["cu"]
            if cu in seen_contig_props:
                tt, circ, qual = seen_contig_props[cu]
                if item["code"] is None or item["code"] != tt or item["circular"] != circ or item["quality"] != qual:
                    item["code_source"] = item["code_source"] if item["code"] == tt and item["code"] is not None else (
                        unann_res[cu][1] if cu in unann_res else "matched_earlier_contig_uid"
                    )
                item["code"], item["circular"], item["quality"] = tt, circ, qual
            else:
                tt, csrc = unann_res[cu]
                check_supported_code(tt)
                seen_contig_props[cu] = (tt, item["circular"], item["quality"])
                item["code"] = tt
                item["code_source"] = csrc
            staged_recs.append(item)

        # 2. Stream O3_FUNGAL mitogenomes
        fungal_raw: list[tuple[str, str, bytes, int, dict]] = []
        to_screen: dict[bytes, tuple[bytes, str]] = {}
        with open(fungal_fa, "rb") as ffh:
            for r in v2.iter_fasta_stream(ffh):
                sid = r.accession
                m = fungal_meta.get(sid, {})
                org = m.get("species") or m.get("genus") or sid
                if not r.seq or not VALID_DNA_RE.match(r.seq):
                    left_out.append((sid, "tower_o3_fungal", "O3_FUNGAL", "mitochondrion", org, len(r.seq), "non_iupac_or_empty"))
                    continue
                cu = ids.contig_uid(r.seq)
                tax = ";".join(m.get(k, "") for k in ("kingdom", "phylum", "class", "order", "family", "genus", "species"))
                if cu not in seen_contig_props and cu not in to_screen:
                    to_screen[cu] = (r.seq, tax)
                fungal_raw.append((sid, r.header, r.seq, r.line_width or 70, m))

        screened_fungal: dict[bytes, dict] = {}
        if to_screen:
            if pool is not None:
                futs = [
                    pool.submit(_fungal_screen_worker, (cu.hex(), sq, tx))
                    for cu, (sq, tx) in sorted(to_screen.items())
                ]
                for fut in futs:
                    chex, res = fut.result()
                    screened_fungal[bytes.fromhex(chex)] = res
            else:
                for cu, (sq, tx) in sorted(to_screen.items()):
                    screened_fungal[cu] = screen_fungal_mt_code(sq, tx)

        for sid, hdr, seq, width, m in fungal_raw:
            cu = ids.contig_uid(seq)
            circ = ("topology=circular" in hdr.lower()) or bool(CIRCULAR_RE.search(hdr))
            qual = "complete" if circ or "completeness=complete" in hdr.lower() else "contig"
            if cu in seen_contig_props:
                tt, circ, qual = seen_contig_props[cu]
                code_src = "matched_earlier_contig_uid"
            else:
                sc = screened_fungal[cu]
                tt = int(sc["genetic_code"])
                code_src = f"screen_{sc['reason']}"
                seen_contig_props[cu] = (tt, circ, qual)
            check_supported_code(tt)
            written = sf.write("fungal", hdr, seq, width, prefix="O3_FUNGAL")
            org = m.get("species") or m.get("genus") or sid
            staged_recs.append({
                "source": "tower_o3_fungal",
                "category": "fungal",
                "subcohort": "O3_FUNGAL",
                "compartment": "mitochondrion",
                "genome": written,
                "accession": m.get("NCBI_accession_number") if m.get("NCBI_accession_number") not in (None, "", "NA") else sid,
                "organism": org,
                "cu": cu,
                "length": len(seq),
                "circular": circ,
                "code": tt,
                "code_source": code_src,
                "quality": qual,
                "n_cds": 0,
                "n_spliced_cds": 0,
                "n_introns": 0,
                "n_homing_rt": 0,
            })
            by_subcohort["O3_FUNGAL"] += 1

        if not staged_recs:
            raise Refused("O4: no qualifying organellar records staged")

        sf.close()
        all_uids = {r["cu"] for r in staged_recs}
        owners, later = _owners_of_uids(
            all_uids, root, root_stages, owner_inputs, later_inputs, sf.work, threads, memory_limit
        )

        listing = make_inputs.Listing()
        genomes: list[dict] = []
        rec_rows: list[tuple] = []
        repeated: list[tuple] = []
        shared: list[tuple] = []
        held_later: list[tuple] = []
        first_of: dict[bytes, str] = {}
        contig_len: dict[bytes, int] = {}

        for ordinal, r in enumerate(staged_recs):
            cu = r["cu"]
            own = owners.get(cu, "")
            lat = "" if own else later.get(cu, "")
            rec_rows.append((
                ordinal, r["source"], r["subcohort"], r["compartment"], r["genome"], r["accession"],
                r["organism"], cu.hex(), r["length"], str(r["circular"]).lower(), r["code"],
                r["code_source"], r["n_cds"], r["n_spliced_cds"], r["n_introns"], r["n_homing_rt"],
                own, lat,
            ))
            if own:
                shared.append((ordinal, r["genome"], cu.hex(), r["length"], own))
            elif lat:
                held_later.append((ordinal, r["genome"], cu.hex(), r["length"], lat))
            contig_len.setdefault(cu, r["length"])
            g = listing.add([cu], r["code"], "prokaryote", r["quality"])
            guid = ids.genome_uid([cu])
            if g is None:
                repeated.append((ordinal, r["genome"], guid.hex(), first_of[guid]))
                continue
            g["circular"] = [cu.hex()] if r["circular"] else []
            genomes.append(g)
            first_of[guid] = r["genome"]
            by_compartment[r["compartment"]] += 1
            by_code[str(r["code"])] += 1

        make_inputs._write_tsv(sf.tmp / "records.tsv", O4_RECORDS_HEAD, rec_rows)
        make_inputs._write_tsv(sf.tmp / "left_out.tsv", O4_LEFT_HEAD, left_out)
        make_inputs._write_tsv(sf.tmp / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
        make_inputs._write_tsv(sf.tmp / "shared_records.tsv", ["ordinal", "genome", "contig_uid", "length", "owner_stage"], shared)
        make_inputs._write_tsv(sf.tmp / "later_records.tsv", ["ordinal", "genome", "contig_uid", "length", "holder"], held_later)
        (sf.tmp / "genomes.json").write_text(json.dumps(genomes))

        used = [k for k in O4_FASTA if any(r["category"] == k for r in staged_recs)]
        source_meta = {
            "ncbi": {
                "source_name": "tower_batch6b_ncbi_o4_org",
                "source_release": "NCBI nuccore protist/cryptophyte organellar & nucleomorph chunks (O4_CRYP, O4_PROT_MT, O4_PROT_PT, O4_WEIRD_MT); staged by c86_o4_gv2_hrom1_inputs.py (DEC-65 (2), DEC-85 (4))",
                "license": "INSDC",
                "redistributable": True,
                "retrieved_at": retrieved_ncbi,
            },
            "fungal": {
                "source_name": "tower_o3_fungal",
                "source_release": "Sarhan et al. 2025 assembled circular fungal mitochondrial genomes (figshare 28750034); genetic codes screened by c86_o4_gv2_hrom1_inputs.py (DEC-80, DEC-85 (4))",
                "license": "CC BY 4.0 as figshare record 28750034 api_metadata.json states it (sources.yaml tower_o3_fungal note: audit wish-list question unresolved)",
                "redistributable": False,
                "retrieved_at": retrieved_fungal,
            },
        }
        _write_objects_json(sf, used, source_meta, Path(inputs_home))
        sha_map, byte_map = sf.commit(used)

        own_c, own_b = _tally_owners(contig_len, owners)
        lat_c, lat_b = _tally_owners(contig_len, later, owners)
        return {
            "stage": O4,
            "records": len(staged_recs),
            "left_out": len(left_out),
            "left_out_by_reason": dict(sorted(Counter(x[6] for x in left_out).items())),
            "records_by_subcohort": dict(sorted(by_subcohort.items())),
            "genomes": len(genomes),
            "repeated_genomes": len(repeated),
            "circular_genomes": sum(1 for g in genomes if g["circular"]),
            "by_compartment": dict(sorted(by_compartment.items())),
            "by_genetic_code": dict(sorted(by_code.items())),
            "contigs_distinct": len(contig_len),
            "bases_distinct": sum(contig_len.values()),
            "held_by": own_c,
            "held_bases_by": own_b,
            "later_held_by": lat_c,
            "later_held_bases_by": lat_b,
            "contigs_new": sum(1 for cu in contig_len if cu not in owners),
            "bases_new": sum(n for cu, n in contig_len.items() if cu not in owners),
            "contigs_new_after_later": sum(1 for cu in contig_len if cu not in owners and cu not in later),
            "bases_new_after_later": sum(n for cu, n in contig_len.items() if cu not in owners and cu not in later),
            "sha256": sha_map,
            "bytes": byte_map,
            **{"listing_" + k: v for k, v in listing.report().items()},
        }
    except BaseException:
        sf.abort()
        raise
    finally:
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)


# ==================================================================================================================
# 5. Stage GV2: Polar, freshwater, alpine-lake, and WWTP giant viruses
# ==================================================================================================================

GV2 = "GV2"
GV2_FASTA = {"gv2": "gv2_giant_viruses.fna.gz"}
GV2_RETRIEVED = "2026-11-26T23:37:56Z"
GV2_SUBCOHORTS = (
    "GV_POLAR_2026",
    "GV_POLAR_PITOT",
    "GV_POLAR_GEVE",
    "GV_GKS",
    "GV_POLAR_MARINE",
    "GV_WWTP",
)
GV2_RECORDS_HEAD = [
    "ordinal", "subcohort", "member", "genome", "record_ordinal", "contig_uid",
    "length", "circular", "dtr_class", "source_contig_uid", "source_length",
    "trimmed_bases", "owner_stage", "later_holder",
]
GV2_LEFT_HEAD = ["subcohort", "member", "accession", "length", "reason"]
GV2_SKIP_MEMBER_HEAD = ["subcohort", "member", "bytes", "reason"]


def _find_in_dirs(roots: list[Path], filename: str) -> Path | None:
    for d in roots:
        p = Path(d) / filename
        if p.is_file():
            return p
    return None


def _all_files_in_dirs(roots: list[Path]) -> list[Path]:
    seen: dict[str, Path] = {}
    for d in roots:
        d = Path(d)
        if not d.is_dir():
            continue
        for p in sorted(d.iterdir()):
            if p.is_file() and p.name != "SHA256SUMS" and not p.name.startswith("."):
                seen.setdefault(p.name, p)
    return [seen[k] for k in sorted(seen)]


def _stem_clean(name: str) -> str:
    base = os.path.basename(name)
    for ext in (".gz", ".fasta", ".fna", ".fa", ".fsa", ".fltr"):
        if base.lower().endswith(ext):
            base = base[: -len(ext)]
    return base


def stage_gv2(
    dirs: dict[str, Path | list[Path]],
    meta_out: Path,
    fasta_out: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    later_inputs: dict[str, Path],
    inputs_home: Path,
    retrieved_at: str = GV2_RETRIEVED,
    work: Path | None = None,
    threads: int = 2,
    memory_limit: str = "4GB",
    exclude_gks_all_mags: bool = False,
    screen_code25: bool = False,
) -> dict:
    norm_dirs: dict[str, list[Path]] = {}
    for sub in GV2_SUBCOHORTS:
        if sub not in dirs:
            raise Refused(f"GV2: missing directory mapping for {sub}")
        v = dirs[sub]
        norm_dirs[sub] = [Path(x) for x in v] if isinstance(v, (list, tuple)) else [Path(v)]

    sf = c78.StageFiles(GV2, Path(meta_out), Path(fasta_out), GV2_FASTA, work)
    try:
        genomes_raw: list[dict] = []
        trims: list[tuple] = []
        left_out: list[tuple] = []
        skipped_members: list[tuple] = []
        dtr_classes: Counter[str] = Counter()
        by_subcohort: Counter[str] = Counter()
        used_genome_names: set[str] = set()
        contig_topology: dict[bytes, bool] = {}  # enforce consistent circularity across shared contigs
        code25_hits: list[dict] = []

        def unique_gname(sub: str, stem: str) -> str:
            cand = stem if stem and stem not in used_genome_names else f"{sub}_{stem}"
            k = 2
            while cand in used_genome_names:
                cand = f"{sub}_{stem}_{k}"
                k += 1
            used_genome_names.add(cand)
            return cand

        def add_bin_genome(sub: str, member: str, stem: str, fobj, *, apply_dtr_if_single: bool = False) -> None:
            raw_recs: list[tuple[str, bytes, int]] = []
            for r in v2.iter_fasta_stream(fobj):
                if not r.seq or not VALID_DNA_RE.match(r.seq):
                    left_out.append((sub, member, r.accession, len(r.seq), "non_iupac_or_empty"))
                    continue
                raw_recs.append((r.header, r.seq, r.line_width or 70))
            if not raw_recs:
                skipped_members.append((sub, member, 0, "no_valid_dna_records"))
                return
            gname = unique_gname(sub, stem)
            contigs_info: list[tuple[bytes, int, int, str, str, bytes, int, bool]] = []
            seqs_for_c25: list[bytes] = []
            for hdr, seq, width in raw_recs:
                src_cu = ids.contig_uid(seq)
                if apply_dtr_if_single and len(raw_recs) == 1:
                    out_seq, k_trim, cls, circ = c78.dec65(seq, flagged=bool(CIRCULAR_RE.search(hdr)))
                    dtr_classes[cls] += 1
                else:
                    out_seq, k_trim, cls, circ = seq, 0, "", False
                cu = ids.contig_uid(out_seq) if k_trim else src_cu
                if cu in contig_topology:
                    circ = contig_topology[cu]
                else:
                    contig_topology[cu] = circ
                first = hdr.split(None, 1)[0] if hdr.strip() else ""
                prefixed_hdr = hdr if first.startswith(gname + "_") or first == gname else f"{gname}|{hdr}"
                written_name = sf.write("gv2", prefixed_hdr, out_seq, width, prefix=gname)
                if k_trim:
                    trims.append((sub, gname, written_name, src_cu.hex(), cu.hex(), len(seq), k_trim, k_trim, len(out_seq)))
                contigs_info.append((cu, len(out_seq), len(seq), written_name, cls, src_cu, k_trim, circ))
                seqs_for_c25.append(out_seq)
            if screen_code25:
                c25 = check_gv_code25(seqs_for_c25)
                if c25["use_c25"]:
                    code25_hits.append({"genome": gname, "subcohort": sub, **c25})
            all_circ = len(contigs_info) == 1 and contigs_info[0][7]
            genomes_raw.append({
                "subcohort": sub,
                "member": member,
                "name": gname,
                "quality": "complete" if all_circ else "contig",
                "circular": [c[0] for c in contigs_info if c[7]],
                "contigs": contigs_info,
            })
            by_subcohort[sub] += 1

        def add_single_record_catalog(sub: str, member: str, fobj) -> None:
            for r in v2.iter_fasta_stream(fobj):
                if not r.seq or not VALID_DNA_RE.match(r.seq):
                    left_out.append((sub, member, r.accession, len(r.seq), "non_iupac_or_empty"))
                    continue
                src_cu = ids.contig_uid(r.seq)
                out_seq, k_trim, cls, circ = c78.dec65(r.seq, flagged=bool(CIRCULAR_RE.search(r.header)))
                dtr_classes[cls] += 1
                cu = ids.contig_uid(out_seq) if k_trim else src_cu
                if cu in contig_topology:
                    circ = contig_topology[cu]
                else:
                    contig_topology[cu] = circ
                gname = unique_gname(sub, r.accession)
                written_name = sf.write("gv2", r.header, out_seq, r.line_width or 70, prefix=sub)
                if k_trim:
                    trims.append((sub, gname, written_name, src_cu.hex(), cu.hex(), len(r.seq), k_trim, k_trim, len(out_seq)))
                if screen_code25:
                    c25 = check_gv_code25([out_seq])
                    if c25["use_c25"]:
                        code25_hits.append({"genome": gname, "subcohort": sub, **c25})
                genomes_raw.append({
                    "subcohort": sub,
                    "member": member,
                    "name": gname,
                    "quality": "complete" if circ else "contig",
                    "circular": [cu] if circ else [],
                    "contigs": [(cu, len(out_seq), len(r.seq), written_name, cls, src_cu, k_trim, circ)],
                })
                by_subcohort[sub] += 1

        # 1. GV_POLAR_2026: GVMAGs_from_polar_regions.zip (skip 1,565 .proteins.faa files!)
        z2026 = _find_in_dirs(norm_dirs["GV_POLAR_2026"], "GVMAGs_from_polar_regions.zip")
        if z2026 is None:
            raise Refused("GV_POLAR_2026: missing GVMAGs_from_polar_regions.zip")
        with zipfile.ZipFile(z2026) as zf:
            for info in sorted(zf.infolist(), key=lambda x: x.filename):
                if info.is_dir():
                    continue
                if not _looks_dna_fasta(info.filename) or info.filename.startswith("__MACOSX/"):
                    why = "protein_fasta" if info.filename.lower().endswith((".faa", ".faa.gz")) else "not_fasta"
                    skipped_members.append(("GV_POLAR_2026", info.filename, info.file_size, why))
                    continue
                with zf.open(info) as raw:
                    add_bin_genome("GV_POLAR_2026", info.filename, _stem_clean(info.filename), raw)

        # 2. GV_POLAR_PITOT: GVMAGs_fna.tar.gz (+ optional tree .txt skipped)
        for p in _all_files_in_dirs(norm_dirs["GV_POLAR_PITOT"]):
            if p.name == "GVMAGs_fna.tar.gz":
                with tarfile.open(p, "r:gz") as tf:
                    members = sorted([m for m in tf.getmembers() if m.isfile()], key=lambda x: x.name)
                    for m in members:
                        if not _looks_dna_fasta(m.name) or os.path.basename(m.name).startswith("."):
                            skipped_members.append(("GV_POLAR_PITOT", m.name, m.size, "not_fasta"))
                            continue
                        ef = tf.extractfile(m)
                        assert ef is not None
                        add_bin_genome("GV_POLAR_PITOT", m.name, _stem_clean(m.name), ef)
            else:
                skipped_members.append(("GV_POLAR_PITOT", p.name, p.stat().st_size, "not_fasta"))

        # 3. GV_POLAR_GEVE: excised confirmed viral regions (*.fna) kept; all_ICEL_viral_contigs.zip host nuclear contigs left out
        for p in _all_files_in_dirs(norm_dirs["GV_POLAR_GEVE"]):
            if p.name == "all_ICEL_viral_contigs.zip":
                with zipfile.ZipFile(p) as zf:
                    for info in sorted(zf.infolist(), key=lambda x: x.filename):
                        if info.is_dir():
                            continue
                        if _looks_dna_fasta(info.filename):
                            with zf.open(info) as raw:
                                for r in v2.iter_fasta_stream(raw):
                                    left_out.append((
                                        "GV_POLAR_GEVE", f"{p.name}:{info.filename}", r.accession,
                                        len(r.seq), "eukaryote_host_nuclear_contig_c60",
                                    ))
                        else:
                            skipped_members.append(("GV_POLAR_GEVE", f"{p.name}:{info.filename}", info.file_size, "not_fasta"))
            elif _looks_dna_fasta(p.name):
                with open(p, "rb") as fh:
                    add_single_record_catalog("GV_POLAR_GEVE", p.name, fh)
            else:
                skipped_members.append(("GV_POLAR_GEVE", p.name, p.stat().st_size, "not_fasta"))

        # 4. GV_GKS: GKS_genomes.tar.gz (tomas_scMAGs kept; all_MAGs excluded if prokaryotic co-binned MAGs)
        zgks = _find_in_dirs(norm_dirs["GV_GKS"], "GKS_genomes.tar.gz")
        if zgks is None:
            raise Refused("GV_GKS: missing GKS_genomes.tar.gz")
        with tarfile.open(zgks, "r:gz") as tf:
            members = sorted([m for m in tf.getmembers() if m.isfile()], key=lambda x: x.name)
            for m in members:
                if not _looks_dna_fasta(m.name):
                    skipped_members.append(("GV_GKS", m.name, m.size, "not_fasta"))
                    continue
                ef = tf.extractfile(m)
                assert ef is not None
                if exclude_gks_all_mags and "/all_MAGs/" in m.name:
                    bp = sum(len(r.seq) for r in v2.iter_fasta_stream(ef))
                    left_out.append(("GV_GKS", m.name, _stem_clean(m.name), bp, "prokaryote_co_binned_mag_not_gv"))
                    continue
                add_bin_genome("GV_GKS", m.name, _stem_clean(m.name), ef)

        # 5. GV_POLAR_MARINE: individual .fa GVMAG files
        for p in _all_files_in_dirs(norm_dirs["GV_POLAR_MARINE"]):
            if _looks_dna_fasta(p.name):
                with open(p, "rb") as fh:
                    add_bin_genome("GV_POLAR_MARINE", p.name, _stem_clean(p.name), fh)
            else:
                skipped_members.append(("GV_POLAR_MARINE", p.name, p.stat().st_size, "not_fasta"))

        # 6. GV_WWTP: genomes.tar.gz (only genomes/ncv/contigs/* and genomes/vph_plv/contigs/* are assemblies!)
        for p in _all_files_in_dirs(norm_dirs["GV_WWTP"]):
            if p.name == "genomes.tar.gz":
                with tarfile.open(p, "r:gz") as tf:
                    members = sorted([m for m in tf.getmembers() if m.isfile()], key=lambda x: x.name)
                    for m in members:
                        if m.name.startswith("genomes/ncv/contigs/") and _looks_dna_fasta(m.name):
                            ef = tf.extractfile(m)
                            assert ef is not None
                            data = ef.read()
                            recs_tmp = list(v2.iter_fasta_stream(io.BytesIO(data)))
                            add_bin_genome("GV_WWTP", m.name, _stem_clean(m.name), io.BytesIO(data),
                                           apply_dtr_if_single=(len(recs_tmp) == 1))
                        elif m.name.startswith("genomes/vph_plv/contigs/") and _looks_dna_fasta(m.name):
                            ef = tf.extractfile(m)
                            assert ef is not None
                            data = ef.read()
                            recs_tmp = list(v2.iter_fasta_stream(io.BytesIO(data)))
                            if len(recs_tmp) == 1:
                                add_bin_genome("GV_WWTP", m.name, _stem_clean(m.name), io.BytesIO(data), apply_dtr_if_single=True)
                            else:
                                add_single_record_catalog("GV_WWTP", m.name, io.BytesIO(data))
                        else:
                            reason = "prokka_or_protein_annotation" if ("/annotations/" in m.name or "/proteins/" in m.name) else "not_genome_contig"
                            skipped_members.append(("GV_WWTP", m.name, m.size, reason))
            else:
                skipped_members.append(("GV_WWTP", p.name, p.stat().st_size, "not_fasta"))

        if not genomes_raw:
            raise Refused("GV2: no viral genomes staged")

        sf.close()
        all_uids = {c[0] for g in genomes_raw for c in g["contigs"]}
        owners, later = _owners_of_uids(
            all_uids, root, root_stages, owner_inputs, later_inputs, sf.work, threads, memory_limit
        )

        listing = make_inputs.Listing()
        out_genomes: list[dict] = []
        rec_rows: list[tuple] = []
        repeated: list[tuple] = []
        shared: list[tuple] = []
        held_later: list[tuple] = []
        first_of: dict[bytes, str] = {}
        contig_len: dict[bytes, int] = {}
        ordinal = 0

        for g in genomes_raw:
            uids = [c[0] for c in g["contigs"]]
            for k, (cu, n, src_n, wname, cls, src_cu, k_trim, circ) in enumerate(g["contigs"]):
                own = owners.get(cu, "")
                lat = "" if own else later.get(cu, "")
                rec_rows.append((
                    ordinal, g["subcohort"], g["member"], g["name"], k, cu.hex(),
                    n, str(circ).lower(), cls, src_cu.hex(), src_n, k_trim, own, lat,
                ))
                if own:
                    shared.append((ordinal, g["name"], cu.hex(), n, own))
                elif lat:
                    held_later.append((ordinal, g["name"], cu.hex(), n, lat))
                contig_len.setdefault(cu, n)
                ordinal += 1
            entry = listing.add(uids, 11, "virus", g["quality"])
            guid = ids.genome_uid(list(dict.fromkeys(uids)))
            if entry is None:
                repeated.append((ordinal - len(uids), g["name"], guid.hex(), first_of[guid]))
                continue
            entry["circular"] = [u.hex() for u in g["circular"]]
            out_genomes.append(entry)
            first_of[guid] = g["name"]

        make_inputs._write_tsv(sf.tmp / "records.tsv", GV2_RECORDS_HEAD, rec_rows)
        make_inputs._write_tsv(sf.tmp / "left_out.tsv", GV2_LEFT_HEAD, left_out)
        make_inputs._write_tsv(sf.tmp / "skipped_members.tsv", GV2_SKIP_MEMBER_HEAD, skipped_members)
        make_inputs._write_tsv(sf.tmp / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
        make_inputs._write_tsv(sf.tmp / "shared_records.tsv", ["ordinal", "genome", "contig_uid", "length", "owner_stage"], shared)
        make_inputs._write_tsv(sf.tmp / "later_records.tsv", ["ordinal", "genome", "contig_uid", "length", "holder"], held_later)
        make_inputs._write_tsv(sf.tmp / "trims.tsv", c62.TRIM_HEAD, trims)
        (sf.tmp / "genomes.json").write_text(json.dumps(out_genomes))

        gv2_source = {
            "gv2": {
                "source_name": "tower_gv2",
                "source_release": "Polar, freshwater, alpine-lake and WWTP giant virus cohorts (GV_POLAR_2026 figshare:32316834, GV_POLAR_PITOT figshare:25000571, GV_POLAR_GEVE figshare:29546846, GV_GKS zenodo:22834836, GV_POLAR_MARINE figshare:26314543, GV_WWTP zenodo:17735769); staged by c86_o4_gv2_hrom1_inputs.py (DEC-65 (1), DEC-80, DEC-88)",
                "license": "CC BY 4.0 / cc-by-4.0 as the job lists state them (sources.yaml tower_gv_*); no separate provider record kept",
                "redistributable": False,
                "retrieved_at": retrieved_at,
            }
        }
        _write_objects_json(sf, ["gv2"], gv2_source, Path(inputs_home))
        sha_map, byte_map = sf.commit(["gv2"])

        own_c, own_b = _tally_owners(contig_len, owners)
        lat_c, lat_b = _tally_owners(contig_len, later, owners)
        return {
            "stage": GV2,
            "genomes_read": len(genomes_raw),
            "genomes_by_subcohort": dict(sorted(by_subcohort.items())),
            "genomes": len(out_genomes),
            "repeated_genomes": len(repeated),
            "circular_genomes": sum(1 for g in out_genomes if g["circular"]),
            "dtr_classes": dict(sorted(dtr_classes.items())),
            "trims": len(trims),
            "trimmed_bases": sum(t[7] for t in trims),
            "left_out": len(left_out),
            "left_out_by_reason": dict(sorted(Counter(x[4] for x in left_out).items())),
            "skipped_members": len(skipped_members),
            "code25_genomes": len(code25_hits),
            "code25_examples": code25_hits[:10],
            "recommended_genecall_mode": "gv-meta-c25" if code25_hits else "gv-meta",
            "records": len(rec_rows),
            "contigs_distinct": len(contig_len),
            "bases_distinct": sum(contig_len.values()),
            "held_by": own_c,
            "held_bases_by": own_b,
            "later_held_by": lat_c,
            "later_held_bases_by": lat_b,
            "contigs_new": sum(1 for cu in contig_len if cu not in owners),
            "bases_new": sum(n for cu, n in contig_len.items() if cu not in owners),
            "contigs_new_after_later": sum(1 for cu in contig_len if cu not in owners and cu not in later),
            "bases_new_after_later": sum(n for cu, n in contig_len.items() if cu not in owners and cu not in later),
            "sha256": sha_map,
            "bytes": byte_map,
            **{"listing_" + k: v for k, v in listing.report().items()},
        }
    except BaseException:
        sf.abort()
        raise


# ==================================================================================================================
# 6. Stage HROM1: Oral/CPR species reps (HROM_ORAL) + code-4 sets (CODE4_GLMC + CODE4_ACTINO1)
# ==================================================================================================================

HROM1 = "HROM1"
HROM1_FASTA = {
    "oral_glmc": "hrom1_oral_glmc.fna.gz",
    "actino1": "hrom1_code4_actino1.fna.gz",
}
HROM1_RETRIEVED_ORAL = "2026-11-27T06:20:21Z"
HROM1_RETRIEVED_ACTINO = "2026-11-26T19:15:16Z"

HROM1_RECORDS_HEAD = [
    "ordinal", "subcohort", "genome", "record_ordinal", "contig_uid", "length",
    "assembly_quality", "genetic_code", "circular", "owner_stage", "later_holder",
]
HROM1_STATUS_HEAD = [
    "genome", "subcohort", "accession", "r232_status", "in_gtdb_r232",
    "is_gtdb_representative", "gtdb_representative", "gtdb_species",
    "in_globdb_r232", "globdb_r232_id", "taxonomy", "completeness",
    "contamination", "species_status", "isolate_status", "contigs", "bp",
]
HROM1_LEFT_HEAD = ["subcohort", "genome", "member", "length", "reason"]


def read_hrom_species_metadata(path: Path) -> dict[str, dict]:
    """Parse HROM-Species-metadata.tsv by positional index because 'species' appears twice in the header."""
    out: dict[str, dict] = {}
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        head = [c.strip() for c in (next(rd, None) or [])]
        if len(head) < 15 or head[0] != "species" or "gtdb_taxonomy" not in head or "species_status" not in head:
            raise Refused(f"{path}: unexpected HROM metadata header {head[:8]!r}")
        i_comp = head.index("completeness")
        i_cont = head.index("contamination")
        i_tax = head.index("gtdb_taxonomy")
        i_stat = head.index("species_status")
        i_iso = head.index("isolate status") if "isolate status" in head else -1
        for row in rd:
            if not row or not row[0].strip():
                continue
            gid = row[0].strip()
            if gid in out:
                raise Refused(f"{path}: duplicate HROM genome id {gid}")
            out[gid] = {
                "completeness": row[i_comp].strip() if i_comp < len(row) else "",
                "contamination": row[i_cont].strip() if i_cont < len(row) else "",
                "taxonomy": row[i_tax].strip() if i_tax < len(row) else "",
                "species_status": row[i_stat].strip() if i_stat < len(row) else "",
                "isolate_status": row[i_iso].strip() if 0 <= i_iso < len(row) else "",
            }
    return out


def read_glmc_metadata_csv(path: Path) -> dict[str, dict]:
    """Parse CODE4_GLMC metadata_mags.csv (UTF-8 with BOM)."""
    out: dict[str, dict] = {}
    with open(path, encoding="utf-8-sig", errors="replace", newline="") as fh:
        rd = csv.DictReader(fh)
        if "Short_name" not in (rd.fieldnames or []):
            raise Refused(f"{path}: missing Short_name column")
        for r in rd:
            gid = (r.get("Short_name") or "").strip()
            if gid:
                out[gid] = {k.strip(): (v or "").strip() for k, v in r.items() if k}
    return out


def read_actino1_manifest(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        fields = rd.fieldnames or []
        col = "accession" if "accession" in fields else ("assembly_accession" if "assembly_accession" in fields else None)
        if col is None:
            raise Refused(f"{path}: missing accession / assembly_accession column")
        for r in rd:
            acc = (r.get(col) or "").strip()
            if acc:
                out[acc] = {k.strip(): (v or "").strip() for k, v in r.items() if k}
    return out


def _quality_from_checkm(comp_s: str, cont_s: str, status_s: str = "") -> str:
    st = status_s.strip().upper()
    if "HQ" in st or st == "NC/HQ":
        return "mag_hq"
    if st == "MQ":
        return "mag_mq"
    try:
        comp = float(comp_s)
        cont = float(cont_s)
        if comp >= 90.0 and cont <= 5.0:
            return "mag_hq"
        if comp >= 50.0 and cont <= 10.0:
            return "mag_mq"
    except ValueError:
        pass
    return "contig"


def _r232_label(x: dict | None, has_acc: bool) -> str:
    if not has_acc or x is None:
        return "not_visible"
    if x["in_gtdb_r232"] and x["in_globdb_r232"]:
        return "in_gtdb_and_globdb_r232"
    if x["in_gtdb_r232"]:
        return "in_gtdb_r232"
    if x["in_globdb_r232"]:
        return "in_globdb_r232"
    return "not_in_r232"


def stage_hrom1(
    hrom_tar: Path,
    hrom_meta: Path,
    glmc_zip: Path,
    glmc_csv: Path,
    actino_dir: Path,
    index_dir: Path,
    meta_out: Path,
    fasta_out: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    later_inputs: dict[str, Path],
    inputs_home: Path,
    hrom_hq_only: bool = False,
    retrieved_oral: str = HROM1_RETRIEVED_ORAL,
    retrieved_actino: str = HROM1_RETRIEVED_ACTINO,
    work: Path | None = None,
    threads: int = 2,
    memory_limit: str = "4GB",
    workers: int = 1,
) -> dict:
    hrom_tar, hrom_meta = Path(hrom_tar), Path(hrom_meta)
    glmc_zip, glmc_csv = Path(glmc_zip), Path(glmc_csv)
    actino_dir, index_dir = Path(actino_dir), Path(index_dir)
    if not (index_dir / "gtdb_r232.parquet").is_file() or not (index_dir / "globdb_r232.parquet").is_file():
        raise Refused(f"{index_dir} is not card C-74's index (missing gtdb_r232.parquet or globdb_r232.parquet)")

    hmeta = read_hrom_species_metadata(hrom_meta)
    gmeta = read_glmc_metadata_csv(glmc_csv)
    act_man_path = actino_dir / "MANIFEST_CODE4_ACTINO1_PRJNA1366607.tsv"
    ameta = read_actino1_manifest(act_man_path)

    sf = c78.StageFiles(HROM1, Path(meta_out), Path(fasta_out), HROM1_FASTA, work)
    pool = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None
    pending: list[tuple[str, object]] = []
    try:
        gc_results: dict[str, dict] = {}
        contigs: dict[str, tuple[bytes, list[int]]] = {}
        circ_uids: dict[str, list[bytes]] = {}
        ginfo: dict[str, dict] = {}
        left_out: list[tuple] = []
        order_gids: list[str] = []

        def submit_screen(gid: str, seqs: list[bytes], tax: str, force_alt: bool) -> None:
            if pool is None:
                gc_results[gid] = screen_prokaryote_code(seqs, tax, force_alt=force_alt)
                return
            pending.append((gid, pool.submit(_prok_screen_worker, (gid, seqs, tax, force_alt))))
            while len(pending) >= 4 * workers:
                s, fut = pending.pop(0)
                gc_results[s] = fut.result()[1]

        def ingest_genome(
            subcohort: str,
            category: str,
            gid: str,
            member: str,
            fobj,
            *,
            tax: str,
            quality: str,
            acc: str,
            comp: str,
            cont: str,
            sp_status: str,
            iso_status: str,
            force_alt: bool,
            allow_complete_circ: bool = False,
        ) -> None:
            recs: list[tuple[str, bytes, int]] = []
            for r in v2.iter_fasta_stream(fobj):
                if not r.seq or not VALID_DNA_RE.match(r.seq):
                    left_out.append((subcohort, gid, f"{member}:{r.accession}", len(r.seq), "non_iupac_or_empty"))
                    continue
                recs.append((r.header, r.seq, r.line_width or 70))
            if not recs:
                left_out.append((subcohort, gid, member, 0, "no_valid_dna_records"))
                return
            if gid in contigs:
                raise Refused(f"HROM1: duplicate genome identifier {gid}")
            cus: list[bytes] = []
            c_circ: list[bytes] = []
            for hdr, seq, width in recs:
                first = hdr.split(None, 1)[0] if hdr.strip() else ""
                prefixed = hdr if first.startswith(gid + "_") or first == gid else f"{gid}|{hdr}"
                sf.write(category, prefixed, seq, width, prefix=gid)
                cu = ids.contig_uid(seq)
                cus.append(cu)
                if allow_complete_circ and (
                    CIRCULAR_RE.search(hdr)
                    or (quality == "complete" and len(recs) == 1)
                    or ("complete genome" in hdr.lower() and "plasmid" not in hdr.lower())
                    or ("complete sequence" in hdr.lower() and "plasmid" in hdr.lower())
                ):
                    c_circ.append(cu)
                if not acc:
                    m_acc = ACC_GCA_RE.search(hdr)
                    if m_acc:
                        acc = m_acc.group(0)
            contigs[gid] = (b"".join(cus), [len(s) for _, s, _ in recs])
            circ_uids[gid] = c_circ
            ginfo[gid] = {
                "subcohort": subcohort,
                "category": category,
                "quality": quality,
                "accession": acc,
                "taxonomy": tax,
                "completeness": comp,
                "contamination": cont,
                "species_status": sp_status,
                "isolate_status": iso_status,
            }
            order_gids.append(gid)
            submit_screen(gid, [s for _, s, _ in recs], tax, force_alt)

        # 1. Stream HROM_ORAL tar.gz
        with tarfile.open(hrom_tar, "r:gz") as tf:
            for m in tf:
                if not m.isfile() or not _looks_fasta(m.name):
                    continue
                gid = _stem_clean(m.name)
                if gid not in hmeta:
                    raise Refused(f"HROM_ORAL: archive member {m.name} ({gid}) missing from {hrom_meta}")
                hm = hmeta[gid]
                if hrom_hq_only and hm["species_status"] != "NC/HQ":
                    ef = tf.extractfile(m)
                    assert ef is not None
                    bp = sum(len(r.seq) for r in v2.iter_fasta_stream(ef))
                    left_out.append(("HROM_ORAL", gid, m.name, bp, f"hrom_non_hq:{hm['species_status']}"))
                    continue
                ef = tf.extractfile(m)
                assert ef is not None
                qual = _quality_from_checkm(hm["completeness"], hm["contamination"], hm["species_status"])
                m_acc = ACC_GCA_RE.search(hm["taxonomy"])
                ingest_genome(
                    "HROM_ORAL", "oral_glmc", gid, m.name, ef,
                    tax=hm["taxonomy"], quality=qual, acc=m_acc.group(0) if m_acc else "",
                    comp=hm["completeness"], cont=hm["contamination"],
                    sp_status=hm["species_status"], iso_status=hm["isolate_status"],
                    force_alt=False, allow_complete_circ=False,
                )

        # 2. Stream CODE4_GLMC merged_MAGs.zip
        with zipfile.ZipFile(glmc_zip) as zf:
            for info in sorted(zf.infolist(), key=lambda x: x.filename):
                if info.is_dir() or not _looks_dna_fasta(info.filename) or info.filename.startswith("__MACOSX/"):
                    continue
                gid = _stem_clean(info.filename)
                gm = gmeta.get(gid, {})
                tax = gm.get("Classification_GTDB_226") or gm.get("Taxonomy", "")
                comp = gm.get("checkm2_completeness") or gm.get("Completeness", "")
                cont = gm.get("checkm2_contamination") or gm.get("Contamination", "")
                mimag = gm.get("MIMAG", "")
                sample = gm.get("Sample", "")
                m_acc = ACC_GCA_RE.search(sample) or ACC_GCA_RE.search(gid)
                qual = _quality_from_checkm(comp, cont, mimag)
                with zf.open(info) as raw:
                    ingest_genome(
                        "CODE4_GLMC", "oral_glmc", gid, info.filename, raw,
                        tax=tax, quality=qual, acc=m_acc.group(0) if m_acc else "",
                        comp=comp, cont=cont, sp_status=mimag, iso_status=sample,
                        force_alt=True, allow_complete_circ=False,
                    )

        # 3. Stream CODE4_ACTINO1 assemblies/*.fna.gz in manifest / sorted order
        asm_dir = actino_dir / "assemblies"
        for acc in sorted(ameta):
            matches = sorted(asm_dir.glob(f"{acc}*_genomic.fna.gz")) + sorted(asm_dir.glob(f"{acc}*.fna.gz"))
            if not matches:
                raise Refused(f"CODE4_ACTINO1: missing assembly FASTA for {acc} in {asm_dir}")
            fna_gz = matches[0]
            am = ameta[acc]
            lvl = am.get("assembly_level", "")
            qual = c58._ASSEMBLY_QUALITY_MAP.get(lvl, "contig")
            org = am.get("organism_name", "")
            with gzip.open(fna_gz, "rb") as gfh:
                ingest_genome(
                    "CODE4_ACTINO1", "actino1", acc, fna_gz.name, gfh,
                    tax=f"d__Bacteria;p__Actinomycetota;{org}", quality=qual, acc=acc,
                    comp="", cont="", sp_status=lvl, iso_status="NCBI",
                    force_alt=True, allow_complete_circ=True,
                )

        while pending:
            s, fut = pending.pop(0)
            gc_results[s] = fut.result()[1]

        if not order_gids:
            raise Refused("HROM1: no genomes staged")

        sf.close()
        kept = sorted(order_gids)
        resolved, genome_uids, conflicts = c58._resolve_within_stage_code_conflicts(kept, contigs, gc_results)

        # Enforce consistent circularity across any contig shared within HROM1
        shared_circ_any: set[bytes] = {cu for gid in kept for cu in circ_uids[gid]}

        all_uids = {
            contigs[g][0][i : i + 16]
            for g in kept
            for i in range(0, len(contigs[g][0]), 16)
        }
        owners, later = _owners_of_uids(
            all_uids, root, root_stages, owner_inputs, later_inputs, sf.work, threads, memory_limit
        )

        # Lookup GTDB R232 & GlobDB r232 membership (recorded, never cut!)
        accs_to_check = sorted({
            r232_index.nover(ginfo[g]["accession"])
            for g in kept
            if ginfo[g]["accession"]
        })
        look_map = {
            r232_index.nover(x["accession"]): x
            for x in (r232_index.lookup(index_dir, accs_to_check) if accs_to_check else [])
        }

        listing = make_inputs.Listing()
        genomes: list[dict] = []
        rec_rows: list[tuple] = []
        gc_rows: list[tuple] = []
        status_rows: list[tuple] = []
        repeated: list[tuple] = []
        shared: list[tuple] = []
        held_later: list[tuple] = []
        first_of: dict[bytes, str] = {}
        contig_len: dict[bytes, int] = {}
        by_subcohort: Counter[str] = Counter()
        status_counts: Counter[str] = Counter()
        ordinal = 0

        for g in kept:
            raw, lens = contigs[g]
            uids = [raw[i : i + 16] for i in range(0, len(raw), 16)]
            info = ginfo[g]
            gc = gc_results[g]
            code = resolved[g]
            gc_rows.append((
                g, HROM1, info["taxonomy"], sum(lens), len(lens), info["quality"], code,
                gc["reason"], str(gc["tested_alt"]).lower(),
                c58._fmt_opt(gc["cd_11"]), c58._fmt_opt(gc["mean_orf_11"]), c58._fmt_opt(gc["n_genes_11"]),
                c58._fmt_opt(gc["cd_25"]), c58._fmt_opt(gc["mean_orf_25"]), c58._fmt_opt(gc["n_genes_25"]),
                c58._fmt_opt(gc["cd_4"]), c58._fmt_opt(gc["mean_orf_4"]), c58._fmt_opt(gc["n_genes_4"]),
                c58._fmt_opt(gc["tga_per_kb"]),
            ))

            acc = info["accession"]
            x = look_map.get(r232_index.nover(acc)) if acc else None
            st_lbl = _r232_label(x, bool(acc))
            status_counts[st_lbl] += 1
            status_rows.append((
                g, info["subcohort"], acc, st_lbl,
                "yes" if (x and x["in_gtdb_r232"]) else ("no" if acc else ""),
                str(x["is_representative"]) if (x and x["in_gtdb_r232"]) else "",
                (x["representative"] or "") if x else "",
                (x["species"] or "") if x else "",
                "yes" if (x and x["in_globdb_r232"]) else ("no" if acc else ""),
                (x["globdb_r232_id"] or "") if x else "",
                info["taxonomy"], info["completeness"], info["contamination"],
                info["species_status"], info["isolate_status"], len(lens), sum(lens),
            ))

            g_circ = [u for u in uids if u in shared_circ_any]
            for k, (cu, n) in enumerate(zip(uids, lens)):
                own = owners.get(cu, "")
                lat = "" if own else later.get(cu, "")
                is_c = cu in shared_circ_any
                rec_rows.append((
                    ordinal, info["subcohort"], g, k, cu.hex(), n,
                    info["quality"], code, str(is_c).lower(), own, lat,
                ))
                if own:
                    shared.append((ordinal, g, k, cu.hex(), n, own))
                elif lat:
                    held_later.append((ordinal, g, k, cu.hex(), n, lat))
                contig_len.setdefault(cu, n)
                ordinal += 1

            e = listing.add(uids, code, "prokaryote", info["quality"])
            if e is None:
                repeated.append((ordinal - len(uids), g, genome_uids[g].hex(), first_of[genome_uids[g]]))
                continue
            e["circular"] = [u.hex() for u in dict.fromkeys(g_circ)]
            genomes.append(e)
            first_of[genome_uids[g]] = g
            by_subcohort[info["subcohort"]] += 1

        make_inputs._write_tsv(sf.tmp / "records.tsv", HROM1_RECORDS_HEAD, rec_rows)
        make_inputs._write_tsv(sf.tmp / "genetic_codes.tsv", c58.GC_HEAD, gc_rows)
        make_inputs._write_tsv(sf.tmp / "r232_status.tsv", HROM1_STATUS_HEAD, status_rows)
        make_inputs._write_tsv(sf.tmp / "left_out.tsv", HROM1_LEFT_HEAD, left_out)
        make_inputs._write_tsv(sf.tmp / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
        make_inputs._write_tsv(
            sf.tmp / "shared_records.tsv",
            ["ordinal", "genome", "record_ordinal", "contig_uid", "length", "owner_stage"],
            shared,
        )
        make_inputs._write_tsv(
            sf.tmp / "later_records.tsv",
            ["ordinal", "genome", "record_ordinal", "contig_uid", "length", "holder"],
            held_later,
        )
        (sf.tmp / "genomes.json").write_text(json.dumps(genomes))

        used = [k for k in HROM1_FASTA if any(ginfo[g]["category"] == k for g in kept)]
        source_meta = {
            "oral_glmc": {
                "source_name": "tower_hrom1_oral_glmc",
                "source_release": "Human Reference Oral Microbiome species reps (decodebiome.org HROM_Genome.All.tar.gz) and leafhopper symbiont/gut MAGs (zenodo:21341769 merged_MAGs.zip); genetic codes 11/4/25 screened by c86_o4_gv2_hrom1_inputs.py (DEC-80, DEC-88)",
                "license": "cc-by-4.0 as the job lists state them (sources.yaml tower_hrom_oral, tower_code4_glmc: no provider record kept)",
                "redistributable": False,
                "retrieved_at": retrieved_oral,
            },
            "actino1": {
                "source_name": "tower_code4_actino1",
                "source_release": "Code-4 Actinomycetota assemblies of BioProject PRJNA1366607 (Equivita/Tapirivita/Gorillivita, 35 INSDC assemblies); staged by c86_o4_gv2_hrom1_inputs.py (DEC-84 (2), DEC-88)",
                "license": "INSDC",
                "redistributable": True,
                "retrieved_at": retrieved_actino,
            },
        }
        _write_objects_json(sf, used, source_meta, Path(inputs_home))
        sha_map, byte_map = sf.commit(used)

        own_c, own_b = _tally_owners(contig_len, owners)
        lat_c, lat_b = _tally_owners(contig_len, later, owners)
        return {
            "stage": HROM1,
            "genomes_input": len(kept),
            "genomes": len(genomes),
            "repeated_genomes": len(repeated),
            "left_out": len(left_out),
            "left_out_by_reason": dict(sorted(Counter(x[4] for x in left_out).items())),
            "genomes_by_subcohort": dict(sorted(by_subcohort.items())),
            "genetic_codes": dict(sorted(Counter(str(resolved[g]) for g in kept).items())),
            "genetic_codes_by_subcohort": {
                sub: dict(sorted(Counter(str(resolved[g]) for g in kept if ginfo[g]["subcohort"] == sub).items()))
                for sub in sorted({ginfo[g]["subcohort"] for g in kept})
            },
            "code_conflicts": len(conflicts),
            "r232_status": dict(sorted(status_counts.items())),
            "contigs_total": sum(len(contigs[g][1]) for g in kept),
            "contigs_distinct": len(contig_len),
            "bases_distinct": sum(contig_len.values()),
            "held_by": own_c,
            "held_bases_by": own_b,
            "later_held_by": lat_c,
            "later_held_bases_by": lat_b,
            "contigs_new": sum(1 for cu in contig_len if cu not in owners),
            "bases_new": sum(n for cu, n in contig_len.items() if cu not in owners),
            "contigs_new_after_later": sum(1 for cu in contig_len if cu not in owners and cu not in later),
            "bases_new_after_later": sum(n for cu, n in contig_len.items() if cu not in owners and cu not in later),
            "sha256": sha_map,
            "bytes": byte_map,
            **{"listing_" + k: v for k, v in listing.report().items()},
        }
    except BaseException:
        sf.abort()
        raise
    finally:
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)


def hrom1_ani(
    meta: Path,
    fasta_paths: list[Path],
    out: Path,
    *,
    whole: dict[str, tuple[Path, list[Path]]],
    by_genus: dict[str, tuple[Path, list[Path]]] | None = None,
    taxonomy_tsv: Path | None = None,
    archaea: tuple[Path, Path] | None = None,
    threads: int = 8,
) -> dict:
    """Record skani >=96% ANI / >=50% larger AF matches against held stages without removing any genome (DEC-84 (2))."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    work = out / "work"
    shutil.rmtree(work, ignore_errors=True)
    q_files = c81.split_genomes(Path(meta), fasta_paths, work / "HROM1")
    (work / "q.txt").write_text("".join(f"{p}\n" for p in sorted(map(str, q_files.values()))))
    query_of = {str(p): g for g, p in q_files.items()}
    genera = {c81.genus_of(r["taxonomy"]) for r in c81._rows_tsv(Path(meta) / "r232_status.tsv")} - {""}

    groups: list[tuple[str, tuple[str, object]]] = [(n, ("whole", v)) for n, v in whole.items()]
    if by_genus:
        groups += [(n, ("genus", v)) for n, v in by_genus.items()]
    if archaea is not None:
        groups.append(("r226_archaea", ("archaea", archaea)))
    keep = c81.genus_genomes(taxonomy_tsv, genera) if (by_genus and taxonomy_tsv is not None) else set()

    res: dict = {
        "queries": len(q_files),
        "genera": len(genera),
        "threads": threads,
        "threshold": {"ani_min": c81.ANI_SAG, "af_min_larger_of_two": c81.AF_SAG},
        "groups": {},
    }
    hits_all: list[dict[str, dict]] = []
    for name, (kind, v) in groups:
        d = work / name
        unlisted: list[str] = []
        if kind == "archaea":
            files = c81.archaea_refs(v[0], v[1], genera, d)
        else:
            files = c81.split_genomes(v[0], v[1], d, keep=keep if kind == "genus" else None, unlisted=unlisted)
        ref_of = {str(p): (name, g) for g, p in files.items()}
        entry: dict = {
            "kind": kind,
            "references": len(files),
            "records_in_no_genome": len(unlisted),
        }
        if files:
            rlist = work / f"r_{name}.txt"
            rlist.write_text("".join(f"{p}\n" for p in sorted(ref_of)))
            table = out / f"hrom1_vs_{name}.tsv"
            t0 = c39_ani.child_cpu()
            subprocess.run(
                ["skani", "dist", "-t", str(threads), "--min-af", str(c81.AF_SAG),
                 "--ql", str(work / "q.txt"), "--rl", str(rlist), "-o", str(table)],
                check=True,
            )
            entry["cpu_s"] = round(c39_ani.child_cpu() - t0, 1)
            h = c81.best_hits(table, query_of, ref_of)
            entry["meet"] = len(h)
            hits_all.append(h)
        shutil.rmtree(d, ignore_errors=True)
        res["groups"][name] = entry

    best = c81.merge_best(hits_all)
    make_inputs._write_tsv(
        out / "hrom1_ani_hits.tsv",
        ["genome", "stage", "ref_genome", "ani", "af"],
        [(g, h["stage"], h["genome"], h["ani"], h["af"]) for g, h in sorted(best.items())],
    )
    res["meet"] = len(best)
    res["meet_by_stage"] = dict(sorted(Counter(h["stage"] for h in best.values()).items()))
    (out / "summary.json").write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    shutil.rmtree(work, ignore_errors=True)
    return res


# ==================================================================================================================
# 7. CLI
# ==================================================================================================================


def _add_common_stage_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--root", type=Path, default=None)
    p.add_argument("--root-stages", default=",".join(ROOT_STAGES))
    p.add_argument("--chain-home", type=Path, default=None, help="r26.11.4 stage inputs home (/tmp/if2453/s2/inputs)")
    p.add_argument("--with-r26-11-5-defaults", action="store_true",
                   help="include all 24 r26.11.5 stages (GRARE1..GBRMGD) in owner_inputs")
    p.add_argument("--owner", action="append", default=[], metavar="NAME=DIR")
    p.add_argument("--later", action="append", default=[], metavar="NAME=DIR")
    p.add_argument("--inputs-home", type=Path, default=Path("/tmp/if2453/s2/inputs"))
    p.add_argument("--work", type=Path, default=None)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--memory-limit", default="4GB")
    p.add_argument("--report-out", type=Path, default=None)


def _resolve_owners_and_later(a) -> tuple[list[str], dict[str, Path], dict[str, Path]]:
    root_stages = [s for s in a.root_stages.split(",") if s]
    owners: dict[str, Path] = {}
    if a.chain_home is not None:
        for s in R26_11_4_STAGES:
            owners[s] = Path(a.chain_home) / s
    if a.with_r26_11_5_defaults:
        for s in R26_11_5_STAGES:
            owners[s] = Path(DEFAULT_R26_11_5_DIRS[s])
    for k, v in c81._pairs(a.owner, "--owner").items():
        owners[k] = Path(v)
    later = {k: Path(v) for k, v in c81._pairs(a.later, "--later").items()}
    return root_stages, owners, later


def _parse_group_sets(items: list[str], what: str) -> dict[str, tuple[Path, list[Path]]]:
    out: dict[str, tuple[Path, list[Path]]] = {}
    for k, v in c81._pairs(items, what).items():
        m, sep, fas = v.partition(":")
        if not sep or not m or not fas:
            raise Refused(f"{what} {k}={v}: expected NAME=META_DIR:FASTA[,FASTA]")
        out[k] = (Path(m), [Path(x) for x in fas.split(",") if x])
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Card C-86: stage O4, GV2, and HROM1 cohorts")
    sub = ap.add_subparsers(dest="cmd", required=True)

    pc = sub.add_parser("census")
    pc.add_argument("paths", nargs="+", type=Path)
    pc.add_argument("--out", type=Path, default=None)

    po4 = sub.add_parser("o4-stage")
    po4.add_argument("--o4-org-dir", type=Path, required=True)
    po4.add_argument("--fungal-fa", type=Path, required=True)
    po4.add_argument("--fungal-tsv", type=Path, required=True)
    po4.add_argument("--meta-out", type=Path, required=True)
    po4.add_argument("--fasta-out", type=Path, required=True)
    _add_common_stage_args(po4)

    pgv = sub.add_parser("gv2-stage")
    pgv.add_argument("--dir", action="append", default=[], metavar="SUB=DIR",
                     help="may be repeated per subcohort (e.g. both /stor/.../dl/SUB and /tmp/.../dl/SUB)")
    pgv.add_argument("--meta-out", type=Path, required=True)
    pgv.add_argument("--fasta-out", type=Path, required=True)
    pgv.add_argument("--exclude-gks-all-mags", action="store_true")
    pgv.add_argument("--screen-code25", action="store_true")
    _add_common_stage_args(pgv)

    phr = sub.add_parser("hrom1-stage")
    phr.add_argument("--hrom-tar", type=Path, required=True)
    phr.add_argument("--hrom-meta", type=Path, required=True)
    phr.add_argument("--glmc-zip", type=Path, required=True)
    phr.add_argument("--glmc-csv", type=Path, required=True)
    phr.add_argument("--actino-dir", type=Path, required=True)
    phr.add_argument("--index-dir", type=Path, required=True)
    phr.add_argument("--meta-out", type=Path, required=True)
    phr.add_argument("--fasta-out", type=Path, required=True)
    phr.add_argument("--hrom-hq-only", action="store_true")
    _add_common_stage_args(phr)

    pha = sub.add_parser("hrom1-ani")
    pha.add_argument("meta", type=Path)
    pha.add_argument("out", type=Path)
    pha.add_argument("--fasta", action="append", required=True, type=Path)
    pha.add_argument("--whole", action="append", default=[], metavar="NAME=META_DIR:FASTA[,FASTA]")
    pha.add_argument("--by-genus", action="append", default=[], metavar="NAME=META_DIR:FASTA[,FASTA]")
    pha.add_argument("--taxonomy", type=Path, default=None)
    pha.add_argument("--archaea", default="", metavar="FASTA:LIST_TSV")
    pha.add_argument("--threads", type=int, default=8)

    a = ap.parse_args(argv)
    try:
        if a.cmd == "census":
            rep = census(a.paths)
            txt = json.dumps(rep, indent=1, sort_keys=True)
            if a.out is not None:
                a.out.parent.mkdir(parents=True, exist_ok=True)
                a.out.write_text(txt + "\n")
            print(txt)
            return 0

        if a.cmd == "hrom1-ani":
            arc = None
            if a.archaea:
                fa, sep, lst = a.archaea.partition(":")
                if not sep or not fa or not lst:
                    raise Refused(f"--archaea {a.archaea}: expected FASTA:LIST_TSV")
                arc = (Path(fa), Path(lst))
            rep = hrom1_ani(
                a.meta, a.fasta, a.out,
                whole=_parse_group_sets(a.whole, "--whole"),
                by_genus=_parse_group_sets(a.by_genus, "--by-genus"),
                taxonomy_tsv=a.taxonomy,
                archaea=arc,
                threads=a.threads,
            )
            print(json.dumps(rep, indent=1, sort_keys=True))
            return 0

        root_stages, owners, later = _resolve_owners_and_later(a)
        if a.cmd == "o4-stage":
            stage_name = O4
            rep = stage_o4(
                a.o4_org_dir, a.fungal_fa, a.fungal_tsv, a.meta_out, a.fasta_out,
                root=a.root, root_stages=root_stages, owner_inputs=owners, later_inputs=later,
                inputs_home=a.inputs_home, work=a.work, threads=a.threads,
                memory_limit=a.memory_limit, workers=a.workers,
            )
        elif a.cmd == "gv2-stage":
            stage_name = GV2
            gv_dirs: dict[str, list[Path]] = {}
            for item in a.dir:
                k, sep, v = item.partition("=")
                if not sep or not k or not v:
                    raise Refused(f"--dir {item!r}: expected SUB=DIR")
                gv_dirs.setdefault(k, []).append(Path(v))
            rep = stage_gv2(
                gv_dirs, a.meta_out, a.fasta_out,
                root=a.root, root_stages=root_stages, owner_inputs=owners, later_inputs=later,
                inputs_home=a.inputs_home, work=a.work, threads=a.threads,
                memory_limit=a.memory_limit, exclude_gks_all_mags=a.exclude_gks_all_mags,
                screen_code25=a.screen_code25,
            )
        else:
            stage_name = HROM1
            rep = stage_hrom1(
                a.hrom_tar, a.hrom_meta, a.glmc_zip, a.glmc_csv, a.actino_dir, a.index_dir,
                a.meta_out, a.fasta_out,
                root=a.root, root_stages=root_stages, owner_inputs=owners, later_inputs=later,
                inputs_home=a.inputs_home, hrom_hq_only=a.hrom_hq_only,
                work=a.work, threads=a.threads, memory_limit=a.memory_limit, workers=a.workers,
            )
    except Refused as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2

    rep["owner_stages"] = [*(root_stages if a.root else []), *owners]
    rep["later_stages"] = [*later]
    rep["s2_check"] = c81._check(stage_name, a.meta_out / stage_name)
    txt = json.dumps(rep, indent=1, sort_keys=True)
    if a.report_out is not None:
        a.report_out.parent.mkdir(parents=True, exist_ok=True)
        a.report_out.write_text(txt + "\n")
    print(txt)
    return 0 if rep["s2_check"]["exit"] == 0 and rep["s2_check"]["problems"] == [] else 1


if __name__ == "__main__":
    raise SystemExit(main())
