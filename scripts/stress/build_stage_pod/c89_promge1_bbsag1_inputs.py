"""Card C-89 (DEC-65, DEC-80, DEC-84 (2), DEC-88): stage tower-landed MGE/plasmid (`PROMGE1`) and microfluidic
single-cell (`BBSAG1`, or `BBSAG1A` + `BBSAG1B` if net-new bases exceed 40 Gbp) cohorts on `wilkcomp02` for `r26.11.6`.

Subcommands:
    python c89_promge1_bbsag1_inputs.py census PATH [PATH ...] [--out JSON]
    python c89_promge1_bbsag1_inputs.py promge1-stage --dir SUB=DIR ... --index-dir DIR
        --meta-out DIR --fasta-out DIR [common options]
    python c89_promge1_bbsag1_inputs.py bbsag1-stage --bbsag20-dir DIR ... --haima-dir DIR ...
        --deathvalley-dir DIR ... --index-dir DIR --meta-out DIR --fasta-out DIR [common options]
    python c89_promge1_bbsag1_inputs.py bbsag1-ani META_DIR OUT_DIR --fasta FASTA ...
        [--whole NAME=META:FASTA[,FASTA]] [--by-genus NAME=META:FASTA[,FASTA]]
        [--taxonomy TSV] [--archaea FASTA:LIST] [--threads 4]
"""

from __future__ import annotations

import argparse
import bz2
import contextlib
import csv
import gzip
import io
import json
import lzma
import os
import re
import shutil
import subprocess
import sys
import tarfile
import warnings
import zipfile
from collections import Counter
from collections.abc import Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import a26_r26_11_4_chain as a26  # noqa: E402
import c39_ani  # noqa: E402
import c48c_virus_inputs as v2  # noqa: E402
import c58_sags_symbionts_niche as c58  # noqa: E402
import c62_metavr_inputs as c62  # noqa: E402
import c77_sag2_inputs as c77  # noqa: E402
import c78_compact_tower_sets as c78  # noqa: E402
import c81_gwv1_sag3_wet1 as c81  # noqa: E402
import make_inputs  # noqa: E402
import r232_index  # noqa: E402
from biobox import fasta, ids  # noqa: E402
from biobox.l1 import genecall  # noqa: E402

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
C86_STAGES: tuple[str, ...] = ("O4", "GV2", "HROM1")

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
DEFAULT_C86_DIRS: dict[str, str] = {
    "O4": "/tmp/if2453/c86/meta/O4",
    "GV2": "/tmp/if2453/c86/meta/GV2",
    "HROM1": "/tmp/if2453/c86/meta/HROM1",
}

STOR_MIN_BYTES = 64 * 1024 * 1024  # AGENTS.md §7.3: >= 67,108,864 B on /stor
TASKS_PER_CHILD = 500  # recycle worker processes well below 280 CPU-min guard
PARQUET_FLUSH = 250_000
BLOCK_BYTES = 4 << 20

IUPAC_BYTES = b"ACGTURYKMSWBDHVNacgturykmswbdhvn"
VALID_DNA_RE = re.compile(rb"^[ACGTURYKMSWBDHVNacgturykmswbdhvn]+$")
GAP_RE = re.compile(rb"[-.]")
PROTEIN_ONLY = set(b"EFILPQZJO*efilpqzjo")
CIRCULAR_RE = re.compile(r"circular", re.IGNORECASE)
PRODIGAL_FFN_HDR_RE = re.compile(r"^[^\s#]+\s+#\s+\d+\s+#\s+\d+\s+#\s+-?1\s+#\s+ID=")
ACC_GCA_RE = re.compile(r"GC[AF]_\d{9}(?:\.\d+)?")

FASTA_EXTS = (
    ".fa", ".fna", ".fasta", ".fas", ".ffn", ".faa", ".prt", ".fsa",
    ".fa.gz", ".fna.gz", ".fasta.gz", ".ffn.gz", ".faa.gz", ".prt.gz",
)
DNA_FASTA_EXTS = (".fa", ".fna", ".fasta", ".fas", ".fsa", ".fa.gz", ".fna.gz", ".fasta.gz")

_CODE25_TAX_RE = re.compile(
    r"(patescibacteria|patescibacteriota|saccharibacteria|saccharimonadia|gracilibacteria|"
    r"absconditabacteria|absconditabacterales|bd1-5|sr1|jaedam01|gn02|tm7|candidatus\s+nanosynbacter)",
    re.IGNORECASE,
)
_CODE4_TAX_RE = re.compile(
    r"(hodgkinia|nasuia|zinderia|mycoplasma|mycoplasmatota|ureaplasma|spiroplasma|mesoplasma|"
    r"entomoplasma|mycoplasmoides|mycoplasmopsis|metamycoplasma|malacoplasma|acholeplasma|"
    r"phytoplasma|equivita|tapirivita|gorillivita|altericodex|stammera)",
    re.IGNORECASE,
)


# ==================================================================================================================
# 1. Fast block-streaming FASTA reader, magic-byte detection, and archive census
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


def _is_macosx_meta(name: str) -> bool:
    norm = name.replace("\\", "/")
    base = os.path.basename(norm.rstrip("/"))
    return norm.startswith("__MACOSX/") or "/__MACOSX/" in norm or base == ".DS_Store" or base.startswith("._")


def _looks_fasta(name: str) -> bool:
    base = os.path.basename(name.lower())
    return not _is_macosx_meta(name) and not base.startswith(".") and base.endswith(FASTA_EXTS)


def _looks_dna_fasta(name: str) -> bool:
    low = name.lower()
    base = os.path.basename(low)
    if _is_macosx_meta(name) or base.startswith("."):
        return False
    if base.endswith((".proteins.fa", ".proteins.fasta", ".faa", ".faa.gz", ".ffn", ".ffn.gz", ".prt", ".prt.gz")):
        return False
    return base.endswith(DNA_FASTA_EXTS)


def iter_fasta_fast(fobj) -> Iterator[tuple[str, str, bytes, int]]:
    """Fast block-streaming FASTA reader yielding `(accession, header, seq_bytes, first_line_width)`."""
    carry = b""
    hdr: str | None = None
    acc: str = ""
    seq_parts: list[bytes] = []
    width: int = 0

    def _emit():
        nonlocal hdr, acc, seq_parts, width
        if hdr is not None:
            seq = b"".join(seq_parts)
            res = (acc, hdr, seq, width or 70)
            hdr = None
            acc = ""
            seq_parts = []
            width = 0
            return res
        return None

    while True:
        chunk = fobj.read(BLOCK_BYTES)
        if not chunk:
            break
        buf = carry + chunk if carry else chunk
        cut = buf.rfind(b"\n") + 1
        if cut == 0:
            carry = buf
            continue
        carry = buf[cut:]
        block = buf[:cut]
        i = 0
        n = len(block)
        while i < n:
            if block[i] == 0x3E:  # '>'
                prev = _emit()
                if prev is not None:
                    yield prev
                nl = block.find(b"\n", i)
                raw_h = block[i + 1 : nl].rstrip(b"\r")
                hdr = raw_h.decode("utf-8")
                acc = hdr.split(None, 1)[0] if hdr.strip() else ""
                i = nl + 1
            else:
                if hdr is None:
                    nl = block.find(b"\n", i)
                    if block[i:nl].strip():
                        raise ValueError("sequence bytes before first FASTA header")
                    i = nl + 1
                    continue
                if width == 0:
                    nl = block.find(b"\n", i)
                    ln = block[i:nl].rstrip(b"\r")
                    if ln:
                        width = len(ln)
                        seq_parts.append(ln)
                    i = nl + 1
                    continue
                nxt = block.find(b"\n>", i)
                if nxt == -1:
                    seq_parts.append(block[i:].translate(None, b"\r\n"))
                    break
                seq_parts.append(block[i : nxt + 1].translate(None, b"\r\n"))
                i = nxt + 1

    if carry:
        for raw_ln in carry.splitlines():
            ln = raw_ln.rstrip(b"\r")
            if not ln:
                continue
            if ln.startswith(b">"):
                prev = _emit()
                if prev is not None:
                    yield prev
                hdr = ln[1:].decode("utf-8")
                acc = hdr.split(None, 1)[0] if hdr.strip() else ""
            else:
                if hdr is None:
                    raise ValueError("sequence bytes before first FASTA header")
                if width == 0:
                    width = len(ln)
                seq_parts.append(ln)
    last = _emit()
    if last is not None:
        yield last


def _seq_alphabet(seq: bytes) -> str:
    if not seq:
        return "empty"
    if not seq.translate(None, IUPAC_BYTES):
        return "dna"
    if GAP_RE.search(seq):
        return "alignment"
    if PROTEIN_ONLY & set(seq):
        return "protein"
    return "non_iupac"


def _fasta_member_stats(fobj, name: str, size: int, max_headers: int = 3) -> dict:
    n = bp = circ = non_iupac_recs = ffn_headers = 0
    mn, mx = None, 0
    heads: list[str] = []
    kinds: set[str] = set()
    try:
        for _, hdr, seq, _ in iter_fasta_fast(fobj):
            n += 1
            bp += len(seq)
            mn = len(seq) if mn is None else min(mn, len(seq))
            mx = max(mx, len(seq))
            circ += bool(CIRCULAR_RE.search(hdr))
            if PRODIGAL_FFN_HDR_RE.match(hdr):
                ffn_headers += 1
            if len(heads) < max_headers:
                heads.append(hdr[:160])
            k = _seq_alphabet(seq)
            kinds.add(k)
            if k != "dna":
                non_iupac_recs += 1
    except (ValueError, UnicodeDecodeError, OSError, EOFError, lzma.LZMAError):
        return {"name": name, "bytes": size, "alphabet": "not_fasta"}
    if n == 0:
        return {"name": name, "bytes": size, "alphabet": "not_fasta"}
    if kinds == {"dna"}:
        alpha = "prodigal_cds_ffn" if ffn_headers == n else "dna"
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
        "prodigal_ffn_headers": ffn_headers,
        "non_dna_records": non_iupac_recs,
        "first_headers": heads,
    }


def _open_maybe_compressed(raw_fobj, name: str):
    low = name.lower()
    if low.endswith(".gz"):
        return gzip.open(raw_fobj, "rb")
    if low.endswith(".bz2"):
        return bz2.open(raw_fobj, "rb")
    if low.endswith(".xz"):
        return lzma.open(raw_fobj, "rb")
    return raw_fobj


def census(paths: Iterable[Path | str]) -> dict:
    out: dict[str, dict] = {}
    for p in map(Path, paths):
        mg = magic_format(p)
        low = p.name.lower()
        members: list[dict] = []
        if mg == "zip":
            with zipfile.ZipFile(p) as zf:
                for info in zf.infolist():
                    if info.is_dir():
                        members.append({"name": info.filename, "bytes": 0, "alphabet": "directory"})
                    elif _is_macosx_meta(info.filename):
                        members.append({"name": info.filename, "bytes": info.file_size, "alphabet": "macosx"})
                    elif _looks_fasta(info.filename):
                        with zf.open(info) as raw:
                            with _open_maybe_compressed(raw, info.filename) as stream:
                                members.append(_fasta_member_stats(stream, info.filename, info.file_size))
                    else:
                        members.append({"name": info.filename, "bytes": info.file_size, "alphabet": "not_fasta"})
        elif mg == "tar" or (mg in ("gzip", "bzip2", "xz") and (".tar." in low or low.endswith((".tgz", ".tbz2", ".txz")))):
            with tarfile.open(p, "r|*") as tf:
                for m in tf:
                    if m.isdir():
                        members.append({"name": m.name.rstrip("/") + "/", "bytes": 0, "alphabet": "directory"})
                    elif _is_macosx_meta(m.name):
                        members.append({"name": m.name, "bytes": m.size, "alphabet": "macosx"})
                    elif m.isfile() and _looks_fasta(m.name):
                        ef = tf.extractfile(m)
                        assert ef is not None
                        with _open_maybe_compressed(ef, m.name) as stream:
                            members.append(_fasta_member_stats(stream, m.name, m.size))
                    else:
                        members.append({"name": m.name, "bytes": m.size, "alphabet": "not_fasta"})
        else:
            with open(p, "rb") as raw:
                if mg == "gzip":
                    stream = gzip.open(raw, "rb")
                elif mg == "bzip2":
                    stream = bz2.open(raw, "rb")
                elif mg == "xz":
                    stream = lzma.open(raw, "rb")
                else:
                    stream = raw
                with stream:
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
# 2. Prokaryotic genetic-code screening (11, 4, 25)
# ==================================================================================================================


def _train_call_stats(seq_strs: list[str], total_bp: int, table: int) -> dict:
    import pyrodigal

    if total_bp <= 0:
        return {"cd": 0.0, "mean_orf": 0.0, "n_genes": 0, "n_cod": 0, "n_tga": 0, "n_tgg": 0}
    finder = pyrodigal.GeneFinder(meta=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        if total_bp < 20_000:
            joined = "".join(seq_strs)
            reps = (20_000 // max(1, len(joined))) + 2
            finder.train(joined * reps, translation_table=table)
        else:
            finder.train(*seq_strs, translation_table=table)

    coding_bp = n_genes = n_cod = n_tga = n_tgg = 0
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
    return {
        "cd": coding_bp / total_bp,
        "mean_orf": (coding_bp / n_genes) if n_genes > 0 else 0.0,
        "n_genes": n_genes,
        "n_cod": n_cod,
        "n_tga": n_tga,
        "n_tgg": n_tgg,
    }


def screen_prokaryote_code(
    seqs: list[bytes],
    taxonomy_hint: str = "",
    *,
    force_alt: bool = False,
    screen_max_bp: int = 500_000,
) -> dict:
    """Screen prokaryotic SAGs/MAGs across genetic codes 11, 4 (UGA->Trp), and 25 (UGA->Gly)."""
    pre_strs, pre_bp, is_full = c58._prefix_seqs_for_screen(seqs, screen_max_bp)
    if pre_bp < 5_000:
        code = 4 if _CODE4_TAX_RE.search(taxonomy_hint) else (25 if _CODE25_TAX_RE.search(taxonomy_hint) else 11)
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
    tax_c25 = bool(_CODE25_TAX_RE.search(taxonomy_hint))
    tax_c4 = bool(_CODE4_TAX_RE.search(taxonomy_hint))
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


# ==================================================================================================================
# 3. Bounded-memory DuckDB ownership & shared file helpers
# ==================================================================================================================


def _sql_str(s: object) -> str:
    return "'" + str(s).replace("'", "''") + "'"


class ContigPartWriter:
    """Spools `(contig_uid, length)` pairs into bounded Parquet chunks on disk so DuckDB computes ownership in O(1) RAM."""

    def __init__(self, out_dir: Path) -> None:
        import pyarrow as pa

        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.schema = pa.schema([("contig_uid", pa.binary(16)), ("length", pa.int64())])
        self.parts: list[str] = []
        self._cu: list[bytes] = []
        self._ln: list[int] = []

    def add(self, cu: bytes, length: int) -> None:
        self._cu.append(cu)
        self._ln.append(length)
        if len(self._cu) >= PARQUET_FLUSH:
            self.flush()

    def flush(self) -> None:
        import pyarrow as pa
        import pyarrow.parquet as pq

        if not self._cu and self.parts:
            return
        p = self.out_dir / f"part-{len(self.parts):05d}.parquet"
        pq.write_table(pa.table({"contig_uid": self._cu, "length": self._ln}, schema=self.schema), p)
        self.parts.append(str(p))
        self._cu.clear()
        self._ln.clear()


def _duck_ownership_from_parts(
    part_files: list[str],
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    later_inputs: dict[str, Path],
    work: Path,
    threads: int,
    memory_limit: str,
) -> dict:
    """Compute exact build-order `contig_uid` ownership in DuckDB with spilling enabled (`AGENTS.md` §4.1)."""
    import duckdb
    from biobox import io as bio_io
    from biobox.l0seq.store import INDEX_DIR

    for name, d in owner_inputs.items():
        if not (Path(d) / "genomes.json").is_file():
            raise Refused(f"--owner {name}: no genomes.json under {d}")
    for name, d in later_inputs.items():
        if not (Path(d) / "genomes.json").is_file():
            raise Refused(f"--later {name}: no genomes.json under {d}")

    duck_tmp = Path(work) / "duck_own"
    shutil.rmtree(duck_tmp, ignore_errors=True)
    duck_tmp.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(
        config={
            "threads": max(1, int(threads)),
            "memory_limit": memory_limit,
            "temp_directory": str(duck_tmp),
            "preserve_insertion_order": False,
        }
    )
    try:
        con.execute(
            f"CREATE TEMP TABLE stage_u AS SELECT contig_uid, max(length)::BIGINT AS length "
            f"FROM read_parquet({a26._lit(part_files)}) GROUP BY contig_uid"
        )
        n_dist, bp_dist = con.execute(
            "SELECT count(*)::BIGINT, coalesce(sum(length), 0)::BIGINT FROM stage_u"
        ).fetchone()

        con.execute("CREATE TEMP TABLE matched_owners (contig_uid BLOB, owner VARCHAR)")
        con.execute("CREATE TEMP TABLE matched_later (contig_uid BLOB, holder VARCHAR)")

        def _load_genomes_into_cur(gj: Path) -> bool:
            txt_head = ""
            with open(gj, "r", encoding="utf-8") as fh:
                txt_head = fh.read(64).strip()
            if txt_head == "[]":
                return False
            con.execute(
                "CREATE OR REPLACE TEMP TABLE cur_g AS SELECT DISTINCT contig_uid FROM "
                f"(SELECT from_hex(unnest(contig_uids)) AS contig_uid FROM read_json({_sql_str(gj)}, format='array', "
                "columns={'contig_uids': 'VARCHAR[]'}, maximum_object_size=1073741824)) "
                "SEMI JOIN stage_u USING (contig_uid)"
            )
            return True

        roots = list(root_stages)
        if root is not None and roots:
            base = Path(root) / INDEX_DIR
            for st in roots:
                sdir = base / f"stage={st}"
                if not sdir.is_dir():
                    raise Refused(f"root {root} has no contig_index for stage {st}")
                files = bio_io.part_files(sdir)
                if not files:
                    continue
                con.execute(
                    f"INSERT INTO matched_owners SELECT DISTINCT p.contig_uid, {_sql_str(st)} "
                    f"FROM read_parquet({a26._lit(files)}, hive_partitioning=false) p "
                    "JOIN stage_u u USING (contig_uid) "
                    "ANTI JOIN matched_owners m USING (contig_uid) "
                    "WHERE p.is_default_variant"
                )

        for st, d in owner_inputs.items():
            if _load_genomes_into_cur(Path(d) / "genomes.json"):
                con.execute(
                    f"INSERT INTO matched_owners SELECT contig_uid, {_sql_str(st)} FROM cur_g "
                    "ANTI JOIN matched_owners USING (contig_uid)"
                )

        for st, d in later_inputs.items():
            if _load_genomes_into_cur(Path(d) / "genomes.json"):
                con.execute(
                    f"INSERT INTO matched_later SELECT contig_uid, {_sql_str(st)} FROM cur_g "
                    "ANTI JOIN matched_owners USING (contig_uid) "
                    "ANTI JOIN matched_later USING (contig_uid)"
                )

        owners_map = {bytes(cu): str(own) for cu, own in con.execute("SELECT contig_uid, owner FROM matched_owners").fetchall()}
        later_map = {bytes(cu): str(h) for cu, h in con.execute("SELECT contig_uid, holder FROM matched_later").fetchall()}

        own_rows = con.execute(
            "SELECT m.owner, count(*)::BIGINT, coalesce(sum(u.length), 0)::BIGINT "
            "FROM matched_owners m JOIN stage_u u USING (contig_uid) GROUP BY m.owner ORDER BY m.owner"
        ).fetchall()
        lat_rows = con.execute(
            "SELECT m.holder, count(*)::BIGINT, coalesce(sum(u.length), 0)::BIGINT "
            "FROM matched_later m JOIN stage_u u USING (contig_uid) GROUP BY m.holder ORDER BY m.holder"
        ).fetchall()

        held_by = {r[0]: int(r[1]) for r in own_rows}
        held_bases_by = {r[0]: int(r[2]) for r in own_rows}
        later_by = {r[0]: int(r[1]) for r in lat_rows}
        later_bases_by = {r[0]: int(r[2]) for r in lat_rows}
        hc = sum(held_by.values())
        hb = sum(held_bases_by.values())
        lc = sum(later_by.values())
        lb = sum(later_bases_by.values())
        return {
            "owners": owners_map,
            "later": later_map,
            "contigs_distinct": int(n_dist),
            "bases_distinct": int(bp_dist),
            "held_by": held_by,
            "held_bases_by": held_bases_by,
            "later_held_by": later_by,
            "later_held_bases_by": later_bases_by,
            "contigs_new": int(n_dist) - hc,
            "bases_new": int(bp_dist) - hb,
            "contigs_new_after_later": int(n_dist) - hc - lc,
            "bases_new_after_later": int(bp_dist) - hb - lb,
        }
    finally:
        con.close()
        shutil.rmtree(duck_tmp, ignore_errors=True)


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


def _commit_with_stor_floor(sf: c78.StageFiles, used_categories: list[str]) -> tuple[dict[str, str], dict[str, int], dict[str, str]]:
    """Commit `StageFiles` and move `records.tsv` to `fasta_dir` if `records.tsv >= 64 MiB` (`AGENTS.md` §7.3)."""
    big_meta_moved: dict[str, Path] = {}
    rec_tmp = sf.tmp / "records.tsv"
    if rec_tmp.is_file() and rec_tmp.stat().st_size >= STOR_MIN_BYTES:
        sf.fasta_dir.mkdir(parents=True, exist_ok=True)
        part = sf.fasta_dir / ".records.tsv.partial"
        dst = sf.fasta_dir / "records.tsv"
        shutil.copyfile(rec_tmp, part)
        os.replace(part, dst)
        rec_tmp.unlink()
        big_meta_moved["records.tsv"] = dst

    sha_map, byte_map = sf.commit(used_categories)
    locations: dict[str, str] = {}
    for k in used_categories:
        fn = sf.fasta_names[k]
        p = sf.fasta_dir / fn if (sf.fasta_dir / fn).is_file() else sf.meta_dir / fn
        locations[fn] = str(p)
    for fn, dst in big_meta_moved.items():
        sha_map[fn] = ids.sha256_file(dst).hex()
        byte_map[fn] = dst.stat().st_size
        locations[fn] = str(dst)
    for p in sorted(sf.meta_dir.iterdir()):
        if p.is_file():
            locations.setdefault(p.name, str(p))
    return sha_map, byte_map, locations


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
            if p.is_file() and p.name not in ("SHA256SUMS", "SOURCE_URLS.tsv") and not p.name.startswith("."):
                seen.setdefault(p.name, p)
    return [seen[k] for k in sorted(seen)]


def _stem_clean(name: str) -> str:
    base = os.path.basename(name)
    changed = True
    while changed:
        changed = False
        for ext in (".gz", ".bz2", ".xz", ".fasta", ".fna", ".fa", ".fas", ".fsa", "_genomic"):
            if base.lower().endswith(ext) and len(base) > len(ext):
                base = base[: -len(ext)]
                changed = True
    return base


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


# ==================================================================================================================
# 4. Stage PROMGE1: 8 MGE / plasmid / Borg / satellite cohorts (`kind="plasmid"`)
# ==================================================================================================================

PROMGE1 = "PROMGE1"
PROMGE1_SUBCOHORTS = (
    "PROMGE2",
    "MFD_PLAS",
    "GSPR",
    "ECE_HODOR",
    "SAT_P4P2",
    "ECE_METH",
    "PICJI_SAT",
    "BORG_ONT",
)
PROMGE1_FASTA = {
    "promge2": "promge1_promge2.fna.gz",
    "open": "promge1_plasmids_eces.fna.gz",
}
PROMGE1_RETRIEVED_PROMGE2 = "2026-11-27T08:28:42Z"
PROMGE1_RETRIEVED_OPEN = "2026-11-27T02:30:35Z"

PROMGE1_RECORDS_HEAD = [
    "ordinal", "subcohort", "genome", "record_ordinal", "contig_uid", "length",
    "circular", "dtr_class", "trimmed_bases", "mge_class", "owner_stage", "later_holder",
]
PROMGE1_LEFT_HEAD = ["subcohort", "member", "accession", "length", "reason"]
PROMGE1_SKIP_HEAD = ["subcohort", "member", "bytes", "reason"]
PROMGE1_STATUS_HEAD = [
    "accession", "subcohort", "role", "r232_status", "in_gtdb_r232",
    "is_gtdb_representative", "gtdb_representative", "gtdb_species",
    "in_globdb_r232", "globdb_r232_id", "records",
]


def _parse_promge2_header(hdr: str) -> tuple[str, dict[str, str], str]:
    tok, _, rest = hdr.partition(" ")
    tags: dict[str, str] = {}
    for kv in rest.split(";"):
        k, sep, v = kv.partition("=")
        if sep:
            tags[k.strip()] = v.strip()
    m_acc = ACC_GCA_RE.search(tok)
    return tok, tags, (m_acc.group(0) if m_acc else "")


def _is_promge2_cellular_only(mge_val: str) -> bool:
    if not mge_val:
        return False
    parts = [p.split(":", 1)[0].strip().lower() for p in mge_val.split(",") if p.strip()]
    return bool(parts) and all(p == "cellular" for p in parts)


_PROMGE2_HOST_ACC_RE = re.compile(r"^MGE_(GC[AF]_\d{9}(?:\.\d+)?)_")


@contextlib.contextmanager
def _open_gzip_fast(path: Path):
    """Stream gzip bytes via `pigz -dc` / `gzip -dc` when available, falling back cleanly to `gzip.open`."""
    gz_bin = shutil.which("pigz") or shutil.which("gzip")
    if gz_bin is None:
        with gzip.open(path, "rb") as gfh:
            yield gfh
        return
    cmd = [gz_bin, "-dc", "-p", "2", str(path)] if os.path.basename(gz_bin) == "pigz" else [gz_bin, "-dc", str(path)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    normal_exit = False
    try:
        assert proc.stdout is not None
        yield proc.stdout
        normal_exit = True
    finally:
        if not normal_exit:
            proc.terminate()
            proc.wait()
        else:
            _, err = proc.communicate()
            if proc.returncode != 0:
                raise Refused(
                    f"{path}: gzip decompression failed (exit {proc.returncode}): {err.decode('utf-8', 'replace')[:200]}"
                )


def _is_promge2_qualifying_member(member_id: str, pg3_rep_accs: set[str]) -> bool:
    """Option 3b: keep all SPIRE metagenomic MAG hosts (non-GCA/GCF) and proGenomes3 `specI` species reps."""
    m = _PROMGE2_HOST_ACC_RE.match(member_id)
    if m is not None:
        return m.group(1) in pg3_rep_accs
    return True


def read_promge2_qualifying_hosts(promge2_dirs: list[Path]) -> tuple[set[str], set[str]]:
    """Return `(pg3_rep_accs, rescued_rep_ids)` from `genomes_taxonomy.tsv.gz` (or `genomes_quality.tsv.gz`)
    and `mge100_clusters.tsv.gz`."""
    tax_p = _find_in_dirs(promge2_dirs, "genomes_taxonomy.tsv.gz")
    qual_p = _find_in_dirs(promge2_dirs, "genomes_quality.tsv.gz")
    if tax_p is None and qual_p is None:
        raise Refused("PROMGE2: --promge2-spire-and-pg3-reps-only requires genomes_taxonomy.tsv.gz or genomes_quality.tsv.gz")

    pg3_rep_accs: set[str] = set()
    if tax_p is not None:
        with gzip.open(tax_p, "rt", encoding="utf-8", errors="replace", newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            if "genome" not in (rd.fieldnames or []) or "representative" not in (rd.fieldnames or []):
                raise Refused(f"{tax_p}: missing genome or representative column")
            for r in rd:
                gid = (r.get("genome") or "").strip()
                if gid.startswith(("GCA_", "GCF_")) and (r.get("representative") or "").strip().lower() == "true":
                    pg3_rep_accs.add(gid)
    else:
        assert qual_p is not None
        with gzip.open(qual_p, "rt", encoding="utf-8", errors="replace", newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            if "genome_id" not in (rd.fieldnames or []) or "representative" not in (rd.fieldnames or []):
                raise Refused(f"{qual_p}: missing genome_id or representative column")
            for r in rd:
                gid = (r.get("genome_id") or "").strip()
                if gid.startswith(("GCA_", "GCF_")) and (r.get("representative") or "").strip().lower() == "true":
                    pg3_rep_accs.add(gid)

    rescued_rep_ids: set[str] = set()
    cl_p = _find_in_dirs(promge2_dirs, "mge100_clusters.tsv.gz")
    if cl_p is not None:
        with _open_gzip_fast(cl_p) as raw, io.TextIOWrapper(raw, encoding="utf-8", errors="replace") as fh:
            hdr = fh.readline()
            if "representative_id" not in hdr or "member_id" not in hdr:
                raise Refused(f"{cl_p}: missing member_id or representative_id column")
            for ln in fh:
                cols = ln.rstrip("\n").split("\t")
                if len(cols) < 3:
                    continue
                mem, rep_id = cols[1], cols[2]
                if rep_id in rescued_rep_ids or _is_promge2_qualifying_member(rep_id, pg3_rep_accs):
                    continue
                if _is_promge2_qualifying_member(mem, pg3_rep_accs):
                    rescued_rep_ids.add(rep_id)
    return pg3_rep_accs, rescued_rep_ids


def read_mfd_plasmid_metadata(path: Path) -> dict[str, str]:
    """Return plasmid_id -> topology string from MFD `pfd_metadata.tsv`."""
    out: dict[str, str] = {}
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        fields = rd.fieldnames or []
        if "topology" not in fields:
            raise Refused(f"{path}: missing topology column")
        for r in rd:
            topo = (r.get("topology") or "").strip()
            for col in ("member_id", "representative"):
                pid = (r.get(col) or "").strip()
                if pid:
                    out.setdefault(pid, topo)
    return out


def read_ncbi_batch6b_manifest(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        fields = rd.fieldnames or []
        col = "accession" if "accession" in fields else ("assembly_accession" if "assembly_accession" in fields else None)
        if col is None:
            raise Refused(f"{path}: missing accession column")
        for r in rd:
            acc = (r.get(col) or "").strip()
            if acc:
                out[acc] = {k.strip(): (v or "").strip() for k, v in r.items() if k}
    return out


def _peek_is_prodigal_ffn(path: Path) -> bool:
    mg = magic_format(path)
    with open(path, "rb") as raw:
        stream = gzip.open(raw, "rb") if mg == "gzip" else raw
        with stream:
            for ln in stream:
                s = ln.strip()
                if not s:
                    continue
                if s.startswith(b">"):
                    return bool(PRODIGAL_FFN_HDR_RE.match(s[1:].decode("utf-8", "replace")))
                return False
    return False


def stage_promge1(
    dirs: dict[str, Path | list[Path]],
    index_dir: Path,
    meta_out: Path,
    fasta_out: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    later_inputs: dict[str, Path],
    inputs_home: Path,
    exclude_promge2_cellular: bool = False,
    promge2_non_nested_only: bool = False,
    promge2_spire_and_pg3_reps_only: bool = False,
    retrieved_promge2: str = PROMGE1_RETRIEVED_PROMGE2,
    retrieved_open: str = PROMGE1_RETRIEVED_OPEN,
    work: Path | None = None,
    threads: int = 2,
    memory_limit: str = "4GB",
) -> dict:
    index_dir = Path(index_dir)
    if not (index_dir / "gtdb_r232.parquet").is_file() or not (index_dir / "globdb_r232.parquet").is_file():
        raise Refused(f"{index_dir} is not card C-74's index (missing gtdb_r232.parquet or globdb_r232.parquet)")

    norm_dirs: dict[str, list[Path]] = {}
    for sub in PROMGE1_SUBCOHORTS:
        if sub not in dirs:
            raise Refused(f"PROMGE1: missing directory mapping for {sub}")
        v = dirs[sub]
        norm_dirs[sub] = [Path(x) for x in v] if isinstance(v, (list, tuple)) else [Path(v)]

    sf = c78.StageFiles(PROMGE1, Path(meta_out), Path(fasta_out), PROMGE1_FASTA, work)
    part_writer = ContigPartWriter(sf.work / "u_parts")
    raw_rec_path = sf.work / "raw_records.tsv"
    genomes_path = sf.tmp / "genomes.json"

    try:
        trims: list[tuple] = []
        left_out: list[tuple] = []
        skipped_members: list[tuple] = []
        repeated: list[tuple] = []
        dtr_classes: Counter[str] = Counter()
        by_subcohort_records: Counter[str] = Counter()
        by_subcohort_genomes: Counter[str] = Counter()
        by_mge_class: Counter[str] = Counter()
        acc_roles: dict[str, tuple[str, str]] = {}
        acc_rec_counts: Counter[str] = Counter()

        used_gnames: set[str] = set()
        seen_contig_circ: dict[bytes, bool] = {}  # enforce consistent circularity across shared contigs
        first_genome_of_guid: dict[bytes, bytes] = {}  # guid_16 -> first genome name UTF-8 bytes
        used_cats: set[str] = set()

        n_genomes = 0
        n_circular_genomes = 0
        n_records = 0
        n_shared_within = 0
        n_repeated_within_genome = 0

        def unique_gname(sub: str, stem: str) -> str:
            cand = stem if stem and stem not in used_gnames else f"{sub}_{stem}"
            k = 2
            while cand in used_gnames:
                cand = f"{sub}_{stem}_{k}"
                k += 1
            used_gnames.add(cand)
            return cand

        raw_rf = open(raw_rec_path, "w", encoding="utf-8", buffering=1 << 20)
        gf = open(genomes_path, "w", encoding="utf-8", buffering=1 << 20)
        gf.write("[")
        first_json = True

        def emit_genome(
            sub: str,
            cat: str,
            gname: str,
            mge_cls: str,
            contigs_data: list[tuple[str, bytes, int, bool, str, int, bytes, int]],
        ) -> None:
            """Write FASTA records, stream raw_records.tsv + Parquet uids, and append genome to genomes.json."""
            nonlocal first_json, n_genomes, n_circular_genomes, n_records, n_shared_within, n_repeated_within_genome
            used_cats.add(cat)
            start_ord = n_records
            uids: list[bytes] = []
            circ_list: list[bytes] = []
            for k, (hdr, out_seq, width, circ, dtr_cls, k_trim, src_cu, src_len) in enumerate(contigs_data):
                cu = ids.contig_uid(out_seq) if k_trim else src_cu
                if cu in seen_contig_circ:
                    circ = seen_contig_circ[cu]
                    n_shared_within += 1
                else:
                    seen_contig_circ[cu] = circ
                first_tok = hdr.split(None, 1)[0] if hdr.strip() else ""
                out_hdr = hdr if (len(contigs_data) == 1 or first_tok.startswith(gname + "_") or first_tok == gname) else f"{gname}|{hdr}"
                written = sf.write(cat, out_hdr, out_seq, width, prefix=sub if len(contigs_data) == 1 else gname)
                if k_trim:
                    trims.append((sub, gname, written, src_cu.hex(), cu.hex(), src_len, k_trim, k_trim, len(out_seq)))
                part_writer.add(cu, len(out_seq))
                raw_rf.write(
                    f"{n_records}\t{sub}\t{gname}\t{k}\t{cu.hex()}\t{len(out_seq)}\t"
                    f"{str(circ).lower()}\t{dtr_cls}\t{k_trim}\t{mge_cls}\n"
                )
                n_records += 1
                by_subcohort_records[sub] += 1
                by_mge_class[mge_cls] += 1
                uids.append(cu)
                if circ:
                    circ_list.append(cu)

            dedup_uids = list(dict.fromkeys(uids))
            n_repeated_within_genome += len(uids) - len(dedup_uids)
            guid = ids.genome_uid(dedup_uids)
            if guid in first_genome_of_guid:
                repeated.append((start_ord, gname, guid.hex(), first_genome_of_guid[guid].decode("utf-8")))
                return
            first_genome_of_guid[guid] = gname.encode("utf-8")
            dedup_circ = [u.hex() for u in dict.fromkeys(circ_list)]
            qual = "complete" if (len(dedup_uids) == 1 and len(dedup_circ) == 1) else "contig"
            entry = {
                "contig_uids": [u.hex() for u in dedup_uids],
                "genetic_code": 11,
                "kind": "plasmid",
                "assembly_quality": qual,
                "circular": dedup_circ,
            }
            if not first_json:
                gf.write(",")
            first_json = False
            gf.write(json.dumps(entry))
            n_genomes += 1
            if dedup_circ:
                n_circular_genomes += 1
            by_subcohort_genomes[sub] += 1

        promge2_excl_counts: Counter[str] = Counter()
        promge2_excl_bases: Counter[str] = Counter()
        promge2_excl_first: dict[str, str] = {}
        pg3_rep_accs: set[str] = set()
        rescued_rep_ids: set[str] = set()
        if promge2_spire_and_pg3_reps_only:
            pg3_rep_accs, rescued_rep_ids = read_promge2_qualifying_hosts(norm_dirs["PROMGE2"])

        def _record_promge2_excl(reason: str, tok: str, seq_len: int) -> None:
            promge2_excl_counts[reason] += 1
            promge2_excl_bases[reason] += seq_len
            promge2_excl_first.setdefault(reason, tok)

        try:
            # 1. PROMGE2: mge100.fna.gz (+ skip tabular .tsv.gz files)
            for p in _all_files_in_dirs(norm_dirs["PROMGE2"]):
                if p.name == "mge100.fna.gz":
                    with _open_gzip_fast(p) as gfh:
                        for acc, hdr, seq, width in iter_fasta_fast(gfh):
                            tok, tags, host_acc = _parse_promge2_header(hdr)
                            mge_val = tags.get("mge", "")
                            mge_type = tags.get("mge_type", "")
                            if not seq or seq.translate(None, IUPAC_BYTES):
                                left_out.append(("PROMGE2", p.name, tok, len(seq), "non_iupac_or_empty"))
                                continue
                            if exclude_promge2_cellular and _is_promge2_cellular_only(mge_val):
                                _record_promge2_excl("promge2_cellular_recombinase_non_mge", tok, len(seq))
                                continue
                            if promge2_non_nested_only and mge_type == "nested":
                                _record_promge2_excl("promge2_nested_composite", tok, len(seq))
                                continue
                            if promge2_spire_and_pg3_reps_only and not (
                                _is_promge2_qualifying_member(tok, pg3_rep_accs) or tok in rescued_rep_ids
                            ):
                                _record_promge2_excl("promge2_pg3_nonrep_host", tok, len(seq))
                                continue
                            if host_acc:
                                acc_roles.setdefault(host_acc, ("PROMGE2", "host_genome"))
                                acc_rec_counts[host_acc] += 1
                            src_cu = ids.contig_uid(seq)
                            mge_cls = f"promge2:{mge_val.split(':', 1)[0]}" if mge_val else "promge2:mge"
                            gname = unique_gname("PROMGE2", tok)
                            emit_genome("PROMGE2", "promge2", gname, mge_cls, [
                                (hdr, seq, width, False, "", 0, src_cu, len(seq))
                            ])
                else:
                    skipped_members.append(("PROMGE2", p.name, p.stat().st_size, "metadata_table_not_fasta"))

            for reason in sorted(promge2_excl_counts):
                left_out.append((
                    "PROMGE2",
                    f"mge100.fna.gz({promge2_excl_counts[reason]}_records)",
                    promge2_excl_first[reason],
                    promge2_excl_bases[reason],
                    reason,
                ))

            # 2. MFD_PLAS: mfd_plasmids_drep_and_members.fasta + pfd_metadata.tsv (skip Prodigal .ffn!)
            mfd_meta_path = _find_in_dirs(norm_dirs["MFD_PLAS"], "pfd_metadata.tsv")
            mfd_topo = read_mfd_plasmid_metadata(mfd_meta_path) if mfd_meta_path is not None else {}
            for p in _all_files_in_dirs(norm_dirs["MFD_PLAS"]):
                if not _looks_dna_fasta(p.name):
                    skipped_members.append(("MFD_PLAS", p.name, p.stat().st_size, "metadata_table_not_fasta"))
                    continue
                if _peek_is_prodigal_ffn(p):
                    skipped_members.append(("MFD_PLAS", p.name, p.stat().st_size, "prodigal_cds_ffn_not_plasmid_contig"))
                    continue
                mg = magic_format(p)
                with open(p, "rb") as raw:
                    stream = gzip.open(raw, "rb") if mg == "gzip" else raw
                    with stream:
                        for acc, hdr, seq, width in iter_fasta_fast(stream):
                            if not seq or seq.translate(None, IUPAC_BYTES):
                                left_out.append(("MFD_PLAS", p.name, acc, len(seq), "non_iupac_or_empty"))
                                continue
                            src_cu = ids.contig_uid(seq)
                            topo = mfd_topo.get(acc, "")
                            if topo in ("Inverted terminal repeat", "Linear"):
                                # Preserve biological TIRs on linear replicons byte-identically (DEC-65)
                                out_seq, k_trim, cls, circ = seq, 0, "linear_or_tir_preserved", False
                            else:
                                flagged = topo in ("Circular", "Direct terminal repeat") or bool(CIRCULAR_RE.search(hdr))
                                out_seq, k_trim, cls, circ = c78.dec65(seq, flagged=flagged)
                                dtr_classes[cls] += 1
                            gname = unique_gname("MFD_PLAS", acc)
                            emit_genome("MFD_PLAS", "open", gname, "mfd_plasmid", [
                                (hdr, out_seq, width, circ, cls, k_trim, src_cu, len(seq))
                            ])

            # 3. GSPR: all_soil_plasmids_derep.fna.gz
            for p in _all_files_in_dirs(norm_dirs["GSPR"]):
                if p.name == "all_soil_plasmids_derep.fna.gz":
                    with open(p, "rb") as raw, gzip.open(raw, "rb") as gfh:
                        for acc, hdr, seq, width in iter_fasta_fast(gfh):
                            if not seq or seq.translate(None, IUPAC_BYTES):
                                left_out.append(("GSPR", p.name, acc, len(seq), "non_iupac_or_empty"))
                                continue
                            src_cu = ids.contig_uid(seq)
                            out_seq, k_trim, cls, circ = c78.dec65(seq, flagged=bool(CIRCULAR_RE.search(hdr)))
                            dtr_classes[cls] += 1
                            gname = unique_gname("GSPR", acc)
                            emit_genome("GSPR", "open", gname, "soil_plasmid", [
                                (hdr, out_seq, width, circ, cls, k_trim, src_cu, len(seq))
                            ])
                else:
                    skipped_members.append(("GSPR", p.name, p.stat().st_size, "not_fasta"))

            # 4. ECE_HODOR: Hodor.38.genomes.fasta, Hodor.prophages.fasta, prophage01..03.all.pure.fa, draft.Hodor.genomes.tar.gz
            hodor_files = _all_files_in_dirs(norm_dirs["ECE_HODOR"])
            hodor_order = sorted(
                hodor_files,
                key=lambda x: (
                    0 if x.name == "Hodor.38.genomes.fasta" else (
                        1 if x.name == "Hodor.prophages.fasta" else (
                            2 if x.name.startswith("prophage") else 3
                        )
                    ),
                    x.name,
                ),
            )
            for p in hodor_order:
                if p.name == "draft.Hodor.genomes.tar.gz":
                    with tarfile.open(p, "r|*") as tf:
                        for m in tf:
                            if m.isdir():
                                continue
                            if _is_macosx_meta(m.name) or not _looks_dna_fasta(m.name):
                                skipped_members.append(("ECE_HODOR", f"{p.name}:{m.name}", m.size, "not_fasta"))
                                continue
                            ef = tf.extractfile(m)
                            assert ef is not None
                            with _open_maybe_compressed(ef, m.name) as stream:
                                recs_bin: list[tuple[str, bytes, int, bytes]] = []
                                for acc, hdr, seq, width in iter_fasta_fast(stream):
                                    if not seq or seq.translate(None, IUPAC_BYTES):
                                        left_out.append(("ECE_HODOR", f"{p.name}:{m.name}", acc, len(seq), "non_iupac_or_empty"))
                                        continue
                                    recs_bin.append((hdr, seq, width, ids.contig_uid(seq)))
                            if not recs_bin:
                                skipped_members.append(("ECE_HODOR", f"{p.name}:{m.name}", m.size, "no_valid_dna_records"))
                                continue
                            gname = unique_gname("ECE_HODOR", _stem_clean(m.name))
                            if len(recs_bin) == 1:
                                hdr, seq, width, src_cu = recs_bin[0]
                                out_seq, k_trim, cls, circ = c78.dec65(seq, flagged=bool(CIRCULAR_RE.search(hdr)))
                                dtr_classes[cls] += 1
                                emit_genome("ECE_HODOR", "open", gname, "hodor_draft", [
                                    (hdr, out_seq, width, circ, cls, k_trim, src_cu, len(seq))
                                ])
                            else:
                                emit_genome("ECE_HODOR", "open", gname, "hodor_draft", [
                                    (hdr, seq, width, False, "", 0, src_cu, len(seq))
                                    for hdr, seq, width, src_cu in recs_bin
                                ])
                elif _looks_dna_fasta(p.name):
                    mcls = "hodor_prophage" if "prophage" in p.name.lower() else "hodor_ece"
                    mg = magic_format(p)
                    with open(p, "rb") as raw:
                        stream = gzip.open(raw, "rb") if mg == "gzip" else raw
                        with stream:
                            for acc, hdr, seq, width in iter_fasta_fast(stream):
                                if not seq or seq.translate(None, IUPAC_BYTES):
                                    left_out.append(("ECE_HODOR", p.name, acc, len(seq), "non_iupac_or_empty"))
                                    continue
                                src_cu = ids.contig_uid(seq)
                                out_seq, k_trim, cls, circ = c78.dec65(seq, flagged=bool(CIRCULAR_RE.search(hdr)))
                                dtr_classes[cls] += 1
                                gname = unique_gname("ECE_HODOR", acc)
                                emit_genome("ECE_HODOR", "open", gname, mcls, [
                                    (hdr, out_seq, width, circ, cls, k_trim, src_cu, len(seq))
                                ])
                else:
                    skipped_members.append(("ECE_HODOR", p.name, p.stat().st_size, "not_fasta"))

            # 5. SAT_P4P2: P4P2_DereplicatedGenomes.zip (Replicons/*.fna only; skip Proteins/*.prt and __MACOSX)
            for p in _all_files_in_dirs(norm_dirs["SAT_P4P2"]):
                if p.name == "P4P2_DereplicatedGenomes.zip":
                    with zipfile.ZipFile(p) as zf:
                        for info in sorted(zf.infolist(), key=lambda x: x.filename):
                            if info.is_dir():
                                continue
                            if _is_macosx_meta(info.filename):
                                skipped_members.append(("SAT_P4P2", info.filename, info.file_size, "macosx"))
                                continue
                            if "/Replicons/" not in info.filename or not _looks_dna_fasta(info.filename):
                                why = "protein_fasta" if ("/Proteins/" in info.filename or info.filename.lower().endswith(".prt")) else "not_replicon_fasta"
                                skipped_members.append(("SAT_P4P2", info.filename, info.file_size, why))
                                continue
                            mcls = "p4_satellite" if "/P4/" in info.filename else "p2_prophage"
                            with zf.open(info) as raw:
                                recs_sat: list[tuple[str, bytes, int, bytes]] = []
                                for acc, hdr, seq, width in iter_fasta_fast(raw):
                                    if not seq or seq.translate(None, IUPAC_BYTES):
                                        left_out.append(("SAT_P4P2", info.filename, acc, len(seq), "non_iupac_or_empty"))
                                        continue
                                    recs_sat.append((hdr, seq, width, ids.contig_uid(seq)))
                            if not recs_sat:
                                skipped_members.append(("SAT_P4P2", info.filename, info.file_size, "no_valid_dna_records"))
                                continue
                            gname = unique_gname("SAT_P4P2", _stem_clean(info.filename))
                            emit_genome("SAT_P4P2", "open", gname, mcls, [
                                (hdr, seq, width, False, "", 0, src_cu, len(seq))
                                for hdr, seq, width, src_cu in recs_sat
                            ])
                else:
                    skipped_members.append(("SAT_P4P2", p.name, p.stat().st_size, "not_fasta"))

            # 6. ECE_METH: Mp_and_ECE_genomes.fa.gz (stage HMp-Plasmid-like / Mp-ICE ECEs; leave out host chromosomes)
            for p in _all_files_in_dirs(norm_dirs["ECE_METH"]):
                if _looks_dna_fasta(p.name):
                    mg = magic_format(p)
                    with open(p, "rb") as raw:
                        stream = gzip.open(raw, "rb") if mg == "gzip" else raw
                        with stream:
                            for acc, hdr, seq, width in iter_fasta_fast(stream):
                                if not seq or seq.translate(None, IUPAC_BYTES):
                                    left_out.append(("ECE_METH", p.name, acc, len(seq), "non_iupac_or_empty"))
                                    continue
                                low_a = acc.lower()
                                if not any(k in low_a for k in ("plasmid", "ece", "ice", "borg")):
                                    left_out.append(("ECE_METH", p.name, acc, len(seq), "host_archaeal_chromosome_not_ece"))
                                    continue
                                src_cu = ids.contig_uid(seq)
                                is_ice = "ice" in low_a
                                flagged = (not is_ice) and ("complete" in low_a or bool(CIRCULAR_RE.search(hdr)))
                                out_seq, k_trim, cls, circ = (seq, 0, "", False) if is_ice else c78.dec65(seq, flagged=flagged)
                                if cls:
                                    dtr_classes[cls] += 1
                                gname = unique_gname("ECE_METH", acc)
                                emit_genome("ECE_METH", "open", gname, "mp_ice" if is_ice else "mp_ece", [
                                    (hdr, out_seq, width, circ, cls, k_trim, src_cu, len(seq))
                                ])
                else:
                    skipped_members.append(("ECE_METH", p.name, p.stat().st_size, "not_fasta"))

            # 7. PICJI_SAT: INSDC BioProject PRJNA1424553 complete assemblies (DEC-65 (2): byte-identical)
            picji_root = norm_dirs["PICJI_SAT"][0]
            picji_mans = sorted(picji_root.glob("MANIFEST_PICJI_SAT_*.tsv"))
            if not picji_mans:
                raise Refused(f"PICJI_SAT: missing MANIFEST_PICJI_SAT_*.tsv in {picji_root}")
            picji_meta = read_ncbi_batch6b_manifest(picji_mans[0])
            picji_asm = picji_root / "assemblies"
            for acc in sorted(picji_meta):
                matches = sorted(picji_asm.glob(f"{acc}*_genomic.fna.gz")) + sorted(picji_asm.glob(f"{acc}*.fna.gz"))
                if not matches:
                    raise Refused(f"PICJI_SAT: missing assembly FASTA for {acc} in {picji_asm}")
                fna_gz = matches[0]
                acc_roles[acc] = ("PICJI_SAT", "insdc_assembly")
                with gzip.open(fna_gz, "rb") as gfh:
                    for racc, hdr, seq, width in iter_fasta_fast(gfh):
                        if not seq or seq.translate(None, IUPAC_BYTES):
                            left_out.append(("PICJI_SAT", fna_gz.name, racc, len(seq), "non_iupac_or_empty"))
                            continue
                        src_cu = ids.contig_uid(seq)
                        low_h = hdr.lower()
                        circ = ("circular" in low_h or "complete sequence" in low_h or "complete genome" in low_h) and "incomplete" not in low_h
                        mcls = "picji_plasmid" if "plasmid" in low_h else "picji_chromosome"
                        gname = unique_gname("PICJI_SAT", f"{acc}_{racc}")
                        acc_rec_counts[acc] += 1
                        emit_genome("PICJI_SAT", "open", gname, mcls, [
                            (hdr, seq, width, circ, "", 0, src_cu, len(seq))
                        ])

            # 8. BORG_ONT: INSDC BioProject PRJNA1119519 Nanopore Borg assembly (DEC-65 (2): linear, byte-identical)
            borg_root = norm_dirs["BORG_ONT"][0]
            borg_mans = sorted(borg_root.glob("MANIFEST_BORG_ONT_*.tsv"))
            if not borg_mans:
                raise Refused(f"BORG_ONT: missing MANIFEST_BORG_ONT_*.tsv in {borg_root}")
            borg_meta = read_ncbi_batch6b_manifest(borg_mans[0])
            borg_asm = borg_root / "assemblies"
            for acc in sorted(borg_meta):
                matches = sorted(borg_asm.glob(f"{acc}*_genomic.fna.gz")) + sorted(borg_asm.glob(f"{acc}*.fna.gz"))
                if not matches:
                    raise Refused(f"BORG_ONT: missing assembly FASTA for {acc} in {borg_asm}")
                fna_gz = matches[0]
                acc_roles[acc] = ("BORG_ONT", "insdc_borg_assembly")
                recs_borg: list[tuple[str, bytes, int, bytes]] = []
                with gzip.open(fna_gz, "rb") as gfh:
                    for racc, hdr, seq, width in iter_fasta_fast(gfh):
                        if not seq or seq.translate(None, IUPAC_BYTES):
                            left_out.append(("BORG_ONT", fna_gz.name, racc, len(seq), "non_iupac_or_empty"))
                            continue
                        recs_borg.append((hdr, seq, width, ids.contig_uid(seq)))
                if not recs_borg:
                    raise Refused(f"BORG_ONT: no valid DNA records in {fna_gz}")
                acc_rec_counts[acc] += len(recs_borg)
                gname = unique_gname("BORG_ONT", acc)
                emit_genome("BORG_ONT", "open", gname, "borg_ont", [
                    (hdr, seq, width, False, "", 0, src_cu, len(seq))
                    for hdr, seq, width, src_cu in recs_borg
                ])
        finally:
            gf.write("]")
            gf.close()
            raw_rf.close()

        if n_genomes == 0:
            raise Refused("PROMGE1: no MGE/plasmid genomes staged")

        sf.close()
        part_writer.flush()
        own_rep = _duck_ownership_from_parts(
            part_writer.parts, root, root_stages, owner_inputs, later_inputs, sf.work, threads, memory_limit
        )
        owners = own_rep["owners"]
        later = own_rep["later"]

        # Stream final records.tsv + shared_records.tsv + later_records.tsv in O(1) RAM
        shared_rows: list[tuple] = []
        later_rows: list[tuple] = []
        with open(raw_rec_path, "r", encoding="utf-8", buffering=1 << 20) as r_in, open(
            sf.tmp / "records.tsv", "w", encoding="utf-8", buffering=1 << 20
        ) as r_out:
            r_out.write("\t".join(PROMGE1_RECORDS_HEAD) + "\n")
            for ln in r_in:
                cols = ln.rstrip("\n").split("\t")
                ordinal_s, sub, gname, k_s, cu_hex, length_s = cols[0], cols[1], cols[2], cols[3], cols[4], cols[5]
                cu_b = bytes.fromhex(cu_hex)
                own = owners.get(cu_b, "")
                lat = "" if own else later.get(cu_b, "")
                r_out.write("\t".join([*cols, own, lat]) + "\n")
                if own:
                    shared_rows.append((int(ordinal_s), gname, cu_hex, int(length_s), own))
                elif lat:
                    later_rows.append((int(ordinal_s), gname, cu_hex, int(length_s), lat))
        raw_rec_path.unlink(missing_ok=True)

        # Lookup INSDC host/assembly accessions in C-74's r232_index
        nover_to_accs: dict[str, list[str]] = {}
        for acc in sorted(acc_roles):
            nover_to_accs.setdefault(r232_index.nover(acc), []).append(acc)
        look_map = {
            r232_index.nover(x["accession"]): x
            for x in (r232_index.lookup(index_dir, sorted(nover_to_accs)) if nover_to_accs else [])
        }
        status_rows: list[tuple] = []
        status_counts: Counter[str] = Counter()
        for acc in sorted(acc_roles):
            sub, role = acc_roles[acc]
            x = look_map.get(r232_index.nover(acc))
            st_lbl = _r232_label(x, True)
            status_counts[st_lbl] += 1
            status_rows.append((
                acc, sub, role, st_lbl,
                "yes" if (x and x["in_gtdb_r232"]) else "no",
                str(x["is_representative"]) if (x and x["in_gtdb_r232"]) else "",
                (x["representative"] or "") if x else "",
                (x["species"] or "") if x else "",
                "yes" if (x and x["in_globdb_r232"]) else "no",
                (x["globdb_r232_id"] or "") if x else "",
                acc_rec_counts[acc],
            ))

        make_inputs._write_tsv(sf.tmp / "left_out.tsv", PROMGE1_LEFT_HEAD, left_out)
        make_inputs._write_tsv(sf.tmp / "skipped_members.tsv", PROMGE1_SKIP_HEAD, skipped_members)
        make_inputs._write_tsv(sf.tmp / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
        make_inputs._write_tsv(sf.tmp / "shared_records.tsv", ["ordinal", "genome", "contig_uid", "length", "owner_stage"], shared_rows)
        make_inputs._write_tsv(sf.tmp / "later_records.tsv", ["ordinal", "genome", "contig_uid", "length", "holder"], later_rows)
        make_inputs._write_tsv(sf.tmp / "trims.tsv", c62.TRIM_HEAD, trims)
        make_inputs._write_tsv(sf.tmp / "r232_status.tsv", PROMGE1_STATUS_HEAD, status_rows)

        used_ordered = [k for k in PROMGE1_FASTA if k in used_cats]
        promge2_desc = (
            "proMGE v2.0 recombinase-delineated prokaryotic MGE catalog (promge.embl.de mge100.fna.gz; non-cellular, non-nested MGEs across SPIRE MAGs and proGenomes3 specI species representatives); staged by c89_promge1_bbsag1_inputs.py (DEC-65, DEC-80 (3), DEC-88)"
            if promge2_spire_and_pg3_reps_only
            else "proMGE v2.0 recombinase-delineated prokaryotic MGE catalog (promge.embl.de mge100.fna.gz); staged by c89_promge1_bbsag1_inputs.py (DEC-65, DEC-80 (3), DEC-88)"
        )
        source_meta = {
            "promge2": {
                "source_name": "tower_promge2",
                "source_release": promge2_desc,
                "license": "CC BY-SA 4.0 as the job list states it (sources.yaml tower_promge2: no provider record of the licence was kept)",
                "redistributable": False,
                "retrieved_at": retrieved_promge2,
            },
            "open": {
                "source_name": "tower_promge1_plasmids_eces",
                "source_release": "Global Soil Plasmidome Resource (zenodo:17136302), Microflora Danica ONT plasmids (zenodo:22944851), Hodor >300 kb megaplasmids & prophages (figshare:31385365), P4/P2 satellites & prophages (figshare:31843870), Methanoperedens ECEs (figshare:32069694), and INSDC PRJNA1424553 & PRJNA1119519; staged by c89_promge1_bbsag1_inputs.py (DEC-65, DEC-80 (3), DEC-88)",
                "license": "CC BY 4.0 / cc-by-4.0 / INSDC (confirmed in sources.yaml tower_gspr, tower_mfd_plas, tower_ece_hodor, tower_sat_p4p2, tower_ece_meth, tower_picji_sat, tower_borg_ont)",
                "redistributable": True,
                "retrieved_at": retrieved_open,
            },
        }
        _write_objects_json(sf, used_ordered, source_meta, Path(inputs_home))
        sha_map, byte_map, locations = _commit_with_stor_floor(sf, used_ordered)

        return {
            "stage": PROMGE1,
            "records": n_records,
            "records_by_subcohort": dict(sorted(by_subcohort_records.items())),
            "genomes": n_genomes,
            "genomes_by_subcohort": dict(sorted(by_subcohort_genomes.items())),
            "repeated_genomes": len(repeated),
            "circular_genomes": n_circular_genomes,
            "by_mge_class": dict(sorted(by_mge_class.items())),
            "dtr_classes": dict(sorted(dtr_classes.items())),
            "trims": len(trims),
            "trimmed_bases": sum(t[7] for t in trims),
            "left_out": len(left_out),
            "left_out_by_reason": dict(sorted(Counter(x[4] for x in left_out).items())),
            "promge2_excluded_records": dict(sorted(promge2_excl_counts.items())),
            "promge2_excluded_bases": dict(sorted(promge2_excl_bases.items())),
            "promge2_rescued_clusters": len(rescued_rep_ids),
            "skipped_members": len(skipped_members),
            "r232_host_accession_status": dict(sorted(status_counts.items())),
            "contigs_distinct": own_rep["contigs_distinct"],
            "bases_distinct": own_rep["bases_distinct"],
            "held_by": own_rep["held_by"],
            "held_bases_by": own_rep["held_bases_by"],
            "later_held_by": own_rep["later_held_by"],
            "later_held_bases_by": own_rep["later_held_bases_by"],
            "contigs_new": own_rep["contigs_new"],
            "bases_new": own_rep["bases_new"],
            "contigs_new_after_later": own_rep["contigs_new_after_later"],
            "bases_new_after_later": own_rep["bases_new_after_later"],
            "shared_within_stage_occurrences": n_shared_within,
            "repeated_within_genome_contigs": n_repeated_within_genome,
            "sha256": sha_map,
            "bytes": byte_map,
            "locations": locations,
        }
    except BaseException:
        sf.abort()
        raise


# ==================================================================================================================
# 5. Stage BBSAG1 (or BBSAG1A + BBSAG1B if net-new bases > max_stage_new_bp)
# ==================================================================================================================

BBSAG1 = "BBSAG1"
BBSAG1_MAX_NEW_BP = 40_000_000_000  # 40 Gbp net-new cap per stage
BBSAG1_RETRIEVED_MAIN = "2026-11-27T05:51:55Z"
BBSAG1_RETRIEVED_DV = "2026-11-27T06:14:46Z"

BBSAG1_RECORDS_HEAD = [
    "ordinal", "subcohort", "genome", "record_ordinal", "contig_uid", "length",
    "assembly_quality", "genetic_code", "circular", "owner_stage", "later_holder",
]
BBSAG1_STATUS_HEAD = [
    "genome", "subcohort", "genome_kind", "accession", "r232_status", "in_gtdb_r232",
    "is_gtdb_representative", "gtdb_representative", "gtdb_species",
    "in_globdb_r232", "globdb_r232_id", "taxonomy", "completeness", "contamination",
    "cluster_id", "contigs", "bp",
]
BBSAG1_LEFT_HEAD = ["subcohort", "genome", "sagmag", "completeness", "contamination", "length", "reason"]
BBSAG1_SKIP_HEAD = ["subcohort", "member", "bytes", "reason"]

_BBSAG20_MEMBER_RE = re.compile(r"^(?:mag/)?([A-Za-z0-9]+)/\1-(\d+)/\1-\2(\.[A-Za-z0-9.]+)$")


def read_bbsag20_summary(path: Path) -> tuple[dict[str, dict], list[tuple]]:
    """Read BBSAG20 `summary.tsv`, keeping MIMAG MQ/HQ (`completeness >= 50.0` and `contamination < 10.0`)
    and recording LQ genomes in `left_out`."""
    kept: dict[str, dict] = {}
    left: list[tuple] = []
    need = ("sampleid", "sagmag", "clusterid", "genomeid", "completeness", "contamination", "totallength", "gtdbr207")
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        miss = [c for c in need if c not in (rd.fieldnames or [])]
        if miss:
            raise Refused(f"{path}: missing column {miss[0]}")
        for r in rd:
            gid = (r.get("genomeid") or "").strip()
            if not gid:
                continue
            if gid in kept:
                raise Refused(f"{path}: duplicate genomeid {gid}")
            sagmag = (r.get("sagmag") or "").strip().upper()
            comp_s = (r.get("completeness") or "").strip()
            cont_s = (r.get("contamination") or "").strip()
            len_s = (r.get("totallength") or "0").strip()
            try:
                comp = float(comp_s)
                cont = float(cont_s)
                tot_len = int(float(len_s))
            except ValueError as exc:
                raise Refused(f"{path}: invalid numeric fields on {gid}") from exc
            if comp >= 50.0 and cont < 10.0:
                tier = "HQ" if (comp >= 90.0 and cont < 5.0) else "MQ"
                qual = ("mag_hq" if tier == "HQ" else "mag_mq") if sagmag == "MAG" else "contig"
                kept[gid] = {
                    "sagmag": sagmag,
                    "tier": tier,
                    "quality": qual,
                    "cluster_id": (r.get("clusterid") or "").strip(),
                    "completeness": comp_s,
                    "contamination": cont_s,
                    "taxonomy": (r.get("gtdbr207") or "").strip(),
                }
            else:
                left.append(("BBSAG20", gid, sagmag, comp_s, cont_s, tot_len, "bbsag20_low_quality"))
    return kept, left


def _bbsag20_member_genomeid(member_name: str, is_mag_tar: bool) -> tuple[str | None, str]:
    m = _BBSAG20_MEMBER_RE.match(member_name)
    if m is None:
        return None, os.path.basename(member_name)
    prefix, num, ext = m.group(1), m.group(2), m.group(3)
    kind_tag = "mag" if is_mag_tar or member_name.startswith("mag/") else "sag"
    return f"{prefix}{kind_tag}-{num}", ext


def _fasta_names_for_bbsag_stage(stage_name: str) -> dict[str, str]:
    low = stage_name.lower()
    return {
        "main": f"{low}_sags_mags.fna.gz",
        "dv": f"{low}_deathvalley.fna.gz",
    }


def _bgzf_batch_worker(chunks: list[bytes]) -> bytes:
    from biobox.l0seq import bgzf

    return b"".join(bgzf._block(c, 6) for c in chunks)


class ParallelBgzfWriter:
    """Byte-identical to `c77.BgzfWriter(path, level=6)`, spreading deflate blocks across `pool`."""

    def __init__(self, path: Path, pool: ProcessPoolExecutor) -> None:
        from biobox.l0seq import bgzf

        self._bgzf = bgzf
        self.fh = open(path, "wb", buffering=1 << 22)
        self.buf = bytearray()
        self.pool = pool
        self._chunks: list[bytes] = []
        self._futs: list[object] = []

    def _submit_batch(self) -> None:
        if not self._chunks:
            return
        batch = self._chunks
        self._chunks = []
        self._futs.append(self.pool.submit(_bgzf_batch_worker, batch))
        while len(self._futs) >= 8:
            self.fh.write(self._futs.pop(0).result())

    def write(self, data: bytes) -> None:
        bs = self._bgzf.BLOCK_DATA
        self.buf += data
        while len(self.buf) >= bs:
            self._chunks.append(bytes(self.buf[:bs]))
            del self.buf[:bs]
            if len(self._chunks) >= 128:
                self._submit_batch()

    def close(self) -> None:
        if self.buf:
            self._chunks.append(bytes(self.buf))
            self.buf.clear()
        self._submit_batch()
        while self._futs:
            self.fh.write(self._futs.pop(0).result())
        self.fh.write(self._bgzf.EOF_BLOCK)
        self.fh.close()


def stage_bbsag1(
    bbsag20_dirs: list[Path] | Path,
    haima_dirs: list[Path] | Path,
    dv_dirs: list[Path] | Path,
    index_dir: Path,
    meta_out: Path,
    fasta_out: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    later_inputs: dict[str, Path],
    inputs_home: Path,
    max_stage_new_bp: int = BBSAG1_MAX_NEW_BP,
    retrieved_main: str = BBSAG1_RETRIEVED_MAIN,
    retrieved_dv: str = BBSAG1_RETRIEVED_DV,
    work: Path | None = None,
    threads: int = 2,
    memory_limit: str = "4GB",
    workers: int = 1,
) -> dict:
    index_dir = Path(index_dir)
    if not (index_dir / "gtdb_r232.parquet").is_file() or not (index_dir / "globdb_r232.parquet").is_file():
        raise Refused(f"{index_dir} is not card C-74's index (missing gtdb_r232.parquet or globdb_r232.parquet)")

    bb_roots = [Path(x) for x in bbsag20_dirs] if isinstance(bbsag20_dirs, (list, tuple)) else [Path(bbsag20_dirs)]
    hm_roots = [Path(x) for x in haima_dirs] if isinstance(haima_dirs, (list, tuple)) else [Path(haima_dirs)]
    dv_roots = [Path(x) for x in dv_dirs] if isinstance(dv_dirs, (list, tuple)) else [Path(dv_dirs)]

    summary_path = _find_in_dirs(bb_roots, "summary.tsv")
    if summary_path is None:
        raise Refused("BBSAG20: missing summary.tsv")
    bb_kept, left_out = read_bbsag20_summary(summary_path)

    work_dir = Path(work) if work is not None else Path(meta_out) / "_bbsag1_work"
    own_work = work is None
    work_dir.mkdir(parents=True, exist_ok=True)
    spool_path = work_dir / "bbsag1_records_spool.bin"
    part_writer = ContigPartWriter(work_dir / "u_parts")

    pool = ProcessPoolExecutor(max_workers=workers, max_tasks_per_child=TASKS_PER_CHILD) if workers > 1 else None
    pending: list[tuple[str, object]] = []

    try:
        gc_results: dict[str, dict] = {}
        contigs: dict[str, tuple[bytes, list[int]]] = {}
        spool_index: dict[str, list[tuple[int, int, int, bool]]] = {}  # gid -> [(offset, byte_len, seq_len, circ)]
        ginfo: dict[str, dict] = {}
        skipped_members: list[tuple] = []
        trims: list[tuple] = []
        dtr_classes: Counter[str] = Counter()
        order_gids: list[str] = []

        def submit_screen(gid: str, seqs: list[bytes], tax: str, force_alt: bool = False) -> None:
            if pool is None:
                gc_results[gid] = screen_prokaryote_code(seqs, tax, force_alt=force_alt)
                return
            pending.append((gid, pool.submit(_prok_screen_worker, (gid, seqs, tax, force_alt))))
            while len(pending) >= 4 * workers:
                s, fut = pending.pop(0)
                gc_results[s] = fut.result()[1]

        with open(spool_path, "wb", buffering=1 << 22) as spool_fh:

            def ingest_genome(
                subcohort: str,
                cat: str,
                gid: str,
                member: str,
                fobj,
                *,
                genome_kind: str,
                quality: str,
                tax: str,
                comp: str,
                cont: str,
                cluster_id: str,
                apply_single_dtr: bool = False,
            ) -> None:
                raw_recs: list[tuple[str, bytes, int]] = []
                for acc, hdr, seq, width in iter_fasta_fast(fobj):
                    if not seq or seq.translate(None, IUPAC_BYTES):
                        left_out.append((subcohort, f"{gid}:{acc}", genome_kind, comp, cont, len(seq), "non_iupac_or_empty"))
                        continue
                    raw_recs.append((hdr, seq, width))
                if not raw_recs:
                    skipped_members.append((subcohort, member, 0, "no_valid_dna_records"))
                    return
                if gid in contigs:
                    raise Refused(f"BBSAG1: duplicate genome identifier {gid}")

                cus: list[bytes] = []
                lens: list[int] = []
                seqs_screen: list[bytes] = []
                sp_entries: list[tuple[int, int, int, bool]] = []
                all_circ = False
                acc_found = ""

                for hdr, seq, width in raw_recs:
                    src_cu = ids.contig_uid(seq)
                    if apply_single_dtr and len(raw_recs) == 1:
                        out_seq, k_trim, cls, circ = c78.dec65(seq, flagged=True)
                        dtr_classes[cls] += 1
                        cu = ids.contig_uid(out_seq) if k_trim else src_cu
                        if k_trim:
                            trims.append((subcohort, gid, hdr.split()[0], src_cu.hex(), cu.hex(), len(seq), k_trim, k_trim, len(out_seq)))
                        if circ:
                            all_circ = True
                    else:
                        out_seq, circ = seq, False
                        cu = src_cu

                    first_tok = hdr.split(None, 1)[0] if hdr.strip() else ""
                    prefixed = hdr if first_tok.startswith(gid + "_") or first_tok == gid else f"{gid}|{hdr}"
                    formatted = b"\n".join(fasta.format_record(prefixed, out_seq, width)) + b"\n"
                    off = spool_fh.tell()
                    spool_fh.write(formatted)
                    sp_entries.append((off, len(formatted), len(out_seq), circ))
                    cus.append(cu)
                    lens.append(len(out_seq))
                    seqs_screen.append(out_seq)
                    part_writer.add(cu, len(out_seq))
                    if not acc_found:
                        m_acc = ACC_GCA_RE.search(hdr)
                        if m_acc:
                            acc_found = m_acc.group(0)

                eff_qual = "complete" if all_circ else quality
                contigs[gid] = (b"".join(cus), lens)
                spool_index[gid] = sp_entries
                ginfo[gid] = {
                    "subcohort": subcohort,
                    "category": cat,
                    "genome_kind": genome_kind,
                    "quality": eff_qual,
                    "accession": acc_found,
                    "taxonomy": tax,
                    "completeness": comp,
                    "contamination": cont,
                    "cluster_id": cluster_id,
                }
                order_gids.append(gid)
                submit_screen(gid, seqs_screen, tax)

            # 1. BBSAG20: *sag.tar.gz + mag.tar.gz
            bb_files = _all_files_in_dirs(bb_roots)
            for p in bb_files:
                if p.name == "summary.tsv":
                    continue
                if p.name == "bbsag20_plasmid.fna.gz":
                    skipped_members.append(("BBSAG20", p.name, p.stat().st_size, "plasmid_subset_of_bbsag20_genomes"))
                    continue
                if not p.name.endswith(".tar.gz"):
                    skipped_members.append(("BBSAG20", p.name, p.stat().st_size, "not_genome_tarball"))
                    continue
                is_mag_tar = p.name == "mag.tar.gz"
                skipped_non_fna = 0
                skipped_non_fna_bytes = 0
                skipped_lq_fna = 0
                skipped_lq_bytes = 0
                with tarfile.open(p, "r|*") as tf:
                    for m in tf:
                        if m.isdir():
                            continue
                        gid, ext = _bbsag20_member_genomeid(m.name, is_mag_tar)
                        if gid is None or ext.lower() != ".fna.gz":
                            skipped_non_fna += 1
                            skipped_non_fna_bytes += m.size
                            continue
                        if gid not in bb_kept:
                            # LQ genome already recorded in left_out from summary.tsv
                            skipped_lq_fna += 1
                            skipped_lq_bytes += m.size
                            continue
                        ef = tf.extractfile(m)
                        assert ef is not None
                        meta_row = bb_kept[gid]
                        with gzip.open(ef, "rb") as gfh:
                            ingest_genome(
                                "BBSAG20", "main", gid, f"{p.name}:{m.name}", gfh,
                                genome_kind=meta_row["sagmag"],
                                quality=meta_row["quality"],
                                tax=meta_row["taxonomy"],
                                comp=meta_row["completeness"],
                                cont=meta_row["contamination"],
                                cluster_id=meta_row["cluster_id"],
                                apply_single_dtr=False,
                            )
                if skipped_non_fna:
                    skipped_members.append((
                        "BBSAG20", f"{p.name}:non_fna_annotations({skipped_non_fna})",
                        skipped_non_fna_bytes, "annotation_faa_ffn_gff_not_genome_fna",
                    ))
                if skipped_lq_fna:
                    skipped_members.append((
                        "BBSAG20", f"{p.name}:lq_fna_members({skipped_lq_fna})",
                        skipped_lq_bytes, "low_quality_excluded_in_left_out",
                    ))

            # 2. HAIMA_LR_SAG: SAGs.zip, long-read-MAGs.zip, short-read-MAGs.zip
            for p in _all_files_in_dirs(hm_roots):
                if p.suffix.lower() != ".zip":
                    skipped_members.append(("HAIMA_LR_SAG", p.name, p.stat().st_size, "not_zip"))
                    continue
                low_p = p.name.lower()
                if "sag" in low_p:
                    gkind, def_qual = "SAG", "contig"
                elif "long-read" in low_p:
                    gkind, def_qual = "LR_MAG", "mag_mq"
                else:
                    gkind, def_qual = "SR_MAG", "mag_mq"
                with zipfile.ZipFile(p) as zf:
                    for info in sorted(zf.infolist(), key=lambda x: x.filename):
                        if info.is_dir():
                            continue
                        if _is_macosx_meta(info.filename) or not _looks_dna_fasta(info.filename):
                            skipped_members.append(("HAIMA_LR_SAG", f"{p.name}:{info.filename}", info.file_size, "not_dna_fasta"))
                            continue
                        stem = _stem_clean(info.filename)
                        gid = stem if stem not in contigs else f"HAIMA_{gkind}_{stem}"
                        is_cbin = (gkind == "LR_MAG") and ("-cbin" in stem.lower() or stem.lower().startswith("cbin"))
                        with zf.open(info) as raw, _open_maybe_compressed(raw, info.filename) as stream:
                            ingest_genome(
                                "HAIMA_LR_SAG", "main", gid, f"{p.name}:{info.filename}", stream,
                                genome_kind="cMAG" if is_cbin else gkind,
                                quality=def_qual,
                                tax="",
                                comp="",
                                cont="",
                                cluster_id="",
                                apply_single_dtr=is_cbin,
                            )

            # 3. DEATHVALLEY_SAG: BLM1-MAGs_SAGs_2011_2021.zip
            for p in _all_files_in_dirs(dv_roots):
                if p.suffix.lower() != ".zip":
                    skipped_members.append(("DEATHVALLEY_SAG", p.name, p.stat().st_size, "not_zip"))
                    continue
                with zipfile.ZipFile(p) as zf:
                    for info in sorted(zf.infolist(), key=lambda x: (x.filename, x.CRC)):
                        if info.is_dir():
                            continue
                        if _is_macosx_meta(info.filename) or not _looks_dna_fasta(info.filename):
                            skipped_members.append(("DEATHVALLEY_SAG", f"{p.name}:{info.filename}", info.file_size, "not_dna_fasta"))
                            continue
                        stem = _stem_clean(info.filename)
                        if stem.endswith("_contigs"):
                            stem = stem[: -len("_contigs")]
                        gid = stem if stem not in contigs else f"DV_{stem}"
                        is_sag = stem.upper().startswith("AM-")
                        with zf.open(info) as raw, _open_maybe_compressed(raw, info.filename) as stream:
                            ingest_genome(
                                "DEATHVALLEY_SAG", "dv", gid, f"{p.name}:{info.filename}", stream,
                                genome_kind="SAG" if is_sag else "MAG",
                                quality="contig" if is_sag else "mag_mq",
                                tax="",
                                comp="",
                                cont="",
                                cluster_id="",
                                apply_single_dtr=False,
                            )

        while pending:
            s, fut = pending.pop(0)
            gc_results[s] = fut.result()[1]

        if not order_gids:
            raise Refused("BBSAG1: no qualifying genomes staged")

        missing_mq_hq = sorted(set(bb_kept) - {g for g in order_gids if ginfo[g]["subcohort"] == "BBSAG20"})
        if missing_mq_hq:
            raise Refused(
                f"BBSAG20: {len(missing_mq_hq)} MQ/HQ genomes in summary.tsv missing from archives (first: {missing_mq_hq[0]})"
            )

        kept_all = sorted(order_gids)
        resolved_codes, genome_uids, code_conflicts = c58._resolve_within_stage_code_conflicts(
            kept_all, contigs, gc_results
        )
        shared_circ_any: set[bytes] = {
            contigs[g][0][k * 16 : (k + 1) * 16]
            for g in kept_all
            for k, (_, _, _, is_c) in enumerate(spool_index[g])
            if is_c
        }

        part_writer.flush()
        own_all = _duck_ownership_from_parts(
            part_writer.parts, root, root_stages, owner_inputs, later_inputs, work_dir, threads, memory_limit
        )
        prior_owners = own_all["owners"]
        later_map = own_all["later"]

        # Lookup INSDC accessions if any exist
        accs_to_check = sorted({
            r232_index.nover(ginfo[g]["accession"])
            for g in kept_all
            if ginfo[g]["accession"]
        })
        look_map = {
            r232_index.nover(x["accession"]): x
            for x in (r232_index.lookup(index_dir, accs_to_check) if accs_to_check else [])
        }

        # Decide single stage (`BBSAG1`) vs deterministic two-stage split (`BBSAG1A`, `BBSAG1B`)
        total_new_bp = own_all["bases_new"]
        if total_new_bp <= max_stage_new_bp:
            stage_slices: list[tuple[str, list[str]]] = [(BBSAG1, kept_all)]
        else:
            # Keep DEATHVALLEY_SAG intact in BBSAG1B; split main genomes deterministically at midpoint <= max_stage_new_bp
            dv_gids = [g for g in kept_all if ginfo[g]["category"] == "dv"]
            main_gids = [g for g in kept_all if ginfo[g]["category"] != "dv"]
            seen_new_cu: set[bytes] = set()
            target_a = min(max_stage_new_bp, total_new_bp // 2)
            cum_a = 0
            gids_a: list[str] = []
            gids_b: list[str] = []
            for g in main_gids:
                raw, lens = contigs[g]
                uids = [raw[i : i + 16] for i in range(0, len(raw), 16)]
                add_bp = sum(n for cu, n in zip(uids, lens) if cu not in prior_owners and cu not in seen_new_cu)
                if gids_a and (cum_a + add_bp > target_a):
                    gids_b.append(g)
                else:
                    gids_a.append(g)
                    cum_a += add_bp
                    for cu in uids:
                        if cu not in prior_owners:
                            seen_new_cu.add(cu)
            gids_b.extend(dv_gids)
            gids_b.sort()
            if not gids_a or not gids_b:
                raise Refused(f"BBSAG1: cannot split {total_new_bp} net-new bp cleanly under cap {max_stage_new_bp}")
            stage_slices = [("BBSAG1A", gids_a), ("BBSAG1B", gids_b)]

        stage_reports: dict[str, dict] = {}
        extra_owners_from_earlier_slice: dict[bytes, str] = {}

        with open(spool_path, "rb") as spool_rd:
            for idx_slice, (st_name, slice_gids) in enumerate(stage_slices):
                fa_names = _fasta_names_for_bbsag_stage(st_name)
                sf = c78.StageFiles(st_name, Path(meta_out), Path(fasta_out), fa_names, work_dir / f"sf_{st_name}")
                if pool is not None:
                    for w in sf.fh.values():
                        if hasattr(w, "fh"):
                            w.fh.close()
                        else:
                            w.close()
                    sf.fh = {
                        k: (ParallelBgzfWriter(p, pool) if p.name.endswith(".gz") else open(p, "wb", buffering=1 << 22))
                        for k, p in sf.partial.items()
                    }
                try:
                    eff_owners = {**extra_owners_from_earlier_slice, **prior_owners}
                    listing = make_inputs.Listing()
                    genomes_out: list[dict] = []
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
                    used_cats: set[str] = set()
                    ordinal = 0

                    for g in slice_gids:
                        info = ginfo[g]
                        cat = info["category"]
                        used_cats.add(cat)
                        raw, lens = contigs[g]
                        uids = [raw[i : i + 16] for i in range(0, len(raw), 16)]
                        code = resolved_codes[g]
                        gc = gc_results[g]

                        # Copy formatted FASTA bytes from spool into StageFiles BGZF writer
                        for off, blen, _, _ in spool_index[g]:
                            spool_rd.seek(off)
                            sf.fh[cat].write(spool_rd.read(blen))

                        gc_rows.append((
                            g, st_name, info["taxonomy"], sum(lens), len(lens), info["quality"], code,
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
                            g, info["subcohort"], info["genome_kind"], acc, st_lbl,
                            "yes" if (x and x["in_gtdb_r232"]) else ("no" if acc else ""),
                            str(x["is_representative"]) if (x and x["in_gtdb_r232"]) else "",
                            (x["representative"] or "") if x else "",
                            (x["species"] or "") if x else "",
                            "yes" if (x and x["in_globdb_r232"]) else ("no" if acc else ""),
                            (x["globdb_r232_id"] or "") if x else "",
                            info["taxonomy"], info["completeness"], info["contamination"],
                            info["cluster_id"], len(lens), sum(lens),
                        ))

                        g_circ = [u for u in uids if u in shared_circ_any]
                        for k, (cu, n) in enumerate(zip(uids, lens)):
                            own = eff_owners.get(cu, "")
                            lat = "" if own else later_map.get(cu, "")
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
                        genomes_out.append(e)
                        first_of[genome_uids[g]] = g
                        by_subcohort[info["subcohort"]] += 1

                    sf.close()
                    slice_left = left_out if idx_slice == 0 else []
                    slice_skip = skipped_members if idx_slice == 0 else []
                    slice_trims = [t for t in trims if t[1] in set(slice_gids)]

                    make_inputs._write_tsv(sf.tmp / "records.tsv", BBSAG1_RECORDS_HEAD, rec_rows)
                    make_inputs._write_tsv(sf.tmp / "genetic_codes.tsv", c58.GC_HEAD, gc_rows)
                    make_inputs._write_tsv(sf.tmp / "r232_status.tsv", BBSAG1_STATUS_HEAD, status_rows)
                    make_inputs._write_tsv(sf.tmp / "left_out.tsv", BBSAG1_LEFT_HEAD, slice_left)
                    make_inputs._write_tsv(sf.tmp / "skipped_members.tsv", BBSAG1_SKIP_HEAD, slice_skip)
                    make_inputs._write_tsv(sf.tmp / "trims.tsv", c62.TRIM_HEAD, slice_trims)
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
                    (sf.tmp / "genomes.json").write_text(json.dumps(genomes_out))

                    used_ordered = [k for k in fa_names if k in used_cats]
                    source_meta = {
                        "main": {
                            "source_name": "tower_bbsag20_haima",
                            "source_release": "Hosokawa et al. human oral & gut SAG-gel microfluidic MQ/HQ single-cell & MAG catalog (figshare:24473008) and Haima cold-seep microfluidic SAGs & long/short-read MAGs (figshare:32084286); genetic codes 11/4/25 screened by c89_promge1_bbsag1_inputs.py (DEC-65, DEC-80 (3), DEC-84 (2), DEC-88)",
                            "license": "CC BY 4.0, MIT (confirmed in sources.yaml tower_bbsag20, tower_haima_lr_sag)",
                            "redistributable": True,
                            "retrieved_at": retrieved_main,
                        },
                        "dv": {
                            "source_name": "tower_deathvalley_sag",
                            "source_release": "Death Valley deep regional groundwater BLM1 MAGs & Bigelow SCGC SAGs (zenodo:21729218 BLM1-MAGs_SAGs_2011_2021.zip); genetic codes 11/4/25 screened by c89_promge1_bbsag1_inputs.py (DEC-80 (2)-(3), DEC-88)",
                            "license": "none stated: zenodo record 21729218 licence field is empty (sources.yaml tower_deathvalley_sag, DEC-80 (2))",
                            "redistributable": False,
                            "retrieved_at": retrieved_dv,
                        },
                    }
                    _write_objects_json(sf, used_ordered, source_meta, Path(inputs_home))
                    sha_map, byte_map, locations = _commit_with_stor_floor(sf, used_ordered)

                    held_c: Counter[str] = Counter()
                    held_b: Counter[str] = Counter()
                    lat_c: Counter[str] = Counter()
                    lat_b: Counter[str] = Counter()
                    for cu, n in contig_len.items():
                        if cu in eff_owners:
                            held_c[eff_owners[cu]] += 1
                            held_b[eff_owners[cu]] += n
                        elif cu in later_map:
                            lat_c[later_map[cu]] += 1
                            lat_b[later_map[cu]] += n

                    for cu in contig_len:
                        if cu not in prior_owners and cu not in extra_owners_from_earlier_slice:
                            extra_owners_from_earlier_slice[cu] = st_name

                    stage_reports[st_name] = {
                        "stage": st_name,
                        "genomes_input": len(slice_gids),
                        "genomes": len(genomes_out),
                        "repeated_genomes": len(repeated),
                        "circular_genomes": sum(1 for g in genomes_out if g["circular"]),
                        "genomes_by_subcohort": dict(sorted(by_subcohort.items())),
                        "genetic_codes": dict(sorted(Counter(str(resolved_codes[g]) for g in slice_gids).items())),
                        "genetic_codes_by_subcohort": {
                            sub: dict(sorted(Counter(str(resolved_codes[g]) for g in slice_gids if ginfo[g]["subcohort"] == sub).items()))
                            for sub in sorted({ginfo[g]["subcohort"] for g in slice_gids})
                        },
                        "r232_status": dict(sorted(status_counts.items())),
                        "left_out": len(slice_left),
                        "skipped_members": len(slice_skip),
                        "trims": len(slice_trims),
                        "contigs_total": sum(len(contigs[g][1]) for g in slice_gids),
                        "contigs_distinct": len(contig_len),
                        "bases_distinct": sum(contig_len.values()),
                        "held_by": dict(sorted(held_c.items())),
                        "held_bases_by": dict(sorted(held_b.items())),
                        "later_held_by": dict(sorted(lat_c.items())),
                        "later_held_bases_by": dict(sorted(lat_b.items())),
                        "contigs_new": sum(1 for cu in contig_len if cu not in eff_owners),
                        "bases_new": sum(n for cu, n in contig_len.items() if cu not in eff_owners),
                        "sha256": sha_map,
                        "bytes": byte_map,
                        "locations": locations,
                        **{"listing_" + k: v for k, v in listing.report().items()},
                    }
                except BaseException:
                    sf.abort()
                    raise

        return {
            "stages_written": [st for st, _ in stage_slices],
            "total_genomes_input": len(kept_all),
            "total_contigs_distinct": own_all["contigs_distinct"],
            "total_bases_distinct": own_all["bases_distinct"],
            "total_bases_new": total_new_bp,
            "code_conflicts": len(code_conflicts),
            "bbsag20_lq_excluded": len([x for x in left_out if x[6] == "bbsag20_low_quality"]),
            "dtr_classes": dict(sorted(dtr_classes.items())),
            "stages": stage_reports,
        }
    finally:
        spool_path.unlink(missing_ok=True)
        shutil.rmtree(work_dir / "u_parts", ignore_errors=True)
        if own_work:
            shutil.rmtree(work_dir, ignore_errors=True)
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)


# ==================================================================================================================
# 6. `skani` >= 96% ANI / >= 50% AF check (recorded, none cut per DEC-84 (2))
# ==================================================================================================================


def _prepare_meta_with_records(meta: Path, work: Path) -> Path:
    """If `records.tsv` was moved to `/stor` (`>= 64 MiB`), link it into a small work meta view for `c81.split_genomes`."""
    meta = Path(meta)
    if (meta / "records.tsv").is_file() or (meta / "records.tsv.gz").is_file():
        return meta
    objs = json.loads((meta / "objects_staging.json").read_text())
    for o in objs:
        cand = Path(o["path"]).parent / "records.tsv"
        if cand.is_file():
            view = work / f"meta_view_{meta.name}"
            view.mkdir(parents=True, exist_ok=True)
            (view / "records.tsv").symlink_to(cand)
            return view
    raise Refused(f"no records.tsv found in {meta} or its staging directory")


def bbsag1_ani(
    meta: Path,
    fasta_paths: list[Path],
    out: Path,
    *,
    whole: dict[str, tuple[Path, list[Path]]],
    by_genus: dict[str, tuple[Path, list[Path]]] | None = None,
    taxonomy_tsv: Path | None = None,
    archaea: tuple[Path, Path] | None = None,
    threads: int = 4,
) -> dict:
    """Record `skani` >= 96% ANI / >= 50% larger AF matches against held stages without removing any genome (DEC-84 (2))."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    work = out / "work"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True, exist_ok=True)

    q_meta = _prepare_meta_with_records(Path(meta), work)
    q_files = c81.split_genomes(q_meta, [Path(x) for x in fasta_paths], work / "Q")
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
        d = work / f"ref_{name}"
        unlisted: list[str] = []
        if kind == "archaea":
            files = c81.archaea_refs(v[0], v[1], genera, d)
        else:
            r_meta = _prepare_meta_with_records(v[0], work) if not (Path(v[0]) / "records.tsv").is_file() else v[0]
            files = c81.split_genomes(r_meta, v[1], d, keep=keep if kind == "genus" else None, unlisted=unlisted)
        ref_of = {str(p): (name, g) for g, p in files.items()}
        entry: dict = {
            "kind": kind,
            "references": len(files),
            "records_in_no_genome": len(unlisted),
        }
        if files:
            rlist = work / f"r_{name}.txt"
            rlist.write_text("".join(f"{p}\n" for p in sorted(ref_of)))
            table = work / f"vs_{name}.tsv"
            t0 = c39_ani.child_cpu()
            subprocess.run(
                [
                    "skani", "dist", "-t", str(threads), "--min-af", str(c81.AF_SAG),
                    "--ql", str(work / "q.txt"), "--rl", str(rlist), "-o", str(table),
                ],
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
        out / "ani_hits.tsv",
        ["genome", "stage", "ref_genome", "ani", "af"],
        [(g, h["stage"], h["genome"], h["ani"], h["af"]) for g, h in sorted(best.items())],
    )
    res["meet"] = len(best)
    res["meet_by_stage"] = dict(sorted(Counter(h["stage"] for h in best.values()).items()))
    (out / "ani_summary.json").write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
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
    p.add_argument("--with-c86-defaults", action="store_true",
                   help="include merged C-86 r26.11.6 stages (O4, GV2, HROM1) in owner_inputs")
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
    if a.with_c86_defaults:
        for s in C86_STAGES:
            owners[s] = Path(DEFAULT_C86_DIRS[s])
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
    ap = argparse.ArgumentParser(description="Card C-89: stage PROMGE1 and BBSAG1 cohorts")
    sub = ap.add_subparsers(dest="cmd", required=True)

    pc = sub.add_parser("census")
    pc.add_argument("paths", nargs="+", type=Path)
    pc.add_argument("--out", type=Path, default=None)

    ppm = sub.add_parser("promge1-stage")
    ppm.add_argument("--dir", action="append", default=[], metavar="SUB=DIR",
                     help="may be repeated per subcohort (e.g. both /stor/.../dl/SUB and /tmp/.../dl/SUB)")
    ppm.add_argument("--index-dir", type=Path, required=True)
    ppm.add_argument("--meta-out", type=Path, required=True)
    ppm.add_argument("--fasta-out", type=Path, required=True)
    ppm.add_argument("--exclude-promge2-cellular", action="store_true")
    ppm.add_argument("--promge2-non-nested-only", action="store_true")
    ppm.add_argument("--promge2-spire-and-pg3-reps-only", action="store_true")
    ppm.add_argument("--retrieved-promge2", default=PROMGE1_RETRIEVED_PROMGE2)
    ppm.add_argument("--retrieved-open", default=PROMGE1_RETRIEVED_OPEN)
    _add_common_stage_args(ppm)

    pbb = sub.add_parser("bbsag1-stage")
    pbb.add_argument("--bbsag20-dir", action="append", required=True, type=Path)
    pbb.add_argument("--haima-dir", action="append", required=True, type=Path)
    pbb.add_argument("--deathvalley-dir", action="append", required=True, type=Path)
    pbb.add_argument("--index-dir", type=Path, required=True)
    pbb.add_argument("--meta-out", type=Path, required=True)
    pbb.add_argument("--fasta-out", type=Path, required=True)
    pbb.add_argument("--max-stage-new-bp", type=int, default=BBSAG1_MAX_NEW_BP)
    pbb.add_argument("--retrieved-main", default=BBSAG1_RETRIEVED_MAIN)
    pbb.add_argument("--retrieved-dv", default=BBSAG1_RETRIEVED_DV)
    _add_common_stage_args(pbb)

    pba = sub.add_parser("bbsag1-ani")
    pba.add_argument("meta", type=Path)
    pba.add_argument("out", type=Path)
    pba.add_argument("--fasta", action="append", required=True, type=Path)
    pba.add_argument("--whole", action="append", default=[], metavar="NAME=META_DIR:FASTA[,FASTA]")
    pba.add_argument("--by-genus", action="append", default=[], metavar="NAME=META_DIR:FASTA[,FASTA]")
    pba.add_argument("--taxonomy", type=Path, default=None)
    pba.add_argument("--archaea", default="", metavar="FASTA:LIST_TSV")
    pba.add_argument("--threads", type=int, default=4)
    pba.add_argument("--report-out", type=Path, default=None)

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

        if a.cmd == "bbsag1-ani":
            arc = None
            if a.archaea:
                fa, sep, lst = a.archaea.partition(":")
                if not sep or not fa or not lst:
                    raise Refused(f"--archaea {a.archaea}: expected FASTA:LIST_TSV")
                arc = (Path(fa), Path(lst))
            rep = bbsag1_ani(
                a.meta, a.fasta, a.out,
                whole=_parse_group_sets(a.whole, "--whole"),
                by_genus=_parse_group_sets(a.by_genus, "--by-genus"),
                taxonomy_tsv=a.taxonomy,
                archaea=arc,
                threads=a.threads,
            )
            txt = json.dumps(rep, indent=1, sort_keys=True)
            if a.report_out is not None:
                a.report_out.parent.mkdir(parents=True, exist_ok=True)
                a.report_out.write_text(txt + "\n")
            print(txt)
            return 0

        root_stages, owners, later = _resolve_owners_and_later(a)
        if a.cmd == "promge1-stage":
            pm_dirs: dict[str, list[Path]] = {}
            for item in a.dir:
                k, sep, v = item.partition("=")
                if not sep or not k or not v:
                    raise Refused(f"--dir {item!r}: expected SUB=DIR")
                pm_dirs.setdefault(k, []).append(Path(v))
            rep = stage_promge1(
                pm_dirs, a.index_dir, a.meta_out, a.fasta_out,
                root=a.root, root_stages=root_stages, owner_inputs=owners, later_inputs=later,
                inputs_home=a.inputs_home,
                exclude_promge2_cellular=a.exclude_promge2_cellular,
                promge2_non_nested_only=a.promge2_non_nested_only,
                promge2_spire_and_pg3_reps_only=a.promge2_spire_and_pg3_reps_only,
                retrieved_promge2=a.retrieved_promge2,
                retrieved_open=a.retrieved_open,
                work=a.work, threads=a.threads, memory_limit=a.memory_limit,
            )
            rep["owner_stages"] = [*(root_stages if a.root else []), *owners]
            rep["later_stages"] = [*later]
            rep["s2_check"] = c81._check(PROMGE1, a.meta_out / PROMGE1)
            ok = rep["s2_check"]["exit"] == 0 and rep["s2_check"]["problems"] == []
        else:
            rep = stage_bbsag1(
                a.bbsag20_dir, a.haima_dir, a.deathvalley_dir, a.index_dir,
                a.meta_out, a.fasta_out,
                root=a.root, root_stages=root_stages, owner_inputs=owners, later_inputs=later,
                inputs_home=a.inputs_home,
                max_stage_new_bp=a.max_stage_new_bp,
                retrieved_main=a.retrieved_main,
                retrieved_dv=a.retrieved_dv,
                work=a.work, threads=a.threads, memory_limit=a.memory_limit, workers=a.workers,
            )
            rep["owner_stages"] = [*(root_stages if a.root else []), *owners]
            rep["later_stages"] = [*later]
            ok = True
            for st_name in rep["stages_written"]:
                chk = c81._check(st_name, a.meta_out / st_name)
                rep["stages"][st_name]["s2_check"] = chk
                if chk["exit"] != 0 or chk["problems"] != []:
                    ok = False

        txt = json.dumps(rep, indent=1, sort_keys=True)
        if a.report_out is not None:
            a.report_out.parent.mkdir(parents=True, exist_ok=True)
            a.report_out.write_text(txt + "\n")
        print(txt)
        return 0 if ok else 1
    except Refused as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
