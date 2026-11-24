"""Card C-59 (DEC-69, DEC-65, DEC-50, DEC-26, DEC-27):
Audit, benchmark, deduplicate, and stage NCBI RefSeq Organellar genomes (`O1`: plastids,
chromatophores, and table-11 mitochondria) from RefSeq release `.gbff.gz` flatfiles, and
record the non-table-11 mitochondrial/organellar remainder in `O1_left_for_c59a.tsv.gz`.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
import re
import shutil
import sys
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # c48c_virus_inputs and make_inputs beside this script

import c48c_virus_inputs as v2  # noqa: E402
import make_inputs  # noqa: E402
from biobox import fasta, ids  # noqa: E402

csv.field_size_limit(1 << 30)
Refused = make_inputs.Refused

VALID_DNA_RE = re.compile(rb"^[ACGTURYKMSWBDHVNacgturykmswbdhvn]+$")
VALID_KINDS = ("plasmid", "prokaryote", "virus")

HOMING_RT_RE = re.compile(
    r"(?:homing\s+endonuclease|laglidadg|giy-?yig|\bhnh\b|maturase|\bmatk\b|\bmatr\b|"
    r"group\s*(?:i|ii|1|2)\s*intron|reverse\s+transcriptase|intron-encoded|intron\s+orf|"
    r"\bi-(?:cre|ceu|chl|sce|ani|dir|ppo|tevi?)\b|\brtl\b)",
    re.IGNORECASE,
)

ORIGIN_WRAP_RE = re.compile(
    r"^(?:complement\()?join\((\d+)\.\.(\d+),(\d+)\.\.(\d+)\)\)?$"
)
_ORIGIN_TRANS = str.maketrans("", "", " 0123456789\t\r\n")

O1_RECORD_HEAD = [
    "ordinal", "source", "compartment", "accession", "organism", "contig_uid",
    "length", "circular", "transl_table", "n_cds", "n_spliced_cds", "n_introns",
    "n_homing_rt", "owner_stage",
]
O1_META_HEAD = [
    "stage", "source", "compartment", "genome_id", "contig_uid", "accession",
    "organism", "length", "circular", "transl_table", "n_cds", "n_spliced_cds",
    "n_introns", "n_homing_rt", "owner_stage",
]
O1_LEFT_HEAD = [
    "accession", "source", "compartment", "organism", "taxonomy", "length",
    "circular", "transl_table", "transl_tables", "n_cds", "n_spliced_cds",
    "n_introns", "n_homing_rt", "reason", "needed_capability",
]


def _base_acc(acc: str) -> str:
    acc = acc.strip()
    return acc.split(".", 1)[0] if "." in acc else acc


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _is_spliced_location(loc: str, seq_len: int, circular: bool) -> bool:
    """Return True if a GenBank CDS location represents multi-exon splicing (not origin wrap)."""
    clean = re.sub(r"\s+", "", loc)
    if "join(" not in clean:
        return False
    if circular and seq_len > 0:
        m = ORIGIN_WRAP_RE.match(clean)
        if m:
            s1, e1, s2, e2 = int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4))
            if (e1 == seq_len and s2 == 1) or (e2 == seq_len and s1 == 1):
                return False
    return True


def _is_homing_or_rt(gene: str, product: str) -> bool:
    text = f"{gene} {product}".strip()
    return bool(text and HOMING_RT_RE.search(text))


def _modal_transl_table(tables: list[int]) -> int | None:
    if not tables:
        return None
    counts = Counter(tables)
    # Highest count, tie-break by smallest translation table
    return min(counts.keys(), key=lambda t: (-counts[t], t))


def parse_refseq_organelle_gbff(gb_fh, *, keep_translations: bool = True) -> Iterator[dict]:
    """Stream RefSeq GenBank flatfile records, extracting organellar metadata, CDSs, introns, and sequence."""
    acc = ""
    ver = ""
    moltype = ""
    circular = False
    defn_lines: list[str] = []
    org = ""
    tax_lines: list[str] = []
    organelle = ""
    cds_list: list[dict] = []
    n_introns = 0
    seq_chunks: list[str] = []
    state = "head"

    # Current feature state
    feat_type = ""
    feat_loc_parts: list[str] = []
    feat_qual_key = ""
    feat_qual_val_parts: list[str] = []
    cur_cds: dict | None = None

    def _finish_qualifier():
        nonlocal feat_qual_key, feat_qual_val_parts, organelle
        if not feat_qual_key:
            return
        val = "".join(feat_qual_val_parts) if feat_qual_key == "translation" else " ".join(feat_qual_val_parts)
        val = val.strip()
        if val.startswith('"') and val.endswith('"') and len(val) >= 2:
            val = val[1:-1]
        if feat_type == "source" and feat_qual_key == "organelle" and not organelle:
            organelle = val
        elif cur_cds is not None:
            if feat_qual_key == "transl_table":
                try:
                    cur_cds["transl_table"] = int(val)
                except ValueError:
                    pass
            elif feat_qual_key == "product":
                cur_cds["product"] = val
            elif feat_qual_key == "gene":
                cur_cds["gene"] = val
            elif feat_qual_key == "translation" and keep_translations:
                cur_cds["translation"] = re.sub(r"\s+", "", val)
        feat_qual_key = ""
        feat_qual_val_parts = []

    def _finish_feature():
        nonlocal feat_type, feat_loc_parts, cur_cds
        _finish_qualifier()
        if cur_cds is not None:
            cur_cds["location"] = "".join(feat_loc_parts)
            cds_list.append(cur_cds)
            cur_cds = None
        feat_type = ""
        feat_loc_parts = []

    def _emit():
        nonlocal acc, ver, moltype, circular, defn_lines, org, tax_lines, organelle
        nonlocal cds_list, n_introns, seq_chunks, state
        _finish_feature()
        if not acc:
            return None
        seq_str = "".join(seq_chunks).upper()
        seq_len = len(seq_str)
        tables: list[int] = []
        n_spliced = 0
        n_homing = 0
        for c in cds_list:
            c["is_spliced"] = _is_spliced_location(c["location"], seq_len, circular)
            c["is_homing_rt"] = _is_homing_or_rt(c.get("gene", ""), c.get("product", ""))
            if c["is_spliced"]:
                n_spliced += 1
            if c["is_homing_rt"]:
                n_homing += 1
            # INSDC Feature Table rule: CDS without /transl_table= uses standard genetic code (1)
            if c.get("transl_table") is None:
                c["transl_table"] = 1
            tables.append(c["transl_table"])
        rec = {
            "accession": ver or acc,
            "base_accession": _base_acc(ver or acc),
            "moltype": moltype,
            "circular": circular,
            "length": seq_len,
            "definition": " ".join(defn_lines),
            "organism": org,
            "taxonomy": " ".join(tax_lines),
            "organelle": organelle,
            "transl_table": _modal_transl_table(tables),
            "transl_tables": sorted(set(tables)),
            "n_cds": len(cds_list),
            "n_spliced_cds": n_spliced,
            "n_introns": n_introns,
            "n_homing_rt": n_homing,
            "cds": cds_list,
            "seq": seq_str.encode("ascii"),
        }
        acc = ""
        ver = ""
        moltype = ""
        circular = False
        defn_lines = []
        org = ""
        tax_lines = []
        organelle = ""
        cds_list = []
        n_introns = 0
        seq_chunks = []
        state = "head"
        return rec

    for raw_line in gb_fh:
        line = raw_line.rstrip("\r\n")
        if line.startswith("//"):
            rec = _emit()
            if rec is not None:
                yield rec
            continue
        if state == "origin":
            seq_chunks.append(line.translate(_ORIGIN_TRANS))
            continue
        if line.startswith("LOCUS"):
            parts = line.split()
            if len(parts) >= 2:
                acc = parts[1]
            low = [p.lower() for p in parts]
            circular = "circular" in low
            for tok in parts[2:]:
                u = tok.upper()
                if "DNA" in u or "RNA" in u:
                    moltype = tok
                    break
            state = "head"
        elif line.startswith("DEFINITION"):
            defn_lines = [line[12:].strip()]
            state = "def"
        elif line.startswith("VERSION"):
            parts = line.split()
            if len(parts) >= 2:
                ver = parts[1]
            state = "head"
        elif line.startswith("  ORGANISM"):
            org = line[12:].strip()
            tax_lines = []
            state = "tax"
        elif line.startswith("FEATURES"):
            state = "features"
        elif line.startswith("ORIGIN"):
            _finish_feature()
            state = "origin"
        elif state == "def":
            if line.startswith("            "):
                defn_lines.append(line.strip())
            else:
                state = "head"
        elif state == "tax":
            if line.startswith("            "):
                tax_lines.append(line.strip())
            else:
                state = "head"
        elif state == "features":
            if not line.startswith(" "):
                _finish_feature()
                state = "head"
                continue
            # New feature key at column 5 (indent 5 spaces, non-space at index 5)
            if len(line) > 5 and line[:5] == "     " and line[5] != " ":
                _finish_feature()
                parts = line[5:].split(None, 1)
                feat_type = parts[0]
                feat_loc_parts = [parts[1].strip()] if len(parts) > 1 else []
                if feat_type == "intron":
                    n_introns += 1
                elif feat_type == "CDS":
                    cur_cds = {
                        "location": "",
                        "gene": "",
                        "product": "",
                        "transl_table": None,
                        "translation": "",
                    }
            elif line.startswith("                     "):
                s = line[21:]
                if s.startswith("/"):
                    _finish_qualifier()
                    if "=" in s:
                        k, v = s[1:].split("=", 1)
                        feat_qual_key = k
                        feat_qual_val_parts = [v]
                    else:
                        feat_qual_key = s[1:]
                        feat_qual_val_parts = ["true"]
                elif feat_qual_key:
                    feat_qual_val_parts.append(s)
                else:
                    # Location continuation line
                    feat_loc_parts.append(s.strip())


def normalize_compartment(organelle: str, definition: str, source_compartment: str) -> str:
    org_l = (organelle or "").lower()
    def_l = (definition or "").lower()
    if "chromatophore" in org_l or "chromatophore" in def_l:
        return "chromatophore"
    if "nucleomorph" in org_l or "nucleomorph" in def_l:
        return "nucleomorph"
    if "plastid" in org_l or "chloroplast" in org_l or "apicoplast" in org_l or "cyanelle" in org_l:
        return "plastid"
    if "mitochondri" in org_l or "kinetoplast" in org_l:
        return "mitochondrion"
    return source_compartment


def partition_organelle_record(rec: dict, *, source_compartment: str) -> dict:
    """Classify one RefSeq organellar record into O1 (table 11) or left_for_c59a."""
    comp = normalize_compartment(rec.get("organelle", ""), rec.get("definition", ""), source_compartment)
    seq = rec.get("seq") or b""
    tt = rec.get("transl_table")
    n_spliced = int(rec.get("n_spliced_cds") or 0)

    if not seq or not VALID_DNA_RE.match(seq):
        return {
            "bucket": "left_for_c59a",
            "compartment": comp,
            "reason": "invalid_or_empty_seq",
            "needed_capability": "manual_review",
        }
    if tt == 11:
        return {
            "bucket": "O1",
            "compartment": comp,
            "reason": "table11_o1",
            "needed_capability": "",
        }
    if tt is None:
        return {
            "bucket": "left_for_c59a",
            "compartment": comp,
            "reason": "no_transl_table",
            "needed_capability": "c59a_transl_table_or_unannotated",
        }
    cap = "c59a_transl_table_and_c60_spliced" if n_spliced > 0 else "c59a_transl_table"
    return {
        "bucket": "left_for_c59a",
        "compartment": comp,
        "reason": f"transl_table_{tt}",
        "needed_capability": cap,
    }


# --------------------------------------------------------------------------------------------------------------
# Pyrodigal recall benchmark against RefSeq annotated CDSs and homing-endonuclease / group-II-intron-RT ORFs


def _aa_match(ref_aa: str, pred_aas: set[str], pred_list: list[str]) -> tuple[bool, bool]:
    """Return (exact_match, overlap_match) of a RefSeq CDS protein against predicted proteins.

    Exact match ignores leading alternative start codon (M vs V/L/I) if full length matches.
    Overlap match checks whether >= 80% of the shorter sequence is contained as a contiguous
    substring (catching start-codon shifts or intron-split exons).
    """
    if not ref_aa:
        return False, False
    if ref_aa in pred_aas:
        return True, True
    ref_tail = ref_aa[1:]
    rlen = len(ref_aa)
    for p in pred_list:
        if len(p) == rlen and p[1:] == ref_tail:
            return True, True
    # Substring / 80% containment check
    min_k = max(15, int(0.8 * rlen))
    for p in pred_list:
        if rlen >= 15 and p in ref_aa and len(p) >= min_k:
            return False, True
        if len(p) >= 15 and ref_aa in p:
            return False, True
        # Shared 20-aa anchor near middle -> check overlap length
        if rlen >= 25 and len(p) >= 25:
            mid = ref_aa[rlen // 3 : rlen // 3 + 20]
            if mid in p:
                return False, True
    return False, False


def _call_proteins_mode(seq_str: str, circular: bool, mode: str, transl_table: int) -> list[str]:
    import pyrodigal
    import pyrodigal_gv
    from biobox.l1.genecall import SINGLE_MODE_MIN_BP, _genes

    if mode == "meta":
        finder = pyrodigal.GeneFinder(meta=True)
    elif mode == "gv-meta":
        finder = pyrodigal_gv.ViralGeneFinder(meta=True)
    elif mode == "single":
        if len(seq_str) < 20_000:
            finder = pyrodigal.GeneFinder(meta=True)
        else:
            finder = pyrodigal.GeneFinder(meta=False)
            finder.train(seq_str, translation_table=transl_table)
    elif mode == "prokaryote_policy":
        if len(seq_str) < SINGLE_MODE_MIN_BP:
            finder = pyrodigal.GeneFinder(meta=True)
        else:
            finder = pyrodigal.GeneFinder(meta=False)
            finder.train(seq_str, translation_table=transl_table)
    else:
        raise ValueError(f"unknown mode {mode!r}")

    return [aa for _, _, _, aa, _, _ in _genes(finder, seq_str, circular)]


def benchmark_pyrodigal_organelles(
    records: list[dict],
    *,
    modes: tuple[str, ...] = ("meta", "single", "prokaryote_policy", "gv-meta"),
) -> dict:
    """Benchmark pyrodigal recall against RefSeq CDS and homing-endonuclease/RT annotations."""
    stats: dict[str, dict[str, dict[str, int]]] = {}
    for rec in records:
        seq_str = rec["seq"].decode("ascii")
        circular = bool(rec["circular"])
        tt = rec.get("transl_table") or 11
        comp = rec.get("compartment") or "unknown"
        stratum = f"{comp}:table_{tt}"
        cds_with_aa = [c for c in rec.get("cds", []) if c.get("translation")]
        if not cds_with_aa:
            continue
        stratum_d = stats.setdefault(stratum, {})
        for mode in modes:
            # Pyrodigal single-mode supports translation tables 1..25 (except 7,8,17,18,19,20)
            try:
                preds = _call_proteins_mode(seq_str, circular, mode, tt)
            except Exception:
                continue
            pred_set = set(preds)
            m = stratum_d.setdefault(mode, {
                "genomes": 0,
                "pred_orfs": 0,
                "cds_total": 0,
                "cds_exact": 0,
                "cds_overlap": 0,
                "unspliced_total": 0,
                "unspliced_exact": 0,
                "unspliced_overlap": 0,
                "spliced_total": 0,
                "spliced_overlap": 0,
                "homing_rt_total": 0,
                "homing_rt_exact": 0,
                "homing_rt_overlap": 0,
            })
            m["genomes"] += 1
            m["pred_orfs"] += len(preds)
            for c in cds_with_aa:
                exact, ov = _aa_match(c["translation"], pred_set, preds)
                m["cds_total"] += 1
                if exact:
                    m["cds_exact"] += 1
                if ov:
                    m["cds_overlap"] += 1
                if c.get("is_spliced"):
                    m["spliced_total"] += 1
                    if ov:
                        m["spliced_overlap"] += 1
                else:
                    m["unspliced_total"] += 1
                    if exact:
                        m["unspliced_exact"] += 1
                    if ov:
                        m["unspliced_overlap"] += 1
                if c.get("is_homing_rt"):
                    m["homing_rt_total"] += 1
                    if exact:
                        m["homing_rt_exact"] += 1
                    if ov:
                        m["homing_rt_overlap"] += 1

    # Add recall percentages
    out: dict[str, dict[str, dict]] = {}
    for stratum, mode_map in sorted(stats.items()):
        out[stratum] = {}
        for mode, m in sorted(mode_map.items()):
            rec_row = dict(m)
            rec_row["cds_exact_pct"] = round(100.0 * m["cds_exact"] / max(1, m["cds_total"]), 2)
            rec_row["cds_overlap_pct"] = round(100.0 * m["cds_overlap"] / max(1, m["cds_total"]), 2)
            rec_row["unspliced_exact_pct"] = round(100.0 * m["unspliced_exact"] / max(1, m["unspliced_total"]), 2)
            rec_row["unspliced_overlap_pct"] = round(100.0 * m["unspliced_overlap"] / max(1, m["unspliced_total"]), 2)
            rec_row["spliced_overlap_pct"] = round(100.0 * m["spliced_overlap"] / max(1, m["spliced_total"]), 2) if m["spliced_total"] else None
            rec_row["homing_rt_exact_pct"] = round(100.0 * m["homing_rt_exact"] / max(1, m["homing_rt_total"]), 2) if m["homing_rt_total"] else None
            rec_row["homing_rt_overlap_pct"] = round(100.0 * m["homing_rt_overlap"] / max(1, m["homing_rt_total"]), 2) if m["homing_rt_total"] else None
            out[stratum][mode] = rec_row
    return out


def audit_gbff_files(
    sources_spec: list[dict],
    *,
    max_per_stratum: int = 30,
    seed: int = 59,
) -> dict:
    """Stream all RefSeq organelle .gbff.gz files, tally compartments/tables/introns, and benchmark pyrodigal."""
    rng = random.Random(seed)
    by_stratum: dict[str, dict[str, int]] = {}
    reservoir: dict[str, list[dict]] = defaultdict(list)
    seen_stratum: Counter[str] = Counter()
    total_records = 0
    total_bases = 0

    for spec in sources_spec:
        src_path = Path(spec["path"])
        src_comp = spec["source_compartment"]
        with _open_text(src_path) as fh:
            for rec in parse_refseq_organelle_gbff(fh, keep_translations=True):
                part = partition_organelle_record(rec, source_compartment=src_comp)
                comp = part["compartment"]
                tt = rec["transl_table"]
                stratum = f"{comp}:table_{tt}"
                total_records += 1
                total_bases += rec["length"]

                st = by_stratum.setdefault(stratum, {
                    "records": 0,
                    "bases": 0,
                    "circular": 0,
                    "ge_100kb": 0,
                    "n_cds": 0,
                    "n_spliced_cds": 0,
                    "n_introns": 0,
                    "n_homing_rt": 0,
                    "records_with_introns": 0,
                    "records_with_homing_rt": 0,
                })
                st["records"] += 1
                st["bases"] += rec["length"]
                if rec["circular"]:
                    st["circular"] += 1
                if rec["length"] >= 100_000:
                    st["ge_100kb"] += 1
                st["n_cds"] += rec["n_cds"]
                st["n_spliced_cds"] += rec["n_spliced_cds"]
                st["n_introns"] += rec["n_introns"]
                st["n_homing_rt"] += rec["n_homing_rt"]
                if rec["n_introns"] > 0 or rec["n_spliced_cds"] > 0:
                    st["records_with_introns"] += 1
                if rec["n_homing_rt"] > 0:
                    st["records_with_homing_rt"] += 1

                # Reservoir sample per stratum for pyrodigal benchmark (bounded memory)
                if rec["seq"] and VALID_DNA_RE.match(rec["seq"]) and rec["n_cds"] > 0:
                    seen_stratum[stratum] += 1
                    k = seen_stratum[stratum]
                    sample_rec = {
                        "accession": rec["accession"],
                        "compartment": comp,
                        "circular": rec["circular"],
                        "transl_table": tt,
                        "cds": rec["cds"],
                        "seq": rec["seq"],
                    }
                    if len(reservoir[stratum]) < max_per_stratum:
                        reservoir[stratum].append(sample_rec)
                    else:
                        j = rng.randrange(k)
                        if j < max_per_stratum:
                            reservoir[stratum][j] = sample_rec

    sampled_records = [r for st in sorted(reservoir) for r in reservoir[st]]
    bench = benchmark_pyrodigal_organelles(sampled_records)
    return {
        "total_records": total_records,
        "total_bases": total_bases,
        "strata": dict(sorted(by_stratum.items())),
        "benchmark": bench,
    }


# --------------------------------------------------------------------------------------------------------------
# Stage O1 (RefSeq Organelles transl_table=11; DEC-65 (2) byte-identical as deposited; redistributable=True)


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


def stage_o1(
    sources_spec: list[dict],
    out: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
    genome_kind: str = "plasmid",
    compress_fasta: bool = False,
    final_dir: Path | None = None,
) -> dict:
    """Stage `O1` inputs from RefSeq organellar `.gbff.gz` files (transl_table=11)."""
    if genome_kind not in VALID_KINDS:
        raise Refused(f"unsupported genome_kind {genome_kind!r}; expected one of {VALID_KINDS}")
    out = Path(out)
    if out.exists():
        raise Refused(f"{out} exists")
    if not sources_spec:
        raise Refused("sources_spec is empty")
    for spec in sources_spec:
        p = Path(spec["path"])
        if not p.is_file():
            raise Refused(f"missing input gbff file: {p}")

    work = out.with_name(f"{out.name}.tmp")
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    obj_dir = Path(final_dir) if final_dir is not None else out

    try:
        rep = _stage_o1(
            sources_spec,
            obj_dir,
            work,
            root,
            list(root_stages),
            owner_inputs,
            retrieved,
            genome_kind,
            compress_fasta,
        )
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise
    work.rename(out)
    rep["sha256"] = {
        p.name: ids.sha256_file(p).hex()
        for p in sorted(out.iterdir())
        if p.is_file()
    }
    rep["bytes"] = {
        p.name: p.stat().st_size
        for p in sorted(out.iterdir())
        if p.is_file()
    }
    return rep


def _open_fasta_out(path: Path, compress: bool):
    if compress:
        raw = open(path, "wb")
        gz = gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=1, mtime=0)
        return raw, gz
    f = open(path, "wb")
    return None, f


def _stage_o1(
    sources_spec: list[dict],
    obj_dir: Path,
    work: Path,
    root: Path | None,
    root_stages: list[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
    genome_kind: str,
    compress_fasta: bool,
) -> dict:
    # Group input files by source_name so plastid.1..3 land in one plastid object FASTA
    by_source: dict[str, list[dict]] = defaultdict(list)
    source_order: list[str] = []
    for spec in sources_spec:
        sn = spec["source_name"]
        if sn not in by_source:
            source_order.append(sn)
        by_source[sn].append(spec)

    o1_records: list[tuple] = []  # (src, comp, acc, org, cu, length, circ, tt, n_cds, n_spliced, n_introns, n_homing)
    left_rows: list[tuple] = []
    objects: list[dict] = []
    by_comp_o1: Counter[str] = Counter()
    by_reason_left: Counter[str] = Counter()
    ext = ".fna.gz" if compress_fasta else ".fna"

    for sn in source_order:
        specs = by_source[sn]
        src_rel = specs[0].get("source_release", "238")
        out_name = f"{sn}_t11{ext}"
        fa_dst = work / out_name
        written = 0
        raw_fh, w = _open_fasta_out(fa_dst, compress_fasta)
        try:
            for spec in specs:
                src_path = Path(spec["path"])
                src_comp = spec["source_compartment"]
                with _open_text(src_path) as fh:
                    for rec in parse_refseq_organelle_gbff(fh, keep_translations=False):
                        part = partition_organelle_record(rec, source_compartment=src_comp)
                        comp = part["compartment"]
                        if part["bucket"] != "O1":
                            by_reason_left[part["reason"]] += 1
                            left_rows.append((
                                rec["accession"],
                                sn,
                                comp,
                                rec["organism"],
                                rec["taxonomy"],
                                rec["length"],
                                str(rec["circular"]).lower(),
                                "" if rec["transl_table"] is None else rec["transl_table"],
                                ",".join(str(t) for t in rec["transl_tables"]),
                                rec["n_cds"],
                                rec["n_spliced_cds"],
                                rec["n_introns"],
                                rec["n_homing_rt"],
                                part["reason"],
                                part["needed_capability"],
                            ))
                            continue

                        seq = rec["seq"]
                        acc = rec["accession"]
                        hdr = f"{acc} {rec['definition']}".strip().encode("ascii", errors="replace").decode("ascii")
                        w.write(b"\n".join(fasta.format_record(hdr, seq, 70)) + b"\n")
                        cu = ids.contig_uid(seq)
                        o1_records.append((
                            sn, comp, acc, rec["organism"], cu, rec["length"],
                            bool(rec["circular"]), 11, rec["n_cds"], rec["n_spliced_cds"],
                            rec["n_introns"], rec["n_homing_rt"],
                        ))
                        by_comp_o1[comp] += 1
                        written += 1
        finally:
            w.close()
            if raw_fh is not None:
                raw_fh.close()

        if written == 0:
            fa_dst.unlink(missing_ok=True)
            continue
        objects.append({
            "path": str(obj_dir / out_name),
            "source_name": sn,
            "source_release": src_rel,
            "license": "INSDC / NCBI RefSeq (public domain)",
            "redistributable": True,
            "retrieved_at": retrieved,
        })

    if not o1_records:
        raise Refused("no table-11 organellar DNA records found for O1")

    all_uids = [r[4] for r in o1_records]
    owners = v2.owners_of(all_uids, root, root_stages, owner_inputs)

    listing = make_inputs.Listing()
    genomes: list[dict] = []
    rec_rows: list[tuple] = []
    meta_rows: list[tuple] = []
    repeated: list[tuple] = []
    shared: list[tuple] = []
    first_of: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}
    contig_comp: dict[bytes, str] = {}

    for ordinal, (sn, comp, acc, org, cu, length, circ, tt, n_cds, n_spl, n_int, n_hom) in enumerate(o1_records):
        own = owners.get(cu, "")
        rec_rows.append((
            ordinal, sn, comp, acc, org, cu.hex(), length, str(circ).lower(),
            tt, n_cds, n_spl, n_int, n_hom, own,
        ))
        if own:
            shared.append((ordinal, acc, cu.hex(), own))
        contig_len.setdefault(cu, length)
        contig_comp.setdefault(cu, comp)
        g = listing.add([cu], 11, genome_kind, "complete")
        guid = ids.genome_uid([cu])
        if g is None:
            repeated.append((ordinal, acc, guid.hex(), first_of[guid]))
            continue
        g["circular"] = [cu.hex()] if circ else []
        genomes.append(g)
        first_of[guid] = acc
        meta_rows.append((
            "O1", sn, comp, guid.hex(), cu.hex(), acc, org, length,
            str(circ).lower(), tt, n_cds, n_spl, n_int, n_hom, own,
        ))

    meta_rows.sort(key=lambda t: (t[2], t[1], t[5], t[4]))
    left_rows.sort(key=lambda t: (t[13], t[2], t[0]))

    make_inputs._write_tsv(work / "records.tsv", O1_RECORD_HEAD, rec_rows)
    make_inputs._write_tsv(work / "o1_metadata.tsv", O1_META_HEAD, meta_rows)
    make_inputs._write_tsv(work / "repeated_records.tsv", ["ordinal", "accession", "genome_uid", "first_accession"], repeated)
    make_inputs._write_tsv(work / "shared_records.tsv", ["ordinal", "accession", "contig_uid", "owner_stage"], shared)

    # Deterministic gzip for O1_left_for_c59a.tsv.gz
    with open(work / "O1_left_for_c59a.tsv.gz", "wb") as raw_l:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_l, compresslevel=1, mtime=0) as gz:
            gz.write(("\t".join(O1_LEFT_HEAD) + "\n").encode("utf-8"))
            for row in left_rows:
                gz.write(("\t".join(str(x) for x in row) + "\n").encode("utf-8"))

    (work / "objects.json").write_text(json.dumps(objects, indent=1))
    (work / "genomes.json").write_text(json.dumps(genomes))

    owner_contigs: Counter[str] = Counter()
    owner_bases: Counter[str] = Counter()
    distinct_by_comp: Counter[str] = Counter()
    bases_by_comp: Counter[str] = Counter()
    for cu, n in contig_len.items():
        c = contig_comp[cu]
        distinct_by_comp[c] += 1
        bases_by_comp[c] += n
        if cu in owners:
            st = owners[cu]
            owner_contigs[st] += 1
            owner_bases[st] += n

    held_contigs = sum(owner_contigs.values())
    held_bases = sum(owner_bases.values())

    return {
        "o1_records": len(o1_records),
        "o1_records_by_compartment": dict(sorted(by_comp_o1.items())),
        "distinct_contigs_by_compartment": dict(sorted(distinct_by_comp.items())),
        "distinct_bases_by_compartment": dict(sorted(bases_by_comp.items())),
        "left_for_c59a": len(left_rows),
        "left_for_c59a_by_reason": dict(sorted(by_reason_left.items())),
        "genome_kind": genome_kind,
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
        "root_stages": root_stages,
        "owner_inputs": list(owner_inputs),
        "build_order": [*root_stages, *owner_inputs, "O1"],
        "objects": len(objects),
        **listing.report(),
    }


def _default_sources(raw_dir: Path) -> list[dict]:
    specs = []
    for i in (1, 2, 3):
        p = raw_dir / f"plastid.{i}.genomic.gbff.gz"
        if p.is_file():
            specs.append({
                "path": str(p),
                "source_name": "ncbi_refseq_plastid_r238",
                "source_release": "238",
                "source_compartment": "plastid",
            })
    mp = raw_dir / "mitochondrion.1.genomic.gbff.gz"
    if mp.is_file():
        specs.append({
            "path": str(mp),
            "source_name": "ncbi_refseq_mitochondrion_r238",
            "source_release": "238",
            "source_compartment": "mitochondrion",
        })
    return specs


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Card C-59: NCBI RefSeq Organellar genomes (O1)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_aud = sub.add_parser("audit", help="Audit transl_table/introns and benchmark pyrodigal recall")
    p_aud.add_argument("--raw-dir", type=Path, required=True)
    p_aud.add_argument("--max-per-stratum", type=int, default=30)
    p_aud.add_argument("--seed", type=int, default=59)
    p_aud.add_argument("--out-json", type=Path, default=None)

    p_st = sub.add_parser("stage", help="Stage O1 inputs for transl_table=11 organellar genomes")
    p_st.add_argument("--raw-dir", type=Path, required=True)
    p_st.add_argument("--out", type=Path, required=True)
    p_st.add_argument("--root", type=Path, default=None)
    p_st.add_argument("--root-stages", nargs="*", default=[])
    p_st.add_argument("--owner-inputs", nargs="*", default=[])
    p_st.add_argument("--retrieved", required=True)
    p_st.add_argument("--genome-kind", choices=VALID_KINDS, default="plasmid")
    p_st.add_argument("--compress-fasta", action="store_true")
    p_st.add_argument("--out-json", type=Path, default=None)

    args = ap.parse_args(argv)
    specs = _default_sources(args.raw_dir)
    if not specs:
        raise Refused(f"no RefSeq organelle .gbff.gz files found in {args.raw_dir}")

    if args.cmd == "audit":
        rep = audit_gbff_files(specs, max_per_stratum=args.max_per_stratum, seed=args.seed)
    else:
        owner_inputs = _parse_owner_inputs(args.owner_inputs)
        rep = stage_o1(
            specs,
            args.out,
            root=args.root,
            root_stages=args.root_stages,
            owner_inputs=owner_inputs,
            retrieved=args.retrieved,
            genome_kind=args.genome_kind,
            compress_fasta=args.compress_fasta,
        )

    txt = json.dumps(rep, indent=2, sort_keys=True)
    if args.out_json is not None:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(txt + "\n")
    print(txt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
