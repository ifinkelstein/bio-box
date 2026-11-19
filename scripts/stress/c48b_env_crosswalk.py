#!/usr/bin/env python3
"""Card C-48b: JGI GOLD and IMG/M environmental metadata crosswalk and DuckDB sidecar tables.

Streaming, bounded-memory builder (1 CPU thread, deterministic Parquet output):
- Streams GOLD v.10 .xlsx sheets (inlineStr and sharedStrings) with PII/contact column stripping.
- Links IMG/M microbiomes (and optional isolate/virus slices) + GOLD Gs/Gb/Gp/Ga IDs to
  SRA runs, SRA experiments, and BioSample accessions.
- Cross-links against on-disk SPIRE study tables (and optional C-43 units / MGnify / INSDC tables)
  to produce the updated non-redundant run take-list stratified by GOLD 5-level ecosystem.
- Builds deterministic Parquet dimension/crosswalk tables and a DuckDB sidecar (`env_sidecar.duckdb`)
  with dimension views so queries join environmental attributes without copying them into fact tables
  (AGENTS.md §3.6).
- Measures environmental metadata coverage across stored stages (S1, P1..P6, H1, PL1) and staged
  inputs (V1, PL2, G1).
"""

from __future__ import annotations

import argparse
import ast
import csv
import gc
import gzip
import hashlib
import json
import os
import re
import tempfile
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
from typing import Iterator, Sequence

import blake3
import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

ROW_GROUP_SIZE = 122_880
BATCH_SIZE = 50_000

_PII_PATTERNS = (
    re.compile(r"\be[-_ ]?mail\b", re.I),
    re.compile(r"\bcontact", re.I),
    re.compile(r"\bsubmitter\b", re.I),
    re.compile(r"\bprincipal\s+investigator\b", re.I),
    re.compile(r"\bpi[-_ ]?name\b", re.I),
    re.compile(r"\bpi[-_ ]?email\b", re.I),
    re.compile(r"\b(?:tele)?phone\b", re.I),
    re.compile(r"\bfax\b", re.I),
    re.compile(r"\bstreet\b", re.I),
    re.compile(r"\bpostal\b", re.I),
)

_SRA_RUN_RE = re.compile(r"\b([EDS]RR\d+)\b")
_SRA_EXP_RE = re.compile(r"\b([EDS]RX\d+)\b")
_BIOSAMPLE_RE = re.compile(r"\b(SAM(?:EA|N|D)\d+|[EDS]RS\d+)\b")
_ASM_ACC_RE = re.compile(r"\b(GC[AF]_\d+)(?:\.\d+)?\b")
_IMGPR_OID_RE = re.compile(r"^IMGPR_plasmid_(\d+)_")
_CELL_COL_RE = re.compile(r"^([A-Z]+)")
_WS_RE = re.compile(r"\s+")

HABITAT_ORDER = ("host-associated", "anthropogenic", "terrestrial", "aquatic", "air")
UNLABELLED = "unlabelled"


def is_pii_column(name: str) -> bool:
    """True when a column header represents contact/PII metadata (AGENTS.md §8)."""
    return any(pat.search(name) for pat in _PII_PATTERNS)


def _norm_header(h: str) -> str:
    return _WS_RE.sub(" ", (h or "").strip())


def _col_ref_to_index(ref: str) -> int:
    m = _CELL_COL_RE.match(ref or "")
    if not m:
        return -1
    idx = 0
    for ch in m.group(1):
        idx = idx * 26 + (ord(ch) - 64)
    return idx - 1


def env_uid(source: str, source_id: str) -> bytes:
    """Deterministic 16-byte BLAKE3 UID for an environmental dimension entity."""
    return blake3.blake3(f"env:{source}:{source_id}".encode("utf-8")).digest(length=16)


def iter_xlsx_rows(xlsx_path: Path | str, sheet_name: str) -> Iterator[dict[str, str]]:
    """Stream rows from an OOXML .xlsx worksheet in bounded memory, stripping PII columns."""
    with zipfile.ZipFile(xlsx_path, "r") as zf:
        wb = ET.fromstring(zf.read("xl/workbook.xml"))
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rmap = {r.attrib["Id"]: r.attrib["Target"] for r in rels}

        target: str | None = None
        for sh in wb.iter():
            if sh.tag.endswith("}sheet") and sh.attrib.get("name") == sheet_name:
                rid = next(v for k, v in sh.attrib.items() if k.endswith("}id"))
                target = rmap[rid]
                if not target.startswith("xl/"):
                    target = "xl/" + target.lstrip("/")
                break
        if target is None:
            return

        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            with zf.open("xl/sharedStrings.xml") as ssf:
                for _, elem in ET.iterparse(ssf, events=("end",)):
                    if elem.tag.endswith("}si"):
                        shared_strings.append("".join(elem.itertext()))
                        elem.clear()

        headers: list[str] | None = None
        keep_indices: list[tuple[int, str]] = []
        with zf.open(target) as sf:
            for _, elem in ET.iterparse(sf, events=("end",)):
                if not elem.tag.endswith("}row"):
                    continue
                vals: dict[int, str] = {}
                seq_idx = 0
                for c in elem:
                    if not c.tag.endswith("}c"):
                        continue
                    c_idx = _col_ref_to_index(c.attrib.get("r", ""))
                    if c_idx < 0:
                        c_idx = seq_idx
                    seq_idx = c_idx + 1
                    t = c.attrib.get("t", "")
                    v = ""
                    for ch in c:
                        if ch.tag.endswith("}v"):
                            v = ch.text or ""
                        elif ch.tag.endswith("}is"):
                            v = "".join(ch.itertext())
                    if t == "s" and v.isdigit():
                        v = shared_strings[int(v)]
                    vals[c_idx] = v.strip()
                elem.clear()

                if headers is None:
                    max_c = max(vals.keys(), default=-1)
                    headers = [_norm_header(vals.get(i, "")) for i in range(max_c + 1)]
                    keep_indices = [
                        (i, h) for i, h in enumerate(headers) if h and not is_pii_column(h)
                    ]
                    continue

                if not vals:
                    continue
                yield {h: vals.get(i, "") for i, h in keep_indices}


def harmonize_gold_ecosystem(
    eco_l1: str = "",
    eco_l2: str = "",
    eco_l3: str = "",
    eco_l4: str = "",
    eco_l5: str = "",
) -> str:
    """Map GOLD 5-level ecosystem path to canonical top-level habitat."""
    l1 = (eco_l1 or "").strip().lower()
    l2 = (eco_l2 or "").strip().lower()
    if l1 == "host-associated":
        return "host-associated"
    if l1 == "engineered":
        return "anthropogenic"
    if l1 == "environmental":
        if l2 == "aquatic":
            return "aquatic"
        if l2 == "terrestrial":
            return "terrestrial"
        if l2 == "air":
            return "air"
    return UNLABELLED


def make_ecosystem_path(l1: str, l2: str, l3: str, l4: str, l5: str) -> str:
    parts = [p.strip() for p in (l1, l2, l3, l4, l5) if (p or "").strip()]
    return "; ".join(parts)


def parse_spire_microntology(cell: str) -> list[str]:
    cell = (cell or "").strip()
    if not cell:
        return []
    v = ast.literal_eval(cell)
    if not isinstance(v, list) or not all(isinstance(x, str) for x in v):
        raise ValueError(f"invalid microntology cell: {cell[:80]!r}")
    return v


def spire_habitat(terms: Sequence[str]) -> tuple[str, str]:
    tops = {t.split(":", 1)[0] for t in terms}
    for top in HABITAT_ORDER:
        if top in tops:
            seconds = sorted(":".join(t.split(":")[:2]) for t in terms if t.startswith(top + ":"))
            return top, (seconds[0] if seconds else top)
    return UNLABELLED, UNLABELLED


def extract_sra_runs(*texts: str) -> list[str]:
    found: set[str] = set()
    for t in texts:
        if t:
            found.update(_SRA_RUN_RE.findall(t))
    return sorted(found)


def extract_sra_exps(*texts: str) -> list[str]:
    found: set[str] = set()
    for t in texts:
        if t:
            found.update(_SRA_EXP_RE.findall(t))
    return sorted(found)


def extract_biosamples(*texts: str) -> list[str]:
    found: set[str] = set()
    for t in texts:
        if t:
            found.update(_BIOSAMPLE_RE.findall(t))
    return sorted(found)


def normalize_asm_acc(text: str) -> str:
    m = _ASM_ACC_RE.search(text or "")
    return m.group(1) if m else ""


def _parse_float(val: str) -> float | None:
    s = (val or "").strip()
    if not s:
        return None
    # Extract leading numeric token if units follow (e.g. "18.5 C")
    m = re.match(r"^([+-]?\d+(?:\.\d+)?)", s)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


def _parse_int(val: str) -> int | None:
    s = (val or "").strip().replace(",", "")
    if not s:
        return None
    try:
        return int(float(s))
    except ValueError:
        return None


def _norm_utilization_status(val: str) -> str:
    s = (val or "").strip()
    sl = s.lower()
    if sl == "unrestricted":
        return "Unrestricted"
    if sl == "restricted":
        return "Restricted"
    return s


def _open_tsv(path: Path | str):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt", encoding="utf-8", newline="")
    return open(p, "rt", encoding="utf-8", newline="")


def _write_arrow_chunks_to_parquet(
    rows_iter: Iterator[dict[str, object]],
    schema: pa.Schema,
    out_path: Path,
    *,
    batch_size: int = BATCH_SIZE,
) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(f".{out_path.name}.tmp")
    total = 0
    names = schema.names
    buf: dict[str, list[object]] = {k: [] for k in names}
    with pq.ParquetWriter(tmp, schema, compression="zstd", write_statistics=True) as writer:
        for row in rows_iter:
            for k in names:
                buf[k].append(row.get(k))
            if len(buf[names[0]]) >= batch_size:
                table = pa.Table.from_pydict(buf, schema=schema)
                writer.write_table(table, row_group_size=ROW_GROUP_SIZE)
                total += table.num_rows
                buf = {k: [] for k in names}
        if buf[names[0]]:
            table = pa.Table.from_pydict(buf, schema=schema)
            writer.write_table(table, row_group_size=ROW_GROUP_SIZE)
            total += table.num_rows
    os.replace(tmp, out_path)
    return total


def _export_duckdb_query_deterministic(
    con: duckdb.DuckDBPyConnection,
    sql: str,
    schema: pa.Schema,
    out_path: Path,
) -> int:
    """Execute sorted DuckDB query and stream batches to a deterministic Parquet file."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(f".{out_path.name}.tmp")
    cur = con.execute(sql)
    names = schema.names
    total = 0
    with pq.ParquetWriter(tmp, schema, compression="zstd", write_statistics=True) as writer:
        while True:
            rows = cur.fetchmany(ROW_GROUP_SIZE)
            if not rows:
                break
            cols = [list(col) for col in zip(*rows)]
            tbl = pa.Table.from_pydict(dict(zip(names, cols)), schema=schema)
            writer.write_table(tbl, row_group_size=ROW_GROUP_SIZE)
            total += tbl.num_rows
    os.replace(tmp, out_path)
    return total


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ------------------------------------------------------------------------------ Schemas

_GOLD_AP_SCHEMA = pa.schema([
    pa.field("ap_gold_id", pa.string()),
    pa.field("ap_type", pa.string()),
    pa.field("img_taxon_oid", pa.string()),
    pa.field("ap_genbank", pa.string()),
    pa.field("ap_sra_runs", pa.string()),
    pa.field("ap_sra_exps", pa.string()),
    pa.field("study_gold_id", pa.string()),
    pa.field("organism_gold_id", pa.string()),
    pa.field("project_gold_ids", pa.string()),
    pa.field("jgi_status", pa.string()),
    pa.field("estimated_size", pa.int64()),
])

_GOLD_SP_SCHEMA = pa.schema([
    pa.field("sp_gold_id", pa.string()),
    pa.field("study_gold_id", pa.string()),
    pa.field("biosample_gold_id", pa.string()),
    pa.field("organism_gold_id", pa.string()),
    pa.field("ncbi_bioproject", pa.string()),
    pa.field("ncbi_biosample", pa.string()),
    pa.field("sra_runs", pa.string()),
    pa.field("sra_exps", pa.string()),
    pa.field("jgi_status", pa.string()),
])

_GOLD_BS_SCHEMA = pa.schema([
    pa.field("biosample_gold_id", pa.string()),
    pa.field("eco_l1", pa.string()),
    pa.field("eco_l2", pa.string()),
    pa.field("eco_l3", pa.string()),
    pa.field("eco_l4", pa.string()),
    pa.field("eco_l5", pa.string()),
    pa.field("geo_location", pa.string()),
    pa.field("latitude", pa.float64()),
    pa.field("longitude", pa.float64()),
])

_GOLD_ORG_SCHEMA = pa.schema([
    pa.field("organism_gold_id", pa.string()),
    pa.field("eco_l1", pa.string()),
    pa.field("eco_l2", pa.string()),
    pa.field("eco_l3", pa.string()),
    pa.field("eco_l4", pa.string()),
    pa.field("eco_l5", pa.string()),
    pa.field("habitat", pa.string()),
    pa.field("ph", pa.float64()),
    pa.field("temperature_c", pa.float64()),
    pa.field("temperature_range", pa.string()),
    pa.field("salinity", pa.string()),
    pa.field("salinity_conc", pa.string()),
    pa.field("oxygen_req", pa.string()),
    pa.field("host_name", pa.string()),
    pa.field("body_site", pa.string()),
    pa.field("body_subsite", pa.string()),
    pa.field("geo_location", pa.string()),
    pa.field("latitude", pa.float64()),
    pa.field("longitude", pa.float64()),
])

_IMG_RAW_SCHEMA = pa.schema([
    pa.field("slice_name", pa.string()),
    pa.field("img_taxon_oid", pa.string()),
    pa.field("ncbi_domain", pa.string()),
    pa.field("ap_type", pa.string()),
    pa.field("ap_gold_id", pa.string()),
    pa.field("sp_gold_id", pa.string()),
    pa.field("study_gold_id", pa.string()),
    pa.field("sequencing_center", pa.string()),
    pa.field("ncbi_asm_acc", pa.string()),
    pa.field("ncbi_bioproject", pa.string()),
    pa.field("ncbi_biosample", pa.string()),
    pa.field("ncbi_genbank", pa.string()),
    pa.field("eco_l1", pa.string()),
    pa.field("eco_l2", pa.string()),
    pa.field("eco_l3", pa.string()),
    pa.field("eco_l4", pa.string()),
    pa.field("eco_l5", pa.string()),
    pa.field("habitat", pa.string()),
    pa.field("ph", pa.float64()),
    pa.field("temperature_c", pa.float64()),
    pa.field("temperature_range", pa.string()),
    pa.field("salinity", pa.string()),
    pa.field("salinity_conc", pa.string()),
    pa.field("oxygen_req", pa.string()),
    pa.field("oxygen_conc", pa.string()),
    pa.field("depth_m", pa.float64()),
    pa.field("elevation_m", pa.float64()),
    pa.field("subsurface_m", pa.float64()),
    pa.field("latitude", pa.float64()),
    pa.field("longitude", pa.float64()),
    pa.field("longhurst_code", pa.string()),
    pa.field("host_name", pa.string()),
    pa.field("body_site", pa.string()),
    pa.field("body_subsite", pa.string()),
    pa.field("geo_location", pa.string()),
    pa.field("asm_bases", pa.int64()),
    pa.field("sra_runs", pa.string()),
    pa.field("sra_exps", pa.string()),
    pa.field("jgi_status", pa.string()),
])

_SPIRE_RAW_SCHEMA = pa.schema([
    pa.field("sample_id", pa.string()),
    pa.field("study", pa.string()),
    pa.field("biosample_acc", pa.string()),
    pa.field("sra_run", pa.string()),
    pa.field("harmonized_habitat", pa.string()),
    pa.field("habitat_l2", pa.string()),
    pa.field("microntology", pa.string()),
    pa.field("latitude", pa.float64()),
    pa.field("longitude", pa.float64()),
])

_C43_REF_SCHEMA = pa.schema([
    pa.field("resource", pa.string()),
    pa.field("uid", pa.string()),
    pa.field("run_acc", pa.string()),
    pa.field("sample_acc", pa.string()),
    pa.field("habitat", pa.string()),
    pa.field("biome", pa.string()),
])

ENV_SAMPLE_DIM_SCHEMA = pa.schema([
    pa.field("env_uid", pa.binary(16), nullable=False),
    pa.field("source", pa.string(), nullable=False),
    pa.field("source_id", pa.string(), nullable=False),
    pa.field("img_taxon_oid", pa.string()),
    pa.field("gold_study_id", pa.string()),
    pa.field("gold_biosample_id", pa.string()),
    pa.field("gold_sp_id", pa.string()),
    pa.field("gold_ap_id", pa.string()),
    pa.field("biosample_acc", pa.string()),
    pa.field("sra_runs", pa.string()),
    pa.field("sra_experiments", pa.string()),
    pa.field("insdc_assembly_acc", pa.string()),
    pa.field("insdc_genbank_acc", pa.string()),
    pa.field("eco_l1", pa.string()),
    pa.field("eco_l2", pa.string()),
    pa.field("eco_l3", pa.string()),
    pa.field("eco_l4", pa.string()),
    pa.field("eco_l5", pa.string()),
    pa.field("ecosystem_path", pa.string()),
    pa.field("harmonized_habitat", pa.string(), nullable=False),
    pa.field("spire_microntology", pa.string()),
    pa.field("mgnify_biome", pa.string()),
    pa.field("habitat", pa.string()),
    pa.field("ph", pa.float64()),
    pa.field("temperature_c", pa.float64()),
    pa.field("temperature_range", pa.string()),
    pa.field("salinity", pa.string()),
    pa.field("salinity_conc", pa.string()),
    pa.field("oxygen_req", pa.string()),
    pa.field("oxygen_conc", pa.string()),
    pa.field("depth_m", pa.float64()),
    pa.field("elevation_m", pa.float64()),
    pa.field("subsurface_m", pa.float64()),
    pa.field("latitude", pa.float64()),
    pa.field("longitude", pa.float64()),
    pa.field("longhurst_code", pa.string()),
    pa.field("host_name", pa.string()),
    pa.field("body_site", pa.string()),
    pa.field("body_subsite", pa.string()),
    pa.field("geo_location", pa.string()),
    pa.field("jgi_utilization_status", pa.string()),
    pa.field("asm_bases", pa.int64()),
    pa.field("redistributable", pa.bool_(), nullable=False),
])

IMG_M_CROSSWALK_SCHEMA = pa.schema([
    pa.field("img_taxon_oid", pa.string(), nullable=False),
    pa.field("gold_ap_id", pa.string()),
    pa.field("gold_sp_id", pa.string()),
    pa.field("gold_study_id", pa.string()),
    pa.field("gold_biosample_id", pa.string()),
    pa.field("biosample_acc", pa.string()),
    pa.field("sra_runs", pa.string()),
    pa.field("sra_experiments", pa.string()),
    pa.field("eco_l1", pa.string()),
    pa.field("eco_l2", pa.string()),
    pa.field("eco_l3", pa.string()),
    pa.field("eco_l4", pa.string()),
    pa.field("eco_l5", pa.string()),
    pa.field("ecosystem_path", pa.string()),
    pa.field("harmonized_habitat", pa.string(), nullable=False),
    pa.field("asm_bases", pa.int64()),
    pa.field("jgi_utilization_status", pa.string()),
    pa.field("is_unrestricted", pa.bool_(), nullable=False),
    pa.field("is_non_restricted", pa.bool_(), nullable=False),
    pa.field("in_spire", pa.bool_(), nullable=False),
    pa.field("in_mgnify", pa.bool_(), nullable=False),
    pa.field("in_insdc", pa.bool_(), nullable=False),
    pa.field("decision", pa.string(), nullable=False),
    pa.field("redistributable", pa.bool_(), nullable=False),
])

GENOME_ENV_LINK_SCHEMA = pa.schema([
    pa.field("stage", pa.string(), nullable=False),
    pa.field("genome_key", pa.string(), nullable=False),
    pa.field("source_accession", pa.string()),
    pa.field("link_method", pa.string(), nullable=False),
    pa.field("env_uid", pa.binary(16)),
    pa.field("env_source", pa.string()),
    pa.field("harmonized_habitat", pa.string(), nullable=False),
    pa.field("ecosystem_path", pa.string()),
    pa.field("redistributable", pa.bool_(), nullable=False),
])


# ------------------------------------------------------------------------------ Staging readers

def _stage_gold_data(gold_data: Path | None, tmp_dir: Path) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    ap_path = tmp_dir / "gold_ap.parquet"
    sp_path = tmp_dir / "gold_sp.parquet"
    bs_path = tmp_dir / "gold_bs.parquet"
    org_path = tmp_dir / "gold_org.parquet"

    def _ap_rows() -> Iterator[dict[str, object]]:
        if not gold_data or not Path(gold_data).is_file():
            return
        for r in iter_xlsx_rows(gold_data, "Analysis Project"):
            jgi_st = ""
            for k, v in r.items():
                if "JGI DATA UTILIZATION STATUS" in k.upper():
                    jgi_st = _norm_utilization_status(v)
                    break
            sra_raw = r.get("AP NCBI SRA IDS", "")
            yield {
                "ap_gold_id": r.get("AP GOLD ID", ""),
                "ap_type": r.get("AP TYPE", ""),
                "img_taxon_oid": r.get("AP IMG TAXON ID", ""),
                "ap_genbank": r.get("AP GENBANK", ""),
                "ap_sra_runs": ";".join(extract_sra_runs(sra_raw)),
                "ap_sra_exps": ";".join(extract_sra_exps(sra_raw)),
                "study_gold_id": r.get("AP STUDY GOLD ID", ""),
                "organism_gold_id": r.get("AP ORGANISM GOLD ID", ""),
                "project_gold_ids": r.get("AP PROJECT GOLD IDS", ""),
                "jgi_status": jgi_st,
                "estimated_size": _parse_int(r.get("AP ESTIMATED SIZE", "")),
            }

    def _sp_rows() -> Iterator[dict[str, object]]:
        if not gold_data or not Path(gold_data).is_file():
            return
        for r in iter_xlsx_rows(gold_data, "Sequencing Project"):
            jgi_st = ""
            for k, v in r.items():
                if "JGI DATA UTILIZATION STATUS" in k.upper():
                    jgi_st = _norm_utilization_status(v)
                    break
            sra_raw = r.get("SRA EXPERIMENT IDS", "")
            yield {
                "sp_gold_id": r.get("PROJECT GOLD ID", ""),
                "study_gold_id": r.get("STUDY GOLD ID", ""),
                "biosample_gold_id": r.get("BIOSAMPLE GOLD ID", ""),
                "organism_gold_id": r.get("ORGANISM GOLD ID", ""),
                "ncbi_bioproject": r.get("NCBI BIOPROJECT ACCESSION", ""),
                "ncbi_biosample": ";".join(extract_biosamples(r.get("NCBI BIOSAMPLE ACCESSION", ""))),
                "sra_runs": ";".join(extract_sra_runs(sra_raw)),
                "sra_exps": ";".join(extract_sra_exps(sra_raw)),
                "jgi_status": jgi_st,
            }

    def _bs_rows() -> Iterator[dict[str, object]]:
        if not gold_data or not Path(gold_data).is_file():
            return
        for r in iter_xlsx_rows(gold_data, "Biosample"):
            yield {
                "biosample_gold_id": r.get("BIOSAMPLE GOLD ID", ""),
                "eco_l1": r.get("BIOSAMPLE ECOSYSTEM", ""),
                "eco_l2": r.get("BIOSAMPLE ECOSYSTEM CATEGORY", ""),
                "eco_l3": r.get("BIOSAMPLE ECOSYSTEM TYPE", ""),
                "eco_l4": r.get("BIOSAMPLE ECOSYSTEM SUBTYPE", ""),
                "eco_l5": r.get("BIOSAMPLE SPECIFIC ECOSYSTEM", ""),
                "geo_location": r.get("BIOSAMPLE GEOGRAPHIC LOCATION", ""),
                "latitude": _parse_float(r.get("BIOSAMPLE LATITUDE", "")),
                "longitude": _parse_float(r.get("BIOSAMPLE LONGITUDE", "")),
            }

    def _org_rows() -> Iterator[dict[str, object]]:
        if not gold_data or not Path(gold_data).is_file():
            return
        for r in iter_xlsx_rows(gold_data, "Organism"):
            yield {
                "organism_gold_id": r.get("ORGANISM GOLD ID", ""),
                "eco_l1": r.get("ORGANISM ECOSYSTEM", ""),
                "eco_l2": r.get("ORGANISM ECOSYSTEM CATEGORY", ""),
                "eco_l3": r.get("ORGANISM ECOSYSTEM TYPE", ""),
                "eco_l4": r.get("ORGANISM ECOSYSTEM SUBTYPE", ""),
                "eco_l5": r.get("ORGANISM SPECIFIC ECOSYSTEM", ""),
                "habitat": r.get("ORGANISM HABITAT", ""),
                "ph": _parse_float(r.get("ORGANISM PH", "")),
                "temperature_c": _parse_float(r.get("ORGANISM SAMPLE COLLECTION TEMPERATURE", "")),
                "temperature_range": r.get("ORGANISM TEMPERATURE RANGE", ""),
                "salinity": r.get("ORGANISM SALINITY", ""),
                "salinity_conc": r.get("ORGANISM SALINITY CONCENTRATION", ""),
                "oxygen_req": r.get("ORGANISM OXYGEN REQUIREMENT", ""),
                "host_name": r.get("ORGANISM ISOLATION HOST NAME", ""),
                "body_site": r.get("ORGANISM HOST BODY SITE", ""),
                "body_subsite": r.get("ORGANISM HOST BODY SUBSITE", ""),
                "geo_location": r.get("ORGANISM GEOGRAPHIC LOCATION", ""),
                "latitude": _parse_float(r.get("ORGANISM LATITUDE", "")),
                "longitude": _parse_float(r.get("ORGANISM LONGITUDE", "")),
            }

    _write_arrow_chunks_to_parquet(_ap_rows(), _GOLD_AP_SCHEMA, ap_path)
    gc.collect()
    _write_arrow_chunks_to_parquet(_sp_rows(), _GOLD_SP_SCHEMA, sp_path)
    gc.collect()
    _write_arrow_chunks_to_parquet(_bs_rows(), _GOLD_BS_SCHEMA, bs_path)
    gc.collect()
    _write_arrow_chunks_to_parquet(_org_rows(), _GOLD_ORG_SCHEMA, org_path)
    gc.collect()
    paths["ap"] = ap_path
    paths["sp"] = sp_path
    paths["bs"] = bs_path
    paths["org"] = org_path
    return paths


def _stage_img_m(img_files: Sequence[tuple[str, Path]], tmp_dir: Path) -> Path:
    out = tmp_dir / "img_raw.parquet"

    def _rows() -> Iterator[dict[str, object]]:
        for slice_name, path in img_files:
            if not path or not Path(path).is_file():
                continue
            with _open_tsv(path) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for r in reader:
                    oid = (r.get("IMG Genome ID") or "").strip()
                    if not oid:
                        continue
                    sra_id = r.get("SRA ID", "")
                    sra_run = r.get("SRA Run", "")
                    yield {
                        "slice_name": slice_name,
                        "img_taxon_oid": oid,
                        "ncbi_domain": (r.get("NCBI Domain") or "").strip(),
                        "ap_type": (r.get("GOLD Analysis Project Type") or "").strip(),
                        "ap_gold_id": (r.get("GOLD Analysis Project ID") or "").strip(),
                        "sp_gold_id": (r.get("GOLD Sequencing Project ID") or "").strip(),
                        "study_gold_id": (r.get("GOLD Study ID") or "").strip(),
                        "sequencing_center": (r.get("Sequencing Center") or "").strip(),
                        "ncbi_asm_acc": normalize_asm_acc(r.get("NCBI Assembly Accession") or ""),
                        "ncbi_bioproject": (r.get("NCBI Bioproject Accession") or "").strip(),
                        "ncbi_biosample": ";".join(extract_biosamples(r.get("NCBI Biosample Accession") or "")),
                        "ncbi_genbank": (r.get("NCBI GenBank ID") or "").strip(),
                        "eco_l1": (r.get("Ecosystem") or "").strip(),
                        "eco_l2": (r.get("Ecosystem Category") or "").strip(),
                        "eco_l3": (r.get("Ecosystem Type") or "").strip(),
                        "eco_l4": (r.get("Ecosystem Subtype") or "").strip(),
                        "eco_l5": (r.get("Specific Ecosystem") or "").strip(),
                        "habitat": (r.get("Habitat") or "").strip(),
                        "ph": _parse_float(r.get("pH") or ""),
                        "temperature_c": _parse_float(r.get("Sample Collection Temperature") or ""),
                        "temperature_range": (r.get("Temperature Range") or "").strip(),
                        "salinity": (r.get("Salinity") or "").strip(),
                        "salinity_conc": (r.get("Salinity Concentration") or "").strip(),
                        "oxygen_req": (r.get("Oxygen Requirement") or "").strip(),
                        "oxygen_conc": (r.get("Oxygen Concentration") or "").strip(),
                        "depth_m": _parse_float(r.get("Depth In Meters") or ""),
                        "elevation_m": _parse_float(r.get("Elevation In Meters") or ""),
                        "subsurface_m": _parse_float(r.get("Subsurface In Meters") or ""),
                        "latitude": _parse_float(r.get("Latitude") or ""),
                        "longitude": _parse_float(r.get("Longitude") or ""),
                        "longhurst_code": (r.get("Longhurst Code") or "").strip(),
                        "host_name": (r.get("Host Name") or "").strip(),
                        "body_site": (r.get("Sample Body Site") or "").strip(),
                        "body_subsite": (r.get("Sample Body Subsite") or "").strip(),
                        "geo_location": (r.get("Geographic Location") or "").strip(),
                        "asm_bases": _parse_int(r.get("Genome Size * assembled") or ""),
                        "sra_runs": ";".join(extract_sra_runs(sra_run, sra_id)),
                        "sra_exps": ";".join(extract_sra_exps(sra_id, sra_run)),
                        "jgi_status": _norm_utilization_status(r.get("JGI Data Utilization Status") or ""),
                    }

    _write_arrow_chunks_to_parquet(_rows(), _IMG_RAW_SCHEMA, out)
    gc.collect()
    return out


def _stage_spire(spire_dir: Path | None, tmp_dir: Path) -> Path:
    out = tmp_dir / "spire_raw.parquet"

    def _rows() -> Iterator[dict[str, object]]:
        if not spire_dir or not Path(spire_dir).is_dir():
            return
        d = Path(spire_dir)
        for name in sorted(os.listdir(d)):
            if not name.endswith(".tsv") or name.startswith("probe_"):
                continue
            study = name[:-4]
            with open(d / name, "rt", encoding="utf-8", newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for r in reader:
                    sid = (r.get("sample_id") or "").strip()
                    if not sid:
                        continue
                    terms = parse_spire_microntology(r.get("microntology") or "")
                    h1, h2 = spire_habitat(terms)
                    runs = extract_sra_runs(sid)
                    bs = extract_biosamples(sid)
                    yield {
                        "sample_id": sid,
                        "study": study,
                        "biosample_acc": ";".join(bs),
                        "sra_run": ";".join(runs),
                        "harmonized_habitat": h1,
                        "habitat_l2": h2,
                        "microntology": ";".join(sorted(set(terms))),
                        "latitude": _parse_float(r.get("lat") or ""),
                        "longitude": _parse_float(r.get("lon") or ""),
                    }

    _write_arrow_chunks_to_parquet(_rows(), _SPIRE_RAW_SCHEMA, out)
    gc.collect()
    return out


def _stage_c43_units(c43_units: Path | None, tmp_dir: Path) -> Path:
    out = tmp_dir / "c43_ref.parquet"

    def _rows() -> Iterator[dict[str, object]]:
        if not c43_units or not Path(c43_units).is_file():
            return
        with _open_tsv(c43_units) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for r in reader:
                res = (r.get("resource") or "").strip()
                uid = (r.get("uid") or "").strip()
                runs = [x for x in (r.get("runs") or "").split(";") if x]
                samples = [x for x in (r.get("samples") or "").split(";") if x]
                hab = (r.get("habitat") or UNLABELLED).strip()
                detail = (r.get("habitat_detail") or "").strip()
                if not runs and not samples:
                    yield {"resource": res, "uid": uid, "run_acc": "", "sample_acc": "", "habitat": hab, "biome": detail}
                else:
                    for run in runs:
                        yield {"resource": res, "uid": uid, "run_acc": run, "sample_acc": "", "habitat": hab, "biome": detail}
                    for samp in samples:
                        yield {"resource": res, "uid": uid, "run_acc": "", "sample_acc": samp, "habitat": hab, "biome": detail}

    _write_arrow_chunks_to_parquet(_rows(), _C43_REF_SCHEMA, out)
    gc.collect()
    return out


def _register_udfs(con: duckdb.DuckDBPyConnection) -> None:
    con.create_function(
        "harmonize_gold_eco",
        lambda l1, l2, l3, l4, l5: harmonize_gold_ecosystem(l1 or "", l2 or "", l3 or "", l4 or "", l5 or ""),
        ["VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR"],
        "VARCHAR",
        null_handling="special",
    )
    con.create_function(
        "make_eco_path",
        lambda l1, l2, l3, l4, l5: make_ecosystem_path(l1 or "", l2 or "", l3 or "", l4 or "", l5 or ""),
        ["VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR", "VARCHAR"],
        "VARCHAR",
        null_handling="special",
    )
    con.create_function(
        "make_env_uid",
        lambda src, sid: env_uid(src or "", sid or ""),
        ["VARCHAR", "VARCHAR"],
        "BLOB",
        null_handling="special",
    )
    con.create_function(
        "merge_semi",
        lambda a, b: ";".join(sorted({x for s in (a or "", b or "") for x in s.split(";") if x})),
        ["VARCHAR", "VARCHAR"],
        "VARCHAR",
        null_handling="special",
    )


def build_crosswalk(
    *,
    gold_data: Path | str | None = None,
    gold_sra_biome: Path | str | None = None,
    gold_eco_paths: Path | str | None = None,
    img_microbiomes: Path | str | None = None,
    img_extra_slices: Sequence[tuple[str, Path | str]] = (),
    spire_dir: Path | str | None = None,
    c43_units: Path | str | None = None,
    out_dir: Path | str,
) -> dict[str, object]:
    """Build bounded-memory, deterministic GOLD/IMG-M crosswalk + DuckDB environmental sidecar."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="c48b_", dir=str(out)) as td:
        tmp = Path(td)
        gold_paths = _stage_gold_data(Path(gold_data) if gold_data else None, tmp)

        img_list: list[tuple[str, Path]] = []
        if img_microbiomes:
            img_list.append(("microbiomes", Path(img_microbiomes)))
        for sname, spath in img_extra_slices:
            img_list.append((sname, Path(spath)))
        img_raw_path = _stage_img_m(img_list, tmp)
        spire_raw_path = _stage_spire(Path(spire_dir) if spire_dir else None, tmp)
        c43_ref_path = _stage_c43_units(Path(c43_units) if c43_units else None, tmp)

        # Optional: validate goldSraBiomeImg.xlsx & GOLDs5levelEcosystemClassificationPaths.xlsx counts
        sra_biome_rows = 0
        if gold_sra_biome and Path(gold_sra_biome).is_file():
            for _ in iter_xlsx_rows(gold_sra_biome, "SRA biomes with IMG annotations"):
                sra_biome_rows += 1
        eco_path_rows = 0
        if gold_eco_paths and Path(gold_eco_paths).is_file():
            for _ in iter_xlsx_rows(gold_eco_paths, "site data"):
                eco_path_rows += 1

        work_db = tmp / "work.duckdb"
        con = duckdb.connect(str(work_db))
        con.execute("SET threads=1")
        con.execute("SET memory_limit='2GB'")
        con.execute("SET preserve_insertion_order=false")
        _register_udfs(con)

        con.execute(f"CREATE VIEW gold_ap AS SELECT * FROM read_parquet('{gold_paths['ap']}')")
        con.execute(f"CREATE VIEW gold_sp AS SELECT * FROM read_parquet('{gold_paths['sp']}')")
        con.execute(f"CREATE VIEW gold_bs AS SELECT * FROM read_parquet('{gold_paths['bs']}')")
        con.execute(f"CREATE VIEW gold_org AS SELECT * FROM read_parquet('{gold_paths['org']}')")
        con.execute(f"CREATE VIEW img_raw AS SELECT * FROM read_parquet('{img_raw_path}')")
        con.execute(f"CREATE VIEW spire_raw AS SELECT * FROM read_parquet('{spire_raw_path}')")
        con.execute(f"CREATE VIEW c43_ref AS SELECT * FROM read_parquet('{c43_ref_path}')")

        # Deduplicate GOLD AP by ap_gold_id deterministically and link first SP -> BS
        con.execute("""
            CREATE TABLE gold_ap_dedup AS
            SELECT
                ap_gold_id,
                min(ap_type) AS ap_type,
                min(nullif(img_taxon_oid, '')) AS img_taxon_oid,
                min(nullif(ap_genbank, '')) AS ap_genbank,
                string_agg(DISTINCT nullif(ap_sra_runs, ''), ';' ORDER BY nullif(ap_sra_runs, '')) AS ap_sra_runs,
                string_agg(DISTINCT nullif(ap_sra_exps, ''), ';' ORDER BY nullif(ap_sra_exps, '')) AS ap_sra_exps,
                min(nullif(study_gold_id, '')) AS study_gold_id,
                min(nullif(organism_gold_id, '')) AS organism_gold_id,
                min(nullif(project_gold_ids, '')) AS project_gold_ids,
                min(nullif(jgi_status, '')) AS jgi_status,
                max(estimated_size) AS estimated_size
            FROM gold_ap
            WHERE ap_gold_id <> ''
            GROUP BY ap_gold_id
        """)

        con.execute("""
            CREATE TABLE gold_sp_dedup AS
            SELECT
                sp_gold_id,
                min(nullif(study_gold_id, '')) AS study_gold_id,
                min(nullif(biosample_gold_id, '')) AS biosample_gold_id,
                min(nullif(organism_gold_id, '')) AS organism_gold_id,
                min(nullif(ncbi_bioproject, '')) AS ncbi_bioproject,
                string_agg(DISTINCT nullif(ncbi_biosample, ''), ';' ORDER BY nullif(ncbi_biosample, '')) AS ncbi_biosample,
                string_agg(DISTINCT nullif(sra_runs, ''), ';' ORDER BY nullif(sra_runs, '')) AS sra_runs,
                string_agg(DISTINCT nullif(sra_exps, ''), ';' ORDER BY nullif(sra_exps, '')) AS sra_exps,
                min(nullif(jgi_status, '')) AS jgi_status
            FROM gold_sp
            WHERE sp_gold_id <> ''
            GROUP BY sp_gold_id
        """)

        # Build covered accession lookup sets from SPIRE and optional C-43 units
        con.execute("""
            CREATE TABLE covered_keys AS
            SELECT 'spire' AS resource, sample_id AS acc FROM spire_raw WHERE sample_id <> ''
            UNION
            SELECT 'spire' AS resource, biosample_acc AS acc FROM spire_raw WHERE biosample_acc <> ''
            UNION
            SELECT 'spire' AS resource, sra_run AS acc FROM spire_raw WHERE sra_run <> ''
            UNION
            SELECT resource, run_acc AS acc FROM c43_ref WHERE run_acc <> ''
            UNION
            SELECT resource, sample_acc AS acc FROM c43_ref WHERE sample_acc <> ''
        """)

        # Enrich IMG/M rows with GOLD AP + SP metadata
        con.execute("""
            CREATE TABLE img_enriched AS
            SELECT
                i.slice_name,
                i.img_taxon_oid,
                coalesce(nullif(i.ap_gold_id, ''), ap.ap_gold_id, '') AS gold_ap_id,
                coalesce(nullif(i.sp_gold_id, ''), split_part(coalesce(ap.project_gold_ids, ''), ',', 1), '') AS gold_sp_id,
                coalesce(nullif(i.study_gold_id, ''), ap.study_gold_id, sp.study_gold_id, '') AS gold_study_id,
                coalesce(sp.biosample_gold_id, '') AS gold_biosample_id,
                merge_semi(i.ncbi_biosample, coalesce(sp.ncbi_biosample, '')) AS biosample_acc,
                merge_semi(i.sra_runs, merge_semi(coalesce(ap.ap_sra_runs, ''), coalesce(sp.sra_runs, ''))) AS sra_runs,
                merge_semi(i.sra_exps, merge_semi(coalesce(ap.ap_sra_exps, ''), coalesce(sp.sra_exps, ''))) AS sra_experiments,
                coalesce(nullif(i.ncbi_asm_acc, ''), '') AS insdc_assembly_acc,
                coalesce(nullif(i.ncbi_genbank, ''), ap.ap_genbank, '') AS insdc_genbank_acc,
                coalesce(nullif(i.eco_l1, ''), bs.eco_l1, '') AS eco_l1,
                coalesce(nullif(i.eco_l2, ''), bs.eco_l2, '') AS eco_l2,
                coalesce(nullif(i.eco_l3, ''), bs.eco_l3, '') AS eco_l3,
                coalesce(nullif(i.eco_l4, ''), bs.eco_l4, '') AS eco_l4,
                coalesce(nullif(i.eco_l5, ''), bs.eco_l5, '') AS eco_l5,
                i.habitat,
                i.ph,
                i.temperature_c,
                i.temperature_range,
                i.salinity,
                i.salinity_conc,
                i.oxygen_req,
                i.oxygen_conc,
                i.depth_m,
                i.elevation_m,
                i.subsurface_m,
                coalesce(i.latitude, bs.latitude) AS latitude,
                coalesce(i.longitude, bs.longitude) AS longitude,
                i.longhurst_code,
                i.host_name,
                i.body_site,
                i.body_subsite,
                coalesce(nullif(i.geo_location, ''), bs.geo_location, '') AS geo_location,
                coalesce(nullif(i.jgi_status, ''), ap.jgi_status, sp.jgi_status, '') AS jgi_utilization_status,
                coalesce(i.asm_bases, ap.estimated_size) AS asm_bases
            FROM (
                SELECT *, row_number() OVER (PARTITION BY img_taxon_oid ORDER BY slice_name) AS rn
                FROM img_raw
            ) i
            LEFT JOIN gold_ap_dedup ap ON (i.ap_gold_id <> '' AND i.ap_gold_id = ap.ap_gold_id)
            LEFT JOIN gold_sp_dedup sp ON (
                coalesce(nullif(i.sp_gold_id, ''), split_part(coalesce(ap.project_gold_ids, ''), ',', 1)) = sp.sp_gold_id
                AND sp.sp_gold_id <> ''
            )
            LEFT JOIN (
                SELECT biosample_gold_id, min(eco_l1) AS eco_l1, min(eco_l2) AS eco_l2,
                       min(eco_l3) AS eco_l3, min(eco_l4) AS eco_l4, min(eco_l5) AS eco_l5,
                       min(geo_location) AS geo_location, min(latitude) AS latitude, min(longitude) AS longitude
                FROM gold_bs WHERE biosample_gold_id <> '' GROUP BY biosample_gold_id
            ) bs ON (sp.biosample_gold_id = bs.biosample_gold_id)
            WHERE i.rn = 1
        """)

        # Compute overlap of IMG/M microbiomes against SPIRE, MGnify, INSDC via runs and biosamples
        con.execute("""
            CREATE TABLE img_m_hits AS
            WITH tokens AS (
                SELECT img_taxon_oid, unnest(string_split(sra_runs, ';')) AS acc
                FROM img_enriched WHERE slice_name = 'microbiomes' AND sra_runs <> ''
                UNION
                SELECT img_taxon_oid, unnest(string_split(biosample_acc, ';')) AS acc
                FROM img_enriched WHERE slice_name = 'microbiomes' AND biosample_acc <> ''
            )
            SELECT
                t.img_taxon_oid,
                bool_or(c.resource = 'spire') AS in_spire,
                bool_or(c.resource = 'mgnify') AS in_mgnify,
                bool_or(c.resource = 'insdc') AS in_insdc
            FROM tokens t
            JOIN covered_keys c ON t.acc = c.acc
            GROUP BY t.img_taxon_oid
        """)

        # Write img_m_crosswalk.parquet (sorted by numeric/lex img_taxon_oid)
        crosswalk_sql = """
            SELECT
                e.img_taxon_oid,
                e.gold_ap_id,
                e.gold_sp_id,
                e.gold_study_id,
                e.gold_biosample_id,
                e.biosample_acc,
                e.sra_runs,
                e.sra_experiments,
                e.eco_l1,
                e.eco_l2,
                e.eco_l3,
                e.eco_l4,
                e.eco_l5,
                make_eco_path(e.eco_l1, e.eco_l2, e.eco_l3, e.eco_l4, e.eco_l5) AS ecosystem_path,
                harmonize_gold_eco(e.eco_l1, e.eco_l2, e.eco_l3, e.eco_l4, e.eco_l5) AS harmonized_habitat,
                e.asm_bases,
                e.jgi_utilization_status,
                (e.jgi_utilization_status = 'Unrestricted') AS is_unrestricted,
                (e.jgi_utilization_status <> 'Restricted') AS is_non_restricted,
                coalesce(h.in_spire, false) AS in_spire,
                coalesce(h.in_mgnify, false) AS in_mgnify,
                coalesce(h.in_insdc, false) AS in_insdc,
                CASE
                    WHEN e.jgi_utilization_status = 'Restricted' THEN 'restricted'
                    WHEN coalesce(h.in_spire, false) OR coalesce(h.in_mgnify, false) OR coalesce(h.in_insdc, false) THEN 'covered'
                    WHEN e.sra_runs <> '' OR e.biosample_acc <> '' OR e.sra_experiments <> '' THEN 'take'
                    ELSE 'take_unlinked'
                END AS decision,
                false AS redistributable
            FROM img_enriched e
            LEFT JOIN img_m_hits h ON e.img_taxon_oid = h.img_taxon_oid
            WHERE e.slice_name = 'microbiomes'
            ORDER BY e.img_taxon_oid
        """
        cw_path = out / "img_m_crosswalk.parquet"
        _export_duckdb_query_deterministic(con, crosswalk_sql, IMG_M_CROSSWALK_SCHEMA, cw_path)

        # Write env_sample_dim.parquet combining IMG/M (redistributable=false) and SPIRE (redistributable=true)
        env_dim_sql = """
            SELECT
                make_env_uid('img_m', e.img_taxon_oid) AS env_uid,
                'img_m' AS source,
                e.img_taxon_oid AS source_id,
                e.img_taxon_oid,
                e.gold_study_id,
                e.gold_biosample_id,
                e.gold_sp_id,
                e.gold_ap_id,
                e.biosample_acc,
                e.sra_runs,
                e.sra_experiments,
                e.insdc_assembly_acc,
                e.insdc_genbank_acc,
                e.eco_l1,
                e.eco_l2,
                e.eco_l3,
                e.eco_l4,
                e.eco_l5,
                make_eco_path(e.eco_l1, e.eco_l2, e.eco_l3, e.eco_l4, e.eco_l5) AS ecosystem_path,
                harmonize_gold_eco(e.eco_l1, e.eco_l2, e.eco_l3, e.eco_l4, e.eco_l5) AS harmonized_habitat,
                '' AS spire_microntology,
                '' AS mgnify_biome,
                e.habitat,
                e.ph,
                e.temperature_c,
                e.temperature_range,
                e.salinity,
                e.salinity_conc,
                e.oxygen_req,
                e.oxygen_conc,
                e.depth_m,
                e.elevation_m,
                e.subsurface_m,
                e.latitude,
                e.longitude,
                e.longhurst_code,
                e.host_name,
                e.body_site,
                e.body_subsite,
                e.geo_location,
                e.jgi_utilization_status,
                e.asm_bases,
                false AS redistributable
            FROM img_enriched e
            UNION ALL
            SELECT
                make_env_uid('spire', s.sample_id) AS env_uid,
                'spire' AS source,
                s.sample_id AS source_id,
                '' AS img_taxon_oid,
                '' AS gold_study_id,
                '' AS gold_biosample_id,
                '' AS gold_sp_id,
                '' AS gold_ap_id,
                s.biosample_acc,
                s.sra_run AS sra_runs,
                '' AS sra_experiments,
                '' AS insdc_assembly_acc,
                '' AS insdc_genbank_acc,
                '' AS eco_l1,
                '' AS eco_l2,
                '' AS eco_l3,
                '' AS eco_l4,
                '' AS eco_l5,
                s.habitat_l2 AS ecosystem_path,
                s.harmonized_habitat,
                s.microntology AS spire_microntology,
                '' AS mgnify_biome,
                s.habitat_l2 AS habitat,
                CAST(NULL AS DOUBLE) AS ph,
                CAST(NULL AS DOUBLE) AS temperature_c,
                '' AS temperature_range,
                '' AS salinity,
                '' AS salinity_conc,
                '' AS oxygen_req,
                '' AS oxygen_conc,
                CAST(NULL AS DOUBLE) AS depth_m,
                CAST(NULL AS DOUBLE) AS elevation_m,
                CAST(NULL AS DOUBLE) AS subsurface_m,
                s.latitude,
                s.longitude,
                '' AS longhurst_code,
                '' AS host_name,
                '' AS body_site,
                '' AS body_subsite,
                '' AS geo_location,
                '' AS jgi_utilization_status,
                CAST(NULL AS BIGINT) AS asm_bases,
                true AS redistributable
            FROM (
                SELECT *, row_number() OVER (PARTITION BY sample_id ORDER BY study) AS rn
                FROM spire_raw
            ) s
            WHERE s.rn = 1
            ORDER BY source, source_id
        """
        env_dim_path = out / "env_sample_dim.parquet"
        _export_duckdb_query_deterministic(con, env_dim_sql, ENV_SAMPLE_DIM_SCHEMA, env_dim_path)

        # Compute summary stats from the written Parquet files
        summary_row = con.execute(f"""
            SELECT
                count(*) AS total_units,
                count(*) FILTER (WHERE is_unrestricted) AS unrestricted_units,
                count(*) FILTER (WHERE is_non_restricted) AS non_restricted_units,
                count(*) FILTER (WHERE decision = 'restricted') AS restricted_units,
                count(*) FILTER (WHERE is_unrestricted AND decision IN ('take', 'take_unlinked')) AS novel_unrestricted_units,
                coalesce(sum(asm_bases) FILTER (WHERE is_unrestricted AND decision IN ('take', 'take_unlinked')), 0) AS novel_unrestricted_bases,
                count(*) FILTER (WHERE is_non_restricted AND decision IN ('take', 'take_unlinked')) AS novel_non_restricted_units,
                coalesce(sum(asm_bases) FILTER (WHERE is_non_restricted AND decision IN ('take', 'take_unlinked')), 0) AS novel_non_restricted_bases,
                count(*) FILTER (WHERE in_spire) AS in_spire_units,
                count(*) FILTER (WHERE in_mgnify) AS in_mgnify_units,
                count(*) FILTER (WHERE in_insdc) AS in_insdc_units,
                count(*) FILTER (WHERE sra_runs <> '') AS with_sra_run,
                count(*) FILTER (WHERE sra_experiments <> '') AS with_sra_exp,
                count(*) FILTER (WHERE biosample_acc <> '') AS with_biosample
            FROM read_parquet('{cw_path}')
        """).fetchone()

        strat_rows = con.execute(f"""
            SELECT
                harmonized_habitat,
                eco_l1,
                eco_l2,
                decision,
                count(*) AS units,
                count(*) FILTER (WHERE is_unrestricted) AS unrestricted_units,
                coalesce(sum(asm_bases), 0) AS asm_bases,
                coalesce(sum(asm_bases) FILTER (WHERE is_unrestricted), 0) AS unrestricted_bases
            FROM read_parquet('{cw_path}')
            GROUP BY 1, 2, 3, 4
            ORDER BY 1, 2, 3, 4
        """).fetchall()
        con.close()

    # Create persistent DuckDB sidecar pointing to the deterministic Parquet files via views
    sidecar_db = out / "env_sidecar.duckdb"
    if sidecar_db.exists():
        sidecar_db.unlink()
    scon = duckdb.connect(str(sidecar_db))
    scon.execute("SET threads=1")
    scon.execute(f"CREATE VIEW env_sample_dim AS SELECT * FROM read_parquet('{env_dim_path.resolve()}')")
    scon.execute(f"CREATE VIEW img_m_crosswalk AS SELECT * FROM read_parquet('{cw_path.resolve()}')")
    scon.close()

    # Write stratified ecosystem summary TSV
    strat_tsv = out / "img_m_by_ecosystem.tsv"
    with open(strat_tsv, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["harmonized_habitat", "eco_l1", "eco_l2", "decision", "units",
                    "unrestricted_units", "asm_bases", "unrestricted_bases"])
        for row in strat_rows:
            w.writerow(row)

    summary = {
        "img_m_total_units": int(summary_row[0]),
        "img_m_unrestricted_units": int(summary_row[1]),
        "img_m_non_restricted_units": int(summary_row[2]),
        "img_m_restricted_units": int(summary_row[3]),
        "img_m_novel_unrestricted_units": int(summary_row[4]),
        "img_m_novel_unrestricted_bases": int(summary_row[5]),
        "img_m_novel_non_restricted_units": int(summary_row[6]),
        "img_m_novel_non_restricted_bases": int(summary_row[7]),
        "img_m_in_spire_units": int(summary_row[8]),
        "img_m_in_mgnify_units": int(summary_row[9]),
        "img_m_in_insdc_units": int(summary_row[10]),
        "img_m_with_sra_run": int(summary_row[11]),
        "img_m_with_sra_exp": int(summary_row[12]),
        "img_m_with_biosample": int(summary_row[13]),
        "gold_sra_biome_rows": sra_biome_rows,
        "gold_eco_path_rows": eco_path_rows,
        "sha256_env_sample_dim": _sha256_file(env_dim_path),
        "sha256_img_m_crosswalk": _sha256_file(cw_path),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def build_stage_coverage(
    *,
    env_dir: Path | str,
    data_root: Path | str | None = None,
    staged_inputs_dir: Path | str | None = None,
    globdb_env_tsv: Path | str | None = None,
) -> dict[str, object]:
    """Compute environmental metadata coverage across stored stages and staged inputs."""
    out = Path(env_dir)
    env_dim_path = out / "env_sample_dim.parquet"
    if not env_dim_path.is_file():
        raise FileNotFoundError(f"missing {env_dim_path}; run build first")

    con = duckdb.connect(":memory:")
    con.execute("SET threads=1")
    con.execute("SET memory_limit='2GB'")
    con.execute("SET preserve_insertion_order=false")
    _register_udfs(con)

    con.execute(f"CREATE VIEW env_dim AS SELECT * FROM read_parquet('{env_dim_path}')")

    # Build lookup tables from env_sample_dim
    con.execute("""
        CREATE TABLE env_by_oid AS
        SELECT img_taxon_oid, min(env_uid) AS env_uid, min(source) AS env_source,
               min(harmonized_habitat) AS harmonized_habitat, min(ecosystem_path) AS ecosystem_path
        FROM env_dim WHERE img_taxon_oid <> '' GROUP BY img_taxon_oid
    """)
    # Assembly accession numeric core (matches GCA_000123456 and GCF_000123456)
    con.execute("""
        CREATE TABLE env_by_asm_core AS
        SELECT regexp_extract(insdc_assembly_acc, 'GC[AF]_(\\d+)', 1) AS asm_core,
               min(env_uid) AS env_uid, min(source) AS env_source,
               min(harmonized_habitat) AS harmonized_habitat, min(ecosystem_path) AS ecosystem_path
        FROM env_dim WHERE insdc_assembly_acc <> '' GROUP BY 1
    """)
    con.execute("""
        CREATE TABLE env_by_gb AS
        SELECT split_part(insdc_genbank_acc, '.', 1) AS gb_core,
               min(env_uid) AS env_uid, min(source) AS env_source,
               min(harmonized_habitat) AS harmonized_habitat, min(ecosystem_path) AS ecosystem_path
        FROM env_dim WHERE insdc_genbank_acc <> '' GROUP BY 1
    """)

    # Optional GlobDB environmental table (L1_* counts -> top habitat)
    globdb_map: dict[str, tuple[str, str]] = {}
    if globdb_env_tsv and Path(globdb_env_tsv).is_file():
        _L1_MAP = {
            "L1_animal_associated": "host-associated",
            "L1_built": "anthropogenic",
            "L1_soil": "terrestrial",
            "L1_sediment": "aquatic",
            "L1_aquatic": "aquatic",
        }
        with open(globdb_env_tsv, "rt", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for r in reader:
                gid = (r.get("ID") or "").strip()
                if not gid:
                    continue
                best_k, best_v = "", 0.0
                for col, hab in _L1_MAP.items():
                    v = _parse_float(r.get(col) or "") or 0.0
                    if v > best_v:
                        best_v = v
                        best_k = col
                if best_k:
                    globdb_map[gid] = (_L1_MAP[best_k], f"globdb:{best_k}")

    coverage: list[dict[str, object]] = []
    link_rows: list[dict[str, object]] = []

    if data_root and Path(data_root).is_dir():
        root = Path(data_root)
        # P1..P6: genome_uid -> l0_manifest.object_path (GCA/GCF)
        p_rows = con.execute(f"""
            WITH gi AS (
                SELECT stage, hex(genome_uid) AS ghex, min(contig_uid) AS cuid
                FROM read_parquet('{root}/l0seq/genome_input/*/*.parquet', hive_partitioning=1)
                WHERE stage IN ('P1','P2','P3','P4','P5','P6')
                GROUP BY stage, genome_uid
            ),
            co AS (
                SELECT stage, contig_uid, min(l0_object_sha256) AS osha
                FROM read_parquet('{root}/l0seq/contig_occurrence/*/*/*.parquet', hive_partitioning=1)
                WHERE stage IN ('P1','P2','P3','P4','P5','P6')
                GROUP BY stage, contig_uid
            ),
            lm AS (
                SELECT stage, object_sha256, min(object_path) AS opath
                FROM read_parquet('{root}/l0/l0_manifest/*/*.parquet', hive_partitioning=1)
                WHERE stage IN ('P1','P2','P3','P4','P5','P6')
                GROUP BY stage, object_sha256
            )
            SELECT
                gi.stage,
                gi.ghex,
                regexp_extract(lm.opath, '(GC[AF]_\\d+)', 1) AS asm_acc,
                e.env_uid,
                e.env_source,
                coalesce(e.harmonized_habitat, 'unlabelled') AS hab,
                coalesce(e.ecosystem_path, '') AS epath
            FROM gi
            JOIN co ON gi.stage = co.stage AND gi.cuid = co.contig_uid
            JOIN lm ON co.stage = lm.stage AND co.osha = lm.object_sha256
            LEFT JOIN env_by_asm_core e ON regexp_extract(lm.opath, 'GC[AF]_(\\d+)', 1) = e.asm_core
            ORDER BY gi.stage, gi.ghex
        """).fetchall()
        by_stage: dict[str, list[int]] = {}
        for st, ghex, asm_acc, euid, esrc, hab, epath in p_rows:
            s = by_stage.setdefault(st, [0, 0])
            s[0] += 1
            if euid is not None:
                s[1] += 1
                link_rows.append({
                    "stage": st, "genome_key": ghex, "source_accession": asm_acc,
                    "link_method": "asm_acc", "env_uid": bytes(euid), "env_source": esrc,
                    "harmonized_habitat": hab, "ecosystem_path": epath, "redistributable": True,
                })
        for st in sorted(by_stage):
            tot, hit = by_stage[st]
            coverage.append({
                "stage": st, "kind": "stored_prokaryote", "total_units": tot,
                "linked_units": hit, "pct_linked": round(100.0 * hit / tot, 2) if tot else 0.0,
                "method": "IMG/M+GOLD assembly_accession",
            })

        # S1, PL1, H1: contig_label.accession -> env_by_gb
        for st, kind in (("S1", "stored_virus"), ("PL1", "stored_plasmid"), ("H1", "stored_virus")):
            st_rows = con.execute(f"""
                SELECT
                    count(*) AS tot,
                    count(e.env_uid) AS hit
                FROM read_parquet('{root}/l0seq/contig_label/stage={st}/*.parquet') c
                LEFT JOIN env_by_gb e ON split_part(c.accession, '.', 1) = e.gb_core
            """).fetchone()
            tot, hit = int(st_rows[0]), int(st_rows[1])
            coverage.append({
                "stage": st, "kind": kind, "total_units": tot,
                "linked_units": hit, "pct_linked": round(100.0 * hit / tot, 2) if tot else 0.0,
                "method": "IMG/M+GOLD genbank_accession",
            })

    if staged_inputs_dir and Path(staged_inputs_dir).is_dir():
        sdir = Path(staged_inputs_dir)
        # V1: UHGV human gut virome
        v1_rec = sdir / "V1" / "records.tsv"
        if v1_rec.is_file():
            with open(v1_rec, "rt", encoding="utf-8") as fh:
                v1_tot = sum(1 for _ in fh) - 1
            coverage.append({
                "stage": "V1", "kind": "staged_virus", "total_units": v1_tot,
                "linked_units": v1_tot, "pct_linked": 100.0 if v1_tot else 0.0,
                "method": "source_provenance (UHGV human gut -> host-associated)",
            })

        # PL2: IMG-PR via IMG taxon OID + gut_plasmids_yu2024 provenance + INSDC GenBank
        pl2_rec = sdir / "PL2" / "records.tsv"
        if pl2_rec.is_file():
            oid_map = {
                r[0]: (bytes(r[1]), r[2], r[3], r[4])
                for r in con.execute("SELECT img_taxon_oid, env_uid, env_source, harmonized_habitat, ecosystem_path FROM env_by_oid").fetchall()
            }
            gb_set = {r[0] for r in con.execute("SELECT gb_core FROM env_by_gb").fetchall()}
            pl2_tot = 0
            pl2_img_hit = 0
            pl2_gut_hit = 0
            pl2_gb_hit = 0
            with open(pl2_rec, "rt", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for r in reader:
                    pl2_tot += 1
                    src = r.get("source") or ""
                    acc = r.get("accession") or ""
                    if src == "plasmidscope_IMG-PR":
                        m = _IMGPR_OID_RE.match(acc)
                        if m and m.group(1) in oid_map:
                            pl2_img_hit += 1
                    elif src == "gut_plasmids_yu2024":
                        pl2_gut_hit += 1
                    else:
                        core = acc.split("_", 1)[-1].split(".", 1)[0]
                        if core in gb_set:
                            pl2_gb_hit += 1
            pl2_hit = pl2_img_hit + pl2_gut_hit + pl2_gb_hit
            coverage.append({
                "stage": "PL2", "kind": "staged_plasmid", "total_units": pl2_tot,
                "linked_units": pl2_hit, "pct_linked": round(100.0 * pl2_hit / pl2_tot, 2) if pl2_tot else 0.0,
                "method": f"IMG-PR OID ({pl2_img_hit}) + gut provenance ({pl2_gut_hit}) + GenBank ({pl2_gb_hit})",
            })

        # G1: GlobDB r226 genomes via objects.json -> globdb high_level_env__value.tsv
        g1_obj = sdir / "G1" / "objects.json"
        if g1_obj.is_file():
            objs = json.loads(g1_obj.read_text(encoding="utf-8"))
            g1_tot = len(objs)
            g1_hit = 0
            for o in objs:
                fname = Path(o["path"]).name
                gid = fname.replace(".fa.gz", "").replace(".fna.gz", "")
                if gid in globdb_map:
                    g1_hit += 1
                    hab, epath = globdb_map[gid]
                    link_rows.append({
                        "stage": "G1", "genome_key": gid, "source_accession": gid,
                        "link_method": "globdb_env", "env_uid": None, "env_source": "globdb",
                        "harmonized_habitat": hab, "ecosystem_path": epath, "redistributable": True,
                    })
            coverage.append({
                "stage": "G1", "kind": "staged_prokaryote", "total_units": g1_tot,
                "linked_units": g1_hit, "pct_linked": round(100.0 * g1_hit / g1_tot, 2) if g1_tot else 0.0,
                "method": "GlobDB r226 high_level_env__value.tsv",
            })

    con.close()

    link_rows.sort(key=lambda r: (str(r["stage"]), str(r["genome_key"])))
    link_path = out / "genome_env_link.parquet"
    _write_arrow_chunks_to_parquet(iter(link_rows), GENOME_ENV_LINK_SCHEMA, link_path)

    sidecar_db = out / "env_sidecar.duckdb"
    if sidecar_db.exists():
        scon = duckdb.connect(str(sidecar_db))
        scon.execute("SET threads=1")
        scon.execute(f"CREATE OR REPLACE VIEW genome_env_link AS SELECT * FROM read_parquet('{link_path.resolve()}')")
        scon.execute("""
            CREATE OR REPLACE VIEW v_genome_env AS
            SELECT
                l.stage,
                l.genome_key,
                l.source_accession,
                l.link_method,
                l.harmonized_habitat,
                l.ecosystem_path,
                l.redistributable,
                d.eco_l1, d.eco_l2, d.eco_l3, d.eco_l4, d.eco_l5,
                d.ph, d.temperature_c, d.salinity, d.oxygen_req,
                d.depth_m, d.elevation_m, d.latitude, d.longitude, d.longhurst_code,
                d.host_name, d.body_site
            FROM genome_env_link l
            LEFT JOIN env_sample_dim d ON l.env_uid = d.env_uid
        """)
        scon.close()

    cov_tsv = out / "stage_coverage.tsv"
    with open(cov_tsv, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["stage", "kind", "total_units", "linked_units", "pct_linked", "method"])
        for row in coverage:
            w.writerow([row["stage"], row["kind"], row["total_units"], row["linked_units"],
                        row["pct_linked"], row["method"]])

    res = {"coverage": coverage, "genome_env_links": len(link_rows), "sha256_genome_env_link": _sha256_file(link_path)}
    (out / "coverage.json").write_text(json.dumps(res, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return res


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_b = sub.add_parser("build", help="Build GOLD + IMG/M crosswalk and DuckDB/Parquet sidecar")
    p_b.add_argument("--gold-data", type=Path, default=None)
    p_b.add_argument("--gold-sra-biome", type=Path, default=None)
    p_b.add_argument("--gold-eco-paths", type=Path, default=None)
    p_b.add_argument("--img-microbiomes", type=Path, required=True)
    p_b.add_argument("--img-extra", action="append", default=[], help="name=path for extra IMG/M slices")
    p_b.add_argument("--spire-dir", type=Path, default=None)
    p_b.add_argument("--c43-units", type=Path, default=None)
    p_b.add_argument("--out-dir", type=Path, required=True)

    p_c = sub.add_parser("coverage", help="Measure environmental coverage across stored stages & staged inputs")
    p_c.add_argument("--env-dir", type=Path, required=True)
    p_c.add_argument("--data-root", type=Path, default=None)
    p_c.add_argument("--staged-inputs-dir", type=Path, default=None)
    p_c.add_argument("--globdb-env-tsv", type=Path, default=None)

    args = ap.parse_args(argv)
    if args.cmd == "build":
        extras = []
        for item in args.img_extra:
            name, path = item.split("=", 1)
            extras.append((name, Path(path)))
        res = build_crosswalk(
            gold_data=args.gold_data,
            gold_sra_biome=args.gold_sra_biome,
            gold_eco_paths=args.gold_eco_paths,
            img_microbiomes=args.img_microbiomes,
            img_extra_slices=extras,
            spire_dir=args.spire_dir,
            c43_units=args.c43_units,
            out_dir=args.out_dir,
        )
        print(json.dumps(res, indent=2, sort_keys=True))
        return 0

    if args.cmd == "coverage":
        res = build_stage_coverage(
            env_dir=args.env_dir,
            data_root=args.data_root,
            staged_inputs_dir=args.staged_inputs_dir,
            globdb_env_tsv=args.globdb_env_tsv,
        )
        print(json.dumps(res, indent=2, sort_keys=True))
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
