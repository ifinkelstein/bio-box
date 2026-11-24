"""Card C-57 (DEC-69, DEC-65, DEC-50, DEC-51, DEC-26, DEC-27): stage curated integrative MGEs, transposons,
IS elements, integrons, and non-redundant prophages (`MGE1` redistributable public, `MGE2` non-redistributable GEMs).

    python c57_integrative_mges.py stage \
        --dl DIR --gems-tar TAR --out-mge1 DIR --out-mge2 DIR \
        --retrieved ISO [--root ROOT] [--root-stages S1,P1,...] \
        [--owner-inputs NAME=DIR ...] [--final-mge1 DIR] [--final-mge2 DIR] \
        [--out-report JSON]

Hard rules enforced:
- No 5 kb length floor: 0.7-3 kb IS elements, MITEs, and integron cassettes stay whole (DEC-69).
- Public curated MGEs are staged byte-identical as deposited (never trimmed, DEC-65).
- Redistributable split (DEC-51, DEC-63 (1)): refuses mixing `redistributable: true` and `redistributable: false`
  objects in one stage; splits cleanly into `MGE1` (public, `redistributable: true`, `kind="plasmid"`) and
  `MGE2` (`GEMs` VirSorter prophages, `redistributable: false`, `kind="virus"`).
- Within-stage deduplication by `contig_uid` (`DEC-27`): first record defines the genome in `genomes.json`,
  repeated sequences are logged in `repeated_records.tsv`.
- Cross-stage ownership (`DEC-26`): build-order `contig_uid` ownership against `--root-stages` (read-only
  `contig_index` default variants) and `--owner-inputs NAME=DIR` (`genomes.json`).
- Deterministic sorted TSV sidecar `mge_metadata.tsv` mapping each staged record to `mge_class`
  (`ICE`, `IME`, `CIME`, `Tn`, `IS`, `integron`, `prophage`), `family`, `host_accession`, and `owner_stage`.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import sys
import tarfile
import zipfile
from collections import Counter
from collections.abc import Iterable
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # c48c_virus_inputs and make_inputs beside this script

import c48c_virus_inputs as v2  # noqa: E402
import make_inputs  # noqa: E402
from biobox import fasta, ids  # noqa: E402

Refused = make_inputs.Refused

MGE_CLASSES = ("ICE", "IME", "CIME", "Tn", "IS", "integron", "prophage")
VALID_DNA_RE = re.compile(rb"^[ACGTURYKMSWBDHVNacgturykmswbdhvn]+$")
INSDC_ACC_RE = re.compile(r"\b((?:NC_|NZ_|AC_|AP_|CP|AE|AL|AM|BA|BX|CM|CR|CT|CU|FM|FN|FO|FP|FQ|FR|HE|HF|HG|LR|LS|LT|OU|OV|OW|OX|OY|OZ|[A-Z]{1,2}\d{5,8})(?:\.\d+)?)\b")

RECORD_HEAD = [
    "ordinal", "source", "accession", "contig_uid", "length", "circular",
    "mge_class", "family", "host_accession", "owner_stage",
]
SIDECAR_HEAD = [
    "stage", "source", "genome_id", "contig_uid", "accession", "length",
    "circular", "mge_class", "family", "host_accession", "owner_stage",
]
TRIM_HEAD = [
    "source", "genome", "accession", "source_contig_uid", "contig_uid",
    "source_length", "repeat_length", "bases_removed", "length",
]
SKIPPED_HEAD = ["source", "accession", "reason", "length"]


# --------------------------------------------------------------------------------------------------------------
# Metadata parsers


def parse_iceberg_header(header: str, default_class: str = "ICE") -> dict[str, str]:
    """Parse ICEberg 3.0 FASTA headers (`>ICEberg|id|name|GenBank|ACC|coords ...` or `>SRS... start-stop`)."""
    tok = header.split(None, 1)[0] if header else ""
    parts = tok.split("|")
    mge_class = default_class if default_class in MGE_CLASSES else "ICE"
    family = ""
    host_acc = ""

    if len(parts) >= 4 and parts[0].lower() == "iceberg":
        # e.g. ICEberg|1010|ICEPae690|GenBank|KY852375|315..86517
        # or   ICEberg|1234|ICE_SXT|NZ_CP012345.1:100..50000
        tag = parts[1].upper()
        if "CIME" in tag:
            mge_class = "CIME"
        elif "IME" in tag:
            mge_class = "IME"
        family = parts[2]
        if len(parts) >= 5 and parts[3].lower() in ("genbank", "refseq", "insdc", "ncbi"):
            host_acc = parts[4].split(":", 1)[0]
        else:
            host_acc = parts[3].split(":", 1)[0]
    else:
        family = "Oral_mICE" if "oral" in default_class.lower() else default_class
        host_acc = tok

    if not host_acc:
        m = INSDC_ACC_RE.search(header)
        if m:
            host_acc = m.group(1)

    return {
        "accession": tok,
        "mge_class": mge_class,
        "family": family or mge_class,
        "host_accession": host_acc,
    }


def load_isfinder_csv(path: Path) -> dict[str, dict[str, str]]:
    """Load ISfinder-sequences `IS.csv` keyed by IS Name."""
    out: dict[str, dict[str, str]] = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            name = (row.get("Name") or "").strip()
            if not name:
                continue
            out[name] = {
                "family": (row.get("Family") or "").strip(),
                "group": (row.get("Group") or "").strip(),
                "origin": (row.get("Origin") or "").strip(),
                "accession": (row.get("Accession Number") or "").strip(),
            }
    return out


def parse_isfinder_header(header: str, is_meta: dict[str, dict[str, str]] | None = None) -> dict[str, str]:
    """Parse ISfinder-sequences `IS.fna` header (`>ISname_Family_Group` or `>ISname_Group_Family`)."""
    tok = header.split(None, 1)[0] if header else ""
    parts = tok.split("_")
    name = parts[0] if parts else tok
    meta = (is_meta or {}).get(name) or {}

    family = meta.get("family") or ""
    if not family:
        # In ISfinder-sequences `IS.fna`, headers are `>IS100_IS21_unknown` (middle is family, last is group)
        # or `>IS605_IS605_IS200/IS605` (last has slash). Pick the most informative non-unknown token.
        if len(parts) >= 3:
            if "/" in parts[-1] or (parts[1] == name and parts[-1].lower() != "unknown"):
                family = parts[-1]
            elif parts[1].lower() not in ("unknown", "na", ""):
                family = parts[1]
            elif parts[-1].lower() not in ("unknown", "na", ""):
                family = parts[-1]
        elif len(parts) == 2:
            family = parts[1]

    host_acc = meta.get("accession") or ""
    if host_acc.upper() in ("NA", "NONE", "UNKNOWN"):
        host_acc = ""

    return {
        "accession": tok,
        "mge_class": "IS",
        "family": family or "IS_unclassified",
        "host_accession": host_acc,
    }


def _map_mgedb_type(raw_type: str) -> str:
    t = raw_type.strip().lower()
    if t in ("is", "insertion sequence", "mite", "mic"):
        return "IS"
    if t in ("tn", "cn", "unit transposon", "composite transposon", "transposon"):
        return "Tn"
    if t == "ice":
        return "ICE"
    if t == "ime":
        return "IME"
    if t == "cime":
        return "CIME"
    if "integron" in t:
        return "integron"
    return "Tn"


def load_mgedb_json(raw_bytes: bytes) -> dict[str, dict[str, str]]:
    """Parse MGEdb `mge_records.json` into name -> {type, family, accession}."""
    data = json.loads(raw_bytes.decode("utf-8", "replace"))
    out: dict[str, dict[str, str]] = {}
    items = data.values() if isinstance(data, dict) else data
    for entry in items:
        if not isinstance(entry, dict):
            continue
        name = str(entry.get("name") or "").strip()
        if not name:
            continue
        acc = ""
        seqs = entry.get("sequences") or []
        if seqs and isinstance(seqs, list) and isinstance(seqs[0], dict):
            acc = str(seqs[0].get("accession") or "").strip()
        out[name] = {
            "type": str(entry.get("type") or "").strip(),
            "family": str(entry.get("family") or "").strip(),
            "accession": acc,
        }
    return out


def parse_mgedb_record(header: str, meta_row: dict[str, str] | None = None) -> dict[str, str]:
    """Parse MGEdb `mge_records.fna` header (`>name|seq_idx|accession`) + metadata row."""
    tok = header.split(None, 1)[0] if header else ""
    parts = tok.split("|")
    name = parts[0] if parts else tok
    host_acc = parts[2] if len(parts) >= 3 else ""
    row = meta_row or {}
    mge_class = _map_mgedb_type(row.get("type") or "")
    family = row.get("family") or name
    if not host_acc:
        host_acc = row.get("accession") or ""
    return {
        "accession": tok,
        "mge_class": mge_class,
        "family": family if family.lower() != "none" else name,
        "host_accession": host_acc,
    }


def parse_tncentral_header(header: str, default_class: str = "Tn") -> dict[str, str]:
    """Parse PlasAnn `tncentral_cleaned.fa` (`>IS1006-NC_004361`, `>Tn7-KX117211`) or `transposon.fasta` (`>Tn1`)."""
    tok = header.split(None, 1)[0] if header else ""
    name = tok
    host_acc = ""
    if "-" in tok:
        left, right = tok.split("-", 1)
        if INSDC_ACC_RE.fullmatch(right):
            name, host_acc = left, right
    low = name.lower()
    if low.startswith("is") or low.startswith("mite"):
        mge_class = "IS"
    elif low.startswith("in") and not low.startswith("int") and len(name) > 2 and name[2].isdigit():
        mge_class = "integron"
    elif low.startswith("inti") or "integron" in low:
        mge_class = "integron"
    elif low.startswith("ice"):
        mge_class = "ICE"
    elif low.startswith("ime"):
        mge_class = "IME"
    elif low.startswith("cime"):
        mge_class = "CIME"
    elif low.startswith("tn") or low.startswith("ctn"):
        mge_class = "Tn"
    else:
        mge_class = default_class if default_class in MGE_CLASSES else "Tn"

    return {
        "accession": tok,
        "mge_class": mge_class,
        "family": name or mge_class,
        "host_accession": host_acc,
    }


def load_parnanen_tax(path: Path) -> dict[str, tuple[str, str]]:
    """Load Parnanen `MGE_tax_table_trim.txt` (record_id -> (raw_class, gene_family))."""
    out: dict[str, tuple[str, str]] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            parts = ln.rstrip("\r\n").split("\t")
            if len(parts) >= 2:
                out[parts[0].strip()] = (parts[1].strip(), parts[2].strip() if len(parts) >= 3 else parts[1].strip())
    return out


def parse_parnanen_record(header: str, tax_row: tuple[str, str] | None = None, default_class: str = "integron") -> dict[str, str] | None:
    """Parse Parnanen MGEdb record. Returns None for pure `plasmid` replicon rows so MGE1 stays integrative MGEs."""
    tok = header.split(None, 1)[0] if header else ""
    raw_cls, fam = tax_row if tax_row else ("", "")
    low_cls = raw_cls.lower()
    low_tok = tok.lower()
    if low_cls == "plasmid":
        return None
    if low_cls in ("integrase", "qacedelta") or "inti" in low_tok or "integron" in low_tok:
        mge_class = "integron"
    elif low_cls.startswith("insertion_element") or low_cls.startswith("is") or low_cls.startswith("ist"):
        mge_class = "IS"
    elif low_cls in ("transposase", "tn916", "tnia", "tnib", "tniab", "tnp-iscr") or low_cls.startswith("tn"):
        mge_class = "Tn"
    else:
        mge_class = default_class if default_class in MGE_CLASSES else "Tn"

    host_acc = ""
    m = INSDC_ACC_RE.search(tok.replace("_", " "))
    if m:
        host_acc = m.group(1)

    return {
        "accession": tok,
        "mge_class": mge_class,
        "family": fam or raw_cls or mge_class,
        "host_accession": host_acc,
    }


def load_gems_virsorter_tsv(fh) -> dict[str, dict[str, str]]:
    """Load GEMs `virsorter.tsv` keyed by `prophage_id`."""
    out: dict[str, dict[str, str]] = {}
    text_fh = (ln.decode("utf-8", "replace") if isinstance(ln, bytes) else ln for ln in fh)
    for row in csv.DictReader(text_fh, delimiter="\t"):
        pid = (row.get("prophage_id") or "").strip()
        if pid:
            out[pid] = row
    return out


def parse_gems_prophage(header: str, virsorter_row: dict[str, str] | None = None) -> dict[str, str]:
    """Parse GEMs `virsorter.fna` header + `virsorter.tsv` row."""
    tok = header.split(None, 1)[0] if header else ""
    row = virsorter_row or {}
    return {
        "accession": tok,
        "mge_class": "prophage",
        "family": "VirSorter",
        "host_accession": (row.get("contig_id") or row.get("bin_id") or "").strip(),
    }


# --------------------------------------------------------------------------------------------------------------
# Staging core


def _open_fasta_in(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return open(path, "rb")


def _parse_record_meta(src: dict, header: str) -> dict[str, str] | None:
    parser = src.get("parser", "generic")
    default_class = src.get("mge_class", "Tn")
    if parser == "iceberg":
        return parse_iceberg_header(header, default_class)
    if parser == "isfinder":
        return parse_isfinder_header(header, src.get("meta_map"))
    if parser == "mgedb":
        tok = header.split(None, 1)[0] if header else ""
        name = tok.split("|", 1)[0]
        return parse_mgedb_record(header, (src.get("meta_map") or {}).get(name))
    if parser == "tncentral":
        return parse_tncentral_header(header, default_class)
    if parser == "parnanen":
        tok = header.split(None, 1)[0] if header else ""
        return parse_parnanen_record(header, (src.get("meta_map") or {}).get(tok), default_class)
    if parser == "gems_prophage":
        tok = header.split(None, 1)[0] if header else ""
        return parse_gems_prophage(header, (src.get("meta_map") or {}).get(tok))
    tok = header.split(None, 1)[0] if header else ""
    return {
        "accession": tok,
        "mge_class": default_class if default_class in MGE_CLASSES else "Tn",
        "family": default_class,
        "host_accession": "",
    }


def _stage_one(
    sources: list[dict],
    out_dir: Path,
    *,
    stage_name: str,
    root: Path | None,
    root_stages: list[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
    final_dir: Path | None = None,
) -> dict:
    if not sources:
        raise Refused(f"{stage_name}: no sources to stage")
    out_dir = Path(out_dir)
    if out_dir.exists():
        raise Refused(f"{out_dir} exists")
    redists = {bool(s["redistributable"]) for s in sources}
    if len(redists) > 1:
        raise Refused(
            f"{stage_name}: refusing to mix redistributable: true and redistributable: false sources in one stage "
            "(DEC-51 / DEC-63 (1)); pass split_by_redistributable=True"
        )

    work = out_dir.with_name(f".{out_dir.name}.tmp")
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    obj_dir = Path(final_dir) if final_dir is not None else out_dir

    try:
        rep = _stage_one_work(
            sources, work, obj_dir, stage_name, root, root_stages, owner_inputs, retrieved
        )
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise

    work.rename(out_dir)
    sum_files = sorted(
        p.name for p in out_dir.iterdir()
        if p.is_file() and p.name not in ("inputs_SUMS", "stage_report.json")
    )
    sha_map = {n: ids.sha256_file(out_dir / n).hex() for n in sum_files}
    bytes_map = {n: (out_dir / n).stat().st_size for n in sum_files}
    with open(out_dir / "inputs_SUMS", "w") as w:
        for n in sum_files:
            w.write(f"{sha_map[n]}  {n}\n")
    rep["sha256"] = sha_map
    rep["bytes"] = bytes_map
    (out_dir / "stage_report.json").write_text(json.dumps(rep, indent=1, sort_keys=True))
    return rep


def _stage_one_work(
    sources: list[dict],
    work: Path,
    obj_dir: Path,
    stage_name: str,
    root: Path | None,
    root_stages: list[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
) -> dict:
    # Pass 1: stream each source FASTA, validate DNA (no 5 kb floor!), write byte-identical staged object
    # (for sources where every record is valid and none skipped, write clean FASTA preserving exact sequence bytes).
    staged_records: list[tuple[str, str, bytes, int, bool, str, str, str, str]] = []
    # (source_name, acc, cu, length, circ, kind, mge_class, family, host_acc)
    skipped: list[tuple[str, str, str, int]] = []
    trims: list[tuple] = []
    objects_json: list[dict] = []
    by_source: dict[str, dict] = {}

    for src in sources:
        sname = src["source_name"]
        out_name = src["output_name"]
        fa_path = Path(src["fasta"])
        if not fa_path.is_file():
            raise Refused(f"{sname}: FASTA file not found: {fa_path}")
        kind = src.get("kind", "plasmid")
        dst_fa = work / out_name
        rec_count = 0
        base_count = 0
        skip_count = 0

        with _open_fasta_in(fa_path) as fh, open(dst_fa, "wb") as w:
            for r in v2.iter_fasta_stream(fh):
                meta = _parse_record_meta(src, r.header)
                if meta is None:
                    skipped.append((sname, r.accession, "excluded_class", len(r.seq)))
                    skip_count += 1
                    continue
                if len(r.seq) == 0:
                    skipped.append((sname, r.accession, "empty_sequence", 0))
                    skip_count += 1
                    continue
                if not VALID_DNA_RE.match(r.seq):
                    skipped.append((sname, r.accession, "non_iupac_dna", len(r.seq)))
                    skip_count += 1
                    continue
                # Public curated MGEs & excised prophages: byte-identical, never trimmed (DEC-65)
                seq = r.seq
                circ = False
                w.write(b"\n".join(fasta.format_record(r.header, seq, r.line_width)) + b"\n")
                cu = ids.contig_uid(seq)
                staged_records.append((
                    sname, meta["accession"], cu, len(seq), circ, kind,
                    meta["mge_class"], meta["family"], meta["host_accession"],
                ))
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
            "redistributable": bool(src["redistributable"]),
            "retrieved_at": retrieved,
        })
        by_source[out_name] = {
            "source_name": sname,
            "records": rec_count,
            "bases": base_count,
            "skipped": skip_count,
        }

    if not staged_records:
        raise Refused(f"{stage_name}: no valid records staged")

    # Compute cross-stage build-order ownership (DEC-26)
    all_uids = [x[2] for x in staged_records]
    owners = v2.owners_of(all_uids, root, root_stages, owner_inputs)

    # Deduplicate within stage (DEC-27): first record of a contig_uid defines the genome
    listing = make_inputs.Listing()
    genomes: list[dict] = []
    rec_rows: list[tuple] = []
    sidecar_rows: list[tuple] = []
    repeated_rows: list[tuple] = []
    shared_rows: list[tuple] = []
    first_of: dict[bytes, str] = {}
    class_counts: Counter[str] = Counter()
    genome_class_counts: Counter[str] = Counter()
    len_buckets: Counter[str] = Counter()

    for ordinal, (sname, acc, cu, length, circ, kind, mge_cls, fam, host_acc) in enumerate(staged_records):
        owner = owners.get(cu, "")
        rec_rows.append((
            ordinal, sname, acc, cu.hex(), length, str(circ).lower(),
            mge_cls, fam, host_acc, owner,
        ))
        if owner:
            shared_rows.append((ordinal, acc, cu.hex(), owner))
        class_counts[mge_cls] += 1

        g = listing.add([cu], 11, kind, "contig")
        guid = ids.genome_uid([cu])
        sidecar_rows.append((
            stage_name, sname, guid.hex(), cu.hex(), acc, length,
            str(circ).lower(), mge_cls, fam, host_acc, owner,
        ))
        if g is None:
            repeated_rows.append((ordinal, acc, guid.hex(), first_of[guid]))
            continue
        g["circular"] = [cu.hex()] if circ else []
        genomes.append(g)
        first_of[guid] = acc
        genome_class_counts[mge_cls] += 1
        if length < 1000:
            len_buckets["<1kb"] += 1
        elif length < 5000:
            len_buckets["1-5kb"] += 1
        elif length < 20000:
            len_buckets["5-20kb"] += 1
        elif length < 100000:
            len_buckets["20-100kb"] += 1
        else:
            len_buckets[">=100kb"] += 1

    make_inputs._write_tsv(work / "records.tsv", RECORD_HEAD, rec_rows)
    make_inputs._write_tsv(work / "mge_metadata.tsv", SIDECAR_HEAD, sidecar_rows)
    make_inputs._write_tsv(work / "repeated_records.tsv", ["ordinal", "accession", "genome_uid", "first_accession"], repeated_rows)
    make_inputs._write_tsv(work / "shared_records.tsv", ["ordinal", "accession", "contig_uid", "owner_stage"], shared_rows)
    make_inputs._write_tsv(work / "trims.tsv", TRIM_HEAD, trims)
    make_inputs._write_tsv(work / "skipped_records.tsv", SKIPPED_HEAD, skipped)

    (work / "objects.json").write_text(json.dumps(objects_json, indent=1))
    (work / "genomes.json").write_text(json.dumps(genomes))

    contig_len: dict[bytes, int] = {}
    for _, _, cu, length, *_ in staged_records:
        contig_len.setdefault(cu, length)

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
        "redistributable": bool(sources[0]["redistributable"]),
        "records": len(staged_records),
        "skipped_records": len(skipped),
        "genomes": len(genomes),
        "contigs": total_contigs,
        "bases": total_bases,
        "repeated_records": len(repeated_rows),
        "trims": len(trims),
        "owned_by_stage": dict(sorted(owner_contigs.items())),
        "owned_bases_by_stage": dict(sorted(owner_bases.items())),
        "held_contigs": held_contigs,
        "held_bases": held_bases,
        "new_contigs": total_contigs - held_contigs,
        "new_bases": total_bases - held_bases,
        "records_by_class": dict(sorted(class_counts.items())),
        "genomes_by_class": dict(sorted(genome_class_counts.items())),
        "genomes_by_length_bucket": dict(sorted(len_buckets.items())),
        "by_source": by_source,
    }


def stage_mges(
    sources: list[dict],
    out_dir: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
    stage_name: str = "MGE1",
    final_dir: Path | None = None,
    split_by_redistributable: bool = False,
    out_restricted: Path | None = None,
    final_restricted: Path | None = None,
    restricted_stage_name: str = "MGE2",
) -> dict[str, dict]:
    """Stage MGE sources into `MGE1` (and `MGE2` if `split_by_redistributable=True`)."""
    rs_list = list(root_stages)
    if not split_by_redistributable:
        rep = _stage_one(
            sources, out_dir, stage_name=stage_name, root=root, root_stages=rs_list,
            owner_inputs=owner_inputs, retrieved=retrieved, final_dir=final_dir,
        )
        return {stage_name: rep}

    pub = [s for s in sources if s["redistributable"]]
    res = [s for s in sources if not s["redistributable"]]
    out: dict[str, dict] = {}
    if pub:
        out[stage_name] = _stage_one(
            pub, out_dir, stage_name=stage_name, root=root, root_stages=rs_list,
            owner_inputs=owner_inputs, retrieved=retrieved, final_dir=final_dir,
        )
    if res:
        if out_restricted is None:
            raise Refused("split_by_redistributable=True with non-redistributable sources requires out_restricted")
        oi2 = dict(owner_inputs)
        if pub:
            oi2[stage_name] = Path(out_dir)
        out[restricted_stage_name] = _stage_one(
            res, out_restricted, stage_name=restricted_stage_name, root=root, root_stages=rs_list,
            owner_inputs=oi2, retrieved=retrieved, final_dir=final_restricted,
        )
    return out


# --------------------------------------------------------------------------------------------------------------
# CLI driver for pod staging


def _extract_archives(dl_dir: Path, gems_tar: Path | None, scratch: Path) -> list[dict]:
    """Prepare source dicts from downloaded archives in `dl_dir` and optional on-pod `gems_tar`."""
    scratch.mkdir(parents=True, exist_ok=True)
    sources: list[dict] = []

    # 1. ICEberg 3.0
    ice_rel = "3.0 (2025-03/06; Wang et al. 2024, NAR 52:D732-D737, doi:10.1093/nar/gkad983)"
    ice_lic = "open access academic database (Wang et al. 2024, NAR, CC BY 4.0); INSDC-derived"
    for fname, mcls, out_name in (
        ("ICE_seq_all.fas", "ICE", "iceberg3_ICE.fna"),
        ("AICE_seq_all.fas", "ICE", "iceberg3_AICE.fna"),
        ("IME_seq_all.fas", "IME", "iceberg3_IME.fna"),
        ("CIME_seq_all.fas", "CIME", "iceberg3_CIME.fna"),
        ("Oral_ICE_seq.fas", "ICE", "iceberg3_Oral_mICE.fna"),
    ):
        p = dl_dir / "iceberg3" / fname
        if p.is_file():
            sources.append({
                "source_name": "iceberg3",
                "source_release": ice_rel,
                "license": ice_lic,
                "redistributable": True,
                "fasta": p,
                "output_name": out_name,
                "parser": "iceberg",
                "kind": "plasmid",
                "mge_class": mcls,
            })

    # 2. ISfinder-sequences
    is_fa = dl_dir / "isfinder_sequences" / "IS.fna"
    is_csv = dl_dir / "isfinder_sequences" / "IS.csv"
    if is_fa.is_file():
        is_meta = load_isfinder_csv(is_csv) if is_csv.is_file() else {}
        sources.append({
            "source_name": "isfinder_sequences",
            "source_release": "commit bedbeb5c (2020-10-06; ISfinder, Siguier et al. 2006, NAR 34:D32-D36)",
            "license": "open community compilation from ISfinder (cite Siguier et al. 2006, doi:10.1093/nar/gkj014)",
            "redistributable": True,
            "fasta": is_fa,
            "output_name": "isfinder_IS.fna",
            "parser": "isfinder",
            "meta_map": is_meta,
            "kind": "plasmid",
            "mge_class": "IS",
        })

    # 3. MGEdb 1.1.1 (MobileElementFinder)
    mgedb_tar = dl_dir / "mgedb" / "MGEdb-1.1.1.tar.gz"
    if mgedb_tar.is_file():
        mge_fna = scratch / "mgedb_records.fna"
        mge_meta: dict[str, dict[str, str]] = {}
        with tarfile.open(mgedb_tar, "r:gz") as tf:
            for member in tf.getmembers():
                if member.name.endswith("mge_records.fna"):
                    with tf.extractfile(member) as src_fh, open(mge_fna, "wb") as dst_fh:
                        shutil.copyfileobj(src_fh, dst_fh)
                elif member.name.endswith("mge_records.json"):
                    mge_meta = load_mgedb_json(tf.extractfile(member).read())
        if mge_fna.is_file():
            sources.append({
                "source_name": "mgedb_mobileelementfinder",
                "source_release": "v1.1.1 (Johansson et al. 2021, JAC 76:101-109, doi:10.1093/jac/dkaa390)",
                "license": "GPLv3",
                "redistributable": True,
                "fasta": mge_fna,
                "output_name": "mgedb_1.1.1.fna",
                "parser": "mgedb",
                "meta_map": mge_meta,
                "kind": "plasmid",
                "mge_class": "Tn",
            })

    # 4. PlasAnn 1.1.6 (TnCentral cleaned + transposon.fasta)
    plasann_zip = dl_dir / "plasann" / "plasann_databases.zip"
    if plasann_zip.is_file():
        with zipfile.ZipFile(plasann_zip) as zf:
            for zmember, out_name in (
                ("Database/tncentral_cleaned.fa", "plasann_tncentral.fna"),
                ("Database/transposon.fasta", "plasann_transposon.fna"),
            ):
                if zmember in zf.namelist():
                    dst = scratch / out_name
                    with zf.open(zmember) as src_fh, open(dst, "wb") as dst_fh:
                        shutil.copyfileobj(src_fh, dst_fh)
                    sources.append({
                        "source_name": "plasann_tncentral",
                        "source_release": "v1.1.6 (Zenodo 10.5281/zenodo.17886383; TnCentral)",
                        "license": "CC BY 4.0",
                        "redistributable": True,
                        "fasta": dst,
                        "output_name": out_name,
                        "parser": "tncentral",
                        "kind": "plasmid",
                        "mge_class": "Tn",
                    })

    # 5. Parnanen MobileGeneticElementDatabase (integrons, transposases, IS)
    parn_tar = dl_dir / "parnanen" / "MGEs_FINAL_99perc_trim.fasta.tar.gz"
    parn_tax = dl_dir / "parnanen" / "MGE_tax_table_trim.txt"
    if parn_tar.is_file():
        parn_fna = scratch / "parnanen_mges.fna"
        with tarfile.open(parn_tar, "r:gz") as tf:
            for member in tf.getmembers():
                if member.isfile() and member.name.endswith(".fasta"):
                    with tf.extractfile(member) as src_fh, open(parn_fna, "wb") as dst_fh:
                        shutil.copyfileobj(src_fh, dst_fh)
                    break
        tax_map = load_parnanen_tax(parn_tax) if parn_tax.is_file() else {}
        if parn_fna.is_file():
            sources.append({
                "source_name": "parnanen_mgedb",
                "source_release": "commit 718f4570 (2021-08-16; Parnanen et al. 2018, Nat. Commun., doi:10.1038/s41467-018-06393-w)",
                "license": "open access academic compilation (Parnanen et al. 2018, Nat. Commun., CC BY 4.0)",
                "redistributable": True,
                "fasta": parn_fna,
                "output_name": "parnanen_mges.fna",
                "parser": "parnanen",
                "meta_map": tax_map,
                "kind": "plasmid",
                "mge_class": "integron",
            })

    # 6. GEMs VirSorter prophages (non-redistributable, on-pod tar)
    if gems_tar is not None and gems_tar.is_file():
        gems_fna = scratch / "gems_virsorter.fna"
        gems_meta: dict[str, dict[str, str]] = {}
        with tarfile.open(gems_tar, "r") as tf:
            for member in tf.getmembers():
                if member.name.endswith("prophages/virsorter.fna"):
                    with tf.extractfile(member) as src_fh, open(gems_fna, "wb") as dst_fh:
                        shutil.copyfileobj(src_fh, dst_fh)
                elif member.name.endswith("prophages/virsorter.tsv"):
                    gems_meta = load_gems_virsorter_tsv(tf.extractfile(member))
        if gems_fna.is_file():
            sources.append({
                "source_name": "jgi_gems_nersc",
                "source_release": "GEMs (Nayfach et al. 2021, Nat. Biotechnol.; virsorter.fna 24,313 prophages)",
                "license": "JGI Data Utilization Policy (GEM README.md)",
                "redistributable": False,
                "fasta": gems_fna,
                "output_name": "gems_virsorter_prophages.fna",
                "parser": "gems_prophage",
                "meta_map": gems_meta,
                "kind": "virus",
                "mge_class": "prophage",
            })

    return sources


def _parse_kv(items: list[str]) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for item in items:
        if "=" not in item:
            raise Refused(f"--owner-inputs expects NAME=DIR, got {item!r}")
        k, v = item.split("=", 1)
        out[k] = Path(v)
    return out


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description="Card C-57: stage integrative MGEs (MGE1) and GEMs prophages (MGE2).")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sp = sub.add_parser("stage")
    sp.add_argument("--dl", type=Path, required=True)
    sp.add_argument("--gems-tar", type=Path, default=None)
    sp.add_argument("--scratch", type=Path, required=True)
    sp.add_argument("--out-mge1", type=Path, required=True)
    sp.add_argument("--out-mge2", type=Path, default=None)
    sp.add_argument("--final-mge1", type=Path, default=None)
    sp.add_argument("--final-mge2", type=Path, default=None)
    sp.add_argument("--retrieved", required=True)
    sp.add_argument("--root", type=Path, default=None)
    sp.add_argument("--root-stages", default="")
    sp.add_argument("--owner-inputs", action="append", default=[])
    sp.add_argument("--out-report", type=Path, default=None)

    args = ap.parse_args(argv)
    root_stages = [s for s in args.root_stages.split(",") if s]
    owner_inputs = _parse_kv(args.owner_inputs)

    sources = _extract_archives(args.dl, args.gems_tar, args.scratch)
    rep = stage_mges(
        sources,
        args.out_mge1,
        root=args.root,
        root_stages=root_stages,
        owner_inputs=owner_inputs,
        retrieved=args.retrieved,
        stage_name="MGE1",
        final_dir=args.final_mge1,
        split_by_redistributable=(args.out_mge2 is not None),
        out_restricted=args.out_mge2,
        final_restricted=args.final_mge2,
        restricted_stage_name="MGE2",
    )
    if args.out_report is not None:
        args.out_report.parent.mkdir(parents=True, exist_ok=True)
        args.out_report.write_text(json.dumps(rep, indent=1, sort_keys=True))
    print(json.dumps(rep, indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except Refused as exc:
        print(f"refused: {exc}", file=sys.stderr)
        sys.exit(2)
