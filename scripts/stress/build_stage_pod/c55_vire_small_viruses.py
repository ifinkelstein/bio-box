"""Card C-55 (DEC-69, DEC-71, DEC-65, DEC-59, DEC-50, DEC-51, DEC-26, DEC-27):
survey, fetch, deduplicate, and stage:
  - `SAT1`: compact public `redistributable: true` viral/satellite stage on `/tmp/if2453/s2/inputs/SAT1/`
    (post-14Apr2025 INPHARED 9Aug2025 + GenBank phage delta incl. PX530545, curated phage satellites
    PICIs/cf-PICIs/SaPIs/PLEs/P4, < 5 kb circular DNA viruses Cressdnaviricota/Microviridae/Tisza, and
    archaeal extremophile viruses; byte-identical as deposited under DEC-65 (2); no 5 kb floor; kind="virus").
  - `V3`: EMBL VIRE (`all_samples_genomes.tar`) 693,656 non-Riboviria DNA vOTU species representatives
    streamed directly into `/stor/scratch/Wilke/if2453/bio-box-staging/V3/vire_v3_dna_votus.fna` with DEC-65 (1)
    exact 20..1,000 bp DTR trimming and metadata in `/tmp/if2453/c55/V3/`.
"""
from __future__ import annotations

import argparse
import csv
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
import urllib.parse
import urllib.request
from collections import Counter
from collections.abc import Iterable, Iterator
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # c48c_virus_inputs and make_inputs beside this script

import c48c_virus_inputs as v2  # noqa: E402
import make_inputs  # noqa: E402
from biobox import fasta, ids  # noqa: E402

csv.field_size_limit(1 << 30)
Refused = make_inputs.Refused

VALID_DNA_RE = re.compile(rb"^[ACGTURYKMSWBDHVNacgturykmswbdhvn]+$")

ARCHAEAL_VIRUS_TAXA = (
    "Adnaviria",
    "Fuselloviridae",
    "Lipothrixviridae",
    "Rudiviridae",
    "Bicaudaviridae",
    "Ampullaviridae",
    "Globuloviridae",
    "Guttaviridae",
    "Pleolipoviridae",
    "Sphaerolipoviridae",
    "Turriviridae",
    "Tristromaviridae",
    "Portogloboviridae",
    "Ovaliviridae",
    "Halspiviridae",
    "Thaspiviridae",
    "Spiraviridae",
    "Clavaviridae",
    "Halomagnusviridae",
    "Pungoviridae",
    "Ungulaviridae",
    "Ahmunviridae",
    "Apasviridae",
    "Hafunviridae",
    "Halomagnusviridae",
    "Saparoviridae",
    "Suoliviridae",
    "Verdandiviridae",
    "Vertoviridae",
    "Yangangviridae",
    "Yanodviridae",
)

ARCHAEAL_HOST_TOKENS = (
    "archaea",
    "sulfolobus",
    "saccharolobus",
    "acidianus",
    "metallosphaera",
    "pyrobaculum",
    "thermoproteus",
    "pyrococcus",
    "thermococcus",
    "aeropyrum",
    "ignicoccus",
    "methanobacterium",
    "methanothermobacter",
    "methanococcus",
    "methanosarcina",
    "halobacterium",
    "haloferax",
    "haloarcula",
    "halorubrum",
    "natrialba",
    "natrinema",
    "halorhabdus",
    "halogranum",
    "haloquadratum",
    "halogeometricum",
    "haloterrigena",
)

SATELLITE_TOKENS = (
    "phage satellite",
    "satellite phage",
    "phage-inducible chromosomal island",
    "pathogenicity island sapi",
    "pici-like element",
    "enterobacteria phage p4",
    "escherichia phage p4",
    "bacteriophage p4",
    "capsid-forming pici",
    "vibrio phage ple",
)

SAT1_RECORD_HEAD = [
    "ordinal", "source", "category", "accession", "contig_uid", "length", "circular", "owner_stage",
]
SAT1_META_HEAD = [
    "stage", "source", "category", "genome_id", "contig_uid", "accession", "length", "circular", "owner_stage",
]
V3_TRIM_HEAD = [
    "source", "genome", "accession", "source_contig_uid", "contig_uid",
    "source_length", "repeat_length", "bases_removed", "length",
]
V3_RECORD_HEAD = [
    "ordinal", "genome", "species_cluster", "realm", "contig_uid", "length", "dtr_class",
    "circular", "source_contig_uid", "source_length", "trimmed_bases", "alignment_to_imgvr",
]
V3_LEFT_HEAD = [
    "genome_id", "biosample_name", "species_cluster", "realm", "contig_length",
    "checkv_checkv_quality", "checkv_completeness", "checkv_contamination",
    "genomad_topology", "alignment_to_imgvr", "reason",
]


def _base_acc(acc: str) -> str:
    acc = acc.strip()
    if "." in acc:
        return acc.split(".", 1)[0]
    return acc


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _open_bin(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return open(path, "rb")


def _read_inphared_accessions(tsv_path: Path) -> set[str]:
    accs: set[str] = set()
    with _open_text(tsv_path) as fh:
        header = fh.readline()
        if not header:
            return accs
        cols = header.rstrip("\r\n").split("\t")
        idx = 0
        for i, c in enumerate(cols):
            if c.strip().lower() == "accession":
                idx = i
                break
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\r\n").split("\t")
            if idx < len(parts) and parts[idx].strip():
                accs.add(_base_acc(parts[idx]))
    return accs


def filter_inphared_delta(new_tsv: Path, s1_tsv: Path) -> set[str]:
    """Return base accessions present in `new_tsv` but absent from `s1_tsv`."""
    s1_accs = _read_inphared_accessions(Path(s1_tsv))
    new_accs = _read_inphared_accessions(Path(new_tsv))
    return new_accs - s1_accs


def parse_genbank_dna_records(gb_fh, *, allow_rna: bool = False) -> Iterator[dict]:
    """Stream GenBank flatfile records, yielding DNA records with uppercase sequence bytes."""
    acc = ""
    ver = ""
    moltype = ""
    circular = False
    defn_lines: list[str] = []
    org = ""
    tax_lines: list[str] = []
    host = ""
    seq_chunks: list[str] = []
    state = "head"

    def _emit():
        nonlocal acc, ver, moltype, circular, defn_lines, org, tax_lines, host, seq_chunks, state
        if not acc:
            return None
        is_rna = "RNA" in moltype.upper()
        seq_str = "".join(seq_chunks).upper()
        rec = None
        if (allow_rna or not is_rna) and seq_str:
            rec = {
                "accession": ver or acc,
                "base_accession": _base_acc(ver or acc),
                "moltype": moltype,
                "circular": circular,
                "length": len(seq_str),
                "definition": " ".join(defn_lines),
                "organism": org,
                "taxonomy": " ".join(tax_lines),
                "host": host,
                "seq": seq_str.encode("ascii"),
            }
        acc = ""
        ver = ""
        moltype = ""
        circular = False
        defn_lines = []
        org = ""
        tax_lines = []
        host = ""
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
            # Strip line numbers and whitespace
            seq_chunks.append("".join(ch for ch in line if ch.isalpha()))
            continue
        if line.startswith("LOCUS"):
            parts = line.split()
            if len(parts) >= 2:
                acc = parts[1]
            # Find moltype and topology tokens in LOCUS line
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
        elif line.startswith("ORIGIN"):
            state = "origin"
        elif line.startswith("FEATURES"):
            state = "features"
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
            s = line.strip()
            if s.startswith("/host=") and not host:
                host = s[6:].strip().strip('"')


def classify_refseq_viral_subset(rec: dict) -> str | None:
    """Classify a RefSeq viral DNA record into a C-55 SAT1 target subset or None."""
    tax = rec.get("taxonomy") or ""
    tax_l = tax.lower()
    defn_l = (rec.get("definition") or "").lower()
    org_l = (rec.get("organism") or "").lower()
    host_l = (rec.get("host") or "").lower()
    length = int(rec.get("length") or 0)
    circular = bool(rec.get("circular"))

    if any(t.lower() in tax_l for t in ARCHAEAL_VIRUS_TAXA):
        return "archaeal_virus"
    if host_l and any(h in host_l for h in ARCHAEAL_HOST_TOKENS):
        return "archaeal_virus"
    if "archaeal virus" in defn_l or "archaeal virus" in org_l:
        return "archaeal_virus"

    combined = f"{defn_l} {org_l} {tax_l}"
    if any(tok in combined for tok in SATELLITE_TOKENS) or "alphasatellitidae" in tax_l or "tolecusatellitidae" in tax_l:
        return "phage_satellite"

    if "microviridae" in tax_l and length <= 8000:
        return "small_circular_dna"
    if length < 5000 and (
        circular
        or "cressdnaviricota" in tax_l
        or "monodnaviria" in tax_l
        or "circoviridae" in tax_l
        or "anelloviridae" in tax_l
        or "geminiviridae" in tax_l
        or "nanoviridae" in tax_l
    ):
        return "small_circular_dna"
    return None


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


# --------------------------------------------------------------------------------------------------------------
# Stage SAT1 (compact public redistributable=True viral/satellite stage; DEC-65 (2) byte-identical, no 5 kb floor)


def stage_sat1(
    sources_spec: list[dict],
    out: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
    final_dir: Path | None = None,
) -> dict:
    """Stage compact public `SAT1` inputs. Refuses any non-redistributable source."""
    out = Path(out)
    if out.exists():
        raise Refused(f"{out} exists")
    if not sources_spec:
        raise Refused("sources_spec is empty")
    for spec in sources_spec:
        if not spec.get("redistributable", False):
            raise Refused(f"SAT1 refuses non-redistributable source: {spec.get('source_name')}")
        fa = Path(spec["fasta"])
        if not fa.is_file():
            raise Refused(f"missing FASTA file: {fa}")

    work = out.with_name(f"{out.name}.tmp")
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    obj_dir = Path(final_dir) if final_dir is not None else out

    try:
        rep = _stage_sat1(sources_spec, obj_dir, work, root, list(root_stages), owner_inputs, retrieved)
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


def _stage_sat1(
    sources_spec: list[dict],
    obj_dir: Path,
    work: Path,
    root: Path | None,
    root_stages: list[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
) -> dict:
    # Pass 1: copy/normalize each source FASTA byte-identically (no DTR trimming per DEC-65 (2)),
    # drop only empty/invalid-alphabet records, collect all contig_uids for build-order ownership lookup.
    records: list[tuple[str, str, str, bytes, int, bool]] = []  # (src, cat, acc, cu, length, circ)
    objects: list[dict] = []
    skipped_invalid = 0
    by_category: Counter[str] = Counter()
    lt_5kb = 0

    for spec in sources_spec:
        src_name = spec["source_name"]
        cat = spec["category"]
        fa_src = Path(spec["fasta"])
        circ_accs = {_base_acc(a) for a in (spec.get("circular_accessions") or ())} | set(spec.get("circular_accessions") or ())
        out_name = spec.get("out_fasta_name") or f"{cat}_{src_name}.fna"
        fa_dst = work / out_name
        written_recs = 0
        with _open_bin(fa_src) as fh, open(fa_dst, "wb") as w:
            for r in v2.iter_fasta_stream(fh):
                if not r.seq or not VALID_DNA_RE.match(r.seq):
                    skipped_invalid += 1
                    continue
                acc = r.accession
                circ = (acc in circ_accs) or (_base_acc(acc) in circ_accs)
                cu = ids.contig_uid(r.seq)
                w.write(b"\n".join(fasta.format_record(r.header, r.seq, r.line_width)) + b"\n")
                records.append((src_name, cat, acc, cu, len(r.seq), circ))
                by_category[cat] += 1
                if len(r.seq) < 5000:
                    lt_5kb += 1
                written_recs += 1
        if written_recs == 0:
            fa_dst.unlink(missing_ok=True)
            continue
        objects.append({
            "path": str(obj_dir / out_name),
            "source_name": src_name,
            "source_release": spec["source_release"],
            "license": spec["license"],
            "redistributable": True,
            "retrieved_at": retrieved,
        })

    if not records:
        raise Refused("no valid DNA records across SAT1 sources")

    all_uids = [r[3] for r in records]
    owners = v2.owners_of(all_uids, root, root_stages, owner_inputs)

    listing = make_inputs.Listing()
    genomes: list[dict] = []
    rec_rows: list[tuple] = []
    meta_rows: list[tuple] = []
    repeated: list[tuple] = []
    shared: list[tuple] = []
    first_of: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}
    contig_cat: dict[bytes, str] = {}

    for ordinal, (src_name, cat, acc, cu, length, circ) in enumerate(records):
        own = owners.get(cu, "")
        rec_rows.append((ordinal, src_name, cat, acc, cu.hex(), length, str(circ).lower(), own))
        if own:
            shared.append((ordinal, acc, cu.hex(), own))
        contig_len.setdefault(cu, length)
        contig_cat.setdefault(cu, cat)
        g = listing.add([cu], 11, "virus", "complete")
        guid = ids.genome_uid([cu])
        if g is None:
            repeated.append((ordinal, acc, guid.hex(), first_of[guid]))
            continue
        g["circular"] = [cu.hex()] if circ else []
        genomes.append(g)
        first_of[guid] = acc
        meta_rows.append(("SAT1", src_name, cat, guid.hex(), cu.hex(), acc, length, str(circ).lower(), own))

    meta_rows.sort(key=lambda t: (t[2], t[1], t[5], t[4]))
    make_inputs._write_tsv(work / "records.tsv", SAT1_RECORD_HEAD, rec_rows)
    make_inputs._write_tsv(work / "sat1_metadata.tsv", SAT1_META_HEAD, meta_rows)
    make_inputs._write_tsv(work / "repeated_records.tsv", ["ordinal", "accession", "genome_uid", "first_accession"], repeated)
    make_inputs._write_tsv(work / "shared_records.tsv", ["ordinal", "accession", "contig_uid", "owner_stage"], shared)
    (work / "objects.json").write_text(json.dumps(objects, indent=1))
    (work / "genomes.json").write_text(json.dumps(genomes))

    owner_contigs: Counter[str] = Counter()
    owner_bases: Counter[str] = Counter()
    distinct_by_cat: Counter[str] = Counter()
    bases_by_cat: Counter[str] = Counter()
    for cu, n in contig_len.items():
        c = contig_cat[cu]
        distinct_by_cat[c] += 1
        bases_by_cat[c] += n
        if cu in owners:
            st = owners[cu]
            owner_contigs[st] += 1
            owner_bases[st] += n

    held_contigs = sum(owner_contigs.values())
    held_bases = sum(owner_bases.values())

    return {
        "records": len(records),
        "skipped_invalid": skipped_invalid,
        "records_by_category": dict(sorted(by_category.items())),
        "distinct_contigs_by_category": dict(sorted(distinct_by_cat.items())),
        "distinct_bases_by_category": dict(sorted(bases_by_cat.items())),
        "lt_5kb_records": lt_5kb,
        "lt_5kb_contigs": sum(1 for n in contig_len.values() if n < 5000),
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
        "build_order": [*root_stages, *owner_inputs, "SAT1"],
        "objects": len(objects),
        **listing.report(),
    }


# --------------------------------------------------------------------------------------------------------------
# Stage V3 (EMBL VIRE DNA vOTU representatives; DEC-65 (1) metagenomic DTR trimming)


def _float_or(val: str, default: float) -> float:
    if val in ("", "None", "NA", "nan", "NaN"):
        return default
    try:
        return float(val)
    except ValueError:
        return default


def select_vire_votu_reps(
    meta_tsv: Path,
    *,
    exclude_riboviria: bool = True,
) -> tuple[dict[str, dict], list[tuple]]:
    """Select one deterministic best representative per DNA species_cluster from VIRE metadata.

    Ranking matches V2/MV1: highest checkv_completeness (capped at 100), lowest checkv_contamination,
    longest contig_length, then genome_id. Returns (selected_by_genome_id, left_for_later_rows).
    """
    best: dict[str, tuple[tuple, dict]] = {}
    all_rows: list[dict] = []
    with _open_text(Path(meta_tsv)) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gid = r["genome_id"].strip()
            sc = r["species_cluster"].strip()
            realm = (r.get("realm") or "").strip()
            length = int(r["contig_length"])
            comp = _float_or(r.get("checkv_completeness", ""), -1.0)
            cont = _float_or(r.get("checkv_contamination", ""), 100.0)
            info = {
                "genome_id": gid,
                "biosample_name": (r.get("biosample_name") or "").strip(),
                "species_cluster": sc,
                "realm": realm,
                "contig_length": length,
                "checkv_checkv_quality": (r.get("checkv_checkv_quality") or "").strip(),
                "checkv_completeness": (r.get("checkv_completeness") or "").strip(),
                "checkv_contamination": (r.get("checkv_contamination") or "").strip(),
                "genomad_topology": (r.get("genomad_topology") or "").strip(),
                "alignment_to_imgvr": (r.get("alignment_to_imgvr") or "").strip(),
            }
            all_rows.append(info)
            if exclude_riboviria and realm == "Riboviria":
                continue
            rank = (min(comp, 100.0), -cont, length, gid)
            cur = best.get(sc)
            if cur is None or rank > cur[0]:
                best[sc] = (rank, info)

    selected = {info["genome_id"]: info for _, info in best.values()}
    left: list[tuple] = []
    for info in all_rows:
        gid = info["genome_id"]
        if gid in selected:
            continue
        reason = "riboviria_deferred_c50" if (exclude_riboviria and info["realm"] == "Riboviria") else "not_representative"
        left.append((
            gid,
            info["biosample_name"],
            info["species_cluster"],
            info["realm"],
            info["contig_length"],
            info["checkv_checkv_quality"],
            info["checkv_completeness"],
            info["checkv_contamination"],
            info["genomad_topology"],
            info["alignment_to_imgvr"],
            reason,
        ))
    left.sort(key=lambda t: (t[10], t[2], t[0]))
    return selected, left


class _HashingReader(io.RawIOBase):
    """Wrap a binary stream to compute sha256 and total bytes while tarfile streams it."""

    def __init__(self, stream):
        self._s = stream
        self.hasher = hashlib.sha256()
        self.nbytes = 0

    def readable(self) -> bool:
        return True

    def read(self, n: int = -1) -> bytes:
        chunk = self._s.read(n)
        if chunk:
            self.hasher.update(chunk)
            self.nbytes += len(chunk)
        return chunk


def stage_v3_from_tar(
    meta_tsv: Path,
    tar_source: Path | str,
    out_meta: Path,
    out_fasta: Path,
    *,
    root: Path | None,
    root_stages: Iterable[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
    final_fasta: Path | None = None,
    left_for_later_path: Path | None = None,
    exclude_riboviria: bool = True,
) -> dict:
    """Stream VIRE `all_samples_genomes.tar` (local file or `-` for stdin) to stage V3 DNA vOTU reps."""
    out_meta = Path(out_meta)
    out_fasta = Path(out_fasta)
    if out_meta.exists():
        raise Refused(f"{out_meta} exists")
    if out_fasta.exists():
        raise Refused(f"{out_fasta} exists")

    selected, left_rows = select_vire_votu_reps(Path(meta_tsv), exclude_riboviria=exclude_riboviria)
    if not selected:
        raise Refused("no VIRE DNA vOTU representatives selected")

    work_meta = out_meta.with_name(f"{out_meta.name}.tmp")
    work_fa = out_fasta.with_name(f"{out_fasta.name}.tmp")
    shutil.rmtree(work_meta, ignore_errors=True)
    work_meta.mkdir(parents=True)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    if work_fa.exists():
        work_fa.unlink()

    try:
        rep = _stage_v3_stream(
            selected,
            left_rows,
            tar_source,
            work_meta,
            work_fa,
            Path(final_fasta) if final_fasta is not None else out_fasta,
            left_for_later_path or (work_meta / "V3_left_for_later.tsv.gz"),
            root,
            list(root_stages),
            owner_inputs,
            retrieved,
        )
    except BaseException:
        shutil.rmtree(work_meta, ignore_errors=True)
        work_fa.unlink(missing_ok=True)
        raise

    work_fa.rename(out_fasta)
    work_meta.rename(out_meta)
    rep["fasta_path"] = str(out_fasta)
    rep["fasta_bytes"] = out_fasta.stat().st_size
    rep["fasta_sha256"] = ids.sha256_file(out_fasta).hex()
    rep["meta_sha256"] = {
        p.name: ids.sha256_file(p).hex()
        for p in sorted(out_meta.iterdir())
        if p.is_file()
    }
    return rep


def _stage_v3_stream(
    selected: dict[str, dict],
    left_rows: list[tuple],
    tar_source: Path | str,
    work_meta: Path,
    work_fa: Path,
    final_fasta: Path,
    left_path: Path,
    root: Path | None,
    root_stages: list[str],
    owner_inputs: dict[str, Path],
    retrieved: str,
) -> dict:
    # Write left_for_later deterministically (mtime=0)
    left_path.parent.mkdir(parents=True, exist_ok=True)
    with open(left_path, "wb") as raw_l:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_l, compresslevel=1, mtime=0) as gz:
            gz.write(("\t".join(V3_LEFT_HEAD) + "\n").encode("utf-8"))
            for row in left_rows:
                gz.write(("\t".join(str(x) for x in row) + "\n").encode("utf-8"))

    # Group selected genome_ids by biosample_name so we can skip tar members that hold 0 selected reps!
    needed_biosamples = {info["biosample_name"] for info in selected.values() if info["biosample_name"]}

    raw_in = sys.stdin.buffer if str(tar_source) == "-" else open(Path(tar_source), "rb")
    hreader = _HashingReader(raw_in)
    dtr_classes: Counter[str] = Counter()
    trims: list[tuple] = []
    staged_recs: list[tuple] = []  # (gid, sc, realm, cu, staged_len, dtr_cls, circ, src_cu, src_len, k, al)
    seen_gids: set[str] = set()
    members_seen = 0
    members_read = 0
    src_bases = 0
    staged_bases = 0

    try:
        with tarfile.open(fileobj=hreader, mode="r|*") as tf, open(work_fa, "wb") as fa_w:
            for member in tf:
                if not member.isfile():
                    continue
                members_seen += 1
                name = member.name
                # Expected: genomes/vire_biosample_<biosample>.fna.gz
                if needed_biosamples and "vire_biosample_" in name:
                    base = name.rsplit("vire_biosample_", 1)[1]
                    bname = base.split(".", 1)[0]
                    if bname not in needed_biosamples:
                        continue
                fobj = tf.extractfile(member)
                if fobj is None:
                    continue
                members_read += 1
                raw_bytes = fobj.read()
                stream = gzip.GzipFile(fileobj=io.BytesIO(raw_bytes)) if name.endswith(".gz") else io.BytesIO(raw_bytes)
                with stream as fh:
                    for r in v2.iter_fasta_stream(fh):
                        gid = r.accession
                        info = selected.get(gid)
                        if info is None or gid in seen_gids:
                            continue
                        seen_gids.add(gid)
                        src_len = len(r.seq)
                        src_cu = ids.contig_uid(r.seq)
                        src_bases += src_len
                        if info["genomad_topology"] == "DTR":
                            e = v2.ends(r.seq)
                            dtr_cls = v2.classify(e)
                            dtr_classes[dtr_cls] += 1
                            if dtr_cls == "exact":
                                trimmed_k = e.exact
                                seq = v2.trim(r.seq, trimmed_k)
                                circ = True
                            elif dtr_cls == "none":
                                trimmed_k = 0
                                seq = r.seq
                                circ = True
                            else:
                                trimmed_k = 0
                                seq = r.seq
                                circ = False
                        else:
                            dtr_cls = "not_dtr"
                            dtr_classes[dtr_cls] += 1
                            trimmed_k = 0
                            seq = r.seq
                            circ = False

                        fa_w.write(b"\n".join(fasta.format_record(r.header, seq, r.line_width)) + b"\n")
                        staged_len = len(seq)
                        staged_bases += staged_len
                        cu = ids.contig_uid(seq) if trimmed_k else src_cu
                        if trimmed_k:
                            trims.append((
                                "vire", gid, gid, src_cu.hex(), cu.hex(),
                                src_len, trimmed_k, trimmed_k, staged_len,
                            ))
                        staged_recs.append((
                            gid, info["species_cluster"], info["realm"], cu, staged_len,
                            dtr_cls, circ, src_cu, src_len, trimmed_k, info["alignment_to_imgvr"],
                        ))
    finally:
        if str(tar_source) != "-":
            raw_in.close()

    if not staged_recs:
        raise Refused("no VIRE vOTU representative sequences extracted from tar")

    trims.sort(key=lambda t: (t[0], t[1], t[2]))
    all_uids = [r[3] for r in staged_recs]
    owners = v2.owners_of(all_uids, root, root_stages, owner_inputs)

    listing = make_inputs.Listing()
    genomes: list[dict] = []
    rec_rows: list[tuple] = []
    repeated: list[tuple] = []
    shared: list[tuple] = []
    first_of: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}

    for ordinal, (gid, sc, realm, cu, staged_len, dtr_cls, circ, src_cu, src_len, k, al) in enumerate(staged_recs):
        rec_rows.append((
            ordinal, gid, sc, realm, cu.hex(), staged_len, dtr_cls,
            str(circ).lower(), src_cu.hex(), src_len, k, al,
        ))
        if cu in owners:
            shared.append((ordinal, gid, cu.hex(), owners[cu]))
        contig_len.setdefault(cu, staged_len)
        g = listing.add([cu], 11, "virus", "vOTU")
        guid = ids.genome_uid([cu])
        if g is None:
            repeated.append((ordinal, gid, guid.hex(), first_of[guid]))
            continue
        g["circular"] = [cu.hex()] if circ else []
        genomes.append(g)
        first_of[guid] = gid

    make_inputs._write_tsv(work_meta / "records.tsv", V3_RECORD_HEAD, rec_rows)
    make_inputs._write_tsv(work_meta / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
    make_inputs._write_tsv(work_meta / "shared_records.tsv", ["ordinal", "genome", "contig_uid", "owner_stage"], shared)
    make_inputs._write_tsv(work_meta / "trims.tsv", V3_TRIM_HEAD, trims)

    objects = [{
        "path": str(final_fasta),
        "source_name": "vire_genomes",
        "source_release": "VIRE all_samples_genomes.tar (2026-01-29; Fullam et al. NAR 2026, doi:10.1093/nar/gkaf1225)",
        "license": "CC BY-SA 4.0",
        "redistributable": True,
        "retrieved_at": retrieved,
    }]
    (work_meta / "objects.json").write_text(json.dumps(objects, indent=1))
    (work_meta / "genomes.json").write_text(json.dumps(genomes))

    owner_contigs: Counter[str] = Counter()
    owner_bases: Counter[str] = Counter()
    for cu, n in contig_len.items():
        if cu in owners:
            st = owners[cu]
            owner_contigs[st] += 1
            owner_bases[st] += n

    held_contigs = sum(owner_contigs.values())
    held_bases = sum(owner_bases.values())

    return {
        "tar_bytes": hreader.nbytes,
        "tar_sha256": hreader.hasher.hexdigest(),
        "tar_members_seen": members_seen,
        "tar_members_read": members_read,
        "selected_reps": len(selected),
        "extracted_reps": len(staged_recs),
        "missing_reps": len(selected) - len(staged_recs),
        "left_for_later_rows": len(left_rows),
        "left_for_later_path": str(left_path),
        "source_bases": src_bases,
        "staged_bases": staged_bases,
        "dtr_classes": dict(sorted(dtr_classes.items())),
        "trims": len(trims),
        "trimmed_bases": sum(t[7] for t in trims),
        "genomes": len(genomes),
        "contigs": len(contig_len),
        "contig_bases": sum(contig_len.values()),
        "lt_5kb_contigs": sum(1 for n in contig_len.values() if n < 5000),
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
        "build_order": [*root_stages, *owner_inputs, "V3"],
        "objects": len(objects),
        **listing.report(),
    }


# --------------------------------------------------------------------------------------------------------------
# Polite public fetchers (User-Agent: curl, <= 2 req/s, stop on 5 consecutive refusals; DEC-50 / AGENTS §8)


def _http_get(url: str, *, min_interval: float = 0.55, state: dict | None = None, timeout: int = 60) -> bytes:
    if state is None:
        state = {"last": 0.0, "refusals": 0}
    wait = min_interval - (time.monotonic() - state["last"])
    if wait > 0:
        time.sleep(wait)
    req = urllib.request.Request(url, headers={"User-Agent": "curl"})
    for attempt in range(3):
        state["last"] = time.monotonic()
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            state["refusals"] = 0
            return data
        except urllib.error.HTTPError as e:
            if e.code in (403, 429, 500, 502, 503, 504):
                state["refusals"] += 1
                if state["refusals"] >= 5:
                    raise Refused(f"stopped after 5 consecutive HTTP refusals (last {e.code} on {url})") from e
                time.sleep(2.0 * (attempt + 1))
                continue
            raise
        except (urllib.error.URLError, TimeoutError) as e:
            state["refusals"] += 1
            if state["refusals"] >= 5:
                raise Refused(f"stopped after 5 consecutive network errors (last {e} on {url})") from e
            time.sleep(2.0 * (attempt + 1))
    raise Refused(f"failed after 3 attempts: {url}")


def _download_file(url: str, dst: Path, *, state: dict) -> tuple[int, str]:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.is_file() and dst.stat().st_size > 0:
        return dst.stat().st_size, ids.sha256_file(dst).hex()
    tmp = dst.with_name(f".{dst.name}.tmp")
    req = urllib.request.Request(url, headers={"User-Agent": "curl"})
    wait = 0.55 - (time.monotonic() - state["last"])
    if wait > 0:
        time.sleep(wait)
    state["last"] = time.monotonic()
    h = hashlib.sha256()
    n = 0
    with urllib.request.urlopen(req, timeout=300) as resp, open(tmp, "wb") as fh:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
            h.update(chunk)
            n += len(chunk)
    tmp.rename(dst)
    state["refusals"] = 0
    return n, h.hexdigest()


def _ncbi_esearch_all(term: str, *, state: dict, batch: int = 10000) -> list[str]:
    uids: list[str] = []
    retstart = 0
    while True:
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore"
            f"&retmax={batch}&retstart={retstart}&retmode=json&term={urllib.parse.quote(term)}"
        )
        d = json.loads(_http_get(url, state=state))["esearchresult"]
        got = d.get("idlist") or []
        uids.extend(got)
        total = int(d.get("count") or 0)
        retstart += len(got)
        if not got or retstart >= total:
            break
    return uids


def _fetch_ncbi_query_fasta(
    name: str,
    term: str,
    out_dir: Path,
    *,
    state: dict,
    exclude_base_accs: set[str] | None = None,
    batch_size: int = 200,
) -> dict:
    """Fetch DNA-only records for an NCBI nuccore query via esearch + esummary + efetch."""
    out_dir.mkdir(parents=True, exist_ok=True)
    fa_path = out_dir / f"{name}.fna"
    meta_path = out_dir / f"{name}_meta.tsv"
    if fa_path.is_file() and meta_path.is_file():
        circ = set()
        n_rows = 0
        with open(meta_path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                n_rows += 1
                if r["circular"] == "true":
                    circ.add(r["accession"])
        return {
            "name": name, "fasta": str(fa_path), "meta": str(meta_path),
            "bytes": fa_path.stat().st_size, "sha256": ids.sha256_file(fa_path).hex(),
            "records": n_rows, "circular_count": len(circ), "cached": True,
        }

    uids = _ncbi_esearch_all(term, state=state)
    exclude = exclude_base_accs or set()
    kept_uids: list[str] = []
    meta_rows: list[tuple[str, str, int, str, str]] = []  # (uid, acc, slen, circ, title)

    for i in range(0, len(uids), batch_size):
        chunk = uids[i:i + batch_size]
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore"
            f"&retmode=json&id={','.join(chunk)}"
        )
        res = json.loads(_http_get(url, state=state)).get("result") or {}
        for uid in chunk:
            info = res.get(uid)
            if not isinstance(info, dict):
                continue
            mol = (info.get("moltype") or "").lower()
            if mol and "dna" not in mol:
                continue
            acc = info.get("accessionversion") or info.get("caption") or ""
            if not acc or _base_acc(acc) in exclude:
                continue
            circ = "true" if (info.get("topology") or "").lower() == "circular" else "false"
            slen = int(info.get("slen") or 0)
            title = (info.get("title") or "").replace("\t", " ").replace("\n", " ")
            kept_uids.append(uid)
            meta_rows.append((uid, acc, slen, circ, title))

    tmp_fa = fa_path.with_name(f".{fa_path.name}.tmp")
    with open(tmp_fa, "wb") as w:
        for i in range(0, len(kept_uids), batch_size):
            chunk = kept_uids[i:i + batch_size]
            url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore"
                f"&rettype=fasta&retmode=text&id={','.join(chunk)}"
            )
            data = _http_get(url, state=state, timeout=120)
            w.write(data)
            if data and not data.endswith(b"\n"):
                w.write(b"\n")
    tmp_fa.rename(fa_path)
    make_inputs._write_tsv(meta_path, ["uid", "accession", "slen", "circular", "title"], meta_rows)
    return {
        "name": name,
        "term": term,
        "esearch_uids": len(uids),
        "kept_uids": len(kept_uids),
        "fasta": str(fa_path),
        "meta": str(meta_path),
        "bytes": fa_path.stat().st_size,
        "sha256": ids.sha256_file(fa_path).hex(),
        "circular_count": sum(1 for r in meta_rows if r[3] == "true"),
        "cached": False,
    }


def _cmd_fetch_sat1(a) -> int:
    dl = Path(a.dl_dir)
    dl.mkdir(parents=True, exist_ok=True)
    state = {"last": 0.0, "refusals": 0}
    report: dict = {"sources": {}}

    # 1. INPHARED 9Aug2025
    inph_dir = dl / "inphared_9aug2025"
    b1, h1 = _download_file(
        "https://millardlab-inphared.s3.climb.ac.uk/9Aug2025_data.tsv.gz",
        inph_dir / "9Aug2025_data.tsv.gz", state=state,
    )
    b2, h2 = _download_file(
        "https://millardlab-inphared.s3.climb.ac.uk/9Aug2025_genomes.fa.gz",
        inph_dir / "9Aug2025_genomes.fa.gz", state=state,
    )
    # Extract S1 14Apr2025_data.tsv.gz from S1 tar if tar given
    s1_tsv = inph_dir / "14Apr2025_data.tsv.gz"
    if not s1_tsv.is_file():
        with tarfile.open(a.s1_tar, "r") as tf:
            fobj = tf.extractfile("14Apr2025_data.tsv.gz")
            if fobj is None:
                raise Refused(f"14Apr2025_data.tsv.gz not found in {a.s1_tar}")
            s1_tsv.write_bytes(fobj.read())
    s1_accs = _read_inphared_accessions(s1_tsv)
    inph_9aug_accs = _read_inphared_accessions(inph_dir / "9Aug2025_data.tsv.gz")
    delta_accs = inph_9aug_accs - s1_accs

    # Extract delta FASTA from 9Aug2025_genomes.fa.gz
    inph_delta_fa = inph_dir / "inphared_9aug2025_delta.fna"
    n_delta = 0
    with gzip.open(inph_dir / "9Aug2025_genomes.fa.gz", "rb") as fh, open(inph_delta_fa, "wb") as w:
        for r in v2.iter_fasta_stream(fh):
            if _base_acc(r.accession) in delta_accs:
                w.write(b"\n".join(fasta.format_record(r.header, r.seq, r.line_width)) + b"\n")
                n_delta += 1
    report["sources"]["inphared_9aug2025"] = {
        "data_tsv_bytes": b1, "data_tsv_sha256": h1,
        "genomes_fa_bytes": b2, "genomes_fa_sha256": h2,
        "s1_accs": len(s1_accs), "new_accs": len(inph_9aug_accs),
        "delta_accs": len(delta_accs), "delta_fasta_records": n_delta,
        "delta_fasta": str(inph_delta_fa),
    }

    # 2. NCBI RefSeq Viral release 238 (viral.1.genomic.gbff.gz)
    rs_dir = dl / "refseq_viral_r238"
    b3, h3 = _download_file(
        "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz",
        rs_dir / "viral.1.genomic.gbff.gz", state=state,
    )
    rs_out = {
        "small_circular_dna": open(rs_dir / "refseq_small_circular_dna.fna", "wb"),
        "archaeal_virus": open(rs_dir / "refseq_archaeal_virus.fna", "wb"),
        "phage_satellite": open(rs_dir / "refseq_phage_satellite.fna", "wb"),
    }
    rs_circ: dict[str, list[str]] = {k: [] for k in rs_out}
    rs_counts: Counter[str] = Counter()
    rs_taken_accs: set[str] = set()
    try:
        with gzip.open(rs_dir / "viral.1.genomic.gbff.gz", "rt", encoding="utf-8", errors="replace") as fh:
            for rec in parse_genbank_dna_records(fh, allow_rna=False):
                cat = classify_refseq_viral_subset(rec)
                if cat is None:
                    continue
                rs_counts[cat] += 1
                rs_taken_accs.add(rec["base_accession"])
                if rec["circular"]:
                    rs_circ[cat].append(rec["accession"])
                hd = f"{rec['accession']} {rec['definition']}"
                rs_out[cat].write(b"\n".join(fasta.format_record(hd, rec["seq"], 70)) + b"\n")
    finally:
        for w in rs_out.values():
            w.close()
    (rs_dir / "refseq_circular.json").write_text(json.dumps(rs_circ, indent=1, sort_keys=True))
    report["sources"]["ncbi_refseq_viral_r238"] = {
        "gbff_bytes": b3, "gbff_sha256": h3,
        "counts": dict(rs_counts),
        "circular_counts": {k: len(v) for k, v in rs_circ.items()},
    }

    # 3. NCBI E-utilities queries (excluding accessions already taken from S1, INPHARED 9Aug2025, or RefSeq r238)
    nc_dir = dl / "ncbi_eutils"
    already_phage = s1_accs | inph_9aug_accs
    queries = [
        (
            "genbank_phage_delta",
            "(phage[Title] OR bacteriophage[Title]) AND complete genome[Title] AND 2025/04/14:2026/12/31[PDAT] AND 1000:1000000[SLEN]",
            already_phage,
        ),
        (
            "genbank_phage_satellites",
            '("phage satellite"[All Fields] OR "phage-inducible chromosomal island"[All Fields] OR "pathogenicity island SaPI"[All Fields] OR "PICI-like element"[All Fields] OR "satellite phage"[All Fields] OR "bacteriophage P4"[Organism] OR "Vibrio phage PLE"[All Fields] OR "capsid-forming PICI"[All Fields] OR "Enterobacteria phage P4"[Organism] OR "Staphylococcus pathogenicity island"[Title]) AND 3000:60000[SLEN]',
            rs_taken_accs,
        ),
        (
            "genbank_archaeal_viruses",
            "(Fuselloviridae[Organism] OR Lipothrixviridae[Organism] OR Rudiviridae[Organism] OR Bicaudaviridae[Organism] OR Ampullaviridae[Organism] OR Globuloviridae[Organism] OR Guttaviridae[Organism] OR Pleolipoviridae[Organism] OR Sphaerolipoviridae[Organism] OR Turriviridae[Organism] OR Tristromaviridae[Organism] OR Portogloboviridae[Organism] OR Ovaliviridae[Organism] OR Halspiviridae[Organism] OR Thaspiviridae[Organism] OR Spiraviridae[Organism] OR Clavaviridae[Organism] OR Adnaviria[Organism] OR archaeal virus[All Fields]) AND complete[Title] AND 1000:300000[SLEN]",
            rs_taken_accs,
        ),
        (
            "genbank_tisza_prjna573942",
            "PRJNA573942[BioProject]",
            rs_taken_accs,
        ),
        (
            "genbank_microviridae",
            "txid10841[Organism:exp] AND complete genome[Title] AND 3000:8000[SLEN]",
            rs_taken_accs,
        ),
        (
            "genbank_cressdnaviricota_lt5k",
            "Cressdnaviricota[Organism] AND complete genome[Title] AND 500:5000[SLEN]",
            rs_taken_accs,
        ),
    ]
    for qname, qterm, excl in queries:
        report["sources"][qname] = _fetch_ncbi_query_fasta(qname, qterm, nc_dir, state=state, exclude_base_accs=excl)

    out_rep = dl / "fetch_sat1_report.json"
    out_rep.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True))
    return 0


def _load_circ_from_meta(meta_tsv: Path) -> set[str]:
    circ: set[str] = set()
    if not meta_tsv.is_file():
        return circ
    with open(meta_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("circular") == "true":
                circ.add(r["accession"])
    return circ


def _cmd_stage_sat1(a) -> int:
    dl = Path(a.dl_dir)
    rs_dir = dl / "refseq_viral_r238"
    nc_dir = dl / "ncbi_eutils"
    rs_circ = json.loads((rs_dir / "refseq_circular.json").read_text())

    specs = [
        {
            "source_name": "inphared_9aug2025",
            "source_release": "INPHARED 9Aug2025 delta vs 14Apr2025",
            "license": "INSDC (public)",
            "redistributable": True,
            "category": "inphared_delta",
            "fasta": dl / "inphared_9aug2025" / "inphared_9aug2025_delta.fna",
            "out_fasta_name": "inphared_9aug2025_delta.fna",
            "circular_accessions": set(),
        },
        {
            "source_name": "ncbi_eutils_c55_delta",
            "source_release": f"NCBI GenBank E-utilities ({a.retrieved[:10]})",
            "license": "INSDC (public)",
            "redistributable": True,
            "category": "genbank_phage_delta",
            "fasta": nc_dir / "genbank_phage_delta.fna",
            "out_fasta_name": "genbank_phage_delta.fna",
            "circular_accessions": _load_circ_from_meta(nc_dir / "genbank_phage_delta_meta.tsv"),
        },
        {
            "source_name": "ncbi_refseq_viral_r238",
            "source_release": "NCBI RefSeq Viral release 238 (2026-11-13)",
            "license": "INSDC / NCBI RefSeq (public domain)",
            "redistributable": True,
            "category": "phage_satellite",
            "fasta": rs_dir / "refseq_phage_satellite.fna",
            "out_fasta_name": "refseq_phage_satellite.fna",
            "circular_accessions": set(rs_circ.get("phage_satellite", [])),
        },
        {
            "source_name": "ncbi_eutils_c55_delta",
            "source_release": f"NCBI GenBank E-utilities ({a.retrieved[:10]})",
            "license": "INSDC (public)",
            "redistributable": True,
            "category": "phage_satellite",
            "fasta": nc_dir / "genbank_phage_satellites.fna",
            "out_fasta_name": "genbank_phage_satellites.fna",
            "circular_accessions": _load_circ_from_meta(nc_dir / "genbank_phage_satellites_meta.tsv"),
        },
        {
            "source_name": "ncbi_refseq_viral_r238",
            "source_release": "NCBI RefSeq Viral release 238 (2026-11-13)",
            "license": "INSDC / NCBI RefSeq (public domain)",
            "redistributable": True,
            "category": "archaeal_virus",
            "fasta": rs_dir / "refseq_archaeal_virus.fna",
            "out_fasta_name": "refseq_archaeal_virus.fna",
            "circular_accessions": set(rs_circ.get("archaeal_virus", [])),
        },
        {
            "source_name": "ncbi_eutils_c55_delta",
            "source_release": f"NCBI GenBank E-utilities ({a.retrieved[:10]})",
            "license": "INSDC (public)",
            "redistributable": True,
            "category": "archaeal_virus",
            "fasta": nc_dir / "genbank_archaeal_viruses.fna",
            "out_fasta_name": "genbank_archaeal_viruses.fna",
            "circular_accessions": _load_circ_from_meta(nc_dir / "genbank_archaeal_viruses_meta.tsv"),
        },
        {
            "source_name": "ncbi_refseq_viral_r238",
            "source_release": "NCBI RefSeq Viral release 238 (2026-11-13)",
            "license": "INSDC / NCBI RefSeq (public domain)",
            "redistributable": True,
            "category": "small_circular_dna",
            "fasta": rs_dir / "refseq_small_circular_dna.fna",
            "out_fasta_name": "refseq_small_circular_dna.fna",
            "circular_accessions": set(rs_circ.get("small_circular_dna", [])),
        },
        {
            "source_name": "ncbi_eutils_c55_delta",
            "source_release": f"NCBI GenBank E-utilities ({a.retrieved[:10]})",
            "license": "INSDC (public)",
            "redistributable": True,
            "category": "small_circular_dna",
            "fasta": nc_dir / "genbank_tisza_prjna573942.fna",
            "out_fasta_name": "genbank_tisza_prjna573942.fna",
            "circular_accessions": _load_circ_from_meta(nc_dir / "genbank_tisza_prjna573942_meta.tsv"),
        },
        {
            "source_name": "ncbi_eutils_c55_delta",
            "source_release": f"NCBI GenBank E-utilities ({a.retrieved[:10]})",
            "license": "INSDC (public)",
            "redistributable": True,
            "category": "small_circular_dna",
            "fasta": nc_dir / "genbank_microviridae.fna",
            "out_fasta_name": "genbank_microviridae.fna",
            "circular_accessions": _load_circ_from_meta(nc_dir / "genbank_microviridae_meta.tsv"),
        },
        {
            "source_name": "ncbi_eutils_c55_delta",
            "source_release": f"NCBI GenBank E-utilities ({a.retrieved[:10]})",
            "license": "INSDC (public)",
            "redistributable": True,
            "category": "small_circular_dna",
            "fasta": nc_dir / "genbank_cressdnaviricota_lt5k.fna",
            "out_fasta_name": "genbank_cressdnaviricota_lt5k.fna",
            "circular_accessions": _load_circ_from_meta(nc_dir / "genbank_cressdnaviricota_lt5k_meta.tsv"),
        },
    ]

    root_stages = [s for s in a.root_stages.split(",") if s]
    owners = _parse_owner_inputs(a.owner_inputs)
    rep = stage_sat1(
        specs,
        a.out_dir,
        root=a.root,
        root_stages=root_stages,
        owner_inputs=owners,
        retrieved=a.retrieved,
        final_dir=a.final_dir,
    )
    if a.out_report is not None:
        a.out_report.parent.mkdir(parents=True, exist_ok=True)
        a.out_report.write_text(json.dumps(rep, indent=1, sort_keys=True) + "\n")
    print(json.dumps(rep, sort_keys=True))
    return 0


def _cmd_stage_v3(a) -> int:
    root_stages = [s for s in a.root_stages.split(",") if s]
    owners = _parse_owner_inputs(a.owner_inputs)
    rep = stage_v3_from_tar(
        a.meta_tsv,
        a.tar_source,
        a.out_meta,
        a.out_fasta,
        root=a.root,
        root_stages=root_stages,
        owner_inputs=owners,
        retrieved=a.retrieved,
        final_fasta=a.final_fasta,
        left_for_later_path=a.left_for_later,
    )
    if a.out_report is not None:
        a.out_report.parent.mkdir(parents=True, exist_ok=True)
        a.out_report.write_text(json.dumps(rep, indent=1, sort_keys=True) + "\n")
    print(json.dumps(rep, sort_keys=True))
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Card C-55: fetch and stage SAT1 and V3 (VIRE).")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_f = sub.add_parser("fetch-sat1")
    p_f.add_argument("--dl-dir", required=True, type=Path)
    p_f.add_argument("--s1-tar", required=True, type=Path)

    p_s = sub.add_parser("stage-sat1")
    p_s.add_argument("--dl-dir", required=True, type=Path)
    p_s.add_argument("--out-dir", required=True, type=Path)
    p_s.add_argument("--retrieved", required=True)
    p_s.add_argument("--root", type=Path, default=None)
    p_s.add_argument("--root-stages", default="")
    p_s.add_argument("--owner-inputs", action="append", default=[])
    p_s.add_argument("--final-dir", type=Path, default=None)
    p_s.add_argument("--out-report", type=Path, default=None)

    p_v = sub.add_parser("stage-v3")
    p_v.add_argument("--meta-tsv", required=True, type=Path)
    p_v.add_argument("--tar-source", required=True)
    p_v.add_argument("--out-meta", required=True, type=Path)
    p_v.add_argument("--out-fasta", required=True, type=Path)
    p_v.add_argument("--retrieved", required=True)
    p_v.add_argument("--root", type=Path, default=None)
    p_v.add_argument("--root-stages", default="")
    p_v.add_argument("--owner-inputs", action="append", default=[])
    p_v.add_argument("--final-fasta", type=Path, default=None)
    p_v.add_argument("--left-for-later", type=Path, default=None)
    p_v.add_argument("--out-report", type=Path, default=None)

    a = ap.parse_args(argv)
    try:
        if a.cmd == "fetch-sat1":
            return _cmd_fetch_sat1(a)
        if a.cmd == "stage-sat1":
            return _cmd_stage_sat1(a)
        if a.cmd == "stage-v3":
            return _cmd_stage_v3(a)
    except Refused as e:
        print(f"c55_vire_small_viruses.py: refused: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
