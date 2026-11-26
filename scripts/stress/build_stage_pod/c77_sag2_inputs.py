"""Card C-77 (DEC-79): stage the GORG-Dark single cells (`SAG2`) from the figshare zip `GORG-Dark_SAGs.zip`.

    python c77_sag2_inputs.py census ZIP
    python c77_sag2_inputs.py stage ZIP ENA_TSV INDEX_DIR META_OUT FASTA_OUT [--root ROOT]
        [--root-stages S1,P1,...] [--chain-home /tmp/if2453/s2/inputs] [--later NAME=DIR ...]
        [--inputs-home /tmp/if2453/s2/inputs] [--retrieved-at T] [--work DIR] [--threads N] [--workers N]
        [--memory-limit 4GB] [--report-out F]

The zip holds one gzip FASTA per single cell (`GORG-Dark SAGs/SCGC_<id>_contigs.fasta.gz`) and the macOS resource
forks (`__MACOSX/...`), which are not data. A member's name carries no accession and no assembly name, so the only
link to card C-58's ENA table (`ena_PRJNA1259050_assemblies.tsv`, ENA's `base_count` per assembly) is the total of
bases: a cell is matched to an assembly when exactly one cell and exactly one row have that total. That is an
inference, labelled so (a tie matches nothing; two cells never take one row).

The accession test (card C-77, DEC-79 (3)) is C-74's `r232_index.lookup`: a matched cell whose accession is in GTDB R232
(in full) or in GlobDB r232 (only the resources AMXMAG, GEM, GOMC and SPIRE carry INSDC accessions there) is STAGED
all the same (the orchestrating session, 2026-11-26: no r232 genome is held here, and the r232 cut of r26.11.6 dedups
against this release); its status is written for every genome in `r232_status.tsv`: in_gtdb_r232 (with representative
and species), in_globdb_r232, in_gtdb_and_globdb_r232, not_in_r232, or not_visible (no ENA row to take an accession
from). There is no ANI test and nothing is removed but what DEC-27 and `s2_check_inputs.py` require.

Staging (the method of card C-58's `stage_sag1`): the cells in id order; the genetic code screened per cell
(`c58.screen_genetic_code`: 11, 25 or 4, the ENA organism name as the taxonomic hint); a contig name that already
carries the cell's id (`AM-100-A01_NODE_1` for `SCGC_AM-100-A01`) is kept, any other is given the member's stem and
`|`; `kind="prokaryote"`, `assembly_quality="contig"`, whole contigs; DEC-27 repeated contig sets listed once
(`repeated_records.tsv`); `shared_records.tsv` the contigs the root's stored stages or the `r26.11.4` stages hold (the
template's `EXPECT=held_by` lines), `later_records.tsv` those that only another stage of the `r26.11.5` chain before
this one holds (`--later`). One BGZF FASTA (written here, no tool), its metadata under META_OUT/SAG2/.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import sys
import zipfile
import zlib
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import a26_r26_11_4_chain as a26  # noqa: E402
import c58_sags_symbionts_niche as c58  # noqa: E402
import make_inputs  # noqa: E402
import r232_index  # noqa: E402
from biobox import fasta, ids  # noqa: E402
from biobox.l0seq import bgzf  # noqa: E402

Refused = make_inputs.Refused

STAGE = "SAG2"
FASTA_NAME = "sag2_gorg_dark.fna.gz"
MEMBER_SUFFIX = "_contigs.fasta.gz"
ROOT_STAGES = ("S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1")
SOURCE = {
    "source_name": "gorg_dark_sags",
    "source_release": "figshare file 54602522, GORG-Dark_SAGs.zip (members dated 2025-05-15)",
    "license": "CC BY 4.0 as the job list states it (sources.yaml tower_sag2); no provider record was read",
    "redistributable": False,  # DEC-80 (3): true only where a recorded licence permits it, none is confirmed
}
ENA_NEED = ("accession", "assembly_name", "base_count", "scientific_name")


# --- the zip -------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Member:
    name: str
    sag: str
    size: int


def _kind(info: zipfile.ZipInfo) -> str:
    if info.is_dir():
        return "directory"
    if info.filename.startswith("__MACOSX/"):
        return "macosx"
    if info.filename.endswith(MEMBER_SUFFIX):
        return "sag"
    return "other"


def census(zip_path: Path) -> dict:
    """What the zip holds, entry by entry: nothing is passed over unnamed."""
    with zipfile.ZipFile(zip_path) as zf:
        infos = zf.infolist()
    kinds = Counter(_kind(i) for i in infos)
    return {"entries": len(infos), "directories": kinds["directory"], "macosx": kinds["macosx"],
            "sags": kinds["sag"], "other": sorted(i.filename for i in infos if _kind(i) == "other")}


def members(zip_path: Path) -> list[Member]:
    """The single-cell members in id order (the id is the member's stem)."""
    with zipfile.ZipFile(zip_path) as zf:
        out = [Member(i.filename, i.filename.rsplit("/", 1)[-1][: -len(MEMBER_SUFFIX)], i.file_size)
               for i in zf.infolist() if _kind(i) == "sag"]
    out.sort(key=lambda m: m.sag)
    dup = [a.sag for a, b in zip(out, out[1:]) if a.sag == b.sag]
    if dup:
        raise Refused(f"two members for single cell {dup[0]}")
    return out


def _records(zf: zipfile.ZipFile, m: Member) -> list[tuple[str, bytes, int]]:
    """(header, sequence, line width) of a member's records, empty sequences dropped."""
    out: list[tuple[str, bytes, int]] = []
    header: str | None = None
    lines: list[bytes] = []
    try:
        with zf.open(m.name) as raw, gzip.open(raw, "rb") as fh:
            for ln in fh:
                ln = ln.rstrip(b"\r\n")
                if not ln:
                    continue
                if ln.startswith(b">"):
                    if header is not None and lines:
                        out.append((header, b"".join(lines), fasta._width(lines)))
                    header, lines = ln[1:].decode("utf-8", "replace"), []
                elif header is None:
                    raise Refused(f"{m.sag}: sequence before the first header")
                else:
                    lines.append(ln)
        if header is not None and lines:
            out.append((header, b"".join(lines), fasta._width(lines)))
    except (OSError, EOFError, zlib.error) as exc:
        raise Refused(f"{m.sag}: member {m.name} is not a readable gzip FASTA: {exc}") from exc
    return out


def scan(zip_path: Path) -> dict[str, dict]:
    """Per cell: total bases, contigs, the contig uids (16 bytes each, in record order) and lengths."""
    out: dict[str, dict] = {}
    with zipfile.ZipFile(zip_path) as zf:
        for m in members(zip_path):
            recs = _records(zf, m)
            out[m.sag] = {"member": m.name, "bp": sum(len(s) for _, s, _ in recs), "contigs": len(recs),
                          "uids": b"".join(ids.contig_uid(s) for _, s, _ in recs),
                          "lengths": [len(s) for _, s, _ in recs]}
    return out


# --- the ENA table -------------------------------------------------------------------------------------------------

def read_ena(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        miss = [c for c in ENA_NEED if c not in (rd.fieldnames or [])]
        if miss:
            raise Refused(f"ENA table {path}: no column {miss[0]}")
        return list(rd)


def match(scanned: dict[str, dict], ena: list[dict]) -> dict[str, dict]:
    """Cell to ENA row by total bases: exactly one cell and exactly one row with that total."""
    cells: dict[int, list[str]] = defaultdict(list)
    for s in sorted(scanned):
        cells[scanned[s]["bp"]].append(s)
    rows: dict[int, list[dict]] = defaultdict(list)
    for r in ena:
        rows[int(r["base_count"])].append(r)
    out = {}
    for s in sorted(scanned):
        bp = scanned[s]["bp"]
        cand = rows.get(bp, [])
        if not cand:
            out[s] = {"accession": None, "basis": "no_ena_row", "row": None}
        elif len(cand) == 1 and len(cells[bp]) == 1:
            out[s] = {"accession": cand[0]["accession"], "basis": "base_count", "row": cand[0]}
        else:
            out[s] = {"accession": None, "basis": "ambiguous_base_count", "row": None}
    return out


# --- the FASTA writer ----------------------------------------------------------------------------------------------

class BgzfWriter:
    """Deterministic BGZF (gzip members of at most 65,280 bytes with the BC field, EOF block last); no tool."""

    def __init__(self, path: Path, level: int = 6):
        self.fh = open(path, "wb", buffering=1 << 22)
        self.buf = bytearray()
        self.level = level

    def write(self, data: bytes) -> None:
        self.buf += data
        while len(self.buf) >= bgzf.BLOCK_DATA:
            self.fh.write(bgzf._block(bytes(self.buf[: bgzf.BLOCK_DATA]), self.level))
            del self.buf[: bgzf.BLOCK_DATA]

    def close(self) -> None:
        if self.buf:
            self.fh.write(bgzf._block(bytes(self.buf), self.level))
            self.buf.clear()
        self.fh.write(bgzf.EOF_BLOCK)
        self.fh.close()


# --- staging -------------------------------------------------------------------------------------------------------

META_FILES = ("genomes.json", "objects.json", "objects_staging.json", "records.tsv", "repeated_records.tsv",
              "shared_records.tsv", "later_records.tsv", "genetic_codes.tsv", "sag_match.tsv", "r232_status.tsv")
MATCH_HEAD = ["sag", "member", "contigs", "bp", "ena_accession", "ena_assembly_name", "ena_base_count",
              "ena_scientific_name", "basis", "seen_by_accession_test", "in_gtdb_r232", "in_globdb_r232", "status",
              "reason"]
STATUS_HEAD = ["genome", "ena_accession", "basis", "r232_status", "in_gtdb_r232", "is_gtdb_representative",
               "gtdb_representative", "gtdb_species", "in_globdb_r232", "globdb_r232_id", "contigs", "bp"]


def r232_status(x: dict | None) -> str:
    """What C-74's `lookup` says of a cell's ENA accession; a cell with none is not visible to the test."""
    if x is None:
        return "not_visible"
    if x["in_gtdb_r232"] and x["in_globdb_r232"]:
        return "in_gtdb_and_globdb_r232"
    if x["in_gtdb_r232"]:
        return "in_gtdb_r232"
    return "in_globdb_r232" if x["in_globdb_r232"] else "not_in_r232"


def _yn(v: bool) -> str:
    return "yes" if v else "no"


def stage(zip_path: Path, ena_tsv: Path, index_dir: Path, meta_out: Path, fasta_out: Path, *, root: Path | None,
          root_stages, owner_inputs: dict[str, Path], later_inputs: dict[str, Path], inputs_home: Path,
          retrieved_at: str = "2026-11-25T23:15:40Z", work: Path | None = None, threads: int = 2,
          memory_limit: str = "4GB", workers: int = 1) -> dict:
    """Stage SAG2; returns the report (module docstring). Nothing is left behind when it refuses."""
    zip_path, ena_tsv, index_dir = Path(zip_path), Path(ena_tsv), Path(index_dir)
    meta_dir, fa_dir = Path(meta_out) / STAGE, Path(fasta_out) / STAGE
    final_fa = fa_dir / FASTA_NAME
    for p in (meta_dir, final_fa):
        if p.exists():
            raise Refused(f"{p} exists")
    if not ena_tsv.is_file():
        raise Refused(f"ENA table {ena_tsv} is missing")
    if not (index_dir / "gtdb_r232.parquet").is_file() or not (index_dir / "globdb_r232.parquet").is_file():
        raise Refused(f"{index_dir} is not card C-74's index (no gtdb_r232.parquet or globdb_r232.parquet)")
    if not members(zip_path):
        raise Refused(f"{zip_path}: no single cell member ({MEMBER_SUFFIX})")
    scanned = scan(zip_path)
    ena = read_ena(ena_tsv)
    matched = match(scanned, ena)

    accs = [r["accession"] for r in matched.values() if r["accession"]]
    look = {x["accession"]: x for x in r232_index.lookup(index_dir, accs)}
    removed: dict[str, dict] = {}  # none: r232 status is recorded, never a reason to drop (coordinator, 2026-11-26)
    kept = sorted(scanned)

    tmp = Path(meta_out) / f"{STAGE}.tmp"
    partial = fa_dir / f".{FASTA_NAME}.partial"
    work_dir = Path(work) if work is not None else Path(meta_out) / f"{STAGE}.work"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    fa_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    try:
        rep = _stage(zip_path, scanned, matched, look, removed, kept, tmp, partial, work_dir, root,
                     list(root_stages), owner_inputs, later_inputs, Path(inputs_home), final_fa, retrieved_at,
                     threads, memory_limit, workers)
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        partial.unlink(missing_ok=True)
        if not any(fa_dir.iterdir()):
            fa_dir.rmdir()
        raise
    if work is None:
        shutil.rmtree(work_dir, ignore_errors=True)
    partial.rename(final_fa)
    tmp.rename(meta_dir)
    files = {n: meta_dir / n for n in META_FILES}
    files[FASTA_NAME] = final_fa
    rep["sha256"] = {n: ids.sha256_file(p).hex() for n, p in files.items()}
    rep["bytes"] = {n: p.stat().st_size for n, p in files.items()}
    return rep


def _stage(zip_path, scanned, matched, look, removed, kept, tmp, partial, work_dir, root, root_stages, owner_inputs,
           later_inputs, inputs_home, final_fa, retrieved_at, threads, memory_limit, workers) -> dict:
    members_by_sag = {m.sag: m for m in members(zip_path)}
    gc_results: dict[str, dict] = {}
    pool = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None
    pending: list[tuple[str, object]] = []

    def submit(sag: str, seqs: list[bytes], org: str) -> None:
        if pool is None:
            gc_results[sag] = c58.screen_genetic_code(seqs, org, screen_max_bp=500_000)
            return
        pending.append((sag, pool.submit(c58._screen_worker, (sag, seqs, org, 500_000))))
        while len(pending) >= 4 * workers:
            s, fut = pending.pop(0)
            gc_results[s] = fut.result()[1]

    names: set[str] = set()
    n_renamed = 0
    w = BgzfWriter(partial)
    try:
        with zipfile.ZipFile(zip_path) as zf:
            for sag in kept:
                recs = _records(zf, members_by_sag[sag])
                cell_id = sag[len("SCGC_"):] if sag.startswith("SCGC_") else sag
                row = matched[sag]["row"]
                for hdr, seq, width in recs:
                    first = hdr.split(None, 1)[0] if hdr.strip() else ""
                    if not first.startswith(cell_id + "_"):
                        hdr = f"{sag}|{hdr}"
                        n_renamed += 1
                    first = hdr.split(None, 1)[0]
                    if first in names:
                        raise Refused(f"contig name {first} occurs twice in the stage")
                    names.add(first)
                    w.write(b"\n".join(fasta.format_record(hdr, seq, width)) + b"\n")
                submit(sag, [s for _, s, _ in recs], row["scientific_name"] if row else "")
        while pending:
            s, fut = pending.pop(0)
            gc_results[s] = fut.result()[1]
    finally:
        w.close()
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)

    contigs = {s: (scanned[s]["uids"], scanned[s]["lengths"]) for s in kept}
    resolved, genome_uids, conflicts = c58._resolve_within_stage_code_conflicts(kept, contigs, gc_results)

    records_path = tmp / "records.tsv"
    ordinal = 0
    with open(records_path, "w", encoding="utf-8", buffering=1 << 20) as rf, open(
        tmp / "genetic_codes.tsv", "w", encoding="utf-8", buffering=1 << 20
    ) as gf:
        rf.write("ordinal\tgenome\trecord_ordinal\tcontig_uid\tlength\tassembly_quality\tgenetic_code\t"
                 "ena_accession\n")
        gf.write("\t".join(c58.GC_HEAD) + "\n")
        for sag in kept:
            raw, lens = contigs[sag]
            gc = gc_results[sag]
            code = resolved[sag]
            row = matched[sag]["row"]
            gf.write("\t".join([
                sag, STAGE, row["scientific_name"] if row else "", str(sum(lens)), str(len(lens)), "contig",
                str(code), gc["reason"], str(gc["tested_alt"]).lower(),
                c58._fmt_opt(gc["cd_11"]), c58._fmt_opt(gc["mean_orf_11"]), c58._fmt_opt(gc["n_genes_11"]),
                c58._fmt_opt(gc["cd_25"]), c58._fmt_opt(gc["mean_orf_25"]), c58._fmt_opt(gc["n_genes_25"]),
                c58._fmt_opt(gc["cd_4"]), c58._fmt_opt(gc["mean_orf_4"]), c58._fmt_opt(gc["n_genes_4"]),
                c58._fmt_opt(gc["tga_per_kb"])]) + "\n")
            for k, n in enumerate(lens):
                rf.write(f"{ordinal}\t{sag}\t{k}\t{raw[k * 16:(k + 1) * 16].hex()}\t{n}\tcontig\t{code}\t"
                         f"{matched[sag]['accession'] or ''}\n")
                ordinal += 1

    owners = c58._owners_from_records(records_path, root, root_stages, owner_inputs, work=work_dir,
                                      threads=threads, memory_limit=memory_limit)
    later = c58._owners_from_records(records_path, None, [], later_inputs, work=work_dir, threads=threads,
                                     memory_limit=memory_limit) if later_inputs else {}

    listing = make_inputs.Listing()
    genomes, repeated, shared, held_later = [], [], [], []
    first_of: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}
    ordinal = 0
    for sag in kept:
        raw, lens = contigs[sag]
        uids = [raw[i: i + 16] for i in range(0, len(raw), 16)]
        for k, (cu, n) in enumerate(zip(uids, lens)):
            contig_len.setdefault(cu, n)
            if cu in owners:
                shared.append((ordinal, sag, k, cu.hex(), n, owners[cu]))
            elif cu in later:
                held_later.append((ordinal, sag, k, cu.hex(), n, later[cu]))
            ordinal += 1
        g = listing.add(uids, resolved[sag], "prokaryote", "contig")
        guid = genome_uids[sag]
        if g is None:
            repeated.append((ordinal - len(uids), sag, guid.hex(), first_of[guid]))
            continue
        genomes.append(g)
        first_of[guid] = sag
    make_inputs._write_tsv(tmp / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
    make_inputs._write_tsv(tmp / "shared_records.tsv",
                           ["ordinal", "genome", "record_ordinal", "contig_uid", "length", "owner_stage"], shared)
    make_inputs._write_tsv(tmp / "later_records.tsv",
                           ["ordinal", "genome", "record_ordinal", "contig_uid", "length", "holder"], held_later)

    match_rows = []
    for sag in sorted(scanned):
        r, row = matched[sag], matched[sag]["row"]
        x = look.get(r["accession"]) if r["accession"] else None
        match_rows.append((
            sag, scanned[sag]["member"], scanned[sag]["contigs"], scanned[sag]["bp"], r["accession"] or "",
            row["assembly_name"] if row else "", row["base_count"] if row else "",
            row["scientific_name"] if row else "", r["basis"], _yn(x is not None),
            _yn(bool(x and x["in_gtdb_r232"])) if x else "", _yn(bool(x and x["in_globdb_r232"])) if x else "",
            "kept", ""))
    make_inputs._write_tsv(tmp / "sag_match.tsv", MATCH_HEAD, match_rows)
    status = {s: r232_status(look.get(matched[s]["accession"]) if matched[s]["accession"] else None) for s in kept}
    make_inputs._write_tsv(tmp / "r232_status.tsv", STATUS_HEAD, [
        (s, matched[s]["accession"] or "", matched[s]["basis"], status[s],
         _yn(bool(look.get(matched[s]["accession"], {}).get("in_gtdb_r232"))) if matched[s]["accession"] else "",
         str(look[matched[s]["accession"]]["is_representative"]) if matched[s]["accession"]
         and look[matched[s]["accession"]]["in_gtdb_r232"] else "",
         (look[matched[s]["accession"]]["representative"] or "") if matched[s]["accession"] else "",
         (look[matched[s]["accession"]]["species"] or "") if matched[s]["accession"] else "",
         _yn(bool(look.get(matched[s]["accession"], {}).get("in_globdb_r232"))) if matched[s]["accession"] else "",
         (look[matched[s]["accession"]]["globdb_r232_id"] or "") if matched[s]["accession"] else "",
         scanned[s]["contigs"], scanned[s]["bp"]) for s in kept])

    obj = {**SOURCE, "retrieved_at": retrieved_at}
    (tmp / "objects.json").write_text(json.dumps([{"path": str(inputs_home / STAGE / FASTA_NAME), **obj}], indent=1))
    (tmp / "objects_staging.json").write_text(json.dumps([{"path": str(final_fa), **obj}], indent=1))
    (tmp / "genomes.json").write_text(json.dumps(genomes))

    def tally(holders: dict[bytes, str], skip: dict[bytes, str] | None = None):
        c, b = Counter(), Counter()
        for cu, n in contig_len.items():
            if cu in holders and not (skip and cu in skip):
                c[holders[cu]] += 1
                b[holders[cu]] += n
        return dict(sorted(c.items())), dict(sorted(b.items()))

    own_c, own_b = tally(owners)
    lat_c, lat_b = tally(later, owners)
    matched_n = sum(1 for r in matched.values() if r["accession"])
    kept_matched = sum(1 for s in kept if matched[s]["accession"])
    codes = Counter(str(resolved[s]) for s in kept)
    return {
        "sags_in_zip": len(scanned), "contigs_in_zip": sum(d["contigs"] for d in scanned.values()),
        "bases_in_zip": sum(d["bp"] for d in scanned.values()),
        "match_basis": dict(sorted(Counter(r["basis"] for r in matched.values()).items())),
        "seen": {"matched_to_ena": matched_n, "not_matched": len(scanned) - matched_n,
                 "kept_after_test": kept_matched},
        "r232_status": dict(sorted(Counter(status.values()).items())),
        "r232_status_bases": {k: sum(scanned[s]["bp"] for s in kept if status[s] == k)
                              for k in sorted(set(status.values()))},
        "gtdb_representatives": sum(1 for s in kept if matched[s]["accession"]
                                    and look[matched[s]["accession"]]["is_representative"]),
        "sags_staged": len(kept), "bases_staged": sum(scanned[s]["bp"] for s in kept),
        "contigs_staged": sum(scanned[s]["contigs"] for s in kept), "contigs_renamed": n_renamed,
        "genomes": len(genomes), "genomes_repeated": len(repeated), "contigs_distinct": len(contig_len),
        "bases_distinct": sum(contig_len.values()),
        "genetic_codes": dict(sorted(codes.items())), "code_conflicts": len(conflicts),
        "held_by": own_c, "held_bases_by": own_b, "later_held_by": lat_c, "later_held_bases_by": lat_b,
        "contigs_new": sum(1 for cu in contig_len if cu not in owners),
        "bases_new": sum(n for cu, n in contig_len.items() if cu not in owners),
        "contigs_new_after_later": sum(1 for cu in contig_len if cu not in owners and cu not in later),
        "bases_new_after_later": sum(n for cu, n in contig_len.items() if cu not in owners and cu not in later),
        **{"listing_" + k: v for k, v in listing.report().items()},
    }


# --- command line --------------------------------------------------------------------------------------------------

def _check(meta: Path) -> dict:
    p = subprocess.run([sys.executable, str(Path(__file__).resolve().parent / "s2_check_inputs.py"), STAGE,
                        str(meta / "genomes.json"), str(meta / "objects_staging.json")],
                       capture_output=True, text=True, check=False)
    doc = json.loads(p.stdout) if p.stdout.strip().startswith("{") else {"stdout": p.stdout}
    return {"exit": p.returncode, "problems": doc.get("problems"), "gene_call_mode": doc.get("gene_call_mode"),
            "genomes": doc.get("genomes"), "genome_input_rows": doc.get("genome_input_rows"),
            "codes": doc.get("codes"), "kinds": doc.get("kinds"), "stderr": p.stderr[-2000:]}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Stage GORG-Dark SAG2 inputs (card C-77, DEC-79)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("census")
    c.add_argument("zip", type=Path)
    st = sub.add_parser("stage")
    for n in ("zip", "ena_tsv", "index_dir", "meta_out", "fasta_out"):
        st.add_argument(n, type=Path)
    st.add_argument("--root", type=Path, default=None)
    st.add_argument("--root-stages", default=",".join(ROOT_STAGES))
    st.add_argument("--chain-home", type=Path, default=None,
                    help="the r26.11.4 stages' inputs, read in a26's chain order as owners")
    st.add_argument("--later", action="append", default=[], metavar="NAME=DIR",
                    help="an r26.11.5-chain stage before SAG2 whose genomes.json is under DIR")
    st.add_argument("--inputs-home", type=Path, default=Path("/tmp/if2453/s2/inputs"))
    st.add_argument("--retrieved-at", default="2026-11-25T23:15:40Z")
    st.add_argument("--work", type=Path, default=None)
    st.add_argument("--threads", type=int, default=2)
    st.add_argument("--workers", type=int, default=1)
    st.add_argument("--memory-limit", default="4GB")
    st.add_argument("--report-out", type=Path, default=None)
    a = ap.parse_args(argv)
    try:
        if a.cmd == "census":
            rep = {**census(a.zip), "member_ids_by_prefix": dict(sorted(Counter(
                m.sag.split("-", 1)[0] for m in members(a.zip)).items()))}
            print(json.dumps(rep, indent=1, sort_keys=True))
            return 0
        owners = {s: a.chain_home / s for s, _ in a26.R26_11_4_CHAIN} if a.chain_home is not None else {}
        later = {}
        for item in a.later:
            k, _, v = item.partition("=")
            if not k or not v:
                raise Refused(f"--later {item!r}: expected NAME=DIR")
            later[k] = Path(v)
        rep = stage(a.zip, a.ena_tsv, a.index_dir, a.meta_out, a.fasta_out, root=a.root,
                    root_stages=[s for s in a.root_stages.split(",") if s], owner_inputs=owners, later_inputs=later,
                    inputs_home=a.inputs_home, retrieved_at=a.retrieved_at, work=a.work, threads=a.threads,
                    memory_limit=a.memory_limit, workers=a.workers)
        rep["owner_stages"] = [*([s for s in a.root_stages.split(",") if s] if a.root else []), *owners]
        rep["later_stages"] = [*later]
        rep["s2_check"] = _check(a.meta_out / STAGE)
        txt = json.dumps(rep, indent=1, sort_keys=True)
        if a.report_out is not None:
            a.report_out.parent.mkdir(parents=True, exist_ok=True)
            a.report_out.write_text(txt + "\n")
        print(txt)
        return 0 if rep["s2_check"]["exit"] == 0 and rep["s2_check"]["problems"] == [] else 1
    except Refused as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
