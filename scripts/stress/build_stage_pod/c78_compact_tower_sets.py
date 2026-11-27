"""Card C-78 (DEC-79, DEC-80): stage the compact sets that landed on `tower` in batch 2: `LRP1`, `MIR2`, `SYM2`,
`DPN1`. `MGE4` is not staged: its three FASTA files are isolate genome assemblies of one conjugation experiment
(the note says why); nothing of it reaches the pod.

    python c78_compact_tower_sets.py census PATH [PATH ...]
    python c78_compact_tower_sets.py lrp1 FASTA META_OUT FASTA_OUT [common options]

Common options: --root ROOT [--root-stages S1,P1,...] [--chain-home /tmp/if2453/s2/inputs] [--later NAME=DIR ...]
[--inputs-home /tmp/if2453/s2/inputs] [--retrieved-at T] [--work DIR] [--threads N] [--memory-limit 4GB]
[--report-out F].

`census` reads each archive (zip, tar or tar.gz) or plain FASTA and names every member: for a member whose name
looks like FASTA it counts records and bases and tells genome DNA (`dna`: every record of IUPAC nucleotide letters)
from `protein`, `alignment` (gap characters) and `not_fasta`. Nothing is passed over unnamed.

Every stage keeps what DEC-27 and `s2_check_inputs.py` keep: a record that is not DNA is skipped and listed; a
single-record genome whose sequence repeats an earlier one is skipped and listed (`skipped_records.tsv`); a multi-record
genome with the same contig set as another is listed once (`repeated_records.tsv`). Ownership: `shared_records.tsv`
the contigs the root's stored stages or the `r26.11.4` stages hold (the template's `EXPECT=held_by` lines),
`later_records.tsv` those only a stage of the `r26.11.5` chain before this one holds (`--later`); both are kept
(DEC-26: the first holder owns the contig's calls). No genome is removed because GTDB R232 or GlobDB r232 lists it.

LRP1 (`tower_lrp1`, Zenodo 13143815 `metagenome_plasmidome_plasmid2.fasta`): every DNA record is one genome,
`kind="plasmid"`, code 11. DEC-65 (1) (metagenome-assembled): an exact 20..1,000 bp terminal repeat is trimmed once and
the record staged circular (quality `complete`); a record whose header says `circular` (the assembler's
`_type_circular`) and has no terminal repeat is circular as it stands; every other record is linear (`contig`).
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import zipfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import a26_r26_11_4_chain as a26  # noqa: E402
import c48c_virus_inputs as v2  # noqa: E402
import c80_kavagutti_inphared as c80  # noqa: E402
import make_inputs  # noqa: E402
from biobox import fasta, ids  # noqa: E402

Refused = make_inputs.Refused

ROOT_STAGES = ("S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1")
VALID_DNA_RE = re.compile(rb"^[ACGTURYKMSWBDHVNacgturykmswbdhvn]+$")
FASTA_EXTS = (".fa", ".fna", ".fasta", ".fas", ".ffn", ".faa", ".fa.gz", ".fna.gz", ".fasta.gz", ".faa.gz")
CIRCULAR_RE = re.compile(r"circular", re.IGNORECASE)
NOT_A_GAP = re.compile(rb"[-.]")
PROTEIN_ONLY = set(b"EFILPQZJO*efilpqzjo")  # letters no nucleotide code uses


# --- census --------------------------------------------------------------------------------------------------------

def _looks_fasta(name: str) -> bool:
    base = os.path.basename(name.lower())
    return not base.startswith(".") and base.endswith(FASTA_EXTS)


def _alphabet_of(recs_seqs: list[bytes]) -> str:
    if any(NOT_A_GAP.search(s) for s in recs_seqs):
        return "alignment"
    if all(VALID_DNA_RE.match(s) for s in recs_seqs if s):
        return "dna"
    return "protein" if any(PROTEIN_ONLY & set(s) for s in recs_seqs) else "other"


def _fasta_stats(fobj, name: str, size: int, max_headers: int = 3) -> dict:
    n = bp = 0
    mn, mx, circ = None, 0, 0
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
                heads.append(r.header)
            kinds.add(_alphabet_of([r.seq]))
    except (ValueError, UnicodeDecodeError):
        return {"name": name, "bytes": size, "alphabet": "not_fasta"}
    if n == 0:
        return {"name": name, "bytes": size, "alphabet": "not_fasta"}
    alpha = "dna" if kinds == {"dna"} else next(k for k in ("alignment", "protein", "other", "dna") if k in kinds)
    return {"name": name, "bytes": size, "alphabet": alpha, "records": n, "bases": bp, "min_len": mn, "max_len": mx,
            "circular_in_header": circ, "first_headers": heads}


def _maybe_gz(fobj, name: str):
    return gzip.open(fobj, "rb") if name.lower().endswith(".gz") else fobj


def census(paths) -> dict:
    """What each archive or FASTA holds, member by member (module docstring)."""
    out: dict[str, dict] = {}
    for p in map(Path, paths):
        members: list[dict] = []
        if zipfile.is_zipfile(p):
            with zipfile.ZipFile(p) as zf:
                for i in zf.infolist():
                    if i.is_dir():
                        members.append({"name": i.filename, "bytes": 0, "alphabet": "directory"})
                    elif _looks_fasta(i.filename) and not i.filename.startswith("__MACOSX/"):
                        with zf.open(i) as raw:
                            members.append(_fasta_stats(_maybe_gz(raw, i.filename), i.filename, i.file_size))
                    else:
                        members.append({"name": i.filename, "bytes": i.file_size, "alphabet": "not_fasta"})
        elif tarfile.is_tarfile(p):
            with tarfile.open(p, "r:*") as tf:
                for m in tf:
                    if m.isdir():
                        members.append({"name": m.name.rstrip("/") + "/", "bytes": 0, "alphabet": "directory"})
                    elif m.isfile() and _looks_fasta(m.name):
                        f = tf.extractfile(m)
                        assert f is not None
                        members.append(_fasta_stats(_maybe_gz(f, m.name), m.name, m.size))
                    else:
                        members.append({"name": m.name, "bytes": m.size, "alphabet": "not_fasta"})
        else:
            with open(p, "rb") as fh:
                members.append(_fasta_stats(_maybe_gz(fh, p.name), p.name, p.stat().st_size))
        kinds = Counter(m["alphabet"] for m in members)
        out[str(p)] = {"members": members, "by_alphabet": dict(sorted(kinds.items())),
                       "dna_records": sum(m.get("records", 0) for m in members if m["alphabet"] == "dna"),
                       "dna_bases": sum(m.get("bases", 0) for m in members if m["alphabet"] == "dna")}
    return out


def extract(archive: Path, out_dir: Path, names: list[str]) -> dict:
    """Copy the named members of a tar or zip to OUT_DIR/<base name>, byte for byte (a table to read); refused, with
    nothing written, when a name is missing, two names share a base name, or a target exists."""
    out_dir = Path(out_dir)
    bases = [os.path.basename(n) for n in names]
    if len(set(bases)) != len(bases):
        raise Refused("two members share a base name")
    for b in bases:
        if (out_dir / b).exists():
            raise Refused(f"{out_dir / b} exists")
    want = set(names)
    got: dict[str, dict] = {}
    out_dir.mkdir(parents=True, exist_ok=True)
    made: list[Path] = []

    def put(name: str, src) -> None:
        dst = out_dir / os.path.basename(name)
        part = dst.with_name(f".{dst.name}.partial")
        made.append(part)
        with open(part, "wb") as w:
            shutil.copyfileobj(src, w, 1 << 22)
        got[name] = {"part": part, "dst": dst}

    try:
        if zipfile.is_zipfile(archive):
            with zipfile.ZipFile(archive) as zf:
                for i in zf.infolist():
                    if i.filename in want:
                        with zf.open(i) as src:
                            put(i.filename, src)
        else:
            with tarfile.open(archive, "r:*") as tf:
                for m in tf:
                    if m.isfile() and m.name in want:
                        f = tf.extractfile(m)
                        assert f is not None
                        put(m.name, f)
        missing = sorted(want - set(got))
        if missing:
            raise Refused(f"{archive}: no member {missing[0]}")
    except BaseException:
        for p in made:
            p.unlink(missing_ok=True)
        raise
    out = {}
    for name, d in sorted(got.items()):
        os.replace(d["part"], d["dst"])
        out[name] = {"path": str(d["dst"]), "sha256": ids.sha256_file(d["dst"]).hex(), "bytes": d["dst"].stat().st_size}
    return out


# --- DEC-65 (1) ----------------------------------------------------------------------------------------------------

def dec65(seq: bytes, *, flagged: bool) -> tuple[bytes, int, str, bool]:
    """DEC-65 (1) on a metagenome-assembled record: an exact 20..1,000 bp terminal repeat is trimmed once and the
    record is circular; a record its source flags circular with no repeat (class `none`) is circular untrimmed; every
    other record (and every other class) is linear and whole. Returns (sequence, bases trimmed, class, circular)."""
    e = v2.ends(seq)
    cls = v2.classify(e)
    if cls == "exact":
        return v2.trim(seq, e.exact), e.exact, cls, True
    return seq, 0, cls, flagged and cls == "none"


# --- the stage writer ----------------------------------------------------------------------------------------------

class StageFiles:
    """One FASTA per category under FASTA_OUT/<STAGE>/ (written as .partial), metadata in META_OUT/<STAGE>.tmp;
    `commit` renames everything into place, `abort` leaves nothing behind."""

    def __init__(self, stage: str, meta_out: Path, fasta_out: Path, fasta_names: dict[str, str], work: Path | None):
        self.stage = stage
        self.meta_dir, self.fasta_dir = Path(meta_out) / stage, Path(fasta_out) / stage
        self.fasta_names = fasta_names
        if self.meta_dir.exists():
            raise Refused(f"{self.meta_dir} exists")
        for n in fasta_names.values():
            if (self.fasta_dir / n).exists():
                raise Refused(f"{self.fasta_dir / n} exists")
        self.tmp = Path(meta_out) / f"{stage}.tmp"
        shutil.rmtree(self.tmp, ignore_errors=True)
        self.tmp.mkdir(parents=True)
        self.fasta_dir.mkdir(parents=True, exist_ok=True)
        self.work = Path(work) if work is not None else self.tmp / "_work"
        self.own_work = work is None
        self.work.mkdir(parents=True, exist_ok=True)
        self.partial = {k: self.fasta_dir / f".{n}.partial" for k, n in fasta_names.items()}
        self.fh = {k: open(p, "wb", buffering=1 << 22) for k, p in self.partial.items()}
        self.names: set[str] = set()

    def write(self, category: str, header: str, seq: bytes, width: int, prefix: str | None = None) -> str:
        """Write one record; a first token already used in the stage is given PREFIX (default: the category)."""
        first = header.split(None, 1)[0] if header.strip() else ""
        if not first or first in self.names:
            header = f"{prefix or category}|{header}"
            first = header.split(None, 1)[0]
        if first in self.names:
            raise Refused(f"contig name {first} occurs twice in stage {self.stage}")
        self.names.add(first)
        self.fh[category].write(b"\n".join(fasta.format_record(header, seq, width or 60)) + b"\n")
        return first

    def close(self) -> None:
        for f in self.fh.values():
            f.close()

    def abort(self) -> None:
        self.close()
        shutil.rmtree(self.tmp, ignore_errors=True)
        for p in self.partial.values():
            p.unlink(missing_ok=True)
        if self.fasta_dir.is_dir() and not any(self.fasta_dir.iterdir()):
            self.fasta_dir.rmdir()

    def commit(self, used: list[str]) -> tuple[dict, dict]:
        self.close()
        if self.own_work:
            shutil.rmtree(self.work, ignore_errors=True)
        parts = c80._meta_parts(self.tmp, self.meta_dir)
        for k, p in self.partial.items():
            if k in used:
                parts[self.fasta_dir / self.fasta_names[k]] = p
            else:
                p.unlink(missing_ok=True)
        c80._commit(parts)
        shutil.rmtree(self.tmp, ignore_errors=True)
        return c80._hashes([*parts])


SKIP_HEAD = ["category", "member", "name", "length", "reason"]


def _finish(sf: StageFiles, genomes: list[dict], trims: list[tuple], skipped: list[tuple], *, root, root_stages,
            owner_inputs, later_inputs, inputs_home, retrieved_at, source: dict, source_label: str, threads: int,
            memory_limit: str, extra: dict | None = None) -> dict:
    """Ownership, the stage's metadata (c80's writer), then everything into place."""
    if not genomes:
        raise Refused(f"{sf.stage}: no genome to stage")
    sf.close()
    all_uids = {c[0] for g in genomes for c in g["contigs"]}
    owners, later = c80._owners(all_uids, root, list(root_stages), owner_inputs, later_inputs, sf.work, threads,
                                memory_limit)
    make_inputs._write_tsv(sf.tmp / "skipped_records.tsv", SKIP_HEAD, skipped)
    used = [k for k in sf.fasta_names if any(g["category"] == k for g in genomes)]
    rep = c80._write_stage_meta(sf.stage, genomes, trims, owners, later, sf.tmp, Path(inputs_home), sf.fasta_dir,
                                retrieved_at, source, {k: sf.fasta_names[k] for k in used}, source_label)
    rep["skipped_records"] = dict(sorted(Counter(s[4].split(" ", 1)[0] for s in skipped).items()))
    rep.update(extra or {})
    rep["sha256"], rep["bytes"] = sf.commit(used)
    return rep


def _single_record_genome(sf: StageFiles, category: str, member: str, header: str, seq: bytes, width: int, *,
                          kind: str, flagged: bool, seen: dict[bytes, str], trims: list, skipped: list,
                          dtr_classes: Counter, quality_linear: str = "contig", sources: list | None = None,
                          refuse=None) -> dict | None:
    """One record as one genome under DEC-65 (1); None when it is skipped (not DNA, an exact repeat, or a contig
    `refuse(src_cu, cu)` names a reason for)."""
    name = header.split(None, 1)[0] if header.strip() else ""
    if not seq or not VALID_DNA_RE.match(seq):
        skipped.append((category, member, name, len(seq), "not_dna"))
        return None
    src_cu = ids.contig_uid(seq)
    out, k, cls, circ = dec65(seq, flagged=flagged)
    cu = ids.contig_uid(out) if k else src_cu
    if cu in seen:
        skipped.append((category, member, name, len(seq), f"exact_repeat_of {seen[cu]}"))
        return None
    why = refuse(src_cu, cu) if refuse is not None else None
    if why:
        skipped.append((category, member, name, len(seq), why))
        return None
    dtr_classes[cls] += 1
    seen[cu] = name
    written = sf.write(category, header, out, width, prefix=Path(member).stem)
    if sources is not None:
        sources.append((written, category, member, len(seq), len(out), str(circ).lower(), cls))
    if k:
        trims.append((category, written, name, src_cu.hex(), cu.hex(), len(seq), k, k, len(out)))
    return {"name": written, "category": category, "kind": kind, "quality": "complete" if circ else quality_linear,
            "circular": [cu] if circ else [], "code": 11, "contigs": [(cu, len(out), len(seq), name, cls, src_cu, k)]}


# --- LRP1 ----------------------------------------------------------------------------------------------------------

LRP1_FASTA = {"plasmidome": "lrp1_urban_water_plasmidome.fna"}
LRP1_SOURCE = {
    "source_name": "tower_lrp1",
    "source_release": "EU urban-water circular plasmidome, Zenodo record 13143815, metagenome_plasmidome_plasmid2.fasta "
                      "(tower batch 2); LRP1 staging by c78_compact_tower_sets.py (DEC-65 (1), DEC-79)",
    "license": "CC BY 4.0 as the job list states it (sources.yaml tower_lrp1); no provider record was read",
    "redistributable": False,  # DEC-80 (3): true only where a recorded licence permits it; none is confirmed
}


def stage_lrp1(fasta_path: Path, meta_out: Path, fasta_out: Path, *, root, root_stages, owner_inputs, later_inputs,
               inputs_home, retrieved_at: str = "2026-11-25T22:09:48Z", work=None, threads: int = 2,
               memory_limit: str = "4GB") -> dict:
    sf = StageFiles("LRP1", meta_out, fasta_out, LRP1_FASTA, work)
    try:
        genomes, trims, skipped = [], [], []
        seen: dict[bytes, str] = {}
        classes: Counter = Counter()
        n = flagged_n = 0
        member = Path(fasta_path).name
        with open(fasta_path, "rb") as fh:
            try:
                for r in v2.iter_fasta_stream(fh):
                    n += 1
                    flagged = bool(CIRCULAR_RE.search(r.header))
                    flagged_n += flagged
                    g = _single_record_genome(sf, "plasmidome", member, r.header, r.seq, r.line_width,
                                              kind="plasmid", flagged=flagged, seen=seen, trims=trims,
                                              skipped=skipped, dtr_classes=classes)
                    if g is not None:
                        genomes.append(g)
            except ValueError as exc:
                raise Refused(f"{fasta_path}: not a FASTA file ({exc})") from exc
        rep = _finish(sf, genomes, trims, skipped, root=root, root_stages=root_stages, owner_inputs=owner_inputs,
                      later_inputs=later_inputs, inputs_home=inputs_home, retrieved_at=retrieved_at,
                      source=LRP1_SOURCE, source_label="tower_lrp1", threads=threads, memory_limit=memory_limit,
                      extra={"records_read": n, "circular_flagged_in_header": flagged_n,
                             "dtr_classes": dict(sorted(classes.items())),
                             "records_under_5kb": sum(1 for g in genomes if g["contigs"][0][1] < 5000)})
    except BaseException:
        sf.abort()
        raise
    return rep


# --- MIR2 ----------------------------------------------------------------------------------------------------------

MIR2_FASTA = {"bulzu": "mir2_bulzu_2025_plv_virophage.fna", "bellas": "mir2_bellas_2026_plv_framework.fna"}
MIR2_SOURCE = {
    "source_name": "tower_mir2",
    "source_release": "Bulzu et al. 2025 cryptophyte PLV and virophage genomes (Zenodo 14540677 "
                      "genomes_polinton_virophages.zip) and Bellas et al. 2026 PLV framework (figshare file 65710359 "
                      "PLV_framework_supplementary_data.tar.gz: fasta_unsorted_PLV/ and annotated_PLVs/); MIR2 staging "
                      "by c78_compact_tower_sets.py (DEC-65 (1), DEC-79)",
    "license": "CC BY 4.0 as the job list states it (sources.yaml tower_mir2); no provider record was read",
    "redistributable": False,  # DEC-80 (3)
}
SOURCES_HEAD = ["genome", "category", "member", "source_length", "length", "circular", "dtr_class"]
MEMBER_HEAD = ["category", "member", "bytes", "reason"]


def _first_alphabet(fobj) -> tuple[str, object]:
    """The alphabet of a FASTA stream's first record, and an iterator over all its records (the first included)."""
    it = v2.iter_fasta_stream(fobj)
    try:
        first = next(it)
    except StopIteration:
        return "not_fasta", iter(())
    except (ValueError, UnicodeDecodeError):
        return "not_fasta", iter(())

    def chain():
        yield first
        yield from it
    return _alphabet_of([first.seq]), chain()


def _stage_members(sf, category: str, stream, genomes, trims, skipped, skipped_members, sources, seen, classes,
                   *, kind: str, flag_of=None, refuse=None) -> None:
    """Every FASTA member of an archive stream whose first record is DNA: one genome per record."""
    for name, size, fobj in stream:
        if fobj is None or not _looks_fasta(name):
            skipped_members.append((category, name, size, "directory" if fobj is None else "not_fasta"))
            continue
        alpha, recs = _first_alphabet(_maybe_gz(fobj, name))
        if alpha != "dna":
            skipped_members.append((category, name, size, alpha))
            continue
        try:
            for r in recs:
                g = _single_record_genome(sf, category, name, r.header, r.seq, r.line_width, kind=kind,
                                          flagged=bool(flag_of and flag_of(r.header)), seen=seen, trims=trims,
                                          skipped=skipped, dtr_classes=classes, sources=sources, refuse=refuse)
                if g is not None:
                    genomes.append(g)
        except (ValueError, UnicodeDecodeError) as exc:
            raise Refused(f"{name}: not a readable FASTA ({exc})") from exc


def _zip_stream(path: Path):
    with zipfile.ZipFile(path) as zf:
        for i in zf.infolist():
            if i.is_dir():
                yield i.filename, 0, None
            elif i.filename.startswith("__MACOSX/"):
                yield i.filename, i.file_size, None
            else:
                with zf.open(i) as raw:
                    yield i.filename, i.file_size, raw


def _tar_stream(path: Path):
    with tarfile.open(path, "r|*") as tf:
        for m in tf:
            if m.isfile():
                yield m.name, m.size, tf.extractfile(m)
            else:
                yield m.name, 0, None


def stage_mir2(genomes_zip: Path, framework_tgz: Path, meta_out: Path, fasta_out: Path, *, root, root_stages,
               owner_inputs, later_inputs, inputs_home, retrieved_at: str = "2026-11-25T22:10:19Z", work=None,
               threads: int = 2, memory_limit: str = "4GB") -> dict:
    """MIR2: the DNA FASTA members of both archives, one `kind="virus"` genome per record (DEC-65 (1)); the protist
    assemblies (`protist_genomes_assembled.tar.gz`) are not read here (card C-60)."""
    sf = StageFiles("MIR2", meta_out, fasta_out, MIR2_FASTA, work)
    try:
        genomes, trims, skipped, skipped_members, sources = [], [], [], [], []
        seen: dict[bytes, str] = {}
        classes: Counter = Counter()
        _stage_members(sf, "bulzu", _zip_stream(Path(genomes_zip)), genomes, trims, skipped, skipped_members,
                       sources, seen, classes, kind="virus")
        _stage_members(sf, "bellas", _tar_stream(Path(framework_tgz)), genomes, trims, skipped, skipped_members,
                       sources, seen, classes, kind="virus")
        make_inputs._write_tsv(sf.tmp / "skipped_members.tsv", MEMBER_HEAD, skipped_members)
        make_inputs._write_tsv(sf.tmp / "record_sources.tsv", SOURCES_HEAD, sources)
        rep = _finish(sf, genomes, trims, skipped, root=root, root_stages=root_stages, owner_inputs=owner_inputs,
                      later_inputs=later_inputs, inputs_home=inputs_home, retrieved_at=retrieved_at,
                      source=MIR2_SOURCE, source_label="tower_mir2", threads=threads, memory_limit=memory_limit,
                      extra={"dtr_classes": dict(sorted(classes.items())),
                             "members_skipped": dict(sorted(Counter(m[3] for m in skipped_members).items())),
                             "by_member": _by_member(sources)})
    except BaseException:
        sf.abort()
        raise
    return rep


def _by_member(sources: list[tuple]) -> dict:
    out: dict[str, dict] = {}
    for genome, cat, member, src_n, n, circ, cls in sources:
        e = out.setdefault(member, {"genomes": 0, "bases": 0})
        e["genomes"] += 1
        e["bases"] += n
    return dict(sorted(out.items()))


# --- SYM2 ----------------------------------------------------------------------------------------------------------

SYM2_FASTA = {"bin": "sym2_protist_symbiont_bins.fna", "virus": "sym2_protist_viruses.fna"}
SYM2_SOURCE = {
    "source_name": "tower_sym2",
    "source_release": "Schulz/Woyke protist symbiont study, NERSC portal ms_protistsymb: bins_fna.tar.gz and "
                      "viruses.fna.tar.gz (assemblies_fna.tar.gz, the host assemblies, not staged: card C-60); SYM2 "
                      "staging by c78_compact_tower_sets.py (DEC-65 (1), DEC-79, DEC-80)",
    "license": "none stated: open-public-NERSC is the job list's access label, not a licence (sources.yaml tower_sym2)",
    "redistributable": False,  # DEC-80 (3): no recorded licence
}
QUALITY_OF_MIMAG = {"HQ": "mag_hq", "MQ": "mag_mq"}
BIN_HEAD = ["genome", "member", "contigs", "bp", "in_metadata", "mimag", "completeness_checkm2",
            "contamination_checkm2", "classification", "symbclass", "host", "genetic_code", "reason"]
VIRUS_HEAD = ["genome", "length", "staged_length", "circular", "in_final_set", "topology", "checkv_quality",
              "taxonomy", "host"]
STATUS_HEAD = ["genome", "accession", "r232_status", "in_gtdb_r232", "is_gtdb_representative", "gtdb_species",
               "in_globdb_r232", "globdb_r232_id", "contigs", "bp"]


def _read_tsv(path: Path, key: str) -> dict[str, dict]:
    import csv

    with open(path, encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if key not in (rd.fieldnames or []):
            raise Refused(f"{path}: no column {key}")
        out: dict[str, dict] = {}
        for r in rd:
            out.setdefault(r[key], r)
    return out


def _screen_default(seqs: list[bytes], organism_name: str = "") -> dict:
    import c58_sags_symbionts_niche as c58

    return c58.screen_genetic_code(seqs, organism_name)


def _screen_job(args):
    name, seqs, org = args
    return name, _screen_default(seqs, org)


def stage_sym2(bins_tgz: Path, viruses_tgz: Path, tables: Path, meta_out: Path, fasta_out: Path, *, root,
               root_stages, owner_inputs, later_inputs, inputs_home, retrieved_at: str = "2026-11-25T22:26:57Z",
               work=None, threads: int = 2, memory_limit: str = "4GB", screen=None, workers: int = 2) -> dict:
    """SYM2: each bin of `bins_fna/` is one `kind="prokaryote"` genome, its code screened (11, 25, 4) as SYM1's with
    the authors' GTDB classification as the hint, quality from MIMAG (HQ, MQ; else `contig`); each record of
    `viruses.fna` is one `kind="virus"` genome under DEC-65 (1) (the authors' `topology` DTR is the source's circular
    flag), except a contig a bin already holds (kinds must agree, `s2_check_inputs.py`). The bins carry no INSDC
    accession: every one is `not_visible` to the r232 accession test."""
    import concurrent.futures

    meta = _read_tsv(Path(tables) / "metadata_mags.tab", "MAGid")
    final = _read_tsv(Path(tables) / "final_virus_df.tsv", "contig_id")
    sf = StageFiles("SYM2", meta_out, fasta_out, SYM2_FASTA, work)
    screen_fn = screen or _screen_default
    pool = concurrent.futures.ProcessPoolExecutor(max_workers=workers) if workers > 0 and screen is None else None
    try:
        genomes, trims, skipped, skipped_members, sources = [], [], [], [], []
        bin_rows, code_rows = [], []
        pending: list = []
        bin_of: dict[bytes, str] = {}
        seen: dict[bytes, str] = {}
        classes: Counter = Counter()

        def settle(g, res):
            g["code"], g["reason"] = int(res["genetic_code"]), res.get("reason", "")
            code_rows.append((g["name"], g["code"], g["reason"], res.get("cd_11"), res.get("mean_orf_11"),
                              res.get("cd_25"), res.get("cd_4")))

        for name, size, fobj in _tar_stream(Path(bins_tgz)):
            if fobj is None or not _looks_fasta(name):
                skipped_members.append(("bin", name, size, "directory" if fobj is None else "not_fasta"))
                continue
            gid = Path(name).name
            for ext in (".fna", ".fasta", ".fa"):
                if gid.endswith(ext):
                    gid = gid[: -len(ext)]
            recs = []
            for r in v2.iter_fasta_stream(fobj):
                if r.seq and VALID_DNA_RE.match(r.seq):
                    recs.append(r)
                else:
                    skipped.append(("bin", name, r.accession, len(r.seq), "not_dna"))
            if not recs:
                skipped_members.append(("bin", name, size, "no_dna_record"))
                continue
            contigs = []
            for r in recs:
                written = sf.write("bin", r.header, r.seq, r.line_width, prefix=gid)
                cu = ids.contig_uid(r.seq)
                bin_of.setdefault(cu, gid)
                contigs.append((cu, len(r.seq), len(r.seq), written, "", cu, 0))
            m = meta.get(gid, {})
            g = {"name": gid, "category": "bin", "kind": "prokaryote",
                 "quality": QUALITY_OF_MIMAG.get(m.get("MIMAG", ""), "contig"), "circular": [], "contigs": contigs,
                 "code": 11, "reason": ""}
            genomes.append(g)
            bin_rows.append([gid, name, len(contigs), sum(c[1] for c in contigs), "yes" if m else "no",
                             m.get("MIMAG", ""), m.get("Completeness_CheckM2", ""), m.get("Contamination_CheckM2", ""),
                             m.get("classification", ""), m.get("symbclass", ""), m.get("Host", "")])
            seqs, org = [r.seq for r in recs], m.get("classification", "")
            if pool is not None:
                pending.append((g, pool.submit(_screen_job, (gid, seqs, org))))
                while len(pending) > 2 * workers:
                    gg, fut = pending.pop(0)
                    settle(gg, fut.result()[1])
            else:
                settle(g, screen_fn(seqs, org))
        while pending:
            gg, fut = pending.pop(0)
            settle(gg, fut.result()[1])
        for row in bin_rows:
            g = next(x for x in genomes if x["name"] == row[0])
            row.extend([g["code"], g["reason"]])

        def in_bin(src_cu, cu):
            hit = bin_of.get(src_cu) or bin_of.get(cu)
            return f"contig_in_bin {hit}" if hit else None

        virus_rows = []
        n_vir = 0
        for name, size, fobj in _tar_stream(Path(viruses_tgz)):
            if fobj is None or not _looks_fasta(name):
                skipped_members.append(("virus", name, size, "directory" if fobj is None else "not_fasta"))
                continue
            for r in v2.iter_fasta_stream(fobj):
                n_vir += 1
                f = final.get(r.accession)
                g = _single_record_genome(sf, "virus", name, r.header, r.seq, r.line_width, kind="virus",
                                          flagged=bool(f and f.get("topology") == "DTR"), seen=seen, trims=trims,
                                          skipped=skipped, dtr_classes=classes, sources=sources, refuse=in_bin)
                if g is None:
                    continue
                genomes.append(g)
                c = g["contigs"][0]
                virus_rows.append((g["name"], c[2], c[1], str(bool(g["circular"])).lower(), "yes" if f else "no",
                                   (f or {}).get("topology", ""), (f or {}).get("checkv_quality", ""),
                                   (f or {}).get("taxonomy", ""), (f or {}).get("host", "")))
        bin_names = {g["name"] for g in genomes if g["category"] == "bin"}
        make_inputs._write_tsv(sf.tmp / "sym2_bins.tsv", BIN_HEAD, [tuple(r) for r in sorted(bin_rows)])
        make_inputs._write_tsv(sf.tmp / "sym2_viruses.tsv", VIRUS_HEAD, virus_rows)
        make_inputs._write_tsv(sf.tmp / "genetic_codes.tsv",
                               ["genome", "genetic_code", "reason", "cd_11", "mean_orf_11", "cd_25", "cd_4"],
                               sorted(tuple("" if v is None else v for v in r) for r in code_rows))
        make_inputs._write_tsv(sf.tmp / "skipped_members.tsv", MEMBER_HEAD, skipped_members)
        make_inputs._write_tsv(sf.tmp / "r232_status.tsv", STATUS_HEAD, [
            (g["name"], "", "not_visible", "", "", "", "", "", len(g["contigs"]), sum(c[1] for c in g["contigs"]))
            for g in sorted(genomes, key=lambda x: x["name"]) if g["category"] == "bin"])
        genomes.sort(key=lambda g: (g["category"], g["name"]))
        in_final = sum(1 for r in virus_rows if r[4] == "yes")
        rep = _finish(sf, genomes, trims, skipped, root=root, root_stages=root_stages, owner_inputs=owner_inputs,
                      later_inputs=later_inputs, inputs_home=inputs_home, retrieved_at=retrieved_at,
                      source=SYM2_SOURCE, source_label="tower_sym2", threads=threads, memory_limit=memory_limit,
                      extra={"bins": len(bin_names), "bins_in_metadata": sum(1 for r in bin_rows if r[4] == "yes"),
                             "bins_by_mimag": dict(sorted(Counter(r[5] or "none" for r in bin_rows).items())),
                             "genetic_codes": dict(sorted(Counter(str(g["code"]) for g in genomes
                                                                  if g["category"] == "bin").items())),
                             "virus_records_read": n_vir, "viruses_staged": len(virus_rows),
                             "viruses_in_final_set": in_final, "final_set_rows": len(final),
                             "dtr_classes": dict(sorted(classes.items())),
                             "r232_status": {"not_visible": len(bin_names)}})
    except BaseException:
        sf.abort()
        raise
    finally:
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)
    return rep


# --- command line --------------------------------------------------------------------------------------------------

def _check(stage_name: str, meta: Path) -> dict:
    return c80._check(stage_name, meta)


def _kv(items: list[str]) -> dict[str, Path]:
    out = c80._kv(items)
    for name, d in out.items():
        if not (Path(d) / "genomes.json").is_file():
            raise Refused(f"--later {name}: no genomes.json in {d}")
    return out


def _common(p: argparse.ArgumentParser) -> None:
    p.add_argument("meta_out", type=Path)
    p.add_argument("fasta_out", type=Path)
    p.add_argument("--root", type=Path, default=None)
    p.add_argument("--root-stages", default=",".join(ROOT_STAGES))
    p.add_argument("--chain-home", type=Path, default=None,
                   help="the r26.11.4 stages' inputs, read in a26's chain order as owners")
    p.add_argument("--later", action="append", default=[], help="NAME=DIR of an r26.11.5 stage before this one")
    p.add_argument("--inputs-home", type=Path, default=Path("/tmp/if2453/s2/inputs"))
    p.add_argument("--retrieved-at", default=None)
    p.add_argument("--work", type=Path, default=None)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--memory-limit", default="4GB")
    p.add_argument("--report-out", type=Path, default=None)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Stage the compact tower sets of card C-78 (DEC-79)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("census")
    c.add_argument("paths", type=Path, nargs="+")
    c.add_argument("--out", type=Path, default=None)
    x = sub.add_parser("extract")
    x.add_argument("archive", type=Path)
    x.add_argument("out_dir", type=Path)
    x.add_argument("members", nargs="+")
    p = sub.add_parser("lrp1")
    p.add_argument("fasta", type=Path)
    _common(p)
    p = sub.add_parser("mir2")
    p.add_argument("genomes_zip", type=Path)
    p.add_argument("framework_tgz", type=Path)
    _common(p)
    p = sub.add_parser("sym2")
    p.add_argument("bins_tgz", type=Path)
    p.add_argument("viruses_tgz", type=Path)
    p.add_argument("tables", type=Path)
    _common(p)
    p.add_argument("--workers", type=int, default=2)
    a = ap.parse_args(argv)
    try:
        if a.cmd == "census":
            txt = json.dumps(census(a.paths), indent=1, sort_keys=True)
            if a.out is not None:
                a.out.parent.mkdir(parents=True, exist_ok=True)
                a.out.write_text(txt + "\n")
            print(txt)
            return 0
        if a.cmd == "extract":
            print(json.dumps(extract(a.archive, a.out_dir, a.members), indent=1, sort_keys=True))
            return 0
        owners = {s: a.chain_home / s for s, _ in a26.R26_11_4_CHAIN} if a.chain_home is not None else {}
        for s, d in owners.items():
            if not (d / "genomes.json").is_file():
                raise Refused(f"--chain-home: no genomes.json for {s} in {d}")
        kw = dict(root=a.root, root_stages=[s for s in a.root_stages.split(",") if s] if a.root else [],
                  owner_inputs=owners, later_inputs=_kv(a.later), inputs_home=a.inputs_home, work=a.work,
                  threads=a.threads, memory_limit=a.memory_limit)
        if a.retrieved_at:
            kw["retrieved_at"] = a.retrieved_at
        if a.cmd == "lrp1":
            stage_name = "LRP1"
            rep = stage_lrp1(a.fasta, a.meta_out, a.fasta_out, **kw)
        elif a.cmd == "mir2":
            stage_name = "MIR2"
            rep = stage_mir2(a.genomes_zip, a.framework_tgz, a.meta_out, a.fasta_out, **kw)
        elif a.cmd == "sym2":
            stage_name = "SYM2"
            rep = stage_sym2(a.bins_tgz, a.viruses_tgz, a.tables, a.meta_out, a.fasta_out, workers=a.workers, **kw)
    except Refused as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2
    rep["owner_stages"] = [*kw["root_stages"], *owners]
    rep["later_stages"] = [*kw["later_inputs"]]
    rep["s2_check"] = _check(stage_name, a.meta_out / stage_name)
    txt = json.dumps(rep, indent=1, sort_keys=True, default=str)
    if a.report_out is not None:
        a.report_out.parent.mkdir(parents=True, exist_ok=True)
        a.report_out.write_text(txt + "\n")
    print(txt)
    return 0 if rep["s2_check"]["exit"] == 0 and rep["s2_check"]["problems"] == [] else 1


if __name__ == "__main__":
    sys.exit(main())
