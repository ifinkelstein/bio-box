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

STOR_MIN_BYTES = 64 * 1024 * 1024  # AGENTS §7.3: no file under 64 MiB goes onto /stor


class StageFiles:
    """One FASTA per category, written first under META_OUT/<STAGE>.fasta.tmp/ (BGZF when its name ends `.gz`),
    metadata in META_OUT/<STAGE>.tmp. `commit` puts a FASTA of STOR_MIN_BYTES and more in FASTA_OUT/<STAGE>/ (copied
    as .partial, then renamed) and a smaller one beside the metadata in META_OUT/<STAGE>/, then renames the metadata
    into place; `abort` leaves nothing behind."""

    def __init__(self, stage: str, meta_out: Path, fasta_out: Path, fasta_names: dict[str, str], work: Path | None,
                 reuse_fasta: bool = False):
        """REUSE_FASTA: a FASTA already at its place on /stor is kept, never overwritten, when the newly written one has
        the same bytes (a rerun that changes only the metadata); different bytes are refused."""
        import c77_sag2_inputs as c77

        self.stage = stage
        self.meta_dir, self.fasta_dir = Path(meta_out) / stage, Path(fasta_out) / stage
        self.fasta_names = fasta_names
        self.reuse_fasta = reuse_fasta
        self.reused: list[str] = []
        if self.meta_dir.exists():
            raise Refused(f"{self.meta_dir} exists")
        for n in fasta_names.values():
            if (self.fasta_dir / n).exists() and not reuse_fasta:
                raise Refused(f"{self.fasta_dir / n} exists")
        self.tmp = Path(meta_out) / f"{stage}.tmp"
        self.ftmp = Path(meta_out) / f"{stage}.fasta.tmp"
        for d in (self.tmp, self.ftmp):
            shutil.rmtree(d, ignore_errors=True)
            d.mkdir(parents=True)
        self.work = Path(work) if work is not None else self.tmp / "_work"
        self.own_work = work is None
        self.work.mkdir(parents=True, exist_ok=True)
        self.partial = {k: self.ftmp / n for k, n in fasta_names.items()}
        self.fh = {k: (c77.BgzfWriter(p) if p.name.endswith(".gz") else open(p, "wb", buffering=1 << 22))
                   for k, p in self.partial.items()}
        self.names: set[str] = set()
        self.closed = False

    def final_path(self, category: str) -> Path:
        """Where a written FASTA goes (call after `close`)."""
        p = self.partial[category]
        big = p.stat().st_size >= STOR_MIN_BYTES
        return (self.fasta_dir if big else self.meta_dir) / self.fasta_names[category]

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
        if not self.closed:
            for f in self.fh.values():
                f.close()
            self.closed = True

    def abort(self) -> None:
        self.close()
        shutil.rmtree(self.tmp, ignore_errors=True)
        shutil.rmtree(self.ftmp, ignore_errors=True)
        for n in self.fasta_names.values():
            (self.fasta_dir / f".{n}.partial").unlink(missing_ok=True)

    def commit(self, used: list[str]) -> tuple[dict, dict]:
        self.close()
        if self.own_work:
            shutil.rmtree(self.work, ignore_errors=True)
        parts = c80._meta_parts(self.tmp, self.meta_dir)
        to_stor: dict[Path, Path] = {}
        for k in used:
            final = self.final_path(k)
            (to_stor if final.parent == self.fasta_dir else parts)[final] = self.partial[k]
        for final in parts:
            if final.exists():
                raise Refused(f"{final} exists")
        for final, src in list(to_stor.items()):
            if not final.exists():
                continue
            if not self.reuse_fasta or ids.sha256_file(final) != ids.sha256_file(src):
                raise Refused(f"{final} exists and differs from the newly written FASTA")
            self.reused.append(final)
            del to_stor[final]
            parts.pop(final, None)
        for final, src in to_stor.items():
            final.parent.mkdir(parents=True, exist_ok=True)
            part = final.with_name(f".{final.name}.partial")
            shutil.copyfile(src, part)
            if ids.sha256_file(part) != ids.sha256_file(src):
                part.unlink()
                raise Refused(f"{final}: the copy to /stor differs from the written file")
        for final in to_stor:
            os.replace(final.with_name(f".{final.name}.partial"), final)
        c80._commit(parts)
        shutil.rmtree(self.tmp, ignore_errors=True)
        shutil.rmtree(self.ftmp, ignore_errors=True)
        return c80._hashes([*parts, *to_stor, *self.reused])


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
    # objects_staging.json names each FASTA where commit puts it (/stor for 64 MiB and more, else beside the metadata)
    where = {sf.fasta_names[k]: str(sf.final_path(k)) for k in used}
    objs = json.loads((sf.tmp / "objects_staging.json").read_text())
    for o in objs:
        o["path"] = where[Path(o["path"]).name]
    (sf.tmp / "objects_staging.json").write_text(json.dumps(objs, indent=1))
    rep["skipped_records"] = dict(sorted(Counter(s[4].split(" ", 1)[0] for s in skipped).items()))
    rep.update(extra or {})
    rep["sha256"], rep["bytes"] = sf.commit(used)
    if sf.reuse_fasta:
        rep["fasta_reused"] = sorted(p.name for p in sf.reused)
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
GIANT_VIRUS_CLASSES = {"nucleocytoviricota"}  # the authors' symbclass of a giant-virus bin (metadata_mags.tab)
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
               work=None, threads: int = 2, memory_limit: str = "4GB", screen=None, workers: int = 2,
               reuse_fasta: bool = False) -> dict:
    """SYM2: each bin of `bins_fna/` is one `kind="prokaryote"` genome, its code screened (11, 25, 4) as SYM1's with
    the authors' GTDB classification as the hint, quality from MIMAG (HQ, MQ; else `contig`) — except a bin the
    authors' `symbclass` calls Nucleocytoviricota, a giant-virus genome: `kind="virus"`, not screened; each record of
    `viruses.fna` is one `kind="virus"` genome under DEC-65 (1) (the authors' `topology` DTR is the source's circular
    flag), except a contig a bin already holds (kinds must agree, `s2_check_inputs.py`). The bins carry no INSDC
    accession: every one is `not_visible` to the r232 accession test."""
    import concurrent.futures

    meta = _read_tsv(Path(tables) / "metadata_mags.tab", "MAGid")
    final = _read_tsv(Path(tables) / "final_virus_df.tsv", "contig_id")
    sf = StageFiles("SYM2", meta_out, fasta_out, SYM2_FASTA, work, reuse_fasta=reuse_fasta)
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
            giant_virus = m.get("symbclass", "").lower() in GIANT_VIRUS_CLASSES
            g = {"name": gid, "category": "bin", "kind": "virus" if giant_virus else "prokaryote",
                 "quality": "contig" if giant_virus else QUALITY_OF_MIMAG.get(m.get("MIMAG", ""), "contig"),
                 "circular": [], "contigs": contigs, "code": 11, "reason": "gv-meta decides" if giant_virus else ""}
            genomes.append(g)
            bin_rows.append([gid, name, len(contigs), sum(c[1] for c in contigs), "yes" if m else "no",
                             m.get("MIMAG", ""), m.get("Completeness_CheckM2", ""), m.get("Contamination_CheckM2", ""),
                             m.get("classification", ""), m.get("symbclass", ""), m.get("Host", "")])
            seqs, org = [r.seq for r in recs], m.get("classification", "")
            if giant_virus:
                continue  # a virus genome: no code screen (pyrodigal-gv calls it by its kind)
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
        prok = {g["name"] for g in genomes if g["category"] == "bin" and g["kind"] == "prokaryote"}
        make_inputs._write_tsv(sf.tmp / "sym2_bins.tsv", BIN_HEAD, [tuple(r) for r in sorted(bin_rows)])
        make_inputs._write_tsv(sf.tmp / "sym2_viruses.tsv", VIRUS_HEAD, virus_rows)
        make_inputs._write_tsv(sf.tmp / "genetic_codes.tsv",
                               ["genome", "genetic_code", "reason", "cd_11", "mean_orf_11", "cd_25", "cd_4"],
                               sorted(tuple("" if v is None else v for v in r) for r in code_rows))
        make_inputs._write_tsv(sf.tmp / "skipped_members.tsv", MEMBER_HEAD, skipped_members)
        make_inputs._write_tsv(sf.tmp / "r232_status.tsv", STATUS_HEAD, [
            (g["name"], "", "not_visible", "", "", "", "", "", len(g["contigs"]), sum(c[1] for c in g["contigs"]))
            for g in sorted(genomes, key=lambda x: x["name"]) if g["name"] in prok and g["category"] == "bin"])
        genomes.sort(key=lambda g: (g["category"], g["name"]))
        in_final = sum(1 for r in virus_rows if r[4] == "yes")
        rep = _finish(sf, genomes, trims, skipped, root=root, root_stages=root_stages, owner_inputs=owner_inputs,
                      later_inputs=later_inputs, inputs_home=inputs_home, retrieved_at=retrieved_at,
                      source=SYM2_SOURCE, source_label="tower_sym2", threads=threads, memory_limit=memory_limit,
                      extra={"bins": len(bin_names), "bins_in_metadata": sum(1 for r in bin_rows if r[4] == "yes"),
                             "bins_by_mimag": dict(sorted(Counter(r[5] or "none" for r in bin_rows).items())),
                             "bins_by_kind": dict(sorted(Counter(g["kind"] for g in genomes
                                                                 if g["category"] == "bin").items())),
                             "genetic_codes": dict(sorted(Counter(str(g["code"]) for g in genomes
                                                                  if g["name"] in prok
                                                                  and g["category"] == "bin").items())),
                             "virus_records_read": n_vir, "viruses_staged": len(virus_rows),
                             "viruses_in_final_set": in_final, "final_set_rows": len(final),
                             "dtr_classes": dict(sorted(classes.items())),
                             "r232_status": {"not_visible": len(prok)}})
    except BaseException:
        sf.abort()
        raise
    finally:
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)
    return rep


# --- DPN1 ----------------------------------------------------------------------------------------------------------
#
# Four steps, so that skani's run is its own unit:
#   dpn1-extract  the genome FASTA of the two archives, one plain file per genome under WORK/genomes/<category>/, and
#                 WORK/genomes.tsv (name, category, file, member, records, bp, accession, CheckM2, taxonomy)
#   dpn1-refs     card C-75's r226 archaea split one genome at a time (in batches, sketched by `skani sketch
#                 --separate-sketches` and deleted), WORK/ref_records.tsv (each record's contig_uid) and
#                 WORK/ref_holders.tsv: the stage that holds each reference (exact contig_uid over the root, the
#                 r26.11.4 stages and the r26.11.5 stages before DPN1); a reference is held when one stage holds at
#                 least half its bases
#   dpn1-skani    `skani dist` of every DPN1 genome against the reference sketches, `skani triangle -E` among them
#   dpn1          the decisions and the stage: a genome whose best hit at ANI >= 96 % and aligned fraction (the larger
#                 of the two) >= 50 % is a HELD reference is removed (the card: deduplicated against C-75's r226
#                 archaea); a hit on a reference no stage holds is recorded and the genome kept; then within the stage,
#                 greedily in rank order (the closed Nanobdellota genomes first, then CheckM2 completeness - 5 x
#                 contamination, then length, then name), a genome within 96 % / 50 % of one kept before it is
#                 removed. The 26 MGE contigs are `kind="virus"`: their metadata table has no class column.

DPN1_FASTA = {"dpann": "dpn1_dpann_genomes.fna.gz", "mge": "dpn1_dpann_mges.fna"}
DPN1_SOURCE = {
    "source_name": "tower_dpn1",
    "source_release": "Nanobdellota genomes of Zenodo 20174424 (Nanobdellota_Zen.zip: genomes/complete_208, "
                      "genomes/noncircular_51, ncbi_all/hq_genomes) and the 3,662 DPANN archaeal genomes and DPANN "
                      "viruses and unclassified MGEs of Zenodo 11215303; DPN1 staging by c78_compact_tower_sets.py "
                      "(skani 96 % ANI against card C-75's r226 archaea and within the stage; DEC-65 (1); DEC-79)",
    "license": "CC BY 4.0 as the job list states it (sources.yaml tower_dpn1); no provider record was read",
    "redistributable": False,  # DEC-80 (3)
}
DPN1_ZIP_PARTS = {"Nanobdellota_Zen/genomes/complete_208/": "nanobdellota_closed",
                  "Nanobdellota_Zen/genomes/noncircular_51/": "nanobdellota_open",
                  "Nanobdellota_Zen/ncbi_all/hq_genomes/": "ncbi_hq"}
DPN1_TAR_PART = ("3662.DPANN.genomes/", "dpann3662")
MANIFEST_HEAD = ["name", "category", "file", "member", "records", "bp", "accession", "completeness",
                 "contamination", "taxonomy"]
ACC_RE = re.compile(r"(GC[AF]_\d{9}(?:\.\d+)?)")
ANI_MIN, AF_MIN = 96.0, 50.0


def _gid(member: str) -> str:
    base = os.path.basename(member)
    for ext in (".fna", ".fasta", ".fa"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return base


def dpn1_extract(nano_zip: Path, dpann_tgz: Path, dpann_meta: Path, work: Path) -> dict:
    """One plain FASTA per genome under WORK/genomes/<category>/ and WORK/genomes.tsv; nothing else of the zip."""
    import csv

    work = Path(work)
    out_dir = work / "genomes"
    if out_dir.exists() or (work / "genomes.tsv").exists():
        raise Refused(f"{out_dir} exists")
    meta: dict[str, dict] = {}
    with open(dpann_meta, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            meta[r["Genome"]] = r
    tmp = work / "genomes.partial"
    shutil.rmtree(tmp, ignore_errors=True)
    rows: list[tuple] = []
    skipped: Counter = Counter()

    def put(category: str, member: str, fobj) -> None:
        gid = _gid(member)
        dst = tmp / category / f"{gid}.fna"
        if dst.exists():
            raise Refused(f"two members give genome {category}/{gid}")
        dst.parent.mkdir(parents=True, exist_ok=True)
        n = bp = 0
        with open(dst, "wb") as w:
            for r in v2.iter_fasta_stream(fobj):
                if not r.seq or not VALID_DNA_RE.match(r.seq):
                    raise Refused(f"{member}: record {r.accession} is not DNA")
                w.write(b"\n".join(fasta.format_record(r.header, r.seq, r.line_width or 60)) + b"\n")
                n += 1
                bp += len(r.seq)
        if not n:
            raise Refused(f"{member}: no record")
        m = meta.get(gid, {})
        acc = ACC_RE.search(gid)
        rows.append((gid, category, str(out_dir / category / f"{gid}.fna"), member, n, bp,
                     acc.group(1) if acc else "", m.get("Completeness (Checkm2)", ""),
                     m.get("Contamination (Checkm2)", ""), m.get("Taxonomy (GTDB r207)", "")))

    try:
        with zipfile.ZipFile(nano_zip) as zf:
            for i in zf.infolist():
                part = next((c for p, c in DPN1_ZIP_PARTS.items() if i.filename.startswith(p)), None)
                if part is None or i.is_dir() or not _looks_fasta(i.filename):
                    skipped["zip_member_not_genome_fasta"] += 0 if i.is_dir() else 1
                    continue
                with zf.open(i) as raw:
                    put(part, i.filename, raw)
        with tarfile.open(dpann_tgz, "r|*") as tf:
            for m in tf:
                if not m.isfile():
                    continue
                if not (m.name.startswith(DPN1_TAR_PART[0]) and _looks_fasta(m.name)):
                    skipped["tar_member_not_genome_fasta"] += 1
                    continue
                put(DPN1_TAR_PART[1], m.name, tf.extractfile(m))
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    rows.sort(key=lambda r: (r[1], r[0]))
    os.replace(tmp, out_dir)
    make_inputs._write_tsv(work / "genomes.tsv", MANIFEST_HEAD, rows)
    (work / "queries.txt").write_text("".join(r[2] + "\n" for r in rows))
    by = Counter(r[1] for r in rows)
    return {"genomes": len(rows), "by_category": dict(sorted(by.items())),
            "bases_by_category": {c: sum(r[5] for r in rows if r[1] == c) for c in sorted(by)},
            "dpann3662_in_metadata": sum(1 for r in rows if r[1] == "dpann3662" and r[7] != ""),
            "metadata_rows": len(meta), "members_skipped": dict(skipped)}


def _run(cmd: list, log: list) -> None:
    import time

    t0 = time.time()
    p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, check=False)
    log.append({"cmd": [str(c) for c in cmd[:2]], "wall_s": round(time.time() - t0, 1), "rc": p.returncode})
    if p.returncode != 0:
        raise Refused(f"{cmd[0]} {cmd[1]} failed (rc={p.returncode}): {p.stderr[-500:]}")


def dpn1_refs(archaea_fna: Path, archaea_list: Path, work: Path, *, root, root_stages, owner_inputs, later_inputs,
              skani: str | None = "skani", threads: int = 8, batch: int = 1000, memory_limit: str = "8GB") -> dict:
    """Split the reference FASTA by its list (records in list order, each genome contiguous), sketch each batch, and
    find the stage that holds each reference (module notes above)."""
    import csv

    import c62_metavr_inputs as c62

    work = Path(work)
    refs = work / "refs"
    if refs.exists():
        raise Refused(f"{refs} exists")
    with open(archaea_list, encoding="utf-8", newline="") as fh:
        lst = sorted(csv.DictReader(fh, delimiter="\t"), key=lambda r: int(r["first_record"]))
    refs.mkdir(parents=True)
    log: list[dict] = []
    sketches: list[str] = []
    rec_path = work / "ref_records.tsv"
    per_genome: dict[str, dict] = {}
    with open(rec_path, "w", encoding="utf-8") as rec_out, gzip.open(archaea_fna, "rb") as fh:
        rec_out.write("genome\tcontig_uid\tlength\n")
        it = v2.iter_fasta_stream(fh)
        for k in range(0, len(lst), batch):
            bdir = refs / f"batch_{k // batch:03d}"
            bdir.mkdir()
            files = []
            for row in lst[k:k + batch]:
                g, n_rec = row["genome"], int(row["records"])
                p = bdir / f"{g}.fna"
                bp = 0
                with open(p, "wb") as w:
                    for _ in range(n_rec):
                        r = next(it)
                        if not r.header.startswith(g):
                            raise Refused(f"reference record {r.accession} is not genome {g}'s")
                        w.write(b"\n".join(fasta.format_record(r.header, r.seq, r.line_width or 60)) + b"\n")
                        rec_out.write(f"{g}\t{ids.contig_uid(r.seq).hex()}\t{len(r.seq)}\n")
                        bp += len(r.seq)
                if bp != int(row["bp"]):
                    raise Refused(f"reference {g}: {bp} bp, the list says {row['bp']}")
                per_genome[g] = {"gtdb": row["gtdb"], "records": n_rec, "bp": bp, "taxonomy": row.get("taxonomy", "")}
                files.append(p)
            if skani:
                (bdir.parent / f"batch_{k // batch:03d}.txt").write_text("".join(f"{p}\n" for p in files))
                sk = refs / f"sk_{k // batch:03d}"
                _run([skani, "sketch", "--separate-sketches", "-t", threads, "-l",
                      bdir.parent / f"batch_{k // batch:03d}.txt", "-o", sk], log)
                for p in files:
                    s = sk / f"{p.name}.sketch"
                    if not s.is_file():
                        raise Refused(f"no sketch {s}")
                    sketches.append(str(s))
                shutil.rmtree(bdir)
        try:
            next(it)
            raise Refused("the reference FASTA has records the list does not name")
        except StopIteration:
            pass
    (work / "ref_sketches.txt").write_text("".join(s + "\n" for s in sketches))
    owners = c62._owners_from_records(rec_path, root, list(root_stages), owner_inputs, work=work / "_own",
                                      threads=2, memory_limit=memory_limit)
    later = c62._owners_from_records(rec_path, None, [], later_inputs, work=work / "_own_later", threads=2,
                                     memory_limit=memory_limit) if later_inputs else {}
    held: dict[str, Counter] = {g: Counter() for g in per_genome}
    with open(rec_path, encoding="utf-8") as fh:
        next(fh)
        for ln in fh:
            g, cu, n = ln.rstrip("\n").split("\t")
            b = bytes.fromhex(cu)
            s = owners.get(b) or later.get(b)
            if s:
                held[g][s] += int(n)
    rows = []
    for g in sorted(per_genome):
        d, h = per_genome[g], held[g]
        holder, hb = (h.most_common(1)[0] if h else ("", 0))
        rows.append((g, d["gtdb"], d["records"], d["bp"], holder, hb, f"{hb / d['bp']:.4f}" if d["bp"] else "0",
                     "yes" if d["bp"] and hb / d["bp"] >= 0.5 else "no", d["taxonomy"]))
    make_inputs._write_tsv(work / "ref_holders.tsv", ["genome", "gtdb", "records", "bp", "holder", "held_bp",
                                                      "held_fraction", "held", "taxonomy"], rows)
    return {"references": len(rows), "held": sum(1 for r in rows if r[7] == "yes"),
            "held_by": dict(sorted(Counter(r[4] for r in rows if r[7] == "yes").items())),
            "held_gtdb": sum(1 for r in rows if r[7] == "yes" and r[1] == "true"),
            "gtdb": sum(1 for r in rows if r[1] == "true"), "sketches": len(sketches), "log": log}


def dpn1_skani(work: Path, *, skani: str = "skani", threads: int = 8) -> dict:
    work = Path(work)
    log: list[dict] = []
    _run([skani, "dist", "-t", threads, "--ql", work / "queries.txt", "--rl", work / "ref_sketches.txt",
          "-o", work / "skani_dist.tsv"], log)
    _run([skani, "triangle", "-E", "-t", threads, "-l", work / "queries.txt", "-o", work / "skani_triangle.tsv"], log)
    return {"log": log}


def parse_skani_pairs(path: Path) -> list[tuple[str, str, float, float]]:
    """(query file, reference file, ANI, the larger aligned fraction) of a skani dist / triangle -E table, the file
    columns as skani printed them (the FASTA path each sketch was made from)."""
    import csv

    out = []
    with open(path, encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        need = ["Ref_file", "Query_file", "ANI", "Align_fraction_ref", "Align_fraction_query"]
        miss = [c for c in need if c not in (rd.fieldnames or [])]
        if miss:
            raise Refused(f"{path}: not skani output (no column {miss[0]})")
        for r in rd:
            out.append((r["Query_file"], r["Ref_file"], float(r["ANI"]),
                        max(float(r["Align_fraction_ref"]), float(r["Align_fraction_query"]))))
    return out


def _score(row: dict) -> tuple:
    """Rank for the within-stage pass: the closed genomes first, then CheckM2 completeness - 5 x contamination
    (unknown last), then length, then name."""
    try:
        q = float(row["completeness"]) - 5 * float(row["contamination"])
        known = 0
    except (TypeError, ValueError):
        q, known = 0.0, 1
    return (0 if row["category"] == "nanobdellota_closed" else 1, known, -q, -int(row["bp"]), row["name"],
            row["category"])


def dpn1_decide(manifest: list[dict], dist: list[tuple], tri: list[tuple], holders: dict[str, dict]) -> dict:
    """Which genomes stay (module notes above). Returns {"kept": [names in rank order], "removed": {name: row},
    "ref_best": {name: (ani, af, ref, holder, held)}}. A genome is keyed by `<category>:<name>` (a name of
    noncircular_51 can also be one of complete_208); skani's rows are matched by file path."""
    by_name: dict[str, dict] = {}
    key_of_file: dict[str, str] = {}
    for r in manifest:
        key = f"{r['category']}:{r['name']}"
        if key in by_name:
            raise Refused(f"genome {key} listed twice")
        by_name[key] = r
        key_of_file[str(r["file"])] = key
        key_of_file.setdefault(str(Path(r["file"]).resolve()), key)
    best: dict[str, tuple] = {}
    best_held: dict[str, tuple] = {}
    for qf, rf, ani, af in dist:
        k = key_of_file.get(qf)
        if k is None:
            raise Refused(f"skani query {qf} is not a DPN1 genome file")
        if ani < ANI_MIN or af < AF_MIN:
            continue
        ref = _gid(rf)
        h = holders.get(ref)
        if h is None:
            raise Refused(f"skani reference {rf} is not in ref_holders.tsv")
        cand = (ani, af, ref)
        if k not in best or cand > best[k][:3]:
            best[k] = (*cand, h.get("holder", ""), h.get("held") == "yes")
        if h.get("held") == "yes" and (k not in best_held or cand > best_held[k][:3]):
            best_held[k] = (*cand, h.get("holder", ""))
    removed: dict[str, dict] = {}
    for k, (ani, af, ref, holder) in sorted(best_held.items()):
        removed[k] = {"reason": "ani_to_held_reference", "match": ref, "ani": ani, "af": af, "holder": holder}
    adj: dict[str, dict[str, tuple]] = {}
    for qf, rf, ani, af in tri:
        a, b = key_of_file.get(qf), key_of_file.get(rf)
        if a is None or b is None:
            raise Refused(f"skani triangle row {qf} / {rf} is not a pair of DPN1 genome files")
        if ani < ANI_MIN or af < AF_MIN or a == b:
            continue
        adj.setdefault(a, {})[b] = (ani, af)
        adj.setdefault(b, {})[a] = (ani, af)
    order = sorted((k for k in by_name if k not in removed), key=lambda k: _score(by_name[k]))
    kept: list[str] = []
    gone: set[str] = set()
    for k in order:
        if k in gone:
            continue
        kept.append(k)
        for nb, (ani, af) in sorted(adj.get(k, {}).items()):
            if nb not in gone and nb not in removed and nb not in kept:
                gone.add(nb)
                removed[nb] = {"reason": "ani_within_stage", "match": by_name[k]["name"], "ani": ani, "af": af,
                               "holder": "DPN1"}
    return {"kept": kept, "removed": removed, "ref_best": best, "by_name": by_name}


DECISION_HEAD = ["genome", "category", "accession", "bp", "decision", "match", "ani", "af", "match_holder",
                 "best_reference", "best_reference_ani", "best_reference_af", "best_reference_holder",
                 "best_reference_held"]


def stage_dpn1(work: Path, mge_fna: Path, index_dir: Path | None, meta_out: Path, fasta_out: Path, *, root,
               root_stages, owner_inputs, later_inputs, inputs_home, retrieved_at: str = "2026-11-25T23:06:42Z",
               threads: int = 2, memory_limit: str = "4GB", screen=None, workers: int = 2,
               stage_work: Path | None = None) -> dict:
    import concurrent.futures
    import csv

    import r232_index

    work = Path(work)
    with open(work / "genomes.tsv", encoding="utf-8", newline="") as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    with open(work / "ref_holders.tsv", encoding="utf-8", newline="") as fh:
        holders = {r["genome"]: r for r in csv.DictReader(fh, delimiter="\t")}
    dist = parse_skani_pairs(work / "skani_dist.tsv")
    tri = parse_skani_pairs(work / "skani_triangle.tsv")
    dec = dpn1_decide(manifest, dist, tri, holders)
    by_name, kept, removed, best = dec["by_name"], dec["kept"], dec["removed"], dec["ref_best"]

    sf = StageFiles("DPN1", meta_out, fasta_out, DPN1_FASTA, stage_work)
    screen_fn = screen or _screen_default
    pool = concurrent.futures.ProcessPoolExecutor(max_workers=workers) if workers > 0 and screen is None else None
    try:
        genomes, trims, skipped, code_rows, pending = [], [], [], [], []
        classes: Counter = Counter()

        def settle(g, res):
            g["code"], g["reason"] = int(res["genetic_code"]), res.get("reason", "")
            code_rows.append((g["name"], g["code"], g["reason"], res.get("cd_11"), res.get("mean_orf_11"),
                              res.get("cd_25"), res.get("cd_4")))

        for k in sorted(kept, key=lambda x: (by_name[x]["category"], by_name[x]["name"])):
            row = by_name[k]
            closed = row["category"] == "nanobdellota_closed"
            contigs, circ, seqs = [], [], []
            with open(row["file"], "rb") as fh:
                for r in v2.iter_fasta_stream(fh):
                    src_cu = ids.contig_uid(r.seq)
                    out, kk, cls, is_c = dec65(r.seq, flagged=closed) if closed else (r.seq, 0, "", False)
                    if closed:
                        classes[cls] += 1
                    cu = ids.contig_uid(out) if kk else src_cu
                    written = sf.write("dpann", r.header, out, r.line_width, prefix=row["name"])
                    if kk:
                        trims.append((row["category"], row["name"], written, src_cu.hex(), cu.hex(), len(r.seq), kk,
                                      kk, len(out)))
                    if is_c:
                        circ.append(cu)
                    contigs.append((cu, len(out), len(r.seq), written, cls, src_cu, kk))
                    seqs.append(out)
            quality = "complete" if closed and len(circ) == len(contigs) else "contig"
            g = {"name": row["name"], "category": "dpann", "kind": "prokaryote", "quality": quality,
                 "circular": circ, "contigs": contigs, "code": 11, "reason": "", "source": row["category"]}
            genomes.append(g)
            org = row.get("taxonomy", "")
            if pool is not None:
                pending.append((g, pool.submit(_screen_job, (row["name"], seqs, org))))
                while len(pending) > 2 * workers:
                    gg, fut = pending.pop(0)
                    settle(gg, fut.result()[1])
            else:
                settle(g, screen_fn(seqs, org))
        while pending:
            gg, fut = pending.pop(0)
            settle(gg, fut.result()[1])

        seen: dict[bytes, str] = {}
        n_mge = 0
        with open(mge_fna, "rb") as fh:
            for r in v2.iter_fasta_stream(fh):
                n_mge += 1
                g = _single_record_genome(sf, "mge", Path(mge_fna).name, r.header, r.seq, r.line_width,
                                          kind="virus", flagged=False, seen=seen, trims=trims, skipped=skipped,
                                          dtr_classes=classes)
                if g is not None:
                    genomes.append(g)

        dec_rows = []
        for k in sorted(by_name, key=lambda x: (by_name[x]["category"], by_name[x]["name"])):
            row, rm, b = by_name[k], removed.get(k), best.get(k)
            dec_rows.append((row["name"], row["category"], row["accession"], row["bp"],
                             rm["reason"] if rm else "kept", rm["match"] if rm else "",
                             rm["ani"] if rm else "", rm["af"] if rm else "", rm["holder"] if rm else "",
                             b[2] if b else "", b[0] if b else "", b[1] if b else "", b[3] if b else "",
                             ("yes" if b[4] else "no") if b else ""))
        make_inputs._write_tsv(sf.tmp / "dpn1_decisions.tsv", DECISION_HEAD, dec_rows)

        staged = [g for g in genomes if g["category"] == "dpann"]
        accs = [by_name[k]["accession"] for k in kept if by_name[k]["accession"]]
        look = {x["accession"]: x for x in r232_index.lookup(index_dir, accs)} if index_dir and accs else {}
        status_rows, status_count = [], Counter()
        acc_of = {by_name[k]["name"]: by_name[k]["accession"] for k in kept}
        for g in sorted(staged, key=lambda x: x["name"]):
            a = acc_of.get(g["name"], "")
            x = look.get(a)
            st = r232_status_of(x) if a else "not_visible"
            status_count[st] += 1
            status_rows.append((g["name"], a, st, _yn(x["in_gtdb_r232"]) if x else "",
                                str(x["is_representative"]) if x and x["in_gtdb_r232"] else "",
                                (x["species"] or "") if x else "", _yn(x["in_globdb_r232"]) if x else "",
                                (x["globdb_r232_id"] or "") if x else "", len(g["contigs"]),
                                sum(c[1] for c in g["contigs"])))
        make_inputs._write_tsv(sf.tmp / "r232_status.tsv", STATUS_HEAD, status_rows)
        make_inputs._write_tsv(sf.tmp / "genetic_codes.tsv",
                               ["genome", "genetic_code", "reason", "cd_11", "mean_orf_11", "cd_25", "cd_4"],
                               sorted(tuple("" if v is None else v for v in r) for r in code_rows))
        genomes.sort(key=lambda g: (g["category"], g["name"]))
        rm_count = Counter(r["reason"] for r in removed.values())
        rm_bp = Counter()
        for k, r in removed.items():
            rm_bp[r["reason"]] += int(by_name[k]["bp"])
        rep = _finish(sf, genomes, trims, skipped, root=root, root_stages=root_stages, owner_inputs=owner_inputs,
                      later_inputs=later_inputs, inputs_home=inputs_home, retrieved_at=retrieved_at,
                      source=DPN1_SOURCE, source_label="tower_dpn1", threads=threads, memory_limit=memory_limit,
                      extra={"candidates": len(by_name),
                             "candidates_by_category": dict(sorted(Counter(r["category"]
                                                                           for r in by_name.values()).items())),
                             "removed": dict(sorted(rm_count.items())), "removed_bp": dict(sorted(rm_bp.items())),
                             "removed_to_held_by": dict(sorted(Counter(r["holder"] for r in removed.values()
                                                                       if r["reason"] == "ani_to_held_reference"
                                                                       ).items())),
                             "kept": len(kept),
                             "kept_by_category": dict(sorted(Counter(by_name[k]["category"] for k in kept).items())),
                             "kept_bp": sum(int(by_name[k]["bp"]) for k in kept),
                             "best_hit_unheld_kept": sum(1 for k in kept if k in best and not best[k][4]),
                             "mge_records": n_mge, "dtr_classes_closed_and_mge": dict(sorted(classes.items())),
                             "genetic_codes": dict(sorted(Counter(str(g["code"]) for g in staged).items())),
                             "r232_status": dict(sorted(status_count.items()))})
    except BaseException:
        sf.abort()
        raise
    finally:
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)
    return rep


def r232_status_of(x: dict | None) -> str:
    if x is None:
        return "not_visible"
    if x["in_gtdb_r232"] and x["in_globdb_r232"]:
        return "in_gtdb_and_globdb_r232"
    if x["in_gtdb_r232"]:
        return "in_gtdb_r232"
    return "in_globdb_r232" if x["in_globdb_r232"] else "not_in_r232"


def _yn(v) -> str:
    return "yes" if v else "no"


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
    p.add_argument("--reuse-fasta", action="store_true",
                   help="keep a FASTA already on /stor when the new one has the same bytes (refuse if not)")
    p = sub.add_parser("dpn1-extract")
    for n in ("nano_zip", "dpann_tgz", "dpann_meta", "work"):
        p.add_argument(n, type=Path)
    p = sub.add_parser("dpn1-refs")
    for n in ("archaea_fna", "archaea_list", "work"):
        p.add_argument(n, type=Path)
    p.add_argument("--root", type=Path, default=None)
    p.add_argument("--root-stages", default=",".join(ROOT_STAGES))
    p.add_argument("--chain-home", type=Path, default=None)
    p.add_argument("--later", action="append", default=[])
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--batch", type=int, default=1000)
    p.add_argument("--skani", default="skani")
    p.add_argument("--memory-limit", default="8GB")
    p = sub.add_parser("dpn1-skani")
    p.add_argument("work", type=Path)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--skani", default="skani")
    p = sub.add_parser("dpn1")
    p.add_argument("dpn1_work", type=Path)
    p.add_argument("mge_fna", type=Path)
    p.add_argument("index", type=Path)
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
        if a.cmd in ("dpn1-extract", "dpn1-skani"):
            rep = (dpn1_extract(a.nano_zip, a.dpann_tgz, a.dpann_meta, a.work) if a.cmd == "dpn1-extract"
                   else dpn1_skani(a.work, skani=a.skani, threads=a.threads))
            txt = json.dumps(rep, indent=1, sort_keys=True)
            (a.work / f"{a.cmd}.json").write_text(txt + "\n")
            print(txt)
            return 0
        owners = {s: a.chain_home / s for s, _ in a26.R26_11_4_CHAIN} if a.chain_home is not None else {}
        if a.cmd == "dpn1-refs":
            rep = dpn1_refs(a.archaea_fna, a.archaea_list, a.work,
                            root=a.root, root_stages=[s for s in a.root_stages.split(",") if s] if a.root else [],
                            owner_inputs=owners, later_inputs=_kv(a.later), skani=a.skani, threads=a.threads,
                            batch=a.batch, memory_limit=a.memory_limit)
            txt = json.dumps(rep, indent=1, sort_keys=True)
            (a.work / "dpn1-refs.json").write_text(txt + "\n")
            print(txt)
            return 0
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
            rep = stage_sym2(a.bins_tgz, a.viruses_tgz, a.tables, a.meta_out, a.fasta_out, workers=a.workers,
                             reuse_fasta=a.reuse_fasta, **kw)
        elif a.cmd == "dpn1":
            stage_name = "DPN1"
            stage_work = kw.pop("work")
            rep = stage_dpn1(a.dpn1_work, a.mge_fna, a.index, a.meta_out, a.fasta_out, workers=a.workers,
                             stage_work=stage_work, **kw)
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
