"""The `--genomes` and `--objects` files of `bio-box build-stage` for card C-10, from the import ladder's inputs.

    python make_inputs.py phage SRC_FASTA_GZ N OUT_DIR [--drop-repeats] [--min-contig-len L]
    python make_inputs.py prok LOCK_JSON SELECTION_TSV OUT_DIR [--source-release datasets-v2-YYYY-MM-DD]
        [--source-name S] [--license L] [--quality-column COL] [--root ROOT] [--drop-repeats] [--min-contig-len L]
    python make_inputs.py plasmid FASTA META_CSV OUT_DIR [--root ROOT] [--source-release R] [--retrieved-at T]
        [--min-contig-len L]
    python make_inputs.py plasmid-sources SPEC_JSON OUT_DIR [--root ROOT] [--min-contig-len L]

--source-release labels every prok object (default datasets-v2-2026-11-09, the ladder's fetch date).
--drop-repeats keeps the rule these inputs were made with before DEC-27: a contig uid already listed by an earlier
genome of the stage is not listed again, and a genome left with none is not listed. The inputs of stages S1, P1 to
P6 and H1 were made with it, and it gives them again byte for byte.

plasmid (card C-26, stage PL1): the PLSDB FASTA, read in place as one source object (source plsdb, release R, default
2024_05_31_v2; retrieved at T, default 2026-11-11T23:14:30Z, the end of its fetch), and PLSDB's nuccore.csv as
META_CSV. One genome per record, keyed by the header's first token = NUCCORE_ACC: kind plasmid, genetic code 11,
circular when NUCCORE_Topology is circular (linear when it is linear or empty; the empty ones are counted),
assembly_quality complete when NUCCORE_Completeness is complete and contig when it is empty. Any other value of
either field, or a record without a metadata row, is refused (exit 2) before anything is written. No record is
dropped for repeating a sequence (DEC-27): a record whose bytes repeat an earlier record's is that record's genome
(genome_uid hashes the contig set), listed once with the first record's metadata, and written to
OUT_DIR/repeated_records.tsv with the first record and whether the two agree on topology and quality; every record
is still ingested as an occurrence of the contig. With --root, the stages on disk in ROOT are read (their
contig_index files, read only): a record whose contig a stage stores as its default variant is still a genome of its
own (that stage owns the contig's calls, DEC-26), and is written to OUT_DIR/shared_records.tsv with the owning stage.
Writes genomes.json and objects.json and prints a JSON report.

plasmid-sources (card C-40, stage PL2): several plasmid sources, each with its own table. SPEC_JSON is
{"sources": [...]}, each source: name (the object's file stem), source_name, source_release, license, redistributable,
retrieved_at, fasta (a list of FASTA files, plain or .gz, read in that order), table (a TSV with a header line; lines
starting with '#' are skipped), id_column (the value equal to a record's header first token), topology_column and
topology (a map from the column's values to circular or linear), completeness_column and completeness (a map to
complete, fragment or unknown), optionally select ({"column", "values"}: a record whose row holds another value is not
taken) and without_row ("refuse", the default, or "skip": a record with no row is not taken). A value missing from a
map, a record without a row under "refuse", an accession twice in a source's FASTA files or twice in its table is
refused (exit 2) before anything is written. The records taken are written, in order, to OUT_DIR/<name>.fasta, the
source's one object in objects.json (header, sequence and line width kept, fasta.format_record). One genome per
distinct sequence over all sources, the first record taken: kind plasmid, genetic code 11, assembly_quality complete
for complete and contig for fragment or unknown, circular only when its topology is circular and it is not a fragment
(a fragment its table calls circular is counted, never circular). No record is dropped for repeating another's bytes
(DEC-27): it is an occurrence of the first record's genome. OUT_DIR/records.tsv lists every record taken (source,
accession, contig_uid, length, the table's topology and completeness values, what was made of them, the first record
of its sequence when it is not the first, and with --root the stage on disk that owns its contig, DEC-26);
OUT_DIR/skipped_records.tsv every record not taken, with the reason (not_selected, no_row). Writes genomes.json and
objects.json and prints a JSON report with the counts per source.

phage: the first N records of the INPHARED FASTA, copied byte for byte (header and sequence lines as in the
source) to OUT_DIR/phage_first<N>.fa; every record its own genome (kind virus, genetic code 11,
linear: FASTA carries no topology), as the ladder's stage S1 was described (step.py, do_genome_input).
prok: one genome per accession of the ladder's lock.json (sorted), its contigs in record order, kind
prokaryote, linear, genetic code from the selection's translation_table (11 when absent), assembly_quality complete by
default or taken from SELECTION_TSV's `--quality-column COL` (one of complete, chromosome, scaffold, contig, mag_hq,
mag_mq; card C-41), as the ladder's stage P1 and MAG stages G1..G41 are described. The source objects are the lock's
files, read in place (default source_name ncbi_genomes and license "INSDC public data", overridden by `--source-name`
and `--license`). Contig uids are ids.contig_uid of each record's sequence. Every record is listed in its genome
(DEC-27, card C-21a): a contig in several genomes is listed in each. A contig repeated inside one genome is listed once
in it, and a genome whose contigs are those of an earlier genome is that genome (genome_uid hashes the contig set) and
is listed once; both are counted in the report, with the contigs listed in more than one genome. With `--root ROOT`,
contigs stored as default variants by stages on disk in ROOT are reported in OUT_DIR/shared_records.tsv (DEC-26) while
kept in their genomes. Any unmapped quality, missing column, duplicate selection row or missing lock row is refused
(exit 2) before anything is written. Writes OUT_DIR/genomes.json and OUT_DIR/objects.json and prints a JSON report.
"""

import csv
import gzip
import json
import sys
from pathlib import Path

from biobox import fasta, ids

USAGE = __doc__.split("\n\n")[1]
HELP = "\n\n".join(__doc__.split("\n\n")[1:])

PROK_QUALITIES = ("complete", "chromosome", "scaffold", "contig", "mag_hq", "mag_mq")

# card C-26: PLSDB's metadata values and what the stage's inputs make of them
PLSDB_TOPOLOGY = {"circular": True, "linear": False, "": False}
PLSDB_QUALITY = {"complete": "complete", "": "contig"}
PLSDB_RELEASE = "2024_05_31_v2"
PLSDB_RETRIEVED = "2026-11-11T23:14:30Z"  # sources.yaml / the S2 inputs note: download_fasta ended 23:14:30 UTC


class Refused(Exception):
    """Inputs the script will not turn into a stage (exit 2, nothing written)."""


def first_records(src: Path, n: int, out: Path) -> int:
    """Copy the lines of the first `n` records of `src` to `out`; return how many were copied."""
    k = 0
    tmp = out.with_name(f".{out.name}.tmp")
    with gzip.open(src, "rb") as fh, open(tmp, "wb") as w:
        for line in fh:
            if line.startswith(b">"):
                k += 1
                if k > n:
                    k = n
                    break
            w.write(line)
    tmp.replace(out)
    return k


def genome(uids: list[bytes], seen: set, code: int, kind: str, quality: str) -> tuple[dict | None, int]:
    """The rule before DEC-27 (`--drop-repeats`; `s2_h1_inputs.py` uses it for stage H1): the genome without the
    uids in `seen`, which it then adds, and how many uids it dropped."""
    kept = [u for u in dict.fromkeys(uids) if u not in seen]
    seen.update(kept)
    dropped = len(uids) - len(kept)
    if not kept:
        return None, dropped
    return {"contig_uids": [u.hex() for u in kept], "genetic_code": code, "kind": kind, "circular": [],
            "assembly_quality": quality}, dropped


class Listing:
    """The rule of DEC-27: every record is listed in its genome. `add` returns the genome's entry, or None for a
    genome already listed (the same contig set), and counts what it did not list."""

    def __init__(self):
        self.genomes: set[bytes] = set()
        self.memberships: dict[bytes, int] = {}  # contig uid -> genomes listing it
        self.in_genome = self.repeated = 0

    def add(self, uids: list[bytes], code: int, kind: str, quality: str) -> dict | None:
        kept = list(dict.fromkeys(uids))
        self.in_genome += len(uids) - len(kept)
        guid = ids.genome_uid(kept)
        if guid in self.genomes:
            self.repeated += 1
            return None
        self.genomes.add(guid)
        for u in kept:
            self.memberships[u] = self.memberships.get(u, 0) + 1
        return {"contig_uids": [u.hex() for u in kept], "genetic_code": code, "kind": kind, "circular": [],
                "assembly_quality": quality}

    def report(self) -> dict:
        return {"memberships": sum(self.memberships.values()),
                "contigs_in_several_genomes": sum(1 for n in self.memberships.values() if n > 1),
                "records_repeated_in_genome": self.in_genome, "genomes_repeated": self.repeated}


def plasmid_meta(path: Path) -> dict[str, dict]:
    """NUCCORE_ACC -> its topology, completeness and STATUS from PLSDB's nuccore.csv (a quoted CSV)."""
    with open(path, newline="") as fh:
        rows = csv.DictReader(fh)
        need = {"NUCCORE_ACC", "NUCCORE_Topology", "NUCCORE_Completeness"}
        if not need <= set(rows.fieldnames or ()):
            raise Refused(f"{path}: no column {sorted(need - set(rows.fieldnames or ()))[0]}")
        meta: dict[str, dict] = {}
        for r in rows:
            acc = r["NUCCORE_ACC"]
            if acc in meta:
                raise Refused(f"{path}: accession {acc} has two rows")
            meta[acc] = {"topology": r["NUCCORE_Topology"], "completeness": r["NUCCORE_Completeness"],
                         "status": r.get("STATUS") or ""}
    return meta


def stage_owners(root: Path) -> dict[bytes, str]:
    """Each contig a stage on disk in `root` stores as its default variant, and that stage (DEC-26: the stage whose
    contig_index holds the default variant owns the contig; `store.stored_by` is the same rule one key at a time)."""
    import pyarrow.parquet as pq

    from biobox import io
    from biobox.l0seq import store

    owners: dict[bytes, str] = {}
    for s in store.stages_on_disk(root):
        for f in io.part_files(Path(root) / store.INDEX_DIR / f"stage={s}"):
            t = pq.read_table(f, columns=["contig_uid", "is_default_variant"])
            for cu, d in zip(t.column("contig_uid").to_pylist(), t.column("is_default_variant").to_pylist()):
                if not d:
                    continue
                if owners.get(cu, s) != s:
                    raise Refused(f"contig {cu.hex()} is a default variant of stages {owners[cu]} and {s}")
                owners[cu] = s
    return owners


def _write_tsv(path: Path, head: list[str], rows: list[tuple]) -> None:
    path.write_text("".join("\t".join(map(str, r)) + "\n" for r in [tuple(head), *rows]))


def plasmid(fa: Path, meta_csv: Path, out: Path, root: Path | None, release: str, retrieved: str,
            *, min_len: int | None = None) -> dict:
    """The plasmid mode (module docstring): reads everything, refuses before writing, then writes the files."""
    meta = plasmid_meta(meta_csv)
    owners = stage_owners(root) if root is not None else None
    listing = Listing()
    genomes: list[dict] = []
    first: dict[bytes, tuple] = {}  # contig uid -> (ordinal, accession, circular, quality, length) of its first record
    seen_acc: set[str] = set()
    topo = {"circular": 0, "linear": 0, "empty": 0}
    qual = {q: 0 for q in PLSDB_QUALITY.values()}
    repeated, shared = [], []
    suppressed = record_bases = records_below = bases_below = 0
    for k, r in enumerate(fasta.iter_records(fa)):
        acc = r.header.split()[0] if r.header.split() else ""
        m = meta.get(acc)
        if m is None:
            raise Refused(f"record {k} of {fa}, {acc!r}: no row in {meta_csv}")
        if acc in seen_acc:
            raise Refused(f"record {k} of {fa}: accession {acc} is a second record")
        seen_acc.add(acc)
        if m["topology"] not in PLSDB_TOPOLOGY:
            raise Refused(f"{acc}: NUCCORE_Topology {m['topology']!r} is not one of {sorted(PLSDB_TOPOLOGY)}")
        if m["completeness"] not in PLSDB_QUALITY:
            raise Refused(f"{acc}: NUCCORE_Completeness {m['completeness']!r} is not one of {sorted(PLSDB_QUALITY)}")
        if min_len is not None and len(r.seq) < min_len:
            records_below += 1
            bases_below += len(r.seq)
            continue
        circ, q = PLSDB_TOPOLOGY[m["topology"]], PLSDB_QUALITY[m["completeness"]]
        topo[m["topology"] or "empty"] += 1
        qual[q] += 1
        suppressed += m["status"] == "suppressed"
        record_bases += len(r.seq)
        cu = ids.contig_uid(r.seq)
        if owners is not None and cu in owners:
            shared.append((k, acc, cu.hex(), owners[cu]))
        if cu in first:
            f = first[cu]
            repeated.append((k, acc, cu.hex(), f[0], f[1], str(circ == f[2]).lower(), str(q == f[3]).lower()))
            listing.add([cu], 11, "plasmid", q)  # counts the genome as repeated; returns None
            continue
        g = listing.add([cu], 11, "plasmid", q)
        g["circular"] = [cu.hex()] if circ else []
        genomes.append(g)
        first[cu] = (k, acc, circ, q, len(r.seq))
    rep = {"records": len(seen_acc), "genomes": len(genomes), "contigs": len(first), "topology": topo,
           "quality": qual, "genomes_circular": sum(1 for g in genomes if g["circular"]),
           "repeated_records": len(repeated),
           "repeats_other_topology": sum(1 for r in repeated if r[5] == "false"),
           "repeats_other_quality": sum(1 for r in repeated if r[6] == "false"),
           "suppressed": suppressed, "record_bases": record_bases,
           "genome_bases": sum(f[4] for f in first.values()), "meta_rows_not_in_fasta": len(set(meta) - seen_acc),
           **listing.report()}
    if owners is not None:
        by_stage = lambda rows: {s: sum(1 for r in rows if r[3] == s) for s in sorted({r[3] for r in rows})}  # noqa: E731
        rep.update({"shared_records": len(shared), "owner_records": by_stage(shared),
                    "owner_contigs": by_stage(list({r[2]: r for r in shared}.values())),
                    "new_contigs": sum(1 for cu in first if cu not in owners),
                    "new_bases": sum(f[4] for cu, f in first.items() if cu not in owners)})
    objects = [{"path": str(fa), "source_name": "plsdb", "source_release": release,
                "license": "INSDC records; PLSDB states no licence (cite doi:10.1093/nar/gkab1111)",
                "redistributable": True, "retrieved_at": retrieved}]
    out.mkdir(parents=True, exist_ok=True)
    _write_tsv(out / "repeated_records.tsv", ["ordinal", "accession", "contig_uid", "first_ordinal",
                                              "first_accession", "same_topology", "same_quality"], repeated)
    if owners is not None:
        _write_tsv(out / "shared_records.tsv", ["ordinal", "accession", "contig_uid", "owner_stage"], shared)
    (out / "objects.json").write_text(json.dumps(objects, indent=1))
    (out / "genomes.json").write_text(json.dumps(genomes))
    rep.update({"objects": len(objects), "codes": {"11": len(genomes)} if genomes else {}})
    if min_len is not None:
        rep.update({"min_contig_len": min_len, "records_below_min_len": records_below,
                    "bases_below_min_len": bases_below})
    return rep


SOURCE_KEYS = {"name", "source_name", "source_release", "license", "redistributable", "retrieved_at", "fasta", "table",
               "id_column", "topology_column", "topology", "completeness_column", "completeness"}
TOPOLOGIES = ("circular", "linear")
COMPLETENESS = {"complete": "complete", "fragment": "contig", "unknown": "contig"}  # -> assembly_quality


def _source_table(src: dict) -> dict[str, tuple[str, str, bool]]:
    """accession -> (topology value, completeness value, selected) from a source's table (module docstring)."""
    sel = src.get("select")
    path = Path(src["table"])
    rows: dict[str, tuple[str, str, bool]] = {}
    with open(path, newline="") as fh:
        lines = (ln for ln in fh if not ln.startswith("#"))
        reader = csv.DictReader(lines, delimiter="\t", quoting=csv.QUOTE_NONE)
        cols = [src["id_column"], src["topology_column"], src["completeness_column"]] + ([sel["column"]] if sel else [])
        missing = [c for c in cols if c not in (reader.fieldnames or ())]
        if missing:
            raise Refused(f"{path}: no column {missing[0]!r}")
        for r in reader:
            acc = r[src["id_column"]]
            if acc in rows:
                raise Refused(f"{path}: accession {acc} has two rows")
            rows[acc] = (r[src["topology_column"]], r[src["completeness_column"]],
                         sel is None or r[sel["column"]] in sel["values"])
    return rows


def _check_spec(spec: dict) -> list[dict]:
    srcs = spec.get("sources") if isinstance(spec, dict) else None
    if not srcs:
        raise Refused("the spec lists no sources")
    names = set()
    for s in srcs:
        missing = sorted(SOURCE_KEYS - set(s))
        if missing:
            raise Refused(f"source {s.get('name')!r}: no {missing[0]}")
        if s["name"] in names or "/" in s["name"]:
            raise Refused(f"source name {s['name']!r} is repeated or not a file stem")
        names.add(s["name"])
        for v in s["topology"].values():
            if v not in TOPOLOGIES:
                raise Refused(f"source {s['name']}: topology maps to {v!r}, not one of {TOPOLOGIES}")
        for v in s["completeness"].values():
            if v not in COMPLETENESS:
                raise Refused(f"source {s['name']}: completeness maps to {v!r}, not one of {sorted(COMPLETENESS)}")
        if s.get("without_row", "refuse") not in ("refuse", "skip"):
            raise Refused(f"source {s['name']}: without_row is {s['without_row']!r}")
    return srcs


def plasmid_sources(spec_path: Path, out: Path, root: Path | None, *, min_len: int | None = None) -> dict:
    """The plasmid-sources mode (module docstring): reads and checks everything into a sibling scratch directory,
    which becomes OUT_DIR only when nothing was refused."""
    import shutil

    srcs = _check_spec(json.loads(Path(spec_path).read_text()))
    if out.exists():
        raise Refused(f"{out} exists")
    owners = stage_owners(root) if root is not None else None
    tmp = out.with_name(f".{out.name}.tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    try:
        rep = _plasmid_sources(srcs, out, tmp, owners, min_len=min_len)
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    tmp.rename(out)
    return rep


def _plasmid_sources(srcs: list[dict], out: Path, tmp: Path, owners: dict[bytes, str] | None,
                     *, min_len: int | None = None) -> dict:
    listing = Listing()
    genomes: list[dict] = []
    first: dict[bytes, tuple] = {}  # contig uid -> (source, accession, length) of its first record taken
    objects, per = [], {}
    n_taken = 0
    owner_records: dict[str, int] = {}
    rec_w = open(tmp / "records.tsv", "w")
    skip_w = open(tmp / "skipped_records.tsv", "w")
    rec_w.write("\t".join(["ordinal", "source", "accession", "contig_uid", "length", "topology_value",
                           "completeness_value", "circular", "assembly_quality", "fragment", "first_source",
                           "first_accession", "owner_stage"]) + "\n")
    skip_w.write("source\taccession\tcontig_uid\tlength\treason\n")
    for src in srcs:
        name = src["name"]
        rows = _source_table(src)
        st = {"fasta_records": 0, "taken": 0, "skipped": {}, "bases_taken": 0, "topology": {}, "completeness": {},
              "genomes": 0, "genome_bases": 0, "genomes_circular": 0, "genomes_complete": 0,
              "circular_but_fragment": 0, "repeats_within_source": 0, "repeats_of_an_earlier_source": 0,
              "owner_records": {}}
        seen: set[str] = set()
        fa_out = tmp / f"{name}.fasta"
        with open(fa_out, "wb") as w:
            for f in src["fasta"]:
                for r in fasta.iter_records(Path(f)):
                    acc = r.accession
                    if acc in seen:
                        raise Refused(f"source {name}: accession {acc} is a second record ({f})")
                    seen.add(acc)
                    st["fasta_records"] += 1
                    row = rows.get(acc)
                    reason = None
                    if row is None:
                        if src.get("without_row", "refuse") == "refuse":
                            raise Refused(f"source {name}: record {acc} of {f} has no row in {src['table']}")
                        reason = "no_row"
                    elif not row[2]:
                        reason = "not_selected"
                    cu = ids.contig_uid(r.seq)
                    if reason:
                        st["skipped"][reason] = st["skipped"].get(reason, 0) + 1
                        skip_w.write(f"{name}\t{acc}\t{cu.hex()}\t{len(r.seq)}\t{reason}\n")
                        continue
                    tv, cv, _ = row
                    if tv not in src["topology"]:
                        raise Refused(f"source {name}: {acc}: {src['topology_column']} {tv!r} is not one of "
                                      f"{sorted(src['topology'])}")
                    if cv not in src["completeness"]:
                        raise Refused(f"source {name}: {acc}: {src['completeness_column']} {cv!r} is not one of "
                                      f"{sorted(src['completeness'])}")
                    if min_len is not None and len(r.seq) < min_len:
                        reason = f"below:{min_len}"
                        st["skipped"][reason] = st["skipped"].get(reason, 0) + 1
                        skip_w.write(f"{name}\t{acc}\t{cu.hex()}\t{len(r.seq)}\t{reason}\n")
                        continue
                    comp = src["completeness"][cv]
                    frag = comp == "fragment"
                    topo_circ = src["topology"][tv] == "circular"
                    circ = topo_circ and not frag
                    q = COMPLETENESS[comp]
                    st["topology"][tv] = st["topology"].get(tv, 0) + 1
                    st["completeness"][cv] = st["completeness"].get(cv, 0) + 1
                    st["circular_but_fragment"] += topo_circ and frag
                    st["taken"] += 1
                    st["bases_taken"] += len(r.seq)
                    w.write(b"\n".join(fasta.format_record(r.header, r.seq, r.line_width)) + b"\n")
                    owner = owners.get(cu, "") if owners is not None else ""
                    if owner:
                        st["owner_records"][owner] = st["owner_records"].get(owner, 0) + 1
                        owner_records[owner] = owner_records.get(owner, 0) + 1
                    f0 = first.get(cu)
                    if f0 is not None:
                        st["repeats_within_source" if f0[0] == name else "repeats_of_an_earlier_source"] += 1
                        listing.add([cu], 11, "plasmid", q)  # counts the genome as repeated; returns None
                    else:
                        g = listing.add([cu], 11, "plasmid", q)
                        g["circular"] = [cu.hex()] if circ else []
                        genomes.append(g)
                        first[cu] = (name, acc, len(r.seq))
                        st["genomes"] += 1
                        st["genome_bases"] += len(r.seq)
                        st["genomes_circular"] += circ
                        st["genomes_complete"] += q == "complete"
                    rec_w.write("\t".join(map(str, [
                        n_taken, name, acc, cu.hex(), len(r.seq), tv, cv, str(circ).lower(), q, str(frag).lower(),
                        "" if f0 is None else f0[0], "" if f0 is None else f0[1], owner])) + "\n")
                    n_taken += 1
        if not st["taken"]:
            raise Refused(f"source {name}: no record taken")
        st["selected_rows_without_record"] = sum(1 for a, v in rows.items() if v[2] and a not in seen)
        per[name] = st
        objects.append({"path": str(out / f"{name}.fasta"), "source_name": src["source_name"],
                        "source_release": src["source_release"], "license": src["license"],
                        "redistributable": bool(src["redistributable"]), "retrieved_at": src["retrieved_at"]})
    rec_w.close()
    skip_w.close()
    (tmp / "objects.json").write_text(json.dumps(objects, indent=1))
    (tmp / "genomes.json").write_text(json.dumps(genomes))
    rep = {"records_taken": n_taken, "genomes": len(genomes), "contigs": len(first),
           "genomes_circular": sum(1 for g in genomes if g["circular"]),
           "genome_bases": sum(f[2] for f in first.values()), "objects": len(objects),
           "codes": {"11": len(genomes)}, "sources": per, **listing.report()}
    if owners is not None:
        rep.update({"owner_records": dict(sorted(owner_records.items())),
                    "new_contigs": sum(1 for cu in first if cu not in owners),
                    "new_bases": sum(f[2] for cu, f in first.items() if cu not in owners)})
    if min_len is not None:
        rep["min_contig_len"] = min_len
    return rep


def _option(argv: list[str], name: str) -> tuple[list[str], str | None]:
    """`argv` without `name VALUE`, and VALUE (None when absent); Refused when VALUE is missing."""
    if name not in argv:
        return argv, None
    k = argv.index(name)
    if k + 1 >= len(argv):
        raise Refused(f"{name} needs a value")
    return argv[:k] + argv[k + 2 :], argv[k + 1]


def _prok_selection(path: Path, qual_col: str | None) -> dict[str, dict]:
    need = ["accession", "translation_table"] + ([qual_col] if qual_col is not None else [])
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = set(reader.fieldnames or ())
        missing = [c for c in need if c not in fields]
        if missing:
            raise Refused(f"{path}: no column {missing[0]!r}")
        meta: dict[str, dict] = {}
        for r in reader:
            acc = r["accession"]
            if acc in meta:
                raise Refused(f"{path}: accession {acc} has two rows")
            if qual_col is not None and r[qual_col] not in PROK_QUALITIES:
                raise Refused(f"{acc}: {qual_col} {r[qual_col]!r} is not one of {PROK_QUALITIES}")
            meta[acc] = r
    return meta


def _parse_min_len(raw: str | None) -> int | None:
    if raw is None:
        return None
    try:
        v = int(raw)
    except ValueError:
        raise Refused(f"--min-contig-len {raw!r}: must be an integer >= 1") from None
    if v < 1:
        raise Refused(f"--min-contig-len {v}: must be >= 1")
    return v


def main(argv: list[str]) -> int:
    if any(a in ("-h", "--help") for a in argv):
        print(HELP)
        return 0
    try:
        argv, release = _option(argv, "--source-release")
        argv, root = _option(argv, "--root")
        argv, retrieved = _option(argv, "--retrieved-at")
        argv, src_name = _option(argv, "--source-name")
        argv, lic = _option(argv, "--license")
        argv, qual_col = _option(argv, "--quality-column")
        argv, min_len_s = _option(argv, "--min-contig-len")
        min_len = _parse_min_len(min_len_s)
    except Refused as e:
        print(f"make_inputs.py: refused: {e}\n{USAGE}", file=sys.stderr)
        return 2
    if argv[:1] == ["plasmid-sources"]:
        if len(argv) != 3 or any(x is not None for x in (release, retrieved, src_name, lic, qual_col)):
            print(USAGE, file=sys.stderr)
            return 2
        try:
            rep = plasmid_sources(Path(argv[1]), Path(argv[2]), None if root is None else Path(root), min_len=min_len)
        except Refused as e:
            print(f"make_inputs.py plasmid-sources: refused: {e}", file=sys.stderr)
            return 2
        print(json.dumps(rep))
        return 0
    if argv[:1] == ["plasmid"]:
        if len(argv) != 4 or "--drop-repeats" in argv or any(x is not None for x in (src_name, lic, qual_col)):
            print(USAGE, file=sys.stderr)
            return 2
        try:
            rep = plasmid(Path(argv[1]), Path(argv[2]), Path(argv[3]), None if root is None else Path(root),
                          release or PLSDB_RELEASE, retrieved or PLSDB_RETRIEVED, min_len=min_len)
        except Refused as e:
            print(f"make_inputs.py plasmid: refused: {e}", file=sys.stderr)
            return 2
        print(json.dumps(rep))
        return 0
    release = release or "datasets-v2-2026-11-09"
    if retrieved is not None:
        print(USAGE, file=sys.stderr)
        return 2
    drop = "--drop-repeats" in argv
    argv = [a for a in argv if a != "--drop-repeats"]
    if len(argv) != 4 or argv[0] not in ("phage", "prok"):
        print(USAGE, file=sys.stderr)
        return 2
    if (argv[0] == "phage" or drop) and any(x is not None for x in (root, src_name, lic, qual_col)):
        print(USAGE, file=sys.stderr)
        return 2
    out = Path(argv[3])
    seen: set = set()
    listing = Listing()
    genomes, dropped, rep = [], 0, {}
    records_below = bases_below = genomes_below = 0

    def add(uids: list[bytes], code: int, kind: str, quality: str = "complete") -> None:
        nonlocal dropped
        if drop:
            g, d = genome(uids, seen, code, kind, quality)
            dropped += d
        else:
            g = listing.add(uids, code, kind, quality)
        genomes.extend([g] if g else [])

    if argv[0] == "phage":
        out.mkdir(parents=True, exist_ok=True)
        n = int(argv[2])
        fa = out / f"phage_first{n}.fa"
        rep["records_copied"] = first_records(Path(argv[1]), n, fa)
        for r in fasta.iter_records(fa):
            if min_len is not None and len(r.seq) < min_len:
                records_below += 1
                bases_below += len(r.seq)
                continue
            add([ids.contig_uid(r.seq)], 11, "virus")
        objects = [{"path": str(fa), "source_name": "inphared", "source_release": "14Apr2025",
                    "license": "AGPL-3.0 (repo); data from INSDC"}]
        shared = []
        owners = None
        contig_len: dict[bytes, int] = {}
    else:
        try:
            lock = json.loads(Path(argv[1]).read_text())["genomes"]
            meta = _prok_selection(Path(argv[2]), qual_col)
            for acc in sorted(lock):
                if acc not in meta:
                    raise Refused(f"lock accession {acc} has no row in {argv[2]}")
            owners = stage_owners(Path(root)) if root is not None else None
            objects, shared = [], []
            contig_len = {}
            for acc in sorted(lock):
                v = lock[acc]
                row = meta[acc]
                tt = row["translation_table"]
                code = int(tt) if tt not in ("", "None") else 11
                q = row[qual_col] if qual_col is not None else "complete"
                uids: list[bytes] = []
                saw_any = False
                for ord_k, r in enumerate(fasta.iter_records(Path(v["file"]))):
                    saw_any = True
                    if min_len is not None and len(r.seq) < min_len:
                        records_below += 1
                        bases_below += len(r.seq)
                        continue
                    cu = ids.contig_uid(r.seq)
                    uids.append(cu)
                    contig_len.setdefault(cu, len(r.seq))
                    if owners is not None and cu in owners:
                        shared.append((acc, ord_k, r.accession, cu.hex(), len(r.seq), owners[cu]))
                if not saw_any:
                    raise Refused(f"{acc}: {v['file']} has no FASTA records")
                if uids:
                    add(uids, code, "prokaryote", q)
                else:
                    genomes_below += 1
                objects.append({"path": v["file"], "source_name": src_name or "ncbi_genomes",
                                "source_release": release, "license": lic or "INSDC public data",
                                "redistributable": True, "retrieved_at": v["retrieved_utc"]})
        except Refused as e:
            print(f"make_inputs.py prok: refused: {e}", file=sys.stderr)
            return 2
        out.mkdir(parents=True, exist_ok=True)
    (out / "genomes.json").write_text(json.dumps(genomes))
    (out / "objects.json").write_text(json.dumps(objects, indent=1))
    if owners is not None:
        _write_tsv(out / "shared_records.tsv",
                   ["accession", "record_ordinal", "record_accession", "contig_uid", "length", "owner_stage"], shared)
    distinct_uids = {bytes.fromhex(u) for g in genomes for u in g["contig_uids"]}
    contigs = len(distinct_uids)
    rep.update({"genomes": len(genomes), "contigs": contigs,
                **({"contigs_not_listed_again": dropped} if drop else listing.report()), "objects": len(objects),
                "codes": {str(c): sum(1 for g in genomes if g["genetic_code"] == c)
                          for c in sorted({g["genetic_code"] for g in genomes})}})
    if min_len is not None:
        rep.update({"min_contig_len": min_len, "records_below_min_len": records_below,
                    "bases_below_min_len": bases_below,
                    **({"genomes_below_min_len": genomes_below} if argv[0] == "prok" else {})})
    if qual_col is not None:
        rep["quality"] = {q: sum(1 for g in genomes if g["assembly_quality"] == q)
                          for q in sorted({g["assembly_quality"] for g in genomes})}
    if owners is not None:
        by_stage = lambda rows: {s: sum(1 for r in rows if r[5] == s) for s in sorted({r[5] for r in rows})}  # noqa: E731
        rep.update({
            "genome_bases": sum(contig_len[cu] for cu in distinct_uids),
            "shared_records": len(shared),
            "owner_records": by_stage(shared),
            "owner_contigs": by_stage(list({r[3]: r for r in shared}.values())),
            "new_contigs": sum(1 for cu in distinct_uids if cu not in owners),
            "new_bases": sum(contig_len[cu] for cu in distinct_uids if cu not in owners),
        })
    print(json.dumps(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
