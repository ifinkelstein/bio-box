"""The `--genomes` and `--objects` files of `bio-box build-stage` for card C-10, from the import ladder's inputs.

    python make_inputs.py phage SRC_FASTA_GZ N OUT_DIR [--drop-repeats]
    python make_inputs.py prok LOCK_JSON SELECTION_TSV OUT_DIR [--source-release datasets-v2-YYYY-MM-DD] [--drop-repeats]
    python make_inputs.py plasmid FASTA META_CSV OUT_DIR [--root ROOT] [--source-release R] [--retrieved-at T]

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

phage: the first N records of the INPHARED FASTA, copied byte for byte (header and sequence lines as in the
source) to OUT_DIR/phage_first<N>.fa; every record its own genome (kind virus, genetic code 11,
linear: FASTA carries no topology), as the ladder's stage S1 was described (step.py, do_genome_input).
prok: one genome per accession of the ladder's lock.json (sorted), its contigs in record order, kind
prokaryote, assembly_quality complete, linear, genetic code from the selection's translation_table (11
when absent), as the ladder's stage P1 was described. The source objects are the ladder's files, read
in place. Contig uids are ids.contig_uid of each record's sequence. Every record is listed in its genome (DEC-27,
card C-21a): a contig in several genomes is listed in each. A contig repeated inside one genome is listed once in
it, and a genome whose contigs are those of an earlier genome is that genome (genome_uid hashes the contig set) and
is listed once; both are counted in the report, with the contigs listed in more than one genome. Writes
OUT_DIR/genomes.json and OUT_DIR/objects.json and prints a JSON report.
"""

import csv
import gzip
import json
import sys
from pathlib import Path

from biobox import fasta, ids

USAGE = __doc__.split("\n\n")[1]
HELP = "\n\n".join(__doc__.split("\n\n")[1:4])

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


def plasmid(fa: Path, meta_csv: Path, out: Path, root: Path | None, release: str, retrieved: str) -> dict:
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
    suppressed = record_bases = 0
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
    return rep


def _option(argv: list[str], name: str) -> tuple[list[str], str | None]:
    """`argv` without `name VALUE`, and VALUE (None when absent); Refused when VALUE is missing."""
    if name not in argv:
        return argv, None
    k = argv.index(name)
    if k + 1 >= len(argv):
        raise Refused(f"{name} needs a value")
    return argv[:k] + argv[k + 2 :], argv[k + 1]


def main(argv: list[str]) -> int:
    if any(a in ("-h", "--help") for a in argv):
        print(HELP)
        return 0
    try:
        argv, release = _option(argv, "--source-release")
        argv, root = _option(argv, "--root")
        argv, retrieved = _option(argv, "--retrieved-at")
    except Refused:
        print(USAGE, file=sys.stderr)
        return 2
    if argv[:1] == ["plasmid"]:
        if len(argv) != 4 or "--drop-repeats" in argv:  # DEC-27: no record is dropped
            print(USAGE, file=sys.stderr)
            return 2
        try:
            rep = plasmid(Path(argv[1]), Path(argv[2]), Path(argv[3]), None if root is None else Path(root),
                          release or PLSDB_RELEASE, retrieved or PLSDB_RETRIEVED)
        except Refused as e:
            print(f"make_inputs.py plasmid: refused: {e}", file=sys.stderr)
            return 2
        print(json.dumps(rep))
        return 0
    release = release or "datasets-v2-2026-11-09"
    if root is not None or retrieved is not None:
        print(USAGE, file=sys.stderr)
        return 2
    drop = "--drop-repeats" in argv
    argv = [a for a in argv if a != "--drop-repeats"]
    if len(argv) != 4 or argv[0] not in ("phage", "prok"):
        print(USAGE, file=sys.stderr)
        return 2
    out = Path(argv[3])
    out.mkdir(parents=True, exist_ok=True)
    seen: set = set()
    listing = Listing()
    genomes, dropped, rep = [], 0, {}

    def add(uids: list[bytes], code: int, kind: str) -> None:
        nonlocal dropped
        if drop:
            g, d = genome(uids, seen, code, kind, "complete")
            dropped += d
        else:
            g = listing.add(uids, code, kind, "complete")
        genomes.extend([g] if g else [])

    if argv[0] == "phage":
        n = int(argv[2])
        fa = out / f"phage_first{n}.fa"
        rep["records_copied"] = first_records(Path(argv[1]), n, fa)
        for r in fasta.iter_records(fa):
            add([ids.contig_uid(r.seq)], 11, "virus")
        objects = [{"path": str(fa), "source_name": "inphared", "source_release": "14Apr2025",
                    "license": "AGPL-3.0 (repo); data from INSDC"}]
    else:
        lock = json.loads(Path(argv[1]).read_text())["genomes"]
        with open(argv[2]) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            meta = {r["accession"]: r for r in (dict(zip(head, ln.rstrip("\n").split("\t"))) for ln in fh)}
        objects = []
        for acc in sorted(lock):
            v = lock[acc]
            tt = meta[acc]["translation_table"]
            code = int(tt) if tt not in ("", "None") else 11
            add([ids.contig_uid(r.seq) for r in fasta.iter_records(Path(v["file"]))], code, "prokaryote")
            objects.append({"path": v["file"], "source_name": "ncbi_genomes",
                            "source_release": release, "license": "INSDC public data",
                            "redistributable": True, "retrieved_at": v["retrieved_utc"]})
    (out / "genomes.json").write_text(json.dumps(genomes))
    (out / "objects.json").write_text(json.dumps(objects, indent=1))
    contigs = len({u for g in genomes for u in g["contig_uids"]})
    rep.update({"genomes": len(genomes), "contigs": contigs,
                **({"contigs_not_listed_again": dropped} if drop else listing.report()), "objects": len(objects),
                "codes": {str(c): sum(1 for g in genomes if g["genetic_code"] == c)
                          for c in sorted({g["genetic_code"] for g in genomes})}})
    print(json.dumps(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
