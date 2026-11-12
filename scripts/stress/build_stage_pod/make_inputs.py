"""The `--genomes` and `--objects` files of `bio-box build-stage` for card C-10, from the import ladder's inputs.

    python make_inputs.py phage SRC_FASTA_GZ N OUT_DIR [--drop-repeats]
    python make_inputs.py prok LOCK_JSON SELECTION_TSV OUT_DIR [--source-release datasets-v2-YYYY-MM-DD] [--drop-repeats]

--source-release labels every prok object (default datasets-v2-2026-11-09, the ladder's fetch date).
--drop-repeats keeps the rule these inputs were made with before DEC-27: a contig uid already listed by an earlier
genome of the stage is not listed again, and a genome left with none is not listed. The inputs of stages S1, P1 to
P6 and H1 were made with it, and it gives them again byte for byte.

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

import gzip
import json
import sys
from pathlib import Path

from biobox import fasta, ids

USAGE = __doc__.split("\n\n")[1]
HELP = "\n\n".join(__doc__.split("\n\n")[1:3])


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


def main(argv: list[str]) -> int:
    if any(a in ("-h", "--help") for a in argv):
        print(HELP)
        return 0
    release = "datasets-v2-2026-11-09"
    if "--source-release" in argv:
        k = argv.index("--source-release")
        if k + 1 >= len(argv):
            print(USAGE, file=sys.stderr)
            return 2
        release = argv[k + 1]
        argv = argv[:k] + argv[k + 2 :]
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
