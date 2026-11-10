"""The `--genomes` and `--objects` files of `bio-box build-stage` for card C-10, from the import ladder's inputs.

    python make_inputs.py phage SRC_FASTA_GZ N OUT_DIR
    python make_inputs.py prok LOCK_JSON SELECTION_TSV OUT_DIR

phage: the first N records of the INPHARED FASTA, copied byte for byte (header and sequence lines as in
the source) to OUT_DIR/phage_first<N>.fa; every record its own genome (kind virus, genetic code 11,
linear: FASTA carries no topology), as the ladder's stage S1 was described (step.py, do_genome_input).
prok: one genome per accession of the ladder's lock.json (sorted), its contigs in record order, kind
prokaryote, assembly_quality complete, linear, genetic code from the selection's translation_table (11
when absent), as the ladder's stage P1 was described. The source objects are the ladder's files, read
in place. Contig uids are ids.contig_uid of each record's sequence; a uid seen before in the stage is
not listed again (counted in the report). Writes OUT_DIR/genomes.json and OUT_DIR/objects.json and
prints a JSON report.
"""

import gzip
import json
import sys
from pathlib import Path

from biobox import fasta, ids

USAGE = __doc__.split("\n\n")[1]


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
    kept = [u for u in dict.fromkeys(uids) if u not in seen]
    seen.update(kept)
    dropped = len(uids) - len(kept)
    if not kept:
        return None, dropped
    return {"contig_uids": [u.hex() for u in kept], "genetic_code": code, "kind": kind, "circular": [],
            "assembly_quality": quality}, dropped


def main(argv: list[str]) -> int:
    if any(a in ("-h", "--help") for a in argv) or len(argv) != 4 or argv[0] not in ("phage", "prok"):
        print(USAGE, file=sys.stderr if len(argv) != 4 else sys.stdout)
        return 0 if any(a in ("-h", "--help") for a in argv) else 2
    out = Path(argv[3])
    out.mkdir(parents=True, exist_ok=True)
    seen: set = set()
    genomes, dropped, rep = [], 0, {}
    if argv[0] == "phage":
        n = int(argv[2])
        fa = out / f"phage_first{n}.fa"
        rep["records_copied"] = first_records(Path(argv[1]), n, fa)
        for r in fasta.iter_records(fa):
            g, d = genome([ids.contig_uid(r.seq)], seen, 11, "virus", "complete")
            dropped += d
            genomes += [g] if g else []
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
            g, d = genome([ids.contig_uid(r.seq) for r in fasta.iter_records(Path(v["file"]))], seen, code,
                          "prokaryote", "complete")
            dropped += d
            genomes += [g] if g else []
            objects.append({"path": v["file"], "source_name": "ncbi_genomes",
                            "source_release": "datasets-v2-2026-11-09", "license": "INSDC public data",
                            "redistributable": True, "retrieved_at": v["retrieved_utc"]})
    (out / "genomes.json").write_text(json.dumps(genomes))
    (out / "objects.json").write_text(json.dumps(objects, indent=1))
    rep.update({"genomes": len(genomes), "contigs": sum(len(g["contig_uids"]) for g in genomes),
                "contigs_not_listed_again": dropped, "objects": len(objects),
                "codes": {str(c): sum(1 for g in genomes if g["genetic_code"] == c)
                          for c in sorted({g["genetic_code"] for g in genomes})}})
    print(json.dumps(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
