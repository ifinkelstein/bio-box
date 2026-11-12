"""The contig uids of the stages a new source could repeat, read-only (S2 inputs, 2026-11-11).

    python s2_stage_uids.py ROOT R50_GENOMES OUT_TSV

Stored stages: every `stage=<S>` partition of `ROOT/l0seq/genome_input`, its `contig_uid` column (the
contigs of the stage's genomes, including one an earlier stage stored first; `contig_index` lists only
the stage's new sequences), from the files `biobox.io.part_files` lists (the code's own reader of a
partition); length -1 (not in that table). Nothing in ROOT is written. Fetched stages: for each `R50_GENOMES/<S>/lock.json`, the
FASTA files it lists (sorted by accession), each record's `ids.contig_uid` (BLAKE3-128 of the
uppercase sequence, `src/biobox/ids.py`) as `fasta.iter_records` reads it, which is what ingest stores.

Writes OUT_TSV (`stage`, `how` store|fasta, `contig_uid` hex, `length`, `accession` or empty), sorted,
and prints a JSON report: rows and distinct uids per (stage, how), and for a stage both stored and
fetched, whether the two sets are equal.
"""

import json
import sys
from pathlib import Path

import pyarrow.parquet as pq

from biobox import fasta, ids, io


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    root, r50, out = Path(argv[0]), Path(argv[1]), Path(argv[2])
    rows: list[tuple[str, str, str, int, str]] = []
    idx = root / "l0seq" / "genome_input"
    for part in sorted(p for p in idx.iterdir() if p.name.startswith("stage=")):
        stage = part.name.split("=", 1)[1]
        for f in io.part_files(part):
            rows += [(stage, "store", u.hex(), -1, "")
                     for u in pq.read_table(f, columns=["contig_uid"]).column("contig_uid").to_pylist()]
    for lock in sorted(r50.glob("*/lock.json")):
        stage = lock.parent.name
        genomes = json.loads(lock.read_text())["genomes"]
        for acc in sorted(genomes):
            for r in fasta.iter_records(Path(genomes[acc]["file"])):
                rows.append((stage, "fasta", ids.contig_uid(r.seq).hex(), len(r.seq), r.accession))
    rows.sort()
    tmp = out.with_name(f".{out.name}.tmp")
    with open(tmp, "w") as w:
        w.write("stage\thow\tcontig_uid\tlength\taccession\n")
        for r in rows:
            w.write("\t".join(map(str, r)) + "\n")
    tmp.replace(out)
    sets: dict[tuple[str, str], set[str]] = {}
    count: dict[tuple[str, str], int] = {}
    for s, h, u, _, _ in rows:
        sets.setdefault((s, h), set()).add(u)
        count[(s, h)] = count.get((s, h), 0) + 1
    rep = {"out": str(out), "sets": {f"{s}/{h}": {"rows": count[(s, h)], "distinct": len(v)}
                                     for (s, h), v in sorted(sets.items())},
           "store_equals_fasta": {s: sets[(s, "store")] == sets[(s, "fasta")]
                                  for s, h in sorted(sets) if h == "store" and (s, "fasta") in sets}}
    print(json.dumps(rep, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
