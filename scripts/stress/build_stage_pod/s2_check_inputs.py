"""Dry check of a stage's `build-stage` inputs with the `inputs` and `ingest` steps' own readers; writes nothing.

    python s2_check_inputs.py STAGE GENOMES_JSON OBJECTS_JSON

`bio-box stage inputs` has no dry-run flag: it parses `--genomes` with `biobox.cli._genomes` (JSON list
to `GenomeInput`), then `l0seq.inputs.write_genome_input` refuses a partition already written, checks
the genomes in `_rows` (no contig twice in a genome, circular contigs inside their genome, no genome
listed twice, and a contig in several genomes only when they agree on its kind, genetic code and
topology: DEC-27, card C-21a) and writes. This runs the parsing and `_rows` only, never `write_genome_input`, so no
root is read or written. It also parses `--objects` as `stage ingest` does (`cli._source_objects`),
checks each object path is a file, checks the stage name (`catalog.check_stages`), each genome's kind
and assembly quality against the `genome_input` schema's enums, and asks `genecall._finder` which
gene-calling mode the first genome of each kind gets. Prints a JSON report; exit 0 when every check
passes, else 1.
"""

import json
import sys
from pathlib import Path

from biobox import catalog, cli
from biobox.l0seq import inputs
from biobox.l1 import genecall
from biobox.registry.schemas import table


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    stage, gpath, opath = argv
    problems = []
    catalog.check_stages([stage])
    genomes = cli._genomes(gpath)
    rows = inputs._rows(genomes)
    objects = cli._source_objects(opath)
    enums = table("genome_input").enums
    for g in genomes:
        for col, val in (("kind", g.kind), ("assembly_quality", g.assembly_quality)):
            if col in enums and val not in enums[col]:
                problems.append(f"genome {g.contig_uids[0].hex()}: {col} {val!r} not in {enums[col]}")
    for o in objects:
        if not Path(o.path).is_file():
            problems.append(f"object {o.path} is not a file")
    modes = {}
    for g in genomes:
        if g.kind not in modes:
            modes[g.kind] = genecall._finder(g, {u: "" for u in g.contig_uids})[1] if g.kind == "virus" else None
    rep = {"stage": stage, "genomes": len(genomes), "genome_input_rows": len(rows),
           "distinct_genome_uids": len({r[1] for r in rows}), "circular_rows": sum(1 for r in rows if r[2]),
           "kinds": sorted({g.kind for g in genomes}), "codes": sorted({g.genetic_code for g in genomes}),
           "assembly_quality": sorted({g.assembly_quality for g in genomes}), "gene_call_mode": modes,
           "objects": [{"path": str(o.path), "source_name": o.source_name, "source_release": o.source_release,
                        "license": o.license, "redistributable": o.redistributable,
                        "retrieved_at": o.retrieved_at.isoformat(), "bytes": Path(o.path).stat().st_size}
                       for o in objects],
           "biobox": cli.__file__, "problems": problems}
    print(json.dumps(rep, indent=1))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
