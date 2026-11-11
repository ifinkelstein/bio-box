"""The `build-stage` inputs of stage H1 (HPGC huge phages), in the form `make_inputs.py phage` gives S1.

    python s2_h1_inputs.py SRC_FASTA_GZ FASTA_OUT OUT_DIR RETRIEVED_AT

SRC_FASTA_GZ is decompressed byte for byte to FASTA_OUT (as `make_inputs.py phage` copies INPHARED;
written under a temporary name and renamed; an existing FASTA_OUT is kept and only hashed). Every record
is its own genome through `make_inputs.genome` itself (imported from this script's directory): kind
virus, genetic code 11, assembly quality complete, no circular contig (FASTA carries no topology).
`kind virus` is what makes gene calling use `gv-meta` (`biobox.l1.genecall._finder`), as for S1.

The one source object: FASTA_OUT, source `hpgc_huge_phages` (the `sources.yaml` entry), release
`figshare:10.6084/m9.figshare.30506864.v1` (the record's DOI and version), licence CC BY 4.0,
redistributable, retrieved at RETRIEVED_AT (ISO 8601, the end of the download).

`make_inputs.genome` does not list a contig uid already listed (a record whose sequence repeats an
earlier record's). DEC-27 says the script must stop doing that; until the code that reads such inputs
exists (C-21a to C-21c) it still does, so every record it drops is listed in OUT_DIR/dropped_records.tsv
(ordinal, accession, contig_uid, first ordinal and accession of that sequence) and counted in the
report. Writes OUT_DIR/{genomes.json, objects.json, dropped_records.tsv} under OUT_DIR.tmp, renamed
when complete; refuses an existing OUT_DIR. Prints a JSON report.
"""

import gzip
import json
import shutil
import sys
from pathlib import Path

from biobox import fasta, ids

import make_inputs

SOURCE = {"source_name": "hpgc_huge_phages", "source_release": "figshare:10.6084/m9.figshare.30506864.v1",
          "license": "CC BY 4.0", "redistributable": True}


def main(argv: list[str]) -> int:
    if len(argv) != 4:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    src, fa, out, retrieved = Path(argv[0]), Path(argv[1]), Path(argv[2]), argv[3]
    if out.exists():
        print(f"refused: {out} exists", file=sys.stderr)
        return 2
    if not fa.exists():
        fa.parent.mkdir(parents=True, exist_ok=True)
        tmp = fa.with_name(f".{fa.name}.tmp")
        with gzip.open(src, "rb") as r, open(tmp, "wb") as w:
            shutil.copyfileobj(r, w, 1 << 22)
        tmp.replace(fa)
    work = out.with_name(f"{out.name}.tmp")
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    seen: set = set()
    first: dict[bytes, tuple[int, str]] = {}
    genomes, dropped = [], ["ordinal\taccession\tcontig_uid\tfirst_ordinal\tfirst_accession"]
    records = bp = 0
    for i, r in enumerate(fasta.iter_records(fa)):
        uid = ids.contig_uid(r.seq)
        g, d = make_inputs.genome([uid], seen, 11, "virus", "complete")
        if g:
            genomes.append(g)
            first[uid] = (i, r.accession)
        else:
            j, a = first[uid]
            dropped.append(f"{i}\t{r.accession}\t{uid.hex()}\t{j}\t{a}")
        records += 1
        bp += len(r.seq)
    objects = [{"path": str(fa), **SOURCE, "retrieved_at": retrieved}]
    (work / "genomes.json").write_text(json.dumps(genomes))
    (work / "objects.json").write_text(json.dumps(objects, indent=1))
    (work / "dropped_records.tsv").write_text("\n".join(dropped) + "\n")
    work.rename(out)
    rep = {"records": records, "bp": bp, "genomes": len(genomes),
           "contigs": sum(len(g["contig_uids"]) for g in genomes), "contigs_not_listed_again": len(dropped) - 1,
           "objects": len(objects), "fasta": str(fa), "fasta_bytes": fa.stat().st_size,
           "fasta_sha256": ids.sha256_file(fa).hex(),
           "codes": {str(c): sum(1 for g in genomes if g["genetic_code"] == c)
                     for c in sorted({g["genetic_code"] for g in genomes})},
           "sha256": {n: ids.sha256_file(out / n).hex() for n in ("genomes.json", "objects.json",
                                                                 "dropped_records.tsv")}}
    print(json.dumps(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
