"""Step 6, laptop side: export every gene slice of the built golden store for an independent check.

    python step6_export.py ROOT OUT_DIR

Writes regions_plus.txt / regions_minus.txt (samtools region strings on the shard's record names),
expected_plus.tsv / expected_minus.tsv (sha256 of the bytes bio-box returns for each region, same
order) and genes.tsv (accession, start, end, strand, 0-based half-open) for the gene-call comparison.
Origin-wrapping genes are listed in genes.tsv but left out of the slice files.
"""

import hashlib
import sys
from pathlib import Path

import pyarrow.parquet as pq

from biobox import io
from biobox.l0seq import SeqStore
from biobox.l1.genecall import GENE_CALL_VERSION

if len(sys.argv) < 3 or sys.argv[1] in ("-h", "--help"):
    print(__doc__)
    sys.exit(0 if sys.argv[1:2] in (["-h"], ["--help"]) else 2)

root, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
store = SeqStore(root)
occ = [r for f in io.part_files(root / f"l1/protein_occurrence/gene_call_version={GENE_CALL_VERSION}/stage=S0")
       for r in pq.read_table(f).to_pylist()]
co = [r for f in io.part_files(root / "l0seq/contig_occurrence", "stage=S0/source=*/part-*.parquet")
      for r in pq.read_table(f).to_pylist()]
acc = {r["contig_uid"]: r["accession"] for r in co}
occ.sort(key=lambda r: (acc[r["contig_uid"]], r["slot_ordinal"]))
files = {s: (open(out / f"regions_{n}.txt", "w"), open(out / f"expected_{n}.tsv", "w")) for s, n in ((1, "plus"), (-1, "minus"))}
n = wraps = 0
with open(out / "genes.tsv", "w") as genes:
    for r in occ:
        cu = r["contig_uid"]
        rd = store.default_digest(cu)
        length = store.length_of(rd)
        genes.write(f"{acc[cu]}\t{r['start']}\t{r['end']}\t{r['strand']}\t{length}\n")
        if r["end"] > length:
            wraps += 1
            continue
        seq = store.get(cu, r["start"], r["end"], r["strand"])
        reg, exp = files[r["strand"]]
        reg.write(f"{rd.hex()}:{r['start'] + 1}-{r['end']}\n")
        exp.write(f"{rd.hex()}:{r['start'] + 1}-{r['end']}\t{hashlib.sha256(seq).hexdigest()}\n")
        n += 1
for a, b in files.values():
    a.close()
    b.close()
print(f"genes {len(occ)}, slices exported {n}, origin-wrapping skipped {wraps}")
