"""Card B-48: the lengths of stage PL1's contigs from its FASTA, and the bases card C-26's rate sample called.

    python b48_pl1_lengths.py FASTA IN OUT_DIR [EVERY]

IN is the stage's inputs directory (`make_inputs.py plasmid ... --root ROOT`). Reads FASTA once, record by record:
- every record's contig uid (`ids.contig_uid`, the store's) and length, written once per contig to
  OUT_DIR/lengths.parquet (contig_uid, contig_len), for the dry plan of the stage's gene calling (PL1 is not
  ingested, so no `contig_index` holds these lengths);
- the records `c26_rate.py` called (every EVERY-th, default 30, of the records whose contig the stage stores: not a
  repeat, not owned by an earlier stage), with their topology from genomes.json: genomes, bases and the bases the
  caller reads (a circular contig with its wrap flanks, `genecall._genes`: n + 2 min(n, WRAP_FLANK)), over all and
  by topology, to OUT_DIR/sample.json. The selection is c26_rate.py's, line for line.
Nothing is written elsewhere.
"""
import csv
import json
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from biobox import fasta, ids
from biobox.l1 import genecall

fa, inp, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
every = int(sys.argv[4]) if len(sys.argv) > 4 else 30
out.mkdir(parents=True, exist_ok=True)


def tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


skip = {r["accession"] for r in tsv(inp / "repeated_records.tsv")} | {r["accession"] for r in tsv(inp / "shared_records.tsv")}
circular = {bytes.fromhex(u) for g in json.loads((inp / "genomes.json").read_text()) for u in g["circular"]}
tot = {k: {"genomes": 0, "bases": 0, "called": 0} for k in ("all", "circular", "linear")}
lengths: dict[bytes, int] = {}
k = records = 0
for r in fasta.iter_records(fa):
    records += 1
    cu = ids.contig_uid(r.seq)
    n = len(r.seq)
    if lengths.setdefault(cu, n) != n:
        sys.exit(f"contig {cu.hex()} has two lengths")
    acc = r.header.split()[0]
    if acc in skip:
        continue
    k += 1
    if k % every:
        continue
    circ = cu in circular
    called = n + (2 * min(n, genecall.WRAP_FLANK) if circ else 0)
    for key in ("all", "circular" if circ else "linear"):
        t = tot[key]
        t["genomes"] += 1
        t["bases"] += n
        t["called"] += called
uids = sorted(lengths)
pq.write_table(pa.table({"contig_uid": pa.array(uids, pa.binary(16)),
                         "contig_len": pa.array([lengths[u] for u in uids], pa.uint64())}), out / "lengths.parquet")
res = {"every": every, "records": records, "contigs": len(lengths), "new_first_records_seen": k,
       "wrap_flank": genecall.WRAP_FLANK, **tot}
(out / "sample.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res))
