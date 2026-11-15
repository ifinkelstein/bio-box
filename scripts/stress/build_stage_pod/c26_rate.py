"""Card C-26 part 1: the gene-calling rate of the PLSDB stage's genomes, measured with the code's own call path.

    python c26_rate.py FASTA IN OUT_JSON [EVERY]

IN is the stage's inputs directory (`make_inputs.py plasmid ... --root ROOT`). Every EVERY-th record (default 30) of
the records whose contig the stage stores (the first record of its sequence, owned by no earlier stage) is called in
this one process and thread with `genecall._call` on its GenomeInput as genomes.json gives it (kind plasmid: meta
mode, `genecall._finder`; a circular contig with the wrap flanks), as the gene-call step calls each genome. The CPU
time of the calls alone (getrusage of this process around each call) is summed with the bases and genes, over all
and by topology, and written to OUT_JSON. Nothing is written elsewhere.
"""
import csv
import json
import resource
import sys
from pathlib import Path

from biobox import fasta, ids
from biobox.l0seq.inputs import GenomeInput
from biobox.l1 import genecall


def cpu() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


fa, inp, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
every = int(sys.argv[4]) if len(sys.argv) > 4 else 30


def tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


skip = {r["accession"] for r in tsv(inp / "repeated_records.tsv")} | {r["accession"] for r in tsv(inp / "shared_records.tsv")}
circular = {bytes.fromhex(u) for g in json.loads((inp / "genomes.json").read_text()) for u in g["circular"]}
tot = {k: {"genomes": 0, "bases": 0, "genes": 0, "cpu_s": 0.0} for k in ("all", "circular", "linear")}
k = 0
for r in fasta.iter_records(fa):
    acc = r.header.split()[0]
    if acc in skip:
        continue
    k += 1
    if k % every:
        continue
    cu = ids.contig_uid(r.seq)
    g = GenomeInput((cu,), 11, kind="plasmid", circular=(cu,) if cu in circular else ())
    seqs = {cu: r.seq.upper().decode("ascii")}
    t0 = cpu()
    calls = genecall._call(g, seqs)
    dt = cpu() - t0
    for key in ("all", "circular" if cu in circular else "linear"):
        t = tot[key]
        t["genomes"] += 1
        t["bases"] += len(r.seq)
        t["genes"] += len(calls)
        t["cpu_s"] += dt
for t in tot.values():
    if t["bases"]:
        t["s_per_mbp"] = round(t["cpu_s"] / t["bases"] * 1e6, 4)
        t["genes_per_mbp"] = round(t["genes"] / t["bases"] * 1e6, 3)
    t["cpu_s"] = round(t["cpu_s"], 3)
res = {"every": every, "new_first_records_seen": k, "modes": sorted({genecall._finder(
    GenomeInput((b"x" * 16,), 11, kind="plasmid"), {b"x" * 16: "ACGT"})[1]}), **tot}
out.write_text(json.dumps(res, indent=1))
print(json.dumps(res))
