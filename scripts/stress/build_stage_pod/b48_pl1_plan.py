"""Card B-48: the estimate and the batch plan of stage PL1's gene calling on its real inputs, without building it.

    python b48_pl1_plan.py ROOT IN LENGTHS_PARQUET OUT_JSON [BUDGET]

PL1 is not ingested, so no `contig_index` of it holds its contigs' lengths; they come from LENGTHS_PARQUET
(`b48_pl1_lengths.py`, every record of the FASTA by `ids.contig_uid`). Everything else is the code's own: each genome of
IN/genomes.json in genome_uid order (the order `iter_genome_input` streams `genome_input` in), a contig owned by an
earlier stage when the default variant of it is in some stage's `contig_index` of ROOT (DEC-26: PL1's ingest stores
only the others), a contig of several PL1 genomes owned by the least of them (DEC-27; the inputs hold none), each
genome's estimate `genecall.genome_cost_us` (0 when it is not called, card B-48 item 5), and the cut
`genecall.cut_plan` at SHARD_MARGIN of BUDGET (default 300) CPU-minutes, as build-stage makes it when the step starts.
The plan (each batch's genomes, bases and estimate; the merge's) is written to OUT_JSON and summarised on stdout.
Reads ROOT and IN only.
"""
import json
import sys
from pathlib import Path

import duckdb
import pyarrow.parquet as pq

from biobox import ids, io
from biobox import stage as bs
from biobox.l0seq.store import INDEX_DIR
from biobox.l1 import genecall

root, inp, lengths_path, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
budget = int(sys.argv[5]) if len(sys.argv) > 5 else 300
t = pq.read_table(lengths_path)
length = dict(zip(t.column("contig_uid").to_pylist(), t.column("contig_len").to_pylist()))
genomes = []
for g in json.loads((inp / "genomes.json").read_text()):
    cus = [bytes.fromhex(u) for u in g["contig_uids"]]
    genomes.append((ids.genome_uid(cus), sorted(cus), {bytes.fromhex(u) for u in g.get("circular", ())}, g["kind"]))
genomes.sort()
con = duckdb.connect(config={"threads": 2, "memory_limit": "8GB"})
files = [f for p in sorted((root / INDEX_DIR).iterdir()) for f in io.part_files(p)]
lits = "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"
owned = {r[0] for r in con.execute(f"SELECT contig_uid FROM read_parquet({lits}, hive_partitioning = false) "
                                   "WHERE is_default_variant").fetchall()}
least: dict[bytes, bytes] = {}
for guid, cus, _, _ in genomes:
    for c in cus:
        least[c] = min(guid, least.get(c, guid))


def items():
    for guid, cus, circ, kind in genomes:
        contigs = [(length[c], c in circ) for c in cus]
        bases = sum(n for n, _ in contigs)
        mine = any(c not in owned and least[c] == guid for c in cus)
        if not mine:
            yield guid, 0, bases, 0.0, False
        else:
            yield guid, genecall.genome_cost_us(kind, contigs), bases, genecall.genome_genes(kind, bases), True


plan = genecall.cut_plan(items, bs.SHARD_MARGIN * budget)
plan.update({"stage": "PL1", "gene_call_version": genecall.GENE_CALL_VERSION, "budget": budget,
             "owned_contigs_in_root": sum(1 for _, cus, _, _ in genomes for c in cus if c in owned),
             "circular_bases_called": sum(genecall.called_bases(length[c], True) for g in genomes for c in g[1]
                                          if c in g[2] and c not in owned)})
out.write_text(json.dumps(plan, indent=1))
print(json.dumps({k: plan[k] for k in ("genomes", "genomes_not_called", "bases", "bases_called", "estimate_cpu_min",
                                       "merge_estimate_cpu_min", "genes_expected", "owned_contigs_in_root",
                                       "largest_genome")}))
for b in plan["batches"]:
    print(json.dumps({k: b[k] for k in ("batch", "first", "genomes", "bases", "estimate_cpu_min")}))
print("refusal:", bs.genecall_refusal(type("P", (), {"cpu_min_per_process": budget})(), plan))
