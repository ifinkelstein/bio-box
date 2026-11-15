"""Card B-48: the gene-calling rate of a built stage's genomes, the calls alone, as card C-26's `c26_rate.py` measured
the PLSDB sample (read-only on the root).

    python b48_call_rate.py ROOT STAGE OUT_JSON [EVERY]

Every EVERY-th genome (default 30) of the stage's `genome_input`, in genome_uid order (`iter_genome_input`), is
called in this one process and thread with `genecall._call` (the step's call path: `_finder`'s mode, the wrap flanks
of a circular contig), its sequences read first from the store of every stage on disk (DEC-26). The CPU time of the
calls alone (getrusage of this process around each call; the reading, the rows, the sorting and the writing of the
step are not in it) is summed with the genomes, bases, called bases (n + 2 min(n, WRAP_FLANK) for a circular contig)
and genes, by mode and topology, and written to OUT_JSON. Nothing is written elsewhere.
"""
import json
import resource
import sys
from pathlib import Path

from biobox.l0seq import store as st
from biobox.l0seq.inputs import iter_genome_input
from biobox.l0seq.store import SeqStore
from biobox.l1 import genecall


def cpu() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


root, stage, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
every = int(sys.argv[4]) if len(sys.argv) > 4 else 30
stages = [stage, *(s for s in st.stages_on_disk(root) if s != stage)]
store = SeqStore(root, stages)
tot: dict[str, dict] = {}
k = 0
for g in iter_genome_input(root, stage):
    k += 1
    if k % every:
        continue
    seqs = genecall._sequences(store, g)
    mode = genecall.mode_of(g.kind, sum(len(s) for s in seqs.values())) if hasattr(genecall, "mode_of") else None
    t0 = cpu()
    calls = genecall._call(g, seqs)
    dt = cpu() - t0
    mode = calls[0]["mode"] if calls else (mode or "none")
    topo = "circular" if g.circular else "linear"
    for key in ("all", f"{mode}/{topo}"):
        t = tot.setdefault(key, {"genomes": 0, "bases": 0, "called": 0, "genes": 0, "cpu_s": 0.0})
        t["genomes"] += 1
        t["bases"] += sum(len(s) for s in seqs.values())
        t["called"] += sum(len(s) + (2 * min(len(s), genecall.WRAP_FLANK) if cu in g.circular else 0)
                           for cu, s in seqs.items())
        t["genes"] += len(calls)
        t["cpu_s"] += dt
store.close()
for t in tot.values():
    t["cpu_s"] = round(t["cpu_s"], 3)
    t["s_per_mbp_called"] = round(t["cpu_s"] / t["called"] * 1e6, 4) if t["called"] else None
    t["genes_per_mbp"] = round(t["genes"] / t["bases"] * 1e6, 3) if t["bases"] else None
res = {"stage": stage, "every": every, "genomes_seen": k, **tot}
out.write_text(json.dumps(res, indent=1))
print(json.dumps(res))
