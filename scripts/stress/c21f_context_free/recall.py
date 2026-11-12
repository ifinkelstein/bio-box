"""Re-call the genes of one stage under gene-call policy p2 or p3, read-only (card C-21f).

    python recall.py ROOT STAGE OUT_DIR {p2|p3} [--tmp DIR] [--shard I/N]

p2 is today's policy, through `biobox.l1.genecall._finder` itself: viral genomes gv-meta, plasmid genomes and
genomes under 100 kb meta, other genomes single mode trained on all their contigs together. Re-calling p2 is the
control: its genes must equal the stored `feature` rows of the stage exactly, or the environment differs from the
one that built the stage and nothing below means anything.

p3 is the design note's context-free policy (`planning/design/C-21_shared_contigs.md` §2b, A3): viral genomes
gv-meta as today; any other contig of at least 100 kb (`genecall.SINGLE_MODE_MIN_BP`) single mode trained on that
contig alone with its genome's genetic code; shorter contigs meta. The genes of a contig are found exactly as the
product does (`genecall._genes`: the circular-contig flank, the 0-based half-open coordinates).

Nothing is written under ROOT (bind it read-only). Writes, in OUT_DIR, one genome at a time:

    genes.parquet    contig_uid, start, end, strand, mode, genetic_code     (one row per gene)
    contigs.parquet  genome_uid, contig_uid, contig_len, is_circular, kind, genome_bases, n_contigs, is_largest,
                     mode, train_cpu_s, find_cpu_s                           (one row per contig)
    done.json        genomes, contigs, genes, CPU seconds by phase and by mode, process rusage

`train_cpu_s` is the thread CPU time (`time.thread_time`) of building and training the finder: per contig under p3
single; under p2 single the genome's training time, on its largest contig's row (0 on the others). `find_cpu_s` is
the thread CPU time of `find_genes` and the coordinate conversion on that contig. A meta finder is built once per
genome and costs nothing to build, as in `genecall._finder`.

`--shard I/N` calls only the genomes whose position in `genome_uid` order is I modulo N (0-based), so a stage can
be split across processes; the comparison reads every shard's directory.
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import sys
import time
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pyrodigal
import pyrodigal_gv

from biobox import ids
from biobox.l0seq.inputs import iter_genome_input
from biobox.l0seq.store import SeqStore
from biobox.l1 import genecall

GENES = pa.schema([("contig_uid", pa.binary(16)), ("start", pa.int64()), ("end", pa.int64()),
                   ("strand", pa.int8()), ("mode", pa.string()), ("genetic_code", pa.uint8())])
CONTIGS = pa.schema([("genome_uid", pa.binary(16)), ("contig_uid", pa.binary(16)), ("contig_len", pa.int64()),
                     ("is_circular", pa.bool_()), ("kind", pa.string()), ("genome_bases", pa.int64()),
                     ("n_contigs", pa.int32()), ("is_largest", pa.bool_()), ("mode", pa.string()),
                     ("train_cpu_s", pa.float64()), ("find_cpu_s", pa.float64())])


def _find(finder, seq: str, circular: bool) -> tuple[list[tuple], float]:
    t = time.thread_time()
    genes = sorted(genecall._genes(finder, seq, circular), key=lambda g: (g[0], g[1], g[2]))
    return genes, time.thread_time() - t


def call_p2(genome, seqs: dict[bytes, str]) -> list[tuple]:
    """(contig_uid, mode, train_s, find_s, genes) per contig, in contig_uid order, through the product's
    `_finder`: one finder per genome, trained on the whole genome in single mode."""
    t = time.thread_time()
    finder, mode = genecall._finder(genome, seqs)
    train = time.thread_time() - t if mode == "single" else 0.0
    largest = max(sorted(seqs), key=lambda cu: len(seqs[cu]))
    out = []
    for cu in sorted(seqs):
        genes, find = _find(finder, seqs[cu], cu in genome.circular)
        out.append((cu, mode, train if cu == largest else 0.0, find, genes))
    return out


def call_p3(genome, seqs: dict[bytes, str]) -> list[tuple]:
    """As `call_p2`, under p3: the mode and the training depend on the contig alone (and its genome's kind and
    code, which genomes sharing a contig must agree on, DEC-27)."""
    if genome.kind not in genecall.KINDS:
        raise ValueError(f"unknown genome kind {genome.kind!r}")
    meta = None
    out = []
    for cu in sorted(seqs):
        seq = seqs[cu]
        train = 0.0
        if genome.kind == "virus":
            if meta is None:
                meta = pyrodigal_gv.ViralGeneFinder(meta=True)
            finder, mode = meta, "gv-meta"
        elif len(seq) >= genecall.SINGLE_MODE_MIN_BP:
            t = time.thread_time()
            finder = pyrodigal.GeneFinder(meta=False)
            finder.train(seq, translation_table=genome.genetic_code)
            train = time.thread_time() - t
            mode = "single-self"
        else:
            if meta is None:
                meta = pyrodigal.GeneFinder(meta=True)
            finder, mode = meta, "meta"
        genes, find = _find(finder, seq, cu in genome.circular)
        out.append((cu, mode, train, find, genes))
    return out


POLICIES = {"p2": call_p2, "p3": call_p3}


def _rusage() -> dict:
    ru = resource.getrusage(resource.RUSAGE_SELF)
    return {"user_s": ru.ru_utime, "sys_s": ru.ru_stime, "maxrss_mib": round(ru.ru_maxrss / 1024.0, 1)}


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root", type=Path)
    ap.add_argument("stage")
    ap.add_argument("out", type=Path)
    ap.add_argument("policy", choices=sorted(POLICIES))
    ap.add_argument("--tmp", type=Path, default=None)
    ap.add_argument("--shard", default="0/1")
    a = ap.parse_args(argv)
    k, n = (int(x) for x in a.shard.split("/"))
    if not 0 <= k < n:
        raise SystemExit(f"--shard {a.shard}: need 0 <= I < N")
    a.out.mkdir(parents=True, exist_ok=False)  # a run directory is written once
    call = POLICIES[a.policy]
    store = SeqStore(a.root, [a.stage])
    gw = pq.ParquetWriter(a.out / "genes.parquet.tmp", GENES, compression="zstd")
    cw = pq.ParquetWriter(a.out / "contigs.parquet.tmp", CONTIGS, compression="zstd")
    tot = {"genomes": 0, "contigs": 0, "genes": 0, "read_cpu_s": 0.0, "train_cpu_s": {}, "find_cpu_s": {},
           "bases": {}}
    t_wall = time.monotonic()
    try:
        for i, g in enumerate(iter_genome_input(a.root, a.stage, tmp=a.tmp)):
            if i % n != k:
                continue
            t = time.thread_time()
            seqs = genecall._sequences(store, g)
            tot["read_cpu_s"] += time.thread_time() - t
            guid = ids.genome_uid(list(g.contig_uids))
            bases = sum(len(s) for s in seqs.values())
            largest = max(sorted(seqs), key=lambda cu: len(seqs[cu]))
            gcols: list[list] = [[] for _ in GENES]
            ccols: list[list] = [[] for _ in CONTIGS]
            for cu, mode, train, find, genes in call(g, seqs):
                for c, v in zip(ccols, (guid, cu, len(seqs[cu]), cu in g.circular, g.kind, bases, len(seqs),
                                        cu == largest, mode, train, find)):
                    c.append(v)
                for start, end, strand, _aa, code, _partial in genes:
                    for c, v in zip(gcols, (cu, start, end, strand, mode, code)):
                        c.append(v)
                tot["train_cpu_s"][mode] = tot["train_cpu_s"].get(mode, 0.0) + train
                tot["find_cpu_s"][mode] = tot["find_cpu_s"].get(mode, 0.0) + find
                tot["bases"][mode] = tot["bases"].get(mode, 0) + len(seqs[cu])
                tot["contigs"] += 1
                tot["genes"] += len(genes)
            gw.write_batch(pa.RecordBatch.from_arrays([pa.array(c, f.type) for c, f in zip(gcols, GENES)],
                                                      schema=GENES))
            cw.write_batch(pa.RecordBatch.from_arrays([pa.array(c, f.type) for c, f in zip(ccols, CONTIGS)],
                                                      schema=CONTIGS))
            tot["genomes"] += 1
            if tot["genomes"] % 25 == 0:
                print(f"recall: {a.policy} {a.shard} {tot['genomes']} genomes, {tot['genes']} genes", flush=True)
    finally:
        gw.close()
        cw.close()
        store.close()
    os.replace(a.out / "genes.parquet.tmp", a.out / "genes.parquet")
    os.replace(a.out / "contigs.parquet.tmp", a.out / "contigs.parquet")
    done = {"policy": a.policy, "stage": a.stage, "shard": a.shard, "wall_s": round(time.monotonic() - t_wall, 3),
            "pyrodigal": pyrodigal.__version__, "pyrodigal_gv": pyrodigal_gv.__version__,
            "host": os.uname().nodename, **tot, "rusage": _rusage()}
    tmp = a.out / "done.json.tmp"
    tmp.write_text(json.dumps(done, indent=1, sort_keys=True) + "\n")
    os.replace(tmp, a.out / "done.json")
    print(f"recall: done {json.dumps(done, sort_keys=True)}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
