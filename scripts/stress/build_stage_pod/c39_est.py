"""Card C-39 step 5: the cost of building stage V1, from the code's own functions and a gene-calling sample.

    python c39_est.py --root ROOT --inputs IN --faa VOTUS_FAA_ZST --out OUT [--every 30]

IN holds what `c39_inputs.py uhgv ... --root ROOT` wrote; VOTUS_FAA_ZST is UHGV's votus_full.faa.zst (the provider's
prodigal-gv proteins of the representatives). Nothing in ROOT is written. Writes OUT/summary.json (and prints it):

- the code's estimates (`biobox.stage`, `biobox.l1.genecall` constants and functions, read from the tree): ingest
  (`stage.ingest_cpu_min` of the FASTA's bytes), gene calling (`genecall.one_process_cpu_min` over every genome's
  bases, circular ones with the wrap flanks, in gv-meta mode, and the batches it takes at 80% of the 300 CPU-minute
  budget), the two CRISPR producers (`stage.crispr_cpu_min` of the bases the stage stores, times
  `CRISPR_NEAR_CAS_COST`), the Pfam scan (`SCAN_CPU_MIN_PER_M` per million residues, shards of
  300 x 0.8 x 10^6 / SCAN_CPU_MIN_PER_M residues);
- a sample measured here: every EVERY-th contig the stage stores (FASTA order) called in this one process and thread
  with `genecall._call` on its GenomeInput as genomes.json gives it, the CPU time of the calls alone (getrusage), genes,
  distinct proteins (`ids.protein_uid`) and their residues, and how many of those proteins the root's protein_seq
  already holds (any stage, the root's one gene-call version);
- the provider's proteins: records, distinct sequences (`ids.protein_uid` after removing a trailing '*') and their
  residues, and how many of the distinct ones the root's protein_seq already holds;
- the bytes each stage's files hold in ROOT beside the bases its contig_index stores (storage per Gbp of S1 and H1,
  the virus stages).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import resource
from collections import Counter
from compression import zstd
from pathlib import Path

import duckdb

from biobox import fasta, ids, io, stage
from biobox.l0seq import store
from biobox.l0seq.inputs import GenomeInput
from biobox.l1 import genecall

ZOPT = {zstd.DecompressionParameter.window_log_max: 31}


def cpu() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def sql_list(files) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def faa_records(path: Path):
    """(header, sequence) of a zstd FASTA of proteins, streamed."""
    head, seq = None, []
    with zstd.open(path, "rt", options=ZOPT) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if ln.startswith(">"):
                if head is not None:
                    yield head, "".join(seq)
                head, seq = ln[1:], []
            elif ln:
                seq.append(ln)
    if head is not None:
        yield head, "".join(seq)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--inputs", required=True, type=Path)
    ap.add_argument("--faa", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--every", type=int, default=30)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    out: dict = {}
    fa = a.inputs / "uhgv_v1.0_votus_full.fna"
    with open(a.inputs / "records.tsv", newline="") as fh:
        recs = list(csv.DictReader(fh, delimiter="\t"))
    with open(a.inputs / "shared_records.tsv", newline="") as fh:
        owned = {r["contig_uid"] for r in csv.DictReader(fh, delimiter="\t")}
    with open(a.inputs / "repeated_records.tsv", newline="") as fh:
        rep = {r["ordinal"] for r in csv.DictReader(fh, delimiter="\t")}
    first = [r for r in recs if r["ordinal"] not in rep]
    called = Counter()  # (mode, circular) -> bases read by the caller, every genome called whole
    for r in first:
        c = r["circular"] == "true"
        called[("gv-meta", c)] += genecall.called_bases(int(r["length"]), c)
    stored_bases = sum(int(r["length"]) for r in first if r["contig_uid"] not in owned)
    gc_min = genecall.one_process_cpu_min(dict(called))
    margin = stage.SHARD_MARGIN * 300
    shard_res = int(300 * stage.SHARD_MARGIN * 1e6 / stage.SCAN_CPU_MIN_PER_M)
    out["code"] = {
        "fasta_bytes": fa.stat().st_size, "ingest_cpu_min": round(stage.ingest_cpu_min(fa.stat().st_size), 2),
        "genome_bases": sum(int(r["length"]) for r in first), "stored_bases": stored_bases,
        "called_bases_by_topology": {str(k[1]): v for k, v in called.items()},
        "genecall_cpu_min_one_process": round(gc_min, 1), "genecall_batches_at_margin": math.ceil(gc_min / margin),
        "genes_expected_code_rate": round(genecall.GENES_PER_MBP["gv-meta"] * stored_bases / 1e6),
        "crispr_default_cpu_min": round(stage.crispr_cpu_min(stored_bases), 2),
        "crispr_filtered_cpu_min": round(stage.crispr_cpu_min(stored_bases) * stage.CRISPR_NEAR_CAS_COST, 2),
        "scan_cpu_min_per_m_residues": stage.SCAN_CPU_MIN_PER_M, "shard_residues": shard_res,
        "shard_margin": stage.SHARD_MARGIN, "call_s_per_mbp": {f"{m} {c}": v for (m, c), v in
                                                               genecall.CALL_S_PER_MBP.items() if m == "gv-meta"},
        "step_overhead": genecall.STEP_OVERHEAD}

    # the sample: every EVERY-th stored contig, called alone
    want = {r["contig_uid"] for i, r in enumerate([r for r in first if r["contig_uid"] not in owned], 1)
            if i % a.every == 0}
    genomes = {g["contig_uids"][0]: g for g in json.loads((a.inputs / "genomes.json").read_text())}
    tot = {k: {"genomes": 0, "bases": 0, "genes": 0, "cpu_s": 0.0} for k in ("all", "circular", "linear")}
    prot: dict[bytes, int] = {}
    for r in fasta.iter_records(fa):
        cu = ids.contig_uid(r.seq)
        if cu.hex() not in want:
            continue
        g = genomes[cu.hex()]
        circ = bool(g["circular"])
        gi = GenomeInput((cu,), g["genetic_code"], kind="virus", circular=(cu,) if circ else (),
                         assembly_quality=g["assembly_quality"])
        t0 = cpu()
        calls = genecall._call(gi, {cu: r.seq.upper().decode("ascii")})
        dt = cpu() - t0
        for c in calls:
            prot.setdefault(ids.protein_uid(c["aa"]), len(c["aa"]))
        for key in ("all", "circular" if circ else "linear"):
            t = tot[key]
            t["genomes"] += 1
            t["bases"] += len(r.seq)
            t["genes"] += len(calls)
            t["cpu_s"] += dt
    for t in tot.values():
        if t["bases"]:
            t["s_per_mbp"] = round(t["cpu_s"] / t["bases"] * 1e6, 4)
            t["genes_per_mbp"] = round(t["genes"] / t["bases"] * 1e6, 2)
        t["cpu_s"] = round(t["cpu_s"], 2)

    con = duckdb.connect(config={"threads": 2, "memory_limit": "6GB", "temp_directory": str(a.out / "duckdb")})
    gcv = sorted(os.listdir(a.root / "l0seq" / "protein_seq"))
    if len(gcv) != 1:
        raise SystemExit(f"gene call versions {gcv}: expected one")
    ps = [f for s in store.stages_on_disk(a.root) for f in io.part_files(
        a.root / "l0seq/protein_seq" / gcv[0] / f"stage={s}", "protein_prefix=*/part-*.parquet")]

    def known(uids) -> int:
        con.execute("CREATE OR REPLACE TABLE q (protein_uid BLOB)")
        con.executemany("INSERT INTO q VALUES (?)", [(u,) for u in uids])
        (n,), = con.execute(f"SELECT count(DISTINCT protein_uid) FROM read_parquet({sql_list(ps)}, "
                            "hive_partitioning = false) WHERE protein_uid IN (SELECT protein_uid FROM q)").fetchall()
        return n

    tot["all"].update({"distinct_proteins": len(prot), "residues": sum(prot.values()),
                       "distinct_in_root": known(sorted(prot))})
    out["sample"] = {"every": a.every, **tot}

    # the provider's proteins
    n_rec = 0
    pp: dict[bytes, int] = {}
    for _, aa in faa_records(a.faa):
        aa = aa.rstrip("*")
        n_rec += 1
        pp.setdefault(ids.protein_uid(aa), len(aa))
    out["provider_proteins"] = {"records": n_rec, "distinct": len(pp), "residues_distinct": sum(pp.values()),
                                "distinct_in_root": known(sorted(pp)), "gene_call_version": gcv[0]}
    del pp

    # bytes per stored base of each stage in the root
    size = Counter()
    for dp, _, fs in os.walk(a.root):
        for f in fs:
            p = os.path.join(dp, f)
            for part in Path(p).relative_to(a.root).parts:
                if part.startswith("stage="):
                    size[part.removeprefix("stage=")] += os.path.getsize(p)
                    break
    stored = {}
    for s in store.stages_on_disk(a.root):
        ci = io.part_files(a.root / store.INDEX_DIR / f"stage={s}")
        (b,), = con.execute(f"SELECT coalesce(sum(contig_len), 0) FROM read_parquet({sql_list(ci)}, "
                            "hive_partitioning = false) WHERE is_default_variant").fetchall()
        stored[s] = {"bytes": size[s], "stored_bases": int(b), "gb_per_gbp": round(size[s] / int(b), 4) if b else None}
    out["root_bytes_by_stage"] = stored
    (a.out / "summary.json").write_text(json.dumps(out, indent=1, sort_keys=True))
    print(json.dumps(out, sort_keys=True))


if __name__ == "__main__":
    main()
