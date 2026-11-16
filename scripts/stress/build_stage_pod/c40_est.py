"""Card C-40: the estimates of building stage PL2 from its inputs, with the code's own functions, and a measured sample.

    python c40_est.py --root ROOT --inputs IN --profile PROFILE --every N --out OUT

IN holds what `make_inputs.py plasmid-sources ... --root ROOT` wrote. Nothing in ROOT is written. Writes OUT/est.json:

- ingest: `stage.ingest_cpu_min` of each object's FASTA bytes (`stage.fasta_bytes`), the form in this tree (bytes read;
  card B-52's form, which takes off bases another stage holds, is not merged in this base);
- gene calling: the items `genecall.cut_plan` reads (one genome per distinct sequence, in genome_uid order; a genome whose
  one contig an earlier stage owns is not called, as `genecall.stage_items` says), each costed by
  `genecall.genome_cost_us` (kind plasmid: meta mode; circular contigs with their wrap flanks), cut at
  `stage.SHARD_MARGIN` of the profile's `cpu_min_per_process`, as build-stage cuts it;
- a sample: every N-th genome the stage would call (no stage owns its contig), called in this process and thread with
  `genecall._call` (the step's call path): CPU-seconds of the calls alone, genes, proteins (`ids.protein_uid` of the
  translations), distinct proteins and their residues, by source and topology; and how many of the sample's distinct
  proteins are already in an earlier stage's `protein_seq` (every stage on disk);
- from the sample: genes, distinct proteins and residues of the stage (the sample's distinct share per gene is an upper
  bound of the stage's, since a larger set repeats more), the Pfam scan (`stage.SCAN_CPU_MIN_PER_M`, shards of
  `stage.shard_residues`), both CRISPR producers (`stage.crispr_cpu_min`, `stage.CRISPR_NEAR_CAS_COST`);
- stage PL1, built: its run records' CPU-minutes by step, its genes, distinct proteins and residues, read from the root,
  for comparison.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import resource
from collections import Counter, defaultdict
from pathlib import Path

import duckdb
import pyarrow.parquet as pq

from biobox import fasta, ids, io, stage
from biobox.l0seq import store
from biobox.l0seq.inputs import GenomeInput
from biobox.l1 import genecall
from biobox.profile import load

csv.field_size_limit(1 << 30)


def cpu() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def tsv(path: Path):
    with open(path, newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)


def sql_list(files) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def pl1_measured(root: Path) -> dict:
    runs = root / "runs" / "stage=PL1"
    by = defaultdict(float)
    for f in sorted(runs.glob("*.json")):
        if f.name.endswith(".publishing.json"):
            continue
        d = json.loads(f.read_text())
        step = d.get("step", f.stem)
        key = ("genecall" if step.startswith("genecall") else "domains-pfam-38.2" if step.startswith("domains")
               else "crispr" if step.startswith("crispr") else step)
        by[key] += float(d.get("cpu_min") or 0)
    gcv = genecall.GENE_CALL_VERSION
    fp = io.part_files(root / "l1/feature_protein" / f"gene_call_version={gcv}" / "stage=PL1")
    ps = io.part_files(root / "l0seq/protein_seq" / f"gene_call_version={gcv}" / "stage=PL1",
                       "protein_prefix=*/part-*.parquet")
    ci = io.part_files(root / store.INDEX_DIR / "stage=PL1")
    con = duckdb.connect(config={"threads": 2, "memory_limit": "4GB"})
    genes = sum(pq.ParquetFile(f).metadata.num_rows for f in fp)
    (nd, res), = con.execute(f"SELECT count(*), sum(length) FROM read_parquet({sql_list(ps)}, "
                             "hive_partitioning = false)").fetchall()
    (nc, nb), = con.execute(f"SELECT count(*), sum(contig_len) FROM read_parquet({sql_list(ci)}, "
                            "hive_partitioning = false) WHERE is_default_variant").fetchall()
    return {"cpu_min_by_step": {k: round(v, 3) for k, v in sorted(by.items())}, "genes": genes,
            "distinct_proteins": nd, "residues": int(res), "stored_contigs": nc, "stored_bases": int(nb),
            "genes_per_mbp": round(genes / nb * 1e6, 2), "distinct_per_gene": round(nd / genes, 4),
            "residues_per_distinct": round(res / nd, 2),
            "genecall_cpu_s_per_mbp_stored": round(by["genecall"] * 60 / nb * 1e6, 4),
            "pfam_cpu_min_per_m_residues": round(by["domains-pfam-38.2"] / res * 1e6, 3)}


def known_proteins(root: Path, uids: set[bytes]) -> int:
    """How many of `uids` are in any stage's protein_seq on disk (default gene-call version)."""
    gcv = genecall.GENE_CALL_VERSION
    files = []
    for s in store.stages_on_disk(root):
        files += io.part_files(root / "l0seq/protein_seq" / f"gene_call_version={gcv}" / f"stage={s}",
                               "protein_prefix=*/part-*.parquet")
    con = duckdb.connect(config={"threads": 2, "memory_limit": "6GB"})
    con.execute("CREATE TABLE q (protein_uid BLOB)")
    con.executemany("INSERT INTO q VALUES (?)", [(u,) for u in sorted(uids)])
    (n,), = con.execute(f"SELECT count(DISTINCT protein_uid) FROM read_parquet({sql_list(files)}, "
                        "hive_partitioning = false) WHERE protein_uid IN (SELECT protein_uid FROM q)").fetchall()
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--inputs", required=True, type=Path)
    ap.add_argument("--profile", required=True, type=Path)
    ap.add_argument("--every", type=int, default=50)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    prof = load(a.profile)
    out: dict = {"profile_cpu_min_per_process": prof.cpu_min_per_process, "shard_margin": stage.SHARD_MARGIN}

    # ingest
    objs = []
    for o in json.loads((a.inputs / "objects.json").read_text()):
        n = stage.fasta_bytes(Path(o["path"]))
        objs.append({"object": Path(o["path"]).name, "fasta_bytes": n, "cpu_min": round(stage.ingest_cpu_min(n), 2)})
    out["ingest"] = {"objects": objs, "fasta_bytes": sum(x["fasta_bytes"] for x in objs),
                     "cpu_min_one_process": round(stage.ingest_cpu_min(sum(x["fasta_bytes"] for x in objs)), 2),
                     "form": "stage.ingest_cpu_min of FASTA bytes read (B-52 not in this base)"}

    # gene calling: one item per genome (first records), not called when an earlier stage owns the contig
    recs = [r for r in tsv(a.inputs / "records.tsv") if not r["first_accession"]]
    items = []
    for r in recs:
        cu = bytes.fromhex(r["contig_uid"])
        n, circ = int(r["length"]), r["circular"] == "true"
        called = not r["owner_stage"]
        guid = ids.genome_uid([cu])
        cost = genecall.genome_cost_us("plasmid", [(n, circ)]) if called else 0
        items.append((guid, cost, n, genecall.genome_genes("plasmid", n) if called else 0.0, called))
    items.sort(key=lambda t: t[0])
    plan = genecall.cut_plan(lambda: iter(items), stage.SHARD_MARGIN * prof.cpu_min_per_process)
    called = {(True,): 0, (False,): 0}
    for r in recs:
        if not r["owner_stage"]:
            called[(r["circular"] == "true",)] += genecall.called_bases(int(r["length"]), r["circular"] == "true")
    out["genecall"] = {"genomes": plan["genomes"], "genomes_not_called": plan["genomes_not_called"],
                       "bases_called": plan["bases_called"], "estimate_cpu_min": plan["estimate_cpu_min"],
                       "batches": len(plan["batches"]),
                       "largest_batch_cpu_min": max(b["estimate_cpu_min"] for b in plan["batches"]),
                       "largest_genome": plan["largest_genome"], "genes_expected": plan["genes_expected"],
                       "merge_estimate_cpu_min": plan["merge_estimate_cpu_min"],
                       "called_bases_with_wrap_flanks": {"circular": called[(True,)], "linear": called[(False,)]},
                       "rates": plan["rates"]}

    # the sample
    new = {r["accession"]: r for r in recs if not r["owner_stage"]}
    k = 0
    tot = defaultdict(lambda: {"genomes": 0, "bases": 0, "genes": 0, "cpu_s": 0.0})
    prot: dict[bytes, int] = {}
    prot_by_src: dict[str, dict[bytes, int]] = defaultdict(dict)
    sample_genes = 0
    for o in json.loads((a.inputs / "objects.json").read_text()):
        src = Path(o["path"]).stem
        for rec in fasta.iter_records(Path(o["path"])):
            r = new.get(rec.accession)
            if r is None or r["source"] != src:
                continue
            k += 1
            if k % a.every:
                continue
            cu = bytes.fromhex(r["contig_uid"])
            circ = r["circular"] == "true"
            g = GenomeInput((cu,), 11, kind="plasmid", circular=(cu,) if circ else ())
            t0 = cpu()
            calls = genecall._call(g, {cu: rec.seq.upper().decode("ascii")})
            dt = cpu() - t0
            sample_genes += len(calls)
            for c in calls:
                aa = c["aa"]
                u = ids.protein_uid(aa)
                prot.setdefault(u, len(aa))
                prot_by_src[src].setdefault(u, len(aa))
            for key in ("all", src, f"{src}|{'circular' if circ else 'linear'}", "circular" if circ else "linear"):
                t = tot[key]
                t["genomes"] += 1
                t["bases"] += len(rec.seq)
                t["genes"] += len(calls)
                t["cpu_s"] += dt
    for t in tot.values():
        if t["bases"]:
            t["s_per_mbp"] = round(t["cpu_s"] / t["bases"] * 1e6, 4)
            t["genes_per_mbp"] = round(t["genes"] / t["bases"] * 1e6, 2)
        t["cpu_s"] = round(t["cpu_s"], 3)
    known = known_proteins(a.root, set(prot)) if prot else 0
    out["sample"] = {"every": a.every, "genomes_to_call_seen": k, "by": dict(tot), "genes": sample_genes,
                     "distinct_proteins": len(prot), "residues_distinct": sum(prot.values()),
                     "distinct_per_gene": round(len(prot) / sample_genes, 4) if sample_genes else None,
                     "residues_per_distinct": round(sum(prot.values()) / len(prot), 2) if prot else None,
                     "distinct_in_earlier_stages": known,
                     "distinct_by_source": {s: len(v) for s, v in sorted(prot_by_src.items())}}

    # derived from the sample and the code's constants
    new_bases = sum(int(r["length"]) for r in new.values())
    s = out["sample"]["by"].get("all", {})
    genes = s.get("genes_per_mbp", 0) * new_bases / 1e6
    dpg = out["sample"]["distinct_per_gene"] or 0
    rpd = out["sample"]["residues_per_distinct"] or 0
    distinct_hi = genes * dpg
    residues_hi = distinct_hi * rpd
    sr = stage.shard_residues(prof)
    out["derived"] = {
        "new_bases": new_bases, "genes": round(genes),
        "distinct_proteins_upper": round(distinct_hi), "residues_upper": round(residues_hi),
        "pfam_cpu_min_upper": round(stage.scan_cpu_min(int(residues_hi)), 1),
        "pfam_shards_upper": math.ceil(residues_hi / sr) if residues_hi else 0, "shard_residues": sr,
        "crispr_default_cpu_min": round(stage.crispr_cpu_min(new_bases), 2),
        "crispr_near_cas_cpu_min": round(stage.crispr_cpu_min(new_bases) * stage.CRISPR_NEAR_CAS_COST, 2),
        "genecall_at_sample_rate_cpu_min": round(s.get("s_per_mbp", 0) * genecall.STEP_OVERHEAD * new_bases / 1e6 / 60, 1),
    }
    out["pl1_measured"] = pl1_measured(a.root)
    (a.out / "est.json").write_text(json.dumps(out, indent=1, sort_keys=True, default=str))
    print(json.dumps(out, sort_keys=True, default=str))


if __name__ == "__main__":
    main()
