"""Card C-26 part 2: the estimates of stage PL1's steps, from the code's own estimate functions, and the check of each
published run record against them (AGENTS §6: a measurement more than 20% from its estimate stops the card; DEC-43:
only steps estimated at 24 CPU-minutes or more are judged). In the runtime image, with the tree's src on PYTHONPATH;
reads only:

    python c26_est.py genecall --root ROOT --genomes GENOMES.JSON --fasta FASTA --stage S [--budget 300] [--out F]
    python c26_est.py before --root ROOT --profile P --inputs IN --stage PL1 [--out-plan F]
    python c26_est.py rates --root ROOT --objects OBJECTS.JSON --stage PL1 [--expect]

`genecall`: the gene-call plan of a stage that is not ingested yet (the form of card B-48's `b48_pl1_plan.py`): each
contig's length from FASTA (`ids.contig_uid` of every record, the store's), each genome of GENOMES.JSON in genome_uid
order, a contig owned by an earlier stage when the default variant of it is in some stage's `contig_index` of ROOT
(DEC-26), each genome's estimate `genecall.genome_cost_us` (0 when no contig is its own: not called, card B-48 item 5),
cut by `genecall.cut_plan` at `SHARD_MARGIN` of the budget, as build-stage cuts it when the step starts; with
`genecall_refusal`. Prints the plan's summary (one JSON document) and writes the whole plan to --out.

`before`: ingest's estimate (`stage.ingest_cpu_min` of the objects' FASTA bytes), the gene-call plan above (budget the
profile's), the bases the stage stores (the contigs no earlier stage owns: the CRISPR steps' input, `crispr_cpu_min`
and `CRISPR_NEAR_CAS_COST`), and the expectations of part 1 for what is measured only after gene calling (genes,
distinct proteins, residues and the scan, all *derived* in `planning/measurements/2026-11-15_plsdb_inputs.md` §4).

`rates`: the stage's run records that exist, each against its estimate: `ingest`; each `genecall.batch-<j>` (the
estimate in the published plan `l1/genecall_plan/stage=<S>/plan.json`) and the merge `genecall` (the plan's
`merge_estimate_cpu_min`); each `domains-pfam-38.2.shard-<k>` (its plan's residues at `SCAN_CPU_MIN_PER_M`); the two
CRISPR steps (`stage_bases` now). With --expect, also the counts part 1 estimated as ranges, once they exist: genes and
distinct proteins (the gene-call record's outputs) and the residues of the scan plan; a count outside its range by more
than 20% of the nearer end is judged like a step. Prints one JSON line; exit 3 when a judged step or count is outside
20% (exit 4 when every such step is a gene-call batch), else 0.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import duckdb
import pyarrow.parquet as pq

from biobox import fasta, ids, io
from biobox import profile as bprofile
from biobox import stage as bs
from biobox.l0seq.store import INDEX_DIR
from biobox.l1 import crispr, genecall
from biobox.l1.genecall import GENE_CALL_VERSION as GCV

LIB, VER = "pfam", "38.2"
TOL = 0.20
# DEC-43 (PI, 2026-11-15): a step estimated under this many CPU-minutes is reported but not judged (a step that small
# costs its start-up and has no rate). Every step estimated at or above it is held to 20% (AGENTS §6).
MIN_JUDGED_CPU_MIN = 24.0
# Part 1's estimates of what is known only after gene calling (planning/measurements/2026-11-15_plsdb_inputs.md §4,
# all *derived* there): genes 5.61 M (PGAP density) to 6.37 M (meta mode on every 30th record); distinct proteins
# 1.38 M (PGAP's distinct proteins, measured lower bound) to 1.76 M; residues 409.7 M to 521.2 M.
EXPECT = {"genes": (5_610_000, 6_370_000), "distinct_proteins": (1_383_442, 1_760_000),
          "residues": (409_706_637, 521_200_000)}


def rec(root: Path, stage: str, step: str) -> dict | None:
    p = root / "runs" / f"stage={stage}" / f"{step}.json"
    return json.loads(p.read_text()) if p.exists() else None


def scan_plan(root: Path, stage: str) -> dict | None:
    p = root / "l1" / "domain_scan" / f"library={LIB}-{VER}" / f"stage={stage}" / "plan.json"
    return json.loads(p.read_text()) if p.exists() else None


def gc_plan(root: Path, stage: str) -> dict | None:
    p = root / genecall.plan_record_rel(stage)
    return json.loads(p.read_text()) if p.exists() else None


def fasta_bytes(objects: Path) -> int:
    return sum(f for _, _, f in bs.object_sizes(objects))


def plan_of(root: Path, genomes_json: Path, fasta_path: Path, stage: str, budget: int) -> dict:
    length: dict[bytes, int] = {}
    records = 0
    for r in fasta.iter_records(fasta_path):
        records += 1
        cu = ids.contig_uid(r.seq)
        if length.setdefault(cu, len(r.seq)) != len(r.seq):
            raise SystemExit(f"contig {cu.hex()} has two lengths")
    genomes = []
    for g in json.loads(genomes_json.read_text()):
        cus = [bytes.fromhex(u) for u in g["contig_uids"]]
        genomes.append((ids.genome_uid(cus), sorted(cus), {bytes.fromhex(u) for u in g.get("circular", ())},
                        g["kind"]))
    genomes.sort()
    files = [f for p in sorted((root / INDEX_DIR).iterdir()) if p.name != f"stage={stage}" for f in io.part_files(p)]
    con = duckdb.connect(config={"threads": 2, "memory_limit": "8GB"})
    lits = "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"
    owned = {r[0] for r in con.execute(f"SELECT contig_uid FROM read_parquet({lits}, hive_partitioning = false) "
                                       "WHERE is_default_variant").fetchall()} if files else set()
    con.close()
    least: dict[bytes, bytes] = {}
    for guid, cus, _, _ in genomes:
        for c in cus:
            least[c] = min(guid, least.get(c, guid))

    def items():
        for guid, cus, circ, kind in genomes:
            contigs = [(length[c], c in circ) for c in cus]
            n = sum(x for x, _ in contigs)
            if not any(c not in owned and least[c] == guid for c in cus):
                yield guid, 0, n, 0.0, False
            else:
                yield guid, genecall.genome_cost_us(kind, contigs), n, genecall.genome_genes(kind, n), True

    plan = genecall.cut_plan(items, bs.SHARD_MARGIN * budget)
    stored = {c for _, cus, _, _ in genomes for c in cus if c not in owned}
    plan.update({"stage": stage, "gene_call_version": GCV, "budget": budget, "fasta_records": records,
                 "contigs_listed": len({c for _, cus, _, _ in genomes for c in cus}),
                 "contigs_owned_by_earlier_stages": sum(1 for _, cus, _, _ in genomes for c in cus if c in owned),
                 "contigs_stored": len(stored), "bases_stored": sum(length[c] for c in stored),
                 "refusal": bs.genecall_refusal(type("P", (), {"cpu_min_per_process": budget})(), plan)})
    return plan


def summary(plan: dict) -> dict:
    keys = ("stage", "budget", "genomes", "genomes_not_called", "bases", "bases_called", "contigs_listed",
            "contigs_owned_by_earlier_stages", "contigs_stored", "bases_stored", "estimate_cpu_min",
            "merge_estimate_cpu_min", "genes_expected", "largest_genome", "refusal")
    return {**{k: plan[k] for k in keys},
            "batches": [{k: b[k] for k in ("batch", "genomes", "bases", "estimate_cpu_min")} for b in plan["batches"]]}


def genecall_mode(a) -> int:
    plan = plan_of(Path(a.root), Path(a.genomes), Path(a.fasta), a.stage, a.budget)
    if a.out:
        Path(a.out).write_text(json.dumps(plan, indent=1, sort_keys=True))
    print(json.dumps(summary(plan), indent=1, sort_keys=True))
    return 0


def before(a) -> int:
    root, prof, inp = Path(a.root), bprofile.load(a.profile), Path(a.inputs)
    n = fasta_bytes(inp / "objects.json")
    fa = [p for p, _, _ in bs.object_sizes(inp / "objects.json")]
    if len(fa) != 1:
        raise SystemExit(f"expected one FASTA object, got {fa}")
    plan = plan_of(root, inp / "genomes.json", Path(fa[0]), a.stage, prof.cpu_min_per_process)
    if a.out_plan:
        Path(a.out_plan).write_text(json.dumps(plan, indent=1, sort_keys=True))
    b = plan["bases_stored"]
    doc = {
        "fasta_bytes": n,
        "ingest_cpu_min": round(bs.ingest_cpu_min(n), 3),
        "ingest_formula": f"{bs.INGEST_FIXED_CPU_MIN} + {n} x {bs.INGEST_CPU_MIN} / {bs.INGEST_BP}",
        "bases_stored": b,
        "crispr_default_cpu_min": round(bs.crispr_cpu_min(b), 3),
        "crispr_near_cpu_min": round(bs.crispr_cpu_min(b) * bs.CRISPR_NEAR_CAS_COST, 3),
        "crispr_formula": f"bases x {bs.CRISPR_CPU_S_PER_GBP} / 1e9 / 60 (x {bs.CRISPR_NEAR_CAS_COST} filtered)",
        "genecall": summary(plan),
        "margin_cpu_min": bs.SHARD_MARGIN * prof.cpu_min_per_process,
        "shard_residues": bs.shard_residues(prof),
        "expected_after_gene_calling": {
            **{k: list(v) for k, v in EXPECT.items()},
            "scan_cpu_min_before_copies": [round(bs.scan_cpu_min(x), 1) for x in EXPECT["residues"]],
            "scan_shards_before_copies": [-(-x // bs.shard_residues(prof)) for x in EXPECT["residues"]],
            "label": "derived, planning/measurements/2026-11-15_plsdb_inputs.md section 4"},
        "labels": "ingest, crispr, genecall: the code's estimate functions (estimate); expected: part 1 (derived)",
    }
    print(json.dumps(doc, indent=1, sort_keys=True))
    return 0


def rates(a) -> int:
    root, st = Path(a.root), a.stage
    rows, bad, counts = [], [], []

    def add(step: str, est: float, why: str) -> None:
        r = rec(root, st, step)
        if r is None:
            return
        ratio = r["cpu_min"] / est if est > 0 else float("inf")
        judged = est >= MIN_JUDGED_CPU_MIN
        rows.append({"step": step, "cpu_min": round(r["cpu_min"], 3), "estimate": round(est, 3),
                     "ratio": round(ratio, 3), "basis": why, "judged": judged, "wall_s": r["wall_s"],
                     "peak_rss_mib": round((r.get("peak_rss_bytes") or 0) / 2**20, 1)})
        if judged and abs(ratio - 1) > TOL:
            bad.append(step)

    def count(name: str, got: int) -> None:
        lo, hi = EXPECT[name]
        judged_bad = got < lo * (1 - TOL) or got > hi * (1 + TOL)
        counts.append({"count": name, "measured": got, "estimate_range": [lo, hi], "outside_20pct": judged_bad})
        if judged_bad:
            bad.append(name)

    add("ingest", bs.ingest_cpu_min(fasta_bytes(Path(a.objects))), "ingest_cpu_min(FASTA bytes)")
    gp = gc_plan(root, st)
    if gp is not None:
        for b in gp["batches"]:
            add(f"genecall.{genecall.batch_name(b['batch'])}", b["estimate_cpu_min"],
                f"plan: {b['genomes']} genomes, {b['bases']} bases called")
        add("genecall", gp["merge_estimate_cpu_min"], f"plan: the merge of {gp['genes_expected']} genes expected")
    p = scan_plan(root, st)
    if p is not None:
        for k, s in enumerate(p["shards"]):
            add(f"domains-{LIB}-{VER}.shard-{k:05d}", bs.scan_cpu_min(s["residues"]),
                f"{s['residues']} residues x {bs.SCAN_CPU_MIN_PER_M} / 1e6")
    if rec(root, st, "ingest") is not None:
        b = bs.stage_bases(root, st)
        add(f"crispr-{crispr.producer()}", bs.crispr_cpu_min(b), f"{b} bases x {bs.CRISPR_CPU_S_PER_GBP} / 1e9 / 60")
        add(f"crispr-{crispr.near_cas_producer(VER)}", bs.crispr_cpu_min(b) * bs.CRISPR_NEAR_CAS_COST,
            f"{b} bases x {bs.CRISPR_CPU_S_PER_GBP} / 1e9 / 60 x {bs.CRISPR_NEAR_CAS_COST}")
    if a.expect:
        g = rec(root, st, "genecall")
        if g is not None:
            outs = [o["path"] for o in g["outputs"] if o["path"].endswith(".parquet")]
            fp = [root / o for o in outs if o.startswith(f"l1/feature_protein/gene_call_version={GCV}/stage={st}/")]
            ps = [root / o for o in outs if o.startswith(f"l0seq/protein_seq/gene_call_version={GCV}/stage={st}/")]
            count("genes", sum(pq.read_metadata(f).num_rows for f in fp))
            count("distinct_proteins", sum(pq.read_metadata(f).num_rows for f in ps))
        if p is not None:
            count("residues", p["residues"])
    print(json.dumps({"stage": st, "checked": len(rows), "outside_20pct": bad, "rows": rows, "counts": counts},
                     sort_keys=True))
    if bad and all(x.startswith("genecall.batch-") for x in bad):
        return 4  # only gene-call batches: the driver lets the other batches publish, then stops before the merge
    return 3 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("genecall", "before", "rates"))
    ap.add_argument("--root", required=True)
    ap.add_argument("--stage", default="PL1")
    ap.add_argument("--objects")
    ap.add_argument("--profile")
    ap.add_argument("--inputs")
    ap.add_argument("--genomes")
    ap.add_argument("--fasta")
    ap.add_argument("--budget", type=int, default=300)
    ap.add_argument("--out")
    ap.add_argument("--out-plan")
    ap.add_argument("--expect", action="store_true")
    a = ap.parse_args()
    return {"genecall": genecall_mode, "before": before, "rates": rates}[a.mode](a)


if __name__ == "__main__":
    sys.exit(main())
