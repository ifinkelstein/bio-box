"""Card C-47: pre-launch CPU estimates and live rate guard for stage C5 (in the form of c26_est.py).

Implements:
- B-52 `ingest_cpu_min(fasta_bytes, held_bases)`
- B-53 `crispr_cpu_min(bases_stored, flank_bases, near_cas=...)`
- DEC-43 (`MIN_JUDGED_CPU_MIN = 24.0`): steps estimated < 24 CPU-min are reported, not judged
- DEC-53: a judged step stops the build ONLY when `ratio > 1.20` (`over_20pct`); a judged step with
  `ratio < 0.80` (`under_20pct`) is recorded and reported without stopping the build.

    python c47_est.py before --root ROOT --profile P --inputs IN --stage C5 [--out-plan F]
    python c47_est.py rates  --root ROOT --objects OBJECTS.JSON --stage C5 [--before BEFORE.JSON]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import duckdb

from biobox import fasta, ids, io
from biobox import profile as bprofile
from biobox import stage as bs
from biobox.l0seq.store import INDEX_DIR
from biobox.l1 import crispr, genecall
from biobox.l1.genecall import GENE_CALL_VERSION as GCV

LIB, VER = "pfam", "38.2"
TOL = 0.20
MIN_JUDGED_CPU_MIN = 24.0  # DEC-43


def rec(root: Path, stage: str, step: str) -> dict | None:
    p = root / "runs" / f"stage={stage}" / f"{step}.json"
    return json.loads(p.read_text()) if p.exists() else None


def scan_plan(root: Path, stage: str) -> dict | None:
    p = root / "l1" / "domain_scan" / f"library={LIB}-{VER}" / f"stage={stage}" / "plan.json"
    return json.loads(p.read_text()) if p.exists() else None


def gc_plan(root: Path, stage: str) -> dict | None:
    p = root / genecall.plan_record_rel(stage)
    return json.loads(p.read_text()) if p.exists() else None


def fasta_bytes_from_objects(objects: Path) -> int:
    return sum(f for _, _, f in bs.object_sizes(objects))


def judge_rate(
    step: str,
    cpu_min: float,
    est: float,
    basis: str,
    *,
    wall_s: float = 0.0,
    peak_rss_mib: float = 0.0,
    min_judged: float = MIN_JUDGED_CPU_MIN,
    tol: float = TOL,
) -> tuple[dict, bool, bool]:
    """Evaluate one step under DEC-43 and DEC-53.

    Returns (row_dict, stop_over_20pct, note_under_20pct).
    Under DEC-53:
      - ratio > 1 + tol on a judged step (est >= min_judged) -> stop_over_20pct = True
      - ratio < 1 - tol on a judged step (est >= min_judged) -> note_under_20pct = True (no stop)
    """
    ratio = cpu_min / est if est > 0 else float("inf")
    judged = est >= min_judged
    over = judged and ratio > 1.0 + tol
    under = judged and ratio < 1.0 - tol
    row = {
        "step": step,
        "cpu_min": round(cpu_min, 3),
        "estimate": round(est, 3),
        "ratio": round(ratio, 3),
        "basis": basis,
        "judged": judged,
        "over_20pct": over,
        "under_20pct": under,
        "wall_s": wall_s,
        "peak_rss_mib": peak_rss_mib,
    }
    return row, over, under


def plan_of(root: Path, genomes_json: Path, objects_json: Path, stage: str, budget: int) -> tuple[dict, int]:
    """Build pre-launch gene-call plan and count total FASTA bytes across all objects in objects.json."""
    length: dict[bytes, int] = {}
    records = 0
    total_fasta_bytes = 0
    for p, _disk_b, fa_b in bs.object_sizes(objects_json):
        total_fasta_bytes += fa_b
        for r in fasta.iter_records(p):
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
    all_contigs = {c for _, cus, _, _ in genomes for c in cus}
    stored = {c for c in all_contigs if c not in owned}
    held = {c for c in all_contigs if c in owned}
    held_bases = sum(length[c] for c in held)
    circ_set = {c for _, _, circ, _ in genomes for c in circ if c in stored}
    flank_bases = sum(min(length[c], 2 * crispr.WRAP_FLANK) for c in circ_set)
    plan.update({
        "stage": stage,
        "gene_call_version": GCV,
        "budget": budget,
        "fasta_records": records,
        "contigs_listed": len(all_contigs),
        "contigs_owned_by_earlier_stages": len(held),
        "held_bases": held_bases,
        "contigs_stored": len(stored),
        "bases_stored": sum(length[c] for c in stored),
        "flank_bases": flank_bases,
        "refusal": bs.genecall_refusal(type("P", (), {"cpu_min_per_process": budget})(), plan),
    })
    return plan, total_fasta_bytes


def summary(plan: dict) -> dict:
    keys = (
        "stage", "budget", "genomes", "genomes_not_called", "bases", "bases_called", "contigs_listed",
        "contigs_owned_by_earlier_stages", "held_bases", "contigs_stored", "bases_stored", "flank_bases",
        "estimate_cpu_min", "merge_estimate_cpu_min", "genes_expected", "largest_genome", "refusal",
    )
    return {
        **{k: plan[k] for k in keys},
        "batches": [{k: b[k] for k in ("batch", "genomes", "bases", "estimate_cpu_min")} for b in plan["batches"]],
    }


def before(a) -> int:
    root, prof, inp = Path(a.root), bprofile.load(a.profile), Path(a.inputs)
    plan, n = plan_of(root, inp / "genomes.json", inp / "objects.json", a.stage, prof.cpu_min_per_process)
    if a.out_plan:
        Path(a.out_plan).write_text(json.dumps(plan, indent=1, sort_keys=True))
    b = plan["bases_stored"]
    hb = plan["held_bases"]
    fb = plan["flank_bases"]
    doc = {
        "fasta_bytes": n,
        "held_bases": hb,
        "ingest_cpu_min": round(bs.ingest_cpu_min(n, hb), 3),
        "ingest_formula": f"{bs.INGEST_FIXED_CPU_MIN} + max(0, {n} - {hb}) x {bs.INGEST_CPU_MIN} / {bs.INGEST_BP}",
        "bases_stored": b,
        "flank_bases": fb,
        "crispr_default_cpu_min": round(bs.crispr_cpu_min(b, fb, near_cas=False), 3),
        "crispr_near_cpu_min": round(bs.crispr_cpu_min(b, fb, near_cas=True), 3),
        "crispr_formula": f"({b} + {fb}) x {bs.CRISPR_CPU_S_PER_GBP} / 1e9 / 60 (x {bs.CRISPR_NEAR_CAS_COST} near-cas)",
        "genecall": summary(plan),
        "margin_cpu_min": bs.SHARD_MARGIN * prof.cpu_min_per_process,
        "shard_residues": bs.shard_residues(prof),
        "guard_policy": "DEC-43 (>= 24 CPU-min judged) + DEC-53 (stop only on ratio > 1.20; record ratio < 0.80)",
    }
    print(json.dumps(doc, indent=1, sort_keys=True))
    return 0


def rates(a) -> int:
    root, st = Path(a.root), a.stage
    rows, over_list, under_list = [], [], []
    before_doc = json.loads(Path(a.before).read_text()) if a.before and Path(a.before).exists() else None

    def add(step: str, est: float, why: str) -> None:
        r = rec(root, st, step)
        if r is None:
            return
        row, over, under = judge_rate(
            step,
            r["cpu_min"],
            est,
            why,
            wall_s=r["wall_s"],
            peak_rss_mib=round((r.get("peak_rss_bytes") or 0) / 2**20, 1),
        )
        rows.append(row)
        if over:
            over_list.append(step)
        if under:
            under_list.append(step)

    if rec(root, st, "ingest") is not None:
        if before_doc is not None and "ingest_cpu_min" in before_doc:
            est_in = float(before_doc["ingest_cpu_min"])
            why_in = "before.json ingest_cpu_min(fasta_bytes, held_bases)"
        else:
            fa_b = fasta_bytes_from_objects(Path(a.objects))
            hb = bs.held_bases(root, st)
            est_in = bs.ingest_cpu_min(fa_b, hb)
            why_in = f"ingest_cpu_min({fa_b}, held={hb})"
        add("ingest", est_in, why_in)

    gp = gc_plan(root, st)
    if gp is not None:
        for b in gp["batches"]:
            add(f"genecall.{genecall.batch_name(b['batch'])}", b["estimate_cpu_min"],
                f"plan: {b['genomes']} genomes, {b['bases']} bases called")
        add("genecall", gp["merge_estimate_cpu_min"], f"plan: merge of {gp['genes_expected']} genes expected")
    elif rec(root, st, "genecall") is not None and before_doc is not None:
        gc_est = before_doc.get("genecall", {}).get("estimate_cpu_min")
        if gc_est is not None:
            add("genecall", float(gc_est), "before.json one-process genecall estimate_cpu_min")

    p = scan_plan(root, st)
    if p is not None:
        for k, s in enumerate(p["shards"]):
            add(f"domains-{LIB}-{VER}.shard-{k:05d}", bs.scan_cpu_min(s["residues"]),
                f"{s['residues']} residues x {bs.SCAN_CPU_MIN_PER_M} / 1e6")

    if rec(root, st, "ingest") is not None:
        b = bs.stage_bases(root, st)
        fb = bs.crispr_flank_bases(root, st)
        add(f"crispr-{crispr.producer()}", bs.crispr_cpu_min(b, fb, near_cas=False),
            f"({b}+{fb}) bases x {bs.CRISPR_CPU_S_PER_GBP} / 1e9 / 60")
        add(f"crispr-{crispr.near_cas_producer(VER)}", bs.crispr_cpu_min(b, fb, near_cas=True),
            f"({b}+{fb}) bases x {bs.CRISPR_CPU_S_PER_GBP} x {bs.CRISPR_NEAR_CAS_COST} / 1e9 / 60")

    print(json.dumps({
        "stage": st,
        "checked": len(rows),
        "over_20pct": over_list,
        "under_20pct": under_list,
        "rows": rows,
    }, sort_keys=True))
    if over_list and all(x.startswith("genecall.batch-") for x in over_list):
        return 4
    return 3 if over_list else 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("before", "rates"))
    ap.add_argument("--root", required=True)
    ap.add_argument("--stage", default="C5")
    ap.add_argument("--objects")
    ap.add_argument("--profile")
    ap.add_argument("--inputs")
    ap.add_argument("--before")
    ap.add_argument("--out-plan")
    a = ap.parse_args()
    return {"before": before, "rates": rates}[a.mode](a)


if __name__ == "__main__":
    sys.exit(main())
