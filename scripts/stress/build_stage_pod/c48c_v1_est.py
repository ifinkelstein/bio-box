"""Card C-48c part 1: pre-launch estimates, the live rate guard and the rehearsal sample of stage V1 (UHGV v1.0 vOTU
representatives), in the form of card C-47's c47_est.py.

    python c48c_v1_est.py before --root ROOT --profile P --inputs IN [--stage V1] [--out-plan F] [--forbid C5]
    python c48c_v1_est.py rates  --root ROOT --objects OBJECTS.JSON [--stage V1] [--before BEFORE.JSON]
    python c48c_v1_est.py sample --inputs IN --out DIR [--every 1000] [--circular 12] [--shared-per-stage 2]

`before` (reads ROOT, writes only --out-plan): the code's own estimates of every judged step (card B-52's
`ingest_cpu_min`, card B-48's gene-call plan at the rates of card B-53, card B-53's `crispr_cpu_min`, the Pfam scan at
`SCAN_CPU_MIN_PER_M` over the residues card C-39 measured on the provider's proteins), the contigs of V1 that each
stage of the root already owns, and the comparison with card C-39's §5 (`compare_c39`). Exit 5 when a stage in
--forbid owns a contig of V1 (stage C5 is parked half built: a release cut without it would miss calls); exit 6 when
an estimate is more than 20% from C-39's and `EXPLAINED` gives no reason (AGENTS §6). The document is printed either
way.

`rates`: DEC-43 (`MIN_JUDGED_CPU_MIN = 24`) and DEC-53 (stop only when a judged step is more than 20% over its
estimate; more than 20% under is recorded, not stopped) over the stage's run records. Exit 3 to stop. There is no
list of accepted records: a stop stays a stop.

`sample`: the rehearsal's inputs, every EVERY-th genome of V1 (genomes.json order, which is the FASTA's) plus the
first CIRCULAR circular genomes not already taken and the first SHARED_PER_STAGE records of each stage that owns a
contig of V1 (`shared_records.tsv`), so that the rehearsal meets DEC-26's skips, card C-38h's owned proteins and
card B-47b's copies. The sample FASTA's records are written as the source holds them; each is checked against its
genome's `contig_uid`.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

LIB, VER = "pfam", "38.2"
TOL = 0.20
MIN_JUDGED_CPU_MIN = 24.0  # DEC-43
FORBID = ("C5",)

# Card C-39's figures (planning/measurements/2026-11-16_virus_votu_inputs.md). Counts fixed by the inputs are
# *measured* there (§3, §4); the rest are the note's *derived* estimates (§5) or the provider's proteins (*measured*,
# prodigal-gv, where our calls are pyrodigal-gv).
C39_EXACT = {
    "genomes": 168_536,  # §3, V1/report.json
    "genome_input_rows": 168_536,  # one contig per genome, §3
    "contigs": 168_536,
    "stored": 168_358,  # §4: 168,536 - 178
    "stored_by_earlier": 178,  # §4: H1 177, S1 1
    "ledger_skipped_earlier": 178,  # DEC-26: one gene-call skip per contig an earlier stage owns
    "bases_stored": 5_178_493_077,  # §4
    "circular": 22_517,  # §3
}
C39_STORED_BY = {"H1": 177, "S1": 1}  # §4
C39_BANDED = {
    "genes": 7_390_000,  # §5, derived: 1,426.36 genes per Mbp (measured sample) x 5,178.49 Mbp
    "distinct_proteins": 5_597_468,  # §5, the provider's distinct proteins (measured)
    "new_distinct_proteins": 5_422_648,  # §5: 5,597,468 - 174,820 already in the root (measured)
    "residues": 1_307_249_946,  # §5, residues of the provider's distinct proteins (measured)
}
C39_ADDED_BYTES = (3_700_000_000, 4_200_000_000)  # §5, derived: 5.178 Gbp x 0.7175 to 0.7316 GB per Gbp, + up to 10%
C39_EST = {  # §5, derived by c39_est.py from the code of 2026-11-16
    "ingest": 23.44,
    "genecall": 687.5,
    "genecall_batches": 3,
    "crispr_default": 15.71,
    "crispr_near": 18.69,
    "scan": 10_406.0,
    "scan_shards": 44,
}
C39_RESIDUES = 1_307_249_946
# Why an estimate of this code may differ from C-39's: the cards merged after C-39 changed the functions.
EXPLAINED = {
    "ingest": "card B-52 (merged 2026-11-16): 0.07 + 4.11 CPU-min per 1e9 FASTA bytes not held by another stage; "
              "C-39 took the earlier rate over every byte",
    "genecall": "card B-53 (merged 2026-11-18): the circular gv-meta rate is 4.3762 CPU-s per Mbp read, fitted on "
                "V1's own sample of 778 circular genomes, where C-39 took the estimated 6.3075; and the 178 genomes "
                "whose contig S1 or H1 owns are not called (card B-48 (5)), where C-39 called every genome",
    "genecall_batches": "follows the gene-call estimate (card B-48's cut at 80% of 300 CPU-minutes)",
    "crispr_default": "card B-53 (merged 2026-11-18): 168.1 CPU-s per Gbp of stored bases plus the circular "
                      "contigs' wrap flanks, where C-39 took the earlier rate over stored bases alone",
    "crispr_near": "as crispr_default, times CRISPR_NEAR_CAS_COST 1.15",
}


# pure functions (tests/r1/test_c48c_v1_est.py) -------------------------------------------------------------------


def judge_rate(step: str, cpu_min: float, est: float, basis: str, *, wall_s: float = 0.0,
               peak_rss_mib: float = 0.0, min_judged: float = MIN_JUDGED_CPU_MIN,
               tol: float = TOL) -> tuple[dict, bool, bool]:
    """One step under DEC-43 and DEC-53: (row, stop, recorded-under). A step estimated under `min_judged` is never
    judged. A judged step stops the build when its ratio is above 1 + tol; below 1 - tol it is recorded only."""
    ratio = cpu_min / est if est > 0 else float("inf")
    judged = est >= min_judged
    over = judged and ratio > 1.0 + tol
    under = judged and ratio < 1.0 - tol
    row = {"step": step, "cpu_min": round(cpu_min, 3), "estimate": round(est, 3), "ratio": round(ratio, 3),
           "basis": basis, "judged": judged, "over_20pct": over, "under_20pct": under, "wall_s": wall_s,
           "peak_rss_mib": peak_rss_mib}
    return row, over, under


def rates_exit(over: list[str]) -> int:
    """The driver stops on 3: any judged step more than 20% over its estimate (DEC-53)."""
    return 3 if over else 0


def owners_problem(by_stage: dict[str, int], forbidden: tuple[str, ...] = FORBID) -> str | None:
    """The launch gate: a message when a stage in `forbidden` owns at least one contig of the new stage."""
    bad = {s: n for s, n in sorted(by_stage.items()) if s in forbidden and n > 0}
    if not bad:
        return None
    return ("contigs of the stage are already owned by " + ", ".join(f"{s} ({n})" for s, n in bad.items())
            + ": do not launch (a release cut without that stage would be missing their calls)")


def band(center: float, tol: float) -> tuple[int, int]:
    return round(center * (1 - tol)), round(center * (1 + tol))


def expected_ranges(tol: float = TOL) -> dict[str, tuple[int, int]]:
    """The range of each count the check holds the built stage to: exact where the inputs fix it, +-tol around the
    C-39 note's estimate otherwise (AGENTS §6: a measurement more than 20% from its estimate is reported)."""
    r = {k: (v, v) for k, v in C39_EXACT.items()}
    r.update({k: band(v, tol) for k, v in C39_BANDED.items()})
    r["added_bytes"] = (round(C39_ADDED_BYTES[0] * (1 - tol)), round(C39_ADDED_BYTES[1] * (1 + tol)))
    return r


def check_counts(observed: dict, expected: dict[str, tuple[int, int]]) -> list[str]:
    """`name: why` for each expected count that is missing or outside its range, in name order."""
    out = []
    for k in sorted(expected):
        lo, hi = expected[k]
        if k not in observed:
            out.append(f"{k}: not observed")
        elif not lo <= observed[k] <= hi:
            out.append(f"{k}: {observed[k]} outside [{lo}, {hi}]")
    return out


def compare_c39(ours: dict[str, float], c39: dict[str, float], explained: dict[str, str] | None = None,
                tol: float = TOL) -> list[dict]:
    """Each estimate both hold: ratio ours / C-39's, whether it is more than `tol` off, and whether that is
    unexplained (no entry in `explained`)."""
    explained = EXPLAINED if explained is None else explained
    rows = []
    for k in sorted(set(ours) & set(c39)):
        ratio = ours[k] / c39[k] if c39[k] else float("inf")
        off = abs(ratio - 1.0) > tol
        rows.append({"name": k, "ours": ours[k], "c39": c39[k], "ratio": round(ratio, 4), "off_20pct": off,
                     "unexplained": off and k not in explained, "why": explained.get(k)})
    return rows


def pick_sample(n: int, circular: list[int], shared: list[int], *, every: int, n_circ: int) -> list[int]:
    """Indices of the rehearsal's genomes: every `every`-th, the first `n_circ` circular ones not already taken, and
    every index in `shared`; sorted, whatever the order of the lists given."""
    pick = set(range(0, n, every))
    extra = [i for i in sorted(set(circular)) if i not in pick][:n_circ]
    return sorted(pick | set(extra) | set(shared))


# the rest reads the root and the inputs, in the image ------------------------------------------------------------


def rec(root: Path, stage: str, step: str) -> dict | None:
    p = root / "runs" / f"stage={stage}" / f"{step}.json"
    return json.loads(p.read_text()) if p.exists() else None


def scan_plan(root: Path, stage: str) -> dict | None:
    p = root / "l1" / "domain_scan" / f"library={LIB}-{VER}" / f"stage={stage}" / "plan.json"
    return json.loads(p.read_text()) if p.exists() else None


def gc_plan(root: Path, stage: str) -> dict | None:
    from biobox.l1 import genecall

    p = root / genecall.plan_record_rel(stage)
    return json.loads(p.read_text()) if p.exists() else None


def owners(root: Path, stage: str) -> dict[bytes, str]:
    """Each contig whose default variant a stage of the root other than `stage` holds, with that stage (the first by
    name when several do)."""
    import duckdb

    from biobox import io
    from biobox.l0seq.store import INDEX_DIR

    out: dict[bytes, str] = {}
    con = duckdb.connect(config={"threads": 2, "memory_limit": "8GB"})
    for d in sorted((root / INDEX_DIR).iterdir()):
        s = d.name.removeprefix("stage=")
        files = io.part_files(d)
        if s == stage or not files:
            continue
        lits = "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"
        for (cu,) in con.execute(f"SELECT contig_uid FROM read_parquet({lits}, hive_partitioning = false) "
                                 "WHERE is_default_variant").fetchall():
            out.setdefault(cu, s)
    con.close()
    return out


def plan_of(root: Path, genomes_json: Path, objects_json: Path, stage: str, budget: int) -> tuple[dict, int]:
    """The gene-call plan build-stage will make (card B-48, as c47_est.py), what the stage stores and holds, and the
    contigs each other stage owns."""
    from biobox import fasta, ids
    from biobox import stage as bs
    from biobox.l1 import crispr, genecall
    from biobox.l1.genecall import GENE_CALL_VERSION as GCV

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
    own = owners(root, stage)
    least: dict[bytes, bytes] = {}
    for guid, cus, _, _ in genomes:
        for c in cus:
            least[c] = min(guid, least.get(c, guid))

    def items():
        for guid, cus, circ, kind in genomes:
            contigs = [(length[c], c in circ) for c in cus]
            n = sum(x for x, _ in contigs)
            if not any(c not in own and least[c] == guid for c in cus):
                yield guid, 0, n, 0.0, False
            else:
                yield guid, genecall.genome_cost_us(kind, contigs), n, genecall.genome_genes(kind, n), True

    plan = genecall.cut_plan(items, bs.SHARD_MARGIN * budget)
    all_contigs = {c for _, cus, _, _ in genomes for c in cus}
    stored = {c for c in all_contigs if c not in own}
    held = {c for c in all_contigs if c in own}
    circ_set = {c for _, _, circ, _ in genomes for c in circ if c in stored}
    plan.update({
        "stage": stage,
        "gene_call_version": GCV,
        "budget": budget,
        "fasta_records": records,
        "contigs_listed": len(all_contigs),
        "contigs_owned_by_earlier_stages": len(held),
        "owned_by_stage": dict(sorted(Counter(own[c] for c in held).items())),
        "held_bases": sum(length[c] for c in held),
        "contigs_stored": len(stored),
        "bases_stored": sum(length[c] for c in stored),
        "circular_stored": len(circ_set),
        # as `stage.crispr_flank_bases` (card B-53): 2 x min(CIRC_FLANK, contig_len) per circular contig stored
        "flank_bases": sum(2 * min(crispr.CIRC_FLANK, length[c]) for c in circ_set),
        "refusal": bs.genecall_refusal(type("P", (), {"cpu_min_per_process": budget})(), plan),
    })
    return plan, total_fasta_bytes


def summary(plan: dict) -> dict:
    keys = ("stage", "budget", "genomes", "genomes_not_called", "bases", "bases_called", "contigs_listed",
            "contigs_owned_by_earlier_stages", "owned_by_stage", "held_bases", "contigs_stored", "bases_stored",
            "circular_stored", "flank_bases", "estimate_cpu_min", "merge_estimate_cpu_min", "genes_expected",
            "largest_genome", "refusal")
    return {**{k: plan[k] for k in keys},
            "batches": [{k: b[k] for k in ("batch", "genomes", "bases", "estimate_cpu_min")} for b in plan["batches"]]}


def before(a) -> int:
    from biobox import profile as bprofile
    from biobox import stage as bs

    root, prof, inp = Path(a.root), bprofile.load(a.profile), Path(a.inputs)
    plan, n = plan_of(root, inp / "genomes.json", inp / "objects.json", a.stage, prof.cpu_min_per_process)
    if a.out_plan:
        Path(a.out_plan).write_text(json.dumps(plan, indent=1, sort_keys=True))
    b, hb, fb = plan["bases_stored"], plan["held_bases"], plan["flank_bases"]
    shard_res = bs.shard_residues(prof)
    gc_total = plan["estimate_cpu_min"]
    ours = {
        "ingest": bs.ingest_cpu_min(n, hb),
        "genecall": gc_total,
        "genecall_batches": len(plan["batches"]) or 1,
        "crispr_default": bs.crispr_cpu_min(b, fb, near_cas=False),
        "crispr_near": bs.crispr_cpu_min(b, fb, near_cas=True),
        "scan": bs.scan_cpu_min(C39_RESIDUES),
        "scan_shards": math.ceil(C39_RESIDUES / shard_res),
    }
    cmp_rows = compare_c39(ours, C39_EST)
    forbid = tuple(x for x in a.forbid.split(",") if x)
    gate = owners_problem(plan["owned_by_stage"], forbid)
    doc = {
        "fasta_bytes": n,
        "held_bases": hb,
        "ingest_cpu_min": round(ours["ingest"], 3),
        "ingest_formula": f"{bs.INGEST_FIXED_CPU_MIN} + max(0, {n} - {hb}) x {bs.INGEST_CPU_MIN} / {bs.INGEST_BP}",
        "bases_stored": b,
        "flank_bases": fb,
        "crispr_default_cpu_min": round(ours["crispr_default"], 3),
        "crispr_near_cpu_min": round(ours["crispr_near"], 3),
        "crispr_formula": f"({b} + {fb}) x {bs.CRISPR_CPU_S_PER_GBP} / 1e9 / 60 (x {bs.CRISPR_NEAR_CAS_COST} near-cas)",
        "genecall": summary(plan),
        "scan_from_c39_residues": {"residues": C39_RESIDUES, "cpu_min": round(ours["scan"], 1),
                                   "shards": ours["scan_shards"], "shard_residues": shard_res,
                                   "rate": f"{bs.SCAN_CPU_MIN_PER_M} CPU-min per 1e6 residues",
                                   "label": "derived from C-39's measured residues of the provider's proteins"},
        "margin_cpu_min": bs.SHARD_MARGIN * prof.cpu_min_per_process,
        "c39_compare": cmp_rows,
        "owners_gate": {"forbidden": list(forbid), "owned_by_stage": plan["owned_by_stage"], "problem": gate},
        "guard_policy": "DEC-43 (>= 24 CPU-min judged) + DEC-53 (stop only on ratio > 1.20; record ratio < 0.80)",
    }
    print(json.dumps(doc, indent=1, sort_keys=True))
    if gate:
        print(f"refused: {gate}", file=sys.stderr)
        return 5
    bad = [r["name"] for r in cmp_rows if r["unexplained"]]
    if bad:
        print(f"refused: estimates more than 20% from C-39's without an explanation: {bad}", file=sys.stderr)
        return 6
    return 0


def rates(a) -> int:
    from biobox import stage as bs
    from biobox.l1 import crispr, genecall

    root, st = Path(a.root), a.stage
    rows, over_list, under_list = [], [], []
    before_doc = json.loads(Path(a.before).read_text()) if a.before and Path(a.before).exists() else None

    def add(step: str, est: float, why: str) -> None:
        r = rec(root, st, step)
        if r is None:
            return
        row, over, under = judge_rate(step, r["cpu_min"], est, why, wall_s=r["wall_s"],
                                      peak_rss_mib=round((r.get("peak_rss_bytes") or 0) / 2**20, 1))
        rows.append(row)
        if over:
            over_list.append(step)
        if under:
            under_list.append(step)

    if rec(root, st, "ingest") is not None:
        if before_doc is not None and "ingest_cpu_min" in before_doc:
            add("ingest", float(before_doc["ingest_cpu_min"]), "before.json ingest_cpu_min(fasta_bytes, held_bases)")
        else:
            fa_b = sum(f for _, _, f in bs.object_sizes(Path(a.objects)))
            hb = bs.held_bases(root, st)
            add("ingest", bs.ingest_cpu_min(fa_b, hb), f"ingest_cpu_min({fa_b}, held={hb})")

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

    print(json.dumps({"stage": st, "checked": len(rows), "over_20pct": over_list, "under_20pct": under_list,
                      "rows": rows}, sort_keys=True))
    return rates_exit(over_list)


def sample(a) -> int:
    from biobox import fasta, ids

    inp, out = Path(a.inputs), Path(a.out)
    out.mkdir(parents=True, exist_ok=False)
    genomes = json.loads((inp / "genomes.json").read_text())
    with open(inp / "records.tsv", newline="") as fh:
        recs = list(csv.DictReader(fh, delimiter="\t"))
    if len(recs) != len(genomes) or any(int(r["ordinal"]) != i for i, r in enumerate(recs)):
        raise SystemExit("records.tsv is not one row per genome in genomes.json order")
    circ = [i for i, g in enumerate(genomes) if g["circular"]]
    with open(inp / "shared_records.tsv", newline="") as fh:
        by_owner: dict[str, list[int]] = {}
        for r in csv.DictReader(fh, delimiter="\t"):
            by_owner.setdefault(r["owner_stage"], []).append(int(r["ordinal"]))
    shared = [i for s in sorted(by_owner) for i in sorted(by_owner[s])[: a.shared_per_stage]]
    pick = pick_sample(len(genomes), circ, shared, every=a.every, n_circ=a.circular)
    want = set(pick)
    objs = json.loads((inp / "objects.json").read_text())
    if len(objs) != 1:
        raise SystemExit("expected one source object")
    fa_out = out / "sample.fna"
    n_bases = 0
    with open(fa_out, "wb") as fh:
        for i, r in enumerate(fasta.iter_records(Path(objs[0]["path"]))):
            if i not in want:
                continue
            if ids.contig_uid(r.seq).hex() != genomes[i]["contig_uids"][0]:
                raise SystemExit(f"record {i}: contig_uid differs from genomes.json")
            fh.write(b"\n".join(fasta.format_record(r.header, r.seq, r.line_width)) + b"\n")
            n_bases += len(r.seq)
    (out / "genomes.json").write_text(json.dumps([genomes[i] for i in pick]))
    (out / "objects.json").write_text(json.dumps([{**objs[0], "path": str(fa_out)}], indent=1))
    rep = {"genomes": len(pick), "bases": n_bases, "circular": sum(1 for i in pick if genomes[i]["circular"]),
           "shared": {s: [i for i in shared if i in by_owner[s]] for s in sorted(by_owner)},
           "every": a.every}
    print(json.dumps(rep, sort_keys=True))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("before", "rates", "sample"))
    ap.add_argument("--root")
    ap.add_argument("--stage", default="V1")
    ap.add_argument("--objects")
    ap.add_argument("--profile")
    ap.add_argument("--inputs")
    ap.add_argument("--before")
    ap.add_argument("--out-plan")
    ap.add_argument("--forbid", default=",".join(FORBID))
    ap.add_argument("--out")
    ap.add_argument("--every", type=int, default=1000)
    ap.add_argument("--circular", type=int, default=12)
    ap.add_argument("--shared-per-stage", type=int, default=2)
    a = ap.parse_args()
    return {"before": before, "rates": rates, "sample": sample}[a.mode](a)


if __name__ == "__main__":
    sys.exit(main())
