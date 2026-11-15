"""Card C-24: the estimates of stage H1's long steps, from the code's own estimate functions, and the check of each
published run record against them (AGENTS §6: a rate more than 20% from its estimate stops the card). In the image,
with the tree's src on PYTHONPATH; reads only:

    python c24_est.py before --root ROOT --profile P --objects OBJECTS.JSON --stage H1 --bases N [--like S1]
    python c24_est.py rates --root ROOT --objects OBJECTS.JSON --stage H1

`before` prints, as one JSON document: ingest's estimate (`stage.ingest_cpu_min` of the objects' FASTA bytes,
`stage.object_sizes`), the CRISPR steps' estimates at N genome bases (`stage.crispr_cpu_min`, times
`CRISPR_NEAR_CAS_COST` for the filtered one; N is the derived figure of the inputs record until the stage is ingested,
then `stage.stage_bases`), and, scaled from the built stage `--like` (its genome bases, distinct proteins and residues,
and its gene calling's and scan's run records; all *derived*: a ratio times H1's bases), the proteins, residues, scan
shards (`stage.shard_residues` of the profile) and scan CPU-minutes (`stage.scan_cpu_min`) to expect, and gene
calling's CPU-minutes at the like stage's rate per Gbp.

`rates` reads the stage's run records that exist and compares each with its estimate: `ingest` (as above), each
`domains-pfam-38.2.shard-<k>` (its plan's residues at `SCAN_CPU_MIN_PER_M`), `crispr-<default producer>` and
`crispr-<filtered producer>` (the stage's `stage_bases` now). Prints one JSON line; exit 3 when a ratio
measured / estimate is outside 1 +- 0.20, else 0.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from biobox import profile as bprofile
from biobox import stage as bs
from biobox.l1 import crispr
from biobox.l1.genecall import GENE_CALL_VERSION

LIB, VER = "pfam", "38.2"
TOL = 0.20


def rec(root: Path, stage: str, step: str) -> dict | None:
    p = root / "runs" / f"stage={stage}" / f"{step}.json"
    return json.loads(p.read_text()) if p.exists() else None


def plan(root: Path, stage: str) -> dict | None:
    p = root / "l1" / "domain_scan" / f"library={LIB}-{VER}" / f"stage={stage}" / "plan.json"
    return json.loads(p.read_text()) if p.exists() else None


def fasta_bytes(objects: Path) -> int:
    return sum(f for _, _, f in bs.object_sizes(objects))


def before(a) -> int:
    root, prof = Path(a.root), bprofile.load(a.profile)
    n = fasta_bytes(Path(a.objects))
    like = a.like
    lb = bs.stage_bases(root, like)
    lres = bs.stage_residues(root, like, GENE_CALL_VERSION)
    lplan = plan(root, like)
    lprot = lplan["proteins"] if lplan else None
    lgc = rec(root, like, "genecall")
    lscan = sum(json.loads(p.read_text())["cpu_min"]
                for p in (root / "runs" / f"stage={like}").glob(f"domains-{LIB}-{VER}.shard-*.json")
                if not p.name.endswith(".publishing.json"))
    k = a.bases / lb
    res = lres * k
    sr = bs.shard_residues(prof)
    doc = {
        "fasta_bytes": n,
        "ingest_cpu_min": round(bs.ingest_cpu_min(n), 3),
        "ingest_formula": f"{bs.INGEST_FIXED_CPU_MIN} + {n} x {bs.INGEST_CPU_MIN} / {bs.INGEST_BP}",
        "bases": a.bases,
        "crispr_default_cpu_min": round(bs.crispr_cpu_min(a.bases), 3),
        "crispr_near_cpu_min": round(bs.crispr_cpu_min(a.bases) * bs.CRISPR_NEAR_CAS_COST, 3),
        "crispr_formula": f"bases x {bs.CRISPR_CPU_S_PER_GBP} / 1e9 / 60 (x {bs.CRISPR_NEAR_CAS_COST} filtered)",
        "like": {"stage": like, "bases": lb, "residues": lres, "proteins": lprot,
                 "genecall_cpu_min": lgc["cpu_min"] if lgc else None, "scan_cpu_min_shards": round(lscan, 2),
                 "scan_rate_cpu_min_per_m": round(lscan / (lres / 1e6), 3) if lres else None},
        "scaled": {"factor": round(k, 4), "residues": round(res), "proteins": round(lprot * k) if lprot else None,
                   "shard_residues": sr, "shards": int(-(-res // sr)), "scan_cpu_min": round(bs.scan_cpu_min(res), 1),
                   "scan_rate": bs.SCAN_CPU_MIN_PER_M,
                   "genecall_cpu_min": round(lgc["cpu_min"] * k, 1) if lgc else None},
        "labels": "ingest, crispr: the code's estimate functions (estimate); scaled: derived, like-stage figure x factor",
    }
    print(json.dumps(doc, indent=1, sort_keys=True))
    return 0


def rates(a) -> int:
    root, st = Path(a.root), a.stage
    rows, bad = [], []

    def add(step: str, est: float, why: str) -> None:
        r = rec(root, st, step)
        if r is None:
            return
        ratio = r["cpu_min"] / est if est > 0 else float("inf")
        rows.append({"step": step, "cpu_min": round(r["cpu_min"], 3), "estimate": round(est, 3),
                     "ratio": round(ratio, 3), "basis": why})
        if abs(ratio - 1) > TOL:
            bad.append(step)

    add("ingest", bs.ingest_cpu_min(fasta_bytes(Path(a.objects))), "ingest_cpu_min(FASTA bytes)")
    p = plan(root, st)
    if p is not None:
        for k, s in enumerate(p["shards"]):
            add(f"domains-{LIB}-{VER}.shard-{k:05d}", bs.scan_cpu_min(s["residues"]),
                f"{s['residues']} residues x {bs.SCAN_CPU_MIN_PER_M} / 1e6")
    if rec(root, st, "ingest") is not None:
        b = bs.stage_bases(root, st)
        add(f"crispr-{crispr.producer()}", bs.crispr_cpu_min(b), f"{b} bases x {bs.CRISPR_CPU_S_PER_GBP} / 1e9 / 60")
        add(f"crispr-{crispr.near_cas_producer(VER)}", bs.crispr_cpu_min(b) * bs.CRISPR_NEAR_CAS_COST,
            f"{b} bases x {bs.CRISPR_CPU_S_PER_GBP} / 1e9 / 60 x {bs.CRISPR_NEAR_CAS_COST}")
    print(json.dumps({"stage": st, "checked": len(rows), "outside_20pct": bad, "rows": rows}, sort_keys=True))
    return 3 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("before", "rates"))
    ap.add_argument("--root", required=True)
    ap.add_argument("--objects", required=True)
    ap.add_argument("--stage", default="H1")
    ap.add_argument("--profile")
    ap.add_argument("--bases", type=int)
    ap.add_argument("--like", default="S1")
    a = ap.parse_args()
    return before(a) if a.mode == "before" else rates(a)


if __name__ == "__main__":
    sys.exit(main())
