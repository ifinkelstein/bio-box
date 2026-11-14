"""Card F-02 session 10: build a full `micro` world (S1-S11) and summarize it, for the before-and-after table.

    python f0210_world_summary.py build OUT SEED [NULL]
    python f0210_world_summary.py summary WORLD OUT_JSON

build: S1-S11 of the `micro` preset (the code of this tree) at SEED into OUT, one stage at a time, and prints one
JSON line per stage with its wall time and the stage digest, then the total and the world digest. NULL (N0-full or
N0-strict) builds that null world instead (same stages, null preset).
summary: counts from the world's files only: genes per genome and the share of each family kind (backbone),
families by kind, plasmids by type, phages and their genes, prophage instances, arrays and spacers, fusion events by
kind, call errors by kind, assembly classes, contigs, sources by kind, the rows of every TSV of every stage (truth
tables included), the stage digests and, when OUT/timings.json exists, the build time per stage.
"""

from __future__ import annotations

import collections
import csv
import json
import sys
import time
from pathlib import Path

USAGE = __doc__.split("\n\n")[1]


def _tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build(out: str, seed: str, null: str | None = None) -> int:
    from biobox.synthworld import assembly, backbone, emit, fusions, intergenic, mobile, nulls, sequences, sources

    o = Path(out)
    params = backbone.load_preset("micro")
    if null:
        params = nulls.null_preset(params, null)
    s = int(seed)
    marker = params.get("null")

    def emit_stage():
        if marker is None:
            return emit.build(o)
        return emit.build(o, null={"world": marker["world"], "null_seed": s, "relabel": [],
                                   "scope": nulls.scope_n0(o, marker)})

    stages = (("S1-S2 backbone", lambda: backbone.build(params, s, o)), ("S3-S4 sequences", lambda: sequences.build(o)),
              ("S5 replicons", lambda: intergenic.build(o)), ("S6 mobile", lambda: mobile.build(o)),
              ("S7 fusions", lambda: fusions.build(o)), ("S8-S9 assembly", lambda: assembly.build(o)),
              ("S10 sources", lambda: sources.build(o)), ("S11 emit", emit_stage))
    times = {}
    total = 0.0
    last = None
    for name, fn in stages:
        t0 = time.perf_counter()
        last = fn()
        dt = time.perf_counter() - t0
        total += dt
        times[name] = round(dt, 2)
        print(json.dumps({"stage": name, "seconds": round(dt, 2), "digest": getattr(last, "digest", "")}), flush=True)
    times["S1-S11"] = round(total, 2)
    (o / "timings.json").write_text(json.dumps(times, indent=1) + "\n")
    print(json.dumps({"total_seconds": round(total, 2), "world_digest": last.world_digest}), flush=True)
    return 0


def summary(world: str, out_json: str) -> int:
    w = Path(world)
    b = w / "backbone"
    fams = {int(r["family"]): r for r in _tsv(b / "families.tsv")}
    genes = _tsv(b / "genes.tsv")
    genomes = _tsv(b / "genomes.tsv")
    reps = _tsv(b / "replicons.tsv")
    n = len(genomes)
    kind = collections.Counter()
    for g in genes:
        k = fams[int(g["family"])]["kind"]
        kind["vehicle" if k in ("vehicle_chassis", "cargo") else "background" if k in ("root", "originated") else k] += 1
    sizes = collections.Counter(g["genome_id"] for g in genes)
    chrom = [int(r["n_genes"]) for r in reps if r["kind"] == "chromosome"]
    doc = {
        "genomes": n,
        "genes_per_genome": {"mean": round(len(genes) / n, 1), "min": min(sizes.values()), "max": max(sizes.values())},
        "chromosome_genes": {"mean": round(sum(chrom) / len(chrom), 1), "min": min(chrom), "max": max(chrom)},
        "genes_per_genome_by_kind": {k: round(v / n, 1) for k, v in sorted(kind.items())},
        "planted_share": round((kind["planted"] + kind["profile_only"] + kind["decoy"]) / len(genes), 4),
        "families_with_extant_copies_by_kind": dict(sorted(collections.Counter(
            r["kind"] for r in fams.values() if int(r["n_extant_copies"]) > 0).items())),
        "plasmid_genomes_by_type": dict(sorted(collections.Counter(r["name"] for r in reps
                                                                   if r["kind"] == "plasmid").items())),
    }
    m = w / "mobile"
    ph = _tsv(m / "phages.tsv")
    doc["phages"] = {"n": len(ph), "genes": sorted(int(r["n_genes"]) for r in ph)}
    doc["prophage_instances"] = len(_tsv(m / "prophages.tsv"))
    doc["arrays"] = len({r["array_id"] for r in _tsv(m / "arrays.tsv")})
    rows = {}
    kinds = {}
    for p in sorted(w.rglob("*.tsv")):
        rel = str(p.relative_to(w))
        if rel.startswith(("zombi2/", "mobile/phages/")):
            continue
        t = _tsv(p)
        rows[rel] = len(t)
        for col in ("kind", "assembly_class"):
            if t and col in t[0] and rel.split("/")[0] in ("fusions", "assembly", "sources", "mobile"):
                kinds[f"{rel}:{col}"] = dict(sorted(collections.Counter(r[col] for r in t).items()))
    doc["kinds"] = kinds
    doc["tsv_rows"] = rows
    man = w / "manifest.json"
    if man.exists():
        mj = json.loads(man.read_text())
        doc["stage_digests"] = {d: v["digest"] for d, v in mj["stages"].items()}
        doc["world_digest"] = mj["world_digest"]
    tj = w / "timings.json"
    if tj.exists():
        doc["seconds"] = json.loads(tj.read_text())
    Path(out_json).write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: doc[k] for k in ("genes_per_genome", "planted_share", "plasmid_genomes_by_type")}))
    return 0


def main(argv: list[str]) -> int:
    if len(argv) >= 3 and argv[0] == "build":
        return build(*argv[1:4])
    if len(argv) == 3 and argv[0] == "summary":
        return summary(argv[1], argv[2])
    raise SystemExit(USAGE)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
