"""Card 0-08a: the report of the sample scan, from its outputs (run ON THE POD, in the image, with the tree's
`src` on PYTHONPATH; the driver runs it when every process is done):

    python rfam_sample_report.py <M> <tree> <out dir> [--partial] [--template <measurement .md>]

Reads, under M (`/tmp/if2453/m08a`): `sample/sample.json`, `sample/plan.json`, `sample/genes.jsonl` (written by
`rfam_sample.py`), and per process k with `out/s<k>.ok`: `runs/s<k>.json` (timed.py), `out/s<k>.hits.jsonl`
(`infernal_igr.py collect`) and `out/s<k>.permodel.json` (`infernal_igr.py permodel`). From the tree: the design
note's tables and rules (`planning/design/0-08_rfam15.1_models.tsv`, `0-08_model_costs.tsv`,
`0-08_model_rules.py`), the frozen set `src/biobox/l1/rfam_sets/prok-1.tsv`, the clans `Rfam-15.1.clanin`, and
`biobox.l1.ncrna.regions` (the producer's intergenic regions).

Writes to the out dir: `report.json` (every number), `results.md` (the tables of the measurement file), and
`0-08a_prok-2_proposed.tsv`. With `--template`, also the measurement file itself, the template's text between
`<!-- results:begin -->` and `<!-- results:end -->` replaced by `results.md`. Without `--partial` every planned
process must be done; with it, the report covers the processes done so far and says so.

Every number is computed from those files; the labels in `results.md` say which are measured and which derived.
Definitions:

- A hit: a line of a process's `--tblout` at or above the model's GA (`--cut_ga` reports nothing below), mapped to
  its contig (`collect`). A family has a hit when its model has one.
- Model classes: (a), (b<k>), (c<k>), (d) by `0-08_model_rules.py classify` over the 8-genome hits of
  `0-08_model_costs.tsv` (exactly as the set prok-1 was made); `prok-1` is membership of the frozen file; (a-) are
  the four eukaryotic rRNA models dropped from (a).
- Inside genes: the bases of a hit covered by the union of the gene calls of its contig (either strand), over the
  hit's length: "wholly" = 1, "mostly" >= 0.5. "Outside the regions": the hit does not lie wholly inside one
  intergenic region with the producer's 200-base margin (`ncrna.regions`), so the producer's intergenic scan
  cannot report it whole (it scans those regions with `--notrunc`).
- A locus is new for a model outside prok-1 when no hit of a prok-1 model overlaps it on the same strand.
- Cost per 10 Gbp of a set (derived, the design note's formula): sum over its models of the model's own CPU line
  in 0-09c's 8-genome runs (intergenic margin 200 `--notrunc` for an intergenic model, whole genome for a
  whole-genome model; `0-08_model_costs.tsv`) / 18,106,231 x 1e10 / 3,600; scaled to the producer's measured
  254.5 for prok-1 by the ratio of the two derived costs.
"""

from __future__ import annotations

import bisect
import csv
import importlib.util
import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

GOLDEN_BASES = 18_106_231  # the 8 genomes of 0-09c, whose per-model CPU lines price a set
PROK1_MEASURED = 254.5  # CPU-hours per 10 Gbp, prok-1 as the producer runs it (2026-11-11_ncrna_producer.md)
FULL_LIBRARY_0909C = 1910  # CPU-hours per 10 Gbp, full library whole genome (0-09c, derived from measured CPU)
PERMUTATIONS = 200
MARGIN = 200


def h10(cpu_s: float, bases: float = GOLDEN_BASES) -> float:
    return cpu_s / bases * 1e10 / 3600


def _tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader((x for x in fh if not x.startswith("#")), delimiter="\t"))


def _jsonl(path: Path) -> list[dict]:
    return [json.loads(x) for x in Path(path).read_text().splitlines() if x.strip()]


def _rules(tree: Path):
    spec = importlib.util.spec_from_file_location("model_rules", tree / "planning/design/0-08_model_rules.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _union_cover(s: int, e: int, ivs: list[tuple[int, int]]) -> int:
    """Bases of [s, e) inside the union of `ivs` (sorted by start)."""
    tot, cur_s, cur_e = 0, None, None
    for a, b in ivs:
        if b <= s or a >= e:
            continue
        a, b = max(a, s), min(b, e)
        if cur_e is None or a > cur_e:
            if cur_e is not None:
                tot += cur_e - cur_s
            cur_s, cur_e = a, b
        else:
            cur_e = max(cur_e, b)
    if cur_e is not None:
        tot += cur_e - cur_s
    return tot


def _fmt(x, nd=1):
    return f"{x:,.{nd}f}" if isinstance(x, float) else f"{x:,}"


def _table(head: list[str], rows: list[list]) -> str:
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


def report(M: Path, tree: Path, outdir: Path, partial: bool, template: Path | None) -> int:
    sys.path.insert(0, str(tree / "src"))
    from biobox.l1 import ncrna

    outdir.mkdir(parents=True, exist_ok=True)
    sample = json.loads((M / "sample/sample.json").read_text())
    plan = json.loads((M / "sample/plan.json").read_text())
    shards = plan["shards"]
    done = [s for s in shards if (M / f"out/s{s['k']:02d}.ok").exists()]
    if len(done) < len(shards) and not partial:
        print(f"only {len(done)} of {len(shards)} processes are done; use --partial", file=sys.stderr)
        return 3
    gen = {g["id"]: g for g in sample["genomes"]}
    done_ids = {gid for s in done for gid in s["genomes"]}
    genomes = [g for g in sample["genomes"] if g["id"] in done_ids]
    prok = [g for g in genomes if g["kind"] == "prokaryote"]
    phage = [g for g in genomes if g["kind"] != "prokaryote"]
    contig = {c["uid"]: (g["id"], c["len"], c["circular"]) for g in genomes for c in g["contigs"]}
    stratum = {g["id"]: f"{g['domain']}:{g['phylum']}" for g in genomes}

    # --- the runs (measured) -------------------------------------------------------------------------------
    runs = []
    for s in done:
        k = f"s{s['k']:02d}"
        r = json.loads((M / f"runs/{k}.json").read_text())
        c = json.loads((M / f"out/{k}.collect.json").read_text())
        mp = json.loads((M / f"fa/{k}.map.json").read_text())
        runs.append({"process": k, "genomes": len(s["genomes"]), "records": mp["records"],
                     "bases": mp["bases_searched"], "cpu_min": r["cpu_min"], "wall_s": r["wall_s"],
                     "peak_rss_mib": r["peak_rss_mib"], "hits": c["hits"], "families": c["families"],
                     "exit_code": r["exit_code"], "started_at": r["started_at"], "host": r["host"],
                     "rate": mp["bases_searched"] / (r["cpu_min"] * 60)})
    cpu_s = sum(r["cpu_min"] for r in runs) * 60
    bases = sum(r["bases"] for r in runs)

    # --- models --------------------------------------------------------------------------------------------
    rules = _rules(tree)
    models = {r["acc"]: r for r in _tsv(tree / "planning/design/0-08_rfam15.1_models.tsv")}
    costs = {r["acc"]: r for r in _tsv(tree / "planning/design/0-08_model_costs.tsv")}
    prok1 = {r["acc"]: int(r["whole_genome"]) for r in _tsv(tree / "src/biobox/l1/rfam_sets/prok-1.tsv")}
    clan = {}
    for line in (tree / "src/biobox/l1/rfam_sets/Rfam-15.1.clanin").read_text().splitlines():
        f = line.split()
        for name in f[1:]:
            clan[name] = f[0]
    cls = {}
    for a, m in models.items():
        grp, rule = rules.classify({"name": m["name"], "desc": m["desc"], "hits": costs[a]["hits"]})
        if a in prok1:
            cat = "prok-1"
        elif grp == "a":
            cat = "a-"
        else:
            cat = grp
        cls[a] = {"grp": grp, "rule": rule, "cat": cat}
    assert {a for a, c in cls.items() if c["grp"] in ("a", "b") and models[a]["name"] not in rules.EUK_RRNA} == set(prok1)

    permodel = defaultdict(float)
    for s in done:
        for m in json.loads((M / f"out/s{s['k']:02d}.permodel.json").read_text()):
            permodel[m["acc"]] += m["cpu_s"] or 0.0
    pm_tot = sum(permodel.values())

    # --- hits ----------------------------------------------------------------------------------------------
    hits = []
    for s in done:
        hits += _jsonl(M / f"out/s{s['k']:02d}.hits.jsonl")
    genes = defaultdict(list)
    for g in _jsonl(M / "sample/genes.jsonl"):
        if g["contig"] in contig:
            genes[g["contig"]].append((g["start"], g["end"], g["strand"], g["pfam"]))
    for v in genes.values():
        v.sort()
    regs = {}
    for cu, (gid, n, circ) in contig.items():
        regs[cu] = ncrna.regions([(a, b) for a, b, _, _ in genes.get(cu, [])], n, circ, MARGIN)
    p1_by = defaultdict(list)
    for h in hits:
        if cls[h["acc"]]["cat"] == "prok-1":
            p1_by[(h["contig"], h["strand"])].append(h)
    for h in hits:
        cu, s0, e0 = h["contig"], h["start"], h["end"]
        n = contig[cu][1]
        L = e0 - s0
        g = genes.get(cu, [])
        h["gene_frac"] = round(_union_cover(s0, e0, [(a, b) for a, b, _, _ in g]) / L, 4) if L else 0.0
        h["genes"] = [{"start": a, "end": b, "strand": st, "pfam": pf, "overlap": min(b, e0) - max(a, s0)}
                      for a, b, st, pf in g if a < e0 and b > s0]
        rs = regs[cu]
        i = bisect.bisect_right(rs, (s0, float("inf"))) - 1
        inside = i >= 0 and rs[i][0] <= s0 and e0 <= rs[i][1]
        if not inside and contig[cu][2]:  # circular: a region across the origin, or a hit past it
            inside = any(a <= s0 + n and e0 + n <= b for a, b in rs) or any(a <= s0 and e0 <= b for a, b in rs)
        h["in_region"] = inside
        h["ga"] = float(models[h["acc"]]["ga"])
        h["genome_kind"] = gen[h["genome"]]["kind"]
        same = [p for p in p1_by[(cu, h["strand"])] if p["start"] < e0 and p["end"] > s0 and p is not h]
        h["prok1_overlap"] = sorted({p["name"] for p in same})
        h["prok1_same_clan"] = sorted({p["name"] for p in same if clan.get(p["name"]) and
                                       clan.get(p["name"]) == clan.get(h["name"])})

    fam_hits = defaultdict(list)
    for h in hits:
        fam_hits[h["acc"]].append(h)
    fam_genomes = {a: {h["genome"] for h in v} for a, v in fam_hits.items()}

    # --- 1. families as genomes are added --------------------------------------------------------------------
    by_genome = defaultdict(set)
    for h in hits:
        by_genome[h["genome"]].add(h["acc"])

    def curve(ids: list[str], grid: list[int]) -> list[dict]:
        out = []
        perms = []
        for i in range(PERMUTATIONS):
            p = sorted(ids)
            random.Random(f"0-08a:curve:{i}").shuffle(p)
            perms.append(p)
        for n in grid:
            if n > len(ids):
                continue
            vals = []
            for p in perms:
                fam = set()
                for gid in p[:n]:
                    fam |= by_genome[gid]
                vals.append(len(fam))
            out.append({"genomes": n, "mean": sum(vals) / len(vals), "min": min(vals), "max": max(vals)})
        return out

    def last_step(ids: list[str]) -> float:
        tot = 0
        for i in range(PERMUTATIONS):
            p = sorted(ids)
            random.Random(f"0-08a:curve:{i}").shuffle(p)
            fam = set()
            for gid in p[:-1]:
                fam |= by_genome[gid]
            tot += len(by_genome[p[-1]] - fam)
        return tot / PERMUTATIONS

    pids, vids = [g["id"] for g in prok], [g["id"] for g in phage]
    curve_p = curve(pids, [8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    curve_v = curve(vids, [25, 50, 100, 150, 200, 250, 300])
    fam_p = set().union(*(by_genome[g] for g in pids)) if pids else set()
    fam_v = set().union(*(by_genome[g] for g in vids)) if vids else set()
    strata = defaultdict(list)
    for g in prok:
        strata[stratum[g["id"]]].append(g["id"])
    phy_rows = []
    for st in sorted(strata):
        fam = set().union(*(by_genome[g] for g in strata[st]))
        other = set().union(set(), *(by_genome[g] for s2 in strata if s2 != st for g in strata[s2]))
        phy_rows.append({"stratum": st, "genomes": len(strata[st]), "families": len(fam),
                         "only_here": sorted(models[a]["name"] for a in fam - other),
                         "d_families": sorted(models[a]["name"] for a in fam if cls[a]["cat"] == "d")})
    golden_fams = {a for a, c in costs.items() if int(c["hits"]) > 0}

    def by_cat(fams):
        c = Counter(cls[a]["cat"] if cls[a]["cat"] != "c" else "c" for a in fams)
        return {k: c.get(k, 0) for k in ("prok-1", "a-", "c", "d")}

    # --- 2. long-tail models that hit ------------------------------------------------------------------------
    def model_row(a):
        v = fam_hits[a]
        return {"acc": a, "name": models[a]["name"], "desc": models[a]["desc"], "rule": cls[a]["rule"],
                "hits": len(v), "genomes": len(fam_genomes[a]),
                "prok_genomes": len({h["genome"] for h in v if h["genome_kind"] == "prokaryote"}),
                "phage_genomes": len({h["genome"] for h in v if h["genome_kind"] != "prokaryote"}),
                "new_loci": sum(1 for h in v if not h["prok1_overlap"]),
                "same_clan_as_prok1_hit": sum(1 for h in v if h["prok1_same_clan"]),
                "clan": clan.get(models[a]["name"]), "max_score": max(h["score"] for h in v), "ga": float(models[a]["ga"]),
                "cpu_s_scan": round(permodel.get(a, 0.0), 2), "cpu_share": permodel.get(a, 0.0) / pm_tot if pm_tot else 0,
                "h10_igr200nt": h10(float(costs[a]["cpu_s_igr200nt"])), "h10_wg": h10(float(costs[a]["cpu_s_wg"])),
                "in_genes_mostly": sum(1 for h in v if h["gene_frac"] >= 0.5),
                "outside_regions": sum(1 for h in v if not h["in_region"]),
                "phyla": sorted({stratum[h["genome"]] for h in v})}

    d_hit = sorted(a for a in fam_hits if cls[a]["cat"] == "d")
    d_rows = [model_row(a) for a in d_hit]

    # --- 3. excluded models that hit -------------------------------------------------------------------------
    c_hits = sorted((h for h in hits if cls[h["acc"]]["cat"] == "c"),
                    key=lambda h: (cls[h["acc"]]["rule"], h["acc"], h["genome"], h["contig"], h["start"]))
    c_rows = [model_row(a) for a in sorted({h["acc"] for h in c_hits})]
    am_rows = [model_row(a) for a in sorted(a for a in fam_hits if cls[a]["cat"] == "a-")]

    # --- 4. rule-(b) models that never hit --------------------------------------------------------------------
    b_models = sorted(a for a, c in cls.items() if c["grp"] == "b")
    b_rule = defaultdict(lambda: {"models": 0, "hit": 0, "never": [], "cpu_s_never": 0.0, "h10_never": 0.0})
    for a in b_models:
        r = b_rule[cls[a]["rule"]]
        r["models"] += 1
        if a in fam_hits:
            r["hit"] += 1
        else:
            r["never"].append(f"{a} {models[a]['name']}")
            r["cpu_s_never"] += permodel.get(a, 0.0)
            r["h10_never"] += h10(float(costs[a]["cpu_s_wg" if prok1[a] else "cpu_s_igr200nt"]))
    a_never = sorted(a for a in prok1 if cls[a]["grp"] == "a" and a not in fam_hits)

    # --- 5. hits inside genes, by family ----------------------------------------------------------------------
    fam5 = []
    for a, v in sorted(fam_hits.items()):
        mostly = sum(1 for h in v if h["gene_frac"] >= 0.5)
        wholly = sum(1 for h in v if h["gene_frac"] >= 1.0)
        out_r = sum(1 for h in v if not h["in_region"])
        if mostly or out_r:
            fam5.append({"acc": a, "name": models[a]["name"], "cat": cls[a]["cat"], "rule": cls[a]["rule"],
                         "whole_genome_prok1": prok1.get(a), "hits": len(v), "wholly": wholly, "mostly": mostly,
                         "outside_regions": out_r, "h10_wg": h10(float(costs[a]["cpu_s_wg"])),
                         "h10_igr200nt": h10(float(costs[a]["cpu_s_igr200nt"]))})
    in_genes_tot = {"hits": len(hits), "wholly": sum(1 for h in hits if h["gene_frac"] >= 1.0),
                    "mostly": sum(1 for h in hits if h["gene_frac"] >= 0.5),
                    "outside_regions": sum(1 for h in hits if not h["in_region"])}
    wg23 = sorted(a for a, w in prok1.items() if w)
    wg23_cover = {a: {"hits": len(fam_hits.get(a, [])), "outside_regions": sum(1 for h in fam_hits.get(a, [])
                                                                                if not h["in_region"])} for a in wg23}

    # --- 6. prok-2 -------------------------------------------------------------------------------------------
    prok2 = dict(prok1)
    for a in d_hit:
        prok2[a] = 0
    newly_wg = []
    for a in sorted(prok2):
        v = fam_hits.get(a, [])
        lost = sum(1 for h in v if not h["in_region"])
        if not prok2[a] and v and 2 * lost > len(v):
            prok2[a] = 1
            newly_wg.append(a)

    def set_cost(s: dict) -> float:
        return h10(sum(float(costs[a]["cpu_s_wg" if w else "cpu_s_igr200nt"]) for a, w in s.items()))

    c1, c2 = set_cost(prok1), set_cost(prok2)
    tsv = outdir / "0-08a_prok-2_proposed.tsv"
    with open(tsv, "w") as fh:
        fh.write("# Rfam 15.1 model set prok-2, PROPOSED (card 0-08a; not frozen, not under src/; for the PI). prok-1's\n"
                 f"# {len(prok1)} models and its whole_genome flags, plus the {len(d_hit)} rule-(d) models with a hit at or above GA\n"
                 f"# in the 0-08a sample ({len(prok)} prokaryotes, {len(phage)} phages); whole_genome 1 also for the\n"
                 f"# {len(newly_wg)} member models more than half of whose sample hits lie outside every intergenic region\n"
                 f"# (margin 200). Written by scripts/stress/annotation_cost/rfam_sample_report.py. {len(prok2)} models;\n"
                 f"# whole_genome 1 marks the {sum(prok2.values())} models scanned on whole contigs.\n")
        fh.write("acc\tname\twhole_genome\n")
        for a in sorted(prok2):
            fh.write(f"{a}\t{models[a]['name']}\t{prok2[a]}\n")

    res = {
        "partial": len(done) < len(shards), "processes_done": len(done), "processes_planned": len(shards),
        "sample": {k: v for k, v in sample.items() if k != "genomes"},
        "genomes_reported": {"prokaryotes": len(prok), "phages": len(phage),
                             "bases": sum(g["bases"] for g in genomes)},
        "runs": runs, "cpu_s": cpu_s, "bases": bases,
        "bases_per_cpu_s": bases / cpu_s if cpu_s else None,
        "h10_full_library": h10(cpu_s, bases) if cpu_s else None,
        "max_process_cpu_min": max((r["cpu_min"] for r in runs), default=None),
        "permodel_cpu_s": pm_tot,
        "hits": len(hits), "families": len(fam_hits), "families_prok": len(fam_p), "families_phage": len(fam_v),
        "families_phage_only": sorted(models[a]["name"] for a in fam_v - fam_p),
        "families_by_cat": by_cat(fam_hits), "families_by_cat_prok": by_cat(fam_p), "families_by_cat_phage": by_cat(fam_v),
        "golden_families": len(golden_fams), "golden_families_hit_here": len(golden_fams & set(fam_hits)),
        "curve_prok": curve_p, "curve_phage": curve_v,
        "last_genome_new_families_prok": last_step(pids) if len(pids) > 1 else None,
        "last_genome_new_families_phage": last_step(vids) if len(vids) > 1 else None,
        "phyla": phy_rows, "d_hit": d_rows, "c_models": c_rows, "a_minus": am_rows,
        "c_hits": [{k: h[k] for k in ("acc", "name", "genome", "genome_kind", "contig", "start", "end", "strand",
                                      "score", "ga", "evalue", "trunc", "gene_frac", "genes", "in_region",
                                      "prok1_overlap", "prok1_same_clan")} | {"rule": cls[h["acc"]]["rule"],
                                                                              "stratum": stratum[h["genome"]]}
                   for h in c_hits],
        "b_rules": {k: v for k, v in sorted(b_rule.items(), key=lambda kv: int(kv[0][1:]))}, "a_in_prok1_never": a_never,
        "in_genes": in_genes_tot, "families_in_genes": fam5, "wg23": wg23_cover,
        "prok2": {"models": len(prok2), "whole_genome": sum(prok2.values()), "added": d_hit, "newly_wg": newly_wg,
                  "h10_derived_prok1": c1, "h10_derived_prok2": c2, "h10_scaled_prok2": PROK1_MEASURED * c2 / c1},
    }
    (outdir / "report.json").write_text(json.dumps(res, indent=1, sort_keys=True, default=str) + "\n")

    md = results_md(res, models, cls, prok1, sample)
    (outdir / "results.md").write_text(md)
    if template is not None and template.exists():
        t = template.read_text()
        b, e = "<!-- results:begin -->", "<!-- results:end -->"
        if b in t and e in t:
            t = t[: t.index(b) + len(b)] + "\n\n" + md + "\n" + t[t.index(e):]
            (outdir / template.name).write_text(t)
    print(json.dumps({k: res[k] for k in ("partial", "processes_done", "hits", "families", "families_by_cat",
                                          "h10_full_library", "max_process_cpu_min")}, default=str))
    print(json.dumps(res["prok2"] | {"added": len(d_hit)}, default=str))
    return 0


def results_md(res: dict, models: dict, cls: dict, prok1: dict, sample: dict) -> str:
    L = []
    w = L.append
    host = res["runs"][0]["host"] if res["runs"] else "?"
    day = min((r["started_at"] for r in res["runs"]), default="?")[:10]
    lab = f"*measured*, {host}, {day}"
    if res["partial"]:
        w(f"**PARTIAL: {res['processes_done']} of {res['processes_planned']} processes done; every number below "
          f"covers only their {res['genomes_reported']['prokaryotes']} prokaryotes and {res['genomes_reported']['phages']} phages.**\n")
    cp, cv = res["curve_prok"], res["curve_phage"]
    p2 = res["prok2"]
    w("### Summary\n")
    w(_table(["quantity", "value", "label"], [
        ["genomes scanned", f"{res['genomes_reported']['prokaryotes']} prokaryotes, {res['genomes_reported']['phages']} phages", "sample rule"],
        ["bases scanned", _fmt(res["bases"]), lab],
        ["CPU (cmsearch, all processes)", f"{res['cpu_s'] / 3600:,.2f} CPU-hours", lab + ", timed.py"],
        ["bases per CPU-second", _fmt(res["bases_per_cpu_s"] or 0, 0), "*derived*: bases / CPU-s"],
        ["full library, CPU-hours per 10 Gbp", _fmt(res["h10_full_library"] or 0, 0),
         f"*derived*: CPU-s / bases x 1e10 / 3600 (0-09c: {FULL_LIBRARY_0909C})"],
        ["whole sample's CPU against the note's 81 CPU-hour *estimate* (M1)",
         f"{res['cpu_s'] / res['bases'] * sample['bases']['all'] / 3600:,.1f} CPU-hours "
         f"({res['cpu_s'] / res['bases'] * sample['bases']['all'] / 3600 / 81 - 1:+.0%})" if res["bases"] else "-",
         "*derived*: CPU-s / bases scanned x " + _fmt(sample["bases"]["all"]) + " sample bases / 3600"
         + (" (projected from the processes done)" if res["partial"] else "")],
        ["largest process", f"{res['max_process_cpu_min']:,.1f} CPU-minutes" if res["max_process_cpu_min"] else "-", lab],
        ["hits at or above GA", _fmt(res["hits"]), lab],
        ["families with a hit (all / prokaryotes / phages)",
         f"{res['families']} / {res['families_prok']} / {res['families_phage']}", lab],
        ["... of them in prok-1 / (a-) / excluded (c) / long tail (d)",
         " / ".join(str(res["families_by_cat"][k]) for k in ("prok-1", "a-", "c", "d")), lab],
        [f"families at {cp[0]['genomes'] if cp else 8} and at {cp[-1]['genomes'] if cp else 100} prokaryotes "
         f"(mean of {PERMUTATIONS} orders)",
         f"{cp[0]['mean']:.1f} and {cp[-1]['mean']:.1f}" if cp else "-", "*derived* from measured hits"],
        ["new families from the last prokaryote added (mean)",
         f"{res['last_genome_new_families_prok']:.2f}" if res["last_genome_new_families_prok"] is not None else "-",
         "*derived*"],
        ["long-tail (d) models that hit: candidates to add", str(len(res["d_hit"])), lab],
        ["excluded (c) models that hit", f"{len(res['c_models'])} models, {len(res['c_hits'])} hits", lab],
        ["rule-(b) models of prok-1 with no hit", str(sum(len(v["never"]) for v in res["b_rules"].values())), lab],
        ["hits mostly inside genes / outside every margin-200 region",
         f"{res['in_genes']['mostly']} / {res['in_genes']['outside_regions']} of {res['in_genes']['hits']}", lab],
        ["proposed prok-2: models (whole genome)", f"{p2['models']} ({p2['whole_genome']})", "rule in item 6"],
        ["prok-2 cost, CPU-hours per 10 Gbp", f"{p2['h10_scaled_prok2']:.1f} (prok-1: {PROK1_MEASURED})",
         f"*derived*: {PROK1_MEASURED} x {p2['h10_derived_prok2']:.1f} / {p2['h10_derived_prok1']:.1f}"],
    ]))
    w("")
    w("### The sample\n")
    w(f"Seed {sample['seed']}; {sample['prokaryotes_in_stages']:,} prokaryote genomes in P1 and P2, GTDB R226 match "
      + ", ".join(f"{k} {v:,}" for k, v in sample["gtdb_match"].items())
      + f"; {sample['phyla_present']} phyla present ({len(sample['archaea_phyla'])} archaeal); phages every "
      f"{sample['phage_stride']}th of S1 by genome_uid. {sample['contigs']:,} contigs ({sample['circular_contigs']} "
      f"circular); bases: prokaryotes {sample['bases']['prokaryotes']:,}, phages {sample['bases']['phages']:,}, all "
      f"{sample['bases']['all']:,} (*measured*: contig_index lengths).\n")
    w(_table(["stratum (GTDB domain:phylum)", "genomes available", "sampled", "contigs", "bases"],
             [[c["stratum"], _fmt(c["available"]), c["sampled"], c["contigs"], _fmt(c["bases"])]
              for c in sample["composition"]]))
    w("")
    w("### The runs\n")
    w(_table(["process", "genomes", "records", "bases", "CPU-min", "wall s", "bases per CPU-s", "peak RSS MiB", "hits",
              "families"],
             [[r["process"], r["genomes"], r["records"], _fmt(r["bases"]), f"{r['cpu_min']:.2f}", f"{r['wall_s']:,.0f}",
               f"{r['rate']:,.0f}", f"{r['peak_rss_mib']:,.0f}", r["hits"], r["families"]] for r in res["runs"]]))
    w(f"\nAll {lab}, `timed.py` records `runs/s<k>.json`; bases per CPU-s *derived* (bases / CPU-s).\n")
    w("### 1. Families with a hit as genomes are added\n")
    w(f"Prokaryotes ({lab}; *derived*: the number of distinct models with a hit in the first n genomes, over "
      f"{PERMUTATIONS} orders seeded `0-08a:curve:<i>`):\n")
    w(_table(["genomes", "families (mean)", "min", "max"],
             [[c["genomes"], f"{c['mean']:.1f}", c["min"], c["max"]] for c in cp]))
    w(f"\nThe 8 golden genomes of 0-09c had {res['golden_families']} families; {res['golden_families_hit_here']} of them hit "
      f"here. By class, prokaryotes: " + ", ".join(f"{k} {v}" for k, v in res["families_by_cat_prok"].items()) + ".\n")
    w("Phages (same method):\n")
    w(_table(["genomes", "families (mean)", "min", "max"],
             [[c["genomes"], f"{c['mean']:.1f}", c["min"], c["max"]] for c in cv]))
    w(f"\nBy class, phages: " + ", ".join(f"{k} {v}" for k, v in res["families_by_cat_phage"].items()) +
      f". Families hit in phages and in no prokaryote of the sample: {', '.join(res['families_phage_only']) or 'none'}.\n")
    w("By GTDB phylum (prokaryotes):\n")
    w(_table(["phylum", "genomes", "families", "only in this phylum", "long-tail (d) families"],
             [[p["stratum"], p["genomes"], p["families"], ", ".join(p["only_here"]) or "-",
               ", ".join(p["d_families"]) or "-"] for p in res["phyla"]]))
    w("")
    w("### 2. Long-tail (d) models that hit: candidates to add\n")
    w("Hits and genomes *measured*; CPU-s is the model's own `# CPU time` line summed over this scan's processes "
      "(*measured*), its share of all models' lines; CPU-h/10 Gbp *derived* from the model's 8-genome lines "
      "(design note formula) for the intergenic (margin 200, `--notrunc`) and whole-genome scans. New loci: hits that "
      "no prok-1 model's hit overlaps on the same strand.\n")
    w(_table(["accession", "name", "description", "hits", "genomes (prok/phage)", "new loci", "max score / GA",
              "CPU-s", "share", "h/10 Gbp igr / wg", "phyla"],
             [[r["acc"], r["name"], r["desc"], r["hits"], f"{r['genomes']} ({r['prok_genomes']}/{r['phage_genomes']})",
               r["new_loci"], f"{r['max_score']:.1f} / {r['ga']:.1f}", f"{r['cpu_s_scan']:.1f}",
               f"{r['cpu_share']:.3%}", f"{r['h10_igr200nt']:.2f} / {r['h10_wg']:.2f}",
               ", ".join(p.split(":", 1)[1] for p in r["phyla"])] for r in res["d_hit"]]))
    tot_share = sum(r["cpu_share"] for r in res["d_hit"])
    w(f"\n{len(res['d_hit'])} models; together {tot_share:.2%} of this scan's per-model CPU (*derived*).\n")
    w("### 3. Excluded (c) models that hit\n")
    w(_table(["accession", "name", "rule", "hits", "genomes (prok/phage)", "new loci", "max score / GA", "share"],
             [[r["acc"], r["name"], r["rule"], r["hits"], f"{r['genomes']} ({r['prok_genomes']}/{r['phage_genomes']})",
               r["new_loci"], f"{r['max_score']:.1f} / {r['ga']:.1f}", f"{r['cpu_share']:.3%}"] for r in res["c_models"]]))
    w("\nEvery hit (*measured*; gene: the gene calls the hit overlaps, either strand, with their Pfam domains; "
      "in genes: fraction of the hit's bases covered by gene calls):\n")
    w(_table(["model", "rule", "genome", "phylum", "contig", "start-end", "strand", "score - GA", "E-value", "in genes",
              "genes overlapped (start-end strand: Pfam)", "prok-1 hits overlapping"],
             [[h["name"], h["rule"], h["genome"][:16], h["stratum"].split(":", 1)[1], h["contig"][:12],
               f"{h['start']:,}-{h['end']:,}", "+" if h["strand"] == 1 else "-", f"{h['score'] - h['ga']:+.1f}",
               f"{h['evalue']:.2g}", f"{h['gene_frac']:.2f}",
               "; ".join(f"{g['start']:,}-{g['end']:,} {'+' if g['strand'] == 1 else '-'}: {g['pfam'] or 'no Pfam'}"
                         for g in h["genes"]) or "-",
               ", ".join(h["prok1_overlap"]) or "-"] for h in res["c_hits"]]))
    w("\nThe eukaryotic rRNA models dropped by (a-), for the record:\n")
    w(_table(["accession", "name", "hits", "new loci (no prok-1 hit on the same strand)"],
             [[r["acc"], r["name"], r["hits"], r["new_loci"]] for r in res["a_minus"]]))
    w("")
    w("### 4. Rule-(b) models of prok-1 that never hit\n")
    w("CPU-s: this scan's per-model lines of the models with no hit (*measured*); h/10 Gbp: their cost in prok-1 as "
      "the producer scans them (*derived*, design-note formula).\n")
    w(_table(["rule", "models", "with a hit", "never hit", "CPU-s of those", "h/10 Gbp of those"],
             [[k, v["models"], v["hit"], len(v["never"]), f"{v['cpu_s_never']:.1f}", f"{v['h10_never']:.2f}"]
              for k, v in res["b_rules"].items()]))
    for k, v in res["b_rules"].items():
        if v["never"]:
            w(f"\n- ({k}) never hit, {len(v['never'])}: " + ", ".join(v["never"]))
    w(f"\nRule-(a) models of prok-1 (a hit in the 8 golden genomes) with no hit in this sample: {len(res['a_in_prok1_never'])}: "
      + ", ".join(f"{a} {models[a]['name']}" for a in res["a_in_prok1_never"]) + "\n")
    w("### 5. Hits inside genes, by family\n")
    w(f"Of {res['in_genes']['hits']:,} hits, {res['in_genes']['wholly']:,} lie wholly and {res['in_genes']['mostly']:,} "
      f"at least half inside gene calls; {res['in_genes']['outside_regions']:,} do not lie wholly inside an intergenic "
      "region with the 200-base margin (*measured* hits against the stage's gene calls; regions by "
      "`biobox.l1.ncrna.regions`). Families with any such hit:\n")
    w(_table(["accession", "name", "class", "whole genome in prok-1", "hits", "wholly in genes", "mostly in genes",
              "outside regions", "h/10 Gbp wg / igr"],
             [[r["acc"], r["name"], {"prok-1": f"prok-1 ({r['rule']})", "a-": "(a-)"}.get(r["cat"], r["rule"]),
               {1: "yes", 0: "no", None: "-"}[r["whole_genome_prok1"]], r["hits"], r["wholly"], r["mostly"],
               r["outside_regions"], f"{r['h10_wg']:.2f} / {r['h10_igr200nt']:.2f}"] for r in res["families_in_genes"]]))
    w("\nThe 23 whole-genome models of prok-1 in this sample:\n")
    w(_table(["accession", "name", "hits", "outside regions"],
             [[a, models[a]["name"], v["hits"], v["outside_regions"]] for a, v in res["wg23"].items()]))
    w("")
    w("### 6. Proposed prok-2\n")
    w(f"`planning/design/0-08a_prok-2_proposed.tsv`: prok-1's {len(prok1)} models with their flags, plus the "
      f"{len(p2['added'])} long-tail models of item 2; a member is scanned whole genome if prok-1 says so or if more "
      f"than half of its hits here lie outside every margin-200 region ({len(p2['newly_wg'])}: "
      + ", ".join(f"{a} {models[a]['name']}" for a in p2["newly_wg"]) + f"). {p2['models']} models, "
      f"{p2['whole_genome']} whole genome.\n")
    w(f"Cost (*derived*): CPU-hours per 10 Gbp = sum over the set of the model's 8-genome CPU line (margin 200 "
      f"`--notrunc` if intergenic, whole genome if whole genome) / 18,106,231 x 1e10 / 3,600 gives "
      f"{p2['h10_derived_prok1']:.1f} for prok-1 (the note's 247.4) and {p2['h10_derived_prok2']:.1f} for prok-2; scaled "
      f"to prok-1's measured {PROK1_MEASURED}: {PROK1_MEASURED} x {p2['h10_derived_prok2']:.1f} / "
      f"{p2['h10_derived_prok1']:.1f} = **{p2['h10_scaled_prok2']:.1f}**.\n")
    return "\n".join(L) + "\n"


def main(argv: list[str]) -> int:
    partial = "--partial" in argv
    argv = [a for a in argv if a != "--partial"]
    template = None
    if "--template" in argv:
        i = argv.index("--template")
        template = Path(argv[i + 1])
        del argv[i:i + 2]
    if len(argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    M, tree, out = Path(argv[0]), Path(argv[1]), Path(argv[2])
    if template is None:
        cands = sorted((tree / "planning/measurements").glob("*_rfam_sample.md"))
        template = cands[-1] if cands else None
    return report(M, tree, out, partial, template)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
