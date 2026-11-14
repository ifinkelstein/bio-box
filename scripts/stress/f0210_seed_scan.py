"""Card F-02 session 10: the `micro` backbone (S1-S2) over many seeds, measured from its emitted tables.

    python f0210_seed_scan.py scan OUT LABEL [--jobs N] [--seeds default,n2 | s1,s2,...] [--set a.b=JSON ...]
                                            [--axes R+R:S/R:S | none] [--crown AGE] [--rule R+R] [--keep]
    python f0210_seed_scan.py table OUT LABEL [LABEL ...]

scan: builds the backbone of the `micro` preset (the generator of this tree, unchanged) at each seed into
OUT/worlds/LABEL/<seed>, N seeds at once (one process each), measures it and writes one row per seed to OUT/LABEL.tsv
and the realized strength grid to OUT/LABEL_grid.tsv. Seeds: ``default`` is the preset's default seed, ``n2`` the 20
N2 seeds (`nulls.n2_seeds` of the default seed); integers are taken as they are. ``--set`` overrides an existing
preset leaf by its path (JSON value). ``--axes`` replaces `lineage_rates.axes` (``origination:1.2/duplication+
transfer+loss:0.8``; ``none`` for no axis). ``--crown`` and ``--rule`` are measurement devices that patch ZOMBI2's
calls in the worker process, never the generator's code: ``--crown`` rescales each drawn species tree to that crown
age (`_patch_crown`), ``--rule`` changes the planted and profile-only groups' rates or recipients (`_patch_rule`:
modbg, modbal, lossT, dup0, excl). The world directory is removed after it is measured unless ``--keep``.
table: prints the rows of the labels as a Markdown table.

Every number is counted here from the backbone's tables (`genes.tsv`, `families.tsv`, `groups.tsv`, `replicons.tsv`,
`truth_units.tsv`, `events.tsv`) and the species tree (`zombi2/species_complete.nwk`); the DEC-15 (c) statistics are
the generator's own (`backbone.json` ``scramble``), with the thresholds of `backbone.check_scrambled`.
"""

from __future__ import annotations

import collections
import csv
import json
import os
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

USAGE = __doc__.split("\n\n")[1]
#: §5.1 genome size: chromosome mean 1,500 +- 30%; the range's ends within 10% (min <= 1,100, max >= 3,600)
MEAN_BAND = (1050.0, 1950.0)
MIN_MAX = (1100, 3600)
#: PLAN.md §4.3 min_detectable_count: about 6-9 vehicle-independent co-occurrences (DS-08); the upper end
MIN_GENOMES = 9


def _tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write_tsv(path: Path, rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _set(params: dict, path: str, value) -> None:
    keys = path.split(".")
    d = params
    for k in keys[:-1]:
        d = d[k]
    d[keys[-1]] = value


RULES = ("modbg", "modbal", "lossT", "dup0", "excl")


def _patch_rule(rule: str | None) -> None:
    """Measurement devices applied to the planted and profile-only groups after the preset is read (the generator's
    code is not changed): ``modbg`` the run's transfer and loss rates; ``modbal`` loss = transfer + the run's
    duplication; ``lossT`` loss = transfer; ``dup0`` no duplication of their own (0 instead of the run's rate);
    ``excl`` entry exclusion as the plasmids have it: a lineage holding the group's first family receives none of
    its families by transfer (ZOMBI2 Recipients weight 0). Several rules are joined by ``+``."""
    if not rule:
        return
    from biobox.synthworld import backbone

    rules = rule.split("+")
    if set(rules) - set(RULES):
        raise ValueError(rule)
    orig = backbone._declared_groups
    if getattr(orig, "f0210_rule", None):
        return

    def patched(p):
        out = orig(p)
        g = p["genomes"]
        for gr in out:
            if gr["kind"] not in ("planted", "profile_only"):
                continue
            if "modbg" in rules:
                gr["transfer"], gr["loss"] = float(g["transfer"]), float(g["loss"])
            if "modbal" in rules:
                gr["loss"] = gr["transfer"] + float(g["duplication"])
            if "lossT" in rules:
                gr["loss"] = gr["transfer"]
            if "dup0" in rules:
                gr["duplication"] = 0.0
        return out

    patched.f0210_rule = rule
    backbone._declared_groups = patched
    if "excl" in rules:
        import zombi2.genomes as zg
        from zombi2.params import Recipients

        fam = zg.family

        def family(name, *a, transfer_to=None, **k):
            if transfer_to is None and name[:3] in ("mod", "pro") and "." in name:
                transfer_to = Recipients().weighted_by(f"genomes:{name.split('.')[0]}.0",
                                                       {"present": 0.0, "absent": 1.0})
            return fam(name, *a, transfer_to=transfer_to, **k)

        zg.family = family


def _patch_crown(crown: float | None) -> None:
    """Rescale every species tree ZOMBI2 draws so that its crown age (first speciation to the present) is
    ``crown``: every time times one factor, before the burst is grafted and the stem set (a measurement device)."""
    if crown is None:
        return
    import types

    import zombi2.species as zs
    from zombi2.tree import rescale

    orig = zs.simulate_species_tree
    if getattr(orig, "f0210_crown", None):
        return

    def patched(*a, **k):
        tree = orig(*a, **k).complete_tree
        present = max(tree.nodes[i].end_time for i in tree.extant_leaves())
        return types.SimpleNamespace(complete_tree=rescale(tree, factor=crown / (present - tree.nodes[tree.root].end_time)))

    patched.f0210_crown = crown
    zs.simulate_species_tree = patched


def _tree(nwk: Path, extant: set[str]) -> dict:
    from biobox.synthworld import backbone

    nodes = backbone.parse_newick(nwk.read_text(encoding="utf-8"))
    depth = []
    for parent, _label, length in nodes:
        depth.append((depth[parent] if parent >= 0 else 0.0) + length)
    stem = nodes[0][2]
    present = max(depth[i] for i, (_, lab, _) in enumerate(nodes) if lab in extant)
    n_internal = len({p for p, _, _ in nodes if p >= 0})
    return {"crown_age": round(present - stem, 4), "tree_height": round(present, 4), "n_nodes": len(nodes),
            "n_speciations": n_internal - 1}


def measure(world: Path, params: dict) -> tuple[dict, list[dict]]:
    b = world / "backbone"
    genomes = _tsv(b / "genomes.tsv")
    gids = [g["genome_id"] for g in genomes]
    fams = {int(r["family"]): r for r in _tsv(b / "families.tsv")}
    genes = _tsv(b / "genes.tsv")
    reps = _tsv(b / "replicons.tsv")
    kind_of = {f: r["kind"] for f, r in fams.items()}
    per_genome = collections.Counter(g["genome_id"] for g in genes)
    by_kind = collections.Counter()
    for g in genes:
        k = kind_of[int(g["family"])]
        by_kind["vehicle" if k in ("vehicle_chassis", "cargo") else
                "background" if k in ("root", "originated") else k] += 1
    total = len(genes)
    chrom = [int(r["n_genes"]) for r in reps if r["kind"] == "chromosome"]
    pls = collections.Counter(r["name"] for r in reps if r["kind"] == "plasmid")
    pls_per = collections.Counter(r["genome_id"] for r in reps if r["kind"] == "plasmid")
    sizes = [per_genome[g] for g in gids]
    tree = _tree(world / "zombi2" / "species_complete.nwk", {g["species_node"] for g in genomes})
    meta = json.loads((b / "backbone.json").read_text(encoding="utf-8"))
    sc = meta["scramble"]
    units = _tsv(b / "truth_units.tsv")
    groups = {r["name"]: r for r in _tsv(b / "groups.tsv")}
    scored = [u for u in units if groups[u["group"]]["scored"] == "true" and int(u["n_members"]) >= 2]
    # modules: co-presence, intact units, copies, transfers and realized co-transfer
    present = collections.defaultdict(set)
    copies = collections.Counter()
    for g in genes:
        present[int(g["family"])].add(g["genome_id"])
        copies[(int(g["family"]), g["genome_id"])] += 1
    tr_rows = collections.Counter()
    ev_fams = collections.defaultdict(set)
    for e in _tsv(b / "events.tsv"):
        if e["kind"].startswith("transfer"):
            tr_rows[int(e["family"])] += 1
            ev_fams[e["event"]].add(int(e["family"]))
    m = params["modules"]
    together, mtr = m["together"], m["member_transfer"]
    grid = []
    n_ok = n_ok6 = 0
    intact_by = collections.defaultdict(set)
    for u in units:
        if u["intact"] == "true":
            intact_by[u["group"]].add(u["genome_id"])
    for name, gr in sorted(groups.items()):
        if gr["kind"] != "planted":
            continue
        i = int(name[3:])
        mem = [int(x) for x in gr["members"].split(";")]
        co = set(gids)
        for f in mem:
            co &= present.get(f, set())
        n_ok += len(co) >= MIN_GENOMES
        n_ok6 += len(co) >= 6
        cps = [copies[(f, g)] for f in mem for g in present.get(f, ())]
        touching = [fs for fs in ev_fams.values() if fs & set(mem)]
        whole = sum(set(mem) <= fs for fs in touching)
        grid.append({"module": name, "size": len(mem), "together": together[i % len(together)],
                     "member_transfer": mtr[(i // len(together)) % len(mtr)],
                     "genomes_all_members": len(co), "genomes_intact_unit": len(intact_by[name]),
                     "copies_per_present_member": round(sum(cps) / len(cps), 3) if cps else 0.0,
                     "transfers_per_member": round(sum(tr_rows[f] for f in mem) / len(mem), 2),
                     "transfer_events": len(touching),
                     "realized_together": round(whole / len(touching), 3) if touching else 0.0})
    row = {"seed": meta["seed"], **tree, "genomes": len(gids),
           "genes_mean": round(total / len(gids), 1), "genes_min": min(sizes), "genes_max": max(sizes),
           "chrom_mean": round(sum(chrom) / len(chrom), 1), "chrom_min": min(chrom), "chrom_max": max(chrom)}
    for k in ("planted", "profile_only", "decoy", "vehicle", "hub", "plasmid", "background"):
        row[f"share_{k}"] = round(by_kind[k] / total, 4)
        row[f"per_genome_{k}"] = round(by_kind[k] / len(gids), 1)
    row["share_planted_all"] = round((by_kind["planted"] + by_kind["profile_only"] + by_kind["decoy"]) / total, 4)
    for k, size in enumerate(params["plasmids"]["sizes"]):
        row[f"pls{k}_genomes"] = pls.get(f"pls{k}", 0)
    row["genomes_0_1_2_plasmids"] = "/".join(str(sum(1 for g in gids if pls_per[g] == n)) for n in (0, 1, 2))
    row["modules_all_members_ge9"] = n_ok
    row["modules_all_members_ge6"] = n_ok6
    row["module_genomes_min"] = min(r["genomes_all_members"] for r in grid) if grid else 0
    row["module_copies_mean"] = round(sum(r["copies_per_present_member"] for r in grid) / len(grid), 3) if grid else 0
    row["plus_fraction"] = sc["plus_fraction"]
    row["initial_adjacency_kept"] = sc["initial_adjacency_kept"]
    row["declared_clustered"] = sc["declared_clustered"]
    row["scored_units"] = len(scored)
    row["scored_units_intact"] = sum(u["intact"] == "true" for u in scored)
    breaks = []
    if not MEAN_BAND[0] <= row["chrom_mean"] <= MEAN_BAND[1]:
        breaks.append("chrom_mean")
    if row["chrom_min"] > MIN_MAX[0]:
        breaks.append("chrom_min>1100")
    if row["chrom_max"] < MIN_MAX[1]:
        breaks.append("chrom_max<3600")
    if row["genomes"] != int(params["species"]["n_extant"]):
        breaks.append("genomes")
    if any(row[f"pls{k}_genomes"] == 0 for k in range(len(params["plasmids"]["sizes"]))):
        breaks.append("plasmid_type_lost")
    if max(pls_per.values(), default=0) > 2:
        breaks.append("more_than_2_plasmids")
    if n_ok < len(grid):
        breaks.append(f"modules_below_{MIN_GENOMES}_genomes:{len(grid) - n_ok}")
    if not (0.4 <= sc["plus_fraction"] <= 0.6 and sc["initial_adjacency_kept"] <= 0.05
            and sc["declared_clustered"] <= 0.5) or row["scored_units_intact"] != row["scored_units"]:
        breaks.append("dec15")
    row["breaks"] = ";".join(breaks)
    return row, grid


def _one(args) -> tuple[dict, list[dict]]:
    out, label, seed, sets, rule, crown, keep = args
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    from biobox.synthworld import backbone

    _patch_rule(rule)
    _patch_crown(crown)
    params = backbone.load_preset("micro")
    for path, value in sets:
        _set(params, path, value)
    world = Path(out) / "worlds" / label / str(seed)
    if world.exists():
        shutil.rmtree(world)
    world.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    backbone.build(params, seed, world)
    dt = time.perf_counter() - t0
    row, grid = measure(world, params)
    row = {"label": label, "seconds": round(dt, 1), **row}
    for r in grid:
        r.update({"label": label, "seed": seed})
    if not keep:
        shutil.rmtree(world)
    print(json.dumps({k: row[k] for k in ("label", "seed", "seconds", "crown_age", "chrom_mean", "share_planted",
                                         "breaks")}), flush=True)
    return row, grid


def _seeds(spec: str) -> list[int]:
    from biobox.synthworld import backbone, nulls

    default = int(backbone.load_preset("micro")["default_seed"])
    out = []
    for s in spec.split(","):
        if s == "default":
            out.append(default)
        elif s == "n2":
            out += nulls.n2_seeds(default, 20)
        elif s.startswith("n2:"):
            out += [nulls.n2_seeds(default, 20)[int(k)] for k in s[3:].split("+")]
        else:
            out.append(int(s))
    return out


def scan(argv: list[str]) -> int:
    out, label = argv[0], argv[1]
    jobs, seed_spec, sets, rule, crown, keep = 4, "default,n2", [], None, None, False
    rest = argv[2:]
    while rest:
        a = rest.pop(0)
        if a == "--jobs":
            jobs = int(rest.pop(0))
        elif a == "--seeds":
            seed_spec = rest.pop(0)
        elif a == "--set":
            k, _, v = rest.pop(0).partition("=")
            sets.append((k, json.loads(v)))
        elif a == "--axes":            # e.g. origination:0.8/duplication+transfer+loss:0.8 ("none" for no axis)
            spec = rest.pop(0)
            axes = [] if spec == "none" else [{"rates": r.split("+"), "sigma": float(s)}
                                              for r, s in (x.split(":") for x in spec.split("/"))]
            sets.append(("lineage_rates.axes", axes))
        elif a == "--rule":
            rule = rest.pop(0)
        elif a == "--crown":
            crown = float(rest.pop(0))
        elif a == "--keep":
            keep = True
        else:
            raise SystemExit(USAGE)
    seeds = _seeds(seed_spec)
    Path(out).mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        results = list(ex.map(_one, [(out, label, s, sets, rule, crown, keep) for s in seeds]))
    rows = [r for r, _ in results]
    grid = [g for _, gs in results for g in gs]
    _write_tsv(Path(out) / f"{label}.tsv", rows)
    _write_tsv(Path(out) / f"{label}_grid.tsv", grid)
    (Path(out) / f"{label}.json").write_text(json.dumps({"sets": sets, "rule": rule, "crown": crown,
                                                         "seeds": seeds}) + "\n")
    return 0


def table(argv: list[str]) -> int:
    out = Path(argv[0])
    cols = ("label", "seed", "crown_age", "tree_height", "genomes", "genes_mean", "genes_min", "genes_max",
            "chrom_mean", "chrom_min", "chrom_max", "share_planted", "share_profile_only", "share_decoy",
            "share_vehicle", "share_hub", "share_plasmid", "share_background", "pls0_genomes", "pls1_genomes",
            "modules_all_members_ge9", "module_copies_mean", "plus_fraction", "initial_adjacency_kept",
            "declared_clustered", "breaks")
    print("| " + " | ".join(cols) + " |")
    print("|" + "---|" * len(cols))
    for label in argv[1:]:
        for r in _tsv(out / f"{label}.tsv"):
            print("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |")
    return 0


def main(argv: list[str]) -> int:
    if len(argv) < 3 or argv[0] not in ("scan", "table"):
        raise SystemExit(USAGE)
    return scan(argv[1:]) if argv[0] == "scan" else table(argv[1:])


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
