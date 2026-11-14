"""Card F-02 sessions 10 and 11: the `micro` backbone (S1-S2) over many seeds, measured from its emitted tables.

    python f0210_seed_scan.py scan OUT LABEL [--jobs N] [--seeds default,n2 | s1,s2,...] [--set a.b=JSON ...]
                                            [--axes R+R:S/R:S | none] [--crown AGE] [--rule R+R] [--n0] [--keep]
    python f0210_seed_scan.py table OUT LABEL [LABEL ...]
    python f0210_seed_scan.py crowns OUT_TSV [--seeds ...]

scan: builds the backbone of the `micro` preset (the generator of this tree) at each seed into
OUT/worlds/LABEL/<seed>, N seeds at once (one process each), measures it and writes one row per seed to OUT/LABEL.tsv
and the realized strength grid to OUT/LABEL_grid.tsv. Seeds: ``default`` is the preset's default seed, ``n2`` the 20
N2 seeds (`nulls.n2_seeds` of the default seed); integers are taken as they are. ``--set`` overrides an existing
preset leaf by its path (JSON value). ``--axes`` replaces `lineage_rates.axes` (``origination:1.2/duplication+
transfer+loss:0.8``; ``none`` for no axis). ``--n0`` also builds the N0-full backbone at the seed and compares its
background genome sizes with the signal's as `tests/r2/test_synthworld_nulls.py` does (geometric-mean ratio, Pearson r).

``--crown`` and ``--rule`` are session 10's measurement devices that patch ZOMBI2's calls in the worker process,
never the generator's code: ``--crown`` rescales each drawn species tree to that crown age measured from the first
speciation (`_patch_crown`), ``--rule`` changes the planted and profile-only groups' rates or recipients
(`_patch_rule`: modbg, modbal, lossT, dup0, excl). Since session 11 the generator has a fixed crown age measured from
the living lineages' MRCA (``species.crown_age``), entry exclusion (``modules.entry_exclusion``,
``profile_only.entry_exclusion``) and its own member duplication (``member_duplication``), so ``--crown`` is
meaningful only with ``--set species.crown_age=null`` and ``excl`` only with ``entry_exclusion`` false (the device
fills a family's recipient rule only where the generator leaves none); the session-10 labels are reproduced with
those three settings put back. The world directory is removed after it is measured unless ``--keep``.
table: prints the rows of the labels as a Markdown table.
crowns: the species trees alone (ZOMBI2's draw at each seed, before any rescaling), one row per seed with the crown
age from the first speciation and from the living lineages' most recent common ancestor (DEC-37 (3)).

Every number is counted here from the backbone's tables (`genes.tsv`, `families.tsv`, `groups.tsv`, `replicons.tsv`,
`truth_units.tsv`, `events.tsv`), ZOMBI2's `genome_events.tsv` (the entry-exclusion replay) and
`lineage_multipliers.tsv`, and the species tree (`zombi2/species_complete.nwk`); the DEC-15 (c) statistics are
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
    tips = [i for i, (_, lab, _) in enumerate(nodes) if lab in extant]
    _, mrca_depth, _ = backbone.extant_crown({i: (p if p >= 0 else None, depth[i]) for i, (p, _, _) in enumerate(nodes)},
                                             tips)
    return {"crown_age": round(present - stem, 4), "tree_height": round(present, 4), "n_nodes": len(nodes),
            "n_speciations": n_internal - 1, "crown_age_mrca": round(present - mrca_depth, 4)}


def _ranks(xs: list[float]) -> list[float]:
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2.0
        i = j + 1
    return r


def _pearson(a: list[float], b: list[float]) -> float:
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    sab = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    saa = sum((x - ma) ** 2 for x in a)
    sbb = sum((y - mb) ** 2 for y in b)
    return sab / (saa * sbb) ** 0.5 if saa > 0 and sbb > 0 else 0.0


def _exclusion_replay(zdir: Path, first_of: dict[int, int], members_of: dict[int, set[int]]) -> tuple[int, int, int]:
    """Replay ZOMBI2's `genome_events.tsv` (speciation rows included) for the module families: per event id, before
    any of its rows is applied, count the transfers of a module family whose recipient lineage holds a live copy of
    the module's first family, and those whose recipient holds a live copy of any of its families. Returns (module
    transfer rows, into a holder of the first family, into a holder of any member)."""
    live: dict[int, collections.Counter] = collections.defaultdict(collections.Counter)   # family -> lineage counts
    rows_ev: list[list[str]] = []
    n = into_first = into_any = 0

    def flush():
        nonlocal n, into_first, into_any
        for kind, fam, parents, kids in rows_ev:
            if kind.startswith("transfer") and fam in first_of:
                rec = kids[1].split("_")[0]
                n += 1
                into_first += live[first_of[fam]][rec] > 0
                into_any += any(live[f][rec] > 0 for f in members_of[fam])
        for kind, fam, parents, kids in rows_ev:
            if fam not in first_of:
                continue
            if parents:
                live[fam][parents.split("_")[0]] -= 1
            for k in kids:
                live[fam][k.split("_")[0]] += 1
        rows_ev.clear()

    last = None
    with open(zdir / "genome_events.tsv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["event"] != last:
                flush()
                last = r["event"]
            rows_ev.append([r["kind"], int(r["family"]), r["parents"], [c for c in r["children"].split(";") if c]])
    flush()
    return n, into_first, into_any


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
    # the core: root families present in at least 40 of the 46 genomes (test_hub_families_have_more_diverse_neighbours
    # needs 100 of them)
    row["root_families_ge40_genomes"] = sum(1 for r in fams.values()
                                            if r["kind"] == "root" and int(r["n_genomes_present"]) >= 40)
    row["scored_units"] = len(scored)
    row["scored_units_intact"] = sum(u["intact"] == "true" for u in scored)
    # session 11 columns: the crown from the living MRCA, the realized transfer axis per level, entry exclusion
    # replayed from ZOMBI2's events, the size spread and its link to the origination axis
    row["crown_age_mrca"] = tree["crown_age_mrca"]
    row["rescale_factor"] = round(float(meta.get("crown", {}).get("rescale_factor", 1.0)), 6)
    levels = sorted(set(mtr))
    med = []
    for k, lv in enumerate(levels):
        v = sorted(r["transfers_per_member"] for r in grid if r["member_transfer"] == lv)
        m_ = (v[(len(v) - 1) // 2] + v[len(v) // 2]) / 2 if v else 0.0
        med.append(m_)
        row[f"tpm_level{k}_median"] = round(m_, 2)
        row[f"tpm_level{k}_range"] = f"{v[0]}-{v[-1]}" if v else ""
        row[f"cells_level{k}_unpowered"] = sum(
            1 for t in together if not any(r["genomes_all_members"] >= MIN_GENOMES for r in grid
                                           if r["member_transfer"] == lv and r["together"] == t))
    row["levels_distinct"] = all(a < b for a, b in zip(med, med[1:]))
    names = {r["name"]: int(f) for f, r in fams.items() if r["name"]}
    first_of, members_of = {}, {}
    for name, gr in groups.items():
        if gr["kind"] in ("planted", "profile_only"):
            mem = {int(x) for x in gr["members"].split(";")}
            for f in mem:
                first_of[f] = names[f"{name}.0"]
                members_of[f] = mem
    n_tr, into_first, into_any = _exclusion_replay(world / "zombi2", first_of, members_of)
    row["module_transfer_rows"] = n_tr
    row["module_transfers_into_first_holder"] = into_first
    row["module_transfers_into_any_holder"] = into_any
    mean_c = sum(chrom) / len(chrom)
    row["chrom_cv"] = round((sum((c - mean_c) ** 2 for c in chrom) / len(chrom)) ** 0.5 / mean_c, 4)
    lm = world / "zombi2" / "lineage_multipliers.tsv"
    rho = ""
    if lm.exists():
        mult = {r["lineage"]: r.get("origination", "") for r in _tsv(lm)} if lm.stat().st_size else {}
        node_of = {g["genome_id"]: g["species_node"] for g in genomes}
        size_of = {r["genome_id"]: int(r["n_genes"]) for r in reps if r["kind"] == "chromosome"}
        pairs = [(size_of[g], float(mult[node_of[g]])) for g in gids if mult.get(node_of[g])]
        if len(pairs) == len(gids):
            rho = round(_pearson(_ranks([a for a, _ in pairs]), _ranks([b for _, b in pairs])), 4)
    row["spearman_chrom_vs_origination_multiplier"] = rho
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
    if not row["levels_distinct"]:
        breaks.append("transfer_levels_not_distinct")
    if params["modules"].get("entry_exclusion") and into_first:
        breaks.append(f"entry_exclusion_violated:{into_first}")
    row["breaks"] = ";".join(breaks)
    return row, grid


def _bg_sizes(world: Path) -> dict[str, int]:
    kind = {int(r["family"]): r["kind"] for r in _tsv(world / "backbone" / "families.tsv")}
    c = collections.Counter(r["genome_id"] for r in _tsv(world / "backbone" / "genes.tsv")
                            if kind[int(r["family"])] in ("root", "originated"))
    return dict(c)


def _one(args) -> tuple[dict, list[dict]]:
    out, label, seed, sets, rule, crown, keep, n0 = args
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    import math

    from biobox.synthworld import backbone, nulls

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
    if n0:       # N0-full's background sizes against the signal's (tests/r2/test_synthworld_nulls.py's comparison)
        w0 = world.with_name(world.name + "_n0full")
        if w0.exists():
            shutil.rmtree(w0)
        backbone.build(nulls.null_preset(params, "N0-full"), seed, w0)
        a, b = _bg_sizes(world), _bg_sizes(w0)
        g = sorted(a)
        ratio = math.exp(sum(math.log(b[x] / a[x]) for x in g) / len(g))
        row["n0_bg_ratio"] = round(ratio, 4)
        row["n0_bg_pearson"] = round(_pearson([a[x] for x in g], [b[x] for x in g]), 4)
        if not keep:
            shutil.rmtree(w0)
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
    jobs, seed_spec, sets, rule, crown, keep, n0 = 4, "default,n2", [], None, None, False, False
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
        elif a == "--n0":
            n0 = True
        else:
            raise SystemExit(USAGE)
    seeds = _seeds(seed_spec)
    Path(out).mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        results = list(ex.map(_one, [(out, label, s, sets, rule, crown, keep, n0) for s in seeds]))
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


def crowns(argv: list[str]) -> int:
    """The species tree of each seed as ZOMBI2 draws it for the backbone (`backbone._run_zombi2`'s call, before any
    rescaling or burst): crown age from the first speciation and from the living lineages' MRCA."""
    import zombi2.species as zs

    from biobox.synthworld import backbone, seeds

    out, seed_spec = argv[0], "default,n2"
    if argv[1:3] and argv[1] == "--seeds":
        seed_spec = argv[2]
    sp = backbone.load_preset("micro")["species"]
    burst = int(sp["burst_size"])
    rows = []
    for s in _seeds(seed_spec):
        tree = zs.simulate_species_tree(float(sp["birth"]), float(sp["death"]),
                                        n_extant=int(sp["n_extant"]) - (burst - 1) if burst > 1 else int(sp["n_extant"]),
                                        seed=seeds.tool_seed(s, "S1/species_tree")).complete_tree
        mrca, t_mrca, present = backbone.extant_crown({i: (n.parent, n.end_time) for i, n in tree.nodes.items()},
                                                      tree.extant_leaves())
        root = tree.nodes[tree.root]
        kids = [c for c in root.children]
        alive = []
        for c in kids:
            stack, n_ext = [c], 0
            while stack:
                x = stack.pop()
                stack.extend(tree.nodes[x].children)
                n_ext += tree.nodes[x].fate == "extant"
            alive.append(n_ext)
        rows.append({"seed": s, "crown_first_speciation": round(present - root.end_time, 4),
                     "crown_mrca": round(present - t_mrca, 4), "mrca_is_root": mrca == tree.root,
                     "extant_under_root_children": "/".join(str(a) for a in alive)})
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    _write_tsv(Path(out), rows)
    v = sorted(r["crown_mrca"] for r in rows)
    print(json.dumps({"n": len(v), "crown_mrca_median": (v[(len(v) - 1) // 2] + v[len(v) // 2]) / 2,
                      "min": v[0], "max": v[-1]}))
    return 0


def main(argv: list[str]) -> int:
    if len(argv) >= 2 and argv[0] == "crowns":
        return crowns(argv[1:])
    if len(argv) < 3 or argv[0] not in ("scan", "table"):
        raise SystemExit(USAGE)
    return scan(argv[1:]) if argv[0] == "scan" else table(argv[1:])


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
