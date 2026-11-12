"""The synthworld backbone: species tree, gene order, strands, families, copies and groups, with truth
(planning/tests/synthworld_spec.md §2.1 S1-S2, §2.2, §5.1, §6.1; card F-02 session 1).

This is the only module that imports ZOMBI2 (pinned at 0.54.0), and it imports it inside `build`, so the
preset guards run without it. ZOMBI2 runs under DEC-15's three conditions:

(a) every adjacency-scored group has ``fates="group"`` (`check_dec15`); truth is in any case computed from
    the realized gene order, never from the declared groups;
(b) ``replacement`` is 0 in every world with groups (`check_dec15`; ZOMBI2 refuses it too);
(c) the root gene order and strands are scrambled by synthworld: the species tree's stem is replaced by a
    stem of `root_scramble.stem_length` on which only inversions and transpositions run (every other rate,
    the declared families' own rates included, is 0 there), and the result is checked
    (`scramble_stats`, `check_scrambled`) on every build.

Shape on disk (``build(params, seed, out)``; ``read(out)`` loads it back)::

    out/zombi2/              ZOMBI2's own files, from OrderedGenomesResult.write(compact=True), kept for audit
    out/backbone/genomes.tsv       genome_id species_node n_genes n_replicons
    out/backbone/genes.tsv         genome_id replicon topology position strand family copy_label group unit
    out/backbone/families.tsv      family family_id kind name group order_index n_extant_copies n_genomes_present
    out/backbone/groups.tsv        group kind members together adjacent fates scored
    out/backbone/truth_units.tsv   genome_id group unit replicon n_members contiguous co_oriented in_order intact strand
    out/backbone/truth_adjacency.tsv genome_id replicon left_position right_position left_family left_strand
                                   right_family right_strand same_unit
    out/backbone/backbone.json     format, preset (resolved), seed, derived tool seeds, versions, root node,
                                   burst tips, scramble statistics, sha256 of every file above, digest

Every table holds **extant genomes only**, sorted by (genome_id, replicon, position). ``genome_id`` is
``g000``... in species-node order; ``replicon`` is the chromosome's index in its genome (ZOMBI2's order);
positions are 0-based; strand is +1/-1; ``family`` is ZOMBI2's integer family id and ``family_id`` its
string label ``F000123``; ``copy_label`` is the extant gene-tree leaf ``<species node>_g<copy>``; in
``genes.tsv`` ``group`` and ``unit`` name the gene's unit and are empty and -1 for a gene in none (a family's
group, adjacent or not, is in ``families.tsv``). A unit is one copy of an adjacent group in one genome,
keyed by (genome_id, unit): ZOMBI2 keeps a unit id across speciations, so one id spans many genomes.
``truth_adjacency`` has one row per adjacent gene pair in the realized order, the wrap pair included on a
circular replicon of three or more genes; ``same_unit`` says both genes are in the same unit.
``digest`` is BLAKE3-128 over the sorted (path, sha256) list of every file but backbone.json.
"""

from __future__ import annotations

import copy
import csv
import hashlib
import json
import shutil
import typing
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

import numpy as np
import yaml
from blake3 import blake3

from biobox.synthworld import seeds

FORMAT = "synthworld-backbone/1"
ZOMBI2_VERSION = "0.54.0"
PRESET_DIR = Path(__file__).resolve().parent / "presets"

#: ZOMBI2 0.54.0 headers the adapter reads (spec §2.2 schema pins, corrected by spike §7 item 5)
PINNED_HEADERS = {
    "gene_order.tsv": "lineage\tchromosome\ttopology\tposition\tstrand\tfamily\tcopy",
    "groups.tsv": "group\tunit\tfamily\tgene",
    "gene_trees.tsv": "family\tcomplete\textant",
    "genome_events.tsv": "time\tkind\tfamily\tparents\tchildren\tchromosome\tstart\tlength\tdest_chromosome"
                         "\tdest_position\tevent",
    "rearrangement_events.tsv": "time\tkind\tlineage\tchromosome\tstart\tlength\tdest_chromosome"
                                "\tdest_position\tflipped\tcuts",
    "initial_genome.tsv": "chromosome\ttopology\tposition\tstrand\tfamily\tcopy",
}
#: scramble thresholds (DEC-15 (c)); derived: on a random circular order of N genes a given pair is
#: adjacent with probability about 2/N, strands are Binomial(N, 1/2), and a declared gene's outside
#: neighbour is declared with probability about (declared genes)/N, near 0.1 in `micro`
SCRAMBLE_MAX_KEPT = 0.05
SCRAMBLE_PLUS = (0.4, 0.6)
SCRAMBLE_MAX_CLUSTERED = 0.5


class Dec15Error(ValueError):
    """A preset or a world breaks one of DEC-15's three conditions."""


class SchemaError(ValueError):
    """A ZOMBI2 output file does not have the pinned header."""


class Genome(NamedTuple):
    genome_id: str
    species_node: str
    n_genes: int
    n_replicons: int


class Gene(NamedTuple):
    genome_id: str
    replicon: int
    topology: str
    position: int
    strand: int
    family: int
    copy_label: str
    group: str
    unit: int


class Family(NamedTuple):
    family: int
    family_id: str
    kind: str          # root | originated | planted | profile_only
    name: str
    group: str
    order_index: int   # position in its group, -1 if none
    n_extant_copies: int
    n_genomes_present: int


class Group(NamedTuple):
    name: str
    kind: str          # planted | profile_only
    members: tuple[int, ...]
    together: float
    adjacent: bool
    fates: str
    scored: bool


class Unit(NamedTuple):
    genome_id: str
    group: str
    unit: int
    replicon: int      # -1 when the members sit on more than one replicon
    n_members: int
    contiguous: bool
    co_oriented: bool
    in_order: bool
    intact: bool
    strand: int        # 0 unless co-oriented


class Adjacency(NamedTuple):
    genome_id: str
    replicon: int
    left_position: int
    right_position: int
    left_family: int
    left_strand: int
    right_family: int
    right_strand: int
    same_unit: bool


@dataclass
class Backbone:
    params: dict
    seed: int
    root_node: str
    burst_tips: list[str]
    genomes: list[Genome]
    genes: list[Gene]
    families: list[Family]
    groups: list[Group]
    units: list[Unit]
    adjacency: list[Adjacency]
    scramble: dict
    digest: str
    out: Path
    meta: dict = field(default_factory=dict, repr=False)

    def group(self, name: str) -> Group:
        for g in self.groups:
            if g.name == name:
                return g
        raise KeyError(name)


# --- presets ---------------------------------------------------------------------------------------------

def _resolve(node):
    if isinstance(node, dict):
        if "v" in node and "tag" in node:
            if node["tag"] not in ("asserted", "fitted", "derived"):
                raise ValueError(f"unknown parameter tag {node['tag']!r}")
            return node["v"]
        return {k: _resolve(v) for k, v in node.items()}
    return node


def _merge(base: dict, over: dict, where: str = "") -> dict:
    out = copy.deepcopy(base)
    for k, v in over.items():
        if k not in out:
            raise KeyError(f"preset has no parameter {where + k!r}")
        out[k] = _merge(out[k], v, f"{where}{k}.") if isinstance(v, dict) else copy.deepcopy(v)
    return out


def load_preset(name: str, overrides: dict | None = None) -> dict:
    """The resolved preset ``presets/<name>.yaml`` (values only), with ``overrides`` merged in."""
    raw = yaml.safe_load((PRESET_DIR / f"{name}.yaml").read_text(encoding="utf-8"))
    params = _resolve(raw)
    return _merge(params, overrides) if overrides else params


def params_digest(params: dict) -> str:
    return blake3(json.dumps(params, sort_keys=True, separators=(",", ":")).encode()).hexdigest(length=16)


# --- DEC-15 ----------------------------------------------------------------------------------------------

def _declared_groups(p: dict) -> list[dict]:
    """The groups a preset declares, in declaration order (module i, then profile-only j)."""
    out = []
    m = p["modules"]
    for i in range(int(m["count"])):
        out.append({"name": f"mod{i:02d}", "kind": "planted", "size": int(m["sizes"][i % len(m["sizes"])]),
                    "together": float(m["together"][i % len(m["together"])]),
                    "transfer": float(m["member_transfer"][(i // len(m["together"])) % len(m["member_transfer"])]),
                    "loss": float(m["member_loss"]), "adjacent": True, "fates": m["fates"],
                    "scored": bool(m["scored"])})
    q = p["profile_only"]
    for j in range(int(q["count"])):
        out.append({"name": f"pro{j:02d}", "kind": "profile_only", "size": int(q["size"]),
                    "together": float(q["together"]), "transfer": float(q["member_transfer"]),
                    "loss": float(q["member_loss"]), "adjacent": False, "fates": "family", "scored": False})
    return out


def check_dec15(p: dict) -> None:
    """Refuse a preset that breaks DEC-15 (planning/DECISIONS.md)."""
    groups = _declared_groups(p)
    for g in groups:
        if g["scored"] and not (g["adjacent"] and g["fates"] == "group"):
            raise Dec15Error(f"DEC-15 (a): group {g['name']} is scored for adjacency, so it needs "
                             f"adjacent=True and fates='group' (got adjacent={g['adjacent']}, "
                             f"fates={g['fates']!r}); set scored false to keep it with truth from the "
                             f"realized order only")
    if groups and float(p["genomes"]["replacement"]) > 0.0:
        raise Dec15Error(f"DEC-15 (b): replacement must be 0 in a world with groups, got "
                         f"{p['genomes']['replacement']}")
    rs = p["root_scramble"]
    if not (float(rs["stem_length"]) > 0 and float(rs["transposition_passes"]) > 0
            and float(rs["inversion_passes"]) >= 0):
        raise Dec15Error(f"DEC-15 (c): the root must be scrambled (stem_length > 0, transposition_passes "
                         f"> 0), got {rs}")


def scramble_stats(initial, root, *, declared: set[int], grouped_pairs: set[frozenset]) -> dict:
    """How scrambled a root layout is against ZOMBI2's initial one. ``initial`` and ``root`` are the
    (family, strand) runs of one circular replicon; ``grouped_pairs`` are within-unit neighbour pairs,
    which the scramble cannot (and must not) break."""
    def pairs(layout):
        fam = [f for f, _ in layout]
        n = len(fam)
        return [frozenset((fam[i], fam[(i + 1) % n])) for i in range(n)] if n >= 3 else []

    before = set(pairs(initial)) - grouped_pairs
    after_list = pairs(root)
    after = set(after_list)
    outside = [q for q in after_list if q & declared and q not in grouped_pairs]
    return {
        "n_genes": len(root),
        "plus_fraction": round(sum(s == 1 for _, s in root) / len(root), 6) if root else 0.0,
        "initial_adjacency_kept": round(len(before & after) / len(before), 6) if before else 0.0,
        "declared_clustered": round(sum(q <= declared for q in outside) / len(outside), 6) if outside else 0.0,
    }


def check_scrambled(stats: dict) -> None:
    lo, hi = SCRAMBLE_PLUS
    if not (lo <= stats["plus_fraction"] <= hi and stats["initial_adjacency_kept"] <= SCRAMBLE_MAX_KEPT
            and stats["declared_clustered"] <= SCRAMBLE_MAX_CLUSTERED):
        raise Dec15Error(f"DEC-15 (c): the root is not scrambled: {stats} (need plus_fraction in "
                         f"{SCRAMBLE_PLUS}, initial_adjacency_kept <= {SCRAMBLE_MAX_KEPT}, "
                         f"declared_clustered <= {SCRAMBLE_MAX_CLUSTERED})")


# --- ZOMBI2 adapter --------------------------------------------------------------------------------------

def check_schema(zdir: Path) -> None:
    zdir = Path(zdir)
    for name, want in PINNED_HEADERS.items():
        with open(zdir / name, encoding="utf-8") as fh:
            got = fh.readline().rstrip("\n")
        if got != want:
            raise SchemaError(f"{name}: header {got!r}, pinned {want!r} (ZOMBI2 {ZOMBI2_VERSION})")


def _graft_burst(tree, size: int, depth: float, rng, Node, Tree):
    """Graft a clonal burst of ``size`` tips onto one extant tip (spec §2.1 S1): the tip's branch ends
    ``depth`` x crown height before the present and a random binary subtree of short branches follows."""
    nodes = {i: Node(n.id, n.parent, n.birth_time, n.end_time, tuple(n.children), n.fate)
             for i, n in tree.nodes.items()}
    extant = sorted(i for i, n in nodes.items() if n.fate == "extant")
    tip = extant[int(rng.integers(len(extant)))]
    present = nodes[tip].end_time
    crown = present - nodes[tree.root].end_time
    h = min(depth * crown, 0.5 * (present - nodes[tip].birth_time))
    split = present - h
    times = sorted(float(t) for t in rng.uniform(split, present, size - 2))
    next_id = max(nodes) + 1

    def fork(parent: int, t: float) -> list[int]:
        nonlocal next_id
        a, b = next_id, next_id + 1
        next_id += 2
        nodes[a] = Node(a, parent, t)
        nodes[b] = Node(b, parent, t)
        nodes[parent].end_time = t
        nodes[parent].children = (a, b)
        nodes[parent].fate = "speciation"
        return [a, b]

    active = fork(tip, split)
    for t in times:
        active.extend(fork(active.pop(int(rng.integers(len(active)))), t))
    for i in active:
        nodes[i].end_time = present
        nodes[i].fate = "extant"
    return Tree(nodes, tree.root), sorted(active)


def _run_zombi2(p: dict, seed: int, zdir: Path) -> tuple[str, list[str], dict]:
    from zombi2 import species as zspecies
    from zombi2.genomes import family as zfamily, group as zgroup, simulate_genomes_ordered
    from zombi2.params import PerCopy, PerLineage
    from zombi2.tree import Node, Tree, with_stem

    sp, g, rs = p["species"], p["genomes"], p["root_scramble"]
    if int(g["chromosomes"]) != 1:
        raise NotImplementedError("more than one replicon (plasmids) is a later F-02 session")
    species_seed = seeds.tool_seed(seed, "S1/species_tree")
    genomes_seed = seeds.tool_seed(seed, "S2/genomes")
    burst = int(sp["burst_size"])
    tree = zspecies.simulate_species_tree(float(sp["birth"]), float(sp["death"]),
                                          n_extant=int(sp["n_extant"]) - (burst - 1) if burst > 1
                                          else int(sp["n_extant"]),
                                          seed=species_seed).complete_tree
    burst_ids: list[int] = []
    if burst > 1:
        tree, burst_ids = _graft_burst(tree, burst, float(sp["burst_depth"]),
                                       seeds.generator(seed, "S1/clonal_burst"), Node, Tree)
    tree = with_stem(tree, float(rs["stem_length"]), mode="set")
    t0, t1 = tree.nodes[tree.root].birth_time, tree.nodes[tree.root].end_time

    def after(rate: float, scope):
        """0 on the scramble stem, ``rate`` from the first speciation on."""
        rate = float(rate)
        return 0.0 if rate == 0.0 else scope(1.0).changing_at({t0: 0.0, t1: rate})

    def stem_and_after(passes: float, rate: float):
        return PerCopy(1.0).changing_at({t0: float(passes) / float(rs["stem_length"]), t1: float(rate)})

    groups = _declared_groups(p)
    families, zgroups = [], []
    for gr in groups:
        members = [f"{gr['name']}.{j}" for j in range(gr["size"])]
        families += [zfamily(m, transfer=after(gr["transfer"], PerCopy), loss=after(gr["loss"], PerCopy))
                     for m in members]
        zgroups.append(zgroup(gr["name"], members, together=gr["together"], adjacent=gr["adjacent"],
                              fates=gr["fates"]))
    res = simulate_genomes_ordered(
        tree,
        duplication=after(g["duplication"], PerCopy), transfer=after(g["transfer"], PerCopy),
        loss=after(g["loss"], PerCopy),
        origination=after(float(g["origination_per_root_family"]) * int(g["initial_families"]), PerLineage),
        inversion=stem_and_after(rs["inversion_passes"], g["inversion"]),
        transposition=stem_and_after(rs["transposition_passes"], g["transposition"]),
        translocation=after(g["translocation"], PerCopy),
        chromosomes=int(g["chromosomes"]), topology=g["topology"],
        duplication_extent=g["duplication_extent"], loss_extent=g["loss_extent"],
        transfer_extent=g["transfer_extent"], inversion_extent=g["inversion_extent"],
        transposition_extent=g["transposition_extent"], translocation_extent=g["translocation_extent"],
        inversion_probability=float(g["inversion_probability"]), transfer_to=g["transfer_to"],
        replacement=float(g["replacement"]), initial_families=int(g["initial_families"]),
        families=families or None, groups=zgroups or None, max_family_size=int(g["max_family_size"]),
        seed=genomes_seed, progress=False)
    res.write(zdir, compact=True)
    labels = tree.labels()
    import zombi2
    meta = {"species_seed": species_seed, "genomes_seed": genomes_seed, "zombi2": zombi2.__version__,
            "numpy": np.__version__, "stem": [t0, t1]}
    return labels[tree.root], [labels[i] for i in burst_ids], meta | {
        "extant": [labels[i] for i in sorted(tree.extant_leaves())]}


def _rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _units_truth(genes: list[Gene], groups: dict[str, Group]) -> list[Unit]:
    n_by = {}
    for x in genes:
        n_by[(x.genome_id, x.replicon)] = n_by.get((x.genome_id, x.replicon), 0) + 1
    members: dict[tuple[str, int], list[Gene]] = {}
    for x in genes:
        if x.unit >= 0:
            members.setdefault((x.genome_id, x.unit), []).append(x)
    order = {f: i for g in groups.values() for i, f in enumerate(g.members)}
    out = []
    for (gid, unit), ms in sorted(members.items()):
        grp = ms[0].group
        reps = {m.replicon for m in ms}
        strands = {m.strand for m in ms}
        one_rep = len(reps) == 1
        contiguous = in_order = False
        if one_rep:
            n = n_by[(gid, ms[0].replicon)]
            circ = ms[0].topology == "circular"
            pos = sorted(m.position for m in ms)
            k = len(pos)
            starts = [s for s in range(k)
                      if all((pos[(s + i) % k] - pos[s]) % n == i if circ else
                             (s == 0 and pos[i] - pos[0] == i) for i in range(k))]
            contiguous = bool(starts)
            if contiguous:
                at = {m.position: m.family for m in ms}
                fams = [at[pos[(starts[0] + i) % k]] for i in range(k)]
                if strands == {-1}:
                    fams.reverse()
                idx = [order[f] for f in fams]
                in_order = idx == sorted(idx) and len(set(idx)) == k
        co = len(strands) == 1
        out.append(Unit(gid, grp, unit, ms[0].replicon if one_rep else -1, len(ms), contiguous, co,
                        in_order, contiguous and co and in_order, next(iter(strands)) if co else 0))
    return out


def _adjacency_truth(genes: list[Gene]) -> list[Adjacency]:
    out = []
    i = 0
    while i < len(genes):
        j = i
        while j < len(genes) and (genes[j].genome_id, genes[j].replicon) == (genes[i].genome_id, genes[i].replicon):
            j += 1
        run = genes[i:j]
        n = len(run)
        last = n if run[0].topology == "circular" and n >= 3 else n - 1
        for k in range(last):
            a, b = run[k], run[(k + 1) % n]
            out.append(Adjacency(a.genome_id, a.replicon, a.position, b.position, a.family, a.strand,
                                 b.family, b.strand, a.unit >= 0 and a.unit == b.unit))
        i = j
    return out


def _tables(zdir: Path, p: dict, extant: list[str]) -> tuple:
    names = json.loads((zdir / "genome_summary.json").read_text(encoding="utf-8"))["families"]["names"]
    declared = _declared_groups(p)
    groups = {}
    family_meta = {}
    for gr in declared:
        ids = tuple(int(names[f"{gr['name']}.{j}"]) for j in range(gr["size"]))
        groups[gr["name"]] = Group(gr["name"], gr["kind"], ids, gr["together"], gr["adjacent"], gr["fates"],
                                   gr["scored"])
        for j, f in enumerate(ids):
            family_meta[f] = (gr["kind"], f"{gr['name']}.{j}", gr["name"], j)
    unit_of = {int(r["gene"]): (r["group"], int(r["unit"])) for r in _rows(zdir / "groups.tsv")}
    genome_of = {node: f"g{k:03d}" for k, node in enumerate(sorted(extant, key=lambda s: int(s[1:])))}
    genes: list[Gene] = []
    replicon_index: dict[str, dict[str, int]] = {}    # node -> ZOMBI2 chromosome id -> index in genome
    for r in _rows(zdir / "gene_order.tsv"):
        node = r["lineage"]
        if node not in genome_of:
            continue
        reps = replicon_index.setdefault(node, {})
        if r["chromosome"] not in reps:
            reps[r["chromosome"]] = len(reps)
        grp, unit = unit_of.get(int(r["copy"][1:]), ("", -1))
        genes.append(Gene(genome_of[node], reps[r["chromosome"]], r["topology"], int(r["position"]),
                          int(r["strand"]), int(r["family"]), f"{node}_{r['copy']}", grp, unit))
    genes.sort(key=lambda x: (x.genome_id, x.replicon, x.position))
    copies: dict[int, int] = {}
    present: dict[int, set[str]] = {}
    n_genes: dict[str, int] = {}
    for x in genes:
        n_genes[x.genome_id] = n_genes.get(x.genome_id, 0) + 1
        copies[x.family] = copies.get(x.family, 0) + 1
        present.setdefault(x.family, set()).add(x.genome_id)
    genomes = [Genome(gid, node, n_genes.get(gid, 0), len(replicon_index.get(node, {})))
               for node, gid in sorted(genome_of.items(), key=lambda kv: kv[1])]
    n_root = int(p["genomes"]["initial_families"])
    families = []
    for f in sorted(set(copies) | set(family_meta)):
        kind, name, grp, idx = family_meta.get(f, ("root" if f < n_root else "originated", "", "", -1))
        families.append(Family(f, f"F{f:06d}", kind, name, grp, idx, copies.get(f, 0), len(present.get(f, ()))))
    return genomes, genes, families, list(groups.values())


def _layouts(zdir: Path, root_node: str):
    init = [(int(r["family"]), int(r["strand"]))
            for r in sorted(_rows(zdir / "initial_genome.tsv"), key=lambda r: int(r["position"]))]
    root = [(int(r["family"]), int(r["strand"]))
            for r in sorted((r for r in _rows(zdir / "gene_order.tsv") if r["lineage"] == root_node),
                            key=lambda r: int(r["position"]))]
    return init, root


# --- writing and reading ---------------------------------------------------------------------------------

_TABLES = {
    "genomes.tsv": Genome, "genes.tsv": Gene, "families.tsv": Family, "groups.tsv": Group,
    "truth_units.tsv": Unit, "truth_adjacency.tsv": Adjacency,
}


def _cell(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, tuple):
        return ";".join(str(x) for x in v)
    return str(v)


def _write_tsv(path: Path, cls, rows) -> None:
    lines = ["\t".join(cls._fields)] + ["\t".join(_cell(v) for v in r) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _parser(cls):
    """One converter per field, from the NamedTuple's resolved type hints."""
    conv = []
    for name in cls._fields:
        typ = typing.get_type_hints(cls)[name]
        if typ is bool:
            conv.append(lambda s: s == "true")
        elif typ in (int, float, str):
            conv.append(typ)
        elif typing.get_origin(typ) is tuple:
            conv.append(lambda s: tuple(int(x) for x in s.split(";")) if s else ())
        else:
            raise TypeError(f"{cls.__name__}.{name}: no TSV reader for {typ!r}")
    return conv


def _parse(cls, row: dict[str, str], conv) -> tuple:
    return cls(*(c(row[name]) for name, c in zip(cls._fields, conv)))


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(params: dict | str, seed: int, out) -> Backbone:
    """Build the backbone of one world into ``out`` (which must not exist) and return it."""
    p = load_preset(params) if isinstance(params, str) else copy.deepcopy(params)
    check_dec15(p)
    seed = seeds.check_seed(seed)
    out = Path(out)
    if out.exists():
        raise FileExistsError(f"{out} exists; a backbone is written once")
    tmp = out.with_name(out.name + ".partial")
    if tmp.exists():
        shutil.rmtree(tmp)
    (tmp / "backbone").mkdir(parents=True)
    zdir = tmp / "zombi2"
    root_node, burst_tips, meta = _run_zombi2(p, seed, zdir)
    check_schema(zdir)
    genomes, genes, families, groups = _tables(zdir, p, meta.pop("extant"))
    units = _units_truth(genes, {g.name: g for g in groups if g.adjacent})
    adjacency = _adjacency_truth(genes)
    init, root = _layouts(zdir, root_node)
    grouped_pairs = {frozenset(q) for g in groups if g.adjacent for q in zip(g.members, g.members[1:])}
    stats = scramble_stats(init, root, declared={f for g in groups for f in g.members},
                           grouped_pairs=grouped_pairs)
    check_scrambled(stats)
    bdir = tmp / "backbone"
    for name, rows in (("genomes.tsv", genomes), ("genes.tsv", genes), ("families.tsv", families),
                       ("groups.tsv", groups), ("truth_units.tsv", units), ("truth_adjacency.tsv", adjacency)):
        _write_tsv(bdir / name, _TABLES[name], rows)
    files = {str(f.relative_to(tmp)): _sha(f) for f in sorted(tmp.rglob("*")) if f.is_file()}
    digest = blake3("".join(f"{k}\t{v}\n" for k, v in sorted(files.items())).encode()).hexdigest(length=16)
    doc = {"format": FORMAT, "preset": p, "params_digest": params_digest(p), "seed": seed,
           "tool_seeds": {"S1/species_tree": meta["species_seed"], "S2/genomes": meta["genomes_seed"]},
           "versions": {"zombi2": meta["zombi2"], "numpy": meta["numpy"]},
           "zombi2_call": "species.simulate_species_tree + genomes.simulate_genomes_ordered (serial engine); "
                          "write(compact=True)",
           "scramble_stem": meta["stem"], "root_node": root_node, "burst_tips": burst_tips,
           "scramble": stats, "files": files, "digest": digest}
    (bdir / "backbone.json").write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    tmp.rename(out)
    return read(out)


def read(out) -> Backbone:
    out = Path(out)
    bdir = out / "backbone"
    doc = json.loads((bdir / "backbone.json").read_text(encoding="utf-8"))
    if doc["format"] != FORMAT:
        raise ValueError(f"{out}: format {doc['format']!r}, expected {FORMAT!r}")
    t = {}
    for name, cls in _TABLES.items():
        conv = _parser(cls)
        t[name] = [_parse(cls, r, conv) for r in _rows(bdir / name)]
    return Backbone(params=doc["preset"], seed=doc["seed"], root_node=doc["root_node"],
                    burst_tips=doc["burst_tips"], genomes=t["genomes.tsv"], genes=t["genes.tsv"],
                    families=t["families.tsv"], groups=t["groups.tsv"], units=t["truth_units.tsv"],
                    adjacency=t["truth_adjacency.tsv"], scramble=doc["scramble"], digest=doc["digest"],
                    out=out, meta=doc)


__all__ = ["Backbone", "Dec15Error", "SchemaError", "build", "read", "load_preset", "check_dec15",
           "check_schema", "scramble_stats", "check_scrambled"]
