"""The synthworld backbone: species tree and synthetic taxonomy, gene order, strands, families, copies,
groups, plasmids, vehicles, hubs and decoys, with truth (planning/tests/synthworld_spec.md §2.1 S1-S2, §2.2,
§5.1, §6.1; card F-02 sessions 1 and 2).

This is the only module that imports ZOMBI2 (pinned at 0.54.0), and it imports it inside `build`, so the
preset guards run without it. ZOMBI2 runs under DEC-15's three conditions:

(a) every adjacency-scored group has ``fates="group"`` (`check_dec15`); truth is in any case computed from
    the realized gene order, never from the declared groups;
(b) ``replacement`` is 0 in every world with groups (`check_dec15`; ZOMBI2 refuses it too);
(c) the root gene order and strands are scrambled by synthworld: the species tree's stem is replaced by a
    stem of `root_scramble.stem_length` on which only inversions and transpositions run (every other rate,
    the declared families' own rates included, is 0 there), and the result is checked
    (`scramble_stats`, `check_scrambled`) on every build.

What the preset declares, beyond ZOMBI2's anonymous root families (session 2 adds the last four):

- ``modules``: planted adjacent groups ``mod00``..., members ``mod00.0``...; ``profile_only``: groups
  ``pro00``... with ``adjacent=False``;
- ``plasmids``: one adjacent group ``pls<i>`` per entry of ``plasmids.sizes`` (members ``pls<i>.<j>``),
  ``fates="group"``, ``together=1``, members carrying the per-plasmid transfer and loss rates divided by the
  plasmid's size, no duplication of their own, and a recipient rule that weighs a lineage already holding
  the plasmid by ``plasmids.present_weight`` (entry exclusion). ZOMBI2 runs one chromosome; synthworld
  **excises** every extant unit of a plasmid group into a circular replicon of its own and joins the
  chromosome's genes around it. (ZOMBI2 deals root genes round-robin over its chromosomes and lands every
  origination and arriving transfer on a uniformly chosen one, so its own replicons converge in size, and
  a chromosome loss can take the chromosome; replicon-specific rates do not exist there.)
- ``vehicles``: adjacent groups ``veh<i>`` = chassis ``veh<i>.0``, cargo ``veh<i>.1``..``veh<i>.<k>``,
  chassis ``veh<i>.<k+1>``, ``fates="group"``, ``together=1``, per-vehicle rates divided by the size;
- ``hubs``: single families ``hub<i>`` transferring at ``transfer_multiplier`` x the run's rate;
- ``decoys``: single families ``dec<i>`` whose transfers are weighted by the presence of module
  ``partner_modules[i]``'s first family (``mod<m>.0``) in the recipient (``absent_weight`` without it);
- ``genomes.transfer_to``: ``distance`` (ZOMBI2 ``Distance(decay=genomes.transfer_decay)``) or ``uniform``;
  a family with a recipient rule of its own (plasmids, decoys) uses that rule instead;
- ``lineage_rates.axes``: each axis one ZOMBI2 ``Random('lineages', Drift(LogNormal(0, sigma)))`` shared by
  the run rates it lists (duplication, transfer, loss, origination), so genome size varies among lineages.

Shape on disk (``build(params, seed, out)``; ``read(out)`` loads it back)::

    out/zombi2/              ZOMBI2's own files, from OrderedGenomesResult.write(compact=True), kept for audit
    out/backbone/genomes.tsv       genome_id species_node n_genes n_replicons
    out/backbone/genes.tsv         genome_id replicon topology position strand family copy_label group unit
    out/backbone/families.tsv      family family_id kind name group order_index n_extant_copies n_genomes_present
                                   partner
    out/backbone/groups.tsv        group kind members together adjacent fates scored
    out/backbone/truth_units.tsv   genome_id group unit replicon n_members contiguous co_oriented in_order intact strand
    out/backbone/truth_adjacency.tsv genome_id replicon left_position right_position left_family left_strand
                                   right_family right_strand same_unit
    out/backbone/replicons.tsv     genome_id replicon kind name unit topology n_genes
    out/backbone/taxonomy.tsv      genome_id species_node tax_lineage
    out/backbone/taxa.tsv          rank taxon node cut_time n_genomes parent
    out/backbone/events.tsv        event time kind family lineage recipient parent children chromosome start
                                   length dest_chromosome dest_position
    out/backbone/gene_tree_nodes.tsv family node parent kind lineage copy_label branch_length extant
    out/backbone/backbone.json     format, preset (resolved), seed, derived tool seeds, versions, root node,
                                   burst tips, scramble statistics, taxonomy cut times, sha256 of every file
                                   above, digest

The first seven tables hold **extant genomes only**, sorted by (genome_id, replicon, position).
``genome_id`` is ``g000``... in species-node order; ``replicon`` 0 is the chromosome and 1, 2, ... are the
genome's plasmids in the order their first gene sits on ZOMBI2's chromosome; positions are 0-based in
the replicon; strand is +1/-1; ``family`` is ZOMBI2's integer family id and ``family_id`` its string label
``F000123``; ``copy_label`` is the extant gene-tree leaf ``<species node>_g<copy>``; in ``genes.tsv``
``group`` and ``unit`` name the gene's unit and are empty and -1 for a gene in none (a family's group,
adjacent or not, is in ``families.tsv``). A unit is one copy of an adjacent group in one genome, keyed by
(genome_id, unit): ZOMBI2 keeps a unit id across speciations, so one id spans many genomes.
``families.kind`` is root, originated, planted, profile_only, plasmid, vehicle_chassis, cargo, hub or decoy;
``partner`` is a decoy's partner family and -1 otherwise. ``groups.kind`` is planted, profile_only,
plasmid or vehicle; only planted groups are scored. ``truth_adjacency`` has one row per adjacent gene pair
in the realized order (plasmids excised), the wrap pair included on a circular replicon of three or more
genes; ``same_unit`` says both genes are in the same unit. ``replicons.name`` is ``chromosome`` or the
plasmid's group, ``unit`` the plasmid's unit (-1 for the chromosome).

``taxonomy.tax_lineage`` is ``d__<domain>;p__Phylum01;c__Class03;...;s__Species12``: a genome's taxon at a
rank is the species-tree branch alive at that rank's cut, ``crown + cut x (present - crown)`` with
``crown`` the first speciation; taxa are numbered per rank in the order of the first genome holding them.
``taxa`` has one row per taxon: the crossing branch (``node``), the cut's absolute time, the number of
genomes and the taxon of the rank above (``parent``).

``events`` is ZOMBI2's ``genome_events.tsv`` typed, in its order, **all lineages** (extinct ones and the
scramble stem included) and **without speciation rows** (the gene trees hold them): ``lineage`` is the
acting lineage (the donor of a transfer), ``recipient`` the transfer's recipient or empty, ``parent`` the
copy that ends and ``children`` the copies that begin (``;``-separated), empty integer cells -1.
``gene_tree_nodes`` is every family's *complete* gene tree from ``gene_trees.tsv``, one row per node in
preorder (``node`` 0 is the root, ``parent`` -1): internal nodes have kind speciation, duplication or
transfer and the lineage of their label ``<event>_<lineage>``; leaves have kind extant (a copy in
``genes.tsv``) or lost (lost copy or extinct lineage) and their ``copy_label``; ``branch_length`` is the
length above the node (the root's reaches back to the family's origination).
``digest`` is BLAKE3-128 over the sorted (path, sha256) list of every file but backbone.json.

``build_phages(params, seed, out, panel_seed=...)`` writes a phage set's backbone in the same shape (card
F-02 session 5; `mobile`): one linear chromosome per phage, no groups, the scramble statistics recorded but
not enforced, and ``role`` and ``panel_seed`` added to backbone.json.

``species_branches(out)`` reads the complete species tree (extinct lineages and the scramble stem included)
from ``out/zombi2/species_complete.nwk``, so later stages need not read ZOMBI2's files (card F-02 session 3).
"""

from __future__ import annotations

import collections
import copy
import csv
import hashlib
import json
import re
import shutil
import typing
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

import numpy as np
import yaml
from blake3 import blake3

from biobox.synthworld import seeds

FORMAT = "synthworld-backbone/2"
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
#: neighbour is declared with probability about (declared genes)/N, near 0.2 in `micro`
SCRAMBLE_MAX_KEPT = 0.05
SCRAMBLE_PLUS = (0.4, 0.6)
SCRAMBLE_MAX_CLUSTERED = 0.5
#: the run rates a lineage axis may vary (ZOMBI2's per-lineage targets that the preset sets)
LINEAGE_RATES = ("duplication", "transfer", "loss", "origination")


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
    kind: str          # root | originated | planted | profile_only | plasmid | vehicle_chassis | cargo | hub | decoy
    name: str
    group: str
    order_index: int   # position in its group, -1 if none
    n_extant_copies: int
    n_genomes_present: int
    partner: int       # a decoy's partner family, -1 otherwise


class Group(NamedTuple):
    name: str
    kind: str          # planted | profile_only | plasmid | vehicle
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


class Replicon(NamedTuple):
    genome_id: str
    replicon: int
    kind: str          # chromosome | plasmid
    name: str          # "chromosome" or the plasmid's group
    unit: int          # the plasmid's unit, -1 for the chromosome
    topology: str
    n_genes: int


class TaxonomyRow(NamedTuple):
    genome_id: str
    species_node: str
    tax_lineage: str


class Taxon(NamedTuple):
    rank: str
    taxon: str
    node: str          # the species-tree branch alive at the rank's cut
    cut_time: float
    n_genomes: int
    parent: str        # the taxon of the rank above (the domain for the first rank)


class GeneEvent(NamedTuple):
    event: int
    time: float
    kind: str
    family: int
    lineage: str       # acting lineage; a transfer's donor
    recipient: str     # a transfer's recipient, else empty
    parent: str
    children: tuple[str, ...]
    chromosome: int
    start: int
    length: int
    dest_chromosome: int
    dest_position: int


class GeneTreeNode(NamedTuple):
    family: int
    node: int
    parent: int
    kind: str          # speciation | duplication | transfer | extant | lost
    lineage: str
    copy_label: str    # leaves only
    branch_length: float
    extant: bool


class SpeciesBranch(NamedTuple):
    lineage: str       # n<id> or e<id>, as in gene-tree labels
    parent: str        # empty for the root (the scramble stem)
    branch_length: float


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
    replicons: list[Replicon] = field(default_factory=list)
    taxonomy: list[TaxonomyRow] = field(default_factory=list)
    taxa: list[Taxon] = field(default_factory=list)
    events: list[GeneEvent] = field(default_factory=list, repr=False)
    gene_tree_nodes: list[GeneTreeNode] = field(default_factory=list, repr=False)
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


def check_preset(p: dict) -> None:
    """Refuse a malformed session-2 parameter (taxonomy, plasmids, vehicles, decoys, transfer rule,
    lineage axes) before anything is simulated."""
    tx = p["taxonomy"]
    cuts = [float(c) for c in tx["cuts"]]
    if len(cuts) != len(tx["ranks"]) or not all(0.0 < a < b < 1.0 for a, b in zip(cuts, cuts[1:])) \
            or not (cuts and 0.0 < cuts[0] and cuts[-1] < 1.0):
        raise ValueError(f"taxonomy cuts must be one per rank, strictly increasing in (0, 1): {tx}")
    for s in p["plasmids"]["sizes"]:
        if int(s) < 2:
            raise ValueError(f"a plasmid is a ZOMBI2 group and needs at least 2 genes, got {s}")
    for k in p["vehicles"]["cargo"]:
        if int(k) < 1:
            raise ValueError(f"a vehicle needs at least one cargo family, got {k}")
    n_mod = int(p["modules"]["count"])
    for m in p["decoys"]["partner_modules"]:
        if not 0 <= int(m) < n_mod:
            raise ValueError(f"decoy partner module {m} is not one of the {n_mod} planted modules")
    if p["genomes"]["transfer_to"] not in ("uniform", "distance"):
        raise ValueError(f"genomes.transfer_to is 'uniform' or 'distance', got {p['genomes']['transfer_to']!r}")
    seen: set[str] = set()
    for ax in p["lineage_rates"]["axes"]:
        if float(ax["sigma"]) < 0:
            raise ValueError(f"lineage_rates: sigma must be >= 0, got {ax}")
        for r in ax["rates"]:
            if r not in LINEAGE_RATES or r in seen:
                raise ValueError(f"lineage_rates: {r!r} is not one of {LINEAGE_RATES} or is in two axes")
            seen.add(r)


# --- DEC-15 ----------------------------------------------------------------------------------------------

def _declared_groups(p: dict) -> list[dict]:
    """The groups a preset declares, in declaration order: module i, profile-only j, plasmid k, vehicle l.
    ``duplication`` None means the run's rate."""
    out = []
    m = p["modules"]
    for i in range(int(m["count"])):
        out.append({"name": f"mod{i:02d}", "kind": "planted", "size": int(m["sizes"][i % len(m["sizes"])]),
                    "together": float(m["together"][i % len(m["together"])]),
                    "transfer": float(m["member_transfer"][(i // len(m["together"])) % len(m["member_transfer"])]),
                    "loss": float(m["member_loss"]), "duplication": None, "adjacent": True,
                    "fates": m["fates"], "scored": bool(m["scored"])})
    q = p["profile_only"]
    for j in range(int(q["count"])):
        out.append({"name": f"pro{j:02d}", "kind": "profile_only", "size": int(q["size"]),
                    "together": float(q["together"]), "transfer": float(q["member_transfer"]),
                    "loss": float(q["member_loss"]), "duplication": None, "adjacent": False,
                    "fates": "family", "scored": False})
    pl = p["plasmids"]
    for k, size in enumerate(int(s) for s in pl["sizes"]):
        out.append({"name": f"pls{k}", "kind": "plasmid", "size": size, "together": 1.0,
                    "transfer": float(pl["transfer"]) / size, "loss": float(pl["loss"]) / size,
                    "duplication": 0.0, "adjacent": True, "fates": "group", "scored": False})
    ve = p["vehicles"]
    for v, cargo in enumerate(int(c) for c in ve["cargo"]):
        size = cargo + 2
        out.append({"name": f"veh{v}", "kind": "vehicle", "size": size, "together": 1.0,
                    "transfer": float(ve["transfer"]) / size, "loss": float(ve["loss"]) / size,
                    "duplication": 0.0, "adjacent": True, "fates": "group", "scored": False})
    return out


def _declared_singles(p: dict) -> list[dict]:
    """The single declared families: hubs, then decoys (with their partner's name)."""
    g, h, d = p["genomes"], p["hubs"], p["decoys"]
    out = [{"name": f"hub{i}", "kind": "hub", "transfer": float(g["transfer"]) * float(h["transfer_multiplier"]),
            "loss": float(h["loss"]), "partner": ""} for i in range(int(h["count"]))]
    out += [{"name": f"dec{i}", "kind": "decoy", "transfer": float(d["transfer"]), "loss": float(d["loss"]),
             "partner": f"mod{int(m):02d}.0"} for i, m in enumerate(d["partner_modules"])]
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


# --- taxonomy (§2.1 S1) and Newick -------------------------------------------------------------------

def cut_taxonomy(nodes: dict, tips, *, crown: float, present: float, cuts) -> dict[str, list[str]]:
    """For each tip, the branch alive at each cut: ``nodes`` maps a label to (parent label or None, start
    time, end time); a cut ``f`` is the time ``crown + f x (present - crown)``; the branch alive at time
    ``c`` is the ancestor (or the tip) with start < c <= end."""
    times = [crown + float(f) * (present - crown) for f in cuts]
    out = {}
    for tip in tips:
        chain = []
        x = tip
        while x is not None:
            chain.append(x)
            x = nodes[x][0]
        out[tip] = [next(y for y in chain if nodes[y][1] < c <= nodes[y][2]) for c in times]
    return out


def name_taxa(per_tip: dict, *, domain: str, ranks, cut_times) -> tuple[dict[str, str], list[Taxon]]:
    """Name the taxa of `cut_taxonomy`'s result: per rank, numbered in the order of the first key (genome)
    holding them, ``<r>__<Rank>NN`` with ``<r>`` the rank's first letter. Returns the lineage string per
    key and the taxa, rank by rank in number order."""
    keys = list(per_tip)
    names: list[dict[str, str]] = [{} for _ in ranks]
    for k in keys:
        for r, node in enumerate(per_tip[k]):
            if node not in names[r]:
                names[r][node] = f"{ranks[r][0]}__{ranks[r].capitalize()}{len(names[r]) + 1:02d}"
    top = f"d__{domain}"
    lineages = {k: ";".join([top] + [names[r][n] for r, n in enumerate(per_tip[k])]) for k in keys}
    taxa = []
    for r, rank in enumerate(ranks):
        for node, name in names[r].items():
            members = [k for k in keys if per_tip[k][r] == node]
            parent = top if r == 0 else names[r - 1][per_tip[members[0]][r - 1]]
            taxa.append(Taxon(rank, name, node, float(cut_times[r]), len(members), parent))
    return lineages, taxa


_NEWICK = re.compile(r"[(),;]|[^(),;]+")


def parse_newick(text: str) -> list[tuple[int, str, float]]:
    """``(parent index, label, branch length)`` per node, in preorder (the order nodes open in the text),
    the root first with parent -1."""
    nodes: list[list] = []
    stack: list[int] = []
    last, closed = -1, False
    for tok in _NEWICK.findall(text.strip()):
        if tok == "(":
            nodes.append([stack[-1] if stack else -1, "", 0.0])
            stack.append(len(nodes) - 1)
            closed = False
        elif tok == ")":
            last = stack.pop()
            closed = True
        elif tok in (",", ";"):
            closed = False
        else:
            label, _, length = tok.partition(":")
            value = float(length) if length else 0.0
            if closed:
                nodes[last][1], nodes[last][2] = label, value
            else:
                nodes.append([stack[-1] if stack else -1, label, value])
            closed = False
    return [(a, b, c) for a, b, c in nodes]


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
    from zombi2.params import Distance, Drift, LogNormal, PerCopy, PerLineage, Random, Recipients
    from zombi2.tree import Node, Tree, with_stem

    sp, g, rs = p["species"], p["genomes"], p["root_scramble"]
    if int(g["chromosomes"]) != 1:
        raise NotImplementedError("synthworld runs one ZOMBI2 chromosome; plasmids are excised groups")
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

    lineage_axis = {}
    for ax in p["lineage_rates"]["axes"]:
        if float(ax["sigma"]) > 0:
            obj = Random("lineages", Drift(LogNormal(0.0, float(ax["sigma"]))))
            lineage_axis.update({r: obj for r in ax["rates"]})

    def after(rate: float, scope, axis: str | None = None):
        """0 on the scramble stem, ``rate`` from the first speciation on; a run rate named in a lineage
        axis also varies among species branches."""
        rate = float(rate)
        if rate == 0.0:
            return 0.0
        r = scope(1.0).changing_at({t0: 0.0, t1: rate})
        return r.varying_among(lineage_axis[axis]) if axis in lineage_axis else r

    def stem_and_after(passes: float, rate: float):
        return PerCopy(1.0).changing_at({t0: float(passes) / float(rs["stem_length"]), t1: float(rate)})

    groups = _declared_groups(p)
    families, zgroups = [], []
    joint = False
    for gr in groups:
        members = [f"{gr['name']}.{j}" for j in range(gr["size"])]
        rule = None
        if gr["kind"] == "plasmid" and float(p["plasmids"]["present_weight"]) != 1.0:
            rule = Recipients().weighted_by(f"genomes:{members[0]}",
                                            {"present": float(p["plasmids"]["present_weight"]), "absent": 1.0})
            joint = True
        own = {"transfer": after(gr["transfer"], PerCopy), "loss": after(gr["loss"], PerCopy)}
        if gr["duplication"] is not None:
            own["duplication"] = after(gr["duplication"], PerCopy)
        families += [zfamily(m, transfer_to=rule, **own) for m in members]
        zgroups.append(zgroup(gr["name"], members, together=gr["together"], adjacent=gr["adjacent"],
                              fates=gr["fates"]))
    aw = float(p["decoys"]["absent_weight"])
    for s in _declared_singles(p):
        rule = None
        if s["partner"]:
            rule = Recipients().weighted_by(f"genomes:{s['partner']}", {"present": 1.0, "absent": aw})
            joint = True
        families.append(zfamily(s["name"], transfer=after(s["transfer"], PerCopy),
                                loss=after(s["loss"], PerCopy), transfer_to=rule))
    transfer_to = Distance(decay=float(g["transfer_decay"])) if g["transfer_to"] == "distance" else "uniform"
    res = simulate_genomes_ordered(
        tree,
        duplication=after(g["duplication"], PerCopy, "duplication"),
        transfer=after(g["transfer"], PerCopy, "transfer"),
        loss=after(g["loss"], PerCopy, "loss"),
        origination=after(float(g["origination_per_root_family"]) * int(g["initial_families"]), PerLineage,
                          "origination"),
        inversion=stem_and_after(rs["inversion_passes"], g["inversion"]),
        transposition=stem_and_after(rs["transposition_passes"], g["transposition"]),
        translocation=after(g["translocation"], PerCopy),
        chromosomes=int(g["chromosomes"]), topology=g["topology"],
        duplication_extent=g["duplication_extent"], loss_extent=g["loss_extent"],
        transfer_extent=g["transfer_extent"], inversion_extent=g["inversion_extent"],
        transposition_extent=g["transposition_extent"], translocation_extent=g["translocation_extent"],
        inversion_probability=float(g["inversion_probability"]), transfer_to=transfer_to,
        replacement=float(g["replacement"]), initial_families=int(g["initial_families"]),
        families=families or None, groups=zgroups or None, joint=joint,
        max_family_size=int(g["max_family_size"]), seed=genomes_seed, progress=False)
    res.write(zdir, compact=True)
    labels = tree.labels()
    nodes = {labels[i]: (labels[n.parent] if n.parent is not None else None, n.birth_time, n.end_time)
             for i, n in tree.nodes.items()}
    extant = [labels[i] for i in sorted(tree.extant_leaves())]
    present = max(tree.nodes[i].end_time for i in tree.extant_leaves())
    cuts = [float(c) for c in p["taxonomy"]["cuts"]]
    import zombi2
    meta = {"species_seed": species_seed, "genomes_seed": genomes_seed, "zombi2": zombi2.__version__,
            "numpy": np.__version__, "stem": [t0, t1], "joint": joint,
            "taxonomy": {"crown": t1, "present": present, "cut_times": [t1 + f * (present - t1) for f in cuts],
                         "per_tip": cut_taxonomy(nodes, extant, crown=t1, present=present, cuts=cuts)}}
    return labels[tree.root], [labels[i] for i in burst_ids], meta | {"extant": extant}


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
    family_meta: dict[int, tuple] = {}
    for gr in declared:
        ids = tuple(int(names[f"{gr['name']}.{j}"]) for j in range(gr["size"]))
        groups[gr["name"]] = Group(gr["name"], gr["kind"], ids, gr["together"], gr["adjacent"], gr["fates"],
                                   gr["scored"])
        for j, f in enumerate(ids):
            kind = gr["kind"]
            if kind == "vehicle":
                kind = "vehicle_chassis" if j in (0, len(ids) - 1) else "cargo"
            family_meta[f] = (kind, f"{gr['name']}.{j}", gr["name"], j, -1)
    for s in _declared_singles(p):
        family_meta[int(names[s["name"]])] = (s["kind"], s["name"], "", -1,
                                              int(names[s["partner"]]) if s["partner"] else -1)
    plasmid_groups = {gr["name"] for gr in declared if gr["kind"] == "plasmid"}
    unit_of = {int(r["gene"]): (r["group"], int(r["unit"])) for r in _rows(zdir / "groups.tsv")}
    genome_of = {node: f"g{k:03d}" for k, node in enumerate(sorted(extant, key=lambda s: int(s[1:])))}
    rows_of: dict[str, list[dict]] = collections.defaultdict(list)
    for r in _rows(zdir / "gene_order.tsv"):
        if r["lineage"] in genome_of:
            rows_of[r["lineage"]].append(r)
    genes: list[Gene] = []
    replicons: list[Replicon] = []
    for node, gid in sorted(genome_of.items(), key=lambda kv: kv[1]):
        rs = sorted(rows_of.get(node, []), key=lambda r: int(r["position"]))
        if len({r["chromosome"] for r in rs}) > 1:
            raise ValueError(f"{node}: more than one ZOMBI2 chromosome; synthworld runs one")
        chrom: list[tuple] = []
        plasmids: dict[tuple[str, int], list[tuple]] = {}
        for r in rs:
            grp, unit = unit_of.get(int(r["copy"][1:]), ("", -1))
            gene = (int(r["strand"]), int(r["family"]), f"{node}_{r['copy']}", grp, unit)
            (plasmids.setdefault((grp, unit), []) if grp in plasmid_groups else chrom).append(gene)
        topology = rs[0]["topology"] if rs else "circular"
        genes += [Gene(gid, 0, topology, i, *x) for i, x in enumerate(chrom)]
        replicons.append(Replicon(gid, 0, "chromosome", "chromosome", -1, topology, len(chrom)))
        for k, ((grp, unit), xs) in enumerate(plasmids.items(), start=1):
            genes += [Gene(gid, k, "circular", i, *x) for i, x in enumerate(xs)]
            replicons.append(Replicon(gid, k, "plasmid", grp, unit, "circular", len(xs)))
    copies: dict[int, int] = {}
    present: dict[int, set[str]] = {}
    n_genes: dict[str, int] = {}
    n_reps: dict[str, int] = {}
    for x in genes:
        n_genes[x.genome_id] = n_genes.get(x.genome_id, 0) + 1
        copies[x.family] = copies.get(x.family, 0) + 1
        present.setdefault(x.family, set()).add(x.genome_id)
    for r in replicons:
        n_reps[r.genome_id] = n_reps.get(r.genome_id, 0) + 1
    genomes = [Genome(gid, node, n_genes.get(gid, 0), n_reps.get(gid, 0))
               for node, gid in sorted(genome_of.items(), key=lambda kv: kv[1])]
    n_root = int(p["genomes"]["initial_families"])
    families = []
    for f in sorted(set(copies) | set(family_meta)):
        kind, name, grp, idx, partner = family_meta.get(f, ("root" if f < n_root else "originated", "", "", -1, -1))
        families.append(Family(f, f"F{f:06d}", kind, name, grp, idx, copies.get(f, 0), len(present.get(f, ())),
                               partner))
    return genomes, genes, families, list(groups.values()), replicons, genome_of


def _events(zdir: Path) -> list[GeneEvent]:
    def num(s: str) -> int:
        return int(s) if s else -1

    out = []
    with open(zdir / "genome_events.tsv", newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        for time, kind, family, parents, children, chrom, start, length, dchrom, dpos, event in reader:
            if kind == "speciation":
                continue
            kids = tuple(c for c in children.split(";") if c)
            out.append(GeneEvent(int(event), float(time), kind, int(family), (parents or kids[0]).split("_")[0],
                                 kids[1].split("_")[0] if kind.startswith("transfer") else "", parents, kids,
                                 num(chrom), num(start), num(length), num(dchrom), num(dpos)))
    return out


def _gene_tree_nodes(zdir: Path, extant_copies: set[str]) -> list[GeneTreeNode]:
    out = []
    for r in _rows(zdir / "gene_trees.tsv"):
        fam = int(r["family"])
        for idx, (parent, label, length) in enumerate(parse_newick(r["complete"])):
            a, _, b = label.partition("_")
            if b.startswith("g"):              # a leaf, "<lineage>_g<copy>"
                ext = label in extant_copies
                out.append(GeneTreeNode(fam, idx, parent, "extant" if ext else "lost", a, label, length, ext))
            else:                              # "<event>_<lineage>"
                out.append(GeneTreeNode(fam, idx, parent, a, b, "", length, False))
    return out


def species_branches(out) -> list[SpeciesBranch]:
    """The complete species tree of a built world (extinct lineages and the scramble stem included), one
    branch per lineage in preorder, the root first: read from ZOMBI2's ``species_complete.nwk`` here so that
    later stages (sequences) need not know ZOMBI2's files (spec §2.2)."""
    nodes = parse_newick((Path(out) / "zombi2" / "species_complete.nwk").read_text(encoding="utf-8"))
    return [SpeciesBranch(label, nodes[parent][1] if parent >= 0 else "", length)
            for parent, label, length in nodes]


def _taxonomy_tables(meta: dict, p: dict, genome_of: dict[str, str]) -> tuple[list[TaxonomyRow], list[Taxon]]:
    tx = meta["taxonomy"]
    by_genome = {genome_of[node]: tx["per_tip"][node] for node in sorted(genome_of, key=genome_of.get)}
    lineages, taxa = name_taxa(by_genome, domain=p["taxonomy"]["domain"], ranks=p["taxonomy"]["ranks"],
                               cut_times=tx["cut_times"])
    node_of = {gid: node for node, gid in genome_of.items()}
    return [TaxonomyRow(gid, node_of[gid], lin) for gid, lin in lineages.items()], taxa


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
    "truth_units.tsv": Unit, "truth_adjacency.tsv": Adjacency, "replicons.tsv": Replicon,
    "taxonomy.tsv": TaxonomyRow, "taxa.tsv": Taxon, "events.tsv": GeneEvent,
    "gene_tree_nodes.tsv": GeneTreeNode,
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
            item = typing.get_args(typ)[0]
            conv.append(lambda s, item=item: tuple(item(x) for x in s.split(";")) if s else ())
        else:
            raise TypeError(f"{cls.__name__}.{name}: no TSV reader for {typ!r}")
    return conv


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(params: dict | str, seed: int, out) -> Backbone:
    """Build the backbone of one world into ``out`` (which must not exist) and return it."""
    p = load_preset(params) if isinstance(params, str) else copy.deepcopy(params)
    check_dec15(p)
    check_preset(p)
    return _build(p, seed, out, None)


def build_phages(params: dict, seed: int, out, *, panel_seed: int) -> Backbone:
    """The backbone of a world's phage set (spec §2.1 S6; card F-02 session 5): the same ZOMBI2 ordered run
    and adapter tables as `build`, from a preset that `mobile.phage_preset` derives from the host's (linear
    topology, no groups and no declared families, no burst, no lineage axes). DEC-15 governs worlds with
    groups, and a phage set has none, so the root scramble is not enforced: its stem only spreads the root's
    genes over both strands (ZOMBI2 lays the root out on +). ``backbone.json`` also records ``role``
    ``phage_set`` and ``panel_seed``, the host world's seed, from which `sequences.build` draws the HMM
    panel, so phage families take their domains from the host's panel."""
    p = copy.deepcopy(params)
    check_preset(p)
    if _declared_groups(p) or _declared_singles(p):
        raise ValueError("a phage set declares no groups and no single families")
    if p["genomes"]["topology"] != "linear" or int(p["genomes"]["chromosomes"]) != 1:
        raise ValueError("a phage set is one linear chromosome per phage")
    return _build(p, seed, out, {"role": "phage_set", "panel_seed": seeds.check_seed(panel_seed)})


def _build(p: dict, seed: int, out, extra: dict | None) -> Backbone:
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
    genomes, genes, families, groups, replicons, genome_of = _tables(zdir, p, meta.pop("extant"))
    units = _units_truth(genes, {g.name: g for g in groups if g.adjacent})
    adjacency = _adjacency_truth(genes)
    taxonomy, taxa = _taxonomy_tables(meta, p, genome_of)
    events = _events(zdir)
    tree_nodes = _gene_tree_nodes(zdir, {x.copy_label for x in genes})
    init, root = _layouts(zdir, root_node)
    grouped_pairs = {frozenset(q) for g in groups if g.adjacent for q in zip(g.members, g.members[1:])}
    stats = scramble_stats(init, root, declared={f.family for f in families if f.kind not in ("root", "originated")},
                           grouped_pairs=grouped_pairs)
    if extra is None:
        check_scrambled(stats)
    bdir = tmp / "backbone"
    for name, rows in (("genomes.tsv", genomes), ("genes.tsv", genes), ("families.tsv", families),
                       ("groups.tsv", groups), ("truth_units.tsv", units), ("truth_adjacency.tsv", adjacency),
                       ("replicons.tsv", replicons), ("taxonomy.tsv", taxonomy), ("taxa.tsv", taxa),
                       ("events.tsv", events), ("gene_tree_nodes.tsv", tree_nodes)):
        _write_tsv(bdir / name, _TABLES[name], rows)
    files = {str(f.relative_to(tmp)): _sha(f) for f in sorted(tmp.rglob("*")) if f.is_file()}
    digest = blake3("".join(f"{k}\t{v}\n" for k, v in sorted(files.items())).encode()).hexdigest(length=16)
    tx = meta["taxonomy"]
    doc = {"format": FORMAT, "preset": p, "params_digest": params_digest(p), "seed": seed,
           "tool_seeds": {"S1/species_tree": meta["species_seed"], "S2/genomes": meta["genomes_seed"]},
           "versions": {"zombi2": meta["zombi2"], "numpy": meta["numpy"]},
           "zombi2_call": "species.simulate_species_tree + genomes.simulate_genomes_ordered (serial engine"
                          f"{', joint=True' if meta['joint'] else ''}); write(compact=True)",
           "scramble_stem": meta["stem"], "root_node": root_node, "burst_tips": burst_tips,
           "scramble": stats, "taxonomy": {"crown": tx["crown"], "present": tx["present"],
                                           "cut_times": tx["cut_times"]},
           "files": files, "digest": digest}
    if extra:
        doc.update(extra)
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
        with open(bdir / name, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            if tuple(header) != cls._fields:
                raise ValueError(f"{bdir / name}: header {header}, expected {list(cls._fields)}")
            t[name] = [cls(*(c(v) for c, v in zip(conv, row))) for row in reader]
    return Backbone(params=doc["preset"], seed=doc["seed"], root_node=doc["root_node"],
                    burst_tips=doc["burst_tips"], genomes=t["genomes.tsv"], genes=t["genes.tsv"],
                    families=t["families.tsv"], groups=t["groups.tsv"], units=t["truth_units.tsv"],
                    adjacency=t["truth_adjacency.tsv"], scramble=doc["scramble"], digest=doc["digest"],
                    out=out, replicons=t["replicons.tsv"], taxonomy=t["taxonomy.tsv"], taxa=t["taxa.tsv"],
                    events=t["events.tsv"], gene_tree_nodes=t["gene_tree_nodes.tsv"], meta=doc)


__all__ = ["Backbone", "Dec15Error", "SchemaError", "SpeciesBranch", "build", "build_phages", "read", "load_preset",
           "check_dec15", "check_preset", "check_schema", "scramble_stats", "check_scrambled", "cut_taxonomy",
           "name_taxa", "parse_newick", "species_branches"]
