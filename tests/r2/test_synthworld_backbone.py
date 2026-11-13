"""Ring R2: the synthworld backbone on the `micro` preset (synthworld_spec §2.1, §2.2, §5.1, §6.1; card F-02).

Each test builds a micro backbone (ZOMBI2 0.54.0, about 10 s) or reads one built by a module fixture. The
oracles here read ZOMBI2's raw files (`<out>/zombi2/`) with their own parsing, never the backbone's.

DEC-15's three conditions each have a test that fails when the condition is broken:
(a) `fates="group"`: every unit of an adjacency-scored module is intact in every extant genome, while the
    same world with `fates="family"` splits units (as the spike measured), and the truth says so;
(b) replacement 0 in a world with groups: the preset guard (tests/r1) and ZOMBI2 itself both refuse it;
(c) the root order and strands are scrambled, read from ZOMBI2's own initial and root layouts.

Session 2 adds a test per planted feature of the rest of the backbone (plasmids, vehicles, hubs, decoys,
distance decay, lineage-varying rates, the synthetic taxonomy) and for the typed event and gene-tree
tables; each recovers its property from the raw files. The determinism tests above cover the new files,
since they compare every file of the world.
"""

from __future__ import annotations

import collections
import csv
import hashlib
import importlib.util
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import pytest

from biobox.synthworld import backbone

pytestmark = [pytest.mark.r2, pytest.mark.covers("F-02")]
needs_zombi2 = pytest.mark.skipif(
    importlib.util.find_spec("zombi2") is None,
    reason="zombi2 0.54.0 is not importable here: CI installs it by `uv sync`; on the pod add its own "
           "directory with BIOBOX_POD_PYDEPS (card F-02)")

REPO = Path(__file__).resolve().parents[2]


def _tree_sha(root: Path) -> dict[str, str]:
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file()}


def _tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def micro(tmp_path_factory):
    if importlib.util.find_spec("zombi2") is None:
        pytest.skip("zombi2 0.54.0 is not importable here")
    out = tmp_path_factory.mktemp("micro") / "world"
    params = backbone.load_preset("micro")
    t0 = time.perf_counter()
    bb = backbone.build(params, params["default_seed"], out)
    return bb, out, time.perf_counter() - t0


@pytest.fixture(scope="module")
def micro_family_fates(tmp_path_factory):
    if importlib.util.find_spec("zombi2") is None:
        pytest.skip("zombi2 0.54.0 is not importable here")
    out = tmp_path_factory.mktemp("micro_family") / "world"
    params = backbone.load_preset("micro", overrides={"modules": {"fates": "family", "scored": False}})
    return backbone.build(params, params["default_seed"], out), out


# --- determinism (§6.1, §6.2; MR-16) -------------------------------------------------------------------

@needs_zombi2
def test_same_seed_gives_byte_identical_backbone_twice(micro, tmp_path):
    bb, out, _ = micro
    again = backbone.build(backbone.load_preset("micro"), bb.seed, tmp_path / "again")
    assert _tree_sha(tmp_path / "again") == _tree_sha(out)
    assert again.digest == bb.digest


@needs_zombi2
@pytest.mark.covers("F-02", "MR-16")
def test_same_seed_in_a_subprocess_with_other_hash_seed_and_threads(micro, tmp_path):
    bb, out, _ = micro
    dst = tmp_path / "sub"
    env = dict(os.environ, PYTHONHASHSEED="4242", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    code = ("import sys; from biobox.synthworld import backbone; "
            "backbone.build(backbone.load_preset('micro'), int(sys.argv[1]), sys.argv[2])")
    subprocess.run([sys.executable, "-c", code, str(bb.seed), str(dst)], check=True, env=env, cwd=tmp_path)
    assert _tree_sha(dst) == _tree_sha(out)


@needs_zombi2
def test_a_different_seed_gives_a_different_backbone(micro, tmp_path):
    bb, out, _ = micro
    other = backbone.build(backbone.load_preset("micro"), bb.seed + 1, tmp_path / "other")
    assert other.digest != bb.digest
    assert (tmp_path / "other/zombi2/species_complete.nwk").read_bytes() != \
        (out / "zombi2/species_complete.nwk").read_bytes()
    assert [g[:6] for g in other.genes[:200]] != [g[:6] for g in bb.genes[:200]]


# --- DEC-15 (a): fates="group" keeps planted units whole; fates="family" does not ------------------------

def _raw_units(out: Path, genomes) -> dict[tuple[str, int], list[tuple[int, int, int, int]]]:
    """(species node, unit) -> [(chromosome, position, strand, family)], read from the raw files."""
    groups = _tsv(out / "zombi2/groups.tsv")
    unit_of = {int(r["gene"]): int(r["unit"]) for r in groups}
    extant = {g.species_node for g in genomes}
    units = collections.defaultdict(list)
    for r in _tsv(out / "zombi2/gene_order.tsv"):
        if r["lineage"] not in extant:
            continue
        gid = int(r["copy"][1:])
        if gid in unit_of:
            units[(r["lineage"], unit_of[gid])].append(
                (int(r["chromosome"]), int(r["position"]), int(r["strand"]), int(r["family"])))
    return units


def _raw_intact(members, n_genes_by_chrom, order) -> bool:
    """Spike §1 rule, written separately: one chromosome, contiguous (ring wrap allowed), one strand,
    declared order (reversed on the minus strand)."""
    chroms = {m[0] for m in members}
    strands = {m[2] for m in members}
    if len(chroms) != 1 or len(strands) != 1:
        return False
    n = n_genes_by_chrom[next(iter(chroms))]
    pos = sorted(m[1] for m in members)
    k = len(pos)
    contiguous = any(all((pos[(s + i) % k] - pos[s]) % n == i for i in range(k)) for s in range(k))
    if not contiguous:
        return False
    by_pos = {m[1]: m[3] for m in members}
    start = next(s for s in range(k) if all((pos[(s + i) % k] - pos[s]) % n == i for i in range(k)))
    fams = [by_pos[pos[(start + i) % k]] for i in range(k)]
    if strands == {-1}:
        fams.reverse()
    idx = [order[f] for f in fams]
    return idx == sorted(idx)


@needs_zombi2
@pytest.mark.covers("F-02", "DEC-15")
def test_dec15_a_every_scored_unit_is_intact_in_every_extant_genome(micro):
    bb, out, _ = micro
    scored = [u for u in bb.units if u.n_members >= 2 and bb.group(u.group).scored]
    assert len(scored) >= 100, f"too few scored units to say anything: {len(scored)}"
    assert all(u.intact for u in scored), [u for u in scored if not u.intact][:5]


@needs_zombi2
@pytest.mark.covers("F-02", "DEC-15")
def test_dec15_a_fates_family_splits_a_unit_and_the_truth_says_so(micro_family_fates):
    bb, out = micro_family_fates
    multi = [u for u in bb.units if u.n_members >= 2]
    broken = [u for u in multi if not u.intact]
    assert broken, "fates='family' left every unit intact; the spike measured 10-14% broken"
    # the realized-order truth agrees with an independent reading of ZOMBI2's raw files
    raw = _raw_units(out, bb.genomes)
    n_by = collections.Counter((r["lineage"], int(r["chromosome"]))
                               for r in _tsv(out / "zombi2/gene_order.tsv"))
    order = {f.family: f.order_index for f in bb.families if f.group}
    node_of = {g.genome_id: g.species_node for g in bb.genomes}
    for u in multi:
        members = raw[(node_of[u.genome_id], u.unit)]
        n_genes = {c: n_by[(node_of[u.genome_id], c)] for c in {m[0] for m in members}}
        assert len(members) == u.n_members
        assert _raw_intact(members, n_genes, order) == u.intact, u
    # and an adjacent pair split by an intruder is not reported as a within-unit adjacency
    same_unit = {(a.genome_id, a.replicon, a.left_position) for a in bb.adjacency if a.same_unit}
    assert sum(u.n_members - 1 for u in multi if u.intact) <= len(same_unit) <= \
        sum(u.n_members - 1 for u in multi)


# --- DEC-15 (b): replacement is refused with groups --------------------------------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "DEC-15")
def test_dec15_b_replacement_with_groups_never_reaches_a_world(micro, tmp_path):
    bb, _, _ = micro
    assert bb.params["genomes"]["replacement"] == 0.0
    p = backbone.load_preset("micro", overrides={"genomes": {"replacement": 0.5}})
    with pytest.raises(backbone.Dec15Error, match=r"\(b\)"):
        backbone.build(p, bb.seed, tmp_path / "replacing")
    assert not (tmp_path / "replacing").exists()


# --- DEC-15 (c): root order and strands scrambled --------------------------------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "DEC-15")
def test_dec15_c_root_order_and_strands_are_scrambled(micro):
    bb, out, _ = micro
    init = _tsv(out / "zombi2/initial_genome.tsv")
    root_label = bb.root_node
    root = [r for r in _tsv(out / "zombi2/gene_order.tsv") if r["lineage"] == root_label]
    assert len(root) == len(init) and sorted(r["family"] for r in root) == sorted(r["family"] for r in init)
    # ZOMBI2's default: every gene on +, declared families in one block at the end of the chromosome
    assert {r["strand"] for r in init} == {"1"}
    # every declared family: group members, and since session 2 the hubs and decoys
    declared = {f.family for f in bb.families if f.kind not in ("root", "originated")}
    tail = [int(r["family"]) for r in init[-len(declared):]]
    assert set(tail) == declared
    # after the synthworld scramble stem
    plus = sum(r["strand"] == "1" for r in root) / len(root)
    assert 0.4 <= plus <= 0.6, plus

    def pairs(rows):
        fam = [int(r["family"]) for r in sorted(rows, key=lambda r: int(r["position"]))]
        return {frozenset((fam[i], fam[(i + 1) % len(fam)])) for i in range(len(fam))}

    unit_pairs = set()
    for g in bb.groups:
        if g.adjacent:
            unit_pairs |= {frozenset(p) for p in zip(g.members, g.members[1:])}
    before = pairs(init) - unit_pairs
    kept = len(before & pairs(root)) / len(before)
    assert kept <= 0.05, kept
    # declared genes no longer sit next to each other more than chance would put them
    fam_root = [int(r["family"]) for r in sorted(root, key=lambda r: int(r["position"]))]
    n = len(fam_root)
    touching = [frozenset((fam_root[i], fam_root[(i + 1) % n])) for i in range(n)]
    with_declared = [p for p in touching if p & declared and p not in unit_pairs]
    clustered = sum(p <= declared for p in with_declared) / len(with_declared)
    assert clustered <= 0.5, clustered


# --- truth from the realized order -----------------------------------------------------------------------

def _plasmid_sizes(params) -> dict[str, int]:
    """Plasmid group name -> declared size, from the preset (naming rule of backbone.py's docstring)."""
    return {f"pls{i}": int(s) for i, s in enumerate(params["plasmids"]["sizes"])}


def _raw_replicons(out: Path, extant: set[str], plasmids: set[str]):
    """lineage -> [(kind, name, topology, [(family, strand, copy label), ...])], read from ZOMBI2's raw
    files and excised here independently: each unit of a plasmid group leaves its chromosome and
    becomes a replicon of its own, in ZOMBI2's order; the chromosome keeps the other genes in order."""
    unit_of = {int(r["gene"]): (r["group"], int(r["unit"])) for r in _tsv(out / "zombi2/groups.tsv")
               if r["group"] in plasmids}
    rows = collections.defaultdict(list)
    for r in _tsv(out / "zombi2/gene_order.tsv"):
        if r["lineage"] in extant:
            rows[(r["lineage"], int(r["chromosome"]))].append(r)
    reps = collections.defaultdict(list)
    for (lineage, _chrom), rs in sorted(rows.items()):
        rs.sort(key=lambda r: int(r["position"]))
        chromosome, units = [], {}
        for r in rs:
            gene = (int(r["family"]), int(r["strand"]), f"{lineage}_{r['copy']}")
            u = unit_of.get(int(r["copy"][1:]))
            (units.setdefault(u, []) if u else chromosome).append(gene)
        reps[lineage].append(("chromosome", "chromosome", rs[0]["topology"], chromosome))
        reps[lineage] += [("plasmid", u[0], "circular", gs) for u, gs in units.items()]
    return reps


@needs_zombi2
def test_truth_adjacency_equals_a_bruteforce_recount(micro):
    bb, out, _ = micro
    extant = {g.species_node: g.genome_id for g in bb.genomes}
    reps = _raw_replicons(out, set(extant), set(_plasmid_sizes(bb.params)))
    want = collections.Counter()
    for lineage, rs in reps.items():
        for _kind, _name, topo, genes in rs:
            n = len(genes)
            last = n if topo == "circular" and n >= 3 else n - 1
            for p in range(last):
                a, b = genes[p][:2], genes[(p + 1) % n][:2]
                want[(extant[lineage], a, b)] += 1
    got = collections.Counter((a.genome_id, (a.left_family, a.left_strand), (a.right_family, a.right_strand))
                              for a in bb.adjacency)
    assert got == want
    assert sum(want.values()) > 40 * 1000


@needs_zombi2
def test_copies_map_one_to_one_to_extant_gene_tree_leaves(micro):
    bb, out, _ = micro
    leaves = set()
    fam_of_leaf = {}
    for r in _tsv(out / "zombi2/gene_trees.tsv"):
        for leaf in re.findall(r"(n\d+_g\d+)", r["extant"]):
            leaves.add(leaf)
            fam_of_leaf[leaf] = int(r["family"])
    labels = [g[6] for g in bb.genes]
    assert len(labels) == len(set(labels)) == len(leaves)
    assert set(labels) == leaves
    assert all(fam_of_leaf[g[6]] == g[5] for g in bb.genes)
    counts = collections.Counter(g[5] for g in bb.genes)
    assert all(f.n_extant_copies == counts.get(f.family, 0) for f in bb.families)


@needs_zombi2
def test_group_membership_matches_raw_groups_file(micro):
    bb, out, _ = micro
    unit_of = {int(r["gene"]): (r["group"], int(r["unit"])) for r in _tsv(out / "zombi2/groups.tsv")}
    for g in bb.genes:
        gid = int(g[6].split("_g")[1])
        assert (g[7], g[8]) == unit_of.get(gid, ("", -1))


@needs_zombi2
def test_schema_pins_refuse_a_changed_zombi2_header(micro, tmp_path):
    _, out, _ = micro
    backbone.check_schema(out / "zombi2")
    bad = tmp_path / "z"
    bad.mkdir()
    for p in (out / "zombi2").iterdir():
        text = p.read_text(encoding="utf-8")
        if p.name == "genome_events.tsv":
            text = text.replace("\tevent\n", "\n", 1)
        (bad / p.name).write_text(text, encoding="utf-8")
    with pytest.raises(backbone.SchemaError, match="genome_events.tsv"):
        backbone.check_schema(bad)


# --- the micro preset's backbone parameters (§5.1) and the runtime budget (§2.6) --------------------------

@needs_zombi2
def test_micro_backbone_parameters(micro):
    bb, _, _ = micro
    assert len(bb.genomes) == 46
    assert len(bb.burst_tips) == 4 and set(bb.burst_tips) <= {g.species_node for g in bb.genomes}
    planted = [g for g in bb.groups if g.kind == "planted"]
    assert len(planted) == 25 and {len(g.members) for g in planted} == {2, 3, 4, 5, 6}
    assert len([g for g in bb.groups if g.kind == "profile_only"]) == 3
    # §5.1 "chromosome mean ≈ 1,500 genes": plasmids are replicons of their own since session 2
    sizes = [r.n_genes for r in bb.replicons if r.kind == "chromosome"]
    assert len(sizes) == len(bb.genomes)
    assert all(g.n_genes == sum(r.n_genes for r in bb.replicons if r.genome_id == g.genome_id)
               for g in bb.genomes)
    mean = sum(sizes) / len(sizes)
    assert 0.7 * 1500 <= mean <= 1.3 * 1500, (mean, min(sizes), max(sizes))


@needs_zombi2
def test_micro_backbone_runtime_within_budget(micro, record_property):
    bb, _, seconds = micro
    record_property("backbone_seconds", round(seconds, 2))
    sizes = [r.n_genes for r in bb.replicons if r.kind == "chromosome"]
    n_pls = collections.Counter(r.genome_id for r in bb.replicons if r.kind == "plasmid")
    print(f"\nF-02 micro backbone: {seconds:.2f} s wall; chromosome genes mean {sum(sizes) / len(sizes):.0f} "
          f"min {min(sizes)} max {max(sizes)}; plasmids per genome "
          f"{sorted(collections.Counter(n_pls[g.genome_id] for g in bb.genomes).items())}; "
          f"seed {bb.seed}; digest {bb.digest}")
    assert seconds < 60.0  # synthworld_spec §1.4 item 1 (pass under 60 s); §2.6 budget: seconds


# --- session 2: the rest of the backbone (§2.1 S1, S2) ----------------------------------------------------
# Each planted feature is recovered from ZOMBI2's raw files by a brute-force reading of its own.

def _names(out: Path) -> dict[str, int]:
    import json
    return {k: int(v) for k, v in
            json.loads((out / "zombi2/genome_summary.json").read_text(encoding="utf-8"))["families"]["names"].items()}


def _newick_times(path: Path, t0: float):
    """Species tree from ZOMBI2's Newick, parsed here: label -> (parent label, start time, end time)."""
    text = path.read_text(encoding="utf-8").strip()
    nodes: dict[str, list] = {}       # label -> [parent, branch length, children]
    opened: list[list[str]] = []      # the children read so far of each open "("
    closed = False                    # the next label names the node a ")" just closed
    for tok in re.findall(r"\(|\)|,|;|[^(),;]+", text):
        if tok == "(":
            opened.append([])
            closed = False
        elif tok == ")":
            closed = True
        elif tok in (",", ";"):
            closed = False
        else:
            label, _, length = tok.partition(":")
            kids = opened.pop() if closed else []
            nodes[label] = [None, float(length or 0.0), kids]
            for k in kids:
                nodes[k][0] = label
            if opened:
                opened[-1].append(label)
            closed = False
    root = next(k for k, v in nodes.items() if v[0] is None)
    times: dict[str, tuple] = {}
    todo = [(root, t0)]
    while todo:
        lab, start = todo.pop()
        parent, length, kids = nodes[lab]
        times[lab] = (parent, start, start + length)
        todo.extend((k, start + length) for k in kids)
    return times


def _t0(out: Path) -> float:
    with open(out / "zombi2/genome_events.tsv", encoding="utf-8") as fh:
        fh.readline()
        return float(fh.readline().split("\t")[0])


@pytest.fixture(scope="module")
def micro_uniform(tmp_path_factory):
    """The micro world with uniform transfer recipients instead of distance decay."""
    if importlib.util.find_spec("zombi2") is None:
        pytest.skip("zombi2 0.54.0 is not importable here")
    out = tmp_path_factory.mktemp("micro_uniform") / "world"
    params = backbone.load_preset("micro", overrides={"genomes": {"transfer_to": "uniform"}})
    return backbone.build(params, params["default_seed"], out), out


@needs_zombi2
@pytest.mark.covers("F-02", "B-SIM")
def test_plasmids_are_whole_replicons_of_their_preset_size(micro):
    """§2.1 S2 "chromosome plus 0-2 plasmid replicons"; §5.1 plasmids of 20-200 genes. Tolerance: none,
    every extant plasmid has exactly its declared genes."""
    bb, out, _ = micro
    sizes = _plasmid_sizes(bb.params)
    node_of = {g.genome_id: g.species_node for g in bb.genomes}
    reps = _raw_replicons(out, set(node_of.values()), set(sizes))
    raw_plasmids = {lin: [(name, {c for _, _, c in gs}) for kind, name, _, gs in rs if kind == "plasmid"]
                    for lin, rs in reps.items()}
    # every extant plasmid unit is whole: its declared number of genes, contiguous in ZOMBI2's order
    for lin, ps in raw_plasmids.items():
        for name, copies in ps:
            assert len(copies) == sizes[name], (lin, name, len(copies))
    # and is excised as a replicon of its own holding exactly those genes
    got = collections.defaultdict(set)
    for g in bb.genes:
        got[(g.genome_id, g.replicon)].add(g.copy_label)
    plasmid_reps = [r for r in bb.replicons if r.kind == "plasmid"]
    for r in plasmid_reps:
        assert r.topology == "circular" and r.n_genes == sizes[r.name] == len(got[(r.genome_id, r.replicon)])
    want = sorted((gid, name, tuple(sorted(c))) for gid, node in node_of.items()
                  for name, c in raw_plasmids[node])
    assert sorted((r.genome_id, r.name, tuple(sorted(got[(r.genome_id, r.replicon)])))
                  for r in plasmid_reps) == want
    # the chromosome keeps the rest, in ZOMBI2's order
    for gid, node in node_of.items():
        chrom = [r for r in bb.replicons if r.genome_id == gid and r.kind == "chromosome"]
        assert len(chrom) == 1 and chrom[0].replicon == 0
        raw_chrom = [c for kind, _, _, gs in reps[node] if kind == "chromosome" for _, _, c in gs]
        ours = [g.copy_label for g in bb.genes if g.genome_id == gid and g.replicon == 0]
        assert ours == raw_chrom
    # both plasmid types are realized, some genomes carry none, none carries more than two (§5.1)
    per = collections.Counter(r.genome_id for r in plasmid_reps)
    assert {r.name for r in plasmid_reps} == set(sizes)
    assert any(per[gid] == 0 for gid in node_of) and max(per.values()) <= 2, sorted(per.values())


@needs_zombi2
@pytest.mark.covers("F-02", "B-SIM")
def test_a_vehicles_cargo_moves_with_its_chassis(micro):
    """§2.1 S2 vehicles (chassis + cargo): every transfer event that moves a vehicle gene moves every gene
    of that vehicle (read off `event` ids), and every extant cargo copy sits between its two chassis."""
    bb, out, _ = micro
    names = _names(out)
    vehicles = {f"veh{i}": [names[f"veh{i}.{j}"] for j in range(int(k) + 2)]
                for i, k in enumerate(bb.params["vehicles"]["cargo"])}
    by_event = collections.defaultdict(set)
    for r in _tsv(out / "zombi2/genome_events.tsv"):
        if r["kind"].startswith("transfer"):
            by_event[r["event"]].add(int(r["family"]))
    moved = collections.Counter()
    for fams in by_event.values():
        for v, members in vehicles.items():
            if fams & set(members):
                assert set(members) <= fams, (v, sorted(fams))
                moved[v] += 1
    assert all(moved[v] >= 3 for v in vehicles), moved
    # in every extant genome, cargo only inside a whole unit with a chassis gene at each end
    unit_of = {int(r["gene"]): (r["group"], int(r["unit"])) for r in _tsv(out / "zombi2/groups.tsv")}
    extant = {g.species_node for g in bb.genomes}
    cargo = {f: v for v, ms in vehicles.items() for f in ms[1:-1]}
    units = collections.defaultdict(list)
    n_cargo = 0
    for r in _tsv(out / "zombi2/gene_order.tsv"):
        if r["lineage"] in extant and int(r["family"]) in cargo:
            n_cargo += 1
        u = unit_of.get(int(r["copy"][1:]))
        if r["lineage"] in extant and u and u[0] in vehicles:
            units[(r["lineage"], u[1])].append((int(r["position"]), int(r["family"]), int(r["strand"])))
    assert n_cargo == sum(len(ms) - 2 for ms in units.values()) > 0
    for (lin, _u), ms in units.items():
        ms.sort()
        fams = [f for _, f, _ in ms]
        v = cargo.get(fams[1]) or cargo.get(fams[-2])
        want = vehicles[v] if ms[0][2] == 1 else vehicles[v][::-1]
        assert fams == want and ms[-1][0] - ms[0][0] == len(ms) - 1, (lin, fams)


@needs_zombi2
@pytest.mark.covers("F-02", "B-SIM")
def test_hub_families_have_more_diverse_neighbours(micro):
    """§2.1 S2 hub families; §4.2 item 4. Neighbour diversity of a family = the number of distinct families
    next to any of its extant copies, read from the raw gene order. Margin (asserted): every hub at least
    3x the median root family present in at least 40 of the 46 genomes. (Per copy the gap is small: a
    transfer of extent 3 that starts on a hub carries its right-hand neighbours with it.)"""
    bb, out, _ = micro
    names = _names(out)
    hubs = [names[f"hub{i}"] for i in range(int(bb.params["hubs"]["count"]))]
    extant = {g.species_node for g in bb.genomes}
    chroms = collections.defaultdict(list)
    for r in _tsv(out / "zombi2/gene_order.tsv"):
        if r["lineage"] in extant:
            chroms[(r["lineage"], r["chromosome"])].append((int(r["position"]), int(r["family"])))
    neighbours = collections.defaultdict(set)
    copies = collections.Counter()
    genomes = collections.defaultdict(set)
    for (lin, _c), genes in chroms.items():
        fam = [f for _, f in sorted(genes)]
        n = len(fam)
        for i, f in enumerate(fam):
            copies[f] += 1
            genomes[f].add(lin)
            neighbours[f] |= {fam[i - 1], fam[(i + 1) % n]}
    diversity = {f: len(neighbours[f]) for f in copies}
    n_root = int(bb.params["genomes"]["initial_families"])
    roots = sorted(diversity[f] for f in copies if f < n_root and len(genomes[f]) >= 40)
    assert len(roots) >= 100
    median = roots[len(roots) // 2]
    got = {h: diversity.get(h, 0) for h in hubs}
    print(f"\nhub neighbour families {got} over copies {[copies[h] for h in hubs]}; root median {median} "
          f"over {len(roots)} families")
    assert all(copies[h] >= 46 for h in hubs), {h: copies[h] for h in hubs}
    assert min(got.values()) >= 3.0 * median, (got, median)


@needs_zombi2
@pytest.mark.covers("F-02", "B-SIM")
def test_decoy_transfers_land_only_where_the_partner_is(micro):
    """§2.1 S2 `Recipients` decoys: a decoy's transfers are weighted by the presence of a planted module's
    family in the recipient, at the micro weight 0 for an absent partner. Replayed from the raw events:
    every decoy transfer lands in a lineage that held a live copy of its partner just before."""
    bb, out, _ = micro
    names = _names(out)
    d = bb.params["decoys"]
    assert float(d["absent_weight"]) == 0.0
    partner_of = {names[f"dec{i}"]: names[f"mod{m:02d}.0"] for i, m in enumerate(d["partner_modules"])}
    times = _newick_times(out / "zombi2/species_complete.nwk", _t0(out))
    live = {p: set() for p in partner_of.values()}
    landed = collections.Counter()
    shares = []
    for r in _tsv(out / "zombi2/genome_events.tsv"):
        fam = int(r["family"])
        kids = [c for c in r["children"].split(";") if c]
        if fam in partner_of and r["kind"].startswith("transfer"):
            recipient = kids[1].split("_")[0]
            carriers = {c.split("_")[0] for c in live[partner_of[fam]]}
            landed[(fam, recipient in carriers)] += 1
            t, donor = float(r["time"]), r["parents"].split("_")[0]
            others = {lab for lab, (_, s, e) in times.items() if s <= t < e and lab != donor}
            shares.append(len(others & carriers) / len(others))
        if fam in live:
            if r["parents"]:
                live[fam].discard(r["parents"])
            live[fam].update(kids)
    mean_share = sum(shares) / len(shares)
    print(f"\ndecoy transfers (family, partner present): {sorted(landed.items())}; partner carriers among "
          f"the other contemporaneous lineages at those moments: mean {mean_share:.3f}")
    assert not any(not ok for _, ok in landed), landed
    assert all(landed[(f, True)] >= 3 for f in partner_of), landed
    # the rule matters: uniform recipients would have missed the carriers about this often
    assert mean_share <= 0.9, mean_share


def _transfer_ranks(out: Path) -> list[tuple[float, float]]:
    """For each transfer event (one per `event` id), among every other lineage alive at that time but the
    donor: the fraction strictly closer to the donor (patristic distance) than the recipient, and the
    fraction strictly farther. With uniform recipients the two have the same expectation (each pair of
    candidates at different distances adds one to each side), ties notwithstanding."""
    times = _newick_times(out / "zombi2/species_complete.nwk", _t0(out))
    seen, ranks = set(), []
    for r in _tsv(out / "zombi2/genome_events.tsv"):
        if not r["kind"].startswith("transfer") or r["event"] in seen:
            continue
        seen.add(r["event"])
        t = float(r["time"])
        donor = r["parents"].split("_")[0]
        recipient = r["children"].split(";")[1].split("_")[0]
        anc = {}
        x = donor
        while x is not None:
            anc[x] = times[x][2]
            x = times[x][0]
        alive = [lab for lab, (_, s, e) in times.items() if s <= t < e and lab != donor]

        def dist(lab):
            y = lab
            while y not in anc:
                y = times[y][0]
            return 2.0 * (t - anc[y])

        d = {lab: dist(lab) for lab in alive}
        if len(alive) < 2:
            continue
        dr = d[recipient]
        closer = sum(v < dr - 1e-12 for v in d.values())
        farther = sum(v > dr + 1e-12 for v in d.values())
        ranks.append((closer / (len(alive) - 1), farther / (len(alive) - 1)))
    return ranks


@needs_zombi2
@pytest.mark.covers("F-02", "B-SIM")
def test_transfers_decay_with_distance_and_not_without_it(micro, micro_uniform):
    """§2.1 S2 distance-decay transfer; §4.2 item 2. Summed over transfer events, the recipient is closer
    to the donor than other contemporaneous lineages at least 2x as often as it is farther under the micro
    decay; under uniform recipients the two sums are within 10% of each other."""
    _, out, _ = micro
    _, out_u = micro_uniform
    decay = _transfer_ranks(out)
    uniform = _transfer_ranks(out_u)

    def close_far(rs):
        return sum(f for _, f in rs), sum(c for c, _ in rs)

    c, f = close_far(decay)
    cu, fu = close_far(uniform)
    print(f"\ntransfer recipients closer/farther than other lineages: decay {c:.0f}/{f:.0f} over "
          f"{len(decay)} events; uniform {cu:.0f}/{fu:.0f} over {len(uniform)}")
    assert len(decay) >= 1000 and len(uniform) >= 1000
    assert c >= 2.0 * f, (c, f)
    assert 0.9 <= cu / fu <= 1.1, (cu, fu)


@needs_zombi2
@pytest.mark.covers("F-02", "B-SIM")
def test_lineage_varying_rates_reach_the_genome_size_range(micro):
    """§5.1 genome size: chromosomes of about 1,000 to 4,000 genes (uniform rates gave 1,393-1,708,
    session 1). Reached when the smallest chromosome is at most 1,100 and the largest at least 3,600
    (within 10% of the range's ends), read from the raw gene order with plasmids excised."""
    bb, out, _ = micro
    extant = {g.species_node for g in bb.genomes}
    reps = _raw_replicons(out, extant, set(_plasmid_sizes(bb.params)))
    sizes = sorted(len(gs) for rs in reps.values() for kind, _, _, gs in rs if kind == "chromosome")
    print(f"\nchromosome genes: min {sizes[0]} median {sizes[len(sizes) // 2]} max {sizes[-1]}")
    assert sizes == sorted(r.n_genes for r in bb.replicons if r.kind == "chromosome")
    assert sizes[0] <= 1100 and sizes[-1] >= 3600, sizes
    # ZOMBI2 recorded a multiplier per species branch for each varying rate, and they differ
    rows = _tsv(out / "zombi2/lineage_multipliers.tsv")
    for axis in bb.params["lineage_rates"]["axes"]:
        for rate in axis["rates"]:
            vals = {float(r[rate]) for r in rows if r[rate]}
            assert len(vals) > 10, (rate, len(vals))


@needs_zombi2
@pytest.mark.covers("F-02", "B-SIM")
def test_taxonomy_is_cut_from_the_species_tree_at_the_preset_depths(micro):
    """§2.1 S1 synthetic taxonomy from depth cuts. Recomputed from the raw Newick: a genome's taxon at a
    rank is the species-tree branch alive at that rank's cut, a fraction of the crown age."""
    bb, out, _ = micro
    tx = bb.params["taxonomy"]
    times = _newick_times(out / "zombi2/species_complete.nwk", _t0(out))
    root = next(lab for lab, v in times.items() if v[0] is None)
    crown = times[root][2]
    present = max(e for _, _, e in times.values())
    taxa = {(t.rank, t.taxon): t for t in bb.taxa}
    lineage_of = {t.genome_id: t.tax_lineage.split(";") for t in bb.taxonomy}
    assert set(lineage_of) == {g.genome_id for g in bb.genomes}
    for g in bb.genomes:
        parts = lineage_of[g.genome_id]
        assert parts[0] == f"d__{tx['domain']}" and len(parts) == 1 + len(tx["ranks"])
        for rank, frac, name in zip(tx["ranks"], tx["cuts"], parts[1:]):
            cut = crown + frac * (present - crown)
            x = g.species_node
            while not (times[x][1] < cut <= times[x][2]):
                x = times[x][0]
            assert taxa[(rank, name)].node == x, (g.genome_id, rank)
    # ranks nest: a taxon has one parent taxon; counts add up; the clonal burst is one species
    for (rank, name), t in taxa.items():
        members = {gid for gid, parts in lineage_of.items() if name in parts}
        assert t.n_genomes == len(members)
        if rank != tx["ranks"][0]:
            assert len({lineage_of[gid][tx["ranks"].index(rank)] for gid in members}) == 1
    burst = {g.genome_id for g in bb.genomes if g.species_node in bb.burst_tips}
    assert len({lineage_of[gid][-1] for gid in burst}) == 1
    n_per_rank = [len({p[i + 1] for p in lineage_of.values()}) for i in range(len(tx["ranks"]))]
    assert n_per_rank == sorted(n_per_rank) and n_per_rank[0] >= 2 and n_per_rank[-1] < len(bb.genomes)


@needs_zombi2
def test_events_table_is_the_raw_genome_events_typed(micro):
    """§2.2: genome_events.tsv as a typed table (speciation rows are left out: the gene trees hold them)."""
    bb, out, _ = micro
    raw = [r for r in _tsv(out / "zombi2/genome_events.tsv") if r["kind"] != "speciation"]
    assert len(bb.events) == len(raw) > 1000
    for e, r in zip(bb.events, raw):
        kids = tuple(c for c in r["children"].split(";") if c)
        assert (e.event, e.kind, e.family, e.parent, e.children) == \
            (int(r["event"]), r["kind"], int(r["family"]), r["parents"], kids)
        assert e.time == float(r["time"])
        acting = (r["parents"] or kids[0]).split("_")[0]
        assert e.lineage == acting
        assert e.recipient == (kids[1].split("_")[0] if r["kind"].startswith("transfer") else "")
    kinds = collections.Counter(e.kind for e in bb.events)
    assert kinds["transfer_additive"] > 0 and kinds["loss"] > 0 and "transfer_replacing" not in kinds


@needs_zombi2
def test_gene_tree_nodes_parse_the_complete_trees(micro):
    """§2.2: gene_trees.tsv as typed nodes. Checks that need a correct parse: every extant copy is one
    extant leaf of its family's tree, and every extant leaf of a tree is the same time from its root
    (all end at the present), within the 7 significant digits ZOMBI2 writes."""
    bb, out, _ = micro
    nodes = collections.defaultdict(list)
    for n in bb.gene_tree_nodes:
        nodes[n.family].append(n)
    raw = {int(r["family"]): r["complete"] for r in _tsv(out / "zombi2/gene_trees.tsv")}
    assert set(nodes) == set(raw)
    extant_leaves = {n.copy_label: n.family for n in bb.gene_tree_nodes if n.extant}
    assert extant_leaves == {g.copy_label: g.family for g in bb.genes}
    for fam, ns in nodes.items():
        labels = re.findall(r"([a-z]+_[ne]\d+|[ne]\d+_g\d+):", raw[fam])
        assert len(ns) == len(labels), fam
        assert [n.node for n in ns] == list(range(len(ns))) and ns[0].parent == -1
        depth = {}
        for n in ns:
            assert n.parent < n.node
            depth[n.node] = (depth[n.parent] if n.parent >= 0 else 0.0) + n.branch_length
        ends = [depth[n.node] for n in ns if n.extant]
        if ends:
            assert max(ends) - min(ends) <= 1e-4 * max(ends), fam
        parents = {n.parent for n in ns}
        for n in ns:
            assert bool(n.copy_label) == (n.node not in parents)
            assert n.kind in (("extant", "lost") if n.copy_label else ("speciation", "duplication", "transfer"))
