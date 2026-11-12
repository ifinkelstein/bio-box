"""Ring R2: the synthworld backbone on the `micro` preset (synthworld_spec §2.1, §2.2, §5.1, §6.1; card F-02).

Each test builds a micro backbone (ZOMBI2 0.54.0, about 10 s) or reads one built by a module fixture. The
oracles here read ZOMBI2's raw files (`<out>/zombi2/`) with their own parsing, never the backbone's.

DEC-15's three conditions each have a test that fails when the condition is broken:
(a) `fates="group"`: every unit of an adjacency-scored module is intact in every extant genome, while the
    same world with `fates="family"` splits units (as the spike measured), and the truth says so;
(b) replacement 0 in a world with groups: the preset guard (tests/r1) and ZOMBI2 itself both refuse it;
(c) the root order and strands are scrambled, read from ZOMBI2's own initial and root layouts.
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
    declared = {f.family for f in bb.families if f.group}
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

@needs_zombi2
def test_truth_adjacency_equals_a_bruteforce_recount(micro):
    bb, out, _ = micro
    extant = {g.species_node: g.genome_id for g in bb.genomes}
    by_chrom = collections.defaultdict(dict)
    for r in _tsv(out / "zombi2/gene_order.tsv"):
        if r["lineage"] in extant:
            by_chrom[(r["lineage"], r["chromosome"], r["topology"])][int(r["position"])] = \
                (int(r["family"]), int(r["strand"]))
    want = collections.Counter()
    for (lineage, _chrom, topo), genes in by_chrom.items():
        n = len(genes)
        last = n if topo == "circular" and n >= 3 else n - 1
        for p in range(last):
            a, b = genes[p], genes[(p + 1) % n]
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
    sizes = [g.n_genes for g in bb.genomes]
    mean = sum(sizes) / len(sizes)
    assert 0.7 * 1500 <= mean <= 1.3 * 1500, (mean, min(sizes), max(sizes))


@needs_zombi2
def test_micro_backbone_runtime_within_budget(micro, record_property):
    bb, _, seconds = micro
    record_property("backbone_seconds", round(seconds, 2))
    sizes = [g.n_genes for g in bb.genomes]
    print(f"\nF-02 micro backbone: {seconds:.2f} s wall; genes/genome mean {sum(sizes) / len(sizes):.0f} "
          f"min {min(sizes)} max {max(sizes)}; seed {bb.seed}; digest {bb.digest}")
    assert seconds < 60.0  # synthworld_spec §1.4 item 1 (pass under 60 s); §2.6 budget: seconds
