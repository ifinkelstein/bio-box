"""Ring R1: the synthworld backbone's pure parts and preset guards, session 2 (synthworld_spec §2.1 S1-S2,
§2.2, §5.1; card F-02). No ZOMBI2 here: the worlds are tested in tests/r2/test_synthworld_backbone.py.
"""

from __future__ import annotations

import copy

import pytest

from biobox.synthworld import backbone

pytestmark = pytest.mark.covers("F-02")  # ring R1 (by directory)


# --- Newick (the gene trees of §2.2) ---------------------------------------------------------------------

def test_parse_newick_numbers_nodes_in_preorder_with_parents_and_lengths():
    text = "((n3_g1:0.5,n4_g2:0.5)speciation_n2:1.25,e5_g3:0.75)transfer_n1:1;"
    got = backbone.parse_newick(text)
    assert got == [
        (-1, "transfer_n1", 1.0),
        (0, "speciation_n2", 1.25),
        (1, "n3_g1", 0.5),
        (1, "n4_g2", 0.5),
        (0, "e5_g3", 0.75),
    ]


def test_parse_newick_reads_a_single_leaf_tree():
    assert backbone.parse_newick("n3_g8132:0.3828482;") == [(-1, "n3_g8132", 0.3828482)]


# --- taxonomy from depth cuts (§2.1 S1) ------------------------------------------------------------------

def test_taxonomy_takes_the_branch_alive_at_each_cut():
    # stem n0 [0, 1); crown at 1, present at 5. n1 [1, 2) -> tips a, b [2, 5); n2 [1, 4.6) -> c, d [4.6, 5)
    nodes = {"n0": (None, 0.0, 1.0), "n1": ("n0", 1.0, 2.0), "n2": ("n0", 1.0, 4.6),
             "a": ("n1", 2.0, 5.0), "b": ("n1", 2.0, 5.0), "c": ("n2", 4.6, 5.0), "d": ("n2", 4.6, 5.0)}
    # cuts at 0.1, 0.5 and 0.95 of the crown age: times 1.4, 3.0 and 4.8
    got = backbone.cut_taxonomy(nodes, ["a", "b", "c", "d"], crown=1.0, present=5.0, cuts=[0.1, 0.5, 0.95])
    assert got == {"a": ["n1", "a", "a"], "b": ["n1", "b", "b"], "c": ["n2", "n2", "c"],
                   "d": ["n2", "n2", "d"]}


def test_taxon_names_number_taxa_by_their_first_genome():
    per_tip = {"g0": ["n7", "n9"], "g1": ["n5", "n5"], "g2": ["n7", "n8"]}
    names, taxa = backbone.name_taxa(per_tip, domain="Bacteria", ranks=["phylum", "genus"],
                                     cut_times=[1.5, 3.0])
    assert names == {"g0": "d__Bacteria;p__Phylum01;g__Genus01",
                     "g1": "d__Bacteria;p__Phylum02;g__Genus02",
                     "g2": "d__Bacteria;p__Phylum01;g__Genus03"}
    assert [(t.rank, t.taxon, t.node, t.n_genomes, t.parent) for t in taxa] == [
        ("phylum", "p__Phylum01", "n7", 2, "d__Bacteria"),
        ("phylum", "p__Phylum02", "n5", 1, "d__Bacteria"),
        ("genus", "g__Genus01", "n9", 1, "p__Phylum01"),
        ("genus", "g__Genus02", "n5", 1, "p__Phylum02"),
        ("genus", "g__Genus03", "n8", 1, "p__Phylum01"),
    ]


# --- the micro preset's session-2 parameters (§5.1) and their guards ------------------------------------

def test_micro_preset_has_the_planted_backbone_features_of_section_5_1():
    p = backbone.load_preset("micro")
    backbone.check_preset(p)
    sizes = p["plasmids"]["sizes"]
    assert 1 <= len(sizes) <= 2 and all(20 <= s <= 200 for s in sizes)          # 0-2 plasmids, 20-200 genes
    assert len(p["vehicles"]["cargo"]) == 3 and all(3 <= k <= 8 for k in p["vehicles"]["cargo"])
    assert p["hubs"]["count"] == 5 and len(p["decoys"]["partner_modules"]) == 5
    assert p["genomes"]["transfer_to"] == "distance" and p["genomes"]["transfer_decay"] > 0
    assert p["taxonomy"]["ranks"] == ["phylum", "class", "order", "family", "genus", "species"]
    assert any(ax["sigma"] > 0 for ax in p["lineage_rates"]["axes"])


@pytest.mark.parametrize("change, match", [
    (lambda p: p["taxonomy"].update(cuts=[0.5, 0.3, 0.6, 0.7, 0.8, 0.9]), "cuts"),
    (lambda p: p["taxonomy"].update(cuts=[0.1, 0.2]), "cuts"),
    (lambda p: p["plasmids"].update(sizes=[1]), "plasmid"),
    (lambda p: p["vehicles"].update(cargo=[0]), "vehicle"),
    (lambda p: p["decoys"].update(partner_modules=[99]), "decoy"),
    (lambda p: p["genomes"].update(transfer_to="nearby"), "transfer_to"),
    (lambda p: p["lineage_rates"].update(axes=[{"rates": ["loss"], "sigma": 0.2},
                                               {"rates": ["loss"], "sigma": 0.2}]), "lineage"),
    (lambda p: p["lineage_rates"].update(axes=[{"rates": ["inversion"], "sigma": 0.2}]), "lineage"),
])
def test_check_preset_refuses_a_malformed_session_2_parameter(change, match):
    p = copy.deepcopy(backbone.load_preset("micro"))
    change(p)
    with pytest.raises(ValueError, match=match):
        backbone.check_preset(p)


def test_decoys_need_the_modules_they_partner():
    p = copy.deepcopy(backbone.load_preset("micro"))
    p["modules"]["count"] = 0
    with pytest.raises(ValueError, match="decoy"):
        backbone.check_preset(p)
    p["decoys"]["partner_modules"] = []
    backbone.check_preset(p)
