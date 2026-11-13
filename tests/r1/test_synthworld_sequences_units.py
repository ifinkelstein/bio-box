"""Ring R1: the pure parts of synthworld's family design and sequence evolver (synthworld_spec §2.1 S3-S4,
§2.3; card F-02 session 3). No ZOMBI2 and no HMMER binary here: the worlds are tested in
tests/r2/test_synthworld_sequences.py and the HMMER parity in tests/r2/test_synthworld_sequences_pfam.py.

Every oracle is written in this file: the genetic code is spelled out codon by codon, the discrete-gamma
category means are checked against sorted Monte Carlo draws, and the profile-F81 closed form is computed
here, not taken from the generator.
"""

from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from biobox.synthworld import backbone, codons, profiles, sequences

pytestmark = pytest.mark.covers("F-02")  # ring R1 (by directory)

AA = "ACDEFGHIKLMNPQRSTVWY"

#: genetic code 11 (NCBI transl_table=11; same amino acids as the standard code), written out per amino acid
CODE11_BY_AA = {
    "F": "TTT TTC", "L": "TTA TTG CTT CTC CTA CTG", "I": "ATT ATC ATA", "M": "ATG", "V": "GTT GTC GTA GTG",
    "S": "TCT TCC TCA TCG AGT AGC", "P": "CCT CCC CCA CCG", "T": "ACT ACC ACA ACG", "A": "GCT GCC GCA GCG",
    "Y": "TAT TAC", "*": "TAA TAG TGA", "H": "CAT CAC", "Q": "CAA CAG", "N": "AAT AAC", "K": "AAA AAG",
    "D": "GAT GAC", "E": "GAA GAG", "C": "TGT TGC", "W": "TGG", "R": "CGT CGC CGA CGG AGA AGG",
    "G": "GGT GGC GGA GGG",
}
CODE11 = {c: aa for aa, cs in CODE11_BY_AA.items() for c in cs.split()}
STARTS11 = {"TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"}


def _cfg(**over) -> dict:
    cfg = copy.deepcopy(backbone.load_preset("micro")["sequences"])
    for path, value in over.items():
        node = cfg
        keys = path.split("__")
        for k in keys[:-1]:
            node = node[k]
        assert keys[-1] in node, path
        node[keys[-1]] = value
    return cfg


# --- genetic code and codon usage (§2.3 "Codons") ---------------------------------------------------------

def test_the_test_table_is_complete():
    assert len(CODE11) == 64


def test_codon_usage_rows_cover_exactly_the_synonymous_codons():
    cfg = _cfg()["codons"]
    tab = codons.usage_table(0.37, cfg["start_codons"], cfg["stop_codons"])
    assert tab.shape == (22, 64)
    assert np.allclose(tab.sum(axis=1), 1.0)
    for i, aa in enumerate(AA):
        support = {codons.CODONS[j] for j in np.nonzero(tab[i] > 0)[0]}
        assert support == set(CODE11_BY_AA[aa].split()), aa
    assert {codons.CODONS[j] for j in np.nonzero(tab[codons.STOP] > 0)[0]} == set(cfg["stop_codons"])
    assert {codons.CODONS[j] for j in np.nonzero(tab[codons.START] > 0)[0]} == set(cfg["start_codons"])
    assert set(cfg["stop_codons"]) <= set(CODE11_BY_AA["*"]) and set(cfg["start_codons"]) <= STARTS11


def test_codon_usage_follows_gc3():
    cfg = _cfg()["codons"]
    for g in (0.2, 0.5, 0.8):
        tab = codons.usage_table(g, cfg["start_codons"], cfg["stop_codons"])
        phe = tab[AA.index("F")]
        assert phe[codons.CODONS.index("TTC")] == pytest.approx(g)          # 2-fold, one GC-ending codon
        ile = tab[AA.index("I")]
        assert ile[codons.CODONS.index("ATC")] == pytest.approx(g / (g + 2 * (1 - g)))
        stop = tab[codons.STOP]
        assert {codons.CODONS[j]: stop[j] for j in np.nonzero(stop)[0]} == pytest.approx(
            {k: v / sum(cfg["stop_codons"].values()) for k, v in cfg["stop_codons"].items()})


def test_translate_agrees_with_the_written_out_code():
    rng = np.random.default_rng(11)
    sense = [c for c, aa in CODE11.items() if aa != "*"]
    for _ in range(200):
        body = [sense[i] for i in rng.integers(len(sense), size=int(rng.integers(1, 60)))]
        start = sorted(STARTS11 & {"ATG", "GTG", "TTG"})[int(rng.integers(3))]
        stop = sorted(CODE11_BY_AA["*"].split())[int(rng.integers(3))]
        cds = start + "".join(body) + stop
        want = "M" + "".join(CODE11[c] for c in body) + "*"
        assert codons.translate(cds) == want


# --- background and profiles (§2.3 "Sites", "Root") -------------------------------------------------------

def test_background_is_hmmers_null_model():
    pyhmmer = pytest.importorskip("pyhmmer")
    bg = np.asarray(pyhmmer.plan7.Background(pyhmmer.easel.Alphabet.amino()).residue_frequencies,
                    dtype=float)[:20]
    assert np.allclose(profiles.BACKGROUND, bg / bg.sum(), atol=1e-6)
    assert profiles.BACKGROUND.sum() == pytest.approx(1.0)


def _toy_profile(M: int, rng, *, mi=0.0, ii=0.0, md=0.0, dd=0.0) -> profiles.Profile:
    match = rng.dirichlet(np.full(20, 0.5), size=M)
    t = np.zeros(M + 1)
    return profiles.Profile(index=0, name="TOY", accession="TOY.1", source="synthetic", match=match,
                            t_mi=t + mi, t_ii=t + ii, t_md=t + md, t_dd=t + dd, ga=None,
                            in_scan_panel=True, promiscuous=False)


def test_root_walk_without_transitions_emits_every_match_state_from_its_row():
    rng = np.random.default_rng(3)
    prof = _toy_profile(5, rng)
    draws = [profiles.sample_domain(prof, np.random.default_rng(1000 + i)) for i in range(4000)]
    assert all(len(aa) == 5 and list(node) == [1, 2, 3, 4, 5] for aa, node in draws)
    counts = np.zeros((5, 20))
    for aa, _ in draws:
        counts[np.arange(5), aa] += 1
    expected = prof.match * len(draws)
    keep = expected >= 5
    chi2 = (((counts - expected) ** 2 / np.where(keep, expected, 1.0)) * keep).sum()
    dof = int(keep.sum()) - 5
    assert chi2 < dof + 5 * math.sqrt(2 * dof), (chi2, dof)   # mean + 5 sd of chi-square


def test_root_walk_inserts_and_deletes_follow_the_transitions():
    rng = np.random.default_rng(4)
    prof = _toy_profile(10, rng)
    t_mi = prof.t_mi.copy()
    t_ii = prof.t_ii.copy()
    t_mi[4], t_ii[4] = 1.0, 0.5            # after match 4 always an insert run, length 1 + Geometric
    t_md = prof.t_md.copy()
    t_dd = prof.t_dd.copy()
    t_md[0], t_dd[1] = 1.0, 1.0            # begin -> D1 -> D2 -> M3 always
    prof = profiles.Profile(**{**prof.__dict__, "t_mi": t_mi, "t_ii": t_ii, "t_md": t_md, "t_dd": t_dd})
    n_ins = []
    for i in range(3000):
        aa, node = profiles.sample_domain(prof, np.random.default_rng(i))
        matches = [int(x) for x in node if x > 0]
        assert matches == list(range(3, 11))
        k = list(node).index(4)
        run = 0
        while node[k + 1 + run] == 0:
            run += 1
        assert run >= 1 and int((node == 0).sum()) == run
        n_ins.append(run)
    # 1 + Geometric(continue 0.5): mean 2, variance 2
    assert abs(np.mean(n_ins) - 2.0) < 5 * math.sqrt(2.0 / len(n_ins))


def test_synthetic_panel_is_deterministic_and_well_formed():
    cfg = _cfg()["profiles"]
    a = profiles.synthetic_panel(5, cfg)
    b = profiles.synthetic_panel(5, cfg)
    c = profiles.synthetic_panel(6, cfg)
    assert len(a) == int(cfg["n_scan"]) + int(cfg["n_dark"])
    assert all(np.array_equal(x.match, y.match) for x, y in zip(a, b))
    assert not np.array_equal(a[0].match, c[0].match)
    assert sum(not p.in_scan_panel for p in a) == int(cfg["n_dark"])
    assert sum(p.promiscuous for p in a) == int(cfg["n_promiscuous"])
    assert all(p.in_scan_panel for p in a if p.promiscuous)
    for p in a:
        assert np.allclose(p.match.sum(axis=1), 1.0) and (p.match >= 0).all()
        lo, hi = cfg["synthetic_match"]["min"], cfg["synthetic_match"]["max"]
        assert lo <= p.M <= hi
        assert len({p.accession for p in a}) == len(a)


# --- rates (§2.3 "Substitution": Γ categories; "Divergence control in closed form") ------------------------

@pytest.mark.parametrize("alpha", [0.5, 1.0, 2.0])
def test_discrete_gamma_category_means_match_sorted_draws(alpha):
    k = 4
    got = sequences.gamma_rates(alpha, k)
    assert got.shape == (k,) and np.all(np.diff(got) > 0)
    assert got.mean() == pytest.approx(1.0, abs=1e-5)
    draws = np.sort(np.random.default_rng(17).gamma(alpha, 1.0 / alpha, size=2_000_000))
    want = draws.reshape(k, -1).mean(axis=1)
    assert np.allclose(got, want, rtol=0.01, atol=0.002), (got, want)


def _closed_form(mu, rates, h, T):
    """§2.3 written out: expected identity mean over sites of e^{-mu r T} + (1 - e^{-mu r T}) h."""
    T = np.atleast_1d(T)[:, None]
    e = np.exp(-mu * rates[None, :] * T)
    return float((e + (1 - e) * h[None, :]).mean())


def test_solve_rate_hits_the_target_identity():
    rng = np.random.default_rng(5)
    rates = sequences.gamma_rates(1.0, 4)[rng.integers(4, size=300)]
    h = rng.uniform(0.06, 0.4, size=300)
    T = rng.uniform(0.2, 3.0, size=40)
    for target in (0.9, 0.6, 0.4):
        mu, expected, saturation = sequences.solve_rate(target, rates, h, T, bounds=(1e-6, 50.0))
        assert expected == pytest.approx(target, abs=1e-4)
        assert _closed_form(mu, rates, h, T) == pytest.approx(target, abs=1e-4)
        assert saturation == pytest.approx(h.mean(), abs=1e-9)
    # a target below the profile's saturation identity is raised to saturation + margin
    mu, expected, saturation = sequences.solve_rate(0.05, rates, h, T, bounds=(1e-6, 50.0))
    assert expected == pytest.approx(h.mean() + sequences.SATURATION_MARGIN, abs=1e-4)


# --- family design and the evolver on a toy tree ------------------------------------------------------------

def _toy_world(M: int, cfg: dict, seed: int = 9):
    rng = np.random.default_rng(seed)
    panel = [profiles.Profile(**{**_toy_profile(M, rng).__dict__, "index": 0})]
    return panel


def _toy_tree(t: float, n_leaves: int = 2):
    nodes = [backbone.GeneTreeNode(0, 0, -1, "speciation", "n0", "", 0.0, False)]
    for i in range(1, n_leaves + 1):
        nodes.append(backbone.GeneTreeNode(0, i, 0, "extant", f"n{i}", f"n{i}_g{i}", t, True))
    branches = [backbone.SpeciesBranch("n0", "", 0.0)] + \
        [backbone.SpeciesBranch(f"n{i}", "n0", t) for i in range(1, n_leaves + 1)]
    return nodes, branches


def test_evolver_identity_matches_the_closed_form_on_a_two_leaf_tree():
    """F81 with site-specific frequencies: two leaves at time t from their root are at path length 2t.
    One long domain (3,000 match sites), no tails, no clock spread; realized identity within 0.03 of the
    closed form (derived: binomial sd over 3,000 sites is below 0.01)."""
    cfg = _cfg(architecture__domain_count_weights=[1.0], architecture__dark_fraction=0.0,
               architecture__promiscuous_fraction=0.0, architecture__n_tail_mean=0,
               architecture__c_tail_mean=0, architecture__linker_mean=0,
               evolution__clock_sigma=0.0, codons__gc3_sigma=0.0)
    panel = _toy_world(3000, cfg)
    design = sequences.design_family(cfg, 21, 0, panel)
    dom = np.isin(design.kind, (sequences.SITE_MATCH, sequences.SITE_INSERT))
    assert dom.sum() == 3000 and len(design.aa) == 3002     # start + 3,000 + stop
    pi = design.pi_table[design.pi_idx[dom]]
    h = (pi ** 2).sum(axis=1)
    r = design.rates[dom]
    for t in (0.05, 0.2, 0.6, 1.5):
        nodes, branches = _toy_tree(t)
        lin = sequences.lineage_model(branches, cfg, seed=21)
        res = sequences.evolve_family(design, nodes, lin, cfg, np.random.default_rng(int(t * 1000)), mu=1.0)
        a, b = (leaf.protein for leaf in res.leaves)
        assert len(a) == len(b) == 3001
        got = np.mean([x == y for x, y in zip(a[1:], b[1:])])
        assert got == pytest.approx(_closed_form(1.0, r, h, 2 * t), abs=0.03), t
        root = res.root_protein
        got_root = np.mean([x == y for x, y in zip(a[1:], root[1:])])
        assert got_root == pytest.approx(_closed_form(1.0, r, h, t), abs=0.03), t
        for leaf in res.leaves:
            assert leaf.subst_per_site == pytest.approx(t)   # mu 1, clock 1


def test_evolver_cds_translates_and_indels_stay_out_of_domains():
    cfg = _cfg(architecture__domain_count_weights=[0.0, 1.0], architecture__dark_fraction=0.0,
               architecture__promiscuous_fraction=0.0, architecture__n_tail_mean=40,
               architecture__c_tail_mean=40, architecture__linker_mean=40, evolution__indel_rate=0.5)
    panel = _toy_world(120, cfg)
    design = sequences.design_family(cfg, 33, 7, panel)
    spans = [(d.start, d.end) for d in design.domains]
    assert len(spans) == 2
    nodes, branches = _toy_tree(1.0, n_leaves=6)
    lin = sequences.lineage_model(branches, cfg, seed=33)
    res = sequences.evolve_family(design, nodes, lin, cfg, np.random.default_rng(1), mu=1.0)
    root_dom = [res.root_protein[s:e] for s, e in spans]
    lengths = set()
    for leaf in res.leaves:
        cds = leaf.cds
        assert len(cds) % 3 == 0 and cds[:3] in STARTS11 and CODE11[cds[-3:]] == "*"
        prot = "M" + "".join(CODE11[cds[i:i + 3]] for i in range(3, len(cds) - 3, 3))
        assert "*" not in prot and prot == leaf.protein
        assert [e - s for s, e in leaf.domain_spans] == [e - s for s, e in spans]
        lengths.add(len(leaf.protein))
        assert leaf.n_indels > 0
        for (s, e), rd in zip(leaf.domain_spans, root_dom):
            assert len(leaf.protein[s:e]) == len(rd)
    assert len(lengths) > 1, "indel rate 0.5 left every tail at the root length"
