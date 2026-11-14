"""Ring R1: the pure parts of synthworld's null worlds (synthworld_spec §4.1; card F-02 session 9). No ZOMBI2 here:
the null worlds of `micro` are built and measured in tests/r2/test_synthworld_nulls.py.

The oracles are written in this file: a recursive walk over two nested presets that lists every leaf that differs,
and a brute-force eligibility loop for N1's profile-preserving relabel (symmetric difference of presence sets,
transfer counts within the tolerance) over a constructed set of families whose matches are known by construction.
"""

from __future__ import annotations

import copy

import numpy as np
import pytest

from biobox.synthworld import backbone, emit, fusions, mobile, nulls

pytestmark = pytest.mark.covers("F-02", "B-NULL")  # ring R1 (by directory)


def _leaf_diff(a, b, where=()):
    """Every path whose leaf differs between two nested dicts (a key on one side only counts as a leaf)."""
    out = set()
    if isinstance(a, dict) and isinstance(b, dict):
        for k in set(a) | set(b):
            if k not in a or k not in b:
                out.add(where + (k,))
            else:
                out |= _leaf_diff(a[k], b[k], where + (k,))
    elif a != b:
        out.add(where)
    return out


# --- N0 presets (§4.1) ----------------------------------------------------------------------------------------

N0_FULL_PATHS = {("modules", "count"), ("profile_only", "count"), ("decoys", "partner_modules"), ("fusions", "nogd"),
                 ("null",)}
N0_STRICT_PATHS = N0_FULL_PATHS | {("vehicles", "cargo"), ("crispr", "arrays")}


@pytest.mark.parametrize("null, paths", [("N0-full", N0_FULL_PATHS), ("N0-strict", N0_STRICT_PATHS)])
def test_n0_preset_differs_from_the_signal_preset_only_where_section_4_1_says(null, paths):
    p = backbone.load_preset("micro")
    q = nulls.null_preset(p, null)
    assert _leaf_diff(p, q) == paths
    assert q["modules"]["count"] == 0 and q["profile_only"]["count"] == 0
    assert q["decoys"]["partner_modules"] == [] and q["fusions"]["nogd"] == 0
    if null == "N0-strict":
        assert q["vehicles"]["cargo"] == [] and q["crispr"]["arrays"] == 0
    else:
        assert q["vehicles"] == p["vehicles"] and q["crispr"] == p["crispr"]
    # the marker names the null world and keeps what the signal preset had, so the truth can say what was removed
    m = q["null"]
    assert m["world"] == null
    assert m["removed"]["modules.count"] == p["modules"]["count"] == 25
    assert m["removed"]["profile_only.count"] == 3 and m["removed"]["fusions.nogd"] == 2
    assert m["removed"]["decoys.partner_modules"] == p["decoys"]["partner_modules"]
    assert ("vehicles.cargo" in m["removed"]) == (null == "N0-strict")
    assert p == backbone.load_preset("micro")          # the signal preset is not touched


@pytest.mark.parametrize("null", ["N0-full", "N0-strict"])
def test_n0_preset_passes_every_stage_guard_and_keeps_the_phage_set(null):
    p = backbone.load_preset("micro")
    q = nulls.null_preset(p, null)
    backbone.check_dec15(q)
    backbone.check_preset(q)
    fusions.check_config(q["fusions"])
    mobile.check_config(q["phages"], q["prophages"], q["crispr"])
    # the phage set is a world of its own from the S1-S5 sections; the null leaves its preset as it was, so S6
    # rebuilds the same phages (spec §4.1: N0 keeps phages and prophages)
    assert mobile.phage_preset(q) == mobile.phage_preset(p)


def test_null_preset_refuses_n1_and_unknown_names_and_a_null_of_a_null():
    p = backbone.load_preset("micro")
    with pytest.raises(ValueError, match="N1"):
        nulls.null_preset(p, "N1")
    with pytest.raises(ValueError):
        nulls.null_preset(p, "N7")
    with pytest.raises(ValueError):
        nulls.null_preset(nulls.null_preset(p, "N0-full"), "N0-strict")


def test_micro_preset_has_the_nulls_section():
    cfg = backbone.load_preset("micro")["nulls"]
    assert cfg["n1"]["max_hamming_fraction"] == 0.1 and cfg["n1"]["transfer_tolerance"] == 0.25
    assert cfg["n1"]["kinds"] == ["root", "originated"] and cfg["n2"]["seeds"] >= 20
    nulls.check_config(cfg)
    for path, value in ((("n1", "max_hamming_fraction"), 1.5), (("n1", "transfer_tolerance"), -0.1),
                        (("n1", "kinds"), ["planted"]), (("n2", "seeds"), 19)):
        bad = copy.deepcopy(cfg)
        bad[path[0]][path[1]] = value
        with pytest.raises(ValueError):
            nulls.check_config(bad)


# --- world ids (§6.2) ------------------------------------------------------------------------------------------

def test_world_id_names_the_null_and_a_null_seed_that_differs_from_the_world_seed():
    assert emit.world_id("micro", 7, None) == "micro-7"
    assert emit.world_id("micro", 7, {"world": "N0-full"}) == "micro-N0-full-7"
    assert emit.world_id("micro", 7, {"world": "N1", "null_seed": 7}) == "micro-N1-7"
    assert emit.world_id("micro", 7, {"world": "N1", "null_seed": 9}) == "micro-N1.9-7"


# --- N2 seeds (§4.1) -------------------------------------------------------------------------------------------

def test_n2_seeds_are_distinct_u64_words_fixed_by_the_master_seed():
    a = nulls.n2_seeds(20261112, 20)
    assert a == nulls.n2_seeds(20261112, 20) and len(set(a)) == 20
    assert all(isinstance(x, int) and 0 <= x < 2**64 for x in a)
    assert 20261112 not in a
    b = nulls.n2_seeds(20261113, 20)
    assert not set(a) & set(b)
    assert nulls.n2_seeds(20261112, 25)[:20] == a


# --- N1 relabel (§4.1, MR-18) ------------------------------------------------------------------------------------

def _eligible(a, profiles, transfers, pool, max_h, tol):
    """The oracle: every free family of the pool within max_h genomes of a's profile and within tol of its
    transfer count (exact integers: |t_x - t_a| <= tol t_a)."""
    out = []
    for x in pool:
        h = len(profiles[a] ^ profiles[x])
        if h <= max_h and abs(transfers[x] - transfers[a]) <= tol * transfers[a] + 1e-9:
            out.append(x)
    return out


def _constructed():
    """20 genomes; module families M0..M5 in two modules; candidates built at a known Hamming distance and
    transfer count from each module family."""
    G = 20
    rng = np.random.default_rng(5)
    profiles, transfers = {}, {}
    modules = [("mod00", ["M0", "M1", "M2"]), ("mod01", ["M3", "M4", "M5"])]
    base = {}
    for k in range(6):
        s = frozenset(int(x) for x in rng.choice(G, size=12, replace=False))
        base[f"M{k}"] = s
        profiles[f"M{k}"] = s
        transfers[f"M{k}"] = 40
    pool = []

    def add(name, src, flips, t):
        s = set(base[src])
        outside = sorted(set(range(G)) - s)
        inside = sorted(s)
        for i in range(flips):           # alternate removing a member and adding a non-member
            if i % 2 == 0:
                s.remove(inside[i // 2])
            else:
                s.add(outside[i // 2])
        profiles[name] = frozenset(s)
        transfers[name] = t
        pool.append(name)

    add("C0a", "M0", 2, 40)      # eligible for M0
    add("C0b", "M0", 2, 50)      # eligible (40 * 1.25 = 50)
    add("C0c", "M0", 3, 40)      # 3 > floor(0.1 * 20) = 2: not eligible
    add("C1a", "M1", 1, 31)      # 31 < 30? no: 40 - 31 = 9 <= 10: eligible
    add("C1b", "M1", 0, 29)      # 11 > 10: not eligible
    add("C2a", "M2", 0, 60)      # not eligible (transfers)
    add("C3a", "M3", 2, 40)
    add("C4a", "M4", 1, 40)
    add("C5a", "M5", 2, 39)
    return G, modules, profiles, transfers, pool


def test_match_relabel_picks_only_eligible_distinct_families_and_marks_the_rest_unmatched():
    G, modules, profiles, transfers, pool = _constructed()
    cfg = {"max_hamming_fraction": 0.1, "transfer_tolerance": 0.25}
    rows = nulls.match_relabel(modules, profiles, transfers, pool, G, cfg, null_seed=11)
    assert [(r.module_id, r.order_index, r.family_id) for r in rows] == \
        [(m, i, f) for m, fs in modules for i, f in enumerate(fs)]
    used = set()
    for r in rows:
        free = [x for x in pool if x not in used]
        ok = _eligible(r.family_id, profiles, transfers, free, 2, 0.25)
        if r.status == "relabeled":
            assert r.replacement_family_id in ok
            used.add(r.replacement_family_id)
            assert r.hamming == len(profiles[r.family_id] ^ profiles[r.replacement_family_id])
            assert (r.n_genomes, r.replacement_n_genomes) == (len(profiles[r.family_id]),
                                                               len(profiles[r.replacement_family_id]))
            assert (r.transfers, r.replacement_transfers) == (transfers[r.family_id],
                                                              transfers[r.replacement_family_id])
        else:
            assert r.status == "unmatched" and not ok and r.replacement_family_id == ""
    by = {r.family_id: r for r in rows}
    assert by["M0"].replacement_family_id in ("C0a", "C0b")
    assert by["M1"].replacement_family_id == "C1a"
    assert by["M2"].status == "unmatched"
    assert {by[f].replacement_family_id for f in ("M3", "M4", "M5")} <= {"C3a", "C4a", "C5a", "C0a", "C0b"}


def test_match_relabel_is_fixed_by_the_null_seed_and_another_seed_can_choose_otherwise():
    G, modules, profiles, transfers, pool = _constructed()
    cfg = {"max_hamming_fraction": 0.1, "transfer_tolerance": 0.25}
    a = nulls.match_relabel(modules, profiles, transfers, pool, G, cfg, null_seed=11)
    assert a == nulls.match_relabel(modules, profiles, transfers, list(reversed(pool)), G, cfg, null_seed=11)
    choices = {nulls.match_relabel(modules, profiles, transfers, pool, G, cfg, null_seed=s)[0].replacement_family_id
               for s in range(1, 40)}
    assert choices == {"C0a", "C0b"}


def test_match_relabel_with_no_transfers_needs_a_family_with_none():
    profiles = {"A": frozenset({0, 1}), "X": frozenset({0, 1}), "Y": frozenset({0, 1})}
    transfers = {"A": 0, "X": 1, "Y": 0}
    rows = nulls.match_relabel([("mod00", ["A"])], profiles, transfers, ["X", "Y"], 10,
                               {"max_hamming_fraction": 0.1, "transfer_tolerance": 0.25}, null_seed=3)
    assert rows[0].replacement_family_id == "Y"
