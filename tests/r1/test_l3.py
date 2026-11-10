"""Ring R1: tapes and pairs against brute-force oracles, and MR-03, MR-04, MR-05 (O-TAPE, O-PAIRS);
the token id layout (card B-09; the token builder itself is tested on micro-S0 in ring R2).

The pair enumeration under test is the pair oracle `pairs_of_contig` (tests/oracles/pairs_oracle.py,
moved from `biobox.l3.pairs` by card B-11), which the catalog's pair macros are compared with in ring R2.

Slots are random (start, end, strand) sets on linear and circular contigs. The oracles below
are written from the definitions in PLAN.md §4.2 and biobox/l3, without the production
helpers: all slot pairs are enumerated and kept by distance, and geometry is recomputed
from coordinates.
"""

from __future__ import annotations

import hashlib
from collections import Counter

import pytest
from hypothesis import assume, given, settings, strategies as st

from biobox.l3.tapes import build_tape, tape_version
from biobox.l3.tokens import token_id, token_level
from oracles.pairs_oracle import FLAG_CENSORED, FLAG_CROSSES_ORIGIN, pairs_of_contig

K = 3


@st.composite
def contigs(draw, circular=None):
    length = draw(st.integers(min_value=200, max_value=5_000))
    circ = draw(st.booleans()) if circular is None else circular
    n = draw(st.integers(min_value=0, max_value=14))
    starts = sorted(draw(st.lists(st.integers(min_value=0, max_value=length - 2), min_size=n, max_size=n, unique=True)))
    slots = []
    for s in starts:
        e = s + draw(st.integers(min_value=1, max_value=400))
        if not circ:
            e = min(e, length)
        else:
            e = min(e, s + length)
        slots.append((s, e, draw(st.sampled_from([1, -1]))))
    if circ:
        assume(all(e <= length for _, e, _ in slots[:-1]))  # only the last slot may wrap
    toks = [draw(st.sampled_from([None, 1, 2, 3])) for _ in slots]
    return slots, length, circ, toks


def _tape(slots, length, circ):
    return {"contig_uid": b"c" * 16, "is_circular": circ, "contig_len": length, **build_tape(slots, length, circ)}


def _tokens(toks):
    return {"family": toks, "nr50": [None] * len(toks), "nr90": [None] * len(toks)}


def _oracle_topology(a, b):
    """a before b along the contig (unwrapped)."""
    if a[2] == b[2]:
        return 1
    if (a[0] <= b[0] and b[1] <= a[1]) or (b[0] <= a[0] and a[1] <= b[1]):
        return 4
    return 2 if a[2] == -1 else 3


def _oracle_pairs(slots, length, circ, toks, k):
    n = len(slots)
    out = []
    for i in range(n):
        for j in range(n):
            if i == j or toks[i] is None or toks[j] is None:
                continue
            def way(d):
                cross = circ and ((d > 0 and j < i) or (d < 0 and j > i))
                left, right = (slots[i], slots[j]) if d > 0 else (slots[j], slots[i])
                if cross:
                    right = (right[0] + length, right[1] + length, right[2])
                return d, cross, left, right, max(left[0], right[0]) - min(left[1], right[1])

            if circ:
                fwd, back = (j - i) % n, (i - j) % n
                if fwd < back:
                    choice = way(fwd)
                elif back < fwd:
                    choice = way(-back)
                else:  # half way round: the smaller gap, forward on a tie
                    a_, b_ = way(fwd), way(-back)
                    choice = b_ if b_[4] < a_[4] else a_
            else:
                choice = way(j - i)
            d, cross, left, right, gap = choice
            if abs(d) > k:
                continue
            censored = not circ and (i < k or n - 1 - i < k)
            out.append((toks[i], toks[j], i, j, d, _oracle_topology(left, right), gap,
                        (FLAG_CENSORED if censored else 0) | (FLAG_CROSSES_ORIGIN if cross else 0)))
    return sorted(out)


def _prod_pairs(slots, length, circ, toks, k):
    rows = pairs_of_contig(_tape(slots, length, circ), _tokens(toks), k)
    return sorted((r["anchor_token"], r["partner_token"], r["anchor_ordinal"], r["partner_ordinal"], r["orf_offset"],
                   r["junction_topology"], r["delta_bp"], r["flags"]) for r in rows)


@settings(max_examples=400, deadline=None)
@given(contigs())
def test_pairs_equal_brute_force(c):
    slots, length, circ, toks = c
    assert _prod_pairs(slots, length, circ, toks, K) == _oracle_pairs(slots, length, circ, toks, K)


@settings(max_examples=400, deadline=None)
@given(contigs())
def test_tape_lists_equal_brute_force(c):
    slots, length, circ, _ = c
    t = build_tape(slots, length, circ)
    n = len(slots)
    assert [len(v) for v in t.values() if v is not None] == [n] * sum(v is not None for v in t.values())  # C-04
    for i in range(n):
        if i == 0 and (not circ or n == 1):
            assert t["delta_bp_prev_slots"][i] is None and t["junction_topology_prevs"][i] is None
            continue
        prev, cur = slots[i - 1], slots[i]
        if i == 0:
            cur = (cur[0] + length, cur[1] + length, cur[2])
        assert t["delta_bp_prev_slots"][i] == max(prev[0], cur[0]) - min(prev[1], cur[1])
        assert t["junction_topology_prevs"][i] == _oracle_topology(prev, cur)
    if circ:
        assert t["orfs_to_contig_starts"] is None and t["orfs_to_contig_ends"] is None
    else:
        assert t["orfs_to_contig_starts"] == list(range(n)) and t["orfs_to_contig_ends"] == list(range(n))[::-1]


@pytest.mark.covers("C-04")
@settings(max_examples=100, deadline=None)
@given(contigs(), st.data())
def test_tape_kinds_are_the_features_kinds(c, data):
    """card B-07: a slot's kind is its feature's kind; the geometry lists do not depend on kinds."""
    slots, length, circ, _ = c
    kinds = data.draw(st.lists(st.sampled_from([1, 2, 3, 4, 15]), min_size=len(slots), max_size=len(slots)))
    t = build_tape(slots, length, circ, kinds)
    assert t["feature_kinds"] == kinds
    assert {k: v for k, v in t.items() if k != "feature_kinds"} == \
        {k: v for k, v in build_tape(slots, length, circ).items() if k != "feature_kinds"}
    assert build_tape(slots, length, circ)["feature_kinds"] == [1] * len(slots)  # default: every slot a CDS
    with pytest.raises(ValueError, match="kinds"):
        build_tape(slots, length, circ, kinds + [1])


def test_tape_version_known_answer():
    # TASKS.md §2: "t-" + the first 12 hex digits of sha256 over the sorted producer names joined by "\n"
    producers = ["pyrodigal-3.6.3-gv-0.3.2-p2", "test-ncrna", "aragorn-1.2.41"]
    want = "t-" + hashlib.sha256("aragorn-1.2.41\npyrodigal-3.6.3-gv-0.3.2-p2\ntest-ncrna".encode("utf-8")).hexdigest()[:12]
    assert tape_version(producers) == want
    assert tape_version(producers[::-1]) == want == tape_version(sorted(producers))
    one = "pyrodigal-3.6.3-gv-0.3.2-p2"
    assert tape_version([one]) == "t-" + hashlib.sha256(one.encode()).hexdigest()[:12]
    assert tape_version([one]) != want
    assert len(want) == 14 and all(ch in "0123456789abcdef" for ch in want[2:])
    for bad in ([], [one, one]):
        with pytest.raises(ValueError):
            tape_version(bad)


def test_token_id_bit_layout_known_answer():
    # PLAN.md §4.1 and card B-09 ("token ids as today"): level (3 bits) | source (5 bits) | local (56 bits)
    assert token_id("family", 1, 0) == 0x2100_0000_0000_0000
    assert token_id("nr90", 0, 5) == 0x6000_0000_0000_0005
    assert token_id("nr50", 31, 2**56 - 1) == 0x5FFF_FFFF_FFFF_FFFF
    assert [token_level(t) for t in (0x2100_0000_0000_0000, 0x4000_0000_0000_0000, 0x6000_0000_0000_0005)] == \
        ["family", "nr50", "nr90"]
    for bad in (("family", 32, 0), ("nr90", 0, 2**56), ("nr50", -1, 0)):
        with pytest.raises(ValueError):
            token_id(*bad)


def _shape(rows, drop_cross=False):
    """Pairs as a multiset, without ordinals (which a transformation relabels)."""
    return Counter((a, b, abs(d), topo, delta, f & ~FLAG_CROSSES_ORIGIN if drop_cross else f)
                   for a, b, _, _, d, topo, delta, f in rows)


@settings(max_examples=300, deadline=None)
@given(contigs())
def test_mr03_reverse_complement_mirrors_pairs(c):
    slots, length, circ, toks = c
    if circ and slots and slots[-1][1] > length:
        return  # a wrapping slot mirrors to one that starts before the origin; MR-04 covers rotation
    # slots are ordered by start; nested slots or shared starts/ends reorder under the mirror
    # (biobox/l3/tapes.py), so MR-03 is stated for slots whose start and end orders agree
    if any(not (a[0] < b[0] and a[1] < b[1]) for a, b in zip(slots, slots[1:])):
        return
    mirrored = sorted(((length - e, length - s, -strand), t) for (s, e, strand), t in zip(slots, toks))
    m_slots, m_toks = [s for s, _ in mirrored], [t for _, t in mirrored]
    before, after = _prod_pairs(slots, length, circ, toks, K), _prod_pairs(m_slots, length, circ, m_toks, K)
    assert _shape(after, drop_cross=True) == _shape(before, drop_cross=True)
    # offsets change sign: partners downstream become upstream
    assert Counter(d for *_, d, _, _, _ in after) == Counter(-d for *_, d, _, _, _ in before) or circ


@settings(max_examples=300, deadline=None)
@given(contigs(circular=True), st.integers(min_value=0, max_value=4_999))
def test_mr04_rotation_keeps_pairs(c, r):
    slots, length, _, toks = c
    r %= length
    rotated = sorted((((s - r) % length, (s - r) % length + (e - s), strand), t) for (s, e, strand), t in zip(slots, toks))
    r_slots, r_toks = [s for s, _ in rotated], [t for _, t in rotated]
    assume(all(e <= length for _, e, _ in r_slots[:-1]))
    before, after = _prod_pairs(slots, length, True, toks, K), _prod_pairs(r_slots, length, True, r_toks, K)
    assert _shape(after, drop_cross=True) == _shape(before, drop_cross=True)


@settings(max_examples=300, deadline=None)
@given(contigs(circular=False), st.data())
def test_mr05_cut_changes_only_windows_across_the_cut(c, data):
    slots, length, _, toks = c
    cuts = [i for i in range(1, len(slots)) if max(e for _, e, _ in slots[:i]) <= slots[i][0]]
    if not cuts:
        return
    cut = data.draw(st.sampled_from(cuts))
    p = slots[cut][0]
    left = _prod_pairs(slots[:cut], p, False, toks[:cut], K)
    right = _prod_pairs([(s - p, e - p, strand) for s, e, strand in slots[cut:]], length - p, False, toks[cut:], K)
    whole = _prod_pairs(slots, length, False, toks, K)
    kept = [r for r in whole if (r[2] < cut) == (r[3] < cut)]
    spans = [r for r in whole if (r[2] < cut) != (r[3] < cut)]
    # only windows crossing the cut change: pairs across it are lost, the rest keep their geometry
    shape = lambda rows: Counter((a, b, d, topo, delta) for a, b, _, _, d, topo, delta, _ in rows)  # noqa: E731
    assert shape(left + right) == shape(kept)
    assert len(whole) == len(kept) + len(spans)
    # and their anchors gain the edge flag when the cut is within K of them
    assert all(f & FLAG_CENSORED for _, _, i, _, _, _, _, f in left if i >= cut - K)
    assert all(f & FLAG_CENSORED for _, _, i, _, _, _, _, f in right if i < K)
