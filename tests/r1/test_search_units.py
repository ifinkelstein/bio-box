"""Ring R1: search units against oracles written from the definitions (S-04, S-09, S-12 lint, S-13).

ORFs: an independent codon table and a per-position definition (an ORF starts after a stop,
or at a linear contig end, and runs to the next stop). Neighborhoods: windows recomputed from
coordinates by enumerating both ways round. Recipes: the paper checks that refuse a recipe.
"""

from __future__ import annotations

import dataclasses
import json
import platform
import resource
import signal
from collections import Counter
from pathlib import Path

import pytest
from hypothesis import assume, given, settings, strategies as st

from biobox import runstep
from biobox.search import export, nbr, orfs, recipe, runner, worker
from biobox.search.recipe import RecipeError

REPO = Path(__file__).resolve().parents[2]

# --- ORFs --------------------------------------------------------------------------------------

_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_B = "TCAG"
TABLE = {a + b + c: _AAS[16 * i + 4 * j + k] for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}


def _table(code):
    t = dict(TABLE)
    if code == 4:
        t["TGA"] = "W"
    elif code == 15:
        t["TAG"] = "Q"
    elif code == 25:
        t["TGA"] = "G"
    return t


def _rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _oracle_orfs(seq, circular, code, min_aa):
    t, L, out = _table(code), len(seq), set()
    for strand in (1, -1):
        s = seq if strand == 1 else _rc(seq)

        def codon(p):
            return "".join(s[(p + x) % L] for x in range(3)) if circular else s[p:p + 3]

        for p in range(L):
            if circular:
                if t[codon((p - 3) % L)] != "*":
                    continue
            elif not (p < 3 or (p >= 3 and t.get(codon(p - 3)) == "*")):
                continue
            aa, q, closed = [], p, False
            while (circular and q - p + 3 <= L) or (not circular and q + 3 <= L):
                c = t[codon(q)]
                if c == "*":
                    closed = True
                    break
                aa.append(c)
                q += 3
            if circular and not closed:
                continue
            if len(aa) < min_aa:
                continue
            lo, hi = p, q + 3 if closed else q
            if strand == 1:
                start, end = lo, hi
            else:
                start, end = L - hi, L - lo
                if circular and start < 0:
                    start, end = start + L, end + L
            partial = not circular and (p < 3 or not closed)  # open at a linear contig end
            out.add((start, end, strand, "".join(aa), partial))
    return out


def _prod_orfs(seq, circular, code, min_aa):
    return {(o.start, o.end, o.strand, o.aa, o.partial)
            for o in orfs.six_frame_orfs(seq.encode(), circular, code, min_aa)}


@settings(max_examples=300, deadline=None)
@given(st.text(alphabet="ACGT", min_size=3, max_size=240), st.booleans(), st.sampled_from([11, 4]))
def test_six_frame_orfs_equal_per_position_oracle(seq, circular, code):
    assert _prod_orfs(seq, circular, code, 3) == _oracle_orfs(seq, circular, code, 3)


@settings(max_examples=200, deadline=None)
@given(st.text(alphabet="ACGT", min_size=30, max_size=300), st.data())
def test_dark_orf_rule_equals_definition(seq, data):
    L = len(seq)
    genes = []
    for _ in range(data.draw(st.integers(0, 3))):
        s = data.draw(st.integers(0, L - 10))
        genes.append((s, data.draw(st.integers(s + 3, L)), data.draw(st.sampled_from([1, -1]))))
    got = {(o.start, o.end, o.strand) for o, _, _ in orfs.dark_orfs(seq.encode(), False, genes, 11, 3, 20, 0.5)}
    covered = set()
    for s, e, _ in genes:
        covered |= set(range(s, e))
    own = {(1, e) for s, e, st_ in genes if st_ == 1} | {(-1, s) for s, e, st_ in genes if st_ == -1}
    want = set()
    for start, end, strand, _, _ in _oracle_orfs(seq, False, 11, 3):
        if (strand, end if strand == 1 else start) in own:
            continue
        frac = len(covered & set(range(start, end))) / (end - start)
        if frac >= 0.5 and not (start < 20 or end > L - 20):
            continue
        want.add((start, end, strand))
    assert got == want


@settings(max_examples=300, deadline=None)
@given(st.text(alphabet="ACGT", min_size=3, max_size=240), st.booleans(), st.sampled_from([4, 11, 15, 25]))
@pytest.mark.covers("W-01")
def test_w01_six_frame_orfs_equal_oracle_for_codes_4_11_15_25(seq, circular, code):
    """Both strands, circular and linear: translation and stop-to-stop boundaries per code."""
    assert _prod_orfs(seq, circular, code, 3) == _oracle_orfs(seq, circular, code, 3)


@settings(max_examples=300, deadline=None)
@given(st.text(alphabet="ACGT", min_size=0, max_size=300), st.sampled_from([4, 11, 15, 25]))
@pytest.mark.covers("W-01")
def test_w01_translate_equals_codon_table(seq, code):
    t = _table(code)
    want = "".join(t[seq[i:i + 3]] for i in range(0, len(seq) - 2, 3))
    assert orfs._translate(seq, code) == want


@pytest.mark.covers("W-01")
def test_w01_tag_is_glutamine_and_not_a_stop_under_code_15():
    body = "GCT" * 40
    seq = "ATG" + body + "TAG" + body + "TAA"
    under11 = {(o.start, o.end, o.strand, o.aa) for o in orfs.six_frame_orfs(seq.encode(), False, 11, 30) if o.strand == 1}
    under15 = {(o.start, o.end, o.strand, o.aa) for o in orfs.six_frame_orfs(seq.encode(), False, 15, 30) if o.strand == 1}
    first = 3 + 120 + 3
    assert (0, first, 1, "M" + "A" * 40) in under11  # the ORF ends at TAG
    assert (first, len(seq), 1, "A" * 40) in under11
    through = (0, len(seq), 1, "M" + "A" * 40 + "Q" + "A" * 40)  # reads through, TAG as Q
    assert through in under15 and through not in under11
    assert (0, first, 1, "M" + "A" * 40) not in under15 and (first, len(seq), 1, "A" * 40) not in under15
    assert orfs._translate("TAGTAATGA", 15) == "Q**"
    assert orfs._translate("TAGTAATGA", 11) == "***"


@pytest.mark.parametrize("code", [0, 7, 26, 99])
@pytest.mark.covers("W-01")
def test_w01_unsupported_code_raises_value_error_naming_it(code):
    with pytest.raises(ValueError, match=str(code)):
        orfs._translate("ATGAAATAG", code)
    with pytest.raises(ValueError, match=str(code)):
        orfs.six_frame_orfs(b"ATGAAATAGATGAAATAG", False, code, 1)
    with pytest.raises(ValueError, match=str(code)):
        orfs.check_code(code)


@pytest.mark.covers("W-01")
def test_w01_supported_codes_include_the_golden_set():
    assert {4, 11, 15, 25} <= set(orfs.supported_codes())


# --- neighborhoods (S-09) ------------------------------------------------------------------------


@st.composite
def contigs(draw, circular=None):
    length = draw(st.integers(min_value=300, max_value=6_000))
    circ = draw(st.booleans()) if circular is None else circular
    n = draw(st.integers(min_value=1, max_value=14))
    starts = sorted(draw(st.lists(st.integers(0, length - 2), min_size=n, max_size=n, unique=True)))
    slots = []
    for s in starts:
        e = min(s + draw(st.integers(1, 500)), length)
        slots.append((s, e, draw(st.sampled_from([1, -1]))))
    return slots, length, circ


def _tape(slots, length, circ):
    return {"contig_uid": b"c" * 16, "contig_len": length, "is_circular": circ,
            "start_bps": [s for s, _, _ in slots], "end_bps": [e for _, e, _ in slots],
            "strands": [x for _, _, x in slots], "protein_uids": [bytes([i]) * 16 for i in range(len(slots))],
            "feature_uids": [bytes([0xF0, i]) * 8 for i in range(len(slots))]}


def _gap(a, b):
    return max(a[0], b[0]) - min(a[1], b[1])


def _oracle_nbr(slots, length, circ, i, p: nbr.Policy):
    n = len(slots)
    a = slots[i]
    if not circ:
        if p.edge == "exclude" and min(a[0], length - a[1]) < p.edge_bp:
            return None
        if p.edge == "require" and (i < p.edge_n or n - 1 - i < p.edge_n):
            return None
    ways = {}
    for j in range(n):
        if j == i:
            continue
        opts = []
        if circ:
            f = (j - i) % n
            for d in (f, f - n):
                if d == 0:
                    continue
                b = slots[j]
                crosses = (d > 0 and j < i) or (d < 0 and j > i)
                if crosses:
                    b = (b[0] + length, b[1] + length, b[2]) if d > 0 else (b[0] - length, b[1] - length, b[2])
                opts.append((d, _gap(a, b), crosses))
        else:
            opts.append((j - i, _gap(a, slots[j]), False))
        if p.window_genes is not None:
            best = min(opts, key=lambda o: (abs(o[0]), -o[0]))
            if abs(best[0]) <= p.window_genes:
                ways[j] = best
        else:
            best = min(opts, key=lambda o: (o[1], -o[0]))
            if best[1] <= p.window_bp:
                ways[j] = best
    operon = None
    if p.operon_gap is not None:
        operon = {i}
        for step in (1, -1):
            cur = i
            while True:
                nxt = (cur + step) % n if circ else cur + step
                if not circ and not 0 <= nxt < n or nxt in operon or slots[nxt][2] != a[2]:
                    break
                x, y = (slots[cur], slots[nxt])
                if circ and ((step == 1 and nxt < cur) or (step == -1 and nxt > cur)):
                    y = (y[0] + step * length, y[1] + step * length, y[2])
                if _gap(x, y) > p.operon_gap:
                    break
                operon.add(nxt)
                cur = nxt
    if circ:
        cens = False
    elif p.window_genes is not None:
        cens = i < p.window_genes or n - 1 - i < p.window_genes
    else:
        cens = a[0] - p.window_bp < 0 or a[1] + p.window_bp > length
    rows = {(i, 0, 0, None, int(cens))}
    for j, (d, g, crosses) in ways.items():
        if p.strand == "same" and slots[j][2] != a[2]:
            continue
        if operon is not None and j not in operon:
            continue
        rows.add((j, d, d * a[2], g, int(cens) | (2 if crosses else 0)))
    return rows


def _prod_nbr(slots, length, circ, i, p):
    got = nbr.neighborhood(_tape(slots, length, circ), i, p)
    return None if got is None else {(r["ordinal"], r["offset"], r["offset_rel"], r["delta_bp"], r["flags"])
                                     for r in got}


POLICIES = st.builds(
    nbr.Policy.parse,
    window=st.sampled_from(["0g", "1g", "3g", "10g", "100bp", "1kb", "3kb"]),
    strand=st.sampled_from(["any", "same"]),
    operon_gap=st.sampled_from([None, 0, 50, 300]),
    edge=st.sampled_from(["truncate", "exclude:200bp", "require:2"]),
    circular=st.sampled_from(["auto", "off"]),
)


@settings(max_examples=500, deadline=None)
@given(contigs(), POLICIES, st.data())
def test_s09_neighborhood_equals_coordinate_oracle(c, policy, data):
    slots, length, circ = c
    i = data.draw(st.integers(0, len(slots) - 1))
    eff_circ = circ and policy.circular == "auto"
    want = _oracle_nbr(slots, length, eff_circ, i, policy)
    assert _prod_nbr(slots, length, circ, i, policy) == want


@settings(max_examples=300, deadline=None)
@given(contigs(circular=False), POLICIES, st.data())
def test_s09_mr03_reverse_complement_keeps_offsets_along_the_anchor(c, policy, data):
    slots, length, _ = c
    assume(all(a[0] < b[0] and a[1] < b[1] for a, b in zip(slots, slots[1:])))
    i = data.draw(st.integers(0, len(slots) - 1))
    m = sorted((length - e, length - s, -x) for s, e, x in slots)
    mi = len(slots) - 1 - i
    before = nbr.neighborhood(_tape(slots, length, False), i, policy)
    after = nbr.neighborhood(_tape(m, length, False), mi, policy)
    key = lambda rows: None if rows is None else Counter(  # noqa: E731
        (r["offset_rel"], r["same_strand"], r["delta_bp"], r["in_operon"], r["flags"]) for r in rows)
    assert key(after) == key(before)


@settings(max_examples=300, deadline=None)
@given(contigs(circular=True), st.integers(0, 5_999), POLICIES, st.data())
def test_s09_mr04_rotation_keeps_circular_neighborhoods(c, r, policy, data):
    slots, length, _ = c
    policy = dataclasses.replace(policy, circular="auto")  # edge policies do not apply on a circle
    r %= length
    i = data.draw(st.integers(0, len(slots) - 1))
    tagged = sorted((((s - r) % length, (s - r) % length + (e - s), x), k) for k, (s, e, x) in enumerate(slots))
    rot = [t for t, _ in tagged]
    assume(all(e <= length for _, e, _ in rot[:-1]))
    ri = [k for _, k in tagged].index(i)
    before = nbr.neighborhood(_tape(slots, length, True), i, policy)
    after = nbr.neighborhood(_tape(rot, length, True), ri, policy)
    n = len(slots)
    # a partner exactly half way round is reached forward; rotation keeps it on the same side only
    # up to that tie, so compare offsets by magnitude there
    key = lambda rows: Counter((abs(x["offset"]) if 2 * abs(x["offset"]) == n else x["offset"], x["delta_bp"],  # noqa: E731
                                x["same_strand"], x["in_operon"]) for x in rows)
    assert key(after) == key(before)


CUT_POLICIES = st.builds(
    nbr.Policy.parse,
    window=st.sampled_from(["1g", "3g", "10g", "100bp", "1kb"]),
    strand=st.sampled_from(["any", "same"]),
    operon_gap=st.sampled_from([None, 0, 50, 300]),
)


@settings(max_examples=300, deadline=None)
@given(contigs(circular=False), CUT_POLICIES, st.data())
def test_s09_mr05_cut_changes_only_windows_across_the_cut(c, policy, data):
    """Cutting a linear contig between genes keeps each side's rows and their geometry, drops the
    rows across the cut, and flags exactly the anchors whose window the cut now reaches."""
    slots, length, _ = c
    cuts = [i for i in range(1, len(slots)) if max(e for _, e, _ in slots[:i]) <= slots[i][0]]
    assume(cuts)
    cut = data.draw(st.sampled_from(cuts))
    p0 = slots[cut][0]
    pieces = [(slots[:cut], p0, 0), ([(s - p0, e - p0, x) for s, e, x in slots[cut:]], length - p0, cut)]
    shape = lambda rows: sorted((r["offset"], r["delta_bp"], r["same_strand"], r["in_operon"]) for r in rows)  # noqa: E731
    for piece, plen, base in pieces:
        for i in range(len(piece)):
            whole = nbr.neighborhood(_tape(slots, length, False), base + i, policy)
            part = nbr.neighborhood(_tape(piece, plen, False), i, policy)
            same_side = [r for r in whole if base <= r["ordinal"] < base + len(piece)]
            assert shape(part) == shape(same_side)
            if policy.window_genes is not None:
                k = policy.window_genes
                assert bool(part[0]["flags"] & nbr.FLAG_CENSORED) == (i < k or len(piece) - 1 - i < k)
            else:
                a = piece[i]
                assert bool(part[0]["flags"] & nbr.FLAG_CENSORED) == (a[0] - policy.window_bp < 0
                                                                      or a[1] + policy.window_bp > plen)


@settings(max_examples=300, deadline=None)
@given(contigs(), POLICIES, st.data())
@pytest.mark.covers("S-09")
def test_b16_neighborhood_rows_name_their_slots_by_feature_uid(c, policy, data):
    """AGENTS §3.3: every row with an ordinal carries the feature_uid at that ordinal and the tape version."""
    slots, length, circ = c
    i = data.draw(st.integers(0, len(slots) - 1))
    tape = {**_tape(slots, length, circ), "tape_version": "t-0123456789ab"}
    rows = nbr.neighborhood(tape, i, policy)
    assume(rows is not None)
    for r in rows:
        assert r["anchor_ordinal"] == i and r["anchor_feature_uid"] == bytes([0xF0, i]) * 8
        assert r["feature_uid"] == bytes([0xF0, r["ordinal"]]) * 8
        assert r["tape_version"] == "t-0123456789ab"
    assert len({r["feature_uid"] for r in rows}) == len(rows)  # no duplicate under the primary key


def test_s09_edge_policies_and_flags_on_a_fixed_contig():
    slots = [(0, 300, 1), (350, 900, 1), (1_000, 1_500, -1), (1_520, 2_000, -1), (5_000, 5_600, 1)]
    t = _tape(slots, 6_000, False)
    p = nbr.Policy.parse(window="2g")
    rows = nbr.neighborhood(t, 0, p)
    assert [r["ordinal"] for r in rows] == [0, 1, 2] and all(r["flags"] & nbr.FLAG_CENSORED for r in rows)
    assert nbr.neighborhood(t, 0, nbr.Policy.parse(window="2g", edge="exclude:3kb")) is None
    assert nbr.neighborhood(t, 2, nbr.Policy.parse(window="2g", edge="require:2")) is not None
    assert nbr.neighborhood(t, 1, nbr.Policy.parse(window="2g", edge="require:2")) is None
    op = nbr.neighborhood(t, 2, nbr.Policy.parse(window="10g", operon_gap=50))
    assert [r["ordinal"] for r in op] == [2, 3] and all(r["in_operon"] for r in op)
    same = nbr.neighborhood(t, 2, nbr.Policy.parse(window="10g", strand="same"))
    # slot 3 lies right of a minus-strand anchor along the contig: upstream of it, so offset_rel -1
    assert [(r["ordinal"], r["offset"], r["offset_rel"]) for r in same] == [(2, 0, 0), (3, 1, -1)]
    circ = nbr.neighborhood(_tape(slots, 6_000, True), 4, nbr.Policy.parse(window="1g"))
    assert {r["ordinal"]: r["flags"] for r in circ} == {3: 0, 4: 0, 0: nbr.FLAG_CROSSES_ORIGIN}
    assert {r["ordinal"]: r["delta_bp"] for r in circ}[0] == 400  # 5_600 -> 6_000 + 0


# --- neighbor summary (W-06) ----------------------------------------------------------------------

C1, C2, G1 = b"\xc1" * 16, b"\xc2" * 16, b"\x91" * 16
P1, P2, P3, P4, P5 = (bytes([0xA0 + i]) * 16 for i in range(1, 6))
A1, A2 = b"\xaa" * 16, b"\xab" * 16
K1, K2 = b"\x31" * 16, b"\x32" * 16


def _nrow(contig, anchor_pu, anchor_ord, anchor_strand, offset, pu, strand, flags=0):
    """One neighborhood row as nbr.fetch writes it; ordinals stand in for slots of a 20-slot tape."""
    fu = lambda o: contig[:8] + bytes([o]) * 8  # noqa: E731
    return {"anchor_protein_uid": anchor_pu, "contig_uid": contig, "anchor_ordinal": anchor_ord,
            "ordinal": anchor_ord + offset, "offset": offset, "offset_rel": offset * anchor_strand, "protein_uid": pu,
            "start": 100 * (anchor_ord + offset), "end": 100 * (anchor_ord + offset) + 90, "strand": strand,
            "delta_bp": None if offset == 0 else 10, "junction_topology": None if offset == 0 else 1,
            "same_strand": strand == anchor_strand, "in_operon": None, "flags": flags, "policy_hash": "p",
            "anchor_feature_uid": fu(anchor_ord), "feature_uid": fu(anchor_ord + offset), "tape_version": "t-x"}


# locus 1 (C1, anchor A1 on +, window not cut): P1 (two domains of F1, one of F2) at +1, P2 (F1) at -2 on the
# other strand, P3 (no domain) at +2; locus 2 (C2, anchor A2 on -, window cut by the contig end): P1 at +1
# along the contig (upstream of the anchor, offset_rel -1), P4 (F3) at -1 along the contig (downstream, +1)
TINY = [_nrow(C1, A1, 5, 1, 0, A1, 1), _nrow(C1, A1, 5, 1, 1, P1, 1), _nrow(C1, A1, 5, 1, -2, P2, -1),
        _nrow(C1, A1, 5, 1, 2, P3, 1),
        _nrow(C2, A2, 1, -1, 0, A2, -1, flags=1), _nrow(C2, A2, 1, -1, 1, P1, -1, flags=1),
        _nrow(C2, A2, 1, -1, -1, P4, 1, flags=1)]
TINY_DOMAINS = [{"protein_uid": P1, "library": "lib-1", "model_acc": "F1", "model_name": "fam1"},
                {"protein_uid": P1, "library": "lib-1", "model_acc": "F1", "model_name": "fam1"},
                {"protein_uid": P1, "library": "lib-1", "model_acc": "F2", "model_name": "fam2"},
                {"protein_uid": P2, "library": "lib-1", "model_acc": "F1", "model_name": "fam1"},
                {"protein_uid": P4, "library": "lib-1", "model_acc": "F3", "model_name": "fam3"},
                {"protein_uid": P5, "library": "lib-1", "model_acc": "F3", "model_name": "fam3"}]
TINY_CLUSTERS = [{"member_uid": u, "level": "nr50", "cluster_uid": c, "rep_uid": r}
                 for u, c, r in ((P1, K1, P2), (P2, K1, P2), (P5, K1, P2), (P3, K2, P3), (A1, K2, P3))]
TINY_LENGTHS = [{"protein_uid": P2, "length": 100}, {"protein_uid": P3, "length": 50}]


def _tiny(by=("pfam", "nr50"), top_n=50, genome_of=None):
    s = nbr.summarize_rows(TINY, by, domains=TINY_DOMAINS, clusters=TINY_CLUSTERS, protein_len=TINY_LENGTHS,
                           genome_of=genome_of, top_n=top_n)
    return s, {(r["level"], r["key"]): r for r in s["rows"]}


@pytest.mark.covers("W-06")
def test_w06_summary_counts_a_multi_domain_neighbor_once_per_family():
    s, rows = _tiny()
    f1, f2, f3 = rows[("pfam", "lib-1:F1")], rows[("pfam", "lib-1:F2")], rows[("pfam", "lib-1:F3")]
    # P1 carries F1 twice and F2 once: one F1 neighbor and one F2 neighbor at each locus it sits in
    assert (f1["n_loci"], f1["n_neighbors"], f1["share_of_loci"]) == (2, 3, 1.0)
    assert (f2["n_loci"], f2["n_neighbors"], f2["share_of_loci"]) == (2, 2, 1.0)
    assert (f3["n_loci"], f3["n_neighbors"], f3["share_of_loci"]) == (1, 1, 0.5)
    assert (f1["library"], f1["model_acc"], f1["model_name"]) == ("lib-1", "F1", "fam1")
    un = rows[("pfam", "")]
    assert un["unannotated"] and (un["n_loci"], un["n_neighbors"], un["offsets"]) == (1, 1, {"2": 1})
    assert [(r["level"], r["key"]) for r in s["rows"] if r["level"] == "pfam"] == [
        ("pfam", "lib-1:F1"), ("pfam", "lib-1:F2"), ("pfam", "lib-1:F3"), ("pfam", "")]


@pytest.mark.covers("W-06")
def test_w06_two_neighbors_of_one_cluster_at_one_locus_are_one_locus_and_two_neighbors():
    s, rows = _tiny()
    k1, k2 = rows[("nr50", K1.hex())], rows[("nr50", K2.hex())]
    assert (k1["n_loci"], k1["n_neighbors"]) == (2, 3)  # P1 and P2 at locus 1, P1 at locus 2
    assert (k1["cluster_uid"], k1["rep_uid"], k1["rep_length"]) == (K1.hex(), P2.hex(), 100)
    assert k1["model_name"] == "fam1"  # members P1 (fam1, fam2), P2 (fam1), P5 (fam3)
    assert (k2["n_loci"], k2["n_neighbors"], k2["rep_length"], k2["model_name"]) == (1, 1, 50, "")
    un = rows[("nr50", "")]  # P4 has no NR50 cluster
    assert (un["n_loci"], un["n_neighbors"], un["offsets"], un["share_of_loci"]) == (1, 1, {"1": 1}, 0.5)
    # every neighbor gene has one cluster or none: the level's neighbors add up
    assert sum(r["n_neighbors"] for r in s["rows"] if r["level"] == "nr50") == s["header"]["n_neighbors"] == 5


@pytest.mark.covers("W-06")
def test_w06_offsets_are_strand_relative_and_the_mode_breaks_ties_toward_the_anchor_then_downstream():
    s, rows = _tiny()
    f1 = rows[("pfam", "lib-1:F1")]
    # locus 1 (+): P1 at +1, P2 at -2; locus 2 (-): P1 one slot right along the contig is upstream, -1
    assert f1["offsets"] == {"-2": 1, "-1": 1, "1": 1} and list(f1["offsets"]) == ["-2", "-1", "1"]
    assert f1["modal_offset"] == 1
    assert f1["same_strand_share"] == round(2 / 3, 4) == 0.6667  # P2 is on the other strand
    f3 = rows[("pfam", "lib-1:F3")]
    assert f3["offsets"] == {"1": 1} and f3["same_strand_share"] == 0.0  # P4 on + beside a - anchor
    k1 = rows[("nr50", K1.hex())]
    assert k1["offsets"] == {"-2": 1, "-1": 1, "1": 1} and k1["modal_offset"] == 1


@pytest.mark.covers("W-06")
def test_w06_top_n_limits_rows_after_totals_and_genomes_come_from_contig_dim():
    s, rows = _tiny(top_n=1)
    h = s["header"]
    assert (h["n_anchor_loci"], h["n_contigs"], h["n_censored_loci"], h["n_anchor_proteins"]) == (2, 2, 1, 2)
    assert h["n_keys"] == {"pfam": 3, "nr50": 2}  # before the cap
    assert [(r["level"], r["key"]) for r in s["rows"]] == [("pfam", "lib-1:F1"), ("pfam", ""),
                                                          ("nr50", K1.hex()), ("nr50", "")]
    assert rows[("pfam", "lib-1:F1")]["share_of_loci"] == 1.0  # shares are of all anchor loci
    assert rows[("pfam", "lib-1:F1")]["n_genomes"] == 2 and h["genome_basis"] == "contig_uid"
    s, rows = _tiny(genome_of=[{"contig_uid": C1, "genome_uid": G1}, {"contig_uid": C2, "genome_uid": G1}])
    assert s["header"]["genome_basis"] == "genome_uid" and rows[("pfam", "lib-1:F1")]["n_genomes"] == 1
    assert all(r["genome_basis"] == "genome_uid" for r in s["rows"])
    # a contig without a dimension row makes the whole summary fall back to contigs
    s, rows = _tiny(genome_of=[{"contig_uid": C1, "genome_uid": G1}])
    assert s["header"]["genome_basis"] == "contig_uid" and rows[("pfam", "lib-1:F1")]["n_genomes"] == 2


_TSV_INT = {"n_loci", "n_neighbors", "n_genomes", "modal_offset", "rep_length"}
_TSV_FLOAT = {"share_of_loci", "same_strand_share"}


def _tsv_cell(col: str, x: str):
    if col in _TSV_INT or col in _TSV_FLOAT:
        return None if x == "" else (int(x) if col in _TSV_INT else float(x))
    if col == "unannotated":
        return {"true": True, "false": False}[x]
    if col == "offsets":
        return {} if x == "" else {o: int(n) for o, n in (p.rsplit(":", 1) for p in x.split(","))}
    return x


def parse_summary_tsv(text: str) -> dict:
    """`bio-box nbr summarize --format tsv` read back: `# key<TAB>json` header lines, then one table per
    level (a column line, then rows), tables separated by a blank line."""
    header, rows, cols = {}, [], None
    for line in text.splitlines():
        if not line:
            cols = None
        elif line.startswith("#"):
            k, tab, v = line[1:].strip(" ").partition("\t")
            if tab:
                header[k] = json.loads(v)
        elif cols is None:
            cols = line.split("\t")
        else:
            rows.append({c: _tsv_cell(c, x) for c, x in zip(cols, line.split("\t"), strict=True)})
    return {"header": header, "rows": rows}


def same_numbers(tsv_rows: list[dict], records: list[dict]) -> bool:
    """Each TSV row equals its NDJSON record on the TSV's columns (an absent name is "" in TSV, null in JSON)."""
    if len(tsv_rows) != len(records):
        return False
    norm = lambda v: "" if v is None else v  # noqa: E731
    return all({k: norm(r[k]) for k in t} == {k: norm(v) for k, v in t.items()} for t, r in zip(tsv_rows, records))


@pytest.mark.covers("W-06")
def test_w06_summary_records_and_tsv_carry_the_same_numbers():
    s, _ = _tiny()
    s["header"].update(release="rX", tape_version="t-x", policy=dataclasses.asdict(nbr.Policy.parse("5g")),
                       policy_hash=nbr.Policy.parse("5g").hash())
    recs = list(export.nbr_summary_records(s))
    assert recs[0]["schema"] == "nbr_summary_header/1" and list(recs[0]) == export.SCHEMAS["nbr_summary_header"]["keys"]
    assert all(list(r) == export.SCHEMAS["nbr_summary"]["keys"] for r in recs[1:]) and len(recs) == 1 + len(s["rows"])
    assert all(r["schema"] == "nbr_summary/1" for r in recs[1:])
    back = parse_summary_tsv("\n".join(export.nbr_summary_tsv(s)))
    assert {k: back["header"][k] for k in ("n_anchor_loci", "n_censored_loci", "release", "tape_version", "policy",
                                           "n_keys")} == {k: recs[0][k] for k in ("n_anchor_loci", "n_censored_loci",
                                                                                 "release", "tape_version", "policy",
                                                                                 "n_keys")}
    assert same_numbers(back["rows"], recs[1:])
    assert {"cluster_uid", "rep_length"} <= set(next(r for r in back["rows"] if r["level"] == "nr50"))
    with pytest.raises(ValueError, match="nr70"):
        nbr.summarize_rows(TINY, ("pfam", "nr70"))


# --- recipes (S-04, S-12 lint) ----------------------------------------------------------------------

BASE = """
recipe: t
written_against: rS0.1
seed: [{sequence: MKVLAAGIVGLLL}]
target: unique
tool: jackhmmer
iterations: {max: 3, inclusion_bits: 27}
controls: {positive: [{sequence: MKVLAAGIVGLLL}]}
"""


def test_s04_recipe_without_positive_control_is_refused():
    with pytest.raises(RecipeError, match="positive control"):
        recipe.load(BASE.replace("controls: {positive: [{sequence: MKVLAAGIVGLLL}]}", ""))
    with pytest.raises(RecipeError, match="positive control"):
        recipe.load(BASE.replace("{positive: [{sequence: MKVLAAGIVGLLL}]}", "{positive: []}"))
    r, _ = recipe.load(BASE)
    assert r["controls"]["negative"] == "shuffled_seed"  # the default negative (MR-09)


@pytest.mark.parametrize("edit,match", [
    (("target: unique", "target: everything"), "target"),
    (("iterations: {max: 3, inclusion_bits: 27}", "iterations: {max: 3, inclusion_evalue: 0.001}"), "bits"),
    (("iterations: {max: 3, inclusion_bits: 27}", "iterations: {max: 3}"), "inclusion_bits"),
    (("tool: jackhmmer", "tool: phmmer"), "one round"),
    (("recipe: t", "recipe: t\nz: 2e7"), "release_default"),
    (("recipe: t", "recipe: t\ncolour: red"), "unknown keys"),
    (("seed: [{sequence: MKVLAAGIVGLLL}]", "seed: [{protein_uid: ABC}]"), "hex"),
    (("seed: [{sequence: MKVLAAGIVGLLL}]", "seed: [{kind: sequence, value: MKV, source: {paper_id: x}}]"), "evidence"),
    (("recipe: t", "recipe: t\ncontext: {edge: sideways}"), "edge"),
])
def test_recipe_paper_checks_refuse(edit, match):
    with pytest.raises(RecipeError, match=match):
        recipe.load(BASE.replace(*edit))


UID32 = "ab" * 16


@pytest.mark.parametrize("kind,value", [("protein_uid", UID32), ("region_uid", UID32), ("sequence", "MKVLAAGIVGLLL"),
                                        ("protein_accession", "MJ0613")])
def test_w08_shorthand_with_label_equals_the_long_form(kind, value):
    s, _ = recipe.load(BASE.replace("seed: [{sequence: MKVLAAGIVGLLL}]", "seed: [{%s: %s, label: x}]" % (kind, value)))
    assert s["seed"][0] == {"kind": kind, "value": value, "label": "x"}
    if kind not in recipe.RUN_POSITIVE_KINDS:  # card W-18: a positive of this kind is refused on paper
        return
    short = BASE.replace("{positive: [{sequence: MKVLAAGIVGLLL}]}",
                         "{positive: [{%s: %s, label: T7 gp1}]}" % (kind, value))
    long = BASE.replace("{positive: [{sequence: MKVLAAGIVGLLL}]}",
                        "{positive: [{kind: %s, value: %s, label: T7 gp1}]}" % (kind, value))
    a, _ = recipe.load(short)
    b, _ = recipe.load(long)
    assert a == b and a["controls"]["positive"][0]["label"] == "T7 gp1"
    assert recipe.recipe_hash(a) == recipe.recipe_hash(b)
    assert s["seed"][0]["value"] == a["controls"]["positive"][0]["value"]


@pytest.mark.covers("W-18")
@pytest.mark.parametrize("form", ["{protein_accession: MJ0613, label: x}", "{kind: protein_accession, value: MJ0613}"])
def test_w18_a_positive_control_the_runner_cannot_resolve_is_refused_by_the_parser(form):
    """Card W-18 (2): the parser and the submit path agree. A `protein_accession` positive control was valid on paper
    and refused at submit ("not executable"): the release records no protein accession to resolve it through, so the
    recipe parser refuses it, naming the kinds that run and the command that finds a protein's uid."""
    text = BASE.replace("{positive: [{sequence: MKVLAAGIVGLLL}]}", "{positive: [%s]}" % form)
    with pytest.raises(RecipeError, match="protein_accession") as e:
        recipe.load(text)
    assert "protein_uid" in str(e.value) and "bio-box find" in str(e.value)
    for kind in recipe.RUN_POSITIVE_KINDS:  # and every kind the runner takes is still accepted
        value = "MKVLAAGIVGLLL" if kind == "sequence" else UID32
        r, _ = recipe.load(BASE.replace("{positive: [{sequence: MKVLAAGIVGLLL}]}", "{positive: [{%s: %s}]}" % (kind, value)))
        assert recipe.executable(r) == []


def test_w08_two_kind_keys_in_one_entry_are_refused():
    with pytest.raises(RecipeError, match="more than one kind"):
        recipe.load(BASE.replace("{positive: [{sequence: MKVLAAGIVGLLL}]}",
                                 "{positive: [{protein_uid: %s, sequence: MKV, label: x}]}" % UID32))
    with pytest.raises(RecipeError, match="unknown keys"):
        recipe.load(BASE.replace("{positive: [{sequence: MKVLAAGIVGLLL}]}",
                                 "{positive: [{sequence: MKVLAAGIVGLLL, colour: red}]}"))


def test_recipe_hash_is_spelling_free_and_gaps_are_listed():
    a, _ = recipe.load(BASE)
    b, _ = recipe.load(BASE.replace("seed: [{sequence: MKVLAAGIVGLLL}]",
                                    "seed:\n  - kind: sequence\n    value: mkvl aagi vglll"))
    assert recipe.recipe_hash(a) == recipe.recipe_hash(b)
    assert recipe.executable(a) == []
    c, _ = recipe.load(BASE.replace("target: unique", "target: igr").replace("tool: jackhmmer", "tool: nhmmer"))
    assert any("target igr" in g for g in recipe.executable(c)) and any("nhmmer" in g for g in recipe.executable(c))


@pytest.mark.covers("W-02")
@pytest.mark.parametrize("kind,extra", [("hmm_file", ""), ("hmm_file", ", model: PF00940"), ("sequences", ""),
                                        ("msa", "")])
def test_w02_file_seed_kinds_parse_on_paper_and_are_executable(kind, extra):
    short = BASE.replace("seed: [{sequence: MKVLAAGIVGLLL}]", "seed: [{%s: /data/x.file%s, label: s}]" % (kind, extra))
    long = BASE.replace("seed: [{sequence: MKVLAAGIVGLLL}]",
                        "seed: [{kind: %s, value: /data/x.file%s, label: s}]" % (kind, extra))
    a, _ = recipe.load(short)
    b, _ = recipe.load(long)
    assert a == b and a["seed"][0]["kind"] == kind and a["seed"][0]["value"] == "/data/x.file"
    assert recipe.executable(a) == []
    if extra:
        assert a["seed"][0]["model"] == "PF00940"


@pytest.mark.covers("W-02")
def test_w02_model_key_belongs_to_hmm_file_only():
    with pytest.raises(RecipeError, match="model"):
        recipe.load(BASE.replace("seed: [{sequence: MKVLAAGIVGLLL}]", "seed: [{sequences: /x.faa, model: PF00940}]"))
    with pytest.raises(RecipeError, match="model"):
        recipe.load(BASE.replace("seed: [{sequence: MKVLAAGIVGLLL}]", "seed: [{hmm: PF00940, model: PF00940}]"))


@pytest.mark.covers("W-02")
@pytest.mark.parametrize("kind", ["hmm_file", "msa", "sequences"])
def test_w02_recipe_hash_covers_file_content_not_path(kind):
    def pinned(path, sha):
        r, _ = recipe.load(BASE.replace("seed: [{sequence: MKVLAAGIVGLLL}]", "seed: [{%s: %s}]" % (kind, path)))
        with pytest.raises(RecipeError, match="sha256"):  # a file seed is hashed only once its content is pinned
            recipe.recipe_hash(r)
        r["seed"][0]["sha256"] = sha
        return r

    a = pinned("/one/place/m.file", "aa" * 32)
    b = pinned("/another/m2.file", "aa" * 32)
    c = pinned("/one/place/m.file", "ab" * 32)
    assert recipe.recipe_hash(a) == recipe.recipe_hash(b) != recipe.recipe_hash(c)
    form = json.dumps(recipe.identity(a))
    assert "/one/place" not in form and "m.file" not in form and "aa" * 32 in form


# --- JSON exports (S-13) ------------------------------------------------------------------------------


def test_s13_json_schemas_match_the_lock():
    locked = json.loads(export.lock_path(REPO).read_text())
    assert export.lock_problems(locked) == []
    removed = {k: {**v, "keys": v["keys"][:-1]} for k, v in export.SCHEMAS.items()}
    assert export.lock_problems(locked, removed)  # a removed key without a bump fails
    bumped = {k: {"version": v["version"] + 1, "keys": v["keys"][:-1]} for k, v in export.SCHEMAS.items()}
    assert export.lock_problems(locked, bumped) == []
    added = {k: {**v, "keys": v["keys"] + ["extra"]} for k, v in export.SCHEMAS.items()}
    assert export.lock_problems(locked, added) == []


ORDINAL_KEYS = {"slot_ordinal", "anchor_ordinal", "ordinal"}


@pytest.mark.covers("S-13")
def test_b16_json_records_with_an_ordinal_carry_feature_uid_and_tape_version():
    """AGENTS §3.3: no exported ordinal without the feature_uid beside it and the tape version it belongs to."""
    tv, cu, anchor, fu, pu = "t-0123456789ab", b"\x01" * 16, b"\x02" * 16, b"\x03" * 16, b"\x04" * 16
    nb = {"anchor_protein_uid": pu, "contig_uid": cu, "anchor_ordinal": 4, "ordinal": 5, "offset": 1, "offset_rel": 1,
          "protein_uid": pu, "start": 10, "end": 20, "strand": 1, "delta_bp": 3, "junction_topology": 1,
          "same_strand": True, "in_operon": None, "flags": 2, "policy_hash": "p", "anchor_feature_uid": anchor,
          "feature_uid": fu, "tape_version": tv}
    [rec] = export.nbr_records([nb])
    assert (rec["anchor_feature_uid"], rec["feature_uid"], rec["tape_version"]) == (anchor.hex(), fu.hex(), tv)
    assert (rec["anchor_ordinal"], rec["ordinal"]) == (4, 5)
    hit = {"target_uid": pu, "target_kind": "protein", "seq_bitscore": 50.0, "bitscore": 40.0}
    occ = {"contig_uid": cu, "ordinal": 5, "feature_uid": fu, "start": 10, "end": 20, "strand": 1}  # card B-17
    [rec] = export.hit_records([hit], {pu: [occ]}, tape_version=tv)
    assert rec["tape_version"] == tv
    assert [(o["ordinal"], o["feature_uid"]) for o in rec["occurrences"]] == [(5, fu.hex())]
    # every record kind whose keys (or nested occurrences) name an ordinal names the tape version too
    for kind, spec in export.SCHEMAS.items():
        keys = set(spec["keys"])
        if keys & ORDINAL_KEYS:
            assert {"feature_uid", "tape_version"} <= keys, kind
            assert "anchor_ordinal" not in keys or "anchor_feature_uid" in keys, kind
        if "occurrences" in keys:
            assert "tape_version" in keys, kind


# --- C-05: the worker's scope, its limits and how an exit is classified ------------------------------


def _found(name: str) -> str:
    return f"/usr/bin/{name}"


def _no_systemd(name: str) -> str | None:
    return None if name == "systemd-run" else f"/usr/bin/{name}"


@pytest.mark.covers("C-05", "S-06")
def test_c05_worker_argv_in_a_scope_and_in_the_fallback():
    job = Path("/w/run/job.json")
    worker_cmd = runner.worker_command(job)
    assert worker_cmd[-1] == str(job)
    lim = runner.Limits(cpus=2, mem_gb=18, nice=15)
    assert runner.worker_argv(job, lim, system="Linux", which=_found) == (
        ["systemd-run", "--user", "--scope", "--quiet", "-p", "MemoryMax=18000000000", "-p", "MemorySwapMax=0",
         "-p", "CPUQuota=200%", "nice", "-n", "15", *worker_cmd], "scope")
    # no memory limit in the run's limits: no MemoryMax (as the fallback sets no RLIMIT_AS then)
    assert runner.worker_argv(job, runner.Limits(cpus=3, nice=15), system="Linux", which=_found) == (
        ["systemd-run", "--user", "--scope", "--quiet", "-p", "MemorySwapMax=0", "-p", "CPUQuota=300%",
         "nice", "-n", "15", *worker_cmd],
        "scope")
    nice = (["nice", "-n", "15", *worker_cmd], "nice")
    assert runner.worker_argv(job, lim, system="Linux", which=_no_systemd) == nice
    assert runner.worker_argv(job, lim, system="Darwin", which=_found) == nice
    assert runner.worker_argv(job, lim, system="Linux", which=lambda name: None) == (worker_cmd, "bare")
    # run-step builds its argv with the same helper
    cmd = ["mmseqs", "easy-search"]
    for system, which in (("Linux", _found), ("Linux", _no_systemd), ("Darwin", _found)):
        assert runstep.build_argv(cmd, 19, "200G", system=system, which=which) == \
            runstep.scope_argv(cmd, 19, mem_max="200G", system=system, which=which)[0]


@pytest.mark.covers("C-05", "S-06")
def test_c05_rlimit_as_only_where_no_scope_caps_memory():
    """RLIMIT_CPU caps each process in both paths; RLIMIT_AS (virtual memory) only without a scope, on Linux."""
    lim = runner.Limits(cpu_min=2, mem_gb=4)
    cpu = {resource.RLIMIT_CPU: (120, 125)}
    assert runner.rlimits(lim, "scope", "Linux") == cpu
    for kind in ("nice", "bare"):
        assert runner.rlimits(lim, kind, "Linux") == {**cpu, resource.RLIMIT_AS: (4 * 10**9, 4 * 10**9)}
        assert runner.rlimits(lim, kind, "Darwin") == cpu
    assert runner.rlimits(runner.Limits(cpu_min=2), "nice", "Linux") == cpu


@pytest.mark.covers("C-05", "S-05")
def test_c05_worker_env_names_the_memory_limit_it_runs_under():
    lim = runner.Limits(mem_gb=4)
    scope = runner.worker_env(lim, "scope", "Linux", {"PATH": "/bin"})
    assert scope["BIOBOX_MEMORY_MAX"] == str(4 * 10**9) and "BIOBOX_RLIMIT_AS" not in scope
    assert scope["PATH"] == "/bin"
    fallback = runner.worker_env(lim, "nice", "Linux", {})
    assert fallback["BIOBOX_RLIMIT_AS"] == str(4 * 10**9) and "BIOBOX_MEMORY_MAX" not in fallback
    for env in (runner.worker_env(lim, "nice", "Darwin", {}), runner.worker_env(runner.Limits(), "scope", "Linux", {})):
        assert "BIOBOX_RLIMIT_AS" not in env and "BIOBOX_MEMORY_MAX" not in env


@pytest.mark.covers("W-04")
def test_w04_runner_caps_the_workers_thread_pools_by_its_cpus():
    want = {"OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "POLARS_MAX_THREADS", "ARROW_IO_THREADS"}
    assert set(worker.THREAD_VARS) == want
    for kind, system in (("scope", "Linux"), ("nice", "Linux"), ("nice", "Darwin"), ("bare", "Linux")):
        env = runner.worker_env(runner.Limits(cpus=3), kind, system, {"OMP_NUM_THREADS": "112"})
        assert {v: env[v] for v in want} == {v: "3" for v in want}
    # the worker's entry sets them itself before it imports anything that starts a pool
    boot = runner.worker_command(Path("/w/job.json"))[2]
    assert boot.index("OMP_NUM_THREADS") < boot.index("import worker")


@pytest.mark.covers("W-09")
def test_w09_every_pyhmmer_search_takes_cpus_from_the_caller(tmp_path, monkeypatch):
    import pyhmmer

    from biobox.search import engine

    calls = []

    def recorder(*args, **kwargs):
        calls.append(kwargs)
        assert kwargs.get("cpus", 0) >= 1, f"pyhmmer search called without a positive cpus: {sorted(kwargs)}"
        return []

    for name in ("hmmsearch", "hmmscan", "phmmer", "jackhmmer", "nhmmer"):
        if hasattr(pyhmmer, name):
            monkeypatch.setattr(pyhmmer, name, recorder)
        if hasattr(pyhmmer.hmmer, name):
            monkeypatch.setattr(pyhmmer.hmmer, name, recorder)
    hmm = engine.build_from_sequence("w09", "MKVLAAGIVGLLLAQPSTEWRHDNCY")
    lib = tmp_path / "lib.hmm"
    with open(lib, "wb") as fh:
        hmm.write(fh)
    target = engine.load_target("protein", iter([("00" * 16, "MKVLAAGIVGLLLAQPSTEWRHDNCY")]))

    engine.families_of("MKVLAAGIVGLLLAQPSTEWRHDNCY", [lib], threads=3)
    engine.search(hmm, target, 25.0, 25.0, threads=5)
    assert [c["cpus"] for c in calls] == [3, 5]
    # no default: a caller that forgets the thread count fails loudly
    with pytest.raises(TypeError):
        engine.families_of("MKVLAAGIVGLLLAQPSTEWRHDNCY", [lib])
    with pytest.raises(TypeError):
        engine.search(hmm, target, 25.0, 25.0)


KILL, TERM, XCPU, SEGV = -signal.SIGKILL, -signal.SIGTERM, -signal.SIGXCPU, -signal.SIGSEGV
OK_RESULT = {"status": "ok", "error": None}


@pytest.mark.covers("C-05", "S-05")
@pytest.mark.parametrize("rc,timed_out,cancelled,cgroup_oom,result,want", [
    (0, False, False, None, OK_RESULT, "ok"),
    (0, False, False, None, {"status": "failed", "error": "positive control not recovered"}, "failed"),
    (0, False, False, None, None, "failed"),                    # exit 0 without a result
    (worker.EXIT_OOM, False, False, None, None, "oom"),          # MemoryError (under RLIMIT_AS in the fallback)
    (KILL, True, False, None, None, "timeout"),                  # the runner's own wall-clock SIGKILL
    (KILL, True, False, True, None, "timeout"),
    (KILL, False, False, True, None, "oom"),                     # the scope's memory.events counts an oom_kill
    (KILL, False, False, False, None, "failed"),                 # the scope says no OOM kill: killed from outside
    (KILL, False, False, None, None, "failed"),                  # W-04: no evidence of OOM: killed from outside
    (worker.EXIT_FAILED, False, False, None, None, "failed"),    # W-04: an uncaught exception in the worker
    (TERM, False, True, None, None, "stopped"),                  # cancel
    (XCPU, False, False, None, None, "timeout"),                 # RLIMIT_CPU
    (SEGV, False, False, None, None, "failed"),
    (1, False, False, None, None, "failed"),
])
def test_c05_exit_classification_from_its_evidence(rc, timed_out, cancelled, cgroup_oom, result, want):
    status, error = runner.classify_exit(rc, timed_out=timed_out, cancelled=cancelled, cgroup_oom=cgroup_oom,
                                         result=result, stderr="Traceback\nRuntimeError: boom\n")
    assert status == want, error
    if want != "ok":
        assert error
    if rc in (1, worker.EXIT_FAILED):
        assert error.endswith("RuntimeError: boom")
    if rc == KILL and not timed_out and not cgroup_oom:
        assert "signal 9" in error and "outside" in error


CG = "/user.slice/user-1000.slice/user@1000.service/app.slice/run-r0123abcd.scope"


def _systemctl(tmp_path: Path, load: str, result: str) -> str:
    """A stand-in for `systemctl` that answers `show` with this LoadState and Result (and ignores the rest)."""
    p = tmp_path / f"systemctl-{load}-{result}"
    p.write_text(f'#!/bin/sh\n[ "$2" = show ] && printf "LoadState={load}\\nResult={result}\\n"\nexit 0\n')
    p.chmod(0o755)
    return str(p)


@pytest.mark.covers("C-05", "S-05")
def test_c05_oom_evidence_from_the_scope(tmp_path):
    work, cg_root = tmp_path / "work", tmp_path / "cgroup"
    work.mkdir()
    assert runner.scope_oom_killed(work, cgroup_root=cg_root, systemctl=None) is None  # the worker left no marker
    (work / runner.CGROUP_FILE).write_text(CG + "\n")
    assert runner.scope_oom_killed(work, cgroup_root=cg_root, systemctl=None) is None  # scope gone, nobody to ask
    events = cg_root / CG.lstrip("/") / "memory.events"
    events.parent.mkdir(parents=True)
    events.write_text("low 0\nhigh 0\nmax 1450\noom 3\noom_kill 1\noom_group_kill 0\n")
    assert runner.scope_oom_killed(work, cgroup_root=cg_root, systemctl=None) is True
    events.write_text("low 0\nhigh 0\nmax 0\noom 0\noom_kill 0\noom_group_kill 0\n")
    assert runner.scope_oom_killed(work, cgroup_root=cg_root, systemctl=None) is False
    events.unlink()  # systemd removed the scope's cgroup: its unit's result says why it ended
    for load, result, want in (("loaded", "oom-kill", True), ("loaded", "success", False),
                               ("loaded", "signal", False), ("not-found", "success", None)):
        got = runner.scope_oom_killed(work, cgroup_root=cg_root, systemctl=_systemctl(tmp_path, load, result))
        assert got is want, (load, result)


# --- A-07: the oom fault reaches its limit on every path; the runner's liveness ---------------------------


def _oom_worker(tmp_path: Path, lim: runner.Limits, kind: str, cap_bytes: int | None = None) -> int:
    """Run the real worker entry with the `oom` fault under the env and rlimits the runner gives `kind`
    on Linux (or on this system for "bare" without mem_gb); its exit code."""
    import os
    import subprocess

    job = tmp_path / kind / "job.json"
    job.parent.mkdir(parents=True)
    job.write_text(json.dumps({"cpus": 1, "fault": "oom"}))
    system = "Linux" if lim.mem_gb is not None else platform.system()
    env = runner.worker_env(lim, kind, system, dict(os.environ), cap_bytes=cap_bytes)
    lims = runner.rlimits(lim, kind, system, cap_bytes)

    def pre():
        for which, value in lims.items():
            resource.setrlimit(which, value)
    return subprocess.run(runner.worker_command(job), env=env, preexec_fn=pre,
                          stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=120).returncode


def _held_by(lim: runner.Limits, kind: str, system: str, cap_bytes: int | None) -> tuple[str, int | None]:
    """What holds the worker's memory on a path, and at how many bytes: ("rlimit_as", n), ("cgroup", n)
    (the scope's MemoryMax, or the inherited unit's cap) or ("none", None)."""
    rl = runner.rlimits(lim, kind, system, cap_bytes)
    if resource.RLIMIT_AS in rl:
        return "rlimit_as", rl[resource.RLIMIT_AS][0]
    if kind == "scope" and lim.mem_gb is not None:
        argv, _ = runner.worker_argv(Path("/w/job.json"), lim, system=system, which=_found)
        return "cgroup", int(next(a for a in argv if a.startswith("MemoryMax=")).split("=")[1])
    if kind == "inherited" and cap_bytes is not None:
        return "cgroup", cap_bytes
    return "none", None


# (kind, system, limits, inherited cap): every way the runner can start a worker
OOM_PATHS = [
    ("scope", "Linux", runner.Limits(mem_gb=4), None),                  # GitHub's ubuntu-24.04 runner
    ("inherited", "Linux", runner.Limits(mem_gb=4), 32 * 2**30),        # the pod under scripts/pod_pytest.sh
    ("inherited", "Linux", runner.Limits(mem_gb=None), 32 * 2**30),     # the unit's cap is the only limit
    ("nice", "Linux", runner.Limits(mem_gb=4), None),
    ("bare", "Linux", runner.Limits(mem_gb=4), None),
    ("nice", "Darwin", runner.Limits(mem_gb=4), None),
    ("scope", "Linux", runner.Limits(mem_gb=None), None),
]


@pytest.mark.covers("A-07", "S-05")
@pytest.mark.parametrize("kind,system,lim,cap", OOM_PATHS)
def test_a07_the_oom_fault_allocates_past_the_limit_that_holds_its_path(kind, system, lim, cap):
    """By construction, not by speed: in one allocation the fault asks for twice the limit the runner named,
    and the limit it names is the one that holds the worker on that path."""
    env = runner.worker_env(lim, kind, system, {}, cap_bytes=cap)
    how, limit = _held_by(lim, kind, system, cap)
    n = worker.oom_bytes(env)
    if how == "none":
        assert n is None
    else:
        assert n == 2 * limit > limit, (how, limit, env)


@pytest.mark.covers("A-07", "S-05")
@pytest.mark.parametrize("kind,system,lim,cap", [p for p in OOM_PATHS if _held_by(p[2], p[0], p[1], p[3])[0] == "cgroup"])
def test_a07_a_cgroup_kill_of_the_oom_fault_is_oom(tmp_path, kind, system, lim, cap):
    """Scope and inherited cap: the kernel answers the allocation with SIGKILL and one more oom_kill in the
    cgroup's memory.events (in an inherited unit, above the count when the worker started)."""
    work, cg_root = tmp_path / "work", tmp_path / "cgroup"
    work.mkdir()
    (work / runner.CGROUP_FILE).write_text(CG + "\n")
    events = cg_root / CG.lstrip("/") / "memory.events"
    events.parent.mkdir(parents=True)
    baseline = 2 if kind == "inherited" else 0  # kills in the shared unit before this worker started
    events.write_text(f"low 0\nhigh 0\nmax 9\noom 9\noom_kill {baseline + 1}\noom_group_kill 0\n")
    evidence = runner.scope_oom_killed(work, cgroup_root=cg_root, systemctl=None, baseline=baseline)
    status, error = runner.classify_exit(KILL, timed_out=False, cancelled=False, cgroup_oom=evidence, result=None)
    assert status == "oom", error
    events.write_text(f"low 0\nhigh 0\nmax 9\noom 9\noom_kill {baseline}\noom_group_kill 0\n")
    evidence = runner.scope_oom_killed(work, cgroup_root=cg_root, systemctl=None, baseline=baseline)
    assert runner.classify_exit(KILL, timed_out=False, cancelled=False, cgroup_oom=evidence, result=None)[0] == "failed"


@pytest.mark.covers("A-07", "S-05")
@pytest.mark.skipif(platform.system() != "Linux", reason="RLIMIT_AS is set only on Linux")
@pytest.mark.parametrize("kind,cap", [("nice", None), ("bare", None), ("inherited", 32 * 2**30)])
def test_a07_the_oom_fault_under_rlimit_as_exits_oom_at_once(tmp_path, kind, cap):
    """RLIMIT_AS paths (no cgroup holds the worker to mem_gb): the real worker's allocation fails in the
    allocator, before any page is touched, and the exit is EXIT_OOM."""
    lim = runner.Limits(mem_gb=4)
    assert _held_by(lim, kind, "Linux", cap)[0] == "rlimit_as"
    rc = _oom_worker(tmp_path, lim, kind, cap)
    assert runner.classify_exit(rc, timed_out=False, cancelled=False, cgroup_oom=None, result=None) == \
        ("oom", "out of memory (MemoryError)"), rc


@pytest.mark.covers("A-07", "S-05")
def test_a07_the_oom_fault_without_a_limit_exits_oom(tmp_path):
    """No limit named (macOS, or no mem_gb): the fault raises the MemoryError an allocation would."""
    lim = runner.Limits(mem_gb=None)
    assert _held_by(lim, "bare", platform.system(), None)[0] == "none"
    rc = _oom_worker(tmp_path, lim, "bare")
    assert runner.classify_exit(rc, timed_out=False, cancelled=False, cgroup_oom=None, result=None)[0] == "oom"


@pytest.mark.covers("A-07", "S-05")
def test_a07_an_oom_fault_that_is_not_stopped_is_failed_not_oom(tmp_path):
    """A named limit that does not hold (a scope without the memory controller, say): the allocation
    succeeds, and the run must not be recorded as `oom`."""
    import os
    import subprocess

    job = tmp_path / "job.json"
    job.write_text(json.dumps({"cpus": 1, "fault": "oom"}))
    env = {**os.environ, "BIOBOX_MEMORY_MAX": str(2**20)}  # 1 MiB, enforced by nothing
    env.pop("BIOBOX_RLIMIT_AS", None)
    p = subprocess.run(runner.worker_command(job), env=env, stdout=subprocess.DEVNULL,
                       stderr=subprocess.PIPE, text=True, timeout=120)
    status, error = runner.classify_exit(p.returncode, timed_out=False, cancelled=False, cgroup_oom=None,
                                         result=None, stderr=p.stderr)
    assert status == "failed" and "not enforced" in error, (p.returncode, p.stderr[-500:])


def _stat_line(pid: int, state: str, flags: int = 0x400100, start: int = 12345) -> str:
    rest = [str(1 + i) for i in range(49)]  # fields 4..52 of proc(5)
    rest[5] = str(flags)                      # field 9
    rest[18] = str(start)                     # field 22
    return f"{pid} (py (thon) x) {state} " + " ".join(rest) + "\n"


def _fake_proc(tmp_path: Path, procs: dict[int, str]) -> Path:
    root = tmp_path / "proc"
    (root / "self").mkdir(parents=True)
    (root / "self" / "stat").write_text(_stat_line(1, "R"))
    for pid, line in procs.items():
        (root / str(pid)).mkdir()
        (root / str(pid) / "stat").write_text(line)
    return root


@pytest.mark.covers("A-07", "A-01")
@pytest.mark.parametrize("line,start,want", [
    (_stat_line(100, "R"), 12345, True),
    (_stat_line(100, "S"), 12345, True),
    (_stat_line(100, "D"), 12345, True),
    (_stat_line(100, "T"), 12345, True),               # stopped is not dead
    (_stat_line(100, "S"), None, True),                # a pending record from before the start time was kept
    (_stat_line(100, "Z"), 12345, False),              # the leader is a zombie (its last thread may hold the lock)
    (_stat_line(100, "X"), 12345, False),
    (_stat_line(100, "R", flags=0x400104), 12345, False),  # PF_EXITING: the kernel is tearing it down
    (_stat_line(100, "S", start=99999), 12345, False),  # the pid now belongs to another process
    (None, 12345, False),                               # the pid is gone
])
def test_a07_runner_process_liveness_from_proc(tmp_path, line, start, want):
    root = _fake_proc(tmp_path, {100: line} if line else {})
    assert runner.runner_process_alive({"pid": 100, "start": start}, proc_root=root) is want


@pytest.mark.covers("A-07")
def test_a07_without_proc_the_lock_alone_decides(tmp_path):
    assert runner.runner_process_alive({"pid": 100, "start": 1}, proc_root=tmp_path / "noproc") is None
    root = _fake_proc(tmp_path, {100: _stat_line(100, "S", start=777)})
    assert runner.proc_start_time(100, proc_root=root) == 777
    assert runner.proc_start_time(101, proc_root=root) is None


@pytest.mark.covers("C-05")
def test_c05_worker_reads_its_cgroup_v2_path():
    assert worker.cgroup_v2_path(f"0::{CG}\n") == CG
    assert worker.cgroup_v2_path(f"12:memory:/x\n1:name=systemd:/y\n0::{CG}\n") == CG
    assert worker.cgroup_v2_path("12:memory:/x\n1:name=systemd:/y\n") is None  # cgroup v1 only: no evidence
