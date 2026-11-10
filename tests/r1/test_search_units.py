"""Ring R1: search units against oracles written from the definitions (S-04, S-09, S-12 lint, S-13).

ORFs: an independent codon table and a per-position definition (an ORF starts after a stop,
or at a linear contig end, and runs to the next stop). Neighborhoods: windows recomputed from
coordinates by enumerating both ways round. Recipes: the paper checks that refuse a recipe.
"""

from __future__ import annotations

import dataclasses
import json
from collections import Counter
from pathlib import Path

import pytest
from hypothesis import assume, given, settings, strategies as st

from biobox.search import export, nbr, orfs, recipe
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
    short = BASE.replace("{positive: [{sequence: MKVLAAGIVGLLL}]}",
                         "{positive: [{%s: %s, label: T7 gp1}]}" % (kind, value))
    long = BASE.replace("{positive: [{sequence: MKVLAAGIVGLLL}]}",
                        "{positive: [{kind: %s, value: %s, label: T7 gp1}]}" % (kind, value))
    a, _ = recipe.load(short)
    b, _ = recipe.load(long)
    assert a == b and a["controls"]["positive"][0]["label"] == "T7 gp1"
    assert recipe.recipe_hash(a) == recipe.recipe_hash(b)
    s, _ = recipe.load(BASE.replace("seed: [{sequence: MKVLAAGIVGLLL}]", "seed: [{%s: %s, label: x}]" % (kind, value)))
    assert s["seed"][0] == {"kind": kind, "value": a["controls"]["positive"][0]["value"], "label": "x"}


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
    occ = {"contig_uid": cu, "slot_ordinal": 5, "feature_uid": fu, "genome_uid": cu, "start": 10, "end": 20,
           "strand": 1, "tier": 1}
    [rec] = export.hit_records([hit], {pu: [occ]}, tape_version=tv)
    assert rec["tape_version"] == tv
    assert [(o["slot_ordinal"], o["feature_uid"]) for o in rec["occurrences"]] == [(5, fu.hex())]
    # every record kind whose keys (or nested occurrences) name an ordinal names the tape version too
    for kind, spec in export.SCHEMAS.items():
        keys = set(spec["keys"])
        if keys & ORDINAL_KEYS:
            assert {"feature_uid", "tape_version"} <= keys, kind
            assert "anchor_ordinal" not in keys or "anchor_feature_uid" in keys, kind
        if "occurrences" in keys:
            assert "tape_version" in keys, kind
