"""Ring R1: the pure parts of synthworld's fusions, splits and NOGD (synthworld_spec §2.1 S7, §3.3 truth_fusion;
card F-02 session 6). No ZOMBI2 here: the micro world is tested in tests/r2/test_synthworld_fusions.py.

The oracles are written in this file: a base-by-base walk that tags every base of an edited record with the
old position it came from, string slicing and reverse complements, and the genetic code spelled out.
"""

from __future__ import annotations

import copy

import numpy as np
import pytest

from biobox.synthworld import backbone, fusions, mobile

pytestmark = pytest.mark.covers("F-02")  # ring R1 (by directory)

CODE11_BY_AA = {
    "F": "TTT TTC", "L": "TTA TTG CTT CTC CTA CTG", "I": "ATT ATC ATA", "M": "ATG", "V": "GTT GTC GTA GTG",
    "S": "TCT TCC TCA TCG AGT AGC", "P": "CCT CCC CCA CCG", "T": "ACT ACC ACA ACG", "A": "GCT GCC GCA GCG",
    "Y": "TAT TAC", "*": "TAA TAG TGA", "H": "CAT CAC", "Q": "CAA CAG", "N": "AAT AAC", "K": "AAA AAG",
    "D": "GAT GAC", "E": "GAA GAG", "C": "TGT TGC", "W": "TGG", "R": "CGT CGC CGA CGG AGA AGG",
    "G": "GGT GGC GGA GGG",
}
CODE11 = {c: aa for aa, cs in CODE11_BY_AA.items() for c in cs.split()}
SENSE = sorted(c for c, a in CODE11.items() if a != "*")


def rc(s: str) -> str:
    return "".join({"A": "T", "C": "G", "G": "C", "T": "A"}[c] for c in reversed(s))


def translate(cds: str) -> str:
    aa = [CODE11[cds[i:i + 3]] for i in range(0, len(cds), 3)]
    assert cds[:3] in ("ATG", "GTG", "TTG") and aa[-1] == "*" and "*" not in aa[:-1], cds
    return "M" + "".join(aa[1:-1])


def cut(seq: str, s: int, e: int) -> str:
    return (seq + seq)[s:e] if e > len(seq) else seq[s:e]


def gene(rng, n_codons: int, start: str = "ATG", stop: str = "TAA") -> str:
    return start + "".join(SENSE[int(i)] for i in rng.integers(len(SENSE), size=n_codons)) + stop


# --- the lift ------------------------------------------------------------------------------------------------

def _tagged(N: int, edits):
    """The edited record as a list of old positions (None for new bases), walked base by base."""
    by_start = {}
    for a, b, l in edits:
        by_start.setdefault(a, []).append((b, l))
    out, i = [], 0
    while i <= N:
        reps = sorted(by_start.get(i, []), key=lambda x: x[0] - i)      # an insertion before a replacement
        jump = None
        for b, l in reps:
            out += [None] * l
            if b > i:
                jump = b
        if jump is not None:
            i = jump
            continue
        if i < N:
            out.append(i)
        i += 1
    return out


def test_edits_lift_matches_a_base_by_base_walk():
    rng = np.random.default_rng(7)
    n_checked = 0
    for trial in range(300):
        N = int(rng.integers(40, 400))
        cuts = sorted(set(int(x) for x in rng.integers(1, N, size=int(rng.integers(0, 10)))))
        edits, k = [], 0
        while k + 1 < len(cuts):
            a = cuts[k]
            b = a if rng.random() < 0.4 else cuts[k + 1]
            edits.append((a, b, int(rng.integers(0 if b > a else 1, 12))))     # an insertion adds bases
            k += 2 if b > a else 1
        if len({a for a, b, _ in edits if a == b}) < sum(a == b for a, b, _ in edits):
            continue
        ed = fusions.Edits(N, edits)
        tags = _tagged(N, edits)
        inv = {x: i for i, x in enumerate(tags) if x is not None}
        assert len(tags) == ed.new_length
        for _ in range(30):
            s = int(rng.integers(0, N))
            e = s + int(rng.integers(1, N))
            if e > N and rng.random() < 0.5:
                e = s + int(rng.integers(1, N - s + 1))
            want = list(range(s, e)) if e <= N else list(range(s, N)) + list(range(0, e - N))
            if ed.touched(s, e):
                with pytest.raises(ValueError):
                    ed.feature(s, e)
                # an old base inside a replaced span is gone, or an insertion splits the feature
                pos = [inv.get(x) for x in want]
                assert None in pos or any((p1 + 1) % len(tags) != p2 for p1, p2 in zip(pos, pos[1:]))
                continue
            s2, e2 = ed.feature(s, e)
            got = [tags[i % len(tags)] for i in range(s2, e2)]
            assert got == want, (edits, s, e)
            n_checked += 1
    assert n_checked > 2000


def test_apply_edits_writes_the_blocks_where_the_lift_says():
    seq = "".join("ACGT"[i % 4] for i in range(100))
    edits = [(10, 20, "NNN"), (50, 50, "XXXXX"), (70, 72, "")]
    new, starts = fusions.apply_edits(seq, edits)
    assert new == seq[:10] + "NNN" + seq[20:50] + "XXXXX" + seq[50:70] + seq[72:]
    assert [new[s:s + len(b)] for s, (_, _, b) in zip(starts, edits)] == ["NNN", "XXXXX", ""]
    ed = fusions.Edits(100, [(a, b, len(x)) for a, b, x in edits])
    assert ed.new_starts() == starts and ed.new_length == len(new)
    assert new[slice(*ed.feature(30, 40))] == seq[30:40]


def test_edits_refuse_overlaps_and_an_insertion_at_zero():
    with pytest.raises(ValueError):
        fusions.Edits(100, [(10, 20, 1), (15, 30, 1)])
    with pytest.raises(ValueError):
        fusions.Edits(100, [(10, 10, 1), (10, 10, 2)])
    with pytest.raises(ValueError):
        fusions.Edits(100, [(0, 0, 3)])
    with pytest.raises(ValueError):
        fusions.Edits(100, [(90, 110, 3)])
    with pytest.raises(ValueError):
        fusions.Edits(100, [(40, 40, 0)])


# --- fusions and splits ----------------------------------------------------------------------------------------

def _pair(rng, gap: int):
    """Two co-directional genes a, b on +: a tight gap (>= 0) or an overlap (-1, -4) made by hand."""
    a = gene(rng, 30, stop="TGA" if gap == -4 else "TAA")
    b = gene(rng, 25)
    if gap == -1:                       # TAATG: the stop's last base is the start's first
        b = "A" + b[1:]
        joined = a + b[1:]
    elif gap == -4:                     # ATGA: TGA overlaps ATG, the A before T of TGA is the start's first
        a = a[:-4] + "A" + "TGA"        # last sense codon ends in A
        while CODE11[a[-6:-3]] == "*":
            a = gene(rng, 30, stop="TGA")
            a = a[:-4] + "A" + "TGA"
        b = "ATG" + "A" + b[4:]         # the second codon starts with A (the stop's last base)
        joined = a[:-4] + b
        assert joined[len(a) - 4:len(a)] == "ATGA"
    else:
        joined = a + "C" * gap + b
    return a, b, joined


@pytest.mark.parametrize("gap", [0, 7, -1, -4])
@pytest.mark.parametrize("strand", [1, -1])
def test_fusion_edit_joins_two_genes_into_one_orf(gap, strand):
    rng = np.random.default_rng(100 + gap)
    a, b, joined = _pair(rng, gap)
    pa, pb = translate(a), translate(b)
    linker = "GGTTCTGCA"                # G S A
    flank_l, flank_r = "TTTTTTTTTT", "AAAAAAAAAA"
    s0 = len(flank_l)
    if strand == 1:
        rec = flank_l + joined + flank_r
        ga = (s0, s0 + len(a), 1)
        gb = (s0 + len(joined) - len(b), s0 + len(joined), 1)
    else:
        rec = flank_l + rc(joined) + flank_r
        ga = (s0 + len(joined) - len(a), s0 + len(joined), -1)
        gb = (s0, s0 + len(b), -1)
    x, y, blk = fusions.fusion_edit(ga, gb, linker)
    new, (b0,) = fusions.apply_edits(rec, [(x, y, blk)])
    fused = new[s0:len(new) - len(flank_r)]
    cds = fused if strand == 1 else rc(fused)
    assert translate(cds) == pa + "GSA" + pb[1:]
    assert new[:s0] == flank_l and new.endswith(flank_r)


@pytest.mark.parametrize("strand", [1, -1])
def test_split_edit_makes_two_orfs_with_the_protein_halves(strand):
    rng = np.random.default_rng(3)
    g = gene(rng, 60, start="GTG")
    prot = translate(g)
    block = "TAA" + "CCCCAGGAGGCCCCCCC" + "ATG"
    rec = "T" * 12 + (g if strand == 1 else rc(g)) + "A" * 9
    s, e = 12, 12 + len(g)
    for p in (1, 17, 59):
        x, y, blk = fusions.split_edit((s, e, strand), p, block)
        new, _ = fusions.apply_edits(rec, [(x, y, blk)])
        body = new[s:e + len(block)]
        body = body if strand == 1 else rc(body)
        h1, h2 = body[:3 * p + 3], body[3 * p + len(block) - 3:]
        assert translate(h1) == prot[:p] and translate(h2) == "M" + prot[p:]
    with pytest.raises(ValueError):
        fusions.split_edit((s, e, strand), 0, block)
    with pytest.raises(ValueError):
        fusions.split_edit((s, e, strand), 61, block)


def test_split_point_lies_between_the_two_domains_and_arch_round_trips():
    arch = fusions.parse_arch("PF00001.1:5-80;SYN0002:91-200;PF:3:1:210-300")
    assert arch == [("PF00001.1", 5, 80), ("SYN0002", 91, 200), ("PF:3:1", 210, 300)]
    assert fusions.format_arch(arch) == "PF00001.1:5-80;SYN0002:91-200;PF:3:1:210-300"
    assert fusions.split_point(arch, 1) == 85 and fusions.split_point(arch, 2) == 205
    assert fusions.split_point([("a", 1, 50), ("b", 50, 90)], 1) == 50
    assert fusions.parse_arch("") == []


# --- presets ---------------------------------------------------------------------------------------------------

def test_micro_preset_has_the_section_5_1_counts_and_check_config_refuses_bad_ones():
    p = backbone.load_preset("micro")
    cfg = p["fusions"]
    assert (cfg["fusions"], cfg["splits"], cfg["nogd"]) == (15, 5, 2)
    fusions.check_config(cfg)
    for key, value in (("clade_size", [0, 4]), ("clade_size", [5, 2]), ("split_spacer", [10, 40]),
                       ("linker_residues", "GSX"), ("nogd", 120), ("fusions", -1)):
        q = copy.deepcopy(cfg)
        q[key] = value
        with pytest.raises(ValueError):
            fusions.check_config(q)


def test_phage_world_preset_leaves_the_later_sections_out():
    """Adding a later stage's section to the preset must not change the phage world (and so S6's digest)."""
    p = backbone.load_preset("micro")
    q = mobile.phage_preset(p)
    assert "fusions" not in q
    before = {k: v for k, v in p.items() if k not in ("phages", "prophages", "crispr", "fusions")}
    assert set(q) == set(before)
