"""Ring R1: the pure parts of synthworld's contamination and duplicates (synthworld_spec §2.1 S10, §3.3 truth_dup and
truth_contam; card F-02 session 7). No ZOMBI2 here: the micro world is tested in tests/r2/test_synthworld_assembly.py.

The oracles are written in this file: the genetic codes spelled out and every single-base change of every codon
tried, and tree path lengths summed edge by edge.
"""

from __future__ import annotations

import numpy as np
import pytest

from biobox.synthworld import mobile, sources

pytestmark = pytest.mark.covers("F-02")  # ring R1 (by directory)

CODE11_BY_AA = {
    "F": "TTT TTC", "L": "TTA TTG CTT CTC CTA CTG", "I": "ATT ATC ATA", "M": "ATG", "V": "GTT GTC GTA GTG",
    "S": "TCT TCC TCA TCG AGT AGC", "P": "CCT CCC CCA CCG", "T": "ACT ACC ACA ACG", "A": "GCT GCC GCA GCG",
    "Y": "TAT TAC", "*": "TAA TAG TGA", "H": "CAT CAC", "Q": "CAA CAG", "N": "AAT AAC", "K": "AAA AAG",
    "D": "GAT GAC", "E": "GAA GAG", "C": "TGT TGC", "W": "TGG", "R": "CGT CGC CGA CGG AGA AGG",
    "G": "GGT GGC GGA GGG",
}
CODE11 = {c: aa for aa, cs in CODE11_BY_AA.items() for c in cs.split()}
CODES = {11: CODE11, 4: dict(CODE11, TGA="W"), 25: dict(CODE11, TGA="G")}
RC = {"A": "T", "C": "G", "G": "C", "T": "A"}


def rc(s):
    return "".join(RC[c] for c in reversed(s))


def test_synonymous_sites_are_exactly_the_protein_keeping_changes():
    rng = np.random.default_rng(31)
    n_sites = 0
    for trial in range(150):
        code = (11, 4, 25)[trial % 3]
        L = int(rng.integers(200, 600))
        seq = "".join("ACGT"[int(i)] for i in rng.integers(4, size=L))
        pieces = []
        for _ in range(int(rng.integers(1, 4))):
            n = int(rng.integers(5, 40))
            s = int(rng.integers(0, L))
            e = s + 3 * n
            if e > L + L // 2:
                continue
            pieces.append((s, e, 1 if rng.random() < 0.5 else -1))
        cut = [bool(rng.random() < 0.3) for _ in pieces]               # genes cut by a contig edge
        whole = [p for p, c in zip(pieces, cut) if not c]
        part = [p for p, c in zip(pieces, cut) if c]
        got, cover = sources.synonymous_sites(seq, whole, code, part)
        # brute force: per base the pieces over it (cut ones too); per whole piece, internal codons alone in one
        # piece, every change tried
        cnt = [0] * L
        for s, e, _ in pieces:
            for x in range(s, e):
                cnt[x % L] += 1
        assert cover.tolist() == cnt
        want = {}
        for s, e, st in whole:
            fwd = (seq + seq)[s:e]
            cds = fwd if st == 1 else rc(fwd)
            for c in range(1, len(cds) // 3 - 1):
                # base i of codon c on the gene strand sits at forward offset 3c + i (+) or len - 1 - 3c - i (-)
                offs = [3 * c + i if st == 1 else len(cds) - 1 - 3 * c - i for i in range(3)]
                if any(cnt[(s + o) % L] != 1 for o in offs):
                    continue
                cod = cds[3 * c:3 * c + 3]
                for i in range(3):
                    keep = []
                    for b in "ACGT":
                        alt = cod[:i] + b + cod[i + 1:]
                        if b != cod[i] and CODES[code][alt] == CODES[code][cod]:
                            keep.append(b if st == 1 else RC[b])
                    if keep:
                        want[(s + offs[i]) % L] = keep
        assert {k: sorted(v) for k, v in got.items()} == {k: sorted(v) for k, v in want.items()}
        n_sites += len(want)
    assert n_sites > 2000


def test_tree_distance_sums_the_branches_between_two_tips():
    # ((a:1,b:2)x:0.5,(c:1.5,(d:0.25,e:0.75)y:0.5)z:1)r:0.1, built as species branches
    class B:
        def __init__(self, lineage, parent, branch_length):
            self.lineage, self.parent, self.branch_length = lineage, parent, branch_length
    edges = {"r": ("", 0.1), "x": ("r", 0.5), "a": ("x", 1.0), "b": ("x", 2.0), "z": ("r", 1.0), "c": ("z", 1.5),
             "y": ("z", 0.5), "d": ("y", 0.25), "e": ("y", 0.75)}
    order = ["r", "x", "a", "b", "z", "c", "y", "d", "e"]
    tree = mobile.species_tree([B(k, *edges[k]) for k in order], {"a", "b", "c", "d", "e"})

    def up(t):
        path = []
        while t:
            path.append(t)
            t = edges[t][0]
        return path

    for x in "abcde":
        for y in "abcde":
            px, py = up(x), up(y)
            common = set(px) & set(py)
            want = sum(edges[n][1] for n in px if n not in common) + sum(edges[n][1] for n in py if n not in common)
            assert sources.tree_distance(tree, x, y) == pytest.approx(want)
