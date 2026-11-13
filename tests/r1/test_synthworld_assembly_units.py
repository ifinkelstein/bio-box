"""Ring R1: the pure parts of synthworld's call errors, genetic codes and assembly (synthworld_spec §2.1 S8-S9,
§3.3 truth_callerr and truth_edges; card F-02 session 7). No ZOMBI2 here: the micro world is tested in
tests/r2/test_synthworld_assembly.py.

The oracles are written in this file: the genetic codes spelled out, a base-by-base walk that tags every base
of an edited record with the old position it came from, contigs cut from a tagged circle by slicing, and a
brute-force ORF scan.
"""

from __future__ import annotations

import numpy as np
import pytest

from biobox.synthworld import assembly

pytestmark = pytest.mark.covers("F-02")  # ring R1 (by directory)

CODE11_BY_AA = {
    "F": "TTT TTC", "L": "TTA TTG CTT CTC CTA CTG", "I": "ATT ATC ATA", "M": "ATG", "V": "GTT GTC GTA GTG",
    "S": "TCT TCC TCA TCG AGT AGC", "P": "CCT CCC CCA CCG", "T": "ACT ACC ACA ACG", "A": "GCT GCC GCA GCG",
    "Y": "TAT TAC", "*": "TAA TAG TGA", "H": "CAT CAC", "Q": "CAA CAG", "N": "AAT AAC", "K": "AAA AAG",
    "D": "GAT GAC", "E": "GAA GAG", "C": "TGT TGC", "W": "TGG", "R": "CGT CGC CGA CGG AGA AGG",
    "G": "GGT GGC GGA GGG",
}
CODE11 = {c: aa for aa, cs in CODE11_BY_AA.items() for c in cs.split()}
CODES = {11: CODE11, 4: dict(CODE11, TGA="W"), 25: dict(CODE11, TGA="G")}   # NCBI tables 4 and 25: TGA only
SENSE = sorted(c for c, a in CODE11.items() if a != "*")
STARTS = ("ATG", "GTG", "TTG")


def rc(s: str) -> str:
    return "".join({"A": "T", "C": "G", "G": "C", "T": "A"}[c] for c in reversed(s))


def protein(cds: str, code: int = 11) -> str:
    t = CODES[code]
    aa = [t[cds[i:i + 3]] for i in range(0, len(cds), 3)]
    assert len(cds) % 3 == 0 and cds[:3] in STARTS and aa[-1] == "*" and "*" not in aa[:-1], cds
    return "M" + "".join(aa[1:-1])


def gene(rng, n_codons: int, start: str = "ATG", stop: str = "TAA", pool=SENSE) -> str:
    return start + "".join(pool[int(i)] for i in rng.integers(len(pool), size=n_codons)) + stop


# --- genetic codes ---------------------------------------------------------------------------------------------

def test_codes_4_and_25_differ_from_11_only_at_tga():
    for code in (11, 4, 25):
        t = assembly.code_table(code)
        assert t == CODES[code]
    rng = np.random.default_rng(3)
    n_tga = 0
    for _ in range(50):
        g = gene(rng, 40, pool=SENSE + ["TGA"] * 2, stop="TAA")
        has = "TGA" in [g[i:i + 3] for i in range(3, len(g) - 3, 3)]
        n_tga += has
        assert (assembly.orf_protein(g, 11) is None) == has
        if not has:
            assert assembly.orf_protein(g, 11) == protein(g, 11)
        for code in (4, 25):
            assert assembly.orf_protein(g, code) == protein(g, code)
    assert 0 < n_tga < 50
    with pytest.raises(ValueError):
        assembly.code_table(2)


# --- the lift ----------------------------------------------------------------------------------------------------

def _tagged(N: int, edits):
    """The edited record as a list of old positions (None for new bases), walked base by base."""
    by_start = {}
    for a, b, l in edits:
        by_start.setdefault(a, []).append((b, l))
    out, i = [], 0
    while i <= N:
        reps = sorted(by_start.get(i, []), key=lambda x: x[0] - i)
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


def _random_edits(rng, N):
    pts = sorted(set(int(x) for x in rng.integers(1, N, size=int(rng.integers(0, 40)))))
    edits, k = [], 0
    while k < len(pts):
        a = pts[k]
        kind = int(rng.integers(3))
        if kind == 0:                                    # insertion
            edits.append((a, a, int(rng.integers(1, 5))))
            k += 1
        elif k + 1 < len(pts):                           # replacement [a, next point)
            edits.append((a, pts[k + 1], int(rng.integers(0, 5))))
            k += 2
        else:
            k += 1
    return edits


def test_lift_touch_test_and_start_match_a_base_by_base_walk():
    rng = np.random.default_rng(11)
    n = 0
    for _ in range(300):
        N = int(rng.integers(50, 400))
        edits = _random_edits(rng, N)
        lift = assembly.Lift(N, edits)
        tag = _tagged(N, edits)
        assert len(tag) == lift.new_length
        where = {p: i for i, p in enumerate(tag) if p is not None}
        for _ in range(30):
            s = int(rng.integers(0, N))
            e = s + int(rng.integers(1, N))
            bases = [x % N for x in range(s, e)]
            brute = any((a < b and any(a <= x < b for x in bases))
                        or (a == b and any(x == a for x in bases[1:])) for a, b, _ in edits)
            if e > N:     # an insertion at N sits between the last base and the first
                brute = brute or any(a == b == N for a, b, _ in edits)
            assert lift.touched(s, e) == brute, (edits, s, e)
            if not brute:
                s2, e2 = lift.feature(s, e)
                assert [tag[x % lift.new_length] for x in range(s2, e2)] == bases
                n += 1
            if s in where:
                assert lift.start(s) == where[s]
    assert n > 500


# --- pieces, placement and contig coordinates ----------------------------------------------------------------

def _contigs_by_slicing(N, ps, orient):
    """Each contig as the list of frame positions of its bases (contig order) and the strand they read on."""
    out = []
    for p, o in zip(ps, orient):
        pos = [x % N for x in range(p.a, p.b)]
        out.append(pos if o == 1 else pos[::-1])
    return out


def test_pieces_tile_the_replicon_and_place_matches_sliced_contigs():
    rng = np.random.default_rng(5)
    n_parts = 0
    for trial in range(400):
        N = int(rng.integers(30, 300))
        circular = bool(rng.random() < 0.8)
        m = int(rng.integers(0, 6))
        cuts = {int(c): ("break_random", "") for c in rng.integers(0 if circular else 1, N, size=m)}
        rot = int(rng.integers(N))
        ps = assembly.pieces(N, circular, cuts, rot)
        # tiling: every base exactly once
        seen = sorted(x % N for p in ps for x in range(p.a, p.b))
        assert seen == list(range(N))
        assert all(p.a < p.b <= p.a + N for p in ps)
        if circular and not cuts:
            assert len(ps) == 1 and ps[0].circular and ps[0].a == rot
        orient = [1 if rng.random() < 0.5 else -1 for _ in ps]
        tags = _contigs_by_slicing(N, ps, orient)
        for _ in range(20):
            s = int(rng.integers(0, N))
            length = int(rng.integers(1, N if circular else N - s + 1))
            e = s + length
            if not circular and e > N:
                continue
            strand = 1 if rng.random() < 0.5 else -1
            fpos = [x % N for x in range(s, e)]
            want = []
            for k, tg in enumerate(tags):
                L = len(tg)
                circ_contig = ps[k].circular
                inside = [i for i, x in enumerate(tg) if x in set(fpos)]
                if not inside:
                    continue
                if circ_contig:
                    # one run, possibly across the contig's origin: start at the base reading the feature's first
                    # (o = 1) or last (o = -1) frame position
                    first = tg.index(fpos[0] if orient[k] == 1 else fpos[-1])
                    want.append((k, first, first + length, strand * orient[k], False, False))
                    continue
                runs, cur = [], [inside[0]]
                for i in inside[1:]:
                    if i == cur[-1] + 1:
                        cur.append(i)
                    else:
                        runs.append(cur)
                        cur = [i]
                runs.append(cur)
                for run in runs:
                    lo, hi = run[0], run[-1] + 1
                    # the feature's 5' end is its frame start on +, its frame end on -
                    fr_lo, fr_hi = tg[lo], tg[hi - 1]
                    if orient[k] == -1:
                        fr_lo, fr_hi = fr_hi, fr_lo
                    clip_left = fr_lo != fpos[0]
                    clip_right = fr_hi != fpos[-1]
                    p5, p3 = (clip_left, clip_right) if strand == 1 else (clip_right, clip_left)
                    want.append((k, lo, hi, strand * orient[k], p5, p3))
            got = []
            for j, x, y, cl, cr in assembly.place(ps, N, circular, s, e):
                L = ps[j].b - ps[j].a
                cs, ce, cst = assembly.to_contig(L, orient[j], x, y, strand)
                p5, p3 = (cl, cr) if strand == 1 else (cr, cl)
                got.append((j, cs, ce, cst, p5, p3))
            assert sorted(got) == sorted(want), (N, circular, cuts, rot, s, e, orient)
            n_parts += len(got)
    assert n_parts > 5000


def test_linear_replicon_ends_and_circular_wrap_piece():
    ps = assembly.pieces(100, False, {30: ("break_repeat", "F1"), 70: ("break_random", "")})
    assert [(p.a, p.b, p.left, p.right, p.fam_left, p.fam_right) for p in ps] == [
        (0, 30, "replicon_end", "break_repeat", "", "F1"), (30, 70, "break_repeat", "break_random", "F1", ""),
        (70, 100, "break_random", "replicon_end", "", "")]
    ps = assembly.pieces(100, True, {30: ("break_random", ""), 70: ("break_repeat", "F2")})
    assert [(p.a, p.b, p.left, p.right) for p in ps] == [(30, 70, "break_random", "break_repeat"),
                                                       (70, 130, "break_repeat", "break_random")]


# --- expected calls ------------------------------------------------------------------------------------------

def _orfs_brute(cds):
    """All maximal stop-to-stop stretches per frame of the gene's own strand: (frame start, stop end)."""
    out = []
    for f in range(3):
        prev = f
        for i in range(f, len(cds) - 2, 3):
            if CODE11[cds[i:i + 3]] == "*":
                out.append((prev, i + 3))
                prev = i + 3
    return out


def test_code11_calls_are_the_5_prime_and_3_prime_orf_pieces():
    rng = np.random.default_rng(17)
    n = 0
    for _ in range(400):
        g = gene(rng, int(rng.integers(60, 200)))
        L = len(g)
        kind = int(rng.integers(3))
        q = int(rng.integers(L // 3, 2 * L // 3))
        if kind == 0:
            g = g[:q] + g[q + 1:]                            # deletion
        elif kind == 1:
            g = g[:q] + "ACGT"[int(rng.integers(4))] + g[q:]  # insertion
        else:
            c = q // 3
            g = g[:3 * c] + "TAA" + g[3 * c + 3:]            # internal stop
        L = len(g)
        orfs = _orfs_brute(g)
        want = []
        # 5' piece: frame 0 from the start to its first stop
        five = next(((a, b) for a, b in orfs if a == 0), None)
        if five and (five[1] - 3) // 3 >= 30:
            want.append(five)
        # 3' piece: the stretch ending at the true stop, from its first start codon
        three = next(((a, b) for a, b in orfs if b == L), None)
        if three:
            st = next((i for i in range(three[0], L - 3, 3) if g[i:i + 3] in STARTS), None)
            if st is not None and (L - 3 - st) // 3 >= 30 and (st, L) not in want:
                want.append((st, L))
        assert assembly.code11_calls(g, 30) == sorted(want), g
        n += len(want)
    assert n > 300


def test_n50_matches_a_brute_force_definition():
    rng = np.random.default_rng(2)
    for _ in range(200):
        xs = [int(x) for x in rng.integers(1, 1000, size=int(rng.integers(1, 30)))]
        tot = sum(xs)
        want = max(x for x in xs if sum(y for y in xs if y >= x) * 2 >= tot)
        assert assembly.n50(xs) == want


# --- recoding for codes 4 and 25 --------------------------------------------------------------------------------

def _rep(seq, genes):
    N = len(seq)
    gs = [{"label": f"x{i}", "s": s, "e": e, "st": st, "rep": 0} for i, (s, e, st) in enumerate(genes)]
    return assembly._Rep(0, "r", seq, N, "chromosome", "circular", assembly._counts(N, gs), np.zeros(N, bool)), gs


def test_recoding_keeps_every_protein_under_codes_4_and_25_and_removes_tga_stops():
    rng = np.random.default_rng(23)
    n_overlap = 0
    for trial in range(200):
        code = 4 if trial % 2 else 25
        pool = SENSE + ["TGG"] * 6 + ["GGA"] * 6
        a = gene(rng, int(rng.integers(20, 60)), stop="TGA", pool=pool)
        b = gene(rng, int(rng.integers(20, 60)), start="ATG", stop=["TAA", "TGA"][int(rng.integers(2))], pool=pool)
        b = b[:3] + "AAA" + b[6:]            # b's second codon starts with A: a's stop TGA is b's T G A
        c = gene(rng, int(rng.integers(20, 60)), stop="TGA", pool=pool)
        # a -4 overlap a|b: a's last sense codon must end in A, then TGA; b starts ATG at that A
        a = a[:-6] + "GCA" + "TGA"
        ab = a[:-4] + b                      # ...GC A TG A..: b's ATG starts at a's last sense base
        assert ab[len(a) - 4:len(a) - 1] == "ATG" and ab[len(a) - 3:len(a)] == "TGA"
        left = "".join("ACGT"[int(i)] for i in rng.integers(4, size=20))
        mid = "".join("ACGT"[int(i)] for i in rng.integers(4, size=30))
        strand = 1 if trial % 3 else -1
        if strand == 1:
            seq = left + ab + mid + c + left
            sa = len(left)
            genes = [(sa, sa + len(a), 1), (sa + len(a) - 4, sa + len(a) - 4 + len(b), 1),
                     (sa + len(ab) + len(mid), sa + len(ab) + len(mid) + len(c), 1)]
        else:
            body = left + ab + mid + c + left
            seq = rc(body)
            M = len(body)
            sa = len(left)
            fw = [(sa, sa + len(a)), (sa + len(a) - 4, sa + len(a) - 4 + len(b)),
                  (sa + len(ab) + len(mid), sa + len(ab) + len(mid) + len(c))]
            genes = [(M - e, M - s, -1) for s, e in fw]
        rp, gs = _rep(seq, genes)
        mods = {}
        assembly._recode(rp, sorted(gs, key=lambda g: g["s"]), code, np.random.default_rng(trial),
                         {"code4_tgg_to_tga": 0.9, "code25_gly_to_tga": 0.5}, mods)
        edits = sorted(rp.edits, key=lambda x: (x[0], x[1]))
        n_overlap += sum(1 for x in edits if x[3] == "recode_overlap")
        new, _ = assembly.fusions.apply_edits(seq, [(x, y, blk) for x, y, blk, _, _ in edits])
        lift = assembly.Lift(len(seq), [(x, y, len(blk)) for x, y, blk, _, _ in edits])
        for g, (s, e, st) in zip(gs, genes):
            old = seq[s:e] if st == 1 else rc(seq[s:e])
            s2 = lift.start(s) + mods.get(g["label"], {}).get("start_delta", 0)
            cds = new[s2:s2 + (e - s)]
            cds = cds if st == 1 else rc(cds)
            assert protein(cds, code) == protein(old, 11)
            assert cds[-3:] != "TGA"
    assert n_overlap == 200
