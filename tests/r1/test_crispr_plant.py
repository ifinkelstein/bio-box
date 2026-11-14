"""Ring R1: contigs with planted CRISPR arrays whose repeats and spacers are known by construction (card E-12).

The constructor is the second oracle of the CRISPR producer (`tests/r2/test_crispr_producer.py`, AGENTS §5.2): the
truth is written down as the sequence is made, never read back from MinCED. A contig is seeded random background
(no 23-base repeat outside the arrays, checked here) with arrays inserted at given positions. An array is `n` copies
of one real repeat with distinct random spacers between them; each spacer's first and last base rotate through the
four bases, so no base follows (or precedes) more than half the repeat copies and a repeat cannot be extended by
the 75% column rule CRT uses (its `getActualRepeatLength`). An array may be inserted reverse-complemented (the "array
on the reverse strand of the sequence"): then its elements, read on the + strand, are the reverse complements in
reverse order. A circular contig is made by rotating a linear one; an element's coordinates follow the producer's
convention (an element wholly past the origin at its position on the contig, one across it with end > n).

The tests here check the constructor itself: every element's bases are the contig's bases at its coordinates, the
background holds no long repeat, and rotation keeps every element's bases.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field

import pytest

pytestmark = pytest.mark.covers("E-12")

# real direct repeats: type I-F (Pectobacterium-like), I-E (E. coli K-12) and II-A (S. pyogenes)
REPEAT_IF = b"GTTCACTGCCGTACAGGCAGCTTAGAAA"  # 28
REPEAT_IE = b"CGGTTTATCCCCGCTGGCGCGGGGAACTC"  # 29
REPEAT_IIA = b"GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC"  # 36
_COMP = bytes.maketrans(b"ACGTN", b"TGCAN")


def revcomp(s: bytes) -> bytes:
    return s.translate(_COMP)[::-1]


@dataclass
class Element:
    kind: str  # repeat or spacer
    ordinal: int
    start: int
    end: int
    seq: bytes


@dataclass
class PlantedArray:
    start: int
    end: int
    elements: list[Element] = field(default_factory=list)


@dataclass
class Contig:
    seq: bytes
    arrays: list[PlantedArray]
    circular: bool = False


def _random(rng: random.Random, n: int) -> bytes:
    return bytes(rng.choice(b"ACGT") for _ in range(n))


def array_block(rng: random.Random, repeat: bytes, n_repeats: int, spacer_len: int = 33) -> tuple[bytes, list[Element]]:
    """The bases of one array (repeat, spacer, ..., repeat) and its elements, coordinates relative to the block."""
    parts, els, at = [], [], 0
    for i in range(n_repeats):
        parts.append(repeat)
        els.append(Element("repeat", 2 * i, at, at + len(repeat), repeat))
        at += len(repeat)
        if i + 1 < n_repeats:
            # first and last base rotate through ACGT, so no base sits beside more than half the copies
            inner = _random(rng, spacer_len + (i % 3) - 2)
            sp = b"ACGT"[i % 4:i % 4 + 1] + inner + b"TGCA"[i % 4:i % 4 + 1]
            parts.append(sp)
            els.append(Element("spacer", 2 * i + 1, at, at + len(sp), sp))
            at += len(sp)
    return b"".join(parts), els


def reversed_block(block: bytes, els: list[Element]) -> tuple[bytes, list[Element]]:
    """The block reverse-complemented, with its elements as they read on the + strand of the result."""
    n = len(block)
    out = [Element(e.kind, 0, n - e.end, n - e.start, revcomp(e.seq)) for e in reversed(els)]
    for k, e in enumerate(out):
        e.ordinal = k
    return revcomp(block), out


def plant(seed: int, length: int, arrays: list[tuple[int, bytes, int, bool]]) -> Contig:
    """A linear contig of `length` bases of seeded random background with arrays (position, repeat, repeats,
    reverse) inserted at increasing positions: each array's block replaces nothing, it is inserted, so the contig
    is longer than `length` by the arrays' bases. The background base beside an array adds one to a column of at
    least four, which stays under CRT's 75% (checked below)."""
    rng = random.Random(seed)
    bg = bytearray(_random(rng, length))
    seq, out, last = b"", [], 0
    for pos, repeat, n, reverse in sorted(arrays, key=lambda a: a[0]):
        block, els = array_block(rng, repeat, n)
        if reverse:
            block, els = reversed_block(block, els)
        seq += bytes(bg[last:pos])
        base = len(seq)
        seq += block
        out.append(PlantedArray(base, base + len(block),
                                [Element(e.kind, e.ordinal, base + e.start, base + e.end, e.seq) for e in els]))
        last = pos
    seq += bytes(bg[last:])
    return Contig(seq, out)


def rotate(c: Contig, k: int) -> Contig:
    """The circular contig `c` read from base k: an element (or array) at [s, e) moves to (s - k) mod n, keeping its
    length, so one across the new origin has end > n."""
    n = len(c.seq)

    def move(s: int, e: int) -> tuple[int, int]:
        s2 = (s - k) % n
        return s2, s2 + (e - s)

    arrays = []
    for a in c.arrays:
        s, e = move(a.start, a.end)
        arrays.append(PlantedArray(s, e, [Element(x.kind, x.ordinal, *move(x.start, x.end), x.seq) for x in a.elements]))
    return Contig(c.seq[k:] + c.seq[:k], arrays, circular=True)


def bases(c: Contig, s: int, e: int) -> bytes:
    """The contig's bases at [s, e), around the origin of a circular contig."""
    if e <= len(c.seq):
        return c.seq[s:e]
    assert c.circular
    return c.seq[s:] + c.seq[:e - len(c.seq)]


def contigs() -> dict[str, Contig]:
    """The planted contigs of the producer's tests (oracle 2): name -> contig."""
    plain = plant(11, 6000, [(2500, REPEAT_IF, 6, False)])
    reverse = plant(12, 6000, [(3100, REPEAT_IE, 5, True)])
    two = plant(13, 9000, [(1800, REPEAT_IF, 4, False), (6200, REPEAT_IIA, 7, False)])
    linear_twin = plant(14, 7000, [(3000, REPEAT_IE, 8, False)])  # the circular one below, before rotation
    a = linear_twin.arrays[0]
    k = (a.elements[5].start + a.elements[5].end) // 2  # inside the third spacer
    across = rotate(linear_twin, k)
    none = plant(15, 5000, [])
    none.circular = True
    return {"plain": plain, "reverse": reverse, "two": two, "linear_twin": linear_twin, "across": across, "none": none}


# --- the constructor's own checks ---------------------------------------------------------------------------------


def _kmers(s: bytes, k: int) -> dict[bytes, int]:
    out: dict[bytes, int] = {}
    for i in range(len(s) - k + 1):
        out[s[i:i + k]] = out.get(s[i:i + k], 0) + 1
    return out


def test_every_element_is_the_contigs_bases_at_its_coordinates():
    for name, c in contigs().items():
        n = len(c.seq)
        for a in c.arrays:  # around the origin of a circular contig, positions agree modulo n
            assert a.elements[0].start == a.start and (a.elements[-1].end - a.end) % n == 0, name
            assert [e.ordinal for e in a.elements] == list(range(len(a.elements))), name
            assert [e.kind for e in a.elements] == ["repeat", "spacer"] * (len(a.elements) // 2) + ["repeat"], name
            for x, y in zip(a.elements, a.elements[1:]):
                assert (x.end - y.start) % n == 0, name
            for e in a.elements:
                assert bases(c, e.start, e.end) == e.seq, (name, e)
                assert 0 <= e.start < len(c.seq) and e.end <= e.start + len(c.seq)


def test_repeats_and_spacers_are_within_minceds_length_ranges_and_spacers_differ():
    for name, c in contigs().items():
        for a in c.arrays:
            reps = {e.seq for e in a.elements if e.kind == "repeat"}
            sps = [e.seq for e in a.elements if e.kind == "spacer"]
            assert len(reps) == 1 and 23 <= len(next(iter(reps))) <= 47, name
            assert all(26 <= len(s) <= 50 for s in sps) and len(set(sps)) == len(sps), name
            # no base follows or precedes more than half of the repeats (no column CRT could extend)
            after = [bases(c, e.end, e.end + 1) for e in a.elements if e.kind == "repeat"]
            before = [bases(c, e.start - 1, e.start) if e.start else c.seq[-1:] for e in a.elements
                      if e.kind == "repeat"]
            for col in (after, before):
                assert max(col.count(b) for b in set(col)) <= len(col) / 2 + 1, (name, col)


def test_the_background_has_no_23_base_repeat_outside_the_arrays():
    for name, c in contigs().items():
        masked = bytearray(c.seq)
        for a in c.arrays:
            for i in range(a.start, a.end):
                masked[i % len(c.seq)] = ord("N")
        ring = bytes(masked) + (bytes(masked[:22]) if c.circular else b"")  # each 23-mer once, also across the origin
        rep = [k for k, v in _kmers(ring, 23).items() if v > 1 and b"N" not in k]
        assert rep == [], name


def test_rotation_moves_every_element_and_one_array_crosses_the_origin():
    cs = contigs()
    lin, circ = cs["linear_twin"], cs["across"]
    n = len(lin.seq)
    assert len(circ.seq) == n and circ.circular
    (a,) = circ.arrays
    assert a.start < n < a.end  # the planted array crosses the origin
    assert [e.seq for e in a.elements] == [e.seq for e in lin.arrays[0].elements]
    assert any(e.start < n < e.end for e in a.elements)  # the third spacer is cut by the origin
    assert any(e.start < a.start for e in a.elements)  # elements past the origin lie at their place on the contig


def test_the_reverse_array_reads_as_reverse_complements_in_reverse_order():
    rng = random.Random(3)
    block, els = array_block(rng, REPEAT_IE, 4)
    rblock, rels = reversed_block(block, els)
    assert [e.seq for e in rels] == [revcomp(e.seq) for e in reversed(els)]
    assert all(rblock[e.start:e.end] == e.seq for e in rels)
    assert {e.seq for e in rels if e.kind == "repeat"} == {revcomp(REPEAT_IE)}
