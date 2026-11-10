"""Slot geometry shared by tapes and pairs (PLAN.md §4.2; data_layer_and_tests.md §2A.2).

A slot is (start, end, strand), 0-based half-open on the contig; on a circular contig a
slot may wrap the origin (`end > contig_len`). Slots on one contig are in tape order
(the ordinal inside one tape version), which is coordinate order.
"""

from __future__ import annotations

# junction_topology codes (PLAN.md §4.2); 0 (intra-ORF fusion) needs segments and is not
# produced until they exist
CO_DIRECTIONAL, DIVERGENT, CONVERGENT, NESTED_ANTISENSE = 1, 2, 3, 4


def junction(left: tuple[int, int, int], right: tuple[int, int, int]) -> int:
    """Topology of two slots, `left` before `right` along the contig (coordinates already unwrapped)."""
    (ls, le, lstr), (rs, re_, rstr) = left, right
    if lstr == rstr:
        return CO_DIRECTIONAL
    if (ls <= rs and re_ <= le) or (rs <= ls and le <= re_):
        return NESTED_ANTISENSE
    return DIVERGENT if lstr == -1 else CONVERGENT


def cyclic_offset(i: int, j: int, n: int) -> int:
    """Signed slot offset from i to j on a circle of n slots: shortest way, +n/2 on a tie."""
    d = (j - i) % n
    return d if d <= n - d else d - n


def unwrap(slot: tuple[int, int, int], shift: int) -> tuple[int, int, int]:
    s, e, strand = slot
    return s + shift, e + shift, strand


def pair_geometry(slots: list[tuple[int, int, int]], i: int, j: int, offset: int,
                  contig_len: int, circular: bool) -> tuple[int, int, bool]:
    """(junction_topology, delta_bp, crosses_origin) of anchor i and partner j at `offset`.

    The pair is ordered left to right along the contig; on a circular contig, a pair whose
    shortest path crosses the origin shifts the right slot by one contig length. delta_bp is
    the gap between the two slots, or minus their overlap: max(starts) - min(ends), which
    reverse complementing leaves unchanged (MR-03) even when slots overlap or nest.
    """
    a, b = slots[i], slots[j]
    if offset > 0:
        left, right = a, b
        crosses = circular and j < i
    else:
        left, right = b, a
        crosses = circular and i < j
    if crosses:
        right = unwrap(right, contig_len)
    return junction(left, right), max(left[0], right[0]) - min(left[1], right[1]), crosses
