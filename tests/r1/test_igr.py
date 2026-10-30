"""Ring R1: IGR construction against a per-base coverage oracle (R-05 logic; Hypothesis)."""

from __future__ import annotations

from hypothesis import given, settings
from hypothesis import strategies as st

from biobox.l1.igr import igrs_of_contig


def oracle(slots, length, circular):
    """Maximal uncovered runs from a per-base mask; on a circular contig the run across the origin merges."""
    covered = [False] * length
    for _, s, e in slots:
        for p in range(s, e):
            covered[p % length] = True
    runs, i = [], 0
    while i < length:
        if covered[i]:
            i += 1
            continue
        j = i
        while j < length and not covered[j]:
            j += 1
        runs.append([i, j])
        i = j
    if circular and len(runs) > 1 and runs[0][0] == 0 and runs[-1][1] == length:
        first = runs.pop(0)
        runs[-1][1] = length + first[1]
    return sorted(tuple(r) for r in runs)


@st.composite
def contigs(draw):
    length = draw(st.integers(1, 400))
    circular = draw(st.booleans())
    n = draw(st.integers(0, 8))
    slots = []
    for k in range(n):
        s = draw(st.integers(0, length - 1))
        max_len = length if circular else length - s
        g = draw(st.integers(1, max(1, min(max_len, 120))))
        e = s + g
        if not circular:
            e = min(e, length)
        slots.append((k, s, e))
    if sum(1 for _, s, e in slots if e > length) > 1:
        slots = [x for x in slots if x[2] <= length] + [next(x for x in slots if x[2] > length)]
    return slots, length, circular


@settings(max_examples=600, deadline=None)
@given(contigs())
def test_igrs_equal_the_uncovered_runs(case):
    slots, length, circular = case
    got = igrs_of_contig(slots, length, circular)
    assert sorted((s, e) for s, e, _, _ in got) == oracle(slots, length, circular)
    ends = {o: e for o, _, e in slots}
    starts = {o: s for o, s, _ in slots}
    for s, e, left, right in got:
        assert 0 <= s < length and s < e and (e <= length or (circular and e <= s + length))
        if left is not None:
            assert ends[left] % length == s % length or ends[left] == s
        if right is not None:
            assert starts[right] == e % length or starts[right] == e


def test_examples():
    # linear: edges have null flanks; overlapping genes leave no IGR
    assert igrs_of_contig([(0, 10, 20), (1, 15, 30), (2, 40, 50)], 60, False) == [
        (0, 10, None, 0), (30, 40, 1, 2), (50, 60, 2, None)]
    # circular, no wrap: the run across the origin has end > length
    assert igrs_of_contig([(0, 10, 20), (1, 40, 50)], 60, True) == [(20, 40, 0, 1), (50, 70, 1, 0)]
    # circular with a wrapping gene covering [55, 65) -> [0, 5) is covered too
    assert igrs_of_contig([(0, 10, 20), (1, 55, 65)], 60, True) == [(5, 10, 1, 0), (20, 55, 0, 1)]
    # no genes
    assert igrs_of_contig([], 60, True) == [(0, 60, None, None)]
