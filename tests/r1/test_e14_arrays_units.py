"""Ring R1, card E-14: the window in bases of a locus and an array's place in it (`biobox.search.arrays`), on
hand-made rows whose answers are worked out here by hand (no store).

The definition (module docstring of `biobox.search.arrays`): a gene window spans from the least start to the
greatest end of the locus's rows, a row reached across a circular origin moved by one contig length to the anchor's
side; a window in bases is [anchor start - W, anchor end + W). An array is within it when it shares at least one base
with the span, partial when part of it lies outside; on a circular contig it is moved by -n, 0 or +n, whichever puts
most of it inside (ties: the smaller gap to the anchor, then no move).
"""

from __future__ import annotations

import pytest

from biobox.search import arrays, nbr

pytestmark = pytest.mark.covers("E-14")  # ring R1 by its directory (tests/conftest.py)

CU, AF = b"\x01" * 16, b"\x02" * 16


def row(offset: int, start: int, end: int, strand: int = 1, flags: int = 0) -> dict:
    return {"contig_uid": CU, "anchor_feature_uid": AF, "offset": offset, "start": start, "end": end,
            "strand": strand, "flags": flags, "tape_version": "t-000000000000"}


def arr(start: int, end: int) -> dict:
    return {"feature_uid": b"\x03" * 16, "producer": "minced-0.4.2-w8-nr3-r23x47-s26x50", "start": start, "end": end,
            "strand": 1, "n_repeats": 3, "n_spacers": 2, "repeat": "ACGT"}


G3 = nbr.Policy.parse(window="3g")
B500 = nbr.Policy.parse(window="500bp")


def test_a_gene_window_spans_its_rows_and_a_base_window_the_anchor_plus_and_minus_w():
    rows = [row(-1, 100, 400), row(0, 500, 900), row(1, 950, 1200)]
    assert arrays.window(rows, G3, 10_000) == (100, 1200)
    assert arrays.window(rows, B500, 10_000) == (0, 1400)


def test_rows_across_the_origin_are_moved_to_the_anchors_side():
    n = 10_000
    # the anchor near the end; the next gene (offset +1) after the origin, the gene at -1 before it
    rows = [row(-1, 9_000, 9_400), row(0, 9_500, 9_900), row(1, 100, 300, flags=nbr.FLAG_CROSSES_ORIGIN)]
    assert arrays.window(rows, G3, n) == (9_000, n + 300)
    # the anchor just after the origin, a gene at -1 before it
    rows = [row(-1, 9_700, 9_950, flags=nbr.FLAG_CROSSES_ORIGIN), row(0, 50, 400)]
    assert arrays.window(rows, G3, n) == (9_700 - n, 400)


@pytest.mark.parametrize("a, inside, partial", [
    ((1200, 1300), False, None),   # starts at the span's end: no shared base
    ((1199, 1300), True, True),    # one base inside
    ((150, 300), True, False),     # wholly inside
    ((50, 101), True, True),       # its last base is the span's first
    ((50, 100), False, None),      # ends where the span starts
])
def test_an_array_is_within_when_it_shares_a_base_and_partial_when_it_runs_out(a, inside, partial):
    rows = [row(-1, 100, 400), row(0, 500, 900), row(1, 950, 1200)]
    lo, hi = arrays.window(rows, G3, 10_000)
    got = arrays.place(arr(*a), rows[1], lo, hi, 10_000, False)
    assert (got is not None) == inside
    if inside:
        assert got["partial"] is partial and got["crosses_origin"] is False


def test_gap_and_side_follow_the_anchors_strand():
    plus, minus = row(0, 500, 900, 1), row(0, 500, 900, -1)
    for anchor, down in ((plus, 1), (minus, -1)):
        after = arrays.place(arr(1000, 1100), anchor, 0, 2000, 10_000, False)
        before = arrays.place(arr(200, 300), anchor, 0, 2000, 10_000, False)
        over = arrays.place(arr(850, 950), anchor, 0, 2000, 10_000, False)
        assert (after["delta_bp"], after["side_rel"]) == (100, down)
        assert (before["delta_bp"], before["side_rel"]) == (200, -down)
        assert (over["delta_bp"], over["side_rel"]) == (-50, 0)


def test_on_a_circle_the_array_is_moved_where_most_of_it_is_inside():
    n = 10_000
    anchor = row(0, 9_500, 9_900)
    # an array just after the origin, in a window that runs past it: moved by +n
    got = arrays.place(arr(100, 300), anchor, 9_000, n + 400, n, True)
    assert got["crosses_origin"] is True and got["delta_bp"] == n + 100 - 9_900 and got["side_rel"] == 1
    assert (got["start"], got["end"]) == (100, 300)  # stored coordinates, not moved ones
    # the same contig read as linear: the array is not in the window
    assert arrays.place(arr(100, 300), anchor, 9_000, n + 400, n, False) is None
    # an array across the origin (end > n) beside an anchor just after it: moved by -n, wholly inside
    got = arrays.place(arr(9_800, n + 50), row(0, 100, 400), 9_700 - n, 800, n, True)
    assert got["crosses_origin"] is True and got["partial"] is False and got["side_rel"] == -1
    # inside without a move: not across
    got = arrays.place(arr(9_920, 9_990), anchor, 9_000, n + 400, n, True)
    assert got["crosses_origin"] is False and got["delta_bp"] == 20


def test_a_release_without_l1_crispr_files_is_not_read():
    class Files:
        def entries(self, name, **want):
            assert name == "l1_crispr", name  # nothing but the check is looked at
            return []

    assert arrays.near([row(0, 1, 2)], G3, files=Files()) is None
