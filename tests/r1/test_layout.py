"""Ring R1: bucket width and part cutting for prefix-bucketed tables (biobox.layout, LY-04, LY-05)."""

from __future__ import annotations

from hypothesis import given, strategies as st

from biobox import layout


def test_prefix_width_scales_with_stage_size():
    t = layout.TARGET_FILE_BYTES
    assert layout.prefix_width(0) == 0
    assert layout.prefix_width(8 * 2**20) == 0           # the golden stage: one file
    assert layout.prefix_width(16 * t - 1) == 0
    assert layout.prefix_width(16 * t) == 1
    assert layout.prefix_width(256 * t) == 2
    assert layout.prefix_width(10**15) == 2


@given(st.binary(min_size=16, max_size=16), st.sampled_from([0, 1, 2]))
def test_a_key_matches_its_own_bucket_and_no_other_at_that_width(key, width):
    b = layout.bucket(key, width)
    assert layout.bucket_matches(b, key)
    if width:
        assert len(b) == width
        others = [f"{i:02x}"[:width] for i in range(256)]
        assert sum(layout.bucket_matches(o, key) for o in set(others)) == 1


@given(st.lists(st.floats(min_value=1, max_value=500), max_size=200), st.integers(min_value=1, max_value=2_000))
def test_cut_respects_the_cap_unless_one_row_exceeds_it(sizes, cap):
    starts = layout.cut(sizes, cap)
    assert starts[0] == 0 and starts == sorted(set(starts))
    for a, b in zip(starts, starts[1:] + [len(sizes)]):
        part = sizes[a:b]
        assert sum(part) <= cap or len(part) == 1
