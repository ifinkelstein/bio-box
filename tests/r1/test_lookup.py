"""Ring R1: `KeyLookup`, a key lookup over sorted Parquet files that reads footers once and one
row group per key (card B-13; AGENTS §4.4; LY-06).

Oracle: a brute-force scan of every row of every file, written here without the lookup, and for
the row groups read, the row groups whose footer min/max holds the key, counted from the footers.
"""

from __future__ import annotations

import random
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from hypothesis import given, settings, strategies as st

from biobox.l0seq.lookup import KeyLookup


def _k(i: int) -> bytes:
    return i.to_bytes(8, "big")


def _write(path: Path, keys: list[bytes], row_group_size: int, tag: str, *, stats=True) -> Path:
    keys = sorted(keys)
    t = pa.table({
        "k": pa.array(keys, pa.binary(8)),
        "v": pa.array([f"{tag}:{i}" for i in range(len(keys))]),
        "n": pa.array(range(len(keys)), pa.int64()),
    })
    pq.write_table(t, path, row_group_size=row_group_size, write_statistics=["k"] if stats else False)
    return path


def _brute(files: list[Path], key: bytes) -> list[tuple]:
    out = []
    for f in files:
        for r in pq.read_table(f).to_pylist():
            if r["k"] == key:
                out.append((r["k"], r["v"], r["n"]))
    return sorted(out)


def _rows(rows: list[dict]) -> list[tuple]:
    return sorted((r["k"], r["v"], r["n"]) for r in rows)


def _candidates(files: list[Path], key: bytes) -> int:
    n = 0
    for f in files:
        md = pq.read_metadata(f)
        for g in range(md.num_row_groups):
            if md.row_group(g).num_rows:
                s = md.row_group(g).column(0).statistics
                n += s.min <= key <= s.max
    return n


@pytest.fixture
def files(tmp_path) -> list[Path]:
    rng = random.Random(20261109)
    # file a: six copies of key 100 at rows 48..53, across the row-group boundary at row 50
    a = [_k(2 * i) for i in range(48)] + [_k(100)] * 6 + [_k(200 + 2 * i) for i in range(146)]
    # file b: other keys, some also in a (several files may hold a key)
    b = sorted({_k(rng.randrange(0, 600)) for _ in range(300)} | {_k(100), _k(250)})
    return [_write(tmp_path / "a.parquet", a, 50, "a"), _write(tmp_path / "b.parquet", b, 37, "b")]


@pytest.mark.covers("LY-06")
def test_lookup_equals_a_brute_force_scan(files):
    lk = KeyLookup(files, "k", ["v", "n"])
    for i in range(700):
        assert _rows(lk.get(_k(i))) == _brute(files, _k(i)), i


@pytest.mark.covers("LY-06")
def test_boundary_keys_and_duplicates_spanning_two_row_groups(files):
    assert pq.read_metadata(files[0]).num_row_groups == 4
    lk = KeyLookup(files, "k", ["v", "n"])
    got = lk.get(_k(100))
    assert sorted(r["n"] for r in got if r["v"].startswith("a:")) == [48, 49, 50, 51, 52, 53]
    assert _rows(got) == _brute(files, _k(100))
    for f in files:  # every row group's min and max key, in every file
        md = pq.read_metadata(f)
        for g in range(md.num_row_groups):
            s = md.row_group(g).column(0).statistics
            for key in (s.min, s.max):
                assert _rows(lk.get(key)) == _brute(files, key)


@pytest.mark.covers("LY-06")
def test_absent_keys(files):
    lk = KeyLookup(files, "k", ["v"])
    assert lk.get(_k(10**9)) == [] and lk.row_groups_read == 0  # above every range: nothing read
    absent = next(i for i in range(1, 96, 2) if not _brute(files, _k(i)))
    assert lk.get(_k(absent)) == []
    assert lk.row_groups_read == _candidates(files, _k(absent))


@pytest.mark.covers("LY-06")
def test_counters_one_row_group_per_lookup_and_cache_hits_read_none(tmp_path):
    f = _write(tmp_path / "one.parquet", [_k(3 * i) for i in range(10_000)], 1_000, "x")
    lk = KeyLookup([f], "k", ["v", "n"], cache_row_groups=4)
    assert lk.footers_read == 1 and lk.row_groups_read == 0  # construction reads no data
    for i in random.Random(5).sample(range(10_000), 40):
        before = lk.row_groups_read
        rows = lk.get(_k(3 * i))
        assert [r["n"] for r in rows] == [i]
        assert lk.row_groups_read - before <= 1
        before = lk.row_groups_read
        assert lk.get(_k(3 * i)) == rows  # a cache hit reads nothing
        assert lk.row_groups_read == before
    assert lk.footers_read == 1
    # the LRU holds 4 row groups: after 5 distinct ones the first is read again, the last is not
    lk = KeyLookup([f], "k", cache_row_groups=4)
    for g in range(5):
        lk.get(_k(3 * (1_000 * g + 500)))
    assert lk.row_groups_read == 5
    lk.get(_k(3 * 500))
    assert lk.row_groups_read == 6
    lk.get(_k(3 * 4_500))
    assert lk.row_groups_read == 6


def test_only_the_asked_columns_are_returned(files):
    lk = KeyLookup(files, "k", ["n"])
    assert lk.get(_k(100)) and all(set(r) == {"k", "n"} for r in lk.get(_k(100)))


def test_row_group_without_statistics_is_refused(tmp_path):
    f = _write(tmp_path / "nostats.parquet", [_k(i) for i in range(10)], 5, "z", stats=False)
    with pytest.raises(ValueError, match="min/max"):
        KeyLookup([f], "k")


def test_file_not_sorted_by_key_is_refused(tmp_path):
    f = tmp_path / "unsorted.parquet"
    pq.write_table(pa.table({"k": pa.array([_k(5), _k(6), _k(1), _k(2)], pa.binary(8))}), f, row_group_size=2)
    with pytest.raises(ValueError, match="not sorted"):
        KeyLookup([f], "k")


@pytest.mark.covers("LY-06")
@settings(max_examples=60, deadline=None)
@given(
    data=st.lists(st.lists(st.integers(0, 60), max_size=80), min_size=1, max_size=3),
    rg=st.integers(1, 9),
    probes=st.lists(st.integers(0, 63), min_size=1, max_size=20),
)
def test_property_lookup_equals_brute_force(tmp_path_factory, data, rg, probes):
    d = tmp_path_factory.mktemp("prop")
    files = [_write(d / f"f{i}.parquet", [_k(x) for x in xs], rg, str(i)) for i, xs in enumerate(data)]
    lk = KeyLookup(files, "k", ["v", "n"], cache_row_groups=2)
    for p in probes:
        assert _rows(lk.get(_k(p))) == _brute(files, _k(p))
        one = KeyLookup(files, "k", cache_row_groups=0)
        one.get(_k(p))
        assert one.row_groups_read == _candidates(files, _k(p))  # only row groups whose range holds the key
