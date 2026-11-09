"""Ring R1: contract checker behaves on tiny in-memory tables."""

from __future__ import annotations

import pyarrow as pa

import pyarrow.parquet as pq
import pytest
from blake3 import blake3
from hypothesis import assume, given, strategies as st

from biobox.ids import KIND, feature_uid, region_uid
from biobox.registry.schemas import TABLES, conformance_problems, table


@pytest.mark.parametrize("name", sorted(TABLES))
def test_registry_schema_survives_parquet_round_trip(name, tmp_path):
    # Parquet cannot hold some Arrow types (e.g. timestamp[s]); a schema that changes on
    # write would fail C-01 on every file we produce.
    spec = TABLES[name]
    path = tmp_path / f"{name}.parquet"
    pq.write_table(spec.schema.empty_table(), path)
    assert pq.read_table(path).schema.equals(spec.schema, check_metadata=False)


def _igr(strands):
    spec = table("igr")
    n = len(strands)
    return spec, pa.table(
        {
            "contig_uid": pa.array([bytes(16)] * n, pa.binary(16)),
            "start": pa.array(range(0, 10 * n, 10), pa.int64()),
            "end": pa.array(range(5, 10 * n + 5, 10), pa.int64()),
            "strand": pa.array(strands, pa.int8()),
            "left_slot": pa.array([None] * n, pa.uint32()),
            "right_slot": pa.array([0] * n, pa.uint32()),
        },
        schema=spec.schema,
    )


def test_conforming_table_has_no_problems():
    spec, t = _igr([1, -1, 1])
    assert conformance_problems(t, spec) == []


def test_enum_violation_is_reported():
    spec, t = _igr([1, 0])
    assert any("strand" in p for p in conformance_problems(t, spec))


def test_duplicate_key_is_reported():
    spec, t = _igr([1])
    assert any("duplicate" in p for p in conformance_problems(pa.concat_tables([t, t]), spec))


def test_wrong_uid_width_is_a_schema_problem():
    spec, t = _igr([1])
    bad = t.set_column(0, pa.field("contig_uid", pa.binary(8), nullable=False), pa.array([bytes(8)], pa.binary(8)))
    assert conformance_problems(bad, spec) == ["igr: schema differs from registry"]


_UID = st.binary(min_size=16, max_size=16)
_COORD = st.integers(min_value=0, max_value=2**40)
_FEATURE = st.tuples(
    _UID, _COORD, st.integers(min_value=1, max_value=2**40), st.sampled_from([1, -1]), st.sampled_from(sorted(KIND.values()))
).map(lambda t: (t[0], t[1], t[1] + t[2], t[3], t[4]))


@pytest.mark.covers("LY-01")
def test_kind_values():
    assert KIND == {"cds": 1, "ncrna": 2, "trna": 3, "crispr_array": 4, "other": 15}


@pytest.mark.covers("LY-01")
def test_feature_uid_known_answer():
    contig = bytes(range(16))
    layout = contig + (100).to_bytes(8, "big") + (400).to_bytes(8, "big") + b"\xff" + b"\x02"
    assert len(layout) == 34
    expected = blake3(layout).digest(length=16)
    assert feature_uid(contig, 100, 400, -1, KIND["ncrna"]) == expected
    assert len(feature_uid(contig, 100, 400, -1, 2)) == 16


@pytest.mark.covers("LY-01")
@given(_FEATURE, _FEATURE)
def test_feature_uid_distinct_inputs_give_distinct_uids(a, b):
    assume(a != b)
    assert feature_uid(*a) != feature_uid(*b)


@pytest.mark.covers("LY-01")
@given(_FEATURE)
def test_feature_uid_strand_kind_and_region_differ(f):
    contig, start, end, strand, kind = f
    assert feature_uid(contig, start, end, strand, kind) != feature_uid(contig, start, end, -strand, kind)
    other_kind = next(k for k in sorted(KIND.values()) if k != kind)
    assert feature_uid(contig, start, end, strand, kind) != feature_uid(contig, start, end, strand, other_kind)
    assert feature_uid(contig, start, end, strand, kind) != region_uid(contig, start, end, strand)


@pytest.mark.covers("LY-01")
def test_feature_uid_allows_origin_wrapping_end():
    assert len(feature_uid(bytes(16), 4_000_000, 4_700_000, 1, 1)) == 16


@pytest.mark.covers("LY-01")
@pytest.mark.parametrize(
    "args",
    [
        (bytes(15), 0, 5, 1, 1),
        (bytes(17), 0, 5, 1, 1),
        (bytes(16), -1, 5, 1, 1),
        (bytes(16), 5, 5, 1, 1),
        (bytes(16), 5, 4, 1, 1),
        (bytes(16), 0, 5, 0, 1),
        (bytes(16), 0, 5, 2, 1),
        (bytes(16), 0, 5, 1, 0),
        (bytes(16), 0, 5, 1, 5),
        (bytes(16), 0, 5, 1, 256),
    ],
)
def test_feature_uid_rejects_invalid_input(args):
    with pytest.raises(ValueError):
        feature_uid(*args)
