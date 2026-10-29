"""Ring R1: contract checker behaves on tiny in-memory tables."""

from __future__ import annotations

import pyarrow as pa

import pyarrow.parquet as pq
import pytest

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
