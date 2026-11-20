"""Ring R1: a multi-part table whose last part is under the small-file floor keeps no such part (card B-57, LY-04).

Oracles: `contracts.check_file_sizes` (LY-04), `check_sort_order` and `check_statistics` on the manifest
form, and the Parquet footers read directly for the row-group caps (AGENTS §4.6).
"""

from __future__ import annotations

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import contracts
from biobox.io import part_files, write_parts
from biobox.registry.schemas import table

TAPES = table("contig_tapes")
EST = 400_000.0  # per row; with MAX_BYTES two rows fill a part
MAX_BYTES = 1_000_000
RG_ROWS, MAX_SLOTS = 2, 8
KW = {"row_group_rows": RG_ROWS, "max_slots": MAX_SLOTS, "slots_column": "n_slots",
      "statistics": ["n_slots", "contig_len"]}


def uid(i: int) -> bytes:
    return i.to_bytes(16, "big")


def tape_row(i: int, n: int) -> dict:
    return {"contig_uid": uid(i), "contig_len": 100 * n + 50, "is_circular": False, "n_slots": n,
            "feature_kinds": [1] * n, "strands": [1] * n, "start_bps": [100 * k for k in range(n)],
            "end_bps": [100 * k + 90 for k in range(n)], "delta_bp_prev_slots": [10] * n,
            "junction_topology_prevs": [1] * n, "orfs_to_contig_starts": list(range(n)),
            "orfs_to_contig_ends": list(range(n))[::-1], "slot_flags": [0] * n,
            "protein_uids": [uid(10_000 + k) for k in range(n)], "feature_uids": [uid(20_000 + k) for k in range(n)]}


def rows(n_rows: int) -> list[dict]:
    return [tape_row(i, 5) for i in range(n_rows)]


def est(b: pa.RecordBatch) -> np.ndarray:
    return np.full(b.num_rows, EST)


def write(rs: list[dict], base, **kw):
    return write_parts([pa.Table.from_pylist(rs, schema=TAPES.schema)], TAPES, base, est, MAX_BYTES, **KW, **kw)


def manifest(root, paths) -> dict:
    return {"root": str(root), "files": [{"table": "contig_tapes", "path": p.relative_to(root).as_posix(),
                                          "partitions": {}, "bytes": p.stat().st_size} for p in paths]}


def read_rows(paths) -> list[dict]:
    return [r for p in paths for r in pq.read_table(p).to_pylist()]


@pytest.mark.covers("LY-04", "LY-02", "LY-03", "LY-08")
def test_a_tail_under_the_floor_is_coalesced_into_the_part_before_it(tmp_path):
    rs = rows(3)  # two rows fill part 0, one row is the tail
    plain = write(rs, tmp_path / "plain", small_file_floor=0)
    assert len(plain) == 2
    tail = plain[-1].stat().st_size
    floor = 500_000
    assert tail < floor < MAX_BYTES
    before = contracts.check_file_sizes(manifest(tmp_path, plain), min_bytes=floor)
    assert before, "the reproduction: a 2-part write with a tail under the floor trips LY-04"

    paths = write(rs, tmp_path / "fixed", small_file_floor=floor)
    assert [p.name for p in paths] == ["part-0.parquet"]
    assert part_files(tmp_path / "fixed") == paths and sorted((tmp_path / "fixed").iterdir()) == paths
    assert contracts.check_file_sizes(manifest(tmp_path, paths), min_bytes=floor) == []
    assert read_rows(paths) == rs
    m = manifest(tmp_path, paths)
    assert contracts.check_sort_order(m) == [] and contracts.check_statistics(m) == []
    md = pq.ParquetFile(paths[0]).metadata
    assert sum(md.row_group(g).num_rows for g in range(md.num_row_groups)) == 3
    for g in range(md.num_row_groups):
        assert md.row_group(g).num_rows <= RG_ROWS
    for b in pq.ParquetFile(paths[0]).iter_batches(batch_size=1):  # every row is 5 slots: one row per group
        assert b.num_rows == 1
    assert md.num_row_groups == 3


@pytest.mark.covers("LY-04")
def test_a_tail_at_or_above_the_floor_and_a_single_part_stay_byte_identical(tmp_path):
    rs = rows(3)
    plain = write(rs, tmp_path / "plain", small_file_floor=0)
    tail = plain[-1].stat().st_size
    kept = write(rs, tmp_path / "at", small_file_floor=tail)  # not under the floor: left alone
    assert [p.read_bytes() for p in kept] == [p.read_bytes() for p in plain]
    one = write(rows(2), tmp_path / "one_a", small_file_floor=10**9)
    one_b = write(rows(2), tmp_path / "one_b", small_file_floor=0)
    assert len(one) == 1 and one[0].read_bytes() == one_b[0].read_bytes()


@pytest.mark.covers("LY-04")
def test_five_parts_only_the_last_two_merge_and_nothing_is_left_behind(tmp_path):
    rs = rows(5)  # parts of 2, 2, 1 rows
    plain = write(rs, tmp_path / "plain", small_file_floor=0)
    assert len(plain) == 3
    paths = write(rs, tmp_path / "fixed", small_file_floor=500_000)
    assert [p.name for p in paths] == ["part-0.parquet", "part-1.parquet"]
    assert paths[0].read_bytes() == plain[0].read_bytes()
    assert read_rows(paths) == rs and len(pq.read_table(paths[1])) == 3
    assert sorted(p.name for p in (tmp_path / "fixed").iterdir()) == ["part-0.parquet", "part-1.parquet"]
    again = write(rs, tmp_path / "again", small_file_floor=500_000)
    assert [p.read_bytes() for p in again] == [p.read_bytes() for p in paths]
