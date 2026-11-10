"""Ring R1: the one Parquet writer (biobox.io, card B-01; LY-01..LY-04, LY-08, F-11, R-09).

Oracles: `contracts._sort_ok` (a row-by-row Python comparison) for the sort check,
`contracts.check_statistics` / `check_sort_order` for the footer, a separate greedy loop for the
slot-bounded row groups, and the registry's `conformance_problems` for the schema.
"""

from __future__ import annotations

import dataclasses
import os
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from hypothesis import given, settings, strategies as st

from biobox import contracts
from biobox.io import duck_batches, write_table
from biobox.registry.schemas import conformance_problems, table

OCC = table("protein_occurrence")
TAPES = table("contig_tapes")
# the tapes' sort key before card B-07, kept here for the writer's nulls-first order on a nullable key
NULL_KEYED = dataclasses.replace(TAPES, sort_key=("tax_order_bucket", "contig_uid"))
DOM = table("l1_domains")
TOK = table("token_dict")


def uid(i: int) -> bytes:
    return i.to_bytes(16, "big")


def occ_row(key: bytes, i: int) -> dict:
    return {"protein_uid": key, "contig_uid": uid(7 * i + 1), "slot_ordinal": i, "genome_uid": uid(3),
            "start": 100 * i, "end": 100 * i + 90, "strand": 1 if i % 2 else -1,
            "tax_order": None if i % 3 == 0 else i, "is_species_rep": i % 2 == 0, "tier": 1,
            "assembly_quality": "complete", "dist_to_contig_end_bp": i, "feature_uid": uid(11 * i + 5)}


def occ_rows(n: int, first: int = 0) -> list[dict]:
    return [occ_row(uid(first + i), first + i) for i in range(n)]


def tape_row(i: int, n: int, bucket: int | None = None) -> dict:
    return {"contig_uid": uid(i), "tax_order_bucket": bucket, "genome_uid": None, "contig_len": 100 * n + 50,
            "is_circular": False, "n_slots": n, "tier": 1, "assembly_quality": "complete",
            "feature_kinds": [1] * n, "strands": [1] * n, "start_bps": [100 * k for k in range(n)],
            "end_bps": [100 * k + 90 for k in range(n)], "delta_bp_prev_slots": [10] * n,
            "junction_topology_prevs": [1] * n, "orfs_to_contig_starts": list(range(n)),
            "orfs_to_contig_ends": list(range(n))[::-1], "slot_flags": [0] * n,
            "protein_uids": [uid(10_000 + k) for k in range(n)], "feature_uids": [uid(20_000 + k) for k in range(n)]}


def dom_row(key: tuple[bytes, int, str]) -> dict:
    p, env, acc = key
    return {"protein_uid": p, "model_acc": acc, "model_name": acc.lower(), "model_len": 200,
            "seq_bitscore": 50.0, "seq_evalue_zref": 1e-10, "dom_bitscore": 40.0, "dom_ievalue_zref": 1e-9,
            "env_start": env, "env_end": env + 30, "ali_start": env, "ali_end": env + 30, "hmm_from": 0,
            "hmm_to": 30, "cutoff": "gathering", "z": 1e8, "dom_z": 1e8}


def inferred(rows: list[dict]) -> pa.Table:
    """Rows as pyarrow infers them: bytes as variable binary, ints as int64 (not the spec types)."""
    return pa.Table.from_pylist(rows)


def batches_of(tbl: pa.Table, cuts) -> list[pa.RecordBatch]:
    bounds = [0, *sorted(cuts), tbl.num_rows]
    out = []
    for a, b in zip(bounds, bounds[1:]):
        part = tbl.slice(a, b - a).combine_chunks().to_batches()
        out.append(part[0] if part else pa.RecordBatch.from_pylist([], schema=tbl.schema))
    return out


def rg_sizes(path: Path) -> list[int]:
    md = pq.ParquetFile(path).metadata
    return [md.row_group(i).num_rows for i in range(md.num_row_groups)]


def manifest(root: Path, name: str, path: Path) -> dict:
    return {"root": str(root), "files": [{"table": name, "path": path.relative_to(root).as_posix()}]}


# --- cast (LY-01) ------------------------------------------------------------------------------


@pytest.mark.covers("LY-01")
def test_inferred_types_are_cast_to_the_spec_including_fixed_size_binary(tmp_path):
    rows = occ_rows(50)
    src = inferred(rows)
    assert src.schema.field("protein_uid").type == pa.binary()
    out = tmp_path / "a/b/part-0.parquet"
    res = write_table([src], OCC, out)
    got = pq.read_table(out)
    assert got.schema.equals(OCC.schema, check_metadata=False)
    assert pq.ParquetFile(out).schema_arrow.field("protein_uid").type == pa.binary(16)
    assert conformance_problems(got, OCC) == []
    assert got.to_pylist() == rows
    assert res["rows"] == 50 and res["row_groups"] == 1 and res["bytes"] == out.stat().st_size


@pytest.mark.covers("LY-01")
def test_list_of_binary_is_cast_to_list_of_uids(tmp_path):
    rows = [tape_row(i, 3) for i in range(5)]
    out = tmp_path / "t.parquet"
    write_table([inferred(rows)], TAPES, out)
    got = pq.read_table(out)
    assert got.schema.equals(TAPES.schema, check_metadata=False)
    assert got.to_pylist() == rows


@pytest.mark.covers("LY-01")
def test_wrong_uid_length_fails_loudly_and_writes_nothing(tmp_path):
    rows = occ_rows(5)
    rows[3]["contig_uid"] = b"\x01" * 15
    with pytest.raises(ValueError, match="contig_uid"):
        write_table([inferred(rows)], OCC, tmp_path / "x.parquet")
    assert list(tmp_path.iterdir()) == []


@pytest.mark.covers("LY-01")
def test_columns_are_reordered_to_the_spec(tmp_path):
    src = inferred(occ_rows(10))
    src = src.select(list(reversed(src.column_names)))
    out = tmp_path / "x.parquet"
    write_table([src], OCC, out)
    assert pq.read_table(out).schema.names == OCC.schema.names


@pytest.mark.covers("LY-01")
@pytest.mark.parametrize("change,name", [("drop", "tier"), ("add", "bogus")])
def test_missing_or_extra_column_is_an_error_naming_it(tmp_path, change, name):
    src = inferred(occ_rows(4))
    src = src.drop_columns([name]) if change == "drop" else src.append_column(name, pa.array([1, 2, 3, 4]))
    with pytest.raises(ValueError, match=name):
        write_table([src], OCC, tmp_path / "x.parquet")
    assert list(tmp_path.iterdir()) == []


@pytest.mark.covers("LY-01")
def test_null_in_a_non_nullable_field_is_an_error(tmp_path):
    rows = occ_rows(4)
    rows[2]["contig_uid"] = None
    with pytest.raises(ValueError, match="contig_uid"):
        write_table([inferred(rows)], OCC, tmp_path / "x.parquet")
    rows = occ_rows(4)
    rows[2]["tax_order"] = None  # nullable: fine
    write_table([inferred(rows)], OCC, tmp_path / "y.parquet")


@pytest.mark.covers("LY-01")
def test_out_of_range_integer_is_an_error(tmp_path):
    rows = occ_rows(4)
    rows[1]["tier"] = 300  # uint8
    with pytest.raises(ValueError, match="tier"):
        write_table([inferred(rows)], OCC, tmp_path / "x.parquet")


# --- sort refusal (LY-02) ----------------------------------------------------------------------


@pytest.mark.covers("LY-02")
def test_unsorted_rows_within_a_batch_are_refused_with_the_row_index(tmp_path):
    rows = occ_rows(10)
    rows[4], rows[5] = rows[5], rows[4]
    with pytest.raises(ValueError, match=r"protein_occurrence.*protein_uid.*row 5"):
        write_table([inferred(rows)], OCC, tmp_path / "x.parquet")
    assert list(tmp_path.iterdir()) == []


@pytest.mark.covers("LY-02")
def test_unsorted_rows_across_a_batch_boundary_are_refused(tmp_path):
    a = inferred(occ_rows(6, first=10))  # keys 10..15
    write_table([a, inferred(occ_rows(6, first=16))], OCC, tmp_path / "ok.parquet")
    write_table([a, inferred(occ_rows(6, first=15))], OCC, tmp_path / "tie.parquet")  # 15 then 15 is allowed
    c = inferred(occ_rows(6, first=12))  # 12 after 15
    with pytest.raises(ValueError, match=r"protein_occurrence.*row 6\b"):
        write_table([a, c], OCC, tmp_path / "x.parquet", row_group_rows=4)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["ok.parquet", "tie.parquet"]


@pytest.mark.covers("LY-02")
def test_unsorted_after_empty_batches_is_still_caught(tmp_path):
    a, c = inferred(occ_rows(3, first=10)), inferred(occ_rows(3, first=5))
    empty = a.slice(0, 0)
    with pytest.raises(ValueError, match=r"row 3\b"):
        write_table([a, empty, empty, c], OCC, tmp_path / "x.parquet")


@pytest.mark.covers("LY-02")
def test_nulls_sort_first(tmp_path):
    ok = [tape_row(5, 1, None), tape_row(9, 1, None), tape_row(1, 1, 0), tape_row(2, 1, 3)]
    write_table([inferred(ok)], NULL_KEYED, tmp_path / "ok.parquet")
    bad = [tape_row(1, 1, 0), tape_row(2, 1, None)]
    with pytest.raises(ValueError, match=r"contig_tapes.*tax_order_bucket.*row 1\b"):
        write_table([inferred(bad)], NULL_KEYED, tmp_path / "bad.parquet")
    across = [inferred([tape_row(1, 1, 4)]), inferred([tape_row(2, 1, None)])]
    with pytest.raises(ValueError, match=r"row 1\b"):
        write_table(across, NULL_KEYED, tmp_path / "across.parquet")
    # within the null bucket the second key decides
    with pytest.raises(ValueError, match=r"row 1\b"):
        write_table([inferred([tape_row(9, 1, None), tape_row(5, 1, None)])], NULL_KEYED, tmp_path / "n.parquet")


@pytest.mark.covers("LY-02")
def test_binary_keys_compare_as_unsigned_bytes(tmp_path):
    lo, hi = b"\x01" + b"\x00" * 15, b"\x80" + b"\x00" * 15
    write_table([inferred([occ_row(lo, 0), occ_row(hi, 1)])], OCC, tmp_path / "ok.parquet")
    with pytest.raises(ValueError, match=r"row 1\b"):
        write_table([inferred([occ_row(hi, 0), occ_row(lo, 1)])], OCC, tmp_path / "bad.parquet")
    assert contracts._sort_ok([(hi,), (lo,)]) == 1  # the oracle agrees


_KEY = st.tuples(st.sampled_from([uid(1), uid(2), b"\x80" + bytes(15), b"\xff" * 16]),
                 st.integers(0, 3), st.sampled_from(["PF1", "PF2", "PFa"]))


@pytest.mark.covers("LY-02")
@settings(max_examples=60, deadline=None)
@given(keys=st.lists(_KEY, max_size=40), sort_it=st.booleans(), cuts=st.lists(st.integers(0, 40), max_size=5),
       rg=st.integers(1, 8))
def test_refusal_agrees_with_the_row_by_row_oracle(keys, sort_it, cuts, rg):
    if sort_it:
        keys = sorted(keys)
    rows = [dom_row(k) for k in keys]
    want = contracts._sort_ok([(r["protein_uid"], r["env_start"], r["model_acc"]) for r in rows])
    src = inferred(rows) if rows else DOM.schema.empty_table()
    with tempfile.TemporaryDirectory() as d:
        out = Path(d) / "x.parquet"
        batches = batches_of(src, [c for c in cuts if c <= len(rows)])
        if want is None:
            assert write_table(batches, DOM, out, row_group_rows=rg)["rows"] == len(rows)
        else:
            with pytest.raises(ValueError, match=rf"row {want}\b"):
                write_table(batches, DOM, out, row_group_rows=rg)
            assert list(Path(d).iterdir()) == []


@pytest.mark.covers("LY-02")
@settings(max_examples=60, deadline=None)
@given(keys=st.lists(st.tuples(st.sampled_from([None, 0, 1, 7]), st.sampled_from([1, 2, 0x80 << 120])), max_size=30),
       sort_it=st.booleans(), cuts=st.lists(st.integers(0, 30), max_size=5))
def test_refusal_with_a_nullable_key_agrees_with_the_oracle(keys, sort_it, cuts):
    if sort_it:
        keys = sorted(keys, key=lambda k: ((0,) if k[0] is None else (1, k[0]), k[1]))
    rows = [tape_row(c, 1, b) for b, c in keys]
    want = contracts._sort_ok([(r["tax_order_bucket"], r["contig_uid"]) for r in rows])
    src = pa.Table.from_pylist(rows, schema=TAPES.schema)
    with tempfile.TemporaryDirectory() as d:
        batches = batches_of(src, [c for c in cuts if c <= len(rows)])
        if want is None:
            write_table(batches, NULL_KEYED, Path(d) / "x.parquet", row_group_rows=4)
        else:
            with pytest.raises(ValueError, match=rf"row {want}\b"):
                write_table(batches, NULL_KEYED, Path(d) / "x.parquet", row_group_rows=4)


# --- row groups and bytes (LY-04, LY-08)-------------------------------------------------------


@pytest.mark.covers("LY-08")
def test_one_batch_many_batches_and_a_table_give_identical_bytes(tmp_path):
    src = inferred(occ_rows(1000))
    cuts = [1, 2, 3, 50, 51, 127, 128, 129, 400, 999]
    inputs = {
        "one": [src.combine_chunks().to_batches()[0]],
        "many": batches_of(src, cuts),
        "table": [pa.Table.from_batches(batches_of(src, cuts))],
        "bare_table": pa.Table.from_batches(batches_of(src, [7, 300])),
        "spec_typed": [pa.Table.from_pylist(occ_rows(1000), schema=OCC.schema)],
    }
    digests = {}
    for name, batches in inputs.items():
        out = tmp_path / f"{name}.parquet"
        write_table(batches, OCC, out, row_group_rows=128)
        digests[name] = out.read_bytes()
    assert len(set(digests.values())) == 1


def _big_occurrence(n: int, seed: int = 1) -> pa.Table:
    """n spec-typed rows sorted by protein_uid, with random (incompressible, dictionary-busting) UIDs."""
    rng = np.random.default_rng(seed)
    keys = np.zeros((n, 16), dtype=np.uint8)
    keys[:, 8:] = np.frombuffer(np.sort(rng.integers(0, 2**63, n)).astype(">u8").tobytes(), np.uint8).reshape(n, 8)
    i = np.arange(n)
    cols = {
        "protein_uid": pa.FixedSizeBinaryArray.from_buffers(pa.binary(16), n, [None, pa.py_buffer(keys.tobytes())]),
        "contig_uid": pa.FixedSizeBinaryArray.from_buffers(pa.binary(16), n, [None, pa.py_buffer(rng.bytes(16 * n))]),
        "slot_ordinal": i % 5000, "genome_uid": pa.array([uid(3)] * n, pa.binary(16)), "start": i * 100,
        "end": i * 100 + 90, "strand": np.where(i % 2, 1, -1), "tax_order": rng.integers(0, 1000, n),
        "is_species_rep": i % 2 == 0, "tier": np.ones(n, np.uint8), "assembly_quality": ["complete"] * n,
        "dist_to_contig_end_bp": i % 7000,
        "feature_uid": pa.FixedSizeBinaryArray.from_buffers(pa.binary(16), n, [None, pa.py_buffer(rng.bytes(16 * n))]),
    }
    return pa.table({f.name: pa.array(cols[f.name]).cast(f.type) for f in OCC.schema}, schema=OCC.schema)


@pytest.mark.covers("LY-08")
def test_bytes_do_not_depend_on_batching_when_pages_and_dictionaries_overflow(tmp_path):
    src = _big_occurrence(200_000)
    write_table([src], OCC, tmp_path / "one.parquet")
    cuts = sorted(np.random.default_rng(2).choice(np.arange(1, 200_000), 150, replace=False).tolist())
    write_table(batches_of(src, cuts), OCC, tmp_path / "split.parquet")
    assert rg_sizes(tmp_path / "one.parquet") == [122_880, 77_120]
    assert (tmp_path / "one.parquet").read_bytes() == (tmp_path / "split.parquet").read_bytes()


@pytest.mark.covers("LY-04")
@pytest.mark.parametrize("n,want", [(1000, [128] * 7 + [104]), (1024, [128] * 8), (5, [5]), (128, [128])])
def test_every_row_group_has_row_group_rows_except_the_last(tmp_path, n, want):
    out = tmp_path / "x.parquet"
    res = write_table(batches_of(inferred(occ_rows(n)), [1, 2, 3, 200, 201]), OCC, out, row_group_rows=128)
    assert rg_sizes(out) == want
    assert res == {"rows": n, "row_groups": len(want), "bytes": out.stat().st_size}


def _greedy_groups(slots: list[int], max_slots: int, max_rows: int) -> list[int]:
    """Oracle: a separate loop over rows. Close before a row would pass the cap; one oversized row alone."""
    sizes, cur_rows, cur_slots = [], 0, 0
    for s in slots:
        if cur_rows and (cur_slots + s > max_slots or cur_rows == max_rows):
            sizes.append(cur_rows)
            cur_rows, cur_slots = 0, 0
        cur_rows += 1
        cur_slots += s
    if cur_rows:
        sizes.append(cur_rows)
    return sizes


@pytest.mark.covers("LY-04")
def test_slot_bounded_row_groups_including_an_oversized_contig(tmp_path):
    slots = [3, 4, 2, 6, 1, 1, 1, 12, 2, 8, 0, 1]
    rows = [tape_row(i, n) for i, n in enumerate(slots)]
    src = inferred(rows)
    out = tmp_path / "one.parquet"
    write_table([src], TAPES, out, max_slots=8, slots_column="n_slots")
    assert rg_sizes(out) == [2, 2, 3, 1, 1, 2, 1]  # [3,4] [2,6] [1,1,1] [12] [2] [8,0] [1]
    assert rg_sizes(out) == _greedy_groups(slots, 8, 10**9)
    many = tmp_path / "many.parquet"
    write_table(batches_of(src, [1, 2, 5, 7, 8, 9]), TAPES, many, max_slots=8, slots_column="n_slots")
    assert many.read_bytes() == out.read_bytes()
    capped = tmp_path / "capped.parquet"
    write_table(batches_of(src, [3, 4]), TAPES, capped, row_group_rows=2, max_slots=8, slots_column="n_slots")
    assert rg_sizes(capped) == _greedy_groups(slots, 8, 2) == [2, 2, 2, 1, 1, 1, 2, 1]


@pytest.mark.covers("LY-04")
@settings(max_examples=60, deadline=None)
@given(slots=st.lists(st.integers(0, 20), max_size=40), cap=st.integers(1, 25), rows=st.integers(1, 6),
       cuts=st.lists(st.integers(0, 40), max_size=6))
def test_slot_bounded_row_groups_match_the_greedy_oracle(slots, cap, rows, cuts):
    src = pa.Table.from_pylist([tape_row(i, n) for i, n in enumerate(slots)], schema=TAPES.schema)
    with tempfile.TemporaryDirectory() as d:
        out = Path(d) / "x.parquet"
        write_table(batches_of(src, [c for c in cuts if c <= len(slots)]), TAPES, out, row_group_rows=rows,
                    max_slots=cap, slots_column="n_slots")
        assert rg_sizes(out) == _greedy_groups(slots, cap, rows)


def tok_row(i: int, n: int) -> dict:
    return {"contig_uid": uid(i), "family_tokens": [None if k % 3 == 0 else k for k in range(n)],
            "nr50_tokens": list(range(n)), "nr90_tokens": [7] * n}


@pytest.mark.covers("LY-04")
def test_a_list_slots_column_counts_its_length_as_the_slots(tmp_path):
    """card B-09: `tape_tokens` has no `n_slots`; the length of one of its lists is the row's slot count."""
    spec = table("tape_tokens")
    slots = [3, 4, 2, 6, 1, 1, 1, 12, 2, 8, 0, 1]
    rows = [tok_row(i, n) for i, n in enumerate(slots)]
    src = inferred(rows)
    out = tmp_path / "one.parquet"
    write_table([src], spec, out, max_slots=8, slots_column="family_tokens")  # nulls inside a list are slots too
    assert rg_sizes(out) == _greedy_groups(slots, 8, 10**9) == [2, 2, 3, 1, 1, 2, 1]
    assert pq.read_table(out).to_pylist() == rows
    many = tmp_path / "many.parquet"
    write_table(batches_of(src, [1, 2, 5, 7, 8, 9]), spec, many, max_slots=8, slots_column="nr50_tokens")
    assert many.read_bytes() == out.read_bytes()
    capped = tmp_path / "capped.parquet"
    write_table(batches_of(src, [3, 4]), spec, capped, row_group_rows=2, max_slots=8, slots_column="nr90_tokens")
    assert rg_sizes(capped) == _greedy_groups(slots, 8, 2)
    # the same cuts as the integer column that holds the lengths
    tapes = tmp_path / "tapes.parquet"
    write_table([inferred([tape_row(i, n) for i, n in enumerate(slots)])], TAPES, tapes, max_slots=8,
                slots_column="n_slots")
    assert rg_sizes(tapes) == rg_sizes(out)


def test_a_slots_column_must_be_an_integer_or_a_non_null_list(tmp_path):
    with pytest.raises(ValueError, match="assembly_quality"):
        write_table([inferred([tape_row(0, 2)])], TAPES, tmp_path / "a.parquet", max_slots=8,
                    slots_column="assembly_quality")
    circular = {**tape_row(1, 2), "orfs_to_contig_starts": None}  # a nullable list column holding a null list
    with pytest.raises(ValueError, match="orfs_to_contig_starts"):
        write_table([inferred([tape_row(0, 2), circular])], TAPES, tmp_path / "b.parquet", max_slots=8,
                    slots_column="orfs_to_contig_starts")
    assert list(tmp_path.iterdir()) == []


@pytest.mark.covers("LY-08")
@settings(max_examples=40, deadline=None)
@given(keys=st.lists(_KEY, max_size=300), cuts=st.lists(st.integers(0, 300), max_size=12),
       rg=st.integers(1, 64))
def test_bytes_do_not_depend_on_how_the_input_is_batched(keys, cuts, rg):
    rows = [dom_row(k) for k in sorted(keys)]
    src = pa.Table.from_pylist(rows, schema=DOM.schema)
    with tempfile.TemporaryDirectory() as d:
        one, split = Path(d) / "one.parquet", Path(d) / "split.parquet"
        write_table([src], DOM, one, row_group_rows=rg)
        write_table(batches_of(inferred(rows) if rows else src, [c for c in cuts if c <= len(rows)]), DOM, split,
                    row_group_rows=rg)
        assert one.read_bytes() == split.read_bytes()
        assert pq.read_table(split).to_pylist() == rows


# --- footer (LY-03, LY-08) ---------------------------------------------------------------------


@pytest.mark.covers("LY-03", "LY-02")
def test_sort_key_columns_have_min_max_in_every_row_group(tmp_path):
    rows = [tape_row(i, 2, None) for i in range(10)] + [tape_row(i, 2, 7) for i in range(10, 20)]
    out = tmp_path / "l3/tapes/tape_version=t-000000000000/stage=S0/part-0.parquet"
    write_table([inferred(rows)], TAPES, out, row_group_rows=4)
    m = manifest(tmp_path, "contig_tapes", out)
    assert contracts.check_statistics(m) == []
    assert contracts.check_sort_order(m) == []
    md = pq.ParquetFile(out).metadata
    leaf = {md.schema.column(i).path: i for i in range(md.num_columns)}
    for rg in range(md.num_row_groups):
        st_ = md.row_group(rg).column(leaf["contig_uid"]).statistics
        assert st_.has_min_max
        assert md.row_group(rg).column(leaf["contig_len"]).statistics is None  # not requested
        assert md.row_group(rg).column(leaf["contig_uid"]).compression == "ZSTD"


@pytest.mark.covers("LY-03")
def test_extra_statistics_resolve_list_columns_to_their_leaf(tmp_path):
    out = tmp_path / "x.parquet"
    write_table([inferred([tape_row(i, 3) for i in range(6)])], TAPES, out,
                statistics=["contig_len", "protein_uids"])
    md = pq.ParquetFile(out).metadata
    leaf = {md.schema.column(i).path: i for i in range(md.num_columns)}
    rg = md.row_group(0)
    assert rg.column(leaf["contig_len"]).statistics.has_min_max
    [path] = [p for p in leaf if p.startswith("protein_uids.")]
    assert rg.column(leaf[path]).statistics.has_min_max
    assert rg.column(leaf["n_slots"]).statistics is None
    with pytest.raises(ValueError, match="nonesuch"):
        write_table([inferred([tape_row(0, 1)])], TAPES, tmp_path / "y.parquet", statistics=["nonesuch"])


@pytest.mark.covers("LY-08")
def test_footer_holds_no_run_dependent_metadata(tmp_path):
    out = tmp_path / "x.parquet"
    write_table([inferred(occ_rows(3))], OCC, out)
    md = pq.ParquetFile(out).metadata
    assert md.created_by.startswith("parquet-cpp-arrow version ")
    assert set(md.metadata) == {b"ARROW:schema"}  # pyarrow's own copy of the spec schema; nothing added
    for needle in (str(tmp_path).encode(), os.uname().nodename.encode()):
        assert needle not in out.read_bytes()


# --- atomic writes (F-11, R-09) ----------------------------------------------------------------


@pytest.mark.covers("F-11")
def test_a_failing_input_leaves_neither_the_file_nor_a_temporary(tmp_path):
    d = tmp_path / "out"
    seen = []

    def source():
        yield inferred(occ_rows(10))
        yield inferred(occ_rows(10, first=10))
        seen.extend(p.name for p in d.iterdir())  # two row groups are on disk under the temporary name
        raise RuntimeError("source died")

    with pytest.raises(RuntimeError, match="source died"):
        write_table(source(), OCC, d / "part-0.parquet", row_group_rows=4)
    assert seen == [f".part-0.parquet.tmp-{os.getpid()}"]
    assert list(d.iterdir()) == []


@pytest.mark.covers("F-11")
def test_a_kill_before_the_rename_leaves_no_final_file(tmp_path):
    out = tmp_path / "part-0.parquet"
    script = textwrap.dedent(f"""
        import os, signal
        import pyarrow as pa
        from biobox.io import write_table
        from biobox.registry.schemas import table
        def source():
            yield pa.Table.from_pylist([{{"token_id": i, "level": "nr90", "source": "s", "key": str(i),
                                          "n_slots": 1}} for i in range(100)])
            os.kill(os.getpid(), signal.SIGKILL)
        write_table(source(), table("token_dict"), {str(out)!r}, row_group_rows=10)
    """)
    proc = subprocess.run([sys.executable, "-c", script], capture_output=True, timeout=120)
    assert proc.returncode == -9, proc.stderr.decode()
    assert not out.exists()
    left = [p.name for p in tmp_path.iterdir()]
    assert len(left) == 1 and left[0].startswith(".part-0.parquet.tmp-") and not left[0].endswith(".parquet")


@pytest.mark.covers("R-09")
def test_an_existing_file_is_never_overwritten(tmp_path):
    out = tmp_path / "x.parquet"
    write_table([inferred(occ_rows(3))], OCC, out)
    before = out.read_bytes()
    with pytest.raises(FileExistsError):
        write_table([inferred(occ_rows(5))], OCC, out)
    assert out.read_bytes() == before
    assert [p.name for p in tmp_path.iterdir()] == ["x.parquet"]


@pytest.mark.covers("LY-01")
@pytest.mark.parametrize("batches", [[], "empty_batches"])
def test_empty_input_writes_a_valid_empty_file(tmp_path, batches):
    if batches == "empty_batches":
        batches = [inferred(occ_rows(2)).slice(0, 0), pa.RecordBatch.from_pylist([], schema=OCC.schema)]
    out = tmp_path / "x.parquet"
    assert write_table(batches, OCC, out)["rows"] == 0
    got = pq.read_table(out)
    assert got.num_rows == 0 and got.schema.equals(OCC.schema, check_metadata=False)
    assert contracts.check_statistics(manifest(tmp_path, "protein_occurrence", out)) == []


# --- DuckDB source -----------------------------------------------------------------------------

_OCC_SQL = """
SELECT CAST(i AS BIGINT) AS dist_to_contig_end_bp,
       from_hex(lpad(to_hex(i), 32, '0')) AS protein_uid,
       from_hex(lpad(to_hex(7 * i + 1), 32, '0')) AS contig_uid,
       i AS slot_ordinal, from_hex(lpad('3', 32, '0')) AS genome_uid,
       from_hex(lpad(to_hex(11 * i + 5), 32, '0')) AS feature_uid,
       100 * i AS start, 100 * i + 90 AS "end", CASE WHEN i % 2 = 1 THEN 1 ELSE -1 END AS strand,
       CASE WHEN i % 3 = 0 THEN NULL ELSE i END AS tax_order, i % 2 = 0 AS is_species_rep, 1 AS tier,
       'complete' AS assembly_quality
FROM range($n) t(i) ORDER BY protein_uid
"""


@pytest.mark.covers("LY-01", "LY-08")
def test_duck_batches_stream_cast_batches_into_write_table(tmp_path):
    con = duckdb.connect()
    batches = list(duck_batches(con, _OCC_SQL, OCC, params={"n": 1000}, batch_rows=128))
    assert len(batches) > 1 and all(b.num_rows <= 128 for b in batches)
    assert all(b.schema.equals(OCC.schema) for b in batches)
    out = tmp_path / "duck.parquet"
    write_table(duck_batches(con, _OCC_SQL, OCC, params={"n": 1000}, batch_rows=100), OCC, out, row_group_rows=256)
    got = pq.read_table(out)
    assert conformance_problems(got, OCC) == []
    assert got.to_pylist() == occ_rows(1000)
    ref = tmp_path / "ref.parquet"
    write_table([inferred(occ_rows(1000))], OCC, ref, row_group_rows=256)
    assert out.read_bytes() == ref.read_bytes()


def test_duck_batches_is_lazy_and_names_a_missing_column():
    con = duckdb.connect()
    gen = duck_batches(con, "SELECT * FROM no_such_table", OCC)  # nothing runs until the first batch
    with pytest.raises(duckdb.Error):
        next(gen)
    with pytest.raises(ValueError, match="assembly_quality"):
        next(duck_batches(con, _OCC_SQL.replace(",\n       'complete' AS assembly_quality", ""), OCC,
                          params={"n": 10}))


def test_duck_batches_cast_blob_lists_to_uid_lists():
    con = duckdb.connect()
    sql = ("SELECT from_hex(lpad(to_hex(i), 32, '0')) AS contig_uid, "
           "[from_hex(lpad(to_hex(i), 16, '0')), from_hex(lpad(to_hex(i + 1), 16, '0'))] AS family_tokens, "
           "[1, 2] AS nr50_tokens, [3, 4] AS nr90_tokens FROM range(3) t(i) ORDER BY 1")
    spec = table("tape_tokens")
    with pytest.raises(ValueError, match="family_tokens"):
        list(duck_batches(con, sql, spec))  # list<blob> is not list<uint64>
    good = sql.replace("[from_hex(lpad(to_hex(i), 16, '0')), from_hex(lpad(to_hex(i + 1), 16, '0'))]", "[i, i]")
    [b] = list(duck_batches(con, good, spec))
    assert b.schema.equals(spec.schema)
    occ_sql = _OCC_SQL.replace("from_hex(lpad('3', 32, '0'))", "from_hex(lpad('3', 30, '0'))")
    with pytest.raises(ValueError, match="genome_uid"):
        list(duck_batches(con, occ_sql, OCC, params={"n": 5}))


def test_token_dict_spec_round_trips(tmp_path):
    rows = [{"token_id": i, "level": "nr90", "source": "mmseqs", "key": f"k{i}", "n_slots": 2} for i in range(20)]
    out = tmp_path / "x.parquet"
    write_table([inferred(rows)], TOK, out)
    assert pq.read_table(out).to_pylist() == rows
