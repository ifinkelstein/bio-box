"""The one Parquet writer (card B-01; AGENTS §3.7, §3.9, §3.10, §4.1, §4.3, §4.6).

`write_table` streams record batches into one zstd Parquet file. Each batch is cast to the
registry schema; rows must be sorted by the spec's sort key (checked within and across batches,
nulls first, bytes unsigned, the order of `contracts._sort_ok`); row groups depend only on the rows,
never on how the input was batched, so the bytes do not either; min/max statistics are written for
the sort-key columns. The file is written under `.<name>.tmp-<pid>` beside its final path and
renamed into place only after the writer closed cleanly. At most one row group is buffered.

`duck_batches` streams a DuckDB query result as batches cast to a spec (DuckDB has no fixed-size
binary type and returns UIDs and digests as BLOB).
"""

from __future__ import annotations

import contextlib
import io
import os
from collections.abc import Iterable, Iterator
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from biobox.registry.schemas import TableSpec

ROW_GROUP_ROWS = 122_880  # AGENTS §4.6
_CAST_ERRORS = (pa.ArrowInvalid, pa.ArrowNotImplementedError, pa.ArrowTypeError)


def conform(batch: pa.RecordBatch, spec: TableSpec) -> pa.RecordBatch:
    """`batch` with the spec's columns, in its order and types. Missing or extra columns, a failed
    cast (a UID of the wrong width, an integer out of range) and a null in a non-nullable field are
    errors naming the column."""
    have = batch.schema.names
    missing = [n for n in spec.schema.names if n not in have]
    extra = [n for n in have if n not in spec.schema.names]
    if missing or extra or len(set(have)) != len(have):
        raise ValueError(f"{spec.name}: columns differ from the registry: missing {missing}, extra {extra}"
                         + ("" if len(set(have)) == len(have) else ", duplicated names"))
    cols = []
    for f in spec.schema:
        col = batch.column(f.name)
        try:
            col = col.cast(f.type)
        except _CAST_ERRORS as e:
            raise ValueError(f"{spec.name}.{f.name}: cannot cast {col.type} to {f.type}: {e}") from None
        if not f.nullable and col.null_count:
            raise ValueError(f"{spec.name}.{f.name}: {col.null_count} nulls in a non-nullable column")
        cols.append(col)
    return pa.RecordBatch.from_arrays(cols, schema=spec.schema)


def _lt_eq(cur: pa.Array, prev: pa.Array) -> tuple[pa.Array, pa.Array]:
    """Elementwise `cur < prev` and `cur == prev`, nulls first (null == null)."""
    lt, eq = pc.less(cur, prev), pc.equal(cur, prev)
    if cur.null_count or prev.null_count:
        cn, pn = pc.is_null(cur), pc.is_null(prev)
        lt = pc.or_(pc.fill_null(lt, False), pc.and_not(cn, pn))
        eq = pc.or_(pc.fill_null(eq, False), pc.and_(cn, pn))
    return lt, eq


def first_unsorted(keys: list[pa.Array]) -> int | None:
    """Index of the first row that sorts before the row above it (lexicographic over `keys`), or None."""
    n = len(keys[0]) if keys else 0
    if n < 2:
        return None
    bad = prefix_eq = None
    for a in keys:
        lt, eq = _lt_eq(a.slice(1), a.slice(0, n - 1))
        step = lt if prefix_eq is None else pc.and_(prefix_eq, lt)
        bad = step if bad is None else pc.or_(bad, step)
        prefix_eq = eq if prefix_eq is None else pc.and_(prefix_eq, eq)
    i = pc.index(bad, True).as_py()
    return None if i < 0 else i + 1


def _leaf_paths(schema: pa.Schema) -> list[str]:
    """Parquet leaf column paths (`protein_uids.list.element`), as `contracts.check_statistics` sees them."""
    buf = io.BytesIO()
    pq.write_table(schema.empty_table(), buf)
    md = pq.ParquetFile(pa.BufferReader(buf.getvalue())).metadata
    return [md.schema.column(i).path for i in range(md.num_columns)]


def _stat_paths(spec: TableSpec, names: list[str]) -> list[str]:
    leaves, out = _leaf_paths(spec.schema), []
    for name in names:
        hit = [p for p in leaves if p == name or p.startswith(name + ".")]
        if not hit:
            raise ValueError(f"{spec.name}: no column {name!r} for statistics")
        out += [p for p in hit if p not in out]
    return out


def _as_batches(batches) -> Iterator[pa.RecordBatch]:
    if isinstance(batches, (pa.Table, pa.RecordBatch)):
        batches = [batches]
    for b in batches:
        if isinstance(b, pa.Table):
            yield from b.to_batches()
        elif isinstance(b, pa.RecordBatch):
            yield b
        else:
            raise TypeError(f"expected pyarrow.RecordBatch or Table, got {type(b).__name__}")


def write_table(batches: Iterable[pa.RecordBatch | pa.Table] | pa.Table, spec: TableSpec, path, *,
                row_group_rows: int = ROW_GROUP_ROWS, max_slots: int | None = None,
                slots_column: str | None = None, statistics: Iterable[str] | None = None) -> dict:
    """Write sorted batches as one Parquet file at `path`, which must not exist yet.

    Every row group holds `row_group_rows` rows except the last. With `max_slots` and `slots_column`,
    a row group is also closed before a row would push its sum of `slots_column` past `max_slots`; a
    single row above `max_slots` is a row group of its own. Returns rows, row groups and file bytes.
    """
    path = Path(path)
    if row_group_rows < 1:
        raise ValueError("row_group_rows must be positive")
    if (max_slots is None) != (slots_column is None):
        raise ValueError("max_slots and slots_column go together")
    if slots_column is not None and not pa.types.is_integer(spec.schema.field(slots_column).type):
        raise ValueError(f"{spec.name}.{slots_column}: slot counts must be an integer column")
    if path.exists():
        raise FileExistsError(f"{path}: stages are append-only; refusing to overwrite")
    key = [c for c in spec.sort_key if c in spec.schema.names]  # partition keys live in the path
    stats = _stat_paths(spec, key + [s for s in (statistics or []) if s not in key])
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")

    pending: list[pa.RecordBatch] = []
    n_pend = s_pend = rows = groups = 0
    last: list[pa.Array] | None = None
    writer = pq.ParquetWriter(tmp, spec.schema, compression="zstd", write_statistics=stats or False)

    def flush() -> None:
        nonlocal n_pend, s_pend, groups
        if n_pend:
            # one contiguous chunk per column, so page cuts and bytes depend only on the rows
            writer.write_table(pa.Table.from_batches(pending, schema=spec.schema).combine_chunks(),
                               row_group_size=n_pend)
            pending.clear()
            n_pend = s_pend = 0
            groups += 1

    try:
        for raw in _as_batches(batches):
            b = conform(raw, spec)
            n = b.num_rows
            if not n:
                continue
            cols = [b.column(c) for c in key]
            at = None
            if last is not None:
                at = 0 if first_unsorted([pa.concat_arrays([p, c.slice(0, 1)]) for p, c in zip(last, cols)]) else None
            if at is None:
                at = first_unsorted(cols)
            if at is not None:
                raise ValueError(f"{spec.name}: rows are not sorted by {key}: row {rows + at} sorts before "
                                 "the row above it")
            last = [c.slice(n - 1, 1) for c in cols]
            if max_slots is not None:
                s = b.column(slots_column)
                if s.null_count or (pc.min(s).as_py() or 0) < 0:
                    raise ValueError(f"{spec.name}.{slots_column}: slot counts must be non-negative")
                cum = np.cumsum(s.to_numpy(zero_copy_only=False), dtype=np.int64)
            off = 0
            while off < n:
                take = min(row_group_rows - n_pend, n - off)
                if max_slots is not None:
                    base = int(cum[off - 1]) if off else 0
                    fit = int(np.searchsorted(cum, base + max_slots - s_pend, side="right")) - off
                    take = min(take, max(fit, 0) if n_pend else max(fit, 1))
                if take:
                    pending.append(b.slice(off, take))
                    n_pend += take
                    if max_slots is not None:
                        s_pend += int(cum[off + take - 1]) - base
                    off += take
                if n_pend == row_group_rows or off < n:
                    flush()
            rows += n
        flush()
        writer.close()
        if path.exists():
            raise FileExistsError(f"{path}: appeared while writing; refusing to overwrite")
        os.replace(tmp, path)
    except BaseException:
        with contextlib.suppress(Exception):
            writer.close()
        tmp.unlink(missing_ok=True)
        raise
    return {"rows": rows, "row_groups": groups, "bytes": path.stat().st_size}


def duck_batches(con, sql: str, spec: TableSpec, params=None, batch_rows: int = ROW_GROUP_ROWS
                 ) -> Iterator[pa.RecordBatch]:
    """Stream the result of `sql` on a DuckDB connection as record batches cast to `spec`.

    Nothing runs until the first batch is asked for; the result is read through an Arrow reader of
    `batch_rows` rows, never materialized.
    """
    rel = con.sql(sql, params=params)
    if rel is None:
        raise ValueError(f"{spec.name}: the statement returned no result")
    reader = rel.to_arrow_reader(batch_rows)
    try:
        for b in reader:
            yield conform(b, spec)
    finally:
        reader.close()
