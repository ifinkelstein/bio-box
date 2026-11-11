"""Look a key up in sorted Parquet files without loading them (card B-13; AGENTS §4.4).

`KeyLookup(files, key_column, columns)` reads each file's footer once and keeps, per row group,
the min and max of `key_column` from its statistics; nothing else is loaded. `get(key)` finds the
row groups whose range holds the key by binary search (a key equal to a boundary may sit in two
adjacent row groups, and several files may hold it), reads only those row groups and only the
asked columns, and finds the key inside a row group by binary search, since every file is sorted by
the key (AGENTS §3.7). An LRU of `cache_row_groups` decoded row groups serves repeats.

Counters: `footers_read`, `row_groups_read` (a cache hit counts nothing).
"""

from __future__ import annotations

import bisect
import threading
from collections import OrderedDict
from collections.abc import Iterable
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


class _Values:
    """A read-only sequence view of an Arrow array, so `bisect` probes O(log n) values."""

    __slots__ = ("_a",)

    def __init__(self, a: pa.Array):
        self._a = a

    def __len__(self) -> int:
        return len(self._a)

    def __getitem__(self, i: int):
        return self._a[i].as_py()


class _File:
    __slots__ = ("path", "metadata", "row_groups", "mins", "maxs")

    def __init__(self, path: Path, metadata: pq.FileMetaData, row_groups: list[int], mins: list, maxs: list):
        self.path, self.metadata, self.row_groups, self.mins, self.maxs = path, metadata, row_groups, mins, maxs


class KeyLookup:
    """Rows of sorted Parquet `files` by the value of `key_column`; see the module docstring."""

    def __init__(self, files: Iterable[Path], key_column: str, columns: Iterable[str] | None = None,
                 cache_row_groups: int = 64):
        self.key_column = key_column
        self.columns = [key_column] + [c for c in dict.fromkeys(columns or []) if c != key_column]
        self.cache_row_groups = cache_row_groups
        self.footers_read = 0
        self.row_groups_read = 0
        self._cache: OrderedDict[tuple[int, int], tuple[pa.Table, _Values]] = OrderedDict()
        self._lock = threading.Lock()
        self._files: list[_File] = []
        for path in files:
            path = Path(path)
            md = pq.read_metadata(path)
            self.footers_read += 1
            leaf = [i for i in range(md.num_columns) if md.schema.column(i).path == key_column]
            if not leaf:
                raise ValueError(f"{path}: no column {key_column!r}")
            rgs, mins, maxs = [], [], []
            for g in range(md.num_row_groups):
                rg = md.row_group(g)
                if rg.num_rows == 0:
                    continue
                st = rg.column(leaf[0]).statistics
                if st is None or not st.has_min_max:
                    raise ValueError(f"{path}: row group {g} has no min/max statistics on {key_column}; "
                                     "a key lookup cannot prune it (LY-03)")
                lo, hi = st.min, st.max
                if lo > hi or (maxs and lo < maxs[-1]):
                    raise ValueError(f"{path}: row group {g} is not sorted by {key_column} (LY-02)")
                rgs.append(g)
                mins.append(lo)
                maxs.append(hi)
            self._files.append(_File(path, md, rgs, mins, maxs))

    @property
    def files(self) -> list[Path]:
        return [f.path for f in self._files]

    def _decoded(self, fi: int, g: int) -> tuple[pa.Table, _Values]:
        hit = self._cache.pop((fi, g), None)
        if hit is not None:
            self._cache[(fi, g)] = hit
            return hit
        f = self._files[fi]
        t = pq.ParquetFile(f.path, metadata=f.metadata).read_row_group(g, columns=self.columns)
        self.row_groups_read += 1
        out = (t, _Values(t.column(self.key_column).combine_chunks()))
        if self.cache_row_groups > 0:
            self._cache[(fi, g)] = out
            if len(self._cache) > self.cache_row_groups:
                self._cache.popitem(last=False)
        return out

    def get(self, key) -> list[dict]:
        """Every row whose `key_column` equals `key`, as dicts of the asked columns (file order)."""
        out: list[dict] = []
        with self._lock:
            for fi, f in enumerate(self._files):
                first = bisect.bisect_left(f.maxs, key)  # first row group whose max >= key
                last = bisect.bisect_right(f.mins, key)  # row groups before `last` have min <= key
                for j in range(first, last):
                    t, keys = self._decoded(fi, f.row_groups[j])
                    a, b = bisect.bisect_left(keys, key), bisect.bisect_right(keys, key)
                    if a < b:
                        out += t.slice(a, b - a).to_pylist()
        return out
