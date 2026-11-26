"""Ring R1, card B-74: stream a stage's own `protein_seq` outside DuckDB in `domains._stage_uid_lengths` and
`domains._measured_uids` when `owned` shared-contig proteins (card C-38h) are present.

Covers:
- `_duck()` sets `DUCK_MEMORY_LIMIT = "1GB"` and `preserve_insertion_order = False`.
- Under a tight `DUCK_MEMORY_LIMIT`, `_stage_uid_lengths` and `_measured_uids` (full and `[lo, hi]`-bounded)
  succeed on a large multi-row-group `protein_seq` stage with `owned` shared contigs, without passing own
  `protein_seq` files to DuckDB.
- Exact parity against an independent Python dict/set oracle across own-only, owned-only, overlapping
  (`min(length)` preserved), multi-owner, empty-files, empty-owned-match, and `[lo, hi]` boundary cases.
- `[lo, hi]` row-group pruning on `own_files` inside `_measured_uids` when `owned` is present.
- Early generator close releases DuckDB connection and temporary scratch directory.
- End-to-end `domains.plan` + `scan_shard` + `merge` byte-for-byte (sha256) parity against `domains.run`
  with both `owned` shared-contig proteins and a copy source that itself has `owned`.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import random
import re
import tempfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pyhmmer
import pytest

from biobox import ids
from biobox import io as bio
from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.registry.schemas import table

pytestmark = pytest.mark.covers("B-74", "C-38h", "MR-08", "MR-12", "D-01", "D-03")

_SIZE_RE = re.compile(r"^\s*([0-9]+(?:\.[0-9]+)?)\s*([A-Za-z]+)\s*$")
_UNITS = {
    "B": 1,
    "KB": 10**3, "KIB": 2**10,
    "MB": 10**6, "MIB": 2**20,
    "GB": 10**9, "GIB": 2**30,
    "TB": 10**12, "TIB": 2**40,
}
_AA = "ACDEFGHIKLMNPQRSTVWY"


def _uid(i: int) -> bytes:
    return i.to_bytes(16, "big")


def _cuid(tag: int) -> bytes:
    return hashlib.blake2b(f"contig-{tag}".encode(), digest_size=16).digest()


def _fuid(contig: bytes, start: int) -> bytes:
    return hashlib.blake2b(contig + start.to_bytes(8, "big"), digest_size=16).digest()


def _write_protein_seq(
    root: Path,
    stage: str,
    rows: list[tuple[bytes, int, str]],
    *,
    prefix: str = "all",
    part: int = 0,
    row_group_rows: int = 16,
) -> Path:
    """Write a sorted `protein_seq` part for `stage`."""
    rows = sorted(rows, key=lambda r: r[0])
    spec = table("protein_seq")
    t = pa.table(
        {
            "protein_uid": pa.array([r[0] for r in rows], type=pa.binary(16)),
            "protein_digest": pa.array(
                [hashlib.sha256(r[2].encode()).digest() for r in rows], type=pa.binary(32)
            ),
            "length": pa.array([r[1] for r in rows], type=pa.uint32()),
            "seq": pa.array([r[2] for r in rows], type=pa.string()),
        },
        schema=spec.schema,
    )
    path = root / f"l0seq/protein_seq/gene_call_version={V}/stage={stage}/protein_prefix={prefix}/part-{part}.parquet"
    bio.write_table(t, spec, path, row_group_rows=row_group_rows)
    return path


def _write_feature_protein(
    root: Path,
    stage: str,
    contig_prots: list[tuple[bytes, bytes]],
    *,
    row_group_rows: int = 64,
) -> Path:
    """Write a sorted `feature_protein` part for `stage` mapping `(contig_uid, protein_uid)`."""
    spec = table("feature_protein")
    ordered = sorted(
        ((cu, i * 300, i * 300 + 270, pu) for i, (cu, pu) in enumerate(contig_prots)),
        key=lambda r: (r[0], r[1], r[2]),
    )
    t = pa.table(
        {
            "contig_uid": pa.array([r[0] for r in ordered], type=pa.binary(16)),
            "feature_uid": pa.array([_fuid(r[0], r[1]) for r in ordered], type=pa.binary(16)),
            "start": pa.array([r[1] for r in ordered], type=pa.uint32()),
            "end": pa.array([r[2] for r in ordered], type=pa.uint32()),
            "strand": pa.array([1] * len(ordered), type=pa.int8()),
            "protein_uid": pa.array([r[3] for r in ordered], type=pa.binary(16)),
            "genetic_code": pa.array([11] * len(ordered), type=pa.uint8()),
            "edit_type": pa.array(["none"] * len(ordered), type=pa.string()),
        },
        schema=spec.schema,
    )
    path = root / f"l1/feature_protein/gene_call_version={V}/stage={stage}/part-0.parquet"
    bio.write_table(t, spec, path, row_group_rows=row_group_rows)
    return path


def _write_genome_input(root: Path, stage: str, contigs: list[bytes]) -> Path:
    spec = table("genome_input")
    gid = b"\x01" * 16
    ordered = sorted(set(contigs))
    t = pa.table(
        {
            "contig_uid": pa.array(ordered, type=pa.binary(16)),
            "genome_uid": pa.array([gid] * len(ordered), type=pa.binary(16)),
            "is_circular": pa.array([False] * len(ordered), type=pa.bool_()),
            "kind": pa.array(["virus"] * len(ordered), type=pa.string()),
            "genetic_code": pa.array([11] * len(ordered), type=pa.uint8()),
            "assembly_quality": pa.array(["high"] * len(ordered), type=pa.string()),
        },
        schema=spec.schema,
    )
    path = root / f"l0seq/genome_input/stage={stage}/part-0.parquet"
    bio.write_table(t, spec, path)
    return path


def _entry(root: Path, path: Path) -> dict:
    return {"path": path.relative_to(root).as_posix(), "bytes": path.stat().st_size}


def _duck_bytes(con) -> float:
    raw = str(con.execute("SELECT current_setting('memory_limit')").fetchone()[0])
    m = _SIZE_RE.match(raw)
    assert m is not None, f"unrecognised DuckDB memory_limit format: {raw!r}"
    return float(m.group(1)) * _UNITS[m.group(2).upper()]


def _spy_duck(monkeypatch) -> list[str]:
    """Wrap `domains._duck` to record every SQL string executed via `.sql(...)` or `.execute(...)`."""
    queries: list[str] = []
    real_duck = domains._duck

    class _Proxy:
        def __init__(self, con):
            self._con = con

        def sql(self, query, *args, **kwargs):
            queries.append(str(query))
            return self._con.sql(query, *args, **kwargs)

        def execute(self, query, *args, **kwargs):
            queries.append(str(query))
            return self._con.execute(query, *args, **kwargs)

        def __getattr__(self, name):
            return getattr(self._con, name)

    @contextlib.contextmanager
    def wrapped():
        with real_duck() as con:
            yield _Proxy(con)

    monkeypatch.setattr(domains, "_duck", wrapped)
    return queries


# --- 1. DuckDB connection config -----------------------------------------------------------------------------


def test_duck_config_sets_1gb_and_disables_preserve_insertion_order():
    assert domains.DUCK_MEMORY_LIMIT == "1GB"
    with domains._duck() as con:
        pio = con.execute("SELECT current_setting('preserve_insertion_order')").fetchone()[0]
        mem = _duck_bytes(con)
    assert pio is False or str(pio).lower() == "false"
    assert 0.9 * 10**9 <= mem <= 1.1 * 2**30


# --- 2. Tight DuckDB memory limit: fails with OOM before fix, passes after -----------------------------------


def test_tight_duck_memory_limit_streams_own_protein_seq_without_oom(tmp_path, monkeypatch):
    """A stage with 300k own proteins across 74 row groups and 12 owned shared-contig proteins must not pass
    its own `protein_seq` files through DuckDB; under a tight 12MB DuckDB cap the old DuckDB union OOMs,
    whereas streaming own files via PyArrow and 2-way merging the small owned stream succeeds."""
    root = tmp_path / "root"
    n_own = 300_000
    # Build 300k own proteins split across two sorted parts using fast PyArrow arrays
    spec = table("protein_seq")
    half = n_own // 2
    own_files: list[Path] = []
    for part_idx, (lo_i, hi_i) in enumerate(((0, half), (half, n_own))):
        idxs = range(lo_i, hi_i)
        uids = [((i * 2) + 100).to_bytes(16, "big") for i in idxs]
        lens = [60 + (i % 40) for i in idxs]
        dig = b"\xab" * 32
        t = pa.table(
            {
                "protein_uid": pa.array(uids, type=pa.binary(16)),
                "protein_digest": pa.array([dig] * (hi_i - lo_i), type=pa.binary(32)),
                "length": pa.array(lens, type=pa.uint32()),
                "seq": pa.array(["M" * 60] * (hi_i - lo_i), type=pa.string()),
            },
            schema=spec.schema,
        )
        p = root / f"l0seq/protein_seq/gene_call_version={V}/stage=B/protein_prefix={part_idx:02x}/part-0.parquet"
        bio.write_table(t, spec, p, row_group_rows=4_096)
        own_files.append(p)

    # Owner stage A has 1 shared contig (12 proteins: 4 overlap with B, 8 unique) + 1 unshared contig (5 proteins)
    c_shared, c_other = _cuid(1), _cuid(2)
    owned_overlap = [(((k * 1000) * 2) + 100, 50) for k in range(4)]  # smaller length than B's >= 60
    owned_unique = [(10 + k, 77 + k) for k in range(4)] + [(n_own * 2 + 500 + k, 88 + k) for k in range(4)]
    unshared = [(n_own * 2 + 900 + k, 99) for k in range(5)]
    a_rows = [(_uid(i), n, "M" * n) for i, n in owned_overlap + owned_unique + unshared]
    a_ps = _write_protein_seq(root, "A", a_rows, row_group_rows=8)
    a_fp = _write_feature_protein(
        root,
        "A",
        [(c_shared, _uid(i)) for i, _ in owned_overlap + owned_unique]
        + [(c_other, _uid(i)) for i, _ in unshared],
    )
    b_gi = _write_genome_input(root, "B", [c_shared])
    owned = {
        "genome_input": [_entry(root, b_gi)],
        "stages": [{"stage": "A", "feature_protein": [_entry(root, a_fp)], "protein_files": [_entry(root, a_ps)]}],
    }
    src = {
        "stage": "B",
        "protein_files": [_entry(root, p) for p in own_files],
        "owned": owned,
    }

    # Build independent oracle
    oracle: dict[bytes, int] = {((i * 2) + 100).to_bytes(16, "big"): 60 + (i % 40) for i in range(n_own)}
    for i, n in owned_overlap + owned_unique:
        u = _uid(i)
        oracle[u] = min(oracle[u], n) if u in oracle else n
    want_items = sorted(oracle.items())

    monkeypatch.setattr(domains, "DUCK_MEMORY_LIMIT", "12MB")
    queries = _spy_duck(monkeypatch)

    # 1. _stage_uid_lengths must succeed under 12MB cap and never reference own_files in SQL
    got_ul = list(domains._stage_uid_lengths(root, own_files, owned))
    assert len(got_ul) == len(want_items)
    assert got_ul[:20] == want_items[:20]
    assert got_ul[-20:] == want_items[-20:]
    # Check min(length) on the 4 overlapping uids
    got_map = dict(got_ul)
    for i, n in owned_overlap:
        assert got_map[_uid(i)] == n
    assert queries, "expected DuckDB to query owned stage A"
    for q in queries:
        for own_p in own_files:
            assert str(own_p) not in q, f"own protein_seq file entered DuckDB SQL: {own_p}"

    # 2. _measured_uids (full) must succeed under 12MB cap and never reference own_files in SQL
    queries.clear()
    got_uids = list(domains._measured_uids(root, src))
    assert len(got_uids) == len(want_items)
    assert got_uids == [u for u, _ in want_items]
    assert queries
    for q in queries:
        for own_p in own_files:
            assert str(own_p) not in q, f"own protein_seq file entered DuckDB SQL: {own_p}"

    # 3. _measured_uids ([lo, hi]-bounded) must succeed under 12MB cap and never reference own_files in SQL
    queries.clear()
    lo, hi = _uid(10), _uid(5_000)
    got_win = list(domains._measured_uids(root, src, lo, hi))
    assert got_win == [u for u, _ in want_items if lo <= u <= hi]
    assert queries
    for q in queries:
        for own_p in own_files:
            assert str(own_p) not in q, f"own protein_seq file entered DuckDB SQL: {own_p}"


# --- 3. Exact parity & edge cases ----------------------------------------------------------------------------


def test_stage_uid_lengths_and_measured_uids_parity_and_edge_cases(tmp_path):
    rng = random.Random(74)
    root = tmp_path / "root"

    c_sh1, c_sh2, c_un1, c_un2, c_miss = (_cuid(k) for k in range(10, 15))
    # Construct overlapping and disjoint uids across B (2 own files), A1, and A2
    pool = sorted(rng.sample(range(100, 10_000), 160))
    b_only = pool[0:50]
    a1_only = pool[50:80]
    a2_only = pool[80:105]
    b_and_a1_b_smaller = pool[105:115]
    b_and_a1_a1_smaller = pool[115:125]
    a1_and_a2_overlap = pool[125:135]
    all_three = pool[135:145]
    a1_unshared = pool[145:153]
    a2_unshared = pool[153:160]

    b_rows_raw: dict[bytes, int] = {}
    for i in b_only:
        b_rows_raw[_uid(i)] = rng.randint(50, 200)
    for i in b_and_a1_b_smaller:
        b_rows_raw[_uid(i)] = 60
    for i in b_and_a1_a1_smaller:
        b_rows_raw[_uid(i)] = 180
    for i in all_three:
        b_rows_raw[_uid(i)] = 150

    b_sorted = sorted(b_rows_raw.items())
    mid = len(b_sorted) // 2
    b_f0 = _write_protein_seq(
        root, "B", [(u, n, "M" * n) for u, n in b_sorted[:mid]], prefix="00", row_group_rows=7
    )
    b_f1 = _write_protein_seq(
        root, "B", [(u, n, "M" * n) for u, n in b_sorted[mid:]], prefix="80", row_group_rows=7
    )
    b_files = [b_f0, b_f1]

    a1_rows_raw: dict[bytes, int] = {}
    for i in a1_only:
        a1_rows_raw[_uid(i)] = rng.randint(50, 200)
    for i in b_and_a1_b_smaller:
        a1_rows_raw[_uid(i)] = 120  # larger than B's 60 -> B wins
    for i in b_and_a1_a1_smaller:
        a1_rows_raw[_uid(i)] = 90  # smaller than B's 180 -> A1 wins
    for i in a1_and_a2_overlap:
        a1_rows_raw[_uid(i)] = 110
    for i in all_three:
        a1_rows_raw[_uid(i)] = 130
    for i in a1_unshared:
        a1_rows_raw[_uid(i)] = 75
    a1_ps = _write_protein_seq(
        root, "A1", [(u, n, "M" * n) for u, n in a1_rows_raw.items()], row_group_rows=9
    )
    a1_shared_uids = [_uid(i) for i in a1_only + b_and_a1_b_smaller + b_and_a1_a1_smaller + a1_and_a2_overlap + all_three]
    a1_fp = _write_feature_protein(
        root, "A1", [(c_sh1, u) for u in a1_shared_uids] + [(c_un1, _uid(i)) for i in a1_unshared]
    )

    a2_rows_raw: dict[bytes, int] = {}
    for i in a2_only:
        a2_rows_raw[_uid(i)] = rng.randint(50, 200)
    for i in a1_and_a2_overlap:
        a2_rows_raw[_uid(i)] = 85  # smaller than A1's 110 -> A2 wins
    for i in all_three:
        a2_rows_raw[_uid(i)] = 95  # smallest of (150, 130, 95) -> A2 wins
    for i in a2_unshared:
        a2_rows_raw[_uid(i)] = 70
    a2_ps = _write_protein_seq(
        root, "A2", [(u, n, "M" * n) for u, n in a2_rows_raw.items()], row_group_rows=9
    )
    a2_shared_uids = [_uid(i) for i in a2_only + a1_and_a2_overlap + all_three]
    a2_fp = _write_feature_protein(
        root, "A2", [(c_sh2, u) for u in a2_shared_uids] + [(c_un2, _uid(i)) for i in a2_unshared]
    )

    b_gi = _write_genome_input(root, "B", [c_sh1, c_sh2])
    owned = {
        "genome_input": [_entry(root, b_gi)],
        "stages": [
            {"stage": "A1", "feature_protein": [_entry(root, a1_fp)], "protein_files": [_entry(root, a1_ps)]},
            {"stage": "A2", "feature_protein": [_entry(root, a2_fp)], "protein_files": [_entry(root, a2_ps)]},
        ],
    }

    # Independent Python oracle
    oracle: dict[bytes, int] = dict(b_rows_raw)
    for u in a1_shared_uids:
        n = a1_rows_raw[u]
        oracle[u] = min(oracle[u], n) if u in oracle else n
    for u in a2_shared_uids:
        n = a2_rows_raw[u]
        oracle[u] = min(oracle[u], n) if u in oracle else n
    want = sorted(oracle.items())
    want_uids = [u for u, _ in want]

    # Full _stage_uid_lengths parity + strict ordering + min(length)
    got_ul = list(domains._stage_uid_lengths(root, b_files, owned))
    assert got_ul == want
    assert all(got_ul[k][0] < got_ul[k + 1][0] for k in range(len(got_ul) - 1))
    assert all(isinstance(u, bytes) and isinstance(n, int) for u, n in got_ul)
    for i in b_and_a1_b_smaller:
        assert oracle[_uid(i)] == 60
    for i in b_and_a1_a1_smaller:
        assert oracle[_uid(i)] == 90
    for i in a1_and_a2_overlap:
        assert oracle[_uid(i)] == 85
    for i in all_three:
        assert oracle[_uid(i)] == 95
    for i in a1_unshared + a2_unshared:
        assert _uid(i) not in oracle

    # _measured_uids full and bounded windows
    src = {"stage": "B", "protein_files": [_entry(root, p) for p in b_files], "owned": owned}
    assert list(domains._measured_uids(root, src)) == want_uids

    windows = [
        (want_uids[0], want_uids[0]),  # single exact lowest uid
        (want_uids[-1], want_uids[-1]),  # single exact highest uid
        (want_uids[10], want_uids[40]),  # middle span across both files and owners
        (_uid(0), _uid(50)),  # empty window below all uids
        (_uid(20_000), _uid(30_000)),  # empty window above all uids
        (_uid(a1_only[0]), _uid(a1_only[0])),  # owned-only exact uid
        (_uid(b_only[0]), _uid(b_only[0])),  # own-only exact uid
    ]
    for lo, hi in windows:
        got_w = list(domains._measured_uids(root, src, lo, hi))
        want_w = [u for u in want_uids if lo <= u <= hi]
        assert got_w == want_w, (lo, hi)

    # Edge case: empty own `files` with non-empty `owned`
    owned_only_oracle: dict[bytes, int] = {}
    for u in a1_shared_uids:
        owned_only_oracle[u] = a1_rows_raw[u]
    for u in a2_shared_uids:
        owned_only_oracle[u] = min(owned_only_oracle[u], a2_rows_raw[u]) if u in owned_only_oracle else a2_rows_raw[u]
    want_owned_only = sorted(owned_only_oracle.items())
    assert list(domains._stage_uid_lengths(root, [], owned)) == want_owned_only
    src_empty_own = {"stage": "B0", "protein_files": [], "owned": owned}
    assert list(domains._measured_uids(root, src_empty_own)) == [u for u, _ in want_owned_only]
    lo_o, hi_o = want_owned_only[2][0], want_owned_only[10][0]
    assert list(domains._measured_uids(root, src_empty_own, lo_o, hi_o)) == [
        u for u, _ in want_owned_only if lo_o <= u <= hi_o
    ]

    # Edge case: non-empty own `files` with `owned` matching 0 rows (gi contig not in A1/A2)
    b_gi_miss = _write_genome_input(root, "B_miss", [c_miss])
    owned_zero = {
        "genome_input": [_entry(root, b_gi_miss)],
        "stages": owned["stages"],
    }
    assert list(domains._stage_uid_lengths(root, b_files, owned_zero)) == b_sorted
    src_zero = {"stage": "B_miss", "protein_files": [_entry(root, p) for p in b_files], "owned": owned_zero}
    assert list(domains._measured_uids(root, src_zero)) == [u for u, _ in b_sorted]

    # Edge case: `owned` with empty `stages` list
    owned_empty_stages = {"genome_input": [_entry(root, b_gi)], "stages": []}
    assert list(domains._stage_uid_lengths(root, b_files, owned_empty_stages)) == b_sorted
    src_empty_stages = {
        "stage": "B_es",
        "protein_files": [_entry(root, p) for p in b_files],
        "owned": owned_empty_stages,
    }
    assert list(domains._measured_uids(root, src_empty_stages)) == [u for u, _ in b_sorted]


# --- 4. Row-group pruning on own_files in _measured_uids when owned is present -------------------------------


def test_measured_uids_prunes_own_files_row_groups_with_owned(tmp_path, monkeypatch):
    root = tmp_path / "root"
    n_own = 200
    own_rows = [(_uid(1_000 + i * 10), 80, "M" * 80) for i in range(n_own)]
    own_p = _write_protein_seq(root, "B", own_rows, row_group_rows=10)
    assert pq.ParquetFile(own_p).metadata.num_row_groups == 20

    c_sh = _cuid(77)
    a_rows = [(_uid(1_000 + 73 * 10), 80, "M" * 80), (_uid(1_000 + 75 * 10 + 5), 90, "M" * 90)]
    a_ps = _write_protein_seq(root, "A", a_rows, row_group_rows=10)
    a_fp = _write_feature_protein(root, "A", [(c_sh, r[0]) for r in a_rows])
    b_gi = _write_genome_input(root, "B", [c_sh])
    src = {
        "stage": "B",
        "protein_files": [_entry(root, own_p)],
        "owned": {
            "genome_input": [_entry(root, b_gi)],
            "stages": [{"stage": "A", "feature_protein": [_entry(root, a_fp)], "protein_files": [_entry(root, a_ps)]}],
        },
    }

    # Spy on ParquetFile.read_row_group and DuckDB SQL
    rg_reads: list[int] = []
    real_read_rg = pq.ParquetFile.read_row_group

    def spy_read_rg(self, i, *args, **kwargs):
        rg_reads.append(i)
        return real_read_rg(self, i, *args, **kwargs)

    monkeypatch.setattr(pq.ParquetFile, "read_row_group", spy_read_rg)
    queries = _spy_duck(monkeypatch)

    # Window strictly inside row group 7 (indices 70..79 -> uids 1700..1790)
    lo, hi = _uid(1_720), _uid(1_760)
    got = list(domains._measured_uids(root, src, lo, hi))
    want = sorted({r[0] for r in own_rows + a_rows if lo <= r[0] <= hi})
    assert got == want
    assert rg_reads == [7], f"expected exactly row group 7 read from own_files, got {rg_reads}"
    assert queries and all(str(own_p) not in q for q in queries)


# --- 5. Early generator close cleans up DuckDB scratch directory ---------------------------------------------


def test_early_close_releases_duckdb_and_scratch(tmp_path, monkeypatch):
    root = tmp_path / "root"
    scratch_base = tmp_path / "duck_tmp"
    scratch_base.mkdir()

    b_p = _write_protein_seq(root, "B", [(_uid(i), 50, "M" * 50) for i in range(10, 40)], row_group_rows=5)
    c_sh = _cuid(99)
    a_rows = [(_uid(i), 50, "M" * 50) for i in range(1, 20)]
    a_ps = _write_protein_seq(root, "A", a_rows, row_group_rows=5)
    a_fp = _write_feature_protein(root, "A", [(c_sh, r[0]) for r in a_rows])
    b_gi = _write_genome_input(root, "B", [c_sh])
    owned = {
        "genome_input": [_entry(root, b_gi)],
        "stages": [{"stage": "A", "feature_protein": [_entry(root, a_fp)], "protein_files": [_entry(root, a_ps)]}],
    }
    src = {"stage": "B", "protein_files": [_entry(root, b_p)], "owned": owned}

    real_mkdtemp = tempfile.mkdtemp
    monkeypatch.setattr(tempfile, "mkdtemp", lambda prefix="tmp", dir=None: real_mkdtemp(prefix=prefix, dir=scratch_base))

    g1 = domains._stage_uid_lengths(root, [b_p], owned)
    assert next(g1) == (_uid(1), 50)
    g1.close()
    assert list(scratch_base.iterdir()) == []

    g2 = domains._measured_uids(root, src)
    assert next(g2) == _uid(1)
    g2.close()
    assert list(scratch_base.iterdir()) == []


# --- 6. End-to-end plan + scan_shard + merge byte-identity with owned + copy sources -------------------------


def _random_protein(rng: random.Random, length: int) -> str:
    return "M" + "".join(rng.choices(_AA, k=length - 1))


def _build_synthetic_hmm(prots: dict[bytes, str], path: Path) -> None:
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet, seed=42)
    background = pyhmmer.plan7.Background(alphabet)
    with open(path, "wb") as fh:
        for i, uid in enumerate(sorted(prots)):
            seq = pyhmmer.easel.TextSequence(name=f"p{i}".encode(), sequence=prots[uid]).digitize(alphabet)
            hmm, _, _ = builder.build(seq, background)
            hmm.name = f"mod{i:02d}".encode()
            hmm.accession = f"B74M{i:04d}".encode()
            hmm.cutoffs.gathering = (25.0, 25.0)
            hmm.write(fh)


def _write_run_record(root: Path, stage: str, step: str, rel_dirs: list[str]) -> None:
    outputs = []
    for d in rel_dirs:
        for p in sorted((root / d).glob("*.parquet")):
            outputs.append({
                "path": p.relative_to(root).as_posix(),
                "bytes": p.stat().st_size,
                "sha256": ids.sha256_file(p).hex(),
            })
    rec = root / f"runs/stage={stage}/{step}.json"
    rec.parent.mkdir(parents=True, exist_ok=True)
    rec.write_text(json.dumps({"stage": stage, "step": step, "outputs": outputs}, sort_keys=True))


def _dir_hashes(base: Path) -> dict[str, str]:
    return {
        p.relative_to(base).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(base.rglob("*"))
        if p.is_file()
    }


def test_e2e_plan_shard_merge_byte_identity_with_owned_and_copy_sources(tmp_path):
    """End-to-end byte-identity parity for:
    - Stage B: own proteins + shared-contig proteins owned by A (`_stage_uid_lengths`, `_shard_proteins`, `merge`)
    - Stage C: own proteins + shared-contig proteins owned by A, copying from B (exercises `_measured_uids` on
      copy source B which itself has `owned`, both full in `plan`/`merge` and `[lo, hi]`-bounded in `scan_shard`).
    Both B and C must produce byte-for-byte identical domain and ledger Parquet files to single-process `domains.run`.
    """
    rng = random.Random(2026)
    root = tmp_path / "root"
    scratch = tmp_path / "scratch"
    scratch.mkdir()

    # Construct 16 distinct synthetic proteins (120..180 aa)
    seqs = {i: _random_protein(rng, 120 + i * 4) for i in range(16)}
    uids = {i: ids.protein_uid(seqs[i]) for i in range(16)}

    # HMM library built from 6 of the proteins spanning A-shared, B-own, and C-own
    hmm = tmp_path / "b74.hmm"
    _build_synthetic_hmm({uids[i]: seqs[i] for i in (0, 2, 6, 8, 12, 14)}, hmm)

    # Stage A (called): shared contig has proteins 0..5, unshared contig has 15
    c_shared, c_a_only = _cuid(101), _cuid(102)
    a_idx = [0, 1, 2, 3, 4, 5, 15]
    _write_protein_seq(root, "A", [(uids[i], len(seqs[i]), seqs[i]) for i in a_idx], row_group_rows=3)
    _write_feature_protein(
        root, "A", [(c_shared, uids[i]) for i in range(6)] + [(c_a_only, uids[15])]
    )

    # Stage B (annotated): own proteins 4, 5 (overlap with c_shared) + 6, 7, 8, 9, 10, 11; genome_input has c_shared
    b_own_idx = [4, 5, 6, 7, 8, 9, 10, 11]
    _write_protein_seq(root, "B", [(uids[i], len(seqs[i]), seqs[i]) for i in b_own_idx], row_group_rows=3)
    _write_genome_input(root, "B", [c_shared])

    # Stage C (annotated): own proteins 8, 9 (overlap with B) + 12, 13, 14; genome_input has c_shared
    c_own_idx = [8, 9, 12, 13, 14]
    _write_protein_seq(root, "C", [(uids[i], len(seqs[i]), seqs[i]) for i in c_own_idx], row_group_rows=2)
    _write_genome_input(root, "C", [c_shared])

    # Plan, scan (small shard_residues so multiple shards run), and merge B
    pb = domains.plan(root, "B", V, "test", "1", hmm, "gathering", shard_residues=400, seq_stages=["B", "A"])
    assert pb["proteins"] == 12  # 0..11
    assert pb["shards"] >= 3
    for k in range(pb["shards"]):
        domains.scan_shard(root, "B", "test", "1", hmm, k, cpus=1, tmp=scratch)
    mb = domains.merge(root, "B", "test", "1", tmp=scratch)
    assert mb["proteins"] == 12 and mb["domains"] > 0
    _write_run_record(root, "B", "domains-test-1", ["l1/domains/library=test-1/stage=B", "ledger/stage=B/step=domains-test-1"])

    # Plan, scan, and merge C with B as copy source (B has `owned`) and A as owner
    pc = domains.plan(root, "C", V, "test", "1", hmm, "gathering", shard_residues=250, seq_stages=["C", "B", "A"])
    assert pc["proteins"] == 11  # c_shared (0..5) + B-overlap (8, 9) + new (12, 13, 14)
    assert pc["copied_proteins"] == 8  # 0..5 + 8, 9 all measured by B
    assert pc["scanned_proteins"] == 3  # 12, 13, 14 only
    assert pc["sources"] == ["B"]
    assert pc["shards"] >= 2
    for k in range(pc["shards"]):
        domains.scan_shard(root, "C", "test", "1", hmm, k, cpus=1, tmp=scratch)
    mc = domains.merge(root, "C", "test", "1", tmp=scratch)
    assert mc["proteins"] == 11 and mc["copied_proteins"] == 8 and mc["domains"] > 0

    # Independent single-process reference builds using `domains.run` on the exact union protein sets
    ref_root = tmp_path / "ref_root"
    _write_protein_seq(ref_root, "B", [(uids[i], len(seqs[i]), seqs[i]) for i in range(12)], row_group_rows=4)
    _write_protein_seq(
        ref_root, "C", [(uids[i], len(seqs[i]), seqs[i]) for i in [0, 1, 2, 3, 4, 5, 8, 9, 12, 13, 14]], row_group_rows=4
    )
    domains.run(ref_root, "B", V, "test", "1", hmm, "gathering", cpus=1)
    domains.run(ref_root, "C", V, "test", "1", hmm, "gathering", cpus=1)

    for st in ("B", "C"):
        assert _dir_hashes(root / f"l1/domains/library=test-1/stage={st}") == _dir_hashes(
            ref_root / f"l1/domains/library=test-1/stage={st}"
        )
        assert _dir_hashes(root / f"ledger/stage={st}/step=domains-test-1") == _dir_hashes(
            ref_root / f"ledger/stage={st}/step=domains-test-1"
        )
