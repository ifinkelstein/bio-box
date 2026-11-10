"""Token postings: one row per slot and token (card B-10; planning/tasks/TASKS.md §2; AGENTS §3.5).

`token_postings` replaces `slot_pairs` (card B-11 removes that table): neighbour questions are a
self-join of the postings on `contig_uid` with a slot distance, so pairs are never materialized.
This card writes the tokens that `tape_tokens` carries: per slot, one row for its canonical family
token, one for its NR50 token and one for its NR90 token, each where the slot has one, all with
`seg_idx` 0. A non-CDS slot has no token and so no posting. Per-domain family postings
(`seg_idx` > 0) come with a later card.

Columns: `token_id`, `contig_uid`, `ordinal` (the 0-based list index in the tape, valid in one tape
version only; `feature_uid` beside it is the stable key, AGENTS §3.3), `seg_idx`, `feature_uid`,
`start`, `end`, `strand` (from the tape lists), `n_slots`, `is_circular`, `contig_len` (from the
tape row).

Writes, through `biobox.io.write_table`:

    l3/token_postings/release=<R>/tape_version=<tv>/token_bucket=<all|h|hh>/part-<k>.parquet

The bucket is the low byte of `token_id` as hex at the width `biobox.layout.prefix_width` picks for
the estimated size, as `slot_pairs.anchor_bucket` (`biobox.l3.pairs.bucket_key`), so a lookup by
token opens one bucket. Within a bucket, rows are sorted by `(token_id, contig_uid, ordinal,
seg_idx)` and cut into parts of `max_file_bytes` estimated bytes; row groups hold 122,880 rows.
With no rows at all, one empty `token_bucket=all/part-0.parquet` is written.

Bounded memory (AGENTS §4.1): one DuckDB connection (`DUCK_MEMORY_LIMIT`, one thread, spilling into
a private scratch directory) unnests the tape lists and the token lists by position, joins them on
`contig_uid`, unions the three levels and sorts by (bucket, sort key); an ordered Arrow reader hands
the rows to the writer, which cuts buckets and parts while streaming.
"""

from __future__ import annotations

import shutil
import tempfile
from collections.abc import Iterator
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa

from biobox import catalog, io, layout
from biobox.l1.genecall import DUCK_MEMORY_LIMIT, GENE_CALL_VERSION
from biobox.l3.tokens import LEVELS, _scan, _tape_files, tape_version_of
from biobox.registry.schemas import table

EST_POSTING_ROW = 36  # zstd bytes per posting row, measured 35.74 (rL.2, pod wilkcomp02, 2026-11-10, scripts/stress/postings_bytes.py; the earlier estimate was 20), rounded up to an int; cuts parts, picks width
_LISTS = ("feature_uids", "start_bps", "end_bps", "strands")
_TOKENS = tuple(f"{lvl}_tokens" for lvl in LEVELS)


class _Peek:
    """(bucket, batch) pieces, consumed from the front; a piece can be split and its rest put back."""

    def __init__(self, it: Iterator[tuple[int, pa.RecordBatch]]):
        self._it, self._head = it, None

    def peek(self) -> tuple[int, pa.RecordBatch] | None:
        if self._head is None:
            self._head = next(self._it, None)
        return self._head

    def take(self, bucket: int, n: int) -> Iterator[pa.RecordBatch]:
        """Up to `n` rows of `bucket` from the front."""
        while n and (h := self.peek()) is not None and h[0] == bucket:
            b = h[1].slice(0, n)
            rest = h[1].slice(b.num_rows)
            self._head = (bucket, rest) if rest.num_rows else None
            n -= b.num_rows
            yield b


def _pieces(reader, names: list[str]) -> Iterator[tuple[int, pa.RecordBatch]]:
    """Each batch of the reader split where the bucket number (`bk`) changes, without `bk`."""
    for b in reader:
        if not b.num_rows:
            continue
        bk = b.column("bk").to_numpy(zero_copy_only=False)
        cuts = (np.flatnonzero(np.diff(bk)) + 1).tolist()
        cols = pa.RecordBatch.from_arrays([b.column(n) for n in names], names=names)
        for a, z in zip([0, *cuts], [*cuts, len(bk)]):
            yield int(bk[a]), cols.slice(a, z - a)


def _bucket_name(bk: int, width: int) -> str:
    """The partition value of bucket number `bk` (the low byte, or its high hex digit at width 1)."""
    return layout.bucket(bytes([bk << 4 if width == 1 else bk]), width)


def run(root: Path, release: str, stages: list[str], gene_call_version: str = GENE_CALL_VERSION,
        producers: list[str] | None = None, prefix_width: int | None = None,
        max_file_bytes: int = layout.MAX_FILE_BYTES, tmp=None, memory_limit: str = DUCK_MEMORY_LIMIT) -> dict:
    """Write the release's postings for the tape version of `producers` (default: the gene-call producer).

    Reads the tapes of that tape version for `stages`, the release's stages (card B-25; a listed stage
    without tapes raises FileNotFoundError, and they must be the stages of the release's `release_def`),
    and its `tape_tokens` of `release` (written by `biobox.l3.tokens.run`). `prefix_width` overrides
    the bucket width; `tmp` holds DuckDB's private scratch directory (system temp by default; never
    the data root); `memory_limit` is DuckDB's cap (a test passes a smaller one, scaled with its store).
    """
    root = Path(root)
    stages = catalog.check_stages(stages)
    if prefix_width not in (None, 0, 1, 2):
        raise ValueError(f"prefix_width {prefix_width}: must be 0, 1 or 2")
    if tmp is not None and Path(tmp).resolve().is_relative_to(root.resolve()):
        raise ValueError(f"{tmp}: the scratch directory may not lie under the data root")
    tv = tape_version_of(gene_call_version, producers)
    tape_files = _tape_files(root, tv, stages)
    tok_dir = root / f"l3/tape_tokens/tape_version={tv}/release={release}"
    tok_files = io.part_files(tok_dir)
    if not tok_files:
        raise FileNotFoundError(f"{tok_dir}: no tape_tokens for tape version {tv} and release {release}; "
                                "run biobox.l3.tokens first")
    base = root / f"l3/token_postings/release={release}/tape_version={tv}"
    if base.exists() and any(base.rglob("*.parquet")):
        raise FileExistsError(f"{base}: postings exist; refusing to overwrite")
    catalog.define_release(root, release, stages)

    tape_src = _scan(tape_files, "contig_uid, contig_len, is_circular, n_slots, " + ", ".join(_LISTS),
                     "NULL::BLOB AS contig_uid, NULL::BIGINT AS contig_len, NULL::BOOLEAN AS is_circular, "
                     "NULL::UINTEGER AS n_slots, NULL::BLOB[] AS feature_uids, NULL::UINTEGER[] AS start_bps, "
                     "NULL::UINTEGER[] AS end_bps, NULL::TINYINT[] AS strands")
    tok_src = _scan(tok_files, "contig_uid, " + ", ".join(_TOKENS), "")
    spec = table("token_postings")
    names = spec.schema.names
    scratch = Path(tempfile.mkdtemp(prefix="biobox-postings-", dir=tmp))
    con = duckdb.connect(config={"threads": 1, "memory_limit": memory_limit,
                                 "temp_directory": str(scratch / "duckdb"), "preserve_insertion_order": True})
    try:
        con.execute(f"CREATE VIEW tapes AS SELECT * FROM {tape_src}")
        con.execute(f"CREATE VIEW toks AS SELECT * FROM {tok_src}")
        _check(con)
        # rows per level: the non-null tokens (every token row has a tape of the same length, by _check)
        counts = con.execute("SELECT " + ", ".join(
            f"coalesce(sum(len(list_filter({c}, lambda x: x IS NOT NULL))), 0)" for c in _TOKENS) + " FROM toks").fetchall()[0]
        by_level = {lvl: int(n) for lvl, n in zip(LEVELS, counts)}
        n_rows = sum(by_level.values())
        width = (layout.prefix_width(EST_POSTING_ROW * n_rows) if prefix_width is None else prefix_width) if n_rows else 0
        bk = {0: "0", 1: "(token_id & 255) >> 4", 2: "token_id & 255"}[width]
        slots = (
            "SELECT t.contig_uid, t.n_slots, t.is_circular, t.contig_len, "
            "unnest(range(len(t.feature_uids))) AS ordinal, "
            + ", ".join(f"unnest(t.{c}) AS {c}" for c in _LISTS) + ", "
            + ", ".join(f"unnest(k.{c}) AS {c}" for c in _TOKENS) + " "
            "FROM tapes t JOIN toks k ON k.contig_uid = t.contig_uid"
        )
        level_rows = " UNION ALL ".join(
            f"SELECT {c} AS token_id, contig_uid, ordinal, 0 AS seg_idx, feature_uids AS feature_uid, "
            f'start_bps AS start, end_bps AS "end", strands AS strand, n_slots, is_circular, contig_len '
            f"FROM s WHERE {c} IS NOT NULL" for c in _TOKENS)
        sql = (f"WITH s AS ({slots}), p AS ({level_rows}) "
               f"SELECT ({bk})::BIGINT AS bk, * FROM p ORDER BY bk, token_id, contig_uid, ordinal, seg_idx")
        per_part = max(1, max_file_bytes // EST_POSTING_ROW)  # layout.cut on rows of equal estimated size
        reader = con.sql(sql).to_arrow_reader(io.ROW_GROUP_ROWS)
        n_files, n_buckets, written = 0, 0, 0
        try:
            src = _Peek(_pieces(reader, names))
            if src.peek() is None:
                written += io.write_table([], spec, base / f"token_bucket={layout.ALL}/part-0.parquet")["rows"]
                n_files = n_buckets = 1
            while (h := src.peek()) is not None:
                b, k = h[0], 0
                d = base / f"token_bucket={_bucket_name(b, width)}"
                while (h := src.peek()) is not None and h[0] == b:
                    written += io.write_table(src.take(b, per_part), spec, d / f"part-{k}.parquet")["rows"]
                    k += 1
                n_files += k
                n_buckets += 1
        finally:
            reader.close()
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)
    if written != n_rows:
        raise RuntimeError(f"token_postings: wrote {written} rows, the tape tokens hold {n_rows}")
    return {"rows": n_rows, "rows_by_level": by_level, "buckets": n_buckets, "files": n_files, "tape_version": tv}


def _check(con) -> None:
    """Refuse a contig whose tape and token lists disagree: a token row without a tape, a tape without a
    token row, a contig twice, or lists of another length than the tape's slots. Names the first contig."""
    lens = " OR ".join([f"len(t.{c}) <> t.n_slots" for c in _LISTS] + [f"len(k.{c}) <> t.n_slots" for c in _TOKENS])
    bad = con.execute(
        "SELECT coalesce(t.contig_uid, k.contig_uid) AS cu, t.contig_uid IS NULL AS no_tape, "
        "k.contig_uid IS NULL AS no_tokens FROM tapes t FULL OUTER JOIN toks k ON k.contig_uid = t.contig_uid "
        f"WHERE t.contig_uid IS NULL OR k.contig_uid IS NULL OR {lens} ORDER BY cu LIMIT 1").fetchall()
    if bad:
        cu, no_tape, no_tokens = bad[0]
        what = ("has tape tokens but no tape" if no_tape else "has a tape but no tape tokens" if no_tokens
                else "has token lists of another length than its tape's slots")
        raise ValueError(f"contig {bytes(cu).hex()} {what}")
    dup = con.execute("SELECT cu FROM (SELECT contig_uid AS cu FROM tapes UNION ALL SELECT contig_uid FROM toks) "
                      "GROUP BY cu HAVING count(*) > 2 ORDER BY cu LIMIT 1").fetchall()
    if dup:
        raise ValueError(f"contig {bytes(dup[0][0]).hex()} has more than one tape or tape-token row")
