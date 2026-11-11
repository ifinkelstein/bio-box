"""Token postings: one row per slot and token (card B-10; planning/tasks/TASKS.md §2; AGENTS §3.5).

`token_postings` replaces `slot_pairs` (card B-11 removes that table): neighbour questions are a
self-join of the postings on `contig_uid` with a slot distance, so pairs are never materialized.
Per slot: one row for its NR50 token and one for its NR90 token from `tape_tokens`, each where the
slot has one, with `seg_idx` 0; and (card B-21) one family row per included domain of the slot's
protein, from `l1_domains`, with `seg_idx` the domain's rank on the protein by envelope start
(ties: library precedence, then model accession; TASKS.md §2), so a protein with two domains puts
both families on its slot and `fusions_of` finds them. A model hitting twice is two rows of one
token. A non-CDS slot has no protein and so no posting.

Family postings are read from the release's family sources: the libraries `token_dict` names at the
family level, whose domain files of the release's stages are read (`biobox.l3.tokens._domain_files`),
and each domain's token is `token_dict`'s (family, library, accession). The canonical family that
`tape_tokens.family_tokens` keeps must be one of the slot's family postings, and a slot with family
postings must have one: tape tokens and domains that disagree (a domain file changed after the
tokens were made) are refused, naming the contig, as is a domain on a slot whose model has no token.

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

Bounded memory (AGENTS §4.1; card B-29): the tape tokens are read one row group at a time (at most
131,072 slots, AGENTS §4.6) beside a cursor on each tape file, both in `contig_uid` order, and each
token row group's postings are exploded with Arrow into a flat Parquet run in a private scratch
directory. The family postings are one more flat run, streamed from one DuckDB statement that ranks
the domains per protein (a window over the domain rows) and joins them to the tapes' slots unnested
in place (a hash join whose build side is the ranked domains; no sort, no list is held whole); the
canonical families of the token row groups are written beside as small runs and checked against it
with two anti-joins. One DuckDB connection (`DUCK_MEMORY_LIMIT`, one thread, spilling into that
directory) checks the inputs, then sorts the flat runs by (bucket, sort key) in one statement; an ordered Arrow
reader hands the rows to the writer, which cuts buckets and parts while streaming. The earlier single
statement (join whole list rows, unnest, union the levels from a materialized CTE, sort) ran out of
its 2 GB at rL.3's 12 M slots: its sort never spilled when the first branch of the union was small
(measured, wilkcomp02, 2026-11-10; planning/measurements/2026-11-10_postings_memory.md).
"""

from __future__ import annotations

import bisect
import shutil
import tempfile
from collections.abc import Iterator
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from biobox import catalog, io, layout, ledger
from biobox.l1.genecall import DUCK_MEMORY_LIMIT, GENE_CALL_VERSION
from biobox.l3.tokens import LEVELS, _domain_files, _scan, _sql_str, _tape_files, tape_version_of
from biobox.registry.schemas import table

EST_POSTING_ROW = 36  # zstd bytes per posting row, measured 35.74 (rL.2, pod wilkcomp02, 2026-11-10, scripts/stress/postings_bytes.py; the earlier estimate was 20), rounded up to an int; cuts parts, picks width
_LISTS = ("feature_uids", "start_bps", "end_bps", "strands")
_TOKENS = tuple(f"{lvl}_tokens" for lvl in LEVELS)
_FAMILY, _CLUSTERS = "family_tokens", tuple(c for c in _TOKENS if c != "family_tokens")
_CANON_SCHEMA = pa.schema([table("tape_tokens").schema.field("contig_uid"), pa.field("ordinal", pa.uint32()),
                           pa.field("token_id", pa.uint64())])  # a canonical family on a slot (card B-21)


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


_TAPE_COLS = ("contig_uid", "contig_len", "is_circular", "n_slots", *_LISTS)
_TAPE_SCHEMA = pa.schema([table("contig_tapes").schema.field(c) for c in _TAPE_COLS])
_TOK_COLS = ("contig_uid", *_TOKENS)
_TOK_SCHEMA = pa.schema([table("tape_tokens").schema.field(c) for c in _TOK_COLS])


class _TapeCursor:
    """One tape file, read one row group at a time and consumed from the front in `contig_uid` order."""

    def __init__(self, path: Path):
        self._pf = pq.ParquetFile(path)
        self._next, self._t, self._uids, self._off = 0, None, [], 0

    def _load(self) -> bool:
        """Whether a row remains (the next row group is read when the current one is used up)."""
        while self._t is None or self._off == self._t.num_rows:
            if self._next == self._pf.metadata.num_row_groups:
                return False
            t = self._pf.read_row_group(self._next, columns=list(_TAPE_COLS))
            self._t = t.select(list(_TAPE_COLS)).cast(_TAPE_SCHEMA)
            self._uids, self._off, self._next = self._t.column("contig_uid").to_pylist(), 0, self._next + 1
        return True

    def through(self, hi: bytes) -> list[pa.Table]:
        """The rows from the front up to and including `contig_uid` `hi`."""
        out = []
        while self._load():
            end = bisect.bisect_right(self._uids, hi, self._off)
            if end > self._off:
                out.append(self._t.slice(self._off, end - self._off))
            self._off = end
            if end < len(self._uids):
                break
        return out

    def first(self) -> bytes | None:
        return self._uids[self._off] if self._load() else None


def _explode_runs(tape_files: list[Path], tok_files: list[Path], spec, out: Path
                  ) -> tuple[list[Path], int, list[Path]]:
    """The cluster postings, unsorted, as flat Parquet runs under `out`, one per tape-token row group, and the
    canonical families of the slots (contig, ordinal, token) as runs beside them; returns the posting runs,
    their row count and the canonical runs (cards B-29, B-21).

    The tape tokens (sorted by `contig_uid` across their parts) are read one row group at a time; each
    tape file is a cursor in `contig_uid` order, so the tape rows of a token row group are the rows up to
    its last contig in every tape file, and each file is read once. A token row group holds at most
    131,072 slots (AGENTS §4.6), so only that and one row group per tape file are held at a time.
    Called after `_check`; tapes and tokens that still do not pair up are a RuntimeError.
    """
    out.mkdir(parents=True)
    cursors = [_TapeCursor(f) for f in tape_files]
    runs, n, canons = [], 0, []
    for f in tok_files:
        pf = pq.ParquetFile(f)
        for i in range(pf.metadata.num_row_groups):
            k = pf.read_row_group(i, columns=list(_TOK_COLS)).select(list(_TOK_COLS)).cast(_TOK_SCHEMA)
            if not k.num_rows:
                continue
            ku = k.column("contig_uid").to_pylist()
            pieces = [p for c in cursors for p in c.through(ku[-1])]
            t = pa.concat_tables(pieces) if pieces else _TAPE_SCHEMA.empty_table()
            tu = t.column("contig_uid").to_pylist()
            order = sorted(range(len(tu)), key=tu.__getitem__)
            if order != list(range(len(tu))):
                t, tu = t.take(pa.array(order, pa.int64())), [tu[j] for j in order]
            if tu != ku:
                at = next((a for a, b in zip(tu, ku) if a != b), max(tu, ku, key=len)[min(len(tu), len(ku))])
                raise RuntimeError(f"contig {at.hex()}: tapes and tape tokens do not pair up in {f.name}, row group {i}")
            posts, canon = _explode(t, k, spec)
            n += posts.num_rows
            if posts.num_rows:
                runs.append(out / f"run-{len(runs):06d}.parquet")
                pq.write_table(posts, runs[-1], compression="zstd", compression_level=1)
            if canon.num_rows:
                canons.append(out / f"canon-{len(canons):06d}.parquet")
                pq.write_table(canon, canons[-1], compression="zstd", compression_level=1)
    left = next((u for c in cursors if (u := c.first()) is not None), None)
    if left is not None:
        raise RuntimeError(f"contig {left.hex()}: a tape without tape tokens")
    return runs, n, canons


def _explode(t: pa.Table, k: pa.Table, spec) -> tuple[pa.Table, pa.Table]:
    """The cluster postings of tape rows `t` and their tape-token rows `k` (the same contigs in the same
    order), level by level in slot order: a row per slot and non-null token, with its 0-based `ordinal`;
    and the canonical family of each slot that has one, as (contig_uid, ordinal, token_id) (card B-21:
    family postings come from the domains, `_family_run`)."""
    col = lambda tbl, c: tbl.column(c).combine_chunks()  # noqa: E731
    lens = pc.list_value_length(col(t, "feature_uids")).to_numpy(zero_copy_only=False).astype(np.int64)
    for tbl, c in [(t, c) for c in _LISTS] + [(k, c) for c in _TOKENS]:
        if not np.array_equal(pc.list_value_length(col(tbl, c)).to_numpy(zero_copy_only=False), lens):
            raise RuntimeError(f"{c}: lists of another length than the tape's feature_uids")
    parent = np.repeat(np.arange(len(lens), dtype=np.int64), lens)
    ordinal = np.arange(len(parent), dtype=np.int64) - np.repeat(np.cumsum(lens) - lens, lens)
    flat = {c: col(t, c).flatten() for c in _LISTS}
    per_contig = {c: col(t, c) for c in ("contig_uid", "n_slots", "is_circular", "contig_len")}

    def valid(c: str) -> tuple[pa.Array, np.ndarray]:
        tok = col(k, c).flatten()
        return tok, np.flatnonzero(tok.is_valid().to_numpy(zero_copy_only=False))

    tok, keep = valid(_FAMILY)
    canon = pa.Table.from_arrays([per_contig["contig_uid"].take(pa.array(parent[keep], pa.int64())),
                                  pa.array(ordinal[keep]), tok.take(pa.array(keep, pa.int64()))],
                                 names=_CANON_SCHEMA.names).cast(_CANON_SCHEMA)
    levels = []
    for c in _CLUSTERS:
        tok, keep = valid(c)
        idx, p = pa.array(keep, pa.int64()), pa.array(parent[keep], pa.int64())
        arrays = {"token_id": tok.take(idx), "contig_uid": per_contig["contig_uid"].take(p),
                  "ordinal": pa.array(ordinal[keep]), "seg_idx": pa.array(np.zeros(len(keep), np.uint16)),
                  "feature_uid": flat["feature_uids"].take(idx), "start": flat["start_bps"].take(idx),
                  "end": flat["end_bps"].take(idx), "strand": flat["strands"].take(idx),
                  "n_slots": per_contig["n_slots"].take(p), "is_circular": per_contig["is_circular"].take(p),
                  "contig_len": per_contig["contig_len"].take(p)}
        levels.append(pa.Table.from_arrays([arrays[n] for n in spec.schema.names], names=spec.schema.names)
                      .cast(spec.schema))
    return pa.concat_tables(levels), canon


def run(root: Path, release: str, stages: list[str], gene_call_version: str = GENE_CALL_VERSION,
        producers: list[str] | None = None, prefix_width: int | None = None,
        max_file_bytes: int = layout.MAX_FILE_BYTES, tmp=None, memory_limit: str = DUCK_MEMORY_LIMIT,
        out: Path | None = None) -> dict:
    """Write the release's postings for the tape version of `producers` (default: the gene-call producer).

    Reads the tapes of that tape version for `stages`, the release's stages (card B-25; a listed stage
    without tapes raises FileNotFoundError, and they must be the stages of the release's `release_def`),
    its `tape_tokens` and `token_dict` of `release` (written by `biobox.l3.tokens.run`), and the domains of
    the release's family libraries in `stages` (card B-21). `prefix_width` overrides
    the bucket width; `tmp` holds DuckDB's private scratch directory (system temp by default; never
    the data root nor `out`); `memory_limit` is DuckDB's cap (a test passes a smaller one, scaled with
    its store).

    Inputs are read from `root`; the postings and (when the root has none) the `release_def` are written
    under `out` (default `root`; card C-08). Refuses, before any work, when the release's postings of this
    tape version already hold a file in `root` or `out`.
    """
    root = Path(root)
    out = root if out is None else Path(out)
    stages = catalog.check_stages(stages)
    if prefix_width not in (None, 0, 1, 2):
        raise ValueError(f"prefix_width {prefix_width}: must be 0, 1 or 2")
    ledger.check_scratch(tmp, root, out)
    tv = tape_version_of(gene_call_version, producers)
    tape_files = _tape_files(root, tv, stages)
    tok_dir = root / f"l3/tape_tokens/tape_version={tv}/release={release}"
    tok_files = io.part_files(tok_dir)
    if not tok_files:
        raise FileNotFoundError(f"{tok_dir}: no tape_tokens for tape version {tv} and release {release}; "
                                "run biobox.l3.tokens first")
    rel_base = f"l3/token_postings/release={release}/tape_version={tv}"
    catalog.verify_release_def(root, release, stages)  # a disagreeing definition is named before a refusal
    ledger.refuse_written(root, out, [rel_base], f"postings of release {release}")
    base = out / rel_base
    catalog.define_release(root, release, stages, out=out)

    tape_src = _scan(tape_files, "contig_uid, contig_len, is_circular, n_slots, protein_uids, " + ", ".join(_LISTS),
                     "NULL::BLOB AS contig_uid, NULL::BIGINT AS contig_len, NULL::BOOLEAN AS is_circular, "
                     "NULL::UINTEGER AS n_slots, NULL::BLOB[] AS protein_uids, NULL::BLOB[] AS feature_uids, "
                     "NULL::UINTEGER[] AS start_bps, NULL::UINTEGER[] AS end_bps, NULL::TINYINT[] AS strands")
    tok_src = _scan(tok_files, "contig_uid, " + ", ".join(_TOKENS), "")
    dict_files = io.part_files(root / f"l3/token_dict/release={release}")
    spec = table("token_postings")
    names = spec.schema.names
    scratch = Path(tempfile.mkdtemp(prefix="biobox-postings-", dir=tmp))
    con = duckdb.connect(config={"threads": 1, "memory_limit": memory_limit,
                                 "temp_directory": str(scratch / "duckdb"), "preserve_insertion_order": True})
    try:
        con.execute(f"CREATE VIEW tapes AS SELECT * FROM {tape_src}")
        con.execute(f"CREATE VIEW toks AS SELECT * FROM {tok_src}")
        _check(con)
        fam_run = scratch / "family.parquet"
        n_fam = _family_run(con, root, release, stages, dict_files, spec, fam_run)
        # rows per level: one family row per included domain on a slot; the non-null cluster tokens (every token
        # row has a tape of the same length, by _check)
        counts = con.execute("SELECT " + ", ".join(
            f"coalesce(sum(len(list_filter({c}, lambda x: x IS NOT NULL))), 0)" for c in _CLUSTERS) + " FROM toks").fetchall()[0]
        by_level = {"family": n_fam, **{c.removesuffix("_tokens"): int(n) for c, n in zip(_CLUSTERS, counts)}}
        n_rows = sum(by_level.values())
        width = (layout.prefix_width(EST_POSTING_ROW * n_rows) if prefix_width is None else prefix_width) if n_rows else 0
        bk = {0: "0", 1: "(token_id & 255) >> 4", 2: "token_id & 255"}[width]
        runs, exploded, canons = _explode_runs(tape_files, tok_files, spec, scratch / "runs")
        if exploded != n_rows - n_fam:
            raise RuntimeError(f"token_postings: exploded {exploded} rows, the tape tokens hold {n_rows - n_fam}")
        _check_canonical(con, canons, fam_run if n_fam else None)
        runs += [fam_run] if n_fam else []
        per_part = max(1, max_file_bytes // EST_POSTING_ROW)  # layout.cut on rows of equal estimated size
        # the flat runs, sorted by (bucket, sort key) in one statement: DuckDB's sort spills
        reader = con.sql(f"SELECT ({bk})::BIGINT AS bk, * FROM read_parquet([{', '.join(_sql_str(r) for r in runs)}], "
                         "hive_partitioning = false) ORDER BY bk, token_id, contig_uid, ordinal, seg_idx"
                         ).to_arrow_reader(io.ROW_GROUP_ROWS) if runs else None
        n_files, n_buckets, written = 0, 0, 0
        try:
            src = _Peek(_pieces(reader or [], names))
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
            if reader is not None:
                reader.close()
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)
    if written != n_rows:
        raise RuntimeError(f"token_postings: wrote {written} rows, the tape tokens and domains hold {n_rows}")
    return {"rows": n_rows, "rows_by_level": by_level, "buckets": n_buckets, "files": n_files, "tape_version": tv}


def _family_run(con, root: Path, release: str, stages: list[str], dict_files: list[Path], spec, path: Path) -> int:
    """Write the family postings, unsorted, to the flat run `path` and return their count (card B-21).

    The family sources are the libraries `token_dict` names at the family level (none without a
    dictionary); their domain files of `stages` are read. Each protein's distinct included domains
    (library, model, envelope start; a domain found in two stages is one) are ranked by envelope start,
    then token id (library precedence, then accession in code point order: `biobox.l3.tokens` ranks keys
    so), from 0: that is `seg_idx`. The slots come from the tapes' lists unnested side by side, so no
    list is repeated per slot. A domain on a slot whose model has no family token raises ValueError.
    """
    if not dict_files:
        return 0
    dicts = _scan(dict_files, "token_id, level, source, key", "")
    libs = [lib for (lib,) in con.execute(f"SELECT DISTINCT source FROM {dicts} WHERE level = 'family' "
                                          "ORDER BY source").fetchall()]
    parts = [f"SELECT {_sql_str(lib)} AS lib, protein_uid, model_acc, env_start FROM read_parquet("
             f"[{', '.join(_sql_str(f) for f in files)}], hive_partitioning = false)"
             for lib in libs if (files := _domain_files(root, lib, stages))]
    if not parts:
        return 0
    sql = (
        f"WITH d AS (SELECT DISTINCT protein_uid, env_start, lib, model_acc FROM ({' UNION ALL '.join(parts)})), "
        f"f AS (SELECT token_id, source, key FROM {dicts} WHERE level = 'family'), "
        "g AS (SELECT d.protein_uid, d.lib, d.model_acc, f.token_id, row_number() OVER (PARTITION BY d.protein_uid "
        "      ORDER BY d.env_start, f.token_id, d.lib, d.model_acc) - 1 AS seg_idx "
        "      FROM d LEFT JOIN f ON f.source = d.lib AND f.key = d.model_acc), "
        "s AS (SELECT contig_uid, n_slots, is_circular, contig_len, unnest(protein_uids) AS pu, "
        "      unnest(range(len(protein_uids))) AS ordinal, unnest(feature_uids) AS feature_uid, "
        "      unnest(start_bps) AS start, unnest(end_bps) AS \"end\", unnest(strands) AS strand FROM tapes) "
        "SELECT g.token_id, s.contig_uid, s.ordinal, g.seg_idx, s.feature_uid, s.start, s.\"end\", s.strand, "
        "s.n_slots, s.is_circular, s.contig_len, g.lib AS _lib, g.model_acc AS _acc FROM s JOIN g ON g.protein_uid = s.pu")
    names = spec.schema.names
    n, writer = 0, None
    reader = con.sql(sql).to_arrow_reader(io.ROW_GROUP_ROWS)
    try:
        for b in reader:
            if not b.num_rows:
                continue
            if b.column("token_id").null_count:
                j = b.column("token_id").is_null().to_numpy(zero_copy_only=False).argmax()
                r = {c: b.column(c)[j].as_py() for c in ("contig_uid", "ordinal", "_lib", "_acc")}
                raise ValueError(f"contig {bytes(r['contig_uid']).hex()}, slot {r['ordinal']}: an included domain of "
                                 f"{r['_lib']} {r['_acc']} has no family token in the token_dict of release {release} "
                                 "(the domains changed after biobox.l3.tokens ran)")
            batch = io.conform(pa.RecordBatch.from_arrays([b.column(c) for c in names], names=names), spec)
            if writer is None:
                writer = pq.ParquetWriter(path, spec.schema, compression="zstd", compression_level=1)
            writer.write_batch(batch)
            n += batch.num_rows
    finally:
        reader.close()
        if writer is not None:
            writer.close()
    return n


def _check_canonical(con, canons: list[Path], fam_run: Path | None) -> None:
    """Refuse tape tokens and family postings that disagree (card B-21): a slot whose canonical family in
    `tape_tokens` is none of its family postings, or a slot with family postings and no canonical family.
    Names the first such contig. Two anti-joins over flat runs (one row per family slot, one per posting)."""
    def src(files: list[Path]) -> str:
        if not files:
            return "(SELECT NULL::BLOB AS contig_uid, NULL::UINTEGER AS ordinal, NULL::UBIGINT AS token_id WHERE false)"
        return (f"(SELECT contig_uid, ordinal, token_id FROM read_parquet([{', '.join(_sql_str(f) for f in files)}], "
                "hive_partitioning = false))")

    canon, fam = src(canons), src([fam_run] if fam_run is not None else [])
    bad = con.execute(f"SELECT c.contig_uid, c.ordinal, c.token_id FROM {canon} c ANTI JOIN {fam} f "
                      "ON f.contig_uid = c.contig_uid AND f.ordinal = c.ordinal AND f.token_id = c.token_id "
                      "ORDER BY 1, 2 LIMIT 1").fetchall()
    if bad:
        cu, i, tid = bad[0]
        raise ValueError(f"contig {bytes(cu).hex()}, slot {i}: its canonical family token {tid} in tape_tokens is on no "
                         "included domain of the slot's protein (the domains changed after biobox.l3.tokens ran)")
    bad = con.execute(f"SELECT f.contig_uid, f.ordinal FROM {fam} f ANTI JOIN {canon} c "
                      "ON c.contig_uid = f.contig_uid AND c.ordinal = f.ordinal ORDER BY 1, 2 LIMIT 1").fetchall()
    if bad:
        cu, i = bad[0]
        raise ValueError(f"contig {bytes(cu).hex()}, slot {i}: family postings from its domains but no canonical "
                         "family in tape_tokens (the domains changed after biobox.l3.tokens ran)")


def _check(con) -> None:
    """Refuse a contig whose tape and token lists disagree: a token row without a tape, a tape without a
    token row, a contig twice, or lists of another length than the tape's slots. Names the first contig.

    The join holds list lengths, not lists (card B-29: joined whole, the token lists were a 340 MB hash
    table at rL.3's 12 M slots, measured)."""
    tape_lists = (*_LISTS, "protein_uids")  # protein_uids: the family postings unnest it beside the others (B-21)
    lens = " OR ".join([f"t.{c} <> t.n_slots" for c in tape_lists] + [f"k.{c} <> t.n_slots" for c in _TOKENS])
    bad = con.execute(
        "SELECT coalesce(t.contig_uid, k.contig_uid) AS cu, t.contig_uid IS NULL AS no_tape, "
        "k.contig_uid IS NULL AS no_tokens "
        f"FROM (SELECT contig_uid, n_slots, {', '.join(f'len({c}) AS {c}' for c in tape_lists)} FROM tapes) t "
        f"FULL OUTER JOIN (SELECT contig_uid, {', '.join(f'len({c}) AS {c}' for c in _TOKENS)} FROM toks) k "
        f"ON k.contig_uid = t.contig_uid WHERE t.contig_uid IS NULL OR k.contig_uid IS NULL OR {lens} "
        "ORDER BY cu LIMIT 1").fetchall()
    if bad:
        cu, no_tape, no_tokens = bad[0]
        what = ("has tape tokens but no tape" if no_tape else "has a tape but no tape tokens" if no_tokens
                else "has token lists of another length than its tape's slots")
        raise ValueError(f"contig {bytes(cu).hex()} {what}")
    dup = con.execute("SELECT cu FROM (SELECT contig_uid AS cu FROM tapes UNION ALL SELECT contig_uid FROM toks) "
                      "GROUP BY cu HAVING count(*) > 2 ORDER BY cu LIMIT 1").fetchall()
    if dup:
        raise ValueError(f"contig {bytes(dup[0][0]).hex()} has more than one tape or tape-token row")
