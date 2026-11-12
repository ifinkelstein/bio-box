"""Release-time tokens for the tapes (PLAN.md §4.1 token namespace; data_layer_and_tests.md §1 L3; card B-09).

Tapes carry stable keys (`protein_uid`); a release maps each slot to tokens at three levels:

    family  the model of each included domain of the protein (card B-21); `tape_tokens.family_tokens`
            holds one per slot, the canonical domain: the first library in `libraries` order with an
            included domain, and in it the domain with the highest bitscore (ties: model
            accession, then envelope start); `tape_tokens.domain_family_tokens` (card B-39) holds per
            slot a list with the token of every distinct included domain (library, model, envelope
            start) in the order of `biobox.l3.postings`' `seg_idx` (envelope start, token id, library,
            accession), so element s is the token of the slot's family posting with seg_idx s; a
            model hitting a protein twice is in it twice, and a slot without a domain has []
    nr50    the protein's NR50 cluster in this release
    nr90    the protein's NR90 cluster in this release

`token_id` is uint64: level (3 bits) | source (5 bits) | local (56 bits). PLAN.md §4.1 gives
32 bits with 24 local and widens to 64 when counts require it; NR90 clusters pass 2^24 by
S5, so the release line starts wide. Source is the library's 1-based position in
`libraries` for family tokens and 0 for clusters; local is the rank of the key within
(level, source) among the keys that occur on at least one slot, in ascending string order of the
key (the model accession for families, the cluster uid as lowercase hex for clusters), so ids are
dense and reproducible but release-local (never a join key outside the release). A family key
occurs on a slot when the slot's protein has an included domain of the model, canonical or not
(card B-21: `biobox.l3.postings` writes one family posting per included domain, and
`fusions_of` finds two families on one slot), so new family keys never move a cluster token.

Writes, under `<root>/l3/`, through `biobox.io.write_table`:

    token_dict/release=<R>/part-<k>.parquet                     token_id -> level, source, key, n_slots
    tape_tokens/tape_version=<tv>/release=<R>/part-<k>.parquet  per contig: one list per level, and the
                                                                list of every slot's domain families

`n_slots` is the number of slots that carry the token: for a cluster, the slots of its members;
for a family, the slots whose protein has at least one included domain of the model.

`token_dict` is sorted by `token_id`, which readers rely on (the reverse lookup token -> key prunes on
it); within each (level, source) the rows are then in key order, so its row groups carry min/max
statistics on `key` as well, and `catalog.TokenKeyLookup` finds a key by footer statistics (card B-34).

and the release's `release_def` (`biobox.catalog.define_release`) if it has none. Only the release's
stages are read (card B-25): tapes and domains of a stage the release does not list are ignored.
The tape version is that of `producers` (default: the gene-call producer alone, card B-07); a
non-CDS slot has no protein and so no token (a null in each list), and a tape without slots has
four empty lists.

Bounded memory (AGENTS §4.1; card B-46): one DuckDB connection (`memory_limit`, default
`DUCK_MEMORY_LIMIT`, one thread, spilling into a private scratch directory) and Python that holds a
batch at a time. DuckDB's hash join does not stay inside its limit when the build side is large: one
LEFT JOIN of rM.1's 48.4 M slots to its 45.8 M NR50 cluster rows ran out of 2 GB alone (measured,
wilkcomp02, 2026-11-12, DuckDB 2.0.0; card B-46's report), while its sorts and hash
aggregates spill. So no join here has a build side that grows with the release:

1. The slots (the tapes' `protein_uids` unnested with their ordinals), the cluster rows and the domain
   rows are split into `k` protein ranges by the hash of the protein uid (`_route`; slots without a
   protein in a range of their own), with `k` a power of two chosen from the input sizes and the
   memory limit (`n_ranges`, `RANGE_ROW_BYTES`).
2. Per range, each slot's canonical family and cluster keys are joined from that range's proteins, and
   the slots of every model with an included domain on them are counted.
3. The family keys (one per model) are ranked in SQL; each cluster level's keys are numbered in one
   pass over its slots sorted by key (`_number`: a key's local id is its rank, its `n_slots` its run),
   which also writes each protein's token into its range.
4. Per range, each slot gets its tokens and its protein's list of domain families.
5. The slots of all ranges, sorted back into (contig, file, row, ordinal) order, are gathered into one
   row of lists per tape in Python (`_tapes_of`) and written through `io.write_parts`; the dictionary is
   written from the family keys and the two levels' numbered runs, already in token order.

A row group of `tape_tokens` holds whole contigs and at most `max_slots` slots (AGENTS §4.6). The
output does not depend on `k` or on the hash: every file is sorted by a unique key.
"""

from __future__ import annotations

import heapq
import re
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
from biobox.l3.tapes import MAX_SLOTS, tape_version
from biobox.registry.schemas import table

LEVELS = {"family": 1, "nr50": 2, "nr90": 3}  # 4 is reserved for subfamilies
CLUSTER_SOURCE = "mmseqs2"
EST_DICT_ROW = 60  # zstd estimate per token_dict row; cuts parts only
EST_TOKEN_ROW, EST_TOKEN_SLOT = 100, 30  # per tape_tokens row and per mean slot; cuts parts only
# card B-46: DuckDB bytes allowed per input row of one protein range (slots, cluster rows or domain rows, whichever
# is largest), so 2 GB gives ranges of at most 2 M rows (derived). An estimate with a margin: one LEFT JOIN of
# rM.1's 48.4 M slots to its 45.8 M NR50 rows did not fit in 2 GB, the three joins and aggregates over rN.3's 10.6 M
# proteins and 21.2 M cluster rows did (measured, wilkcomp02, 2026-11-12, DuckDB 2.0.0)
RANGE_ROW_BYTES = 1_000
ROUTE_HELD_ROWS = 1 << 20  # rows the router holds over all its files before writing; bounds Python memory
ROUTE_OPEN_FILES = 512  # files the router writes at once; more ranges are routed in several passes


def token_id(level: str, source: int, local: int) -> int:
    if not 0 <= source < 32 or not 0 <= local < 2**56:
        raise ValueError(f"token ({level}, {source}, {local}) out of range")
    return LEVELS[level] << 61 | source << 56 | local


def token_level(tid: int) -> str:
    return {v: k for k, v in LEVELS.items()}[tid >> 61]


def tape_version_of(gene_call_version: str = GENE_CALL_VERSION, producers: list[str] | None = None) -> str:
    """The tape version of `producers`, by default the gene-call producer alone."""
    return tape_version([gene_call_version] if producers is None else producers)


def _tape_files(root: Path, tv: str, stages: list[str]) -> list[Path]:
    """The tape files of `tv` for `stages`, stage by stage in the order given and each stage's parts in
    part order (card B-25). A listed stage without tapes of `tv` raises FileNotFoundError."""
    files = []
    for s in catalog.check_stages(stages):
        found = io.part_files(Path(root) / f"l3/tapes/tape_version={tv}/stage={s}")
        if not found:
            raise FileNotFoundError(f"stage {s}: no tapes of tape version {tv} under {root}")
        files += found
    return files


def read_tapes(root: Path, gene_call_version: str = GENE_CALL_VERSION, columns: list[str] | None = None,
               producers: list[str] | None = None, *, stages: list[str]) -> list[dict]:
    """The tapes of the tape version for `stages` as a list (a test and `pairs` helper; builders use
    `iter_tapes`)."""
    tv = tape_version_of(gene_call_version, producers)
    return [r for f in _tape_files(root, tv, stages) for r in pq.read_table(f, columns=columns).to_pylist()]


def iter_tapes(root: Path, tape_version: str, stages: list[str], columns: list[str] | None = None) -> Iterator[dict]:
    """The tapes of `tape_version` one row at a time in `contig_uid` order, over the files of `stages`.

    Each file is read one row group at a time and the files are merged on `contig_uid` (equal uids
    in the order of `stages`, then part order), so one row group per file is held.
    """
    def rows(f: Path) -> Iterator[dict]:
        pf = pq.ParquetFile(f)
        cols = None if columns is None else list(dict.fromkeys(["contig_uid", *columns]))
        for i in range(pf.metadata.num_row_groups):
            yield from pf.read_row_group(i, columns=cols).to_pylist()

    keep = None if columns is None else set(columns)
    for r in heapq.merge(*(rows(f) for f in _tape_files(root, tape_version, stages)), key=lambda r: r["contig_uid"]):
        yield r if keep is None else {k: v for k, v in r.items() if k in keep}


def _domain_files(root: Path, library: str, stages: list[str]) -> list[Path]:
    """The domain files of `library` for `stages`, in the order of `stages`. A library need not be
    scanned on every stage (a stage without its files adds none)."""
    return [f for s in stages for f in io.part_files(Path(root) / f"l1/domains/library={library}/stage={s}")]


def libraries_in(root: Path, stages: list[str]) -> list[str]:
    """The libraries with domain files in at least one of `stages`, by name (card B-25: a library that
    only an unlisted stage has does not enter the release, so it cannot renumber the family sources)."""
    base = Path(root) / "l1/domains"
    names = sorted(p.name.removeprefix("library=") for p in base.iterdir() if p.name.startswith("library=")) \
        if base.exists() else []
    return [n for n in names if _domain_files(root, n, catalog.check_stages(stages))]


def _sql_str(s: str) -> str:
    return "'" + str(s).replace("'", "''") + "'"


def _scan(files: list[Path], cols: str, empty: str, numbered: bool = False) -> str:
    """A subquery over `files` (with `fidx`, the file's position in `files`, and `rn`, the row in the
    file, when `numbered`), or a typed empty one when there are none (read_parquet needs a file)."""
    if not files:
        return f"(SELECT {empty} WHERE false)"
    if not numbered:
        paths = "[" + ", ".join(_sql_str(f) for f in files) + "]"
        return f"(SELECT {cols} FROM read_parquet({paths}, hive_partitioning = false))"
    return "(" + " UNION ALL ".join(
        f"SELECT {i}::BIGINT AS fidx, file_row_number AS rn, {cols} "
        f"FROM read_parquet({_sql_str(f)}, hive_partitioning = false, file_row_number = true)"
        for i, f in enumerate(files)) + ")"


_UNITS = {"b": 1, "kb": 10**3, "mb": 10**6, "gb": 10**9, "tb": 10**12,
          "kib": 2**10, "mib": 2**20, "gib": 2**30, "tib": 2**40}


def limit_bytes(limit: str) -> int:
    """DuckDB's `memory_limit` string ("2GB", "115MB", "1.5GiB"; KB is 1000 bytes, KiB 1024, as DuckDB reads
    them) in bytes."""
    m = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*([kmgt]?i?b)\s*", str(limit), re.IGNORECASE)
    if not m or m[2].lower() not in _UNITS:
        raise ValueError(f"memory limit {limit!r}: give a number and a unit (B, KB, MB, GB, TB, KiB, MiB, GiB, TiB)")
    return int(float(m[1]) * _UNITS[m[2].lower()])


def n_ranges(n_rows: int, memory_limit: str) -> int:
    """The number of protein ranges for inputs of at most `n_rows` rows each (slots, cluster rows, domain rows):
    the least power of two that leaves at most `memory_limit` / `RANGE_ROW_BYTES` rows of each in a range."""
    per = max(1, limit_bytes(memory_limit) // RANGE_ROW_BYTES)
    k = 1
    while k * per < n_rows:
        k *= 2
    return k


def _route(con, sql: str, paths: list[Path]) -> None:
    """Stream `sql` into the scratch files `paths` by its column `_r` (a row goes to `paths[_r]`, without `_r`).
    Every file is written, an empty one too, and rows keep their order within a file. At most `ROUTE_OPEN_FILES`
    files are open at once (the pod's soft limit is 4,096 descriptors): with more paths, `sql` runs once per group
    of that many. At most about `ROUTE_HELD_ROWS` rows are held besides the batch being split."""
    for lo in range(0, len(paths), ROUTE_OPEN_FILES):
        group = paths[lo:lo + ROUTE_OPEN_FILES]
        reader = con.sql(sql).to_arrow_reader(io.ROW_GROUP_ROWS)
        schema = pa.schema([f for f in reader.schema if f.name != "_r"])
        writers, held, n_held = [], [[] for _ in group], [0] * len(group)
        per = max(1024, ROUTE_HELD_ROWS // len(group))

        def flush(i: int) -> None:
            if held[i]:
                writers[i].write_table(pa.Table.from_batches(held[i], schema=schema))
                held[i].clear()
                n_held[i] = 0

        try:
            writers += [pq.ParquetWriter(p, schema, compression="zstd", compression_level=1) for p in group]
            for b in reader:
                if not b.num_rows:
                    continue
                r = b.column("_r").to_numpy(zero_copy_only=False).astype(np.int64) - lo
                order = np.argsort(r, kind="stable")
                cut = np.searchsorted(r[order], np.arange(len(group) + 1))
                if cut[0] == cut[-1]:
                    continue
                body = b.drop_columns(["_r"]).take(pa.array(order[cut[0]:cut[-1]]))
                for i in np.flatnonzero(np.diff(cut)):
                    held[i].append(body.slice(cut[i] - cut[0], cut[i + 1] - cut[i]))
                    n_held[i] += int(cut[i + 1] - cut[i])
                    if n_held[i] >= per:
                        flush(i)
            for i in range(len(group)):
                flush(i)
        finally:
            reader.close()
            for w in writers:
                w.close()


def _number(con, sql: str, level: int, dict_path: Path, token_path: Path) -> int:
    """Number the cluster keys of one level from `sql`, the key of every slot that has one with its protein and the
    protein's range (`key`, `pu`, `_r`), in (key, protein) order; return how many keys there are.

    A key's local id is its rank among the distinct keys in key order and its `n_slots` the slots in its run, as
    the key ranking of the release asks (no join to a key table: the runs of the sorted stream are the keys).
    Writes the keys (token_id, key, n_slots) to `dict_path` in token order, and each protein's token once,
    (_r, pu, token_id), to `token_path`."""
    reader = con.sql(sql).to_arrow_reader(io.ROW_GROUP_ROWS)
    dict_schema = pa.schema([("token_id", pa.uint64()), ("key", pa.string()), ("n_slots", pa.uint64())])
    tok_schema = pa.schema([("_r", pa.uint64()), ("pu", pa.binary()), ("token_id", pa.uint64())])
    base = level << 61
    nxt, last_key, last_pu, run = 0, None, None, 0  # keys numbered so far; the open run's key, protein, length
    dict_w = pq.ParquetWriter(dict_path, dict_schema, compression="zstd", compression_level=1)
    tok_w = pq.ParquetWriter(token_path, tok_schema, compression="zstd", compression_level=1)

    def keys_out(keys: list, locs: np.ndarray, n: np.ndarray) -> None:
        if len(keys):
            dict_w.write_table(pa.table({"token_id": pa.array(base + locs.astype(np.uint64), pa.uint64()),
                                         "key": pa.array(keys, pa.string()), "n_slots": pa.array(n, pa.uint64())}))

    try:
        for b in reader:
            n = b.num_rows
            if not n:
                continue
            key, pu = b.column("key"), b.column("pu")
            new = np.empty(n, bool)
            new[0] = last_key is None or key[0].as_py() != last_key
            new[1:] = pc.not_equal(key.slice(1), key.slice(0, n - 1)).to_numpy(zero_copy_only=False)
            same_pu = np.empty(n, bool)
            same_pu[0] = not new[0] and pu[0].as_py() == last_pu
            same_pu[1:] = pc.equal(pu.slice(1), pu.slice(0, n - 1)).to_numpy(zero_copy_only=False)
            local = nxt - 1 + np.cumsum(new)
            starts = np.flatnonzero(new)
            ends = np.r_[starts[1:], n]
            if last_key is not None and (not len(starts) or starts[0] > 0):  # the open run goes on into this batch
                run += int(starts[0]) if len(starts) else n
            if len(starts):
                if last_key is not None:
                    keys_out([last_key], np.array([nxt - 1]), np.array([run]))
                done = starts[:-1]
                keys_out(key.take(pa.array(done)).to_pylist(), local[done], ends[:-1] - done)
                run = int(ends[-1] - starts[-1])
                last_key = key[int(starts[-1])].as_py()
                nxt += len(starts)
            if nxt > 2**56:
                raise ValueError("more than 2^56 keys at one (level, source): widen the token local field")
            last_pu = pu[n - 1].as_py()
            keep = pa.array(np.flatnonzero(new | ~same_pu))  # a protein on several slots: its token once
            tok_w.write_table(pa.table({"_r": b.column("_r").take(keep).cast(pa.uint64()),
                                        "pu": pu.take(keep).cast(pa.binary()),
                                        "token_id": pa.array(base + local[keep.to_numpy()].astype(np.uint64))}))
        if last_key is not None:
            keys_out([last_key], np.array([nxt - 1]), np.array([run]))
    finally:
        reader.close()
        dict_w.close()
        tok_w.close()
    return nxt


_TOK_LISTS = (("family_tokens", "f"), ("nr50_tokens", "a"), ("nr90_tokens", "b"), ("domain_family_tokens", "d"))


def _tape_lists(rows: pa.RecordBatch) -> pa.RecordBatch:
    """The tape-token rows of `rows`, slot rows of whole tapes in (contig_uid, fidx, rn, ord) order: one row per
    tape with its tokens as lists in slot order (a tape without slots has one row with a null `ord`: empty lists)."""
    fidx = rows.column("fidx").to_numpy(zero_copy_only=False)
    rn = rows.column("rn").to_numpy(zero_copy_only=False)
    starts = np.r_[0, np.flatnonzero((fidx[1:] != fidx[:-1]) | (rn[1:] != rn[:-1])) + 1]
    slot = rows.column("ord").is_valid()
    counts = np.add.reduceat(slot.to_numpy(zero_copy_only=False).astype(np.int64), starts)
    offsets = pa.array(np.r_[0, np.cumsum(counts)], pa.int32())
    vals = rows.filter(slot)
    return pa.RecordBatch.from_arrays(
        [rows.column("contig_uid").take(pa.array(starts))] +
        [pa.ListArray.from_arrays(offsets, vals.column(c)) for _, c in _TOK_LISTS],
        names=["contig_uid"] + [n for n, _ in _TOK_LISTS])


def _tapes_of(reader) -> Iterator[pa.RecordBatch]:
    """Tape-token batches from a reader of slot rows in (contig_uid, fidx, rn, ord) order; a tape cut by a batch
    boundary is held until its last slot is read."""
    held = None
    for b in reader:
        if held is not None:
            b = pa.Table.from_batches([held, b]).combine_chunks().to_batches()[0] if b.num_rows else held
        if not b.num_rows:
            continue
        fidx = b.column("fidx").to_numpy(zero_copy_only=False)
        rn = b.column("rn").to_numpy(zero_copy_only=False)
        change = np.flatnonzero((fidx[1:] != fidx[:-1]) | (rn[1:] != rn[:-1])) + 1
        cut = int(change[-1]) if len(change) else 0
        held = b.slice(cut)
        if cut:
            yield _tape_lists(b.slice(0, cut))
    if held is not None and held.num_rows:
        yield _tape_lists(held)


def _dict_rows(con, fam_keys: str, names: dict[int, str], runs: list[tuple[str, Path]]) -> Iterator[pa.RecordBatch]:
    """The `token_dict` rows in token order: the family keys (by source, then key), then each cluster level's run
    file (already in key order) with its level and source."""
    spec = table("token_dict")
    if names:  # no library, no family key
        yield from io.duck_batches(con, f"SELECT k.token_id, 'family' AS level, s.name AS source, k.key, k.n_slots "
                                        f"FROM read_parquet({fam_keys}) k JOIN (VALUES "
                                        + ", ".join(f"({i}, {_sql_str(n)})" for i, n in names.items()) +
                                        ") s(src, name) ON s.src = k.src ORDER BY k.token_id", spec)
    for level, path in runs:
        pf = pq.ParquetFile(path)
        for i in range(pf.metadata.num_row_groups):
            t = pf.read_row_group(i)
            yield pa.RecordBatch.from_arrays(
                [t.column("token_id").combine_chunks(), pa.array([level] * t.num_rows, pa.string()),
                 pa.array([CLUSTER_SOURCE] * t.num_rows, pa.string()), t.column("key").combine_chunks(),
                 t.column("n_slots").combine_chunks()], names=spec.schema.names)


def run(root: Path, release: str, stages: list[str], gene_call_version: str = GENE_CALL_VERSION,
        libraries: list[str] | None = None, producers: list[str] | None = None,
        max_file_bytes: int = layout.MAX_FILE_BYTES, max_slots: int = MAX_SLOTS, tmp=None,
        out: Path | None = None, memory_limit: str = DUCK_MEMORY_LIMIT) -> dict:
    """Write the release's `token_dict` and the `tape_tokens` of the tape version of `producers`.

    `stages`: the release's stages (card B-25); only their tapes and domains are read, a listed stage
    without tapes raises FileNotFoundError, and the release's `release_def` is written, or checked to
    list the same stages (`catalog.define_release`). `libraries`: the domain libraries in precedence
    order (default: every library with domains in `stages`, by name). `tmp` holds DuckDB's private
    scratch directory and the builder's scratch runs (system temp by default; never the data root nor
    `out`). `memory_limit` is DuckDB's cap, which also sets the number of protein ranges (a test passes a
    smaller one, scaled with its store; card B-46).

    Inputs (tapes, domains, the release's clusters) are read from `root`; both tables and (when the
    root has none) the `release_def` are written under `out` (default `root`; card C-08). Refuses,
    before any work, when either table of this release already holds a file in `root` or `out`.
    """
    root = Path(root)
    out = root if out is None else Path(out)
    stages = catalog.check_stages(stages)
    libraries = libraries_in(root, stages) if libraries is None else list(libraries)
    if len(libraries) > 31:
        raise ValueError("more than 31 libraries: widen the token source field")
    if max_slots < 1:
        raise ValueError("max_slots must be positive")
    ledger.check_scratch(tmp, root, out)
    tv = tape_version_of(gene_call_version, producers)
    catalog.verify_release_def(root, release, stages)  # a disagreeing definition is named before a refusal
    ledger.refuse_written(root, out, [f"l3/token_dict/release={release}", f"l3/tape_tokens/tape_version={tv}/release={release}"],
                          f"tokens of release {release}")
    tape_files = _tape_files(root, tv, stages)
    catalog.define_release(root, release, stages, out=out)

    # the stages' tapes, numbered by file and row: rows with equal contig_uid keep that order
    tape_src = _scan(tape_files, "contig_uid, protein_uids",
                     "NULL::BIGINT AS fidx, NULL::BIGINT AS rn, NULL::BLOB AS contig_uid, "
                     "NULL::BLOB[] AS protein_uids", numbered=True)
    dom_parts = []
    for rank, lib in enumerate(libraries):  # a library listed twice: its first position wins, as its rank is smaller
        files = _domain_files(root, lib, stages)
        if files:
            dom_parts.append(_scan(files, f"{rank} AS rank, {_sql_str(lib)} AS lib, protein_uid, model_acc, "
                                          "dom_bitscore, env_start", ""))
    doms = " UNION ALL ".join(f"SELECT * FROM {p}" for p in dom_parts) or (
        "SELECT NULL::INTEGER AS rank, NULL::VARCHAR AS lib, NULL::BLOB AS protein_uid, NULL::VARCHAR AS model_acc, "
        "NULL::FLOAT AS dom_bitscore, NULL::UINTEGER AS env_start WHERE false")
    cdir = root / f"l2/cluster_member/release={release}"
    clusters = _scan(io.part_files(cdir), "member_uid, level, cluster_uid",
                     "NULL::BIGINT AS fidx, NULL::BIGINT AS rn, NULL::BLOB AS member_uid, NULL::VARCHAR AS level, "
                     "NULL::BLOB AS cluster_uid", numbered=True)
    scratch = Path(tempfile.mkdtemp(prefix="biobox-tokens-", dir=tmp))
    con = duckdb.connect(config={"threads": 1, "memory_limit": memory_limit,
                                 "temp_directory": str(scratch / "duckdb"), "preserve_insertion_order": True})

    def at(name: str, r: int | None = None) -> Path:
        return scratch / (f"{name}.parquet" if r is None else f"{name}-{r:05d}.parquet")

    def sql(path: Path | list[Path]) -> str:
        return f"read_parquet([{', '.join(_sql_str(p) for p in ([path] if isinstance(path, Path) else path))}])"

    try:
        con.execute(f"CREATE VIEW tapes AS SELECT * FROM {tape_src}")
        (n_tapes, n_slots), = con.execute("SELECT count(*), coalesce(sum(len(protein_uids)), 0) FROM tapes").fetchall()
        n_rows = max(n_slots, sum(pq.ParquetFile(f).metadata.num_rows for f in io.part_files(cdir)),
                     sum(pq.ParquetFile(f).metadata.num_rows for lib in libraries for f in _domain_files(root, lib, stages)))
        k = n_ranges(n_rows, memory_limit)
        ranges = range(k)
        # (1) the slots (ord 1-based; a tape without slots has one row with a null ord), the cluster rows and the
        # domain rows, each split into k ranges by the hash of the protein uid; slots without a protein go to range k
        # (an unnest in the projection streams; a lateral `unnest(...) WITH ORDINALITY` did not, at 50 MB)
        _route(con, f"SELECT CASE WHEN pu IS NULL THEN {k} ELSE hash(pu) % {k} END AS _r, * FROM ("
                    "  SELECT fidx, rn, contig_uid, CASE WHEN n > 0 THEN generate_subscripts(pus, 1) END AS ord, "
                    "  unnest(pus) AS pu FROM (SELECT fidx, rn, contig_uid, len(protein_uids) AS n, "
                    "    CASE WHEN len(protein_uids) = 0 THEN [NULL::BLOB] ELSE protein_uids END AS pus FROM tapes))",
               [at("slots", r) for r in range(k + 1)])
        _route(con, f"SELECT hash(member_uid) % {k} AS _r, member_uid, level, cluster_uid, fidx, rn FROM {clusters} "
                    "WHERE level IN ('nr50', 'nr90')", [at("members", r) for r in ranges])
        _route(con, f"SELECT hash(protein_uid) % {k} AS _r, * FROM ({doms})", [at("domains", r) for r in ranges])
        # (2) per range: each slot's canonical family and cluster keys, joined to the range's proteins only, and the
        # slots of every model with an included domain on them (card B-21), counted
        for r in ranges:
            dom, mem = sql(at("domains", r)), sql(at("members", r))
            cl = ("(SELECT member_uid, arg_max(cluster_uid, {{'f': fidx, 'r': rn}}) AS c FROM {mem} WHERE level = '{lv}' "
                  "GROUP BY member_uid)")  # one row by the primary key; were there more, the last in file order
            con.execute(
                "COPY (SELECT s.fidx, s.rn, s.contig_uid, s.ord, s.pu, f.f.src AS fam_src, f.f.acc AS fam_key, "
                f"lower(hex(c50.c)) AS nr50_key, lower(hex(c90.c)) AS nr90_key FROM {sql(at('slots', r))} s "
                # family of a protein: the smallest (library rank, -bitscore, accession, envelope start)
                "LEFT JOIN (SELECT protein_uid, arg_min({'src': rank + 1, 'acc': model_acc}, "
                f"  {{'r': rank, 'b': -dom_bitscore, 'a': model_acc, 'e': env_start}}) AS f FROM {dom} "
                "  GROUP BY protein_uid) f ON f.protein_uid = s.pu "
                f"LEFT JOIN {cl.format(mem=mem, lv='nr50')} c50 ON c50.member_uid = s.pu "
                f"LEFT JOIN {cl.format(mem=mem, lv='nr90')} c90 ON c90.member_uid = s.pu) "
                f"TO {_sql_str(at('base', r))} (FORMAT parquet, COMPRESSION zstd)")
            # every model with an included domain on a protein, once per protein; a library listed twice takes its
            # first position, as the canonical family does
            con.execute(
                f"COPY (SELECT d.src, d.model_acc AS key, count(*) AS n FROM {sql(at('slots', r))} s JOIN ("
                f"  SELECT protein_uid, min(rank) + 1 AS src, model_acc FROM {dom} GROUP BY protein_uid, lib, model_acc"
                f") d ON d.protein_uid = s.pu GROUP BY d.src, d.model_acc) TO {_sql_str(at('famcount', r))} "
                "(FORMAT parquet, COMPRESSION zstd)")
        # (3) the keys that occur on a slot, ranked within (level, source) in key order: the families by their
        # counts (a few per model and range), each cluster level by numbering the runs of its slots sorted by key
        fam_keys = _sql_str(at("famkeys"))
        con.execute(f"COPY (SELECT src, key, n_slots, local, ({LEVELS['family']}::UBIGINT << 61) | "
                    "(src::UBIGINT << 56) | local::UBIGINT AS token_id FROM ("
                    "  SELECT src, key, n_slots, row_number() OVER (PARTITION BY src ORDER BY key) - 1 AS local FROM ("
                    f"    SELECT src, key, sum(n)::BIGINT AS n_slots FROM {sql([at('famcount', r) for r in ranges])} "
                    f"    GROUP BY src, key))) TO {fam_keys} (FORMAT parquet, COMPRESSION zstd)")
        (n_fam, too_many), = con.execute(f"SELECT count(*), count(*) FILTER (WHERE local >= {2**56}) "
                                         f"FROM read_parquet({fam_keys})").fetchall()
        if too_many:
            raise ValueError("more than 2^56 keys at one (level, source): widen the token local field")
        counts = {"family": n_fam}
        bases = sql([at("base", r) for r in ranges])
        for level in ("nr50", "nr90"):
            counts[level] = _number(con, f"SELECT {level}_key AS key, pu, hash(pu) % {k} AS _r FROM {bases} "
                                         f"WHERE {level}_key IS NOT NULL ORDER BY key, pu", LEVELS[level],
                                    at(f"dict-{level}"), at(f"tokrun-{level}"))
            _route(con, f"SELECT * FROM {sql(at(f'tokrun-{level}'))}", [at(f"tok-{level}", r) for r in ranges])
        # (4) per range: every slot's tokens. Card B-39: the family tokens of each protein's distinct included
        # domains in the order of the postings' seg_idx (biobox.l3.postings._family_run: envelope start, token id,
        # library, accession), one list per protein
        for r in ranges:
            con.execute(
                "COPY (SELECT s.contig_uid, s.fidx, s.rn, s.ord, kf.token_id AS f, t50.token_id AS a, "
                f"t90.token_id AS b, coalesce(pf.toks, []::UBIGINT[]) AS d FROM {sql(at('base', r))} s "
                "LEFT JOIN (SELECT protein_uid, list_transform(list_sort(list({'e': env_start, 't': token_id, 'l': lib, "
                "  'a': model_acc})), lambda x: x.t) AS toks FROM ("
                "    SELECT d.protein_uid, d.env_start, d.lib, d.model_acc, k.token_id FROM ("
                "      SELECT protein_uid, env_start, lib, model_acc, min(rank) + 1 AS src "
                f"     FROM {sql(at('domains', r))} GROUP BY protein_uid, env_start, lib, model_acc) d "
                f"    JOIN read_parquet({fam_keys}) k ON k.src = d.src AND k.key = d.model_acc) "
                "  GROUP BY protein_uid) pf ON pf.protein_uid = s.pu "
                f"LEFT JOIN read_parquet({fam_keys}) kf ON kf.src = s.fam_src AND kf.key = s.fam_key "
                f"LEFT JOIN {sql(at('tok-nr50', r))} t50 ON t50.pu = s.pu "
                f"LEFT JOIN {sql(at('tok-nr90', r))} t90 ON t90.pu = s.pu) "
                f"TO {_sql_str(at('slot-tokens', r))} (FORMAT parquet, COMPRESSION zstd)")
        con.execute("COPY (SELECT contig_uid, fidx, rn, ord, NULL::UBIGINT AS f, NULL::UBIGINT AS a, NULL::UBIGINT AS b, "
                    f"[]::UBIGINT[] AS d FROM {sql(at('slots', k))}) TO {_sql_str(at('slot-tokens', k))} "
                    "(FORMAT parquet, COMPRESSION zstd)")  # slots without a protein have no token

        # (5) the dictionary in token order, and one row of lists per tape from its slots sorted back into place
        spec = table("token_dict")
        n_dict = len(io.write_parts(
            _dict_rows(con, fam_keys, {i + 1: lib for i, lib in enumerate(libraries)},
                       [(level, at(f"dict-{level}")) for level in ("nr50", "nr90")]),
            spec, out / f"l3/token_dict/release={release}", EST_DICT_ROW, max_file_bytes,
            statistics=["n_slots", "key"]))  # key: token_of's point lookup by (level, key), card B-34
        est_tok = EST_TOKEN_ROW + EST_TOKEN_SLOT * max(1, n_slots // max(1, n_tapes))
        reader = con.sql(f"SELECT * FROM {sql([at('slot-tokens', r) for r in range(k + 1)])} "
                         "ORDER BY contig_uid, fidx, rn, ord").to_arrow_reader(io.ROW_GROUP_ROWS)
        try:
            n_tok = len(io.write_parts(_tapes_of(reader), table("tape_tokens"),
                                       out / f"l3/tape_tokens/tape_version={tv}/release={release}", est_tok,
                                       max_file_bytes, max_slots=max_slots, slots_column="nr90_tokens"))
        finally:
            reader.close()
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)
    return {"tokens": {lvl: counts.get(lvl, 0) for lvl in LEVELS}, "contigs": n_tapes, "files": n_dict + n_tok}
