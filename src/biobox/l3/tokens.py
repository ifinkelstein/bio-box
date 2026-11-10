"""Release-time tokens for the tapes (PLAN.md §4.1 token namespace; data_layer_and_tests.md §1 L3; card B-09).

Tapes carry stable keys (`protein_uid`); a release maps each slot to tokens at three levels:

    family  the canonical domain of the protein: the first library in `libraries` order with
            an included domain, and in it the domain with the highest bitscore (ties: model
            accession, then envelope start)
    nr50    the protein's NR50 cluster in this release
    nr90    the protein's NR90 cluster in this release

`token_id` is uint64: level (3 bits) | source (5 bits) | local (56 bits). PLAN.md §4.1 gives
32 bits with 24 local and widens to 64 when counts require it; NR90 clusters pass 2^24 by
S5, so the release line starts wide. Source is the library's 1-based position in
`libraries` for family tokens and 0 for clusters; local is the rank of the key within
(level, source) among the keys that occur on at least one slot, in ascending string order of the
key (the model accession for families, the cluster uid as lowercase hex for clusters), so ids are
dense and reproducible but release-local (never a join key outside the release).

Writes, under `<root>/l3/`, through `biobox.io.write_table`:

    token_dict/release=<R>/part-<k>.parquet                     token_id -> level, source, key, n_slots
    tape_tokens/tape_version=<tv>/release=<R>/part-<k>.parquet  per contig: one list per level

The tape version is that of `producers` (default: the gene-call producer alone, card B-07); a
non-CDS slot has no protein and so no token (a null in each list), and a tape without slots has
three empty lists.

Bounded memory (AGENTS §4.1): one DuckDB connection (`DUCK_MEMORY_LIMIT`, one thread, spilling
into a private scratch directory) unnests the tapes' `protein_uids` with their ordinals, joins them
to the canonical families and the clusters, ranks the keys and aggregates the tokens back into one
list per tape; the results stream through `io.duck_batches` into `io.write_table`. A row group of
`tape_tokens` holds whole contigs and at most `max_slots` slots (AGENTS §4.6). The list aggregate
collects (ordinal, tokens) structs and sorts them per tape: DuckDB spills an unordered `list()` but
not `list(... ORDER BY ...)`.
"""

from __future__ import annotations

import heapq
import shutil
import tempfile
from collections.abc import Iterator
from pathlib import Path

import duckdb
import pyarrow.parquet as pq

from biobox import io, layout
from biobox.l1.genecall import DUCK_MEMORY_LIMIT, GENE_CALL_VERSION
from biobox.l3.tapes import MAX_SLOTS, tape_version
from biobox.registry.schemas import table

LEVELS = {"family": 1, "nr50": 2, "nr90": 3}  # 4 is reserved for subfamilies
CLUSTER_SOURCE = "mmseqs2"
EST_DICT_ROW = 60  # zstd estimate per token_dict row; cuts parts only
EST_TOKEN_ROW, EST_TOKEN_SLOT = 100, 30  # per tape_tokens row and per mean slot; cuts parts only


def token_id(level: str, source: int, local: int) -> int:
    if not 0 <= source < 32 or not 0 <= local < 2**56:
        raise ValueError(f"token ({level}, {source}, {local}) out of range")
    return LEVELS[level] << 61 | source << 56 | local


def token_level(tid: int) -> str:
    return {v: k for k, v in LEVELS.items()}[tid >> 61]


def tape_version_of(gene_call_version: str = GENE_CALL_VERSION, producers: list[str] | None = None) -> str:
    """The tape version of `producers`, by default the gene-call producer alone."""
    return tape_version([gene_call_version] if producers is None else producers)


def _tape_files(root: Path, tv: str) -> list[Path]:
    return io.part_files(Path(root) / "l3/tapes", f"tape_version={tv}/stage=*/part-*.parquet")


def read_tapes(root: Path, gene_call_version: str = GENE_CALL_VERSION, columns: list[str] | None = None,
               producers: list[str] | None = None) -> list[dict]:
    """Every tape of the tape version as a list (a test and `pairs` helper; builders use `iter_tapes`)."""
    tv = tape_version_of(gene_call_version, producers)
    return [r for f in _tape_files(root, tv) for r in pq.read_table(f, columns=columns).to_pylist()]


def iter_tapes(root: Path, tape_version: str, columns: list[str] | None = None) -> Iterator[dict]:
    """The tapes of `tape_version` one row at a time in `contig_uid` order, over every stage's files.

    Each file is read one row group at a time and the files are merged on `contig_uid` (equal uids
    in file path order), so one row group per file is held.
    """
    def rows(f: Path) -> Iterator[dict]:
        pf = pq.ParquetFile(f)
        cols = None if columns is None else list(dict.fromkeys(["contig_uid", *columns]))
        for i in range(pf.metadata.num_row_groups):
            yield from pf.read_row_group(i, columns=cols).to_pylist()

    keep = None if columns is None else set(columns)
    for r in heapq.merge(*(rows(f) for f in _tape_files(root, tape_version)), key=lambda r: r["contig_uid"]):
        yield r if keep is None else {k: v for k, v in r.items() if k in keep}


def libraries_in(root: Path) -> list[str]:
    return sorted({p.name.removeprefix("library=") for p in (Path(root) / "l1/domains").glob("library=*")})


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


def _write_parts(con, sql: str, spec_name: str, base: Path, est_row: int, max_file_bytes: int, **kw) -> int:
    """Stream `sql` into `base/part-<k>.parquet` through `io.write_parts`, rows of `est_row` estimated
    bytes each, so every part holds max(1, max_file_bytes // est_row) rows. At least part-0 is written."""
    spec = table(spec_name)
    return len(io.write_parts(io.duck_batches(con, sql, spec), spec, base, est_row, max_file_bytes, **kw))


def run(root: Path, release: str, gene_call_version: str = GENE_CALL_VERSION,
        libraries: list[str] | None = None, producers: list[str] | None = None,
        max_file_bytes: int = layout.MAX_FILE_BYTES, max_slots: int = MAX_SLOTS, tmp=None) -> dict:
    """Write the release's `token_dict` and the `tape_tokens` of the tape version of `producers`.

    `libraries`: the domain libraries in precedence order (default: every library under `l1/domains`,
    by name). `tmp` holds DuckDB's private scratch directory (system temp by default; never the data root).
    """
    root = Path(root)
    libraries = libraries_in(root) if libraries is None else list(libraries)
    if len(libraries) > 31:
        raise ValueError("more than 31 libraries: widen the token source field")
    if max_slots < 1:
        raise ValueError("max_slots must be positive")
    if tmp is not None and Path(tmp).resolve().is_relative_to(root.resolve()):
        raise ValueError(f"{tmp}: the scratch directory may not lie under the data root")
    tv = tape_version_of(gene_call_version, producers)

    # every stage's tapes, numbered by file and row: rows with equal contig_uid keep that order
    tape_src = _scan(_tape_files(root, tv), "contig_uid, protein_uids",
                     "NULL::BIGINT AS fidx, NULL::BIGINT AS rn, NULL::BLOB AS contig_uid, "
                     "NULL::BLOB[] AS protein_uids", numbered=True)
    dom_parts = []
    for rank, lib in enumerate(libraries):  # a library listed twice: its first position wins, as its rank is smaller
        files = io.part_files(root / "l1/domains", f"library={lib}/stage=*/part-*.parquet")
        if files:
            dom_parts.append(_scan(files, f"{rank} AS rank, protein_uid, model_acc, dom_bitscore, env_start", ""))
    doms = " UNION ALL ".join(f"SELECT * FROM {p}" for p in dom_parts) or (
        "SELECT NULL::INTEGER AS rank, NULL::BLOB AS protein_uid, NULL::VARCHAR AS model_acc, "
        "NULL::FLOAT AS dom_bitscore, NULL::UINTEGER AS env_start WHERE false")
    cdir = root / f"l2/cluster_member/release={release}"
    clusters = _scan(io.part_files(cdir), "member_uid, level, cluster_uid",
                     "NULL::BIGINT AS fidx, NULL::BIGINT AS rn, NULL::BLOB AS member_uid, NULL::VARCHAR AS level, "
                     "NULL::BLOB AS cluster_uid", numbered=True)
    names = ", ".join(f"({s}, {_sql_str(n)})" for s, n in [(0, CLUSTER_SOURCE)] +
                      [(i + 1, lib) for i, lib in enumerate(libraries)])
    level_names = ", ".join(f"({v}, {_sql_str(k)})" for k, v in LEVELS.items())
    lv_fam, lv_50, lv_90 = LEVELS["family"], LEVELS["nr50"], LEVELS["nr90"]

    scratch = Path(tempfile.mkdtemp(prefix="biobox-tokens-", dir=tmp))
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb"), "preserve_insertion_order": True})
    slots, keys = _sql_str(scratch / "slots.parquet"), _sql_str(scratch / "keys.parquet")
    try:
        con.execute(f"CREATE VIEW tapes AS SELECT * FROM {tape_src}")
        # family of a protein: the smallest (library rank, -bitscore, accession, envelope start)
        con.execute(f"CREATE VIEW fam AS SELECT protein_uid, arg_min({{'src': rank + 1, 'acc': model_acc}}, "
                    f"{{'r': rank, 'b': -dom_bitscore, 'a': model_acc, 'e': env_start}}) AS f FROM ({doms}) "
                    "GROUP BY protein_uid")
        # cluster of a protein at a level (one row by the primary key; were there more, the last in file order)
        con.execute(f"CREATE VIEW cl AS SELECT member_uid, level, arg_max(cluster_uid, {{'f': fidx, 'r': rn}}) AS c "
                    f"FROM {clusters} WHERE level IN ('nr50', 'nr90') GROUP BY member_uid, level")
        # one row per slot (ord 1-based), or one with a null ord for a tape without slots
        con.execute(
            "COPY (SELECT t.fidx, t.rn, t.contig_uid, CASE WHEN t.n > 0 THEN u.ord END AS ord, "
            "f.f.src AS fam_src, f.f.acc AS fam_key, lower(hex(c50.c)) AS nr50_key, lower(hex(c90.c)) AS nr90_key "
            "FROM (SELECT fidx, rn, contig_uid, len(protein_uids) AS n, "
            "      CASE WHEN len(protein_uids) = 0 THEN [NULL::BLOB] ELSE protein_uids END AS pus FROM tapes) t "
            "CROSS JOIN unnest(t.pus) WITH ORDINALITY AS u(pu, ord) "
            "LEFT JOIN fam f ON f.protein_uid = u.pu "
            "LEFT JOIN cl c50 ON c50.member_uid = u.pu AND c50.level = 'nr50' "
            "LEFT JOIN cl c90 ON c90.member_uid = u.pu AND c90.level = 'nr90') "
            f"TO {slots} (FORMAT parquet, COMPRESSION zstd)")
        # the keys that occur on a slot, ranked within (level, source)
        con.execute(
            "COPY (SELECT lv, src, key, n_slots, local, (lv::UBIGINT << 61) | (src::UBIGINT << 56) | local::UBIGINT "
            "AS token_id FROM ("
            "  SELECT lv, src, key, count(*) AS n_slots, "
            "  row_number() OVER (PARTITION BY lv, src ORDER BY key) - 1 AS local FROM ("
            f"    SELECT {lv_fam} AS lv, fam_src AS src, fam_key AS key FROM read_parquet({slots}) WHERE fam_key IS NOT NULL"
            f"    UNION ALL SELECT {lv_50}, 0, nr50_key FROM read_parquet({slots}) WHERE nr50_key IS NOT NULL"
            f"    UNION ALL SELECT {lv_90}, 0, nr90_key FROM read_parquet({slots}) WHERE nr90_key IS NOT NULL) "
            f"  GROUP BY lv, src, key)) TO {keys} (FORMAT parquet, COMPRESSION zstd)")
        (too_many,), = con.execute(f"SELECT count(*) FROM read_parquet({keys}) WHERE local >= {2**56}").fetchall()
        if too_many:
            raise ValueError("more than 2^56 keys at one (level, source): widen the token local field")
        counts = dict(con.execute(f"SELECT lv, count(*) FROM read_parquet({keys}) GROUP BY lv").fetchall())
        (n_tapes, n_slots), = con.execute("SELECT count(*), coalesce(sum(len(protein_uids)), 0) FROM tapes").fetchall()

        n_dict = _write_parts(
            con, f"SELECT k.token_id, l.level, s.name AS source, k.key, k.n_slots FROM read_parquet({keys}) k "
                 f"JOIN (VALUES {level_names}) l(lv, level) ON l.lv = k.lv "
                 f"JOIN (VALUES {names}) s(src, name) ON s.src = k.src ORDER BY k.token_id",
            "token_dict", root / f"l3/token_dict/release={release}", EST_DICT_ROW, max_file_bytes,
            statistics=["n_slots"])
        est_tok = EST_TOKEN_ROW + EST_TOKEN_SLOT * max(1, n_slots // max(1, n_tapes))
        n_tok = _write_parts(
            con, "SELECT contig_uid, "
                 "coalesce(list_transform(l, lambda x: x.f), []::UBIGINT[]) AS family_tokens, "
                 "coalesce(list_transform(l, lambda x: x.a), []::UBIGINT[]) AS nr50_tokens, "
                 "coalesce(list_transform(l, lambda x: x.b), []::UBIGINT[]) AS nr90_tokens FROM ("
                 "  SELECT s.fidx, s.rn, s.contig_uid, list_sort(list({'o': s.ord, 'f': kf.token_id, "
                 "    'a': k50.token_id, 'b': k90.token_id}) FILTER (WHERE s.ord IS NOT NULL)) AS l "
                 f"  FROM read_parquet({slots}) s "
                 f"  LEFT JOIN read_parquet({keys}) kf ON kf.lv = {lv_fam} AND kf.src = s.fam_src AND kf.key = s.fam_key "
                 f"  LEFT JOIN read_parquet({keys}) k50 ON k50.lv = {lv_50} AND k50.key = s.nr50_key "
                 f"  LEFT JOIN read_parquet({keys}) k90 ON k90.lv = {lv_90} AND k90.key = s.nr90_key "
                 "  GROUP BY s.fidx, s.rn, s.contig_uid) "
                 "ORDER BY contig_uid, fidx, rn",
            "tape_tokens", root / f"l3/tape_tokens/tape_version={tv}/release={release}", est_tok, max_file_bytes,
            max_slots=max_slots, slots_column="nr90_tokens")
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)
    return {"tokens": {lvl: counts.get(lv, 0) for lvl, lv in LEVELS.items()}, "contigs": n_tapes,
            "files": n_dict + n_tok}
