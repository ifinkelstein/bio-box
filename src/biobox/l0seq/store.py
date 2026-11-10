"""Write and read the L0-seq store (operations O-INGEST, O-SLICE, O-EXTRACT).

Layout under `<root>/` (planning/data_layer_and_tests.md §2A.1, §3.1):

    l0/l0_manifest.parquet
    l0seq/store_params.json
    l0seq/seq_shards/stage=<S>/shard-<NNNN>.fa.gz   (+ .fa.gz.fai, .fa.gz.gzi)
    l0seq/contig_occurrence/stage=<S>/source=<name>/part-<k>.parquet
    l0seq/contig_index/stage=<S>/part-<k>.parquet
    l0seq/seq_locator/stage=<S>/part-<k>.parquet

Each table is cut into parts at `max_file_bytes` of estimated size (`biobox.io.write_parts`, card
B-19); a small stage is one `part-0.parquet`. Readers take every part, in part order.

Sequences are bgzip FASTA shards (DEC-17, after K-13): each stored byte sequence once,
named by its `raw_digest` hex. One ingest batch writes one shard (card B-15, TASKS §2): a
batch is a maximal run of source objects, in object sha256 order, whose new sequences total
at most `shard_bases` bases (one larger object is a batch of its own); records are sorted by
`raw_digest` inside a shard and there is no order across shards. Ingest holds one batch's
sequences at a time. samtools, pysam and MMseqs2 read the shards as they are.

A stage never rewrites an earlier stage's files (AGENTS §3.10): a sequence already stored
adds only occurrence rows, and a stage's `seq_locator` file holds only the sequences first
stored in that stage. `l0_manifest` and `store_params.json` are small metadata files
rewritten per stage.

`SeqStore` loads neither index (card B-13): `seq_locator` and `contig_index` are looked up
through `KeyLookup` (footers once, one row group per key) and each shard's `.fai` is bisected
on disk (`BgzfFasta`).
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from biobox import fasta, ids, io, layout, ledger, provenance
from biobox.l0seq.bgzf import BgzfFasta, write_bgzf_fasta
from biobox.l0seq.lookup import KeyLookup
from biobox.registry.schemas import table

STORE_FORMAT = "bgzf-fasta-1"
SHARD_BASES = 2**30  # ~0.33 GB on disk at 0.31 bytes/bp: inside the 256 MB-1 GB band (LY-04)
LINE_WIDTH = 80
BGZF_LEVEL = 6
DUCK_MEMORY_LIMIT = "1GB"  # ingest's DuckDB cap; larger sorts spill into the private scratch directory
SCRATCH_FLUSH_ROWS = 65_536  # pass-1 rows buffered before a scratch row group is written
# Estimated zstd bytes per row; they only cut parts at `max_file_bytes`. Estimates from the column
# widths (random 16- and 32-byte keys do not compress); none of these tables has a measurement yet
# (planning/measurements/2026-11-09_design_stress_test.md §3 and §6a list none).
EST_OCCURRENCE_ROW = 120  # estimate: contig_occurrence (accession and header strings, raw_digest, contig_uid)
EST_INDEX_ROW = 56  # estimate: contig_index (contig_uid 16 + raw_digest 32 + length and flag)
EST_LOCATOR_ROW = 40  # estimate: seq_locator (raw_digest 32 + dictionary-coded shard path + length)

_COMP = bytes.maketrans(b"ACGTRYKMBVDHNacgtrykmbvdhn", b"TGCAYRMKVBHDNtgcayrmkvbhdn")


def revcomp(seq: bytes) -> bytes:
    return seq.translate(_COMP)[::-1]


@dataclass(frozen=True)
class SourceObject:
    path: Path
    source_name: str
    source_release: str
    license: str = "unknown"
    redistributable: bool = False
    precedence: int = 100  # lower wins when choosing a contig's default variant
    retrieved_at: datetime = field(default_factory=lambda: datetime(1970, 1, 1, tzinfo=timezone.utc))


def _write(tbl: pa.Table, path: Path, row_group_size: int, statistics: list[str] | bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(tbl, path, compression="zstd", row_group_size=row_group_size, write_statistics=statistics)


LOCATOR_DIR = "l0seq/seq_locator"
INDEX_DIR = "l0seq/contig_index"


def locator_files(root: Path) -> list[Path]:
    """The store's `seq_locator` files (the parts of each stage that stored a sequence), in stage order
    and then part order."""
    return io.part_files(Path(root) / LOCATOR_DIR, "stage=*/part-*.parquet")


def index_files(root: Path) -> list[Path]:
    return io.part_files(Path(root) / INDEX_DIR, "stage=*/part-*.parquet")


def has_sequences(root: Path) -> bool:
    """True if `root` holds a sequence store (at least one stage's `seq_locator`)."""
    return bool(locator_files(root))


def _sql_path(p: Path) -> str:
    return "'" + str(p).replace("'", "''") + "'"


def _sql_files(files: list[Path]) -> str:
    return "[" + ", ".join(_sql_path(f) for f in files) + "]"


# pass 1's private scratch rows, one per record: the occurrence columns plus what DuckDB decides on
_SEEN = pa.schema([
    ("source", pa.string()),  # source_name, the contig_occurrence partition
    ("source_uid", pa.uint64()),
    ("source_release", pa.string()),
    ("l0_object_sha256", pa.binary(32)),
    ("record_ordinal", pa.uint64()),
    ("accession", pa.string()),
    ("original_header", pa.string()),
    ("raw_digest", pa.binary(32)),
    ("contig_uid", pa.binary(16)),
    ("line_width", pa.uint32()),
    ("contig_len", pa.int64()),
    ("precedence", pa.int64()),
    ("obj_idx", pa.uint32()),  # rank of the object in sha256 order
])


class _Scratch:
    """An unsorted Parquet file written in row groups of `SCRATCH_FLUSH_ROWS` rows."""

    def __init__(self, path: Path, schema: pa.Schema):
        self._schema = schema
        self._cols: list[list] = [[] for _ in schema]
        self._writer = pq.ParquetWriter(path, schema, compression="zstd", write_statistics=False)

    def add(self, row: tuple) -> None:
        for col, v in zip(self._cols, row):
            col.append(v)
        if len(self._cols[0]) >= SCRATCH_FLUSH_ROWS:
            self.flush()

    def flush(self) -> None:
        if self._cols[0]:
            arrays = [pa.array(c, type=f.type) for c, f in zip(self._cols, self._schema)]
            self._writer.write_batch(pa.RecordBatch.from_arrays(arrays, schema=self._schema))
            self._cols = [[] for _ in self._schema]

    def close(self) -> None:
        self.flush()
        self._writer.close()


def _object_order(objects: list[SourceObject]) -> list[tuple[bytes, SourceObject]]:
    """(sha256, object) in ascending sha256; ties (the same bytes listed twice) by the object's fields."""
    keyed = [(ids.sha256_file(o.path), o) for o in objects]
    return sorted(keyed, key=lambda t: (t[0], t[1].source_name, t[1].source_release, t[1].precedence,
                                        Path(t[1].path).name, str(t[1].path)))


def _batches(per_object: list[tuple[int, int]], shard_bases: int) -> list[list[int]]:
    """Objects (index, new bases), in sha256 order, cut into maximal runs of at most `shard_bases`
    new bases; an object above `shard_bases` is a batch of its own. Objects that store nothing new
    are in no batch (TASKS §2, card B-15)."""
    batches: list[list[int]] = []
    held = 0
    for i, bases in per_object:
        if not batches or held + bases > shard_bases:
            batches.append([])
            held = 0
        batches[-1].append(i)
        held += bases
    return batches


def _wanted(con, batch_rows: int = io.ROW_GROUP_ROWS):
    """(obj_idx, record_ordinal, raw_digest, contig_len) of every first occurrence, in object then
    record order."""
    reader = con.sql("SELECT obj_idx, record_ordinal, raw_digest, contig_len FROM new_seq "
                     "ORDER BY obj_idx, record_ordinal").to_arrow_reader(batch_rows)
    try:
        for b in reader:
            yield from zip(*(b.column(i).to_pylist() for i in range(4)))
    finally:
        reader.close()


def _stat(path: Path) -> tuple[int, int]:
    st = Path(path).stat()
    return st.st_size, st.st_mtime_ns


def stage_dirs(stage: str) -> list[str]:
    """The root-relative directories that hold one stage's ingest files (the ledger's included)."""
    return [f"l0seq/seq_shards/stage={stage}", f"l0seq/contig_occurrence/stage={stage}", f"{INDEX_DIR}/stage={stage}",
            f"{LOCATOR_DIR}/stage={stage}", f"ledger/stage={stage}/step=ingest"]


def stage_files(base: Path, stage: str) -> list[Path]:
    """Every file under `base` in the stage's ingest directories, a writer's temporary files included."""
    return sorted(p for d in stage_dirs(stage) if (Path(base) / d).is_dir()
                  for p in (Path(base) / d).rglob("*") if not p.is_dir())


def _no_hook(point: str) -> None:
    pass


def ingest(
    objects: list[SourceObject],
    root: Path,
    stage: str,
    shard_bases: int = SHARD_BASES,
    line_width: int = LINE_WIDTH,
    tmp: Path | None = None,
    max_file_bytes: int = layout.MAX_FILE_BYTES,
    out: Path | None = None,
    _hook: Callable[[str], None] | None = None,
) -> dict:
    """Ingest FASTA source objects for one stage. Returns counts for the ledger.

    Batch ingest (card B-15; AGENTS §4.1, §4.3). The objects are taken in ascending sha256, so no
    output depends on the order of `objects`. Pass 1 streams every record once into a private
    scratch table (sequences are hashed and dropped). One DuckDB connection then decides which
    digests are new (an anti-join against the earlier stages' `seq_locator` files), the first
    occurrence of each by (object sha256, record ordinal), and the stage's `contig_index` rows.
    Pass 2 re-reads the objects batch by batch (`_batches`) and writes one shard per batch, its new
    sequences sorted by `raw_digest`; only one batch's sequences are held at a time. Pass 2 refuses a
    source that changed since pass 1: an object whose size or mtime differs, and a kept record whose
    length or `raw_digest` differs (card B-23). `tmp` holds the scratch directory (system temp by default; never under the data root
    or `out`). The tables are cut into parts at `max_file_bytes` of estimated size (`io.write_parts`).

    Where files go, and a rerun after a crash (card B-23). The earlier stages are read from `root`;
    this stage's files are written under `out` (default `root`) at their root-relative paths. Under
    `bio-box run-step` pass `out=$BIOBOX_STAGING` and `root=$BIOBOX_ROOT`: run-step publishes the
    stage from staging only when ingest succeeded and empties staging before the next attempt, so a
    killed ingest is rerun from the start into a clean staging directory and gives the same bytes.
    A direct call (`out` is the root) needs a clean stage directory: ingest refuses, before it reads
    a source, when any file (a temporary one included) is in the stage's directories (`stage_dirs`)
    under `root` or `out`. After a direct ingest died, remove those directories and run it again; the
    `l0_manifest` is written last, so a stage that died earlier has not touched it. Nothing here
    overwrites a file: `io.write_table` still refuses to. `store_params.json` is written under `out`
    only when the root has none; when it has one, the parameters must equal it.

    `_hook` (tests only) is called with each point at which a kill leaves the stage half written:
    "pass1", "shard:<j>" after each shard, "contig_occurrence", "contig_index", "seq_locator",
    "ledger", "l0_manifest"; an exception it raises is not caught.
    """
    root = Path(root)
    out = root if out is None else Path(out)
    hook = _hook or _no_hook
    if tmp is not None:
        for base, what in ((root, "the data root"), (out, "the output directory")):
            if Path(tmp).resolve().is_relative_to(base.resolve()):
                raise ValueError(f"{tmp}: the scratch directory may not lie under {what}")
    for base in dict.fromkeys([root.resolve(), out.resolve()]):
        left = stage_files(base, stage)
        if left:
            raise FileExistsError(
                f"{left[0]}: stage {stage} already has {len(left)} ingest file(s) under {base}; ingest "
                f"writes a stage once. A direct call needs a clean stage directory: if an earlier ingest "
                f"died, remove {', '.join(stage_dirs(stage))} and run it again, or ingest under run-step")
    params_path = root / "l0seq" / "store_params.json"
    params = {"format": STORE_FORMAT, "shard_bases": shard_bases, "line_width": line_width, "bgzf_level": BGZF_LEVEL}
    if params_path.exists():
        if json.loads(params_path.read_text()) != params:
            raise ValueError("store parameters differ from the existing store")
    else:
        (out / "l0seq").mkdir(parents=True, exist_ok=True)
        (out / "l0seq" / "store_params.json").write_text(json.dumps(params, sort_keys=True))

    ordered = _object_order(objects)
    old_loc, old_idx = locator_files(root), index_files(root)  # the earlier stages', read before writing
    scratch = Path(tempfile.mkdtemp(prefix="biobox-ingest-", dir=tmp))
    try:
        # pass 1: every record once; occurrence and decision columns to scratch, sequences dropped
        seen_path = scratch / "seen.parquet"
        seen = _Scratch(seen_path, _SEEN)
        n_in = [0] * len(ordered)
        stats = [_stat(obj.path) for _, obj in ordered]
        try:
            for i, (sha, obj) in enumerate(ordered):
                src_uid = ids.source_uid(obj.source_name, obj.source_release)
                for k, rec in enumerate(fasta.iter_records(Path(obj.path))):
                    seen.add((obj.source_name, src_uid, obj.source_release, sha, k, rec.accession, rec.header,
                              ids.raw_digest(rec.seq), ids.contig_uid(rec.seq), rec.line_width, len(rec.seq),
                              obj.precedence, i))
                    n_in[i] = k + 1
                rec = None
        finally:
            seen.close()
        hook("pass1")
        con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                     "temp_directory": str(scratch / "duckdb"), "preserve_insertion_order": True})
        try:
            return _decide_and_write(con, root, out, stage, ordered, n_in, stats, seen_path, old_loc, old_idx,
                                     scratch, shard_bases, line_width, max_file_bytes, hook)
        finally:
            con.close()
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def _decide_and_write(con, root: Path, out: Path, stage: str, ordered: list[tuple[bytes, SourceObject]],
                      n_in: list[int], stats: list[tuple[int, int]], seen_path: Path, old_loc: list[Path],
                      old_idx: list[Path], scratch: Path, shard_bases: int, line_width: int, max_file_bytes: int,
                      hook: Callable[[str], None]) -> dict:
    con.execute(f"CREATE VIEW seen AS SELECT * FROM read_parquet({_sql_path(seen_path)})")
    con.execute("CREATE VIEW old_loc AS " + (
        f"SELECT raw_digest FROM read_parquet({_sql_files(old_loc)}, union_by_name = true)" if old_loc
        else "SELECT NULL::BLOB AS raw_digest WHERE false"))
    con.execute("CREATE VIEW old_idx AS " + (
        f"SELECT contig_uid, raw_digest, is_default_variant FROM read_parquet({_sql_files(old_idx)}, "
        "union_by_name = true)" if old_idx
        else "SELECT NULL::BLOB AS contig_uid, NULL::BLOB AS raw_digest, NULL::BOOLEAN AS is_default_variant WHERE false"))
    # new digests: not in any earlier stage, each at its first occurrence by (object sha256, ordinal)
    new_path = scratch / "new_seq.parquet"
    con.execute(f"""
        COPY (SELECT raw_digest, contig_len, obj_idx, record_ordinal FROM (
                SELECT s.raw_digest, s.contig_len, s.obj_idx, s.record_ordinal,
                       row_number() OVER (PARTITION BY s.raw_digest ORDER BY s.obj_idx, s.record_ordinal) AS rn
                FROM seen s ANTI JOIN old_loc l ON s.raw_digest = l.raw_digest)
              WHERE rn = 1 ORDER BY obj_idx, record_ordinal)
        TO {_sql_path(new_path)} (FORMAT parquet)""")
    con.execute(f"CREATE VIEW new_seq AS SELECT * FROM read_parquet({_sql_path(new_path)})")
    kept = [0] * len(ordered)
    per_object = []
    for i, n, bases in con.execute("SELECT obj_idx, count(*), sum(contig_len) FROM new_seq "
                                   "GROUP BY obj_idx ORDER BY obj_idx").fetchall():
        kept[i] = n
        per_object.append((i, int(bases)))
    batches = _batches(per_object, shard_bases)
    shard_of = [f"l0seq/seq_shards/stage={stage}/shard-{j:04d}.fa.gz" for j in range(len(batches))]
    con.register("batch_of", pa.table({
        "obj_idx": pa.array([i for b in batches for i in b], pa.uint32()),
        "shard": pa.array([shard_of[j] for j, b in enumerate(batches) for _ in b], pa.string())}))

    # pass 2: one batch's first occurrences at a time, one shard per batch
    wanted = _wanted(con)
    nxt = next(wanted, None)
    try:
        for j, batch in enumerate(batches):
            seqs: dict[bytes, bytes] = {}
            for i in batch:
                path = Path(ordered[i][1].path)
                if _stat(path) != stats[i]:
                    raise ValueError(f"{path}: changed since the first pass (size or mtime)")
                for k, rec in enumerate(fasta.iter_records(path)):
                    if nxt is None or nxt[0] != i:
                        break  # nothing more to keep from this object
                    if nxt[1] == k:
                        # stored bytes are identified by their digest (AGENTS §3.2): hash again (about
                        # 3 s per 240 s at 1,000x, card B-23), after the cheap length check
                        if len(rec.seq) != nxt[3] or ids.raw_digest(rec.seq) != nxt[2]:
                            raise ValueError(f"{path}: record {k} changed since the first pass")
                        seqs[nxt[2]] = rec.seq
                        nxt = next(wanted, None)
                rec = None
                if nxt is not None and nxt[0] == i:
                    raise ValueError(f"{path}: fewer records than in the first pass")
            write_bgzf_fasta(((rd.hex(), seqs[rd]) for rd in sorted(seqs)), out / shard_of[j], line_width, BGZF_LEVEL)
            del seqs
            hook(f"shard:{j}")
    finally:
        wanted.close()
    n_new = sum(kept)

    # contig_occurrence: partitioned by source, sorted by (object sha, ordinal)
    occ_spec = table("contig_occurrence")
    occ_cols = ", ".join(f'"{c}"' for c in occ_spec.schema.names)
    for (src,) in con.execute("SELECT DISTINCT source FROM seen ORDER BY source").fetchall():
        sql = (f"SELECT {occ_cols} FROM seen WHERE source = {_sql_path(src)} "
               "ORDER BY l0_object_sha256, record_ordinal, obj_idx")
        io.write_parts(io.duck_batches(con, sql, occ_spec), occ_spec,
                       out / f"l0seq/contig_occurrence/stage={stage}/source={src}", EST_OCCURRENCE_ROW, max_file_bytes)
    hook("contig_occurrence")

    # contig_index: new (contig_uid, raw_digest) pairs; the default is the variant with the lowest
    # (precedence, object sha256, record ordinal) of the stage, unless an earlier stage gave one
    idx_spec = table("contig_index")
    idx_sql = """
        WITH pairs AS (SELECT contig_uid, raw_digest, min(contig_len) AS contig_len FROM seen
                       GROUP BY contig_uid, raw_digest),
             ranked AS (SELECT contig_uid, raw_digest, row_number() OVER (PARTITION BY contig_uid
                            ORDER BY precedence, l0_object_sha256, record_ordinal, obj_idx) AS rn FROM seen),
             best AS (SELECT contig_uid, raw_digest AS best_digest FROM ranked WHERE rn = 1),
             old_default AS (SELECT DISTINCT contig_uid FROM old_idx WHERE is_default_variant)
        SELECT p.contig_uid, p.raw_digest, p.contig_len,
               (d.contig_uid IS NULL AND b.best_digest = p.raw_digest) AS is_default_variant
        FROM pairs p JOIN best b ON b.contig_uid = p.contig_uid
             LEFT JOIN old_default d ON d.contig_uid = p.contig_uid
        WHERE NOT EXISTS (SELECT 1 FROM old_idx o WHERE o.contig_uid = p.contig_uid AND o.raw_digest = p.raw_digest)
        ORDER BY p.contig_uid, p.raw_digest"""
    (n_idx,), = con.execute(f"SELECT count(*) FROM ({idx_sql})").fetchall()
    if n_idx:
        io.write_parts(io.duck_batches(con, idx_sql, idx_spec), idx_spec,
                       out / f"{INDEX_DIR}/stage={stage}", EST_INDEX_ROW, max_file_bytes)
    hook("contig_index")

    # seq_locator: this stage's new sequences only, sorted by raw_digest; never rewritten (AGENTS §3.10)
    if n_new:
        loc_spec = table("seq_locator")
        io.write_parts(io.duck_batches(con, "SELECT n.raw_digest, b.shard, n.contig_len FROM new_seq n "
                                            "JOIN batch_of b USING (obj_idx) ORDER BY n.raw_digest", loc_spec),
                       loc_spec, out / f"{LOCATOR_DIR}/stage={stage}", EST_LOCATOR_ROW, max_file_bytes)

    hook("seq_locator")

    # ledger: per object, records in; kept = sequences first stored from it; merged = the rest. The
    # objects are in sha256 order, the ledger's sort key; rows stream to the writer
    tools = provenance.tool_versions("pyarrow")
    ledger.write(out, stage, "ingest", (ledger.row(sha, "l0_object", n_in[i], kept[i], n_in[i] - kept[i], tools=tools)
                                        for i, (sha, _) in enumerate(ordered)))
    hook("ledger")

    # l0_manifest: small metadata, rewritten with all stages; written last and renamed into place, so
    # a stage that died earlier has not touched it (card B-23)
    man_spec = table("l0_manifest")
    man_path = root / "l0/l0_manifest.parquet"
    man = pa.Table.from_pylist([{
        "source_uid": ids.source_uid(obj.source_name, obj.source_release),
        "source_name": obj.source_name,
        "source_release": obj.source_release,
        "object_path": Path(obj.path).name,
        "object_sha256": sha,
        "bytes": Path(obj.path).stat().st_size,
        "license": obj.license,
        "redistributable": obj.redistributable,
        "retrieved_at": obj.retrieved_at,
    } for sha, obj in ordered], schema=man_spec.schema)
    if man_path.exists():
        man = pa.concat_tables([pq.read_table(man_path), man])
    man = man.sort_by([("source_uid", "ascending"), ("object_path", "ascending")])
    man_out = out / "l0/l0_manifest.parquet"
    man_tmp = man_out.with_name(f".{man_out.name}.tmp-{os.getpid()}")
    _write(man, man_tmp, 122_880)
    os.replace(man_tmp, man_out)
    hook("l0_manifest")
    return {"objects": len(ordered), "records": sum(n_in), "new_sequences": n_new, "index_rows": n_idx,
            "shards": len(batches)}


class SeqStore:
    """Reader: slices by (contig_uid, start, end, strand), 0-based half-open.

    `seq_locator` maps a raw_digest to its shard; the shard's `.fai` and `.gzi` find the
    blocks. `counters()` sums blocks decompressed and bytes read over all shards (LY-06).

    Nothing is loaded on open but the Parquet footers of the two indexes: `locator` (by
    `raw_digest`) and `index` (`contig_index` by `contig_uid`) are `KeyLookup`s that read one
    row group per key, with an LRU of `cache_row_groups` decoded row groups each. An unknown
    key raises KeyError.
    """

    def __init__(self, root: Path, cache_blocks: int = 64, cache_row_groups: int = 64):
        self.root = Path(root)
        self.locator_files = locator_files(self.root)
        if not self.locator_files:
            raise FileNotFoundError(f"{self.root / LOCATOR_DIR}: no seq_locator files")
        self.index_files = index_files(self.root)
        self.locator = KeyLookup(self.locator_files, "raw_digest", ["shard", "contig_len"], cache_row_groups)
        self.index = KeyLookup(self.index_files, "contig_uid", ["raw_digest", "is_default_variant"], cache_row_groups)
        self._cache_blocks = cache_blocks
        self._shards: dict[str, BgzfFasta] = {}

    def _where(self, rd: bytes) -> tuple[str, int]:
        rows = self.locator.get(rd)
        if not rows:
            raise KeyError(rd)
        return rows[0]["shard"], rows[0]["contig_len"]

    def default_digest(self, contig_uid: bytes) -> bytes:
        """raw_digest of the contig's default variant."""
        for r in self.index.get(contig_uid):
            if r["is_default_variant"]:
                return r["raw_digest"]
        raise KeyError(contig_uid)

    def _shard(self, rel: str) -> BgzfFasta:
        if rel not in self._shards:
            self._shards[rel] = BgzfFasta(self.root / rel, cache_blocks=self._cache_blocks)
        return self._shards[rel]

    def shard_of(self, rd: bytes) -> str:
        return self._where(rd)[0]

    def counters(self) -> dict[str, int]:
        s = self._shards.values()
        return {"blocks_read": sum(b.blocks_read for b in s),
                "compressed_bytes_read": sum(b.compressed_bytes_read for b in s),
                "uncompressed_bytes_decoded": sum(b.uncompressed_bytes_decoded for b in s)}

    def close(self) -> None:
        for b in self._shards.values():
            b.close()

    def length_of(self, rd: bytes) -> int:
        return self._where(rd)[1]

    def locate(self, rd: bytes) -> tuple[str, int]:
        """(shard path under the root, contig_len) of a stored variant, in one locator lookup."""
        return self._where(rd)

    def read(self, shard: str, rd: bytes, start: int, end: int) -> bytes:
        """Forward bytes of a variant already located (`locate`), with no index lookup; batch
        readers (`catalog.bb_seq`) use it to read in shard order."""
        return self._shard(shard).fetch(rd.hex(), start, end)

    def raw(self, rd: bytes, start: int = 0, end: int | None = None) -> bytes:
        """Forward-strand bytes of one stored variant."""
        rel, length = self._where(rd)
        end = length if end is None else end
        if not 0 <= start <= end <= length:
            raise ValueError(f"slice {start}:{end} outside 0:{length}")
        return self._shard(rel).fetch(rd.hex(), start, end)

    def get(
        self, contig_uid: bytes, start: int = 0, end: int | None = None, strand: int = 1, circular: bool = False
    ) -> bytes:
        """Slice the default variant. On a circular contig `end` may exceed the length (origin wrap)."""
        rd = self.default_digest(contig_uid)
        n = self.length_of(rd)
        if circular and end is not None and end > n:
            if not 0 <= start < n < end <= start + n:
                raise ValueError(f"wrapping slice {start}:{end} invalid for length {n}")
            seq = self.raw(rd, start, n) + self.raw(rd, 0, end - n)
        else:
            seq = self.raw(rd, start, end)
        if strand == -1:
            return revcomp(seq)
        if strand != 1:
            raise ValueError("strand must be +1 or -1")
        return seq

    def reconstruct_object(self, l0_object_sha256: bytes) -> list[bytes]:
        """FASTA lines of one source object, from occurrence rows and chunks (R-02)."""
        occ = pa.concat_tables(
            pq.read_table(f) for f in io.part_files(self.root / "l0seq/contig_occurrence", "stage=*/source=*/part-*.parquet")
        )
        rows = sorted(
            (r for r in occ.to_pylist() if r["l0_object_sha256"] == l0_object_sha256),
            key=lambda r: r["record_ordinal"],
        )
        lines: list[bytes] = []
        for r in rows:
            lines.extend(fasta.format_record(r["original_header"], self.raw(r["raw_digest"]), r["line_width"]))
        return lines
