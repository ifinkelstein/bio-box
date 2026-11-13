"""Write and read the L0-seq store (operations O-INGEST, O-SLICE, O-EXTRACT).

Layout under `<root>/` (planning/data_layer_and_tests.md §2A.1, §3.1):

    l0/l0_manifest/stage=<S>/part-<k>.parquet
    l0seq/store_params.json
    l0seq/seq_shards/stage=<S>/shard-<NNNN>.fa.gz   (+ .fa.gz.fai, .fa.gz.gzi)
    l0seq/contig_occurrence/stage=<S>/source=<name>/part-<k>.parquet
    l0seq/contig_label/stage=<S>/part-<k>.parquet
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

`contig_label` (card B-27) is `contig_occurrence` by contig, for readers that look a contig's
accession and header up by `contig_uid` (`contig_occurrence` is sorted by object and ordinal, so
such a lookup cannot prune it): one row per contig of the stage, sorted by `contig_uid`, with its
first record by (`l0_object_sha256`, `record_ordinal`) and its number of records. It is made from
one DuckDB sort of the stage's records by (`contig_uid`, `l0_object_sha256`, `record_ordinal`),
which spills into the private scratch directory, read a batch at a time (`_label_batches`).

A stage never rewrites an earlier stage's files (AGENTS §3.10): a sequence already stored
adds only occurrence rows, a stage's `seq_locator` file holds only the sequences first
stored in that stage, and its `l0_manifest` partition lists only the objects it ingested
(card B-22). `store_params.json` (the `store_params` artifact) is written once, by the first
stage; every later stage checks it and writes none.

A stage too large for one process (AGENTS §7.2) is ingested by batch processes (card B-30): pass-1 shards,
a plan, one process per batch (one shard each, DEC-21) and a merge, whose files are `ingest`'s byte for byte;
their records are the `ingest_plan` artifact, `l0seq/ingest_plan/stage=<S>/` (section below `_write_tables`).

`SeqStore` loads neither index (card B-13): `seq_locator` and `contig_index` are looked up
through `KeyLookup` (footers once, one row group per key) and each shard's `.fai` is bisected
on disk (`BgzfFasta`). It reads only the stages it is opened for (card B-28): `SeqStore(root, ["S0"])`,
or `SeqStore(root, catalog.ReleaseFiles.of(root, "rS0.1"))` for a release's files.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import tempfile
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from biobox import fasta, ids, io, layout, ledger, provenance
from biobox.l0seq.bgzf import BgzfFasta, write_bgzf_fasta
from biobox.l0seq.lookup import KeyLookup
from biobox.registry.schemas import ARTIFACTS, table

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
# contig_label: 47.9 bytes per row, measured on wilkcomp02, 2026-11-10, rL.3's three stages (32,815 contigs;
# planning/measurements/2026-11-10_label_and_membership_lookups.md); it only cuts parts
EST_LABEL_ROW = 48
EST_INDEX_ROW = 56  # estimate: contig_index (contig_uid 16 + raw_digest 32 + length and flag)
EST_LOCATOR_ROW = 40  # estimate: seq_locator (raw_digest 32 + dictionary-coded shard path + length)
EST_MANIFEST_ROW = 160  # estimate: l0_manifest (sha256 32 + object name, source and license strings + 3 scalars)

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


LOCATOR_DIR = "l0seq/seq_locator"
INDEX_DIR = "l0seq/contig_index"
LABEL_DIR = table("contig_label").location
MANIFEST_DIR = table("l0_manifest").location
PARAMS_PATH = ARTIFACTS["store_params"].location


def locator_files(root: Path) -> list[Path]:
    """Every stage's `seq_locator` files on disk, by stage directory and then part order. Only ingest
    reads these (a new stage's anti-join must see every sequence already stored); a reader takes its
    stages' files (`SeqStore`, card B-28)."""
    return io.part_files(Path(root) / LOCATOR_DIR, "stage=*/part-*.parquet")


def index_files(root: Path) -> list[Path]:
    """Every stage's `contig_index` files on disk; ingest only, as `locator_files`."""
    return io.part_files(Path(root) / INDEX_DIR, "stage=*/part-*.parquet")


# ledger reason of a contig that another stage owns (DEC-26), followed by that stage's name, or that another genome of
# the same stage owns (DEC-27), followed by `<stage>/<genome_uid hex>` (`ContigOwners.owner`)
OWNED_BY = "owned_by:"


def stages_on_disk(root: Path) -> list[str]:
    """The stages whose `seq_locator` directory exists in the root, sorted by name: the stages a new stage's ingest
    compares its sequences against (`locator_files`). `build-stage` passes them to gene calling (card C-20), whose
    genomes may name a contig that one of them stored first; nothing is read from a stage that holds none."""
    base = Path(root) / LOCATOR_DIR
    return sorted(p.name.removeprefix("stage=") for p in base.glob("stage=*") if p.is_dir()) if base.is_dir() else []


def stored_by(root: Path, contig_uids, skip=()) -> dict[bytes, str]:
    """For each of `contig_uids`, the stage on disk (not in `skip`) whose `contig_index` holds its default variant,
    where one does. For error messages that name a contig's owning stage (DEC-26); a reader never uses it to choose
    what to read (AGENTS §3.9)."""
    out: dict[bytes, str] = {}
    want = list(dict.fromkeys(contig_uids))
    for s in stages_on_disk(root):
        if s in skip:
            continue
        look = KeyLookup(io.part_files(Path(root) / INDEX_DIR / f"stage={s}"), "contig_uid", ["is_default_variant"])
        for cu in want:
            if cu not in out and any(r["is_default_variant"] for r in look.get(cu)):
                out[cu] = s
    return out


class ContigOwners:
    """Which stage owns each contig (DEC-26, card C-20): the stage whose `contig_index` holds the contig's default
    variant. Ingest writes a sequence once, in the first stage that sees it, and a later stage's `contig_index`
    holds no default variant of a contig an earlier stage stored, so that stage is unique; it owns the contig's
    gene calls, IGRs and tape. `stages` are the stages to look in, `stage` (the caller's own) first; one key
    lookup per stage at most. A contig owned by none of them raises KeyError naming the stage on disk that owns
    it, if any.

    Within the caller's stage, a contig in several of its genomes belongs to the genome with the least
    `genome_uid` (DEC-27, card C-21a): `shared` maps each such contig to that genome (`inputs.shared_contigs`)."""

    def __init__(self, root: Path, stage: str, stages, shared: dict[bytes, bytes] | None = None):
        self.root, self.stage = Path(root), stage
        self.stages = list(dict.fromkeys([stage, *stages]))
        self.shared = {} if shared is None else dict(shared)
        self._index = {s: KeyLookup(io.part_files(self.root / INDEX_DIR / f"stage={s}"), "contig_uid",
                                    ["is_default_variant"]) for s in self.stages}

    def owner(self, contig_uid: bytes, genome_uid: bytes | None = None) -> str:
        """The owning stage; with `genome_uid`, a contig of the caller's stage that another of its genomes owns is
        `<stage>/<that genome_uid in hex>` (`OWNER_SEP`). Either way, the owner is the caller exactly when the
        result is `self.stage`, and the result is what follows `OWNED_BY` in a skip's ledger reason."""
        for s in self.stages:
            if any(r["is_default_variant"] for r in self._index[s].get(contig_uid)):
                g = self.shared.get(contig_uid) if s == self.stage and genome_uid is not None else None
                return s if g is None or g == genome_uid else f"{s}{OWNER_SEP}{g.hex()}"
        raise KeyError(missing_contig(self.root, contig_uid, self.stages))


OWNER_SEP = "/"  # between the stage and the genome in an owner within a stage (`ContigOwners.owner`)


def skipped_entry(contig_uid: bytes, owner: str) -> dict:
    """A producer's result entry for a contig it skipped (C-20): the contig, the owning stage and, for a contig that
    another genome of the producer's own stage owns (C-21a), that genome."""
    stage, _, genome = owner.partition(OWNER_SEP)
    return {"contig_uid": contig_uid.hex(), "owner_stage": stage, **({"owner_genome": genome} if genome else {})}


def missing_contig(root: Path, contig_uid: bytes, stages) -> str:
    """The message for a contig that none of `stages` holds: the stage on disk that owns it, if one does."""
    found = stored_by(root, [contig_uid], skip=set(stages)).get(contig_uid)
    if found is None:
        return f"contig_uid {contig_uid.hex()} is in none of the stages {list(stages)} read, nor in any stage on disk"
    return (f"contig_uid {contig_uid.hex()} is in none of the stages {list(stages)} read: stage {found} stored it "
            f"first and owns it (DEC-26); read that stage too (seq_stages, --seq-stages)")


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


def _order_key(t: tuple[bytes, SourceObject]) -> tuple:
    return (t[0], t[1].source_name, t[1].source_release, t[1].precedence, Path(t[1].path).name, str(t[1].path))


def _object_order(objects: list[SourceObject]) -> list[tuple[bytes, SourceObject]]:
    """(sha256, object) in ascending sha256; ties (the same bytes listed twice) by the object's fields."""
    keyed = [(ids.sha256_file(o.path), o) for o in objects]
    return sorted(keyed, key=_order_key)


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


def _wanted(con, batch_rows: int = io.ROW_GROUP_ROWS, objects: tuple[int, int] | None = None):
    """(obj_idx, record_ordinal, raw_digest, contig_len) of every first occurrence, in object then
    record order; with `objects` (first, last), of the objects in that range only (one batch, card B-30)."""
    where = "" if objects is None else f"WHERE obj_idx BETWEEN {int(objects[0])} AND {int(objects[1])} "
    reader = con.sql("SELECT obj_idx, record_ordinal, raw_digest, contig_len FROM new_seq "
                     f"{where}ORDER BY obj_idx, record_ordinal").to_arrow_reader(batch_rows)
    try:
        for b in reader:
            yield from zip(*(b.column(i).to_pylist() for i in range(4)))
    finally:
        reader.close()


_LABEL_COLUMNS = ("contig_uid", "accession", "original_header", "l0_object_sha256", "record_ordinal")


def _label_batches(con, spec, batch_rows: int = io.ROW_GROUP_ROWS):
    """The stage's `contig_label` rows in contig_uid order (card B-27): from one DuckDB sort of every record of
    `seen` by (contig_uid, l0_object_sha256, record_ordinal), which spills, each contig's first row with
    `n_records`, its number of rows. The sorted rows come `batch_rows` at a time; a contig is cut where the
    contig_uid changes, and the last contig of a batch is held (one row and a count) until the next batch shows
    whether it goes on. Ties (one object listed twice) are broken by `obj_idx`; their rows are equal."""
    reader = con.sql(f"SELECT {', '.join(_LABEL_COLUMNS)} FROM seen "
                     "ORDER BY contig_uid, l0_object_sha256, record_ordinal, obj_idx").to_arrow_reader(batch_rows)
    held: tuple[pa.RecordBatch, int] | None = None  # the last contig seen: its first row and its rows so far

    def labelled(rows: pa.RecordBatch, counts) -> pa.RecordBatch:
        n = pa.array(np.asarray(counts, dtype=np.uint32), pa.uint32())
        return io.conform(pa.RecordBatch.from_arrays([*rows.columns, n], names=[*_LABEL_COLUMNS, "n_records"]), spec)

    try:
        for b in reader:
            n = b.num_rows
            if not n:
                continue
            uids = b.column(0)
            cut = np.flatnonzero(pc.not_equal(uids.slice(1), uids.slice(0, n - 1)).to_numpy(zero_copy_only=False)) + 1
            starts = np.concatenate([[0], cut]).astype(np.int64)
            counts = np.diff(np.concatenate([starts, [n]]))
            first = 0
            if held is not None:
                if held[0].column(0)[0].as_py() == uids[0].as_py():  # the held contig goes on into this batch
                    held = (held[0], held[1] + int(counts[0]))
                    first = 1
                if first == 0 or len(starts) > 1:
                    yield labelled(held[0], [held[1]])
                    held = None
            if held is None:
                done = starts[first:-1]
                if len(done):
                    yield labelled(b.take(pa.array(done)), counts[first:-1])
                held = (b.slice(int(starts[-1]), 1), int(counts[-1]))
        if held is not None:
            yield labelled(held[0], [held[1]])
    finally:
        reader.close()


def _stat(path: Path) -> tuple[int, int]:
    st = Path(path).stat()
    return st.st_size, st.st_mtime_ns


def stage_dirs(stage: str) -> list[str]:
    """The root-relative directories that hold one stage's ingest files (the ledger's included)."""
    return [f"l0seq/seq_shards/stage={stage}", f"l0seq/contig_occurrence/stage={stage}", f"{LABEL_DIR}/stage={stage}",
            f"{INDEX_DIR}/stage={stage}", f"{LOCATOR_DIR}/stage={stage}", f"ledger/stage={stage}/step=ingest",
            f"{MANIFEST_DIR}/stage={stage}"]


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
    under `root` or `out`. After a direct ingest died, remove those directories and run it again.
    Nothing here overwrites a file: `io.write_table` still refuses to.

    Store parameters and the object list (card B-22). `store_params.json` (the `store_params`
    artifact) is checked before any source is read: when the root has one, the parameters must equal
    it (ValueError otherwise) and nothing is written; when it has none, this is the store's first
    stage and the file is written under `out`, after the stage's tables, so run-step publishes it with
    them. The stage's `l0_manifest` partition (`l0/l0_manifest/stage=<S>/`) lists the objects of this
    stage only and is written last; no earlier stage's row is read or rewritten.

    `_hook` (tests only) is called with each point at which a kill leaves the stage half written:
    "pass1", "shard:<j>" after each shard, "contig_occurrence", "contig_label", "contig_index",
    "seq_locator", "ledger", "l0_manifest"; an exception it raises is not caught.
    """
    root = Path(root)
    out = root if out is None else Path(out)
    hook = _hook or _no_hook
    if tmp is not None:
        for base, what in ((root, "the data root"), (out, "the output directory")):
            if Path(tmp).resolve().is_relative_to(base.resolve()):
                raise ValueError(f"{tmp}: the scratch directory may not lie under {what}")
    _refuse_stage_written(root, out, stage, stage_dirs(stage))
    params = _params(shard_bases, line_width)
    new_params = _check_params(root, out, params)

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
            report = _decide_and_write(con, root, out, stage, ordered, n_in, stats, seen_path, old_loc, old_idx,
                                       scratch, shard_bases, line_width, max_file_bytes, hook)
        finally:
            con.close()
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    # the store's parameters, by its first stage only, after every table of the stage; then the
    # stage's l0_manifest partition, last (card B-22)
    if new_params is not None:
        _write_new(out / PARAMS_PATH, new_params)
    _write_l0_manifest(out, stage, ordered, stats, max_file_bytes)
    hook("l0_manifest")
    return report


def _refuse_stage_written(root: Path, out: Path, stage: str, dirs: list[str]) -> None:
    """FileExistsError when any file (a temporary one included) is in `dirs` of the stage under `root` or `out`."""
    for base in dict.fromkeys([Path(root).resolve(), Path(out).resolve()]):
        left = sorted(p for d in dirs if (base / d).is_dir() for p in (base / d).rglob("*") if not p.is_dir())
        if left:
            raise FileExistsError(
                f"{left[0]}: stage {stage} already has {len(left)} ingest file(s) under {base}; ingest "
                f"writes a stage once. A direct call needs a clean stage directory: if an earlier ingest "
                f"died, remove {', '.join(stage_dirs(stage))} and run it again, or ingest under run-step")


def _params(shard_bases: int, line_width: int) -> dict:
    return {"format": STORE_FORMAT, "shard_bases": shard_bases, "line_width": line_width, "bgzf_level": BGZF_LEVEL}


def _params_bytes(params: dict) -> bytes:
    return json.dumps(params, sort_keys=True).encode()


def _check_params(root: Path, out: Path, params: dict) -> bytes | None:
    """The `store_params` artifact (card B-22), checked before any source is read. The root has it:
    the parameters must equal it (ValueError otherwise) and None is returned, since a later stage
    writes none. The root has none: this is the store's first stage, and the bytes to write under
    `out` are returned (FileExistsError if `out` already holds a copy)."""
    have, mine = root / PARAMS_PATH, out / PARAMS_PATH
    if have.exists():
        stored = json.loads(have.read_text())
        if stored != params:
            raise ValueError(f"{have}: the store's parameters {stored} differ from this ingest's {params}; "
                             "store_params is written once, by a store's first stage, and every later stage "
                             "must use the same parameters")
        return None
    if mine.exists():
        raise FileExistsError(f"{mine}: already exists; store_params is written once, by the store's first stage")
    return _params_bytes(params)


def _write_new(path: Path, data: bytes) -> None:
    """`data` at `path`, which must not exist, under a temporary name renamed into place (AGENTS §3.9)."""
    if path.exists():
        raise FileExistsError(f"{path}: already exists; refusing to overwrite")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def _write_l0_manifest(out: Path, stage: str, ordered: list[tuple[bytes, SourceObject]],
                       stats: list[tuple[int, int]], max_file_bytes: int) -> None:
    """The stage's own `l0_manifest` partition: one row per object it ingested, sorted by
    (source_uid, object_path), ties by sha256 order; written once and never rewritten (AGENTS §3.10,
    card B-22). Its rows are as many as the objects ingest was handed."""
    spec = table("l0_manifest")
    src = [ids.source_uid(o.source_name, o.source_release) for _, o in ordered]
    name = [Path(o.path).name for _, o in ordered]
    order = sorted(range(len(ordered)), key=lambda i: (src[i], name[i], i))
    objs = [ordered[i][1] for i in order]
    tbl = pa.table({
        "source_uid": pa.array([src[i] for i in order], pa.uint64()),
        "source_name": pa.array([o.source_name for o in objs], pa.string()),
        "source_release": pa.array([o.source_release for o in objs], pa.string()),
        "object_path": pa.array([name[i] for i in order], pa.string()),
        "object_sha256": pa.array([ordered[i][0] for i in order], pa.binary(32)),
        "bytes": pa.array([stats[i][0] for i in order], pa.int64()),  # the size pass 2 checked
        "license": pa.array([o.license for o in objs], pa.string()),
        "redistributable": pa.array([o.redistributable for o in objs], pa.bool_()),
        "retrieved_at": pa.array([o.retrieved_at for o in objs], pa.timestamp("ms", tz="UTC")),
    }, schema=spec.schema)
    io.write_parts(tbl, spec, out / MANIFEST_DIR / f"stage={stage}", EST_MANIFEST_ROW, max_file_bytes)


def _decide_and_write(con, root: Path, out: Path, stage: str, ordered: list[tuple[bytes, SourceObject]],
                      n_in: list[int], stats: list[tuple[int, int]], seen_path: Path, old_loc: list[Path],
                      old_idx: list[Path], scratch: Path, shard_bases: int, line_width: int, max_file_bytes: int,
                      hook: Callable[[str], None]) -> dict:
    con.execute(f"CREATE VIEW seen AS SELECT * FROM read_parquet({_sql_path(seen_path)})")
    _old_views(con, old_loc, old_idx)
    kept, batches, _ = _decide(con, scratch / "new_seq.parquet", len(ordered), shard_bases)
    shard_of = _shard_paths(stage, len(batches))

    # pass 2: one batch's first occurrences at a time, one shard per batch
    def unchanged(i: int, path: Path) -> None:
        if _stat(path) != stats[i]:
            raise ValueError(f"{path}: changed since the first pass (size or mtime)")

    wanted = _wanted(con)
    nxt = next(wanted, None)
    try:
        for j, batch in enumerate(batches):
            nxt = _write_shard(out / shard_of[j], [(i, Path(ordered[i][1].path)) for i in batch], wanted, nxt,
                               unchanged, line_width)
            hook(f"shard:{j}")
    finally:
        wanted.close()
    return _write_tables(con, out, stage, ordered, n_in, kept, batches, max_file_bytes, hook)


def _old_views(con, old_loc: list[Path], old_idx: list[Path]) -> None:
    """Views `old_loc` and `old_idx` over the earlier stages' `seq_locator` and `contig_index` files."""
    con.execute("CREATE VIEW old_loc AS " + (
        f"SELECT raw_digest FROM read_parquet({_sql_files(old_loc)}, union_by_name = true)" if old_loc
        else "SELECT NULL::BLOB AS raw_digest WHERE false"))
    con.execute("CREATE VIEW old_idx AS " + (
        f"SELECT contig_uid, raw_digest, is_default_variant FROM read_parquet({_sql_files(old_idx)}, "
        "union_by_name = true)" if old_idx
        else "SELECT NULL::BLOB AS contig_uid, NULL::BLOB AS raw_digest, NULL::BOOLEAN AS is_default_variant WHERE false"))


def _shard_paths(stage: str, n: int) -> list[str]:
    return [f"l0seq/seq_shards/stage={stage}/shard-{j:04d}.fa.gz" for j in range(n)]


def _decide(con, new_path: Path, n_objects: int, shard_bases: int) -> tuple[list[int], list[list[int]], list[int]]:
    """Over the views `seen` and `old_loc`: write the stage's new sequences to `new_path` (the view `new_seq`)
    and cut the objects into batches (`_batches`). Returns (sequences kept per object, the batches, new bases
    per object)."""
    # new digests: not in any earlier stage, each at its first occurrence by (object sha256, ordinal)
    con.execute(f"""
        COPY (SELECT raw_digest, contig_len, obj_idx, record_ordinal FROM (
                SELECT s.raw_digest, s.contig_len, s.obj_idx, s.record_ordinal,
                       row_number() OVER (PARTITION BY s.raw_digest ORDER BY s.obj_idx, s.record_ordinal) AS rn
                FROM seen s ANTI JOIN old_loc l ON s.raw_digest = l.raw_digest)
              WHERE rn = 1 ORDER BY obj_idx, record_ordinal)
        TO {_sql_path(new_path)} (FORMAT parquet)""")
    con.execute(f"CREATE VIEW new_seq AS SELECT * FROM read_parquet({_sql_path(new_path)})")
    kept, new_bases = [0] * n_objects, [0] * n_objects
    per_object = []
    for i, n, bases in con.execute("SELECT obj_idx, count(*), sum(contig_len) FROM new_seq "
                                   "GROUP BY obj_idx ORDER BY obj_idx").fetchall():
        kept[i], new_bases[i] = n, int(bases)
        per_object.append((i, int(bases)))
    return kept, _batches(per_object, shard_bases), new_bases


def _write_shard(path: Path, members: list[tuple[int, Path]], wanted, nxt, check: Callable[[int, Path], None],
                 line_width: int):
    """Pass 2 for one batch: read its objects (`members`, (obj_idx, path) in order; `check` is called on each
    before it is read) for their first occurrences, which `wanted` yields in object then record order (`nxt` is
    the next one, read already), and write them as one shard sorted by `raw_digest`. Only this batch's
    sequences are held. Returns the next first occurrence after the batch's (None at the end)."""
    seqs: dict[bytes, bytes] = {}
    for i, src in members:
        check(i, src)
        for k, rec in enumerate(fasta.iter_records(src)):
            if nxt is None or nxt[0] != i:
                break  # nothing more to keep from this object
            if nxt[1] == k:
                # stored bytes are identified by their digest (AGENTS §3.2): hash again (about
                # 3 s per 240 s at 1,000x, card B-23), after the cheap length check
                if len(rec.seq) != nxt[3] or ids.raw_digest(rec.seq) != nxt[2]:
                    raise ValueError(f"{src}: record {k} changed since the first pass")
                seqs[nxt[2]] = rec.seq
                nxt = next(wanted, None)
        rec = None
        if nxt is not None and nxt[0] == i:
            raise ValueError(f"{src}: fewer records than in the first pass")
    write_bgzf_fasta(((rd.hex(), seqs[rd]) for rd in sorted(seqs)), path, line_width, BGZF_LEVEL)
    del seqs
    return nxt


def _write_tables(con, out: Path, stage: str, ordered: list[tuple[bytes, SourceObject]], n_in: list[int],
                  kept: list[int], batches: list[list[int]], max_file_bytes: int, hook: Callable[[str], None]) -> dict:
    """The stage's tables but `l0_manifest`, over the views `seen`, `old_idx` and `new_seq`, from one connection:
    `contig_occurrence`, `contig_label`, `contig_index`, `seq_locator` (each sequence's shard is its batch's) and
    the ledger, in that order, `hook` called after each. Returns ingest's counts."""
    shard_of = _shard_paths(stage, len(batches))
    con.register("batch_of", pa.table({
        "obj_idx": pa.array([i for b in batches for i in b], pa.uint32()),
        "shard": pa.array([shard_of[j] for j, b in enumerate(batches) for _ in b], pa.string())}))
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

    # contig_label: contig_occurrence by contig (card B-27), every source of the stage together; a stage without
    # records has neither table
    (n_records,), = con.execute("SELECT count(*) FROM seen").fetchall()
    if n_records:
        lab_spec = table("contig_label")
        io.write_parts(_label_batches(con, lab_spec), lab_spec, out / f"{LABEL_DIR}/stage={stage}", EST_LABEL_ROW,
                       max_file_bytes)
    hook("contig_label")

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
    return {"objects": len(ordered), "records": sum(n_in), "new_sequences": n_new, "index_rows": n_idx,
            "shards": len(batches)}


# --- ingest by batch processes (card B-30) -----------------------------------------------------------------
#
# The single-process `ingest` above is cut into processes, each under the CPU-minute budget (AGENTS §7.2), that
# give its files byte for byte:
#
# - `ingest_pass1`, one process per pass-1 shard: the objects in `pass1_order` (by name and path: no read is
#   needed to sort them) cut into runs of at most `pass1_bytes` bytes of file (`pass1_shards`; an object above
#   that is a shard of its own). Each streams its objects' records into a scratch run of the rows pass 1 writes
#   (`_SEEN`, with the object's pass-1 index for `obj_idx`) and publishes its record `pass1-<k>.json` (the
#   objects' names, sizes, sha256 and record counts, and the run's sha256).
# - `ingest_plan`: one DuckDB connection over the pass-1 runs, each object numbered by its rank in sha256 order
#   as `ingest` numbers it, makes `ingest`'s decisions (`_decide`): the new sequences (the scratch decision
#   table `new_seq.parquet`) and the batches, a pure function of the sorted inputs and `shard_bases` (DEC-21).
#   It publishes `plan.json`: the objects in sha256 order, each batch's objects and shard, the earlier stages'
#   index files it decided against, and the decision table's sha256.
# - `ingest_batch`, one process per batch: its objects read again (each checked against the plan's size and
#   sha256) for its new sequences, written as its shard (`_write_shard`, as `ingest` writes it).
# - `ingest_merge`, under the step name `ingest`: the tables, ledger, `store_params` and `l0_manifest`
#   (`_write_tables`, as `ingest` writes them) over the same views, once every batch's shard is published.
#
# The pass-1 cut decides only which process reads which object: no output depends on it, nor on the order of
# the objects given. The scratch runs live in `ingest_runs_dir` (under `tmp`, never in the root); when scratch
# is lost after a record was published (a reboot clears /tmp), `rescan_pass1` and `rescan_plan` make a run
# again, and it must have the sha256 its record holds (`ingest_run_problems` names what is missing).

PLAN_DIR = ARTIFACTS["ingest_plan"].location
PLAN_FORMAT = 1
RUNS_DIR = "biobox-ingest"  # under the scratch directory: one directory per root and stage
NEW_SEQ = "new_seq.parquet"  # the plan's decision table, in the stage's runs directory


def plan_dir(stage: str) -> str:
    """The root-relative directory of the stage's `ingest_plan` records."""
    return f"{PLAN_DIR}/stage={stage}"


def _pass1_rel(stage: str, k: int) -> str:
    return f"{plan_dir(stage)}/pass1-{k:05d}.json"


def _plan_rel(stage: str) -> str:
    return f"{plan_dir(stage)}/plan.json"


def _table_dirs(stage: str) -> list[str]:
    """The stage's ingest directories but the shards': what only the merge writes."""
    return [d for d in stage_dirs(stage) if not d.startswith("l0seq/seq_shards/")]


def _pass1_key(o: SourceObject) -> tuple:
    return (Path(o.path).name, str(o.path), o.source_name, o.source_release, o.precedence, o.license,
            o.redistributable, o.retrieved_at.isoformat())


def pass1_order(objects: list[SourceObject]) -> list[SourceObject]:
    """The objects in pass-1 order: by file name, then path and the other fields. Nothing is read to sort them."""
    return sorted(objects, key=_pass1_key)


def pass1_shards(objects: list[SourceObject], pass1_bytes: int) -> list[list[int]]:
    """Indices into `pass1_order(objects)`, cut into maximal runs of at most `pass1_bytes` bytes of file (the size
    on disk); an object above that is a shard of its own. A stage without objects is one empty shard."""
    if isinstance(pass1_bytes, bool) or not isinstance(pass1_bytes, int) or pass1_bytes < 1:
        raise ValueError(f"pass1_bytes {pass1_bytes!r}: must be a positive integer")
    order = pass1_order(objects)
    return _batches([(i, Path(o.path).stat().st_size) for i, o in enumerate(order)], pass1_bytes) or [[]]


def _runs_base(tmp: Path | None, runs: Path | None) -> Path | None:
    """Where the batch steps' shared runs live: `runs` when given, else `tmp` (one directory for both)."""
    return tmp if runs is None else runs


def ingest_runs_dir(tmp: Path | None, root: Path, stage: str) -> Path:
    """The scratch directory of the stage's pass-1 runs and decision table: under `tmp` (else the system's),
    named by the root and the stage, so two roots or stages never mix their runs."""
    key = hashlib.sha256(f"{Path(root).resolve()}\n{stage}".encode()).hexdigest()[:16]
    return Path(tmp if tmp is not None else tempfile.gettempdir()) / RUNS_DIR / f"{stage}-{key}"


def _json_bytes(doc: dict) -> bytes:
    return json.dumps(doc, sort_keys=True, indent=1).encode() + b"\n"


def _find(root: Path, out: Path | None, rel: str) -> Path | None:
    """`rel` in the root, else in `out` (written by a direct call into another output root)."""
    for base in [Path(root)] + ([] if out is None else [Path(out)]):
        if (base / rel).is_file():
            return base / rel
    return None


def _doc(root: Path, out: Path | None, rel: str, what: str) -> tuple[bytes, dict]:
    p = _find(root, out, rel)
    if p is None:
        raise FileNotFoundError(f"{Path(root) / rel}: no {what}")
    raw = p.read_bytes()
    return raw, json.loads(raw)


def _write_once(root: Path, out: Path, rel: str, data: bytes, what: str) -> None:
    for base in dict.fromkeys([Path(root).resolve(), Path(out).resolve()]):
        if (base / rel).exists():
            raise FileExistsError(f"{base / rel}: {what} is written once; refusing to overwrite")
    _write_new(Path(out) / rel, data)


def read_pass1(root: Path, stage: str, k: int, out: Path | None = None) -> dict:
    """The published record of pass-1 shard k (from the root, else from `out`)."""
    return _doc(root, out, _pass1_rel(stage, k), f"record of pass-1 shard {k} of stage {stage}")[1]


def read_ingest_plan(root: Path, stage: str, out: Path | None = None) -> dict:
    """The stage's published ingest plan (from the root, else from `out`)."""
    return _doc(root, out, _plan_rel(stage), f"ingest plan of stage {stage}; write it with ingest_plan")[1]


def _hint(stage: str, what: str) -> str:
    return (f"its record is published, so make it again into this scratch directory: python -m biobox.l0seq.store "
            f"rescan --stage {stage} {what} (build-stage rescans lost runs itself)")


def _run_why(path: Path, sha256: str) -> str | None:
    if not path.is_file():
        return "is missing"
    if ids.sha256_file(path).hex() != sha256:
        return "differs from the sha256 its record holds"
    return None


def _checked_run(path: Path, sha256: str, what: str, stage: str, how: str) -> Path:
    """`path`, a scratch run whose record holds `sha256`; FileNotFoundError or ValueError naming the rescan."""
    why = _run_why(path, sha256)
    if why == "is missing":
        raise FileNotFoundError(f"{path}: the scratch run of {what} of stage {stage} is missing; {_hint(stage, how)}")
    if why is not None:
        raise ValueError(f"{path}: the scratch run of {what} of stage {stage} {why}; {_hint(stage, how)}")
    return path


def _duck(scratch: Path):
    return duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                  "temp_directory": str(Path(scratch) / "duckdb"), "preserve_insertion_order": True})


def _pass1_run(order: list[SourceObject], members: list[int], path: Path) -> list[dict]:
    """Pass 1 over the objects `members` (indices into `order`): their records into the scratch run `path`, as
    `ingest` writes them, `obj_idx` being the pass-1 index. Returns each object's name, size, sha256, records and
    bases."""
    seen = _Scratch(path, _SEEN)
    docs = []
    try:
        for ci in members:
            obj = order[ci]
            src = Path(obj.path)
            size, sha = src.stat().st_size, ids.sha256_file(src)
            src_uid = ids.source_uid(obj.source_name, obj.source_release)
            n = bases = 0
            for k, rec in enumerate(fasta.iter_records(src)):
                seen.add((obj.source_name, src_uid, obj.source_release, sha, k, rec.accession, rec.header,
                          ids.raw_digest(rec.seq), ids.contig_uid(rec.seq), rec.line_width, len(rec.seq),
                          obj.precedence, ci))
                n, bases = k + 1, bases + len(rec.seq)
            rec = None
            docs.append({"index": ci, "name": src.name, "bytes": size, "sha256": sha.hex(), "records": n,
                         "bases": bases})
    finally:
        seen.close()
    return docs


def ingest_pass1(objects: list[SourceObject], root: Path, stage: str, k: int, pass1_bytes: int,
                 shard_bases: int = SHARD_BASES, line_width: int = LINE_WIDTH, tmp: Path | None = None,
                 out: Path | None = None, runs_tmp: Path | None = None) -> dict:
    """Pass 1 of shard `k` of `pass1_shards(objects, pass1_bytes)`: its run into scratch (`ingest_runs_dir` under
    `runs_tmp`, the scratch the batch steps share; default `tmp`, this step's own) and its record under `out`
    (default `root`). Refuses, before it reads a source, a stage whose ingest files exist, a shard whose record or
    a plan exists (in `root` or `out`), and store parameters other than the store's."""
    root = Path(root)
    out = root if out is None else Path(out)
    ledger.check_scratch(tmp, root, out)
    ledger.check_scratch(runs_tmp, root, out)
    tmp = _runs_base(tmp, runs_tmp)
    _refuse_stage_written(root, out, stage, stage_dirs(stage))
    shards = pass1_shards(objects, pass1_bytes)
    if isinstance(k, bool) or not isinstance(k, int) or not 0 <= k < len(shards):
        raise ValueError(f"pass-1 shard {k!r}: stage {stage} has {len(shards)} at {pass1_bytes} bytes per shard")
    rel = _pass1_rel(stage, k)
    for r, what in ((rel, f"pass-1 shard {k}"), (_plan_rel(stage), "the plan")):
        if _find(root, out, r) is not None:
            raise FileExistsError(f"{r}: {what} of stage {stage} is written already; refusing to run it again")
    params = _params(shard_bases, line_width)
    _check_params(root, out, params)
    path = ingest_runs_dir(tmp, root, stage) / f"pass1-{k:05d}.parquet"
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        docs = _pass1_run(pass1_order(objects), shards[k], part)
    except BaseException:
        part.unlink(missing_ok=True)
        raise
    os.replace(part, path)  # a run left by an attempt that died before its record is replaced
    rec = {"format": PLAN_FORMAT, "stage": stage, "shard": k, "shards": len(shards), "pass1_bytes": pass1_bytes,
           "store_params": params, "objects": docs, "records": sum(d["records"] for d in docs),
           "run_sha256": ids.sha256_file(path).hex(), "run_bytes": path.stat().st_size}
    _write_once(root, out, rel, _json_bytes(rec), f"the record of pass-1 shard {k}")
    return {"shard": k, "shards": len(shards), "objects": len(docs), "records": rec["records"]}


def _pass1_records(objects: list[SourceObject], root: Path, out: Path | None,
                   stage: str) -> tuple[list[SourceObject], list[tuple[bytes, dict]]]:
    """(pass1_order(objects), every pass-1 record as (bytes, doc)), each checked to be of these objects."""
    _, r0 = _doc(root, out, _pass1_rel(stage, 0), f"record of pass-1 shard 0 of stage {stage}")
    shards = pass1_shards(objects, r0["pass1_bytes"])
    if len(shards) != r0["shards"]:
        raise ValueError(f"stage {stage}: the objects make {len(shards)} pass-1 shards, the records {r0['shards']}; "
                         "the records are of another object list")
    order = pass1_order(objects)
    recs = []
    for k, members in enumerate(shards):
        raw, r = _doc(root, out, _pass1_rel(stage, k), f"record of pass-1 shard {k} of stage {stage}; run that shard")
        if (r["shard"], r["shards"], r["pass1_bytes"]) != (k, r0["shards"], r0["pass1_bytes"]) or \
                [(d["index"], d["name"]) for d in r["objects"]] != [(i, Path(order[i].path).name) for i in members]:
            raise ValueError(f"the record of pass-1 shard {k} of stage {stage} is not of these objects")
        recs.append((raw, r))
    return order, recs


def _file_docs(root: Path, files: list[Path]) -> list[dict]:
    return [{"path": f.relative_to(root).as_posix(), "bytes": f.stat().st_size} for f in files]


def _earlier(root: Path) -> dict:
    """The earlier stages' index files a stage's decisions depend on, root-relative with their sizes."""
    return {"seq_locator": _file_docs(root, locator_files(root)), "contig_index": _file_docs(root, index_files(root))}


def _seen_view(con, runs: list[Path], obj_of: list[int]) -> None:
    """The view `seen` over the pass-1 runs, each object numbered by its rank in sha256 order (`obj_of`, by
    pass-1 index): the rows `ingest`'s own pass 1 writes."""
    con.register("obj_map", pa.table({"pass1_idx": pa.array(range(len(obj_of)), pa.uint32()),
                                      "obj_idx": pa.array(obj_of, pa.uint32())}))
    cols = ", ".join(f's."{n}"' for n in _SEEN.names if n != "obj_idx")
    con.execute(f"CREATE VIEW seen AS SELECT {cols}, m.obj_idx FROM read_parquet({_sql_files(runs)}) s "
                "JOIN obj_map m ON s.obj_idx = m.pass1_idx")


def _sha_order(order: list[SourceObject], recs: list[tuple[bytes, dict]]) -> list[dict]:
    """The pass-1 records' objects in `ingest`'s order (`_object_order`: sha256, then the object's fields; then
    the pass-1 index)."""
    docs = {d["index"]: d for _, r in recs for d in r["objects"]}
    return [docs[ci] for ci in sorted(docs, key=lambda ci: (_order_key((bytes.fromhex(docs[ci]["sha256"]), order[ci])),
                                                            ci))]


def _make_plan(order: list[SourceObject], recs: list[tuple[bytes, dict]], runs: list[Path], root: Path, stage: str,
               params: dict, dest: Path, tmp: Path | None) -> dict:
    """`ingest`'s decisions over the pass-1 runs: the decision table written to `dest`, the plan returned."""
    docs = _sha_order(order, recs)
    obj_of = [0] * len(order)
    for rank, d in enumerate(docs):
        obj_of[d["index"]] = rank
    earlier = _earlier(root)
    scratch = Path(tempfile.mkdtemp(prefix="biobox-ingest-plan-", dir=tmp))
    try:
        con = _duck(scratch)
        try:
            _seen_view(con, runs, obj_of)
            _old_views(con, locator_files(root), index_files(root))
            kept, batches, new_bases = _decide(con, dest, len(docs), params["shard_bases"])
        finally:
            con.close()
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    shard_of = _shard_paths(stage, len(batches))
    return {"format": PLAN_FORMAT, "stage": stage, "store_params": params, "pass1_bytes": recs[0][1]["pass1_bytes"],
            "pass1_sha256": [hashlib.sha256(raw).hexdigest() for raw, _ in recs], "earlier": earlier,
            "objects": [{**d, "kept": kept[i], "kept_bases": new_bases[i]} for i, d in enumerate(docs)],
            "batches": [{"shard": shard_of[j], "objects": b, "new_sequences": sum(kept[i] for i in b),
                         "new_bases": sum(new_bases[i] for i in b), "bases_read": sum(docs[i]["bases"] for i in b)}
                        for j, b in enumerate(batches)],
            "records": sum(d["records"] for d in docs), "new_sequences": sum(kept),
            "new_seq_sha256": ids.sha256_file(dest).hex(), "new_seq_bytes": dest.stat().st_size}


def _runs(recs: list[tuple[bytes, dict]], tmp: Path | None, root: Path, stage: str) -> list[Path]:
    d = ingest_runs_dir(tmp, root, stage)
    return [_checked_run(d / f"pass1-{k:05d}.parquet", r["run_sha256"], f"pass-1 shard {k}", stage, f"--pass1 {k}")
            for k, (_, r) in enumerate(recs)]


def ingest_plan(objects: list[SourceObject], root: Path, stage: str, shard_bases: int = SHARD_BASES,
                line_width: int = LINE_WIDTH, tmp: Path | None = None, out: Path | None = None,
                runs_tmp: Path | None = None) -> dict:
    """The stage's decisions from its pass-1 runs: the decision table into scratch (both under `runs_tmp`, default
    `tmp`; DuckDB spills into `tmp`), the plan under `out` (default `root`). Needs every pass-1 record and a run
    with its sha256. Refuses a stage whose ingest files or plan exist, and store parameters other than the store's
    or the pass-1 records'."""
    root = Path(root)
    out = root if out is None else Path(out)
    ledger.check_scratch(tmp, root, out)
    ledger.check_scratch(runs_tmp, root, out)
    shared = _runs_base(tmp, runs_tmp)
    _refuse_stage_written(root, out, stage, stage_dirs(stage))
    rel = _plan_rel(stage)
    if _find(root, out, rel) is not None:
        raise FileExistsError(f"{rel}: the plan of stage {stage} is written once; refusing to overwrite")
    params = _params(shard_bases, line_width)
    _check_params(root, out, params)
    order, recs = _pass1_records(objects, root, out, stage)
    if any(r["store_params"] != params for _, r in recs):
        raise ValueError(f"stage {stage}: the pass-1 records were made with other store parameters than {params}")
    runs = _runs(recs, shared, root, stage)
    dest = ingest_runs_dir(shared, root, stage) / NEW_SEQ
    part = dest.with_name(f".{dest.name}.tmp-{os.getpid()}")
    try:
        doc = _make_plan(order, recs, runs, root, stage, params, part, tmp)
    except BaseException:
        part.unlink(missing_ok=True)
        raise
    os.replace(part, dest)
    _write_once(root, out, rel, _json_bytes(doc), "the ingest plan")
    return {"objects": len(doc["objects"]), "records": doc["records"], "new_sequences": doc["new_sequences"],
            "batches": len(doc["batches"])}


def _check_plan_objects(order: list[SourceObject], plan: dict, stage: str) -> None:
    if sorted(d["index"] for d in plan["objects"]) != list(range(len(order))) or \
            any(Path(order[d["index"]].path).name != d["name"] for d in plan["objects"]):
        raise ValueError(f"the plan of stage {stage} is of another object list")


def ingest_batch(objects: list[SourceObject], root: Path, stage: str, j: int, tmp: Path | None = None,
                 out: Path | None = None, runs_tmp: Path | None = None) -> dict:
    """Batch `j` of the stage's plan: its objects read again for their new sequences, written as its shard under
    `out` (default `root`), as `ingest` writes it. Each object must have the size and sha256 of pass 1. Refuses a
    stage whose tables exist and a shard that exists, in `root` or `out`; needs the plan's decision table in
    scratch (under `runs_tmp`, default `tmp`) with its sha256."""
    root = Path(root)
    out = root if out is None else Path(out)
    ledger.check_scratch(tmp, root, out)
    ledger.check_scratch(runs_tmp, root, out)
    _refuse_stage_written(root, out, stage, _table_dirs(stage))
    _, plan = _doc(root, out, _plan_rel(stage), f"ingest plan of stage {stage}")
    if isinstance(j, bool) or not isinstance(j, int) or not 0 <= j < len(plan["batches"]):
        raise ValueError(f"batch {j!r}: the plan of stage {stage} has {len(plan['batches'])} batches")
    b = plan["batches"][j]
    for base in dict.fromkeys([root.resolve(), out.resolve()]):
        d = base / b["shard"]
        if d.parent.is_dir() and any(p.name.startswith((d.name, f".{d.name}")) for p in d.parent.iterdir()):
            raise FileExistsError(f"{d}: the shard of batch {j} of stage {stage} is written already (or a "
                                  "temporary file of it); refusing to write it again")
    params = plan["store_params"]
    _check_params(root, out, params)
    order = pass1_order(objects)
    _check_plan_objects(order, plan, stage)
    new = _checked_run(ingest_runs_dir(_runs_base(tmp, runs_tmp), root, stage) / NEW_SEQ, plan["new_seq_sha256"],
                       "the plan's decision table", stage, "--plan")
    docs = plan["objects"]

    def unchanged(i: int, src: Path) -> None:
        if src.stat().st_size != docs[i]["bytes"] or ids.sha256_file(src).hex() != docs[i]["sha256"]:
            raise ValueError(f"{src}: changed since pass 1 (its size or sha256 is not the plan's); the stage must "
                             "be ingested again from the start")

    scratch = Path(tempfile.mkdtemp(prefix="biobox-ingest-batch-", dir=tmp))
    try:
        con = _duck(scratch)
        try:
            con.execute(f"CREATE VIEW new_seq AS SELECT * FROM read_parquet({_sql_path(new)})")
            wanted = _wanted(con, objects=(b["objects"][0], b["objects"][-1]))
            try:
                nxt = _write_shard(out / b["shard"], [(i, Path(order[docs[i]["index"]].path)) for i in b["objects"]],
                                   wanted, next(wanted, None), unchanged, params["line_width"])
            finally:
                wanted.close()
        finally:
            con.close()
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    if nxt is not None:
        raise RuntimeError(f"batch {j} of stage {stage}: the decision table holds sequences outside its objects")
    return {"batch": j, "shard": b["shard"], "objects": len(b["objects"]), "new_sequences": b["new_sequences"],
            "new_bases": b["new_bases"]}


def ingest_merge(objects: list[SourceObject], root: Path, stage: str, tmp: Path | None = None,
                 max_file_bytes: int = layout.MAX_FILE_BYTES, out: Path | None = None,
                 runs_tmp: Path | None = None, _hook: Callable[[str], None] | None = None) -> dict:
    """The stage's tables, ledger, `store_params` (its first stage only) and `l0_manifest` under `out` (default
    `root`), as `ingest` writes them, from the plan and the pass-1 runs. Needs every batch's shard published (in
    `root` or `out`), every scratch run (under `runs_tmp`, default `tmp`) with the sha256 its record holds, and the
    earlier stages' index files the plan decided against. `_hook` (tests only) is called at `ingest`'s points from
    "contig_occurrence" on."""
    root = Path(root)
    out = root if out is None else Path(out)
    hook = _hook or _no_hook
    ledger.check_scratch(tmp, root, out)
    ledger.check_scratch(runs_tmp, root, out)
    shared = _runs_base(tmp, runs_tmp)
    _refuse_stage_written(root, out, stage, _table_dirs(stage))
    raw, plan = _doc(root, out, _plan_rel(stage), f"ingest plan of stage {stage}")
    params = plan["store_params"]
    new_params = _check_params(root, out, params)
    missing = [b["shard"] for b in plan["batches"]
               if any(_find(root, out, b["shard"] + s) is None for s in ("", ".fai", ".gzi"))]
    if missing:
        raise FileNotFoundError(f"stage {stage}: {', '.join(missing)} not published; run those batches first")
    if _earlier(root) != plan["earlier"]:
        raise ValueError(f"stage {stage}: the earlier stages' seq_locator and contig_index files are not those its "
                         "plan decided against; ingest the stage again from the start")
    order, recs = _pass1_records(objects, root, out, stage)
    _check_plan_objects(order, plan, stage)
    if [hashlib.sha256(r).hexdigest() for r, _ in recs] != plan["pass1_sha256"]:
        raise ValueError(f"stage {stage}: the pass-1 records are not those the plan was made from")
    runs = _runs(recs, shared, root, stage)
    new = _checked_run(ingest_runs_dir(shared, root, stage) / NEW_SEQ, plan["new_seq_sha256"],
                       "the plan's decision table", stage, "--plan")
    docs = plan["objects"]
    ordered = [(bytes.fromhex(d["sha256"]), order[d["index"]]) for d in docs]
    obj_of = [0] * len(order)
    for rank, d in enumerate(docs):
        obj_of[d["index"]] = rank
    scratch = Path(tempfile.mkdtemp(prefix="biobox-ingest-merge-", dir=tmp))
    try:
        con = _duck(scratch)
        try:
            _seen_view(con, runs, obj_of)
            _old_views(con, locator_files(root), index_files(root))
            con.execute(f"CREATE VIEW new_seq AS SELECT * FROM read_parquet({_sql_path(new)})")
            report = _write_tables(con, out, stage, ordered, [d["records"] for d in docs], [d["kept"] for d in docs],
                                   [b["objects"] for b in plan["batches"]], max_file_bytes, hook)
        finally:
            con.close()
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    if new_params is not None:
        _write_new(out / PARAMS_PATH, new_params)
    _write_l0_manifest(out, stage, ordered, [(d["bytes"], 0) for d in docs], max_file_bytes)
    hook("l0_manifest")
    return report


def ingest_run_problems(root: Path, stage: str, tmp: Path | None) -> list[tuple[str, int | None, str]]:
    """The scratch runs whose record is published in the root but which cannot be used as they are, as (what,
    shard, why): ("pass1", k, ...) for a pass-1 run, then ("plan", None, ...) for the decision table; why is "is
    missing" or "differs from the sha256 its record holds". Empty when every published run is in place."""
    d = ingest_runs_dir(tmp, root, stage)
    found: list[tuple[str, int | None, str]] = []
    p0 = _find(root, None, _pass1_rel(stage, 0))
    for k in range(json.loads(p0.read_bytes())["shards"] if p0 is not None else 0):
        p = _find(root, None, _pass1_rel(stage, k))
        why = None if p is None else _run_why(d / f"pass1-{k:05d}.parquet", json.loads(p.read_bytes())["run_sha256"])
        if why is not None:
            found.append(("pass1", k, why))
    p = _find(root, None, _plan_rel(stage))
    why = None if p is None else _run_why(d / NEW_SEQ, json.loads(p.read_bytes())["new_seq_sha256"])
    if why is not None:
        found.append(("plan", None, why))
    return found


def rescan_pass1(objects: list[SourceObject], root: Path, stage: str, k: int, tmp: Path | None = None) -> dict:
    """Pass 1 of shard k again into scratch, for a shard whose record is published in the root but whose run is
    gone or damaged. The record is not written again: the new run must reproduce it (its objects and its sha256),
    else it is deleted, the run in scratch is left as it was, and a ValueError names the shard."""
    root = Path(root)
    ledger.check_scratch(tmp, root)
    _refuse_stage_written(root, root, stage, _table_dirs(stage))
    name = f"pass-1 shard {k} of stage {stage}"
    _, rec = _doc(root, None, _pass1_rel(stage, k), f"record of {name} to scan again for")
    shards = pass1_shards(objects, rec["pass1_bytes"])
    order = pass1_order(objects)
    if len(shards) != rec["shards"] or shards[k] != [d["index"] for d in rec["objects"]]:
        raise ValueError(f"{name}: these objects make other pass-1 shards than its record's")
    for d in rec["objects"]:
        src = Path(order[d["index"]].path)
        if src.name != d["name"] or src.stat().st_size != d["bytes"]:
            raise ValueError(f"{name}: {src} is not the object its record lists (name or size); scanned again, it "
                             "cannot reproduce the published record, so the plan and the merge cannot run")
    path = ingest_runs_dir(tmp, root, stage) / f"pass1-{k:05d}.parquet"
    path.parent.mkdir(parents=True, exist_ok=True)
    new = path.with_name(f".rescan-{os.getpid()}-{path.name}")
    try:
        docs = _pass1_run(order, shards[k], new)
    except BaseException:
        new.unlink(missing_ok=True)
        raise
    have = ids.sha256_file(new).hex()
    if docs != rec["objects"] or have != rec["run_sha256"]:
        new.unlink()
        raise ValueError(f"{name}: scanned again, its run has sha256 {have}, the record holds {rec['run_sha256']}; "
                         "the published record is not the pass 1 of these inputs, so the plan and the merge cannot "
                         "run (the run in scratch is left as it was)")
    os.replace(new, path)
    return {"pass1": k, "objects": len(docs), "run_sha256": have}


def rescan_plan(objects: list[SourceObject], root: Path, stage: str, tmp: Path | None = None) -> dict:
    """The plan's decision table again into scratch, for a plan published in the root whose table is gone or
    damaged; every pass-1 run must be in place (rescan those first). The new table and plan must equal the
    published ones, else the table is deleted and a ValueError names the plan."""
    root = Path(root)
    ledger.check_scratch(tmp, root)
    _refuse_stage_written(root, root, stage, _table_dirs(stage))
    _, plan = _doc(root, None, _plan_rel(stage), f"ingest plan of stage {stage} to make again")
    order, recs = _pass1_records(objects, root, None, stage)
    if [hashlib.sha256(r).hexdigest() for r, _ in recs] != plan["pass1_sha256"]:
        raise ValueError(f"the plan of stage {stage}: the pass-1 records are not those it was made from")
    runs = _runs(recs, tmp, root, stage)
    dest = ingest_runs_dir(tmp, root, stage) / NEW_SEQ
    new = dest.with_name(f".rescan-{os.getpid()}-{dest.name}")
    try:
        doc = _make_plan(order, recs, runs, root, stage, plan["store_params"], new, tmp)
    except BaseException:
        new.unlink(missing_ok=True)
        raise
    if doc != plan:
        new.unlink()
        raise ValueError(f"the plan of stage {stage}: made again, its decision table has sha256 "
                         f"{doc['new_seq_sha256']}, the plan holds {plan['new_seq_sha256']} (or the plan differs); "
                         "the published plan is not the one of these inputs, so the batches and the merge cannot run")
    os.replace(new, dest)
    return {"plan": True, "new_sequences": doc["new_sequences"], "new_seq_sha256": doc["new_seq_sha256"]}


def main(argv: list[str] | None = None) -> int:
    """`python -m biobox.l0seq.store rescan --root R --stage S --objects F --tmp D (--pass1 K | --plan)`:
    `rescan_pass1` or `rescan_plan` in its own process, as build-stage runs it (card B-30). It writes nothing into
    the root, so it is not a run-step. Prints one JSON line; exit 2 with the reason on a refusal."""
    ap = argparse.ArgumentParser(prog="python -m biobox.l0seq.store")
    ap.add_argument("action", choices=["rescan"])
    for name in ("--root", "--stage", "--objects", "--tmp"):
        ap.add_argument(name, required=True)
    what = ap.add_mutually_exclusive_group(required=True)
    what.add_argument("--pass1", type=int, help="the pass-1 shard to scan again")
    what.add_argument("--plan", action="store_true", help="make the plan's decision table again")
    a = ap.parse_args(argv)
    from biobox.cli import _source_objects

    try:
        objs = _source_objects(a.objects)
        if a.plan:
            res = rescan_plan(objs, Path(a.root), a.stage, tmp=Path(a.tmp))
        else:
            res = rescan_pass1(objs, Path(a.root), a.stage, a.pass1, tmp=Path(a.tmp))
    except (ValueError, FileNotFoundError, FileExistsError) as e:
        print(f"ingest rescan: {e}", file=sys.stderr)
        return 2
    print(json.dumps({"action": a.action, "result": res, "stage": a.stage}, sort_keys=True))
    return 0


class SeqStore:
    """Reader: slices by (contig_uid, start, end, strand), 0-based half-open.

    `seq_locator` maps a raw_digest to its shard; the shard's `.fai` and `.gzi` find the
    blocks. `counters()` sums blocks decompressed and bytes read over all shards (LY-06).

    Nothing is loaded on open but the Parquet footers of the two indexes: `locator` (by
    `raw_digest`) and `index` (`contig_index` by `contig_uid`) are `KeyLookup`s that read one
    row group per key, with an LRU of `cache_row_groups` decoded row groups each. An unknown
    key raises KeyError.

    The store reads the stages it is opened for, and no other (card B-28; AGENTS §3.9). `stages` is
    either a list of stage names, whose own index directories are listed (a per-stage builder passes its
    stage and the earlier stages its inputs may name), or a `catalog.ReleaseFiles`, whose file lists are
    taken as given (a release's reader: its manifest's files, else those of its `release_def` stages).
    There is no default: a store never lists every stage on disk. A contig or sequence of a stage not
    opened raises KeyError naming the stages.
    """

    def __init__(self, root: Path, stages, cache_blocks: int = 64, cache_row_groups: int = 64):
        self.root = Path(root)
        if isinstance(stages, (str, bytes)) or stages is None or isinstance(stages, int):
            raise TypeError(f"SeqStore needs the stages it reads (a list of names or catalog.ReleaseFiles), "
                            f"got {stages!r}")
        if hasattr(stages, "table"):  # catalog.ReleaseFiles
            self.stages = list(stages.stages)
            self.locator_files = stages.table("seq_locator")
            self.index_files = stages.table("contig_index")
            self._occurrence_files = lambda: stages.table("contig_occurrence")
        else:
            from biobox.catalog import check_stages

            self.stages = check_stages(stages)
            self.locator_files = [f for s in self.stages for f in io.part_files(self.root / LOCATOR_DIR / f"stage={s}")]
            self.index_files = [f for s in self.stages for f in io.part_files(self.root / INDEX_DIR / f"stage={s}")]
            occ = self.root / table("contig_occurrence").location
            self._occurrence_files = lambda: [f for s in self.stages
                                              for f in io.part_files(occ / f"stage={s}", "source=*/part-*.parquet")]
        if not self.locator_files:
            raise FileNotFoundError(f"{self.root / LOCATOR_DIR}: no seq_locator files in stages {self.stages}")
        self.locator = KeyLookup(self.locator_files, "raw_digest", ["shard", "contig_len"], cache_row_groups)
        self.index = KeyLookup(self.index_files, "contig_uid", ["raw_digest", "is_default_variant"], cache_row_groups)
        self._cache_blocks = cache_blocks
        self._shards: dict[str, BgzfFasta] = {}

    def _where(self, rd: bytes) -> tuple[str, int]:
        rows = self.locator.get(rd)
        if not rows:
            raise KeyError(f"raw_digest {rd.hex()} is in none of the stages {self.stages} this store reads")
        return rows[0]["shard"], rows[0]["contig_len"]

    def default_digest(self, contig_uid: bytes) -> bytes:
        """raw_digest of the contig's default variant."""
        for r in self.index.get(contig_uid):
            if r["is_default_variant"]:
                return r["raw_digest"]
        raise KeyError(f"contig_uid {contig_uid.hex()} is in none of the stages {self.stages} this store reads")

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
        """FASTA lines of one source object, from occurrence rows and chunks (R-02); the occurrence rows are
        those of the store's stages."""
        occ = pa.concat_tables(pq.read_table(f) for f in self._occurrence_files())
        rows = sorted(
            (r for r in occ.to_pylist() if r["l0_object_sha256"] == l0_object_sha256),
            key=lambda r: r["record_ordinal"],
        )
        lines: list[bytes] = []
        for r in rows:
            lines.extend(fasta.format_record(r["original_header"], self.raw(r["raw_digest"]), r["line_width"]))
        return lines


if __name__ == "__main__":
    sys.exit(main())
