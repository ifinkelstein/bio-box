"""Per-release search pack: build and verify (SEARCH_PLAN SD-06, SD-07; data_layer_and_tests.md §2B.2; C-12).

    <pack dir>/release=<R>/<pack_id>/
      pack_manifest.json         shards, meta files, counts, sources, rules, builder versions, and the
                                 tape_version the region flanks (feature uids) refer to
      unique/part-<k>.faa.zst    every exact-unique protein, named by protein_uid hex
      igr/part-<k>.fna.zst       IGRs >= 30 nt (+ strand), named by region_uid hex
      igr/igr_meta.parquet       search_region rows for the IGRs
      dark_orf/part-<k>.faa.zst  dark ORFs (biobox.search.orfs), named by region_uid hex
      dark_orf/dark_orf_meta.parquet
      models/<library>.hmm       the release's HMM libraries, re-written as plain HMMER3 text

S0 targets are `unique`, `igr`, `dark_orf` and `models` (ST0: no reduction). `spacers` waits
for CRISPR calls (O-ARRAY) and is listed under `missing_targets`; `nr50`, `nr90`, `segments`
and `igr90` come from S1, and the DIAMOND index of `unique` when DIAMOND enters the image.

A shard is cut once it holds `shard_text_bytes` of uncompressed FASTA (2 GiB: about 1 GiB of
protein and 0.6 GiB of nucleotide zstd), so every shard but the last of a target is over
256 MB. `pack_id` is the BLAKE3-128 of `pack_manifest.json`; the manifest holds no
timestamps, so the same inputs give the same pack_id (D-01).

The pack lives outside the data root (card W-03), in a pack directory (`bio-box pack build
--pack-dir DIR`; by default `<profile tmp>/packs` with a host profile, else beside the search work
directory, `runner.default_pack_base`); a pack directory inside the root is refused. It is a derived,
rebuildable cache of a release, not a release file: the same inputs give the same bytes and the same
`pack_id`, its own manifest records what it was built from, and `verify` (C-12) checks it. It is kept
off the data root because on shared storage every release's pack has files under the 64 MiB floor
(LY-04; planning/measurements/2026-11-10_search_pack_size.md). `build` writes the pack under a
temporary name in `release=<R>/` and renames it to its `pack_id`; it never overwrites a pack, and it is a
no-op, before any work, when a pack built from the same inputs (`_fingerprint`) is there and verifies.

A pack holds the release's files only (card B-28): every input is listed by `catalog.ReleaseFiles`,
from the release manifest when there is one, else from the stages of the release's `release_def`, so
a stage on disk the release does not list never enters its pack.

Bounded memory (card W-16; AGENTS §4.1). Every target's records are in uid order across the whole release
(protein_uid for `unique`, region_uid for `igr` and `dark_orf`), and region uids are hashes, so the order is
a global sort, and the bytes of the pack are those of that order. No target is held whole: it goes through
an external merge sort (`_Sorter`): up to `RUN_ROWS` rows are held, sorted and written as a run to a private
scratch directory; runs are merged `FAN_IN` at a time, each read back `MERGE_BATCH` rows at a time, until
one merge is left, which streams into the FASTA shards and the meta file (row groups of
`META_ROW_GROUP_ROWS`, each written from one contiguous chunk, as the whole-table write before W-16 cut them,
so the bytes are the same). The sort is stable, so rows with equal keys keep the order they were read in,
and duplicates resolve as before. Regions are made one contig at a time. A contig's genes come from
`protein_occurrence` (sorted by protein) regrouped by contig through the same sort; its tape (length,
topology, feature strands), its IGRs and its genetic code are read from files that are each sorted by
`contig_uid`, merged on it a batch per file at a time (`_merged`).

So what is held at a time is: one sorter's buffer (`RUN_ROWS` rows: about 550 bytes per region with its
sequence, about 400 per protein plus its length, so tens of MB), `FAN_IN` runs' read batches, the meta
file's pending row group (at most `META_ROW_GROUP_ROWS` rows, about 12 MB), one contig's sequence and its
ORFs, and one batch of each input file with Arrow's read buffer of its row group. None of these grows with
the number of regions, genes or proteins of the release; the largest contig and the number of input files
bound the last two. `verify` streams too: shard records are counted while decompressing and distinct
proteins through a sorter.

The scratch directory is `release=<R>/.scratch-<pid>-<hex>` beside the build's `.building-<pid>-<hex>`
directory (both outside the data root, both ignored by `find`). Progress goes to stderr, one timestamped
line (`runstep.log`, as `run-step` writes its own: card C-13) per part as it starts and ends, per scratch run
and per shard written. A build stopped by SIGTERM, SIGHUP or SIGINT removes its `.building-*` and scratch
directories, says so on stderr and exits 143, 129 or 130 (SystemExit, or KeyboardInterrupt for SIGINT): the
next call starts again from nothing (a part is not resumed). A build killed by SIGKILL or the kernel's
OOM killer cannot clean up; its two directories are left for removal by exact path.
"""

from __future__ import annotations

import contextlib
import heapq
import itertools
import json
import os
import shutil
import signal
import threading
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Iterator
from operator import itemgetter
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from blake3 import blake3

from biobox import ids, io, layout, provenance, runstep
from biobox.catalog import ReleaseFiles
from biobox.l0seq.store import SeqStore
from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.l1.igr import igr_version
from biobox.l3.geometry import junction
from biobox.l3.tapes import tape_version
from biobox.registry.schemas import UID, table
from biobox.search import fastz, orfs

PACK_FORMAT = "biobox-pack-1"
SHARD_TEXT_BYTES = 2 * 2**30
IGR_MIN_NT = 30
Z = {"z_prot": domains.Z_PROT, "dom_z": domains.DOM_Z, "z_nt_mb": 1e5}  # SD-04, DEC-04
TARGET_KIND = {"unique": "protein", "igr": "igr", "dark_orf": "dark_orf"}
MISSING_AT_S0 = {"spacers": "no CRISPR spacer calls in this release (O-ARRAY not built)"}
# card W-16, the external sort (`_Sorter`): rows held before a sorted run is written (tens of MB: module docstring),
# runs merged at once, and rows per row group of a run and per batch read back from one
RUN_ROWS = 65_536
FAN_IN = 64
MERGE_BATCH = 1_024
META_ROW_GROUP_ROWS = io.ROW_GROUP_ROWS  # row groups of <target>_meta.parquet (AGENTS §4.6; as before W-16)
BATCH_ROWS = 8_192  # rows per batch read from an input file, and per batch of meta rows
TAPE_BATCH_ROWS = 64  # tape rows (contigs, each with its slot lists) per batch read from a tape file
_STOP_CODES = {signal.SIGTERM: 143, signal.SIGHUP: 129}  # exit codes of a build stopped by a signal (128 + signum)

MANIFEST = "pack_manifest.json"
# the manifest keys that, with each target's `sources`, say what a pack was built from (`_fingerprint`)
FINGERPRINT_KEYS = ("pack_format", "release_id", "gene_call_version", "igr_version", "tape_version", "z", "rules",
                    "builder", "missing_targets")
_REGION = table("search_region").schema
# a region row in the sorter: its sort key (region_uid, sequence), then the other search_region columns
_REGION_SORT = pa.schema([_REGION.field("region_uid"), pa.field("seq", pa.string()), *list(_REGION)[1:]])
_PROTEIN_SORT = pa.schema([pa.field("protein_uid", UID), pa.field("seq", pa.string())])
_GENE_SORT = pa.schema([table("protein_occurrence").schema.field(c) for c in ("contig_uid", "start", "end", "strand")])


class NoPack(FileNotFoundError):
    """No pack of the release in the pack directory; the message names the command that builds one."""


class Stopped(SystemExit):
    """A build ended by SIGTERM or SIGHUP, after its directories were removed; the exit code is 128 + signum."""


def _log(msg: str) -> None:
    runstep.log("pack build", msg)


def release_dir(packs: Path, release: str) -> Path:
    """Where the packs of a release live in the pack directory `packs`."""
    return Path(packs) / f"release={release}"


def pack_dir(packs: Path, release: str, pack_id: bytes | str) -> Path:
    """The pack of `release` with `pack_id` in the pack directory `packs`: `<packs>/release=<R>/<pack_id>`."""
    return release_dir(packs, release) / (pack_id if isinstance(pack_id, str) else pack_id.hex())


def outside_root(root: Path, path: Path, what: str) -> Path:
    """`path`, after refusing (ValueError) a path that is the data root or lies under it (card W-03: nothing in
    the search path writes under the root)."""
    r, p = Path(root).resolve(), Path(path).resolve()
    if p == r or p.is_relative_to(r):
        raise ValueError(f"{what} {path} lies under the data root {root}; the search path writes nothing under the "
                         "root: give a directory outside it")
    return Path(path)


def build_command(root: Path, release: str, packs: Path) -> str:
    """The `pack build` command a missing pack's refusal names (card W-16: `--model` is optional; it is needed
    only for the release's own domain libraries, and an `hmm_file` or `sequences` seed needs none)."""
    return (f"bio-box pack build --root {root} --release {release} --pack-dir {packs} "
            "(add --model <library>-<version>=<HMM file> for each domain library of the release, if it has any)")


def find(packs: Path, release: str, pack_id: bytes | str | None = None, *, root: Path | None = None) -> Path:
    """The pack directory of `release` in `packs` (of `pack_id`, when given). NoPack, naming the command that
    builds it, when there is none; ValueError when several packs of the release are there and none is named."""
    if pack_id is not None:
        d = pack_dir(packs, release, pack_id)
        found = [d] if (d / MANIFEST).is_file() else []
    else:
        base = release_dir(packs, release)
        found = sorted(p for p in base.iterdir() if not p.name.startswith(".") and (p / MANIFEST).is_file()) \
            if base.is_dir() else []
    if not found:
        what = "search pack " if pack_id is None else \
            f"search pack {pack_id if isinstance(pack_id, str) else pack_id.hex()} "
        raise NoPack(f"no {what}of release {release} in the pack directory {packs}; build it with: "
                     f"{build_command(root or '<data root>', release, packs)}")
    if len(found) > 1:
        raise ValueError(f"{len(found)} search packs of release {release} in {release_dir(packs, release)} "
                         f"({', '.join(p.name for p in found)}); keep the one to search and remove the others")
    return found[0]


def _rel(root: Path, p: Path) -> str:
    return p.relative_to(root).as_posix()


def _source_hash(root: Path, paths: list[Path]) -> dict:
    h = blake3()
    for p in sorted(paths):
        h.update(f"{_rel(root, p)}:{ids.sha256_file(p).hex()}\n".encode())
    return {"files": len(paths), "blake3": h.hexdigest()}


def _protein_files(rf: ReleaseFiles, gcv: str) -> list[Path]:
    return rf.table("protein_seq", gene_call_version=gcv)


def _occurrence_files(rf: ReleaseFiles, gcv: str) -> list[Path]:
    return rf.table("protein_occurrence", gene_call_version=gcv)


def _igr_files(rf: ReleaseFiles, iv: str) -> list[Path]:
    """The IGR files of igr version iv (`<tape_version>-igr1`, card B-20)."""
    return rf.table("igr", igr_version=iv)


def _tape_files(rf: ReleaseFiles, tv: str) -> list[Path]:
    return rf.table("contig_tapes", tape_version=tv)


def manifest_tape_version(manifest: dict) -> str:
    """The tape version a pack's flanks refer to, and so the one its runs read (card B-16); a manifest
    from before B-16 names none and was built from the gene-call producer's."""
    return manifest.get("tape_version") or tape_version([manifest["gene_call_version"]])


# --- streaming readers (card W-16) ------------------------------------------------------------------------------


def _file_rows(path: Path, cols: list[str], batch_rows: int) -> Iterator[tuple]:
    """The rows of one Parquet file as tuples of `cols`, `batch_rows` at a time; the file must be sorted by its
    first column (ValueError otherwise, since `_merged` relies on it)."""
    prev = None
    for b in pq.ParquetFile(path).iter_batches(batch_size=batch_rows, columns=cols):
        lists = [b.column(c).to_pylist() for c in cols]
        del b
        keys = lists[0]
        if keys and ((prev is not None and keys[0] < prev) or any(a > b for a, b in zip(keys, keys[1:]))):
            raise ValueError(f"{path}: rows are not sorted by {cols[0]}")
        if keys:
            prev = keys[-1]
        yield from zip(*lists)


def _merged(files: list[Path], cols: list[str], batch_rows: int | None = None) -> Iterator[tuple]:
    """The rows of `files` (each sorted by `cols[0]`) as tuples of `cols`, merged in `cols[0]` order; equal keys
    come in file order, then row order, as reading the files one after another gives them. One batch of
    `batch_rows` rows (default `BATCH_ROWS`) of each file is held at a time (and Arrow's read buffer of its row
    group)."""
    n = batch_rows or BATCH_ROWS
    return heapq.merge(*(_file_rows(f, cols, n) for f in files), key=itemgetter(0))


def _grouped(rows: Iterator[tuple]) -> Iterator[tuple[bytes, list[tuple]]]:
    """(key, its rows) for rows that come in key order (the key is each row's first field)."""
    cur, held = None, []
    for r in rows:
        if r[0] != cur:
            if held:
                yield cur, held
            cur, held = r[0], []
        held.append(r)
    if held:
        yield cur, held


class _Keyed:
    """A stream of (key, value) in key order, read forward: `take(key)` is the value of `key`, or `default`
    when the stream has none (keys before it are skipped)."""

    def __init__(self, it: Iterator[tuple], default=None):
        self._it, self._default = iter(it), default
        self._head = next(self._it, None)

    def take(self, key):
        while self._head is not None and self._head[0] < key:
            self._head = next(self._it, None)
        if self._head is not None and self._head[0] == key:
            v = self._head[1]
            self._head = next(self._it, None)
            return v
        return self._default


def _contig_tapes(files: list[Path], with_strands: bool) -> Iterator[tuple[bytes, tuple[int, bool, dict | None]]]:
    """(contig_uid, (contig_len, is_circular, strand of each feature uid or None)) per contig of the tapes, in
    contig_uid order; a contig on several tape rows takes the last row's length, is circular when any row says
    so, and has the strands of all its rows (as the whole-release dictionaries before W-16 had them)."""
    cols = ["contig_uid", "contig_len", "is_circular"] + (["feature_uids", "strands"] if with_strands else [])
    for cu, rows in _grouped(_merged(files, cols, TAPE_BATCH_ROWS)):
        strands = None
        if with_strands:
            strands = {}
            for r in rows:
                strands.update(zip(r[3], r[4]))
        yield cu, (rows[-1][1], any(r[2] for r in rows), strands)


def _contig_codes(files: list[Path]) -> Iterator[tuple[bytes, int]]:
    """(contig_uid, genetic code) in contig_uid order: the commonest code among the contig's CDS features (ties:
    the code met first in file order, then row order, as a Counter fed in that order gives it)."""
    for cu, rows in _grouped(_merged(files, ["contig_uid", "genetic_code"])):
        yield cu, Counter(r[1] for r in rows).most_common(1)[0][0]


def _check_codes(rf: ReleaseFiles, gcv: str) -> None:
    """Fail fast, before anything is written, on a genetic code the dark-ORF stage cannot translate."""
    ok: set[int] = set()
    for cu, code in _contig_codes(rf.table("feature_protein", gene_call_version=gcv)):
        if code in ok:
            continue
        try:
            orfs.check_code(code)
        except ValueError as e:
            raise ValueError(f"{e}; first contig using it: {cu.hex()}") from None
        ok.add(code)


def _chunks(rows: Iterable[tuple], n: int) -> Iterator[list[tuple]]:
    it = iter(rows)
    while chunk := list(itertools.islice(it, n)):
        yield chunk


def _table(rows: list[tuple], schema: pa.Schema) -> pa.Table:
    """Rows (tuples in the order of `schema`) as one table, each column one contiguous array."""
    cols = list(zip(*rows)) if rows else [()] * len(schema)
    return pa.Table.from_arrays([pa.array(c, type=f.type) for c, f in zip(cols, schema)], schema=schema)


class _Sorter:
    """An external merge sort (card W-16) of rows (tuples in the order of `schema`) by `key`, stable: rows with
    equal keys come out in the order they were added. Up to `RUN_ROWS` rows are held; then they are sorted and
    written as a run (Parquet, row groups of `MERGE_BATCH` rows) under `directory`. `rows()` merges the runs
    `FAN_IN` at a time into longer runs until at most `FAN_IN` are left and streams the last merge, reading
    each run `MERGE_BATCH` rows at a time; with no run written, the held rows are sorted in memory. Progress
    lines name `what` (none when it is None)."""

    def __init__(self, directory: Path, schema: pa.Schema, key: Callable | None, what: str | None):
        self.dir, self.schema, self.key, self.what = directory, schema, key, what
        self.n = 0
        self._held: list[tuple] = []
        self._runs: list[Path] = []
        self._k = 0

    def add(self, row: tuple) -> None:
        self._held.append(row)
        self.n += 1
        if len(self._held) >= RUN_ROWS:
            self._spill()

    def _write(self, rows: Iterable[tuple]) -> Path:
        self.dir.mkdir(parents=True, exist_ok=True)
        path = self.dir / f"run-{self._k:06d}.parquet"
        self._k += 1
        with pq.ParquetWriter(path, self.schema, compression="zstd", compression_level=1) as w:
            for chunk in _chunks(rows, MERGE_BATCH):
                w.write_table(_table(chunk, self.schema), row_group_size=MERGE_BATCH)
        return path

    def _spill(self) -> None:
        self._held.sort(key=self.key)
        self._runs.append(self._write(self._held))
        self._held = []
        if self.what and len(self._runs) % 16 == 1:
            _log(f"{self.what}: {self.n:,} rows read, sorted run {len(self._runs)} written")

    def _read(self, path: Path) -> Iterator[tuple]:
        names = self.schema.names
        for b in pq.ParquetFile(path).iter_batches(batch_size=MERGE_BATCH):
            yield from zip(*(b.column(n).to_pylist() for n in names))

    def _merge(self, paths: list[Path]) -> Iterator[tuple]:
        return heapq.merge(*(self._read(p) for p in paths), key=self.key)

    def rows(self) -> Iterator[tuple]:
        if not self._runs:
            self._held.sort(key=self.key)
            held, self._held = self._held, []
            yield from held
            return
        if self._held:
            self._spill()
        runs = self._runs
        while len(runs) > FAN_IN:
            if self.what:
                _log(f"{self.what}: merging {len(runs)} sorted runs, {FAN_IN} at a time")
            merged = []
            for i in range(0, len(runs), FAN_IN):
                group = runs[i:i + FAN_IN]
                if len(group) == 1:
                    merged.append(group[0])
                    continue
                merged.append(self._write(self._merge(group)))
                for p in group:
                    p.unlink()
            runs = merged
        self._runs = runs
        yield from self._merge(runs)

    def close(self) -> None:
        shutil.rmtree(self.dir, ignore_errors=True)


def _gc(seq: bytes) -> float:
    s = seq.upper()
    acgt = sum(s.count(b) for b in (b"A", b"C", b"G", b"T"))
    return (s.count(b"G") + s.count(b"C")) / acgt if acgt else 0.0


def _region(kind: str, cu: bytes, start: int, end: int, strand: int, n_residues: int, **kw) -> dict:
    row = {"region_uid": ids.region_uid(cu, start, end, strand), "target_kind": kind, "contig_uid": cu,
           "start": start, "end": end, "strand": strand, "n_residues": n_residues,
           "orientation": None, "gc": None, "dist_to_contig_end_bp": None,
           "partial": False, "at_contig_end": False, "gene_overlap_frac": None, "genetic_code": None,
           "left_feature_uid": None, "right_feature_uid": None}
    row.update(kw)
    return row


# --- writers ---------------------------------------------------------------------------------------------------


class _Shards(fastz.ShardWriter):
    """`fastz.ShardWriter` that logs each shard as it is closed."""

    def __init__(self, directory: Path, suffix: str, shard_text_bytes: int, target: str):
        super().__init__(directory, suffix, shard_text_bytes)
        self.target = target

    def _closed(self) -> None:
        s = self.shards[-1]
        _log(f"{self.target}: shard {s['path'].name} written: {s['n_seq']:,} sequences, {s['n_residues']:,} "
             f"residues, {s['path'].stat().st_size:,} bytes")

    def close(self) -> list[dict]:
        was_open = self._out is not None
        out = super().close()
        if was_open:
            self._closed()
        return out

    def done(self) -> list[dict]:
        """Close, as `fastz.write_records` ends: an empty target is one empty shard."""
        if not self.shards:
            self._open()
        return self.close()


class _Meta:
    """`<target>_meta.parquet`, written a row group at a time: the same bytes as `pq.write_table` of the whole
    table with `row_group_size=META_ROW_GROUP_ROWS` (each row group from one contiguous chunk per column)."""

    def __init__(self, path: Path):
        self.path = path
        self._w = pq.ParquetWriter(path, _REGION, compression="zstd")
        self._pending: list[pa.RecordBatch] = []
        self._n = self.groups = self.rows = 0

    def _flush(self) -> None:
        t = pa.Table.from_batches(self._pending, schema=_REGION).combine_chunks()
        self._w.write_table(t, row_group_size=META_ROW_GROUP_ROWS)
        self._pending, self._n = [], 0
        self.groups += 1

    def add(self, batch: pa.RecordBatch) -> None:
        off = 0
        while off < batch.num_rows:
            take = min(META_ROW_GROUP_ROWS - self._n, batch.num_rows - off)
            self._pending.append(batch.slice(off, take))
            self._n += take
            self.rows += take
            off += take
            if self._n == META_ROW_GROUP_ROWS:
                self._flush()

    def close(self) -> None:
        if self._n or not self.groups:  # no rows at all: the empty table, as pq.write_table writes it
            self._flush()
        self._w.close()


def _region_row(row: dict, seq: str) -> tuple:
    """A region as a row of `_REGION_SORT`."""
    return (row["region_uid"], seq, *(row[n] for n in _REGION.names[1:]))


def _write_regions(regions: _Sorter, out: Path, target: str, suffix: str, shard_bytes: int) -> list[dict]:
    """The sorted regions into the target's FASTA shards and its meta file; the shards' list. A region on two
    tape rows (an IGR of a contig on two stages) is the same row twice, as it was before W-16."""
    _log(f"{target}: sorting {regions.n:,} regions")
    w = _Shards(out / target, suffix, shard_bytes, target)
    meta = _Meta(out / f"{target}/{target}_meta.parquet")
    try:
        for chunk in _chunks(regions.rows(), BATCH_ROWS):
            for r in chunk:
                w.write(r[0].hex(), r[1])
            meta.add(_table([(r[0], *r[2:]) for r in chunk], _REGION).to_batches()[0])
        shards = w.done()
    finally:
        meta.close()
        regions.close()
    _log(f"{target}: {meta.path.name} written: {meta.rows:,} rows in {meta.groups} row groups")
    return shards


# --- the targets -------------------------------------------------------------------------------------------------


def _build_unique(rf: ReleaseFiles, out: Path, gcv: str, scratch: Path, shard_bytes: int) -> list[dict]:
    """Every distinct protein in protein_uid order; a protein in several files takes the last file's sequence
    (file order, then row order), as one dictionary updated file by file gave it before W-16: the sort is
    stable and rows go in in that order, so the last row of each protein_uid is kept."""
    files = _protein_files(rf, gcv)
    _log(f"unique: start, {len(files)} protein_seq files, {_n_rows(files):,} rows")
    proteins = _Sorter(scratch / "unique", _PROTEIN_SORT, itemgetter(0), "unique")
    w = _Shards(out / "unique", ".faa.zst", shard_bytes, "unique")
    try:
        for f in files:
            for b in pq.ParquetFile(f).iter_batches(batch_size=BATCH_ROWS, columns=["protein_uid", "seq"]):
                for row in zip(b.column("protein_uid").to_pylist(), b.column("seq").to_pylist()):
                    proteins.add(row)
        _log(f"unique: sorting {proteins.n:,} rows")
        last = None
        for row in proteins.rows():
            if last is not None and row[0] != last[0]:
                w.write(last[0].hex(), last[1].replace("*", ""))
            last = row
        if last is not None:
            w.write(last[0].hex(), last[1].replace("*", ""))
        shards = w.done()
    finally:
        proteins.close()
    _log(f"unique: done, {sum(s['n_seq'] for s in shards):,} proteins in {len(shards)} shards")
    return shards


def _slice(contig: bytes, start: int, end: int, circular: bool) -> bytes:
    """`SeqStore.get(cu, start, end, 1, circular)` of a contig whose whole forward sequence is `contig`."""
    n = len(contig)
    if circular and end > n:
        if not 0 <= start < n < end <= start + n:
            raise ValueError(f"wrapping slice {start}:{end} invalid for length {n}")
        return contig[start:] + contig[:end - n]
    if not 0 <= start <= end <= n:
        raise ValueError(f"slice {start}:{end} outside 0:{n}")
    return contig[start:end]


def _n_rows(files: list[Path]) -> int:
    return sum(pq.ParquetFile(f).metadata.num_rows for f in files)


class _Progress:
    """A progress line each time another twentieth of `total` contigs is done."""

    def __init__(self, what: str, total: int, sorter: _Sorter):
        self.what, self.total, self.sorter, self.done, self._next = what, total, sorter, 0, 1

    def tick(self) -> None:
        self.done += 1
        if self.total and self.done * 20 >= self._next * self.total:
            self._next = self.done * 20 // self.total + 1
            _log(f"{self.what}: {self.done:,} of {self.total:,} contigs, {self.sorter.n:,} regions so far")


def _build_igr(rf: ReleaseFiles, out: Path, iv: str, tv: str, store: SeqStore | None, scratch: Path,
               shard_bytes: int) -> list[dict]:
    files = _igr_files(rf, iv)
    tape_files = _tape_files(rf, tv) if store is not None else []
    _log(f"igr: start, {len(files)} IGR files ({_n_rows(files):,} rows), {len(tape_files)} tape files")
    regions = _Sorter(scratch / "igr", _REGION_SORT, itemgetter(0, 1), "igr")
    try:
        tapes = _Keyed(_contig_tapes(tape_files, True))
        cols = ["contig_uid", "start", "end", "left_feature_uid", "right_feature_uid"]
        for cu, rows in _grouped(_merged(files, cols)):
            rows = [r for r in rows if r[2] - r[1] >= IGR_MIN_NT]
            if not rows:
                continue
            tape = tapes.take(cu)
            if tape is None:
                raise ValueError(f"IGR on contig {cu.hex()}, which has no tape of tape version {tv} in the release")
            L, circ, strands = tape
            contig = store.get(cu)  # the whole contig once; each IGR is a slice of it
            for _, s, e, lf, rt in rows:
                seq = _slice(contig, s, e, circ)
                # flanks of an IGR never nest, so only their strands matter
                orient = junction((0, 1, strands[lf]), (2, 3, strands[rt])) if lf and rt else None
                row = _region("igr", cu, s, e, 1, e - s, left_feature_uid=lf, right_feature_uid=rt,
                              orientation=orient, gc=_gc(seq), dist_to_contig_end_bp=None if circ else min(s, L - e),
                              partial=not circ and (s == 0 or e == L))
                regions.add(_region_row(row, seq.decode("ascii")))
        shards = _write_regions(regions, out, "igr", ".fna.zst", shard_bytes)
    finally:
        regions.close()
    _log(f"igr: done, {sum(s['n_seq'] for s in shards):,} IGRs in {len(shards)} shards")
    return shards


def _contig_genes(occ: list[Path], scratch: Path) -> _Sorter:
    """The genes (contig_uid, start, end, strand) of `protein_occurrence`, which is sorted by protein, in a sorter
    that gives them back in (contig_uid, start, end, strand) order."""
    genes = _Sorter(scratch / "genes", _GENE_SORT, None, "dark_orf genes")
    for f in occ:
        for b in pq.ParquetFile(f).iter_batches(batch_size=BATCH_ROWS, columns=list(_GENE_SORT.names)):
            for row in zip(*(b.column(n).to_pylist() for n in _GENE_SORT.names)):
                genes.add(row)
    return genes


def _build_dark_orf(rf: ReleaseFiles, out: Path, gcv: str, tv: str, store: SeqStore | None, scratch: Path,
                    shard_bytes: int) -> list[dict]:
    tape_files = _tape_files(rf, tv) if store is not None else []
    occ = _occurrence_files(rf, gcv)
    n_contigs = _n_rows(tape_files)
    _log(f"dark_orf: start, {len(tape_files)} tape files ({n_contigs:,} contigs), {len(occ)} protein_occurrence files")
    regions = _Sorter(scratch / "dark_orf", _REGION_SORT, itemgetter(0, 1), "dark_orf")
    gene_rows = _contig_genes(occ, scratch) if tape_files else _Sorter(scratch / "genes", _GENE_SORT, None, "genes")
    try:
        _log(f"dark_orf: {gene_rows.n:,} genes read; finding dark ORFs contig by contig")
        genes = _Keyed(((cu, [r[1:] for r in rows]) for cu, rows in _grouped(gene_rows.rows())), [])
        codes = _Keyed(_contig_codes(rf.table("feature_protein", gene_call_version=gcv)), 11)
        progress = _Progress("dark_orf", n_contigs, regions)
        for cu, (_, circ, _) in _contig_tapes(tape_files, False):
            code = codes.take(cu)
            for o, frac, at_end in orfs.dark_orfs(store.get(cu), circ, genes.take(cu), code):
                row = _region("dark_orf", cu, o.start, o.end, o.strand, len(o.aa), partial=o.partial,
                              at_contig_end=at_end, gene_overlap_frac=frac, genetic_code=code)
                regions.add(_region_row(row, o.aa))
            progress.tick()
        shards = _write_regions(regions, out, "dark_orf", ".faa.zst", shard_bytes)
    finally:
        gene_rows.close()
        regions.close()
    _log(f"dark_orf: done, {sum(s['n_seq'] for s in shards):,} dark ORFs in {len(shards)} shards")
    return shards


def _build_models(out: Path, models: dict[str, Path]) -> list[dict]:
    shards = []
    (out / "models").mkdir(parents=True, exist_ok=True)
    _log(f"models: start, {len(models)} libraries")
    for lib in sorted(models):
        path = out / f"models/{lib}.hmm"
        n, m = 0, 0
        with open(path, "wb") as fh:
            for hmm in domains.iter_hmms(models[lib]):
                hmm.write(fh)
                n += 1
                m += hmm.M
        shards.append({"path": path, "n_seq": n, "n_residues": m})
        _log(f"models: {path.name} written: {n:,} models")
    return shards


def _fingerprint(manifest: dict) -> dict:
    """What a pack was built from: the build's rules, versions and tools, and each target's sources (input
    files by path and sha256, library files by sha256). The pack's bytes are a function of these (D-01)."""
    return {**{k: manifest.get(k) for k in FINGERPRINT_KEYS},
            "sources": {t: info.get("sources") for t, info in sorted(manifest.get("targets", {}).items())}}


@contextlib.contextmanager
def _stop_on_signals():
    """While building: SIGTERM and SIGHUP raise `Stopped` (SIGINT raises KeyboardInterrupt, as by default), so the
    build's cleanup runs; a second signal during the cleanup is ignored. Main thread only (elsewhere, nothing)."""
    if threading.current_thread() is not threading.main_thread():
        yield
        return

    def stop(signum, frame):
        for s in _STOP_CODES:
            signal.signal(s, signal.SIG_IGN)
        raise Stopped(_STOP_CODES[signum])

    previous = {s: signal.signal(s, stop) for s in _STOP_CODES}
    try:
        yield
    finally:
        for s, h in previous.items():
            signal.signal(s, h)


def build(root: Path, release: str, gene_call_version: str = GENE_CALL_VERSION,
          models: dict[str, Path] | None = None, shard_text_bytes: int = SHARD_TEXT_BYTES,
          producers: list[str] | None = None, release_manifest: dict | None = None, *, packs: Path) -> dict:
    """Build the release's pack from the data root into the pack directory `packs` (card W-03), at
    `<packs>/release=<R>/<pack_id>/`; returns the manifest (with `pack_id`).

    `packs` may not lie under the root (ValueError, before anything is read). When a pack of the release
    built from the same inputs is there and verifies, nothing is built or written and its manifest is
    returned. Otherwise the pack is written under a temporary name in `release=<R>/` and renamed to its
    pack_id; if a pack of that id is there already (and does not verify), FileExistsError: a pack is never
    overwritten. A build stopped by SIGTERM or SIGHUP removes its temporary and scratch directories and raises
    `Stopped` (a SystemExit with code 128 + signum; card W-16); so does SIGINT, with KeyboardInterrupt.

    `producers` (default: the gene-call producer alone) name the feature set: the IGRs are those of
    its igr version and the region flanks are features of its tape version, which the manifest
    records as `tape_version` (card B-16). Every input file is the release's (card B-28): from
    `release_manifest` (default: `manifests/release=<R>.json` when it exists), else from the stages of
    the release's `release_def` (FileNotFoundError when it has neither)."""
    root = Path(root).resolve()
    packs = outside_root(root, packs, "the pack directory")
    with _stop_on_signals():
        return _build(root, release, gene_call_version, models, shard_text_bytes, producers, release_manifest, packs)


def _build(root: Path, release: str, gcv: str, models: dict[str, Path] | None, shard_text_bytes: int,
           producers: list[str] | None, release_manifest: dict | None, packs: Path) -> dict:
    rf = ReleaseFiles.of(root, release, release_manifest)
    _log(f"release {release}: checking genetic codes")
    _check_codes(rf, gcv)
    producers = sorted(producers) if producers is not None else [gcv]
    tv, iv = tape_version(producers), igr_version(producers)
    seq_files = rf.artifact("seq_shards")
    _log(f"release {release}: hashing the input files")
    seq_hash = _source_hash(root, seq_files)
    sources = {
        "unique": {"protein_seq": _source_hash(root, _protein_files(rf, gcv))},
        "igr": {"igr": _source_hash(root, _igr_files(rf, iv)), "seq_shards": seq_hash},
        "dark_orf": {"protein_occurrence": _source_hash(root, _occurrence_files(rf, gcv)),
                     "contig_tapes": _source_hash(root, _tape_files(rf, tv)),
                     "seq_shards": seq_hash},
        "models": {lib: {"files": 1, "sha256": ids.sha256_file(p).hex()} for lib, p in sorted((models or {}).items())},
    }
    head = {
        "pack_format": PACK_FORMAT,
        "release_id": release,
        "gene_call_version": gcv,
        "igr_version": iv,
        "tape_version": tv,  # the tape version the region flanks refer to, and the one a run reads
        "z": Z,
        "rules": {"igr_min_nt": IGR_MIN_NT, "dark_orf_min_aa": orfs.MIN_AA, "dark_orf_end_bp": orfs.END_BP,
                  "dark_orf_max_gene_overlap": orfs.MAX_GENE_OVERLAP, "shard_text_bytes": shard_text_bytes},
        "builder": provenance.tool_versions("bio-box", "pyhmmer", "pyarrow"),
        "missing_targets": MISSING_AT_S0,
    }
    want = _fingerprint({**head, "targets": {t: {"sources": s} for t, s in sources.items()}})
    base = release_dir(packs, release)
    for d in sorted(p for p in base.iterdir() if not p.name.startswith(".")) if base.is_dir() else []:
        if (d / MANIFEST).is_file():
            have, pid = read_manifest(d)
            if _fingerprint(have) == want and not verify_dir(root, release, d, release_manifest=release_manifest):
                _log(f"release {release}: pack {pid.hex()} is there, built from the same inputs: nothing to do")
                return {**have, "pack_id": pid.hex()}  # built from the same inputs: a no-op
    base.mkdir(parents=True, exist_ok=True)
    tag = f"{os.getpid()}-{os.urandom(4).hex()}"
    out, scratch = base / f".building-{tag}", base / f".scratch-{tag}"
    out.mkdir()
    try:
        scratch.mkdir()
        _log(f"release {release}: building in {out} (scratch {scratch})")
        manifest, data = _build_into(out, root, rf, head, sources, gcv, iv, tv, models or {}, shard_text_bytes,
                                     scratch)
        final = pack_dir(packs, release, pack_id_of(data))
        if final.exists():
            if not verify_dir(root, release, final, release_manifest=release_manifest):
                shutil.rmtree(out)
                _log(f"release {release}: pack {final.name} is there already and verifies; this build removed")
                return {**manifest, "pack_id": final.name}
            raise FileExistsError(f"{final}: a pack with this pack_id exists and does not verify; a pack is never "
                                  "overwritten: remove it, then build again")
        os.rename(out, final)
    except BaseException as e:
        shutil.rmtree(out, ignore_errors=True)
        shutil.rmtree(scratch, ignore_errors=True)
        if isinstance(e, (Stopped, KeyboardInterrupt)):
            why = "SIGINT" if isinstance(e, KeyboardInterrupt) else f"a signal (exit {e.code})"
            _log(f"release {release}: stopped by {why}; removed {out} and {scratch}; nothing is kept, the next "
                 "call starts again")
        raise
    shutil.rmtree(scratch, ignore_errors=True)
    _log(f"release {release}: pack {final.name} written: {len(manifest['shards']) + len(manifest['meta']) + 1} "
         f"files, {sum(s['bytes'] for s in manifest['shards'] + manifest['meta']):,} bytes in shards and meta")
    return {**manifest, "pack_id": final.name}


def _build_into(out: Path, root: Path, rf: ReleaseFiles, head: dict, sources: dict, gcv: str, iv: str, tv: str,
                models: dict[str, Path], shard_text_bytes: int, scratch: Path) -> tuple[dict, bytes]:
    """Write the pack's files and its manifest into the empty directory `out`, with scratch files in the empty
    directory `scratch`; (manifest, manifest bytes)."""
    # a protein-only root (a synthetic world) has no sequence store: its nucleotide targets are empty
    store = SeqStore(root, rf) if rf.table("seq_locator") else None
    try:
        built = {
            "unique": _build_unique(rf, out, gcv, scratch, shard_text_bytes),
            "igr": _build_igr(rf, out, iv, tv, store, scratch, shard_text_bytes),
            "dark_orf": _build_dark_orf(rf, out, gcv, tv, store, scratch, shard_text_bytes),
            "models": _build_models(out, models),
        }
    finally:
        if store is not None:
            store.close()
    manifest = {
        **head,
        "targets": {t: {"kind": TARGET_KIND.get(t, t), "n_seq": sum(s["n_seq"] for s in sh),
                        "n_residues": sum(s["n_residues"] for s in sh), "sources": sources[t]}
                    for t, sh in built.items()},
        "shards": [{"target": t, "path": _rel(out, s["path"]), "sha256": ids.sha256_file(s["path"]).hex(),
                    "bytes": s["path"].stat().st_size, "n_seq": s["n_seq"], "n_residues": s["n_residues"]}
                   for t, sh in built.items() for s in sh],
        "meta": [{"target": t, "path": f"{t}/{t}_meta.parquet",
                  "sha256": ids.sha256_file(out / f"{t}/{t}_meta.parquet").hex(),
                  "bytes": (out / f"{t}/{t}_meta.parquet").stat().st_size} for t in ("igr", "dark_orf")],
    }
    data = json.dumps(manifest, sort_keys=True, indent=1).encode()
    (out / MANIFEST).write_bytes(data)
    return manifest, data


def pack_id_of(manifest_bytes: bytes) -> bytes:
    return blake3(manifest_bytes).digest(length=16)


def read_manifest(pdir: Path) -> tuple[dict, bytes]:
    """(manifest, pack_id) of the pack in the directory `pdir` (`find`, `pack_dir`)."""
    data = (Path(pdir) / MANIFEST).read_bytes()
    return json.loads(data), pack_id_of(data)


def _count_records(path: Path) -> int:
    """The FASTA records of a zstd shard, counted while decompressing (a header is a line that starts with `>`)."""
    n, at_line_start = 0, True
    with pa.CompressedInputStream(pa.OSFile(str(path)), "zstd") as fh:
        while chunk := fh.read(1 << 24):
            n += chunk.count(b"\n>") + (1 if at_line_start and chunk[:1] == b">" else 0)
            at_line_start = chunk.endswith(b"\n")
    return n


def _count_source(rf: ReleaseFiles, manifest: dict, target: str, out: Path) -> int | None:
    gcv = manifest["gene_call_version"]
    if target == "unique":
        # distinct protein uids through a sorter, in a scratch directory beside the pack (dot-named, so `find`
        # and `build` pass it by), removed after
        uids = _Sorter(out.parent / f".verify-{os.getpid()}-{os.urandom(4).hex()}",
                       pa.schema([pa.field("protein_uid", UID)]), None, None)
        try:
            for f in _protein_files(rf, gcv):
                for b in pq.ParquetFile(f).iter_batches(batch_size=BATCH_ROWS, columns=["protein_uid"]):
                    for u in b.column("protein_uid").to_pylist():
                        uids.add((u,))
            return sum(1 for _ in itertools.groupby(uids.rows()))
        finally:
            uids.close()
    if target == "igr":
        n, low = 0, manifest["rules"]["igr_min_nt"]
        for f in _igr_files(rf, manifest.get("igr_version") or igr_version([gcv])):
            pf = pq.ParquetFile(f)
            for i in range(pf.metadata.num_row_groups):
                t = pf.read_row_group(i, columns=["start", "end"])
                n += pc.sum(pc.greater_equal(pc.subtract(t.column("end"), t.column("start")), low)).as_py() or 0
        return n
    if target == "dark_orf":
        return pq.read_metadata(out / "dark_orf/dark_orf_meta.parquet").num_rows
    return None


def verify(root: Path, release: str, pack_id: bytes | None = None,
           min_shard_bytes: int = layout.TARGET_FILE_BYTES, release_manifest: dict | None = None, *,
           packs: Path) -> list[str]:
    """C-12: problems with the release's pack in the pack directory `packs` (the one of `pack_id`, when
    given; `find`); empty means it verifies. Sequence counts are checked against the release's files under
    `root` (`release_manifest`, as for `build`)."""
    try:
        out = find(packs, release, pack_id, root=root)
    except (FileNotFoundError, ValueError) as e:
        return [str(e)]
    return verify_dir(root, release, out, pack_id, min_shard_bytes, release_manifest)


def verify_dir(root: Path, release: str, out: Path, pack_id: bytes | None = None,
               min_shard_bytes: int = layout.TARGET_FILE_BYTES, release_manifest: dict | None = None) -> list[str]:
    """`verify` of the pack in the directory `out`: also, the directory is named by the manifest's hash."""
    root, out = Path(root).resolve(), Path(out)
    if not (out / MANIFEST).exists():
        return [f"{out}: no {MANIFEST}"]
    try:
        rf = ReleaseFiles.of(root, release, release_manifest)
    except FileNotFoundError as e:
        return [f"release {release}: {e}"]
    manifest, pid = read_manifest(out)
    problems = []
    if pack_id is not None and pack_id != pid:
        problems.append(f"pack_id {pack_id.hex()} differs from the manifest hash {pid.hex()}")
    if out.name != pid.hex():
        problems.append(f"{out}: the directory is not named by the manifest hash {pid.hex()}")
    if manifest.get("release_id") != release:
        problems.append(f"{out}: a pack of release {manifest.get('release_id')}, not {release}")
    listed = {MANIFEST}
    by_target: dict[str, list[dict]] = defaultdict(list)
    for entry in manifest["shards"] + manifest["meta"]:
        listed.add(entry["path"])
        p = out / entry["path"]
        if not p.exists():
            problems.append(f"{entry['path']}: listed but missing")
            continue
        if p.stat().st_size != entry["bytes"]:
            problems.append(f"{entry['path']}: {p.stat().st_size} bytes, manifest says {entry['bytes']}")
        if ids.sha256_file(p).hex() != entry["sha256"]:
            problems.append(f"{entry['path']}: sha256 differs from the manifest")
    for entry in manifest["shards"]:
        by_target[entry["target"]].append(entry)
    for p in sorted(x for x in out.rglob("*") if x.is_file()):
        if _rel(out, p) not in listed:
            problems.append(f"{_rel(out, p)}: on disk but not in the manifest")
    for target, shards in sorted(by_target.items()):
        info = manifest["targets"][target]
        if sum(s["n_seq"] for s in shards) != info["n_seq"]:
            problems.append(f"{target}: shard counts do not sum to the target's n_seq")
        # models are one file per library (small by nature, SEARCH_PLAN §3.2), not size-cut shards
        for s in (sorted(shards, key=lambda s: s["path"])[:-1] if target != "models" else []):
            if s["bytes"] < min_shard_bytes:
                problems.append(f"{s['path']}: {s['bytes']} bytes, under {min_shard_bytes} and not the last shard")
        if target == "models":
            for s in shards:
                p = out / s["path"]
                if p.exists() and sum(1 for _ in domains.iter_hmms(p)) != s["n_seq"]:
                    problems.append(f"{s['path']}: model count differs from the manifest")
            continue
        for s in shards:
            p = out / s["path"]
            if p.exists():
                n = _count_records(p)
                if n != s["n_seq"]:
                    problems.append(f"{s['path']}: {n} records, manifest says {s['n_seq']}")
        expect = _count_source(rf, manifest, target, out)
        if expect is not None and expect != info["n_seq"]:
            problems.append(f"{target}: {info['n_seq']} sequences, its source has {expect}")
    return problems
