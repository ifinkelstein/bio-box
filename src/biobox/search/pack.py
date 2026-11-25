"""Per-release search pack: build and verify (SEARCH_PLAN SD-06, SD-07; data_layer_and_tests.md §2B.2; C-12).

    <pack dir>/release=<R>/<pack_id>/
      pack_manifest.json         shards, meta files, counts, sources, rules, builder versions, and the
                                 tape_version the region flanks (feature uids) refer to
      pack_verified.json         the build's record that the pack verified (card W-19), naming the manifest
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
directory, `runner.default_pack_base`); a pack directory inside the root is refused. The host profile's
`<tmp>/packs` is the host's one pack directory (card W-22): `pack`, `search`, `nbr` and `releases` find a
release's pack there with only `--profile`. It is a derived,
rebuildable cache of a release, not a release file: the same inputs give the same bytes and the same
`pack_id`, its own manifest records what it was built from, and `verify` (C-12) checks it. It is kept
off the data root because on shared storage every release's pack has files under the 64 MiB floor
(LY-04; planning/measurements/2026-11-10_search_pack_size.md). `build` writes the pack under a
temporary name in `release=<R>/` and renames it to its `pack_id`; it never overwrites a pack, and it is a
no-op, before any work, when a pack built from the same inputs (`_fingerprint`) is there and verifies.

Verified once, at build (card W-19). Before it renames the pack into place, `build` runs the checks of `verify`
on it (every file's size and sha256, each shard's record count, each target's count against its sources; not the
shard-size rule, which `pack verify` keeps) and writes `pack_verified.json`: the pack id and the full BLAKE3 of
the manifest's bytes. The record sits beside the manifest, not in it, so the manifest, and the pack id that is its
hash, are what they were before W-19: a rebuild of the same inputs gives the same pack id. A search checks only
the manifest, the record and the files' sizes (`search_check`), reads no shard to do it and writes nothing into
the pack directory. A pack built before W-19 has no record: a search verifies it in full for that run, with its
scratch in the run's work directory (`verify_dir(scratch=...)`), and `build` with the same inputs, a no-op that
verifies the pack in full, adds the record without touching the pack's files. `verify` accepts the record and
checks that it names the manifest.

Nothing but a build writes into the pack directory (card W-22), so it may be read-only for everyone but its
builder. `verify`, and the full verify of an existing pack that a no-op `build` runs, sort the release's protein
uids under `scratch` (the CLI gives the profile's per-user work directory, `<tmp>/search/<user>/work`), else
under the system temporary directory (`tempfile.gettempdir()`), never beside the pack. A no-op `build` of a pack
without a record, in a directory its caller cannot write, adds no record and says so (`_add_record`). A build that
makes a pack writes only into `release=<R>/` of the pack directory: its `.building-*` and `.scratch-*` directories
and the renamed pack.

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

The estimate (card W-35; AGENTS §7.2). One process is about 86 CPU-minutes at 12.3 Gbp (measured, `rN.1`), so a
50 Gbp release passes the 300 CPU-minute budget. Before any work, `pack build` counts the release's `protein_seq`
rows (file metadata) and genome bases (the `contig_len` of its tapes, once per contig: `release_counts`), applies the
rates below (`part_cpu_min`, `one_process_cpu_min`) and prints the estimate. Without `--parts` a build whose
one-process estimate is over `stage.SHARD_MARGIN` of the profile's `cpu_min_per_process` is refused (exit 2) with the
`--parts` that fits (`one_process_refusal`); with `--parts N` each part's estimate must fit (`parts_refusal`). The
one-process build itself (`build`) is unchanged: the same code, the same bytes and pack ids.

The build in parts (`build_in_parts`, `bio-box pack build --parts N --jobs J`). A driver process, then parts, each
its own process (`python -m biobox.search.pack part`), then one last process (`... assemble`):

- The driver checks the genetic codes, hashes the inputs (the manifest's `sources`, and so the fingerprint a
  pack is matched by, as `build` does), is a no-op when a pack of the same fingerprint is there and verifies, cuts
  the release (below) and writes the plan, `<work>/plan-of-<N>.json`. It launches the parts that have no output yet,
  at most J at once, one thread each, within the allowance of the hour (`stage._launch`, as `build-stage` launches
  scan shards), each under the memory cap it inherits (or a `MemoryMax` scope: `runstep.memory_cap`) at the
  profile's `nice_monolithic`; `STOP` in the data root is checked before each launch (exit 75). Then it launches
  the last process the same way, and removes the work directory once the pack is in place.
- A part computes one target's records over one range of the release with the readers of the one-process build,
  sorts them with the same `_Sorter` and writes them as one sorted run, `<work>/<name>/run.parquet`, with
  `done.json` (rows, its CPU-minutes, wall time and peak RSS, its estimate). The parts are `unique` (every
  `protein_seq` row: one part), and `igr-<k>-of-<N>` and `dark_orf-<k>-of-<N>` for k < N, each over the contigs whose
  uid lies in its range [lo, hi). The cut (`_cut`): the release's contigs in uid order (its tapes), a new part
  starting at the first contig whose genome bases before it reach k/N of the release's (the dark-ORF work is linear
  in bases, measured below); a contig is never split, and a part with no contig writes an empty run. Each reader is
  given the range: files sorted by contig uid skip row groups by their statistics and keep the rows in range,
  `protein_occurrence` (sorted by protein) keeps the rows of the range's contigs.
- The last process merges each target's runs in part order (`heapq.merge` on the target's sort key), writes the
  shards and meta files with the writers of the one-process build (`_write_unique`, `_write_region_rows`), the
  models (`_build_models`) and the manifest, verifies the pack (card W-19), writes the record and renames the pack
  into `<pack dir>/release=<R>/<pack id>`.

Where. `<work>` is `<work base>/release=<R>/<key>`, `key` a hash of the fingerprint, under the profile's tmp
(`<tmp>/scratch/pack` unless `--work` names another base), never under the data root and never in the pack
directory: the one thing a build in parts writes there is the last rename, so the work directory must be on the pack
directory's file system (checked before any part). A part writes into `.<name>.tmp-<pid>-<hex>` and renames it to
`<name>` once its run and `done.json` are written, and a reader takes only a part with its `done.json`: a part
stopped by SIGTERM, SIGHUP or SIGINT removes its temporary directory (exit 143, 129, 130), one killed by SIGKILL
leaves it, and the next driver, which holds the work directory's lock (`.lock`, one driver at a time), removes it
before it launches anything. The last process builds the pack in `<work>/.building-<pid>-<hex>` the same way. A part
whose output is there is not run again, so a build stopped between parts resumes at the parts without output.

Why the bytes cannot depend on N, on J or on which process wrote a run. The one-process bytes are those of one
stable sort of a target's rows in reading order: IGRs and dark ORFs are made contig by contig in contig uid order,
and the shards are cut by `shard_text_bytes` of the sorted stream and the meta row groups by `META_ROW_GROUP_ROWS`
rows of it, never by run or process. A part reads the rows of its contigs in that same order, and the ranges are
consecutive in uid order, so the parts in part order concatenate to the one-process reading order; a stable sort
of each part followed by a merge that takes equal keys in part order (`heapq.merge` is stable across its inputs) is
the stable sort of the whole. Equal keys (region uid and sequence) have one contig, since the region uid hashes it,
so they meet inside one part anyway. `unique` is one part reading the files in the one-process order, and keeps the
last row of a protein as it does. A run holds the rows in the sorter's own schema, through which the one-process
sorter's runs also pass; the float columns are cast to the meta schema's float32 in both paths. So the files, the
manifest and the pack id (its hash) are the same, whatever N, J or the order the parts ran in.
"""

from __future__ import annotations

import argparse
import bisect
import contextlib
import fcntl
import heapq
import itertools
import json
import math
import os
import re
import resource
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass
from datetime import datetime
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
# card W-19: the build's record that it verified the pack (C-12), beside the manifest and not in it, so the manifest
# and the pack id (its hash) are what they were before the record existed (module docstring)
VERIFIED = "pack_verified.json"
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


_WHO = threading.local()  # who the progress lines name: "pack build", or "pack part <name>" in a part (card W-35)


def _log(msg: str) -> None:
    runstep.log(getattr(_WHO, "name", "pack build"), msg)


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


@dataclass(frozen=True)
class _Range:
    """The contig uids [lo, hi) of a part (card W-35); None is unbounded on that side; `empty` holds none."""
    lo: bytes | None = None
    hi: bytes | None = None
    empty: bool = False

    @classmethod
    def of(cls, part: dict) -> "_Range":
        """The range of a part of the plan (`_cut`), whose bounds are hex or None."""
        return cls(None if part["lo"] is None else bytes.fromhex(part["lo"]),
                   None if part["hi"] is None else bytes.fromhex(part["hi"]), part["empty"])

    def has(self, key: bytes) -> bool:
        return not self.empty and (self.lo is None or key >= self.lo) and (self.hi is None or key < self.hi)

    def may_hold(self, bounds: tuple[bytes, bytes] | None) -> bool:
        """Whether a row group whose keys lie in `bounds` (min, max; None when unknown) may hold a key in range."""
        if self.empty:
            return False
        if bounds is None:
            return True
        return (self.lo is None or bounds[1] >= self.lo) and (self.hi is None or bounds[0] < self.hi)

    def cut(self, keys: list) -> tuple[int, int]:
        """The slice of `keys` (sorted) in range."""
        if self.empty:
            return 0, 0
        a = 0 if self.lo is None else bisect.bisect_left(keys, self.lo)
        z = len(keys) if self.hi is None else bisect.bisect_left(keys, self.hi)
        return a, max(a, z)


def _group_bounds(md: pq.FileMetaData, i: int, name: str) -> tuple[bytes, bytes] | None:
    """(min, max) of the column `name` in row group `i` from its statistics; None when they are missing."""
    rg = md.row_group(i)
    for j in range(rg.num_columns):
        c = rg.column(j)
        if c.path_in_schema == name:
            s = c.statistics
            if s is not None and s.has_min_max and isinstance(s.min, bytes) and isinstance(s.max, bytes):
                return s.min, s.max
            return None
    return None


def _batches(path: Path, cols: list[str], batch_rows: int, rng: _Range | None = None) -> Iterator[pa.RecordBatch]:
    """Batches of at most `batch_rows` rows of `cols` of one Parquet file, read one row group at a time. Not
    `ParquetFile.iter_batches`: it reads ahead, a whole file at a time (64 readers of 26 MB files held 1.03 GB
    of Arrow memory, against 28 MB read a row group each: measured, wilkcomp02, 2026-11-10, pyarrow of image
    bio-box-test_74b3e05). With `rng` (card W-35), a row group whose statistics on `cols[0]` say it holds no key
    in range is not read (the rows read are not filtered here)."""
    pf = pq.ParquetFile(path)
    for i in range(pf.metadata.num_row_groups):
        if rng is not None and not rng.may_hold(_group_bounds(pf.metadata, i, cols[0])):
            continue
        yield from pf.read_row_group(i, columns=cols, use_threads=False).to_batches(max_chunksize=batch_rows)


def _file_rows(path: Path, cols: list[str], batch_rows: int, rng: _Range | None = None) -> Iterator[tuple]:
    """The rows of one Parquet file as tuples of `cols`, `batch_rows` at a time; the file must be sorted by its
    first column (ValueError otherwise, since `_merged` relies on it). With `rng`, only the rows whose first column
    is in range (card W-35)."""
    prev = None
    for b in _batches(path, cols, batch_rows, rng):
        lists = [b.column(c).to_pylist() for c in cols]
        del b
        keys = lists[0]
        if keys and ((prev is not None and keys[0] < prev) or any(a > b for a, b in zip(keys, keys[1:]))):
            raise ValueError(f"{path}: rows are not sorted by {cols[0]}")
        if keys:
            prev = keys[-1]
        if rng is not None:
            a, z = rng.cut(keys)
            if (a, z) != (0, len(keys)):
                lists = [x[a:z] for x in lists]
        yield from zip(*lists)


def _merged(files: list[Path], cols: list[str], batch_rows: int | None = None,
            rng: _Range | None = None) -> Iterator[tuple]:
    """The rows of `files` (each sorted by `cols[0]`) as tuples of `cols`, merged in `cols[0]` order; equal keys
    come in file order, then row order, as reading the files one after another gives them. One batch of
    `batch_rows` rows (default `BATCH_ROWS`) of each file is held at a time (and Arrow's read buffer of its row
    group). With `rng`, the rows whose key is in range (card W-35), in the same order."""
    n = batch_rows or BATCH_ROWS
    return heapq.merge(*(_file_rows(f, cols, n, rng) for f in files), key=itemgetter(0))


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


def _contig_tapes(files: list[Path], with_strands: bool,
                  rng: _Range | None = None) -> Iterator[tuple[bytes, tuple[int, bool, dict | None]]]:
    """(contig_uid, (contig_len, is_circular, strand of each feature uid or None)) per contig of the tapes, in
    contig_uid order; a contig on several tape rows takes the last row's length, is circular when any row says
    so, and has the strands of all its rows (as the whole-release dictionaries before W-16 had them). With `rng`,
    the contigs in range (card W-35)."""
    cols = ["contig_uid", "contig_len", "is_circular"] + (["feature_uids", "strands"] if with_strands else [])
    for cu, rows in _grouped(_merged(files, cols, TAPE_BATCH_ROWS, rng)):
        strands = None
        if with_strands:
            strands = {}
            for r in rows:
                strands.update(zip(r[3], r[4]))
        yield cu, (rows[-1][1], any(r[2] for r in rows), strands)


def _contig_codes(files: list[Path], rng: _Range | None = None) -> Iterator[tuple[bytes, int]]:
    """(contig_uid, genetic code) in contig_uid order: the commonest code among the contig's CDS features (ties:
    the code met first in file order, then row order, as a Counter fed in that order gives it). With `rng`, the
    contigs in range (card W-35)."""
    for cu, rows in _grouped(_merged(files, ["contig_uid", "genetic_code"], rng=rng)):
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


def _write_run(path: Path, schema: pa.Schema, rows: Iterable[tuple]) -> int:
    """Rows (tuples in the order of `schema`) as a sorted run at `path`, row groups of `MERGE_BATCH` rows; the rows
    written. A sorter's scratch run, and a part's output (card W-35)."""
    n = 0
    # no statistics: a merged run has thousands of row groups, and min/max of sequences would fill its footer
    with pq.ParquetWriter(path, schema, compression="zstd", compression_level=1, write_statistics=False) as w:
        for chunk in _chunks(rows, MERGE_BATCH):
            w.write_table(_table(chunk, schema), row_group_size=MERGE_BATCH)
            n += len(chunk)
    return n


def _read_run(path: Path, schema: pa.Schema) -> Iterator[tuple]:
    """The rows of a run (`_write_run`) as tuples, `MERGE_BATCH` rows read at a time."""
    names = schema.names
    for b in _batches(path, names, MERGE_BATCH):
        yield from zip(*(b.column(n).to_pylist() for n in names))


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
        _write_run(path, self.schema, rows)
        return path

    def _spill(self) -> None:
        self._held.sort(key=self.key)
        self._runs.append(self._write(self._held))
        self._held = []
        if self.what and len(self._runs) % 16 == 1:
            _log(f"{self.what}: {self.n:,} rows read, sorted run {len(self._runs)} written")

    def _read(self, path: Path) -> Iterator[tuple]:
        return _read_run(path, self.schema)

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
    try:
        return _write_region_rows(regions.rows(), out, target, suffix, shard_bytes)
    finally:
        regions.close()


def _write_region_rows(rows: Iterator[tuple], out: Path, target: str, suffix: str, shard_bytes: int) -> list[dict]:
    """Region rows of `_REGION_SORT` in pack order into the target's FASTA shards and its meta file; the shards'
    list. The one writer of both builds (card W-35): the one-process build gives it its sorter's rows, the build in
    parts the merge of its parts' runs."""
    w = _Shards(out / target, suffix, shard_bytes, target)
    meta = _Meta(out / f"{target}/{target}_meta.parquet")
    try:
        for chunk in _chunks(rows, BATCH_ROWS):
            for r in chunk:
                w.write(r[0].hex(), r[1])
            meta.add(_table([(r[0], *r[2:]) for r in chunk], _REGION).to_batches()[0])
        shards = w.done()
    finally:
        meta.close()
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
    try:
        _unique_rows(rf, gcv, proteins)
        _log(f"unique: sorting {proteins.n:,} rows")
        return _write_unique(proteins.rows(), out, shard_bytes)
    finally:
        proteins.close()


def _unique_rows(rf: ReleaseFiles, gcv: str, proteins: _Sorter) -> None:
    """Every `protein_seq` row of the release into the sorter, in file order, then row order."""
    for f in _protein_files(rf, gcv):
        for b in _batches(f, ["protein_uid", "seq"], BATCH_ROWS):
            for row in zip(b.column("protein_uid").to_pylist(), b.column("seq").to_pylist()):
                proteins.add(row)


def _write_unique(rows: Iterator[tuple], out: Path, shard_bytes: int) -> list[dict]:
    """Protein rows in protein_uid order, equal uids in reading order, into the `unique` shards, the last row of each
    uid kept; the shards' list. The one writer of both builds (card W-35)."""
    w = _Shards(out / "unique", ".faa.zst", shard_bytes, "unique")
    last = None
    for row in rows:
        if last is not None and row[0] != last[0]:
            w.write(last[0].hex(), last[1].replace("*", ""))
        last = row
    if last is not None:
        w.write(last[0].hex(), last[1].replace("*", ""))
    shards = w.done()
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
        _igr_regions(rf, iv, tv, store, regions)
        shards = _write_regions(regions, out, "igr", ".fna.zst", shard_bytes)
    finally:
        regions.close()
    _log(f"igr: done, {sum(s['n_seq'] for s in shards):,} IGRs in {len(shards)} shards")
    return shards


def _igr_regions(rf: ReleaseFiles, iv: str, tv: str, store: SeqStore | None, regions: _Sorter,
                 rng: _Range | None = None) -> None:
    """The IGRs of at least `IGR_MIN_NT` of the release (of the contigs in `rng`, card W-35) into the sorter, contig by
    contig in contig_uid order."""
    files = _igr_files(rf, iv)
    tape_files = _tape_files(rf, tv) if store is not None else []
    tapes = _Keyed(_contig_tapes(tape_files, True, rng))
    cols = ["contig_uid", "start", "end", "left_feature_uid", "right_feature_uid"]
    for cu, rows in _grouped(_merged(files, cols, rng=rng)):
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


def _contig_genes(occ: list[Path], scratch: Path, rng: _Range | None = None) -> _Sorter:
    """The genes (contig_uid, start, end, strand) of `protein_occurrence`, which is sorted by protein, in a sorter
    that gives them back in (contig_uid, start, end, strand) order. With `rng`, the genes of the contigs in range
    (card W-35): each batch's contig column is read first and only those rows taken."""
    genes = _Sorter(scratch / "genes", _GENE_SORT, None, "dark_orf genes")
    for f in occ:
        for b in _batches(f, list(_GENE_SORT.names), BATCH_ROWS):
            if rng is not None:
                keep = [i for i, c in enumerate(b.column("contig_uid").to_pylist()) if rng.has(c)]
                if not keep:
                    continue
                if len(keep) < b.num_rows:
                    b = b.take(pa.array(keep, pa.int64()))
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
    try:
        _dark_orf_regions(rf, gcv, tv, store, scratch, regions, n_contigs)
        shards = _write_regions(regions, out, "dark_orf", ".faa.zst", shard_bytes)
    finally:
        regions.close()
    _log(f"dark_orf: done, {sum(s['n_seq'] for s in shards):,} dark ORFs in {len(shards)} shards")
    return shards


def _dark_orf_regions(rf: ReleaseFiles, gcv: str, tv: str, store: SeqStore | None, scratch: Path, regions: _Sorter,
                      n_contigs: int, rng: _Range | None = None) -> None:
    """The dark ORFs of the release's contigs (those in `rng`, card W-35) into the sorter, contig by contig in
    contig_uid order; the genes are sorted by contig under `scratch` first."""
    tape_files = _tape_files(rf, tv) if store is not None else []
    occ = _occurrence_files(rf, gcv)
    gene_rows = _contig_genes(occ, scratch, rng) if tape_files else _Sorter(scratch / "genes", _GENE_SORT, None,
                                                                            "genes")
    try:
        _log(f"dark_orf: {gene_rows.n:,} genes read; finding dark ORFs contig by contig")
        genes = _Keyed(((cu, [r[1:] for r in rows]) for cu, rows in _grouped(gene_rows.rows())), [])
        codes = _Keyed(_contig_codes(rf.table("feature_protein", gene_call_version=gcv), rng), 11)
        progress = _Progress("dark_orf", n_contigs, regions)
        for cu, (_, circ, _) in _contig_tapes(tape_files, False, rng):
            code = codes.take(cu)
            for o, frac, at_end in orfs.dark_orfs(store.get(cu), circ, genes.take(cu), code):
                row = _region("dark_orf", cu, o.start, o.end, o.strand, len(o.aa), partial=o.partial,
                              at_contig_end=at_end, gene_overlap_frac=frac, genetic_code=code)
                regions.add(_region_row(row, o.aa))
            progress.tick()
    finally:
        gene_rows.close()


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
          producers: list[str] | None = None, release_manifest: dict | None = None, *, packs: Path,
          scratch: Path | None = None) -> dict:
    """Build the release's pack from the data root into the pack directory `packs` (card W-03), at
    `<packs>/release=<R>/<pack_id>/`; returns the manifest (with `pack_id`).

    `packs` may not lie under the root (ValueError, before anything is read). When a pack of the release
    built from the same inputs is there and verifies, nothing is built or written and its manifest is
    returned; that verify sorts under `scratch` (default: the system temporary directory), never in the pack
    directory (card W-22). Otherwise the pack is written under a temporary name in `release=<R>/` and renamed
    to its pack_id; if a pack of that id is there already (and does not verify), FileExistsError: a pack is never
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
        return _build(root, release, gene_call_version, models, shard_text_bytes, producers, release_manifest, packs,
                      scratch)


def _build(root: Path, release: str, gcv: str, models: dict[str, Path] | None, shard_text_bytes: int,
           producers: list[str] | None, release_manifest: dict | None, packs: Path, verify_scratch: Path | None
           ) -> dict:
    rf = ReleaseFiles.of(root, release, release_manifest)
    _log(f"release {release}: checking genetic codes")
    _check_codes(rf, gcv)
    producers = sorted(producers) if producers is not None else [gcv]
    tv, iv = tape_version(producers), igr_version(producers)
    head, sources = _head_and_sources(root, rf, release, gcv, iv, tv, models, shard_text_bytes)
    want = _fingerprint({**head, "targets": {t: {"sources": s} for t, s in sources.items()}})
    have = _existing(root, release, packs, want, release_manifest, verify_scratch)
    if have is not None:
        return have  # built from the same inputs: a no-op
    base = release_dir(packs, release)
    base.mkdir(parents=True, exist_ok=True)
    tag = f"{os.getpid()}-{os.urandom(4).hex()}"
    out, scratch = base / f".building-{tag}", base / f".scratch-{tag}"
    out.mkdir()
    try:
        scratch.mkdir()
        _log(f"release {release}: building in {out} (scratch {scratch})")
        manifest, data = _build_into(out, root, rf, head, sources, gcv, iv, tv, models or {}, shard_text_bytes,
                                     scratch)
        final, placed = _place(root, rf, release, out, manifest, data, scratch, packs, release_manifest,
                               verify_scratch)
    except BaseException as e:
        shutil.rmtree(out, ignore_errors=True)
        shutil.rmtree(scratch, ignore_errors=True)
        if isinstance(e, (Stopped, KeyboardInterrupt)):
            why = "SIGINT" if isinstance(e, KeyboardInterrupt) else f"a signal (exit {e.code})"
            _log(f"release {release}: stopped by {why}; removed {out} and {scratch}; nothing is kept, the next "
                 "call starts again")
        raise
    shutil.rmtree(scratch, ignore_errors=True)
    if placed:
        _log(f"release {release}: pack {final.name} written: {len(manifest['shards']) + len(manifest['meta']) + 2} "
             f"files, {sum(s['bytes'] for s in manifest['shards'] + manifest['meta']):,} bytes in shards and meta")
    return {**manifest, "pack_id": final.name}


def _head_and_sources(root: Path, rf: ReleaseFiles, release: str, gcv: str, iv: str, tv: str,
                      models: dict[str, Path] | None, shard_text_bytes: int) -> tuple[dict, dict]:
    """The manifest's head (rules, versions, tools) and each target's `sources` (every input file by path and
    sha256, the libraries by sha256): what the pack is built from (`_fingerprint`)."""
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
    return head, sources


def _existing(root: Path, release: str, packs: Path, want: dict, release_manifest: dict | None,
              verify_scratch: Path | None) -> dict | None:
    """The manifest (with `pack_id`) of a pack of the release in `packs` built from the same inputs (`want`, a
    `_fingerprint`) that verifies in full, its record added when it has none; None when there is no such pack."""
    base = release_dir(packs, release)
    for d in sorted(p for p in base.iterdir() if not p.name.startswith(".")) if base.is_dir() else []:
        if (d / MANIFEST).is_file():
            have, pid = read_manifest(d)
            if _fingerprint(have) == want and not verify_dir(root, release, d, release_manifest=release_manifest,
                                                             scratch=verify_scratch):
                _log(f"release {release}: pack {pid.hex()} is there, built from the same inputs: nothing to do")
                _add_record(d, release)  # a pack built before W-19 has none; it has just verified in full
                return {**have, "pack_id": pid.hex()}
    return None


def _place(root: Path, rf: ReleaseFiles, release: str, out: Path, manifest: dict, data: bytes, scratch: Path,
           packs: Path, release_manifest: dict | None, verify_scratch: Path | None) -> tuple[Path, bool]:
    """Verify the pack just written in `out` (card W-19), write its record and rename it to `<packs>/release=<R>/<pack
    id>`; (that directory, whether this build put it there). A pack of that id already there that verifies is kept and
    `out` removed; one that does not verify is never overwritten (FileExistsError)."""
    pid = pack_id_of(data)
    # card W-19: verify once, here, and record it beside the manifest; a search then checks the record and the
    # sizes only. The shard-size rule is not a property of the files and is left to `pack verify`.
    _log(f"release {release}: verifying the pack (C-12)")
    problems = _problems(rf, release, out, manifest, 0, scratch)
    if problems:
        raise ValueError(f"the pack built in {out} does not verify (C-12): {problems[:3]}")
    _write_record(out, pid, release, data)
    _log(f"release {release}: the pack verifies; record {VERIFIED} written")
    final = pack_dir(packs, release, pid)
    if final.exists():
        if not verify_dir(root, release, final, release_manifest=release_manifest, scratch=verify_scratch):
            shutil.rmtree(out)
            _log(f"release {release}: pack {final.name} is there already and verifies; this build removed")
            _add_record(final, release)
            return final, False
        raise FileExistsError(f"{final}: a pack with this pack_id exists and does not verify; a pack is never "
                              "overwritten: remove it, then build again")
    os.rename(out, final)
    return final, True


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
    return _write_manifest(out, head, sources, built)


def _write_manifest(out: Path, head: dict, sources: dict, built: dict[str, list[dict]]) -> tuple[dict, bytes]:
    """The manifest of the pack whose files `built` lists (per target, in the order unique, igr, dark_orf, models),
    written into `out`; (manifest, its bytes). Both builds end here (card W-35)."""
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


# --- the build's verification record and the search's check (card W-19) -------------------------------------------


def _record_bytes(pid: bytes, release: str, manifest_bytes: bytes) -> bytes:
    """The record a build writes once the pack verifies (`VERIFIED`): which pack, by its id and the full BLAKE3 of
    its manifest's bytes. No timestamp or host: the same inputs give the same bytes (D-01)."""
    return json.dumps({"pack_id": pid.hex(), "release_id": release,
                       "manifest_blake3": blake3(manifest_bytes).hexdigest(),
                       "check": "pack verify (C-12) at build: every listed file's size and sha256, each shard's "
                                "record count, each target's count against its sources", "problems": []},
                      sort_keys=True, indent=1).encode()


def read_record(pdir: Path) -> dict | None:
    """The build's verification record of the pack in `pdir`, None when there is none (a pack built before W-19)
    or it is not JSON."""
    try:
        return json.loads((Path(pdir) / VERIFIED).read_bytes())
    except FileNotFoundError:
        return None
    except ValueError:
        return {}


def _record_matches(rec: dict, pid: bytes, release: str, manifest_bytes: bytes) -> bool:
    return (rec.get("pack_id") == pid.hex() and rec.get("release_id") == release and rec.get("problems") == []
            and rec.get("manifest_blake3") == blake3(manifest_bytes).hexdigest())


def _write_record(d: Path, pid: bytes, release: str, manifest_bytes: bytes) -> None:
    """Write the record into the pack directory `d`, under a temporary name renamed into place."""
    tmp = d / f".{VERIFIED}.tmp-{os.getpid()}"
    tmp.write_bytes(_record_bytes(pid, release, manifest_bytes))
    os.replace(tmp, d / VERIFIED)


def _add_record(d: Path, release: str) -> None:
    """Add the record to a pack that has just verified in full and has none (built before W-19): the pack's own
    files are not touched. A pack directory this user cannot write keeps no record; the build says so."""
    if (d / VERIFIED).exists():
        return
    data = (d / MANIFEST).read_bytes()
    try:
        _write_record(d, pack_id_of(data), release, data)
    except OSError as e:
        _log(f"release {release}: pack {d.name} verified, but its record could not be written ({e}); every search "
             "of it verifies it in full")
        return
    _log(f"release {release}: pack {d.name} verified; its record {VERIFIED} added (the pack's files are unchanged)")


def search_check(pdir: Path, release: str) -> tuple[list[str], bool]:
    """What a search checks of its pack (card W-19), and whether the build recorded that the pack verified:
    (problems, recorded). The checks: the manifest's hash is the directory's name (the pack id), the pack is of
    `release`, every file the manifest lists is there with the size the manifest gives, and a record, when there
    is one, names this manifest (its pack id and full BLAKE3). It reads the manifest and the record and stats the
    files; it reads no shard and writes nothing. `recorded` is False for a pack built before W-19 (no record): the
    caller verifies it in full (`verify_dir`). `pack verify` stays the full check."""
    pdir = Path(pdir)
    data = (pdir / MANIFEST).read_bytes()
    manifest, pid = json.loads(data), pack_id_of(data)
    problems = []
    if pdir.name != pid.hex():
        problems.append(f"{pdir}: the directory is not named by the manifest hash {pid.hex()}")
    if manifest.get("release_id") != release:
        problems.append(f"{pdir}: a pack of release {manifest.get('release_id')}, not {release}")
    for entry in manifest["shards"] + manifest["meta"]:
        try:
            size = (pdir / entry["path"]).stat().st_size
        except FileNotFoundError:
            problems.append(f"{entry['path']}: listed but missing")
            continue
        if size != entry["bytes"]:
            problems.append(f"{entry['path']}: {size} bytes, manifest says {entry['bytes']}")
    rec = read_record(pdir)
    recorded = rec is not None and _record_matches(rec, pid, release, data)
    if rec is not None and not recorded:
        problems.append(f"{VERIFIED}: does not name this manifest (pack id {pid.hex()})")
    return problems, recorded


def _count_records(path: Path) -> int:
    """The FASTA records of a zstd shard, counted while decompressing (a header is a line that starts with `>`)."""
    n, at_line_start = 0, True
    with pa.CompressedInputStream(pa.OSFile(str(path)), "zstd") as fh:
        while chunk := fh.read(1 << 24):
            n += chunk.count(b"\n>") + (1 if at_line_start and chunk[:1] == b">" else 0)
            at_line_start = chunk.endswith(b"\n")
    return n


def _count_source(rf: ReleaseFiles, manifest: dict, target: str, out: Path, scratch: Path | None = None) -> int | None:
    gcv = manifest["gene_call_version"]
    if target == "unique":
        # distinct protein uids through a sorter, in a scratch directory under `scratch` (a search's: its work
        # directory, card W-19; `pack verify`'s: the profile's work directory, card W-22), else under the system
        # temporary directory; never beside the pack, so nothing but a build writes into the pack directory (card
        # W-22); removed after
        base = Path(scratch) if scratch is not None else Path(tempfile.gettempdir())
        uids = _Sorter(base / f".verify-{os.getpid()}-{os.urandom(4).hex()}",
                       pa.schema([pa.field("protein_uid", UID)]), None, None)
        try:
            for f in _protein_files(rf, gcv):
                for b in _batches(f, ["protein_uid"], BATCH_ROWS):
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
           packs: Path, scratch: Path | None = None) -> list[str]:
    """C-12: problems with the release's pack in the pack directory `packs` (the one of `pack_id`, when
    given; `find`); empty means it verifies. Sequence counts are checked against the release's files under
    `root` (`release_manifest`, as for `build`). Nothing is written into `packs` (card W-22): the count of
    distinct proteins sorts under `scratch` (`verify_dir`)."""
    try:
        out = find(packs, release, pack_id, root=root)
    except (FileNotFoundError, ValueError) as e:
        return [str(e)]
    return verify_dir(root, release, out, pack_id, min_shard_bytes, release_manifest, scratch=scratch)


def verify_dir(root: Path, release: str, out: Path, pack_id: bytes | None = None,
               min_shard_bytes: int = layout.TARGET_FILE_BYTES, release_manifest: dict | None = None, *,
               scratch: Path | None = None) -> list[str]:
    """`verify` of the pack in the directory `out`: also, the directory is named by the manifest's hash, and the
    build's verification record, when there is one, names this manifest. The distinct proteins are counted through
    a sorter under `scratch` (default: the system temporary directory; never beside the pack, card W-22)."""
    root, out = Path(root).resolve(), Path(out)
    if not (out / MANIFEST).exists():
        return [f"{out}: no {MANIFEST}"]
    try:
        rf = ReleaseFiles.of(root, release, release_manifest)
    except FileNotFoundError as e:
        return [f"release {release}: {e}"]
    data = (out / MANIFEST).read_bytes()
    manifest, pid = json.loads(data), pack_id_of(data)
    problems = []
    if pack_id is not None and pack_id != pid:
        problems.append(f"pack_id {pack_id.hex()} differs from the manifest hash {pid.hex()}")
    if out.name != pid.hex():
        problems.append(f"{out}: the directory is not named by the manifest hash {pid.hex()}")
    rec = read_record(out)
    if rec is not None and not _record_matches(rec, pid, release, data):
        problems.append(f"{VERIFIED}: does not name this manifest (pack id {pid.hex()})")
    return problems + _problems(rf, release, out, manifest, min_shard_bytes, scratch)


_SHARD_PART_RE = re.compile(r"^(.*(?:^|/))part-(\d+)(\.[^/]*)?$")


def _shard_sort_key(path: str) -> tuple[int, str, int, str]:
    """Numeric sort key for `.../part-<k>.<ext>` shard paths (card W-41); deterministic fallback for other paths."""
    m = _SHARD_PART_RE.match(path)
    if m:
        return (0, m.group(1), int(m.group(2)), m.group(3) or "")
    return (1, path, 0, "")


def _problems(rf: ReleaseFiles, release: str, out: Path, manifest: dict, min_shard_bytes: int,
              scratch: Path | None) -> list[str]:
    """The checks of `verify_dir` on the pack's contents, whatever its directory is called (the build runs them on
    its temporary directory before renaming it: card W-19)."""
    problems = []
    if manifest.get("release_id") != release:
        problems.append(f"{out}: a pack of release {manifest.get('release_id')}, not {release}")
    listed = {MANIFEST, VERIFIED}
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
        for s in (sorted(shards, key=lambda s: _shard_sort_key(s["path"]))[:-1] if target != "models" else []):
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
        expect = _count_source(rf, manifest, target, out, scratch)
        if expect is not None and expect != info["n_seq"]:
            problems.append(f"{target}: {info['n_seq']} sequences, its source has {expect}")
    return problems


# --- the estimate (card W-35) --------------------------------------------------------------------------------------

# CPU-minutes per unit of each process's work (card W-35; module docstring), fitted to the run records of the pack of
# `rN.3` built in 4 parts (measured, wilkcomp02, 2026-11-14, unit `w35-rn3`, each part's `done.json`;
# planning/measurements/2026-11-14_pack_build_in_parts.md; 10,634,028 protein_seq rows, 12,276,916,069 genome bases,
# one library of 30,134 models) and rounded up, so that no part's estimate is below its record
# (`tests/r2/test_w35_pack_parts.py` holds the records and checks each estimate within [record, 1.2 x record]):
# - unique: 2.58 CPU-minutes for the 10,634,028 rows (0.238 per million, derived); 0.25.
# - igr: 1.97 to 2.13 for each quarter of the bases, about 3.07 Gbp (0.68 per Gbp at the slowest, derived); 0.70.
# - dark_orf: 15.69 to 16.65 for each quarter (5.41 per Gbp at the slowest, derived); 5.45.
# - the last process: 13.19, of which the models about 2.6 (160 s of wall time in its log), the rest split by its log
#   into the unique shards and half the verify (about 155 s: 0.24 per million rows, derived) and the region shards, the
#   manifest's hashes and half the verify (about 512 s: 0.70 per Gbp, derived); 0.25 and 0.70.
# - the driver: 0.95 (the unit's 89.68 CPU-minutes less the parts' 88.73, derived); 0.075 per Gbp.
# The rates of the one-process estimate are the same: 93.0 for rN.3 against its 79.5 measured in one process (+17%),
# 90.4 for rN.1 (no library) against 86.0 (+5%) (derived; measured in planning/measurements/2026-11-11_build_r50.md and
# 2026-11-11_pack_build_rN1.md). Before the measurement the rates were derived from the one-process log of rN.3
# (dark_orf 4.81 per Gbp: 11% under the slowest part's record; the last process 0.40 per million rows and 0.80 per
# Gbp: 27% over its record).
UNIQUE_CPU_MIN_PER_M_ROWS = 0.25
IGR_CPU_MIN_PER_GBP = 0.70
DARK_ORF_CPU_MIN_PER_GBP = 5.45
ASSEMBLE_CPU_MIN_PER_M_ROWS = 0.25
ASSEMBLE_CPU_MIN_PER_GBP = 0.70
MODELS_CPU_MIN_PER_LIBRARY = 2.6  # Pfam-A 38.2, 30,134 models (rN.3's one library)
DRIVER_CPU_MIN_PER_GBP = 0.075
PART_START_CPU_MIN = 0.05  # a process's start: imports, opening the release (estimate)
# the budget without a host profile: AGENTS §7.2's 300 CPU-minutes per process
DEFAULT_CPU_MIN_PER_PROCESS = 300


@dataclass(frozen=True)
class Counts:
    """What the estimate is made from: the release's `protein_seq` rows (a protein stored by two stages counts twice),
    its genome bases and contigs (the `contig_len` of its tapes, once per contig) and the domain libraries packed."""
    proteins: int
    bases: int
    contigs: int
    libraries: int = 0


def _margin() -> float:
    from biobox.stage import SHARD_MARGIN

    return SHARD_MARGIN


def _counts(rf: ReleaseFiles, gcv: str, tv: str, libraries: int) -> Counts:
    proteins = _n_rows(_protein_files(rf, gcv))
    bases = contigs = 0
    if rf.table("seq_locator"):  # a protein-only root has no nucleotide targets
        for _, (length, _, _) in _contig_tapes(_tape_files(rf, tv), False):
            bases += length
            contigs += 1
    return Counts(proteins, bases, contigs, libraries)


def release_counts(root: Path, release: str, gene_call_version: str = GENE_CALL_VERSION,
                   producers: list[str] | None = None, release_manifest: dict | None = None,
                   libraries: int = 0) -> Counts:
    """The release's counts for the estimate (`Counts`), read from the files' metadata and the tapes' lengths, one batch
    of contigs at a time."""
    rf = ReleaseFiles.of(Path(root).resolve(), release, release_manifest)
    gcv = gene_call_version
    return _counts(rf, gcv, tape_version(sorted(producers) if producers else [gcv]), libraries)


def part_cpu_min(target: str, *, rows: int = 0, bases: int = 0, libraries: int = 0) -> float:
    """The estimate of one process of the build: a part of `target` (`unique` over `rows` protein_seq rows, `igr` or
    `dark_orf` over `bases` genome bases), the last process (`assemble`) or the driver (`driver`) of a release of
    `rows` rows, `bases` bases and `libraries` domain libraries (rates above)."""
    m, g = rows / 1e6, bases / 1e9
    if target == "unique":
        return PART_START_CPU_MIN + UNIQUE_CPU_MIN_PER_M_ROWS * m
    if target == "igr":
        return PART_START_CPU_MIN + IGR_CPU_MIN_PER_GBP * g
    if target == "dark_orf":
        return PART_START_CPU_MIN + DARK_ORF_CPU_MIN_PER_GBP * g
    if target == "assemble":
        return (PART_START_CPU_MIN + ASSEMBLE_CPU_MIN_PER_M_ROWS * m + ASSEMBLE_CPU_MIN_PER_GBP * g
                + MODELS_CPU_MIN_PER_LIBRARY * libraries)
    if target == "driver":
        return PART_START_CPU_MIN + DRIVER_CPU_MIN_PER_GBP * g
    raise ValueError(f"no estimate for {target!r}")


def one_process_cpu_min(c: Counts) -> float:
    """The estimate of the build in one process: the work of every part, the last process and the driver."""
    return (part_cpu_min("unique", rows=c.proteins) + part_cpu_min("igr", bases=c.bases)
            + part_cpu_min("dark_orf", bases=c.bases)
            + part_cpu_min("assemble", rows=c.proteins, bases=c.bases, libraries=c.libraries)
            + part_cpu_min("driver", bases=c.bases) - 4 * PART_START_CPU_MIN)


def _rates() -> str:
    return (f"{UNIQUE_CPU_MIN_PER_M_ROWS} CPU-minutes per million protein_seq rows for unique, {IGR_CPU_MIN_PER_GBP} and "
            f"{DARK_ORF_CPU_MIN_PER_GBP} per Gbp for igr and dark_orf, {ASSEMBLE_CPU_MIN_PER_M_ROWS} per million rows, "
            f"{ASSEMBLE_CPU_MIN_PER_GBP} per Gbp and {MODELS_CPU_MIN_PER_LIBRARY} per library for the last process: "
            "the rates in biobox.search.pack, card W-35")


def estimate_line(release: str, c: Counts, budget: int) -> str:
    """The estimate as `pack build` prints it before any work."""
    return (f"release {release}: estimate {one_process_cpu_min(c):,.1f} CPU-minutes in one process for {c.bases:,} "
            f"genome bases in {c.contigs:,} contigs, {c.proteins:,} protein_seq rows and {c.libraries} domain "
            f"libraries ({_rates()}); the budget of a process is {_margin() * budget:g} CPU-minutes "
            f"({_margin():.0%} of cpu_min_per_process {budget}, AGENTS §7.2)")


def parts_needed(c: Counts, budget: int) -> int:
    """The fewest parts whose even share of the bases puts each region part within `SHARD_MARGIN` of `budget`."""
    room = _margin() * budget - PART_START_CPU_MIN
    rate = max(IGR_CPU_MIN_PER_GBP, DARK_ORF_CPU_MIN_PER_GBP) * c.bases / 1e9
    return max(1, math.ceil(rate / room)) if room > 0 else 1


def one_process_refusal(c: Counts, budget: int) -> str | None:
    """Why the build may not run in one process, or None: its estimate is over `SHARD_MARGIN` of `budget`, the
    profile's `cpu_min_per_process` (as `stage` refuses a step with no shards)."""
    est, lim = one_process_cpu_min(c), _margin() * budget
    if est <= lim:
        return None
    return (f"the pack of {c.bases:,} genome bases and {c.proteins:,} protein_seq rows is estimated at {est:,.1f} "
            f"CPU-minutes in one process ({_rates()}), over {_margin():.0%} of the profile's cpu_min_per_process of "
            f"{budget} (AGENTS §7.2): build it in parts, with --parts {parts_needed(c, budget)} or more and --jobs J "
            "(bio-box pack build --help)")


def parts_refusal(parts: list[dict], budget: int, n: int) -> str | None:
    """Why the build in `n` parts may not run, or None: a process estimated over `SHARD_MARGIN` of `budget`. A region
    part over it needs more parts (or holds one contig too large for any); `unique` and the last process cannot be cut."""
    lim = _margin() * budget
    over = [p for p in parts if p["estimate_cpu_min"] > lim]
    if not over:
        return None
    say = ", ".join(f"{p['name']} {p['estimate_cpu_min']:,.1f} CPU-minutes ({p['contigs']:,} contigs, {p['bases']:,} "
                    f"bases, {p.get('rows', 0):,} rows)" for p in over)
    cut = [p for p in over if p["target"] in ("igr", "dark_orf")]
    if len(cut) < len(over):
        fix = "unique and the last process are not cut into parts: split the release"
    elif all(p["contigs"] == 1 for p in cut):
        fix = "each holds one contig, which no number of parts cuts"
    else:
        worst = max(p["estimate_cpu_min"] for p in cut)
        fix = f"give more parts (--parts {max(n + 1, math.ceil(n * worst / lim))} or more)"
    return (f"in {n} parts, over {_margin():.0%} of the profile's cpu_min_per_process of {budget} (AGENTS §7.2): "
            f"{say}; {fix}")


# --- the build in parts (card W-35; module docstring) ----------------------------------------------------------------

_MODULE = "biobox.search.pack"  # what `python -m` runs for a part and for the last process
PARTS_WORK = "scratch/pack"  # the default work base, under the profile's tmp
RUN, DONE, ASSEMBLED = "run.parquet", "done.json", "assembled.json"
REGION_TARGETS = ("igr", "dark_orf")


def _cut(rf: ReleaseFiles, tv: str, n: int) -> list[dict]:
    """The release's contigs (its tapes, in uid order) cut into `n` consecutive ranges: contig c starts range k at the
    first contig whose midpoint, in genome bases from the first contig, is at k/n of the release's bases or past it
    (`2 n before + n len >= 2 k total`, exact in integers). Each range as {lo, hi (uid hex, None unbounded), empty,
    contigs, bases}; a range no contig starts holds none. Two passes over the tapes' lengths, one batch at a time."""
    files = _tape_files(rf, tv) if rf.table("seq_locator") else []
    total = sum(length for _, (length, _, _) in _contig_tapes(files, False))
    parts = [{"lo": None, "hi": None, "empty": False, "contigs": 0, "bases": 0} for _ in range(n)]
    k, before = 0, 0
    for cu, (length, _, _) in _contig_tapes(files, False):
        while k + 1 < n and 2 * n * before + n * length >= 2 * (k + 1) * total:
            k += 1
            parts[k]["lo"] = parts[k - 1]["hi"] = cu.hex()
        parts[k]["contigs"] += 1
        parts[k]["bases"] += length
        before += length
    for p in parts[k + 1:]:
        p["empty"] = True
    return parts


def _done(wdir: Path, part: dict) -> dict | None:
    """The record of a part whose output is in the work directory `wdir` and is of its range (its `done.json` and its
    run, of the size the record gives), else None."""
    try:
        rec = json.loads((wdir / part["name"] / DONE).read_text())
    except (FileNotFoundError, ValueError):
        return None
    run = wdir / part["name"] / RUN
    if any(rec.get(k) != part[k] for k in ("lo", "hi", "empty", "target")) or not run.is_file() \
            or run.stat().st_size != rec.get("run_bytes"):
        return None
    return rec


def _usage(t0: float) -> dict:
    """This process's CPU-minutes, wall seconds since `t0` and peak RSS (KiB: `ru_maxrss` on Linux)."""
    ru = resource.getrusage(resource.RUSAGE_SELF)
    return {"cpu_min": round((ru.ru_utime + ru.ru_stime) / 60, 3), "wall_s": round(time.monotonic() - t0, 1),
            "peak_rss_kib": ru.ru_maxrss}


def _usage_line(rec: dict, est: float) -> str:
    return (f"{rec['cpu_min']:,.2f} CPU-minutes (estimate {est:,.2f}), wall {rec['wall_s']:,.1f} s, peak RSS "
            f"{rec['peak_rss_kib'] / 1024:,.1f} MiB")


def _allowance(prof) -> Callable[[int], int]:
    """The threads work in flight may total now, for a launch of `threads` (`runstep.allowance`, read again at each
    call, as `build-stage` reads it before each shard): said once and again when it changes."""
    last: list = [None]

    def allow(threads: int) -> int:
        a = runstep.allowance(prof, datetime.now(), runstep.read_cpu(), threads=threads)
        if a.threads != last[0]:
            _log(a.line)
            last[0] = a.threads
        return a.threads

    return allow


def _part_runner(cap: dict) -> Callable[..., int]:
    """Run a part (`stage._launch`'s runner) as its own process in this process's group, under the memory cap `cap`
    (`runstep.memory_cap`: the one inherited, or a `MemoryMax` scope it opens) at the profile's `nice_monolithic`,
    with `threads` in the thread variables, as `build-stage` runs a shard rescan. A signal gives 128 + its number."""

    def run(prof, root: Path, stage: str, step: str, cmd, *, threads: int) -> int:
        argv = runstep.build_argv(list(cmd), prof.nice_monolithic, prof.mem_scope,
                                  inherited=cap if cap["how"] == "inherited" else None)
        env = dict(os.environ, BIOBOX_THREADS=str(threads), **{v: str(threads) for v in runstep.THREAD_VARS})
        code = subprocess.run(argv, env=env, check=False).returncode
        return code if code >= 0 else 128 - code

    return run


def build_in_parts(root: Path, release: str, gene_call_version: str = GENE_CALL_VERSION,
                   models: dict[str, Path] | None = None, shard_text_bytes: int = SHARD_TEXT_BYTES,
                   producers: list[str] | None = None, release_manifest: dict | None = None, *, packs: Path, prof,
                   parts: int, jobs: int | None = None, work: Path | None = None, scratch: Path | None = None,
                   run: Callable[..., int] | None = None,
                   allow: Callable[[int], int] | None = None) -> tuple[int, dict | None]:
    """Build the release's pack in parts (card W-35, module docstring), each its own process: `unique`, `parts`
    ranges of contigs for each of `igr` and `dark_orf`, then the last process; (0, the manifest with `pack_id`), or
    (the exit code of a part that failed or was stopped, or 75 for `STOP`, None). The pack is that of `build` with the
    same arguments, byte for byte.

    `prof` is the host profile: its `cpu_min_per_process` and `SHARD_MARGIN` bound every process's estimate
    (ValueError before any work otherwise: `parts_refusal`), its memory cap and nice are each part's, its tmp holds the
    work base (`<tmp>/scratch/pack`, or `work`). At most `jobs` parts run at once (default: the allowance at the start),
    one thread each, within `allow` (default: `runstep.allowance` read before each launch). `run` (tests) stands in for
    the runner of a part, with `stage._launch`'s signature. `scratch` is where the full verify of a pack already there
    sorts (as for `build`). ValueError, before any part, for a work directory under the root or on another file system
    than the pack directory, or one another driver holds."""
    root = Path(root).resolve()
    packs = outside_root(root, packs, "the pack directory")
    if isinstance(parts, bool) or not isinstance(parts, int) or parts < 1:
        raise ValueError(f"--parts {parts}: must be a whole number, at least 1")
    if jobs is not None and (isinstance(jobs, bool) or jobs < 1):
        raise ValueError(f"--jobs {jobs}: must be at least 1")
    base_work = outside_root(root, Path(work) if work is not None else Path(prof.tmp) / PARTS_WORK,
                             "the work directory")
    with _stop_on_signals():
        return _build_in_parts(root, release, gene_call_version, models, shard_text_bytes, producers, release_manifest,
                               packs, prof, parts, jobs, base_work, scratch, run, allow)


def _build_in_parts(root: Path, release: str, gcv: str, models: dict[str, Path] | None, shard_text_bytes: int,
                    producers: list[str] | None, release_manifest: dict | None, packs: Path, prof, n: int,
                    jobs: int | None, base_work: Path, verify_scratch: Path | None, run, allow
                    ) -> tuple[int, dict | None]:
    from biobox import stage

    rf = ReleaseFiles.of(root, release, release_manifest)
    producers = sorted(producers) if producers is not None else [gcv]
    tv, iv = tape_version(producers), igr_version(producers)
    models = {lib: Path(p).resolve() for lib, p in sorted((models or {}).items())}
    # the estimate of every process, before any work (card W-35)
    cut = _cut(rf, tv, n)
    c = Counts(_n_rows(_protein_files(rf, gcv)), sum(p["bases"] for p in cut), sum(p["contigs"] for p in cut),
               len(models))
    plan_parts = [{"name": "unique", "target": "unique", "lo": None, "hi": None, "empty": False, "contigs": c.contigs,
                   "bases": c.bases, "rows": c.proteins}]
    for t in REGION_TARGETS:
        plan_parts += [{**p, "name": f"{t}-{k}-of-{n}", "target": t, "rows": 0} for k, p in enumerate(cut)]
    for p in plan_parts:
        p["estimate_cpu_min"] = part_cpu_min(p["target"], rows=p["rows"], bases=0 if p["target"] == "unique" else
                                             p["bases"])
    last = {"name": "assemble", "target": "assemble", "contigs": c.contigs, "bases": c.bases, "rows": c.proteins,
            "estimate_cpu_min": part_cpu_min("assemble", rows=c.proteins, bases=c.bases, libraries=c.libraries)}
    budget = prof.cpu_min_per_process
    _log(f"release {release}: in {n} parts: " + "; ".join(
        f"{p['name']} {p['estimate_cpu_min']:,.2f} CPU-minutes ({p['contigs']:,} contigs, {p['bases']:,} bases)"
        for p in [*plan_parts, last]) + f"; budget {_margin() * budget:g} CPU-minutes per process")
    why = parts_refusal([*plan_parts, last], budget, n)
    if why is not None:
        raise ValueError(why)
    cap = None
    if run is None:
        try:
            cap = runstep.memory_cap(prof.mem_scope, require=prof.require_mem_cap)
        except runstep.NoCap as e:
            raise ValueError(f"no part can start: {e}") from None
        run = _part_runner(cap)
    _log(f"release {release}: checking genetic codes")
    _check_codes(rf, gcv)
    head, sources = _head_and_sources(root, rf, release, gcv, iv, tv, models, shard_text_bytes)
    want = _fingerprint({**head, "targets": {t: {"sources": s} for t, s in sources.items()}})
    wdir = base_work / f"release={release}" / blake3(json.dumps(want, sort_keys=True).encode()).hexdigest()[:32]
    have = _existing(root, release, packs, want, release_manifest, verify_scratch)
    if have is not None:
        return 0, have
    base = release_dir(packs, release)
    base.mkdir(parents=True, exist_ok=True)
    wdir.mkdir(parents=True, exist_ok=True)
    if os.stat(wdir).st_dev != os.stat(base).st_dev:
        raise ValueError(f"the work directory {wdir} is on another file system than the pack directory {packs}; the "
                         "pack is renamed from one into the other: give a --work on the pack directory's file system")
    lock = open(wdir / ".lock", "a")
    try:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError(f"another pack build of release {release} in parts holds {wdir / '.lock'}") from None
        for p in sorted(wdir.iterdir()):
            if p.name.startswith(".") and p.is_dir():  # a killed part's or last process's (SIGKILL): never taken
                shutil.rmtree(p)
                _log(f"release {release}: removed {p}, left by a process that was killed")
        plan = {"card": "W-35", "root": str(root), "release": release, "release_manifest": release_manifest,
                "gene_call_version": gcv, "igr_version": iv, "tape_version": tv, "packs": str(Path(packs).resolve()),
                "models": {lib: str(p) for lib, p in models.items()}, "head": head, "sources": sources,
                "parts": plan_parts, "assemble": last}
        plan_path = wdir / f"plan-of-{n}.json"
        tmp = wdir / f".{plan_path.name}.tmp-{os.getpid()}"
        tmp.write_text(json.dumps(plan, sort_keys=True, indent=1))
        os.replace(tmp, plan_path)
        todo = [p for p in plan_parts if _done(wdir, p) is None]
        for p in plan_parts:
            if p not in todo:
                _log(f"part {p['name']}: its output is there from an earlier call; not run again")
        allow = allow if allow is not None else _allowance(prof)
        jobs = jobs if jobs is not None else max(1, allow(1))
        _log(f"release {release}: {len(todo)} of {len(plan_parts)} parts to run, at most {jobs} at once, one thread "
             f"each; work directory {wdir}")
        steps = [stage.Step(release, p["name"], ("part", "--plan", str(plan_path), "--name", p["name"]), module=_MODULE)
                 for p in todo]
        code = stage._launch(prof, root, steps, jobs, 1, run, None, allow) if steps else None
        for p in plan_parts:
            rec = _done(wdir, p)
            if rec is not None:
                _log(f"part {p['name']} done: {rec['rows']:,} rows, {_usage_line(rec, p['estimate_cpu_min'])}"
                     + ("" if p in todo else " (an earlier call's)"))
        if code is not None:
            _log(f"release {release}: stopped before the last process (exit {code}); the parts with output stay in "
                 f"{wdir}, and the next call runs only the others")
            return code, None
        code = stage._launch(prof, root, [stage.Step(release, "assemble", ("assemble", "--plan", str(plan_path)),
                                                     module=_MODULE)], 1, 1, run, None, allow)
        if code is not None:
            _log(f"release {release}: the last process did not finish (exit {code}); the parts stay in {wdir}")
            return code, None
        rec = json.loads((wdir / ASSEMBLED).read_text())
        _log(f"part assemble done: pack {rec['pack_id']}, {_usage_line(rec, last['estimate_cpu_min'])}")
        man, pid = read_manifest(pack_dir(packs, release, rec["pack_id"]))
        shutil.rmtree(wdir)  # the parts' runs, the plan and the lock: scratch, once the pack is in place
    finally:
        lock.close()
    _log(f"release {release}: pack {pid.hex()} built in {n} parts; work directory {wdir} removed")
    return 0, {**man, "pack_id": pid.hex()}


def _run_part(plan: dict, wdir: Path, name: str) -> None:
    """One part (module docstring): its target's rows over its range, sorted into one run, written with its record into
    a temporary directory renamed to `<wdir>/<name>` at the end."""
    t0 = time.monotonic()
    part = next((p for p in plan["parts"] if p["name"] == name), None)
    if part is None:
        raise ValueError(f"no part {name} in the plan")
    if _done(wdir, part) is not None:
        _log("its output is there already; nothing to do")
        return
    tmp = wdir / f".{name}.tmp-{os.getpid()}-{os.urandom(4).hex()}"
    tmp.mkdir()
    try:
        root, target = Path(plan["root"]), part["target"]
        rf = ReleaseFiles.of(root, plan["release"], plan["release_manifest"])
        gcv, rng = plan["gene_call_version"], _Range.of(part)
        _log(f"{target}: start, {part['contigs']:,} contigs, {part['bases']:,} bases (estimate "
             f"{part['estimate_cpu_min']:,.2f} CPU-minutes), in {tmp}")
        if target == "unique":
            schema = _PROTEIN_SORT
            sorter = _Sorter(tmp / "sort", schema, itemgetter(0), "unique")
            _unique_rows(rf, gcv, sorter)
        else:
            schema = _REGION_SORT
            sorter = _Sorter(tmp / "sort", schema, itemgetter(0, 1), target)
            store = SeqStore(root, rf) if rf.table("seq_locator") else None
            try:
                if target == "igr":
                    _igr_regions(rf, plan["igr_version"], plan["tape_version"], store, sorter, rng)
                else:
                    _dark_orf_regions(rf, gcv, plan["tape_version"], store, tmp, sorter, part["contigs"], rng)
            finally:
                if store is not None:
                    store.close()
        try:
            _log(f"{target}: sorting {sorter.n:,} rows into one run")
            rows = _write_run(tmp / RUN, schema, sorter.rows())
        finally:
            sorter.close()
        rec = {**{k: part[k] for k in ("name", "target", "lo", "hi", "empty", "contigs", "bases")}, "rows": rows,
               "run_bytes": (tmp / RUN).stat().st_size, "estimate_cpu_min": part["estimate_cpu_min"], **_usage(t0)}
        (tmp / DONE).write_text(json.dumps(rec, sort_keys=True, indent=1))
        os.rename(tmp, wdir / name)
    except BaseException as e:
        shutil.rmtree(tmp, ignore_errors=True)
        if isinstance(e, (Stopped, KeyboardInterrupt)):
            why = "SIGINT" if isinstance(e, KeyboardInterrupt) else f"a signal (exit {e.code})"
            _log(f"stopped by {why}; removed {tmp}; the next call runs this part again")
        raise
    _log(f"done: {rows:,} rows, {_usage_line(rec, part['estimate_cpu_min'])}")


def _assemble(plan: dict, wdir: Path) -> None:
    """The last process (module docstring): every target's runs merged in part order into the shards and meta files,
    the models, the manifest; the pack verified, recorded and renamed into the pack directory."""
    t0 = time.monotonic()
    root, release, packs = Path(plan["root"]), plan["release"], Path(plan["packs"])
    rf = ReleaseFiles.of(root, release, plan["release_manifest"])
    by = {t: [p for p in plan["parts"] if p["target"] == t] for t in ("unique", *REGION_TARGETS)}
    recs = {p["name"]: _done(wdir, p) for ps in by.values() for p in ps}
    missing = [k for k, v in recs.items() if v is None]
    if missing:
        raise ValueError(f"parts without output in {wdir}: {', '.join(missing)}")
    shard_bytes = plan["head"]["rules"]["shard_text_bytes"]
    tag = f"{os.getpid()}-{os.urandom(4).hex()}"
    out, scratch = wdir / f".building-{tag}", wdir / f".scratch-{tag}"
    out.mkdir()

    def merged(target: str, schema: pa.Schema, key: Callable) -> Iterator[tuple]:
        ps = by[target]
        _log(f"{target}: merging {len(ps)} part runs, {sum(recs[p['name']]['rows'] for p in ps):,} rows")
        return heapq.merge(*(_read_run(wdir / p["name"] / RUN, schema) for p in ps), key=key)

    try:
        scratch.mkdir()
        _log(f"release {release}: building in {out} (scratch {scratch})")
        built = {
            "unique": _write_unique(merged("unique", _PROTEIN_SORT, itemgetter(0)), out, shard_bytes),
            "igr": _write_region_rows(merged("igr", _REGION_SORT, itemgetter(0, 1)), out, "igr", ".fna.zst",
                                      shard_bytes),
            "dark_orf": _write_region_rows(merged("dark_orf", _REGION_SORT, itemgetter(0, 1)), out, "dark_orf",
                                           ".faa.zst", shard_bytes),
            "models": _build_models(out, {lib: Path(p) for lib, p in plan["models"].items()}),
        }
        manifest, data = _write_manifest(out, plan["head"], plan["sources"], built)
        final, placed = _place(root, rf, release, out, manifest, data, scratch, packs, plan["release_manifest"], wdir)
        rec = {"pack_id": final.name, "placed": placed, **_usage(t0)}
        tmp = wdir / f".{ASSEMBLED}.tmp-{os.getpid()}"
        tmp.write_text(json.dumps(rec, sort_keys=True, indent=1))
        os.replace(tmp, wdir / ASSEMBLED)
    except BaseException as e:
        shutil.rmtree(out, ignore_errors=True)
        shutil.rmtree(scratch, ignore_errors=True)
        if isinstance(e, (Stopped, KeyboardInterrupt)):
            why = "SIGINT" if isinstance(e, KeyboardInterrupt) else f"a signal (exit {e.code})"
            _log(f"stopped by {why}; removed {out} and {scratch}; the parts stay, and the next call runs this again")
        raise
    shutil.rmtree(scratch, ignore_errors=True)
    if placed:
        _log(f"release {release}: pack {final.name} written: {len(manifest['shards']) + len(manifest['meta']) + 2} "
             f"files, {sum(s['bytes'] for s in manifest['shards'] + manifest['meta']):,} bytes in shards and meta")
    _log(f"done: {_usage_line(rec, plan['assemble']['estimate_cpu_min'])}")


def main(argv: list[str] | None = None) -> int:
    """`python -m biobox.search.pack part --plan P --name N` (one part) or `... assemble --plan P` (the last process):
    the processes `bio-box pack build --parts` launches (card W-35). 0, or 128 + the signal that stopped it."""
    ap = argparse.ArgumentParser(prog="python -m biobox.search.pack",
                                 description="one process of `bio-box pack build --parts` (card W-35)")
    ap.add_argument("action", choices=["part", "assemble"])
    ap.add_argument("--plan", required=True, help="the plan the driver wrote, <work>/plan-of-<N>.json")
    ap.add_argument("--name", help="part: its name in the plan")
    a = ap.parse_args(argv)
    if a.action == "part" and not a.name:
        ap.error("part needs --name")
    plan_path = Path(a.plan)
    plan = json.loads(plan_path.read_text())
    before = getattr(_WHO, "name", None)
    _WHO.name = f"pack part {a.name if a.action == 'part' else 'assemble'}"
    try:
        with _stop_on_signals():
            if a.action == "part":
                _run_part(plan, plan_path.parent, a.name)
            else:
                _assemble(plan, plan_path.parent)
    except Stopped as e:
        return e.code
    except KeyboardInterrupt:
        return 130
    finally:
        if before is None:
            del _WHO.name
        else:
            _WHO.name = before
    return 0


if __name__ == "__main__":
    sys.exit(main())
