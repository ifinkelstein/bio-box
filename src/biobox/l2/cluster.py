"""Protein clustering with MMseqs2 (operation O-CLUSTER).

NR90 clusters all distinct proteins of one gene-calling version; NR50 clusters
the NR90 representatives, and `mergeclusters` maps every protein to its NR50
representative. Each member is then re-aligned to its representative at each
level (`mmseqs align`, thresholds off) so `cluster_member_aln` has one row per
protein per level, representatives included as self rows.

Writes under `<root>/`:

    l2/cluster_member/release=<R>/part-<k>.parquet
    l2/cluster_member_aln/release=<R>/part-<k>.parquet
    l2/member_cluster/release=<R>/part-<k>.parquet
    ledger/release=<R>/step=cluster/part-0.parquet

Two ways to cluster a release, with byte-identical files (card B-41): `run`, the two `cluster`
workflow calls each one process; or in parts (`plan`, `step`, `finish`, `run_split`; section
"clustering in parts" below), each module call its own process and every per-query one in parts,
which also writes its steps' records, `l2/cluster_plan/release=<R>/<step>.json`.

The tables go through `biobox.io.write_parts`, cut at `max_file_bytes` of estimated size (card
B-19). `member_cluster` (card B-27) is `cluster_member` sorted by (level, member_uid), for readers
that look a protein's clusters up by member (`cluster_member` is sorted by cluster, so such a lookup
cannot prune it): a member reads one row group per level, and a lookup at one level reads only that
level's rows; `write_member_cluster` makes it from the `cluster_member` files just written, with
one DuckDB sort that spills.

The proteins are those of the release's stages only (card B-25), so a stage written later or
while the release is built never enters its clusters.

Inputs are written to FASTA sorted by protein_uid and MMseqs2 runs with
`--shuffle 0` and explicit `--threads`, so output does not depend on input order
or on the host's core count.

Bounded memory (card B-24; AGENTS §4.1). No list or dict of proteins, members or alignments is
held; the Python side holds a batch of rows and one cluster's members at a time.

- Each stage's `protein_seq` files are one run sorted by protein_uid; the runs are merged into the
  FASTA `SCRATCH_ROWS` proteins at a time, so the sequences are never sorted or held.
- MMseqs2's two TSV outputs per level are converted to flat Parquet a block at a time; an alignment
  row carries its line number.
- Every other step is a plain DuckDB sort of those flat files (`DUCK_MEMORY_LIMIT`, one thread,
  spilling into the private working directory; no join, window or hash aggregate, whose state does
  not spill as a sort does) followed by a merge in Python: the members sorted by uid against the
  proteins (every protein is a member exactly once), and against the alignment lines sorted by
  (member, line), which keeps each member's first line whose query is its representative; the
  members sorted by (representative, member), whose cluster uid is hashed one cluster at a time.
  The rows found go to flat Parquet, and one more sort per table gives the writer its sort order.

The MMseqs2 command lines, parameters and thread counts are those of the builder before card B-24,
so cluster membership is unchanged. Before the card the Python side held every row in lists: 9.57 GB
at 2.05 M proteins (measured, planning/measurements/2026-11-10_cluster_memory.md).
"""

from __future__ import annotations

import bisect
import hashlib
import heapq
import itertools
import json
import math
import mmap
import os
import re
import shlex
import shutil
import subprocess
import tempfile
from collections.abc import Callable, Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq
from blake3 import blake3

from biobox import catalog, ids, io, layout, ledger, provenance
from biobox.registry.schemas import ARTIFACTS, table

# Estimated zstd bytes per row; they only cut parts. Measured on the golden set 36.5 (`cluster_member`)
# and 26.2 (`cluster_member_aln`) bytes per protein and level, 31.0 and 31.0 on the S1 pilot
# (planning/measurements/2026-11-09_design_stress_test.md §3, §6a).
EST_MEMBER_ROW = 37
EST_ALN_ROW = 31
# member_cluster: 40.5 bytes per row, measured on wilkcomp02, 2026-11-10, rL.3's 21,241,368 rows
# (planning/measurements/2026-11-10_label_and_membership_lookups.md)
EST_MEMBER_CLUSTER_ROW = 40

# UniRef-style: member covered >= 80% by the alignment to its representative (cov-mode 1,
# coverage of the target = member), greedy by length (cluster-mode 2, CD-HIT-like).
LEVELS = {
    "nr90": {"min_seq_id": 0.9, "coverage": 0.8},
    "nr50": {"min_seq_id": 0.5, "coverage": 0.8},
}
COMMON = ["--cov-mode", "1", "--cluster-mode", "2"]

DUCK_MEMORY_LIMIT = "2GB"  # the sorting connection's cap, as genecall's and postings'; larger sorts spill
SCRATCH_ROWS = 16_384  # rows per batch read, written or handed on by the Python side
_CSV_BLOCK = 1 << 22  # bytes per block when converting an MMseqs2 TSV
_ALN_COLUMNS = {"query": pa.string(), "target": pa.string(), "qstart": pa.int64(), "qend": pa.int64(),
                "tstart": pa.int64(), "tend": pa.int64(), "fident": pa.float64(), "alnlen": pa.int64(),
                "cigar": pa.string()}
_ALN_FORMAT = ",".join(_ALN_COLUMNS)


def params_hash(level: str) -> str:
    p = {"level": level, **LEVELS[level], "common": COMMON, "tool": "mmseqs2"}
    return hashlib.sha256(json.dumps(p, sort_keys=True).encode()).hexdigest()[:16]


def cluster_uid(member_uids: list[bytes], phash: str) -> bytes:
    return blake3(b"".join(sorted(member_uids)) + phash.encode()).digest(length=16)


class MMseqs:
    def __init__(self, workdir: Path, threads: int, exe: str = "mmseqs"):
        if shutil.which(exe) is None:
            raise FileNotFoundError(f"{exe} not on PATH")
        self.exe, self.threads, self.workdir = exe, threads, Path(workdir)
        self.log = open(self.workdir / "mmseqs.log", "a")

    def __call__(self, *args: str, threaded: bool = True) -> None:
        cmd = [self.exe, *args, "-v", "1"]
        if threaded:
            cmd += ["--threads", str(self.threads)]
        self.log.write(" ".join(cmd) + "\n")
        self.log.flush()
        subprocess.run(cmd, check=True, stdout=self.log, stderr=subprocess.STDOUT, cwd=self.workdir)


def _protein_files(root: Path, gene_call_version: str, stages: list[str]) -> list[list[Path]]:
    """The `protein_seq` files of each of `stages`, in bucket and part order (card B-25: no other stage
    is read; a listed stage without `protein_seq` files raises FileNotFoundError)."""
    out = []
    for s in stages:
        found = io.part_files(root / f"l0seq/protein_seq/gene_call_version={gene_call_version}/stage={s}",
                              "protein_prefix=*/part-*.parquet")
        if not found:
            raise FileNotFoundError(f"stage {s}: no protein_seq files for gene_call_version={gene_call_version}")
        out.append(found)
    return out


def _stage_proteins(files: list[Path]) -> Iterator[tuple[bytes, str]]:
    """(protein_uid, seq) of one stage's `protein_seq` files, in protein_uid order, `SCRATCH_ROWS` at a
    time. The files are sorted by protein_uid and bucketed by its prefix (`protein_seq`'s sort key), so in
    bucket and part order they are one sorted run; a uid that does not follow the one before is a
    RuntimeError, so no unsorted input is ever written as if it were sorted."""
    last = None
    for f in files:
        for b in pq.ParquetFile(f).iter_batches(batch_size=SCRATCH_ROWS, columns=["protein_uid", "seq"]):
            uids = b.column(0).to_pylist()
            if uids and last is not None and uids[0] <= last or any(x >= y for x, y in zip(uids, uids[1:])):
                raise RuntimeError(f"{f}: protein_seq rows are not in strictly increasing protein_uid order")
            if uids:
                last = uids[-1]
            yield from zip(uids, b.column(1).to_pylist())


def _sql_path(p: Path) -> str:
    return "'" + str(p).replace("'", "''") + "'"


class _Flat:
    """Rows appended one at a time and written to a flat Parquet file `SCRATCH_ROWS` at a time."""

    def __init__(self, path: Path, schema: pa.Schema):
        self.schema = schema
        self.cols: list[list] = [[] for _ in schema]
        self.w = pq.ParquetWriter(path, schema, compression="zstd", compression_level=1)

    def add(self, *row) -> None:
        for c, v in zip(self.cols, row):
            c.append(v)
        if len(self.cols[0]) >= SCRATCH_ROWS:
            self.flush()

    def flush(self) -> None:
        if self.cols[0]:
            self.w.write_table(pa.table([pa.array(c, f.type) for c, f in zip(self.cols, self.schema)],
                                        schema=self.schema))
            for c in self.cols:
                c.clear()

    def close(self) -> None:
        self.flush()
        self.w.close()


def _schema(**cols: pa.DataType) -> pa.Schema:
    return pa.schema([pa.field(n, t, nullable=False) for n, t in cols.items()])


def _write_fasta(stage_files: list[list[Path]], faa: Path, uids: Path) -> int:
    """The distinct proteins of the stages' files as FASTA records `>{uid hex}` sorted by protein_uid:
    each stage's files are one sorted run (`_stage_proteins`) and the runs are merged, so no sort holds
    the sequences. A protein in two stages is written once (its uid is the hash of its sequence). The
    uids, as hex in the same order, go to the flat file `uids`. Returns their number."""
    n, last = 0, None
    lines: list[str] = []
    out = _Flat(uids, _schema(uid=pa.string()))
    try:
        with open(faa, "w") as fh:
            for u, s in heapq.merge(*(_stage_proteins(f) for f in stage_files), key=lambda r: r[0]):
                if u == last:
                    continue
                last, h = u, u.hex()
                out.add(h)
                lines.append(f">{h}\n{s}\n")
                n += 1
                if len(lines) == SCRATCH_ROWS:
                    fh.write("".join(lines))
                    lines.clear()
            fh.write("".join(lines))
    finally:
        out.close()
    return n


def _tsv_to_parquet(src: Path, dst: Path, columns: dict[str, pa.DataType], line_numbers: bool = False) -> None:
    """An MMseqs2 TSV (no header, no quoting) as flat Parquet, read a block at a time in file order; with
    `line_numbers`, a column `ln` holds each row's 0-based line number among the non-empty lines."""
    names = list(columns)
    reader = pacsv.open_csv(src, read_options=pacsv.ReadOptions(column_names=names, use_threads=False,
                                                                  block_size=_CSV_BLOCK),
                            parse_options=pacsv.ParseOptions(delimiter="\t", quote_char=False, escape_char=False),
                            convert_options=pacsv.ConvertOptions(column_types=columns, null_values=[],
                                                                 strings_can_be_null=False))
    schema = _schema(**columns, **({"ln": pa.int64()} if line_numbers else {}))
    n = 0
    with pq.ParquetWriter(dst, schema, compression="zstd", compression_level=1) as w:
        for b in reader:
            cols = [b.column(c) for c in names]
            if line_numbers:
                cols.append(pa.array(np.arange(n, n + b.num_rows, dtype=np.int64)))
            n += b.num_rows
            w.write_batch(pa.RecordBatch.from_arrays(cols, schema=schema))


def _sorted(con, sql: str, path: Path) -> Path:
    """The result of `sql` (which ends in its ORDER BY) as a flat Parquet file in that order: one sort at
    a time, so a sort never shares the memory limit with another one being read."""
    con.execute(f"COPY ({sql}) TO {_sql_path(path)} (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE {SCRATCH_ROWS})")
    return path


def _flat_rows(path: Path, key: int = 0) -> Iterator[tuple]:
    """The rows of a flat Parquet file, `SCRATCH_ROWS` at a time; with `key`, a row whose first `key`
    values sort before the row above's is a RuntimeError (a merge must never see unsorted input)."""
    last = None
    for b in pq.ParquetFile(path).iter_batches(batch_size=SCRATCH_ROWS):
        for row in zip(*(b.column(i).to_pylist() for i in range(b.num_columns))):
            if key:
                k = row[:key]
                if last is not None and k < last:
                    raise RuntimeError(f"{path.name}: rows are not sorted on their first {key} columns")
                last = k
            yield row


def _pick_alignments(con, w: Path, level: str, prot: Path | None = None) -> int:
    """Check that every protein is a member of exactly one cluster at `level`, and write to
    `pick_<level>.parquet` each member's first alignment line whose query is the member's representative
    (as the builder before card B-24 kept it: lines in file order, a line of another query passed over).
    Returns the number of members with an alignment. Three sorted streams are merged: the proteins
    (`prot`, default `w/prot.parquet`, uid order), the members (by uid) and the alignment lines (by target,
    then line)."""
    prot = w / "prot.parquet" if prot is None else prot
    clu, aln = _sql_path(w / f"clu_{level}.parquet"), _sql_path(w / f"aln_{level}.parquet")
    members = _flat_rows(_sorted(con, f"SELECT member, rep FROM read_parquet({clu}) ORDER BY member",
                                 w / f"by_member_{level}.parquet"), key=1)
    lines = _flat_rows(_sorted(con, "SELECT target, ln, query, qstart, qend, tstart, tend, fident, alnlen, cigar "
                                    f"FROM read_parquet({aln}) ORDER BY target, ln", w / f"by_target_{level}.parquet"),
                       key=2)
    out = _Flat(w / f"pick_{level}.parquet", _schema(**{c: t for c, t in _ALN_COLUMNS.items()}))
    picked = 0
    try:
        head = next(lines, None)
        uids = (u for (u,) in _flat_rows(prot, key=1))
        for mr, u in itertools.zip_longest(members, uids, fillvalue=None):
            m, r = (None, None) if mr is None else mr
            if m != u:
                lost = m is None or (u is not None and u < m)
                raise RuntimeError(f"{level}: clustering does not cover every protein once: "
                                   + (f"protein {u} is in no cluster" if lost else
                                      f"member {m} is listed twice or is not a protein"))
            while head is not None and head[0] < m:  # a line whose target is no member: passed over
                head = next(lines, None)
            got = False
            while head is not None and head[0] == m:
                if not got and head[2] == r:
                    out.add(head[2], head[0], *head[3:])
                    got = True
                    picked += 1
                head = next(lines, None)
    finally:
        members.close()
        lines.close()
        out.close()
    return picked


def _cluster_members(con, w: Path, level: str, phash: str) -> int:
    """Write `mem_<level>.parquet`, every member with its representative and cluster uid; returns the
    number of clusters. The members arrive sorted by (rep, member), so each cluster's sorted members are
    hashed when its last one has passed (the digest of `cluster_uid`)."""
    rows = _flat_rows(_sorted(con, f"SELECT rep, member FROM read_parquet({_sql_path(w / f'clu_{level}.parquet')}) "
                                   "ORDER BY rep, member", w / f"by_rep_{level}.parquet"), key=2)
    out = _Flat(w / f"mem_{level}.parquet", _schema(member=pa.string(), rep=pa.string(), cluster_uid=pa.binary(16)))
    tail = phash.encode()
    n, cur, ms = 0, None, []

    def close() -> None:
        h = blake3()
        for m in ms:
            h.update(bytes.fromhex(m))
        h.update(tail)
        cu = h.digest(length=16)
        for m in ms:
            out.add(m, cur, cu)

    try:
        for r, m in rows:
            if r != cur:
                if ms:
                    close()
                    n += 1
                cur, ms = r, []
            ms.append(m)
        if ms:
            close()
            n += 1
    finally:
        rows.close()
        out.close()
    return n


def write_member_cluster(member_files: list[Path], dest: Path, *, tmp: Path | None = None,
                         max_file_bytes: int = layout.MAX_FILE_BYTES, memory_limit: str | None = None) -> dict:
    """Write `member_cluster` (card B-27) under `dest` (`l2/member_cluster/release=<R>`): the rows of the
    release's `cluster_member` files `member_files` (member_uid, level, cluster_uid, rep_uid), sorted by
    (level, member_uid) with one DuckDB sort on one thread, capped at `memory_limit` (default
    `DUCK_MEMORY_LIMIT`) and spilling into a private directory under `tmp` (system temp by default), streamed
    to the writer a batch at a time. Refuses (FileExistsError) when `dest` already holds parts, as every
    writer does. Returns rows and parts."""
    spec = table("member_cluster")
    if not member_files:
        raise FileNotFoundError(f"{dest}: no cluster_member files to sort into member_cluster")
    have = io.part_files(dest)
    if have:
        raise FileExistsError(f"{have[0]}: member_cluster is written once; refusing to overwrite")
    listing = "[" + ", ".join(_sql_path(Path(f)) for f in member_files) + "]"
    with tempfile.TemporaryDirectory(dir=tmp) as d:
        con = duckdb.connect(config={"threads": 1, "temp_directory": str(Path(d) / "duckdb"),
                                     "memory_limit": DUCK_MEMORY_LIMIT if memory_limit is None else memory_limit})
        try:
            sql = (f"SELECT {', '.join(spec.schema.names)} FROM read_parquet({listing}, hive_partitioning = false) "
                   "ORDER BY level, member_uid")
            paths = io.write_parts(io.duck_batches(con, sql, spec), spec, dest, EST_MEMBER_CLUSTER_ROW, max_file_bytes)
        finally:
            con.close()
    return {"rows": sum(pq.read_metadata(p).num_rows for p in paths), "parts": len(paths)}


def run(root: Path, release: str, gene_call_version: str, stages: list[str], threads: int = 1,
        tmp: Path | None = None, max_file_bytes: int = layout.MAX_FILE_BYTES, out: Path | None = None,
        memory_limit: str | None = None) -> dict:
    """Cluster the distinct proteins of `stages`, the release's stages (card B-25), and write the
    release's cluster tables; the release's `release_def` is written, or checked to list the same
    stages (`catalog.define_release`).

    The proteins are read from `root`; the cluster tables, the release's clustering ledger and (when
    the root has none) its `release_def` are written under `out` (default `root`; card C-08). Refuses,
    before any work, when those tables or that ledger already hold a file in `root` or `out`; with
    `out` given, `tmp` may lie under neither. `memory_limit` is DuckDB's cap (default
    `DUCK_MEMORY_LIMIT`, read when called)."""
    root = Path(root)
    if out is not None:
        ledger.check_scratch(tmp, root, out)
    out = root if out is None else Path(out)
    stages = catalog.check_stages(stages)
    catalog.verify_release_def(root, release, stages)  # a disagreeing definition is named before a refusal
    ledger.refuse_written(root, out, [f"l2/cluster_member/release={release}", f"l2/cluster_member_aln/release={release}",
                                      f"l2/member_cluster/release={release}",
                                      ledger.release_path("", release, "cluster").parent.as_posix()],
                          f"clustering of release {release}")
    files = _protein_files(root, gene_call_version, stages)
    limit = DUCK_MEMORY_LIMIT if memory_limit is None else memory_limit
    with tempfile.TemporaryDirectory(dir=tmp) as d:
        w = Path(d)
        n_prot = _write_fasta(files, w / "in.faa", w / "prot.parquet")
        if not n_prot:
            raise FileNotFoundError(f"release {release}: no proteins in stages {stages}")
        catalog.define_release(root, release, stages, out=out)
        mm = MMseqs(w, threads)
        mm("createdb", "in.faa", "db", "--dbtype", "1", "--shuffle", "0", "--write-lookup", "1", threaded=False)
        mm("cluster", "db", "clu90", "tmp90", "--min-seq-id", str(LEVELS["nr90"]["min_seq_id"]),
           "-c", str(LEVELS["nr90"]["coverage"]), *COMMON)
        mm("createsubdb", "clu90", "db", "rep90", threaded=False)
        mm("cluster", "rep90", "clu50r", "tmp50", "--min-seq-id", str(LEVELS["nr50"]["min_seq_id"]),
           "-c", str(LEVELS["nr50"]["coverage"]), *COMMON)
        mm("mergeclusters", "db", "clu50", "clu90", "clu50r")
        return _write_tables(mm, w, w, out, release, n_prot, limit, max_file_bytes, member_align=True)


def _write_tables(mm: MMseqs, w: Path, f: Path, out: Path, release: str, n_prot: int, limit: str,
                  max_file_bytes: int, member_align: bool) -> dict:
    """From the clusterings `clu90` and `clu50` in the MMseqs2 working directory `w` (and `w/prot.parquet`): each
    level's members' TSV, the members aligned to their representatives (`aln_<level>`, made here when
    `member_align`, else made before, in parts, card B-41) and that alignment's TSV; then the release's cluster
    tables, its `member_cluster` and its ledger under `out`. The TSVs and the Python side's flat files go to `f`
    (`w` itself, or a directory under it). Closes `mm`'s log. Returns `run`'s report."""
    rel = "" if f == w else f.relative_to(w).as_posix() + "/"
    for level, clu in (("nr90", "clu90"), ("nr50", "clu50")):
        mm("createtsv", "db", "db", clu, f"{rel}{clu}.tsv")
        _tsv_to_parquet(f / f"{clu}.tsv", f / f"clu_{level}.parquet", {"rep": pa.string(), "member": pa.string()})
        if member_align:
            mm("align", "db", "db", clu, f"aln_{level}", "-a", "-e", "inf", "--min-seq-id", "0", "-c", "0",
               "--alignment-mode", "3")
        mm("convertalis", "db", "db", f"aln_{level}", f"{rel}aln_{level}.tsv", "--format-output", _ALN_FORMAT)
        _tsv_to_parquet(f / f"aln_{level}.tsv", f / f"aln_{level}.parquet", _ALN_COLUMNS, line_numbers=True)
    mm.log.close()

    con = duckdb.connect(config={"threads": 1, "memory_limit": limit, "temp_directory": str(f / "duckdb"),
                                 "preserve_insertion_order": True})
    try:
        aligned = {lvl: _pick_alignments(con, f, lvl, prot=w / "prot.parquet") for lvl in LEVELS}
        counts = {lvl: _cluster_members(con, f, lvl, params_hash(lvl)) for lvl in LEVELS}
        mem_spec, aln_spec = table("cluster_member"), table("cluster_member_aln")
        order = sorted(LEVELS)  # `level` leads both sort keys

        def members() -> Iterator[pa.RecordBatch]:
            for lvl in order:
                yield from io.duck_batches(con, (
                    f"SELECT unhex(member) AS member_uid, '{lvl}' AS level, cluster_uid, unhex(rep) AS rep_uid, "
                    f"'{params_hash(lvl)}' AS cluster_params_hash FROM read_parquet({_sql_path(f / f'mem_{lvl}.parquet')}) "
                    "ORDER BY cluster_uid, member_uid"), mem_spec)

        def alignments() -> Iterator[pa.RecordBatch]:
            # convertalis: query = representative, target = member; coordinates 1-based inclusive
            for lvl in order:
                yield from io.duck_batches(con, (
                    f"SELECT unhex(target) AS member_uid, '{lvl}' AS level, unhex(query) AS rep_uid, "
                    "tstart - 1 AS qstart, tend AS qend, qstart - 1 AS tstart, qend AS tend, fident AS pident, "
                    f"alnlen, encode(cigar) AS cigar FROM read_parquet({_sql_path(f / f'pick_{lvl}.parquet')}) "
                    "ORDER BY rep_uid, member_uid"), aln_spec)

        io.write_parts(members(), mem_spec, out / f"l2/cluster_member/release={release}", EST_MEMBER_ROW,
                       max_file_bytes)
        io.write_parts(alignments(), aln_spec, out / f"l2/cluster_member_aln/release={release}", EST_ALN_ROW,
                       max_file_bytes)
    finally:
        con.close()
    # the lookup by member (card B-27), from the cluster_member files just written
    write_member_cluster(io.part_files(out / f"l2/cluster_member/release={release}"),
                         out / f"l2/member_cluster/release={release}", tmp=f, max_file_bytes=max_file_bytes,
                         memory_limit=limit)

    # the release's own ledger (ledger/release=<R>/step=cluster, card B-28): NR90 sees every protein, NR50
    # the NR90 representatives
    n_in = {"nr90": n_prot, "nr50": counts["nr90"]}
    tools = provenance.tool_versions("mmseqs2", "pyarrow")
    ledger.write_release(out, release, "cluster",
                 sorted((ledger.row(lvl.encode(), "level", n_in[lvl], counts[lvl], n_in[lvl] - counts[lvl], tools=tools)
                         for lvl in LEVELS), key=lambda r: (r["object_key"], r["reason_code"])))
    return {"proteins": n_prot, "clusters": counts, "aligned": aligned}


# --- clustering in parts (card B-41; design planning/design/B-37_split_clustering.md) -------------------------
#
# `run` makes two `cluster` workflow calls, each one process. In parts, the module calls those workflow calls
# make (`workflow`, from their `-v 3` log: 15 at NR90, 23 at NR50; a test checks the list against the image's
# MMseqs2) are replayed one process each, with the logged arguments, and so are `run`'s own calls between and
# after them (`calls`). Every per-query module (`prefilter`, `rescorediagonal`, `align`, the member `align`) can
# run as N processes, each over a consecutive range of the sorted keys of its split input (prefilter's query
# database, the others' result database; the other databases whole), joined in key order (`join_dbs`). The
# parts of a prefilter ask for the self hit (`--add-self-matches 1`): a part's query is not the target by path,
# so MMseqs2 does not search it as a database against itself, and the self hit of a short or masked protein is
# lost otherwise (measured in the design). The target is never split, and a prefilter process whose log does
# not say `Process prefiltering step 1 of 1` is refused (`check_prefilter_log`): MMseqs2's own target split
# changes results. With these, every database of the pipeline and the published tables equal `run`'s at any N.
#
# Steps, each its own process (build-stage runs each as a `run-step`; `run_split` runs them in this process):
#
#   cluster.plan                   the FASTA and `createdb` (`plan`), into the work directory
#   <call>                         a call that runs whole: one-process modules, and a per-query call in one part
#   <call>.cut                     a per-query call's parts (`cut`): from its input's keys and the budget
#   <call>.part-<p>, <call>        its parts, then their join under the call's own name
#   cluster                        `finish`: TSVs, the Python side and the tables, as `run` writes them
#
# Each step but the last writes one record, `l2/cluster_plan/release=<R>/<step>.json` (the `cluster_plan`
# artifact): the digests of its inputs (from their producers' records) and of its outputs, each a sha256 over
# the database's type and (key, entry) pairs in key order (`db_digest`), not over files (the `.0`, `.1` layout
# depends on the thread count). The MMseqs2 databases stay in the work directory, `<tmp>/biobox-cluster/<R>-<key>`
# (`work_dir`, `key` from the root and the plan's inputs); the published files are `run`'s only. A part's files
# stay until its join is published (`clean_parts`); the work directory is removed after `cluster` publishes.
#
# The parts are a pure function of the input's sorted keys, their sizes and the plan's cut, never of threads or
# the hour: a forced count (`parts`, for tests and trials), or what fits `part_cpu_min` (80% of the budget in
# build-stage) at the call's estimate (`CALL_CPU_S`). Lost scratch: `lost` names every recorded step whose
# outputs are missing or are not its record's, and `redo` makes each again, which must reproduce its record.

PLAN_FORMAT = 1
WORK_DIR = "biobox-cluster"  # under the scratch directory: one work directory per root, release and plan
PLAN_STEP = "cluster.plan"
FINAL_STEP = "cluster"
PLAN_DIR = ARTIFACTS["cluster_plan"].location
# per-query modules: the positional argument that is split (prefilter's query database; the result database of
# the others), consecutive key ranges of it, one part per process
PER_QUERY = {"prefilter": 0, "rescorediagonal": 2, "align": 2}
UNTHREADED = {"createsubdb"}  # as `run` calls it and as the B-37 replay did
PREFILTER_ONE_STEP = "Process prefiltering step 1 of 1"
FINISH_INPUTS = ("db", "db_h", "prot.parquet", "clu90", "clu50", "aln_nr90", "aln_nr50")
_SUB_MAT = "aa:blosum62.out,nucl:nucleotide.out"
_PCA, _PCB = "substitution:1.100,context:1.400", "substitution:4.100,context:5.800"
_INDEX_BLOCK = 1 << 22  # bytes per block when reading an MMseqs2 index

# CPU-seconds of each per-query call, measured on wilkcomp02, 2026-11-11, MMseqs2 18.8cc5c, 2 threads per process,
# for 100,000 proteins of stage S1 replayed one module call per process (the B-37 design §2, run `s1_100k_s1`,
# `times.jsonl`), and the run's total over all its 48 calls, 489.8. A call's estimate for a release is its share
# of the release's estimate in one process (`biobox.stage.cluster_cpu_min`, card B-36), *derived*: the shares were
# measured at 100,000 proteins only (the one-process modules' share was the same at 25,000 and 50,000).
CALL_CPU_S = {
    "cluster.nr90.01.rescorediagonal": 0.297, "cluster.nr90.06.rescorediagonal": 0.334,
    "cluster.nr90.07.align": 9.495, "cluster.nr90.11.prefilter": 26.062, "cluster.nr90.12.align": 139.51,
    "cluster.nr50.01.rescorediagonal": 0.256, "cluster.nr50.06.rescorediagonal": 0.315,
    "cluster.nr50.07.align": 14.593, "cluster.nr50.11.prefilter": 19.9, "cluster.nr50.12.align": 51.403,
    "cluster.nr50.15.prefilter": 28.588, "cluster.nr50.16.align": 42.644, "cluster.nr50.19.prefilter": 58.334,
    "cluster.nr50.20.align": 51.046, "cluster.nr90.members": 13.856, "cluster.nr50.members": 23.167,
}
MEASURED_TOTAL_CPU_S = 489.8
# Each prefilter part builds the k-mer index of the whole target again: 12.44 CPU-seconds per extra part at 59,000
# to 78,316 targets (*derived* in the B-37 design §3 from the 4-part run: (282.2 - 132.9) / 12, measured, same
# host and day), taken per target at the smaller count (the larger rate, *derived*: 12.44 / 59,000).
PREFILTER_INDEX_CPU_S_PER_TARGET = 12.44 / 59_000


def _kmermatcher(db: str, out: str, min_id: str) -> tuple[str, ...]:
    return ("kmermatcher", db, out, "--sub-mat", _SUB_MAT, "--alph-size", "aa:13,nucl:5", "--min-seq-id", min_id,
            "--kmer-per-seq", "21", "--spaced-kmer-mode", "1", "--kmer-per-seq-scale", "aa:0.000,nucl:0.200",
            "--adjust-kmer-len", "0", "--mask", "0", "--mask-prob", "0.9", "--mask-lower-case", "0",
            "--mask-n-repeat", "0", "--cov-mode", "1", "-k", "0", "-c", "0.8", "--max-seq-len", "65535",
            "--hash-shift", "67", "--split-memory-limit", "0", "--include-only-extendable", "0",
            "--ignore-multi-kmer", "0", "--compressed", "0", "--cluster-weight-threshold", "0.9")


def _rescore(q: str, res: str, out: str, mode: str, filter_hits: str, min_id: str) -> tuple[str, ...]:
    return ("rescorediagonal", q, q, res, out, "--sub-mat", _SUB_MAT, "--rescore-mode", mode, "--wrapped-scoring",
            "0", "--filter-hits", filter_hits, "-e", "0.001", "-c", "0.8", "-a", "0", "--cov-mode", "1",
            "--min-seq-id", min_id, "--min-aln-len", "0", "--seq-id-mode", "0", "--add-self-matches", "0",
            "--sort-results", "0", "--db-load-mode", "0", "--compressed", "0")


def _clust(db: str, res: str, out: str) -> tuple[str, ...]:
    return ("clust", db, res, out, "--cluster-mode", "2", "--max-iterations", "1000", "--similarity-type", "2",
            "--compressed", "0", "--cluster-weight-threshold", "0.9")


def _subdb(ids_: str, db: str, out: str) -> tuple[str, ...]:
    return ("createsubdb", ids_, db, out, "--subdb-mode", "1")


def _align(q: str, res: str, out: str, min_id: str, comp_bias: str) -> tuple[str, ...]:
    return ("align", q, q, res, out, "--sub-mat", _SUB_MAT, "-a", "0", "--alignment-mode", "3",
            "--alignment-output-mode", "0", "--wrapped-scoring", "0", "-e", "0.001", "--min-seq-id", min_id,
            "--min-aln-len", "0", "--seq-id-mode", "0", "--alt-ali", "0", "-c", "0.8", "--cov-mode", "1",
            "--max-seq-len", "65535", "--comp-bias-corr", comp_bias, "--comp-bias-corr-scale", "1",
            "--max-rejected", "2147483647", "--max-accept", "2147483647", "--add-self-matches", "0",
            "--db-load-mode", "0", "--pca", _PCA, "--pcb", _PCB, "--score-bias", "0", "--realign", "0",
            "--realign-score-bias", "-0.2", "--realign-max-seqs", "2147483647", "--corr-score-weight", "0",
            "--gap-open", "aa:11,nucl:5", "--gap-extend", "aa:1,nucl:2", "--zdrop", "40", "--compressed", "0")


def _prefilter(q: str, out: str, s: str, comp_bias: str, diag: str, min_ungapped: str) -> tuple[str, ...]:
    return ("prefilter", q, q, out, "--sub-mat", _SUB_MAT, "--seed-sub-mat", "aa:VTML80.out,nucl:nucleotide.out",
            "-s", s, "-k", "0", "--target-search-mode", "0", "--k-score", "seq:2147483647,prof:2147483647",
            "--alph-size", "aa:21,nucl:5", "--max-seq-len", "65535", "--max-seqs", "20", "--split", "0",
            "--split-mode", "2", "--split-memory-limit", "0", "-c", "0.8", "--cov-mode", "1", "--comp-bias-corr",
            comp_bias, "--comp-bias-corr-scale", "1", "--diag-score", diag, "--exact-kmer-matching", "0", "--mask",
            "1", "--mask-prob", "0.9", "--mask-lower-case", "0", "--mask-n-repeat", "0", "--min-ungapped-score",
            min_ungapped, "--add-self-matches", "0", "--spaced-kmer-mode", "1", "--db-load-mode", "0", "--pca", _PCA,
            "--pcb", _PCB, "--compressed", "0")


def workflow(level: str) -> list[tuple[str, ...]]:
    """The module calls of `run`'s `cluster` workflow call at `level`, as MMseqs2 18.8cc5c makes them (its `-v 3`
    log on the fixture, B-37 design §1), without `--threads` and `-v`; paths are relative to the work directory,
    the workflow's temporary directory being `tmp90/w` or `tmp50/w` and its linclust directory `.../linclust/l`.
    The cascade without reassignment: linclust (redundancy), then one prefilter-align-clust step at NR90 (`-s 1`)
    and three at NR50 (`-s 1`, `2`, `3`)."""
    if level == "nr90":
        db, out, min_id, w, lin_bias, steps = "db", "clu90", "0.9", "tmp90/w", "0", [("1", "0", "0", "0")]
    elif level == "nr50":
        db, out, min_id, w, lin_bias = "rep90", "clu50r", "0.5", "tmp50/w", "1"
        steps = [("1", "0", "0", "0"), ("2", "1", "1", "15"), ("3", "1", "1", "15")]
    else:
        raise ValueError(f"level {level!r}: one of {sorted(LEVELS)}")
    lin = f"{w}/linclust/l"
    red = f"{lin}/input_step_redundancy"
    calls_ = [
        _kmermatcher(db, f"{lin}/pref", min_id),
        _rescore(db, f"{lin}/pref", f"{lin}/pref_rescore1", "0", "0", min_id),
        _clust(db, f"{lin}/pref_rescore1", f"{lin}/pre_clust"),
        _subdb(f"{lin}/order_redundancy", db, red),
        _subdb(f"{lin}/order_redundancy", f"{lin}/pref", f"{lin}/pref_filter1"),
        ("filterdb", f"{lin}/pref_filter1", f"{lin}/pref_filter2", "--filter-file", f"{lin}/order_redundancy",
         "--compressed", "0"),
        _rescore(red, f"{lin}/pref_filter2", f"{lin}/pref_rescore2", "1", "1", min_id),
        _align(red, f"{lin}/pref_rescore2", f"{lin}/aln", min_id, lin_bias),
        _clust(red, f"{lin}/aln", f"{lin}/clust"),
        ("mergeclusters", db, f"{w}/clu_redundancy", f"{lin}/pre_clust", f"{lin}/clust", "--compressed", "0"),
        _subdb(f"{w}/clu_redundancy", db, f"{w}/input_step_redundancy"),
    ]
    inp = f"{w}/input_step_redundancy"
    for k, (s, bias, diag, ungapped) in enumerate(steps):
        if k:
            nxt = f"{w}/input_step{k}"
            calls_.append(_subdb(f"{w}/clu_step{k - 1}", inp, nxt))
            inp = nxt
        calls_ += [_prefilter(inp, f"{w}/pref_step{k}", s, bias, diag, ungapped),
                   _align(inp, f"{w}/pref_step{k}", f"{w}/aln_step{k}", min_id, bias),
                   _clust(inp, f"{w}/aln_step{k}", f"{w}/clu_step{k}")]
    calls_.append(("mergeclusters", db, out, f"{w}/clu_redundancy", *(f"{w}/clu_step{k}" for k in range(len(steps))),
                   "--compressed", "0"))
    return calls_


@dataclass(frozen=True)
class Call:
    """One module call of clustering in parts: its step name and its MMseqs2 arguments (paths relative to the work
    directory; no `--threads` or `-v`)."""
    name: str
    argv: tuple[str, ...]

    @property
    def module(self) -> str:
        return self.argv[0]

    @property
    def split(self) -> bool:
        """A per-query call: it can run in parts."""
        return self.module in PER_QUERY

    @property
    def positional(self) -> list[str]:
        return list(itertools.takewhile(lambda a: not a.startswith("-"), self.argv[1:]))

    @property
    def out_index(self) -> int:
        """The output's index among the positional arguments (mergeclusters: its second; else its last)."""
        return 1 if self.module == "mergeclusters" else len(self.positional) - 1

    @property
    def output(self) -> str:
        return self.positional[self.out_index]

    @property
    def order_file(self) -> str | None:
        """The linclust `clust` call's companion: `linclust.sh` writes the first column of `pre_clust.index` to
        `order_redundancy` (awk) right after it; here the same step writes it."""
        return f"{Path(self.output).parent.as_posix()}/order_redundancy" if self.output.endswith("/pre_clust") else None

    @property
    def inputs(self) -> list[str]:
        p = self.positional
        flagged = [self.argv[i + 1] for i, a in enumerate(self.argv) if a == "--filter-file"]
        return list(dict.fromkeys(p[:self.out_index] + p[self.out_index + 1:] + flagged))


def calls() -> list[Call]:
    """Every module call after the plan, in order (module docstring): the NR90 workflow call's, the NR90
    representatives, the NR50 workflow call's, the merge of the two levels, and each level's member alignment."""
    out = []
    for level in ("nr90", "nr50"):
        out += [Call(f"cluster.{level}.{i:02d}.{a[0]}", a) for i, a in enumerate(workflow(level))]
        if level == "nr90":
            out.append(Call("cluster.rep90", ("createsubdb", "clu90", "db", "rep90")))
    out.append(Call("cluster.merge", ("mergeclusters", "db", "clu50", "clu90", "clu50r")))
    for level, clu in (("nr90", "clu90"), ("nr50", "clu50")):
        out.append(Call(f"cluster.{level}.members", ("align", "db", "db", clu, f"aln_{level}", "-a", "-e", "inf",
                                                     "--min-seq-id", "0", "-c", "0", "--alignment-mode", "3")))
    return out


def part_name(call: str, k: int) -> str:
    return f"{call}.part-{k:05d}"


_PART = re.compile(r"(?P<call>.+)\.part-(?P<k>\d{5})")


def _parse(name: str) -> tuple[Call, str, int | None]:
    """(the call, "call" | "cut" | "part", the part) of a step name; ValueError for a name of no step."""
    by = {c.name: c for c in calls()}
    if name.endswith(".cut") and name[:-4] in by and by[name[:-4]].split:
        return by[name[:-4]], "cut", None
    m = _PART.fullmatch(name)
    if m and m["call"] in by and by[m["call"]].split:
        return by[m["call"]], "part", int(m["k"])
    if name in by:
        return by[name], "call", None
    raise ValueError(f"{name}: not a step of clustering in parts (cluster.calls(), a .cut or .part-<p> of a "
                     f"per-query call, or {PLAN_STEP})")


# --- MMseqs2 databases: read, digest, join -------------------------------------------------------------------


def _data_files(base: Path) -> list[Path]:
    """The data file of database `base`, or its numbered files `.0`, `.1`, ... (one per thread that wrote it)."""
    base = Path(base)
    if base.is_file():
        return [base]
    out = []
    while (p := base.with_name(f"{base.name}.{len(out)}")).is_file():
        out.append(p)
    if not out:
        raise FileNotFoundError(f"{base}: no MMseqs2 data file (neither it nor {base.name}.0 exists)")
    return out


def _index_blocks(index: Path) -> Iterator[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """(keys, offsets, lengths) of an MMseqs2 `.index` file, a block at a time, in file order."""
    if index.stat().st_size == 0:
        return
    cols = {"k": pa.uint64(), "o": pa.uint64(), "n": pa.uint64()}
    reader = pacsv.open_csv(index, read_options=pacsv.ReadOptions(column_names=list(cols), use_threads=False,
                                                                    block_size=_INDEX_BLOCK),
                            parse_options=pacsv.ParseOptions(delimiter="\t", quote_char=False, escape_char=False),
                            convert_options=pacsv.ConvertOptions(column_types=cols))
    for b in reader:
        yield tuple(b.column(i).to_numpy() for i in range(3))


def _index_sorted(index: Path) -> bool:
    last = None
    for k, _, _ in _index_blocks(index):
        if len(k) and ((last is not None and k[0] <= last) or bool(np.any(k[1:] <= k[:-1]))):
            return False
        if len(k):
            last = k[-1]
    return True


def index_rows(base: Path) -> Iterator[tuple[int, int, int]]:
    """(key, offset, length) of every entry of database `base` in key order, streamed: in file order when its
    index is sorted (as MMseqs2 writes it), else through one DuckDB sort that spills to a private directory."""
    index = Path(f"{base}.index")
    if _index_sorted(index):
        for k, o, n in _index_blocks(index):
            yield from zip(k.tolist(), o.tolist(), n.tolist())
        return
    with tempfile.TemporaryDirectory(prefix=".sort-", dir=index.parent) as d:
        con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT, "temp_directory": d})
        try:
            res = con.execute(f"SELECT k, o, n FROM read_csv({_sql_path(index)}, delim = '\t', header = false, "
                              "columns = {'k': 'UBIGINT', 'o': 'UBIGINT', 'n': 'UBIGINT'}) ORDER BY k")
            while rows := res.fetchmany(SCRATCH_ROWS):
                yield from rows
        finally:
            con.close()


class _Data:
    """The data files of a database as one run of bytes (MMseqs2's offsets run across `.0`, `.1`, ...), mapped."""

    def __init__(self, base: Path):
        self.starts, self.maps, self.files, pos = [], [], [], 0
        for f in _data_files(base):
            fh = open(f, "rb")
            size = os.fstat(fh.fileno()).st_size
            self.files.append(fh)
            self.maps.append(mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) if size else None)
            self.starts.append(pos)
            pos += size

    def read(self, off: int, n: int) -> bytes:
        i = bisect.bisect_right(self.starts, off) - 1
        m = self.maps[i] if i >= 0 else None
        lo = off - self.starts[i] if i >= 0 else 0
        if m is None or lo + n > len(m):
            return b""
        return m[lo:lo + n]

    def close(self) -> None:
        for m in self.maps:
            if m is not None:
                m.close()
        for fh in self.files:
            fh.close()


def db_entries(base: Path) -> Iterator[tuple[int, bytes]]:
    """(key, entry) of every entry of database `base`, in key order, each entry with its closing NUL; an entry
    that does not end in a NUL where its index says is a ValueError."""
    data = _Data(base)
    try:
        for k, o, n in index_rows(base):
            e = data.read(o, n)
            if len(e) != n or not e.endswith(b"\0"):
                raise ValueError(f"{base}: entry {k} at {o}+{n} is not a NUL-terminated entry of its data")
            yield k, e
    finally:
        data.close()


def db_digest(base: Path) -> str:
    """sha256 over the database's type (its `.dbtype` file) and its (key, entry) pairs in key order: what its
    readers see, whatever data files and index order it was written with."""
    h = hashlib.sha256(b"dbtype\0" + Path(f"{base}.dbtype").read_bytes())
    for k, e in db_entries(base):
        h.update(k.to_bytes(8, "little"))
        h.update(len(e).to_bytes(8, "little"))
        h.update(e)
    return h.hexdigest()


def _digest(path: Path) -> str | None:
    """A database's `db_digest`, a plain file's sha256, or None when it is not there."""
    path = Path(path)
    try:
        if Path(f"{path}.dbtype").is_file():
            return db_digest(path)
        return ids.sha256_file(path).hex() if path.is_file() else None
    except (FileNotFoundError, ValueError):
        return None


_DB_FILE = r"(_h)?(\.(\d+|index|dbtype|lookup|source))?"


def _clear(w: Path, rel: str) -> None:
    """Remove database `rel` of the work directory (data, numbered data, index, type, lookup, source, headers),
    as links where they are links, before a step writes it again; nothing else."""
    base = Path(w) / rel
    if not base.parent.is_dir():
        base.parent.mkdir(parents=True)
        return
    pat = re.compile(re.escape(base.name) + _DB_FILE)
    for p in base.parent.iterdir():
        if pat.fullmatch(p.name) and (p.is_symlink() or p.is_file()):
            p.unlink()


def join_dbs(parts: Sequence[Path], out: Path) -> None:
    """Write database `out` from `parts`, each a consecutive key range after the one before: their entries in key
    order, one data file. ValueError when a part's key does not follow the one before or their types differ."""
    out = Path(out)
    types = {Path(f"{p}.dbtype").read_bytes() for p in parts}
    if len(types) != 1:
        raise ValueError(f"{out}: the parts are of {len(types)} database types")
    data, index = out.with_name(f".{out.name}.join-{os.getpid()}"), out.with_name(f".{out.name}.join-{os.getpid()}.index")
    try:
        last, off = None, 0
        with open(data, "wb") as d, open(index, "w") as ix:
            for p in parts:
                for k, e in db_entries(p):
                    if last is not None and k <= last:
                        raise ValueError(f"{p}: key {k} does not follow key {last} of the part before; the parts are "
                                         "not consecutive key ranges in order")
                    d.write(e)
                    ix.write(f"{k}\t{off}\t{len(e)}\n")
                    off, last = off + len(e), k
        _clear(out.parent, out.name)
        os.replace(data, out)
        os.replace(index, Path(f"{out}.index"))
        Path(f"{out}.dbtype").write_bytes(types.pop())
    finally:
        data.unlink(missing_ok=True)
        index.unlink(missing_ok=True)


def check_prefilter_log(text: str, where: str) -> None:
    """ValueError unless the prefilter log `text` says `Process prefiltering step 1 of 1` and no other step: with
    `--split 0` MMseqs2 picks its own split from the memory it sees, and a target split changes results (B-37
    design §3: 3,081 of 78,316 entries)."""
    steps = re.findall(r"Process prefiltering step (\d+) of (\d+)", text)
    if steps != [("1", "1")]:
        raise ValueError(f"{where}: its prefilter log must say `{PREFILTER_ONE_STEP}` once, but says "
                         f"{['step %s of %s' % s for s in steps] or 'no step'}: MMseqs2 split the search by its own "
                         "memory estimate, which changes results when it splits the target; refused")


# --- the cut ------------------------------------------------------------------------------------------------


def cut(rows: Callable[[], Iterable[tuple[int, int]]], est_cpu_min: float, *, part_cpu_min: float | None = None,
        overhead_cpu_min: float = 0.0, parts: int | None = None) -> list[dict]:
    """A per-query call's parts: consecutive ranges of `rows()` (key, units: residues or entry bytes; strictly
    increasing keys, read twice), each part the keys whose units before them fall in its share of the total, so
    the parts are about equal in work (a key larger than a share moves the next cut on, but every part gets at
    least one key). Their number is `parts` when given (forced), else the least that keeps
    each part's estimate, its share of `est_cpu_min` plus `overhead_cpu_min` (what every part pays whatever its
    size), within `part_cpu_min`; never more than the keys. A pure function of the rows and these numbers. Each
    part: first and last key, keys, units, and its estimate `cpu_min`. No rows: one empty part."""
    n = total = 0
    last = None
    for k, u in rows():
        if last is not None and k <= last:
            raise ValueError(f"key {k} after {last}: the keys are not in strictly increasing order")
        last, n, total = k, n + 1, total + u
    if parts is not None:
        if isinstance(parts, bool) or not isinstance(parts, int) or parts < 1:
            raise ValueError(f"parts {parts!r}: must be a positive integer")
        want = parts
    else:
        if part_cpu_min is None or part_cpu_min <= 0:
            raise ValueError(f"part_cpu_min {part_cpu_min!r}: must be positive (or give parts)")
        room = part_cpu_min - overhead_cpu_min
        if room <= 0:
            raise ValueError(f"the per-part overhead of {overhead_cpu_min:.3f} CPU-minutes alone fills the part budget "
                             f"of {part_cpu_min:.3f}: no number of parts fits it")
        want = math.ceil(est_cpu_min / room)
    if not n:
        return [{"first": None, "last": None, "keys": 0, "units": 0, "cpu_min": round(est_cpu_min + overhead_cpu_min, 6)}]
    m = max(1, min(want, n))
    out: list[dict] = []
    cum = prev = 0
    for i, (k, u) in enumerate(rows()):
        p = min(m - 1, cum * m // total) if total else min(m - 1, i * m // n)
        # every part gets a key: no part is passed over, and the keys left can fill the parts left
        p = 0 if i == 0 else min(prev + 1, max(p, prev, m - (n - i)))
        prev = p
        if not out or out[-1]["p"] != p:
            out.append({"p": p, "first": k, "last": k, "keys": 0, "units": 0})
        cur = out[-1]
        cur["last"], cur["keys"], cur["units"] = k, cur["keys"] + 1, cur["units"] + u
        cum += u
    for d in out:
        del d["p"]
        share = d["units"] / total if total else d["keys"] / n
        d["cpu_min"] = round(est_cpu_min * share + overhead_cpu_min, 6)
    return out


# --- records and the work directory ----------------------------------------------------------------------------


def plan_dir(release: str) -> str:
    """The root-relative directory of the release's `cluster_plan` records."""
    return f"{PLAN_DIR}/release={release}"


def _rel(release: str, name: str) -> str:
    return f"{plan_dir(release)}/{name}.json"


def _json_bytes(doc: dict) -> bytes:
    return json.dumps(doc, sort_keys=True, indent=1).encode() + b"\n"


def _find(root: Path, out: Path | None, rel: str) -> Path | None:
    """`rel` in the root, else in `out` (written by a direct call into another output root)."""
    for base in [Path(root)] + ([] if out is None else [Path(out)]):
        if (base / rel).is_file():
            return base / rel
    return None


def _record_raw(root: Path, release: str, name: str, out: Path | None = None) -> bytes | None:
    p = _find(root, out, _rel(release, name))
    return None if p is None else p.read_bytes()


def read_record(root: Path, release: str, name: str, out: Path | None = None) -> dict | None:
    """Step `name`'s published record (from the root, else from `out`), or None."""
    raw = _record_raw(root, release, name, out)
    return None if raw is None else json.loads(raw)


def _write_record(root: Path, out: Path, release: str, name: str, doc: dict) -> None:
    rel = _rel(release, name)
    for base in dict.fromkeys([Path(root).resolve(), Path(out).resolve()]):
        if (base / rel).exists():
            raise FileExistsError(f"{base / rel}: the record of {name} is written once; refusing to overwrite")
    path = Path(out) / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    tmp.write_bytes(_json_bytes(doc))
    os.replace(tmp, path)


def _table_dirs(release: str) -> list[str]:
    return [f"l2/cluster_member/release={release}", f"l2/cluster_member_aln/release={release}",
            f"l2/member_cluster/release={release}", ledger.release_path("", release, "cluster").parent.as_posix()]


_INPUT_KEYS = ("format", "release", "stages", "gene_call_version", "mmseqs2", "files", "cut")


def _work_dir(tmp: Path, root: Path, inputs: dict) -> Path:
    key = hashlib.sha256(f"{Path(root).resolve()}\n".encode() + _json_bytes(inputs)).hexdigest()[:16]
    return Path(tmp) / WORK_DIR / f"{inputs['release']}-{key}"


def work_dir(tmp: Path, root: Path, release: str, out: Path | None = None) -> Path:
    """The release's MMseqs2 work directory under the scratch directory `tmp`: named by the root and the plan's
    inputs, so the work of another plan, release or root never mixes in."""
    plan_ = read_record(root, release, PLAN_STEP, out)
    if plan_ is None:
        raise FileNotFoundError(f"release {release}: no plan of clustering in parts ({_rel(release, PLAN_STEP)}); "
                                "write it with bio-box release cluster plan")
    return _work_dir(tmp, root, {k: plan_[k] for k in _INPUT_KEYS})


def _mm(w: Path, log: Path, argv: Sequence[str], threads: int | None) -> None:
    """One MMseqs2 module call in the work directory at `-v 3`, its output appended to `log`."""
    exe = shutil.which("mmseqs")
    if exe is None:
        raise FileNotFoundError("mmseqs not on PATH")
    cmd = [exe, *argv, "-v", "3", *(["--threads", str(threads)] if threads is not None else [])]
    with open(log, "a") as fh:
        fh.write("### " + shlex.join(cmd[1:]) + "\n")
        fh.flush()
        subprocess.run(cmd, check=True, cwd=w, stdout=fh, stderr=subprocess.STDOUT)


def _make_db(root: Path, w: Path, plan_: dict) -> tuple[int, dict[str, str]]:
    """The plan's work: the proteins of its files as FASTA sorted by protein_uid (`_write_fasta`, as `run`) and
    `createdb` (as `run`), into the work directory; the FASTA is removed after. Returns the proteins and the
    outputs' digests."""
    groups = []
    for group in plan_["files"]:
        paths = []
        for f in group:
            p = Path(root) / f["path"]
            if p.stat().st_size != f["bytes"]:
                raise ValueError(f"{p}: {p.stat().st_size} bytes, the plan of release {plan_['release']} lists "
                                 f"{f['bytes']}")
            paths.append(p)
        groups.append(paths)
    for rel in ("db", "db_h"):
        _clear(w, rel)
    (w / "logs").mkdir(parents=True, exist_ok=True)
    n = _write_fasta(groups, w / "in.faa", w / "prot.parquet")
    log = w / "logs" / f"{PLAN_STEP}.log"
    log.write_text("")
    _mm(w, log, ("createdb", "in.faa", "db", "--dbtype", "1", "--shuffle", "0", "--write-lookup", "1"), None)
    (w / "in.faa").unlink()
    return n, {p: _digest(w / p) for p in ("db", "db_h", "prot.parquet")}


def plan(root: Path, release: str, gene_call_version: str, stages: list[str], tmp: Path, *,
         parts: int | None = None, part_cpu_min: float | None = None, out: Path | None = None) -> dict:
    """Step `cluster.plan`: the release's proteins (the `protein_seq` files of `stages`, card B-25) into the MMseqs2
    database of a fresh work directory, and the plan record under `out` (default `root`): the files, the MMseqs2
    version, the cut (`parts`, a forced count, or `part_cpu_min`, the CPU-minutes a part may be estimated at), the
    protein count, the release's estimate in one process and the outputs' digests. Refuses when the release's
    cluster tables or the plan exist in `root` or `out`."""
    root = Path(root)
    out_root = root if out is None else Path(out)
    if tmp is None:
        raise ValueError("clustering in parts needs a scratch directory (--tmp)")
    ledger.check_scratch(tmp, root, out_root)
    if (parts is None) == (part_cpu_min is None):
        raise ValueError("give exactly one of parts (a forced count) and part_cpu_min (the budget of a part)")
    if parts is not None and (isinstance(parts, bool) or not isinstance(parts, int) or parts < 1):
        raise ValueError(f"parts {parts!r}: must be a positive integer")
    if part_cpu_min is not None and not part_cpu_min > 0:
        raise ValueError(f"part_cpu_min {part_cpu_min!r}: must be positive")
    stages = catalog.check_stages(stages)
    catalog.verify_release_def(root, release, stages)
    ledger.refuse_written(root, out_root, _table_dirs(release), f"clustering of release {release}")
    if _find(root, out, _rel(release, PLAN_STEP)) is not None:
        raise FileExistsError(f"{_rel(release, PLAN_STEP)}: the plan of release {release} is written once")
    files = _protein_files(root, gene_call_version, stages)
    inputs = {"format": PLAN_FORMAT, "release": release, "stages": stages, "gene_call_version": gene_call_version,
              "mmseqs2": provenance.tool_versions("mmseqs2")["mmseqs2"],
              "files": [[{"path": p.relative_to(root).as_posix(), "bytes": p.stat().st_size} for p in g] for g in files],
              "cut": {"parts": parts} if parts is not None else {"part_cpu_min": part_cpu_min}}
    w = _work_dir(tmp, root, inputs)
    if w.exists():
        shutil.rmtree(w)  # an earlier attempt of this plan that left no record
    w.mkdir(parents=True)
    n, outputs = _make_db(root, w, inputs)
    if not n:
        raise FileNotFoundError(f"release {release}: no proteins in stages {stages}")
    from biobox.stage import cluster_cpu_min  # the release's estimate in one process (card B-36)

    doc = {**inputs, "proteins": n, "release_cpu_min": round(cluster_cpu_min(n), 6), "outputs": outputs}
    _write_record(root, out_root, release, PLAN_STEP, doc)
    return {"proteins": n, "release_cpu_min": doc["release_cpu_min"], "cut": doc["cut"]}


@dataclass(frozen=True)
class _Ctx:
    root: Path
    out: Path | None
    release: str
    plan: dict
    plan_sha256: str
    w: Path
    threads: int


def _ctx(root: Path, release: str, tmp: Path, out: Path | None, threads: int) -> _Ctx:
    raw = _record_raw(root, release, PLAN_STEP, out)
    if raw is None:
        raise FileNotFoundError(f"release {release}: no plan of clustering in parts; run {PLAN_STEP} first")
    p = json.loads(raw)
    have = provenance.tool_versions("mmseqs2")["mmseqs2"]
    if have != p["mmseqs2"]:
        raise ValueError(f"release {release}: MMseqs2 is {have}, the plan was made with {p['mmseqs2']}; a build "
                         "is made with one version")
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError(f"threads {threads!r}: give the thread count explicitly (at least 1)")
    ledger.check_scratch(tmp, root, out)
    w = _work_dir(tmp, root, {k: p[k] for k in _INPUT_KEYS})
    return _Ctx(Path(root), None if out is None else Path(out), release, p, hashlib.sha256(raw).hexdigest(), w, threads)


def _producers() -> dict[str, str]:
    """Which step writes each path of the work directory."""
    out = {p: PLAN_STEP for p in ("db", "db_h", "prot.parquet")}
    for c in calls():
        out[c.output] = c.name
        if c.order_file:
            out[c.order_file] = c.name
    return out


def _input_digests(ctx: _Ctx, paths: Iterable[str]) -> dict[str, str]:
    """The digests of `paths` as their producers' records hold them (no input is read again)."""
    prods, got = _producers(), {}
    for p in paths:
        rec = read_record(ctx.root, ctx.release, prods[p], ctx.out)
        if rec is None or p not in rec.get("outputs", {}):
            raise FileNotFoundError(f"release {ctx.release}: {p} has no record yet; run {prods[p]} first")
        got[p] = rec["outputs"][p]
    return got


def _outputs(w: Path, paths: Iterable[str]) -> dict[str, str]:
    got = {}
    for p in paths:
        d = _digest(w / p)
        if d is None:
            raise FileNotFoundError(f"{w / p}: the step wrote no such database or file")
        got[p] = d
    return got


def _cut_doc(ctx: _Ctx, c: Call) -> dict:
    src = c.positional[PER_QUERY[c.module]]
    seq = c.module == "prefilter"
    est = ctx.plan["release_cpu_min"] * CALL_CPU_S[c.name] / MEASURED_TOTAL_CPU_S
    overhead = 0.0
    if seq:  # every part builds the k-mer index of the whole target
        with open(ctx.w / f"{c.positional[1]}.index", "rb") as fh:
            targets = sum(1 for line in fh if line.strip())
        overhead = PREFILTER_INDEX_CPU_S_PER_TARGET * targets / 60

    def rows() -> Iterator[tuple[int, int]]:
        for k, _, n in index_rows(ctx.w / src):
            yield k, (max(0, n - 2) if seq else n)  # a sequence entry is its residues, a newline and a NUL

    got = cut(rows, est, overhead_cpu_min=overhead, **ctx.plan["cut"])
    return {"split": src, "units": "residues" if seq else "entry_bytes", "cpu_min": round(est, 6),
            "overhead_cpu_min": round(overhead, 6), "parts": got}


def _run_call(ctx: _Ctx, c: Call, log: Path, name: str) -> tuple[list[str], dict[str, str]]:
    """Call `c` whole, as the workflow made it; its outputs' digests."""
    _clear(ctx.w, c.output)
    _mm(ctx.w, log, c.argv, None if c.module in UNTHREADED else ctx.threads)
    if c.module == "prefilter":
        check_prefilter_log(log.read_text(), name)
    paths = [c.output] + ([f"{c.output}_h"] if Path(ctx.w / f"{c.output}_h.dbtype").is_file() else [])
    if c.order_file:  # linclust.sh: awk '{ print $1 }' pre_clust.index > order_redundancy
        order = ctx.w / c.order_file
        tmp = order.with_name(f".{order.name}.tmp-{os.getpid()}")
        with open(ctx.w / f"{c.output}.index") as src, open(tmp, "w") as dst:
            for line in src:
                dst.write(line.split(None, 1)[0] + "\n" if line.strip() else "\n")
        os.replace(tmp, order)
        paths.append(c.order_file)
    return list(c.argv), _outputs(ctx.w, paths)


def _parts_rel(c: Call) -> str:
    return f"parts/{c.name}"


def _run_part(ctx: _Ctx, c: Call, k: int, cut_doc: dict, log: Path, name: str,
              self_matches: bool) -> tuple[list[str], dict[str, str]]:
    """Part `k` of call `c`: its keys of the split input as a sub-database (linked for a sequence database, copied
    for a result database, as the B-37 replay made them), the module over it with the other databases whole (a
    prefilter part asking for the self hit), into `parts/<call>/out-<k>`."""
    pos = PER_QUERY[c.module]
    src, p = c.positional[pos], cut_doc["parts"][k]
    rel = _parts_rel(c)
    (ctx.w / rel).mkdir(parents=True, exist_ok=True)
    ids_path = ctx.w / rel / f"ids-{k:05d}"
    n = 0
    with open(ids_path, "w") as fh:
        for key, _, _ in index_rows(ctx.w / src):
            if p["first"] <= key <= p["last"]:
                fh.write(f"{key}\n")
                n += 1
    if n != p["keys"]:
        raise ValueError(f"{name}: {n} keys of {src} in {p['first']}..{p['last']}, its cut lists {p['keys']}")
    pin, pout = f"{rel}/in-{k:05d}", f"{rel}/out-{k:05d}"
    _clear(ctx.w, pin)
    _clear(ctx.w, pout)
    _mm(ctx.w, log, ("createsubdb", f"{rel}/ids-{k:05d}", src, pin, "--subdb-mode",
                     "1" if c.module == "prefilter" else "0"), None)
    argv = list(c.argv)
    argv[1 + pos], argv[1 + c.out_index] = pin, pout
    if c.module == "prefilter" and self_matches:
        i = argv.index("--add-self-matches")
        argv[i + 1] = "1"
    _mm(ctx.w, log, argv, ctx.threads)
    if c.module == "prefilter":
        check_prefilter_log(log.read_text(), name)
    return argv, _outputs(ctx.w, [pout])


def _execute(ctx: _Ctx, name: str, self_matches: bool = True) -> dict:
    """Do step `name` in the work directory; return its record."""
    c, kind, k = _parse(name)
    (ctx.w / "logs").mkdir(parents=True, exist_ok=True)
    log = ctx.w / "logs" / f"{name}.log"
    log.write_text("")
    doc = {"format": PLAN_FORMAT, "release": ctx.release, "step": name, "plan_sha256": ctx.plan_sha256,
           "module": c.module}
    if kind == "cut":
        src = c.positional[PER_QUERY[c.module]]
        return {**doc, "inputs": _input_digests(ctx, dict.fromkeys([src, *c.positional[1:2]] if c.module == "prefilter"
                                                                   else [src])), **_cut_doc(ctx, c)}
    cut_doc = read_record(ctx.root, ctx.release, f"{c.name}.cut", ctx.out) if c.split else None
    if c.split and cut_doc is None:
        raise FileNotFoundError(f"{name}: {c.name}.cut has no record; run it first")
    if kind == "part":
        if not 0 <= k < len(cut_doc["parts"]) or len(cut_doc["parts"]) < 2:
            raise ValueError(f"{name}: {c.name} is cut into {len(cut_doc['parts'])} part(s)")
        argv, outputs = _run_part(ctx, c, k, cut_doc, log, name, self_matches)
        p = cut_doc["parts"][k]
        return {**doc, "part": k, "first": p["first"], "last": p["last"], "keys": p["keys"], "argv": argv,
                "inputs": _input_digests(ctx, c.inputs), "outputs": outputs}
    if c.split and len(cut_doc["parts"]) > 1:  # the join of the parts, under the call's own name
        rel, n = _parts_rel(c), len(cut_doc["parts"])
        inputs = {}
        for j in range(n):
            rec = read_record(ctx.root, ctx.release, part_name(c.name, j), ctx.out)
            pout = f"{rel}/out-{j:05d}"
            if rec is None:
                raise FileNotFoundError(f"{name}: {part_name(c.name, j)} has no record; run it first")
            have = _digest(ctx.w / pout)
            if have != rec["outputs"][pout]:
                raise ValueError(f"{name}: {ctx.w / pout} {'is missing' if have is None else 'differs from its record'}; "
                                 f"make it again with {part_name(c.name, j)}")
            inputs[pout] = have
        join_dbs([ctx.w / f"{rel}/out-{j:05d}" for j in range(n)], ctx.w / c.output)
        return {**doc, "parts": n, "inputs": inputs, "outputs": _outputs(ctx.w, [c.output])}
    argv, outputs = _run_call(ctx, c, log, name)
    return {**doc, "argv": argv, "inputs": _input_digests(ctx, c.inputs), "outputs": outputs}


def step(root: Path, release: str, name: str, tmp: Path, threads: int, *, out: Path | None = None,
         _self_matches: bool = True) -> dict:
    """Step `name` of clustering in parts (module docstring) with `threads` MMseqs2 threads: its work in the work
    directory under `tmp`, its record under `out` (default `root`). Refuses a step whose record exists, a release
    whose cluster tables are written, and a step whose inputs have no record. `_self_matches` False (tests only)
    runs prefilter parts as the naive split did, without asking for the self hit."""
    root = Path(root)
    out_root = root if out is None else Path(out)
    ctx = _ctx(root, release, tmp, out, threads)
    ledger.refuse_written(root, out_root, _table_dirs(release), f"clustering of release {release}")
    if _find(root, out, _rel(release, name)) is not None:
        raise FileExistsError(f"{_rel(release, name)}: step {name} of release {release} is done once")
    doc = _execute(ctx, name, _self_matches)
    _write_record(root, out_root, release, name, doc)
    return {"step": name, "outputs": sorted(doc.get("outputs", {})),
            **({"parts": len(doc["parts"])} if name.endswith(".cut") else {})}


def sequence(root: Path, release: str, out: Path | None = None) -> list[str]:
    """The step names in order, as far as they are known now: the plan, then each call (a per-query call's cut,
    then, once the cut is recorded, its parts and its join when it has more than one part); the list stops after a
    cut without a record. `cluster` (`finish`) follows the last call."""
    names = [PLAN_STEP]
    for c in calls():
        if c.split:
            names.append(f"{c.name}.cut")
            rec = read_record(root, release, f"{c.name}.cut", out)
            if rec is None:
                return names
            if len(rec["parts"]) > 1:
                names += [part_name(c.name, k) for k in range(len(rec["parts"]))]
        names.append(c.name)
    return names


def clean_parts(root: Path, release: str, name: str, tmp: Path, out: Path | None = None) -> None:
    """Remove the parts of call `name` from the work directory, once its join is published."""
    if read_record(root, release, name, out) is None:
        raise FileNotFoundError(f"{name}: its join is not published; its parts are still its input")
    shutil.rmtree(work_dir(tmp, root, release, out) / f"parts/{name}", ignore_errors=True)


def lost(root: Path, release: str, tmp: Path, out: Path | None = None) -> list[str]:
    """Every step with a record whose outputs are not in the work directory as its record holds them (missing, or
    another digest), in step order: what `redo` must make again before the build goes on (a reboot clears /tmp).
    A cut has no outputs; a part whose join is published and in place is not needed again."""
    if read_record(root, release, PLAN_STEP, out) is None:
        return []
    w, good, bad = work_dir(tmp, root, release, out), {}, []

    def ok(name: str) -> bool | None:
        if name not in good:
            rec = read_record(root, release, name, out)
            good[name] = None if rec is None else all(_digest(w / p) == d for p, d in rec["outputs"].items())
        return good[name]

    for name in sequence(root, release, out):
        if read_record(root, release, name, out) is None:
            break
        if name.endswith(".cut"):
            continue
        m = _PART.fullmatch(name)
        if m and ok(m["call"]):
            continue
        if not ok(name):
            bad.append(name)
    return bad


def redo(root: Path, release: str, name: str, tmp: Path, threads: int, *, out: Path | None = None) -> dict:
    """Make step `name` again in the work directory for a step whose record is published but whose outputs were
    lost or damaged (`lost`). Its record is not written again: the step done again must reproduce it, else a
    ValueError names the step (the record is then not the work of these inputs: tell the PI)."""
    root = Path(root)
    rec = read_record(root, release, name, out)
    if rec is None:
        raise FileNotFoundError(f"{name} of release {release}: no published record to make it again for; run it as "
                                "a step")
    if name == PLAN_STEP:
        w = work_dir(tmp, root, release, out)
        w.mkdir(parents=True, exist_ok=True)
        n, outputs = _make_db(root, w, rec)
        doc = {**rec, "proteins": n, "outputs": outputs}
    else:
        doc = _execute(_ctx(root, release, tmp, out, threads), name)
    if doc != rec:
        diff = sorted(k for k in set(doc) | set(rec) if doc.get(k) != rec.get(k))
        raise ValueError(f"{name} of release {release}: made again, its record would differ in {diff}; the published "
                         "record is not this step on these inputs, so the build cannot go on (tell the PI)")
    return {"step": name, "outputs": sorted(doc.get("outputs", {})), "verified": True}


def finish(root: Path, release: str, tmp: Path, threads: int, *, out: Path | None = None,
           max_file_bytes: int = layout.MAX_FILE_BYTES, memory_limit: str | None = None) -> dict:
    """Step `cluster` of clustering in parts: from the work directory's clusterings and member alignments, the
    TSVs, the Python side and the tables, `member_cluster`, the ledger and (when the root has none) the release's
    `release_def`, as `run` writes them, under `out` (default `root`). Needs every step's record and its inputs in
    the work directory with their recorded digests."""
    root = Path(root)
    out_root = root if out is None else Path(out)
    ctx = _ctx(root, release, tmp, out, threads)
    stages = catalog.verify_release_def(root, release, ctx.plan["stages"])
    ledger.refuse_written(root, out_root, _table_dirs(release), f"clustering of release {release}")
    names = sequence(root, release, out)
    todo = [n for n in names if read_record(root, release, n, out) is None]
    if todo or names[-1] != calls()[-1].name:
        raise FileNotFoundError(f"release {release}: step {todo[0] if todo else names[-1]} and those after it have no "
                                "record; the clustering is not done")
    want = _input_digests(ctx, FINISH_INPUTS)
    for p, d in want.items():
        if _digest(ctx.w / p) != d:
            raise ValueError(f"{ctx.w / p}: missing or not its record's ({_producers()[p]}); make it again with bio-box "
                             "release cluster redo (build-stage does this itself)")
    catalog.define_release(root, release, stages, out=out_root)
    f = ctx.w / "finish"
    if f.exists():
        shutil.rmtree(f)  # an earlier attempt's
    f.mkdir()
    limit = DUCK_MEMORY_LIMIT if memory_limit is None else memory_limit
    try:
        return _write_tables(MMseqs(ctx.w, threads), ctx.w, f, out_root, release, ctx.plan["proteins"], limit,
                             max_file_bytes, member_align=False)
    finally:
        shutil.rmtree(f, ignore_errors=True)


def run_split(root: Path, release: str, gene_call_version: str, stages: list[str], *, tmp: Path,
              parts: int | None = None, part_cpu_min: float | None = None, threads: int = 1, out: Path | None = None,
              until: str | None = None, max_file_bytes: int = layout.MAX_FILE_BYTES, memory_limit: str | None = None,
              _self_matches: bool = True) -> dict | None:
    """Clustering in parts in this process: the plan, every step without a record in order (each join's parts
    removed once it is recorded), `finish`, then the work directory removed. A step with a record is passed over,
    so a call resumes; `until` stops before that step (`FINAL_STEP`: before `finish`). The files are `run`'s."""
    root, tmp = Path(root), Path(tmp)
    tmp.mkdir(parents=True, exist_ok=True)
    if read_record(root, release, PLAN_STEP, out) is None:
        if until == PLAN_STEP:
            return None
        plan(root, release, gene_call_version, stages, tmp, parts=parts, part_cpu_min=part_cpu_min, out=out)
    while True:
        todo = [n for n in sequence(root, release, out) if read_record(root, release, n, out) is None]
        if not todo:
            break
        if todo[0] == until:
            return None
        step(root, release, todo[0], tmp, threads, out=out, _self_matches=_self_matches)
        c, kind, _ = _parse(todo[0])
        if kind == "call" and c.split and len(read_record(root, release, f"{c.name}.cut", out)["parts"]) > 1:
            clean_parts(root, release, c.name, tmp, out)
    if until == FINAL_STEP:
        return None
    res = finish(root, release, tmp, threads, out=out, max_file_bytes=max_file_bytes, memory_limit=memory_limit)
    shutil.rmtree(work_dir(tmp, root, release, out), ignore_errors=True)
    return res
