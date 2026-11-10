"""Protein clustering with MMseqs2 (operation O-CLUSTER).

NR90 clusters all distinct proteins of one gene-calling version; NR50 clusters
the NR90 representatives, and `mergeclusters` maps every protein to its NR50
representative. Each member is then re-aligned to its representative at each
level (`mmseqs align`, thresholds off) so `cluster_member_aln` has one row per
protein per level, representatives included as self rows.

Writes under `<root>/`:

    l2/cluster_member/release=<R>/part-<k>.parquet
    l2/cluster_member_aln/release=<R>/part-<k>.parquet
    ledger/release=<R>/step=cluster/part-0.parquet

Both tables go through `biobox.io.write_parts`, cut at `max_file_bytes` of estimated size (card
B-19).

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

import hashlib
import heapq
import itertools
import json
import shutil
import subprocess
import tempfile
from collections.abc import Iterator
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq
from blake3 import blake3

from biobox import catalog, io, layout, ledger, provenance
from biobox.registry.schemas import table

# Estimated zstd bytes per row; they only cut parts. Measured on the golden set 36.5 (`cluster_member`)
# and 26.2 (`cluster_member_aln`) bytes per protein and level, 31.0 and 31.0 on the S1 pilot
# (planning/measurements/2026-11-09_design_stress_test.md §3, §6a).
EST_MEMBER_ROW = 37
EST_ALN_ROW = 31

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


def _pick_alignments(con, w: Path, level: str) -> int:
    """Check that every protein is a member of exactly one cluster at `level`, and write to
    `pick_<level>.parquet` each member's first alignment line whose query is the member's representative
    (as the builder before card B-24 kept it: lines in file order, a line of another query passed over).
    Returns the number of members with an alignment. Three sorted streams are merged: the proteins
    (`prot.parquet`, uid order), the members (by uid) and the alignment lines (by target, then line)."""
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
        prot = (u for (u,) in _flat_rows(w / "prot.parquet", key=1))
        for mr, u in itertools.zip_longest(members, prot, fillvalue=None):
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
        for level, clu in (("nr90", "clu90"), ("nr50", "clu50")):
            mm("createtsv", "db", "db", clu, f"{clu}.tsv")
            _tsv_to_parquet(w / f"{clu}.tsv", w / f"clu_{level}.parquet", {"rep": pa.string(), "member": pa.string()})
            mm("align", "db", "db", clu, f"aln_{level}", "-a", "-e", "inf", "--min-seq-id", "0", "-c", "0",
               "--alignment-mode", "3")
            mm("convertalis", "db", "db", f"aln_{level}", f"aln_{level}.tsv", "--format-output", _ALN_FORMAT)
            _tsv_to_parquet(w / f"aln_{level}.tsv", w / f"aln_{level}.parquet", _ALN_COLUMNS, line_numbers=True)
        mm.log.close()

        con = duckdb.connect(config={"threads": 1, "memory_limit": limit, "temp_directory": str(w / "duckdb"),
                                     "preserve_insertion_order": True})
        try:
            aligned = {lvl: _pick_alignments(con, w, lvl) for lvl in LEVELS}
            counts = {lvl: _cluster_members(con, w, lvl, params_hash(lvl)) for lvl in LEVELS}
            mem_spec, aln_spec = table("cluster_member"), table("cluster_member_aln")
            order = sorted(LEVELS)  # `level` leads both sort keys

            def members() -> Iterator[pa.RecordBatch]:
                for lvl in order:
                    yield from io.duck_batches(con, (
                        f"SELECT unhex(member) AS member_uid, '{lvl}' AS level, cluster_uid, unhex(rep) AS rep_uid, "
                        f"'{params_hash(lvl)}' AS cluster_params_hash FROM read_parquet({_sql_path(w / f'mem_{lvl}.parquet')}) "
                        "ORDER BY cluster_uid, member_uid"), mem_spec)

            def alignments() -> Iterator[pa.RecordBatch]:
                # convertalis: query = representative, target = member; coordinates 1-based inclusive
                for lvl in order:
                    yield from io.duck_batches(con, (
                        f"SELECT unhex(target) AS member_uid, '{lvl}' AS level, unhex(query) AS rep_uid, "
                        "tstart - 1 AS qstart, tend AS qend, qstart - 1 AS tstart, qend AS tend, fident AS pident, "
                        f"alnlen, encode(cigar) AS cigar FROM read_parquet({_sql_path(w / f'pick_{lvl}.parquet')}) "
                        "ORDER BY rep_uid, member_uid"), aln_spec)

            io.write_parts(members(), mem_spec, out / f"l2/cluster_member/release={release}", EST_MEMBER_ROW,
                           max_file_bytes)
            io.write_parts(alignments(), aln_spec, out / f"l2/cluster_member_aln/release={release}", EST_ALN_ROW,
                           max_file_bytes)
        finally:
            con.close()

    # the release's own ledger (ledger/release=<R>/step=cluster, card B-28): NR90 sees every protein, NR50
    # the NR90 representatives
    n_in = {"nr90": n_prot, "nr50": counts["nr90"]}
    tools = provenance.tool_versions("mmseqs2", "pyarrow")
    ledger.write_release(out, release, "cluster",
                 sorted((ledger.row(lvl.encode(), "level", n_in[lvl], counts[lvl], n_in[lvl] - counts[lvl], tools=tools)
                         for lvl in LEVELS), key=lambda r: (r["object_key"], r["reason_code"])))
    return {"proteins": n_prot, "clusters": counts, "aligned": aligned}
