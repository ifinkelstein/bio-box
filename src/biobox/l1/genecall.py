"""Gene calling (operation O-CALL) and the protein tables it feeds.

Policy p2 (approved 2026-10-29, PLAN.md §3.1; field practice in
planning/DECISIONS.md DEC-13):

- viral contigs: pyrodigal-gv in meta mode (prodigal-gv models for giant
  viruses and alternative codes, including code 15), code taken per gene;
- standalone plasmids: meta mode;
- prokaryotic genomes >= 100 kb: single mode trained on all contigs of the
  genome together with the genome's code; below 100 kb: meta mode;
- circular contigs: genes may wrap across the origin. Genes are called on the
  contig with `WRAP_FLANK` bp of the other end added on both sides and kept
  when their start lies on the contig, so each gene appears once with full
  context; a wrapping gene has `end > contig_len`.

Writes, per gene-calling version and stage (paths under `<root>/`), all through `biobox.io.write_table`:

    l1/feature/producer=<v>/stage=<S>/part-<k>.parquet                       (card B-05; kind cds)
    l1/feature_protein/gene_call_version=<v>/stage=<S>/part-<k>.parquet      (card B-05)
    l1/protein_occurrence/gene_call_version=<v>/stage=<S>/part-<k>.parquet
    l0seq/protein_seq/gene_call_version=<v>/stage=<S>/protein_prefix=<all|h|hh>/part-<k>.parquet  (biobox.layout)
    ledger/stage=<S>/step=genecall/part-0.parquet

Every table is cut into parts at `max_file_bytes` of estimated size (`io.write_parts`, card B-19); a
small stage is one `part-0.parquet` per table.

Bounded memory (AGENTS §4.1): consecutive genomes are grouped into batches of at most
`batch_bases` bases (a larger genome is a batch alone; lengths from `contig_index`), and a batch is
called genome by genome, in the calling process or, with `workers` > 1, in one of `workers`
worker processes, at most `2 * workers` batches in flight (card B-32). Each batch comes back as
two Arrow record batches that are appended to unsorted scratch Parquet files in a private
temporary directory; one DuckDB connection (`DUCK_MEMORY_LIMIT`, spilling into that directory,
one thread so that the output order is fixed) then sorts each table into its file, so the output
depends neither on the batches nor on the worker count. Nothing is published if a genome fails.

Genes of a contig are ordered by (start, end, strand); every table names a gene by its `feature_uid`
and coordinates, never by its position in that order (AGENTS §3.3). Genome, tier and assembly quality
are not copied onto the gene rows: they live in `contig_dim` (AGENTS §3.6).

Coordinates are 0-based half-open (pyrodigal reports 1-based inclusive).
`dist_to_contig_end_bp` is 0 for genes running off a linear contig end and
equals the contig length on circular contigs, which have no edge.
"""

from __future__ import annotations

import pickle
import queue
import shutil
import subprocess
import sys
import tempfile
import threading
from collections.abc import Iterable, Iterator
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pyrodigal
import pyrodigal_gv

from biobox import ids, io, layout, ledger, provenance
from biobox.l0seq.inputs import GenomeInput, iter_genome_input
from biobox.l0seq.lookup import KeyLookup
from biobox.l0seq.store import SeqStore
from biobox.registry.schemas import table

__all__ = ["GENE_CALL_VERSION", "GenomeInput", "call_genome", "run"]

PROTEIN_BYTES_PER_AA = 0.6  # zstd estimate (data_layer_and_tests.md §2A.4); sizes buckets only
SINGLE_MODE_MIN_BP = 100_000  # Prodigal: below this, use anonymous/meta mode
WRAP_FLANK = 50_000
GENE_CALL_VERSION = f"pyrodigal-{pyrodigal.__version__}-gv-{pyrodigal_gv.__version__}-p2"
KINDS = ("prokaryote", "virus", "plasmid")
SCRATCH_FLUSH_ROWS = 65_536  # rows a scratch writer buffers before it writes a row group
# bases of consecutive genomes in one worker task (card B-32): chosen, not tuned. Derived: about 2.4 s
# of meta-mode calling per task (one thread called 61.95 Mbp of short contigs in 146 s, measured on
# wilkcomp02, planning/measurements/2026-11-10_genecall_short_contigs.md), so the per-task cost is
# small, and a 100,000-contig slice still makes 62 tasks for 8 workers
BATCH_BASES = 1_000_000
DUCK_MEMORY_LIMIT = "2GB"  # the sorting connection's cap; larger sorts spill into the scratch directory
# Estimated zstd bytes per row; they only cut parts at `max_file_bytes` (card B-19). Measured values
# are from planning/measurements/2026-11-09_design_stress_test.md (§3 golden set, §6a S1 pilot), taken
# before `feature_uid` (16 random bytes, which do not compress) was added to these tables.
EST_FEATURE_ROW = 24  # estimate: feature_uid 16 + delta-coded coordinates; not measured
EST_FEATURE_PROTEIN_ROW = 40  # estimate: feature_uid and protein_uid, 16 random bytes each; not measured
# estimate: 32.2 B/gene measured (§3; 30.5 in §6a) + 16 for feature_uid, before card B-17 removed the
# ordinal, genome and dimension columns; not re-measured
EST_OCCURRENCE_ROW = 48

# unsorted per-gene scratch rows, appended genome by genome
_GENE_SCRATCH = pa.schema([
    ("contig_uid", pa.binary(16)), ("feature_uid", pa.binary(16)),
    ("start", pa.int64()), ("end", pa.int64()), ("strand", pa.int8()), ("protein_uid", pa.binary(16)),
    ("genetic_code", pa.uint8()), ("dist_to_contig_end_bp", pa.uint32()),
    ("mode", pa.string()), ("wraps_origin", pa.bool_()),
])
_CANDIDATE_SCRATCH = pa.schema([("protein_uid", pa.binary(16)), ("aa", pa.string())])
_CONTIG_SCRATCH = pa.schema([("contig_uid", pa.binary(16))])  # every contig, also those without genes
_SCRATCH_FILES = {"genes": "genes.parquet", "candidates": "candidates.parquet", "contigs": "contigs.parquet"}


def _finder(genome: GenomeInput, seqs: dict[bytes, str]):
    if genome.kind not in KINDS:
        raise ValueError(f"unknown genome kind {genome.kind!r}")
    if genome.kind == "virus":
        return pyrodigal_gv.ViralGeneFinder(meta=True), "gv-meta"
    if genome.kind == "plasmid" or sum(len(s) for s in seqs.values()) < SINGLE_MODE_MIN_BP:
        return pyrodigal.GeneFinder(meta=True), "meta"
    finder = pyrodigal.GeneFinder(meta=False)
    finder.train(*seqs.values(), translation_table=genome.genetic_code)
    return finder, "single"


def _genes(finder, seq: str, circular: bool):
    """Yield (start, end, strand, aa, code, partial) with 0-based half-open coordinates."""
    n = len(seq)
    if not circular:
        for g in finder.find_genes(seq):
            yield g.begin - 1, g.end, g.strand, g.translate(include_stop=False), g.translation_table, bool(
                g.partial_begin or g.partial_end
            )
        return
    flank = min(n, WRAP_FLANK)
    ext = seq[n - flank :] + seq + seq[:flank]
    for g in finder.find_genes(ext):
        start = g.begin - 1 - flank
        if not 0 <= start < n or g.partial_begin or g.partial_end:
            continue
        end = g.end - flank
        yield start, end, g.strand, g.translate(include_stop=False), g.translation_table, False


def _sequences(store: SeqStore, genome: GenomeInput) -> dict[bytes, str]:
    return {cu: store.get(cu).upper().decode("ascii") for cu in genome.contig_uids}


def _call(genome: GenomeInput, seqs: dict[bytes, str]) -> list[dict]:
    finder, mode = _finder(genome, seqs)
    calls: list[dict] = []
    for cu in sorted(seqs):
        seq = seqs[cu]
        circular = cu in genome.circular
        genes = sorted(_genes(finder, seq, circular), key=lambda t: (t[0], t[1], t[2]))
        for start, end, strand, aa, code, partial in genes:
            if circular:
                dist = len(seq)
            else:
                dist = 0 if partial else min(start, len(seq) - end)
            calls.append({
                "contig_uid": cu,
                "start": start,
                "end": end,
                "strand": strand,
                "aa": aa,
                "genetic_code": code,
                "mode": mode,
                "wraps_origin": end > len(seq),
                "dist_to_contig_end_bp": dist,
            })
    return calls


def call_genome(store: SeqStore, genome: GenomeInput) -> list[dict]:
    return _call(genome, _sequences(store, genome))


def _call_batch(store: SeqStore, batch: list[GenomeInput]) -> tuple[pa.RecordBatch, pa.RecordBatch, list[str], int]:
    """Call the genes of a batch of genomes, one genome at a time, in input order: the gene and the
    protein-candidate scratch rows as two record batches, the modes used and the number of genes that
    wrap the origin. Runs in a worker process (or in the caller's when `workers` is 1)."""
    cds = ids.KIND["cds"]
    gcols: list[list] = [[] for _ in _GENE_SCRATCH]
    pus: list[bytes] = []
    aas: list[str] = []
    modes: set[str] = set()
    n_wrap = 0
    for genome in batch:
        for c in call_genome(store, genome):
            pu = ids.protein_uid(c["aa"])
            fu = ids.feature_uid(c["contig_uid"], c["start"], c["end"], c["strand"], cds)
            row = (c["contig_uid"], fu, c["start"], c["end"], c["strand"], pu, c["genetic_code"],
                   c["dist_to_contig_end_bp"], c["mode"], c["wraps_origin"])
            for col, v in zip(gcols, row):
                col.append(v)
            pus.append(pu)
            aas.append(c["aa"])
            modes.add(c["mode"])
            n_wrap += c["wraps_origin"]
    genes = pa.RecordBatch.from_arrays([pa.array(c, type=f.type) for c, f in zip(gcols, _GENE_SCRATCH)],
                                       schema=_GENE_SCRATCH)
    cands = pa.RecordBatch.from_arrays([pa.array(pus, type=pa.binary(16)), pa.array(aas, type=pa.string())],
                                       schema=_CANDIDATE_SCRATCH)
    return genes, cands, sorted(modes), n_wrap


def _worker_main(root: str) -> None:
    """A worker process: open the store once, then read `(seq, batch)` pickles from stdin until it
    closes and answer each with `("ok", seq, result)` or `("error", seq, exception)` on stdout."""
    store = SeqStore(Path(root))
    inp, out = sys.stdin.buffer, sys.stdout.buffer
    sys.stdout = sys.stderr  # a stray print must not corrupt the result stream
    try:
        while True:
            try:
                seq, batch = pickle.load(inp)
            except EOFError:
                return
            try:
                msg = ("ok", seq, _call_batch(store, batch))
            except Exception as e:  # noqa: BLE001 - raised again in the parent
                msg = ("error", seq, e)
            try:
                data = pickle.dumps(msg, protocol=pickle.HIGHEST_PROTOCOL)
            except Exception as e:  # noqa: BLE001 - an exception that does not pickle
                data = pickle.dumps(("error", seq, RuntimeError(f"{type(msg[2]).__name__}: {msg[2]} ({e})")))
            out.write(data)
            out.flush()
    finally:
        store.close()


class _Workers:
    """`n` worker processes, each `python -c` running `_worker_main`. Not `multiprocessing`: its spawn
    and forkserver start methods run the caller's main script again in every worker (the stress
    drivers have no main guard), and fork is unsafe in a process that has threads. Tasks go to a
    worker's stdin and results come back on its stdout as pickles; one reader thread per worker
    drains stdout into a shared queue, so neither side can block on a full pipe."""

    def __init__(self, n: int, root: Path):
        src = str(Path(__file__).resolve().parents[2])  # the directory that holds this `biobox`
        boot = (f"import sys; sys.path.insert(0, {src!r}); "
                f"from biobox.l1.genecall import _worker_main; _worker_main({str(root)!r})")
        self._q: queue.SimpleQueue = queue.SimpleQueue()
        self._procs = [subprocess.Popen([sys.executable, "-c", boot], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                       for _ in range(n)]
        self._readers = [threading.Thread(target=self._drain, args=(i, p), daemon=True, name=f"genecall-reader-{i}")
                         for i, p in enumerate(self._procs)]
        for t in self._readers:
            t.start()

    def _drain(self, i: int, p: subprocess.Popen) -> None:
        try:
            while True:
                self._q.put((i, pickle.load(p.stdout)))
        except BaseException as e:  # noqa: BLE001 - EOF or a torn stream: the worker is gone
            self._q.put((i, ("eof", None, e)))

    def send(self, i: int, seq: int, batch: list[GenomeInput]) -> None:
        p = self._procs[i]
        try:
            pickle.dump((seq, batch), p.stdin, protocol=pickle.HIGHEST_PROTOCOL)
            p.stdin.flush()
        except BrokenPipeError:
            raise RuntimeError(f"gene-call worker {i} exited with code {p.wait()}") from None

    def receive(self) -> tuple[int, int, tuple]:
        """(worker, seq, result) of the next batch any worker finishes; a worker's exception is raised."""
        i, (kind, seq, value) = self._q.get()
        if kind == "ok":
            return i, seq, value
        if kind == "error":
            raise value
        raise RuntimeError(f"gene-call worker {i} exited with code {self._procs[i].wait()} before answering")

    def close(self, kill: bool) -> None:
        for p in self._procs:
            if kill:
                p.kill()
            try:
                p.stdin.close()
            except OSError:
                pass
        for p in self._procs:
            p.wait()
        for t in self._readers:
            t.join()
        for p in self._procs:
            p.stdout.close()


def _genome_bases(lengths: KeyLookup, genome: GenomeInput) -> int:
    """Bases of a genome's default contig variants, from `contig_index` (one key lookup per contig)."""
    n = 0
    for cu in genome.contig_uids:
        rows = [r for r in lengths.get(cu) if r["is_default_variant"]]
        if not rows:
            raise KeyError(cu)
        n += rows[0]["contig_len"]
    return n


def _batches(genomes: Iterable[GenomeInput], bases_of, batch_bases: int) -> Iterator[list[GenomeInput]]:
    """Consecutive genomes in input order, cut before the genome that would take a batch past
    `batch_bases` bases; a genome larger than that is a batch alone."""
    batch: list[GenomeInput] = []
    held = 0
    for g in genomes:
        n = bases_of(g)
        if batch and held + n > batch_bases:
            yield batch
            batch, held = [], 0
        batch.append(g)
        held += n
    if batch:
        yield batch


def _calls_in_order(store: SeqStore, root: Path, batches: Iterable[list[GenomeInput]], workers: int
                    ) -> Iterator[tuple[list[GenomeInput], tuple]]:
    """(batch, `_call_batch` result) in input order. With `workers` > 1, batches are called in worker
    processes (`_Workers`), each with its own `SeqStore` (card B-32: in threads, sequence reads and
    per-gene Python work held the GIL and the workers starved on short contigs). A worker takes one
    batch at a time; batches are read from `batches` lazily, so that at most `2 * workers` are in
    flight or done and waiting for their turn. A worker's exception is raised here, and the workers
    are then killed."""
    if workers <= 1:
        for b in batches:
            yield b, _call_batch(store, b)
        return
    pool = _Workers(workers, root)
    finished = False
    try:
        it = iter(batches)
        sent: dict[int, list[GenomeInput]] = {}
        done: dict[int, tuple] = {}
        idle = list(range(workers))
        n_sent = n_out = 0
        more = True
        while True:
            while more and idle and n_sent - n_out < 2 * workers:
                b = next(it, None)
                if b is None:
                    more = False
                    break
                pool.send(idle.pop(), n_sent, b)
                sent[n_sent] = b
                n_sent += 1
            if n_out in done:
                yield sent.pop(n_out), done.pop(n_out)
                n_out += 1
            elif n_out == n_sent:
                break
            else:
                i, seq, res = pool.receive()
                idle.append(i)
                done[seq] = res
        finished = True
    finally:
        pool.close(kill=not finished)


class _Scratch:
    """An unsorted Parquet file written in row groups of `SCRATCH_FLUSH_ROWS` rows."""

    def __init__(self, path: Path, schema: pa.Schema):
        self.path, self._schema = path, schema
        self._held: list[pa.RecordBatch] = []
        self._n = 0
        self._writer = pq.ParquetWriter(path, schema, compression="zstd", write_statistics=False)

    def add(self, b: pa.RecordBatch) -> None:
        if b.num_rows:
            self._held.append(b)
            self._n += b.num_rows
        while self._n >= SCRATCH_FLUSH_ROWS:
            t = pa.Table.from_batches(self._held, schema=self._schema)
            self._writer.write_table(t.slice(0, SCRATCH_FLUSH_ROWS), row_group_size=SCRATCH_FLUSH_ROWS)
            rest = t.slice(SCRATCH_FLUSH_ROWS)
            self._held, self._n = rest.to_batches(), rest.num_rows

    def close(self) -> None:
        if self._n:
            self._writer.write_table(pa.Table.from_batches(self._held, schema=self._schema),
                                     row_group_size=SCRATCH_FLUSH_ROWS)
            self._held, self._n = [], 0
        self._writer.close()


def _sql_path(p: Path) -> str:
    return "'" + str(p).replace("'", "''") + "'"


class _Cursor:
    """Sorted record batches, consumed from the front in slices."""

    def __init__(self, batches: Iterator[pa.RecordBatch]):
        self._it, self._b, self._off = batches, None, 0

    def head(self) -> pa.RecordBatch | None:
        while self._b is None or self._off >= self._b.num_rows:
            self._b, self._off = next(self._it, None), 0
            if self._b is None:
                return None
        return self._b.slice(self._off)

    def take(self, n: int) -> None:
        self._off += n


def _buckets(b: pa.RecordBatch, width: int) -> np.ndarray:
    """`layout.bucket` of every row as an integer: 0 at width 0, the first hex digit or byte otherwise."""
    if width == 0:
        return np.zeros(b.num_rows, dtype=np.int64)
    col = b.column("protein_uid")
    first = np.frombuffer(col.buffers()[1], dtype=np.uint8)[col.offset * 16 : (col.offset + len(col)) * 16 : 16]
    return (first >> 4 if width == 1 else first).astype(np.int64)


def _bucket_name(code: int, width: int) -> str:
    return layout.ALL if width == 0 else bytes([code if width == 2 else code << 4]).hex()[:width]


def _part(rows: _Cursor, code: int, width: int, max_file_bytes: int) -> Iterator[pa.RecordBatch]:
    """The rows of one `protein_seq` file: bucket `code`, cut where `layout.cut` cuts (same float sums)."""
    held = 0.0
    while (b := rows.head()) is not None:
        codes = _buckets(b, width)
        same = int(np.argmax(codes != code)) if (codes != code).any() else b.num_rows
        est = 40 + PROTEIN_BYTES_PER_AA * b.column("length").to_numpy()[:same].astype(np.float64)
        cs = np.cumsum(np.concatenate([[held], est]))  # cs[i] = bytes held before row i
        cut = np.flatnonzero((cs[:-1] != 0) & (cs[1:] > max_file_bytes))
        n = int(cut[0]) if len(cut) else same
        if n:
            yield b.slice(0, n)
            rows.take(n)
        if n < b.num_rows:
            return
        held = float(cs[n])


def _write_protein_seq(con, sql: str, base: Path, width: int, max_file_bytes: int) -> None:
    spec = table("protein_seq")
    rows = _Cursor(io.duck_batches(con, sql, spec))
    while (b := rows.head()) is not None:
        code = int(_buckets(b.slice(0, 1), width)[0])
        k = 0
        while (b := rows.head()) is not None and int(_buckets(b.slice(0, 1), width)[0]) == code:
            io.write_table(_part(rows, code, width, max_file_bytes), spec,
                           base / f"protein_prefix={_bucket_name(code, width)}/part-{k}.parquet")
            k += 1


def run(root: Path, stage: str, genomes: Iterable[GenomeInput] | None = None, version: str = GENE_CALL_VERSION,
        prefix_width: int | None = None, max_file_bytes: int = layout.MAX_FILE_BYTES, tmp=None,
        workers: int = 1, batch_bases: int = BATCH_BASES) -> dict:
    """Call genes and write the stage's tables. `prefix_width` (0, 1, 2) overrides the size rule.

    `genomes`: any iterable, consumed once; None reads the stage's `genome_input` table. `tmp`: the
    directory that holds the private scratch directory (system temp by default; never the data root).
    `workers`: worker processes calling genes at once; `batch_bases`: bases of consecutive genomes
    sent to a worker as one task (a larger genome is a task alone). The output bytes depend on
    neither. `max_file_bytes`: every table is cut into parts of at most this many estimated bytes
    (one row above it is a part alone).
    """
    root = Path(root)
    if workers < 1:
        raise ValueError(f"workers must be at least 1, got {workers}")
    if batch_bases < 1:
        raise ValueError(f"batch_bases must be at least 1, got {batch_bases}")
    if tmp is not None and Path(tmp).resolve().is_relative_to(root.resolve()):
        raise ValueError(f"{tmp}: the scratch directory may not lie under the data root")
    store = SeqStore(root)
    scratch = Path(tempfile.mkdtemp(prefix="biobox-genecall-", dir=tmp))
    try:
        if genomes is None:
            genomes = iter_genome_input(root, stage, tmp=scratch)
        gs = _Scratch(scratch / _SCRATCH_FILES["genes"], _GENE_SCRATCH)
        cs = _Scratch(scratch / _SCRATCH_FILES["candidates"], _CANDIDATE_SCRATCH)
        ks = _Scratch(scratch / _SCRATCH_FILES["contigs"], _CONTIG_SCRATCH)
        n_genes = n_wrap = 0
        modes: set[str] = set()
        lengths = KeyLookup(store.index_files, "contig_uid", ["contig_len", "is_default_variant"])
        batches = _batches(genomes, lambda g: _genome_bases(lengths, g), batch_bases)
        calls_iter = _calls_in_order(store, root, batches, workers)
        try:
            for batch, (genes, cands, batch_modes, batch_wrap) in calls_iter:
                ks.add(pa.RecordBatch.from_arrays(
                    [pa.array([cu for g in batch for cu in g.contig_uids], type=pa.binary(16))], schema=_CONTIG_SCRATCH))
                gs.add(genes)
                cs.add(cands)
                modes.update(batch_modes)
                n_wrap += batch_wrap
                n_genes += genes.num_rows
        finally:
            calls_iter.close()
            for s in (gs, cs, ks):
                s.close()
        n_prot = _publish(root, stage, version, scratch, prefix_width, max_file_bytes)
        return {"genes": n_genes, "distinct_proteins": n_prot, "modes": sorted(modes), "wrapping_genes": n_wrap}
    finally:
        store.close()
        shutil.rmtree(scratch, ignore_errors=True)


def _publish(root: Path, stage: str, version: str, scratch: Path, prefix_width: int | None,
             max_file_bytes: int) -> int:
    """Sort the scratch rows into the stage's files, each table cut into parts at `max_file_bytes`;
    returns the number of distinct proteins. The ledger rows stream from DuckDB to the writer."""
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb"), "preserve_insertion_order": True})
    try:
        src = {k: f"read_parquet({_sql_path(scratch / name)})" for k, name in _SCRATCH_FILES.items()}
        con.execute(f"CREATE VIEW genes AS SELECT * FROM {src['genes']}")
        con.execute(f"CREATE VIEW candidates AS SELECT * FROM {src['candidates']}")
        con.execute(f"CREATE VIEW contigs AS SELECT DISTINCT contig_uid FROM {src['contigs']}")
        gcv = f"gene_call_version={version}/stage={stage}"
        outputs = [
            ("feature", f"l1/feature/producer={version}/stage={stage}", EST_FEATURE_ROW,
             f'SELECT feature_uid, contig_uid, start, "end", strand, {ids.KIND["cds"]}::UTINYINT AS kind FROM genes '
             'ORDER BY contig_uid, start, "end", strand'),
            ("feature_protein", f"l1/feature_protein/{gcv}", EST_FEATURE_PROTEIN_ROW,
             'SELECT contig_uid, feature_uid, start, "end", strand, protein_uid, genetic_code, \'none\' AS edit_type '
             'FROM genes ORDER BY contig_uid, start, "end", strand'),
            ("protein_occurrence", f"l1/protein_occurrence/{gcv}", EST_OCCURRENCE_ROW,
             'SELECT protein_uid, contig_uid, start, "end", strand, dist_to_contig_end_bp, feature_uid '
             'FROM genes '
             'ORDER BY protein_uid, contig_uid, start, "end", strand'),
        ]
        for name, rel, est, sql in outputs:
            io.write_parts(io.duck_batches(con, sql, table(name)), table(name), root / rel, est, max_file_bytes)

        distinct = "(SELECT DISTINCT protein_uid, aa FROM candidates)"
        (n_prot, n_aa), = con.execute(f"SELECT count(*), coalesce(sum(length(aa)), 0) FROM {distinct}").fetchall()
        # the stage's estimated protein bytes, sum over proteins of 40 + PROTEIN_BYTES_PER_AA * length
        est = 40 * n_prot + PROTEIN_BYTES_PER_AA * n_aa
        width = layout.prefix_width(est) if prefix_width is None else prefix_width
        _write_protein_seq(
            con, f"SELECT protein_uid, unhex(sha256(aa)) AS protein_digest, length(aa) AS length, aa AS seq "
                 f"FROM {distinct} ORDER BY protein_uid",
            root / f"l0seq/protein_seq/{gcv}", width, max_file_bytes)

        # ledger: per contig, genes in; kept = first gene of its protein sequence in (contig, start, end,
        # strand) order; rows stream from DuckDB in contig_uid order (the ledger's sort key) to the writer
        tools = provenance.tool_versions("pyrodigal", "pyrodigal-gv", "pyarrow")

        def ledger_rows() -> Iterator[dict]:
            reader = con.sql(
                "WITH g AS (SELECT contig_uid, row_number() OVER (PARTITION BY protein_uid ORDER BY contig_uid, "
                'start, "end", strand) = 1 AS first FROM genes) '
                "SELECT c.contig_uid, count(g.contig_uid) AS n, count(*) FILTER (WHERE g.first) AS k "
                "FROM contigs c LEFT JOIN g USING (contig_uid) GROUP BY c.contig_uid ORDER BY c.contig_uid"
            ).to_arrow_reader(io.ROW_GROUP_ROWS)
            try:
                for b in reader:
                    for cu, n, k in zip(*(b.column(i).to_pylist() for i in range(3))):
                        yield ledger.row(cu, "contig", n, k, n - k, tools=tools)
            finally:
                reader.close()

        ledger.write(root, stage, "genecall", ledger_rows())
    finally:
        con.close()
    return n_prot
