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

Bounded memory (AGENTS §4.1): a genome of at least `THREAD_MIN_BASES` bases is a batch alone, and
shorter consecutive genomes are grouped into batches of at most `batch_bases` bases (lengths from
`contig_index`). A batch is called genome by genome, in the calling process or, with `workers` > 1,
a long genome in a worker thread and a batch of short ones in a worker process, at most `workers`
batches at once (card B-32). Each batch comes back as two Arrow record batches that are appended,
in the order the batches finish, to unsorted scratch Parquet files in a private
temporary directory; one DuckDB connection (`DUCK_MEMORY_LIMIT`, spilling into that directory,
one thread so that the output order is fixed) then sorts each table into its file, so the output
depends neither on the batches nor on the worker count. Nothing is published if a genome fails.

Genes of a contig are ordered by (start, end, strand); every table names a gene by its `feature_uid`
and coordinates, never by its position in that order (AGENTS §3.3). Genome, tier and assembly quality
are not copied onto the gene rows: they live in `contig_dim` (AGENTS §3.6).

A contig shared with a genome of an earlier stage (DEC-26, card C-20): ingest stored its bytes in the stage
that saw them first, which owns its genes, IGRs and tape. This stage reads it from that stage (`seq_stages`),
so a genome is still trained and called whole, but writes no row for it and records it as skipped (the
result's `skipped_contigs`, a ledger row with reason `owned_by:<stage>`). Which genome's calls a shared contig
carries therefore depends on which stage was built first; no other contig's rows do.

A contig in several genomes of this stage (DEC-27, card C-21a): the genome with the least `genome_uid` writes its
rows. Every genome is still called whole; the other genomes' rows for the contig are dropped before they reach the
scratch files, and the skip is recorded (`skipped_contigs` with `owner_genome`; a ledger row with reason
`owned_by:<stage>/<genome_uid>` and the counts of the contig's `ok` row). Within a stage the output does not depend
on the order of the genomes.

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
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pyrodigal
import pyrodigal_gv

from biobox import ids, io, layout, ledger, provenance
from biobox.l0seq.inputs import GenomeInput, stage_genomes
from biobox.l0seq.lookup import KeyLookup
from biobox.l0seq.store import OWNED_BY, ContigOwners, SeqStore, missing_contig, skipped_entry
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
# a genome at least this long is called alone in a worker thread, shorter ones in batches in worker
# processes (card B-32). Chosen; the reasoning is an estimate from the one-thread profile of the slice
# (measured, wilkcomp02): a 620 bp contig cost 0.93 ms in `find_genes`, which releases the GIL, and
# about 0.5 ms of Python that holds it, mostly reading the sequence. `find_genes` grows with length
# and the Python part barely does, so at 10 kb the Python part is about 5% (derived: 0.93 ms x 16
# against about 0.75 ms) and threads scale, while a process costs seconds of imports to start
THREAD_MIN_BASES = 10_000
DUCK_MEMORY_LIMIT = "2GB"  # the sorting connection's cap; larger sorts spill into the scratch directory
SKIP_REASON = OWNED_BY  # ledger reason of a contig another stage owns (DEC-26), followed by that stage
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


def _rows(calls_per_genome: Iterable[tuple[bytes, list[dict]]]
          ) -> tuple[pa.RecordBatch, pa.RecordBatch, list[str], int, pa.Array]:
    """The gene and the protein-candidate scratch rows of some genomes' calls, given as (genome_uid, calls), as two
    record batches, the modes used, the number of genes that wrap the origin, and each gene row's genome_uid."""
    cds = ids.KIND["cds"]
    gcols: list[list] = [[] for _ in _GENE_SCRATCH]
    pus: list[bytes] = []
    aas: list[str] = []
    guids: list[bytes] = []
    modes: set[str] = set()
    n_wrap = 0
    for guid, calls in calls_per_genome:
        guids += [guid] * len(calls)
        for c in calls:
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
    return genes, cands, sorted(modes), n_wrap, pa.array(guids, type=pa.binary(16))


def _guid(genome: GenomeInput) -> bytes:
    return ids.genome_uid(list(genome.contig_uids))


def _call_batch(store: SeqStore, batch: list[GenomeInput]) -> tuple:
    """`_rows` of a batch of genomes, called one genome at a time in input order. Runs in a worker
    process, or in the caller's when `workers` is 1."""
    return _rows((_guid(g), call_genome(store, g)) for g in batch)


def _call_read(genome: GenomeInput, seqs: dict[bytes, str]) -> tuple:
    """`_rows` of one genome whose sequences the calling thread has read. Runs in a worker thread."""
    return _rows([(_guid(genome), _call(genome, seqs))])


def _worker_main(root: str, stages: list[str]) -> None:
    """A worker process: open the store (for the caller's `stages`, card B-28) once, then read
    `(seq, batch)` pickles from stdin until it closes and answer each with `("ok", seq, result)` or
    `("error", seq, exception)` on stdout."""
    store = SeqStore(Path(root), stages)
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
    """Up to `n` worker processes, each `python -c` running `_worker_main`, started when a batch finds
    none idle (a process costs about 2.8 s of imports, measured: B-32 measurement file). Not
    `multiprocessing`: its spawn and forkserver start methods run the caller's main script again in
    every worker (the stress drivers have no main guard), and fork is unsafe in a process that has
    threads. A batch goes to a worker's stdin and its result comes back on stdout as a pickle; one
    reader thread per worker drains stdout into the shared queue `q` as `(i, message)`, so neither
    side can block on a full pipe."""

    def __init__(self, n: int, root: Path, stages: list[str], q: queue.SimpleQueue):
        src = str(Path(__file__).resolve().parents[2])  # the directory that holds this `biobox`
        self._boot = (f"import sys; sys.path.insert(0, {src!r}); from biobox.l1.genecall import _worker_main; "
                      f"_worker_main({str(root)!r}, {[str(s) for s in stages]!r})")
        self._n, self._q = n, q
        self._procs: list[subprocess.Popen] = []
        self._readers: list[threading.Thread] = []
        self._idle: list[int] = []

    def _drain(self, i: int, p: subprocess.Popen) -> None:
        try:
            while True:
                self._q.put((i, pickle.load(p.stdout)))
        except BaseException as e:  # noqa: BLE001 - EOF or a torn stream: the worker is gone
            self._q.put((i, ("eof", None, e)))

    def send(self, seq: int, batch: list[GenomeInput]) -> None:
        """Give a batch to an idle worker, starting one if none is idle (at most `n`)."""
        if not self._idle:
            if len(self._procs) >= self._n:
                raise RuntimeError("no idle gene-call worker")  # the caller keeps at most n busy
            i = len(self._procs)
            p = subprocess.Popen([sys.executable, "-c", self._boot], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            self._procs.append(p)
            t = threading.Thread(target=self._drain, args=(i, p), daemon=True, name=f"genecall-reader-{i}")
            self._readers.append(t)
            t.start()
            self._idle.append(i)
        i = self._idle.pop()
        p = self._procs[i]
        try:
            pickle.dump((seq, batch), p.stdin, protocol=pickle.HIGHEST_PROTOCOL)
            p.stdin.flush()
        except BrokenPipeError:
            raise RuntimeError(f"gene-call worker {i} exited with code {p.wait()}") from None

    def answered(self, i: int, kind: str) -> None:
        """Worker `i` sent message `kind`: it is idle again, or gone."""
        if kind == "eof":
            raise RuntimeError(f"gene-call worker {i} exited with code {self._procs[i].wait()} before answering")
        self._idle.append(i)

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


def _genome_bases(store: SeqStore, lengths: KeyLookup, genome: GenomeInput) -> int:
    """Bases of a genome's default contig variants, from the store's `contig_index` files (one key
    lookup per contig). A contig the store does not hold raises KeyError naming the stage on disk that
    owns it, if any (DEC-26)."""
    n = 0
    for cu in genome.contig_uids:
        rows = [r for r in lengths.get(cu) if r["is_default_variant"]]
        if not rows:
            raise KeyError(missing_contig(store.root, cu, store.stages))
        n += rows[0]["contig_len"]
    return n


def _owned(genomes: Iterable[GenomeInput], owners: ContigOwners, skipped: dict[bytes, str],
           drop: set[tuple[bytes, bytes]]) -> Iterator[GenomeInput]:
    """The genomes, unchanged; each contig that another stage (DEC-26) or another genome of this stage (DEC-27) owns
    is added to `skipped` with its owner, and (genome_uid, contig_uid) to `drop`, before its genome is handed on."""
    for g in genomes:
        guid = _guid(g)
        for cu in g.contig_uids:
            s = owners.owner(cu, guid)
            if s != owners.stage:
                skipped[cu] = s
                drop.add((guid, cu))
        yield g


def _drop(genes: pa.RecordBatch, cands: pa.RecordBatch, modes: list[str], n_wrap: int, guids: pa.Array,
          skip: set[tuple[bytes, bytes]]) -> tuple[pa.RecordBatch, pa.RecordBatch, list[str], int]:
    """A batch's rows without the genes that `skip` names by (genome_uid, contig_uid): a contig another stage owns
    (DEC-26) or, in a genome other than its owner, a contig another genome of this stage owns (DEC-27). `guids` is
    each gene row's genome; the gene and candidate rows are parallel, one per gene."""
    if not skip:
        return genes, cands, modes, n_wrap
    keep = pa.array([k not in skip for k in zip(guids.to_pylist(), genes.column("contig_uid").to_pylist())],
                    pa.bool_())  # only for a batch that holds a skipped contig
    genes, cands = genes.filter(keep), cands.filter(keep)
    return (genes, cands, sorted(set(genes.column("mode").to_pylist())),
            int(pc.sum(genes.column("wraps_origin").cast(pa.int64())).as_py() or 0))


def _batches(genomes: Iterable[GenomeInput], bases_of, batch_bases: int, thread_min_bases: int
             ) -> Iterator[tuple[bool, list[GenomeInput]]]:
    """`(long, batch)` in input order. A genome of at least `thread_min_bases` bases is a batch alone
    (`long`); shorter consecutive genomes are grouped, cut before the genome that would take a batch
    past `batch_bases` bases (a genome larger than that is a batch alone)."""
    batch: list[GenomeInput] = []
    held = 0
    for g in genomes:
        n = bases_of(g)
        if n >= thread_min_bases or (batch and held + n > batch_bases):
            if batch:
                yield False, batch
            batch, held = [], 0
        if n >= thread_min_bases:
            yield True, [g]
            continue
        batch.append(g)
        held += n
    if batch:
        yield False, batch


def _called(store: SeqStore, root: Path, batches: Iterable[tuple[bool, list[GenomeInput]]], workers: int
            ) -> Iterator[tuple[list[GenomeInput], tuple]]:
    """(batch, `_rows` result); at most `workers` batches are called at once, and batches are read
    from `batches` lazily, one for each that finishes.

    Card B-32: a long genome (`long`) is called in a worker thread, its sequences read here in the
    calling thread so that the store is never shared between threads: pyrodigal releases the GIL in
    `find_genes`, which dominates a long genome's cost. Batches of short genomes go to worker
    processes (`_Workers`), each with its own `SeqStore`: in threads, sequence reads and per-gene
    Python work held the GIL and the workers starved. With workers, results come in the order they
    finish, so that a long genome does not hold back the batches behind it; the output does not
    depend on that order, because `_publish` sorts every table on a key that is unique when each
    contig's rows come from one genome (`_drop` keeps only its owner's; MR-01). A worker's exception is raised
    here; the processes are then killed and the threads finish their genome."""
    if workers <= 1:
        for _, b in batches:
            yield b, _call_batch(store, b)
        return
    q: queue.SimpleQueue = queue.SimpleQueue()  # (worker index, or None for a thread; message)
    procs = _Workers(workers, root, store.stages, q)
    threads: ThreadPoolExecutor | None = None

    def threaded(seq: int):
        def put(f) -> None:
            if not f.cancelled():
                e = f.exception()
                q.put((None, ("ok", seq, f.result()) if e is None else ("error", seq, e)))
        return put

    finished = False
    try:
        it = iter(batches)
        sent: dict[int, list[GenomeInput]] = {}
        n_sent = 0
        more = True
        while True:
            while more and len(sent) < workers:
                nxt = next(it, None)
                if nxt is None:
                    more = False
                    break
                long, b = nxt
                if long:
                    if threads is None:
                        threads = ThreadPoolExecutor(max_workers=workers, thread_name_prefix="genecall")
                    threads.submit(_call_read, b[0], _sequences(store, b[0])).add_done_callback(threaded(n_sent))
                else:
                    procs.send(n_sent, b)
                sent[n_sent] = b
                n_sent += 1
            if not sent:
                break
            i, (kind, seq, value) = q.get()
            if i is not None:
                procs.answered(i, kind)
            if kind == "error":
                raise value
            yield sent.pop(seq), value
        finished = True
    finally:
        procs.close(kill=not finished)
        if threads is not None:
            threads.shutdown(wait=True, cancel_futures=True)


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
        workers: int = 1, seq_stages: list[str] | None = None, out: Path | None = None,
        batch_bases: int = BATCH_BASES) -> dict:
    """Call genes and write the stage's tables. `prefix_width` (0, 1, 2) overrides the size rule.

    `out` (card C-08): the output root, default `root`. Sequences and `genome_input` are read from
    `root`; every table and the ledger are written under `out` (`$BIOBOX_STAGING` under run-step).
    Refuses, before any work, when one of the stage's gene-call partitions already holds a file in
    `root` or `out`; `tmp` may lie under neither.

    `genomes`: any iterable, consumed once; None reads the stage's `genome_input` table. `tmp`: the
    directory that holds the private scratch directory (system temp by default; never the data root).
    `workers`: batches called at once, in worker threads or processes; `batch_bases`: bases of
    consecutive short genomes sent to a worker process as one task (a larger genome is a task alone).
    The output bytes depend on neither. `max_file_bytes`: every table is cut into parts of at most
    this many estimated bytes (one row above it is a part alone).
    `seq_stages`: the stages whose sequences the genomes may name (card B-28): the stage itself (the
    default) and, when its contigs repeat sequences stored before it, the earlier stages that hold them.
    A contig in none of them raises KeyError naming the stages read and the stage on disk that owns it.

    A contig stored first by another stage (DEC-26, card C-20): that stage owns its genes. Its genome is
    still called whole (single mode trains on every contig of the genome), but no row is written for the
    contig; it is listed in the result (`skipped_contigs`: contig uid and owning stage, by contig uid)
    and in the ledger, one row per contig with reason `owned_by:<stage>` and every count 0.

    A contig in several genomes of this stage (DEC-27, card C-21a): the genome with the least `genome_uid`
    writes its rows (`inputs.shared_contigs`; given `genomes`, from them). The others are called whole too,
    and their rows for it are dropped; it is listed in `skipped_contigs` with `owner_genome`, and has a
    second ledger row, reason `owned_by:<stage>/<genome_uid>`, with the counts of its `ok` row (L-01).
    """
    root = Path(root)
    out = root if out is None else Path(out)
    if workers < 1:
        raise ValueError(f"workers must be at least 1, got {workers}")
    if batch_bases < 1:
        raise ValueError(f"batch_bases must be at least 1, got {batch_bases}")
    ledger.check_scratch(tmp, root, out)
    gcv = f"gene_call_version={version}/stage={stage}"
    ledger.refuse_written(root, out, [f"l1/feature/producer={version}/stage={stage}", f"l1/feature_protein/{gcv}",
                                      f"l1/protein_occurrence/{gcv}", f"l0seq/protein_seq/{gcv}",
                                      f"ledger/stage={stage}/step=genecall"], f"gene calls of stage {stage}")
    stages = list(dict.fromkeys([stage, *([] if seq_stages is None else seq_stages)]))
    store = SeqStore(root, stages)
    scratch = Path(tempfile.mkdtemp(prefix="biobox-genecall-", dir=tmp))
    try:
        genomes, shared = stage_genomes(root, stage, genomes, tmp=scratch)
        gs = _Scratch(scratch / _SCRATCH_FILES["genes"], _GENE_SCRATCH)
        cs = _Scratch(scratch / _SCRATCH_FILES["candidates"], _CANDIDATE_SCRATCH)
        ks = _Scratch(scratch / _SCRATCH_FILES["contigs"], _CONTIG_SCRATCH)
        n_genes = n_wrap = 0
        modes: set[str] = set()
        # contig -> owner (`ContigOwners.owner`); few: 4 of 54,265 contigs at 50 Gbp were shared at all (C-20)
        skipped: dict[bytes, str] = {}
        drop: set[tuple[bytes, bytes]] = set()  # (genome_uid, contig_uid) whose rows are not written
        lengths = KeyLookup(store.index_files, "contig_uid", ["contig_len", "is_default_variant"])
        owned = _owned(genomes, ContigOwners(root, stage, stages, shared), skipped, drop)
        batches = _batches(owned, lambda g: _genome_bases(store, lengths, g), batch_bases, THREAD_MIN_BASES)
        calls_iter = _called(store, root, batches, workers)
        try:
            for batch, (genes, cands, batch_modes, batch_wrap, guids) in calls_iter:
                pairs = [(_guid(g), cu) for g in batch for cu in g.contig_uids]
                skip = {p for p in pairs if p in drop}
                genes, cands, batch_modes, batch_wrap = _drop(genes, cands, batch_modes, batch_wrap, guids, skip)
                ks.add(pa.RecordBatch.from_arrays(
                    [pa.array([p[1] for p in pairs if p not in skip], type=pa.binary(16))],
                    schema=_CONTIG_SCRATCH))
                gs.add(genes)
                cs.add(cands)
                modes.update(batch_modes)
                n_wrap += batch_wrap
                n_genes += genes.num_rows
        finally:
            calls_iter.close()
            for s in (gs, cs, ks):
                s.close()
        n_prot = _publish(out, stage, version, scratch, prefix_width, max_file_bytes, skipped)
        return {"genes": n_genes, "distinct_proteins": n_prot, "modes": sorted(modes), "wrapping_genes": n_wrap,
                "skipped_contigs": [skipped_entry(cu, s) for cu, s in sorted(skipped.items())]}
    finally:
        store.close()
        shutil.rmtree(scratch, ignore_errors=True)


def _publish(root: Path, stage: str, version: str, scratch: Path, prefix_width: int | None,
             max_file_bytes: int, skipped: dict[bytes, str] | None = None) -> int:
    """Sort the scratch rows into the stage's files, each table cut into parts at `max_file_bytes`;
    returns the number of distinct proteins. The ledger rows stream from DuckDB to the writer; a contig
    in `skipped` (owned by another stage, DEC-26) has one ledger row, reason `owned_by:<stage>`, counts 0;
    one owned by another genome of this stage (DEC-27) has, beside its `ok` row, a row with reason
    `owned_by:<stage>/<genome_uid>` and the same counts."""
    skipped = skipped or {}
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb"), "preserve_insertion_order": True})
    try:
        con.register("skipped", pa.table({"contig_uid": pa.array(list(skipped), pa.binary(16)),
                                          "reason": pa.array([f"{SKIP_REASON}{s}" for s in skipped.values()],
                                                             pa.string())}))
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
                ", ok AS MATERIALIZED (SELECT c.contig_uid, count(g.contig_uid) AS n, "
                "count(*) FILTER (WHERE g.first) AS k FROM contigs c LEFT JOIN g USING (contig_uid) "
                "GROUP BY c.contig_uid) "
                "SELECT contig_uid, n, k, reason FROM (SELECT contig_uid, n, k, 'ok' AS reason FROM ok "
                # a skip within the stage repeats its contig's counts (L-01); another stage's contig has none here
                "UNION ALL SELECT s.contig_uid, coalesce(ok.n, 0), coalesce(ok.k, 0), s.reason FROM skipped s "
                "LEFT JOIN ok USING (contig_uid)) ORDER BY contig_uid, reason"
            ).to_arrow_reader(io.ROW_GROUP_ROWS)
            try:
                for b in reader:
                    for cu, n, k, why in zip(*(b.column(i).to_pylist() for i in range(4))):
                        yield ledger.row(cu, "contig", n, k, n - k, reason=why, tools=tools)
            finally:
                reader.close()

        ledger.write(root, stage, "genecall", ledger_rows())
    finally:
        con.close()
    return n_prot
