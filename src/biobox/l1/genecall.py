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

    l1/feature/producer=<v>/stage=<S>/part-0.parquet                       (card B-05; kind cds)
    l1/feature_protein/gene_call_version=<v>/stage=<S>/part-0.parquet      (card B-05)
    l1/protein_occurrence/gene_call_version=<v>/stage=<S>/part-0.parquet
    l1/slot_protein/gene_call_version=<v>/stage=<S>/part-0.parquet         (legacy; B-17 removes it)
    l0seq/protein_seq/gene_call_version=<v>/stage=<S>/protein_prefix=<all|h|hh>/part-<k>.parquet  (biobox.layout)
    ledger/stage=<S>/step=genecall/part-0.parquet

Bounded memory (AGENTS §4.1): genomes are called one at a time (or `workers` at a time, at most
`2 * workers` in flight) and their genes appended to unsorted scratch Parquet files in a private
temporary directory; one DuckDB connection (`DUCK_MEMORY_LIMIT`, spilling into that directory,
one thread so that the output order is fixed) then sorts each table into its file. Nothing is
published if a genome fails.

Coordinates are 0-based half-open (pyrodigal reports 1-based inclusive).
`dist_to_contig_end_bp` is 0 for genes running off a linear contig end and
equals the contig length on circular contigs, which have no edge.
"""

from __future__ import annotations

import shutil
import tempfile
from collections import deque
from collections.abc import Iterable, Iterator
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pyrodigal
import pyrodigal_gv

from biobox import ids, io, layout, ledger, provenance
from biobox.l0seq.inputs import GenomeInput, iter_genome_input
from biobox.l0seq.store import SeqStore
from biobox.registry.schemas import table

__all__ = ["GENE_CALL_VERSION", "GenomeInput", "call_genome", "run"]

PROTEIN_BYTES_PER_AA = 0.6  # zstd estimate (data_layer_and_tests.md §2A.4); sizes buckets only
SINGLE_MODE_MIN_BP = 100_000  # Prodigal: below this, use anonymous/meta mode
WRAP_FLANK = 50_000
GENE_CALL_VERSION = f"pyrodigal-{pyrodigal.__version__}-gv-{pyrodigal_gv.__version__}-p2"
KINDS = ("prokaryote", "virus", "plasmid")
SCRATCH_FLUSH_ROWS = 65_536  # rows a scratch writer buffers before it writes a row group
DUCK_MEMORY_LIMIT = "2GB"  # the sorting connection's cap; larger sorts spill into the scratch directory

# unsorted per-gene scratch rows, appended genome by genome
_GENE_SCRATCH = pa.schema([
    ("contig_uid", pa.binary(16)), ("slot_ordinal", pa.uint32()), ("feature_uid", pa.binary(16)),
    ("start", pa.int64()), ("end", pa.int64()), ("strand", pa.int8()), ("protein_uid", pa.binary(16)),
    ("genetic_code", pa.uint8()), ("genome_uid", pa.binary(16)), ("dist_to_contig_end_bp", pa.uint32()),
    ("assembly_quality", pa.string()), ("tier", pa.uint8()), ("is_species_rep", pa.bool_()),
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
    guid = ids.genome_uid(list(genome.contig_uids))
    calls: list[dict] = []
    for cu in sorted(seqs):
        seq = seqs[cu]
        circular = cu in genome.circular
        genes = sorted(_genes(finder, seq, circular), key=lambda t: (t[0], t[1], t[2]))
        for ordinal, (start, end, strand, aa, code, partial) in enumerate(genes):
            if circular:
                dist = len(seq)
            else:
                dist = 0 if partial else min(start, len(seq) - end)
            calls.append({
                "contig_uid": cu,
                "slot_ordinal": ordinal,
                "genome_uid": guid,
                "start": start,
                "end": end,
                "strand": strand,
                "aa": aa,
                "genetic_code": code,
                "mode": mode,
                "wraps_origin": end > len(seq),
                "dist_to_contig_end_bp": dist,
                "assembly_quality": genome.assembly_quality,
                "tier": genome.tier,
                "is_species_rep": genome.is_species_rep,
            })
    return calls


def call_genome(store: SeqStore, genome: GenomeInput) -> list[dict]:
    return _call(genome, _sequences(store, genome))


def _calls_in_order(store: SeqStore, genomes: Iterable[GenomeInput], workers: int
                    ) -> Iterator[tuple[GenomeInput, list[dict]]]:
    """(genome, calls) in input order. With `workers` > 1, genomes are called in threads (pyrodigal
    releases the GIL in `find_genes`), submitted lazily so that at most `2 * workers` are in flight
    or waiting. Sequences are read here, in the calling thread, so the store is never shared
    between threads. A worker's exception is raised here; the genomes not yet started are cancelled."""
    if workers <= 1:
        for g in genomes:
            yield g, call_genome(store, g)
        return
    ex = ThreadPoolExecutor(max_workers=workers, thread_name_prefix="genecall")
    pending: deque = deque()
    try:
        for g in genomes:
            pending.append((g, ex.submit(_call, g, _sequences(store, g))))
            if len(pending) >= 2 * workers:
                g0, fut = pending.popleft()
                yield g0, fut.result()
        while pending:
            g0, fut = pending.popleft()
            yield g0, fut.result()
    finally:
        ex.shutdown(wait=True, cancel_futures=True)


class _Scratch:
    """An unsorted Parquet file written in row groups of `SCRATCH_FLUSH_ROWS` rows."""

    def __init__(self, path: Path, schema: pa.Schema):
        self.path, self._schema = path, schema
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
        workers: int = 1) -> dict:
    """Call genes and write the stage's tables. `prefix_width` (0, 1, 2) overrides the size rule.

    `genomes`: any iterable, consumed once; None reads the stage's `genome_input` table. `tmp`: the
    directory that holds the private scratch directory (system temp by default; never the data root).
    `workers`: genomes called at once; the output bytes do not depend on it.
    """
    root = Path(root)
    if workers < 1:
        raise ValueError(f"workers must be at least 1, got {workers}")
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
        cds = ids.KIND["cds"]
        calls_iter = _calls_in_order(store, genomes, workers)
        try:
            for g, calls in calls_iter:
                for cu in g.contig_uids:
                    ks.add((cu,))
                for c in calls:
                    pu = ids.protein_uid(c["aa"])
                    fu = ids.feature_uid(c["contig_uid"], c["start"], c["end"], c["strand"], cds)
                    gs.add((c["contig_uid"], c["slot_ordinal"], fu, c["start"], c["end"], c["strand"], pu,
                            c["genetic_code"], c["genome_uid"], c["dist_to_contig_end_bp"], c["assembly_quality"],
                            c["tier"], c["is_species_rep"], c["mode"], c["wraps_origin"]))
                    cs.add((pu, c["aa"]))
                    modes.add(c["mode"])
                    n_wrap += c["wraps_origin"]
                n_genes += len(calls)
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
    """Sort the scratch rows into the stage's files; returns the number of distinct proteins. Each
    table is one `part-0.parquet` for now; splitting a large stage's tables into parts is a later card."""
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb"), "preserve_insertion_order": True})
    try:
        src = {k: f"read_parquet({_sql_path(scratch / name)})" for k, name in _SCRATCH_FILES.items()}
        con.execute(f"CREATE VIEW genes AS SELECT * FROM {src['genes']}")
        con.execute(f"CREATE VIEW candidates AS SELECT * FROM {src['candidates']}")
        con.execute(f"CREATE VIEW contigs AS SELECT DISTINCT contig_uid FROM {src['contigs']}")
        gcv = f"gene_call_version={version}/stage={stage}"
        outputs = [
            ("feature", f"l1/feature/producer={version}/stage={stage}",
             f'SELECT feature_uid, contig_uid, start, "end", strand, {ids.KIND["cds"]}::UTINYINT AS kind FROM genes '
             'ORDER BY contig_uid, start, "end", strand'),
            ("feature_protein", f"l1/feature_protein/{gcv}",
             'SELECT contig_uid, feature_uid, start, "end", strand, protein_uid, genetic_code, \'none\' AS edit_type '
             'FROM genes ORDER BY contig_uid, start, "end", strand'),
            ("protein_occurrence", f"l1/protein_occurrence/{gcv}",
             'SELECT protein_uid, contig_uid, slot_ordinal, genome_uid, start, "end", strand, '
             'NULL::UINTEGER AS tax_order, is_species_rep, tier, assembly_quality, dist_to_contig_end_bp, feature_uid '
             'FROM genes '
             'ORDER BY protein_uid, contig_uid, slot_ordinal'),
            ("slot_protein", f"l1/slot_protein/{gcv}",
             "SELECT contig_uid, slot_ordinal, protein_uid, genetic_code, 'none' AS edit_type FROM genes "
             "ORDER BY contig_uid, slot_ordinal"),
        ]
        for name, rel, sql in outputs:
            io.write_table(io.duck_batches(con, sql, table(name)), table(name), root / rel / "part-0.parquet")

        distinct = "(SELECT DISTINCT protein_uid, aa FROM candidates)"
        (n_prot, n_aa), = con.execute(f"SELECT count(*), coalesce(sum(length(aa)), 0) FROM {distinct}").fetchall()
        # the stage's estimated protein bytes, sum over proteins of 40 + PROTEIN_BYTES_PER_AA * length
        est = 40 * n_prot + PROTEIN_BYTES_PER_AA * n_aa
        width = layout.prefix_width(est) if prefix_width is None else prefix_width
        _write_protein_seq(
            con, f"SELECT protein_uid, unhex(sha256(aa)) AS protein_digest, length(aa) AS length, aa AS seq "
                 f"FROM {distinct} ORDER BY protein_uid",
            root / f"l0seq/protein_seq/{gcv}", width, max_file_bytes)

        # ledger: per contig, genes in; kept = first gene of its protein sequence in (contig, slot) order
        reader = con.sql(
            "WITH g AS (SELECT contig_uid, row_number() OVER (PARTITION BY protein_uid ORDER BY contig_uid, slot_ordinal)"
            " = 1 AS first FROM genes) "
            "SELECT c.contig_uid, count(g.contig_uid) AS n, count(*) FILTER (WHERE g.first) AS k "
            "FROM contigs c LEFT JOIN g USING (contig_uid) GROUP BY c.contig_uid ORDER BY c.contig_uid"
        ).to_arrow_reader(io.ROW_GROUP_ROWS)
        tools = provenance.tool_versions("pyrodigal", "pyrodigal-gv", "pyarrow")
        # ledger.write takes a list: one small row per contig is held here
        rows = [ledger.row(cu, "contig", n, k, n - k, tools=tools)
                for b in reader for cu, n, k in zip(*(b.column(i).to_pylist() for i in range(3)))]
        ledger.write(root, stage, "genecall", rows)
    finally:
        con.close()
    return n_prot
