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

A genome all of whose contigs another stage (DEC-26) or another genome of this stage (DEC-27) owns is not called
(card B-48, item 5): every row of its calls would be dropped, so calling it changes no byte. Its contigs are skipped
and recorded exactly as before (each in `skipped_contigs`; the ledger rows `owned_by:<stage>` with counts 0, or
`owned_by:<stage>/<genome_uid>` with the counts of the owner's `ok` row), and nothing else of it was ever written:
not a gene row, a contig row, a mode or a wrapping gene. On stage PL1 that is 7,046 genomes (*derived*, card C-26).

Coordinates are 0-based half-open (pyrodigal reports 1-based inclusive).
`dist_to_contig_end_bp` is 0 for genes running off a linear contig end and
equals the contig length on circular contigs, which have no edge.

Gene calling in batches (card B-48). The step is estimated before it starts (`plan`, `genome_cost_us`): each genome
at the rate of the mode the code calls it in (`mode_of`, the rule of `_finder`) and the topology of each contig, per
base the caller reads (a circular contig with its wrap flanks, `called_bases`), plus each process's fixed cost and the
writing's cost per gene; the rates and their sources are the constants below. `build-stage` runs the step in one
process, as before and with the same argv, when that estimate is within `stage.SHARD_MARGIN` of the budget; else in
batches: the genomes in `genome_uid` order (the order `iter_genome_input` streams them in) cut into the fewest
consecutive runs of equal estimated cost each of which fits the margin (`cut_plan`), one process per batch
(`bio-box stage genecall batch`), then one merge (`... merge`) that is the `run-step` step named `genecall`.

What the one-process call does that could depend on the whole stage, and why batches do not change it:

- Training. Single mode trains on the contigs of one genome (`_finder`); meta and gv-meta train nothing. No step
  trains on more than one genome, and a batch holds whole genomes, so each genome's calls are those of one process.
- Which genome writes a contig. DEC-26: the stage whose `contig_index` holds it (the same lookup in every batch).
  DEC-27: the least `genome_uid` among the stage's genomes that list it, from `shared_contigs` over the whole stage's
  `genome_input`, which every batch reads whole. So the rows dropped and the skips recorded do not depend on the cut.
- Row order. Every table is sorted by `_publish` on a key that is unique when one genome writes each contig (the
  scratch rows' order already varies with the worker count). The merge runs `_publish` once over the scratch rows of
  every batch: the same rows, so the same order.
- File cuts (`max_file_bytes`) and the protein bucket width. Both are computed by `_publish` from the sorted stream
  and from the distinct proteins of the whole stage, which the merge sees whole; `io.write_parts` cuts by the rows'
  estimated sizes alone, never by how its input was batched.
- Ids. `feature_uid`, `protein_uid` and `genome_uid` are functions of a row's own fields.
- The ledger. Its `n_kept` (the first gene of each protein over the stage, in (contig, start, end, strand) order) is a
  window over every gene of the stage, computed by the merge over all batches; its skip rows come from the union of
  the batches' skips, and a contig's owner is the same in every batch. Rows are ordered by (contig, reason).
- The result. `genes` and `wrapping_genes` are sums over the batches, `modes` their union (each the modes of the rows
  written), `distinct_proteins` the merge's count, `skipped_contigs` the union sorted by contig uid.

So the merge writes the one-process call's files, ledger and result, byte for byte, whatever the number of batches,
the number at once, the threads of each and the order they finish in. A batch's rows go only into the plan's work
directory in scratch, never into the root: its scratch rows, its skips and `done.json` (counts, its CPU-minutes, wall
time and peak RSS, its estimate) into `.batch-<j>.tmp-<pid>-<hex>`, renamed to `batch-<j>` at the end. The merge takes
only a batch directory whose `done.json` names the plan's digest of its genomes; a batch stopped by SIGTERM, SIGHUP or
SIGINT removes its temporary directory, one killed by SIGKILL leaves it dot-named, and the next driver removes it
before it launches anything. A batch with output is not run again.

Where a batch's rows live once its step is published (card B-50, decided on the card: its "Amended 2026-11-15"). In
scratch, and nowhere else. Each piece of the step is a `run-step` step of the stage, recorded in `runs/stage=<S>/` like
ingest's batches: `genecall.plan` (`stage genecall plan --plan P`: the plan build-stage made, checked against the
stage's `genome_input`, published as it is) and `genecall.batch-<j>` (`stage genecall batch`: the batch called into
scratch as above). What each stages, because `run-step` publishes only a file a registered table or artifact claims and
refuses an empty staging, is one small JSON of the artifact `genecall_plan`, at
`l1/genecall_plan/stage=<S>/plan.json` (the plan file's bytes) and `.../batch-<j>.json` (`batch_record`: the fields
of its `done.json` that do not depend on the run, which are its counts, modes, skips, genome digest and estimate, and the
plan's sha256; its CPU-minutes, wall time and peak memory are the run record's). The artifact is not in any manifest
and is `small_ok`, as `ingest_plan` is. So the root gains per batch two small files, the record and its run record (370 to
393 and 1,401 bytes on card B-48's fixture: *measured*, wilkcomp02, 2026-11-15, card B-50), and per stage the plan's
two; no gene-call row is ever written there twice. The rows
are intermediate: the merge reads them from scratch and supersedes them with the stage's tables. Copying them into the
root, as ingest's batches publish their shards, would add every batch's rows beside the merged tables (about the
gene-call bytes of the stage again) for files nothing reads after the merge.
That satisfies AGENTS §3.9 and §3.10: every root file is written once under a temporary name and renamed into place by
`run-step`; nothing published is ever rewritten (the records are written once, `_write_record`); the merge reads the
batches named by the published plan, and checks each one's scratch output against its published record. Scratch lost
after a batch was published (a reboot clears /tmp) is made again: `build-stage` runs that batch again outside
`run-step` (`stage genecall batch --check-record`), which must reproduce the published record (its counts and genome
digest: the rows' order in scratch depends on the workers, so not their sha256), else it stops before the merge.
A direct call of `stage genecall plan|batch` with neither `--out` nor `$BIOBOX_STAGING` writes no record: only scratch,
as card B-48's commands did.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import pickle
import queue
import resource
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

import duckdb
from blake3 import blake3
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pyrodigal
import pyrodigal_gv

from biobox import ids, io, layout, ledger, provenance
from biobox.l0seq.inputs import GenomeInput, iter_genome_input, stage_genomes
from biobox.l0seq.lookup import KeyLookup
from biobox.l0seq.store import OWNED_BY, ContigOwners, SeqStore, missing_contig, skipped_entry
from biobox.registry.schemas import ARTIFACTS, table

__all__ = ["GENE_CALL_VERSION", "GenomeInput", "call_genome", "merge", "plan", "publish_plan", "run", "run_batch"]

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

# --- the estimate of the step (cards B-48, B-53; planning/measurements/2026-11-15_genecall_batches.md,
# planning/measurements/2026-11-17_crispr_genecall_estimates.md) ------------------------------------------------------
# The calls alone, CPU-seconds per Mbp of the bases the caller reads (`called_bases`: a circular contig with its wrap
# flanks), by mode and topology (`mode_of`; True = circular). Each measured rate is *measured* on wilkcomp02, one
# thread, with the step's own call path (`_call`) timed alone, on every 30th genome or stored contig: single mode
# 0.7941 on stage P6 (72 genomes, 351.38 Mbp, 2026-11-15; P2's 0.7838 is the lower of the two), gv-meta linear 4.2019
# on stage S1 (908 genomes, 55.45 Mbp, 2026-11-15; H1's 3.7795 and V1's 3.7222 on 4,833 linear virus genomes, 137.47
# Mbp, 2026-11-16, are lower), gv-meta circular 4.3762 on stage V1 (778 circular UHGV v1.0 virus genomes,
# `/tmp/if2453/s2/inputs/V1/`, 34,564,134 contig bp, 88,404,110 bp read with 50 kb wrap flanks, 386.87 CPU-s = 4.37615
# CPU-s per Mbp read, 11.1929 per contig Mbp; 2026-11-16, card C-39 `c39_est.py`, `V1/est.json`, paired with
# `called_bases` in card B-53's note), meta 2.9450 on PLSDB's linear records and 4.4205 on its circular ones (180 and
# 1,874 records of stage PL1, 17.14 and 295.46 Mbp called, 2026-11-15; card C-26's sample, `c26_rate.py`). Each is the
# largest measurement, and each is rounded up. Only circular single mode is not measured (no built stage has a circular
# contig in single mode): it is the linear rate times meta mode's circular-to-linear ratio, 1.5011 (*derived*:
# 4.4205 / 2.9449; an *estimate* for circular single mode).
CALL_S_PER_MBP = {("single", False): 0.7941, ("single", True): 0.7941 * 1.5011,
                  ("gv-meta", False): 4.2019, ("gv-meta", True): 4.3762,
                  ("meta", False): 2.9450, ("meta", True): 4.4205}
# The whole step over its calls alone: the reading of the sequences, the rows, the scratch files, the sorting and
# writing of every table, and the workers' contention, as one factor on the calls. *Derived*, the largest of the eight
# built stages' run records (wilkcomp02, 2026-11-10 to 15; card B-48's note): each record's cpu_min over its calls at
# the rates above, 1.0635 (H1) to 1.2252 (S1), rounded up to 1.23. With it every record's estimate is at it or above
# it, by at most 15.8% (H1).
STEP_OVERHEAD = 1.23
# Every process's own cost: interpreter and imports, the store's footers, `genome_input` streamed and its shared
# contigs. *Estimate* (C-10 measured 0.07 CPU-minutes for a step that does almost nothing; a batch streams the whole
# stage's genome_input once more), checked against the batches of stage S1 (card B-48's note).
PROCESS_S = 15.0
# The merge of a step in batches: `_publish` over every batch's rows, per gene written. *Derived* from the merge of stage
# S1 in 4 batches (run record, wilkcomp02, 2026-11-15, card B-48's note: 0.432 CPU-minutes for 2,398,087 genes, so
# (25.95 s - PROCESS_S) / 2,398,087 = 4.57 us a gene), raised by a fifth and rounded up.
PUBLISH_S_PER_GENE = 0.0000055
# Genes per Mbp of contig, by mode, for the merge's estimate (the genes are known only after calling). *Measured*:
# gv-meta 1,441.0 on S1 (2,398,087 genes, 1,664.11 Mbp: its run record's feature rows), single 917.4 on P5 (9,144,001
# genes, 9,967.90 Mbp), meta 1,079.0 on PLSDB's circular records (card C-26's sample). The largest of each mode.
GENES_PER_MBP = {"single": 917.4, "gv-meta": 1441.1, "meta": 1079.1}

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


def mode_of(kind: str, bases: int) -> str:
    """The mode a genome of `kind` and `bases` bases is called in (policy p2, module docstring)."""
    if kind not in KINDS:
        raise ValueError(f"unknown genome kind {kind!r}")
    if kind == "virus":
        return "gv-meta"
    if kind == "plasmid" or bases < SINGLE_MODE_MIN_BP:
        return "meta"
    return "single"


def _finder(genome: GenomeInput, seqs: dict[bytes, str]):
    mode = mode_of(genome.kind, sum(len(s) for s in seqs.values()))
    if mode == "gv-meta":
        return pyrodigal_gv.ViralGeneFinder(meta=True), "gv-meta"
    if mode == "meta":
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
    is added to `skipped` with its owner, and (genome_uid, contig_uid) to `drop`, before its genome is handed on. A
    genome all of whose contigs are so owned is not handed on (card B-48, item 5: every row of its calls would be
    dropped); its skips are recorded all the same."""
    for g in genomes:
        guid = _guid(g)
        mine = False
        for cu in g.contig_uids:
            s = owners.owner(cu, guid)
            if s != owners.stage:
                skipped[cu] = s
                drop.add((guid, cu))
            else:
                mine = True
        if mine:
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
    _refuse_written(root, out, stage, version)
    stages = list(dict.fromkeys([stage, *([] if seq_stages is None else seq_stages)]))
    store = SeqStore(root, stages)
    scratch = Path(tempfile.mkdtemp(prefix="biobox-genecall-", dir=tmp))
    try:
        genomes, shared = stage_genomes(root, stage, genomes, tmp=scratch)
        counts, skipped = _call_into(store, root, ContigOwners(root, stage, stages, shared), genomes, workers,
                                     batch_bases, scratch)
        files = {k: [scratch / name] for k, name in _SCRATCH_FILES.items()}
        n_prot = _publish(out, stage, version, files, scratch, prefix_width, max_file_bytes, skipped)
        return _result(counts, n_prot, skipped)
    finally:
        store.close()
        shutil.rmtree(scratch, ignore_errors=True)


def _refuse_written(root: Path, out: Path, stage: str, version: str) -> None:
    gcv = f"gene_call_version={version}/stage={stage}"
    ledger.refuse_written(root, out, [f"l1/feature/producer={version}/stage={stage}", f"l1/feature_protein/{gcv}",
                                      f"l1/protein_occurrence/{gcv}", f"l0seq/protein_seq/{gcv}",
                                      f"ledger/stage={stage}/step=genecall"], f"gene calls of stage {stage}")


def _result(counts: dict, n_prot: int, skipped: dict[bytes, str]) -> dict:
    """The step's result (the run's printed `result`): the same keys and order whether in one process or merged."""
    return {"genes": counts["genes"], "distinct_proteins": n_prot, "modes": sorted(counts["modes"]),
            "wrapping_genes": counts["wrapping_genes"],
            "skipped_contigs": [skipped_entry(cu, s) for cu, s in sorted(skipped.items())]}


def _call_into(store: SeqStore, root: Path, owners: ContigOwners, genomes: Iterable[GenomeInput], workers: int,
               batch_bases: int, scratch: Path) -> tuple[dict, dict[bytes, str]]:
    """Call `genomes` and append their rows to the three unsorted scratch files in `scratch` (`_SCRATCH_FILES`);
    return the counts (genes, wrapping genes and modes of the rows written) and the skipped contigs with their owners.
    The calling loop of `run`, and of a batch (card B-48)."""
    gs = _Scratch(scratch / _SCRATCH_FILES["genes"], _GENE_SCRATCH)
    cs = _Scratch(scratch / _SCRATCH_FILES["candidates"], _CANDIDATE_SCRATCH)
    ks = _Scratch(scratch / _SCRATCH_FILES["contigs"], _CONTIG_SCRATCH)
    n_genes = n_wrap = 0
    modes: set[str] = set()
    # contig -> owner (`ContigOwners.owner`); few: 4 of 54,265 contigs at 50 Gbp were shared at all (C-20)
    skipped: dict[bytes, str] = {}
    drop: set[tuple[bytes, bytes]] = set()  # (genome_uid, contig_uid) whose rows are not written
    lengths = KeyLookup(store.index_files, "contig_uid", ["contig_len", "is_default_variant"])
    owned = _owned(genomes, owners, skipped, drop)
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
    return {"genes": n_genes, "wrapping_genes": n_wrap, "modes": sorted(modes)}, skipped


def _publish(root: Path, stage: str, version: str, files: dict[str, list[Path]], scratch: Path,
             prefix_width: int | None, max_file_bytes: int, skipped: dict[bytes, str] | None = None) -> int:
    """Sort the scratch rows of `files` (each of `_SCRATCH_FILES`' keys: its scratch files, one in one process, one
    per batch in a merge, card B-48) into the stage's files, each table cut into parts at `max_file_bytes`;
    returns the number of distinct proteins. DuckDB spills under `scratch`. The ledger rows stream from DuckDB to the
    writer; a contig in `skipped` (owned by another stage, DEC-26) has one ledger row, reason `owned_by:<stage>`,
    counts 0; one owned by another genome of this stage (DEC-27) has, beside its `ok` row, a row with reason
    `owned_by:<stage>/<genome_uid>` and the same counts."""
    skipped = skipped or {}
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb"), "preserve_insertion_order": True})
    try:
        con.register("skipped", pa.table({"contig_uid": pa.array(list(skipped), pa.binary(16)),
                                          "reason": pa.array([f"{SKIP_REASON}{s}" for s in skipped.values()],
                                                             pa.string())}))
        src = {k: "read_parquet([" + ", ".join(_sql_path(p) for p in files[k]) + "], hive_partitioning = false)"
               for k in _SCRATCH_FILES}
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


# --- the estimate, the plan, the batches and the merge (card B-48; module docstring) --------------------------------

PLAN_FILE, DONE_FILE, SKIPPED_FILE = "plan.json", "done.json", "skipped.parquet"
WORK_BASE = "genecall"  # under the scratch directory: <scratch>/genecall/<stage>.<root key>/<plan key>/
_SKIPPED_SCHEMA = pa.schema([("contig_uid", pa.binary(16)), ("owner", pa.string())])
_STOP_CODES = {signal.SIGTERM: 143, signal.SIGHUP: 129}  # a batch stopped by a signal exits 128 + signum


class Stopped(SystemExit):
    """A batch ended by SIGTERM or SIGHUP, after its temporary directory was removed; the code is 128 + signum."""


def _log(msg: str) -> None:
    print(f"{datetime.now().astimezone().isoformat(timespec='seconds')} genecall: {msg}", file=sys.stderr, flush=True)


def called_bases(n: int, circular: bool) -> int:
    """Bases the caller reads for a contig of `n` bases: a circular one with `WRAP_FLANK` (or all of it) on each side
    (`_genes`)."""
    return n + 2 * min(n, WRAP_FLANK) if circular else n


def calls_cpu_s(mode: str, contigs: Iterable[tuple[int, bool]]) -> float:
    """CPU-seconds of the calls alone of a genome in `mode` with these (length, circular) contigs, at
    `CALL_S_PER_MBP`: what card C-26's sample and this card's samples measured."""
    return sum(CALL_S_PER_MBP[(mode, c)] * called_bases(n, c) for n, c in contigs) / 1e6


def genome_cost_us(kind: str, contigs: Iterable[tuple[int, bool]]) -> int:
    """The estimate of one genome in the step, in CPU-microseconds (rounded up): its calls times `STEP_OVERHEAD`."""
    contigs = list(contigs)
    mode = mode_of(kind, sum(n for n, _ in contigs))
    return int(-(-STEP_OVERHEAD * calls_cpu_s(mode, contigs) * 1e6 // 1))


def one_process_cpu_min(called: dict[tuple[str, bool], int]) -> float:
    """The one-process estimate of a stage whose called genomes read `called` bases (with the wrap flanks) by (mode,
    circular): a process's own cost and `STEP_OVERHEAD` times the calls (genome costs are linear in their bases)."""
    return (PROCESS_S + STEP_OVERHEAD * sum(CALL_S_PER_MBP[k] * b for k, b in called.items()) / 1e6) / 60


def genome_genes(kind: str, bases: int) -> float:
    """Genes a genome of `kind` and `bases` bases is expected to have (`GENES_PER_MBP`), for the merge's estimate."""
    return GENES_PER_MBP[mode_of(kind, bases)] * bases / 1e6


def merge_cpu_min(genes: float) -> float:
    """The merge's estimate: a process's own cost and `_publish` over `genes` genes."""
    return (PROCESS_S + PUBLISH_S_PER_GENE * genes) / 60


def _sql_files(files: Iterable[Path]) -> str:
    return "[" + ", ".join(_sql_path(Path(f)) for f in files) + "]"


def stage_items(root: Path, stage: str, stages: list[str], tmp: Path | None = None
                ) -> Callable[[], Iterator[tuple[bytes, int, int, float, bool]]]:
    """A maker of the stage's genomes as the plan reads them, one pass each call: (genome_uid, cost in CPU-us, bases,
    expected genes, called) in genome_uid order (the order `iter_genome_input` streams them in, and so the order a
    batch counts them in). One DuckDB query a pass: `genome_input` joined to the default variants of `stages`'
    `contig_index` (the stage that holds one owns the contig, DEC-26) and to the least genome of each contig within
    the stage (DEC-27); a genome none of whose contigs is its own is not called and costs nothing (item 5). A contig
    in none of `stages` raises KeyError naming the stage on disk that owns it, as the step does."""
    from biobox.l0seq.store import INDEX_DIR

    root = Path(root)
    gi = io.part_files(root / table("genome_input").location / f"stage={stage}")
    if not gi:
        raise FileNotFoundError(f"{root / table('genome_input').location}/stage={stage}: no genome_input for stage "
                                f"{stage}")
    ci = []
    for s in stages:
        files = io.part_files(root / INDEX_DIR / f"stage={s}")
        if files:
            lit = "'" + s.replace("'", "''") + "'"
            ci.append(f"SELECT contig_uid, contig_len, {lit} AS st FROM read_parquet({_sql_files(files)}, "
                      "hive_partitioning = false) WHERE is_default_variant")
    own = "'" + stage.replace("'", "''") + "'"
    sql = (f"WITH gi AS (SELECT genome_uid, contig_uid, is_circular, kind FROM read_parquet({_sql_files(gi)}, "
           "hive_partitioning = false)), "
           f"ci AS ({' UNION ALL '.join(ci) if ci else 'SELECT NULL::BLOB AS contig_uid, 0 AS contig_len, NULL AS st'}), "
           "least AS (SELECT contig_uid, min(genome_uid) AS g FROM gi GROUP BY contig_uid) "
           f"SELECT gi.genome_uid, gi.contig_uid, gi.is_circular, gi.kind, ci.contig_len, "
           f"coalesce(ci.st = {own} AND least.g = gi.genome_uid, false) AS mine FROM gi "
           "LEFT JOIN ci USING (contig_uid) JOIN least USING (contig_uid) ORDER BY gi.genome_uid, gi.contig_uid")

    def items() -> Iterator[tuple[bytes, int, int, float, bool]]:
        d = Path(tempfile.mkdtemp(prefix="biobox-genecall-plan-", dir=tmp))
        con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT, "temp_directory": str(d)})
        try:
            reader = con.sql(sql).to_arrow_reader(io.ROW_GROUP_ROWS)
            guid, rows = None, []
            for b in reader:
                for r in zip(*(b.column(i).to_pylist() for i in range(6))):
                    if r[0] != guid:
                        if rows:
                            yield _item(root, stages, guid, rows)
                        guid, rows = r[0], []
                    rows.append(r)
            if rows:
                yield _item(root, stages, guid, rows)
        finally:
            con.close()
            shutil.rmtree(d, ignore_errors=True)

    return items


def _item(root: Path, stages: list[str], guid: bytes, rows: list[tuple]) -> tuple[bytes, int, int, float, bool]:
    for _, cu, _, _, n, _ in rows:
        if n is None:
            raise KeyError(missing_contig(root, cu, stages))
    kind = rows[0][3]
    contigs = [(n, bool(c)) for _, _, c, _, n, _ in rows]
    bases = sum(n for n, _ in contigs)
    if not any(r[5] for r in rows):
        mode_of(kind, bases)  # an unknown kind is refused here as in the step
        return guid, 0, bases, 0.0, False
    return guid, genome_cost_us(kind, contigs), bases, genome_genes(kind, bases), True


def cut_plan(items: Callable[[], Iterable[tuple[bytes, int, int, float, bool]]], batch_cpu_min: float,
             batches: int | None = None) -> dict:
    """The batches of a step whose genomes `items()` gives, in order, as (genome_uid, cost in CPU-us, bases, expected
    genes, called), one pass per call (`stage_items`).

    The cut: with N batches, a genome starts batch k at the first genome whose midpoint, in cost from the first
    genome, is at k/N of the total or past it (`2 N before + N cost >= 2 k total`, exact in integers, as card W-35
    cuts contigs); a genome is never split. N is `batches` when given, else the fewest from ceil(total / limit) up
    whose every batch fits `limit` (`batch_cpu_min` less a process's own cost, `PROCESS_S`), one pass over the items
    per N tried; when none fits (a genome alone is over it), the first N tried. The plan's `bases` are every genome's,
    `bases_called` those of the genomes called. Each batch: its first genome's index, its genomes, the bases
    (`bases`) and expected genes of the ones it calls, its estimate (CPU-minutes, its own process
    included) and the BLAKE3-128 digest of its genome uids in order, which the batch checks. A pure function of the
    items, the budget and `batches`."""
    total = n = not_called = bases = bases_called = 0
    genes = 0.0
    big, big_uid = 0, None
    for guid, cost, b, ge, called in items():
        total, n, bases, genes = total + cost, n + 1, bases + b, genes + ge
        not_called += not called
        bases_called += b if called else 0
        if cost > big:
            big, big_uid = cost, guid
    limit = int(batch_cpu_min * 60e6 - PROCESS_S * 1e6)
    if batches is not None:
        if isinstance(batches, bool) or batches < 1:
            raise ValueError(f"batches must be at least 1, got {batches}")
        tries = [batches]
    else:
        first = max(1, -(-total // limit)) if limit > 0 else 1
        # a batch of the cut holds at most total / N + the largest genome (its genomes' midpoints lie within total / N),
        # so N = ceil(total / (limit - largest)) fits for certain: the search ends there; no N fits a genome over it
        tries = range(first, max(first, -(-total // (limit - big))) + 1) if big < limit else [first]
    cut = None
    for k in tries:
        got = _cut(items, k, total)
        if cut is None:
            cut = got
        if max(c["cost_us"] for c in got) <= limit:
            cut = got
            break
    out = []
    for j, c in enumerate(cut):
        out.append({"batch": j, "first": c["first"], "genomes": c["genomes"], "bases": c["bases"],
                    "genes_expected": round(c["genes"], 1),
                    "estimate_cpu_min": round((PROCESS_S * 1e6 + c["cost_us"]) / 60e6, 3),
                    "digest": c["digest"].hexdigest(length=16)})
    one = (PROCESS_S * 1e6 + total) / 60e6
    return {"card": "B-48", "genomes": n, "genomes_not_called": not_called, "bases": bases,
            "bases_called": bases_called, "genes_expected": round(genes, 1), "estimate_cpu_min": round(one, 3),
            "largest_genome": None if big_uid is None else {"genome_uid": big_uid.hex(),
                                                            "estimate_cpu_min": round(big / 60e6, 3)},
            "batch_cpu_min": batch_cpu_min, "batches": out, "merge_estimate_cpu_min": round(merge_cpu_min(genes), 3),
            "rates": {"call_s_per_mbp": {f"{m}/{'circular' if c else 'linear'}": r
                                         for (m, c), r in sorted(CALL_S_PER_MBP.items())},
                      "step_overhead": STEP_OVERHEAD, "process_s": PROCESS_S,
                      "publish_s_per_gene": PUBLISH_S_PER_GENE}}


def _cut(items, k: int, total: int) -> list[dict]:
    cut = [{"first": 0, "genomes": 0, "bases": 0, "genes": 0.0, "cost_us": 0, "digest": blake3()} for _ in range(k)]
    j, before, i = 0, 0, 0
    for guid, cost, b, ge, called in items():
        while j + 1 < k and 2 * k * before + k * cost >= 2 * (j + 1) * total:
            j += 1
            cut[j]["first"] = i
        c = cut[j]
        c["genomes"] += 1
        c["cost_us"] += cost
        if called:
            c["bases"] += b
            c["genes"] += ge
        c["digest"].update(guid)
        before += cost
        i += 1
    for x in cut[j + 1:]:
        x["first"] = i
    return cut


def plan(root: Path, stage: str, *, version: str = GENE_CALL_VERSION, seq_stages: list[str] | None = None,
         batch_cpu_min: float, batches: int | None = None, tmp=None) -> dict:
    """The estimate and the batches of the stage's gene calling (`cut_plan` over `stage_items`), with the stage, the
    gene-call version and the stages read (`seq_stages`, the stage first). Reads `genome_input` and `contig_index`
    only; writes nothing (DuckDB spills under `tmp`)."""
    stages = list(dict.fromkeys([stage, *([] if seq_stages is None else seq_stages)]))
    p = cut_plan(stage_items(Path(root), stage, stages, tmp), batch_cpu_min, batches)
    return {"stage": stage, "gene_call_version": version, "seq_stages": stages, **p}


def work_dir(scratch: Path, root: Path, stage: str, plan_: dict) -> Path:
    """Where a plan's batches live: `<scratch>/genecall/<stage>.<root key>/<plan key>/`, the keys BLAKE3 of the root's
    absolute path and of the plan, so a plan of another root, or another plan of this stage, never shares it."""
    rk = blake3(str(Path(root).resolve()).encode()).hexdigest(length=8)
    pk = blake3(json.dumps(plan_, sort_keys=True).encode()).hexdigest(length=16)
    return Path(scratch) / WORK_BASE / f"{stage}.{rk}" / pk


def write_plan(plan_: dict, work: Path) -> Path:
    """`plan_` as `<work>/plan.json` (sorted-key JSON), written under a temporary name and renamed; the same file is
    left as it is."""
    work.mkdir(parents=True, exist_ok=True)
    path = work / PLAN_FILE
    data = json.dumps(plan_, sort_keys=True, indent=1).encode()
    if not path.is_file() or path.read_bytes() != data:
        t = work / f".{PLAN_FILE}.tmp-{os.getpid()}"
        t.write_bytes(data)
        os.replace(t, path)
    return path


def read_plan(path: Path) -> dict:
    p = json.loads(Path(path).read_text())
    if p.get("card") != "B-48" or not isinstance(p.get("batches"), list):
        raise ValueError(f"{path}: not a gene-call plan (card B-48)")
    return p


def batch_name(j: int) -> str:
    return f"batch-{j:05d}"


def batch_done(work: Path, plan_: dict, j: int) -> dict | None:
    """The record of batch `j` whose output is in `work` and is of the plan's batch (its `done.json` names the plan's
    digest and genomes, and its files have the recorded sizes), else None."""
    d = Path(work) / batch_name(j)
    try:
        rec = json.loads((d / DONE_FILE).read_text())
    except (FileNotFoundError, ValueError):
        return None
    want = plan_["batches"][j]
    if rec.get("digest") != want["digest"] or rec.get("genomes") != want["genomes"] or rec.get("first") != want["first"]:
        return None
    for name, size in rec.get("files", {}).items():
        if not (d / name).is_file() or (d / name).stat().st_size != size:
            return None
    return rec


# --- the published records of a step in batches (card B-50; module docstring) -------------------------------------

RECORDS_DIR = ARTIFACTS["genecall_plan"].location  # l1/genecall_plan
RECORD_CARD = "B-50"
# the fields of a batch's done.json that do not depend on the run (its usage is the run record's), in its record
RECORD_FIELDS = ("batch", "first", "genomes", "digest", "genes", "wrapping_genes", "modes", "skipped", "estimate_cpu_min")


def records_dir(stage: str) -> str:
    """The root-relative directory of a stage's `genecall_plan` records."""
    return f"{RECORDS_DIR}/stage={stage}"


def plan_record_rel(stage: str) -> str:
    return f"{records_dir(stage)}/{PLAN_FILE}"


def batch_record_rel(stage: str, j: int) -> str:
    return f"{records_dir(stage)}/{batch_name(j)}.json"


def published_plan(root: Path, stage: str) -> tuple[bytes, dict] | None:
    """The stage's published plan record (its bytes and the plan), or None."""
    p = Path(root) / plan_record_rel(stage)
    if not p.is_file():
        return None
    raw = p.read_bytes()
    return raw, read_plan(p)


def published_batch(root: Path, stage: str, j: int) -> dict | None:
    """The published record of batch `j` of the stage, or None."""
    p = Path(root) / batch_record_rel(stage, j)
    return json.loads(p.read_text()) if p.is_file() else None


def batch_record(plan_: dict, plan_sha256: str, done: dict) -> dict:
    """What batch `done` publishes (`genecall_plan`'s batch-<j>.json): the plan's stage and version, the sha256 of the
    plan file's bytes, and `RECORD_FIELDS` of its done.json."""
    return {"card": RECORD_CARD, "stage": plan_["stage"], "gene_call_version": plan_["gene_call_version"],
            "plan_sha256": plan_sha256, **{k: done[k] for k in RECORD_FIELDS}}


def record_bytes(doc: dict) -> bytes:
    return json.dumps(doc, sort_keys=True, indent=1).encode() + b"\n"


def record_problem(published: dict, mine: dict) -> str | None:
    """How a batch's record made now differs from its published record (each differing field, both values), or None."""
    diff = sorted(k for k in set(published) | set(mine) if published.get(k) != mine.get(k))
    if not diff:
        return None
    return "; ".join(f"{k}: published {published.get(k)!r}, now {mine.get(k)!r}" for k in diff)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_record(root: Path, out: Path, rel: str, data: bytes, what: str) -> Path:
    """`data` as `out / rel`, under a temporary name renamed into place; refused (FileExistsError) when `root` or `out`
    already holds `rel`: a record is written once."""
    for base in dict.fromkeys([Path(root).resolve(), Path(out).resolve()]):
        if (base / rel).exists():
            raise FileExistsError(f"{base / rel}: {what} is written once; refusing to overwrite")
    dst = Path(out) / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    t = dst.with_name(f".{dst.name}.tmp-{os.getpid()}")
    t.write_bytes(data)
    os.replace(t, dst)
    return dst


def _check_digests(root: Path, stage: str, p: dict, plan_path: Path, tmp) -> None:
    """Every batch's digest of its genome uids is that of its genomes in the stage's `genome_input` (streamed once, in
    the order the batches count them), and the plan counts every genome; else ValueError naming the first batch that
    differs."""
    bat = p["batches"]
    ds = [blake3() for _ in bat]
    j = n = 0
    gen = iter_genome_input(root, stage, tmp=tmp)
    try:
        for g in gen:
            while j < len(bat) and n >= bat[j]["first"] + bat[j]["genomes"]:
                j += 1
            if j < len(bat):
                ds[j].update(_guid(g))
            n += 1
    finally:
        gen.close()
    if n != p["genomes"]:
        raise ValueError(f"{plan_path}: a plan of {p['genomes']:,} genomes, but stage {stage}'s genome_input holds "
                         f"{n:,}: the plan was made from other inputs")
    bad = [k for k, (b, d) in enumerate(zip(bat, ds)) if d.hexdigest(length=16) != b["digest"]]
    if bad:
        raise ValueError(f"{plan_path}: batch {bad[0]}'s digest of its genomes is not that of genomes "
                         f"{bat[bad[0]]['first']} to {bat[bad[0]]['first'] + bat[bad[0]]['genomes']} of stage {stage}'s "
                         "genome_input: the plan was made from other inputs")


def publish_plan(root: Path, stage: str, plan_path: Path, out: Path, *, version: str = GENE_CALL_VERSION,
                 tmp=None) -> dict:
    """The step `genecall.plan` (card B-50): the plan file at `plan_path` (build-stage's, in its work directory),
    checked to be of the stage, the version and the stage's `genome_input` (`_check_digests`), written as it is, byte
    for byte, to `out` as the record `l1/genecall_plan/stage=<S>/plan.json`. Refuses a record that `root` or `out`
    already holds. Writes nothing else; DuckDB spills under `tmp`."""
    root, plan_path = Path(root), Path(plan_path)
    ledger.check_scratch(tmp, root, out)
    raw = plan_path.read_bytes()
    p = read_plan(plan_path)
    if p["stage"] != stage or p["gene_call_version"] != version:
        raise ValueError(f"{plan_path}: a plan of stage {p['stage']} and version {p['gene_call_version']}, not "
                         f"{stage} and {version}")
    _check_digests(root, stage, p, plan_path, tmp)
    _write_record(root, out, plan_record_rel(stage), raw, f"the gene-call plan of stage {stage}")
    return {"record": plan_record_rel(stage), "plan_sha256": _sha256(raw), "batches": len(p["batches"]),
            "genomes": p["genomes"], "estimate_cpu_min": p["estimate_cpu_min"],
            "batch_estimates_cpu_min": [x["estimate_cpu_min"] for x in p["batches"]],
            "merge_estimate_cpu_min": p["merge_estimate_cpu_min"]}


def publish_batch(root: Path, stage: str, plan_path: Path, done: dict, out: Path) -> Path:
    """The record of the batch whose done.json is `done` (`batch_record`), written once to `out`."""
    raw = Path(plan_path).read_bytes()
    rec = batch_record(read_plan(plan_path), _sha256(raw), done)
    return _write_record(root, out, batch_record_rel(stage, done["batch"]), record_bytes(rec),
                         f"the record of gene-call batch {done['batch']} of stage {stage}")


@contextlib.contextmanager
def _stop_on_signals():
    """While a batch runs: SIGTERM and SIGHUP raise `Stopped` (SIGINT raises KeyboardInterrupt, as by default), so its
    cleanup runs; a second signal during the cleanup is ignored. Main thread only (elsewhere, nothing)."""
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


def _usage(t0: float) -> dict:
    """This process's CPU-minutes with its waited-for children (the worker processes), wall seconds since `t0` and
    the larger peak RSS of the two (KiB, `ru_maxrss` on Linux)."""
    me, kids = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return {"cpu_min": round((me.ru_utime + me.ru_stime + kids.ru_utime + kids.ru_stime) / 60, 3),
            "wall_s": round(time.monotonic() - t0, 1), "peak_rss_kib": max(me.ru_maxrss, kids.ru_maxrss)}


def _slice(genomes: Iterable[GenomeInput], first: int, count: int, digest) -> Iterator[GenomeInput]:
    """The genomes `first` to `first + count` of `genomes`, the digest updated with each one's uid; the rest are not
    read."""
    if count == 0:
        return
    for i, g in enumerate(genomes):
        if i < first:
            continue
        digest.update(_guid(g))
        yield g
        if i + 1 >= first + count:
            return


def run_batch(root: Path, stage: str, plan_path: Path, j: int, *, version: str = GENE_CALL_VERSION, workers: int = 1,
              batch_bases: int = BATCH_BASES, check_record: bool = False) -> dict:
    """Batch `j` of the plan at `plan_path` (module docstring): its genomes called as `run` calls them (`_call_into`,
    `workers` at once), their rows, skips and record written into `<work>/.batch-<j>.tmp-<pid>-<hex>`, renamed to
    `<work>/batch-<j>` at the end. Writes nothing into the root. A batch whose output is there is not run again. Before
    any calling, the plan's digest of the batch's genomes is checked against the stage's `genome_input` (ValueError).
    A signal (SIGTERM, SIGHUP) removes the temporary directory and exits 128 + signum. With `check_record` (card B-50:
    a batch made again after its scratch was lost) its output must reproduce the batch's record published in the root
    (`record_problem`), else it is removed and ValueError names the batch and the fields."""
    rec = _run_batch(root, stage, plan_path, j, version=version, workers=workers, batch_bases=batch_bases)
    if check_record:
        published = published_batch(root, stage, j)
        if published is None:
            raise FileNotFoundError(f"{Path(root) / batch_record_rel(stage, j)}: batch {j} of stage {stage} has no "
                                    "published record to check against")
        why = record_problem(published, batch_record(read_plan(plan_path), _sha256(Path(plan_path).read_bytes()), rec))
        if why is not None:
            shutil.rmtree(Path(plan_path).parent / batch_name(j), ignore_errors=True)
            raise ValueError(f"batch {j} of stage {stage}: made again into scratch, it does not reproduce its published "
                             f"record {batch_record_rel(stage, j)} ({why}); its output is removed and the merge does not "
                             "run. A record a batch does not reproduce is not the gene calling of these inputs: tell "
                             "the PI")
        _log(f"batch {j}: reproduces its published record {batch_record_rel(stage, j)}")
    return rec


def _run_batch(root: Path, stage: str, plan_path: Path, j: int, *, version: str, workers: int,
               batch_bases: int) -> dict:
    t0 = time.monotonic()
    root, plan_path = Path(root), Path(plan_path)
    work = plan_path.parent
    p = read_plan(plan_path)
    if p["stage"] != stage or p["gene_call_version"] != version:
        raise ValueError(f"{plan_path}: a plan of stage {p['stage']} and version {p['gene_call_version']}, not "
                         f"{stage} and {version}")
    if not 0 <= j < len(p["batches"]):
        raise ValueError(f"{plan_path}: no batch {j} (the plan has {len(p['batches'])})")
    if workers < 1:
        raise ValueError(f"workers must be at least 1, got {workers}")
    ledger.check_scratch(work, root)
    want = p["batches"][j]
    rec = batch_done(work, p, j)
    if rec is not None:
        _log(f"batch {j}: its output is in {work / batch_name(j)}; not run again")
        return rec
    stages = p["seq_stages"]
    tmpd = work / f".{batch_name(j)}.tmp-{os.getpid()}-{os.urandom(4).hex()}"
    tmpd.mkdir()
    with _stop_on_signals():
        try:
            scratch = tmpd / "scratch"
            scratch.mkdir()
            d, gen = blake3(), iter_genome_input(root, stage, tmp=scratch)
            try:
                for _ in _slice(gen, want["first"], want["genomes"], d):
                    pass
            finally:
                gen.close()
            if d.hexdigest(length=16) != want["digest"]:
                raise ValueError(f"{plan_path}: batch {j}'s digest of its genomes is not that of genomes "
                                 f"{want['first']} to {want['first'] + want['genomes']} of stage {stage}'s "
                                 "genome_input: the plan was made from other inputs")
            _log(f"batch {j}: calling {want['genomes']:,} genomes ({want['bases']:,} bases to call, estimate "
                 f"{want['estimate_cpu_min']:,.2f} CPU-minutes) with {workers} workers, in {tmpd}")
            store = SeqStore(root, stages)
            genomes, shared = stage_genomes(root, stage, None, tmp=scratch)
            try:
                owners = ContigOwners(root, stage, stages, shared)
                counts, skipped = _call_into(store, root, owners, _slice(genomes, want["first"], want["genomes"],
                                                                         blake3()), workers, batch_bases, tmpd)
            finally:
                genomes.close()  # the stream stopped after the batch's last genome: close its DuckDB connection
                store.close()
            _log(f"batch {j}: writing {counts['genes']:,} genes' rows and {len(skipped):,} skips")
            pq.write_table(pa.table({"contig_uid": pa.array(sorted(skipped), pa.binary(16)),
                                     "owner": pa.array([skipped[c] for c in sorted(skipped)], pa.string())},
                                    schema=_SKIPPED_SCHEMA), tmpd / SKIPPED_FILE)
            shutil.rmtree(scratch)
            files = {name: (tmpd / name).stat().st_size for name in (*_SCRATCH_FILES.values(), SKIPPED_FILE)}
            rec = {"batch": j, "first": want["first"], "genomes": want["genomes"], "digest": want["digest"],
                   **counts, "skipped": len(skipped), "files": files, "estimate_cpu_min": want["estimate_cpu_min"],
                   **_usage(t0)}
            (tmpd / DONE_FILE).write_text(json.dumps(rec, sort_keys=True))
            os.rename(tmpd, work / batch_name(j))
        except BaseException:
            shutil.rmtree(tmpd, ignore_errors=True)
            raise
    _log(f"batch {j} done: {rec['genes']:,} genes, {rec['cpu_min']:,.2f} CPU-minutes (estimate "
         f"{rec['estimate_cpu_min']:,.2f}), wall {rec['wall_s']:,.1f} s, peak RSS {rec['peak_rss_kib'] / 1024:,.1f} MiB")
    return rec


def _check_published(root: Path, stage: str, plan_path: Path, p: dict, recs: list[dict]) -> None:
    """When the root holds the stage's published plan (card B-50: build-stage published every piece), the plan file is
    that record, byte for byte, and each batch's output in scratch is its published record (`record_problem`); else
    ValueError, or FileNotFoundError for a batch without a record. A plan never published (a direct call, card B-48's
    commands) is not checked."""
    pub = published_plan(root, stage)
    if pub is None:
        return
    raw = Path(plan_path).read_bytes()
    if pub[0] != raw:
        raise ValueError(f"{plan_path}: not the plan published in {Path(root) / plan_record_rel(stage)}; the batches of "
                         "another plan are not merged")
    for j, r in enumerate(recs):
        published = published_batch(root, stage, j)
        if published is None:
            raise FileNotFoundError(f"{Path(root) / batch_record_rel(stage, j)}: batch {j} of stage {stage} has no "
                                    "published record; run its step first")
        why = record_problem(published, batch_record(p, _sha256(raw), r))
        if why is not None:
            raise ValueError(f"batch {j} of stage {stage}: its output in scratch is not its published record "
                             f"{batch_record_rel(stage, j)} ({why})")


def merge(root: Path, stage: str, plan_path: Path, *, version: str = GENE_CALL_VERSION, out: Path | None = None,
          tmp=None, prefix_width: int | None = None, max_file_bytes: int = layout.MAX_FILE_BYTES) -> dict:
    """The step `genecall` of a stage called in batches: `_publish` once over every batch's rows, so the files, the
    ledger and the result are those of `run` on the stage (module docstring). Refuses, before any work, a stage whose
    gene-call partitions hold a file (as `run`), a plan of another stage or version, and a plan any batch of which has
    no output (FileNotFoundError naming them); and, when the root holds the stage's published plan (card B-50), a plan
    file other than it or a batch whose output is not its published record (`_check_published`)."""
    t0 = time.monotonic()
    root = Path(root)
    out = root if out is None else Path(out)
    ledger.check_scratch(tmp, root, out)
    _refuse_written(root, out, stage, version)
    plan_path = Path(plan_path)
    work = plan_path.parent
    p = read_plan(plan_path)
    if p["stage"] != stage or p["gene_call_version"] != version:
        raise ValueError(f"{plan_path}: a plan of stage {p['stage']} and version {p['gene_call_version']}, not "
                         f"{stage} and {version}")
    recs = [batch_done(work, p, j) for j in range(len(p["batches"]))]
    missing = [j for j, r in enumerate(recs) if r is None]
    if missing:
        raise FileNotFoundError(f"{work}: no output of batch {', batch '.join(map(str, missing))} of "
                                f"{len(p['batches'])}; run them first (bio-box stage genecall batch)")
    _check_published(root, stage, plan_path, p, recs)
    skipped: dict[bytes, str] = {}
    for j in range(len(recs)):
        t = pq.read_table(work / batch_name(j) / SKIPPED_FILE)
        for cu, s in zip(t.column("contig_uid").to_pylist(), t.column("owner").to_pylist()):
            if skipped.setdefault(cu, s) != s:
                raise ValueError(f"contig {cu.hex()}: two batches name two owners, {skipped[cu]} and {s}")
    counts = {"genes": sum(r["genes"] for r in recs), "wrapping_genes": sum(r["wrapping_genes"] for r in recs),
              "modes": sorted({m for r in recs for m in r["modes"]})}
    files = {k: [work / batch_name(j) / name for j in range(len(recs))] for k, name in _SCRATCH_FILES.items()}
    scratch = Path(tempfile.mkdtemp(prefix="biobox-genecall-merge-", dir=tmp))
    try:
        _log(f"merging {len(recs)} batches: {counts['genes']:,} genes, {len(skipped):,} skipped contigs")
        n_prot = _publish(out, stage, version, files, scratch, prefix_width, max_file_bytes, skipped)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    u = _usage(t0)
    _log(f"merge done: {n_prot:,} distinct proteins, {u['cpu_min']:,.2f} CPU-minutes (estimate "
         f"{p['merge_estimate_cpu_min']:,.2f}), wall {u['wall_s']:,.1f} s, peak RSS {u['peak_rss_kib'] / 1024:,.1f} MiB")
    return _result(counts, n_prot, skipped)
