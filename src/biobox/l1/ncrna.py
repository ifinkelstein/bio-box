"""Non-coding RNA features from Rfam covariance models with Infernal `cmsearch` (operation O-SCAN, card E-01).

Every contig of a stage is searched, both strands, with every model of one Rfam CM library, using
Infernal's own filter pipeline (00-spec.md §3.5, D-08), each model's GA cutoff (`--cut_ga`) and the
fixed nucleotide search space `Z_nt = 1e5` Mb (DEC-04), so a stored E-value does not change as the
corpus grows and a contig's hits do not depend on the other contigs searched. Every hit that
`cmsearch --cut_ga` reports (score at or above the model's GA) is kept: neither the spec nor the
decision rows say how overlapping hits of different models (clan members) are resolved, so none is
removed here. A hit's locus is a `feature` row of kind ncrna; two models that hit the same
coordinates and strand share one feature (and its `feature_uid`) and keep one `l1_ncrna` row each.

The producer is `infernal-<infernal version>-rfam-<library version>`. Infernal's version is read from
the binary (`cmsearch -h`, `provenance.tool_versions`); the library's is the `refdb.yaml` version whose
sha256 is the file's (`Rfam.cm.gz` of Rfam 15.1 gives `15.1`), else `sha` + the first 12 hex digits of
the file's sha256 (a subset such as the test fixture's), so neither is written in the code.

Writes, under the output root `out` (the data root by default; `$BIOBOX_STAGING` under run-step, card
C-08), all through `biobox.io.write_parts`:

    l1/feature/producer=<p>/stage=<S>/part-<k>.parquet     kind ncrna, sorted by (contig_uid, start, end, strand)
    l1/ncrna/producer=<p>/stage=<S>/part-<k>.parquet       one row per hit (registry `l1_ncrna`)
    ledger/stage=<S>/step=ncrna-<p>/part-0.parquet          per contig: hits in, kept

Coordinates are 0-based half-open on the contig (cmsearch reports 1-based inclusive, `from > to` on
the minus strand), strand +1/-1. Circular contigs: a contig of length n is searched as one record
holding F bases of its end, the contig, and F bases of its start (repeated when n < F), where F is the
longest hit any model of the library can report (the largest of the models' `W` and their filter
HMMs' `MAXL`). A hit is kept when its leftmost base lies on the contig itself, so each hit, also one
that crosses the origin (then `end > contig_len`), is found once with its full context; a hit in the
flanks is the copy of one kept elsewhere, and truncated hits at the record's ends lie in the flanks.
A hit longer than the circle is impossible on the molecule and is counted in the ledger as rejected
(`exceeds_circle`). On a linear contig cmsearch may report hits truncated by the contig's ends, which
keep their truncation flags.

Bounded memory (AGENTS §4.1): contigs are read one at a time and appended to a FASTA batch file in a
private scratch directory; a batch of at least `batch_bases` bases (or the last one) is searched by one
`cmsearch` call with `cpus` worker threads, and its hits are appended to an unsorted scratch Parquet
file. One DuckDB connection (one thread, spilling into the scratch directory) then sorts the hits into
the output files, so the bytes depend neither on the batches, nor on the order of the genomes, nor on
the thread count. Nothing is published if a batch fails.
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import tempfile
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from biobox import ids, io, layout, ledger, provenance
from biobox.l0seq.inputs import GenomeInput, iter_genome_input
from biobox.l0seq.store import SeqStore
from biobox.registry.schemas import table

__all__ = ["Z_NT_MB", "library_version", "producer", "read_models", "run"]

Z_NT_MB = 1e5  # the nucleotide search space in Mb for every cmsearch of the release line (DEC-04)
KIND_NCRNA = ids.KIND["ncrna"]
LIBRARY = "rfam"  # the refdb.yaml library this producer reads
# bases of contigs searched by one cmsearch call: chosen, not tuned. Each call reads and configures every
# model of the library again, so a batch should be large beside that; the scratch FASTA is about this size
BATCH_BASES = 50_000_000
DUCK_MEMORY_LIMIT = "2GB"  # the sorting connection's cap; larger sorts spill into the scratch directory
# estimated zstd bytes per row; they only cut parts at `max_file_bytes` (card B-19). Not measured
EST_FEATURE_ROW = 24  # estimate: feature_uid 16 + delta-coded coordinates (as genecall's)
EST_NCRNA_ROW = 48  # estimate: feature_uid 16 + coordinates, scores, a dictionary-coded model name
REJECT_CIRCLE = "exceeds_circle"

_HIT_SCRATCH = pa.schema([f for f in table("l1_ncrna").schema])
_CONTIG_SCRATCH = pa.schema([("contig_uid", pa.binary(16)), ("n_in", pa.int64()), ("n_rejected", pa.int64())])


@dataclass(frozen=True)
class Model:
    """One covariance model's header values (the CM file's NAME, ACC, CLEN, W, GA and its filter HMM's MAXL)."""
    name: str
    acc: str
    clen: int
    w: int
    maxl: int
    ga: float | None


def _open_text(path: Path):
    path = Path(path)
    return gzip.open(path, "rt") if path.name.endswith(".gz") else open(path)


def read_models(cm_path: Path) -> dict[str, Model]:
    """The models of a CM file keyed by accession (by name for a model without one), read line by line.
    Each Infernal 1.1 model is a CM header (from `INFERNAL1/a` to the `CM` line) followed by its filter
    HMM (from `HMMER3/f` to the `HMM` line); only those headers are read."""
    out: dict[str, Model] = {}
    head: dict[str, str] = {}
    where = None  # "cm" or "hmm" while in a header
    with _open_text(cm_path) as fh:
        for line in fh:
            if line.startswith("INFERNAL1/"):
                if head:
                    raise ValueError(f"{cm_path}: a CM begins before the filter HMM of {head.get('NAME')}")
                where, head = "cm", {}
                continue
            if line.startswith("HMMER3/"):
                if where is None and not head:
                    raise ValueError(f"{cm_path}: an HMM outside a CM; not an Infernal 1.1 CM file")
                where = "hmm"
                continue
            if where is None:
                continue
            parts = line.split(None, 1)
            if not parts:
                continue
            key = parts[0]
            if where == "cm":
                if key == "CM":
                    where = None
                elif key in ("NAME", "ACC", "CLEN", "W", "GA"):
                    head[key] = parts[1].split()[0]
            else:
                if key == "MAXL":
                    head["MAXL"] = parts[1].split()[0]
                elif key == "HMM":
                    m = Model(name=head["NAME"], acc=head.get("ACC", head["NAME"]), clen=int(head["CLEN"]),
                              w=int(head["W"]), maxl=int(head.get("MAXL", 0)),
                              ga=float(head["GA"]) if "GA" in head else None)
                    if m.acc in out:
                        raise ValueError(f"{cm_path}: model {m.acc} appears twice")
                    out[m.acc] = m
                    where, head = None, {}
    if head or where is not None:
        raise ValueError(f"{cm_path}: the file ends inside the model {head.get('NAME')}")
    if not out:
        raise ValueError(f"{cm_path}: no covariance model")
    return out


def library_version(cm_path: Path) -> str:
    """The Rfam library's version: the `refdb.yaml` version of the entry whose sha256 is the file's, else
    `sha` + the first 12 hex digits of the file's sha256 (a library not frozen in refdb.yaml)."""
    digest = ids.sha256_file(Path(cm_path)).hex()
    for db in provenance.refdb_libraries():
        if db["name"] == LIBRARY and db.get("sha256") == digest:
            return db["version"]
    return "sha" + digest[:12]


def infernal_version() -> str:
    return provenance.tool_versions("infernal")["infernal"]


def producer(cm_path: Path) -> str:
    """`infernal-<version>-rfam-<version>`, both read from the binary and the library file."""
    return f"infernal-{infernal_version()}-{LIBRARY}-{library_version(cm_path)}"


def _paths(stage: str, prod: str) -> tuple[str, str, str]:
    return (f"l1/feature/producer={prod}/stage={stage}", f"{table('l1_ncrna').location}/producer={prod}/stage={stage}",
            f"ledger/stage={stage}/step=ncrna-{prod}")


def _sql_path(p: Path) -> str:
    return "'" + str(p).replace("'", "''") + "'"


def _extended(seq: bytes, flank: int) -> bytes:
    """`flank` bases of the circle's end, the contig, and `flank` bases of its start (the contig repeated
    as often as a flank longer than it needs)."""
    if flank == 0:
        return seq
    ring = seq * -(-flank // len(seq))
    return ring[len(ring) - flank:] + seq + ring[:flank]


class _Batch:
    """One FASTA scratch file of contigs, records named by contig_uid hex, and what each record is."""

    def __init__(self, path: Path):
        self.path = path
        self.fh = open(path, "wb")
        self.contigs: dict[str, tuple[bytes, int, int]] = {}  # name -> (contig_uid, contig_len, flank)
        self.bases = 0

    def add(self, cu: bytes, seq: bytes, flank: int | None) -> None:
        name = cu.hex()
        if name in self.contigs:
            return  # a contig named twice in the input is searched once
        rec = seq if flank is None else _extended(seq, flank)
        self.fh.write(b">" + name.encode() + b"\n" + rec + b"\n")
        self.contigs[name] = (cu, len(seq), 0 if flank is None else flank)
        self.bases += len(rec)

    def close(self) -> None:
        self.fh.close()


def _cmsearch_argv(exe: str, cm: Path, fna: Path, tbl: Path, cpus: int) -> list[str]:
    return [exe, "--cpu", str(cpus), "-Z", f"{Z_NT_MB:g}", "--cut_ga", "--noali", "--tformat", "fasta",
            "--tblout", str(tbl), "-o", os.devnull, str(cm), str(fna)]


def _parse_tblout(path: Path) -> Iterator[dict]:
    """The hits of a `cmsearch --tblout` file (format 1): one dict per data line, fields as printed."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split(None, 17)
            if len(f) < 17:
                raise ValueError(f"{path}: a tblout line with {len(f)} fields: {line!r}")
            yield {"target": f[0], "query": f[2], "query_acc": f[3], "mdl": f[4], "mdl_from": int(f[5]),
                   "mdl_to": int(f[6]), "seq_from": int(f[7]), "seq_to": int(f[8]), "strand": f[9], "trunc": f[10],
                   "gc": float(f[12]), "score": float(f[14]), "evalue": float(f[15]), "inc": f[16]}


def _search(exe: str, cm: Path, models: dict[str, Model], batch: _Batch, cpus: int, scratch: Path
            ) -> tuple[pa.RecordBatch, pa.RecordBatch]:
    """Search one batch file; its hits on the contigs as rows of `l1_ncrna`, and per contig the hits in and
    rejected."""
    tbl = scratch / "batch.tbl"
    tbl.unlink(missing_ok=True)
    done = subprocess.run(_cmsearch_argv(exe, cm, batch.path, tbl, cpus), stdin=subprocess.DEVNULL,
                          stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if done.returncode != 0:
        raise RuntimeError(f"cmsearch exited with code {done.returncode}: {done.stderr.decode(errors='replace')[-2000:]}")
    rows: list[dict] = []
    n_in: dict[str, int] = dict.fromkeys(batch.contigs, 0)
    n_rej: dict[str, int] = dict.fromkeys(batch.contigs, 0)
    for h in _parse_tblout(tbl):
        cu, n, flank = batch.contigs[h["target"]]
        left, right = min(h["seq_from"], h["seq_to"]) - 1, max(h["seq_from"], h["seq_to"])
        start, end = left - flank, right - flank
        if flank and not 0 <= start < n:
            continue  # the copy, in a flank, of a hit kept where its leftmost base lies on the contig
        n_in[h["target"]] += 1
        if flank and end - start > n:
            n_rej[h["target"]] += 1
            continue
        m = models.get(h["query_acc"] if h["query_acc"] != "-" else h["query"])
        if m is None:
            raise ValueError(f"cmsearch reported model {h['query_acc']} ({h['query']}), which {cm} does not hold")
        if h["strand"] not in ("+", "-"):
            raise ValueError(f"cmsearch reported strand {h['strand']!r}")
        strand = 1 if h["strand"] == "+" else -1
        rows.append({
            "feature_uid": ids.feature_uid(cu, start, end, strand, KIND_NCRNA), "contig_uid": cu,
            "start": start, "end": end, "strand": strand, "model_acc": m.acc, "model_name": h["query"],
            "model_len": m.clen, "mdl_start": h["mdl_from"] - 1, "mdl_end": h["mdl_to"],
            "mdl_kind": h["mdl"], "bitscore": h["score"], "evalue_zref": h["evalue"], "ga_bits": m.ga,
            "z_mb": Z_NT_MB, "trunc_5p": "5'" in h["trunc"], "trunc_3p": "3'" in h["trunc"], "gc": h["gc"],
            "contig_len": n,
        })
    names = sorted(batch.contigs)
    contigs = pa.RecordBatch.from_arrays(
        [pa.array([batch.contigs[k][0] for k in names], type=pa.binary(16)),
         pa.array([n_in[k] for k in names], type=pa.int64()), pa.array([n_rej[k] for k in names], type=pa.int64())],
        schema=_CONTIG_SCRATCH)
    return pa.RecordBatch.from_pylist(rows, schema=_HIT_SCRATCH), contigs


def _contigs(store: SeqStore, genomes: Iterable[GenomeInput]) -> Iterator[tuple[bytes, bytes, bool]]:
    """(contig_uid, upper-case sequence, circular) of every contig of the genomes, one at a time, in input order."""
    for g in genomes:
        circ = set(g.circular)
        for cu in g.contig_uids:
            yield cu, store.get(cu).upper(), cu in circ


def run(root: Path, stage: str, cm_path: Path, cpus: int, genomes: Iterable[GenomeInput] | None = None, *,
        tmp=None, out: Path | None = None, seq_stages: list[str] | None = None, batch_bases: int = BATCH_BASES,
        max_file_bytes: int = layout.MAX_FILE_BYTES) -> dict:
    """Search the stage's contigs with the CM library `cm_path` and write its ncRNA features and hits.

    `cpus`: cmsearch worker threads (`--cpu`), given explicitly; the output bytes do not depend on it.
    `genomes`: any iterable, consumed once; None reads the stage's `genome_input` table (which also says
    which contigs are circular). `seq_stages`: the stages whose sequences the genomes may name (card
    B-28), default the stage itself. `tmp` holds the private scratch directory (system temp by default)
    and may lie under neither the root nor `out`. `batch_bases`: bases per cmsearch call; the output
    does not depend on it. `max_file_bytes`: each table is cut into parts of at most this many
    estimated bytes.

    Inputs are read from `root`; every output and the ledger are written under `out` (default `root`;
    card C-08). Refuses, before reading a sequence, when one of the producer's partitions of this stage
    already holds a file in `root` or `out`.
    """
    root = Path(root)
    out = root if out is None else Path(out)
    if not isinstance(cpus, int) or isinstance(cpus, bool) or cpus < 1:
        raise ValueError(f"cpus {cpus!r}: give cmsearch's thread count explicitly (at least 1)")
    if batch_bases < 1:
        raise ValueError(f"batch_bases must be at least 1, got {batch_bases}")
    ledger.check_scratch(tmp, root, out)
    exe = shutil.which("cmsearch")
    if exe is None:
        raise FileNotFoundError("cmsearch (Infernal) not on PATH")
    cm_path = Path(cm_path)
    prod = producer(cm_path)
    feat_dir, hit_dir, step_dir = _paths(stage, prod)
    ledger.refuse_written(root, out, [feat_dir, hit_dir, step_dir], f"ncRNA features {prod} of stage {stage}")
    models = read_models(cm_path)
    no_ga = sorted(a for a, m in models.items() if m.ga is None)
    if no_ga:
        raise ValueError(f"{cm_path}: {len(no_ga)} model(s) without a GA cutoff, first {no_ga[0]}; --cut_ga needs one")
    flank = max(max(m.w, m.maxl) for m in models.values())
    tools = {**provenance.tool_versions("infernal", "pyarrow"), LIBRARY: library_version(cm_path)}
    store = SeqStore(root, [stage] if seq_stages is None else list(seq_stages))
    scratch = Path(tempfile.mkdtemp(prefix="biobox-ncrna-", dir=tmp))
    try:
        cm = cm_path
        if cm_path.name.endswith(".gz"):  # cmsearch reads a plain CM file
            cm = scratch / "library.cm"
            with gzip.open(cm_path, "rb") as src, open(cm, "wb") as dst:
                shutil.copyfileobj(src, dst, 1 << 20)
        if genomes is None:
            genomes = iter_genome_input(root, stage, tmp=scratch)
        hits_w = pq.ParquetWriter(scratch / "hits.parquet", _HIT_SCRATCH, compression="zstd", write_statistics=False)
        contigs_w = pq.ParquetWriter(scratch / "contigs.parquet", _CONTIG_SCRATCH, compression="zstd",
                                     write_statistics=False)
        stats = {"contigs": 0, "bases": 0, "batches": 0, "hits": 0}
        try:
            batch = None

            def flush(b: _Batch) -> None:
                b.close()
                hits, contigs = _search(exe, cm, models, b, cpus, scratch)
                hits_w.write_batch(hits)
                contigs_w.write_batch(contigs)
                stats["batches"] += 1
                stats["hits"] += hits.num_rows
                b.path.unlink()

            for cu, seq, circ in _contigs(store, genomes):
                if batch is None:
                    batch = _Batch(scratch / "batch.fna")
                batch.add(cu, seq, flank if circ else None)
                stats["contigs"] += 1
                stats["bases"] += len(seq)
                if batch.bases >= batch_bases:
                    flush(batch)
                    batch = None
            if batch is not None:
                flush(batch)
        finally:
            hits_w.close()
            contigs_w.close()
        _publish(out, stage, prod, scratch, tools, max_file_bytes)
        return {**stats, "producer": prod, "models": len(models), "flank": flank}
    finally:
        store.close()
        shutil.rmtree(scratch, ignore_errors=True)


def _publish(out: Path, stage: str, prod: str, scratch: Path, tools: dict[str, str], max_file_bytes: int) -> None:
    """Sort the scratch hits into the producer's `feature` and `l1_ncrna` files and write the ledger."""
    feat_dir, hit_dir, _ = _paths(stage, prod)
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb"), "preserve_insertion_order": True})
    try:
        con.execute(f"CREATE VIEW hits AS SELECT DISTINCT * FROM read_parquet({_sql_path(scratch / 'hits.parquet')})")
        con.execute("CREATE VIEW contigs AS SELECT contig_uid, sum(n_in) AS n_in, sum(n_rejected) AS n_rejected "
                    f"FROM read_parquet({_sql_path(scratch / 'contigs.parquet')}) GROUP BY contig_uid")
        spec = table("feature")
        io.write_parts(io.duck_batches(
            con, f'SELECT DISTINCT feature_uid, contig_uid, start, "end", strand, {KIND_NCRNA}::UTINYINT AS kind '
                 'FROM hits ORDER BY contig_uid, start, "end", strand', spec),
            spec, out / feat_dir, EST_FEATURE_ROW, max_file_bytes)
        spec = table("l1_ncrna")
        cols = ", ".join(f'"{c}"' for c in spec.schema.names)
        io.write_parts(io.duck_batches(
            con, f'SELECT {cols} FROM hits ORDER BY contig_uid, start, "end", strand, model_acc', spec),
            spec, out / hit_dir, EST_NCRNA_ROW, max_file_bytes)

        def ledger_rows() -> Iterator[dict]:
            """Per contig (every contig searched, also one without a hit), in contig_uid order: hits in, kept,
            and the hits longer than a circular contig rejected."""
            reader = con.sql("SELECT contig_uid, n_in, n_rejected FROM contigs ORDER BY contig_uid"
                             ).to_arrow_reader(io.ROW_GROUP_ROWS)
            try:
                for b in reader:
                    for cu, n, rej in zip(*(b.column(i).to_pylist() for i in range(3))):
                        if rej:
                            yield ledger.row(cu, "contig", n, n - rej, n_rejected=rej, reason=REJECT_CIRCLE,
                                             tools=tools)
                        yield ledger.row(cu, "contig", n, n - rej, tools=tools)
            finally:
                reader.close()

        ledger.write(out, stage, f"ncrna-{prod}", ledger_rows())
    finally:
        con.close()
