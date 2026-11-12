"""Non-coding RNA features from Rfam covariance models with Infernal `cmsearch` (operation O-SCAN, cards E-01, E-07).

What is scanned (DEC-28, the specification `planning/design/0-08_ncrna_models.md` §2, §4, §6). A model set file
(`rfam_sets/<name>.tsv`: accession, name, and `whole_genome` 0 or 1; `prok-1` by default, 588 models of Rfam 15.1)
names the models; the producer fetches them from the Rfam library (`cmfetch -f`) into scratch, as two CM files:

- the set's intergenic models (`whole_genome` 0) search the contig's intergenic regions with a margin of `MARGIN`
  (200) bases into the flanking genes, with `cmsearch --notrunc` (no re-search of record ends for truncated hits);
- the set's whole-genome models (`whole_genome` 1: RNA inside coding sequence, 23 in prok-1) search whole contigs
  with Infernal's default truncated-hit handling, as card E-01 searched every model.

So no model searches a base twice: the two model sets are disjoint and the regions of a contig are disjoint.
Both use Infernal's own filter pipeline, each model's GA cutoff (`--cut_ga`) and the fixed nucleotide search space
`Z_nt = 1e5` Mb (DEC-04), so a stored E-value does not change as the corpus grows and a hit's GA decision does not
depend on the other models or contigs searched. Every hit at or above GA is kept (clan competition is card E-08).
A hit's locus is a `feature` row of kind ncrna; two models that hit the same coordinates and strand share one
feature (and its `feature_uid`) and keep one `l1_ncrna` row each.

Regions (`regions`): the intergenic bases of a contig are those no gene call of the stage's gene-call version
covers (either strand; the runs `biobox.l1.igr.igrs_of_contig` gives); a region is a maximal run of bases that lie
within `MARGIN` bases of an intergenic base: each intergenic run widened by the margin on both sides, clipped at
the ends of a linear contig, taken around the origin of a circular one, overlapping or touching runs merged (also
across the origin). A region that crosses the origin of a circular contig is one region `[s, e)` with `e > n`,
cut from the doubled contig; a region that covers the whole circle is the whole contig (a contig without genes)
and is searched as a circular whole contig is (below). A contig whose genes cover it has no region. Each region
is one record named `<contig_uid hex>_<s>`, cut from the + strand (cmsearch searches both), and a hit at record
offsets `[a, b)` lies at `[s + a, s + b)` on the contig, less n when that start is past the origin, so a hit
crossing the origin has `end > contig_len` (AGENTS §3.1).

Whole contigs: a linear contig is one record. A circular contig of length n is searched as one record holding
F bases of its end, the contig, and F bases of its start (repeated when n < F), where F is the longest hit a
model of that CM can report (the largest of the models' `W` and their filter HMMs' `MAXL`); a hit is kept when its
leftmost base lies on the contig itself, so each hit, also one across the origin, is found once with its full
context; a hit longer than the circle is impossible on the molecule and is counted in the ledger as rejected
(`exceeds_circle`). On a linear contig cmsearch may report hits truncated by the contig's ends, which keep their
truncation flags.

The producer is `infernal-<infernal version>-rfam-<library version>-<set version>-<REGION_RULE>`. Infernal's version
is read from the binary (`cmsearch -h`); the library's is the `refdb.yaml` version whose sha256 is the file's
(`Rfam.cm.gz` of Rfam 15.1 gives `15.1`), else `sha` + the first 12 hex digits of its sha256; the set's is its name
without hyphens when the file is a frozen set of `FROZEN_SETS` (`prok-1` gives `prok1`), else `sha` + 12 hex
digits of the file's sha256; `REGION_RULE` (`igr200nt`) names the margin and `--notrunc`. A change of any of them
is a new producer, so a new tape version.

Writes, under the output root `out` (the data root by default; `$BIOBOX_STAGING` under run-step, card C-08), all
through `biobox.io.write_parts`:

    l1/feature/producer=<p>/stage=<S>/part-<k>.parquet     kind ncrna, sorted by (contig_uid, start, end, strand)
    l1/ncrna/producer=<p>/stage=<S>/part-<k>.parquet       one row per hit (registry `l1_ncrna`)
    ledger/stage=<S>/step=ncrna-<p>/part-0.parquet          per contig: hits in, kept

A contig that another stage stored first (DEC-26, card C-20) is that stage's: it is not searched here and has one
ledger row, reason `owned_by:<stage>`, counts 0.

Two ways to scan a stage, with byte-identical output (card E-07, as card B-24 shards the domain scan):

- `run`: every contig in one process.
- In shards, each its own process under the 300 CPU-minute limit (AGENTS §7.2): `plan` cuts the stage's contigs,
  sorted by `contig_uid` (its `genome_input` rows), into shards of at most `scan_bases` genome bases (each contig's
  `contig_len` in the stage's `contig_index`; a contig another stage owns counts 0; a contig longer than the limit
  is a shard of its own) and writes the plan, `l1/domain_scan/library=<p>/stage=<S>/plan.json` (the `domain_scan`
  artifact of the registry, its library partition naming this producer); `scan_shard` scans shard k into two
  sorted scratch runs, `<tmp>/biobox-ncrna-scan/<plan>/shard-<k>.hits.parquet` and `.contigs.parquet` (`<plan>`: the
  first 16 hex of the plan's sha256), and writes its record `.../shard-<k>.json` beside the plan (the runs' sha256
  and counts); `merge` checks every run against its record and writes the runs as `run` writes its files (one
  DuckDB sort of the same rows by the same statement). A run lost from scratch after its record was published is
  scanned again by `rescan_shard` (`python -m biobox.l1.ncrna rescan`, as card C-12 does for domains), which must
  reproduce the record's sha256. Shards are a pure function of the sorted contigs and the limit.

Bounded memory (AGENTS §4.1): contigs and the gene calls are streamed in `contig_uid` order (one contig's genes at a
time) and each contig's records are appended to a FASTA batch file in a private scratch directory; a batch of at
least `batch_bases` bases (or the last one) is searched by one `cmsearch` call with `cpus` worker threads, and its
hits are appended to an unsorted scratch Parquet file. One DuckDB connection (one thread, spilling into the scratch
directory) then sorts the hits into the output files, so the bytes depend neither on the batches, nor on the shards,
nor on the thread count. Nothing is published if a batch fails.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from biobox import ids, io, layout, ledger, provenance
from biobox.l0seq.inputs import GenomeInput
from biobox.l0seq.store import INDEX_DIR, ContigOwners, SeqStore
from biobox.l0seq.store import OWNED_BY as SKIP_REASON
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.l1.igr import igrs_of_contig
from biobox.registry.schemas import ARTIFACTS, table

__all__ = ["DEFAULT_SET", "MARGIN", "REGION_RULE", "Z_NT_MB", "library_version", "merge", "plan", "producer",
           "read_models", "read_set", "regions", "rescan_shard", "run", "scan_shard", "set_version"]

Z_NT_MB = 1e5  # the nucleotide search space in Mb for every cmsearch of the release line (DEC-04)
KIND_NCRNA = ids.KIND["ncrna"]
LIBRARY = "rfam"  # the refdb.yaml library this producer reads
MARGIN = 200  # bases of each flanking gene added to an intergenic region (DEC-28)
REGION_RULE = f"igr{MARGIN}nt"  # intergenic regions at this margin, searched with --notrunc
SET_DIR = Path(__file__).resolve().parent / "rfam_sets"
DEFAULT_SET = "prok-1"
# Frozen model sets: name -> sha256 of the file. A frozen file is never regenerated; a new set is a new file
# with a new name (DEC-28). prok-1: planning/design/0-08_ncrna_models.md, 588 models, 23 of them whole-genome.
FROZEN_SETS = {"prok-1": "092660962def1ec76466de3816558451f478637ffb43d141ef684273f68d6179"}
# bases of records searched by one cmsearch call: chosen, not tuned. Each call reads and configures every
# model of its CM file again, so a batch should be large beside that; the scratch FASTA is about this size
BATCH_BASES = 50_000_000
DUCK_MEMORY_LIMIT = "2GB"  # the sorting connection's cap; larger sorts spill into the scratch directory
# estimated zstd bytes per row; they only cut parts at `max_file_bytes` (card B-19). Not measured
EST_FEATURE_ROW = 24  # estimate: feature_uid 16 + delta-coded coordinates (as genecall's)
EST_NCRNA_ROW = 48  # estimate: feature_uid 16 + coordinates, scores, a dictionary-coded model name
REJECT_CIRCLE = "exceeds_circle"
PLAN_FORMAT = 1
SCAN_DIR = "biobox-ncrna-scan"  # under the scratch directory: the shards' runs, one directory per plan
_READ_ROWS = 65_536  # rows per read of genome_input, contig_index and the gene calls

_HIT_SCRATCH = pa.schema([f for f in table("l1_ncrna").schema])
_CONTIG_SCRATCH = pa.schema([("contig_uid", pa.binary(16)), ("n_in", pa.int64()), ("n_rejected", pa.int64()),
                             ("reason", pa.string())])


@dataclass(frozen=True)
class Model:
    """One covariance model's header values (the CM file's NAME, ACC, CLEN, W, GA and its filter HMM's MAXL)."""
    name: str
    acc: str
    clen: int
    w: int
    maxl: int
    ga: float | None


@dataclass(frozen=True)
class SetModel:
    """One row of a model set file: the model and whether it searches whole contigs (else intergenic regions)."""
    acc: str
    name: str
    whole_genome: bool


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
    digest = _sha256(Path(cm_path))
    for db in provenance.refdb_libraries():
        if db["name"] == LIBRARY and db.get("sha256") == digest:
            return db["version"]
    return "sha" + digest[:12]


@lru_cache(maxsize=16)
def _sha256_cached(path: str, size: int, mtime_ns: int) -> str:
    return ids.sha256_file(Path(path)).hex()


def _sha256(path: Path) -> str:
    """sha256 hex of a file, computed once per (path, size, mtime) in a process (a library is read many times)."""
    st = Path(path).stat()
    return _sha256_cached(str(Path(path).resolve()), st.st_size, st.st_mtime_ns)


def infernal_version() -> str:
    return provenance.tool_versions("infernal")["infernal"]


# --- the model set ---------------------------------------------------------------------------------------------


def set_path(model_set: str | Path) -> Path:
    """A model set file: `model_set` as a path when such a file exists, else `rfam_sets/<model_set>.tsv`."""
    p = Path(model_set)
    if p.is_file():
        return p
    q = SET_DIR / f"{model_set}.tsv"
    if q.is_file():
        return q
    raise FileNotFoundError(f"model set {model_set!r}: neither a file nor one of {sorted(x.stem for x in SET_DIR.glob('*.tsv'))}")


def set_version(model_set: str | Path) -> str:
    """The set's version in the producer name: a frozen set's name without hyphens (`prok-1` gives `prok1`), else
    `sha` + the first 12 hex digits of the file's sha256."""
    digest = _sha256(set_path(model_set))
    for name, sha in sorted(FROZEN_SETS.items()):
        if sha == digest:
            return name.replace("-", "")
    return "sha" + digest[:12]


def read_set(model_set: str | Path) -> list[SetModel]:
    """The models of a set file, in file order: `#` lines are comments, then the header `acc name whole_genome`
    (tab-separated) and one row per model. A file of `rfam_sets/` named as a frozen set must have its sha256."""
    path = set_path(model_set)
    frozen = FROZEN_SETS.get(path.stem)
    if frozen is not None and path.resolve().parent == SET_DIR and _sha256(path) != frozen:
        raise ValueError(f"{path}: the frozen model set {path.stem} has another sha256 than {frozen}; a frozen set is "
                         "never changed (a new set is a new file with a new name)")
    rows = [ln.rstrip("\n").split("\t") for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    if not rows or rows[0] != ["acc", "name", "whole_genome"]:
        raise ValueError(f"{path}: the first line that is not a comment must be the header acc, name, whole_genome")
    out, seen = [], set()
    for k, r in enumerate(rows[1:], 2):
        if len(r) != 3 or r[2] not in ("0", "1") or not r[0]:
            raise ValueError(f"{path}: row {k} is not <acc> <name> <0|1>: {r!r}")
        if r[0] in seen:
            raise ValueError(f"{path}: model {r[0]} is listed twice")
        seen.add(r[0])
        out.append(SetModel(r[0], r[1], r[2] == "1"))
    if not out:
        raise ValueError(f"{path}: no model")
    return out


def producer(cm_path: Path, model_set: str | Path = DEFAULT_SET) -> str:
    """`infernal-<version>-rfam-<version>-<set version>-<region rule>`, read from the binary, the library file and
    the set file (module docstring)."""
    return f"infernal-{infernal_version()}-{LIBRARY}-{library_version(cm_path)}-{set_version(model_set)}-{REGION_RULE}"


def _paths(stage: str, prod: str) -> tuple[str, str, str]:
    return (f"l1/feature/producer={prod}/stage={stage}", f"{table('l1_ncrna').location}/producer={prod}/stage={stage}",
            f"ledger/stage={stage}/step=ncrna-{prod}")


def _sql_path(p: Path) -> str:
    return "'" + str(p).replace("'", "''") + "'"


def _sql_files(files: Iterable[Path]) -> str:
    return "[" + ", ".join(_sql_path(f) for f in files) + "]"


# --- regions ---------------------------------------------------------------------------------------------------


def regions(genes: Iterable[tuple[int, int]], n: int, circular: bool, margin: int = MARGIN) -> list[tuple[int, int]]:
    """The intergenic regions of one contig of length n with gene calls `genes` ((start, end), either strand; on a
    circular contig a gene across the origin has end > n), widened by `margin` (module docstring): sorted
    (start, end), 0 <= start < n; on a circular contig a region across the origin has end > n, and one that covers
    the whole circle is (0, n). Empty when the genes cover every base."""
    if n < 1:
        return []
    runs = [(s, e) for s, e, _, _ in igrs_of_contig([(i, s, e) for i, (s, e) in enumerate(genes)], n, circular)]
    if not runs:
        return []
    if not circular:
        w = sorted((max(0, s - margin), min(n, e + margin)) for s, e in runs)
    else:
        w = []
        for s, e in runs:
            s, e = s - margin, e + margin
            if e - s >= n:
                return [(0, n)]
            w.append((s % n, s % n + (e - s)))
        w.sort()
    out = [list(w[0])]
    for s, e in w[1:]:
        if s <= out[-1][1]:  # overlapping or touching
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    if circular:
        # only the last run can pass the origin; it may reach into the first runs (their coordinates + n)
        while len(out) > 1 and out[-1][1] >= out[0][0] + n:
            first = out.pop(0)
            out[-1][1] = max(out[-1][1], first[1] + n)
        if out[-1][1] - out[-1][0] >= n:
            return [(0, n)]
    return [(s, e) for s, e in out]


def _extended(seq: bytes, flank: int) -> bytes:
    """`flank` bases of the circle's end, the contig, and `flank` bases of its start (the contig repeated
    as often as a flank longer than it needs)."""
    if flank == 0:
        return seq
    ring = seq * -(-flank // len(seq))
    return ring[len(ring) - flank:] + seq + ring[:flank]


# --- the scan --------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class _Rec:
    """What one FASTA record is: its contig, the contig's length, the contig offset of the record's first base
    (`base`), the flank of a circular whole-contig record (leftmost-base rule) and whether the contig is circular."""
    contig: bytes
    n: int
    base: int
    flank: int
    circular: bool


class _Batch:
    """One FASTA scratch file of records for one CM file, and what each record is."""

    def __init__(self, path: Path):
        self.path = path
        self.fh = open(path, "wb")
        self.recs: dict[str, _Rec] = {}
        self.bases = 0

    def add(self, name: str, rec: bytes, what: _Rec) -> None:
        if name in self.recs:
            return  # a contig named twice in the input is searched once
        self.fh.write(b">" + name.encode() + b"\n" + rec + b"\n")
        self.recs[name] = what
        self.bases += len(rec)

    def close(self) -> None:
        self.fh.close()


def _cmsearch_argv(exe: str, cm: Path, fna: Path, tbl: Path, cpus: int, notrunc: bool = False) -> list[str]:
    return [exe, "--cpu", str(cpus), "-Z", f"{Z_NT_MB:g}", "--cut_ga", *(["--notrunc"] if notrunc else []), "--noali",
            "--tformat", "fasta", "--tblout", str(tbl), "-o", os.devnull, str(cm), str(fna)]


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


def _search(exe: str, cm: Path, models: dict[str, Model], batch: _Batch, cpus: int, scratch: Path, notrunc: bool
            ) -> tuple[pa.RecordBatch, pa.RecordBatch]:
    """Search one batch file; its hits on the contigs as rows of `l1_ncrna`, and per contig the hits in and
    rejected."""
    tbl = scratch / "batch.tbl"
    tbl.unlink(missing_ok=True)
    done = subprocess.run(_cmsearch_argv(exe, cm, batch.path, tbl, cpus, notrunc), stdin=subprocess.DEVNULL,
                          stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if done.returncode != 0:
        raise RuntimeError(f"cmsearch exited with code {done.returncode}: {done.stderr.decode(errors='replace')[-2000:]}")
    rows: list[dict] = []
    n_in: dict[bytes, int] = {}
    n_rej: dict[bytes, int] = {}
    for h in _parse_tblout(tbl):
        r = batch.recs[h["target"]]
        n = r.n
        left, right = min(h["seq_from"], h["seq_to"]) - 1, max(h["seq_from"], h["seq_to"])
        if r.flank:
            start, end = left - r.flank, right - r.flank
            if not 0 <= start < n:
                continue  # the copy, in a flank, of a hit kept where its leftmost base lies on the contig
        else:
            start, end = r.base + left, r.base + right
            if r.circular and start >= n:  # a region across the origin: past it, the start is taken modulo n
                start, end = start - n, end - n
        n_in[r.contig] = n_in.get(r.contig, 0) + 1
        if r.flank and end - start > n:
            n_rej[r.contig] = n_rej.get(r.contig, 0) + 1
            continue
        m = models.get(h["query_acc"] if h["query_acc"] != "-" else h["query"])
        if m is None:
            raise ValueError(f"cmsearch reported model {h['query_acc']} ({h['query']}), which {cm} does not hold")
        if h["strand"] not in ("+", "-"):
            raise ValueError(f"cmsearch reported strand {h['strand']!r}")
        strand = 1 if h["strand"] == "+" else -1
        rows.append({
            "feature_uid": ids.feature_uid(r.contig, start, end, strand, KIND_NCRNA), "contig_uid": r.contig,
            "start": start, "end": end, "strand": strand, "model_acc": m.acc, "model_name": h["query"],
            "model_len": m.clen, "mdl_start": h["mdl_from"] - 1, "mdl_end": h["mdl_to"],
            "mdl_kind": h["mdl"], "bitscore": h["score"], "evalue_zref": h["evalue"], "ga_bits": m.ga,
            "z_mb": Z_NT_MB, "trunc_5p": "5'" in h["trunc"], "trunc_3p": "3'" in h["trunc"], "gc": h["gc"],
            "contig_len": n,
        })
    keys = sorted(n_in)
    contigs = pa.RecordBatch.from_arrays(
        [pa.array(keys, type=pa.binary(16)), pa.array([n_in[k] for k in keys], type=pa.int64()),
         pa.array([n_rej.get(k, 0) for k in keys], type=pa.int64()), pa.array([None] * len(keys), type=pa.string())],
        schema=_CONTIG_SCRATCH)
    return pa.RecordBatch.from_pylist(rows, schema=_HIT_SCRATCH), contigs


@dataclass
class _Library:
    """The set's models fetched from the library into scratch: one CM file per kind (None when the set has no
    model of that kind), their header values, and the flank of each kind's circular whole-contig records."""
    whole_cm: Path | None
    igr_cm: Path | None
    whole: dict[str, Model]
    igr: dict[str, Model]

    @property
    def flank_whole(self) -> int:
        return max((max(m.w, m.maxl) for m in self.whole.values()), default=0)

    @property
    def flank_igr(self) -> int:
        return max((max(m.w, m.maxl) for m in self.igr.values()), default=0)


def _check_set(cm_path: Path, model_set) -> tuple[list[SetModel], dict[str, Model]]:
    """The set's models and the library's; ValueError when a set model is not in the library or has no GA."""
    members = read_set(model_set)
    lib = read_models(cm_path)
    missing = [m.acc for m in members if m.acc not in lib]
    if missing:
        raise ValueError(f"{cm_path}: {len(missing)} model(s) of the set {set_path(model_set).name} are not in the "
                         f"library, first {missing[0]}")
    no_ga = [m.acc for m in members if lib[m.acc].ga is None]
    if no_ga:
        raise ValueError(f"{cm_path}: {len(no_ga)} model(s) of the set without a GA cutoff, first {no_ga[0]}; --cut_ga "
                         "needs one")
    return members, lib


def _fetch(cm_path: Path, model_set, scratch: Path) -> _Library:
    """`cmfetch -f` the set's models of each kind from the library (decompressed into scratch when gzipped) into
    `scratch/whole.cm` and `scratch/igr.cm`, and check that each holds exactly its models."""
    members, _ = _check_set(cm_path, model_set)
    exe = shutil.which("cmfetch")
    if exe is None:
        raise FileNotFoundError("cmfetch (Infernal) not on PATH")
    plain = Path(cm_path)
    if plain.name.endswith(".gz"):  # cmfetch and cmsearch read a plain CM file
        plain = scratch / "library.cm"
        with gzip.open(cm_path, "rb") as src, open(plain, "wb") as dst:
            shutil.copyfileobj(src, dst, 1 << 20)
    got: dict[str, tuple[Path | None, dict[str, Model]]] = {}
    for kind, want in (("whole", sorted(m.acc for m in members if m.whole_genome)),
                       ("igr", sorted(m.acc for m in members if not m.whole_genome))):
        if not want:
            got[kind] = (None, {})
            continue
        keys, cm = scratch / f"{kind}.keys", scratch / f"{kind}.cm"
        keys.write_text("".join(a + "\n" for a in want))
        done = subprocess.run([exe, "-f", "-o", str(cm), str(plain), str(keys)], stdin=subprocess.DEVNULL,
                              stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        if done.returncode != 0:
            raise RuntimeError(f"cmfetch exited with code {done.returncode}: {done.stderr.decode(errors='replace')[-2000:]}")
        models = read_models(cm)
        if sorted(models) != want:
            raise RuntimeError(f"cmfetch gave {len(models)} model(s) for the {len(want)} {kind} models of the set")
        got[kind] = (cm, models)
    if plain != Path(cm_path):
        plain.unlink()
    return _Library(got["whole"][0], got["igr"][0], got["whole"][1], got["igr"][1])


def _range_rows(files: list[Path], cols: list[str], lo: bytes | None = None, hi: bytes | None = None
                ) -> Iterator[tuple]:
    """Rows (`cols`, the first being `contig_uid`) of files sorted by contig_uid, in order, with lo <= contig_uid
    <= hi when given; reads only the row groups whose contig_uid range meets [lo, hi]. A row that sorts before the
    one above it is a RuntimeError."""
    last = None
    for f in files:
        pf = pq.ParquetFile(f)
        col = pf.schema_arrow.get_field_index(cols[0])
        for g in range(pf.metadata.num_row_groups):
            st = pf.metadata.row_group(g).column(col).statistics
            if st is not None and st.has_min_max and ((lo is not None and bytes(st.max) < lo)
                                                      or (hi is not None and bytes(st.min) > hi)):
                continue
            for b in pf.iter_batches(batch_size=_READ_ROWS, row_groups=[g], columns=cols):
                for row in zip(*(b.column(i).to_pylist() for i in range(len(cols)))):
                    cu = row[0]
                    if (lo is not None and cu < lo) or (hi is not None and cu > hi):
                        continue
                    if last is not None and cu < last:
                        raise RuntimeError(f"{f}: rows are not in contig_uid order")
                    last = cu
                    yield row


def _with_genes(contigs: Iterable[tuple[bytes, bool]], feats: Iterator[tuple] | None
                ) -> Iterator[tuple[bytes, bool, list[tuple[int, int]]]]:
    """(contig_uid, circular, its genes) for each of `contigs` (sorted by contig_uid), the genes from the sorted
    gene-call rows `feats` (contig_uid, start, end): a merge join, one contig's genes held at a time."""
    cur = None if feats is None else next(feats, None)
    for cu, circ in contigs:
        while cur is not None and cur[0] < cu:
            cur = next(feats, None)
        genes = []
        while cur is not None and cur[0] == cu:
            genes.append((cur[1], cur[2]))
            cur = next(feats, None)
        yield cu, circ, genes


def _gene_files(root: Path, stage: str, gene_call_version: str) -> list[Path]:
    return io.part_files(Path(root) / f"l1/feature/producer={gene_call_version}/stage={stage}")


def _genome_input_files(root: Path, stage: str) -> list[Path]:
    return io.part_files(Path(root) / table("genome_input").location / f"stage={stage}")


class _Scan:
    """The scan of a run of contigs into unsorted scratch Parquet files: the hits (`l1_ncrna` rows) and, per contig,
    hits in and rejected (a row of zeros for each contig searched, the counts per batch, the reason for a contig
    another stage owns)."""

    def __init__(self, exe: str, lib: _Library, cpus: int, scratch: Path, batch_bases: int, margin: int):
        self.exe, self.lib, self.cpus, self.scratch = exe, lib, cpus, scratch
        self.batch_bases, self.margin = batch_bases, margin
        self.hits_path, self.contigs_path = scratch / "hits.parquet", scratch / "contigs.parquet"
        self.hits_w = pq.ParquetWriter(self.hits_path, _HIT_SCRATCH, compression="zstd", write_statistics=False)
        self.contigs_w = pq.ParquetWriter(self.contigs_path, _CONTIG_SCRATCH, compression="zstd", write_statistics=False)
        self.batches: dict[str, _Batch | None] = {"whole": None, "igr": None}
        self.stats = {"contigs": 0, "bases": 0, "batches": 0, "hits": 0, "regions": 0, "region_bases": 0}

    def _flush(self, kind: str) -> None:
        b = self.batches[kind]
        if b is None:
            return
        self.batches[kind] = None
        b.close()
        cm, models = (self.lib.whole_cm, self.lib.whole) if kind == "whole" else (self.lib.igr_cm, self.lib.igr)
        hits, contigs = _search(self.exe, cm, models, b, self.cpus, self.scratch, notrunc=kind == "igr")
        self.hits_w.write_batch(hits)
        self.contigs_w.write_batch(contigs)
        self.stats["batches"] += 1
        self.stats["hits"] += hits.num_rows
        b.path.unlink()

    def _add(self, kind: str, name: str, rec: bytes, what: _Rec) -> None:
        if self.batches[kind] is None:
            self.batches[kind] = _Batch(self.scratch / f"batch-{kind}.fna")
        b = self.batches[kind]
        b.add(name, rec, what)
        if b.bases >= self.batch_bases:
            self._flush(kind)

    def skip(self, cu: bytes, owner: str) -> None:
        self.contigs_w.write_batch(pa.RecordBatch.from_arrays(
            [pa.array([cu], pa.binary(16)), pa.array([0], pa.int64()), pa.array([0], pa.int64()),
             pa.array([f"{SKIP_REASON}{owner}"], pa.string())], schema=_CONTIG_SCRATCH))

    def contig(self, cu: bytes, seq: bytes, circular: bool, genes: list[tuple[int, int]]) -> None:
        n = len(seq)
        self.stats["contigs"] += 1
        self.stats["bases"] += n
        self.contigs_w.write_batch(pa.RecordBatch.from_arrays(
            [pa.array([cu], pa.binary(16)), pa.array([0], pa.int64()), pa.array([0], pa.int64()),
             pa.array([None], pa.string())], schema=_CONTIG_SCRATCH))
        if self.lib.whole:
            f = self.lib.flank_whole if circular else 0
            self._add("whole", cu.hex(), _extended(seq, f), _Rec(cu, n, 0, f, circular))
        if self.lib.igr:
            ring = seq + seq if circular else seq
            for s, e in regions(genes, n, circular, self.margin):
                self.stats["regions"] += 1
                self.stats["region_bases"] += e - s
                if circular and e - s == n:  # the whole circle: as a circular whole contig
                    f = self.lib.flank_igr
                    self._add("igr", f"{cu.hex()}_{s}", _extended(seq, f), _Rec(cu, n, 0, f, True))
                else:
                    self._add("igr", f"{cu.hex()}_{s}", ring[s:e], _Rec(cu, n, s, 0, circular))

    def close(self, ok: bool) -> None:
        try:
            if ok:
                for kind in ("whole", "igr"):
                    self._flush(kind)
        finally:
            for b in self.batches.values():
                if b is not None:
                    b.close()
            self.hits_w.close()
            self.contigs_w.close()


def _scan(root: Path, stage: str, lib: _Library, contigs: Iterator[tuple[bytes, bool, list]], stages: list[str],
          cpus: int, scratch: Path, batch_bases: int) -> tuple[_Scan, dict[bytes, str]]:
    """Search `contigs` ((contig_uid, circular, genes), sorted by contig_uid) into scratch; a contig another stage
    owns is not searched (DEC-26) and is returned with its owner."""
    exe = shutil.which("cmsearch")
    if exe is None:
        raise FileNotFoundError("cmsearch (Infernal) not on PATH")
    store = SeqStore(root, stages)
    owners = ContigOwners(root, stage, stages)
    skipped: dict[bytes, str] = {}
    sc = _Scan(exe, lib, cpus, scratch, batch_bases, MARGIN)
    ok = False
    try:
        for cu, circ, genes in contigs:
            s = owners.owner(cu)
            if s != stage:
                skipped[cu] = s
                sc.skip(cu, s)
                continue
            sc.contig(cu, store.get(cu).upper(), circ, genes)
        ok = True
    finally:
        sc.close(ok)
        store.close()
    return sc, skipped


def _mkscratch(prefix: str, tmp) -> Path:
    """A private scratch directory in `tmp` (made when it does not exist yet), else in the system's."""
    if tmp is not None:
        Path(tmp).mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix=prefix, dir=tmp))


def _check_args(cpus: int, batch_bases: int) -> None:
    if not isinstance(cpus, int) or isinstance(cpus, bool) or cpus < 1:
        raise ValueError(f"cpus {cpus!r}: give cmsearch's thread count explicitly (at least 1)")
    if batch_bases < 1:
        raise ValueError(f"batch_bases must be at least 1, got {batch_bases}")


def _tools(cm_path: Path) -> dict[str, str]:
    return {**provenance.tool_versions("infernal", "pyarrow"), LIBRARY: library_version(cm_path)}


def run(root: Path, stage: str, cm_path: Path, cpus: int, genomes: Iterable[GenomeInput] | None = None, *,
        model_set: str | Path = DEFAULT_SET, gene_call_version: str = GENE_CALL_VERSION, tmp=None,
        out: Path | None = None, seq_stages: list[str] | None = None, batch_bases: int = BATCH_BASES,
        max_file_bytes: int = layout.MAX_FILE_BYTES) -> dict:
    """Search the stage's contigs with the set `model_set` of the CM library `cm_path` and write its ncRNA features
    and hits (module docstring).

    `cpus`: cmsearch worker threads (`--cpu`), given explicitly; the output bytes do not depend on it.
    `genomes`: any iterable, consumed once; None reads the stage's `genome_input` table (which also says which
    contigs are circular). Contigs are searched in contig_uid order. `gene_call_version`: the gene calls whose
    intergenic regions the set's intergenic models search (needed only when the set has such models).
    `seq_stages`: the stages whose sequences the genomes may name (card B-28), default the stage itself. A contig
    that another stage stored first (DEC-26, card C-20) is that stage's: it is not searched here, and is listed in
    the result (`skipped_contigs`) and in the ledger (reason `owned_by:<stage>`, counts 0). `tmp` holds the private
    scratch directory (system temp by default) and may lie under neither the root nor `out`. `batch_bases`: bases
    per cmsearch call; the output does not depend on it. `max_file_bytes`: each table is cut into parts of at most
    this many estimated bytes.

    Inputs are read from `root`; every output and the ledger are written under `out` (default `root`; card C-08).
    Refuses, before reading a sequence, when one of the producer's partitions of this stage already holds a file in
    `root` or `out`.
    """
    root = Path(root)
    out = root if out is None else Path(out)
    _check_args(cpus, batch_bases)
    ledger.check_scratch(tmp, root, out)
    cm_path = Path(cm_path)
    prod = producer(cm_path, model_set)
    feat_dir, hit_dir, step_dir = _paths(stage, prod)
    ledger.refuse_written(root, out, [feat_dir, hit_dir, step_dir], f"ncRNA features {prod} of stage {stage}")
    members, _ = _check_set(cm_path, model_set)
    gene_files = _gene_files(root, stage, gene_call_version)
    if any(not m.whole_genome for m in members) and not gene_files:
        raise FileNotFoundError(f"stage {stage}: no gene calls of {gene_call_version}, whose intergenic regions the "
                                "set's intergenic models search; run gene calling first")
    tools = _tools(cm_path)
    stages = list(dict.fromkeys([stage, *([] if seq_stages is None else seq_stages)]))
    scratch = _mkscratch("biobox-ncrna-", tmp)
    try:
        lib = _fetch(cm_path, model_set, scratch)
        if genomes is None:
            gi = _genome_input_files(root, stage)
            if not gi:
                raise FileNotFoundError(f"stage {stage}: no genome_input")
            contigs = ((cu, bool(c)) for cu, c in _range_rows(gi, ["contig_uid", "is_circular"]))
        else:
            contigs = _distinct(iter(sorted({(cu, cu in set(g.circular)) for g in genomes for cu in g.contig_uids})))
        feats = _range_rows(gene_files, ["contig_uid", "start", "end"]) if lib.igr else None
        sc, skipped = _scan(root, stage, lib, _with_genes(contigs, feats), stages, cpus, scratch, batch_bases)
        _publish(out, stage, prod, [sc.hits_path], [sc.contigs_path], scratch, tools, max_file_bytes)
        return {**sc.stats, "producer": prod, "models": len(lib.whole) + len(lib.igr), "models_whole": len(lib.whole),
                "models_intergenic": len(lib.igr), "flank": lib.flank_whole, "flank_igr": lib.flank_igr,
                "skipped_contigs": [{"contig_uid": cu.hex(), "owner_stage": s} for cu, s in sorted(skipped.items())]}
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def _publish(out: Path, stage: str, prod: str, hit_files: list[Path], contig_files: list[Path], scratch: Path,
             tools: dict[str, str], max_file_bytes: int) -> None:
    """Sort the scratch hits into the producer's `feature` and `l1_ncrna` files and write the ledger: per contig
    searched, hits in and kept (and the hits longer than a circular contig rejected); a contig owned by another
    stage (DEC-26) has one row, reason `owned_by:<stage>`, counts 0. The one writer of `run` and `merge`, so their
    files are the same."""
    feat_dir, hit_dir, _ = _paths(stage, prod)
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb"), "preserve_insertion_order": True})
    try:
        con.execute(f"CREATE VIEW hits AS SELECT DISTINCT * FROM read_parquet({_sql_files(hit_files)})")
        con.execute("CREATE VIEW contigs AS SELECT contig_uid, sum(n_in) AS n_in, sum(n_rejected) AS n_rejected, "
                    f"max(reason) AS reason FROM read_parquet({_sql_files(contig_files)}) GROUP BY contig_uid")
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
            reader = con.sql("SELECT contig_uid, n_in, n_rejected, reason FROM contigs ORDER BY contig_uid"
                             ).to_arrow_reader(io.ROW_GROUP_ROWS)
            try:
                for b in reader:
                    for cu, n, rej, why in zip(*(b.column(i).to_pylist() for i in range(4))):
                        if why is not None:  # a contig another stage owns: not searched here
                            yield ledger.row(cu, "contig", 0, 0, reason=why, tools=tools)
                            continue
                        if rej:
                            yield ledger.row(cu, "contig", n, n - rej, n_rejected=rej, reason=REJECT_CIRCLE,
                                             tools=tools)
                        yield ledger.row(cu, "contig", n, n - rej, tools=tools)
            finally:
                reader.close()

        ledger.write(out, stage, f"ncrna-{prod}", ledger_rows())
    finally:
        con.close()


# --- the scan in shards (card E-07) ----------------------------------------------------------------------------


def scan_dir(stage: str, prod: str) -> str:
    """The root-relative directory of the shard plan and records of one producer and stage (the registry's
    `domain_scan` artifact, its library partition naming the producer)."""
    return f"{ARTIFACTS['domain_scan'].location}/library={prod}/stage={stage}"


def _plan_rel(stage: str, prod: str) -> str:
    return f"{scan_dir(stage, prod)}/plan.json"


def _shard_rel(stage: str, prod: str, k: int) -> str:
    return f"{scan_dir(stage, prod)}/shard-{k:05d}.json"


def _json_bytes(doc: dict) -> bytes:
    return json.dumps(doc, sort_keys=True, indent=1).encode() + b"\n"


def _write_once(root: Path, out: Path, rel: str, data: bytes, what: str) -> None:
    for base in dict.fromkeys([Path(root).resolve(), Path(out).resolve()]):
        if (base / rel).exists():
            raise FileExistsError(f"{base / rel}: {what} is written once; refusing to overwrite")
    path = Path(out) / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def _find(root: Path, out: Path | None, rel: str) -> Path | None:
    """`rel` in the root, else in `out` (written earlier by a direct call into another output root)."""
    for base in [Path(root)] + ([] if out is None else [Path(out)]):
        if (base / rel).is_file():
            return base / rel
    return None


def _plan_bytes(root: Path, stage: str, prod: str, out: Path | None = None) -> bytes:
    p = _find(root, out, _plan_rel(stage, prod))
    if p is None:
        raise FileNotFoundError(f"{Path(root) / _plan_rel(stage, prod)}: no shard plan for ncrna-{prod} of stage "
                                f"{stage}; write it with ncrna.plan")
    return p.read_bytes()


def read_plan(root: Path, stage: str, prod: str, out: Path | None = None) -> dict:
    """The stage's shard plan for one producer (from the root, else from `out`)."""
    return json.loads(_plan_bytes(root, stage, prod, out))


def run_paths(tmp: Path | None, root: Path, stage: str, prod: str, k: int, out: Path | None = None
              ) -> tuple[Path, Path]:
    """Where `scan_shard` writes shard k's sorted runs (hits, contigs): under the scratch directory (`tmp`, else the
    system's), in a directory named by the plan's sha256, so runs of another plan never mix in."""
    digest = hashlib.sha256(_plan_bytes(root, stage, prod, out)).hexdigest()[:16]
    d = Path(tmp if tmp is not None else tempfile.gettempdir()) / SCAN_DIR / digest
    return d / f"shard-{k:05d}.hits.parquet", d / f"shard-{k:05d}.contigs.parquet"


def _file_docs(root: Path, files: list[Path]) -> list[dict]:
    return [{"path": f.relative_to(root).as_posix(), "bytes": f.stat().st_size} for f in files]


def _cut(rows: Iterator[tuple[bytes, int]], scan_bases: int) -> Iterator[dict]:
    """Consecutive shards of at most `scan_bases` bases; a contig that would not fit closes the shard before it (a
    contig above the limit is a shard of its own)."""
    cur = None
    for cu, n in rows:
        if cur is not None and cur["bases"] + n > scan_bases:
            yield cur
            cur = None
        if cur is None:
            cur = {"first": cu.hex(), "last": cu.hex(), "contigs": 0, "bases": 0}
        cur["last"] = cu.hex()
        cur["contigs"] += 1
        cur["bases"] += n
    if cur is not None:
        yield cur


def _contig_bases(root: Path, gi: list[Path], idx: list[Path], scratch: Path) -> Iterator[tuple[bytes, int]]:
    """(contig_uid, bases to search) of the stage's `genome_input` contigs in contig_uid order: the contig's length
    in the stage's `contig_index` (default variant), 0 for a contig another stage stored (DEC-26)."""
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb")})
    try:
        sql = (f"SELECT g.contig_uid, coalesce(c.contig_len, 0) AS n FROM (SELECT DISTINCT contig_uid FROM "
               f"read_parquet({_sql_files(gi)}, hive_partitioning = false)) g LEFT JOIN (SELECT contig_uid, "
               f"max(contig_len) AS contig_len FROM read_parquet({_sql_files(idx)}, hive_partitioning = false) WHERE "
               "is_default_variant GROUP BY contig_uid) c ON c.contig_uid = g.contig_uid ORDER BY g.contig_uid")
        reader = con.sql(sql).to_arrow_reader(_READ_ROWS)
        try:
            for b in reader:
                yield from zip(b.column(0).to_pylist(), b.column(1).to_pylist())
        finally:
            reader.close()
    finally:
        con.close()


def plan(root: Path, stage: str, cm_path: Path, scan_bases: int, *, model_set: str | Path = DEFAULT_SET,
         gene_call_version: str = GENE_CALL_VERSION, tmp=None, out: Path | None = None) -> dict:
    """Write the stage's shard plan under `out` (default `root`): its `genome_input` contigs sorted by contig_uid,
    cut into shards of at most `scan_bases` genome bases. Refuses when the stage's features of this producer are
    written already, or a plan exists, in `root` or `out`. Returns the shard count and the totals."""
    root = Path(root)
    out = root if out is None else Path(out)
    if not isinstance(scan_bases, int) or isinstance(scan_bases, bool) or scan_bases < 1:
        raise ValueError(f"scan_bases {scan_bases!r}: must be a positive integer")
    ledger.check_scratch(tmp, root, out)
    cm_path = Path(cm_path)
    prod = producer(cm_path, model_set)
    ledger.refuse_written(root, out, list(_paths(stage, prod)), f"ncRNA features {prod} of stage {stage}")
    rel = _plan_rel(stage, prod)
    for base in dict.fromkeys([root.resolve(), out.resolve()]):
        if (base / rel).exists():
            raise FileExistsError(f"{base / rel}: the shard plan is written once; refusing to overwrite")
    members, _ = _check_set(cm_path, model_set)
    gi = _genome_input_files(root, stage)
    idx = io.part_files(root / INDEX_DIR / f"stage={stage}")
    genes = _gene_files(root, stage, gene_call_version) if any(not m.whole_genome for m in members) else []
    if not gi or not idx:
        raise FileNotFoundError(f"stage {stage}: no genome_input or contig_index files")
    if any(not m.whole_genome for m in members) and not genes:
        raise FileNotFoundError(f"stage {stage}: no gene calls of {gene_call_version}; run gene calling first")
    scratch = _mkscratch("biobox-ncrna-plan-", tmp)
    try:
        shards = list(_cut(_contig_bases(root, gi, idx, scratch), scan_bases))
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    doc = {"format": PLAN_FORMAT, "stage": stage, "producer": prod, "gene_call_version": gene_call_version,
           "cm_sha256": _sha256(cm_path), "set_sha256": _sha256(set_path(model_set)), "margin": MARGIN,
           "scan_bases": scan_bases, "files": {"genome_input": _file_docs(root, gi), "contig_index": _file_docs(root, idx),
                                               "genes": _file_docs(root, genes)},
           "contigs": sum(s["contigs"] for s in shards), "bases": sum(s["bases"] for s in shards), "shards": shards}
    _write_once(root, out, rel, _json_bytes(doc), "the shard plan")
    return {"shards": len(shards), "contigs": doc["contigs"], "bases": doc["bases"], "producer": prod}


def _checked_files(root: Path, docs: list[dict]) -> list[Path]:
    out = []
    for f in docs:
        path = Path(root) / f["path"]
        if path.stat().st_size != f["bytes"]:
            raise ValueError(f"{path}: {path.stat().st_size} bytes, the shard plan lists {f['bytes']}")
        out.append(path)
    return out


def _shard_scan(root: Path, stage: str, cm_path: Path, model_set, k: int, cpus: int, seq_stages, scratch: Path,
                out: Path | None) -> tuple[bytes, dict, _Scan, dict[bytes, str], str]:
    """Scan shard k of the plan into `scratch`, checked: k is one of the plan's shards, the library and the set are
    the plan's, and the contigs and bases read are the shard's."""
    cm_path = Path(cm_path)
    prod = producer(cm_path, model_set)
    raw = _plan_bytes(root, stage, prod, out)
    p = json.loads(raw)
    if not 0 <= k < len(p["shards"]):
        raise ValueError(f"shard {k}: the plan of ncrna-{prod} for stage {stage} has {len(p['shards'])} shards")
    for what, have, want in (("library", _sha256(cm_path), p["cm_sha256"]),
                             ("model set", _sha256(set_path(model_set)), p["set_sha256"])):
        if have != want:
            raise ValueError(f"the {what}'s sha256 {have} is not the plan's ({want})")
    s = p["shards"][k]
    lo, hi = bytes.fromhex(s["first"]), bytes.fromhex(s["last"])
    gi = _checked_files(root, p["files"]["genome_input"])
    genes = _checked_files(root, p["files"]["genes"])
    lib = _fetch(cm_path, model_set, scratch)
    read = {"contigs": 0}

    def contigs():
        for cu, c in _range_rows(gi, ["contig_uid", "is_circular"], lo, hi):
            read["contigs"] += 1
            yield cu, bool(c)

    feats = _range_rows(genes, ["contig_uid", "start", "end"], lo, hi) if lib.igr else None
    stages = list(dict.fromkeys([stage, *([] if seq_stages is None else seq_stages)]))
    sc, skipped = _scan(root, stage, lib, _with_genes(_distinct(contigs()), feats), stages, cpus, scratch,
                        BATCH_BASES)
    if read["contigs"] < s["contigs"] or sc.stats["contigs"] + len(skipped) != s["contigs"] \
            or sc.stats["bases"] != s["bases"]:
        raise ValueError(f"shard {k}: {sc.stats['contigs'] + len(skipped)} contigs ({sc.stats['bases']} bases) read, "
                         f"the plan lists {s['contigs']} ({s['bases']} bases)")
    return raw, s, sc, skipped, prod


def _distinct(contigs: Iterator[tuple[bytes, bool]]) -> Iterator[tuple[bytes, bool]]:
    last = None
    for cu, c in contigs:
        if cu != last:
            yield cu, c
        last = cu


def _sorted_runs(sc: _Scan, scratch: Path, hits_path: Path, contigs_path: Path) -> None:
    """The scan's scratch files as two sorted runs, so their bytes depend on the rows only (not on the batches or
    the threads): the distinct hits by the table's key, the contig counts summed per contig."""
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb"), "preserve_insertion_order": True})
    try:
        cols = ", ".join(f'"{c}"' for c in _HIT_SCRATCH.names)
        for sql, schema, path in (
                (f"SELECT DISTINCT {cols} FROM read_parquet({_sql_path(sc.hits_path)}) "
                 'ORDER BY contig_uid, start, "end", strand, model_acc, bitscore', _HIT_SCRATCH, hits_path),
                (f"SELECT contig_uid, sum(n_in)::BIGINT, sum(n_rejected)::BIGINT, max(reason) FROM "
                 f"read_parquet({_sql_path(sc.contigs_path)}) GROUP BY contig_uid ORDER BY contig_uid",
                 _CONTIG_SCRATCH, contigs_path)):
            tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
            w = pq.ParquetWriter(tmp, schema, compression="zstd", write_statistics=False)
            reader = con.sql(sql).to_arrow_reader(io.ROW_GROUP_ROWS)
            try:
                for b in reader:  # DuckDB returns UIDs as BLOB: cast back to the scratch schema
                    w.write_batch(pa.RecordBatch.from_arrays(
                        [c.cast(f.type) for c, f in zip(b.columns, schema)], schema=schema))
            finally:
                reader.close()
                w.close()
            os.replace(tmp, path)
    finally:
        con.close()


def scan_shard(root: Path, stage: str, cm_path: Path, k: int, cpus: int, *, model_set: str | Path = DEFAULT_SET,
               seq_stages: list[str] | None = None, tmp=None, out: Path | None = None) -> dict:
    """Scan shard `k` of the stage's plan with `cpus` threads: write its sorted runs under scratch (`run_paths`) and
    its record (the runs' sha256 and counts) under `out` (default `root`). Refuses a shard whose record exists in
    `root` or `out`, a library or set other than the plan's, and contigs that do not match the plan's counts."""
    root = Path(root)
    out = root if out is None else Path(out)
    _check_args(cpus, BATCH_BASES)
    ledger.check_scratch(tmp, root, out)
    prod = producer(Path(cm_path), model_set)
    ledger.refuse_written(root, out, list(_paths(stage, prod)), f"ncRNA features {prod} of stage {stage}")
    rel = _shard_rel(stage, prod, k)
    for base in dict.fromkeys([root.resolve(), out.resolve()]):
        if (base / rel).exists():
            raise FileExistsError(f"{base / rel}: shard {k} is scanned once; refusing to overwrite")
    scratch = _mkscratch("biobox-ncrna-", tmp)
    try:
        raw, s, sc, skipped, prod = _shard_scan(root, stage, cm_path, model_set, k, cpus, seq_stages, scratch, out)
        hits_path, contigs_path = run_paths(tmp, root, stage, prod, k, out)
        hits_path.parent.mkdir(parents=True, exist_ok=True)
        _sorted_runs(sc, scratch, hits_path, contigs_path)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    rec = {"format": PLAN_FORMAT, "plan_sha256": hashlib.sha256(raw).hexdigest(), "shard": k, "first": s["first"],
           "last": s["last"], "contigs": s["contigs"], "bases": s["bases"], "hits": sc.stats["hits"],
           "skipped": len(skipped), "hits_sha256": ids.sha256_file(hits_path).hex(),
           "contigs_sha256": ids.sha256_file(contigs_path).hex()}
    _write_once(root, out, rel, _json_bytes(rec), f"the record of shard {k}")
    return {"shard": k, "contigs": s["contigs"], "bases": s["bases"], "hits": sc.stats["hits"],
            "regions": sc.stats["regions"], "batches": sc.stats["batches"]}


def shard_record(root: Path, stage: str, prod: str, k: int) -> dict | None:
    """The published record of shard k (in the root), or None when there is none."""
    p = _find(root, None, _shard_rel(stage, prod, k))
    return None if p is None else json.loads(p.read_bytes())


def run_problem(root: Path, stage: str, prod: str, k: int, tmp: Path | None) -> str | None:
    """Why the scratch runs of shard k, whose record is published, cannot be merged as they are, or None: "is
    missing" or "differs from the sha256 its record holds"."""
    rec = shard_record(root, stage, prod, k)
    if rec is None:
        raise FileNotFoundError(f"shard {k} of ncrna-{prod}, stage {stage}: no record")
    for path, key in zip(run_paths(tmp, root, stage, prod, k), ("hits_sha256", "contigs_sha256")):
        if not path.is_file():
            return "is missing"
        if ids.sha256_file(path).hex() != rec[key]:
            return "differs from the sha256 its record holds"
    return None


def rescan_shard(root: Path, stage: str, cm_path: Path, k: int, cpus: int, *, model_set: str | Path = DEFAULT_SET,
                 seq_stages: list[str] | None = None, tmp=None) -> dict:
    """Scan shard k again into scratch (`run_paths`) for a shard whose record is published in `root` but whose runs
    are gone or damaged, e.g. after a reboot cleared /tmp (as card C-12 does for domains). The record is not written
    again: the new runs must have the sha256 it holds. On a mismatch the new runs are deleted, the runs in scratch
    are left as they were, and a ValueError names the shard."""
    root = Path(root)
    ledger.check_scratch(tmp, root)
    _check_args(cpus, BATCH_BASES)
    prod = producer(Path(cm_path), model_set)
    ledger.refuse_written(root, root, list(_paths(stage, prod)), f"ncRNA features {prod} of stage {stage}")
    rec = shard_record(root, stage, prod, k)
    if rec is None:
        raise FileNotFoundError(f"shard {k} of ncrna-{prod}, stage {stage}: no published record to scan it again for; "
                                "scan it with ncrna.scan_shard")
    scratch = _mkscratch("biobox-ncrna-", tmp)
    try:
        raw, _, sc, _, _ = _shard_scan(root, stage, cm_path, model_set, k, cpus, seq_stages, scratch, None)
        if rec["plan_sha256"] != hashlib.sha256(raw).hexdigest() or rec["shard"] != k:
            raise ValueError(f"shard {k} of ncrna-{prod}, stage {stage}: its record is of another plan or shard")
        paths = run_paths(tmp, root, stage, prod, k)
        paths[0].parent.mkdir(parents=True, exist_ok=True)
        new = [p.with_name(f".rescan-{os.getpid()}-{p.name}") for p in paths]
        _sorted_runs(sc, scratch, *new)
        have = [ids.sha256_file(p).hex() for p in new]
        if have != [rec["hits_sha256"], rec["contigs_sha256"]]:
            for p in new:
                p.unlink()
            raise ValueError(f"shard {k} of ncrna-{prod}, stage {stage}: scanned again, its runs have sha256 {have}, "
                             "the record holds others; the published record is not the scan of these inputs, so the "
                             "merge cannot run (the runs in scratch are left as they were)")
        for a, b in zip(new, paths):
            os.replace(a, b)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    return {"shard": k, "hits": sc.stats["hits"], "hits_sha256": have[0], "contigs_sha256": have[1]}


def merge(root: Path, stage: str, cm_path: Path, *, model_set: str | Path = DEFAULT_SET, tmp=None,
          max_file_bytes: int = layout.MAX_FILE_BYTES, out: Path | None = None) -> dict:
    """Write the shards' runs as the stage's ncRNA files and ledger under `out` (default `root`): the files `run`
    writes. Every shard needs its record (in `root` or `out`) and runs whose sha256 are the record's; otherwise
    nothing is written."""
    root = Path(root)
    out_root = root if out is None else Path(out)
    ledger.check_scratch(tmp, root, out_root)
    cm_path = Path(cm_path)
    prod = producer(cm_path, model_set)
    ledger.refuse_written(root, out_root, list(_paths(stage, prod)), f"ncRNA features {prod} of stage {stage}")
    raw = _plan_bytes(root, stage, prod, out)
    p, digest = json.loads(raw), hashlib.sha256(raw).hexdigest()
    missing = [k for k in range(len(p["shards"])) if _find(root, out, _shard_rel(stage, prod, k)) is None]
    if missing:
        raise FileNotFoundError(f"ncrna-{prod} of stage {stage}: no record of shard {', '.join(map(str, missing))} "
                                f"of {len(p['shards'])}; scan it with ncrna.scan_shard")
    hits, contigs, n_hits = [], [], 0
    for k in range(len(p["shards"])):
        rec = json.loads(_find(root, out, _shard_rel(stage, prod, k)).read_bytes())
        if rec["plan_sha256"] != digest or rec["shard"] != k:
            raise ValueError(f"shard {k}: its record is of another plan or shard")
        for path, key in zip(run_paths(tmp, root, stage, prod, k, out), ("hits_sha256", "contigs_sha256")):
            if not path.is_file():
                raise FileNotFoundError(f"{path}: the run of shard {k} is not in the scratch directory; its record is "
                                        "published, so scan the shard again into this scratch directory")
            if ids.sha256_file(path).hex() != rec[key]:
                raise ValueError(f"{path}: sha256 differs from the record of shard {k}")
        h, c = run_paths(tmp, root, stage, prod, k, out)
        hits.append(h)
        contigs.append(c)
        n_hits += rec["hits"]
    scratch = _mkscratch("biobox-ncrna-merge-", tmp)
    try:
        _publish(out_root, stage, prod, hits, contigs, scratch, _tools(cm_path), max_file_bytes)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    return {"producer": prod, "shards": len(p["shards"]), "contigs": p["contigs"], "bases": p["bases"], "hits": n_hits}


def main(argv: list[str] | None = None) -> int:
    """`python -m biobox.l1.ncrna rescan --root R --stage S --cm F [--ncrna-set SET] --shard K --cpus N --tmp D
    [--seq-stages S0,S1]`: `rescan_shard` in its own process, as build-stage runs it. It writes nothing into the root,
    so it is not a run-step. Prints one JSON line; exit 2 with the reason on a refusal."""
    ap = argparse.ArgumentParser(prog="python -m biobox.l1.ncrna")
    ap.add_argument("action", choices=["rescan"])
    for name in ("--root", "--stage", "--cm", "--tmp"):
        ap.add_argument(name, required=True)
    ap.add_argument("--ncrna-set", default=DEFAULT_SET)
    ap.add_argument("--seq-stages")
    ap.add_argument("--shard", type=int, required=True)
    ap.add_argument("--cpus", type=int, required=True)
    a = ap.parse_args(argv)
    try:
        res = rescan_shard(Path(a.root), a.stage, Path(a.cm), a.shard, a.cpus, model_set=a.ncrna_set,
                           seq_stages=None if a.seq_stages is None else a.seq_stages.split(","), tmp=Path(a.tmp))
    except (ValueError, FileNotFoundError, FileExistsError) as e:
        print(f"ncrna rescan: {e}", file=sys.stderr)
        return 2
    print(json.dumps({"action": a.action, "result": res, "stage": a.stage}, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
