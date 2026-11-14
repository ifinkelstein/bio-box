"""CRISPR arrays as features, called by MinCED (card E-12, `PLAN.md` step 4; operation O-ARRAY).

The tool. MinCED 0.4.2 (pinned as `minced=0.4.2` in `container/environment.yml`), a Java program derived from CRT
(Bland et al. 2007): it scans each FASTA record for runs of a direct repeat separated by spacers. It has no threads,
no random seed and no notion of topology. It is run as `java <JVM_OPTIONS> -jar minced.jar <MINCED_ARGS> -gffFull
batch.fna batch.txt batch.gff`; the jar is the one the conda wrapper `minced` runs (the wrapper's own directory,
symlinks resolved), and `java` is `$JAVA_HOME/bin/java` when it exists, else the one on PATH, as the wrapper chooses.
Calling java directly is what lets the JVM's settings be fixed (below); the wrapper takes none.

Parameters (`MINCED_ARGS`, given explicitly, so a change of MinCED's defaults cannot change the calls silently; they
are MinCED's and CRT's own defaults, which Prokka also runs):

- `-minNR 3`: at least three repeats (two spacers). Two copies of a 23 to 47 base direct repeat at spacer distance
  occur by chance and in repeated genes; a two-repeat array is not called by MinCED's default, and lowering it is a
  change of producer that needs its false-positive rate measured first (later card).
- `-minRL 23 -maxRL 47`: repeat length range, CRT's range from the repeats of CRISPRdb.
- `-minSL 26 -maxSL 50`: spacer length range, CRT's. Spacers of type III systems up to about 70 bases fall outside it
  (stated as a limit; a later card may widen it with the false-positive rate measured).
- `-searchWL 8`: the length of the k-mer search window that seeds repeat candidates (range 6 to 9).

The producer is `minced-<version>-w8-nr3-r23x47-s26x50` (`PARAM_TAG`): MinCED's version, read from the jar
(`--version`), and the parameters; a change of either, or of the origin rule below, is a new producer, so new
partitions and, should arrays go on a tape, a new tape version.

Output format. With `-gffFull` MinCED writes a GFF3 file: per array a `repeat_region` line (1-based inclusive range,
score = repeats, `ID=CRISPR<k>`, `rpt_unit_seq` the consensus repeat) followed by one `repeat_unit` line per repeat
(`Parent=CRISPR<k>`); `##gff-version 3` is repeated before each record's lines. That file is parsed (`parse_gff`):
it has the record name on every line, no timing and no layout. The text report (the first output file, a scratch file
deleted after; MinCED deletes an existing output file before it writes, so it cannot be `/dev/null`) is not parsed:
with `-gffFull` it carries no record names at all, and without it, it carries the time each record took (`Time to
find repeats: 456 ms`), so it is neither machine-readable nor byte-stable. Spacers are the gaps
between consecutive repeat units, which is how CRT defines them (its text report prints the same bases, as the test
oracle checks); repeat and spacer sequences are cut from the bases the producer gave MinCED.

Coordinates. MinCED reports 1-based inclusive positions; a range `[a, b]` is `[a - 1, b)` here (0-based half-open,
AGENTS §3.1), checked on a hand-made contig (`tests/r1/test_crispr_parse.py`).

Circular contigs (MinCED knows no topology). A linear contig is one record. A circular contig of length n is one
record holding F bases of its end, the contig, and F bases of its start, F = min(n, `CIRC_FLANK`), as the ncRNA
producer searches a circular whole contig; an array is kept when its leftmost base lies on the contig itself (record
offsets [F, F + n)), so each array, also one across the origin, is found once with its whole context. An array across
the origin has start < n < end (AGENTS §3.1); an element (repeat or spacer) of it that lies wholly past the origin is
stored at its position on the contig (start - n), and one that crosses the origin has end > n, as an ncRNA hit is.
An array longer than the circle is impossible on the molecule and is counted in the ledger as rejected
(`exceeds_circle`). An array kept on a circular record whose end comes within one repeat and one spacer of the
record's end may have been cut by the record; when F = n it is an array as long as the circle (rejected as above),
when F < n the run fails (RuntimeError): no array of `CIRC_FLANK` bases is known (*estimate*: the longest arrays
reported hold a few hundred repeats, under 50 kb), so it is a finding, not a case to absorb. CRT's scan works left to
right on each record and finds an array from its own bases and its neighbours', so an array's call does not depend on
where in the record it lies or on the other records of the batch (tested: a contig alone and in a batch, a planted
array on a linear contig and rotated across the origin of a circular one give the same elements).

Strand. MinCED reports repeats as they read on the + strand of the record and does not say which way the array is
transcribed. Every array is stored with strand +1 (the strand its repeat and spacer sequences are read on), and so
is each element. The orientation (CRISPRDirection, PLAN.md step 4) is a later card; it would be a new nullable column
of `l1_crispr` (LY-09), not a change of strand, so a `feature_uid` (which hashes the strand, DEC-18) does not change
when the orientation becomes known. A spacer as a search target is the + strand sequence and its reverse complement.

Determinism (AGENTS §4.3). MinCED's GFF depends only on the records (tested: byte-identical across runs). The JVM is
given a fixed heap (`-Xms` = `-Xmx` = `JAVA_HEAP`), the serial collector, one active processor, no perf-data file
(which it would write under /tmp) and the scratch directory as its temporary directory (`JVM_OPTIONS`); none of
these changes MinCED's output, they fix its footprint. Up to `cpus` MinCED processes run at once, each on its own
batch file; the rows go to unsorted scratch Parquet files and one DuckDB connection (one thread) sorts them into the
outputs, so the bytes depend neither on `cpus`, nor on the batches, nor on the genome order.

Writes, under the output root `out` (the data root by default; `$BIOBOX_STAGING` under run-step, card C-08), through
`biobox.io.write_parts`:

    l1/feature/producer=<p>/stage=<S>/part-<k>.parquet   kind crispr_array, one row per array, sorted by
                                                          (contig_uid, start, end, strand)
    l1/crispr/producer=<p>/stage=<S>/part-<k>.parquet    one row per repeat and per spacer (registry
                                                          `l1_crispr`), sorted by (contig_uid, feature_uid, ordinal)
    ledger/stage=<S>/step=crispr-<p>/part-0.parquet       per contig: arrays in, kept, rejected (longer than a
                                                          circular contig); tool versions: minced, java, pyarrow
                                                          and the arguments (`minced_args`)

`l1_crispr` rows: the array's `feature_uid` (the stable key), the element's `region_uid` (`ids.region_uid` of its
contig, start, end and strand: the key a spacer will have as a search target, spec invariant 16), `element`
(`repeat` or `spacer`), `ordinal` (its position in the array from 0, repeats even and spacers odd; valid only beside
the array's `feature_uid`, never a key, AGENTS §3.3), start, end, strand, `seq` (the element's bases, uppercase, on
the stored strand) and `contig_len`. The sequence is stored, not only its coordinates or a digest: a spacer is 26 to
50 bases, so a row is small (*estimate* 80 bytes zstd), and the spacer corpus and any protospacer search read the
spacers without opening the sequence store; the test checks every `seq` against the contig's bases.

A contig that another stage stored first (DEC-26, card C-20) is that stage's: it is not searched here and has one
ledger row, reason `owned_by:<stage>`, counts 0. A contig in several genomes of the stage (`genome_input` lists it
once per genome, DEC-27, card C-21a) is searched once, for the genome with the least `genome_uid`; beside its `ok` row
it has a ledger row with reason `owned_by:<stage>/<genome_uid>` and the same counts, and the result lists it with
`owner_genome`. As the ncRNA producer does (its helpers are used).

Tapes. Arrays are `feature` rows of their own producer and go on a tape only when the tape step is given that
producer (`--producers <gene caller>,<this producer>`), which makes a tape version of its own, as the ncRNA producer's
features do (tested). `build-stage --crispr` does not pass it yet: the C-08 provenance check
(`biobox.provenance._producer_problems`) knows only gene-call and ncRNA producer names, so a release whose tape
version lists this producer would fail it; that and what arrays as slots mean for `nbr` windows are a later card.

Bounded memory (AGENTS §4.1): contigs are streamed in `contig_uid` order from the stage's `genome_input`, one contig's
sequence at a time, appended to a FASTA batch file in a private scratch directory; a batch of at least `batch_bases`
bases (or the last) is one MinCED process with the fixed heap. Only a batch's parsed rows are held (its arrays, a
few hundred rows), then appended to the scratch Parquet files. Nothing is published if a batch fails.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from collections.abc import Iterable, Iterator
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from biobox import ids, io, layout, ledger, provenance
from biobox.l0seq.inputs import GenomeInput
from biobox.l0seq.store import OWNER_SEP, ContigOwners, SeqStore, skipped_entry
from biobox.l0seq.store import OWNED_BY as SKIP_REASON
from biobox.l1.ncrna import _distinct, _genome_input_files, _mkscratch, _range_rows, _sql_files
from biobox.registry.schemas import table

__all__ = ["CIRC_FLANK", "JAVA_HEAP", "MINCED_ARGS", "PARAM_TAG", "Array", "elements", "minced_version", "parse_gff",
           "place", "producer", "run"]

KIND_CRISPR = ids.KIND["crispr_array"]
STRAND = 1  # every array and element: the strand its sequences are read on (module docstring)
# MinCED's parameters, in the order they are passed (module docstring). Its defaults, given explicitly
MINCED_ARGS = (("-searchWL", 8), ("-minNR", 3), ("-minRL", 23), ("-maxRL", 47), ("-minSL", 26), ("-maxSL", 50))
PARAM_TAG = "w8-nr3-r23x47-s26x50"  # the parameters above in the producer name
MAX_UNIT = 47 + 50  # the longest repeat plus the longest spacer MinCED reports with MINCED_ARGS
# bases of each circular contig's end and start around it in its record (module docstring); an estimate of a bound,
# not measured: no reported CRISPR array is half as long
CIRC_FLANK = 100_000
JAVA_HEAP = "512m"  # MinCED's fixed heap. A 13.4 Mbp genome ran in 256 MiB (147 MB RSS; measured, wilkcomp02,
# 2026-11-13, GCA_037157765.1 of stage P2), and the GFF was byte-identical to a run with 1 GiB; 512 MiB leaves room
JVM_OPTIONS = (f"-Xms{JAVA_HEAP}", f"-Xmx{JAVA_HEAP}", "-XX:+UseSerialGC", "-XX:ActiveProcessorCount=1",
               "-XX:-UsePerfData")
# bases of records per MinCED process: chosen, not tuned. A process costs a JVM start (about 0.5 CPU-s, measured in
# the note), so a batch should be large beside that; one record is held at a time, whatever the batch
BATCH_BASES = 50_000_000
DUCK_MEMORY_LIMIT = "2GB"  # the sorting connection's cap; larger sorts spill into the scratch directory
# estimated zstd bytes per row; they only cut parts at `max_file_bytes` (card B-19). Not measured
EST_FEATURE_ROW = 24  # estimate: feature_uid 16 + delta-coded coordinates (as genecall's)
EST_ELEMENT_ROW = 80  # estimate: two UIDs 32 + coordinates + a 26 to 50 base sequence
REJECT_CIRCLE = "exceeds_circle"

_ELEMENT_SCRATCH = pa.schema([f for f in table("l1_crispr").schema])
_ARRAY_SCRATCH = pa.schema([("feature_uid", pa.binary(16)), ("contig_uid", pa.binary(16)), ("start", pa.int64()),
                            ("end", pa.int64()), ("strand", pa.int8())])
_CONTIG_SCRATCH = pa.schema([("contig_uid", pa.binary(16)), ("n_in", pa.int64()), ("n_rejected", pa.int64()),
                             ("reason", pa.string())])
_VERSION_RE = re.compile(r"^minced (\S+)", re.M)
_JAVA_RE = re.compile(r'version "([^"]+)"')


# --- the tool ---------------------------------------------------------------------------------------------------


def _tool() -> tuple[str, Path]:
    """(java, minced.jar): the jar beside the conda wrapper `minced` (symlinks resolved), and `$JAVA_HOME/bin/java`
    when it exists, else `java` on PATH, as the wrapper picks them."""
    wrapper = shutil.which("minced")
    if wrapper is None:
        raise FileNotFoundError("minced not on PATH (it is in the bio-box image, container/environment.yml)")
    jar = Path(os.path.realpath(wrapper)).parent / "minced.jar"
    if not jar.is_file():
        raise FileNotFoundError(f"{jar}: no minced.jar beside the minced wrapper {wrapper}")
    home = os.environ.get("JAVA_HOME")
    java = str(Path(home) / "bin" / "java") if home and (Path(home) / "bin" / "java").exists() else shutil.which("java")
    if java is None:
        raise FileNotFoundError("java not on PATH (MinCED is a Java program)")
    return java, jar


@lru_cache(maxsize=4)
def _versions(java: str, jar: str) -> tuple[str, str]:
    """(MinCED's version from `--version`, the JVM's from `java -version`), once per process."""
    m = subprocess.run([java, *JVM_OPTIONS, "-jar", jar, "--version"], capture_output=True, text=True,
                       stdin=subprocess.DEVNULL)
    v = _VERSION_RE.search(m.stdout + m.stderr)
    if v is None:
        raise ValueError(f"{jar}: no `minced <version>` line in the output of --version")
    j = subprocess.run([java, "-version"], capture_output=True, text=True, stdin=subprocess.DEVNULL)
    jv = _JAVA_RE.search(j.stderr + j.stdout)
    return v[1], (jv[1] if jv else "unknown")


def minced_version() -> str:
    java, jar = _tool()
    return _versions(java, str(jar))[0]


def producer() -> str:
    """`minced-<version>-<PARAM_TAG>` (module docstring)."""
    return f"minced-{minced_version()}-{PARAM_TAG}"


def minced_args() -> list[str]:
    return [str(x) for kv in MINCED_ARGS for x in kv]


def _argv(java: str, jar: Path, fna: Path, txt: Path, gff: Path, scratch: Path) -> list[str]:
    return [java, *JVM_OPTIONS, f"-Djava.io.tmpdir={scratch}", "-jar", str(jar), *minced_args(), "-gffFull",
            str(fna), str(txt), str(gff)]


def _tools(java: str, jar: Path) -> dict[str, str]:
    """The ledger's tool versions: MinCED's, the JVM's, pyarrow's, and MinCED's arguments."""
    mv, jv = _versions(java, str(jar))
    return {"minced": mv, "java": jv, **provenance.tool_versions("pyarrow"), "minced_args": " ".join(minced_args())}


# --- parsing and placing ----------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Array:
    """One array as MinCED reports it on a record: 0-based half-open record offsets of the array and of each repeat
    unit, in order."""
    record: str
    start: int
    end: int
    units: tuple[tuple[int, int], ...]


def _attrs(field: str) -> dict[str, str]:
    return dict(kv.split("=", 1) for kv in field.split(";") if "=" in kv)


def parse_gff(lines: Iterable[str]) -> Iterator[Array]:
    """The arrays of a MinCED `-gffFull` file, in file order. Each `repeat_region` line opens an array and the
    `repeat_unit` lines after it (`Parent=` its ID) are its repeats; positions 1-based inclusive become 0-based
    half-open. ValueError when a unit names another parent, units overlap or are out of order, the region's range is
    not from its first unit's start to its last unit's end, or its score is not its number of units."""
    cur: dict | None = None

    def close() -> Array:
        units = tuple(cur["units"])
        if not units:
            raise ValueError(f"array {cur['id']} on {cur['record']}: no repeat_unit")
        if len(units) != cur["score"]:
            raise ValueError(f"array {cur['id']} on {cur['record']}: score {cur['score']}, {len(units)} repeat units")
        if (units[0][0], units[-1][1]) != (cur["start"], cur["end"]):
            raise ValueError(f"array {cur['id']} on {cur['record']}: range [{cur['start']}, {cur['end']}) is not its "
                             f"units' [{units[0][0]}, {units[-1][1]})")
        for (a, b), (c, _) in zip(units, units[1:]):
            if c <= b:
                raise ValueError(f"array {cur['id']} on {cur['record']}: repeat units [{a}, {b}) and one at {c} "
                                 "overlap, touch or are out of order")
        return Array(cur["record"], cur["start"], cur["end"], units)

    for k, line in enumerate(lines, 1):
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) != 9:
            raise ValueError(f"GFF line {k}: {len(f)} fields, not 9: {line!r}")
        rec, kind, a, b, attrs = f[0], f[2], int(f[3]), int(f[4]), _attrs(f[8])
        if a < 1 or b < a:
            raise ValueError(f"GFF line {k}: range {a}..{b}")
        if kind == "repeat_region":
            if cur is not None:
                yield close()
            cur = {"record": rec, "id": attrs.get("ID"), "start": a - 1, "end": b, "score": int(f[5]), "units": []}
        elif kind == "repeat_unit":
            if cur is None or attrs.get("Parent") != cur["id"] or rec != cur["record"]:
                raise ValueError(f"GFF line {k}: a repeat_unit outside its repeat_region: {line!r}")
            cur["units"].append((a - 1, b))
        else:
            raise ValueError(f"GFF line {k}: feature type {kind!r}, not repeat_region or repeat_unit")
    if cur is not None:
        yield close()


def elements(units: tuple[tuple[int, int], ...]) -> list[tuple[str, int, int, int]]:
    """(element, ordinal, start, end) of an array's repeats and the spacers between them, in order: repeat i has
    ordinal 2i, the spacer after it 2i + 1; a spacer is the gap from one repeat's end to the next one's start."""
    out = []
    for i, (a, b) in enumerate(units):
        out.append(("repeat", 2 * i, a, b))
        if i + 1 < len(units):
            out.append(("spacer", 2 * i + 1, b, units[i + 1][0]))
    return out


@dataclass(frozen=True)
class Rec:
    """What one FASTA record is: its contig, the contig's length, the flank of a circular contig's record (the
    contig starts at record offset `flank`) and whether the contig is circular."""
    contig: bytes
    n: int
    flank: int
    circular: bool
    length: int  # bases in the record


def place(arr: Array, rec: Rec) -> tuple[int, int] | str | None:
    """Where an array of a record lies on its contig: (start, end) 0-based half-open, start < n, end > n for one
    across a circular contig's origin; None for a copy in a circular record's flank (the leftmost-base rule);
    `REJECT_CIRCLE` for an array longer than the circle. RuntimeError for an array a circular record may have cut
    (module docstring)."""
    start, end = arr.start - rec.flank, arr.end - rec.flank
    if not rec.circular:
        return start, end
    if not 0 <= start < rec.n:
        return None
    if end - start > rec.n or arr.end > rec.length - MAX_UNIT:
        if end - start > rec.n or rec.flank == rec.n:
            return REJECT_CIRCLE
        raise RuntimeError(f"contig {rec.contig.hex()}: an array at [{start}, {end}) reaches within {MAX_UNIT} bases "
                           f"of its record's end, so the record ({rec.flank} bases each side, CIRC_FLANK) may have cut "
                           "it; an array this long is not known, raise CIRC_FLANK in a new producer")
    return start, end


def extended(seq: bytes, flank: int) -> bytes:
    """`flank` bases of the circle's end, the contig, and `flank` bases of its start (flank <= len(seq))."""
    return seq if flank == 0 else seq[len(seq) - flank:] + seq + seq[:flank]


# --- one batch ---------------------------------------------------------------------------------------------------


class _Batch:
    """One FASTA scratch file of records, and what each record is."""

    def __init__(self, path: Path):
        self.path = path
        self.fh = open(path, "wb")
        self.recs: dict[str, Rec] = {}
        self.seqs: dict[str, tuple[int, int]] = {}  # record -> (byte offset of its sequence line, length)
        self.bases = 0
        self.at = 0

    def add(self, name: str, rec: bytes, what: Rec) -> None:
        if name in self.recs:
            return  # a contig named twice in the input is searched once
        head = b">" + name.encode() + b"\n"
        self.fh.write(head + rec + b"\n")
        self.seqs[name] = (self.at + len(head), len(rec))
        self.at += len(head) + len(rec) + 1
        self.recs[name] = what
        self.bases += len(rec)

    def close(self) -> None:
        self.fh.close()


def _search(java: str, jar: Path, batch: _Batch, scratch: Path) -> tuple[pa.RecordBatch, pa.RecordBatch,
                                                                          pa.RecordBatch]:
    """MinCED on one batch file: its arrays (feature rows), their elements (`l1_crispr` rows) and per contig with an
    array the arrays in and rejected. Deletes the batch's files."""
    gff, txt = batch.path.with_suffix(".gff"), batch.path.with_suffix(".txt")
    try:
        done = subprocess.run(_argv(java, jar, batch.path, txt, gff, scratch), stdin=subprocess.DEVNULL,
                              stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, cwd=scratch)
        if done.returncode != 0:
            raise RuntimeError(f"minced exited with code {done.returncode}: "
                               f"{done.stderr.decode(errors='replace')[-2000:]}")
        with open(gff) as fh:
            arrays = list(parse_gff(fh))
        feats: list[dict] = []
        rows: list[dict] = []
        n_in: dict[bytes, int] = {}
        n_rej: dict[bytes, int] = {}
        with open(batch.path, "rb") as fna:
            for arr in arrays:
                rec = batch.recs.get(arr.record)
                if rec is None:
                    raise ValueError(f"minced reported an array on record {arr.record!r}, not in the batch")
                where = place(arr, rec)
                if where is None:
                    continue  # the copy, in a flank, of an array kept where its leftmost base lies on the contig
                n_in[rec.contig] = n_in.get(rec.contig, 0) + 1
                if where == REJECT_CIRCLE:
                    n_rej[rec.contig] = n_rej.get(rec.contig, 0) + 1
                    continue
                start, end = where
                off, length = batch.seqs[arr.record]
                if arr.end > length:
                    raise ValueError(f"minced reported [{arr.start}, {arr.end}) on a record of {length} bases")
                fna.seek(off + arr.start)
                bases = fna.read(arr.end - arr.start)
                fu = ids.feature_uid(rec.contig, start, end, STRAND, KIND_CRISPR)
                feats.append({"feature_uid": fu, "contig_uid": rec.contig, "start": start, "end": end,
                              "strand": STRAND})
                for kind, ordinal, a, b in elements(arr.units):
                    s, e = a - rec.flank, b - rec.flank
                    if rec.circular and s >= rec.n:  # wholly past the origin: its position on the contig
                        s, e = s - rec.n, e - rec.n
                    rows.append({"feature_uid": fu, "contig_uid": rec.contig, "region_uid": ids.region_uid(
                        rec.contig, s, e, STRAND), "element": kind, "ordinal": ordinal, "start": s, "end": e,
                        "strand": STRAND, "seq": bases[a - arr.start:b - arr.start].decode("ascii"),
                        "contig_len": rec.n})
        keys = sorted(n_in)
        contigs = pa.RecordBatch.from_arrays(
            [pa.array(keys, type=pa.binary(16)), pa.array([n_in[k] for k in keys], type=pa.int64()),
             pa.array([n_rej.get(k, 0) for k in keys], type=pa.int64()), pa.array([None] * len(keys), pa.string())],
            schema=_CONTIG_SCRATCH)
        return (pa.RecordBatch.from_pylist(feats, schema=_ARRAY_SCRATCH),
                pa.RecordBatch.from_pylist(rows, schema=_ELEMENT_SCRATCH), contigs)
    finally:
        for f in (batch.path, gff, txt):
            f.unlink(missing_ok=True)


class _Scan:
    """The scan of a run of contigs into unsorted scratch Parquet files: arrays, elements and, per contig, arrays in
    and rejected (a row of zeros for each contig searched, the counts per batch, the reason for a contig another stage
    owns). Up to `cpus` MinCED processes at once; only this thread writes the scratch files."""

    def __init__(self, java: str, jar: Path, cpus: int, scratch: Path, batch_bases: int):
        self.java, self.jar, self.cpus, self.scratch, self.batch_bases = java, jar, cpus, scratch, batch_bases
        self.paths = {k: scratch / f"{k}.parquet" for k in ("arrays", "elements", "contigs")}
        self.writers = {k: pq.ParquetWriter(self.paths[k], s, compression="zstd", write_statistics=False)
                        for k, s in (("arrays", _ARRAY_SCRATCH), ("elements", _ELEMENT_SCRATCH),
                                     ("contigs", _CONTIG_SCRATCH))}
        self.batch: _Batch | None = None
        self.pool = ThreadPoolExecutor(max_workers=cpus)
        self.running: set[Future] = set()
        self.stats = {"contigs": 0, "bases": 0, "batches": 0, "arrays": 0, "repeats": 0, "spacers": 0}
        self.shared: dict[bytes, str] = {}  # contig in several genomes of the stage -> its owner (DEC-27)
        self.k = 0

    def _collect(self, fut: Future) -> None:
        arrays, elements_, contigs = fut.result()
        self.writers["arrays"].write_batch(arrays)
        self.writers["elements"].write_batch(elements_)
        self.writers["contigs"].write_batch(contigs)
        self.stats["batches"] += 1
        self.stats["arrays"] += arrays.num_rows
        kinds = elements_.column("element").to_pylist()
        self.stats["repeats"] += kinds.count("repeat")
        self.stats["spacers"] += kinds.count("spacer")

    def _drain(self, keep: int) -> None:
        """Collect finished batches until at most `keep` run."""
        while len(self.running) > keep:
            done, _ = wait(self.running, return_when=FIRST_COMPLETED)
            for fut in done:  # in any order: the scratch rows are sorted before they are published
                self.running.discard(fut)
                self._collect(fut)

    def _flush(self) -> None:
        b, self.batch = self.batch, None
        if b is None:
            return
        b.close()
        self._drain(self.cpus - 1)
        self.running.add(self.pool.submit(_search, self.java, self.jar, b, self.scratch))

    def skip(self, cu: bytes, owner: str) -> None:
        self.writers["contigs"].write_batch(pa.RecordBatch.from_arrays(
            [pa.array([cu], pa.binary(16)), pa.array([0], pa.int64()), pa.array([0], pa.int64()),
             pa.array([f"{SKIP_REASON}{owner}"], pa.string())], schema=_CONTIG_SCRATCH))

    def contig(self, cu: bytes, seq: bytes, circular: bool) -> None:
        n = len(seq)
        self.stats["contigs"] += 1
        self.stats["bases"] += n
        self.writers["contigs"].write_batch(pa.RecordBatch.from_arrays(
            [pa.array([cu], pa.binary(16)), pa.array([0], pa.int64()), pa.array([0], pa.int64()),
             pa.array([None], pa.string())], schema=_CONTIG_SCRATCH))
        if n == 0:
            return
        f = min(CIRC_FLANK, n) if circular else 0
        rec = extended(seq, f)
        if self.batch is None:
            self.batch = _Batch(self.scratch / f"batch-{self.k:06d}.fna")
            self.k += 1
        self.batch.add(cu.hex(), rec, Rec(cu, n, f, circular, len(rec)))
        if self.batch.bases >= self.batch_bases:
            self._flush()

    def close(self, ok: bool) -> None:
        try:
            if ok:
                self._flush()
                self._drain(0)
        finally:
            for fut in self.running:
                fut.cancel()
            self.pool.shutdown(wait=True)
            if self.batch is not None:
                self.batch.close()
            for w in self.writers.values():
                w.close()


def _scan(root: Path, stage: str, contigs: Iterator[tuple[bytes, bool, bytes | None]], stages: list[str], cpus: int,
          scratch: Path, batch_bases: int, java: str, jar: Path) -> tuple[_Scan, dict[bytes, str]]:
    """Search `contigs` ((contig_uid, circular, owner genome), sorted by contig_uid) into scratch; a contig another
    stage owns is not searched (DEC-26) and is returned with its owner. A contig in several genomes of this stage
    (DEC-27) is searched once, for its owner genome, and `sc.shared` holds it with that owner (`<stage>/<genome_uid>`),
    which its contig row carries as the reason of the other genomes' skip."""
    store = SeqStore(root, stages)
    owners = ContigOwners(root, stage, stages)
    skipped: dict[bytes, str] = {}
    sc = _Scan(java, jar, cpus, scratch, batch_bases)
    ok = False
    try:
        for cu, circ, genome in contigs:
            s = owners.owner(cu)
            if s != stage:
                skipped[cu] = s
                sc.skip(cu, s)
                continue
            if genome is not None:
                sc.shared[cu] = f"{stage}{OWNER_SEP}{genome.hex()}"
                sc.skip(cu, sc.shared[cu])
            sc.contig(cu, store.get(cu).upper(), circ)
        ok = True
    finally:
        sc.close(ok)
        store.close()
    return sc, skipped


def _paths(stage: str, prod: str) -> tuple[str, str, str]:
    return (f"l1/feature/producer={prod}/stage={stage}", f"{table('l1_crispr').location}/producer={prod}/stage={stage}",
            f"ledger/stage={stage}/step=crispr-{prod}")


def _check_args(cpus: int, batch_bases: int) -> None:
    if not isinstance(cpus, int) or isinstance(cpus, bool) or cpus < 1:
        raise ValueError(f"cpus {cpus!r}: give the number of MinCED processes at once explicitly (at least 1)")
    if not isinstance(batch_bases, int) or batch_bases < 1:
        raise ValueError(f"batch_bases must be a positive integer, got {batch_bases!r}")


def run(root: Path, stage: str, cpus: int, genomes: Iterable[GenomeInput] | None = None, *, tmp=None,
        out: Path | None = None, seq_stages: list[str] | None = None, batch_bases: int = BATCH_BASES,
        max_file_bytes: int = layout.MAX_FILE_BYTES) -> dict:
    """Call the CRISPR arrays of the stage's contigs with MinCED and write its features, elements and ledger (module
    docstring).

    `cpus`: MinCED processes at once, given explicitly; the output bytes do not depend on it. `genomes`: any iterable,
    consumed once; None reads the stage's `genome_input` table (which also says which contigs are circular). Contigs
    are searched in contig_uid order. `seq_stages`: the stages whose sequences the genomes may name (card B-28),
    default the stage itself; a contig another stage stored first (DEC-26) is skipped and listed in the result
    (`skipped_contigs`) and the ledger, and one in several genomes of the stage is searched once (module docstring).
    `tmp` holds the private scratch directory (system temp by default) and may lie under neither the root nor `out`.
    `batch_bases`: bases per MinCED process; the output does not depend on it. `max_file_bytes`: each table is cut
    into parts of at most this many estimated bytes.

    Inputs are read from `root`; every output and the ledger are written under `out` (default `root`; card C-08).
    Refuses, before reading a sequence, when one of the producer's partitions of this stage already holds a file in
    `root` or `out`."""
    root = Path(root)
    out = root if out is None else Path(out)
    _check_args(cpus, batch_bases)
    ledger.check_scratch(tmp, root, out)
    java, jar = _tool()
    prod = producer()
    ledger.refuse_written(root, out, list(_paths(stage, prod)), f"CRISPR arrays {prod} of stage {stage}")
    tools = _tools(java, jar)
    stages = list(dict.fromkeys([stage, *([] if seq_stages is None else seq_stages)]))
    scratch = _mkscratch("biobox-crispr-", tmp)
    try:
        if genomes is None:
            gi = _genome_input_files(root, stage)
            if not gi:
                raise FileNotFoundError(f"stage {stage}: no genome_input")
            contigs = _distinct(_range_rows(gi, ["contig_uid", "is_circular", "genome_uid"]))
        else:
            contigs = _distinct(iter(sorted({(cu, cu in set(g.circular), ids.genome_uid(list(g.contig_uids)))
                                             for g in genomes for cu in g.contig_uids})))
        sc, skipped = _scan(root, stage, contigs, stages, cpus, scratch, batch_bases, java, jar)
        rejected = _publish(out, stage, prod, sc, scratch, tools, max_file_bytes)
        return {**sc.stats, "producer": prod, "rejected": rejected, "minced_args": tools["minced_args"],
                "java": tools["java"],
                "skipped_contigs": [skipped_entry(cu, s) for cu, s in sorted({**skipped, **sc.shared}.items())]}
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def _sql_path(p: Path) -> str:
    return "'" + str(p).replace("'", "''") + "'"


def _publish(out: Path, stage: str, prod: str, sc: _Scan, scratch: Path, tools: dict[str, str],
             max_file_bytes: int) -> int:
    """Sort the scratch rows into the producer's `feature` and `l1_crispr` files and write the ledger: per contig
    searched, arrays in, kept and rejected (longer than a circular contig); a contig owned by another stage (DEC-26)
    has one row, reason `owned_by:<stage>`, counts 0; one in several genomes of the stage (DEC-27) has a row, reason
    `owned_by:<stage>/<genome_uid>`, after its `ok` row and with its counts. Returns the arrays rejected."""
    feat_dir, el_dir, _ = _paths(stage, prod)
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb"), "preserve_insertion_order": True})
    try:
        con.execute(f"CREATE VIEW arrays AS SELECT DISTINCT * FROM read_parquet({_sql_path(sc.paths['arrays'])})")
        con.execute(f"CREATE VIEW elements AS SELECT DISTINCT * FROM read_parquet({_sql_path(sc.paths['elements'])})")
        con.execute("CREATE VIEW contigs AS SELECT contig_uid, sum(n_in) AS n_in, sum(n_rejected) AS n_rejected, "
                    f"max(reason) AS reason FROM read_parquet({_sql_files([sc.paths['contigs']])}) GROUP BY contig_uid")
        spec = table("feature")
        io.write_parts(io.duck_batches(
            con, f'SELECT feature_uid, contig_uid, start, "end", strand, {KIND_CRISPR}::UTINYINT AS kind FROM arrays '
                 'ORDER BY contig_uid, start, "end", strand', spec),
            spec, out / feat_dir, EST_FEATURE_ROW, max_file_bytes)
        spec = table("l1_crispr")
        cols = ", ".join(f'"{c}"' for c in spec.schema.names)
        io.write_parts(io.duck_batches(
            con, f"SELECT {cols} FROM elements ORDER BY contig_uid, feature_uid, ordinal", spec),
            spec, out / el_dir, EST_ELEMENT_ROW, max_file_bytes)
        rejected = con.execute("SELECT coalesce(sum(n_rejected), 0) FROM contigs").fetchall()[0][0]

        def ledger_rows() -> Iterator[dict]:
            reader = con.sql("SELECT contig_uid, n_in, n_rejected, reason FROM contigs ORDER BY contig_uid"
                             ).to_arrow_reader(io.ROW_GROUP_ROWS)
            try:
                for b in reader:
                    for cu, n, rej, why in zip(*(b.column(i).to_pylist() for i in range(4))):
                        if why is not None and OWNER_SEP not in why.removeprefix(SKIP_REASON):
                            yield ledger.row(cu, "contig", 0, 0, reason=why, tools=tools)  # another stage's
                            continue
                        kept = n - rej
                        if rej:
                            yield ledger.row(cu, "contig", n, kept, n_rejected=rej, reason=REJECT_CIRCLE, tools=tools)
                        yield ledger.row(cu, "contig", n, kept, tools=tools)
                        if why is not None:  # searched once for the genome that owns it (DEC-27); its reason
                            # sorts after `ok` and repeats the contig's counts (L-01)
                            yield ledger.row(cu, "contig", n, kept, reason=why, tools=tools)
            finally:
                reader.close()

        ledger.write(out, stage, f"crispr-{prod}", ledger_rows())
        return int(rejected)
    finally:
        con.close()
