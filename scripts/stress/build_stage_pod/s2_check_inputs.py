"""Dry check of a stage's `build-stage` inputs with the `inputs` and `ingest` steps' own readers; writes nothing.

    python s2_check_inputs.py [--tables-only] STAGE GENOMES_JSON OBJECTS_JSON

`bio-box stage inputs` has no dry-run flag: it parses `--genomes` with `biobox.cli._genomes` (JSON list
to `GenomeInput`), then `l0seq.inputs.write_genome_input` refuses a partition already written, checks
the genomes in `_rows` (no contig twice in a genome, circular contigs inside their genome, no genome
listed twice, and a contig in several genomes only when they agree on its kind, genetic code and
topology: DEC-27, card C-21a) and writes. This runs the parsing and `_rows` only, never `write_genome_input`, so no
root is read or written. It also parses `--objects` as `stage ingest` does (`cli._source_objects`),
checks each object path is a file, checks the stage name (`catalog.check_stages`), each genome's kind
and assembly quality against the `genome_input` schema's enums, and asks `genecall._finder` which
gene-calling mode the first genome of each kind gets. Prints a JSON report; exit 0 when every check
passes, else 1; 2 on a usage error.

The byte pass (card R-21; skipped with `--tables-only`, for a re-run after a change of metadata only). The checks
above compare tables; `LR2` passed them twice with 691 gzip files written into its FASTA as sequence (card R-19).
So each object is read once, streaming, with card R-19's reader (`r19_fasta_scan._Scan`: a line piece of at most
1 MiB, records split as `biobox.fasta.iter_records` splits them), and each record's `contig_uid` is hashed as it is
read (`biobox.ids.contig_uid`: blake3-128 of the upper-cased sequence, here over buffered pieces). Named problems,
one line per object and kind, with the count and the first record: `header_not_utf8` (ingest decodes every
header); `gap_byte` (`-` or `.` in a sequence); `non_iupac_byte` (any other byte outside `ACGTURYSWKMBDHVN`, either
case); `empty_record`; `record_not_in_genomes` (its `contig_uid` is in no genome of GENOMES_JSON);
`sequence_before_first_header`, `header_too_long` (over 1 MiB), `read_error` and `compression_not_as_named` (the
first two bytes say gzip and the name does not end in `.gz`, or the reverse: ingest opens by the name), so that the
check agrees with ingest; and `contigs_in_no_object` once for the stage. A UTF-8 header that is not ASCII is counted,
not a problem (ingest accepts it). Memory: the stage's contig uids as a dict over the parsed genomes' own uid objects,
and one line piece and a 1 MiB hash buffer; no sequence is held.
"""

import gzip
import io
import json
import resource
import sys
import time
import zlib
from pathlib import Path

from blake3 import blake3

from biobox import catalog, cli
from biobox.l0seq import inputs
from biobox.l1 import genecall
from biobox.registry.schemas import table

_HERE = str(Path(__file__).resolve().parent)
if _HERE not in sys.path:  # absent under `python -I` or when imported from elsewhere
    sys.path.append(_HERE)
import r19_fasta_scan as r19  # noqa: E402  (card R-19's reader, in this directory)

USAGE = __doc__.split("\n\n")[1]
GAPS = b"-."
HASH_BUF = 1 << 20


class _UidScan(r19._Scan):
    """Card R-19's scan, with each record's contig uid hashed as its pieces are read and looked up in `want`."""

    def __init__(self, want: dict) -> None:
        super().__init__()
        self.want = want
        self.buf = bytearray()
        self.hasher = blake3()
        self.rec_gap = 0
        self.rec_other = 0
        self.gap_bytes = {"-": 0, ".": 0}
        self.records_with_gap = 0
        self.records_with_other = 0
        self.other_bytes = 0
        self.not_in_genomes = 0
        self.first: dict[str, tuple[int, str]] = {}  # problem name -> (record, header shown)

    def open_record(self, first_piece: bytes) -> None:
        super().open_record(first_piece)  # closes the record before it through close_record
        self.buf = bytearray()
        self.hasher = blake3()
        self.rec_gap = 0
        self.rec_other = 0

    def seq_piece(self, piece: bytes) -> None:
        before = self.bad_bytes
        super().seq_piece(piece)
        if self.bad_bytes != before:
            dash, dot = piece.count(b"-"), piece.count(b".")
            self.gap_bytes["-"] += dash
            self.gap_bytes["."] += dot
            self.rec_gap += dash + dot
            self.rec_other += self.bad_bytes - before - dash - dot
        self.buf += piece
        if len(self.buf) >= HASH_BUF:
            self.hasher.update(self.buf.upper())
            self.buf = bytearray()

    def _first(self, name: str, n: int, hdr: bytes) -> None:
        if name not in self.first:
            self.first[name] = (n, r19._show(hdr))

    def close_record(self) -> None:
        if self.records == 0 or self.header is None:
            return
        n, hdr = self.records, bytes(self.header)
        if self.buf:
            self.hasher.update(self.buf.upper())
            self.buf = bytearray()
        uid = self.hasher.digest(length=16)
        if uid in self.want:
            self.want[uid] = True
        else:
            self.not_in_genomes += 1
            self._first("record_not_in_genomes", n, hdr)
        if self.rec_gap:
            self.records_with_gap += 1
            self._first("gap_byte", n, hdr)
        if self.rec_other:
            self.records_with_other += 1
            self.other_bytes += self.rec_other
            self._first("non_iupac_byte", n, hdr)
        super().close_record()
        for reason in ("header_not_utf8", "empty_record", "header_too_long"):  # R-19's tests, its first record
            if reason in self.first_record_of:
                self._first(reason, n, hdr)


def _scan_object(path: Path, want: dict) -> tuple[dict, list[str]]:
    t0, r0 = time.monotonic(), resource.getrusage(resource.RUSAGE_SELF)
    sc = _UidScan(want)
    comp, read_error = None, None
    try:
        with open(path, "rb", buffering=0) as raw:
            buffered = io.BufferedReader(raw, buffer_size=r19.PIECE)
            comp = r19._compression(buffered)
            try:
                if comp == "gzip":
                    with gzip.GzipFile(fileobj=buffered, mode="rb") as gz:
                        sc.feed(gz)
                else:
                    sc.feed(buffered)
            except (OSError, EOFError, zlib.error) as e:  # as r19_fasta_scan.scan_file catches them
                read_error = f"{type(e).__name__}: {e}"
                sc.close_record()
    except OSError as e:
        read_error = f"{type(e).__name__}: {e}"
    r1 = resource.getrusage(resource.RUSAGE_SELF)
    c = sc.counts
    rep = {"path": str(path), "compression": comp, "records": sc.records, "bases": sc.bases,
           "headers_not_utf8": c["headers_not_utf8"],
           # UTF-8 but not printable ASCII (a header that is not UTF-8 is never printable ASCII): accepted by ingest
           "headers_not_ascii": c["headers_not_printable_ascii"] - c["headers_not_utf8"],
           "records_with_gap_bytes": sc.records_with_gap, "gap_bytes": sc.gap_bytes,
           "records_with_other_non_iupac": sc.records_with_other, "other_non_iupac_bytes": sc.other_bytes,
           "other_non_iupac_byte_values": {f"0x{b:02x}": k for b, k in sorted(
               ((b, k) for b, k in sc.bad_values.items() if b not in GAPS), key=lambda x: (-x[1], x[0]))[:10]},
           "empty_records": c["empty_records"], "lines_before_first_header": c["lines_before_first_header"],
           "headers_too_long": c["headers_too_long"], "records_not_in_genomes": sc.not_in_genomes,
           "read_error": read_error, "wall_s": round(time.monotonic() - t0, 3),
           "cpu_s": round(r1.ru_utime - r0.ru_utime + r1.ru_stime - r0.ru_stime, 3)}

    def at(name: str) -> str:
        k, h = sc.first[name]
        return f"first record {k} ({h!r})"

    probs = []
    named_gz = str(path).endswith(".gz")  # biobox.fasta._open
    if comp is not None and (comp == "gzip") != named_gz:
        probs.append(f"object {path}: compression_not_as_named: its first two bytes say {comp} and its name "
                     f"{'ends' if named_gz else 'does not end'} in .gz; ingest opens it by the name")
    if read_error:
        probs.append(f"object {path}: read_error: {read_error}")
    if c["lines_before_first_header"]:
        probs.append(f"object {path}: sequence_before_first_header: {c['lines_before_first_header']} line(s); "
                     "ingest refuses the file")
    for name, n, what in (
        ("header_not_utf8", c["headers_not_utf8"], "record(s) whose header is not UTF-8 (ingest stops at the first)"),
        ("header_too_long", c["headers_too_long"], "record(s) whose header is over 1 MiB"),
        ("gap_byte", sc.records_with_gap,
         f"record(s) hold {sc.gap_bytes['-']} '-' and {sc.gap_bytes['.']} '.' bytes"),
        ("non_iupac_byte", sc.records_with_other,
         f"record(s) hold {sc.other_bytes} byte(s) outside the IUPAC letters other than '-' and '.' "
         f"(most frequent: {rep['other_non_iupac_byte_values']})"),
        ("empty_record", c["empty_records"], "record(s) with no sequence"),
        ("record_not_in_genomes", sc.not_in_genomes, "record(s) whose contig_uid is in no genome of genomes.json"),
    ):
        if n:
            probs.append(f"object {path}: {name}: {n} {what}, {at(name)}")
    if sc.records == 0 and not read_error:
        probs.append(f"object {path}: no_records: the file holds no FASTA record")
    return rep, probs


def _peak_rss_bytes() -> int:
    try:
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith("VmHWM:"):
                return int(line.split()[1]) * 1024
    except OSError:
        pass
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss if sys.platform == "darwin" else rss * 1024


def byte_pass(genomes, objects) -> tuple[dict, list[str]]:
    """Read every object once (see the module docstring); the report and the named problems."""
    t0, r0 = time.monotonic(), resource.getrusage(resource.RUSAGE_SELF)
    want = {u: False for g in genomes for u in g.contig_uids}
    reps, problems = [], []
    for o in objects:
        if not Path(o.path).is_file():
            continue  # named by the table checks
        rep, probs = _scan_object(Path(o.path), want)
        reps.append(rep)
        problems += probs
    missing = [u for u, seen in want.items() if not seen]
    if missing:
        problems.append(f"contigs_in_no_object: {len(missing)} contig(s) of genomes.json are in no object "
                        f"(first: {missing[0].hex()})")
    r1 = resource.getrusage(resource.RUSAGE_SELF)
    return {"objects": reps, "records": sum(r["records"] for r in reps), "bases": sum(r["bases"] for r in reps),
            "contigs": len(want), "contigs_in_no_object": len(missing),
            "wall_s": round(time.monotonic() - t0, 3),
            "cpu_s": round(r1.ru_utime - r0.ru_utime + r1.ru_stime - r0.ru_stime, 3),
            "peak_rss_bytes": _peak_rss_bytes()}, problems


def main(argv: list[str]) -> int:
    tables_only = "--tables-only" in argv
    argv = [a for a in argv if a != "--tables-only"]
    if len(argv) != 3 or any(a.startswith("--") for a in argv):
        print(USAGE, file=sys.stderr)
        return 2
    stage, gpath, opath = argv
    problems = []
    catalog.check_stages([stage])
    genomes = cli._genomes(gpath)
    rows = inputs._rows(genomes)
    objects = cli._source_objects(opath)
    enums = table("genome_input").enums
    for g in genomes:
        for col, val in (("kind", g.kind), ("assembly_quality", g.assembly_quality)):
            if col in enums and val not in enums[col]:
                problems.append(f"genome {g.contig_uids[0].hex()}: {col} {val!r} not in {enums[col]}")
    for o in objects:
        if not Path(o.path).is_file():
            problems.append(f"object {o.path} is not a file")
    modes = {}
    for g in genomes:
        if g.kind not in modes:
            modes[g.kind] = genecall._finder(g, {u: "" for u in g.contig_uids})[1] if g.kind == "virus" else None
    bp = None
    if not tables_only:
        bp, bprobs = byte_pass(genomes, objects)
        problems += bprobs
    rep = {"stage": stage, "genomes": len(genomes), "genome_input_rows": len(rows),
           "distinct_genome_uids": len({r[1] for r in rows}), "circular_rows": sum(1 for r in rows if r[2]),
           "kinds": sorted({g.kind for g in genomes}), "codes": sorted({g.genetic_code for g in genomes}),
           "assembly_quality": sorted({g.assembly_quality for g in genomes}), "gene_call_mode": modes,
           "objects": [{"path": str(o.path), "source_name": o.source_name, "source_release": o.source_release,
                        "license": o.license, "redistributable": o.redistributable,
                        "retrieved_at": o.retrieved_at.isoformat(), "bytes": Path(o.path).stat().st_size}
                       for o in objects],
           "biobox": cli.__file__, "tables_only": tables_only, "byte_pass": bp, "problems": problems}
    print(json.dumps(rep, indent=1))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
