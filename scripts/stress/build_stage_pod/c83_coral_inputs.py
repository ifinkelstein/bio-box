"""Card C-83 (DEC-85): stage the two coral-reef sets on `tower` for `r26.11.5`: the Great Barrier Reef genomes
(`GBR_MGD`, stage `GBRMGD`) and Tara Pacific's Reef Microbiomics Database genomes (`RMD1`).

    python c83_coral_inputs.py census TAR OUT_TSV [--report-out F]

`census` reads an archive once, streaming (`tarfile` mode `r|*`: the outer gzip or bzip2 is taken from the file's
bytes), and writes one row per member in archive order: its size in the tar, its own compression told by its first
bytes (gzip `1f 8b`, bzip2 `BZh`, xz `fd 37 7a 58 5a 00`, zstd `28 b5 2f fd`; never by its name), and its content read
from the decompressed bytes by `parse_fasta`: blocks of 4 MiB cut at their last newline, a line that starts with `>`
opens a record (its header without the trailing CR), and the sequence between two headers is its bytes without CR and
LF, counted in C (`bytes.translate`), never line by line (card R-19's line reader ran at about 5 MB/s here). A
space or tab in a sequence line is counted as an other byte; `s2_check_inputs.py` (card R-21) reads the staged output
with R-19's own reader. The kind of a member:

- `dna`: FASTA whose sequence bytes are all IUPAC nucleotide letters (either case);
- `dna_with_gaps`: FASTA of nucleotides holding `-` or `.` (an alignment; ingest refuses a gap byte, DEC-84 (1));
- `dna_with_other_bytes`: FASTA, mostly nucleotides, holding another byte (refused by ingest likewise);
- `protein`: FASTA whose sequence is less than half `ACGTUN` (amino acids);
- `empty`: no byte; `tar`: an archive inside the archive; `text`: UTF-8 that is not FASTA (a table, a README);
  `binary`: anything else; `zstd`: a zstd member (named, not read: the standard library has no reader).

Memory: counters and the current header only; no sequence is held.
"""

from __future__ import annotations

import argparse
import bz2
import csv
import gzip
import io
import json
import lzma
import re
import shutil
import subprocess
import sys
import tarfile
import zlib
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import a26_r26_11_4_chain as a26  # noqa: E402
import c39_ani  # noqa: E402
import c58_sags_symbionts_niche as c58  # noqa: E402
import c77_sag2_inputs as c77  # noqa: E402
import c81_gwv1_sag3_wet1 as c81  # noqa: E402
import make_inputs  # noqa: E402
import r232_index  # noqa: E402
from biobox import fasta, ids  # noqa: E402

Refused = make_inputs.Refused

PIECE = 1 << 20
BLOCK = 4 << 20
MAX_LINE = 256 << 20  # a header line longer than this is refused
MAGIC = ((b"\x1f\x8b", "gzip"), (b"BZh", "bzip2"), (b"\xfd7zXZ\x00", "xz"), (b"\x28\xb5\x2f\xfd", "zstd"))
ACGTN = b"ACGTUNacgtun"
IUPAC = b"ACGTURYSWKMBDHVNacgturyswkmbdhvn"
GAPS = b"-."
CRLF = b"\r\n"
SHOW = 200
TASKS_PER_CHILD = 500
SUFFIXES = (".gz", ".bz2", ".xz", ".zst", ".fasta", ".fna", ".fa", ".fas", ".ffn", ".faa", ".fsa", ".seq", ".txt")
CENSUS_HEAD = ["member", "tar_bytes", "compression", "kind", "records", "bases", "acgtn_bases", "gap_bytes",
               "other_bytes", "headers_not_utf8", "empty_records", "min_record", "max_record", "first_header",
               "stem"]


# === reading a member by its bytes ==================================================================================

def compression_of(head: bytes) -> str:
    for magic, name in MAGIC:
        if head.startswith(magic):
            return name
    return "none"


def open_member(fh) -> tuple[str, io.BufferedIOBase | None]:
    """(compression, decompressed stream) of a binary file object, the compression told by its first bytes; for
    `zstd` the stream is None (named, not read)."""
    buf = fh if isinstance(fh, io.BufferedReader) else io.BufferedReader(_Raw(fh), buffer_size=PIECE)
    comp = compression_of(buf.peek(6)[:6])
    if comp == "gzip":
        return comp, gzip.GzipFile(fileobj=buf, mode="rb")
    if comp == "bzip2":
        return comp, bz2.BZ2File(buf, mode="rb")
    if comp == "xz":
        return comp, lzma.LZMAFile(buf, mode="rb")
    if comp == "zstd":
        return comp, None
    return comp, buf


class _Raw(io.RawIOBase):
    """Any object with `read` as a raw stream, so that it can be buffered and peeked."""

    def __init__(self, fh):
        self._fh = fh

    def readable(self) -> bool:
        return True

    def readinto(self, b) -> int:
        data = self._fh.read(len(b))
        n = len(data)
        b[:n] = data
        return n


def stem_of(member: str) -> str:
    """`bins/A.fna.gz` -> `A`: the base name without the compression and FASTA suffixes it carries."""
    name = member.rsplit("/", 1)[-1]
    changed = True
    while changed:
        changed = False
        for s in SUFFIXES:
            if name.lower().endswith(s) and len(name) > len(s):
                name = name[: -len(s)]
                changed = True
    return name


def parse_fasta(fh, on_header, on_seq) -> None:
    """Read FASTA from a binary stream in blocks: `on_header(header bytes)` for each line that starts with `>`,
    `on_seq(bytes)` for the sequence bytes after it (CR and LF removed; any number of calls per record), and
    `on_seq` also for what comes before the first header (the caller tells)."""
    carry = b""
    while True:
        block = fh.read(BLOCK)
        if not block:
            break
        t = carry + block if carry else block
        cut = t.rfind(b"\n") + 1
        if cut == 0:
            carry = t
            if len(carry) > MAX_LINE:
                raise Refused(f"a line of more than {MAX_LINE} bytes")
            continue
        carry = t[cut:]
        _lines(t[:cut] if cut < len(t) else t, on_header, on_seq)
    if carry:
        _lines(carry + b"\n", on_header, on_seq)


def _lines(t: bytes, on_header, on_seq) -> None:
    """T is whole lines (it ends with LF) and starts at a line's start."""
    i, n = 0, len(t)
    while i < n:
        if t[i] == 0x3E:  # '>'
            e = t.find(b"\n", i)
            on_header(t[i + 1:e].rstrip(b"\r"))
            i = e + 1
            continue
        j = t.find(b"\n>", i)
        if j == -1:
            on_seq(t[i:].translate(None, CRLF))
            return
        on_seq(t[i:j + 1].translate(None, CRLF))
        i = j + 1


class _Census:
    """Counts over one member's FASTA: records, bases (sequence bytes), ACGTUN bases, gap bytes, other bytes (neither
    IUPAC nor a gap), headers that are not UTF-8, empty records, the least and largest record."""

    def __init__(self) -> None:
        self.records = self.bases = self.acgtn = self.gaps = self.other = 0
        self.not_utf8 = self.empty = self.before = 0
        self.min_len: int | None = None
        self.max_len = 0
        self.first_header = ""
        self.cur: int | None = None

    def header(self, h: bytes) -> None:
        self._close()
        self.records += 1
        try:
            h.decode("utf-8")
        except UnicodeDecodeError:
            self.not_utf8 += 1
        if self.records == 1:
            self.first_header = h.decode("utf-8", "backslashreplace")[:SHOW]
        self.cur = 0

    def seq(self, s: bytes) -> None:
        if self.cur is None:
            if s.strip():
                self.before += 1
            return
        n = len(s)
        self.cur += n
        self.bases += n
        self.acgtn += n - len(s.translate(None, ACGTN))
        self.gaps += s.count(b"-") + s.count(b".")
        self.other += len(s.translate(None, IUPAC + GAPS))

    def _close(self) -> None:
        if self.cur is None:
            return
        if self.cur == 0:
            self.empty += 1
        self.min_len = self.cur if self.min_len is None else min(self.min_len, self.cur)
        self.max_len = max(self.max_len, self.cur)
        self.cur = None

    def finish(self) -> None:
        self._close()


def kind_of(sc: _Census | None, head: bytes, nbytes: int) -> str:
    if nbytes == 0:
        return "empty"
    if sc is not None and sc.records > 0 and sc.before == 0:
        if sc.bases and sc.acgtn < 0.5 * sc.bases:
            return "protein"
        if sc.gaps:
            return "dna_with_gaps"
        if sc.other:
            return "dna_with_other_bytes"
        return "dna"
    if len(head) >= 262 and head[257:262] == b"ustar":
        return "tar"
    try:
        head.decode("utf-8")
        return "text"
    except UnicodeDecodeError:
        # a piece cut inside a character still counts as text
        try:
            head[:-4].decode("utf-8")
            return "text"
        except UnicodeDecodeError:
            return "binary"


class _Count(io.RawIOBase):
    def __init__(self, fh):
        self._fh = fh
        self.n = 0
        self.head = b""

    def readable(self) -> bool:
        return True

    def readinto(self, b) -> int:
        data = self._fh.read(len(b))
        k = len(data)
        b[:k] = data
        if len(self.head) < 65536:
            self.head += data[: 65536 - len(self.head)]
        self.n += k
        return k


def scan_member(fh) -> dict:
    """The census row of one member's bytes (without its name): compression, kind and counts."""
    comp, stream = open_member(fh)
    if stream is None:
        return {"compression": comp, "kind": comp, "records": 0, "bases": 0, "acgtn_bases": 0, "gap_bytes": 0,
                "other_bytes": 0, "headers_not_utf8": 0, "empty_records": 0, "min_record": "", "max_record": "",
                "first_header": ""}
    counted = _Count(stream)
    buf = io.BufferedReader(counted, buffer_size=PIECE)
    sc = None
    if buf.peek(PIECE).lstrip()[:1] == b">":
        sc = _Census()
        parse_fasta(buf, sc.header, sc.seq)
        sc.finish()
    while buf.read(BLOCK):
        pass
    kind = kind_of(sc, counted.head, counted.n)
    if sc is None or kind not in ("dna", "dna_with_gaps", "dna_with_other_bytes", "protein"):
        return {"compression": comp, "kind": kind, "records": 0, "bases": 0, "acgtn_bases": 0, "gap_bytes": 0,
                "other_bytes": 0, "headers_not_utf8": 0, "empty_records": 0, "min_record": "", "max_record": "",
                "first_header": ""}
    return {"compression": comp, "kind": kind, "records": sc.records, "bases": sc.bases, "acgtn_bases": sc.acgtn,
            "gap_bytes": sc.gaps, "other_bytes": sc.other, "headers_not_utf8": sc.not_utf8,
            "empty_records": sc.empty, "min_record": sc.min_len if sc.min_len is not None else "",
            "max_record": sc.max_len, "first_header": sc.first_header.replace("\t", " ")}


def iter_tar(tar_path: Path):
    """(TarInfo, file object or None) of each entry, in archive order, streaming; a refusal for an unreadable one."""
    try:
        with tarfile.open(tar_path, "r|*") as tf:
            for m in tf:
                yield m, (tf.extractfile(m) if m.isfile() else None)
    except (tarfile.TarError, OSError, EOFError, zlib.error, lzma.LZMAError) as exc:
        raise Refused(f"{tar_path}: not a readable archive: {exc}") from exc


def census(tar_path: Path, out_tsv: Path) -> dict:
    """Write OUT_TSV (CENSUS_HEAD, one row per member) under a temporary name and rename it; return the summary."""
    tar_path, out_tsv = Path(tar_path), Path(out_tsv)
    tmp = out_tsv.with_name(f".{out_tsv.name}.partial")
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    dirs, others = 0, Counter()
    by_kind: dict[str, dict[str, int]] = {}
    by_comp: Counter = Counter()
    n = 0
    try:
        with open(tmp, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(CENSUS_HEAD)
            for m, f in iter_tar(tar_path):
                if m.isdir():
                    dirs += 1
                    continue
                if f is None:
                    others[m.type.decode("ascii", "replace") if isinstance(m.type, bytes) else str(m.type)] += 1
                    continue
                try:
                    row = scan_member(f)
                except (OSError, EOFError, zlib.error, lzma.LZMAError) as exc:
                    raise Refused(f"{tar_path}: member {m.name} could not be read: {exc}") from exc
                n += 1
                by_comp[row["compression"]] += 1
                k = by_kind.setdefault(row["kind"], {"members": 0, "records": 0, "bases": 0})
                k["members"] += 1
                k["records"] += row["records"]
                k["bases"] += row["bases"]
                w.writerow([m.name, m.size, *(row[c] for c in CENSUS_HEAD[2:-1]), stem_of(m.name)])
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    tmp.replace(out_tsv)
    return {"archive": str(tar_path), "members": n, "directories": dirs, "other_entries": dict(sorted(others.items())),
            "by_kind": dict(sorted(by_kind.items())), "by_compression": dict(sorted(by_comp.items()))}


# === staging (both sets) ===========================================================================================

def read_census(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if rows and set(CENSUS_HEAD) - set(rows[0]):
        raise Refused(f"{path}: not a census table (no column {sorted(set(CENSUS_HEAD) - set(rows[0]))[0]})")
    return rows


def read_records(stream, name: str) -> list[tuple[str, bytes]]:
    """(header, sequence) of each record of one member; refuses a header that is not UTF-8, sequence before the first
    header, an empty record and a sequence byte that is not an IUPAC nucleotide letter (a gap byte among them)."""
    recs: list[list] = []

    def on_header(h: bytes) -> None:
        try:
            recs.append([h.decode("utf-8"), []])
        except UnicodeDecodeError as exc:
            raise Refused(f"{name}: header {h[:80]!r} is not UTF-8") from exc

    def on_seq(s: bytes) -> None:
        if not recs:
            if s.strip():
                raise Refused(f"{name}: sequence before the first header")
            return
        recs[-1][1].append(s)

    parse_fasta(stream, on_header, on_seq)
    out = []
    for hdr, parts in recs:
        seq = b"".join(parts)
        if not seq:
            raise Refused(f"{name}: record {hdr[:80]!r} is empty")
        bad = seq.translate(None, IUPAC)
        if bad:
            raise Refused(f"{name}: record {hdr[:80]!r} holds the byte 0x{bad[0]:02x}, not a nucleotide letter")
        out.append((hdr, seq))
    return out


def quality_of_mag(completeness: str, contamination: str) -> str:
    """`mag_hq` above 90 % complete and under 5 % contaminated (MIMAG's numbers, as cards C-66 and C-76), else
    `mag_mq`."""
    try:
        comp, cont = float(completeness), float(contamination)
    except ValueError:
        return "mag_mq"
    return "mag_hq" if comp > 90.0 and cont < 5.0 else "mag_mq"


def stage_set(stage: str, tar_path: Path, plan: dict[str, dict], meta_out: Path, fasta_out: Path, *,
              fasta_name: str, source: dict, retrieved_at: str, status_head: list[str], record_extra: str,
              left_out: list[tuple], root: Path | None, root_stages, owner_inputs: dict[str, Path],
              later_inputs: dict[str, Path], inputs_home: Path, work: Path | None = None, threads: int = 2,
              memory_limit: str = "4GB", workers: int = 1) -> dict:
    """Stage the members of TAR_PATH that PLAN names (member -> {genome, kind, quality, hint, extra, status}: the
    genome's id, `genome_input` kind and assembly quality, the taxonomic hint of the code screen, the value of the
    records table's last column RECORD_EXTRA, and its row of r232_status.tsv under STATUS_HEAD); every other member is
    in LEFT_OUT (member, kind, records, bases, reason), written as left_out.tsv. Whole contigs; a contig name that does
    not begin with its genome's id and `_` is given the id and `|`; the genetic code screened per genome (card C-58's
    screen, 11, 25 or 4); one BGZF FASTA in archive order, metadata in genome order. Nothing is left behind when it
    refuses."""
    tar_path = Path(tar_path)
    meta_dir, fa_dir = Path(meta_out) / stage, Path(fasta_out) / stage
    final_fa = fa_dir / fasta_name
    for p in (meta_dir, final_fa):
        if p.exists():
            raise Refused(f"{p} exists")
    if not plan:
        raise Refused(f"{stage}: nothing to stage")
    dup = [g for g, n in Counter(v["genome"] for v in plan.values()).items() if n > 1]
    if dup:
        raise Refused(f"{stage}: two members would stage as genome {sorted(dup)[0]}")
    c81._check_later(later_inputs)
    tmp = Path(meta_out) / f"{stage}.tmp"
    partial = fa_dir / f".{fasta_name}.partial"
    work_dir = Path(work) if work is not None else Path(meta_out) / f"{stage}.work"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    fa_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    try:
        rep = _stage_set(stage, tar_path, plan, tmp, partial, work_dir, final_fa, fasta_name, source, retrieved_at,
                         status_head, record_extra, left_out, root, list(root_stages), owner_inputs, later_inputs,
                         Path(inputs_home), threads, memory_limit, workers)
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        partial.unlink(missing_ok=True)
        if fa_dir.exists() and not any(fa_dir.iterdir()):
            fa_dir.rmdir()
        raise
    if work is None:
        shutil.rmtree(work_dir, ignore_errors=True)
    partial.rename(final_fa)
    tmp.rename(meta_dir)
    files = {p.name: p for p in sorted(meta_dir.iterdir())}
    files[fasta_name] = final_fa
    rep["sha256"] = {n: ids.sha256_file(p).hex() for n, p in files.items()}
    rep["bytes"] = {n: p.stat().st_size for n, p in files.items()}
    return rep


def _stage_set(stage, tar_path, plan, tmp, partial, work_dir, final_fa, fasta_name, source, retrieved_at,
               status_head, record_extra, left_out, root, root_stages, owner_inputs, later_inputs, inputs_home,
               threads, memory_limit, workers) -> dict:
    gc_results: dict[str, dict] = {}
    # each screen worker is replaced after TASKS_PER_CHILD genomes, so that no process passes the pod's 300
    # CPU-minutes however large the set (RMD1's 13,549 genomes would give each of 3 workers more, *derived* from
    # GBRMGD's stage: about 260 CPU-minutes for 5,283 bins over 4 processes)
    pool = (ProcessPoolExecutor(max_workers=workers, max_tasks_per_child=TASKS_PER_CHILD) if workers > 1 else None)
    pending: list[tuple[str, object]] = []
    by_genome = {v["genome"]: v for v in plan.values()}

    def submit(g: str, seqs: list[bytes]) -> None:
        hint = by_genome[g]["hint"]
        if pool is None:
            gc_results[g] = c58.screen_genetic_code(seqs, hint, screen_max_bp=500_000)
            return
        pending.append((g, pool.submit(c58._screen_worker, (g, seqs, hint, 500_000))))
        while len(pending) >= 4 * workers:
            s, fut = pending.pop(0)
            gc_results[s] = fut.result()[1]

    contigs: dict[str, tuple[bytes, list[int]]] = {}
    member_comp: Counter = Counter()
    names: set[str] = set()
    n_renamed = 0
    w = c77.BgzfWriter(partial)
    try:
        for m, f in iter_tar(tar_path):
            if f is None or m.name not in plan:
                continue
            g = plan[m.name]["genome"]
            comp, stream = open_member(f)
            if stream is None:
                raise Refused(f"{m.name}: a {comp} member, which this script does not read")
            member_comp[comp] += 1
            try:
                recs = read_records(stream, m.name)
            except (OSError, EOFError, zlib.error, lzma.LZMAError) as exc:
                raise Refused(f"{tar_path}: member {m.name} could not be read: {exc}") from exc
            if not recs:
                raise Refused(f"{m.name} holds no sequence")
            for hdr, seq in recs:
                first = hdr.split(None, 1)[0] if hdr.strip() else ""
                if not first.startswith(g + "_"):
                    hdr = f"{g}|{hdr}"
                    n_renamed += 1
                first = hdr.split(None, 1)[0]
                if first in names:
                    raise Refused(f"contig name {first} occurs twice in the stage")
                names.add(first)
                w.write(b"\n".join(fasta.format_record(hdr, seq, 80)) + b"\n")
            contigs[g] = (b"".join(ids.contig_uid(s) for _, s in recs), [len(s) for _, s in recs])
            submit(g, [s for _, s in recs])
        while pending:
            s, fut = pending.pop(0)
            gc_results[s] = fut.result()[1]
    finally:
        w.close()
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)
    lost = sorted(set(by_genome) - set(contigs))
    if lost:
        raise Refused(f"{len(lost)} members were planned but not read; first {lost[0]}")
    kept = sorted(contigs)

    resolved, genome_uids, conflicts = c58._resolve_within_stage_code_conflicts(kept, contigs, gc_results)
    records_path = tmp / "records.tsv"
    ordinal = 0
    with open(records_path, "w", encoding="utf-8", buffering=1 << 20) as rf, open(
            tmp / "genetic_codes.tsv", "w", encoding="utf-8", buffering=1 << 20) as gf:
        rf.write(f"ordinal\tgenome\trecord_ordinal\tcontig_uid\tlength\tassembly_quality\tgenetic_code\t"
                 f"{record_extra}\n")
        gf.write("\t".join(c58.GC_HEAD) + "\n")
        for g in kept:
            raw, lens = contigs[g]
            gc, v = gc_results[g], by_genome[g]
            gf.write("\t".join([
                g, stage, v["hint"], str(sum(lens)), str(len(lens)), v["quality"], str(resolved[g]), gc["reason"],
                str(gc["tested_alt"]).lower(),
                c58._fmt_opt(gc["cd_11"]), c58._fmt_opt(gc["mean_orf_11"]), c58._fmt_opt(gc["n_genes_11"]),
                c58._fmt_opt(gc["cd_25"]), c58._fmt_opt(gc["mean_orf_25"]), c58._fmt_opt(gc["n_genes_25"]),
                c58._fmt_opt(gc["cd_4"]), c58._fmt_opt(gc["mean_orf_4"]), c58._fmt_opt(gc["n_genes_4"]),
                c58._fmt_opt(gc["tga_per_kb"])]) + "\n")
            for k, n in enumerate(lens):
                rf.write(f"{ordinal}\t{g}\t{k}\t{raw[k * 16:(k + 1) * 16].hex()}\t{n}\t{v['quality']}\t"
                         f"{resolved[g]}\t{v['extra']}\n")
                ordinal += 1
    owners, later = c81._owners(records_path, root, root_stages, owner_inputs, later_inputs, work_dir, threads,
                                memory_limit)

    listing = make_inputs.Listing()
    genomes, repeated, shared, held_later = [], [], [], []
    first_of: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}
    ordinal = 0
    for g in kept:
        raw, lens = contigs[g]
        uids = [raw[i: i + 16] for i in range(0, len(raw), 16)]
        for k, (cu, n) in enumerate(zip(uids, lens)):
            contig_len.setdefault(cu, n)
            if cu in owners:
                shared.append((ordinal, g, k, cu.hex(), n, owners[cu]))
            elif cu in later:
                held_later.append((ordinal, g, k, cu.hex(), n, later[cu]))
            ordinal += 1
        e = listing.add(uids, resolved[g], by_genome[g]["kind"], by_genome[g]["quality"])
        if e is None:
            repeated.append((ordinal - len(uids), g, genome_uids[g].hex(), first_of[genome_uids[g]]))
            continue
        genomes.append(e)
        first_of[genome_uids[g]] = g
    c81._write_tsv(tmp / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
    c81._write_tsv(tmp / "shared_records.tsv",
                   ["ordinal", "genome", "record_ordinal", "contig_uid", "length", "owner_stage"], shared)
    c81._write_tsv(tmp / "later_records.tsv",
                   ["ordinal", "genome", "record_ordinal", "contig_uid", "length", "holder"], held_later)
    c81._write_tsv(tmp / "left_out.tsv", ["member", "kind", "records", "bases", "reason"], left_out)
    c81._write_tsv(tmp / "r232_status.tsv", status_head,
                   [(*by_genome[g]["status"], len(contigs[g][1]), sum(contigs[g][1])) for g in kept])

    obj = {**source, "retrieved_at": retrieved_at}
    (tmp / "objects.json").write_text(json.dumps([{"path": str(inputs_home / stage / fasta_name), **obj}], indent=1))
    (tmp / "objects_staging.json").write_text(json.dumps([{"path": str(final_fa), **obj}], indent=1))
    (tmp / "genomes.json").write_text(json.dumps(genomes))

    own_c, own_b = c81._tally(contig_len, owners)
    lat_c, lat_b = c81._tally(contig_len, later, owners)
    bp = {g: sum(contigs[g][1]) for g in kept}
    status_col = status_head.index("r232_status")
    status = {g: by_genome[g]["status"][status_col] for g in kept}
    return {
        "genomes_read": len(kept), "contigs": sum(len(contigs[g][1]) for g in kept), "bases": sum(bp.values()),
        "contigs_renamed": n_renamed, "member_compression": dict(sorted(member_comp.items())),
        "by_kind": dict(sorted(Counter(by_genome[g]["kind"] for g in kept).items())),
        "by_quality": dict(sorted(Counter(by_genome[g]["quality"] for g in kept).items())),
        "bases_by_quality": {q: sum(bp[g] for g in kept if by_genome[g]["quality"] == q)
                             for q in sorted({by_genome[g]["quality"] for g in kept})},
        "r232_status": dict(sorted(Counter(status.values()).items())),
        "r232_status_bases": {k: sum(bp[g] for g in kept if status[g] == k) for k in sorted(set(status.values()))},
        "genetic_codes": dict(sorted(Counter(str(resolved[g]) for g in kept).items())),
        "code_conflicts": len(conflicts),
        "left_out": len(left_out),
        "left_out_by_reason": dict(sorted(Counter(r[-1] for r in left_out).items())),
        "genomes": len(genomes), "genomes_repeated": len(repeated), "contigs_distinct": len(contig_len),
        "bases_distinct": sum(contig_len.values()),
        "held_by": own_c, "held_bases_by": own_b, "later_held_by": lat_c, "later_held_bases_by": lat_b,
        "contigs_new": sum(1 for cu in contig_len if cu not in owners),
        "bases_new": sum(n for cu, n in contig_len.items() if cu not in owners),
        "contigs_new_after_later": sum(1 for cu in contig_len if cu not in owners and cu not in later),
        "bases_new_after_later": sum(n for cu, n in contig_len.items() if cu not in owners and cu not in later),
        **{"listing_" + k: v for k, v in listing.report().items()},
    }


LEAVE_REASON = {
    "protein": "protein FASTA (not genome sequence)",
    "dna_with_gaps": "an alignment (gap bytes; not a genome)",
    "dna_with_other_bytes": "FASTA holding bytes that are not nucleotide letters (refused by ingest, DEC-84 (1))",
    "text": "a table or text (not sequence)",
    "binary": "not sequence",
    "empty": "an empty member",
    "tar": "an archive inside the archive (not read)",
    "zstd": "a zstd member (not read)",
}


def _left(r: dict, reason: str) -> tuple:
    return (r["member"], r["kind"], r["records"], r["bases"], reason)


# === GBRMGD ========================================================================================================

GBRMGD = "GBRMGD"
GBRMGD_FASTA = "gbrmgd_bins.fna.gz"
GBRMGD_DIR = "prokaryote_bins/"
GBRMGD_SOURCE = {
    "source_name": "tower_gbr_mgd",
    "source_release": "Great Barrier Reef Microbial Genomes Database (GBR_MGD), zenodo 20503966: "
                      "prokaryote_bins_5283.tar.bz2; staged by c83_coral_inputs.py",
    "license": "none stated: the job list copied the Zenodo record's licence field as the literal {} (DEC-80 (2))",
    "redistributable": False,  # DEC-80 (3): no recorded licence permits it
}
GBRMGD_RETRIEVED = "2026-11-26T07:15:49Z"  # tower's SOURCE_URLS.tsv
GBRMGD_STATUS_HEAD = ["genome", "site", "binner", "r232_status", "gtdb_classification", "contigs", "bp"]


def plan_gbrmgd(census_rows: list[dict]) -> tuple[dict[str, dict], list[tuple]]:
    """Every DNA member under prokaryote_bins/ is a prokaryote bin, genome id its stem (`<site>__<binner and bin>`);
    the archive carries no quality table, so `mag_mq` (the set's floor is not stated in what was fetched) and no
    taxonomic hint; no accession, so `not_visible` to the r232 test. Anything else is left out with its reason."""
    plan: dict[str, dict] = {}
    left: list[tuple] = []
    for r in census_rows:
        m = r["member"]
        if r["kind"] == "dna" and m.startswith(GBRMGD_DIR) and "/" not in m[len(GBRMGD_DIR):]:
            g = r["stem"]
            site, _, binner = g.partition("__")
            plan[m] = {"genome": g, "kind": "prokaryote", "quality": "mag_mq", "hint": "", "extra": site,
                       "status": (g, site, binner, "not_visible", "")}
        elif r["kind"] == "dna":
            left.append(_left(r, f"DNA outside {GBRMGD_DIR} (not one of the set's prokaryote bins)"))
        else:
            left.append(_left(r, LEAVE_REASON.get(r["kind"], r["kind"])))
    return plan, left


# === RMD1 ==========================================================================================================

RMD1 = "RMD1"
RMD1_FASTA = "rmd1_reef_genomes.fna.gz"
RMD1_SOURCE = {
    "source_name": "tower_batch6b_coral_rmd1",
    "source_release": "Reef Microbiomics Database v1.0 (Tara Pacific), RMDv1-genomes.tar.gz from sunagawalab.ethz.ch, "
                      "described by genomes-aggregated-summary.tsv of zenodo 17844029 (tower_rmd1); staged by "
                      "c83_coral_inputs.py",
    "license": "none stated for the tarball: the job list's 'open-academic' is the worker's label, not a licence; "
               "the tables' Zenodo record states cc-by-4.0 (DEC-80 (2))",
    "redistributable": False,  # DEC-80 (3): no recorded licence covers the genome tarball
}
RMD1_RETRIEVED = "2026-11-27T06:14:12Z"  # tower's SOURCE_URLS.tsv, RMDv1-genomes.tar.gz
RMD1_NEED = ("data_type", "biosample", "genome", "GTDBTK_TAXONOMY", "COMPLETENESS", "CONTAMINATION", "#SCAFFOLDS",
             "GENOME_SIZE", "secondary_cluster_representative")
RMD1_KIND = {"TPAC-CORAL-MAGS": "mag", "EXTERNAL-MAGS": "mag", "EXTERNAL-ISOG": "isolate"}
RMD1_STATUS_HEAD = ["genome", "data_type", "genome_kind", "biosample", "gtdb_r232_accession", "via", "r232_status",
                    "in_gtdb_r232", "is_gtdb_representative", "gtdb_representative", "gtdb_species",
                    "in_globdb_r232", "globdb_r232_id", "gtdb_r232_genomes_of_biosample", "gtdb_classification",
                    "completeness", "contamination", "secondary_cluster_representative", "contigs", "bp"]


def rmd1_table(path: Path) -> dict[str, dict]:
    """genomes-aggregated-summary.tsv (a quoted TSV): genome id -> its row; refuses a missing column, a repeated
    genome and a data type it cannot name as MAG or isolate."""
    with open(path, newline="", encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t", quotechar='"')
        miss = [c for c in RMD1_NEED if c not in (rd.fieldnames or [])]
        if miss:
            raise Refused(f"{path}: no column {miss[0]!r}")
        out: dict[str, dict] = {}
        for r in rd:
            g = r["genome"].strip()
            if g in out:
                raise Refused(f"{path}: genome {g} has two rows")
            if r["data_type"] not in RMD1_KIND:
                raise Refused(f"{path}: genome {g} has the data type {r['data_type']!r}, neither MAG nor isolate")
            out[g] = {k: r[k].strip() for k in RMD1_NEED}
    return out


def plan_rmd1(census_rows: list[dict], table: dict[str, dict], index_dir: Path, gtdb_meta_dir: Path,
              ) -> tuple[dict[str, dict], list[tuple], dict]:
    """Every DNA member whose stem is a genome of the table is staged (`kind="prokaryote"`: the table is GTDB-Tk's
    and CheckM's, prokaryotes only); a MAG's quality from its completeness and contamination, an isolate's
    `scaffold` (a draft: the table gives no assembly level). The r232 test: an isolate's BioSample through GTDB
    R232's metadata (`ncbi_biosample`) to its accession, then card C-74's `lookup`; a MAG's BioSample is its
    metagenome's, shared by every MAG of that sample, so a MAG is `not_visible` and the number of GTDB R232 genomes of
    that BioSample is recorded beside it. Members that are not DNA, or not in the table, are left out."""
    plan: dict[str, dict] = {}
    left: list[tuple] = []
    for r in census_rows:
        if r["kind"] != "dna":
            left.append(_left(r, LEAVE_REASON.get(r["kind"], r["kind"])))
        elif r["stem"] not in table:
            left.append(_left(r, "DNA whose name is no genome of genomes-aggregated-summary.tsv"))
        else:
            plan[r["member"]] = {"genome": r["stem"]}
    rows = {v["genome"]: table[v["genome"]] for v in plan.values()}
    bios = {d["biosample"] for d in rows.values() if d["biosample"]}
    by_bs, _ = c81.gtdb_crosswalk(gtdb_meta_dir, bios, set())
    via = {g: by_bs[d["biosample"]][0] for g, d in rows.items()
           if RMD1_KIND[d["data_type"]] == "isolate" and by_bs.get(d["biosample"])}
    look = {x["accession"]: x for x in r232_index.lookup(index_dir, sorted(set(via.values())))}

    def yn(v: bool) -> str:
        return "yes" if v else "no"

    for m, v in plan.items():
        g = v["genome"]
        d = rows[g]
        gk = RMD1_KIND[d["data_type"]]
        x = look.get(via[g]) if g in via else None
        status = c81.r232_status(x, gk == "isolate" and bool(d["biosample"]))
        v.update({"kind": "prokaryote", "hint": d["GTDBTK_TAXONOMY"], "extra": d["data_type"],
                  "quality": quality_of_mag(d["COMPLETENESS"], d["CONTAMINATION"]) if gk == "mag" else "scaffold",
                  "genome_kind": gk,
                  "status": (g, d["data_type"], gk, d["biosample"], via.get(g, ""), "biosample" if g in via else "",
                             status, yn(bool(x and x["in_gtdb_r232"])),
                             str(x["is_representative"]) if x and x["in_gtdb_r232"] else "",
                             (x["representative"] or "") if x else "", (x["species"] or "") if x else "",
                             yn(bool(x and x["in_globdb_r232"])), (x["globdb_r232_id"] or "") if x else "",
                             len(by_bs.get(d["biosample"], [])), d["GTDBTK_TAXONOMY"], d["COMPLETENESS"],
                             d["CONTAMINATION"], d["secondary_cluster_representative"])})
    info = {"table_rows": len(table), "table_rows_by_type": dict(sorted(Counter(d["data_type"] for d in
                                                                                table.values()).items())),
            "table_rows_not_in_archive": len(set(table) - {v["genome"] for v in plan.values()}),
            "isolates_found_by_biosample": len(via),
            "mags_whose_biosample_gtdb_r232_holds": sum(1 for g, d in rows.items() if RMD1_KIND[d["data_type"]] ==
                                                         "mag" and by_bs.get(d["biosample"]))}
    return plan, left, info


# === ANI against held sets (recorded, not cut) and exact overlap with a stage staged later ========================

def genomes_where(meta: Path, column: str, values: set[str]) -> set[str]:
    """The genomes of a staged set's records.tsv whose COLUMN is one of VALUES (to take a set's MAGs and not its
    unbinned contigs or viral contigs)."""
    out: set[str] = set()
    with open(c81._records_tsv(Path(meta)), newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if column not in (rd.fieldnames or []):
            raise Refused(f"{meta}: records.tsv has no column {column}")
        gcol = "genome" if "genome" in (rd.fieldnames or []) else "genome_id"
        for r in rd:
            if r[column] in values:
                out.add(r[gcol])
    return out


def all_archaea_refs(fna: Path, list_tsv: Path, out_dir: Path) -> dict[str, Path]:
    """Every genome of card C-75's r226 archaea (one FASTA in list order, `first_record` and `records` per genome),
    one file each; the sets staged here carry no genus for every genome, so no genus is used to choose."""
    want: dict[int, str] = {}
    with open(list_tsv, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for i in range(int(r["first_record"]), int(r["first_record"]) + int(r["records"])):
                want[i] = r["genome"]
    out_dir.mkdir(parents=True, exist_ok=True)
    files: dict[str, Path] = {}
    cur, fh = None, None
    try:
        for i, r in enumerate(fasta.iter_records(Path(fna))):
            g = want.get(i)
            if g is None:
                continue
            if g != cur:  # a genome's records are consecutive in the FASTA (the list's first_record and records)
                if fh is not None:
                    fh.close()
                if g in files:
                    raise Refused(f"{list_tsv}: the records of {g} are not consecutive")
                files[g] = out_dir / (re.sub(r"[^A-Za-z0-9_.-]", "_", g) + ".fna")
                cur, fh = g, open(files[g], "wb")
            fh.write(b">" + r.accession.encode() + b"\n" + r.seq + b"\n")
    finally:
        if fh is not None:
            fh.close()
    return files


def ani(meta: Path, fasta_paths: list[Path], out: Path, *, refs: dict[str, tuple[Path, list[Path], tuple | None]],
        archaea: tuple[Path, Path] | None, threads: int = 8, label: str = "Q") -> dict:
    """The staged set's genomes against each reference group in its own `skani dist -t THREADS --min-af 50` process
    (skani's default sketch, so that no process holds every sketch): REFS name -> (meta dir, FASTA files, optional
    (column, values) that picks the group's genomes), and card C-75's r226 archaea whole. A genome meets a reference
    at ANI >= 96 over an aligned fraction >= 50 (the larger of the two, card C-66's rule); the best hit over groups
    is kept. Recorded, never cut."""
    out = Path(out)
    work = out / "work"
    shutil.rmtree(work, ignore_errors=True)
    q_files = c81.split_genomes(Path(meta), list(fasta_paths), work / label)
    (work / "q.txt").write_text("".join(f"{p}\n" for p in sorted(map(str, q_files.values()))))
    query_of = {str(p): g for g, p in q_files.items()}
    groups: list[tuple[str, object]] = [(n, v) for n, v in refs.items()]
    if archaea is not None:
        groups.append(("r226_archaea", None))
    res: dict = {"queries": len(q_files), "threads": threads, "groups": {},
                 "threshold": {"ani_min": c81.ANI_SAG, "af_min_larger_of_two": c81.AF_SAG}}
    hits_all = []
    for name, v in groups:
        d = work / f"ref_{name}"
        unlisted: list[str] = []
        if v is None:
            files = all_archaea_refs(archaea[0], archaea[1], d)
        else:
            m, fas, pick = v
            keep = genomes_where(m, pick[0], set(pick[1])) if pick else None
            files = c81.split_genomes(m, fas, d, keep=keep, unlisted=unlisted)
        ref_of = {str(p): (name, g) for g, p in files.items()}
        entry = {"references": len(files), "records_in_no_genome": len(unlisted),
                 "picked_by": list(v[2]) if v is not None and v[2] else None}
        if files:
            (work / f"r_{name}.txt").write_text("".join(f"{p}\n" for p in sorted(ref_of)))
            table = out / f"vs_{name}.tsv"
            t0 = c39_ani.child_cpu()
            subprocess.run(["skani", "dist", "-t", str(threads), "--min-af", str(c81.AF_SAG), "--ql",
                            str(work / "q.txt"), "--rl", str(work / f"r_{name}.txt"), "-o", str(table)], check=True)
            entry["cpu_s"] = round(c39_ani.child_cpu() - t0, 1)
            h = c81.best_hits(table, query_of, ref_of)
            entry["meet"] = len(h)
            hits_all.append(h)
        shutil.rmtree(d, ignore_errors=True)
        res["groups"][name] = entry
    best = c81.merge_best(hits_all)
    c81._write_tsv(out / "ani_hits.tsv", ["genome", "stage", "ref_genome", "ani", "af"],
                   [(g, h["stage"], h["genome"], h["ani"], h["af"]) for g, h in sorted(best.items())])
    res["meet"] = len(best)
    res["meet_by_stage"] = dict(sorted(Counter(h["stage"] for h in best.values()).items()))
    (out / "summary.json").write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    shutil.rmtree(work, ignore_errors=True)
    return res


def overlap(meta: Path, others: dict[str, Path], work: Path, threads: int = 2, memory_limit: str = "4GB") -> dict:
    """Contigs (and bases) of a staged set's records.tsv that each of OTHERS (name -> dir with genomes.json) also
    lists, by exact contig_uid, each counted against every other set (not by ownership): for a stage staged after
    this one and earlier in the chain, which card A-29 orders."""
    rec = c81._records_tsv(Path(meta))
    length: dict[str, int] = {}
    with open(rec, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            length.setdefault(r["contig_uid"], int(r["length"]))
    res = {}
    for name, d in others.items():
        held = c58._owners_from_records(rec, None, [], {name: Path(d)}, work=Path(work), threads=threads,
                                        memory_limit=memory_limit)
        res[name] = {"contigs": len(held), "bases": sum(length[cu.hex()] for cu in held)}
    return res


# === command line ==================================================================================================

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Stage the coral-reef sets GBRMGD and RMD1 (card C-83, DEC-85)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("census")
    c.add_argument("tar", type=Path)
    c.add_argument("out_tsv", type=Path)
    c.add_argument("--report-out", type=Path, default=None)
    g = sub.add_parser("gbrmgd-stage")
    for n in ("tar", "census_tsv", "meta_out", "fasta_out"):
        g.add_argument(n, type=Path)
    c81._stage_args(g, GBRMGD_RETRIEVED)
    g.add_argument("--workers", type=int, default=1)
    n = sub.add_parser("ani")
    n.add_argument("meta", type=Path)
    n.add_argument("fasta", help="FASTA[,FASTA]")
    n.add_argument("out", type=Path)
    n.add_argument("--ref", action="append", default=[], metavar="NAME=META_DIR:FASTA[,FASTA][@COLUMN=V1|V2]")
    n.add_argument("--archaea", default="", metavar="FASTA:LIST_TSV")
    n.add_argument("--threads", type=int, default=8)
    n.add_argument("--label", default="Q")
    n.add_argument("--report-out", type=Path, default=None)
    for cmd in ("rmd1-plan", "rmd1-stage"):
        r = sub.add_parser(cmd)
        for n_ in ("tar", "census_tsv", "summary_tsv", "index_dir", "gtdb_meta_dir"):
            r.add_argument(n_, type=Path)
        if cmd == "rmd1-stage":
            r.add_argument("meta_out", type=Path)
            r.add_argument("fasta_out", type=Path)
            c81._stage_args(r, RMD1_RETRIEVED)
            r.add_argument("--workers", type=int, default=1)
        else:
            r.add_argument("--report-out", type=Path, default=None)
    o = sub.add_parser("overlap")
    o.add_argument("meta", type=Path)
    o.add_argument("--with", dest="others", action="append", default=[], metavar="NAME=DIR")
    o.add_argument("--work", type=Path, required=True)
    o.add_argument("--report-out", type=Path, default=None)
    a = ap.parse_args(argv)

    def emit(rep: dict) -> None:
        txt = json.dumps(rep, indent=1, sort_keys=True)
        if a.report_out is not None:
            a.report_out.parent.mkdir(parents=True, exist_ok=True)
            a.report_out.write_text(txt + "\n")
        print(txt)

    try:
        if a.cmd == "census":
            emit(census(a.tar, a.out_tsv))
            return 0
        if a.cmd == "ani":
            refs = {}
            for k, v in c81._pairs(a.ref, "--ref").items():
                v, _, pick = v.partition("@")
                m, _, fas = v.partition(":")
                if not m or not fas:
                    raise Refused(f"--ref {k}={v}: expected NAME=META_DIR:FASTA[,FASTA][@COLUMN=V1|V2]")
                p = None
                if pick:
                    col, _, vals = pick.partition("=")
                    if not col or not vals:
                        raise Refused(f"--ref {k}: expected @COLUMN=V1|V2, not @{pick}")
                    p = (col, tuple(vals.split("|")))
                refs[k] = (Path(m), [Path(x) for x in fas.split(",") if x], p)
            arc = None
            if a.archaea:
                fa, _, lst = a.archaea.partition(":")
                if not fa or not lst:
                    raise Refused(f"--archaea {a.archaea}: expected FASTA:LIST_TSV")
                arc = (Path(fa), Path(lst))
            emit(ani(a.meta, [Path(x) for x in a.fasta.split(",") if x], a.out, refs=refs, archaea=arc,
                     threads=a.threads, label=a.label))
            return 0
        if a.cmd == "overlap":
            others = {k: Path(v) for k, v in c81._pairs(a.others, "--with").items()}
            emit(overlap(a.meta, others, a.work))
            return 0
        if a.cmd == "rmd1-plan":
            rows = read_census(a.census_tsv)
            plan, left, info = plan_rmd1(rows, rmd1_table(a.summary_tsv), a.index_dir, a.gtdb_meta_dir)
            by_member = {r["member"]: r for r in rows}
            kinds = Counter(v["genome_kind"] for v in plan.values())
            emit({**info, "planned_members": len(plan), "planned_by_genome_kind": dict(sorted(kinds.items())),
                  "planned_bases": sum(int(by_member[m]["bases"]) for m in plan),
                  "planned_bases_by_genome_kind": {k: sum(int(by_member[m]["bases"]) for m, v in plan.items()
                                                          if v["genome_kind"] == k) for k in sorted(kinds)},
                  "planned_by_quality": dict(sorted(Counter(v["quality"] for v in plan.values()).items())),
                  "planned_r232_status": dict(sorted(Counter(v["status"][6] for v in plan.values()).items())),
                  "left_out": len(left), "left_out_by_reason": dict(sorted(Counter(x[-1] for x in left).items())),
                  "left_out_bases": sum(int(x[3] or 0) for x in left)})
            return 0
        owners = {st: a.chain_home / st for st, _ in a26.R26_11_4_CHAIN} if a.chain_home is not None else {}
        later = {k: Path(v) for k, v in c81._pairs(a.later, "--later").items()}
        roots = [x for x in a.root_stages.split(",") if x]
        common = dict(root=a.root, root_stages=roots, owner_inputs=owners, later_inputs=later,
                      inputs_home=a.inputs_home, retrieved_at=a.retrieved_at, work=a.work, threads=a.threads,
                      memory_limit=a.memory_limit, workers=a.workers)
        if a.cmd == "gbrmgd-stage":
            name = GBRMGD
            plan, left = plan_gbrmgd(read_census(a.census_tsv))
            rep = stage_set(GBRMGD, a.tar, plan, a.meta_out, a.fasta_out, fasta_name=GBRMGD_FASTA,
                            source=GBRMGD_SOURCE, status_head=GBRMGD_STATUS_HEAD, record_extra="site",
                            left_out=left, **common)
        elif a.cmd == "rmd1-stage":
            name = RMD1
            plan, left, info = plan_rmd1(read_census(a.census_tsv), rmd1_table(a.summary_tsv), a.index_dir,
                                         a.gtdb_meta_dir)
            rep = stage_set(RMD1, a.tar, plan, a.meta_out, a.fasta_out, fasta_name=RMD1_FASTA, source=RMD1_SOURCE,
                            status_head=RMD1_STATUS_HEAD, record_extra="data_type", left_out=left, **common)
            rep.update(info)
            kinds = {v["genome"]: v["genome_kind"] for v in plan.values()}
            rep["by_genome_kind"] = dict(sorted(Counter(kinds.values()).items()))
        else:
            return 2
        rep["owner_stages"] = [*(roots if a.root else []), *owners]
        rep["later_stages"] = [*later]
        rep["s2_check"] = c81._check(name, a.meta_out / name)
        emit(rep)
        return 0 if rep["s2_check"]["exit"] == 0 and rep["s2_check"]["problems"] == [] else 1
    except Refused as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
