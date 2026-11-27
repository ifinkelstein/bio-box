"""Card R-19: scan a staged FASTA (plain or gzip) for bytes that are not sequence.

Standard library only, one file per process, bounded memory (a line is read in pieces of at most 1 MiB; nothing but
counters and the current header is held). Records are split as `biobox.fasta.iter_records` splits them: a line that
starts with '>' opens a record, trailing CR/LF is stripped, whitespace-only lines are skipped.

Reported per file: records, bases (sequence bytes), headers that are not UTF-8, headers that are not printable ASCII
(bytes 0x20 to 0x7e and tab), records holding a byte outside the IUPAC nucleotide letters `ACGTURYSWKMBDHVN` (either
case), the count of such bytes, empty records, lines before the first header, a read error (truncated gzip, trailing
bytes that are not gzip), the first record of each kind, and the first offending record with its header, the header
before it, its reasons and its first bad byte. The compression is taken from the file's first two bytes, not its name.

Subcommands:
    scan PATH OUT_JSON [--stage S]          one file; OUT_JSON written under a temporary name and renamed;
                                            exit 0 clean, 1 damaged, 2 unreadable
    scan-stage OBJECTS_JSON OUT_DIR --stage S
                                            every path of the stage's objects.json, each `scan` in its own process,
                                            one after another; a file whose JSON exists is not scanned again; then
                                            OUT_DIR/<S>.stage.json (exit 0 once written)
    table OUT_DIR                           one Markdown row per scanned file of every <S>.stage.json in OUT_DIR
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import io
import json
import os
import re
import resource
import socket
import subprocess
import sys
import time
import zlib
import gzip
from pathlib import Path

PIECE = 1 << 20
HEADER_CAP = 1 << 20
IUPAC = b"ACGTURYSWKMBDHVNacgturyswkmbdhvn"
PRINTABLE = bytes(range(0x20, 0x7F)) + b"\t"
_NON_IUPAC_RE = re.compile(b"[^" + re.escape(IUPAC) + b"]")
_SHOW = 200

REASONS = (
    "sequence_before_first_header",
    "header_not_utf8",
    "header_not_printable_ascii",
    "header_too_long",
    "non_iupac_byte",
    "empty_record",
)


class _HashingRaw(io.RawIOBase):
    """Reads a file and hashes every byte read, so the sha256 costs no second pass."""

    def __init__(self, fh):
        self._fh = fh
        self.sha = hashlib.sha256()
        self.nbytes = 0

    def readable(self) -> bool:
        return True

    def readinto(self, b) -> int:
        n = self._fh.readinto(b)
        if n:
            self.sha.update(memoryview(b)[:n])
            self.nbytes += n
        return n or 0


def _show(header: bytes) -> str:
    return header.decode("utf-8", "backslashreplace")[:_SHOW]


class _Scan:
    def __init__(self) -> None:
        self.records = 0
        self.bases = 0
        self.counts = {
            "headers_not_utf8": 0,
            "headers_not_printable_ascii": 0,
            "headers_too_long": 0,
            "records_with_non_iupac": 0,
            "non_iupac_bytes": 0,
            "empty_records": 0,
            "lines_before_first_header": 0,
        }
        self.offending_records = 0
        self.bad_values: dict[int, int] = {}  # byte value -> count, over every sequence byte outside IUPAC
        self.first_record_of: dict[str, int] = {}
        self.first_offending: dict | None = None
        self.prev_header: str | None = None
        # the open record
        self.header = bytearray()
        self.header_long = False
        self.seq_len = 0
        self.bad_bytes = 0
        self.first_bad: tuple[int, int] | None = None  # (byte, offset)

    def _mark(self, reason: str, record: int) -> None:
        self.first_record_of.setdefault(reason, record)

    def before_header(self) -> None:
        self.counts["lines_before_first_header"] += 1
        self._mark("sequence_before_first_header", 0)
        if self.first_offending is None:
            self.first_offending = {
                "record": 0,
                "header": None,
                "previous_header": None,
                "reasons": ["sequence_before_first_header"],
            }

    def open_record(self, first_piece: bytes) -> None:
        self.close_record()
        self.records += 1
        self.header = bytearray()
        self.header_long = False
        self.seq_len = 0
        self.bad_bytes = 0
        self.first_bad = None
        self.header_piece(first_piece)

    def header_piece(self, piece: bytes) -> None:
        room = HEADER_CAP - len(self.header)
        if len(piece) > room:
            self.header_long = True
            piece = piece[: max(room, 0)]
        self.header += piece

    def seq_piece(self, piece: bytes) -> None:
        if not piece:
            return
        bad = piece.translate(None, IUPAC)
        if bad:
            self.bad_bytes += len(bad)
            for b in set(bad):
                self.bad_values[b] = self.bad_values.get(b, 0) + bad.count(b)
            if self.first_bad is None:
                m = _NON_IUPAC_RE.search(piece)
                off = m.start() if m else 0
                self.first_bad = (piece[off], self.seq_len + off)
        self.seq_len += len(piece)

    def close_record(self) -> None:
        if self.records == 0 or self.header is None:
            return
        n = self.records
        hdr = bytes(self.header)
        reasons: list[str] = []
        try:
            hdr.decode("utf-8")
        except UnicodeDecodeError:
            reasons.append("header_not_utf8")
            self.counts["headers_not_utf8"] += 1
        if hdr.translate(None, PRINTABLE):
            reasons.append("header_not_printable_ascii")
            self.counts["headers_not_printable_ascii"] += 1
        if self.header_long:
            reasons.append("header_too_long")
            self.counts["headers_too_long"] += 1
        if self.bad_bytes:
            reasons.append("non_iupac_byte")
            self.counts["records_with_non_iupac"] += 1
            self.counts["non_iupac_bytes"] += self.bad_bytes
        if self.seq_len == 0:
            reasons.append("empty_record")
            self.counts["empty_records"] += 1
        self.bases += self.seq_len
        shown = _show(hdr)
        if reasons:
            self.offending_records += 1
            for r in reasons:
                self._mark(r, n)
            if self.first_offending is None:
                first = {"record": n, "header": shown, "previous_header": self.prev_header, "reasons": reasons}
                if self.first_bad is not None:
                    first["first_bad_byte"] = f"0x{self.first_bad[0]:02x}"
                    first["first_bad_offset"] = self.first_bad[1]
                first["record_bases"] = self.seq_len
                self.first_offending = first
        self.prev_header = shown
        self.header = None  # closed

    def feed(self, fh) -> None:
        at_start = True
        in_header = False
        carry = b""  # trailing CRs of a piece that did not end its line
        self.header = None
        while True:
            piece = fh.readline(PIECE)
            if not piece:
                break
            ends = piece.endswith(b"\n")
            content = piece[:-1] if ends else piece
            if carry:
                content = carry + content
                carry = b""
            if ends:
                content = content.rstrip(b"\r")
            elif content.endswith(b"\r"):
                stripped = content.rstrip(b"\r")
                carry = content[len(stripped):]
                content = stripped
            if at_start:
                if ends and not content.strip():
                    continue  # a blank line, skipped as biobox.fasta skips it
                if content.startswith(b">"):
                    self.open_record(content[1:])
                    in_header = True
                elif self.records == 0:
                    self.before_header()
                    in_header = False
                else:
                    in_header = False
                    self.seq_piece(content)
            elif self.records == 0:
                pass  # the rest of a line before the first header, counted once
            elif in_header:
                self.header_piece(content)
            else:
                self.seq_piece(content)
            at_start = ends
        self.close_record()


def _compression(buffered: io.BufferedReader) -> str:
    head = buffered.peek(2)[:2]
    return "gzip" if head == b"\x1f\x8b" else "none"


def scan_file(path: Path | str) -> dict:
    """Scan one FASTA, plain or gzip, and return its report (see the module docstring)."""
    path = Path(path)
    t0 = time.monotonic()
    r0 = resource.getrusage(resource.RUSAGE_SELF)
    st = path.stat()
    sc = _Scan()
    read_error = None
    with open(path, "rb", buffering=0) as raw_fh:
        raw = _HashingRaw(raw_fh)
        buffered = io.BufferedReader(raw, buffer_size=PIECE)
        comp = _compression(buffered)
        try:
            if comp == "gzip":
                with gzip.GzipFile(fileobj=buffered, mode="rb") as gz:
                    sc.feed(gz)
            else:
                sc.feed(buffered)
        except (OSError, EOFError, zlib.error) as e:  # gzip.BadGzipFile is an OSError
            read_error = f"{type(e).__name__}: {e}"
            if sc.header is not None:
                sc.close_record()
        while buffered.read(PIECE):  # the rest of the bytes, for the sha256
            pass
    r1 = resource.getrusage(resource.RUSAGE_SELF)
    c = sc.counts
    damaged = bool(
        read_error
        or sc.records == 0
        or c["headers_not_utf8"]
        or c["headers_not_printable_ascii"]
        or c["headers_too_long"]
        or c["records_with_non_iupac"]
        or c["empty_records"]
        or c["lines_before_first_header"]
    )
    rep = {
        "path": str(path),
        "file_bytes": st.st_size,
        "sha256": raw.sha.hexdigest(),
        "bytes_hashed": raw.nbytes,
        "compression": comp,
        "status": "damaged" if damaged else "clean",
        "records": sc.records,
        "bases": sc.bases,
        **c,
        "offending_records": sc.offending_records,
        # the most frequent byte values outside IUPAC (at most 10), so a gap '-' tells itself from compressed data
        "non_iupac_byte_values": {f"0x{b:02x}": n for b, n in sorted(sc.bad_values.items(), key=lambda x: (-x[1], x[0]))[:10]},
        "distinct_non_iupac_byte_values": len(sc.bad_values),
        "read_error": read_error,
        "first_record_of": sc.first_record_of,
        "first_offending": sc.first_offending,
        "wall_s": round(time.monotonic() - t0, 3),
        "cpu_s": round((r1.ru_utime - r0.ru_utime) + (r1.ru_stime - r0.ru_stime), 3),
    }
    return rep


def _write_json_atomic(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    tmp.write_text(json.dumps(obj, indent=1, sort_keys=False) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def _script_sha256() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def object_paths(objects_json: Path | str) -> list[Path]:
    """The paths an objects.json names, absolute (a relative path is taken from the file's directory)."""
    objects_json = Path(objects_json)
    entries = json.loads(objects_json.read_text(encoding="utf-8"))
    out = []
    for e in entries:
        p = Path(e["path"])
        out.append(p if p.is_absolute() else objects_json.parent / p)
    return out


def cmd_scan(args: argparse.Namespace) -> int:
    path = Path(args.path)
    out = Path(args.out_json)
    meta = {
        "card": "R-19",
        "stage": args.stage,
        "host": socket.gethostname(),
        "python": sys.version.split()[0],
        "script_sha256": _script_sha256(),
        "started_at": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    try:
        rep = scan_file(path)
    except OSError as e:
        rep = {"path": str(path), "status": "unreadable", "read_error": f"{type(e).__name__}: {e}"}
        _write_json_atomic(out, {**meta, **rep})
        return 2
    rep = {**meta, **rep, "finished_at": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
    _write_json_atomic(out, rep)
    print(json.dumps({k: rep[k] for k in ("path", "status", "records", "bases", "wall_s")}), flush=True)
    return 0 if rep["status"] == "clean" else 1


_ROW_KEYS = (
    "status", "records", "bases", "headers_not_utf8", "headers_not_printable_ascii", "records_with_non_iupac",
    "non_iupac_bytes", "empty_records", "lines_before_first_header", "read_error", "first_offending",
    "sha256", "file_bytes", "compression", "wall_s", "cpu_s",
)


def cmd_scan_stage(args: argparse.Namespace) -> int:
    stage = args.stage
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = object_paths(args.objects_json)
    names = [p.name for p in paths]
    if len(set(names)) != len(names):
        print(f"{args.objects_json}: two objects share a file name", file=sys.stderr)
        return 2
    per_file = []
    for p in paths:
        out = out_dir / f"{stage}__{p.name}.json"
        if not out.is_file():
            subprocess.run(
                [sys.executable, "-I", str(Path(__file__).resolve()), "scan", str(p), str(out), "--stage", stage],
                check=False,
            )
        if out.is_file():
            rep = json.loads(out.read_text(encoding="utf-8"))
        else:
            rep = {"status": "no_report"}
        per_file.append({"file": p.name, "path": str(p), **{k: rep.get(k) for k in _ROW_KEYS}})
    summary = {
        "card": "R-19",
        "stage": stage,
        "objects_json": str(Path(args.objects_json)),
        "files": len(per_file),
        "damaged": [f["file"] for f in per_file if f["status"] != "clean"],
        "records": sum(f["records"] or 0 for f in per_file),
        "bases": sum(f["bases"] or 0 for f in per_file),
        "per_file": per_file,
    }
    _write_json_atomic(out_dir / f"{stage}.stage.json", summary)
    print(json.dumps({"stage": stage, "files": summary["files"], "damaged": summary["damaged"]}), flush=True)
    return 0


def cmd_table(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    print("| Stage | File | Status | Records | Bases | Headers not UTF-8 | Headers not printable ASCII "
          "| Records with non-IUPAC | Non-IUPAC bytes | Empty records | Lines before first header | Read error |")
    print("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for sj in sorted(out_dir.glob("*.stage.json")):
        s = json.loads(sj.read_text(encoding="utf-8"))
        for f in s["per_file"]:
            print("| " + " | ".join([
                f"`{s['stage']}`", f"`{f['file']}`", str(f["status"]),
                *(f"{f.get(k) or 0:,}" for k in (
                    "records", "bases", "headers_not_utf8", "headers_not_printable_ascii",
                    "records_with_non_iupac", "non_iupac_bytes", "empty_records", "lines_before_first_header",
                )),
                str(f.get("read_error") or "none"),
            ]) + " |")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Card R-19: scan staged FASTA for bytes that are not sequence")
    sub = p.add_subparsers(dest="cmd", required=True)
    s1 = sub.add_parser("scan")
    s1.add_argument("path")
    s1.add_argument("out_json")
    s1.add_argument("--stage", default="")
    s1.set_defaults(func=cmd_scan)
    s2 = sub.add_parser("scan-stage")
    s2.add_argument("objects_json")
    s2.add_argument("out_dir")
    s2.add_argument("--stage", required=True)
    s2.set_defaults(func=cmd_scan_stage)
    s3 = sub.add_parser("table")
    s3.add_argument("out_dir")
    s3.set_defaults(func=cmd_table)
    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
