"""Exact-byte FASTA reading and writing for the sequence store.

Records keep their header line and sequence bytes exactly as in the source
(case and IUPAC codes preserved). `line_width` is the wrap width of the record,
or 0 when the wrapping is irregular, in which case round trips compare the
unwrapped sequence.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Record:
    header: str  # without the leading '>'
    seq: bytes
    line_width: int

    @property
    def accession(self) -> str:
        return self.header.split(None, 1)[0] if self.header else ""


def _open(path: Path):
    return gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")


def normalized_lines(path: Path) -> list[bytes]:
    """Decompressed lines with CRLF folded to LF and blank lines dropped."""
    with _open(path) as fh:
        return [ln.rstrip(b"\r\n") for ln in fh if ln.strip()]


def _width(lines: list[bytes]) -> int:
    if not lines:
        return 0
    w = len(lines[0])
    if all(len(ln) == w for ln in lines[:-1]) and len(lines[-1]) <= w:
        return w
    return 0


def iter_records(path: Path) -> Iterator[Record]:
    """The records of `path` one at a time, reading the file as it goes (card B-15).

    Only the record being assembled is held. Lines are taken as `normalized_lines` gives them (CRLF
    folded to LF, blank lines dropped); `read` is this, listed.
    """
    header: str | None = None
    seq_lines: list[bytes] = []
    with _open(path) as fh:
        for raw in fh:
            if not raw.strip():
                continue
            ln = raw.rstrip(b"\r\n")
            if ln.startswith(b">"):
                if header is not None:
                    rec = Record(header, b"".join(seq_lines), _width(seq_lines))
                    seq_lines = []
                    yield rec
                header = ln[1:].decode()
                seq_lines = []
            else:
                if header is None:
                    raise ValueError(f"{path}: sequence before the first header")
                seq_lines.append(ln)
    if header is not None:
        yield Record(header, b"".join(seq_lines), _width(seq_lines))


def read(path: Path) -> list[Record]:
    return list(iter_records(path))


def format_record(header: str, seq: bytes, line_width: int) -> list[bytes]:
    lines = [b">" + header.encode()]
    if line_width <= 0:
        lines.append(seq)
    else:
        lines.extend(seq[i : i + line_width] for i in range(0, len(seq), line_width))
    return lines
