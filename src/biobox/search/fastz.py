"""zstd-compressed FASTA for search-pack shards (`part-<k>.faa.zst`, `part-<k>.fna.zst`).

Records are `>name` plus the sequence on one line. Shards are standard zstd frames
(pyarrow's codec), so `zstd -dc` reads them.

Readers that do not hold a shard (card W-19). A pack's shard is up to 2 GiB of text (`pack.SHARD_TEXT_BYTES`), and
a search reads a target of many shards, so nothing in the search path reads a whole shard into memory:
`blocks` hands the decompressed stream to Easel's FASTA parser (pyhmmer's `SequenceFile`, digital mode) and
yields digitized blocks of about a given number of residues; `lookup` finds the records of a few names in one
shard, a chunk (`CHUNK` bytes of text) at a time, stopping at the place of the last name sought, since a pack's
shards hold their records in name order (`biobox.search.pack`, module docstring); `first_name` reads a shard's
first header only. `read_records` reads a whole shard and is for small shards and tests.
"""

from __future__ import annotations

import io
from pathlib import Path
from typing import Iterable, Iterator

import pyarrow as pa

CHUNK = 1 << 24  # bytes of decompressed text read at a time by `lookup` and `blocks`


class ShardWriter:
    """Writes records to `<dir>/part-<k><suffix>`, starting a new shard once the current one
    holds `shard_text_bytes` of uncompressed text. Returns per-shard (path, n_seq, n_residues)."""

    def __init__(self, directory: Path, suffix: str, shard_text_bytes: int):
        self.dir, self.suffix, self.limit = Path(directory), suffix, shard_text_bytes
        self.dir.mkdir(parents=True, exist_ok=True)
        self.shards: list[dict] = []
        self._out = None
        self._held = 0

    def _open(self) -> None:
        path = self.dir / f"part-{len(self.shards)}{self.suffix}"
        self._out = pa.CompressedOutputStream(str(path), "zstd")
        self.shards.append({"path": path, "n_seq": 0, "n_residues": 0})
        self._held = 0

    def write(self, name: str, seq: str) -> None:
        if self._out is None or self._held >= self.limit:
            self.close()
            self._open()
        text = f">{name}\n{seq}\n".encode("ascii")
        self._out.write(text)
        self._held += len(text)
        self.shards[-1]["n_seq"] += 1
        self.shards[-1]["n_residues"] += len(seq)

    def close(self) -> list[dict]:
        if self._out is not None:
            self._out.close()
            self._out = None
        return self.shards


def write_records(directory: Path, suffix: str, records: Iterable[tuple[str, str]], shard_text_bytes: int) -> list[dict]:
    w = ShardWriter(directory, suffix, shard_text_bytes)
    for name, seq in records:
        w.write(name, seq)
    if not w.shards:  # an empty target is one empty shard, so it is still listed and hashed
        w._open()
    return w.close()


def read_records(path: Path) -> Iterator[tuple[str, str]]:
    """Every record of a shard; the whole shard is read into memory first (small shards and tests only)."""
    with pa.CompressedInputStream(pa.OSFile(str(path)), "zstd") as fh:
        data = fh.read().decode("ascii")
    name = None
    for line in data.splitlines():
        if line.startswith(">"):
            name = line[1:]
        elif name is not None:
            yield name, line
            name = None


def zst(data: bytes) -> bytes:
    """A zstd frame of `data` (search_profile blobs); read back with `unzst` or `zstd -d`."""
    sink = pa.BufferOutputStream()
    with pa.CompressedOutputStream(sink, "zstd") as out:
        out.write(data)
    return sink.getvalue().to_pybytes()


def unzst(blob: bytes) -> bytes:
    with pa.CompressedInputStream(pa.BufferReader(blob), "zstd") as fh:
        return fh.read()


def _open(path: Path):
    return pa.CompressedInputStream(pa.OSFile(str(path)), "zstd")


def first_name(path: Path) -> str | None:
    """The name of a shard's first record, None for an empty shard; only the start of the shard is decompressed."""
    head = b""
    with _open(path) as fh:
        while b"\n" not in head:
            chunk = fh.read(1 << 16)
            if not chunk:
                break
            head += chunk
    if not head.startswith(b">"):
        return None
    end = head.find(b"\n")
    return head[1:end if end >= 0 else len(head)].decode("ascii")


def lookup(path: Path, names: Iterable[str]) -> dict[str, str]:
    """{name: sequence} of those of `names` that are records of the shard `path`, whose records are in name order
    (every pack shard). The shard is read `CHUNK` bytes at a time, and the read stops once the next record's name
    sorts after every name still sought; at most a chunk and one partial record are held."""
    want = sorted(set(names))
    found: dict[str, str] = {}
    if not want:
        return found
    keys = {n: b">" + n.encode("ascii") + b"\n" for n in want}
    buf = b""
    with _open(path) as fh:
        while want:
            chunk = fh.read(CHUNK)
            eof = not chunk
            buf += chunk
            # the complete records: up to the start of the last header, which may still be cut (all of it at eof)
            cut = len(buf) if eof else buf.rfind(b"\n>") + 1
            if cut <= 0:
                continue  # one record longer than the chunk: read on
            done, buf = buf[:cut], buf[cut:]  # `done` and `buf` each start at a record
            for n in list(want):
                key = keys[n]
                at = done.find(b"\n" + key)
                start = at + 1 + len(key) if at >= 0 else (len(key) if done.startswith(key) else -1)
                if start < 0:
                    continue
                end = done.find(b"\n", start)
                found[n] = done[start:end if end >= 0 else len(done)].decode("ascii")
                want.remove(n)
            if eof or not want:
                break
            nl = buf.find(b"\n")
            if nl > 0 and buf[1:nl].decode("ascii") > want[-1]:
                break  # every name still sought sorts before the next record: none of them is in the shard
    return found


class _Stream(io.RawIOBase):
    """A pyarrow input stream (with its first bytes already read) as a raw binary file, for pyhmmer."""

    def __init__(self, stream, head: bytes):
        self._s, self._head = stream, head

    def readable(self) -> bool:
        return True

    def readinto(self, b) -> int:
        if self._head:
            n = min(len(b), len(self._head))
            b[:n], self._head = self._head[:n], self._head[n:]
            return n
        data = self._s.read(len(b))
        b[:len(data)] = data
        return len(data)


def blocks(path: Path, alphabet, residues: int) -> Iterator:
    """The records of a shard as digitized `pyhmmer.easel.DigitalSequenceBlock`s of at least `residues` residues
    each (the last one may hold fewer), read by Easel's FASTA parser from the decompressed stream: one block is held
    at a time (card W-19). Nothing for an empty shard."""
    import pyhmmer

    with _open(path) as raw:
        head = raw.read(1 << 16)
        if not head:
            return
        with pyhmmer.easel.SequenceFile(_Stream(raw, head), "fasta", digital=True, alphabet=alphabet) as sf:
            while True:
                block = sf.read_block(residues=max(1, int(residues)))
                if not len(block):
                    return
                yield block
