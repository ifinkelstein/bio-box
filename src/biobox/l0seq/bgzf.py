"""bgzip FASTA with `.fai` and `.gzi`: the L0-seq storage format (DEC-17, K-13).

`write_bgzf_fasta` writes the BGZF layout samtools/htslib read: gzip members of at most
65,280 uncompressed bytes with the `BC` extra field, a 28-byte EOF block, a `.fai`
(name, length, offset, line bases, line width) and a `.gzi` (block starts as
(compressed, uncompressed) uint64 pairs, first block omitted). Output is deterministic
(MTIME 0, fixed level).

`BgzfFasta` reads a slice back through the two indexes and counts the blocks it
decompresses and the bytes it reads, which is how LY-06 scores a sequence query.
"""

from __future__ import annotations

import bisect
import os
import struct
import threading
import zlib
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path

BLOCK_DATA = 0xFF00  # bgzip's uncompressed block payload
EOF_BLOCK = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def _block(data: bytes, level: int) -> bytes:
    c = zlib.compressobj(level, zlib.DEFLATED, -15)
    payload = c.compress(data) + c.flush()
    bsize = 18 + len(payload) + 8  # header + deflate + CRC32/ISIZE
    header = struct.pack("<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, bsize - 1)
    return header + payload + struct.pack("<II", zlib.crc32(data), len(data))


@dataclass(frozen=True)
class FaiEntry:
    name: str
    length: int
    offset: int
    line_bases: int
    line_width: int

    def byte_of(self, pos: int) -> int:
        """Uncompressed file offset of sequence position `pos` (0-based)."""
        return self.offset + (pos // self.line_bases) * self.line_width + pos % self.line_bases


def write_bgzf_fasta(records, out: Path, line_width: int = 80, level: int = 6) -> dict:
    """Write (name, seq bytes) records as bgzip FASTA plus `.fai` and `.gzi`. Returns sizes.

    The three files are written under temporary names and renamed into place, the `.fa.gz`
    last, so a reader never sees a shard without its indexes.
    """
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = {s: out.with_name(f".{out.name}{s}.tmp-{os.getpid()}") for s in ("", ".fai", ".gzi")}
    fai: list[FaiEntry] = []
    gzi: list[tuple[int, int]] = []
    buf = bytearray()
    upos = cpos = 0  # uncompressed bytes emitted so far, compressed bytes written so far

    with open(tmp[""], "wb") as fh:
        def flush(final: bool = False):
            nonlocal buf, upos, cpos
            while len(buf) >= BLOCK_DATA or (final and buf):
                piece = bytes(buf[:BLOCK_DATA])
                del buf[:BLOCK_DATA]
                if upos:
                    gzi.append((cpos, upos))
                blk = _block(piece, level)
                fh.write(blk)
                cpos += len(blk)
                upos += len(piece)

        written = 0  # uncompressed offset of the next byte appended to buf
        for name, seq in records:
            header = b">" + name.encode() + b"\n"
            buf += header
            written += len(header)
            fai.append(FaiEntry(name, len(seq), written, line_width, line_width + 1))
            for i in range(0, len(seq), line_width):
                line = seq[i : i + line_width] + b"\n"
                buf += line
                written += len(line)
                if len(buf) >= BLOCK_DATA:
                    flush()
        flush(final=True)
        fh.write(EOF_BLOCK)
        cpos += len(EOF_BLOCK)

    tmp[".fai"].write_text(
        "".join(f"{e.name}\t{e.length}\t{e.offset}\t{e.line_bases}\t{e.line_width}\n" for e in fai))
    with open(tmp[".gzi"], "wb") as fh:
        fh.write(struct.pack("<Q", len(gzi)))
        for c, u in gzi:
            fh.write(struct.pack("<QQ", c, u))
    for s in (".fai", ".gzi", ""):
        os.replace(tmp[s], str(out) + s)
    return {"bytes": cpos, "uncompressed_bytes": upos, "blocks": len(gzi) + 1, "records": len(fai)}


def read_fai(path: Path) -> list[FaiEntry]:
    out = []
    for line in Path(path).read_text().splitlines():
        n, length, off, lb, lw = line.split("\t")[:5]
        out.append(FaiEntry(n, int(length), int(off), int(lb), int(lw)))
    return out


def read_gzi(path: Path) -> list[tuple[int, int]]:
    raw = Path(path).read_bytes()
    (n,) = struct.unpack_from("<Q", raw)
    if len(raw) != 8 + 16 * n:
        raise ValueError(f"{path}: {len(raw)} bytes for {n} entries")
    return [struct.unpack_from("<QQ", raw, 8 + 16 * i) for i in range(n)]


def block_starts(path: Path) -> list[tuple[int, int]]:
    """(compressed, uncompressed) start of every data block, from the block headers alone.

    Reads 18 bytes and the 4-byte ISIZE per block; the EOF block is not listed. Raises on a
    block without the BGZF `BC` field or on a missing EOF block.
    """
    out = []
    size = Path(path).stat().st_size
    c = u = 0
    with open(path, "rb") as fh:
        while c < size:
            fh.seek(c)
            h = fh.read(18)
            if len(h) < 18 or h[:4] != b"\x1f\x8b\x08\x04" or h[12:14] != b"BC":
                raise ValueError(f"{path}: no BGZF block header at {c}")
            bsize = struct.unpack_from("<H", h, 16)[0] + 1
            fh.seek(c + bsize - 4)
            (isize,) = struct.unpack("<I", fh.read(4))
            if isize == 0 and c + bsize == size:
                return out
            out.append((c, u))
            c += bsize
            u += isize
    raise ValueError(f"{path}: no EOF block")


class BgzfFasta:
    """Slice a bgzip FASTA through its `.fai` and `.gzi`.

    Counters: `blocks_read` (blocks decompressed), `compressed_bytes_read`,
    `uncompressed_bytes_decoded`. With `cache_blocks` > 0 a small LRU of decoded blocks
    serves repeated slices; a cache hit counts nothing. Reads use pread, so one reader is
    safe to share between threads.
    """

    def __init__(self, path: Path, cache_blocks: int = 0):
        self._fd: int | None = None
        self.path = Path(path)
        self.fai: dict[str, FaiEntry] = {e.name: e for e in read_fai(str(path) + ".fai")}
        pairs = read_gzi(str(path) + ".gzi")
        self._cstart = [0] + [c for c, _ in pairs]
        self._ustart = [0] + [u for _, u in pairs]
        self._size = self.path.stat().st_size
        self._cache: OrderedDict[int, bytes] = OrderedDict()
        self._cache_max = cache_blocks
        self._lock = threading.Lock()
        self.blocks_read = 0
        self.compressed_bytes_read = 0
        self.uncompressed_bytes_decoded = 0

    def close(self) -> None:
        if self._fd is not None:
            os.close(self._fd)
            self._fd = None

    def __del__(self):
        self.close()

    def _blocks_for(self, ulo: int, uhi: int) -> range:
        first = bisect.bisect_right(self._ustart, ulo) - 1
        last = bisect.bisect_right(self._ustart, uhi - 1) - 1
        return range(first, last + 1)

    def units_for(self, name: str, start: int, end: int) -> int:
        """Blocks a slice touches, from the indexes alone (no I/O)."""
        if start == end:
            return 0
        e = self.fai[name]
        return len(self._blocks_for(e.byte_of(start), e.byte_of(end - 1) + 1))

    def _decoded(self, b: int) -> bytes:
        hit = self._cache.pop(b, None)
        if hit is not None:
            self._cache[b] = hit
            return hit
        if self._fd is None:
            self._fd = os.open(self.path, os.O_RDONLY)
        cstart = self._cstart[b]
        cend = self._cstart[b + 1] if b + 1 < len(self._cstart) else self._size - len(EOF_BLOCK)
        blk = os.pread(self._fd, cend - cstart, cstart)
        data = zlib.decompress(blk[18:-8], -15)
        self.blocks_read += 1
        self.compressed_bytes_read += len(blk)
        self.uncompressed_bytes_decoded += len(data)
        if self._cache_max:
            self._cache[b] = data
            if len(self._cache) > self._cache_max:
                self._cache.popitem(last=False)
        return data

    def fetch(self, name: str, start: int, end: int) -> bytes:
        e = self.fai[name]
        if not 0 <= start <= end <= e.length:
            raise ValueError(f"slice {start}:{end} outside 0:{e.length}")
        if start == end:
            return b""
        ulo, uhi = e.byte_of(start), e.byte_of(end - 1) + 1
        blocks = self._blocks_for(ulo, uhi)
        with self._lock:
            data = b"".join(self._decoded(b) for b in blocks)
        base = self._ustart[blocks[0]]
        return data[ulo - base : uhi - base].replace(b"\n", b"")
