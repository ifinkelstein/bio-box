"""bgzip FASTA export with `.fai` and `.gzi` (operation O-EXTRACT, data_layer_and_tests.md §2A.1).

`write_bgzf_fasta` writes the BGZF layout samtools/htslib read: gzip members of at most
65,280 uncompressed bytes with the `BC` extra field, a 28-byte EOF block, a `.fai`
(name, length, offset, line bases, line width) and a `.gzi` (block starts as
(compressed, uncompressed) uint64 pairs, first block omitted). Output is deterministic
(MTIME 0, fixed level).

`BgzfFasta` reads a slice back through the two indexes and reports which blocks it
decompressed, so the K-13 bake-off can count BGZF blocks the way LY-06 counts row groups.
"""

from __future__ import annotations

import bisect
import struct
import zlib
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
    """Write (name, seq bytes) records as bgzip FASTA plus `.fai` and `.gzi`. Returns sizes."""
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fai: list[FaiEntry] = []
    gzi: list[tuple[int, int]] = []
    buf = bytearray()
    upos = cpos = 0  # uncompressed bytes emitted so far, compressed bytes written so far

    with open(out, "wb") as fh:
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

    Path(str(out) + ".fai").write_text(
        "".join(f"{e.name}\t{e.length}\t{e.offset}\t{e.line_bases}\t{e.line_width}\n" for e in fai))
    with open(str(out) + ".gzi", "wb") as fh:
        fh.write(struct.pack("<Q", len(gzi)))
        for c, u in gzi:
            fh.write(struct.pack("<QQ", c, u))
    return {"bytes": cpos, "uncompressed_bytes": upos, "blocks": len(gzi) + 1, "records": len(fai)}


class BgzfFasta:
    """Slice a bgzip FASTA through its `.fai` and `.gzi`; counts blocks decompressed."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.fai: dict[str, FaiEntry] = {}
        for line in Path(str(path) + ".fai").read_text().splitlines():
            n, length, off, lb, lw = line.split("\t")[:5]
            self.fai[n] = FaiEntry(n, int(length), int(off), int(lb), int(lw))
        raw = Path(str(path) + ".gzi").read_bytes()
        (n,) = struct.unpack_from("<Q", raw)
        pairs = [struct.unpack_from("<QQ", raw, 8 + 16 * i) for i in range(n)]
        self._cstart = [0] + [c for c, _ in pairs]
        self._ustart = [0] + [u for _, u in pairs]
        self._size = self.path.stat().st_size
        self.blocks_read = 0
        self.compressed_bytes_read = 0
        self.uncompressed_bytes_decoded = 0

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

    def fetch(self, name: str, start: int, end: int) -> bytes:
        e = self.fai[name]
        if not 0 <= start <= end <= e.length:
            raise ValueError(f"slice {start}:{end} outside 0:{e.length}")
        if start == end:
            return b""
        ulo, uhi = e.byte_of(start), e.byte_of(end - 1) + 1
        blocks = self._blocks_for(ulo, uhi)
        data = bytearray()
        with open(self.path, "rb") as fh:
            for b in blocks:
                cstart = self._cstart[b]
                cend = self._cstart[b + 1] if b + 1 < len(self._cstart) else self._size - len(EOF_BLOCK)
                fh.seek(cstart)
                blk = fh.read(cend - cstart)
                data += zlib.decompress(blk[18:-8], -15)
                self.blocks_read += 1
                self.compressed_bytes_read += len(blk)
        self.uncompressed_bytes_decoded += len(data)
        text = bytes(data[ulo - self._ustart[blocks[0]] : uhi - self._ustart[blocks[0]]])
        return text.replace(b"\n", b"")


def extract(store, contigs, out: Path, line_width: int = 80, level: int = 6) -> dict:
    """`bio-box seq extract --bgzip --faidx`: default variant of each contig, named by contig_uid hex."""
    return write_bgzf_fasta(((cu.hex(), store.get(cu)) for cu in contigs), out, line_width, level)
