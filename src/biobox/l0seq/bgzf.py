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
import mmap
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


def _fai_entry(line: bytes) -> FaiEntry:
    n, length, off, lb, lw = line.split(b"\t")[:5]
    return FaiEntry(n.decode(), int(length), int(off), int(lb), int(lw))


def fai_first_unsorted(path: Path) -> int | None:
    """0-based line number of the first `.fai` line whose name is not above the previous one, or None.

    Streams the file; names must be strictly increasing (bytes order, which is str order for the
    ASCII hex names the store writes).
    """
    prev = None
    with open(path, "rb") as fh:
        for i, line in enumerate(fh):
            name = line.split(b"\t", 1)[0]
            if prev is not None and name <= prev:
                return i
            prev = name
    return None


class BgzfFasta:
    """Slice a bgzip FASTA through its `.fai` and `.gzi`.

    The `.fai` is not loaded. Its records are sorted by name (the store writes them sorted by
    `raw_digest` hex), so `entry(name)` bisects the memory-mapped file on line starts and parses
    about log2(n) lines (`fai_lines_read` counts them); a small LRU keeps parsed entries. Opening
    reads no line of it. Sortedness is checked lazily: a lookup whose bisection lands between
    lines out of order raises; `check_sorted()` reads the whole file once (contracts). An export
    whose records are in caller order (`extract`) is opened with `name_sorted=False`, which reads
    the whole `.fai` into a dict, as the export's caller chose its size.

    The `.gzi` block table is loaded whole: 16 bytes on disk per 65,280-byte block, about 15,500
    blocks and 250 KB per Gbp of FASTA text, about 1.1 MB per Gbp as two Python int lists
    (*derived*: 1e9 x 81/80 / 65,280 blocks; 2 x (8-byte pointer + 28-byte int) per block).

    Counters: `blocks_read` (blocks decompressed), `compressed_bytes_read`,
    `uncompressed_bytes_decoded`, `fai_lines_read`. With `cache_blocks` > 0 a small LRU of
    decoded blocks serves repeated slices; a cache hit counts nothing. Reads use pread and the
    caches sit behind one lock, so one reader is safe to share between threads.
    """

    FAI_CACHE = 1024

    def __init__(self, path: Path, cache_blocks: int = 0, name_sorted: bool = True):
        self._fd: int | None = None
        self.path = Path(path)
        self.fai_path = Path(str(path) + ".fai")
        self._fai_all: dict[str, FaiEntry] | None = (
            None if name_sorted else {e.name: e for e in read_fai(self.fai_path)})
        self._fai_mm: mmap.mmap | None = None
        self._fai_size: int | None = None
        self._fai_cache: OrderedDict[str, FaiEntry] = OrderedDict()
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
        self.fai_lines_read = 0

    def close(self) -> None:
        if self._fd is not None:
            os.close(self._fd)
            self._fd = None
        if getattr(self, "_fai_mm", None) is not None:
            self._fai_mm.close()
            self._fai_mm = None
            self._fai_size = None

    def __del__(self):
        try:
            self.close()
        except Exception:  # at interpreter exit module globals (os) may already be gone
            pass

    def check_sorted(self) -> bool:
        """True if the `.fai` names are strictly increasing (reads the whole file once)."""
        return fai_first_unsorted(self.fai_path) is None

    def _fai_map(self) -> tuple[mmap.mmap | None, int]:
        if self._fai_size is None:
            size = self.fai_path.stat().st_size
            if size:
                with open(self.fai_path, "rb") as fh:
                    self._fai_mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
            self._fai_size = size
        return self._fai_mm, self._fai_size

    def _bisect(self, name: str) -> FaiEntry:
        mm, size = self._fai_map()
        if mm is None:
            raise KeyError(name)
        target = name.encode()
        seen: dict[int, tuple[bytes, int]] = {}  # line start -> (line without newline, next line start)

        def line_at(s: int) -> tuple[bytes, int]:
            if s not in seen:
                e = mm.find(b"\n", s)
                e = size if e < 0 else e
                self.fai_lines_read += 1
                seen[s] = (mm[s:e], e + 1)
            return seen[s]

        def name_at(s: int) -> bytes:
            return line_at(s)[0].split(b"\t", 1)[0]

        # lo is a line start and every line starting before lo has a name < target; every line
        # starting at or after hi has a name >= target. Ends with lo the first line >= target.
        lo, hi = 0, size
        while lo < hi:
            mid = (lo + hi) // 2
            if mid == 0 or mm[mid - 1] == 0x0A:
                s = mid
            else:
                s = mm.find(b"\n", mid)
                s = size if s < 0 else s + 1
            if s >= hi:  # no line starts in [mid, hi)
                hi = mid
            elif name_at(s) < target:
                lo = line_at(s)[1]
            else:
                hi = s

        def unsorted(a: bytes, b: bytes) -> ValueError:
            return ValueError(f"{self.fai_path}: the index is not sorted by name ({a.decode()!r} comes "
                              f"before {b.decode()!r}); a lookup cannot bisect it (check_sorted)")

        # lazy sortedness check on the neighbours of the landing point
        if lo > 0:
            prev = name_at(mm.rfind(b"\n", 0, lo - 1) + 1)
            if prev >= target:
                raise unsorted(prev, target)
        if lo >= size:
            raise KeyError(name)
        cur, nxt = line_at(lo)
        cname = cur.split(b"\t", 1)[0]
        if nxt < size and name_at(nxt) <= cname:
            raise unsorted(cname, name_at(nxt))
        if cname != target:
            raise KeyError(name)
        return _fai_entry(cur)

    def entry(self, name: str) -> FaiEntry:
        """The `.fai` record of `name`; KeyError if the shard has none."""
        if self._fai_all is not None:
            return self._fai_all[name]
        with self._lock:
            hit = self._fai_cache.pop(name, None)
            if hit is None:
                hit = self._bisect(name)
            self._fai_cache[name] = hit
            if len(self._fai_cache) > self.FAI_CACHE:
                self._fai_cache.popitem(last=False)
            return hit

    def _blocks_for(self, ulo: int, uhi: int) -> range:
        first = bisect.bisect_right(self._ustart, ulo) - 1
        last = bisect.bisect_right(self._ustart, uhi - 1) - 1
        return range(first, last + 1)

    def units_for(self, name: str, start: int, end: int) -> int:
        """Blocks a slice touches, from the indexes alone (no I/O)."""
        if start == end:
            return 0
        e = self.entry(name)
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
        e = self.entry(name)
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
