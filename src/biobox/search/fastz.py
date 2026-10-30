"""zstd-compressed FASTA for search-pack shards (`part-<k>.faa.zst`, `part-<k>.fna.zst`).

Records are `>name` plus the sequence on one line. Shards are standard zstd frames
(pyarrow's codec), so `zstd -dc` reads them.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Iterator

import pyarrow as pa


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
    with pa.CompressedInputStream(pa.OSFile(str(path)), "zstd") as fh:
        data = fh.read().decode("ascii")
    name = None
    for line in data.splitlines():
        if line.startswith(">"):
            name = line[1:]
        elif name is not None:
            yield name, line
            name = None
