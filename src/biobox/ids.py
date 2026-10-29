"""Content identifiers (planning/data_layer_and_tests.md §2)."""

from __future__ import annotations

import hashlib

from blake3 import blake3


def contig_uid(seq: bytes) -> bytes:
    """BLAKE3-128 of the canonical uppercase sequence."""
    return blake3(seq.upper()).digest(length=16)


def raw_digest(seq: bytes) -> bytes:
    """sha256 of the exact stored bytes (case and IUPAC preserved)."""
    return hashlib.sha256(seq).digest()


def source_uid(source_name: str, source_release: str) -> int:
    """First 8 bytes of sha256(name, release) as uint64."""
    h = hashlib.sha256(f"{source_name}\0{source_release}".encode()).digest()
    return int.from_bytes(h[:8], "big")


def sha256_file(path) -> bytes:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.digest()
