"""Content identifiers (planning/data_layer_and_tests.md §2)."""

from __future__ import annotations

import hashlib

from blake3 import blake3


def contig_uid(seq: bytes) -> bytes:
    """BLAKE3-128 of the canonical uppercase sequence."""
    return blake3(seq.upper()).digest(length=16)


def protein_uid(aa: str) -> bytes:
    """BLAKE3-128 of the translated amino-acid sequence (no stop character)."""
    return blake3(aa.encode("ascii")).digest(length=16)


def genome_uid(contig_uids: list[bytes]) -> bytes:
    """BLAKE3-128 over the sorted contig UIDs, so case variants give the same genome."""
    return blake3(b"".join(sorted(contig_uids))).digest(length=16)


def region_uid(contig_uid: bytes, start: int, end: int, strand: int) -> bytes:
    """BLAKE3-128 of (contig_uid, start, end, strand): a nucleotide search target (IGR, dark ORF)."""
    return blake3(contig_uid + start.to_bytes(8, "big") + end.to_bytes(8, "big")
                  + strand.to_bytes(1, "big", signed=True)).digest(length=16)


KIND: dict[str, int] = {"cds": 1, "ncrna": 2, "trna": 3, "crispr_array": 4, "other": 15}


def feature_uid(contig_uid: bytes, start: int, end: int, strand: int, kind: int) -> bytes:
    """BLAKE3-128 of (contig_uid, start, end, strand, kind): the stable key of a feature.

    An origin-wrapping feature on a circular contig has end > contig_len; that is allowed.
    """
    if len(contig_uid) != 16:
        raise ValueError(f"contig_uid must be 16 bytes, got {len(contig_uid)}")
    if start < 0 or end <= start:
        raise ValueError(f"need 0 <= start < end, got start={start} end={end}")
    if strand not in (1, -1):
        raise ValueError(f"strand must be +1 or -1, got {strand}")
    if kind not in KIND.values():
        raise ValueError(f"kind must be one of {sorted(KIND.values())}, got {kind}")
    return blake3(contig_uid + start.to_bytes(8, "big") + end.to_bytes(8, "big")
                  + strand.to_bytes(1, "big", signed=True) + kind.to_bytes(1, "big")).digest(length=16)


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
