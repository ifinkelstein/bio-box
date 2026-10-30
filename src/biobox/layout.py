"""File sizing for prefix-bucketed tables (LY-04, data_layer_and_tests.md §3.2).

A table bucketed by key prefix (`protein_seq` by `protein_prefix`) picks the bucket width
per stage from the stage's size, so a small stage is one file and a large one is spread
over 16 or 256 buckets of at least `TARGET_FILE_BYTES`:

    width 0   protein_prefix=all   the whole stage
    width 1   protein_prefix=<h>   first hex digit of the key
    width 2   protein_prefix=<hh>  first byte of the key

Within a bucket, rows (sorted by key) are cut into `part-<k>.parquet` files of at most
`MAX_FILE_BYTES` (estimated). Both choices depend only on the stage's content (D-01, D-04).
A reader finds a key's bucket in every stage with `protein_prefix IN ('all', h, hh)`.
"""

from __future__ import annotations

TARGET_FILE_BYTES = 256 * 2**20
MAX_FILE_BYTES = 2**30
ALL = "all"


def prefix_width(est_bytes: float, target: int = TARGET_FILE_BYTES) -> int:
    """Widest of 2, 1, 0 hex digits whose buckets still average at least `target` bytes."""
    for width in (2, 1):
        if est_bytes / 16**width >= target:
            return width
    return 0


def bucket(key: bytes, width: int) -> str:
    return ALL if width == 0 else key[:1].hex()[:width]


def bucket_matches(value: str, key: bytes) -> bool:
    """LY-05: a row with `key` may live under partition `value`."""
    return value == ALL or (len(value) in (1, 2) and key[:1].hex().startswith(value))


def cut(sizes: list[float], max_bytes: int = MAX_FILE_BYTES) -> list[int]:
    """Start index of each part when rows of these estimated sizes are cut at `max_bytes`."""
    starts, held = [0], 0.0
    for i, n in enumerate(sizes):
        if held and held + n > max_bytes:
            starts.append(i)
            held = 0.0
        held += n
    return starts
