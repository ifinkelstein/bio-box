"""bgzip FASTA export (operation O-EXTRACT, data_layer_and_tests.md §2A.1).

The store's shards are already bgzip FASTA named by `raw_digest`; `extract` writes a
chosen set of contigs (default variants) named by `contig_uid` for a tool run.
"""

from __future__ import annotations

from pathlib import Path

from biobox.l0seq.bgzf import EOF_BLOCK, BgzfFasta, FaiEntry, write_bgzf_fasta  # noqa: F401


def extract(store, contigs, out: Path, line_width: int = 80, level: int = 6) -> dict:
    """`bio-box seq extract --bgzip --faidx`: default variant of each contig, named by contig_uid hex."""
    return write_bgzf_fasta(((cu.hex(), store.get(cu)) for cu in contigs), out, line_width, level)
