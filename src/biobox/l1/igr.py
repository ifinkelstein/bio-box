"""Intergenic regions as coordinates (data_layer_and_tests.md §2A.1; R-05).

An IGR is a maximal run of bases covered by no slot of one gene-calling version, so
overlapping and nested genes produce none between them. Every run of >= 1 nt is stored;
the >= 30 nt filter (Weinberg convention, SEARCH_PLAN SD-07) is applied when the search
pack is built. Coordinates are 0-based half-open on the + strand (`strand` = 1);
orientation is derived from the flanking slots. On linear contigs the end runs have a
null flank; on circular contigs the run across the origin has `end > contig_len`
(C-05), and a contig with no genes is one IGR with both flanks null.

Writes `<root>/l0seq/igr/igr_version=<v>/stage=<S>/part-0.parquet`.
"""

from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.registry.schemas import table

IGR_RULE = "igr1"


def igr_version(gene_call_version: str = GENE_CALL_VERSION) -> str:
    return f"{gene_call_version}-{IGR_RULE}"


def igrs_of_contig(slots: list[tuple[int, int, int]], length: int, circular: bool) -> list[tuple[int, int, int | None, int | None]]:
    """(start, end, left_slot, right_slot) for one contig; slots are (slot_ordinal, start, end)."""
    if not slots:
        return [(0, length, None, None)] if length else []
    genes = sorted(slots, key=lambda s: (s[1], s[2], s[0]))
    out = []
    wraps = [(e - length, o) for o, s, e in genes if e > length] if circular else []
    if wraps:
        # an origin-wrapping gene also covers [0, end - length); scan from its tail, and the
        # run before the first gene closes the circle, so no closing run is added
        reach, reach_slot = max(wraps)
        rest = genes
    else:
        reach, reach_slot = genes[0][2], genes[0][0]
        rest = genes[1:]
        if not circular and genes[0][1] > 0:
            out.append((0, genes[0][1], None, genes[0][0]))
    for ordinal, s, e in rest:
        if s > reach:
            out.append((reach, s, reach_slot, ordinal))
        if e > reach:
            reach, reach_slot = e, ordinal
    if circular and not wraps:
        first_ordinal, first_start = genes[0][0], genes[0][1]
        if reach < first_start + length:
            start, end = reach, first_start + length if first_start else length
            if start >= length:  # the last gene ends exactly at the origin
                start, end = start - length, end - length
            out.append((start, end, reach_slot, first_ordinal))
    elif not circular and reach < length:
        out.append((reach, length, reach_slot, None))
    return sorted(out)


def run(root: Path, stage: str, circular: set[bytes] | frozenset = frozenset(),
        gene_call_version: str = GENE_CALL_VERSION) -> dict:
    """Build IGRs for every contig of a stage from protein_occurrence and contig_index."""
    root = Path(root)
    occ = pq.read_table(root / f"l1/protein_occurrence/gene_call_version={gene_call_version}/stage={stage}/part-0.parquet",
                        columns=["contig_uid", "slot_ordinal", "start", "end"]).to_pylist()
    idx = pq.read_table(root / f"l0seq/contig_index/stage={stage}/part-0.parquet").to_pylist()
    lengths = {r["contig_uid"]: r["contig_len"] for r in idx if r["is_default_variant"]}
    by_contig: dict[bytes, list] = {cu: [] for cu in lengths}
    for r in occ:
        by_contig.setdefault(r["contig_uid"], []).append((r["slot_ordinal"], r["start"], r["end"]))
    rows = []
    for cu in sorted(by_contig):
        for s, e, left, right in igrs_of_contig(by_contig[cu], lengths[cu], cu in circular):
            rows.append({"contig_uid": cu, "start": s, "end": e, "strand": 1, "left_slot": left, "right_slot": right})
    v = igr_version(gene_call_version)
    out = root / f"l0seq/igr/igr_version={v}/stage={stage}/part-0.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=table("igr").schema), out, compression="zstd", row_group_size=122_880)
    return {"igrs": len(rows), "igrs_ge_30": sum(1 for r in rows if r["end"] - r["start"] >= 30),
            "bp": sum(r["end"] - r["start"] for r in rows)}
