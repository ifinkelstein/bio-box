"""Contig tapes, coordinates only (operation O-TAPE; PLAN.md §4.2, data_layer_and_tests.md §1 L3).

One row per contig of a stage (every default variant in `contig_index`, as for IGRs), with
parallel lists holding one entry per gene slot; the list index is `slot_ordinal`, so
`(contig_uid, slot_ordinal)` is the slot key everywhere. Tapes are append-only per stage and
carry stable keys only (`protein_uid`); tokens are a release-time join (`biobox.l3.tokens`).

S0 fills the columns L1 can supply: kinds, strands, coordinates, the junction to the
previous slot, censoring distances and protein UIDs. Columns that need later sidecars
(segments, directons, regimes, taxonomy, structure) are added as nullable columns when
their producers exist (LY-09).

Junctions and distances on a linear contig: the first slot has no previous junction (null),
and `orfs_to_contig_starts/ends` count the slots between a slot and the contig ends,
saturating at 65,535 with `slot_flags` bit 0 set (C-07). On a circular contig the first
slot's previous slot is the last one across the origin, and the distance lists are null
(no edge).

Slots are in start order (slot_ordinal), so where genes nest or share a start, reverse
complementing a contig reorders them and their slot offsets change; MR-03 holds exactly for
slots whose start and end orders agree (tests/r1/test_l3.py).

Writes `<root>/l3/tapes/gene_call_version=<v>/stage=<S>/part-<k>.parquet`, sorted by
(tax_order_bucket, contig_uid), cut at 1 GB (estimated).
"""

from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from biobox import layout, ledger, provenance
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.l3.geometry import pair_geometry
from biobox.registry.schemas import table

KIND_CDS = 1
U16_MAX = 65_535
FLAG_DIST_CLIPPED = 1
ROW_GROUP_CONTIGS = 8_192  # data_layer_and_tests.md §3.2


def build_tape(slots: list[tuple[int, int, int]], contig_len: int, circular: bool) -> dict:
    """List columns for one contig from its slots (start, end, strand) in slot_ordinal order."""
    n = len(slots)
    delta: list[int | None] = []
    junc: list[int | None] = []
    for i in range(n):
        if i == 0 and (not circular or n == 1):
            delta.append(None)
            junc.append(None)
            continue
        t, d, _ = pair_geometry(slots, i, (i - 1) % n, -1, contig_len, circular)
        junc.append(t)
        delta.append(d)
    flags = [0] * n
    to_start = to_end = None
    if not circular:
        to_start, to_end = [], []
        for i in range(n):
            a, b = i, n - 1 - i
            if a > U16_MAX or b > U16_MAX:
                flags[i] |= FLAG_DIST_CLIPPED
            to_start.append(min(a, U16_MAX))
            to_end.append(min(b, U16_MAX))
    return {
        "feature_kinds": [KIND_CDS] * n,
        "strands": [s[2] for s in slots],
        "start_bps": [s[0] for s in slots],
        "end_bps": [s[1] for s in slots],
        "delta_bp_prev_slots": delta,
        "junction_topology_prevs": junc,
        "orfs_to_contig_starts": to_start,
        "orfs_to_contig_ends": to_end,
        "slot_flags": flags,
    }


def _occurrence(root: Path, gene_call_version: str, stage: str) -> list[dict]:
    p = root / f"l1/protein_occurrence/gene_call_version={gene_call_version}/stage={stage}/part-0.parquet"
    if not p.exists():
        return []
    cols = ["protein_uid", "contig_uid", "slot_ordinal", "genome_uid", "start", "end", "strand", "tier",
            "assembly_quality"]
    return pq.read_table(p, columns=cols).to_pylist()


def run(root: Path, stage: str, circular: set[bytes] | frozenset = frozenset(),
        gene_call_version: str = GENE_CALL_VERSION, max_file_bytes: int = layout.MAX_FILE_BYTES) -> dict:
    """Write the stage's coordinate tapes from protein_occurrence and contig_index."""
    root = Path(root)
    idx = pq.read_table(root / f"l0seq/contig_index/stage={stage}/part-0.parquet",
                        columns=["contig_uid", "contig_len", "is_default_variant"]).to_pylist()
    lengths = {r["contig_uid"]: r["contig_len"] for r in idx if r["is_default_variant"]}
    by_contig: dict[bytes, list[dict]] = {cu: [] for cu in lengths}
    for r in _occurrence(root, gene_call_version, stage):
        by_contig.setdefault(r["contig_uid"], []).append(r)
    rows, ledger_rows = [], []
    for cu in sorted(by_contig):
        occ = sorted(by_contig[cu], key=lambda r: r["slot_ordinal"])
        if [r["slot_ordinal"] for r in occ] != list(range(len(occ))):
            raise ValueError(f"{cu.hex()}: slot ordinals are not 0..n-1")
        slots = [(r["start"], r["end"], r["strand"]) for r in occ]
        first = occ[0] if occ else {}
        rows.append({
            "contig_uid": cu,
            "tax_order_bucket": None,
            "genome_uid": first.get("genome_uid"),
            "contig_len": lengths[cu],
            "is_circular": cu in circular,
            "n_slots": len(occ),
            "tier": first.get("tier"),
            "assembly_quality": first.get("assembly_quality"),
            **build_tape(slots, lengths[cu], cu in circular),
            "protein_uids": [r["protein_uid"] for r in occ],
        })
        ledger_rows.append(ledger.row(cu, "contig", len(occ), len(occ), tools=provenance.tool_versions("pyarrow")))
    spec = table("contig_tapes")
    est = [200 + 30 * r["n_slots"] for r in rows]
    starts = layout.cut(est, max_file_bytes)
    base = root / f"l3/tapes/gene_call_version={gene_call_version}/stage={stage}"
    for k, (a, b) in enumerate(zip(starts, starts[1:] + [len(rows)])):
        base.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.Table.from_pylist(rows[a:b], schema=spec.schema), base / f"part-{k}.parquet",
                       compression="zstd", row_group_size=ROW_GROUP_CONTIGS,
                       write_statistics=["contig_uid", "tax_order_bucket", "n_slots", "contig_len"])
    ledger.write(root, stage, "tape", ledger_rows)
    return {"contigs": len(rows), "slots": sum(r["n_slots"] for r in rows), "files": len(starts)}
