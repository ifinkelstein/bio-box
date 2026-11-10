"""Slot pairs within +-k ORFs (operation O-PAIRS; PLAN.md §4.2 `slot_pairs`).

For every slot with a token at a level, one row per other slot within `K_MAX` slots that has
a token at the same level, in both directions (anchor -> partner and back). Queries narrow
to their own k <= K_MAX with `abs(orf_offset) <= k`. On a linear contig the window stops at
the contig ends and an anchor whose window is cut is flagged; on a circular contig offsets
are the shortest way round (a partner exactly half way round takes the way with the smaller
bp gap, forward on a tie), and a pair whose path crosses the origin is flagged. Until
regulon loci exist the locus is the whole contig.

Columns beyond the plan's minimum that S0 cannot fill (tax_order, vehicle, BioProject,
partner badge) are nullable or come later (LY-09).

Writes `<root>/l3/slot_pairs/release=<R>/anchor_bucket=<all|h|hh>/part-<k>.parquet`, sorted
by (anchor_token, partner_token, contig_uid, anchor_ordinal, partner_ordinal). The bucket
is the anchor token's low byte as hex at the width `biobox.layout` picks for the release's
size, so a lookup by token opens one bucket.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from biobox import layout
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.l3.geometry import cyclic_offset, pair_geometry
from biobox.l3.tokens import LEVELS, read_tapes, tape_version_of
from biobox.registry.schemas import table

K_MAX = 10
FLAG_CENSORED = 1  # the anchor's +-K_MAX window is cut by a linear contig end
FLAG_CROSSES_ORIGIN = 2
BYTES_PER_ROW = 20  # PLAN.md §4.2 size model: 15-25 B per row


def bucket_key(tid: int) -> bytes:
    return bytes([tid & 0xFF])


def window(i: int, n: int, circular: bool, k: int) -> list[tuple[int, int]]:
    """(partner ordinal, signed offset) for the slots within k of slot i, nearest first."""
    out = []
    for d in sorted(range(-k, k + 1), key=lambda d: (abs(d), d)):
        if d == 0:
            continue
        if circular:
            j = (i + d) % n
            if j == i or cyclic_offset(i, j, n) != d:
                continue
        else:
            j = i + d
            if not 0 <= j < n:
                continue
        out.append((j, d))
    return out


def pairs_of_contig(tape: dict, tokens: dict[str, list[int | None]], k: int = K_MAX) -> list[dict]:
    slots = list(zip(tape["start_bps"], tape["end_bps"], tape["strands"]))
    n, circ, length = len(slots), tape["is_circular"], tape["contig_len"]
    rows = []
    for i in range(n):
        flags_i = FLAG_CENSORED if not circ and (i < k or n - 1 - i < k) else 0
        for j, d in window(i, n, circ, k):
            topo, delta, crosses = pair_geometry(slots, i, j, d, length, circ)
            if circ and 2 * abs(d) == n:
                # equally near both ways round: take the smaller gap (reverse complement keeps it, MR-03)
                other = pair_geometry(slots, i, j, -d, length, circ)
                if other[1] < delta:
                    d, (topo, delta, crosses) = -d, other
            for lvl in LEVELS:
                a, b = tokens[lvl][i], tokens[lvl][j]
                if a is None or b is None:
                    continue
                rows.append({"anchor_token": a, "partner_token": b, "contig_uid": tape["contig_uid"],
                             "anchor_ordinal": i, "partner_ordinal": j, "orf_offset": d,
                             "junction_topology": topo, "delta_bp": delta,
                             "flags": flags_i | (FLAG_CROSSES_ORIGIN if crosses else 0),
                             "tier": tape["tier"], "tax_order": None})
    return rows


def run(root: Path, release: str, gene_call_version: str = GENE_CALL_VERSION, k: int = K_MAX,
        prefix_width: int | None = None, max_file_bytes: int = layout.MAX_FILE_BYTES,
        producers: list[str] | None = None) -> dict:
    """`producers`: the tape version's producers, by default the gene-call producer alone (card B-07)."""
    root = Path(root)
    tv = tape_version_of(gene_call_version, producers)
    tapes = {t["contig_uid"]: t for t in read_tapes(root, gene_call_version, producers=producers)}
    tok_files = sorted((root / f"l3/tape_tokens/tape_version={tv}/release={release}").glob("part-*.parquet"))
    rows = []
    for f in tok_files:
        for tok in pq.read_table(f).to_pylist():
            t = tapes[tok["contig_uid"]]
            rows += pairs_of_contig(t, {lvl: tok[f"{lvl}_tokens"] for lvl in LEVELS}, k)
    width = layout.prefix_width(BYTES_PER_ROW * len(rows)) if prefix_width is None else prefix_width
    buckets: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        buckets[layout.bucket(bucket_key(r["anchor_token"]), width)].append(r)
    spec = table("slot_pairs")
    n_files = 0
    for b in sorted(buckets):
        rs = sorted(buckets[b], key=lambda r: (r["anchor_token"], r["partner_token"], r["contig_uid"],
                                               r["anchor_ordinal"], r["partner_ordinal"]))
        starts = layout.cut([BYTES_PER_ROW] * len(rs), max_file_bytes)
        base = root / f"l3/slot_pairs/release={release}/anchor_bucket={b}"
        base.mkdir(parents=True, exist_ok=True)
        for part, (a, z) in enumerate(zip(starts, starts[1:] + [len(rs)])):
            pq.write_table(pa.Table.from_pylist(rs[a:z], schema=spec.schema), base / f"part-{part}.parquet",
                           compression="zstd", row_group_size=122_880,
                           write_statistics=["anchor_token", "partner_token", "contig_uid", "anchor_ordinal",
                                             "partner_ordinal"])
            n_files += 1
    return {"pairs": len(rows), "buckets": len(buckets), "files": n_files, "k": k}
