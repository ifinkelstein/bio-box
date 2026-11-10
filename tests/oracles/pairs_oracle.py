"""The slot-pair oracle: `pairs_of_contig`, `window`, the flags and `K_MAX`, moved unchanged from
`biobox.l3.pairs` when card B-11 removed the materialized `slot_pairs` table (AGENTS §3.5).

The catalog's pair macros (`pairs_of`, `q1_count`, `partners_of`) answer neighbour questions from
`token_postings` at query time; the tests compare them with this loop over the tapes.

For every slot with a token at a level, one row per other slot within `k` slots that has a token at
the same level, in both directions (anchor -> partner and back). On a linear contig the window stops
at the contig ends and an anchor whose window is cut is flagged; on a circular contig offsets are the
shortest way round (a partner exactly half way round takes the way with the smaller bp gap, forward
on a tie), and a pair whose path crosses the origin is flagged.
"""

from __future__ import annotations

from biobox.l3.geometry import cyclic_offset, pair_geometry
from biobox.l3.tokens import LEVELS

K_MAX = 10
FLAG_CENSORED = 1  # the anchor's +-K_MAX window is cut by a linear contig end
FLAG_CROSSES_ORIGIN = 2


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
