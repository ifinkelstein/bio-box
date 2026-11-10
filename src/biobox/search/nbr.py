"""Neighborhood fetch with explicit semantics (operation O-NBR; SEARCH_PLAN SD-08; S-09).

For each anchor slot, the slots in its window, with the anchor itself at offset 0:

    window      `10g`: slots within 10 slot offsets; `20kb` / `500bp`: slots whose gap to the
                anchor (max(starts) - min(ends), negative for overlaps) is at most that many bp
    strand      `any`, or `same`: only slots on the anchor's strand
    operon_gap  bp; when set, only slots joined to the anchor by a run of same-strand slots
                whose consecutive gaps are all <= operon_gap (`in_operon`)
    edge        linear contigs only. `truncate`: keep the anchor, cut its window at the
                contig end and set flag bit 0 when the cut removes window positions;
                `exclude:<bp>`: drop anchors within <bp> of a contig end (e.g. exclude:3kb);
                `require:<n>`: drop anchors with fewer than n slots on either side
    circular    `auto`: wrap on contigs whose tape says circular; `off`: treat every contig as linear

On a circular contig a slot is reached the shorter way round (by slot offset for gene
windows, forward on a tie; by gap for bp windows, forward on a tie), and flag bit 1 marks a
neighbor reached across the origin. `offset` is along the contig; `offset_rel` is along the
anchor's strand, positive downstream of it.
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.l3.geometry import cyclic_offset, pair_geometry
from biobox.l3.tapes import tape_version
from biobox.registry.schemas import table

FLAG_CENSORED = 1
FLAG_CROSSES_ORIGIN = 2
_BP = re.compile(r"^(\d+(?:\.\d+)?)(bp|kb|mb)$")
_UNITS = {"bp": 1, "kb": 1_000, "mb": 1_000_000}


class PolicyError(ValueError):
    pass


def parse_bp(text: str) -> int:
    m = _BP.match(str(text).strip().lower())
    if not m:
        raise PolicyError(f"{text!r}: expected a length like 500bp, 3kb")
    return int(float(m[1]) * _UNITS[m[2]])


@dataclass(frozen=True)
class Policy:
    window_genes: int | None = 10
    window_bp: int | None = None
    strand: str = "any"
    operon_gap: int | None = None
    edge: str = "truncate"  # truncate | exclude | require
    edge_bp: int = 0
    edge_n: int = 0
    circular: str = "auto"

    @classmethod
    def parse(cls, window: str = "10g", strand: str = "any", operon_gap: int | str | None = None,
              edge: str = "truncate", circular: str = "auto") -> "Policy":
        w = str(window).strip().lower()
        wg = wb = None
        if re.fullmatch(r"\d+g", w):
            wg = int(w[:-1])
        else:
            wb = parse_bp(w)
        if strand not in ("any", "same"):
            raise PolicyError(f"strand {strand!r}: expected any or same")
        if circular not in ("auto", "off"):
            raise PolicyError(f"circular {circular!r}: expected auto or off")
        og = None if operon_gap is None else (int(operon_gap) if str(operon_gap).isdigit() else parse_bp(operon_gap))
        e, ebp, en = edge, 0, 0
        if edge.startswith("exclude:"):
            e, ebp = "exclude", parse_bp(edge.split(":", 1)[1])
        elif edge.startswith("require:"):
            e, en = "require", int(edge.split(":", 1)[1])
        elif edge != "truncate":
            raise PolicyError(f"edge {edge!r}: expected truncate, exclude:<bp> or require:<n>")
        return cls(wg, wb, strand, og, e, ebp, en, circular)

    def hash(self) -> str:
        return hashlib.sha256(json.dumps(asdict(self), sort_keys=True).encode()).hexdigest()[:16]


def _operon(slots, i: int, n: int, circ: bool, gap: int, length: int) -> set[int]:
    """Slots joined to i by same-strand runs with every consecutive gap <= gap."""
    members = {i}
    strand = slots[i][2]
    for step in (1, -1):
        cur = i
        for _ in range(n - 1):
            nxt = cur + step
            if circ:
                nxt %= n
            elif not 0 <= nxt < n:
                break
            if nxt in members or slots[nxt][2] != strand:
                break
            _, d, _ = pair_geometry(slots, cur, nxt, step, length, circ)
            if d > gap:
                break
            members.add(nxt)
            cur = nxt
    return members


def _candidates(slots, i: int, n: int, circ: bool, p: Policy, length: int):
    """(j, offset, topo, delta, crosses) for every slot in i's window, before strand/operon filters."""
    out = []
    if p.window_genes is not None:
        k = p.window_genes
        for j in range(n):
            if j == i:
                continue
            d = cyclic_offset(i, j, n) if circ else j - i
            if abs(d) <= k:
                out.append((j, d, *pair_geometry(slots, i, j, d, length, circ)))
        return out
    for j in range(n):
        if j == i:
            continue
        if circ:
            fwd, back = (j - i) % n, -((i - j) % n)
            a = pair_geometry(slots, i, j, fwd, length, circ)
            b = pair_geometry(slots, i, j, back, length, circ)
            d, g = (back, b) if b[1] < a[1] else (fwd, a)
        else:
            d = j - i
            g = pair_geometry(slots, i, j, d, length, circ)
        if g[1] <= p.window_bp:
            out.append((j, d, *g))
    return out


def neighborhood(tape: dict, i: int, policy: Policy) -> list[dict] | None:
    """Rows for anchor slot i of one tape; None when the edge policy drops the anchor."""
    slots = list(zip(tape["start_bps"], tape["end_bps"], tape["strands"]))
    n, length = len(slots), tape["contig_len"]
    circ = bool(tape["is_circular"]) and policy.circular == "auto"
    s, e, strand = slots[i]
    censored = False
    if not circ:
        if policy.edge == "exclude" and min(s, length - e) < policy.edge_bp:
            return None
        if policy.edge == "require" and (i < policy.edge_n or n - 1 - i < policy.edge_n):
            return None
        if policy.window_genes is not None:
            censored = i < policy.window_genes or n - 1 - i < policy.window_genes
        else:
            censored = s - policy.window_bp < 0 or e + policy.window_bp > length
    operon = _operon(slots, i, n, circ, policy.operon_gap, length) if policy.operon_gap is not None else None
    base = FLAG_CENSORED if censored else 0
    rows = [_row(tape, i, i, 0, strand, None, None, operon, base, policy)]
    for j, d, topo, delta, crosses in _candidates(slots, i, n, circ, policy, length):
        if policy.strand == "same" and slots[j][2] != strand:
            continue
        if operon is not None and j not in operon:
            continue
        rows.append(_row(tape, i, j, d, strand, topo, delta, operon,
                         base | (FLAG_CROSSES_ORIGIN if crosses else 0), policy))
    rows.sort(key=lambda r: r["offset"])
    return rows


def _row(tape, i, j, d, anchor_strand, topo, delta, operon, flags, policy) -> dict:
    return {"anchor_protein_uid": tape["protein_uids"][i], "contig_uid": tape["contig_uid"], "anchor_ordinal": i,
            "ordinal": j, "offset": d, "offset_rel": d * anchor_strand, "protein_uid": tape["protein_uids"][j],
            "start": tape["start_bps"][j], "end": tape["end_bps"][j], "strand": tape["strands"][j],
            "delta_bp": delta, "junction_topology": topo, "same_strand": tape["strands"][j] == anchor_strand,
            "in_operon": None if operon is None else j in operon, "flags": flags, "policy_hash": policy.hash()}


def _tapes(root: Path, gene_call_version: str, contigs: set[bytes]) -> dict[bytes, dict]:
    out = {}
    tv = tape_version([gene_call_version])  # the gene-call producer's tapes (card B-07)
    for f in sorted((Path(root) / "l3/tapes").glob(f"tape_version={tv}/stage=*/part-*.parquet")):
        t = pq.read_table(f, filters=[("contig_uid", "in", sorted(contigs))]) if contigs else None
        for r in (t.to_pylist() if t is not None else []):
            out[r["contig_uid"]] = r
    return out


def occurrences(root: Path, protein_uids: set[bytes], gene_call_version: str = GENE_CALL_VERSION) -> list[dict]:
    files = sorted((Path(root) / "l1/protein_occurrence").glob(f"gene_call_version={gene_call_version}/stage=*/part-*.parquet"))
    if not protein_uids:
        return []
    return [r for f in files for r in pq.read_table(f, filters=[("protein_uid", "in", sorted(protein_uids))]).to_pylist()]


def fetch(root: Path, protein_uids: set[bytes], policy: Policy,
          gene_call_version: str = GENE_CALL_VERSION) -> tuple[list[dict], dict]:
    """Neighborhoods of every occurrence of the anchor proteins; (rows, counts)."""
    occ = occurrences(root, protein_uids, gene_call_version)
    tapes = _tapes(root, gene_call_version, {r["contig_uid"] for r in occ})
    rows, dropped = [], 0
    for r in sorted(occ, key=lambda r: (r["contig_uid"], r["slot_ordinal"])):
        got = neighborhood(tapes[r["contig_uid"]], r["slot_ordinal"], policy)
        if got is None:
            dropped += 1
        else:
            rows += got
    rows.sort(key=lambda r: (r["anchor_protein_uid"], r["contig_uid"], r["anchor_ordinal"], r["offset"]))
    return rows, {"anchors": len(occ), "anchors_dropped_by_edge": dropped, "rows": len(rows)}


def write(rows: list[dict], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=table("neighborhoods").schema), path, compression="zstd")
    return path
