"""Step 4: what changes when one ncRNA slot is added and one gene is split into two segment slots.

    python step4_keys.py ROOT

Uses the built golden store (E. coli K-12 chromosome) and the repo's own igr/tape/pair functions on
the modified slot list. Rows are compared twice: under the current key (contig_uid, slot_ordinal)
and under a coordinate key (contig_uid, start, end, strand, kind).
"""

import json
import sys
from pathlib import Path

import pyarrow.parquet as pq

from biobox import io
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.l1.igr import igrs_of_contig
from biobox.l3.tapes import build_tape

# the pair enumeration lives on as the test oracle since card B-11 removed the materialized pairs
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))
from oracles.pairs_oracle import pairs_of_contig  # noqa: E402

root = Path(sys.argv[1])
occ = [r for f in io.part_files(root / f"l1/protein_occurrence/gene_call_version={GENE_CALL_VERSION}/stage=S0")
       for r in pq.read_table(f).to_pylist()]
idx = {r["contig_uid"]: r["contig_len"] for f in io.part_files(root / "l0seq/contig_index/stage=S0")
       for r in pq.read_table(f).to_pylist()}
ecoli = next(cu for cu, n in idx.items() if n == 4_641_652)
L = idx[ecoli]
base = sorted((r for r in occ if r["contig_uid"] == ecoli), key=lambda r: r["slot_ordinal"])
# features as (start, end, strand, kind, token-key); kind 1 CDS, 2 ncRNA, 3 segment
feats0 = [(r["start"], r["end"], r["strand"], 1, r["protein_uid"].hex()) for r in base]

# (a) one 377-nt ncRNA in the first IGR >= 500 nt at or after 3.27 Mb
igr0 = igrs_of_contig([(i, f[0], f[1]) for i, f in enumerate(feats0)], L, True)
s, e, _, _ = next(g for g in igr0 if g[0] >= 3_270_000 and g[1] - g[0] >= 500)
nc = (s + 60, s + 60 + 377, -1, 2, "ncRNA:new")
# (b) the longest protein split at its midpoint (in-frame) into two segment slots
longest = max(range(len(feats0)), key=lambda i: feats0[i][1] - feats0[i][0])
gs, ge, gstr, _, gkey = feats0[longest]
mid = gs + ((ge - gs) // 6) * 3
segs = [(gs, mid, gstr, 3, gkey + ":seg0"), (mid, ge, gstr, 3, gkey + ":seg1")]


def world(feats):
    feats = sorted(feats, key=lambda f: (f[0], f[1], f[2]))
    slots = [(f[0], f[1], f[2]) for f in feats]
    tape = build_tape(slots, L, True)
    tape.update(contig_uid=ecoli, is_circular=True, contig_len=L, tier=1)
    toks = {k: i for i, k in enumerate(sorted({f[4] for f in feats}))}
    t = [toks[f[4]] for f in feats]
    tk = {"family": [None] * len(feats), "nr50": [None] * len(feats), "nr90": t}  # one level; x3 at most
    pairs = pairs_of_contig(tape, tk)
    igrs = igrs_of_contig([(i, f[0], f[1]) for i, f in enumerate(feats)], L, True)
    return feats, pairs, igrs


def diff(a: set, b: set) -> dict:
    return {"before": len(a), "after": len(b), "unchanged": len(a & b), "removed_or_rewritten": len(a - b), "new": len(b - a)}


def tables(feats, pairs, igrs, coord_key: bool):
    def sk(i):  # slot key
        f = feats[i]
        return (f[0], f[1], f[2], f[3]) if coord_key else i

    slot_rows = {(sk(i), f[4]) for i, f in enumerate(feats)}  # slot_protein / protein_occurrence / tape entry
    igr_rows = {(g[0], g[1], None if g[2] is None else sk(g[2]), None if g[3] is None else sk(g[3])) for g in igrs}
    pair_rows = {(sk(p["anchor_ordinal"]), sk(p["partner_ordinal"]), p["orf_offset"], p["junction_topology"], p["delta_bp"],
                  p["flags"]) for p in pairs}
    pair_keys = {(sk(p["anchor_ordinal"]), sk(p["partner_ordinal"])) for p in pairs}
    return {"slot rows (slot_protein, protein_occurrence, tape list entries)": slot_rows, "igr rows (left_slot/right_slot)": igr_rows,
            "slot_pairs rows, one token level": pair_rows, "slot_pairs primary keys, one token level": pair_keys}


w0 = world(feats0)
out = {"contig": "E. coli K-12 MG1655 chromosome", "slots_before": len(feats0), "ncRNA_at": [nc[0], nc[1]],
       "ncRNA_ordinal": sorted(feats0 + [nc], key=lambda f: (f[0], f[1], f[2])).index(nc),
       "split_gene_ordinal": longest, "split_gene_len_bp": ge - gs}
for name, feats1 in (("add_ncRNA", feats0 + [nc]),
                     ("split_gene", [f for i, f in enumerate(feats0) if i != longest] + segs),
                     ("both", [f for i, f in enumerate(feats0) if i != longest] + segs + [nc])):
    w1 = world(feats1)
    out[name] = {}
    for label, ck in (("ordinal_key (current)", False), ("coordinate_key (proposed)", True)):
        t0, t1 = tables(*w0, ck), tables(*w1, ck)
        out[name][label] = {k: diff(t0[k], t1[k]) for k in t0}
print(json.dumps(out, indent=1))
