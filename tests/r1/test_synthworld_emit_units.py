"""Ring R1: the pure parts of synthworld's emit stage (synthworld_spec §2.1 S11, §3.1; card F-02 session 8). No ZOMBI2
here: the micro world is emitted in tests/r2/test_synthworld_emit.py.

The oracles are written in this file: a record carried through S6 insertions, S7 and S8 edits and S9 cuts as an
array of each base's S5 position (provenance), so a feature's place on a contig is found base by base, never by the
lift arithmetic under test; reverse complements spelled out; stable keys hashed here with blake3.
"""

from __future__ import annotations

import collections
import csv
import random
from pathlib import Path

import pytest
from blake3 import blake3

from biobox.synthworld import backbone, emit

pytestmark = pytest.mark.covers("F-02")  # ring R1 (by directory)
RC = {"A": "T", "C": "G", "G": "C", "T": "A"}


def rc(s: str) -> str:
    return "".join(RC[c] for c in reversed(s))


def _write(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["record"], delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def test_stop_bp_and_assembly_quality_follow_the_spec():
    # §3.1: just past the stop codon on +, the stop codon's lowest base on -
    assert emit.stop_bp(100, 400, 1, False) == 400
    assert emit.stop_bp(100, 400, -1, False) == 100
    assert emit.stop_bp(100, 400, 1, True) is None
    assert emit.stop_bp(4990, 5020, 1, False) == 5020          # across a circular origin, end > length
    assert emit.assembly_quality("complete", 1.0) == "complete"
    assert emit.assembly_quality("draft", 0.99) == "contig"
    assert emit.assembly_quality("mag", 0.93) == "mag_hq" and emit.assembly_quality("mag", 0.89) == "mag_mq"


def test_check_config_accepts_micro_and_refuses_bad_values():
    p = backbone.load_preset("micro")
    emit.check_config(p["emit"])
    for bad in ({"mislabel_rate": 1.5}, {"biomes": []}, {"biomes": ["soil", "soil"]}, {"biome_rank": "kingdom"}):
        with pytest.raises(ValueError):
            emit.check_config(dict(p["emit"], **bad))


def test_window_places_every_feature_base_where_the_window_puts_it():
    """A feature on a contig, cut by a contaminant window (forward or reverse-complemented): the parts returned
    hold exactly the feature's bases inside the window, each at the window position the cut gives it."""
    rnd = random.Random(8)
    n = 0
    for _ in range(400):
        L = rnd.randint(30, 120)
        seq = "".join(rnd.choice("ACGT") for _ in range(L))
        circ = rnd.random() < 0.5
        cs = rnd.randrange(L)
        ln = rnd.randint(1, L if circ else L - cs)
        ce = cs + ln
        st = rnd.choice((1, -1))
        p5, p3 = rnd.random() < 0.2, rnd.random() < 0.2
        x0 = rnd.randrange(L)
        y0 = rnd.randint(x0 + 1, L)
        o = rnd.choice((1, -1))
        win = seq[x0:y0] if o == 1 else rc(seq[x0:y0])
        Wn = y0 - x0
        # brute force: contig positions of the feature, then their window positions
        pos = [(cs + i) % L for i in range(ln)]
        inside = [p for p in pos if x0 <= p < y0]
        want_pos = sorted((p - x0) if o == 1 else (Wn - 1 - (p - x0)) for p in inside)
        got = emit.window(cs, ce, st, p5, p3, L, x0, y0, o)
        got_pos = sorted(q for a, b, *_ in got for q in range(a, b))
        assert got_pos == want_pos
        for a, b, wst, q5, q3 in got:
            assert wst == st * o
            # the part's bases, read on its strand, are a run of the feature's bases read on the feature's strand
            fbases = "".join(seq[p] for p in pos)
            fbases = fbases if st == 1 else rc(fbases)
            part = win[a:b] if wst == 1 else rc(win[a:b])
            assert part in fbases
            # 5' end kept iff the feature's 5' base is in the window and was not cut before
            five = pos[0] if st == 1 else pos[-1]
            three = pos[-1] if st == 1 else pos[0]
            k5 = [p for p in range(a, b)]
            five_w = (five - x0) if o == 1 else (Wn - 1 - (five - x0))
            three_w = (three - x0) if o == 1 else (Wn - 1 - (three - x0))
            assert q5 == (p5 or not (x0 <= five < y0 and five_w in k5))
            assert q3 == (p3 or not (x0 <= three < y0 and three_w in k5))
            n += 1
    assert n > 150          # parts checked (a window misses many features)


# --- the lift and placement chain against a provenance oracle ------------------------------------------------------

def _apply(seq: str, prov: list, edits):
    """Replacements [a, b) -> block (a = b an insertion; a point 0 of an insertion list means N): the new record and
    its provenance (each base's S5 position, None for a new base)."""
    out, pv, prev = [], [], 0
    for a, b, blk in sorted(edits, key=lambda x: (x[0], x[1])):
        out.append(seq[prev:a])
        pv += prov[prev:a]
        out.append(blk)
        pv += [None] * len(blk)
        prev = b
    out.append(seq[prev:])
    pv += prov[prev:]
    return "".join(out), pv


def _edits(rnd, N, k, insert_only=False):
    pts = sorted(rnd.sample(range(1, N), 2 * k))
    out = []
    for i in range(k):
        a, b = pts[2 * i], pts[2 * i + 1]
        if insert_only or rnd.random() < 0.4:
            b = a
        b = min(b, a + rnd.randint(0, 4))
        blk = "".join(rnd.choice("ACGT") for _ in range(rnd.randint(0 if a < b else 1, 5)))
        out.append((a, b, blk))
    return out


def _world(tmp: Path, rnd):
    """One circular record through S5 -> S9, the tables the Placer reads, and the oracle's arrays."""
    rec, gid = "g000_r0", "g000"
    N5 = rnd.randint(1500, 2500)
    s5 = "".join(rnd.choice("ACGT") for _ in range(N5))
    ins = [(a, a, "".join(rnd.choice("ACGT") for _ in range(rnd.randint(20, 60)))) for a in
           sorted(rnd.sample(range(1, N5 + 1), 3))]
    s6, p6 = _apply(s5, list(range(N5)), [(a, a, b) if a < N5 else (N5, N5, b) for a, _, b in ins])
    e7 = _edits(rnd, len(s6), 3)
    s7, p7 = _apply(s6, p6, e7)
    e8 = _edits(rnd, len(s7), 8)
    s8, p8 = _apply(s7, p7, e8)
    N8 = len(s8)
    cuts = sorted(k * N8 // 4 + rnd.randrange(N8 // 8) for k in range(4))      # contigs of at least N8 / 8
    spans = [(c, cuts[i + 1] if i + 1 < len(cuts) else cuts[0] + N8) for i, c in enumerate(cuts)]
    contigs, rows = {}, []
    for k, (a, b) in enumerate(spans):
        o = rnd.choice((1, -1))
        fwd = (s8 + s8)[a:b]
        s = fwd if o == 1 else rc(fwd)
        cid = f"{gid}_c{k:04d}"
        uid = blake3(s.encode()).hexdigest(length=16)
        contigs[cid] = s
        rows.append({"contig_id": cid, "genome_id": gid, "replicon": 0, "record": rec, "length": len(s),
                     "replicon_start": a, "replicon_end": b, "orientation": o, "rotation_offset": 0,
                     "is_circular": "false", "emitted": "false" if k == 2 else "true", "contig_uid": uid})
    _write(tmp / "assembly" / "contigs.tsv", rows)
    _write(tmp / "assembly" / "replicons.tsv", [{"record": rec, "genome_id": gid, "kind": "chromosome",
                                                  "topology": "circular", "length": N8}])
    _write(tmp / "replicons" / "replicons.tsv", [{"record": rec, "length": N5}])
    _write(tmp / "mobile" / "replicons.tsv", [{"record": rec, "length": len(s6)}])
    _write(tmp / "fusions" / "replicons.tsv", [{"record": rec, "length": len(s7)}])
    _write(tmp / "mobile" / "insertions.tsv", [{"record": rec, "at": a, "length": len(b)} for a, _, b in ins])
    _write(tmp / "fusions" / "edits.tsv", [{"record": rec, "start": a, "end": b, "length": len(x)} for a, b, x in e7])
    _write(tmp / "assembly" / "edits.tsv", [{"record": rec, "start": a, "end": b, "length": len(x)} for a, b, x in e8])
    # a near duplicate of contig 0 and a reverse-complemented contaminant window of contig 1
    c1 = contigs[f"{gid}_c0001"]
    x0 = rnd.randrange(len(c1) // 3)
    y0 = rnd.randint(x0 + 50, len(c1))
    win = rc(c1[x0:y0])
    cuid = blake3(win.encode()).hexdigest(length=16)
    contam = [{"contig_id": "g009_x0", "donor_contig_id": f"{gid}_c0001", "donor_start": x0, "donor_end": y0,
               "orientation": -1, "host_genome_id": "g009", "contig_uid": cuid}]
    dup = [{"kind": "near", "rep_contig_id": f"{gid}_c0000", "member_contig_uid": "ab" * 16,
            "member_genome_id": gid, "member_contig_id": f"{gid}_c0000_src_c"}]
    occ_seq = {r["contig_uid"]: contigs[r["contig_id"]] for r in rows if r["emitted"] == "true"}
    occ_seq[cuid] = win
    occ_seq["ab" * 16] = contigs[f"{gid}_c0000"]
    return {"rec": rec, "s5": s5, "s6": s6, "s7": s7, "p6": p6, "p7": p7, "p8": p8, "s8": s8, "rows": rows,
            "contigs": contigs, "contam": contam, "dup": dup, "occ_seq": occ_seq}


def _oracle(w, stage: int, s: int, e: int):
    """Where S<stage>-frame feature [s, e) is on the S8 frame, from provenance: (s8, e8), or None when a later
    edit touched it (a base replaced, or a new base inside)."""
    src = {5: None, 6: w["p6"], 7: w["p7"]}[stage]
    N = {5: len(w["s5"]), 6: len(w["s6"]), 7: len(w["s7"])}[stage]
    ids_ = [((p if src is None else src[p]), p) for p in ((s + i) % N for i in range(e - s))]
    if any(x is None for x, _ in ids_):
        return None
    # the S8 positions of those S5 positions, in order, must be one circular run of the same length
    where = {}
    for i, x in enumerate(w["p8"]):
        if x is not None:
            where[x] = i
    if stage == 5:
        want = [p for _, p in ids_]
    else:
        want = [x for x, _ in ids_]
    if any(x not in where for x in want):
        return None
    pos = [where[x] for x in want]
    N8 = len(w["p8"])
    if any((pos[i + 1] - pos[i]) % N8 != 1 for i in range(len(pos) - 1)):
        return None                      # a base new in a later stage lies inside the run
    return pos[0], pos[0] + len(pos)


def test_placer_lifts_and_places_features_as_the_provenance_oracle_does(tmp_path):
    rnd = random.Random(20261113)
    n = collections.Counter()
    for trial in range(12):
        d = tmp_path / f"w{trial}"
        w = _world(d, rnd)
        placer = emit.Placer(d, w["contam"], w["dup"])
        rec = w["rec"]
        for _ in range(250):
            stage = rnd.choice((5, 6, 7))
            src = {5: w["s5"], 6: w["s6"], 7: w["s7"]}[stage]
            N = len(src)
            s = rnd.randrange(N)
            e = s + rnd.randint(1, 300)
            st = rnd.choice((1, -1))
            if stage > 5 and any({6: w["p6"], 7: w["p7"]}[stage][p % N] is None for p in range(s, e)):
                continue           # a feature made of a later stage's own bases has no S5 provenance to follow
            want = _oracle(w, stage, s, e)
            got = placer.lift(rec, stage, s, e)
            if want is None:
                assert isinstance(got, str), (trial, stage, s, e, got)
                n["touched"] += 1
                continue
            N8 = len(w["s8"])
            assert not isinstance(got, str), (trial, stage, s, e, got)
            assert (got[0] % N8, got[1] - got[0]) == (want[0] % N8, want[1] - want[0])
            fb = (src + src)[s:e]
            fb = fb if st == 1 else rc(fb)
            status, pls = placer.place(rec, got[0], got[1], st)
            # brute force: which emitted contig bases hold the feature's S8 bases
            hold = collections.Counter()
            whole_on = set()
            for r in w["rows"]:
                a, b, o, L = int(r["replicon_start"]), int(r["replicon_end"]), int(r["orientation"]), int(r["length"])
                k = sum(1 for i in range(want[0], want[1]) if a <= i % N8 < b or a <= i % N8 + N8 < b)
                if k and r["emitted"] == "true":
                    hold[r["contig_uid"]] += k
                    hold["ab" * 16] += k if r["contig_id"].endswith("c0000") else 0
                if k == want[1] - want[0]:
                    whole_on.add(r["contig_id"])
            if not hold:
                assert status == "dropped" and not pls
                n["dropped"] += 1
                continue
            assert status == ("whole" if whole_on else "cut")
            got_k = collections.Counter()
            for pl in pls:
                seq = w["occ_seq"][pl.occ.uid]
                part = (seq + seq)[pl.start:pl.end] if pl.end > len(seq) else seq[pl.start:pl.end]
                part = part if pl.strand == 1 else rc(part)
                if pl.p5 or pl.p3:
                    assert part in fb
                else:
                    assert part == fb
                if pl.occ.origin != "contaminant":
                    got_k[pl.occ.uid] += pl.end - pl.start
            assert {k: v for k, v in got_k.items() if v} == {k: v for k, v in hold.items() if v}
            n[status] += 1
    print(f"\nplacements {dict(n)}")
    assert n["whole"] > 100 and n["cut"] > 5 and n["touched"] > 20 and n["dropped"] > 5
