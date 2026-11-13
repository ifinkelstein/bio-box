"""Ring R2: synthworld call errors, genetic codes 4 and 25, and assembly on the `micro` preset (synthworld_spec
§2.1 S8-S9, §3.3 truth_callerr and truth_edges, §3.2 truth_contig, §5.1, §6.1, §7.3; card F-02 session 7).

A module fixture builds the micro world's S1-S7 (ZOMBI2 0.54.0), timing each stage, records every file the
earlier stages wrote, and builds S8-S9 once. The oracles are written here: FASTA is parsed, reverse complements
and the genetic codes are spelled out, the S8 frame is rebuilt from the contigs by slicing, the lifts (S5 -> S6
by insertions.tsv, S6 -> S7 and S7 -> S8 by edits.tsv) are recomputed by this file's own formulas, ORFs are
scanned by brute force, clades are read from the species tree's Newick, and pyrodigal is the external gene
caller. What they take from the generator are model inputs (the preset, the S7 world), never its realized
statistics.
"""

from __future__ import annotations

import bisect
import collections
import csv
import hashlib
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest
from blake3 import blake3

from biobox.synthworld import assembly, backbone, fusions, intergenic, mobile, sequences

pytestmark = [pytest.mark.r2, pytest.mark.covers("F-02")]
needs_zombi2 = pytest.mark.skipif(
    importlib.util.find_spec("zombi2") is None,
    reason="zombi2 0.54.0 is not importable here: CI installs it by `uv sync`; on the pod the image has it")

#: the digests of sessions 1 to 6's stages for the micro default seed (measured 2026-11-13; planning/measurements/
#: 2026-11-13_synthworld_fusions.md): those stages must not change
DIGESTS = {"backbone": "025fdd0d27fdc65e678782216ecd21fb", "sequences": "34ebce1324413a261acd874d7ee4fe3e",
           "replicons": "55bda97155d75148e5cdca7b90a1d7be", "mobile": "0bd7f0db5691de58ea5e1a0df1e6e570",
           "fusions": "db09abf3fe79ae34a73e3ab0664a9b13"}
EARLIER = ("backbone", "zombi2", "sequences", "replicons", "mobile", "fusions")

CODE11_BY_AA = {
    "F": "TTT TTC", "L": "TTA TTG CTT CTC CTA CTG", "I": "ATT ATC ATA", "M": "ATG", "V": "GTT GTC GTA GTG",
    "S": "TCT TCC TCA TCG AGT AGC", "P": "CCT CCC CCA CCG", "T": "ACT ACC ACA ACG", "A": "GCT GCC GCA GCG",
    "Y": "TAT TAC", "*": "TAA TAG TGA", "H": "CAT CAC", "Q": "CAA CAG", "N": "AAT AAC", "K": "AAA AAG",
    "D": "GAT GAC", "E": "GAA GAG", "C": "TGT TGC", "W": "TGG", "R": "CGT CGC CGA CGG AGA AGG",
    "G": "GGT GGC GGA GGG",
}
CODE11 = {c: aa for aa, cs in CODE11_BY_AA.items() for c in cs.split()}
CODE = {11: CODE11, 4: dict(CODE11, TGA="W"), 25: dict(CODE11, TGA="G")}   # NCBI tables 4 and 25 reassign TGA
STARTS = {"ATG", "GTG", "TTG"}
_RC = {"A": "T", "C": "G", "G": "C", "T": "A"}


def rc(s: str) -> str:
    return "".join(_RC[c] for c in reversed(s))


def codons_of(cds: str) -> list[str]:
    return [cds[i:i + 3] for i in range(0, len(cds) - len(cds) % 3, 3)]


def translate(cds: str, code: int = 11) -> str:
    aa = [CODE[code][c] for c in codons_of(cds)]
    assert len(cds) % 3 == 0 and cds[:3] in STARTS and aa[-1] == "*" and "*" not in aa[:-1], cds[:12]
    return "M" + "".join(aa[1:-1])


def uid(s: str) -> str:
    return blake3(s.encode()).hexdigest(length=16)


def _tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _fasta(path: Path) -> dict[str, str]:
    out, name, chunks = {}, None, []
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(chunks)
            name, chunks = line[1:].split()[0], []
        else:
            assert len(line) <= 80
            chunks.append(line)
    if name is not None:
        out[name] = "".join(chunks)
    return out


def _tree_sha(root: Path) -> dict[str, str]:
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file()}


def _cut(seq: str, start: int, end: int, strand: int = 1) -> str:
    n = len(seq)
    assert 0 <= start <= end <= start + n, (start, end, n)
    s = (seq + seq)[start:end] if end > n else seq[start:end]
    return s if strand == 1 else rc(s)


def _stage_digests(out: Path) -> dict[str, str]:
    return {d: json.loads((out / d / f"{d}.json").read_text())["digest"] for d in DIGESTS}


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    if importlib.util.find_spec("zombi2") is None:
        pytest.skip("zombi2 0.54.0 is not importable here")
    out = tmp_path_factory.mktemp("micro_assembly") / "world"
    params = backbone.load_preset("micro")
    times = {}
    for name, fn in (("S1-S2 backbone", lambda: backbone.build(params, params["default_seed"], out)),
                     ("S3-S4 sequences", lambda: sequences.build(out)),
                     ("S5 replicons", lambda: intergenic.build(out)),
                     ("S6 mobile", lambda: mobile.build(out)),
                     ("S7 fusions", lambda: fusions.build(out))):
        t0 = time.perf_counter()
        fn()
        times[name] = time.perf_counter() - t0
    before = {d: _tree_sha(out / d) for d in EARLIER}
    digests = _stage_digests(out)
    t0 = time.perf_counter()
    aw = assembly.build(out)
    times["S8-S9 assembly"] = time.perf_counter() - t0
    return aw, out, times, before, digests


@pytest.fixture(scope="module")
def t(world):
    """The tables, and the S8 frame rebuilt from the contigs (emitted and dropped) by slicing."""
    _, out, _, _, _ = world
    ad = out / "assembly"
    contigs = _tsv(ad / "contigs.tsv")
    reps = {r["record"]: r for r in _tsv(ad / "replicons.tsv")}
    cseq = {}
    for d in ("contigs", "dropped"):
        for f in sorted((ad / d).glob("*.fna")):
            cseq.update(_fasta(f))
    frame_buf = {rec: bytearray(b"-" * int(r["length"])) for rec, r in reps.items()}
    cover = {rec: np.zeros(int(r["length"]), dtype=np.int16) for rec, r in reps.items()}
    for c in contigs:
        rec, a, b = c["record"], int(c["replicon_start"]), int(c["replicon_end"])
        s = cseq[c["contig_id"]]
        assert len(s) == int(c["length"]) == b - a
        fwd = (s if c["orientation"] == "1" else rc(s)).encode()
        N = len(frame_buf[rec])
        if b <= N:
            frame_buf[rec][a:b] = fwd
            cover[rec][a:b] += 1
        else:
            frame_buf[rec][a:N] = fwd[:N - a]
            frame_buf[rec][0:b - N] = fwd[N - a:]
            cover[rec][a:N] += 1
            cover[rec][0:b - N] += 1
    frames = {rec: bytes(buf).decode() for rec, buf in frame_buf.items()}
    return {"contigs": contigs, "cseq": cseq, "reps": reps, "frames": frames, "cover": cover,
            "genes": _tsv(ad / "genes.tsv"), "edits": _tsv(ad / "edits.tsv"), "callerr": _tsv(ad / "callerr.tsv"),
            "pieces": _tsv(ad / "gene_contigs.tsv"), "genomes": _tsv(ad / "genomes.tsv"),
            "repeats": _tsv(ad / "repeats.tsv"), "g7": _tsv(out / "fusions" / "genes.tsv"),
            "s7": {k: v for f in sorted((out / "fusions" / "fasta").glob("*.fna")) for k, v in _fasta(f).items()}}


def _species(out: Path):
    """child -> parent of the complete species tree, from its Newick (parsed here)."""
    text = (out / "zombi2" / "species_complete.nwk").read_text(encoding="utf-8").strip()
    tokens, buf = [], ""
    for ch in text:
        if ch in "(),;":
            if buf:
                tokens.append(buf)
                buf = ""
            tokens.append(ch)
        else:
            buf += ch
    parent, stack, nodes, last, label_next = {}, [], [], None, False
    for tok in tokens:
        if tok == "(":
            nodes.append(None)
            stack.append(len(nodes) - 1)
            label_next = False
        elif tok == ")":
            last = stack.pop()
            label_next = True
        elif tok in ",;":
            label_next = False
        else:
            name = tok.split(":")[0]
            if label_next:
                nodes[last] = name
                k = last
            else:
                nodes.append(name)
                k = len(nodes) - 1
            parent[k] = stack[-1] if stack else None
            label_next = False
    return {nodes[k]: (nodes[p] if p is not None else None) for k, p in parent.items()}


def _below(par: dict, node: str, tips) -> set[str]:
    out = set()
    for tip in tips:
        x = tip
        while x is not None:
            if x == node:
                out.add(tip)
                break
            x = par[x]
    return out


# --- this file's lifts and placement -----------------------------------------------------------------------------

class _EditLift:
    """Replacements [a, b) -> l bases (a = b an insertion) on a record of N bases, as position masks: a feature is
    touched when a replaced base lies in it or an insertion point lies strictly inside it, and otherwise moves by
    the summed length change of the edits ending at or before its start."""

    def __init__(self, N: int, edits):
        self.N = N
        self.replaced = np.zeros(N + 1, dtype=bool)
        self.point = np.zeros(N + 1, dtype=bool)
        delta = np.zeros(N + 2, dtype=np.int64)
        for a, b, l in edits:
            if a == b:
                self.point[a] = True
            else:
                self.replaced[a:b] = True
            delta[b] += l - (b - a)
        self.before = np.cumsum(delta)          # before[x]: change of the edits with b <= x
        self.total = int(delta.sum())

    def _touched(self, s, e):
        if self.replaced[s:e].any() or self.point[s + 1:e].any():
            return True
        return False

    def lift(self, s: int, e: int):
        N = self.N
        if e <= N:
            if self._touched(s, e):
                return None
            d = int(self.before[s])
            return s + d, e + d
        if self._touched(s, N) or self.point[N] or self.point[0] or self._touched(0, e - N):
            return None
        return s + int(self.before[s]), N + self.total + (e - N)


def _ins_lift(points, N, s, e):
    """S5 -> S6: blocks at points a (N for the origin) of lengths l."""
    s2 = s + sum(l for a, l in points if a <= s)
    if e <= N:
        return s2, e + sum(l for a, l in points if a < e)
    return s2, N + sum(l for _, l in points) + (e - N)


class _Contigs:
    """Where a frame interval lies on the contigs of one record (written here: interval arithmetic on the frame
    spans of contigs.tsv)."""

    def __init__(self, rows, N, circular):
        self.rows = sorted(rows, key=lambda r: int(r["replicon_start"]))
        self.starts = [int(r["replicon_start"]) for r in self.rows]
        self.N, self.circular = N, circular

    def parts(self, s, e, strand):
        """[(contig row, cs, ce, cstrand, clipped 5', clipped 3')]."""
        N, out = self.N, []
        if len(self.rows) == 1 and self.rows[0]["is_circular"] == "true":
            r = self.rows[0]
            a, L = int(r["replicon_start"]), int(r["length"])
            x = (s - a) % N
            if r["orientation"] == "1":
                return [(r, x, x + (e - s), strand, False, False)]
            cs = (L - x - (e - s)) % L
            return [(r, cs, cs + (e - s), -strand, False, False)]
        for shift in ((-N, 0, N) if self.circular else (0,)):
            fs, fe = s + shift, e + shift
            for r in self.rows:
                a, b = int(r["replicon_start"]), int(r["replicon_end"])
                if a < fe and b > fs:
                    x, y = max(a, fs) - a, min(b, fe) - a
                    cl, cr = fs < a, fe > b
                    L = b - a
                    if r["orientation"] == "1":
                        cs, ce, st = x, y, strand
                    else:
                        cs, ce, st = L - y, L - x, -strand
                    p5, p3 = (cl, cr) if strand == 1 else (cr, cl)
                    out.append((r, cs, ce, st, p5, p3))
        return out


def _orf_pieces(cds: str, min_codons: int):
    """Brute force: per frame the stop-to-stop stretches; the 5' piece is frame 0's first, the 3' piece the one
    ending at the true stop, from its first start codon."""
    L = len(cds)
    stretches = []
    for f in range(3):
        prev = f
        for i in range(f, L - 2, 3):
            if CODE11[cds[i:i + 3]] == "*":
                stretches.append((prev, i + 3))
                prev = i + 3
    out = []
    five = next(((a, b) for a, b in stretches if a == 0), None)
    if five and (five[1] - 3) // 3 >= min_codons:
        out.append(five)
    three = next(((a, b) for a, b in stretches if b == L), None)
    if three:
        st = next((i for i in range(three[0], L - 3, 3) if cds[i:i + 3] in STARTS), None)
        if st is not None and (L - 3 - st) // 3 >= min_codons and (st, L) not in out:
            out.append((st, L))
    return sorted(out)


# --- the earlier stages do not change (sessions 1 to 6) ------------------------------------------------------------

@needs_zombi2
def test_earlier_stages_keep_every_file_and_their_digests(world):
    aw, out, times, before, digests = world
    for d in EARLIER:
        assert _tree_sha(out / d) == before[d], d
    assert _stage_digests(out) == digests == DIGESTS
    assert {d: aw.meta[f"{d}_digest"] for d in DIGESTS} == DIGESTS
    print(f"\nearlier stage digests {digests}")


# --- determinism (§6.1, §6.2; MR-16) ---------------------------------------------------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "MR-16")
def test_same_seed_gives_byte_identical_assembly_stage_in_a_subprocess(world, tmp_path, record_property):
    _, out, _, _, _ = world
    dst = tmp_path / "copy"
    dst.mkdir()
    for d in EARLIER:
        shutil.copytree(out / d, dst / d)
    env = dict(os.environ, PYTHONHASHSEED="4242", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    code = ("import os, sys, time; from biobox.synthworld import assembly; t0 = time.perf_counter(); "
            "assembly.build(sys.argv[1]); dt = time.perf_counter() - t0; p = '/proc/self/status'; "
            "hwm = next((l.split()[1] for l in open(p) if l.startswith('VmHWM')), -1) if os.path.exists(p) else -1; "
            "print(hwm, round(dt, 2))")
    proc = subprocess.run([sys.executable, "-c", code, str(dst)], check=True, env=env, cwd=tmp_path,
                          capture_output=True, text=True)
    assert _tree_sha(dst / "assembly") == _tree_sha(out / "assembly")
    hwm_kb, secs = proc.stdout.split()[-2:]
    print(f"\nassembly stage in a child: {secs} s wall, VmHWM {int(hwm_kb) / 1024:.0f} MB")
    record_property("assembly_child_seconds", float(secs))
    if int(hwm_kb) > 0:
        record_property("assembly_child_maxrss_mb", round(int(hwm_kb) / 1024, 1))
        assert int(hwm_kb) / 1024 < 1024      # one genome at a time (AGENTS §4.1)


# --- S9: the contigs tile the S8 frame; classes, breaks, edges (truth_edges, truth_contig) ----------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "MR-04", "MR-05")
def test_contigs_tile_the_s8_frame_and_their_edges_are_true(world, t):
    _, out, _, _, _ = world
    cfg = backbone.load_preset("micro")["assembly"]
    for rec, cov in t["cover"].items():
        assert int(cov.min()) == 1 and int(cov.max()) == 1, rec          # every base in exactly one contig
        assert uid(t["frames"][rec]) == t["reps"][rec]["contig_uid"]
        assert len(t["frames"][rec]) == int(t["reps"][rec]["length"])
    gen = {g["genome_id"]: g for g in t["genomes"]}
    n = len(gen)
    want = [int(round(cfg["class_weights"][0] * n)), int(round(cfg["class_weights"][1] * n))]
    classes = collections.Counter(g["assembly_class"] for g in gen.values())
    assert [classes["complete"], classes["draft"]] == want and classes["mag"] == n - sum(want)
    emitted_files = {}
    for f in sorted((out / "assembly" / "contigs").glob("*.fna")):
        emitted_files.update(_fasta(f))
    by_genome = collections.defaultdict(list)
    cuts = collections.defaultdict(dict)        # record -> frame position -> (edge, family)
    for c in t["contigs"]:
        by_genome[c["genome_id"]].append(c)
        assert uid(t["cseq"][c["contig_id"]]) == c["contig_uid"]
        assert (c["contig_id"] in emitted_files) == (c["emitted"] == "true")
        L = int(c["length"])
        klass = gen[c["genome_id"]]["assembly_class"]
        topo = t["reps"][c["record"]]["topology"]
        assert (c["is_circular"] == "true") == (klass == "complete" and topo == "circular")
        if c["drop_reason"] == "short":
            assert L < cfg["min_contig"] and c["emitted"] == "false"
        elif c["emitted"] == "true":
            assert L >= cfg["min_contig"] and c["drop_reason"] == ""
        else:
            assert c["drop_reason"] == "mag_dropout" and klass == "mag"
        # edges in frame direction: the edge at replicon_start and the one at replicon_end
        lo_edge, hi_edge = (c["left_edge"], c["right_edge"]) if c["orientation"] == "1" else \
            (c["right_edge"], c["left_edge"])
        lo_fam, hi_fam = (c["repeat_family_left"], c["repeat_family_right"]) if c["orientation"] == "1" else \
            (c["repeat_family_right"], c["repeat_family_left"])
        if c["is_circular"] == "true":
            assert lo_edge == hi_edge == "none" and L == int(t["reps"][c["record"]]["length"])
            assert int(c["rotation_offset"]) == int(c["replicon_start"])
            continue
        N = int(t["reps"][c["record"]]["length"])
        for pos, edge, fam in ((int(c["replicon_start"]), lo_edge, lo_fam), (int(c["replicon_end"]) % N, hi_edge, hi_fam)):
            if edge == "replicon_end":
                assert topo == "linear" and pos == 0
            else:
                assert edge in ("break_repeat", "break_random") and (fam != "") == (edge == "break_repeat")
                assert cuts[c["record"]].setdefault(pos, (edge, fam)) == (edge, fam)
    # each cut is a contig's start and another's end; a repeat break is at the midpoint of a repeat copy
    starts = collections.Counter((c["record"], int(c["replicon_start"])) for c in t["contigs"]
                                 if c["is_circular"] == "false")
    rep_mid = {(r["genome_id"], int(r["replicon"]), (int(r["start"]) + int(r["end"])) // 2 % len(t["frames"][
        f"{r['genome_id']}_r{r['replicon']}"])): r["family_id"] for r in t["repeats"]}
    n_rep = 0
    for rec, cs in cuts.items():
        gid, rep = rec.rsplit("_r", 1)
        for pos, (edge, fam) in cs.items():
            assert starts[(rec, pos)] == 1
            if edge == "break_repeat":
                assert rep_mid[(gid, int(rep), pos)] == fam
                n_rep += 1
    # N50, completeness and the MAG dropout rule
    ratios = []
    for gid, g in gen.items():
        cs = by_genome[gid]
        em = sorted((int(c["length"]) for c in cs if c["emitted"] == "true"), reverse=True)
        G = int(g["length"])
        assert G == sum(len(t["frames"][r]) for r in {c["record"] for c in cs})
        acc, n50 = 0, 0
        for x in em:
            acc += x
            if 2 * acc >= sum(em):
                n50 = x
                break
        assert int(g["n50"]) == n50 and int(g["emitted_length"]) == sum(em)
        if g["assembly_class"] == "mag":
            c = float(g["target_completeness"])
            lo, hi = cfg["mag_completeness"]
            assert lo <= c <= hi and sum(em) / G >= c and (sum(em) - min(em)) / G < c
        if g["assembly_class"] != "complete":
            ratios.append(n50 / int(g["target_n50"]))
            assert int(g["target_n50"]) in (cfg["draft_n50"] + [cfg["mag_n50"]])
    gm = math.exp(sum(math.log(x) for x in ratios) / len(ratios))
    assert 0.5 < gm < 2.0
    # every complete circular contig has a gene across its origin (MR-04)
    circ = {c["contig_id"] for c in t["contigs"] if c["is_circular"] == "true"}
    spanning = {p["contig_id"] for p in t["pieces"] if p["contig_id"] in circ and p["origin_spanning"] == "true"}
    assert spanning == circ and len(circ) >= 5
    # orientation is a fair coin
    k = sum(c["orientation"] == "-1" for c in t["contigs"])
    m = len(t["contigs"])
    assert abs(k - m / 2) < 4 * math.sqrt(m / 4)
    print(f"\nclasses {dict(classes)}; contigs {m} ({sum(c['emitted'] == 'true' for c in t['contigs'])} emitted, "
          f"{sum(c['drop_reason'] == 'short' for c in t['contigs'])} short, "
          f"{sum(c['drop_reason'] == 'mag_dropout' for c in t['contigs'])} MAG dropout); repeat breaks {n_rep}; "
          f"N50/target geometric mean {gm:.3f}, range {min(ratios):.3f}-{max(ratios):.3f}; "
          f"circular contigs {len(circ)}; reversed {k}")


# --- S8: the frame is the S7 replicons with the listed edits ------------------------------------------------------

@needs_zombi2
def test_s8_frame_is_the_s7_replicons_with_the_listed_edits(t):
    by_rec = collections.defaultdict(list)
    for r in t["edits"]:
        by_rec[r["record"]].append(r)
    kinds = collections.Counter(r["kind"] for r in t["edits"])
    for rec, s7 in t["s7"].items():
        ed = sorted(by_rec.get(rec, []), key=lambda r: (int(r["start"]), int(r["end"])))
        parts, prev, shift = [], 0, 0
        for r in ed:
            a, b, l = int(r["start"]), int(r["end"]), int(r["length"])
            assert a >= prev
            parts.append(s7[prev:a])
            new = t["frames"][rec][a + shift:a + shift + l]
            assert (int(r["new_start"]), int(r["new_end"])) == (a + shift, a + shift + l)
            assert new != s7[a:b] or a == b
            parts.append(new)
            shift += l - (b - a)
            prev = b
        parts.append(s7[prev:])
        assert "".join(parts) == t["frames"][rec], rec
        assert int(t["reps"][rec]["s7_length"]) == len(s7)
    print(f"\nedits by kind {dict(kinds)}")


# --- genetic codes 4 and 25 --------------------------------------------------------------------------------------

@needs_zombi2
def test_one_clade_uses_code_4_and_one_genome_code_25(world, t):
    _, out, _, _, _ = world
    cfg = backbone.load_preset("micro")["callerr"]
    code = {g["genome_id"]: int(g["genetic_code"]) for g in t["genomes"]}
    node_gid = {r["species_node"]: r["genome_id"] for r in _tsv(out / "backbone" / "genomes.tsv")}
    par = _species(out)
    four = {g for g, c in code.items() if c == 4}
    assert collections.Counter(code.values())[25] == 1
    lo, hi = cfg["code4_clade_size"]
    clades = [x for x in par if {node_gid[y] for y in _below(par, x, node_gid)} == four]
    assert clades and lo <= len(four) <= hi
    # no gene of a code-4 or code-25 genome stops with TGA; internal TGA codons carry Trp (4) or Gly (25)
    n_w = collections.Counter()
    alt_rows = {r["true_gene_id"]: r for r in t["callerr"] if r["kind"] == "alt_code"}
    with_tga = set()
    for g in t["genes"]:
        c = code[g["genome_id"]]
        if c == 11 or g["feature_kind"] != "CDS":
            continue
        cds = _cut(t["frames"][f"{g['genome_id']}_r{g['replicon']}"], int(g["start"]), int(g["end"]), int(g["strand"]))
        cods = codons_of(cds)
        assert cods[-1] in ("TAA", "TAG")
        if "TGA" in cods[1:-1]:
            with_tga.add(g["copy_label"])
            first = cods.index("TGA", 1)
            r = alt_rows[g["copy_label"]]
            pos = int(r["error_pos"])
            N = len(t["frames"][r["record"]])
            want = (int(g["start"]) + 3 * first) % N if g["strand"] == "1" else (int(g["end"]) - 3 * first - 3) % N
            assert pos == want
        if c == 4:
            n_w["TGA"] += cods[1:-1].count("TGA")
            n_w["TGG"] += cods[1:-1].count("TGG")
    assert with_tga == set(alt_rows)
    share = n_w["TGA"] / (n_w["TGA"] + n_w["TGG"])
    assert share > 0.8                    # code4_tgg_to_tga 0.9 less the codons shared with another gene
    print(f"\ncode 4: {sorted(four)}; code 25: {[g for g, c in code.items() if c == 25]}; alt_code genes "
          f"{len(alt_rows)}; Trp as TGA in code 4 {share:.3f}")


# --- every gene reads as its model ------------------------------------------------------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "B-SIM")
def test_every_gene_translates_under_its_genome_code_to_its_protein_uid(t):
    code = {g["genome_id"]: int(g["genetic_code"]) for g in t["genomes"]}
    err = {r["callerr_id"]: r for r in t["callerr"]}
    g7 = {g["copy_label"]: g for g in t["g7"]}
    n = collections.Counter()
    labels = set()
    for g in t["genes"]:
        labels.add(g["copy_label"])
        seq = t["frames"][f"{g['genome_id']}_r{g['replicon']}"]
        s, e, st = int(g["start"]), int(g["end"]), int(g["strand"])
        assert 0 <= s < len(seq) and s < e < s + len(seq)
        assert (g["origin_spanning"] == "true") == (e > len(seq))
        assert int(g["genetic_code"]) == code[g["genome_id"]]
        cds = _cut(seq, s, e, st)
        kind = err[g["callerr_id"]]["kind"] if g["callerr_id"] else ""
        if g["feature_kind"] == "pseudogene":
            assert kind == "internal_stop" and g["protein_uid"] == ""
            assert "*" in [CODE11[c] for c in codons_of(cds)[:-1]]
        elif kind.startswith("frameshift"):
            old = _cut(t["s7"][f"{g['genome_id']}_r{g['replicon']}"], int(g7[g["copy_label"]]["start"]),
                       int(g7[g["copy_label"]]["end"]), st)
            m = min(len(old), len(cds))
            p = next((i for i in range(m) if old[i] != cds[i]), m)            # the first base that differs
            if kind == "frameshift_del":
                assert len(cds) == len(old) - 1 and old[:p] + old[p + 1:] == cds
            else:
                assert len(cds) == len(old) + 1 and cds[:p] + cds[p + 1:] == old
            assert uid(translate(old)) == g["protein_uid"]
        else:
            assert uid(translate(cds, code[g["genome_id"]])) == g["protein_uid"], g["copy_label"]
        n[g["feature_kind"] + ("/" + kind if kind else "")] += 1
    smorfs = {g["copy_label"] for g in t["genes"] if g["feature_kind"] == "smORF"}
    assert labels == set(g7) | smorfs
    assert sum(g["origin_spanning"] == "true" for g in t["genes"]) >= 5       # §5.1
    print(f"\ngenes by kind {dict(n)}")


# --- call errors are planted as designed (truth_callerr) -------------------------------------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "F-08")
def test_call_errors_are_planted_as_designed_and_expected_calls_are_the_code11_orfs(t):
    cfg = backbone.load_preset("micro")["callerr"]
    code = {g["genome_id"]: int(g["genetic_code"]) for g in t["genomes"]}
    genes = {g["copy_label"]: g for g in t["genes"]}
    g7 = {g["copy_label"]: g for g in t["g7"]}
    by_rec = collections.defaultdict(list)
    for g in t["genes"]:
        by_rec[f"{g['genome_id']}_r{g['replicon']}"].append(g)
    status = collections.defaultdict(list)
    for p in t["pieces"]:
        status[p["copy_label"]].append(p)
    kinds = collections.Counter()
    for r in t["callerr"]:
        kinds[r["kind"]] += 1
        g = genes[r["true_gene_id"]]
        seq = t["frames"][r["record"]]
        N = len(seq)
        s, e, st = int(r["true_start"]), int(r["true_end"]), int(r["strand"])
        assert (s, e, st) == (int(g["start"]), int(g["end"]), int(g["strand"]))
        cds = _cut(seq, s, e, st)
        k = r["kind"]
        if k != "alt_code":
            assert code[r["genome_id"]] == 11 and g["callerr_id"] == r["callerr_id"]
        pos = int(r["error_pos"])
        q = (pos - s) % N if st == 1 else (e - 1 - pos) % N          # the error's base index on the gene strand
        n_cod = len(cds) // 3
        if k in ("frameshift_ins", "frameshift_del", "internal_stop"):
            assert n_cod >= cfg["min_gene_codons"] - 1 and n_cod // 3 - 1 <= q // 3 <= 2 * n_cod // 3 + 1
            if k == "internal_stop":
                assert CODE11[cds[3 * (q // 3):3 * (q // 3) + 3]] == "*" and r["true_protein_uid"] == ""
            want = _orf_pieces(cds, cfg["min_call_codons"])
        elif k == "alt_start_downstream":
            assert cds[:3] == "TTG" and q // 3 == 0
            old = _cut(t["s7"][r["record"]], int(g7[g["copy_label"]]["start"]), int(g7[g["copy_label"]]["end"]), st)
            assert old[:3] in ("ATG", "GTG") and translate(old) == translate(cds)
            a = int(r["alt_start"])
            j = (a - s) % N // 3 if st == 1 else (e - 3 - a) % N // 3
            lo, hi = cfg["alt_downstream_codons"]
            assert lo <= j <= hi and cds[3 * j:3 * j + 3] == "ATG"
            want = [(0, len(cds))]
        elif k == "alt_start_upstream":
            a = int(r["alt_start"])
            kk = (s - a) // 3 if st == 1 else (a + 3 - e) // 3
            lo, hi = cfg["alt_upstream_codons"]
            assert lo <= kk <= hi and a == pos
            up = _cut(seq, a, s, 1) if st == 1 else _cut(seq, e, a + 3, -1)
            assert len(up) == 3 * kk and up[:3] == "ATG"
            assert all(CODE11[c] != "*" for c in codons_of(up))
            for h in by_rec[r["record"]]:                                # no gene in the window
                hs, he = int(h["start"]), int(h["end"])
                w0, w1 = (a, s) if st == 1 else (e, a + 3)
                assert he <= w0 or hs >= w1 or h is g
            old = _cut(t["s7"][r["record"]], int(g7[g["copy_label"]]["start"]), int(g7[g["copy_label"]]["end"]), st)
            assert old == cds
            want = [(0, len(cds))]
        elif k == "smorf":
            lo, hi = cfg["smorf_aa"]
            assert g["family_id"] == "smORF" and lo <= n_cod - 1 <= hi
            assert translate(cds)
            pad = cfg["smorf_pad"]
            up = _cut(seq, s - pad - 6 - 10, s, 1) if st == 1 else _cut(seq, e, e + pad + 6 + 10, -1)
            assert "AGGAGG" in up
            for h in by_rec[r["record"]]:
                assert int(h["end"]) <= s or int(h["start"]) >= e or h is g
            want = [(0, len(cds))]
        else:
            assert code[r["genome_id"]] in (4, 25)
            want = _orf_pieces(cds, cfg["min_call_codons"])
        got = []
        for item in filter(None, r["expected_calls"].split(";")):
            span, strand = item.split(":")
            x, y = (int(v) for v in span.split("-"))
            assert int(strand) == st
            # back to offsets on the gene strand
            a = (x - s) % N if st == 1 else (e - y) % N
            got.append((a, a + (y - x)))
        assert sorted(got) == want, (r["callerr_id"], got, want)
        ps = status[r["true_gene_id"]]
        em = [p for p in ps if p["emitted"] == "true"]
        want_status = "dropped" if not em else "whole" if len(ps) == 1 and ps[0]["partial_5"] == ps[0]["partial_3"] \
            == "false" else "partial"
        assert r["status"] == want_status
    for kind, n in cfg["counts"].items():
        assert kinds[kind] == n, kind
    print(f"\ncall errors by kind {dict(kinds)}; status {dict(collections.Counter(r['status'] for r in t['callerr']))}")


# --- the external gene caller sees the errors -----------------------------------------------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "F-08")
def test_pyrodigal_sees_the_planted_errors_and_the_alternative_code(t, record_property):
    pyrodigal = pytest.importorskip("pyrodigal")
    code = {g["genome_id"]: int(g["genetic_code"]) for g in t["genomes"]}
    errs = collections.defaultdict(list)
    for r in t["callerr"]:
        errs[r["genome_id"]].append(r)
    g11 = max((g for g in code if code[g] == 11), key=lambda g: (len(errs[g]), g))
    g4 = sorted(g for g in code if code[g] == 4)[0]

    def calls(gid, table):
        recs = sorted(rec for rec in t["frames"] if rec.startswith(gid + "_r"))
        gf = pyrodigal.GeneFinder(meta=False)
        gf.train(*[t["frames"][rec].encode() for rec in recs], translation_table=table)
        out = {}
        for rec in recs:
            out[rec] = {(x.begin - 1, x.end, x.strand) for x in gf.find_genes(t["frames"][rec].encode())}
        return out

    c11 = calls(g11, 11)
    res = collections.Counter()
    for r in errs[g11]:
        s, e, st = int(r["true_start"]), int(r["true_end"]), int(r["strand"])
        got = c11[r["record"]]
        stop_match = any(x[2] == st and (x[1] == e if st == 1 else x[0] == s) for x in got)
        exact = (s, e, st) in got
        if r["kind"] in ("frameshift_ins", "frameshift_del", "internal_stop"):
            assert not exact, r["callerr_id"]           # a broken gene is never called whole
        res[(r["kind"], "exact")] += exact
        res[(r["kind"], "stop")] += stop_match
        res[(r["kind"], "n")] += 1
    # the code-4 genome: stop-anchored recall of its genes with an internal TGA, code 11 against code 4
    alt = [r for r in errs[g4] if r["kind"] == "alt_code"]
    rec4 = {}
    for table in (11, 4):
        cc = calls(g4, table)
        hit = sum(any(x[2] == int(r["strand"]) and (x[1] == int(r["true_end"]) if r["strand"] == "1"
                                                     else x[0] == int(r["true_start"])) for x in cc[r["record"]])
                  for r in alt)
        rec4[table] = hit / len(alt)
    assert len(alt) > 100 and rec4[4] > rec4[11] + 0.5
    summary = {f"{k}_{m}": v for (k, m), v in sorted(res.items())}
    print(f"\npyrodigal on {g11} (code 11): {summary}\npyrodigal on {g4}: stop-anchored recall of {len(alt)} "
          f"genes with an internal TGA, table 11 {rec4[11]:.3f}, table 4 {rec4[4]:.3f}")
    record_property("pyrodigal_code4_recall_t11", round(rec4[11], 3))
    record_property("pyrodigal_code4_recall_t4", round(rec4[4], 3))


# --- every gene piece and every earlier feature maps onto the contigs ---------------------------------------------------

def _contig_index(t):
    by_rec = collections.defaultdict(list)
    for c in t["contigs"]:
        by_rec[c["record"]].append(c)
    return {rec: _Contigs(rows, len(t["frames"][rec]), t["reps"][rec]["topology"] == "circular")
            for rec, rows in by_rec.items()}


@needs_zombi2
@pytest.mark.covers("F-02", "MR-05")
def test_every_gene_piece_reads_the_frame_bases_on_its_contig(t):
    idx = _contig_index(t)
    rows = collections.defaultdict(list)
    for p in t["pieces"]:
        rows[p["copy_label"]].append((p["contig_id"], int(p["start"]), int(p["end"]), int(p["strand"]),
                                      p["partial_5"] == "true", p["partial_3"] == "true"))
    n = collections.Counter()
    for g in t["genes"]:
        rec = f"{g['genome_id']}_r{g['replicon']}"
        s, e, st = int(g["start"]), int(g["end"]), int(g["strand"])
        parts = idx[rec].parts(s, e, st)
        want = sorted((r["contig_id"], cs, ce, cst, p5, p3) for r, cs, ce, cst, p5, p3 in parts)
        assert sorted(rows[g["copy_label"]]) == want, g["copy_label"]
        whole = _cut(t["frames"][rec], s, e, st)
        got = ""
        for r, cs, ce, cst, p5, p3 in sorted(parts, key=lambda x: (x[4], not x[5])):
            piece = _cut(t["cseq"][r["contig_id"]], cs, ce, cst)
            assert piece in whole
            got += piece
        assert len(got) == len(whole)
        n["whole" if len(parts) == 1 and not (parts[0][4] or parts[0][5]) else "cut"] += 1
        n["on_dropped"] += any(r["emitted"] == "false" for r, *_ in parts)
    print(f"\ngenes on contigs {dict(n)}")
    assert n["cut"] > 0 and n["on_dropped"] > 0


@needs_zombi2
def test_every_earlier_feature_maps_onto_the_contigs_or_is_accounted_for(world, t):
    """S5 elements, motifs and IGRs, S6 insertions, prophages, arrays and spacer slots, and S7 genes, lifted
    through S6's insertions and S7's and S8's edits, read their own stage's bases on the contigs (whole, or in
    parts across a break), or an edit of their stage's successors touched them: S7's (the fused pairs, already
    tested by the S7 tests) or S8's, which touches a gene only by its own call error or by recoding, an IGR only
    by an upstream alternative start or a smORF, and never an element, a motif or a spacer."""
    _, out, _, _, _ = world
    idx = _contig_index(t)
    code = {g["genome_id"]: int(g["genetic_code"]) for g in t["genomes"]}
    s5 = {k: v for f in sorted((out / "replicons" / "fasta").glob("*.fna")) for k, v in _fasta(f).items()}
    s6 = {k: v for f in sorted((out / "mobile" / "fasta").glob("*.fna")) for k, v in _fasta(f).items()}
    ins = collections.defaultdict(list)
    for r in _tsv(out / "mobile" / "insertions.tsv"):
        ins[r["record"]].append((int(r["at"]), int(r["length"])))
    e7 = collections.defaultdict(list)
    for r in _tsv(out / "fusions" / "edits.tsv"):
        e7[r["record"]].append((int(r["start"]), int(r["end"]), int(r["length"])))
    e8 = collections.defaultdict(list)
    e8_kind = collections.defaultdict(list)
    for r in t["edits"]:
        e8[r["record"]].append((int(r["start"]), int(r["end"]), int(r["length"])))
        e8_kind[r["record"]].append(r)
    errors_by_gene = {r["true_gene_id"]: r["kind"] for r in t["callerr"] if r["kind"] != "alt_code"}
    features = collections.defaultdict(list)      # record -> (stage, kind, s, e, strand, name)
    for g in t["g7"]:
        features[f"{g['genome_id']}_r{g['replicon']}"].append((7, "gene", int(g["start"]), int(g["end"]),
                                                               int(g["strand"]), g["copy_label"]))
    for name, cols in (("insertions.tsv", ("start", "end")), ("prophages.tsv", ("start", "end")),
                       ("arrays.tsv", ("locus_start", "locus_end")),
                       ("spacer_targets.tsv", ("spacer_start", "spacer_end"))):
        for r in _tsv(out / "mobile" / name):
            if name == "spacer_targets.tsv" and r["target_kind"] not in ("phage", "none", "prophage"):
                continue
            features[r["record"]].append((6, name.split(".")[0], int(r[cols[0]]), int(r[cols[1]]), 1, ""))
    for name in ("elements.tsv", "rna.tsv", "junctions.tsv"):
        for r in _tsv(out / "replicons" / name):
            s, e = int(r["start"]), int(r["end"])
            if e > s:
                features[f"{r['genome_id']}_r{r['replicon']}"].append((5, name.split(".")[0], s, e, 1, ""))
    n = collections.Counter()
    for rec in sorted(features):
        N5, N6, N7 = len(s5[rec]), len(s6[rec]), len(t["s7"][rec])
        l7 = _EditLift(N6, e7[rec])
        l8 = _EditLift(N7, e8[rec])
        own8 = None
        for stage, kind, s, e, st, name in features[rec]:
            src = {5: s5, 6: s6, 7: t["s7"]}[stage][rec]
            bases = _cut(src, s, e)
            x = (s, e)
            if stage == 5:
                x = _ins_lift(ins[rec], N5, *x)
            if stage <= 6:
                x = l7.lift(*x)
                if x is None:
                    n[f"{kind}_s7_event"] += 1
                    continue
            y = l8.lift(*x)
            if y is None:
                if own8 is None:
                    own8 = [(int(r["start"]), int(r["end"]), r["kind"], r["callerr_id"]) for r in e8_kind[rec]]
                hit = {k for a, b, k, _ in own8 if (a < b and a < x[1] and b > x[0]) or (a == b and x[0] < a < x[1])
                       or (x[1] > N7 and ((a < b and a < x[1] - N7) or (a == b and a < x[1] - N7) or a == N7))}
                if kind == "gene":
                    assert hit <= {errors_by_gene.get(name)} or (code[rec.rsplit("_r", 1)[0]] != 11 and all(
                        h.startswith("recode") for h in hit)), (name, hit)
                elif kind == "junctions":
                    assert hit <= {"alt_start_upstream", "smorf", "recode_overlap"}, (rec, s, e, hit)
                elif kind in ("insertions", "prophages", "arrays"):
                    assert code[rec.rsplit("_r", 1)[0]] != 11 and all(h.startswith("recode") for h in hit)
                else:
                    raise AssertionError(f"{kind} {rec} [{s}, {e}) touched by S8 edits {hit}")
                n[f"{kind}_s8_edit"] += 1
                continue
            parts = idx[rec].parts(y[0], y[1], 1)
            got = "".join(_cut(t["cseq"][r["contig_id"]], cs, ce, cst) for r, cs, ce, cst, _, _ in
                          sorted(parts, key=lambda p: (p[4], not p[5])))
            if len(parts) == 1:
                assert got == bases, (kind, rec, s, e)
                n[f"{kind}_whole"] += 1
            else:
                assert sorted(got) == sorted(bases) and len(got) == len(bases)
                for r, cs, ce, cst, _, _ in parts:
                    assert _cut(t["cseq"][r["contig_id"]], cs, ce, cst) in bases
                n[f"{kind}_cut"] += 1
            n[f"{kind}_on_dropped"] += any(r["emitted"] == "false" for r, *_ in parts)
    print(f"\nearlier features on the contigs {dict(sorted(n.items()))}")
    assert n["gene_cut"] > 0 and n["elements_whole"] > 100_000


# --- runtime (§2.6) -----------------------------------------------------------------------------------------------

@needs_zombi2
def test_micro_stage_times(world, record_property):
    aw, _, times, _, _ = world
    for k, v in times.items():
        record_property(k.split()[0] + "_seconds", round(v, 2))
    total = sum(times.values())
    print("\n" + "; ".join(f"{k} {v:.2f} s" for k, v in times.items()) + f"; S1-S9 {total:.2f} s; "
          f"{aw.meta['counts']}; digest {aw.digest}")
    assert times["S8-S9 assembly"] < 120.0
