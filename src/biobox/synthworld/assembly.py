"""synthworld gene-call errors, genetic codes 4 and 25, assembly classes and fragmentation
(planning/tests/synthworld_spec.md §2.1 S8 and S9, §3.3 ``truth_callerr`` and ``truth_edges``, §3.2
``truth_contig``, §5.1; card F-02 session 7). Reads a world whose S1-S7 are built and writes
``world/assembly/``: the S7 host replicons with S8's edits applied and S9's cuts made, in **one** rewrite.

**One new frame.** S8's edits are applied to every replicon of ``fusions/`` in one pass, giving the *S8 frame*
(``edits.tsv`` is the lift from the S7 frame, in the shape of ``fusions/edits.tsv``). S9 cuts the S8 frame into
contigs and writes only the contigs: every base of the S8 frame is in exactly one contig, emitted
(``contigs/``) or not (``dropped/``), and ``contigs.tsv`` places each contig on the frame, so the S8 frame is
the "pre-rotation, forward replicon frame" of §3.3 ``truth_edges``. The S8 frame itself is not written
(``replicons.tsv`` gives its lengths, GC and BLAKE3-128).

**S8, call errors** (``callerr.counts``; ids ``CE0000`` ... in the order of ``ERROR_KINDS``). Each error goes to
a genome drawn uniformly among the code-11 genomes (stream ``S8/place``), then to a gene or gap of that genome
(stream ``S8/error/<id>``: the candidates in sorted order, permuted, the first one free). Only host genes of S7
(``source`` host, not across an origin) are used, at most one error per gene, and no two errors' footprints
(the gene, or the window or point they edit) overlap. An error's *gene strand* coordinate q is a base index of
the gene read 5'->3'.

- ``frameshift_del`` / ``frameshift_ins``: one base of a codon in the gene's middle third (codon index in
  [n/3, 2n/3), n codons, the codon shared with no other gene) is deleted, or a base (uniform over ACGT) is
  inserted before it. The intended model (``true_start``, ``true_end``) is the gene with the indel, its
  ``true_protein_uid`` the S7 protein (a sequencing indel in a real gene); ``expected_arch`` not_split.
- ``internal_stop``: one base of a middle-third codon is changed so that the codon reads as a stop (uniform over
  the codon and base choices that make one); the gene becomes a pseudogene (``feature_kind`` pseudogene, no
  protein_uid); ``expected_arch`` not_split.
- ``alt_start_downstream``: a gene starting ATG or GTG (its start codon shared with no other gene) with an ATG
  at codon index j in ``alt_downstream_codons`` has its start recoded to TTG (still initiator M: the protein is
  unchanged); ``alt_start`` is the competing ATG.
- ``alt_start_upstream``: k codons (uniform on ``alt_upstream_codons``) before the true start, in frame, an ATG
  is written, and every in-frame stop between has its first base changed T -> C, so the ORF reaches the planted
  start; the window holds no gene and no edited base lies on a protected base (an S5 element or motif, an S6
  insertion, an S7 edit or split, lifted to the S7 frame). The gene and its protein are unchanged; ``alt_start``
  is the planted ATG.
- ``smorf``: a block, pad + RBS (AGGAGG) + r bases + ATG + (a - 1) codons + TAA + pad (``smorf_pad`` background
  bases from the genome's intergenic composition, r uniform on ``intergenic.cassettes.rbs_spacer``, a uniform
  on ``smorf_aa``, residues uniform over the 20 amino acids, codons from the genome's usage), on a strand drawn
  with probability 1/2, is inserted at a point of an intergenic gap of at least ``smorf_min_gap`` bases, at least
  ``smorf_pad`` from its ends, next to no protected base. The smORF is a new gene ``<genome>_<id>`` of family
  ``smORF``.

**Genetic codes 4 and 25** (stream ``S8/codes``): the code-4 clade is a species branch with ``code4_clade_size``
extant genomes below it, the code-25 genome a tip outside it; every other genome is code 11 and only code-11
genomes get call errors. In a code-4 or code-25 genome every gene (prophage, Cas, fused, split and NOGD genes
included) is recoded so that it encodes the same protein under the genome's code (stream ``S8/recode/<genome>``,
genes in record and start order):
- a TGA stop becomes TAA (its middle base); when that base is also in the next gene (a -4 overlap ``ATGA``, the
  downstream start Z+TG), ``AAT`` is inserted between the stop's T and G (``ZTGA`` -> ``ZTAATGA``): the
  upstream gene now stops with TAA, the downstream one starts with ATG three bases later, and both proteins are
  unchanged (the junction becomes a -1 overlap);
- in code 4 an internal TGG (Trp) becomes TGA with probability ``code4_tgg_to_tga``; in code 25 an internal Gly
  codon becomes TGA with probability ``code25_gly_to_tga``; a codon with a base in another gene is left alone.
A gene that then holds an in-frame TGA has a ``callerr.tsv`` row of kind ``alt_code`` (``error_pos`` the first
TGA). All substitutions are single-base edits.

**expected_calls** are, in the S8 frame, what a code-11 reading of the error's gene gives: for alt starts and
smORFs the true model; for frameshifts, internal stops and alt-code genes the 5' piece (true start to the first
code-11 stop in its frame) and the 3' piece (from the first start codon after the previous in-frame stop to the
true stop, within the gene), each listed when it has at least ``min_call_codons`` sense codons (middle pieces
are not listed).

**S9, assembly** (stream ``S9/classes``: a permutation of the genomes, the first round(0.4 n) complete, the next
round(0.4 n) draft, the rest MAG, by ``assembly.class_weights``; then stream ``S9/genome/<genome>``):
- *complete*: one contig per replicon; a circular replicon is rotated to an origin inside one of its genes (a
  gene uniform, an offset uniform over its interior bases), so every complete circular contig has a gene across
  its origin (MR-04, §5.1 at least 5), and is circular (``is_circular``);
- *draft* (target N50 uniform over ``draft_n50``) and *MAG* (``mag_n50``): n = max(circular replicons,
  round(G x ``n50_factor`` / N50)) breaks over the genome's G bases: first one at the midpoint of each repeat
  copy (vehicle chassis families with two or more copies in the genome; copies of a family of equal length >=
  ``repeat_min_length`` at >= ``repeat_min_identity`` identity), a random subset of them if they exceed n, then
  uniform random breaks for the rest, and one more on any circular replicon still unbroken; a circular
  replicon with cuts c1 < ... < cm gives pieces [ci, ci+1) and [cm, c1 + N) (end > N), a linear one [0, c1) ...
  [cm, N);
- pieces under ``min_contig`` bases are not emitted (``drop_reason`` short); in a MAG, pieces are then dropped in
  a random order while the emitted share of G stays >= c (c uniform on ``mag_completeness``; ``drop_reason``
  mag_dropout);
- every contig is reverse-complemented with probability ``reverse_probability`` (``orientation`` -1).

Shape on disk (``build(world)``)::

    contigs/<genome_id>.fna, .fai   emitted contigs, records <genome_id>_c<kkkk> in contig order
    dropped/<genome_id>.fna, .fai   contigs not emitted
    genomes.tsv     genome_id genetic_code assembly_class target_n50 n50 target_completeness completeness length
                    emitted_length n_contigs n_emitted n_breaks n_repeat_breaks
    replicons.tsv   genome_id replicon record kind topology length s7_length n_edits gc contig_uid
    edits.tsv       edit_id callerr_id kind genome_id replicon record start end length new_start new_end
    genes.tsv       genome_id replicon copy_label family_id source feature_kind start end strand origin_spanning
                    protein_uid genetic_code callerr_id
    callerr.tsv     callerr_id kind genome_id replicon record true_gene_id true_start true_end strand error_pos
                    alt_start true_protein_uid expected_calls expected_arch genetic_code status
    contigs.tsv     contig_id genome_id replicon record length replicon_start replicon_end orientation
                    rotation_offset is_circular left_edge right_edge repeat_family_left repeat_family_right
                    n_partial_genes emitted drop_reason gc contig_uid
    gene_contigs.tsv genome_id copy_label contig_id start end strand partial_5 partial_3 origin_spanning emitted
    repeats.tsv     genome_id replicon copy_label family_id start end strand reason
    assembly.json   format, seed, the callerr and assembly presets and their digest, upstream digests, versions,
                    stream paths, counts, sha256 of every file, digest

Coordinates: ``edits``, ``genes``, ``callerr`` and ``repeats`` are on the S8 frame (0-based half-open on the +
strand of the replicon record; a feature across the origin has end > length); ``edits.start``/``end`` are on
the S7 frame (the lift: as ``fusions/edits.tsv``). ``contigs.replicon_start``/``replicon_end`` place a contig on
the S8 frame (end > length when it crosses the frame's origin; a complete circular contig spans [r, r + N) with
``rotation_offset`` r): the contig is ``frame[replicon_start:replicon_end]`` (circularly), reverse-complemented
when ``orientation`` is -1. ``left_edge`` is the edge at the contig's coordinate 0 and ``right_edge`` the edge
at its end (contig orientation): ``replicon_end``, ``break_repeat``, ``break_random`` or ``none`` (a complete
circular contig); ``repeat_family_*`` the repeat's family at a repeat break. ``gene_contigs`` has one row per
piece of a gene on a contig, in contig coordinates (strand on the contig; ``partial_5``/``partial_3`` the gene's
5'/3' end cut off by a contig edge; ``origin_spanning`` end > contig length on a circular contig; ``emitted``
the contig's), so a gene cut by a break has two rows and a gene on a dropped contig says so. ``callerr.status``:
whole (one emitted piece, not partial), partial, or dropped (no emitted piece). ``expected_calls`` items are
``start-end:strand``. ``genes.protein_uid`` is BLAKE3-128 of the protein under the genome's code (empty for a
pseudogene; a frameshift's is its intended protein). Floats at 6 significant digits. One genome is in memory
at a time. ``digest`` is BLAKE3-128 over the sorted (path, sha256) list of every file but assembly.json.
"""

from __future__ import annotations

import bisect
import collections
import json
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

import numpy as np
from blake3 import blake3

from biobox import ids
from biobox.synthworld import backbone, codons, fusions, intergenic, mobile, seeds, sequences

FORMAT = "synthworld-assembly/1"
ERROR_KINDS = ("frameshift_ins", "frameshift_del", "alt_start_upstream", "alt_start_downstream", "smorf",
               "internal_stop")
CLASSES = ("complete", "draft", "mag")
STARTS = ("ATG", "GTG", "TTG")
STOPS11 = ("TAA", "TAG", "TGA")
CODES = (11, 4, 25)
RBS = intergenic.RBS
q6 = sequences.q6
revcomp = intergenic.revcomp
_COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


# --- genetic codes -----------------------------------------------------------------------------------------

def code_table(code: int) -> dict[str, str]:
    """NCBI transl_table 11, 4 (TGA = W) or 25 (TGA = G)."""
    t = dict(codons.CODE11)
    if code == 4:
        t["TGA"] = "W"
    elif code == 25:
        t["TGA"] = "G"
    elif code != 11:
        raise ValueError(f"genetic code {code} is not one of {CODES}")
    return t


_TABLE = {c: code_table(c) for c in CODES}


def translate(cds: str, code: int = 11) -> str:
    """Codon by codon under ``code``, the first codon read as M when it is ATG, GTG or TTG; stops read ``*``."""
    t = _TABLE[code]
    out = [t[cds[i:i + 3]] for i in range(0, len(cds) - len(cds) % 3, 3)]
    if out and cds[:3] in STARTS:
        out[0] = "M"
    return "".join(out)


def orf_protein(cds: str, code: int = 11) -> str | None:
    """The protein (initiator M, no stop) when ``cds`` reads as one ORF under ``code``, else None."""
    if len(cds) % 3 or len(cds) < 6 or cds[:3] not in STARTS:
        return None
    p = translate(cds, code)
    if not p.endswith("*") or "*" in p[:-1]:
        return None
    return p[:-1]


def cut(seq: str, s: int, e: int, strand: int = 1) -> str:
    """Bases [s, e) of a circular record (e > len wraps), reverse-complemented on strand -1."""
    n = len(seq)
    x = seq[s:] + seq[:e - n] if e > n else seq[s:e]
    return x if strand == 1 else revcomp(x)


def code11_calls(cds: str, min_codons: int) -> list[tuple[int, int]]:
    """The 5' and 3' code-11 ORF pieces of a gene read 5'->3' (module docstring), as offsets in ``cds``."""
    L = len(cds)
    out = []
    for i in range(0, L - 2, 3):
        if cds[i:i + 3] in STOPS11:
            if i // 3 >= min_codons:
                out.append((0, i + 3))
            break
    q = None
    for j in range(L - 6, -1, -3):
        if cds[j:j + 3] in STOPS11:
            q = j
            break
    lo = q + 3 if q is not None else (L - 3) % 3
    for j in range(lo, L - 3, 3):
        if cds[j:j + 3] in STARTS:
            if (L - 3 - j) // 3 >= min_codons and (j, L) not in out:
                out.append((j, L))
            break
    return sorted(out)


def _piece_coords(s: int, e: int, strand: int, a: int, b: int, N: int) -> tuple[int, int]:
    x, y = (s + a, s + b) if strand == 1 else (e - b, e - a)
    if x >= N:
        x, y = x - N, y - N
    return x, y


# --- the lift from the S7 frame ------------------------------------------------------------------------------

class Lift(fusions.Edits):
    """`fusions.Edits` (replacements [a, b) -> l bases, a = b an insertion) with a vectorized touch test."""

    def _range(self, x: int, y: int) -> tuple[int, int]:
        return int(np.searchsorted(self.b, x, side="right")), int(np.searchsorted(self.a, y, side="left"))

    def touching(self, s: int, e: int) -> list[int]:
        """Indices (sorted edit order) of the edits that touch [s, e) (an edit [a, b) with a < e and b > s, or
        an insertion strictly inside)."""
        s, e = int(s), int(e)
        parts = [(s, e)] if e <= self.N else [(s, self.N + 1), (-1, e - self.N)]
        out: list[int] = []
        for x, y in parts:
            i0, i1 = self._range(x, y)
            out += range(i0, max(i0, i1))
        return out

    def touched(self, s: int, e: int) -> bool:
        return bool(self.touching(s, e))

    def start(self, x: int) -> int:
        """Where position x of the old record is (x not inside a replacement that changes length); an insertion
        at x puts its bases before x."""
        return int(x) + int(self.cum[np.searchsorted(self.b, int(x), side="right")])


# --- placing frame features on contigs ----------------------------------------------------------------------

class Piece(NamedTuple):
    replicon: int
    a: int                 # [a, b) on the S8 frame, b - a <= N, b > N across the origin
    b: int
    circular: bool         # a complete circular contig
    left: str              # edge at a (frame direction)
    right: str             # edge at b
    fam_left: str
    fam_right: str


def pieces(N: int, circular: bool, cuts: dict[int, tuple[str, str]], rotation: int | None = None) -> list[Piece]:
    """The pieces of a replicon of N bases, frame order: ``cuts`` maps a cut position to (edge label, repeat
    family); a circular replicon without cuts is one circular piece from ``rotation``."""
    pos = sorted(cuts)
    if circular:
        if not pos:
            r = int(rotation or 0)
            return [Piece(0, r, r + N, True, "none", "none", "", "")]
        out = []
        for i, c in enumerate(pos):
            d = pos[i + 1] if i + 1 < len(pos) else pos[0] + N
            out.append(Piece(0, c, d, False, cuts[c][0], cuts[d % N][0], cuts[c][1], cuts[d % N][1]))
        return out
    bounds = [0] + [c for c in pos if 0 < c < N] + [N]
    lab = dict(cuts)
    lab[0] = lab[N] = ("replicon_end", "")
    return [Piece(0, x, y, False, lab[x][0], lab[y][0], lab[x][1], lab[y][1]) for x, y in zip(bounds, bounds[1:])]


def place(ps: list[Piece], N: int, circular: bool, s: int, e: int) -> list[tuple[int, int, int, bool, bool]]:
    """The parts of frame feature [s, e) on pieces ``ps`` (frame order): (piece index, x, y, clipped at x,
    clipped at y), x and y offsets from the piece's a in frame direction (y > piece length on a circular
    piece, across its origin)."""
    if len(ps) == 1 and ps[0].circular:
        x = (s - ps[0].a) % N
        return [(0, x, x + (e - s), False, False)]
    starts = [p.a for p in ps]
    out = []
    for k in ((-N, 0, N) if circular else (0,)):
        fs, fe = s + k, e + k
        i = max(0, bisect.bisect_right(starts, fs) - 1)
        while i < len(ps) and ps[i].a < fe:
            p = ps[i]
            if p.b > fs:
                x, y = max(fs, p.a), min(fe, p.b)
                out.append((i, x - p.a, y - p.a, fs < p.a, fe > p.b))
            i += 1
    return sorted(out)


def to_contig(L: int, orientation: int, x: int, y: int, strand: int) -> tuple[int, int, int]:
    """Offsets [x, y) of a piece of length L (frame direction) in contig coordinates, with the strand."""
    if orientation == 1:
        return x, y, strand
    s = (L - y) % L
    return s, s + (y - x), -strand


# --- tables --------------------------------------------------------------------------------------------------

class GenomeRow(NamedTuple):
    genome_id: str
    genetic_code: int
    assembly_class: str
    target_n50: int
    n50: int
    target_completeness: float
    completeness: float
    length: int
    emitted_length: int
    n_contigs: int
    n_emitted: int
    n_breaks: int
    n_repeat_breaks: int


class RepliconRow(NamedTuple):
    genome_id: str
    replicon: int
    record: str
    kind: str
    topology: str
    length: int
    s7_length: int
    n_edits: int
    gc: float
    contig_uid: str


class EditRow(NamedTuple):
    edit_id: str
    callerr_id: str
    kind: str
    genome_id: str
    replicon: int
    record: str
    start: int
    end: int
    length: int
    new_start: int
    new_end: int


class GeneRow(NamedTuple):
    genome_id: str
    replicon: int
    copy_label: str
    family_id: str
    source: str
    feature_kind: str
    start: int
    end: int
    strand: int
    origin_spanning: bool
    protein_uid: str
    genetic_code: int
    callerr_id: str


class CallerrRow(NamedTuple):
    callerr_id: str
    kind: str
    genome_id: str
    replicon: int
    record: str
    true_gene_id: str
    true_start: int
    true_end: int
    strand: int
    error_pos: int
    alt_start: int
    true_protein_uid: str
    expected_calls: tuple[str, ...]
    expected_arch: str
    genetic_code: int
    status: str


class ContigRow(NamedTuple):
    contig_id: str
    genome_id: str
    replicon: int
    record: str
    length: int
    replicon_start: int
    replicon_end: int
    orientation: int
    rotation_offset: int
    is_circular: bool
    left_edge: str
    right_edge: str
    repeat_family_left: str
    repeat_family_right: str
    n_partial_genes: int
    emitted: bool
    drop_reason: str
    gc: float
    contig_uid: str


class GenePieceRow(NamedTuple):
    genome_id: str
    copy_label: str
    contig_id: str
    start: int
    end: int
    strand: int
    partial_5: bool
    partial_3: bool
    origin_spanning: bool
    emitted: bool


class RepeatRow(NamedTuple):
    genome_id: str
    replicon: int
    copy_label: str
    family_id: str
    start: int
    end: int
    strand: int
    reason: str


_TABLES = {"genomes.tsv": GenomeRow, "replicons.tsv": RepliconRow, "edits.tsv": EditRow, "genes.tsv": GeneRow,
           "callerr.tsv": CallerrRow, "contigs.tsv": ContigRow, "gene_contigs.tsv": GenePieceRow,
           "repeats.tsv": RepeatRow}


def check_config(ce: dict, asm: dict) -> None:
    if set(ce["counts"]) != set(ERROR_KINDS) or min(int(v) for v in ce["counts"].values()) < 0:
        raise ValueError(f"callerr.counts gives a count >= 0 for each of {ERROR_KINDS}")
    if sum(int(v) for v in ce["counts"].values()) > 9999:
        raise ValueError("at most 9999 call errors (ids CE0000 ...)")
    for k in ("alt_upstream_codons", "alt_downstream_codons", "smorf_aa", "code4_clade_size"):
        lo, hi = (int(x) for x in ce[k])
        if not 1 <= lo <= hi:
            raise ValueError(f"callerr.{k} must be 1 <= lo <= hi, got {ce[k]}")
    if int(ce["alt_upstream_codons"][0]) < 6:
        raise ValueError("callerr.alt_upstream_codons starts at 6, beyond an RBS and its spacer")
    if int(ce["smorf_aa"][0]) < 2 or int(ce["min_gene_codons"]) < 9:
        raise ValueError("callerr.smorf_aa starts at 2 and min_gene_codons at 9")
    for k in ("code4_tgg_to_tga", "code25_gly_to_tga"):
        if not 0.0 <= float(ce[k]) <= 1.0:
            raise ValueError(f"callerr.{k} is a probability")
    w = [float(x) for x in asm["class_weights"]]
    if len(w) != 3 or min(w) < 0 or abs(sum(w) - 1.0) > 1e-9:
        raise ValueError("assembly.class_weights: three weights (complete, draft, mag) summing to 1")
    lo, hi = (float(x) for x in asm["mag_completeness"])
    if not 0.0 < lo <= hi <= 1.0:
        raise ValueError("assembly.mag_completeness must lie in (0, 1]")
    if min(int(x) for x in asm["draft_n50"]) < 1 or int(asm["mag_n50"]) < 1 or int(asm["min_contig"]) < 1:
        raise ValueError("N50 targets and min_contig are positive")


@dataclass
class AssemblyWorld:
    seed: int
    genomes: list[GenomeRow]
    replicons: list[RepliconRow]
    digest: str
    out: Path
    meta: dict = field(default_factory=dict, repr=False)


# --- S8 planning, one genome -----------------------------------------------------------------------------------

def _codon_pos(g: dict, c: int, N: int) -> tuple[int, int, int]:
    """+ strand positions of codon c of gene g, in reading order."""
    s, e, st = g["s"], g["e"], g["st"]
    if st == 1:
        return tuple((s + 3 * c + i) % N for i in range(3))   # type: ignore[return-value]
    return tuple((e - 1 - 3 * c - i) % N for i in range(3))   # type: ignore[return-value]


def _read(seq: str, pos, st: int) -> str:
    x = "".join(seq[p] for p in pos)
    return x if st == 1 else "".join(_COMP[b] for b in x)


@dataclass
class _Rep:
    rep: int
    record: str
    seq: str
    N: int
    kind: str
    topology: str
    counts: np.ndarray            # genes covering each base
    used: np.ndarray              # footprints of planned errors
    prot: np.ndarray | None = None
    edits: list = field(default_factory=list)    # (a, b, block, kind, owner)


def _counts(N: int, genes: list[dict]) -> np.ndarray:
    d = np.zeros(N + 1, dtype=np.int32)
    for g in genes:
        s, e = g["s"], g["e"]
        if e <= N:
            d[s] += 1
            d[e] -= 1
        else:
            d[s] += 1
            d[N] -= 1
            d[0] += 1
            d[e - N] -= 1
    return np.cumsum(d)[:N]


def _subst(rp: _Rep, pos, old: str, new: str, st: int, kind: str, owner: str) -> None:
    """Single-base edits turning codon ``old`` (read on strand st at ``pos``) into ``new``."""
    for p, x, y in zip(pos, old, new):
        if x != y:
            rp.edits.append((p, p + 1, y if st == 1 else _COMP[y], kind, owner))


def _middle_codons(g: dict, rp: _Rep) -> list[int]:
    n = (g["e"] - g["s"]) // 3
    return [c for c in range(max(1, n // 3), 2 * n // 3)
            if all(rp.counts[p] == 1 for p in _codon_pos(g, c, rp.N))]


def _free(rp: _Rep, s: int, e: int) -> bool:
    return not rp.used[s:e].any()


def _plan_error(eid: str, kind: str, rng, reps: dict[int, _Rep], cands: list[dict], cfg: dict, ctx: dict) -> dict:
    """Plan one error (module docstring); marks its footprint and appends its edits."""
    if kind == "smorf":
        return _plan_smorf(eid, rng, reps, cfg, ctx)
    k_up = int(rng.integers(int(cfg["alt_upstream_codons"][0]), int(cfg["alt_upstream_codons"][1]) + 1)) \
        if kind == "alt_start_upstream" else 0
    for i in rng.permutation(len(cands)):
        g = cands[int(i)]
        rp = reps[g["rep"]]
        s, e, st, N = g["s"], g["e"], g["st"], rp.N
        if g["label"] in ctx["used_genes"] or not _free(rp, s, e):
            continue
        n = (e - s) // 3
        if kind in ("frameshift_ins", "frameshift_del", "internal_stop"):
            if n < int(cfg["min_gene_codons"]):
                continue
            mid = _middle_codons(g, rp)
            if not mid:
                continue
            if kind == "internal_stop":
                opts = []
                for c in mid:
                    pos = _codon_pos(g, c, N)
                    old = _read(rp.seq, pos, st)
                    for j in range(3):
                        for b in "ACGT":
                            new = old[:j] + b + old[j + 1:]
                            if b != old[j] and new in STOPS11:
                                opts.append((c, j, b))
                c, j, b = opts[int(rng.integers(len(opts)))]
                pos = _codon_pos(g, c, N)
                _subst(rp, pos, _read(rp.seq, pos, st), _read(rp.seq, pos, st)[:j] + b
                       + _read(rp.seq, pos, st)[j + 1:], st, kind, eid)
                info = {"q": 3 * c + j, "err": pos[j]}
            else:
                c = mid[int(rng.integers(len(mid)))]
                j = int(rng.integers(3))
                q = 3 * c + j
                if kind == "frameshift_del":
                    p = s + q if st == 1 else e - 1 - q
                    rp.edits.append((p, p + 1, "", kind, eid))
                    info = {"q": q, "err": p, "len_delta": -1}
                else:
                    b = "ACGT"[int(rng.integers(4))]
                    p = s + q if st == 1 else e - q
                    rp.edits.append((p, p, b if st == 1 else _COMP[b], kind, eid))
                    info = {"q": q, "err": p, "len_delta": 1, "base": b}
            rp.used[s:e] = True
        elif kind == "alt_start_downstream":
            p0 = _codon_pos(g, 0, N)
            old = _read(rp.seq, p0, st)
            if old not in ("ATG", "GTG") or any(rp.counts[p] != 1 for p in p0):
                continue
            lo, hi = (int(x) for x in cfg["alt_downstream_codons"])
            js = [j for j in range(lo, min(hi, n - 2) + 1) if _read(rp.seq, _codon_pos(g, j, N), st) == "ATG"]
            if not js:
                continue
            j = js[0]
            _subst(rp, p0, old, "TTG", st, kind, eid)
            alt = _codon_pos(g, j, N)
            info = {"err": min(p0), "alt": min(alt), "j": j}
            rp.used[s:e] = True
        elif kind == "alt_start_upstream":
            lo, hi = (s - 3 * k_up, s) if st == 1 else (e, e + 3 * k_up)
            if lo < 1 or hi > N - 1 or rp.counts[lo:hi].any() or rp.used[lo:hi].any():
                continue
            # codon m of the window, m = 0 the planted start (5' end), m = k_up the true start
            def cpos(m: int) -> tuple[int, int, int]:
                if st == 1:
                    return (lo + 3 * m, lo + 3 * m + 1, lo + 3 * m + 2)
                return (hi - 1 - 3 * m, hi - 2 - 3 * m, hi - 3 - 3 * m)
            plan = [(cpos(0), _read(rp.seq, cpos(0), st), "ATG")]
            for m in range(1, k_up):
                cod = _read(rp.seq, cpos(m), st)
                if cod in STOPS11:
                    plan.append((cpos(m), cod, "C" + cod[1:]))
            touched = [p for pos, old, new in plan for p, x, y in zip(pos, old, new) if x != y]
            if rp.prot is not None and any(rp.prot[p] for p in touched):
                continue
            for pos, old, new in plan:
                _subst(rp, pos, old, new, st, kind, eid)
            info = {"err": min(cpos(0)), "alt": min(cpos(0)), "k": k_up}
            rp.used[lo:hi] = True
            rp.used[s:e] = True
        else:
            raise ValueError(kind)
        ctx["used_genes"].add(g["label"])
        return dict(info, eid=eid, kind=kind, rep=g["rep"], gene=g)
    raise ValueError(f"{eid} ({kind}): no free candidate in {ctx['gid']}")


def _plan_smorf(eid: str, rng, reps: dict[int, _Rep], cfg: dict, ctx: dict) -> dict:
    a_lo, a_hi = (int(x) for x in cfg["smorf_aa"])
    n_aa = int(rng.integers(a_lo, a_hi + 1))
    st = 1 if rng.random() < 0.5 else -1
    r = int(rng.integers(ctx["rbs"][0], ctx["rbs"][1] + 1))
    pad, min_gap = int(cfg["smorf_pad"]), int(cfg["smorf_min_gap"])
    cod, prot = fusions._draw_codons(rng, codons.AA, n_aa - 1, ctx["usage"])
    comp = ctx["comp"]
    block = (fusions._background(rng, pad, comp) + RBS + fusions._background(rng, r, comp) + "ATG" + cod + "TAA"
             + fusions._background(rng, pad, comp))
    off, glen = pad + len(RBS) + r, 3 * (n_aa + 1)
    gaps = []
    for rep in sorted(reps):
        rp = reps[rep]
        z = np.concatenate([[1], (rp.counts == 0).astype(np.int8) ^ 1, [1]])
        d = np.diff(z)
        for x, y in zip(np.flatnonzero(d == -1), np.flatnonzero(d == 1)):
            if y - x >= min_gap:
                gaps.append((rep, int(x), int(y)))
    for i in rng.permutation(len(gaps)):
        rep, x, y = gaps[int(i)]
        rp = reps[rep]
        lo, hi = x + pad, y - pad
        pts = np.arange(lo, hi + 1)
        bad = rp.used[pts - 1] | rp.used[np.minimum(pts, rp.N - 1)]
        if rp.prot is not None:
            bad |= rp.prot[pts - 1] | rp.prot[np.minimum(pts, rp.N - 1)]
        ok = pts[~bad]
        if not len(ok):
            continue
        p = int(ok[int(rng.integers(len(ok)))])
        blk = block if st == 1 else revcomp(block)
        rp.edits.append((p, p, blk, "smorf", eid))
        rp.used[max(0, p - 1):p + 1] = True
        g_off = off if st == 1 else len(block) - off - glen
        return {"eid": eid, "kind": "smorf", "rep": rep, "gene": None, "point": p, "strand": st,
                "g_off": g_off, "glen": glen, "protein": "M" + prot, "err": p}
    raise ValueError(f"{eid} (smorf): no gap in {ctx['gid']}")


def _recode(rp: _Rep, genes: list[dict], code: int, rng, cfg: dict, mods: dict) -> None:
    """Recode every gene of a replicon of a code-4 or code-25 genome (module docstring)."""
    N, seq = rp.N, rp.seq
    f = float(cfg["code4_tgg_to_tga"] if code == 4 else cfg["code25_gly_to_tga"])
    by_start = {(g["s"] % N, g["st"]): g for g in genes}
    for g in genes:
        st, n = g["st"], (g["e"] - g["s"]) // 3
        # internal codons
        for c in range(1, n - 1):
            pos = _codon_pos(g, c, N)
            cod = _read(seq, pos, st)
            hit = cod == "TGG" if code == 4 else cod[:2] == "GG"
            if not hit:
                continue
            u = rng.random()
            if u < f and all(rp.counts[p] == 1 for p in pos):
                _subst(rp, pos, cod, "TGA", st, "recode_codon", g["label"])
        # the stop
        pos = _codon_pos(g, n - 1, N)
        cod = _read(seq, pos, st)
        if cod != "TGA":
            continue
        if rp.counts[pos[1]] == 1:
            _subst(rp, pos, cod, "TAA", st, "recode_stop", g["label"])
            continue
        # a -4 overlap: the downstream gene starts at the stop's first base - 1 (on its strand)
        if st == 1:
            x = (g["e"] - 4) % N
            h = by_start.get((x, 1))
            point = (x + 2) % N or N
            if h is None:
                raise AssertionError(f"{g['label']}: a shared TGA stop that is not a -4 overlap")
            rp.edits.append((point, point, "AAT", "recode_overlap", g["label"]))
            mods.setdefault(h["label"], {})["start_delta"] = 3
        else:
            # downstream gene h lies to the left, h.end = g.start + 4
            h = next((y for y in genes if y["st"] == -1 and y["e"] % N == (g["s"] + 4) % N and y is not g), None)
            if h is None:
                raise AssertionError(f"{g['label']}: a shared TGA stop that is not a -4 overlap")
            point = (g["s"] + 2) % N or N
            rp.edits.append((point, point, "ATT", "recode_overlap", g["label"]))
            mods.setdefault(g["label"], {})["start_delta"] = 3


# --- protected bases (S5 elements, S6 insertions, S7 edits) on the S7 frame -------------------------------------

def _protected(rp: _Rep, s5_len: int, s6_len: int, elements: list[dict], ins6: list[dict], ed7: list[dict],
               spans7: list[tuple[int, int]]) -> np.ndarray:
    N = rp.N
    mask = np.zeros(N, dtype=bool)
    lift6 = mobile.Lift(s5_len, [(int(r["at"]), int(r["length"])) for r in ins6])
    lift7 = Lift(s6_len, [(int(r["start"]), int(r["end"]), int(r["length"])) for r in ed7])
    spans = list(spans7)
    for r in ed7:
        spans.append((int(r["new_start"]), int(r["new_end"])))
    for r in ins6:
        s, e = int(r["start"]), int(r["end"])
        try:
            spans.append(lift7.feature(s, e))
        except ValueError:
            pass
    for r in elements:
        s, e = int(r["start"]), int(r["end"])
        try:
            spans.append(lift7.feature(*lift6.feature(s, e)))
        except ValueError:
            pass
    for s, e in spans:
        if e <= N:
            mask[s:e] = True
        else:
            mask[s:] = True
            mask[:e - N] = True
    return mask


class _Stream:
    """A TSV sorted by genome_id, read forward one genome at a time on request."""

    def __init__(self, path: Path):
        self.it = intergenic._genome_rows(path)
        self.pending = next(self.it, None)

    def rows(self, gid: str) -> list[dict]:
        while self.pending is not None and self.pending[0] < gid:
            self.pending = next(self.it, None)
        if self.pending is not None and self.pending[0] == gid:
            out = self.pending[1]
            self.pending = next(self.it, None)
            return out
        return []


# --- S9, one genome -------------------------------------------------------------------------------------------

def n50(lengths) -> int:
    xs = sorted((int(x) for x in lengths), reverse=True)
    half, acc = sum(xs) / 2.0, 0
    for x in xs:
        acc += x
        if acc >= half:
            return x
    return 0


def _repeat_loci(genes: list[GeneRow], seqs: dict[int, str], fam_kind: dict[str, str], cfg: dict) -> list[tuple]:
    by_fam = collections.defaultdict(list)
    for g in genes:
        if g.family_id and g.source != "smorf":
            by_fam[g.family_id].append(g)
    out = {}
    L, ident = int(cfg["repeat_min_length"]), float(cfg["repeat_min_identity"])
    for fam in sorted(by_fam):
        gs = by_fam[fam]
        if len(gs) < 2:
            continue
        if fam_kind.get(fam) == "vehicle_chassis":
            for g in gs:
                out[g.copy_label] = (g, "chassis")
            continue
        arr = {}
        for g in gs:
            if g.end - g.start >= L:
                arr[g.copy_label] = np.frombuffer(cut(seqs[g.replicon], g.start, g.end, g.strand).encode(), np.uint8)
        labs = sorted(arr)
        for i, x in enumerate(labs):
            for y in labs[i + 1:]:
                if len(arr[x]) == len(arr[y]) and float((arr[x] == arr[y]).mean()) >= ident:
                    for lab in (x, y):
                        g = next(z for z in gs if z.copy_label == lab)
                        out.setdefault(lab, (g, "duplicate"))
    return [out[k] for k in sorted(out, key=lambda k: (out[k][0].replicon, out[k][0].start, k))]


# --- build -------------------------------------------------------------------------------------------------------

def build(world) -> AssemblyWorld:
    """Write ``world/assembly/`` (module docstring) from a world whose S1-S7 are built."""
    world = Path(world)
    docs = {}
    for d, n, mod in (("backbone", "backbone", backbone), ("sequences", "sequences", sequences),
                      ("replicons", "replicons", intergenic), ("mobile", "mobile", mobile),
                      ("fusions", "fusions", fusions)):
        docs[d] = json.loads((world / d / f"{n}.json").read_text(encoding="utf-8"))
        if docs[d]["format"] != mod.FORMAT:
            raise ValueError(f"{world}/{d}: format {docs[d]['format']!r}, expected {mod.FORMAT!r}")
    bb = docs["backbone"]
    preset = bb["preset"]
    ce, asm = preset["callerr"], preset["assembly"]
    check_config(ce, asm)
    seed = int(bb["seed"])
    out = world / "assembly"
    if out.exists():
        raise FileExistsError(f"{out} exists; the assembly stage is written once")
    tmp = world / "assembly.partial"
    if tmp.exists():
        shutil.rmtree(tmp)
    (tmp / "contigs").mkdir(parents=True)
    (tmp / "dropped").mkdir()

    branches = backbone.species_branches(world)
    host_genomes = mobile._tsv(world / "backbone" / "genomes.tsv")
    gid_of = {r["species_node"]: r["genome_id"] for r in host_genomes}
    node_of = {g: n for n, g in gid_of.items()}
    all_gids = sorted(node_of)
    tree = mobile.species_tree(branches, gid_of)
    lin_seq = sequences.lineage_model(branches, docs["sequences"]["sequences_preset"], seed)
    lin_igr = intergenic.lineages(world, preset["intergenic"]["composition"])
    fam_kind = {r["family_id"]: r["kind"] for r in mobile._tsv(world / "backbone" / "families.tsv")}
    rbs = tuple(int(x) for x in preset["intergenic"]["cassettes"]["rbs_spacer"])

    # --- genetic codes (S8/codes)
    rng = seeds.generator(seed, "S8/codes")
    lo, hi = (int(x) for x in ce["code4_clade_size"])
    cand = [x for x in tree.labels[1:] if lo <= len(tree.below[x]) <= hi]
    if not cand:
        raise ValueError(f"no species branch has {lo} to {hi} extant genomes below it")
    node4 = cand[int(rng.integers(len(cand)))]
    code_of = {g: 11 for g in all_gids}
    for tip in tree.below[node4]:
        code_of[gid_of[tip]] = 4
    rest = [g for g in all_gids if code_of[g] == 11]
    gid25 = rest[int(rng.integers(len(rest)))]
    code_of[gid25] = 25

    # --- error placement (S8/place)
    rng = seeds.generator(seed, "S8/place")
    eligible = [g for g in all_gids if code_of[g] == 11]
    errors_of: dict[str, list[tuple[str, str]]] = collections.defaultdict(list)
    k = 0
    for kind in ERROR_KINDS:
        for _ in range(int(ce["counts"][kind])):
            errors_of[eligible[int(rng.integers(len(eligible)))]].append((f"CE{k:04d}", kind))
            k += 1

    # --- assembly classes (S9/classes)
    rng = seeds.generator(seed, "S9/classes")
    w = [float(x) for x in asm["class_weights"]]
    n_g = len(all_gids)
    n_c, n_d = int(round(w[0] * n_g)), int(round(w[1] * n_g))
    perm = rng.permutation(n_g)
    class_of = {}
    for rank, i in enumerate(perm):
        class_of[all_gids[int(i)]] = "complete" if rank < n_c else "draft" if rank < n_c + n_d else "mag"

    # --- the earlier frames' small tables
    ins6 = collections.defaultdict(list)
    for r in mobile._tsv(world / "mobile" / "insertions.tsv"):
        ins6[r["record"]].append(r)
    ed7 = collections.defaultdict(list)
    for r in mobile._tsv(world / "fusions" / "edits.tsv"):
        ed7[r["record"]].append(r)
    spans7 = collections.defaultdict(list)
    for r in mobile._tsv(world / "fusions" / "fusions.tsv"):
        spans7[r["record"]].append((int(r["start"]), int(r["end"])))
    len5 = {r["record"]: int(r["length"]) for r in mobile._tsv(world / "replicons" / "replicons.tsv")}
    len6 = {r["record"]: int(r["length"]) for r in mobile._tsv(world / "mobile" / "replicons.tsv")}
    reps7 = collections.defaultdict(list)
    for r in mobile._tsv(world / "fusions" / "replicons.tsv"):
        reps7[r["genome_id"]].append(r)
    elements = _Stream(world / "replicons" / "elements.tsv")

    names = ("genomes.tsv", "replicons.tsv", "edits.tsv", "genes.tsv", "callerr.tsv", "contigs.tsv",
             "gene_contigs.tsv", "repeats.tsv")
    wr = {n: mobile._Rows(tmp / n, _TABLES[n]) for n in names}
    counts = collections.Counter()
    genome_rows, rep_rows = [], []
    for gid, rows in mobile._by_genome(world / "fusions" / "genes.tsv"):
        code, node = code_of[gid], node_of[gid]
        recs = dict(intergenic._fasta_records(world / "fusions" / "fasta" / f"{gid}.fna"))
        reps: dict[int, _Rep] = {}
        genes_by_rep = collections.defaultdict(list)
        for r in rows:
            genes_by_rep[int(r["replicon"])].append({
                "row": r, "rep": int(r["replicon"]), "label": r["copy_label"], "s": int(r["start"]),
                "e": int(r["end"]), "st": int(r["strand"])})
        for rr in sorted(reps7[gid], key=lambda r: int(r["replicon"])):
            rep, rec = int(rr["replicon"]), rr["record"]
            seq = recs[rec]
            N = len(seq)
            reps[rep] = _Rep(rep, rec, seq, N, rr["kind"], rr["topology"], _counts(N, genes_by_rep[rep]),
                             np.zeros(N, dtype=bool))
        my_errors = errors_of.get(gid, [])
        el_rows = elements.rows(gid)
        mods: dict[str, dict] = {}
        planned: list[dict] = []
        if any(k in ("alt_start_upstream", "smorf") for _, k in my_errors):
            el_by = collections.defaultdict(list)
            for r in el_rows:
                el_by[f"{gid}_r{r['replicon']}"].append(r)
            for rp in reps.values():
                rp.prot = _protected(rp, len5[rp.record], len6[rp.record], el_by[rp.record], ins6[rp.record],
                                     ed7[rp.record], spans7[rp.record])
        if my_errors:
            cands = [g for rep in sorted(genes_by_rep) for g in genes_by_rep[rep]
                     if g["row"]["source"] == "host" and g["row"]["origin_spanning"] == "false"]
            cands.sort(key=lambda g: (g["rep"], g["s"], g["label"]))
            li = lin_seq.index[node]
            ctx = {"gid": gid, "used_genes": set(), "rbs": rbs, "usage": lin_seq.cum[li],
                   "comp": lin_igr.comp_cum[lin_igr.index[node]]}
            for eid, kind in my_errors:
                pe = _plan_error(eid, kind, seeds.generator(seed, f"S8/error/{eid}"), reps, cands, ce, ctx)
                planned.append(pe)
                if pe["gene"] is not None:
                    m = mods.setdefault(pe["gene"]["label"], {})
                    m["callerr"] = eid
                    m["len_delta"] = pe.get("len_delta", 0)
                    if kind == "internal_stop":
                        m["kind"] = "pseudogene"
        if code != 11:
            rr_ = seeds.generator(seed, f"S8/recode/{gid}")
            for rep in sorted(reps):
                gs = sorted(genes_by_rep[rep], key=lambda g: (g["s"], g["label"]))
                _recode(reps[rep], gs, code, rr_, ce, mods)

        # --- apply S8's edits: the S8 frame
        seqs8: dict[int, str] = {}
        genes8: list[GeneRow] = []
        callerr: list[CallerrRow] = []
        for rep in sorted(reps):
            rp = reps[rep]
            ed = sorted(rp.edits, key=lambda x: (x[0], x[1]))
            new, starts = fusions.apply_edits(rp.seq, [(a, b, blk) for a, b, blk, _, _ in ed])
            lift = Lift(rp.N, [(a, b, len(blk)) for a, b, blk, _, _ in ed])
            if len(new) != lift.new_length:
                raise AssertionError(f"{rp.record}: edited length disagrees with the lift")
            N8 = len(new)
            for i, ((a, b, blk, kind, owner), b0) in enumerate(zip(ed, starts)):
                wr["edits.tsv"].row(EditRow(f"{rp.record}_e{i:05d}", owner if owner.startswith("CE") else "",
                                            kind, gid, rep, rp.record, a, b, len(blk), b0, b0 + len(blk)))
            touched_by = collections.defaultdict(set)
            for i, (a, b, blk, kind, owner) in enumerate(ed):
                touched_by[owner].add(i)
            out_genes = []
            for g in genes_by_rep[rep]:
                r, s, e = g["row"], g["s"], g["e"]
                m = mods.get(g["label"], {})
                own = set(lift.touching(s, e))
                if not own and not m:
                    s2, e2 = lift.feature(s, e)
                else:
                    s2 = lift.start(s) + m.get("start_delta", 0)
                    e2 = s2 + (e - s) + m.get("len_delta", 0)
                    if s2 >= N8:
                        s2, e2 = s2 - N8, e2 - N8
                kind_ = m.get("kind", "CDS")
                uid = r["protein_uid"] if kind_ == "CDS" else ""
                if own or m:
                    _check_gene(new, rp.seq, g, s2, e2, code, m, r)
                out_genes.append(GeneRow(gid, rep, g["label"], r["family_id"], r["source"], kind_, s2, e2, g["st"],
                                         e2 > N8, uid, code, m.get("callerr", "")))
            for pe in planned:
                if pe["rep"] != rep or pe["kind"] != "smorf":
                    continue
                i = next(i for i, x in enumerate(ed) if x[4] == pe["eid"])
                s2 = starts[i] + pe["g_off"]
                e2 = s2 + pe["glen"]
                cds = cut(new, s2, e2, pe["strand"])
                if orf_protein(cds) != pe["protein"]:
                    raise AssertionError(f"{pe['eid']}: the smORF does not encode its protein")
                lab = f"{gid}_{pe['eid']}"
                out_genes.append(GeneRow(gid, rep, lab, "smORF", "smorf", "smORF", s2, e2, pe["strand"], False,
                                         _uid(pe["protein"]), code, pe["eid"]))
            out_genes.sort(key=lambda x: (x.start, x.end, x.copy_label))
            by_label = {x.copy_label: x for x in out_genes}
            for pe in planned:
                if pe["rep"] == rep:
                    callerr.append(_callerr_row(pe, gid, rp, new, lift, by_label, ce, code))
            if code != 11:
                for x in out_genes:
                    callerr += _alt_code_rows(x, new, rp.record, ce)
            genes8 += out_genes
            seqs8[rep] = new
            gc = q6((new.count("G") + new.count("C")) / N8)
            rrow = RepliconRow(gid, rep, rp.record, rp.kind, rp.topology, N8, rp.N, len(ed), gc,
                               ids.contig_uid(new.encode("ascii")).hex())
            rep_rows.append(rrow)
            wr["replicons.tsv"].row(rrow)
            counts["edits"] += len(ed)
        for x in genes8:
            wr["genes.tsv"].row(x)
        counts["genes"] += len(genes8)

        # --- S9: classes, breaks, contigs
        grow, contig_rows, piece_rows, rep_loci = _assemble(gid, code, class_of[gid], reps, seqs8, genes8, fam_kind,
                                                            asm, seed, tmp)
        status = _status(piece_rows)
        for x in sorted(callerr, key=lambda c: c.callerr_id):
            wr["callerr.tsv"].row(x._replace(status=status.get(x.true_gene_id, "dropped")))
            counts["callerr_" + x.kind] += 1
        for x in contig_rows:
            wr["contigs.tsv"].row(x)
        for x in piece_rows:
            wr["gene_contigs.tsv"].row(x)
        for x in rep_loci:
            wr["repeats.tsv"].row(x)
        wr["genomes.tsv"].row(grow)
        genome_rows.append(grow)
        counts["genomes"] += 1
        counts["contigs"] += len(contig_rows)
        counts["contigs_emitted"] += sum(c.emitted for c in contig_rows)
        counts["genomes_" + grow.assembly_class] += 1
    for x in wr.values():
        x.close()

    files = {str(f.relative_to(tmp)): sequences._sha(f) for f in sorted(tmp.rglob("*")) if f.is_file()}
    digest = blake3("".join(f"{k}\t{v}\n" for k, v in sorted(files.items())).encode()).hexdigest(length=16)
    cfg = {"callerr": ce, "assembly": asm}
    doc = {"format": FORMAT, "seed": seed, "preset": cfg, "params_digest": backbone.params_digest(cfg),
           **{f"{d}_digest": docs[d]["digest"] for d in ("backbone", "sequences", "replicons", "mobile", "fusions")},
           "versions": {"numpy": np.__version__},
           "codes": {"code4_clade_node": node4, "code4_genomes": sorted(g for g in all_gids if code_of[g] == 4),
                     "code25_genome": gid25},
           "stream_paths": ["S8/codes", "S8/place", "S8/error/<id>", "S8/recode/<genome>", "S9/classes",
                            "S9/genome/<genome>"],
           "counts": dict(sorted(counts.items())), "files": files, "digest": digest}
    (tmp / "assembly.json").write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    tmp.rename(out)
    return read(world)


def _uid(prot: str) -> str:
    return blake3(prot.encode("ascii")).hexdigest(length=16)


def _check_gene(new: str, old: str, g: dict, s2: int, e2: int, code: int, m: dict, r: dict) -> None:
    """An edited gene reads as designed (the build refuses otherwise)."""
    st = g["st"]
    cds_old = cut(old, g["s"], g["e"], st)
    cds_new = cut(new, s2, e2, st)
    want = orf_protein(cds_old)
    if want is None or _uid(want) != r["protein_uid"]:
        raise AssertionError(f"{g['label']}: the S7 gene does not read as its protein")
    if "callerr" in m and m.get("len_delta"):
        if len(cds_new) != len(cds_old) + m["len_delta"]:
            raise AssertionError(f"{g['label']}: the frameshift does not change the length by one")
        return
    if m.get("kind") == "pseudogene":
        diff = [i for i, (x, y) in enumerate(zip(cds_old, cds_new)) if x != y]
        if len(cds_new) != len(cds_old) or len(diff) != 1 or translate(cds_new, 11)[diff[0] // 3] != "*":
            raise AssertionError(f"{g['label']}: the internal stop is not one base making a stop")
        return
    if orf_protein(cds_new, code) != want:
        raise AssertionError(f"{g['label']}: an edited gene no longer encodes its protein under code {code}")


def _callerr_row(pe: dict, gid: str, rp: _Rep, new: str, lift: Lift, genes: dict, cfg: dict,
                 code: int) -> CallerrRow:
    N8 = len(new)
    if pe["kind"] == "smorf":
        g = genes[f"{gid}_{pe['eid']}"]
        err = g.start
        alt = -1
    else:
        g = genes[pe["gene"]["label"]]
        if pe["kind"] == "frameshift_ins":
            err = lift.start(pe["err"]) - 1
        else:
            err = lift.start(pe["err"])
        alt = lift.start(pe["alt"]) if "alt" in pe else -1
    err %= N8
    if pe["kind"] in ("frameshift_ins", "frameshift_del", "internal_stop"):
        cds = cut(new, g.start, g.end, g.strand)
        calls = [_piece_coords(g.start, g.end, g.strand, a, b, N8)
                 for a, b in code11_calls(cds, int(cfg["min_call_codons"]))]
        arch = "not_split"
    else:
        calls = [(g.start, g.end)]
        arch = "single"
    uid = "" if pe["kind"] == "internal_stop" else g.protein_uid
    return CallerrRow(pe["eid"], pe["kind"], gid, rp.rep, rp.record, g.copy_label, g.start, g.end, g.strand, err,
                      alt, uid, tuple(f"{a}-{b}:{g.strand}" for a, b in calls), arch, code, "")


def _alt_code_rows(g: GeneRow, new: str, record: str, cfg: dict) -> list[CallerrRow]:
    N8 = len(new)
    cds = cut(new, g.start, g.end, g.strand)
    first = next((c for c in range(1, len(cds) // 3 - 1) if cds[3 * c:3 * c + 3] == "TGA"), None)
    if first is None:
        return []
    err = (g.start + 3 * first) if g.strand == 1 else (g.end - 3 * first - 3)
    calls = [_piece_coords(g.start, g.end, g.strand, a, b, N8) for a, b in code11_calls(cds, int(cfg["min_call_codons"]))]
    return [CallerrRow(f"AC_{g.copy_label}", "alt_code", g.genome_id, g.replicon, record, g.copy_label, g.start,
                       g.end, g.strand, err % N8, -1, g.protein_uid, tuple(f"{a}-{b}:{g.strand}" for a, b in calls),
                       "not_split", g.genetic_code, "")]


def _status(rows: list[GenePieceRow]) -> dict[str, str]:
    by = collections.defaultdict(list)
    for r in rows:
        by[r.copy_label].append(r)
    out = {}
    for lab, rs in by.items():
        em = [r for r in rs if r.emitted]
        if not em:
            out[lab] = "dropped"
        elif len(rs) == 1 and not (rs[0].partial_5 or rs[0].partial_3):
            out[lab] = "whole"
        else:
            out[lab] = "partial"
    return out


def _assemble(gid: str, code: int, klass: str, reps: dict[int, _Rep], seqs8: dict[int, str], genes8: list[GeneRow],
              fam_kind: dict, cfg: dict, seed: int, tmp: Path):
    rng = seeds.generator(seed, f"S9/genome/{gid}")
    order = sorted(seqs8)
    lens = {r: len(seqs8[r]) for r in order}
    circ = {r: reps[r].topology == "circular" for r in order}
    G = sum(lens.values())
    target_n50, target_c = 0, 1.0
    if klass == "draft":
        target_n50 = int(cfg["draft_n50"][int(rng.integers(len(cfg["draft_n50"])))])
    elif klass == "mag":
        target_n50 = int(cfg["mag_n50"])
        lo, hi = (float(x) for x in cfg["mag_completeness"])
        target_c = q6(lo + (hi - lo) * float(rng.random()))
    loci = _repeat_loci(genes8, seqs8, fam_kind, cfg)
    rep_rows = [RepeatRow(gid, g.replicon, g.copy_label, g.family_id, g.start, g.end, g.strand, why)
                for g, why in loci]
    genes_by_rep = collections.defaultdict(list)
    for g in genes8:
        genes_by_rep[g.replicon].append(g)
    cuts: dict[int, dict[int, tuple[str, str]]] = {r: {} for r in order}
    rot: dict[int, int] = {}
    n_rep_breaks = 0
    if klass == "complete":
        for r in order:
            if circ[r]:
                gs = [g for g in genes_by_rep[r] if g.end - g.start >= 3]
                if gs:
                    g = gs[int(rng.integers(len(gs)))]
                    rot[r] = (g.start + 1 + int(rng.integers(g.end - g.start - 1))) % lens[r]
                else:
                    rot[r] = int(rng.integers(lens[r]))
    else:
        n_target = max(sum(circ.values()), int(round(G * float(cfg["n50_factor"]) / target_n50)))
        use = list(range(len(loci)))
        if len(use) > n_target:
            use = sorted(int(i) for i in rng.permutation(len(loci))[:n_target])
        for i in use:
            g, _ = loci[i]
            p = ((g.start + g.end) // 2) % lens[g.replicon]
            if (circ[g.replicon] or 0 < p) and p not in cuts[g.replicon]:
                cuts[g.replicon][p] = ("break_repeat", g.family_id)
                n_rep_breaks += 1
        n_rand = max(0, n_target - n_rep_breaks)
        offs = np.cumsum([0] + [lens[r] for r in order])
        for x in rng.integers(0, G, size=n_rand):
            j = int(np.searchsorted(offs, int(x), side="right")) - 1
            r, p = order[j], int(x) - int(offs[j])
            if (circ[r] or p > 0) and p not in cuts[r]:
                cuts[r][p] = ("break_random", "")
        for r in order:
            if circ[r] and not cuts[r]:
                cuts[r][int(rng.integers(lens[r]))] = ("break_random", "")
    # pieces, frame order, replicon by replicon
    allp: list[tuple[int, Piece]] = []
    for r in order:
        for p in pieces(lens[r], circ[r], cuts[r], rot.get(r)):
            allp.append((r, p._replace(replicon=r)))
    L = [p.b - p.a for _, p in allp]
    drop = ["short" if x < int(cfg["min_contig"]) else "" for x in L]
    if klass == "mag":
        emitted = sum(x for x, d in zip(L, drop) if not d)
        for i in rng.permutation(len(allp)):
            i = int(i)
            if drop[i]:
                continue
            if (emitted - L[i]) / G >= target_c:
                drop[i] = "mag_dropout"
                emitted -= L[i]
    orient = [-1 if u < float(cfg["reverse_probability"]) else 1 for u in rng.random(len(allp))]
    # gene pieces
    per_piece = collections.defaultdict(list)
    piece_rows: list[GenePieceRow] = []
    by_rep_idx = collections.defaultdict(list)
    for k, (r, p) in enumerate(allp):
        by_rep_idx[r].append(k)
    cids = [f"{gid}_c{k:04d}" for k in range(len(allp))]
    for g in genes8:
        idx = by_rep_idx[g.replicon]
        ps = [allp[k][1] for k in idx]
        for j, x, y, cl, cr in place(ps, lens[g.replicon], circ[g.replicon], g.start, g.end):
            k = idx[j]
            p = allp[k][1]
            Lk = p.b - p.a
            cs, ce, cst = to_contig(Lk, orient[k], x, y, g.strand)
            p5, p3 = (cl, cr) if g.strand == 1 else (cr, cl)
            row = GenePieceRow(gid, g.copy_label, cids[k], cs, ce, cst, p5, p3, ce > Lk, not drop[k])
            piece_rows.append(row)
            per_piece[k].append(row)
    # contigs
    fa = sequences._Fasta(tmp / "contigs" / f"{gid}.fna")
    fd = sequences._Fasta(tmp / "dropped" / f"{gid}.fna")
    contig_rows = []
    for k, (r, p) in enumerate(allp):
        s = cut(seqs8[r], p.a, p.b, orient[k])
        left, right = (p.left, p.right) if orient[k] == 1 else (p.right, p.left)
        fl, fr = (p.fam_left, p.fam_right) if orient[k] == 1 else (p.fam_right, p.fam_left)
        (fd if drop[k] else fa).write(cids[k], s)
        npart = sum(1 for x in per_piece[k] if x.partial_5 or x.partial_3)
        contig_rows.append(ContigRow(cids[k], gid, r, reps[r].record, len(s), p.a, p.b, orient[k],
                                     p.a if p.circular else 0, p.circular, left, right, fl, fr, npart, not drop[k],
                                     drop[k], q6((s.count("G") + s.count("C")) / len(s)),
                                     ids.contig_uid(s.encode("ascii")).hex()))
    fa.close()
    fd.close()
    em = [x for x, d in zip(L, drop) if not d]
    grow = GenomeRow(gid, code, klass, target_n50, n50(em), target_c, q6(sum(em) / G), G, sum(em), len(allp),
                     len(em), sum(len(c) for c in cuts.values()), n_rep_breaks)
    return grow, contig_rows, piece_rows, rep_rows


def read(world) -> AssemblyWorld:
    """Load ``world/assembly/``'s genome and replicon tables."""
    world = Path(world)
    adir = world / "assembly"
    doc = json.loads((adir / "assembly.json").read_text(encoding="utf-8"))
    if doc["format"] != FORMAT:
        raise ValueError(f"{adir}: format {doc['format']!r}, expected {FORMAT!r}")
    t = {}
    for name in ("genomes.tsv", "replicons.tsv"):
        cls = _TABLES[name]
        conv = backbone._parser(cls)
        t[name] = [cls(*(c(r[f]) for c, f in zip(conv, cls._fields))) for r in mobile._tsv(adir / name)]
    return AssemblyWorld(doc["seed"], t["genomes.tsv"], t["replicons.tsv"], doc["digest"], world, doc)


__all__ = ["FORMAT", "ERROR_KINDS", "CLASSES", "code_table", "translate", "orf_protein", "cut", "code11_calls",
           "Lift", "Piece", "pieces", "place", "to_contig", "n50", "check_config", "build", "read", "AssemblyWorld"]
