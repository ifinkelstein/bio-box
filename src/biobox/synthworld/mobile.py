"""synthworld phages, prophages and CRISPR arrays (planning/tests/synthworld_spec.md §2.1 S6, §3.3
``truth_spacer_target``, §5.1, §5.3; card F-02 session 5). Reads a world whose S1-S5 are built (``backbone``,
``sequences``, ``replicons``) and writes ``world/mobile/``.

**Phage set.** A world of its own under ``mobile/phages/``, built by the S1-S5 code: `backbone.build_phages`
runs ZOMBI2 a second time (``phages.count`` extant phages on a tree of ``phages.birth``/``death``; one linear
chromosome per phage; the ``phages.genomes`` rates; a stem with only inversions, so the root has blocks on
both strands; no groups), then `sequences.build` and `intergenic.build` realize it. Its preset is the host's
with those values (`phage_preset`): no RNA motifs, the host's other sections as they stand. Its master seed is
derived from stream ``S6/phages/world``, so none of its streams coincide with the host's; its HMM panel is the
host's (``panel_seed``). Phage ``g<k>`` of that world is ``phage_<k>``, its families ``P<family>`` and its gene
copies ``ph_<copy_label>``; each phage is one linear record (``phages.fna``), end flanks included.

**Prophages.** ``prophages.count`` phages (stream ``S6/placement``: a permutation, the first ones, in phage
order) are integrated, ``prophage_00`` ... . Each is integrated on a host species-tree branch chosen among
those with ``clade_size`` extant genomes below them, at a time uniform along the branch, on a strand drawn with
probability 1/2, and is inherited by every extant genome below (the clade). Its site is the 3' end of an
**anchor** gene: a copy of a root or originated family that has exactly one copy in every clade genome, on the
chromosome, with no overlap at its 3' end (stream ``S6/anchor/<item>`` draws it among those families; the
placement stream ranks ``ANCHOR_TRIES`` branches per item, and the first with a free candidate is used; anchors
are distinct over the world). The prophage is inserted right after the anchor's stop codon (+ strand anchor:
at its end; - strand: at its start), so no gene, element or IGR of S5 is cut. From the integration point it
evolves down the clade (stream ``S6/prophage/<id>``): on a branch of length t in lineage b, every base outside
the phage's genes has a substitution event with probability 1 - exp(-r c_b t) (r = ``prophages.rate``, c_b the
lineage clock of S4) and takes a base from b's intergenic composition (S5); every internal codon of a gene
(start and stop codons, and codons sharing a base with another gene, excluded) has a synonymous event with the
same probability and is redrawn from b's codon usage for its amino acid. A prophage's genes therefore encode
the phage's proteins exactly.

**CRISPR-Cas loci.** ``crispr.arrays`` arrays ``CR00`` ... are seeded like prophages (branch, time, strand,
anchor; ``crispr.clade_size``). Each is one inserted locus: a Cas operon of ``cas_roles`` genes, then a leader
of ``leader_length`` bases, then the array (repeat, spacer, repeat, ..., repeat), all on the locus strand, so
the leader end of the array (5', spacer index 0) faces the operon. **As built, the Cas operon is not a ZOMBI2
group** (§2.1 S2 lists Cas groups; no earlier session built one, and adding a group to the host's ZOMBI2 run
would change every world): it is planted with its array and inherited down the clade only (no transfer).
Cas gene k of array a is family ``CAS_FAMILY_BASE + 10 a + k``, designed and evolved by the S3-S4 code on the
clade's tree (streams ``S3/family/<id>``, ``S4/family/<id>``), with S5 cassettes (``S5/family/<id>``) and the
S5 junction rules, the operon being one planted unit (tight junctions, overlaps by protein-neutral edits). Its
end flanks are ``locus_flank`` bases of the first gene's U and the last gene's D. The repeat (``repeat_length``,
its first ``repeat_stem`` bases the reverse complement of its last) and the leader (the seeding lineage's
intergenic composition) are drawn once per array from stream ``S6/crispr/<array>`` and do not change.

**Spacers** (stream ``S6/crispr/<array>``, after the repeat and leader). ``initial_spacers`` are acquired at the
seeding point; then along every branch of the clade, events in a Poisson process: an acquisition at rate
lambda = ``acquisitions_mean`` / T, inserted at the leader end (index 0), and a deletion starting at each spacer
but the leader-proximal one at rate delta = ``deletions_per_spacer`` / T (T the time from the seeding point to the
present), removing 1 + geometric(``deletion_mean_length`` - 1) spacers but never taking the array below
``min_spacers``; when an acquisition takes the array past ``max_spacers`` the trailer-end spacers go. A spacer is
acquired once and inherited unchanged. It is random (probability ``decoy_fraction``) or a slice of a phage: a
phage uniform among the set, a strand with probability 1/2, a protospacer of ``spacer_length`` uniform among the
positions followed by a PAM (``pam`` NGG read 5'->3' on the protospacer's strand, directly 3' of it; type II-A)
or, with probability ``no_pam_fraction``, among those followed by three bases that are not one, then
n_mismatches (``mismatch_weights`` over 0, 1, ...) positions changed to another base.

**Targets** (§3.3 ``truth_spacer_target``: one row per spacer x target). A spacer's targets are every window
of the phage records and of the integrated prophages (both strands) within ``target_max_mismatches`` of it,
found by an exact-block index (pigeonhole over ``target_max_mismatches`` + 1 blocks, then Hamming), the
designed protospacer among them; a spacer without one has a single row of kind ``none``. Host sequence
outside prophages is not scanned: a random 30-mer lies within 5 mismatches of a given one with probability
3.0e-11 (derived: sum over k <= 5 of C(30, k) 3^k / 4^30), about 6e-3 per spacer over a micro host (2e8
positions, derived).

**Host replicons.** Every host replicon of S5 is written again with its insertions (``fasta/<genome>.fna``,
records ``<genome>_r<replicon>`` as in S5). ``insertions.tsv`` is the lift from the S5 frame: with insertions
of lengths l_i at points a_i of a replicon of length N (a point 0 counts as N), a feature [s, e) of S5 moves to
s + sum(l_i : a_i <= s) and e + sum(l_i : a_i < e); one across the origin (e > N) ends at N' + (e - N), N' the
new length. No S5 feature straddles a point (`Lift.feature` checks). ``genes.tsv`` lists every gene of the new
replicons: the host's lifted (``source`` host) and the inserted ones (``prophage``, ``cas``).

Shape on disk (``build(world, pfam_hmm=None, rfam_cm=None)``)::

    phages/                  the phage set's world (backbone/, zombi2/, sequences/, replicons/)
    phages.fna, .fai         one record per phage, phage_<k>
    phages.tsv      phage_id genome_id species_node record length n_genes gc contig_uid prophage_id
                    n_prophage_instances
    phage_genes.tsv phage_id copy_label family_id start end strand cds_edit protein_uid
    fasta/<genome_id>.fna, .fai   host replicons with their insertions
    replicons.tsv   genome_id replicon record kind topology length s5_length n_genes n_insertions gc contig_uid
    genes.tsv       genome_id replicon copy_label family_id source insertion_id start end strand
                    origin_spanning cds_edit protein_uid
    insertions.tsv  insertion_id kind element genome_id replicon record at length start end strand anchor_family
                    anchor_copy
    prophages.tsv   prophage_id phage_id genome_id record start end strand clade_node insertion_time
                    anchor_family n_subst identity
    cas.tsv         array_id cas_module_id role order_index family family_id dark n_domains domain_hmms
                    target_identity rate_multiplier
    arrays.tsv      array_id genome_id record cas_module_id strand locus_start locus_end leader_start leader_end
                    array_start array_end n_spacers repeat clade_node insertion_time anchor_family spacer_ids
    spacers.tsv     spacer_id array_id seq acquisition_node acquisition_time source phage_id proto_start proto_end
                    proto_strand n_mismatches pam_seq pam_present
    spacer_targets.tsv array_id genome_id record array_start array_end array_strand repeat_consensus
                    cas_module_id spacer_id spacer_index spacer_start spacer_end spacer_seq acquisition_node
                    target_kind target_id target_record proto_start proto_end proto_strand n_mismatches pam_seq
                    pam_present designed
    mobile.json     format, seed, the phages/prophages/crispr presets and their digest, upstream digests, the
                    phage world's digests, versions, stream paths, counts, sha256 of every file, digest

Coordinates are 0-based half-open on the + strand of their record; strand +1/-1; a host feature across the
origin has end > length. ``prophages.identity`` is 1 - n_subst / phage length (the instance read on the phage's
strand against the phage record). A spacer's ``seq`` is 5'->3' on its array's strand, and equals its
protospacer (read 5'->3' on ``proto_strand``) but at ``n_mismatches`` positions; ``pam_seq`` is the three bases
3' of the protospacer on that strand (empty past a phage end). ``spacer_targets.target_kind`` is ``phage``,
``prophage`` or ``none``; ``target_id`` the phage or ``<prophage_id>@<genome_id>``; ``designed`` marks the
protospacer the spacer was cut from. ``acquisition_node`` is the species-tree branch of the acquisition.
Rows are in genome, replicon, start order (``genes``), in element then genome order (``insertions``, instances),
in id order (``spacers``). Floats are written at 6 significant digits. One host genome is in memory at a time;
the inserted loci (a prophage per clade genome, a locus per clade genome) are held while the host is
rewritten. ``digest`` is BLAKE3-128 over the sorted (path, sha256) list of every file but mobile.json.
"""

from __future__ import annotations

import collections
import copy
import csv
import json
import math
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

import numpy as np
from blake3 import blake3

from biobox import ids
from biobox.synthworld import backbone, codons, intergenic, profiles, seeds, sequences

FORMAT = "synthworld-mobile/1"
CAS_FAMILY_BASE = 900_000
ANCHOR_TRIES = 8
ANCHOR_KINDS = ("root", "originated")
q6 = sequences.q6
revcomp = intergenic.revcomp

_ACGT = np.frombuffer(b"ACGT", dtype=np.uint8)
_CODE = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
_TCAG = np.zeros(256, dtype=np.int64)
for _i, _b in enumerate(b"TCAG"):
    _TCAG[_b] = _i
_CMP = np.zeros(256, dtype=np.uint8)
for _a, _b in zip(b"ACGT", b"TGCA"):
    _CMP[_a] = _b


def to_codes(seq) -> np.ndarray:
    """ACGT (str, bytes or uint8 array) -> 0..3."""
    b = np.frombuffer(seq.encode("ascii"), dtype=np.uint8) if isinstance(seq, str) else np.asarray(seq, np.uint8)
    c = _CODE[b]
    if (c == 255).any():
        raise ValueError("a sequence with letters outside ACGT")
    return c


def pam_matches(pam: str, rule: str) -> bool:
    """Whether ``pam`` reads as ``rule`` (IUPAC N matches any base)."""
    return len(pam) == len(rule) and all(r == "N" or a == r for a, r in zip(pam, rule))


# --- the species tree and clades --------------------------------------------------------------------------

@dataclass
class SpeciesTree:
    labels: list[str]                     # preorder, the stem first
    parent: dict[str, str]
    children: dict[str, list[str]]
    start: dict[str, float]
    end: dict[str, float]
    below: dict[str, tuple[str, ...]]     # extant tips in the subtree, preorder
    present: float


def species_tree(branches, extant) -> SpeciesTree:
    extant = set(extant)
    labels = [b.lineage for b in branches]
    parent = {b.lineage: b.parent for b in branches}
    children: dict[str, list[str]] = {x: [] for x in labels}
    start, end = {}, {}
    for b in branches:
        if b.parent:
            children[b.parent].append(b.lineage)
        start[b.lineage] = end[b.parent] if b.parent else 0.0
        end[b.lineage] = start[b.lineage] + float(b.branch_length)
    below: dict[str, tuple[str, ...]] = {}
    for x in reversed(labels):
        below[x] = ((x,) if x in extant else ()) + tuple(t for c in children[x] for t in below[c])
    present = max(end[x] for x in extant)
    return SpeciesTree(labels, parent, children, start, end, below, present)


class CladeNode(NamedTuple):
    lineage: str
    parent: int            # index in the clade, -1 for the root
    branch_length: float   # the root's: from the insertion time to the end of its branch
    start: float           # absolute time at the top of the (part) branch
    tip: bool


@dataclass
class Clade:
    node: str
    t_insert: float
    nodes: list[CladeNode]  # preorder over the branches with an extant descendant

    @property
    def tips(self) -> list[str]:
        return [n.lineage for n in self.nodes if n.tip]


def clade(tree: SpeciesTree, node: str, t_insert: float) -> Clade:
    if not tree.start[node] <= t_insert <= tree.end[node]:
        raise ValueError(f"time {t_insert} is not on branch {node}")
    out: list[CladeNode] = []

    def walk(x: str, par: int) -> None:
        if not tree.below[x]:
            return
        root = par == -1
        bl = tree.end[x] - (t_insert if root else tree.start[x])
        out.append(CladeNode(x, par, bl, t_insert if root else tree.start[x], not tree.children[x]))
        me = len(out) - 1
        for c in tree.children[x]:
            walk(c, me)

    walk(node, -1)
    return Clade(node, float(t_insert), out)


def gene_tree(cl: Clade, family: int, suffix: str) -> list[backbone.GeneTreeNode]:
    """The clade as a gene tree of one family (`backbone.GeneTreeNode` rows), leaves ``<lineage>_<suffix>``."""
    return [backbone.GeneTreeNode(int(family), i, n.parent, "extant" if n.tip else "speciation", n.lineage,
                                  f"{n.lineage}_{suffix}" if n.tip else "", n.branch_length, n.tip)
            for i, n in enumerate(cl.nodes)]


def _down(cl: Clade, state0, branch, copy_state):
    """Evolve ``state0`` from the clade's root down its nodes in preorder: ``branch(state, node)`` in place;
    returns the tips' states by lineage."""
    n_children = [0] * len(cl.nodes)
    for n in cl.nodes:
        if n.parent >= 0:
            n_children[n.parent] += 1
    states: dict[int, object] = {}
    tips = {}
    for i, n in enumerate(cl.nodes):
        if n.parent == -1:
            st = copy_state(state0)
        else:
            src = states[n.parent]
            n_children[n.parent] -= 1
            st = copy_state(src) if n_children[n.parent] > 0 else src
            if n_children[n.parent] == 0:
                del states[n.parent]
        branch(st, n)
        if n.tip:
            tips[n.lineage] = st
        else:
            states[i] = st
    return tips


# --- the lift from the S5 frame -----------------------------------------------------------------------------

class Lift:
    """Insertions (point a in [0, N], length l; a point 0 counts as N) into a record of length N."""

    def __init__(self, N: int, points):
        self.N = int(N)
        pts = sorted((self.N if int(a) == 0 else int(a), int(l)) for a, l in points)
        self.a = np.array([a for a, _ in pts], dtype=np.int64)
        self.l = np.array([l for _, l in pts], dtype=np.int64)
        if len(self.a) and (self.a.min() < 0 or self.a.max() > self.N):
            raise ValueError(f"insertion points must lie in [0, {self.N}]")
        self.cum = np.concatenate([[0], np.cumsum(self.l)])
        self.new_length = self.N + int(self.l.sum())

    def start(self, x: int) -> int:
        return int(x) + int(self.cum[np.searchsorted(self.a, x, side="right")])

    def end(self, x: int) -> int:
        return int(x) + int(self.cum[np.searchsorted(self.a, x, side="left")])

    def feature(self, s: int, e: int) -> tuple[int, int]:
        """[s, e) of the old record in the new one (e > N across the origin); refuses a feature that
        straddles an insertion point."""
        s, e = int(s), int(e)
        for a in self.a:
            if s < a < e or (e > self.N and (s < a + self.N < e or a == self.N)):
                raise ValueError(f"feature [{s}, {e}) straddles an insertion at {a}")
        if e <= self.N:
            return self.start(s), self.end(e)
        return self.start(s), self.new_length + (e - self.N)

    def blocks(self) -> list[int]:
        """Where each insertion starts in the new record, in point order."""
        return [int(a) + int(self.cum[k]) for k, a in enumerate(self.a)]


def insert_blocks(seq: str, points) -> tuple[str, list[int]]:
    """``seq`` with blocks inserted (``points``: (a, block) with a in [0, N], 0 counting as N, sorted by a and
    then in the given order); returns the new sequence and each block's start, in the order given."""
    N = len(seq)
    order = sorted(range(len(points)), key=lambda k: (N if points[k][0] == 0 else points[k][0], k))
    parts, starts, prev, shift = [], [0] * len(points), 0, 0
    for k in order:
        a = N if points[k][0] == 0 else int(points[k][0])
        parts.append(seq[prev:a])
        starts[k] = a + shift
        parts.append(points[k][1])
        shift += len(points[k][1])
        prev = a
    parts.append(seq[prev:])
    return "".join(parts), starts


def flip(span: tuple[int, int], strand: int, length: int, block_strand: int) -> tuple[int, int, int]:
    """A feature of a block laid out on + (span, strand) once the block is inserted on ``block_strand``."""
    s, e = span
    return (s, e, strand) if block_strand == 1 else (length - e, length - s, -strand)


# --- prophage evolution -------------------------------------------------------------------------------------

def coding_sites(L: int, genes) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """For a linear record of length L with ``genes`` (start, end, strand): the bases outside every gene, and
    each gene's internal codons (start and stop codons left out, and any codon sharing a base with another
    gene) as their leftmost record position and strand."""
    cover = np.zeros(L, dtype=np.int32)
    for s, e, _ in genes:
        cover[s:e] += 1
    pos, strand = [], []
    for s, e, st in genes:
        n = (e - s) // 3
        k = np.arange(1, n - 1, dtype=np.int64)
        p = s + 3 * k if st == 1 else e - 3 * (k + 1)
        ok = (cover[p] + cover[p + 1] + cover[p + 2]) == 3
        pos.append(p[ok])
        strand.append(np.full(int(ok.sum()), st, dtype=np.int8))
    return (np.flatnonzero(cover == 0), np.concatenate(pos) if pos else np.zeros(0, np.int64),
            np.concatenate(strand) if strand else np.zeros(0, np.int8))


def read_codons(seq: np.ndarray, pos: np.ndarray, strand: np.ndarray) -> np.ndarray:
    """Codon indices (`codons.CODONS` order) of the codons at ``pos`` read on ``strand`` (ASCII array)."""
    b0, b1, b2 = seq[pos], seq[pos + 1], seq[pos + 2]
    minus = strand == -1
    c0 = np.where(minus, _CMP[b2], b0)
    c1 = np.where(minus, _CMP[b1], b1)
    c2 = np.where(minus, _CMP[b0], b2)
    return 16 * _TCAG[c0] + 4 * _TCAG[c1] + _TCAG[c2]


def write_codons(seq: np.ndarray, pos: np.ndarray, strand: np.ndarray, idx: np.ndarray) -> None:
    c = codons.CODON_BYTES[idx]
    minus = strand == -1
    seq[pos] = np.where(minus, _CMP[c[:, 2]], c[:, 0])
    seq[pos + 1] = np.where(minus, _CMP[c[:, 1]], c[:, 1])
    seq[pos + 2] = np.where(minus, _CMP[c[:, 0]], c[:, 2])


def evolve_prophage(seq0: str, genes, cl: Clade, clock, usage_cum, comp_cum, index: dict[str, int], rate: float,
                    rng: np.random.Generator) -> dict[str, str]:
    """The prophage at each clade tip (module docstring): ``clock``, ``usage_cum`` (n, 22, 64) and
    ``comp_cum`` (n, 4) per lineage index (``index``)."""
    igr, pos, strand = coding_sites(len(seq0), genes)
    base = np.frombuffer(seq0.encode("ascii"), dtype=np.uint8).copy()

    def branch(st: np.ndarray, n: CladeNode) -> None:
        li = index[n.lineage]
        p = -math.expm1(-float(rate) * float(clock[li]) * float(n.branch_length))
        if p <= 0.0:
            return
        hit = igr[rng.random(len(igr)) < p]
        if hit.size:
            st[hit] = _ACGT[np.minimum(np.count_nonzero(comp_cum[li][None, :] <= rng.random(hit.size)[:, None],
                                                        axis=1), 3)]
        ch = np.flatnonzero(rng.random(len(pos)) < p)
        if ch.size:
            old = read_codons(st, pos[ch], strand[ch])
            res = codons.CODE11_RESIDUE[old].astype(np.int64)
            if (res >= 20).any():
                raise AssertionError("an internal codon of a prophage gene is a stop codon")
            write_codons(st, pos[ch], strand[ch], codons.draw(usage_cum[li], res, rng.random(ch.size)))

    tips = _down(cl, base, branch, lambda s: s.copy())
    return {k: v.tobytes().decode("ascii") for k, v in tips.items()}


# --- spacers ------------------------------------------------------------------------------------------------

@dataclass
class Spacer:
    spacer_id: str
    array_id: str
    seq: str
    node: str
    time: float
    source: str            # phage | decoy
    phage: int = -1
    start: int = -1
    end: int = -1
    strand: int = 0
    n_mismatches: int = 0
    pam_seq: str = ""
    pam_present: bool = False


class PhageSites:
    """Per phage record, the protospacer starts on each strand followed (or not) by a PAM."""

    def __init__(self, seqs: list[str], length: int, pam: str):
        self.seqs = seqs
        self.length = int(length)
        self.pam = pam
        self.sites = []
        k = len(pam)
        for s in seqs:
            L, n = len(s), self.length
            arr = np.frombuffer(s.encode("ascii"), dtype=np.uint8)
            plus = np.arange(0, max(0, L - n - k + 1))
            minus = np.arange(k, max(k, L - n + 1))
            ok_p = np.ones(len(plus), dtype=bool)
            ok_m = np.ones(len(minus), dtype=bool)
            for j, r in enumerate(pam):
                if r != "N":
                    ok_p &= arr[plus + n + j] == ord(r)
                    ok_m &= _CMP[arr[minus - 1 - j]] == ord(r)     # the PAM on - is revcomp(s[w - k:w])
            self.sites.append({(1, True): plus[ok_p], (1, False): plus[~ok_p],
                               (-1, True): minus[ok_m], (-1, False): minus[~ok_m]})

    def protospacer(self, j: int, w: int, strand: int) -> tuple[str, str]:
        """The protospacer at [w, w + length) of phage j read on ``strand``, and its PAM bases."""
        s, n, k = self.seqs[j], self.length, len(self.pam)
        if strand == 1:
            return s[w:w + n], s[w + n:w + n + k]
        return revcomp(s[w:w + n]), revcomp(s[max(0, w - k):w])


def mutate(seq: str, m: int, rng: np.random.Generator) -> str:
    """``seq`` with ``m`` distinct positions changed to another base."""
    if m == 0:
        return seq
    b = bytearray(seq.encode("ascii"))
    for p in sorted(int(x) for x in rng.choice(len(b), size=m, replace=False)):
        b[p] = _ACGT[(int(_CODE[b[p]]) + int(rng.integers(1, 4))) % 4]
    return b.decode("ascii")


def draw_spacer(cfg: dict, sites: PhageSites, rng: np.random.Generator, sid: str, aid: str, node: str,
                time: float) -> Spacer:
    n = int(cfg["spacer_length"])
    if rng.random() < float(cfg["decoy_fraction"]):
        return Spacer(sid, aid, _ACGT[rng.integers(0, 4, size=n)].tobytes().decode("ascii"), node, time, "decoy")
    j = int(rng.integers(len(sites.seqs)))
    strand = 1 if rng.random() < 0.5 else -1
    pam = not (rng.random() < float(cfg["no_pam_fraction"]))
    cand = sites.sites[j][(strand, pam)]
    if not len(cand):
        raise ValueError(f"phage {j} has no protospacer site of kind {(strand, pam)}")
    w = int(cand[int(rng.integers(len(cand)))])
    proto, pam_seq = sites.protospacer(j, w, strand)
    wts = np.cumsum(np.asarray(cfg["mismatch_weights"], dtype=np.float64))
    m = min(int(np.searchsorted(wts / wts[-1], rng.random(), side="right")), len(wts) - 1)
    return Spacer(sid, aid, mutate(proto, m, rng), node, time, "phage", j, w, w + n, strand, m, pam_seq,
                  pam_matches(pam_seq, sites.pam))


def evolve_array(cl: Clade, cfg: dict, acquire, rng: np.random.Generator) -> dict[str, list[int]]:
    """Spacer ids per clade tip, index 0 leader-proximal (module docstring). ``acquire(node, time)`` registers
    a new spacer and returns its id (ids increase in acquisition order)."""
    T = max(1e-9, max(n.start + n.branch_length for n in cl.nodes) - cl.t_insert)
    lam = float(cfg["acquisitions_mean"]) / T
    delta = float(cfg["deletions_per_spacer"]) / T
    nmin, nmax = int(cfg["min_spacers"]), int(cfg["max_spacers"])
    dm = float(cfg["deletion_mean_length"])
    first: list[int] = []
    for _ in range(int(cfg["initial_spacers"])):
        first.insert(0, acquire(cl.node, cl.t_insert))
    del first[nmax:]

    def branch(sp: list[int], n: CladeNode) -> None:
        t = 0.0
        while True:
            k = len(sp) - 1 if len(sp) > nmin else 0
            rate = lam + delta * k
            t += float(rng.exponential(1.0 / rate))
            if t >= n.branch_length:
                return
            if rng.random() * rate < lam:
                sp.insert(0, acquire(n.lineage, n.start + t))
                del sp[nmax:]
            else:
                at = 1 + int(rng.integers(k))
                ell = 1 + (int(rng.geometric(1.0 / dm)) - 1 if dm > 1.0 else 0)
                ell = min(ell, len(sp) - at, len(sp) - nmin)
                del sp[at:at + ell]

    return _down(cl, first, branch, list)


# --- the target scan ----------------------------------------------------------------------------------------

class KmerIndex:
    """Exact k-mer positions of one sequence (codes 0..3)."""

    def __init__(self, codes: np.ndarray, k: int):
        self.codes = np.asarray(codes, dtype=np.uint8)
        self.k = int(k)
        n = len(self.codes) - self.k + 1
        km = np.zeros(max(0, n), dtype=np.int64)
        for i in range(self.k):
            km = km * 4 + self.codes[i:i + max(0, n)]
        self.order = np.argsort(km, kind="stable")
        self.sorted = km[self.order]

    def positions(self, code: int) -> np.ndarray:
        lo, hi = np.searchsorted(self.sorted, code, "left"), np.searchsorted(self.sorted, code, "right")
        return self.order[lo:hi]


def block_length(n: int, max_mm: int) -> int:
    """The length of the max_mm + 1 disjoint blocks of an n-mer (pigeonhole: one is exact)."""
    return max(1, n // (int(max_mm) + 1))


def scan(index: KmerIndex, Q: np.ndarray, max_mm: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Every window of ``index``'s sequence within ``max_mm`` mismatches of a query (rows of ``Q``, codes, one
    length n), on the sequence's own strand: (query row, window start, mismatches). Candidates are the windows
    holding one of the max_mm + 1 disjoint blocks of a query exactly, so none is missed (pigeonhole) when
    ``index.k`` is at most n // (max_mm + 1)."""
    Q = np.atleast_2d(np.asarray(Q, dtype=np.uint8))
    m, n = Q.shape
    k, L = index.k, len(index.codes)
    empty = np.zeros(0, np.int64)
    if L < n or m == 0:
        return empty, empty, empty
    if k * (int(max_mm) + 1) > n:
        raise ValueError(f"blocks of {k} cannot hold {max_mm} + 1 disjoint exact blocks in {n}")
    pw = 4 ** np.arange(k - 1, -1, -1, dtype=np.int64)
    qs, ws = [], []
    for j in range(int(max_mm) + 1):
        off = j * k
        code = Q[:, off:off + k].astype(np.int64) @ pw
        lo = np.searchsorted(index.sorted, code, "left")
        cnt = np.searchsorted(index.sorted, code, "right") - lo
        tot = int(cnt.sum())
        if not tot:
            continue
        qi = np.repeat(np.arange(m), cnt)
        at = np.repeat(lo, cnt) + (np.arange(tot) - np.repeat(np.cumsum(cnt) - cnt, cnt))
        w = index.order[at] - off
        ok = (w >= 0) & (w <= L - n)
        qs.append(qi[ok])
        ws.append(w[ok])
    if not qs:
        return empty, empty, empty
    key = np.unique(np.concatenate(qs).astype(np.int64) * (L + 1) + np.concatenate(ws))
    qi, w = key // (L + 1), key % (L + 1)
    mm = (index.codes[w[:, None] + np.arange(n)] != Q[qi]).sum(axis=1)
    keep = mm <= int(max_mm)
    return qi[keep], w[keep], mm[keep]


def pam_at(seq: str, start: int, n: int, strand: int, k: int, circular: bool) -> str:
    """The k bases 3' of the window [start, start + n) read on ``strand`` (empty past a linear end)."""
    N = len(seq)
    a, b = (start + n, start + n + k) if strand == 1 else (start - k, start)
    if circular:
        p = "".join(seq[i % N] for i in range(a, b))
    elif a < 0 or b > N:
        return ""
    else:
        p = seq[a:b]
    return p if strand == 1 else revcomp(p)


def find_targets(seq: str, queries: list[str], max_mm: int, pam: str, *, circular: bool = False, lo: int = 0,
                 hi: int | None = None) -> list[tuple[int, int, int, int, int, str]]:
    """The windows of ``seq[lo:hi]`` within ``max_mm`` of each query on both strands, as (query index, start,
    end, strand, mismatches, PAM bases) in ``seq`` coordinates; the PAM is read from ``seq`` (wrapping when
    ``circular``). A window on -1 reads as the query on the minus strand."""
    if not queries:
        return []
    hi = len(seq) if hi is None else hi
    n = len(queries[0])
    if any(len(q) != n for q in queries):
        raise ValueError("queries of one length only")
    ix = KmerIndex(to_codes(seq[lo:hi]), block_length(n, max_mm))
    out = []
    for strand, qs in ((1, queries), (-1, [revcomp(q) for q in queries])):
        qi, w, mm = scan(ix, np.stack([to_codes(q) for q in qs]), max_mm)
        for i, s, x in zip(qi.tolist(), w.tolist(), mm.tolist()):
            out.append((i, lo + s, lo + s + n, strand, x, pam_at(seq, lo + s, n, strand, len(pam), circular)))
    return out


class Target(NamedTuple):
    kind: str              # phage | prophage
    target_id: str
    record: str
    start: int
    end: int
    strand: int
    n_mismatches: int
    pam_seq: str
    pam_present: bool


# --- the phage set ------------------------------------------------------------------------------------------

def phage_preset(preset: dict) -> dict:
    """The phage set's preset: the host's, with the ``phages`` values, no groups or declared families, no burst,
    no lineage axes and no RNA motifs."""
    ph = preset["phages"]
    p = {k: copy.deepcopy(v) for k, v in preset.items() if k not in ("phages", "prophages", "crispr")}
    p["species"] = dict(p["species"], n_extant=int(ph["count"]), birth=float(ph["birth"]),
                        death=float(ph["death"]), burst_size=1)
    p["genomes"] = dict(p["genomes"], **copy.deepcopy(ph["genomes"]))
    p["root_scramble"] = copy.deepcopy(ph["root_scramble"])
    p["modules"] = dict(p["modules"], count=0)
    p["profile_only"] = dict(p["profile_only"], count=0)
    p["plasmids"] = dict(p["plasmids"], sizes=[])
    p["vehicles"] = dict(p["vehicles"], cargo=[])
    p["hubs"] = dict(p["hubs"], count=0)
    p["decoys"] = dict(p["decoys"], partner_modules=[])
    p["lineage_rates"] = dict(p["lineage_rates"], axes=[])
    p["rna_motifs"] = dict(p["rna_motifs"], source="synthetic", n_models=0, rfam_models=[])
    return p


def phage_seed(seed: int) -> int:
    w = seeds.seed_sequence(seed, "S6/phages/world").generate_state(2, dtype=np.uint32)
    return (int(w[0]) << 32) | int(w[1])


# --- I/O ----------------------------------------------------------------------------------------------------

class PhageRow(NamedTuple):
    phage_id: str
    genome_id: str
    species_node: str
    record: str
    length: int
    n_genes: int
    gc: float
    contig_uid: str
    prophage_id: str
    n_prophage_instances: int


class PhageGeneRow(NamedTuple):
    phage_id: str
    copy_label: str
    family_id: str
    start: int
    end: int
    strand: int
    cds_edit: str
    protein_uid: str


class RepliconRow(NamedTuple):
    genome_id: str
    replicon: int
    record: str
    kind: str
    topology: str
    length: int
    s5_length: int
    n_genes: int
    n_insertions: int
    gc: float
    contig_uid: str


class GeneRow(NamedTuple):
    genome_id: str
    replicon: int
    copy_label: str
    family_id: str
    source: str
    insertion_id: str
    start: int
    end: int
    strand: int
    origin_spanning: bool
    cds_edit: str
    protein_uid: str


class InsertionRow(NamedTuple):
    insertion_id: str
    kind: str
    element: str
    genome_id: str
    replicon: int
    record: str
    at: int
    length: int
    start: int
    end: int
    strand: int
    anchor_family: str
    anchor_copy: str


class ProphageRow(NamedTuple):
    prophage_id: str
    phage_id: str
    genome_id: str
    record: str
    start: int
    end: int
    strand: int
    clade_node: str
    insertion_time: float
    anchor_family: str
    n_subst: int
    identity: float


class CasRow(NamedTuple):
    array_id: str
    cas_module_id: str
    role: str
    order_index: int
    family: int
    family_id: str
    dark: bool
    n_domains: int
    domain_hmms: tuple[str, ...]
    target_identity: float
    rate_multiplier: float


class ArrayRow(NamedTuple):
    array_id: str
    genome_id: str
    record: str
    cas_module_id: str
    strand: int
    locus_start: int
    locus_end: int
    leader_start: int
    leader_end: int
    array_start: int
    array_end: int
    n_spacers: int
    repeat: str
    clade_node: str
    insertion_time: float
    anchor_family: str
    spacer_ids: tuple[str, ...]


class SpacerRow(NamedTuple):
    spacer_id: str
    array_id: str
    seq: str
    acquisition_node: str
    acquisition_time: float
    source: str
    phage_id: str
    proto_start: int
    proto_end: int
    proto_strand: int
    n_mismatches: int
    pam_seq: str
    pam_present: bool


class SpacerTargetRow(NamedTuple):
    array_id: str
    genome_id: str
    record: str
    array_start: int
    array_end: int
    array_strand: int
    repeat_consensus: str
    cas_module_id: str
    spacer_id: str
    spacer_index: int
    spacer_start: int
    spacer_end: int
    spacer_seq: str
    acquisition_node: str
    target_kind: str
    target_id: str
    target_record: str
    proto_start: int
    proto_end: int
    proto_strand: int
    n_mismatches: int
    pam_seq: str
    pam_present: bool
    designed: bool


_TABLES = {"phages.tsv": PhageRow, "phage_genes.tsv": PhageGeneRow, "replicons.tsv": RepliconRow,
           "genes.tsv": GeneRow, "insertions.tsv": InsertionRow, "prophages.tsv": ProphageRow, "cas.tsv": CasRow,
           "arrays.tsv": ArrayRow, "spacers.tsv": SpacerRow, "spacer_targets.tsv": SpacerTargetRow}


class _Rows:
    def __init__(self, path: Path, cls):
        self.fh = open(path, "w", encoding="utf-8", newline="")
        self.fh.write("\t".join(cls._fields) + "\n")

    def row(self, r) -> None:
        self.fh.write("\t".join(backbone._cell(v) for v in r) + "\n")

    def close(self) -> None:
        self.fh.close()


def check_config(ph: dict, pro: dict, cr: dict) -> None:
    if int(pro["count"]) > int(ph["count"]):
        raise ValueError(f"{pro['count']} prophages need as many phages, the set has {ph['count']}")
    for name, cfg in (("prophages", pro), ("crispr", cr)):
        lo, hi = (int(x) for x in cfg["clade_size"])
        if not 1 <= lo <= hi:
            raise ValueError(f"{name}.clade_size must be 1 <= lo <= hi, got {cfg['clade_size']}")
    if int(cr["arrays"]) > 100 or len(cr["cas_roles"]) > 10 or not cr["cas_roles"]:
        raise ValueError("at most 100 arrays and 1 to 10 Cas roles (family numbers 900000 + 10 a + k)")
    if len(cr["mismatch_weights"]) > int(cr["target_max_mismatches"]) + 1:
        raise ValueError("a spacer's designed mismatches cannot exceed target_max_mismatches")
    if not 0 <= int(cr["min_spacers"]) <= int(cr["initial_spacers"]) <= int(cr["max_spacers"]):
        raise ValueError("need min_spacers <= initial_spacers <= max_spacers")
    if not 2 * int(cr["repeat_stem"]) <= int(cr["repeat_length"]):
        raise ValueError("repeat_stem is at most half the repeat")
    if set(str(cr["pam"])) - set("ACGTN"):
        raise ValueError(f"pam is IUPAC A, C, G, T or N, got {cr['pam']!r}")


def _tsv(path: Path) -> list[dict[str, str]]:
    return intergenic._tsv(path)


def _by_genome(path: Path):
    """A TSV sorted by genome_id, one genome's rows at a time."""
    return intergenic._genome_rows(path)


@dataclass
class Item:
    """One insertion in a clade: a prophage or a CRISPR-Cas locus."""
    item_id: str
    kind: str              # prophage | crispr
    tries: list[tuple[str, float]]   # ranked (branch, fraction of the branch)
    strand: int
    phage: int = -1
    node: str = ""
    t_insert: float = 0.0
    clade: Clade | None = None
    anchor: int = -1
    tips: dict = field(default_factory=dict)       # lineage -> per-genome content


@dataclass
class MobileWorld:
    seed: int
    phages: list[PhageRow]
    replicons: list[RepliconRow]
    insertions: list[InsertionRow]
    prophages: list[ProphageRow]
    cas: list[CasRow]
    arrays: list[ArrayRow]
    spacers: list[SpacerRow]
    digest: str
    out: Path
    meta: dict = field(default_factory=dict, repr=False)


def _phage_world(tmp: Path, preset: dict, seed: int, pfam_hmm, rfam_cm) -> Path:
    pw = tmp / "phages"
    backbone.build_phages(phage_preset(preset), phage_seed(seed), pw, panel_seed=seed)
    sequences.build(pw, pfam_hmm=pfam_hmm)
    intergenic.build(pw, rfam_cm=rfam_cm)
    return pw


def build(world, *, pfam_hmm=None, rfam_cm=None) -> MobileWorld:
    """Write ``world/mobile/`` (module docstring) from a world whose S1-S5 are built."""
    world = Path(world)
    bb = json.loads((world / "backbone" / "backbone.json").read_text(encoding="utf-8"))
    sq = json.loads((world / "sequences" / "sequences.json").read_text(encoding="utf-8"))
    rp = json.loads((world / "replicons" / "replicons.json").read_text(encoding="utf-8"))
    for doc, want in ((bb, backbone.FORMAT), (sq, sequences.FORMAT), (rp, intergenic.FORMAT)):
        if doc["format"] != want:
            raise ValueError(f"{world}: format {doc['format']!r}, expected {want!r}")
    preset = bb["preset"]
    ph, pro, cr = preset["phages"], preset["prophages"], preset["crispr"]
    check_config(ph, pro, cr)
    seed = int(bb["seed"])
    out = world / "mobile"
    if out.exists():
        raise FileExistsError(f"{out} exists; the mobile stage is written once")
    tmp = world / "mobile.partial"
    if tmp.exists():
        shutil.rmtree(tmp)
    (tmp / "fasta").mkdir(parents=True)

    # --- the phage set
    pw = _phage_world(tmp, preset, seed, pfam_hmm, rfam_cm)
    pgen = {r["species_node"]: r["genome_id"] for r in _tsv(pw / "backbone" / "genomes.tsv")}
    phage_ids = {gid: f"phage_{int(gid[1:]):02d}" for gid in sorted(pgen.values())}
    precs = {}
    for f in sorted((pw / "replicons" / "fasta").glob("*.fna")):
        for name, s in intergenic._fasta_records(f):
            precs[name.split("_r")[0]] = s
    pgids = sorted(phage_ids)
    pseqs = [precs[g] for g in pgids]
    puid = {r["copy_label"]: r["protein_uid"] for r in _tsv(pw / "sequences" / "genes.tsv")}
    pgenes: dict[str, list[dict]] = collections.defaultdict(list)
    for r in _tsv(pw / "replicons" / "genes.tsv"):
        pgenes[r["genome_id"]].append(r)

    # --- host model inputs
    branches = backbone.species_branches(world)
    host_genomes = _tsv(world / "backbone" / "genomes.tsv")
    gid_of = {r["species_node"]: r["genome_id"] for r in host_genomes}
    tree = species_tree(branches, gid_of)
    seq_cfg = sq["sequences_preset"]
    igr_cfg = dict(preset["intergenic"], _motif_rate_factor=float(preset["rna_motifs"]["rate_factor"]))
    lin_seq = sequences.lineage_model(branches, seq_cfg, seed)
    lin_igr = intergenic.lineages(world, preset["intergenic"]["composition"])
    kappa = float(seq_cfg["evolution"]["kappa_syn"])
    kinds = {int(r["family"]): r["kind"] for r in _tsv(world / "backbone" / "families.tsv")}

    # --- placement (stream S6/placement): prophages, then arrays
    rng = seeds.generator(seed, "S6/placement")
    items: list[Item] = []
    n_pro = int(pro["count"])
    chosen = sorted(int(x) for x in rng.permutation(len(pgids))[:n_pro])
    for spec, n_items in ((pro, n_pro), (cr, int(cr["arrays"]))):
        lo, hi = (int(x) for x in spec["clade_size"])
        cand = [x for x in tree.labels[1:] if lo <= len(tree.below[x]) <= hi]
        if not cand:
            raise ValueError(f"no species branch has {lo} to {hi} extant genomes below it")
        for i in range(n_items):
            perm = rng.permutation(len(cand))[:ANCHOR_TRIES]
            tries = [(cand[int(k)], float(rng.random())) for k in perm]
            strand = 1 if rng.random() < 0.5 else -1
            if spec is pro:
                items.append(Item(f"prophage_{i:02d}", "prophage", tries, strand, phage=chosen[i]))
            else:
                items.append(Item(f"CR{i:02d}", "crispr", tries, strand))

    # --- pass 1: per (item, try), the families with one copy per clade genome that can anchor
    gen_tries: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    for ii, it in enumerate(items):
        for ti, (node, _) in enumerate(it.tries):
            for tip in tree.below[node]:
                gen_tries[gid_of[tip]].append((ii, ti))
    stats: dict[tuple[int, int], collections.Counter] = collections.defaultdict(collections.Counter)
    juncs = _by_genome(world / "replicons" / "junctions.tsv")
    for gid, rows in _by_genome(world / "replicons" / "genes.tsv"):
        jg, jrows = next(juncs)
        if jg != gid:
            raise ValueError(f"genes.tsv and junctions.tsv disagree on genome order at {gid}")
        if gid not in gen_tries:
            continue
        ok = _anchor_families(rows, jrows, kinds)
        for key in gen_tries[gid]:
            stats[key].update(ok)
    used: set[int] = set()
    for ii, it in enumerate(items):
        for ti, (node, frac) in enumerate(it.tries):
            k = len(tree.below[node])
            fams = sorted(f for f, c in stats[(ii, ti)].items() if c == k and f not in used)
            if fams:
                a = seeds.generator(seed, f"S6/anchor/{it.item_id}")
                it.anchor = fams[int(a.integers(len(fams)))]
                it.node = node
                it.t_insert = q6(tree.start[node] + frac * (tree.end[node] - tree.start[node]))
                it.t_insert = min(max(it.t_insert, tree.start[node]), tree.end[node])
                it.clade = clade(tree, node, it.t_insert)
                used.add(it.anchor)
                break
        if it.clade is None:
            raise ValueError(f"{it.item_id}: none of {len(it.tries)} branches has a free anchor family")

    # --- the inserted contents
    panel = profiles.panel(seq_cfg["profiles"], seed, pfam_hmm)
    cas_rows: list[CasRow] = []
    spacers: dict[str, list[Spacer]] = {}
    sites = PhageSites(pseqs, int(cr["spacer_length"]), str(cr["pam"]))
    for it in items:
        if it.kind == "prophage":
            g = pgids[it.phage]
            genes = [(int(r["start"]), int(r["end"]), int(r["strand"])) for r in pgenes[g]]
            it.tips = evolve_prophage(pseqs[it.phage], genes, it.clade, lin_seq.clock, lin_seq.cum, lin_igr.comp_cum,
                                      lin_seq.index, float(pro["rate"]),
                                      seeds.generator(seed, f"S6/prophage/{it.item_id}"))
        else:
            a = int(it.item_id[2:])
            it.tips, rows, spacers[it.item_id] = _crispr_locus(it, a, cr, seq_cfg, igr_cfg, kappa, seed, panel,
                                                               lin_seq, lin_igr, sites)
            cas_rows += rows

    # --- phage records and tables
    fa = sequences._Fasta(tmp / "phages.fna")
    t_ph = _Rows(tmp / "phages.tsv", PhageRow)
    t_pg = _Rows(tmp / "phage_genes.tsv", PhageGeneRow)
    pro_of = {it.phage: it for it in items if it.kind == "prophage"}
    phage_rows = []
    for j, g in enumerate(pgids):
        pid, s = phage_ids[g], pseqs[j]
        fa.write(pid, s)
        it = pro_of.get(j)
        node = next(n for n, x in pgen.items() if x == g)
        row = PhageRow(pid, g, node, pid, len(s), len(pgenes[g]), q6((s.count("G") + s.count("C")) / len(s)),
                       ids.contig_uid(s.encode("ascii")).hex(), it.item_id if it else "",
                       len(it.tips) if it else 0)
        t_ph.row(row)
        phage_rows.append(row)
        for r in pgenes[g]:
            t_pg.row(PhageGeneRow(pid, f"ph_{r['copy_label']}", "P" + r["family_id"][1:], int(r["start"]),
                                  int(r["end"]), int(r["strand"]), r["cds_edit"], puid[r["copy_label"]]))
    fa.close()
    t_ph.close()
    t_pg.close()

    # --- pass 2: one host genome at a time
    max_mm = int(cr["target_max_mismatches"])
    pam = str(cr["pam"])
    distinct = [s for aid in sorted(spacers) for s in spacers[aid]]
    queries = [s.seq for s in distinct]
    targets: dict[str, list[Target]] = {s.spacer_id: [] for s in distinct}
    for j, g in enumerate(pgids):
        for qi, a, b, st, mm, p in find_targets(pseqs[j], queries, max_mm, pam):
            targets[distinct[qi].spacer_id].append(Target("phage", phage_ids[g], phage_ids[g], a, b, st, mm, p,
                                                          pam_matches(p, pam)))
    reps_s5 = collections.defaultdict(list)
    for r in _tsv(world / "replicons" / "replicons.tsv"):
        reps_s5[r["genome_id"]].append(r)
    by_genome_items: dict[str, list[Item]] = collections.defaultdict(list)
    for it in items:
        for tip in it.clade.tips:
            by_genome_items[gid_of[tip]].append(it)
    writers = {n: _Rows(tmp / n, c) for n, c in _TABLES.items()
               if n in ("replicons.tsv", "genes.tsv", "insertions.tsv")}
    instances_pro: list[ProphageRow] = []
    instances_arr: list[tuple[ArrayRow, list[tuple[int, int]]]] = []
    ins_rows: list[InsertionRow] = []
    rep_rows: list[RepliconRow] = []
    counts = collections.Counter()
    for gid, rows in _by_genome(world / "replicons" / "genes.tsv"):
        recs = dict(intergenic._fasta_records(world / "replicons" / "fasta" / f"{gid}.fna"))
        node = next(n for n, x in gid_of.items() if x == gid)
        my = by_genome_items.get(gid, [])
        # each insertion's point, from its anchor copy on the chromosome
        planned: dict[int, list] = collections.defaultdict(list)
        for it in my:
            anc = [r for r in rows if int(r["family_id"][1:]) == it.anchor]
            if len(anc) != 1 or int(anc[0]["replicon"]) != 0:
                raise AssertionError(f"{it.item_id}: anchor F{it.anchor:06d} is not single-copy in {gid}")
            r = anc[0]
            rec = f"{gid}_r0"
            N = len(recs[rec])
            s, e = int(r["start"]), int(r["end"])
            at = (e if e <= N else e - N) if int(r["strand"]) == 1 else s
            content = it.tips[node]
            block = content if it.kind == "prophage" else content["seq"]
            if it.strand == -1:
                block = revcomp(block)
            planned[0].append((at, block, it, r["copy_label"]))
        by_rep = collections.defaultdict(list)
        for r in rows:
            by_rep[int(r["replicon"])].append(r)
        fa = sequences._Fasta(tmp / "fasta" / f"{gid}.fna")
        for rr in sorted(reps_s5[gid], key=lambda r: int(r["replicon"])):
            rep = int(rr["replicon"])
            rec = rr["record"]
            old = recs[rec]
            N = len(old)
            pl = planned.get(rep, [])
            new, starts = insert_blocks(old, [(a, b) for a, b, _, _ in pl])
            lift = Lift(N, [(a, len(b)) for a, b, _, _ in pl])
            if len(new) != lift.new_length:
                raise AssertionError(f"{rec}: inserted length disagrees with the lift")
            fa.write(rec, new)
            out_genes = []
            for r in by_rep[rep]:
                s2, e2 = lift.feature(int(r["start"]), int(r["end"]))
                out_genes.append(GeneRow(gid, rep, r["copy_label"], r["family_id"], "host", "", s2, e2,
                                         int(r["strand"]), e2 > len(new), r["cds_edit"], ""))
            for (a, block, it, anchor_copy), b0 in zip(pl, starts):
                iid = f"{it.item_id}@{gid}"
                L = len(block)
                row = InsertionRow(iid, it.kind, it.item_id, gid, rep, rec, a if a else N, L, b0, b0 + L, it.strand,
                                   f"F{it.anchor:06d}", anchor_copy)
                writers["insertions.tsv"].row(row)
                ins_rows.append(row)
                if it.kind == "prophage":
                    g = pgids[it.phage]
                    for pr in pgenes[g]:
                        s, e, st = flip((int(pr["start"]), int(pr["end"])), int(pr["strand"]), L, it.strand)
                        out_genes.append(GeneRow(gid, rep, f"ph_{pr['copy_label']}@{gid}", "P" + pr["family_id"][1:],
                                                 "prophage", iid, b0 + s, b0 + e, st, False, pr["cds_edit"],
                                                 puid[pr["copy_label"]]))
                    inst = it.tips[node]
                    n_sub = int(np.count_nonzero(np.frombuffer(inst.encode("ascii"), np.uint8)
                                                 != np.frombuffer(pseqs[it.phage].encode("ascii"), np.uint8)))
                    instances_pro.append(ProphageRow(it.item_id, phage_ids[g], gid, rec, b0, b0 + L, it.strand,
                                                     it.node, it.t_insert, f"F{it.anchor:06d}", n_sub,
                                                     q6(1.0 - n_sub / L)))
                else:
                    c = it.tips[node]
                    for lab, fid, span, edit, uid in c["genes"]:
                        s, e, st = flip(span, 1, L, it.strand)
                        out_genes.append(GeneRow(gid, rep, lab, fid, "cas", iid, b0 + s, b0 + e, st, False, edit, uid))
                    ls, le, _ = flip(c["leader"], 1, L, it.strand)
                    as_, ae, _ = flip(c["array"], 1, L, it.strand)
                    sp_spans = [flip(x, 1, L, it.strand)[:2] for x in c["spacers"]]
                    instances_arr.append((ArrayRow(it.item_id, gid, rec, f"cas_{it.item_id}", it.strand, b0, b0 + L,
                                                   b0 + ls, b0 + le, b0 + as_, b0 + ae, len(c["ids"]), c["repeat"],
                                                   it.node, it.t_insert, f"F{it.anchor:06d}",
                                                   tuple(c["ids"])), [(b0 + x, b0 + y) for x, y in sp_spans]))
                # targets inside this prophage instance
                if it.kind == "prophage":
                    for qi, a2, b2, st, mm, p in find_targets(new, queries, max_mm, pam,
                                                              circular=rr["topology"] == "circular", lo=b0,
                                                              hi=b0 + L):
                        targets[distinct[qi].spacer_id].append(Target("prophage", iid, rec, a2, b2, st, mm, p,
                                                                      pam_matches(p, pam)))
            out_genes.sort(key=lambda x: (x.start, x.end, x.copy_label))
            for x in out_genes:
                writers["genes.tsv"].row(x)
            rrow = RepliconRow(gid, rep, rec, rr["kind"], rr["topology"], len(new), N, len(out_genes), len(pl),
                               q6((new.count("G") + new.count("C")) / len(new)), ids.contig_uid(new.encode("ascii")).hex())
            writers["replicons.tsv"].row(rrow)
            rep_rows.append(rrow)
            counts["genes"] += len(out_genes)
            counts["replicons"] += 1
        fa.close()
        counts["genomes"] += 1
    for w in writers.values():
        w.close()

    # --- truth of prophages, Cas, arrays, spacers and targets
    t = _Rows(tmp / "prophages.tsv", ProphageRow)
    instances_pro.sort(key=lambda r: (r.prophage_id, r.genome_id))
    for r in instances_pro:
        t.row(r)
    t.close()
    t = _Rows(tmp / "cas.tsv", CasRow)
    for r in cas_rows:
        t.row(r)
    t.close()
    instances_arr.sort(key=lambda x: (x[0].array_id, x[0].genome_id))
    t = _Rows(tmp / "arrays.tsv", ArrayRow)
    for r, _ in instances_arr:
        t.row(r)
    t.close()
    spacer_rows = []
    t = _Rows(tmp / "spacers.tsv", SpacerRow)
    for s in distinct:
        row = SpacerRow(s.spacer_id, s.array_id, s.seq, s.node, q6(s.time), s.source,
                        phage_ids[pgids[s.phage]] if s.phage >= 0 else "", s.start, s.end, s.strand,
                        s.n_mismatches, s.pam_seq, s.pam_present)
        t.row(row)
        spacer_rows.append(row)
    t.close()
    by_id = {s.spacer_id: s for s in distinct}
    for s in distinct:
        tg = sorted(set(targets[s.spacer_id]), key=lambda x: (x.kind != "phage", x.target_id, x.start, x.strand))
        targets[s.spacer_id] = tg
        if s.source == "phage":
            want = (phage_ids[pgids[s.phage]], s.start, s.end, s.strand, s.n_mismatches)
            if not any((x.target_id, x.start, x.end, x.strand, x.n_mismatches) == want for x in tg):
                raise AssertionError(f"{s.spacer_id}: the scan misses its designed protospacer {want}")
    t = _Rows(tmp / "spacer_targets.tsv", SpacerTargetRow)
    n_target_rows = 0
    for ar, spans in instances_arr:
        for i, (sid, (ss, se)) in enumerate(zip(ar.spacer_ids, spans)):
            s = by_id[sid]
            head = (ar.array_id, ar.genome_id, ar.record, ar.array_start, ar.array_end, ar.strand, ar.repeat,
                    ar.cas_module_id, sid, i, ss, se, s.seq, s.node)
            if not targets[sid]:
                t.row(SpacerTargetRow(*head, "none", "", "", -1, -1, 0, -1, "", False, False))
                n_target_rows += 1
            for x in targets[sid]:
                designed = (s.source == "phage" and x.kind == "phage" and x.target_id == phage_ids[pgids[s.phage]]
                            and (x.start, x.strand) == (s.start, s.strand))
                t.row(SpacerTargetRow(*head, x.kind, x.target_id, x.record, x.start, x.end, x.strand,
                                      x.n_mismatches, x.pam_seq, x.pam_present, designed))
                n_target_rows += 1
    t.close()
    counts.update(phages=len(pgids), prophage_instances=len(instances_pro), array_instances=len(instances_arr),
                  spacers=len(distinct), spacer_target_rows=n_target_rows, insertions=len(ins_rows))

    files = {str(f.relative_to(tmp)): sequences._sha(f) for f in sorted(tmp.rglob("*")) if f.is_file()}
    digest = blake3("".join(f"{k}\t{v}\n" for k, v in sorted(files.items())).encode()).hexdigest(length=16)
    pre = {"phages": ph, "prophages": pro, "crispr": cr}
    pdocs = {d: json.loads((pw / d / f"{n}.json").read_text(encoding="utf-8"))["digest"]
             for d, n in (("backbone", "backbone"), ("sequences", "sequences"), ("replicons", "replicons"))}
    doc = {"format": FORMAT, "seed": seed, "preset": pre, "params_digest": backbone.params_digest(pre),
           "backbone_digest": bb["digest"], "sequences_digest": sq["digest"], "replicons_digest": rp["digest"],
           "phage_world": {"seed": phage_seed(seed), "digests": pdocs}, "versions": {"numpy": np.__version__},
           "stream_paths": ["S6/phages/world", "S6/placement", "S6/anchor/<item>", "S6/prophage/<prophage_id>",
                            "S6/crispr/<array_id>", "S3/family/F9<aa>0<k>", "S4/family/F9<aa>0<k>",
                            "S5/family/F9<aa>0<k>"],
           "items": {it.item_id: {"clade_node": it.node, "insertion_time": it.t_insert,
                                  "anchor_family": f"F{it.anchor:06d}", "strand": it.strand,
                                  "genomes": [gid_of[x] for x in it.clade.tips]} for it in items},
           "counts": dict(sorted(counts.items())), "files": files, "digest": digest}
    (tmp / "mobile.json").write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    tmp.rename(out)
    return read(world)


def _anchor_families(rows: list[dict], jrows: list[dict], kinds: dict[int, str]) -> list[int]:
    """The families of one genome that can anchor: one copy in the genome, on the chromosome, of a root or
    originated family, with no overlap at its 3' end (S5 junctions)."""
    n = collections.Counter(int(r["family_id"][1:]) for r in rows)
    over_after = {}     # (replicon, left_position) -> overlap
    over_before = {}    # (replicon, right_position) -> overlap
    for j in jrows:
        over_after[(int(j["replicon"]), int(j["left_position"]))] = int(j["overlap_bp"])
        over_before[(int(j["replicon"]), int(j["right_position"]))] = int(j["overlap_bp"])
    out = []
    for r in rows:
        f = int(r["family_id"][1:])
        if n[f] != 1 or int(r["replicon"]) != 0 or kinds.get(f) not in ANCHOR_KINDS:
            continue
        key = (0, int(r["position"]))
        o = over_after.get(key) if int(r["strand"]) == 1 else over_before.get(key)
        if o == 0:
            out.append(f)
    return out


def _crispr_locus(it: Item, a: int, cr: dict, seq_cfg: dict, igr_cfg: dict, kappa: float, seed: int, panel,
                  lin_seq, lin_igr, sites: PhageSites):
    """One array's Cas operon, leader, repeat and spacers at every clade tip (module docstring)."""
    cl = it.clade
    roles = list(cr["cas_roles"])
    J = igr_cfg["junctions"]
    leaves_by_role = []
    rows = []
    for k, role in enumerate(roles):
        fam = CAS_FAMILY_BASE + 10 * a + k
        d = sequences.design_family(seq_cfg, seed, fam, panel)
        nodes = gene_tree(cl, fam, f"{it.item_id}.{role}")
        res = sequences.evolve_family(d, nodes, lin_seq, seq_cfg, seeds.generator(seed, f"S4/family/{d.family_id}"))
        cas = intergenic.evolve_family(d.family_id, nodes, res.mu, lin_igr, igr_cfg, kappa, seed)
        leaves = {lf.lineage: lf for lf in res.leaves}
        cass = {c.lineage: c for c in cas}
        leaves_by_role.append((fam, d, leaves, cass))
        rows.append(CasRow(it.item_id, f"cas_{it.item_id}", role, k, fam, d.family_id, d.dark, d.n_domains,
                           tuple(x.accession for x in d.domains), d.target_identity, res.mu))
    rng = seeds.generator(seed, f"S6/crispr/{it.item_id}")
    R, S = int(cr["repeat_length"]), int(cr["spacer_length"])
    k = int(cr["repeat_stem"])
    stem = _ACGT[rng.integers(0, 4, size=k)].tobytes().decode("ascii")
    repeat = stem + _ACGT[rng.integers(0, 4, size=R - 2 * k)].tobytes().decode("ascii") + revcomp(stem)
    li = lin_igr.index[cl.node]
    leader = _ACGT[np.minimum(np.count_nonzero(lin_igr.comp_cum[li][None, :]
                                               <= rng.random(int(cr["leader_length"]))[:, None], axis=1), 3)]
    leader = leader.tobytes().decode("ascii")
    reg: list[Spacer] = []

    def acquire(node: str, time: float) -> int:
        reg.append(draw_spacer(cr, sites, rng, f"{it.item_id}.s{len(reg):03d}", it.item_id, node, time))
        return len(reg) - 1

    lists = evolve_array(cl, cr, acquire, rng)
    F = int(cr["locus_flank"])
    tips = {}
    for tip in cl.tips:
        genes = []
        for fam, d, leaves, cass in leaves_by_role:
            lf, c = leaves[tip], cass[tip]
            genes.append(intergenic.GeneIn(lf.copy_label, fam, 1, f"cas/{it.item_id}", lf.cds, c))
        juncs = [intergenic.junction(genes[i], genes[i + 1], J, seed) for i in range(len(genes) - 1)]
        flanks = (intergenic.end_flank(genes[0], F, "left"), intergenic.end_flank(genes[-1], F, "right"))
        asm = intergenic.assemble(genes, juncs, False, 0, flanks)
        sp = [reg[x] for x in lists[tip]]
        arr = "".join(repeat + s.seq for s in sp) + repeat
        Lc = len(asm.seq)
        a0 = Lc + len(leader)
        gene_info = []
        for (fam, d, leaves, cass), g, span, edit in zip(leaves_by_role, genes, asm.gene_spans, asm.edits):
            prot = leaves[tip].protein
            gene_info.append((g.copy_label, d.family_id, span, edit, blake3(prot.encode("ascii")).hexdigest(length=16)))
        tips[tip] = {"seq": asm.seq + leader + arr, "genes": gene_info, "leader": (Lc, a0),
                     "array": (a0, a0 + len(arr)), "spacers": [(a0 + R + i * (R + S), a0 + R + i * (R + S) + S)
                                                                for i in range(len(sp))],
                     "ids": [s.spacer_id for s in sp], "repeat": repeat}
    return tips, rows, reg


def read(world) -> MobileWorld:
    """Load ``world/mobile/``'s tables but genes and spacer_targets (which grow with the world)."""
    world = Path(world)
    mdir = world / "mobile"
    doc = json.loads((mdir / "mobile.json").read_text(encoding="utf-8"))
    if doc["format"] != FORMAT:
        raise ValueError(f"{mdir}: format {doc['format']!r}, expected {FORMAT!r}")
    t = {}
    for name in ("phages.tsv", "replicons.tsv", "insertions.tsv", "prophages.tsv", "cas.tsv", "arrays.tsv",
                 "spacers.tsv"):
        cls = _TABLES[name]
        conv = backbone._parser(cls)
        with open(mdir / name, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            if tuple(header) != cls._fields:
                raise ValueError(f"{mdir / name}: header {header}, expected {list(cls._fields)}")
            t[name] = [cls(*(c(v) for c, v in zip(conv, row))) for row in reader]
    return MobileWorld(doc["seed"], t["phages.tsv"], t["replicons.tsv"], t["insertions.tsv"], t["prophages.tsv"],
                       t["cas.tsv"], t["arrays.tsv"], t["spacers.tsv"], doc["digest"], world, doc)


__all__ = ["FORMAT", "CAS_FAMILY_BASE", "Lift", "KmerIndex", "PhageSites", "Spacer", "SpeciesTree", "Clade",
           "MobileWorld", "build", "read", "phage_preset", "phage_seed", "species_tree", "clade", "gene_tree",
           "insert_blocks", "flip", "coding_sites", "read_codons", "write_codons", "evolve_prophage", "mutate",
           "draw_spacer", "evolve_array", "scan", "find_targets", "pam_at", "block_length", "pam_matches",
           "to_codes", "check_config"]
