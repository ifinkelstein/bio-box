"""synthworld intergenic regions, regulatory elements, RNA motifs and replicon sequences
(planning/tests/synthworld_spec.md §2.1 S5, §2.4, §2.5; card F-02 session 4). Reads a built backbone
(``synthworld-backbone/2``) and its sequences (``synthworld-sequences/1``), never ZOMBI2.

**Cassettes (§2.4).** Every gene copy carries an upstream cassette U and a downstream cassette D, laid out on
its coding strand. U (``cassettes.upstream_length`` bases, ending at the start codon) holds, 5'->3': background,
a sigma70-like -35 box (TTGACA), a spacer of 15-19, a -10 box (TATAAT), the 5' UTR's first part, the RNA
motif of a linked family (§2.5; empty otherwise), the UTR's second part, an RBS (AGGAGG) and a spacer of
``rbs_spacer`` bases before the start codon. D (``downstream_length`` bases, starting after the stop codon)
holds a gap, a terminator (a GC-rich stem, a loop, the closing stem and a T tract) and background. A
family's layout (spacer, UTR parts, RBS spacer, gap, stem and loop lengths) is drawn once, at its root, from
stream ``S5/family/<family_id>``; with a motif of m bases, U's first m bases are dropped, so U keeps its
length and the motif pushes the promoter upstream. Each copy also carries three length scores z_U, z_T,
z_D (standard normal at the root).

**Evolution along the gene tree (§2.4).** The family's complete gene tree (``backbone/gene_tree_nodes.tsv``)
pruned to the ancestors of extant copies, as in S4, from the same stream: on a branch of length t in
lineage b (clock c_b, ``sequences/lineages.tsv``), with nu = ``intergene_speed`` x kappa_syn x mu_f
(mu_f the family's rate, ``sequences/families.tsv``; so a cassette site's event rate is ``intergene_speed``
times a codon's synonymous event rate), a background site has a substitution event with probability
1 - exp(-nu c_b t) and takes a base from lineage b's intergenic composition; an element site (box, RBS, T
tract) has rate ``element_rate_factor`` x nu and takes a base from its column (the consensus base with
probability ``element_consensus``, each other base a third of the rest); a terminator stem pair has the
element rate and is replaced by a canonical pair drawn from ``pair_weights``. The length scores change at
length events, Poisson at ``length_event_rate`` per unit (unscaled) time; at each, z <- rho z + sqrt(1 -
rho^2) N(0, 1) with rho = ``length_event_rho``, so a copy's preferred lengths are inherited unchanged
between events and with small multiplicative noise at them, and stay N(0, 1) at every node. Lineage b's intergenic GC is expit(``gc_slope`` logit(GC3_b) + ``gc_offset``),
GC3_b from S4; a site's base is A, C, G, T with probabilities (1 - g)/2, g/2, g/2, (1 - g)/2.

**RNA motifs (§2.5).** ``rna_motifs.n_models`` models (`rna`: synthetic, or Rfam CMs through cmemit) are
linked to ``families_per_model`` families each, drawn from stream ``S5/motif/links`` among the families of
``family_kinds`` present in at least ``min_genomes`` genomes (so no planted, plasmid, vehicle, hub or decoy
family carries a motif). Each linked family's tree is cut at height h = ``clade_height`` x the species crown
age into clades (each branch crossing that height opens one; a family younger than h is one clade, opened
at its origination); clade k of link l gets emission ``<l>.<k>`` (preorder of the crossing branches), from
``cmemit`` with tool seed ``S5/emit/<link>`` or from the synthetic model with stream ``S5/emit/<link>``. On the
crossing branch the cassette evolves to the cut, the emission is installed in U, and the motif then evolves
with the cassette at ``rna_motifs.rate_factor`` x nu (`rna.evolve_motif`: unpaired positions take the lineage's
composition, pairs change only to canonical pairs). Motifs sit in the 5' UTR of their gene, on its strand;
none is planted inside a coding sequence (§2.5 names none there).

**Junctions (§2.4)**, between neighbouring genes on a replicon in the realized order of ``backbone/genes.tsv``
(plasmids excised; on a circular replicon also the last gene and the first): by strands co-directional,
divergent (<- ->) or convergent (-> <-). A co-directional junction is ``codir_tight`` when both genes are in
one unit of a planted group (``backbone/groups.tsv`` kind planted; plasmid and vehicle units take the
background rule), else when the transcriptionally downstream gene carries a motif it is ``codir_loose``
(a motif marks an operon's first gene), else it is tight when u < ``p_tight`` with u from BLAKE3 of
``<seed>/S5/tight/<upstream family>/<downstream family>`` (the families in transcription order, so an
inverted block keeps its classes and the same adjacency gets the same class in every lineage), else loose.
Lengths come from the junction's distribution, L = clip(round(offset + scale exp(sigma z)), min, max) with
z one of the copies' scores:
- ``codir_tight``: u = Phi(z_T of the downstream gene); u below ``tight_overlap.minus4`` tries a 4-base
  overlap (``ATGA``: the upstream stop TGA overlaps the downstream start), u below ``minus4 + minus1`` a
  1-base overlap (``TAATG``: the stop's last base is the start's first), else L from ``tight`` at
  z = Phi^-1((u - minus4 - minus1) / (1 - minus4 - minus1)). An overlap is realized with protein-neutral edits
  only: the upstream stop codon becomes another stop, the downstream start codon another initiation codon
  (-1: TAA/TGA before ATG, TAG before GTG, a TTG start becomes ATG; -4: the stop becomes TGA and the start
  Z+TG, Z the third base of the upstream gene's last sense codon). When that cannot hold (-4 with Z = C, or a
  downstream second codon not starting with A) the junction takes L from ``tight`` at z = Phi^-1(u / (minus4
  + minus1)) (its row keeps the tried overlap as ``target_length``). The sequence is the last L bases of the downstream U (RBS part).
- ``codir_loose``: L from ``loose`` at the downstream z_U; the sequence is the downstream U's last L bases,
  preceded by the upstream D's terminator span when L also holds the downstream promoter (§2.1 S5: an
  operon's last gene may carry a terminator).
- ``divergent``: L from ``divergent`` at (z_U + z_U')/sqrt 2, split between revcomp(U_left) and U_right in
  proportion to their ``loose`` lengths; ``convergent``: L from ``convergent`` at (z_D + z_D')/sqrt 2, split
  between D_left and revcomp(D_right) in proportion to their ``convergent`` lengths.
A part never exceeds its cassette; a motif-carrying gene's part is raised to hold its motif (rule
``motif_min``). A linear replicon has an end flank of ``linear_end_flank`` bases at each end (the end gene's
U or D facing it). Overlaps exist only in ``codir_tight`` junctions.

**Replicons.** One sequence per replicon: genes (the CDS of ``sequences/cds.fna`` with the edits above,
reverse-complemented on strand -1) and junctions in order. A circular replicon is cut at a random origin
(stream ``S5/origin/<genome_id>/<replicon>``, uniform over its bases), so a feature across the origin has
end > length (AGENTS §3.1); topology is the backbone's (``backbone/replicons.tsv``).

Shape on disk (``build(world, rfam_cm=None)`` writes ``world/replicons/``; ``read(world)`` loads the small
tables)::

    fasta/<genome_id>.fna, .fai  one record per replicon, ``<genome_id>_r<replicon>``, 80 columns
    replicons.tsv   genome_id replicon record kind topology length n_genes origin_offset gc contig_uid
    genes.tsv       genome_id replicon position copy_label family_id start end strand origin_spanning
                    operon_id start_codon stop_codon cds_edit
    junctions.tsv   genome_id replicon left_position right_position left_copy right_copy junction rule
                    target_length length overlap_bp start end igr_uid gc
    elements.tsv    genome_id replicon left_position kind start end strand copy_label rna_id
    rna.tsv         rna_id link model_acc model_name source emission_id exp local_mode attached_copy
                    family_id genome_id replicon start end strand length dist_to_start
                    n_subst_from_emission n_compensatory structure
    emissions.tsv   emission_id link model_acc clade_node height sequence structure
    links.tsv       link model_index model_acc model_name family family_id n_clades clade_height
    models.tsv      index accession name source length ga n_pairs structure
    lineages.tsv    lineage gc3 gc_igr
    replicons.json  format, seed, the ``intergenic`` and ``rna_motifs`` presets and their digest, the
                    backbone and sequences digests, the motif source (file name, sha256, Rfam release),
                    versions, stream paths, sha256 of every file above, digest

Rows are sorted by (genome_id, replicon, position); ``emissions`` and ``links`` by link. Coordinates are
0-based half-open on the + strand of the replicon record; a feature across the origin of a circular
replicon has end > length. ``origin_offset`` is where coordinate 0 lies, counted from the start of the gene
at position 0 in the uncut circle. ``junctions``: one row per junction (``left_position`` the gene before
it in replicon order, ``right_position`` the one after, -1 at a linear end, whose rows have junction
``left_end`` or ``right_end``); an overlap has start = end = the left gene's end and ``overlap_bp`` > 0;
``igr_uid`` is BLAKE3-128 of the IGR's + strand bases (``biobox.ids.contig_uid``), empty when it has none;
``target_length`` is the drawn length before a motif raised it (-1 or -4 for a tried overlap, realized when
``overlap_bp`` > 0); ``rule`` says why the class: ``group`` (one planted unit), ``motif`` (a motif-carrying
downstream gene), ``drawn`` (the hash rule, or the strands), ``motif_min`` (a motif raised the length) or
``end`` (a linear end flank).
``elements``: the planted elements that lie wholly inside an IGR (an element a short IGR trims away is not
there), by junction; strand is the owning gene's. ``rna``: one row per motif instance, ``dist_to_start`` the
bases from the motif's end to the start codon, ``n_subst_from_emission`` and ``n_compensatory`` the positions
that differ from the clade's emission and the pairs whose two bases both differ, ``structure`` the
instance's pairs (dot-bracket, 5'->3' on its strand). ``cds_edit`` lists ``start:<from>><to>`` and
``stop:<from>><to>`` edits (empty when none). ``operon_id`` is ``<genome_id>_r<replicon>_o<position>``,
the position of the operon's first gene in replicon order: an operon is a maximal chain of genes joined by
``codir_tight`` junctions. Floats are written at 6 significant digits. The per-genome working files live in
``replicons.partial/scratch`` while building and are removed before publishing; at any time one family
(evolving) or one genome (assembling) is in memory. ``digest`` is BLAKE3-128 over the sorted (path, sha256)
list of every file but replicons.json.
"""

from __future__ import annotations

import csv
import json
import math
import shutil
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

import numpy as np
from blake3 import blake3

from biobox import ids
from biobox.synthworld import backbone, rna, seeds, sequences

FORMAT = "synthworld-replicons/1"
q6 = sequences.q6
ELEMENT_KINDS = ("promoter_m35", "promoter_m10", "rbs", "terminator", "rna_motif")
JUNCTIONS = ("codir_tight", "codir_loose", "divergent", "convergent", "left_end", "right_end")
M35, M10, RBS = "TTGACA", "TATAAT", "AGGAGG"
STOPS = ("TAA", "TAG", "TGA")
_STARTS_OK = ("ATG", "GTG", "TTG")
_COMP = str.maketrans("ACGT", "TGCA")
_ND = statistics.NormalDist()
_FLUSH_LINES = 4096


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def phi(z: float) -> float:
    return _ND.cdf(float(z))


def phi_inv(u: float) -> float:
    return _ND.inv_cdf(min(max(float(u), 1e-12), 1.0 - 1e-12))


def length_of(dist: dict, z: float) -> int:
    """clip(round(offset + scale exp(sigma z)), min, max), the value rounded at 6 significant digits first."""
    x = q6(float(dist["offset"]) + float(dist["scale"]) * math.exp(float(dist["sigma"]) * float(z)))
    return int(min(max(math.floor(x + 0.5), int(dist["min"])), int(dist["max"])))


def gc_igr(gc3: float, comp: dict) -> float:
    g = float(gc3)
    x = float(comp["gc_slope"]) * math.log(g / (1.0 - g)) + float(comp["gc_offset"])
    return q6(1.0 / (1.0 + math.exp(-x)))


def composition(g: float) -> np.ndarray:
    return np.array([(1.0 - g) / 2.0, g / 2.0, g / 2.0, (1.0 - g) / 2.0])


def tight_u(seed: int, up: int, down: int) -> float:
    """The background tight-junction draw of an ordered family pair (transcription order), in [0, 1)."""
    d = blake3(f"{int(seed)}/S5/tight/{int(up)}/{int(down)}".encode("ascii")).digest(length=8)
    return int.from_bytes(d, "big") / 2.0**64


# --- cassette layout --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Layout:
    u_len: int
    d_len: int
    sp1: int           # -35 box end to -10 box start
    utr_a: int         # -10 box end to the motif
    utr_b: int         # the motif to the RBS
    sp2: int           # RBS end to the start codon
    term_start: int    # stop codon end to the terminator
    stem: int
    loop: int
    t_tract: int

    # positions in U without its motif (m = 0); `u_elements` shifts them for a motif of m bases
    @property
    def rbs_end(self) -> int:
        return self.u_len - self.sp2

    @property
    def rbs_start(self) -> int:
        return self.rbs_end - len(RBS)

    @property
    def mi(self) -> int:
        """Where the motif is inserted into U."""
        return self.rbs_start - self.utr_b

    @property
    def m10_start(self) -> int:
        return self.mi - self.utr_a - len(M10)

    @property
    def m35_start(self) -> int:
        return self.m10_start - self.sp1 - len(M35)

    @property
    def term_end(self) -> int:
        return self.term_start + 2 * self.stem + self.loop + self.t_tract

    def u_elements(self, m: int) -> list[tuple[str, int, int]]:
        """(kind, start, end) in the realized U (its first m bases dropped after the motif's insertion)."""
        out = [("promoter_m35", self.m35_start - m, self.m35_start - m + len(M35)),
               ("promoter_m10", self.m10_start - m, self.m10_start - m + len(M10))]
        if m:
            out.append(("rna_motif", self.mi - m, self.mi))
        out.append(("rbs", self.rbs_start, self.rbs_end))
        return out

    def d_elements(self) -> list[tuple[str, int, int]]:
        return [("terminator", self.term_start, self.term_end)]

    def need_motif(self, m: int) -> int:
        """Bases of U, from the start codon back, that hold the motif (0 without one)."""
        return self.u_len - (self.mi - m) if m else 0

    def need_promoter(self, m: int) -> int:
        return self.u_len - (self.m35_start - m)


def _geometric(rng, mean: float) -> int:
    m = float(mean)
    return int(rng.geometric(1.0 / (m + 1.0))) - 1 if m > 0 else 0


def _choice(rng, values, weights) -> int:
    w = np.cumsum(np.asarray(weights, dtype=np.float64))
    return int(values[min(int(np.searchsorted(w / w[-1], rng.random(), side="right")), len(values) - 1)])


def draw_layout(cas: dict, rng) -> Layout:
    lo, hi = (int(x) for x in cas["rbs_spacer"])
    utr_max = int(cas["utr_max"])
    lay = Layout(
        u_len=int(cas["upstream_length"]), d_len=int(cas["downstream_length"]),
        sp1=_choice(rng, cas["promoter_spacer"]["values"], cas["promoter_spacer"]["weights"]),
        utr_a=min(int(cas["utr_min"]) + _geometric(rng, cas["utr_mean"]), utr_max),
        utr_b=min(int(cas["utr_min"]) + _geometric(rng, cas["utr_mean"]), utr_max),
        sp2=int(rng.integers(lo, hi + 1)),
        term_start=min(_geometric(rng, cas["terminator_gap_mean"]), int(cas["terminator_gap_max"])),
        stem=int(rng.integers(int(cas["stem"][0]), int(cas["stem"][1]) + 1)),
        loop=int(rng.integers(int(cas["loop"][0]), int(cas["loop"][1]) + 1)),
        t_tract=int(cas["t_tract"]))
    if lay.term_end > lay.d_len:
        raise ValueError(f"downstream_length {lay.d_len} cannot hold a terminator ending at {lay.term_end}")
    return lay


# --- position weight columns --------------------------------------------------------------------------------

def element_columns(cas: dict) -> np.ndarray:
    """(1 + 19, 4) cumulative columns: row 0 a placeholder for the lineage's composition, rows 1-6 the -35
    box, 7-12 the -10 box, 13-18 the RBS, 19 the T tract."""
    p = float(cas["element_consensus"])
    rows = [np.full(4, 0.25)]
    for cons in (M35, M10, RBS, "T"):
        for b in cons:
            r = np.full(4, (1.0 - p) / 3.0)
            r[rna.BASES.index(b)] = p
            rows.append(r)
    tab = np.cumsum(np.vstack(rows), axis=1)
    tab[:, -1] = 1.0
    return tab


_ROW_M35, _ROW_M10, _ROW_RBS, _ROW_TT = 1, 7, 13, 19


def site_model(lay: Layout, element_rate: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Over U (without motif) then D: each site's column row (0 background) and rate factor (1 background,
    ``element_rate`` element, 0 a stem base, which changes with its pair), and the stem pairs (k, 2)."""
    n = lay.u_len + lay.d_len
    row = np.zeros(n, dtype=np.int64)
    rate = np.ones(n)
    for start, first in ((lay.m35_start, _ROW_M35), (lay.m10_start, _ROW_M10), (lay.rbs_start, _ROW_RBS)):
        row[start:start + 6] = np.arange(first, first + 6)
        rate[start:start + 6] = element_rate
    d0 = lay.u_len + lay.term_start
    s2 = d0 + lay.stem + lay.loop
    tt = s2 + lay.stem
    row[tt:tt + lay.t_tract] = _ROW_TT
    rate[tt:tt + lay.t_tract] = element_rate
    pairs = np.array([(d0 + i, s2 + lay.stem - 1 - i) for i in range(lay.stem)], dtype=np.int64).reshape(-1, 2)
    rate[pairs.ravel()] = 0.0
    if lay.m35_start < 0:
        raise ValueError(f"upstream_length {lay.u_len} is too short for the layout {lay}")
    return row, rate, pairs


# --- evolution of one family's cassettes --------------------------------------------------------------------

@dataclass
class Lineages:
    index: dict[str, int]
    clock: np.ndarray
    gc3: np.ndarray
    gc_igr: np.ndarray
    comp_cum: np.ndarray       # (n, 4)


def lineages(world: Path, comp: dict) -> Lineages:
    rows = _tsv(world / "sequences" / "lineages.tsv")
    gcs = [gc_igr(float(r["gc3"]), comp) for r in rows]
    cum = np.cumsum(np.stack([composition(g) for g in gcs]), axis=1)
    cum[:, -1] = 1.0
    return Lineages({r["lineage"]: i for i, r in enumerate(rows)}, np.array([float(r["clock"]) for r in rows]),
                    np.array([float(r["gc3"]) for r in rows]), np.array(gcs), cum)


@dataclass
class Link:
    link: str
    model: rna.MotifModel
    family: int
    family_id: str
    n_clades: int = 0


@dataclass
class CopyCassette:
    copy_label: str
    lineage: str
    u: str
    d: str
    z: tuple[float, float, float]      # z_U, z_T, z_D
    layout: Layout
    m: int
    link: str = ""
    emission: str = ""
    n_subst: int = 0
    n_comp: int = 0
    structure: str = ""


@dataclass
class _Node:
    nt: np.ndarray
    z: np.ndarray
    motif: np.ndarray | None = None
    emission: int = -1


@dataclass
class EmissionRow:
    emission_id: str
    link: str
    model_acc: str
    clade_node: int
    height: float
    emission: rna.Emission


def _draw(cum_rows: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Per row, the first base whose cumulative probability exceeds u (rows over ACGT, last column 1)."""
    return (cum_rows[:, :3] <= u[:, None]).sum(axis=1).astype(np.uint8)


def _cell(v) -> str:
    return "true" if v is True else "false" if v is False else str(v)


class _Rows:
    """A TSV table written row by row (`backbone._cell`'s format for the types used here)."""

    def __init__(self, path: Path, cls):
        self.fh = open(path, "w", encoding="utf-8", newline="")
        self.fh.write("\t".join(cls._fields) + "\n")

    def row(self, r) -> None:
        self.fh.write("\t".join(map(_cell, r)) + "\n")

    def close(self) -> None:
        self.fh.close()


def evolve_family(fid: str, nodes, mu: float, lin: Lineages, cfg: dict, kappa: float, seed: int,
                  link: Link | None = None, *, clade_height: float = 0.0, emit=None,
                  emissions_out: list | None = None) -> list[CopyCassette]:
    """The extant copies' cassettes of one family (module docstring), from stream ``S5/family/<fid>``.
    ``emit(n)`` gives a linked family's n clade emissions (`rna.Emission`)."""
    cas, ev = cfg["cassettes"], cfg["evolution"]
    rng = seeds.generator(seed, f"S5/family/{fid}")
    lay = draw_layout(cas, rng)
    efac = float(ev["element_rate_factor"])
    row, rate, stem_pairs = site_model(lay, efac)
    cols = element_columns(cas)
    pw = cas["pair_weights"]
    pair_cum = np.cumsum([float(pw[p]) for p in rna.CANONICAL])
    pair_cum = pair_cum / pair_cum[-1]
    pair_cum[-1] = 1.0
    ev_rate = float(cfg["junctions"]["length_event_rate"])
    jump_rho = float(cfg["junctions"]["length_event_rho"])
    nu = float(ev["intergene_speed"]) * float(kappa) * float(mu)
    mfac = float(cfg["_motif_rate_factor"])

    by_id = {n.node: n for n in nodes}
    keep: set[int] = set()
    for n in nodes:
        if n.extant:
            x = n.node
            while x != -1 and x not in keep:
                keep.add(x)
                x = by_id[x].parent
    kept = sorted(keep)
    if not kept:
        raise ValueError(f"{fid}: no extant copy")
    n_children = {x: 0 for x in kept}
    for x in kept:
        if by_id[x].parent != -1:
            n_children[by_id[x].parent] += 1
    t_end: dict[int, float] = {}
    for x in kept:
        par = by_id[x].parent
        t_end[x] = (t_end[par] if par != -1 else 0.0) + float(by_id[x].branch_length)
    present = max(t_end[x] for x in kept if by_id[x].extant)

    # clades (§2.5 step 1): the branches crossing height h, in preorder; a family younger than h is one clade
    cut: dict[int, float] = {}            # node -> time above the cut on its branch
    if link is not None:
        h = float(clade_height)
        root = kept[0]
        if present - (t_end[root] - float(by_id[root].branch_length)) <= h:
            cut[root] = 0.0
        else:
            for x in kept:
                lo = present - t_end[x]
                hi = lo + float(by_id[x].branch_length)
                if lo <= h < hi:
                    cut[x] = hi - h
        em = emit(len(cut))
        link.n_clades = len(cut)
        clade_of = {x: k for k, x in enumerate(sorted(cut))}
        for x, k in clade_of.items():
            emissions_out.append(EmissionRow(f"{link.link}.{k}", link.link, link.model.accession, x,
                                             q6(present - t_end[x] + float(by_id[x].branch_length) - cut[x]),
                                             em[k]))
        em_codes = [rna.to_codes(e.seq) for e in em]
        em_pairs = [np.array(e.pairs, dtype=np.int64).reshape(-1, 2) for e in em]
        em_unp = []
        for e in em:
            paired = np.zeros(len(e.seq), dtype=bool)
            if e.pairs:
                paired[np.array(e.pairs).ravel()] = True
            em_unp.append(np.nonzero(~paired)[0])

    # the root's cassettes
    li0 = lin.index[by_id[kept[0]].lineage]
    S = lay.u_len + lay.d_len
    tab0 = cols.copy()
    tab0[0] = lin.comp_cum[li0]
    nt = _draw(tab0[row], rng.random(S))
    if len(stem_pairs):
        k = np.minimum(np.count_nonzero(pair_cum[None, :] <= rng.random(len(stem_pairs))[:, None], axis=1), 5)
        nt[stem_pairs[:, 0]] = rna._PAIR_A[k]
        nt[stem_pairs[:, 1]] = rna._PAIR_B[k]
    root = _Node(nt, rng.standard_normal(3))

    n_pairs = len(stem_pairs)
    tables: dict[int, np.ndarray] = {}

    def branch(st: _Node, t: float, li: int) -> None:
        if t <= 0.0:
            return
        x = nu * float(lin.clock[li]) * t
        hit = np.flatnonzero(rng.random(S) < -np.expm1(-x * rate))
        if hit.size:
            tab = tables.get(li)
            if tab is None:
                tab = tables[li] = cols.copy()
                tab[0] = lin.comp_cum[li]
            st.nt[hit] = _draw(tab[row[hit]], rng.random(hit.size))
        if n_pairs:
            ph = np.flatnonzero(rng.random(n_pairs) < -math.expm1(-x * efac))
            if ph.size:
                c = np.minimum(np.count_nonzero(pair_cum[None, :] <= rng.random(ph.size)[:, None], axis=1), 5)
                st.nt[stem_pairs[ph, 0]] = rna._PAIR_A[c]
                st.nt[stem_pairs[ph, 1]] = rna._PAIR_B[c]
        if st.motif is not None:
            rna.evolve_motif(st.motif, em_pairs[st.emission], em_unp[st.emission], x * mfac, lin.comp_cum[li],
                             pair_cum, rng)
        k = int(rng.poisson(ev_rate * t))
        if k:
            r = jump_rho ** k                      # k AR(1) steps compose to one with rho^k
            st.z = r * st.z + math.sqrt(1.0 - r * r) * rng.standard_normal(3)

    states: dict[int, _Node] = {}
    out: list[CopyCassette] = []
    leaves = {}
    for x in kept:
        n = by_id[x]
        par = n.parent
        src = root if par == -1 else states[par]
        if par == -1 or n_children[par] > 1:
            st = _Node(src.nt.copy(), src.z.copy(), None if src.motif is None else src.motif.copy(), src.emission)
        else:
            st = src
        if par != -1:
            n_children[par] -= 1
            if n_children[par] == 0:
                del states[par]
        li = lin.index[n.lineage]
        bl = float(n.branch_length)
        if x in cut:
            branch(st, cut[x], li)
            k = clade_of[x]
            st.motif, st.emission = em_codes[k].copy(), k
            branch(st, bl - cut[x], li)
        else:
            branch(st, bl, li)
        if n_children[x]:
            states[x] = st
        if n.extant:
            leaves[x] = st
    for x in sorted(leaves, key=lambda y: by_id[y].copy_label):
        st = leaves[x]
        u = rna.to_str(st.nt[:lay.u_len])
        d = rna.to_str(st.nt[lay.u_len:])
        m = 0
        rec = {}
        if st.motif is not None:
            ms = rna.to_str(st.motif)
            m = len(ms)
            u = (u[:lay.mi] + ms + u[lay.mi:])[m:]
            e = em[st.emission]
            diff = st.motif != em_codes[st.emission]
            prs = em_pairs[st.emission]
            comp_ = int(np.count_nonzero(diff[prs[:, 0]] & diff[prs[:, 1]])) if len(prs) else 0
            rec = dict(link=link.link, emission=f"{link.link}.{st.emission}", n_subst=int(diff.sum()),
                       n_comp=comp_, structure=rna.dot_bracket(m, e.pairs))
        z = tuple(q6(v) for v in st.z)
        out.append(CopyCassette(by_id[x].copy_label, by_id[x].lineage, u, d, z, lay, m, **rec))
    return out


# --- junctions ----------------------------------------------------------------------------------------------

class GeneIn(NamedTuple):
    """One gene of a replicon as the junction rules see it."""
    copy_label: str
    family: int
    strand: int
    unit_key: str          # "<group>/<unit>" of a planted unit, else ""
    cds: str               # coding strand, start to stop
    cas: CopyCassette


@dataclass
class Junction:
    kind: str
    rule: str              # why the class and length: group | drawn | motif | motif_min | end
    target: int
    length: int
    overlap: int
    seq: str               # + strand IGR bases
    elements: list[tuple[str, int, int, int, str]]   # kind, start, end (in seq), strand, copy_label
    edits: dict[str, dict[str, str]] = field(default_factory=dict)   # copy_label -> {"start"|"stop": codon}


def _u_part(g: GeneIn, ell: int) -> tuple[str, list]:
    """The last ``ell`` bases of g's U and the elements wholly inside them (part coordinates)."""
    lay = g.cas.layout
    off = lay.u_len - ell
    els = [(k, s - off, e - off) for k, s, e in lay.u_elements(g.cas.m) if s >= off]
    return g.cas.u[off:], els


def _d_part(g: GeneIn, ell: int) -> tuple[str, list]:
    els = [(k, s, e) for k, s, e in g.cas.layout.d_elements() if e <= ell]
    return g.cas.d[:ell], els


class _Builder:
    def __init__(self):
        self.parts: list[str] = []
        self.elements: list = []
        self.n = 0

    def add(self, g: GeneIn, part: tuple[str, list], reverse: bool) -> None:
        seq, els = part
        ell = len(seq)
        if reverse:
            self.parts.append(revcomp(seq))
            self.elements += [(k, self.n + ell - e, self.n + ell - s, g.strand, g.copy_label) for k, s, e in els]
        else:
            self.parts.append(seq)
            self.elements += [(k, self.n + s, self.n + e, g.strand, g.copy_label) for k, s, e in els]
        self.n += ell

    def seq(self) -> str:
        return "".join(self.parts)


def _split(L: int, wa: float, wb: float, cap_a: int, cap_b: int, need_a: int, need_b: int) -> tuple[int, int, bool]:
    la = int(math.floor(L * wa / (wa + wb) + 0.5))
    la = min(max(la, L - cap_b), cap_a)
    lb = L - la
    raised = False
    if la < need_a:
        la, raised = need_a, True
    if lb < need_b:
        lb, raised = need_b, True
    return la, lb, raised


def overlap_edits(up_cds: str, down_cds: str, o: int) -> dict[str, str] | None:
    """The protein-neutral edits that make the upstream stop and the downstream start overlap by ``o`` (1 or
    4) bases (module docstring), as {"stop": codon, "start": codon} (unchanged codons included), or None."""
    stop, start = up_cds[-3:], down_cds[:3]
    if o == 1:
        if start == "TTG":
            start = "ATG"
        if start == "ATG":
            stop = stop if stop in ("TAA", "TGA") else "TAA"
        elif start == "GTG":
            stop = "TAG"
        else:
            return None
        return {"stop": stop, "start": start}
    if o == 4:
        z = up_cds[-4]
        if z == "C" or len(down_cds) < 7 or down_cds[3] != "A":
            return None
        return {"stop": "TGA", "start": z + "TG"}
    raise ValueError(f"overlap {o}: only 1 and 4 are planted")


def junction(a: GeneIn, b: GeneIn, J: dict, seed: int, *, same: bool = False) -> Junction:
    """The junction between ``a`` (left in replicon order) and ``b`` (right); ``same`` when they are one gene
    (a circular replicon of one gene). Module docstring."""
    la_ = a.cas.layout
    if a.strand == b.strand:
        up, down = (a, b) if a.strand == 1 else (b, a)
        rev = a.strand == -1
        if not same and a.unit_key and a.unit_key == b.unit_key:
            kind, rule = "codir_tight", "group"
        elif down.cas.m:
            kind, rule = "codir_loose", "motif"
        elif not same and tight_u(seed, up.family, down.family) < float(J["p_tight"]):
            kind, rule = "codir_tight", "drawn"
        else:
            kind, rule = "codir_loose", "drawn"
        bld = _Builder()
        if kind == "codir_tight":
            ov = J["tight_overlap"]
            p4, p1 = float(ov["minus4"]), float(ov["minus1"])
            u = phi(down.cas.z[1])
            o = 4 if u < p4 else 1 if u < p4 + p1 else 0
            if o:
                ed = overlap_edits(up.cds, down.cds, o)
                if ed is not None:
                    edits = {up.copy_label: {"stop": ed["stop"]}, down.copy_label: {"start": ed["start"]}}
                    return Junction(kind, rule, -o, 0, o, "", [], edits)
                L = length_of(J["tight"], phi_inv(u / (p4 + p1)))       # the overlap cannot hold
                target = -o
            else:
                target = L = length_of(J["tight"], phi_inv((u - p4 - p1) / (1.0 - p4 - p1)))
            parts = [(down, _u_part(down, L))]
        else:
            target = L = length_of(J["loose"], down.cas.z[0])
            need = down.cas.layout.need_motif(down.cas.m)
            if L < need:
                L, rule = need, "motif_min"
            k = up.cas.layout.term_end
            if L >= k + down.cas.layout.need_promoter(down.cas.m) and L - k <= down.cas.layout.u_len:
                parts = [(up, _d_part(up, k)), (down, _u_part(down, L - k))]
            else:
                parts = [(down, _u_part(down, L))]
        for g, part in (reversed(parts) if rev else parts):
            bld.add(g, part, rev)
        return Junction(kind, rule, target, L, 0, bld.seq(), bld.elements)
    if a.strand == -1:      # divergent: <- a | b ->
        z = (a.cas.z[0] + b.cas.z[0]) / math.sqrt(2.0)
        target = L = length_of(J["divergent"], z)
        L = min(L, la_.u_len + b.cas.layout.u_len)
        wa, wb = length_of(J["loose"], a.cas.z[0]), length_of(J["loose"], b.cas.z[0])
        x, y, raised = _split(L, wa, wb, la_.u_len, b.cas.layout.u_len, la_.need_motif(a.cas.m),
                              b.cas.layout.need_motif(b.cas.m))
        bld = _Builder()
        bld.add(a, _u_part(a, x), True)
        bld.add(b, _u_part(b, y), False)
        return Junction("divergent", "motif_min" if raised else "drawn", target, x + y, 0, bld.seq(), bld.elements)
    z = (a.cas.z[2] + b.cas.z[2]) / math.sqrt(2.0)        # convergent: a -> | <- b
    target = L = length_of(J["convergent"], z)
    L = min(L, la_.d_len + b.cas.layout.d_len)
    wa, wb = length_of(J["convergent"], a.cas.z[2]), length_of(J["convergent"], b.cas.z[2])
    x, y, _ = _split(L, max(wa, 1), max(wb, 1), la_.d_len, b.cas.layout.d_len, 0, 0)
    bld = _Builder()
    bld.add(a, _d_part(a, x), False)
    bld.add(b, _d_part(b, y), True)
    return Junction("convergent", "drawn", target, x + y, 0, bld.seq(), bld.elements)


def end_flank(g: GeneIn, F: int, side: str) -> Junction:
    """A linear replicon's end flank: the end gene's U or D that faces the end, ``F`` bases (more if its
    motif needs them)."""
    bld = _Builder()
    if (side == "left") == (g.strand == 1):     # the gene's 5' end faces the end: its U
        ell = min(max(F, g.cas.layout.need_motif(g.cas.m)), g.cas.layout.u_len)
        bld.add(g, _u_part(g, ell), g.strand == -1)
    else:
        ell = min(F, g.cas.layout.d_len)
        bld.add(g, _d_part(g, ell), g.strand == -1)
    return Junction(f"{side}_end", "end", ell, ell, 0, bld.seq(), bld.elements)


def apply_edits(cds: str, edit: dict[str, str]) -> tuple[str, str]:
    """The CDS with its start and stop codons replaced as ``edit`` says, and the ``cds_edit`` text."""
    notes = []
    if "start" in edit and edit["start"] != cds[:3]:
        notes.append(f"start:{cds[:3]}>{edit['start']}")
        cds = edit["start"] + cds[3:]
    if "stop" in edit and edit["stop"] != cds[-3:]:
        notes.append(f"stop:{cds[-3:]}>{edit['stop']}")
        cds = cds[:-3] + edit["stop"]
    return cds, ";".join(notes)


# --- replicon assembly --------------------------------------------------------------------------------------

@dataclass
class Assembled:
    seq: str
    gene_spans: list[tuple[int, int]]          # per gene, start and end in the cut frame (end > N across origin)
    junction_spans: list[tuple[int, int]]      # per junction (and end flank), start and end of its bases
    origin_offset: int
    edits: list[str]                           # per gene, its cds_edit text
    gene_seqs: list[str]                       # per gene, its + strand bases as laid out


def assemble(genes: list[GeneIn], juncs: list[Junction], circular: bool, origin: int = 0,
             flanks: tuple[Junction, Junction] | None = None) -> Assembled:
    """Lay genes and junctions out (circular: junction j follows gene j, the last one closing the circle;
    linear: n - 1 junctions between the genes and the two end flanks), then cut a circle at ``origin`` bases
    after the first gene's start. Checks that overlapping bases agree."""
    n = len(genes)
    edit: dict[str, dict] = {}
    for j in juncs:
        for c, e in j.edits.items():
            edit.setdefault(c, {}).update(e)
    plus, notes = [], []
    for g in genes:
        e = edit.get(g.copy_label, {})
        if e.get("start", "ATG") not in _STARTS_OK or e.get("stop", "TAA") not in STOPS:
            raise AssertionError(f"{g.copy_label}: an overlap edit {e} is not protein-neutral")
        cds, note = apply_edits(g.cds, e)
        plus.append(cds if g.strand == 1 else revcomp(cds))
        notes.append(note)
    pos = 0
    pieces: list[tuple[int, str]] = []
    gspan, jspan = [], []
    if not circular:
        if flanks is None:
            raise ValueError("a linear replicon needs its two end flanks")
        pieces.append((0, flanks[0].seq))
        jspan.append((0, flanks[0].length))
        pos = flanks[0].length
    n_j = n if circular else n - 1
    if len(juncs) != n_j:
        raise ValueError(f"{n} genes need {n_j} junctions, got {len(juncs)}")
    for i in range(n):
        gspan.append((pos, pos + len(plus[i])))
        pieces.append((pos, plus[i]))
        end = pos + len(plus[i])
        if i < n_j:
            j = juncs[i]
            if j.overlap:
                jspan.append((end, end))
                pos = end - j.overlap
            else:
                jspan.append((end, end + j.length))
                pieces.append((end, j.seq))
                pos = end + j.length
        else:
            pos = end
    if not circular:
        pieces.append((pos, flanks[1].seq))
        jspan.append((pos, pos + flanks[1].length))
        pos += flanks[1].length
    N = pos
    for i in range(n_j):
        o_i = juncs[i].overlap
        if o_i and plus[i][-o_i:] != plus[(i + 1) % n][:o_i]:
            raise AssertionError(f"junction {i}: the overlapping bases of its genes disagree")
    if sum(len(s) for _, s in pieces) - sum(j.overlap for j in juncs) != N:
        raise AssertionError("genes and junctions do not tile the replicon")
    buf = bytearray(N)
    for at, s in pieces:
        b = s.encode("ascii")
        a = at % N
        if a + len(b) <= N:
            buf[a:a + len(b)] = b
        else:
            buf[a:] = b[:N - a]
            buf[:len(b) - (N - a)] = b[N - a:]
    o = int(origin) % N if circular else 0
    seq = (bytes(buf[o:]) + bytes(buf[:o])).decode("ascii")

    def cut(s: int, e: int) -> tuple[int, int]:
        s2 = (s - o) % N
        return s2, s2 + (e - s)
    return Assembled(seq, [cut(*x) for x in gspan], [cut(*x) for x in jspan], o, notes, plus)


def operons(kinds: list[str], n: int, circular: bool) -> list[int]:
    """Per gene, the position of its operon's first gene: genes are joined across a ``codir_tight``
    junction (junction j follows gene j)."""
    tight = [k == "codir_tight" for k in kinds]
    first = list(range(n))
    if circular and n and all(tight[:n]):
        return [0] * n
    start = 0
    if circular:
        start = next(i for i in range(n) if not tight[i - 1])   # a gene whose left junction breaks
    for k in range(n):
        i = (start + k) % n
        left = (i - 1) % n
        joined = tight[left] if (circular or i > 0) else False
        first[i] = first[left] if joined and k > 0 else i
    return first


# --- I/O ----------------------------------------------------------------------------------------------------

def _tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class RepliconRow(NamedTuple):
    genome_id: str
    replicon: int
    record: str
    kind: str
    topology: str
    length: int
    n_genes: int
    origin_offset: int
    gc: float
    contig_uid: str


class GeneRow(NamedTuple):
    genome_id: str
    replicon: int
    position: int
    copy_label: str
    family_id: str
    start: int
    end: int
    strand: int
    origin_spanning: bool
    operon_id: str
    start_codon: str
    stop_codon: str
    cds_edit: str


class JunctionRow(NamedTuple):
    genome_id: str
    replicon: int
    left_position: int
    right_position: int
    left_copy: str
    right_copy: str
    junction: str
    rule: str
    target_length: int
    length: int
    overlap_bp: int
    start: int
    end: int
    igr_uid: str
    gc: float


class ElementRow(NamedTuple):
    genome_id: str
    replicon: int
    left_position: int
    kind: str
    start: int
    end: int
    strand: int
    copy_label: str
    rna_id: str


class RnaRow(NamedTuple):
    rna_id: str
    link: str
    model_acc: str
    model_name: str
    source: str
    emission_id: str
    exp: float
    local_mode: bool
    attached_copy: str
    family_id: str
    genome_id: str
    replicon: int
    start: int
    end: int
    strand: int
    length: int
    dist_to_start: int
    n_subst_from_emission: int
    n_compensatory: int
    structure: str


class EmissionTsv(NamedTuple):
    emission_id: str
    link: str
    model_acc: str
    clade_node: int
    height: float
    sequence: str
    structure: str


class LinkRow(NamedTuple):
    link: str
    model_index: int
    model_acc: str
    model_name: str
    family: int
    family_id: str
    n_clades: int
    clade_height: float


class ModelRow(NamedTuple):
    index: int
    accession: str
    name: str
    source: str
    length: int
    ga: float
    n_pairs: int
    structure: str


class LineageRow(NamedTuple):
    lineage: str
    gc3: float
    gc_igr: float


_TABLES = {"replicons.tsv": RepliconRow, "genes.tsv": GeneRow, "junctions.tsv": JunctionRow,
           "elements.tsv": ElementRow, "rna.tsv": RnaRow, "emissions.tsv": EmissionTsv, "links.tsv": LinkRow,
           "models.tsv": ModelRow, "lineages.tsv": LineageRow}
_SMALL = ("replicons.tsv", "emissions.tsv", "links.tsv", "models.tsv", "lineages.tsv")

_CAS_FIELDS = ("copy_label", "family", "cds", "u", "d", "z", "layout", "m", "link", "emission", "n_subst",
               "n_comp", "structure")


def _cas_line(family: int, cds: str, c: CopyCassette) -> str:
    lay = c.layout
    return "\t".join([c.copy_label, str(family), cds, c.u, c.d, ",".join(repr(v) for v in c.z),
                      ",".join(str(getattr(lay, f)) for f in Layout.__dataclass_fields__), str(c.m), c.link,
                      c.emission, str(c.n_subst), str(c.n_comp), c.structure]) + "\n"


def _cas_parse(line: str) -> tuple[int, str, CopyCassette]:
    f = line.rstrip("\n").split("\t")
    lay = Layout(*(int(v) for v in f[6].split(",")))
    z = tuple(float(v) for v in f[5].split(","))
    return int(f[1]), f[2], CopyCassette(f[0], f[0].split("_")[0], f[3], f[4], z, lay, int(f[7]), f[8], f[9],
                                         int(f[10]), int(f[11]), f[12])


def check_config(cfg: dict, motifs: dict) -> None:
    J = cfg["junctions"]
    if not 0.0 <= float(J["p_tight"]) <= 1.0:
        raise ValueError(f"p_tight must lie in [0, 1], got {J['p_tight']}")
    ov = J["tight_overlap"]
    if not (0.0 <= float(ov["minus4"]) and 0.0 <= float(ov["minus1"]) and float(ov["minus4"]) + float(ov["minus1"]) < 1.0):
        raise ValueError(f"tight_overlap: two probabilities summing below 1, got {ov}")
    cas = cfg["cassettes"]
    for kind in ("tight", "loose", "divergent", "convergent"):
        d = J[kind]
        if not (0 <= int(d["min"]) <= int(d["max"]) and float(d["scale"]) > 0 and float(d["sigma"]) >= 0):
            raise ValueError(f"junctions.{kind}: need 0 <= min <= max, scale > 0, sigma >= 0, got {d}")
    if int(J["loose"]["max"]) > int(cas["upstream_length"]) or int(J["tight"]["max"]) > int(cas["upstream_length"]):
        raise ValueError("junction lengths taken from one U cannot exceed cassettes.upstream_length")
    if motifs["source"] not in ("synthetic", "rfam"):
        raise ValueError(f"rna_motifs.source is 'synthetic' or 'rfam', got {motifs['source']!r}")


def _sha(path: Path) -> str:
    return sequences._sha(path)


def _fasta_records(path: Path):
    name, chunks = None, []
    with open(path, encoding="ascii") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def _models(motifs: dict, cas: dict, seed: int, rfam_cm, scratch: Path) -> tuple[list[rna.MotifModel], dict]:
    if motifs["source"] == "synthetic":
        return rna.synthetic_models(motifs, cas["pair_weights"], seed), {"source": "synthetic", "file": None,
                                                                        "sha256": None, "release": None}
    if rfam_cm is None:
        raise ValueError("rna_motifs.source 'rfam' needs the Rfam CM library (rfam_cm=)")
    accs = list(motifs["rfam_models"])
    got = rna.extract_rfam(Path(rfam_cm), accs, scratch / "models")
    out = []
    for k, acc in enumerate(accs):
        f = scratch / "models" / f"{acc}.cm"
        out.append(rna.MotifModel(k, acc, got[acc]["name"], "rfam", got[acc]["clen"], got[acc]["ga"], f))
    return out, {"source": "rfam", "file": Path(rfam_cm).name, "sha256": _sha(Path(rfam_cm)),
                 "release": str(motifs["rfam_release"])}


def _links(world: Path, motifs: dict, models: list[rna.MotifModel], seed: int, have: set[int]) -> dict[int, Link]:
    fams = _tsv(world / "backbone" / "families.tsv")
    kinds = set(motifs["family_kinds"])
    elig = sorted(int(r["family"]) for r in fams if r["kind"] in kinds
                  and int(r["n_genomes_present"]) >= int(motifs["min_genomes"]) and int(r["family"]) in have)
    per = int(motifs["families_per_model"])
    need = per * len(models)
    if len(elig) < need:
        raise ValueError(f"{need} motif links need as many eligible families, the world has {len(elig)}")
    rng = seeds.generator(seed, "S5/motif/links")
    pick = [elig[int(i)] for i in rng.permutation(len(elig))[:need]]
    return {f: Link(f"L{k:02d}", models[k // per], f, f"F{f:06d}") for k, f in enumerate(pick)}


@dataclass
class RepWorld:
    seed: int
    replicons: list[RepliconRow]
    emissions: list[EmissionTsv]
    links: list[LinkRow]
    models: list[ModelRow]
    lineages: list[LineageRow]
    digest: str
    out: Path
    meta: dict = field(default_factory=dict, repr=False)


def build(world, *, rfam_cm=None) -> RepWorld:
    """Write ``world/replicons/`` from ``world/backbone/`` and ``world/sequences/`` (module docstring)."""
    world = Path(world)
    bb = json.loads((world / "backbone" / "backbone.json").read_text(encoding="utf-8"))
    if bb["format"] != backbone.FORMAT:
        raise ValueError(f"{world}: backbone format {bb['format']!r}, expected {backbone.FORMAT!r}")
    sq = json.loads((world / "sequences" / "sequences.json").read_text(encoding="utf-8"))
    if sq["format"] != sequences.FORMAT:
        raise ValueError(f"{world}: sequences format {sq['format']!r}, expected {sequences.FORMAT!r}")
    preset = bb["preset"]
    cfg, motifs = preset["intergenic"], preset["rna_motifs"]
    check_config(cfg, motifs)
    seed = int(bb["seed"])
    out = world / "replicons"
    if out.exists():
        raise FileExistsError(f"{out} exists; replicons are written once")
    tmp = world / "replicons.partial"
    if tmp.exists():
        shutil.rmtree(tmp)
    scratch = tmp / "scratch"
    (tmp / "fasta").mkdir(parents=True)
    scratch.mkdir()
    lin = lineages(world, cfg["composition"])
    kappa = float(sq["sequences_preset"]["evolution"]["kappa_syn"])
    cfg = dict(cfg, _motif_rate_factor=float(motifs["rate_factor"]))
    fam_rows = {int(r["family"]): r for r in _tsv(world / "sequences" / "families.tsv")}
    models, src = _models(motifs, cfg["cassettes"], seed, rfam_cm, scratch)
    links = _links(world, motifs, models, seed, set(fam_rows))
    tx = bb["taxonomy"]
    h = float(motifs["clade_height"]) * (float(tx["present"]) - float(tx["crown"]))
    genome_of = {r["species_node"]: r["genome_id"] for r in _tsv(world / "backbone" / "genomes.tsv")}
    versions = {"numpy": np.__version__}
    if motifs["source"] == "rfam":
        versions["infernal"] = rna.infernal_version()

    # stage 1: cassettes, one family at a time, into per-genome working files
    buffers: dict[str, list[str]] = {}
    n_buf = [0]

    def flush() -> None:
        for gid, lines in sorted(buffers.items()):
            with open(scratch / f"{gid}.cas", "a", encoding="ascii") as fh:
                fh.writelines(lines)
        buffers.clear()
        n_buf[0] = 0

    emissions: list[EmissionRow] = []
    cds_iter = _fasta_records(world / "sequences" / "cds.fna")
    for family, nodes in sequences._stream_trees(world / "backbone" / "gene_tree_nodes.tsv"):
        if not any(n.extant for n in nodes):
            continue
        fr = fam_rows[family]
        link = links.get(family)
        emit = None
        if link is not None:
            if link.model.source == "rfam":
                def emit(n, link=link):
                    return rna.cmemit(link.model.cm_file, n, seeds.tool_seed(seed, f"S5/emit/{link.link}"),
                                      exp=float(motifs["exp"]), local=bool(motifs["local"]))
            else:
                def emit(n, link=link):
                    return rna.synthetic_emit(link.model, n, seeds.generator(seed, f"S5/emit/{link.link}"))
        cas = evolve_family(fr["family_id"], nodes, float(fr["rate_multiplier"]), lin, cfg, kappa, seed, link,
                            clade_height=h, emit=emit, emissions_out=emissions)
        for c in cas:
            name, cds = next(cds_iter)
            if name != c.copy_label:
                raise ValueError(f"cds.fna holds {name} where {c.copy_label} was expected (family order)")
            buffers.setdefault(genome_of[c.lineage], []).append(_cas_line(family, cds, c))
            n_buf[0] += 1
        if n_buf[0] >= _FLUSH_LINES:
            flush()
    flush()
    if next(cds_iter, None) is not None:
        raise ValueError("cds.fna holds more records than the gene trees' extant copies")

    # stage 2: one genome at a time
    planted = {r["name"] for r in _tsv(world / "backbone" / "groups.tsv") if r["kind"] == "planted"}
    reps = {}
    for r in _tsv(world / "backbone" / "replicons.tsv"):
        reps[(r["genome_id"], int(r["replicon"]))] = r
    J = cfg["junctions"]
    writers = {name: _Rows(tmp / name, cls) for name, cls in _TABLES.items()
               if name in ("replicons.tsv", "genes.tsv", "junctions.tsv", "elements.tsv", "rna.tsv")}
    link_by_id = {lk.link: lk for lk in links.values()}
    counts = {"genomes": 0, "replicons": 0, "genes": 0, "junctions": 0, "rna": 0}
    for gid, rows in _genome_rows(world / "backbone" / "genes.tsv"):
        cas_of: dict[str, tuple[int, str, CopyCassette]] = {}
        path = scratch / f"{gid}.cas"
        with open(path, encoding="ascii") as fh:
            for line in fh:
                fam, cds, c = _cas_parse(line)
                cas_of[c.copy_label] = (fam, cds, c)
        fa = sequences._Fasta(tmp / "fasta" / f"{gid}.fna")
        by_rep: dict[int, list[dict]] = {}
        for r in rows:
            by_rep.setdefault(int(r["replicon"]), []).append(r)
        for rep in sorted(by_rep):
            rr = sorted(by_rep[rep], key=lambda r: int(r["position"]))
            genes = []
            for r in rr:
                fam, cds, c = cas_of[r["copy_label"]]
                key = f"{r['group']}/{r['unit']}" if r["group"] in planted else ""
                genes.append(GeneIn(r["copy_label"], fam, int(r["strand"]), key, cds, c))
            topo = reps[(gid, rep)]["topology"]
            circular = topo == "circular"
            n = len(genes)
            if circular:
                juncs = [junction(genes[i], genes[(i + 1) % n], J, seed, same=n == 1) for i in range(n)]
                flanks = None
                N_guess = sum(len(g.cds) for g in genes) + sum(j.length - j.overlap for j in juncs)
                origin = int(seeds.generator(seed, f"S5/origin/{gid}/{rep}").integers(N_guess))
            else:
                juncs = [junction(genes[i], genes[i + 1], J, seed) for i in range(n - 1)]
                F = int(J["linear_end_flank"])
                flanks = (end_flank(genes[0], F, "left"), end_flank(genes[-1], F, "right"))
                origin = 0
            asm = assemble(genes, juncs, circular, origin, flanks)
            N = len(asm.seq)
            record = f"{gid}_r{rep}"
            fa.write(record, asm.seq)
            seqb = asm.seq.encode("ascii")
            gc = q6((asm.seq.count("G") + asm.seq.count("C")) / N)
            writers["replicons.tsv"].row(RepliconRow(gid, rep, record, reps[(gid, rep)]["kind"], topo, N, n,
                                                     asm.origin_offset, gc, ids.contig_uid(seqb).hex()))
            first = operons([j.kind for j in juncs], n, circular)
            for i, g in enumerate(genes):
                s, e = asm.gene_spans[i]
                cds_real = asm.gene_seqs[i] if g.strand == 1 else revcomp(asm.gene_seqs[i])
                writers["genes.tsv"].row(GeneRow(gid, rep, int(rr[i]["position"]), g.copy_label,
                                                 f"F{g.family:06d}", s, e, g.strand, e > N,
                                                 f"{gid}_r{rep}_o{int(rr[first[i]]['position'])}",
                                                 cds_real[:3], cds_real[-3:], asm.edits[i]))
            spans = asm.junction_spans
            all_j = ([flanks[0]] + juncs + [flanks[1]]) if flanks else juncs
            for k, j in enumerate(all_j):
                s, e = spans[k]
                if flanks:
                    li, ri = k - 1, k
                    li = -1 if k == 0 else li
                    ri = -1 if k == len(all_j) - 1 else ri
                else:
                    li, ri = k, (k + 1) % n
                lp = int(rr[li]["position"]) if li >= 0 else -1
                rp = int(rr[ri]["position"]) if ri >= 0 else -1
                igr = j.seq
                writers["junctions.tsv"].row(JunctionRow(
                    gid, rep, lp, rp, genes[li].copy_label if li >= 0 else "",
                    genes[ri].copy_label if ri >= 0 else "", j.kind, j.rule, j.target, j.length, j.overlap, s, e,
                    ids.contig_uid(igr.encode("ascii")).hex() if igr else "",
                    q6((igr.count("G") + igr.count("C")) / len(igr)) if igr else float("nan")))
                for kind, es, ee, strand, owner in j.elements:
                    a = s + es
                    if a >= N:
                        a -= N
                    rid = f"{owner}.rna" if kind == "rna_motif" else ""
                    writers["elements.tsv"].row(ElementRow(gid, rep, lp, kind, a, a + (ee - es), strand, owner, rid))
                    if kind == "rna_motif":
                        fam, _, c = cas_of[owner]
                        lk = link_by_id[c.link]
                        writers["rna.tsv"].row(RnaRow(
                            rid, lk.link, lk.model.accession, lk.model.name, lk.model.source, c.emission,
                            float(motifs["exp"]) if lk.model.source == "rfam" else float("nan"),
                            bool(motifs["local"]) and lk.model.source == "rfam", owner, f"F{fam:06d}", gid, rep,
                            a, a + (ee - es), strand, ee - es,
                            c.layout.utr_b + len(RBS) + c.layout.sp2, c.n_subst, c.n_comp, c.structure))
                        counts["rna"] += 1
            counts["replicons"] += 1
            counts["genes"] += n
            counts["junctions"] += len(all_j)
        fa.close()
        path.unlink()
        counts["genomes"] += 1
    for w in writers.values():
        w.close()
    shutil.rmtree(scratch)
    motif_carriers = sum(int(r["n_extant_copies"]) for r in _tsv(world / "backbone" / "families.tsv")
                         if int(r["family"]) in links)
    if counts["rna"] != motif_carriers:
        raise AssertionError(f"{counts['rna']} motif instances for {motif_carriers} copies of linked families")

    t = _Rows(tmp / "emissions.tsv", EmissionTsv)
    for e in sorted(emissions, key=lambda e: (e.link, int(e.emission_id.split(".")[1]))):
        t.row(EmissionTsv(e.emission_id, e.link, e.model_acc, e.clade_node, e.height, e.emission.seq,
                          rna.dot_bracket(len(e.emission.seq), e.emission.pairs)))
    t.close()
    t = _Rows(tmp / "links.tsv", LinkRow)
    for lk in sorted(links.values(), key=lambda x: x.link):
        t.row(LinkRow(lk.link, lk.model.index, lk.model.accession, lk.model.name, lk.family, lk.family_id,
                      lk.n_clades, q6(h)))
    t.close()
    t = _Rows(tmp / "models.tsv", ModelRow)
    for m in models:
        t.row(ModelRow(m.index, m.accession, m.name, m.source, m.length, m.ga, len(m.pairs), m.structure))
    t.close()
    t = _Rows(tmp / "lineages.tsv", LineageRow)
    for lab, i in sorted(lin.index.items(), key=lambda kv: kv[1]):
        t.row(LineageRow(lab, float(lin.gc3[i]), float(lin.gc_igr[i])))
    t.close()
    files = {str(f.relative_to(tmp)): _sha(f) for f in sorted(tmp.rglob("*")) if f.is_file()}
    digest = blake3("".join(f"{k}\t{v}\n" for k, v in sorted(files.items())).encode()).hexdigest(length=16)
    pre = {"intergenic": preset["intergenic"], "rna_motifs": preset["rna_motifs"]}
    doc = {"format": FORMAT, "seed": seed, "preset": pre, "params_digest": backbone.params_digest(pre),
           "backbone_digest": bb["digest"], "sequences_digest": sq["digest"], "motif_source": src,
           "versions": versions, "counts": counts,
           "stream_paths": ["S5/family/<family_id>", "S5/motif/panel", "S5/motif/links", "S5/emit/<link>",
                            "S5/origin/<genome_id>/<replicon>", "S5/tight/<up>/<down> (BLAKE3 hash)"],
           "files": files, "digest": digest}
    (tmp / "replicons.json").write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    tmp.rename(out)
    return read(world)


def _genome_rows(path: Path):
    """``backbone/genes.tsv`` (sorted by genome) one genome at a time."""
    gid, rows = None, []
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["genome_id"] != gid and rows:
                yield gid, rows
                rows = []
            gid = r["genome_id"]
            rows.append(r)
    if rows:
        yield gid, rows


def read(world) -> RepWorld:
    """Load ``world/replicons/``'s small tables (the per-gene, per-junction, per-element and per-motif tables
    grow with the world and are read by the caller as needed)."""
    world = Path(world)
    rdir = world / "replicons"
    doc = json.loads((rdir / "replicons.json").read_text(encoding="utf-8"))
    if doc["format"] != FORMAT:
        raise ValueError(f"{rdir}: format {doc['format']!r}, expected {FORMAT!r}")
    t = {}
    for name in _SMALL:
        cls = _TABLES[name]
        conv = backbone._parser(cls)
        with open(rdir / name, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            if tuple(header) != cls._fields:
                raise ValueError(f"{rdir / name}: header {header}, expected {list(cls._fields)}")
            t[name] = [cls(*(c(v) for c, v in zip(conv, row))) for row in reader]
    return RepWorld(doc["seed"], t["replicons.tsv"], t["emissions.tsv"], t["links.tsv"], t["models.tsv"],
                    t["lineages.tsv"], doc["digest"], world, doc)


__all__ = ["FORMAT", "Layout", "GeneIn", "Junction", "Assembled", "RepWorld", "build", "read", "junction",
           "end_flank", "assemble", "overlap_edits", "apply_edits", "operons", "length_of", "gc_igr",
           "tight_u", "revcomp", "evolve_family", "draw_layout", "site_model", "element_columns"]
