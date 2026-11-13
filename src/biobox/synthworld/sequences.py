"""synthworld family design and protein-coding sequences (planning/tests/synthworld_spec.md §2.1 S3-S4, §2.3;
card F-02 session 3). Reads a built backbone (`backbone`, format ``synthworld-backbone/2``) and never ZOMBI2.

**Family design (S3).** Every family with an extant copy gets, from its own stream ``S3/family/<family_id>``
and the world's HMM panel (`profiles`): a dark flag (``architecture.dark_fraction``; a dark family draws
its domains from the profiles held out of the scan panel), a domain count (``domain_count_weights`` for
1, 2, 3, ...), a promiscuous-domain trap prepended with probability ``promiscuous_fraction`` (not to dark
families), its domains (uniform over the pool, with replacement), N-tail, linker and C-tail lengths
(geometric on 0, 1, ... with the preset means), a difficulty class and a target identity (uniform within the
class's range; classes taken in name order, weights ``evolution.difficulty``), and its root protein:
the initiator, the N-tail, each domain drawn by `profiles.sample_domain` (a walk of the profile's core model),
the linkers, the C-tail; tail, linker and insert-state residues from the world background. Each root site
gets a discrete-gamma rate category (``gamma_alpha``, ``gamma_categories``; Yang's category means).

**Evolution (S4)**, from the family's stream ``S4/family/<family_id>``, down the family's complete gene tree
(``backbone/gene_tree_nodes.tsv``) pruned to the nodes ancestral to an extant copy (a branch with no extant
descendant cannot change an extant sequence, so pruning leaves the extant copies' distribution unchanged).
On a branch of length t in lineage b (the branch above a gene-tree node lies in that node's lineage):
- site s has a substitution event with probability 1 - exp(-mu_f r_s c_b t); the new residue is drawn from
  the site's equilibrium, its match-state emission row (domain match sites) or the background (all other
  sites): F81 with site-specific frequencies. The new residue's codon is drawn from lineage b's usage;
- every codon (the initiation and stop codons included) has a synonymous event with probability
  1 - exp(-kappa_syn mu_f c_b t) and is redrawn from lineage b's usage for its residue;
- each tail and linker segment of current length n has Poisson(indel_rate mu_f c_b t (n + 1)) indel events,
  insertion or deletion with probability 1/2, length uniform on 1..``indel_max_length``, inserted residues
  from the background with new rate categories; domains have no indels (§2.3 v1).
  Within a branch, substitutions are applied before indels.
mu_f solves the closed form of §2.3 for the family's target: the mean over pairs of extant copies (all pairs,
or ``target_pairs`` drawn) of mean_s [h_s + (1 - h_s) exp(-mu r_s T)], T the path length in clock-scaled time
and h_s the site's sum of squared frequencies, over the root's tail, linker and domain sites; a family with
one copy targets its identity to the root. A target below the saturation identity mean_s h_s is raised to
it plus `SATURATION_MARGIN`; mu_f is bounded by ``rate_bounds``.

**Lineages.** On the complete species tree (`backbone.species_branches`, the stem included), in preorder,
stream ``S4/lineages`` draws per branch b: logit(GC3_b) = logit(GC3 of the parent branch, or ``gc3_root``
at the stem) + N(0, gc3_sigma^2 t_b), and the relaxed-clock multiplier c_b = exp(N(-sigma^2 / 2, sigma^2)),
sigma = ``clock_sigma`` (uncorrelated, mean 1). Each branch's GC3 is constant along it and sets its codon
usage (`codons.usage_table`), so a transferred gene arrives with donor codons and ameliorates by
synonymous and non-synonymous events in its new lineage.

Shape on disk (``build(world)`` writes ``world/sequences/``; ``read(world)`` loads it)::

    cds.fna, cds.fna.fai         coding sequences, initiation codon to stop codon inclusive, 5'->3' on the
                                 coding strand, one record per extant copy named by its copy_label
    proteins.faa, .fai           their translations under code 11 (initiator M, no stop character)
    roots.faa, .fai              each family's root protein, named by family_id
    families.tsv   family family_id dark promiscuous root_source n_domains domain_hmms n_tail c_tail linkers
                   root_length difficulty target_identity rate_multiplier expected_identity
                   saturation_identity n_target_pairs realized_mean_identity realized_min_identity
                   n_extant_copies
    domains.tsv    family family_id index profile accession name aa_start aa_end n_match n_insert n_delete
                   promiscuous in_scan_panel
    genes.tsv      copy_label family family_id genome_id species_node genetic_code cds_length
                   protein_length start_codon stop_codon protein_uid gc gc3 path_time subst_per_site
                   realized_subst_per_site root_identity hgt_on_lineage n_indels domain_arch
    lineages.tsv   lineage parent branch_length clock gc3
    panel.tsv      index name accession source match_length in_scan_panel promiscuous ga_seq ga_dom
                   saturation
    sequences.json format, seed, the ``sequences`` preset and its digest, the backbone digest, the panel
                   source (file name, sha256, release), versions, sha256 of every file above, digest

FASTA is uncompressed, 80 columns, with a samtools-style ``.fai`` (name, length, offset, 80, 81). Records
and rows are in family order (the order of ``backbone/gene_tree_nodes.tsv``) and, within a family, by
copy_label; ``domains`` by family and index. A gene is placed by joining ``backbone/genes.tsv`` on
copy_label (genome, replicon, position, strand; `iter_genome` does it one genome at a time); a CDS on
strand -1 is laid out reverse-complemented. Coordinates are 0-based half-open: ``domains.aa_start/aa_end``
on the root protein, ``genes.domain_arch`` (``accession:start-end;...``, in domain order) on the copy's
protein. ``root_identity`` and ``realized_subst_per_site`` count the root's tail, linker and domain sites
still present in the copy (inserted residues are not counted); ``subst_per_site`` is the expected number of
substitution events per site from the family root, mu_f x the clock-scaled path length; ``path_time`` the
same path in unscaled time; ``hgt_on_lineage`` the transfers on the root-to-copy path that the path enters
on the recipient side. Identities over pairs (``realized_*_identity``) are over the root's sites present in
both copies, over all pairs when a family has at most ``exact_identity_max_copies`` copies (else
that many pairs drawn), empty cells ``nan`` below two copies. Floats are written at 6 significant digits.
The genetic code is 11 everywhere (codes 4 and 25 are S8's). ``digest`` is BLAKE3-128 over the sorted
(path, sha256) list of every file but sequences.json.
"""

from __future__ import annotations

import csv
import functools
import hashlib
import json
import math
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

import numpy as np
from blake3 import blake3

from biobox.synthworld import backbone, codons, profiles, seeds

FORMAT = "synthworld-sequences/1"
SITE_START, SITE_TAIL, SITE_MATCH, SITE_INSERT, SITE_STOP = range(5)
SATURATION_MARGIN = 0.01
_INDEL_KINDS = ("ntail", "linker", "ctail")


def q6(x) -> float:
    """A float at 6 significant digits (outputs and model constants derived from floats), so that a last-bit
    difference in a libm or SIMD exp cannot reach a written value."""
    return float(f"{float(x):.6g}")


# --- discrete gamma (§2.3 "r_s a Γ rate category") ------------------------------------------------------

def _gammap(a: float, x: float) -> float:
    """Regularized lower incomplete gamma P(a, x): series below a + 1, Lentz continued fraction above."""
    if x <= 0.0:
        return 0.0
    log_front = -x + a * math.log(x) - math.lgamma(a)
    if x < a + 1.0:
        ap, term = a, 1.0 / a
        total = term
        for _ in range(100_000):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * 1e-16:
                break
        return min(1.0, total * math.exp(log_front))
    tiny = 1e-300
    b = x + 1.0 - a
    c, d = 1.0 / tiny, 1.0 / b
    h = d
    for i in range(1, 100_000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        d = tiny if abs(d) < tiny else d
        c = b + an / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-16:
            break
    return max(0.0, 1.0 - math.exp(log_front) * h)


def gamma_rates(alpha: float, k: int) -> np.ndarray:
    """Yang's (1994) discrete gamma: the means of the k equal-probability categories of Gamma(alpha, rate
    alpha), mean 1. Category j's mean is k [P(alpha+1, alpha q_j) - P(alpha+1, alpha q_{j-1})] with q_j the
    j/k quantile."""
    return np.array(_gamma_rates(float(alpha), int(k)))


@functools.lru_cache(maxsize=16)
def _gamma_rates(alpha: float, k: int) -> tuple[float, ...]:
    a = float(alpha)
    if a <= 0 or int(k) < 1:
        raise ValueError(f"gamma_alpha must be > 0 and gamma_categories >= 1, got {alpha}, {k}")
    if k == 1:
        return (1.0,)
    qs = [0.0]
    for j in range(1, k):
        lo, hi = 0.0, 1.0
        while _gammap(a, a * hi) < j / k:
            hi *= 2.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if _gammap(a, a * mid) < j / k:
                lo = mid
            else:
                hi = mid
        qs.append(0.5 * (lo + hi))
    cdf1 = [_gammap(a + 1.0, a * q) for q in qs] + [1.0]
    return tuple(q6(k * (cdf1[j + 1] - cdf1[j])) for j in range(k))


# --- divergence control in closed form (§2.3) ------------------------------------------------------------

def _collapse(rates: np.ndarray, h: np.ndarray):
    """Sites grouped by rate: (rate values, weight, mean h) per group (expected identity is linear in h)."""
    vals, inv = np.unique(rates, return_inverse=True)
    w = np.bincount(inv, minlength=len(vals)) / len(rates)
    hbar = np.bincount(inv, weights=h, minlength=len(vals)) / np.maximum(1, np.bincount(inv, minlength=len(vals)))
    return vals, w, hbar


def _expected(mu: float, vals, w, hbar, T: np.ndarray) -> float:
    e = np.exp(-mu * np.outer(vals, T)).mean(axis=1)
    return float((w * (hbar + (1.0 - hbar) * e)).sum())


def solve_rate(target: float, rates: np.ndarray, h: np.ndarray, T, *, bounds) -> tuple[float, float, float]:
    """mu such that the mean over path lengths ``T`` and sites of h_s + (1 - h_s) exp(-mu r_s T) equals
    ``target`` (raised to the saturation identity mean(h) + `SATURATION_MARGIN` when below it), within
    ``bounds``. Returns (mu at 6 significant digits, the expected identity at that mu, the saturation)."""
    T = np.asarray(T, dtype=np.float64)
    vals, w, hbar = _collapse(np.asarray(rates, dtype=np.float64), np.asarray(h, dtype=np.float64))
    sat = float((w * hbar).sum())
    goal = max(float(target), sat + SATURATION_MARGIN)
    lo, hi = (float(b) for b in bounds)
    if goal >= _expected(lo, vals, w, hbar, T):
        mu = lo
    elif goal <= _expected(hi, vals, w, hbar, T):
        mu = hi
    else:
        a, b = math.log(lo), math.log(hi)
        for _ in range(64):                # the bracket shrinks to below 1e-17 in log mu, far under q6
            m = 0.5 * (a + b)
            if _expected(math.exp(m), vals, w, hbar, T) > goal:
                a = m
            else:
                b = m
        mu = math.exp(0.5 * (a + b))
    mu = q6(mu)
    return mu, _expected(mu, vals, w, hbar, T), sat


# --- family design (S3) -----------------------------------------------------------------------------------

@dataclass
class Domain:
    profile: int
    accession: str
    name: str
    start: int        # root protein coordinates, 0-based half-open
    end: int
    n_match: int
    n_insert: int
    n_delete: int
    promiscuous: bool
    in_scan_panel: bool


@dataclass
class Design:
    family: int
    family_id: str
    dark: bool
    promiscuous: bool
    n_domains: int                 # drawn count, the promiscuous trap not included
    domains: list[Domain]
    n_tail: int
    c_tail: int
    linkers: list[int]
    difficulty: str
    target_identity: float
    aa: np.ndarray                 # (L,) residue codes: START, 0..19, STOP
    kind: np.ndarray               # (L,) SITE_*
    segment: np.ndarray            # (L,) segment index, increasing along the protein
    seg_kind: list[str]            # start, ntail, domain, linker, ctail, stop
    pi_table: np.ndarray           # (P, 20) equilibrium rows; row 0 the background
    pi_idx: np.ndarray             # (L,) row of pi_table per site
    rates: np.ndarray              # (L,) rate category means; 0 at START and STOP
    rate_values: np.ndarray        # (k,) the category means

    @property
    def root_length(self) -> int:
        """Protein length (the stop codon has no residue)."""
        return int(len(self.aa) - 1)


def _geometric(rng, mean: float) -> int:
    m = float(mean)
    return int(rng.geometric(1.0 / (m + 1.0))) - 1 if m > 0 else 0


def design_family(cfg: dict, seed: int, family: int, panel: list[profiles.Profile]) -> Design:
    """The family's architecture, target and root protein (module docstring), from ``S3/family/<id>``."""
    fid = f"F{int(family):06d}"
    arch, ev = cfg["architecture"], cfg["evolution"]
    if arch["root_source"] != "hmm_sample":
        raise NotImplementedError(f"root_source {arch['root_source']!r}: only 'hmm_sample' is implemented "
                                  f"(real_seed_member needs Pfam's seed alignments)")
    rng = seeds.generator(seed, f"S3/family/{fid}")
    dark = bool(rng.random() < float(arch["dark_fraction"]))
    w = np.asarray(arch["domain_count_weights"], dtype=np.float64)
    cw = np.cumsum(w / w.sum())
    n_dom = 1 + min(int(np.searchsorted(cw, rng.random(), side="right")), len(w) - 1)
    promiscuous = bool(rng.random() < float(arch["promiscuous_fraction"])) and not dark
    scan = [p for p in panel if p.in_scan_panel and not p.promiscuous]
    held = [p for p in panel if not p.in_scan_panel]
    traps = [p for p in panel if p.promiscuous]
    pool = held if dark else scan
    if not pool:
        raise ValueError(f"{fid}: the panel has no {'dark' if dark else 'scan'} profiles")
    chosen = [pool[int(i)] for i in rng.integers(len(pool), size=n_dom)]
    trap = traps[int(rng.integers(len(traps)))] if traps else None
    if promiscuous and trap is not None:
        chosen = [trap] + chosen
    else:
        promiscuous = False
    n_tail = _geometric(rng, arch["n_tail_mean"])
    linkers = [_geometric(rng, arch["linker_mean"]) for _ in range(len(chosen) - 1)]
    c_tail = _geometric(rng, arch["c_tail_mean"])
    classes = sorted(ev["difficulty"])
    cw = np.cumsum([float(ev["difficulty"][c]["weight"]) for c in classes])
    cls = classes[min(int(np.searchsorted(cw / cw[-1], rng.random(), side="right")), len(classes) - 1)]
    lo, hi = (float(x) for x in ev["difficulty"][cls]["identity"])
    target = q6(lo + (hi - lo) * rng.random())

    bg_cum = codons.cumulative(profiles.BACKGROUND[None, :])
    aa_parts, kind_parts, seg_parts, pi_parts = [], [], [], []
    seg_kind: list[str] = []
    pi_rows = [profiles.BACKGROUND]
    n_pi = 1
    domains: list[Domain] = []
    pos = 0

    def add(kind_name: str, aa: np.ndarray, kinds: np.ndarray, pi: np.ndarray) -> None:
        nonlocal pos
        seg_kind.append(kind_name)
        aa_parts.append(aa.astype(np.int8))
        kind_parts.append(kinds.astype(np.int8))
        seg_parts.append(np.full(len(aa), len(seg_kind) - 1, dtype=np.int16))
        pi_parts.append(pi.astype(np.int32))
        pos += len(aa)

    def tail(kind_name: str, n: int) -> None:
        a = codons.draw(bg_cum, np.zeros(n, dtype=np.int64), rng.random(n))
        add(kind_name, a, np.full(n, SITE_TAIL), np.zeros(n))

    add("start", np.array([codons.START]), np.array([SITE_START]), np.zeros(1))
    tail("ntail", n_tail)
    for i, prof in enumerate(chosen):
        if i:
            tail("linker", linkers[i - 1])
        a, node = profiles.sample_domain(prof, rng)
        is_match = node > 0
        start = pos
        add("domain", a, np.where(is_match, SITE_MATCH, SITE_INSERT),
            np.where(is_match, n_pi + node - 1, 0))
        pi_rows.append(prof.match)
        n_pi += prof.M
        n_match = int(is_match.sum())
        domains.append(Domain(prof.index, prof.accession, prof.name, start, pos, n_match,
                              int(len(node) - n_match), prof.M - n_match, prof.promiscuous, prof.in_scan_panel))
    tail("ctail", c_tail)
    add("stop", np.array([codons.STOP]), np.array([SITE_STOP]), np.zeros(1))
    aa = np.concatenate(aa_parts)
    kind = np.concatenate(kind_parts)
    rate_values = gamma_rates(ev["gamma_alpha"], ev["gamma_categories"])
    rates = rate_values[rng.integers(len(rate_values), size=len(aa))]
    rates[(kind == SITE_START) | (kind == SITE_STOP)] = 0.0
    return Design(int(family), fid, dark, promiscuous, n_dom, domains, n_tail, c_tail, linkers, cls, target,
                  aa, kind, np.concatenate(seg_parts), seg_kind, np.vstack(pi_rows), np.concatenate(pi_parts),
                  rates, rate_values)


# --- lineages ---------------------------------------------------------------------------------------------

@dataclass
class LineageModel:
    labels: list[str]
    index: dict[str, int]
    parent: list[str]
    branch_length: list[float]
    clock: np.ndarray
    gc3: np.ndarray
    cum: np.ndarray                # (n, 22, 64) cumulative codon usage per lineage


def lineage_model(branches, cfg: dict, seed: int) -> LineageModel:
    """Per species branch, its GC3 and clock multiplier (module docstring), from ``S4/lineages``."""
    ev, cod = cfg["evolution"], cfg["codons"]
    rng = seeds.generator(seed, "S4/lineages")
    g0 = float(cod["gc3_root"])
    s_gc, s_clock = float(cod["gc3_sigma"]), float(ev["clock_sigma"])
    x: dict[str, float] = {}
    labels, parents, lengths, clock, gc3 = [], [], [], [], []
    for b in branches:
        z_gc, z_clock = rng.standard_normal(2)
        base = x[b.parent] if b.parent else math.log(g0 / (1.0 - g0))
        x[b.lineage] = base + s_gc * math.sqrt(max(0.0, float(b.branch_length))) * float(z_gc)
        labels.append(b.lineage)
        parents.append(b.parent)
        lengths.append(float(b.branch_length))
        clock.append(q6(math.exp(s_clock * float(z_clock) - 0.5 * s_clock * s_clock)))
        gc3.append(q6(1.0 / (1.0 + math.exp(-x[b.lineage]))))
    cum = np.stack([codons.cumulative(codons.usage_table(g, cod["start_codons"], cod["stop_codons"]))
                    for g in gc3])
    return LineageModel(labels, {lab: i for i, lab in enumerate(labels)}, parents, lengths,
                        np.array(clock), np.array(gc3), cum)


# --- evolution (S4) ---------------------------------------------------------------------------------------

class _State:
    """One sequence along the tree: per site its residue, codon, root site of origin (or a new id), number of
    substitution events since the root, equilibrium row, rate, segment and kind."""
    __slots__ = ("aa", "codon", "origin", "events", "pi_idx", "rates", "segment", "kind")

    def __init__(self, **arrays):
        for k in self.__slots__:
            setattr(self, k, arrays[k])

    def copy(self) -> "_State":
        return _State(**{k: getattr(self, k).copy() for k in self.__slots__})

    def insert(self, at: int, **arrays) -> None:
        for k in self.__slots__:
            a = getattr(self, k)
            setattr(self, k, np.concatenate([a[:at], arrays[k].astype(a.dtype), a[at:]]))

    def delete(self, a: int, b: int) -> None:
        for k in self.__slots__:
            x = getattr(self, k)
            setattr(self, k, np.concatenate([x[:a], x[b:]]))


@dataclass
class Leaf:
    copy_label: str
    lineage: str
    protein: str
    cds: str
    path_time: float
    subst_per_site: float
    realized_subst: float
    root_identity: float
    hgt: int
    n_indels: int
    domain_spans: list[tuple[int, int]]
    homolog: np.ndarray = field(repr=False)   # (root length + 1,) residue at each root site, -1 if absent


@dataclass
class FamilyResult:
    mu: float
    expected_identity: float
    saturation: float
    n_pairs: int
    root_protein: str
    root_cds: str
    leaves: list[Leaf]
    realized_mean_identity: float = float("nan")
    realized_min_identity: float = float("nan")


def _protein(state: _State) -> str:
    return codons.LETTERS[state.aa[:-1].astype(np.int64)].tobytes().decode("ascii")


def _cds(state: _State) -> str:
    return codons.CODON_BYTES[state.codon.astype(np.int64)].tobytes().decode("ascii")


def _branch(state: _State, e: float, tab: np.ndarray, mu: float, cfg: dict, rng, pi_cum: np.ndarray,
            bg_cum: np.ndarray, rate_values: np.ndarray, indel_segs: np.ndarray, next_origin: list[int]) -> int:
    """Evolve ``state`` in place along a branch of clock-scaled length ``e`` (module docstring); returns the
    number of indel events applied."""
    ev = cfg["evolution"]
    L = len(state.aa)
    p = -np.expm1(-mu * state.rates * e)
    hit = np.nonzero(rng.random(L) < p)[0]
    if hit.size:
        new = codons.draw(pi_cum, state.pi_idx[hit], rng.random(hit.size))
        state.aa[hit] = new
        state.codon[hit] = codons.draw(tab, new, rng.random(hit.size))
        state.events[hit] += 1
    ps = -math.expm1(-float(ev["kappa_syn"]) * mu * e)
    syn = np.nonzero(rng.random(L) < ps)[0]
    if syn.size:
        state.codon[syn] = codons.draw(tab, state.aa[syn].astype(np.int64), rng.random(syn.size))
    if not len(indel_segs):
        return 0
    # Poisson(rate x (n_s + 1)) per segment, drawn as one Poisson total split multinomially (the same law)
    n_indel_sites = int(np.count_nonzero((state.kind == SITE_TAIL))) + len(indel_segs)
    total = int(rng.poisson(float(ev["indel_rate"]) * mu * e * n_indel_sites))
    if total == 0:
        return 0
    lo = np.searchsorted(state.segment, indel_segs, side="left")
    hi = np.searchsorted(state.segment, indel_segs, side="right")
    sizes = (hi - lo + 1).astype(np.float64)
    counts = rng.multinomial(total, sizes / sizes.sum())
    applied = 0
    max_len = int(ev["indel_max_length"])
    for s, n_ev in zip(indel_segs, counts):
        for _ in range(int(n_ev)):
            a = int(np.searchsorted(state.segment, s, side="left"))
            b = int(np.searchsorted(state.segment, s, side="right"))
            insertion = rng.random() < 0.5
            ell = int(rng.integers(1, max_len + 1))
            if insertion:
                at = a + int(rng.integers(0, b - a + 1))
                new = codons.draw(bg_cum, np.zeros(ell, dtype=np.int64), rng.random(ell))
                state.insert(at, aa=new, codon=codons.draw(tab, new, rng.random(ell)),
                             origin=np.arange(next_origin[0], next_origin[0] + ell),
                             events=np.zeros(ell), pi_idx=np.zeros(ell),
                             rates=rate_values[rng.integers(len(rate_values), size=ell)],
                             segment=np.full(ell, s), kind=np.full(ell, SITE_TAIL))
                next_origin[0] += ell
                applied += 1
            elif b > a:
                at = a + int(rng.integers(0, b - a))
                state.delete(at, min(at + ell, b))
                applied += 1
    return applied


def evolve_family(design: Design, nodes, lin: LineageModel, cfg: dict, rng: np.random.Generator,
                  mu: float | None = None) -> FamilyResult:
    """Evolve the family down its gene tree (``nodes``: the family's `backbone.GeneTreeNode` rows in
    preorder). ``mu`` given skips the closed-form solve (tests)."""
    ev = cfg["evolution"]
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
        raise ValueError(f"{design.family_id}: no extant copy")
    n_children = {x: 0 for x in kept}
    for x in kept:
        if by_id[x].parent != -1:
            n_children[by_id[x].parent] += 1
    depth: dict[int, float] = {}
    time_: dict[int, float] = {}
    hgt: dict[int, int] = {}
    lin_of: dict[int, int] = {}
    for x in kept:
        n = by_id[x]
        li = lin.index[n.lineage]
        lin_of[x] = li
        par = n.parent
        e = float(lin.clock[li]) * float(n.branch_length)
        depth[x] = (depth[par] if par != -1 else 0.0) + e
        time_[x] = (time_[par] if par != -1 else 0.0) + float(n.branch_length)
        recipient = par != -1 and by_id[par].kind == "transfer" and n.lineage != by_id[par].lineage
        hgt[x] = (hgt[par] if par != -1 else 0) + int(recipient)
    leaves = sorted((x for x in kept if by_id[x].extant), key=lambda x: by_id[x].copy_label)

    # path lengths for the closed-form target: pairs of extant copies, or the copy to the root
    if len(leaves) >= 2:
        ia, ib = np.triu_indices(len(leaves), 1)
        m = int(ev["target_pairs"])
        if len(ia) > m:
            pick = np.sort(rng.choice(len(ia), size=m, replace=False))
            ia, ib = ia[pick], ib[pick]
        chains = {}
        for x in leaves:
            c, y = [], x
            while y != -1:
                c.append(y)
                y = by_id[y].parent
            chains[x] = c
        T = []
        for i, j in zip(ia, ib):
            a, b = leaves[int(i)], leaves[int(j)]
            sa = set(chains[a])
            mrca = next(y for y in chains[b] if y in sa)
            T.append(depth[a] + depth[b] - 2.0 * depth[mrca])
        T = np.array(T)
    else:
        T = np.array([depth[leaves[0]]])
    sites = np.isin(design.kind, (SITE_TAIL, SITE_MATCH, SITE_INSERT))
    h = (design.pi_table[design.pi_idx[sites]] ** 2).sum(axis=1)
    if mu is None:
        mu, expected, sat = solve_rate(design.target_identity, design.rates[sites], h, T,
                                       bounds=tuple(ev["rate_bounds"]))
    else:
        vals, w, hbar = _collapse(design.rates[sites], h)
        expected, sat = _expected(float(mu), vals, w, hbar, T), float(h.mean())

    pi_cum = codons.cumulative(design.pi_table)
    bg_cum = codons.cumulative(profiles.BACKGROUND[None, :])
    L0 = len(design.aa)
    root_tab = lin.cum[lin_of[kept[0]]]
    root = _State(aa=design.aa.copy(), codon=codons.draw(root_tab, design.aa.astype(np.int64), rng.random(L0)),
                  origin=np.arange(L0, dtype=np.int64), events=np.zeros(L0, dtype=np.int32),
                  pi_idx=design.pi_idx.copy(), rates=design.rates.copy(), segment=design.segment.copy(),
                  kind=design.kind.copy())
    indel_segs = np.array([i for i, k in enumerate(design.seg_kind) if k in _INDEL_KINDS], dtype=np.int16)
    dom_segs = np.array([i for i, k in enumerate(design.seg_kind) if k == "domain"], dtype=np.int16)
    next_origin = [L0]
    states: dict[int, tuple[_State, int]] = {}
    out_leaves: dict[int, Leaf] = {}
    root_protein, root_cds = _protein(root), _cds(root)
    for x in kept:
        n = by_id[x]
        par = n.parent
        src, n_indel = (root, 0) if par == -1 else states[par]
        # the parent's last child to be evolved takes the parent's state over instead of copying it
        st = src.copy() if par == -1 or n_children[par] > 1 else src
        if par != -1:
            n_children[par] -= 1
            if n_children[par] == 0:
                del states[par]
        n_indel += _branch(st, float(lin.clock[lin_of[x]]) * float(n.branch_length), lin.cum[lin_of[x]],
                           float(mu), cfg, rng, pi_cum, bg_cum, design.rate_values, indel_segs, next_origin)
        if n_children[x]:
            states[x] = (st, n_indel)
        if n.extant:
            if not np.array_equal(codons.CODE11_RESIDUE[st.codon[1:].astype(np.int64)], st.aa[1:]):
                raise AssertionError(f"{n.copy_label}: a codon does not encode its residue under code 11")
            mask = (st.origin < L0) & (st.kind >= SITE_TAIL) & (st.kind <= SITE_INSERT)
            homolog = np.full(L0, -1, dtype=np.int8)
            homolog[st.origin[mask]] = st.aa[mask]
            lo = np.searchsorted(st.segment, dom_segs, side="left")
            hi = np.searchsorted(st.segment, dom_segs, side="right")
            out_leaves[x] = Leaf(
                n.copy_label, n.lineage, _protein(st), _cds(st), q6(time_[x]), q6(float(mu) * depth[x]),
                q6(st.events[mask].mean()), q6((st.aa[mask] == design.aa[st.origin[mask]]).mean()), hgt[x],
                n_indel, [(int(a), int(b)) for a, b in zip(lo, hi)], homolog)
    res = FamilyResult(float(mu), q6(expected), q6(sat), len(T), root_protein, root_cds,
                       [out_leaves[x] for x in leaves])
    _realized_identity(res, int(ev["exact_identity_max_copies"]), rng)
    return res


def _realized_identity(res: FamilyResult, max_exact: int, rng) -> None:
    n = len(res.leaves)
    if n < 2:
        return
    H = np.stack([leaf.homolog for leaf in res.leaves])
    valid = H >= 0
    vals = []
    if n <= max_exact:
        for i in range(n - 1):
            both = valid[i + 1:] & valid[i]
            eq = (H[i + 1:] == H[i]) & both
            vals.append(eq.sum(axis=1) / np.maximum(1, both.sum(axis=1)))
        ident = np.concatenate(vals)
    else:
        ia, ib = np.triu_indices(n, 1)
        pick = rng.choice(len(ia), size=max_exact, replace=False)
        both = valid[ia[pick]] & valid[ib[pick]]
        ident = ((H[ia[pick]] == H[ib[pick]]) & both).sum(axis=1) / np.maximum(1, both.sum(axis=1))
    res.realized_mean_identity = q6(ident.mean())
    res.realized_min_identity = q6(ident.min())


# --- tables -----------------------------------------------------------------------------------------------

class SeqFamily(NamedTuple):
    family: int
    family_id: str
    dark: bool
    promiscuous: bool
    root_source: str
    n_domains: int
    domain_hmms: tuple[str, ...]
    n_tail: int
    c_tail: int
    linkers: tuple[int, ...]
    root_length: int
    difficulty: str
    target_identity: float
    rate_multiplier: float
    expected_identity: float
    saturation_identity: float
    n_target_pairs: int
    realized_mean_identity: float
    realized_min_identity: float
    n_extant_copies: int


class SeqDomain(NamedTuple):
    family: int
    family_id: str
    index: int
    profile: int
    accession: str
    name: str
    aa_start: int
    aa_end: int
    n_match: int
    n_insert: int
    n_delete: int
    promiscuous: bool
    in_scan_panel: bool


class SeqGene(NamedTuple):
    copy_label: str
    family: int
    family_id: str
    genome_id: str
    species_node: str
    genetic_code: int
    cds_length: int
    protein_length: int
    start_codon: str
    stop_codon: str
    protein_uid: str
    gc: float
    gc3: float
    path_time: float
    subst_per_site: float
    realized_subst_per_site: float
    root_identity: float
    hgt_on_lineage: int
    n_indels: int
    domain_arch: str


class SeqLineage(NamedTuple):
    lineage: str
    parent: str
    branch_length: float
    clock: float
    gc3: float


class PanelRow(NamedTuple):
    index: int
    name: str
    accession: str
    source: str
    match_length: int
    in_scan_panel: bool
    promiscuous: bool
    ga_seq: float
    ga_dom: float
    saturation: float


_TABLES = {"families.tsv": SeqFamily, "domains.tsv": SeqDomain, "genes.tsv": SeqGene,
           "lineages.tsv": SeqLineage, "panel.tsv": PanelRow}


class _Tsv:
    def __init__(self, path: Path, cls):
        self.fh = open(path, "w", encoding="utf-8", newline="")
        self.fh.write("\t".join(cls._fields) + "\n")

    def row(self, r) -> None:
        self.fh.write("\t".join(backbone._cell(v) for v in r) + "\n")

    def close(self) -> None:
        self.fh.close()


class _Fasta:
    """FASTA at 80 columns with a samtools-style .fai written alongside, record by record."""

    def __init__(self, path: Path):
        self.fh = open(path, "wb")
        self.fai = open(str(path) + ".fai", "w", encoding="ascii", newline="")
        self.offset = 0

    def write(self, name: str, seq: str) -> None:
        head = f">{name}\n".encode("ascii")
        body = b"".join(seq[i:i + 80].encode("ascii") + b"\n" for i in range(0, len(seq), 80))
        self.fh.write(head)
        self.offset += len(head)
        self.fai.write(f"{name}\t{len(seq)}\t{self.offset}\t80\t81\n")
        self.fh.write(body)
        self.offset += len(body)

    def close(self) -> None:
        self.fh.close()
        self.fai.close()


def _gc(cds: str) -> tuple[float, float]:
    """GC of the whole CDS, and GC3 over its sense codons (initiation and stop codons left out)."""
    third = cds[5:-3:3]
    gc = (cds.count("G") + cds.count("C")) / len(cds)
    gc3 = (third.count("G") + third.count("C")) / len(third) if third else 0.0
    return q6(gc), q6(gc3)


def _stream_trees(path: Path):
    """The backbone's gene_tree_nodes.tsv, one family at a time, in file order."""
    conv = backbone._parser(backbone.GeneTreeNode)
    fam, rows = None, []
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        if tuple(header) != backbone.GeneTreeNode._fields:
            raise ValueError(f"{path}: header {header}")
        for raw in reader:
            n = backbone.GeneTreeNode(*(c(v) for c, v in zip(conv, raw)))
            if n.family != fam and rows:
                yield fam, rows
                rows = []
            fam = n.family
            rows.append(n)
    if rows:
        yield fam, rows


def check_config(cfg: dict) -> None:
    """Refuse a ``sequences`` preset section this session does not implement or cannot use."""
    if int(cfg["codons"]["genetic_code"]) != 11:
        raise NotImplementedError("only genetic code 11 is implemented (codes 4 and 25 are S8, spec §2.1)")
    for cls, spec in cfg["evolution"]["difficulty"].items():
        lo, hi = (float(x) for x in spec["identity"])
        if not 0.0 < lo <= hi < 1.0 or float(spec["weight"]) < 0:
            raise ValueError(f"difficulty {cls}: identity range must lie in (0, 1), got {spec}")
    if not cfg["architecture"]["domain_count_weights"] or min(cfg["architecture"]["domain_count_weights"]) < 0:
        raise ValueError("architecture.domain_count_weights must be non-negative and non-empty")


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@dataclass
class SeqWorld:
    seed: int
    families: list[SeqFamily]
    domains: list[SeqDomain]
    genes: list[SeqGene]
    lineages: list[SeqLineage]
    panel: list[PanelRow]
    digest: str
    out: Path
    meta: dict = field(default_factory=dict, repr=False)


def build(world, *, pfam_hmm=None) -> SeqWorld:
    """Write ``world/sequences/`` from the backbone in ``world`` (module docstring) and return it."""
    world = Path(world)
    bb = json.loads((world / "backbone" / "backbone.json").read_text(encoding="utf-8"))
    if bb["format"] != backbone.FORMAT:
        raise ValueError(f"{world}: backbone format {bb['format']!r}, expected {backbone.FORMAT!r}")
    cfg, seed = bb["preset"]["sequences"], int(bb["seed"])
    check_config(cfg)
    out = world / "sequences"
    if out.exists():
        raise FileExistsError(f"{out} exists; sequences are written once")
    tmp = world / "sequences.partial"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir()
    panel = profiles.panel(cfg["profiles"], seed, pfam_hmm)
    lin = lineage_model(backbone.species_branches(world), cfg, seed)
    genome_of = {}
    with open(world / "backbone" / "genomes.tsv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genome_of[r["species_node"]] = r["genome_id"]
    with open(world / "backbone" / "families.tsv", newline="", encoding="utf-8") as fh:
        expected_genes = sum(int(r["n_extant_copies"]) for r in csv.DictReader(fh, delimiter="\t"))

    fa_cds, fa_prot, fa_root = _Fasta(tmp / "cds.fna"), _Fasta(tmp / "proteins.faa"), _Fasta(tmp / "roots.faa")
    t_fam, t_dom, t_gene = (_Tsv(tmp / "families.tsv", SeqFamily), _Tsv(tmp / "domains.tsv", SeqDomain),
                            _Tsv(tmp / "genes.tsv", SeqGene))
    n_genes = 0
    for family, nodes in _stream_trees(world / "backbone" / "gene_tree_nodes.tsv"):
        if not any(n.extant for n in nodes):
            continue
        d = design_family(cfg, seed, family, panel)
        res = evolve_family(d, nodes, lin, cfg, seeds.generator(seed, f"S4/family/{d.family_id}"))
        fa_root.write(d.family_id, res.root_protein)
        t_fam.row(SeqFamily(family, d.family_id, d.dark, d.promiscuous, cfg["architecture"]["root_source"],
                            d.n_domains, tuple(x.accession for x in d.domains), d.n_tail, d.c_tail,
                            tuple(d.linkers), d.root_length, d.difficulty, d.target_identity, res.mu,
                            res.expected_identity, res.saturation, res.n_pairs, res.realized_mean_identity,
                            res.realized_min_identity, len(res.leaves)))
        for i, x in enumerate(d.domains):
            t_dom.row(SeqDomain(family, d.family_id, i, x.profile, x.accession, x.name, x.start, x.end,
                                x.n_match, x.n_insert, x.n_delete, x.promiscuous, x.in_scan_panel))
        for leaf in res.leaves:
            fa_cds.write(leaf.copy_label, leaf.cds)
            fa_prot.write(leaf.copy_label, leaf.protein)
            gc, gc3 = _gc(leaf.cds)
            arch = ";".join(f"{x.accession}:{a}-{b}" for x, (a, b) in zip(d.domains, leaf.domain_spans))
            t_gene.row(SeqGene(leaf.copy_label, family, d.family_id, genome_of[leaf.lineage], leaf.lineage, 11,
                               len(leaf.cds), len(leaf.protein), leaf.cds[:3], leaf.cds[-3:],
                               blake3(leaf.protein.encode("ascii")).hexdigest(length=16), gc, gc3,
                               leaf.path_time, leaf.subst_per_site, leaf.realized_subst, leaf.root_identity,
                               leaf.hgt, leaf.n_indels, arch))
            n_genes += 1
    for w in (fa_cds, fa_prot, fa_root, t_fam, t_dom, t_gene):
        w.close()
    if n_genes != expected_genes:
        raise ValueError(f"wrote {n_genes} genes, the backbone has {expected_genes} extant copies")
    t = _Tsv(tmp / "lineages.tsv", SeqLineage)
    for i, lab in enumerate(lin.labels):
        t.row(SeqLineage(lab, lin.parent[i], q6(lin.branch_length[i]), float(lin.clock[i]), float(lin.gc3[i])))
    t.close()
    t = _Tsv(tmp / "panel.tsv", PanelRow)
    for p in panel:
        ga = p.ga or (float("nan"), float("nan"))
        t.row(PanelRow(p.index, p.name, p.accession, p.source, p.M, p.in_scan_panel, p.promiscuous, q6(ga[0]),
                       q6(ga[1]), q6(p.saturation)))
    t.close()
    files = {f.name: _sha(f) for f in sorted(tmp.iterdir()) if f.is_file()}
    digest = blake3("".join(f"{k}\t{v}\n" for k, v in sorted(files.items())).encode()).hexdigest(length=16)
    versions = {"numpy": np.__version__}
    src = {"source": cfg["profiles"]["source"], "file": None, "sha256": None, "release": None}
    if cfg["profiles"]["source"] == "pfam":
        import pyhmmer
        versions["pyhmmer"] = pyhmmer.__version__
        src.update(file=Path(pfam_hmm).name, sha256=_sha(Path(pfam_hmm)),
                   release=str(cfg["profiles"]["pfam_release"]))
    doc = {"format": FORMAT, "seed": seed, "sequences_preset": cfg,
           "params_digest": backbone.params_digest(cfg), "backbone_digest": bb["digest"], "panel": src,
           "versions": versions, "stream_paths": ["S3/panel", "S3/panel/<NAME>", "S3/family/<family_id>",
                                                  "S4/lineages", "S4/family/<family_id>"],
           "n_genes": n_genes, "files": files, "digest": digest}
    (tmp / "sequences.json").write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    tmp.rename(out)
    return read(world, genes=False)


def read(world, *, genes: bool = True) -> SeqWorld:
    """Load ``world/sequences/``; ``genes=False`` leaves the per-gene table out (one row per extant copy,
    the only table that grows with the world's gene count), as `build` does to stay bounded."""
    world = Path(world)
    sdir = world / "sequences"
    doc = json.loads((sdir / "sequences.json").read_text(encoding="utf-8"))
    if doc["format"] != FORMAT:
        raise ValueError(f"{sdir}: format {doc['format']!r}, expected {FORMAT!r}")
    t = {}
    for name, cls in _TABLES.items():
        if name == "genes.tsv" and not genes:
            t[name] = []
            continue
        conv = backbone._parser(cls)
        with open(sdir / name, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            if tuple(header) != cls._fields:
                raise ValueError(f"{sdir / name}: header {header}, expected {list(cls._fields)}")
            t[name] = [cls(*(c(v) for c, v in zip(conv, row))) for row in reader]
    return SeqWorld(doc["seed"], t["families.tsv"], t["domains.tsv"], t["genes.tsv"], t["lineages.tsv"],
                    t["panel.tsv"], doc["digest"], world, doc)


@dataclass
class SequenceModel:
    """A built world's sequence model, to regenerate a family's design (its profile rows, rates and root)."""
    cfg: dict
    seed: int
    panel: list[profiles.Profile]

    def design(self, family: int) -> Design:
        return design_family(self.cfg, self.seed, family, self.panel)


def model(world, *, pfam_hmm=None) -> SequenceModel:
    doc = json.loads((Path(world) / "sequences" / "sequences.json").read_text(encoding="utf-8"))
    cfg, seed = doc["sequences_preset"], int(doc["seed"])
    return SequenceModel(cfg, seed, profiles.panel(cfg["profiles"], seed, pfam_hmm))


class Placed(NamedTuple):
    replicon: int
    position: int
    strand: int
    copy_label: str
    cds: str


def iter_genome(world, genome_id: str):
    """One genome's genes in realized order (``backbone/genes.tsv``) with strand and CDS (coding strand),
    the CDS read through ``cds.fna.fai``."""
    world = Path(world)
    fai = {}
    with open(world / "sequences" / "cds.fna.fai", encoding="ascii") as fh:
        for line in fh:
            name, length, offset, lb, lw = line.rstrip("\n").split("\t")
            fai[name] = (int(length), int(offset), int(lb), int(lw))
    with open(world / "sequences" / "cds.fna", "rb") as seq, \
            open(world / "backbone" / "genes.tsv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["genome_id"] != genome_id:
                continue
            n, off, lb, lw = fai[r["copy_label"]]
            seq.seek(off)
            raw = seq.read(n + n // lb + 1).replace(b"\n", b"")[:n]
            yield Placed(int(r["replicon"]), int(r["position"]), int(r["strand"]), r["copy_label"],
                         raw.decode("ascii"))


__all__ = ["FORMAT", "Design", "Domain", "FamilyResult", "Leaf", "LineageModel", "SeqWorld", "SequenceModel",
           "build", "read", "model", "iter_genome", "design_family", "evolve_family", "lineage_model",
           "gamma_rates", "solve_rate", "check_config", "q6"]
