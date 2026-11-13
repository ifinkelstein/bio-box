"""Group a search's loci by their flanks (card W-25; `planning/design/W-25_flanks.md`).

Each locus (an anchor gene and its window, `nbr.fetch` with edge `truncate`) is two ordered token lists, upstream and
downstream, oriented to the anchor's strand and read outward from it. A token is (id, side, strand relative to the
anchor); the id is the gene's protein's NR50 or NR90 cluster uid, or at `pfam` its domain architecture (as `hits.tsv`
writes it). A gene with no id at the level gets a token of its own in its locus, never shared.

Distances between two loci, on lists truncated to the positions both observe on each side (a side cut by the end of a
linear contig observes only its genes): J = 1 - Jaccard of the token sets; A = 1 - Jaccard of the sets of adjacent
pairs (anchor, t1), (t1, t2), ... per side (the adjacency index of BiG-SCAPE); combined = 1 - (w (1 - J) + (1 - w)
(1 - A)) with `W_COMBINED`, chosen on the planted test of tests/r2/test_w25_flanks.py. The method is reimplemented from
the papers (Domainator, BiG-SCAPE, Flankophile); no code is copied.

Loci are grouped by the anchor's own cluster; within a group, loci with equal lists and cuts collapse to one leaf (named
by the least feature_uid of its loci). Per group, at most `max_leaves` leaves (the least names) go into the pairwise
step: the J, A and positions matrices, a UPGMA tree weighted by each leaf's loci (so it equals UPGMA on the loci), and
context types from cutting it at `cut`. Only one group's matrices are held at a time (112 MB at 2,000 leaves, derived:
two float64 and one int32 K x K matrix and the UPGMA working copy).
"""

from __future__ import annotations

import atexit
import colorsys
import json
import math
import os
import re
import shutil
import sys
import tempfile
import urllib.parse
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from biobox.registry.schemas import table
from biobox.search import export, nbr

LEVELS = ("nr50", "nr90", "pfam")
GROUP_BY = ("anchor-nr90", "anchor-nr50", "none")
DISTANCES = ("j", "a", "combined")
DEFAULT_WINDOW = "5g"
DEFAULT_LEVEL = "nr50"
DEFAULT_GROUP_BY = "anchor-nr90"
DEFAULT_MAX_LEAVES = 2000
DEFAULT_CUT = 0.5
# the largest w on 0, 0.05, ..., 1 at which every planted between-type distance exceeds every within-type distance
# (design note; test_w25_w_is_the_largest_weight_on_j_that_still_separates_the_planted_types)
W_COMBINED = 0.2
CUT_TOLERANCE = 1e-9  # a merge at the cut joins a type even when its average carries rounding
ANCHOR = "^"
_SAFE = re.compile(r"[^A-Za-z0-9._-]")

README = """\
Files written by `bio-box nbr flanks` (README.txt; card W-25, planning/design/W-25_flanks.md)

A locus is one anchor gene (a gene of a hit protein) and its window of N genes on each side. Each flanking gene is a
token: its id at the level (NR50 or NR90 cluster uid of the release, or the protein's Pfam architecture), its side
(u: upstream, d: downstream, along the anchor's strand) and its strand relative to the anchor (+ same, - opposite).
A gene with no id at the level is a token of its own (~ and its feature_uid): it is never shared, so unannotated
neighbours never make two loci similar. A side is cut when a linear contig ends inside the window.

groups.tsv          # key<TAB>json facts of the run, then one row per group: the anchor cluster (or all), its loci,
                    its leaves (loci with equal token lists and cuts collapse into one leaf), the leaves in the tree,
                    sampled = more leaves than --max-loci (the leaves with the least names are in the tree), and the
                    number of context types.
context_types.tsv   one row per locus: feature_uid (the anchor gene), group, context_type (empty when the locus's
                    leaf is not in the tree), leaf (the least feature_uid of the leaf's loci), in_tree, the anchor
                    gene (contig_uid, accession, start, end, strand, 0-based half-open), anchor_protein_uid, the
                    anchor's nr90 and nr50, stage, kind, genome_uid, n_genomes and genome_uids (every genome that
                    holds the contig), tax_order (empty: no release has taxonomy yet), cut_up, cut_down (a linear
                    contig end cut that side), n_up, n_down (flanking genes seen on each side).
conservation.tsv    per group, context type and neighbour id (key): n_loci_type (the type's loci), n_loci (loci with
                    the id anywhere in their window), share_of_loci, n_neighbors (genes), n_genomes (distinct genomes
                    of those loci's contigs), offsets (offset:count along the anchor's strand, negative upstream),
                    modal_offset (ties: nearer the anchor, then downstream), same_strand_share. One row per type
                    with unannotated=true gathers the genes with no id.
groups/<group>/     leaves.tsv: one row per leaf (n_loci, in_tree, tree_order, context_type, cuts, tokens outward
                    from the anchor, its loci's feature_uids); distance_j.tsv, distance_a.tsv: J = 1 - Jaccard of the
                    token sets, A = 1 - Jaccard of the adjacent token pairs (anchor, t1), (t1, t2), ... per side, each
                    pair of leaves compared over the positions both observe on each side (4 decimals; 1 when nothing
                    is compared); positions.tsv: those positions per pair; tree.nwk: UPGMA (average linkage, leaves
                    weighted by their loci) on the distance --distance names, node heights = merge distance / 2.
                    Context types are the clades merged at a distance <= --cut, numbered by loci, then leaf name.
gff/                (--export gff) one GFF3 per locus, locus_<accession>_<feature_uid>.gff3: a CDS row per gene of
                    the window, 1-based inclusive, ID=feature_uid, protein_uid, token, offset_rel, anchor; and
                    protein_groups.tsv, protein_uid to its id at the level (its own uid when it has none).
groups/<group>/figure.svg, figure.pdf (--figure; card W-27) the group's synteny figure. One row per leaf in the tree,
                    in tree order, with the UPGMA tree at the left (the dashed line is --cut). Each row draws its
                    leaf's first locus (least feature_uid) to scale: kb from the anchor's centre along the anchor's
                    strand, upstream to the left; a gene is an arrow pointing right when it is on the anchor's strand;
                    the anchor is black and centred; // marks the end of a linear contig that cut the window. A
                    family (id at the level) has one colour in every group of the run; in a group, families held by
                    fewer than 10% of the group's loci are grey, genes with no id white. The legend names the coloured
                    families (at nr50 or nr90: the cluster uid's first 8 hex digits and its genes' commonest Pfam
                    architecture) with the share of the group's loci holding each. Columns: the drawn locus's
                    accession, the loci collapsed into the leaf, context type, genome kind, stage, taxonomy
                    (tax_order; - where the release has none) and the anchor's NR90 cluster (8 hex digits); a value
                    that differs between the leaf's loci is given with its count. At most --figure-rows rows (default
                    60): above it the figure draws the leaf with the most loci of each context type in type order,
                    then the leaves with the most loci (ties: the least name), restricts the tree to them and says so
                    in its title. No date is written in either file: the same inputs give the same bytes.

No ordinal is written: genes are named by feature_uid. Ids of clusters are of this release only.
"""


@dataclass
class Gene:
    offset_rel: int
    feature_uid: bytes
    protein_uid: bytes | None
    start: int
    end: int
    strand: int
    key: str | None  # the id at the level; None when the gene has none
    same: bool  # on the anchor's strand
    arch: str | None = None  # the protein's Pfam architecture (the figure's family names); None without a library


@dataclass
class Locus:
    feature_uid: bytes
    contig_uid: bytes
    anchor_protein_uid: bytes
    start: int
    end: int
    strand: int
    anchor_key: str | None
    nr90: bytes | None
    nr50: bytes | None
    up: list[Gene] = field(default_factory=list)
    down: list[Gene] = field(default_factory=list)
    cut_up: bool = False
    cut_down: bool = False
    accession: str | None = None
    contig_len: int | None = None
    circular: bool = False
    stage: str | None = None
    kind: str | None = None
    genome_uid: bytes | None = None
    genome_uids: tuple[bytes, ...] = ()
    tax_order: int | None = None

    def tokens(self, side: str) -> list[str]:
        genes = self.up if side == "u" else self.down
        return [token(g, side, self.feature_uid) for g in genes]


def token(g: Gene, side: str, anchor_fu: bytes) -> str:
    sign = "+" if g.same else "-"
    if g.key is None:
        return f"~{anchor_fu.hex()}{g.feature_uid.hex()}|{side}|{sign}"
    return f"{g.key}|{side}|{sign}"


def shown(g: Gene) -> str:
    """A gene's token as leaves.tsv writes it: id and relative strand, or ~feature_uid when it has no id."""
    return ("~" + g.feature_uid.hex()) if g.key is None else f"{g.key}:{'+' if g.same else '-'}"


@dataclass
class Leaf:
    name: bytes
    loci: list[Locus]
    up: list[str]
    down: list[str]
    cut_up: bool
    cut_down: bool
    window: int

    def obs(self, side: str) -> int:
        lst, cut = (self.up, self.cut_up) if side == "u" else (self.down, self.cut_down)
        return len(lst) if cut else self.window


class FlanksError(ValueError):
    pass


def parse_window(text: str) -> int:
    m = re.fullmatch(r"(\d+)g", str(text).strip().lower())
    if not m or int(m[1]) < 1:
        raise FlanksError(f"--window {text}: nbr flanks compares gene positions, so it takes a gene window such as 5g "
                          "(at least 1 gene on each side), not a length in bases")
    return int(m[1])


# --- loading the loci -----------------------------------------------------------------------------------------------


def load_loci(root: Path, anchors: set[bytes], window: int, *, release: str, files, level: str,
              library: str | None, gene_call_version: str, tape_version: str | None, threads: int | None = None,
              tmp=None) -> tuple[list[Locus], dict]:
    """The anchors' loci with their tokens at `level` and the facts of each anchor; (loci by feature_uid, facts)."""
    policy = nbr.Policy.parse(f"{window}g", "any", None, "truncate", "auto")
    rows, counts = nbr.fetch(root, anchors, policy, gene_call_version, tape_version=tape_version, files=files)
    tv = counts["tape_version"]
    with export.open_store(root, release=release, gene_call_version=gene_call_version, tape_version=tv,
                           manifest=files, threads=threads, tmp=tmp,
                           relations=("occ", "prot", "dom", "clus", "lab", "cidx", "dim", "gin")) as st:
        con = st.con
        lib = st.library(library)
        if level == "pfam" and lib is None:
            raise FlanksError(f"release {release} has no domain library, so --level pfam has no ids; use --level nr50 "
                              "or nr90")
        con.register("_nb_in", pa.Table.from_pylist(rows, schema=table("neighborhoods").schema))
        con.execute("CREATE TEMP TABLE _nb AS SELECT * FROM _nb_in")
        con.unregister("_nb_in")
        con.execute("CREATE TEMP TABLE _p AS SELECT DISTINCT protein_uid FROM _nb WHERE protein_uid IS NOT NULL")
        con.execute("CREATE TEMP TABLE _c AS SELECT DISTINCT contig_uid FROM _nb")
        export._annotate_contigs(con, with_ordinals=False, genomes=True)
        export._annotate_proteins(con, lib)
        # the anchor gene's stage: the stage of its protein_occurrence file (as hits.tsv); its kind there (genome_input)
        con.execute("""CREATE TEMP TABLE _anc AS
                       SELECT o.contig_uid, o.feature_uid, min(o.stage) AS stage FROM occ o
                       WHERE o.protein_uid IN (SELECT DISTINCT anchor_protein_uid FROM _nb)
                       GROUP BY o.contig_uid, o.feature_uid""")
        con.execute("CREATE TEMP TABLE _kind AS SELECT contig_uid, stage, min(kind) AS kind FROM gin "
                    "WHERE contig_uid IN (SELECT contig_uid FROM _c) GROUP BY contig_uid, stage")
        contig = {r[0]: r[1:] for r in con.execute(
            """SELECT c.contig_uid, l.accession, d.genome_uid, g.genome_uids FROM _c c
               LEFT JOIN _lab l USING (contig_uid) LEFT JOIN _dim d USING (contig_uid)
               LEFT JOIN _gens g USING (contig_uid) ORDER BY c.contig_uid""").fetchall()}
        anchor_facts = {(r[0], r[1]): r[2:] for r in con.execute(
            """SELECT a.contig_uid, a.feature_uid, a.stage, k.kind FROM _anc a
               LEFT JOIN _kind k ON k.contig_uid = a.contig_uid AND k.stage = a.stage
               ORDER BY a.contig_uid, a.feature_uid""").fetchall()}
        genes = export._stream(con, """
            SELECT n.contig_uid, n.anchor_feature_uid, n.anchor_protein_uid, n."offset", n.offset_rel, n.feature_uid,
                   n.start, n."end", n.strand, n.protein_uid, cl.nr90, cl.nr50, a.arch
            FROM _nb n LEFT JOIN _cl cl ON cl.protein_uid = n.protein_uid
            LEFT JOIN _arch a ON a.protein_uid = n.protein_uid
            ORDER BY n.contig_uid, n.anchor_feature_uid, n."offset" """)
        loci: dict[tuple[bytes, bytes], Locus] = {}
        neighbours: dict[tuple[bytes, bytes], list[dict]] = {}
        for r in genes:
            k = (r["contig_uid"], r["anchor_feature_uid"])
            if r["offset"] == 0:
                loci[k] = Locus(r["anchor_feature_uid"], r["contig_uid"], r["anchor_protein_uid"], r["start"],
                                r["end"], r["strand"], _key(r, level), r["nr90"], r["nr50"])
            else:
                neighbours.setdefault(k, []).append(r)
    tapes = nbr._tapes(files, tv, {cu for cu, _ in loci})
    tax = _tax_orders(files, release, {cu for cu, _ in loci})
    out = []
    for k in sorted(loci, key=lambda k: k[1]):
        loc = loci[k]
        for r in neighbours.get(k, []):
            g = Gene(r["offset_rel"], r["feature_uid"], r["protein_uid"], r["start"], r["end"], r["strand"],
                     _key(r, level), r["strand"] == loc.strand, r["arch"])
            (loc.up if g.offset_rel < 0 else loc.down).append(g)
        loc.up.sort(key=lambda g: -g.offset_rel)
        loc.down.sort(key=lambda g: g.offset_rel)
        tape = tapes[loc.contig_uid]
        loc.circular, loc.contig_len = bool(tape["is_circular"]), tape["contig_len"]
        loc.cut_up = not loc.circular and len(loc.up) < window
        loc.cut_down = not loc.circular and len(loc.down) < window
        acc, genome, genome_uids = contig.get(loc.contig_uid, (None, None, None))
        loc.accession, loc.genome_uid = acc, genome
        loc.genome_uids = tuple(bytes.fromhex(x) for x in genome_uids.split(",")) if genome_uids else ()
        loc.stage, loc.kind = anchor_facts.get(k, (None, None))
        loc.tax_order = tax.get(loc.contig_uid)
        out.append(loc)
    return out, {"tape_version": tv, "library": lib, "policy": policy, "anchors_dropped_by_edge":
                 counts["anchors_dropped_by_edge"]}


def _key(r: dict, level: str) -> str | None:
    if r["protein_uid"] is None:
        return None
    if level == "pfam":
        return r["arch"]  # None: no domain of the library
    v = r[level]
    return v.hex() if v is not None else None


def _tax_orders(files, release: str, contigs: set[bytes]) -> dict[bytes, int | None]:
    out = {}
    for path in files.table("contig_dim", release=release):
        t = pq.read_table(path, columns=["contig_uid", "tax_order"], partitioning=None,
                          filters=[("contig_uid", "in", sorted(contigs))]) if contigs else None
        for r in (t.to_pylist() if t is not None else []):
            out.setdefault(r["contig_uid"], r["tax_order"])
    return out


# --- leaves and distances -------------------------------------------------------------------------------------------


def group_of(loc: Locus, group_by: str) -> str:
    if group_by == "none":
        return "all"
    v = loc.nr90 if group_by == "anchor-nr90" else loc.nr50
    return v.hex() if v is not None else "none"


def leaves_of(loci: list[Locus], window: int) -> list[Leaf]:
    """Loci with equal token lists and cuts as one leaf each, sorted by name (the least feature_uid of its loci)."""
    by: dict[tuple, list[Locus]] = {}
    for loc in loci:
        key = (tuple(loc.tokens("u")), tuple(loc.tokens("d")), loc.cut_up, loc.cut_down)
        by.setdefault(key, []).append(loc)
    out = []
    for (up, down, cu, cd), ls in by.items():
        ls = sorted(ls, key=lambda x: x.feature_uid)
        out.append(Leaf(ls[0].feature_uid, ls, list(up), list(down), cu, cd, window))
    return sorted(out, key=lambda x: x.name)


class _Bits:
    """Each leaf's token and adjacent-pair sets as prefix bitsets per side, over one group's universe."""

    def __init__(self, leaves: list[Leaf]):
        index: dict = {}
        self.tok: list[dict[str, list[int]]] = []
        self.pair: list[dict[str, list[int]]] = []
        for leaf in leaves:
            tk, pr = {}, {}
            for side, lst in (("u", leaf.up), ("d", leaf.down)):
                t_pref, p_pref = [0], [0]
                prev = ANCHOR
                for t in lst:
                    t_bit = index.setdefault(("t", t), len(index))
                    p_bit = index.setdefault(("p", side, prev, t), len(index))
                    t_pref.append(t_pref[-1] | (1 << t_bit))
                    p_pref.append(p_pref[-1] | (1 << p_bit))
                    prev = t
                tk[side], pr[side] = t_pref, p_pref
            self.tok.append(tk)
            self.pair.append(pr)


def _jac(x: int, y: int) -> float:
    union = (x | y).bit_count()
    return 1.0 if union == 0 else 1.0 - (x & y).bit_count() / union


def matrices(leaves: list[Leaf]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(J, A, positions) for every pair of leaves."""
    n = len(leaves)
    bits = _Bits(leaves)
    mj = np.zeros((n, n), dtype=np.float64)
    ma = np.zeros((n, n), dtype=np.float64)
    mp = np.zeros((n, n), dtype=np.int32)
    obs = [(leaf.obs("u"), leaf.obs("d"), len(leaf.up), len(leaf.down)) for leaf in leaves]
    for i in range(n):
        ou_i, od_i, lu_i, ld_i = obs[i]
        ti, pi = bits.tok[i], bits.pair[i]
        mp[i, i] = ou_i + od_i
        for j in range(i + 1, n):
            ou_j, od_j, lu_j, ld_j = obs[j]
            mu, md = min(ou_i, ou_j), min(od_i, od_j)
            tj, pj = bits.tok[j], bits.pair[j]
            si = ti["u"][min(mu, lu_i)] | ti["d"][min(md, ld_i)]
            sj = tj["u"][min(mu, lu_j)] | tj["d"][min(md, ld_j)]
            qi = pi["u"][min(mu, lu_i)] | pi["d"][min(md, ld_i)]
            qj = pj["u"][min(mu, lu_j)] | pj["d"][min(md, ld_j)]
            mj[i, j] = mj[j, i] = _jac(si, sj)
            ma[i, j] = ma[j, i] = _jac(qi, qj)
            mp[i, j] = mp[j, i] = mu + md
    return mj, ma, mp


def combined(mj: np.ndarray, ma: np.ndarray, w: float = W_COMBINED) -> np.ndarray:
    return 1.0 - (w * (1.0 - mj) + (1.0 - w) * (1.0 - ma))


def upgma(d: np.ndarray, weights: list[int]) -> list[tuple[int, int, float]]:
    """Average linkage on `d` with each item weighted by `weights`: the merges (kept slot, merged slot, distance) in
    order. Each step merges the active pair with the least distance; ties go to the lexicographically least (k, l)
    with k < l, where a cluster's slot is its least member's index (the merged cluster keeps slot k)."""
    n = len(d)
    if n < 2:
        return []
    m = np.array(d, dtype=np.float64, copy=True)
    np.fill_diagonal(m, np.inf)
    size = np.array(weights, dtype=np.float64)
    active = np.ones(n, dtype=bool)
    rowmin = np.full(n, np.inf)
    rowarg = np.full(n, -1, dtype=np.int64)
    idx = np.arange(n)

    def recompute(r: int) -> None:
        row = m[r, r + 1:]
        if row.size:
            a = int(np.argmin(row))
            rowmin[r], rowarg[r] = row[a], r + 1 + a
        else:
            rowmin[r], rowarg[r] = np.inf, -1

    for r in range(n):
        recompute(r)
    merges = []
    for _ in range(n - 1):
        k = int(np.argmin(rowmin))
        l, h = int(rowarg[k]), float(rowmin[k])
        merges.append((k, l, h))
        new = (size[k] * m[k] + size[l] * m[l]) / (size[k] + size[l])
        m[k, :] = new
        m[:, k] = new
        m[k, k] = np.inf
        m[l, :] = np.inf
        m[:, l] = np.inf
        size[k] += size[l]
        active[l] = False
        rowmin[l], rowarg[l] = np.inf, -1
        recompute(k)
        rows = active & (idx < l) & (idx != k)
        stale = rows & ((rowarg == k) | (rowarg == l))
        for r in np.nonzero(stale)[0]:
            recompute(int(r))
        upd = rows & ~stale & (idx < k)
        col = m[:, k]
        better = upd & ((col < rowmin) | ((col == rowmin) & (k < rowarg)))
        rowmin[better] = col[better]
        rowarg[better] = k
    return merges


def tree(names: list[str], merges: list[tuple[int, int, float]]) -> tuple[str, list[int]]:
    """(Newick, leaf order) of the merges: node height = distance / 2, the left child the kept (lower) slot."""
    sub = {i: n for i, n in enumerate(names)}
    height = {i: 0.0 for i in range(len(names))}
    order = {i: [i] for i in range(len(names))}
    for k, l, h in merges:
        top = h / 2
        sub[k] = f"({sub[k]}:{top - height[k]:.6f},{sub[l]}:{top - height[l]:.6f})"
        height[k] = top
        order[k] = order[k] + order[l]
        del sub[l], height[l], order[l]
    (root,) = sub
    return sub[root] + ";", order[root]


def cut_types(n: int, merges: list[tuple[int, int, float]], cut: float) -> list[int]:
    """The cluster slot of each item after the merges at a distance <= cut."""
    held = {i: [i] for i in range(n)}
    for k, l, h in merges:
        if h <= cut + CUT_TOLERANCE:
            held[k] += held.pop(l)
    of = [0] * n
    for s, items in held.items():
        for i in items:
            of[i] = s
    return of


# --- conservation ---------------------------------------------------------------------------------------------------


def conservation(loci: list[Locus]) -> list[dict]:
    """One row per neighbour id over `loci` (one context type), then one row for the genes with no id."""
    stats: dict[str | None, dict] = {}
    for loc in loci:
        for g in loc.up + loc.down:
            s = stats.setdefault(g.key, {"loci": set(), "n": 0, "offsets": Counter(), "same": 0, "genomes": set()})
            s["loci"].add(loc.feature_uid)
            s["n"] += 1
            s["offsets"][g.offset_rel] += 1
            s["same"] += g.same
            s["genomes"].update(loc.genome_uids or ((loc.genome_uid,) if loc.genome_uid else (loc.contig_uid,)))
    rows = []
    for key, s in stats.items():
        offs = s["offsets"]
        modal = max(offs, key=lambda o: (offs[o], -abs(o), o))
        rows.append({"key": key or "", "unannotated": key is None, "n_loci_type": len(loci), "n_loci": len(s["loci"]),
                     "share_of_loci": round(len(s["loci"]) / len(loci), 4), "n_neighbors": s["n"],
                     "n_genomes": len(s["genomes"]), "offsets": ",".join(f"{o}:{offs[o]}" for o in sorted(offs)),
                     "modal_offset": modal, "same_strand_share": round(s["same"] / s["n"], 4)})
    return sorted(rows, key=lambda r: (r["unannotated"], -r["n_loci"], r["key"]))


# --- the figure (card W-27) -----------------------------------------------------------------------------------------
#
# One figure per group, groups/<group>/figure.svg and figure.pdf: a row per leaf in the tree, in tree order, with the
# UPGMA tree at the left and columns at the right. Each row draws its leaf's first locus (the least feature_uid, the
# leaf's name) to scale: x is bp from the anchor's centre along the anchor's strand (upstream left), a gene is an arrow
# pointing right when it is on the anchor's strand. Colours: a family (the gene's id at the level) has one colour in
# every group of a run (ranked over the run by loci holding it, then id); in a group, a family held by fewer than
# FIGURE_MIN_SHARE of the group's loci is grey; a gene with no id is white; the anchor is black. Above `figure_rows`
# leaves the figure draws the leaf with the most loci of each context type (in type order), then the leaves with the
# most loci (ties: the least name), and restricts the tree to them. The figure is planned as data (`figure_plan`, what
# the tests check) and drawn by `render_figure` with matplotlib's SVG and PDF backends only (no pyplot, no display),
# the bundled DejaVu fonts named by file, a fixed `svg.hashsalt` and no date in either file, so the bytes depend only
# on the inputs and the locked matplotlib.

DEFAULT_FIGURE_ROWS = 60
FIGURE_MIN_SHARE = 0.1
FIGURE_ANCHOR = "#222222"
FIGURE_GREY = "#bdbdbd"
FIGURE_NONE = "#ffffff"
FIGURE_FILES = ("figure.svg", "figure.pdf")
# 20 distinct colours, none grey or black (those mean rare, no id and anchor); further families take golden-angle hues
FIGURE_PALETTE = ("#e69f00", "#56b4e9", "#009e73", "#f0e442", "#0072b2", "#d55e00", "#cc79a7", "#8c564b", "#17becf",
                  "#bcbd22", "#9467bd", "#ff9da7", "#2ca02c", "#aec7e8", "#ffbb78", "#98df8a", "#c5b0d5", "#c49c94",
                  "#9edae5", "#f7b6d2")
FIGURE_COLUMNS = ("accession", "loci", "type", "kind", "stage", "taxon", "anchor NR90")
_COLUMN_IN = (0.95, 0.38, 0.36, 1.2, 0.85, 0.5, 0.75)  # column widths, inches
_COLUMN_CHARS = (16, 7, 5, 22, 15, 8, 8)  # longest text per column before it is cut with "~"
_SVG_HASHSALT = "bio-box W-27"
_ROW_IN, _TREE_IN, _GENES_IN, _LEGEND_ROW_IN = 0.2, 1.7, 7.4, 0.14


@dataclass
class FigGene:
    feature_uid: str
    key: str | None
    x0: float
    x1: float
    forward: bool  # on the anchor's strand: the arrow points right
    colour: str
    anchor: bool = False


@dataclass
class FigRow:
    leaf: str
    drawn: str  # the feature_uid of the locus drawn (the leaf's first)
    n_loci: int
    context_type: int | None
    genes: list[FigGene]
    breaks: list[float]  # x of a linear contig's end that cut a side
    columns: list[str]


@dataclass
class FigurePlan:
    group: str
    title: list[str]
    rows: list[FigRow]
    segments: list[tuple[float, float, float, float]]  # the tree: (distance, row, distance, row) line ends
    tree_max: float
    cut: float
    half_width: float  # bp on each side of the anchor's centre
    legend: list[tuple[str, str]]  # (colour, label)
    columns: list[str] = field(default_factory=lambda: list(FIGURE_COLUMNS))
    xlabel: str = ""
    distance: str = "j"


def figure_unavailable() -> str | None:
    """Why --figure cannot run here (matplotlib not importable), or None."""
    try:
        _mpl()
    except ImportError as e:
        return (f"--figure needs matplotlib, which cannot be imported here ({e}); it is a locked dependency (uv.lock), "
                "so install the locked environment or run without --figure")
    return None


def _mpl():
    """matplotlib, imported on first use; its cache goes to a private temporary directory unless MPLCONFIGDIR is set
    (never the home directory: on the pod that is on /stor, AGENTS §7.3)."""
    if "matplotlib" not in sys.modules and "MPLCONFIGDIR" not in os.environ:
        d = tempfile.mkdtemp(prefix="biobox-mpl-")
        os.environ["MPLCONFIGDIR"] = d
        atexit.register(shutil.rmtree, d, True)
    import matplotlib
    from matplotlib.figure import Figure
    from matplotlib.font_manager import FontProperties
    from matplotlib.patches import Patch, Polygon

    return matplotlib, Figure, FontProperties, Patch, Polygon


def figure_colour(i: int) -> str:
    """The run's i-th family colour: the palette, then golden-angle hues at alternating saturation and value."""
    if i < len(FIGURE_PALETTE):
        return FIGURE_PALETTE[i]
    k = i - len(FIGURE_PALETTE)
    r, g, b = colorsys.hsv_to_rgb((k * 0.6180339887) % 1.0, 0.45 + 0.3 * (k % 2), 0.95 - 0.2 * ((k // 2) % 2))
    return "#%02x%02x%02x" % (round(r * 255), round(g * 255), round(b * 255))


def group_shares(loci: list[Locus]) -> dict[str, float]:
    """Per family (id at the level) in the flanks of `loci`: the share of the loci that hold it."""
    held: Counter = Counter()
    for loc in loci:
        held.update({g.key for g in loc.up + loc.down if g.key is not None})
    return {k: held[k] / len(loci) for k in sorted(held)}


def run_colours(groups: dict[str, list[Locus]]) -> dict[str, str]:
    """One colour per family coloured in any group of the run (share >= FIGURE_MIN_SHARE there), ranked by the loci
    of the run holding it (descending), then id: the same family has the same colour in every group's figure."""
    held: Counter = Counter()
    coloured: set[str] = set()
    for gid in sorted(groups):
        shares = group_shares(groups[gid])
        coloured.update(k for k, s in shares.items() if s >= FIGURE_MIN_SHARE)
        for loc in groups[gid]:
            held.update({g.key for g in loc.up + loc.down if g.key is not None})
    ranked = sorted(coloured, key=lambda k: (-held[k], k))
    return {k: figure_colour(i) for i, k in enumerate(ranked)}


def family_names(loci: list[Locus], level: str) -> dict[str, str]:
    """Each family's legend name: at pfam its architecture; at nr50 or nr90 the cluster uid's first 8 hex digits and
    the most common Pfam architecture of its genes in the run (ties: the least), when the release has a library."""
    archs: dict[str, Counter] = {}
    for loc in loci:
        for g in loc.up + loc.down:
            if g.key is not None:
                c = archs.setdefault(g.key, Counter())
                if g.arch:
                    c[g.arch] += 1
    out = {}
    for k in sorted(archs):
        if level == "pfam":
            out[k] = k
            continue
        c = archs[k]
        best = min(c, key=lambda a: (-c[a], a)) if c else None
        out[k] = f"{k[:8]} {best}" if best else k[:8]
    return out


def _cut_text(s: str, n: int) -> str:
    s = str(s).replace("$", "\\$")  # never mathtext
    return s if len(s) <= n else s[:n - 1] + "~"


def _agg(values: list) -> str:
    """One value, or each value with its count (most first, then by name) when the leaf's loci differ."""
    c = Counter("-" if v is None else str(v) for v in values)
    if len(c) == 1:
        return next(iter(c))
    return " ".join(f"{v}:{c[v]}" for v in sorted(c, key=lambda v: (-c[v], v)))


def _x_of(loc: Locus, g: Gene, side: str) -> tuple[float, float]:
    """A gene's span in bp from the anchor's centre along the anchor's strand; on a circular contig the copy of it on
    its side of the anchor nearest the anchor."""
    centre, s = (loc.start + loc.end) / 2, loc.strand
    best = None
    for shift in ((0, -loc.contig_len, loc.contig_len) if loc.circular and loc.contig_len else (0,)):
        a, b = sorted(((g.start + shift - centre) * s, (g.end + shift - centre) * s))
        mid = (a + b) / 2
        if (mid < 0) == (side == "u") and (best is None or abs(mid) < abs(sum(best) / 2)):
            best = (a, b)
    if best is None:
        best = tuple(sorted(((g.start - centre) * s, (g.end - centre) * s)))
    return best


def _breaks(loc: Locus) -> list[float]:
    """x of each contig end that cut a side: upstream lies at lower coordinates for a + anchor."""
    centre, s = (loc.start + loc.end) / 2, loc.strand
    out = []
    if loc.cut_up:
        out.append(((0 if s == 1 else loc.contig_len) - centre) * s)
    if loc.cut_down:
        out.append(((loc.contig_len if s == 1 else 0) - centre) * s)
    return out


def figure_leaves(kept: list[Leaf], ctype: dict[int, int], tree_order: list[int], rows: int) -> list[int]:
    """The leaves (indices into `kept`) the figure draws, in tree order: all when at most `rows`; else the leaf with
    the most loci of each context type in type order (ties: the least name), then the leaves with the most loci."""
    rank = {i: r for r, i in enumerate(tree_order)}
    if len(kept) <= rows:
        return sorted(range(len(kept)), key=rank.get)
    best = lambda ix: sorted(ix, key=lambda i: (-len(kept[i].loci), kept[i].name))  # noqa: E731
    by_type: dict[int, list[int]] = {}
    for i in range(len(kept)):
        by_type.setdefault(ctype[i], []).append(i)
    chosen: list[int] = []
    for t in sorted(by_type):
        if len(chosen) < rows:
            chosen.append(best(by_type[t])[0])
    taken = set(chosen)
    chosen += [i for i in best(range(len(kept))) if i not in taken][:rows - len(chosen)]
    return sorted(chosen, key=rank.get)


def _segments(n: int, merges: list[tuple[int, int, float]], y: dict[int, int]) -> list[tuple]:
    """The tree restricted to the drawn leaves (`y`: index -> row): each merge of two drawn subtrees is a junction at
    its distance, each subtree a line from its own height to it. A subtree without drawn leaves vanishes; one that
    meets such a subtree keeps its height."""
    node: dict[int, tuple[float, float] | None] = {i: ((float(y[i]), 0.0) if i in y else None) for i in range(n)}
    out = []
    for k, l, h in merges:
        a, b = node[k], node[l]
        if a and b:
            for yy, xx in (a, b):
                if xx != h:
                    out.append((xx, yy, h, yy))
            out.append((h, a[0], h, b[0]))
            node[k] = ((a[0] + b[0]) / 2, h)
        else:
            node[k] = a or b
        del node[l]
    return out


def figure_plan(gid: str, group_loci: list[Locus], kept: list[Leaf], merges, tree_order: list[int],
                ctype: dict[int, int], colours: dict[str, str], names: dict[str, str], *, rows: int, facts: dict,
                window: int, distance: str, cut: float, n_leaves: int) -> FigurePlan:
    """What one group's figure draws (no matplotlib here)."""
    shares = group_shares(group_loci)
    coloured = {k for k, s in shares.items() if s >= FIGURE_MIN_SHARE}
    shown = figure_leaves(kept, ctype, tree_order, rows)
    y = {i: r for r, i in enumerate(shown)}
    used: set[str | None] = set()
    out_rows = []
    for i in shown:
        leaf = kept[i]
        loc = leaf.loci[0]
        genes = [FigGene(loc.feature_uid.hex(), loc.anchor_key, -(loc.end - loc.start) / 2, (loc.end - loc.start) / 2,
                         True, FIGURE_ANCHOR, True)]
        for side, lst in (("u", loc.up), ("d", loc.down)):
            for g in lst:
                x0, x1 = _x_of(loc, g, side)
                colour = FIGURE_NONE if g.key is None else (colours[g.key] if g.key in coloured else FIGURE_GREY)
                used.add(g.key if g.key in coloured or g.key is None else "")
                genes.append(FigGene(g.feature_uid.hex(), g.key, x0, x1, g.same, colour))
        cols = [loc.accession or loc.contig_uid.hex()[:12], str(len(leaf.loci)), str(ctype[i]),
                _agg([x.kind for x in leaf.loci]), _agg([x.stage for x in leaf.loci]),
                _agg([x.tax_order for x in leaf.loci]), _agg([x.nr90.hex()[:8] if x.nr90 else None for x in leaf.loci])]
        out_rows.append(FigRow(leaf.name.hex(), loc.feature_uid.hex(), len(leaf.loci), ctype[i], genes, _breaks(loc),
                               cols))
    segments = _segments(len(kept), merges, y)
    reach = [abs(v) for r in out_rows for g in r.genes for v in (g.x0, g.x1)] + \
            [abs(v) for r in out_rows for v in r.breaks]
    half = max(reach) * 1.03 if reach else 1.0
    legend = [(FIGURE_ANCHOR, "anchor (the search's hit), centred")]
    for k in sorted((k for k in used if k), key=lambda k: (-shares[k], colours[k], k)):
        legend.append((colours[k], f"{names.get(k, k[:8])}  {round(100 * shares[k])}%"))
    if "" in used:
        legend.append((FIGURE_GREY, f"other families (each in < {FIGURE_MIN_SHARE:.0%} of the group's loci)"))
    if None in used:
        legend.append((FIGURE_NONE, f"no id at level {facts['level']} (never shared)"))
    n_loci = len(group_loci)
    title = [f"nbr flanks: group {gid}, {n_loci} loci in {n_leaves} leaves ({len(kept)} in the tree, "
             f"{len(set(ctype.values()))} context types)",
             f"release {facts['release']}, level {facts['level']}, window {window}g, distance {distance}, cut {cut:g}; "
             f"a family's colour is the same in every group of the run; % = share of the group's loci holding it"]
    if len(shown) == len(kept):
        title.append(f"rows: all {len(kept)} leaves in the tree, in tree order; each row draws its leaf's first locus "
                     "(least feature_uid); loci = the loci collapsed into the leaf")
    else:
        title.append(f"rows: {len(shown)} of {len(kept)} leaves (--figure-rows {rows}): the leaf with the most loci of "
                     "each context type, then the leaves with the most loci; the tree is restricted to them")
    if n_leaves > len(kept):
        title.append(f"{n_leaves - len(kept)} leaves beyond --max-loci are not in the tree and are not drawn")
    return FigurePlan(gid, title, out_rows, segments, max((s[2] for s in segments), default=0.0), cut, half, legend,
                      xlabel="kb from the anchor's centre along its strand (upstream left); // = end of a linear "
                             "contig", distance=distance)


def _ticks(lo: float, hi: float, target: int = 8) -> list[float]:
    raw = (hi - lo) / target
    mag = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1.0
    step = next(m * mag for m in (1, 2, 5, 10) if m * mag >= raw)
    first = math.ceil(lo / step)
    return [k * step for k in range(first, math.floor(hi / step) + 1)]


def _fmt_tick(v: float, scale: float) -> str:
    s = f"{v / scale:.6g}"
    return "0" if s in ("-0", "0") else s


def render_figure(plan: FigurePlan, gdir: Path) -> None:
    """Draw `plan` into gdir/figure.svg and figure.pdf (each under a temporary name, then renamed); one figure at a
    time, cleared after writing."""
    matplotlib, Figure, FontProperties, Patch, Polygon = _mpl()
    fonts = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
    regular, bold = str(fonts / "DejaVuSans.ttf"), str(fonts / "DejaVuSans-Bold.ttf")
    fp = lambda size, path=regular: FontProperties(fname=path, size=size)  # noqa: E731
    n = max(len(plan.rows), 1)
    ncol = 3
    legend_rows = math.ceil(len(plan.legend) / ncol)
    left, gap, cols_in = 0.3, 0.12, sum(_COLUMN_IN)
    width = left + _TREE_IN + gap + _GENES_IN + 0.2 + cols_in + 0.25
    head_in = 0.22 + 0.16 * len(plan.title) + 0.28
    rows_in = n * _ROW_IN
    axis_in = 0.5
    legend_in = legend_rows * _LEGEND_ROW_IN + 0.1
    height = head_in + rows_in + axis_in + legend_in + 0.15
    rc = {"svg.hashsalt": _SVG_HASHSALT, "svg.fonttype": "none", "pdf.fonttype": 3, "pdf.compression": 6,
          "font.family": "DejaVu Sans", "axes.unicode_minus": False, "path.simplify": False,
          "figure.dpi": 72, "savefig.dpi": 72, "savefig.bbox": None}
    with matplotlib.rc_context(rc):
        fig = Figure(figsize=(width, height))

        def box(x_in, w_in):
            return fig.add_axes((x_in / width, (axis_in + legend_in + 0.15) / height, w_in / width, rows_in / height))

        for k, line in enumerate(plan.title):
            fig.text(left / width, 1 - (0.2 + 0.16 * k) / height, _cut_text(line, 230), va="top", ha="left",
                     fontproperties=fp(8 if k == 0 else 6.5, bold if k == 0 else regular))
        tree_ax = box(left, _TREE_IN)
        gene_ax = box(left + _TREE_IN + gap, _GENES_IN)
        col_ax = box(left + _TREE_IN + gap + _GENES_IN + 0.2, cols_in)
        for ax in (tree_ax, gene_ax, col_ax):
            ax.set_ylim(n - 0.5, -0.5)
            for side in ("top", "right", "left"):
                ax.spines[side].set_visible(False)
            ax.set_yticks([])
        # tree: distance from 0 (the leaves, at the right) to the root (left)
        tmax = max(plan.tree_max, plan.cut if plan.cut <= 1 else 0, 1e-6) * 1.04
        tree_ax.set_xlim(tmax, 0)
        for x0, y0, x1, y1 in plan.segments:
            tree_ax.plot([x0, x1], [y0, y1], color="#333333", linewidth=0.7, solid_capstyle="butt")
        tree_ax.axvline(plan.cut, color="#b2182b", linewidth=0.6, linestyle=(0, (2, 2)))
        tree_ax.text(plan.cut, -0.75, f"cut {plan.cut:g}", ha="center", va="bottom", color="#b2182b",
                     fontproperties=fp(6))
        tt = [t for t in _ticks(0, tmax, 4) if t <= tmax]
        tree_ax.set_xticks(tt)
        tree_ax.set_xticklabels([_fmt_tick(t, 1) for t in tt], fontproperties=fp(6))
        tree_ax.set_xlabel(f"merge distance ({plan.distance})", fontproperties=fp(6.5))
        # genes
        half = plan.half_width
        gene_ax.set_xlim(-half, half)
        head = 0.022 * 2 * half
        gene_ax.axvline(0, color="#999999", linewidth=0.4, linestyle=(0, (1, 2)), zorder=0)
        for r, row in enumerate(plan.rows):
            if row.context_type is not None and row.context_type % 2 == 0:
                gene_ax.axhspan(r - 0.5, r + 0.5, color="#f2f2f2", linewidth=0, zorder=0)
            xs = [v for g in row.genes for v in (g.x0, g.x1)] + row.breaks
            gene_ax.plot([min(xs), max(xs)], [r, r], color="#8c8c8c", linewidth=0.6, zorder=1)
            for g in row.genes:
                hl = min(head, 0.5 * (g.x1 - g.x0))
                h = 0.34 if g.anchor else 0.3
                if g.forward:
                    pts = [(g.x0, r - h), (g.x1 - hl, r - h), (g.x1, r), (g.x1 - hl, r + h), (g.x0, r + h)]
                else:
                    pts = [(g.x1, r - h), (g.x0 + hl, r - h), (g.x0, r), (g.x0 + hl, r + h), (g.x1, r + h)]
                gene_ax.add_patch(Polygon(pts, closed=True, facecolor=g.colour, edgecolor="#000000" if g.anchor
                                          else "#4d4d4d", linewidth=0.8 if g.anchor else 0.4, zorder=2))
            for b in row.breaks:
                d = 0.006 * 2 * half
                for off in (-d * 0.6, d * 0.6):
                    gene_ax.plot([b + off - d / 2, b + off + d / 2], [r + 0.38, r - 0.38], color="#b2182b",
                                 linewidth=1.0, zorder=3)
        ticks = _ticks(-half, half)
        gene_ax.set_xticks(ticks)
        gene_ax.set_xticklabels([_fmt_tick(t, 1000) for t in ticks], fontproperties=fp(6))
        gene_ax.set_xlabel(plan.xlabel, fontproperties=fp(6.5))
        # columns
        col_ax.set_xlim(0, cols_in)
        col_ax.set_xticks([])
        col_ax.spines["bottom"].set_visible(False)
        xs = [sum(_COLUMN_IN[:k]) for k in range(len(_COLUMN_IN))]
        for name, x in zip(plan.columns, xs):
            col_ax.text(x, -0.75, name, ha="left", va="bottom", fontproperties=fp(6.5, bold))
        for r, row in enumerate(plan.rows):
            for k, (v, x) in enumerate(zip(row.columns, xs)):
                col_ax.text(x, r, _cut_text(v, _COLUMN_CHARS[k]), ha="left", va="center", fontproperties=fp(6.5))
        # legend, below everything
        lax = fig.add_axes((left / width, 0.1 / height, (width - left - 0.25) / width, legend_in / height))
        lax.set_axis_off()
        handles = [Patch(facecolor=c, edgecolor="#4d4d4d", linewidth=0.4) for c, _ in plan.legend]
        lax.legend(handles, [_cut_text(lab, 70) for _, lab in plan.legend], loc="upper left", ncol=ncol,
                   prop=fp(6.5), frameon=False, borderaxespad=0, handlelength=1.6, handleheight=0.9,
                   columnspacing=1.5, labelspacing=0.35)
        for name, fmt, meta in (("figure.svg", "svg", {"Date": None}),
                                ("figure.pdf", "pdf", {"CreationDate": None, "ModDate": None})):
            tmp = gdir / f".{name}.tmp-{os.getpid()}"
            try:
                fig.savefig(tmp, format=fmt, metadata=meta)
                os.replace(tmp, gdir / name)
            finally:
                if tmp.exists():
                    tmp.unlink()
        fig.clear()
    del fig


# --- writing --------------------------------------------------------------------------------------------------------


def _cell(v) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, bytes):
        return v.hex()
    return str(v).replace("\t", " ").replace("\n", " ")


def _tsv(cols: list[str], rows) -> Iterator[str]:
    yield "\t".join(cols)
    for r in rows:
        yield "\t".join(_cell(r[c]) for c in cols)


CONTEXT_COLUMNS = ["feature_uid", "group", "context_type", "leaf", "in_tree", "contig_uid", "accession", "start", "end",
                   "strand", "anchor_protein_uid", "anchor_nr90", "anchor_nr50", "stage", "kind", "genome_uid",
                   "n_genomes", "genome_uids", "tax_order", "cut_up", "cut_down", "n_up", "n_down"]
CONSERVATION_COLUMNS = ["group", "context_type", "key", "unannotated", "n_loci_type", "n_loci", "share_of_loci",
                        "n_neighbors", "n_genomes", "offsets", "modal_offset", "same_strand_share"]
GROUP_COLUMNS = ["group", "n_loci", "n_leaves", "n_leaves_in_tree", "sampled", "n_context_types"]
LEAF_COLUMNS = ["leaf", "n_loci", "in_tree", "tree_order", "context_type", "cut_up", "cut_down", "n_up", "n_down",
                "tokens_up", "tokens_down", "feature_uids"]


def _matrix(names: list[str], m: np.ndarray, fmt: str) -> Iterator[str]:
    yield "\t".join(["leaf", *names])
    for i, name in enumerate(names):
        yield name + "\t" + "\t".join(fmt % v for v in m[i])


def _gff(loc: Locus) -> list[str]:
    seqid = loc.accession or loc.contig_uid.hex()
    q = lambda s: urllib.parse.quote(str(s), safe="+-._:|~ ")  # noqa: E731  (GFF3 escapes ; = & , % tab)
    lines = ["##gff-version 3", f"##sequence-region {q(seqid)} 1 {loc.contig_len}"]
    if loc.circular:
        lines.append(f"{q(seqid)}\tbio-box\tregion\t1\t{loc.contig_len}\t.\t+\t.\tID={loc.contig_uid.hex()};"
                     "Is_circular=true")
    anchor = Gene(0, loc.feature_uid, loc.anchor_protein_uid, loc.start, loc.end, loc.strand, loc.anchor_key, True)
    for g in sorted(loc.up + [anchor] + loc.down, key=lambda g: (g.start, g.end, g.feature_uid)):
        tok = ("~" + g.feature_uid.hex()) if g.key is None else g.key
        attrs = [f"ID={g.feature_uid.hex()}"]
        if g.protein_uid is not None:
            attrs += [f"Name={g.protein_uid.hex()}", f"protein_uid={g.protein_uid.hex()}"]
        attrs += [f"token={q(tok)}", f"offset_rel={g.offset_rel}", f"anchor={'true' if g is anchor else 'false'}"]
        lines.append(f"{q(seqid)}\tbio-box\tCDS\t{g.start + 1}\t{g.end}\t.\t{'+' if g.strand == 1 else '-'}\t0\t"
                     + ";".join(attrs))
    return lines


def write(out_dir: Path, loci: list[Locus], facts: dict, *, window: int, group_by: str, distance: str, cut: float,
          max_leaves: int, export_gff: bool, notes: list[str], figure: bool = False,
          figure_rows: int = DEFAULT_FIGURE_ROWS) -> dict:
    """Write every output into `out_dir` (which exists and is empty); one group's matrices, and with `figure` one
    group's figure, at a time. Returns counts for the stderr lines; appends to `notes` one line per sampled group."""
    groups: dict[str, list[Locus]] = {}
    for loc in loci:
        groups.setdefault(group_of(loc, group_by), []).append(loc)
    order = sorted(groups, key=lambda g: (-len(groups[g]), g))
    ctx_rows, cons_rows, group_rows = [], [], []
    peak_leaves = 0
    n_figures = 0
    if figure:
        colours, fam_names = run_colours(groups), family_names(loci, facts["level"])
    for gid in order:
        leaves = leaves_of(groups[gid], window)
        kept = leaves[:max_leaves]
        sampled = len(leaves) > len(kept)
        if sampled:
            notes.append(f"note: group {gid}: {len(leaves)} leaves after collapsing, above --max-loci {max_leaves}: the "
                         f"{max_leaves} with the least names are in the tree, the other {len(leaves) - len(kept)} "
                         "(their loci) have in_tree=false and no context type")
        peak_leaves = max(peak_leaves, len(kept))
        names = [leaf.name.hex() for leaf in kept]
        mj, ma, mp = matrices(kept)
        d = {"j": mj, "a": ma, "combined": None}[distance]
        if d is None:
            d = combined(mj, ma)
        merges = upgma(d, [len(leaf.loci) for leaf in kept])
        gdir = out_dir / "groups" / gid
        gdir.mkdir(parents=True)
        export.write_lines(_matrix(names, mj, "%.4f"), gdir / "distance_j.tsv")
        export.write_lines(_matrix(names, ma, "%.4f"), gdir / "distance_a.tsv")
        export.write_lines(_matrix(names, mp, "%d"), gdir / "positions.tsv")
        newick, tree_order = tree(names, merges)
        export.write_lines([newick], gdir / "tree.nwk")
        del mj, ma, mp, d
        slot = cut_types(len(kept), merges, cut)
        members: dict[int, list[int]] = {}
        for i, s in enumerate(slot):
            members.setdefault(s, []).append(i)
        ranked = sorted(members.values(), key=lambda ix: (-sum(len(kept[i].loci) for i in ix), kept[ix[0]].name))
        ctype = {i: t + 1 for t, ix in enumerate(ranked) for i in ix}
        rank = {i: r for r, i in enumerate(tree_order)}
        leaf_rows = []
        for i, leaf in enumerate(leaves):
            in_tree = i < len(kept)
            rep = leaf.loci[0]
            leaf_rows.append({"leaf": leaf.name, "n_loci": len(leaf.loci), "in_tree": in_tree,
                              "tree_order": rank.get(i) if in_tree else None,
                              "context_type": ctype.get(i) if in_tree else None, "cut_up": leaf.cut_up,
                              "cut_down": leaf.cut_down, "n_up": len(leaf.up), "n_down": len(leaf.down),
                              "tokens_up": ",".join(shown(g) for g in rep.up),
                              "tokens_down": ",".join(shown(g) for g in rep.down),
                              "feature_uids": ",".join(x.feature_uid.hex() for x in leaf.loci)})
            for loc in leaf.loci:
                ctx_rows.append({"feature_uid": loc.feature_uid, "group": gid,
                                 "context_type": ctype.get(i) if in_tree else None, "leaf": leaf.name,
                                 "in_tree": in_tree, "contig_uid": loc.contig_uid, "accession": loc.accession,
                                 "start": loc.start, "end": loc.end, "strand": f"{loc.strand:+d}",
                                 "anchor_protein_uid": loc.anchor_protein_uid, "anchor_nr90": loc.nr90,
                                 "anchor_nr50": loc.nr50, "stage": loc.stage, "kind": loc.kind,
                                 "genome_uid": loc.genome_uid, "n_genomes": len(loc.genome_uids) or None,
                                 "genome_uids": ",".join(x.hex() for x in loc.genome_uids) or None,
                                 "tax_order": loc.tax_order, "cut_up": loc.cut_up, "cut_down": loc.cut_down,
                                 "n_up": len(loc.up), "n_down": len(loc.down)})
        export.write_lines(_tsv(LEAF_COLUMNS, leaf_rows), gdir / "leaves.tsv")
        if figure:
            render_figure(figure_plan(gid, groups[gid], kept, merges, tree_order, ctype, colours, fam_names,
                                      rows=figure_rows, facts=facts, window=window, distance=distance, cut=cut,
                                      n_leaves=len(leaves)), gdir)
            n_figures += 1
        for t, ix in enumerate(ranked):
            for r in conservation([loc for i in ix for loc in kept[i].loci]):
                cons_rows.append({"group": gid, "context_type": t + 1, **r})
        group_rows.append({"group": gid, "n_loci": len(groups[gid]), "n_leaves": len(leaves),
                           "n_leaves_in_tree": len(kept), "sampled": sampled, "n_context_types": len(ranked)})
    ctx_rows.sort(key=lambda r: (order.index(r["group"]), r["context_type"] is None, r["context_type"] or 0,
                                 r["feature_uid"]))
    export.write_lines(_tsv(CONTEXT_COLUMNS, ctx_rows), out_dir / "context_types.tsv")
    export.write_lines(_tsv(CONSERVATION_COLUMNS, cons_rows), out_dir / "conservation.tsv")
    head = {"release": facts["release"], "tape_version": facts["tape_version"], "level": facts["level"],
            "library": facts["library"], "window": f"{window}g", "group_by": group_by, "distance": distance,
            "w": W_COMBINED if distance == "combined" else None, "cut": cut, "max_loci": max_leaves,
            "n_anchor_proteins_given": facts["n_anchor_proteins_given"], "n_loci": len(loci),
            "anchors_dropped_by_edge": facts["anchors_dropped_by_edge"], "n_groups": len(order),
            "policy_hash": facts["policy"].hash()}
    export.write_lines([f"# {k}\t{json.dumps(v, separators=(',', ':'))}" for k, v in head.items()]
                       + list(_tsv(GROUP_COLUMNS, group_rows)), out_dir / "groups.tsv")
    (out_dir / "README.txt").write_text(README, encoding="utf-8")
    n_gff = 0
    if export_gff:
        gff = out_dir / "gff"
        gff.mkdir()
        proteins: dict[bytes, str] = {}
        for loc in sorted(loci, key=lambda x: x.feature_uid):
            name = f"locus_{_SAFE.sub('_', loc.accession or loc.contig_uid.hex())}_{loc.feature_uid.hex()}.gff3"
            export.write_lines(_gff(loc), gff / name)
            n_gff += 1
            anchor = Gene(0, loc.feature_uid, loc.anchor_protein_uid, loc.start, loc.end, loc.strand, loc.anchor_key,
                          True)
            for g in loc.up + [anchor] + loc.down:
                if g.protein_uid is not None:
                    proteins.setdefault(g.protein_uid, g.key if g.key is not None else g.protein_uid.hex())
        export.write_lines(_tsv(["protein_uid", "group"], ({"protein_uid": u, "group": proteins[u]}
                                                          for u in sorted(proteins))), gff / "protein_groups.tsv")
    return {"groups": len(order), "loci": len(loci), "context_rows": len(ctx_rows), "conservation_rows": len(cons_rows),
            "gff": n_gff, "peak_leaves": peak_leaves, "figures": n_figures,
            "leaves": sum(r["n_leaves"] for r in group_rows)}


def run(root: Path, anchors: set[bytes], out: Path, *, release: str, files, window: int, level: str, group_by: str,
        distance: str, cut: float, max_leaves: int, export_gff: bool, library: str | None, gene_call_version: str,
        tape_version: str | None, threads: int | None = None, tmp=None, figure: bool = False,
        figure_rows: int = DEFAULT_FIGURE_ROWS) -> tuple[dict, list[str]]:
    """Load, group and write into `out` (new or empty, outside the root): filled under a temporary name beside it and
    renamed into place. Returns (counts, notes)."""
    loci, facts = load_loci(root, anchors, window, release=release, files=files, level=level, library=library,
                            gene_call_version=gene_call_version, tape_version=tape_version, threads=threads, tmp=tmp)
    facts.update({"release": release, "level": level, "n_anchor_proteins_given": len(anchors)})
    staging = out.parent / f".{out.name}.tmp-{os.getpid()}"
    shutil.rmtree(staging, ignore_errors=True)
    staging.mkdir(parents=True)
    notes: list[str] = []
    try:
        counts = write(staging, loci, facts, window=window, group_by=group_by, distance=distance, cut=cut,
                       max_leaves=max_leaves, export_gff=export_gff, notes=notes, figure=figure,
                       figure_rows=figure_rows)
        os.replace(staging, out)  # an empty `out` is replaced; one that gained files meanwhile is not (OSError)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    return counts, notes
