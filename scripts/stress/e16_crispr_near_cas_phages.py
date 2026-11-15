#!/usr/bin/env python3
"""Card E-16: the chance rate of the filtered CRISPR producer (DEC-40, card E-15) on phages (stage S1).

Usage:
    python scripts/stress/e16_crispr_near_cas_phages.py measure --root ROOT --calls E13_VIEW --kept E15_VIEW
        --tree TREE --out DIR [--stages S1,P1,P2] [--seq-stages S0,S1] [--draws 10000] [--seed 16] [--e13-extra FILE]

Runs on the pod, in the test or runtime image, one thread; reads only, writes only into `--out` (a new directory).
Every input is read by an explicit path (AGENTS §3.9; never a glob over `l1/feature` or `l1/crispr`):

- the calls of MinCED at 2 repeats (`minced-0.4.2-w8-nr2-r23x47-s26x50`, card E-13) and the default producer's
  arrays (`...-nr3-...`) from E13_VIEW (`/tmp/if2453/e13/view`); the filtered producer's kept calls
  (`...-nr2-...-cas1kb-caspfam1-pfam-38.2`, card E-15) from E15_VIEW (`/tmp/if2453/e15/view`);
- from ROOT: the stage's gene calls (`l1/feature_protein/gene_call_version=pyrodigal-3.7.1-gv-0.3.2-p2`), its Pfam
  scan (`l1/domains/library=pfam-38.2`), `l0seq/genome_input` and `l0seq/contig_index`, and for S1 the record names
  (`l0seq/contig_occurrence/stage=S1/source=inphared`);
- the 31 Cas families with their roles from the worked example's SQL (TREE/examples/crispr_cas_neighbourhoods/sql/
  q01_loci_by_domains.sql, every role whose name starts with `cas`), not from the producer's list; checked equal to
  card E-13's `CAS_FAMILIES`.

The three rules (the producer's, `biobox.l1.crispr` docstring), at a call's span [s, e) on its contig: (a) it shares
no base with an array of the default producer on the contig (d_s < e and s < d_e for none), (b) the gap
max(0, g_s - e, s - g_e) to some Cas gene (a called gene whose protein has a domain of one of the 31 families in the
stage's pfam-38.2 scan) is at most 1,000 bases, (c) no called gene has g_s <= s and e <= g_e. On a circular contig every
gene and default array is also taken shifted by -n and +n. `rules_at` evaluates them directly from the coordinates
(sorted starts with a running maximum of ends: the first rule holds at [s, e) when some gene starting at or before
e + 1,000 ends at or after s - 1,000, and so on); applied to the real calls it must keep exactly the producer's calls
(the first oracle check).

The null (`exact_keep`, `simulate`). Randomized: the start of each call. Held fixed: the call's length L and its
genome (so the number of calls per genome and their lengths), the genome's contigs with their lengths and topology,
its called genes, which of them are Cas genes and of which families, and the default producer's arrays. A call's
start is uniform over every place where it fits in its genome: on a linear contig of n bases the n - L + 1 starts,
on a circular one the n starts (the call runs across the origin when s + L > n), a contig weighted by its number
of starts. At that place the three rules are applied as above. p_i, the probability that call i is kept, is computed
exactly from the geometry (`exact_keep`: in the space of starts, the union of the windows [g_s - L - 1000, g_e + 1000]
of the Cas genes minus the starts that overlap a default array, [d_s - L + 1, d_e - 1], and those inside a gene,
[g_s, g_e - L], counted as integer intervals), and estimated by Monte Carlo (`simulate`: `--draws` places per call,
the rules applied by `rules_at`), so two implementations of one number check each other. Each call has its own
random stream (numpy PCG64, seeded from `--seed` and the call's feature_uid), so the result does not depend on the
order of the calls. A draw is one placement of every call; K_d, the calls kept in draw d, has expectation sum p_i
and is a sum of independent Bernoulli variables (a Poisson-binomial, `tail_at_least`).

Three readings of the null, all reported:
- N2 (the card's): every call at 2 repeats is placed, the three rules applied at the new place.
- N3: only the calls that pass rule (a) at their real place (the "extra calls" of card E-13, which are not a
  default array called again) are placed, the three rules applied at the new place.
- E13: card E-13's quantity behind "at most 2.96 per Gbp" (its rule A5 bound): the sum over the extra calls that are
  not wholly inside a gene at their real place of f, the share of the starts of the call's genome, every contig read
  as linear (n - L + 1 starts), that lie within 1 kb of a Cas gene (`e13_window`). Rules (a) and (c) are not applied
  at the new place, which is why it is an upper bound on what the rule keeps by chance. Reproduced on P1 + P2 (the
  second oracle check) exactly from the geometry and by Monte Carlo, and compared per call with card E-13's stored
  f_1_kb (`--e13-extra`, its `extra_arrays.tsv`).

Contigs a stage lists but does not store (card E-17, DEC-26): a contig that `genome_input` lists and an earlier stage
stores (its default variant is in that stage's `contig_index`) is that stage's: both CRISPR producers and gene calling
skip it in this stage. The null takes the stage's own stored contigs; `--seq-stages S0,S1` (as `bio-box stage crispr`
takes it) names the stages to look in, and the report gives `contigs_listed`, `contigs_left_out` and
`contigs_left_out_by_stage`. A listed contig that no named stage stores is an error.

Also, for scale, the arrays of the default producer (3 repeats) placed in the same way: observed and expected within
1 kb of a Cas gene, and within 10 genes by card E-13's own `Windows` and `placement_share` (its 32.9-fold).

Per Cas family: the expected number of calls kept by chance beside a gene of that family (within 1 kb), the exact
share restricted to the windows of that family's genes (a call beside two families counts for both), with the
observed kept calls beside it.

Writes into `--out`: `summary.json` (every number), `kept_s1.tsv` (the producer's kept calls on S1 with their phage,
Cas genes and group), `families.tsv` and `roles.tsv` (per stage), `calls_<stage>.tsv` (per call: its p_i and
E-13's f), `draws_<stage>.npz` (K_d of N2, N3 and E13 per draw).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import resource
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

GAP = 1000
GENE_CALL = "pyrodigal-3.7.1-gv-0.3.2-p2"
LIBRARY = "pfam-38.2"
NR2 = "minced-0.4.2-w8-nr2-r23x47-s26x50"
NR3 = "minced-0.4.2-w8-nr3-r23x47-s26x50"
NEAR = f"{NR2}-cas1kb-caspfam1-{LIBRARY}"
STRESS = Path(__file__).resolve().parent


# --- one contig's geometry ---------------------------------------------------------------------------------------


def _shifts(n: int, circ: bool) -> tuple[int, ...]:
    return (-n, 0, n) if circ else (0,)


def _sorted_cummax(starts: np.ndarray, ends: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    o = np.argsort(starts, kind="stable")
    s, e = starts[o], ends[o]
    return s, (np.maximum.accumulate(e) if len(e) else e)


@dataclass
class Contig:
    """A contig of n bases: its genes (starts, ends, sorted by start), which are Cas genes and of which families, and
    the default producer's arrays."""

    n: int
    circ: bool
    gs: np.ndarray
    ge: np.ndarray
    cas: np.ndarray  # indices into gs/ge
    cas_fams: list  # per Cas gene, a frozenset of Pfam accessions
    ds: np.ndarray
    de: np.ndarray
    uid: bytes = b""
    _rules: dict = field(default_factory=dict, repr=False)

    def copies(self, s: np.ndarray, e: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        sh = _shifts(self.n, self.circ)
        return (np.concatenate([s + k for k in sh]).astype(np.int64),
                np.concatenate([e + k for k in sh]).astype(np.int64))

    def rule(self, name: str):
        """Sorted starts and running maximum of ends of the copies of: `genes`, `cas`, `default`, or `fam:<PF>`."""
        if name not in self._rules:
            if name == "genes":
                s, e = self.gs, self.ge
            elif name == "cas":
                s, e = self.gs[self.cas], self.ge[self.cas]
            elif name == "default":
                s, e = self.ds, self.de
            else:
                pf = name.split(":", 1)[1]
                idx = np.array([i for i, f in zip(self.cas, self.cas_fams) if pf in f], dtype=np.int64)
                s, e = self.gs[idx], self.ge[idx]
            self._rules[name] = _sorted_cummax(*self.copies(np.asarray(s, np.int64), np.asarray(e, np.int64)))
        return self._rules[name]

    def families(self) -> list[str]:
        return sorted(set().union(*self.cas_fams)) if self.cas_fams else []


def make_contig(n: int, circular: bool, genes=(), cas=None, defaults=(), uid: bytes = b"") -> Contig:
    """A contig from plain lists: genes [(start, end)], cas {gene index: {Pfam accessions}}, defaults [(start, end)]."""
    genes = sorted(genes)
    cas = cas or {}
    gs = np.array([g[0] for g in genes], dtype=np.int64)
    ge = np.array([g[1] for g in genes], dtype=np.int64)
    idx = sorted(cas)
    return Contig(n, circular, gs, ge, np.array(idx, dtype=np.int64), [frozenset(cas[i]) for i in idx],
                  np.array([d[0] for d in defaults], dtype=np.int64), np.array([d[1] for d in defaults], dtype=np.int64),
                  uid)


def _near(c: Contig, name: str, s: np.ndarray, e: np.ndarray, gap: int = GAP) -> np.ndarray:
    st, cm = c.rule(name)
    out = np.zeros(len(s), dtype=bool)
    if not len(st):
        return out
    i = np.searchsorted(st, e + gap, side="right")  # copies with start <= e + gap
    ok = i > 0
    out[ok] = cm[i[ok] - 1] >= s[ok] - gap
    return out


def rules_at(c: Contig, s: np.ndarray, e: np.ndarray, gap: int = GAP):
    """(overlaps a default array, within `gap` of a Cas gene, wholly inside a gene) at spans [s, e), from the
    coordinates."""
    s, e = np.asarray(s, np.int64), np.asarray(e, np.int64)
    ov = np.zeros(len(s), dtype=bool)
    st, cm = c.rule("default")
    if len(st):
        i = np.searchsorted(st, e, side="left")  # copies with start < e
        ok = i > 0
        ov[ok] = cm[i[ok] - 1] > s[ok]
    near = _near(c, "cas", s, e, gap)
    ins = np.zeros(len(s), dtype=bool)
    st, cm = c.rule("genes")
    if len(st):
        i = np.searchsorted(st, s, side="right")  # copies with start <= s
        ok = i > 0
        ins[ok] = cm[i[ok] - 1] >= e[ok]
    return ov, near, ins


# --- exact shares in the space of starts ------------------------------------------------------------------------


def _merge(lo, hi) -> list[tuple[int, int]]:
    """Closed integer intervals [lo, hi] (empty ones dropped), merged and sorted."""
    lo, hi = np.asarray(lo, np.int64), np.asarray(hi, np.int64)
    keep = lo <= hi
    lo, hi = lo[keep], hi[keep]
    if not len(lo):
        return []
    o = np.argsort(lo, kind="stable")
    lo, hi = lo[o].tolist(), hi[o].tolist()
    out, cl, ch = [], lo[0], hi[0]
    for a, b in zip(lo[1:], hi[1:]):
        if a > ch + 1:
            out.append((cl, ch))
            cl, ch = a, b
        elif b > ch:
            ch = b
    out.append((cl, ch))
    return out


def _size(m: list[tuple[int, int]]) -> int:
    return sum(b - a + 1 for a, b in m)


def _inter_size(m1: list[tuple[int, int]], m2: list[tuple[int, int]]) -> int:
    i = j = tot = 0
    while i < len(m1) and j < len(m2):
        lo, hi = max(m1[i][0], m2[j][0]), min(m1[i][1], m2[j][1])
        if lo <= hi:
            tot += hi - lo + 1
        if m1[i][1] < m2[j][1]:
            i += 1
        else:
            j += 1
    return tot


def _domain(c: Contig, length: int, linear: bool = False) -> int:
    if c.circ and not linear:
        return c.n if length <= c.n else 0
    return max(0, c.n - length + 1)


def _clip(lo: np.ndarray, hi: np.ndarray, top: int) -> tuple[np.ndarray, np.ndarray]:
    return np.maximum(lo, 0), np.minimum(hi, top)


def _windows(c: Contig, length: int, top: int, name: str = "cas", gap: int = GAP) -> list[tuple[int, int]]:
    st, _ = c.rule(name)
    if not len(st):
        return []
    # the copies' own ends (not the running maximum): recompute from the copies
    if name == "cas":
        s, e = c.copies(c.gs[c.cas], c.ge[c.cas])
    else:
        pf = name.split(":", 1)[1]
        idx = np.array([i for i, f in zip(c.cas, c.cas_fams) if pf in f], dtype=np.int64)
        s, e = c.copies(c.gs[idx], c.ge[idx])
    return _merge(*_clip(s - length - gap, e + gap, top))


def _removed(c: Contig, length: int, top: int, span: tuple[int, int]) -> list[tuple[int, int]]:
    """Starts in `span` at which a call overlaps a default array or lies wholly inside a gene."""
    a, b = span
    lo, hi = [], []
    if len(c.ds):
        s, e = c.copies(c.ds, c.de)
        lo.append(s - length + 1)
        hi.append(e - 1)
    if len(c.gs):
        s, e = c.copies(c.gs, c.ge)
        lo.append(s)
        hi.append(e - length)
    if not lo:
        return []
    lo, hi = np.concatenate(lo), np.concatenate(hi)
    keep = (hi >= a) & (lo <= b) & (lo <= hi)
    return _merge(*_clip(lo[keep], hi[keep], top))


def exact_keep(genome: list[Contig], length: int, gap: int = GAP) -> dict:
    """The starts of a call of `length` bases in its genome (`domain`), those within `gap` of a Cas gene (`window`),
    those at which the three rules keep it (`keep`), and per family the kept starts within `gap` of a gene of that
    family (`family_keep`)."""
    dom = win = keep = 0
    fam: Counter = Counter()
    for c in genome:
        d = _domain(c, length)
        dom += d
        if not d or not len(c.cas):
            continue
        top = d - 1
        w = _windows(c, length, top, "cas", gap)
        if not w:
            continue
        r = _removed(c, length, top, (w[0][0], w[-1][1]))
        win += _size(w)
        keep += _size(w) - _inter_size(w, r)
        for pf in c.families():
            wf = _windows(c, length, top, f"fam:{pf}", gap)
            k = _size(wf) - _inter_size(wf, r)
            if k:
                fam[pf] += k
    return {"domain": dom, "window": win, "keep": keep, "family_keep": dict(fam)}


def e13_window(genome: list[Contig], length: int, gap: int = GAP) -> tuple[int, int]:
    """Card E-13's window share as (numerator, denominator): every contig read as linear, n - L + 1 starts; the
    starts within `gap` of a Cas gene (its copies at -n, 0, +n on a circular contig), clipped to [0, n - L]."""
    num = den = 0
    for c in genome:
        d = _domain(c, length, linear=True)
        den += d
        if d and len(c.cas):
            num += _size(_windows(c, length, d - 1, "cas", gap))
    return num, den


# --- Monte Carlo ------------------------------------------------------------------------------------------------


def simulate(genome: list[Contig], length: int, draws: int, rng: np.random.Generator, linear: bool = False,
             families: bool = False) -> dict:
    """`draws` uniform places of a call of `length` bases in its genome (every contig linear with `linear`, as card
    E-13 placed them); per place whether it lies within 1 kb of a Cas gene (`near`) and whether the three rules keep
    it (`kept`); with `families`, per family the number of kept places beside a gene of it."""
    w = np.array([_domain(c, length, linear) for c in genome], dtype=np.int64)
    tot = int(w.sum())
    if tot == 0:
        z = np.zeros(draws, dtype=bool)
        return {"near": z, "kept": z.copy(), "family": {}}
    ci = np.zeros(draws, dtype=np.int64) if len(w) == 1 else rng.choice(len(w), size=draws, p=w / tot)
    s = rng.integers(0, w[ci])
    e = s + length
    near = np.zeros(draws, dtype=bool)
    kept = np.zeros(draws, dtype=bool)
    fam: Counter = Counter()
    for k, c in enumerate(genome):
        m = np.flatnonzero(ci == k)
        if not len(m) or not len(c.cas):
            continue
        ov, nr, ins = rules_at(c, s[m], e[m])
        near[m] = nr
        kp = nr & ~ov & ~ins
        kept[m] = kp
        if families and kp.any():
            mk = m[kp]
            for pf in c.families():
                fam[pf] += int(_near(c, f"fam:{pf}", s[mk], e[mk]).sum())
    return {"near": near, "kept": kept, "family": dict(fam)}


def tail_at_least(ps, k: int) -> float:
    """P(K >= k) for K a sum of independent Bernoulli(p) (the Poisson-binomial), by dynamic programming."""
    pmf = np.zeros(1)
    pmf[0] = 1.0
    for p in ps:
        if p <= 0:
            continue
        nxt = np.zeros(len(pmf) + 1)
        nxt[:-1] += pmf * (1 - p)
        nxt[1:] += pmf * p
        pmf = nxt
    return float(pmf[k:].sum()) if k < len(pmf) else 0.0


def call_rng(seed: int, stage: str, feature_uid: bytes, salt: int) -> np.random.Generator:
    words = np.frombuffer(feature_uid, dtype=np.uint32).tolist()
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(
        [seed, salt, *[ord(ch) for ch in stage], *words])))


# --- reading ----------------------------------------------------------------------------------------------------


def _files(d: Path) -> list[Path]:
    """The part files of one partition directory, by explicit name (the directory is never globbed above it)."""
    fs = sorted((p for p in d.iterdir() if re.fullmatch(r"part-\d+\.parquet", p.name)),
                key=lambda p: int(p.stem.split("-")[1]))
    if not fs:
        raise SystemExit(f"{d}: no part files")
    return fs


def _read(d: Path, cols: list[str]):
    import pyarrow as pa
    import pyarrow.parquet as pq

    return pa.concat_tables([pq.read_table(f, columns=cols, partitioning=None) for f in _files(d)])


def cas_families_from_example(tree: Path) -> dict[str, str]:
    sql = (tree / "examples/crispr_cas_neighbourhoods/sql/q01_loci_by_domains.sql").read_text()
    out = {pf: role for pf, role in re.findall(r"\('(PF\d{5})', '([A-Za-z0-9]+)'\)", sql) if role.startswith("cas")}
    if len(out) != 31:
        raise SystemExit(f"{len(out)} Cas families in the example's SQL, not 31")
    return out


def _load_e13():
    spec = importlib.util.spec_from_file_location("e13_crispr_min_repeats", STRESS / "e13_crispr_min_repeats.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@dataclass
class Stage:
    name: str
    contigs: dict  # contig_uid -> Contig
    genome_of: dict  # contig_uid -> genome_uid
    genome_contigs: dict  # genome_uid -> [contig_uid]
    bases: int
    calls: list  # nr2 calls: dicts
    nr3: list
    kept_producer: set
    cas_genes: list  # (contig_uid, start, end, strand, protein_uid, frozenset of PF)
    fam_names: dict
    cas_gene_total: Counter
    cas_genome_total: dict
    listed: int = 0  # contigs of genome_input
    skipped: dict = field(default_factory=dict)  # contig_uid -> the stage that stores a listed contig left out


def load_stage(root: Path, calls_view: Path, kept_view: Path, st: str, fams: dict[str, str],
               seq_stages=()) -> Stage:
    import pyarrow as pa
    import pyarrow.compute as pc

    pa.set_cpu_count(1)
    pa.set_io_thread_count(1)
    gi = _read(root / f"l0seq/genome_input/stage={st}", ["contig_uid", "genome_uid", "is_circular"]).to_pylist()
    genome_of: dict = {}
    circ: dict = {}
    for r in gi:
        cu = r["contig_uid"]
        genome_of[cu] = min(genome_of.get(cu, r["genome_uid"]), r["genome_uid"])  # DEC-27: the least genome
        circ[cu] = circ.get(cu, False) or bool(r["is_circular"])
    ci = _read(root / f"l0seq/contig_index/stage={st}", ["contig_uid", "contig_len", "is_default_variant"]).to_pylist()
    nlen: dict = {}
    for r in ci:
        # the contigs the stage stores: it holds their default variant (DEC-26, `ContigOwners`)
        if r["is_default_variant"] and r["contig_uid"] in genome_of:
            if nlen.setdefault(r["contig_uid"], r["contig_len"]) != r["contig_len"]:
                raise SystemExit(f"{st}: contig {r['contig_uid'].hex()} has two lengths")
    # a listed contig another stage stores is that stage's (DEC-26): left out, as both CRISPR producers and gene
    # calling skip it; one that no named stage stores is an error
    left = sorted(set(genome_of) - set(nlen))
    skipped: dict = {}  # contig_uid -> the stage that stores it
    for other in dict.fromkeys(s for s in seq_stages if s != st):
        todo = [cu for cu in left if cu not in skipped]
        if not todo:
            break
        if not (root / f"l0seq/contig_index/stage={other}").is_dir():
            continue  # a named stage with no contig_index stores nothing
        t = _read(root / f"l0seq/contig_index/stage={other}", ["contig_uid", "is_default_variant"])
        t = t.filter(pc.and_(t.column("is_default_variant"),
                             pc.is_in(t.column("contig_uid"), value_set=pa.array(todo, pa.binary(16)))))
        for cu in t.column("contig_uid").to_pylist():
            skipped.setdefault(cu, other)
    lost = [cu for cu in left if cu not in skipped]
    if lost:
        raise SystemExit(f"{st}: {len(lost)} contigs of genome_input are stored by none of the stages "
                         f"{list(dict.fromkeys([st, *seq_stages]))} (first {lost[0].hex()}); name the stage that "
                         f"stores them with --seq-stages")
    listed = len(genome_of)
    for cu in skipped:
        del genome_of[cu]
        del circ[cu]
    genome_contigs: dict = defaultdict(list)
    for cu in sorted(genome_of):
        genome_contigs[genome_of[cu]].append(cu)

    def feats(view: Path, prod: str) -> list[dict]:
        t = _read(view / f"l1/feature/producer={prod}/stage={st}", ["feature_uid", "contig_uid", "start", "end"])
        return sorted(t.to_pylist(), key=lambda r: (r["contig_uid"], r["start"], r["end"]))

    calls, nr3 = feats(calls_view, NR2), feats(calls_view, NR3)
    kept = {(r["contig_uid"], r["start"], r["end"]) for r in feats(kept_view, NEAR)}
    el = _read(calls_view / f"l1/crispr/producer={NR2}/stage={st}", ["feature_uid", "contig_len"])
    lens = dict(zip(el.column("feature_uid").to_pylist(), el.column("contig_len").to_pylist()))
    for c in calls:
        if lens[c["feature_uid"]] != nlen[c["contig_uid"]]:
            raise SystemExit(f"{st}: a call's contig_len differs from contig_index")

    # Cas proteins: the stage's pfam-38.2 rows of the 31 families
    dom = _read(root / f"l1/domains/library={LIBRARY}/stage={st}", ["protein_uid", "model_acc", "model_name"])
    acc = pc.utf8_slice_codeunits(dom.column("model_acc"), 0, 7)
    m = pc.is_in(acc, value_set=pa.array(sorted(fams)))
    cas_prot: dict = defaultdict(set)
    fam_names: dict = {}
    for pu, a, nm in zip(dom.column("protein_uid").filter(m).to_pylist(), acc.filter(m).to_pylist(),
                         dom.column("model_name").filter(m).to_pylist()):
        cas_prot[pu].add(a)
        fam_names[a] = nm
    del dom, acc, m

    # genes of the contigs whose genome holds a call (nr2 or nr3); Cas genes everywhere (for the family table)
    with_call = {genome_of[c["contig_uid"]] for c in calls + nr3}
    want = sorted(cu for g in with_call for cu in genome_contigs[g])
    want_arr = pa.array(want, pa.binary(16))
    cas_arr = pa.array(sorted(cas_prot), pa.binary(16))
    gene_parts = []
    cas_genes = []
    import pyarrow.parquet as pq

    for f in _files(root / f"l1/feature_protein/gene_call_version={GENE_CALL}/stage={st}"):
        pf = pq.ParquetFile(f)
        for g in range(pf.metadata.num_row_groups):
            t = pf.read_row_group(g, columns=["contig_uid", "start", "end", "strand", "protein_uid"])
            isc = pc.is_in(t.column("protein_uid"), value_set=cas_arr)
            for r in t.filter(isc).to_pylist():
                cas_genes.append((r["contig_uid"], r["start"], r["end"], r["strand"], r["protein_uid"],
                                  frozenset(cas_prot[r["protein_uid"]])))
            idx = pc.index_in(t.column("contig_uid"), value_set=want_arr)
            keep = pc.is_valid(idx)
            gene_parts.append((np.asarray(idx.filter(keep).to_numpy(zero_copy_only=False), dtype=np.int64),
                               t.column("start").filter(keep).to_numpy().astype(np.int64),
                               t.column("end").filter(keep).to_numpy().astype(np.int64),
                               np.asarray(isc.filter(keep).to_numpy(zero_copy_only=False), dtype=bool)))
    gi_ = np.concatenate([p[0] for p in gene_parts]) if gene_parts else np.zeros(0, np.int64)
    gs_ = np.concatenate([p[1] for p in gene_parts]) if gene_parts else np.zeros(0, np.int64)
    ge_ = np.concatenate([p[2] for p in gene_parts]) if gene_parts else np.zeros(0, np.int64)
    del gene_parts
    o = np.lexsort((ge_, gs_, gi_))
    gi_, gs_, ge_ = gi_[o], gs_[o], ge_[o]
    bounds = np.searchsorted(gi_, np.arange(len(want) + 1))
    cas_by_contig: dict = defaultdict(list)
    for cu, s, e, _strand, _pu, fs in cas_genes:
        cas_by_contig[cu].append((s, e, fs))
    d_by: dict = defaultdict(list)
    for a in nr3:
        d_by[a["contig_uid"]].append((a["start"], a["end"]))
    contigs: dict = {}
    for k, cu in enumerate(want):
        gs, ge = gs_[bounds[k]:bounds[k + 1]], ge_[bounds[k]:bounds[k + 1]]
        cas_spans = defaultdict(set)
        for s, e, fs in cas_by_contig.get(cu, []):
            cas_spans[(s, e)] |= fs
        idx, fl = [], []
        for i, (s, e) in enumerate(zip(gs.tolist(), ge.tolist())):
            if (s, e) in cas_spans:
                idx.append(i)
                fl.append(frozenset(cas_spans[(s, e)]))
        if len({(s, e) for s, e in zip(gs.tolist(), ge.tolist()) if (s, e) in cas_spans}) != len(cas_spans):
            raise SystemExit(f"{st}: a Cas gene of {cu.hex()} is not among its genes")
        ds = np.array([d[0] for d in d_by.get(cu, [])], dtype=np.int64)
        de = np.array([d[1] for d in d_by.get(cu, [])], dtype=np.int64)
        contigs[cu] = Contig(int(nlen[cu]), bool(circ[cu]), gs, ge, np.array(idx, dtype=np.int64), fl, ds, de, cu)
    tot: Counter = Counter()
    gen: dict = defaultdict(set)
    for cu, _s, _e, _strand, _pu, fs in cas_genes:
        for pf in fs:
            tot[pf] += 1
            gen[pf].add(genome_of.get(cu))
    return Stage(st, contigs, genome_of, dict(genome_contigs), int(sum(nlen.values())), calls, nr3, kept, cas_genes,
                 fam_names, tot, {k: len(v) for k, v in gen.items()}, listed, skipped)


# --- the measurement ----------------------------------------------------------------------------------------------


def _ci(x: np.ndarray) -> dict:
    m, sd = float(x.mean()), float(x.std(ddof=1))
    se = sd / len(x) ** 0.5
    q = np.percentile(x, [2.5, 50, 97.5], method="lower")
    return {"mean": round(m, 4), "sd": round(sd, 4), "se": round(se, 5),
            "mean_ci95": [round(m - 1.96 * se, 4), round(m + 1.96 * se, 4)],
            "draws_p2_5": int(q[0]), "draws_median": int(q[1]), "draws_p97_5": int(q[2]), "draws_max": int(x.max())}


def measure_stage(stg: Stage, draws: int, seed: int, fams: dict[str, str], e13mod) -> tuple[dict, list[dict], dict]:
    st = stg.name
    genome = {g: [stg.contigs[cu] for cu in cus if cu in stg.contigs] for g, cus in stg.genome_contigs.items()}
    rows = []
    tot = {"N2": np.zeros(draws, np.int64), "N3": np.zeros(draws, np.int64), "E13": np.zeros(draws, np.int64)}
    fam_exact: dict = {"N2": Counter(), "N3": Counter()}
    fam_mc: Counter = Counter()
    mine_kept = set()
    t0 = time.process_time()
    for k, c in enumerate(stg.calls):
        cu, s, e = c["contig_uid"], c["start"], c["end"]
        L = e - s
        cg = stg.contigs[cu]
        g = genome[stg.genome_of[cu]]
        ov, nr, ins = (x[0] for x in rules_at(cg, np.array([s]), np.array([e])))
        if not ov and nr and not ins:
            mine_kept.add((cu, s, e))
        # E-13's "extra" (overlap with an array at 3, linear, as its `_overlaps`) and "inside a gene" (with copies)
        e13_extra = not any(d0 < e and s < d1 for d0, d1 in zip(cg.ds.tolist(), cg.de.tolist()))
        ex = exact_keep(g, L)
        num, den = e13_window(g, L)
        p = ex["keep"] / ex["domain"] if ex["domain"] else 0.0
        f = num / den if den else 0.0
        in_e13 = e13_extra and not ins
        rng = call_rng(seed, st, c["feature_uid"], 0)
        any_cas = any(len(x.cas) for x in g)
        if any_cas:
            sim = simulate(g, L, draws, rng, families=True)
            kd = sim["kept"]
            tot["N2"] += kd
            if not ov:
                tot["N3"] += kd
            for pf, v in sim["family"].items():
                fam_mc[pf] += v
            if in_e13:
                tot["E13"] += simulate(g, L, draws, call_rng(seed, st, c["feature_uid"], 1), linear=True)["near"]
        for pf, v in ex["family_keep"].items():
            fam_exact["N2"][pf] += v / ex["domain"]
            if not ov:
                fam_exact["N3"][pf] += v / ex["domain"]
        rows.append({"stage": st, "feature_uid": c["feature_uid"].hex(), "contig_uid": cu.hex(), "start": s, "end": e,
                     "length": L, "overlaps_default": bool(ov), "near_1_kb": bool(nr), "inside_gene": bool(ins),
                     "kept_rules": bool(not ov and nr and not ins), "e13_extra": e13_extra, "in_e13_set": in_e13,
                     "genome_has_cas": any_cas, "p_keep": p, "f_e13_1kb": f, "keep_starts": ex["keep"],
                     "window_starts": ex["window"], "domain_starts": ex["domain"], "e13_num": num, "e13_den": den})
        if k and k % 5000 == 0:
            print(f"{st}: {k} of {len(stg.calls)} calls, {time.process_time() - t0:.0f} CPU-s", file=sys.stderr,
                  flush=True)
    gbp = stg.bases / 1e9
    ps_all = [r["p_keep"] for r in rows]
    ps_n3 = [r["p_keep"] for r in rows if not r["overlaps_default"]]
    obs = len(stg.kept_producer)
    res = {"stage": st, "bases": stg.bases, "gbp": gbp, "calls_nr2": len(stg.calls), "arrays_nr3": len(stg.nr3),
           "contigs_listed": stg.listed, "contigs_left_out": len(stg.skipped),
           "contigs_left_out_by_stage": dict(sorted(Counter(stg.skipped.values()).items())),
           "contigs": len(stg.genome_of), "genomes": len(stg.genome_contigs),
           "circular_contigs": sum(1 for x in stg.contigs.values() if x.circ),
           "cas_genes": len(stg.cas_genes), "genomes_with_cas_gene": len({stg.genome_of.get(x[0]) for x in stg.cas_genes}),
           "kept_producer": obs, "kept_rules_at_real_place": len(mine_kept),
           "oracle_kept_equal": mine_kept == stg.kept_producer,
           "only_mine": sorted([cu.hex(), s, e] for cu, s, e in mine_kept - stg.kept_producer),
           "only_producer": sorted([cu.hex(), s, e] for cu, s, e in stg.kept_producer - mine_kept),
           "rejected_rule_a": sum(r["overlaps_default"] for r in rows),
           "rejected_rule_b": sum(not r["overlaps_default"] and not r["near_1_kb"] for r in rows),
           "rejected_rule_c": sum(not r["overlaps_default"] and r["near_1_kb"] and r["inside_gene"] for r in rows),
           "e13_extra": sum(r["e13_extra"] for r in rows), "e13_set": sum(r["in_e13_set"] for r in rows),
           "calls_in_genomes_with_cas": sum(r["genome_has_cas"] for r in rows), "draws": draws, "seed": seed}
    nulls = {}
    for name, ps, arr in (("N2", ps_all, tot["N2"]), ("N3", ps_n3, tot["N3"]),
                          ("E13", [r["f_e13_1kb"] for r in rows if r["in_e13_set"]], tot["E13"])):
        exp = float(sum(ps))
        nulls[name] = {"calls": len(ps), "expected_exact": round(exp, 4),
                       "sd_exact": round(float(np.sqrt(sum(p * (1 - p) for p in ps))), 4),
                       "expected_per_gbp": round(exp / gbp, 4), "monte_carlo": _ci(arr),
                       "monte_carlo_per_gbp_ci95": [round(x / gbp, 4) for x in _ci(arr)["mean_ci95"]]}
        if name != "E13":
            nulls[name]["observed_kept"] = obs
            nulls[name]["observed_per_gbp"] = round(obs / gbp, 4)
            nulls[name]["enrichment"] = round(obs / exp, 3) if exp else None
            nulls[name]["excess"] = round(obs - exp, 3)
            nulls[name]["p_at_least_observed_exact"] = tail_at_least(ps, obs)
            nulls[name]["p_at_least_observed_mc"] = float((arr >= obs).mean())
    res["nulls"] = nulls

    # the default producer's arrays placed in the same way, for scale
    obs1 = obs10 = 0
    exp1 = exp1_e13 = exp10 = 0.0
    e13w: dict = {}
    for a in stg.nr3:
        cu, s, e = a["contig_uid"], a["start"], a["end"]
        cg = stg.contigs[cu]
        g = genome[stg.genome_of[cu]]
        _, nr, _ = rules_at(cg, np.array([s]), np.array([e]))
        obs1 += bool(nr[0])
        num, den = e13_window(g, e - s)
        exp1 += num / den if den else 0.0
        win = []
        for x in g:
            if len(x.cas):
                if x.uid not in e13w:
                    e13w[x.uid] = (e13mod.Windows(x.gs, x.ge, x.n, x.circ, x.cas),
                                   e13mod.Windows(x.gs, x.ge, x.n, x.circ, x.cas, k=1))
                win.append((x.n, *e13w[x.uid]))
            else:
                win.append((x.n, None, None))
        f = e13mod.placement_share(e - s, win)
        exp10 += f["10_genes"]
        exp1_e13 += f["1_kb"]
        if len(cg.cas):
            obs10 += e13w[cu][0].near_genes(s, e)
    res["default_arrays"] = {"arrays": len(stg.nr3), "within_1_kb_observed": obs1,
                             "within_1_kb_expected": round(exp1, 3), "within_1_kb_expected_e13_function": round(exp1_e13, 3),
                             "within_1_kb_enrichment": round(obs1 / exp1, 2) if exp1 else None,
                             "within_10_genes_observed": int(obs10), "within_10_genes_expected": round(exp10, 3),
                             "within_10_genes_enrichment": round(obs10 / exp10, 2) if exp10 else None}
    fam_rows = []
    kept_beside: Counter = Counter()
    for cu, s, e in stg.kept_producer:
        cg = stg.contigs[cu]
        for pf in cg.families():
            if _near(cg, f"fam:{pf}", np.array([s]), np.array([e]))[0]:
                kept_beside[pf] += 1
    for pf in sorted(fams):
        fam_rows.append({"stage": st, "family": pf, "name": stg.fam_names.get(pf, ""), "role": fams[pf],
                         "cas_genes": stg.cas_gene_total.get(pf, 0), "genomes": stg.cas_genome_total.get(pf, 0),
                         "expected_kept_by_chance_N2": round(fam_exact["N2"].get(pf, 0.0), 4),
                         "expected_kept_by_chance_N3": round(fam_exact["N3"].get(pf, 0.0), 4),
                         "monte_carlo_N2": round(fam_mc.get(pf, 0) / draws, 4),
                         "observed_kept_beside": kept_beside.get(pf, 0)})
    res["cpu_s"] = round(time.process_time() - t0, 1)
    return res, rows, {"families": fam_rows, "draws": tot}


def _consensus(reps: list[str]) -> str:
    if len({len(r) for r in reps}) != 1:
        return Counter(reps).most_common(1)[0][0]
    out = []
    for col in zip(*reps):
        cnt = Counter(col)
        best = max(cnt.values())
        out.append(next(ch for ch in col if cnt[ch] == best))  # a tie goes to the earliest repeat
    return "".join(out)


HOST_RE = re.compile(r"^\S+\s+(?:UNVERIFIED:\s+)?(.+?)\s+(?:phage|virus|bacteriophage)\b", re.I)


def list_kept(root: Path, kept_view: Path, stg: Stage, fams: dict[str, str], e13mod) -> list[dict]:
    """The producer's kept calls of the stage with their contig, record name, elements and Cas genes."""
    import pyarrow as pa
    import pyarrow.compute as pc

    kept = sorted(stg.kept_producer)
    el = _read(kept_view / f"l1/crispr/producer={NEAR}/stage={stg.name}",
               ["feature_uid", "contig_uid", "element", "ordinal", "start", "end", "seq"]).to_pylist()
    f = _read(kept_view / f"l1/feature/producer={NEAR}/stage={stg.name}",
              ["feature_uid", "contig_uid", "start", "end"]).to_pylist()
    fu_of = {(r["contig_uid"], r["start"], r["end"]): r["feature_uid"] for r in f}
    els: dict = defaultdict(list)
    for r in el:
        els[r["feature_uid"]].append(r)
    names: dict = {}
    if stg.name == "S1":
        co = _read(root / "l0seq/contig_occurrence/stage=S1/source=inphared",
                   ["contig_uid", "l0_object_sha256", "record_ordinal", "accession", "original_header"])
        co = co.filter(pc.is_in(co.column("contig_uid"), value_set=pa.array(sorted({k[0] for k in kept}),
                                                                           pa.binary(16))))
        for r in sorted(co.to_pylist(), key=lambda r: (r["l0_object_sha256"], r["record_ordinal"])):
            names.setdefault(r["contig_uid"], (r["accession"], r["original_header"]))
    cas_detail: dict = defaultdict(list)
    for cu, s, e, strand, pu, fs in stg.cas_genes:
        cas_detail[cu].append((s, e, strand, pu, fs))
    out = []
    for cu, s, e in kept:
        cg = stg.contigs[cu]
        reps = [r["seq"] for r in sorted(els[fu_of[(cu, s, e)]], key=lambda r: r["ordinal"]) if r["element"] == "repeat"]
        nsp = sum(1 for r in els[fu_of[(cu, s, e)]] if r["element"] == "spacer")
        best = None
        beside, within10 = [], []
        w10 = e13mod.Windows(cg.gs, cg.ge, cg.n, cg.circ, cg.cas)
        for j, (gs, ge, strand, pu, fs) in enumerate(sorted(cas_detail[cu])):
            gap = min(max(0, gs + k - e, s - (ge + k)) for k in _shifts(cg.n, cg.circ))
            if best is None or gap < best[0]:
                best = (gap, gs, ge, strand, pu, fs)
            if gap <= GAP:
                beside.append((gs, ge, fs))
            gi = int(np.flatnonzero((cg.gs == gs) & (cg.ge == ge))[0])
            one = e13mod.Windows(cg.gs, cg.ge, cg.n, cg.circ, np.array([gi]))
            if one.near_genes(s, e):
                within10.append((gs, ge, fs))
        roles10 = sorted({fams[pf] for _, _, fs in within10 for pf in fs})
        genes10 = len({(a, b) for a, b, _ in within10})
        roles_beside = sorted({fams[pf] for _, _, fs in beside for pf in fs})
        if len(roles10) >= 2 and genes10 >= 2:
            group = "a"
        elif roles_beside and set(roles_beside) <= {"cas3", "cas4"}:
            group = "b"
        else:
            group = "c"
        acc, head = names.get(cu, ("", ""))
        m = HOST_RE.match(head or "")
        gap, gs, ge, strand, pu, fs = best
        dflt = [(d0, d1) for d0, d1 in zip(cg.ds.tolist(), cg.de.tolist())]
        out.append({"stage": stg.name, "accession": acc, "header": head, "host_or_name": m.group(1) if m else "",
                    "contig_uid": cu.hex(), "contig_len": cg.n, "topology": "circular" if cg.circ else "linear",
                    "start": s, "end": e, "repeats": len(reps), "spacers": nsp, "repeat_len": len(reps[0]),
                    "consensus": _consensus(reps),
                    "nearest_cas_families": ",".join(sorted(fs)),
                    "nearest_cas_names": ",".join(stg.fam_names.get(pf, "") for pf in sorted(fs)),
                    "nearest_cas_role": ",".join(sorted({fams[pf] for pf in fs})),
                    "nearest_cas_protein_uid": pu.hex(), "nearest_cas_gene": f"{gs}-{ge}",
                    "nearest_cas_strand": strand, "gap_bp": gap,
                    "cas_roles_within_1kb": ",".join(roles_beside),
                    "cas_roles_within_10_genes": ",".join(roles10), "cas_genes_within_10_genes": genes10,
                    "group": group, "default_arrays_on_contig": len(dflt),
                    "icp1": "ICP1" in (head or "")})
    out.sort(key=lambda r: (r["group"], r["accession"], r["start"]))
    return out


def _write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    cols = list(rows[0])

    def cell(v):
        if isinstance(v, bool):
            return "1" if v else "0"
        if isinstance(v, float):
            return f"{v:.6g}"
        return str(v)

    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(cell(r[c]) for c in cols) + "\n")


def compare_e13(rows: list[dict], path: Path) -> dict:
    """Per call, this script's E-13 set and f against card E-13's `extra_arrays.tsv` (its P1 and P2 extra calls)."""
    theirs = {}
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            if r["inside_gene"] == "0":
                theirs[(r["stage"], r["contig_uid"], int(r["start"]), int(r["end"]))] = float(r["f_1_kb"])
    mine = {(r["stage"], r["contig_uid"], r["start"], r["end"]): r["f_e13_1kb"] for r in rows if r["in_e13_set"]}
    both = set(mine) & set(theirs)
    diffs = [abs(round(mine[k], 6) - theirs[k]) for k in both]
    return {"e13_rows_not_inside_gene": len(theirs), "mine": len(mine), "both": len(both),
            "only_mine": len(set(mine) - set(theirs)), "only_e13": len(set(theirs) - set(mine)),
            "f_equal_at_6_decimals": sum(d < 5e-7 for d in diffs), "f_max_abs_diff": max(diffs) if diffs else None,
            "sum_f_mine": round(sum(mine.values()), 4), "sum_f_e13_rounded": round(sum(theirs.values()), 4)}


def cmd_measure(a) -> int:
    t0, c0 = time.time(), time.process_time()
    root, out, tree = Path(a.root), Path(a.out), Path(a.tree)
    out.mkdir(parents=True, exist_ok=False)
    fams = cas_families_from_example(tree)
    e13 = _load_e13()
    if fams != e13.CAS_FAMILIES:
        raise SystemExit("the example's Cas families differ from card E-13's CAS_FAMILIES")
    summary: dict = {"families": len(fams), "gap": GAP, "producers": {"nr2": NR2, "nr3": NR3, "near": NEAR},
                     "gene_call_version": GENE_CALL, "library": LIBRARY, "stages": {}}
    seq = [x for x in a.seq_stages.split(",") if x]
    all_rows, fam_rows = [], []
    for st in a.stages.split(","):
        ts = time.time()
        stg = load_stage(root, Path(a.calls), Path(a.kept), st, fams, seq)
        print(f"{st}: loaded in {time.time() - ts:.0f} s", file=sys.stderr, flush=True)
        res, rows, extra = measure_stage(stg, a.draws, a.seed, fams, e13)
        summary["stages"][st] = res
        print(json.dumps({st: {k: v for k, v in res.items() if k in ("kept_producer", "kept_rules_at_real_place",
                                                                    "oracle_kept_equal", "nulls")}}),
              file=sys.stderr, flush=True)
        _write_tsv(out / f"calls_{st}.tsv", rows)
        np.savez_compressed(out / f"draws_{st}.npz", **extra["draws"])
        fam_rows += extra["families"]
        all_rows += rows
        if st == "S1":
            kept = list_kept(root, Path(a.kept), stg, fams, e13)
            _write_tsv(out / "kept_s1.tsv", kept)
            summary["kept_s1_groups"] = dict(sorted(Counter(r["group"] for r in kept).items()))
            summary["kept_s1_icp1"] = sum(r["icp1"] for r in kept)
            summary["kept_s1_by_nearest_role"] = dict(sorted(Counter(r["nearest_cas_role"] for r in kept).items()))
        del stg
    prok = [s for s in ("P1", "P2") if s in summary["stages"]]
    if prok:
        sub = [summary["stages"][s] for s in prok]
        gbp = sum(x["gbp"] for x in sub)
        comb = {"gbp": gbp, "kept_producer": sum(x["kept_producer"] for x in sub)}
        for nm in ("N2", "N3", "E13"):
            exp = sum(x["nulls"][nm]["expected_exact"] for x in sub)
            dr = sum(np.load(out / f"draws_{s}.npz")[nm] for s in prok)
            comb[nm] = {"expected_exact": round(exp, 4), "expected_per_gbp": round(exp / gbp, 4),
                        "monte_carlo": _ci(dr), "monte_carlo_per_gbp_ci95": [round(x / gbp, 4) for x in
                                                                             _ci(dr)["mean_ci95"]]}
        summary["p1_p2"] = comb
        if a.e13_extra:
            summary["p1_p2"]["e13_per_call"] = compare_e13([r for r in all_rows if r["stage"] in prok],
                                                           Path(a.e13_extra))
    roles: dict = {}
    for r in fam_rows:
        k = (r["stage"], r["role"])
        x = roles.setdefault(k, {"stage": r["stage"], "role": r["role"], "families": 0, "cas_genes": 0,
                                 "expected_kept_by_chance_N2_sum_of_families": 0.0, "observed_kept_beside_sum": 0})
        x["families"] += 1
        x["cas_genes"] += r["cas_genes"]
        x["expected_kept_by_chance_N2_sum_of_families"] = round(
            x["expected_kept_by_chance_N2_sum_of_families"] + r["expected_kept_by_chance_N2"], 4)
        x["observed_kept_beside_sum"] += r["observed_kept_beside"]
    _write_tsv(out / "families.tsv", fam_rows)
    _write_tsv(out / "roles.tsv", list(roles.values()))
    summary["cost"] = {"wall_s": round(time.time() - t0, 1), "cpu_s": round(time.process_time() - c0, 1),
                       "peak_rss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)}
    (out / "summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True, default=str) + "\n")
    print(json.dumps(summary["cost"]))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("measure", help="the null, the oracle checks and the tables")
    m.add_argument("--root", required=True)
    m.add_argument("--calls", required=True, help="card E-13's view: the calls at 2 and the arrays at 3")
    m.add_argument("--kept", required=True, help="card E-15's view: the filtered producer's kept calls")
    m.add_argument("--tree", required=True, help="a checkout (the example's SQL and card E-13's script)")
    m.add_argument("--out", required=True)
    m.add_argument("--stages", default="S1,P1,P2")
    m.add_argument("--seq-stages", default="",
                   help="S0,S1: the stages whose sequences a stage's genomes may name, as the producers take it "
                        "(bio-box stage crispr --seq-stages): a contig that genome_input lists and one of them stores "
                        "is left out of the null and counted (DEC-26); a listed contig that none stores is an error "
                        "(default: none)")
    m.add_argument("--draws", type=int, default=10_000)
    m.add_argument("--seed", type=int, default=16)
    m.add_argument("--e13-extra", help="card E-13's extra_arrays.tsv, to compare f per call")
    a = ap.parse_args(argv)
    return {"measure": cmd_measure}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
