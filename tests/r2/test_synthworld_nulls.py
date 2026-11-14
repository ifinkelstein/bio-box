"""Ring R2: synthworld's null worlds on the `micro` preset (synthworld_spec §4.1, §4.3, §6.1, §6.2; card F-02 session
9; MR-18, B-NULL).

A module fixture builds, each in a child process of its own and several at once: the signal world (S1-S11) at the
default seed; N0-full and N0-strict at the default seed twice (another PYTHONHASHSEED and thread count) and once at
another seed; then N1 of the signal world twice and once with another null seed. Every child reports its wall time
and peak RSS (VmHWM). The oracles are written here: presets are compared leaf by leaf, files by sha256, phyletic
profiles and transfer counts are counted from `assembly/genes.tsv` and `backbone/events.tsv`, adjacency and flank
contexts are read off the emitted truth genes sorted along each emitted contig, the naive co-occurrence test of §4.3 is
a hypergeometric tail summed here from log-binomials, and permutation and binomial bounds are computed here. What they
take from the generator are model inputs and its truth tables, never the statistic under test.
"""

from __future__ import annotations

import collections
import csv
import hashlib
import importlib.util
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pytest

from biobox.synthworld import backbone, emit, nulls

pytestmark = [pytest.mark.r2, pytest.mark.covers("F-02", "B-NULL")]
needs_zombi2 = pytest.mark.skipif(
    importlib.util.find_spec("zombi2") is None,
    reason="zombi2 0.54.0 is not importable here: CI installs it by `uv sync`; on the pod the image has it")

#: the signal world of the micro default seed, every stage (measured 2026-11-13, sessions 1 to 8): unchanged
SIGNAL = {"backbone": "025fdd0d27fdc65e678782216ecd21fb", "sequences": "34ebce1324413a261acd874d7ee4fe3e",
          "replicons": "55bda97155d75148e5cdca7b90a1d7be", "mobile": "0bd7f0db5691de58ea5e1a0df1e6e570",
          "fusions": "db09abf3fe79ae34a73e3ab0664a9b13", "assembly": "2f425158dc38f2302c1d02e9ebed30da",
          "sources": "770007bde8615f41d9dd5f2eb912e753", "emit": "a2f0e448484dab63d41281a5ef6a00c0"}
SIGNAL_WORLD = "f0262027e4ccdcbc72cd7d90f23fbb5c"
STAGES = ("backbone", "sequences", "replicons", "mobile", "fusions", "assembly", "sources", "emit")
SEED = 20261112
SEED2 = SEED + 1
#: worlds built at once; each child is one process (numpy at its thread count)
JOBS = max(1, min(7, os.cpu_count() or 1))

CHILD = r'''
import json, os, sys, time
from biobox.synthworld import backbone, nulls
kind, seed, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
t0 = time.perf_counter()
p = backbone.load_preset("micro")
if kind == "signal":
    ew = nulls.build_stages(p, seed, out)
elif kind.startswith("N0"):
    ew = nulls.build(kind, p, seed, out)
else:
    _, sig, ns = kind.split("|")
    ew = nulls.build_n1(sig, out, null_seed=int(ns))
dt = time.perf_counter() - t0
hwm = -1
if os.path.exists("/proc/self/status"):
    hwm = int(next(l.split()[1] for l in open("/proc/self/status") if l.startswith("VmHWM")))
print(json.dumps({"seconds": round(dt, 2), "vmhwm_kb": hwm, "world_digest": ew.world_digest}))
'''


def _env(variant: str) -> dict:
    threads = "1" if variant == "a" else "4"
    return dict(os.environ, PYTHONHASHSEED="919" if variant == "a" else "0", OMP_NUM_THREADS=threads,
                OPENBLAS_NUM_THREADS=threads, MKL_NUM_THREADS=threads, NUMEXPR_NUM_THREADS=threads)


def _run_jobs(jobs, cwd: Path) -> dict:
    """Run [(name, kind, seed, out, variant)] as child processes, at most JOBS at once; their reports by name."""
    pending, running, done = list(jobs), [], {}
    while pending or running:
        while pending and len(running) < JOBS:
            name, kind, seed, out, variant = pending.pop(0)
            p = subprocess.Popen([sys.executable, "-c", CHILD, kind, str(seed), str(out)], env=_env(variant),
                                 cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            running.append((name, p))
        name, p = running.pop(0)
        out_, err = p.communicate()
        if p.returncode:
            raise AssertionError(f"{name} failed:\n{err[-4000:]}")
        done[name] = json.loads(out_.strip().splitlines()[-1])
    return done


@pytest.fixture(scope="module")
def worlds(tmp_path_factory):
    if importlib.util.find_spec("zombi2") is None:
        pytest.skip("zombi2 0.54.0 is not importable here")
    base = tmp_path_factory.mktemp("nulls")
    W = {}
    jobs = []
    for name, kind, seed, variant in (("signal", "signal", SEED, "a"),
                                      ("n0full", "N0-full", SEED, "a"), ("n0full_b", "N0-full", SEED, "b"),
                                      ("n0full_s2", "N0-full", SEED2, "a"),
                                      ("n0strict", "N0-strict", SEED, "a"), ("n0strict_b", "N0-strict", SEED, "b"),
                                      ("n0strict_s2", "N0-strict", SEED2, "a")):
        W[name] = base / name
        jobs.append((name, kind, seed, W[name], variant))
    t0 = time.perf_counter()
    reports = _run_jobs(jobs, base)
    for name, ns, variant in (("n1", SEED, "a"), ("n1_b", SEED, "b"), ("n1_s2", SEED2, "a")):
        W[name] = base / name
    reports.update(_run_jobs([(name, f"N1|{W['signal']}|{ns}", SEED, W[name], variant)
                              for name, ns, variant in (("n1", SEED, "a"), ("n1_b", SEED, "b"),
                                                        ("n1_s2", SEED2, "a"))], base))
    reports["_wall"] = round(time.perf_counter() - t0, 1)
    return W, reports


def _tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _tree_sha(root: Path) -> dict[str, str]:
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file()}


def _digests(w: Path) -> dict[str, str]:
    man = json.loads((w / "manifest.json").read_text())
    return {d: v["digest"] for d, v in man["stages"].items()}


def _leaf_diff(a, b, where=()):
    out = set()
    if isinstance(a, dict) and isinstance(b, dict):
        for k in set(a) | set(b):
            out |= {where + (k,)} if (k not in a or k not in b) else _leaf_diff(a[k], b[k], where + (k,))
    elif a != b:
        out.add(where)
    return out


def _truth(w: Path, name: str) -> list[dict[str, str]]:
    return _tsv(w / "emit" / "truth" / f"{name}.tsv")


def _preset(w: Path) -> dict:
    return json.loads((w / "backbone" / "backbone.json").read_text())["preset"]


def _presence(w: Path) -> dict[str, frozenset]:
    """Per family, the host genomes holding a copy in the S8 genomes (a fused gene counts for both families)."""
    pres = collections.defaultdict(set)
    for r in _tsv(w / "assembly" / "genes.tsv"):
        for f in r["family_id"].split("+"):
            pres[f].add(r["genome_id"])
    return {f: frozenset(s) for f, s in pres.items()}


def _transfers(w: Path) -> collections.Counter:
    out = collections.Counter()
    for r in _tsv(w / "backbone" / "events.tsv"):
        if r["kind"].startswith("transfer"):
            out[f"F{int(r['family']):06d}"] += 1
    return out


def _adjacency(w: Path) -> tuple[dict, dict]:
    """Per unordered family pair: the genomes where the two are neighbours on an emitted assembly contig, and the
    distinct (left flank, right flank) families around the pair (None at a contig end), read in the pair's order."""
    by = collections.defaultdict(list)
    for r in _truth(w, "truth_gene"):
        if r["occurrence"] == "assembly":
            by[r["contig_uid"]].append((int(r["start"]), int(r["end"]), r["family_id"], r["genome_id"]))
    genomes, contexts = collections.defaultdict(set), collections.defaultdict(set)
    for gs in by.values():
        gs.sort()
        f = [x[2] for x in gs]
        for i in range(len(f) - 1):
            x = f[i - 1] if i > 0 else None
            y = f[i + 2] if i + 2 < len(f) else None
            key, ctx = ((f[i], f[i + 1]), (x, y)) if f[i] <= f[i + 1] else ((f[i + 1], f[i]), (y, x))
            genomes[key].add(gs[i][3])
            contexts[key].add(ctx)
    return genomes, contexts


def _linked_sets(w: Path) -> list[set]:
    """The world's own linked structures other than planted modules: each vehicle, each plasmid and each CRISPR-Cas
    locus's Cas families (vehicles of their own, which the traps score)."""
    fid = {int(r["family"]): r["family_id"] for r in _tsv(w / "backbone" / "families.tsv")}
    sets = [{fid[int(x)] for x in g["members"].split(";")} for g in _tsv(w / "backbone" / "groups.tsv")
            if g["kind"] in ("vehicle", "plasmid")]
    cas = collections.defaultdict(set)
    for r in _tsv(w / "mobile" / "cas.tsv"):
        cas[r["cas_module_id"]].add(r["family_id"])
    return sets + list(cas.values())


def _modules(w: Path, kind: str = "planted") -> dict[str, list[str]]:
    out = collections.defaultdict(list)
    for r in _truth(w, "truth_module"):
        if r["kind"] == kind and r["member_role"] == "core":
            out[r["module_id"]].append(r["family_id"])
    return dict(out)


def _within(pair, sets) -> bool:
    return any(pair[0] in s and pair[1] in s for s in sets)


def _binom_sf(k: int, n: int, p: float) -> float:
    """P(X >= k), X ~ Binomial(n, p), summed exactly."""
    return float(sum(math.comb(n, j) * p ** j * (1 - p) ** (n - j) for j in range(k, n + 1)))


def _excess_test(k0: int, n0: int, k1: int, n1: int) -> float:
    """One-sided exact conditional test that rate k0/n0 exceeds k1/n1: given k0 + k1, k0 ~ Binomial(k0 + k1,
    n0 / (n0 + n1)) when the two rates are equal; the p-value P(X >= k0)."""
    return _binom_sf(k0, k0 + k1, n0 / (n0 + n1))


def _codir_classes(w: Path) -> dict[tuple[str, str], str]:
    """The class of each co-directional junction drawn by the background rule (rule ``drawn``), keyed by its two
    families in transcription order; a key met twice must have one class."""
    g = {r["copy_label"]: (r["family_id"], r["strand"]) for r in _tsv(w / "replicons" / "genes.tsv")}
    out = {}
    for j in _tsv(w / "replicons" / "junctions.tsv"):
        if j["rule"] != "drawn" or not j["junction"].startswith("codir"):
            continue
        a, b = g[j["left_copy"]], g[j["right_copy"]]
        key = (a[0], b[0]) if a[1] == "1" else (b[0], a[0])
        assert out.setdefault(key, j["junction"]) == j["junction"], key
    return out


def _host_kind(w: Path) -> dict[str, str]:
    return {r["family_id"]: r["kind"] for r in _tsv(w / "backbone" / "families.tsv")}


def _hypergeom_sf_table(G: int) -> np.ndarray:
    """T[K, n, k] = P(X >= k), X hypergeometric: k of n draws from G genomes of which K carry the first family."""
    def lc(a, b):
        return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)
    T = np.ones((G + 1, G + 1, G + 2))
    for K in range(G + 1):
        for n in range(G + 1):
            pm = np.zeros(G + 2)
            for k in range(max(0, n - (G - K)), min(K, n) + 1):
                pm[k] = math.exp(lc(K, k) + lc(G - K, n - k) - lc(G, n))
            T[K, n, :] = np.cumsum(pm[::-1])[::-1]
    return T


# --- the signal world does not change; a null world is a world of its own --------------------------------------------

@needs_zombi2
def test_signal_world_keeps_every_stage_digest(worlds):
    W, rep = worlds
    assert _digests(W["signal"]) == SIGNAL
    assert rep["signal"]["world_digest"] == SIGNAL_WORLD
    man = json.loads((W["signal"] / "manifest.json").read_text())
    assert man["world_id"] == f"micro-{SEED}" and "null" not in man
    assert not (W["signal"] / "emit" / "truth" / "truth_null.tsv").exists()


@needs_zombi2
@pytest.mark.covers("F-02", "B-NULL", "MR-16")
def test_every_null_world_is_byte_identical_in_two_processes_with_other_hash_seeds_and_threads(worlds):
    from blake3 import blake3

    W, rep = worlds
    for a, b in (("n0full", "n0full_b"), ("n0strict", "n0strict_b"), ("n1", "n1_b")):
        assert _tree_sha(W[a]) == _tree_sha(W[b]), a
        assert rep[a]["world_digest"] == rep[b]["world_digest"]
    names = ("signal", "n0full", "n0strict", "n1", "n0full_s2", "n0strict_s2", "n1_s2")
    assert len({rep[n]["world_digest"] for n in names}) == len(names)
    for n, null in (("n0full", "N0-full"), ("n0strict", "N0-strict"), ("n1", "N1")):
        man = json.loads((W[n] / "manifest.json").read_text())
        assert man["world_id"] == f"micro-{null}-{SEED}" and man["null"] == {"world": null, "null_seed": SEED}
        lines = "".join(f"{k}\t{v}\n" for k, v in sorted(_digests(W[n]).items()))
        assert man["world_digest"] == rep[n]["world_digest"] == blake3(lines.encode()).hexdigest(length=16)
        assert {r["world_id"] for r in _truth(W[n], "truth_module")} <= {man["world_id"]}
    assert json.loads((W["n1_s2"] / "manifest.json").read_text())["world_id"] == f"micro-N1.{SEED2}-{SEED}"
    assert json.loads((W["n0full_s2"] / "manifest.json").read_text())["world_id"] == f"micro-N0-full-{SEED2}"
    print("\nworld digests: " + "; ".join(f"{n} {rep[n]['world_digest']}" for n in names))


# --- N0: the preset, what is kept exactly, what is kept in distribution ---------------------------------------------

@needs_zombi2
def test_n0_recorded_presets_differ_from_the_signal_only_by_the_removed_groups(worlds):
    W, _ = worlds
    sig = _preset(W["signal"])
    full = {("modules", "count"), ("profile_only", "count"), ("decoys", "partner_modules"), ("fusions", "nogd"),
            ("null",)}
    for n, want in (("n0full", full), ("n0strict", full | {("vehicles", "cargo"), ("crispr", "arrays")})):
        p = _preset(W[n])
        assert _leaf_diff(sig, p) == want, n
        assert p == _preset(W[n + "_s2"])         # another seed, the same null preset
        assert json.loads((W[n + "_s2"] / "backbone" / "backbone.json").read_text())["seed"] == SEED2


@needs_zombi2
@pytest.mark.covers("F-02", "B-NULL")
def test_n0_keeps_the_species_tree_lineages_phages_and_every_per_genome_draw_exactly(worlds):
    W, _ = worlds
    s = W["signal"]
    same_files = ("zombi2/species_complete.nwk", "backbone/taxonomy.tsv", "backbone/taxa.tsv", "sequences/lineages.tsv",
                  "replicons/lineages.tsv", "mobile/phages.fna", "mobile/phage_genes.tsv")
    for n in ("n0full", "n0strict"):
        w = W[n]
        for f in same_files:
            assert (w / f).read_bytes() == (s / f).read_bytes(), (n, f)
        assert _tree_sha(w / "mobile" / "phages") == _tree_sha(s / "mobile" / "phages"), n
        bs, bn = (json.loads((x / "backbone" / "backbone.json").read_text()) for x in (s, w))
        assert bn["burst_tips"] == bs["burst_tips"] and bn["tool_seeds"] == bs["tool_seeds"]
        assert [r["species_node"] for r in _tsv(w / "backbone" / "genomes.tsv")] == \
            [r["species_node"] for r in _tsv(s / "backbone" / "genomes.tsv")]
        # S6 prophage placement, S8 codes, S9 classes and targets, S10 roles and donors, S11 labels
        key = ("prophage_id", "phage_id", "clade_node", "genome_id", "strand")
        assert [tuple(r[k] for k in key) for r in _tsv(w / "mobile" / "prophages.tsv")] == \
            [tuple(r[k] for k in key) for r in _tsv(s / "mobile" / "prophages.tsv")]
        key = ("genome_id", "genetic_code", "assembly_class", "target_n50", "target_completeness")
        assert [tuple(r[k] for k in key) for r in _tsv(w / "assembly" / "genomes.tsv")] == \
            [tuple(r[k] for k in key) for r in _tsv(s / "assembly" / "genomes.tsv")]
        plan = lambda x: sorted((r["callerr_id"], r["kind"], r["genome_id"])          # noqa: E731
                                for r in _tsv(x / "assembly" / "callerr.tsv") if r["kind"] != "alt_code")
        assert plan(w) == plan(s)
        key = ("source_label", "genome_id", "kind", "rep_genome_id", "accession", "bioproject")
        assert [tuple(r[k] for k in key) for r in _tsv(w / "sources" / "sources.tsv")] == \
            [tuple(r[k] for k in key) for r in _tsv(s / "sources" / "sources.tsv")]
        contam = lambda x: collections.Counter((r["host_genome_id"], r["donor_genome_id"])   # noqa: E731
                                               for r in _tsv(x / "sources" / "contam.tsv"))
        assert contam(w) == contam(s)
        dup = lambda x: {(r["kind"], r["rep_genome_id"], r["member_genome_id"], r["member_source_label"])  # noqa: E731
                         for r in _tsv(x / "sources" / "dup.tsv")}     # one row per contig: contigs are content
        assert dup(w) == dup(s)
        drop = {"world_id", "release_tag", "genome_uid", "n_contigs", "length", "gc_true", "completeness_true"}
        tg = lambda x: [{k: v for k, v in r.items() if k not in drop} for r in _truth(x, "truth_genome")]  # noqa: E731
        assert tg(w) == tg(s)
        # §2.4: the background tight rule gives a pair of families the same class in both worlds
        cn, cs = _codir_classes(w), _codir_classes(s)
        common = set(cn) & set(cs)
        assert len(common) >= 100 and all(cn[k] == cs[k] for k in common)
        print(f"\n{n}: identical species tree, taxonomy, lineage tables, phage set, prophage clades, codes, classes, "
              f"call-error plan, sources, contamination and duplicate roles, labels; {len(common)} shared "
              f"background co-directional family pairs, all of one class in both worlds")


@needs_zombi2
@pytest.mark.covers("F-02", "B-NULL")
def test_n0_keeps_the_background_process_in_distribution(worlds, record_property):
    """Gene content is a new draw of the same process (ZOMBI2's gene-level draws change with the groups), so the
    background is compared in distribution, each with its stated test:

    - background genome size: per genome, the root and originated genes; geometric mean of the per-genome ratio
      null / signal within [0.9, 1.1] (a third of §5.1's +-30% size band) and Pearson r >= 0.9 over the 46 genomes
      (size is set by the per-lineage rates, which ZOMBI2 draws with the tree);
    - phyletic patterns: the genomes-present counts of the 1,000 root families, two-sample Kolmogorov-Smirnov
      D below the alpha = 0.01 critical value 1.628 sqrt(2 / 1000);
    - strands: the + share of genes in [0.4, 0.6] (backbone.SCRAMBLE_PLUS, as the signal world);
    - operons: the tight share of distinct background co-directional family pairs within 4 binomial standard errors
      of intergenic.junctions.p_tight (the hash rule is a fixed law);
    - composition: genome GC within 0.01 of the same genome's in the signal world (codon usage and intergenic
      composition are the lineage's, identical; gene content moves it only through the family mix).
    """
    W, _ = worlds
    s = W["signal"]

    def bg_sizes(w):
        kind = {int(r["family"]): r["kind"] for r in _tsv(w / "backbone" / "families.tsv")}
        c = collections.Counter(r["genome_id"] for r in _tsv(w / "backbone" / "genes.tsv")
                                if kind[int(r["family"])] in ("root", "originated"))
        return np.array([c[g] for g in sorted(c)], dtype=float)

    def ks(a, b):
        a, b = np.sort(a), np.sort(b)
        x = np.union1d(a, b)
        return float(np.max(np.abs(np.searchsorted(a, x, "right") / len(a) - np.searchsorted(b, x, "right") / len(b))))

    def gc(w):
        return np.array([float(r["gc_true"]) for r in _truth(w, "truth_genome") if r["is_phage"] == "false"])

    p_tight = float(_preset(s)["intergenic"]["junctions"]["p_tight"])
    root_s = [int(r["n_genomes_present"]) for r in _tsv(s / "backbone" / "families.tsv") if r["kind"] == "root"]
    sig_sizes, sig_gc = bg_sizes(s), gc(s)
    for n in ("n0full", "n0strict"):
        w = W[n]
        sz = bg_sizes(w)
        ratio = float(np.exp(np.mean(np.log(sz / sig_sizes))))
        r = float(np.corrcoef(sz, sig_sizes)[0, 1])
        assert 0.9 <= ratio <= 1.1 and r >= 0.9, (n, ratio, r)
        root = [int(x["n_genomes_present"]) for x in _tsv(w / "backbone" / "families.tsv") if x["kind"] == "root"]
        d = ks(np.array(root), np.array(root_s))
        assert d < 1.628 * math.sqrt((len(root) + len(root_s)) / (len(root) * len(root_s))), (n, d)
        genes = _tsv(w / "backbone" / "genes.tsv")
        plus = sum(g["strand"] == "1" for g in genes) / len(genes)
        assert backbone.SCRAMBLE_PLUS[0] <= plus <= backbone.SCRAMBLE_PLUS[1]
        kinds = _host_kind(w)
        cl = {k: v for k, v in _codir_classes(w).items() if kinds.get(k[0]) in ("root", "originated")
              and kinds.get(k[1]) in ("root", "originated")}
        share = sum(v == "codir_tight" for v in cl.values()) / len(cl)
        assert abs(share - p_tight) <= 4 * math.sqrt(p_tight * (1 - p_tight) / len(cl)), (n, share, len(cl))
        dgc = np.abs(gc(w) - sig_gc)
        assert dgc.max() <= 0.01, (n, dgc.max())
        record_property(f"{n}_bg_size_ratio", round(ratio, 4))
        print(f"\n{n}: background genes per genome {sz.mean():.1f} vs signal {sig_sizes.mean():.1f} (geometric mean "
              f"ratio {ratio:.3f}, r {r:.3f}); root presence KS D {d:.4f}; + strand {plus:.3f}; tight share {share:.4f} "
              f"over {len(cl)} distinct pairs (p_tight {p_tight}); GC |diff| mean {dgc.mean():.5f} max {dgc.max():.5f}")


@needs_zombi2
@pytest.mark.covers("F-02", "B-NULL")
def test_n0_has_no_planted_signal(worlds, record_property):
    """Truth: no planted or profile_only rows, no decoy family, no NOGD, no planted junction. Genomes: two brute-force
    measures of every adjacent family pair on the emitted contigs, the genomes in which it is adjacent and its distinct
    flank contexts. tau is the 99.9th percentile of the signal world's background pairs (not planted, not within one
    vehicle, plasmid or Cas locus). Positive control: planted pairs exceed tau far more often than background pairs.
    Then N0's pairs (the same exclusions) exceed tau no more often than the signal background, by a one-sided exact
    conditional binomial test at alpha = 0.001 (two measures x two worlds); leaked modules (the signal's 111 planted
    pairs, 34 to 61 above tau) would fail it."""
    W, _ = worlds
    s = W["signal"]
    gs, cs = _adjacency(s)
    planted = {k for k in gs if _within(k, [set(m) for m in _modules(s).values()])}
    linked = _linked_sets(s)
    bgk = [k for k in gs if k not in planted and not _within(k, linked)]
    taus = {}
    for name, M in (("genomes", gs), ("contexts", cs)):
        tau = float(np.percentile([len(M[k]) for k in bgk], 99.9))
        kb, kp = sum(len(M[k]) > tau for k in bgk), sum(len(M[k]) > tau for k in planted)
        assert _excess_test(kp, len(planted), kb, len(bgk)) < 1e-6, (name, kp, kb)
        taus[name] = (tau, kb, len(bgk), kp, len(planted))
    for n in ("n0full", "n0strict"):
        w = W[n]
        assert not _modules(w, "planted") and not _modules(w, "profile_only")
        assert not [r for r in _truth(w, "truth_module") if r["member_role"].startswith("nogd")]
        kinds = collections.Counter(_host_kind(w).values())
        assert not kinds["planted"] and not kinds["profile_only"] and not kinds["decoy"]
        assert not [r for r in _tsv(w / "fusions" / "events.tsv") if r["kind"] == "nogd"]
        assert not [j for j in _tsv(w / "replicons" / "junctions.tsv") if j["rule"] == "group"]
        assert not [t for t in _truth(w, "truth_traps") if t["kind"] == "decoy_profile"]
        g0, c0 = _adjacency(w)
        lk = _linked_sets(w)
        k0 = [k for k in g0 if not _within(k, lk)]
        for name, M in (("genomes", g0), ("contexts", c0)):
            tau, kb, nb, kp, npl = taus[name]
            ka = sum(len(M[k]) > tau for k in k0)
            p = _excess_test(ka, len(k0), kb, nb)
            assert p >= 1e-3, (n, name, ka, len(k0), kb, nb, p)
            record_property(f"{n}_{name}_above_tau", ka)
            print(f"\n{n} {name}: tau {tau:g}; signal planted {kp}/{npl}, signal background {kb}/{nb}, null {ka}/"
                  f"{len(k0)} above it (one-sided p {p:.3g})")


@needs_zombi2
@pytest.mark.covers("F-02", "B-NULL")
def test_n0_strict_has_no_vehicles_or_cas_loci_and_keeps_plasmids_and_prophages(worlds):
    W, _ = worlds
    p = _preset(W["signal"])
    for n in ("n0strict", "n0strict_s2"):
        w = W[n]
        kinds = collections.Counter(_host_kind(w).values())
        assert not kinds["vehicle_chassis"] and not kinds["cargo"]
        assert not _tsv(w / "mobile" / "cas.tsv") and not _tsv(w / "mobile" / "arrays.tsv")
        assert not _truth(w, "truth_spacer_target") and not _truth(w, "truth_module")
        assert not [t for t in _truth(w, "truth_traps") if t["kind"] == "mge_cargo"]
        # plasmids stay (§5.3); a plasmid type can die out at another seed (loss 0.7 > transfer 0.3 per unit time)
        reps = [r for r in _tsv(w / "backbone" / "replicons.tsv") if r["kind"] == "plasmid"]
        assert (reps or n.endswith("_s2")) and {int(r["n_genes"]) for r in reps} <= set(p["plasmids"]["sizes"])
        assert len({r["prophage_id"] for r in _tsv(w / "mobile" / "prophages.tsv")}) == p["prophages"]["count"]
    for n in ("n0full", "n0full_s2"):
        w = W[n]
        assert len({r["array_id"] for r in _tsv(w / "mobile" / "arrays.tsv")}) == p["crispr"]["arrays"]
        units = collections.defaultdict(int)
        for u in _tsv(w / "backbone" / "truth_units.tsv"):
            if u["group"].startswith("veh"):
                assert u["intact"] == "true"
                units[u["group"]] = int(u["n_members"])
        assert set(units.values()) <= {c + 2 for c in p["vehicles"]["cargo"]}


@needs_zombi2
@pytest.mark.covers("F-02", "B-NULL")
def test_another_seed_gives_another_null_world_with_the_same_structure(worlds):
    W, rep = worlds
    p = _preset(W["signal"])
    cnt = p["callerr"]["counts"]
    for n in ("n0full", "n0strict"):
        a, b = W[n], W[n + "_s2"]
        assert rep[n]["world_digest"] != rep[n + "_s2"]["world_digest"]
        assert (a / "zombi2" / "species_complete.nwk").read_bytes() != (b / "zombi2" / "species_complete.nwk").read_bytes()
        for w in (a, b):
            assert len(_tsv(w / "backbone" / "genomes.tsv")) == p["species"]["n_extant"]
            assert len(_tsv(w / "mobile" / "phages.tsv")) == p["phages"]["count"]
            ev = collections.Counter(r["kind"] for r in _tsv(w / "fusions" / "events.tsv"))
            assert (ev["fusion"], ev["split"], ev["nogd"]) == (p["fusions"]["fusions"], p["fusions"]["splits"], 0)
            plan = collections.Counter(r["kind"] for r in _tsv(w / "assembly" / "callerr.tsv") if r["kind"] != "alt_code")
            assert dict(plan) == dict(cnt)
            src = collections.Counter(r["kind"] for r in _tsv(w / "sources" / "sources.tsv"))
            assert src["exact_dup"] == p["duplicates"]["exact"] and src["near_dup"] == p["duplicates"]["near"]
            assert len({r["host_genome_id"] for r in _tsv(w / "sources" / "contam.tsv")}) == p["contamination"]["genomes"]
            assert not _modules(w, "planted") and not _modules(w, "profile_only")
            kinds = _host_kind(w)
            cl = {k: v for k, v in _codir_classes(w).items() if kinds.get(k[0]) in ("root", "originated")
                  and kinds.get(k[1]) in ("root", "originated")}
            share = sum(v == "codir_tight" for v in cl.values()) / len(cl)
            q = float(p["intergenic"]["junctions"]["p_tight"])
            assert abs(share - q) <= 4 * math.sqrt(q * (1 - q) / len(cl))
        strip = lambda w: [{k: v for k, v in r.items() if k not in ("world_id", "release_tag")}   # noqa: E731
                           for r in _truth(w, "truth_null")]
        assert strip(a) == strip(b)


@needs_zombi2
@pytest.mark.covers("F-02", "B-NULL")
def test_n0_truth_null_says_what_was_removed_and_how_pairs_are_evaluated(worlds):
    W, _ = worlds
    p = _preset(W["signal"])
    fid_of = lambda w: {int(r["family"]): r["family_id"] for r in _tsv(w / "backbone" / "families.tsv")}  # noqa: E731
    for n in ("n0full", "n0strict"):
        w = W[n]
        rows = _truth(w, "truth_null")
        assert not _truth(w, "truth_null_relabel")
        by = collections.defaultdict(list)
        for r in rows:
            assert r["null"] == {"n0full": "N0-full", "n0strict": "N0-strict"}[n]
            by[(r["kind"], r["action"])].append(r)
        assert [r["item_id"] for r in by[("planted_module", "removed")]] == \
            [f"mod{i:02d}" for i in range(p["modules"]["count"])]
        assert len(by[("profile_only", "removed")]) == p["profile_only"]["count"]
        assert len(by[("decoy", "removed")]) == len(p["decoys"]["partner_modules"])
        assert len(by[("nogd", "removed")]) == p["fusions"]["nogd"]
        fid = fid_of(w)
        groups = _tsv(w / "backbone" / "groups.tsv")
        veh = {g["name"]: [fid[int(x)] for x in g["members"].split(";")] for g in groups if g["kind"] == "vehicle"}
        pls = {g["name"]: [fid[int(x)] for x in g["members"].split(";")] for g in groups if g["kind"] == "plasmid"}
        cas = collections.defaultdict(list)
        for r in _tsv(w / "mobile" / "cas.tsv"):
            cas[r["cas_module_id"]].append(r["family_id"])
        assert {r["item_id"]: r["family_ids"].split(";") for r in by[("vehicle_chassis", "excluded")]} == \
            {v: [m[0], m[-1]] for v, m in veh.items()}
        assert {r["item_id"]: sorted(r["family_ids"].split(";")) for r in by[("cas_operon", "excluded")]} == \
            {k: sorted(v) for k, v in cas.items()}
        assert {r["item_id"]: r["family_ids"].split(";") for r in by[("plasmid", "null")]} == pls
        assert [r["family_ids"] for r in by[("all_pairs", "null")]] == ["*"]
        if n == "n0strict":
            assert len(by[("vehicle", "removed")]) == len(p["vehicles"]["cargo"])
            assert len(by[("crispr_cas_locus", "removed")]) == p["crispr"]["arrays"]
            assert not veh and not cas
        else:
            assert veh and len(cas) == p["crispr"]["arrays"]
            assert {r["kind"] for r in _truth(w, "truth_module")} == {"vehicle_chassis", "cas_operon"}


@needs_zombi2
@pytest.mark.covers("F-02", "B-NULL")
def test_n0_naive_cooccurrence_statistic_is_miscalibrated_section_4_3(worlds, record_property):
    """§4.3 positive control for confounding: on N0 a hypergeometric co-occurrence test over the S8 genomes, with no
    clade conditioning and no dedup, must show genomic inflation lambda >= 1.5 or >= 5x the expected share of pairs at
    p < 1e-3. The tail is summed here from log-binomials; families in 2 to G - 1 genomes. (The test is discrete and
    conservative, so lambda from the median p-value sits below 1 even under structure; the excess is the criterion
    that can show it.)"""
    W, _ = worlds
    for n in ("n0full", "n0strict"):
        pres = _presence(W[n])
        G = len({g for s in pres.values() for g in s})
        gi = {g: i for i, g in enumerate(sorted({g for s in pres.values() for g in s}))}
        fams = sorted(f for f, s in pres.items() if 2 <= len(s) <= G - 1)
        P = np.zeros((len(fams), G), dtype=np.int32)
        for i, f in enumerate(fams):
            P[i, [gi[g] for g in pres[f]]] = 1
        C, k = P @ P.T, P.sum(1)
        iu = np.triu_indices(len(fams), 1)
        p = _hypergeom_sf_table(G)[k[iu[0]], k[iu[1]], C[iu]]
        lam = NormalDist().inv_cdf(1 - float(np.median(p)) / 2) ** 2 / 0.454936
        excess = float(np.mean(p < 1e-3)) / 1e-3
        assert lam >= 1.5 or excess >= 5.0, (n, lam, excess)
        record_property(f"{n}_naive_excess", round(excess, 2))
        print(f"\n{n}: naive co-occurrence over {len(fams)} families, {len(p)} pairs: lambda {lam:.3f}, "
              f"share at p < 1e-3 {excess / 1000:.4f} ({excess:.1f}x the 0.001 a calibrated test gives)")


# --- N1 -------------------------------------------------------------------------------------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "MR-18")
def test_n1_shares_s1_to_s10_and_the_emitted_contigs_with_its_signal_world(worlds):
    W, _ = worlds
    s = W["signal"]
    for n in ("n1", "n1_s2"):
        w = W[n]
        for d in nulls.SHARED:
            assert _tree_sha(w / d) == _tree_sha(s / d), (n, d)
        assert {d: v for d, v in _digests(w).items() if d != "emit"} == {d: v for d, v in SIGNAL.items() if d != "emit"}
        assert _tree_sha(w / "emit" / "fasta") == _tree_sha(s / "emit" / "fasta")
        for f in ("ingest.tsv", "genomes.json", "truth.gff"):
            assert (w / "emit" / f).read_bytes() == (s / "emit" / f).read_bytes(), f
        strip = lambda x, t, drop=(): [{k: v for k, v in r.items() if k not in ("world_id",) + drop}  # noqa: E731
                                       for r in _truth(x, t)]
        for t in emit.TABLES:
            if t in ("truth_module", "truth_gene", "truth_fusion"):
                continue
            assert strip(w, t) == strip(s, t), (n, t)
        assert strip(w, "truth_gene", ("module_id",)) == strip(s, "truth_gene", ("module_id",))
        assert strip(w, "truth_fusion", ("in_module",)) == strip(s, "truth_fusion", ("in_module",))


@needs_zombi2
@pytest.mark.covers("F-02", "MR-18")
def test_n1_relabel_keeps_profiles_and_transfers_within_the_bounds(worlds, record_property):
    """Each replacement, counted here from the S8 genomes and the backbone events: a root or originated family, used
    once, within floor(0.1 x 46) = 4 genomes of the module family's profile and within 25% of its transfers; an
    unmatched member has no free family that qualifies (brute force in the same order). Kept co-occurrence: for every
    evaluated label pair the shared genomes change by at most the sum of the two Hamming distances (exact)."""
    W, _ = worlds
    s = W["signal"]
    pres, tr, kind = _presence(s), _transfers(s), _host_kind(s)
    cfg = _preset(s)["nulls"]["n1"]
    G = len(_tsv(s / "backbone" / "genomes.tsv"))
    max_h = math.floor(cfg["max_hamming_fraction"] * G + 1e-9)
    tol = cfg["transfer_tolerance"]
    pool = sorted(f for f, k in kind.items() if k in cfg["kinds"])
    sig_mods = _modules(s)
    for n in ("n1", "n1_s2"):
        w = W[n]
        rows = _truth(w, "truth_null_relabel")
        assert [(r["module_id"], r["family_id"]) for r in rows] == \
            [(m, f) for m, fs in sig_mods.items() for f in fs]
        used = set()
        for r in rows:
            a = r["family_id"]
            free = [x for x in pool if x not in used]
            ok = [x for x in free if len(pres.get(a, frozenset()) ^ pres.get(x, frozenset())) <= max_h
                  and abs(tr[x] - tr[a]) <= tol * tr[a] + 1e-9]
            if r["status"] == "relabeled":
                x = r["replacement_family_id"]
                assert x in ok and kind[x] in cfg["kinds"]
                assert int(r["hamming"]) == len(pres[a] ^ pres.get(x, frozenset())) <= max_h
                assert (int(r["transfers"]), int(r["replacement_transfers"])) == (tr[a], tr[x])
                used.add(x)
            else:
                assert r["status"] == "unmatched" and not ok
        mapping = {r["family_id"]: r["replacement_family_id"] for r in rows if r["status"] == "relabeled"}
        ham = {r["family_id"]: int(r["hamming"]) for r in rows if r["status"] == "relabeled"}
        mods = _modules(w)
        assert mods.keys() == sig_mods.keys()
        n_pairs = 0
        for m, fs in sig_mods.items():
            assert mods[m] == [mapping.get(f, f) for f in fs]
            for i, a in enumerate(fs):
                for b in fs[i + 1:]:
                    if a not in mapping and b not in mapping:
                        continue
                    a2, b2 = mapping.get(a, a), mapping.get(b, b)
                    before = len(pres[a] & pres[b])
                    after = len(pres.get(a2, frozenset()) & pres.get(b2, frozenset()))
                    assert abs(after - before) <= ham.get(a, 0) + ham.get(b, 0)
                    n_pairs += 1
        assert not [r for r in _truth(w, "truth_module") if r["member_role"].startswith("nogd")]
        # module_id follows the labels: genes of a replacement carry its module, genes of a replaced family none
        mod_of = {f: m for m, fs in mods.items() for f in fs}
        for g in _truth(w, "truth_gene"):
            fams = g["family_id"].split("+")
            if any(f in mod_of and f in mapping.values() for f in fams):
                assert g["module_id"] == next(mod_of[f] for f in fams if f in mod_of and f in mapping.values())
            elif any(f in mapping for f in fams):
                assert g["module_id"] == ""
        record_property(f"{n}_relabeled", len(mapping))
        record_property(f"{n}_evaluated_pairs", n_pairs)
        print(f"\n{n}: {len(mapping)} of {len(rows)} module members relabeled ({sum(1 for m, fs in sig_mods.items() if any(f in mapping for f in fs))} "
              f"of {len(sig_mods)} modules), {n_pairs} evaluated label pairs; unmatched module families' transfers "
              f"{sorted(int(r['transfers']) for r in rows if r['status'] == 'unmatched')}")
    a = {r["family_id"]: r["replacement_family_id"] for r in _truth(W["n1"], "truth_null_relabel")}
    b = {r["family_id"]: r["replacement_family_id"] for r in _truth(W["n1_s2"], "truth_null_relabel")}
    # another null seed draws other replacements (and so frees or takes other candidates for later members)
    assert a != b and any(a.values()) and any(b.values())


@needs_zombi2
@pytest.mark.covers("F-02", "MR-18")
def test_n1_label_pairs_are_no_more_adjacent_than_matched_random_family_pairs(worlds, record_property):
    """Destroyed: the adjacency of module labels. The mean number of genomes in which a label pair is adjacent on the
    emitted contigs is compared with 2,000 draws of matched random pairs (each label replaced by a root or originated
    family present in the same number of genomes +-2, fixed seed); one-sided permutation p must be >= 0.01. Positive
    control: the signal world's module pairs under the same procedure have p <= 1/2001."""
    W, _ = worlds
    s = W["signal"]
    gs, _ = _adjacency(s)
    pres, kind = _presence(s), _host_kind(s)
    pool = sorted(f for f, k in kind.items() if k in ("root", "originated") and f in pres)
    by_n = collections.defaultdict(list)
    for f in pool:
        by_n[len(pres[f])].append(f)
    rng = np.random.default_rng(20261113)

    def adj(a, b):
        return len(gs.get((a, b) if a <= b else (b, a), ()))

    def perm_p(pairs):
        obs = float(np.mean([adj(a, b) for a, b in pairs]))
        cands = {}
        for f in {x for p in pairs for x in p}:
            k = len(pres.get(f, ()))
            cands[f] = [x for d in range(-2, 3) for x in by_n.get(k + d, ())]
        null = np.empty(2000)
        for i in range(2000):
            null[i] = np.mean([adj(cands[a][rng.integers(len(cands[a]))], cands[b][rng.integers(len(cands[b]))])
                               for a, b in pairs])
        return obs, float(null.mean()), (1 + int(np.sum(null >= obs))) / 2001

    sig_pairs = [(a, b) for fs in _modules(s).values() for i, a in enumerate(fs) for b in fs[i + 1:]]
    obs_s, null_s, p_s = perm_p(sig_pairs)
    assert p_s <= 1 / 2001 + 1e-12, (obs_s, null_s)
    for n in ("n1", "n1_s2"):
        w = W[n]
        moved = {r["replacement_family_id"] for r in _truth(w, "truth_null_relabel") if r["status"] == "relabeled"}
        pairs = [(a, b) for fs in _modules(w).values() for i, a in enumerate(fs) for b in fs[i + 1:]
                 if a in moved or b in moved]
        obs, null_mean, p = perm_p(pairs)
        assert p >= 0.01, (n, obs, null_mean, p)
        record_property(f"{n}_adjacency_p", round(p, 4))
        print(f"\n{n}: {len(pairs)} evaluated label pairs adjacent in {obs:.3f} genomes on average, matched random "
              f"pairs {null_mean:.3f}, one-sided p {p:.3f}; signal module pairs {obs_s:.2f} against {null_s:.3f} "
              f"(p {p_s:.4f})")


# --- time and memory ---------------------------------------------------------------------------------------------------

@needs_zombi2
def test_null_world_times_and_peak_memory(worlds, record_property):
    _, rep = worlds
    for n, r in rep.items():
        if n == "_wall":
            continue
        record_property(f"{n}_seconds", r["seconds"])
        if r["vmhwm_kb"] > 0:
            record_property(f"{n}_vmhwm_mb", round(r["vmhwm_kb"] / 1024, 1))
            assert r["vmhwm_kb"] / 1024 < 1024, n
    record_property("fixture_wall_seconds", rep["_wall"])
    print("\n" + "; ".join(f"{n} {r['seconds']} s, VmHWM {r['vmhwm_kb'] / 1024:.0f} MB" for n, r in rep.items()
                           if n != "_wall") + f"; fixture wall {rep['_wall']} s with {JOBS} at once")
