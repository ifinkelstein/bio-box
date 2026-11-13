"""Ring R2: synthworld family design and protein-coding sequences on the `micro` preset (synthworld_spec
§2.1 S3-S4, §2.3, §5.1, §5.3, §6.1; card F-02 session 3).

A module fixture builds the micro backbone (ZOMBI2 0.54.0) and its sequences once. The oracles are written
here: the genetic code is spelled out codon by codon, FASTA and the gene trees are parsed by this file,
identities are counted residue by residue, and the profile-F81 expectation is computed from its closed form.
What the oracles take from the generator are model *inputs* (a family's rate, the lineage clock, a
family's site profile and rate categories), never its realized sequences or statistics.

The micro preset's default profile source is `synthetic` (no Pfam file in CI); the Pfam-seeded world and
the HMMER parity checks are in tests/r2/test_synthworld_sequences_pfam.py.
"""

from __future__ import annotations

import collections
import csv
import hashlib
import importlib.util
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest

from biobox import ids
from biobox.synthworld import backbone, sequences

pytestmark = [pytest.mark.r2, pytest.mark.covers("F-02")]
needs_zombi2 = pytest.mark.skipif(
    importlib.util.find_spec("zombi2") is None,
    reason="zombi2 0.54.0 is not importable here: CI installs it by `uv sync`; on the pod add its own "
           "directory with BIOBOX_POD_PYDEPS (card F-02)")

AA = "ACDEFGHIKLMNPQRSTVWY"
CODE11_BY_AA = {
    "F": "TTT TTC", "L": "TTA TTG CTT CTC CTA CTG", "I": "ATT ATC ATA", "M": "ATG", "V": "GTT GTC GTA GTG",
    "S": "TCT TCC TCA TCG AGT AGC", "P": "CCT CCC CCA CCG", "T": "ACT ACC ACA ACG", "A": "GCT GCC GCA GCG",
    "Y": "TAT TAC", "*": "TAA TAG TGA", "H": "CAT CAC", "Q": "CAA CAG", "N": "AAT AAC", "K": "AAA AAG",
    "D": "GAT GAC", "E": "GAA GAG", "C": "TGT TGC", "W": "TGG", "R": "CGT CGC CGA CGG AGA AGG",
    "G": "GGT GGC GGA GGG",
}
CODE11 = {c: aa for aa, cs in CODE11_BY_AA.items() for c in cs.split()}
STARTS11 = {"TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"}   # NCBI transl_table=11 initiation codons


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


def _copy_backbone(src: Path, dst: Path) -> Path:
    dst.mkdir(parents=True)
    for d in ("backbone", "zombi2"):
        shutil.copytree(src / d, dst / d)
    return dst


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    if importlib.util.find_spec("zombi2") is None:
        pytest.skip("zombi2 0.54.0 is not importable here")
    out = tmp_path_factory.mktemp("micro_seq") / "world"
    params = backbone.load_preset("micro")
    backbone.build(params, params["default_seed"], out)
    t0 = time.perf_counter()
    sw = sequences.build(out)
    return sw, out, time.perf_counter() - t0


# --- determinism (§6.1, §6.2; MR-16) -----------------------------------------------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "MR-16")
def test_same_seed_gives_byte_identical_sequences_in_a_subprocess(world, tmp_path, record_property):
    sw, out, _ = world
    dst = _copy_backbone(out, tmp_path / "copy")
    env = dict(os.environ, PYTHONHASHSEED="4242", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    code = ("import resource, sys; from biobox.synthworld import sequences; sequences.build(sys.argv[1]); "
            "print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)")
    proc = subprocess.run([sys.executable, "-c", code, str(dst)], check=True, env=env, cwd=tmp_path,
                          capture_output=True, text=True)
    assert _tree_sha(dst / "sequences") == _tree_sha(out / "sequences")
    # bounded memory: the sequence stage holds one family at a time (the child's own peak RSS, Linux KiB)
    rss_mb = int(proc.stdout.split()[-1]) / 1024
    record_property("sequences_child_maxrss_mb", round(rss_mb, 1))
    print(f"\nsequence stage child peak RSS {rss_mb:.0f} MB")
    assert rss_mb < 1024


@needs_zombi2
def test_another_seed_gives_other_sequences(world, tmp_path):
    sw, out, _ = world
    params = backbone.load_preset("micro")
    other = tmp_path / "other"
    backbone.build(params, params["default_seed"] + 1, other)
    sw2 = sequences.build(other)
    assert sw2.digest != sw.digest
    a, b = _fasta(out / "sequences/roots.faa"), _fasta(other / "sequences/roots.faa")
    common = sorted(set(a) & set(b))
    # family ids common to both worlds are mostly root families with an extant copy in both
    # (391 at seeds 20261112 and 20261113, measured on wilkcomp02 2026-11-12)
    assert len(common) >= 200
    assert sum(a[k] == b[k] for k in common) == 0


# --- coding sequences (§2.3 "Codons"; §3.2 truth_gene; §7.3 O-TRANSLATE) ---------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "O-TRANSLATE")
def test_every_cds_translates_to_its_protein_with_code_11(world):
    sw, out, _ = world
    cds = _fasta(out / "sequences/cds.fna")
    prot = _fasta(out / "sequences/proteins.faa")
    genes = {r["copy_label"]: r for r in _tsv(out / "sequences/genes.tsv")}
    raw = {r["copy_label"] for r in _tsv(out / "backbone/genes.tsv")}
    assert set(cds) == set(prot) == set(genes) == raw and len(raw) > 40_000
    for name, s in cds.items():
        assert len(s) % 3 == 0 and set(s) <= set("ACGT"), name
        codons_ = [s[i:i + 3] for i in range(0, len(s), 3)]
        assert codons_[0] in STARTS11, name
        assert CODE11[codons_[-1]] == "*", name
        body = [CODE11[c] for c in codons_[1:-1]]
        assert "*" not in body, name
        assert prot[name] == "M" + "".join(body), name
        g = genes[name]
        assert int(g["cds_length"]) == len(s) and int(g["protein_length"]) == len(prot[name])
        assert g["start_codon"] == codons_[0] and g["stop_codon"] == codons_[-1]
        assert g["protein_uid"] == ids.protein_uid(prot[name]).hex()
        assert int(g["genetic_code"]) == 11


@needs_zombi2
def test_fasta_index_points_at_each_record(world):
    _, out, _ = world
    for fa in ("cds.fna", "proteins.faa"):
        path = out / "sequences" / fa
        data = path.read_bytes()
        recs = _fasta(path)
        rows = [line.split("\t") for line in (out / "sequences" / (fa + ".fai")).read_text().splitlines()]
        assert [r[0] for r in rows] == list(recs)
        for name, length, offset, linebases, linewidth in rows[::97]:
            n, off, lb, lw = int(length), int(offset), int(linebases), int(linewidth)
            got = b"".join(data[off + i * lw: off + i * lw + min(lb, n - i * lb)] for i in range(-(-n // lb)))
            assert got.decode() == recs[name]


@needs_zombi2
def test_genes_are_placed_on_the_realized_order_and_strand(world):
    _, out, _ = world
    cds = _fasta(out / "sequences/cds.fna")
    raw = _tsv(out / "backbone/genes.tsv")
    for gid in ("g000", "g017"):
        want = [(int(r["replicon"]), int(r["position"]), int(r["strand"]), r["copy_label"])
                for r in raw if r["genome_id"] == gid]
        got = list(sequences.iter_genome(out, gid))
        assert [(g.replicon, g.position, g.strand, g.copy_label) for g in got] == want
        assert all(g.cds == cds[g.copy_label] for g in got)


# --- divergence along the gene tree (§2.3 "Substitution", "Divergence control") ---------------------------

def _gene_trees(out: Path) -> dict[int, list[tuple[int, str, float]]]:
    """family -> [(parent index, label, branch length)] per node in text order, parsed here."""
    trees = {}
    for r in _tsv(out / "zombi2/gene_trees.tsv"):
        nodes, stack, last, closed = [], [], -1, False
        for tok in re.findall(r"\(|\)|,|;|[^(),;]+", r["complete"].strip()):
            if tok == "(":
                nodes.append([stack[-1] if stack else -1, "", 0.0])
                stack.append(len(nodes) - 1)
                closed = False
            elif tok == ")":
                last, closed = stack.pop(), True
            elif tok in (",", ";"):
                closed = False
            else:
                label, _, length = tok.partition(":")
                if closed:
                    nodes[last][1:] = [label, float(length or 0)]
                else:
                    nodes.append([stack[-1] if stack else -1, label, float(length or 0)])
                closed = False
        trees[int(r["family"])] = [tuple(n) for n in nodes]
    return trees


def _lineage(label: str) -> str:
    """'<lineage>_g<copy>' for a leaf, '<event>_<lineage>' for an internal node."""
    a, _, b = label.partition("_")
    return a if b.startswith("g") else b


def _path_length(nodes, leaf_a: int, leaf_b: int, clock: dict[str, float]) -> float:
    def up(x):
        chain = []
        while x != -1:
            chain.append(x)
            x = nodes[x][0]
        return chain
    ca, cb = up(leaf_a), up(leaf_b)
    common = set(ca) & set(cb)
    return sum(clock[_lineage(nodes[x][1])] * nodes[x][2] for x in ca + cb if x not in common)


def _spearman(x, y) -> float:
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    return float(np.corrcoef(rx, ry)[0, 1])


def _domain_string(protein: str, arch: str) -> str:
    out = []
    for part in arch.split(";"):
        _acc, span = part.rsplit(":", 1)
        s, e = (int(v) for v in span.split("-"))
        out.append(protein[s:e])
    return "".join(out)


@needs_zombi2
@pytest.mark.covers("F-02", "B-SIM")
def test_pairwise_identity_falls_with_tree_distance_as_the_model_expects(world):
    """Two copies of a family at path length T (gene-tree branch lengths times the lineage clock, from
    ZOMBI2's raw trees) share a residue at a domain site with probability e^{-mu r T} + (1 - e^{-mu r T}) h,
    h the site's sum of squared frequencies (§2.3). Up to 3 pairs per family; pairs binned into quintiles
    of the expectation; tolerance 0.02 per bin on the mean (derived: a pair's identity has sd at most about
    0.05 over its domain sites, and a bin holds at least 200 pairs, so 0.02 is more than 5 standard errors).
    Identity falls with distance: Spearman of mu T against identity at most -0.7."""
    sw, out, _ = world
    prot = _fasta(out / "sequences/proteins.faa")
    genes = _tsv(out / "sequences/genes.tsv")
    fams = {int(r["family"]): r for r in _tsv(out / "sequences/families.tsv")}
    clock = {r["lineage"]: float(r["clock"]) for r in _tsv(out / "sequences/lineages.tsv")}
    trees = _gene_trees(out)
    model = sequences.model(out)
    by_fam = collections.defaultdict(list)
    for g in genes:
        by_fam[int(g["family"])].append(g)
    rng = np.random.default_rng(7)
    rows = []
    for fam in sorted(by_fam):
        copies = by_fam[fam]
        if len(copies) < 2:
            continue
        design = model.design(fam)
        dom = np.isin(design.kind, (sequences.SITE_MATCH, sequences.SITE_INSERT))
        pi = design.pi_table[design.pi_idx[dom]]
        h, r = (pi ** 2).sum(axis=1), design.rates[dom]
        mu = float(fams[fam]["rate_multiplier"])
        nodes = trees[fam]
        leaf_index = {lab: i for i, (_, lab, _) in enumerate(nodes)}
        for _ in range(min(3, len(copies) * (len(copies) - 1) // 2)):
            i, j = rng.choice(len(copies), size=2, replace=False)
            a, b = copies[i], copies[j]
            sa = _domain_string(prot[a["copy_label"]], a["domain_arch"])
            sb = _domain_string(prot[b["copy_label"]], b["domain_arch"])
            assert len(sa) == len(sb) == int(dom.sum()), fam
            T = _path_length(nodes, leaf_index[a["copy_label"]], leaf_index[b["copy_label"]], clock)
            e = np.exp(-mu * r * T)
            rows.append((mu * T, float((e + (1 - e) * h).mean()),
                         sum(x == y for x, y in zip(sa, sb)) / len(sa)))
    rows = np.array(rows)
    assert len(rows) >= 1500, len(rows)
    order = np.argsort(rows[:, 1])
    report = []
    for chunk in np.array_split(order, 5):
        exp_mean, got_mean = rows[chunk, 1].mean(), rows[chunk, 2].mean()
        report.append((round(exp_mean, 3), round(got_mean, 3), len(chunk)))
        assert abs(exp_mean - got_mean) <= 0.02, report
    rho = _spearman(rows[:, 0], rows[:, 2])
    print(f"\nidentity by expectation quintile (expected, realized, pairs): {report}; spearman(mu T, id) "
          f"{rho:.3f} over {len(rows)} pairs")
    assert rho <= -0.7


@needs_zombi2
@pytest.mark.covers("F-02", "B-SIM")
def test_copies_of_unrelated_families_are_no_more_similar_than_chance(world):
    """Pairs of proteins from families that share no profile: identity at aligned positions (from position
    1, so the initiator M is not counted) exceeds the chance identity of the two compositions by at most
    0.01 on average over 2,000 pairs (derived: chance identity is about 0.06 and the per-pair sd about
    0.015, so the mean's standard error is below 0.001), and no pair shares more than 5% of its 5-mers
    (two unrelated 300-residue proteins share about 300^2/20^5, under 0.03, 5-mers in all)."""
    _, out, _ = world
    prot = _fasta(out / "sequences/proteins.faa")
    fams = {int(r["family"]): set(r["domain_hmms"].split(";")) for r in _tsv(out / "sequences/families.tsv")}
    one_copy = {}
    for g in _tsv(out / "sequences/genes.tsv"):
        one_copy.setdefault(int(g["family"]), g["copy_label"])
    keys = sorted(one_copy)
    rng = np.random.default_rng(13)
    excess, shared = [], []
    while len(excess) < 2000:
        fa, fb = (keys[i] for i in rng.choice(len(keys), size=2, replace=False))
        if fams[fa] & fams[fb]:
            continue
        a, b = prot[one_copy[fa]][1:], prot[one_copy[fb]][1:]
        n = min(len(a), len(b))
        ident = sum(x == y for x, y in zip(a[:n], b[:n])) / n
        ca = np.array([a.count(c) for c in AA], float) / len(a)
        cb = np.array([b.count(c) for c in AA], float) / len(b)
        excess.append(ident - float(ca @ cb))
        ka = {a[i:i + 5] for i in range(len(a) - 4)}
        kb = {b[i:i + 5] for i in range(len(b) - 4)}
        shared.append(len(ka & kb) / max(1, len(ka)))
    print(f"\nunrelated pairs: mean identity excess {np.mean(excess):.4f}, max shared 5-mer fraction "
          f"{max(shared):.4f}")
    assert abs(np.mean(excess)) <= 0.01
    assert max(shared) <= 0.05


# --- family design against the preset (§2.1 S3, §5.1, §5.3) -----------------------------------------------

def _within(count: int, n: int, p: float, k: float = 4.0) -> bool:
    return abs(count - n * p) <= k * math.sqrt(n * p * (1 - p)) + 1


@needs_zombi2
@pytest.mark.covers("F-02", "B-SIM")
def test_family_architecture_lengths_and_composition_follow_the_preset(world):
    """Domain counts, dark and promiscuous fractions within 4 binomial sd of the preset; tail and linker
    lengths (geometric, mean m, variance m(m+1)) within 4 sd of the mean; the synthetic panel's domain
    lengths log-normal around the preset median (within 15%) and spread (within 20%); residues outside
    domains drawn from the background (Jensen-Shannon divergence at most 0.002 nats); extant protein length
    median in the asserted band 200-400 aa (§7.1 fits it to the golden genomes later)."""
    _, out, _ = world
    cfg = backbone.load_preset("micro")["sequences"]
    arch, prof = cfg["architecture"], cfg["profiles"]
    fams = _tsv(out / "sequences/families.tsv")
    n = len(fams)
    w = arch["domain_count_weights"]
    for k, p in enumerate(w, start=1):
        assert _within(sum(int(f["n_domains"]) == k for f in fams), n, p), (k, p)
    dark = [f for f in fams if f["dark"] == "true"]
    assert _within(len(dark), n, arch["dark_fraction"])
    lit = [f for f in fams if f["dark"] == "false"]
    assert not any(f["promiscuous"] == "true" for f in dark)
    assert _within(sum(f["promiscuous"] == "true" for f in lit), len(lit), arch["promiscuous_fraction"])
    for col, key in (("n_tail", "n_tail_mean"), ("c_tail", "c_tail_mean")):
        m = float(arch[key])
        xs = [int(f[col]) for f in fams]
        assert abs(np.mean(xs) - m) <= 4 * math.sqrt(m * (m + 1) / len(xs)), (col, np.mean(xs))
    links = [int(x) for f in fams for x in f["linkers"].split(";") if x]
    m = float(arch["linker_mean"])
    assert abs(np.mean(links) - m) <= 4 * math.sqrt(m * (m + 1) / len(links))
    # root length is the sum of its parts
    doms = collections.defaultdict(list)
    for d in _tsv(out / "sequences/domains.tsv"):
        doms[int(d["family"])].append(int(d["aa_end"]) - int(d["aa_start"]))
    for f in fams:
        parts = 1 + int(f["n_tail"]) + int(f["c_tail"]) + sum(int(x) for x in f["linkers"].split(";") if x)
        assert int(f["root_length"]) == parts + sum(doms[int(f["family"])])
    # the synthetic panel's match lengths
    panel = _tsv(out / "sequences/panel.tsv")
    M = np.log([int(p["match_length"]) for p in panel])
    sm = prof["synthetic_match"]
    assert abs(np.exp(np.median(M)) / sm["median"] - 1) <= 0.15
    assert abs(M.std() / sm["sigma"] - 1) <= 0.2
    # composition outside domains is the background; protein lengths
    prot = _fasta(out / "sequences/proteins.faa")
    counts = collections.Counter()
    lengths = []
    for g in _tsv(out / "sequences/genes.tsv"):
        p = prot[g["copy_label"]]
        lengths.append(len(p))
        mask = np.ones(len(p), bool)
        mask[0] = False
        for part in g["domain_arch"].split(";"):
            s, e = (int(v) for v in part.rsplit(":", 1)[1].split("-"))
            mask[s:e] = False
        counts.update(c for c, keep in zip(p, mask) if keep)
    f = np.array([counts[c] for c in AA], float)
    f /= f.sum()
    from biobox.synthworld import profiles
    bg = profiles.BACKGROUND
    mid = (f + bg) / 2
    js = 0.5 * float((f * np.log(f / mid)).sum() + (bg * np.log(bg / mid)).sum())
    med = float(np.median(lengths))
    print(f"\nnon-domain composition JS {js:.5f} nats over {sum(counts.values())} residues; protein length "
          f"median {med:.0f}, 5% {np.quantile(lengths, 0.05):.0f}, 95% {np.quantile(lengths, 0.95):.0f}")
    assert js <= 0.002
    assert 200 <= med <= 400


@needs_zombi2
@pytest.mark.covers("F-02", "B-SIM")
def test_gc3_follows_the_lineage_codon_usage(world):
    """§2.3 "GC and codon usage follow a Brownian process on the species tree": a genome's mean GC3 over
    its genes, counted here from the CDS, tracks its species branch's GC3 (Pearson r at least 0.6, slope
    positive), and the branches spread by at least 0.1."""
    _, out, _ = world
    cds = _fasta(out / "sequences/cds.fna")
    target = {r["lineage"]: float(r["gc3"]) for r in _tsv(out / "sequences/lineages.tsv")}
    node_of = {r["genome_id"]: r["species_node"] for r in _tsv(out / "backbone/genomes.tsv")}
    gc3 = collections.defaultdict(list)
    for r in _tsv(out / "backbone/genes.tsv"):
        s = cds[r["copy_label"]]
        thirds = s[5:-3:3]
        gc3[r["genome_id"]].append(sum(c in "GC" for c in thirds) / len(thirds))
    gids = sorted(gc3)
    x = np.array([target[node_of[g]] for g in gids])
    y = np.array([np.mean(gc3[g]) for g in gids])
    r = float(np.corrcoef(x, y)[0, 1])
    slope = float(np.polyfit(x, y, 1)[0])
    print(f"\ngenome GC3 vs lineage GC3: r {r:.3f}, slope {slope:.3f}, lineage range {x.min():.3f}-{x.max():.3f}, "
          f"genome range {y.min():.3f}-{y.max():.3f}")
    assert x.max() - x.min() >= 0.1
    assert r >= 0.6 and slope > 0


# --- runtime (§2.6) ---------------------------------------------------------------------------------------

@needs_zombi2
def test_micro_sequences_runtime_within_budget(world, record_property):
    sw, _, seconds = world
    record_property("sequences_seconds", round(seconds, 2))
    print(f"\nF-02 micro sequences: {seconds:.2f} s wall; {sw.meta['n_genes']} genes in {len(sw.families)} families; "
          f"digest {sw.digest}")
    assert seconds < 120.0  # synthworld_spec §2.6: the S4 evolver under 2 min single core
