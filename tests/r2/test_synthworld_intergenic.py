"""Ring R2: synthworld IGRs, regulatory elements, RNA motifs and replicon sequences on the `micro` preset
(synthworld_spec §2.1 S5, §2.4, §2.5, §2.6, §5.1, §6.1; card F-02 session 4).

A module fixture builds the micro backbone (ZOMBI2 0.54.0), its sequences and its replicons once. The oracles
are written here: FASTA is parsed, reverse complements, translation, base coverage, junction types from
strands, IGR GC, PWM matches, terminator pairing and the motif scan are computed by this file; the IGR
runs come from the store's own code (`biobox.l1.igr.igrs_of_contig`, `biobox.l1.ncrna.regions`), which
knows nothing of synthworld. What the oracles take from the generator are model inputs (the preset, a
clade's emission, a lineage's intergenic GC), never its realized statistics.

The micro preset's motif source is `synthetic` (CI has no Infernal or Rfam); the Rfam CMs and the cmsearch
parity are in tests/r2/test_synthworld_rna_rfam.py.
"""

from __future__ import annotations

import collections
import csv
import hashlib
import importlib.util
import json
import math
import os
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest
from blake3 import blake3

from biobox.l1 import igr as l1_igr
from biobox.l1 import ncrna as l1_ncrna
from biobox.synthworld import backbone, intergenic, sequences

pytestmark = [pytest.mark.r2, pytest.mark.covers("F-02")]
needs_zombi2 = pytest.mark.skipif(
    importlib.util.find_spec("zombi2") is None,
    reason="zombi2 0.54.0 is not importable here: CI installs it by `uv sync`; on the pod add its own "
           "directory with BIOBOX_POD_PYDEPS (card F-02)")

CODE11_BY_AA = {
    "F": "TTT TTC", "L": "TTA TTG CTT CTC CTA CTG", "I": "ATT ATC ATA", "M": "ATG", "V": "GTT GTC GTA GTG",
    "S": "TCT TCC TCA TCG AGT AGC", "P": "CCT CCC CCA CCG", "T": "ACT ACC ACA ACG", "A": "GCT GCC GCA GCG",
    "Y": "TAT TAC", "*": "TAA TAG TGA", "H": "CAT CAC", "Q": "CAA CAG", "N": "AAT AAC", "K": "AAA AAG",
    "D": "GAT GAC", "E": "GAA GAG", "C": "TGT TGC", "W": "TGG", "R": "CGT CGC CGA CGG AGA AGG",
    "G": "GGT GGC GGA GGG",
}
CODE11 = {c: aa for aa, cs in CODE11_BY_AA.items() for c in cs.split()}
STARTS = {"ATG", "GTG", "TTG"}
STOPS = {"TAA", "TAG", "TGA"}
PAIRS_OK = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G"), ("G", "T"), ("T", "G")}
_RC = {"A": "T", "C": "G", "G": "C", "T": "A"}


def rc(s: str) -> str:
    return "".join(_RC[c] for c in reversed(s))


def translate(cds: str) -> str:
    aa = [CODE11[cds[i:i + 3]] for i in range(0, len(cds), 3)]
    assert aa[-1] == "*" and "*" not in aa[:-1]
    return "M" + "".join(aa[1:-1])


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


def _cut(seq: str, start: int, end: int, strand: int) -> str:
    """Bases [start, end) of a circular or linear record (end > length across the origin), 5'->3' on strand."""
    n = len(seq)
    assert 0 <= start <= end <= start + n
    s = (seq + seq)[start:end] if end > n else seq[start:end]
    return s if strand == 1 else rc(s)


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    if importlib.util.find_spec("zombi2") is None:
        pytest.skip("zombi2 0.54.0 is not importable here")
    out = tmp_path_factory.mktemp("micro_igr") / "world"
    params = backbone.load_preset("micro")
    backbone.build(params, params["default_seed"], out)
    sequences.build(out)
    t0 = time.perf_counter()
    rw = intergenic.build(out)
    return rw, out, time.perf_counter() - t0


@pytest.fixture(scope="module")
def tables(world):
    return load_tables(world)


def load_tables(world) -> dict:
    rw, out, _ = world
    rd = out / "replicons"
    seqs = {}
    for f in sorted((rd / "fasta").glob("*.fna")):
        seqs.update(_fasta(f))
    reps = {(r["genome_id"], int(r["replicon"])): r for r in _tsv(rd / "replicons.tsv")}
    genes = _tsv(rd / "genes.tsv")
    juncs = _tsv(rd / "junctions.tsv")
    elems = _tsv(rd / "elements.tsv")
    rnas = _tsv(rd / "rna.tsv")
    return {"seqs": seqs, "reps": reps, "genes": genes, "juncs": juncs, "elems": elems, "rna": rnas,
            "seq_of": {k: seqs[r["record"]] for k, r in reps.items()}}


# --- determinism (§6.1, §6.2; MR-16) -----------------------------------------------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "MR-16")
def test_same_seed_gives_byte_identical_replicons_in_a_subprocess(world, tmp_path, record_property):
    rw, out, _ = world
    dst = tmp_path / "copy"
    dst.mkdir()
    for d in ("backbone", "zombi2", "sequences"):
        shutil.copytree(out / d, dst / d)
    env = dict(os.environ, PYTHONHASHSEED="4242", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    code = ("import os, sys; from biobox.synthworld import intergenic; intergenic.build(sys.argv[1]); "
            "p = '/proc/self/status'; "
            "print(next((l.split()[1] for l in open(p) if l.startswith('VmHWM')), -1) if os.path.exists(p) else -1)")
    proc = subprocess.run([sys.executable, "-c", code, str(dst)], check=True, env=env, cwd=tmp_path,
                          capture_output=True, text=True)
    assert _tree_sha(dst / "replicons") == _tree_sha(out / "replicons")
    hwm_kb = int(proc.stdout.split()[-1])
    if hwm_kb < 0:
        return
    rss_mb = hwm_kb / 1024
    record_property("replicons_child_maxrss_mb", round(rss_mb, 1))
    print(f"\nreplicon stage child peak RSS {rss_mb:.0f} MB")
    assert rss_mb < 1024      # one family or one genome at a time (AGENTS §4.1)


@needs_zombi2
def test_another_seed_gives_other_replicons(world, tmp_path):
    rw, out, _ = world
    params = backbone.load_preset("micro")
    other = tmp_path / "other"
    backbone.build(params, params["default_seed"] + 1, other)
    sequences.build(other)
    rw2 = intergenic.build(other)
    assert rw2.digest != rw.digest
    assert not {r.contig_uid for r in rw.replicons} & {r.contig_uid for r in rw2.replicons}
    assert [r.structure for r in rw.models] != [r.structure for r in rw2.models]


# --- genes on the replicons (§2.1 S5; AGENTS §3.1; §7.3 O-SLICE, O-TRANSLATE) -------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "O-SLICE", "O-TRANSLATE")
def test_every_gene_slices_to_its_coding_sequence(world, tables):
    _, out, _ = world
    cds = _fasta(out / "sequences/cds.fna")
    prot = _fasta(out / "sequences/proteins.faa")
    bb = {r["copy_label"]: r for r in _tsv(out / "backbone/genes.tsv")}
    assert {g["copy_label"] for g in tables["genes"]} == set(bb) == set(cds)
    n_edit = 0
    for g in tables["genes"]:
        seq = tables["seq_of"][(g["genome_id"], int(g["replicon"]))]
        s, e, strand = int(g["start"]), int(g["end"]), int(g["strand"])
        b = bb[g["copy_label"]]
        assert (b["genome_id"], b["replicon"], b["position"], b["strand"]) == \
            (g["genome_id"], g["replicon"], g["position"], g["strand"])
        assert 0 <= s < len(seq) and e - s == len(cds[g["copy_label"]])
        assert (g["origin_spanning"] == "true") == (e > len(seq))
        got = _cut(seq, s, e, strand)
        want = cds[g["copy_label"]]
        if g["cds_edit"]:
            n_edit += 1
            for part in g["cds_edit"].split(";"):
                where, change = part.split(":")
                old, new = change.split(">")
                assert old != new
                if where == "start":
                    assert want[:3] == old and new in STARTS
                    want = new + want[3:]
                else:
                    assert where == "stop" and want[-3:] == old and new in STOPS
                    want = want[:-3] + new
        assert got == want, g["copy_label"]
        assert translate(got) == prot[g["copy_label"]]
        assert (g["start_codon"], g["stop_codon"]) == (got[:3], got[-3:])
    print(f"\n{len(tables['genes'])} genes, {n_edit} with protein-neutral overlap edits")


@needs_zombi2
def test_genes_and_igrs_tile_every_replicon_and_overlap_only_where_planted(world, tables):
    _, out, _ = world
    bb_reps = {(r["genome_id"], int(r["replicon"])): r for r in _tsv(out / "backbone/replicons.tsv")}
    assert set(tables["reps"]) == set(bb_reps)
    genes_by = collections.defaultdict(list)
    for g in tables["genes"]:
        genes_by[(g["genome_id"], int(g["replicon"]))].append(g)
    juncs_by = collections.defaultdict(list)
    for j in tables["juncs"]:
        juncs_by[(j["genome_id"], int(j["replicon"]))].append(j)
    elems_by = collections.defaultdict(list)
    for x in tables["elems"]:
        elems_by[(x["genome_id"], int(x["replicon"]))].append(x)
    n_over = collections.Counter()
    for key, rep in tables["reps"].items():
        seq = tables["seq_of"][key]
        N = len(seq)
        assert int(rep["length"]) == N and rep["topology"] == bb_reps[key]["topology"] == "circular"
        assert int(rep["n_genes"]) == int(bb_reps[key]["n_genes"]) == len(genes_by[key])
        assert rep["contig_uid"] == blake3(seq.encode()).hexdigest(length=16)
        cover = np.zeros(N, dtype=np.int64)
        gene_cover = np.zeros(N, dtype=np.int64)
        for g in genes_by[key]:
            s, e = int(g["start"]), int(g["end"])
            idx = np.arange(s, e) % N
            cover[idx] += 1
            gene_cover[idx] += 1
        over = 0
        for j in juncs_by[key]:
            s, e, L, o = int(j["start"]), int(j["end"]), int(j["length"]), int(j["overlap_bp"])
            assert e - s == L and 0 <= s < N
            igr = _cut(seq, s, e, 1) if L else ""
            if L:
                cover[np.arange(s, e) % N] += 1
                assert j["igr_uid"] == blake3(igr.encode()).hexdigest(length=16)
                assert abs(float(j["gc"]) - (igr.count("G") + igr.count("C")) / L) < 1e-5
            else:
                assert j["igr_uid"] == ""
            if o:
                assert j["junction"] == "codir_tight" and o in (1, 4) and int(j["target_length"]) == -o
                n_over[o] += 1
                over += o
        # every base is a gene's or an IGR's, exactly once, but the planted overlaps' bases (two genes)
        assert set(np.unique(cover)) <= {1, 2} and int((cover == 2).sum()) == over
        assert int((gene_cover == 2).sum()) == over
        # every element lies inside an IGR, never inside a gene
        for x in elems_by[key]:
            s, e = int(x["start"]), int(x["end"])
            assert 0 <= s < N and gene_cover[np.arange(s, e) % N].max() == 0, x
    print(f"\nplanted overlaps by length: {dict(n_over)}")
    assert n_over[4] > 0 and n_over[1] > 0


@needs_zombi2
def test_overlaps_read_their_templates_and_are_realized_exactly_when_the_codons_allow(world, tables):
    """A -4 overlap reads Z+TGA on the coding strand (the upstream stop TGA inside the downstream start ZTG
    and its next base A); a -1 overlap shares the stop's last base with the start's first. A tried -4 is
    realized iff the upstream gene's last sense codon does not end in C and the downstream gene's second
    codon starts with A (the edits change only start and stop codons); a tried -1 always is, with these
    starts."""
    _, out, _ = world
    cds = _fasta(out / "sequences/cds.fna")
    gene = {g["copy_label"]: g for g in tables["genes"]}
    tried = collections.Counter()
    done = collections.Counter()
    n_tight = 0
    for j in tables["juncs"]:
        if j["junction"] != "codir_tight":
            continue
        n_tight += 1
        a, b = gene[j["left_copy"]], gene[j["right_copy"]]
        up, down = (a, b) if int(a["strand"]) == 1 else (b, a)
        t, o = int(j["target_length"]), int(j["overlap_bp"])
        if t >= 0:
            assert o == 0
            continue
        tried[-t] += 1
        cu, cd = cds[up["copy_label"]], cds[down["copy_label"]]
        possible = (cu[-4] != "C" and cd[3] == "A") if t == -4 else cd[:3] in STARTS
        assert (o == -t) == possible, (j, cu[-6:], cd[:6])
        if o:
            done[o] += 1
            seq = tables["seq_of"][(j["genome_id"], int(j["replicon"]))]
            ds, de = int(down["start"]), int(down["end"])
            head = _cut(seq, ds, ds + 4, 1) if int(down["strand"]) == 1 else _cut(seq, de - 4, de, -1)
            tail_up = _cut(seq, int(up["start"]), int(up["end"]), int(up["strand"]))
            assert tail_up[-o:] == head[:o]
            if o == 4:
                assert head[1:] == "TGA" and tail_up[-3:] == "TGA"
    p4 = backbone.load_preset("micro")["intergenic"]["junctions"]["tight_overlap"]
    f4, f1 = tried[4] / n_tight, tried[1] / n_tight
    print(f"\n{n_tight} tight junctions: tried -4 {f4:.3f}, -1 {f1:.3f}; realized -4 {done[4] / n_tight:.3f}, "
          f"-1 {done[1] / n_tight:.3f}")
    # the tries follow Phi(z_T), z_T N(0, 1) per downstream copy (copies of one family share it closely, so
    # the tolerance allows for about 1,000 independent draws: 4 sd of 0.6 is 0.06)
    assert abs(f4 - p4["minus4"]) <= 0.06 and abs(f1 - p4["minus1"]) <= 0.04


@needs_zombi2
def test_features_across_the_origin_follow_the_convention(tables):
    n_cross = collections.Counter()
    for name, rows in (("genes", tables["genes"]), ("juncs", tables["juncs"]), ("elems", tables["elems"]),
                       ("rna", tables["rna"])):
        for r in rows:
            N = len(tables["seq_of"][(r["genome_id"], int(r["replicon"]))])
            s, e = int(r["start"]), int(r["end"])
            assert 0 <= s < N and s <= e < s + N, (name, r)
            n_cross[name] += e > N
    print(f"\nfeatures across the origin: {dict(n_cross)}")
    assert n_cross["genes"] >= 5        # §5.1: at least 5 origin-spanning genes


# --- junctions (§2.4; spike §7 item 8) ---------------------------------------------------------------------

@needs_zombi2
def test_junction_types_and_classes_follow_the_rules(world, tables):
    """Types from the strands (written here); a planted unit's inner pairs tight; a motif-carrying gene's
    upstream junction loose; every other co-directional class a function of the ordered family pair, tight
    for a share of distinct pairs within 4 sd of p_tight; declared families of different groups next to each
    other are tight no more often than background pairs (spike §7 item 8: no inherited module adjacency)."""
    _, out, _ = world
    p = backbone.load_preset("micro")
    bb = {r["copy_label"]: r for r in _tsv(out / "backbone/genes.tsv")}
    fams = {int(r["family"]): r for r in _tsv(out / "backbone/families.tsv")}
    planted = {r["name"] for r in _tsv(out / "backbone/groups.tsv") if r["kind"] == "planted"}
    motif = {r["attached_copy"] for r in tables["rna"]}
    cls_of: dict[tuple[int, int], set] = collections.defaultdict(set)
    declared_mixed = set()
    n_group = 0
    for j in tables["juncs"]:
        a, b = bb[j["left_copy"]], bb[j["right_copy"]]
        sa, sb = int(a["strand"]), int(b["strand"])
        want = ("codir" if sa == sb else "divergent" if sa == -1 else "convergent")
        assert j["junction"].startswith(want), j
        if sa != sb:
            continue
        up, down = (a, b) if sa == 1 else (b, a)
        same_unit = a["group"] in planted and a["group"] == b["group"] and a["unit"] == b["unit"]
        if same_unit:
            assert j["junction"] == "codir_tight" and j["rule"] == "group"
            n_group += 1
            continue
        if down["copy_label"] in motif:
            assert j["junction"] == "codir_loose" and j["rule"] in ("motif", "motif_min")
            continue
        assert j["rule"] == "drawn"
        key = (int(up["family"]), int(down["family"]))
        cls_of[key].add(j["junction"])
        fa, fb = fams[key[0]], fams[key[1]]
        if "" not in (fa["name"], fb["name"]) and (fa["group"] != fb["group"] or fa["group"] == ""):
            declared_mixed.add(key)      # two declared families of different groups (or singles) side by side
    assert all(len(v) == 1 for v in cls_of.values())
    n = len(cls_of)
    p_t = p["intergenic"]["junctions"]["p_tight"]
    frac = sum("codir_tight" in v for v in cls_of.values()) / n
    sd = math.sqrt(p_t * (1 - p_t) / n)
    m = len(declared_mixed)
    mixed = sum("codir_tight" in cls_of[k] for k in declared_mixed) / max(1, m)
    print(f"\n{n_group} within-unit junctions; {n} distinct background co-directional pairs, tight {frac:.4f}; "
          f"{m} distinct pairs of declared families of different groups, tight {mixed:.3f}")
    assert n_group > 0 and abs(frac - p_t) <= 4 * sd
    assert m >= 10 and mixed <= p_t + 4 * math.sqrt(p_t * (1 - p_t) / m)


def _quantile(dist: dict, q: float) -> float:
    z = statistics.NormalDist().inv_cdf(q)
    return min(max(round(dist["offset"] + dist["scale"] * math.exp(dist["sigma"] * z)), dist["min"]), dist["max"])


@needs_zombi2
def test_junction_lengths_follow_the_preset_by_type(tables):
    """Realized lengths by junction type (drawn lengths: no motif raised them; tight ones without a tried
    overlap) against the preset's quantiles at 0.1 ... 0.9, within max(3 bases, 12%) (derived: copies of one
    family carry closely inherited length scores, so the effective sample is about the number of distinct
    family pairs, over a thousand per type; a quantile's standard error is then below 5% of its value)."""
    J = backbone.load_preset("micro")["intergenic"]["junctions"]
    by = collections.defaultdict(list)
    for j in tables["juncs"]:
        if j["rule"] not in ("drawn", "group"):
            continue
        t = int(j["target_length"])
        if j["junction"] == "codir_tight" and t < 0:
            continue
        by[j["junction"]].append(int(j["length"]))
    report = {}
    for kind, dist in (("codir_tight", J["tight"]), ("codir_loose", J["loose"]), ("divergent", J["divergent"]),
                       ("convergent", J["convergent"])):
        xs = np.array(by[kind])
        assert len(xs) >= 500, (kind, len(xs))
        rows = []
        for q in (0.1, 0.25, 0.5, 0.75, 0.9):
            got, want = float(np.quantile(xs, q)), _quantile(dist, q)
            rows.append((q, got, want))
            assert abs(got - want) <= max(3.0, 0.12 * want), (kind, rows)
        report[kind] = (len(xs), [(g, w) for _, g, w in rows])
    print(f"\njunction length quantiles (realized, preset) by type: {report}")


# --- composition, elements, conservation (§2.4) ------------------------------------------------------------

@needs_zombi2
def test_igr_composition_follows_the_lineage(world, tables):
    _, out, _ = world
    node_of = {r["genome_id"]: r["species_node"] for r in _tsv(out / "backbone/genomes.tsv")}
    gc_igr = {r["lineage"]: float(r["gc_igr"]) for r in _tsv(out / "replicons/lineages.tsv")}
    gc, at = collections.Counter(), collections.Counter()
    for j in tables["juncs"]:
        L = int(j["length"])
        if not L:
            continue
        s = _cut(tables["seq_of"][(j["genome_id"], int(j["replicon"]))], int(j["start"]), int(j["end"]), 1)
        gc[j["genome_id"]] += s.count("G") + s.count("C")
        at[j["genome_id"]] += L
    gids = sorted(at)
    x = np.array([gc_igr[node_of[g]] for g in gids])
    y = np.array([gc[g] / at[g] for g in gids])
    r = float(np.corrcoef(x, y)[0, 1])
    print(f"\ngenome IGR GC vs lineage intergenic GC: r {r:.3f}, lineage {x.min():.3f}-{x.max():.3f}, "
          f"genome {y.min():.3f}-{y.max():.3f}")
    assert x.max() - x.min() >= 0.05 and r >= 0.6


def _terminator_ok(t: str) -> bool:
    """Some stem k in 6..10 and loop in 3..6 with the 7-base T tract after the closing stem, every stem pair
    canonical."""
    for k in range(6, 11):
        loop = len(t) - 7 - 2 * k
        if 3 <= loop <= 6 and all((t[i], t[2 * k + loop - 1 - i]) in PAIRS_OK for i in range(k)):
            return True
    return False


@needs_zombi2
def test_planted_elements_read_as_their_motifs(tables):
    """Boxes and RBSs on their strand match their consensus at about the column probability (0.7), far above
    the IGR background; every terminator has a canonical stem, a loop and a T tract."""
    p = backbone.load_preset("micro")["intergenic"]["cassettes"]["element_consensus"]
    cons = {"promoter_m35": "TTGACA", "promoter_m10": "TATAAT", "rbs": "AGGAGG"}
    match = collections.defaultdict(list)
    n_term, t_frac = 0, []
    for x in tables["elems"]:
        seq = tables["seq_of"][(x["genome_id"], int(x["replicon"]))]
        s = _cut(seq, int(x["start"]), int(x["end"]), int(x["strand"]))
        if x["kind"] in cons:
            assert len(s) == 6
            match[x["kind"]].append(sum(a == b for a, b in zip(s, cons[x["kind"]])) / 6)
        elif x["kind"] == "terminator":
            n_term += 1
            assert _terminator_ok(s), s
            t_frac.append(s[-7:].count("T") / 7)
    # background 6-mers: every 50th junction's first 6 bases, on +
    bg = collections.defaultdict(list)
    for j in tables["juncs"][::50]:
        if int(j["length"]) >= 6:
            s = _cut(tables["seq_of"][(j["genome_id"], int(j["replicon"]))], int(j["start"]), int(j["start"]) + 6, 1)
            for k, c in cons.items():
                bg[k].append(sum(a == b for a, b in zip(s, c)) / 6)
    report = {k: (len(v), round(float(np.mean(v)), 3), round(float(np.mean(bg[k])), 3)) for k, v in match.items()}
    print(f"\nelement consensus match (n, planted, background): {report}; terminators {n_term}, "
          f"T tract {np.mean(t_frac):.3f}")
    for k, v in match.items():
        assert len(v) >= 1000 and abs(np.mean(v) - p) <= 0.03 and np.mean(bg[k]) <= 0.4
    assert n_term >= 1000 and abs(np.mean(t_frac) - p) <= 0.03


@needs_zombi2
def test_orthologous_igrs_are_near_identical_in_close_relatives(world, tables):
    """The same junction (left family, right family, type) in the clonal burst's genomes against genomes
    of different phyla: equal-length IGRs share >= 0.9 of their bases within the burst, and at least 0.15
    more than across phyla."""
    _, out, _ = world
    bbj = json.loads((out / "backbone/backbone.json").read_text(encoding="utf-8"))
    node_gid = {r["species_node"]: r["genome_id"] for r in _tsv(out / "backbone/genomes.tsv")}
    burst = {node_gid[t] for t in bbj["burst_tips"]}
    phylum = {r["genome_id"]: r["tax_lineage"].split(";")[1] for r in _tsv(out / "backbone/taxonomy.tsv")}
    fam = {g["copy_label"]: g["family_id"] for g in tables["genes"]}
    by_key = collections.defaultdict(dict)
    for j in tables["juncs"]:
        if int(j["length"]) < 20:
            continue
        key = (fam[j["left_copy"]], fam[j["right_copy"]], j["junction"])
        by_key[key].setdefault(j["genome_id"], j)

    def ident(j1, j2):
        a = _cut(tables["seq_of"][(j1["genome_id"], int(j1["replicon"]))], int(j1["start"]), int(j1["end"]), 1)
        b = _cut(tables["seq_of"][(j2["genome_id"], int(j2["replicon"]))], int(j2["start"]), int(j2["end"]), 1)
        return sum(x == y for x, y in zip(a, b)) / len(a) if len(a) == len(b) else None

    close, far = [], []
    for key, per in by_key.items():
        gs = sorted(per)
        bs = [g for g in gs if g in burst]
        if len(bs) >= 2:
            v = ident(per[bs[0]], per[bs[1]])
            if v is not None:
                close.append(v)
        g1 = gs[0]
        g2 = next((g for g in gs if phylum[g] != phylum[g1]), None)
        if g2 is not None:
            v = ident(per[g1], per[g2])
            if v is not None:
                far.append(v)
    print(f"\northologous IGR identity: burst {np.mean(close):.3f} over {len(close)}, across phyla "
          f"{np.mean(far):.3f} over {len(far)}")
    assert len(close) >= 200 and len(far) >= 200
    assert np.mean(close) >= 0.9 and np.mean(close) >= np.mean(far) + 0.15


# --- RNA motifs (§2.5) ------------------------------------------------------------------------------------

def _hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


@needs_zombi2
@pytest.mark.covers("F-02", "B-SIM")
def test_every_planted_motif_is_recovered_at_its_truth_coordinates_by_an_independent_scan(world, tables):
    """Each instance cut at its coordinates differs from its clade's emission at exactly its
    n_subst_from_emission positions, n_compensatory pairs differ at both bases, every pair is canonical;
    dist_to_start is its distance to its gene's start codon. Scanning both strands of the store's own
    intergenic region around it (`ncrna.regions`, IGR runs widened by 200) for the emission, by Hamming
    distance, finds the instance at its coordinates as the best window."""
    _, out, _ = world
    em = {r["emission_id"]: r for r in _tsv(out / "replicons/emissions.tsv")}
    links = {r["link"]: r for r in _tsv(out / "replicons/links.tsv")}
    gene = {g["copy_label"]: g for g in tables["genes"]}
    carriers = collections.Counter(g["family_id"] for g in tables["genes"])
    per_link = collections.Counter(r["link"] for r in tables["rna"])
    for lk, row in links.items():
        assert per_link[lk] == carriers[row["family_id"]] > 0
    genes_by = collections.defaultdict(list)
    for g in tables["genes"]:
        genes_by[(g["genome_id"], int(g["replicon"]))].append((int(g["start"]), int(g["end"])))
    regions = {k: l1_ncrna.regions(v, len(tables["seq_of"][k]), True) for k, v in genes_by.items()}
    best, n_div = 0, 0
    divs = []
    for r in tables["rna"]:
        key = (r["genome_id"], int(r["replicon"]))
        seq = tables["seq_of"][key]
        N = len(seq)
        s, e, strand = int(r["start"]), int(r["end"]), int(r["strand"])
        inst = _cut(seq, s, e, strand)
        ref = em[r["emission_id"]]
        assert len(inst) == len(ref["sequence"]) == int(r["length"])
        assert _hamming(inst, ref["sequence"]) == int(r["n_subst_from_emission"])
        pairs = [(i, j) for i, j in _pairs(ref["structure"])]
        assert r["structure"] == ref["structure"]
        assert sum(inst[i] != ref["sequence"][i] and inst[j] != ref["sequence"][j] for i, j in pairs) == \
            int(r["n_compensatory"])
        assert all((inst[i], inst[j]) in PAIRS_OK for i, j in pairs if (ref["sequence"][i], ref["sequence"][j]) in PAIRS_OK)
        g = gene[r["attached_copy"]]
        assert int(g["strand"]) == strand and g["family_id"] == r["family_id"]
        gs, ge = int(g["start"]), int(g["end"])
        dist = (gs - e) % N if strand == 1 else (s - ge) % N
        assert dist == int(r["dist_to_start"])
        divs.append(int(r["n_subst_from_emission"]) / len(inst))
        # the independent scan, both strands of the region that holds the instance
        reg = next((a, b) for a, b in regions[key] if a <= s and e <= b or a <= s + N and e + N <= b)
        a, b = reg
        span = _cut(seq, a, b, 1)
        L = len(inst)
        want = np.frombuffer(ref["sequence"].encode(), np.uint8)
        fwd = np.frombuffer(span.encode(), np.uint8)
        rev = np.frombuffer(rc(span).encode(), np.uint8)
        hf = (np.lib.stride_tricks.sliding_window_view(fwd, L) != want).sum(axis=1)
        hr = (np.lib.stride_tricks.sliding_window_view(rev, L) != want).sum(axis=1)
        low = int(min(hf.min(), hr.min()))
        # a window k of the reverse complement covers [len - k - L, len - k) of the span
        winners = [((a + int(k)) % N, 1) for k in np.nonzero(hf == low)[0]]
        winners += [((a + len(span) - int(k) - L) % N, -1) for k in np.nonzero(hr == low)[0]]
        assert (s, strand) in winners or divs[-1] > 0.3, (r, low, winners[:3])
        best += winners == [(s, strand)]
        n_div += divs[-1] > 0.3
    n = len(tables["rna"])
    print(f"\n{n} motif instances in {len(links)} links: the truth window is the unique best scan window for "
          f"{best}; divergence from the emission median {np.median(divs):.3f}, max {max(divs):.3f}, "
          f"{n_div} above 0.3")
    assert n >= 200 and best >= 0.95 * n


def _pairs(db: str):
    st, out = [], []
    for i, c in enumerate(db):
        if c == "(":
            st.append(i)
        elif c == ")":
            out.append((st.pop(), i))
    assert not st
    return sorted(out)


@needs_zombi2
@pytest.mark.covers("F-02", "R-05")
def test_the_store_igr_code_finds_the_planted_igrs_and_motifs(tables):
    """`biobox.l1.igr.igrs_of_contig` on the genes' coordinates gives exactly the junctions of positive
    length, and every planted motif lies inside one of those runs."""
    genes_by = collections.defaultdict(list)
    for g in tables["genes"]:
        genes_by[(g["genome_id"], int(g["replicon"]))].append((int(g["start"]), int(g["end"])))
    juncs_by = collections.defaultdict(set)
    for j in tables["juncs"]:
        if int(j["length"]):
            juncs_by[(j["genome_id"], int(j["replicon"]))].add((int(j["start"]), int(j["end"])))
    runs_by = {}
    for key, gs in genes_by.items():
        N = len(tables["seq_of"][key])
        runs = l1_igr.igrs_of_contig([(i, s, e) for i, (s, e) in enumerate(gs)], N, True)
        runs_by[key] = [(s, e) for s, e, _, _ in runs]
        assert set(runs_by[key]) == juncs_by[key], key
    for r in tables["rna"]:
        key = (r["genome_id"], int(r["replicon"]))
        N = len(tables["seq_of"][key])
        s, e = int(r["start"]), int(r["end"])
        assert any(a <= s and e <= b or a <= s + N and e + N <= b for a, b in runs_by[key]), r


# --- runtime (§2.6) ---------------------------------------------------------------------------------------

@needs_zombi2
def test_micro_replicons_runtime_within_budget(world, record_property):
    rw, _, seconds = world
    record_property("replicons_seconds", round(seconds, 2))
    c = rw.meta["counts"]
    print(f"\nF-02 micro replicons: {seconds:.2f} s wall; {c}; digest {rw.digest}")
    # synthworld_spec §2.6: all of synthworld within 3 min; backbone and sequences took about 115 s
    # (measured, session 3), so this stage gets 60 s
    assert seconds < 60.0
