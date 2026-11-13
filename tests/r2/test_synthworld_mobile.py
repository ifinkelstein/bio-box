"""Ring R2: synthworld phages, prophages and CRISPR arrays on the `micro` preset (synthworld_spec §2.1 S6,
§3.3 truth_spacer_target, §5.1, §6.1, §7.3; card F-02 session 5).

A module fixture builds the micro world's S1-S5 (ZOMBI2 0.54.0) and its mobile stage once. The oracles are
written here: FASTA is parsed, reverse complements and translation are spelled out, the lift of S5
coordinates is recomputed from insertions.tsv by this file's own formula, clades are read from the species
tree's Newick, repeats are counted by string search, and spacer targets are found by a brute-force Hamming
scan over every window. What they take from the generator are model inputs (the preset, the S5 world, a
spacer's sequence), never its realized statistics.
"""

from __future__ import annotations

import collections
import csv
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest
from blake3 import blake3

from biobox.synthworld import backbone, intergenic, mobile, sequences

pytestmark = [pytest.mark.r2, pytest.mark.covers("F-02")]
needs_zombi2 = pytest.mark.skipif(
    importlib.util.find_spec("zombi2") is None,
    reason="zombi2 0.54.0 is not importable here: CI installs it by `uv sync`; on the pod the image has it")

CODE11_BY_AA = {
    "F": "TTT TTC", "L": "TTA TTG CTT CTC CTA CTG", "I": "ATT ATC ATA", "M": "ATG", "V": "GTT GTC GTA GTG",
    "S": "TCT TCC TCA TCG AGT AGC", "P": "CCT CCC CCA CCG", "T": "ACT ACC ACA ACG", "A": "GCT GCC GCA GCG",
    "Y": "TAT TAC", "*": "TAA TAG TGA", "H": "CAT CAC", "Q": "CAA CAG", "N": "AAT AAC", "K": "AAA AAG",
    "D": "GAT GAC", "E": "GAA GAG", "C": "TGT TGC", "W": "TGG", "R": "CGT CGC CGA CGG AGA AGG",
    "G": "GGT GGC GGA GGG",
}
CODE11 = {c: aa for aa, cs in CODE11_BY_AA.items() for c in cs.split()}
STARTS = {"ATG", "GTG", "TTG"}
_RC = {"A": "T", "C": "G", "G": "C", "T": "A"}


def rc(s: str) -> str:
    return "".join(_RC[c] for c in reversed(s))


def translate(cds: str) -> str:
    aa = [CODE11[cds[i:i + 3]] for i in range(0, len(cds), 3)]
    assert cds[:3] in STARTS and aa[-1] == "*" and "*" not in aa[:-1], cds[:12]
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


def _cut(seq: str, start: int, end: int, strand: int = 1) -> str:
    n = len(seq)
    assert 0 <= start <= end <= start + n
    s = (seq + seq)[start:end] if end > n else seq[start:end]
    return s if strand == 1 else rc(s)


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    if importlib.util.find_spec("zombi2") is None:
        pytest.skip("zombi2 0.54.0 is not importable here")
    out = tmp_path_factory.mktemp("micro_mobile") / "world"
    params = backbone.load_preset("micro")
    backbone.build(params, params["default_seed"], out)
    sequences.build(out)
    intergenic.build(out)
    t0 = time.perf_counter()
    mw = mobile.build(out)
    return mw, out, time.perf_counter() - t0


@pytest.fixture(scope="module")
def t(world):
    _, out, _ = world
    md = out / "mobile"
    seqs = {}
    for f in sorted((md / "fasta").glob("*.fna")):
        seqs.update(_fasta(f))
    s5 = {}
    for f in sorted((out / "replicons" / "fasta").glob("*.fna")):
        s5.update(_fasta(f))
    return {"seqs": seqs, "s5": s5, "phages": _fasta(md / "phages.fna"),
            "genes": _tsv(md / "genes.tsv"), "ins": _tsv(md / "insertions.tsv"),
            "pro": _tsv(md / "prophages.tsv"), "arrays": _tsv(md / "arrays.tsv"), "spacers": _tsv(md / "spacers.tsv"),
            "targets": _tsv(md / "spacer_targets.tsv"), "cas": _tsv(md / "cas.tsv"),
            "phage_rows": _tsv(md / "phages.tsv"), "phage_genes": _tsv(md / "phage_genes.tsv"),
            "reps": _tsv(md / "replicons.tsv")}


def _species(out: Path):
    """Parent, children and extant tips of the complete species tree, from its Newick (parsed here)."""
    text = (out / "zombi2" / "species_complete.nwk").read_text(encoding="utf-8").strip()
    parent, stack, last, label_next = {}, [], None, False
    nodes: list[str | None] = []
    tokens = []
    buf = ""
    for ch in text:
        if ch in "(),;":
            if buf:
                tokens.append(buf)
                buf = ""
            tokens.append(ch)
        else:
            buf += ch
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
            if label_next:                 # the label of the clade just closed
                nodes[last] = name
                k = last
            else:                          # a leaf
                nodes.append(name)
                k = len(nodes) - 1
            parent[k] = stack[-1] if stack else None
            label_next = False
    return {nodes[k]: (nodes[p] if p is not None else None) for k, p in parent.items()}


def _below(par: dict, node: str, tips: set[str]) -> set[str]:
    out = set()
    for tip in tips:
        x = tip
        while x is not None:
            if x == node:
                out.add(tip)
                break
            x = par[x]
    return out


def _ancestors(par: dict, tip: str) -> set[str]:
    out, x = set(), tip
    while x is not None:
        out.add(x)
        x = par[x]
    return out


# --- determinism (§6.1, §6.2; MR-16) -----------------------------------------------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "MR-16")
def test_same_seed_gives_byte_identical_mobile_stage_in_a_subprocess(world, tmp_path, record_property):
    mw, out, _ = world
    dst = tmp_path / "copy"
    dst.mkdir()
    for d in ("backbone", "zombi2", "sequences", "replicons"):
        shutil.copytree(out / d, dst / d)
    env = dict(os.environ, PYTHONHASHSEED="4242", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    code = ("import os, sys, time; from biobox.synthworld import mobile; t0 = time.perf_counter(); "
            "mobile.build(sys.argv[1]); dt = time.perf_counter() - t0; p = '/proc/self/status'; "
            "hwm = next((l.split()[1] for l in open(p) if l.startswith('VmHWM')), -1) if os.path.exists(p) else -1; "
            "print(hwm, round(dt, 2))")
    proc = subprocess.run([sys.executable, "-c", code, str(dst)], check=True, env=env, cwd=tmp_path,
                          capture_output=True, text=True)
    assert _tree_sha(dst / "mobile") == _tree_sha(out / "mobile")
    hwm_kb, secs = proc.stdout.split()[-2:]
    print(f"\nmobile stage in a child: {secs} s wall, VmHWM {int(hwm_kb) / 1024:.0f} MB")
    record_property("mobile_child_seconds", float(secs))
    if int(hwm_kb) < 0:
        return
    rss_mb = int(hwm_kb) / 1024
    record_property("mobile_child_maxrss_mb", round(rss_mb, 1))
    assert rss_mb < 1024      # one host genome at a time (AGENTS §4.1)


# --- the phage set (§2.1 S6, §5.1) -------------------------------------------------------------------------

@needs_zombi2
def test_phage_set_is_twelve_linear_phages_whose_genes_encode_their_proteins(world, t):
    _, out, _ = world
    pw = out / "mobile" / "phages"
    p = backbone.load_preset("micro")
    assert len(t["phage_rows"]) == p["phages"]["count"] == len(t["phages"])
    assert {r["topology"] for r in _tsv(pw / "backbone" / "replicons.tsv")} == {"linear"}
    assert {r["topology"] for r in _tsv(pw / "replicons" / "replicons.tsv")} == {"linear"}
    prot = _fasta(pw / "sequences" / "proteins.faa")
    sizes = collections.Counter(r["phage_id"] for r in t["phage_genes"])
    fams = collections.defaultdict(set)
    for g in t["phage_genes"]:
        seq = t["phages"][g["phage_id"]]
        s, e, strand = int(g["start"]), int(g["end"]), int(g["strand"])
        assert 0 <= s < e <= len(seq)
        assert translate(_cut(seq, s, e, strand)) == prot[g["copy_label"][3:]]
        assert g["protein_uid"] == blake3(prot[g["copy_label"][3:]].encode()).hexdigest(length=16)
        fams[g["family_id"]].add(g["phage_id"])
    for r in t["phage_rows"]:
        seq = t["phages"][r["record"]]
        assert int(r["length"]) == len(seq) and r["contig_uid"] == blake3(seq.encode()).hexdigest(length=16)
        assert int(r["n_genes"]) == sizes[r["phage_id"]]
    shared = sum(len(v) >= 2 for v in fams.values()) / len(fams)
    print(f"\nphage genes per phage {sorted(sizes.values())}; lengths "
          f"{sorted(int(r['length']) for r in t['phage_rows'])}; {len(fams)} families, {shared:.2f} in >= 2 phages")
    assert all(40 <= n <= 80 for n in sizes.values())        # §5.1: 40-80 genes each
    assert shared >= 0.2                                       # a related set
    # the phage set takes its domains from the host's HMM panel (one panel per world)
    assert (pw / "sequences" / "panel.tsv").read_bytes() == (out / "sequences" / "panel.tsv").read_bytes()


# --- the host replicons after the insertions -----------------------------------------------------------------

def _lift_rows(ins_rows):
    by = collections.defaultdict(list)
    for r in ins_rows:
        by[r["record"]].append((int(r["at"]), int(r["length"])))
    return by


def _lift(points, N, s, e):
    """This file's lift: blocks at points a (N for an insertion at the origin) of lengths l."""
    s2 = s + sum(l for a, l in points if a <= s)
    if e <= N:
        return s2, e + sum(l for a, l in points if a < e)
    return s2, N + sum(l for _, l in points) + (e - N)


@needs_zombi2
def test_every_host_feature_keeps_its_bases_after_the_insertions(world, t):
    """Every host gene of S5 (joined by copy_label), every S5 element and motif, and every S5 junction lifted
    by insertions.tsv reads the same bases in the new replicon as in S5; host genes keep their order."""
    _, out, _ = world
    pts = _lift_rows(t["ins"])
    s5_genes = {g["copy_label"]: g for g in _tsv(out / "replicons" / "genes.tsv")}
    host = [g for g in t["genes"] if g["source"] == "host"]
    assert {g["copy_label"] for g in host} == set(s5_genes)
    reps = {r["record"]: r for r in t["reps"]}
    for g in host:
        old = s5_genes[g["copy_label"]]
        rec = f"{g['genome_id']}_r{g['replicon']}"
        N = len(t["s5"][rec])
        assert int(reps[rec]["s5_length"]) == N
        want = _lift(pts[rec], N, int(old["start"]), int(old["end"]))
        assert (int(g["start"]), int(g["end"])) == want
        assert _cut(t["seqs"][rec], *want) == _cut(t["s5"][rec], int(old["start"]), int(old["end"]))
        assert (g["origin_spanning"] == "true") == (want[1] > len(t["seqs"][rec]))
    n = 0
    for name in ("elements.tsv", "rna.tsv", "junctions.tsv"):
        for r in _tsv(out / "replicons" / name):
            rec = f"{r['genome_id']}_r{r['replicon']}"
            s, e = int(r["start"]), int(r["end"])
            if e == s:
                continue
            N = len(t["s5"][rec])
            assert _cut(t["seqs"][rec], *_lift(pts[rec], N, s, e)) == _cut(t["s5"][rec], s, e)
            n += 1
    for rec, seq in t["seqs"].items():
        added = sum(l for _, l in pts[rec])
        assert len(seq) == len(t["s5"][rec]) + added
        assert reps[rec]["contig_uid"] == blake3(seq.encode()).hexdigest(length=16)
    print(f"\n{len(host)} host genes and {n} S5 elements, motifs and IGRs read the same after "
          f"{len(t['ins'])} insertions")


@needs_zombi2
def test_genes_across_the_origin_follow_the_convention(t):
    n = 0
    for g in t["genes"]:
        N = len(t["seqs"][f"{g['genome_id']}_r{g['replicon']}"])
        s, e = int(g["start"]), int(g["end"])
        assert 0 <= s < N and s < e < s + N
        n += e > N
    assert n >= 5       # §5.1: at least 5 origin-spanning genes


# --- prophages ---------------------------------------------------------------------------------------------

@needs_zombi2
def test_prophages_sit_in_exactly_their_clade_after_their_anchor_and_encode_the_phage(world, t):
    _, out, _ = world
    par = _species(out)
    node_gid = {r["species_node"]: r["genome_id"] for r in _tsv(out / "backbone" / "genomes.tsv")}
    tips = set(node_gid)
    p = backbone.load_preset("micro")
    pw = out / "mobile" / "phages"
    prot = _fasta(pw / "sequences" / "proteins.faa")
    pgenes = collections.defaultdict(list)
    for g in t["phage_genes"]:
        pgenes[g["phage_id"]].append(g)
    genes_by_ins = collections.defaultdict(list)
    for g in t["genes"]:
        if g["insertion_id"]:
            genes_by_ins[g["insertion_id"]].append(g)
    gene_at = {g["copy_label"]: g for g in t["genes"]}
    by_pro = collections.defaultdict(list)
    for r in t["pro"]:
        by_pro[r["prophage_id"]].append(r)
    assert len(by_pro) == p["prophages"]["count"]
    idents = []
    for pid, rows in by_pro.items():
        node = rows[0]["clade_node"]
        want = {node_gid[x] for x in _below(par, node, tips)}
        assert {r["genome_id"] for r in rows} == want
        lo, hi = p["prophages"]["clade_size"]
        assert lo <= len(want) <= hi
        phage = t["phages"][rows[0]["phage_id"]]
        for r in rows:
            seq = t["seqs"][r["record"]]
            s, e, strand = int(r["start"]), int(r["end"]), int(r["strand"])
            inst = _cut(seq, s, e, strand)
            assert len(inst) == len(phage)
            ident = sum(a == b for a, b in zip(inst, phage)) / len(phage)
            assert abs(ident - float(r["identity"])) < 1e-5 and ident >= 0.9
            idents.append(ident)
            # right after the anchor's stop codon
            ins = next(x for x in t["ins"] if x["insertion_id"] == f"{pid}@{r['genome_id']}")
            anc = gene_at[ins["anchor_copy"]]
            assert anc["family_id"] == ins["anchor_family"] == r["anchor_family"]
            if int(anc["strand"]) == 1:
                assert int(anc["end"]) % len(seq) == s
            else:
                assert int(anc["start"]) == e % len(seq)
            # every phage gene, at its place in the instance, encodes the phage's protein
            got = genes_by_ins[ins["insertion_id"]]
            assert len(got) == len(pgenes[r["phage_id"]])
            for g in got:
                gs, ge = int(g["start"]), int(g["end"])
                assert s <= gs < ge <= e
                label = g["copy_label"].split("@")[0][3:]
                assert translate(_cut(seq, gs, ge, int(g["strand"]))) == prot[label]
    print(f"\n{len(t['pro'])} prophage instances of {len(by_pro)} prophages; identity to the phage "
          f"{min(idents):.4f}-{max(idents):.4f}")
    assert min(idents) < 1.0      # prophages diverge down their clade


# --- CRISPR-Cas loci ------------------------------------------------------------------------------------------

def _count(seq: str, word: str) -> int:
    n, i = 0, seq.find(word)
    while i >= 0:
        n += 1
        i = seq.find(word, i + 1)
    return n


@needs_zombi2
def test_crispr_loci_read_cas_genes_leader_and_repeat_spacer_units_in_their_clade(world, t):
    _, out, _ = world
    par = _species(out)
    node_gid = {r["species_node"]: r["genome_id"] for r in _tsv(out / "backbone" / "genomes.tsv")}
    tips = set(node_gid)
    cr = backbone.load_preset("micro")["crispr"]
    spacer = {r["spacer_id"]: r for r in t["spacers"]}
    genes_by_ins = collections.defaultdict(list)
    for g in t["genes"]:
        if g["insertion_id"]:
            genes_by_ins[g["insertion_id"]].append(g)
    by_arr = collections.defaultdict(list)
    for r in t["arrays"]:
        by_arr[r["array_id"]].append(r)
    assert len(by_arr) == cr["arrays"]
    n_sp = []
    for aid, rows in by_arr.items():
        want = {node_gid[x] for x in _below(par, rows[0]["clade_node"], tips)}
        assert {r["genome_id"] for r in rows} == want
        assert cr["clade_size"][0] <= len(want) <= cr["clade_size"][1]
        for r in rows:
            seq = t["seqs"][r["record"]]
            strand = int(r["strand"])
            R = r["repeat"]
            k = cr["repeat_stem"]
            assert len(R) == cr["repeat_length"] and R[:k] == rc(R[-k:])
            ids = r["spacer_ids"].split(";")
            n = int(r["n_spacers"])
            assert len(ids) == n and cr["min_spacers"] <= n <= cr["max_spacers"]
            n_sp.append(n)
            arr = _cut(seq, int(r["array_start"]), int(r["array_end"]), strand)
            assert arr == "".join(R + spacer[x]["seq"] for x in ids) + R
            # the repeat occurs in the replicon only in this array (either strand), n + 1 times
            assert _count(seq, R) + _count(seq, rc(R)) == n + 1
            # the leader sits on the array's 5' side, between the Cas operon and the array
            ls, le, as_, ae = (int(r[c]) for c in ("leader_start", "leader_end", "array_start", "array_end"))
            assert le - ls == cr["leader_length"]
            assert (le == as_) if strand == 1 else (ae == ls)
            cas = sorted(genes_by_ins[f"{aid}@{r['genome_id']}"], key=lambda g: int(g["start"]) * strand)
            assert [g["copy_label"].split(".")[-1] for g in cas] == list(cr["cas_roles"])
            for g in cas:
                gs, ge = int(g["start"]), int(g["end"])
                assert int(g["strand"]) == strand and int(r["locus_start"]) <= gs < ge <= int(r["locus_end"])
                assert (ge <= ls) if strand == 1 else (gs >= le)
                p = translate(_cut(seq, gs, ge, strand))
                assert g["protein_uid"] == blake3(p.encode()).hexdigest(length=16)
    print(f"\n{len(t['arrays'])} array instances of {len(by_arr)} arrays; spacers per instance {min(n_sp)}-{max(n_sp)}, "
          f"median {int(np.median(n_sp))}")


@needs_zombi2
def test_spacers_are_acquired_at_the_leader_end_down_the_clade(world, t):
    """Within an instance spacer ids (acquisition order) fall from the leader; a spacer's acquisition branch is
    on the instance genome's own path below the seeding branch; instances order shared spacers alike."""
    _, out, _ = world
    par = _species(out)
    gid_node = {r["genome_id"]: r["species_node"] for r in _tsv(out / "backbone" / "genomes.tsv")}
    spacer = {r["spacer_id"]: r for r in t["spacers"]}
    by_arr = collections.defaultdict(list)
    for r in t["arrays"]:
        by_arr[r["array_id"]].append(r)
    n_branch = 0
    for aid, rows in by_arr.items():
        lists = {}
        for r in rows:
            ids = r["spacer_ids"].split(";")
            k = [int(x.split(".s")[1]) for x in ids]
            assert all(a > b for a, b in zip(k, k[1:]))
            anc = _ancestors(par, gid_node[r["genome_id"]])
            assert r["clade_node"] in anc
            for x in ids:
                node = spacer[x]["acquisition_node"]
                assert node in anc
                n_branch += node != r["clade_node"]
            times = [float(spacer[x]["acquisition_time"]) for x in ids]
            assert all(a >= b for a, b in zip(times, times[1:]))
            lists[r["genome_id"]] = ids
        for a in lists.values():
            for b in lists.values():
                sh = set(a) & set(b)
                assert [x for x in a if x in sh] == [x for x in b if x in sh]
    assert n_branch > 0


# --- spacer targets (§3.3, §7.3) -----------------------------------------------------------------------------

def _windows(seq: str, q: str, max_mm: int):
    """Every window of seq within max_mm of q on either strand, by brute force: (start, strand, mismatches)."""
    n = len(q)
    if len(seq) < n:
        return set()
    a = np.frombuffer(seq.encode(), np.uint8)
    v = np.lib.stride_tricks.sliding_window_view(a, n)
    out = set()
    for strand, want in ((1, q), (-1, rc(q))):
        d = (v != np.frombuffer(want.encode(), np.uint8)).sum(axis=1)
        out |= {(int(w), strand, int(d[w])) for w in np.flatnonzero(d <= max_mm)}
    return out


@needs_zombi2
@pytest.mark.covers("F-02", "B-SIM")
def test_designed_spacers_match_their_protospacers_and_a_brute_force_scan_finds_every_target(world, t):
    cr = backbone.load_preset("micro")["crispr"]
    max_mm = cr["target_max_mismatches"]
    by_spacer = collections.defaultdict(set)
    designed = {}
    for r in t["targets"]:
        if r["target_kind"] == "none":
            by_spacer[r["spacer_id"]] |= set()
            continue
        key = (r["target_kind"], r["target_id"], int(r["proto_start"]), int(r["proto_strand"]), int(r["n_mismatches"]))
        by_spacer[r["spacer_id"]].add(key)
        if r["designed"] == "true":
            designed[r["spacer_id"]] = key
        # each target reads as the spacer but at n_mismatches positions, and its PAM is what follows it
        seq = t["phages"][r["target_record"]] if r["target_kind"] == "phage" else t["seqs"][r["target_record"]]
        s, e, st = int(r["proto_start"]), int(r["proto_end"]), int(r["proto_strand"])
        proto = _cut(seq, s, e, st)
        assert sum(a != b for a, b in zip(proto, r["spacer_seq"])) == int(r["n_mismatches"]) <= max_mm
        if r["target_kind"] == "phage":
            pam = seq[e:e + 3] if st == 1 else (rc(seq[s - 3:s]) if s >= 3 else "")
            if st == 1 and e + 3 > len(seq):
                pam = ""
        else:
            pam = _cut(seq + seq, e, e + 3, 1) if st == 1 else rc(_cut(seq + seq, s - 3 + len(seq), s + len(seq), 1))
        assert r["pam_seq"] == pam and (r["pam_present"] == "true") == (len(pam) == 3 and pam[1:] == "GG")
    in_arrays = {x for r in t["arrays"] for x in r["spacer_ids"].split(";")}
    assert {r["spacer_id"] for r in t["targets"]} == in_arrays     # spacers.tsv also keeps the ones lost
    for sp in t["spacers"]:
        if sp["source"] == "phage" and sp["spacer_id"] in in_arrays:
            d = designed[sp["spacer_id"]]
            assert d[:2] == ("phage", sp["phage_id"]) and d[2:] == (int(sp["proto_start"]), int(sp["proto_strand"]),
                                                                   int(sp["n_mismatches"]))
    # a brute-force scan of the phages and the prophage instances, for about 60 of the spacers in an array
    present = sorted({r["spacer_id"] for r in t["targets"]})
    seqs = {r["spacer_id"]: r["spacer_seq"] for r in t["targets"]}
    blocks = [(r["prophage_id"] + "@" + r["genome_id"], r["record"], int(r["start"]), int(r["end"])) for r in t["pro"]]
    n_checked = n_targets = 0
    for sid in present[::max(1, len(present) // 60)]:
        q = seqs[sid]
        want = set()
        for pid, seq in t["phages"].items():
            want |= {("phage", pid, w, st, d) for w, st, d in _windows(seq, q, max_mm)}
        for iid, rec, s, e in blocks:
            want |= {("prophage", iid, s + w, st, d) for w, st, d in _windows(t["seqs"][rec][s:e], q, max_mm)}
        assert by_spacer[sid] == want, sid
        n_checked += 1
        n_targets += len(want)
    kinds = collections.Counter(r["target_kind"] for r in t["targets"])
    print(f"\n{len(t['spacers'])} spacers, {len(t['targets'])} spacer x target rows {dict(kinds)}; brute force "
          f"agrees on {n_checked} spacers ({n_targets} targets)")
    assert n_checked >= min(50, len(present)) and kinds["phage"] > 0 and kinds["prophage"] > 0


@needs_zombi2
def test_spacer_mix_follows_the_preset(t):
    cr = backbone.load_preset("micro")["crispr"]
    sp = t["spacers"]
    n = len(sp)
    dec = sum(r["source"] == "decoy" for r in sp)
    ph = [r for r in sp if r["source"] == "phage"]
    pam = sum(r["pam_present"] == "true" for r in ph) / len(ph)
    mm = collections.Counter(int(r["n_mismatches"]) for r in ph)
    print(f"\n{n} spacers: decoys {dec}, PAM present {pam:.3f}, mismatches {dict(sorted(mm.items()))}")
    assert abs(dec / n - cr["decoy_fraction"]) <= 4 * (cr["decoy_fraction"] * (1 - cr["decoy_fraction"]) / n) ** 0.5
    q = 1 - cr["no_pam_fraction"]
    assert abs(pam - q) <= 4 * (q * (1 - q) / len(ph)) ** 0.5
    assert max(mm) <= cr["target_max_mismatches"] and mm[0] > mm[5]


# --- runtime (§2.6) ---------------------------------------------------------------------------------------

@needs_zombi2
def test_micro_mobile_runtime_within_budget(world, record_property):
    mw, _, seconds = world
    record_property("mobile_seconds", round(seconds, 2))
    print(f"\nF-02 micro mobile stage: {seconds:.2f} s wall; {mw.meta['counts']}; digest {mw.digest}")
    # synthworld_spec §2.6: S1-S5 took about 157 s of the 3 minutes (measured, session 4), leaving about 23 s
    # for S6-S11; this stage is allowed 15 s of them
    assert seconds < 15.0
