"""Ring R2: synthworld fusions, splits and NOGD on the `micro` preset (synthworld_spec §2.1 S7, §3.3 truth_fusion
and truth_module, §5.1, §6.1, §7.3; card F-02 session 6).

A module fixture builds the micro world's S1-S6 (ZOMBI2 0.54.0), records every file the earlier stages wrote,
and builds S7 once. The oracles are written here: FASTA is parsed, reverse complements and translation are
spelled out, the lifts (S5 -> S6 by insertions.tsv, S6 -> S7 by edits.tsv) are recomputed by this file's own
formulas, clades are read from the species tree's Newick, proteins come from sequences/proteins.faa, and every
gene's protein_uid is checked against the translation of its bases. What they take from the generator are model
inputs (the preset, the S6 world), never its realized statistics.
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

import pytest
from blake3 import blake3

from biobox.synthworld import backbone, fusions, intergenic, mobile, sequences

pytestmark = [pytest.mark.r2, pytest.mark.covers("F-02")]
needs_zombi2 = pytest.mark.skipif(
    importlib.util.find_spec("zombi2") is None,
    reason="zombi2 0.54.0 is not importable here: CI installs it by `uv sync`; on the pod the image has it")

#: the S6 digest of the micro default seed as session 5 measured it (planning/measurements/
#: 2026-11-13_synthworld_mobile.md): sessions 1 to 5's world must not change
S6_DIGEST = "0bd7f0db5691de58ea5e1a0df1e6e570"
EARLIER = ("backbone", "zombi2", "sequences", "replicons", "mobile")

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
    assert len(cds) % 3 == 0 and cds[:3] in STARTS and aa[-1] == "*" and "*" not in aa[:-1], cds[:12]
    return "M" + "".join(aa[1:-1])


def uid(s: str) -> str:
    return blake3(s.encode()).hexdigest(length=16)


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


def _stage_digests(out: Path) -> dict[str, str]:
    return {d: json.loads((out / d / f"{n}.json").read_text())["digest"]
            for d, n in (("backbone", "backbone"), ("sequences", "sequences"), ("replicons", "replicons"),
                         ("mobile", "mobile"))}


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    if importlib.util.find_spec("zombi2") is None:
        pytest.skip("zombi2 0.54.0 is not importable here")
    out = tmp_path_factory.mktemp("micro_fusions") / "world"
    params = backbone.load_preset("micro")
    backbone.build(params, params["default_seed"], out)
    sequences.build(out)
    intergenic.build(out)
    mobile.build(out)
    before = {d: _tree_sha(out / d) for d in EARLIER}
    digests = _stage_digests(out)
    t0 = time.perf_counter()
    fw = fusions.build(out)
    return fw, out, time.perf_counter() - t0, before, digests


@pytest.fixture(scope="module")
def t(world):
    _, out, _, _, _ = world
    fd = out / "fusions"
    seqs, s6, s5 = {}, {}, {}
    for d, into in ((fd, seqs), (out / "mobile", s6), (out / "replicons", s5)):
        for f in sorted((d / "fasta").glob("*.fna")):
            into.update(_fasta(f))
    return {"seqs": seqs, "s6": s6, "s5": s5, "genes": _tsv(fd / "genes.tsv"), "edits": _tsv(fd / "edits.tsv"),
            "events": _tsv(fd / "events.tsv"), "fus": _tsv(fd / "fusions.tsv"), "nogd": _tsv(fd / "nogd.tsv"),
            "reps": _tsv(fd / "replicons.tsv"), "g6": _tsv(out / "mobile" / "genes.tsv"),
            "prot": _fasta(out / "sequences" / "proteins.faa")}


def _species(out: Path):
    """child -> parent of the complete species tree, from its Newick (parsed here)."""
    text = (out / "zombi2" / "species_complete.nwk").read_text(encoding="utf-8").strip()
    tokens, buf = [], ""
    for ch in text:
        if ch in "(),;":
            if buf:
                tokens.append(buf)
                buf = ""
            tokens.append(ch)
        else:
            buf += ch
    parent, stack, nodes, last, label_next = {}, [], [], None, False
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
            if label_next:
                nodes[last] = name
                k = last
            else:
                nodes.append(name)
                k = len(nodes) - 1
            parent[k] = stack[-1] if stack else None
            label_next = False
    return {nodes[k]: (nodes[p] if p is not None else None) for k, p in parent.items()}


def _below(par: dict, node: str, tips) -> set[str]:
    out = set()
    for tip in tips:
        x = tip
        while x is not None:
            if x == node:
                out.add(tip)
                break
            x = par[x]
    return out


# --- this file's lifts --------------------------------------------------------------------------------------

def _ins_lift(points, N, s, e):
    """S5 -> S6: blocks at points a (N for the origin) of lengths l."""
    s2 = s + sum(l for a, l in points if a <= s)
    if e <= N:
        return s2, e + sum(l for a, l in points if a < e)
    return s2, N + sum(l for _, l in points) + (e - N)


def _touches(edits, N, s, e):
    parts = [(s, e)] if e <= N else [(s, N + 1), (-1, e - N)]
    return any((a < b and a < y and b > x) or (a == b and x < a < y) for a, b, _ in edits for x, y in parts)


def _edit_lift(edits, N, s, e):
    """S6 -> S7: replacements [a, b) -> l bases; None when an edit touches the feature."""
    if _touches(edits, N, s, e):
        return None
    d = sum(l - (b - a) for a, b, l in edits if b <= s)
    if e <= N:
        return s + d, e + d
    return s + d, N + sum(l - (b - a) for a, b, l in edits) + (e - N)


# --- the earlier stages do not change (sessions 1 to 5) -------------------------------------------------------

@needs_zombi2
def test_earlier_stages_keep_every_file_and_their_digests(world):
    fw, out, _, before, digests = world
    for d in EARLIER:
        assert _tree_sha(out / d) == before[d], d
    assert _stage_digests(out) == digests
    assert digests["mobile"] == S6_DIGEST
    doc = fw.meta
    assert (doc["backbone_digest"], doc["sequences_digest"], doc["mobile_digest"]) == \
        (digests["backbone"], digests["sequences"], digests["mobile"])
    print(f"\nearlier stage digests {digests}")


# --- determinism (§6.1, §6.2; MR-16) -----------------------------------------------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "MR-16")
def test_same_seed_gives_byte_identical_fusions_stage_in_a_subprocess(world, tmp_path, record_property):
    _, out, _, _, _ = world
    dst = tmp_path / "copy"
    dst.mkdir()
    for d in EARLIER:
        shutil.copytree(out / d, dst / d)
    env = dict(os.environ, PYTHONHASHSEED="4242", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    code = ("import os, sys, time; from biobox.synthworld import fusions; t0 = time.perf_counter(); "
            "fusions.build(sys.argv[1]); dt = time.perf_counter() - t0; p = '/proc/self/status'; "
            "hwm = next((l.split()[1] for l in open(p) if l.startswith('VmHWM')), -1) if os.path.exists(p) else -1; "
            "print(hwm, round(dt, 2))")
    proc = subprocess.run([sys.executable, "-c", code, str(dst)], check=True, env=env, cwd=tmp_path,
                          capture_output=True, text=True)
    assert _tree_sha(dst / "fusions") == _tree_sha(out / "fusions")
    hwm_kb, secs = proc.stdout.split()[-2:]
    print(f"\nfusions stage in a child: {secs} s wall, VmHWM {int(hwm_kb) / 1024:.0f} MB")
    record_property("fusions_child_seconds", float(secs))
    if int(hwm_kb) > 0:
        record_property("fusions_child_maxrss_mb", round(int(hwm_kb) / 1024, 1))
        assert int(hwm_kb) / 1024 < 1024      # one host genome at a time (AGENTS §4.1)


# --- events sit in exactly their clades --------------------------------------------------------------------

@needs_zombi2
def test_events_follow_the_preset_and_sit_in_exactly_their_clade(world, t):
    _, out, _, _, _ = world
    cfg = backbone.load_preset("micro")["fusions"]
    par = _species(out)
    node_gid = {r["species_node"]: r["genome_id"] for r in _tsv(out / "backbone" / "genomes.tsv")}
    kinds = collections.Counter(r["kind"] for r in t["events"])
    assert kinds == {"fusion": cfg["fusions"], "split": cfg["splits"], "nogd": cfg["nogd"]}
    fams = []
    inst = collections.defaultdict(set)
    for r in t["fus"]:
        inst[r["fusion_id"]].add(r["genome_id"])
    for r in t["nogd"]:
        inst[r["event_id"]].add(r["genome_id"])
    for r in t["events"]:
        want = {node_gid[x] for x in _below(par, r["clade_node"], node_gid)}
        assert set(r["genomes"].split(";")) == want == inst[r["event_id"]], r["event_id"]
        assert cfg["clade_size"][0] <= len(want) <= cfg["clade_size"][1]
        fams += [r["family_a"]] + ([r["family_b"]] if r["kind"] == "fusion" else [])
    assert len(fams) == len(set(fams))             # no family in two events
    print(f"\nevents {dict(kinds)}; instances {sum(len(v) for v in inst.values())}; "
          f"clade sizes {sorted(len(v) for v in inst.values())}")


# --- fused ORFs, split halves and NOGD genes read as designed ---------------------------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "MR-11")
def test_fused_orfs_encode_the_two_proteins_joined_by_the_linker(t):
    ev = {r["event_id"]: r for r in t["events"]}
    g6 = {g["copy_label"]: g for g in t["g6"]}
    genes = {g["copy_label"]: g for g in t["genes"]}
    by_rec = collections.defaultdict(list)
    for g in t["g6"]:
        by_rec[(g["genome_id"], g["replicon"])].append(g)
    n = 0
    for r in t["fus"]:
        if r["kind"] != "fusion":
            continue
        (label,) = r["gene_ids"].split(";")
        la, lb = label.split("+")
        a, b = g6[la], g6[lb]
        assert (a["family_id"], b["family_id"]) == (r["family_a"], r["family_b"])
        # in S6 a and b were neighbours, co-directional, a upstream
        rows = by_rec[(a["genome_id"], a["replicon"])]
        i, j = rows.index(a), rows.index(b)
        st = int(a["strand"])
        assert int(b["strand"]) == st == int(r["strand"]) and j - i == st
        seq = t["seqs"][r["record"]]
        got = translate(_cut(seq, int(r["start"]), int(r["end"]), st))
        linker = ev[r["fusion_id"]]["linker"]
        lp = translate("ATG" + linker + "TAA")[1:]
        want = t["prot"][la] + lp + t["prot"][lb][1:]
        assert got == want
        assert int(r["junction_aa"]) == len(t["prot"][la]) and int(r["linker_aa"]) == len(lp)
        assert genes[label]["protein_uid"] == uid(got) and r["contig_uid"] == uid(seq)
        n += 1
    assert n == sum(1 for r in t["events"] if r["kind"] == "fusion" for _ in r["genomes"].split(";"))


@needs_zombi2
@pytest.mark.covers("F-02", "MR-11")
def test_split_halves_encode_the_protein_cut_between_two_domains(world, t):
    _, out, _, _, _ = world
    arch = {r["copy_label"]: r["domain_arch"] for r in _tsv(out / "sequences" / "genes.tsv")}
    ev = {r["event_id"]: r for r in t["events"]}
    genes = {g["copy_label"]: g for g in t["genes"]}
    n = 0
    for r in t["fus"]:
        if r["kind"] != "split":
            continue
        h1, h2 = (genes[x] for x in r["gene_ids"].split(";"))
        parent = h1["parents"]
        prot = t["prot"][parent]
        k = int(ev[r["fusion_id"]]["split_domain_index"])
        assert k == int(r["split_domain_index"]) >= 1
        doms = [x.rsplit(":", 1)[1].split("-") for x in arch[parent].split(";")]
        p = int(r["junction_aa"])
        assert int(doms[k - 1][1]) <= p <= int(doms[k][0])              # between domain k - 1 and domain k
        seq = t["seqs"][r["record"]]
        st = int(r["strand"])
        p1 = translate(_cut(seq, int(h1["start"]), int(h1["end"]), st))
        p2 = translate(_cut(seq, int(h2["start"]), int(h2["end"]), st))
        assert p1 == prot[:p] and p2 == "M" + prot[p:]
        assert (h1["protein_uid"], h2["protein_uid"]) == (uid(p1), uid(p2))
        # the first half is upstream, the spacer between them holds an RBS
        gap = (int(h1["end"]), int(h2["start"])) if st == 1 else (int(h2["end"]), int(h1["start"]))
        spacer = _cut(seq, *gap, st)
        assert spacer == ev[r["fusion_id"]]["spacer"] and "AGGAGG" in spacer
        assert (int(r["start"]), int(r["end"])) == (min(int(h1["start"]), int(h2["start"])),
                                                    max(int(h1["end"]), int(h2["end"])))
        n += 1
    assert n == sum(len(r["genomes"].split(";")) for r in t["events"] if r["kind"] == "split")


@needs_zombi2
def test_nogd_puts_an_unrelated_family_in_place_of_a_module_member(world, t):
    _, out, _, _, _ = world
    bgenes = {r["copy_label"]: r for r in _tsv(out / "backbone" / "genes.tsv")}
    groups = {r["name"]: r for r in _tsv(out / "backbone" / "groups.tsv")}
    seq_fams = {r["family_id"]: r for r in _tsv(out / "sequences" / "families.tsv")}
    nfam = {r["family_id"]: r for r in _tsv(out / "fusions" / "nogd_families.tsv")}
    g6 = {g["copy_label"]: g for g in t["g6"]}
    by_rec7 = collections.defaultdict(list)
    for g in t["genes"]:
        by_rec7[(g["genome_id"], g["replicon"])].append(g)
    by_rec6 = collections.defaultdict(list)
    for g in t["g6"]:
        by_rec6[(g["genome_id"], g["replicon"])].append(g)
    prots = collections.defaultdict(list)
    for r in t["nogd"]:
        old = g6[r["replaced_copy"]]
        assert old["family_id"] == r["family_replaced"]
        b = bgenes[r["replaced_copy"]]
        assert b["group"] == r["module_id"] and groups[b["group"]]["kind"] == "planted"
        seq = t["seqs"][r["record"]]
        p = translate(_cut(seq, int(r["start"]), int(r["end"]), int(r["strand"])))
        assert uid(p) == r["protein_uid"] and int(r["strand"]) == int(old["strand"])
        prots[r["event_id"]].append(p)
        # the same neighbours on both sides as in S6, the replaced gene's place taken
        key = (old["genome_id"], old["replicon"])
        names6 = [g["copy_label"] for g in by_rec6[key]]
        names7 = [g["copy_label"] for g in by_rec7[key]]
        i = names6.index(r["replaced_copy"])
        j = names7.index(r["copy_label"])
        assert names7[j - 1] == names6[i - 1] and names7[j + 1] == names6[i + 1]
        # unrelated: no domain in common with the family it displaces
        assert not set(nfam[r["family_replacement"]]["domain_hmms"].split(";")) \
            & set(seq_fams[r["family_replaced"]]["domain_hmms"].split(";"))
    assert len(prots) == backbone.load_preset("micro")["fusions"]["nogd"]


# --- every gene and every earlier feature maps into the new replicons ------------------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "B-SIM")
def test_every_gene_translates_to_its_protein_uid(t):
    n = collections.Counter()
    for g in t["genes"]:
        seq = t["seqs"][f"{g['genome_id']}_r{g['replicon']}"]
        s, e = int(g["start"]), int(g["end"])
        assert 0 <= s < len(seq) and s < e < s + len(seq)
        assert (g["origin_spanning"] == "true") == (e > len(seq))
        assert uid(translate(_cut(seq, s, e, int(g["strand"])))) == g["protein_uid"], g["copy_label"]
        n[g["source"]] += 1
    print(f"\ngenes by source {dict(n)}")
    assert sum(g["origin_spanning"] == "true" for g in t["genes"]) >= 5      # §5.1


@needs_zombi2
def test_every_earlier_feature_reads_the_same_bases_or_belongs_to_its_event(world, t):
    """S6 genes, insertions, prophages, arrays and spacers, and S5 elements, motifs and IGRs (through S6's
    insertions), lifted by edits.tsv, read the same bases in S7; a feature an edit touches is one of that event's
    genes or lies between a fused pair."""
    _, out, _, _, _ = world
    ed = collections.defaultdict(list)
    ev_of = collections.defaultdict(list)
    for r in t["edits"]:
        ed[r["record"]].append((int(r["start"]), int(r["end"]), int(r["length"])))
        ev_of[r["record"]].append(r)
    for rec, x in ed.items():
        x.sort()
    ins = collections.defaultdict(list)
    for r in _tsv(out / "mobile" / "insertions.tsv"):
        ins[r["record"]].append((int(r["at"]), int(r["length"])))
    reps = {r["record"]: r for r in t["reps"]}
    # edited spans per record: parents of the event genes, and the fused pairs' copies
    parents = {}
    for g in t["genes"]:
        for p in filter(None, g["parents"].split(";")):
            parents[p] = g["event_id"]
    fused_pairs = {tuple(g["parents"].split(";")) for g in t["genes"] if g["source"] == "fused"}
    fused_juncs = {(a, b) for a, b in fused_pairs} | {(b, a) for a, b in fused_pairs}

    def check(rec, s, e, owner=None):
        N6 = len(t["s6"][rec])
        got = _edit_lift(ed.get(rec, []), N6, s, e)
        if got is None:
            assert owner is not None, (rec, s, e)
            return 0
        assert _cut(t["seqs"][rec], *got) == _cut(t["s6"][rec], s, e)
        return 1

    n = collections.Counter()
    g7 = {g["copy_label"]: g for g in t["genes"]}
    for g in t["g6"]:
        rec = f"{g['genome_id']}_r{g['replicon']}"
        if g["copy_label"] in parents:
            assert g["copy_label"] not in g7
            n["event_parent"] += 1
            continue
        n["gene"] += check(rec, int(g["start"]), int(g["end"]))
        assert (int(g7[g["copy_label"]]["start"]), int(g7[g["copy_label"]]["end"])) == \
            _edit_lift(ed.get(rec, []), len(t["s6"][rec]), int(g["start"]), int(g["end"]))
    for name, cols in (("insertions.tsv", ("start", "end")), ("prophages.tsv", ("start", "end")),
                       ("arrays.tsv", ("locus_start", "locus_end")), ("spacer_targets.tsv", ("spacer_start", "spacer_end"))):
        for r in _tsv(out / "mobile" / name):
            n[name] += check(r["record"], int(r[cols[0]]), int(r[cols[1]]))
    # S5 features through both lifts
    s5_genes = {g["copy_label"]: g for g in _tsv(out / "replicons" / "genes.tsv")}
    owner_by_pos = {(g["genome_id"], g["replicon"], g["position"]): g["copy_label"] for g in s5_genes.values()}
    for name in ("elements.tsv", "rna.tsv", "junctions.tsv"):
        for r in _tsv(out / "replicons" / name):
            rec = f"{r['genome_id']}_r{r['replicon']}"
            s, e = int(r["start"]), int(r["end"])
            if e == s:
                continue
            s6, e6 = _ins_lift(ins[rec], len(t["s5"][rec]), s, e)
            owner = None
            if name == "junctions.tsv":
                pair = (r["left_copy"], r["right_copy"])
                owner = pair if pair in fused_juncs else None
            elif name == "elements.tsv":
                lc = owner_by_pos.get((r["genome_id"], r["replicon"], r["left_position"]))
                owner = r["copy_label"] if r["copy_label"] in parents or lc in parents else None
            else:
                owner = r["attached_copy"] if r["attached_copy"] in parents else None
            got = check(rec, s6, e6, owner)
            n[name] += got
            n[name + "_in_event"] += 1 - got
    for rec, seq in t["seqs"].items():
        assert int(reps[rec]["s6_length"]) == len(t["s6"][rec])
        assert len(seq) == len(t["s6"][rec]) + sum(l - (b - a) for a, b, l in ed.get(rec, []))
        assert reps[rec]["contig_uid"] == uid(seq)
    print(f"\nfeatures read the same or belong to an event: {dict(n)}")
    assert n["junctions.tsv_in_event"] > 0 and n["elements.tsv_in_event"] > 0


# --- runtime (§2.6) ---------------------------------------------------------------------------------------

@needs_zombi2
def test_micro_fusions_runtime_within_budget(world, record_property):
    fw, _, seconds, _, _ = world
    record_property("fusions_seconds", round(seconds, 2))
    print(f"\nF-02 micro fusions stage: {seconds:.2f} s wall; {fw.meta['counts']}; digest {fw.digest}")
    # synthworld_spec §2.6: S1-S6 take about 168 s of the 3 minutes (measured, sessions 4 and 5), leaving about
    # 12 s for S7-S11 (derived); this bound is S7's share, measured at session 6
    assert seconds < 6.0
