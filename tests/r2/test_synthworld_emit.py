"""Ring R2: synthworld's emit stage on the `micro` preset (synthworld_spec §2.1 S11, §3, §6.2, §7.3; card F-02
session 8).

A module fixture builds the micro world's S1-S10 (ZOMBI2 0.54.0), timing each stage, records every file the earlier
stages wrote, and builds S11 once. The oracles are written here: FASTA is parsed and checked at 80 columns, reverse
complements and the genetic codes are spelled out, uids are hashed with blake3 in this file (with their own byte
packing for feature_uid and region_uid), the S8 frame is rebuilt from the S9 contigs by slicing, a feature's expected
bases are read from its own stage's records (S5 replicons, S6 and S7 host records, the S8 frame), and a contig's
orientation comes from the S9 and S10 tables, never from the emitted truth.
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

from biobox.synthworld import assembly, backbone, emit, fusions, intergenic, mobile, sequences, sources

pytestmark = [pytest.mark.r2, pytest.mark.covers("F-02")]
needs_zombi2 = pytest.mark.skipif(
    importlib.util.find_spec("zombi2") is None,
    reason="zombi2 0.54.0 is not importable here: CI installs it by `uv sync`; on the pod the image has it")

#: the stage digests of sessions 1 to 7 for the micro default seed (measured 2026-11-13): emit must not change them
DIGESTS = {"backbone": "025fdd0d27fdc65e678782216ecd21fb", "sequences": "34ebce1324413a261acd874d7ee4fe3e",
           "replicons": "55bda97155d75148e5cdca7b90a1d7be", "mobile": "0bd7f0db5691de58ea5e1a0df1e6e570",
           "fusions": "db09abf3fe79ae34a73e3ab0664a9b13", "assembly": "2f425158dc38f2302c1d02e9ebed30da",
           "sources": "770007bde8615f41d9dd5f2eb912e753"}
EARLIER = ("backbone", "zombi2", "sequences", "replicons", "mobile", "fusions", "assembly", "sources")
#: the S11 digest and the world digest of the micro default seed, measured at this session's first complete emit
#: (wilkcomp02, 2026-11-13; the same in a separate world build in the runtime image)
EMIT_DIGEST = "a2f0e448484dab63d41281a5ef6a00c0"
WORLD_DIGEST = "f0262027e4ccdcbc72cd7d90f23fbb5c"

CODE11_BY_AA = {
    "F": "TTT TTC", "L": "TTA TTG CTT CTC CTA CTG", "I": "ATT ATC ATA", "M": "ATG", "V": "GTT GTC GTA GTG",
    "S": "TCT TCC TCA TCG AGT AGC", "P": "CCT CCC CCA CCG", "T": "ACT ACC ACA ACG", "A": "GCT GCC GCA GCG",
    "Y": "TAT TAC", "*": "TAA TAG TGA", "H": "CAT CAC", "Q": "CAA CAG", "N": "AAT AAC", "K": "AAA AAG",
    "D": "GAT GAC", "E": "GAA GAG", "C": "TGT TGC", "W": "TGG", "R": "CGT CGC CGA CGG AGA AGG",
    "G": "GGT GGC GGA GGG",
}
CODE11 = {c: aa for aa, cs in CODE11_BY_AA.items() for c in cs.split()}
CODE = {11: CODE11, 4: dict(CODE11, TGA="W"), 25: dict(CODE11, TGA="G")}
STARTS = {"ATG", "GTG", "TTG"}
_RC = {"A": "T", "C": "G", "G": "C", "T": "A"}
KIND = {"cds": 1, "ncrna": 2, "crispr_array": 4}          # biobox.ids.KIND, spelled out


def rc(s: str) -> str:
    return "".join(_RC[c] for c in reversed(s))


def uid(s: str) -> str:
    return blake3(s.encode()).hexdigest(length=16)


def feature_uid(cuid: str, s: int, e: int, st: int, kind: int) -> str:
    return blake3(bytes.fromhex(cuid) + s.to_bytes(8, "big") + e.to_bytes(8, "big") + st.to_bytes(1, "big", signed=True)
                  + kind.to_bytes(1, "big")).hexdigest(length=16)


def region_uid(cuid: str, s: int, e: int) -> str:
    return blake3(bytes.fromhex(cuid) + s.to_bytes(8, "big") + e.to_bytes(8, "big")
                  + (1).to_bytes(1, "big", signed=True)).hexdigest(length=16)


def translate(cds: str, code: int) -> str | None:
    if len(cds) % 3 or cds[:3] not in STARTS:
        return None
    aa = [CODE[code][cds[i:i + 3]] for i in range(0, len(cds), 3)]
    if aa[-1] != "*" or "*" in aa[:-1]:
        return None
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


def _cut(seq: str, s: int, e: int, st: int = 1) -> str:
    n = len(seq)
    assert 0 <= s <= e <= s + n, (s, e, n)
    x = (seq + seq)[s:e] if e > n else seq[s:e]
    return x if st == 1 else rc(x)


def _tree_sha(root: Path) -> dict[str, str]:
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file()}


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    if importlib.util.find_spec("zombi2") is None:
        pytest.skip("zombi2 0.54.0 is not importable here")
    out = tmp_path_factory.mktemp("micro_emit") / "world"
    params = backbone.load_preset("micro")
    times = {}
    for name, fn in (("S1-S2 backbone", lambda: backbone.build(params, params["default_seed"], out)),
                     ("S3-S4 sequences", lambda: sequences.build(out)),
                     ("S5 replicons", lambda: intergenic.build(out)),
                     ("S6 mobile", lambda: mobile.build(out)),
                     ("S7 fusions", lambda: fusions.build(out)),
                     ("S8-S9 assembly", lambda: assembly.build(out)),
                     ("S10 sources", lambda: sources.build(out))):
        t0 = time.perf_counter()
        fn()
        times[name] = time.perf_counter() - t0
    before = {d: _tree_sha(out / d) for d in EARLIER}
    t0 = time.perf_counter()
    ew = emit.build(out)
    times["S11 emit"] = time.perf_counter() - t0
    return ew, out, times, before


@pytest.fixture(scope="module")
def t(world):
    _, out, _, _ = world
    ed = out / "emit"
    ingest = _tsv(ed / "ingest.tsv")
    fa = {}
    for rel in sorted({r["file"] for r in ingest}):
        fa.update(_fasta(ed / rel))
    seq_of_uid = {}
    for r in ingest:
        seq_of_uid[r["contig_uid"]] = fa[r["contig_accession"]]
    # the S8 frame, rebuilt from the S9 contigs (emitted and dropped)
    ad = out / "assembly"
    contigs = _tsv(ad / "contigs.tsv")
    reps = {r["record"]: r for r in _tsv(ad / "replicons.tsv")}
    cseq = {}
    for d in ("contigs", "dropped"):
        for f in sorted((ad / d).glob("*.fna")):
            cseq.update(_fasta(f))
    buf = {rec: bytearray(b"-" * int(r["length"])) for rec, r in reps.items()}
    for c in contigs:
        s = cseq[c["contig_id"]]
        fwd = (s if c["orientation"] == "1" else rc(s)).encode()
        a, b, N = int(c["replicon_start"]), int(c["replicon_end"]), len(buf[c["record"]])
        if b <= N:
            buf[c["record"]][a:b] = fwd
        else:
            buf[c["record"]][a:N] = fwd[:N - a]
            buf[c["record"]][0:b - N] = fwd[N - a:]
    frames = {rec: bytes(v).decode() for rec, v in buf.items()}
    # each emitted contig's orientation against the frame, from the S9 and S10 tables
    by_id = {c["contig_id"]: c for c in contigs}
    orient = {c["contig_uid"]: int(c["orientation"]) for c in contigs if c["emitted"] == "true"}
    dup = _tsv(out / "sources" / "dup.tsv")
    contam = _tsv(out / "sources" / "contam.tsv")
    for d in dup:
        orient[d["member_contig_uid"]] = int(by_id[d["rep_contig_id"]]["orientation"])
    for c in contam:
        orient[c["contig_uid"]] = int(by_id[c["donor_contig_id"]]["orientation"]) * int(c["orientation"])
    muts = collections.defaultdict(set)
    member_uid = {d["member_contig_id"]: d["member_contig_uid"] for d in dup}
    for m in _tsv(out / "sources" / "mutations.tsv"):
        muts[member_uid[m["member_contig_id"]]].add(int(m["position"]))
    truth = {name: _tsv(ed / "truth" / f"{name}.tsv") for name in emit.TABLES}
    return {"ingest": ingest, "fa": fa, "seq": seq_of_uid, "frames": frames, "contigs": contigs, "orient": orient,
            "muts": muts, "dup": dup, "contam": contam, "truth": truth, "genes8": _tsv(ad / "genes.tsv"),
            "pieces": _tsv(ad / "gene_contigs.tsv"), "genomes": _tsv(ad / "genomes.tsv")}


def _matches(t, cu: str, s: int, e: int, st: int, want: str) -> bool:
    """The emitted bases [s, e) of contig ``cu`` read on strand ``st`` equal ``want``, or (a near duplicate's copy)
    differ from it only at the copy's listed mutation positions."""
    seq = t["seq"][cu]
    got = _cut(seq, s, e, 1)
    w = want if st == 1 else rc(want)
    if got == w:
        return True
    if len(got) != len(w) or cu not in t["muts"]:
        return False
    return all((s + i) % len(seq) in t["muts"][cu] for i, (x, y) in enumerate(zip(got, w)) if x != y)


# --- the earlier stages do not change ------------------------------------------------------------------------------

@needs_zombi2
def test_earlier_stages_keep_every_file_and_their_digests(world):
    ew, out, _, before = world
    for d in EARLIER:
        assert _tree_sha(out / d) == before[d], d
    got = {d: json.loads((out / d / f"{d}.json").read_text())["digest"] for d in DIGESTS}
    assert got == DIGESTS
    assert {d: ew.meta[f"{d}_digest"] for d in DIGESTS} == DIGESTS
    assert (ew.digest, ew.world_digest) == (EMIT_DIGEST, WORLD_DIGEST)
    print(f"\nemit digest {ew.digest}; world digest {ew.world_digest}; release_tag {ew.release_tag}")


@needs_zombi2
def test_manifest_digests_are_recomputable(world):
    ew, out, _, _ = world
    man = json.loads((out / "manifest.json").read_text())
    stages = {d: v["digest"] for d, v in man["stages"].items()}
    assert {d: stages[d] for d in DIGESTS} == DIGESTS and stages["emit"] == ew.digest
    lines = "".join(f"{k}\t{v}\n" for k, v in sorted(stages.items()))
    assert man["world_digest"] == blake3(lines.encode()).hexdigest(length=16)
    tag = "".join(f"{k}\t{v}\n" for k, v in sorted(DIGESTS.items()))
    assert man["release_tag"] == "sw-" + blake3(tag.encode()).hexdigest(length=8)
    doc = json.loads((out / "emit" / "emit.json").read_text())
    files = _tree_sha(out / "emit")
    files.pop("emit.json")
    assert files == doc["files"]
    lines = "".join(f"{k}\t{v}\n" for k, v in sorted(files.items()))
    assert doc["digest"] == blake3(lines.encode()).hexdigest(length=16)
    assert man["versions"]["zombi2"] == "0.54.0"


# --- determinism (§6.1, §6.2; MR-16) ---------------------------------------------------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "MR-16")
def test_same_seed_gives_byte_identical_emit_stage_in_a_subprocess(world, tmp_path, record_property):
    _, out, _, _ = world
    dst = tmp_path / "copy"
    dst.mkdir()
    for d in EARLIER:
        shutil.copytree(out / d, dst / d)
    env = dict(os.environ, PYTHONHASHSEED="919", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    code = ("import os, sys, time; from biobox.synthworld import emit; t0 = time.perf_counter(); "
            "emit.build(sys.argv[1]); dt = time.perf_counter() - t0; p = '/proc/self/status'; "
            "hwm = next((l.split()[1] for l in open(p) if l.startswith('VmHWM')), -1) if os.path.exists(p) else -1; "
            "print(hwm, round(dt, 2))")
    proc = subprocess.run([sys.executable, "-c", code, str(dst)], check=True, env=env, cwd=tmp_path,
                          capture_output=True, text=True)
    assert _tree_sha(dst / "emit") == _tree_sha(out / "emit")
    assert (dst / "manifest.json").read_bytes() == (out / "manifest.json").read_bytes()
    hwm_kb, secs = proc.stdout.split()[-2:]
    print(f"\nemit stage in a child: {secs} s wall, VmHWM {int(hwm_kb) / 1024:.0f} MB")
    record_property("emit_child_seconds", float(secs))
    if int(hwm_kb) > 0:
        record_property("emit_child_maxrss_mb", round(int(hwm_kb) / 1024, 1))
        assert int(hwm_kb) / 1024 < 1024


# --- sources, the ingest manifest and the build-stage inputs --------------------------------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "MR-06")
def test_every_source_file_holds_its_contigs_as_the_earlier_stages_wrote_them(world, t):
    _, out, _, _ = world
    occ = {r["contig_accession"]: r for r in _tsv(out / "sources" / "source_contigs.tsv")}
    phages = {r["phage_id"]: r for r in _tsv(out / "mobile" / "phages.tsv")}
    ph_seq = _fasta(out / "mobile" / "phages.fna")
    srcs = collections.defaultdict(list)
    for r in t["ingest"]:
        srcs[r["file"]].append(r)
        s = t["fa"][r["contig_accession"]]
        assert len(s) == int(r["length"]) and uid(s) == r["contig_uid"]
        if r["contig_accession"] not in occ:
            p = phages[r["genome_id"]]
            assert s == ph_seq[p["record"]] and r["kind"] == "virus" and r["replicon_class"] == "virus"
    # every host contig occurrence of S10 is emitted once, under its accession, in its source's file
    for acc, o in occ.items():
        assert acc in t["fa"]
    files_seen = {}
    for o in occ.values():
        if o["file"] not in files_seen:
            files_seen[o["file"]] = _fasta(out / o["file"])
        assert t["fa"][o["contig_accession"]] == files_seen[o["file"]][o["record"]]
    # one file per source, records in source_contigs order
    for f, rows in srcs.items():
        assert list(_fasta(out / "emit" / f)) == [r["contig_accession"] for r in rows]
        assert len({(r["source_label"], r["genome_id"]) for r in rows}) == 1
    # exact duplicates share contig_uid; near duplicates do not
    by_uid = collections.defaultdict(set)
    for r in t["ingest"]:
        by_uid[r["contig_uid"]].add(r["source_label"])
    for d in t["dup"]:
        # a near duplicate's contig with no site to change (S10: wholly inside a gene cut at both edges) keeps it
        same = d["kind"] == "exact" or int(d["n_mutations"]) == 0
        assert (d["rep_contig_uid"] == d["member_contig_uid"]) == same
        if d["kind"] == "exact":
            assert {"src_a", d["member_source_label"]} <= by_uid[d["rep_contig_uid"]]
    gen = {g["genome_id"]: g for g in t["genomes"]}
    for r in t["ingest"]:
        if r["kind"] == "prokaryote":
            g = gen[r["genome_id"]]
            assert int(r["genetic_code"]) == int(g["genetic_code"])
            assert r["source_class"] == {"mag": "MAG"}.get(g["assembly_class"], "isolate")
            assert r["molecule_type"] == "DNA"
    print(f"\nsources {len(srcs)}; contig occurrences {len(t['ingest'])}; distinct contigs {len(t['seq'])}")


@needs_zombi2
@pytest.mark.covers("F-02", "MR-06")
def test_build_stage_inputs_are_accepted_by_the_builders_readers(world, t):
    """The builders' own readers are the consumers here (the `--genomes`/`--objects` parsers and genome_input's
    refusals, DEC-27): the emitted files must pass them."""
    from biobox import cli
    from biobox.l0seq import inputs

    _, out, _, _ = world
    ed = out / "emit"
    genomes = cli._genomes(str(ed / "genomes.json"))
    objects = cli._source_objects(str(ed / "objects.json"))
    rows = inputs._rows(genomes)
    assert {r[0].hex() for r in rows} == set(t["seq"])
    assert all(o.path.exists() for o in objects) and len(objects) == len({r["file"] for r in t["ingest"]})
    n_exact = sum(1 for d in t["dup"] if d["kind"] == "exact")
    n_src = len({(r["source_label"], r["genome_id"]) for r in t["ingest"]})
    exact_sources = len({(d["member_source_label"], d["member_genome_id"]) for d in t["dup"] if d["kind"] == "exact"})
    assert len(genomes) == n_src - exact_sources and n_exact > 0
    circ = {r["contig_uid"] for r in t["ingest"] if r["topology"] == "circular"}
    assert {c.hex() for g in genomes for c in g.circular} == circ
    complete = {c["contig_uid"] for c in t["contigs"] if c["is_circular"] == "true" and c["emitted"] == "true"}
    assert complete <= circ
    # labels: a mislabelled source names another genus; a biome is shared by a class
    tg = t["truth"]["truth_genome"]
    lab = collections.Counter(r["mislabeled"] for r in tg if r["is_phage"] == "false")
    for r in tg:
        if r["mislabeled"] == "true":
            g = lambda x: next(y for y in x.split(";") if y.startswith("g__"))   # noqa: E731
            assert g(r["tax_lineage"]) != g(r["tax_lineage_true"])
    biome = collections.defaultdict(set)
    for r in tg:
        if r["is_phage"] == "false":
            biome[next(y for y in r["tax_lineage_true"].split(";") if y.startswith("c__"))].add(r["biome"])
    assert all(len(v) == 1 for v in biome.values())
    print(f"\ngenomes listed {len(genomes)} ({n_src} sources, {exact_sources} exact duplicates listed once); "
          f"mislabelled {lab['true']} of {sum(lab.values())}; biomes by class {len(biome)}")


# --- truth in emitted-contig coordinates ----------------------------------------------------------------------------

@needs_zombi2
@pytest.mark.covers("F-02", "B-SIM", "MR-05")
def test_every_truth_gene_reads_its_s8_gene_on_the_emitted_contig(t):
    g8 = {g["copy_label"]: g for g in t["genes8"]}
    frameshift = {r["true_gene_id"] for r in _tsv_callerr(t) if r["kind"].startswith("frameshift")}
    code = {g["genome_id"]: int(g["genetic_code"]) for g in t["genomes"]}
    n = collections.Counter()
    keys = set()
    for r in t["truth"]["truth_gene"]:
        s, e, st, cu = int(r["start"]), int(r["end"]), int(r["strand"]), r["contig_uid"]
        assert r["feature_uid"] == feature_uid(cu, s, e, st, KIND["cds"])
        assert r["feature_uid"] not in keys
        keys.add(r["feature_uid"])
        seq = t["seq"][cu]
        got = _cut(seq, s, e, st)
        assert (r["origin_spanning"] == "true") == (e > len(seq))
        want_stop = "" if r["partial_3"] == "true" else str(e if st == 1 else s)
        assert r["stop_bp"] == want_stop
        if r["occurrence"] == "phage":
            assert translate(got, 11) is not None and uid(translate(got, 11)) == r["protein_uid"]
            n["phage"] += 1
            continue
        g = g8[r["gene_id"]]
        whole = _cut(t["frames"][f"{g['genome_id']}_r{g['replicon']}"], int(g["start"]), int(g["end"]),
                     int(g["strand"]))
        status = "cut" if r["partial_5"] == "true" or r["partial_3"] == "true" else "whole"
        assert r["status"] == status
        if status == "whole":
            assert _matches(t, cu, s, e, st, whole), r["gene_id"]
            if r["occurrence"] == "near_dup":     # its changes are synonymous (S10)
                tab = CODE[code[g["genome_id"]]]
                assert [tab[got[i:i + 3]] for i in range(0, len(got) - 2, 3)] == \
                    [tab[whole[i:i + 3]] for i in range(0, len(whole) - 2, 3)]
            else:
                assert got == whole
            if r["feature_kind"] in ("CDS", "smORF") and r["gene_id"] not in frameshift:
                assert uid(translate(got, code[g["genome_id"]])) == r["protein_uid"]
        else:
            assert len(got) < len(whole)
            if r["occurrence"] != "near_dup":
                assert got in whole
            p5, p3 = r["partial_5"] == "true", r["partial_3"] == "true"
            if r["occurrence"] != "near_dup" and not (p5 and p3):
                assert whole.startswith(got) if p3 else whole.endswith(got)
        n[f"{r['occurrence']}_{status}"] += 1
    # every emitted piece of S9 is a row of the assembly occurrence
    em = sum(1 for p in t["pieces"] if p["emitted"] == "true")
    assert n["assembly_whole"] + n["assembly_cut"] == em
    print(f"\ntruth_gene rows by occurrence and status {dict(n)}")


def _tsv_callerr(t):
    return t["truth"]["truth_callerr"]


@needs_zombi2
@pytest.mark.covers("F-02", "F-08")
def test_call_errors_sit_on_their_gene_on_the_emitted_contig(t):
    g8 = {g["copy_label"]: g for g in t["genes8"]}
    n = collections.Counter()
    for r in t["truth"]["truth_callerr"]:
        n[r["status"]] += 1
        if r["status"] != "whole":
            assert r["contig_uid"] == "" and r["status"] in ("partial", "dropped")
            continue
        s, e, st, cu = int(r["true_start"]), int(r["true_end"]), int(r["strand"]), r["contig_uid"]
        g = g8[r["true_gene_id"]]
        whole = _cut(t["frames"][f"{g['genome_id']}_r{g['replicon']}"], int(g["start"]), int(g["end"]),
                     int(g["strand"]))
        assert _matches(t, cu, s, e, st, whole), r["callerr_id"]
        assert r["feature_uid"] == feature_uid(cu, s, e, st, KIND["cds"])
        L = len(t["seq"][cu])               # S9 may rotate a complete circular contig's origin into the gene
        if r["kind"] != "alt_start_upstream":
            assert r["error_pos"] != "" and (int(r["error_pos"]) - s) % L < e - s
        else:
            a = int(r["alt_start"])
            up = _cut(t["seq"][cu], a, s) if st == 1 else _cut(t["seq"][cu], e, a + 3, -1)
            assert up[:3] == "ATG" and len(up) % 3 == 0
        for item in filter(None, r["expected_calls"].split(";")):
            span, cst = item.split(":")
            x, y = (int(v) for v in span.split("-"))
            ox = (x - s) % L
            # a call lies inside the gene, on its strand, and keeps its true start or its true stop
            assert int(cst) == st and ox + (y - x) <= e - s and (ox == 0 or ox + (y - x) == e - s)
    print(f"\ntruth_callerr by status {dict(n)}")
    assert n["whole"] > 0


@needs_zombi2
def test_s5_s6_s7_features_read_their_own_stage_bases_on_the_emitted_contigs(world, t):
    """IGRs and RNA motifs (S5 records), spacers and protospacers (S6 spacer table), fused ORFs (S7 records):
    every whole row's bases equal the feature's own bases (a near duplicate may differ only at its listed
    mutations); the rest say why not."""
    _, out, _, _ = world
    s5 = {k: v for f in sorted((out / "replicons" / "fasta").glob("*.fna")) for k, v in _fasta(f).items()}
    s7 = {k: v for f in sorted((out / "fusions" / "fasta").glob("*.fna")) for k, v in _fasta(f).items()}
    n = collections.Counter()
    for r in t["truth"]["truth_igr"]:
        n[f"igr_{r['status']}"] += 1
        if r["status"] not in ("whole", "cut"):
            assert r["status"] in ("dropped", "s6_insertion", "s7_event", "s8_edit") and r["contig_uid"] == ""
            continue
        s, e = int(r["start"]), int(r["end"])
        if s == e:
            assert r["igr_uid"] == ""
            continue
        rec = f"{r['genome_id']}_r{r['replicon']}"
        own = _cut(s5[rec], int(r["s5_start"]), int(r["s5_end"]))
        cu = r["contig_uid"]
        got = _cut(t["seq"][cu], s, e)
        assert r["igr_uid"] == uid(got) and r["region_uid"] == region_uid(cu, s, e)
        if r["status"] == "whole":
            assert _matches(t, cu, s, e, t["orient"][cu], own), (rec, r["s5_start"], cu, s, e)
        else:
            assert len(got) < len(own)
    rna = {r["rna_id"]: r for r in _tsv(out / "replicons" / "rna.tsv")}
    for r in t["truth"]["truth_rna"]:
        n[f"rna_{r['status']}"] += 1
        if r["status"] != "whole":
            continue
        x = rna[r["rna_id"]]
        own = _cut(s5[f"{x['genome_id']}_r{x['replicon']}"], int(x["start"]), int(x["end"]), int(x["strand"]))
        cu, s, e, st = r["contig_uid"], int(r["start"]), int(r["end"]), int(r["strand"])
        assert _matches(t, cu, s, e, st, own), r["rna_id"]
        assert r["feature_uid"] == feature_uid(cu, s, e, st, KIND["ncrna"])
    for r in t["truth"]["truth_spacer_target"]:
        n[f"spacer_{r['status']}_{r['target_status']}"] += 1
        if r["status"] == "whole":
            cu = r["contig_uid"]
            assert _matches(t, cu, int(r["spacer_start"]), int(r["spacer_end"]), int(r["array_strand"]),
                            r["spacer_seq"]), r["spacer_id"]
            assert r["array_feature_uid"] == feature_uid(cu, int(r["array_start"]), int(r["array_end"]),
                                                         int(r["array_strand"]), KIND["crispr_array"])
        if r["target_status"] == "whole" and r["target_contig_uid"]:
            tu, ps, pe, pst = r["target_contig_uid"], int(r["proto_start"]), int(r["proto_end"]), int(r["proto_strand"])
            pro = _cut(t["seq"][tu], ps, pe, 1)
            sp = r["spacer_seq"] if pst == 1 else rc(r["spacer_seq"])
            diff = sum(1 for x, y in zip(pro, sp) if x != y)
            m_w = sum(1 for i in range(len(pro)) if (ps + i) % len(t["seq"][tu]) in t["muts"].get(tu, ()))
            assert len(pro) == len(sp)
            # exactly n_mismatches, but on a near duplicate's copy each S10 mutation inside may add or remove one
            assert abs(diff - int(r["n_mismatches"])) <= m_w, (r["spacer_id"], tu)
            n["protospacer_with_s10_mutations"] += m_w > 0
    fus = {(x["fusion_id"], x["genome_id"]): x for x in _tsv(out / "fusions" / "fusions.tsv")}
    for r in t["truth"]["truth_fusion"]:
        n[f"fusion_{r['status']}"] += 1
        if r["status"] != "whole":
            continue
        src = fus[(r["fusion_id"], r["genome_id"])]
        own = _cut(s7[src["record"]], int(src["start"]), int(src["end"]), int(src["strand"]))
        assert _matches(t, r["contig_uid"], int(r["start"]), int(r["end"]), int(r["strand"]), own), r["fusion_id"]
    print(f"\nS5-S7 features on the emitted contigs {dict(sorted(n.items()))}")
    assert n["igr_whole"] > 50_000 and n["rna_whole"] > 500


@needs_zombi2
@pytest.mark.covers("F-02", "MR-04", "MR-05")
def test_contig_and_edge_tables_cover_every_emitted_contig(world, t):
    _, out, _, _ = world
    tc = t["truth"]["truth_contig"]
    emitted = [r for r in tc if r["is_dropped"] == "false"]
    assert sorted(r["accession"] for r in emitted) == sorted(r["contig_accession"] for r in t["ingest"])
    for r in emitted:
        assert uid(t["fa"][r["accession"]]) == r["contig_uid"] and int(r["length"]) == len(t["fa"][r["accession"]])
    dropped = {r["contig_uid"] for r in tc if r["is_dropped"] == "true"}
    assert dropped == {c["contig_uid"] for c in t["contigs"] if c["emitted"] == "false"}
    edges = t["truth"]["truth_edges"]
    assert sorted(r["contig_uid"] for r in edges) == sorted(t["seq"])
    for r in edges:
        assert int(r["replicon_end"]) - int(r["replicon_start"]) == len(t["seq"][r["contig_uid"]])
    circ = [r for r in emitted if r["is_circular"] == "true"]
    assert len(circ) >= 5 and all(r["origin"] in ("assembly", "exact_dup", "near_dup") for r in circ)
    assert sum(r["is_contaminant"] == "true" for r in emitted) == len(t["contam"])
    print(f"\ncontigs emitted {len(emitted)} ({len(t['seq'])} distinct), dropped {len(dropped)}, circular "
          f"{len(circ)}, prophage-bearing {sum(r['prophage_bearing'] == 'true' for r in emitted)}")


@needs_zombi2
def test_gff_lists_every_gene_on_every_accession(world, t):
    _, out, _, _ = world
    accs = collections.Counter(r["contig_uid"] for r in t["ingest"])
    want = sum(accs[r["contig_uid"]] for r in t["truth"]["truth_gene"])
    lines = (out / "emit" / "truth.gff").read_text().splitlines()
    assert lines[0] == "##gff-version 3"
    got = sum(1 for x in lines if not x.startswith("#") and x.split("\t")[2] in ("CDS", "pseudogene"))
    assert got == want


@needs_zombi2
def test_family_module_and_trap_tables_are_complete(world, t):
    _, out, _, _ = world
    fam = {r["family_id"] for r in t["truth"]["truth_family"]}
    genes = {r["family_id"] for r in t["truth"]["truth_gene"]}
    missing = {f for f in genes - fam if f != "smORF" and "+" not in f}
    assert not missing, sorted(missing)[:5]
    groups = [g for g in _tsv(out / "backbone" / "groups.tsv") if g["kind"] in ("planted", "profile_only", "vehicle")]
    core = [r for r in t["truth"]["truth_module"] if r["member_role"] == "core" and r["kind"] != "cas_operon"]
    assert len(core) == sum(len(g["members"].split(";")) for g in groups)
    kinds = collections.Counter(r["kind"] for r in t["truth"]["truth_traps"])
    assert {"promiscuous_domain", "mge_cargo", "hub_family", "decoy_profile", "clonal_burst"} <= set(kinds)
    print(f"\nfamilies {len(fam)}; module rows {len(t['truth']['truth_module'])}; traps {dict(kinds)}")


# --- runtime (§2.6) -----------------------------------------------------------------------------------------------

@needs_zombi2
def test_micro_stage_times(world, record_property):
    ew, _, times, _ = world
    for k, v in times.items():
        record_property(k.split()[0] + "_seconds", round(v, 2))
    total = sum(times.values())
    record_property("S1-S11_seconds", round(total, 2))
    print("\n" + "; ".join(f"{k} {v:.2f} s" for k, v in times.items()) + f"; S1-S11 {total:.2f} s; "
          f"counts {ew.meta['counts']}")
    assert times["S11 emit"] < 120.0
