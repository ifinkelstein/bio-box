"""Ring R1: unit tests for card C-47 (DEC-32 / DEC-35 NR90 representative Pfam inference evaluation
and DEC-43 / DEC-53 build-driver rate guard).

Tests written FIRST before implementing `scripts/stress/c47_nr90_pfam_eval.py`.
Uses only `scripts/` (copied into the test image by `container/Dockerfile`, `tests/r0/test_image_files.py`).
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-47")  # ring R1 by directory (tests/conftest.py)

ROOT = Path(__file__).resolve().parents[2]
EVAL_SCRIPT = ROOT / "scripts" / "stress" / "c47_nr90_pfam_eval.py"
EST_SCRIPT = ROOT / "scripts" / "stress" / "build_stage_pod" / "c47_est.py"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None, f"missing script: {path}"
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _brute_force_map(cigar: bytes, qs: int, ts: int, env_start: int, env_end: int) -> tuple[int, int] | None:
    """Independent residue-by-residue oracle for B-47 §3 envelope mapping (query=rep: M both, I rep, D member)."""
    import re

    r_to_m: dict[int, int] = {}
    rp, mp = ts, qs
    for num, op in re.findall(rb"(\d+)([MID])", cigar):
        n = int(num)
        if op == b"M":
            for k in range(n):
                r_to_m[rp + k] = mp + k
            rp += n
            mp += n
        elif op == b"I":
            rp += n
        else:
            mp += n
    aligned = [r_to_m[r] for r in range(env_start, env_end) if r in r_to_m]
    if not aligned:
        return None
    return aligned[0], aligned[-1] + 1


def test_c47_dec35_length_and_coverage_classification():
    m = _load("c47_nr90_pfam_eval", EVAL_SCRIPT)
    # Exact 20% boundary is NOT scanned for length (|L_m - L_r| > 0.2 * L_r is strict)
    assert m.length_mismatch(80, 100) is False
    assert m.length_mismatch(120, 100) is False
    assert m.length_mismatch(79, 100) is True
    assert m.length_mismatch(121, 100) is True

    # Envelope coverage f inside representative aligned span [tstart, tend)
    assert m.envelope_coverage(10, 60, 0, 100) == pytest.approx(1.0)
    assert m.envelope_coverage(10, 60, 20, 100) == pytest.approx(0.8)
    assert m.envelope_coverage(10, 60, 21, 100) == pytest.approx(0.78)
    assert m.envelope_coverage(10, 60, 60, 100) == pytest.approx(0.0)

    # DEC-35 (4): if ANY domain row has f < 0.8, the member is scanned (role='unmapped'), not inferred
    assert m.classify_member(is_rep=True, lm=100, lr=100, rep_envs=[(10, 60)], tstart=0, tend=100) == "rep"
    assert m.classify_member(is_rep=False, lm=75, lr=100, rep_envs=[(10, 60)], tstart=0, tend=100) == "length"
    assert m.classify_member(is_rep=False, lm=95, lr=100, rep_envs=[(10, 60), (60, 100)], tstart=0, tend=90) == "unmapped"
    assert m.classify_member(is_rep=False, lm=95, lr=100, rep_envs=[(10, 60), (60, 100)], tstart=0, tend=92) == "inferred"
    # Representative with 0 domains -> inferred (no unmapped domains)
    assert m.classify_member(is_rep=False, lm=95, lr=100, rep_envs=[], tstart=0, tend=95) == "inferred"
    # Missing alignment for non-rep -> scanned as unmapped
    assert m.classify_member(is_rep=False, lm=95, lr=100, rep_envs=[(10, 60)], tstart=None, tend=None) == "unmapped"


def test_c47_cigar_envelope_mapping_matches_brute_force_oracle():
    m = _load("c47_nr90_pfam_eval", EVAL_SCRIPT)
    cases = [
        # (cigar, qs, qe, ts, te, env_start, env_end)
        (b"50M", 5, 55, 10, 60, 15, 45),
        # Insertion in rep (I) splits envelope start/end -> moves inward to next aligned residue
        (b"10M5I20M", 0, 30, 0, 35, 12, 30),
        (b"10M5I20M", 0, 30, 0, 35, 2, 13),
        # Deletion in rep / insertion in member (D) inside envelope is included
        (b"10M4D20M", 3, 37, 1, 31, 5, 25),
        # Clipped start/end when f >= 0.8 (envelope slightly outside [ts, te))
        (b"30M", 10, 40, 20, 50, 18, 48),
        (b"12M2I8M3D10M", 4, 37, 6, 38, 10, 36),
    ]
    for cigar, qs, qe, ts, te, es, ee in cases:
        got = m.map_envelope(cigar, qs, qe, ts, te, es, ee)
        want = _brute_force_map(cigar, qs, ts, es, ee)
        assert got == want, f"mismatch for {cigar!r} [{es},{ee}): got {got}, oracle {want}"


def test_c47_pair_evaluation_phylum_b47b_and_stop_gate():
    m = _load("c47_nr90_pfam_eval", EVAL_SCRIPT)
    # 5 proteins in 2 clusters across 2 phyla:
    # p1 (rep, len 100, p__A): {PF00001, PF00002} (copied via B-47b)
    # p2 (member of p1, len 98, p__A, f=1.0): true {PF00001} -> inferred {PF00001, PF00002} => 1 right, 1 false
    # p3 (member of p1, len 70, p__A, length mismatch): true {PF00001, PF00003} -> scanned => 2 right, 0 false, 0 missed
    # p4 (rep, len 200, p__B): {PF00010}
    # p5 (member of p4, len 195, p__B, f=0.5 < 0.8 -> unmapped scanned): true {PF00010, PF00011} -> scanned => 2 right
    # p6 (member of p4, len 200, p__B, f=0.9 >= 0.8 -> inferred): true {PF00010, PF00012} -> inferred {PF00010} => 1 right, 1 missed
    proteins = {
        "p1": {"len": 100, "phylum": "p__A", "copied": True},
        "p2": {"len": 98, "phylum": "p__A", "copied": False},
        "p3": {"len": 70, "phylum": "p__A", "copied": False},
        "p4": {"len": 200, "phylum": "p__B", "copied": False},
        "p5": {"len": 195, "phylum": "p__B", "copied": False},
        "p6": {"len": 200, "phylum": "p__B", "copied": False},
    }
    clusters = {
        "p1": ("p1", None, None),
        "p2": ("p1", 0, 100),
        "p3": ("p1", 0, 70),
        "p4": ("p4", None, None),
        "p5": ("p4", 0, 100),   # rep domain [50, 150) has overlap 50/100 = 0.5 < 0.8 -> unmapped scanned
        "p6": ("p4", 60, 200),  # rep domain [50, 150) has overlap 90/100 = 0.9 >= 0.8 -> inferred
    }
    domains = {
        "p1": [("PF00001", 10, 50), ("PF00002", 60, 90)],
        "p2": [("PF00001", 10, 50)],
        "p3": [("PF00001", 10, 50), ("PF00003", 55, 68)],
        "p4": [("PF00010", 50, 150)],
        "p5": [("PF00010", 50, 150), ("PF00011", 155, 190)],
        "p6": [("PF00010", 50, 150), ("PF00012", 155, 195)],
    }
    res = m.evaluate_clusters(proteins, clusters, domains)
    ov = res["overall"]
    # True pairs total: 2 + 1 + 2 + 1 + 2 + 2 = 10
    # Right pairs: p1(2) + p2(1) + p3(2) + p4(1) + p5(2) + p6(1) = 9
    # False pairs: p2(1) = 1; Missed pairs: p6(1) = 1
    assert ov["true_pairs"] == 10
    assert ov["pairs_right"] == 9
    assert ov["pairs_false"] == 1
    assert ov["pairs_missed"] == 1
    assert ov["wrong_rate"] == pytest.approx(1 / 10)
    assert ov["missed_rate"] == pytest.approx(1 / 10)
    assert ov["stop_triggered"] is True  # both 10% > 1%
    # Scanned (before B-47b): p1(100) + p3(70) + p4(200) + p5(195) = 565 of 863 residues
    assert ov["scanned_residues"] == 565
    assert ov["total_residues"] == 863
    # After B-47b copy (p1 is copied, not scanned): 565 - 100 = 465 scanned
    assert ov["scanned_residues_after_b47b"] == 465
    assert ov["copied_residues"] == 100
    assert set(res["by_phylum"]) == {"p__A", "p__B"}
    assert res["by_phylum"]["p__A"]["pairs_false"] == 1
    assert res["by_phylum"]["p__B"]["pairs_missed"] == 1


def _uid(name: str) -> bytes:
    return name.encode().ljust(16, b"\0")


def test_c47_eval_main_on_a_root_with_the_registry_columns(tmp_path):
    """`main` end to end on a small root whose tables carry the registry's column names (l0_manifest's object_path and
    object_sha256, contig_occurrence's l0_object_sha256, feature_protein, protein_seq, l1_domains, cluster_member,
    cluster_member_aln), each read by its run record's file list. The scenario is the one of the in-memory test above;
    its counts were worked out by hand there."""
    import json

    import pyarrow as pa
    import pyarrow.parquet as pq

    m = _load("c47_nr90_pfam_eval", EVAL_SCRIPT)
    from biobox.l1.genecall import GENE_CALL_VERSION as GCV

    root = tmp_path / "root"
    U = pa.binary(16)

    def put(rel: str, table: pa.Table) -> str:
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(table, p)
        return rel

    def record(stage: str, step: str, rels: list[str]) -> None:
        p = root / "runs" / f"stage={stage}" / f"{step}.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps({"outputs": [{"path": r} for r in rels]}))

    lens = {"p1": 100, "p2": 98, "p3": 70, "p4": 200, "p5": 195, "p6": 200}
    names = sorted(lens)
    ps = put(f"l0seq/protein_seq/gene_call_version={GCV}/stage=C5/protein_prefix=00/part-0.parquet", pa.table({
        "protein_uid": pa.array([_uid(n) for n in names], U), "protein_digest": pa.array([b"\0" * 32] * 6, pa.binary(32)),
        "length": pa.array([lens[n] for n in names], pa.uint32()), "seq": ["M" * lens[n] for n in names]}))
    contig = {"p1": "cA", "p2": "cA", "p3": "cA", "p4": "cB", "p5": "cB", "p6": "cB"}
    fp = put(f"l1/feature_protein/gene_call_version={GCV}/stage=C5/part-0.parquet", pa.table({
        "contig_uid": pa.array([_uid(contig[n]) for n in names], U), "feature_uid": pa.array([_uid("f" + n) for n in names], U),
        "start": pa.array([10 * k for k in range(6)], pa.uint64()), "end": pa.array([10 * k + 9 for k in range(6)], pa.uint64()),
        "strand": pa.array([1] * 6, pa.int8()), "protein_uid": pa.array([_uid(n) for n in names], U),
        "genetic_code": pa.array([11] * 6, pa.uint8()), "edit_type": ["none"] * 6}))
    record("C5", "genecall", [ps, fp])
    sha = {"GCA_A": b"\x01" * 32, "GCA_B": b"\x02" * 32}
    man = put("l0/l0_manifest/stage=C5/part-0.parquet", pa.table({
        "object_path": ["GCA_A.fna.gz", "GCA_B.fna.gz"], "object_sha256": pa.array([sha["GCA_A"], sha["GCA_B"]], pa.binary(32))}))
    occ = put("l0seq/contig_occurrence/stage=C5/source=ncbi_genomes/part-0.parquet", pa.table({
        "l0_object_sha256": pa.array([sha["GCA_A"], sha["GCA_B"]], pa.binary(32)),
        "record_ordinal": pa.array([0, 0], pa.uint64()), "contig_uid": pa.array([_uid("cA"), _uid("cB")], U)}))
    record("C5", "ingest", [man, occ])
    doms = [("p1", "PF00001", 10, 50), ("p1", "PF00002", 60, 90), ("p2", "PF00001", 10, 50), ("p3", "PF00001", 10, 50),
            ("p3", "PF00003", 55, 68), ("p4", "PF00010", 50, 150), ("p5", "PF00010", 50, 150), ("p5", "PF00011", 155, 190),
            ("p6", "PF00010", 50, 150), ("p6", "PF00012", 155, 195)]
    dom = put("l1/domains/library=pfam-38.2/stage=C5/part-0.parquet", pa.table({
        "protein_uid": pa.array([_uid(d[0]) for d in doms], U), "model_acc": [d[1] for d in doms],
        "env_start": pa.array([d[2] for d in doms], pa.uint32()), "env_end": pa.array([d[3] for d in doms], pa.uint32())}))
    record("C5", "domains-pfam-38.2", [dom])
    src = put(f"l0seq/protein_seq/gene_call_version={GCV}/stage=P1/protein_prefix=00/part-0.parquet",
              pa.table({"protein_uid": pa.array([_uid("p1")], U)}))
    plan = root / "l1/domain_scan/library=pfam-38.2/stage=C5/plan.json"
    plan.parent.mkdir(parents=True, exist_ok=True)
    plan.write_text(json.dumps({"copy": {"sources": [{"stage": "P1", "protein_files": [{"path": src}]}]}}))
    rep = {"p1": "p1", "p2": "p1", "p3": "p1", "p4": "p4", "p5": "p4", "p6": "p4"}
    cm = tmp_path / "cm.parquet"
    pq.write_table(pa.table({"member_uid": pa.array([_uid(n) for n in names], U), "level": ["nr90"] * 6,
                             "cluster_uid": pa.array([_uid("c" + rep[n]) for n in names], U),
                             "rep_uid": pa.array([_uid(rep[n]) for n in names], U), "cluster_params_hash": ["h"] * 6}), cm)
    # q* member, t* representative; the CIGAR's query is the representative (I: representative only, D: member only)
    aln_rows = [("p1", 0, 100, 0, 100, b"100M"), ("p2", 0, 98, 0, 100, b"50M2I48M"), ("p3", 0, 70, 0, 70, b"70M"),
                ("p4", 0, 200, 0, 200, b"200M"), ("p5", 0, 100, 0, 100, b"100M"), ("p6", 60, 200, 60, 200, b"140M")]
    aln = tmp_path / "aln.parquet"
    pq.write_table(pa.table({
        "member_uid": pa.array([_uid(r[0]) for r in aln_rows], U), "level": ["nr90"] * 6,
        "rep_uid": pa.array([_uid(rep[r[0]]) for r in aln_rows], U),
        **{k: pa.array([r[i] for r in aln_rows], pa.uint32()) for i, k in enumerate(("qstart", "qend", "tstart", "tend"), 1)},
        "pident": pa.array([1.0] * 6, pa.float32()), "alnlen": pa.array([100] * 6, pa.uint32()),
        "cigar": pa.array([r[5] for r in aln_rows], pa.binary())}), aln)
    sel = tmp_path / "selection.tsv"
    sel.write_text("accession\tphylum\nGCA_A\tp__A\nGCA_B\tp__B\n")
    objs = tmp_path / "objects.json"
    objs.write_text(json.dumps([{"path": "/x/GCA_A.fna.gz"}, {"path": "/x/GCA_B.fna.gz"}]))
    out = tmp_path / "out.json"
    assert m.main(["eval", "--root", str(root), "--stage", "C5", "--stage-cm", str(cm), "--stage-aln", str(aln),
                   "--selection", str(sel), "--objects", str(objs), "--out", str(out), "--tmp", str(tmp_path / "t")]) == 0
    doc = json.loads(out.read_text())
    ov = doc["stage_local"]["overall"]
    assert (ov["true_pairs"], ov["pairs_right"], ov["pairs_false"], ov["pairs_missed"]) == (10, 9, 1, 1)
    assert (ov["scanned_residues"], ov["total_residues"], ov["scanned_residues_after_b47b"]) == (565, 863, 465)
    assert (ov["reps"], ov["length_scanned"], ov["unmapped_scanned"], ov["inferred_members"]) == (2, 1, 1, 2)
    assert doc["dec35_6_stop_triggered"] is True
    ph = doc["stage_local"]["by_phylum"]
    assert set(ph) == {"p__A", "p__B"}
    assert (ph["p__A"]["pairs_false"], ph["p__A"]["pairs_missed"]) == (1, 0)
    assert (ph["p__B"]["pairs_false"], ph["p__B"]["pairs_missed"]) == (0, 1)
    # p2's PF00001 maps onto itself; p6's PF00010 [50, 150) maps from the alignment's first column (60), 10 from its own
    co = doc["stage_local"]["coordinates"]
    assert co["n_pairs"] == 2 and co["both_within_10"] == 1.0 and co["both_within_5"] == 0.5


def test_c47_est_dec53_rate_guard():
    est = _load("c47_est", EST_SCRIPT)
    # DEC-43: < 24 CPU-min is never judged even if ratio is 2.0
    row_small, over, under = est.judge_rate("crispr-default", 20.0, 10.0, "test")
    assert row_small["judged"] is False
    assert over is False and under is False

    # DEC-53: >= 24 CPU-min and ratio < 0.80 -> under_20pct=True, over_20pct=False (no stop!)
    row_under, over, under = est.judge_rate("ingest", 18.0, 30.0, "test")
    assert row_under["judged"] is True
    assert over is False and under is True

    # DEC-53: >= 24 CPU-min and ratio > 1.20 -> over_20pct=True (STOP!)
    row_over, over, under = est.judge_rate("domains.shard-00000", 290.0, 240.0, "test")
    assert row_over["judged"] is True
    assert over is True and under is False

    # Within +-20% -> neither
    row_ok, over, under = est.judge_rate("genecall", 150.0, 152.8, "test")
    assert row_ok["judged"] is True
    assert over is False and under is False
