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
