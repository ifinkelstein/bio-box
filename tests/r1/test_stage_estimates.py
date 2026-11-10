"""Card C-13 (5), (6): the CPU-minute estimates of `build-stage` against the first real build (card C-10).

The oracle is the run records and scan shard records of that build, as its measurements file gives them
(`planning/measurements/2026-11-10_build_stage_pod.md`, *measured* on wilkcomp02, 2026-11-10), not the formulas in
`stage.py`. AGENTS §6: an estimate more than 20% from a measurement is a finding.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from biobox import profile
from biobox import stage as bs

pytestmark = pytest.mark.covers("C-13")

# The P1 Pfam scan of C-10: 19 shard processes, 2,343.74 CPU-minutes in all, 294,516,869 residues (shards.json).
P1_SCAN_CPU_MIN, P1_RESIDUES = 2343.74, 294_516_869
# Ingest of C-10's two stages: bytes of FASTA read, CPU-minutes of the ingest run record.
S1A_BYTES, S1A_INGEST_CPU_MIN = 63_357_113, 0.36
P1_BYTES, P1_INGEST_CPU_MIN = 1_029_384_352, 4.44
# The import ladder's P2 ingest, the measurement the per-byte rate comes from (42.9 CPU-min, 9,596,180,328 bp).
P2_BYTES, P2_INGEST_CPU_MIN = 9_596_180_328, 42.9


def _profile() -> profile.Profile:
    return profile.load(Path(__file__).resolve().parents[2] / "profiles" / "wilkcomp02.yaml")


def test_c13_a_shard_planned_for_the_budget_fits_its_margin_at_the_measured_pfam_rate():
    """A shard of `shard_residues` residues, scanned at the rate the P1 scan measured (7.96 CPU-min per million
    residues, derived: 2,343.74 / 294.52), stays within `SHARD_MARGIN` (80%) of the 300 CPU-minute budget; and the
    margin leaves room for a library a further 25% slower (derived: 1 / 0.8 = 1.25) before the budget is passed."""
    prof = _profile()
    measured = P1_SCAN_CPU_MIN / (P1_RESIDUES / 1e6)
    assert bs.SCAN_CPU_MIN_PER_M >= measured  # the constant is the larger measured rate
    shard = bs.shard_residues(prof) * measured / 1e6
    assert shard <= bs.SHARD_MARGIN * prof.cpu_min_per_process
    assert shard * 1.25 <= prof.cpu_min_per_process


def test_c13_the_scan_estimate_of_the_p1_scan_is_within_20_percent_of_its_measurement():
    est = bs.scan_cpu_min(P1_RESIDUES)
    assert abs(est - P1_SCAN_CPU_MIN) / P1_SCAN_CPU_MIN <= 0.20
    assert est >= P1_SCAN_CPU_MIN * 0.999  # it errs high, if at all (the rounding of the rate)


@pytest.mark.parametrize("n_bytes,measured", [(S1A_BYTES, S1A_INGEST_CPU_MIN), (P1_BYTES, P1_INGEST_CPU_MIN),
                                              (P2_BYTES, P2_INGEST_CPU_MIN)], ids=["S1a", "P1", "P2"])
def test_c13_the_ingest_estimate_is_within_20_percent_on_a_small_and_a_large_stage(n_bytes, measured):
    """C-10 found the estimate 21% under the measurement on a 63 MB stage (0.28 against 0.36 CPU-min): a step's fixed
    start-up cost dominates it. With the measured constant term the estimate is within 20% at every size measured."""
    est = bs.ingest_cpu_min(n_bytes)
    assert abs(est - measured) / measured <= 0.20, (n_bytes, est, measured)


def test_c13_the_ingest_estimate_has_a_constant_term():
    """The start-up cost of a step that does almost nothing: C-10's `inputs` step of S1a, 0.07 CPU-min."""
    assert bs.ingest_cpu_min(0) == bs.INGEST_FIXED_CPU_MIN == 0.07
