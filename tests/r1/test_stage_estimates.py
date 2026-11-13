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


# --- clustering (card B-36, DEC-25) -----------------------------------------------------------------------------
# The clustering run records of the 10 Gbp rebuild, 16 threads each (*measured*, wilkcomp02, 2026-11-10, branch
# `build/r10`, planning/measurements/2026-11-10_rebuild_10gbp.md, "Per-step numbers, phase 1"): proteins clustered,
# CPU-minutes of the run record. `rN.1a` is phages only (stage S1).
RN1A = (1_134_550, 84.75)
RN1B = (2_053_075, 284.10)
RN1 = (10_620_684, 2195.91)


@pytest.mark.covers("B-36")
@pytest.mark.parametrize("n,measured", [RN1B, RN1], ids=["rN.1b", "rN.1"])
def test_b36_the_cluster_estimate_is_within_20_percent_of_the_two_upper_points(n, measured):
    est = bs.cluster_cpu_min(n)
    assert abs(est - measured) / measured <= 0.20, (n, est, measured)


@pytest.mark.covers("B-36")
@pytest.mark.parametrize("n,measured", [RN1A, RN1B, RN1], ids=["rN.1a", "rN.1b", "rN.1"])
def test_b36_the_cluster_estimate_is_below_none_of_the_three_measured_points(n, measured):
    """A refusal gate errs high: the estimate is at least each measurement (to the rounding of the records)."""
    assert bs.cluster_cpu_min(n) >= measured * 0.999, (n, bs.cluster_cpu_min(n), measured)


@pytest.mark.covers("B-36")
def test_b36_the_shipped_profile_caps_clustering_at_20000_cpu_minutes():
    assert _profile().cluster_cpu_min_max == 20000


def _write_profile(tmp_path: Path, **change) -> Path:
    import yaml

    d = yaml.safe_load((Path(__file__).resolve().parents[2] / "profiles" / "wilkcomp02.yaml").read_text())
    d.update(change)
    for k, v in list(d.items()):
        if v is ...:
            del d[k]
    p = tmp_path / "p.yaml"
    p.write_text(yaml.safe_dump(d))
    return p


@pytest.mark.covers("B-36")
def test_b36_the_cluster_cap_is_optional_may_be_null_and_must_be_an_int(tmp_path):
    assert profile.load(_write_profile(tmp_path, cluster_cpu_min_max=...)).cluster_cpu_min_max is None
    assert profile.load(_write_profile(tmp_path, cluster_cpu_min_max=None)).cluster_cpu_min_max is None
    for bad in ("12000", True, 1.5):
        with pytest.raises(ValueError, match="cluster_cpu_min_max"):
            profile.load(_write_profile(tmp_path, cluster_cpu_min_max=bad))
    with pytest.raises(ValueError, match="cpu_min_per_process"):  # null is allowed for the optional key only
        profile.load(_write_profile(tmp_path, cpu_min_per_process=None))


@pytest.mark.covers("B-36", "B-42")
def test_b36_with_the_cap_only_a_call_run_as_one_process_may_pass_300_up_to_the_cap():
    """Card B-42 (DEC-29 narrows DEC-25): the cap no longer lets the whole release cluster in one process over 300
    CPU-minutes; that is refused whatever the cap, and the release clusters in parts. Over 300 only a call that no
    number of parts fits may run, as one process, up to the cap, and above it is refused, naming the cap."""
    import dataclasses

    prof = _profile()
    assert bs.cluster_cpu_min(RN1[0]) > prof.cpu_min_per_process
    for cap in (20000, 1000):
        why = bs.cluster_refusal(dataclasses.replace(prof, cluster_cpu_min_max=cap), "rN.1", RN1[0])
        assert why is not None and "10,620,684" in why and "cpu_min_per_process" in why  # in one process: refused
        assert "DEC-29" in why and "in parts" in why and "--cluster-parts" in why  # names the split path
    assert bs.cluster_refusal(prof, "rN.1", RN1[0], one_process=[]) is None  # in parts: nothing passes 300
    call = "cluster.nr90.11.prefilter"
    assert bs.cluster_refusal(prof, "rN.1", RN1[0], one_process=[(call, RN1[1])]) is None  # one call, under the cap
    capped = dataclasses.replace(prof, cluster_cpu_min_max=1000)
    why = bs.cluster_refusal(capped, "rN.1", RN1[0], one_process=[(call, RN1[1])])
    assert why is not None and "cluster_cpu_min_max" in why and "1000" in why and call in why


@pytest.mark.covers("B-36", "B-42")
def test_b36_without_the_cap_clustering_over_the_budget_in_one_process_is_refused():
    import dataclasses

    prof = dataclasses.replace(_profile(), cluster_cpu_min_max=None)
    why = bs.cluster_refusal(prof, "rN.1", RN1[0])
    assert why is not None and "cpu_min_per_process" in why and "cluster_cpu_min_max" in why and "DEC-29" in why
    small = 1000  # proteins: far under the budget at any fitted rate
    assert bs.cluster_cpu_min(small) < prof.cpu_min_per_process
    assert bs.cluster_refusal(prof, "rN.1", small) is None
