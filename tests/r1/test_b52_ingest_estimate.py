"""Ring R1: the estimate of the ingest step against the run records of every built stage (card B-52, DEC-47 (2);
AGENTS §6).

The oracle is measured, never the formulas of `stage`: the ingest run records of the nine stages built into the data
root `/tmp/if2453/r10/root` (wilkcomp02, 2026-11-10 to 16, `runs/stage=<S>/ingest.json`, `cpu_min`), beside what each
stage read and what the root held before it, read on 2026-11-16 by the script of
`planning/measurements/2026-11-16_ingest_estimate.md` (read-only, unit `b52-table`): the bytes of FASTA of the record's
`--objects`, and the bases of the contigs the stage's `genome_input` lists that a stage ingested before it holds as a
default variant in its `contig_index`. Both are known before the step starts (the inputs step publishes
`genome_input` before ingest). No stage was ingested in batches.

The guard of a build judges a step by |record - estimate| / estimate (DEC-47: 26.612 against 33.267 is "20.005%
under"), so that is the measure here.
"""

from __future__ import annotations

import pytest

from biobox import stage as bs

pytestmark = pytest.mark.covers("B-52")

# stage: (bytes of FASTA read, bases of its listed contigs another stage held before it, the run record's cpu_min)
RECORDS = {
    "S1": (2_121_124_779, 0, 7.540384),  # 34,076 records of 27,258 distinct contigs (repeats are read, not stored)
    "P1": (1_029_384_352, 0, 4.611177),
    "P2": (9_716_599_273, 0, 38.33493),
    "P3": (10_074_413_081, 6_963, 43.889413),
    "P4": (10_101_497_437, 7_053, 42.979174),
    "P5": (10_092_827_773, 185_683, 43.854939),
    "P6": (10_087_447_605, 17_198, 41.485983),
    "H1": (1_962_116_403, 141_683_767, 7.644566),  # 514 contigs S1 holds
    "PL1": (7_425_668_755, 1_177_838_866, 26.61184),  # 7,046 contigs S1 and P1 to P6 hold (DEC-47)
}


def _off(stage: str) -> float:
    n_bytes, held, cpu_min = RECORDS[stage]
    est = bs.ingest_cpu_min(n_bytes, held_bases=held)
    return (cpu_min - est) / est


@pytest.mark.parametrize("stage", list(RECORDS))
def test_b52_each_built_stage_s_ingest_record_is_within_20_percent_of_its_estimate(stage):
    n_bytes, held, cpu_min = RECORDS[stage]
    est = bs.ingest_cpu_min(n_bytes, held_bases=held)
    assert abs(cpu_min - est) <= 0.20 * est, (stage, round(est, 3), cpu_min, f"{_off(stage):+.1%}")


def test_b52_eight_of_the_nine_records_are_within_10_percent_and_only_s1_is_not():
    """S1's FASTA repeats 6,818 records of contigs it lists (420,240,502 bases read twice, *derived*: 2,084,350,440 -
    1,664,109,938 from the table of the note), which ingest reads but does not store; which records repeat is known
    only after hashing, so its estimate counts them as stored."""
    assert [s for s in RECORDS if abs(_off(s)) > 0.10] == ["S1"]


def test_b52_without_bases_held_elsewhere_the_estimate_is_the_bytes_at_the_rate():
    """The form before card B-52 is kept for a caller that knows nothing of the root (`_refusal`'s one object, an
    ingest batch's bases read): every byte counted as stored, which errs high."""
    assert bs.ingest_cpu_min(0) == bs.INGEST_FIXED_CPU_MIN
    assert bs.ingest_cpu_min(10**9) == pytest.approx(bs.INGEST_FIXED_CPU_MIN + bs.INGEST_CPU_MIN * 10**9 / bs.INGEST_BP)
    assert bs.ingest_cpu_min(10**9, held_bases=4 * 10**8) == pytest.approx(bs.ingest_cpu_min(6 * 10**8))
    assert bs.ingest_cpu_min(10**6, held_bases=10**7) == bs.INGEST_FIXED_CPU_MIN  # never below the start-up cost
