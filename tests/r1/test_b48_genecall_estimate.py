"""Ring R1: the estimate of the gene-call step against measurements (card B-48; AGENTS §6).

The oracle is measured, never the formulas of `genecall`: the run records of the eight stages gene-called so far in
the data root `/tmp/if2453/r10/root` (wilkcomp02, 2026-11-10 to 15, `runs/stage=<S>/genecall.json`, `cpu_min`) beside
what each stage asked of gene calling (`scripts/stress/build_stage_pod/b48_stage_inputs.py`: every genome's contigs,
their lengths and topology, and the mode the code calls it in; every built stage is linear, one mode a stage), and the
calls-only samples of card C-26 (stage PL1's records, `c26_rate.py`) and of this card (`b48_call_rate.py`, every 30th
genome of S1, H1, P2 and P6). The numbers are those of `planning/measurements/2026-11-15_genecall_batches.md`.

A gene-call estimate is a plan's: a batch cut by it must not overrun it, so each estimate is at least its record and,
by AGENTS §6, not more than 20% above it (card W-35's test, `test_w35_each_estimate_is_at_or_above_its_run_record`).
"""

from __future__ import annotations

import pytest

from biobox.l1 import genecall

pytestmark = pytest.mark.covers("B-48")

# stage: (mode, bases the caller read, the run record's cpu_min). H1's 1,961,515,105 include the 514 genomes
# (141,683,767 bases) all of whose contigs another stage stores: the code of 2026-11-15 called them, and dropped the
# rows (DEC-26).
RECORDS = {
    "S1": ("gv-meta", 1_664_109_938, 142.785628),
    "H1": ("gv-meta", 1_961_515_105, 146.09093),
    "P1": ("single", 1_016_625_803, 15.391635),
    "P2": ("single", 9_596_180_328, 143.287698),
    "P3": ("single", 9_949_084_485, 153.677083),
    "P4": ("single", 9_976_464_003, 151.405728),
    "P5": ("single", 9_967_897_315, 155.692349),
    "P6": ("single", 9_962_553_301, 157.196096),
}
# calls-only samples, one thread: (mode, circular, bases called, CPU-seconds of the calls)
SAMPLES = {
    "PL1 circular": ("meta", True, 295_462_404, 1306.076),
    "PL1 linear": ("meta", False, 17_141_185, 50.48),
    "S1": ("gv-meta", False, 55_446_997, 232.98),
    "H1": ("gv-meta", False, 65_344_071, 246.97),
    "P2": ("single", False, 304_636_972, 238.786),
    "P6": ("single", False, 351_376_483, 279.027),
}


@pytest.mark.parametrize("stage", sorted(RECORDS))
def test_b48_the_one_process_estimate_of_each_built_stage_is_at_its_record_or_within_20_percent_above(stage):
    mode, called, cpu_min = RECORDS[stage]
    est = genecall.one_process_cpu_min({(mode, False): called})
    assert cpu_min <= est <= 1.2 * cpu_min, (stage, round(est, 2), cpu_min)


@pytest.mark.parametrize("name", sorted(SAMPLES))
def test_b48_the_calls_rate_of_each_mode_is_at_its_samples_and_within_20_percent(name):
    mode, circular, called, cpu_s = SAMPLES[name]
    est = genecall.CALL_S_PER_MBP[(mode, circular)] * called / 1e6
    assert cpu_s <= est <= 1.2 * cpu_s, (name, round(est, 2), cpu_s)


# Stage S1 gene-called again in 4 batches by build-stage in a scratch view of the root (wilkcomp02, 2026-11-15, unit
# b48-s1-batches, 4 workers each, card B-48's note): each batch's bases read and its `done.json` CPU-minutes (its own
# process and its worker processes), and the merge's run record with the genes it wrote.
S1_BATCHES = [(416_023_553, 30.86), (416_029_336, 30.70), (416_075_372, 31.87), (415_981_677, 30.78)]
S1_MERGE = (1_664_109_938, 2_398_087, 0.432473)


@pytest.mark.parametrize("j", range(len(S1_BATCHES)))
def test_b48_each_batch_of_s1_is_at_its_estimate_or_within_20_percent_under_it(j):
    called, cpu_min = S1_BATCHES[j]
    est = genecall.one_process_cpu_min({("gv-meta", False): called})  # a batch's estimate: its process and its calls
    assert cpu_min <= est <= 1.2 * cpu_min, (j, round(est, 2), cpu_min)


def test_b48_the_merge_of_s1_is_at_its_estimate_or_within_20_percent_under_it():
    bases, genes, cpu_min = S1_MERGE
    assert abs(genecall.genome_genes("virus", bases) - genes) / genes < 0.001  # the expected genes of the plan
    est = genecall.merge_cpu_min(genecall.genome_genes("virus", bases))
    assert cpu_min <= est <= 1.2 * cpu_min, (round(est, 3), cpu_min)


def test_b48_the_pl1_sample_as_a_whole_is_within_20_percent():
    """Card C-26's sample of PL1, 2,054 records in meta mode (circular and linear together): 1,356.556 CPU-seconds."""
    est = sum(genecall.CALL_S_PER_MBP[("meta", c)] * b / 1e6 for c, b in ((True, 295_462_404), (False, 17_141_185)))
    assert 1356.556 <= est <= 1.2 * 1356.556


def test_b48_the_mode_and_the_bases_read_follow_the_callers_rule():
    """The estimate's mode is `_finder`'s (virus gv-meta; plasmid, or a genome under SINGLE_MODE_MIN_BP, meta; else
    single), and a circular contig is read with its wrap flanks (`_genes`)."""
    assert genecall.mode_of("virus", 10**7) == "gv-meta"
    assert genecall.mode_of("plasmid", 10**7) == "meta"
    assert genecall.mode_of("prokaryote", genecall.SINGLE_MODE_MIN_BP - 1) == "meta"
    assert genecall.mode_of("prokaryote", genecall.SINGLE_MODE_MIN_BP) == "single"
    with pytest.raises(ValueError, match="unknown genome kind"):
        genecall.mode_of("bogus", 1)
    assert genecall.called_bases(30_000, True) == 90_000  # the whole contig on each side
    assert genecall.called_bases(200_000, True) == 300_000  # WRAP_FLANK on each side
    assert genecall.called_bases(200_000, False) == 200_000
    assert genecall.WRAP_FLANK == 50_000
