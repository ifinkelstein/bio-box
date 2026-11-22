"""Ring R1: meta-mode long-linear rate and short-batch accumulation (card B-63; AGENTS §6).

Oracles are measured on wilkcomp02 (2026-11-21, `scripts/stress/b63_meta_genecall_arms.py`,
`planning/measurements/2026-11-21_b63_meta_genecall.md`):
- Arm A (`bio-box_faec1f4.sif`, 1 thread): meta-mode CPU-s/Mbp increases 2.4x-3.9x with contig length
  (1-5 kb: 1.38-2.01 s/Mbp; 10-30 kb: 1.89-3.20 s/Mbp; >=50 kb: 3.27-5.76 s/Mbp); card C-40's 1-in-50
  IMG-PR linear sample ran at 3.9980 s/Mbp (`2026-11-16_pl2_inputs.md` §5).
- Arm B (4 workers, `stage genecall` on 1-in-200 IMG-PR + mMGE = 3,742 genomes, 82.58 Mbp):
  - Linear IMG-PR + mMGE (38 circular, 11.07 Mbp <10 kb, 71.51 Mbp >=10 kb): 7.637 CPU-min baseline.
  - DEC-63 circular-DTR IMG-PR + mMGE (684 circular): 8.643 CPU-min baseline.
"""

from __future__ import annotations

import pytest

from biobox.l0seq.inputs import GenomeInput
from biobox.l1 import genecall

pytestmark = pytest.mark.covers("B-63")

# Card C-40 1-in-50 IMG-PR linear plasmid sample (wilkcomp02, 2026-11-17, `2026-11-16_pl2_inputs.md` §5):
# 13,115 linear genomes, 313.52 Mbp, mean 23,905 bp, 3.9980 CPU-s/Mbp.
IMGPR_SAMPLE_MEAN_BP = 23_905
IMGPR_SAMPLE_CALL_S_PER_MBP = 3.9980


def test_b63_long_linear_meta_rate_covers_imgpr_sample_within_20_percent():
    """A linear meta-mode contig >= 10 kb must be priced at or above IMG-PR's measured 3.9980 s/Mbp
    and within 20% of it (AGENTS §6)."""
    est_s_per_mbp = genecall.calls_cpu_s("meta", [(IMGPR_SAMPLE_MEAN_BP, False)]) / (IMGPR_SAMPLE_MEAN_BP / 1e6)
    assert IMGPR_SAMPLE_CALL_S_PER_MBP <= est_s_per_mbp <= 1.20 * IMGPR_SAMPLE_CALL_S_PER_MBP, (
        round(est_s_per_mbp, 4),
        IMGPR_SAMPLE_CALL_S_PER_MBP,
    )


def test_b63_long_linear_meta_factor_scales_with_monkeypatched_call_s_per_mbp(monkeypatch):
    """Monkeypatching `genecall.CALL_S_PER_MBP` (as `test_b48_genecall_batches.py` and `b50_batched_tree.py` do)
    must scale long linear meta contigs proportionally."""
    base_calls = genecall.calls_cpu_s("meta", [(50_000, False)])
    base_cost = genecall.genome_cost_us("plasmid", [(50_000, False)])
    monkeypatch.setattr(genecall, "CALL_S_PER_MBP", {m: r * 10.0 for m, r in genecall.CALL_S_PER_MBP.items()})
    assert genecall.calls_cpu_s("meta", [(50_000, False)]) == pytest.approx(base_calls * 10.0)
    assert genecall.genome_cost_us("plasmid", [(50_000, False)]) == pytest.approx(base_cost * 10.0, abs=10)


def test_b63_short_linear_meta_circular_meta_and_other_modes_unchanged():
    """Linear meta contigs below `META_LINEAR_LONG_BP`, circular meta contigs, `gv-meta`, and `single`
    keep their exact `CALL_S_PER_MBP` rates."""
    # Short linear meta (< 10 kb)
    n_short = genecall.META_LINEAR_LONG_BP - 1
    assert genecall.calls_cpu_s("meta", [(n_short, False)]) == pytest.approx(
        genecall.CALL_S_PER_MBP[("meta", False)] * n_short / 1e6
    )
    # Circular meta (any length)
    for n in (3_000, 25_000, 120_000):
        assert genecall.calls_cpu_s("meta", [(n, True)]) == pytest.approx(
            genecall.CALL_S_PER_MBP[("meta", True)] * genecall.called_bases(n, True) / 1e6
        )
    # gv-meta and single (linear and circular)
    for mode in ("gv-meta", "single"):
        for circ in (False, True):
            n = 150_000
            assert genecall.calls_cpu_s(mode, [(n, circ)]) == pytest.approx(
                genecall.CALL_S_PER_MBP[(mode, circ)] * genecall.called_bases(n, circ) / 1e6
            )


def test_b63_pl1_batch_records_stay_within_20_percent_with_long_linear_meta_rate():
    """Even if every linear contig in PL1's four batches were >= 10 kb, each batch's estimate remains
    at or above its measured run record and within 20% (AGENTS §6)."""
    pl1_batches = [
        (110_736_578, 2_332_913_588, 213.891242),
        (113_188_587, 2_331_391_450, 212.302178),
        (121_049_374, 2_326_070_477, 209.892586),
        (126_487_304, 2_322_450_973, 213.041484),
    ]
    for lin_bp, circ_bp, cpu_min in pl1_batches:
        # Represent the batch's linear bases as >=10 kb contigs and circular bases as called_bases
        calls_s = (
            genecall.calls_cpu_s("meta", [(lin_bp, False)])
            + genecall.CALL_S_PER_MBP[("meta", True)] * circ_bp / 1e6
        )
        est = (genecall.PROCESS_S + genecall.STEP_OVERHEAD * calls_s) / 60.0
        assert cpu_min <= est <= 1.20 * cpu_min, (round(est, 3), cpu_min)


def test_b63_interleaved_long_genomes_do_not_fragment_short_batches():
    """Interleaving >= thread_min_bases genomes with < thread_min_bases genomes must NOT flush the
    short-genome batch early; short genomes accumulate up to `batch_bases`."""
    short_gs = [GenomeInput((i.to_bytes(16, "big"),), 11, kind="plasmid") for i in range(1, 201)]
    long_gs = [GenomeInput(((1000 + i).to_bytes(16, "big"),), 11, kind="plasmid") for i in range(1, 201)]
    lengths = {g.contig_uids[0]: 3_000 for g in short_gs}
    lengths.update({g.contig_uids[0]: 20_000 for g in long_gs})
    interleaved = [g for pair in zip(short_gs, long_gs) for g in pair]

    out = list(genecall._batches(interleaved, lambda g: lengths[g.contig_uids[0]], 500_000, 10_000))
    long_tasks = [b for is_long, b in out if is_long]
    short_batches = [b for is_long, b in out if not is_long]

    assert len(long_tasks) == 200
    assert [b[0] for b in long_tasks] == long_gs
    # 200 * 3,000 = 600,000 bp at batch_bases=500,000 -> 2 batches (166 + 34 genomes), not 200 singletons
    assert len(short_batches) == 2, f"expected 2 short batches, got {len(short_batches)}"
    assert [g for b in short_batches for g in b] == short_gs
    assert all(sum(lengths[g.contig_uids[0]] for g in b) <= 500_000 for b in short_batches)


def test_b63_cut_plan_records_long_linear_meta_parameters_in_rates():
    """`cut_plan` records `meta_linear_long_bp` and `meta_linear_long_factor` under `rates` while
    leaving `rates['call_s_per_mbp']` unchanged."""
    items = lambda: [(b"\x01" * 16, genecall.genome_cost_us("plasmid", [(25_000, False)]), 25_000, 27.0, True)]
    p = genecall.cut_plan(items, batch_cpu_min=240.0)
    assert p["rates"]["call_s_per_mbp"]["meta/linear"] == genecall.CALL_S_PER_MBP[("meta", False)]
    assert p["rates"]["meta_linear_long_bp"] == genecall.META_LINEAR_LONG_BP
    assert p["rates"]["meta_linear_long_factor"] == genecall.META_LINEAR_LONG_FACTOR
