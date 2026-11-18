"""Ring R1: refitted estimates of the CRISPR steps on circular plasmids and gene calling on circular viruses
(card B-53; AGENTS §6).

The oracles are measured, never the formulas under test:

- the CRISPR run records of all nine stages built into `/tmp/if2453/r10/root` (wilkcomp02, 2026-11-15 to 16,
  `runs/stage=<S>/crispr-*.json`, `cpu_min`), for both the default producer (`minced-0.4.2-w8-nr3-r23x47-s26x50`) and
  the filtered two-repeat producer (`...-nr2-...-cas1kb-caspfam1-pfam-38.2`), beside each stage's stored genome bases
  (`contig_index` default variants) and the circular wrap-flank bases MinCED scans around its circular contigs
  (`2 * min(CIRC_FLANK, contig_len)`, read on 2026-11-17 by `b53_table.py` in
  `planning/measurements/2026-11-17_crispr_genecall_estimates.md`);
- the gene-calling batch run records of stage `PL1` (`runs/stage=PL1/genecall.batch-0000*.json`, wilkcomp02,
  2026-11-16, card C-26 part 2; the eight one-process stages are held by `test_b48_genecall_estimate.py`);
- the one-thread `genecall._call` sample of stage `V1` (every 30th stored contig of UHGV v1.0,
  `/tmp/if2453/s2/inputs/V1/`, measured on wilkcomp02, 2026-11-16, card C-39, `c39_est.py`,
  `planning/measurements/2026-11-16_virus_votu_inputs.md` §5 and `V1/est.json`), paired with the sample's `called_bases`
  (with 50 kb wrap flanks) read on 2026-11-17 (`b53_table.py`).

As in `test_b52_ingest_estimate.py`, a step's guard judges |record - estimate| / estimate <= 0.20 (DEC-43, DEC-47);
a gene-call rate or batch estimate must also sit at or above its measurement so batches do not overrun the 300
CPU-minute cap (`test_b48_genecall_estimate.py`).
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from biobox import fasta, ids, profile
from biobox import stage as bs
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l1 import crispr, genecall

pytestmark = pytest.mark.covers("B-53")

# stage: (stored_bases, circular_flank_bases, default_cpu_min, near_cas_cpu_min)
# Read from /tmp/if2453/r10/root on wilkcomp02, 2026-11-17 (`b53_table.py`). Stages S1, P1..P6 and H1 list every
# contig linear; PL1 stores 61,647 contigs (5,924,552,875 bp), 56,548 of them circular (5,453,091,032 bp), adding
# 6,205,263,778 wrap-flank bases at `crispr.CIRC_FLANK` = 100,000 (12,129,816,653 bases scanned by MinCED).
CRISPR_RECORDS = {
    "S1": (1_664_109_938, 0, 5.085492, 5.486575),
    "P1": (1_016_625_803, 0, 3.084777, 3.382703),
    "P2": (9_596_180_328, 0, 26.48908, 31.794481),
    "P3": (9_949_077_522, 0, 29.194508, 33.385812),
    "P4": (9_976_456_950, 0, 27.585047, 32.606163),
    "P5": (9_967_711_632, 0, 27.505385, 32.360891),
    "P6": (9_962_536_103, 0, 27.937717, 32.796186),
    "H1": (1_819_831_338, 0, 4.924614, 5.876108),
    "PL1": (5_924_552_875, 6_205_263_778, 30.938485, 35.513282),
}

# Stage PL1's four gene-call batches (`runs/stage=PL1/genecall.batch-0000*.json`, wilkcomp02, 2026-11-16):
# (linear_called_bases, circular_called_bases, run_record_cpu_min), all in `meta` mode.
PL1_BATCHES = [
    (110_736_578, 2_332_913_588, 213.891242),
    (113_188_587, 2_331_391_450, 212.302178),
    (121_049_374, 2_326_070_477, 209.892586),
    (126_487_304, 2_322_450_973, 213.041484),
]

# Card C-39's gene-calling sample of stage V1 (`V1/est.json`, wilkcomp02, 2026-11-16, one thread, `genecall._call`
# timed alone on every 30th stored contig) paired with `called_bases` (`b53_table.py`, 2026-11-17):
# 778 circular virus genomes (34,564,134 contig bp, 88,404,110 bp read with 50 kb wrap flanks, 386.87 CPU-s) and
# 4,833 linear virus genomes (137,465,952 bp, 511.67 CPU-s).
V1_CIRC_CONTIG_BP, V1_CIRC_CALLED_BP, V1_CIRC_CPU_S = 34_564_134, 88_404_110, 386.87
V1_LIN_CONTIG_BP, V1_LIN_CALLED_BP, V1_LIN_CPU_S = 137_465_952, 137_465_952, 511.67
# Whole V1 inputs (168,536 virus genomes: 146,019 linear, 4,233,918,303 bp; 22,517 circular, 988,370,473 contig bp,
# 2,494,857,527 bp read with 50 kb wrap flanks).
V1_ALL_LIN_CALLED_BP = 4_233_918_303
V1_ALL_CIRC_CONTIG_BP, V1_ALL_CIRC_CALLED_BP = 988_370_473, 2_494_857_527


def _crispr_off(stage: str, near_cas: bool) -> float:
    bases, flank, def_min, near_min = CRISPR_RECORDS[stage]
    cpu_min = near_min if near_cas else def_min
    est = bs.crispr_cpu_min(bases, flank_bases=flank, near_cas=near_cas)
    return (cpu_min - est) / est


@pytest.mark.parametrize("near_cas", [False, True], ids=["default", "near_cas"])
@pytest.mark.parametrize("stage", ["PL1", "S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1"])
def test_b53_each_built_stage_s_crispr_record_is_within_20_percent_of_its_estimate(stage, near_cas):
    bases, flank, def_min, near_min = CRISPR_RECORDS[stage]
    cpu_min = near_min if near_cas else def_min
    est = bs.crispr_cpu_min(bases, flank_bases=flank, near_cas=near_cas)
    assert abs(cpu_min - est) <= 0.20 * est, (stage, near_cas, round(est, 3), cpu_min, f"{_crispr_off(stage, near_cas):+.1%}")


def test_b53_all_eighteen_crispr_records_are_within_10_percent():
    bad = [(s, near, f"{_crispr_off(s, near):+.2%}")
           for s in CRISPR_RECORDS for near in (False, True) if abs(_crispr_off(s, near)) > 0.10]
    assert bad == []


def test_b53_without_flank_bases_crispr_cpu_min_keeps_its_one_argument_form():
    """Callers outside `stage.py` (`c24_est.py`, `c26_est.py`, `c39_est.py`, `test_crispr_parse.py`) pass `bases`
    alone: every base scanned is counted at `CRISPR_CPU_S_PER_GBP`."""
    assert bs.crispr_cpu_min(0) == 0.0
    assert bs.crispr_cpu_min(10**9) == pytest.approx(bs.CRISPR_CPU_S_PER_GBP / 60)
    assert bs.crispr_cpu_min(6 * 10**8, flank_bases=4 * 10**8) == pytest.approx(bs.crispr_cpu_min(10**9))
    assert bs.crispr_cpu_min(10**9, near_cas=True) == pytest.approx(bs.crispr_cpu_min(10**9) * bs.CRISPR_NEAR_CAS_COST)


def _stage_with_plasmids(tmp_path: Path) -> tuple[Path, dict[str, bytes]]:
    """Stage S0 owns circular contig `shared`; stage PL stores linear `lin`, short circular `short` (< CIRC_FLANK,
    listed in two genomes of PL under DEC-27) and long circular `long` (> CIRC_FLANK), and also lists `shared`
    (DEC-26)."""
    import json

    root = tmp_path / "root"
    raw = {
        "shared": b"ACGT" * 5_000,   # 20 kb, owned by S0: skipped by PL's CRISPR scan
        "lin": b"TGCA" * 7_500,      # 30 kb, linear in PL: 0 flank bases
        "short": b"CGTA" * 10_000,   # 40 kb, circular in PL (< 100 kb): 2 * 40 kb = 80 kb flank bases (counted once)
        "long": b"GTAC" * 62_500,    # 250 kb, circular in PL (> 100 kb): 2 * 100 kb = 200 kb flank bases
    }
    cu = {k: ids.contig_uid(v) for k, v in raw.items()}
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    for st, names, genomes in [
        ("S0", ["shared"], [GenomeInput((cu["shared"],), 11, kind="plasmid", circular=(cu["shared"],))]),
        ("PL", ["shared", "lin", "short", "long"], [
            GenomeInput((cu["shared"],), 11, kind="plasmid", circular=(cu["shared"],)),
            GenomeInput((cu["lin"],), 11, kind="plasmid", circular=()),
            GenomeInput((cu["short"],), 11, kind="plasmid", circular=(cu["short"],)),
            GenomeInput((cu["short"], cu["long"]), 11, kind="plasmid", circular=(cu["short"], cu["long"])),
        ]),
    ]:
        fna = in_dir / f"{st}.fna"
        fna.write_bytes(b"\n".join(ln for n in names for ln in fasta.format_record(n, raw[n], 80)) + b"\n")
        write_genome_input(root, st, genomes)
        ingest([SourceObject(fna, "b53_test", "1", "CC0", True)], root, st)
        runs = root / "runs" / f"stage={st}"
        runs.mkdir(parents=True, exist_ok=True)
        gi_files = sorted((root / "l0seq" / "genome_input" / f"stage={st}").glob("part-*.parquet"))
        (runs / "inputs.json").write_text(json.dumps({
            "outputs": [{"path": p.relative_to(root).as_posix()} for p in gi_files]
        }))
        ci_files = sorted((root / "l0seq" / "contig_index" / f"stage={st}").glob("part-*.parquet"))
        (runs / "ingest.json").write_text(json.dumps({
            "outputs": [{"path": p.relative_to(root).as_posix()} for p in ci_files]
        }))
    return root, raw


def test_b53_crispr_flank_bases_and_refusal_count_circular_wrap_flanks_of_stored_contigs(tmp_path):
    """`crispr_flank_bases` equals a brute-force sum of `2 * min(CIRC_FLANK, len)` over the distinct circular contigs
    the stage stores (skipping contigs owned by an earlier stage and counting a contig shared within the stage once),
    and `crispr_refusal` refuses a circular plasmid stage whose scanned bases pass `SHARD_MARGIN` even when its stored
    bases alone would fit."""
    root, raw = _stage_with_plasmids(tmp_path)
    stored_expected = len(raw["lin"]) + len(raw["short"]) + len(raw["long"])  # 320,000 bp
    flank_expected = 2 * min(crispr.CIRC_FLANK, len(raw["short"])) + 2 * min(crispr.CIRC_FLANK, len(raw["long"]))  # 280,000 bp
    assert bs.stage_bases(root, "PL") == stored_expected == 320_000

    prof = profile.load(Path(__file__).resolve().parents[2] / "profiles" / "wilkcomp02.yaml")
    # 0.8 * 0.0016 = 0.00128 CPU-min: 320 kb stored alone fits (~0.00090–0.00097), while 600 kb scanned (~0.00168) exceeds it
    tight = dataclasses.replace(prof, cpu_min_per_process=0.0016)
    b = bs.Build(root=root, stage="PL", release="rB53.1", genomes=tmp_path / "g.json", objects=tmp_path / "o.json",
                 stages=("S0", "PL"), crispr=True)
    step = bs.Step("PL", "crispr-minced", ("stage", "crispr"), check="crispr")
    why = bs.crispr_refusal(tight, b, root, step)
    assert why is not None and "600,000" in why and "280,000" in why and "no shards" in why, why
    assert bs.crispr_flank_bases(root, "PL") == flank_expected == 280_000


@pytest.mark.parametrize("j", range(len(PL1_BATCHES)))
def test_b53_each_genecall_batch_of_pl1_is_within_20_percent_of_its_estimate(j):
    lin_bp, circ_bp, cpu_min = PL1_BATCHES[j]
    est = genecall.one_process_cpu_min({("meta", False): lin_bp, ("meta", True): circ_bp})
    assert cpu_min <= est <= 1.20 * cpu_min, (j, round(est, 3), cpu_min)
    assert abs(cpu_min - est) <= 0.20 * est, (j, round(est, 3), cpu_min)


def test_b53_circular_gv_meta_rate_matches_the_v1_sample_and_holds_v1_within_20_percent():
    """Before card B-53, `CALL_S_PER_MBP[("gv-meta", True)]` was an estimate (6.3075 s per Mbp read = 4.2019 * 1.5011),
    giving 557.61 CPU-s against 386.87 CPU-s (+44.1%) on card C-39's 778 circular virus genomes of stage V1. With the
    measured rate (4.3762 s per Mbp read), both the sample's calls and V1's projected one-process cost sit at or within
    20% above the measurement."""
    est_circ_s = genecall.CALL_S_PER_MBP[("gv-meta", True)] * V1_CIRC_CALLED_BP / 1e6
    assert V1_CIRC_CPU_S <= est_circ_s <= 1.20 * V1_CIRC_CPU_S, (round(est_circ_s, 2), V1_CIRC_CPU_S)

    est_lin_s = genecall.CALL_S_PER_MBP[("gv-meta", False)] * V1_LIN_CALLED_BP / 1e6
    assert V1_LIN_CPU_S <= est_lin_s <= 1.20 * V1_LIN_CPU_S, (round(est_lin_s, 2), V1_LIN_CPU_S)

    # Whole-stage V1 projection from C-39's sample rates vs `one_process_cpu_min`
    sample_lin_rate = V1_LIN_CPU_S / (V1_LIN_CONTIG_BP / 1e6)
    sample_circ_contig_rate = V1_CIRC_CPU_S / (V1_CIRC_CONTIG_BP / 1e6)
    projected = (genecall.PROCESS_S + genecall.STEP_OVERHEAD * (
        sample_lin_rate * V1_ALL_LIN_CALLED_BP / 1e6 + sample_circ_contig_rate * V1_ALL_CIRC_CONTIG_BP / 1e6
    )) / 60
    est_v1 = genecall.one_process_cpu_min({("gv-meta", False): V1_ALL_LIN_CALLED_BP,
                                           ("gv-meta", True): V1_ALL_CIRC_CALLED_BP})
    assert projected <= est_v1 <= 1.20 * projected, (round(est_v1, 2), round(projected, 2))
    assert abs(projected - est_v1) <= 0.20 * est_v1, (round(est_v1, 2), round(projected, 2))
