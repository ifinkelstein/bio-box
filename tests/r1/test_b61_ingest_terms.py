"""Ring R1: card B-61's terms of the ingest estimate: bytes out of `.gz` objects, and bytes in records under the least
contig length (read and dropped), with every stage built from plain FASTA without a least length unchanged to the
byte (AGENTS §6; `planning/measurements/2026-11-21_ingest_terms.md`).

The oracles are measured, never the formulas of `stage`: the run records of the built stages (`runs/stage=<S>/
ingest.json` in `/tmp/if2453/r10/root`, read 2026-11-21) and the arms of cards B-58 and B-61 (each a `build-stage
--depth stored` in a scratch root on wilkcomp02, `cpu_min` of the ingest step's run record, the mean of three
interleaved runs; `/tmp/if2453/b58/results.jsonl` and `/tmp/if2453/b61/results.jsonl`, `scripts/stress/
b58_ingest_arms.py`), with what each read counted by that script's own loop over the decompressed lines. The
counter of the records' lengths is held against a loop over the lines written here and against `fasta.iter_records`.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from biobox import fasta
from biobox import stage as bs

pytestmark = pytest.mark.covers("B-61")

# --- every stage built from plain FASTA without a least length: the estimate is card B-52's, to the byte -----------

# stage: (bytes of FASTA read, bases of its listed contigs another stage held before it, card B-52's estimate as the
# code computed it before card B-61: `repr(stage.ingest_cpu_min(n, held))`, wilkcomp02, 2026-11-21, the tree of
# integration f60ed52 in the image). Bytes and held bases: tests/r1/test_b52_ingest_estimate.py (S1 to PL1) and
# planning/measurements/2026-11-20_stage_V1/estimate/before.json (V1). None of their objects is `.gz` and none was
# built with --min-contig-len (*measured*: the root's l0_manifest and run records, 2026-11-21).
PLAIN = {
    "S1": (2_121_124_779, 0, 8.787822841690001),
    "P1": (1_029_384_352, 0, 4.300769686720001),
    "P2": (9_716_599_273, 0, 40.00522301203001),
    "P3": (10_074_413_081, 6_963, 41.475809144980005),
    "P4": (10_101_497_437, 7_053, 41.587125478240004),
    "P5": (10_092_827_773, 185_683, 41.550758989900004),
    "P6": (10_087_447_605, 17_198, 41.529338972770006),
    "H1": (1_962_116_403, 141_683_767, 7.551978133960001),
    "PL1": (7_425_668_755, 1_177_838_866, 25.748580843790002),
    "V1": (5_228_187_536, 43_795_699, 21.377850450070003),
}


@pytest.mark.parametrize("stage", list(PLAIN))
def test_b61_a_plain_stage_without_a_least_length_keeps_its_estimate_to_the_byte(stage):
    n, held, before = PLAIN[stage]
    assert bs.ingest_cpu_min(n, held) == before
    assert bs.ingest_cpu_min(n, held, gz_bytes=0, below_bytes=0) == before


# --- C5: moved by the gzip term only -------------------------------------------------------------------------------

# C5 (planning/measurements/2026-11-20_stage_C5/before.json; runs/stage=C5/ingest.json): 2,032 `.fna.gz` objects,
# 9,509,853,070 bytes of FASTA, 267,273,180 bases held by earlier stages, no least length; record 45.803255 CPU-min.
C5 = (9_509_853_070, 267_273_180, 45.803255)
C5_BEFORE = 38.0570033479  # card B-52's estimate, computed as PLAIN's (DEC-56 accepted 45.803 against it)


def test_b61_c5_moves_by_the_gzip_term_alone_and_its_record_is_within_the_guard():
    n, held, rec = C5
    est = bs.ingest_cpu_min(n, held, gz_bytes=n)
    assert est - C5_BEFORE == pytest.approx(bs.INGEST_GZ_CPU_MIN * n / bs.INGEST_BP, rel=1e-12)
    assert bs.ingest_cpu_min(n, held) == C5_BEFORE
    # the gzip rate is fitted on the arms, never on C5: its ratio after is what the arms give, and under DEC-53's 1.20
    assert 1.0 < rec / est <= 1.20, round(rec / est, 4)


# --- every arm of cards B-58 and B-61 ----------------------------------------------------------------------------

# arm: (n bytes of FASTA, of them out of .gz, bases held, bytes in records under 5,000 bp, mean cpu_min of the ingest
# record over its runs). B-58's runs: /tmp/if2453/b58/results.jsonl (planning/measurements/2026-11-21_c5_ingest_cost.md
# §2: gz is all six runs); B-61's: /tmp/if2453/b61/results.jsonl. What each read: /tmp/if2453/b61/report.json
# (`b58_ingest_arms.py groups-inputs`, `c5-counts`, `g1-inputs`; `line_counts`).
C5S, C5X3, S58 = 944_466_060, 3_159_336_047, 1_297_869_659
ARMS = {
    "gz": (C5S, C5S, 0, 0, 4.390360),
    "plain": (C5S, 0, 0, 0, 4.144285),
    "few": (C5S, 0, 0, 0, 4.074595),
    "view": (C5S, C5S, 0, 0, 4.112413),
    "held A": (C5S, C5S, 0, 0, 4.238483),
    "held B": (C5S, C5S, 932_761_780, 0, 0.362013),
    "scale": (C5X3, C5X3, 0, 0, 14.672343),
    "spx": (S58, S58, 0, 0, 6.475237),
    "spx5000": (S58, S58, 0, 1_023_387_001, 1.719564),
    # card B-61's groups of SPX1 (2 files each); `<g> B` is the same objects again as stage B of the root (all held)
    "spxg1": (3_129_093, 3_129_093, 0, 0, 0.067698),
    "spxg1 B": (3_129_093, 3_129_093, 2_874_761, 0, 0.053197),
    "spxg1m": (3_129_093, 3_129_093, 0, 2_932_981, 0.054333),
    "spxg2": (134_563_454, 134_563_454, 0, 0, 0.699910),
    "spxg2 B": (134_563_454, 134_563_454, 128_631_463, 0, 0.119550),
    "spxg2m": (134_563_454, 134_563_454, 0, 105_088_104, 0.228299),
    "spxg3": (450_359_726, 450_359_726, 0, 0, 2.274661),
    "spxg3 B": (450_359_726, 450_359_726, 420_250_690, 0, 0.360313),
    "spxg3m": (450_359_726, 450_359_726, 0, 427_508_895, 0.348061),
    "spxg4": (963_086_613, 963_086_613, 0, 0, 4.709720),
    "spxg4 B": (963_086_613, 963_086_613, 912_348_458, 0, 0.549665),
    "spxg4m": (963_086_613, 963_086_613, 0, 841_176_755, 0.907762),  # six runs: three more beside g1's
    "spxg5": (2_124_973_510, 2_124_973_510, 0, 0, 10.167476),
    "spxg5 B": (2_124_973_510, 2_124_973_510, 2_035_093_705, 0, 1.050529),
    "spxg5m": (2_124_973_510, 2_124_973_510, 0, 1_366_543_629, 3.763481),
    "spxg6": (5_563_406_950, 5_563_406_950, 0, 0, 27.678208),
    "spxg6 B": (5_563_406_950, 5_563_406_950, 5_267_907_736, 0, 3.374714),
    "spxg6m": (5_563_406_950, 5_563_406_950, 0, 4_901_131_525, 4.860855),
    # every 10th object of G1 (714 GlobDB MAGs as .fa.gz, 121,666 contigs), every contig kept
    "g1": (996_211_156, 996_211_156, 0, 0, 4.239964),
}
# Not judged at 15% (planning/measurements/2026-11-21_ingest_terms.md §3): spxg1's three arms read 3 MB and are the
# start-up of a process (0.05 to 0.07 CPU-min; DEC-43: a step that small has no rate), held to 0.025 CPU-min instead;
# spxg3 B and spxg6 B hold 0.7 and 6.9 million short contigs (633 and 806 bytes of FASTA each on average) and pay a cost per
# record the estimate cannot know before ingest without counting records (about 4 to 6 microseconds each, *derived*),
# held to 25% and named here.
START_UP = ("spxg1", "spxg1 B", "spxg1m")
SHORT_HELD = ("spxg3 B", "spxg6 B")


def _est(arm: str) -> float:
    n, z, h, d, _ = ARMS[arm]
    return bs.ingest_cpu_min(n, h, gz_bytes=z, below_bytes=d)


@pytest.mark.parametrize("arm", [a for a in ARMS if a not in START_UP + SHORT_HELD])
def test_b61_each_arm_is_within_15_percent_of_its_estimate(arm):
    rec = ARMS[arm][4]
    assert abs(rec - _est(arm)) <= 0.15 * _est(arm), (arm, round(_est(arm), 4), rec, f"{rec / _est(arm):.3f}")


@pytest.mark.parametrize("arm", START_UP)
def test_b61_an_arm_of_3_mb_is_its_start_up_within_a_fortieth_of_a_minute(arm):
    assert abs(ARMS[arm][4] - _est(arm)) <= 0.025, (arm, round(_est(arm), 4), ARMS[arm][4])


@pytest.mark.parametrize("arm", SHORT_HELD)
def test_b61_millions_of_short_held_contigs_stay_within_25_percent(arm):
    rec = ARMS[arm][4]
    assert 1.0 < rec / _est(arm) <= 1.25, (arm, round(_est(arm), 4), rec, f"{rec / _est(arm):.3f}")


# G1 (planning/measurements/2026-11-18_globdb_inputs/G1/est.json, `ingest`): 7,139 `.fa.gz` objects, 10,087,932,724
# bytes of FASTA, 5,699,697 bases held when its inputs were staged, no least length; card B-52's estimate 41.508.
G1 = (10_087_932_724, 5_699_697)


def test_b61_g1_predicted_from_its_own_sample_is_within_15_percent_of_the_new_estimate():
    """*Derived*: the g1 arm's cost per byte (every 10th object, every contig kept; the mean of three runs) carried to
    the whole stage. Card B-52's estimate stays where it was (41.508); the new one adds the gzip term."""
    n, held = G1
    old = bs.ingest_cpu_min(n, held)
    new = bs.ingest_cpu_min(n, held, gz_bytes=n)
    assert round(old, 3) == 41.508
    sn, _, _, _, rec = ARMS["g1"]
    pred = bs.INGEST_FIXED_CPU_MIN + (rec - bs.INGEST_FIXED_CPU_MIN) / sn * (n - held)
    assert abs(pred - new) <= 0.15 * new, (round(pred, 2), round(new, 2), round(pred / new, 3))
    assert pred / old < 1.20  # nor would it have stopped under the old estimate, at the arms' load


# SPX1, counted by the new code read-only on its staged inputs (`b58_ingest_arms.py spx1-plan`: `object_sizes`,
# `below_bytes` at --min-contig-len 5000, `gz_bytes`; wilkcomp02, 2026-11-21, /tmp/if2453/b61/plan.json): 304 `.fa.gz`
# objects, 75,267,693,104 bytes on disk (the largest 4,938,946,579), 256,948,810,775 bytes of FASTA, 220,918,826,783 of
# them in records under 5,000 bp; no contig held (card C-48e: 0 shared with the root).
SPX1 = dict(on_disk=75_267_693_104, largest=4_938_946_579, fasta=256_948_810_775, below=220_918_826_783)
B58_SPX1 = 264.0  # card B-58's derivation from its two spx arms (planning/measurements/2026-11-21_c5_ingest_cost.md §4)


def test_b61_spx1_is_estimated_near_card_b58_s_figure_and_its_first_pass_is_two_shards():
    s = SPX1
    est = bs.ingest_cpu_min(s["fasta"], 0, gz_bytes=s["fasta"], below_bytes=s["below"])
    assert round(bs.ingest_cpu_min(s["fasta"]), 1) == 1056.1  # card B-52's: every byte stored
    assert abs(est - B58_SPX1) <= 0.15 * B58_SPX1, round(est, 1)
    assert round(est, 1) == 245.9
    prof = type("P", (), {"cpu_min_per_process": 300})()
    per_shard = bs._pass1_bytes(prof, s["on_disk"], s["fasta"], s["fasta"], s["below"])
    before = bs._pass1_bytes(prof, s["on_disk"], s["fasta"])
    assert -(-s["on_disk"] // before) >= 5  # card B-52's cut: five shards or more
    # maximal runs of at most per_shard bytes: the first holds more than per_shard less the largest object, so the
    # rest fits one more shard when this holds
    assert s["on_disk"] - (per_shard - s["largest"]) <= per_shard
    assert per_shard == 73_464_662_516


def test_b61_the_fitted_rates_are_the_note_s():
    """The constants and where they came from: the joint fit of the note (§3) gave 0.258 and 0.131 per Gbyte (0.252 to
    0.259 and 0.113 to 0.138 with one arm left out); the rate under the least length is taken at the top, rounded up."""
    assert bs.INGEST_GZ_CPU_MIN == 0.26 and bs.INGEST_BELOW_CPU_MIN == 0.14


# --- the counter of the records under the least length ---------------------------------------------------------------


def _oracle(data: bytes, min_len: int) -> tuple[int, int]:
    """(bytes, bytes in records of fewer than min_len bases) by a loop over the lines."""
    below, cur = 0, None
    for ln in data.splitlines(keepends=True):
        if ln.startswith(b">"):
            if cur is not None and cur[1] < min_len:
                below += cur[0]
            cur = [len(ln), 0]
        elif cur is not None:
            cur[0] += len(ln)
            cur[1] += len(ln.rstrip(b"\r\n"))
    if cur is not None and cur[1] < min_len:
        below += cur[0]
    return len(data), below


CASES = {
    "one-line records": b">a x\nACGTACGTAC\n>b\nACG\n>c long header with > inside\nACGTACGTACGTACGTACGT\n",
    "wrapped lines": b">a\nACGTA\nCGTAC\nGT\n>b\nAC\nGT\n>c\nACGTACGTAC\nACGTACGTAC\nA\n",
    "CRLF and blank lines": b">a\r\nACGTACG\r\n\r\nTAC\r\n>b\r\nACG\r\n\n>c\r\n\r\n",
    "no final newline": b">a\nACGTACGTACGT\n>b\nACGTA",
    "text before the first header": b"junk\n>a\nACGTACGTACGT\n>b\nAC\n",
    "empty record last": b">a\nACGTACGTACGT\n>b\n",
    "empty": b"",
}


@pytest.mark.parametrize("case", list(CASES))
@pytest.mark.parametrize("read", [1, 2, 3, 5, 7, 64, 1 << 20])
@pytest.mark.parametrize("gz", [False, True], ids=["plain", "gz"])
def test_b61_the_counter_finds_the_bytes_in_records_under_the_least_length(tmp_path, monkeypatch, case, read, gz):
    """Every chunk size, so that headers, `\\n>` and records cross the edges of the reads."""
    data = CASES[case]
    p = tmp_path / ("x.fa.gz" if gz else "x.fa")
    p.write_bytes(gzip.compress(data, mtime=0) if gz else data)
    monkeypatch.setattr(bs, "_READ", read)
    for min_len in (1, 4, 10, 11, 12, 13, 21):
        assert bs.fasta_counts(p, min_len) == _oracle(data, min_len), min_len
    assert bs.fasta_counts(p, 5)[0] == bs.fasta_bytes(p)


def test_b61_the_counter_drops_what_ingest_drops(tmp_path, monkeypatch):
    """The records the counter puts under the least length are those `fasta.iter_records` reads as shorter."""
    monkeypatch.setattr(bs, "_READ", 3)
    data = b"".join(CASES[c] for c in ("one-line records", "wrapped lines", "CRLF and blank lines"))
    p = tmp_path / "all.fa"
    p.write_bytes(data)
    for min_len in range(1, 25):
        lens = [len(r.seq) for r in fasta.iter_records(p)]
        recs = data.split(b"\n>")
        sizes = [len(r) + (0 if i == len(recs) - 1 else 1) + (1 if i else 0) for i, r in enumerate(recs)]
        assert len(sizes) == len(lens)
        want = sum(s for s, n in zip(sizes, lens) if n < min_len)
        assert bs.fasta_counts(p, min_len) == (len(data), want), min_len


def test_b61_a_line_of_spaces_counts_as_bases_so_the_estimate_errs_high(tmp_path):
    """`fasta.iter_records` drops a line of only spaces; the counter counts them, so the record looks longer than
    ingest finds it (kept, never dropped where ingest keeps it)."""
    p = tmp_path / "s.fa"
    p.write_bytes(b">a\nACG\n      \n")
    assert [len(r.seq) for r in fasta.iter_records(p)] == [3]
    assert bs.fasta_counts(p, 5) == (p.stat().st_size, 0)  # 9 "bases": kept by the counter, dropped by ingest
    assert bs.fasta_counts(p, 10) == (p.stat().st_size, p.stat().st_size)


def test_b61_object_sizes_reads_once_and_below_bytes_follows_the_least_length(tmp_path, monkeypatch):
    data = CASES["wrapped lines"]
    paths = []
    for i, gz in enumerate((True, False)):
        p = tmp_path / (f"o{i}.fa.gz" if gz else f"o{i}.fa")
        p.write_bytes(gzip.compress(data, mtime=0) if gz else data)
        paths.append(p)
    objs = tmp_path / "objects.json"
    objs.write_text("[" + ",".join(f'{{"path": "{p}"}}' for p in paths) + "]")
    calls = []
    real = bs.fasta_counts
    monkeypatch.setattr(bs, "fasta_counts", lambda p, n: calls.append(p) or real(p, n))
    monkeypatch.setattr(bs, "_SIZES", {})
    monkeypatch.setattr(bs, "_BELOW", {})
    sizes = bs.object_sizes(objs, 11)
    assert [(p, f) for p, _, f in sizes] == [(p, len(data)) for p in paths]
    assert bs.gz_bytes(sizes) == len(data)
    assert bs.below_bytes(objs, 11) == 2 * _oracle(data, 11)[1]
    assert bs.below_bytes(objs, None) == 0
    assert len(calls) == 2  # one read of each object for the sizes and the counts together
