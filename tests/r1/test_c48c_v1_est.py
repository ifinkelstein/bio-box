"""Ring R1: the pure functions of card C-48c's stage-V1 driver estimates (`scripts/stress/build_stage_pod/c48c_v1_est.py`).

Written first (AGENTS §5.1). Each function is checked against values written out here by hand or against a
brute-force loop written separately in this file, never against the function itself:

- `judge_rate`: DEC-43 (a step is judged only when estimated at 24 CPU-minutes or more) and DEC-53 (a judged step
  stops the build only when more than 20% over its estimate; more than 20% under is recorded, not stopped);
- `rates_exit`: the exit code the driver stops on;
- `owners_problem`: the launch gate on contigs of V1 that stage C5 (parked, half built) already owns;
- `band` and `expected_ranges`: the ranges the check holds the built stage's counts to, from card C-39's note
  (`planning/measurements/2026-11-16_virus_votu_inputs.md`), exact where the inputs fix a count, +-20% where the
  note's figure is an estimate (AGENTS §6);
- `check_counts`: which observed counts fall outside their range;
- `compare_c39`: the pre-launch estimates against the note's §5, flagging a difference over 20% that the code does
  not explain;
- `pick_sample`: the rehearsal's sample of genomes (deterministic, holds the shared and circular genomes asked for).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-48c")  # ring R1 by directory (tests/conftest.py)

ROOT = Path(__file__).resolve().parents[2]
EST_SCRIPT = ROOT / "scripts" / "stress" / "build_stage_pod" / "c48c_v1_est.py"


def _load():
    spec = importlib.util.spec_from_file_location("c48c_v1_est", EST_SCRIPT)
    assert spec is not None and spec.loader is not None, f"missing script: {EST_SCRIPT}"
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c48c_v1_est"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def est():
    return _load()


# DEC-43 and DEC-53 -----------------------------------------------------------------------------------------------


def test_small_step_is_reported_not_judged(est):
    row, over, under = est.judge_rate("crispr-x", 20.0, 10.0, "test")  # twice its estimate, but estimated under 24
    assert (row["judged"], over, under) == (False, False, False)
    assert row["ratio"] == 2.0


def test_judged_from_24_cpu_minutes(est):
    row, over, under = est.judge_rate("s", 30.0, 24.0, "test")  # 1.25 on a step estimated at exactly 24
    assert (row["judged"], over, under) == (True, True, False)
    row, over, under = est.judge_rate("s", 30.0, 23.999, "test")
    assert (row["judged"], over, under) == (False, False, False)


def test_over_20_percent_stops_under_is_recorded(est):
    _, over, under = est.judge_rate("domains-pfam-38.2.shard-00003", 290.0, 240.0, "test")  # 1.208
    assert (over, under) == (True, False)
    _, over, under = est.judge_rate("genecall.batch-00001", 150.0, 200.0, "test")  # 0.75: recorded, no stop
    assert (over, under) == (False, True)
    _, over, under = est.judge_rate("genecall.batch-00002", 192.0, 200.0, "test")  # 0.96
    assert (over, under) == (False, False)


def test_bounds_are_strict(est):
    _, over, under = est.judge_rate("s", 120.0, 100.0, "test")  # exactly 1.20: not over
    assert (over, under) == (False, False)
    _, over, under = est.judge_rate("s", 80.0, 100.0, "test")  # exactly 0.80: not under
    assert (over, under) == (False, False)
    _, over, _ = est.judge_rate("s", 120.01, 100.0, "test")
    assert over


def test_rates_exit(est):
    assert est.rates_exit([]) == 0
    assert est.rates_exit(["genecall.batch-00000"]) == 3
    assert est.rates_exit(["domains-pfam-38.2.shard-00001", "ingest"]) == 3


# the C5 gate -----------------------------------------------------------------------------------------------------


def test_owners_problem(est):
    assert est.owners_problem({"H1": 177, "S1": 1}) is None
    assert est.owners_problem({"H1": 177, "S1": 1, "C5": 0}) is None
    msg = est.owners_problem({"H1": 177, "C5": 1})
    assert msg is not None and "C5" in msg and "1" in msg
    assert est.owners_problem({"X9": 2}, forbidden=("X9",)) is not None


# expected ranges --------------------------------------------------------------------------------------------------


def test_band(est):
    assert est.band(100, 0.2) == (80, 120)
    assert est.band(7_390_000, 0.2) == (5_912_000, 8_868_000)
    assert est.band(5, 0.0) == (5, 5)


def test_expected_ranges_from_the_c39_note(est):
    r = est.expected_ranges()
    # fixed by the inputs: values copied by hand from the C-39 note (§3, §4), all measured there
    exact = {
        "genomes": 168_536,
        "genome_input_rows": 168_536,
        "contigs": 168_536,
        "stored": 168_358,
        "stored_by_earlier": 178,
        "ledger_skipped_earlier": 178,
        "bases_stored": 5_178_493_077,
        "circular": 22_517,
    }
    for k, v in exact.items():
        assert r[k] == (v, v), k
    # estimates of the note's §5, +-20% (AGENTS §6)
    assert r["genes"] == (5_912_000, 8_868_000)  # about 7.39 M genes
    assert r["distinct_proteins"] == (round(5_597_468 * 0.8), round(5_597_468 * 1.2))
    assert r["residues"] == (round(1_307_249_946 * 0.8), round(1_307_249_946 * 1.2))


def _brute_problems(obs: dict, exp: dict) -> list[str]:
    out = []
    for k in sorted(exp):
        lo, hi = exp[k]
        if k not in obs:
            out.append(k)
        elif not (lo <= obs[k] <= hi):
            out.append(k)
    return out


def test_check_counts(est):
    exp = {"a": (10, 10), "b": (80, 120), "c": (0, 5)}
    cases = [
        {"a": 10, "b": 80, "c": 5},
        {"a": 10, "b": 120, "c": 0},
        {"a": 11, "b": 100, "c": 3},
        {"a": 10, "b": 79, "c": 3},
        {"a": 10, "b": 121, "c": 6},
        {"b": 100, "c": 3},
    ]
    for obs in cases:
        got = est.check_counts(obs, exp)
        assert sorted(p.split(":")[0] for p in got) == _brute_problems(obs, exp), obs


# against the C-39 note --------------------------------------------------------------------------------------------


def test_compare_c39_flags_only_unexplained(est):
    c39 = {"ingest": 23.44, "genecall": 687.5, "crispr_default": 15.71}
    ours = {"ingest": 21.38, "genecall": 500.0, "crispr_default": 30.0}  # 0.912, 0.727, 1.910
    rows = est.compare_c39(ours, c39, explained={"genecall": "a reason"})
    by = {r["name"]: r for r in rows}
    assert by["ingest"]["off_20pct"] is False and by["ingest"]["unexplained"] is False
    assert by["genecall"]["off_20pct"] is True and by["genecall"]["unexplained"] is False
    assert by["crispr_default"]["off_20pct"] is True and by["crispr_default"]["unexplained"] is True
    assert round(by["genecall"]["ratio"], 3) == 0.727


# the rehearsal's sample ---------------------------------------------------------------------------------------------


def _brute_sample(n, circular, shared, every, n_circ):
    pick = set(range(0, n, every))
    extra = [i for i in sorted(circular) if i not in pick][:n_circ]
    return sorted(pick | set(extra) | set(shared))


def test_pick_sample(est):
    n = 1000
    circular = [3, 17, 200, 201, 500, 999]
    shared = [42, 900]
    got = est.pick_sample(n, circular, shared, every=100, n_circ=3)
    assert got == _brute_sample(n, circular, shared, 100, 3)
    assert got == sorted(set(got))
    assert {42, 900} <= set(got)
    assert est.pick_sample(n, list(reversed(circular)), list(reversed(shared)), every=100, n_circ=3) == got
