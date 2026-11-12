"""Ring R1: a clustering call that cannot be split under the budget runs as one process (card B-45, DEC-29).

Without MMseqs2: the cut is made from hand-written `.index` files (as `test_cluster_call_estimates` does), and the
gate of `build-stage` (`stage.cluster_refusal`) is a pure function of the profile and the estimates.

- A prefilter whose parts each rebuild the target's k-mer table (no shared index with `--diag-score 0`, card B-44),
  when that table alone fills the part budget: with the plan's one-process cap (the profile's `cluster_cpu_min_max`,
  recorded in the plan's `cut`) it is one part, its record saying why and its estimate (work plus table); above the cap
  it is refused; without the cap it is refused as before card B-45, the text naming DEC-29.
- Every other call is cut as before: the cap changes no other cut.
- The overhead can be forced by an explicit plan parameter (`rebuild_cpu_min`, recorded in `cut`), and the real
  rates reach the same case when the part budget is below the table's estimate.
- `release cluster plan` shows which calls may run as one process and why (`cluster.one_process_outlook`).
- The release-level gate: clustering in parts is checked call by call (each call the plan runs as one process
  against the cap), not by the one-process estimate of the whole release; in one process, as before.
"""

from __future__ import annotations

import dataclasses
import math
from pathlib import Path

import pytest

from biobox import profile
from biobox import stage as bs
from biobox.l2 import cluster

pytestmark = pytest.mark.covers("MR-16", "B-36")  # ring R1 (by directory)

PART = 240.0  # 80% of the 300 CPU-minute budget, build-stage's part budget
FAST = ("cluster.nr90.11.prefilter", "cluster.nr50.11.prefilter")  # their parts rebuild the table (B-44)
SENSITIVE = ("cluster.nr50.15.prefilter", "cluster.nr50.19.prefilter")  # their parts share an index (B-43)
BY = {c.name: c for c in cluster.calls()}


def cut_doc(tmp_path: Path, name: str, units: int, keys: int, plan_cut: dict, proteins: int = 10_620_684,
            shared_index: bool = True) -> dict:
    """The cut record the code makes for call `name` whose split input holds `keys` entries and `units` residues or
    entry bytes in all (for a prefilter its query database is its target, so `keys` targets)."""
    c = BY[name]
    prefilter = c.module == "prefilter"
    d = tmp_path / name
    src = d / f"{c.positional[cluster.PER_QUERY[c.module]]}.index"
    src.parent.mkdir(parents=True, exist_ok=True)
    pad = 2 if prefilter else 0  # a sequence entry is its residues, a newline and a NUL
    src.write_text("".join(f"{k}\t{k * 10}\t{units // keys + pad}\n" for k in range(keys)))
    ctx = cluster._Ctx(root=d, out=None, release="rX", plan={"proteins": proteins, "cut": plan_cut}, plan_sha256="",
                       w=d, threads=1, shared_index=shared_index)
    return cluster._cut_doc(ctx, c)


def test_b45_a_call_whose_rebuilt_table_fills_the_part_budget_is_one_part_under_the_cap(tmp_path):
    units, keys = 400_000_000, 1_000
    for name in FAST:
        got = cut_doc(tmp_path, name, units, keys,
                      {"part_cpu_min": PART, "rebuild_cpu_min": 300.0, "one_process_cpu_min_max": 20000.0})
        work = cluster.call_cpu_min(name, 10_620_684, units)
        assert got["cpu_min"] == pytest.approx(work, abs=1e-5)
        assert got["overhead_cpu_min"] == pytest.approx(300.0)
        assert len(got["parts"]) == 1 and got["parts"][0]["keys"] == keys
        assert got["parts"][0]["cpu_min"] == pytest.approx(work + 300.0, abs=1e-4)
        one = got["one_process"]
        assert one["decision"] == "DEC-29" and one["cpu_min_max"] == 20000.0
        assert one["cpu_min"] == pytest.approx(work + 300.0, abs=1e-4)
        assert "rebuild" in one["why"] and "B-44" in one["why"]
        assert "index_cpu_min" not in got and "index" not in got


def test_b45_above_the_cap_the_one_process_call_is_refused(tmp_path):
    units = 400_000_000
    work = cluster.call_cpu_min(FAST[0], 10_620_684, units)
    with pytest.raises(ValueError, match=r"DEC-29.*over the one-process cap of 300"):
        cut_doc(tmp_path, FAST[0], units, 1_000,
                {"part_cpu_min": PART, "rebuild_cpu_min": 300.0, "one_process_cpu_min_max": 300.0})
    # at the cap exactly it runs
    got = cut_doc(tmp_path, FAST[0], units, 1_000,
                  {"part_cpu_min": PART, "rebuild_cpu_min": 300.0, "one_process_cpu_min_max": work + 300.0})
    assert got["one_process"]["decision"] == "DEC-29"


def test_b45_without_the_cap_the_refusal_is_todays_and_names_dec29(tmp_path):
    for name in FAST:
        with pytest.raises(ValueError) as e:
            cut_doc(tmp_path, name, 400_000_000, 1_000, {"part_cpu_min": PART, "rebuild_cpu_min": 300.0})
        text = str(e.value)
        assert "alone fills the part budget of 240.000: no number of parts fits it" in text  # today's words
        assert "DEC-29" in text and "cluster_cpu_min_max" in text and name in text


def test_b45_the_real_table_rate_reaches_the_same_case(tmp_path):
    """No forced overhead: a part budget below the real table estimate of 1,000 targets."""
    keys = 1_000
    table = cluster.part_overhead_cpu_min(BY[FAST[0]], keys, shared=False)
    got = cut_doc(tmp_path, FAST[0], 1_000_000, keys, {"part_cpu_min": table * 0.9, "one_process_cpu_min_max": 1.0})
    assert len(got["parts"]) == 1 and got["one_process"]["decision"] == "DEC-29"
    assert got["overhead_cpu_min"] == pytest.approx(table, abs=1e-6)
    with pytest.raises(ValueError, match="DEC-29"):
        cut_doc(tmp_path, FAST[0], 1_000_000, keys, {"part_cpu_min": table * 0.9})
    # just above the table, the call is cut as before, into parts, without the cap mattering
    a = cut_doc(tmp_path, FAST[0], 1_000_000, keys, {"part_cpu_min": table * 1.1, "one_process_cpu_min_max": 1.0})
    b = cut_doc(tmp_path, FAST[0], 1_000_000, keys, {"part_cpu_min": table * 1.1})
    assert a == b and "one_process" not in a


@pytest.mark.parametrize("name", sorted(c.name for c in cluster.calls() if c.split))
def test_b45_every_other_call_is_cut_as_before(tmp_path, name):
    """The cap in the plan changes no cut that is not DEC-29's case: realistic sizes at 10.6 M proteins."""
    c = BY[name]
    keys = 2_000 if c.module == "prefilter" else 4
    units = 3_000_000_000
    with_cap = cut_doc(tmp_path / "a", name, units, keys, {"part_cpu_min": PART, "one_process_cpu_min_max": 20000.0})
    without = cut_doc(tmp_path / "b", name, units, keys, {"part_cpu_min": PART})
    assert with_cap == without and "one_process" not in with_cap
    forced = cut_doc(tmp_path / "c", name, units, keys, {"parts": 3, "one_process_cpu_min_max": 20000.0})
    assert len(forced["parts"]) == 3 and "one_process" not in forced  # a forced count is never DEC-29's case


def test_b45_a_call_that_shares_an_index_is_not_dec29s_case(tmp_path):
    """The sensitive prefilters share an index: their parts load it, so a forced rebuild only stops them from being
    one part (they are two, sharing the index); with the index off (tests only) their parts rebuild, and the case is
    DEC-29's."""
    for name in SENSITIVE:
        got = cut_doc(tmp_path / name, name, 100_000_000, 1_000, {"part_cpu_min": PART, "rebuild_cpu_min": 300.0})
        assert len(got["parts"]) == 2 and "index_cpu_min" in got and "one_process" not in got
        off = cut_doc(tmp_path / f"{name}-off", name, 100_000_000, 1_000,
                      {"part_cpu_min": PART, "rebuild_cpu_min": 300.0, "one_process_cpu_min_max": 20000.0},
                      shared_index=False)
        assert len(off["parts"]) == 1 and off["one_process"]["decision"] == "DEC-29"


def test_b45_the_plan_shows_which_calls_may_run_as_one_process_and_why():
    rate = cluster.PREFILTER_INDEX_CPU_S_PER_TARGET
    from_targets = math.ceil(PART * 60 / rate)  # the table of a part alone is 240 CPU-minutes: about 68.3 M targets
    assert 68_000_000 < from_targets < 68_500_000
    small = cluster.one_process_outlook(10_620_684, {"part_cpu_min": PART, "one_process_cpu_min_max": 20000.0})
    assert sorted(small["calls"]) == sorted(FAST)
    assert all(v == {"from_targets": from_targets, "possible": False} for v in small["calls"].values())
    assert small["cpu_min_max"] == 20000.0 and "B-44" in small["why"] and "DEC-29" in small["why"]
    big = cluster.one_process_outlook(95_000_000, {"part_cpu_min": PART})
    assert all(v["possible"] for v in big["calls"].values()) and big["cpu_min_max"] is None
    assert "refused" in big["why"]  # without the cap such a call is refused
    forced = cluster.one_process_outlook(1_000, {"part_cpu_min": PART, "rebuild_cpu_min": 300.0})
    assert all(v == {"from_targets": 0, "possible": True} for v in forced["calls"].values())
    assert cluster.one_process_outlook(95_000_000, {"parts": 4})["calls"] == {}  # a forced count never is


@pytest.mark.parametrize("kw", [{"one_process_cpu_min_max": -1.0}, {"one_process_cpu_min_max": float("nan")},
                                {"rebuild_cpu_min": -0.5}, {"rebuild_cpu_min": float("inf")},
                                {"one_process_cpu_min_max": True}])
def test_b45_the_plans_new_parameters_are_checked_before_any_work(tmp_path, kw):
    with pytest.raises(ValueError, match=next(iter(kw))):
        cluster.plan(tmp_path / "root", "rX", "gv", ["S0"], tmp_path / "tmp", part_cpu_min=PART, **kw)
    assert not (tmp_path / "tmp").exists()


# --- the release-level gate (stage.cluster_refusal) ---------------------------------------------------------------


def _profile(**kw) -> profile.Profile:
    p = profile.load(Path(__file__).resolve().parents[2] / "profiles" / "wilkcomp02.yaml")
    return dataclasses.replace(p, **kw)


def test_b45_in_parts_the_gate_checks_the_one_process_calls_not_the_whole_release():
    n = 95_000_000
    prof = _profile(cluster_cpu_min_max=20000)
    assert bs.cluster_cpu_min(n) > 20000
    # in one process (the unsplit mode): today's check, refused
    why = bs.cluster_refusal(prof, "rX", n)
    assert why is not None and "cluster_cpu_min_max" in why and "DEC-25" in why
    # in parts: no call yet planned as one process, or each under the cap: not refused
    assert bs.cluster_refusal(prof, "rX", n, one_process=[]) is None
    assert bs.cluster_refusal(prof, "rX", n, one_process=[("cluster.nr90.11.prefilter", 3945.4)]) is None
    # a call over the cap is refused, by name
    why = bs.cluster_refusal(prof, "rX", n, one_process=[("cluster.nr90.11.prefilter", 3945.4),
                                                         ("cluster.nr50.11.prefilter", 20000.5)])
    assert why is not None and "cluster.nr50.11.prefilter" in why and "DEC-29" in why
    assert "cluster.nr90.11.prefilter" not in why and "cluster_cpu_min_max of 20000" in why


def test_b45_without_the_cap_the_gate_is_unchanged_in_either_mode():
    prof = _profile(cluster_cpu_min_max=None)
    n = 95_000_000
    for one in (None, [], [("cluster.nr90.11.prefilter", 1.0)]):
        why = bs.cluster_refusal(prof, "rX", n, one_process=one)
        assert why is not None and "no cluster_cpu_min_max" in why and "DEC-25" in why
    assert bs.cluster_refusal(prof, "rX", 1_000, one_process=[]) is None


def test_b45_a_call_is_checked_against_the_profile_before_its_launch():
    prof = _profile(cluster_cpu_min_max=20000)
    assert bs.one_process_refusal(prof, "rX", "cluster.nr90.11.prefilter", 3945.4) is None
    why = bs.one_process_refusal(prof, "rX", "cluster.nr90.11.prefilter", 20001.0)
    assert why is not None and "DEC-29" in why and "cluster_cpu_min_max of 20000" in why
    why = bs.one_process_refusal(_profile(cluster_cpu_min_max=None), "rX", "cluster.nr90.11.prefilter", 1.0)
    assert why is not None and "no cluster_cpu_min_max" in why and "DEC-29" in why
    assert bs.DEC29_EXEMPT == "exempt: DEC-29"
