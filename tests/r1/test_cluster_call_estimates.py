"""Ring R1: the per-call CPU estimates of clustering in parts against the measured run records (card C-22).

`c22_cluster_calls_measured.tsv` is the measured sample (wilkcomp02, 2026-11-11; its header says which files): every
per-query call, and every part of a call cut into parts, of the clustering of release `rN.3` (10,620,684 proteins) and,
as a second sample at a fifth of the size, of `rC.1b` (2,053,075), with the CPU-minutes the run record measured, the
size of the call's own input (the queries and the residues or entry bytes the cut splits) and, for a prefilter, its
targets. The estimate is read from the code (`cluster.call_cpu_min`, `cluster.part_overhead_cpu_min`) the way the cut
uses it: a part's share of its call's work by units, plus what every part pays.

- Every `rN.3` part is within 20% of its record (AGENTS §6); the model was fitted on that sample.
- The second sample is out of that fit: every call that took 2 CPU-minutes or more is within 20% there too.
- The whole-release estimate of card B-36 (`stage.cluster_cpu_min`, the DEC-25 gate) and the sum of the per-call
  estimates plus the measured one-process rest agree within 20% at both measured sizes.
- The part counts that follow from the estimates keep every `rN.3` call's measured CPU-minutes per part under the
  300 CPU-minute budget of AGENTS §7.2 (the old constants left one call over it).
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import pytest

from biobox import stage
from biobox.l2 import cluster

pytestmark = pytest.mark.covers("MR-16")  # ring R1 (by directory)

TSV = Path(__file__).with_name("c22_cluster_calls_measured.tsv")
BUDGET, PART_BUDGET = 300.0, 240.0  # AGENTS §7.2 and 80% of it (build-stage's cut)


def sample() -> list[dict]:
    with open(TSV) as fh:
        rows = list(csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"))
    return [r for r in rows if r["step"] != "one_process_rest"]


def rest(release: str) -> float:
    with open(TSV) as fh:
        rows = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
        return next(float(r["cpu_min"]) for r in rows if r["step"] == "one_process_rest" and r["sample"] == release)


def estimates(release: str) -> list[tuple[dict, float, float]]:
    """(row, the part's estimate, its measured CPU-minutes) for every row of `release`: the call's work from its
    total units, shared out by the part's units, plus the per-part overhead as the cut adds it (the table's rebuild
    for a prefilter: both samples ran before card B-43's shared index)."""
    by = {c.name: c for c in cluster.calls()}
    groups = defaultdict(list)
    for r in sample():
        if r["sample"] == release:
            groups[r["step"]].append(r)
    out = []
    for name, parts in groups.items():
        total = sum(int(p["units"]) for p in parts)
        work = cluster.call_cpu_min(name, int(parts[0]["proteins"]), total)
        for p in parts:
            oh = cluster.part_overhead_cpu_min(by[name], int(p["targets"]), shared=False) if p["targets"] else 0.0
            out.append((p, work * int(p["units"]) / total + oh, float(p["cpu_min"])))
    return out


def test_c22_the_sample_is_the_two_releases_every_per_query_call():
    names = {c.name for c in cluster.calls() if c.split}
    for release, proteins in (("rN.3", 10_620_684), ("rC.1b", 2_053_075)):
        rows = [r for r in sample() if r["sample"] == release]
        assert {r["step"] for r in rows} == names
        assert {int(r["proteins"]) for r in rows} == {proteins}
    assert [r["parts"] for r in sample() if r["sample"] == "rN.3" and r["step"] == "cluster.nr90.12.align"] == ["3"] * 3


def test_c22_each_rN3_part_is_estimated_within_20_percent_of_its_record():
    bad = [(p["step"], p["part"], round(est, 2), meas, f"{est / meas - 1:+.0%}")
           for p, est, meas in estimates("rN.3") if abs(est / meas - 1) > 0.20]
    assert not bad, bad


def test_c22_the_second_sample_is_within_20_percent_for_every_call_of_2_cpu_minutes_or_more():
    bad = [(p["step"], round(est, 2), meas, f"{est / meas - 1:+.0%}")
           for p, est, meas in estimates("rC.1b") if meas >= 2.0 and abs(est / meas - 1) > 0.20]
    assert not bad, bad


@pytest.mark.parametrize("release", ["rN.3", "rC.1b"])
def test_c22_the_whole_release_estimate_agrees_with_the_sum_of_the_calls(release):
    n = int(next(r["proteins"] for r in sample() if r["sample"] == release))
    calls = sum(est for _, est, _ in estimates(release)) + rest(release)
    whole = stage.cluster_cpu_min(n)
    assert abs(calls / whole - 1) <= 0.20, (release, round(calls, 1), round(whole, 1))


def test_c22_the_cut_the_estimates_give_keeps_every_rN3_part_under_the_budget():
    """Cut each call at 80% of the budget from the estimate of its whole work; the measured work of the call, shared
    out evenly, plus the measured per-part overhead (the rebuild's estimate) stays under 300."""
    by = {c.name: c for c in cluster.calls()}
    groups = defaultdict(list)
    for p, est, meas in estimates("rN.3"):
        groups[p["step"]].append((p, est, meas))
    over = []
    for name, parts in groups.items():
        oh = (cluster.part_overhead_cpu_min(by[name], int(parts[0][0]["targets"]), shared=False)
              if parts[0][0]["targets"] else 0.0)
        work_est = sum(est - oh for _, est, _ in parts)
        k = max(1, math.ceil(work_est / (PART_BUDGET - oh)))
        work_meas = sum(meas - oh for _, _, meas in parts)
        if work_meas / k + oh > BUDGET:
            over.append((name, k, round(work_meas / k + oh, 1)))
    assert not over, over


def cut_doc(tmp_path, name: str, units: int, targets: int, plan_cut: dict, proteins: int = 10_620_684) -> dict:
    """The cut record the code makes for call `name` whose split input holds `units` residues or entry bytes in all
    (four entries; for a prefilter, whose query database is its target, `targets` entries); only the databases'
    `.index` files exist."""
    c = {x.name: x for x in cluster.calls()}[name]
    prefilter = c.module == "prefilter"
    src = tmp_path / f"{c.positional[cluster.PER_QUERY[c.module]]}.index"
    src.parent.mkdir(parents=True, exist_ok=True)
    n = targets if prefilter else 4
    pad = 2 if prefilter else 0  # a sequence entry is its residues, a newline and a NUL
    src.write_text("".join(f"{k}\t{k * 10}\t{units // n + pad}\n" for k in range(n)))
    if prefilter:
        assert c.positional[1] == c.positional[0]  # the cascade searches its input against itself
    ctx = cluster._Ctx(root=tmp_path, out=None, release="rX", plan={"proteins": proteins, "cut": plan_cut},
                       plan_sha256="", w=tmp_path, threads=1)
    return cluster._cut_doc(ctx, c)


def test_c22_a_cut_estimates_from_its_own_input_not_the_release_total(tmp_path):
    name, units = "cluster.nr90.12.align", 1_000_000_000
    got = cut_doc(tmp_path, name, units, 0, {"part_cpu_min": PART_BUDGET})
    assert got["cpu_min"] == pytest.approx(cluster.call_cpu_min(name, 10_620_684, units), abs=1e-5)
    # half the input at the same release size: half the work, and the parts follow it
    half = cut_doc(tmp_path, name, units // 2, 0, {"part_cpu_min": PART_BUDGET})
    assert half["cpu_min"] == pytest.approx(got["cpu_min"] / 2, rel=1e-6)
    assert len(got["parts"]) == math.ceil(got["cpu_min"] / PART_BUDGET) == 2
    assert all(p["cpu_min"] <= PART_BUDGET for p in got["parts"])


def test_c22_a_one_part_prefilter_that_shares_an_index_is_estimated_with_its_own_table_build(tmp_path):
    """A sensitive prefilter call cut into one part builds the table itself (no shared index is made for one part), so
    its estimate carries the rebuild, not the index's load; it stays one part only if that fits the part budget,
    else it is two parts that share the index and carry the load."""
    by = {x.name: x for x in cluster.calls()}
    name, targets = "cluster.nr50.19.prefilter", 400_000
    rebuild = cluster.part_overhead_cpu_min(by[name], targets, shared=False)
    load = cluster.part_overhead_cpu_min(by[name], targets, shared=True)
    assert load < rebuild / 4
    one = cut_doc(tmp_path, name, 1_000_000_000, targets, {"part_cpu_min": PART_BUDGET})
    assert len(one["parts"]) == 1 and "index_cpu_min" not in one
    assert one["overhead_cpu_min"] == pytest.approx(rebuild, abs=1e-5)
    assert one["parts"][0]["cpu_min"] == pytest.approx(one["cpu_min"] + rebuild, abs=1e-4)
    # just over: the work alone fits 240, the work and the rebuild do not
    units = int((PART_BUDGET - rebuild / 2) / 163.11e-9)
    two = cut_doc(tmp_path, name, units, targets, {"part_cpu_min": PART_BUDGET})
    assert two["cpu_min"] < PART_BUDGET < two["cpu_min"] + rebuild
    assert len(two["parts"]) == 2 and two["overhead_cpu_min"] == pytest.approx(load, abs=1e-5)
    assert two["index_cpu_min"] == pytest.approx(cluster.index_cpu_min(targets), abs=1e-5)
    # a forced single part has no index either
    forced = cut_doc(tmp_path, name, units, targets, {"parts": 1})
    assert len(forced["parts"]) == 1 and forced["overhead_cpu_min"] == pytest.approx(rebuild, abs=1e-5)
    # a call that cannot share one pays the rebuild in every part
    fast = cut_doc(tmp_path, "cluster.nr90.11.prefilter", units, targets, {"part_cpu_min": PART_BUDGET})
    assert fast["overhead_cpu_min"] == pytest.approx(
        cluster.part_overhead_cpu_min(by["cluster.nr90.11.prefilter"], targets, False), abs=1e-5)
