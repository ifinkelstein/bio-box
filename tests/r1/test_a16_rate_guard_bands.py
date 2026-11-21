"""Ring R1: card A-16, the rate guards' two lines (DEC-64, amending DEC-43 and DEC-53).

A judged step (estimated at 24 CPU-minutes or more) more than 1.50 times its estimate stops the build, as a step over
1.20 did before; one over 1.20 and at or under 1.50 is listed under `recorded_over_20pct`, stays in `outside_20pct`,
and the build goes on. The line that stops is compared on the reported three-decimal ratio, as card A-13's bounds are.
Both guards: the release train's `train_est.py rates` and the stage driver's `stage_est.py rates` (with its `JUDGE`
figures and conf `ACCEPTED` lines), and `stage_check.py`, which fails a stage's check on what the guard stops on.

The roots are constructed here (no data root, no MinCED, MMseqs2 or hmmsearch: the CRISPR producer's name comes from
a fixed MinCED version set with `monkeypatch`, as tests/r1/test_stage_driver.py does).
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("A-16", "A-13", "A-12", "C-46")  # ring R1 by its directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
TRAIN = BSP / "train_est.py"
STAGE = BSP / "stage_est.py"
DECISIONS = REPO / "planning" / "DECISIONS.md"

CALL = "cluster.nr50.20.align"
P0, P1 = f"{CALL}.part-00000", f"{CALL}.part-00001"
DEC61 = ("r26.11.2", f"{CALL}.part-*", 1.50, "DEC-61")


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod  # a dataclass looks its module up by name
    spec.loader.exec_module(mod)
    return mod


def _te():
    return _load("train_est", TRAIN)


def _se():
    return _load("stage_est", STAGE)


@pytest.fixture
def no_minced(monkeypatch):
    """The producer names without the MinCED binary: a fixed version."""
    from biobox.l1 import crispr

    monkeypatch.setattr(crispr, "minced_version", lambda: "0.4.2")


# --- the constants ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("path", [TRAIN, STAGE], ids=["train_est", "stage_est"])
def test_two_named_lines_with_the_decision_row_where_they_are_defined(path):
    mod = _load(path.stem, path)
    assert mod.RECORD_ABOVE == 1.20 and mod.STOP_ABOVE == 1.50
    assert mod.RECORD_ABOVE == 1 + mod.TOL  # the line DEC-53 stopped on is the line that now records
    src = path.read_text()
    for name, value in (("RECORD_ABOVE", "1.20"), ("STOP_ABOVE", "1.50")):
        (line,) = [x for x in src.splitlines() if x.startswith(f"{name} = ")]
        assert line.startswith(f"{name} = {value}") and "DEC-64" in line, line
    # unchanged: DEC-43's least judged estimate and the line under which a step is recorded (1 - TOL = 0.80)
    assert mod.MIN_JUDGED_CPU_MIN == 24.0 and mod.TOL == 0.20


# --- train_est.py rates ----------------------------------------------------------------------------------------


def _record(root: Path, rel: str, step: str, cpu_min: float) -> None:
    p = root / "runs" / f"stage={rel}" / f"{step}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"step": step, "cpu_min": cpu_min, "wall_s": 10.0, "peak_rss_bytes": 100_000_000}))


def _train_root(tmp_path: Path, rel: str, p0: float, *, p1: float = 100.0, extra: dict | None = None):
    """Two parts of one call, each estimated at exactly 100 CPU-minutes (cut 200 x factor 0.5): the ratio is cpu/100."""
    root = tmp_path / "root"
    plan = root / "l2" / "cluster_plan" / f"release={rel}" / f"{CALL}.cut.json"
    plan.parent.mkdir(parents=True)
    plan.write_text(json.dumps({"parts": [{"cpu_min": 200.0}] * 2}))
    _record(root, rel, P0, p0)
    _record(root, rel, P1, p1)
    steps = {}
    for step, (cpu, e) in (extra or {}).items():
        _record(root, rel, step, cpu)
        steps[step] = e
    est = tmp_path / "before.json"
    est.write_text(json.dumps({"projection": {"factors": {CALL: 0.5}, "steps": steps}}))
    return root, est


def _rates(te, root, est, capsys, rel):
    rc = te.rates(argparse.Namespace(root=str(root), release=rel, est=str(est), out=None))
    cap = capsys.readouterr()
    return rc, json.loads(cap.out)


def _steps(rows):
    return [r["step"] for r in rows]


# cpu of part-00000 (estimate 100): (reported ratio, recorded in the band, stops)
TRAIN_LINES = {
    "1.200": (120.0, 1.2, False, False),     # at the line that reports: inside 20%
    "1.201": (120.1, 1.201, True, False),    # over it: recorded, goes on
    "1.215": (121.5, 1.215, True, False),    # DEC-61's part, under the rule alone
    "1.500": (150.0, 1.5, True, False),      # at the line that stops: recorded, goes on
    "1.5004": (150.04, 1.5, True, False),    # reported 1.500: the line that stops reads the reported ratio
    "1.501": (150.1, 1.501, False, True),    # over the line that stops
    "1.600": (160.0, 1.6, False, True),
}


@pytest.mark.parametrize("case", list(TRAIN_LINES))
def test_train_both_sides_of_both_lines(tmp_path, monkeypatch, capsys, case):
    te = _te()
    monkeypatch.setattr(te, "ACCEPTED", ())
    cpu, ratio, recorded, stops = TRAIN_LINES[case]
    root, est = _train_root(tmp_path, "rA16.0", cpu)
    rc, doc = _rates(te, root, est, capsys, "rA16.0")
    (row,) = [r for r in doc["outside_20pct"] if r["step"] == P0] or [None]
    assert rc == (3 if stops else 0), doc
    assert _steps(doc["over_20pct"]) == ([P0] if stops else [])
    assert _steps(doc["recorded_over_20pct"]) == ([P0] if recorded else [])
    if recorded or stops:
        assert row is not None and row["ratio"] == ratio and row["estimate"] == 100.0
    else:
        assert row is None, "1.200 is not outside 20%"
    assert doc["under_20pct"] == []
    assert "accepted_over_20pct" not in doc  # no entry names this release


def test_train_recorded_row_names_its_step_and_ratio_and_worst_over_is_as_before(tmp_path, monkeypatch, capsys):
    te = _te()
    monkeypatch.setattr(te, "ACCEPTED", ())
    root, est = _train_root(tmp_path, "rA16.0", 130.0, p1=110.0)
    rc, doc = _rates(te, root, est, capsys, "rA16.0")
    assert rc == 0
    (rec,) = doc["recorded_over_20pct"]
    assert rec == {"step": P0, "cpu_min": 130.0, "wall_s": 10.0, "peak_rss_gb": 0.1, "estimate": 100.0,
                   "ratio": 1.3, "judged": True, "basis": "cut estimate 200.000 x factor 0.5"}
    assert doc["outside_20pct"] == [rec] and doc["over_20pct"] == []
    assert doc["worst_over"] == rec  # the worst judged ratio, as before card A-16


def test_train_a_step_just_over_1_20_unrounded_is_recorded_as_before(tmp_path, monkeypatch, capsys):
    """The line that reports is DEC-53's comparison, unchanged: 1.2004 (reported 1.200) was outside 20% and stopped;
    now it is recorded and the train goes on."""
    te = _te()
    monkeypatch.setattr(te, "ACCEPTED", ())
    root, est = _train_root(tmp_path, "rA16.0", 120.04)
    rc, doc = _rates(te, root, est, capsys, "rA16.0")
    assert rc == 0
    assert _steps(doc["recorded_over_20pct"]) == [P0] and _steps(doc["outside_20pct"]) == [P0]


@pytest.mark.parametrize("cpu,code", [(135.0, 0), (150.0, 0), (150.1, 3), (240.0, 3)])
def test_train_script_as_the_train_runs_it(tmp_path, cpu, code):
    """`release_train.sh` runs `python3 -I train_est.py rates ...` and acts on the exit status alone."""
    root, est = _train_root(tmp_path, "rA16.0", cpu)
    res = subprocess.run([sys.executable, "-I", str(TRAIN), "rates", "--root", str(root), "--release", "rA16.0",
                          "--est", str(est)], capture_output=True, text=True, check=False)
    assert res.returncode == code, res.stderr
    doc = json.loads(res.stdout)
    assert _steps(doc["over_20pct"]) == ([P0] if code else [])
    assert _steps(doc["recorded_over_20pct"]) == ([] if code else [P0])
    assert res.stderr == ""


@pytest.mark.parametrize("cpu", [79.9, 50.0, 10.0])
def test_train_a_step_under_0_80_is_recorded_and_never_stops(tmp_path, monkeypatch, capsys, cpu):
    te = _te()
    monkeypatch.setattr(te, "ACCEPTED", ())
    root, est = _train_root(tmp_path, "rA16.0", cpu)
    rc, doc = _rates(te, root, est, capsys, "rA16.0")
    assert rc == 0
    assert _steps(doc["under_20pct"]) == [P0] and _steps(doc["outside_20pct"]) == [P0]
    assert doc["over_20pct"] == [] and doc["recorded_over_20pct"] == []


def test_train_a_step_estimated_under_24_is_not_judged(tmp_path, monkeypatch, capsys):
    te = _te()
    monkeypatch.setattr(te, "ACCEPTED", ())
    root, est = _train_root(tmp_path, "rA16.0", 100.0, extra={"tokens": (230.0, 23.9), "dim": (1.0, 10.0)})
    rc, doc = _rates(te, root, est, capsys, "rA16.0")
    assert rc == 0, "tokens at 9.6 times an estimate of 23.9 is not judged (DEC-43)"
    assert doc["judged"] == 2
    assert doc["outside_20pct"] == [] and doc["recorded_over_20pct"] == [] and doc["over_20pct"] == []


# --- ACCEPTED: an entry matters only above 1.50 ---------------------------------------------------------------


def test_the_dec61_entry_is_kept_unchanged_and_valid():
    te = _te()
    assert DEC61 in [tuple(e) for e in te.ACCEPTED]
    assert te.accepted_table_problems(te.ACCEPTED) == []
    assert "| DEC-61 |" in DECISIONS.read_text()


def test_the_allowance_for_entries_made_before_dec64_is_dec61s_entry_alone():
    te = _te()
    assert [tuple(e) for e in te.KEPT_BEFORE_DEC64] == [DEC61]
    for e in te.KEPT_BEFORE_DEC64:
        assert int(e[3].removeprefix("DEC-")) < 64, "only a row older than DEC-64 can be kept"
        assert te.RECORD_ABOVE < e[2] <= te.STOP_ABOVE


NEW_ENTRIES_REFUSED = {
    "bound at 1.50": ("r26.11.3", f"{CALL}.part-*", 1.50, "DEC-65"),
    "bound 1.30": ("r26.11.3", f"{CALL}.part-*", 1.30, "DEC-65"),
    "bound 1.201": ("r26.11.3", f"{CALL}.part-*", 1.201, "DEC-65"),
    "DEC-61's bound for another release": ("r26.11.3", f"{CALL}.part-*", 1.50, "DEC-61"),
    "DEC-61's bound for another pattern": ("r26.11.2", "cluster.nr50.19.prefilter.part-*", 1.50, "DEC-61"),
    "DEC-61's entry under another row": ("r26.11.2", f"{CALL}.part-*", 1.50, "DEC-62"),
    "DEC-61's entry with a lower bound": ("r26.11.2", f"{CALL}.part-*", 1.40, "DEC-61"),
    "an older row with a band bound": ("r26.11.1", f"{CALL}.part-*", 1.30, "DEC-55"),
}


@pytest.mark.parametrize("why", sorted(NEW_ENTRIES_REFUSED))
def test_a_new_entry_at_or_under_1_50_is_refused_and_stops_the_train(tmp_path, monkeypatch, capsys, why):
    te = _te()
    entry = NEW_ENTRIES_REFUSED[why]
    problems = te.accepted_table_problems((entry,))
    assert len(problems) == 1 and "1.50" in problems[0] and "DEC-64" in problems[0], problems
    monkeypatch.setattr(te, "ACCEPTED", (DEC61, entry) if entry[:2] != DEC61[:2] else (entry,))
    root, est = _train_root(tmp_path, entry[0], 100.0)  # nothing outside 20%: only the table can stop it
    rc, doc = _rates(te, root, est, capsys, entry[0])
    assert rc == 3 and "DEC-64" in doc["accepted_table_error"]


@pytest.mark.parametrize("bound", [1.501, 1.55, 1.60])
def test_a_new_entry_over_1_50_and_at_most_the_ceiling_is_valid(bound):
    te = _te()
    assert te.ACCEPT_CEILING == 1.60
    assert te.accepted_table_problems((("r26.11.3", f"{CALL}.part-*", bound, "DEC-65"),)) == []
    assert te.accepted_table_problems((("r26.11.3", f"{CALL}.part-*", 1.601, "DEC-65"),)) != []


def test_the_dec61_entry_accepts_nothing_the_rule_does_not_already_pass(tmp_path, monkeypatch, capsys):
    """With and without DEC-61's entry, every ratio of its own steps stops or goes on alike."""
    te = _te()
    kept = te.ACCEPTED
    for k, cpu in enumerate((119.9, 120.0, 120.1, 121.5, 135.0, 149.9, 150.0, 150.04, 150.1, 151.0, 160.0, 200.0)):
        out = {}
        for name, table in (("kept", kept), ("empty", ())):
            monkeypatch.setattr(te, "ACCEPTED", table)
            root, est = _train_root(tmp_path / f"{k}-{name}", "r26.11.2", cpu)
            rc, doc = _rates(te, root, est, capsys, "r26.11.2")
            out[name] = (rc, _steps(doc["over_20pct"]), _steps(doc["recorded_over_20pct"]),
                         _steps(doc["outside_20pct"]), _steps(doc["under_20pct"]))
        assert out["kept"] == out["empty"], cpu


def test_dec61s_part_at_1_215_is_recorded_by_the_rule_and_matched_by_its_entry(tmp_path, monkeypatch, capsys):
    """The stop of r26.11.2: part-00001 at 190.462 against 156.713 (1.215). Under DEC-64 it is in the band; DEC-61's
    entry still matches it, so it is listed under both keys, and the train goes on as it does without the entry."""
    te = _te()
    root, est = _train_root(tmp_path, "r26.11.2", 100.0, p1=121.5)
    rc, doc = _rates(te, root, est, capsys, "r26.11.2")
    assert rc == 0
    (rec,) = doc["recorded_over_20pct"]
    (acc,) = doc["accepted_over_20pct"]
    assert rec == acc and rec["step"] == P1 and rec["ratio"] == 1.215 and rec["accepted"] == "DEC-61"
    assert _steps(doc["outside_20pct"]) == [P1] and doc["over_20pct"] == []
    assert doc["worst_over"]["step"] == P0  # as card A-13 made it: an accepted record is not the worst over


def test_dec61s_parts_above_1_50_still_stop_naming_the_entry(tmp_path, capsys):
    te = _te()
    root, est = _train_root(tmp_path, "r26.11.2", 100.0, p1=150.1)
    rc, doc = _rates(te, root, est, capsys, "r26.11.2")
    assert rc == 3
    (over,) = doc["over_20pct"]
    assert over["step"] == P1 and over["accepted_bound_exceeded"] == {"decision": "DEC-61", "step": f"{CALL}.part-*",
                                                                     "max_ratio": 1.5}
    assert doc["recorded_over_20pct"] == [] and doc["accepted_over_20pct"] == []


@pytest.mark.parametrize("cpu,code", [(155.0, 0), (155.1, 3)])
def test_a_new_entry_accepts_above_1_50_up_to_its_bound(tmp_path, monkeypatch, capsys, cpu, code):
    te = _te()
    monkeypatch.setattr(te, "ACCEPTED", (te.Accept("r26.11.3", f"{CALL}.part-*", 1.55, "DEC-65"),))
    root, est = _train_root(tmp_path, "r26.11.3", cpu)
    rc, doc = _rates(te, root, est, capsys, "r26.11.3")
    assert rc == code
    assert doc["recorded_over_20pct"] == [], "a step over 1.50 is never in the band"
    assert _steps(doc["accepted_over_20pct"]) == ([] if code else [P0])
    assert _steps(doc["over_20pct"]) == ([P0] if code else [])


def test_a_malformed_table_still_stops_the_train(tmp_path, monkeypatch, capsys):
    te = _te()
    monkeypatch.setattr(te, "ACCEPTED", (("r26.11.3", "*", 1.55, "DEC-65"),))
    root, est = _train_root(tmp_path, "r26.11.3", 100.0)
    rc, doc = _rates(te, root, est, capsys, "r26.11.3")
    assert rc == 3 and doc["accepted_table_error"]


# --- stage_est.py: judge_items, the JUDGE figures, the conf's ACCEPTED lines ----------------------------------


def _item(name, cpu, est, complete=True):
    return {"name": name, "cpu_min": cpu, "estimate": est, "basis": "test", "complete": complete, "records": 1}


STAGE_LINES = {  # cpu against an estimate of 100: (recorded, stops)
    "1.200": (120.0, False, False),
    "1.201": (120.1, True, False),
    "1.500": (150.0, True, False),
    "1.5004": (150.04, True, False),
    "1.501": (150.1, False, True),
}


@pytest.mark.parametrize("complete", [True, False], ids=["complete", "running"])
@pytest.mark.parametrize("name", ["ingest", "genecall.note"], ids=["code", "judge-figure"])
@pytest.mark.parametrize("case", list(STAGE_LINES))
def test_stage_both_sides_of_both_lines(case, name, complete):
    """The code's estimates and a conf's JUDGE figure (the group `<step>.note`) follow the same two lines; a group
    still running is judged on its partial sum, as before."""
    se = _se()
    cpu, recorded, stops = STAGE_LINES[case]
    doc = se.judge_items([_item(name, cpu, 100.0, complete)], {})
    assert doc["over_20pct"] == ([name] if stops else [])
    assert doc["recorded_over_20pct"] == ([name] if recorded else [])
    assert se.rates_exit(doc) == (3 if stops else 0)
    (row,) = doc["rows"]
    assert row["over_20pct"] is stops and row["recorded_over_20pct"] is recorded


def test_stage_under_0_80_is_recorded_and_never_stops_and_under_24_is_not_judged():
    se = _se()
    doc = se.judge_items([_item("crispr-x", 10.0, 100.0), _item("tapes", 230.0, 23.9)], {})
    assert doc["under_20pct"] == ["crispr-x"] and doc["over_20pct"] == [] and doc["recorded_over_20pct"] == []
    assert se.rates_exit(doc) == 0
    assert [r["judged"] for r in doc["rows"]] == [True, False]


def test_stage_an_accepted_line_matters_only_above_1_50():
    se = _se()
    acc = {("ingest", 160.0, 100.0): "DEC-65"}
    doc = se.judge_items([_item("ingest", 160.0, 100.0)], acc)
    assert doc["over_20pct"] == [] and doc["accepted_over_20pct"] == {"ingest": "DEC-65"}
    assert doc["recorded_over_20pct"] == [] and se.rates_exit(doc) == 0
    other = se.judge_items([_item("ingest", 160.1, 100.0)], acc)
    assert other["over_20pct"] == ["ingest"] and se.rates_exit(other) == 3


def test_stage_dec56s_record_is_in_the_band_and_its_line_changes_nothing():
    """DEC-56's line (C5's ingest, 45.803 against 38.057: 1.204) is a record of what was accepted before DEC-64: with
    or without it the step is recorded and the build goes on; with it the line is still reported."""
    se = _se()
    item = _item("ingest", 45.803, 38.057)
    with_line = se.judge_items([item], {("ingest", 45.803, 38.057): "DEC-56"})
    without = se.judge_items([item], {})
    for doc in (with_line, without):
        assert doc["over_20pct"] == [] and doc["recorded_over_20pct"] == ["ingest"] and se.rates_exit(doc) == 0
    assert with_line["accepted_over_20pct"] == {"ingest": "DEC-56"} and without["accepted_over_20pct"] == {}


def test_stage_guard_policy_names_both_lines():
    src = STAGE.read_text()
    assert src.count("DEC-64 (record above 1.20; stop above 1.50)") == 2  # rates and before


# --- stage_est.py rates and stage_check.py on a constructed root ----------------------------------------------


def _driver_tests():
    return _load("a16_stage_driver_tests", REPO / "tests" / "r1" / "test_stage_driver.py")


def _batch_root(tmp_path, cpu: float, estimate: float = 100.0, judge: str = ""):
    from biobox.l1 import genecall

    td = _driver_tests()
    batch = f"genecall.{genecall.batch_name(0)}"

    def records(S):
        S["records"]["ZQ"][batch] = []
        S["cpu"][batch] = cpu
        S["gcplan"]["batches"][0]["estimate_cpu_min"] = estimate

    r = td.build_root(tmp_path, records)
    if judge:
        r["conf"].write_text(r["conf"].read_text() + judge)
    return td, r, batch


@pytest.mark.parametrize("cpu,recorded,stops", [(120.0, False, False), (130.0, True, False), (150.0, True, False),
                                                (150.1, False, True)])
def test_stage_rates_doc_on_a_published_plan(tmp_path, no_minced, cpu, recorded, stops):
    from biobox import profile as bprofile

    se = _se()
    _, r, batch = _batch_root(tmp_path, cpu)
    doc = se.rates_doc(se.read_conf(r["conf"], "ZQ"), r["root"], bprofile.load(r["profile"]), tmp_path / "none",
                       str(r["before_json"]), tmp_path / "t")
    assert doc["over_20pct"] == ([batch] if stops else [])
    assert doc["recorded_over_20pct"] == ([batch] if recorded else [])
    assert se.rates_exit(doc) == (3 if stops else 0)


@pytest.mark.parametrize("cpu,stops", [(37.5, False), (37.6, True)])
def test_stage_rates_doc_judges_a_judge_figure_by_the_same_lines(tmp_path, no_minced, cpu, stops):
    """JUDGE=genecall 25: the group genecall.note is every genecall record (the batch and the merge's 0.1) against the
    note's 25; the batch alone, against the plan's 40, is inside 20%."""
    from biobox import profile as bprofile

    se = _se()
    judge = "JUDGE=genecall 25 planning/measurements/2026-11-20_stage_V1.md\n"
    _, r, batch = _batch_root(tmp_path, cpu - 0.1, estimate=40.0, judge=judge)
    doc = se.rates_doc(se.read_conf(r["conf"], "ZQ"), r["root"], bprofile.load(r["profile"]), tmp_path / "none",
                       str(r["before_json"]), tmp_path / "t")
    by = {x["step"]: x for x in doc["rows"]}
    assert by["genecall.note"]["cpu_min"] == pytest.approx(cpu)
    assert by["genecall.note"]["ratio"] == (1.504 if stops else 1.5)
    assert doc["over_20pct"] == (["genecall.note"] if stops else [])
    assert doc["recorded_over_20pct"] == ([] if stops else ["genecall.note"])
    assert se.rates_exit(doc) == (3 if stops else 0)


def test_stage_check_passes_a_stage_with_a_step_recorded_in_the_band(tmp_path, monkeypatch, no_minced):
    td, r, batch = _batch_root(tmp_path, 130.0)
    rc, s = td.run_check(r, monkeypatch)
    assert rc == 0, s["failed"]
    assert s["steps"]["rates"]["recorded_over_20pct"] == [batch]


def test_stage_check_fails_a_stage_with_a_step_over_1_50(tmp_path, monkeypatch, no_minced):
    td, r, batch = _batch_root(tmp_path, 150.1)
    rc, s = td.run_check(r, monkeypatch)
    assert rc == 1
    assert any(batch in f and "rates over" in f for f in s["failed"]), s["failed"]


# --- what the drivers read, and the guards that are not in these scripts --------------------------------------


DRIVERS = [BSP / "release_train.sh", BSP / "stage_driver.sh", BSP / "stage_queue.sh",
           REPO / "scripts" / "orchestrator" / "rel_wait.sh"]
KEYS = ("over_20pct", "recorded_over_20pct", "outside_20pct", "under_20pct", "accepted_over_20pct", "worst_over",
        "exit_code")


@pytest.mark.parametrize("path", DRIVERS, ids=lambda p: p.name)
def test_the_drivers_read_the_rate_guard_by_its_exit_status_alone(path):
    src = path.read_text()
    assert not [k for k in KEYS if k in src], "a driver that reads a key of the guard's output"
    if path.name == "release_train.sh":
        assert re.search(r'\[ "\$rc" = 0 \] \|\| stop_unit "\$unit" "rate guard', src)
        assert re.search(r'\[ "\$rc" = 0 \] \|\| fail "rate guard at end', src)
    if path.name == "stage_driver.sh":
        assert '[ "$RATES_RC" = 0 ] || stop_unit' in src and '[ "$RATES_RC" = 0 ] || fail' in src


def test_the_process_and_space_guards_are_elsewhere_and_unchanged():
    for est in (TRAIN, STAGE):
        src = est.read_text()
        assert "r50_watch" not in src and "CPU_MAX" not in src and "disk_usage" not in src and "statvfs" not in src
    rt = (BSP / "release_train.sh").read_text()
    sd = (BSP / "stage_driver.sh").read_text()
    for line in ("CPU_MAX=${CPU_MAX:-280}", "TMP_MIN_GB=${TMP_MIN_GB:-150}", "SCRATCH_MAX_GB=${SCRATCH_MAX_GB:-250}"):
        assert line in rt, line
    for line in ("CPU_MAX=${CPU_MAX:-280}", "TMP_MIN_GB=${TMP_MIN_GB:-200}", "SCRATCH_MAX_GB=${SCRATCH_MAX_GB:-100}",
                 '[ "$CPU_MAX" -le 280 ] || refuse'):
        assert line in sd, line
    for src in (rt, sd):
        assert "r50_watch.py" in src and '--cpu-max "$CPU_MAX"' in src
