"""Ring R1: card A-13, the release train's table of records accepted by a decision row (DEC-43, DEC-53).

`train_est.py rates` stops the train (exit 3) on a judged step more than 20% over its estimate. Its table `ACCEPTED`
carries records the PI accepted, each by a decision row, as `c35_est.py` carried DEC-55 and `c47_est.py` DEC-56.
The table is empty in card A-13's commit. These tests inject entries with `monkeypatch` and drive `rates` on a small
root of run records constructed here (no data root, no MinCED, MMseqs2 or hmmsearch).
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("A-13", "C-46")  # ring R1 by its directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
POD_DIR = REPO / "scripts" / "stress" / "build_stage_pod"
SCRIPT = POD_DIR / "train_est.py"
DECISIONS = REPO / "planning" / "DECISIONS.md"

REL = "r26.11.2"
CALL = "cluster.nr50.20.align"
P0, P1, P2 = (f"{CALL}.part-{k:05d}" for k in range(3))
BASIS = "cut estimate 200.000 x factor 0.5"


def _te():
    if str(POD_DIR) not in sys.path:
        sys.path.insert(0, str(POD_DIR))
    return __import__("train_est")


def _record(root: Path, rel: str, step: str, cpu_min: float) -> None:
    p = root / "runs" / f"stage={rel}" / f"{step}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"step": step, "cpu_min": cpu_min, "wall_s": 10.0, "peak_rss_bytes": 100_000_000}))


def _root(tmp_path: Path, rel: str = REL, *, postings: float = 50.0) -> tuple[Path, Path]:
    """Three parts of one call (estimates 100 each: cut 200 x factor 0.5), two earlier-record steps, one unknown."""
    root = tmp_path / "root"
    plan = root / "l2" / "cluster_plan" / f"release={rel}" / f"{CALL}.cut.json"
    plan.parent.mkdir(parents=True)
    plan.write_text(json.dumps({"parts": [{"cpu_min": 200.0}] * 3}))
    _record(root, rel, P0, 100.0)    # ratio 1.000
    _record(root, rel, P1, 121.5)    # ratio 1.215: over 20%, judged
    _record(root, rel, P2, 70.0)     # ratio 0.700: under 20%, judged
    _record(root, rel, "postings", postings)
    _record(root, rel, "tokens", 30.0)   # estimate 10: never judged (DEC-43)
    _record(root, rel, "catalog", 1.0)   # no estimate
    est = tmp_path / "before.json"
    est.write_text(json.dumps({"projection": {"factors": {CALL: 0.5}, "steps": {"postings": 50.0, "tokens": 10.0}}}))
    return root, est


def _row(step: str, cpu: float, est: float, ratio: float, basis: str = BASIS, judged: bool = True) -> dict:
    return {"step": step, "cpu_min": cpu, "wall_s": 10.0, "peak_rss_gb": 0.1,
            "estimate": est, "ratio": ratio, "judged": judged, "basis": basis}


def _expected_empty_table(rel: str = REL) -> dict:
    """What `rates` printed before card A-13 on the root of `_root`, worked out by hand from the guard's rules."""
    p1 = _row(P1, 121.5, 100.0, 1.215)
    p2 = _row(P2, 70.0, 100.0, 0.7)
    return {
        "release": rel,
        "records": 6,
        "judged": 4,                        # three parts and postings; tokens is under 24 CPU-minutes
        "cpu_min": 372.5,                   # 1 + 100 + 121.5 + 70 + 50 + 30
        "estimate_of_recorded": 360.0,      # 0 + 100 + 100 + 100 + 50 + 10
        "worst_over": p1,
        "largest": p1,
        "over_20pct": [p1],
        "under_20pct": [p2],
        "outside_20pct": [p1, p2],
        "no_estimate": ["catalog"],
        "no_estimate_n": 1,
    }


def _rates(te, root: Path, est: Path, capsys, rel: str = REL) -> tuple[int, dict, str]:
    rc = te.rates(argparse.Namespace(root=str(root), release=rel, est=str(est), out=None))
    cap = capsys.readouterr()
    return rc, json.loads(cap.out), cap.err


# --- The empty table changes nothing ---------------------------------------------------------------------------


def test_empty_table_output_is_unchanged_to_the_byte(tmp_path, monkeypatch, capsys):
    te = _te()
    monkeypatch.setattr(te, "ACCEPTED", ())
    root, est = _root(tmp_path)
    rc = te.rates(argparse.Namespace(root=str(root), release=REL, est=str(est), out=None))
    assert rc == 3
    assert capsys.readouterr().out == json.dumps(_expected_empty_table(), sort_keys=True) + "\n"


def test_script_as_the_train_runs_it_is_unchanged_for_a_release_no_entry_names(tmp_path):
    """`release_train.sh` runs `python3 -I train_est.py rates --root R --release R --est E` and reads the exit status."""
    rel = "rA13.0"  # a name no entry of the repository's table can carry: the release-name pattern is r<YY>.<MM>.<n>
    root, est = _root(tmp_path, rel)
    res = subprocess.run([sys.executable, "-I", str(SCRIPT), "rates", "--root", str(root), "--release", rel,
                          "--est", str(est)], capture_output=True, text=True, check=False)
    assert res.returncode == 3, res.stderr
    assert res.stdout == json.dumps(_expected_empty_table(rel), sort_keys=True) + "\n"
    assert res.stderr == ""


# --- An entry accepts at or under its bound, and nothing else --------------------------------------------------


@pytest.mark.parametrize("bound", [1.215, 1.30, 1.60])
def test_accepted_at_and_under_the_bound(tmp_path, monkeypatch, capsys, bound):
    te = _te()
    monkeypatch.setattr(te, "ACCEPTED", (te.Accept(REL, f"{CALL}.part-*", bound, "DEC-55"),))
    root, est = _root(tmp_path)
    rc, doc, err = _rates(te, root, est, capsys)
    assert rc == 0, doc
    assert err == ""
    assert doc["over_20pct"] == []
    (acc,) = doc["accepted_over_20pct"]
    assert acc["step"] == P1 and acc["ratio"] == 1.215 and acc["accepted"] == "DEC-55"
    # Still outside 20% and still in the measurement note's inputs; the under-estimate is untouched.
    assert [r["step"] for r in doc["outside_20pct"]] == [P1, P2]
    assert [r["step"] for r in doc["under_20pct"]] == [P2]
    assert doc["largest"]["step"] == P1
    # worst_over holds only what is not accepted: the worst of the rest is a ratio of 1.0.
    assert doc["worst_over"]["step"] != P1 and doc["worst_over"]["ratio"] == 1.0


def test_stopped_above_the_bound_and_says_the_entry_was_exceeded(tmp_path, monkeypatch, capsys):
    te = _te()
    monkeypatch.setattr(te, "ACCEPTED", (te.Accept(REL, f"{CALL}.part-*", 1.21, "DEC-55"),))
    root, est = _root(tmp_path)
    rc, doc, _ = _rates(te, root, est, capsys)
    assert rc == 3
    assert doc["accepted_over_20pct"] == []
    (over,) = doc["over_20pct"]
    assert over["step"] == P1
    assert over["accepted_bound_exceeded"] == {"decision": "DEC-55", "step": f"{CALL}.part-*", "max_ratio": 1.21}
    assert doc["worst_over"]["step"] == P1


def test_stopped_for_another_release(tmp_path, monkeypatch, capsys):
    te = _te()
    monkeypatch.setattr(te, "ACCEPTED", (te.Accept("r26.11.3", f"{CALL}.part-*", 1.60, "DEC-55"),
                                         te.Accept("r26.11.1", f"{CALL}.part-*", 1.60, "DEC-55")))
    root, est = _root(tmp_path)
    rc, doc, _ = _rates(te, root, est, capsys)
    assert rc == 3
    assert doc == _expected_empty_table()  # no entry names this release: the output is as with the empty table


def test_stopped_for_another_step_at_the_same_ratio(tmp_path, monkeypatch, capsys):
    te = _te()
    monkeypatch.setattr(te, "ACCEPTED", (te.Accept(REL, f"{CALL}.part-*", 1.30, "DEC-55"),))
    root, est = _root(tmp_path, postings=60.75)  # 60.75 / 50 = 1.215, as part-00001
    rc, doc, _ = _rates(te, root, est, capsys)
    assert rc == 3
    assert [r["step"] for r in doc["over_20pct"]] == ["postings"]
    assert "accepted_bound_exceeded" not in doc["over_20pct"][0]
    assert [r["step"] for r in doc["accepted_over_20pct"]] == [P1]
    assert doc["worst_over"]["step"] == "postings"


def test_step_name_is_matched_whole_and_case_sensitively(tmp_path, monkeypatch, capsys):
    te = _te()
    for k, pattern in enumerate((f"{CALL}.part-00000", f"{CALL}", f"{CALL.upper()}.PART-*", "nr50.20.align.part-*")):
        monkeypatch.setattr(te, "ACCEPTED", (te.Accept(REL, pattern, 1.60, "DEC-55"),))
        root, est = _root(tmp_path / f"p{k}")
        rc, doc, _ = _rates(te, root, est, capsys)
        assert rc == 3, pattern
        assert doc["accepted_over_20pct"] == [], pattern


def test_step_under_0_80_is_unaffected(tmp_path, monkeypatch, capsys):
    te = _te()
    monkeypatch.setattr(te, "ACCEPTED", (te.Accept(REL, f"{CALL}.part-0000[12]", 1.30, "DEC-55"),))
    root, est = _root(tmp_path)
    rc, doc, _ = _rates(te, root, est, capsys)
    assert rc == 0
    (under,) = doc["under_20pct"]
    assert under == _row(P2, 70.0, 100.0, 0.7)  # no `accepted` field on a record under its estimate
    assert [r["step"] for r in doc["accepted_over_20pct"]] == [P1]


def test_unjudged_step_is_unaffected(tmp_path, monkeypatch, capsys):
    te = _te()
    monkeypatch.setattr(te, "ACCEPTED", (te.Accept(REL, "tokens", 1.30, "DEC-55"),
                                         te.Accept(REL, f"{CALL}.part-*", 1.30, "DEC-55")))
    root, est = _root(tmp_path)
    rc, doc, _ = _rates(te, root, est, capsys)
    assert rc == 0
    assert [r["step"] for r in doc["accepted_over_20pct"]] == [P1]  # tokens (ratio 3.0, estimate 10) is not judged


def test_two_matching_entries_the_tighter_bound_decides(tmp_path, monkeypatch, capsys):
    te = _te()
    wide = te.Accept(REL, f"{CALL}.part-*", 1.30, "DEC-55")
    tight = te.Accept(REL, f"{CALL}.part-00001", 1.21, "DEC-56")
    for k, table in enumerate(((wide, tight), (tight, wide))):  # order in the table does not matter
        monkeypatch.setattr(te, "ACCEPTED", table)
        root, est = _root(tmp_path / f"t{k}")
        rc, doc, _ = _rates(te, root, est, capsys)
        assert rc == 3
        assert doc["over_20pct"][0]["accepted_bound_exceeded"]["decision"] == "DEC-56"
    tight_ok = te.Accept(REL, f"{CALL}.part-0000?", 1.25, "DEC-56")
    monkeypatch.setattr(te, "ACCEPTED", (wide, tight_ok))
    root, est = _root(tmp_path / "ok")
    rc, doc, _ = _rates(te, root, est, capsys)
    assert rc == 0
    assert doc["accepted_over_20pct"][0]["accepted"] == "DEC-56"


# --- A malformed table stops the train with a message, never ignored ------------------------------------------


MALFORMED = {
    "bound at 1.20": (REL, f"{CALL}.part-*", 1.20, "DEC-55"),
    "bound under 1.20": (REL, f"{CALL}.part-*", 1.10, "DEC-55"),
    "bound above the ceiling": (REL, f"{CALL}.part-*", 1.601, "DEC-55"),
    "bound far above": (REL, f"{CALL}.part-*", 5.0, "DEC-55"),
    "bound not a number": (REL, f"{CALL}.part-*", "1.3", "DEC-55"),
    "bound a bool": (REL, f"{CALL}.part-*", True, "DEC-55"),
    "bound nan": (REL, f"{CALL}.part-*", float("nan"), "DEC-55"),
    "empty pattern": (REL, "", 1.30, "DEC-55"),
    "star alone": (REL, "*", 1.30, "DEC-55"),
    "wildcards alone": (REL, "*?*", 1.30, "DEC-55"),
    "decision without number": (REL, f"{CALL}.part-*", 1.30, "DEC-"),
    "decision lower case": (REL, f"{CALL}.part-*", 1.30, "dec-55"),
    "decision with suffix": (REL, f"{CALL}.part-*", 1.30, "DEC-55a"),
    "decision bare number": (REL, f"{CALL}.part-*", 1.30, "55"),
    "decision leading zero": (REL, f"{CALL}.part-*", 1.30, "DEC-055"),
    "empty release": ("", f"{CALL}.part-*", 1.30, "DEC-55"),
    "release pattern": ("r26.11.*", f"{CALL}.part-*", 1.30, "DEC-55"),
    "three fields": (REL, f"{CALL}.part-*", 1.30),
}


@pytest.mark.parametrize("why", sorted(MALFORMED))
def test_malformed_entry_stops_with_a_message(tmp_path, monkeypatch, capsys, why):
    te = _te()
    monkeypatch.setattr(te, "ACCEPTED", (MALFORMED[why],))
    root, est = _root(tmp_path)
    for f in (root / "runs" / f"stage={REL}").glob("*.json"):  # nothing over 20%: only the table can stop it
        if f.stem == P1:
            f.unlink()
    rc, doc, err = _rates(te, root, est, capsys)
    assert rc == 3  # the one status on which release_train.sh stops the unit
    assert "ACCEPTED" in err and doc["accepted_table_error"]
    with pytest.raises(te.AcceptedTableError):
        te.check_rates(root, REL, est)


def test_the_same_release_and_pattern_twice_is_refused(tmp_path, monkeypatch, capsys):
    te = _te()
    monkeypatch.setattr(te, "ACCEPTED", (te.Accept(REL, f"{CALL}.part-*", 1.30, "DEC-55"),
                                         te.Accept(REL, f"{CALL}.part-*", 1.25, "DEC-56")))
    root, est = _root(tmp_path)
    rc, doc, err = _rates(te, root, est, capsys)
    assert rc == 3
    assert "twice" in err


def test_the_ceiling_itself_is_accepted_as_a_bound():
    te = _te()
    assert te.ACCEPT_CEILING == 1.60  # DEC-55, the largest ratio the PI has accepted
    assert te.accepted_table_problems((te.Accept(REL, f"{CALL}.part-*", 1.60, "DEC-55"),)) == []


# --- The repository's own table --------------------------------------------------------------------------------


def _decision_rows() -> set[str]:
    return set(re.findall(r"^\| (DEC-[0-9]+) \|", DECISIONS.read_text(), flags=re.MULTILINE))


def _entries_without_a_row(table) -> list[str]:
    rows = _decision_rows()
    return [e[3] for e in table if e[3] not in rows]


def test_every_entry_has_a_decision_row():
    te = _te()
    assert {"DEC-43", "DEC-53", "DEC-55", "DEC-56"} <= _decision_rows()  # the reader of DECISIONS.md works
    assert te.accepted_table_problems(te.ACCEPTED) == []
    assert _entries_without_a_row(te.ACCEPTED) == []
    # An entry naming a row that does not exist is caught.
    assert _entries_without_a_row((te.Accept(REL, f"{CALL}.part-*", 1.30, "DEC-99999"),)) == ["DEC-99999"]


def test_nothing_outside_the_script_can_add_an_entry():
    src = SCRIPT.read_text()
    assert "environ" not in src and "getenv" not in src
    res = subprocess.run([sys.executable, "-I", str(SCRIPT), "--help"], capture_output=True, text=True, check=True)
    assert "accept" not in res.stdout.lower()
