"""Ring R1: test_run rows and the stage gate query (data_layer_and_tests.md §4.12)."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

import pyarrow.parquet as pq
import pytest

from biobox import gate, testrun
from biobox.registry.schemas import conformance_problems, table


@pytest.mark.parametrize("name,ids", [
    ("test_r03_random_slices", ["R-03"]),
    ("test_k02_k03_spills", ["K-02", "K-03"]),
    ("test_ly06_slice", ["LY-06"]),
    ("test_mr14_soft_mask", ["MR-14"]),
    ("test_c10_occurrence", ["C-10"]),
    ("test_tables_conform", []),
    ("test_r04", ["R-04"]),
])
def test_ids_from_name(name, ids):
    assert testrun.ids_from_name(name) == ids


def test_expand_ranges():
    assert gate.expand("C-01..C-03") == ["C-01", "C-02", "C-03"]
    assert gate.expand("Q1..Q3") == ["Q1", "Q2", "Q3"]
    assert gate.expand("LY-09..LY-10") == ["LY-09", "LY-10"]
    assert gate.expand("S-02d") == ["S-02d"]


T0 = datetime(2026, 10, 30, 12, tzinfo=timezone.utc)


def _rows(outcomes: dict[str, str], ring="R2", t=T0):
    return [testrun.make_row(tid, ring, out, 0.1, t, f"tests/x.py::{tid}") for tid, out in outcomes.items()]


def test_session_file_conforms_and_is_indexed(tmp_path):
    p1 = testrun.write_session(_rows({"C-01": "pass"}), tmp_path, day=T0)
    p2 = testrun.write_session(_rows({"C-02": "fail"}), tmp_path, day=T0)
    assert [p1.name, p2.name] == ["2026-10-30-0.parquet", "2026-10-30-1.parquet"]
    assert conformance_problems(pq.read_table(p1), table("test_run")) == []
    assert gate._files(tmp_path) == [str(p1), str(p2)]


def test_gate_reports_missing_failing_and_rings(tmp_path, monkeypatch):
    specs = {"S9": {"B": ["C-01", "C-02", "C-03"], "W": ["R-08"], "rings": ["R2", "R3"], "milestones": []}}
    monkeypatch.setattr(gate, "load_gates", lambda: specs)
    testrun.write_session(_rows({"C-01": "pass", "C-02": "pass"}) + _rows({"G-01": "skip"}, ring="R3"), tmp_path, day=T0)
    rep = gate.gate_report(tmp_path, "S9")
    assert rep["B"] == {"passing": ["C-01", "C-02"], "failing": [], "missing": ["C-03"]}
    assert rep["rings"]["R2"]["green"] and not rep["rings"]["R3"]["green"]  # R3 only skipped
    assert not rep["passes"]
    # a later failing run of C-02 makes it failing; the latest R2 session is no longer green
    testrun.write_session(_rows({"C-02": "fail", "C-03": "pass"}, t=T0 + timedelta(hours=1)), tmp_path, day=T0)
    rep = gate.gate_report(tmp_path, "S9")
    assert rep["B"]["failing"] == ["C-02"] and rep["B"]["missing"] == []
    assert not rep["rings"]["R2"]["green"]


def test_skipped_rows_count_as_not_run(tmp_path, monkeypatch):
    specs = {"S9": {"B": ["C-01", "C-02"], "W": [], "rings": [], "milestones": []}}
    monkeypatch.setattr(gate, "load_gates", lambda: specs)
    testrun.write_session(_rows({"C-01": "skip"}) + _rows({"C-02": "skip"}) + _rows({"C-02": "pass"}), tmp_path, day=T0)
    rep = gate.gate_report(tmp_path, "S9")
    assert rep["B"] == {"passing": ["C-02"], "failing": [], "missing": ["C-01"]}


def test_gate_filters_by_image(tmp_path, monkeypatch):
    specs = {"S9": {"B": ["C-01"], "W": [], "rings": [], "milestones": []}}
    monkeypatch.setattr(gate, "load_gates", lambda: specs)
    monkeypatch.setenv("BIOBOX_IMAGE_DIGEST", "sha256:aaa")
    testrun.write_session(_rows({"C-01": "pass"}), tmp_path, day=T0)
    assert gate.gate_report(tmp_path, "S9", image="sha256:aaa")["passes"]
    assert gate.gate_report(tmp_path, "S9", image="sha256:bbb")["B"]["missing"] == ["C-01"]
