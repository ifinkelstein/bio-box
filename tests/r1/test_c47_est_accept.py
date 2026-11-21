"""Ring R1: card C-47's rate guard accepts exactly one record, by name and value (DEC-56).

DEC-56 (PI, 2026-11-20): stage C5's ingest, 45.803 CPU-minutes against the code's estimate of 38.057 (ratio 1.204),
is accepted and reported; every other step, and any other ingest record, is still judged as DEC-43 and DEC-53 say.
The tests drive `c47_est.py rates` end to end on a small root of run records written here.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-47")  # ring R1 by directory (tests/conftest.py)

EST_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "stress" / "build_stage_pod" / "c47_est.py"
ST = "C5"


def _load():
    spec = importlib.util.spec_from_file_location("c47_est_accept", EST_SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c47_est_accept"] = mod
    spec.loader.exec_module(mod)
    return mod


def _record(root: Path, step: str, cpu_min: float) -> None:
    p = root / "runs" / f"stage={ST}" / f"{step}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"cpu_min": cpu_min, "wall_s": 60.0, "peak_rss_bytes": 2**20, "outputs": []}))


def _rates(est, monkeypatch, tmp_path: Path, *, ingest: float, ingest_est: float,
           merge: tuple[float, float] | None = None, batch: tuple[float, float] | None = None):
    root = tmp_path / "root"
    _record(root, "ingest", ingest)
    before = tmp_path / "before.json"
    before.write_text(json.dumps({"ingest_cpu_min": ingest_est}))
    if merge is not None or batch is not None:
        b_cpu, b_est = batch if batch is not None else (100.0, 100.0)
        m_cpu, m_est = merge if merge is not None else (1.0, 1.0)
        plan = {"batches": [{"batch": 0, "genomes": 10, "bases": 1000, "estimate_cpu_min": b_est}],
                "merge_estimate_cpu_min": m_est, "genes_expected": 1000}
        pp = root / est.genecall.plan_record_rel(ST)
        pp.parent.mkdir(parents=True, exist_ok=True)
        pp.write_text(json.dumps(plan))
        _record(root, f"genecall.{est.genecall.batch_name(0)}", b_cpu)
        _record(root, "genecall", m_cpu)
    # The CRISPR estimates read the stage's stored bases; no CRISPR record exists here, so they are never judged.
    monkeypatch.setattr(est.bs, "stage_bases", lambda root, st: 0)
    monkeypatch.setattr(est.bs, "crispr_flank_bases", lambda root, st: 0)
    # `rates` names the CRISPR steps by their producer, which asks the MinCED binary for its version. The guard's
    # logic does not depend on it, and GitHub Actions R1 has no MinCED (card E-12): a fixed version stands in.
    monkeypatch.setattr(est.crispr, "minced_version", lambda: "0.0.0")
    a = argparse.Namespace(root=str(root), stage=ST, before=str(before), objects=None)
    return est, a


def _run(est, a, capsys) -> tuple[int, dict]:
    rc = est.rates(a)
    return rc, json.loads(capsys.readouterr().out)


def test_dec56_accepted_ingest_does_not_stop(tmp_path, monkeypatch, capsys):
    est, a = _rates(_load(), monkeypatch, tmp_path, ingest=45.803, ingest_est=38.057)
    rc, doc = _run(est, a, capsys)
    assert rc == 0
    assert doc["over_20pct"] == []
    assert doc["accepted_over_20pct"] == {"ingest": "DEC-56"}
    (row,) = [r for r in doc["rows"] if r["step"] == "ingest"]
    assert row["accepted"] == "DEC-56"
    assert row["ratio"] == pytest.approx(1.204)
    assert row["over_20pct"] is False


@pytest.mark.parametrize("ingest,ingest_est", [(45.804, 38.057), (45.802, 38.057), (45.803, 38.058),
                                               (45.803, 38.056), (50.0, 38.057)])
def test_dec56_other_ingest_values_still_stop(tmp_path, monkeypatch, capsys, ingest, ingest_est):
    est, a = _rates(_load(), monkeypatch, tmp_path, ingest=ingest, ingest_est=ingest_est)
    rc, doc = _run(est, a, capsys)
    assert rc == 3
    assert doc["over_20pct"] == ["ingest"]
    assert doc["accepted_over_20pct"] == {}


def test_dec56_another_step_at_the_same_ratio_stops(tmp_path, monkeypatch, capsys):
    est, a = _rates(_load(), monkeypatch, tmp_path, ingest=45.803, ingest_est=38.057, merge=(45.803, 38.057))
    rc, doc = _run(est, a, capsys)
    assert rc == 3
    assert doc["over_20pct"] == ["genecall"]
    assert doc["accepted_over_20pct"] == {"ingest": "DEC-56"}


def test_dec56_genecall_batch_overrun_still_returns_4(tmp_path, monkeypatch, capsys):
    est, a = _rates(_load(), monkeypatch, tmp_path, ingest=45.803, ingest_est=38.057, batch=(130.0, 100.0))
    rc, doc = _run(est, a, capsys)
    assert rc == 4
    assert doc["over_20pct"] == [f"genecall.{est.genecall.batch_name(0)}"]
    assert doc["accepted_over_20pct"] == {"ingest": "DEC-56"}


def test_dec56_accepts_only_stage_c5(tmp_path, monkeypatch, capsys):
    est, a = _rates(_load(), monkeypatch, tmp_path, ingest=45.803, ingest_est=38.057)
    src = Path(a.root) / "runs" / f"stage={ST}"
    src.rename(src.with_name("stage=C6"))
    a.stage = "C6"
    rc, doc = _run(est, a, capsys)
    assert rc == 3
    assert doc["over_20pct"] == ["ingest"]
