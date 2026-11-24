"""Ring R1: card B-68 — refit shared-index `.cut` (`mmseqs createindex`) estimates in `train_est.py`.

Oracles (*measured* on wilkcomp02, read-only from `/tmp/if2453/r10/root`, 18 threads):
- `r26.11.1`: `nr50.15` (`10,064,063` targets, `2,852,219,962` residues, `-s 2`, `34.457` CPU-min)
  and `nr50.19` (`9,421,654` targets, `2,698,588,280` residues, `-s 3`, `42.162` CPU-min).
- `r26.11.2`: `nr50.15` (`10,927,842` targets, `3,029,785,935` residues, `-s 2`, `36.452` CPU-min)
  and `nr50.19` (`10,235,527` targets, `2,868,261,896` residues, `-s 3`, `44.629` CPU-min).
- `r26.11.3`: `nr50.15` (`24,298,242` targets, `6,476,355,769` residues, `-s 2`, `113.827` CPU-min).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog, provenance
from biobox.l1.genecall import GENE_CALL_VERSION as GCV

pytestmark = pytest.mark.covers("B-68", "B-64", "C-46")

REPO = Path(__file__).resolve().parents[2]
POD_DIR = REPO / "scripts" / "stress" / "build_stage_pod"
TRAIN = POD_DIR / "train_est.py"
PROFILE = REPO / "profiles" / "wilkcomp02.yaml"


def _te():
    spec = importlib.util.spec_from_file_location("train_est_b68", TRAIN)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["train_est_b68"] = mod
    spec.loader.exec_module(mod)
    return mod


# 5 measured 18-thread shared-index `createindex` `.cut` records:
#   (release, call, sensitivity, targets, residues, n_parts, measured_cpu_min)
MEASURED_CREATEINDEX = [
    ("r26.11.1", "cluster.nr50.15.prefilter", "2", 10_064_063, 2_852_219_962, 24, 34.457),
    ("r26.11.1", "cluster.nr50.19.prefilter", "3", 9_421_654, 2_698_588_280, 24, 42.162),
    ("r26.11.2", "cluster.nr50.15.prefilter", "2", 10_927_842, 3_029_785_935, 3, 36.452),
    ("r26.11.2", "cluster.nr50.19.prefilter", "3", 10_235_527, 2_868_261_896, 6, 44.629),
    ("r26.11.3", "cluster.nr50.15.prefilter", "2", 24_298_242, 6_476_355_769, 9, 113.827),
]


def _write_cut_and_record(
    root: Path,
    release: str,
    call: str,
    *,
    targets: int,
    residues: int,
    n_parts: int,
    measured_cpu_min: float,
    sensitivity: str | None = None,
    shared_index: bool = True,
    include_index_doc: bool = True,
) -> None:
    pdir = root / "l2" / "cluster_plan" / f"release={release}"
    rdir = root / "runs" / f"stage={release}"
    pdir.mkdir(parents=True, exist_ok=True)
    rdir.mkdir(parents=True, exist_ok=True)

    parts = []
    for p in range(n_parts):
        pk = targets // n_parts + (targets % n_parts if p == n_parts - 1 else 0)
        pu = residues // n_parts + (residues % n_parts if p == n_parts - 1 else 0)
        parts.append({"first": p, "last": p, "keys": pk, "units": pu, "cpu_min": 150.0})

    cut_doc: dict = {
        "format": 1,
        "release": release,
        "step": f"{call}.cut",
        "cpu_min": 150.0 * n_parts,
        "overhead_cpu_min": 0.8 if shared_index else 0.0,
        "parts": parts,
    }
    if shared_index and n_parts > 1:
        cut_doc["index_cpu_min"] = round(targets * 1e-6, 6)
        if include_index_doc:
            argv = ["createindex", f"parts/{call}/target", f"parts/{call}/index-tmp"]
            if sensitivity is not None:
                argv.extend(["-s", sensitivity])
            argv.extend(["-k", "0", "--spaced-kmer-mode", "1"])
            cut_doc["index"] = {"path": f"parts/{call}/target.idx", "id": "0" * 64, "argv": argv}

    (pdir / f"{call}.cut.json").write_text(json.dumps(cut_doc))
    (rdir / f"{call}.cut.json").write_text(json.dumps({
        "stage": release,
        "step": f"{call}.cut",
        "exit_code": 0,
        "cpu_min": measured_cpu_min,
        "wall_s": 651.4,
        "peak_rss_bytes": 55_980_000_000,
    }))


def test_b68_check_rates_r26_11_3_stale_before_json_uses_cut_json(tmp_path):
    """Live r26.11.3 stop regression: even when before.json holds stale pre-B-68 steps
    (as on an A-15 resume), check_rates prices shared-index .cut steps from <call>.cut.json."""
    root = tmp_path / "root"
    _write_cut_and_record(
        root, "r26.11.3", "cluster.nr50.15.prefilter",
        targets=24_298_242, residues=6_476_355_769, n_parts=9,
        measured_cpu_min=113.827, sensitivity="2",
    )
    est = tmp_path / "before.json"
    est.write_text(json.dumps({
        "projection": {
            "factors": {},
            "steps": {
                "cluster.nr50.15.prefilter.cut": 64.562,
                "cluster.nr50.19.prefilter.cut": 79.05,
            },
        }
    }))

    te = _te()
    doc = te.check_rates(root, "r26.11.3", est)
    assert doc["exit_code"] == 0, doc
    assert doc["over_20pct"] == [], doc["over_20pct"]
    assert doc["outside_20pct"] == [], doc["outside_20pct"]
    assert doc["worst_over"]["step"] == "cluster.nr50.15.prefilter.cut"
    assert 0.90 <= doc["worst_over"]["ratio"] <= 1.10, doc["worst_over"]
    assert "cut createindex (-s 2, 24,298,242 targets)" == doc["worst_over"]["basis"]


@pytest.mark.parametrize(
    "release,call,sens,targets,residues,n_parts,measured",
    MEASURED_CREATEINDEX,
    ids=[f"{r}-{c}" for r, c, *_ in MEASURED_CREATEINDEX],
)
def test_b68_all_five_measured_createindex_records_within_10_percent(
    tmp_path, release, call, sens, targets, residues, n_parts, measured
):
    """All 5 measured 18-thread createindex .cut records across r26.11.1, r26.11.2, and r26.11.3
    sit within ±10% (0.90–1.10x) of createindex_cpu_min and pass check_rates with stale before.json."""
    te = _te()
    est_val = te.createindex_cpu_min(residues, targets, sens)
    ratio = measured / est_val
    assert 0.90 <= ratio <= 1.10, (release, call, measured, est_val, ratio)

    root = tmp_path / "root"
    _write_cut_and_record(
        root, release, call,
        targets=targets, residues=residues, n_parts=n_parts,
        measured_cpu_min=measured, sensitivity=sens,
    )
    est = tmp_path / "before.json"
    est.write_text(json.dumps({"projection": {"factors": {}, "steps": {f"{call}.cut": 50.0}}}))
    doc = te.check_rates(root, release, est)
    assert doc["exit_code"] == 0, doc
    assert doc["over_20pct"] == []
    assert doc["outside_20pct"] == []
    assert 0.90 <= doc["worst_over"]["ratio"] <= 1.10


def test_b68_createindex_cpu_min_edge_cases_and_fallback():
    te = _te()
    assert te.CREATEINDEX_T0 == 10_000_000
    assert te.CREATEINDEX_GROWTH == 0.48
    assert te.CREATEINDEX_RATE_BY_SENS == {"2": 11.80, "3": 15.75}
    assert te.createindex_cpu_min(0, 10_000_000, "2") == 0.0
    assert te.createindex_cpu_min(1_000_000_000, 0, "2") == 0.0
    assert te.createindex_cpu_min(-1, 10_000_000, "2") == 0.0
    # At T0 and 1e9 residues, equals the base rate
    assert te.createindex_cpu_min(1_000_000_000, 10_000_000, "2") == pytest.approx(11.80, rel=1e-6)
    assert te.createindex_cpu_min(1_000_000_000, 10_000_000, "3") == pytest.approx(15.75, rel=1e-6)
    # Unknown sensitivity falls back cleanly to a positive finite estimate
    other = te.createindex_cpu_min(1_000_000_000, 10_000_000, "4")
    assert other > 0.0


def test_b68_check_rates_fallback_sensitivity_when_index_argv_omits_s(tmp_path):
    """When <call>.cut.json has index_cpu_min without index.argv -s, infer -s 3 for .19. and -s 2 otherwise."""
    root = tmp_path / "root"
    _write_cut_and_record(
        root, "r26.11.2", "cluster.nr50.19.prefilter",
        targets=10_235_527, residues=2_868_261_896, n_parts=6,
        measured_cpu_min=44.629, sensitivity=None, shared_index=True, include_index_doc=False,
    )
    est = tmp_path / "before.json"
    est.write_text(json.dumps({"projection": {"factors": {}, "steps": {"cluster.nr50.19.prefilter.cut": 10.0}}}))
    te = _te()
    doc = te.check_rates(root, "r26.11.2", est)
    assert doc["exit_code"] == 0, doc
    assert doc["outside_20pct"] == []
    assert doc["worst_over"]["basis"] == "cut createindex (-s 3, 10,235,527 targets)"


def test_b68_non_indexed_cut_steps_keep_existing_estimate_path(tmp_path):
    """Non-indexed .cut steps (no index / index_cpu_min, or len(parts) == 1) fall through to steps[name]."""
    root = tmp_path / "root"
    # Multi-part align.cut with no shared index
    _write_cut_and_record(
        root, "r26.11.3", "cluster.nr50.16.align",
        targets=24_298_242, residues=6_476_355_769, n_parts=9,
        measured_cpu_min=1.2, shared_index=False,
    )
    # Multi-part nr90.11.prefilter.cut with no shared index
    _write_cut_and_record(
        root, "r26.11.3", "cluster.nr90.11.prefilter",
        targets=60_000_000, residues=18_000_000_000, n_parts=1,
        measured_cpu_min=2.5, shared_index=False,
    )
    # Single-part nr50.15.prefilter.cut (len(parts) == 1 => no shared index built)
    _write_cut_and_record(
        root, "r26.11.3", "cluster.nr50.15.prefilter",
        targets=2_000_000, residues=600_000_000, n_parts=1,
        measured_cpu_min=0.8, shared_index=False,
    )
    est = tmp_path / "before.json"
    est.write_text(json.dumps({
        "projection": {
            "factors": {},
            "steps": {
                "cluster.nr50.16.align.cut": 1.1,
                "cluster.nr90.11.prefilter.cut": 2.4,
                "cluster.nr50.15.prefilter.cut": 0.75,
            },
        }
    }))

    out_path = tmp_path / "rates.json"
    te = _te()
    doc = te.check_rates(root, "r26.11.3", est, out_path=out_path)
    assert doc["exit_code"] == 0
    assert doc["judged"] == 0
    rows = {r["step"]: r for r in json.loads(out_path.read_text())["rows"]}
    for step, expected_est in (
        ("cluster.nr50.16.align.cut", 1.1),
        ("cluster.nr90.11.prefilter.cut", 2.4),
        ("cluster.nr50.15.prefilter.cut", 0.75),
    ):
        assert rows[step]["estimate"] == expected_est
        assert rows[step]["basis"] == "earlier record x ratio"


def test_b68_rates_cli_runs_stdlib_only_under_python_I(tmp_path):
    """Host invocation `python3 -I train_est.py rates ...` must succeed with stdlib only."""
    root = tmp_path / "root"
    _write_cut_and_record(
        root, "r26.11.3", "cluster.nr50.15.prefilter",
        targets=24_298_242, residues=6_476_355_769, n_parts=9,
        measured_cpu_min=113.827, sensitivity="2",
    )
    est = tmp_path / "before.json"
    est.write_text(json.dumps({"projection": {"factors": {}, "steps": {"cluster.nr50.15.prefilter.cut": 64.562}}}))
    res = subprocess.run(
        [sys.executable, "-I", str(TRAIN), "rates", "--root", str(root), "--release", "r26.11.3", "--est", str(est)],
        capture_output=True, text=True, check=False,
    )
    assert res.returncode == 0, res.stderr
    doc = json.loads(res.stdout)
    assert doc["over_20pct"] == []
    assert doc["outside_20pct"] == []
    assert 0.90 <= doc["worst_over"]["ratio"] <= 1.10


def _write_stage_stubs(root: Path, stage: str) -> None:
    ps_rel = f"l0seq/protein_seq/gene_call_version={GCV}/stage={stage}/protein_prefix=00/part-0.parquet"
    fp_rel = f"l1/feature_protein/gene_call_version={GCV}/stage={stage}/part-0.parquet"
    for rel, col in ((ps_rel, "protein_uid"), (fp_rel, "feature_uid")):
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.table({col: [b"\x00" * 16]}), p)
    rec = {
        "stage": stage, "step": "genecall", "exit_code": 0,
        "outputs": [
            {"path": ps_rel, "bytes": (root / ps_rel).stat().st_size},
            {"path": fp_rel, "bytes": (root / fp_rel).stat().st_size},
        ],
    }
    rd = root / "runs" / f"stage={stage}"
    rd.mkdir(parents=True, exist_ok=True)
    (rd / "genecall.json").write_text(json.dumps(rec))


def _run_before(
    tmp_path: Path,
    monkeypatch,
    capsys,
    *,
    root: Path,
    earlier: str,
    release: str,
    stages: list[str],
    earlier_proteins: int,
    release_proteins: int,
    parts: int | None = None,
) -> dict:
    if shutil.which("mmseqs") is None:
        monkeypatch.setattr(provenance, "tool_versions", lambda *names: {n: "stub" for n in names})

    import duckdb

    counts = iter([release_proteins, earlier_proteins])

    class _FakeCon:
        def execute(self, _sql):
            val = next(counts)
            return type("_Res", (), {"fetchone": lambda self: (val,)})()

        def close(self):
            pass

    monkeypatch.setattr(duckdb, "connect", lambda *a, **kw: _FakeCon())

    te = _te()
    ns = argparse.Namespace(
        root=str(root), release=release, stages=",".join(stages),
        earlier=earlier, parts=parts, profile=str(PROFILE),
        jobs=2, cpus=18, tmp=str(tmp_path / "tmp"), duck_threads=1,
    )
    assert te.before(ns) == 0
    return json.loads(capsys.readouterr().out)


def test_b68_before_projects_shared_index_cut_including_k_old_1_to_k_gt_1(tmp_path, monkeypatch, capsys):
    """`before` projects shared-index `.cut` steps via `createindex_cpu_min` when projected k > 1,
    even when the earlier release had k_old == 1 (no shared index built), while keeping k == 1 and
    non-indexed `.cut` steps on their small non-indexed estimate."""
    root = tmp_path / "root"
    for s in ("S1", "S2"):
        _write_stage_stubs(root, s)

    earlier = "r26.11.0"
    catalog.define_release(root, earlier, ["S1"])
    pdir = root / "l2" / "cluster_plan" / f"release={earlier}"
    rdir = root / "runs" / f"stage={earlier}"
    pdir.mkdir(parents=True, exist_ok=True)
    rdir.mkdir(parents=True, exist_ok=True)
    (pdir / "cluster.plan.json").write_text(json.dumps({
        "format": 1, "release": earlier, "stages": ["S1"], "proteins": 15_000_000,
        "cut": {"part_cpu_min": 240.0, "one_process_cpu_min_max": 20000.0},
    }))

    # Earlier release had k_old == 1 for nr50.15.prefilter (small non-indexed .cut = 0.4 CPU-min),
    # k_old == 2 for nr50.19.prefilter (shared-index .cut = 15.0 CPU-min), and k_old == 2 for nr50.16.align.
    u15, t15 = 1_500_000_000, 5_000_000
    u19, t19 = 1_400_000_000, 4_700_000
    u16, t16 = 1_600_000_000, 5_000_000
    for call, k_old, u, t, oh, cut_meas in (
        ("cluster.nr50.15.prefilter", 1, u15, t15, 0.4, 0.4),
        ("cluster.nr50.19.prefilter", 2, u19, t19, 0.4, 15.0),
        ("cluster.nr50.16.align", 2, u16, t16, 0.0, 0.6),
    ):
        parts = [{"first": i, "last": i, "keys": t // k_old, "units": u // k_old, "cpu_min": 180.0} for i in range(k_old)]
        (pdir / f"{call}.cut.json").write_text(json.dumps({
            "format": 1, "release": earlier, "step": f"{call}.cut",
            "cpu_min": 180.0 * k_old, "overhead_cpu_min": oh, "parts": parts,
        }))
        for i in range(k_old):
            sname = f"{call}.part-{i:05d}" if k_old > 1 else call
            (rdir / f"{sname}.json").write_text(json.dumps({"stage": earlier, "step": sname, "cpu_min": 180.0}))
        (rdir / f"{call}.cut.json").write_text(json.dumps({"stage": earlier, "step": f"{call}.cut", "cpu_min": cut_meas}))

    # Grow 4x to 60M proteins: both nr50.15 and nr50.19 project to k > 1
    doc_grown = _run_before(
        tmp_path, monkeypatch, capsys,
        root=root, earlier=earlier, release="r26.11.9", stages=["S1", "S2"],
        earlier_proteins=15_000_000, release_proteins=60_000_000,
    )
    te = _te()
    steps = doc_grown["projection"]["steps"]
    r = 60_000_000 / 15_000_000
    assert steps["cluster.nr50.15.prefilter.cut"] == round(te.createindex_cpu_min(round(u15 * r), round(t15 * r), "2"), 3)
    assert steps["cluster.nr50.19.prefilter.cut"] == round(te.createindex_cpu_min(round(u19 * r), round(t19 * r), "3"), 3)
    # Non-indexed align.cut keeps its small linear estimate
    assert steps["cluster.nr50.16.align.cut"] == round(0.6 * r, 3)

    # With --parts 1 (k == 1), no shared index is built so nr50.15.prefilter.cut stays on its small estimate
    doc_one = _run_before(
        tmp_path, monkeypatch, capsys,
        root=root, earlier=earlier, release="r26.11.9", stages=["S1", "S2"],
        earlier_proteins=15_000_000, release_proteins=60_000_000, parts=1,
    )
    steps_one = doc_one["projection"]["steps"]
    assert steps_one["cluster.nr50.15.prefilter.cut"] == round(0.4 * r, 3)
    assert steps_one["cluster.nr50.19.prefilter.cut"] < 5.0
