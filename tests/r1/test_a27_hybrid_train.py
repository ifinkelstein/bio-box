"""Ring R1: card A-27 — release-train driver updates for hybrid clustering (`DEC-77`).

Covers:
1. `train_est.py`: `--cluster-method` (default `hybrid-b1`), recording `cluster_method` in `before.json`,
   filtering `calib["calls"]` and `calib["others"]` to active `l2.calls(mode)` across a `cascade -> hybrid-b1`
   transition, `B-72` level-split finish steps (`cluster.finish.nr90`, `cluster.finish.nr50`) in `before` and
   `rates`, and parity between `train_est.createindex_cpu_min` and `biobox.l2.cluster.createindex_cpu_min`
   (while keeping `python3 -I train_est.py rates` stdlib-only).
2. `r50_watch.py`: `DEC-77` one-process exemption on shared-index `.cut` (`createindex`) steps for both
   published records (`cpu_budget_reason == "exempt: DEC-77"` backed by `<call>.cut.json`) and live
   `release cluster step|redo` processes (reading `one_process_cpu_min_max` from `cluster.plan.json` while
   `<call>.cut.json` is not yet written), while strictly capping all non-exempt steps (`part-*`,
   non-shared-index `.cut`, `cluster.finish.nr90`/`nr50`, `cluster`) at `--cpu-max` (280 CPU-min).
3. `release_train.sh` & `train_plan.py`: defaults (`CLUSTER_METHOD=hybrid-b1`, `SCRATCH_MAX_GB=350`,
   `TMP_MIN_GB=100`), passing `--cluster-method` to `train_plan.py`, `train_est.py before`, and `build-stage`,
   backwards-compatible resume of a pre-A-27 train without `cluster_method`, and B-49 reuse detection under
   `--cluster-method`.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog, provenance
from biobox.l1.genecall import GENE_CALL_VERSION as GCV
from biobox.l2 import cluster as l2

pytestmark = pytest.mark.covers("A-27", "DEC-77", "DEC-64")

REPO = Path(__file__).resolve().parents[2]
POD_DIR = REPO / "scripts" / "stress" / "build_stage_pod"
TRAIN_EST = POD_DIR / "train_est.py"
TRAIN_PLAN = POD_DIR / "train_plan.py"
WATCH = POD_DIR / "r50_watch.py"
RELEASE_TRAIN = POD_DIR / "release_train.sh"
PROFILE = REPO / "profiles" / "wilkcomp02.yaml"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _te():
    return _load("train_est_a27", TRAIN_EST)


def _watch():
    return _load("r50_watch_a27", WATCH)


# --- 1. train_est.py: hybrid mode filtering, level-split finish, createindex parity -----------------------------


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


def _make_cascade_earlier_release(root: Path, earlier: str = "r26.11.3", *, split_finish: bool = False) -> None:
    """Create a realistic cascade earlier release (like r26.11.3) with dropped and kept calls."""
    for s in ("S1", "S2"):
        _write_stage_stubs(root, s)
    catalog.define_release(root, earlier, ["S1"])
    pdir = root / "l2" / "cluster_plan" / f"release={earlier}"
    rdir = root / "runs" / f"stage={earlier}"
    pdir.mkdir(parents=True, exist_ok=True)
    rdir.mkdir(parents=True, exist_ok=True)
    (pdir / "cluster.plan.json").write_text(json.dumps({
        "format": 1, "release": earlier, "stages": ["S1"], "proteins": 100_000_000,
        "cut": {"part_cpu_min": 240.0, "one_process_cpu_min_max": 20000.0},
    }))

    # Split calls: nr90.11/12 (dropped in hybrid-a/b1), nr50.15/16 (kept in hybrid-a/b1),
    # nr50.19/20 (dropped in hybrid-b1, kept in hybrid-a/cascade)
    split_calls = [
        ("cluster.nr90.11.prefilter", 1, 30_000_000_000, 70_000_000, 350.0, 5000.0, 2.0),
        ("cluster.nr90.12.align", 8, 25_000_000_000, 70_000_000, 0.0, 3200.0, 1.5),
        ("cluster.nr50.15.prefilter", 9, 6_500_000_000, 24_000_000, 0.8, 1200.0, 114.0),
        ("cluster.nr50.16.align", 4, 6_000_000_000, 24_000_000, 0.0, 600.0, 0.9),
        ("cluster.nr50.19.prefilter", 20, 6_000_000_000, 22_000_000, 0.8, 4000.0, 145.0),
        ("cluster.nr50.20.align", 6, 5_500_000_000, 22_000_000, 0.0, 900.0, 0.8),
    ]
    for call, k, u, t, oh, total_meas, cut_meas in split_calls:
        parts = [
            {"first": i, "last": i, "keys": t // k, "units": u // k, "cpu_min": total_meas / k}
            for i in range(k)
        ]
        (pdir / f"{call}.cut.json").write_text(json.dumps({
            "format": 1, "release": earlier, "step": f"{call}.cut",
            "cpu_min": total_meas, "overhead_cpu_min": oh, "parts": parts,
        }))
        for i in range(k):
            sname = f"{call}.part-{i:05d}" if k > 1 else call
            (rdir / f"{sname}.json").write_text(json.dumps({
                "stage": earlier, "step": sname, "cpu_min": total_meas / k, "peak_rss_bytes": 40_000_000_000,
            }))
        (rdir / f"{call}.cut.json").write_text(json.dumps({
            "stage": earlier, "step": f"{call}.cut", "cpu_min": cut_meas,
        }))

    # Non-split calls (kept and dropped)
    others = {
        "cluster.plan": 80.0,
        "cluster.nr90.00.kmermatcher": 60.0,
        "cluster.nr90.10.createsubdb": 15.0,    # dropped in hybrid-a/b1
        "cluster.nr90.13.clust": 25.0,          # dropped in hybrid-a/b1
        "cluster.nr90.14.mergeclusters": 10.0,  # dropped in hybrid-a/b1
        "cluster.nr50.00.kmermatcher": 45.0,
        "cluster.nr50.18.createsubdb": 12.0,    # dropped in hybrid-b1
        "cluster.nr50.21.clust": 20.0,          # dropped in hybrid-b1
        "tokens": 30.0,
    }
    if split_finish:
        others["cluster.finish.nr90"] = 85.0
        others["cluster.finish.nr50"] = 80.0
        others["cluster"] = 8.0
    else:
        others["cluster"] = 170.0

    for step, cpu in others.items():
        (rdir / f"{step}.json").write_text(json.dumps({
            "stage": earlier, "step": step, "cpu_min": cpu,
            "argv": ["--threads", "18"],
        }))


def _run_before_cli(
    tmp_path: Path,
    monkeypatch,
    capsys,
    *,
    root: Path,
    earlier: str = "r26.11.3",
    release: str = "r26.11.4",
    stages: list[str] = ("S1", "S2"),
    earlier_proteins: int = 100_000_000,
    release_proteins: int = 200_000_000,
    cluster_method: str | None = None,
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
    argv = [
        "before",
        "--root", str(root),
        "--release", release,
        "--stages", ",".join(stages),
        "--earlier", earlier,
        "--profile", str(PROFILE),
        "--jobs", "2",
        "--cpus", "18",
        "--tmp", str(tmp_path / "tmp"),
        "--duck-threads", "1",
    ]
    if cluster_method is not None:
        argv.extend(["--cluster-method", cluster_method])
    assert te.main(argv) == 0
    return json.loads(capsys.readouterr().out)


def test_train_est_before_filters_dropped_calls_for_hybrid_b1_and_records_method(tmp_path, monkeypatch, capsys):
    """Transitioning cascade r26.11.3 -> hybrid-b1 r26.11.4 must drop nr90.11/12 and nr50.19/20
    (and their .cut / non-split cascade steps) so they do not inflate projected CPU."""
    root = tmp_path / "root"
    _make_cascade_earlier_release(root)

    # Default CLI --cluster-method is hybrid-b1
    doc_b1 = _run_before_cli(tmp_path, monkeypatch, capsys, root=root)
    assert doc_b1["cluster_method"] == "hybrid-b1"
    b1_calls = {c["call"] for c in doc_b1["projection"]["calls"]}
    b1_steps = set(doc_b1["projection"]["steps"])
    b1_factors = set(doc_b1["projection"]["factors"])

    # Kept in hybrid-b1
    assert {"cluster.nr50.15.prefilter", "cluster.nr50.16.align"} <= b1_calls
    assert "cluster.nr50.15.prefilter.cut" in b1_steps
    assert "cluster.nr90.00.kmermatcher" in b1_steps
    assert "tokens" in b1_steps

    # Dropped in hybrid-b1
    dropped = {
        "cluster.nr90.11.prefilter", "cluster.nr90.12.align",
        "cluster.nr50.19.prefilter", "cluster.nr50.20.align",
    }
    assert b1_calls.isdisjoint(dropped), b1_calls
    assert b1_factors.isdisjoint(dropped), b1_factors
    for s in (
        "cluster.nr90.11.prefilter.cut", "cluster.nr90.12.align.cut",
        "cluster.nr90.10.createsubdb", "cluster.nr90.13.clust", "cluster.nr90.14.mergeclusters",
        "cluster.nr50.19.prefilter.cut", "cluster.nr50.20.align.cut",
        "cluster.nr50.18.createsubdb", "cluster.nr50.21.clust",
    ):
        assert s not in b1_steps, f"dropped step {s} leaked into projection.steps"

    # Compare against cascade and hybrid-a: hybrid-b1 < hybrid-a < cascade
    doc_a = _run_before_cli(tmp_path, monkeypatch, capsys, root=root, cluster_method="hybrid-a")
    doc_cas = _run_before_cli(tmp_path, monkeypatch, capsys, root=root, cluster_method="cascade")
    assert doc_a["cluster_method"] == "hybrid-a"
    assert doc_cas["cluster_method"] == "cascade"
    a_calls = {c["call"] for c in doc_a["projection"]["calls"]}
    assert "cluster.nr50.19.prefilter" in a_calls and "cluster.nr90.11.prefilter" not in a_calls
    assert doc_b1["projection"]["cluster_cpu_min"] < doc_a["projection"]["cluster_cpu_min"]
    assert doc_a["projection"]["cluster_cpu_min"] < doc_cas["projection"]["cluster_cpu_min"]


def test_train_est_before_and_rates_handle_level_split_finish_transition(tmp_path, monkeypatch, capsys):
    """Transitioning from an unsplit-finish earlier release (`cluster` = 170 CPU-min) to a B-72
    level-split release (`cluster.finish.nr90`, `cluster.finish.nr50`, lightweight `cluster`)
    projects nr90/nr50 cleanly and passes check_rates without false stops."""
    root = tmp_path / "root"
    _make_cascade_earlier_release(root, split_finish=False)

    # 2x growth: old cluster = 170 -> total 340 -> 170 each for nr90 and nr50
    doc = _run_before_cli(tmp_path, monkeypatch, capsys, root=root, cluster_method="hybrid-b1")
    steps = doc["projection"]["steps"]
    assert steps["cluster.finish.nr90"] == pytest.approx(170.0, rel=1e-3)
    assert steps["cluster.finish.nr50"] == pytest.approx(170.0, rel=1e-3)
    assert "cluster" not in steps  # unsplit earlier `cluster` must not judge lightweight B-72 `cluster`

    # Now publish B-72 level-split finish records on r26.11.4 and run check_rates
    rdir = root / "runs" / "stage=r26.11.4"
    rdir.mkdir(parents=True, exist_ok=True)
    for name, cpu in (("cluster.finish.nr90", 185.0), ("cluster.finish.nr50", 160.0), ("cluster", 12.0)):
        (rdir / f"{name}.json").write_text(json.dumps({
            "stage": "r26.11.4", "step": name, "cpu_min": cpu, "wall_s": 60.0, "peak_rss_bytes": 10_000_000_000,
        }))

    est_path = tmp_path / "before.json"
    est_path.write_text(json.dumps(doc))
    te = _te()
    res = te.check_rates(root, "r26.11.4", est_path)
    assert res["exit_code"] == 0, res
    assert res["over_20pct"] == []
    assert res["outside_20pct"] == []

    # Also test legacy/stale before.json that only has unsplit `steps["cluster"] = 340.0`
    legacy_est = tmp_path / "legacy_before.json"
    legacy_est.write_text(json.dumps({"projection": {"factors": {}, "steps": {"cluster": 340.0}}}))
    res_legacy = te.check_rates(root, "r26.11.4", legacy_est)
    assert res_legacy["exit_code"] == 0, res_legacy
    assert res_legacy["over_20pct"] == []
    assert res_legacy["outside_20pct"] == []


def test_train_est_before_preserves_direct_scaling_when_earlier_already_level_split(tmp_path, monkeypatch, capsys):
    """When --earlier already has cluster.finish.nr90, cluster.finish.nr50, and lightweight cluster,
    each scales directly by r."""
    root = tmp_path / "root"
    _make_cascade_earlier_release(root, split_finish=True)
    doc = _run_before_cli(tmp_path, monkeypatch, capsys, root=root, cluster_method="hybrid-b1")
    steps = doc["projection"]["steps"]
    assert steps["cluster.finish.nr90"] == pytest.approx(170.0, rel=1e-3)
    assert steps["cluster.finish.nr50"] == pytest.approx(160.0, rel=1e-3)
    assert steps["cluster"] == pytest.approx(16.0, rel=1e-3)


def test_train_est_createindex_cpu_min_matches_l2_cluster_and_rates_stays_stdlib_only(tmp_path):
    """train_est.createindex_cpu_min must stay identical to biobox.l2.cluster.createindex_cpu_min,
    and `python3 -I train_est.py rates` must run stdlib-only on the host."""
    te = _te()
    for res, tgt, sens in (
        (6_476_355_769, 24_298_242, "2"),
        (2_868_261_896, 10_235_527, "3"),
        (18_500_000_000, 70_000_000, "2"),
        (17_000_000_000, 66_000_000, "3"),
        (1_000_000_000, 10_000_000, "4"),
        (0, 10_000_000, "2"),
    ):
        assert te.createindex_cpu_min(res, tgt, sens) == pytest.approx(
            l2.createindex_cpu_min(res, tgt, sens), rel=1e-12
        )

    # Host stdlib-only check
    root = tmp_path / "root"
    rdir = root / "runs" / "stage=r26.11.4"
    rdir.mkdir(parents=True)
    (rdir / "cluster.finish.nr90.json").write_text(json.dumps({"step": "cluster.finish.nr90", "cpu_min": 180.0}))
    (rdir / "cluster.json").write_text(json.dumps({"step": "cluster", "cpu_min": 10.0}))
    est = tmp_path / "before.json"
    est.write_text(json.dumps({"projection": {"factors": {}, "steps": {"cluster": 340.0}}}))
    proc = subprocess.run(
        [sys.executable, "-I", str(TRAIN_EST), "rates", "--root", str(root), "--release", "r26.11.4", "--est", str(est)],
        capture_output=True, text=True, check=False,
    )
    assert proc.returncode == 0, proc.stderr + proc.stdout
    assert json.loads(proc.stdout)["over_20pct"] == []


# --- 2. r50_watch.py: DEC-77 createindex .cut exemption ---------------------------------------------------------


REL_W = "r26.11.4"
CUT_15 = "cluster.nr50.15.prefilter.cut"
CUT_19 = "cluster.nr50.19.prefilter.cut"
DEC77_REASON = "exempt: DEC-77"


def _setup_watch_env(tmp_path: Path) -> tuple[Path, Path, Path, Path, float]:
    root = tmp_path / "root"
    scratch = tmp_path / "scratch"
    state = tmp_path / "state"
    root.mkdir()
    scratch.mkdir()
    state.mkdir()
    since = tmp_path / "since"
    since.write_text("mark\n")
    t0 = since.stat().st_mtime
    os.utime(since, (t0 - 60, t0 - 60))
    return root, scratch, since, state, t0


def _write_plan(root: Path, rel: str = REL_W, doc: dict | str | None = None) -> Path:
    path = root / "l2" / "cluster_plan" / f"release={rel}" / "cluster.plan.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    if doc is None:
        doc = {
            "format": 2,
            "release": rel,
            "method": "hybrid-b1",
            "cut": {"part_cpu_min": 240.0, "one_process_cpu_min_max": 20000.0},
        }
    path.write_text(doc if isinstance(doc, str) else json.dumps(doc))
    return path


def _write_dec77_cut(root: Path, rel: str = REL_W, step: str = CUT_15, doc: dict | str | None = None) -> Path:
    path = root / "l2" / "cluster_plan" / f"release={rel}" / f"{step}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    if doc is None:
        doc = {
            "parts": [
                {"first": 0, "last": 4, "keys": 35_000_000, "units": 9_000_000_000, "cpu_min": 200.0},
                {"first": 5, "last": 9, "keys": 35_000_000, "units": 9_000_000_000, "cpu_min": 200.0},
            ],
            "index_cpu_min": 70.0,
            "createindex_cpu_min": 560.0,
            "one_process": {
                "decision": "DEC-77",
                "cpu_min": 560.0,
                "cpu_min_max": 20000.0,
                "why": "shared index createindex exceeds part budget",
            },
        }
    path.write_text(doc if isinstance(doc, str) else json.dumps(doc))
    return path


def _write_run_rec(root: Path, rel: str, name: str, doc: dict, *, mtime: float) -> Path:
    path = root / "runs" / f"stage={rel}" / f"{name}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc))
    os.utime(path, (mtime, mtime))
    return path


def _run_watch(monkeypatch, capsys, tmp_path: Path, *, live_rows: list[dict] | None = None, cpu_max: float = 280.0):
    mod = _watch()
    argv = [
        str(WATCH),
        "--root", str(tmp_path / "root"),
        "--scratch", str(tmp_path / "scratch"),
        "--since", str(tmp_path / "since"),
        "--state", str(tmp_path / "state"),
        "--tmp-min-gb", "0",
        "--scratch-max-gb", "1000",
        "--cpu-max", str(cpu_max),
    ]
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(mod, "live", lambda _prefix: list(live_rows or []))
    rc = mod.main()
    return rc, json.loads(capsys.readouterr().out.strip())


@pytest.mark.parametrize("cut_step", [CUT_15, CUT_19])
def test_watch_published_dec77_cut_passes_under_cap_and_fails_above_cap(monkeypatch, capsys, tmp_path, cut_step):
    root, scratch, since, state, t0 = _setup_watch_env(tmp_path)
    _write_dec77_cut(root, REL_W, cut_step)
    _write_run_rec(root, REL_W, cut_step, {"cpu_min": 560.0, "cpu_budget_reason": DEC77_REASON}, mtime=t0)

    rc, out = _run_watch(monkeypatch, capsys, tmp_path)
    assert rc == 0, out
    assert out["reasons"] == []

    # Host stdlib-only invocation (`python3 -I`)
    proc = subprocess.run(
        [
            sys.executable, "-I", str(WATCH),
            "--root", str(root), "--scratch", str(scratch), "--since", str(since), "--state", str(state),
            "--tmp-min-gb", "0", "--scratch-max-gb", "1000", "--cpu-max", "280",
        ],
        capture_output=True, text=True, check=False,
    )
    assert proc.returncode == 0, proc.stderr + proc.stdout
    assert json.loads(proc.stdout)["reasons"] == []

    # Above cpu_min_max must still fail
    _write_run_rec(root, REL_W, cut_step, {"cpu_min": 20000.5, "cpu_budget_reason": DEC77_REASON}, mtime=t0)
    rc, out = _run_watch(monkeypatch, capsys, tmp_path)
    assert rc == 3, out
    assert any("20000" in r for r in out["reasons"]), out["reasons"]


@pytest.mark.parametrize("action", ["step", "redo"])
@pytest.mark.parametrize("cut_step", [CUT_15, CUT_19])
def test_watch_live_dec77_cut_uses_cluster_plan_while_running_and_enforces_cap(
    monkeypatch, capsys, tmp_path, action, cut_step
):
    """While <call>.cut is running, <call>.cut.json is NOT yet on disk; r50_watch.py reads
    one_process_cpu_min_max from cluster.plan.json ONLY for shared-index prefilter .cut steps."""
    root, _, _, _, _ = _setup_watch_env(tmp_path)
    _write_plan(root, REL_W)
    step = f"release cluster {action} release={REL_W} name={cut_step}"

    rc, out = _run_watch(monkeypatch, capsys, tmp_path, live_rows=[{"pid": 555, "step": step, "cpu_min": 560.0}])
    assert rc == 0, out
    assert out["reasons"] == []

    rc, out = _run_watch(monkeypatch, capsys, tmp_path, live_rows=[{"pid": 555, "step": step, "cpu_min": 20000.5}])
    assert rc == 3, out
    assert any("20000" in r for r in out["reasons"]), out["reasons"]


def test_watch_live_dec77_cut_fails_without_plan_cap_or_with_invalid_existing_cut_json(monkeypatch, capsys, tmp_path):
    root, _, _, _, _ = _setup_watch_env(tmp_path)
    step = f"release cluster step release={REL_W} name={CUT_15}"

    # No cluster.plan.json -> fails at 280
    rc, out = _run_watch(monkeypatch, capsys, tmp_path, live_rows=[{"pid": 555, "step": step, "cpu_min": 560.0}])
    assert rc == 3 and any("280" in r for r in out["reasons"]), out

    # Plan without part_cpu_min (forced --parts) -> fails at 280
    _write_plan(root, REL_W, {"format": 2, "cut": {"parts": 4, "one_process_cpu_min_max": 20000.0}})
    rc, out = _run_watch(monkeypatch, capsys, tmp_path, live_rows=[{"pid": 555, "step": step, "cpu_min": 560.0}])
    assert rc == 3 and any("280" in r for r in out["reasons"]), out

    # Valid plan, but <call>.cut.json already exists on disk WITHOUT valid DEC-77 one_process -> fails at 280
    _write_plan(root, REL_W)
    _write_dec77_cut(root, REL_W, CUT_15, {"parts": [{"first": 0}, {"first": 1}]})
    rc, out = _run_watch(monkeypatch, capsys, tmp_path, live_rows=[{"pid": 555, "step": step, "cpu_min": 560.0}])
    assert rc == 3 and any("280" in r for r in out["reasons"]), out


@pytest.mark.parametrize(
    "record_name,live_step",
    [
        ("cluster.nr90.11.prefilter.cut", f"release cluster step release={REL_W} name=cluster.nr90.11.prefilter.cut"),
        ("cluster.nr50.16.align.cut", f"release cluster step release={REL_W} name=cluster.nr50.16.align.cut"),
        ("cluster.nr50.15.prefilter.part-00000", f"release cluster step release={REL_W} name=cluster.nr50.15.prefilter.part-00000"),
        ("cluster.finish.nr90", f"release cluster step release={REL_W} name=cluster.finish.nr90"),
        ("cluster.finish.nr50", f"release cluster step release={REL_W} name=cluster.finish.nr50"),
        ("cluster", f"release cluster finish release={REL_W}"),
        ("cluster.nr50.15.prefilter", f"release cluster step release={REL_W} name=cluster.nr50.15.prefilter"),
    ],
    ids=[
        "non-shared-prefilter-cut",
        "align-cut",
        "split-part",
        "finish-nr90",
        "finish-nr50",
        "cluster-publish",
        "prefilter-join",
    ],
)
def test_watch_non_exempt_steps_stay_strictly_capped_at_280(monkeypatch, capsys, tmp_path, record_name, live_step):
    root, _, _, _, t0 = _setup_watch_env(tmp_path)
    _write_plan(root, REL_W)
    _write_dec77_cut(root, REL_W, CUT_15)

    # Published record fails at 300 > 280 even if claiming exempt: DEC-77
    _write_run_rec(root, REL_W, record_name, {"cpu_min": 300.0, "cpu_budget_reason": DEC77_REASON}, mtime=t0)
    rc, out = _run_watch(monkeypatch, capsys, tmp_path)
    assert rc == 3 and any("280" in r for r in out["reasons"]), out

    # Live process fails at 300 > 280
    (root / "runs" / f"stage={REL_W}" / f"{record_name}.json").unlink()
    rc, out = _run_watch(monkeypatch, capsys, tmp_path, live_rows=[{"pid": 777, "step": live_step, "cpu_min": 300.0}])
    assert rc == 3 and any("280" in r for r in out["reasons"]), out


def test_watch_published_dec77_cut_requires_reason_and_valid_cut_json(monkeypatch, capsys, tmp_path):
    root, _, _, _, t0 = _setup_watch_env(tmp_path)
    _write_dec77_cut(root, REL_W, CUT_15)

    # Missing cpu_budget_reason -> fails at 280
    _write_run_rec(root, REL_W, CUT_15, {"cpu_min": 560.0}, mtime=t0)
    rc, out = _run_watch(monkeypatch, capsys, tmp_path)
    assert rc == 3 and any("280" in r for r in out["reasons"]), out

    # Wrong decision in cut.json (DEC-29 instead of DEC-77) -> fails at 280
    _write_dec77_cut(root, REL_W, CUT_15, {
        "parts": [{"first": 0}, {"first": 1}],
        "one_process": {"decision": "DEC-29", "cpu_min": 560.0, "cpu_min_max": 20000.0},
    })
    _write_run_rec(root, REL_W, CUT_15, {"cpu_min": 560.0, "cpu_budget_reason": DEC77_REASON}, mtime=t0)
    rc, out = _run_watch(monkeypatch, capsys, tmp_path)
    assert rc == 3 and any("280" in r for r in out["reasons"]), out

    # Single-part cut.json (len(parts) == 1) -> fails at 280
    _write_dec77_cut(root, REL_W, CUT_15, {
        "parts": [{"first": 0}],
        "one_process": {"decision": "DEC-77", "cpu_min": 560.0, "cpu_min_max": 20000.0},
    })
    rc, out = _run_watch(monkeypatch, capsys, tmp_path)
    assert rc == 3 and any("280" in r for r in out["reasons"]), out


# --- 3. release_train.sh & train_plan.py: defaults, --cluster-method wiring, backwards-compat resume -----------


def test_release_train_defaults_match_dec77():
    src = RELEASE_TRAIN.read_text()
    for line in (
        "CLUSTER_METHOD=${CLUSTER_METHOD:-hybrid-b1}",
        "SCRATCH_MAX_GB=${SCRATCH_MAX_GB:-350}",
        "TMP_MIN_GB=${TMP_MIN_GB:-100}",
        "CPU_MAX=${CPU_MAX:-280}",
    ):
        assert line in src, f"missing default {line!r} in release_train.sh"


def _make_stub_train_world(tmp_path: Path, *, plan_cluster_method: str | None = "hybrid-b1"):
    pod = tmp_path / "if2453"
    fake = tmp_path / "fakebin"
    bin_dir = tmp_path / "bin"
    tree = tmp_path / "tree"
    root = pod / "r10" / "root"
    work = pod / "c46" / "train"
    host_tmp = pod / "bio-box"
    capture = tmp_path / "capture.log"
    for d in (fake, bin_dir, tree / "scripts", tree / "profiles", root / "manifests", host_tmp):
        d.mkdir(parents=True, exist_ok=True)
    (root / "manifests" / "release=r26.11.3.json").write_text("{}\n")

    plan: dict = {
        "release": "r26.11.4",
        "stages": ["S1", "P1"],
        "new_stages": ["P1"],
        "latest_release": "r26.11.3",
        "nothing_new": False,
        "earlier_manifest_releases": ["r26.11.3"],
        "last_stage": {
            "stage": "P1", "depth": "called", "genomes": "g.json", "objects": "o.json", "crispr": False,
        },
    }
    if plan_cluster_method is not None:
        plan["cluster_method"] = plan_cluster_method
    plan_file = tmp_path / "plan.json"
    plan_file.write_text(json.dumps(plan))

    prof = tree / "profiles" / "wilkcomp02.yaml"
    prof.write_text(
        f"max_threads_day: 40\nmax_threads_night: 80\nmem_scope: 750G\n"
        f"tmp: {host_tmp}\nday_start_hour: 8\nday_end_hour: 20\n"
    )

    stubs = {
        "systemd-run": (
            '#!/bin/sh\nwhile [ $# -gt 0 ]; do\n'
            '  case "$1" in --user|--scope|--quiet|--collect|--unit=*) shift ;; -p) shift 2 ;; *) break ;; esac\n'
            'done\nexec "$@"\n'
        ),
        "nice": '#!/bin/sh\n[ "$1" = -n ] && shift 2\nexec "$@"\n',
        "sleep": "#!/bin/sh\nexit 0\n",
        "systemctl": (
            '#!/bin/sh\ncase "$*" in\n  *is-active*) exit 3 ;;\n'
            '  *show*-p*ActiveState*) echo inactive; exit 0 ;;\nesac\nexit 0\n'
        ),
        "singularity": (
            '#!/bin/sh\n'
            'echo "singularity $*" >> "$CAPTURE_LOG"\n'
            'case "$*" in\n'
            '  *" release verify "*) echo "release R: verified, manifest hash abc123"; exit 0 ;;\n'
            '  *"cli releases "*) echo "r26.11.4 listed"; exit 0 ;;\n'
            '  *train_plan.py*--json*) cat "$FAKE_PLAN"; exit 0 ;;\n'
            '  *train_plan.py*) echo "Action: cut"; exit 0 ;;\n'
            '  *train_est.py*) echo "{}"; exit 0 ;;\n'
            '  *" pack verify "*) echo "pack verified"; exit 0 ;;\n'
            'esac\necho "RESULT: fake"\nexit 0\n'
        ),
        "python3": (
            '#!/bin/sh\ncase "${2:-}" in -c|-) exec "$REAL_PYTHON" "$@" ;; esac\n'
            'case "$*" in\n  *train_check.py*)\n'
            '    out=""; prev=""\n    for a in "$@"; do [ "$prev" = "--out" ] && out=$a; prev=$a; done\n'
            '    mkdir -p "$out"\n'
            '    printf "path\\tbytes\\tsha256\\nx\\t1\\tab\\n" > "$out/root_added.tsv"\n'
            '    echo "# announcement" > "$out/announcement.md"\n    echo "check ok" ;;\n'
            '  *) echo "{}" ;;\nesac\nexit 0\n'
        ),
    }
    for name, body in stubs.items():
        p = fake / name
        p.write_text(body)
        p.chmod(0o755)

    (tree / "scripts" / "pod_run.sh").write_text(
        '#!/bin/sh\nlog=""; unit=""\n'
        'while [ $# -gt 0 ]; do\n  case "$1" in\n'
        '    --log) log=$2; shift 2 ;;\n    --*) shift 2 ;;\n'
        '    *) unit=$1; shift; break ;;\n  esac\ndone\n'
        'echo "pod_run unit=$unit cmd=$*" >> "$CAPTURE_LOG"\n'
        'rel=""; prev=""\nfor a in "$@"; do [ "$prev" = "--release" ] && rel=$a; prev=$a; done\n'
        'mkdir -p "$(dirname "$log")"\n'
        '{\n'
        '  [ "$1" = "build-stage" ] && echo "build-stage: stage P1 and release $rel: all 9 steps done"\n'
        '  [ "$1" = "pack" ] && echo "pack 02dec8312fff239cfed68b255e02f60e built"\n'
        '  echo "pod_run: $unit exited 0"\n'
        '} >> "$log"\nexit 0\n'
    )
    (tree / "scripts" / "release_examples.sh").write_text("#!/bin/sh\nexit 0\n")

    drv = bin_dir / "release_train.sh"
    drv.write_text(RELEASE_TRAIN.read_text().replace("/tmp/if2453", str(pod)))
    for other in POD_DIR.iterdir():
        link = bin_dir / other.name
        if other.name != "release_train.sh" and not link.exists():
            link.symlink_to(other)
    capture.write_text("")

    env = {
        **{k: v for k, v in os.environ.items() if k not in ("CLUSTER_METHOD", "BIOBOX_TRAIN_NOW")},
        "PATH": f"{fake}:{os.environ['PATH']}",
        "CAPTURE_LOG": str(capture),
        "FAKE_PLAN": str(plan_file),
        "REAL_PYTHON": sys.executable,
        "TREE": str(tree),
        "ROOT": str(root),
        "WORK": str(work),
        "REAL_ROOT": str(tmp_path / "other_root"),
        "START_MIN_GB": "0",
        "TMP_MIN_GB": "0",
        "PACK_PARTS": "2",
        "PACK_JOBS": "2",
        "PACK_MEM": "8G",
        "EX_THREADS": "2",
        "EX_MEM": "8G",
    }
    return drv, capture, root, work, env


def test_release_train_wires_cluster_method_default_and_override(tmp_path):
    # Default: hybrid-b1 passed to train_plan.py, train_est.py before, and build-stage
    drv, capture, _, _, env = _make_stub_train_world(tmp_path / "d1", plan_cluster_method="hybrid-b1")
    got = subprocess.run(["bash", str(drv), "run", "--no-live-cpu"], env=env, capture_output=True, text=True, timeout=60)
    assert got.returncode == 0, got.stdout + got.stderr
    lines = capture.read_text().splitlines()
    assert any("train_plan.py" in ln and "--cluster-method hybrid-b1" in ln for ln in lines), lines
    assert any("train_est.py before" in ln and "--cluster-method hybrid-b1" in ln for ln in lines), lines
    assert any("cmd=build-stage " in ln and "--cluster-method hybrid-b1" in ln for ln in lines), lines

    # Explicit override: CLUSTER_METHOD=cascade
    drv2, capture2, _, _, env2 = _make_stub_train_world(tmp_path / "d2", plan_cluster_method="cascade")
    env2["CLUSTER_METHOD"] = "cascade"
    got2 = subprocess.run(["bash", str(drv2), "run", "--no-live-cpu"], env=env2, capture_output=True, text=True, timeout=60)
    assert got2.returncode == 0, got2.stdout + got2.stderr
    lines2 = capture2.read_text().splitlines()
    assert any("train_plan.py" in ln and "--cluster-method cascade" in ln for ln in lines2), lines2
    assert any("train_est.py before" in ln and "--cluster-method cascade" in ln for ln in lines2), lines2
    assert any("cmd=build-stage " in ln and "--cluster-method cascade" in ln for ln in lines2), lines2


def test_release_train_backwards_compatible_resume_without_cluster_method(tmp_path):
    """Resuming a pre-A-27 $WORK whose plan.json has no cluster_method and whose driver.log
    build-stage command has no --cluster-method must succeed without command-mismatch refusal."""
    drv, capture, root, work, env = _make_stub_train_world(tmp_path, plan_cluster_method=None)

    # Pre-populate $WORK as a pre-A-27 unfinished launch of r26.11.4
    (work / "before").mkdir(parents=True, exist_ok=True)
    (work / "estimate").mkdir(parents=True, exist_ok=True)
    (work / "logs").mkdir(parents=True, exist_ok=True)
    (work / "run" / "c46-train-rel-1").mkdir(parents=True, exist_ok=True)
    (work / "plan.json").write_text(Path(env["FAKE_PLAN"]).read_text())
    (work / "before" / "listing.tsv").write_text("manifests/release=r26.11.3.json\t3\t1.0\n")
    (work / "before" / "verify-r26.11.3.log").write_text("release r26.11.3: verified, manifest hash abc123\n")
    (work / "estimate" / "before.json").write_text("{}\n")
    prof_src = (Path(env["TREE"]) / "profiles" / "wilkcomp02.yaml").read_text()
    (work / "profile.yaml").write_text(
        "\n".join(f"tmp: {work}/tmp" if ln.startswith("tmp:") else ln for ln in prof_src.splitlines()) + "\n"
    )
    # Mark r26.11.4 as started (one run record) and not yet finished (no manifest)
    (root / "runs" / "stage=r26.11.4").mkdir(parents=True, exist_ok=True)
    (root / "runs" / "stage=r26.11.4" / "cluster.plan.json").write_text('{"step":"cluster.plan"}\n')
    first_cmd = (
        f"build-stage --profile {work}/profile.yaml --root {root} "
        f"--stage P1 --release r26.11.4 --stages S1,P1 --depth called --genomes g.json --objects o.json "
        f"--threads 18 --scan-jobs 2 --scan-cpus 18"
    )
    (work / "driver.log").write_text(f"2026-11-24 22:00:00 CST launch c46-train-rel-1: {first_cmd}\n")

    # Dry-run must accept the resume (exit 0) without adding --cluster-method
    dry = subprocess.run(["bash", str(drv), "--dry-run", "--no-live-cpu"], env=env, capture_output=True, text=True, timeout=60)
    assert dry.returncode == 0, dry.stdout + dry.stderr
    assert "Action: resume r26.11.4" in dry.stdout, dry.stdout
    assert "the same as the first launch's" in dry.stdout, dry.stdout

    # Run must resume cleanly and launch the exact same command (no --cluster-method)
    got = subprocess.run(["bash", str(drv), "run", "--no-live-cpu"], env=env, capture_output=True, text=True, timeout=60)
    assert got.returncode == 0, got.stdout + got.stderr
    bs_lines = [ln for ln in capture.read_text().splitlines() if "cmd=build-stage " in ln]
    assert len(bs_lines) == 1 and "--cluster-method" not in bs_lines[0], bs_lines


def test_train_plan_accepts_cluster_method_and_detects_b49_reuse(tmp_path, monkeypatch):
    """train_plan.py --cluster-method hybrid-b1 detects B-49 reuse of a hybrid-b1 release
    and does not cross-reuse a cascade release."""
    if shutil.which("mmseqs") is None:
        monkeypatch.setattr(provenance, "tool_versions", lambda *names: {n: "stub" for n in names})

    tp = _load("train_plan_a27", TRAIN_PLAN)
    root = tmp_path / "root"
    _write_stage_stubs(root, "S1")
    from biobox import profile as bprofile

    prof = bprofile.load(PROFILE)
    cut = l2.cut_of(2, None, float(prof.cluster_cpu_min_max))
    inputs_b1 = l2.plan_inputs(root, "r26.11.2", GCV, ["S1"], cut, mode="hybrid-b1")
    pdir = root / "l2" / "cluster_plan" / "release=r26.11.2"
    pdir.mkdir(parents=True, exist_ok=True)
    (pdir / "cluster.plan.json").write_text(json.dumps({**inputs_b1, "proteins": 1, "outputs": []}))
    # Write a minimal cluster.json so _source_problem sees r26.11.2 as complete for cluster tables check
    # (or check reuse_reasons directly for method mismatch vs no method mismatch)
    est_b1 = tp.estimate_and_reuse(root, "r26.11.3", ["S1"], prof, parts=2, cluster_method="hybrid-b1")
    est_cas = tp.estimate_and_reuse(root, "r26.11.3", ["S1"], prof, parts=2, cluster_method="cascade")
    assert est_b1["cluster_method"] == "hybrid-b1"
    assert est_cas["cluster_method"] == "cascade"
    assert "method" not in (est_b1["reuse_reasons"].get("r26.11.2") or "")
    assert "method" in (est_cas["reuse_reasons"].get("r26.11.2") or "")
