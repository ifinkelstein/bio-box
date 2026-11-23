"""Ring R1: card A-22, r50_watch.py honours DEC-29's one_process clustering exemption (DEC-38)."""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("A-22", "DEC-29", "DEC-38")  # ring R1 by its directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
WATCH = REPO / "scripts" / "stress" / "build_stage_pod" / "r50_watch.py"

REL = "r26.11.3"
CALL = "cluster.nr90.11.prefilter"
OTHER_CALL = "cluster.nr90.12.align"
DEC29_REASON = "exempt: DEC-29"


def _load_watch():
    spec = importlib.util.spec_from_file_location("r50_watch_a22", WATCH)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _write_cut(root: Path, rel: str = REL, call: str = CALL, doc: dict | str | None = None) -> Path:
    path = root / "l2" / "cluster_plan" / f"release={rel}" / f"{call}.cut.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    if doc is None:
        doc = {
            "parts": [{"first": 0, "last": 9, "keys": 10, "units": 100, "cpu_min": 4200.0}],
            "one_process": {
                "decision": "DEC-29",
                "cpu_min": 4200.0,
                "cpu_min_max": 20000.0,
                "why": "per-part k-mer table rebuild exceeds part budget",
            },
        }
    if isinstance(doc, str):
        path.write_text(doc)
    else:
        path.write_text(json.dumps(doc))
    return path


def _write_run_record(root: Path, rel: str, name: str, doc: dict, *, mtime: float) -> Path:
    path = root / "runs" / f"stage={rel}" / f"{name}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc))
    os.utime(path, (mtime, mtime))
    return path


def _setup_env(tmp_path: Path) -> tuple[Path, Path, Path, Path, float]:
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


def _run_watch(
    monkeypatch,
    capsys,
    tmp_path: Path,
    *,
    live_rows: list[dict] | None = None,
    cpu_max: float = 280.0,
    cluster_cpu_max: float | None = None,
    stages: str | None = None,
) -> tuple[int, dict]:
    mod = _load_watch()
    root = tmp_path / "root"
    scratch = tmp_path / "scratch"
    since = tmp_path / "since"
    state = tmp_path / "state"
    argv = [
        str(WATCH),
        "--root", str(root),
        "--scratch", str(scratch),
        "--since", str(since),
        "--state", str(state),
        "--tmp-min-gb", "0",
        "--scratch-max-gb", "1000",
        "--cpu-max", str(cpu_max),
    ]
    if cluster_cpu_max is not None:
        argv.extend(["--cluster-cpu-max", str(cluster_cpu_max)])
    if stages is not None:
        argv.extend(["--stages", stages])
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(mod, "live", lambda _prefix: list(live_rows or []))
    rc = mod.main()
    captured = capsys.readouterr().out.strip()
    return rc, json.loads(captured)


# --- DEC-29 live process exemption -----------------------------------------------------------------------------


@pytest.mark.parametrize("action", ["step", "redo"])
def test_live_dec29_one_process_passes_under_cap_and_fails_above_cap(monkeypatch, capsys, tmp_path, action):
    root, _, _, _, _ = _setup_env(tmp_path)
    _write_cut(root)
    step = f"release cluster {action} release={REL} name={CALL}"

    rc, out = _run_watch(monkeypatch, capsys, tmp_path, live_rows=[{"pid": 1234, "step": step, "cpu_min": 3500.0}])
    assert rc == 0, out
    assert out["reasons"] == []

    rc, out = _run_watch(monkeypatch, capsys, tmp_path, live_rows=[{"pid": 1234, "step": step, "cpu_min": 20000.5}])
    assert rc == 3, out
    assert any("20000" in r for r in out["reasons"]), out["reasons"]


def test_step_of_formats_release_cluster_step_and_redo():
    mod = _load_watch()
    for action in ("step", "redo"):
        cmd = [
            "python", "-m", "biobox.cli", "release", "cluster", action,
            "--root", "/tmp/root", "--release", REL, "--name", CALL,
        ]
        assert mod.step_of(cmd) == f"release cluster {action} release={REL} name={CALL}"


# --- DEC-29 published run record exemption ---------------------------------------------------------------------


def test_published_dec29_one_process_passes_under_cap_and_fails_above_cap(monkeypatch, capsys, tmp_path):
    root, scratch, since, state, t0 = _setup_env(tmp_path)
    _write_cut(root)
    _write_run_record(root, REL, CALL, {"cpu_min": 3500.0, "cpu_budget_reason": DEC29_REASON}, mtime=t0)

    rc, out = _run_watch(monkeypatch, capsys, tmp_path)
    assert rc == 0, out
    assert out["reasons"] == []

    # Also verify isolated stdlib-only invocation (`python3 -I`).
    proc = subprocess.run(
        [
            sys.executable, "-I", str(WATCH),
            "--root", str(root),
            "--scratch", str(scratch),
            "--since", str(since),
            "--state", str(state),
            "--tmp-min-gb", "0",
            "--scratch-max-gb", "1000",
            "--cpu-max", "280",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr + proc.stdout
    assert json.loads(proc.stdout)["reasons"] == []

    _write_run_record(root, REL, CALL, {"cpu_min": 20000.5, "cpu_budget_reason": DEC29_REASON}, mtime=t0)
    rc, out = _run_watch(monkeypatch, capsys, tmp_path)
    assert rc == 3, out
    assert any("20000" in r for r in out["reasons"]), out["reasons"]


# --- Strictness / security checks ------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "record_doc",
    [
        {"cpu_min": 3500.0},
        {"cpu_min": 3500.0, "cpu_budget_reason": "exempt: DEC-25"},
        {"cpu_min": 3500.0, "cpu_budget_reason": ""},
    ],
    ids=["missing-reason", "wrong-reason", "empty-reason"],
)
def test_published_record_without_dec29_reason_still_fails(monkeypatch, capsys, tmp_path, record_doc):
    root, _, _, _, t0 = _setup_env(tmp_path)
    _write_cut(root)
    _write_run_record(root, REL, CALL, record_doc, mtime=t0)

    rc, out = _run_watch(monkeypatch, capsys, tmp_path)
    assert rc == 3, out
    assert any("280" in r for r in out["reasons"]), out["reasons"]


@pytest.mark.parametrize(
    "cut_doc",
    [
        None,  # missing cut.json
        "{not json",
        {
            "parts": [{"first": 0, "last": 9}],
            "one_process": {"decision": "DEC-25", "cpu_min": 4200.0, "cpu_min_max": 20000.0},
        },
        {
            "parts": [{"first": 0, "last": 4}, {"first": 5, "last": 9}],
            "one_process": {"decision": "DEC-29", "cpu_min": 4200.0, "cpu_min_max": 20000.0},
        },
        {
            "parts": [],
            "one_process": {"decision": "DEC-29", "cpu_min": 4200.0, "cpu_min_max": 20000.0},
        },
        {
            "parts": [{"first": 0, "last": 9}],
            "one_process": {"decision": "DEC-29", "cpu_min": 4200.0, "cpu_min_max": 0.0},
        },
        {
            "parts": [{"first": 0, "last": 9}],
            "one_process": {"decision": "DEC-29", "cpu_min": 4200.0, "cpu_min_max": -100.0},
        },
        {
            "parts": [{"first": 0, "last": 9}],
            "one_process": {"decision": "DEC-29", "cpu_min": 4200.0, "cpu_min_max": True},
        },
        {
            "parts": [{"first": 0, "last": 9}],
            "one_process": {"decision": "DEC-29", "cpu_min": 4200.0, "cpu_min_max": "20000"},
        },
        '{"parts": [{"first": 0, "last": 9}], "one_process": {"decision": "DEC-29", "cpu_min": 4200.0, "cpu_min_max": Infinity}}',
        '{"parts": [{"first": 0, "last": 9}], "one_process": {"decision": "DEC-29", "cpu_min": 4200.0, "cpu_min_max": NaN}}',
    ],
    ids=[
        "missing-cut",
        "invalid-json",
        "wrong-decision",
        "two-parts",
        "zero-parts",
        "zero-cap",
        "negative-cap",
        "bool-cap",
        "str-cap",
        "inf-cap",
        "nan-cap",
    ],
)
def test_invalid_or_missing_cut_json_rejects_dec29_for_published_and_live(monkeypatch, capsys, tmp_path, cut_doc):
    root, _, _, _, t0 = _setup_env(tmp_path)
    if cut_doc is not None:
        _write_cut(root, doc=cut_doc)

    # Published record with DEC-29 reason must still fail at 3500 > 280.
    _write_run_record(root, REL, CALL, {"cpu_min": 3500.0, "cpu_budget_reason": DEC29_REASON}, mtime=t0)
    rc, out = _run_watch(monkeypatch, capsys, tmp_path)
    assert rc == 3, out
    assert any("280" in r for r in out["reasons"]), out["reasons"]

    # Live process must also fail at 3500 > 280.
    (root / "runs" / f"stage={REL}" / f"{CALL}.json").unlink()
    step = f"release cluster step release={REL} name={CALL}"
    rc, out = _run_watch(monkeypatch, capsys, tmp_path, live_rows=[{"pid": 1234, "step": step, "cpu_min": 3500.0}])
    assert rc == 3, out
    assert any("280" in r for r in out["reasons"]), out["reasons"]


@pytest.mark.parametrize(
    ("record_name", "live_step"),
    [
        (f"{CALL}.part-00000", f"release cluster step release={REL} name={CALL}.part-00000"),
        (f"{CALL}.cut", f"release cluster step release={REL} name={CALL}.cut"),
        ("cluster", f"release cluster finish release={REL}"),
        (OTHER_CALL, f"release cluster step release={REL} name={OTHER_CALL}"),
    ],
    ids=["split-part", "cut-step", "finish", "other-call"],
)
def test_non_exempt_clustering_steps_still_fail_at_280(monkeypatch, capsys, tmp_path, record_name, live_step):
    root, _, _, _, t0 = _setup_env(tmp_path)
    _write_cut(root)  # valid DEC-29 cut exists for CALL only

    # Published record fails even if it falsely claims DEC-29.
    _write_run_record(root, REL, record_name, {"cpu_min": 3500.0, "cpu_budget_reason": DEC29_REASON}, mtime=t0)
    rc, out = _run_watch(monkeypatch, capsys, tmp_path)
    assert rc == 3, out
    assert any("280" in r for r in out["reasons"]), out["reasons"]

    # Live process also fails at 3500 > 280.
    (root / "runs" / f"stage={REL}" / f"{record_name}.json").unlink()
    rc, out = _run_watch(monkeypatch, capsys, tmp_path, live_rows=[{"pid": 2000, "step": live_step, "cpu_min": 3500.0}])
    assert rc == 3, out
    assert any("280" in r for r in out["reasons"]), out["reasons"]


def test_dec29_cap_rejects_traversal_and_cap_below_default(tmp_path):
    mod = _load_watch()
    root = tmp_path / "root"
    _write_cut(
        root,
        doc={
            "parts": [{"first": 0, "last": 9}],
            "one_process": {"decision": "DEC-29", "cpu_min": 150.0, "cpu_min_max": 200.0},
        },
    )
    assert mod.dec29_cap(root, REL, CALL, 280.0) == 280.0
    assert mod.dec29_cap(root, f"../{REL}", CALL, 280.0) == 280.0
    assert mod.dec29_cap(root, REL, f"cluster.foo/{CALL}", 280.0) == 280.0

