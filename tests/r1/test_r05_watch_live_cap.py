"""Ring R1: card R-05, the watch raises a running cut's CPU cap only where the plan grants DEC-77's exemption.

`build-stage` (`src/biobox/stage.py`, `_cluster_in_parts`) runs a shared-index prefilter's `<call>.cut` step with
`budget_reason="exempt: DEC-77"` only when `l2.cluster.cut_createindex_estimate` (the shared `createindex_cpu_min`
over the call's target in the plan's own work directory) is over the plan's part budget (`cut.part_cpu_min`), and
refuses to launch it when that estimate is over the profile's `cluster_cpu_min_max`. `r50_watch.py` must apply the
same test to a live `.cut` process whose cut record is not yet written: the exempt cap (the plan's
`cut.one_process_cpu_min_max`) where stage.py grants it, `--cpu-max` everywhere else (no plan, an unreadable plan, a
plan without its input fields, a step of another clustering method, no target index, an unsafe name).

The expected decision is the library's own (`l2.cut_createindex_estimate` on the same fixture), so a drift of the
watch's stdlib mirror from `l2/cluster.py` fails here. No MMseqs2 is needed: the target is an MMseqs2 `.index` file
written by the test.
"""
from __future__ import annotations

import importlib.util
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

from biobox.l2 import cluster as l2

pytestmark = pytest.mark.covers("R-05", "A-27", "DEC-77", "DEC-29", "DEC-38")  # ring R1 by its directory

REPO = Path(__file__).resolve().parents[2]
WATCH = REPO / "scripts" / "stress" / "build_stage_pod" / "r50_watch.py"

REL = "r26.11.4"
CUT_15 = "cluster.nr50.15.prefilter.cut"
CUT_19 = "cluster.nr50.19.prefilter.cut"
BUDGET = 240.0  # SHARD_MARGIN (0.8) of the 300 CPU-minute budget, as build-stage gives it
CAP = 20000.0  # profiles/wilkcomp02.yaml cluster_cpu_min_max, as build-stage gives it to the plan
TARGETS = 100


def _load_watch():
    spec = importlib.util.spec_from_file_location("r50_watch_r05", WATCH)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _call(step: str, mode: str) -> l2.Call:
    return {c.name: c for c in l2.calls(mode)}[step[: -len(".cut")]]


def _plan_doc(rel: str = REL, *, method: str | None = "hybrid-b1", fmt: int = 2, cut: dict | None = None) -> dict:
    doc = {
        "format": fmt,
        "release": rel,
        "stages": ["S1"],
        "gene_call_version": "gc-r05",
        "mmseqs2": "18.8cc5c",
        "files": {"l1/protein_seq/stage=S1/part-00000.parquet": "0" * 64},
        "cut": {"part_cpu_min": BUDGET, "one_process_cpu_min_max": CAP} if cut is None else cut,
        "proteins": 10_000,
        "release_cpu_min": 1.0,
        "outputs": {},
    }
    if method is not None:
        doc["method"] = method
    return doc


def _write_plan(root: Path, doc: dict | str, rel: str = REL) -> Path:
    path = root / "l2" / "cluster_plan" / f"release={rel}" / "cluster.plan.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(doc if isinstance(doc, str) else json.dumps(doc, sort_keys=True, indent=1) + "\n")
    return path


def _per_row(est: float, sens: str) -> int:
    """The entry length (residues + newline + NUL) of each of TARGETS entries whose createindex estimate is `est`."""
    rate = l2.CREATEINDEX_RATE_BY_SENS[sens]
    residues = est / (rate * 1e-9 * (TARGETS / l2.CREATEINDEX_T0) ** l2.CREATEINDEX_GROWTH)
    return math.ceil(residues / TARGETS) + 2


def _write_index(tmp: Path, root: Path, plan: dict, step: str, est: float) -> Path:
    """The call's target `.index` in the plan's work directory under `tmp`, sized so that its createindex estimate is
    `est` CPU-minutes (just above), as an MMseqs2 index: key, offset, length; keys sorted."""
    mode = plan.get("method", l2.DEFAULT_CLUSTERING_MODE)
    c = _call(step, mode)
    sens = l2._opt(c.argv, "-s")
    w = l2._work_dir(tmp, root, l2._plan_input_dict(plan))
    index = Path(f"{w / c.positional[1]}.index")
    index.parent.mkdir(parents=True, exist_ok=True)
    n = _per_row(est, sens)
    index.write_text("".join(f"{k}\t{k * n}\t{n}\n" for k in range(TARGETS)))
    return index


def _setup(tmp_path: Path) -> tuple[Path, Path]:
    root, scratch, state = tmp_path / "root", tmp_path / "scratch", tmp_path / "state"
    for d in (root, scratch, state):
        d.mkdir()
    since = tmp_path / "since"
    since.write_text("mark\n")
    t0 = since.stat().st_mtime
    os.utime(since, (t0 - 60, t0 - 60))
    return root, scratch


def _run_watch(monkeypatch, capsys, tmp_path: Path, live_rows: list[dict], cpu_max: float = 280.0):
    mod = _load_watch()
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
    monkeypatch.setattr(mod, "live", lambda _prefix: [dict(r) for r in live_rows])
    rc = mod.main()
    return rc, json.loads(capsys.readouterr().out.strip())


def _live(step: str, cpu_min: float, action: str = "step", **extra) -> dict:
    return {"pid": 4242, "step": f"release cluster {action} release={REL} name={step}", "cpu_min": cpu_min, **extra}


def _library_grants(root: Path, tmp: Path, step: str, mode: str) -> float | None:
    """stage.py's test: the estimate for which it passes `exempt: DEC-77`, or None (`cut_createindex_estimate`)."""
    return l2.cut_createindex_estimate(root, REL, _call(step, mode), tmp)


CASES = [(CUT_15, "hybrid-b1"), (CUT_15, "cascade"), (CUT_19, "cascade"), (CUT_19, "hybrid-a")]
CASE_IDS = [f"{s.split('.')[2]}-{m}" for s, m in CASES]


# --- (1) under the part budget: the default cap ----------------------------------------------------------------


@pytest.mark.parametrize("action", ["step", "redo"])
@pytest.mark.parametrize(("step", "mode"), CASES, ids=CASE_IDS)
@pytest.mark.parametrize("est", [100.0, 239.0])
def test_live_cut_under_part_budget_is_capped_at_default(monkeypatch, capsys, tmp_path, step, mode, est, action):
    root, scratch = _setup(tmp_path)
    plan = _plan_doc(method=mode)
    _write_plan(root, plan)
    _write_index(scratch, root, plan, step, est)
    assert _library_grants(root, scratch, step, mode) is None  # stage.py runs it without exempt: DEC-77

    rc, out = _run_watch(monkeypatch, capsys, tmp_path, [_live(step, 300.0, action)])
    assert rc == 3, out
    assert any("over 280" in r for r in out["reasons"]), out["reasons"]

    rc, out = _run_watch(monkeypatch, capsys, tmp_path, [_live(step, 279.0, action)])
    assert rc == 0 and out["reasons"] == [], out


# --- (2) over the part budget and within cluster_cpu_min_max: the exempt cap ------------------------------------


@pytest.mark.parametrize(("step", "mode"), CASES, ids=CASE_IDS)
@pytest.mark.parametrize("est", [241.0, 600.0, 19_000.0])
def test_live_cut_over_part_budget_within_cap_gets_exempt_cap(monkeypatch, capsys, tmp_path, step, mode, est):
    root, scratch = _setup(tmp_path)
    plan = _plan_doc(method=mode)
    _write_plan(root, plan)
    _write_index(scratch, root, plan, step, est)
    got = _library_grants(root, scratch, step, mode)
    assert got is not None and BUDGET < got <= CAP  # stage.py grants exempt: DEC-77

    rc, out = _run_watch(monkeypatch, capsys, tmp_path, [_live(step, 3_500.0)])
    assert rc == 0 and out["reasons"] == [], out

    rc, out = _run_watch(monkeypatch, capsys, tmp_path, [_live(step, CAP + 0.5)])
    assert rc == 3, out
    assert any("over 20000" in r for r in out["reasons"]), out["reasons"]


# --- an estimate over cluster_cpu_min_max: stage.py refuses to launch it; the watch gives it no exemption ---------


@pytest.mark.parametrize(("step", "mode"), CASES, ids=CASE_IDS)
def test_live_cut_estimated_over_cluster_cap_gets_default(monkeypatch, capsys, tmp_path, step, mode):
    root, scratch = _setup(tmp_path)
    plan = _plan_doc(method=mode)
    _write_plan(root, plan)
    _write_index(scratch, root, plan, step, 30_000.0)
    got = _library_grants(root, scratch, step, mode)
    assert got is not None and got > CAP  # stage.py: one_process_refusal, never launched

    rc, out = _run_watch(monkeypatch, capsys, tmp_path, [_live(step, 300.0)])
    assert rc == 3 and any("over 280" in r for r in out["reasons"]), out


# --- the watch's decision is stage.py's, across the budget -------------------------------------------------------


@pytest.mark.parametrize(("step", "mode"), CASES, ids=CASE_IDS)
def test_live_cut_decision_equals_stage_test(tmp_path, step, mode):
    mod = _load_watch()
    root, scratch = _setup(tmp_path)
    for i, est in enumerate((1.0, 120.0, 239.9, 240.1, 241.0, 2_000.0, 19_999.0, 20_001.0, 1e6)):
        plan = _plan_doc(method=mode, cut={"part_cpu_min": BUDGET, "one_process_cpu_min_max": CAP,
                                          "rebuild_cpu_min": float(i)})  # a new work directory each time
        _write_plan(root, plan)
        _write_index(scratch, root, plan, step, est)
        lib = _library_grants(root, scratch, step, mode)
        want = CAP if lib is not None and lib <= CAP else 280.0
        got = mod.dec77_cut_cap(root, REL, step, 280.0, live_proc=True, tmp=scratch)
        assert got == want, (est, lib, got)


# --- no plan, an unreadable plan, an unsafe name, or nothing to estimate: the default ----------------------------


def _bad_plans() -> list[tuple[str, dict | str]]:
    no_inputs = _plan_doc(method="cascade")
    del no_inputs["files"]
    return [
        ("invalid-json", "{not json"),
        ("not-an-object", "[1, 2]"),
        ("no-input-fields", no_inputs),
        ("other-method-step", _plan_doc(method="hybrid-b1")),  # CUT_19 is no step of hybrid-b1
        ("linclust", _plan_doc(method="linclust")),
        ("unknown-method", _plan_doc(method="hybrid-z")),
        ("forced-parts", _plan_doc(method="cascade", cut={"parts": 4, "one_process_cpu_min_max": CAP})),
        ("no-cap", _plan_doc(method="cascade", cut={"part_cpu_min": BUDGET})),
        ("bool-budget", _plan_doc(method="cascade", cut={"part_cpu_min": True, "one_process_cpu_min_max": CAP})),
        ("str-cap", _plan_doc(method="cascade", cut={"part_cpu_min": BUDGET, "one_process_cpu_min_max": "20000"})),
        ("unsafe-plan-release", _plan_doc(rel="../r26.11.4", method="cascade")),
    ]


@pytest.mark.parametrize("which", [name for name, _ in _bad_plans()])
def test_live_cut_with_unusable_plan_gets_default(monkeypatch, capsys, tmp_path, which):
    root, scratch = _setup(tmp_path)
    doc = dict(_bad_plans())[which]
    _write_plan(root, doc)
    # a target sized over the part budget, in the work directory a cascade plan of this release would have
    _write_index(scratch, root, _plan_doc(method="cascade"), CUT_19, 600.0)
    if isinstance(doc, dict) and "files" in doc and doc.get("method") in l2.CLUSTERING_MODES:
        try:
            _write_index(scratch, root, doc, CUT_19, 600.0)
        except (KeyError, ValueError):
            pass  # no such step under this method
    rc, out = _run_watch(monkeypatch, capsys, tmp_path, [_live(CUT_19, 300.0)])
    assert rc == 3 and any("over 280" in r for r in out["reasons"]), out


def test_live_cut_without_plan_or_target_index_gets_default(monkeypatch, capsys, tmp_path):
    root, scratch = _setup(tmp_path)
    rc, out = _run_watch(monkeypatch, capsys, tmp_path, [_live(CUT_15, 300.0)])  # no plan
    assert rc == 3 and any("over 280" in r for r in out["reasons"]), out

    plan = _plan_doc()
    _write_plan(root, plan)
    rc, out = _run_watch(monkeypatch, capsys, tmp_path, [_live(CUT_15, 300.0)])  # no target index in scratch
    assert rc == 3 and any("over 280" in r for r in out["reasons"]), out

    index = _write_index(scratch, root, plan, CUT_15, 600.0)
    index.write_text("0\t0\tnot-a-number\n")  # unreadable index
    rc, out = _run_watch(monkeypatch, capsys, tmp_path, [_live(CUT_15, 300.0)])
    assert rc == 3 and any("over 280" in r for r in out["reasons"]), out


def test_unsafe_names_get_default(tmp_path):
    mod = _load_watch()
    root, scratch = _setup(tmp_path)
    plan = _plan_doc()
    _write_plan(root, plan)
    _write_index(scratch, root, plan, CUT_15, 600.0)
    assert mod.dec77_cut_cap(root, REL, CUT_15, 280.0, live_proc=True, tmp=scratch) == CAP  # the safe control
    for rel, step in ((f"../{REL}", CUT_15), (f"{REL}/..", CUT_15), (REL, f"x/{CUT_15}"), (REL, f"..{CUT_15}"),
                      (REL, "cluster.nr90.11.prefilter.cut"), (REL, "cluster.nr50.16.align.cut"), ("", CUT_15)):
        assert mod.dec77_cut_cap(root, rel, step, 280.0, live_proc=True, tmp=scratch) == 280.0, (rel, step)
    assert mod.dec77_cut_cap(root, REL, CUT_15, 280.0, live_proc=True, tmp=None) == 280.0  # no scratch known
    assert mod.dec77_cut_cap(root, REL, CUT_15, 280.0, live_proc=False, tmp=scratch) == 280.0  # not live


# --- a format-1 plan (no method: cascade) hashes its work directory without `method`, and grants no DEC-77 ------
# Card R-02: a PLAN_FORMAT 1 cut keeps the pre-B-71 record and refusal, so `cut_createindex_estimate` is None for it
# and build-stage never runs its cut with `exempt: DEC-77`; the watch gives it the default cap.


@pytest.mark.parametrize("step", [CUT_15, CUT_19])
@pytest.mark.parametrize("est", [100.0, 600.0, 19_000.0])
def test_format1_plan_work_dir_and_no_exemption(monkeypatch, capsys, tmp_path, step, est):
    mod = _load_watch()
    root, scratch = _setup(tmp_path)
    plan = _plan_doc(method=None, fmt=1)
    _write_plan(root, plan)
    assert mod.work_dir(scratch, root, plan) == l2._work_dir(scratch, root, l2._plan_input_dict(plan))
    plan2 = _plan_doc(method="cascade")
    assert mod.work_dir(scratch, root, plan2) == l2._work_dir(scratch, root, l2._plan_input_dict(plan2))
    assert mod.work_dir(scratch, root, plan) != mod.work_dir(scratch, root, plan2)

    index = _write_index(scratch, root, plan, step, est)  # in the format-1 plan's own work directory
    assert index.parent.parent.parent == mod.work_dir(scratch, root, plan)
    assert _library_grants(root, scratch, step, "cascade") is None  # card R-02: no DEC-77 for a format-1 plan
    assert mod.dec77_cut_cap(root, REL, step, 280.0, live_proc=True, tmp=scratch) == 280.0
    rc, out = _run_watch(monkeypatch, capsys, tmp_path, [_live(step, 300.0)])
    assert rc == 3 and any("over 280" in r for r in out["reasons"]), out

    # the same target under a format-2 cascade plan: exempt where the library grants it (format 2 unchanged)
    _write_plan(root, plan2)
    _write_index(scratch, root, plan2, step, est)
    lib = _library_grants(root, scratch, step, "cascade")
    want = CAP if lib is not None and lib <= CAP else 280.0
    assert (lib is not None) == (est > BUDGET)
    assert mod.dec77_cut_cap(root, REL, step, 280.0, live_proc=True, tmp=scratch) == want


# --- the process's own --root and --tmp locate its work directory ------------------------------------------------


def test_live_rows_carry_the_process_root_and_tmp(monkeypatch):
    mod = _load_watch()
    cg = "0::/user.slice/user-1.slice/user@1.service/app.slice/rt-rel-x.service\n"
    cmd = ["python", "-m", "biobox.cli", "release", "cluster", "step", "--root", "/r/root", "--gene-call-version",
           "g", "--release", REL, "--name", CUT_15, "--threads", "18", "--tmp", "/w/tmp/scratch"]
    procs = {
        10: {"ppid": 1, "ticks": 0, "cmd": ["bash", "pod_run.sh"], "cgroup": cg},
        11: {"ppid": 10, "ticks": 6000 * mod.TICK, "cmd": cmd, "cgroup": cg},
        12: {"ppid": 11, "ticks": 0, "cmd": ["mmseqs", "createindex"], "cgroup": cg},
    }
    monkeypatch.setattr(mod, "procs", lambda: procs)
    rows = mod.live("rt-rel-")
    assert rows == [{"pid": 11, "step": f"release cluster step release={REL} name={CUT_15}", "cpu_min": 100.0,
                     "root": "/r/root", "tmp": "/w/tmp/scratch"}]


def test_live_cut_uses_the_process_tmp_and_root_over_the_watch_args(monkeypatch, capsys, tmp_path):
    root, scratch = _setup(tmp_path)
    plan = _plan_doc()
    _write_plan(root, plan)
    own = tmp_path / "own-tmp"
    _write_index(own, root, plan, CUT_15, 600.0)  # the process's scratch, not the watch's --scratch
    rc, out = _run_watch(monkeypatch, capsys, tmp_path, [_live(CUT_15, 3_500.0, tmp=str(own), root=str(root))])
    assert rc == 0 and out["reasons"] == [], out
    assert all(set(r) == {"pid", "step", "cpu_min"} for r in out["live"]), out["live"]  # the line stays short

    rc, out = _run_watch(monkeypatch, capsys, tmp_path, [_live(CUT_15, 3_500.0)])  # the watch's --scratch: no index
    assert rc == 3 and any("over 280" in r for r in out["reasons"]), out

    _write_index(scratch, root, plan, CUT_15, 100.0)  # under budget where the process works; over in --scratch: own
    rc, out = _run_watch(monkeypatch, capsys, tmp_path, [_live(CUT_15, 300.0, tmp=str(scratch), root=str(root))])
    assert rc == 3 and any("over 280" in r for r in out["reasons"]), out


# --- the counts of a target index are kept in the state directory, and not reused once the index changes ---------


def test_index_counts_are_cached_and_invalidated(monkeypatch, capsys, tmp_path):
    mod = _load_watch()
    root, scratch = _setup(tmp_path)
    plan = _plan_doc()
    _write_plan(root, plan)
    index = _write_index(scratch, root, plan, CUT_15, 600.0)
    cache = tmp_path / "state" / mod.INDEX_COUNTS
    n = _per_row(600.0, "2")
    assert mod.index_counts(index, cache) == (TARGETS, TARGETS * (n - 2))
    assert cache.is_file()
    assert mod.index_counts(index, cache) == (TARGETS, TARGETS * (n - 2))

    index.write_text("".join(f"{k}\t0\t10\n" for k in range(TARGETS)))  # the target changed: count again
    os.utime(index, ns=(1, 1))
    assert mod.index_counts(index, cache) == (TARGETS, TARGETS * 8)
    rc, out = _run_watch(monkeypatch, capsys, tmp_path, [_live(CUT_15, 300.0)])
    assert rc == 3 and any("over 280" in r for r in out["reasons"]), out

    cache.write_text("{broken")  # a broken cache is counted again, never trusted
    index.write_text("".join(f"{k}\t{k * n}\t{n}\n" for k in range(TARGETS)))
    assert mod.index_counts(index, cache) == (TARGETS, TARGETS * (n - 2))
    rc, out = _run_watch(monkeypatch, capsys, tmp_path, [_live(CUT_15, 3_500.0)])
    assert rc == 0 and out["reasons"] == [], out


# --- the stdlib mirrors in the watch equal the library ----------------------------------------------------------


def test_watch_mirrors_equal_the_library():
    mod = _load_watch()
    assert mod.DEFAULT_CLUSTERING_MODE == l2.DEFAULT_CLUSTERING_MODE
    assert mod.WORK_DIR == l2.WORK_DIR
    assert mod.PLAN_DIR == l2.PLAN_DIR
    assert mod.PLAN_INPUT_KEYS_V1 == l2._INPUT_KEYS_V1 and mod.PLAN_INPUT_KEYS == l2._INPUT_KEYS
    assert set(mod.SHARED_INDEX_CUTS) == set(l2.CLUSTERING_MODES)
    for mode in l2.CLUSTERING_MODES:
        want = {f"{c.name}.cut": (c.positional[1], l2._opt(c.argv, "-s"))
                for c in l2.calls(mode) if c.split and c.shares_index}
        assert mod.SHARED_INDEX_CUTS[mode] == want, mode
    assert mod.SHARED_INDEX_CUT_STEPS == frozenset(s for m in mod.SHARED_INDEX_CUTS.values() for s in m)
    for res, tgt, sens in ((6_476_355_769, 24_298_242, "2"), (2_868_261_896, 10_235_527, "3"),
                           (18_500_000_000, 70_000_000, "2"), (1_000_000_000, 10_000_000, "4"),
                           (0, 10_000_000, "2"), (10, 0, "3"), (10**12, 100, "x")):
        assert mod.createindex_cpu_min(res, tgt, sens) == l2.createindex_cpu_min(res, tgt, sens)


def test_watch_runs_stdlib_only(tmp_path):
    root, scratch = _setup(tmp_path)
    plan = _plan_doc()
    _write_plan(root, plan)
    _write_index(scratch, root, plan, CUT_15, 600.0)
    proc = subprocess.run(
        [sys.executable, "-I", str(WATCH), "--root", str(root), "--scratch", str(scratch),
         "--since", str(tmp_path / "since"), "--state", str(tmp_path / "state"), "--tmp-min-gb", "0",
         "--scratch-max-gb", "1000", "--cpu-max", "280", "--unit-prefix", "r05-no-such-unit-"],
        capture_output=True, text=True, check=False,
    )
    assert proc.returncode == 0, proc.stderr + proc.stdout
    assert json.loads(proc.stdout)["reasons"] == []
