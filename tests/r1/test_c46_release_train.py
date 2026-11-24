"""Ring R1: unit tests for card C-46 (the release train: naming, stage selection, dry-run, rate guard, check)."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

import pytest

from biobox import catalog
from biobox.l1 import crispr
from biobox.l1.genecall import GENE_CALL_VERSION as GCV
from biobox.profile import Profile

pytestmark = pytest.mark.covers("C-46")  # ring R1 by its directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
POD_DIR = REPO / "scripts" / "stress" / "build_stage_pod"


@pytest.fixture(autouse=True)
def _mock_minced_if_absent(monkeypatch):
    if shutil.which("minced") is None:  # GitHub Actions R1 has no MinCED binary (card E-12)
        monkeypatch.setattr(crispr, "minced_version", lambda: "0.4.2")


def _import_pod_module(name: str):
    if str(POD_DIR) not in sys.path:
        sys.path.insert(0, str(POD_DIR))
    return __import__(name)


def _write_record(root: Path, stage: str, step: str, *, started_at: str = "2026-11-18T00:00:00+00:00",
                  exit_code: int = 0, outputs: list[dict] | None = None, argv: list[str] | None = None,
                  cpu_min: float = 1.0, wall_s: float = 10.0) -> Path:
    d = root / "runs" / f"stage={stage}"
    d.mkdir(parents=True, exist_ok=True)
    rec = {
        "stage": stage,
        "step": step,
        "exit_code": exit_code,
        "started_at": started_at,
        "cpu_min": cpu_min,
        "wall_s": wall_s,
        "peak_rss_bytes": 100_000_000,
        "outputs": outputs or [],
        "argv": argv or [],
    }
    p = d / f"{step}.json"
    p.write_text(json.dumps(rec, indent=1, sort_keys=True))
    return p


def _make_file(root: Path, rel: str, data: bytes = b"x") -> dict:
    p = root / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(data)
    return {"path": rel, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def _make_stage(root: Path, stage: str, *, depth: str = "called", started_at: str = "2026-11-18T00:00:00+00:00",
                libraries: list[tuple[str, str]] | None = None, crispr_enabled: bool = False,
                missing_step: str | None = None, bad_exit_step: str | None = None,
                intent_step: str | None = None, corrupt_output_step: str | None = None) -> None:
    """Create minimal valid stage run records and output files for `stage` at `depth`."""
    steps_to_write: list[tuple[str, list[str]]] = [
        ("inputs", ["python", "-m", "biobox.cli", "stage", "inputs", "--stage", stage,
                    "--genomes", f"/in/{stage}/genomes.json", "--gene-call-version", GCV]),
        ("ingest", ["python", "-m", "biobox.cli", "stage", "ingest", "--stage", stage,
                    "--objects", f"/in/{stage}/objects.json", "--gene-call-version", GCV]),
    ]
    if depth in ("called", "annotated"):
        steps_to_write.extend([
            ("genecall", ["python", "-m", "biobox.cli", "stage", "genecall", "--stage", stage]),
            ("igr", ["python", "-m", "biobox.cli", "stage", "igr", "--stage", stage]),
            ("tapes", ["python", "-m", "biobox.cli", "stage", "tapes", "--stage", stage]),
        ])
        if crispr_enabled:
            steps_to_write.append(
                (f"crispr-{crispr.producer()}",
                 ["python", "-m", "biobox.cli", "stage", "crispr", "--stage", stage])
            )
    if depth == "annotated":
        for lib_name, lib_ver in (libraries or [("pfam", "38.2")]):
            steps_to_write.append(
                (f"domains-{lib_name}-{lib_ver}",
                 ["python", "-m", "biobox.cli", "stage", "domains", "--stage", stage,
                  "--library", lib_name, "--library-version", lib_ver,
                  "--hmm", f"/ref/{lib_name}-{lib_ver}.hmm", "--cutoff", "gathering"])
            )
            if crispr_enabled and lib_name == crispr.CAS_LIBRARY:
                steps_to_write.append(
                    (f"crispr-{crispr.near_cas_producer(lib_ver)}",
                     ["python", "-m", "biobox.cli", "stage", "crispr", "--near-cas",
                      "--library", lib_name, "--library-version", lib_ver, "--stage", stage])
                )

    for step_name, argv in steps_to_write:
        if step_name == missing_step:
            continue
        out = _make_file(root, f"dummy/{stage}/{step_name}.dat", data=f"{stage}:{step_name}".encode())
        if step_name == corrupt_output_step:
            (root / out["path"]).write_bytes(b"corrupted-size")
        ec = 1 if step_name == bad_exit_step else 0
        _write_record(root, stage, step_name, started_at=started_at, exit_code=ec, outputs=[out], argv=argv)

    if intent_step is not None:
        (root / "runs" / f"stage={stage}" / f"{intent_step}.publishing.json").write_text('{"step": "x"}')


def _define_release(root: Path, release: str, stages: list[str], *, with_manifest: bool = True) -> None:
    catalog.define_release(root, release, stages)
    if with_manifest:
        m = root / "manifests" / f"release={release}.json"
        m.parent.mkdir(parents=True, exist_ok=True)
        m.write_text(json.dumps({"release": release, "stages": stages, "files": [], "artifacts": []}))


def _profile() -> Profile:
    return Profile(
        max_threads_day=40, max_threads_night=56, mem_scope="200G",
        tmp="/tmp/if2453/bio-box", durable="/stor/work/test", bulk="/stor/scratch/test",
        nice_monolithic=19, nice_fanout=15, day_start_hour=8, day_end_hour=20,
        cpu_min_per_process=300, small_file_floor_bytes=1000,
    )


# --- Rule (1): release naming (DEC-52) ----------------------------------------------------------------------


def test_next_release_name_fresh_and_incrementing(tmp_path):
    tp = _import_pod_module("train_plan")
    assert tp.next_release_name(tmp_path, date(2026, 11, 19)) == "r26.11.1"

    _define_release(tmp_path, "rN.3", ["S0"])
    _define_release(tmp_path, "rM.2", ["S0"])
    _define_release(tmp_path, "rS0.1", ["S0"])
    assert tp.next_release_name(tmp_path, date(2026, 11, 19)) == "r26.11.1"

    _define_release(tmp_path, "r26.11.1", ["S0"])
    assert tp.next_release_name(tmp_path, date(2026, 11, 19)) == "r26.11.2"

    # Numeric increment beyond 9 (not lexicographical)
    for k in range(2, 10):
        _define_release(tmp_path, f"r26.11.{k}", ["S0"])
    assert tp.next_release_name(tmp_path, date(2026, 11, 19)) == "r26.11.10"
    _define_release(tmp_path, "r26.11.10", ["S0"])
    assert tp.next_release_name(tmp_path, date(2026, 11, 19)) == "r26.11.11"

    # New month resets n to 1
    assert tp.next_release_name(tmp_path, date(2026, 12, 1)) == "r26.12.1"


def test_explicit_release_name_refused_if_exists_or_latest(tmp_path):
    tp = _import_pod_module("train_plan")
    _define_release(tmp_path, "r26.11.1", ["S0"], with_manifest=False)
    m_only = tmp_path / "manifests" / "release=r26.11.2.json"
    m_only.parent.mkdir(parents=True, exist_ok=True)
    m_only.write_text("{}")

    assert tp.resolve_release_name(tmp_path, date(2026, 11, 19), explicit="r26.11.9") == "r26.11.9"
    with pytest.raises(tp.ReleaseTrainError):
        tp.resolve_release_name(tmp_path, date(2026, 11, 19), explicit="r26.11.1")
    with pytest.raises(tp.ReleaseTrainError):
        tp.resolve_release_name(tmp_path, date(2026, 11, 19), explicit="r26.11.2")
    with pytest.raises(tp.ReleaseTrainError):
        tp.resolve_release_name(tmp_path, date(2026, 11, 19), explicit="latest")


# --- Rule (2): stage selection & completeness (DEC-52, DEC-54, C-38a) ---------------------------------------


def test_stage_selection_preserves_prior_order_and_appends_by_inputs_time(tmp_path):
    tp = _import_pod_module("train_plan")
    # Prior stages in release_def order [P2, P1] even though P1 has an earlier timestamp
    _make_stage(tmp_path, "P1", depth="called", started_at="2026-11-01T00:00:00+00:00")
    _make_stage(tmp_path, "P2", depth="called", started_at="2026-11-02T00:00:00+00:00")
    _define_release(tmp_path, "rM.2", ["P2"])
    _define_release(tmp_path, "r26.11.1", ["P2", "P1"])
    # Also add release run dirs to ensure they are never mistaken for stages
    _write_record(tmp_path, "r26.11.1", "catalog")
    _write_record(tmp_path, "rM.1", "cluster.plan")

    # Newly complete stages: B started before A, plus C tied with A (tiebreak by name: A before C)
    _make_stage(tmp_path, "B", depth="called", started_at="2026-11-10T01:00:00+00:00")
    _make_stage(tmp_path, "C", depth="annotated", started_at="2026-11-12T05:00:00+00:00",
                libraries=[("pfam", "38.2")], crispr_enabled=True)
    _make_stage(tmp_path, "A", depth="called", started_at="2026-11-12T05:00:00+00:00")

    sel = tp.select_stages(tmp_path)
    assert sel["latest_release"] == "r26.11.1"
    assert sel["prior_stages"] == ["P2", "P1"]
    assert sel["new_stages"] == ["B", "A", "C"]
    assert sel["stages"] == ["P2", "P1", "B", "A", "C"]
    assert sel["excluded"] == []
    assert sel["nothing_new"] is False


def test_stage_selection_excludes_stored_and_incomplete_stages_with_reasons(tmp_path):
    tp = _import_pod_module("train_plan")
    _make_stage(tmp_path, "S0", depth="called", started_at="2026-11-01T00:00:00+00:00")
    _define_release(tmp_path, "r26.11.1", ["S0"])

    # Valid new stage
    _make_stage(tmp_path, "S1", depth="annotated", started_at="2026-11-02T00:00:00+00:00",
                libraries=[("pfam", "38.2")], crispr_enabled=True)
    # Stored-only stage (DEC-54: excluded)
    _make_stage(tmp_path, "S_STORED", depth="stored", started_at="2026-11-03T00:00:00+00:00")
    # Missing tapes step
    _make_stage(tmp_path, "S_NO_TAPES", depth="called", started_at="2026-11-04T00:00:00+00:00",
                missing_step="tapes")
    # Failed genecall step
    _make_stage(tmp_path, "S_BAD_EXIT", depth="called", started_at="2026-11-05T00:00:00+00:00",
                bad_exit_step="genecall")
    # Unfinished publish intent file
    _make_stage(tmp_path, "S_INTENT", depth="called", started_at="2026-11-06T00:00:00+00:00",
                intent_step="tapes")
    # Size-mismatched output
    _make_stage(tmp_path, "S_CORRUPT", depth="called", started_at="2026-11-07T00:00:00+00:00",
                corrupt_output_step="igr")
    # Annotated stage with Pfam plan started but domains-pfam-38.2.json not yet published
    _make_stage(tmp_path, "S_SCAN_WIP", depth="called", started_at="2026-11-08T00:00:00+00:00")
    _write_record(tmp_path, "S_SCAN_WIP", "domains-pfam-38.2.plan",
                  argv=["python", "-m", "biobox.cli", "stage", "domains", "--library", "pfam",
                        "--library-version", "38.2"])
    # Annotated + CRISPR stage missing the filtered near-cas CRISPR step
    _make_stage(tmp_path, "S_NO_NEAR_CAS", depth="annotated", started_at="2026-11-09T00:00:00+00:00",
                libraries=[("pfam", "38.2")], crispr_enabled=True,
                missing_step=f"crispr-{crispr.near_cas_producer('38.2')}")

    sel = tp.select_stages(tmp_path)
    assert sel["stages"] == ["S0", "S1"]
    assert sel["new_stages"] == ["S1"]
    excl = {e["stage"]: e for e in sel["excluded"]}
    assert set(excl) == {"S_STORED", "S_NO_TAPES", "S_BAD_EXIT", "S_INTENT",
                         "S_CORRUPT", "S_SCAN_WIP", "S_NO_NEAR_CAS"}
    assert "stored" in excl["S_STORED"]["reason"].lower()
    assert any("tapes" in p for p in excl["S_NO_TAPES"]["problems"])
    assert any("exit_code" in p for p in excl["S_BAD_EXIT"]["problems"])
    assert any("publishing.json" in p for p in excl["S_INTENT"]["problems"])
    assert any("missing or size-mismatched" in p for p in excl["S_CORRUPT"]["problems"])
    assert any("domains-pfam-38.2" in p for p in excl["S_SCAN_WIP"]["problems"])
    assert any("near" in p or "cas1kb" in p for p in excl["S_NO_NEAR_CAS"]["problems"])


def test_zero_new_stages_and_fresh_root(tmp_path):
    tp = _import_pod_module("train_plan")
    # Fresh root with one complete stage and no prior release
    _make_stage(tmp_path, "S0", depth="called", started_at="2026-11-01T00:00:00+00:00")
    sel_fresh = tp.select_stages(tmp_path)
    assert sel_fresh["latest_release"] is None
    assert sel_fresh["stages"] == ["S0"]
    assert sel_fresh["new_stages"] == ["S0"]
    assert sel_fresh["nothing_new"] is False

    # After defining a release over S0, zero new stages
    _define_release(tmp_path, "r26.11.1", ["S0"])
    sel_same = tp.select_stages(tmp_path)
    assert sel_same["latest_release"] == "r26.11.1"
    assert sel_same["stages"] == ["S0"]
    assert sel_same["new_stages"] == []
    assert sel_same["nothing_new"] is True


def test_dry_run_cli_prints_plan_and_writes_nothing(tmp_path):
    _make_stage(tmp_path, "S0", depth="called", started_at="2026-11-01T00:00:00+00:00")
    _define_release(tmp_path, "r26.11.1", ["S0"])
    _make_stage(tmp_path, "S1", depth="called", started_at="2026-11-02T00:00:00+00:00")
    _make_stage(tmp_path, "S2", depth="stored", started_at="2026-11-03T00:00:00+00:00")

    before = sorted(p.relative_to(tmp_path).as_posix() for p in tmp_path.rglob("*"))
    prof_path = REPO / "profiles" / "wilkcomp02.yaml"
    res = subprocess.run(
        [sys.executable, str(POD_DIR / "train_plan.py"), "--root", str(tmp_path),
         "--profile", str(prof_path), "--date", "2026-11-19", "--no-live-cpu"],
        capture_output=True, text=True, check=False,
    )
    assert res.returncode == 0, res.stderr
    assert "r26.11.2" in res.stdout
    assert "S0,S1" in res.stdout
    assert "S2" in res.stdout
    assert "40" in res.stdout and "80" in res.stdout  # thread plan by the hour (DEC-74)
    after = sorted(p.relative_to(tmp_path).as_posix() for p in tmp_path.rglob("*"))
    assert before == after  # dry-run writes nothing


def test_dry_run_cli_refuses_existing_explicit_release(tmp_path):
    _make_stage(tmp_path, "S0", depth="called", started_at="2026-11-01T00:00:00+00:00")
    _define_release(tmp_path, "r26.11.1", ["S0"])
    prof_path = REPO / "profiles" / "wilkcomp02.yaml"
    res = subprocess.run(
        [sys.executable, str(POD_DIR / "train_plan.py"), "--root", str(tmp_path),
         "--profile", str(prof_path), "--date", "2026-11-19", "--release", "r26.11.1", "--no-live-cpu"],
        capture_output=True, text=True, check=False,
    )
    assert res.returncode == 2
    assert "already exists" in res.stderr


# --- Rate guard (DEC-43, DEC-53) ----------------------------------------------------------------------------


def test_rate_guard_dec43_and_dec53(tmp_path):
    te = _import_pod_module("train_est")
    rel = "r26.11.2"
    est_path = tmp_path / "est.json"
    est_path.write_text(json.dumps({
        "projection": {
            "factors": {},
            "steps": {
                "tokens": 10.0,       # < 24 CPU-min -> not judged (DEC-43)
                "postings": 30.0,     # >= 24 CPU-min -> judged
                "cluster.plan": 40.0, # >= 24 CPU-min -> judged
            },
        }
    }))
    # Small step off by 3x -> not judged -> exit 0
    _write_record(tmp_path, rel, "tokens", cpu_min=30.0)
    # Judged step >20% UNDER estimate -> logged in under_20pct, exit 0 (DEC-53)
    _write_record(tmp_path, rel, "postings", cpu_min=15.0)
    res_under = te.check_rates(tmp_path, rel, est_path)
    assert res_under["over_20pct"] == []
    assert len(res_under["under_20pct"]) == 1
    assert res_under["under_20pct"][0]["step"] == "postings"
    assert res_under["exit_code"] == 0

    # Judged step over 1.20 and at or under 1.50 of its estimate (55 / 40 = 1.375) -> recorded, exit 0 (DEC-64;
    # before DEC-64 this record stopped the train, and this test asserted exit 3 for it)
    _write_record(tmp_path, rel, "cluster.plan", cpu_min=55.0)
    res_band = te.check_rates(tmp_path, rel, est_path)
    assert res_band["over_20pct"] == []
    assert [x["step"] for x in res_band["recorded_over_20pct"]] == ["cluster.plan"]
    assert res_band["exit_code"] == 0

    # Judged step more than 1.50 times its estimate (61 / 40 = 1.525) -> exit 3 (DEC-43, DEC-64)
    _write_record(tmp_path, rel, "cluster.plan", cpu_min=61.0)
    res_over = te.check_rates(tmp_path, rel, est_path)
    assert len(res_over["over_20pct"]) == 1
    assert res_over["over_20pct"][0]["step"] == "cluster.plan"
    assert res_over["exit_code"] == 3


# --- Post-release check & root_added.tsv + announcement draft -----------------------------------------------


def test_train_check_writes_root_added_and_announcement(tmp_path):
    tc = _import_pod_module("train_check")
    root = tmp_path / "root"
    root.mkdir()
    _make_file(root, "manifests/release=r26.11.1.json", b'{"release":"r26.11.1","stages":["S0"],"files":[],"artifacts":[]}')
    before = tmp_path / "before.tsv"
    before.write_text("manifests/release=r26.11.1.json\t62\t1700000000.0\n")

    # Add files for r26.11.2
    _make_file(root, "manifests/release=r26.11.2.json",
               json.dumps({"release": "r26.11.2", "stages": ["S0", "S1"], "files": [], "artifacts": []}).encode())
    _make_file(root, "runs/stage=r26.11.2/catalog.json", b'{"step":"catalog"}')
    after = tmp_path / "after.tsv"
    after.write_text(
        "manifests/release=r26.11.1.json\t62\t1700000000.0\n"
        f"manifests/release=r26.11.2.json\t{(root / 'manifests/release=r26.11.2.json').stat().st_size}\t1700000100.0\n"
        f"runs/stage=r26.11.2/catalog.json\t{(root / 'runs/stage=r26.11.2/catalog.json').stat().st_size}\t1700000100.0\n"
    )
    out_dir = tmp_path / "out"
    rc = tc.run_check(root=root, before=before, after=after, release="r26.11.2", stages=["S0", "S1"],
                      earlier="r26.11.1", out=out_dir)
    assert rc == 0
    tsv = (out_dir / "root_added.tsv").read_text().splitlines()
    assert tsv[0] == "path\tbytes\tsha256"
    assert len(tsv) == 3
    ann = (out_dir / "announcement.md").read_text()
    assert "r26.11.2" in ann and "S0,S1" in ann

    # Stray file not naming r26.11.2 must fail
    _make_file(root, "manifests/release=r26.11.1.json", b"modified")
    after_bad = tmp_path / "after_bad.tsv"
    after_bad.write_text(
        "manifests/release=r26.11.1.json\t8\t1700000200.0\n"
        f"manifests/release=r26.11.2.json\t{(root / 'manifests/release=r26.11.2.json').stat().st_size}\t1700000100.0\n"
    )
    rc_bad = tc.run_check(root=root, before=before, after=after_bad, release="r26.11.2", stages=["S0", "S1"],
                          earlier="r26.11.1", out=tmp_path / "out_bad")
    assert rc_bad == 1


# --- Card A-25: post-manifest parallelism (DEC-75, DEC-73, DEC-74, C-28) ------------------------------------


def _a25_make_train_world(tmp_path: Path, *, mem_scope: str = "750G", max_day: int = 40,
                          max_night: int = 80, live_allowance: str | None = None) -> tuple[Path, Path, Path, dict[str, str]]:
    """Isolated stubbed release_train.sh environment capturing pod_run.sh and release_examples.sh args/env."""
    import os

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
    (root / "manifests" / "release=r26.11.1.json").write_text("{}\n")

    plan = {
        "release": "r26.11.2",
        "stages": ["S1", "P1"],
        "new_stages": ["P1"],
        "latest_release": "r26.11.1",
        "nothing_new": False,
        "earlier_manifest_releases": ["r26.11.1"],
        "last_stage": {
            "stage": "P1", "depth": "called", "genomes": "g.json", "objects": "o.json", "crispr": False,
        },
    }
    plan_file = tmp_path / "plan.json"
    plan_file.write_text(json.dumps(plan))

    prof = tree / "profiles" / "wilkcomp02.yaml"
    prof.write_text(
        f"max_threads_day: {max_day}\n"
        f"max_threads_night: {max_night}\n"
        f"mem_scope: {mem_scope}\n"
        f"tmp: {host_tmp}\n"
        f"day_start_hour: 8\n"
        f"day_end_hour: 20\n"
    )

    runstep_out = "" if live_allowance is None else live_allowance
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
            '#!/bin/sh\ncase "$*" in\n'
            '  *" release verify "*) echo "release R: verified, manifest hash abc123"; exit 0 ;;\n'
            '  *"cli releases "*) echo "r26.11.2 listed"; exit 0 ;;\n'
            '  *train_plan.py*--json*) cat "$FAKE_PLAN"; exit 0 ;;\n'
            '  *train_plan.py*) echo "Action: cut"; exit 0 ;;\n'
            '  *train_est.py*) echo "{}"; exit 0 ;;\n'
            '  *" pack verify "*) echo "pack verified"; exit 0 ;;\n'
            f'  *"biobox.runstep"*) echo "{runstep_out}"; exit 0 ;;\n'
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
        '#!/bin/sh\nlog=""; unit=""; mem=""; threads=""\n'
        'while [ $# -gt 0 ]; do\n  case "$1" in\n'
        '    --log) log=$2; shift 2 ;;\n    --mem) mem=$2; shift 2 ;;\n'
        '    --threads) threads=$2; shift 2 ;;\n    --*) shift 2 ;;\n'
        '    *) unit=$1; shift; break ;;\n  esac\ndone\n'
        'echo "pod_run unit=$unit mem=$mem threads=$threads cmd=$*" >> "$CAPTURE_LOG"\n'
        'rel=""; prev=""\nfor a in "$@"; do [ "$prev" = "--release" ] && rel=$a; prev=$a; done\n'
        'mkdir -p "$(dirname "$log")"\n'
        '{\n'
        '  [ "$1" = "build-stage" ] && echo "build-stage: stage P1 and release $rel: all 9 steps done"\n'
        '  [ "$1" = "pack" ] && echo "pack 02dec8312fff239cfed68b255e02f60e built"\n'
        '  echo "pod_run: $unit exited 0"\n'
        '} >> "$log"\nexit 0\n'
    )
    (tree / "scripts" / "release_examples.sh").write_text(
        '#!/bin/sh\necho "release_examples rel=$1 THREADS=${THREADS:-UNSET} MEM=${MEM:-UNSET}" >> "$CAPTURE_LOG"\nexit 0\n'
    )

    drv = bin_dir / "release_train.sh"
    drv.write_text((POD_DIR / "release_train.sh").read_text().replace("/tmp/if2453", str(pod)))
    for other in POD_DIR.iterdir():
        link = bin_dir / other.name
        if other.name != "release_train.sh" and not link.exists():
            link.symlink_to(other)
    capture.write_text("")

    env = {
        **{k: v for k, v in os.environ.items()
           if k not in ("PACK_PARTS", "PACK_JOBS", "PACK_MEM", "EX_THREADS", "EX_MEM", "BIOBOX_TRAIN_NOW")},
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
    }
    return drv, capture, work, env


def _a25_parse_capture(capture: Path) -> tuple[str, str]:
    lines = capture.read_text().splitlines()
    pack_line = next(ln for ln in lines if "cmd=pack build " in ln)
    ex_line = next(ln for ln in lines if ln.startswith("release_examples "))
    return pack_line, ex_line


@pytest.mark.covers("A-25", "C-46", "C-28")
@pytest.mark.parametrize("now_iso,expected_threads", [
    ("2026-11-24T14:00:00", 16),  # Tuesday day window: min(16, 40) = 16
    ("2026-11-24T22:00:00", 24),  # Tuesday night window: min(24, 80) = 24
    ("2026-11-28T14:00:00", 24),  # Saturday weekend window: min(24, 80) = 24
])
def test_a25_train_defaults_day_and_night_on_750g_profile(tmp_path, now_iso, expected_threads):
    drv, capture, _, env = _a25_make_train_world(tmp_path, mem_scope="750G", max_day=40, max_night=80)
    env["BIOBOX_TRAIN_NOW"] = now_iso
    got = subprocess.run(["bash", str(drv), "run", "--no-live-cpu"], env=env, capture_output=True, text=True, timeout=60)
    assert got.returncode == 0, got.stdout + got.stderr
    pack_line, ex_line = _a25_parse_capture(capture)
    assert f"mem=128G threads={expected_threads} " in pack_line, pack_line
    assert f"--parts 24 --jobs {expected_threads}" in pack_line, pack_line
    assert f"THREADS={expected_threads} MEM=64G" in ex_line, ex_line


@pytest.mark.covers("A-25", "C-46", "C-28")
@pytest.mark.parametrize("live_allowance,now_iso,expected_threads", [
    ("10", "2026-11-24T22:00:00", 10),  # night, allowance 10 -> min(24, 10) = 10
    ("0", "2026-11-24T14:00:00", 1),    # day, allowance 0 -> at least 1
])
def test_a25_train_defaults_respect_live_allowance(tmp_path, live_allowance, now_iso, expected_threads):
    drv, capture, _, env = _a25_make_train_world(tmp_path, mem_scope="750G", live_allowance=live_allowance)
    env["BIOBOX_TRAIN_NOW"] = now_iso
    got = subprocess.run(["bash", str(drv), "run"], env=env, capture_output=True, text=True, timeout=60)
    assert got.returncode == 0, got.stdout + got.stderr
    pack_line, ex_line = _a25_parse_capture(capture)
    assert f"mem=128G threads={expected_threads} " in pack_line, pack_line
    assert f"--jobs {expected_threads}" in pack_line, pack_line
    assert f"THREADS={expected_threads} MEM=64G" in ex_line, ex_line


@pytest.mark.covers("A-25", "C-46")
@pytest.mark.parametrize("mem_scope,expected_pack_mem,expected_ex_mem", [
    ("48G", "48G", "48G"),
    ("16G", "16G", "16G"),
    ("80G", "80G", "64G"),
])
def test_a25_train_defaults_clamp_mem_to_profile_mem_scope(tmp_path, mem_scope, expected_pack_mem, expected_ex_mem):
    drv, capture, _, env = _a25_make_train_world(tmp_path, mem_scope=mem_scope)
    env["BIOBOX_TRAIN_NOW"] = "2026-11-24T14:00:00"
    got = subprocess.run(["bash", str(drv), "run", "--no-live-cpu"], env=env, capture_output=True, text=True, timeout=60)
    assert got.returncode == 0, got.stdout + got.stderr
    pack_line, ex_line = _a25_parse_capture(capture)
    assert f"mem={expected_pack_mem} " in pack_line, pack_line
    assert f"MEM={expected_ex_mem}" in ex_line, ex_line


@pytest.mark.covers("A-25", "C-46")
def test_a25_train_explicit_overrides_take_precedence(tmp_path):
    drv, capture, _, env = _a25_make_train_world(tmp_path, mem_scope="750G", max_day=40, max_night=80)
    env.update({
        "BIOBOX_TRAIN_NOW": "2026-11-24T22:00:00",
        "PACK_PARTS": "2",
        "PACK_JOBS": "2",
        "PACK_MEM": "8G",
        "EX_THREADS": "2",
        "EX_MEM": "12G",
    })
    got = subprocess.run(["bash", str(drv), "run", "--no-live-cpu"], env=env, capture_output=True, text=True, timeout=60)
    assert got.returncode == 0, got.stdout + got.stderr
    pack_line, ex_line = _a25_parse_capture(capture)
    assert "mem=8G threads=2 " in pack_line, pack_line
    assert "--parts 2 --jobs 2" in pack_line, pack_line
    assert "THREADS=2 MEM=12G" in ex_line, ex_line


@pytest.mark.covers("A-25", "C-28")
@pytest.mark.parametrize("env_extra,expected_cap", [
    ({"LAUNCH_DOW": "2", "LAUNCH_HOUR": "14", "LAUNCH_UID": "1000", "LAUNCH_PS_OUT": ""}, "40"),
    ({"LAUNCH_DOW": "2", "LAUNCH_HOUR": "22", "LAUNCH_UID": "1000", "LAUNCH_PS_OUT": ""}, "80"),
    ({"LAUNCH_DOW": "6", "LAUNCH_HOUR": "14", "LAUNCH_UID": "1000", "LAUNCH_PS_OUT": ""}, "80"),
    # root (>200%), own user (>3000%), and another login user at exactly 0.5 CPUs (50.0%) -> still 80
    (
        {
            "LAUNCH_DOW": "2", "LAUNCH_HOUR": "22", "LAUNCH_UID": "1000",
            "LAUNCH_PS_OUT": "0 200.0\n999 80.0\n1000 3000.0\n2000 50.0\n",
        },
        "80",
    ),
    # another login user > 0.5 CPUs (50.1%) at night -> drops to 40 (C-28)
    (
        {
            "LAUNCH_DOW": "2", "LAUNCH_HOUR": "22", "LAUNCH_UID": "1000",
            "LAUNCH_PS_OUT": "1000 100.0\n2000 50.1\n",
        },
        "40",
    ),
    # two other login users summing > 0.5 CPUs (55.0%) on weekend -> drops to 40 (C-28)
    (
        {
            "LAUNCH_DOW": "7", "LAUNCH_HOUR": "12", "LAUNCH_UID": "1000",
            "LAUNCH_PS_OUT": "2000 30.0\n2001 25.0\n",
        },
        "40",
    ),
    # explicit CAP override always wins
    ({"CAP": "12", "LAUNCH_DOW": "2", "LAUNCH_HOUR": "22", "LAUNCH_UID": "1000", "LAUNCH_PS_OUT": ""}, "12"),
])
def test_a25_launch_cap_day_night_weekend_and_others_busy(env_extra, expected_cap):
    import os

    launch = REPO / "examples" / "_shared" / "launch.sh"
    env = {k: v for k, v in os.environ.items()
           if k not in ("CAP", "LAUNCH_DOW", "LAUNCH_HOUR", "LAUNCH_UID", "LAUNCH_PS_OUT")}
    env.update(env_extra)
    got = subprocess.run(
        ["bash", "-c", f'. "{launch}" && cap'],
        env=env, capture_output=True, text=True, timeout=30,
    )
    assert got.returncode == 0, got.stderr
    assert got.stdout.strip() == expected_cap, f"got {got.stdout!r}, stderr={got.stderr!r}"
