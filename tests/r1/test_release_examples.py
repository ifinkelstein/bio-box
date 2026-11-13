"""`scripts/release_examples.sh`: every worked example run on one release, and one summary (card E-11, AGENTS §5.7).

Nothing here touches the pod. `--dry-run` prints what the script would run; `--check-only` runs only its checking
step, on `observed.tsv` files already in OUT (here fabricated: the examples' committed `results/observed.tsv`, or
hand-written ones), with `LOCAL_PYTHON` in place of the image. A planted failure is an expectation edited in a
scratch copy of `examples/`; the script must then exit non-zero and say which check failed.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("E-11")  # ring R1 by its directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "release_examples.sh"
EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
NAMES = ["crispr_cas_neighbourhoods", "rnap"]


def _env(tree: Path, out: Path, **extra: str) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if k not in ("RELEASE", "OUT", "TREE", "ROOT")}
    env.update({"TREE": str(tree), "OUT": str(out), "ROOT": "/nonexistent/root", "PROFILE": "/nonexistent/p.yaml",
                "IMAGE": "/nonexistent/x.sif", "DIGEST": "sha256:" + "0" * 64, "RELEASE": "rN.3"})
    env.update(extra)
    return env


def _run(args: list[str], env: dict[str, str]) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", str(SCRIPT), *args], env=env, capture_output=True, text=True, timeout=120)


def _scratch_tree(tmp_path: Path) -> Path:
    tree = tmp_path / "tree"
    shutil.copytree(EXAMPLES, tree / "examples", ignore=shutil.ignore_patterns("figures", "__pycache__"))
    return tree


def _fabricate(out: Path, tree: Path) -> None:
    """OUT as a run would leave it: each example's observed.tsv (its committed recorded run's values)."""
    for name in NAMES:
        (out / name).mkdir(parents=True)
        shutil.copy(tree / "examples" / name / "results" / "observed.tsv", out / name / "observed.tsv")


def _summary(out: Path) -> list[list[str]]:
    return [line.split("\t") for line in (out / "summary.tsv").read_text().splitlines()]


def test_the_dry_run_lists_every_example_and_writes_nothing(tmp_path):
    out = tmp_path / "out"
    p = _run(["--dry-run", "rN.3"], _env(REPO, out))
    assert p.returncode == 0, p.stdout + p.stderr
    for name in NAMES:
        assert f"examples/{name}/run.sh" in p.stdout, p.stdout
        assert f"{out}/{name}" in p.stdout  # each example's own OUT
    assert "_shared" not in p.stdout.replace("examples/_shared/check.py", "")  # the shared directory is not an example
    assert "check.py" in p.stdout  # the checking step is listed too
    assert not out.exists()


def test_the_dry_run_refuses_a_release_that_is_not_the_environments(tmp_path):
    p = _run(["--dry-run", "rN.2"], _env(REPO, tmp_path / "out"))
    assert p.returncode == 2 and "RELEASE" in p.stderr


def test_without_a_release_it_prints_its_usage_and_exits_2(tmp_path):
    p = _run([], _env(REPO, tmp_path / "out"))
    assert p.returncode == 2 and "usage" in p.stderr


def test_the_checking_step_passes_on_the_recorded_runs_and_writes_one_summary_and_a_note(tmp_path):
    tree = _scratch_tree(tmp_path)
    out = tmp_path / "out"
    _fabricate(out, tree)
    p = _run(["--check-only", "rN.3"], _env(tree, out, LOCAL_PYTHON=sys.executable))
    assert p.returncode == 0, p.stdout + p.stderr
    rows = _summary(out)
    assert rows[0] == ["example", "check", "kind", "expected", "observed", "verdict"]
    assert {r[0] for r in rows[1:]} == set(NAMES)
    assert all(len(r) == 6 for r in rows)
    assert {r[5] for r in rows[1:]} == {"ok"}
    note = (out / "note.md").read_text()
    assert "rN.3" in note and "scripts/release_examples.sh" in note
    for name in NAMES:
        assert name in note


@pytest.mark.parametrize("example,check_id,old,new", [
    ("rnap", "t7_control_bits", ">= 30", ">= 400"),
    ("crispr_cas_neighbourhoods", "tniq_cas7f_loci_with_cas1_cas2_cas3", "== 0", "== 1"),
])
def test_a_planted_failure_in_a_scratch_copy_makes_the_checking_step_exit_non_zero(tmp_path, example, check_id, old,
                                                                                     new):
    tree = _scratch_tree(tmp_path)
    expected = tree / "examples" / example / "expected.tsv"
    lines = expected.read_text().splitlines(keepends=True)
    hit = [k for k, line in enumerate(lines) if line.startswith(check_id + "\t")]
    assert len(hit) == 1 and lines[hit[0]].rstrip("\n").endswith("\t" + old)
    lines[hit[0]] = lines[hit[0]].rstrip("\n")[: -len(old)] + new + "\n"
    expected.write_text("".join(lines))
    out = tmp_path / "out"
    _fabricate(out, tree)
    p = _run(["--check-only", "rN.3"], _env(tree, out, LOCAL_PYTHON=sys.executable))
    assert p.returncode == 1, p.stdout + p.stderr
    failed = [r for r in _summary(out)[1:] if r[5] != "ok"]
    assert [(r[0], r[1], r[2], r[3], r[5]) for r in failed] == [(example, check_id, "invariant", new, "FAIL")]
    assert "1 invariant" in p.stderr and example in p.stderr


def test_an_example_without_an_observed_file_fails_its_invariants_as_missing(tmp_path):
    tree = _scratch_tree(tmp_path)
    out = tmp_path / "out"
    _fabricate(out, tree)
    (out / "rnap" / "observed.tsv").unlink()  # as when its run.sh stopped before the end
    p = _run(["--check-only", "rN.3"], _env(tree, out, LOCAL_PYTHON=sys.executable))
    assert p.returncode == 1
    rnap = [r for r in _summary(out)[1:] if r[0] == "rnap"]
    assert rnap and {r[5] for r in rnap} == {"missing"}
    assert {r[5] for r in _summary(out)[1:] if r[0] != "rnap"} == {"ok"}


def test_a_malformed_expected_file_exits_2(tmp_path):
    tree = _scratch_tree(tmp_path)
    expected = tree / "examples" / "rnap" / "expected.tsv"
    expected.write_text(expected.read_text().replace("\tinvariant\t>= 30", "\tinvariant\t=> 30", 1))
    out = tmp_path / "out"
    _fabricate(out, tree)
    p = _run(["--check-only", "rN.3"], _env(tree, out, LOCAL_PYTHON=sys.executable))
    assert p.returncode == 2, p.stdout + p.stderr


def test_the_checking_step_refuses_an_out_that_already_has_a_summary(tmp_path):
    tree = _scratch_tree(tmp_path)
    out = tmp_path / "out"
    _fabricate(out, tree)
    env = _env(tree, out, LOCAL_PYTHON=sys.executable)
    assert _run(["--check-only", "rN.3"], env).returncode == 0
    p = _run(["--check-only", "rN.3"], env)
    assert p.returncode == 2 and "summary" in p.stderr
