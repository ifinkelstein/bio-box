"""Ring R2, card E-12: `build-stage --crispr` runs the CRISPR array producer as a step after gene calling.

On the C-04 fixture (lambda and M. genitalium, `test_builders_out.direct`), as card E-07's test does for ncRNA:

- the step list with `--crispr` is the list without it plus `crispr-<producer>` after the `features` slot, off by
  default; the IGR, tape, token and posting steps get no `--producers` (the arrays are not on the tape: the tape
  version is the one without `--crispr`);
- a build with `--crispr` ends with every file of the build without it byte-identical (the manifest and catalog
  apart, which list the new files), plus the producer's files, which equal the direct producer's on the same stage;
- above 80% of the budget at the measured rate the step is refused (exit 2) before it runs, since it has no shards.
"""

from __future__ import annotations

import dataclasses
import shutil
import sys
from pathlib import Path

import pytest

from biobox import cli, runstep
from biobox import stage as bs

sys.path.insert(0, str(Path(__file__).parent))
import test_build_stage as c04  # noqa: E402
import test_builders_out as c08  # noqa: E402

direct = c08.direct  # the module-scoped direct build of S0 and release c08.REL (a fixture)

pytestmark = [pytest.mark.r2, pytest.mark.covers("E-12", "D-01", "MR-20"),
              pytest.mark.skipif(shutil.which("minced") is None, reason="minced (MinCED) not on PATH; it is in the "
                                                                       "bio-box image, not in GitHub Actions"),
              pytest.mark.skipif(c08.MMSEQS is None, reason="build-stage clusters with MMseqs2 (in the image)")]


def name() -> str:
    from biobox.l1 import crispr

    return f"crispr-{crispr.producer()}"


def with_crispr(direct, root: Path) -> bs.Build:
    return dataclasses.replace(c04.build_of(direct, root), crispr=True)


@pytest.fixture(scope="module")
def with_arrays(direct, tmp_path_factory):
    root = tmp_path_factory.mktemp("e12-build") / "root"
    assert cli.main(c04.argv_of(direct, root) + ["--crispr"]) == 0
    return root


def test_the_step_follows_the_features_slot_and_passes_no_producers(direct, tmp_path):
    prof = direct["prof"]
    got = bs.steps(prof, with_crispr(direct, tmp_path / "root"))
    names = [s.name for s in got]
    k = c04.NAMES.index("features") + 1
    assert names == c04.NAMES[:k] + [name()] + c04.NAMES[k:]
    st = got[k]
    assert st.check == "crispr" and st.argv[:2] == ("stage", "crispr") and st.scratch is not None
    assert st.argv[st.argv.index("--cpus") + 1] == "1"  # at most --threads (1 here)
    assert all("--producers" not in s.argv for s in got if s.argv is not None)
    assert [s.name for s in bs.steps(prof, c04.build_of(direct, tmp_path / "root"))] == c04.NAMES  # off by default
    assert got == bs.steps(prof, with_crispr(direct, tmp_path / "root"))


def test_with_crispr_the_build_adds_the_producers_files_and_changes_no_other(direct, with_arrays, tmp_path):
    without, withc = c08.data(c08.tree(direct["root"])), c08.data(c08.tree(with_arrays))
    added = {k: v for k, v in withc.items() if k not in without}
    assert {k: withc[k] for k in without} == without
    p = name().removeprefix("crispr-")
    assert sorted(added) == sorted([f"l1/feature/producer={p}/stage=S0/part-0.parquet",
                                    f"l1/crispr/producer={p}/stage=S0/part-0.parquet",
                                    f"ledger/stage=S0/step=crispr-{p}/part-0.parquet"])
    assert (with_arrays / "runs/stage=S0" / f"{name()}.json").is_file()
    from biobox.l1 import crispr

    out = tmp_path / "out"
    crispr.run(direct["root"], "S0", 1, tmp=tmp_path / "scratch", out=out)
    assert c08.tree(out) == added


def test_a_step_estimated_above_the_budget_is_refused_before_it_runs(direct, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(bs, "CRISPR_CPU_S_PER_GBP", 1e9)  # the fixture's 0.6 Mbp far above 80% of the budget
    root = tmp_path / "root"
    assert bs.build_stage(direct["prof"], with_crispr(direct, root), run=c04.runner) == runstep.EXIT_REFUSED
    err = capsys.readouterr().err
    assert f"refused before S0.{name()}" in err and "no shards" in err
    assert not (root / "l1/crispr").exists()
    assert "stage=S0/genecall" in c04.record_names(root) and not any(
        n.startswith("stage=S0/crispr-") for n in c04.record_names(root))
