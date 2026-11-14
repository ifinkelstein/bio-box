"""Ring R2, card E-15 (DEC-40): `build-stage --crispr` runs the filtered CRISPR producer after the stage's Pfam scan.

On the C-04 fixture (lambda and M. genitalium, `test_builders_out.direct`), as card E-12's test does for the default
producer:

- with `--crispr` and a domain library named `pfam`, the step list is card E-12's list plus `crispr-<filtered
  producer>` right after the last domain scan step (`stage crispr --near-cas --library pfam --library-version V`); with
  `--crispr` and no `pfam` library (the C-04 fixture's library is `test-1`) the list is card E-12's, unchanged, and
  the build says in one line that only the default producer runs;
- a build with `--crispr --library pfam 38.2 ...` ends with the default producer's files byte-identical to card
  E-12's build and the filtered producer's files equal to the direct call on the same stage; a second call runs no
  step's command again and changes no file;
- the filtered step's estimate is the default's at `CRISPR_NEAR_CAS_COST` times its rate; above 80% of the budget it is
  refused (exit 2) before it runs, after the default step ran.
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
from test_crispr_build_stage import direct, name, with_arrays  # noqa: E402,F401  (the E-12 fixtures)

pytestmark = [pytest.mark.r2, pytest.mark.covers("E-15", "D-01", "MR-20"),
              pytest.mark.skipif(shutil.which("minced") is None, reason="minced (MinCED) not on PATH; it is in the "
                                                                       "bio-box image, not in GitHub Actions"),
              pytest.mark.skipif(c08.MMSEQS is None, reason="build-stage clusters with MMseqs2 (in the image)")]

PFAM_VERSION = "38.2"


def near_name() -> str:
    from biobox.l1 import crispr

    return f"crispr-{crispr.near_cas_producer(PFAM_VERSION)}"


def with_pfam(direct, root: Path) -> bs.Build:
    w = direct["work"]
    return dataclasses.replace(c04.build_of(direct, root), crispr=True,
                               libraries=(bs.Library("pfam", PFAM_VERSION, w.lib, "gathering"),))


def argv_pfam(direct, root: Path) -> list[str]:
    argv = c04.argv_of(direct, root)
    k = argv.index("--library")
    return argv[:k + 1] + ["pfam", PFAM_VERSION] + argv[k + 3:] + ["--crispr"]


@pytest.fixture(scope="module")
def built(direct, tmp_path_factory):
    root = tmp_path_factory.mktemp("e15-build") / "root"
    assert cli.main(argv_pfam(direct, root)) == 0
    return root


def test_the_filtered_step_follows_the_domain_scan_only_with_a_pfam_library(direct, tmp_path):
    prof = direct["prof"]
    got = bs.steps(prof, with_pfam(direct, tmp_path / "root"))
    names = [s.name for s in got]
    k = c04.NAMES.index("features") + 1
    scan = "domains-pfam-38.2"
    e12 = c04.NAMES[:k] + [name()] + c04.NAMES[k:]
    j = e12.index("domains-test-1")
    assert names == e12[:j] + [scan, near_name()] + e12[j + 1:]
    st = got[names.index(near_name())]
    assert st.check == "crispr-cas" and st.argv[:2] == ("stage", "crispr") and st.scratch is not None
    a = list(st.argv)
    assert "--near-cas" in a and a[a.index("--library") + 1] == "pfam" and a[a.index("--library-version") + 1] == "38.2"
    assert a[a.index("--cpus") + 1] == "1"
    assert all("--producers" not in s.argv for s in got if s.argv is not None)
    # without a pfam library (the fixture's own library, test-1): card E-12's list, unchanged
    e12_build = dataclasses.replace(c04.build_of(direct, tmp_path / "root"), crispr=True)
    assert [s.name for s in bs.steps(prof, e12_build)] == e12
    # a pfam library without --crispr: no CRISPR step at all
    no_crispr = dataclasses.replace(with_pfam(direct, tmp_path / "root"), crispr=False)
    assert not any(s.name.startswith("crispr-") for s in bs.steps(prof, no_crispr))
    assert got == bs.steps(prof, with_pfam(direct, tmp_path / "root"))


def test_without_a_pfam_library_the_build_says_only_the_default_runs(direct, tmp_path, capsys):
    e12_build = dataclasses.replace(c04.build_of(direct, tmp_path / "root"), crispr=True)
    assert bs.build_stage(direct["prof"], e12_build, run=lambda *a, **k: 9) == 9  # stops at its first step
    err = capsys.readouterr().err
    lines = [ln for ln in err.splitlines() if "only the default CRISPR producer" in ln]
    assert len(lines) == 1, err
    with_lib = with_pfam(direct, tmp_path / "root2")
    assert bs.build_stage(direct["prof"], with_lib, run=lambda *a, **k: 9) == 9
    assert "only the default CRISPR producer" not in capsys.readouterr().err


def test_the_build_adds_both_producers_files_and_keeps_the_defaults_bytes(direct, with_arrays, built, tmp_path):
    from biobox.l1 import crispr

    p3, p2 = name().removeprefix("crispr-"), near_name().removeprefix("crispr-")
    mine, e12 = c08.data(c08.tree(built)), c08.data(c08.tree(with_arrays))
    for k in (f"l1/feature/producer={p3}/stage=S0/part-0.parquet", f"l1/crispr/producer={p3}/stage=S0/part-0.parquet",
              f"ledger/stage=S0/step=crispr-{p3}/part-0.parquet"):
        assert mine[k] == e12[k], k  # the default producer: card E-12's bytes
    added = {k: v for k, v in mine.items() if f"producer={p2}" in k or f"step=crispr-{p2}" in k}
    assert sorted(added) == sorted([f"l1/feature/producer={p2}/stage=S0/part-0.parquet",
                                    f"l1/crispr/producer={p2}/stage=S0/part-0.parquet",
                                    f"ledger/stage=S0/step=crispr-{p2}/part-0.parquet"])
    assert (built / "runs/stage=S0" / f"{near_name()}.json").is_file()
    # the direct call on a copy of the built stage without the filtered producer's files gives the same files
    root = tmp_path / "root"
    shutil.copytree(built, root)
    for rel in (f"l1/feature/producer={p2}", f"l1/crispr/producer={p2}", f"ledger/stage=S0/step=crispr-{p2}"):
        shutil.rmtree(root / rel)
    out = tmp_path / "out"
    crispr.run_near_cas(root, "S0", 1, library_version=PFAM_VERSION, tmp=tmp_path / "scratch", out=out,
                        seq_stages=["S0"])
    assert c08.tree(out) == added


def test_a_second_call_runs_no_command_and_changes_no_file(direct, built, tmp_path):
    before = c08.tree(built)
    seen = []

    def counting(prof, root, stage, step, cmd, **kw):
        seen.append(step)
        return c04.runner(prof, root, stage, step, cmd, **kw)

    assert bs.build_stage(direct["prof"], with_pfam(direct, built), run=counting) == 0
    assert c08.tree(built) == before
    assert near_name() in seen  # handed to run-step, which skips it: its run record and outputs exist


def test_the_filtered_step_above_the_budget_is_refused_after_the_default_ran(direct, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(bs, "CRISPR_NEAR_CAS_COST", 1e9)  # the default's estimate as before, the filtered one huge
    root = tmp_path / "root"
    assert bs.build_stage(direct["prof"], with_pfam(direct, root), run=c04.runner) == runstep.EXIT_REFUSED
    err = capsys.readouterr().err
    assert f"refused before S0.{near_name()}" in err and "no shards" in err
    recs = c04.record_names(root)
    assert f"stage=S0/{name()}" in recs and "stage=S0/domains-pfam-38.2" in recs
    assert f"stage=S0/{near_name()}" not in recs
    p2 = near_name().removeprefix("crispr-")
    assert not (root / f"l1/crispr/producer={p2}").exists()
