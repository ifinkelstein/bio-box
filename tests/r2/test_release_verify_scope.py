"""C-19: an older release still verifies after a library is scanned on its stages.

Found on wilkcomp02 (2026-11-11, `planning/measurements/2026-11-11_build_r50.md`): releases built before Pfam was
scanned on their stages reported the scan's domain and ledger files as `not in manifest`, and `release verify` exited
1, although nothing they list had changed. Adding a library's files to a stage is allowed (AGENTS §3.10 forbids
rewriting, not adding).

The store, built by `build-stage` (as `tests/r2/test_build_stage.py` runs it, with the C-08 fixture's inputs):

1. stage S0 and release A over it, with no library;
2. the library scanned on S0 and release B over S0, with it;
3. stage S1 and release C over S0 and S1, with no library given: C's tokens hold the library's families of S0 only;
4. the library scanned on S1 and release D over S0 and S1.

Then: A and C verify with no problem, and so do B and D; C's verify warns, exit unchanged, that the library has domain
files for S0 but not S1; every check a release owns still fails as before (a listed file changed or missing, a second
part beside a listed one, a file no table claims), and a stray in the library's partition is reported by the release
that lists that partition, not by the one that does not.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

from biobox import catalog, cli, runstep
from biobox import stage as bs

sys.path.insert(0, str(Path(__file__).parent))
import test_builders_out as c08  # noqa: E402

pytestmark = [
    pytest.mark.r2,
    pytest.mark.skipif(c08.MMSEQS is None, reason="build-stage clusters with MMseqs2 (in the image)"),
]

A, B, C, D = "rS0.c19a", "rS0.c19b", "rS1.c19c", "rS1.c19d"
LIB = "test-1"  # the library's partition value: <name>-<version>


def runner(prof, root, stage, step, cmd, **kw) -> int:
    """run-step at a fixed hour and load, as `test_build_stage.runner` runs it."""
    return runstep.run_step(prof, root, stage, step, cmd, clock=lambda: c08.NOON_MONDAY, load=lambda: 0.0, **kw)


@pytest.fixture(scope="module")
def store(tmp_path_factory):
    d = tmp_path_factory.mktemp("c19")
    w = c08.Work(d / "work")
    w.prepare()
    prof = c08.profile_at(d / "prof")
    root = d / "root"

    def build(stage: str, release: str, stages: tuple[str, ...], library: bool) -> None:
        libs = (bs.Library("test", "1", w.lib, "gathering"),) if library else ()
        b = bs.Build(root=root, stage=stage, release=release, stages=stages, genomes=w.genomes(stage),
                     objects=w.objects(stage), libraries=libs, threads=1)
        assert bs.build_stage(prof, b, run=runner, clock=lambda: c08.NOON_MONDAY, load=lambda: 0.0) == 0

    build("S0", A, ("S0",), library=False)
    c08.write_library(root, w.lib)  # HMMs of S0's longest lambda proteins, from the stage just called
    build("S0", B, ("S0",), library=True)
    build("S1", C, ("S0", "S1"), library=False)
    build("S1", D, ("S0", "S1"), library=True)
    return root


def manifest(root: Path, release: str) -> dict:
    return catalog.read_manifest(catalog.manifest_path(root, release))


def run_verify(root: Path, release: str, capsys) -> tuple[int, str, str]:
    code = cli.main(["release", "verify", "--root", str(root), "--release", release])
    got = capsys.readouterr()
    return code, got.out, got.err


def domain_stages(man: dict) -> list[tuple[str, str]]:
    return sorted({(f["partitions"]["library"], f["partitions"]["stage"]) for f in man["files"]
                   if f["table"] == "l1_domains"})


def test_c19_the_store_is_as_described(store):
    assert domain_stages(manifest(store, A)) == []
    assert domain_stages(manifest(store, B)) == [(LIB, "S0")]
    assert domain_stages(manifest(store, C)) == [(LIB, "S0")]
    assert domain_stages(manifest(store, D)) == [(LIB, "S0"), (LIB, "S1")]
    for s in ("S0", "S1"):  # what the later scans added under the earlier releases' stages
        assert (store / f"l1/domains/library={LIB}/stage={s}/part-0.parquet").is_file()
        assert (store / f"ledger/stage={s}/step=domains-{LIB}/part-0.parquet").is_file()


@pytest.mark.covers("LY-10")
def test_c19_a_release_built_before_a_library_was_scanned_on_its_stage_verifies(store, capsys):
    for release in (A, B, C, D):
        assert catalog.verify_manifest(manifest(store, release), store) == [], release
        code, out, _ = run_verify(store, release, capsys)
        assert (code, "not in manifest" in out) == (0, False), (release, out)
        assert f"release {release}: verified" in out


@pytest.mark.covers("LY-10")
def test_c19_a_release_built_before_one_of_its_stages_was_scanned_warns_naming_the_stages(store, capsys):
    code, _, err = run_verify(store, C, capsys)
    assert code == 0
    lines = [ln for ln in err.splitlines() if f"library {LIB}" in ln]
    assert len(lines) == 1 and lines[0].startswith("warning: ")
    assert "S0" in lines[0] and "S1" in lines[0], lines[0]
    assert catalog.partial_libraries(manifest(store, C)) == [(LIB, ["S0"], ["S1"])]
    for release in (A, B, D):  # no library, or every stage scanned: no such warning
        code, _, err = run_verify(store, release, capsys)
        assert code == 0 and "library " not in err, (release, err)
        assert catalog.partial_libraries(manifest(store, release)) == []


@pytest.mark.covers("LY-10")
def test_c19_what_an_older_release_owns_is_still_checked(store, tmp_path):
    c = tmp_path / "c"
    shutil.copytree(store, c, symlinks=True)
    m_a, m_b = manifest(c, A), manifest(c, B)
    # a second part beside a file A lists, in the stage the library was scanned on
    occ = next(f["path"] for f in m_a["files"] if f["table"] == "protein_occurrence")
    second = occ.replace("part-0.parquet", "part-9.parquet")
    shutil.copy(c / occ, c / second)
    # a file no table claims, under a directory A lists
    unclaimed = occ.replace("part-0.parquet", "extra/part-0.parquet")
    (c / unclaimed).parent.mkdir()
    shutil.copy(c / occ, c / unclaimed)
    # a stray in the library's partition: B lists that partition, A does not
    dom = f"l1/domains/library={LIB}/stage=S0/part-0.parquet"
    lib_stray = dom.replace("part-0.parquet", "part-9.parquet")
    shutil.copy(c / dom, c / lib_stray)
    want = sorted([second, unclaimed])
    assert catalog.verify_manifest(m_a, c) == [f"not in manifest: {p}" for p in want]
    assert catalog.verify_manifest(m_b, c) == [f"not in manifest: {p}" for p in sorted(want + [lib_stray])]
    # a listed file changed, and one removed
    changed = next(f["path"] for f in m_a["files"] if f["table"] == "contig_index")
    data = bytearray((c / changed).read_bytes())
    data[len(data) // 2] ^= 0xFF
    (c / changed).write_bytes(bytes(data))
    gone = next(f["path"] for f in m_a["files"] if f["table"] == "token_dict")
    (c / gone).unlink()
    got = catalog.verify_manifest(m_a, c)
    assert f"sha256 differs: {changed}" in got and f"missing: {gone}" in got
    assert cli.main(["release", "verify", "--root", str(c), "--release", A]) == 1
