"""Card C-34: what `release verify` of a published release says when one of its stages gains `contig_label` files.

Stages S1, P1 and P2 of the data root were ingested before card B-27 and have no `contig_label` file; `bio-box stage
labels` (card C-34) writes `l0seq/contig_label/stage=<S>/` for them as new files. A release published before that
(`rN.1`, `rN.2`, `rN.3`, `rM.2`) does not list those files. Before the card `verify_manifest` reported each one as
`not in manifest: l0seq/contig_label/stage=<S>/part-0.parquet` in every release that lists the stage: `contig_label`
has no partition but its stage, so card C-19's `_owner_key` gave it none and every release listing the stage owned it.

The rule (card C-34, `catalog.STAGE_ADDED_TABLES`): a `contig_label` file of stage S is a release's concern exactly
when the release's manifest lists a `contig_label` file of S. What tells the two cases apart is the manifest alone:

- it lists none for S (the release was built before S had labels; its readers take S's labels from S's
  `contig_occurrence`, card W-36): every `contig_label` file of S was added after the release and is ignored;
- it lists one or more for S: those are checked as every listed file is (missing, size, sha256), and any other file of
  `contig_label/stage=S/` is `not in manifest`, exactly as before the card.

Nothing else changes: a changed or missing `contig_occurrence` (or any listed file) of S is reported as before, and a
path no table claims is reported. The general case (any stage table added to a stage after a release) is not ruled:
the code has no notion of it beyond C-19's partitions, so the rule names `contig_label` alone.

The world here is small and built in the test: three stages ingested from FASTA written by the test (S0 with its
label file, S1 and S2 ingested with the label writer left out, as ingest worked before card B-27), release A over S0
and S1 (lists S0's labels only) and release B over S1 and S2 (lists none). The label files then given to S1 and S2 are
the ones a re-ingest of the same inputs writes (the oracle of `tests/r2/test_c34_stage_labels.py` shows that
`stage labels` writes the same bytes), copied into place: this test needs no code of the card to place them.
"""

from __future__ import annotations

import random
import shutil
from pathlib import Path

import pytest

from biobox import catalog, cli, fasta, io
from biobox.l0seq import SourceObject, ingest

pytestmark = pytest.mark.covers("C-34", "LY-10")  # ring R1 by its directory (tests/conftest.py)

A, B, C = "rC34.a", "rC34.b", "rC34.c"
LABEL = "l0seq/contig_label"


def _fasta(path: Path, seed: int, n: int) -> Path:
    rng = random.Random(seed)
    lines = []
    for i in range(n):
        seq = "".join(rng.choice("ACGT") for _ in range(400 + 37 * i)).encode()
        lines += fasta.format_record(f"c34_{seed}_{i} a test contig", seq, 70)
    path.write_bytes(b"\n".join(lines) + b"\n")
    return path


def _ingest(root: Path, stage: str, objects, *, labels: bool) -> None:
    """Ingest `objects`; without `labels` the `contig_label` writer is left out (as ingest worked before card B-27)."""
    mp = pytest.MonkeyPatch()
    real = io.write_parts

    def write_parts(batches, spec, base_dir, *a, **kw):
        if spec.name == "contig_label" and not labels:
            return []
        return real(batches, spec, base_dir, *a, **kw)

    mp.setattr(io, "write_parts", write_parts)
    try:
        ingest(objects, root, stage)
    finally:
        mp.undo()


def _objects(src: Path, stage: str) -> list[SourceObject]:
    seed = {"S0": 1, "S1": 2, "S2": 3}[stage]
    return [SourceObject(_fasta(src / f"{stage}.fna", seed, 3), "c34", "1", "test", True)]


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    d = tmp_path_factory.mktemp("c34r1")
    src = d / "src"
    src.mkdir()
    root = d / "root"
    for s in ("S0", "S1", "S2"):
        _ingest(root, s, _objects(src, s), labels=(s == "S0"))
    for rel, stages in ((A, ["S0", "S1"]), (B, ["S1", "S2"])):
        catalog.define_release(root, rel, stages)
        catalog.write_manifest(catalog.build_manifest(root, rel, stages))
    # the label files a re-ingest of the same inputs writes, for the two stages ingested without them
    oracle = d / "oracle"
    for s in ("S1", "S2"):
        ingest(_objects(src, s), oracle / s, s)
    return {"root": root, "oracle": oracle}


def _manifest(root: Path, rel: str) -> dict:
    return catalog.read_manifest(catalog.manifest_path(root, rel))


def _label_stages(man: dict) -> set[str]:
    return {f["partitions"]["stage"] for f in man["files"] if f["table"] == "contig_label"}


def _give_labels(root: Path, oracle: Path, stage: str) -> list[str]:
    """Copy the oracle's label files of `stage` into the root (new files only); their root-relative paths."""
    out = []
    for p in io.part_files(oracle / stage / f"{LABEL}/stage={stage}"):
        dst = root / f"{LABEL}/stage={stage}" / p.name
        assert not dst.exists()
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(p, dst)
        out.append(dst.relative_to(root).as_posix())
    return out


@pytest.fixture
def labelled(world, tmp_path):
    """A copy of the root in which S1 and S2 have their label files (added after A and B were published)."""
    c = tmp_path / "root"
    shutil.copytree(world["root"], c, symlinks=True)
    added = _give_labels(c, world["oracle"], "S1") + _give_labels(c, world["oracle"], "S2")
    return c, added


def test_c34_the_world_is_as_described(world):
    root = world["root"]
    assert _label_stages(_manifest(root, A)) == {"S0"}
    assert _label_stages(_manifest(root, B)) == set()
    for s in ("S1", "S2"):
        assert not io.part_files(root / f"{LABEL}/stage={s}")
        assert io.part_files(world["oracle"] / s / f"{LABEL}/stage={s}")
    for rel in (A, B):
        assert catalog.verify_manifest(_manifest(root, rel), root) == [], rel


def test_c34_a_release_ignores_the_label_files_a_stage_gained_after_it(labelled, capsys):
    """The new label files of S1 and S2 are not the concern of A (which lists S0's labels, not S1's) nor of B (which
    lists none): both verify, by `verify_manifest` and by `release verify`, as before the files existed."""
    c, added = labelled
    assert sorted(added) == [f"{LABEL}/stage=S1/part-0.parquet", f"{LABEL}/stage=S2/part-0.parquet"]
    for rel in (A, B):
        assert catalog.verify_manifest(_manifest(c, rel), c) == [], rel
        capsys.readouterr()
        assert cli.main(["release", "verify", "--root", str(c), "--release", rel]) == 0, rel
        out = capsys.readouterr().out
        assert f"release {rel}: verified" in out and "not in manifest" not in out, out


def test_c34_what_a_release_lists_of_a_stage_is_still_checked(labelled):
    """The rule is narrow: A lists S0's label file, so a second part beside it, a change to it or its loss is reported
    exactly as before; a changed `contig_occurrence` of S1 and a path no table claims under S1's label directory are
    reported as before too."""
    c, _ = labelled
    m_a, m_b = _manifest(c, A), _manifest(c, B)
    s0 = f"{LABEL}/stage=S0/part-0.parquet"
    assert s0 in {f["path"] for f in m_a["files"]}
    second = f"{LABEL}/stage=S0/part-9.parquet"
    shutil.copy(c / s0, c / second)
    unclaimed = f"{LABEL}/stage=S1/extra/part-0.parquet"
    (c / unclaimed).parent.mkdir()
    shutil.copy(c / s0, c / unclaimed)
    assert catalog.verify_manifest(m_a, c) == [f"not in manifest: {p}" for p in sorted([second, unclaimed])]
    assert catalog.verify_manifest(m_b, c) == [f"not in manifest: {unclaimed}"]  # B lists no S0: not its file
    occ = next(f["path"] for f in m_a["files"] if f["table"] == "contig_occurrence" and f["partitions"]["stage"] == "S1")
    data = bytearray((c / occ).read_bytes())
    data[len(data) // 2] ^= 0xFF
    (c / occ).write_bytes(bytes(data))
    (c / s0).unlink()
    got = catalog.verify_manifest(m_a, c)
    assert f"sha256 differs: {occ}" in got and f"missing: {s0}" in got, got
    assert f"sha256 differs: {occ}" in catalog.verify_manifest(m_b, c)
    assert cli.main(["release", "verify", "--root", str(c), "--release", A]) == 1


def test_c34_a_release_that_lists_a_stage_s_labels_still_owns_them(labelled):
    """Release C, defined after S1 gained its labels, lists them: for C a second part beside S1's label file and the
    loss of that file are reported; for A (published before) neither is its concern."""
    c, _ = labelled
    catalog.define_release(c, C, ["S0", "S1"])
    catalog.write_manifest(catalog.build_manifest(c, C, ["S0", "S1"]))
    m_c, m_a = _manifest(c, C), _manifest(c, A)
    assert _label_stages(m_c) == {"S0", "S1"}
    assert catalog.verify_manifest(m_c, c) == [] and catalog.verify_manifest(m_a, c) == []
    s1 = f"{LABEL}/stage=S1/part-0.parquet"
    second = f"{LABEL}/stage=S1/part-1.parquet"
    shutil.copy(c / s1, c / second)
    assert catalog.verify_manifest(m_c, c) == [f"not in manifest: {second}"]
    assert catalog.verify_manifest(m_a, c) == []
    (c / s1).unlink()
    assert catalog.verify_manifest(m_c, c) == [f"missing: {s1}", f"not in manifest: {second}"]
    assert catalog.verify_manifest(m_a, c) == []
