"""Ring R2: end-to-end retire, verify, reader refusals, and tar restore on the micro-S0 fixture (card C-54)."""

from __future__ import annotations

import hashlib
import shutil
import sys
import tarfile
from pathlib import Path

import pytest

from biobox import catalog, cli, runstep
from biobox import stage as bs
from biobox.l2 import cluster

sys.path.insert(0, str(Path(__file__).parent))
import test_builders_out as c08  # noqa: E402

pytestmark = [
    pytest.mark.r2,
    pytest.mark.covers("C-54", "LY-10", "R-09"),
    pytest.mark.skipif(c08.MMSEQS is None, reason="build-stage clusters with MMseqs2 (in the image)"),
]


def _runner(prof, root, stage, step, cmd, **kw) -> int:
    return runstep.run_step(prof, root, stage, step, cmd, clock=lambda: c08.NOON_MONDAY, load=lambda: 0.0, **kw)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _tar_release(root: Path, release: str, tar_path: Path) -> str:
    man = catalog.read_manifest(catalog.manifest_path(root, release))
    with tarfile.open(tar_path, "w") as tf:
        for f in man["files"]:
            if f["partitions"].get("release") == release:
                tf.add(root / f["path"], arcname=f["path"])
        pdir = root / "l2" / "cluster_plan" / f"release={release}"
        if pdir.is_dir():
            for p in sorted(pdir.glob("*.json")):
                tf.add(p, arcname=p.relative_to(root).as_posix())
    return _sha256(tar_path)


@pytest.fixture(scope="module")
def fixture_store(tmp_path_factory):
    d = tmp_path_factory.mktemp("c54_fix")
    w = c08.Work(d / "work")
    w.prepare()
    prof = c08.profile_at(d / "prof")
    root = d / "root"

    def build(stage: str, release: str, stages: tuple[str, ...]) -> None:
        b = bs.Build(
            root=root,
            stage=stage,
            release=release,
            stages=stages,
            genomes=w.genomes(stage),
            objects=w.objects(stage),
            libraries=(),
            threads=1,
        )
        assert bs.build_stage(prof, b, run=_runner, clock=lambda: c08.NOON_MONDAY, load=lambda: 0.0) == 0

    build("S0", "r26.11.1", ("S0",))
    build("S1", "r26.11.2", ("S0", "S1"))
    build("S1", "r26.11.3", ("S0", "S1"))
    return root, d / "prof" / "profile.yaml"


def test_c54_full_retire_refusals_and_restore_round_trip(fixture_store, tmp_path, capsys):
    src_root, prof_path = fixture_store
    root = tmp_path / "root"
    shutil.copytree(src_root, root, symlinks=True)

    # Regenerate catalog for copied root so catalog_info.root points to this root
    for rel in ("r26.11.1", "r26.11.2", "r26.11.3"):
        man = catalog.read_manifest(catalog.manifest_path(root, rel))
        cat_path = catalog.catalog_path(root, rel)
        catalog.generate_catalog(man, out_path=cat_path, root=root)

    # Grab a real protein_uid and accession from r26.11.1 before retiring
    con = catalog.open_catalog(catalog.catalog_path(root, "r26.11.1"))
    try:
        puid_hex = con.execute("SELECT lower(hex(protein_uid)) FROM protein_seq LIMIT 1").fetchone()[0]
        acc = con.execute("SELECT accession FROM contig_label LIMIT 1").fetchone()[0]
    finally:
        con.close()

    # Build a search pack for r26.11.1 and r26.11.3
    packs = tmp_path / "packs"
    assert cli.main(["pack", "build", "--root", str(root), "--release", "r26.11.1", "--pack-dir", str(packs)]) == 0
    assert cli.main(["pack", "build", "--root", str(root), "--release", "r26.11.3", "--pack-dir", str(packs)]) == 0
    capsys.readouterr()

    # Backup r26.11.1 to a tar
    tar_path = tmp_path / "added_r26_11_1.tar"
    tar_sha = _tar_release(root, "r26.11.1", tar_path)
    h_before = catalog.manifest_hash(catalog.read_manifest(catalog.manifest_path(root, "r26.11.1")))

    # Retire r26.11.1
    assert cli.main([
        "release", "retire", "--root", str(root), "--release", "r26.11.1",
        "--tar", f"{tar_path}={tar_sha}",
        "--copy", f"tower:/mnt/user/research-bulk/bio-box/root_backups/added_r26_11_1.tar={tar_sha}",
    ]) == 0
    capsys.readouterr()

    # 1. Verify r26.11.1 (archived) and r26.11.2/3 (--previous r26.11.1)
    assert cli.main(["release", "verify", "--root", str(root), "--release", "r26.11.1"]) == 0
    out_v = capsys.readouterr().out
    assert f"release r26.11.1: verified (archived:" in out_v
    assert f"manifest hash {h_before}" in out_v
    assert cli.main(["release", "verify", "--root", str(root), "--release", "r26.11.2", "--previous", "r26.11.1"]) == 0
    assert cli.main(["release", "verify", "--root", str(root), "--release", "r26.11.3", "--previous", "r26.11.2"]) == 0
    capsys.readouterr()

    # 2. B-49 skips archived release
    assert cluster._source_problem(root, "r26.11.4", "r26.11.1", {}) == "release is archived"

    # 3. Every reader command on archived r26.11.1 exits 2 cleanly naming archived/release=r26.11.1.json and tar -C
    recipe = tmp_path / "recipe.yaml"
    recipe.write_text(
        f"recipe: test_c54\nwritten_against: r26.11.1\nseed:\n  - protein_uid: {puid_hex}\n"
        f"target: unique\ntool: hmmsearch\niterations:\n  max: 1\n  inclusion_bits: 25.0\n"
        f"context:\n  window: 5g\ncontrols:\n  positive:\n    - protein_uid: {puid_hex}\n"
    )
    reader_cmds = [
        ["sql", "--root", str(root), "--release", "r26.11.1", "-e", "SELECT 1"],
        ["find", "--root", str(root), "--release", "r26.11.1", "--accession", acc],
        ["seq", "--root", str(root), "--release", "r26.11.1", "--protein", puid_hex],
        ["export", "proteins", "--root", str(root), "--release", "r26.11.1"],
        ["search", "run", str(recipe), "--root", str(root), "--release", "r26.11.1", "--wait",
         "--pack-dir", str(packs), "--runs", str(tmp_path / "runs"), "--work", str(tmp_path / "swork")],
        ["nbr", "--root", str(root), "--release", "r26.11.1", "--anchors", puid_hex],
    ]
    for argv in reader_cmds:
        rc = cli.main(argv)
        err = capsys.readouterr().err
        assert rc == 2, (argv, rc, err)
        assert "archived/release=r26.11.1.json" in err, (argv, err)
        assert "tar -C" in err, (argv, err)

    # 4. Full restore from tar + remove archived/release=r26.11.1.json -> unarchived verify, sql, seq work again!
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(root)
    (root / "archived" / "release=r26.11.1.json").unlink()

    assert cli.main(["release", "verify", "--root", str(root), "--release", "r26.11.1"]) == 0
    out_restored = capsys.readouterr().out
    assert "archived:" not in out_restored
    assert f"manifest hash {h_before}" in out_restored

    assert cli.main(["sql", "--root", str(root), "--release", "r26.11.1", "-e", "SELECT count(*) FROM contig_dim"]) == 0
    assert cli.main(["seq", "--root", str(root), "--release", "r26.11.1", "--protein", puid_hex]) == 0
    capsys.readouterr()
