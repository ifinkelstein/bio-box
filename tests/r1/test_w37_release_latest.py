"""Card W-37 (DEC-52): `--release latest` names, at each call, the release with the greatest name that has a manifest
(and, for a command that searches, a verified search pack). Ring R1: roots made here with manifest files written by
hand, so no release is built.

The order is DEC-52's, written out by hand below as the oracle: names of the scheme `r<YY>.<MM>.<n>` compare by year,
month and number as numbers (`r26.11.10` after `r26.11.9`); every other name (`rN.*`, `rM.*`, `rS<k>.<n>`, a test's
own) sorts before all of them, among themselves by name, since a manifest records no order of building. A build
command refuses `latest`, and a release named `latest` cannot be defined. Resolving writes nothing (AGENTS §3.10: no
pointer file)."""

from __future__ import annotations

import hashlib
import json
import random
from pathlib import Path

import pytest

from biobox import catalog, cli

pytestmark = pytest.mark.covers("W-37")

# DEC-52's order, by hand: the older names by name, then the new scheme by (year, month, number) as numbers
OLD = ["rM.1", "rM.2", "rN.1", "rN.1a", "rN.2", "rN.3", "rS0.1", "rW37.a", "r26.11", "r26.11.1a"]
OLD_BY_NAME = ["r26.11", "r26.11.1a", "rM.1", "rM.2", "rN.1", "rN.1a", "rN.2", "rN.3", "rS0.1", "rW37.a"]
NEW = ["r26.11.1", "r26.11.2", "r26.11.9", "r26.11.10", "r26.12.1", "r27.01.1", "r27.1.2", "r27.10.3"]


def _tree(root: Path) -> dict[str, str]:
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else "dir"
            for p in sorted(root.rglob("*"))}


def _manifest(root: Path, name: str) -> None:
    """A manifest file of a release with no file: enough for the resolution and for `bio-box releases`."""
    p = catalog.manifest_path(root, name)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"release_id": name, "stages": ["S0"], "files": [], "artifacts": []}))


def _root(tmp_path: Path, manifests, defined=()) -> Path:
    root = tmp_path / "root"
    root.mkdir()
    for name in manifests:
        _manifest(root, name)
    for name in defined:  # a release_def and no manifest
        catalog.define_release(root, name, ["S0"])
    return root


def _cli(capsys, *argv) -> tuple[int, str, str]:
    code = cli.main(list(argv))
    got = capsys.readouterr()
    return code, got.out, got.err


def test_w37_release_order_is_dec52s_names_by_number_after_every_other_name():
    names = OLD + NEW
    random.Random(37).shuffle(names)
    assert sorted(names, key=catalog.release_order) == OLD_BY_NAME + NEW
    # numbers, not text: r26.11.10 after r26.11.9, r26.12.1 after both, r27.1.2 after r27.01.1
    for lo, hi in (("r26.11.9", "r26.11.10"), ("r26.11.10", "r26.12.1"), ("r27.01.1", "r27.1.2"), ("rN.3", "r26.11.1")):
        assert catalog.release_order(lo) < catalog.release_order(hi), (lo, hi)


def test_w37_latest_is_the_greatest_name_with_a_manifest_and_resolving_writes_nothing(tmp_path):
    root = _root(tmp_path, ["rN.3", "rM.2", "r26.11.9", "r26.11.10"], defined=["r26.12.1"])
    before = _tree(root)
    assert catalog.latest_release(root) == "r26.11.10"  # r26.12.1 has a release_def only: not a candidate
    # a condition that only some releases meet: the greatest of those
    assert catalog.latest_release(root, lambda r: None if r in ("r26.11.9", "rM.2") else "no pack") == "r26.11.9"
    assert catalog.latest_release(root, lambda r: None if r in ("rN.3", "rM.2") else "no pack") == "rN.3"
    with pytest.raises(catalog.NoLatestRelease) as e:
        catalog.latest_release(root, lambda r: f"{r} has no pack")
    for name in ("rN.3", "rM.2", "r26.11.9", "r26.11.10"):
        assert f"{name} has no pack" in str(e.value)
    with pytest.raises(catalog.NoLatestRelease, match="manifest"):
        catalog.latest_release(tmp_path / "empty")
    assert _tree(root) == before  # no pointer file, nothing else (AGENTS §3.10)


def test_w37_a_release_named_latest_cannot_be_defined(tmp_path, capsys):
    root = _root(tmp_path, [])
    with pytest.raises(ValueError, match="latest"):
        catalog.define_release(root, "latest", ["S0"])
    assert not catalog.release_def_path(root, "latest").exists()


BUILDS = [
    ["release", "define", "--stages", "S0"],
    ["release", "cluster"],
    ["release", "cluster", "plan", "--tmp", "{tmp}"],
    ["release", "tokens", "--without-clusters"],
    ["release", "postings"],
    ["release", "dim"],
    ["release", "manifest"],
    ["release", "catalog"],
    ["pack", "build"],
    ["build-stage", "--profile", "{tmp}/profile.yaml", "--stage", "S0", "--genomes", "{tmp}/g.json", "--objects",
     "{tmp}/o.json"],
]


@pytest.mark.parametrize("argv", BUILDS, ids=lambda a: " ".join(a[:3]))
def test_w37_a_build_command_refuses_latest(tmp_path, capsys, monkeypatch, argv):
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)
    monkeypatch.delenv("BIOBOX_STAGING", raising=False)
    root = _root(tmp_path, ["rN.3", "r26.11.9"])
    before = _tree(root)
    args = [x.format(tmp=tmp_path) for x in argv]
    code, out, err = _cli(capsys, *args, "--root", str(root), "--release", "latest")
    assert code == 2, (args, out, err)
    assert "--release latest" in err and "build command" in err, err
    assert "r26.11.9" not in out
    assert _tree(root) == before


READS = [
    ["find", "--accession", "NC_001604"],
    ["seq", "--protein", "00" * 16],
    ["export", "proteins"],
    ["sql", "-e", "SELECT 1 AS x"],
    ["nbr", "summarize", "00" * 16],
    ["pack", "verify"],
    ["release", "verify"],
]


@pytest.mark.parametrize("argv", READS, ids=lambda a: " ".join(a[:2]))
def test_w37_a_reading_command_with_no_release_that_has_a_manifest_exits_2_with_the_reason(tmp_path, capsys,
                                                                                           monkeypatch, argv):
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)
    root = _root(tmp_path, [], defined=["r26.11.9"])  # defined, but no manifest
    before = _tree(root)
    extra = ["--pack-dir", str(tmp_path / "packs")] if argv[0] == "pack" else []
    code, out, err = _cli(capsys, *argv, "--root", str(root), "--release", "latest", *extra)
    assert code == 2, (argv, out, err)
    # the refusal of the resolution itself (before W-37 some of these commands refused with another message that
    # happened to show `--release latest` in a suggested command: run 44104739971)
    assert "refused: --release latest:" in err and "manifest" in err, err
    assert out == ""
    assert _tree(root) == before


def test_w37_a_search_with_no_verified_pack_exits_2_naming_each_release_and_writes_nothing(tmp_path, capsys,
                                                                                          monkeypatch):
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)
    root = _root(tmp_path, ["rN.3", "r26.11.9"])
    packs = tmp_path / "packs"
    packs.mkdir()
    recipe = tmp_path / "r.yaml"
    recipe.write_text("recipe: w37\n")
    before = _tree(root)
    where = ["--root", str(root), "--release", "latest", "--pack-dir", str(packs), "--runs", str(tmp_path / "runs"),
             "--work", str(tmp_path / "work")]
    for argv in (["search", "run", str(recipe), *where, "--wait"],
                 ["search", "hmm", *where, "--hmm", str(tmp_path / "x.hmm"), "--bits", "25", "--positive",
                  "NC_001604:3170-5822", "--out", str(tmp_path / "out")]):
        code, out, err = _cli(capsys, *argv)
        assert code == 2, (argv, out, err)
        assert "--release latest" in err and "rN.3" in err and "r26.11.9" in err and "pack" in err, err
    assert _tree(root) == before
    for d in ("runs", "work", "out"):
        assert not (tmp_path / d).exists(), d


def test_w37_releases_marks_the_release_latest_names(tmp_path, capsys, monkeypatch):
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)
    root = _root(tmp_path, ["rN.3", "rM.2", "r26.11.9", "r26.11.10"], defined=["r26.12.1"])
    before = _tree(root)
    code, out, err = _cli(capsys, "releases", "--root", str(root), "--json")
    assert code == 0, err
    marks = {r["release"]: r["latest"] for r in map(json.loads, out.splitlines())}
    assert marks == {"rN.3": [], "rM.2": [], "r26.11.9": [], "r26.11.10": ["read"], "r26.12.1": []}
    code, out, err = _cli(capsys, "releases", "--root", str(root))
    assert code == 0, err
    lines = [x.split("\t") for x in out.splitlines()]
    shown = {r[0]: dict(zip(lines[0], r))["latest"] for r in lines[1:]}
    assert shown == {"rN.3": "-", "rM.2": "-", "r26.11.9": "-", "r26.11.10": "read", "r26.12.1": "-"}
    assert _tree(root) == before
