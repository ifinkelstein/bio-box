"""Card W-37 (DEC-52): `--release latest` names, at each call, the release with the greatest name that has a manifest
(and, for a command that searches, a verified search pack). Ring R1: roots made here with manifest files written by
hand, so no release is built.

The order is DEC-52's, written out by hand below as the oracle (card amended 2026-11-17): any name of neither scheme
(`rS<k>.<n>`, a test's own) first, by name; then the seven names of the old scheme in the order they were built
(DEC-52: `rN.1a`, `rN.1b`, `rN.1`, `rN.2`, `rN.3`, `rM.1`, `rM.2`); then names of the scheme `r<YY>.<MM>.<n>` by
year, month and number as numbers (`r26.11.10` after `r26.11.9`). A build command refuses `latest`, and a release named
`latest` cannot be defined. Resolving writes nothing (AGENTS §3.10: no pointer file)."""

from __future__ import annotations

import hashlib
import json
import random
from pathlib import Path

import pytest

from biobox import catalog, cli

pytestmark = pytest.mark.covers("W-37")

# DEC-52's order, by hand (card amended 2026-11-17): other names by name, the old scheme as built, then the new scheme
# by (year, month, number) as numbers
OTHER = ["r26.11", "r26.11.1a", "rN.4", "rS0.1", "rW37.a"]  # by name; `rN.4` is no release of the fixed list
BUILT = ["rN.1a", "rN.1b", "rN.1", "rN.2", "rN.3", "rM.1", "rM.2"]  # DEC-52 (1) and STATUS: the order of building
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
    names = OTHER + BUILT + NEW
    random.Random(37).shuffle(names)
    assert sorted(names, key=catalog.release_order) == OTHER + BUILT + NEW
    # numbers, not text: r26.11.10 after r26.11.9, r26.12.1 after both, r27.1.2 after r27.01.1; the old scheme as built
    for lo, hi in (("r26.11.9", "r26.11.10"), ("r26.11.10", "r26.12.1"), ("r27.01.1", "r27.1.2"), ("rN.3", "rM.1"),
                   ("rM.1", "rM.2"), ("rM.2", "r26.11.1"), ("rN.1a", "rN.1"), ("rS0.1", "rN.1a")):
        assert catalog.release_order(lo) < catalog.release_order(hi), (lo, hi)


def _pack(packs: Path, release: str) -> None:
    """A search pack of `release` with no shard, its directory named by its manifest's hash and with the build's record
    that it verified, as `pack.build` writes them: enough for the check a search makes of its pack."""
    from biobox.search import pack

    data = json.dumps({"release_id": release, "shards": [], "meta": []}, sort_keys=True).encode()
    pid = pack.pack_id_of(data)
    d = pack.pack_dir(packs, release, pid)
    d.mkdir(parents=True)
    (d / pack.MANIFEST).write_bytes(data)
    (d / pack.VERIFIED).write_bytes(pack._record_bytes(pid, release, data))


def test_w37_latest_of_the_old_releases_is_the_last_built_rM2_not_rN3(tmp_path, capsys, monkeypatch):
    """The pod's data root as of 2026-11-17: rN.3 and rM.2 both have a manifest and a verified pack, and the first
    release of the new scheme has a release_def only. `latest` is rM.2, built after rN.3, for reading and searching."""
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)
    root = _root(tmp_path, BUILT, defined=["r26.11.1"])
    packs = tmp_path / "packs"
    for r in ("rN.3", "rM.2"):
        _pack(packs, r)
    assert catalog.latest_release(root) == "rM.2"
    assert catalog.latest_release(root, lambda r: None if r in ("rN.3", "rM.2") else "no pack") == "rM.2"
    code, out, err = _cli(capsys, "releases", "--root", str(root), "--pack-dir", str(packs), "--json")
    assert code == 0, err
    marks = {r["release"]: r["latest"] for r in map(json.loads, out.splitlines())}
    assert marks == {**{r: [] for r in BUILT}, "rM.2": ["read", "search"], "r26.11.1": []}, marks


def test_w37_latest_is_the_greatest_name_with_a_manifest_and_resolving_writes_nothing(tmp_path):
    root = _root(tmp_path, ["rN.3", "rM.2", "r26.11.9", "r26.11.10"], defined=["r26.12.1"])
    before = _tree(root)
    assert catalog.latest_release(root) == "r26.11.10"  # r26.12.1 has a release_def only: not a candidate
    # a condition that only some releases meet: the greatest of those
    assert catalog.latest_release(root, lambda r: None if r in ("r26.11.9", "rM.2") else "no pack") == "r26.11.9"
    assert catalog.latest_release(root, lambda r: None if r in ("rN.3", "rM.2") else "no pack") == "rM.2"
    assert catalog.latest_release(root, lambda r: None if r == "rN.3" else "no pack") == "rN.3"
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
