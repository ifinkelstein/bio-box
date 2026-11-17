"""Card W-37 (DEC-52): `--release latest` on built releases (ring R2).

The micro-S0 search world (`test_micro_s0_search.world`): release REL (`rS0.l3`, a name outside the scheme
`r<YY>.<MM>.<n>`) with its manifest, catalog and a verified search pack. A copy then gains releases over the same
stage named by the scheme: `r26.11.9` with a manifest and no pack, and `r26.11.10` with only a `release_def`. So
`latest` is `r26.11.9` for a command that reads a release and REL for one that searches; with no new release it is
REL for both.

The oracle of each output is the same command with the release's own name, run here; the names that `latest` must
give are written by hand from DEC-52's rule. Every output and run record must name the release, never `latest`.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

import pyarrow.parquet as pq
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from test_micro_s0_search import (REL, _copy, _packs_of, _recipe, _seed, _w06_anchors, _w07_hmm,  # noqa: E402
                                  _w07_spec, world)  # noqa: F401  (the module-scoped search world)
from test_w23_commands import _s0_genome_input  # noqa: E402

from biobox import catalog, cli, dim  # noqa: E402
from biobox.search import pack, runner  # noqa: E402

pytestmark = [pytest.mark.r2, pytest.mark.covers("W-37")]

READ_NEW = "r26.11.9"  # a manifest, no pack
DEF_ONLY = "r26.11.10"  # a release_def, no manifest


@pytest.fixture(scope="module")
def base(world, tmp_path_factory):
    """A copy of the world with S0's `genome_input`, REL's `contig_dim`, manifest and catalog (as card W-23's `rel`),
    and the world's pack beside it."""
    d = tmp_path_factory.mktemp("w37")
    root = _copy(world, d)
    _s0_genome_input(root)
    (d / "dimtmp").mkdir()
    dim.run(root, REL, ["S0"], tmp=d / "dimtmp")
    man = catalog.build_manifest(root, REL, ["S0"])
    catalog.write_manifest(man)
    catalog.generate_catalog(man, root=root)
    return root


def _clone(base: Path, tmp_path: Path, *, new: bool) -> Path:
    """A copy of `base` and its pack directory; with `new`, also READ_NEW (manifest) and DEF_ONLY (release_def)."""
    root = tmp_path / "a"
    shutil.copytree(base, root)
    shutil.copytree(_packs_of(base), _packs_of(root))
    catalog.define_release(root, DEF_ONLY, ["S0"])
    if new:
        catalog.define_release(root, READ_NEW, ["S0"])
        catalog.write_manifest(catalog.build_manifest(root, READ_NEW, ["S0"]))
    return root


def _tree(root: Path) -> dict[str, str]:
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else "dir"
            for p in sorted(root.rglob("*"))}


def _cli(capsys, argv) -> tuple[int, str, str]:
    capsys.readouterr()
    code = cli.main([str(x) for x in argv])
    got = capsys.readouterr()
    return code, got.out, got.err


def _no_latest_in(d: Path) -> None:
    """No file under `d` says `latest`."""
    for p in sorted(d.rglob("*")):
        if p.is_file():
            assert b"latest" not in p.read_bytes(), p


def _lambda_accession(root: Path) -> str:
    return next(r["accession"] for f in sorted((root / "l0seq/contig_label").rglob("*.parquet"))
                for r in pq.read_table(f, partitioning=None).to_pylist() if r["accession"].startswith("NC_001416"))


def test_w37_without_a_new_release_latest_is_the_old_release_for_every_reading_command(base, tmp_path, capsys,
                                                                                      monkeypatch):
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)
    root = _clone(base, tmp_path, new=False)  # DEF_ONLY has no manifest: not a candidate
    before = _tree(root)
    prot = sorted(_w06_anchors(root))[0].hex()
    uids = ",".join(sorted(u.hex() for u in _w06_anchors(root)))
    r = ["--root", root]
    commands = {
        "sql": ["sql", *r, "-e", "SELECT value FROM catalog_info WHERE key = 'release_id'"],
        "seq": ["seq", *r, "--protein", prot],
        "find": ["find", *r, "--accession", _lambda_accession(root)],
        "export": ["export", "proteins", *r],
        "nbr": ["nbr", "summarize", uids, *r, "--window", "3g", "--json"],
        "pack verify": ["pack", "verify", *r, "--pack-dir", _packs_of(root)],
        "release verify": ["release", "verify", *r],
    }
    for what, argv in commands.items():
        code, out, err = _cli(capsys, [*argv, "--release", "latest"])
        assert code == 0, (what, err)
        assert f"--release latest: {REL}" in err, (what, err)
        code2, named, _ = _cli(capsys, [*argv, "--release", REL])
        assert code2 == 0 and out == named and out, what
        assert "latest" not in out, what
    code, out, _ = _cli(capsys, commands["sql"] + ["--release", "latest"])
    assert out.splitlines() == ["value", REL]
    assert _tree(root) == before


def test_w37_reading_takes_the_newest_manifest_and_searching_the_newest_verified_pack(base, world, tmp_path, capsys,
                                                                                     monkeypatch):
    """The test's name keeps `latest` out of its tmp_path (pytest names it after the test), which the run record and
    recipe.yaml hold as paths: so a `latest` found below is one the command wrote."""
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)
    root = _clone(base, tmp_path, new=True)
    packs = _packs_of(root)
    before = _tree(root)
    # reading: READ_NEW, the greatest name with a manifest; its own name in the output
    prot = sorted(_w06_anchors(root))[0].hex()
    code, out, err = _cli(capsys, ["seq", "--root", root, "--release", "latest", "--protein", prot])
    assert code == 0 and f"--release latest: {READ_NEW}" in err, err
    assert out.startswith(f">{prot} release={READ_NEW} ")
    assert _cli(capsys, ["seq", "--root", root, "--release", READ_NEW, "--protein", prot])[1] == out
    f_latest, f_named = tmp_path / "latest.faa", tmp_path / "named.faa"
    assert _cli(capsys, ["export", "proteins", "--root", root, "--release", "latest", "--out", f_latest])[0] == 0
    assert _cli(capsys, ["export", "proteins", "--root", root, "--release", READ_NEW, "--out", f_named])[0] == 0
    assert f_latest.read_bytes() == f_named.read_bytes()
    # `releases` marks both
    code, out, err = _cli(capsys, ["releases", "--root", root, "--pack-dir", packs, "--json"])
    assert code == 0, err
    marks = {r["release"]: r["latest"] for r in map(json.loads, out.splitlines())}
    assert marks == {REL: ["search"], READ_NEW: ["read"], DEF_ONLY: []}, marks
    code, out, _ = _cli(capsys, ["releases", "--root", root, "--pack-dir", packs])
    lines = [x.split("\t") for x in out.splitlines()]
    assert {x[0]: dict(zip(lines[0], x))["latest"] for x in lines[1:]} == {REL: "search", READ_NEW: "read",
                                                                           DEF_ONLY: "-"}
    # searching: REL, the only release with a verified pack; the run and DIR name it
    hmm, pos = _w07_hmm(world, tmp_path)
    out_dir, runs, work = tmp_path / "out", tmp_path / "runs", tmp_path / "work"
    code, out, err = _cli(capsys, ["search", "hmm", "--root", root, "--release", "latest", "--hmm", hmm, "--bits", "25",
                                   "--neighbors", "3", "--positive", _w07_spec(root, pos), "--out", out_dir,
                                   "--cpus", "1", "--pack-dir", packs, "--runs", runs, "--work", work])
    assert code == 0, err
    assert f"--release latest: {REL}" in err, err
    rid = out.splitlines()[0].split("\t")[0]
    assert (out_dir / "RUN").read_text() == f"run\t{rid}\nrelease\t{REL}\n"
    assert f"written_against: {REL}\n" in (out_dir / "recipe.yaml").read_text()
    rec = runner.read_run(runner.run_dir(runs, rid))
    assert rec["ran_against"] == REL
    command = json.loads(rec["limits_json"])["command"]
    assert command[command.index("--release") + 1] == REL and "latest" not in command, command
    _no_latest_in(out_dir)
    _no_latest_in(runner.run_dir(runs, rid))
    # `search run --wait`, the other way a search starts
    hit, _ = _seed(root)
    recipe = tmp_path / "recipe.yaml"
    recipe.write_text(_recipe(hit, hit))
    one = tmp_path / "one"
    code, out, err = _cli(capsys, ["search", "run", recipe, "--root", root, "--release=latest", "--wait", "--cpus", "1",
                                   "--pack-dir", packs, "--work", work, "--out", one])
    assert code == 0, err
    assert f"--release latest: {REL}" in err, err
    rec = runner.read_run(one)
    assert rec["ran_against"] == REL
    command = json.loads(rec["limits_json"])["command"]
    assert f"--release={REL}" in command and not any("latest" in x for x in command), command
    _no_latest_in(one)
    assert _tree(root) == before  # nothing written under the root: no pointer file


def test_w37_a_search_takes_only_a_pack_whose_build_recorded_that_it_verified(base, world, tmp_path, capsys,
                                                                            monkeypatch):
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)
    root = _clone(base, tmp_path, new=True)
    packs = _packs_of(root)
    (pack.find(packs, REL) / pack.VERIFIED).unlink()  # a pack built before card W-19: no record
    hmm, pos = _w07_hmm(world, tmp_path)
    before = _tree(root)
    argv = ["search", "hmm", "--root", root, "--release", "latest", "--hmm", hmm, "--bits", "25",
            "--positive", _w07_spec(root, pos), "--out", tmp_path / "out", "--cpus", "1", "--pack-dir", packs,
            "--runs", tmp_path / "runs", "--work", tmp_path / "work"]
    code, out, err = _cli(capsys, argv)
    assert code == 2 and out == "", err
    assert "--release latest" in err and pack.VERIFIED in err and REL in err and READ_NEW in err, err
    code, out, _ = _cli(capsys, ["releases", "--root", root, "--pack-dir", packs, "--json"])
    assert {r["release"]: r["latest"] for r in map(json.loads, out.splitlines())} == \
        {REL: [], READ_NEW: ["read"], DEF_ONLY: []}
    shutil.rmtree(packs)  # no pack directory at all
    packs.mkdir()
    code, out, err = _cli(capsys, argv)
    assert code == 2 and "--release latest" in err and "no search pack" in err, err
    for d in ("out", "runs", "work"):
        assert not (tmp_path / d).exists(), d
    assert _tree(root) == before
