"""Card W-29: `results archive RUN --under DIR` finds the output directories of a run by the `RUN` file each export
writes (key<TAB>value lines: `run`, `release`), refuses a directory whose `RUN` names another run, and gives the same
archive bytes for the directories it found as for the same directories named with `--with` (ring R1, constructed
runs; the real exports are in tests/r2/test_w29_run_file_search.py). The `RUN` text here is written by hand from the
card's definition, not by the function under test."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from test_results import _archive, _make_out, _make_run, _profile, _rows, _tree, _visible  # noqa: E402

pytestmark = pytest.mark.covers("W-29", "W-28")


def _run_file(d: Path, n: int, release: str = "rX.1") -> None:
    d.mkdir(parents=True, exist_ok=True)
    (d / "RUN").write_text(f"run\t{bytes([n]).hex() * 16}\nrelease\t{release}\n", encoding="utf-8")


def _spool_tar(base: Path, n: int) -> bytes:
    return (base / "tmp/results/spool" / f"{bytes([n]).hex() * 16}.tar").read_bytes()


def test_w29_under_finds_the_directories_that_name_the_run_and_gives_the_bytes_of_naming_them(tmp_path, capsys):
    run = _make_run(tmp_path / "runs", 1)
    a, b = _make_out(tmp_path / "work" / "hmm", 1), _make_out(tmp_path / "work" / "deeper" / "flanks", 5)
    _run_file(a, 1)
    _run_file(b, 1)
    (tmp_path / "work" / "unrelated").mkdir()  # no RUN: not found
    (tmp_path / "work" / "unrelated" / "x.tsv").write_text("x\n")
    found, named = tmp_path / "found", tmp_path / "named"
    code, out, err = _archive(capsys, _profile(found), run, "--under", str(tmp_path / "work"))
    assert code == 0, err
    code, out2, err = _archive(capsys, _profile(named), run, "--with", str(a), "--with", str(b))
    assert code == 0, err
    assert _spool_tar(found, 1) == _spool_tar(named, 1)
    assert _rows(out)[0]["outputs"] == "flanks,hmm" == _rows(out2)[0]["outputs"]


def test_w29_under_skips_a_directory_that_names_another_run_and_one_below_the_depth_bound(tmp_path, capsys):
    from biobox.search import results

    run = _make_run(tmp_path / "runs", 2)
    mine = _make_out(tmp_path / "w" / "mine", 2)
    _run_file(mine, 2)
    other = _make_out(tmp_path / "w" / "other", 7)
    _run_file(other, 7)  # another run's export beside it: left alone by --under
    deep = tmp_path / "w"
    for k in range(results.UNDER_DEPTH + 1):
        deep = deep / f"d{k}"
    _make_out(deep / "too_deep", 2)
    _run_file(deep / "too_deep", 2)
    code, out, err = _archive(capsys, _profile(tmp_path / "p"), run, "--under", str(tmp_path / "w"))
    assert code == 0, err
    assert _rows(out)[0]["outputs"] == "mine"
    # a directory exactly at the bound is found
    ok = tmp_path / "v"
    for k in range(results.UNDER_DEPTH):
        ok = ok / f"d{k}"
    _make_out(ok, 2)
    _run_file(ok, 2)
    code, out, err = _archive(capsys, _profile(tmp_path / "p2"), run, "--under", str(tmp_path / "v"))
    assert code == 0 and _rows(out)[0]["outputs"] == f"d{results.UNDER_DEPTH - 1}", err


def test_w29_a_with_directory_whose_run_file_names_another_run_is_refused_naming_both(tmp_path, capsys):
    prof = _profile(tmp_path / "p")
    run = _make_run(tmp_path / "runs", 3)
    other = _make_out(tmp_path / "other", 3)
    _run_file(other, 9)
    code, out, err = _archive(capsys, prof, run, "--with", str(other))
    assert code == 2 and "refused" in err, err
    assert ("09" * 16) in err and ("03" * 16) in err, err
    assert not out and _visible(tmp_path / "p" / "tmp/results/spool") == []
    # the same run id but another release: refused too
    same = _make_out(tmp_path / "same", 3)
    _run_file(same, 3, release="rY.9")
    code, out, err = _archive(capsys, prof, run, "--with", str(same))
    assert code == 2 and "rY.9" in err and "rX.1" in err, err
    # a RUN file that is not two key<TAB>value lines
    bad = _make_out(tmp_path / "bad", 3)
    (bad / "RUN").write_text("this is not a run file\n")
    code, out, err = _archive(capsys, prof, run, "--with", str(bad))
    assert code == 2 and "RUN" in err, err


def test_w29_a_directory_without_a_run_file_is_accepted_through_with_as_before(tmp_path, capsys):
    prof = _profile(tmp_path / "p")
    run = _make_run(tmp_path / "runs", 4)
    out = _make_out(tmp_path / "plain", 4)
    code, stdout, err = _archive(capsys, prof, run, "--with", str(out))
    assert code == 0, err
    # and the old recipe check still refuses another run's `search hmm` directory
    run5 = _make_run(tmp_path / "runs", 5, recipe_yaml="recipe: another\n")
    code, _, err = _archive(capsys, prof, run5, "--with", str(out))
    assert code == 2 and "not the recipe run" in err, err
    # --under never takes a directory without a RUN file
    code, stdout, err = _archive(capsys, _profile(tmp_path / "p3"), run, "--under", str(tmp_path / "plain"))
    assert code == 0 and _rows(stdout)[0]["outputs"] == "-", stdout


def test_w29_a_found_directory_inside_a_found_directory_is_kept_once_in_the_outer(tmp_path, capsys):
    run = _make_run(tmp_path / "runs", 6)
    outer = _make_out(tmp_path / "w" / "hmm", 6)
    _run_file(outer, 6)
    inner = _make_out(outer / "flanks", 8)
    _run_file(inner, 6)
    found, named = tmp_path / "found", tmp_path / "named"
    code, _, err = _archive(capsys, _profile(found), run, "--under", str(tmp_path / "w"))
    assert code == 0, err
    code, _, err = _archive(capsys, _profile(named), run, "--with", str(outer))
    assert code == 0, err
    assert _spool_tar(found, 6) == _spool_tar(named, 6)


def test_w29_under_may_be_repeated_and_a_directory_given_both_ways_is_kept_once(tmp_path, capsys):
    run = _make_run(tmp_path / "runs", 10)
    a = _make_out(tmp_path / "x" / "a", 10)
    b = _make_out(tmp_path / "y" / "b", 11)
    _run_file(a, 10)
    _run_file(b, 10)
    both, named = tmp_path / "both", tmp_path / "named"
    code, _, err = _archive(capsys, _profile(both), run, "--under", str(tmp_path / "x"), "--under",
                            str(tmp_path / "y"), "--with", str(a))
    assert code == 0, err
    code, _, err = _archive(capsys, _profile(named), run, "--with", str(a), "--with", str(b))
    assert code == 0, err
    assert _spool_tar(both, 10) == _spool_tar(named, 10)
    assert _tree(a)["RUN"]  # the RUN file is archived with its directory
