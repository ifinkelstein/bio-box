"""Ring R2: card W-29 on the micro-S0 fixture. A run's record keeps the user's command line (in `limits_json`), and
every directory export written from a run (`search hmm`, `nbr --format tsv --out DIR`, `nbr flanks --out DIR`) holds a
`RUN` file: the run id and the release, one `key<TAB>value` per line, nothing else. `results archive RUN --under DIR`
finds those directories and gives the bytes of naming them with `--with`."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest
import yaml

from biobox import cli
from biobox.search import runner

sys.path.insert(0, str(Path(__file__).parent))
import test_micro_s0_search as s0  # noqa: E402
from test_micro_s0_search import world  # noqa: E402,F401  (the module fixture: micro-S0, its pack)

pytestmark = [pytest.mark.r2, pytest.mark.covers("W-29", "W-28", "W-07")]


def _profile(base: Path) -> Path:
    (base / "durable").mkdir(parents=True)
    p = base / "profile.yaml"
    p.write_text(yaml.safe_dump({
        "max_threads_day": 2, "max_threads_night": 2, "mem_scope": "8G", "tmp": str(base / "tmp"),
        "durable": str(base / "durable"), "bulk": str(base / "bulk"), "nice_monolithic": 19, "nice_fanout": 15,
        "day_start_hour": 8, "day_end_hour": 20, "cpu_min_per_process": 300, "small_file_floor_bytes": 1000}))
    return p


def _hmm_argv(world, tmp_path, out: Path, runs: Path) -> list[str]:
    hmm, pos = s0._w07_hmm(world, tmp_path)
    return ["search", "hmm", "--root", str(world[0]), "--release", s0.REL, "--hmm", str(hmm), "--bits", "25",
            "--positive", s0._w07_spec(world[0], pos), "--out", str(out), "--cpus", "1",
            "--pack-dir", str(s0._packs(world)), "--runs", str(runs), "--work", str(tmp_path / "work")]


def test_w29_search_hmm_records_the_users_command_and_writes_a_run_file(world, tmp_path, capsys):
    out, runs = tmp_path / "out", tmp_path / "runs"
    argv = _hmm_argv(world, tmp_path, out, runs)
    capsys.readouterr()
    code = cli.main(argv)
    got = capsys.readouterr()
    assert code == 0, got.err
    rid = got.out.splitlines()[0].split("\t")[0]
    # the file: the run id and the release, one key<TAB>value per line, nothing else (no path, host or date)
    assert (out / "RUN").read_bytes() == f"run\t{rid}\nrelease\t{s0.REL}\n".encode()
    assert not (out / "loci" / "RUN").exists()  # the tables of loci/ are as `nbr` alone writes them, one per locus
    # the record: the user's command line, not the worker's (`argv` is the worker's)
    rec = runner.read_run(runner.run_dir(runs, rid))
    lim = json.loads(rec["limits_json"])
    assert lim["command"][0] == "bio-box" and lim["command"][1:] == argv
    assert rec["argv"] and rec["argv"] != lim["command"]  # the existing field is as it was


def test_w29_nbr_tables_and_flanks_of_a_run_carry_its_run_file_and_a_uid_list_carries_none(world, tmp_path, capsys):
    out, runs = tmp_path / "out", tmp_path / "runs"
    capsys.readouterr()
    assert cli.main(_hmm_argv(world, tmp_path, out, runs)) == 0, capsys.readouterr().err
    rid = (out / "RUN").read_text().split("\t")[1].split("\n")[0]
    where = ["--root", str(world[0]), "--pack-dir", str(s0._packs(world)), "--runs", str(runs)]
    want = f"run\t{rid}\nrelease\t{s0.REL}\n"
    loci = tmp_path / "loci"
    assert cli.main(["nbr", rid, *where, "--window", "5g", "--format", "tsv", "--out", str(loci)]) == 0
    assert (loci / "RUN").read_text() == want
    flk = tmp_path / "flanks"
    capsys.readouterr()
    code = cli.main(["nbr", "flanks", rid, *where, "--window", "5g", "--out", str(flk)])
    assert code == 0, capsys.readouterr().err
    assert (flk / "RUN").read_text() == want and (flk / "README.txt").is_file()
    # the exports of one run are found by --under, and archive to the bytes of --with
    found, named = tmp_path / "found", tmp_path / "named"
    prof_f, prof_n = _profile(found), _profile(named)
    keep = tmp_path / "keep"
    keep.mkdir()
    (keep / "x.txt").write_text("not an export\n")
    capsys.readouterr()
    assert cli.main(["results", "archive", rid, "--runs", str(runs), "--profile", str(prof_f), "--under",
                     str(tmp_path)]) == 0, capsys.readouterr().err
    line = capsys.readouterr().out.splitlines()
    assert line[1].split("\t")[7] == "flanks,loci,out"
    assert cli.main(["results", "archive", rid, "--runs", str(runs), "--profile", str(prof_n), "--with", str(out),
                     "--with", str(loci), "--with", str(flk)]) == 0
    spool = Path("tmp/results/spool") / f"{rid}.tar"
    assert (found / spool).read_bytes() == (named / spool).read_bytes()
    # protein uids given directly name no run: no RUN file
    hits = runner.read_hits(runner.run_dir(runs, rid))
    anchors = ",".join(sorted({h["target_uid"].hex() for h in hits if h["target_kind"] == "protein"}))
    bare = tmp_path / "bare"
    capsys.readouterr()
    code = cli.main(["nbr", "--anchors", anchors, "--release", s0.REL, "--root", str(world[0]), "--pack-dir",
                     str(s0._packs(world)), "--window", "5g", "--format", "tsv", "--out", str(bare)])
    assert code == 0, capsys.readouterr().err
    assert not (bare / "RUN").exists()


def test_w29_a_run_file_of_another_run_in_an_nbr_output_directory_is_refused(world, tmp_path, capsys):
    out, runs = tmp_path / "out", tmp_path / "runs"
    capsys.readouterr()
    assert cli.main(_hmm_argv(world, tmp_path, out, runs)) == 0, capsys.readouterr().err
    rid = (out / "RUN").read_text().split("\t")[1].split("\n")[0]
    clash = tmp_path / "clash"
    clash.mkdir()
    (clash / "RUN").write_text(f"run\t{'ab' * 16}\nrelease\t{s0.REL}\n")
    where = ["--root", str(world[0]), "--pack-dir", str(s0._packs(world)), "--runs", str(runs)]
    capsys.readouterr()
    assert cli.main(["nbr", rid, *where, "--window", "5g", "--format", "tsv", "--out", str(clash)]) == 2
    err = capsys.readouterr().err
    assert "ab" * 16 in err and rid in err, err
    assert not list(clash.glob("locus_*"))  # nothing written


def test_w29_a_submitted_run_keeps_the_command_the_user_typed(world, tmp_path, capsys):
    root = world[0]
    runs, work = tmp_path / "runs", tmp_path / "work"
    hmm, pos = s0._w07_hmm(world, tmp_path)
    rec = tmp_path / "recipe.yaml"
    rec.write_text(s0._recipe(pos, pos))
    argv = ["search", "run", str(rec), "--root", str(root), "--release", s0.REL, "--cpus", "1", "--pack-dir",
            str(s0._packs(world)), "--runs", str(runs), "--work", str(work)]
    capsys.readouterr()
    code = cli.main(argv)
    got = capsys.readouterr()
    assert code == 0, got.err
    rid = got.out.strip()
    import time

    d = runner.run_dir(runs, rid)
    deadline = time.monotonic() + 120
    while not (d / "run.parquet").exists() and time.monotonic() < deadline:
        time.sleep(0.5)
    row = runner.read_run(d)
    lim = json.loads(row["limits_json"])
    assert lim["command"] == ["bio-box", *argv]  # not the detached runner's `--wait --submitted --run-id ...`
