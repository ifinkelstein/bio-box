"""Ring R2: card W-28 (DEC-33) on the fixture: a real `search hmm` run on micro-S0 archived with its outputs,
bundled into a durable area of the test's own (never /stor) at a lowered floor, listed and restored equal."""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

from biobox import cli
from biobox.search import runner

sys.path.insert(0, str(Path(__file__).parent))
import test_micro_s0_search as s0  # noqa: E402
from test_micro_s0_search import world  # noqa: E402,F401  (the module fixture: micro-S0, its pack)

pytestmark = [pytest.mark.r2, pytest.mark.covers("W-28", "W-07")]


def _profile(base: Path, floor: int) -> Path:
    (base / "durable").mkdir(parents=True)
    p = base / "profile.yaml"
    p.write_text(yaml.safe_dump({
        "max_threads_day": 2, "max_threads_night": 2, "mem_scope": "8G", "tmp": str(base / "tmp"),
        "durable": str(base / "durable"), "bulk": str(base / "bulk"), "nice_monolithic": 19, "nice_fanout": 15,
        "day_start_hour": 8, "day_end_hour": 20, "cpu_min_per_process": 300, "small_file_floor_bytes": floor}))
    return p


def _tree(d: Path) -> dict[str, bytes]:
    return {p.relative_to(d).as_posix(): p.read_bytes() for p in sorted(Path(d).rglob("*")) if p.is_file()}


def test_w28_a_fixture_search_is_archived_with_its_outputs_bundled_listed_and_restored_equal(world, tmp_path, capsys,
                                                                                              monkeypatch):
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)
    hmm, pos = s0._w07_hmm(world, tmp_path)
    out, runs = tmp_path / "out", tmp_path / "runs"
    argv = ["search", "hmm", "--root", str(world[0]), "--release", s0.REL, "--hmm", str(hmm), "--bits", "25",
            "--positive", s0._w07_spec(world[0], pos), "--out", str(out), "--cpus", "1",
            "--pack-dir", str(s0._packs(world)), "--runs", str(runs), "--work", str(tmp_path / "work")]
    capsys.readouterr()
    code = cli.main(argv)
    got = capsys.readouterr()
    assert code == 0, got.err
    rid = got.out.splitlines()[0].split("\t")[0]
    rdir = runner.run_dir(runs, rid)
    hits = tmp_path / "extra" / "hits_again.tsv"  # an output written apart from the run, given with --with
    assert cli.main(["search", "hits", rid, "--root", str(world[0]), "--runs", str(runs), "--pack-dir",
                     str(s0._packs(world)), "--format", "tsv", "--out", str(hits)]) == 0
    # two spools on two profiles: the same run gives the same archive
    spooled = []
    for name in ("one", "two"):
        prof = _profile(tmp_path / name, floor=1000)  # lowered: the fixture run is far under 64 MiB
        capsys.readouterr()
        assert cli.main(["results", "archive", rid, "--runs", str(runs), "--profile", str(prof),
                         "--with", str(out), "--with", str(hits)]) == 0
        got = capsys.readouterr()
        line = dict(zip(got.out.splitlines()[0].split("\t"), got.out.splitlines()[1].split("\t")))
        assert line["outputs"] == "hits_again.tsv,out" and line["release"] == s0.REL and line["status"] == "ok"
        spooled.append((tmp_path / name / "tmp/results/spool" / f"{rid}.tar").read_bytes())
    assert spooled[0] == spooled[1]
    # a `search hmm` directory holds its run's recipe: given with another run, it is refused
    other = tmp_path / "other_run"
    shutil.copytree(rdir, other)
    pq.write_table(pq.read_table(other / "run.parquet").set_column(
        2, "recipe_yaml", pa.array(["recipe: another\n"], pa.string())), other / "run.parquet")
    capsys.readouterr()
    assert cli.main(["results", "archive", str(other), "--profile", str(tmp_path / "two" / "profile.yaml"),
                     "--with", str(out)]) == 2
    assert "not the recipe run" in capsys.readouterr().err
    prof = tmp_path / "one" / "profile.yaml"
    assert cli.main(["results", "flush", "--profile", str(prof)]) == 0
    bundle = Path(capsys.readouterr().out.split("\t")[0])
    assert bundle.parent == tmp_path / "one" / "durable" / "results" and bundle.is_file()
    assert sorted(n for n in os.listdir(bundle.parent) if not n.startswith(".")) == [bundle.name]
    assert cli.main(["results", "list", "--profile", str(prof)]) == 0
    rows = capsys.readouterr().out.splitlines()
    assert len(rows) == 2 and rows[1].startswith(rid) and rows[1].endswith(str(bundle))
    back = tmp_path / "back"
    assert cli.main(["results", "get", rid, "--profile", str(prof), "--out", str(back)]) == 0
    restored = _tree(back)
    assert restored.pop("FILES.tsv")
    want = {f"run/{k}": v for k, v in _tree(rdir).items()}
    want.update({f"outputs/out/{k}": v for k, v in _tree(out).items()})
    want["outputs/hits_again.tsv"] = hits.read_bytes()
    assert restored == want
    assert {"outputs/out/hits.tsv", "outputs/out/hits.faa", "outputs/out/neighbors.tsv", "run/run.parquet",
            "run/hits.parquet"} <= set(restored)
