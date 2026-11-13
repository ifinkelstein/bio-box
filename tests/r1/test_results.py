"""Card W-28 (DEC-33): a run's results archived, bundled onto durable storage, listed and restored (ring R1 by
directory). Runs here are constructed run directories (a `search_run` record and files of seeded random bytes);
the durable storage is a directory of the test's own (never /stor); the floor is lowered by the test profile's
`small_file_floor_bytes`. The fixture run of a real search is in tests/r2/test_results_fixture.py."""

from __future__ import annotations

import hashlib
import io
import json
import os
import random
import shutil
import subprocess
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

from biobox import cli
from biobox.registry.schemas import table

pytestmark = pytest.mark.covers("W-28")

WHEN = datetime(2026, 11, 3, 9, 30, 5, tzinfo=timezone.utc)


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)


def _profile(base: Path, floor: int = 10**6) -> Path:
    """A host profile whose tmp, durable and bulk are directories of the test; durable exists."""
    (base / "durable").mkdir(parents=True, exist_ok=True)
    data = {"max_threads_day": 2, "max_threads_night": 2, "mem_scope": "8G", "tmp": str(base / "tmp"),
            "durable": str(base / "durable"), "bulk": str(base / "bulk"), "nice_monolithic": 19, "nice_fanout": 15,
            "day_start_hour": 8, "day_end_hour": 20, "cpu_min_per_process": 300, "small_file_floor_bytes": floor}
    p = base / "profile.yaml"
    p.write_text(yaml.safe_dump(data))
    return p


def _record(rid: bytes, **over) -> pa.Table:
    """One search_run row: given fields, null where nullable, an empty value of the type elsewhere."""
    schema = table("search_run").schema
    given = {"run_id": rid, "ran_against": "rX.1", "written_against": "rX.1", "status": "ok", "submitted_at": WHEN,
             "host": "somehost", "argv": ["cli.py", "search", "run", "recipe.yaml"], "recipe_yaml": "recipe: x\n",
             **over}
    row = {}
    for f in schema:
        if f.name in given:
            row[f.name] = given[f.name]
        elif f.nullable:
            row[f.name] = None
        elif pa.types.is_fixed_size_binary(f.type):
            row[f.name] = b"\0" * f.type.byte_width
        elif pa.types.is_list(f.type):
            row[f.name] = []
        elif pa.types.is_string(f.type):
            row[f.name] = ""
        elif pa.types.is_timestamp(f.type):
            row[f.name] = WHEN
        else:
            row[f.name] = 0
    return pa.Table.from_pylist([row], schema=schema)


def _bytes(seed: int, n: int) -> bytes:
    return random.Random(seed).randbytes(n)


def _make_run(runs: Path, n: int, size: int = 50_000, reverse: bool = False, **over) -> Path:
    """Run directory `run=<id>` in `runs`: its record and three tables of seeded random bytes. `reverse` writes the
    files in the other order."""
    rid = bytes([n]) * 16
    d = runs / f"run={rid.hex()}"
    d.mkdir(parents=True)
    files = [("hits.parquet", _bytes(n, size)), ("profiles.parquet", _bytes(n + 100, 1000)),
             ("neighborhoods.parquet", _bytes(n + 200, 3 * size))]
    for name, data in (reversed(files) if reverse else files):
        (d / name).write_bytes(data)
    pq.write_table(_record(rid, **over), d / "run.parquet")
    return d


def _make_out(base: Path, n: int, reverse: bool = False) -> Path:
    """An output directory as `search hmm` writes one: a recipe, tables and a loci/ directory."""
    out = base
    files = [("recipe.yaml", b"recipe: x\n"), ("hits.tsv", _bytes(n + 1, 2000)), ("hits.faa", _bytes(n + 2, 3000)),
             ("loci/locus_A_1.tsv", _bytes(n + 3, 500)), ("loci/locus_B_2.tsv", _bytes(n + 4, 600))]
    for name, data in (reversed(files) if reverse else files):
        (out / name).parent.mkdir(parents=True, exist_ok=True)
        (out / name).write_bytes(data)
    return out


def _tree(d: Path) -> dict[str, bytes]:
    return {p.relative_to(d).as_posix(): p.read_bytes() for p in sorted(Path(d).rglob("*")) if p.is_file()}


def _main(capsys, *argv) -> tuple[int, str, str]:
    capsys.readouterr()
    code = cli.main(list(argv))
    got = capsys.readouterr()
    return code, got.out, got.err


def _archive(capsys, prof: Path, run: Path, *extra) -> tuple[int, str, str]:
    return _main(capsys, "results", "archive", str(run), "--profile", str(prof), *extra)


def _rows(out: str) -> list[dict]:
    lines = out.splitlines()
    cols = lines[0].split("\t")
    return [dict(zip(cols, x.split("\t"))) for x in lines[1:]]


def _visible(d: Path) -> list[str]:
    return sorted(n for n in os.listdir(d) if not n.startswith(".")) if d.is_dir() else []


# --- determinism ------------------------------------------------------------------------------------------------


def test_w28_the_same_run_gives_the_same_archive_bytes_wherever_and_however_its_files_were_written(tmp_path,
                                                                                                    capsys):
    from biobox.search import results

    a = _make_run(tmp_path / "a" / "runs", 1)
    out_a = _make_out(tmp_path / "a" / "rnap", 1)
    b = _make_run(tmp_path / "elsewhere" / "deeper" / "runs", 1, reverse=True)
    out_b = _make_out(tmp_path / "elsewhere" / "rnap", 1, reverse=True)
    for p in [*b.iterdir(), *out_b.rglob("*")]:  # other times and modes
        os.utime(p, (1_600_000_000, 1_600_000_000))
        if p.is_file():
            p.chmod(0o600)
    la = results.write_archive(a, [out_a], tmp_path / "a.tar")
    lb = results.write_archive(b, [out_b], tmp_path / "b.tar")
    data = (tmp_path / "a.tar").read_bytes()
    assert data == (tmp_path / "b.tar").read_bytes() and la == lb
    assert la["sha256"] == hashlib.sha256(data).hexdigest() and la["size"] == len(data)
    assert la["run_date"] == "2026-11-03T09:30:05Z"  # the run's own date, not the time of archiving
    assert (la["run_id"], la["release"], la["status"], la["outputs"]) == ("01" * 16, "rX.1", "ok", "rnap")
    assert str(tmp_path).encode() not in data  # no host path in the archive
    with tarfile.open(tmp_path / "a.tar") as tf:
        members = tf.getmembers()
    names = [m.name for m in members]
    assert names[0] == "FILES.tsv" and names[1:] == sorted(names[1:])
    assert set(names[1:]) == {"run/run.parquet", "run/hits.parquet", "run/profiles.parquet",
                              "run/neighborhoods.parquet", "outputs/rnap/recipe.yaml", "outputs/rnap/hits.tsv",
                              "outputs/rnap/hits.faa", "outputs/rnap/loci/locus_A_1.tsv",
                              "outputs/rnap/loci/locus_B_2.tsv"}
    assert la["n_files"] == 9
    for m in members:
        assert m.isfile() and (m.mtime, m.uid, m.gid, m.uname, m.gname, m.mode) == (0, 0, 0, "", "", 0o644), m.name


def test_w28_two_processes_with_other_hash_seeds_spool_the_same_bytes(tmp_path, capsys):
    run = _make_run(tmp_path / "runs", 2)
    out = _make_out(tmp_path / "rnap", 2)
    got = []
    for seed in ("0", "4242"):
        prof = _profile(tmp_path / f"seed{seed}")
        env = {**os.environ, "PYTHONHASHSEED": seed}
        env.pop("BIOBOX_PROFILE", None)
        p = subprocess.run([sys.executable, "-m", "biobox.cli", "results", "archive", str(run), "--profile",
                            str(prof), "--with", str(out)], env=env, capture_output=True, text=True, timeout=120)
        assert p.returncode == 0, p.stderr
        got.append((tmp_path / f"seed{seed}" / "tmp/results/spool" / f"{'02' * 16}.tar").read_bytes())
    assert got[0] == got[1]


# --- archive: no-op, refusal, crash -----------------------------------------------------------------------------


def test_w28_archiving_a_run_twice_is_a_no_op_and_a_changed_run_is_refused_naming_both_sha256(tmp_path, capsys):
    prof = _profile(tmp_path)
    run = _make_run(tmp_path / "runs", 3)
    out = _make_out(tmp_path / "rnap", 3)
    code, stdout, err = _archive(capsys, prof, run, "--with", str(out))
    assert code == 0, err
    first = _rows(stdout)[0]
    spool = tmp_path / "tmp/results/spool"
    assert _visible(spool) == [f"{'03' * 16}.tar", f"{'03' * 16}.tsv"]
    before = _tree(spool)
    code, stdout, err = _archive(capsys, prof, run, "--with", str(out))
    assert code == 0 and "already kept" in err and "nothing written" in err, err
    assert _rows(stdout)[0]["sha256"] == first["sha256"] and _tree(spool) == before
    (out / "hits.tsv").write_bytes(b"changed\n")
    code, stdout, err = _archive(capsys, prof, run, "--with", str(out))
    assert code == 2 and "refused" in err and first["sha256"] in err and "never overwritten" in err, err
    other = [w for w in err.split() if len(w.strip(";.,()")) == 64 and w.strip(";.,()") != first["sha256"]]
    assert other, err  # the new archive's sha256 is named too
    assert _tree(spool) == before and not stdout


def test_w28_the_run_must_be_closed_and_a_profile_is_required(tmp_path, capsys):
    prof = _profile(tmp_path)
    d = tmp_path / "runs" / f"run={'04' * 16}"
    d.mkdir(parents=True)
    code, _, err = _archive(capsys, prof, d)
    assert code == 2 and "not closed" in err, err
    code, _, err = _main(capsys, "results", "list")
    assert code == 2 and "--profile" in err, err


def test_w28_a_crash_before_the_rename_into_the_spool_leaves_nothing_visible(tmp_path, capsys):
    from biobox import profile
    from biobox.search import results

    prof = _profile(tmp_path)
    where = results.places(profile.load(prof))
    run = _make_run(tmp_path / "runs", 5)

    def crash(src, dst):
        raise KeyboardInterrupt("killed before the rename")

    with pytest.raises(KeyboardInterrupt):
        results.archive(run, [], where, publish=crash)
    assert _visible(where.spool) == []
    code, stdout, _ = _main(capsys, "results", "list", "--profile", str(prof))
    assert code == 0 and _rows(stdout) == []


def test_w28_outputs_are_kept_by_their_names_and_a_search_hmm_directory_of_another_run_is_refused(tmp_path, capsys):
    prof = _profile(tmp_path)
    out = _make_out(tmp_path / "rnap", 6)
    run = _make_run(tmp_path / "runs", 6, recipe_yaml="recipe: x\n")
    code, stdout, err = _archive(capsys, prof, run)  # no outputs: archived, and the note says how to add them
    assert code == 0 and "without outputs" in err and _rows(stdout)[0]["outputs"] == "-", err
    run6b = _make_run(tmp_path / "runs2", 6, recipe_yaml="recipe: x\n")
    prof2 = _profile(tmp_path / "p2")
    code, stdout, err = _archive(capsys, prof2, run6b, "--with", str(out), "--with", str(out / "hits.tsv"))
    assert code == 2 and "is inside the output" in err and not stdout, err  # it would be kept twice
    code, stdout, err = _archive(capsys, prof2, run6b, "--with", str(out))
    assert code == 0, err
    assert _rows(stdout)[0]["outputs"] == "rnap" and _rows(stdout)[0]["n_files"] == "9"
    with tarfile.open(tmp_path / "p2/tmp/results/spool" / f"{'06' * 16}.tar") as tf:
        assert "outputs/rnap/loci/locus_B_2.tsv" in tf.getnames()
    # a run whose recipe is not the one in that directory: the directory is another run's, refused
    run7 = _make_run(tmp_path / "runs", 7, recipe_yaml="recipe: other\n")
    code, stdout, err = _archive(capsys, prof, run7, "--with", str(out))
    assert code == 2 and "not the recipe run" in err and not stdout, err
    # outputs given twice under one name are refused; a symbolic link is refused
    run8 = _make_run(tmp_path / "runs", 8)
    twin = _make_out(tmp_path / "x" / "rnap", 8)
    code, _, err = _archive(capsys, prof, run8, "--with", str(out), "--with", str(twin))
    assert code == 2 and "two outputs are named 'rnap'" in err, err
    (tmp_path / "linked").mkdir()
    (tmp_path / "linked" / "hits.tsv").symlink_to(out / "hits.tsv")
    code, _, err = _archive(capsys, prof, run8, "--with", str(tmp_path / "linked"))
    assert code == 2 and "symbolic link" in err, err


# --- where ------------------------------------------------------------------------------------------------------


def test_w28_a_results_area_not_on_the_profiles_durable_storage_is_refused(tmp_path, capsys, monkeypatch):
    prof = _profile(tmp_path)
    durable = tmp_path / "durable"
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    code, _, err = _main(capsys, "results", "list", "--profile", str(prof), "--results", str(elsewhere))
    assert code == 2 and "not inside the profile's durable storage" in err, err
    # spelled inside durable, but a link leads out of it: the resolved path decides
    (durable / "away").symlink_to(elsewhere)
    code, _, err = _main(capsys, "results", "list", "--profile", str(prof), "--results", str(durable / "away" / "r"))
    assert code == 2 and "not inside" in err, err
    code, _, err = _main(capsys, "results", "list", "--profile", str(prof), "--results",
                         str(durable / ".." / "elsewhere"))
    assert code == 2 and "not inside" in err, err
    # inside durable, but another file system (a mount): refused by the device, whatever the path says
    from biobox.search import results

    real = os.stat
    mount = durable / "mounted"
    mount.mkdir()

    def stat(p, *a, **k):
        s = real(p, *a, **k)
        if Path(p) == mount.resolve():
            return os.stat_result((s.st_mode, s.st_ino, s.st_dev + 1, s.st_nlink, s.st_uid, s.st_gid, s.st_size,
                                   int(s.st_atime), int(s.st_mtime), int(s.st_ctime)))
        return s

    monkeypatch.setattr(results.os, "stat", stat)
    code, _, err = _main(capsys, "results", "list", "--profile", str(prof), "--results", str(mount / "r"))
    assert code == 2 and "another file system" in err, err
    monkeypatch.setattr(results.os, "stat", real)
    # inside durable: accepted, and the default is <durable>/results
    assert _main(capsys, "results", "list", "--profile", str(prof), "--results", str(durable / "mine"))[0] == 0
    assert _main(capsys, "results", "list", "--profile", str(prof))[0] == 0
    # a profile whose durable storage does not exist: refused, nothing made
    shutil.rmtree(durable)
    code, _, err = _main(capsys, "results", "list", "--profile", str(prof))
    assert code == 2 and "does not exist" in err and not durable.exists(), err


# --- the size rule and flush ------------------------------------------------------------------------------------


def test_w28_a_small_run_never_reaches_the_durable_area_and_a_bundle_over_the_floor_does(tmp_path, capsys):
    floor = 600_000  # lowered for the test; the pod's profile has 67108864
    prof = _profile(tmp_path, floor=floor)
    area = tmp_path / "durable" / "results"
    spool = tmp_path / "tmp/results/spool"
    runs = [_make_run(tmp_path / "runs", n, size=50_000) for n in (11, 12, 13)]  # about 200 kB each
    assert _archive(capsys, prof, runs[0])[0] == 0
    code, stdout, err = _main(capsys, "results", "flush", "--profile", str(prof))
    assert code == 0 and not stdout and "under the floor of 600000" in err and "nothing written" in err, err
    assert not area.exists() or os.listdir(area) == []  # nothing on durable storage, not even a temporary file
    assert len(_visible(spool)) == 2
    for r in runs[1:]:
        assert _archive(capsys, prof, r)[0] == 0
    sizes = [int(_rows(_archive(capsys, prof, r)[1])[0]["size"]) for r in runs]  # each a no-op now
    assert all(s < floor for s in sizes) and sum(sizes) >= floor
    code, stdout, err = _main(capsys, "results", "flush", "--profile", str(prof))
    assert code == 0, err
    bundle = Path(stdout.split("\t")[0])
    assert _visible(area) == [bundle.name] and bundle.stat().st_size >= floor
    assert bundle.name.startswith("results-20261103-20261103-") and "3 runs" in stdout
    assert hashlib.sha256(bundle.read_bytes()).hexdigest() in stdout
    assert _visible(spool) == []  # the spool's entries go once the bundle is in place
    with tarfile.open(bundle) as tf:
        assert tf.getnames() == ["MANIFEST.tsv", *(f"{bytes([n]).hex() * 16}.tar" for n in (11, 12, 13))]
    # flushing again: nothing to do; archiving a bundled run again: a no-op naming the bundle
    code, _, err = _main(capsys, "results", "flush", "--profile", str(prof))
    assert code == 0 and "holds no runs" in err
    code, stdout, err = _archive(capsys, prof, runs[0])
    assert code == 0 and str(bundle) in err and _visible(spool) == [], err


def test_w28_below_floor_writes_a_small_bundle_only_when_asked_and_says_it_breaks_the_rule(tmp_path, capsys):
    prof = _profile(tmp_path, floor=10**9)
    assert _archive(capsys, prof, _make_run(tmp_path / "runs", 14))[0] == 0
    code, stdout, err = _main(capsys, "results", "flush", "--profile", str(prof), "--below-floor")
    assert code == 0 and "PI's agreement" in err and "1 runs" in stdout, err
    assert len(_visible(tmp_path / "durable" / "results")) == 1


def test_w28_a_flush_killed_before_its_rename_leaves_nothing_visible_and_the_runs_spooled(tmp_path, capsys):
    from biobox import profile
    from biobox.search import results

    prof = _profile(tmp_path, floor=1000)
    where = results.places(profile.load(prof))
    run = _make_run(tmp_path / "runs", 15)
    assert _archive(capsys, prof, run)[0] == 0

    def crash(src, dst):
        raise KeyboardInterrupt("killed before the rename")

    with pytest.raises(KeyboardInterrupt):
        results.flush(where, publish=crash)
    assert _visible(where.area) == []  # nothing a reader takes
    partial = [n for n in os.listdir(where.area) if n.startswith(".")]
    assert len(partial) == 1 and ".part-" in partial[0]
    code, stdout, err = _main(capsys, "results", "list", "--profile", str(prof))
    assert code == 0 and [r["where"] for r in _rows(stdout)] == [str(where.spool / f"{'0f' * 16}.tar")]
    assert partial[0] in err and "ask the PI" in err, err
    code, stdout, err = _main(capsys, "results", "flush", "--profile", str(prof))  # again: written; partial kept
    assert code == 0 and len(_visible(where.area)) == 1, err
    assert partial[0] in os.listdir(where.area)


def test_w28_a_bundle_is_never_overwritten(tmp_path, capsys):
    prof = _profile(tmp_path, floor=1000)
    run = _make_run(tmp_path / "runs", 16)
    assert _archive(capsys, prof, run)[0] == 0
    spool = tmp_path / "tmp/results/spool"
    saved = _tree(spool)
    code, stdout, _ = _main(capsys, "results", "flush", "--profile", str(prof))
    bundle = Path(stdout.split("\t")[0])
    good = bundle.read_bytes()
    # the same spool again (as after a flush that stopped before removing it): the bundle is there, the same: no-op
    for name, data in saved.items():
        (spool / name).write_bytes(data)
    code, _, err = _main(capsys, "results", "flush", "--profile", str(prof))
    assert code == 0 and "in a bundle already" in err and _visible(spool) == [], err
    # a file of that name with other bytes (here another run's bundle): refused, and it is left as it is
    other = _profile(tmp_path / "other", floor=1000)
    assert _archive(capsys, other, _make_run(tmp_path / "other" / "runs", 17))[0] == 0
    code, stdout, _ = _main(capsys, "results", "flush", "--profile", str(other))
    foreign = Path(stdout.split("\t")[0]).read_bytes()
    for name, data in saved.items():
        (spool / name).write_bytes(data)
    (bundle.parent / "keep").mkdir()
    shutil.move(bundle, bundle.parent / "keep" / bundle.name)
    bundle.write_bytes(foreign)
    code, _, err = _main(capsys, "results", "flush", "--profile", str(prof))
    assert code == 2 and "never overwritten" in err and bundle.read_bytes() == foreign, err
    assert (bundle.parent / "keep" / bundle.name).read_bytes() == good
    assert len(_visible(spool)) == 2  # the run stays spooled
    # a file named as a bundle that is not one is reported, not read as one
    bundle.write_bytes(b"not a tar")
    code, _, err = _main(capsys, "results", "list", "--profile", str(prof))
    assert code == 1 and "not one" in err, err


# --- list and get -----------------------------------------------------------------------------------------------


def test_w28_list_reads_only_manifests_and_get_restores_every_file_with_its_sha256_checked(tmp_path, capsys,
                                                                                         monkeypatch):
    prof = _profile(tmp_path, floor=1000)
    a, b = _make_run(tmp_path / "runs", 21), _make_run(tmp_path / "runs", 22)
    out = _make_out(tmp_path / "flanks_nr50", 21)
    assert _archive(capsys, prof, a, "--with", str(out))[0] == 0
    assert _main(capsys, "results", "flush", "--profile", str(prof))[0] == 0
    assert _archive(capsys, prof, b)[0] == 0  # b stays spooled
    read: list[str] = []
    real = tarfile.TarFile.extractfile

    def spy(self, member):
        read.append(member if isinstance(member, str) else member.name)
        return real(self, member)

    monkeypatch.setattr(tarfile.TarFile, "extractfile", spy)
    code, stdout, err = _main(capsys, "results", "list", "--profile", str(prof))
    monkeypatch.setattr(tarfile.TarFile, "extractfile", real)
    assert code == 0, err
    assert read == ["MANIFEST.tsv"]  # the bundle's manifest; no run's archive is read
    rows = _rows(stdout)
    assert [r["run_id"] for r in rows] == ["15" * 16, "16" * 16]
    assert Path(rows[0]["where"]).parent == tmp_path / "durable" / "results"
    assert rows[1]["where"] == str(tmp_path / "tmp/results/spool" / f"{'16' * 16}.tar")
    assert rows[0]["outputs"] == "flanks_nr50" and rows[0]["run_date"] == "2026-11-03T09:30:05Z"
    assert "1 runs" in err  # the spool's summary
    code, stdout, _ = _main(capsys, "results", "list", "--profile", str(prof), "--json")
    assert [json.loads(x)["kept"] for x in stdout.splitlines()] == ["bundle", "spool"]
    for rid, run, outs in (("15" * 16, a, {"flanks_nr50": out}), ("16" * 16, b, {})):
        dest = tmp_path / "back" / rid
        code, stdout, err = _main(capsys, "results", "get", rid, "--profile", str(prof), "--out", str(dest))
        assert code == 0, err
        got = _tree(dest)
        assert got.pop("FILES.tsv").startswith(b"# bio-box results archive 1\n")
        want = {f"run/{k}": v for k, v in _tree(run).items()}
        for name, d in outs.items():
            want.update({f"outputs/{name}/{k}": v for k, v in _tree(d).items()})
        assert got == want
        assert not [n for n in os.listdir(dest.parent) if n.startswith(".")]  # no temporary left beside it
    code, _, err = _main(capsys, "results", "get", "ab" * 16, "--profile", str(prof), "--out", str(tmp_path / "n"))
    assert code == 2 and "in no bundle" in err and not (tmp_path / "n").exists()
    code, _, err = _main(capsys, "results", "get", "15" * 16, "--profile", str(prof), "--out", str(tmp_path / "back"))
    assert code == 2 and "not empty" in err


def test_w28_get_refuses_an_archive_whose_bytes_are_not_its_manifests(tmp_path, capsys):
    prof = _profile(tmp_path)
    assert _archive(capsys, prof, _make_run(tmp_path / "runs", 23))[0] == 0
    tar = tmp_path / "tmp/results/spool" / f"{'17' * 16}.tar"
    data = bytearray(tar.read_bytes())
    data[3000] ^= 0xFF  # inside a member's bytes
    tar.write_bytes(bytes(data))
    dest = tmp_path / "back"
    code, _, err = _main(capsys, "results", "get", "17" * 16, "--profile", str(prof), "--out", str(dest))
    assert code == 1 and "corrupt" in err and "sha256" in err, err
    assert not dest.exists() and not [n for n in os.listdir(tmp_path) if n.startswith(".back")]
    # a FILES.tsv check on its own: an archive with a right sha256 in its line but a member that FILES.tsv does
    # not describe is refused too
    from biobox.search import results

    bad = tmp_path / "bad.tar"
    with tarfile.open(bad, "w", format=tarfile.PAX_FORMAT) as tf:
        files = b"# bio-box results archive 1\npath\tsize\tsha256\nrun/x\t1\t" + b"0" * 64 + b"\n"
        for name, payload in (("FILES.tsv", files), ("run/x", b"y")):
            ti = tarfile.TarInfo(name)
            ti.size = len(payload)
            tf.addfile(ti, io.BytesIO(payload))
    (tmp_path / "into").mkdir()
    with pytest.raises(results.Corrupt, match="sha256"):
        results._restore(bad, tmp_path / "into")
