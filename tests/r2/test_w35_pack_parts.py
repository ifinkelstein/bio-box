"""Ring R2: the search pack built in parts under the CPU budget (card W-35; AGENTS §7.2, §4.1, §4.3).

The oracle of every byte-identity test is the one-process build (`pack.build`) of the same release: its files by
sha256 and its pack id, which the W-16 tests check against an in-memory oracle of their own. Two releases: the
micro-S0 world (three contigs, a domain library, built as in test_micro_s0_search) and the synthetic store of
`scripts/stress/synth_scale.py` at 10 copies (30 contigs), so that 2, 3 and 5 parts cut it in different places.
Parts run as their own processes (the default runner), except in the test that lowers the sort's caps, where they
run in this process so the caps reach them. The allowance of the hour is fixed here (`ALLOW`): the pod's real
reading counts other agents' tests, which would hold a launch.

The estimate is checked against brute-force counts of the fixture's tables (its genome bases and protein rows) and,
for its rates, against the run records of the build of `rN.3` in parts (measured on wilkcomp02, the numbers in
`planning/measurements/2026-11-14_pack_build_in_parts.md`).
"""

from __future__ import annotations

import contextlib
import hashlib
import io as stdio
import json
import runpy
import signal
import subprocess
import sys
import time
from pathlib import Path

import pyarrow.parquet as pq
import pytest
import yaml

from biobox import catalog, cli, profile
from biobox.l1 import igr
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.search import pack

sys.path.insert(0, str(Path(__file__).parent))
import test_micro_s0_l3 as micro  # noqa: E402

pytestmark = [pytest.mark.r2, pytest.mark.covers("W-35")]

REPO = Path(__file__).resolve().parents[2]
REL = micro.REL
SYN = "rS0.syn"  # synth_scale's release
FIX = REPO / "tests/fixtures/micro_s0"


def ALLOW(threads: int) -> int:
    """A fixed allowance of the hour (module docstring)."""
    return 16


def _profile(d: Path, cpu_min: int = 300) -> Path:
    data = {"max_threads_day": 4, "max_threads_night": 4, "mem_scope": "4G", "tmp": str(d / "ptmp"),
            "durable": str(d / "durable"), "bulk": str(d / "bulk"), "nice_monolithic": 19, "nice_fanout": 15,
            "day_start_hour": 8, "day_end_hour": 20, "cpu_min_per_process": cpu_min,
            "small_file_floor_bytes": 1000, "require_mem_cap": False}
    path = d / "profile.yaml"
    path.write_text(json.dumps(data))  # JSON is YAML
    return path


def _files(pdir: Path) -> dict[str, str]:
    """Every file of a pack directory by its path in the pack: sha256."""
    return {p.relative_to(pdir).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(pdir.rglob("*")) if p.is_file()}


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """The micro-S0 world with its IGRs, and its pack built in one process: (root, hmm, packs, one-process files)."""
    d = tmp_path_factory.mktemp("w35")
    root = d / "a"
    circ = micro.build(root)
    igr.run(root, "S0", circ)
    hmm = d / "a.hmm"
    man = pack.build(root, REL, models={"test-1": hmm}, packs=d / "one")
    return root, hmm, man["pack_id"], _files(pack.pack_dir(d / "one", REL, man["pack_id"]))


@pytest.fixture(scope="module")
def syn(tmp_path_factory):
    """`synth_scale.py` at 10 copies (as the W-16 memory test runs it) with its IGRs, and its one-process pack."""
    d = tmp_path_factory.mktemp("w35syn")
    script = REPO / "scripts/stress/synth_scale.py"
    argv = sys.argv
    sys.argv = [str(script), str(REPO), str(d), "10", "--shard-bases", str(2**26), "--seed", "1"]
    try:
        with contextlib.redirect_stdout(stdio.StringIO()):
            g = runpy.run_path(str(script), run_name="__main__")
    finally:
        sys.argv = argv
    g["con"].close()
    root = Path(g["root"])
    igr.run(root, "S0")
    man = pack.build(root, SYN, packs=d / "one")
    return root, None, man["pack_id"], _files(pack.pack_dir(d / "one", SYN, man["pack_id"]))


def _in_parts(store, tmp_path: Path, n: int, jobs: int, release: str, **kw):
    root, hmm, _, _ = store
    prof = profile.load(_profile(tmp_path))
    models = {"test-1": hmm} if hmm is not None else None
    return pack.build_in_parts(root, release, models=models, packs=tmp_path / "packs", prof=prof, parts=n,
                               jobs=jobs, allow=ALLOW, **kw)


@pytest.mark.parametrize("which,n,jobs", [("world", 1, 1), ("world", 2, 1), ("world", 3, 3), ("world", 5, 3),
                                          ("syn", 2, 1), ("syn", 3, 3), ("syn", 5, 3), ("syn", 5, 1)])
def test_w35_a_pack_built_in_parts_is_the_one_process_pack_byte_for_byte(world, syn, tmp_path, which, n, jobs):
    """Every file of the pack (shards, meta files, models, manifest, the verification record) has the sha256 of the
    one-process build's, so the pack id is the same; nothing is left in the pack directory but the pack, and the work
    directory is gone."""
    store = world if which == "world" else syn
    release = REL if which == "world" else SYN
    code, man = _in_parts(store, tmp_path, n, jobs, release)
    assert code == 0, code
    assert man["pack_id"] == store[2]
    pdir = pack.pack_dir(tmp_path / "packs", release, man["pack_id"])
    assert _files(pdir) == store[3]
    assert [p.name for p in pack.release_dir(tmp_path / "packs", release).iterdir()] == [pdir.name]
    work = tmp_path / "ptmp" / "scratch" / "pack" / f"release={release}"
    assert not work.exists() or list(work.iterdir()) == [], list(work.iterdir())


def test_w35_the_cut_puts_each_contig_in_exactly_one_part_and_follows_the_bases(syn):
    """The 30 contigs of the synthetic store are split into 2, 3 and 5 consecutive ranges of contig uids: every contig
    of the tapes is in exactly one part (a brute-force loop over the tape files), and parts hold bases near 1/N of the
    release's each (within the largest contig)."""
    root = syn[0]
    rf = catalog.ReleaseFiles.of(root, SYN)
    lengths = {}
    for f in rf.table("contig_tapes", tape_version=pack.tape_version([V])):
        for r in pq.read_table(f, columns=["contig_uid", "contig_len"]).to_pylist():
            lengths[r["contig_uid"]] = r["contig_len"]
    total, big = sum(lengths.values()), max(lengths.values())
    assert len(lengths) == 30
    for n in (2, 3, 5):
        parts = pack._cut(rf, pack.tape_version([V]), n)
        assert len(parts) == n
        seen = []
        for p in parts:
            rng = pack._Range.of(p)
            mine = sorted(c for c in lengths if rng.has(c))
            assert p["contigs"] == len(mine) and p["bases"] == sum(lengths[c] for c in mine)
            seen += mine
            assert abs(p["bases"] - total / n) <= big, (n, p["bases"], total / n)
        assert sorted(seen) == sorted(lengths)  # each contig once


def test_w35_many_sorted_runs_merge_passes_shards_and_row_groups_in_parts(world, tmp_path, monkeypatch):
    """With the sort's caps far below the fixture's size (as the W-16 test sets them), each part goes through several
    sorted runs and merge passes and the last process cuts several shards and meta row groups: the pack in 3 parts,
    run in this process so the caps reach them, equals the one-process pack built under the same caps."""
    root, hmm, _, _ = world
    for k, v in (("RUN_ROWS", 50), ("FAN_IN", 3), ("MERGE_BATCH", 16), ("META_ROW_GROUP_ROWS", 64), ("BATCH_ROWS", 7),
                 ("TAPE_BATCH_ROWS", 1)):
        monkeypatch.setattr(pack, k, v)
    one = pack.build(root, REL, models={"test-1": hmm}, shard_text_bytes=20_000, packs=tmp_path / "one")
    want = _files(pack.pack_dir(tmp_path / "one", REL, one["pack_id"]))
    assert sum(1 for s in one["shards"] if s["target"] == "dark_orf") > 1

    def here(prof, root, stage, step, cmd, *, threads):  # the part in this process: its argv after `-m <module>`
        return pack.main(list(cmd[3:]))

    prof = profile.load(_profile(tmp_path))
    code, man = pack.build_in_parts(root, REL, models={"test-1": hmm}, shard_text_bytes=20_000, packs=tmp_path / "p",
                                    prof=prof, parts=3, jobs=2, allow=ALLOW, run=here)
    assert code == 0 and man["pack_id"] == one["pack_id"]
    assert _files(pack.pack_dir(tmp_path / "p", REL, man["pack_id"])) == want


# --- the estimate ----------------------------------------------------------------------------------------------------


def _brute_counts(root: Path, release: str) -> tuple[int, int, int]:
    """(protein_seq rows, genome bases, contigs) of the release by a loop over its tables."""
    rf = catalog.ReleaseFiles.of(root, release)
    rows = sum(pq.read_table(f).num_rows for f in rf.table("protein_seq", gene_call_version=V))
    lengths = {}
    for f in rf.table("contig_tapes", tape_version=pack.tape_version([V])):
        for r in pq.read_table(f, columns=["contig_uid", "contig_len"]).to_pylist():
            lengths[r["contig_uid"]] = r["contig_len"]
    return rows, sum(lengths.values()), len(lengths)


def test_w35_the_counts_of_the_estimate_equal_a_brute_force_loop(world, syn):
    """The estimate's counts: protein_seq rows and genome bases (once per contig), from a loop over the tables; on
    the world, the bases are also the lengths `sources.yaml` lists for the three genomes."""
    for store, release in ((world, REL), (syn, SYN)):
        c = pack.release_counts(store[0], release)
        assert (c.proteins, c.bases, c.contigs) == _brute_counts(store[0], release)
    src = yaml.safe_load((FIX / "sources.yaml").read_text())
    assert pack.release_counts(world[0], REL).bases == sum(o["length"] for o in src["objects"]) == 797_481


def test_w35_pack_build_prints_its_estimate_before_any_work(world, tmp_path, capsys):
    root, hmm, pid, _ = world
    prof = _profile(tmp_path)
    capsys.readouterr()
    assert cli.main(["pack", "build", "--root", str(root), "--release", REL, "--model", f"test-1={hmm}",
                     "--pack-dir", str(tmp_path / "packs"), "--profile", str(prof)]) == 0
    got = capsys.readouterr()
    assert {json.loads(x)["pack_id"] for x in got.out.splitlines()} == {pid}
    msgs = [x.split(" pack build: ", 1)[1] for x in got.err.splitlines() if " pack build: " in x]
    est = [m for m in msgs if m.startswith(f"release {REL}: estimate")]
    assert len(est) == 1 and "in one process" in est[0] and "797,481 genome bases" in est[0], est
    assert msgs.index(est[0]) < next(i for i, m in enumerate(msgs) if "checking genetic codes" in m)


# the dark-ORF rate that puts the world over the budget in one process but leaves its largest contig, M. genitalium
# (580,076 of 797,481 bases), under it: 203 of 240 CPU-minutes alone, 279 in all (derived)
HIGH_RATE = 3.5e5


def test_w35_over_the_margin_without_parts_is_refused_and_parts_that_fit_build_the_same_pack(world, tmp_path, capsys,
                                                                                         monkeypatch):
    """With the dark-ORF rate raised so that the world's estimate passes 80% of a 300 CPU-minute budget: `pack
    build` without `--parts` exits 2 naming `--parts`, before any work (nothing in the pack directory); with `--parts
    1` the one dark-ORF part is over the margin and refused, naming its estimate; with `--parts 3` (each contig then in
    a part of its own or beside a small one, by the cut's rule) it builds the one-process pack."""
    root, hmm, pid, files = world
    monkeypatch.setattr(pack, "DARK_ORF_CPU_MIN_PER_GBP", HIGH_RATE)
    monkeypatch.setattr(pack, "_allowance", lambda prof: ALLOW)
    prof = _profile(tmp_path)
    packs = tmp_path / "packs"
    base = ["pack", "build", "--root", str(root), "--release", REL, "--model", f"test-1={hmm}", "--pack-dir",
            str(packs), "--profile", str(prof)]
    capsys.readouterr()
    assert cli.main(base) == 2
    err = capsys.readouterr().err
    assert "refused" in err and "--parts 2" in err and "cpu_min_per_process of 300" in err, err
    assert "checking genetic codes" not in err and not packs.exists()
    assert cli.main([*base, "--parts", "1", "--jobs", "1"]) == 2
    err = capsys.readouterr().err
    assert "refused" in err and "dark_orf-0-of-1" in err, err
    assert _empty(pack.release_dir(packs, REL))
    assert cli.main([*base, "--parts", "3", "--jobs", "2"]) == 0
    got = capsys.readouterr()
    assert {json.loads(x)["pack_id"] for x in got.out.splitlines()} == {pid}, got.err
    assert _files(pack.pack_dir(packs, REL, pid)) == files


def test_w35_parts_without_a_profile_or_with_work_under_the_root_are_refused(world, tmp_path, capsys, monkeypatch):
    root, hmm, _, _ = world
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)
    base = ["pack", "build", "--root", str(root), "--release", REL, "--pack-dir", str(tmp_path / "packs"),
            "--parts", "2"]
    assert cli.main(base) == 2
    assert "--profile" in capsys.readouterr().err
    assert cli.main([*base, "--profile", str(_profile(tmp_path)), "--work", str(root / "w")]) == 2
    assert "under the data root" in capsys.readouterr().err
    assert not (root / "w").exists()


# --- resuming, STOP and killed parts ---------------------------------------------------------------------------------


def _empty(d: Path) -> bool:
    return not d.exists() or not any(d.iterdir())


def _copy(root: Path, to: Path) -> Path:
    import shutil

    shutil.copytree(root, to)
    return to


def _subprocess_runner(launched: list[str], after=None):
    """A runner like the default one (each part its own process, here without nice or scope), that lists the steps it
    launches and calls `after(step)` once each has exited 0."""

    def run(prof, root, stage, step, cmd, *, threads):
        launched.append(step)
        code = subprocess.run(list(cmd), check=False).returncode
        code = code if code >= 0 else 128 - code
        if code == 0 and after is not None:
            after(step)
        return code

    return run


def test_w35_stop_between_parts_then_the_next_call_runs_only_the_parts_without_output(world, tmp_path):
    """`STOP` in the data root after the first part: the driver launches no more and exits 75, with that part's output
    kept and nothing in the pack directory; without `STOP` the next call launches only the other parts and the last
    process, and gives the one-process pack."""
    root = _copy(world[0], tmp_path / "root")
    store = (root, *world[1:])
    launched: list[str] = []
    code, man = _in_parts(store, tmp_path, 3, 1, REL,
                          run=_subprocess_runner(launched, after=lambda s: (root / "STOP").write_text("")))
    assert code == 75 and man is None
    assert launched == ["unique"]
    assert _empty(pack.release_dir(tmp_path / "packs", REL))
    (root / "STOP").unlink()
    again: list[str] = []
    code, man = _in_parts(store, tmp_path, 3, 2, REL, run=_subprocess_runner(again))
    assert code == 0 and man["pack_id"] == world[2]
    assert "unique" not in again and sorted(again) == sorted(
        [f"{t}-{k}-of-3" for t in ("igr", "dark_orf") for k in range(3)] + ["assemble"])
    assert _files(pack.pack_dir(tmp_path / "packs", REL, man["pack_id"])) == world[3]


_HELD = r"""
import signal, sys
from pathlib import Path
from biobox.search import pack
flag, point = sys.argv[-2:]
log = pack._log


def held(msg):  # the part stops here until a signal comes
    log(msg)
    if point in msg:
        Path(flag).write_text(msg)
        signal.pause()


pack._log = held
sys.exit(pack.main(sys.argv[1:-2]))
"""


def _poll(pred, timeout: float = 120.0) -> bool:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if pred():
            return True
        time.sleep(0.05)
    return False


@pytest.mark.parametrize("victim,point,sig", [("dark_orf-1-of-3", "genes read", signal.SIGTERM),
                                              ("dark_orf-1-of-3", "genes read", signal.SIGKILL),
                                              ("igr-2-of-3", "igr: start", signal.SIGTERM),
                                              ("assemble", "dark_orf: merging", signal.SIGTERM)])
def test_w35_a_killed_part_leaves_nothing_a_reader_takes_and_the_next_call_finishes_the_pack(
        world, tmp_path, victim, point, sig):
    """A part (or the last process) is held at a point of its work, then sent a signal: the driver exits with the
    part's code and launches nothing more; the part's output directory does not exist (by SIGTERM its temporary
    directory is removed too; by SIGKILL it is left, dot-named) and nothing is in the pack directory. The next call
    removes what SIGKILL left, runs that part again and gives the one-process pack."""
    root, hmm, pid, files = world
    flag = tmp_path / "held"
    work: list[Path] = []

    def run(prof, root_, stage, step, cmd, *, threads):
        if step != victim:
            return _subprocess_runner([])(prof, root_, stage, step, cmd, threads=threads)
        plan = Path(cmd[cmd.index("--plan") + 1])
        work.append(plan.parent)
        p = subprocess.Popen([sys.executable, "-c", _HELD, *cmd[3:], str(flag), point])
        try:
            assert _poll(flag.exists), f"{step} never reached {point}"
            tmp = [x for x in plan.parent.iterdir() if x.name.startswith((f".{step}.tmp-", ".building-"))]
            assert tmp, list(plan.parent.iterdir())  # it was writing when it was stopped
            p.send_signal(sig)
            code = p.wait(timeout=120)
        finally:
            if p.poll() is None:
                p.kill()
        return code if code >= 0 else 128 - code

    code, man = _in_parts(world, tmp_path, 3, 1, REL, run=run)
    assert man is None and code == (143 if sig == signal.SIGTERM else 137), code
    wd = work[0]
    assert not (wd / victim).exists()
    left = [x.name for x in wd.iterdir() if x.name.startswith((f".{victim}.tmp-", ".building-", ".scratch-"))]
    assert (left != []) == (sig == signal.SIGKILL), left
    assert _empty(pack.release_dir(tmp_path / "packs", REL))
    code, man = _in_parts(world, tmp_path, 3, 2, REL)
    assert code == 0 and man["pack_id"] == pid
    assert _files(pack.pack_dir(tmp_path / "packs", REL, pid)) == files
    assert not wd.exists()


def test_w35_a_part_records_its_cpu_wall_time_peak_and_estimate(world, tmp_path, capsys):
    """Each part's line on stderr gives its CPU-minutes, wall time, peak RSS and estimate (card W-35: the pod run
    reports them for every part)."""
    capsys.readouterr()
    code, _ = _in_parts(world, tmp_path, 2, 1, REL)
    assert code == 0
    err = capsys.readouterr().err
    for name in ("unique", "igr-0-of-2", "igr-1-of-2", "dark_orf-0-of-2", "dark_orf-1-of-2", "assemble"):
        line = [x for x in err.splitlines() if f"part {name} done:" in x]
        assert len(line) == 1, (name, err[-3000:])
        assert "CPU-minutes" in line[0] and "wall" in line[0] and "peak RSS" in line[0] and "estimate" in line[0]
