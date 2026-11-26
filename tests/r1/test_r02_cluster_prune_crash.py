"""Ring R1: clustering scratch pruning survives a kill, a format-1 plan can be redone, and the method is recorded
(card R-02; DEC-77, DEC-24; audit of 2026-11-25, §2.2).

(a), (b) `clean_consumed` records a prune before it removes anything, and the resume path of build-stage
(`biobox.stage._cluster_in_parts`: `lost()` and its redo, then `clean_consumed` again for every call) finishes a
recorded prune and prunes again what a redo wrote again. The resume path is stage.py's own; only the run-step of
each step is stood in for (`World.one`: a step with a run record is passed over, as run-step passes it, else it runs
in this process), and the redo launcher `_launch` (`World.launch`), whose launches are counted: a resume after a kill
must launch none. The kill is an exception raised from a patched `shutil.rmtree`, `_clear` or `_write_pruned` while
`clean_consumed` runs for the step under test. Without MMseqs2 the steps are made by the test (`World.fake_step`: each
call's outputs are small files and its record their digests, as `step` records them), so the ordering tests run on
GitHub Actions; with MMseqs2 they are the real steps, and the tables are compared with the one-process `run` (an
independent driver: `cluster.run`).

(c) A format-1 work directory is written with the format-1 code path from its plan on (the plan record as the code
before card B-71 wrote it, every step run under it), and is made again from nothing to its records. The cut record
of a call sharing an index must hold the fields of the code before card B-71 (`git show b741c28`, `_execute` and
`_cut_doc`), not `createindex_cpu_min`.

(d) DEC-77 (2): the clustering method is in the plan record and in the release's clustering ledger (`method_rows`);
a cascade ledger is the ledger of before card B-71.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import random
import shutil
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import ids, io, ledger, profile
from biobox import stage as bs
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l2 import cluster
from biobox.registry.schemas import table

pytestmark = pytest.mark.covers("R-02", "B-72", "B-71", "DEC-77", "DEC-24")

NEED_MMSEQS = pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs not on PATH")
REL = "rR02"
PROFILE = Path(__file__).resolve().parents[2] / "profiles" / "wilkcomp02.yaml"
# the fields of a shared-index cut record cut into more than one part, as the code before card B-71 wrote them
# (git show b741c28:src/biobox/l2/cluster.py: `_execute` for a cut, `_cut_doc`, `_make_index`)
PRE_B71_SHARED_CUT_KEYS = {"format", "release", "step", "plan_sha256", "module", "inputs", "split", "units", "cpu_min",
                           "overhead_cpu_min", "parts", "index_cpu_min", "index"}
PRE_B71_PLAN_KEYS = {"format", "release", "stages", "gene_call_version", "mmseqs2", "files", "cut", "proteins",
                     "release_cpu_min", "outputs"}


class Killed(Exception):
    """The process dies here."""


def _triggers(mode: str) -> list[str]:
    """Every step whose `clean_consumed` prunes something under `mode`, in step order."""
    names = [c.name for c in cluster.calls(mode)] + list(cluster.FINISH_LEVEL_STEPS)
    return [n for n in names if cluster._prune_spec(n, mode) is not None]


KILLS = ("write_pruned", "first", "last")


def _kills(mode: str) -> list[tuple[str, str, str]]:
    """(trigger, where the kill is, mode) for every prune of `mode`: as `_write_pruned` is entered; in the removal of
    its first target, before it is removed (`shutil.rmtree` for a directory, `_clear` for a database); in the removal
    of its last target, after it is removed (the worst case for a record written last)."""
    return [(t, where, mode) for t in _triggers(mode) for where in KILLS]


def _snapshot(w: Path, skip_logs: bool = False) -> dict[str, str]:
    return {p.relative_to(w).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(Path(w).rglob("*")) if p.is_file() and not p.is_symlink()
            and not (skip_logs and (p.name == "mmseqs.log" or p.relative_to(w).parts[0] == "logs"))}


def _tables(base: Path) -> dict[str, str]:
    base = Path(base)
    dirs = [f"l2/cluster_member/release={REL}", f"l2/cluster_member_aln/release={REL}",
            f"l2/member_cluster/release={REL}", f"ledger/release={REL}/step=cluster"]
    got = {}
    for d in dirs:
        for p in sorted((base / d).glob("*.parquet")):
            got[p.relative_to(base).as_posix()] = hashlib.sha256(p.read_bytes()).hexdigest()
    return got


class World:
    """One release clustered in parts by build-stage's `_cluster_in_parts`, its run-step stood in for (`one`)."""

    def __init__(self, tmp_path: Path, mode: str, *, real: bool = False, proteins: Path | None = None):
        self.root = (tmp_path / "root").resolve()
        if proteins is not None:
            shutil.copytree(proteins, self.root)
        else:
            self.root.mkdir(parents=True)
        self.scratch = (tmp_path / "scratch").resolve()
        self.mode, self.real = mode, real
        self.prof = dataclasses.replace(profile.load(PROFILE), tmp=str(tmp_path / "prof"))
        self.b = bs.Build(root=self.root, stage="S0", genomes=tmp_path / "g.json", objects=tmp_path / "o.json",
                          release=REL, stages=("S0",), tmp=self.scratch, threads=1, cluster_parts=1,
                          cluster_method=mode)
        self.by = {c.name: c for c in cluster.calls(mode)}
        self.executed: list[str] = []
        self.launched: list[str] = []
        self.final: dict[str, str] | None = None

    @property
    def w(self) -> Path:
        return cluster.work_dir(self.scratch, self.root, REL)

    # --- the stand-ins ---------------------------------------------------------------------------------------

    def one(self, step: bs.Step, **kw) -> int | None:
        rec = bs._record(self.root, step)
        if rec.exists():  # run-step passes a published step over
            return None
        self.executed.append(step.name)
        (self.real_step if self.real else self.fake_step)(step.name)
        rec.parent.mkdir(parents=True, exist_ok=True)
        rec.write_text(json.dumps({"exit_code": 0, "stage": REL, "step": step.name}))
        return None

    def launch(self, prof, root, todo, jobs, cpus, run, hook, allow):  # `bs._launch` of `_redo_lost_cluster`
        for s in todo:
            self.launched.append(s.name)
            name = s.name.removesuffix(".redo")
            if self.real:
                assert cluster.redo(self.root, REL, name, self.scratch, 1)["verified"]
            else:
                self.fake_outputs(name)
            if hook is not None:
                hook(f"after:{s.name}")
        return None

    def go(self, hook=None) -> int | None:
        def no_parts(*a, **k):
            raise AssertionError("cluster_parts=1: no part step is run")

        return bs._cluster_in_parts(self.prof, self.b, self.root, self.one, hook, 1, 1, no_parts, lambda t: t)

    # --- the steps, real ---------------------------------------------------------------------------------------

    def real_step(self, name: str) -> None:
        if name == cluster.PLAN_STEP:
            cluster.plan(self.root, REL, V, ["S0"], self.scratch, parts=1, mode=self.mode)
        elif name == cluster.FINAL_STEP:
            self.final = _snapshot(self.w, skip_logs=True)
            cluster.finish(self.root, REL, self.scratch, 1)
        else:
            cluster.step(self.root, REL, name, self.scratch, 1)

    # --- the steps, made by the test (no MMseqs2) ---------------------------------------------------------------

    def _plan_doc(self) -> dict:
        return {"format": cluster.PLAN_FORMAT, "release": REL, "stages": ["S0"], "gene_call_version": V,
                "mmseqs2": "fake", "files": [], "cut": {"parts": 1}, "method": self.mode}

    def _paths(self, name: str) -> list[str]:
        if name == cluster.PLAN_STEP:
            return ["db", "db_h", "prot.parquet"]
        if name in cluster.FINISH_LEVEL_STEPS:
            return list(cluster._finish_level_outputs(cluster._FINISH_LEVEL_MAP[name]))
        c = self.by[name]
        return [c.output] + ([c.order_file] if c.order_file else [])

    def fake_outputs(self, name: str) -> dict[str, str]:
        """Write step `name`'s outputs (the same bytes each time, as a redo reproduces them); their digests."""
        got = {}
        for rel in self._paths(name):
            p = self.w / rel if name != cluster.PLAN_STEP else self._w0 / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(f"{name} {rel}\n".encode())
            got[rel] = hashlib.sha256(p.read_bytes()).hexdigest()
        return got

    def fake_step(self, name: str) -> None:
        if name == cluster.PLAN_STEP:
            doc = self._plan_doc()
            self._w0 = cluster._work_dir(self.scratch, self.root, cluster._plan_input_dict(doc))
            self._w0.mkdir(parents=True)
            doc |= {"proteins": 3, "release_cpu_min": 0.0, "outputs": self.fake_outputs(name)}
            cluster._write_record(self.root, self.root, REL, name, doc)
            return
        if name == cluster.FINAL_STEP:
            self.final = _snapshot(self.w)
            return
        plan_sha = hashlib.sha256(cluster._record_raw(self.root, REL, cluster.PLAN_STEP)).hexdigest()
        doc = {"format": cluster.PLAN_FORMAT, "release": REL, "step": name, "plan_sha256": plan_sha}
        if name.endswith(".cut"):
            doc |= {"module": self.by[name[:-4]].module,
                    "parts": [{"first": 0, "last": 0, "keys": 1, "units": 1, "cpu_min": 0.0}]}
        else:
            doc |= {"module": "finish" if name in cluster.FINISH_LEVEL_STEPS else self.by[name].module,
                    "outputs": self.fake_outputs(name)}
        cluster._write_record(self.root, self.root, REL, name, doc)

    def records(self) -> dict[str, bytes]:
        d = self.root / cluster.plan_dir(REL)
        return {p.name: p.read_bytes() for p in sorted(d.glob("*.json"))}


def _arm(monkeypatch, trigger: str, where: str, mode: str) -> None:
    """Raise Killed once, inside `clean_consumed(trigger)`, at `where` (`KILLS`)."""
    targets = cluster._prune_spec(trigger, mode)[2]
    kind, rel = targets[0] if where == "first" else targets[-1]
    state = {"in": False, "fired": False}
    real_clean, real_rmtree, real_clear, real_write = (cluster.clean_consumed, shutil.rmtree, cluster._clear,
                                                       cluster._write_pruned)

    def clean(root, release, step_name, tmp, out=None):
        state["in"] = step_name == trigger and not state["fired"]
        try:
            return real_clean(root, release, step_name, tmp, out)
        finally:
            state["in"] = False

    def fire() -> None:
        state["fired"], state["in"] = True, False
        raise Killed(f"{where} in clean_consumed({trigger})")

    def removal(hit: bool, remove) -> None:
        if hit and state["in"] and where == "first":
            fire()
        remove()
        if hit and state["in"] and where == "last":
            fire()

    def rmtree(path, *a, **k):
        removal(kind == "dir" and Path(path).as_posix().endswith("/" + rel), lambda: real_rmtree(path, *a, **k))

    def clear(w, rel_):
        removal(kind == "db" and rel_ == rel, lambda: real_clear(w, rel_))

    def write(w, entries):
        if state["in"] and where == "write_pruned":
            fire()
        real_write(w, entries)

    monkeypatch.setattr(cluster, "clean_consumed", clean)
    monkeypatch.setattr(shutil, "rmtree", rmtree)
    monkeypatch.setattr(cluster, "_clear", clear)
    monkeypatch.setattr(cluster, "_write_pruned", write)


def _crash_and_resume(tmp_path, monkeypatch, trigger, where, mode, *, ref: World | None = None, real=False,
                      proteins=None) -> tuple[World, World]:
    """An uninterrupted clustering (`ref`, made here when not given), and one killed inside `clean_consumed(trigger)`
    then resumed by build-stage's path, which must make nothing again; returns both."""
    if ref is None:
        ref = World(tmp_path / "ref", mode, real=real, proteins=proteins)
        monkeypatch.setattr(bs, "_launch", ref.launch)
        assert ref.go() is None
        assert ref.launched == []

    got = World(tmp_path / "got", mode, real=real, proteins=proteins)
    monkeypatch.setattr(bs, "_launch", got.launch)
    targets = cluster._prune_spec(trigger, mode)[2]
    with monkeypatch.context() as m:
        _arm(m, trigger, where, mode)
        with pytest.raises(Killed):
            got.go()
    w = got.w
    mine = [e for e in cluster._read_pruned(w) if e["step"] == trigger]
    if where == "write_pruned":  # killed before the record: nothing removed yet
        assert mine == [] and all(cluster._present(w, k, rel) for k, rel in targets), \
            "a target was removed before the prune was recorded"
    else:  # killed while removing: the prune is recorded
        assert len(mine) == 1, f"{trigger}: a target was removed and no prune is recorded"
    if where == "first":
        assert all(cluster._present(w, k, rel) for k, rel in targets)
    assert cluster.lost(got.root, REL, got.scratch) == []  # nothing to make again

    assert got.go() is None  # the resume
    assert got.launched == [], f"the resume made published steps again: {got.launched}"
    assert sorted(got.executed) == sorted(ref.executed)
    assert len(got.executed) == len(set(got.executed))  # no step ran twice
    return ref, got


# --- (a) a kill inside clean_consumed, without MMseqs2 ---------------------------------------------------------------


@pytest.mark.parametrize("trigger,where,mode", _kills("hybrid-b1") + _kills("cascade"))
def test_r02_a_kill_inside_clean_consumed_resumes_without_redo(tmp_path, monkeypatch, trigger, where, mode):
    """A kill before the prune's record or while removing its targets, at every prune of both modes: no step is made
    again on the resume, every prune is finished, and the scratch before the tables and every record are those of an
    uninterrupted clustering."""
    ref, got = _crash_and_resume(tmp_path, monkeypatch, trigger, where, mode)
    assert got.final == ref.final  # the scratch before the tables: every prune done, .pruned.json the same
    assert got.records() == ref.records()
    for t in _triggers(mode):
        for _, rel in cluster._prune_spec(t, mode)[2]:
            assert rel not in got.final and not any(k.startswith(rel + "/") for k in got.final), f"{rel} left"
    assert not got.w.exists()  # removed once the tables are written


def test_r02_the_prune_record_keeps_the_b72_form(tmp_path, monkeypatch):
    """`.pruned.json` keeps the form card B-72 wrote (one entry per prune: step, milestones, steps; sorted keys, indent
    1), so a work directory pruned by the code before this card is read as before; an entry written in that form by
    hand is honoured by `lost()` and by `clean_consumed`."""
    mode = "hybrid-b1"
    world = World(tmp_path, mode)
    monkeypatch.setattr(bs, "_launch", world.launch)

    def stop(point: str) -> None:
        if point == "after:cluster.rep90":
            raise Killed(point)

    with pytest.raises(Killed):
        world.go(stop)
    ms, steps, _ = cluster._prune_spec("cluster.rep90", mode)
    want = json.dumps({"entries": [{"step": "cluster.rep90", "milestones": ms, "steps": steps}]}, sort_keys=True,
                      indent=1).encode() + b"\n"
    assert (world.w / ".pruned.json").read_bytes() == want
    # a B-72 entry, written by hand over the same scratch, for a prune whose targets are gone
    (world.w / ".pruned.json").write_bytes(want)
    assert cluster.lost(world.root, REL, world.scratch) == []
    assert cluster.clean_consumed(world.root, REL, "cluster.rep90", world.scratch) == []


# --- (b) a redo after a lost milestone prunes again, without MMseqs2 -------------------------------------------------


def _lost_milestone(world: World, stop_after: str, milestone: str) -> list[str]:
    def stop(point: str) -> None:
        if point == f"after:{stop_after}":
            raise Killed(point)

    with pytest.raises(Killed):
        world.go(stop)
    gone = [p for p in world.w.glob(milestone + "*") if p.is_file()]
    assert gone
    for p in gone:
        p.unlink()
    return cluster.lost(world.root, REL, world.scratch)


def test_r02_a_redo_after_a_lost_milestone_prunes_again(tmp_path, monkeypatch):
    """Pruned at `cluster.rep90` (hybrid-b1: `tmp90/w`), then its milestone `clu90` is lost after `cluster.merge`:
    the resume makes the NR90 steps again (which write `tmp90/w` again) and prunes `tmp90/w` again; the scratch before
    the tables is an uninterrupted clustering's."""
    mode = "hybrid-b1"
    ref = World(tmp_path / "ref", mode)
    monkeypatch.setattr(bs, "_launch", ref.launch)
    assert ref.go() is None

    world = World(tmp_path / "got", mode)
    monkeypatch.setattr(bs, "_launch", world.launch)
    names = _lost_milestone(world, "cluster.merge", "clu90")
    nr90 = [c.name for c in cluster.calls(mode) if c.name.startswith("cluster.nr90.") and c.name != "cluster.nr90.members"]
    assert names == nr90
    assert world.go() is None
    assert world.launched == [f"{n}.redo" for n in nr90]
    assert not any(k.startswith("tmp90/") for k in world.final), "a redo's intermediates were not pruned again"
    assert world.final == ref.final
    assert world.records() == ref.records()


# --- (c) PLAN_FORMAT 1 cut records, without MMseqs2 ------------------------------------------------------------------


def _cut_doc(tmp_path: Path, fmt: int, units: int, keys: int, plan_cut: dict) -> dict:
    c = {c.name: c for c in cluster.calls("cascade")}["cluster.nr50.15.prefilter"]
    d = tmp_path / f"f{fmt}"
    src = d / f"{c.positional[cluster.PER_QUERY[c.module]]}.index"
    src.parent.mkdir(parents=True, exist_ok=True)
    src.write_text("".join(f"{k}\t{k * 10}\t{units // keys + 2}\n" for k in range(keys)))
    plan = {"format": fmt, "proteins": 10_620_684, "cut": plan_cut} | ({} if fmt == 1 else {"method": "cascade"})
    ctx = cluster._Ctx(root=d, out=None, release="rX", plan=plan, plan_sha256="", w=d, threads=1)
    return cluster._cut_doc(ctx, c)


def test_r02_a_format_1_cut_record_has_no_createindex_estimate(tmp_path):
    """A shared-index cut into two parts: under a format-1 plan its record holds what the code before card B-71 wrote
    (`index_cpu_min`, no `createindex_cpu_min` and no DEC-77 `one_process`); under format 2, as card B-71 made it."""
    cut1 = _cut_doc(tmp_path, 1, 500_000_000, 2_000, {"parts": 2, "one_process_cpu_min_max": 20_000.0})
    assert set(cut1) == {"split", "units", "cpu_min", "overhead_cpu_min", "parts", "index_cpu_min"}
    cut2 = _cut_doc(tmp_path, 2, 500_000_000, 2_000, {"parts": 2, "one_process_cpu_min_max": 20_000.0})
    assert set(cut2) == set(cut1) | {"createindex_cpu_min"}
    assert {k: v for k, v in cut2.items() if k != "createindex_cpu_min"} == cut1
    # over the part budget by the DEC-77 estimate: format 2 runs it under DEC-77, format 1 has no such exemption
    big = {"part_cpu_min": 240.0, "one_process_cpu_min_max": 20_000.0}
    assert _cut_doc(tmp_path / "big", 2, 4_000_000_000_000, 1_000, big)["one_process"]["decision"] == "DEC-77"
    assert "one_process" not in _cut_doc(tmp_path / "big", 1, 4_000_000_000_000, 1_000, big)


def test_r02_a_format_1_plan_gets_no_dec77_createindex_exemption_in_build_stage(tmp_path):
    """`cut_createindex_estimate` (build-stage's DEC-77 check before a shared-index cut) is None for a format-1 plan."""
    root, scratch = tmp_path / "root", tmp_path / "scratch"
    c = {c.name: c for c in cluster.calls("cascade")}["cluster.nr50.15.prefilter"]
    for fmt in (1, 2):
        rel = f"r{fmt}"
        doc = {"format": fmt, "release": rel, "stages": ["S0"], "gene_call_version": V, "mmseqs2": "fake",
               "files": [], "cut": {"part_cpu_min": 240.0}} | ({} if fmt == 1 else {"method": "cascade"})
        cluster._write_record(root, root, rel, cluster.PLAN_STEP, doc)
        idx = cluster.work_dir(scratch, root, rel) / f"{c.positional[1]}.index"
        idx.parent.mkdir(parents=True)
        idx.write_text("".join(f"{k}\t0\t{4_000_000_000}\n" for k in range(1_000)))
    assert cluster.cut_createindex_estimate(root, "r2", c, scratch) > 240.0
    assert cluster.cut_createindex_estimate(root, "r1", c, scratch) is None


# --- (d) the method in the ledger, without MMseqs2 --------------------------------------------------------------------


def test_r02_the_clustering_ledger_records_a_method_other_than_cascade(monkeypatch):
    tools = {"mmseqs2": "fake", "pyarrow": "fake"}
    monkeypatch.setattr(cluster.provenance, "tool_versions", lambda *names: {n: tools[n] for n in names})
    counts = {"nr90": 70, "nr50": 20}
    # cascade: the two rows of before card B-71, as `_publish_prepared` built them then
    old = sorted((ledger.row(lvl.encode(), "level", n_in, counts[lvl], n_in - counts[lvl], tools=tools)
                  for lvl, n_in in (("nr90", 100), ("nr50", 70))), key=lambda r: (r["object_key"], r["reason_code"]))
    assert cluster._cluster_ledger_rows(100, counts, "cascade") == old
    assert cluster.ledger_method(old) == "cascade"
    for mode in ("hybrid-b1", "hybrid-a", "linclust"):
        rows = cluster._cluster_ledger_rows(100, counts, mode)
        assert [r for r in rows if r["reason_code"] == "ok"] == old
        meth = [r for r in rows if r["reason_code"] != "ok"]
        assert [r["reason_code"] for r in meth] == [f"method:{mode}"] * 2
        assert [{**r, "reason_code": "ok"} for r in meth] == old  # same counts and tools: L-01 and L-03 hold
        assert cluster.ledger_method(rows) == mode
    with pytest.raises(ValueError, match="2 clustering methods"):
        cluster.ledger_method(cluster.method_rows(old, "hybrid-a") + cluster.method_rows(old, "hybrid-b1"))


# --- with MMseqs2: the real steps ---------------------------------------------------------------------------------


def _synth_root(root: Path, bases: int = 25, seed: int = 202) -> Path:
    """Synthetic proteins with NR90 and NR50 families (as tests/r1/test_cluster_pruning.py makes them)."""
    rng = random.Random(seed)
    aa = "ACDEFGHIKLMNPQRSTVWY"

    def mutate(s: str, rate: float) -> str:
        return "".join(rng.choice(aa.replace(c, "")) if rng.random() < rate else c for c in s)

    seqs = set()
    for _ in range(bases):
        base = "M" + "".join(rng.choice(aa) for _ in range(rng.randint(120, 260)))
        for _ in range(4):
            var = mutate(base, 0.25)
            for _ in range(3):
                seqs.add(mutate(var, 0.04))
    rows = sorted((ids.protein_uid(s), s) for s in seqs)
    t = pa.table({"protein_uid": [u for u, _ in rows],
                  "protein_digest": [hashlib.sha256(s.encode()).digest() for _, s in rows],
                  "length": [len(s) for _, s in rows], "seq": [s for _, s in rows]})
    io.write_table(t, table("protein_seq"),
                   root / f"l0seq/protein_seq/gene_call_version={V}/stage=S0/protein_prefix=all/part-0.parquet")
    return root


@pytest.fixture(scope="module")
def proteins(tmp_path_factory):
    return _synth_root(tmp_path_factory.mktemp("r02-proteins") / "root")


@pytest.fixture(scope="module")
def oracle(proteins, tmp_path_factory):
    """The one-process clustering (`cluster.run`) of the proteins, by mode: the tables an interrupted clustering in
    parts must end with."""
    got = {}

    def of(mode: str) -> Path:
        if mode not in got:
            base = tmp_path_factory.mktemp(f"r02-oracle-{mode}")
            (base / "tmp").mkdir()
            cluster.run(proteins, REL, V, ["S0"], threads=1, tmp=base / "tmp", out=base / "out", mode=mode)
            got[mode] = base / "out"
        return got[mode]

    return of


@pytest.fixture(scope="module")
def reference(proteins, tmp_path_factory):
    """An uninterrupted clustering in parts by build-stage's path, by mode (made once)."""
    got = {}

    def of(mode: str) -> World:
        if mode not in got:
            ref = World(tmp_path_factory.mktemp(f"r02-ref-{mode}"), mode, real=True, proteins=proteins)
            assert ref.go() is None  # a first clustering: no plan published, so nothing is made again
            got[mode] = ref
        return got[mode]

    return of


REAL_KILLS = [("cluster.rep90", "last", "hybrid-b1"), ("cluster.nr50.14.createsubdb", "last", "hybrid-b1"),
              ("cluster.merge", "write_pruned", "hybrid-b1"), ("cluster.finish.nr90", "first", "hybrid-b1"),
              ("cluster.nr90.13.clust", "last", "cascade")]


@NEED_MMSEQS
@pytest.mark.parametrize("trigger,where,mode", REAL_KILLS)
def test_r02_mmseqs_a_kill_inside_clean_consumed_resumes_without_redo_byte_identical(
        tmp_path, monkeypatch, proteins, oracle, reference, trigger, where, mode):
    """The real steps: a kill inside `clean_consumed` and build-stage's resume make no step again, and the tables are
    the one-process clustering's (`cluster.run`), byte for byte."""
    ref = reference(mode)
    ref, got = _crash_and_resume(tmp_path, monkeypatch, trigger, where, mode, ref=ref, real=True, proteins=proteins)
    want = _tables(oracle(mode))
    assert want and _tables(ref.root) == want
    assert _tables(got.root) == want
    assert got.records() == ref.records()
    assert sorted(got.final) == sorted(ref.final)  # the same files in scratch before the tables (logs aside)
    assert got.final[".pruned.json"] == ref.final[".pruned.json"]


@NEED_MMSEQS
def test_r02_mmseqs_a_redo_after_a_lost_milestone_prunes_again_byte_identical(tmp_path, monkeypatch, proteins,
                                                                             oracle):
    mode = "hybrid-b1"
    world = World(tmp_path, mode, real=True, proteins=proteins)
    monkeypatch.setattr(bs, "_launch", world.launch)
    names = _lost_milestone(world, "cluster.merge", "clu90")
    assert "cluster.nr90.00.kmermatcher" in names and "cluster.nr90.09.mergeclusters" in names
    assert world.go() is None
    assert world.launched == [f"{n}.redo" for n in names]
    assert not any(k.startswith("tmp90/") for k in world.final), "a redo's intermediates were not pruned again"
    assert _tables(world.root) == _tables(oracle(mode))


@NEED_MMSEQS
def test_r02_mmseqs_a_format_1_work_directory_is_redone_to_its_records(tmp_path, proteins, oracle):
    """A format-1 clustering (cascade, every call sharing an index cut into two parts), its plan written as the code
    before card B-71 wrote it and every step run under it; its scratch is then lost (a reboot), and every step is made
    again to its record, after which the tables are the one-process clustering's."""
    root, tmp = tmp_path / "root", tmp_path / "scratch"
    shutil.copytree(proteins, root)
    inputs = cluster.plan_inputs(root, REL, V, ["S0"], cluster.cut_of(parts=2), mode="cascade")
    v1 = {k: inputs[k] for k in cluster._INPUT_KEYS_V1} | {"format": 1}
    w = cluster._work_dir(tmp, root, v1)
    w.mkdir(parents=True)
    n, outputs = cluster._make_db(root, w, v1)  # the plan step's work (`plan`), into the format-1 work directory
    from biobox.stage import cluster_cpu_min

    cluster._write_record(root, root, REL, cluster.PLAN_STEP,
                          {**v1, "proteins": n, "release_cpu_min": round(cluster_cpu_min(n), 6), "outputs": outputs})
    assert set(cluster.read_record(root, REL, cluster.PLAN_STEP)) == PRE_B71_PLAN_KEYS
    assert cluster.work_dir(tmp, root, REL) == w
    cluster.run_split(root, REL, V, ["S0"], tmp=tmp, threads=1, until=cluster.FINAL_STEP)

    recs = {p.name: p.read_bytes() for p in sorted((root / cluster.plan_dir(REL)).glob("*.json"))}
    assert {json.loads(b)["format"] for b in recs.values()} == {1}
    shared = [c.name for c in cluster.calls("cascade") if c.shares_index]
    cut = {s: cluster.read_record(root, REL, f"{s}.cut") for s in shared}
    split = [s for s in shared if len(cut[s]["parts"]) > 1]
    assert split, "no call sharing an index was cut into parts"
    for s in split:
        assert set(cut[s]) == PRE_B71_SHARED_CUT_KEYS, f"{s}.cut: {sorted(set(cut[s]) ^ PRE_B71_SHARED_CUT_KEYS)}"

    shutil.rmtree(w)  # the scratch is lost
    names = cluster.lost(root, REL, tmp)
    assert names[0] == cluster.PLAN_STEP and all(f"{s}.cut" in names for s in split)
    for name in names:
        assert cluster.redo(root, REL, name, tmp, 1)["verified"]
    assert cluster.lost(root, REL, tmp) == []
    assert {p.name: p.read_bytes() for p in sorted((root / cluster.plan_dir(REL)).glob("*.json"))} == recs
    cluster.finish(root, REL, tmp, 1)
    assert _tables(root) == _tables(oracle("cascade"))


def _ledger(base: Path, release: str = REL) -> list[dict]:
    return pq.read_table(Path(base) / f"ledger/release={release}/step=cluster/part-0.parquet").to_pylist()


def _fake_run_record(root: Path, release: str) -> None:
    """The run record of step `cluster` that run-step would write (as tests/r1/test_cluster_hybrid.py makes it)."""
    outputs = []
    for sub in (f"l2/cluster_member/release={release}", f"l2/cluster_member_aln/release={release}",
                f"l2/member_cluster/release={release}", f"ledger/release={release}/step=cluster"):
        for p in sorted((root / sub).glob("*.parquet")):
            outputs.append({"path": p.relative_to(root).as_posix(), "bytes": p.stat().st_size,
                            "sha256": ids.sha256_file(p).hex()})
    rec = root / f"runs/stage={release}/cluster.json"
    rec.parent.mkdir(parents=True, exist_ok=True)
    rec.write_text(json.dumps({"stage": release, "step": "cluster", "exit_code": 0, "outputs": outputs}))


@NEED_MMSEQS
def test_r02_mmseqs_the_method_is_in_the_ledger_and_the_plan(tmp_path, proteins, oracle):
    """hybrid-b1, in one process and in parts: the clustering ledger names the method (and L-01 holds), the plan record
    names it, and a reuse of that clustering keeps it; a cascade ledger has only its two rows of before card B-71."""
    one = _ledger(oracle("hybrid-b1"))
    assert cluster.ledger_method(one) == "hybrid-b1"
    assert sorted(r["reason_code"] for r in one) == ["method:hybrid-b1", "method:hybrid-b1", "ok", "ok"]
    assert ledger.check_balance(oracle("hybrid-b1"), "S0", ["cluster"], REL) == []
    cas = _ledger(oracle("cascade"))
    assert [r["reason_code"] for r in cas] == ["ok", "ok"] and cluster.ledger_method(cas) == "cascade"

    root = tmp_path / "root"
    shutil.copytree(proteins, root)
    cluster.run_split(root, "r1", V, ["S0"], tmp=tmp_path / "tmp1", parts=1, threads=1, mode="hybrid-b1",
                      prune=True, split_finish=True)
    assert cluster.read_record(root, "r1", cluster.PLAN_STEP)["method"] == "hybrid-b1"
    assert _ledger(root, "r1") == one  # the same rows in parts as in one process
    _fake_run_record(root, "r1")
    assert cluster.reuse(root, "r2", V, ["S0"], tmp_path / "tmp2", parts=1, mode="hybrid-b1")["reused_from"] == "r1"
    reused = _ledger(root, "r2")
    assert cluster.ledger_method(reused) == "hybrid-b1"
    assert [r for r in reused if not r["reason_code"].startswith("reused:")] == one
    assert ledger.check_balance(root, "S0", ["cluster"], "r2") == []
