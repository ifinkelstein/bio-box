"""Ring R2: gene calling in batch steps under the CPU budget, byte for byte the one-process call (card B-48).

The oracle of every byte-identity test is the one-process call (`genecall.run`, the step `build-stage` runs today):
its files by sha256, its ledger, and its result. The world, built from the micro-S0 genomes as new FASTA objects:

    S0 (stored first): G0 = (M. genitalium M circular, plasmid P = T4[:8000]) prokaryote code 4, single mode
                       GL = (phage lambda L) virus, gv-meta
    S1 (the stage gene-called in batches; 14 genomes):
        GA = (reverse complement of M circular, P)   prokaryote code 4, single mode on two contigs; P is S0's (DEC-26)
        GB = (L)                                     virus; its one contig is S0's: not called at all (item 5)
        GC = (T4)                                    virus, gv-meta, 169 kb
        GD = (Q circular) plasmid, meta mode with the wrap flanks (Q = T4[8000:38000], whole contig as its flank)
        GE = (M[100000:190000]) prokaryote code 4 under 100 kb: meta; its proteins are also GA's (the ledger's
             first-of-protein window spans batches)
        GF = (X, Y), GG = (X) a contig X in two genomes of S1 (DEC-27; both prokaryote, code 11, linear, so they
             agree on its call inputs): the least genome_uid writes it
        H0..H4 = five 8 kb pieces of T4, viruses under THREAD_MIN_BASES (worker processes); their proteins are GC's
        K0, K1 = two 20 kb pieces of M, prokaryote code 4 (meta)

Batches are run as the commands build-stage runs (`python -m biobox.cli stage genecall batch|merge`), each its own
process; `build-stage` is run with its steps as real run-step subprocesses, and with the gene-calling rates raised in
this process (where the plan is made) so that S1 is over `SHARD_MARGIN` of the budget: the batches follow the plan's
file, so the rates do not reach them. The allowance of the hour is fixed (`IDLE`): the pod's real reading counts
other agents' processes, which would hold a launch.
"""

from __future__ import annotations

import hashlib
import json
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest

from biobox import catalog, fasta, ids, ledger, runstep
from biobox import stage as bs
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l0seq.store import revcomp
from biobox.l1 import genecall
from biobox.l1.genecall import GENE_CALL_VERSION as V

sys.path.insert(0, str(Path(__file__).parent))
import test_builders_out as c08  # noqa: E402

pytestmark = [pytest.mark.r2, pytest.mark.covers("B-48", "D-01")]

FIX = c08.FIX
NOON_MONDAY = c08.NOON_MONDAY
IDLE = runstep.CpuUse(window_s=5.0)  # card C-28: nothing else computes


def _seq(acc: str) -> bytes:
    return fasta.read(FIX / f"{acc}.fna.gz")[0].seq


M, L, T4 = _seq("NC_000908.2"), _seq("NC_001416.1"), _seq("NC_000866.4")
SEQS = {"M": M, "L": L, "T4": T4, "P": T4[:8000], "M_RC": revcomp(M), "Q": T4[8000:38000], "E": M[100_000:190_000],
        "X": T4[60_000:66_000], "Y": M[300_000:330_000],
        **{f"H{i}": T4[70_000 + 15_000 * i: 78_000 + 15_000 * i] for i in range(5)},
        "K0": M[400_000:420_000], "K1": M[450_000:470_000]}
CU = {n: ids.contig_uid(s) for n, s in SEQS.items()}
# per stage: objects (file, [(header, contig name)]) and genomes (contig names, code, kind, circular names)
STAGES = {
    "S0": {"objects": [("b48_s0_host.fna", [("hostA chromosome", "M"), ("hostA plasmid", "P")]),
                       ("b48_s0_lambda.fna", [("phage lambda", "L")])],
           "genomes": [(("M", "P"), 4, "prokaryote", ("M",)), (("L",), 11, "virus", ())]},
    "S1": {"objects": [("b48_s1_a.fna", [("hostB chromosome", "M_RC"), ("hostB plasmid", "P"), ("lambda again", "L"),
                                         ("T4", "T4"), ("plasmid Q", "Q"), ("fragment E", "E")]),
                       ("b48_s1_b.fna", [("X", "X"), ("Y", "Y"), *((f"piece {i}", f"H{i}") for i in range(5)),
                                         ("K0", "K0"), ("K1", "K1")])],
           "genomes": [(("M_RC", "P"), 4, "prokaryote", ("M_RC",)), (("L",), 11, "virus", ()),
                       (("T4",), 11, "virus", ()), (("Q",), 11, "plasmid", ("Q",)), (("E",), 4, "prokaryote", ()),
                       (("X", "Y"), 11, "prokaryote", ()), (("X",), 11, "prokaryote", ()),
                       *(((f"H{i}",), 11, "virus", ()) for i in range(5)),
                       (("K0",), 4, "prokaryote", ()), (("K1",), 4, "prokaryote", ())]},
}
SEQ_STAGES = ["S1", "S0"]
BUDGET = 300  # the profile's cpu_min_per_process (c08.profile_at)


def genome_inputs(stage: str) -> list[GenomeInput]:
    return [GenomeInput(tuple(CU[n] for n in names), code, kind=kind, circular=tuple(CU[n] for n in circ))
            for names, code, kind, circ in STAGES[stage]["genomes"]]


def write_inputs(base: Path) -> dict[str, tuple[Path, Path]]:
    """The FASTA objects and the CLI's genomes and objects JSON files of each stage."""
    base.mkdir(parents=True, exist_ok=True)
    out = {}
    for stage, spec in STAGES.items():
        objs = []
        for name, recs in spec["objects"]:
            lines = [ln for h, n in recs for ln in fasta.format_record(h, SEQS[n], 80)]
            (base / name).write_bytes(b"\n".join(lines) + b"\n")
            objs.append({"path": str(base / name), "source_name": "b48_fixture", "source_release": "1",
                         "license": "INSDC (public)", "redistributable": True})
        g = base / f"genomes-{stage}.json"
        g.write_text(json.dumps([{"contig_uids": [u.hex() for u in x.contig_uids], "genetic_code": x.genetic_code,
                                  "kind": x.kind, "circular": [u.hex() for u in x.circular],
                                  "assembly_quality": x.assembly_quality} for x in genome_inputs(stage)]))
        o = base / f"objects-{stage}.json"
        o.write_text(json.dumps(objs))
        out[stage] = (g, o)
    return out


def source_objects(base: Path, stage: str) -> list[SourceObject]:
    return [SourceObject(base / name, "b48_fixture", "1", "INSDC (public)", True) for name, _ in STAGES[stage]["objects"]]


def tree(root: Path) -> dict[str, str]:
    """sha256 of every file under `root`, except run-step's `runs/`."""
    root = Path(root)
    if not root.exists():
        return {}
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file() and p.relative_to(root).parts[0] != "runs"}


GENECALL_DIRS = (f"l1/feature/producer={V}/stage=S1", f"l1/feature_protein/gene_call_version={V}/stage=S1",
                 f"l1/protein_occurrence/gene_call_version={V}/stage=S1", f"l0seq/protein_seq/gene_call_version={V}/stage=S1",
                 "ledger/stage=S1/step=genecall")


def genecall_files(t: dict[str, str]) -> dict[str, str]:
    return {k: v for k, v in t.items() if k.startswith(GENECALL_DIRS)}


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """S0 and S1 ingested with their genome inputs (no gene calls), and the one-process call of S1 into `one`."""
    d = tmp_path_factory.mktemp("b48")
    inputs = write_inputs(d / "in")
    root = d / "root"
    for s in ("S0", "S1"):
        write_genome_input(root, s, genome_inputs(s))
        ingest(source_objects(d / "in", s), root, s)
    (d / "scratch").mkdir()
    result = genecall.run(root, "S1", seq_stages=SEQ_STAGES, out=d / "one", tmp=d / "scratch")
    return {"d": d, "root": root, "inputs": inputs, "result": result, "one": tree(d / "one")}


def cli(*args: str, check: bool = True) -> subprocess.CompletedProcess:
    p = subprocess.run([sys.executable, "-m", "biobox.cli", *args], capture_output=True, text=True)
    if check:
        assert p.returncode == 0, p.stderr[-4000:]
    return p


def batch_args(root: Path, plan: Path, j: int, workers: int = 1) -> list[str]:
    return ["stage", "genecall", "batch", "--root", str(root), "--stage", "S1", "--plan", str(plan), "--batch", str(j),
            "--workers", str(workers)]


def merge_args(root: Path, plan: Path, out: Path, tmp: Path) -> list[str]:
    return ["stage", "genecall", "merge", "--root", str(root), "--stage", "S1", "--plan", str(plan), "--out", str(out),
            "--tmp", str(tmp)]


def make_plan(world, tmp_path: Path, n: int) -> Path:
    plan = genecall.plan(world["root"], "S1", seq_stages=SEQ_STAGES, batch_cpu_min=bs.SHARD_MARGIN * BUDGET, batches=n,
                         tmp=tmp_path)
    return genecall.write_plan(plan, genecall.work_dir(tmp_path / "work", world["root"], "S1", plan))


def run_batches(world, plan: Path, n: int, jobs: int, workers: int = 1) -> None:
    """The plan's batches, `jobs` processes at once, each `stage genecall batch`."""
    pending = list(range(n))
    while pending:
        procs = [subprocess.Popen([sys.executable, "-m", "biobox.cli", *batch_args(world["root"], plan, j, workers)],
                                  stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
                 for j in pending[:jobs]]
        for p in procs:
            _, err = p.communicate()
            assert p.returncode == 0, err[-4000:]
        pending = pending[jobs:]


# --- the batches and the merge equal the one-process call ------------------------------------------------------------


@pytest.mark.parametrize("n,jobs,workers", [(1, 1, 1), (2, 1, 1), (3, 3, 1), (3, 1, 3), (5, 1, 1), (5, 3, 2)])
def test_b48_batches_and_the_merge_are_the_one_process_call_byte_for_byte(world, tmp_path, n, jobs, workers):
    """Every gene-call file and the ledger have the one-process call's sha256, and the merge prints its result, at 1,
    2, 3 and 5 batches, 1 and 3 batches at once, 1 to 3 threads each; nothing is written into the root."""
    before = tree(world["root"])
    plan = make_plan(world, tmp_path, n)
    assert len(genecall.read_plan(plan)["batches"]) == n
    run_batches(world, plan, n, jobs, workers)
    p = cli(*merge_args(world["root"], plan, tmp_path / "out", tmp_path))
    assert json.loads(p.stdout)["result"] == json.loads(json.dumps(world["result"]))
    assert tree(tmp_path / "out") == world["one"]
    assert len(genecall_files(world["one"])) == len(world["one"]) == 5  # feature, protein tables, ledger
    assert tree(world["root"]) == before


def test_b48_the_plan_holds_every_genome_once_in_uid_order_cut_by_cost(world, tmp_path):
    """The plan's batches are consecutive runs of the stage's genomes in genome_uid order (a brute-force loop over the
    fixture's genome list and FASTA lengths): every genome once, the bases of each batch theirs, and each batch's
    estimate under the budget given; a genome whose every contig another stage or genome owns is counted and costs
    nothing."""
    guid = {g: ids.genome_uid(list(g.contig_uids)) for g in genome_inputs("S1")}
    gs = sorted(guid, key=guid.get)
    length = {CU[n]: len(s) for n, s in SEQS.items()}
    # the rule, written out: S0 stored P and L; X is in GF and GG, and the least genome_uid of the two writes it
    owner_of_x = min(guid[g] for g in gs if CU["X"] in g.contig_uids)
    called = {g: any(c not in (CU["P"], CU["L"]) and (c != CU["X"] or guid[g] == owner_of_x) for c in g.contig_uids)
              for g in gs}
    for n in (1, 2, 3, 5, 14):
        plan = genecall.plan(world["root"], "S1", seq_stages=SEQ_STAGES, batch_cpu_min=240.0, batches=n, tmp=tmp_path)
        bat = plan["batches"]
        assert [b["first"] for b in bat] == [sum(x["genomes"] for x in bat[:k]) for k in range(n)]
        assert sum(b["genomes"] for b in bat) == plan["genomes"] == len(gs) == 14
        for b in bat:
            mine = gs[b["first"]: b["first"] + b["genomes"]]
            assert b["bases"] == sum(length[c] for g in mine if called[g] for c in g.contig_uids)
            assert b["estimate_cpu_min"] <= 240.0
        assert plan["bases"] == sum(length[c] for g in gs for c in g.contig_uids)
        assert plan["bases_called"] == sum(b["bases"] for b in bat)
        # GB always; GG = (X) too when GF writes X
        assert plan["genomes_not_called"] == sum(not c for c in called.values()) == 1 + (
            owner_of_x != ids.genome_uid([CU["X"]]))
    auto = genecall.plan(world["root"], "S1", seq_stages=SEQ_STAGES, batch_cpu_min=240.0, tmp=tmp_path)
    assert len(auto["batches"]) == 1  # the fixture is far under the budget


def test_b48_the_merge_refuses_until_every_batch_is_there_and_a_batch_refuses_another_plan(world, tmp_path):
    plan = make_plan(world, tmp_path, 3)
    run_batches(world, plan, 2, 1)  # batches 0 and 1 only
    p = cli(*merge_args(world["root"], plan, tmp_path / "out", tmp_path), check=False)
    assert p.returncode == 2 and "batch 2" in p.stderr, p.stderr
    assert not (tmp_path / "out").exists()
    other = json.loads(plan.read_text())
    other["batches"][0]["digest"] = "0" * 32  # a plan made from other genome inputs
    bad = plan.parent.parent / "other" / "plan.json"
    bad.parent.mkdir()
    bad.write_text(json.dumps(other))
    p = cli(*batch_args(world["root"], bad, 0), check=False)
    assert p.returncode == 2 and "digest" in p.stderr, p.stderr
    assert not any(x.name.startswith(("batch-", ".batch-")) for x in bad.parent.iterdir())


def test_b48_a_batch_with_output_is_not_run_again(world, tmp_path):
    plan = make_plan(world, tmp_path, 2)
    run_batches(world, plan, 2, 1)
    done = (plan.parent / "batch-00000" / "done.json").read_bytes()
    p = cli(*batch_args(world["root"], plan, 0))
    assert "not run again" in p.stderr
    assert (plan.parent / "batch-00000" / "done.json").read_bytes() == done


def test_b48_plan_action_prints_the_estimate_and_writes_nothing_in_the_root(world, tmp_path):
    before = tree(world["root"])
    p = cli("stage", "genecall", "plan", "--root", str(world["root"]), "--stage", "S1", "--seq-stages", "S1,S0",
            "--runs-tmp", str(tmp_path / "w"), "--batch-cpu-min", "240", "--batches", "3")
    res = json.loads(p.stdout)["result"]
    assert res["batches"] == 3 and res["genomes"] == 14 and Path(res["plan"]).is_file()
    assert Path(res["plan"]).is_relative_to(tmp_path / "w")
    assert res["estimate_cpu_min"] > 0 and len(res["batch_estimates_cpu_min"]) == 3
    p = cli("stage", "genecall", "plan", "--root", str(world["root"]), "--stage", "S1", "--runs-tmp",
            str(world["root"] / "w"), "--batch-cpu-min", "240", check=False)
    assert p.returncode == 2 and "under the data root" in p.stderr, p.stderr
    assert tree(world["root"]) == before


# --- item 5: a genome all of whose contigs another stage owns is not called ---------------------------------------


def test_b48_a_genome_whose_contigs_another_stage_owns_is_not_called_and_its_skip_is_recorded(world, tmp_path,
                                                                                            monkeypatch):
    """GB = (lambda) of S1: S0 stored lambda. It is never called (counted through `_call_batch` and `_call_read`, the
    in-process and thread paths), its contig is skipped as before (`owned_by:S0`, counts 0), and the stage's files are
    those of the same stage without GB (an independent construction: S1's genome_input without that genome), but for
    that one ledger row and its entry in `skipped_contigs`."""
    called: list[bytes] = []
    real_batch, real_read = genecall._call_batch, genecall._call_read

    def batch(store, gs):
        called.extend(ids.genome_uid(list(g.contig_uids)) for g in gs)
        return real_batch(store, gs)

    def read(g, seqs):
        called.append(ids.genome_uid(list(g.contig_uids)))
        return real_read(g, seqs)

    monkeypatch.setattr(genecall, "_call_batch", batch)
    monkeypatch.setattr(genecall, "_call_read", read)
    res = genecall.run(world["root"], "S1", seq_stages=SEQ_STAGES, out=tmp_path / "a", tmp=tmp_path)
    assert ids.genome_uid([CU["L"]]) not in called and len(called) >= 12
    assert res == world["result"]
    rows = [r for r in ledger.read(tmp_path / "a", "S1", "genecall") if r["object_key"] == CU["L"]]
    assert [(r["reason_code"], r["n_in"], r["n_kept"]) for r in rows] == [("owned_by:S0", 0, 0)]
    # the same stage without GB, in a root of its own
    root2 = tmp_path / "root2"
    d = world["d"]
    write_genome_input(root2, "S0", genome_inputs("S0"))
    ingest(source_objects(d / "in", "S0"), root2, "S0")
    write_genome_input(root2, "S1", [g for g in genome_inputs("S1") if g.contig_uids != (CU["L"],)])
    ingest(source_objects(d / "in", "S1"), root2, "S1")
    res2 = genecall.run(root2, "S1", seq_stages=SEQ_STAGES, out=tmp_path / "b", tmp=tmp_path)
    a, b = tree(tmp_path / "a"), tree(tmp_path / "b")
    assert {k: v for k, v in a.items() if "ledger" not in k} == {k: v for k, v in b.items() if "ledger" not in k}
    la = [r for r in ledger.read(tmp_path / "a", "S1", "genecall") if r["object_key"] != CU["L"]]
    assert la == ledger.read(tmp_path / "b", "S1", "genecall")
    assert res["skipped_contigs"] == sorted([*res2["skipped_contigs"], {"contig_uid": CU["L"].hex(),
                                                                       "owner_stage": "S0"}],
                                            key=lambda e: e["contig_uid"])


# --- build-stage ------------------------------------------------------------------------------------------------------

REL = "rB48.1"


class Killed(Exception):
    """Raised by the injected hook: the build stops after gene calling (no MMseqs2 needed)."""


def stop_after_genecall(seen: list[str] | None = None):
    def hook(point: str) -> None:
        if seen is not None:
            seen.append(point)
        if point == "after:genecall":
            raise Killed(point)
    return hook


def runner(prof, root, stage, step, cmd, **kw) -> int:
    return runstep.run_step(prof, root, stage, step, cmd, clock=lambda: NOON_MONDAY, load=lambda: 0.0, **kw)


def s0_root(world, root: Path) -> Path:
    """A root holding S0 ingested with its genome inputs, as S1's build finds it."""
    write_genome_input(root, "S0", genome_inputs("S0"))
    ingest(source_objects(world["d"] / "in", "S0"), root, "S0")
    return root


def build_of(world, root: Path, tmp: Path, **kw) -> bs.Build:
    g, o = world["inputs"]["S1"]
    return bs.Build(root=root, stage="S1", release=REL, stages=("S0", "S1"), genomes=g, objects=o, threads=1, tmp=tmp,
                    **kw)


def build(prof, b: bs.Build, hook, **kw) -> int | None:
    try:
        return bs.build_stage(prof, b, run=runner, clock=lambda: NOON_MONDAY, load=lambda: 0.0, cpu=lambda: IDLE,
                              _hook=hook, **kw)
    except Killed:
        return None


@pytest.fixture
def prof(tmp_path):
    return c08.profile_at(tmp_path / "prof")


def scale_rates(monkeypatch, k: float) -> None:
    """Every gene-calling rate times `k`, in this process (where build-stage makes the plan)."""
    monkeypatch.setattr(genecall, "CALL_S_PER_MBP", {m: r * k for m, r in genecall.CALL_S_PER_MBP.items()})


def rate_factor(largest: float) -> float:
    """The factor on every rate that makes S1's largest genome alone `largest` of the margin (its cost from the
    fixture's own contig lengths and kinds, `genecall.genome_cost_us`): the stage is then over the margin in one
    process, and fits in batches when `largest` < 1."""
    big = max(genecall.genome_cost_us(kind, [(len(SEQS[n]), n in circ) for n in names])
              for names, _, kind, circ in STAGES["S1"]["genomes"]) / 60e6
    return largest * bs.SHARD_MARGIN * BUDGET / big


def over_budget(world, monkeypatch, tmp_path, largest: float = 0.9) -> float:
    k = rate_factor(largest)
    scale_rates(monkeypatch, k)
    return k


@pytest.fixture(scope="module")
def one_process_build(world, tmp_path_factory):
    """S1 built by build-stage at the code's rates (one process, as before card B-48), stopped after gene calling."""
    d = tmp_path_factory.mktemp("b48-one")
    prof = c08.profile_at(d / "prof")
    root = s0_root(world, d / "root")
    assert build(prof, build_of(world, root, d / "scratch"), stop_after_genecall()) is None
    return root


def record(root: Path, name: str = "genecall") -> dict:
    return json.loads((root / f"runs/stage=S1/{name}.json").read_text())


def test_b48_within_the_margin_build_stage_calls_genes_in_one_process_as_before(world, one_process_build):
    """At the code's rates the fixture is far under the margin: the step's argv is the step list's (`stage genecall
    --root ... --workers N --seq-stages ...`, no action), its files are the one-process call's, and no batch ran."""
    rec = record(one_process_build)
    argv = rec["argv"][rec["argv"].index("biobox.cli") + 1:]
    assert argv[:3] == ["stage", "genecall", "--root"] and "--workers" in argv and "--seq-stages" in argv
    assert argv[argv.index("--seq-stages") + 1] == "S1,S0"
    assert genecall_files(tree(one_process_build)) == world["one"]
    assert sorted(p.name for p in (one_process_build / "runs/stage=S1").iterdir()) == [
        "genecall.json", "ingest.json", "inputs.json"]


@pytest.mark.parametrize("jobs", [1, 3])
def test_b48_over_the_margin_build_stage_calls_genes_in_batches_and_merges_the_one_process_files(
        world, one_process_build, prof, tmp_path, monkeypatch, jobs):
    """With the rates raised so that S1 is over the margin, build-stage plans batches, launches at most `jobs` at once
    (each its own process, `--scan-cpus` threads), and the merge is the run-step `genecall`: its record lists exactly
    the one-process record's outputs, with their sha256, and the root holds the same files; nothing in `runs/` names a
    batch, and the plan's work directory is removed."""
    over_budget(world, monkeypatch, tmp_path)
    root = s0_root(world, tmp_path / "root")
    live, peak, names = [0], [0], []
    lock = threading.Lock()

    def batch_run(prof_, root_, stage, step, cmd, *, threads):
        with lock:
            live[0] += 1
            peak[0] = max(peak[0], live[0])
            names.append(step)
        try:
            return bs._batch_process(prof_, root_, stage, step, cmd, threads=threads)
        finally:
            with lock:
                live[0] -= 1

    seen: list[str] = []
    b = build_of(world, root, tmp_path / "scratch", scan_jobs=jobs, scan_cpus=1)
    assert build(prof, b, stop_after_genecall(seen), batch_run=batch_run) is None
    n = len(names)
    assert n >= 2 and peak[0] <= jobs, (n, peak)
    assert sorted(names) == [f"genecall.batch-{j:05d}" for j in range(n)]
    assert [p for p in seen if p.startswith("after:genecall")][-1] == "after:genecall"
    assert sorted(p for p in seen if p.startswith("after:genecall.batch-")) == [f"after:{x}" for x in sorted(names)]
    rec, one = record(root), record(one_process_build)
    assert rec["outputs"] == one["outputs"]
    assert rec["argv"][rec["argv"].index("genecall") + 1] == "merge"
    assert genecall_files(tree(root)) == genecall_files(tree(one_process_build))
    assert tree(root) == tree(one_process_build)
    assert sorted(p.name for p in (root / "runs/stage=S1").iterdir()) == ["genecall.json", "ingest.json", "inputs.json"]
    assert not list((tmp_path / "scratch").glob("genecall*/**/plan.json"))


def test_b48_stop_between_batches_then_the_next_call_runs_only_the_batches_without_output(
        world, one_process_build, prof, tmp_path, monkeypatch):
    over_budget(world, monkeypatch, tmp_path)
    root = s0_root(world, tmp_path / "root")
    b = build_of(world, root, tmp_path / "scratch", scan_jobs=1, scan_cpus=1)

    def stop(point: str) -> None:
        if point == "after:genecall.batch-00000":
            (root / "STOP").write_text("")
        stop_after_genecall()(point)

    assert build(prof, b, stop) == runstep.EXIT_STOP
    assert not (root / "runs/stage=S1/genecall.json").exists() and not genecall_files(tree(root))
    plans = list((tmp_path / "scratch").glob("genecall/*/*/plan.json"))
    assert len(plans) == 1 and sorted(p.name for p in plans[0].parent.iterdir() if p.name.startswith("batch-")) == [
        "batch-00000"]
    (root / "STOP").unlink()
    again: list[str] = []

    def batch_run(prof_, root_, stage, step, cmd, *, threads):
        again.append(step)
        return bs._batch_process(prof_, root_, stage, step, cmd, threads=threads)

    assert build(prof, b, stop_after_genecall(), batch_run=batch_run) is None
    assert "genecall.batch-00000" not in again and again
    assert tree(root) == tree(one_process_build)


_HELD = r"""
import signal, sys
from pathlib import Path
from biobox import cli
from biobox.l1 import genecall
flag, point = sys.argv[-2:]
log = genecall._log


def held(msg):  # the batch stops here until a signal comes
    log(msg)
    if point in msg:
        Path(flag).write_text(msg)
        signal.pause()


genecall._log = held
sys.exit(cli.main(sys.argv[1:-2]))
"""


def _poll(pred, timeout: float = 120.0) -> bool:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if pred():
            return True
        time.sleep(0.05)
    return False


@pytest.mark.parametrize("sig,point", [(signal.SIGTERM, "calling"), (signal.SIGKILL, "calling"),
                                       (signal.SIGTERM, "writing")])
def test_b48_a_killed_batch_leaves_nothing_a_reader_takes_and_the_next_call_gives_the_same_bytes(
        world, one_process_build, prof, tmp_path, monkeypatch, sig, point):
    """Batch 1 is held at a point of its work, then sent a signal: build-stage returns the batch's code and launches
    no more; the batch's directory does not exist (by SIGTERM its temporary directory is removed too; by SIGKILL it is
    left, dot-named) and the root has no gene-call file. The next call removes what SIGKILL left, runs that batch
    again and the merge, and the root holds the one-process files."""
    over_budget(world, monkeypatch, tmp_path)
    root = s0_root(world, tmp_path / "root")
    b = build_of(world, root, tmp_path / "scratch", scan_jobs=1, scan_cpus=1)
    flag, work = tmp_path / "held", []

    def batch_run(prof_, root_, stage, step, cmd, *, threads):
        if step != "genecall.batch-00001":
            return bs._batch_process(prof_, root_, stage, step, cmd, threads=threads)
        plan = Path(cmd[cmd.index("--plan") + 1])
        work.append(plan.parent)
        p = subprocess.Popen([sys.executable, "-c", _HELD, *cmd[3:], str(flag), point])
        try:
            assert _poll(flag.exists), f"{step} never reached {point}"
            assert [x for x in plan.parent.iterdir() if x.name.startswith(".batch-00001.tmp-")]
            p.send_signal(sig)
            code = p.wait(timeout=120)
        finally:
            if p.poll() is None:
                p.kill()
        return code if code >= 0 else 128 - code

    code = build(prof, b, stop_after_genecall(), batch_run=batch_run)
    assert code == (143 if sig == signal.SIGTERM else 137), code
    wd = work[0]
    assert not (wd / "batch-00001").exists()
    left = [x.name for x in wd.iterdir() if x.name.startswith(".batch-00001.tmp-")]
    assert (left != []) == (sig == signal.SIGKILL), left
    assert not genecall_files(tree(root)) and not (root / "runs/stage=S1/genecall.json").exists()
    assert build(prof, b, stop_after_genecall()) is None
    assert tree(root) == tree(one_process_build)
    assert not wd.exists()


def test_b48_a_genome_over_the_margin_alone_is_refused_before_any_batch(world, prof, tmp_path, monkeypatch, capsys):
    """Batches hold whole genomes (single mode trains on the whole genome), so a genome estimated over the margin
    alone cannot be split: the step is refused (exit 2) before any batch, naming the genome and card B-48."""
    over_budget(world, monkeypatch, tmp_path, largest=1.2)
    root = s0_root(world, tmp_path / "root")
    capsys.readouterr()
    code = build(prof, build_of(world, root, tmp_path / "scratch", scan_jobs=1, scan_cpus=1), stop_after_genecall())
    assert code == runstep.EXIT_REFUSED
    err = capsys.readouterr().err
    assert "B-48" in err and "genome" in err and "refused" in err, err[-3000:]
    assert not genecall_files(tree(root)) and not (root / "runs/stage=S1/genecall.json").exists()


def test_b48_build_stage_says_the_estimate_before_gene_calling(world, prof, tmp_path, capsys):
    root = s0_root(world, tmp_path / "root")
    capsys.readouterr()
    assert build(prof, build_of(world, root, tmp_path / "scratch"), stop_after_genecall()) is None
    err = capsys.readouterr().err
    line = [x for x in err.splitlines() if "S1.genecall: " in x and "estimated at" in x]
    assert len(line) == 1 and "in one process" in line[0] and "14 genomes" in line[0], err[-3000:]


# --- a release built on a stage gene-called in batches -----------------------------------------------------------------


@pytest.mark.skipif(c08.MMSEQS is None, reason="build-stage clusters with MMseqs2 (in the image)")
def test_b48_a_release_on_a_stage_called_in_batches_is_the_release_on_the_one_process_stage(
        world, prof, tmp_path, monkeypatch):
    """S0 then S1 with their releases, every step through build-stage: once with S1 gene-called in one process, once in
    batches. Every file of the root (the later tables: IGRs, tapes, clusters, tokens, postings, dimensions) is the
    same, but the manifests and catalogs, which hold the root's path: their file lists are compared instead. The gene-
    call run record is read as before by the cluster gate and the protein count (`release_proteins`)."""
    k = rate_factor(0.9)
    roots = {}
    for how in ("one", "batches"):
        root = tmp_path / how / "root"
        for stage, release, stages in (("S0", "rB48.0", ("S0",)), ("S1", REL, ("S0", "S1"))):
            g, o = world["inputs"][stage]
            b = bs.Build(root=root, stage=stage, release=release, stages=stages, genomes=g, objects=o, threads=1,
                         tmp=tmp_path / how / "scratch", scan_jobs=1, scan_cpus=1)
            seen: list[str] = []
            with monkeypatch.context() as m:
                if how == "batches" and stage == "S1":
                    scale_rates(m, k)
                code = bs.build_stage(prof, b, run=runner, clock=lambda: NOON_MONDAY, load=lambda: 0.0,
                                      cpu=lambda: IDLE, _hook=seen.append)
            assert code == 0, (how, stage)
            batched = any(p.startswith("after:genecall.batch-") for p in seen)
            assert batched == (how == "batches" and stage == "S1"), (how, stage, seen)
        roots[how] = root
    a, b = (c08.data(tree(roots[h])) for h in ("one", "batches"))
    assert a == b
    assert len(genecall_files(a)) == 5 and any(x.startswith("l3/") for x in a)
    for rel in ("rB48.0", REL):
        ma, mb = (catalog.read_manifest(catalog.manifest_path(roots[h], rel)) for h in ("one", "batches"))
        assert ma["files"] == mb["files"]
    assert bs.release_proteins(roots["one"], ["S0", "S1"], V) == bs.release_proteins(roots["batches"], ["S0", "S1"], V)
    assert bs.stage_residues(roots["one"], "S1", V) == bs.stage_residues(roots["batches"], "S1", V)
    assert "merge" in json.loads((roots["batches"] / "runs/stage=S1/genecall.json").read_text())["argv"]
