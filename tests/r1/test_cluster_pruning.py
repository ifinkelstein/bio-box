"""Ring R1: dead-intermediate clustering scratch pruning (`clean_consumed`), milestone-aware
`lost()`/`redo()`, and level-split `cluster finish` (`cluster.finish.nr90`, `cluster.finish.nr50`,
lightweight `cluster`) (card B-72, DEC-77).
"""

from __future__ import annotations

import hashlib
import json
import random
import shutil
from pathlib import Path

import pyarrow as pa
import pytest

from biobox import ids, io
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l2 import cluster
from biobox.registry.schemas import table

pytestmark = pytest.mark.covers("B-72", "DEC-77", "B-41", "B-49")

NEED_MMSEQS = pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs not on PATH")
PLAN_RECORDS = "l2/cluster_plan/"


def _synth_root(root: Path, bases: int = 80, seed: int = 72) -> Path:
    """Synthetic protein corpus with NR90 and NR50 families large enough for clear scratch measurement."""
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
    t = pa.table({
        "protein_uid": [u for u, _ in rows],
        "protein_digest": [hashlib.sha256(s.encode()).digest() for _, s in rows],
        "length": [len(s) for _, s in rows],
        "seq": [s for _, s in rows],
    })
    io.write_table(
        t,
        table("protein_seq"),
        root / f"l0seq/protein_seq/gene_call_version={V}/stage=S0/protein_prefix=all/part-0.parquet",
    )
    return root


def _dir_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file() and not p.is_symlink())


def _tree(root: Path, skip: str | None = PLAN_RECORDS) -> dict[str, str]:
    root = Path(root)
    return {
        p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file() and not (skip and p.relative_to(root).as_posix().startswith(skip))
    }


def test_b72_finish_level_step_constants():
    assert cluster.FINISH_LEVEL_STEPS == ("cluster.finish.nr90", "cluster.finish.nr50")
    assert hasattr(cluster, "clean_consumed")


@NEED_MMSEQS
@pytest.mark.parametrize("mode", ["cascade", "hybrid-b1"])
def test_b72_clean_consumed_rules_and_peak_scratch_reduction(tmp_path, mode):
    """`clean_consumed` unlinks dead intermediates only after verified milestones, never touches
    `clu_redundancy` or `clu_step*` before final `nr50.*.mergeclusters`, and cuts peak scratch by >= 40%.
    Uses `parts=1` so the fixed ~900 MB MMseqs2 `createindex` hash table does not dominate a ~1k-protein fixture."""
    root = _synth_root(tmp_path / "root")

    # 1. Measure peak scratch WITHOUT pruning (unpruned baseline)
    tmp_no, out_no = tmp_path / f"tmp-no-{mode}", tmp_path / f"out-no-{mode}"
    cluster.plan(root, "rNo", V, ["S0"], tmp_no, parts=1, out=out_no, mode=mode)
    w_no = cluster.work_dir(tmp_no, root, "rNo", out=out_no)
    peak_unpruned = _dir_bytes(w_no)
    while True:
        todo = [n for n in cluster.sequence(root, "rNo", out_no) if cluster.read_record(root, "rNo", n, out_no) is None]
        if not todo:
            break
        name = todo[0]
        cluster.step(root, "rNo", name, tmp_no, threads=2, out=out_no)
        peak_unpruned = max(peak_unpruned, _dir_bytes(w_no))

    # 2. Measure peak scratch WITH `clean_consumed` and verify exact safety invariants
    tmp_pr, out_pr = tmp_path / f"tmp-pr-{mode}", tmp_path / f"out-pr-{mode}"
    cluster.plan(root, "rPr", V, ["S0"], tmp_pr, parts=1, out=out_pr, mode=mode)
    w_pr = cluster.work_dir(tmp_pr, root, "rPr", out=out_pr)
    peak_pruned = _dir_bytes(w_pr)
    final_nr50_merge = "cluster.nr50.22.mergeclusters" if mode == "cascade" else "cluster.nr50.18.mergeclusters"
    while True:
        todo = [n for n in cluster.sequence(root, "rPr", out_pr) if cluster.read_record(root, "rPr", n, out_pr) is None]
        if not todo:
            break
        name = todo[0]
        cluster.step(root, "rPr", name, tmp_pr, threads=2, out=out_pr)
        c, kind, _ = cluster._parse(name, mode)
        if kind == "call":
            cluster.clean_consumed(root, "rPr", c.name, tmp_pr, out=out_pr)
            if c.name == "cluster.rep90":
                assert not (w_pr / "tmp90/w").exists()
                assert (w_pr / "clu90.dbtype").is_file()
                assert (w_pr / "rep90.dbtype").is_file()
            elif c.name == "cluster.nr50.09.mergeclusters":
                assert not (w_pr / "tmp50/w/linclust").exists()
                assert (w_pr / "tmp50/w/clu_redundancy.dbtype").is_file()
            elif c.name == "cluster.nr50.13.clust":
                assert not (w_pr / "tmp50/w/pref_step0.dbtype").exists()
                assert not (w_pr / "tmp50/w/aln_step0.dbtype").exists()
                # clu_redundancy and clu_step0 MUST remain intact for final mergeclusters!
                assert (w_pr / "tmp50/w/clu_redundancy.dbtype").is_file()
                assert (w_pr / "tmp50/w/clu_step0.dbtype").is_file()
                # input_step_redundancy is still needed by nr50.14.createsubdb!
                assert (w_pr / "tmp50/w/input_step_redundancy.dbtype").is_file()
            elif c.name == "cluster.nr50.14.createsubdb":
                assert not (w_pr / "tmp50/w/input_step_redundancy.dbtype").exists()
                assert (w_pr / "tmp50/w/clu_redundancy.dbtype").is_file()
                assert (w_pr / "tmp50/w/clu_step0.dbtype").is_file()
            elif c.name == "cluster.nr50.17.clust":
                assert not (w_pr / "tmp50/w/pref_step1.dbtype").exists()
                assert not (w_pr / "tmp50/w/aln_step1.dbtype").exists()
                assert (w_pr / "tmp50/w/clu_redundancy.dbtype").is_file()
                assert (w_pr / "tmp50/w/clu_step0.dbtype").is_file()
                assert (w_pr / "tmp50/w/clu_step1.dbtype").is_file()
            elif c.name == final_nr50_merge:
                assert not (w_pr / "tmp50/w").exists()
                assert (w_pr / "clu50r.dbtype").is_file()
        peak_pruned = max(peak_pruned, _dir_bytes(w_pr))

    reduction = 1.0 - (peak_pruned / peak_unpruned)
    print(f"\n[B-72 scratch {mode}] unpruned_peak={peak_unpruned:,} B  pruned_peak={peak_pruned:,} B  reduction={reduction:.2%}")
    assert reduction >= 0.40, (
        f"mode={mode}: peak scratch reduction {reduction:.1%} < 40% "
        f"(unpruned={peak_unpruned:,} B, pruned={peak_pruned:,} B)"
    )

    # Pruned intermediates are absent, milestones intact -> lost() MUST be empty!
    assert cluster.lost(root, "rPr", tmp_pr, out=out_pr) == []


@NEED_MMSEQS
def test_b72_clean_consumed_refuses_unrecorded_or_corrupted_milestone(tmp_path):
    root = _synth_root(tmp_path / "root", bases=20, seed=73)
    tmp, out = tmp_path / "tmp", tmp_path / "out"
    cluster.run_split(root, "rRef", V, ["S0"], tmp=tmp, parts=1, threads=1, out=out,
                      until="cluster.rep90", mode="hybrid-b1")
    with pytest.raises(FileNotFoundError):
        cluster.clean_consumed(root, "rRef", "cluster.rep90", tmp, out=out)
    cluster.step(root, "rRef", "cluster.rep90", tmp, threads=1, out=out)
    w = cluster.work_dir(tmp, root, "rRef", out=out)
    # Corrupt rep90 -> clean_consumed must refuse and NOT delete tmp90/w
    (w / "rep90.dbtype").unlink()
    with pytest.raises(ValueError):
        cluster.clean_consumed(root, "rRef", "cluster.rep90", tmp, out=out)
    assert (w / "tmp90/w").is_dir()


@NEED_MMSEQS
@pytest.mark.parametrize("mode", ["cascade", "hybrid-b1"])
def test_b72_level_split_finish_byte_identity_and_scratch_cleanup(tmp_path, mode):
    """Level-split `cluster.finish.nr90` -> `cluster.finish.nr50` -> `cluster` with `clean_consumed`
    produces 100% byte-identical `cluster_member`, `cluster_member_aln`, `member_cluster`, and ledger
    against `cluster.run()` and unsplit `run_split()`."""
    root = _synth_root(tmp_path / "root", bases=35, seed=74)
    out_run = tmp_path / f"out-run-{mode}"
    out_unsplit = tmp_path / f"out-unsplit-{mode}"
    out_split = tmp_path / f"out-split-{mode}"
    (tmp_path / f"tmp-run-{mode}").mkdir()

    cluster.run(root, "rB72", V, ["S0"], threads=2, tmp=tmp_path / f"tmp-run-{mode}", out=out_run, mode=mode)
    cluster.run_split(root, "rB72", V, ["S0"], tmp=tmp_path / f"tmp-unsplit-{mode}", parts=2, threads=2,
                      out=out_unsplit, mode=mode, prune=False, split_finish=False)
    cluster.run_split(root, "rB72", V, ["S0"], tmp=tmp_path / f"tmp-split-{mode}", parts=2, threads=2,
                      out=out_split, mode=mode, prune=True, split_finish=True)

    want = _tree(out_run)
    assert _tree(out_unsplit) == want
    assert _tree(out_split) == want

    # Level-split finish wrote records for both finish levels
    for fs in cluster.FINISH_LEVEL_STEPS:
        rec = cluster.read_record(root, "rB72", fs, out=out_split)
        assert rec is not None and rec["module"] == "finish"


@NEED_MMSEQS
def test_b72_lost_and_redo_milestone_recovery(tmp_path):
    """When pruned intermediates are absent and milestones are intact, `lost()` is empty; when a
    milestone (`clu90`, `clu50r`, or `finish/nr90/*`) is deleted, `lost()`/`redo()` rebuilds the
    required upstream steps and reproduces exact records and byte-identical tables."""
    root = _synth_root(tmp_path / "root", bases=25, seed=75)
    oracle = tmp_path / "oracle"
    (tmp_path / "tmp-oracle").mkdir()
    cluster.run(root, "rRec", V, ["S0"], threads=2, tmp=tmp_path / "tmp-oracle", out=oracle, mode="hybrid-b1")

    tmp, out = tmp_path / "tmp", tmp_path / "out"
    # Run with pruning up to cluster.finish.nr50 (so cluster.finish.nr90 is also published)
    cluster.run_split(root, "rRec", V, ["S0"], tmp=tmp, parts=2, threads=2, out=out,
                      until="cluster.finish.nr50", mode="hybrid-b1", prune=True, split_finish=True)
    w = cluster.work_dir(tmp, root, "rRec", out=out)

    # Intermediates pruned, milestones intact -> lost() is empty
    assert not (w / "tmp90/w").exists()
    assert not (w / "tmp50/w").exists()
    assert not (w / "finish/tmp_nr90").exists()
    assert cluster.lost(root, "rRec", tmp, out=out) == []

    # 1. Delete finish/nr90/cluster_member.parquet -> must rebuild aln_nr90 (if pruned) + cluster.finish.nr90
    (w / "finish/nr90/cluster_member.parquet").unlink()
    lost_fin = cluster.lost(root, "rRec", tmp, out=out)
    assert "cluster.finish.nr90" in lost_fin
    for name in lost_fin:
        assert cluster.redo(root, "rRec", name, tmp, threads=1, out=out)["verified"]
    assert cluster.lost(root, "rRec", tmp, out=out) == []

    # 2. Delete clu90 milestone -> must rebuild nr90 linclust steps + clu90 producer
    for p in w.glob("clu90*"):
        if p.is_file():
            p.unlink()
    lost_90 = cluster.lost(root, "rRec", tmp, out=out)
    assert "cluster.nr90.00.kmermatcher" in lost_90
    assert "cluster.nr90.09.mergeclusters" in lost_90
    for name in lost_90:
        assert cluster.redo(root, "rRec", name, tmp, threads=1, out=out)["verified"]
    assert cluster.lost(root, "rRec", tmp, out=out) == []

    # Finish the build and verify 100% byte-identity against oracle
    cluster.run_split(root, "rRec", V, ["S0"], tmp=tmp, parts=2, threads=2, out=out,
                      mode="hybrid-b1", prune=True, split_finish=True)
    assert _tree(out) == _tree(oracle)


@NEED_MMSEQS
@pytest.mark.parametrize("stop_at", ["cluster.nr50.00.kmermatcher", "cluster.finish.nr50"])
def test_b72_crash_resume_mid_run_byte_identical(tmp_path, stop_at):
    """Stopping mid-run after `cluster.rep90` (before `nr50.00`) or after `cluster.finish.nr90`
    (before `cluster.finish.nr50`) and resuming produces 100% byte-identical final tables."""
    root = _synth_root(tmp_path / "root", bases=25, seed=76)
    oracle = tmp_path / "oracle"
    (tmp_path / "tmp-oracle").mkdir()
    cluster.run(root, "rRes", V, ["S0"], threads=2, tmp=tmp_path / "tmp-oracle", out=oracle, mode="hybrid-b1")

    tmp, out = tmp_path / "tmp", tmp_path / "out"
    cluster.run_split(root, "rRes", V, ["S0"], tmp=tmp, parts=2, threads=2, out=out,
                      until=stop_at, mode="hybrid-b1", prune=True, split_finish=True)
    w = cluster.work_dir(tmp, root, "rRes", out=out)
    assert not (w / "tmp90/w").exists()
    assert cluster.lost(root, "rRes", tmp, out=out) == []
    if stop_at == "cluster.finish.nr50":
        assert cluster.read_record(root, "rRes", "cluster.finish.nr90", out=out) is not None
        assert not (w / "finish/tmp_nr90").exists()
        assert (w / "finish/nr90/cluster_member.parquet").is_file()

    # Resume to completion
    cluster.run_split(root, "rRes", V, ["S0"], tmp=tmp, parts=2, threads=2, out=out,
                      mode="hybrid-b1", prune=True, split_finish=True)
    assert _tree(out) == _tree(oracle)
