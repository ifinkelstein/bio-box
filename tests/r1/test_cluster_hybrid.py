"""Ring R1: hybrid clustering modes (hybrid-a / hybrid-b1 / linclust), PLAN_FORMAT = 2,
B-49 format-1 compatibility and cross-method refusal, and createindex DEC-77 exemption (card B-71).
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

from biobox import catalog, cli, ids, io, profile, runstep
from biobox import stage as bs
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l2 import cluster
from biobox.registry.schemas import table

pytestmark = pytest.mark.covers("B-71", "DEC-77", "B-41", "B-49")

NEED_MMSEQS = pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs not on PATH")
PART = 240.0
PLAN_RECORDS = "l2/cluster_plan/"


def _profile(**kw) -> profile.Profile:
    p = profile.load(Path(__file__).resolve().parents[2] / "profiles" / "wilkcomp02.yaml")
    return dataclasses.replace(p, **kw)


def _synth_root(root: Path, bases: int = 40, seed: int = 71) -> Path:
    """Small synthetic protein corpus with NR90 and NR50 families for fast deterministic tests."""
    rng = random.Random(seed)
    aa = "ACDEFGHIKLMNPQRSTVWY"

    def mutate(s: str, rate: float) -> str:
        return "".join(rng.choice(aa.replace(c, "")) if rng.random() < rate else c for c in s)

    seqs = set()
    for _ in range(bases):
        base = "M" + "".join(rng.choice(aa) for _ in range(rng.randint(79, 220)))
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


def _tree(root: Path, skip: str | None = PLAN_RECORDS) -> dict[str, str]:
    root = Path(root)
    return {
        p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file() and not (skip and p.relative_to(root).as_posix().startswith(skip))
    }


def _cut_doc(
    tmp_path: Path,
    name: str,
    units: int,
    keys: int,
    plan_cut: dict,
    *,
    proteins: int = 10_620_684,
    mode: str = "cascade",
    shared_index: bool = True,
) -> dict:
    by = {c.name: c for c in cluster.calls(mode)}
    c = by[name]
    prefilter = c.module == "prefilter"
    d = tmp_path / f"{mode}-{name}"
    src = d / f"{c.positional[cluster.PER_QUERY[c.module]]}.index"
    src.parent.mkdir(parents=True, exist_ok=True)
    pad = 2 if prefilter else 0
    src.write_text("".join(f"{k}\t{k * 10}\t{units // keys + pad}\n" for k in range(keys)))
    ctx = cluster._Ctx(
        root=d,
        out=None,
        release="rX",
        plan={"format": cluster.PLAN_FORMAT, "method": mode, "proteins": proteins, "cut": plan_cut},
        plan_sha256="",
        w=d,
        threads=1,
        shared_index=shared_index,
    )
    return cluster._cut_doc(ctx, c)


def test_b71_clustering_modes_and_validation(tmp_path):
    assert cluster.CLUSTERING_MODES == ("cascade", "hybrid-a", "hybrid-b1", "linclust")
    assert cluster.DEFAULT_CLUSTERING_MODE == "cascade"
    with pytest.raises(ValueError, match="bogus"):
        cluster.workflow("nr90", "bogus")
    with pytest.raises(ValueError, match="bogus"):
        cluster.calls("bogus")
    with pytest.raises(ValueError, match="bogus"):
        cluster.plan_inputs(tmp_path, "rX", V, ["S0"], cluster.cut_of(parts=1), mode="bogus")
    with pytest.raises(ValueError, match="bogus"):
        cluster.plan(tmp_path / "root", "rX", V, ["S0"], tmp_path / "tmp", parts=1, mode="bogus")
    with pytest.raises(ValueError, match="bogus"):
        cluster.run(tmp_path / "root", "rX", V, ["S0"], tmp=tmp_path / "tmp", mode="bogus")


def test_b71_workflow_and_calls_for_all_modes():
    # cascade: 100% unchanged
    assert cluster.workflow("nr90") == cluster.workflow("nr90", "cascade")
    assert cluster.workflow("nr50") == cluster.workflow("nr50", "cascade")
    assert cluster.calls() == cluster.calls("cascade")
    assert len(cluster.workflow("nr90", "cascade")) == 15
    assert len(cluster.workflow("nr50", "cascade")) == 23
    assert len(cluster.calls("cascade")) == 42

    # hybrid-a: 10-call linclust NR90 + 23-call -s 1,2,3 NR50
    w90_a = cluster.workflow("nr90", "hybrid-a")
    w50_a = cluster.workflow("nr50", "hybrid-a")
    assert len(w90_a) == 10
    assert w90_a[0][0] == "kmermatcher"
    assert w90_a[9][0] == "mergeclusters" and w90_a[9][2] == "clu90"
    assert w50_a == cluster.workflow("nr50", "cascade")
    calls_a = cluster.calls("hybrid-a")
    assert len(calls_a) == 37
    names_a = [c.name for c in calls_a]
    assert "cluster.nr90.09.mergeclusters" in names_a
    assert "cluster.nr90.11.prefilter" not in names_a
    assert "cluster.nr50.15.prefilter" in names_a
    assert "cluster.nr50.19.prefilter" in names_a
    assert names_a[-3:] == ["cluster.merge", "cluster.nr90.members", "cluster.nr50.members"]

    # hybrid-b1: 10-call linclust NR90 + 19-call -s 1,2 NR50 (drops 18..21)
    w90_b1 = cluster.workflow("nr90", "hybrid-b1")
    w50_b1 = cluster.workflow("nr50", "hybrid-b1")
    assert w90_b1 == w90_a
    assert len(w50_b1) == 19
    assert w50_b1[:18] == w50_a[:18]
    assert w50_b1[18][0] == "mergeclusters" and w50_b1[18][2] == "clu50r"
    assert "tmp50/w/clu_step1" in w50_b1[18] and "tmp50/w/clu_step2" not in w50_b1[18]
    calls_b1 = cluster.calls("hybrid-b1")
    assert len(calls_b1) == 33
    names_b1 = [c.name for c in calls_b1]
    assert "cluster.nr50.15.prefilter" in names_b1
    assert "cluster.nr50.18.mergeclusters" in names_b1
    assert "cluster.nr50.19.prefilter" not in names_b1
    assert [c.name for c in calls_b1 if c.shares_index] == ["cluster.nr50.15.prefilter"]

    # linclust: 10-call linclust for both NR90 and NR50
    w90_l = cluster.workflow("nr90", "linclust")
    w50_l = cluster.workflow("nr50", "linclust")
    assert w90_l == w90_a
    assert len(w50_l) == 10
    assert w50_l[9][0] == "mergeclusters" and w50_l[9][2] == "clu50r"
    calls_l = cluster.calls("linclust")
    assert len(calls_l) == 24
    assert not [c for c in calls_l if c.module == "prefilter"]


def test_b71_plan_format_2_and_input_keys(tmp_path):
    assert cluster.PLAN_FORMAT == 2
    assert "method" in cluster._INPUT_KEYS
    assert cluster.plan_method({"format": 1}) == "cascade"
    assert cluster.plan_method({"format": 2, "method": "hybrid-b1"}) == "hybrid-b1"

    base = {
        "release": "rX",
        "stages": ["S0"],
        "gene_call_version": V,
        "mmseqs2": "18.8cc5c",
        "files": [],
        "cut": {"parts": 1},
    }
    w_v1 = cluster._work_dir(tmp_path, tmp_path / "root", {"format": 1, **base})
    w_cas = cluster._work_dir(tmp_path, tmp_path / "root", {"format": 2, "method": "cascade", **base})
    w_b1 = cluster._work_dir(tmp_path, tmp_path / "root", {"format": 2, "method": "hybrid-b1", **base})
    assert w_cas != w_b1
    assert w_v1 != w_b1


def test_b71_params_hash_is_semantic_and_unchanged():
    h90 = cluster.params_hash("nr90")
    h50 = cluster.params_hash("nr50")
    assert len(h90) == 16 and len(h50) == 16 and h90 != h50
    u1, u2 = ids.protein_uid("MAAAA"), ids.protein_uid("MAAAC")
    assert cluster.cluster_uid([u1, u2], h90) == cluster.cluster_uid([u2, u1], h90)


def test_b71_createindex_cpu_min_and_dec77_exemption(tmp_path):
    assert cluster.CREATEINDEX_T0 == 10_000_000
    assert cluster.CREATEINDEX_GROWTH == pytest.approx(0.48)
    assert cluster.CREATEINDEX_RATE_BY_SENS == {"2": 11.80, "3": 15.75}
    assert cluster.createindex_cpu_min(0, 1_000, "2") == 0.0
    assert cluster.createindex_cpu_min(1_000, 0, "2") == 0.0

    # r26.11.3 measured -s 2 point: 24,298,242 targets, 6.476B residues -> ~116.8 CPU-min
    est_r3 = cluster.createindex_cpu_min(6_476_000_000, 24_298_242, "2")
    want_r3 = 11.80 * 1e-9 * 6_476_000_000 * (24_298_242 / 10_000_000) ** 0.48
    assert est_r3 == pytest.approx(want_r3, rel=1e-9)

    # Small input: createindex under budget -> no one_process on cut record
    small_units, small_keys = 500_000_000, 2_000
    small = _cut_doc(
        tmp_path / "small",
        "cluster.nr50.15.prefilter",
        small_units,
        small_keys,
        {"part_cpu_min": PART, "one_process_cpu_min_max": 20_000.0},
        mode="hybrid-b1",
    )
    if len(small["parts"]) == 1:
        small = _cut_doc(
            tmp_path / "small2",
            "cluster.nr50.15.prefilter",
            small_units,
            small_keys,
            {"parts": 2, "one_process_cpu_min_max": 20_000.0},
            mode="hybrid-b1",
        )
    assert len(small["parts"]) > 1
    assert "index_cpu_min" in small
    assert small["createindex_cpu_min"] == pytest.approx(
        cluster.createindex_cpu_min(small_units, small_keys, "2"), abs=1e-5
    )
    assert "one_process" not in small

    # Large residues so createindex_cpu_min > part_cpu_min (240) and <= one_process_cpu_min_max (20000)
    # 1_000 keys, 4e12 residues -> createindex_cpu_min(4e12, 1000, "2") = 11.80 * 4000 * (1e-4)^0.48 = ~567.8 CPU-min
    big_units, big_keys = 4_000_000_000_000, 1_000
    ci_est = cluster.createindex_cpu_min(big_units, big_keys, "2")
    assert 300.0 < ci_est < 20_000.0
    exempt = _cut_doc(
        tmp_path / "exempt",
        "cluster.nr50.15.prefilter",
        big_units,
        big_keys,
        {"part_cpu_min": PART, "one_process_cpu_min_max": 20_000.0},
        mode="hybrid-b1",
    )
    assert len(exempt["parts"]) > 1
    assert exempt["createindex_cpu_min"] == pytest.approx(ci_est, abs=1e-5)
    assert exempt["one_process"]["decision"] == "DEC-77"
    assert exempt["one_process"]["cpu_min"] == pytest.approx(ci_est, abs=1e-5)
    assert exempt["one_process"]["cpu_min_max"] == 20_000.0

    # Without cap: refused naming DEC-77
    with pytest.raises(ValueError, match="DEC-77"):
        _cut_doc(
            tmp_path / "nocap",
            "cluster.nr50.15.prefilter",
            big_units,
            big_keys,
            {"part_cpu_min": PART},
            mode="hybrid-b1",
        )

    # Above cap: refused naming DEC-77 and the cap
    with pytest.raises(ValueError, match=r"over the one-process cap of 300.*DEC-77"):
        _cut_doc(
            tmp_path / "overcap",
            "cluster.nr50.15.prefilter",
            big_units,
            big_keys,
            {"part_cpu_min": PART, "one_process_cpu_min_max": 300.0},
            mode="hybrid-b1",
        )


def test_b71_stage_and_cli_cluster_method_wiring(tmp_path):
    prof = _profile(cluster_cpu_min_max=20_000)
    assert bs.DEC77_EXEMPT == "exempt: DEC-77"
    b_def = bs.Build(root=tmp_path / "r", stage="S0", release="r1", genomes=tmp_path / "g.json",
                     objects=tmp_path / "o.json")
    assert b_def.cluster_method == "cascade"
    b_b1 = dataclasses.replace(b_def, cluster_method="hybrid-b1")
    make_def = bs.cluster_steps(prof, b_def)
    make_b1 = bs.cluster_steps(prof, b_b1)
    assert "--cluster-method" not in make_def("plan").argv
    assert "--cluster-method" in make_b1("plan").argv
    i = make_b1("plan").argv.index("--cluster-method")
    assert make_b1("plan").argv[i + 1] == "hybrid-b1"
    with pytest.raises(ValueError, match="bogus"):
        bs.steps(prof, dataclasses.replace(b_def, cluster_method="bogus"))


@NEED_MMSEQS
@pytest.mark.parametrize("mode", ["hybrid-b1", "hybrid-a", "linclust"])
def test_b71_deterministic_byte_identity_across_threads_and_split(tmp_path, mode):
    root = _synth_root(tmp_path / "root")
    out_t1 = tmp_path / "out-t1"
    out_t4 = tmp_path / "out-t4"
    out_sp = tmp_path / "out-split"
    (tmp_path / "tmp-t1").mkdir()
    (tmp_path / "tmp-t4").mkdir()
    cluster.run(root, "rB71", V, ["S0"], threads=1, tmp=tmp_path / "tmp-t1", out=out_t1, mode=mode)
    cluster.run(root, "rB71", V, ["S0"], threads=4, tmp=tmp_path / "tmp-t4", out=out_t4, mode=mode)
    cluster.run_split(root, "rB71", V, ["S0"], tmp=tmp_path / "tmp-sp", parts=2, threads=2, out=out_sp, mode=mode)

    t1 = _tree(out_t1)
    assert _tree(out_t4) == t1
    assert _tree(out_sp) == t1

    plan_doc = cluster.read_record(root, "rB71", cluster.PLAN_STEP, out=out_sp)
    assert plan_doc["format"] == 2
    assert plan_doc["method"] == mode
    if mode == "hybrid-b1":
        assert cluster.read_record(root, "rB71", "cluster.nr90.09.mergeclusters", out=out_sp) is not None
        assert cluster.read_record(root, "rB71", "cluster.nr90.11.prefilter", out=out_sp) is None
        assert cluster.read_record(root, "rB71", "cluster.nr50.15.prefilter", out=out_sp) is not None
        assert cluster.read_record(root, "rB71", "cluster.nr50.18.mergeclusters", out=out_sp) is not None
        assert cluster.read_record(root, "rB71", "cluster.nr50.19.prefilter", out=out_sp) is None


@NEED_MMSEQS
def test_b71_plan_format_1_compat_and_no_cross_method_reuse(tmp_path):
    root = _synth_root(tmp_path / "root", bases=25, seed=72)

    def _fake_cluster_run_record(release: str) -> None:
        outputs = []
        for sub in (
            f"l2/cluster_member/release={release}",
            f"l2/cluster_member_aln/release={release}",
            f"l2/member_cluster/release={release}",
            f"ledger/release={release}/step=cluster",
        ):
            for p in sorted((root / sub).glob("*.parquet")):
                outputs.append({
                    "path": p.relative_to(root).as_posix(),
                    "bytes": p.stat().st_size,
                    "sha256": ids.sha256_file(p).hex(),
                })
        rec_path = root / f"runs/stage={release}/cluster.json"
        rec_path.parent.mkdir(parents=True, exist_ok=True)
        rec_path.write_text(json.dumps({
            "stage": release, "step": "cluster", "exit_code": 0, "outputs": outputs,
        }))

    # 1. Build R1 in cascade mode, then downgrade its cluster.plan.json to PLAN_FORMAT=1 (no "method")
    cluster.run_split(root, "r1", V, ["S0"], tmp=tmp_path / "tmp-r1", parts=1, threads=1, mode="cascade")
    _fake_cluster_run_record("r1")
    plan_path = root / cluster.plan_dir("r1") / f"{cluster.PLAN_STEP}.json"
    p1 = json.loads(plan_path.read_bytes())
    p1_v1 = {k: v for k, v in p1.items() if k != "method"} | {"format": 1}
    raw_v1 = cluster._json_bytes(p1_v1)
    sha_v1 = hashlib.sha256(raw_v1).hexdigest()
    plan_path.write_bytes(raw_v1)
    for rec_file in (root / cluster.plan_dir("r1")).glob("*.json"):
        if rec_file.name == f"{cluster.PLAN_STEP}.json":
            continue
        d = json.loads(rec_file.read_bytes())
        d["format"] = 1
        d["plan_sha256"] = sha_v1
        rec_file.write_bytes(cluster._json_bytes(d))

    # Format-1 plan still reads cleanly as "cascade" and reuses for a format-2 "cascade" release
    cut = cluster.cut_of(parts=1)
    assert cluster.reuse_source(root, "r2", V, ["S0"], cut, tmp=tmp_path / "tmp-u", mode="cascade") == "r1"
    res = cluster.reuse(root, "r2", V, ["S0"], tmp_path / "tmp-r2", parts=1, mode="cascade")
    assert res["reused_from"] == "r1"
    p2 = cluster.read_record(root, "r2", cluster.PLAN_STEP)
    assert p2["format"] == 2 and p2["method"] == "cascade"

    # Cross-method reuse between cascade (r1) and hybrid-b1 MUST be refused
    assert cluster.reuse_source(root, "r3", V, ["S0"], cut, tmp=tmp_path / "tmp-u", mode="hybrid-b1") is None
    with pytest.raises(ValueError, match="method"):
        cluster.reuse(root, "r3", V, ["S0"], tmp_path / "tmp-r3", source="r1", parts=1, mode="hybrid-b1")

    # Build R3 in hybrid-b1 mode -> R4 in hybrid-b1 reuses R3, R4 in cascade does NOT
    cluster.run_split(root, "r3", V, ["S0"], tmp=tmp_path / "tmp-r3b", parts=1, threads=1, mode="hybrid-b1")
    _fake_cluster_run_record("r3")
    assert cluster.reuse_source(root, "r4", V, ["S0"], cut, tmp=tmp_path / "tmp-u", mode="hybrid-b1") == "r3"
    assert cluster.reuse_source(root, "r4", V, ["S0"], cut, tmp=tmp_path / "tmp-u", mode="hybrid-a") is None
    res4 = cluster.reuse(root, "r4", V, ["S0"], tmp_path / "tmp-r4", parts=1, mode="hybrid-b1")
    assert res4["reused_from"] == "r3"
