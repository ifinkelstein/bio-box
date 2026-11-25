"""Ring R1: unit tests for card C-38c (linclust alone vs cascade on r26.11.1's proteins).

Tests written FIRST before implementing `scripts/stress/c38c_linclust_vs_cascade.py`
(`AGENTS.md` §5.1). Covers:
- Extraction of the exact redundancy-phase linclust module calls from `biobox.l2.cluster.workflow`
  at NR90 and NR50 (`--shuffle 0`, release coverage/identity parameters, no prefilter calls).
- Constructed-truth cluster comparison metrics in bounded DuckDB (clusters & singletons per level,
  cascade NR50 -> linclust cluster split distribution, missed-join identity-band histogram).
- `N19` / `N20` 20%-rule comparison helpers (`AGENTS.md` §6).
- Extrapolation math to 350,000,000 and 1,000,000,000 proteins.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox.l2 import cluster
from biobox.registry.schemas import table

pytestmark = pytest.mark.covers("C-38c")  # ring R1 by directory (tests/conftest.py)

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "stress" / "c38c_linclust_vs_cascade.py"


def _load():
    spec = importlib.util.spec_from_file_location("c38c_linclust_vs_cascade", SCRIPT)
    assert spec is not None and spec.loader is not None, f"missing script: {SCRIPT}"
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c38c_linclust_vs_cascade"] = mod
    spec.loader.exec_module(mod)
    return mod


def _uid(n: int) -> bytes:
    return n.to_bytes(16, "big")


# ---------------------------------------------------------------------------
# 1. Extraction of redundancy-phase linclust calls from cluster.workflow
# ---------------------------------------------------------------------------


def test_c38c_extracts_exact_redundancy_phase_calls_from_workflow():
    mod = _load()
    for level, min_id in (("nr90", "0.9"), ("nr50", "0.5")):
        wf = cluster.workflow(level)
        got = mod.linclust_redundancy_calls(level)
        # Exactly the 10 calls of the redundancy phase (up to mergeclusters -> clu_redundancy)
        assert len(got) == 10
        assert [c.argv for c in got] == wf[:10]
        modules = [c.module for c in got]
        assert "prefilter" not in modules
        assert modules == [
            "kmermatcher",
            "rescorediagonal",
            "clust",
            "createsubdb",
            "createsubdb",
            "filterdb",
            "rescorediagonal",
            "align",
            "clust",
            "mergeclusters",
        ]
        # Release identity and coverage parameters preserved verbatim
        km = got[0].argv
        assert km[km.index("--min-seq-id") + 1] == min_id
        assert km[km.index("-c") + 1] == "0.8"
        assert km[km.index("--cov-mode") + 1] == "1"
        al = got[7].argv
        assert al[al.index("--min-seq-id") + 1] == min_id
        assert al[al.index("-c") + 1] == "0.8"
        assert al[al.index("--cov-mode") + 1] == "1"
        cl = got[8].argv
        assert cl[cl.index("--cluster-mode") + 1] == "2"


def test_c38c_pipeline_calls_and_createdb_shuffle_zero():
    mod = _load()
    cdb = mod.createdb_argv("in.faa", "db")
    assert cdb[0] == "createdb"
    assert cdb[cdb.index("--shuffle") + 1] == "0"
    assert cdb[cdb.index("--dbtype") + 1] == "1"

    pipe = mod.linclust_pipeline_calls()
    names = [c.name for c in pipe]
    assert len(pipe) == 22  # 10 nr90 + rep90 + 10 nr50 + merge
    assert all(c.module != "prefilter" for c in pipe)
    assert "cluster.rep90" in names
    assert "cluster.merge" in names
    assert pipe[9].output == "clu90"
    assert pipe[10].argv == ("createsubdb", "clu90", "db", "rep90")
    assert pipe[20].output == "clu50r"
    assert pipe[21].argv == ("mergeclusters", "db", "clu50", "clu90", "clu50r")


# ---------------------------------------------------------------------------
# 2. Constructed-truth cluster comparison metrics (bounded DuckDB)
# ---------------------------------------------------------------------------


def _write_constructed_world(tmp_path: Path):
    """Build a small hand-verifiable clustering world of 12 proteins (uids 1..12).

    Cascade NR90:
      C90_1 = {1, 2} (rep 1; pident 0.95)
      C90_2 = {3, 4} (rep 3; pident 0.92)
      singletons: {5}, {6}, {7}, {8}, {9}, {10}, {11}, {12}
      -> 10 clusters, 8 singletons

    Linclust NR90:
      L90_1 = {1, 2} (rep 1)
      singletons: {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}, {12}
      -> 11 clusters, 10 singletons (missed join: 4 from rep 3 at pident 0.92 -> >=90%)

    Cascade NR50 (3 multi-member clusters + 2 singletons = 5 clusters):
      C50_A (size 5, rep 1): {1, 2(0.95), 3(0.75), 4(0.65), 5(0.55)}
      C50_B (size 3, rep 6): {6, 7(0.85), 8(0.45)}
      C50_C (size 2, rep 9): {9, 10(0.62)}
      C50_D (size 1, rep 11): {11}
      C50_E (size 1, rep 12): {12}
      -> 5 clusters, 2 singletons

    Linclust NR50:
      C50_A splits into 3 linclust clusters: {1, 2}, {3, 4}, {5}
        - rep 1 is in {1, 2}: member 2 joined; members 3 (0.75, 70-80%), 4 (0.65, 60-70%), 5 (0.55, 50-60%) missed;
          of those, only 5 is a linclust singleton.
      C50_B stays intact as 1 linclust cluster: {6, 7, 8} -> 0 missed
      C50_C splits into 2 singletons: {9}, {10} -> member 10 (0.62, 60-70%) missed & singleton
      C50_D stays {11}
      C50_E stays {12}
      -> 8 clusters, 4 singletons ({5}, {9}, {10}, {11}, {12} -> wait: 5 singletons!)
      Let's count linclust NR50 clusters: {1,2}, {3,4}, {5}, {6,7,8}, {9}, {10}, {11}, {12} -> 8 clusters, 5 singletons!
    """
    mem_rows = []
    aln_rows = []

    def add_cluster(level: str, cid: int, rep: int, members_pident: list[tuple[int, float]]):
        cu = _uid(cid)
        ru = _uid(rep)
        for m, pid in members_pident:
            mu = _uid(m)
            mem_rows.append({
                "member_uid": mu,
                "level": level,
                "cluster_uid": cu,
                "rep_uid": ru,
                "cluster_params_hash": f"test-{level}",
            })
            aln_rows.append({
                "member_uid": mu,
                "level": level,
                "rep_uid": ru,
                "qstart": 0,
                "qend": 100,
                "tstart": 0,
                "tend": 100,
                "pident": float(pid),
                "alnlen": 100,
                "cigar": b"100M",
            })

    # Cascade NR90
    add_cluster("nr90", 901, 1, [(1, 1.0), (2, 0.95)])
    add_cluster("nr90", 902, 3, [(3, 1.0), (4, 0.92)])
    for i in range(5, 13):
        add_cluster("nr90", 900 + i, i, [(i, 1.0)])

    # Cascade NR50
    add_cluster("nr50", 501, 1, [(1, 1.0), (2, 0.95), (3, 0.75), (4, 0.65), (5, 0.55)])
    add_cluster("nr50", 502, 6, [(6, 1.0), (7, 0.85), (8, 0.45)])
    add_cluster("nr50", 503, 9, [(9, 1.0), (10, 0.62)])
    add_cluster("nr50", 504, 11, [(11, 1.0)])
    add_cluster("nr50", 505, 12, [(12, 1.0)])

    mem_path = tmp_path / "cluster_member.parquet"
    aln_path = tmp_path / "cluster_member_aln.parquet"
    pq.write_table(pa.Table.from_pylist(mem_rows, schema=table("cluster_member").schema), mem_path)
    pq.write_table(pa.Table.from_pylist(aln_rows, schema=table("cluster_member_aln").schema), aln_path)

    # Linclust TSVs (rep_hex \t member_hex)
    def write_tsv(path: Path, pairs: list[tuple[int, int]]):
        path.write_text("".join(f"{_uid(r).hex()}\t{_uid(m).hex()}\n" for r, m in pairs))

    lin90 = tmp_path / "clu90.tsv"
    write_tsv(lin90, [(1, 1), (1, 2)] + [(i, i) for i in range(3, 13)])

    lin50 = tmp_path / "clu50.tsv"
    write_tsv(
        lin50,
        [(1, 1), (1, 2), (3, 3), (3, 4), (5, 5), (6, 6), (6, 7), (6, 8), (9, 9), (10, 10), (11, 11), (12, 12)],
    )
    return [mem_path], [aln_path], {"nr90": lin90, "nr50": lin50}


def test_c38c_constructed_truth_comparison_metrics(tmp_path: Path):
    mod = _load()
    mem_paths, aln_paths, lin_tsvs = _write_constructed_world(tmp_path)
    res = mod.compare_clusterings_duckdb(mem_paths, aln_paths, lin_tsvs, memory_limit="512MB", tmp_dir=tmp_path)

    # Level summaries
    nr90 = res["levels"]["nr90"]
    assert nr90["n_proteins"] == 12
    assert nr90["cascade_clusters"] == 10
    assert nr90["cascade_singletons"] == 8
    assert nr90["linclust_clusters"] == 11
    assert nr90["linclust_singletons"] == 10

    nr50 = res["levels"]["nr50"]
    assert nr50["n_proteins"] == 12
    assert nr50["cascade_clusters"] == 5
    assert nr50["cascade_singletons"] == 2
    assert nr50["linclust_clusters"] == 8
    assert nr50["linclust_singletons"] == 5

    # Cascade NR50 -> linclust NR50 split distribution
    splits_all = res["nr50_splits"]["all"]
    assert splits_all["n_cascade_clusters"] == 5
    assert splits_all["unsplit_clusters"] == 3  # B, D, E
    assert splits_all["split_clusters"] == 2  # A (3 parts), C (2 parts)
    assert splits_all["max_linclust_clusters"] == 3
    assert math.isclose(splits_all["mean_linclust_clusters"], (3 + 1 + 2 + 1 + 1) / 5, rel_tol=1e-9)

    splits_ge2 = res["nr50_splits"]["size>=2"]
    assert splits_ge2["n_cascade_clusters"] == 3  # A, B, C
    assert splits_ge2["unsplit_clusters"] == 1  # B
    assert splits_ge2["split_clusters"] == 2  # A, C
    assert math.isclose(splits_ge2["mean_linclust_clusters"], (3 + 1 + 2) / 3, rel_tol=1e-9)

    # Missed joins by identity band at NR50
    bands50 = {r["band"]: r for r in res["missed_joins_by_band"]["nr50"]}
    # Total non-rep cascade NR50 members: 2(0.95), 3(0.75), 4(0.65), 5(0.55), 7(0.85), 8(0.45), 10(0.62) = 7
    assert bands50["<50%"]["cascade_nonrep"] == 1
    assert bands50["<50%"]["missed_by_linclust"] == 0
    assert bands50["50-60%"]["cascade_nonrep"] == 1
    assert bands50["50-60%"]["missed_by_linclust"] == 1
    assert bands50["50-60%"]["missed_to_singleton"] == 1
    assert bands50["60-70%"]["cascade_nonrep"] == 2
    assert bands50["60-70%"]["missed_by_linclust"] == 2
    assert bands50["60-70%"]["missed_to_singleton"] == 1  # member 10 is singleton; member 4 joined rep 3
    assert bands50["70-80%"]["cascade_nonrep"] == 1
    assert bands50["70-80%"]["missed_by_linclust"] == 1
    assert bands50["70-80%"]["missed_to_singleton"] == 0  # member 3 is rep of {3, 4}
    assert bands50["80-90%"]["cascade_nonrep"] == 1
    assert bands50["80-90%"]["missed_by_linclust"] == 0
    assert bands50[">=90%"]["cascade_nonrep"] == 1
    assert bands50[">=90%"]["missed_by_linclust"] == 0

    # Missed joins at NR90: member 4 at 0.92 (>=90%)
    bands90 = {r["band"]: r for r in res["missed_joins_by_band"]["nr90"]}
    assert bands90[">=90%"]["cascade_nonrep"] == 2
    assert bands90[">=90%"]["missed_by_linclust"] == 1
    assert bands90[">=90%"]["missed_to_singleton"] == 1


# ---------------------------------------------------------------------------
# 3. N19 / N20 20%-rule helpers and extrapolation math
# ---------------------------------------------------------------------------


def test_c38c_twenty_percent_rule_and_n19_n20_helpers():
    mod = _load()
    within = mod.check_twenty_percent_rule(measured=380.0, estimate=336.0)
    assert not within["exceeds_20pct"]
    assert math.isclose(within["pct_diff"], (380.0 - 336.0) / 336.0 * 100.0, rel_tol=1e-9)

    over = mod.check_twenty_percent_rule(measured=500.0, estimate=336.0)
    assert over["exceeds_20pct"]

    under = mod.check_twenty_percent_rule(measured=250.0, estimate=336.0)
    assert under["exceeds_20pct"]
    assert under["pct_diff"] < -20.0

    # Design note constants N19, N20
    assert math.isclose(mod.N19_NR90_SHARE_EST, 0.074, rel_tol=1e-9)
    assert math.isclose(mod.N19_BOTH_SHARE_EST, 0.06, rel_tol=1e-9)
    assert math.isclose(mod.N19_CPU_PER_M_EST, 14.0, rel_tol=1e-9)
    assert math.isclose(mod.N20_BYTES_PER_SEQ_EST, 336.0, rel_tol=1e-9)


def test_c38c_extrapolation_to_350m_and_1b():
    mod = _load()
    n_meas = 49_102_041
    cpu_meas = 500.0
    rss_meas = 20 * (1024**3)  # 20 GiB
    rows = mod.extrapolate_linclust(n_meas, cpu_meas, rss_meas, target_ns=(350_000_000, 1_000_000_000))
    assert [r["target_proteins"] for r in rows] == [350_000_000, 1_000_000_000]

    r350, r1b = rows
    scale350 = 350_000_000 / n_meas
    nlogn350 = (350_000_000 * math.log2(350_000_000)) / (n_meas * math.log2(n_meas))
    assert math.isclose(r350["linclust_cpu_min_linear"], cpu_meas * scale350, rel_tol=1e-9)
    assert math.isclose(r350["linclust_cpu_min_nlogn"], cpu_meas * nlogn350, rel_tol=1e-9)
    assert r350["linclust_cpu_min_nlogn"] > r350["linclust_cpu_min_linear"]
    assert math.isclose(r350["linclust_peak_rss_gib_linear"], (rss_meas / 1024**3) * scale350, rel_tol=1e-9)
    assert r350["label"] == "derived (extrapolation)"

    # Cascade N17 power law: 284.10 * (n / 2,053,075)^1.244
    expected_cas_1b = 284.10 * (1_000_000_000 / 2_053_075) ** 1.244
    assert math.isclose(r1b["cascade_cpu_min_n17"], expected_cas_1b, rel_tol=1e-9)
    assert r1b["speedup_vs_cascade_linear"] > 10.0
