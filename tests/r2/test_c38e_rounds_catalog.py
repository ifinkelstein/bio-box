"""Ring R2: card C-38e — §3.8 test (5) O(1) row-group pruning on merged round tables + old-manifest verification.

Covers:
  - §3.8 test (5): a release with ten hand-written rounds and merged `release=` copies of `round_domains`
    and `round_coverage` sorted by `protein_uid` reads O(1) row groups (1 of N) on a lookup by `protein_uid`
    (both literal and joined key) as measured by `catalog.profile_query` (AGENTS §4.4, LY-06).
  - Pre-C-38e published-release manifest compatibility: a manifest whose `tables` map predates C-38e
    (no `round_domains`, `round_coverage`, `release_round`) still verifies unchanged (`verify_manifest == []`,
    identical `manifest_hash`) after ten rounds are written into the same root.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog, ids, io
from biobox.registry.schemas import TABLES

pytestmark = [pytest.mark.r2, pytest.mark.covers("C-38e", "LY-06", "LY-10")]

LIB = "pfam-38.2"
N_ROUNDS = 10
PROTS_PER_ROUND = 200
RG_SIZE = 100  # 2,000 proteins -> 20 row groups in merged round_coverage / round_domains


def _puid(i: int) -> bytes:
    return ids.protein_uid(f"M{'K' * 12}{i:07d}")


def _write_minimal_stage(root: Path, stage: str = "S0") -> None:
    spec = TABLES["contig_index"]
    cu = ids.contig_uid(f"GATTACA{stage}".encode() * 30)
    rd = ids.raw_digest(f"GATTACA{stage}".encode() * 30)
    t = pa.table(
        {
            "contig_uid": pa.array([cu], pa.binary(16)),
            "raw_digest": pa.array([rd], pa.binary(32)),
            "contig_len": pa.array([len(f"GATTACA{stage}") * 30], pa.int64()),
            "is_default_variant": pa.array([True], pa.bool_()),
        },
        schema=spec.schema,
    )
    io.write_table(t, spec, root / spec.location / f"stage={stage}" / "part-0.parquet")


def _dom_table(uids: list[bytes]) -> pa.Table:
    spec = TABLES["round_domains"]
    n = len(uids)
    return pa.table(
        {
            "protein_uid": pa.array(uids, pa.binary(16)),
            "model_acc": pa.array([f"PF{i % 50:05d}.1" for i in range(n)], pa.string()),
            "model_name": pa.array([f"fam_{i % 50}" for i in range(n)], pa.string()),
            "model_len": pa.array([120] * n, pa.uint32()),
            "seq_bitscore": pa.array([60.0] * n, pa.float32()),
            "seq_evalue_zref": pa.array([1e-12] * n, pa.float64()),
            "dom_bitscore": pa.array([55.0] * n, pa.float32()),
            "dom_ievalue_zref": pa.array([1e-11] * n, pa.float64()),
            "env_start": pa.array([2] * n, pa.uint32()),
            "env_end": pa.array([15] * n, pa.uint32()),
            "ali_start": pa.array([2] * n, pa.uint32()),
            "ali_end": pa.array([15] * n, pa.uint32()),
            "hmm_from": pa.array([0] * n, pa.uint32()),
            "hmm_to": pa.array([13] * n, pa.uint32()),
            "cutoff": pa.array(["gathering"] * n, pa.string()),
            "z": pa.array([1e8] * n, pa.float64()),
            "dom_z": pa.array([1e8] * n, pa.float64()),
        },
        schema=spec.schema,
    )


def _cov_table(uids: list[bytes], n_rows: list[int]) -> pa.Table:
    spec = TABLES["round_coverage"]
    return pa.table(
        {
            "protein_uid": pa.array(uids, pa.binary(16)),
            "n_rows": pa.array(n_rows, pa.uint32()),
        },
        schema=spec.schema,
    )


@pytest.fixture(scope="module")
def ten_round_store(tmp_path_factory):
    d = tmp_path_factory.mktemp("c38e_10rounds")
    root = d / "root"
    _write_minimal_stage(root, "S0")

    # Release built before any round (simulate pre-C-38e manifest with no round tables in manifest["tables"])
    catalog.define_release(root, "rOld.1", ["S0"])
    m_old = catalog.build_manifest(root, "rOld.1", ["S0"])
    m_old_pre_c38e = copy.deepcopy(m_old)
    for tname in ("round_domains", "round_coverage", "release_round"):
        m_old_pre_c38e["tables"].pop(tname, None)
    catalog.write_manifest(m_old_pre_c38e)

    # Write 10 complete rounds
    dspec, cspec = TABLES["round_domains"], TABLES["round_coverage"]
    all_dom_tables = []
    all_cov_tables = []
    all_uids: list[bytes] = []
    for r in range(N_ROUNDS):
        rid = f"r-{r:016x}"
        uids = sorted(_puid(r * PROTS_PER_ROUND + j) for j in range(PROTS_PER_ROUND))
        all_uids.extend(uids)
        dt = _dom_table(uids)
        ct = _cov_table(uids, [1] * len(uids))
        io.write_table(dt, dspec, root / dspec.location / f"round={rid}" / f"library={LIB}" / "part-0.parquet")
        io.write_table(ct, cspec, root / cspec.location / f"round={rid}" / f"library={LIB}" / "part-0.parquet")
        all_dom_tables.append(dt)
        all_cov_tables.append(ct)

    # Define release taking all 10 rounds and write merged sorted release= copies with multiple row groups
    catalog.define_release(root, "rTen.1", ["S0"])
    merged_dom = pa.concat_tables(all_dom_tables).sort_by([(c, "ascending") for c in dspec.sort_key])
    merged_cov = pa.concat_tables(all_cov_tables).sort_by([(c, "ascending") for c in cspec.sort_key])
    dom_merged_path = root / dspec.location / "release=rTen.1" / f"library={LIB}" / "part-0.parquet"
    cov_merged_path = root / cspec.location / "release=rTen.1" / f"library={LIB}" / "part-0.parquet"
    io.write_table(merged_dom, dspec, dom_merged_path, row_group_rows=RG_SIZE)
    io.write_table(merged_cov, cspec, cov_merged_path, row_group_rows=RG_SIZE)

    man = catalog.build_manifest(root, "rTen.1", ["S0"])
    catalog.write_manifest(man)
    cat_path = catalog.generate_catalog(man, out_path=d / "rTen.1.duckdb", tape_version="t-test")
    return {
        "root": root,
        "m_old": m_old_pre_c38e,
        "m_ten": man,
        "cat_path": cat_path,
        "dom_merged_path": dom_merged_path,
        "cov_merged_path": cov_merged_path,
        "all_uids": sorted(all_uids),
    }


def test_c38e_published_pre_round_manifest_verifies_unchanged_with_ten_rounds(ten_round_store):
    root, m_old = ten_round_store["root"], ten_round_store["m_old"]
    assert "round_domains" not in m_old["tables"]
    assert catalog.verify_manifest(m_old, root) == []
    assert catalog.verify_manifest(ten_round_store["m_ten"], root) == []
    assert catalog.release_rounds(root, "rTen.1") == [f"r-{r:016x}" for r in range(N_ROUNDS)]


def test_c38e_lookup_by_protein_in_ten_round_release_reads_one_row_group(ten_round_store):
    """§3.8 test (5): `profile_query` of a lookup by `protein_uid` on the merged `release=` round files
    reads 1 of 20 row groups (both literal key and joined key)."""
    dom_rgs = pq.read_metadata(ten_round_store["dom_merged_path"]).num_row_groups
    cov_rgs = pq.read_metadata(ten_round_store["cov_merged_path"]).num_row_groups
    assert dom_rgs == (N_ROUNDS * PROTS_PER_ROUND) // RG_SIZE == 20
    assert cov_rgs == 20

    target = ten_round_store["all_uids"][len(ten_round_store["all_uids"]) // 2]
    con = catalog.open_catalog(ten_round_store["cat_path"])
    try:
        # 1. Literal lookup on merged round_domains
        p_dom = catalog.profile_query(
            con,
            "SELECT protein_uid, model_acc FROM round_domains "
            "WHERE release = catalog_release() AND library = ? AND protein_uid = ?",
            [LIB, target],
        )
        assert len(p_dom["result"]) == 1 and p_dom["result"][0][0] == target
        dom_scans = [s for s in p_dom["scans"] if any("release=rTen.1" in f for f in s["files"])]
        assert len(dom_scans) == 1
        assert dom_scans[0]["row_groups_total"] == 20
        assert dom_scans[0]["row_groups_scanned"] == 1

        # 2. Literal lookup on merged round_coverage
        p_cov = catalog.profile_query(
            con,
            "SELECT protein_uid, n_rows FROM round_coverage "
            "WHERE release = catalog_release() AND library = ? AND protein_uid = ?",
            [LIB, target],
        )
        assert p_cov["result"] == [(target, 1)]
        cov_scans = [s for s in p_cov["scans"] if any("release=rTen.1" in f for f in s["files"])]
        assert len(cov_scans) == 1
        assert cov_scans[0]["row_groups_total"] == 20
        assert cov_scans[0]["row_groups_scanned"] == 1

        # 3. Joined-key lookup on merged round_domains (AGENTS §4.4)
        p_join = catalog.profile_query(
            con,
            "WITH q AS (SELECT ?::BLOB AS protein_uid) "
            "SELECT d.protein_uid, d.model_acc FROM round_domains d "
            "JOIN q ON d.protein_uid = q.protein_uid "
            "WHERE d.release = catalog_release() AND d.library = ?",
            [target, LIB],
        )
        assert len(p_join["result"]) == 1 and p_join["result"][0][0] == target
        join_scans = [s for s in p_join["scans"] if any("release=rTen.1" in f for f in s["files"])]
        assert len(join_scans) == 1
        assert join_scans[0]["row_groups_total"] == 20
        assert join_scans[0]["row_groups_scanned"] == 1
    finally:
        con.close()
