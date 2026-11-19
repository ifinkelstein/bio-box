"""Ring R1: card C-38e — round tables, scope rule (e), release_round, and mixed-width ledger views (DEC-54).

Covers:
  - Registry TableSpecs (`round_domains`, `round_coverage`, `release_round`), `ledger` third key set
    `("round", "step")`, and `round_scan` ArtifactSpec (`in_manifest=False`, `small_ok=True`).
  - `define_release` writes `release_round` once beside `release_def` listing every complete round on disk
    when rounds exist (and writes none when 0 complete rounds exist, leaving pre-round releases untouched).
  - §3.8 test (4): scope rule (e) in `_scope`, `_in_release`, `ReleaseFiles`, `build_manifest`, and
    `verify_manifest`: a release built before any round ignores round files and keeps its manifest hash
    unchanged; a release listing rounds verifies; a round file no release lists is not reported as stray
    by any release's verify; a stray part inside a listed round IS reported by the release that lists it.
  - Handover from C-38d: `_view_sql` / `generate_catalog` tracks trailing nullable columns per
    `(table, partition_key_set)` branch so an 8-column pre-C-38d `stage=` ledger unions cleanly with
    10-column `release=` and `round=` ledgers (NULL on the trailing base columns for the 8-column rows).
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog, ids, io
from biobox.registry.schemas import ARTIFACTS, TABLES, conformance_problems

pytestmark = pytest.mark.covers("C-38e", "C-01", "LY-05", "LY-09", "LY-10")

LIB = "pfam-38.2"


def _uid(n: int) -> bytes:
    return ids.protein_uid(f"M{'A' * 10}{n:06d}")


def _write_minimal_stage(root: Path, stage: str = "S0") -> None:
    """Minimal stage files (`contig_index` + pre-C-38d or full `ledger`) so `build_manifest` accepts `stage`."""
    spec = TABLES["contig_index"]
    cu = ids.contig_uid(f"ACGTACGT{stage}".encode() * 20)
    rd = ids.raw_digest(f"ACGTACGT{stage}".encode() * 20)
    t = pa.table(
        {
            "contig_uid": pa.array([cu], pa.binary(16)),
            "raw_digest": pa.array([rd], pa.binary(32)),
            "contig_len": pa.array([160 + len(stage) * 20], pa.int64()),
            "is_default_variant": pa.array([True], pa.bool_()),
        },
        schema=spec.schema,
    )
    io.write_table(t, spec, root / spec.location / f"stage={stage}" / "part-0.parquet")


def _write_round(
    root: Path,
    round_id: str,
    protein_indices: list[int],
    *,
    library: str = LIB,
    write_domains: bool = True,
    write_coverage: bool = True,
    write_ledger: bool = True,
    write_scan_plan: bool = True,
) -> None:
    """Write hand-crafted round Parquet files and optional `round_scan` plan under `root`."""
    uids = sorted(_uid(i) for i in protein_indices)
    if write_domains:
        dspec = TABLES["round_domains"]
        # Give even-indexed proteins a hit; odd-indexed proteins have 0 rows (covered with n_rows=0)
        hit_uids = [u for idx, u in enumerate(uids) if idx % 2 == 0]
        dt = pa.table(
            {
                "protein_uid": pa.array(hit_uids, pa.binary(16)),
                "model_acc": pa.array([f"PF{i:05d}.1" for i in range(len(hit_uids))], pa.string()),
                "model_name": pa.array([f"dom_{i}" for i in range(len(hit_uids))], pa.string()),
                "model_len": pa.array([100] * len(hit_uids), pa.uint32()),
                "seq_bitscore": pa.array([50.0] * len(hit_uids), pa.float32()),
                "seq_evalue_zref": pa.array([1e-10] * len(hit_uids), pa.float64()),
                "dom_bitscore": pa.array([48.0] * len(hit_uids), pa.float32()),
                "dom_ievalue_zref": pa.array([1e-9] * len(hit_uids), pa.float64()),
                "env_start": pa.array([1] * len(hit_uids), pa.uint32()),
                "env_end": pa.array([10] * len(hit_uids), pa.uint32()),
                "ali_start": pa.array([1] * len(hit_uids), pa.uint32()),
                "ali_end": pa.array([10] * len(hit_uids), pa.uint32()),
                "hmm_from": pa.array([0] * len(hit_uids), pa.uint32()),
                "hmm_to": pa.array([9] * len(hit_uids), pa.uint32()),
                "cutoff": pa.array(["gathering"] * len(hit_uids), pa.string()),
                "z": pa.array([1e8] * len(hit_uids), pa.float64()),
                "dom_z": pa.array([1e8] * len(hit_uids), pa.float64()),
            },
            schema=dspec.schema,
        )
        io.write_table(dt, dspec, root / dspec.location / f"round={round_id}" / f"library={library}" / "part-0.parquet")
    if write_coverage:
        cspec = TABLES["round_coverage"]
        n_rows = [1 if idx % 2 == 0 else 0 for idx in range(len(uids))]
        ct = pa.table(
            {
                "protein_uid": pa.array(uids, pa.binary(16)),
                "n_rows": pa.array(n_rows, pa.uint32()),
            },
            schema=cspec.schema,
        )
        io.write_table(ct, cspec, root / cspec.location / f"round={round_id}" / f"library={library}" / "part-0.parquet")
    if write_ledger:
        lspec = TABLES["ledger"]
        lt = pa.table(
            {
                "object_key": pa.array([library.encode()], pa.binary()),
                "object_kind": pa.array(["library"], pa.string()),
                "n_in": pa.array([len(uids)], pa.uint64()),
                "n_kept": pa.array([sum(1 for idx in range(len(uids)) if idx % 2 == 0)], pa.uint64()),
                "n_merged": pa.array([0], pa.uint64()),
                "n_rejected": pa.array([sum(1 for idx in range(len(uids)) if idx % 2 == 1)], pa.uint64()),
                "reason_code": pa.array(["no_hit"], pa.string()),
                "tool_versions_json": pa.array(['{"pyhmmer":"0.12.0"}'], pa.string()),
                "n_bases_in": pa.array([None], pa.uint64()),
                "n_bases_rejected": pa.array([None], pa.uint64()),
            },
            schema=lspec.schema,
        )
        io.write_table(lt, lspec, root / lspec.location / f"round={round_id}" / f"step=domains-{library}" / "part-0.parquet")
    if write_scan_plan:
        plan_path = root / ARTIFACTS["round_scan"].location / f"round={round_id}" / f"library={library}" / "plan.json"
        plan_path.parent.mkdir(parents=True, exist_ok=True)
        plan_path.write_text(json.dumps({"round": round_id, "library": library, "n_requested": len(uids)}, sort_keys=True))


# --- (1) Registry specs & classification ---------------------------------------------------------------


def test_c38e_registry_specs_and_classification():
    rd = TABLES["round_domains"]
    assert rd.location == "l1/round_domains"
    assert rd.layer == "L1"
    assert rd.partition_keys == ("round", "library")
    assert rd.alt_partition_keys == (("release", "library"),)
    assert rd.primary_key == ("protein_uid", "model_acc", "env_start", "round", "library")
    assert rd.sort_key == ("protein_uid", "env_start", "model_acc")
    assert rd.schema.equals(TABLES["l1_domains"].schema, check_metadata=False)
    assert rd.enums == TABLES["l1_domains"].enums
    assert not rd.small_ok

    rc = TABLES["round_coverage"]
    assert rc.location == "l1/round_coverage"
    assert rc.layer == "L1"
    assert rc.partition_keys == ("round", "library")
    assert rc.alt_partition_keys == (("release", "library"),)
    assert rc.primary_key == ("protein_uid", "round", "library")
    assert rc.sort_key == ("protein_uid",)
    assert rc.schema.names == ["protein_uid", "n_rows"]
    assert rc.schema.field("protein_uid").type == pa.binary(16) and not rc.schema.field("protein_uid").nullable
    assert rc.schema.field("n_rows").type == pa.uint32() and not rc.schema.field("n_rows").nullable
    assert not rc.small_ok

    rr = TABLES["release_round"]
    assert rr.location == "release_round"
    assert rr.layer == "release"
    assert rr.partition_keys == ("release",)
    assert rr.primary_key == ("release", "round")
    assert rr.sort_key == ("position",)
    assert rr.schema.names == ["position", "round"]
    assert rr.schema.field("position").type == pa.uint16() and not rr.schema.field("position").nullable
    assert rr.schema.field("round").type == pa.string() and not rr.schema.field("round").nullable
    assert rr.small_ok

    led = TABLES["ledger"]
    assert ("round", "step") in led.alt_partition_keys

    rs = ARTIFACTS["round_scan"]
    assert rs.location == "l1/round_scan"
    assert rs.partition_keys == ("round", "library")
    assert not rs.in_manifest
    assert rs.small_ok

    # Classification & scope rule (e)
    assert catalog.classify("l1/round_domains/round=r-001/library=pfam-38.2/part-0.parquet") == (
        rd, {"round": "r-001", "library": LIB}
    )
    assert catalog.classify("l1/round_domains/release=r26.11.2/library=pfam-38.2/part-0.parquet") == (
        rd, {"release": "r26.11.2", "library": LIB}
    )
    assert catalog.classify("l1/round_coverage/round=r-001/library=pfam-38.2/part-0.parquet") == (
        rc, {"round": "r-001", "library": LIB}
    )
    assert catalog.classify("l1/round_coverage/release=r26.11.2/library=pfam-38.2/part-0.parquet") == (
        rc, {"release": "r26.11.2", "library": LIB}
    )
    assert catalog.classify("release_round/release=r26.11.2/part-0.parquet") == (rr, {"release": "r26.11.2"})
    assert catalog.classify("ledger/round=r-001/step=domains-pfam-38.2/part-0.parquet") == (
        led, {"round": "r-001", "step": f"domains-{LIB}"}
    )
    assert catalog.classify_artifact("l1/round_scan/round=r-001/library=pfam-38.2/plan.json") == (
        rs, {"round": "r-001", "library": LIB}, ""
    )
    assert catalog.classify_artifact("l1/round_scan/round=r-001/library=pfam-38.2/shard-00000.json") == (
        rs, {"round": "r-001", "library": LIB}, ""
    )
    assert catalog._scope("round_domains", {"round": "r-001", "library": LIB}) == ("round", "r-001")
    assert catalog._scope("round_domains", {"release": "r26.11.2", "library": LIB}) == ("release", "r26.11.2")
    assert catalog._scope("ledger", {"round": "r-001", "step": f"domains-{LIB}"}) == ("round", "r-001")


# --- (2) define_release & complete_rounds --------------------------------------------------------------


def test_c38e_define_release_writes_release_round_only_when_complete_rounds_exist(tmp_path):
    root = tmp_path / "root"
    _write_minimal_stage(root, "S0")

    # 0 rounds on disk: define_release writes release_def only, not release_round
    assert catalog.complete_rounds(root) == []
    catalog.define_release(root, "rPre.1", ["S0"])
    assert catalog.release_def_path(root, "rPre.1").is_file()
    assert not catalog.release_round_path(root, "rPre.1").exists()
    assert catalog.release_rounds(root, "rPre.1") == []

    # Incomplete rounds (plan only, domains only, coverage only) are NOT complete
    _write_round(root, "r-plan-only", [1], write_domains=False, write_coverage=False, write_ledger=False, write_scan_plan=True)
    _write_round(root, "r-dom-only", [2], write_domains=True, write_coverage=False, write_ledger=False)
    _write_round(root, "r-cov-only", [3], write_domains=False, write_coverage=True, write_ledger=False)
    assert catalog.complete_rounds(root) == []

    # Re-calling define_release on already-defined rPre.1 is a no-op even after complete rounds appear
    _write_round(root, "r-002", [10, 11])
    _write_round(root, "r-001", [20, 21])
    assert catalog.complete_rounds(root) == ["r-001", "r-002"]
    catalog.define_release(root, "rPre.1", ["S0"])
    assert not catalog.release_round_path(root, "rPre.1").exists()

    # New release defined after r-001 and r-002 completed: writes release_round once
    catalog.define_release(root, "rWith.1", ["S0"])
    rr_path = catalog.release_round_path(root, "rWith.1")
    assert rr_path.is_file()
    t = pq.read_table(rr_path, partitioning=None)
    assert conformance_problems(t, TABLES["release_round"]) == []
    assert t.to_pylist() == [{"position": 0, "round": "r-001"}, {"position": 1, "round": "r-002"}]
    assert catalog.release_rounds(root, "rWith.1") == ["r-001", "r-002"]

    # Adding r-003 later does NOT alter rWith.1's release_round when define_release is called again
    before = (rr_path.read_bytes(), rr_path.stat().st_mtime_ns)
    _write_round(root, "r-003", [30, 31])
    assert catalog.complete_rounds(root) == ["r-001", "r-002", "r-003"]
    catalog.define_release(root, "rWith.1", ["S0"])
    assert (rr_path.read_bytes(), rr_path.stat().st_mtime_ns) == before
    assert catalog.release_rounds(root, "rWith.1") == ["r-001", "r-002"]

    # Staging `out` receives release_round beside release_def when defining a fresh release
    staging = tmp_path / "staging"
    catalog.define_release(root, "rStage.1", ["S0"], out=staging)
    assert catalog.release_def_path(staging, "rStage.1").is_file()
    assert catalog.release_round_path(staging, "rStage.1").is_file()
    assert not catalog.release_round_path(root, "rStage.1").exists()
    assert catalog.release_rounds(staging, "rStage.1") == ["r-001", "r-002", "r-003"]


# --- (3) §3.8 test (4): scope rule (e), verify_manifest, ReleaseFiles ---------------------------------


def test_c38e_scope_rule_e_verify_and_release_files(tmp_path):
    root = tmp_path / "root"
    _write_minimal_stage(root, "S0")

    # 1. Release R_PRE built before any round exists
    catalog.define_release(root, "rPre.1", ["S0"])
    m_pre = catalog.build_manifest(root, "rPre.1", ["S0"])
    catalog.write_manifest(m_pre)
    h_pre = catalog.manifest_hash(m_pre)
    assert catalog.verify_manifest(m_pre, root) == []

    # 2. Rounds r-001 and r-002 arrive; R_PRE still verifies with hash unchanged and ignores them
    _write_round(root, "r-001", [1, 2, 3, 4])
    _write_round(root, "r-002", [5, 6, 7, 8])
    assert catalog.verify_manifest(m_pre, root) == []
    assert catalog.manifest_hash(m_pre) == h_pre
    assert catalog.ReleaseFiles.of(root, "rPre.1").entries("round_domains") == []
    assert catalog.ReleaseFiles(root, "rPre.1", ["S0"]).entries("round_domains") == []

    # 3. Release R_WITH defined after r-001 and r-002; lists both rounds + merged release= files
    catalog.define_release(root, "rWith.1", ["S0"])
    # Also write merged release= files (alt_partition_keys)
    dspec, cspec = TABLES["round_domains"], TABLES["round_coverage"]
    merged_dom = pa.concat_tables([
        pq.read_table(root / dspec.location / f"round={r}" / f"library={LIB}" / "part-0.parquet")
        for r in ("r-001", "r-002")
    ]).sort_by([(c, "ascending") for c in dspec.sort_key])
    merged_cov = pa.concat_tables([
        pq.read_table(root / cspec.location / f"round={r}" / f"library={LIB}" / "part-0.parquet")
        for r in ("r-001", "r-002")
    ]).sort_by([(c, "ascending") for c in cspec.sort_key])
    io.write_table(merged_dom, dspec, root / dspec.location / "release=rWith.1" / f"library={LIB}" / "part-0.parquet")
    io.write_table(merged_cov, cspec, root / cspec.location / "release=rWith.1" / f"library={LIB}" / "part-0.parquet")

    m_with = catalog.build_manifest(root, "rWith.1", ["S0"])
    catalog.write_manifest(m_with)
    assert catalog.verify_manifest(m_with, root) == []
    assert catalog.verify_manifest(m_pre, root) == []

    listed_tables = {f["table"] for f in m_with["files"]}
    assert {"release_round", "round_domains", "round_coverage", "ledger"} <= listed_tables
    listed_rounds = {f["partitions"]["round"] for f in m_with["files"] if "round" in f["partitions"]}
    assert listed_rounds == {"r-001", "r-002"}
    # round_scan has in_manifest=False: never listed as an artifact
    assert not any(a["artifact"] == "round_scan" for a in m_with["artifacts"])

    # 4. Round r-003 (complete, unlisted by any release) and r-unlisted (domains only, unlisted) arrive
    _write_round(root, "r-003", [100, 101])
    _write_round(root, "r-unlisted", [200], write_coverage=False, write_ledger=False)
    assert catalog.verify_manifest(m_pre, root) == []
    assert catalog.verify_manifest(m_with, root) == []

    # ReleaseFiles without manifest (from release_def + release_round) never opens r-003 or r-unlisted
    rf_no_man = catalog.ReleaseFiles(root, "rWith.1", ["S0"])
    got_rounds = {parts["round"] for _, parts in rf_no_man.entries("round_domains") if "round" in parts}
    assert got_rounds == {"r-001", "r-002"}
    assert rf_no_man.entries("round_domains", round="r-003") == []
    got_merged = [p for p, parts in rf_no_man.entries("round_domains") if parts.get("release") == "rWith.1"]
    assert len(got_merged) == 1

    # 5. Stray part inside a listed round (r-001): reported by R_WITH, ignored by R_PRE
    src = root / f"l1/round_domains/round=r-001/library={LIB}/part-0.parquet"
    stray = root / f"l1/round_domains/round=r-001/library={LIB}/part-9.parquet"
    shutil.copy(src, stray)
    assert catalog.verify_manifest(m_pre, root) == []
    assert catalog.verify_manifest(m_with, root) == [
        f"not in manifest: l1/round_domains/round=r-001/library={LIB}/part-9.parquet"
    ]
    stray.unlink()

    # 6. Unclaimed path under an unlisted round: reported by BOTH releases
    unclaimed = root / f"l1/round_domains/round=r-003/extra/part-0.parquet"
    unclaimed.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(src, unclaimed)
    want_unclaimed = ["not in manifest: l1/round_domains/round=r-003/extra/part-0.parquet"]
    assert catalog.verify_manifest(m_pre, root) == want_unclaimed
    assert catalog.verify_manifest(m_with, root) == want_unclaimed


# --- (4) Handover from C-38d: 8-col stage= ledger + 10-col release= / round= ledgers -------------------


def test_c38e_catalog_ledger_unions_8col_stage_with_10col_release_and_round(tmp_path):
    """Pre-C-38d 8-column `stage=` ledger files beside 10-column `release=` and `round=` ledger files must
    generate a catalog whose `ledger` view queries cleanly with NULL trailing base columns on the 8-col rows."""
    root = tmp_path / "root"
    _write_minimal_stage(root, "S0")

    full_spec = TABLES["ledger"]
    old_schema = pa.schema(list(full_spec.schema)[:8])
    assert len(old_schema) == 8 and len(full_spec.schema) == 10

    # 8-column pre-C-38d stage= ledger
    t_stage = pa.table(
        {
            "object_key": pa.array([b"obj_stage"], pa.binary()),
            "object_kind": pa.array(["l0_object"], pa.string()),
            "n_in": pa.array([10], pa.uint64()),
            "n_kept": pa.array([10], pa.uint64()),
            "n_merged": pa.array([0], pa.uint64()),
            "n_rejected": pa.array([0], pa.uint64()),
            "reason_code": pa.array(["kept"], pa.string()),
            "tool_versions_json": pa.array(['{"biobox":"0.1"}'], pa.string()),
        },
        schema=old_schema,
    )
    p_stage = root / "ledger/stage=S0/step=ingest/part-0.parquet"
    p_stage.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(t_stage, p_stage)

    # 10-column round= ledger (via a complete round)
    _write_round(root, "r-led1", [1, 2], write_ledger=False)
    # Write the round ledger with explicit non-null n_bases_in / n_bases_rejected to verify propagation
    t_round = pa.table(
        {
            "object_key": pa.array([b"pfam-38.2"], pa.binary()),
            "object_kind": pa.array(["library"], pa.string()),
            "n_in": pa.array([2], pa.uint64()),
            "n_kept": pa.array([1], pa.uint64()),
            "n_merged": pa.array([0], pa.uint64()),
            "n_rejected": pa.array([1], pa.uint64()),
            "reason_code": pa.array(["no_hit"], pa.string()),
            "tool_versions_json": pa.array(['{"pyhmmer":"0.12.0"}'], pa.string()),
            "n_bases_in": pa.array([500], pa.uint64()),
            "n_bases_rejected": pa.array([250], pa.uint64()),
        },
        schema=full_spec.schema,
    )
    io.write_table(t_round, full_spec, root / f"ledger/round=r-led1/step=domains-{LIB}/part-0.parquet")

    # 10-column release= ledger
    catalog.define_release(root, "rMix.1", ["S0"])
    t_rel = pa.table(
        {
            "object_key": pa.array([b"nr90"], pa.binary()),
            "object_kind": pa.array(["level"], pa.string()),
            "n_in": pa.array([5], pa.uint64()),
            "n_kept": pa.array([5], pa.uint64()),
            "n_merged": pa.array([0], pa.uint64()),
            "n_rejected": pa.array([0], pa.uint64()),
            "reason_code": pa.array(["kept"], pa.string()),
            "tool_versions_json": pa.array(['{"mmseqs":"15"}'], pa.string()),
            "n_bases_in": pa.array([1234], pa.uint64()),
            "n_bases_rejected": pa.array([0], pa.uint64()),
        },
        schema=full_spec.schema,
    )
    io.write_table(t_rel, full_spec, root / "ledger/release=rMix.1/step=cluster/part-0.parquet")

    man = catalog.build_manifest(root, "rMix.1", ["S0"])
    n_cols_by_path = {f["path"]: f["n_columns"] for f in man["files"] if f["table"] == "ledger"}
    assert n_cols_by_path == {
        "ledger/stage=S0/step=ingest/part-0.parquet": 8,
        "ledger/release=rMix.1/step=cluster/part-0.parquet": 10,
        f"ledger/round=r-led1/step=domains-{LIB}/part-0.parquet": 10,
    }

    cat_path = catalog.generate_catalog(man, out_path=tmp_path / "rMix.1.duckdb", tape_version="t-test")
    con = catalog.open_catalog(cat_path)
    try:
        rows = con.execute(
            "SELECT stage, release, round, step, n_in, n_bases_in, n_bases_rejected "
            "FROM ledger ORDER BY step"
        ).fetchall()
    finally:
        con.close()

    assert rows == [
        (None, "rMix.1", None, "cluster", 5, 1234, 0),
        (None, None, "r-led1", f"domains-{LIB}", 2, 500, 250),
        ("S0", None, None, "ingest", 10, None, None),
    ]
