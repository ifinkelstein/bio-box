"""Ring R1: unit tests for retiring old releases' large release= files (DEC-66 S-4, card C-54)."""

from __future__ import annotations

import hashlib
import json
import subprocess
import tarfile
from pathlib import Path

import pyarrow as pa
import pytest

from biobox import catalog, cli, contracts, ids, io, ledger
from biobox.l2 import cluster
from biobox.registry.schemas import ARTIFACTS, TABLES

pytestmark = [pytest.mark.covers("C-54", "LY-10")]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_minimal_stage(root: Path, stage: str = "S0") -> None:
    """Write minimal valid stage Parquet tables (contig_index, contig_tapes, protein_seq, ledger)."""
    cuid = b"\x01" * 16
    puid = b"\x02" * 16
    fuid = b"\x03" * 16
    rd = b"\x04" * 32

    ci_spec = TABLES["contig_index"]
    io.write_table(
        pa.table(
            {
                "contig_uid": pa.array([cuid], pa.binary(16)),
                "raw_digest": pa.array([rd], pa.binary(32)),
                "contig_len": pa.array([1000], pa.int64()),
                "is_default_variant": pa.array([True], pa.bool_()),
            },
            schema=ci_spec.schema,
        ),
        ci_spec,
        root / ci_spec.location / f"stage={stage}" / "part-0.parquet",
    )

    ps_spec = TABLES["protein_seq"]
    io.write_table(
        pa.table(
            {
                "protein_uid": pa.array([puid], pa.binary(16)),
                "protein_digest": pa.array([rd], pa.binary(32)),
                "length": pa.array([10], pa.uint32()),
                "seq": pa.array(["MAAAAAAAAA"], pa.string()),
            },
            schema=ps_spec.schema,
        ),
        ps_spec,
        root / ps_spec.location / "gene_call_version=pyrodigal-3.6.3" / f"stage={stage}" / "protein_prefix=all" / "part-0.parquet",
    )

    tv = "t-0123456789ab"
    tv_spec = TABLES["tape_version_def"]
    tv_path = root / tv_spec.location / f"tape_version={tv}" / "part-0.parquet"
    if not tv_path.exists():
        io.write_table(
            pa.table({"producer": pa.array(["pyrodigal-3.6.3"], pa.string())}, schema=tv_spec.schema),
            tv_spec,
            tv_path,
        )

    ct_spec = TABLES["contig_tapes"]
    io.write_table(
        pa.table(
            {
                "contig_uid": pa.array([cuid], pa.binary(16)),
                "contig_len": pa.array([1000], pa.int64()),
                "is_circular": pa.array([False], pa.bool_()),
                "n_slots": pa.array([1], pa.uint32()),
                "feature_kinds": pa.array([[1]], pa.list_(pa.uint8())),
                "strands": pa.array([[1]], pa.list_(pa.int8())),
                "start_bps": pa.array([[0]], pa.list_(pa.uint32())),
                "end_bps": pa.array([[30]], pa.list_(pa.uint32())),
                "delta_bp_prev_slots": pa.array([[0]], pa.list_(pa.int32())),
                "junction_topology_prevs": pa.array([[0]], pa.list_(pa.uint8())),
                "orfs_to_contig_starts": pa.array([[0]], pa.list_(pa.uint16())),
                "orfs_to_contig_ends": pa.array([[0]], pa.list_(pa.uint16())),
                "slot_flags": pa.array([[0]], pa.list_(pa.uint8())),
                "protein_uids": pa.array([[puid]], pa.list_(pa.binary(16))),
                "feature_uids": pa.array([[fuid]], pa.list_(pa.binary(16))),
            },
            schema=ct_spec.schema,
        ),
        ct_spec,
        root / ct_spec.location / f"tape_version={tv}" / f"stage={stage}" / "part-0.parquet",
    )


def _write_minimal_release(root: Path, release: str, stages: list[str] = ["S0"], with_plan: bool = True) -> dict:
    """Write minimal release= tables, ledger, release_def, cluster_plan, manifest, and catalog."""
    cuid = b"\x01" * 16
    puid = b"\x02" * 16
    fuid = b"\x03" * 16
    guid = b"\x05" * 16
    cluid = b"\x06" * 16
    tv = "t-0123456789ab"

    catalog.define_release(root, release, stages)

    # cluster_member
    cm_spec = TABLES["cluster_member"]
    io.write_table(
        pa.table(
            {
                "member_uid": pa.array([puid], pa.binary(16)),
                "level": pa.array(["nr90"], pa.string()),
                "cluster_uid": pa.array([cluid], pa.binary(16)),
                "rep_uid": pa.array([puid], pa.binary(16)),
                "cluster_params_hash": pa.array(["abcd" * 4], pa.string()),
            },
            schema=cm_spec.schema,
        ),
        cm_spec,
        root / cm_spec.location / f"release={release}" / "part-0.parquet",
    )

    # member_cluster
    mc_spec = TABLES["member_cluster"]
    io.write_table(
        pa.table(
            {
                "member_uid": pa.array([puid], pa.binary(16)),
                "level": pa.array(["nr90"], pa.string()),
                "cluster_uid": pa.array([cluid], pa.binary(16)),
                "rep_uid": pa.array([puid], pa.binary(16)),
            },
            schema=mc_spec.schema,
        ),
        mc_spec,
        root / mc_spec.location / f"release={release}" / "part-0.parquet",
    )

    # token_dict
    td_spec = TABLES["token_dict"]
    io.write_table(
        pa.table(
            {
                "token_id": pa.array([1], pa.uint64()),
                "level": pa.array(["nr90"], pa.string()),
                "source": pa.array(["mmseqs2"], pa.string()),
                "key": pa.array([cluid.hex()], pa.string()),
                "n_slots": pa.array([1], pa.uint64()),
            },
            schema=td_spec.schema,
        ),
        td_spec,
        root / td_spec.location / f"release={release}" / "part-0.parquet",
    )

    # tape_tokens
    tt_spec = TABLES["tape_tokens"]
    io.write_table(
        pa.table(
            {
                "contig_uid": pa.array([cuid], pa.binary(16)),
                "family_tokens": pa.array([[0]], pa.list_(pa.uint64())),
                "nr50_tokens": pa.array([[1]], pa.list_(pa.uint64())),
                "nr90_tokens": pa.array([[1]], pa.list_(pa.uint64())),
                "domain_family_tokens": pa.array([[[]]], pa.list_(pa.list_(pa.uint64()))),
            },
            schema=tt_spec.schema,
        ),
        tt_spec,
        root / tt_spec.location / f"tape_version={tv}" / f"release={release}" / "part-0.parquet",
    )

    # token_postings
    tp_spec = TABLES["token_postings"]
    io.write_table(
        pa.table(
            {
                "token_id": pa.array([1], pa.uint64()),
                "contig_uid": pa.array([cuid], pa.binary(16)),
                "ordinal": pa.array([0], pa.uint32()),
                "seg_idx": pa.array([0], pa.uint16()),
                "feature_uid": pa.array([fuid], pa.binary(16)),
                "start": pa.array([0], pa.int64()),
                "end": pa.array([30], pa.int64()),
                "strand": pa.array([1], pa.int8()),
                "n_slots": pa.array([1], pa.uint32()),
                "is_circular": pa.array([False], pa.bool_()),
                "contig_len": pa.array([1000], pa.int64()),
            },
            schema=tp_spec.schema,
        ),
        tp_spec,
        root / tp_spec.location / f"release={release}" / f"tape_version={tv}" / "token_bucket=all" / "part-0.parquet",
    )

    # contig_dim
    cd_spec = TABLES["contig_dim"]
    io.write_table(
        pa.table(
            {
                "contig_uid": pa.array([cuid], pa.binary(16)),
                "genome_uid": pa.array([guid], pa.binary(16)),
                "contig_len": pa.array([1000], pa.int64()),
                "is_circular": pa.array([False], pa.bool_()),
                "tier": pa.array([None], pa.uint8()),
                "assembly_quality": pa.array(["complete"], pa.string()),
                "tax_order": pa.array([None], pa.uint32()),
                "is_species_rep": pa.array([True], pa.bool_()),
                "n_genomes": pa.array([1], pa.uint32()),
            },
            schema=cd_spec.schema,
        ),
        cd_spec,
        root / cd_spec.location / f"release={release}" / "part-0.parquet",
    )

    # genome_dim
    gd_spec = TABLES["genome_dim"]
    io.write_table(
        pa.table(
            {
                "genome_uid": pa.array([guid], pa.binary(16)),
                "n_contigs": pa.array([1], pa.uint32()),
                "genome_len": pa.array([1000], pa.int64()),
                "kind": pa.array(["virus"], pa.string()),
                "genetic_code": pa.array([11], pa.uint8()),
                "assembly_quality": pa.array(["complete"], pa.string()),
                "tier": pa.array([None], pa.uint8()),
                "tax_order": pa.array([None], pa.uint32()),
                "is_species_rep": pa.array([True], pa.bool_()),
            },
            schema=gd_spec.schema,
        ),
        gd_spec,
        root / gd_spec.location / f"release={release}" / "part-0.parquet",
    )

    # release ledger (must STAY on disk)
    ledger.write_release(
        root,
        release,
        "cluster",
        [
            {
                "object_key": b"nr90",
                "object_kind": "level",
                "n_in": 1,
                "n_kept": 1,
                "n_merged": 0,
                "n_rejected": 0,
                "reason_code": "ok",
                "tool_versions_json": "{}",
            }
        ],
    )

    if with_plan:
        pdir = root / "l2" / "cluster_plan" / f"release={release}"
        pdir.mkdir(parents=True, exist_ok=True)
        (pdir / "cluster.plan.json").write_text(json.dumps({"release": release}, sort_keys=True) + "\n")

    man = catalog.build_manifest(root, release, stages)
    catalog.write_manifest(man)
    catalog.generate_catalog(man)
    return man


def _make_backup_tar(root: Path, release: str, tar_path: Path) -> str:
    """Create a real tar of the release's retirable files and return its sha256."""
    man = catalog.read_manifest(catalog.manifest_path(root, release))
    with tarfile.open(tar_path, "w") as tf:
        for f in man["files"]:
            if f["partitions"].get("release") == release:
                tf.add(root / f["path"], arcname=f["path"])
        pdir = root / "l2" / "cluster_plan" / f"release={release}"
        if pdir.is_dir():
            for p in sorted(pdir.glob("*.json")):
                tf.add(p, arcname=p.relative_to(root).as_posix())
    return _sha256(tar_path)


def _copy_arg(sha: str, name: str = "backup.tar") -> str:
    """A `--copy <LOCATION>=<SHA256>` value (card R-04): the tower copy of a tar, with the tar's sha256."""
    return f"tower:/mnt/user/research-bulk/bio-box/root_backups/{name}={sha}"


def test_c54_release_archived_artifact_registered():
    assert "release_archived" in ARTIFACTS
    spec = ARTIFACTS["release_archived"]
    assert spec.location == "archived"
    assert spec.partition_keys == ()
    assert spec.pattern == r"release=[^/=]+\.json"
    assert spec.in_manifest is False
    assert spec.small_ok is True
    assert contracts.check_file_sizes(
        {
            "files": [],
            "artifacts": [
                {"artifact": "release_archived", "partitions": {}, "path": "archived/release=r26.11.1.json", "bytes": 100},
                {"artifact": "release_archived", "partitions": {}, "path": "archived/release=r26.11.2.json", "bytes": 100},
            ],
        }
    ) == []


def test_c54_retire_refusals_and_keep_latest_policy(tmp_path, capsys):
    root = tmp_path / "root"
    _write_minimal_stage(root, "S0")
    for rel in ("rS0.1", "rN.3", "rM.2", "r26.11.1", "r26.11.2", "r26.11.3"):
        _write_minimal_release(root, rel)

    tar_path = tmp_path / "backup.tar"
    tar_sha = _make_backup_tar(root, "r26.11.1", tar_path)
    tar_arg = f"{tar_path}={tar_sha}"
    cp = ["--copy", _copy_arg(tar_sha)]

    # Refuse latest
    assert cli.main(["release", "retire", "--root", str(root), "--release", "latest", "--tar", tar_arg, *cp]) == 2
    # Refuse signed milestone rS0.1
    assert cli.main(["release", "retire", "--root", str(root), "--release", "rS0.1", "--tar", tar_arg, *cp]) == 2
    assert not (root / "archived" / "release=rS0.1.json").exists()
    # Refuse missing manifest
    assert cli.main(["release", "retire", "--root", str(root), "--release", "r26.11.99", "--tar", tar_arg, *cp]) == 2
    # Refuse --keep-latest 0
    assert cli.main(["release", "retire", "--root", str(root), "--release", "r26.11.1", "--keep-latest", "0", "--tar", tar_arg, *cp]) == 2
    # Refuse empty --tar, malformed --tar, missing tar file, wrong tar sha256
    assert cli.main(["release", "retire", "--root", str(root), "--release", "r26.11.1", *cp]) == 2
    assert cli.main(["release", "retire", "--root", str(root), "--release", "r26.11.1", "--tar", "not-valid", *cp]) == 2
    assert cli.main(["release", "retire", "--root", str(root), "--release", "r26.11.1", "--tar", f"{tmp_path / 'nope.tar'}={tar_sha}", *cp]) == 2
    assert cli.main(["release", "retire", "--root", str(root), "--release", "r26.11.1", "--tar", f"{tar_path}={'0' * 64}",
                     "--copy", _copy_arg("0" * 64)]) == 2
    assert not (root / "archived" / "release=r26.11.1.json").exists()

    # Refuse newest 2 rolling releases (r26.11.2, r26.11.3) under default --keep-latest 2
    assert cli.main(["release", "retire", "--root", str(root), "--release", "r26.11.2", "--tar", tar_arg, *cp]) == 2
    assert cli.main(["release", "retire", "--root", str(root), "--release", "r26.11.3", "--tar", tar_arg, *cp]) == 2
    assert not (root / "archived" / "release=r26.11.2.json").exists()
    assert not (root / "archived" / "release=r26.11.3.json").exists()

    # Old-scheme releases and r26.11.1 are eligible with --keep-latest 2
    for rel in ("rN.3", "rM.2", "r26.11.1"):
        t = tmp_path / f"{rel}.tar"
        h = _make_backup_tar(root, rel, t)
        assert cli.main(["release", "retire", "--root", str(root), "--release", rel, "--tar", f"{t}={h}",
                         "--copy", _copy_arg(h, f"{rel}.tar")]) == 0
        assert (root / "archived" / f"release={rel}.json").is_file()

    # Now only r26.11.2 and r26.11.3 remain non-archived rolling: retiring r26.11.2 with --keep-latest 2 MUST be refused
    t2 = tmp_path / "r26.11.2.tar"
    h2 = _make_backup_tar(root, "r26.11.2", t2)
    cp2 = ["--copy", _copy_arg(h2, "r26.11.2.tar")]
    assert cli.main(["release", "retire", "--root", str(root), "--release", "r26.11.2", "--keep-latest", "2", "--tar", f"{t2}={h2}", *cp2]) == 2
    assert not (root / "archived" / "release=r26.11.2.json").exists()

    # With --keep-latest 1, retiring r26.11.2 succeeds (leaving r26.11.3)
    assert cli.main(["release", "retire", "--root", str(root), "--release", "r26.11.2", "--keep-latest", "1", "--tar", f"{t2}={h2}", *cp2]) == 0
    assert (root / "archived" / "release=r26.11.2.json").is_file()
    capsys.readouterr()


def test_c54_retire_idempotency_and_tamper_detection(tmp_path, capsys):
    root = tmp_path / "root"
    _write_minimal_stage(root, "S0")
    for rel in ("r26.11.1", "r26.11.2", "r26.11.3"):
        _write_minimal_release(root, rel)

    # Set up a minimal host profile and dummy search pack for r26.11.1
    prof_tmp = tmp_path / "prof_tmp"
    pack_dir = prof_tmp / "packs" / "release=r26.11.1" / ("a" * 32)
    pack_dir.mkdir(parents=True)
    (pack_dir / "pack_manifest.json").write_text("{}\n")
    prof_path = tmp_path / "profile.yaml"
    prof_path.write_text(
        f"max_threads_day: 2\nmax_threads_night: 2\nday_start_hour: 8\nday_end_hour: 20\n"
        f"mem_scope: 8G\nrequire_mem_cap: false\nnice_monolithic: 19\nnice_fanout: 15\n"
        f"cpu_min_per_process: 300\nsmall_file_floor_bytes: 1000\ntmp: {prof_tmp}\n"
        f"durable: {tmp_path / 'dur'}\nbulk: {tmp_path / 'bulk'}\n"
    )

    tar_path = tmp_path / "r26.11.1.tar"
    tar_sha = _make_backup_tar(root, "r26.11.1", tar_path)
    tar_arg = f"{tar_path}={tar_sha}"
    copy_arg = f"tower:/mnt/user/research-bulk/bio-box/root_backups/r26.11.1.tar={tar_sha}"

    man_before = catalog.read_manifest(catalog.manifest_path(root, "r26.11.1"))
    h_before = catalog.manifest_hash(man_before)

    # Retire
    code = cli.main([
        "release", "retire", "--root", str(root), "--release", "r26.11.1",
        "--tar", tar_arg, "--copy", copy_arg, "--profile", str(prof_path),
    ])
    assert code == 0

    # Check archived/release=r26.11.1.json
    arch_path = root / "archived" / "release=r26.11.1.json"
    assert arch_path.is_file()
    raw = arch_path.read_text()
    rec = json.loads(raw)
    assert raw == json.dumps(rec, indent=2, sort_keys=True) + "\n"
    assert rec["schema_version"] == 1
    assert rec["release"] == "r26.11.1"
    assert rec["manifest_hash"] == h_before
    assert rec["tars"] == [{"path": str(tar_path), "sha256": tar_sha}]
    # card R-04: a copy is recorded as its location and sha256, the operator's statement
    assert rec["copies"] == [{
        "location": "tower:/mnt/user/research-bulk/bio-box/root_backups/r26.11.1.tar",
        "sha256": tar_sha,
        "stated_by": "operator",
    }]
    assert len(rec["removed_files"]) >= 7
    assert len(rec["removed_plan_files"]) == 1

    # Check kept files stay on disk
    assert (root / "manifests" / "release=r26.11.1.json").is_file()
    assert (root / "catalog" / "release=r26.11.1.duckdb").is_file()
    assert (root / "release_def" / "release=r26.11.1" / "part-0.parquet").is_file()
    assert (root / "ledger" / "release=r26.11.1" / "step=cluster" / "part-0.parquet").is_file()
    assert (root / "l0seq" / "contig_index" / "stage=S0" / "part-0.parquet").is_file()

    # Check retired files and search pack are gone
    for entry in rec["removed_files"] + rec["removed_plan_files"]:
        assert not (root / entry["path"]).exists()
    assert not (prof_tmp / "packs" / "release=r26.11.1").exists()

    # Idempotent re-run with same tar
    assert cli.main(["release", "retire", "--root", str(root), "--release", "r26.11.1", "--tar", tar_arg, "--copy", copy_arg]) == 0

    # Conflict with different tar
    other_tar = tmp_path / "other.tar"
    other_tar.write_bytes(b"different")
    other_sha = _sha256(other_tar)
    assert cli.main(["release", "retire", "--root", str(root), "--release", "r26.11.1", "--tar", f"{other_tar}={other_sha}",
                     "--copy", _copy_arg(other_sha, "other.tar")]) == 2

    # Crash recovery: straggler file still on disk -> verify fails, re-run retire cleans it up -> verify passes
    straggler = root / rec["removed_files"][0]["path"]
    straggler.parent.mkdir(parents=True, exist_ok=True)
    straggler.write_bytes(b"leftover")
    assert catalog.verify_manifest(man_before, root)
    assert cli.main(["release", "verify", "--root", str(root), "--release", "r26.11.1"]) == 1
    assert cli.main(["release", "retire", "--root", str(root), "--release", "r26.11.1", "--tar", tar_arg, "--copy", copy_arg]) == 0
    assert not straggler.exists()
    assert catalog.verify_manifest(man_before, root) == []

    # Tamper detection on archived/release=r26.11.1.json
    # 1. Wrong manifest_hash
    arch_path.write_text(json.dumps({**rec, "manifest_hash": "0" * 32}, indent=2, sort_keys=True) + "\n")
    assert catalog.verify_manifest(man_before, root)
    assert cli.main(["release", "verify", "--root", str(root), "--release", "r26.11.1"]) == 1

    # 2. Forged sha256 in removed_files
    bad_files = [{**rec["removed_files"][0], "sha256": "f" * 64}, *rec["removed_files"][1:]]
    arch_path.write_text(json.dumps({**rec, "removed_files": bad_files}, indent=2, sort_keys=True) + "\n")
    assert catalog.verify_manifest(man_before, root)

    # 3. Listing a stage= file in removed_files
    stage_f = next(f for f in man_before["files"] if "stage" in f["partitions"])
    bad_stage = [*rec["removed_files"], {"path": stage_f["path"], "size_bytes": stage_f["bytes"], "sha256": stage_f["sha256"]}]
    arch_path.write_text(json.dumps({**rec, "removed_files": bad_stage}, indent=2, sort_keys=True) + "\n")
    assert catalog.verify_manifest(man_before, root)

    # 4. Restore valid record, then delete an unarchived file (release_def) -> verify fails
    arch_path.write_text(raw)
    assert catalog.verify_manifest(man_before, root) == []
    rdef = root / "release_def" / "release=r26.11.1" / "part-0.parquet"
    rdef_bytes = rdef.read_bytes()
    rdef.unlink()
    assert catalog.verify_manifest(man_before, root)
    assert cli.main(["release", "verify", "--root", str(root), "--release", "r26.11.1"]) == 1
    rdef.write_bytes(rdef_bytes)
    assert catalog.verify_manifest(man_before, root) == []
    capsys.readouterr()


def test_c54_verify_output_matches_stage_driver_and_release_train_hash_of(tmp_path, capsys):
    root = tmp_path / "root"
    _write_minimal_stage(root, "S0")
    for rel in ("r26.11.1", "r26.11.2", "r26.11.3"):
        _write_minimal_release(root, rel)

    assert cli.main(["release", "verify", "--root", str(root), "--release", "r26.11.1"]) == 0
    out_before = capsys.readouterr().out

    tar_path = tmp_path / "r26.11.1.tar"
    tar_sha = _make_backup_tar(root, "r26.11.1", tar_path)
    assert cli.main(["release", "retire", "--root", str(root), "--release", "r26.11.1", "--tar", f"{tar_path}={tar_sha}",
                     "--copy", _copy_arg(tar_sha)]) == 0
    capsys.readouterr()

    assert cli.main(["release", "verify", "--root", str(root), "--release", "r26.11.1"]) == 0
    out_after = capsys.readouterr().out
    assert "verified (archived: " in out_after
    assert "release files, stages S0, manifest hash " in out_after

    # Exact sed command used by stage_driver.sh and release_train.sh hash_of()
    log_before, log_after = tmp_path / "before.log", tmp_path / "after.log"
    log_before.write_text(out_before)
    log_after.write_text(out_after)
    sed_expr = r"s/.*manifest hash \([0-9a-f]*\).*/\1/p"
    h_before = subprocess.check_output(["sed", "-n", sed_expr, str(log_before)], text=True).strip().splitlines()[-1]
    h_after = subprocess.check_output(["sed", "-n", sed_expr, str(log_after)], text=True).strip().splitlines()[-1]
    assert len(h_before) == 32 and h_before == h_after


def test_c54_latest_and_releases_and_b49_skip_archived(tmp_path, capsys):
    root = tmp_path / "root"
    _write_minimal_stage(root, "S0")
    for rel in ("r26.11.1", "r26.11.2", "r26.11.3"):
        _write_minimal_release(root, rel)

    tar_path = tmp_path / "r26.11.1.tar"
    tar_sha = _make_backup_tar(root, "r26.11.1", tar_path)
    assert cli.main(["release", "retire", "--root", str(root), "--release", "r26.11.1", "--tar", f"{tar_path}={tar_sha}",
                     "--copy", _copy_arg(tar_sha)]) == 0
    capsys.readouterr()

    # latest_release skips archived
    assert catalog.latest_release(root) == "r26.11.3"
    # Even if r26.11.1 were the only one left, latest_release refuses with "archived"
    r_only = tmp_path / "only_arch"
    _write_minimal_stage(r_only, "S0")
    for rel in ("r26.11.1", "r26.11.2"):
        _write_minimal_release(r_only, rel)
    t1 = tmp_path / "o1.tar"
    h1 = _make_backup_tar(r_only, "r26.11.1", t1)
    assert cli.main(["release", "retire", "--root", str(r_only), "--release", "r26.11.1", "--keep-latest", "1", "--tar", f"{t1}={h1}",
                     "--copy", _copy_arg(h1, "o1.tar")]) == 0
    capsys.readouterr()
    (r_only / "manifests" / "release=r26.11.2.json").unlink()
    with pytest.raises(catalog.NoLatestRelease, match="archived"):
        catalog.latest_release(r_only)

    # bio-box releases TSV & JSON
    packs = tmp_path / "packs"
    packs.mkdir()
    assert cli.main(["releases", "--root", str(root), "--pack-dir", str(packs)]) == 0
    tsv = capsys.readouterr().out.splitlines()
    header = tsv[0].split("\t")
    assert tuple(header) == cli.RELEASES_COLUMNS + ("pack", "latest")
    rows = {ln.split("\t")[0]: dict(zip(header, ln.split("\t"))) for ln in tsv[1:]}
    assert rows["r26.11.1"]["manifest"] == "archived"
    assert rows["r26.11.1"]["genomes"] == "-"
    assert rows["r26.11.1"]["contigs"] == "-"
    assert rows["r26.11.1"]["pack"] == "-"
    assert rows["r26.11.1"]["latest"] == "-"
    assert rows["r26.11.3"]["manifest"] == "yes"
    assert rows["r26.11.3"]["latest"] == "read"

    assert cli.main(["releases", "--root", str(root), "--pack-dir", str(packs), "--json"]) == 0
    jrows = {r["release"]: r for r in (json.loads(ln) for ln in capsys.readouterr().out.splitlines())}
    assert jrows["r26.11.1"]["archived"] is True
    assert jrows["r26.11.1"]["genomes"] is None
    assert jrows["r26.11.1"]["contigs"] is None
    assert jrows["r26.11.1"]["pack"] is None

    # B-49 _source_problem skips archived release
    assert cluster._source_problem(root, "r26.11.4", "r26.11.1", {}) == "release is archived"
