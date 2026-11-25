"""Card W-41: pack._problems must sort shards numerically by part-<k>, not lexicographically.

When a target has >= 10 shards (part-0 .. part-16 for unique, part-0 .. part-24 for dark_orf),
lexicographic sorting places "part-9" last and "part-16" / "part-24" in the middle, falsely
flagging the real last shard as under the size floor and hiding an undersized part-9.
"""

from __future__ import annotations

import json
import random
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog, ids
from biobox.catalog import ReleaseFiles
from biobox.search import pack

pytestmark = pytest.mark.covers("W-41", "C-12", "W-35")

_AA = "ACDEFGHIKLMNPQRSTVWY"


def _write_zstd_fasta(path: Path, records: list[tuple[str, str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = "".join(f">{name}\n{seq}\n" for name, seq in records).encode("ascii")
    with pa.CompressedOutputStream(pa.OSFile(str(path), "wb"), "zstd") as fh:
        fh.write(payload)
    return path.stat().st_size


def _make_pack(
    tmp_path: Path,
    target: str,
    ext: str,
    n_shards: int,
    short_indices: set[int],
    *,
    release: str = "rW41.1",
) -> tuple[Path, ReleaseFiles, dict, int]:
    """Build a real, self-consistent minimal pack with `n_shards` shards on disk."""
    rng = random.Random(41)
    root = tmp_path / "root"
    root.mkdir(parents=True, exist_ok=True)
    packs = tmp_path / "packs"
    staging = tmp_path / "staging"
    staging.mkdir(parents=True, exist_ok=True)

    shards: list[dict] = []
    meta: list[dict] = []
    all_uids: list[bytes] = []
    full_sizes: list[int] = []
    short_sizes: list[int] = []

    seq_idx = 0
    for k in range(n_shards):
        rel_path = f"{target}/part-{k}{ext}"
        full_path = staging / rel_path
        uid = seq_idx.to_bytes(16, "big")
        all_uids.append(uid)
        seq_idx += 1
        if k in short_indices:
            seq = "M" + "".join(rng.choice(_AA) for _ in range(19))
        else:
            seq = "M" + "".join(rng.choice(_AA) for _ in range(1199))
        size = _write_zstd_fasta(full_path, [(uid.hex(), seq)])
        if k in short_indices:
            short_sizes.append(size)
        else:
            full_sizes.append(size)
        shards.append({
            "target": target,
            "path": rel_path,
            "bytes": size,
            "sha256": ids.sha256_file(full_path).hex(),
            "n_seq": 1,
        })

    assert full_sizes and short_sizes, "test setup requires both full and short shards"
    assert max(short_sizes) < min(full_sizes), (short_sizes, full_sizes)
    min_shard_bytes = (max(short_sizes) + min(full_sizes)) // 2

    gcv = "prodigal-test"
    manifest_files: list[dict] = []
    if target == "unique":
        pseq_rel = f"l1/protein_seq/gene_call_version={gcv}/stage=S0/part-0.parquet"
        pseq_path = root / pseq_rel
        pseq_path.parent.mkdir(parents=True, exist_ok=True)
        tbl = pa.table({"protein_uid": pa.array(all_uids, type=pa.binary(16))})
        pq.write_table(tbl, pseq_path)
        manifest_files.append({
            "table": "protein_seq",
            "path": pseq_rel,
            "partitions": {"gene_call_version": gcv, "stage": "S0"},
        })
    elif target == "dark_orf":
        meta_rel = "dark_orf/dark_orf_meta.parquet"
        meta_path = staging / meta_rel
        meta_path.parent.mkdir(parents=True, exist_ok=True)
        tbl = pa.table({"region_uid": pa.array(all_uids, type=pa.binary(16))})
        pq.write_table(tbl, meta_path)
        meta.append({
            "target": "dark_orf",
            "path": meta_rel,
            "bytes": meta_path.stat().st_size,
            "sha256": ids.sha256_file(meta_path).hex(),
        })

    rel_manifest_path = catalog.manifest_path(root, release)
    rel_manifest_path.parent.mkdir(parents=True, exist_ok=True)
    rel_manifest_path.write_text(
        json.dumps({"release_id": release, "stages": ["S0"], "files": manifest_files, "artifacts": []})
    )
    rf = ReleaseFiles.of(root, release)

    manifest = {
        "release_id": release,
        "gene_call_version": gcv,
        "shards": shards,
        "meta": meta,
        "targets": {target: {"n_seq": len(all_uids)}},
    }
    data = json.dumps(manifest, sort_keys=True, indent=1).encode()
    pid = pack.pack_id_of(data)
    out = pack.pack_dir(packs, release, pid)
    out.parent.mkdir(parents=True, exist_ok=True)
    staging.rename(out)
    (out / pack.MANIFEST).write_bytes(data)
    return out, rf, manifest, min_shard_bytes


def test_w41_unique_17_shards_only_last_under_floor_verifies(tmp_path: Path):
    """Positive: 17 shards (part-0..part-16) with only part-16 under min_shard_bytes -> 0 problems."""
    out, rf, manifest, min_bytes = _make_pack(tmp_path, "unique", ".faa.zst", 17, {16})
    problems = pack._problems(rf, "rW41.1", out, manifest, min_bytes, tmp_path / "scratch")
    assert problems == []


def test_w41_unique_17_shards_r26_11_3_exact_byte_floor(tmp_path: Path):
    """Exact r26.11.3 byte counts: 17 shards, part-16 = 124_053_559 < 268_435_456."""
    min_bytes = 268_435_456
    root = tmp_path / "root"
    root.mkdir()
    rel_manifest_path = catalog.manifest_path(root, "r26.11.3")
    rel_manifest_path.parent.mkdir(parents=True, exist_ok=True)
    rel_manifest_path.write_text(
        json.dumps({"release_id": "r26.11.3", "stages": ["S0"], "files": [], "artifacts": []})
    )
    rf = ReleaseFiles.of(root, "r26.11.3")
    out = tmp_path / "pack"
    out.mkdir()

    def size_floor_problems(short_idx: int) -> list[str]:
        shards = [
            {
                "target": "unique",
                "path": f"unique/part-{k}.faa.zst",
                "bytes": 124_053_559 if k == short_idx else min_bytes,
                "sha256": "0" * 64,
                "n_seq": 1,
            }
            for k in range(17)
        ]
        manifest = {
            "release_id": "r26.11.3",
            "gene_call_version": "prodigal-test",
            "shards": shards,
            "meta": [],
            "targets": {"unique": {"n_seq": 17}},
        }
        probs = pack._problems(rf, "r26.11.3", out, manifest, min_bytes, tmp_path / "scratch")
        return [p for p in probs if "and not the last shard" in p]

    assert size_floor_problems(16) == []
    assert size_floor_problems(9) == [
        f"unique/part-9.faa.zst: 124053559 bytes, under {min_bytes} and not the last shard"
    ]


def test_w41_unique_17_shards_part9_under_floor_flagged(tmp_path: Path):
    """Negative 1: 17 shards where part-9 (middle shard, not last) is under min_shard_bytes."""
    out, rf, manifest, min_bytes = _make_pack(tmp_path, "unique", ".faa.zst", 17, {9})
    problems = pack._problems(rf, "rW41.1", out, manifest, min_bytes, tmp_path / "scratch")
    part9_bytes = manifest["shards"][9]["bytes"]
    assert problems == [
        f"unique/part-9.faa.zst: {part9_bytes} bytes, under {min_bytes} and not the last shard"
    ]


def test_w41_unique_17_shards_part2_and_part16_under_floor_flags_only_part2(tmp_path: Path):
    """Negative 2: part-2 and part-16 both under min_shard_bytes -> flags part-2, not part-16."""
    out, rf, manifest, min_bytes = _make_pack(tmp_path, "unique", ".faa.zst", 17, {2, 16})
    problems = pack._problems(rf, "rW41.1", out, manifest, min_bytes, tmp_path / "scratch")
    part2_bytes = manifest["shards"][2]["bytes"]
    assert problems == [
        f"unique/part-2.faa.zst: {part2_bytes} bytes, under {min_bytes} and not the last shard"
    ]


def test_w41_dark_orf_25_shards_order(tmp_path: Path):
    """25 dark_orf shards (part-0..part-24): only part-24 may be under the size floor."""
    out_ok, rf_ok, man_ok, min_ok = _make_pack(tmp_path / "ok", "dark_orf", ".faa.zst", 25, {24})
    assert pack._problems(rf_ok, "rW41.1", out_ok, man_ok, min_ok, tmp_path / "s1") == []

    out_bad, rf_bad, man_bad, min_bad = _make_pack(tmp_path / "bad", "dark_orf", ".faa.zst", 25, {9, 24})
    part9_bytes = man_bad["shards"][9]["bytes"]
    assert pack._problems(rf_bad, "rW41.1", out_bad, man_bad, min_bad, tmp_path / "s2") == [
        f"dark_orf/part-9.faa.zst: {part9_bytes} bytes, under {min_bad} and not the last shard"
    ]
