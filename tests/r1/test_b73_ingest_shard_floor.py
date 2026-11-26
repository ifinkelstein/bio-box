"""Ring R1: coalesce undersized seq_shards batches against LY-04 (card B-73; LY-04, D-01, D-03).

Oracles:
- `contracts.check_file_sizes` (LY-04) on the resulting `seq_shards` manifest entries.
- Single-process `store.ingest` vs batched `ingest_pass1`/`ingest_plan`/`ingest_batch`/`ingest_merge`
  byte identity (D-01, D-03).
"""

from __future__ import annotations

import hashlib
import random
from pathlib import Path

import pytest

from biobox import contracts
from biobox.l0seq import SourceObject, store

pytestmark = [pytest.mark.covers("LY-04", "D-01", "D-03")]

PLAN_RECORDS = "l0seq/ingest_plan/"


def _is_contiguous_partition(batches: list[list[int]], per_object: list[tuple[int, int]]) -> bool:
    flat = [i for b in batches for i in b]
    if flat != [i for i, _ in per_object]:
        return False
    return all(len(b) > 0 for b in batches)


def _batch_totals(batches: list[list[int]], per_object: list[tuple[int, int]]) -> list[int]:
    bases_of = dict(per_object)
    return [sum(bases_of[i] for i in b) for b in batches]


def _write_fasta(path: Path, n_bases: int, seed: int, name: str) -> Path:
    rng = random.Random(seed)
    seq = "".join(rng.choices("ACGT", k=n_bases))
    lines = [f">{name}\n"] + [seq[i : i + 80] + "\n" for i in range(0, len(seq), 80)]
    path.write_text("".join(lines))
    return path


def _tree(root: Path) -> dict[str, str]:
    return {
        p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def _shard_manifest(root: Path, stage: str) -> dict:
    shard_dir = root / f"l0seq/seq_shards/stage={stage}"
    artifacts = [
        {
            "artifact": "seq_shards",
            "path": p.relative_to(root).as_posix(),
            "partitions": {"stage": stage},
            "bytes": p.stat().st_size,
        }
        for p in sorted(shard_dir.glob("shard-*.fa.gz"))
    ]
    return {"root": str(root), "files": [], "artifacts": artifacts}


def test_leading_small_plus_large_v2_shape_coalesces():
    # V2 incident shape: jgi_inovirus (~50 Mbp) + jgi_img_vr (~8.15 Gbp)
    per_object = [(0, 50_000_000), (1, 8_150_000_000)]
    got = store._batches(per_object, store.SHARD_BASES, store.MIN_SHARD_BASES)
    assert got == [[0, 1]]


def test_large_plus_trailing_small_o1_shape_coalesces():
    per_object = [(0, 2_500_000_000), (1, 45_000_000)]
    got = store._batches(per_object, store.SHARD_BASES, store.MIN_SHARD_BASES)
    assert got == [[0, 1]]


def test_interior_small_batch_merges_into_smaller_adjacent_neighbor():
    # Three greedy batches: [0] (2.0G), [1] (50M, undersized), [2] (1.2G) -> merges right into [2]
    per_object = [(0, 2_000_000_000), (1, 50_000_000), (2, 1_200_000_000)]
    got = store._batches(per_object, store.SHARD_BASES, store.MIN_SHARD_BASES)
    assert got == [[0], [1, 2]]
    assert _is_contiguous_partition(got, per_object)

    # Mirror: left neighbor is smaller -> merges left into [0]
    per_object_left = [(0, 1_200_000_000), (1, 50_000_000), (2, 2_000_000_000)]
    got_left = store._batches(per_object_left, store.SHARD_BASES, store.MIN_SHARD_BASES)
    assert got_left == [[0, 1], [2]]
    assert _is_contiguous_partition(got_left, per_object_left)

    # Tie between left and right neighbors -> deterministic tie-break to left neighbor
    per_object_tie = [(0, 1_500_000_000), (1, 50_000_000), (2, 1_500_000_000)]
    got_tie = store._batches(per_object_tie, store.SHARD_BASES, store.MIN_SHARD_BASES)
    assert got_tie == [[0, 1], [2]]
    assert _is_contiguous_partition(got_tie, per_object_tie)


def test_multi_object_pl2_shape_meets_floor_and_preserves_order():
    per_object = [
        (0, 30_000_000),
        (1, 20_000_000),
        (2, 1_400_000_000),
        (3, 80_000_000),
        (4, 1_800_000_000),
        (5, 40_000_000),
        (6, 60_000_000),
        (7, 1_100_000_000),
        (8, 25_000_000),
    ]
    got = store._batches(per_object, store.SHARD_BASES, store.MIN_SHARD_BASES)
    assert _is_contiguous_partition(got, per_object)
    assert len(got) > 1
    totals = _batch_totals(got, per_object)
    assert all(t >= store.MIN_SHARD_BASES for t in totals)


def test_guard_leaves_greedy_cut_unchanged_when_floor_disabled_or_above_shard_bases():
    per_object = [(0, 50_000), (1, 180_000)]
    # Default min_shard_bases=0: greedy cut unchanged (two batches)
    assert store._batches(per_object, 200_000) == [[0], [1]]
    assert store._batches(per_object, 200_000, 0) == [[0], [1]]
    # min_shard_bases >= shard_bases (e.g. small-fixture SHARD_BASES=200_000 with default MIN_SHARD_BASES):
    # guard disables coalescing so small-fixture tests are untouched
    assert store._batches(per_object, 200_000, store.MIN_SHARD_BASES) == [[0], [1]]
    assert store._batches(per_object, 200_000, 200_000) == [[0], [1]]


def test_ingest_and_batched_ingest_byte_identical_with_coalescing_and_pass_ly04(tmp_path, monkeypatch):
    src_dir = tmp_path / "src"
    src_dir.mkdir()
    # One small object (~12 kb) and two large objects (~160 kb each); at shard_bases=120_000,
    # greedy cutting makes 3 shards (one tiny, two large), while min_shard_bases=40_000
    # coalesces the small object into its adjacent large shard -> 2 shards, both above floor.
    paths = [
        _write_fasta(src_dir / "small.fa", 12_000, seed=101, name="small_contig"),
        _write_fasta(src_dir / "large_a.fa", 160_000, seed=202, name="large_a_contig"),
        _write_fasta(src_dir / "large_b.fa", 160_000, seed=303, name="large_b_contig"),
    ]
    objects = [SourceObject(p, "synth", "1", "public", True) for p in paths]
    shard_bases = 120_000
    min_bytes = 20_000  # between the small bgzipped shard (~4 KB) and the large shards (~50 KB)

    # First verify the uncoalesced baseline trips LY-04
    monkeypatch.setattr(store, "MIN_SHARD_BASES", 0)
    raw_root = tmp_path / "uncoalesced"
    store.ingest(objects, raw_root, "S0", shard_bases=shard_bases)
    raw_manifest = _shard_manifest(raw_root, "S0")
    assert len(raw_manifest["artifacts"]) == 3
    assert contracts.check_file_sizes(raw_manifest, min_bytes=min_bytes)

    # Now enable coalescing and compare single-process ingest vs batched ingest
    monkeypatch.setattr(store, "MIN_SHARD_BASES", 40_000)
    single_root = tmp_path / "single"
    store.ingest(objects, single_root, "S0", shard_bases=shard_bases)

    fixed_manifest = _shard_manifest(single_root, "S0")
    assert len(fixed_manifest["artifacts"]) == 2
    assert contracts.check_file_sizes(fixed_manifest, min_bytes=min_bytes) == []

    batched_root = tmp_path / "batched"
    scratch = tmp_path / "scratch"
    pass1_bytes = 1  # force one pass-1 shard per object
    n_pass1 = len(store.pass1_shards(objects, pass1_bytes))
    assert n_pass1 == 3
    for k in range(n_pass1):
        store.ingest_pass1(objects, batched_root, "S0", k, pass1_bytes, shard_bases=shard_bases, tmp=scratch)
    store.ingest_plan(objects, batched_root, "S0", shard_bases=shard_bases, tmp=scratch)
    plan = store.read_ingest_plan(batched_root, "S0")
    assert len(plan["batches"]) == 2
    for j in range(len(plan["batches"])):
        store.ingest_batch(objects, batched_root, "S0", j, tmp=scratch)
    store.ingest_merge(objects, batched_root, "S0", tmp=scratch)

    batched_data = {k: v for k, v in _tree(batched_root).items() if not k.startswith(PLAN_RECORDS)}
    assert batched_data == _tree(single_root)
