"""Ring R2: bounded, row-group-pruning `export.protein_seqs` and `hits_faa` at multi-stage scale (card B-65).

Oracles are independent brute-force loops over the Parquet files written here.
"""

from __future__ import annotations

import hashlib
import shutil
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

sys.path.insert(0, str(Path(__file__).parent))
import test_release_readers as rr  # noqa: E402
from test_release_readers import R2, V, store  # noqa: E402,F401

from biobox import catalog  # noqa: E402
from biobox.l3 import tapes  # noqa: E402
from biobox.registry.schemas import table  # noqa: E402
from biobox.search import export  # noqa: E402

pytestmark = [pytest.mark.r2, pytest.mark.covers("B-65", "W-23", "B-27")]

REPO = Path(__file__).resolve().parents[2]
GCV = V
REL = "rS0.b65"
STAGES = ("S0", "S1", "S2", "S3")


def _write_protein_seq(root: Path, stage: str, rows: list[tuple[bytes, int, str]], *,
                       row_group_size: int = 4, write_statistics: bool = True,
                       prefix: str = "00") -> Path:
    schema = table("protein_seq").schema
    path = root / f"l0seq/protein_seq/gene_call_version={GCV}/stage={stage}/protein_prefix={prefix}/part-00000.parquet"
    path.parent.mkdir(parents=True, exist_ok=True)
    tbl = pa.table({
        "protein_uid": pa.array([u for u, _, _ in rows], type=pa.binary(16)),
        "protein_digest": pa.array([hashlib.sha256(s.encode()).digest() for _, _, s in rows], type=pa.binary(32)),
        "length": pa.array([n for _, n, _ in rows], type=pa.uint32()),
        "seq": pa.array([s for _, _, s in rows], type=pa.string()),
    }, schema=schema)
    pq.write_table(tbl, path, row_group_size=row_group_size, write_statistics=write_statistics)
    return path


def _oracle_protein_seqs(root: Path, uids: list[bytes]) -> dict[bytes, tuple[int, str]]:
    want = set(uids)
    out: dict[bytes, tuple[int, str]] = {}
    for p in sorted((root / "l0seq/protein_seq").rglob("*.parquet")):
        for r in pq.read_table(p, columns=["protein_uid", "length", "seq"], partitioning=None).to_pylist():
            u = r["protein_uid"]
            if u not in want:
                continue
            n, s = int(r["length"]), str(r["seq"])
            if u not in out:
                out[u] = (n, s)
            else:
                out[u] = (max(out[u][0], n), min(out[u][1], s))
    return out


def _build_clustered_store(root: Path) -> tuple[list[bytes], int]:
    """4 stages x 15 row groups (row_group_size=4) = 60 row groups total.

    Each stage holds:
      - 16 low-key rows in row groups 0..3 (leading byte 0x00)
      - 40 middle-key rows in row groups 4..13 (leading bytes 0x10..0x9f)
      - 4 high-key rows in row group 14 (leading byte 0xff)
    Selecting all 64 low keys plus 1 high key gives 65 unique keys (> 60 total row groups, exceeding DuckDB's
    dynamic_or_filter_threshold = max(50, 60)) whose [min, max] spans every row group of every stage, while
    only 17 of 60 row groups actually intersect the keys.
    """
    low_uids: list[bytes] = []
    high_uid: bytes | None = None
    for si, stage in enumerate(STAGES):
        rows: list[tuple[bytes, int, str]] = []
        for i in range(16):
            u = bytes([0x00, si, i >> 8, i & 0xFF]) + b"\x00" * 12
            s = f"M{stage}L{i:03d}" + ("A" * (20 + i))
            rows.append((u, len(s), s))
            low_uids.append(u)
        for i in range(40):
            u = bytes([0x10 + i, si, 0, i]) + b"\x00" * 12
            s = f"M{stage}M{i:03d}" + ("C" * (25 + i))
            rows.append((u, len(s), s))
        for i in range(4):
            u = bytes([0xFF, si, 0, i]) + b"\x00" * 12
            s = f"M{stage}H{i:03d}" + ("G" * (30 + i))
            rows.append((u, len(s), s))
            if si == 0 and i == 0:
                high_uid = u
        rows.sort(key=lambda r: r[0])
        _write_protein_seq(root, stage, rows, row_group_size=4)
    catalog.define_release(root, REL, list(STAGES))
    assert high_uid is not None
    return low_uids + [high_uid], 60


def test_b65_protein_seqs_prunes_above_duckdb_or_filter_threshold(tmp_path):
    """When len(uids) > total row groups (exceeding DuckDB's dynamic_or_filter_threshold) and [min, max] spans
    every row group, protein_seqs still reads only the row groups that intersect the keys."""
    root = tmp_path / "clustered"
    uids, total_rgs = _build_clustered_store(root)
    assert len(uids) == 65 > total_rgs

    got, scans = export.protein_seqs(root, uids, release=REL, gene_call_version=GCV)
    assert got == _oracle_protein_seqs(root, uids)
    assert len(got) == 65

    rgs_scanned = sum(s["row_groups_scanned"] for s in scans)
    rgs_total = sum(s["row_groups_total"] for s in scans)
    assert rgs_total == total_rgs
    assert rgs_scanned == 17, f"expected exact 17 intersecting row groups out of {total_rgs}, got {rgs_scanned}: {scans}"


def test_b65_protein_seqs_cross_stage_empty_missing_and_undefined_release(tmp_path):
    """Cross-stage aggregation matches max(length), min(seq); empty/missing/duplicate inputs and undefined
    release behave identically to the contract."""
    root = tmp_path / "cross"
    shared_u = b"\x05" * 16
    only_s0 = b"\x01" * 16
    only_s1 = b"\x09" * 16
    missing = b"\x77" * 16
    _write_protein_seq(root, "S0", [
        (only_s0, 5, "MAAAA"),
        (shared_u, 6, "MBBBBB"),
    ], row_group_size=2)
    _write_protein_seq(root, "S1", [
        (shared_u, 8, "MAAAAA"),
        (only_s1, 7, "MCCCCCC"),
    ], row_group_size=2)
    catalog.define_release(root, REL, ["S0", "S1"])

    got, scans = export.protein_seqs(root, [shared_u, only_s0, missing, shared_u], release=REL, gene_call_version=GCV)
    assert got == _oracle_protein_seqs(root, [shared_u, only_s0, missing])
    assert got[shared_u] == (8, "MAAAAA")
    assert missing not in got
    assert sum(s["row_groups_scanned"] for s in scans) == 2

    empty_got, empty_scans = export.protein_seqs(root, [], release=REL, gene_call_version=GCV)
    assert empty_got == {}
    assert sum(s["row_groups_scanned"] for s in empty_scans) == 0
    assert sum(s["files_read"] for s in empty_scans) == 0

    with pytest.raises(ValueError, match="rNoSuch"):
        export.protein_seqs(root, [only_s0], release="rNoSuch", gene_call_version=GCV)


def test_b65_protein_seqs_refuses_unsorted_or_statless_parquet(tmp_path):
    """LY-02 / LY-03: protein_seqs refuses protein_seq Parquet files lacking min/max stats or not sorted."""
    root_nostat = tmp_path / "nostat"
    u0, u1 = b"\x01" * 16, b"\x02" * 16
    _write_protein_seq(root_nostat, "S0", [(u0, 4, "MAAA"), (u1, 4, "MCCC")],
                       row_group_size=1, write_statistics=False)
    catalog.define_release(root_nostat, REL, ["S0"])
    with pytest.raises(ValueError, match="LY-03"):
        export.protein_seqs(root_nostat, [u0], release=REL, gene_call_version=GCV)

    root_unsorted = tmp_path / "unsorted"
    _write_protein_seq(root_unsorted, "S0", [(u1, 4, "MCCC"), (u0, 4, "MAAA")],
                       row_group_size=1, write_statistics=True)
    catalog.define_release(root_unsorted, REL, ["S0"])
    with pytest.raises(ValueError, match="LY-02"):
        export.protein_seqs(root_unsorted, [u0], release=REL, gene_call_version=GCV)


def test_b65_open_store_disables_preserve_insertion_order(tmp_path):
    """open_store sets preserve_insertion_order=false so DuckDB does not buffer unordered Parquet chunks."""
    root = tmp_path / "opt"
    _write_protein_seq(root, "S0", [(b"\x01" * 16, 4, "MAAA")], row_group_size=1)
    catalog.define_release(root, REL, ["S0"])
    with export.open_store(root, release=REL, gene_call_version=GCV, relations=("prot",)) as st:
        val = st.con.execute("SELECT current_setting('preserve_insertion_order')").fetchone()[0]
    assert val is False or str(val).lower() == "false"


def test_b65_hits_faa_streams_in_batches_matching_oracle(store, tmp_path, monkeypatch):
    """hits_faa streams chunks through protein_seqs across multiple stages (including shared proteins) and
    matches an independent Parquet oracle byte-for-byte, raising ValueError on a missing protein_seq row."""
    root = tmp_path / "r2_copy"
    shutil.copytree(store["root"], root)
    for f in sorted((root / "l0seq/protein_seq").rglob("*.parquet")):
        pq.write_table(pq.read_table(f, partitioning=None), f, row_group_size=4)

    all_prots = _oracle_protein_seqs(root, list(store["s0_proteins"] | store["s1_only"]))
    sample = sorted(all_prots)[:15]
    if store["p"] not in sample:
        sample.append(store["p"])  # shared across S0 and S1
    if store["p_mut"] not in sample:
        sample.append(store["p_mut"])  # S1 only
    sample = sorted(set(sample))
    hits = [
        {"target_kind": "protein", "target_uid": u, "seq_bitscore": 100.0 + i, "seq_evalue_zref": 1e-10 / (i + 1)}
        for i, u in enumerate(sample)
    ]
    tv = tapes.tape_version([V])

    # Baseline with default batch size vs forced small batch size (_BATCH_ROWS = 4)
    monkeypatch.setattr(export, "_BATCH_ROWS", 4)
    lines = list(export.hits_faa(root, hits, release=R2, tape_version=tv, gene_call_version=V))
    assert len(lines) == 2 * len(sample)
    got_seqs = {}
    for hdr, seq in zip(lines[0::2], lines[1::2]):
        assert hdr.startswith(">")
        u = bytes.fromhex(hdr[1:].split(" ", 1)[0])
        got_seqs[u] = seq
        assert f"length={all_prots[u][0]}" in hdr
    assert list(got_seqs) == sample
    assert got_seqs == {u: all_prots[u][1] for u in sample}

    # Missing protein_seq row raises the exact ValueError
    ghost = b"\xde\xad" * 8
    bad_hits = [{"target_kind": "protein", "target_uid": ghost, "seq_bitscore": 50.0, "seq_evalue_zref": 1e-5}]
    with pytest.raises(ValueError, match=f"hit protein {ghost.hex()} has no protein_seq row"):
        list(export.hits_faa(root, bad_hits, release=R2, tape_version=tv, gene_call_version=V))


def test_b65_example_scripts_pass_threads_to_seq_oracle():
    """Both worked example scripts pass --threads \"$THREADS\" to `bb seq-oracle seq`."""
    for rel_path in ("examples/crispr_cas_neighbourhoods/run.sh", "examples/rnap/run.sh"):
        text = (REPO / rel_path).read_text(encoding="utf-8")
        seq_lines = [ln for ln in text.splitlines() if "bb seq-oracle seq " in ln]
        assert len(seq_lines) == 1, f"{rel_path}: expected 1 seq-oracle line, found {seq_lines}"
        assert '--threads "$THREADS"' in seq_lines[0], f"{rel_path}: missing --threads \"$THREADS\" in {seq_lines[0]}"
