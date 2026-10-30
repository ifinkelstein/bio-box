"""Ring R2: micro-S0 fixture through O-INGEST and O-SLICE.

Covers R-01, R-02, R-03 (data_layer_and_tests.md §2A.5), C-01/C-02/C-06 via the
registry contract checker, LY-07 and LY-08 (TEST_ARCHITECTURE §7), and
MR-01, MR-02, MR-03, MR-14 (TEST_ARCHITECTURE §4).
"""

from __future__ import annotations

import gzip
import hashlib
import random
import shutil
from pathlib import Path

import duckdb
import pyarrow.parquet as pq
import pytest
import yaml

from biobox import fasta, ids
from biobox.l0seq import SeqStore, SourceObject, ingest, revcomp
from biobox.registry.schemas import conformance_problems, table

pytestmark = pytest.mark.r2

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())


def _objects(paths):
    return [
        SourceObject(Path(p), CONTRACT["source_name"], CONTRACT["source_release"], CONTRACT["license"], True)
        for p in paths
    ]


def _fixture_paths():
    return [FIX / o["file"] for o in CONTRACT["objects"]]


def _files(root: Path) -> dict[str, bytes]:
    return {str(p.relative_to(root)): p.read_bytes() for p in sorted(root.rglob("*.parquet"))}


@pytest.fixture(scope="module")
def store_root(tmp_path_factory):
    root = tmp_path_factory.mktemp("micro_s0")
    # small chunks and row groups so slices cross chunk and row-group boundaries
    ingest(_objects(_fixture_paths()), root, stage="S0", chunk_size=4096, row_group_chunks=16)
    return root


# --- contract -----------------------------------------------------------------


def test_fixture_matches_contract():
    for o in CONTRACT["objects"]:
        p = FIX / o["file"]
        assert ids.sha256_file(p).hex() == o["sha256_gz"]
        with gzip.open(p, "rb") as fh:
            assert hashlib.sha256(fh.read()).hexdigest() == o["sha256_fna"]
        recs = fasta.read(p)
        assert [r.accession for r in recs] == [o["accession"]]
        assert len(recs[0].seq) == o["length"]


@pytest.mark.covers("C-01", "C-02", "C-06", "LY-01")
@pytest.mark.parametrize(
    "name,glob",
    [
        ("seq_chunks", "l0seq/seq_chunks/stage=*/digest_prefix=*/part-0.parquet"),
        ("contig_occurrence", "l0seq/contig_occurrence/stage=*/source=*/part-0.parquet"),
        ("contig_index", "l0seq/contig_index/stage=*/part-0.parquet"),
        ("seq_locator", "l0seq/seq_locator.parquet"),
        ("l0_manifest", "l0/l0_manifest.parquet"),
    ],
)
def test_written_tables_conform(store_root, name, glob):
    files = sorted(store_root.glob(glob))
    assert files, f"no files for {name}"
    for f in files:
        assert conformance_problems(pq.read_table(f), table(name)) == [], f


# --- R-01, R-02, R-03 -----------------------------------------------------------


def test_r01_record_digests(store_root):
    store = SeqStore(store_root)
    for p in _fixture_paths():
        for rec in fasta.read(p):
            rd = ids.raw_digest(rec.seq)
            assert ids.raw_digest(store.raw(rd)) == rd
            assert store.get(ids.contig_uid(rec.seq)) == rec.seq


def test_r02_object_round_trip(store_root):
    store = SeqStore(store_root)
    for p in _fixture_paths():
        expected = hashlib.sha256(b"\n".join(fasta.normalized_lines(p)) + b"\n").digest()
        got = store.reconstruct_object(ids.sha256_file(p))
        assert hashlib.sha256(b"\n".join(got) + b"\n").digest() == expected, p.name


def test_r03_random_slices_both_strands(store_root):
    store = SeqStore(store_root)
    rng = random.Random(20261029)
    for p in _fixture_paths():
        seq = fasta.read(p)[0].seq
        uid = ids.contig_uid(seq)
        for _ in range(200):
            a, b = sorted(rng.randrange(0, len(seq) + 1) for _ in range(2))
            assert store.get(uid, a, b, 1) == seq[a:b]
            assert store.get(uid, a, b, -1) == revcomp(seq[a:b])
        assert store.get(uid, len(seq) - 1, len(seq)) == seq[-1:]


def test_slice_out_of_bounds_raises(store_root):
    store = SeqStore(store_root)
    seq = fasta.read(_fixture_paths()[0])[0].seq
    with pytest.raises(ValueError):
        store.get(ids.contig_uid(seq), 0, len(seq) + 1)


# --- layout tests ----------------------------------------------------------------


def test_ly07_duckdb_reads_the_same_rows(store_root):
    files = sorted(str(f) for f in store_root.glob("l0seq/seq_chunks/stage=*/digest_prefix=*/part-0.parquet"))
    con = duckdb.connect()
    n_duck, bytes_duck = con.execute(
        "SELECT count(*), sum(length(seq)) FROM read_parquet(?)", [files]
    ).fetchone()
    arrow = [pq.read_table(f) for f in files]
    assert n_duck == sum(t.num_rows for t in arrow)
    assert bytes_duck == sum(sum(len(s) for s in t.column("seq").to_pylist()) for t in arrow)
    typ = con.execute("SELECT typeof(raw_digest) FROM read_parquet(?) LIMIT 1", [files]).fetchone()[0]
    assert typ == "BLOB"


def test_ly08_rewrite_is_byte_identical(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    for r in (a, b):
        ingest(_objects(_fixture_paths()), r, stage="S0", chunk_size=4096, row_group_chunks=16)
    assert _files(a) == _files(b)


# --- metamorphic relations -------------------------------------------------------


def test_mr01_input_order_does_not_matter(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    ingest(_objects(_fixture_paths()), a, stage="S0", chunk_size=4096, row_group_chunks=16)
    ingest(_objects(list(reversed(_fixture_paths()))), b, stage="S0", chunk_size=4096, row_group_chunks=16)
    assert _files(a) == _files(b)


def _write_fasta(path: Path, header: str, seq: bytes, width: int = 70) -> Path:
    path.write_bytes(b"\n".join(fasta.format_record(header, seq, width)) + b"\n")
    return path


def test_mr02_renamed_header_keeps_uids_and_chunks(tmp_path):
    lam = fasta.read(_fixture_paths()[0])[0]
    src = tmp_path / "src"
    src.mkdir()
    renamed = _write_fasta(src / "renamed.fna", "XX_000001.1 renamed lambda", lam.seq)
    a, b = tmp_path / "a", tmp_path / "b"
    ingest(_objects([_fixture_paths()[0]]), a, stage="S0", chunk_size=4096, row_group_chunks=16)
    ingest(_objects([renamed]), b, stage="S0", chunk_size=4096, row_group_chunks=16)
    fa, fb = _files(a), _files(b)
    for k in fa:
        if "seq_chunks" in k or "contig_index" in k or "seq_locator" in k:
            assert fa[k] == fb[k], k
    occ_a = pq.read_table(next(a.glob("l0seq/contig_occurrence/**/*.parquet"))).to_pylist()[0]
    occ_b = pq.read_table(next(b.glob("l0seq/contig_occurrence/**/*.parquet"))).to_pylist()[0]
    assert occ_a["contig_uid"] == occ_b["contig_uid"]
    assert occ_a["accession"] != occ_b["accession"]


def test_mr03_reverse_complement_slices_mirror(tmp_path):
    lam = fasta.read(_fixture_paths()[0])[0]
    src = tmp_path / "src"
    src.mkdir()
    rc_path = _write_fasta(src / "lambda_rc.fna", "lambda_rc", revcomp(lam.seq))
    root = tmp_path / "s"
    ingest(_objects([_fixture_paths()[0], rc_path]), root, stage="S0", chunk_size=4096, row_group_chunks=16)
    store = SeqStore(root)
    fwd, rc = ids.contig_uid(lam.seq), ids.contig_uid(revcomp(lam.seq))
    n = len(lam.seq)
    rng = random.Random(3)
    for _ in range(200):
        a, b = sorted(rng.randrange(0, n + 1) for _ in range(2))
        assert store.get(rc, n - b, n - a, -1) == store.get(fwd, a, b, 1)


def test_mr14_soft_mask_keeps_contig_uid(tmp_path):
    lam = fasta.read(_fixture_paths()[0])[0]
    masked = lam.seq[:1000].lower() + lam.seq[1000:]
    src = tmp_path / "src"
    src.mkdir()
    mpath = _write_fasta(src / "lambda_masked.fna", "lambda_masked", masked)
    root = tmp_path / "s"
    objs = _objects([_fixture_paths()[0]]) + [
        SourceObject(mpath, "other_source", "1", "INSDC (public)", True, precedence=200)
    ]
    ingest(objs, root, stage="S0", chunk_size=4096, row_group_chunks=16)
    store = SeqStore(root)
    assert ids.contig_uid(masked) == ids.contig_uid(lam.seq)
    assert ids.raw_digest(masked) != ids.raw_digest(lam.seq)
    assert store.raw(ids.raw_digest(masked)) == masked
    assert store.get(ids.contig_uid(lam.seq)) == lam.seq  # default variant: higher-precedence source


def test_identical_sequence_stored_once_across_sources_and_stages(tmp_path):
    root = tmp_path / "s"
    ingest(_objects([_fixture_paths()[0]]), root, stage="S0", chunk_size=4096, row_group_chunks=16)
    before = _files(root)
    copy_dir = tmp_path / "copy"
    copy_dir.mkdir()
    dup = copy_dir / "lambda_again.fna.gz"
    shutil.copy(_fixture_paths()[0], dup)
    report = ingest(
        [SourceObject(dup, "second_source", "1", "INSDC (public)", True)], root, stage="S1", chunk_size=4096, row_group_chunks=16
    )
    assert report["new_sequences"] == 0 and report["index_rows"] == 0
    after = _files(root)
    for k, v in before.items():
        if "seq_chunks" in k or "contig_index" in k or "contig_occurrence" in k:
            assert after[k] == v, f"S1 rewrote an S0 file: {k}"
    assert not list(root.glob("l0seq/seq_chunks/stage=S1/**/*.parquet"))
    assert list(root.glob("l0seq/contig_occurrence/stage=S1/source=second_source/part-0.parquet"))
