"""Ring R1, card B-55: `crispr._Scan` buffers scratch Parquet writes to `io.ROW_GROUP_ROWS`.

Before B-55, `crispr._Scan` called `pq.ParquetWriter.write_batch` with a 1-row `RecordBatch` on every `contig()` and
`skip()` call, and once per MinCED batch in `_collect()`, without passing `row_group_size`. On multi-million-contig
stages (3.93M contigs in the C-42 pilot), `contigs.parquet` held 3.93M 1-row row groups (~2.19 GB of Parquet footer
metadata). These tests verify that every scratch Parquet file (`contigs.parquet`, `arrays.parquet`, `elements.parquet`,
`topology.parquet`) has at most `max(1, ceil(num_rows / io.ROW_GROUP_ROWS))` row groups and preserves every row.
"""

from __future__ import annotations

import math
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import ids, io
from biobox.l1 import crispr

pytestmark = pytest.mark.covers("B-55", "E-12", "E-15")


def _stub_search_factory(every_n: int = 50):
    """Stub for `crispr._search` that emits one synthetic array (3 repeats, 2 spacers) on every `every_n`-th record
    of each batch, unlinks the batch FASTA, and returns the three scratch RecordBatches."""

    def _stub_search(java: str, jar: Path, batch: crispr._Batch, scratch: Path,
                     min_repeats: int = crispr.MIN_REPEATS):
        try:
            feats: list[dict] = []
            els: list[dict] = []
            n_in: dict[bytes, int] = {}
            n_rej: dict[bytes, int] = {}
            for idx, rec in enumerate(batch.recs.values()):
                if idx % every_n != 0:
                    continue
                cu = rec.contig
                n_in[cu] = n_in.get(cu, 0) + 1
                start, end = 10, 170
                fu = ids.feature_uid(cu, start, end, crispr.STRAND, crispr.KIND_CRISPR)
                feats.append({"feature_uid": fu, "contig_uid": cu, "start": start, "end": end,
                              "strand": crispr.STRAND})
                units = ((10, 38), (71, 99), (142, 170))
                for kind, ordinal, a, b in crispr.elements(units):
                    els.append({"feature_uid": fu, "contig_uid": cu,
                                "region_uid": ids.region_uid(cu, a, b, crispr.STRAND),
                                "element": kind, "ordinal": ordinal, "start": a, "end": b,
                                "strand": crispr.STRAND, "seq": "A" * (b - a), "contig_len": rec.n})
            keys = sorted(n_in)
            contigs = pa.RecordBatch.from_arrays(
                [pa.array(keys, type=pa.binary(16)), pa.array([n_in[k] for k in keys], type=pa.int64()),
                 pa.array([n_rej.get(k, 0) for k in keys], type=pa.int64()),
                 pa.array([None] * len(keys), pa.string())],
                schema=crispr._CONTIG_SCRATCH)
            return (pa.RecordBatch.from_pylist(feats, schema=crispr._ARRAY_SCRATCH),
                    pa.RecordBatch.from_pylist(els, schema=crispr._ELEMENT_SCRATCH),
                    contigs)
        finally:
            batch.path.unlink(missing_ok=True)

    return _stub_search


def _assert_bounded_row_groups(paths: dict[str, Path]) -> dict[str, tuple[int, int]]:
    out: dict[str, tuple[int, int]] = {}
    for name, path in paths.items():
        pf = pq.ParquetFile(path)
        n_rows = pf.metadata.num_rows
        n_rg = pf.num_row_groups
        bound = max(1, math.ceil(n_rows / io.ROW_GROUP_ROWS))
        assert n_rg <= bound, (f"scratch {name}.parquet has {n_rg} row groups for {n_rows} rows "
                               f"(bound {bound} at ROW_GROUP_ROWS={io.ROW_GROUP_ROWS})")
        out[name] = (n_rows, n_rg)
    return out


def test_scan_buffers_10k_contigs_to_row_group_rows(tmp_path, monkeypatch):
    monkeypatch.setattr(crispr, "_search", _stub_search_factory(every_n=50))
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    sc = crispr._Scan("java", tmp_path / "minced.jar", cpus=2, scratch=scratch, batch_bases=100_000, topology=True)

    n_skip = 500
    n_contigs = 10_000
    seq = b"ACGT" * 50  # 200 bp
    for i in range(n_skip):
        sc.skip(i.to_bytes(16, "big"), "S0")
    for i in range(n_contigs):
        cu = (n_skip + i).to_bytes(16, "big")
        sc.contig(cu, b"" if i == 0 else seq, circular=(i % 7 == 0))
    sc.close(ok=True)

    summary = _assert_bounded_row_groups(sc.paths)
    assert sc.stats["contigs"] == n_contigs
    assert sc.stats["batches"] > 1
    assert sc.stats["arrays"] > 0
    assert sc.stats["repeats"] == 3 * sc.stats["arrays"]
    assert sc.stats["spacers"] == 2 * sc.stats["arrays"]
    assert summary["topology"][0] == n_contigs
    assert summary["arrays"][0] == sc.stats["arrays"]
    assert summary["elements"][0] == 5 * sc.stats["arrays"]
    # contigs.parquet holds n_skip + n_contigs base rows + 1 row per contig with an array from _collect
    assert summary["contigs"][0] == n_skip + n_contigs + sc.stats["arrays"]


def test_scan_flushes_multi_row_group_buffers_at_exact_boundaries(tmp_path, monkeypatch):
    monkeypatch.setattr(crispr, "_search", _stub_search_factory(every_n=1000))
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    sc = crispr._Scan("java", tmp_path / "minced.jar", cpus=2, scratch=scratch, batch_bases=2_000_000, topology=True)

    n_skip = 1_000
    n_contigs = io.ROW_GROUP_ROWS + 7_120  # 130_000 -> crosses one ROW_GROUP_ROWS boundary
    seq = b"ACGT" * 50  # 200 bp
    for i in range(n_skip):
        sc.skip(i.to_bytes(16, "big"), "S0")
    for i in range(n_contigs):
        cu = (n_skip + i).to_bytes(16, "big")
        sc.contig(cu, seq, circular=(i % 11 == 0))
    sc.close(ok=True)

    summary = _assert_bounded_row_groups(sc.paths)
    assert summary["topology"] == (n_contigs, 2)
    assert pq.ParquetFile(sc.paths["topology"]).metadata.row_group(0).num_rows == io.ROW_GROUP_ROWS
    assert summary["contigs"][0] == n_skip + n_contigs + sc.stats["arrays"]
    assert summary["contigs"][1] == math.ceil(summary["contigs"][0] / io.ROW_GROUP_ROWS)
    assert pq.ParquetFile(sc.paths["contigs"]).metadata.row_group(0).num_rows == io.ROW_GROUP_ROWS
    assert summary["arrays"][0] == sc.stats["arrays"]
    assert summary["elements"][0] == 5 * sc.stats["arrays"]


def test_scan_close_not_ok_releases_cleanly(tmp_path):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    sc = crispr._Scan("java", tmp_path / "minced.jar", cpus=1, scratch=scratch, batch_bases=100_000, topology=True)
    sc.skip((1).to_bytes(16, "big"), "S0")
    sc.contig((2).to_bytes(16, "big"), b"ACGT" * 50, circular=False)
    sc.close(ok=False)
    assert all(p.exists() for p in sc.paths.values())
