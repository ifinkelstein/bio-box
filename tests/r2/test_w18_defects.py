"""Card W-18: small defects from the trial search on rN.1 (catalog items 1 and 6).

(1) `open_catalog` must not default DuckDB's temporary directory to a path under the data root.
(6) On a release built before card B-27 (no `contig_label` files) the catalog's `contig_label` view reads
    `contig_occurrence` instead of being empty.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest
import yaml

from biobox import catalog, fasta
from biobox.l0seq import SourceObject, ingest

pytestmark = [pytest.mark.r2, pytest.mark.covers("W-18")]

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("w18")
    lam_path = FIX / CONTRACT["objects"][0]["file"]
    lam = fasta.read(lam_path)[0]
    src = tmp / "src"
    src.mkdir()
    dup = src / "lambda_again.fna"
    dup.write_bytes(b"\n".join(fasta.format_record("lambda_again", lam.seq, 70)) + b"\n")
    root = tmp / "s"
    ingest([SourceObject(lam_path, CONTRACT["source_name"], CONTRACT["source_release"], CONTRACT["license"], True),
            SourceObject(dup, "other_source", "1", "INSDC (public)", True, precedence=200)],
           root, "S0", shard_bases=200_000)
    return root, tmp


def test_w18_open_catalog_temp_directory_is_not_under_the_data_root(world):
    root, tmp = world
    man = catalog.build_manifest(root, "x", ["S0"])
    cat = catalog.generate_catalog(man, tmp / "c1.duckdb")
    c = catalog.open_catalog(cat)
    try:
        temp = Path(c.execute("SELECT current_setting('temp_directory')").fetchone()[0]).resolve()
    finally:
        c.close()
    assert not temp.is_relative_to(root.resolve()), temp
    assert temp != Path(f"{cat.resolve()}.tmp"), temp  # the old default: beside the catalog file, under a root
    assert temp.is_relative_to(Path(tempfile.gettempdir()).resolve()), temp


def test_w18_open_catalog_honours_an_explicit_temp_directory(world):
    root, tmp = world
    cat = catalog.generate_catalog(catalog.build_manifest(root, "x", ["S0"]), tmp / "c2.duckdb")
    want = tmp / "spill"
    c = catalog.open_catalog(cat, temp_directory=want)
    try:
        assert Path(c.execute("SELECT current_setting('temp_directory')").fetchone()[0]) == want
    finally:
        c.close()


def test_w18_contig_label_view_without_label_files_reads_contig_occurrence(world):
    """A release whose manifest has no `contig_label` file (built before B-27): the view has the rows the
    registry table would have, taken from `contig_occurrence`, one per contig, first record by
    (l0_object_sha256, record_ordinal), with the number of records."""
    root, tmp = world
    man = catalog.build_manifest(root, "x", ["S0"])
    assert any(f["table"] == "contig_label" for f in man["files"])  # the fixture has them: this is the oracle
    full = catalog.open_catalog(catalog.generate_catalog(man, tmp / "full.duckdb"))
    want = full.execute("SELECT contig_uid, accession, original_header, n_records, l0_object_sha256, record_ordinal, "
                        "n_occurrences FROM contig_label ORDER BY contig_uid").fetchall()
    full.close()
    assert want
    old = dict(man, files=[f for f in man["files"] if f["table"] != "contig_label"])
    c = catalog.open_catalog(catalog.generate_catalog(old, tmp / "old.duckdb"))
    try:
        got = c.execute("SELECT contig_uid, accession, original_header, n_records, l0_object_sha256, record_ordinal, "
                        "n_occurrences FROM contig_label ORDER BY contig_uid").fetchall()
        assert got == want
        assert got[0][3] == 2 and got[0][6] == 2
        assert c.execute("SELECT count(*) FROM contig_label WHERE contig_uid = ?", [want[0][0]]).fetchone()[0] == 1
    finally:
        c.close()
