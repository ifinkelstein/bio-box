"""Ring R2: IGRs on micro-S0 (R-05; C-03 and C-05 for the igr table; LY-08)."""

from __future__ import annotations

from pathlib import Path

import pyarrow.parquet as pq
import pytest
import yaml

from biobox import catalog, contracts, fasta, ids
from biobox.l0seq import SeqStore, SourceObject, ingest
from biobox.l1 import igr
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, run as call_genes

pytestmark = pytest.mark.r2

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())


def _build(root: Path) -> set[bytes]:
    ingest([SourceObject(FIX / o["file"], "s", "1") for o in CONTRACT["objects"]], root, "S0", chunk_size=4096, row_group_chunks=16)
    genomes, circular = [], set()
    for o in CONTRACT["objects"]:
        cu = ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)
        if o["topology"] == "circular":
            circular.add(cu)
        genomes.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"], circular=(cu,) if cu in circular else ()))
    call_genes(root, "S0", genomes)
    igr.run(root, "S0", circular)
    return circular


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = tmp_path_factory.mktemp("igr")
    circular = _build(root)
    rows = pq.read_table(root / f"l0seq/igr/igr_version={igr.igr_version()}/stage=S0/part-0.parquet").to_pylist()
    return root, circular, rows


def _slots(root):
    occ = pq.read_table(root / f"l1/protein_occurrence/gene_call_version={GENE_CALL_VERSION}/stage=S0/part-0.parquet")
    return {(r["contig_uid"], r["slot_ordinal"]): r for r in occ.to_pylist()}


def test_r05_igrs_lie_between_their_flanking_slots(built):
    root, circular, rows = built
    slots = _slots(root)
    store = SeqStore(root)
    assert rows
    for r in rows:
        n = store.length_of(store._default[r["contig_uid"]])
        if r["left_slot"] is not None:
            assert slots[(r["contig_uid"], r["left_slot"])]["end"] % n == r["start"] % n
        if r["right_slot"] is not None:
            assert slots[(r["contig_uid"], r["right_slot"])]["start"] == r["end"] % n or \
                slots[(r["contig_uid"], r["right_slot"])]["start"] == r["end"]
    # no slot overlaps any IGR interior
    by_contig = {}
    for (cu, _), s in slots.items():
        by_contig.setdefault(cu, []).append(s)
    for r in rows:
        n = store.length_of(store._default[r["contig_uid"]])
        span = {p % n for p in range(r["start"], r["end"])}
        for s in by_contig.get(r["contig_uid"], []):
            assert not span & {p % n for p in range(s["start"], s["end"])}


def test_r05_igr_slices_equal_the_source(built, tmp_path):
    root, circular, rows = built
    store = SeqStore(root)
    con = catalog.open_catalog(catalog.generate_catalog(catalog.build_manifest(root, "x"), tmp_path / "c.duckdb"))
    src = {ids.contig_uid(fasta.read(FIX / o["file"])[0].seq): fasta.read(FIX / o["file"])[0].seq for o in CONTRACT["objects"]}
    for r in rows[::7]:
        seq = src[r["contig_uid"]]
        want = (seq + seq)[r["start"]:r["end"]]
        assert store.get(r["contig_uid"], r["start"], r["end"], 1, circular=r["contig_uid"] in circular) == want
        got = con.execute("SELECT seq_get_circular(?, ?, ?, 1)", [r["contig_uid"], r["start"], r["end"]]).fetchall()[0][0]
        assert got.encode() == want
    con.close()


def test_igr_table_passes_contracts(built, tmp_path):
    root, _, _ = built
    man = catalog.build_manifest(root, "x")
    con = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / "c.duckdb"))
    assert contracts.check_foreign_keys(con) == []
    assert contracts.check_coordinates(con) == []
    con.close()
    assert contracts.check_sort_order(man) == []


def test_ly08_igr_rewrite_is_byte_identical(built, tmp_path):
    root, circular, _ = built
    out = root / f"l0seq/igr/igr_version={igr.igr_version()}/stage=S0/part-0.parquet"
    before = out.read_bytes()
    igr.run(root, "S0", circular)
    assert out.read_bytes() == before
