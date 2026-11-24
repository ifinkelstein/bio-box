"""Ring R1: curated integrative MGEs, transposons, IS elements, integrons, and prophages (card C-57, DEC-69).

Tests constructed inputs for `scripts/stress/build_stage_pod/c57_integrative_mges.py`:
- Byte-identical staging with NO 5 kb floor (0.7–3 kb IS elements and integron cassettes stay whole)
- DEC-65 topology handling (public curated MGEs never trimmed, staged byte-identical; optional metagenomic DTR trim)
- DEC-51 / DEC-63 (1) redistributable split (refuses mixed redistributable in one stage; clean MGE1/MGE2 split)
- DEC-26 / DEC-27 contig_uid deduplication within-stage and ownership against earlier root + staged inputs
- Deterministic MGE class / family / host accession TSV sidecar and byte-identical reruns
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pyarrow as pa
import pytest

from biobox import fasta, ids, io
from biobox.registry.schemas import table

pytestmark = pytest.mark.covers("C-57")  # ring R1 by directory

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"


def _mod():
    path = BSP / "c57_integrative_mges.py"
    spec = importlib.util.spec_from_file_location("c57_integrative_mges", path)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c57_integrative_mges"] = mod
    spec.loader.exec_module(mod)
    return mod


def _seq(seed: str, length: int) -> bytes:
    """Deterministic non-repeating nucleotide sequence of exact length."""
    alphabet = b"ACGT"
    out = bytearray()
    h = hashlib.sha256(seed.encode()).digest()
    i = 0
    while len(out) < length:
        if i % 32 == 0:
            h = hashlib.sha256(seed.encode() + i.to_bytes(4, "little")).digest()
        b = h[i % 32]
        out.append(alphabet[b & 3])
        out.append(alphabet[(b >> 2) & 3])
        out.append(alphabet[(b >> 4) & 3])
        out.append(alphabet[(b >> 6) & 3])
        i += 1
    return bytes(out[:length])


def _write_fasta(path: Path, records: list[tuple[str, bytes]]) -> None:
    with open(path, "wb") as fh:
        for header, seq in records:
            fh.write(b"\n".join(fasta.format_record(header, seq, 80)) + b"\n")


def _write_root_stage(root: Path, stage: str, seqs: list[bytes]) -> None:
    spec = table("contig_index")
    uids = [ids.contig_uid(s) for s in seqs]
    digests = [hashlib.sha256(s.upper()).digest() for s in seqs]
    lens = [len(s) for s in seqs]
    tbl = pa.table(
        {
            "contig_uid": pa.array(uids, type=pa.binary(16)),
            "raw_digest": pa.array(digests, type=pa.binary(32)),
            "contig_len": pa.array(lens, type=pa.int64()),
            "is_default_variant": pa.array([True] * len(seqs), type=pa.bool_()),
        },
        schema=spec.schema,
    )
    out = root / "l0seq" / "contig_index" / f"stage={stage}" / "part-0.parquet"
    io.write_table(tbl, spec, out)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))


def test_header_parsers_extract_class_family_and_host():
    m = _mod()

    # ICEberg 3.0
    ice = m.parse_iceberg_header("ICEberg|1234|ICE_SXT|NZ_CP012345.1:100..50000", "ICE")
    assert ice["mge_class"] == "ICE"
    assert ice["host_accession"] == "NZ_CP012345.1"

    ime = m.parse_iceberg_header("ICEberg|99|IME_foo|NC_000913.3", "IME")
    assert ime["mge_class"] == "IME"
    assert ime["host_accession"] == "NC_000913.3"

    # ISfinder-sequences: >(IS name)_(IS group)_(IS family)
    isf = m.parse_isfinder_header(
        "IS605_IS605_IS200/IS605",
        {"IS605": {"family": "IS200/IS605", "origin": "Helicobacter pylori", "accession": "U60177"}},
    )
    assert isf["mge_class"] == "IS"
    assert isf["family"] == "IS200/IS605"
    assert isf["host_accession"] == "U60177"

    # MGEdb (MobileElementFinder)
    mge_tn = m.parse_mgedb_record("Tn7", {"type": "unit transposon", "family": "Tn7", "accession": "KX117211"})
    assert mge_tn["mge_class"] == "Tn"
    assert mge_tn["family"] == "Tn7"
    assert mge_tn["host_accession"] == "KX117211"

    mge_is = m.parse_mgedb_record("IS110", {"type": "insertion sequence", "family": "IS110", "accession": "Y00434"})
    assert mge_is["mge_class"] == "IS"
    assert mge_is["family"] == "IS110"

    # GEMs prophage
    pro = m.parse_gems_prophage(
        "2013843002_3_DCKB1_C1262",
        {"bin_id": "2013843002_3", "contig_id": "2013843002.a:DCKB1_C1262", "contig_start": "3756", "contig_stop": "14626"},
    )
    assert pro["mge_class"] == "prophage"
    assert pro["host_accession"] == "2013843002.a:DCKB1_C1262"


def test_no_5kb_floor_and_public_byte_identical_staging(tmp_path):
    """Short 0.7-3 kb IS elements and integron cassettes stay whole; public records are never trimmed."""
    m = _mod()
    # 720 bp IS element with an exact 40 bp terminal repeat (must NOT be trimmed for public curated MGE!)
    rep = _seq("rep40", 40)
    is_seq = rep + _seq("is_core", 640) + rep  # 720 bp
    cassette_seq = _seq("integron_cassette", 1150)  # 1.15 kb
    ice_seq = _seq("ice_large", 12000)  # 12 kb

    is_fa = tmp_path / "is.fna"
    int_fa = tmp_path / "integrons.fna"
    ice_fa = tmp_path / "ice.fna"
    _write_fasta(is_fa, [("IS110_IS110_IS110", is_seq)])
    _write_fasta(int_fa, [("intI1_cassette_1", cassette_seq)])
    _write_fasta(ice_fa, [("ICEberg|1|ICE_test|NC_000913.3", ice_seq)])

    out = tmp_path / "MGE1"
    sources = [
        {
            "source_name": "isfinder_sequences",
            "source_release": "2020-10-06",
            "license": "open",
            "redistributable": True,
            "fasta": is_fa,
            "output_name": "isfinder_IS.fna",
            "parser": "isfinder",
            "kind": "plasmid",
            "mge_class": "IS",
        },
        {
            "source_name": "parnanen_mgedb",
            "source_release": "2021-08-16",
            "license": "CC BY 4.0",
            "redistributable": True,
            "fasta": int_fa,
            "output_name": "parnanen_integrons.fna",
            "parser": "parnanen",
            "kind": "plasmid",
            "mge_class": "integron",
        },
        {
            "source_name": "iceberg3",
            "source_release": "3.0",
            "license": "CC BY 4.0",
            "redistributable": True,
            "fasta": ice_fa,
            "output_name": "iceberg3_ICE.fna",
            "parser": "iceberg",
            "kind": "plasmid",
            "mge_class": "ICE",
        },
    ]

    rep_doc = m.stage_mges(
        sources,
        out,
        root=None,
        root_stages=[],
        owner_inputs={},
        retrieved="2026-11-23T00:00:00Z",
        stage_name="MGE1",
    )

    # All 3 records kept (including 720 bp IS and 1150 bp integron cassette: NO 5 kb floor)
    genomes = json.loads((out / "genomes.json").read_text())
    assert len(genomes) == 3
    assert rep_doc["MGE1"]["genomes"] == 3
    assert rep_doc["MGE1"]["trims"] == 0

    # Public IS record was NOT trimmed despite having a 40 bp terminal repeat
    staged_is = list(fasta.iter_records(out / "isfinder_IS.fna"))
    assert len(staged_is) == 1
    assert staged_is[0].seq == is_seq
    assert len(staged_is[0].seq) == 720
    assert genomes[0]["circular"] == []

    # Sidecar TSV has deterministic rows with class/family/host
    sidecar = _read_tsv(out / "mge_metadata.tsv")
    assert len(sidecar) == 3
    assert [r["mge_class"] for r in sidecar] == ["IS", "integron", "ICE"]
    assert sidecar[0]["family"] == "IS110"
    assert sidecar[2]["host_accession"] == "NC_000913.3"
    assert [int(r["length"]) for r in sidecar] == [720, 1150, 12000]


def test_refuses_mixed_redistributable_without_split_and_splits_cleanly(tmp_path):
    """DEC-51 / DEC-63 (1): refuse mixing redistributable flags in one stage; split cleanly into MGE1 + MGE2."""
    m = _mod()
    pub_seq = _seq("pub_tn", 2400)
    jgi_seq = _seq("jgi_prophage", 15000)

    pub_fa = tmp_path / "pub.fna"
    jgi_fa = tmp_path / "jgi.fna"
    _write_fasta(pub_fa, [("Tn7_test", pub_seq)])
    _write_fasta(jgi_fa, [("2013843002_3_DCKB1_C1262", jgi_seq)])

    sources = [
        {
            "source_name": "plasann_tncentral",
            "source_release": "1.1.6",
            "license": "CC BY 4.0",
            "redistributable": True,
            "fasta": pub_fa,
            "output_name": "tncentral.fna",
            "parser": "tncentral",
            "kind": "plasmid",
            "mge_class": "Tn",
        },
        {
            "source_name": "jgi_gems_nersc",
            "source_release": "2021-01",
            "license": "JGI Data Utilization Policy",
            "redistributable": False,
            "fasta": jgi_fa,
            "output_name": "gems_virsorter_prophages.fna",
            "parser": "gems_prophage",
            "kind": "virus",
            "mge_class": "prophage",
        },
    ]

    with pytest.raises(m.Refused, match="redistributable"):
        m.stage_mges(
            sources,
            tmp_path / "mixed",
            root=None,
            root_stages=[],
            owner_inputs={},
            retrieved="2026-11-23T00:00:00Z",
            stage_name="MGE1",
            split_by_redistributable=False,
        )

    out1 = tmp_path / "MGE1"
    out2 = tmp_path / "MGE2"
    rep = m.stage_mges(
        sources,
        out1,
        root=None,
        root_stages=[],
        owner_inputs={},
        retrieved="2026-11-23T00:00:00Z",
        stage_name="MGE1",
        split_by_redistributable=True,
        out_restricted=out2,
        restricted_stage_name="MGE2",
    )

    o1 = json.loads((out1 / "objects.json").read_text())
    o2 = json.loads((out2 / "objects.json").read_text())
    assert len(o1) == 1 and o1[0]["redistributable"] is True
    assert len(o2) == 1 and o2[0]["redistributable"] is False
    assert rep["MGE1"]["genomes"] == 1
    assert rep["MGE2"]["genomes"] == 1
    g2 = json.loads((out2 / "genomes.json").read_text())
    assert g2[0]["kind"] == "virus"


def test_dedup_within_stage_and_build_order_ownership_and_determinism(tmp_path):
    """DEC-26 ownership against root + staged inputs, DEC-27 within-stage dedup, and byte-identical rerun."""
    m = _mod()
    s_root = _seq("held_by_P1", 1800)
    s_staged = _seq("held_by_PL2", 2200)
    s_new = _seq("novel_ice", 4500)

    root = tmp_path / "root"
    _write_root_stage(root, "P1", [s_root])

    pl2_dir = tmp_path / "PL2"
    pl2_dir.mkdir()
    (pl2_dir / "genomes.json").write_text(
        json.dumps([{"contig_uids": [ids.contig_uid(s_staged).hex()], "genetic_code": 11, "kind": "plasmid", "circular": [], "assembly_quality": "complete"}])
    )

    fa = tmp_path / "mge.fna"
    # s_new appears twice to verify DEC-27 within-stage dedup
    _write_fasta(
        fa,
        [
            ("ICEberg|1|ICE_root|NC_1", s_root),
            ("ICEberg|2|ICE_pl2|NC_2", s_staged),
            ("ICEberg|3|ICE_new|NC_3", s_new),
            ("ICEberg|4|ICE_new_dup|NC_4", s_new),
        ],
    )

    sources = [
        {
            "source_name": "iceberg3",
            "source_release": "3.0",
            "license": "CC BY 4.0",
            "redistributable": True,
            "fasta": fa,
            "output_name": "iceberg3_ICE.fna",
            "parser": "iceberg",
            "kind": "plasmid",
            "mge_class": "ICE",
        }
    ]

    out_a = tmp_path / "run_a"
    out_b = tmp_path / "run_b"
    final_dir = tmp_path / "final_MGE1"

    rep_a = m.stage_mges(
        sources,
        out_a,
        root=root,
        root_stages=["P1"],
        owner_inputs={"PL2": pl2_dir},
        retrieved="2026-11-23T00:00:00Z",
        stage_name="MGE1",
        final_dir=final_dir,
    )
    rep_b = m.stage_mges(
        sources,
        out_b,
        root=root,
        root_stages=["P1"],
        owner_inputs={"PL2": pl2_dir},
        retrieved="2026-11-23T00:00:00Z",
        stage_name="MGE1",
        final_dir=final_dir,
    )

    # 4 records -> 3 distinct genomes (1 repeated)
    assert rep_a["MGE1"]["records"] == 4
    assert rep_a["MGE1"]["genomes"] == 3
    assert rep_a["MGE1"]["repeated_records"] == 1
    assert rep_a["MGE1"]["owned_by_stage"] == {"P1": 1, "PL2": 1}
    assert rep_a["MGE1"]["new_contigs"] == 1
    assert rep_a["MGE1"]["new_bases"] == 4500

    # Byte-identical determinism across reruns
    assert rep_a["MGE1"]["sha256"] == rep_b["MGE1"]["sha256"]
    for name in ("genomes.json", "objects.json", "records.tsv", "repeated_records.tsv", "shared_records.tsv", "trims.tsv", "mge_metadata.tsv"):
        assert (out_a / name).read_bytes() == (out_b / name).read_bytes()
