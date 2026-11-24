"""Ring R1 unit tests for card C-59: NCBI RefSeq Organellar genomes (O1)."""
from __future__ import annotations

import gzip
import importlib.util
import io
import json
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-59")

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"


def _mod():
    name = "c59_refseq_organelles"
    if name in sys.modules:
        return sys.modules[name]
    path = BSP / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None, f"missing {path}"
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _make_gbff_record(
    *,
    acc: str,
    seq: str,
    circular: bool = True,
    organelle: str = "plastid:chloroplast",
    organism: str = "Arabidopsis thaliana",
    taxonomy: str = "Eukaryota; Viridiplantae; Streptophyta.",
    definition: str = "Arabidopsis thaliana chloroplast, complete genome.",
    cds_entries: list[dict] | None = None,
    n_introns: int = 0,
) -> str:
    topo = "circular" if circular else "linear"
    lines = [
        f"LOCUS       {acc.split('.')[0]:<16} {len(seq)} bp    DNA     {topo}   PLN 12-NOV-2026",
        f"DEFINITION  {definition}",
        f"ACCESSION   {acc.split('.')[0]}",
        f"VERSION     {acc}",
        f"  ORGANISM  {organism}",
        f"            {taxonomy}",
        "FEATURES             Location/Qualifiers",
        f"     source          1..{len(seq)}",
        f'                     /organism="{organism}"',
    ]
    if organelle:
        lines.append(f'                     /organelle="{organelle}"')
    for i in range(n_introns):
        lines.append(f"     intron          {10 + i * 5}..{12 + i * 5}")
    for cds in cds_entries or []:
        lines.append(f"     CDS             {cds['location']}")
        if cds.get("gene"):
            lines.append(f'                     /gene="{cds["gene"]}"')
        if cds.get("transl_table") is not None:
            lines.append(f"                     /transl_table={cds['transl_table']}")
        if cds.get("product"):
            lines.append(f'                     /product="{cds["product"]}"')
        if cds.get("translation"):
            aa = cds["translation"]
            lines.append(f'                     /translation="{aa}"')
    lines.append("ORIGIN")
    for i in range(0, len(seq), 60):
        chunk = seq[i : i + 60].lower()
        groups = " ".join(chunk[j : j + 10] for j in range(0, len(chunk), 10))
        lines.append(f"{i + 1:>9} {groups}")
    lines.append("//")
    return "\n".join(lines) + "\n"


def _orf_dna(n_codons: int = 120) -> str:
    """Construct a valid table-11 ORF with a Shine-Dalgarno-like upstream flank."""
    # ATG + repeating GCT (Ala) / GAA (Glu) / CTG (Leu) / AAA (Lys) + TAA
    body = "GCTGAACTGAAA" * (n_codons // 4)
    return "AAGGAGGTGATATA" + "ATG" + body + "TAA" + "ACGT" * 20


def test_parse_refseq_organelle_gbff_extracts_transl_table_introns_and_homing_rt():
    m = _mod()
    seq_p = "ATG" + "GCT" * 50 + "TAA" + "C" * 50
    seq_m = "ATG" + "TGA" * 10 + "GCT" * 40 + "TAA" + "G" * 50
    gbff = _make_gbff_record(
        acc="NC_000932.1",
        seq=seq_p,
        circular=True,
        organelle="plastid:chloroplast",
        cds_entries=[
            {
                "location": "1..156",
                "gene": "matK",
                "transl_table": 11,
                "product": "maturase K",
                "translation": "M" + "A" * 50,
            },
            {
                "location": "join(10..50,90..150)",
                "gene": "rps12",
                "transl_table": 11,
                "product": "ribosomal protein S12",
                "translation": "M" + "A" * 30,
            },
        ],
        n_introns=2,
    ) + _make_gbff_record(
        acc="NC_001224.1",
        seq=seq_m,
        circular=True,
        organelle="mitochondrion",
        organism="Saccharomyces cerevisiae",
        definition="Saccharomyces cerevisiae mitochondrion, complete genome.",
        cds_entries=[
            {
                "location": "1..156",
                "gene": "ai1",
                "transl_table": 3,
                "product": "group II intron-encoded reverse transcriptase",
                "translation": "M" + "W" * 10 + "A" * 40,
            },
        ],
        n_introns=1,
    )
    recs = list(m.parse_refseq_organelle_gbff(io.StringIO(gbff)))
    assert len(recs) == 2

    p, mt = recs
    assert p["accession"] == "NC_000932.1"
    assert p["circular"] is True
    assert p["transl_table"] == 11
    assert p["organelle"] == "plastid:chloroplast"
    assert p["n_cds"] == 2
    assert p["n_spliced_cds"] == 1
    assert p["n_introns"] == 2
    assert p["n_homing_rt"] == 1
    assert p["seq"] == seq_p.encode("ascii")

    assert mt["accession"] == "NC_001224.1"
    assert mt["transl_table"] == 3
    assert mt["n_homing_rt"] == 1
    assert mt["n_spliced_cds"] == 0
    assert mt["n_introns"] == 1


def test_origin_wrapping_join_is_not_counted_as_spliced_intron():
    m = _mod()
    seq = "A" * 300
    gbff = _make_gbff_record(
        acc="NC_999999.1",
        seq=seq,
        circular=True,
        cds_entries=[
            {
                "location": "join(250..300,1..50)",
                "gene": "psbA",
                "transl_table": 11,
                "product": "photosystem II protein D1",
                "translation": "M" + "K" * 30,
            },
        ],
    )
    recs = list(m.parse_refseq_organelle_gbff(io.StringIO(gbff)))
    assert len(recs) == 1
    assert recs[0]["n_spliced_cds"] == 0


def test_partition_routes_table11_to_o1_and_non_table11_to_c59a():
    m = _mod()
    seq = _orf_dna(80)
    r11 = {
        "accession": "NC_001.1",
        "circular": True,
        "length": len(seq),
        "organelle": "plastid:chloroplast",
        "definition": "Chloroplast complete genome",
        "transl_table": 11,
        "n_cds": 5,
        "n_spliced_cds": 0,
        "n_introns": 0,
        "n_homing_rt": 1,
        "seq": seq.encode("ascii"),
    }
    r4 = dict(r11, accession="NC_002.1", organelle="mitochondrion", transl_table=4, n_spliced_cds=2)
    r_none = dict(r11, accession="NC_003.1", organelle="mitochondrion", transl_table=None, n_cds=0)

    b11 = m.partition_organelle_record(r11, source_compartment="plastid")
    b4 = m.partition_organelle_record(r4, source_compartment="mitochondrion")
    b_none = m.partition_organelle_record(r_none, source_compartment="mitochondrion")

    assert b11["bucket"] == "O1" and b11["compartment"] == "plastid"
    assert b4["bucket"] == "left_for_c59a" and b4["reason"] == "transl_table_4"
    assert "c59a" in b4["needed_capability"]
    assert b_none["bucket"] == "left_for_c59a" and b_none["reason"] == "no_transl_table"


def test_stage_o1_byte_identical_dedup_circular_and_deterministic(tmp_path):
    from biobox import ids

    m = _mod()
    seq_p1 = _orf_dna(100)
    seq_p2 = _orf_dna(120)
    seq_m11 = _orf_dna(90)
    seq_m2 = _orf_dna(85)

    plastid_gb = tmp_path / "plastid.1.genomic.gbff.gz"
    with gzip.open(plastid_gb, "wt", encoding="utf-8") as fh:
        fh.write(_make_gbff_record(
            acc="NC_100001.1", seq=seq_p1, circular=True, organelle="plastid:chloroplast",
            cds_entries=[{"location": "15..377", "transl_table": 11, "product": "maturase K", "translation": "M" + "A" * 50}],
        ))
        # Repeat of seq_p1 under a second accession (DEC-27)
        fh.write(_make_gbff_record(
            acc="NC_100002.1", seq=seq_p1, circular=True, organelle="plastid:chloroplast",
            cds_entries=[{"location": "15..377", "transl_table": 11, "product": "maturase K", "translation": "M" + "A" * 50}],
        ))
        fh.write(_make_gbff_record(
            acc="NC_100003.1", seq=seq_p2, circular=False, organelle="chromatophore",
            definition="Paulinella chromatophora chromatophore, complete genome",
            cds_entries=[{"location": "15..437", "transl_table": 11, "product": "hypothetical protein", "translation": "M" + "A" * 60}],
        ))

    mito_gb = tmp_path / "mitochondrion.1.genomic.gbff.gz"
    with gzip.open(mito_gb, "wt", encoding="utf-8") as fh:
        fh.write(_make_gbff_record(
            acc="NC_200001.1", seq=seq_m11, circular=True, organelle="mitochondrion",
            organism="Reclinomonas americana",
            definition="Reclinomonas americana mitochondrion, complete genome",
            cds_entries=[{"location": "15..347", "transl_table": 11, "product": "LAGLIDADG homing endonuclease", "translation": "M" + "A" * 40}],
        ))
        fh.write(_make_gbff_record(
            acc="NC_200002.1", seq=seq_m2, circular=True, organelle="mitochondrion",
            organism="Homo sapiens",
            definition="Homo sapiens mitochondrion, complete genome",
            cds_entries=[{"location": "15..332", "transl_table": 2, "product": "cytochrome c oxidase subunit I", "translation": "M" + "A" * 40}],
        ))

    # Fake earlier stage owning seq_p2
    prev_dir = tmp_path / "prev_PL1"
    prev_dir.mkdir()
    cu_p2 = ids.contig_uid(seq_p2.encode("ascii")).hex()
    (prev_dir / "genomes.json").write_text(json.dumps([{"contig_uids": [cu_p2], "genetic_code": 11, "kind": "plasmid"}]))

    sources_spec = [
        {"path": str(plastid_gb), "source_name": "ncbi_refseq_plastid_r238", "source_release": "238", "source_compartment": "plastid"},
        {"path": str(mito_gb), "source_name": "ncbi_refseq_mitochondrion_r238", "source_release": "238", "source_compartment": "mitochondrion"},
    ]

    out1 = tmp_path / "O1_a"
    out2 = tmp_path / "O1_b"
    rep1 = m.stage_o1(
        sources_spec, out1, root=None, root_stages=(), owner_inputs={"PL1": prev_dir},
        retrieved="2026-11-24T00:00:00Z", genome_kind="plasmid",
    )
    rep2 = m.stage_o1(
        sources_spec, out2, root=None, root_stages=(), owner_inputs={"PL1": prev_dir},
        retrieved="2026-11-24T00:00:00Z", genome_kind="plasmid",
    )

    # 3 distinct table-11 genomes (4 records, 1 repeated); 1 table-2 record left for C-59a
    assert rep1["o1_records"] == 4
    assert rep1["genomes"] == 3
    assert rep1["repeated_genomes"] == 1
    assert rep1["circular_genomes"] == 2
    assert rep1["held_contigs"] == 1
    assert rep1["new_contigs"] == 2
    assert rep1["left_for_c59a"] == 1

    # Byte-identical FASTA & deterministic staging outputs (DEC-65 (2) & AGENTS §4.3)
    for fname in ("genomes.json", "o1_metadata.tsv", "O1_left_for_c59a.tsv.gz", "repeated_records.tsv", "shared_records.tsv"):
        assert (out1 / fname).read_bytes() == (out2 / fname).read_bytes(), f"non-deterministic {fname}"

    # Verify objects.json redistributable=True and genomes.json genetic_code=11, kind="plasmid"
    objs = json.loads((out1 / "objects.json").read_text())
    assert all(o["redistributable"] is True for o in objs)
    genomes = json.loads((out1 / "genomes.json").read_text())
    assert all(g["genetic_code"] == 11 and g["kind"] == "plasmid" for g in genomes)

    # Verify left_for_c59a.tsv.gz records NC_200002.1 with transl_table=2
    with gzip.open(out1 / "O1_left_for_c59a.tsv.gz", "rt", encoding="utf-8") as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh]
    assert len(rows) == 2
    assert rows[1][0] == "NC_200002.1"
    assert "transl_table_2" in rows[1]


def test_benchmark_pyrodigal_organelles_recovers_table11_orf():
    m = _mod()
    seq = _orf_dna(160)
    aa = "M" + "AELK" * 40
    rec = {
        "accession": "NC_TEST.1",
        "compartment": "plastid",
        "circular": False,
        "transl_table": 11,
        "cds": [
            {
                "location": "15..497",
                "gene": "matK",
                "product": "maturase K",
                "transl_table": 11,
                "translation": aa,
                "is_spliced": False,
                "is_homing_rt": True,
            },
        ],
        "seq": seq.encode("ascii"),
    }
    bench = m.benchmark_pyrodigal_organelles([rec], modes=("meta", "gv-meta"))
    assert "plastid:table_11" in bench
    assert bench["plastid:table_11"]["meta"]["cds_total"] == 1
    assert bench["plastid:table_11"]["meta"]["homing_rt_total"] == 1
    assert bench["plastid:table_11"]["meta"]["cds_overlap_pct"] == 100.0
    assert bench["plastid:table_11"]["meta"]["homing_rt_overlap_pct"] == 100.0
