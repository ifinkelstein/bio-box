"""Ring R1 tests for card C-58: single-cell SAGs (GORG), genome-reduced endosymbionts (SYM1),
post-R226 prokaryotes, and niche-habitat MAG compendia (DEC-69).
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
import random
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-58")

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"


def _load_c58():
    path = BSP / "c58_sags_symbionts_niche.py"
    spec = importlib.util.spec_from_file_location("c58_sags_symbionts_niche", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c58_sags_symbionts_niche"] = mod
    spec.loader.exec_module(mod)
    return mod


_SENSE_CODONS = [
    a + b + d for a in "ACGT" for b in "ACGT" for d in "ACGT" if (a + b + d) not in ("TAA", "TAG", "TGA")
]


def _make_coding_genome(seed: int, n_genes: int = 120, gene_aa: int = 280, *, internal_tga: bool = False) -> bytes:
    """Build a synthetic bacterial contig with balanced sense codons, optionally with internal TGA codons."""
    rng = random.Random(seed)
    parts: list[str] = ["A" * 40]
    for _ in range(n_genes):
        spacer = "".join(rng.choice("ACGT") for _ in range(7))
        codons = [rng.choice(_SENSE_CODONS) for _ in range(gene_aa - 2)]
        if internal_tga:
            step = max(1, len(codons) // 9)
            for k in range(1, 9):
                codons[min(len(codons) - 1, k * step)] = "TGA"
        intergenic = "".join(rng.choice("AT") for _ in range(25)) + "AGGAGG" + spacer
        parts.append(intergenic + "ATG" + "".join(codons) + "TAA")
    parts.append("T" * 40)
    return "".join(parts).encode("ascii")


def _write_ncbi_summary(path: Path, rows: list[dict[str, str]]) -> None:
    cols = [
        "assembly_accession", "bioproject", "biosample", "wgs_master", "refseq_category",
        "taxid", "species_taxid", "organism_name", "infraspecific_name", "isolate",
        "version_status", "assembly_level", "release_type", "genome_rep", "seq_rel_date",
        "asm_name", "asm_submitter", "gbrs_paired_asm", "paired_asm_comp", "ftp_path",
        "excluded_from_refseq", "relation_to_type_material", "asm_not_live_date",
        "assembly_type", "group", "genome_size", "genome_size_ungapped", "gc_percent",
        "replicon_count", "scaffold_count", "contig_count", "annotation_provider",
        "annotation_name", "annotation_date", "total_gene_count", "protein_coding_gene_count",
        "non_coding_gene_count", "pubmed_id",
    ]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("#   See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt\n")
        fh.write("#" + "\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "na")) for c in cols) + "\n")


def test_parse_and_diff_ncbi_assembly_summary_separates_symbionts_and_post_r226(tmp_path):
    c58 = _load_c58()
    refseq = tmp_path / "refseq_bac.txt"
    genbank = tmp_path / "genbank_bac.txt"
    _write_ncbi_summary(refseq, [
        # Already in GTDB R226 (by paired GCA)
        {
            "assembly_accession": "GCF_000005845.2", "taxid": "511145", "species_taxid": "562",
            "organism_name": "Escherichia coli str. K-12", "version_status": "latest",
            "assembly_level": "Complete Genome", "genome_rep": "Full", "seq_rel_date": "2013/09/26",
            "gbrs_paired_asm": "GCA_000005845.2", "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2",
            "group": "bacteria", "genome_size": "4641652", "contig_count": "1",
        },
        # Ultra-reduced insect endosymbiont (< 200 kb, absent from GTDB R226)
        {
            "assembly_accession": "GCF_000010525.1", "taxid": "387662", "species_taxid": "387662",
            "organism_name": "Candidatus Carsonella ruddii PV", "version_status": "latest",
            "assembly_level": "Complete Genome", "genome_rep": "Full", "seq_rel_date": "2006/10/12",
            "gbrs_paired_asm": "GCA_000010525.1", "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/010/525/GCF_000010525.1_ASM1052v1",
            "group": "bacteria", "genome_size": "159662", "contig_count": "1",
        },
        # Second strain of same species_taxid (lower quality / older -> deduplicated)
        {
            "assembly_accession": "GCF_000999995.1", "taxid": "387662", "species_taxid": "387662",
            "organism_name": "Candidatus Carsonella ruddii str. B", "version_status": "latest",
            "assembly_level": "Contig", "genome_rep": "Full", "seq_rel_date": "2020/01/01",
            "gbrs_paired_asm": "GCA_000999995.1", "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/999/995/GCF_000999995.1_x",
            "group": "bacteria", "genome_size": "155000", "contig_count": "4",
        },
        # Post-R226 standard prokaryote (>= 1 Mbp, new species)
        {
            "assembly_accession": "GCF_999000015.1", "taxid": "9990001", "species_taxid": "9990001",
            "organism_name": "Novelibacter postfreeze", "version_status": "latest",
            "assembly_level": "Complete Genome", "genome_rep": "Full", "seq_rel_date": "2025/08/15",
            "gbrs_paired_asm": "GCA_999000015.1", "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/999/000/015/GCF_999000015.1_x",
            "group": "bacteria", "genome_size": "3450000", "contig_count": "1",
        },
    ])
    _write_ncbi_summary(genbank, [
        # Duplicate of GCF_000010525.1 via paired GCA -> must not double-count
        {
            "assembly_accession": "GCA_000010525.1", "taxid": "387662", "species_taxid": "387662",
            "organism_name": "Candidatus Carsonella ruddii PV", "version_status": "latest",
            "assembly_level": "Complete Genome", "genome_rep": "Full", "seq_rel_date": "2006/10/12",
            "gbrs_paired_asm": "GCF_000010525.1", "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/010/525/GCA_000010525.1_ASM1052v1",
            "group": "bacteria", "genome_size": "159662", "contig_count": "1",
        },
        # Excluded from GTDB: < 1 Mbp obligate endosymbiont in GenBank only
        {
            "assembly_accession": "GCA_000022225.1", "taxid": "555111", "species_taxid": "555111",
            "organism_name": "Candidatus Hodgkinia cicadicola DSEM", "version_status": "latest",
            "assembly_level": "Complete Genome", "genome_rep": "Full", "seq_rel_date": "2009/07/01",
            "gbrs_paired_asm": "na", "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/022/225/GCA_000022225.1_ASM2222v1",
            "group": "bacteria", "genome_size": "143795", "contig_count": "1",
        },
        # Suppressed / non-latest assembly -> ignored
        {
            "assembly_accession": "GCA_000033335.1", "taxid": "555222", "species_taxid": "555222",
            "organism_name": "Candidatus Tremblaya princeps", "version_status": "suppressed",
            "assembly_level": "Complete Genome", "genome_rep": "Full", "seq_rel_date": "2011/01/01",
            "gbrs_paired_asm": "na", "ftp_path": "na",
            "group": "bacteria", "genome_size": "138927", "contig_count": "1",
        },
    ])

    gtdb_accs = {"GCA_000005845", "GCF_000005845"}
    globdb_accs = set()
    res = c58.diff_ncbi_summaries(
        [("refseq", "bacteria", refseq), ("genbank", "bacteria", genbank)],
        gtdb_base_accs=gtdb_accs,
        globdb_base_accs=globdb_accs,
    )
    sym_accs = [r["assembly_accession"] for r in res["sym1_candidates"]]
    post_accs = [r["assembly_accession"] for r in res["post_r226_species_reps"]]
    assert sym_accs == ["GCF_000010525.1", "GCA_000022225.1"]
    assert post_accs == ["GCF_999000015.1"]
    assert res["counts"]["in_gtdb_or_globdb"] == 1
    assert res["counts"]["suppressed_skipped"] == 1


def test_screen_genetic_code_distinguishes_code11_code4_and_code25():
    c58 = _load_c58()
    seq_11 = _make_coding_genome(11, n_genes=130, gene_aa=280, internal_tga=False)
    seq_uga = _make_coding_genome(17, n_genes=130, gene_aa=280, internal_tga=True)

    r11 = c58.screen_genetic_code([seq_11], organism_name="Buchnera aphidicola")
    assert r11["genetic_code"] == 11

    r4 = c58.screen_genetic_code([seq_uga], organism_name="Candidatus Hodgkinia cicadicola")
    assert r4["genetic_code"] == 4
    assert r4["tested_alt"] is True
    assert r4["tga_per_kb"] > 0

    r25 = c58.screen_genetic_code([seq_uga], organism_name="Candidatus Gracilibacteria bacterium")
    assert r25["genetic_code"] == 25
    assert r25["tested_alt"] is True

    # Tiny < 100 kb genome must not crash pyrodigal training
    tiny = _make_coding_genome(23, n_genes=30, gene_aa=200, internal_tga=True)
    rtiny = c58.screen_genetic_code([tiny], organism_name="Candidatus Nasuia deltocephalinicola")
    assert rtiny["genetic_code"] in (4, 11, 25)


def test_parse_gorg_sag_id_and_stage_sag1_deterministic_and_dedup(tmp_path):
    c58 = _load_c58()
    assert c58.parse_gorg_sag_id("AG-313-A02_NODE_1_length_12345_cov_10.2") == "AG-313-A02"
    assert c58.parse_gorg_sag_id("AH-135-C04_contig_7") == "AH-135-C04"

    s1 = _make_coding_genome(101, n_genes=125, gene_aa=270, internal_tga=False)
    s2 = _make_coding_genome(102, n_genes=125, gene_aa=270, internal_tga=False)
    s3 = _make_coding_genome(103, n_genes=125, gene_aa=270, internal_tga=True)

    src_fa = tmp_path / "GORG_v1.fasta.gz"
    with gzip.open(src_fa, "wb") as fh:
        fh.write(b">AG-313-A02_NODE_1\n" + s1 + b"\n")
        fh.write(b">AG-313-A02_NODE_2\n" + s2 + b"\n")
        fh.write(b">AG-313-A03_NODE_1\n" + s3 + b"\n")
        # Identical genome to AG-313-A03 -> repeated_records.tsv under DEC-27
        fh.write(b">AG-313-A04_NODE_1\n" + s3 + b"\n")

    # Earlier stage owns s1
    owner_dir = tmp_path / "owner_P1"
    owner_dir.mkdir()
    from biobox import ids
    (owner_dir / "genomes.json").write_text(json.dumps([
        {"contig_uids": [ids.contig_uid(s1).hex()], "genetic_code": 11, "kind": "prokaryote",
         "circular": [], "assembly_quality": "complete"}
    ]))

    out1 = tmp_path / "out1"
    out2 = tmp_path / "out2"
    rep1 = c58.stage_sag1(
        src_fa, out1,
        owner_inputs={"P1": owner_dir},
        retrieved_at="2026-11-24T00:00:00Z",
    )
    rep2 = c58.stage_sag1(
        src_fa, out2,
        owner_inputs={"P1": owner_dir},
        retrieved_at="2026-11-24T00:00:00Z",
    )

    assert rep1["genomes_distinct"] == 2
    assert rep1["genomes_repeated"] == 1
    assert rep1["contigs_held_earlier"] == 1
    assert rep1["contigs_new"] == 2
    assert rep1["genetic_codes"]["11"] == 1
    assert rep1["genetic_codes"].get("4", 0) + rep1["genetic_codes"].get("25", 0) == 2

    # Byte-identical outputs across runs (AGENTS §4.3; objects.json holds the output path)
    for fname in ("sag1_gorg.fna.gz", "genomes.json", "records.tsv",
                  "repeated_records.tsv", "shared_records.tsv", "genetic_codes.tsv"):
        h1 = hashlib.sha256((out1 / "SAG1" / fname).read_bytes()).hexdigest()
        h2 = hashlib.sha256((out2 / "SAG1" / fname).read_bytes()).hexdigest()
        assert h1 == h2, f"{fname} not deterministic"

    objs = json.loads((out1 / "SAG1" / "objects.json").read_text())
    assert objs[0]["redistributable"] is False


def test_stage_sym1_assigns_verified_codes_and_preserves_determinism(tmp_path):
    c58 = _load_c58()
    fa_dir = tmp_path / "sym_fna"
    fa_dir.mkdir()
    s_buch = _make_coding_genome(201, n_genes=130, gene_aa=280, internal_tga=False)
    s_hodg = _make_coding_genome(202, n_genes=130, gene_aa=280, internal_tga=True)
    with gzip.open(fa_dir / "GCF_000010525.1.fna.gz", "wb") as fh:
        fh.write(b">NC_008512.1 Candidatus Carsonella ruddii PV, complete genome\n" + s_buch + b"\n")
    with gzip.open(fa_dir / "GCA_000022225.1.fna.gz", "wb") as fh:
        fh.write(b">CP001226.1 Candidatus Hodgkinia cicadicola DSEM, complete genome\n" + s_hodg + b"\n")

    selected = [
        {
            "assembly_accession": "GCF_000010525.1", "species_taxid": "387662", "taxid": "387662",
            "organism_name": "Candidatus Carsonella ruddii PV", "assembly_level": "Complete Genome",
            "group": "bacteria", "genome_size": len(s_buch), "seq_rel_date": "2006/10/12",
            "fasta_path": str(fa_dir / "GCF_000010525.1.fna.gz"),
        },
        {
            "assembly_accession": "GCA_000022225.1", "species_taxid": "555111", "taxid": "555111",
            "organism_name": "Candidatus Hodgkinia cicadicola DSEM", "assembly_level": "Complete Genome",
            "group": "bacteria", "genome_size": len(s_hodg), "seq_rel_date": "2009/07/01",
            "fasta_path": str(fa_dir / "GCA_000022225.1.fna.gz"),
        },
    ]

    out1 = tmp_path / "sym_out1"
    out2 = tmp_path / "sym_out2"
    rep1 = c58.stage_sym1(selected, out1, retrieved_at="2026-11-24T00:00:00Z")
    rep2 = c58.stage_sym1(selected, out2, retrieved_at="2026-11-24T00:00:00Z")

    assert rep1["genomes_distinct"] == 2
    assert rep1["genetic_codes"] == {"11": 1, "4": 1}
    genomes = json.loads((out1 / "SYM1" / "genomes.json").read_text())
    assert sorted(g["genetic_code"] for g in genomes) == [4, 11]
    assert all(g["kind"] == "prokaryote" for g in genomes)
    objs = json.loads((out1 / "SYM1" / "objects.json").read_text())
    assert objs[0]["redistributable"] is True

    for fname in ("sym1_endosymbionts.fna.gz", "genomes.json", "sym1_metadata.tsv", "genetic_codes.tsv"):
        assert (out1 / "SYM1" / fname).read_bytes() == (out2 / "SYM1" / fname).read_bytes()


def test_niche_compendia_register_aligns_with_c66():
    c58 = _load_c58()
    reg = c58.niche_compendia_register()
    stages = [r["stage"] for r in reg]
    assert stages[0] == "NX1"
    assert len(stages) >= 4
    # NX1 must align with C-66 (MEER Hadal + eHOMD + mats/serpentinizing/sponges/lichens/rhizosphere)
    assert "MEER" in reg[0]["compendia"]
    assert "eHOMD" in reg[0]["compendia"]
    # None of the registered niche compendia may duplicate GlobDB r226's already-absorbed catalogues as sole target
    assert all(r["globdb_overlap"] != "100% absorbed" for r in reg)
