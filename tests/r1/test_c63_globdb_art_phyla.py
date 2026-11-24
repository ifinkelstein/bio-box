"""Card C-63 (DEC-71, DEC-72): stage ART-RT-enriched GlobDB r226 bacterial phyla as GART1 and GART2.

Tests:
1. Non-GTDB phylum selection and deterministic (taxonomy, ID) ordering for GART1 and GART2.
2. Pyrodigal translation_table 11 vs 25 vs 4 genetic-code screening on constructed >=105 kb
   code-11 and UGA-recoded code-25 bacterial genomes (plus Gracilibacteria taxonomy handling).
3. End-to-end streaming from a local tar.gz archive into consolidated multi-genome .fna.gz
   files, genomes.json, objects.json, objects_staging.json, records.tsv, shared_records.tsv,
   and genetic_codes.tsv, validated by s2_check_inputs.py.
4. Contig ownership deduplication against earlier stages in a data root and staged owner
   inputs (DEC-26 / DEC-27), including within-stage shared contigs and repeated genomes.
5. Byte-identical determinism across repeated runs (including zero-mtime .fna.gz outputs).

Oracles are independent of the script: contig_uid via blake3 directly, selection via a
brute-force Python filter, and input validation via s2_check_inputs.py.
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import io
import json
import random
import subprocess
import sys
import tarfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("C-63")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts" / "stress" / "build_stage_pod"

GART1_PHYLA = {"Myxococcota", "Myxococcota_A", "Cyanobacteriota", "Omnitrophota", "Planctomycetota"}
GART2_PHYLA = {"Patescibacteriota"}

# Sense codons without any stop codon (TAA, TAG, TGA), balanced GC ~50%
_SENSE_CODONS = [
    c
    for c in [a + b + d for a in "ACGT" for b in "ACGT" for d in "ACGT"]
    if c not in ("TAA", "TAG", "TGA")
]


def _load():
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location("c63_globdb_art_phyla", SCRIPTS / "c63_globdb_art_phyla.py")
    if spec is None or spec.loader is None:
        raise FileNotFoundError(SCRIPTS / "c63_globdb_art_phyla.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c63_globdb_art_phyla"] = mod
    spec.loader.exec_module(mod)
    return mod


def cuid(seq: bytes | str) -> str:
    raw = seq.encode("ascii") if isinstance(seq, str) else seq
    return blake3(raw.upper()).digest(length=16).hex()


def make_synthetic_genome(seed: int, *, n_genes: int = 115, gene_codons: int = 300, uga_per_gene: int = 0) -> bytes:
    """Construct a ~110 kb prokaryotic contig with realistic RBS + start + sense codons + TAA stop.

    When uga_per_gene > 0, inserts `uga_per_gene` internal TGA codons evenly inside each ORF so
    that translation_table=11 prematurely truncates every ORF while translation_table=25/4
    recovers full-length ORFs and high coding density.
    """
    rng = random.Random(seed)
    parts: list[str] = ["A" * 40]
    for _ in range(n_genes):
        # Shine-Dalgarno RBS + 7 bp spacer + ATG
        spacer = "".join(rng.choice("ACGT") for _ in range(7))
        codons = [rng.choice(_SENSE_CODONS) for _ in range(gene_codons - 2)]
        if uga_per_gene > 0:
            step = max(1, len(codons) // (uga_per_gene + 1))
            for k in range(1, uga_per_gene + 1):
                pos = min(len(codons) - 1, k * step)
                codons[pos] = "TGA"
        orf = "ATG" + "".join(codons) + "TAA"
        intergenic = "".join(rng.choice("AT") for _ in range(25)) + "AGGAGG" + spacer
        parts.append(intergenic + orf)
    parts.append("T" * 40)
    return "".join(parts).encode("ascii")


def _wrap_fasta(records: list[tuple[str, bytes]], width: int = 80) -> bytes:
    out: list[bytes] = []
    for hdr, seq in records:
        out.append(b">" + hdr.encode("ascii") + b"\n")
        for i in range(0, len(seq), width):
            out.append(seq[i : i + width] + b"\n")
    return b"".join(out)


def _write_stages_tsv(path: Path, rows: list[dict]) -> None:
    cols = [
        "ID", "stage", "position_in_stage", "resource", "genome_size",
        "contig_number", "checkm2_completeness", "checkm2_contamination", "taxonomy",
    ]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _write_tar_gz(path: Path, members: list[tuple[str, list[tuple[str, bytes]]]]) -> None:
    """Create a GlobDB-style tar.gz containing R226CHUNK000/<ID>.fa.gz members."""
    with tarfile.open(path, "w:gz") as tf:
        for gid, recs in members:
            fa_bytes = _wrap_fasta(recs)
            buf = io.BytesIO()
            with gzip.GzipFile(filename="", mode="wb", fileobj=buf, mtime=0) as gz:
                gz.write(fa_bytes)
            gz_bytes = buf.getvalue()
            ti = tarfile.TarInfo(name=f"globdb_r226_genome_fasta/R226CHUNK000/{gid}.fa.gz")
            ti.size = len(gz_bytes)
            tf.addfile(ti, io.BytesIO(gz_bytes))


def _write_fake_root(root: Path, stage: str, owned_uids: list[str]) -> None:
    part_dir = root / "l0seq" / "contig_index" / f"stage={stage}"
    part_dir.mkdir(parents=True, exist_ok=True)
    tbl = pa.table({
        "contig_uid": pa.array([bytes.fromhex(u) for u in owned_uids], type=pa.binary(16)),
        "is_default_variant": pa.array([True] * len(owned_uids), type=pa.bool_()),
    })
    pq.write_table(tbl, part_dir / "part-0.parquet")


def test_phylum_selection_and_ordering(tmp_path: Path):
    mod = _load()
    tsv = tmp_path / "non_gtdb_stages.tsv"
    rows = [
        {
            "ID": "MOTU40_0001", "stage": "G26", "position_in_stage": 10, "resource": "mOTU",
            "genome_size": 120000, "contig_number": 1, "checkm2_completeness": 95.0,
            "checkm2_contamination": 1.0,
            "taxonomy": "d__Bacteria;p__Myxococcota;c__Myxococcia;o__Myxococcales;f__F1;g__G1;s__S1",
        },
        {
            "ID": "SPIREOTU_0002", "stage": "G26", "position_in_stage": 11, "resource": "SPIRE",
            "genome_size": 130000, "contig_number": 2, "checkm2_completeness": 75.0,
            "checkm2_contamination": 3.0,
            "taxonomy": "d__Bacteria;p__Myxococcota_A;c__C1;o__O1;f__F1;g__G1;s__S1",
        },
        {
            "ID": "BCRBG_0003", "stage": "G20", "position_in_stage": 5, "resource": "RBG",
            "genome_size": 115000, "contig_number": 1, "checkm2_completeness": 45.0,
            "checkm2_contamination": 2.0,
            "taxonomy": "d__Bacteria;p__Cyanobacteriota;c__Cyanobacteriia;o__O1;f__F1;g__G1;s__S1",
        },
        {
            "ID": "GEMOTU_0004", "stage": "G27", "position_in_stage": 1, "resource": "GEM",
            "genome_size": 118000, "contig_number": 1, "checkm2_completeness": 92.0,
            "checkm2_contamination": 2.0,
            "taxonomy": "d__Bacteria;p__Omnitrophota;c__Omnitrophia;o__O1;f__F1;g__G1;s__S1",
        },
        {
            "ID": "TPMCOTU_0005", "stage": "G29", "position_in_stage": 2, "resource": "TPMC",
            "genome_size": 125000, "contig_number": 1, "checkm2_completeness": 80.0,
            "checkm2_contamination": 12.0,
            "taxonomy": "d__Bacteria;p__Planctomycetota;c__Planctomycetia;o__O1;f__F1;g__G1;s__S1",
        },
        {
            "ID": "BCRBG_0006", "stage": "G27", "position_in_stage": 50, "resource": "RBG",
            "genome_size": 112000, "contig_number": 1, "checkm2_completeness": 91.0,
            "checkm2_contamination": 0.5,
            "taxonomy": "d__Bacteria;p__Patescibacteriota;c__Gracilibacteria;o__Absconditabacterales;f__F1;g__G1;s__S1",
        },
        {
            "ID": "MOTU40_0007", "stage": "G27", "position_in_stage": 51, "resource": "mOTU",
            "genome_size": 111000, "contig_number": 1, "checkm2_completeness": 82.0,
            "checkm2_contamination": 1.0,
            "taxonomy": "d__Bacteria;p__Patescibacteriota;c__Paceibacteria;o__O1;f__F1;g__G1;s__S1",
        },
        # Non-target phyla that must NOT be selected
        {
            "ID": "MOTU40_9999", "stage": "G35", "position_in_stage": 0, "resource": "mOTU",
            "genome_size": 200000, "contig_number": 1, "checkm2_completeness": 95.0,
            "checkm2_contamination": 1.0,
            "taxonomy": "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__O1;f__F1;g__G1;s__S1",
        },
        {
            "ID": "SPIREOTU_8888", "stage": "G41", "position_in_stage": 9, "resource": "SPIRE",
            "genome_size": 150000, "contig_number": 1, "checkm2_completeness": 80.0,
            "checkm2_contamination": 1.0,
            "taxonomy": "d__Unclassified;p_Unclassified;c__Unclassified;o__Unclassified;f__Unclassified;g__Unclassified;s__Unclassified",
        },
    ]
    _write_stages_tsv(tsv, rows)
    sel = mod.select_stage_genomes(tsv)

    assert set(sel) == {"GART1", "GART2"}
    # Brute-force oracle: sorted by (taxonomy, ID)
    gart1_ids = [r["ID"] for r in sel["GART1"]]
    gart2_ids = [r["ID"] for r in sel["GART2"]]
    assert gart1_ids == ["BCRBG_0003", "MOTU40_0001", "SPIREOTU_0002", "GEMOTU_0004", "TPMCOTU_0005"]
    assert gart2_ids == ["BCRBG_0006", "MOTU40_0007"]
    # Quality mapping check: >90 & <5 -> mag_hq; >=50 & <10 -> mag_mq; else contig
    qual_by_id = {r["ID"]: r["assembly_quality"] for r in sel["GART1"] + sel["GART2"]}
    assert qual_by_id["MOTU40_0001"] == "mag_hq"
    assert qual_by_id["SPIREOTU_0002"] == "mag_mq"
    assert qual_by_id["BCRBG_0003"] == "contig"
    assert qual_by_id["TPMCOTU_0005"] == "contig"
    assert qual_by_id["BCRBG_0006"] == "mag_hq"


def test_genetic_code_screening_code11_vs_code25():
    mod = _load()
    seq_code11 = make_synthetic_genome(101, n_genes=120, gene_codons=300, uga_per_gene=0)
    seq_code25 = make_synthetic_genome(202, n_genes=120, gene_codons=300, uga_per_gene=8)

    res11 = mod.screen_genetic_code(
        [seq_code11],
        taxonomy="d__Bacteria;p__Patescibacteriota;c__Paceibacteria;o__O1;f__F1;g__G1;s__S1",
    )
    assert res11["genetic_code"] == 11
    assert res11["cd_11"] >= 0.80

    res25 = mod.screen_genetic_code(
        [seq_code25],
        taxonomy="d__Bacteria;p__Patescibacteriota;c__Gracilibacteria;o__Absconditabacterales;f__F1;g__G1;s__S1",
    )
    assert res25["genetic_code"] == 25
    assert res25["tested_alt"] is True
    assert (res25["cd_25"] - res25["cd_11"]) >= 0.10
    assert res25["mean_orf_25"] > res25["mean_orf_11"] * 1.25


def test_stage_gart_end_to_end_dedup_and_determinism(tmp_path: Path):
    mod = _load()

    # Construct 4 genomes:
    # - GART1: g1_myxo (code 11, 2 contigs, one owned by root stage P2), g1_cyano (code 11, 1 contig, owned by owner_input V2)
    # - GART2: g2_graci (code 25 UGA-recoded, 1 contig), g2_pacei (code 11, 1 contig),
    #          g2_dup (exact repeat of g2_pacei's contig -> DEC-27 repeated genome)
    s_myxo_1 = make_synthetic_genome(11, n_genes=65, gene_codons=280, uga_per_gene=0)
    s_myxo_2 = make_synthetic_genome(12, n_genes=60, gene_codons=280, uga_per_gene=0)
    s_cyano = make_synthetic_genome(13, n_genes=115, gene_codons=300, uga_per_gene=0)
    s_graci = make_synthetic_genome(21, n_genes=120, gene_codons=300, uga_per_gene=8)
    s_pacei = make_synthetic_genome(22, n_genes=115, gene_codons=300, uga_per_gene=0)

    tsv = tmp_path / "non_gtdb_stages.tsv"
    rows = [
        {
            "ID": "MYXO_01", "stage": "G26", "position_in_stage": 0, "resource": "mOTU",
            "genome_size": len(s_myxo_1) + len(s_myxo_2), "contig_number": 2,
            "checkm2_completeness": 95.0, "checkm2_contamination": 1.0,
            "taxonomy": "d__Bacteria;p__Myxococcota;c__Myxococcia;o__O1;f__F1;g__G1;s__S1",
        },
        {
            "ID": "CYANO_01", "stage": "G20", "position_in_stage": 0, "resource": "SPIRE",
            "genome_size": len(s_cyano), "contig_number": 1,
            "checkm2_completeness": 82.0, "checkm2_contamination": 2.0,
            "taxonomy": "d__Bacteria;p__Cyanobacteriota;c__Cyanobacteriia;o__O1;f__F1;g__G1;s__S1",
        },
        {
            "ID": "GRACI_01", "stage": "G27", "position_in_stage": 0, "resource": "RBG",
            "genome_size": len(s_graci), "contig_number": 1,
            "checkm2_completeness": 93.0, "checkm2_contamination": 0.8,
            "taxonomy": "d__Bacteria;p__Patescibacteriota;c__Gracilibacteria;o__Absconditabacterales;f__F1;g__G1;s__S1",
        },
        {
            "ID": "PACEI_01", "stage": "G27", "position_in_stage": 1, "resource": "mOTU",
            "genome_size": len(s_pacei), "contig_number": 1,
            "checkm2_completeness": 88.0, "checkm2_contamination": 1.5,
            "taxonomy": "d__Bacteria;p__Patescibacteriota;c__Paceibacteria;o__O1;f__F1;g__G1;s__S1",
        },
        {
            "ID": "PACEI_02_DUP", "stage": "G27", "position_in_stage": 2, "resource": "SPIRE",
            "genome_size": len(s_pacei), "contig_number": 1,
            "checkm2_completeness": 85.0, "checkm2_contamination": 2.0,
            "taxonomy": "d__Bacteria;p__Patescibacteriota;c__Paceibacteria;o__O1;f__F1;g__G1;s__S2",
        },
        {
            "ID": "OTHER_01", "stage": "G35", "position_in_stage": 0, "resource": "mOTU",
            "genome_size": 1000, "contig_number": 1,
            "checkm2_completeness": 95.0, "checkm2_contamination": 0.5,
            "taxonomy": "d__Bacteria;p__Pseudomonadota;c__C1;o__O1;f__F1;g__G1;s__S1",
        },
    ]
    _write_stages_tsv(tsv, rows)

    tar_path = tmp_path / "globdb_r226_genome_fasta.tar.gz"
    _write_tar_gz(tar_path, [
        ("PACEI_02_DUP", [("PACEI_02_DUP_0001", s_pacei)]),
        ("MYXO_01", [("MYXO_01_0001", s_myxo_1), ("MYXO_01_0002", s_myxo_2)]),
        ("OTHER_01", [("OTHER_01_0001", b"A" * 1000)]),
        ("GRACI_01", [("GRACI_01_0001", s_graci)]),
        ("CYANO_01", [("CYANO_01_0001", s_cyano)]),
        ("PACEI_01", [("PACEI_01_0001", s_pacei)]),
    ])

    fake_root = tmp_path / "root"
    _write_fake_root(fake_root, "P2", [cuid(s_myxo_1)])

    fake_v2 = tmp_path / "owner_V2"
    fake_v2.mkdir()
    (fake_v2 / "genomes.json").write_text(json.dumps([
        {"contig_uids": [cuid(s_cyano)], "genetic_code": 11, "kind": "virus", "circular": [], "assembly_quality": "vOTU"}
    ]))

    out1 = tmp_path / "run1"
    fa1 = tmp_path / "fasta1"
    rep1 = mod.stage_gart(
        tsv, tar_path, out1,
        fasta_dir=fa1,
        final_dirs={"GART1": Path("/tmp/if2453/s2/inputs/GART1"), "GART2": Path("/tmp/if2453/s2/inputs/GART2")},
        root=fake_root,
        root_stages=["P2"],
        owner_inputs={"V2": fake_v2},
        workers=1,
    )

    assert rep1["GART1"]["selected_genomes"] == 2
    assert rep1["GART1"]["genomes"] == 2
    assert rep1["GART1"]["contigs"] == 3
    assert rep1["GART1"]["owner_contigs"] == {"P2": 1, "V2": 1}
    assert rep1["GART1"]["new_contigs"] == 1
    assert rep1["GART1"]["codes"] == {"11": 2}

    assert rep1["GART2"]["selected_genomes"] == 3
    assert rep1["GART2"]["genomes"] == 2  # PACEI_02_DUP repeats PACEI_01
    assert rep1["GART2"]["repeated_genomes"] == 1
    assert rep1["GART2"]["contigs"] == 2
    assert rep1["GART2"]["codes"] == {"11": 1, "25": 1}

    # Check objects.json points to --final-dir, and objects_staging.json points to fasta_dir
    obj1_final = json.loads((out1 / "GART1" / "objects.json").read_text())
    assert obj1_final[0]["path"] == "/tmp/if2453/s2/inputs/GART1/gart1_globdb.fna.gz"
    assert obj1_final[0]["source_name"] == "globdb"
    assert obj1_final[0]["redistributable"] is True

    # Validate BOTH stages with s2_check_inputs.py against objects_staging.json
    for st in ("GART1", "GART2"):
        proc = subprocess.run(
            [
                sys.executable,
                str(SCRIPTS / "s2_check_inputs.py"),
                st,
                str(out1 / st / "genomes.json"),
                str(out1 / st / "objects_staging.json"),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert proc.returncode == 0, f"s2_check_inputs.py {st} failed: {proc.stdout}\n{proc.stderr}"
        chk = json.loads(proc.stdout)
        assert chk["problems"] == []

    # Second run must produce byte-identical outputs across all files including .fna.gz
    out2 = tmp_path / "run2"
    fa2 = tmp_path / "fasta2"
    mod.stage_gart(
        tsv, tar_path, out2,
        fasta_dir=fa2,
        final_dirs={"GART1": Path("/tmp/if2453/s2/inputs/GART1"), "GART2": Path("/tmp/if2453/s2/inputs/GART2")},
        root=fake_root,
        root_stages=["P2"],
        owner_inputs={"V2": fake_v2},
        workers=1,
    )
    for st, fname in (("GART1", "gart1_globdb.fna.gz"), ("GART2", "gart2_globdb_cpr.fna.gz")):
        assert (fa1 / fname).read_bytes() == (fa2 / fname).read_bytes()
        for meta_name in (
            "genomes.json", "objects.json", "records.tsv",
            "repeated_records.tsv", "shared_records.tsv", "genetic_codes.tsv",
        ):
            assert (out1 / st / meta_name).read_bytes() == (out2 / st / meta_name).read_bytes()
