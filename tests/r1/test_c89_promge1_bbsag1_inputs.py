"""Card C-89 (DEC-65, DEC-80, DEC-84 (2), DEC-88): stage tower-landed MGE/plasmid (`PROMGE1`) and microfluidic
single-cell (`BBSAG1`) cohorts for `r26.11.6`.

Independent oracles in this file:
- `contig_uid` and `genome_uid` computed directly via `blake3`;
- synthetic coding sequences generated from a sense-codon pool with planted internal `TGA` codons
  (`code 4` Trp vs `code 25` Gly) and standard `code 11` controls;
- hand-constructed `.tar.gz`, `.tar.xz`, `.zip`, `.fna.gz`, Prodigal `-d` `.ffn` single-gene FASTAs,
  macOS `__MACOSX` resource forks, duplicate-entry `.xlsx` zips, and Parquet `r232_index` tables;
- every test runs under both normal `PATH` and tools-hidden `PATH=/usr/bin:/bin` (skipping only live `skani` when absent).
"""

from __future__ import annotations

import bz2
import csv
import gzip
import importlib.util
import io
import json
import lzma
import random
import shutil
import subprocess
import sys
import tarfile
import warnings
import zipfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("C-89")

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
TEMPLATES = BSP / "stages" / "templates"
NOTE = "planning/measurements/2026-11-27_c89_promge1_bbsag1_inputs.md"


def _load(name: str, file: str):
    sys.path.insert(0, str(BSP))
    spec = importlib.util.spec_from_file_location(name, BSP / file)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(BSP / file)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _mod():
    return _load("c89_promge1_bbsag1_inputs", "c89_promge1_bbsag1_inputs.py")


@pytest.fixture(params=["normal-path", "tools-hidden"], autouse=True)
def path_mode(request, monkeypatch):
    if request.param == "tools-hidden":
        monkeypatch.setenv("PATH", "/usr/bin:/bin")
    return request.param


def cuid(seq: str | bytes) -> str:
    b = seq.encode("ascii") if isinstance(seq, str) else seq
    return blake3(b.upper()).digest(length=16).hex()


def guid(seqs: list[str | bytes]) -> str:
    raws = sorted(bytes.fromhex(cuid(s)) for s in dict.fromkeys(seqs))
    return blake3(b"".join(raws)).digest(length=16).hex()


def _fa(records: list[tuple[str, str]], width: int = 70) -> str:
    out = io.StringIO()
    for hdr, s in records:
        out.write(f">{hdr}\n")
        for i in range(0, len(s), width):
            out.write(s[i : i + width] + "\n")
    return out.getvalue()


def _gz_bytes(text: str) -> bytes:
    return gzip.compress(text.encode("utf-8"), mtime=0)


def _rows(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ------------------------------------------------------------------------------------------------------------------
# Independent synthetic coding-sequence oracle (for pyrodigal single-mode genetic-code screen 11 / 4 / 25)
# ------------------------------------------------------------------------------------------------------------------

BASE_SENSE_CODONS = [
    "GCT", "GCC", "GCA", "CGT", "CGC", "AAT", "AAC", "GAT", "GAC", "TGT",
    "CAA", "CAG", "GAA", "GAG", "GGT", "GGC", "GGA", "CAT", "CAC", "ATT",
    "ATC", "CTT", "CTC", "TTA", "TTG", "AAA", "AAG", "TTT", "TTC", "CCT",
    "CCC", "TCT", "TCC", "ACT", "ACC", "TAT", "TAC", "GTT", "GTC", "GTA",
]


def synth_coding_genome(
    seed: int,
    *,
    n_genes: int = 24,
    codons_per_gene: int = 330,
    tga_every: int = 0,
    tgg_every: int = 85,
) -> str:
    rng = random.Random(seed)
    spacer = "TAATTTTAAAGGAGATATACAT"
    chunks: list[str] = ["AATTCCGGTTAA"]
    for _ in range(n_genes):
        orf = ["ATG"]
        for k in range(1, codons_per_gene):
            if tga_every and k % tga_every == 0:
                orf.append("TGA")
            elif tgg_every and k % tgg_every == 0:
                orf.append("TGG")
            else:
                orf.append(rng.choice(BASE_SENSE_CODONS))
        orf.append("TAA")
        chunks.append("".join(orf))
        chunks.append(spacer)
    return "".join(chunks)


def _rc(seq: str) -> str:
    tr = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(tr)[::-1]


def _write_r232_index(d: Path) -> None:
    d.mkdir(parents=True, exist_ok=True)
    gt = [
        ("GCF_058350605", "GCA_058350605", True, "GCF_058350605", "s__Escherichia coli"),
        ("GCA_040223915", "GCA_040223915", False, "GCA_040223911", "s__Methanoperedens sp040223915"),
    ]
    pq.write_table(
        pa.table({
            "accession": [g[0] for g in gt],
            "genbank_accession": [g[1] for g in gt],
            "is_representative": [g[2] for g in gt],
            "representative": [g[3] for g in gt],
            "species": [g[4] for g in gt],
        }),
        d / "gtdb_r232.parquet",
    )
    pq.write_table(
        pa.table({
            "globdb_id": ["GCA_040223915"],
            "insdc_accession": ["GCA_040223915"],
            "original_id": [None],
            "cluster_id": [None],
        }),
        d / "globdb_r232.parquet",
    )


# ==================================================================================================================
# 1. Magic-byte format detection & archive member census
# ==================================================================================================================


def test_magic_bytes_and_archive_census_flags_prodigal_ffn_proteins_and_macosx(tmp_path: Path):
    mod = _mod()
    gz_fna = tmp_path / "mge.fna.gz"
    gz_fna.write_bytes(_gz_bytes(_fa([("c1 circular", "ACGTACGT" * 25)])))
    bz_file = tmp_path / "sample.bz2"
    bz_file.write_bytes(bz2.compress(b"hello"))
    xz_tar = tmp_path / "sample.tar.xz"
    tbuf = io.BytesIO()
    with tarfile.open(fileobj=tbuf, mode="w") as tf:
        b = _fa([("m1", "ACGT" * 60)]).encode()
        ti = tarfile.TarInfo("bins/m1.fa")
        ti.size = len(b)
        tf.addfile(ti, io.BytesIO(b))
    xz_tar.write_bytes(lzma.compress(tbuf.getvalue()))

    # Prodigal -d single-gene nucleotide dump (like MFD01716.bin.1.47.fasta)
    ffn_fake = tmp_path / "MFD01716.bin.1.47.fasta"
    ffn_fake.write_text(
        ">contig_14181_1 # 3952 # 4143 # -1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=AGxAG\n"
        "ATGTATTACTATTACATCTCACTTATACATGTATAA\n"
        ">contig_14181_2 # 6199 # 6438 # 1 # ID=1_2;partial=00;start_type=ATG;rbs_motif=None\n"
        "ATGGGAGATGCAGTATCAACAACAAGTAAATAGTATTAA\n"
    )

    zpath = tmp_path / "sat_and_sag.zip"
    with zipfile.ZipFile(zpath, "w") as zf:
        zf.writestr("1.P4P2/", "")
        zf.writestr("__MACOSX/1.P4P2/._foo.fna", b"\x00\x05\x16\x07mac")
        zf.writestr("1.P4P2/.DS_Store", b"\x00\x00\x00\x01Bud1")
        zf.writestr("1.P4P2/P4/Proteins/sat1.prt", ">p1\nMEFILPQZ*\n")
        zf.writestr("1.P4P2/P4/Replicons/sat1.fna", _fa([("sat1", "ACGTURYKMSWBDHVN" * 12)]))
        zf.writestr("aln.fasta", ">a1\nACGT--ACGT\n")
        zf.writestr("bad.fna", ">b1\nACGT1234ACGT\n")
        # Duplicate zip entry name (as in BLM1-MAGs_SAGs_2011_2021.zip)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            for _ in range(2):
                zf.writestr("information regarding the MAGs and SAGs of BLM1.xlsx", b"PK_xlsx_stub")

    assert mod.magic_format(gz_fna) == "gzip"
    assert mod.magic_format(bz_file) == "bzip2"
    assert mod.magic_format(xz_tar) == "xz"
    assert mod.magic_format(zpath) == "zip"
    assert mod.magic_format(ffn_fake) == "fasta"

    rep = mod.census([gz_fna, xz_tar, ffn_fake, zpath])
    assert rep[str(gz_fna)]["dna_records"] == 1
    assert rep[str(xz_tar)]["dna_records"] == 1 and rep[str(xz_tar)]["dna_bases"] == 240
    assert rep[str(ffn_fake)]["members"][0]["alphabet"] == "prodigal_cds_ffn"

    zmem = [(m["name"], m["alphabet"]) for m in rep[str(zpath)]["members"]]
    assert ("1.P4P2/", "directory") in zmem
    assert ("__MACOSX/1.P4P2/._foo.fna", "macosx") in zmem
    assert ("1.P4P2/.DS_Store", "macosx") in zmem
    assert ("1.P4P2/P4/Proteins/sat1.prt", "protein") in zmem
    assert ("1.P4P2/P4/Replicons/sat1.fna", "dna") in zmem
    assert ("aln.fasta", "alignment") in zmem
    assert ("bad.fna", "non_iupac") in zmem
    assert sum(1 for n, a in zmem if n.endswith(".xlsx") and a == "not_fasta") == 2


# ==================================================================================================================
# 2. Genetic code screen (11, 4, 25) & 59-stage prefix definition
# ==================================================================================================================


def test_genetic_code_screens_11_4_25_and_59_stage_prefix_defaults():
    mod = _mod()
    # 13 root + 19 r26.11.4 + 24 r26.11.5 + 3 C-86 stages (O4, GV2, HROM1) = 59 prior stages
    assert len(mod.ROOT_STAGES) == 13
    assert len(mod.R26_11_4_STAGES) == 19
    assert len(mod.R26_11_5_STAGES) == 24
    assert mod.C86_STAGES == ("O4", "GV2", "HROM1")
    assert len(set(mod.ROOT_STAGES) | set(mod.R26_11_4_STAGES) | set(mod.R26_11_5_STAGES) | set(mod.C86_STAGES)) == 59

    seq_11 = synth_coding_genome(411, tga_every=0, tgg_every=80).encode()
    seq_c4 = synth_coding_genome(404, tga_every=75, tgg_every=0).encode()
    seq_c25 = synth_coding_genome(425, tga_every=16, tgg_every=90).encode()

    assert mod.screen_prokaryote_code([seq_11], "d__Bacteria;p__Bacteroidota;g__Prevotella")["genetic_code"] == 11
    assert mod.screen_prokaryote_code([seq_c4], "d__Bacteria;p__Firmicutes_A;g__Unclassified")["genetic_code"] == 4
    assert mod.screen_prokaryote_code([seq_c25], "d__Bacteria;p__Patescibacteriota;c__Saccharimonadia")["genetic_code"] == 25
    # Degenerate / high-N assembly (>5 kb of Ns) and fragmented multi-contig SAG (>500 kb with short chaff) must not crash
    assert mod.screen_prokaryote_code([b"N" * 12_000], "")["genetic_code"] == 11
    chaff = [b"ACGT" * 15 for _ in range(400)]
    assert mod.screen_prokaryote_code([seq_c25 * 22, *chaff], "d__Bacteria;p__Patescibacteriota")["genetic_code"] == 25


# ==================================================================================================================
# 3. Stage PROMGE1 (8 MGE / plasmid / Borg / satellite subcohorts)
# ==================================================================================================================


def _build_promge1_fixture(tmp_path: Path) -> dict:
    dl = tmp_path / "dl_promge1"
    for sub in ("PROMGE2", "GSPR", "MFD_PLAS", "ECE_HODOR", "SAT_P4P2", "ECE_METH",
                "PICJI_SAT/assemblies", "BORG_ONT/assemblies"):
        (dl / sub).mkdir(parents=True)

    s_pm_is = synth_coding_genome(501, n_genes=4, codons_per_gene=140)
    s_pm_phg = synth_coding_genome(502, n_genes=5, codons_per_gene=150)
    s_pm_cell = synth_coding_genome(503, n_genes=3, codons_per_gene=120)
    s_pm_rescued = synth_coding_genome(519, n_genes=4, codons_per_gene=135)
    s_pm_nonrep = synth_coding_genome(520, n_genes=4, codons_per_gene=125)
    s_pm_nested = synth_coding_genome(521, n_genes=5, codons_per_gene=145)
    s_gspr_core = synth_coding_genome(504, n_genes=5, codons_per_gene=160)
    s_gspr_dtr = s_gspr_core + s_gspr_core[:42]  # 42 bp exact DTR -> trimmed & circular (DEC-65 (1))
    s_mfd_circ = synth_coding_genome(505, n_genes=5, codons_per_gene=150)
    tir_arm = "ACGTTGCATGCAACGTTGCAACGTTGCATGCAACGT"  # 36 bp TIR
    s_mfd_tir = tir_arm + synth_coding_genome(506, n_genes=4, codons_per_gene=140) + _rc(tir_arm)
    s_hodor_38 = synth_coding_genome(507, n_genes=6, codons_per_gene=160)
    s_hodor_mag1 = synth_coding_genome(508, n_genes=4, codons_per_gene=140)
    s_hodor_mag2 = synth_coding_genome(509, n_genes=4, codons_per_gene=140)
    s_hodor_pro = synth_coding_genome(510, n_genes=4, codons_per_gene=140)
    s_sat_p4 = synth_coding_genome(511, n_genes=4, codons_per_gene=140)
    s_sat_p2 = synth_coding_genome(512, n_genes=5, codons_per_gene=150)
    s_meth_ece = synth_coding_genome(513, n_genes=5, codons_per_gene=150)
    s_meth_host = synth_coding_genome(514, n_genes=6, codons_per_gene=160)
    s_picji_chr = synth_coding_genome(515, n_genes=6, codons_per_gene=160)
    s_picji_pl = synth_coding_genome(516, n_genes=4, codons_per_gene=130)
    s_borg1 = synth_coding_genome(517, n_genes=5, codons_per_gene=150)
    s_borg2 = synth_coding_genome(518, n_genes=5, codons_per_gene=150)

    # 1. PROMGE2 (includes pg3 specI rep, SPIRE MAG, cluster-rescued pg3 non-rep, pure pg3 non-rep, nested, and cellular)
    (dl / "PROMGE2" / "mge100.fna.gz").write_bytes(_gz_bytes(_fa([
        ("MGE_GCA_014491485.1_1.SAMN1.CTG1:100-2000 mge=is_tn:1;genome_type=ACC;mge_type=non-nested;size=1901;n_genes=4;mgeR=rve:1;mge100=proMGE100_2.0_0000001_is_tn_1_1", s_pm_is),
        ("MGE_SAMEA5618361.psa_megahit.psb_metabat2.00024_k141_95:500-9000 mge=phage:1;genome_type=ACC;mge_type=non-nested;size=8501;n_genes=5;mgeR=int_p2:1;mge100=proMGE100_2.0_0000002_phage_1_1", s_pm_phg),
        ("MGE_GCA_088888885.1_1.SAMN2.CTG1:50-3000 mge=ce:1;genome_type=ACC;mge_type=non-nested;size=2951;n_genes=4;mgeR=int_sxt:1;mge100=proMGE100_2.0_0000004_ce_2_1", s_pm_rescued),
        ("MGE_GCA_099999995.1_1.SAMN3.CTG1:50-2500 mge=is_tn:1;genome_type=ACC;mge_type=non-nested;size=2451;n_genes=4;mgeR=rve:1;mge100=proMGE100_2.0_0000005_is_tn_1_1", s_pm_nonrep),
        ("MGE_GCA_014491485.1_1.SAMN1.CTG4:10-4000 mge=phage:1,is_tn:1;genome_type=ACC;mge_type=nested;size=3991;n_genes=5;mgeR=int_p2:1;mge100=proMGE100_2.0_0000006_nested_1_1", s_pm_nested),
        ("MGE_GCA_014491485.1_1.SAMN1.CTG3:10-1500 mge=cellular:1;genome_type=COR;mge_type=non-nested;size=1491;n_genes=3;mgeR=xer:1;mge100=proMGE100_2.0_0000003_cellular_1_1", s_pm_cell),
    ])))
    (dl / "PROMGE2" / "genomes_taxonomy.tsv.gz").write_bytes(_gz_bytes(
        "genome\tsource\tspeci\tkingdom\tphylum\tclass\ttax_order\tfamily\tgenus\tspecies\trepresentative\n"
        "GCA_014491485.1\tpg3\tspecI_v4_00001\tBacteria\tProteobacteria\tGamma\tEntero\tEntero\tSalmonella\tS. enterica\tTrue\n"
        "GCA_088888885.1\tpg3\tspecI_v4_00001\tBacteria\tProteobacteria\tGamma\tEntero\tEntero\tSalmonella\tS. enterica\tFalse\n"
        "GCA_099999995.1\tpg3\tspecI_v4_00001\tBacteria\tProteobacteria\tGamma\tEntero\tEntero\tSalmonella\tS. enterica\tFalse\n"
        "spire_mag_00000001\tspire\tspecI_v4_00002\tBacteria\tBacteroidota\tBacteroidia\tBacteroidales\tPrevotellaceae\tPrevotella\tPrevotella sp.\tFalse\n"
    ))
    (dl / "PROMGE2" / "mge100_clusters.tsv.gz").write_bytes(_gz_bytes(
        "mge100_cluster_id\tmember_id\trepresentative_id\tn_members\tn_species\n"
        "proMGE100_1.0_0000001_is_tn_1_1\tMGE_GCA_014491485.1_1.SAMN1.CTG1:100-2000\tMGE_GCA_014491485.1_1.SAMN1.CTG1:100-2000\t1\t1\n"
        "proMGE100_1.0_0000002_phage_1_1\tMGE_SAMEA5618361.psa_megahit.psb_metabat2.00024_k141_95:500-9000\tMGE_SAMEA5618361.psa_megahit.psb_metabat2.00024_k141_95:500-9000\t1\t1\n"
        "proMGE100_1.0_0000004_ce_2_1\tMGE_GCA_088888885.1_1.SAMN2.CTG1:50-3000\tMGE_GCA_088888885.1_1.SAMN2.CTG1:50-3000\t2\t1\n"
        "proMGE100_1.0_0000004_ce_2_1\tMGE_SAMEA5618361.psa_megahit.psb_metabat2.00024_k141_99:50-3000\tMGE_GCA_088888885.1_1.SAMN2.CTG1:50-3000\t2\t1\n"
        "proMGE100_1.0_0000005_is_tn_1_1\tMGE_GCA_099999995.1_1.SAMN3.CTG1:50-2500\tMGE_GCA_099999995.1_1.SAMN3.CTG1:50-2500\t1\t1\n"
        "proMGE100_1.0_0000003_cellular_1_1\tMGE_GCA_014491485.1_1.SAMN1.CTG3:10-1500\tMGE_GCA_014491485.1_1.SAMN1.CTG3:10-1500\t5\t1\n"
    ))

    # 2. GSPR (includes 1 exact repeat of s_mfd_circ to test DEC-27)
    (dl / "GSPR" / "all_soil_plasmids_derep.fna.gz").write_bytes(_gz_bytes(_fa([
        ("2088090014|GPIPI_16944593", s_gspr_dtr),
        ("2088090015|GPICI_9089551", s_mfd_circ),
    ])))
    (dl / "GSPR" / "accessions.txt").write_text("3300012987\tIMG\n")

    # 3. MFD_PLAS: mfd_plasmids_drep_and_members.fasta + pfd_metadata.tsv + Prodigal .ffn MFD01716.bin.1.47.fasta (must skip!)
    (dl / "MFD_PLAS" / "mfd_plasmids_drep_and_members.fasta").write_text(_fa([
        ("contig_100057__MFD09222", s_mfd_circ),
        ("contig_200001__MFD01111", s_mfd_tir),
    ]))
    (dl / "MFD_PLAS" / "pfd_metadata.tsv").write_text(
        "representative\tmember_id\tptu\tdptu\trps\tdb\tsource\ttechnology\ttopology\tlength\n"
        f"contig_100057__MFD09222\tcontig_100057__MFD09222\tptu1\tdPTU_1\tNA\tMFD\tMFD\tONT\tCircular\t{len(s_mfd_circ)}\n"
        f"contig_200001__MFD01111\tcontig_200001__MFD01111\tptu2\tdPTU_2\tNA\tMFD\tMFD\tONT\tInverted terminal repeat\t{len(s_mfd_tir)}\n"
    )
    (dl / "MFD_PLAS" / "MFD01716.bin.1.47.fasta").write_text(
        ">contig_14181_1 # 3952 # 4143 # -1 # ID=1_1;partial=00;start_type=ATG\nATGTATTACTATTACATAA\n"
    )

    # 4. ECE_HODOR
    (dl / "ECE_HODOR" / "Hodor.38.genomes.fasta").write_text(_fa([
        ("ERR11473975_NODE_51_curated_circular", s_hodor_38),
    ]))
    (dl / "ECE_HODOR" / "Hodor.prophages.fasta").write_text(_fa([
        ("huge_phage_1554_curated_complete_prophage", s_hodor_pro),
    ]))
    (dl / "ECE_HODOR" / "prophage01.all.pure.fa").write_text(_fa([
        ("GENOME083466_1_1", s_hodor_pro),  # exact repeat of Hodor.prophages record -> DEC-27 collapse
    ]))
    with tarfile.open(dl / "ECE_HODOR" / "draft.Hodor.genomes.tar.gz", "w:gz") as tf:
        b = _fa([("d1", s_hodor_mag1), ("d2", s_hodor_mag2)]).encode()
        ti = tarfile.TarInfo("2.all.pebblescout/DRR127535.hodor.MAG.fa")
        ti.size = len(b)
        tf.addfile(ti, io.BytesIO(b))

    # 5. SAT_P4P2: Replicons/*.fna staged, Proteins/*.prt and __MACOSX skipped
    with zipfile.ZipFile(dl / "SAT_P4P2" / "P4P2_DereplicatedGenomes.zip", "w") as zf:
        zf.writestr("__MACOSX/1.P4P2_DereplicatedGenomes/P4/Replicons/._p4.fna", b"mac")
        zf.writestr("1.P4P2_DereplicatedGenomes/P4/Proteins/Contig.1.P4.prt", ">p1\nMEFILPQ*\n")
        zf.writestr("1.P4P2_DereplicatedGenomes/P4/Replicons/Contig.1.P4.fna", _fa([("c1_p4", s_sat_p4)]))
        zf.writestr("1.P4P2_DereplicatedGenomes/P2/Replicons/ESCO.Prophage.1.P2.fna", _fa([("c1_p2", s_sat_p2)]))

    # 6. ECE_METH: circular HMp-Plasmid-like ECE + Methanoperedens host chromosome (host chromosome left out or flagged)
    (dl / "ECE_METH" / "Mp_and_ECE_genomes.fa.gz").write_bytes(_gz_bytes(_fa([
        ('FINAL_77kb_HMp-Plasmid-like_SR-VP_PACBIO-HIFI_complete id=1 bin="FINAL_77kb"', s_meth_ece),
        ('FINAL_2-76MB_Methanoperedens_SR-VP_PACBIO-HIFI_complete id=2 bin="FINAL_2-76MB"', s_meth_host),
    ])))

    # 7. PICJI_SAT (INSDC complete genome + plasmid; DEC-65 (2) byte-identical, redistributable: true)
    (dl / "PICJI_SAT" / "MANIFEST_PICJI_SAT_PRJNA1424553.tsv").write_text(
        "stage\tbioproject\taccession\torganism_name\tbase_count\tassembly_level\turl\n"
        f"PICJI_SAT\tPRJNA1424553\tGCF_058350605.1\tEscherichia coli\t{len(s_picji_chr) + len(s_picji_pl)}\tComplete Genome\thttps://ftp.ncbi.nlm.nih.gov/x\n"
    )
    (dl / "PICJI_SAT" / "assemblies" / "GCF_058350605.1_genomic.fna.gz").write_bytes(_gz_bytes(_fa([
        ("CM177969.1 Escherichia coli strain M24FEP-TDT chromosome, complete sequence, whole genome shotgun sequence", s_picji_chr),
        ("JBVKFH010000002.1 Escherichia coli strain M24FEP-TDT plasmid pB25FEP-TDT, whole genome shotgun sequence", s_picji_pl),
    ])))

    # 8. BORG_ONT (INSDC Nanopore Borg assembly; linear, DEC-65 (2) byte-identical, redistributable: true)
    (dl / "BORG_ONT" / "MANIFEST_BORG_ONT_PRJNA1119519.tsv").write_text(
        "stage\tbioproject\taccession\torganism_name\tbase_count\tassembly_level\turl\n"
        f"BORG_ONT\tPRJNA1119519\tGCA_040223915.1\tCandidatus Methanoperedens sp.\t{len(s_borg1) + len(s_borg2)}\tContig\thttps://ftp.ncbi.nlm.nih.gov/y\n"
    )
    (dl / "BORG_ONT" / "assemblies" / "GCA_040223915.1_genomic.fna.gz").write_bytes(_gz_bytes(_fa([
        ("JBEHDZ010000001.1 MAG: Candidatus Methanoperedens sp. isolate bMp part_1, whole genome shotgun sequence", s_borg1),
        ("JBEHDZ010000002.1 MAG: Candidatus Methanoperedens sp. isolate bMp part_2, whole genome shotgun sequence", s_borg2),
    ])))

    idx = tmp_path / "index_p"
    _write_r232_index(idx)

    owner_pl2 = tmp_path / "owners_p" / "PL2"
    owner_pl2.mkdir(parents=True)
    (owner_pl2 / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(s_mfd_circ)]}]))
    owner_hrom1 = tmp_path / "owners_p" / "HROM1"
    owner_hrom1.mkdir(parents=True)
    (owner_hrom1 / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(s_pm_is)]}]))

    return {
        "dirs": {
            "PROMGE2": [dl / "PROMGE2"],
            "GSPR": [dl / "GSPR"],
            "MFD_PLAS": [dl / "MFD_PLAS"],
            "ECE_HODOR": [dl / "ECE_HODOR"],
            "SAT_P4P2": [dl / "SAT_P4P2"],
            "ECE_METH": [dl / "ECE_METH"],
            "PICJI_SAT": [dl / "PICJI_SAT"],
            "BORG_ONT": [dl / "BORG_ONT"],
        },
        "index_dir": idx,
        "owners": {"PL2": owner_pl2, "HROM1": owner_hrom1},
        "seqs": {
            "pm_is": s_pm_is, "pm_phg": s_pm_phg, "pm_cell": s_pm_cell,
            "pm_rescued": s_pm_rescued, "pm_nonrep": s_pm_nonrep, "pm_nested": s_pm_nested,
            "gspr_core": s_gspr_core, "mfd_circ": s_mfd_circ, "mfd_tir": s_mfd_tir,
            "hodor_38": s_hodor_38, "hodor_mag1": s_hodor_mag1, "hodor_mag2": s_hodor_mag2,
            "hodor_pro": s_hodor_pro, "sat_p4": s_sat_p4, "sat_p2": s_sat_p2,
            "meth_ece": s_meth_ece, "meth_host": s_meth_host,
            "picji_chr": s_picji_chr, "picji_pl": s_picji_pl,
            "borg1": s_borg1, "borg2": s_borg2,
        },
    }


def test_promge1_stage_eight_subcohorts_dec65_dtr_tir_repeat_collapse_and_s2_check(tmp_path: Path):
    mod = _mod()
    fx = _build_promge1_fixture(tmp_path)
    rep1 = mod.stage_promge1(
        fx["dirs"], fx["index_dir"], tmp_path / "meta_a", tmp_path / "fasta_a",
        root=None, root_stages=(), owner_inputs=fx["owners"], later_inputs={},
        inputs_home=Path("/tmp/if2453/s2/inputs"), work=tmp_path / "work_a",
        exclude_promge2_cellular=True, promge2_non_nested_only=True, promge2_spire_and_pg3_reps_only=True,
    )
    rep2 = mod.stage_promge1(
        fx["dirs"], fx["index_dir"], tmp_path / "meta_b", tmp_path / "fasta_b",
        root=None, root_stages=(), owner_inputs=fx["owners"], later_inputs={},
        inputs_home=Path("/tmp/if2453/s2/inputs"), work=tmp_path / "work_b",
        exclude_promge2_cellular=True, promge2_non_nested_only=True, promge2_spire_and_pg3_reps_only=True,
    )

    meta = tmp_path / "meta_a" / "PROMGE1"
    # Under 64 MiB: no small files written onto /stor (fasta_a/PROMGE1)
    fa_dir = tmp_path / "fasta_a" / "PROMGE1"
    assert not fa_dir.exists() or list(fa_dir.iterdir()) == []

    genomes = json.loads((meta / "genomes.json").read_text())
    assert {g["kind"] for g in genomes} == {"plasmid"}
    assert {g["genetic_code"] for g in genomes} == {11}

    s = fx["seqs"]
    by_uids = {tuple(g["contig_uids"]): g for g in genomes}

    # GSPR 42 bp exact DTR trimmed once & marked circular (DEC-65 (1))
    assert (cuid(s["gspr_core"]),) in by_uids
    assert by_uids[(cuid(s["gspr_core"]),)]["circular"] == [cuid(s["gspr_core"])]
    trims = _rows(meta / "trims.tsv")
    assert any(int(t["bases_removed"]) == 42 for t in trims)

    # MFD_PLAS: Circular plasmid is circular; Inverted terminal repeat plasmid stays linear & untrimmed!
    assert by_uids[(cuid(s["mfd_circ"]),)]["circular"] == [cuid(s["mfd_circ"])]
    assert (cuid(s["mfd_tir"]),) in by_uids
    assert by_uids[(cuid(s["mfd_tir"]),)]["circular"] == []

    # BORG_ONT: multi-contig Borg MAG stays linear & untrimmed under DEC-65 (2)
    borg_key = (cuid(s["borg1"]), cuid(s["borg2"]))
    assert borg_key in by_uids and by_uids[borg_key]["circular"] == []

    # Option 3b PROMGE2 filters: pg3 rep + SPIRE MAG + cluster-rescued ce kept; cellular, nested, and pure pg3 non-rep excluded
    all_cus = {u for g in genomes for u in g["contig_uids"]}
    assert {cuid(s["pm_is"]), cuid(s["pm_phg"]), cuid(s["pm_rescued"])} <= all_cus
    assert cuid(s["pm_cell"]) not in all_cus
    assert cuid(s["pm_nested"]) not in all_cus
    assert cuid(s["pm_nonrep"]) not in all_cus
    assert rep1["promge2_rescued_clusters"] == 1
    assert rep1["promge2_excluded_records"] == {
        "promge2_cellular_recombinase_non_mge": 1,
        "promge2_nested_composite": 1,
        "promge2_pg3_nonrep_host": 1,
    }
    left = _rows(meta / "left_out.tsv")
    left_reasons = {r["reason"] for r in left}
    assert {
        "promge2_cellular_recombinase_non_mge",
        "promge2_nested_composite",
        "promge2_pg3_nonrep_host",
        "host_archaeal_chromosome_not_ece",
    } <= left_reasons

    # Prodigal .ffn MFD01716.bin.1.47.fasta and SAT_P4P2 Proteins/*.prt skipped in skipped_members.tsv
    skipped = _rows(meta / "skipped_members.tsv")
    skipped_names = {r["member"] for r in skipped}
    assert any("MFD01716.bin.1.47.fasta" in n for n in skipped_names)
    assert any("Proteins" in n for n in skipped_names)

    # DEC-27 exact repeats recorded (GSPR s_mfd_circ duplicate + prophage01 s_hodor_pro duplicate)
    assert rep1["repeated_genomes"] >= 2

    # Ownership against PL2 and C-86's HROM1
    assert rep1["held_by"] == {"HROM1": 1, "PL2": 1}

    # Redistributable provenance split: PROMGE2 is false (DEC-80 (3)), open CC-BY / INSDC is true
    objs = {Path(o["path"]).name: o for o in json.loads((meta / "objects.json").read_text())}
    assert objs["promge1_promge2.fna.gz"]["redistributable"] is False
    assert objs["promge1_plasmids_eces.fna.gz"]["redistributable"] is True

    # Byte-for-byte determinism
    for f in ("genomes.json", "objects.json", "records.tsv", "left_out.tsv", "skipped_members.tsv",
              "repeated_records.tsv", "shared_records.tsv", "trims.tsv", "r232_status.tsv",
              "promge1_promge2.fna.gz", "promge1_plasmids_eces.fna.gz"):
        assert (meta / f).read_bytes() == (tmp_path / "meta_b" / "PROMGE1" / f).read_bytes(), f

    chk = subprocess.run(
        [sys.executable, str(BSP / "s2_check_inputs.py"), "PROMGE1",
         str(meta / "genomes.json"), str(meta / "objects_staging.json")],
        capture_output=True, text=True, check=False,
    )
    assert chk.returncode == 0, chk.stdout + chk.stderr
    assert json.loads(chk.stdout)["problems"] == []


# ==================================================================================================================
# 4. Stage BBSAG1 (BBSAG20 MQ/HQ SAGs+MAGs + HAIMA_LR_SAG + DEATHVALLEY_SAG, codes 11/4/25, split <= 40 Gbp)
# ==================================================================================================================


def _build_bbsag1_fixture(tmp_path: Path) -> dict:
    dl = tmp_path / "dl_bbsag1"
    for sub in ("BBSAG20", "HAIMA_LR_SAG", "DEATHVALLEY_SAG"):
        (dl / sub).mkdir(parents=True)

    s_sag_hq11 = synth_coding_genome(601, n_genes=24, codons_per_gene=320)
    s_sag_mq25 = synth_coding_genome(602, n_genes=24, codons_per_gene=320, tga_every=16, tgg_every=90)
    s_sag_lq = synth_coding_genome(603, n_genes=12, codons_per_gene=200)
    s_mag_hq4 = synth_coding_genome(604, n_genes=24, codons_per_gene=320, tga_every=70, tgg_every=0)
    s_mag_lq = synth_coding_genome(605, n_genes=12, codons_per_gene=200)
    s_haima_sag = synth_coding_genome(606, n_genes=24, codons_per_gene=320)
    s_haima_cbin_core = synth_coding_genome(607, n_genes=24, codons_per_gene=320)
    s_haima_cbin_dtr = s_haima_cbin_core + s_haima_cbin_core[:50]  # 50 bp exact DTR on cbin cMAG
    s_haima_sr = synth_coding_genome(608, n_genes=24, codons_per_gene=320)
    s_dv_sag = synth_coding_genome(609, n_genes=24, codons_per_gene=320)

    # BBSAG20 summary.tsv: note QLF001sag-00001 AND QLF001mag-00001 share member stem QLF001-00001!
    head = [
        "sampleid", "sample", "sagmag", "clusterid", "genomeid", "completeness", "contamination",
        "quality", "rrna5S", "rrna16S", "rrna23S", "trnaaa", "contig", "totallength", "n50",
        "gc", "cds", "gtdbr207", "plasmids", "phages", "arg",
    ]
    rows = [
        ["QLX001", "feces", "SAG", "cl_1", "QLF001sag-00001", "94.2", "1.1", "88.7", "1", "1", "1", "18", "1",
         str(len(s_sag_hq11)), str(len(s_sag_hq11)), "45.0", "24", "d__Bacteria;p__Bacteroidota;g__Prevotella;s__", "0", "0", "0"],
        ["QLX001", "saliva", "SAG", "cl_2", "QLF001sag-00002", "68.5", "2.0", "58.5", "0", "1", "0", "14", "1",
         str(len(s_sag_mq25)), str(len(s_sag_mq25)), "32.0", "24", "d__Bacteria;p__Patescibacteriota;c__Gracilibacteria;s__", "0", "0", "0"],
        ["QLX001", "saliva", "SAG", "cl_3", "QLF001sag-00003", "35.0", "1.0", "30.0", "0", "0", "0", "5", "1",
         str(len(s_sag_lq)), str(len(s_sag_lq)), "40.0", "12", "d__Bacteria;p__Firmicutes;s__", "0", "0", "0"],
        ["QLX001", "feces", "MAG", "cl_4", "QLF001mag-00001", "96.0", "0.8", "92.0", "1", "1", "1", "20", "1",
         str(len(s_mag_hq4)), str(len(s_mag_hq4)), "29.0", "24", "d__Bacteria;p__Mycoplasmatota;g__Mycoplasma;s__", "1", "0", "0"],
        ["QLX001", "feces", "MAG", "cl_5", "QLF001mag-00002", "85.0", "14.5", "12.5", "1", "0", "0", "12", "1",
         str(len(s_mag_lq)), str(len(s_mag_lq)), "44.0", "12", "d__Bacteria;p__Bacteroidota;s__", "0", "0", "0"],
    ]
    (dl / "BBSAG20" / "summary.tsv").write_text("\n".join("\t".join(r) for r in [head, *rows]) + "\n")

    with tarfile.open(dl / "BBSAG20" / "QLF001sag.tar.gz", "w:gz") as tf:
        for num, seq in (("00001", s_sag_hq11), ("00002", s_sag_mq25), ("00003", s_sag_lq)):
            for ext, payload in (
                (".fna.gz", _gz_bytes(_fa([(f"ctg_{num}", seq)]))),
                (".faa.gz", _gz_bytes(">p1\nMEFILPQ*\n")),
                (".ffn.gz", _gz_bytes(">c1 # 1 # 30 # 1 # ID=1_1\nATGAAATAA\n")),
                (".gff.gz", _gz_bytes("##gff-version 3\n")),
            ):
                ti = tarfile.TarInfo(f"QLF001/QLF001-{num}/QLF001-{num}{ext}")
                ti.size = len(payload)
                tf.addfile(ti, io.BytesIO(payload))

    with tarfile.open(dl / "BBSAG20" / "mag.tar.gz", "w:gz") as tf:
        for num, seq in (("00001", s_mag_hq4), ("00002", s_mag_lq)):
            for ext, payload in (
                (".fna.gz", _gz_bytes(_fa([(f"mctg_{num}", seq)]))),
                (".faa.gz", _gz_bytes(">p1\nMEFILPQ*\n")),
                (".gff.gz", _gz_bytes("##gff-version 3\n")),
            ):
                ti = tarfile.TarInfo(f"mag/QLF001/QLF001-{num}/QLF001-{num}{ext}")
                ti.size = len(payload)
                tf.addfile(ti, io.BytesIO(payload))

    # bbsag20_plasmid.fna.gz contains plasmid contigs extracted from MAGs (skipped as subset of MAG contigs)
    (dl / "BBSAG20" / "bbsag20_plasmid.fna.gz").write_bytes(_gz_bytes(_fa([("QLF001mag-00001_1", s_mag_hq4)])))

    # HAIMA_LR_SAG
    with zipfile.ZipFile(dl / "HAIMA_LR_SAG" / "SAGs.zip", "w") as zf:
        zf.writestr("SAGs/", "")
        zf.writestr("SAGs/HM_SAG_0001.fa", _fa([("hm_s1", s_haima_sag)]))
    with zipfile.ZipFile(dl / "HAIMA_LR_SAG" / "long-read-MAGs.zip", "w") as zf:
        zf.writestr("CP11-cbin2.fasta", _fa([("cbin2_ctg1", s_haima_cbin_dtr)]))
    with zipfile.ZipFile(dl / "HAIMA_LR_SAG" / "short-read-MAGs.zip", "w") as zf:
        zf.writestr("SR_bin1.fa", _fa([("sr1", s_haima_sr)]))

    # DEATHVALLEY_SAG (includes duplicate .xlsx entry and 1 exact repeat of s_haima_sag)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        with zipfile.ZipFile(dl / "DEATHVALLEY_SAG" / "BLM1-MAGs_SAGs_2011_2021.zip", "w") as zf:
            zf.writestr("AM-294-E03_contigs.fa", _fa([("am_ctg1", s_dv_sag)]))
            zf.writestr("AM-DUP_contigs.fa", _fa([("dup_ctg1", s_haima_sag)]))
            zf.writestr("information regarding the MAGs and SAGs of BLM1.xlsx", b"xlsx1")
            zf.writestr("information regarding the MAGs and SAGs of BLM1.xlsx", b"xlsx1")

    idx = tmp_path / "index_b"
    _write_r232_index(idx)

    owner_sag1 = tmp_path / "owners_b" / "SAG1"
    owner_sag1.mkdir(parents=True)
    (owner_sag1 / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(s_sag_hq11)]}]))

    return {
        "bbsag20_dirs": [dl / "BBSAG20"],
        "haima_dirs": [dl / "HAIMA_LR_SAG"],
        "dv_dirs": [dl / "DEATHVALLEY_SAG"],
        "index_dir": idx,
        "owners": {"SAG1": owner_sag1},
        "seqs": {
            "sag_hq11": s_sag_hq11, "sag_mq25": s_sag_mq25, "sag_lq": s_sag_lq,
            "mag_hq4": s_mag_hq4, "mag_lq": s_mag_lq,
            "haima_sag": s_haima_sag, "haima_cbin_core": s_haima_cbin_core,
            "haima_sr": s_haima_sr, "dv_sag": s_dv_sag,
        },
    }


def test_bbsag1_stage_filters_lq_screens_codes_handles_id_collision_splits_and_passes_s2_check(tmp_path: Path):
    mod = _mod()
    fx = _build_bbsag1_fixture(tmp_path)

    # Case A: single stage BBSAG1 (net-new bases <= max_stage_new_bp)
    rep = mod.stage_bbsag1(
        fx["bbsag20_dirs"], fx["haima_dirs"], fx["dv_dirs"], fx["index_dir"],
        tmp_path / "meta_single", tmp_path / "fasta_single",
        root=None, root_stages=(), owner_inputs=fx["owners"], later_inputs={},
        inputs_home=Path("/tmp/if2453/s2/inputs"), work=tmp_path / "work_single",
        max_stage_new_bp=40_000_000_000,
    )
    assert rep["stages_written"] == ["BBSAG1"]
    meta = tmp_path / "meta_single" / "BBSAG1"
    genomes = json.loads((meta / "genomes.json").read_text())
    # Kept: 2 MQ/HQ SAGs + 1 HQ MAG (BBSAG20) + 3 HAIMA + 1 unique DV (1 DV exact repeat of HAIMA SAG) = 7 distinct genomes
    assert len(genomes) == 7
    assert {g["kind"] for g in genomes} == {"prokaryote"}

    s = fx["seqs"]
    by_cu = {g["contig_uids"][0]: g for g in genomes}
    # Both QLF001sag-00001 and QLF001mag-00001 survived without stem collision, with codes 11, 25, 4
    assert by_cu[cuid(s["sag_hq11"])]["genetic_code"] == 11
    assert by_cu[cuid(s["sag_mq25"])]["genetic_code"] == 25
    assert by_cu[cuid(s["mag_hq4"])]["genetic_code"] == 4
    # LQ SAG and LQ MAG excluded into left_out.tsv
    assert cuid(s["sag_lq"]) not in by_cu
    assert cuid(s["mag_lq"]) not in by_cu
    left = _rows(meta / "left_out.tsv")
    assert {r["genome"] for r in left} == {"QLF001sag-00003", "QLF001mag-00002"}

    # HAIMA cbin2 50 bp exact DTR trimmed & circular
    assert by_cu[cuid(s["haima_cbin_core"])]["circular"] == [cuid(s["haima_cbin_core"])]
    assert by_cu[cuid(s["haima_cbin_core"])]["assembly_quality"] == "complete"

    # Redistributable split: BBSAG20+HAIMA (true) vs DEATHVALLEY (false)
    objs = {Path(o["path"]).name: o for o in json.loads((meta / "objects.json").read_text())}
    assert objs["bbsag1_sags_mags.fna.gz"]["redistributable"] is True
    assert objs["bbsag1_deathvalley.fna.gz"]["redistributable"] is False

    chk = subprocess.run(
        [sys.executable, str(BSP / "s2_check_inputs.py"), "BBSAG1",
         str(meta / "genomes.json"), str(meta / "objects_staging.json")],
        capture_output=True, text=True, check=False,
    )
    assert chk.returncode == 0, chk.stdout + chk.stderr
    assert json.loads(chk.stdout)["problems"] == []

    # Case B: deterministic split into BBSAG1A + BBSAG1B when net-new bases exceed cap (with workers=2 ParallelBgzfWriter)
    rep_split = mod.stage_bbsag1(
        fx["bbsag20_dirs"], fx["haima_dirs"], fx["dv_dirs"], fx["index_dir"],
        tmp_path / "meta_split", tmp_path / "fasta_split",
        root=None, root_stages=(), owner_inputs=fx["owners"], later_inputs={},
        inputs_home=Path("/tmp/if2453/s2/inputs"), work=tmp_path / "work_split",
        max_stage_new_bp=90_000, workers=2,
    )
    rep_split_1w = mod.stage_bbsag1(
        fx["bbsag20_dirs"], fx["haima_dirs"], fx["dv_dirs"], fx["index_dir"],
        tmp_path / "meta_split_1w", tmp_path / "fasta_split_1w",
        root=None, root_stages=(), owner_inputs=fx["owners"], later_inputs={},
        inputs_home=Path("/tmp/if2453/s2/inputs"), work=tmp_path / "work_split_1w",
        max_stage_new_bp=90_000, workers=1,
    )
    assert rep_split["stages_written"] == ["BBSAG1A", "BBSAG1B"]
    for st_name in ("BBSAG1A", "BBSAG1B"):
        s_2w = {k: v for k, v in rep_split["stages"][st_name]["sha256"].items() if k != "objects_staging.json"}
        s_1w = {k: v for k, v in rep_split_1w["stages"][st_name]["sha256"].items() if k != "objects_staging.json"}
        assert s_2w == s_1w
    for st_name in ("BBSAG1A", "BBSAG1B"):
        m_st = tmp_path / "meta_split" / st_name
        chk_st = subprocess.run(
            [sys.executable, str(BSP / "s2_check_inputs.py"), st_name,
             str(m_st / "genomes.json"), str(m_st / "objects_staging.json")],
            capture_output=True, text=True, check=False,
        )
        assert chk_st.returncode == 0, chk_st.stdout + chk_st.stderr
        assert json.loads(chk_st.stdout)["problems"] == []


def test_bbsag1_ani_records_96pct_matches_without_cutting_including_closed_cmags(tmp_path: Path):
    if shutil.which("skani") is None:
        pytest.skip("skani is not on PATH (GitHub Actions, or the tools-hidden run)")
    mod = _mod()
    fx = _build_bbsag1_fixture(tmp_path)
    mod.stage_bbsag1(
        fx["bbsag20_dirs"], fx["haima_dirs"], fx["dv_dirs"], fx["index_dir"],
        tmp_path / "meta_ani", tmp_path / "fasta_ani",
        root=None, root_stages=(), owner_inputs={}, later_inputs={},
        inputs_home=Path("/tmp/if2453/s2/inputs"), work=tmp_path / "work_ani",
        max_stage_new_bp=40_000_000_000,
    )
    meta = tmp_path / "meta_ani" / "BBSAG1"
    ref_meta = tmp_path / "ref_HROM1"
    ref_meta.mkdir()
    ref_fa = tmp_path / "hrom1_ref.fna"
    cbin_seq = fx["seqs"]["haima_cbin_core"]
    ref_fa.write_text(_fa([("hrom1_match", cbin_seq)]))
    (ref_meta / "records.tsv").write_text(f"genome\tcontig_uid\nHROM_REF_1\t{cuid(cbin_seq)}\n")

    staged_fas = [Path(o["path"]) for o in json.loads((meta / "objects_staging.json").read_text())]
    res = mod.bbsag1_ani(meta, staged_fas, tmp_path / "ani_out", whole={"HROM1": (ref_meta, [ref_fa])}, threads=1)
    assert res["meet"] == 1 and res["meet_by_stage"] == {"HROM1": 1}
    hits = _rows(tmp_path / "ani_out" / "ani_hits.tsv")
    assert len(hits) == 1 and hits[0]["genome"] == "CP11-cbin2"
    # Closed cMAG never removed on ANI match (DEC-84 (2))
    assert len(json.loads((meta / "genomes.json").read_text())) == 7


# ==================================================================================================================
# 5. Conf templates (PROMGE1.conf and BBSAG1.conf / BBSAG1A.conf + BBSAG1B.conf)
# ==================================================================================================================


def test_templates_pass_stage_est_conf_and_match_measurement_note():
    if not (REPO / NOTE).is_file():
        pytest.skip(f"{NOTE} not yet written (pre-staging fixture test pass)")
    est = _load("stage_est", "stage_est.py")
    bbsag_stages = ["BBSAG1"] if (TEMPLATES / "BBSAG1.conf").is_file() else ["BBSAG1A", "BBSAG1B"]
    stages = ["PROMGE1", *bbsag_stages]
    note_text = (REPO / NOTE).read_text(encoding="utf-8")
    for st in stages:
        tpl = TEMPLATES / f"{st}.conf"
        assert tpl.is_file(), f"missing template {tpl}"
        c = est.read_conf(tpl, st, tree=REPO)
        assert c.crispr is True
        assert str(c.inputs) == f"/tmp/if2453/s2/inputs/{st}"
        assert c.min_contig_len is None
        assert NOTE in c.notes
        assert c.genecall_mode == "default"
        assert (c.depth == "annotated" and c.library == ("pfam", "38.2", "gathering")) or (
            c.depth == "called" and c.library is None
        )
        sum_names = {n for _, n in c.sums}
        assert {"genomes.json", "objects.json"} < sum_names
        for h, _ in c.sums:
            assert h in note_text, f"{st}: SUM {h} missing from {NOTE}"
        for _, (n, _) in c.expect.items():
            assert n > 0
