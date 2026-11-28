"""Card C-86 (DEC-85 (4), DEC-80, DEC-88): stage light tower-landed r26.11.6 cohorts O4, GV2, and HROM1.

Independent oracles in this file:
- contig_uid and genome_uid computed directly via blake3;
- synthetic coding sequences built from sense-codon pools with planted internal TGA (code 4 Trp vs code 25 Gly)
  and internal TAG (code 16 Leu) codons;
- hand-written GenBank flatfiles (.gbff), tar.gz and zip archives, duplicate-column TSVs, BOM-prefixed CSVs, and
  Parquet r232_index tables;
- every test runs under both normal PATH and tools-hidden PATH=/usr/bin:/bin (skipping only live skani when absent).
"""

from __future__ import annotations

import bz2
import csv
import gzip
import importlib.util
import io
import json
import random
import shutil
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("C-86")

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
TEMPLATES = BSP / "stages" / "templates"
NOTE = "planning/measurements/2026-11-27_c86_o4_gv2_hrom1_inputs.md"


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
    return _load("c86_o4_gv2_hrom1_inputs", "c86_o4_gv2_hrom1_inputs.py")


@pytest.fixture(params=["normal-path", "tools-hidden"], autouse=True)
def path_mode(request, monkeypatch):
    if request.param == "tools-hidden":
        monkeypatch.setenv("PATH", "/usr/bin:/bin")
    return request.param


def cuid(seq: str | bytes) -> str:
    b = seq.encode("ascii") if isinstance(seq, str) else seq
    return blake3(b.upper()).digest(length=16).hex()


def guid(seqs: list[str | bytes]) -> str:
    raws = sorted(bytes.fromhex(cuid(s)) for s in seqs)
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


def _gz_lines(path: Path) -> list[str]:
    return gzip.decompress(path.read_bytes()).decode("utf-8").splitlines()


def _rows(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ------------------------------------------------------------------------------------------------------------------
# Synthetic coding-sequence generator (independent oracle for pyrodigal / pyrodigal-gv code screens)
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
    tag_every: int = 0,
    tgg_every: int = 85,
) -> str:
    """Build a synthetic ~25 kb genome of dense ORFs separated by short RBS-bearing intergenic spacers."""
    rng = random.Random(seed)
    spacer = "TAATTTTAAAGGAGATATACAT"  # stop + AT-rich spacer + Shine-Dalgarno
    chunks: list[str] = ["AATTCCGGTTAA"]
    for _ in range(n_genes):
        orf = ["ATG"]
        for k in range(1, codons_per_gene):
            if tga_every and k % tga_every == 0:
                orf.append("TGA")
            elif tag_every and k % tag_every == 0:
                orf.append("TAG")
            elif tgg_every and k % tgg_every == 0:
                orf.append("TGG")
            else:
                orf.append(rng.choice(BASE_SENSE_CODONS))
        orf.append("TAA")
        chunks.append("".join(orf))
        chunks.append(spacer)
    return "".join(chunks)


def format_gbff_record(
    accession: str,
    seq: str,
    *,
    moltype: str = "DNA",
    circular: bool = True,
    definition: str = "Protist organelle complete genome.",
    organism: str = "Protista sp.",
    taxonomy: str = "Eukaryota; Sar; Stramenopiles.",
    organelle: str | None = "mitochondrion",
    transl_table: int | None = 4,
    spliced_cds: int = 0,
    unspliced_cds: int = 2,
    introns: int = 0,
) -> str:
    topo = "circular" if circular else "linear"
    lines = [
        f"LOCUS       {accession.split('.')[0]:<16} {len(seq)} bp    {moltype:<7} {topo:<8} PLN 27-NOV-2026",
        f"DEFINITION  {definition}",
        f"ACCESSION   {accession.split('.')[0]}",
        f"VERSION     {accession}",
        f"SOURCE      {organism}",
        f"  ORGANISM  {organism}",
        f"            {taxonomy}",
        "FEATURES             Location/Qualifiers",
        f"     source          1..{len(seq)}",
        f'                     /organism="{organism}"',
    ]
    if organelle is not None:
        lines.append(f'                     /organelle="{organelle}"')
    for i in range(introns):
        s = 50 + i * 120
        lines.append(f"     intron          {s}..{s + 40}")
    for i in range(spliced_cds):
        s = 100 + i * 300
        lines.append(f"     CDS             join({s}..{s + 90},{s + 150}..{s + 240})")
        if transl_table is not None:
            lines.append(f"                     /transl_table={transl_table}")
        lines.append('                     /product="spliced protein"')
    for i in range(unspliced_cds):
        s = 100 + (spliced_cds + i) * 300
        e = min(len(seq), s + 210)
        lines.append(f"     CDS             {s}..{e}")
        if transl_table is not None:
            lines.append(f"                     /transl_table={transl_table}")
        lines.append('                     /product="cytochrome oxidase"')
    lines.append("ORIGIN")
    low = seq.lower()
    for pos in range(0, len(low), 60):
        row = low[pos : pos + 60]
        groups = " ".join(row[j : j + 10] for j in range(0, len(row), 10))
        lines.append(f"{pos + 1:9d} {groups}")
    lines.append("//")
    return "\n".join(lines) + "\n"


# ==================================================================================================================
# 1. Magic-byte format detection & archive census
# ==================================================================================================================


def test_magic_bytes_and_archive_census(tmp_path: Path):
    mod = _mod()
    gz_fna = tmp_path / "sample.fna.gz"
    gz_fna.write_bytes(_gz_bytes(_fa([("c1 circular", "ACGTACGT" * 20)])))
    bz_file = tmp_path / "sample.bz2"
    bz_file.write_bytes(bz2.compress(b"hello"))
    gbff_file = tmp_path / "sample.gbff"
    gbff_file.write_text(format_gbff_record("KX000001.1", "ACGT" * 200))

    zpath = tmp_path / "mixed.zip"
    with zipfile.ZipFile(zpath, "w") as zf:
        zf.writestr("bins/", "")
        zf.writestr("bins/good.fna", _fa([("g1_1", "ACGTURYKMSWBDHVN" * 10)]))
        zf.writestr("bins/prot.faa", ">p1\nMEFILPQZ*\n")
        zf.writestr("bins/aln.fasta", ">a1\nACGT--ACGT\n")
        zf.writestr("bins/bad.fna", ">b1\nACGT1234ACGT\n")
        zf.writestr("meta.xlsx", b"PK_fake_sheet")

    tpath = tmp_path / "mixed.tar.gz"
    with tarfile.open(tpath, "w:gz") as tf:
        data = _fa([("v1", "ACGT" * 100)]).encode()
        ti = tarfile.TarInfo("viruses/v1.fna")
        ti.size = len(data)
        tf.addfile(ti, io.BytesIO(data))

    assert mod.magic_format(gz_fna) == "gzip"
    assert mod.magic_format(bz_file) == "bzip2"
    assert mod.magic_format(zpath) == "zip"
    assert mod.magic_format(tpath) == "gzip"
    assert mod.magic_format(gbff_file) == "gbff"

    rep = mod.census([gz_fna, zpath, tpath, gbff_file])
    assert rep[str(gz_fna)]["dna_records"] == 1
    assert rep[str(gz_fna)]["members"][0]["circular_in_header"] == 1
    zmem = {m["name"]: m["alphabet"] for m in rep[str(zpath)]["members"]}
    assert zmem == {
        "bins/": "directory",
        "bins/good.fna": "dna",
        "bins/prot.faa": "protein",
        "bins/aln.fasta": "alignment",
        "bins/bad.fna": "non_iupac",
        "meta.xlsx": "not_fasta",
    }
    assert rep[str(tpath)]["dna_records"] == 1 and rep[str(tpath)]["dna_bases"] == 400
    assert rep[str(gbff_file)]["gbff_records"] == 1


# ==================================================================================================================
# 2. Genetic code screens (1, 3, 4, 11, 16, 22, 25)
# ==================================================================================================================


def test_genetic_code_screens_across_standard_organellar_and_recoded_tables():
    mod = _mod()
    seq_std = synth_coding_genome(11, tga_every=0, tag_every=0, tgg_every=80).encode()
    # Code 4: internal TGA ~1.25% (every 80th codon), TGG absent -> UGA=Trp
    seq_c4 = synth_coding_genome(14, tga_every=80, tag_every=0, tgg_every=0).encode()
    # Code 25: internal TGA ~6.25% (every 16th codon = Gly), TGG ~1.1% (every 90th codon = Trp) -> UGA=Gly
    seq_c25 = synth_coding_genome(25, tga_every=16, tag_every=0, tgg_every=90).encode()
    # Code 16: internal TAG ~5% (every 20th codon = Leu), no internal TGA -> UAG=Leu
    seq_c16 = synth_coding_genome(16, tga_every=0, tag_every=20, tgg_every=80).encode()

    for code in (1, 3, 4, 11, 16, 22, 25):
        mod.check_supported_code(code)

    # Prokaryote screen (HROM1): catches code 11, code 4 even on Equivita (no tax regex match), and code 25 on CPR
    r11 = mod.screen_prokaryote_code([seq_std], "d__Bacteria;p__Bacteroidota;g__Prevotella;s__")
    assert r11["genetic_code"] == 11
    r4 = mod.screen_prokaryote_code([seq_c4], "d__Bacteria;p__Actinomycetota;g__Equivita;s__Equivita altericodex")
    assert r4["genetic_code"] == 4
    r25 = mod.screen_prokaryote_code([seq_c25], "d__Bacteria;p__Patescibacteriota;c__Gracilibacteria;s__")
    assert r25["genetic_code"] == 25

    # Fungal mitogenome screen (O4 O3_FUNGAL): distinguishes standard (1), mold mt (4), and TAG->Leu mt (16)
    f1 = mod.screen_fungal_mt_code(seq_std, "Ascomycota")
    f4 = mod.screen_fungal_mt_code(seq_c4, "Ascomycota")
    f16 = mod.screen_fungal_mt_code(seq_c16, "Chytridiomycota")
    assert f1["genetic_code"] in (1, 4)  # standard ORFs have no internal stops under either 1 or 4; 1 preferred when tied
    assert f1["genetic_code"] == 1
    assert f4["genetic_code"] == 4
    assert f16["genetic_code"] == 16

    # Giant virus code-25 check (GV2): rescues true code-25 (TGA=Gly, TGG=Trp), keeps code-4 and code-11 as gv-meta
    gv_std = mod.check_gv_code25([seq_std])
    gv_c4 = mod.check_gv_code25([seq_c4])
    gv_c25 = mod.check_gv_code25([seq_c25])
    assert gv_std["use_c25"] is False
    assert gv_c4["use_c25"] is False
    assert gv_c25["use_c25"] is True


# ==================================================================================================================
# 3. O4 staging (O4_ORG .gbff/.fna chunks + O3_FUNGAL .fasta/.tsv)
# ==================================================================================================================


def _build_o4_fixture(tmp_path: Path) -> dict:
    dl = tmp_path / "dl_o4"
    for sub in ("O4_CRYP", "O4_PROT_MT", "O4_PROT_PT", "O4_WEIRD_MT", "O3_FUNGAL"):
        (dl / sub).mkdir(parents=True)

    s_pt11 = synth_coding_genome(101, n_genes=6, codons_per_gene=200)
    s_nm1 = synth_coding_genome(102, n_genes=6, codons_per_gene=200)
    s_mt4 = synth_coding_genome(103, n_genes=24, codons_per_gene=320, tga_every=60, tgg_every=0)
    s_mt16 = synth_coding_genome(104, n_genes=8, codons_per_gene=200, tag_every=25)
    s_mt22 = synth_coding_genome(105, n_genes=8, codons_per_gene=200)
    s_nuc = synth_coding_genome(106, n_genes=8, codons_per_gene=200)
    s_rna = "ACGT" * 300
    s_bad = "ACGT" * 100 + "XXXX" + "ACGT" * 100
    s_fung4 = synth_coding_genome(107, n_genes=24, codons_per_gene=320, tga_every=50, tgg_every=0)
    s_fung16 = synth_coding_genome(108, n_genes=24, codons_per_gene=320, tag_every=20, tgg_every=80)

    # O4_CRYP: 1 circular plastid (table 11), 1 linear nucleomorph (table 1), 1 spliced nuclear chromosome (left out)
    (dl / "O4_CRYP" / "O4_CRYP_chunk_0001.fna").write_text(
        _fa([
            ("M81884.1 Epifagus virginiana chloroplast, complete genome", s_pt11),
            ("CP000001.1 Guillardia theta nucleomorph chromosome 1", s_nm1),
            ("NC_999999.1 Cryptomonas nuclear chromosome 1", s_nuc),
        ])
    )
    (dl / "O4_CRYP" / "O4_CRYP_chunk_0001.gbff").write_text(
        format_gbff_record("M81884.1", s_pt11, circular=True, definition="Epifagus virginiana chloroplast",
                           organism="Epifagus virginiana", organelle="plastid", transl_table=11)
        + format_gbff_record("CP000001.1", s_nm1, circular=False, definition="Guillardia theta nucleomorph chr 1",
                             organism="Guillardia theta", organelle="nucleomorph", transl_table=1)
        + format_gbff_record("NC_999999.1", s_nuc, circular=False, definition="Cryptomonas nuclear chromosome 1",
                             organism="Cryptomonas paramecium", organelle=None, transl_table=1,
                             spliced_cds=5, unspliced_cds=0, introns=5)
    )

    # O4_PROT_MT: 1 circular table-4 mtDNA, 1 RNA record (left out), 1 non-IUPAC record (left out)
    (dl / "O4_PROT_MT" / "O4_PROT_MT_chunk_0001.fna").write_text(
        _fa([
            ("KT000001.1 Cafeteria roenbergensis mitochondrion, complete genome", s_mt4),
            ("XM_000001.1 Protist mRNA", s_rna),
            ("BAD00001.1 Bad bases mitochondrion", s_bad),
        ])
    )
    (dl / "O4_PROT_MT" / "O4_PROT_MT_chunk_0001.gbff").write_text(
        format_gbff_record("KT000001.1", s_mt4, circular=True, organelle="mitochondrion", transl_table=4)
        + format_gbff_record("XM_000001.1", s_rna, moltype="mRNA", circular=False, organelle=None, transl_table=1)
        + format_gbff_record("BAD00001.1", s_bad, circular=False, organelle="mitochondrion", transl_table=4)
    )

    # O4_PROT_PT: 1 linear table-16 mitochondrion/plastid held by O2
    (dl / "O4_PROT_PT" / "O4_PROT_PT_chunk_0001.fna").write_text(
        _fa([("JQ000016.1 Pedinomonas mitochondrion", s_mt16)])
    )
    (dl / "O4_PROT_PT" / "O4_PROT_PT_chunk_0001.gbff").write_text(
        format_gbff_record("JQ000016.1", s_mt16, circular=False, organelle="mitochondrion", transl_table=16)
    )

    # O4_WEIRD_MT: 1 circular table-22 mitochondrion held by RMD1
    (dl / "O4_WEIRD_MT" / "O4_WEIRD_MT_chunk_0001.fna").write_text(
        _fa([("HM000022.1 Scenedesmus obliquus mitochondrion", s_mt22)])
    )
    (dl / "O4_WEIRD_MT" / "O4_WEIRD_MT_chunk_0001.gbff").write_text(
        format_gbff_record("HM000022.1", s_mt22, circular=True, organelle="mitochondrion", transl_table=22)
    )

    # O3_FUNGAL: 3 records (1 code-4 new, 1 exact repeat of KT000001.1, 1 code-16 new)
    fungal_fa = dl / "O3_FUNGAL" / "Assembled_fungal_mitochondrial_genomes.fasta"
    fungal_fa.write_text(
        _fa([
            ("SRR100001 [organism=Aspergillus_sp] [Location=mitochondrion] [topology=circular] [completeness=complete]",
             s_fung4),
            ("SRR100002 [organism=Dup_sp] [Location=mitochondrion] [topology=circular] [completeness=complete]",
             s_mt4),
            ("SRR100003 [organism=Spizellomyces_sp] [Location=mitochondrion] [topology=circular]",
             s_fung16),
        ])
    )
    fungal_tsv = dl / "O3_FUNGAL" / "All_fungal_mitochondrial_genomes_accessions.tsv"
    fungal_tsv.write_text(
        "seq_id\ttax_id\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tgenome_size\tGC\tsource\tNCBI_accession_number\n"
        f"SRR100001\t5052\tFungi\tAscomycota\tEurotiomycetes\tEurotiales\tAspergillaceae\tAspergillus\tAspergillus sp.\t{len(s_fung4)}\t0.35\tSRA\tNA\n"
        f"SRR100002\t5052\tFungi\tAscomycota\tEurotiomycetes\tEurotiales\tAspergillaceae\tAspergillus\tDup sp.\t{len(s_mt4)}\t0.35\tNCBI\tKT000001\n"
        f"SRR100003\t4805\tFungi\tChytridiomycota\tSpizellomycetes\tSpizellomycetales\tSpizellomycetaceae\tSpizellomyces\tSpizellomyces sp.\t{len(s_fung16)}\t0.35\tSRA\tNA\n"
    )

    owner_o2 = tmp_path / "owners" / "O2"
    owner_o2.mkdir(parents=True)
    (owner_o2 / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(s_mt16)]}]))
    owner_rmd1 = tmp_path / "owners" / "RMD1"
    owner_rmd1.mkdir(parents=True)
    (owner_rmd1 / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(s_mt22)]}]))

    return {
        "o4_org_dir": dl,
        "fungal_fa": fungal_fa,
        "fungal_tsv": fungal_tsv,
        "owners": {"O2": owner_o2, "RMD1": owner_rmd1},
        "seqs": {
            "pt11": s_pt11, "nm1": s_nm1, "mt4": s_mt4, "mt16": s_mt16, "mt22": s_mt22,
            "fung4": s_fung4, "fung16": s_fung16,
        },
        "tmp": tmp_path,
    }


def test_o4_stage_routes_organellar_codes_drops_nuclear_and_rna_and_passes_s2_check(tmp_path: Path):
    mod = _mod()
    fx = _build_o4_fixture(tmp_path)
    rep1 = mod.stage_o4(
        fx["o4_org_dir"], fx["fungal_fa"], fx["fungal_tsv"],
        tmp_path / "meta_a", tmp_path / "fasta_a",
        root=None, root_stages=(), owner_inputs=fx["owners"], later_inputs={},
        inputs_home=Path("/tmp/if2453/s2/inputs"), work=tmp_path / "work_a",
    )
    rep2 = mod.stage_o4(
        fx["o4_org_dir"], fx["fungal_fa"], fx["fungal_tsv"],
        tmp_path / "meta_b", tmp_path / "fasta_b",
        root=None, root_stages=(), owner_inputs=fx["owners"], later_inputs={},
        inputs_home=Path("/tmp/if2453/s2/inputs"), work=tmp_path / "work_b",
    )

    meta = tmp_path / "meta_a" / "O4"
    # Under 64 MiB: both FASTA objects stay in meta_out/O4/, nothing under fasta_out/O4/
    fa_dir = tmp_path / "fasta_a" / "O4"
    assert not fa_dir.exists() or list(fa_dir.iterdir()) == []
    assert (meta / "o4_ncbi_organelles.fna.gz").is_file()
    assert (meta / "o4_fungal_mitogenomes.fna.gz").is_file()

    genomes = json.loads((meta / "genomes.json").read_text())
    # 5 kept from O4_ORG + 2 unique kept from O3_FUNGAL (1 repeat SRR100002 listed in repeated_records.tsv) = 7 genomes
    assert len(genomes) == 7 and rep1["genomes"] == 7 and rep1["repeated_genomes"] == 1
    by_cu = {g["contig_uids"][0]: g for g in genomes}
    s = fx["seqs"]
    assert by_cu[cuid(s["pt11"])]["genetic_code"] == 11 and by_cu[cuid(s["pt11"])]["circular"] == [cuid(s["pt11"])]
    assert by_cu[cuid(s["nm1"])]["genetic_code"] == 1 and by_cu[cuid(s["nm1"])]["circular"] == []
    assert by_cu[cuid(s["mt4"])]["genetic_code"] == 4 and by_cu[cuid(s["mt4"])]["circular"] == [cuid(s["mt4"])]
    assert by_cu[cuid(s["mt16"])]["genetic_code"] == 16 and by_cu[cuid(s["mt16"])]["circular"] == []
    assert by_cu[cuid(s["mt22"])]["genetic_code"] == 22 and by_cu[cuid(s["mt22"])]["circular"] == [cuid(s["mt22"])]
    assert by_cu[cuid(s["fung4"])]["genetic_code"] == 4
    assert by_cu[cuid(s["fung16"])]["genetic_code"] == 16

    left = _rows(meta / "left_out.tsv")
    left_accs = {r["accession"]: r["reason"] for r in left}
    assert set(left_accs) == {"NC_999999.1", "XM_000001.1", "BAD00001.1"}
    assert "nuclear" in left_accs["NC_999999.1"]
    assert "moltype" in left_accs["XM_000001.1"]
    assert "non_iupac" in left_accs["BAD00001.1"]

    rep_rows = _rows(meta / "repeated_records.tsv")
    assert len(rep_rows) == 1 and rep_rows[0]["genome"] == "SRR100002"

    shared = {(r["genome"], r["owner_stage"]) for r in _rows(meta / "shared_records.tsv")}
    assert shared == {("JQ000016.1", "O2"), ("HM000022.1", "RMD1")}
    assert rep1["held_by"] == {"O2": 1, "RMD1": 1}

    objs = {Path(o["path"]).name: o for o in json.loads((meta / "objects.json").read_text())}
    assert objs["o4_ncbi_organelles.fna.gz"]["redistributable"] is True
    assert objs["o4_fungal_mitogenomes.fna.gz"]["redistributable"] is False

    # Byte-for-byte determinism across runs
    for f in ("genomes.json", "objects.json", "records.tsv", "left_out.tsv", "repeated_records.tsv",
              "shared_records.tsv", "o4_ncbi_organelles.fna.gz", "o4_fungal_mitogenomes.fna.gz"):
        assert (meta / f).read_bytes() == (tmp_path / "meta_b" / "O4" / f).read_bytes(), f
    assert {k: v for k, v in rep1["sha256"].items() if k != "objects_staging.json"} == {
        k: v for k, v in rep2["sha256"].items() if k != "objects_staging.json"
    }

    chk = subprocess.run(
        [sys.executable, str(BSP / "s2_check_inputs.py"), "O4",
         str(meta / "genomes.json"), str(meta / "objects_staging.json")],
        capture_output=True, text=True, check=False,
    )
    assert chk.returncode == 0, chk.stdout + chk.stderr
    assert json.loads(chk.stdout)["problems"] == []


# ==================================================================================================================
# 4. GV2 staging (6 polar/freshwater/WWTP giant-virus dirs)
# ==================================================================================================================


def _build_gv2_fixture(tmp_path: Path) -> dict:
    dl = tmp_path / "dl_gv2"
    for sub in ("GV_POLAR_2026", "GV_POLAR_PITOT", "GV_POLAR_GEVE", "GV_GKS", "GV_POLAR_MARINE", "GV_WWTP"):
        (dl / sub).mkdir(parents=True)

    c_ant1 = synth_coding_genome(201, n_genes=4, codons_per_gene=150)
    c_ant2 = synth_coding_genome(202, n_genes=4, codons_per_gene=150)
    c_pat = synth_coding_genome(203, n_genes=4, codons_per_gene=150)
    c_pitot = synth_coding_genome(204, n_genes=4, codons_per_gene=150)
    c_gks_gv = synth_coding_genome(205, n_genes=4, codons_per_gene=150)
    c_gks_bac = synth_coding_genome(206, n_genes=4, codons_per_gene=150)
    c_vr_core = synth_coding_genome(207, n_genes=4, codons_per_gene=150)
    c_vr_dtr = c_vr_core + c_vr_core[:35]  # 35 bp exact terminal repeat
    c_icel_vr = synth_coding_genome(208, n_genes=4, codons_per_gene=150)
    c_icel_host = synth_coding_genome(209, n_genes=6, codons_per_gene=200)
    c_ncv = synth_coding_genome(210, n_genes=4, codons_per_gene=150)
    c_plv = synth_coding_genome(211, n_genes=4, codons_per_gene=150)

    # 1. GV_POLAR_2026
    with zipfile.ZipFile(dl / "GV_POLAR_2026" / "GVMAGs_from_polar_regions.zip", "w") as zf:
        zf.writestr("GVMAGs_from_polar_regions/Antarctic_seawater/ANT_P3_gi-53.fna",
                    _fa([("c1", c_ant1), ("c2", c_ant2)]))
        zf.writestr("GVMAGs_from_polar_regions/Source_data.xlsx", b"xlsx_bytes")

    # 2. GV_POLAR_PITOT
    with tarfile.open(dl / "GV_POLAR_PITOT" / "GVMAGs_fna.tar.gz", "w:gz") as tf:
        b = _fa([("pitot_c1", c_pitot)]).encode()
        ti = tarfile.TarInfo("GVMAGs_fna/Pitot_ARC_1.fna")
        ti.size = len(b)
        tf.addfile(ti, io.BytesIO(b))

    # 3. GV_GKS: tomas_scMAGs (viral) + all_MAGs (host/prokaryotic MAGs to exclude when exclude_gks_all_mags=True)
    with tarfile.open(dl / "GV_GKS" / "GKS_genomes.tar.gz", "w:gz") as tf:
        for mname, seq in (("GKS_genomes/tomas_scMAGs/GKS1_GV.fna", c_gks_gv),
                           ("GKS_genomes/all_MAGs/GKS1_MB2.1.fltr.fa", c_gks_bac)):
            b = _fa([("ctg1", seq)]).encode()
            ti = tarfile.TarInfo(mname)
            ti.size = len(b)
            tf.addfile(ti, io.BytesIO(b))

    # 4. GV_POLAR_GEVE: confirmed excised VRs + host chromosome zip (host zip left out for C-60)
    (dl / "GV_POLAR_GEVE" / "fcylindrus_vrconfirmed.fna").write_text(_fa([("fcyl_vr1", c_vr_dtr)]))
    (dl / "GV_POLAR_GEVE" / "ICEL_confirmed_VR_noflank.fna").write_text(_fa([("ICEL_vr1", c_icel_vr)]))
    with zipfile.ZipFile(dl / "GV_POLAR_GEVE" / "all_ICEL_viral_contigs.zip", "w") as zf:
        zf.writestr("all_ICEL_viral_contigs.fna", _fa([("JACBWV010000610.1", c_icel_host)]))

    # 5. GV_POLAR_MARINE: exact repeat of ANT_P3_gi-53 + PAT_7_gi-38 held by GV1
    (dl / "GV_POLAR_MARINE" / "ANT_P3_gi-53.fa").write_text(_fa([("c1", c_ant1), ("c2", c_ant2)]))
    (dl / "GV_POLAR_MARINE" / "PAT_7_gi-38.fa").write_text(_fa([("pat_c1", c_pat)]))

    # 6. GV_WWTP: ncv/contigs + vph_plv/contigs kept; annotations/prokka and proteins skipped
    with tarfile.open(dl / "GV_WWTP" / "genomes.tar.gz", "w:gz") as tf:
        items = [
            ("genomes/ncv/contigs/NCV_bin1.fasta", _fa([("ncv_c1", c_ncv)]).encode()),
            ("genomes/vph_plv/contigs/PLV_1.fna", _fa([("plv_c1", c_plv)]).encode()),
            ("genomes/vph_plv/annotations/prokka/PLV_1/PLV_1.ffn", _fa([("cds1", "ATGAAATAA")]).encode()),
            ("genomes/vph_plv/annotations/prokka/PLV_1/PLV_1.fna", _fa([("plv_c1", c_plv)]).encode()),
            ("genomes/ncv/proteins/NCV_bin1.faa", b">p1\nMEFILPQ*\n"),
        ]
        for mname, b in items:
            ti = tarfile.TarInfo(mname)
            ti.size = len(b)
            tf.addfile(ti, io.BytesIO(b))

    owner_gv1 = tmp_path / "owners" / "GV1"
    owner_gv1.mkdir(parents=True)
    (owner_gv1 / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(c_pat)]}]))

    return {
        "dirs": {sub: dl / sub for sub in (
            "GV_POLAR_2026", "GV_POLAR_PITOT", "GV_POLAR_GEVE", "GV_GKS", "GV_POLAR_MARINE", "GV_WWTP"
        )},
        "owners": {"GV1": owner_gv1},
        "seqs": {"vr_core": c_vr_core, "pat": c_pat, "icel_host": c_icel_host, "gks_bac": c_gks_bac},
        "tmp": tmp_path,
    }


def test_gv2_stage_filters_prokka_and_host_contigs_deduplicates_and_passes_s2_check(tmp_path: Path):
    mod = _mod()
    fx = _build_gv2_fixture(tmp_path)
    rep1 = mod.stage_gv2(
        fx["dirs"], tmp_path / "meta_a", tmp_path / "fasta_a",
        root=None, root_stages=(), owner_inputs=fx["owners"], later_inputs={},
        inputs_home=Path("/tmp/if2453/s2/inputs"), work=tmp_path / "work_a",
        exclude_gks_all_mags=True,
    )
    rep2 = mod.stage_gv2(
        fx["dirs"], tmp_path / "meta_b", tmp_path / "fasta_b",
        root=None, root_stages=(), owner_inputs=fx["owners"], later_inputs={},
        inputs_home=Path("/tmp/if2453/s2/inputs"), work=tmp_path / "work_b",
        exclude_gks_all_mags=True,
    )
    meta = tmp_path / "meta_a" / "GV2"
    genomes = json.loads((meta / "genomes.json").read_text())
    # Staged: ANT_P3_gi-53 (1), Pitot_ARC_1 (1), GKS1_GV (1), fcyl_vr1 (1), ICEL_vr1 (1), PAT_7_gi-38 (1),
    # NCV_bin1 (1), PLV_1 (1) = 8 distinct viral genomes; ANT_P3_gi-53 in GV_POLAR_MARINE is 1 exact repeat
    assert len(genomes) == 8 and rep1["genomes"] == 8 and rep1["repeated_genomes"] == 1
    assert {g["kind"] for g in genomes} == {"virus"}

    # Verify 35 bp exact DTR on fcyl_vr1 trimmed & circular under DEC-65 (1)
    vr_cu = cuid(fx["seqs"]["vr_core"])
    vr_g = next(g for g in genomes if g["contig_uids"] == [vr_cu])
    assert vr_g["circular"] == [vr_cu]

    # Host nuclear contig & prokaryotic GKS all_MAGs recorded in left_out.tsv and absent from genomes.json
    all_cus = {u for g in genomes for u in g["contig_uids"]}
    assert cuid(fx["seqs"]["icel_host"]) not in all_cus
    assert cuid(fx["seqs"]["gks_bac"]) not in all_cus
    left = _rows(meta / "left_out.tsv")
    reasons = {r["reason"] for r in left}
    assert "eukaryote_host_nuclear_contig_c60" in reasons
    assert "prokaryote_co_binned_mag_not_gv" in reasons

    # Prokka annotations and proteins skipped in skipped_members.tsv
    skipped_m = _rows(meta / "skipped_members.tsv")
    assert any("prokka" in r["member"] for r in skipped_m)

    assert rep1["held_by"] == {"GV1": 1}
    objs = json.loads((meta / "objects.json").read_text())
    assert all(o["redistributable"] is False for o in objs)

    for f in ("genomes.json", "objects.json", "records.tsv", "left_out.tsv", "skipped_members.tsv",
              "repeated_records.tsv", "shared_records.tsv", "trims.tsv"):
        assert (meta / f).read_bytes() == (tmp_path / "meta_b" / "GV2" / f).read_bytes(), f
    assert {k: v for k, v in rep1["sha256"].items() if k != "objects_staging.json"} == {
        k: v for k, v in rep2["sha256"].items() if k != "objects_staging.json"
    }

    chk = subprocess.run(
        [sys.executable, str(BSP / "s2_check_inputs.py"), "GV2",
         str(meta / "genomes.json"), str(meta / "objects_staging.json")],
        capture_output=True, text=True, check=False,
    )
    assert chk.returncode == 0, chk.stdout + chk.stderr
    assert json.loads(chk.stdout)["problems"] == []


# ==================================================================================================================
# 5. HROM1 staging (HROM_ORAL + CODE4_GLMC + CODE4_ACTINO1, codes 11/4/25, r232_status, skani)
# ==================================================================================================================


def _write_r232_index(d: Path) -> None:
    d.mkdir(parents=True)
    gt = [
        ("GCA_056388515", "GCA_056388515", True, "GCA_056388515", "s__Equivita altericodex"),
        ("GCA_000999995", "GCA_000999995", False, "GCA_000999991", "s__Zinderia insecticola"),
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
            "globdb_id": ["GCA_000999995"],
            "insdc_accession": ["GCA_000999995"],
            "original_id": [None],
            "cluster_id": [None],
        }),
        d / "globdb_r232.parquet",
    )


def _build_hrom1_fixture(tmp_path: Path) -> dict:
    dl = tmp_path / "dl_hrom1"
    for sub in ("HROM_ORAL", "CODE4_GLMC", "CODE4_ACTINO1/assemblies"):
        (dl / sub).mkdir(parents=True)

    s_oral_std = synth_coding_genome(301, n_genes=24, codons_per_gene=320)
    s_oral_cpr25 = synth_coding_genome(302, n_genes=24, codons_per_gene=320, tga_every=16, tgg_every=90)
    s_oral_mq = synth_coding_genome(303, n_genes=24, codons_per_gene=320)
    s_glmc_c4 = synth_coding_genome(304, n_genes=24, codons_per_gene=320, tga_every=60, tgg_every=0)
    s_actino_c4 = synth_coding_genome(305, n_genes=24, codons_per_gene=320, tga_every=60, tgg_every=0)

    # 1. HROM_ORAL: note duplicate "species" column in header (col 0 and col 12)!
    hrom_tar = dl / "HROM_ORAL" / "HROM_Genome.All.tar.gz"
    with tarfile.open(hrom_tar, "w:gz") as tf:
        for gid, seq in (
            ("HROM_Genome_0001", s_oral_std),
            ("HROM_Genome_0002", s_oral_cpr25),
            ("HROM_Genome_0003", s_oral_mq),
        ):
            b = _fa([(f"{gid}_ctg1", seq)]).encode()
            ti = tarfile.TarInfo(f"HROM_representative_Genomes/{gid}.fna")
            ti.size = len(b)
            tf.addfile(ti, io.BytesIO(b))

    hrom_meta = dl / "HROM_ORAL" / "HROM-Species-metadata.tsv"
    hrom_head = [
        "species", "completeness", "contamination", "N50", "genome length", "gtdb_taxonomy",
        "domain", "phylum", "class", "order", "family", "genus", "species", "species_status",
        "isolate status", "gut-oral shared status", "count_MQ", "count_NC", "count_total",
        "continent", "country", "Isolation source",
    ]
    hrom_rows = [
        ["HROM_Genome_0001", "98.5", "0.5", "50000", str(len(s_oral_std)),
         "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella oralis",
         "d__Bacteria", "p__Bacteroidota", "c__Bacteroidia", "o__Bacteroidales", "f__Bacteroidaceae", "g__Prevotella",
         "s__Prevotella oralis/HROM_Genome_0001", "NC/HQ", "Isolate", "not-shared", "0", "1", "1", "USA", "North America", "Oral"],
        ["HROM_Genome_0002", "92.0", "1.1", "40000", str(len(s_oral_cpr25)),
         "d__Bacteria;p__Patescibacteriota;c__Gracilibacteria;o__Absconditabacterales;f__SR1;g__SR1;s__",
         "d__Bacteria", "p__Patescibacteriota", "c__Gracilibacteria", "o__Absconditabacterales", "f__SR1", "g__SR1",
         "s__HROM_Genome_0002/HROM_Genome_0002", "NC/HQ", "No", "not-shared", "0", "1", "1", "China", "Asia", "Dental"],
        ["HROM_Genome_0003", "75.0", "3.0", "25000", str(len(s_oral_mq)),
         "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__",
         "d__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Lactobacillales", "f__Streptococcaceae", "g__Streptococcus",
         "s__HROM_Genome_0003/HROM_Genome_0003", "MQ", "No", "not-shared", "1", "0", "1", "UK", "Europe", "Saliva"],
    ]
    hrom_meta.write_text("\n".join("\t".join(r) for r in [hrom_head, *hrom_rows]) + "\n")

    # 2. CODE4_GLMC: zip + UTF-8 BOM metadata_mags.csv (real header columns)
    glmc_zip = dl / "CODE4_GLMC" / "merged_MAGs.zip"
    with zipfile.ZipFile(glmc_zip, "w") as zf:
        zf.writestr("merged/", "")
        zf.writestr("merged/GLMC_Bin012.fasta", _fa([("zind_ctg1", s_glmc_c4)]))
    glmc_csv = dl / "CODE4_GLMC" / "metadata_mags.csv"
    glmc_csv.write_bytes(
        (
            "\ufeffShort_name,Sample,Host,Subfamily,Classification_GTDB_226,Genome_Size(bp),Contig_Count,N50,CDS,GC(%),tRNA,rRNA,checkm_completeness,checkm_contamination,checkm2_completeness,checkm2_contamination,MIMAG\n"
            f"GLMC_Bin012,GCA_000999995.1,Dalbulus maidis,Deltocephalinae,d__Bacteria;p__Pseudomonadota;g__Zinderia;s__Zinderia insecticola,{len(s_glmc_c4)},1,{len(s_glmc_c4)},98,18.6,29,2,80.0,0.1,96.0,0.2,HQ\n"
        ).encode("utf-8")
    )

    # 3. CODE4_ACTINO1: manifest (real column 'accession') + .fna.gz (complete closed chromosome -> never cut by ANI under DEC-84 (2))
    act_fna = dl / "CODE4_ACTINO1" / "assemblies" / "GCA_056388515.1_genomic.fna.gz"
    act_fna.write_bytes(_gz_bytes(_fa([("CP100001.1 Equivita altericodex chromosome, complete genome", s_actino_c4)])))
    act_man = dl / "CODE4_ACTINO1" / "MANIFEST_CODE4_ACTINO1_PRJNA1366607.tsv"
    act_man.write_text(
        "stage\tbioproject\taccession\torganism_name\tbase_count\tassembly_level\turl\n"
        f"CODE4_ACTINO1\tPRJNA1366607\tGCA_056388515.1\tEquivita altericodex\t{len(s_actino_c4)}\tComplete Genome\thttps://ftp.ncbi.nlm.nih.gov/x\n"
    )

    idx = tmp_path / "index"
    _write_r232_index(idx)

    owner_dpn1 = tmp_path / "owners" / "DPN1"
    owner_dpn1.mkdir(parents=True)
    (owner_dpn1 / "genomes.json").write_text(json.dumps([{"contig_uids": [cuid(s_oral_std)]}]))

    return {
        "hrom_tar": hrom_tar,
        "hrom_meta": hrom_meta,
        "glmc_zip": glmc_zip,
        "glmc_csv": glmc_csv,
        "actino_dir": dl / "CODE4_ACTINO1",
        "index_dir": idx,
        "owners": {"DPN1": owner_dpn1},
        "seqs": {
            "oral_std": s_oral_std, "oral_cpr25": s_oral_cpr25, "oral_mq": s_oral_mq,
            "glmc_c4": s_glmc_c4, "actino_c4": s_actino_c4,
        },
        "tmp": tmp_path,
    }


def test_hrom1_stage_screens_codes_4_11_25_records_r232_status_and_passes_s2_check(tmp_path: Path):
    mod = _mod()
    fx = _build_hrom1_fixture(tmp_path)
    rep1 = mod.stage_hrom1(
        fx["hrom_tar"], fx["hrom_meta"], fx["glmc_zip"], fx["glmc_csv"], fx["actino_dir"], fx["index_dir"],
        tmp_path / "meta_a", tmp_path / "fasta_a",
        root=None, root_stages=(), owner_inputs=fx["owners"], later_inputs={},
        inputs_home=Path("/tmp/if2453/s2/inputs"), work=tmp_path / "work_a",
    )
    rep2 = mod.stage_hrom1(
        fx["hrom_tar"], fx["hrom_meta"], fx["glmc_zip"], fx["glmc_csv"], fx["actino_dir"], fx["index_dir"],
        tmp_path / "meta_b", tmp_path / "fasta_b",
        root=None, root_stages=(), owner_inputs=fx["owners"], later_inputs={},
        inputs_home=Path("/tmp/if2453/s2/inputs"), work=tmp_path / "work_b",
    )

    meta = tmp_path / "meta_a" / "HROM1"
    genomes = json.loads((meta / "genomes.json").read_text())
    # All genomes kept (none removed for GTDB R232 or GlobDB r232 membership!)
    assert len(genomes) == 5 and rep1["genomes"] == 5
    by_cu = {g["contig_uids"][0]: g for g in genomes}
    s = fx["seqs"]
    assert by_cu[cuid(s["oral_std"])]["genetic_code"] == 11
    assert by_cu[cuid(s["oral_cpr25"])]["genetic_code"] == 25
    assert by_cu[cuid(s["glmc_c4"])]["genetic_code"] == 4
    assert by_cu[cuid(s["actino_c4"])]["genetic_code"] == 4
    assert by_cu[cuid(s["actino_c4"])]["assembly_quality"] == "complete"
    assert by_cu[cuid(s["oral_std"])]["assembly_quality"] == "mag_hq"
    assert by_cu[cuid(s["oral_mq"])]["assembly_quality"] == "mag_mq"

    # r232_status.tsv records GTDB R232 & GlobDB r232 status without cutting any genome
    st = {r["genome"]: r for r in _rows(meta / "r232_status.tsv")}
    assert st["GCA_056388515.1"]["r232_status"] == "in_gtdb_r232"
    assert st["GLMC_Bin012"]["r232_status"] == "in_gtdb_and_globdb_r232"
    assert st["HROM_Genome_0001"]["r232_status"] == "not_visible"

    assert rep1["held_by"] == {"DPN1": 1}
    objs = {Path(o["path"]).name: o for o in json.loads((meta / "objects.json").read_text())}
    assert objs["hrom1_code4_actino1.fna.gz"]["redistributable"] is True
    assert objs["hrom1_oral_glmc.fna.gz"]["redistributable"] is False

    for f in ("genomes.json", "objects.json", "records.tsv", "genetic_codes.tsv", "r232_status.tsv",
              "repeated_records.tsv", "shared_records.tsv"):
        assert (meta / f).read_bytes() == (tmp_path / "meta_b" / "HROM1" / f).read_bytes(), f

    chk = subprocess.run(
        [sys.executable, str(BSP / "s2_check_inputs.py"), "HROM1",
         str(meta / "genomes.json"), str(meta / "objects_staging.json")],
        capture_output=True, text=True, check=False,
    )
    assert chk.returncode == 0, chk.stdout + chk.stderr
    assert json.loads(chk.stdout)["problems"] == []


def test_hrom1_ani_records_96pct_matches_without_cutting_including_closed_genomes(tmp_path: Path):
    if shutil.which("skani") is None:
        pytest.skip("skani is not on PATH (GitHub Actions, or the tools-hidden run)")
    mod = _mod()
    fx = _build_hrom1_fixture(tmp_path)
    mod.stage_hrom1(
        fx["hrom_tar"], fx["hrom_meta"], fx["glmc_zip"], fx["glmc_csv"], fx["actino_dir"], fx["index_dir"],
        tmp_path / "meta_s", tmp_path / "fasta_s",
        root=None, root_stages=(), owner_inputs={}, later_inputs={},
        inputs_home=Path("/tmp/if2453/s2/inputs"), work=tmp_path / "work_s",
    )
    meta = tmp_path / "meta_s" / "HROM1"
    ref_meta = tmp_path / "ref_NX1"
    ref_meta.mkdir()
    ref_fa = tmp_path / "nx1.fna"
    # Plant a reference matching the closed Equivita genome (GCA_056388515.1) at 100% ANI
    act_seq = fx["seqs"]["actino_c4"]
    ref_fa.write_text(_fa([("nx1_equivita_copy", act_seq)]))
    (ref_meta / "records.tsv").write_text(f"genome\tcontig_uid\nNX1_BIN_1\t{cuid(act_seq)}\n")

    staged_fas = [Path(o["path"]) for o in json.loads((meta / "objects_staging.json").read_text())]
    res = mod.hrom1_ani(
        meta, staged_fas, tmp_path / "ani_out",
        whole={"NX1": (ref_meta, [ref_fa])}, threads=1,
    )
    assert res["meet"] == 1 and res["meet_by_stage"] == {"NX1": 1}
    hits = _rows(tmp_path / "ani_out" / "hrom1_ani_hits.tsv")
    assert len(hits) == 1 and hits[0]["genome"] == "GCA_056388515.1"
    # Confirm genomes.json still holds all 5 genomes (closed genome never cut per DEC-84 (2))
    assert len(json.loads((meta / "genomes.json").read_text())) == 5


# ==================================================================================================================
# 6. Conf templates (O4.conf, GV2.conf, HROM1.conf)
# ==================================================================================================================


@pytest.mark.parametrize("stage,after", [("O4", ("GBRMGD",)), ("GV2", ("O4",)), ("HROM1", ("GV2",))])
def test_templates_pass_stage_est_conf_and_match_measurement_note(stage: str, after: tuple[str, ...]):
    est = _load("stage_est", "stage_est.py")
    c = est.read_conf(TEMPLATES / f"{stage}.conf", stage, tree=REPO)
    assert c.after == after
    assert c.crispr is True
    assert str(c.inputs) == f"/tmp/if2453/s2/inputs/{stage}"
    assert c.min_contig_len is None
    assert NOTE in c.notes
    if stage == "O4":
        assert c.depth == "annotated" and c.library == ("pfam", "38.2", "gathering")
    elif stage == "GV2":
        assert (c.depth == "annotated" and c.library == ("pfam", "38.2", "gathering")) or (
            c.depth == "called" and c.library is None
        )
    else:
        assert c.depth == "called" and c.library is None
    if stage == "GV2":
        assert c.genecall_mode in ("gv-meta", "gv-meta-c25")
    else:
        assert c.genecall_mode == "default"
    sum_names = {n for _, n in c.sums}
    assert {"genomes.json", "objects.json"} < sum_names
    note_text = (REPO / NOTE).read_text(encoding="utf-8")
    for h, _ in c.sums:
        assert h in note_text, f"{stage}: SUM {h} missing from {NOTE}"
    for owner, (n, _) in c.expect.items():
        assert n > 0
