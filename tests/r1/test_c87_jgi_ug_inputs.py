"""Card C-87 (DEC-80, DEC-88): stream-census and stage pod-resident JGI IUIG and IMG-UMAG prokaryote archives.

Every test constructs independent ground truth in `tmp_path` (BLAKE3 digests directly from `blake3`, synthetic
prokaryotic coding sequences with standard and alternative UGA-readthrough genetic codes `11`, `25`, `4`, corrupt
unselected archive payloads proving zero extraction of unselected entries, and independent brute-force filters).
Each test runs under both the normal `PATH` and `PATH=/usr/bin:/bin` so GitHub Actions (without pod binaries) passes.
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
import tracemalloc
import zipfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("C-87")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
TEMPLATES = BSP / "stages" / "templates"
NOTE = "planning/measurements/2026-11-27_c87_jgi_ug_inputs.md"


def _load(name: str, file: str):
    if str(BSP) not in sys.path:
        sys.path.insert(0, str(BSP))
    target = BSP / file
    spec = importlib.util.spec_from_file_location(name, target)
    if spec is None or spec.loader is None or not target.is_file():
        raise FileNotFoundError(target)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _mod():
    return _load("c87_jgi_ug_inputs", "c87_jgi_ug_inputs.py")


@pytest.fixture(params=["normal-path", "tools-hidden"], autouse=True)
def path_mode(request, monkeypatch):
    if request.param == "tools-hidden":
        monkeypatch.setenv("PATH", "/usr/bin:/bin")
    return request.param


def cuid(seq: bytes) -> str:
    return blake3(seq.upper()).digest(length=16).hex()


def guid(seq_list: list[bytes]) -> str:
    seen: list[bytes] = []
    for s in seq_list:
        u = blake3(s.upper()).digest(length=16)
        if u not in seen:
            seen.append(u)
    h = blake3()
    for u in seen:
        h.update(u)
    return h.digest(length=16).hex()


_SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in ("TAA", "TAG", "TGA")]


def synth_contig(
    seed: int,
    *,
    n_genes: int = 115,
    codons: int = 300,
    uga_per_gene: int = 0,
    uga_codon_bias: str = "gly",
) -> bytes:
    """Synthetic prokaryotic contig (~108 kb by default) with RBS + ATG + sense codons + TAA.

    When `uga_per_gene > 0`, inserts in-frame internal `TGA` codons inside every ORF so code 11 truncates ORFs
    while code 25 (`TGA` -> Gly) and code 4 (`TGA` -> Trp) read through full-length genes.
    """
    rng = random.Random(seed)
    parts = ["A" * 40]
    for _ in range(n_genes):
        spacer = "".join(rng.choice("ACGT") for _ in range(7))
        cods = [rng.choice(_SENSE) for _ in range(codons - 2)]
        if uga_per_gene:
            step = max(1, len(cods) // (uga_per_gene + 1))
            for k in range(1, uga_per_gene + 1):
                idx = min(len(cods) - 1, k * step)
                cods[idx] = "TGA"
                if uga_codon_bias == "trp" and idx + 1 < len(cods):
                    cods[idx + 1] = "TGG"
        parts.append("".join(rng.choice("AT") for _ in range(25)) + "AGGAGG" + spacer + "ATG" + "".join(cods) + "TAA")
    parts.append("T" * 40)
    return "".join(parts).encode("ascii")


L = 80 + 115 * (38 + 3 * 300)


def _fa_bytes(records: list[tuple[str, bytes]], width: int = 80) -> bytes:
    out: list[bytes] = []
    for hdr, seq in records:
        out.append(b">" + hdr.encode("utf-8") + b"\n")
        for i in range(0, len(seq), width):
            out.append(seq[i : i + width] + b"\n")
    return b"".join(out)


def _gz_bytes(raw: bytes | str) -> bytes:
    if isinstance(raw, str):
        raw = raw.encode("utf-8")
    buf = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buf, mtime=0) as gz:
        gz.write(raw)
    return buf.getvalue()


def _inner_tgz(files: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz", format=tarfile.GNU_FORMAT) as tf:
        for name, data in files.items():
            ti = tarfile.TarInfo(name=name)
            ti.size = len(data)
            ti.mtime = 0
            tf.addfile(ti, io.BytesIO(data))
    return buf.getvalue()


def _nested_jgi_zip(zip_path: Path, inner_tar_name: str, members: list[tuple[str, bytes]]) -> None:
    tar_buf = io.BytesIO()
    with tarfile.open(fileobj=tar_buf, mode="w", format=tarfile.GNU_FORMAT) as tf:
        for mname, data in members:
            ti = tarfile.TarInfo(name=mname)
            ti.size = len(data)
            ti.mtime = 0
            tf.addfile(ti, io.BytesIO(data))
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("Download_0000001_File_Manifest.csv", "filename,size\n")
        zf.writestr(inner_tar_name, tar_buf.getvalue())


def _build_r232_fixture(tmp_path: Path) -> Path:
    r232 = _load("r232_index", "r232_index.py")
    gt = tmp_path / "GTDB_R232_META"
    gl = tmp_path / "GLOBDB_R232_TABLES"
    r6 = tmp_path / "r226"
    idx = tmp_path / "r232_index"
    for d in (gt, gl, r6):
        d.mkdir(parents=True, exist_ok=True)

    meta_head = "\t".join([
        "accession",
        "gtdb_genome_representative",
        "gtdb_representative",
        "gtdb_taxonomy",
        "ncbi_genbank_assembly_accession",
        "ncbi_biosample",
        "ncbi_wgs_master",
        "ncbi_taxid",
        "ncbi_organism_name",
    ])
    ar_rows = [
        meta_head,
        "GB_GCA_000000011.1\tGB_GCA_000000011.1\tt\td__Archaea;p__Thermoproteota;c__C;o__O;f__F;g__G;s__Archaeum known\tGCA_000000011.1\tSAMN00000011\tAAAA01000000\t111111\tArchaeum known",
    ]
    bac_rows = [
        meta_head,
        "RS_GCF_000000033.2\tRS_GCF_000000033.2\tt\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\tGCA_000000033.2\tSAMN00000033\tBBBB01000000\t511145\tEscherichia coli K12",
    ]
    (gt / "ar53_metadata_r232.tsv.gz").write_bytes(_gz_bytes("\n".join(ar_rows) + "\n"))
    (gt / "bac120_metadata_r232.tsv.gz").write_bytes(_gz_bytes("\n".join(bac_rows) + "\n"))

    glob_tax = "\n".join([
        "GCA_000000011\td__Archaea;p__Thermoproteota;c__C;o__O;f__F;g__G;s__Archaeum known",
        "GCF_000000033\td__Bacteria;p__Pseudomonadota;c__C;o__O;f__F;g__Escherichia;s__Escherichia coli",
        "GEMOTU_00001\td__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__F;g__G;s__Bacteroidota gem1",
    ]) + "\n"
    (gl / "globdb_r232_taxonomy.tsv.gz").write_bytes(_gz_bytes(glob_tax))
    ck = "\n".join([
        "Name\tCompleteness\tContamination\tGenome_Size\tTotal_Contigs",
        "GCA_000000011\t95.0\t0.5\t1800000\t5",
        "GCF_000000033\t99.0\t0.1\t4600000\t1",
        "GEMOTU_00001\t92.0\t1.0\t2400000\t20",
    ]) + "\n"
    (gl / "globdb_r232_checkm2.tsv.gz").write_bytes(_gz_bytes(ck))
    (gl / "globdb_r232_dataset_list.tsv").write_text(
        "Name\tIdentifier\tTitle\tPaper link\tSpecies reps in dataset\tNumber of genomes in GlobDB r232\n"
        "GTDB\tGCA & GCF\tt\tl\t2\t2\n"
        "GEM\tGEMOTU\tt\tl\t2\t1\n"
    )
    dict_buf = io.BytesIO()
    with tarfile.open(fileobj=dict_buf, mode="w:gz") as tf:
        gem_dict = (
            "GlobDB_ID\totu_id\tgenome_id\n"
            "GEMOTU_00001\tOTU-1\t3300000001_1\n"
            "GEMOTU_00002\tOTU-2\t3300000001_2\n"
        ).encode("utf-8")
        ti = tarfile.TarInfo("dicts/GEMOTU_dictionary.tsv")
        ti.size = len(gem_dict)
        tf.addfile(ti, io.BytesIO(gem_dict))
    (gl / "globdb_r232_dictionaries.tar.gz").write_bytes(dict_buf.getvalue())

    (r6 / "globdb_r226_taxonomy.tsv.gz").write_bytes(_gz_bytes("GEMOTU_00001\td__Bacteria;p__Bacteroidota;s__x\n"))
    (r6 / "globdb_r226_dictionaries.tar.gz").write_bytes(dict_buf.getvalue())
    r232.build(gl, gt, r6, idx)
    return idx


def _read_tsv(path: Path) -> list[dict]:
    opener = (lambda p: gzip.open(p, "rt", encoding="utf-8", newline="")) if path.suffix == ".gz" else (
        lambda p: open(p, encoding="utf-8", newline="")
    )
    with opener(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_zip_infolist_and_o1_seek_without_extracting_unselected_members(tmp_path: Path):
    m = _mod()
    s_keep = synth_contig(11, n_genes=10)

    # 1. Direct zip mode: unselected member has corrupt gzip bytes; O(1) ZipFile.open(info) must only open selected!
    direct_zip = tmp_path / "direct.zip"
    with zipfile.ZipFile(direct_zip, "w", compression=zipfile.ZIP_STORED) as zf:
        zf.writestr("UNSELECTED_CORRUPT.fna.gz", b"\x1f\x8b\x08\x00CORRUPT_PAYLOAD_DO_NOT_READ")
        zf.writestr("KEEP_001.fna.gz", _gz_bytes(_fa_bytes([("KEEP_001_c1", s_keep)])))

    info_summary = m.inspect_zip_entries(direct_zip)
    assert info_summary["layout"] == "direct_members"
    assert info_summary["genome_member_count"] == 2

    got_direct = list(m.iter_selected_from_zip(direct_zip, {"KEEP_001"}))
    assert len(got_direct) == 1
    gid, recs = got_direct[0]
    assert gid == "KEEP_001"
    assert len(recs) == 1 and recs[0][1] == s_keep

    # 2. Nested JGI zip -> outer .tar -> inner <id>.tar.gz mode: unselected inner .tar.gz is corrupt!
    nested_zip = tmp_path / "nested.zip"
    good_tgz = _inner_tgz({"637000099/637000099.fna": _fa_bytes([("scaf99", s_keep)])})
    corrupt_tgz = b"NOT_GZIP_UNSELECTED_INNER_TAR_GZ_MUST_BE_SKIPPED"
    _nested_jgi_zip(
        nested_zip,
        "Custom_MPIIUIG/2023-11-14_1.0/aggregate_tarball.tar",
        [("637000001.tar.gz", corrupt_tgz), ("637000099.tar.gz", good_tgz)],
    )
    nested_summary = m.inspect_zip_entries(nested_zip)
    assert nested_summary["layout"] == "nested_tar"
    assert nested_summary["inner_tar"] == "Custom_MPIIUIG/2023-11-14_1.0/aggregate_tarball.tar"

    got_nested = list(m.iter_selected_from_zip(nested_zip, {"637000099"}))
    assert len(got_nested) == 1
    assert got_nested[0][0] == "637000099"
    assert got_nested[0][1][0][1] == s_keep


def test_metadata_crossjoin_r232_index_and_novelty_prioritization(tmp_path: Path):
    m = _mod()
    idx = _build_r232_fixture(tmp_path)

    # Metadata zip: IUIG + uMAG
    meta_zip = tmp_path / "jgi_inovirus_umag_iuig_meta_20261118.zip"
    iuig_meta = "\n".join([
        "taxon_oid\ttaxon_display_name\tdomain\tphylum\tir_class\tir_order\tfamily\tgenus\tspecies\tseq_status\tseq_center\tanalysis_project_id\tncbi_taxon_id\tis_public\tadd_date",
        # 1: Known GTDB R232 bacterium via GCA crosswalk -> recorded as in_gtdb_and_globdb_r232, not picked for JGIUG1
        "637000001\tEscherichia coli K12\tBacteria\tPseudomonadota\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae\tEscherichia\tEscherichia coli\tFinished\tJGI\tGa0001\t511145\tYes\t01-JAN-20",
        # 2: Novel HQ Archaeon -> top priority for JGIUG1!
        "280000001\tUncultivated Asgard archaeon\tArchaea\tCandidatus Asgardarchaeota\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified\tPermanent Draft\tJGI\tGa0002\t\tYes\t01-JAN-23",
        # 3: Novel MQ common-phylum bacterium -> spillover into JGIUG_left_for_later.tsv.gz when cap binds!
        "280000002\tUncultivated Pseudomonadota sp\tBacteria\tPseudomonadota\tAlphaproteobacteria\tunclassified\tunclassified\tunclassified\tunclassified\tPermanent Draft\tJGI\tGa0003\t\tYes\t01-JAN-23",
        # 4: Virus row -> non-prokaryote skipped
        "280000003\tUncultivated phage\tViruses\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified\tPermanent Draft\tJGI\tGa0004\t\tYes\t01-JAN-23",
        # 5: Eukaryote row -> non-prokaryote deferred to C-60
        "280000004\tUncultivated alga\tEukaryota\tChlorophyta\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified\tPermanent Draft\tJGI\tGa0005\t\tYes\t01-JAN-23",
    ]) + "\n"
    umag_meta = "\n".join([
        "taxon_oid\ttaxon_display_name\tdomain\tphylum\tir_class\tir_order\tfamily\tgenus\tspecies\tseq_status\tseq_center\tanalysis_project_id\tncbi_taxon_id\tis_public\tadd_date",
        "3300000001\tMarine vent metagenome\t*Microbiome\tEnvironmental\tAquatic\tThermal springs\tVent\tUnclassified\tVent\tPermanent Draft\tJGI\tGa0100\t\tYes\t01-JAN-23",
        "3300000002\tSubsurface aquifer metagenome\t*Microbiome\tEnvironmental\tTerrestrial\tDeep subsurface\tAquifer\tUnclassified\tAquifer\tPermanent Draft\tJGI\tGa0200\t\tYes\t01-JAN-24",
    ]) + "\n"
    with zipfile.ZipFile(meta_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("Custom_MPIIUIG/2023-11-14_1.0/img_genomes_tarballs_metadata_111023.txt", iuig_meta)
        zf.writestr("CusIMGUMAG/IMG_Data/unrestricted_metagenomes_for_bin_tarballs.02062024.txt", umag_meta)

    # IMG/M 70-col tar mapping 637000001 -> GCF_000000033.2
    img_m_tar = tmp_path / "img_m_meta.tar"
    img_70 = _gz_bytes(
        "IMG Genome ID\tNCBI Assembly Accession\tNCBI Biosample Accession\tUncultured Type\tGenome Size * assembled\tScaffold Count * assembled\tEcosystem\tEcosystem Category\n"
        "637000001\tGCF_000000033.2\tSAMN00000033\t\t4600000\t1\tHost-associated\tHuman\n"
        "280000001\t\t\tSingle Cell\t120000\t1\tEnvironmental\tAquatic\n"
        "280000002\t\t\tMetagenomic\t120000\t1\tEnvironmental\tSoil\n"
    )
    with tarfile.open(img_m_tar, "w") as tf:
        ti = tarfile.TarInfo("img_m/img_m_prok_70col.tsv.gz")
        ti.size = len(img_70)
        tf.addfile(ti, io.BytesIO(img_70))

    # GEMs OTU + CheckM tar: 3300000001_1 is GEMOTU_00001 in GlobDB r232; 3300000001_2 is OTU-2 collapsed in r232 dict;
    # 280000001 has HQ CheckM; 280000002 has MQ CheckM.
    gems_tar = tmp_path / "gems_meta.tar"
    g2o = (
        "otu_id\tgenome_id\ttype\tstudy\n"
        "OTU-1\t3300000001_1\tMAG\tGEM\n"
        "OTU-2\t3300000001_2\tMAG\tGEM\n"
    ).encode("utf-8")
    g2c = (
        "genome_id\tcompleteness\tcontamination\tquality_score\tpass_qc\n"
        "280000001\t94.0\t1.0\t89.0\t1\n"
        "280000002\t65.0\t3.0\t50.0\t1\n"
    ).encode("utf-8")
    with tarfile.open(gems_tar, "w") as tf:
        for nm, raw in [("gems/otus/genome_to_otu.tsv", g2o), ("gems/otus/genome_to_checkm.tsv", g2c)]:
            ti = tarfile.TarInfo(nm)
            ti.size = len(raw)
            tf.addfile(ti, io.BytesIO(raw))

    s_arch = synth_contig(21, n_genes=20)
    s_cpr = synth_contig(22, n_genes=20, uga_per_gene=6)
    s_pseu = synth_contig(23, n_genes=20)
    s_gem = synth_contig(24, n_genes=20)

    iuig_zip = tmp_path / "iuig.zip"
    _nested_jgi_zip(
        iuig_zip,
        "Custom_MPIIUIG/2023-11-14_1.0/aggregate_tarball.tar",
        [
            ("637000001.tar.gz", b"CORRUPT_KNOWN_GTDB_SKIPPED"),
            ("280000001.tar.gz", _inner_tgz({"280000001/280000001.fna": _fa_bytes([("a1", s_arch)])})),
            ("280000002.tar.gz", _inner_tgz({"280000002/280000002.fna": _fa_bytes([("p1", s_pseu)])})),
            ("280000003.tar.gz", b"CORRUPT_VIRUS_SKIPPED"),
            ("280000004.tar.gz", b"CORRUPT_EUK_SKIPPED"),
        ],
    )

    mbin1 = (
        "# IMG TAXON ID: 3300000001\n"
        "IMG Bin ID\tBin Quality\tBin Lineage\tGTDB-TK lineage\tBin Completeness\tBin Contamintation\tTotal Number of Bases\tNumber of genes\tNum of 5s rRNA\tNum of 16s rRNA\tNum of 23s rRNA\tNum of tRNA\tIMG scaffold IDs of members\n"
        f"3300000001_1\tHQ\tBacteria;Bacteroidetes\td__Bacteria;p__Bacteroidota;c__Bacteroidia;o__O;f__F;g__G;s__Bacteroidota gem1\t92.0\t1.0\t{len(s_gem)}\t20\t1\t1\t1\t20\ts1\n"
        f"3300000001_2\tMQ\tBacteria;Bacteroidetes\td__Bacteria;p__Bacteroidota;c__Bacteroidia;o__O;f__F;g__G;s__\t70.0\t2.0\t{len(s_gem)}\t20\t0\t1\t0\t15\ts2\n"
    ).encode("utf-8")
    mbin2 = (
        "# IMG TAXON ID: 3300000002\n"
        "IMG Bin ID\tBin Quality\tBin Lineage\tGTDB-TK lineage\tBin Completeness\tBin Contamintation\tTotal Number of Bases\tNumber of genes\tNum of 5s rRNA\tNum of 16s rRNA\tNum of 23s rRNA\tNum of tRNA\tIMG scaffold IDs of members\n"
        f"3300000002_1\tHQ\tBacteria;Candidatus Pacebacteria\td__Bacteria;p__Patescibacteriota;c__Paceibacteria;o__O;f__F;g__G;s__\t91.5\t0.5\t{len(s_cpr)}\t20\t1\t1\t1\t22\tcpr1\n"
        "3300000002_2\tLQ\tBacteria;Proteobacteria\td__Bacteria;p__Pseudomonadota;c__C;o__O;f__F;g__G;s__\t32.0\t1.0\t50000\t10\t0\t0\t0\t4\tlq1\n"
    ).encode("utf-8")

    umag_zip = tmp_path / "umag.zip"
    _nested_jgi_zip(
        umag_zip,
        "CusIMGUMAG/IMG_Data/IMG_prokaryote_binning_files.tar",
        [
            ("./3300000001_1.tar.gz", _inner_tgz({"3300000001_1.fna": _fa_bytes([("s1", s_gem)]), "mbin_datafile_3300000001.txt": mbin1})),
            ("./3300000001_2.tar.gz", b"CACHED_MBIN_NO_NEED_TO_DECOMPRESS_DURING_CENSUS"),
            ("./3300000002_1.tar.gz", _inner_tgz({"3300000002_1.fna": _fa_bytes([("cpr1", s_cpr)]), "mbin_datafile_3300000002.txt": mbin2})),
            ("./3300000002_2.tar.gz", b"LQ_CACHED_MBIN_SKIPPED"),
        ],
    )

    census_dir = tmp_path / "census"
    # Cap JGIUG1 at 2 genomes' bases so HQ Archaea (280000001) + HQ CPR (3300000002_1) enter JGIUG1,
    # while MQ Pseudomonadota (280000002) spills over into JGIUG_left_for_later.tsv.gz!
    cap_bp = len(s_arch) + len(s_cpr)
    rep = m.census_and_select(
        iuig_zip=iuig_zip,
        umag_zip=umag_zip,
        meta_zip=meta_zip,
        r232_index_dir=idx,
        gtdb_meta_dir=tmp_path / "GTDB_R232_META",
        gems_meta_tar=gems_tar,
        img_m_tar=img_m_tar,
        out_dir=census_dir,
        max_stage_bp=cap_bp,
    )

    assert rep["iuig"]["non_prokaryote_rows"] == {"Eukaryota": 1, "Viruses": 1}
    assert rep["selected"]["genomes"] == 2
    sel_rows = _read_tsv(census_dir / "selection.tsv")
    # Archaea HQ ranked #1, CPR/Patescibacteriota HQ ranked #2!
    assert [r["genome_id"] for r in sel_rows] == ["280000001", "3300000002_1"]

    left_rows = _read_tsv(census_dir / "JGIUG_left_for_later.tsv.gz")
    left_by_id = {r["genome_id"]: r for r in left_rows}
    # 280000002 is spillover_cap; 637000001, 3300000001_1, 3300000001_2 are r232-covered HQ/MQ reserves
    assert "280000002" in left_by_id
    assert left_by_id["280000002"]["defer_reason"] == "spillover_stage_cap"
    assert left_by_id["3300000001_1"]["r232_status"] == "in_globdb_r232"
    assert left_by_id["637000001"]["r232_status"] == "in_gtdb_and_globdb_r232"
    assert "3300000002_2" not in left_by_id  # LQ (<50% completeness) excluded from HQ/MQ reserve

    # Verify reselect_census preserves deterministic selection and enforces strict not_in_r232
    rep_res = m.reselect_census(census_dir, max_stage_bp=cap_bp, min_unscored_bp=10_000)
    assert rep_res["selected"]["genomes"] == 2
    assert [r["genome_id"] for r in _read_tsv(census_dir / "selection.tsv")] == ["280000001", "3300000002_1"]


def test_genetic_code_screening_iupac_validation_dec27_dedup_and_s2_check(tmp_path: Path):
    m = _mod()
    idx = _build_r232_fixture(tmp_path)

    s_arch_shared = synth_contig(101)  # shared with root P2
    s_arch_new = synth_contig(102)
    s_cpr_25 = synth_contig(103, uga_per_gene=8)  # UGA readthrough -> code 25
    s_later_held = synth_contig(104)  # held by RMD1 (later r26.11.5 stage)
    s_bad_iupac = b"ACGT" * 200 + b"---" + b"ACGT" * 200  # gap bytes -> refused by DEC-84 (1)

    meta_zip = tmp_path / "meta.zip"
    iuig_meta = "\n".join([
        "taxon_oid\ttaxon_display_name\tdomain\tphylum\tir_class\tir_order\tfamily\tgenus\tspecies\tseq_status\tseq_center\tanalysis_project_id\tncbi_taxon_id\tis_public\tadd_date",
        "290000001\tArch 1\tArchaea\tThermoproteota\tC\tO\tF\tG\tunclassified\tPermanent Draft\tJGI\tGa1\t\tYes\t01-JAN-23",
        "290000002\tArch 1 exact repeat\tArchaea\tThermoproteota\tC\tO\tF\tG\tunclassified\tPermanent Draft\tJGI\tGa2\t\tYes\t01-JAN-23",
        "290000003\tBad gap genome\tArchaea\tThermoproteota\tC\tO\tF\tG\tunclassified\tPermanent Draft\tJGI\tGa3\t\tYes\t01-JAN-23",
    ]) + "\n"
    umag_meta = (
        "taxon_oid\ttaxon_display_name\tdomain\tphylum\tir_class\tir_order\tfamily\tgenus\tspecies\tseq_status\tseq_center\tanalysis_project_id\tncbi_taxon_id\tis_public\tadd_date\n"
        "3300000099\tDeep aquifer\t*Microbiome\tEnvironmental\tAquatic\tGroundwater\tSubsurface\tUnclassified\tSubsurface\tPermanent Draft\tJGI\tGa99\t\tYes\t01-JAN-24\n"
    )
    with zipfile.ZipFile(meta_zip, "w") as zf:
        zf.writestr("Custom_MPIIUIG/2023-11-14_1.0/img_genomes_tarballs_metadata_111023.txt", iuig_meta)
        zf.writestr("CusIMGUMAG/IMG_Data/unrestricted_metagenomes_for_bin_tarballs.02062024.txt", umag_meta)

    iuig_zip = tmp_path / "iuig.zip"
    _nested_jgi_zip(
        iuig_zip,
        "Custom_MPIIUIG/2023-11-14_1.0/aggregate_tarball.tar",
        [
            ("290000001.tar.gz", _inner_tgz({"290000001/290000001.fna": _fa_bytes([("c_sh", s_arch_shared), ("c_new", s_arch_new)])})),
            ("290000002.tar.gz", _inner_tgz({"290000002/290000002.fna": _fa_bytes([("c_sh2", s_arch_shared), ("c_new2", s_arch_new)])})),
            ("290000003.tar.gz", _inner_tgz({"290000003/290000003.fna": _fa_bytes([("c_bad", s_bad_iupac)])})),
        ],
    )

    mbin99 = (
        "# IMG TAXON ID: 3300000099\n"
        "IMG Bin ID\tBin Quality\tBin Lineage\tGTDB-TK lineage\tBin Completeness\tBin Contamintation\tTotal Number of Bases\tNumber of genes\tNum of 5s rRNA\tNum of 16s rRNA\tNum of 23s rRNA\tNum of tRNA\tIMG scaffold IDs of members\n"
        f"3300000099_1\tHQ\tBacteria;Candidatus Gracilibacteria\td__Bacteria;p__Patescibacteriota;c__Gracilibacteria;o__Absconditabacterales;f__F;g__G;s__\t93.0\t0.8\t{len(s_cpr_25) + len(s_later_held)}\t120\t1\t1\t1\t24\ts_cpr,s_rmd\n"
    ).encode("utf-8")
    umag_zip = tmp_path / "umag.zip"
    _nested_jgi_zip(
        umag_zip,
        "CusIMGUMAG/IMG_Data/IMG_prokaryote_binning_files.tar",
        [
            (
                "./3300000099_1.tar.gz",
                _inner_tgz({
                    "3300000099_1.fna": _fa_bytes([("s_cpr", s_cpr_25), ("s_rmd", s_later_held)]),
                    "mbin_datafile_3300000099.txt": mbin99,
                }),
            ),
        ],
    )

    # Root stage P2 holds s_arch_shared; r26.11.4 stage V2 holds nothing; r26.11.5 stage RMD1 holds s_later_held
    root = tmp_path / "root"
    p2_dir = root / "l0seq" / "contig_index" / "stage=P2"
    p2_dir.mkdir(parents=True)
    pq.write_table(
        pa.table({
            "contig_uid": pa.array([bytes.fromhex(cuid(s_arch_shared))], type=pa.binary(16)),
            "is_default_variant": pa.array([True]),
        }),
        p2_dir / "part-0.parquet",
    )
    v2_dir = tmp_path / "chain" / "V2"
    v2_dir.mkdir(parents=True)
    (v2_dir / "genomes.json").write_text("[]")
    rmd1_dir = tmp_path / "later" / "RMD1"
    rmd1_dir.mkdir(parents=True)
    (rmd1_dir / "genomes.json").write_text(
        json.dumps([{"contig_uids": [cuid(s_later_held)], "genetic_code": 11, "kind": "prokaryote", "circular": [], "assembly_quality": "mag_hq"}])
    )

    def _run_once(tag: str, workers: int) -> dict:
        sub = tmp_path / tag
        cdir = sub / "census"
        mdir = sub / "meta"
        fdir = sub / "fasta"
        m.census_and_select(
            iuig_zip=iuig_zip,
            umag_zip=umag_zip,
            meta_zip=meta_zip,
            r232_index_dir=idx,
            gtdb_meta_dir=tmp_path / "GTDB_R232_META",
            gems_meta_tar=None,
            img_m_tar=None,
            out_dir=cdir,
            max_stage_bp=10 * L,
        )
        rep_stage = m.stage_jgiug1(
            census_dir=cdir,
            iuig_zip=iuig_zip,
            umag_zip=umag_zip,
            meta_out=mdir,
            fasta_out=fdir,
            root=root,
            root_stages=["P2"],
            owner_inputs={"V2": v2_dir},
            later_inputs={"RMD1": rmd1_dir},
            inputs_home=Path("/tmp/if2453/s2/inputs"),
            workers=workers,
            max_tasks_per_child=2,
        )
        return {"census": cdir, "meta": mdir / "JGIUG1", "fasta": fdir / "JGIUG1", "rep": rep_stage}

    r1 = _run_once("w1", workers=1)
    r2 = _run_once("w2", workers=2)

    rep = r1["rep"]
    assert rep["s2_check"]["exit"] == 0 and rep["s2_check"]["problems"] == [], rep["s2_check"]
    # 290000003 dropped to left_out.tsv for gap bytes; 290000002 is DEC-27 repeat of 290000001 -> 2 distinct genomes!
    assert rep["genomes"] == 2
    assert rep["genomes_repeated"] == 1
    assert rep["non_iupac_dropped_genomes"] == 1
    assert rep["codes"] == {"11": 1, "25": 1}
    assert rep["held_by"] == {"P2": 1}
    assert rep["later_held_by"] == {"RMD1": 1}
    assert rep["contigs_distinct"] == 4
    assert rep["contigs_new"] == 3  # after root + r26.11.4
    assert rep["contigs_new_after_later"] == 2  # after RMD1 too

    left_out = _read_tsv(r1["meta"] / "left_out.tsv")
    assert [x["genome_id"] for x in left_out] == ["290000003"]
    assert "gap" in left_out[0]["reason"] or "iupac" in left_out[0]["reason"]

    repeated = _read_tsv(r1["meta"] / "repeated_records.tsv")
    assert len(repeated) == 1
    assert repeated[0]["genome"] == "290000002"
    assert repeated[0]["first_genome"] == "290000001"
    assert repeated[0]["genome_uid"] == guid([s_arch_shared, s_arch_new])

    objs = json.loads((r1["meta"] / "objects.json").read_text())
    assert objs[0]["redistributable"] is False
    assert objs[0]["path"] == "/tmp/if2453/s2/inputs/JGIUG1/jgiug1_genomes.fna.gz"

    # Verify reconcile_stage_codes reads records.tsv + genetic_codes.tsv cleanly and passes s2_check
    rep_rec = m.reconcile_stage_codes(r1["meta"], r1["meta"] / "records.tsv")
    assert rep_rec["s2_check"]["exit"] == 0 and rep_rec["s2_check"]["problems"] == []
    assert rep_rec["sha256"]["genomes.json"] == r2["rep"]["sha256"]["genomes.json"]

    # Determinism across workers=1 and workers=2 (AGENTS.md §4.3)
    for fname in (
        "genomes.json",
        "objects.json",
        "records.tsv",
        "repeated_records.tsv",
        "shared_records.tsv",
        "later_records.tsv",
        "genetic_codes.tsv",
        "r232_status.tsv",
        "left_out.tsv",
    ):
        assert (r1["meta"] / fname).read_bytes() == (r2["meta"] / fname).read_bytes(), fname
    assert (r1["fasta"] / "jgiug1_genomes.fna.gz").read_bytes() == (r2["fasta"] / "jgiug1_genomes.fna.gz").read_bytes()


def test_bounded_memory_streaming_100x_synthetic(tmp_path: Path):
    """AGENTS.md §4.1 & §4.2: streaming 100 synthetic genomes keeps peak Python heap < 16 MiB."""
    m = _mod()
    direct_zip = tmp_path / "synth100.zip"
    want = set()
    with zipfile.ZipFile(direct_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for i in range(100):
            gid = f"SYNTH_{i:04d}"
            want.add(gid)
            seq = (("ACGT" * 250 + ("AAAA" if i % 2 else "TTTT") * 50) * 35).encode("ascii")
            zf.writestr(f"{gid}.fna.gz", _gz_bytes(_fa_bytes([(f"{gid}_c0", seq)])))

    tracemalloc.start()
    try:
        n_seen = 0
        total_bp = 0
        for _gid, recs in m.iter_selected_from_zip(direct_zip, want):
            n_seen += 1
            total_bp += sum(len(s) for _, s, _ in recs)
        _, peak_bytes = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    assert n_seen == 100
    assert total_bp > 4_000_000
    assert peak_bytes < 16 * 1024 * 1024, f"peak traced heap {peak_bytes} B exceeded 16 MiB bound"


def test_write_conf_template_passes_stage_est(tmp_path: Path):
    m = _mod()
    est = _load("stage_est", "stage_est.py")

    fake_tree = tmp_path / "tree"
    note_rel = "planning/measurements/2026-11-27_fixture_note.md"
    (fake_tree / "planning" / "measurements").mkdir(parents=True)
    (fake_tree / "planning").mkdir(parents=True, exist_ok=True)
    (fake_tree / "planning" / "DECISIONS.md").write_text("| DEC-88 | 2026-11-27 | ok |\n")
    sums = {
        "genomes.json": "a" * 64,
        "objects.json": "b" * 64,
        "jgiug1_genomes.fna.gz": "c" * 64,
    }
    (fake_tree / note_rel).write_text("\n".join(sums.values()) + "\n")
    tpl = tmp_path / "JGIUG1.conf"
    m.write_conf_template(
        tpl,
        stage="JGIUG1",
        note=note_rel,
        after="GBRMGD",
        sums=sums,
        expect_held_by={"P2": 3, "GART1": 12},
        summary_comment="Fixture JGIUG1 template",
    )

    conf = est.read_conf(tpl, "JGIUG1", tree=fake_tree)
    assert (conf.depth, conf.crispr, conf.library, conf.after, conf.genecall_mode, conf.min_contig_len) == (
        "called",
        True,
        None,
        ("GBRMGD",),
        "default",
        None,
    )
    assert conf.expect == {"P2": (3, note_rel), "GART1": (12, note_rel)}

    proc = subprocess.run(
        [sys.executable, "-I", str(BSP / "stage_est.py"), "conf", str(tpl), "--stage", "JGIUG1", "--tree", str(fake_tree)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert "STAGE=JGIUG1" in proc.stdout and "DEPTH=called" in proc.stdout and "CRISPR=yes" in proc.stdout


def test_committed_jgiug1_template_passes_stage_est():
    est = _load("stage_est", "stage_est.py")
    a26 = _load("a26_r26_11_4_chain", "a26_r26_11_4_chain.py")
    a29 = _load("a29_r26_11_5_chain", "a29_r26_11_5_chain.py")
    allowed = set(a29.ROOT_STAGES) | {s for s, _ in a26.R26_11_4_CHAIN} | set(a29.R26_11_5_CHAIN) | {"RMD1", "GBRMGD"}
    note_text = (REPO / NOTE).read_text(encoding="utf-8")
    tpl = TEMPLATES / "JGIUG1.conf"
    c = est.read_conf(tpl, "JGIUG1", tree=REPO)
    assert (c.depth, c.crispr, c.library, c.genecall_mode, c.min_contig_len) == ("called", True, None, "default", None)
    assert str(c.inputs) == "/tmp/if2453/s2/inputs/JGIUG1"
    assert NOTE in c.notes
    assert {name for _, name in c.sums} == {"genomes.json", "objects.json", "jgiug1_genomes.fna.gz"}
    for h, _ in c.sums:
        assert h in note_text, f"JGIUG1: SUM {h} missing from {NOTE}"
    assert set(c.expect) <= allowed
    for _owner, (k, _) in c.expect.items():
        assert k > 0


def test_transitive_shared_contig_code_component_resolution():
    """Regression test: a 3-hop shared-contig chain C(11) - B(25) - A(25) with guid(C) < guid(A) < guid(B)
    must assign a uniform genetic_code across the entire connected component (DEC-27).
    """
    m = _mod()
    u_ab = b"\x01" * 16
    u_bc = b"\x02" * 16
    u_a_only = b"\x03" * 16
    u_c_only = b"\x04" * 16
    contigs = {
        "gA": [(u_ab, 5000), (u_a_only, 5000)],
        "gB": [(u_ab, 5000), (u_bc, 5000)],
        "gC": [(u_bc, 5000), (u_c_only, 5000)],
    }
    # Force guid(gC) < guid(gA) < guid(gB) so a naive 1-pass greedy loop visits gC(11), then gA(25), then gB(25->11),
    # leaving u_ab in gA=25 vs gB=11!
    guid_map = {"gC": b"\x00" * 16, "gA": b"\x10" * 16, "gB": b"\x20" * 16}
    initial = {"gA": 25, "gB": 25, "gC": 11}
    code_map, conflicts = m._resolve_component_codes(
        ["gA", "gB", "gC"],
        lambda g: contigs[g],
        guid_map,
        initial,
    )
    assert code_map == {"gA": 11, "gB": 11, "gC": 11}
    assert len(conflicts) == 2
