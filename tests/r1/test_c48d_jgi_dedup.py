"""Card C-48d: unit tests for streaming deduplication and staging of GEMs, IMG-UMAG
(prokaryote & eukaryote MAGs), and IUIG after GTDB R226 and GlobDB r226.

All fixtures are synthetic and built in `tmp_path` (no top-level repo files opened,
satisfying `tests/r0/test_image_files.py`).
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import io
import json
import resource
import shutil
import sys
import tarfile
import tracemalloc
import zipfile
from pathlib import Path

import pytest

from biobox import ids

pytestmark = pytest.mark.covers("C-48d")  # ring R1 by directory (tests/conftest.py)

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = ROOT / "scripts" / "stress" / "build_stage_pod"
META_SCRIPT = SCRIPTS_DIR / "c48d_metadata.py"
STREAM_SCRIPT = SCRIPTS_DIR / "c48d_stream.py"
STAGE_SCRIPT = SCRIPTS_DIR / "c48d_skani_stage.py"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None, f"missing script: {path}"
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _make_fa_bytes(records: list[tuple[str, str]]) -> bytes:
    lines = []
    for name, seq in records:
        lines.append(f">{name}\n{seq}\n")
    return "".join(lines).encode("ascii")


def _make_inner_tgz(files: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz", format=tarfile.GNU_FORMAT) as tf:
        for fname, data in files.items():
            ti = tarfile.TarInfo(name=fname)
            ti.size = len(data)
            ti.mtime = 0
            tf.addfile(ti, io.BytesIO(data))
    return buf.getvalue()


def _make_nested_zip(zip_path: Path, inner_tar_name: str, members: list[tuple[str, bytes]]) -> None:
    tar_buf = io.BytesIO()
    with tarfile.open(fileobj=tar_buf, mode="w", format=tarfile.GNU_FORMAT) as tf:
        for mname, data in members:
            ti = tarfile.TarInfo(name=mname)
            ti.size = len(data)
            ti.mtime = 0
            tf.addfile(ti, io.BytesIO(data))
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(inner_tar_name, tar_buf.getvalue())


def _make_gems_meta_tar(tar_path: Path) -> None:
    genome_meta = "\n".join([
        "genome_id\tmetagenome_id\tgenome_length\tnum_contigs\tn50\tnum_16s\tnum_5s\tnum_23s\tnum_trna\tcompleteness\tcontamination\tquality_score\tmimag_quality\totu_id\tecosystem\tecosystem_category\tecosystem_type\thabitat\tlongitude\tlatitude",
        "3300000001_1\t3300000001\t2500000\t40\t80000\t1\t1\t1\t20\t95.0\t1.0\t90.0\tHQ\tOTU-0001\tEnvironmental\tAquatic\tMarine\tOcean\t0\t0",
        "3300000001_2\t3300000001\t2100000\t60\t45000\t0\t1\t0\t15\t75.0\t2.0\t65.0\tMQ\tOTU-0002\tEnvironmental\tTerrestrial\tSoil\tForest\t0\t0",
        "3300000002_3\t3300000002\t1900000\t90\t22000\t0\t0\t0\t10\t45.0\t1.0\t40.0\tLQ\tOTU-0003\tEnvironmental\tTerrestrial\tSoil\tGrass\t0\t0",
        "3300000002_4\t3300000002\t3100000\t50\t70000\t1\t1\t1\t22\t88.0\t12.5\t25.5\tMQ\tOTU-0004\tHost-associated\tMammals\tDigestive system\tGut\t0\t0",
        "3300000003_5\t3300000003\t2800000\t35\t95000\t1\t1\t1\t24\t92.0\t0.5\t89.5\tHQ\tOTU-0005\tEngineered\tWastewater\tActivated sludge\tSludge\t0\t0",
    ]) + "\n"
    otu_reps = "\n".join([
        "otu_id\tgenome_id",
        "OTU-0001\t3300000001_1",
        "OTU-0002\t3300000001_2",
        "OTU-0003\t3300000002_3",
        "OTU-0004\t3300000002_4",
        "OTU-0005\t3300000003_5",
        "OTU-0006\tGCA_000001005.1",
    ]) + "\n"
    otu_tax = "\n".join([
        "otu_id\ttaxonomy\trank\tt_domain\tt_phylum\tt_class\tt_order\tt_family\tt_genus\tt_species",
        "OTU-0001\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__X;f__Y;g__Z;s__\tgenus\tBacteria\tPseudomonadota\tGammaproteobacteria\tX\tY\tZ\t",
        "OTU-0002\td__Bacteria;p__Acidobacteriota;c__A;o__B;f__C;g__D;s__\tgenus\tBacteria\tAcidobacteriota\tA\tB\tC\tD\t",
        "OTU-0003\td__Bacteria;p__Bacillota;c__A;o__B;f__C;g__D;s__\tgenus\tBacteria\tBacillota\tA\tB\tC\tD\t",
        "OTU-0004\td__Bacteria;p__Bacteroidota;c__A;o__B;f__C;g__D;s__\tgenus\tBacteria\tBacteroidota\tA\tB\tC\tD\t",
        "OTU-0005\td__Archaea;p__Thermoproteota;c__A;o__B;f__C;g__D;s__\tgenus\tArchaea\tThermoproteota\tA\tB\tC\tD\t",
    ]) + "\n"
    virsorter = "\n".join([
        "bin_id\tprophage_id\tcontig_id\tcontig_start\tcontig_stop",
        "3300000001_1\tVP_1\tGa0001_1001\t100\t35100",
        "3300000001_2\tVP_2\tGa0001_2001\t500\t45500",
        "3300000003_5\tVP_3\tGa0003_5001\t1000\t26000",
    ]) + "\n"
    with tarfile.open(tar_path, "w", format=tarfile.GNU_FORMAT) as tf:
        for name, text in [
            ("gems/genomes/genome_metadata.tsv", genome_meta),
            ("gems/otus/otu_representative_genomes.tsv", otu_reps),
            ("gems/otus/otu_taxonomy.tsv", otu_tax),
            ("gems/prophages/virsorter.tsv", virsorter),
        ]:
            raw = text.encode("utf-8")
            ti = tarfile.TarInfo(name=name)
            ti.size = len(raw)
            ti.mtime = 0
            tf.addfile(ti, io.BytesIO(raw))


def test_metadata_gems_and_virsorter_dedup_against_globdb_and_gtdb(tmp_path: Path):
    m = _load("c48d_metadata", META_SCRIPT)
    gems_tar = tmp_path / "gems_meta.tar"
    _make_gems_meta_tar(gems_tar)

    # GlobDB already holds OTU-0001 (3300000001_1) as a GEMOTU rep; OTU-0003 is LQ (<50% comp);
    # OTU-0004 has >10% contam; so only OTU-0002 and OTU-0005 survive as HQ/MQ novel GEM candidates.
    globdb_refs = {
        "globdb_gem_ids": {"3300000001_1"},
        "globdb_all_ids": {"3300000001_1", "GEMOTU_00001"},
        "gtdb_accessions": {"GCA_000001005.1"},
        "gtdb_stems": {"GCA_000001005"},
    }
    res = m.analyze_gems_metadata(gems_tar, globdb_refs)
    assert res["total_mags"] == 5
    assert res["hq_mq_mags"] == 3  # 1, 2, 5
    assert res["already_in_globdb_or_gtdb"] == 1  # 3300000001_1
    cand_ids = [r["genome_id"] for r in res["candidates"]]
    # Sorted by (taxonomy, genome_id): d__Archaea (3300000003_5) before d__Bacteria (3300000001_2)
    assert cand_ids == ["3300000003_5", "3300000001_2"]
    # VirSorter prophage summary
    vp = res["virsorter"]
    assert vp["total_prophages"] == 3
    assert vp["unique_bins"] == 3
    assert vp["on_globdb_bins"] == 1
    assert vp["on_novel_candidate_bins"] == 2
    assert vp["median_len_bp"] == 35000


def test_metadata_iuig_and_umag_filtering(tmp_path: Path):
    m = _load("c48d_metadata", META_SCRIPT)
    meta_zip = tmp_path / "umag_iuig_meta.zip"
    iuig_txt = "\n".join([
        "taxon_oid\ttaxon_display_name\tdomain\tphylum\tir_class\tir_order\tfamily\tgenus\tspecies\tseq_status\tseq_center\tanalysis_project_id\tncbi_taxon_id\tis_public\tadd_date",
        "637000001\tEscherichia coli K12\tBacteria\tPseudomonadota\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae\tEscherichia\tEscherichia coli\tFinished\tJGI\tGa0000001\t511145\tYes\t01-JAN-06",
        "637000002\tNovel thermophile sp1\tArchaea\tThermoproteota\tNitrososphaeria\tNitrososphaerales\tNitrososphaeraceae\tNovel\tNovel sp1\tPermanent Draft\tJGI\tGa0000002\t9999991\tYes\t01-JAN-22",
        "637000003\tSaccharomyces cerevisiae\tEukaryota\tAscomycota\tSaccharomycetes\tSaccharomycetales\tSaccharomycetaceae\tSaccharomyces\tSaccharomyces cerevisiae\tFinished\tJGI\tGa0000003\t4932\tYes\t01-JAN-07",
    ]) + "\n"
    umag_txt = "\n".join([
        "taxon_oid\ttaxon_display_name\tdomain\tphylum\tir_class\tir_order\tfamily\tgenus\tspecies\tseq_status\tseq_center\tanalysis_project_id\tncbi_taxon_id\tis_public\tadd_date",
        "3300000001\tMarine sample 1\t*Microbiome\tEnvironmental\tAquatic\tMarine\tOceanic\tUnclassified\tMarine\tPermanent Draft\tJGI\tGa0100001\t\tYes\t01-JAN-20",
        "3300000009\tSoil sample 2\t*Microbiome\tEnvironmental\tTerrestrial\tSoil\tForest\tUnclassified\tSoil\tPermanent Draft\tJGI\tGa0100002\t\tYes\t01-JAN-21",
    ]) + "\n"
    with zipfile.ZipFile(meta_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("Custom_MPIIUIG/2023-11-14_1.0/img_genomes_tarballs_metadata_111023.txt", iuig_txt)
        zf.writestr("CusIMGUMAG/IMG_Data/unrestricted_metagenomes_for_bin_tarballs.02062024.txt", umag_txt)

    # Mark 637000001 as already in GTDB via ncbi_taxon_id or IMG->GCA mapping
    refs = {
        "gtdb_ncbi_taxids": {"511145"},
        "gtdb_species_names": {"Escherichia coli"},
        "img_oid_to_gca": {"637000001": "GCF_000005845.2"},
        "gtdb_stems": {"GCF_000005845"},
        "gems_metagenome_ids": {"3300000001"},
    }
    res = m.analyze_iuig_umag_metadata(meta_zip, refs)
    assert res["iuig"]["total_rows"] == 3
    assert res["iuig"]["by_domain"]["Bacteria"] == 1
    assert res["iuig"]["by_domain"]["Archaea"] == 1
    assert res["iuig"]["by_domain"]["Eukaryota"] == 1
    assert res["iuig"]["prok_matched_gtdb_or_globdb"] == 1
    assert res["iuig"]["prok_novel_candidate_oids"] == ["637000002"]
    assert res["umag"]["total_metagenomes"] == 2
    assert res["umag"]["shared_with_gems_metagenomes"] == 1


def test_stream_gems_and_nested_bins_filter_hash_and_overlap(tmp_path: Path):
    s = _load("c48d_stream", STREAM_SCRIPT)

    seq_shared = "ATG" * 400
    seq_novel_1 = "ACGT" * 350
    seq_novel_2 = "GCTA" * 300
    shared_uid = ids.contig_uid(seq_shared.encode("ascii"))
    known_uids = {shared_uid: "P1"}

    # 1) GEMs fna tar: 2 members, only one in candidate_ids; the skipped member has corrupted gzip
    # to prove the streamer does NOT decompress skipped members.
    gems_fna_tar = tmp_path / "gems_fna.tar"
    good_gz = gzip.compress(_make_fa_bytes([("c1", seq_shared), ("c2", seq_novel_1)]))
    bad_gz = b"NOT_VALID_GZIP_MUST_NOT_BE_READ"
    with tarfile.open(gems_fna_tar, "w", format=tarfile.GNU_FORMAT) as tf:
        for name, data in [("fna/3300000001_1.fna.gz", bad_gz), ("fna/3300000001_2.fna.gz", good_gz)]:
            ti = tarfile.TarInfo(name=name)
            ti.size = len(data)
            ti.mtime = 0
            tf.addfile(ti, io.BytesIO(data))

    out_extract = tmp_path / "extract_gems"
    g_res = s.stream_gems_fna_tar(
        gems_fna_tar,
        candidate_ids={"3300000001_2"},
        known_contig_uids=known_uids,
        extract_dir=out_extract,
    )
    assert g_res["members_seen"] == 2
    assert g_res["members_skipped"] == 1
    assert g_res["members_processed"] == 1
    assert g_res["shared_contigs"] == 1
    assert g_res["shared_by_stage"] == {"P1": 1}
    assert g_res["new_contigs"] == 1
    assert (out_extract / "3300000001_2.fna.gz").is_file()

    # 2) Nested IMG-UMAG eukaryote zip
    euk_zip = tmp_path / "euk.zip"
    gff_text = "\n".join([
        "scaf1\tProdigal v2.6.3\tCDS\t1\t300\t10.0\t+\t0\tID=g1;translation_table=11",
        "scaf1\tGeneMark.hmm-2 v1.05\tCDS\t400\t760\t12.0\t-\t0\tID=g2;translation_table=11",
    ]) + "\n"
    euk_tgz = _make_inner_tgz({
        "3300037603_4.fna": _make_fa_bytes([("scaf1", seq_novel_2)]),
        "3300037603_4.gff": gff_text.encode("utf-8"),
    })
    _make_nested_zip(euk_zip, "CusIMGUMAG/IMG_Data/IMG_eukaryote_binning_files.tar", [("./3300037603_4.tar.gz", euk_tgz)])

    e_res = s.stream_nested_bin_zip(
        euk_zip,
        inner_tar_name="CusIMGUMAG/IMG_Data/IMG_eukaryote_binning_files.tar",
        mode="umag_euk",
        candidate_ids=None,
        known_contig_uids=known_uids,
    )
    assert e_res["members_processed"] == 1
    row = e_res["records"][0]
    assert row["bin_id"] == "3300037603_4"
    assert row["bp"] == len(seq_novel_2)
    assert row["gff_cds_count"] == 2
    assert "Prodigal v2.6.3" in row["gff_tools"]

    # 3) Nested IMG-UMAG prokaryote zip with HQ/MQ filtering from mbin_datafile
    prok_zip = tmp_path / "prok.zip"
    mbin_txt = "\n".join([
        "# IMG TAXON ID: 3300006461",
        "IMG Bin ID\tBin Quality\tBin Lineage\tGTDB-TK lineage\tBin Completeness\tBin Contamintation\tTotal Number of Bases\tNumber of genes\tNum of 5s rRNA\tNum of 16s rRNA\tNum of 23s rRNA\tNum of tRNA\tIMG scaffold IDs of members",
        "3300006461_13\tHQ\tBacteria;Firmicutes\tBacteria;Bacillota;Bacilli;O;F;G;S\t94.5\t1.2\t1200\t10\t1\t1\t1\t20\ts1",
        "3300006461_99\tLQ\tBacteria;Firmicutes\tBacteria;Bacillota;Bacilli;O;F;G;S\t35.0\t1.0\t1200\t10\t0\t0\t0\t5\ts1",
    ]) + "\n"
    hq_tgz = _make_inner_tgz({
        "3300006461_13.fna": _make_fa_bytes([("s1", seq_novel_1)]),
        "mbin_datafile_3300006461.txt": mbin_txt.encode("utf-8"),
    })
    lq_tgz = _make_inner_tgz({
        "3300006461_99.fna": _make_fa_bytes([("s1", seq_novel_2)]),
        "mbin_datafile_3300006461.txt": mbin_txt.encode("utf-8"),
    })
    _make_nested_zip(
        prok_zip,
        "CusIMGUMAG/IMG_Data/IMG_prokaryote_binning_files.tar",
        [("./3300006461_13.tar.gz", hq_tgz), ("./3300006461_99.tar.gz", lq_tgz)],
    )
    p_res = s.stream_nested_bin_zip(
        prok_zip,
        inner_tar_name="CusIMGUMAG/IMG_Data/IMG_prokaryote_binning_files.tar",
        mode="umag_prok",
        candidate_ids=None,
        known_contig_uids=known_uids,
        filter_hq_mq=True,
    )
    assert p_res["members_seen"] == 2
    assert p_res["hq_mq_kept"] == 1
    assert p_res["lq_filtered"] == 1
    assert p_res["records"][0]["bin_id"] == "3300006461_13"
    assert p_res["records"][0]["completeness"] == 94.5


def test_streaming_bounded_rss_on_100x_synthetic_archive(tmp_path: Path):
    """AGENTS.md §4.1 & §4.2: streaming 100 synthetic MAGs keeps peak Python heap bounded."""
    s = _load("c48d_stream", STREAM_SCRIPT)
    tar_path = tmp_path / "synth_100x.tar"
    # 100 synthetic MAGs x ~50 kbp = ~5 Mbp total; peak traced allocation must stay well below total
    with tarfile.open(tar_path, "w", format=tarfile.GNU_FORMAT) as tf:
        for i in range(100):
            seq = ("ACGT" * 250 + f"{'AAAA' if i % 2 else 'TTTT'}" * 50) * 40  # ~48 kbp
            gz = gzip.compress(_make_fa_bytes([(f"c_{i}_0", seq)]))
            ti = tarfile.TarInfo(name=f"fna/SYNTH_{i:04d}.fna.gz")
            ti.size = len(gz)
            ti.mtime = 0
            tf.addfile(ti, io.BytesIO(gz))

    tracemalloc.start()
    try:
        res = s.stream_gems_fna_tar(tar_path, candidate_ids=None, known_contig_uids={})
        _, peak_bytes = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    assert res["members_processed"] == 100
    assert res["new_contigs"] == 100
    # Traced Python heap peak must stay < 16 MiB (strictly bounded, never holding all 100 genomes)
    assert peak_bytes < 16 * 1024 * 1024, f"peak traced heap {peak_bytes} B too high"
    assert rss_mb > 0


def test_skani_species_rep_selection_and_stage_proposals(tmp_path: Path):
    st = _load("c48d_skani_stage", STAGE_SCRIPT)

    candidates = [
        {"genome_id": "G_A", "source": "gems", "completeness": 96.0, "contamination": 1.0, "bp": 4_000_000, "n50": 90_000, "taxonomy": "d__Bacteria;p__P1;c__C1;o__O1;f__F1;g__Gen1;s__", "ecosystem": "Aquatic"},
        {"genome_id": "G_B", "source": "gems", "completeness": 80.0, "contamination": 2.0, "bp": 3_800_000, "n50": 40_000, "taxonomy": "d__Bacteria;p__P1;c__C1;o__O1;f__F1;g__Gen1;s__", "ecosystem": "Aquatic"},
        {"genome_id": "G_C", "source": "umag_prok", "completeness": 91.0, "contamination": 0.5, "bp": 5_000_000, "n50": 75_000, "taxonomy": "d__Bacteria;p__P2;c__C2;o__O2;f__F2;g__Gen2;s__", "ecosystem": "Soil"},
        {"genome_id": "G_D", "source": "iuig", "completeness": 99.0, "contamination": 0.1, "bp": 4_500_000, "n50": 120_000, "taxonomy": "d__Bacteria;p__P3;c__C3;o__O3;f__F3;g__Gen3;s__", "ecosystem": "Marine"},
    ]
    # G_D matches reference REF_GTDB_1 at 97.2% ANI / 80% AF -> dropped as redundant with GTDB/GlobDB.
    # G_A and G_B match each other at 96.1% ANI / 70% AF -> G_A wins (higher quality score).
    # G_C has no hit >= 95% ANI -> survives as novel species rep.
    edges = [
        ("REF_GTDB_1", "G_D", 97.2, 80.0, 82.0),
        ("G_A", "G_B", 96.1, 70.0, 72.0),
        ("G_A", "G_C", 84.0, 60.0, 60.0),
    ]
    sel = st.select_species_reps(
        candidates,
        skani_edges=edges,
        ani_threshold=95.0,
        min_af=50.0,
        ref_ids={"REF_GTDB_1"},
    )
    assert [r["genome_id"] for r in sel["novel_reps"]] == ["G_A", "G_C"]
    assert sel["dropped_by_ref"] == ["G_D"]
    assert sel["clustered_into_novel_rep"] == {"G_B": "G_A"}

    # Stage proposal: prokaryote ~10 Gbp (test with 6 Mbp target to force 2 stages)
    prok_stages = st.propose_stages(sel["novel_reps"], prefix="J", target_bp=6_000_000, redistributable=False)
    assert len(prok_stages) == 2
    assert all(stage["redistributable"] is False for stage in prok_stages)
    assert prok_stages[0]["stage"] == "J1"
    assert prok_stages[1]["stage"] == "J2"

    # Eukaryote stage proposal: separate E1.. stages with explicit intron-aware gene-calling note
    euk_bins = [
        {"genome_id": "3300037603_4", "source": "umag_euk", "bp": 15_000_000, "taxonomy": "d__Eukaryota;p__Fungi", "ecosystem": "Soil"},
        {"genome_id": "3300055075_55", "source": "umag_euk", "bp": 12_000_000, "taxonomy": "d__Eukaryota;p__Chlorophyta", "ecosystem": "Marine"},
    ]
    euk_stages = st.propose_stages(
        euk_bins,
        prefix="E",
        target_bp=20_000_000,
        redistributable=False,
        kind="eukaryote",
        depth_recommendation="stored",
        gene_caller_note="intron-bearing eukaryote bins require eukaryotic gene caller (MetaEuk/AUGUSTUS/Helixer); do not run pyrodigal",
    )
    assert euk_stages[0]["stage"] == "E1"
    assert euk_stages[0]["kind"] == "eukaryote"
    assert euk_stages[0]["depth_recommendation"] == "stored"
    assert "intron" in euk_stages[0]["gene_caller_note"].lower()

    # Locked pilot input manifest (deterministic, redistributable=False)
    fa_dir = tmp_path / "pilot_fa"
    fa_dir.mkdir()
    p_a = fa_dir / "G_A.fna.gz"
    p_c = fa_dir / "G_C.fna.gz"
    p_a.write_bytes(gzip.compress(_make_fa_bytes([("a1", "ACGT" * 200)])))
    p_c.write_bytes(gzip.compress(_make_fa_bytes([("c1", "TGCA" * 200)])))

    out1 = tmp_path / "pilot1"
    out2 = tmp_path / "pilot2"
    m1 = st.write_pilot_manifest(
        sel["novel_reps"],
        {"G_A": p_a, "G_C": p_c},
        out1,
        source_name="jgi_gems_umag_iuig",
        source_release="2026-11-18",
        license_str="JGI Data Utilization Policy",
        redistributable=False,
    )
    m2 = st.write_pilot_manifest(
        list(reversed(sel["novel_reps"])),
        {"G_A": p_a, "G_C": p_c},
        out2,
        source_name="jgi_gems_umag_iuig",
        source_release="2026-11-18",
        license_str="JGI Data Utilization Policy",
        redistributable=False,
    )
    assert m1["genomes"] == 2
    assert (out1 / "genomes.json").read_bytes() == (out2 / "genomes.json").read_bytes()
    assert (out1 / "objects.json").read_bytes() == (out2 / "objects.json").read_bytes()
    assert (out1 / "inputs_SUMS").read_bytes() == (out2 / "inputs_SUMS").read_bytes()
    objs = json.loads((out1 / "objects.json").read_text())
    assert all(o["redistributable"] is False for o in objs)
    # Verify inputs_SUMS hashes match actual files
    for line in (out1 / "inputs_SUMS").read_text().splitlines():
        digest, fname = line.split()
        assert hashlib.sha256((out1 / fname).read_bytes()).hexdigest() == digest


@pytest.mark.skipif(shutil.which("skani") is None, reason="skani not on PATH")
def test_skani_real_cli_dedup(tmp_path: Path):
    """Oracle check against real skani binary inside the pinned container image."""
    st = _load("c48d_skani_stage", STAGE_SCRIPT)
    import random

    rng = random.Random(42)
    base_seq = "".join(rng.choices("ACGT", k=60_000))
    # Mutate 1% -> ~99% ANI copy
    seq_near = list(base_seq)
    for idx in range(0, len(seq_near), 100):
        seq_near[idx] = "T" if seq_near[idx] != "T" else "A"
    seq_near_s = "".join(seq_near)
    # Independent sequence -> <80% ANI
    seq_far = "".join(rng.choices("ACGT", k=60_000))

    f_ref = tmp_path / "REF1.fna.gz"
    f_dup = tmp_path / "CAND_DUP.fna.gz"
    f_nov = tmp_path / "CAND_NOV.fna.gz"
    f_ref.write_bytes(gzip.compress(_make_fa_bytes([("r1", base_seq)])))
    f_dup.write_bytes(gzip.compress(_make_fa_bytes([("d1", seq_near_s)])))
    f_nov.write_bytes(gzip.compress(_make_fa_bytes([("n1", seq_far)])))

    res = st.run_skani_dedup(
        candidate_fastas={"CAND_DUP": f_dup, "CAND_NOV": f_nov},
        ref_fastas={"REF1": f_ref},
        work_dir=tmp_path / "skani_work",
        threads=1,
        ani_threshold=95.0,
        min_af=50.0,
    )
    assert any(e[0] == "REF1" and e[1] == "CAND_DUP" and e[2] >= 95.0 for e in res["edges"])
    assert not any(e[1] == "CAND_NOV" and e[2] >= 95.0 for e in res["edges"])
