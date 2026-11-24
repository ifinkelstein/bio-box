"""Ring R1 tests for card C-66: post-GlobDB r226 long-read (LR2) and niche/holobiont (NX1)
RT-hotspot compendia (DEC-69, DEC-71, DEC-65, DEC-59, DEC-50, DEC-51, DEC-26, DEC-27).
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import random
import subprocess
import sys
import tarfile
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-66")

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"


def _mod():
    spec = importlib.util.spec_from_file_location("c66_rt_hotspot_compendia", BSP / "c66_rt_hotspot_compendia.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c66_rt_hotspot_compendia"] = mod
    spec.loader.exec_module(mod)
    return mod


def _make_seq(seed: int, length: int) -> bytes:
    rng = random.Random(seed)
    return bytes(rng.choice(b"ACGT") for _ in range(length))


def test_c66_emerge_v2_selects_novel_species_reps_and_excludes_gtdb_globdb(tmp_path):
    m = _mod()
    tsv = tmp_path / "metadata_MAGs_v2_EMERGE.tsv"
    tsv.write_text(
        "genome\tcompleteness\tcontamination\tgtdb_taxonomy\tis_representative\tspecies_cluster\n"
        "MAG_A1\t95.0\t1.2\td__Bacteria;p__Acidobacteriota;s__GCA_000111222.1\tTrue\tclust_1\n"
        "MAG_A2\t98.0\t0.5\td__Bacteria;p__Myxococcota;s__\tTrue\tclust_2\n"
        "MAG_A3\t85.0\t2.0\td__Bacteria;p__Myxococcota;s__\tFalse\tclust_2\n"
        "MAG_A4\t45.0\t1.0\td__Bacteria;p__Planctomycetota;s__\tTrue\tclust_3\n"
        "MAG_A5\t91.0\t12.0\td__Bacteria;p__Planctomycetota;s__\tTrue\tclust_4\n"
        "MAG_A6\t72.0\t3.0\td__Bacteria;p__Patescibacteriota;s__\tTrue\tclust_5\n"
    )
    selected, deferred = m.filter_emerge_v2_metadata(
        tsv,
        held_species={"s__GCA_000111222.1"},
        min_completeness=50.0,
        max_contamination=10.0,
    )
    assert [r["genome"] for r in selected] == ["MAG_A2", "MAG_A6"]
    reasons = {r["genome"]: r["reason"] for r in deferred}
    assert reasons["MAG_A1"] == "species_in_gtdb_or_globdb"
    assert reasons["MAG_A3"] == "not_representative"
    assert reasons["MAG_A4"] == "below_mimag_mq"
    assert reasons["MAG_A5"] == "below_mimag_mq"


def test_c66_parse_skani_hits_deduplicates_at_96_ani_and_50_af(tmp_path):
    m = _mod()
    tsv = tmp_path / "skani_hits.tsv"
    tsv.write_text(
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        "/ref/G1.fa\t/q/mag_dup.fa\t97.2\t65.0\t60.0\tG1\tmag_dup\n"
        "/ref/G2.fa\t/q/mag_low_ani.fa\t94.5\t80.0\t78.0\tG2\tmag_low_ani\n"
        "/ref/G3.fa\t/q/mag_low_af.fa\t98.1\t30.0\t42.0\tG3\tmag_low_af\n"
        "/ref/G4.fa\t/q/mag_edge.fa\t96.0\t50.0\t49.0\tG4\tmag_edge\n"
    )
    matched = m.parse_skani_matches(tsv, min_ani=96.0, min_af=50.0)
    assert set(matched.keys()) == {"mag_dup", "mag_edge"}
    assert matched["mag_dup"]["ani"] == pytest.approx(97.2)


def test_c66_diff_ehomd_and_ena_against_gtdb_and_globdb(tmp_path):
    m = _mod()
    ehomd_tsv = tmp_path / "GCA_ID_info.txt"
    ehomd_tsv.write_text(
        "GCA-ID\tHMT-ID\tGenus\tSpecies\tStrain\tContigs\tCombined_Size\tStatus\tSequence_Source\n"
        "GCA_000001111.1\tHMT-001\tStreptococcus\toralis\ts1\t1\t2000000\tComplete Genome\tftp://x\n"
        "GCA_000002222.2\tHMT-002\tTM7x\tsp.\ts2\t5\t700000\tContig\tftp://y\n"
        "GCA_000003333.1\tHMT-003\tCapnocytophaga\tsp.\tERR_MAG\t20\t2200000\tContig\tftp://z\n"
    )
    novel, skipped = m.filter_ehomd_novel(
        ehomd_tsv,
        gtdb_accs={"GCA_000001111"},
        globdb_accs={"GCA_000002222"},
    )
    assert [r["gca_id"] for r in novel] == ["GCA_000003333.1"]
    assert len(skipped) == 2


def test_c66_stage_lr2_trims_exact_dtr_filters_unbinned_5kb_and_passes_s2_check(tmp_path):
    m = _mod()
    # Construct a binned MAG (2 contigs) and an unbinned HiFi assembly (1 circular >=5kb with 40bp exact DTR,
    # 1 linear >=5kb, 1 short <5kb dropped, 1 contig already in the binned MAG dropped from unbinned).
    c_bin1 = _make_seq(101, 6_000)
    c_bin2 = _make_seq(102, 5_500)
    dtr = _make_seq(201, 40)
    c_circ_core = dtr + _make_seq(202, 6_000)
    c_circ_raw = c_circ_core + dtr  # exact 40 bp terminal repeat
    c_lin = _make_seq(203, 7_000)
    c_short = _make_seq(204, 3_000)

    mag_dir = tmp_path / "mags"
    mag_dir.mkdir()
    mag_fa = mag_dir / "rumen_mag1.fa"
    mag_fa.write_bytes(b">b1\n" + c_bin1 + b"\n>b2\n" + c_bin2 + b"\n")

    asm_fa = tmp_path / "hifi_asm.fa"
    asm_fa.write_bytes(
        b">contig_circ circular=true\n" + c_circ_raw + b"\n"
        b">contig_lin\n" + c_lin + b"\n"
        b">contig_short\n" + c_short + b"\n"
        b">contig_dup_bin\n" + c_bin1 + b"\n"
    )

    mag_entries = [
        {
            "genome_id": "bickhart:rumen_mag1",
            "source_name": "bickhart_2022_rumen_hifi",
            "path": str(mag_fa),
            "assembly_quality": "mag_hq",
            "taxonomy": "d__Bacteria;p__Bacillota",
            "genetic_code": 11,
        }
    ]
    unbinned_entries = [
        {
            "source_name": "bickhart_2022_rumen_hifi",
            "path": str(asm_fa),
        }
    ]

    meta_dir = tmp_path / "meta_LR2"
    fasta_dir = tmp_path / "stor_LR2"
    rep = m.stage_compendium(
        stage="LR2",
        mag_entries=mag_entries,
        unbinned_entries=unbinned_entries,
        out_meta_dir=meta_dir,
        out_fasta_dir=fasta_dir,
        final_fasta_dir=tmp_path / "final_LR2",
        fasta_name="lr2_compendia.fna.gz",
        source_meta={
            "source_name": "c66_lr2_compendia",
            "source_release": "2026-11-24",
            "license": "CC BY 4.0",
            "redistributable": True,
        },
        min_unbinned_len=5_000,
        screen_codes=False,
    )

    # 1 MAG + 2 unbinned >= 5 kb contigs (short dropped, dup_bin already in MAG skipped)
    assert rep["mag_genomes"] == 1
    assert rep["unbinned_contigs_kept"] == 2
    assert rep["unbinned_short_dropped"] == 1
    assert rep["unbinned_already_binned_skipped"] == 1
    assert rep["dtr_trimmed"] == 1
    assert rep["circular_contigs"] == 1

    # Verify s2_check_inputs.py passes on objects_staging.json
    proc = subprocess.run(
        [
            sys.executable,
            str(BSP / "s2_check_inputs.py"),
            "LR2",
            str(meta_dir / "genomes.json"),
            str(meta_dir / "objects_staging.json"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    chk = json.loads(proc.stdout)
    assert chk["problems"] == []
    assert chk["genomes"] == 3
    assert chk["circular_rows"] == 1


def test_c66_stage_nx1_handles_owner_overlap_and_passes_s2_check(tmp_path):
    m = _mod()
    from biobox import ids

    c1 = _make_seq(301, 8_000)
    c2 = _make_seq(302, 9_000)
    c_shared = _make_seq(303, 7_500)

    # Fake earlier stage holding c_shared
    owner_dir = tmp_path / "owner_S1"
    owner_dir.mkdir()
    (owner_dir / "genomes.json").write_text(
        json.dumps([
            {
                "contig_uids": [ids.contig_uid(c_shared).hex()],
                "genetic_code": 11,
                "kind": "prokaryote",
                "circular": [],
                "assembly_quality": "mag_hq",
            }
        ])
    )

    mag_fa = tmp_path / "lost_city_1.fa"
    mag_fa.write_bytes(b">c1\n" + c1 + b"\n>c2\n" + c2 + b"\n>cs\n" + c_shared + b"\n")

    meta_dir = tmp_path / "meta_NX1"
    fasta_dir = tmp_path / "stor_NX1"
    rep = m.stage_compendium(
        stage="NX1",
        mag_entries=[
            {
                "genome_id": "lost_city:mag1",
                "source_name": "lost_city_serpentinizing_mags",
                "path": str(mag_fa),
                "assembly_quality": "mag_hq",
                "taxonomy": "d__Archaea;p__Methanobacteriota",
                "genetic_code": 11,
            }
        ],
        unbinned_entries=[],
        out_meta_dir=meta_dir,
        out_fasta_dir=fasta_dir,
        final_fasta_dir=tmp_path / "final_NX1",
        fasta_name="nx1_compendia.fna.gz",
        source_meta={
            "source_name": "c66_nx1_compendia",
            "source_release": "2026-11-24",
            "license": "CC BY 4.0 / INSDC",
            "redistributable": True,
        },
        owner_inputs={"S1": owner_dir},
        screen_codes=False,
    )

    assert rep["distinct_contigs"] == 3
    assert rep["shared_contigs"] == 1
    assert rep["net_new_contigs"] == 2

    proc = subprocess.run(
        [
            sys.executable,
            str(BSP / "s2_check_inputs.py"),
            "NX1",
            str(meta_dir / "genomes.json"),
            str(meta_dir / "objects_staging.json"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert json.loads(proc.stdout)["problems"] == []
