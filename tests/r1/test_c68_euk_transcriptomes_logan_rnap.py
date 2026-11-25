"""Ring R1 tests for card C-68: Spliced protist & plant transcriptomes (`ET1`: `MMETSP` + `1KP`)
and targeted `Logan` `ssRNAP` baiting (`LGNRNAP1`, W-40 `ssrnap` panel).
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
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-68")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
STAGES = BSP / "stages" / "templates"


def _mod():
    spec = importlib.util.spec_from_file_location(
        "c68_euk_transcriptomes_logan_rnap", BSP / "c68_euk_transcriptomes_logan_rnap.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c68_euk_transcriptomes_logan_rnap"] = mod
    spec.loader.exec_module(mod)
    return mod


def _stage_est():
    spec = importlib.util.spec_from_file_location("stage_est", BSP / "stage_est.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["stage_est"] = mod
    spec.loader.exec_module(mod)
    return mod


def _make_seq(seed: int, length: int) -> bytes:
    rng = random.Random(seed)
    return bytes(rng.choice(b"ACGT") for _ in range(length))


def test_c68_cdna_length_filter_rna_normalization_and_isoform_dedup():
    m = _mod()
    s_long_1 = _make_seq(1, 2_400)
    s_long_2 = _make_seq(2, 1_800)
    s_rna = s_long_2.replace(b"T", b"U")  # identical after U -> T normalization
    s_exact_1500 = _make_seq(3, 1_500)
    s_short_1499 = _make_seq(4, 1_499)

    records = [
        ("TRINITY_DN10_c0_g1_i1 len=2400", s_long_1),
        ("TRINITY_DN10_c0_g1_i2 len=1800", s_long_2),
        ("TRINITY_DN10_c0_g1_i3 len=1800", s_rna),
        ("scaffold-ABCD-100001", s_exact_1500),
        ("TRINITY_DN99_c0_g1_i1 len=1499", s_short_1499),
        ("BAD_SEQ", b"ACGT12345"),
    ]

    # Without isoform dedup: exact duplicate (s_rna == s_long_2) collapsed, < 1500 and invalid logged
    kept_all, left_all, st_all = m.filter_and_dedup_cdna_records(
        records,
        sample_id="MMETSP0001",
        source_name="mmetsp_reassemblies",
        organism="Naegleria gruberi",
        min_len=1_500,
        dedup_isoforms=False,
    )
    assert len(kept_all) == 3
    assert st_all["within_sample_exact_dup"] == 1
    assert st_all["below_min_len"] == 1
    assert st_all["below_min_len_bp"] == 1_499
    assert st_all["invalid_seq"] == 1
    reasons_all = {r["transcript_id"]: r["reason"] for r in left_all}
    assert reasons_all["TRINITY_DN99_c0_g1_i1"] == "below_min_len_1500nt"
    assert reasons_all["BAD_SEQ"] == "invalid_or_empty_seq"

    # With isoform dedup: only longest isoform (i1, 2400 nt) survives for locus TRINITY_DN10_c0_g1
    kept_iso, left_iso, st_iso = m.filter_and_dedup_cdna_records(
        records,
        sample_id="MMETSP0001",
        source_name="mmetsp_reassemblies",
        organism="Naegleria gruberi",
        min_len=1_500,
        dedup_isoforms=True,
    )
    assert [r["transcript_id"] for r in kept_iso] == ["TRINITY_DN10_c0_g1_i1", "scaffold-ABCD-100001"]
    assert st_iso["shorter_isoform_dropped"] == 1
    iso_reasons = {r["transcript_id"]: r["reason"] for r in left_iso}
    assert iso_reasons["TRINITY_DN10_c0_g1_i2"] == "shorter_trinity_isoform"


def test_c68_taxon_translation_table_routing_and_orfs_check_code_compat():
    m = _mod()
    from biobox.search import orfs

    cases = [
        ("Arabidopsis thaliana", "Eukaryota; Viridiplantae; Streptophyta", 1),
        ("Trypanosoma brucei", "Eukaryota; Discoba; Euglenozoa; Kinetoplastea", 1),
        ("Stentor coeruleus", "Eukaryota; Sar; Alveolata; Ciliophora; Heterotrichea", 1),
        ("Blepharisma japonicum", "Eukaryota; Sar; Alveolata; Ciliophora; Heterotrichea", 1),
        ("Tetrahymena thermophila", "Eukaryota; Sar; Alveolata; Ciliophora; Oligohymenophorea", 6),
        ("Paramecium tetraurelia", "Eukaryota; Sar; Alveolata; Ciliophora; Oligohymenophorea", 6),
        ("Oxytricha trifallax", "Eukaryota; Sar; Alveolata; Ciliophora; Spirotrichea", 6),
        ("Acetabularia acetabulum", "Eukaryota; Viridiplantae; Chlorophyta; Ulvophyceae; Dasycladales", 6),
        ("Hexamita inflata", "Eukaryota; Metamonada; Fornicata; Hexamitidae", 6),
        ("Euplotes crassus", "Eukaryota; Sar; Alveolata; Ciliophora; Spirotrichea; Euplotida", 10),
        ("Mesodinium pulex", "Eukaryota; Sar; Alveolata; Ciliophora; Litostomatea", 29),
        ("Vorticella campanula", "Eukaryota; Sar; Alveolata; Ciliophora; Oligohymenophorea; Peritrichia", 30),
    ]

    for org, tax, expected_bio in cases:
        res_defer = m.resolve_taxon_genetic_code(org, tax, unsupported_code_policy="defer")
        assert res_defer["bio_transl_table"] == expected_bio, f"{org}: got {res_defer}"
        if expected_bio in (1, 6, 10):
            assert res_defer["orfs_supported"] is True
            assert res_defer["action"] == "stage"
            assert res_defer["effective_transl_table"] == expected_bio
            orfs.check_code(res_defer["effective_transl_table"])
        else:
            # 29 and 30: either natively supported if orfs.py is later extended, or cleanly deferred
            if res_defer["orfs_supported"]:
                assert res_defer["action"] == "stage"
                orfs.check_code(res_defer["effective_transl_table"])
            else:
                assert res_defer["action"] == "defer"
                assert res_defer["effective_transl_table"] is None
                assert f"transl_table_{expected_bio}" in res_defer["defer_reason"]

        # fallback_table6 policy must ALWAYS yield a code accepted by orfs.check_code
        res_fb = m.resolve_taxon_genetic_code(org, tax, unsupported_code_policy="fallback_table6")
        assert res_fb["action"] == "stage"
        assert res_fb["effective_transl_table"] is not None
        orfs.check_code(res_fb["effective_transl_table"])
        if expected_bio in (29, 30) and not res_fb["orfs_supported"]:
            assert res_fb["effective_transl_table"] == 6


def test_c68_stage_et1_determinism_left_for_later_and_s2_check_inputs(tmp_path):
    m = _mod()
    from biobox.search import orfs

    s_plant_1 = _make_seq(201, 2_800)
    s_plant_short = _make_seq(202, 900)
    s_ciliate_1 = _make_seq(203, 3_100)
    s_euplotes_1 = _make_seq(204, 2_200)
    s_mesodinium_1 = _make_seq(205, 2_600)

    samples = [
        {
            "sample_id": "MMETSP0100",
            "source_name": "mmetsp_reassemblies",
            "source_release": "zenodo-3247846",
            "license": "CC BY 4.0",
            "organism": "Tetrahymena thermophila",
            "taxonomy": "Eukaryota; Sar; Alveolata; Ciliophora",
            "records": [
                ("TRINITY_DN1_c0_g1_i1", s_ciliate_1),
                # Same sequence as plant sample with different genetic_code -> cross-code conflict safely dropped
                ("TRINITY_DN2_c0_g1_i1", s_plant_1),
            ],
        },
        {
            "sample_id": "MMETSP0200",
            "source_name": "mmetsp_reassemblies",
            "source_release": "zenodo-3247846",
            "license": "CC BY 4.0",
            "organism": "Euplotes focardii",
            "taxonomy": "Eukaryota; Sar; Alveolata; Ciliophora; Euplotida",
            "records": [("TRINITY_DN5_c0_g1_i1", s_euplotes_1)],
        },
        {
            "sample_id": "MMETSP0300",
            "source_name": "mmetsp_reassemblies",
            "source_release": "zenodo-3247846",
            "license": "CC BY 4.0",
            "organism": "Mesodinium pulex",
            "taxonomy": "Eukaryota; Sar; Alveolata; Ciliophora",
            "records": [("TRINITY_DN9_c0_g1_i1", s_mesodinium_1)],
        },
        {
            "sample_id": "ONEKP_ABCD",
            "source_name": "onekp_transcriptomes",
            "source_release": "gigadb-100910",
            "license": "CC0 1.0",
            "organism": "Arabidopsis thaliana",
            "taxonomy": "Eukaryota; Viridiplantae",
            "records": [
                ("scaffold-ABCD-1", s_plant_1),
                ("scaffold-ABCD-2", s_plant_short),
            ],
        },
    ]

    # Run 1 (forward order) vs Run 2 (reversed sample & record order) -> byte-identical outputs
    run1_meta, run1_fa = tmp_path / "run1_meta", tmp_path / "run1_fa"
    run2_meta, run2_fa = tmp_path / "run2_meta", tmp_path / "run2_fa"

    rep1 = m.stage_et1(
        samples=samples,
        out_meta_dir=run1_meta,
        out_fasta_dir=run1_fa,
        min_len=1_500,
        unsupported_code_policy="defer",
    )
    rev_samples = [{**s, "records": list(reversed(s["records"]))} for s in reversed(samples)]
    rep2 = m.stage_et1(
        samples=rev_samples,
        out_meta_dir=run2_meta,
        out_fasta_dir=run2_fa,
        min_len=1_500,
        unsupported_code_policy="defer",
    )

    assert rep1["fasta_sha256"] == rep2["fasta_sha256"]
    assert rep1["genomes_json_sha256"] == rep2["genomes_json_sha256"]
    assert rep1["et1_metadata_tsv_sha256"] == rep2["et1_metadata_tsv_sha256"]
    assert rep1["left_for_later_tsv_gz_sha256"] == rep2["left_for_later_tsv_gz_sha256"]

    # Check ET1_left_for_later.tsv.gz records both <1500 nt transcript and deferred table-29 sample
    with gzip.open(run1_meta / "ET1_left_for_later.tsv.gz", "rt", encoding="utf-8") as fh:
        left_rows = list(csv.DictReader(fh, delimiter="\t"))
    reasons = {r["reason"] for r in left_rows}
    assert "below_min_len_1500nt" in reasons
    if not m.check_orfs_code_supported(29):
        assert "transl_table_29_unsupported_by_orfs_py" in reasons

    # Verify s2_check_inputs.py passes with 0 problems and every emitted code passes orfs.check_code
    proc = subprocess.run(
        [
            sys.executable,
            str(BSP / "s2_check_inputs.py"),
            "ET1",
            str(run1_meta / "genomes.json"),
            str(run1_meta / "objects_staging.json"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, f"s2_check_inputs failed: {proc.stdout}\n{proc.stderr}"
    chk = json.loads(proc.stdout)
    assert chk["problems"] == []
    for code in chk["codes"]:
        orfs.check_code(code)


def test_c68_load_ssrnap_cutoffs_and_stage_lgnrnap1(tmp_path):
    m = _mod()
    from biobox import ids

    cutoffs = m.load_ssrnap_cutoffs()
    assert "palm:mito_invertron_rpo" in cutoffs
    assert "basal_protist_rpot" in cutoffs
    assert "dpbb_single:crass_vrnap_dpbb" in cutoffs
    assert cutoffs["palm:mito_invertron_rpo"]["seq_bits"] == 30.0
    assert cutoffs["palm:mito_invertron_rpo"]["dom_bits"] == 20.0

    s_new1 = _make_seq(301, 3_500)
    s_new2 = _make_seq(302, 4_200)
    s_shared = _make_seq(303, 2_900)

    earlier_dir = tmp_path / "earlier_ET1"
    earlier_dir.mkdir()
    (earlier_dir / "genomes.json").write_text(
        json.dumps([{
            "contig_uids": [ids.contig_uid(s_shared).hex()],
            "genetic_code": 1,
            "kind": "prokaryote",
            "circular": [],
            "assembly_quality": "contig",
        }]),
        encoding="utf-8",
    )

    baited = [
        {
            "run_acc": "ERR1234567",
            "contig_id": "ERR1234567_10",
            "seq": s_new1,
            "circular": False,
            "trimmed_k": 0,
            "raw_len": len(s_new1),
            "bait_models": ["palm:mito_invertron_rpo"],
            "cluster_ids": ["ERR1234567_10_1"],
            "biome": "marine metagenome",
        },
        {
            "run_acc": "SRR7654321",
            "contig_id": "SRR7654321_22",
            "seq": s_new2,
            "circular": True,
            "trimmed_k": 31,
            "raw_len": len(s_new2) + 31,
            "bait_models": ["dpbb_single:crass_vrnap_dpbb"],
            "cluster_ids": ["SRR7654321_22_3"],
            "biome": "freshwater metatranscriptome [METATRANSCRIPTOMIC]",
        },
        {
            "run_acc": "DRR0000100",
            "contig_id": "DRR0000100_5",
            "seq": s_shared,
            "circular": False,
            "trimmed_k": 0,
            "raw_len": len(s_shared),
            "bait_models": ["palm:basal_protist_rpot"],
            "cluster_ids": ["DRR0000100_5_1"],
            "biome": "soil metagenome",
        },
    ]

    meta_dir = tmp_path / "meta_LGNRNAP1"
    fa_dir = tmp_path / "fa_LGNRNAP1"
    rep = m.stage_lgnrnap1(
        baited_contigs=baited,
        out_meta_dir=meta_dir,
        out_fasta_dir=fa_dir,
        owner_inputs={"ET1": earlier_dir},
    )
    assert rep["stage"] == "LGNRNAP1"
    assert rep["distinct_contigs"] == 3
    assert rep["shared_contigs"] == 1
    assert rep["net_new_contigs"] == 2
    assert rep["shared_by_stage"] == {"ET1": 1}
    assert rep["circular_contigs"] == 1

    proc = subprocess.run(
        [
            sys.executable,
            str(BSP / "s2_check_inputs.py"),
            "LGNRNAP1",
            str(meta_dir / "genomes.json"),
            str(meta_dir / "objects_staging.json"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, f"s2_check_inputs failed: {proc.stdout}\n{proc.stderr}"
    chk = json.loads(proc.stdout)
    assert chk["problems"] == []


def test_c68_conf_templates_pass_stage_est_validation():
    est = _stage_est()

    c_et1 = est.read_conf(STAGES / "ET1.conf", "ET1", tree=REPO)
    assert c_et1.stage == "ET1"
    assert c_et1.depth == "called"
    assert c_et1.min_contig_len == 1500
    assert c_et1.genecall_mode == "default"
    assert c_et1.crispr is False
    assert c_et1.after == ("LGRT1",)

    c_lgn = est.read_conf(STAGES / "LGNRNAP1.conf", "LGNRNAP1", tree=REPO)
    assert c_lgn.stage == "LGNRNAP1"
    assert c_lgn.depth == "called"
    assert c_lgn.min_contig_len == 2000
    assert c_lgn.genecall_mode == "gv-meta"
    assert c_lgn.crispr is True
    assert c_lgn.after == ("ET1",)

    # Every SUM digest in both confs must appear in its named NOTE file
    for conf in (c_et1, c_lgn):
        note_text = "\n".join((REPO / n).read_text(encoding="utf-8") for n in conf.notes)
        for h, name in conf.sums:
            assert h in note_text, f"{conf.stage}.conf SUM for {name} missing from note"
