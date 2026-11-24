"""Ring R1: unit tests for `bio-box search panel` (card W-39; SD-04, S-03, S-04, MR-12)."""

from __future__ import annotations

import io
from pathlib import Path

import pyhmmer
import pytest

from biobox import ids
from biobox.search import engine, panel

pytestmark = [pytest.mark.covers("W-39", "SD-04", "S-03", "S-04", "MR-12")]


def _make_hmm_bytes(name: str, seq: str, *, acc: str | None = None, ga: tuple[float, float] | None = None,
                    tc: tuple[float, float] | None = None) -> bytes:
    hmm = engine.build_from_sequence(name, seq)
    if acc is not None:
        hmm.accession = acc.encode()
    if ga is not None:
        hmm.cutoffs.gathering = ga
    if tc is not None:
        hmm.cutoffs.trusted = tc
    buf = io.BytesIO()
    hmm.write(buf)
    return buf.getvalue()


_SEQ_A = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"
_SEQ_B = "MSALVNVLGGVAGKTTLSRLIAQALRGQGFSVLVIDGDPVSRNLAEALGLPVVEVATGDLAEALAALRDAGVDVLIVDTPGIGDLAELARLLGVRVVVVDAADLAELARALGLPVVEVATGDLAEALAALRDAGVDVLIVDTPGIGDLAELARLLGVRVVVVDAADLAELAR"
_SEQ_C = "MGEKLYVLDVDSGAIQALRRAGIPVVEAGSGQRLADALRAAGVDVVVTDPGIGDLAELARLLGVRVVVVDAADLAELARALGLPVVEVATGDLAEALAALRDAGVDVLIVDTPGIGDLAELARLLGVRVVVV"


def test_model_key_assignment_and_duplicate_refusal(tmp_path):
    """model_key is [NS:]ACC when ACC is present and unique, else [NS:]NAME; duplicate keys refuse."""
    f1 = tmp_path / "m1.hmm"
    f1.write_bytes(
        _make_hmm_bytes("ModelA", _SEQ_A, acc="PF00001.1", ga=(25.0, 25.0))
        + _make_hmm_bytes("ModelB", _SEQ_B, acc=None, ga=(20.0, 20.0))
    )
    loaded = panel.load_panel_specs(
        hmm_specs=[f"rt:{f1}"],
        library_specs=[],
        decoy_hmm_specs=[],
        decoy_library_specs=[],
        cutoff="gathering",
    )
    keys = [m.model_key for m in loaded]
    assert keys == ["rt:PF00001.1", "rt:ModelB"]
    assert all(not m.is_decoy for m in loaded)
    assert all(len(m.hmm_sha256) == 64 for m in loaded)

    # Same content in a different path produces identical hmm_sha256 (DATE stripped)
    f2 = tmp_path / "copy.hmm"
    f2.write_bytes(f1.read_bytes())
    loaded2 = panel.load_panel_specs(
        hmm_specs=[f"rt:{f2}"],
        library_specs=[],
        decoy_hmm_specs=[],
        decoy_library_specs=[],
        cutoff="gathering",
    )
    assert [m.hmm_sha256 for m in loaded2] == [m.hmm_sha256 for m in loaded]

    # Duplicate keys across files raise PanelError
    with pytest.raises(panel.PanelError, match="duplicate"):
        panel.load_panel_specs(
            hmm_specs=[f"rt:{f1}", f"rt:{f2}"],
            library_specs=[],
            decoy_hmm_specs=[],
            decoy_library_specs=[],
            cutoff="gathering",
        )


def test_cutoffs_tsv_and_heterogeneous_cutoffs(tmp_path):
    """Mix GA, fallback-bits, and cutoffs-tsv (with i_evalue_max, hmm_cov_min, min_dom_frac); missing cutoff refuses."""
    f1 = tmp_path / "het.hmm"
    f1.write_bytes(
        _make_hmm_bytes("HasGA", _SEQ_A, acc="PF10001", ga=(27.5, 22.0))
        + _make_hmm_bytes("NoGA", _SEQ_B, acc=None, ga=None)
        + _make_hmm_bytes("Custom", _SEQ_C, acc="CUST1", ga=None)
    )
    # Without fallback or TSV, NoGA and Custom must be refused
    with pytest.raises(panel.PanelError, match="cutoff"):
        panel.load_panel_specs(
            hmm_specs=[str(f1)],
            library_specs=[],
            decoy_hmm_specs=[],
            decoy_library_specs=[],
            cutoff="gathering",
        )

    tsv = tmp_path / "cutoffs.tsv"
    tsv.write_text(
        "# model_key\tseq_bits\tdom_bits\ti_evalue_max\thmm_cov_min\tmin_dom_frac\n"
        "CUST1\t32.0\t20.0\t1e-5\t0.5\t0.45\n",
        encoding="utf-8",
    )
    loaded = panel.load_panel_specs(
        hmm_specs=[str(f1)],
        library_specs=[],
        decoy_hmm_specs=[],
        decoy_library_specs=[],
        cutoff="gathering",
        fallback_bits=30.0,
        cutoffs_tsv=tsv,
    )
    by_key = {m.model_key: m for m in loaded}
    assert by_key["PF10001"].seq_bits == pytest.approx(27.5)
    assert by_key["PF10001"].dom_bits == pytest.approx(22.0)
    assert by_key["NoGA"].seq_bits == pytest.approx(30.0)
    assert by_key["CUST1"].seq_bits == pytest.approx(32.0)
    assert by_key["CUST1"].dom_bits == pytest.approx(20.0)
    assert by_key["CUST1"].i_evalue_max == pytest.approx(1e-5)
    assert by_key["CUST1"].hmm_cov_min == pytest.approx(0.5)
    assert by_key["CUST1"].min_dom_frac == pytest.approx(0.45)


def test_min_dom_frac_rejects_repeat_stacking():
    """Repeat-stacked multi-domain hits (low max_dom_bits / seq_bits) are rejected by min_dom_frac."""
    # Repeat-stacked protein: 8 weak 6-bit domains summing to 48 seq bits -> dom_seq_frac = 0.125
    repeat_hit = {
        "protein_uid": b"\x01" * 16,
        "model_key": "art:DISTORF",
        "seq_bitscore": 48.0,
        "max_dom_bits": 6.0,
        "n_model_domains": 8,
    }
    # Genuine single-domain hit: 45 seq bits, 45 dom bits -> dom_seq_frac = 1.0
    real_hit = {
        "protein_uid": b"\x02" * 16,
        "model_key": "art:DISTORF",
        "seq_bitscore": 45.0,
        "max_dom_bits": 45.0,
        "n_model_domains": 1,
    }
    assert panel.passes_min_dom_frac(repeat_hit["max_dom_bits"], repeat_hit["seq_bitscore"], 0.4) is False
    assert panel.passes_min_dom_frac(real_hit["max_dom_bits"], real_hit["seq_bitscore"], 0.4) is True
    assert panel.passes_min_dom_frac(repeat_hit["max_dom_bits"], repeat_hit["seq_bitscore"], None) is True


def test_overlap_competition_decoy_veto_and_fusion_retention():
    """Overlap competition sets best_on_overlap, best_in_namespace, and vetoed_by_decoy while keeping fusions."""
    p1 = b"\x11" * 16
    p2 = b"\x22" * 16
    p3 = b"\x33" * 16
    rows = [
        # p1: competing targets in two namespaces
        {"protein_uid": p1, "model_key": "pfam:RVT_1", "namespace": "pfam", "is_decoy": False,
         "bitscore": 80.0, "i_evalue_zref": 1e-20, "env_start": 10, "env_end": 110},
        {"protein_uid": p1, "model_key": "myrt:RT_A", "namespace": "myrt", "is_decoy": False,
         "bitscore": 60.0, "i_evalue_zref": 1e-15, "env_start": 15, "env_end": 105},
        # p2: decoy outscores target on overlapping region -> vetoed_by_decoy=True
        {"protein_uid": p2, "model_key": "decoy:RdRP_1", "namespace": "decoy", "is_decoy": True,
         "bitscore": 95.0, "i_evalue_zref": 1e-25, "env_start": 20, "env_end": 120},
        {"protein_uid": p2, "model_key": "pfam:RVT_1", "namespace": "pfam", "is_decoy": False,
         "bitscore": 42.0, "i_evalue_zref": 1e-9, "env_start": 25, "env_end": 115},
        # p3: 3-domain fusion (non-overlapping envelopes) -> all 3 retain best_on_overlap=True
        {"protein_uid": p3, "model_key": "art:NTERM", "namespace": "art", "is_decoy": False,
         "bitscore": 70.0, "i_evalue_zref": 1e-18, "env_start": 10, "env_end": 100},
        {"protein_uid": p3, "model_key": "pfam:RVT_1", "namespace": "pfam", "is_decoy": False,
         "bitscore": 85.0, "i_evalue_zref": 1e-22, "env_start": 120, "env_end": 320},
        {"protein_uid": p3, "model_key": "art:CTERM", "namespace": "art", "is_decoy": False,
         "bitscore": 55.0, "i_evalue_zref": 1e-12, "env_start": 340, "env_end": 410},
    ]
    panel.annotate_overlap_competition(rows)
    by = {(r["protein_uid"], r["model_key"]): r for r in rows}

    assert by[(p1, "pfam:RVT_1")]["best_on_overlap"] is True
    assert by[(p1, "myrt:RT_A")]["best_on_overlap"] is False
    assert by[(p1, "myrt:RT_A")]["best_in_namespace"] is True
    assert by[(p1, "myrt:RT_A")]["vetoed_by_decoy"] is False

    assert by[(p2, "decoy:RdRP_1")]["best_on_overlap"] is True
    assert by[(p2, "decoy:RdRP_1")]["vetoed_by_decoy"] is False
    assert by[(p2, "pfam:RVT_1")]["best_on_overlap"] is False
    assert by[(p2, "pfam:RVT_1")]["vetoed_by_decoy"] is True

    for mk in ("art:NTERM", "pfam:RVT_1", "art:CTERM"):
        assert by[(p3, mk)]["best_on_overlap"] is True
        assert by[(p3, mk)]["vetoed_by_decoy"] is False


def test_circular_origin_spanning_operon_and_linear_censoring():
    """Origin-spanning hits on a circular contig merge into one cluster; linear contigs censor at ends."""
    cu_circ = ids.contig_uid(b"ACGT" * 5000)
    cu_lin = ids.contig_uid(b"TGCA" * 5000)
    # Circular contig: 10 slots, length 10000; hits at slot 9 (9100..9900) and slot 0 (100..900) -> 0 intervening genes, 200 bp gap
    f9 = ids.feature_uid(cu_circ, 9100, 9900, 1, ids.KIND["cds"])
    f0 = ids.feature_uid(cu_circ, 100, 900, 1, ids.KIND["cds"])
    tape_circ = {
        "contig_uid": cu_circ, "contig_len": 10000, "is_circular": True, "n_slots": 10,
        "feature_uids": [f0] + [ids.feature_uid(cu_circ, i * 1000, i * 1000 + 800, 1, ids.KIND["cds"]) for i in range(1, 9)] + [f9],
        "start_bps": [100] + [i * 1000 for i in range(1, 9)] + [9100],
        "end_bps": [900] + [i * 1000 + 800 for i in range(1, 9)] + [9900],
        "strands": [1] * 10,
        "protein_uids": [bytes([i]) * 16 for i in range(10)],
    }
    orf_hits_circ = [
        {"contig_uid": cu_circ, "feature_uid": f0, "ordinal": 0, "start": 100, "end": 900, "strand": 1,
         "protein_uid": b"\x00" * 16, "model_key": "Cas7f", "model_name": "Cas7f", "model_len": 200,
         "bitscore": 80.0, "hmm_from": 0, "hmm_to": 200, "is_decoy": False, "best_on_overlap": True,
         "contig_accession": "CIRC1", "stage": "S0", "genome_kind": "plasmid"},
        {"contig_uid": cu_circ, "feature_uid": f9, "ordinal": 9, "start": 9100, "end": 9900, "strand": 1,
         "protein_uid": b"\x09" * 16, "model_key": "TniQ", "model_name": "TniQ", "model_len": 200,
         "bitscore": 90.0, "hmm_from": 0, "hmm_to": 200, "is_decoy": False, "best_on_overlap": True,
         "contig_accession": "CIRC1", "stage": "S0", "genome_kind": "plasmid"},
    ]
    clusters = panel.call_loci_for_contig(
        tape_circ, orf_hits_circ, max_gap_genes=5, max_gap_bp=10000, min_models=2, neighbors=5,
    )
    assert len(clusters) == 1
    c = clusters[0]
    assert c["censored"] is False
    assert c["n_models"] == 2
    assert c["n_panel_orfs"] == 2
    assert c["models_ordered"] == "TniQ(+)>Cas7f(+)"
    assert c["span_bp"] == 1800

    # Linear contig with same slot positions: does NOT wrap, so 2 separate 1-model clusters -> 0 clusters at min_models=2
    tape_lin = {**tape_circ, "contig_uid": cu_lin, "is_circular": False}
    orf_hits_lin = [{**h, "contig_uid": cu_lin} for h in orf_hits_circ]
    clusters_lin = panel.call_loci_for_contig(
        tape_lin, orf_hits_lin, max_gap_genes=5, max_gap_bp=10000, min_models=2, neighbors=5,
    )
    assert len(clusters_lin) == 0
    # At min_models=1, both edge clusters are flagged censored=True
    clusters_lin1 = panel.call_loci_for_contig(
        tape_lin, orf_hits_lin, max_gap_genes=5, max_gap_bp=10000, min_models=1, neighbors=5,
    )
    assert len(clusters_lin1) == 2
    assert all(x["censored"] is True for x in clusters_lin1)


def test_cluster_filters_min_models_distinct_orfs_require_forbid():
    """--distinct-orfs, --require, and --forbid filter clusters accurately."""
    cu = ids.contig_uid(b"ACTG" * 2000)
    f0 = ids.feature_uid(cu, 1000, 2000, 1, ids.KIND["cds"])
    f1 = ids.feature_uid(cu, 2100, 3000, 1, ids.KIND["cds"])
    tape = {
        "contig_uid": cu, "contig_len": 8000, "is_circular": False, "n_slots": 2,
        "feature_uids": [f0, f1], "start_bps": [1000, 2100], "end_bps": [2000, 3000],
        "strands": [1, 1], "protein_uids": [b"\x01" * 16, b"\x02" * 16],
    }
    # Fusion ORF f0 hits both M1 and M2
    fusion_hits = [
        {"contig_uid": cu, "feature_uid": f0, "ordinal": 0, "start": 1000, "end": 2000, "strand": 1,
         "protein_uid": b"\x01" * 16, "model_key": "M1", "model_name": "M1", "model_len": 150,
         "bitscore": 70.0, "hmm_from": 0, "hmm_to": 150, "is_decoy": False, "best_on_overlap": True,
         "contig_accession": "LIN1", "stage": "S0", "genome_kind": "virus"},
        {"contig_uid": cu, "feature_uid": f0, "ordinal": 0, "start": 1000, "end": 2000, "strand": 1,
         "protein_uid": b"\x01" * 16, "model_key": "M2", "model_name": "M2", "model_len": 150,
         "bitscore": 65.0, "hmm_from": 0, "hmm_to": 150, "is_decoy": False, "best_on_overlap": True,
         "contig_accession": "LIN1", "stage": "S0", "genome_kind": "virus"},
    ]
    assert len(panel.call_loci_for_contig(tape, fusion_hits, min_models=2, distinct_orfs=False)) == 1
    assert len(panel.call_loci_for_contig(tape, fusion_hits, min_models=2, distinct_orfs=True)) == 0
    assert len(panel.call_loci_for_contig(tape, fusion_hits, min_models=2, require={"M1"})) == 1
    assert len(panel.call_loci_for_contig(tape, fusion_hits, min_models=2, require={"M3"})) == 0
    assert len(panel.call_loci_for_contig(tape, fusion_hits, min_models=2, forbid={"M2"})) == 0


def test_four_locus_diagnostics():
    """Verify max_adjacent_igr_bp, bridged_noncoding_bp, cas_less_crispr, and split_model_pair."""
    cu = ids.contig_uid(b"GCTA" * 4000)
    # 5 slots:
    # slot 0: hit N-half of RVT (1000..1500, +1)
    # slot 1: hit C-half of RVT (1520..2000, +1) -> split_model_pair=True
    # IGR (2000..2600, 600 bp)
    # slot 2: unannotated smORF (2600..2780, 60 aa < 100 aa, no domains)
    # IGR (2780..3580, 800 bp) -> bridged_noncoding_bp = 3580 - 2000 = 1580 bp!
    # slot 3: normal protein (3580..4500, 306 aa)
    f = [
        ids.feature_uid(cu, 1000, 1500, 1, ids.KIND["cds"]),
        ids.feature_uid(cu, 1520, 2000, 1, ids.KIND["cds"]),
        ids.feature_uid(cu, 2600, 2780, 1, ids.KIND["cds"]),
        ids.feature_uid(cu, 3580, 4500, 1, ids.KIND["cds"]),
    ]
    p = [bytes([i + 1]) * 16 for i in range(4)]
    tape = {
        "contig_uid": cu, "contig_len": 16000, "is_circular": False, "n_slots": 4,
        "feature_uids": f, "start_bps": [1000, 1520, 2600, 3580], "end_bps": [1500, 2000, 2780, 4500],
        "strands": [1, 1, 1, 1], "protein_uids": p,
    }
    hits = [
        {"contig_uid": cu, "feature_uid": f[0], "ordinal": 0, "start": 1000, "end": 1500, "strand": 1,
         "protein_uid": p[0], "model_key": "pfam:PF00078", "model_name": "RVT_1", "model_len": 240,
         "bitscore": 55.0, "hmm_from": 5, "hmm_to": 115, "is_decoy": False, "best_on_overlap": True,
         "contig_accession": "DIAG1", "stage": "S0", "genome_kind": "virus"},
        {"contig_uid": cu, "feature_uid": f[1], "ordinal": 1, "start": 1520, "end": 2000, "strand": 1,
         "protein_uid": p[1], "model_key": "pfam:PF00078", "model_name": "RVT_1", "model_len": 240,
         "bitscore": 60.0, "hmm_from": 125, "hmm_to": 235, "is_decoy": False, "best_on_overlap": True,
         "contig_accession": "DIAG1", "stage": "S0", "genome_kind": "virus"},
    ]
    igrs = [
        {"contig_uid": cu, "start": 2000, "end": 2600, "strand": 1, "left_feature_uid": f[1], "right_feature_uid": f[2]},
        {"contig_uid": cu, "start": 2780, "end": 3580, "strand": 1, "left_feature_uid": f[2], "right_feature_uid": f[3]},
    ]
    prot_lens = {p[0]: 166, p[1]: 160, p[2]: 60, p[3]: 306}
    prot_doms_no_cas = {p[3]: {"SomeOtherDom"}}
    crispr_arrays = [{"start": 2100, "end": 2500}]

    res = panel.call_loci_for_contig(
        tape, hits, min_models=1, igrs=igrs, protein_lengths=prot_lens,
        protein_domains=prot_doms_no_cas, crispr_arrays=crispr_arrays,
    )
    assert len(res) == 1
    c = res[0]
    assert c["split_model_pair"] is True
    assert c["max_adjacent_igr_bp"] == 800
    assert c["bridged_noncoding_bp"] == 1580
    assert c["n_crispr_arrays"] == 1
    assert c["cas_less_crispr"] is True

    # Adding a Cas domain within +-10 genes flips cas_less_crispr to False
    prot_doms_with_cas = {p[3]: {"Cas_Csy3"}}
    res2 = panel.call_loci_for_contig(
        tape, hits, min_models=1, igrs=igrs, protein_lengths=prot_lens,
        protein_domains=prot_doms_with_cas, crispr_arrays=crispr_arrays,
    )
    assert res2[0]["cas_less_crispr"] is False


def test_substep_planner_cpu_min_bound():
    """plan_tasks splits (shard_slice, model_batch) so every sub-step stays <= 200 CPU-min."""
    shards = [{"path": f"unique/part-{i}.faa.zst", "n_residues": 5_000_000_000} for i in range(4)]
    models = [{"model_key": f"M{i}", "model_len": 250} for i in range(40)]
    tasks = panel.plan_tasks(shards, models, max_cpu_min=200.0)
    assert len(tasks) > 1
    assert all(t["est_cpu_min"] <= 200.0 + 1e-6 for t in tasks)
    # Every (shard, model) pair covered exactly once
    covered = {(s, m) for t in tasks for s in t["shard_indices"] for m in t["model_indices"]}
    expected = {(s, m) for s in range(len(shards)) for m in range(len(models))}
    assert covered == expected


def test_sql_panel_name_validation():
    """Panel names in `bio-box sql --panel NAME=DIR` must be valid identifiers and not collide with TABLES."""
    assert panel.validate_panel_name("rt") == "rt"
    assert panel.validate_panel_name("prok_rt_1") == "prok_rt_1"
    with pytest.raises(ValueError):
        panel.validate_panel_name("1bad")
    with pytest.raises(ValueError):
        panel.validate_panel_name("bad-name")
    with pytest.raises(ValueError):
        panel.validate_panel_name("l1_domains")
