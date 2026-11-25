"""Ring R1: unit tests for card W-38a (`scripts/stress/w38a_plm_bakeoff.py`; SD-21, MR-01, MR-03, S-11)."""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import pytest

pytestmark = [pytest.mark.covers("W-38a", "SD-21", "MR-01", "MR-03", "S-11")]

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "stress" / "w38a_plm_bakeoff.py"


def _load_mod():
    name = "w38a_plm_bakeoff"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, SCRIPT)
    assert spec is not None and spec.loader is not None, f"missing script: {SCRIPT}"
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_default_release_is_live_and_manifest_backed(tmp_path):
    """Default release avoids DEC-66-retiring r26.11.1; open_release_files refuses missing manifest without globbing."""
    m = _load_mod()
    assert m.DEFAULT_RELEASE != "r26.11.1"
    parser = m.build_parser()
    ns = parser.parse_args(["build-holdout", "--positives", "p.jsonl", "--out", "out"])
    assert ns.release == m.DEFAULT_RELEASE

    # Missing release under an empty root must fail closed via ReleaseFiles.of (no silent glob fallback)
    with pytest.raises(FileNotFoundError):
        m.open_release_files(tmp_path, "r26.11.2")


def test_strict_holdout_leakage_filter_boundaries_and_stratification():
    """Holdout requires strict <20% AA ID AND hmmsearch E > 10, fails closed on NaN/missing, and tags fusions."""
    m = _load_mod()
    candidates = [
        # Pass: remote homolog, monodomain
        {
            "feature_uid": "f01", "protein_uid": "p01", "system_id": "DS-1",
            "aa_len": 350, "n_domains": 1,
            "max_train_aa_id": 0.199, "min_train_hmm_evalue": 10.01, "min_train_foldseek_evalue": 0.05,
        },
        # Pass: novel-fold multi-domain fusion > 1,022 aa
        {
            "feature_uid": "f02", "protein_uid": "p02", "system_id": "DS-2",
            "aa_len": 1250, "n_domains": 3,
            "max_train_aa_id": 0.0, "min_train_hmm_evalue": math.inf, "min_train_foldseek_evalue": 4.2,
        },
        # Reject: exact 0.20 AA ID boundary
        {
            "feature_uid": "f03", "protein_uid": "p03", "system_id": "DS-3",
            "aa_len": 400, "n_domains": 1,
            "max_train_aa_id": 0.20, "min_train_hmm_evalue": 50.0, "min_train_foldseek_evalue": 5.0,
        },
        # Reject: exact E = 10.0 boundary
        {
            "feature_uid": "f04", "protein_uid": "p04", "system_id": "DS-4",
            "aa_len": 400, "n_domains": 1,
            "max_train_aa_id": 0.12, "min_train_hmm_evalue": 10.0, "min_train_foldseek_evalue": 5.0,
        },
        # Reject: NaN homology value (fail-closed)
        {
            "feature_uid": "f05", "protein_uid": "p05", "system_id": "DS-5",
            "aa_len": 700, "n_domains": 2,
            "max_train_aa_id": float("nan"), "min_train_hmm_evalue": 100.0,
        },
        # Reject: missing hmm E-value (fail-closed)
        {
            "feature_uid": "f06", "protein_uid": "p06", "system_id": "DS-6",
            "aa_len": 700, "n_domains": 2,
            "max_train_aa_id": 0.10, "min_train_hmm_evalue": None,
        },
    ]

    kept, rejected = m.filter_holdout_leakage(candidates)
    assert [r["feature_uid"] for r in kept] == ["f01", "f02"]
    assert [r["feature_uid"] for r in rejected] == ["f03", "f04", "f05", "f06"]

    by_f = {r["feature_uid"]: r for r in kept}
    assert by_f["f01"]["holdout_stratum"] == "remote_homolog"
    assert by_f["f01"]["is_fusion"] is False
    assert by_f["f01"]["exceeds_esm2_cap"] is False

    assert by_f["f02"]["holdout_stratum"] == "novel_fold"
    assert by_f["f02"]["is_fusion"] is True
    assert by_f["f02"]["exceeds_esm2_cap"] is True


def test_circular_origin_wrapping_spire_truncation_and_provenance():
    """Origin-spanning circular operon slices cleanly into [1500, 5000] bp linear fragment keeping feature_uid."""
    m = _load_mod()
    # 10,000 bp circular plasmid with an operon wrapping origin:
    # slot A: 9300..9900 (target), slot B: 9950..10400 (origin-wrapping end > contig_len, target),
    # slot C: 500..1200 (neighbor), slot D: 4500..5200 (distant, excluded from 3000 bp slice)
    locus = {
        "contig_uid": "aa" * 16,
        "contig_len": 10000,
        "is_circular": True,
        "anchor_feature_uid": "f_b",
        "target_feature_uids": ["f_a", "f_b"],
        "slots": [
            {"feature_uid": "f_c", "protein_uid": "p_c", "start": 500, "end": 1200, "strand": 1, "kind": 1},
            {"feature_uid": "f_d", "protein_uid": "p_d", "start": 4500, "end": 5200, "strand": -1, "kind": 1},
            {"feature_uid": "f_a", "protein_uid": "p_a", "start": 9300, "end": 9900, "strand": 1, "kind": 1},
            {"feature_uid": "f_b", "protein_uid": "p_b", "start": 9950, "end": 10400, "strand": 1, "kind": 1},
        ],
    }

    frag = m.truncate_locus_to_spire(locus, seed=42, target_len_bp=3000)
    assert frag["is_circular"] is False
    assert frag["parent_is_circular"] is True
    assert frag["contig_len"] == 3000
    assert frag["crosses_origin_in_parent"] is True
    assert frag["n_target_slots_retained"] == 2

    fuids = [s["parent_feature_uid"] for s in frag["slots"]]
    assert "f_a" in fuids and "f_b" in fuids and "f_c" in fuids
    assert "f_d" not in fuids

    # Check local coordinates are strictly inside [0, 3000], ordered, and edge distances are valid
    for idx, slot in enumerate(frag["slots"]):
        assert 0 <= slot["local_start"] < slot["local_end"] <= 3000
        assert slot["dist_to_contig_end_bp"] == min(slot["local_start"], 3000 - slot["local_end"])
        assert slot["orfs_to_contig_start"] == idx
        assert slot["orfs_to_contig_end"] == len(frag["slots"]) - 1 - idx

    # Linear contig near start clamps to [0, L] without negative coordinates and flags partial boundary ORF
    linear_locus = {
        "contig_uid": "bb" * 16,
        "contig_len": 20000,
        "is_circular": False,
        "anchor_feature_uid": "lin_1",
        "target_feature_uids": ["lin_1"],
        "slots": [
            {"feature_uid": "lin_1", "protein_uid": "lp1", "start": 100, "end": 900, "strand": 1, "kind": 1},
            {"feature_uid": "lin_2", "protein_uid": "lp2", "start": 1700, "end": 2600, "strand": 1, "kind": 1},
        ],
    }
    lfrag = m.truncate_locus_to_spire(linear_locus, seed=42, target_len_bp=2000)
    assert lfrag["window_start_bp"] == 0
    assert lfrag["window_end_bp"] == 2000
    by_f = {s["parent_feature_uid"]: s for s in lfrag["slots"]}
    assert by_f["lin_1"]["is_partial"] is False
    assert by_f["lin_2"]["is_partial"] is True
    assert by_f["lin_2"]["local_end"] == 2000


def test_determinism_and_input_order_invariance():
    """Holdout filtering, background hash-reservoir sampling, and locus truncation are input-order invariant."""
    m = _load_mod()
    bg_rows = [
        {
            "feature_uid": f"bg_{i:04d}",
            "protein_uid": f"pu_{i:04d}",
            "aa_len": 300 + (i % 50),
            "trap_cohort": ["ntpase", "methyltransferase", "transposase", "hydrolase", "general"][i % 5],
        }
        for i in range(200)
    ]
    s1 = m.sample_background_reservoir(bg_rows, total_size=40, seed=12345)
    s2 = m.sample_background_reservoir(list(reversed(bg_rows)), total_size=40, seed=12345)
    assert s1 == s2
    assert len(s1) == 40

    # Verify trap cohort quotas respected (15/15/15/10/45% of 40 = 6/6/6/4/18)
    counts = {}
    for r in s1:
        counts[r["trap_cohort"]] = counts.get(r["trap_cohort"], 0) + 1
    assert counts == {
        "general": 18,
        "hydrolase": 4,
        "methyltransferase": 6,
        "ntpase": 6,
        "transposase": 6,
    }


def _make_planted_holdout():
    """Construct a deterministic synthetic holdout with 20 positives across all strata + 1000 negatives."""
    rows = []
    # 10 remote_homolog positives (5 monodomain, 5 fusion including 2 > 1,022 aa)
    for i in range(10):
        is_fus = i >= 5
        rows.append({
            "feature_uid": f"pos_rem_{i}",
            "protein_uid": f"pp_rem_{i}",
            "system_id": f"SYS_REM_{i}",
            "label": 1,
            "holdout_stratum": "remote_homolog",
            "is_fusion": is_fus,
            "exceeds_esm2_cap": i >= 8,
            "trap_cohort": "positive",
        })
    # 10 novel_fold positives (5 monodomain, 5 fusion)
    for i in range(10):
        rows.append({
            "feature_uid": f"pos_nov_{i}",
            "protein_uid": f"pp_nov_{i}",
            "system_id": f"SYS_NOV_{i}",
            "label": 1,
            "holdout_stratum": "novel_fold",
            "is_fusion": i >= 5,
            "exceeds_esm2_cap": False,
            "trap_cohort": "positive",
        })
    # 1000 background negatives across all 5 cohorts
    cohorts = ["ntpase", "methyltransferase", "transposase", "hydrolase", "general"]
    for i in range(1000):
        rows.append({
            "feature_uid": f"neg_{i:04d}",
            "protein_uid": f"np_{i:04d}",
            "system_id": "",
            "label": 0,
            "holdout_stratum": "background",
            "is_fusion": False,
            "exceeds_esm2_cap": False,
            "trap_cohort": cohorts[i % len(cohorts)],
        })
    return rows


def test_metrics_oracle_and_gates_pass_and_fail_scenarios(tmp_path):
    """Exact PR-AUC / Recall@FPR against hand-calculated values, plus planted Gate 1-3 pass & fail checks."""
    m = _load_mod()

    # Exact hand-calculated oracle: y=[1, 0, 1, 0], scores=[0.9, 0.8, 0.7, 0.1]
    # Sorted: (1 @ rank 1 -> P=1.0), (0 @ rank 2), (1 @ rank 3 -> P=2/3), (0 @ rank 4)
    # AP = (1.0 + 2/3) / 2 = 5/6
    assert m.average_precision([1, 0, 1, 0], [0.9, 0.8, 0.7, 0.1]) == pytest.approx(5.0 / 6.0)
    assert m.recall_at_k([1, 0, 1, 0], [0.9, 0.8, 0.7, 0.1], 2) == pytest.approx(0.5)

    holdout = _make_planted_holdout()
    base_neg = {f"neg_{i:04d}": 0.1 + (i / 10000.0) for i in range(1000)}  # max neg = 0.1999

    def _scores(high_fuids: set[str]) -> dict[str, float]:
        s = dict(base_neg)
        for r in holdout:
            if r["label"] == 1:
                s[r["feature_uid"]] = 0.95 if r["feature_uid"] in high_fuids else 0.05
        return s

    # --- Scenario A: All 3 gates PASS; ProstT5 within 2 pp of ESM-2 3B 3Di -> ProstT5 selected ---
    # remote_and_fusion has 15 positives (10 remote + 5 novel fusions)
    jack_hits = {f"pos_rem_{i}" for i in range(3)}  # 3/15 = 20%
    prostt5_hits = {f"pos_rem_{i}" for i in range(10)} | {"pos_nov_5", "pos_nov_6"}  # 12/15 = 80% (+60 pp)
    esm2_3b_hits = prostt5_hits  # same recall -> ProstT5 preferred (within 2 pp)
    hhblits_hits = {f"pos_rem_{i}" for i in range(5)}
    # Option B recovers 4 marginal novel systems (4/20 = 20% >= 10%) with zero fusion loss
    esmc_hits = prostt5_hits | {"pos_nov_0", "pos_nov_1", "pos_nov_2", "pos_nov_3", "pos_nov_7", "pos_nov_8", "pos_nov_9"}
    # Gate 3: nbr_gba gets 2/10 novel_fold; geneclr gets 8/10 (4x >= 2x), 10% truncation drop (<15%), 0% NR90 drop
    nbr_hits = {"pos_nov_0", "pos_nov_1"}
    geneclr_hits = {f"pos_nov_{i}" for i in range(8)} | {f"pos_rem_{i}" for i in range(6)}
    geneclr_trunc = {f"pos_nov_{i}" for i in range(7)} | {f"pos_rem_{i}" for i in range(6)}  # 13 vs 14 -> 7.1% drop

    eval_pass = {
        "methods": {
            "jackhmmer_3iter": m.evaluate_method(holdout, _scores(jack_hits)),
            "hhblits_3iter": m.evaluate_method(holdout, _scores(hhblits_hits)),
            "prostt5_foldseek": m.evaluate_method(holdout, _scores(prostt5_hits), latency_ratio_vs_afdb=1.1),
            "esm2_3b_3di_foldseek": m.evaluate_method(holdout, _scores(esm2_3b_hits), latency_ratio_vs_afdb=1.3),
            "esmc_300m_segments_ivfpq": m.evaluate_method(holdout, _scores(esmc_hits)),
            "nbr_gba": m.evaluate_method(holdout, _scores(nbr_hits)),
            "esm_df_geneclr_df": m.evaluate_method(
                holdout, _scores(geneclr_hits),
                truncated_scores=_scores(geneclr_trunc),
                nr90_scores=_scores(geneclr_hits),
            ),
        }
    }
    verdict_pass = m.check_gates(eval_pass)
    assert verdict_pass["gate1_3di"]["pass"] is True
    assert verdict_pass["gate1_3di"]["recommended_encoder"] == "prostt5_foldseek"
    assert verdict_pass["gate2_dense_ann"]["pass"] is True
    assert verdict_pass["gate3_tape_context"]["pass"] is True
    assert verdict_pass["gate3_tape_context"]["action"] == "BUILD_OPTION_C_AND_D"

    # --- Scenario B: All 3 gates FAIL ---
    # Gate 1 fails (+6.7 pp < +20 pp); Gate 2 fails (0 marginal systems -> REJECT_OPTION_B_PERMANENTLY);
    # Gate 3 fails (50% truncation drop >= 15% -> KEEP_OPTION_D_ONLY)
    weak_3di = jack_hits | {"pos_rem_3"}
    geneclr_bad_trunc = {f"pos_nov_{i}" for i in range(3)} | {f"pos_rem_{i}" for i in range(4)}  # 50% drop
    eval_fail = {
        "methods": {
            "jackhmmer_3iter": m.evaluate_method(holdout, _scores(jack_hits)),
            "hhblits_3iter": m.evaluate_method(holdout, _scores(hhblits_hits)),
            "prostt5_foldseek": m.evaluate_method(holdout, _scores(weak_3di), latency_ratio_vs_afdb=1.1),
            "esmc_300m_full_fp16": m.evaluate_method(holdout, _scores(weak_3di)),
            "nbr_gba": m.evaluate_method(holdout, _scores(nbr_hits)),
            "defense_predictor_fixed": m.evaluate_method(
                holdout, _scores(geneclr_hits),
                truncated_scores=_scores(geneclr_bad_trunc),
                nr90_scores=_scores(geneclr_hits),
            ),
        }
    }
    verdict_fail = m.check_gates(eval_fail)
    assert verdict_fail["gate1_3di"]["pass"] is False
    assert verdict_fail["gate1_3di"]["action"] == "HOLD_3DI"
    assert verdict_fail["gate2_dense_ann"]["pass"] is False
    assert verdict_fail["gate2_dense_ann"]["action"] == "REJECT_OPTION_B_PERMANENTLY"
    assert verdict_fail["gate3_tape_context"]["pass"] is False
    assert verdict_fail["gate3_tape_context"]["action"] == "KEEP_OPTION_D_ONLY"

    # --- End-to-end CLI round-trip (`evaluate` -> `check-gates`) in tmp_path ---
    holdout_path = tmp_path / "holdout.jsonl"
    holdout_path.write_text("\n".join(json.dumps(r) for r in holdout) + "\n", encoding="utf-8")
    (tmp_path / "jack.json").write_text(json.dumps(_scores(jack_hits)), encoding="utf-8")
    (tmp_path / "prostt5.json").write_text(json.dumps(_scores(prostt5_hits)), encoding="utf-8")
    pred_manifest = {
        "jackhmmer_3iter": {"scores": "jack.json"},
        "prostt5_foldseek": {"scores": "prostt5.json", "latency_ratio_vs_afdb": 1.05},
    }
    manifest_path = tmp_path / "preds.json"
    manifest_path.write_text(json.dumps(pred_manifest), encoding="utf-8")
    eval_out = tmp_path / "eval.json"
    gates_out = tmp_path / "gates.json"

    assert m.main(["evaluate", "--holdout", str(holdout_path), "--predictions", str(manifest_path), "--out", str(eval_out)]) == 0
    assert m.main(["check-gates", "--eval-json", str(eval_out), "--out", str(gates_out)]) == 0
    cli_gates = json.loads(gates_out.read_text(encoding="utf-8"))
    assert cli_gates["gate1_3di"]["pass"] is True
    assert cli_gates["gate1_3di"]["recommended_encoder"] == "prostt5_foldseek"
