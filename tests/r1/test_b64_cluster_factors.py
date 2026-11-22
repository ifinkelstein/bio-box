"""Ring R1: card B-64 — refit release clustering calibration factors on two releases (`r26.11.1` and `r26.11.2`).

Oracles (*measured* on wilkcomp02, 2026-11-21, read-only from `/tmp/if2453/r10/root`):
- Published `l2/cluster_plan/release=<R>/cluster.plan.json` and `*.cut.json` and `runs/stage=<R>/*.json`
  for `rM.2` (`45,815,843` proteins, legacy pre-B-54 cuts, 24 parts), `r26.11.1` (`49,102,041` proteins,
  legacy pre-B-54 cuts, 24 parts), and `r26.11.2` (`57,780,176` proteins, post-B-54 cuts at `part_cpu_min = 240.0`).
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog, provenance
from biobox.l1.genecall import GENE_CALL_VERSION as GCV
from biobox.l2 import cluster as l2

pytestmark = pytest.mark.covers("B-64", "B-54", "C-46")

REPO = Path(__file__).resolve().parents[2]
POD_DIR = REPO / "scripts" / "stress" / "build_stage_pod"
PROFILE = REPO / "profiles" / "wilkcomp02.yaml"


def _te():
    if str(POD_DIR) not in sys.path:
        sys.path.insert(0, str(POD_DIR))
    return __import__("train_est")


# Measured split-call oracles from /tmp/if2453/r10/root (wilkcomp02, 2026-11-21):
#   call: (n_parts, total_keys, total_units, overhead_cpu_min, stored_work_cpu_min, meas_sum, lo_part_meas, hi_part_meas)
RM2_CUTS = {
    "cluster.nr50.01.rescorediagonal": (24, 31_183_876, 4_840_240_391, 0.0, 14.230, 43.092, 1.677, 1.959),
    "cluster.nr50.06.rescorediagonal": (24, 17_348_804, 1_767_125_978, 0.0, 7.422, 14.749, 0.592, 0.650),
    "cluster.nr50.07.align": (24, 17_348_804, 698_116_275, 0.0, 309.014, 279.402, 10.479, 12.412),
    "cluster.nr50.11.prefilter": (24, 13_520_716, 4_070_572_296, 47.513, 114.852, 1307.653, 53.262, 56.853),
    "cluster.nr50.12.align": (24, 13_520_716, 3_038_129_206, 0.0, 992.040, 1259.101, 48.659, 54.859),
    "cluster.nr50.15.prefilter": (24, 8_915_865, 2_601_223_781, 0.701, 298.080, 421.930, 17.077, 18.256),
    "cluster.nr50.16.align": (24, 8_915_865, 2_717_279_156, 0.0, 785.742, 707.665, 27.332, 31.417),
    "cluster.nr50.19.prefilter": (24, 8_312_546, 2_454_500_276, 0.654, 708.012, 957.957, 37.821, 41.582),
    "cluster.nr50.20.align": (24, 8_312_546, 2_611_996_035, 0.0, 880.178, 685.875, 26.692, 29.623),
    "cluster.nr50.members": (24, 8_192_065, 409_423_542, 0.0, 360.995, 332.320, 11.977, 14.904),
    "cluster.nr90.01.rescorediagonal": (24, 45_815_843, 3_177_657_566, 0.0, 17.382, 49.190, 1.889, 2.203),
    "cluster.nr90.06.rescorediagonal": (24, 33_085_999, 1_705_376_160, 0.0, 13.677, 22.746, 0.886, 1.028),
    "cluster.nr90.07.align": (24, 33_085_999, 1_356_134_142, 0.0, 365.383, 340.227, 13.259, 15.418),
    "cluster.nr90.11.prefilter": (24, 31_920_514, 10_194_781_853, 112.173, 690.781, 3598.107, 145.797, 156.429),
    "cluster.nr90.12.align": (24, 31_920_514, 8_616_200_048, 0.0, 3161.456, 3500.653, 141.469, 148.995),
    "cluster.nr90.members": (24, 31_183_876, 432_415_353, 0.0, 132.674, 148.455, 4.861, 9.346),
}

R26_11_1_CUTS = {
    "cluster.nr50.01.rescorediagonal": (24, 33_117_773, 5_122_699_229, 0.0, 15.061, 44.221, 1.708, 2.050),
    "cluster.nr50.06.rescorediagonal": (24, 18_876_326, 1_929_507_605, 0.0, 8.104, 16.208, 0.631, 0.716),
    "cluster.nr50.07.align": (24, 18_876_326, 758_111_957, 0.0, 338.841, 289.749, 10.925, 12.984),
    "cluster.nr50.11.prefilter": (24, 14_912_484, 4_385_758_277, 52.404, 127.398, 1433.588, 57.972, 61.754),
    "cluster.nr50.12.align": (24, 14_912_484, 3_333_344_123, 0.0, 1088.437, 1337.461, 51.441, 57.807),
    "cluster.nr50.15.prefilter": (24, 10_064_063, 2_852_219_962, 0.791, 335.328, 466.314, 18.727, 20.858),
    "cluster.nr50.16.align": (24, 10_064_063, 3_073_734_693, 0.0, 896.236, 755.212, 29.591, 33.357),
    "cluster.nr50.19.prefilter": (24, 9_421_654, 2_698_588_280, 0.741, 799.736, 1373.837, 55.185, 60.154),
    "cluster.nr50.20.align": (24, 9_421_654, 2_964_472_779, 0.0, 1013.591, 743.502, 29.228, 32.596),
    "cluster.nr50.members": (24, 9_291_312, 440_098_571, 0.0, 390.468, 356.505, 13.333, 15.872),
    "cluster.nr90.01.rescorediagonal": (24, 49_102_041, 3_347_121_799, 0.0, 18.309, 50.854, 2.018, 2.264),
    "cluster.nr90.06.rescorediagonal": (24, 35_153_162, 1_767_760_409, 0.0, 14.177, 22.817, 0.889, 1.053),
    "cluster.nr90.07.align": (24, 35_153_162, 1_412_500_763, 0.0, 385.612, 340.406, 13.429, 15.109),
    "cluster.nr90.11.prefilter": (24, 33_882_143, 10_664_275_036, 119.066, 755.872, 3832.621, 156.004, 165.790),
    "cluster.nr90.12.align": (24, 33_882_143, 9_125_396_776, 0.0, 3348.291, 3656.605, 148.257, 156.607),
    "cluster.nr90.members": (24, 33_117_773, 463_925_032, 0.0, 142.341, 162.484, 5.221, 10.204),
}

R26_11_2_CUTS = {
    "cluster.nr50.01.rescorediagonal": (1, 35_475_280, 5_488_078_645, 0.0, 16.135, 17.273, 17.273, 17.273),
    "cluster.nr50.06.rescorediagonal": (1, 20_249_969, 2_078_959_494, 0.0, 8.732, 8.695, 8.695, 8.695),
    "cluster.nr50.07.align": (2, 20_249_969, 821_248_164, 0.0, 320.051, 311.919, 151.696, 160.223),
    "cluster.nr50.11.prefilter": (1, 16_085_986, 4_633_835_937, 56.528, 181.444, 232.881, 232.881, 232.881),
    "cluster.nr50.12.align": (6, 16_085_986, 3_597_818_614, 0.0, 1429.454, 1409.127, 222.245, 242.589),
    "cluster.nr50.15.prefilter": (3, 10_927_842, 3_029_785_935, 0.859, 499.261, 495.159, 159.513, 169.014),
    "cluster.nr50.16.align": (4, 10_927_842, 3_343_819_070, 0.0, 857.450, 796.609, 192.055, 204.092),
    "cluster.nr50.19.prefilter": (6, 10_235_527, 2_868_261_896, 0.805, 1366.864, 1463.418, 241.588, 246.920),
    "cluster.nr50.20.align": (4, 10_235_527, 3_225_662_877, 0.0, 854.567, 770.781, 185.898, 198.289),
    "cluster.nr50.members": (2, 10_091_519, 519_001_993, 0.0, 411.789, 405.570, 194.813, 210.757),
    "cluster.nr90.01.rescorediagonal": (1, 57_780_176, 4_040_497_318, 0.0, 22.102, 18.424, 18.424, 18.424),
    "cluster.nr90.06.rescorediagonal": (1, 37_975_640, 1_892_665_885, 0.0, 15.179, 11.165, 11.165, 11.165),
    "cluster.nr90.07.align": (2, 37_975_640, 1_519_158_769, 0.0, 386.267, 351.672, 172.026, 179.645),
    "cluster.nr90.11.prefilter": (11, 36_380_963, 11_254_104_592, 127.847, 1153.040, 2314.891, 204.757, 221.745),
    "cluster.nr90.12.align": (17, 36_380_963, 9_811_385_661, 0.0, 3947.101, 3977.816, 228.887, 237.195),
    "cluster.nr90.members": (1, 35_475_280, 544_385_754, 0.0, 182.122, 180.476, 180.476, 180.476),
}

# Measured non-split clustering steps: step: (cpu_min, threads_or_None)
RM2_OTHERS = {
    "cluster.plan": (11.410, None),
    "cluster.nr90.00.kmermatcher": (28.782, 32),
    "cluster.nr50.00.kmermatcher": (29.305, 48),
    "cluster.nr50.15.prefilter.cut": (32.207, 48),
    "cluster.nr50.19.prefilter.cut": (39.469, 48),
    "cluster": (80.439, 48),
}

R26_11_1_OTHERS = {
    "cluster.plan": (10.920, None),
    "cluster.nr90.00.kmermatcher": (30.520, 18),
    "cluster.nr50.00.kmermatcher": (19.093, 18),
    "cluster.nr50.15.prefilter.cut": (34.457, 18),
    "cluster.nr50.19.prefilter.cut": (42.162, 18),
    "cluster": (75.077, 18),
}

R26_11_2_OTHERS = {
    "cluster.plan": (13.040, None),
    "cluster.nr90.00.kmermatcher": (35.564, 18),
    "cluster.nr50.00.kmermatcher": (20.354, 18),
    "cluster.nr50.15.prefilter.cut": (36.452, 18),
    "cluster.nr50.19.prefilter.cut": (44.629, 18),
    "cluster": (94.653, 18),
}


def _write_stage_stubs(root: Path, stage: str) -> None:
    ps_rel = f"l0seq/protein_seq/gene_call_version={GCV}/stage={stage}/protein_prefix=00/part-0.parquet"
    fp_rel = f"l1/feature_protein/gene_call_version={GCV}/stage={stage}/part-0.parquet"
    for rel, col in ((ps_rel, "protein_uid"), (fp_rel, "feature_uid")):
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.table({col: [b"\x00" * 16]}), p)
    rec = {
        "stage": stage, "step": "genecall", "exit_code": 0,
        "outputs": [{"path": ps_rel, "bytes": (root / ps_rel).stat().st_size},
                    {"path": fp_rel, "bytes": (root / fp_rel).stat().st_size}],
    }
    rd = root / "runs" / f"stage={stage}"
    rd.mkdir(parents=True, exist_ok=True)
    (rd / "genecall.json").write_text(json.dumps(rec))


def _write_release(
    root: Path,
    release: str,
    stages: list[str],
    proteins: int,
    cuts: dict[str, tuple],
    others: dict[str, tuple[float, int | None]],
    *,
    recompute_b54_cuts: bool = False,
) -> None:
    catalog.define_release(root, release, stages)
    pdir = root / "l2" / "cluster_plan" / f"release={release}"
    rdir = root / "runs" / f"stage={release}"
    pdir.mkdir(parents=True, exist_ok=True)
    rdir.mkdir(parents=True, exist_ok=True)
    (pdir / "cluster.plan.json").write_text(json.dumps({
        "format": 1, "release": release, "stages": stages, "proteins": proteins,
        "cut": {"part_cpu_min": 240.0, "one_process_cpu_min_max": 20000.0},
    }))
    for call, (k, total_keys, total_units, oh, stored_w, meas_sum, lo_m, hi_m) in cuts.items():
        work = l2.call_cpu_min(call, proteins, total_units, targets=total_keys) if recompute_b54_cuts else stored_w
        parts = []
        mid_m = (meas_sum - lo_m - hi_m) / (k - 2) if k > 2 else 0.0
        for p in range(k):
            pk = total_keys // k + (total_keys % k if p == k - 1 else 0)
            pu = total_units // k + (total_units % k if p == k - 1 else 0)
            share = pu / total_units if total_units else 1.0 / k
            parts.append({"first": p, "last": p, "keys": pk, "units": pu, "cpu_min": round(work * share + oh, 6)})
            if k == 1:
                m = meas_sum
            elif p == 0:
                m = lo_m
            elif p == k - 1:
                m = hi_m
            else:
                m = mid_m
            step_name = f"{call}.part-{p:05d}" if k > 1 else call
            (rdir / f"{step_name}.json").write_text(json.dumps({
                "stage": release, "step": step_name, "exit_code": 0,
                "cpu_min": round(m, 3), "wall_s": 100.0, "peak_rss_bytes": 10_000_000_000,
            }))
        (pdir / f"{call}.cut.json").write_text(json.dumps({
            "format": 1, "release": release, "step": f"{call}.cut",
            "cpu_min": round(work, 6), "overhead_cpu_min": oh, "parts": parts,
        }))
    for step, (cpu_min, threads) in others.items():
        argv = ["mmseqs", "kmermatcher", "--threads", str(threads)] if threads is not None else []
        (rdir / f"{step}.json").write_text(json.dumps({
            "stage": release, "step": step, "exit_code": 0,
            "cpu_min": cpu_min, "wall_s": 100.0, "peak_rss_bytes": 10_000_000_000, "argv": argv,
        }))


def _run_before_and_rates(
    tmp_path: Path,
    monkeypatch,
    capsys,
    *,
    earlier: str,
    earlier_stages: list[str],
    earlier_proteins: int,
    release: str,
    release_stages: list[str],
    release_proteins: int,
    root: Path,
    cpus: int = 18,
    parts: int | None = None,
) -> tuple[dict, dict]:
    if shutil.which("mmseqs") is None:
        monkeypatch.setattr(provenance, "tool_versions", lambda *names: {n: "stub" for n in names})

    # Stub duckdb.connect so distinct(stages) returns the exact measured release protein counts
    import duckdb

    counts = iter([release_proteins, earlier_proteins])

    class _FakeCon:
        def execute(self, _sql):
            val = next(counts)
            return type("_Res", (), {"fetchone": lambda self: (val,)})()

        def close(self):
            pass

    monkeypatch.setattr(duckdb, "connect", lambda *a, **kw: _FakeCon())

    te = _te()
    ns = argparse.Namespace(
        root=str(root), release=release, stages=",".join(release_stages),
        earlier=earlier, parts=parts, profile=str(PROFILE),
        jobs=2, cpus=cpus, tmp=str(tmp_path / "tmp"), duck_threads=1,
    )
    rc = te.before(ns)
    assert rc == 0
    before_doc = json.loads(capsys.readouterr().out)
    est_path = tmp_path / f"before_{release}.json"
    est_path.write_text(json.dumps(before_doc))
    rates_doc = te.check_rates(root, release, est_path)
    return before_doc, rates_doc


def test_b64_r26_11_1_earlier_to_r26_11_2_all_judged_steps_within_20_percent(tmp_path, monkeypatch, capsys):
    """Legacy pre-B-54 `earlier = r26.11.1` calibrated for post-B-54 `r26.11.2`: every judged clustering
    step and part (>= 24 CPU-min, DEC-43) of `r26.11.2` must lie within 20% (0.80 <= ratio <= 1.20)."""
    root = tmp_path / "root"
    for s in ("S1", "C5"):
        _write_stage_stubs(root, s)
    _write_release(root, "r26.11.1", ["S1"], 49_102_041, R26_11_1_CUTS, R26_11_1_OTHERS, recompute_b54_cuts=False)
    _write_release(root, "r26.11.2", ["S1", "C5"], 57_780_176, R26_11_2_CUTS, R26_11_2_OTHERS, recompute_b54_cuts=False)

    before_doc, rates_doc = _run_before_and_rates(
        tmp_path, monkeypatch, capsys,
        earlier="r26.11.1", earlier_stages=["S1"], earlier_proteins=49_102_041,
        release="r26.11.2", release_stages=["S1", "C5"], release_proteins=57_780_176,
        root=root, cpus=18,
    )
    assert rates_doc["exit_code"] == 0, rates_doc["over_20pct"]
    assert rates_doc["over_20pct"] == [], rates_doc["over_20pct"]
    assert rates_doc["recorded_over_20pct"] == [], rates_doc["recorded_over_20pct"]
    assert rates_doc["under_20pct"] == [], rates_doc["under_20pct"]
    assert rates_doc["outside_20pct"] == [], rates_doc["outside_20pct"]

    # Unjudged rescorediagonal calls must not be inflated above 24 CPU-min by 24-part earlier factors
    factors = before_doc["projection"]["factors"]
    for name in ("cluster.nr50.01.rescorediagonal", "cluster.nr90.01.rescorediagonal", "cluster.nr90.06.rescorediagonal"):
        assert factors[name] == 1.0, (name, factors[name])

    # Pre-launch projection for target-scaled prefilters must respect part_cpu_min=240 rather than 24 fixed parts
    by_call = {c["call"]: c for c in before_doc["projection"]["calls"]}
    assert by_call["cluster.nr50.11.prefilter"]["projected"] < 400.0
    assert by_call["cluster.nr90.11.prefilter"]["projected"] < 3500.0


def test_b64_rM2_earlier_to_r26_11_1_recomputed_b54_all_judged_steps_within_20_percent(tmp_path, monkeypatch, capsys):
    """Legacy `earlier = rM.2` (with 48-thread nr50 kmermatcher) calibrated for `r26.11.1` (recomputed under B-54,
    18-thread kmermatcher): every judged clustering step and part of `r26.11.1` must lie within 20%."""
    root = tmp_path / "root"
    for s in ("S1", "H1"):
        _write_stage_stubs(root, s)
    _write_release(root, "rM.2", ["S1"], 45_815_843, RM2_CUTS, RM2_OTHERS, recompute_b54_cuts=False)
    _write_release(root, "r26.11.1", ["S1", "H1"], 49_102_041, R26_11_1_CUTS, R26_11_1_OTHERS, recompute_b54_cuts=True)

    before_doc, rates_doc = _run_before_and_rates(
        tmp_path, monkeypatch, capsys,
        earlier="rM.2", earlier_stages=["S1"], earlier_proteins=45_815_843,
        release="r26.11.1", release_stages=["S1", "H1"], release_proteins=49_102_041,
        root=root, cpus=18,
    )
    assert rates_doc["exit_code"] == 0, rates_doc["over_20pct"]
    assert rates_doc["over_20pct"] == [], rates_doc["over_20pct"]
    assert rates_doc["recorded_over_20pct"] == [], rates_doc["recorded_over_20pct"]
    assert rates_doc["under_20pct"] == [], rates_doc["under_20pct"]
    assert rates_doc["outside_20pct"] == [], rates_doc["outside_20pct"]
    # 48t -> 18t kmermatcher thread scaling keeps nr50.00.kmermatcher accurate (~20.1 vs 19.093)
    assert before_doc["projection"]["steps"]["cluster.nr50.00.kmermatcher"] < 24.0


def test_b64_r26_11_2_earlier_to_next_release_all_judged_steps_within_20_percent(tmp_path, monkeypatch, capsys):
    """Post-B-54 `earlier = r26.11.2` calibrated for a post-B-54 next release: all judged steps within 20%."""
    root = tmp_path / "root"
    for s in ("S1", "C6"):
        _write_stage_stubs(root, s)
    _write_release(root, "r26.11.2", ["S1"], 57_780_176, R26_11_2_CUTS, R26_11_2_OTHERS, recompute_b54_cuts=False)
    _write_release(root, "r26.11.3", ["S1", "C6"], 57_780_176, R26_11_2_CUTS, R26_11_2_OTHERS, recompute_b54_cuts=False)

    before_doc, rates_doc = _run_before_and_rates(
        tmp_path, monkeypatch, capsys,
        earlier="r26.11.2", earlier_stages=["S1"], earlier_proteins=57_780_176,
        release="r26.11.3", release_stages=["S1", "C6"], release_proteins=57_780_176,
        root=root, cpus=18,
    )
    assert rates_doc["exit_code"] == 0
    assert rates_doc["outside_20pct"] == [], rates_doc["outside_20pct"]
    assert all(v == 1.0 for v in before_doc["projection"]["factors"].values())


def test_b64_systematic_drift_beyond_deadband_is_still_calibrated(tmp_path, monkeypatch, capsys):
    """If a future release genuinely drifts >15% from B-54 (e.g. 1.35x on nr50.20.align), calibration tracks it."""
    root = tmp_path / "root"
    for s in ("S1", "C6"):
        _write_stage_stubs(root, s)
    # Construct earlier and next with nr50.20.align systematically 1.35x its B-54 estimate (854.567 -> 1153.665)
    drifted = dict(R26_11_2_CUTS)
    k, t, u, oh, w, _, _, _ = drifted["cluster.nr50.20.align"]
    drifted_sum = round(w * 1.35, 3)
    drifted["cluster.nr50.20.align"] = (k, t, u, oh, w, drifted_sum, round(drifted_sum / k, 3), round(drifted_sum / k, 3))

    _write_release(root, "r26.11.2", ["S1"], 57_780_176, drifted, R26_11_2_OTHERS, recompute_b54_cuts=False)
    _write_release(root, "r26.11.3", ["S1", "C6"], 57_780_176, drifted, R26_11_2_OTHERS, recompute_b54_cuts=False)

    before_doc, rates_doc = _run_before_and_rates(
        tmp_path, monkeypatch, capsys,
        earlier="r26.11.2", earlier_stages=["S1"], earlier_proteins=57_780_176,
        release="r26.11.3", release_stages=["S1", "C6"], release_proteins=57_780_176,
        root=root, cpus=18,
    )
    assert before_doc["projection"]["factors"]["cluster.nr50.20.align"] == pytest.approx(1.35, rel=1e-3)
    assert rates_doc["exit_code"] == 0
    assert rates_doc["outside_20pct"] == [], rates_doc["outside_20pct"]
