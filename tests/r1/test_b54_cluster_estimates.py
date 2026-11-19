"""Ring R1: refitted intermediate cascade part-cut estimates by surviving target DB size and kmermatcher thread-scaling
estimates (card B-54, DEC-53, DEC-55; AGENTS §6).

The oracles are measured, never the formulas under test:

- the published cluster plan cut records (`l2/cluster_plan/release=<R>/*.cut.json`) and run records
  (`runs/stage=<R>/cluster.*.json`) of `rN.3` (`10,620,684` proteins), `rM.1` (`45,815,843`), `rM.2` (`45,815,843`),
  and `r26.11.1` (`49,102,041`) under `/tmp/if2453/r10/root` (read-only on wilkcomp02, 2026-11-19), plus `rC.1b`
  (`2,053,075` proteins, stages `S1,P1`, the same stages as `rN.1b`, from `tests/r1/c22_cluster_calls_measured.tsv`);
- the eight `kmermatcher` run records (`cluster.nr90.00.kmermatcher` and `cluster.nr50.00.kmermatcher`) across `rN.3`
  (16 threads), `rM.1` (32 threads), `rM.2` (32 and 48 threads), and `r26.11.1` (18 threads).

As in `test_b52_ingest_estimate.py` and `test_b53_step_estimates.py`, the guard judges every step estimated at 24
CPU-minutes or more (DEC-43) by |record - estimate| / estimate <= 0.20 (AGENTS §6).
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from biobox import ledger, provenance
from biobox.l2 import cluster

pytestmark = pytest.mark.covers("B-54")  # ring R1 (by directory)

# Per-release, per-call measured totals for all 12 substantive split calls (the four tiny rescorediagonal calls are
# 0.35 to 2.63 CPU-min per part, never judged under DEC-43):
#   (release, call, proteins, targets_or_keys, total_units, n_parts, shared_index, measured_cpu_min)
# Read from /tmp/if2453/r10/root on wilkcomp02 (2026-11-19) and c22_cluster_calls_measured.tsv (`rC.1b`).
CASCADE_CALLS = [
    # rC.1b (2,053,075 proteins; S1,P1 = rN.1b stages; 1 part, pre-B-43 rebuild)
    ("rC.1b", "cluster.nr90.07.align", 2_053_075, 1_444_105, 30_164_641, 1, False, 4.476),
    ("rC.1b", "cluster.nr90.11.prefilter", 2_053_075, 1_406_903, 393_789_271, 1, False, 8.506),
    ("rC.1b", "cluster.nr90.12.align", 2_053_075, 1_406_903, 238_594_813, 1, False, 91.564),
    ("rC.1b", "cluster.nr90.members", 2_053_075, 1_395_712, 16_709_202, 1, False, 5.532),
    ("rC.1b", "cluster.nr50.07.align", 2_053_075, 1_048_925, 23_786_783, 1, False, 6.802),
    ("rC.1b", "cluster.nr50.11.prefilter", 2_053_075, 921_613, 249_281_723, 1, False, 5.162),
    ("rC.1b", "cluster.nr50.12.align", 2_053_075, 921_613, 109_549_428, 1, False, 39.446),
    ("rC.1b", "cluster.nr50.15.prefilter", 2_053_075, 784_957, 210_104_758, 1, False, 10.351),
    ("rC.1b", "cluster.nr50.16.align", 2_053_075, 784_957, 158_337_507, 1, False, 31.677),
    ("rC.1b", "cluster.nr50.19.prefilter", 2_053_075, 767_512, 206_770_502, 1, False, 20.458),
    ("rC.1b", "cluster.nr50.20.align", 2_053_075, 767_512, 202_861_569, 1, False, 35.396),
    ("rC.1b", "cluster.nr50.members", 2_053_075, 761_833, 16_075_323, 1, False, 10.666),
    # rN.3 (10,620,684 proteins; S1,P1,P2; pre-B-43 rebuild)
    ("rN.3", "cluster.nr90.07.align", 10_620_684, 8_546_571, 259_562_849, 1, False, 52.975),
    ("rN.3", "cluster.nr90.11.prefilter", 10_620_684, 8_363_309, 2_639_280_400, 1, False, 98.543),
    ("rN.3", "cluster.nr90.12.align", 10_620_684, 8_363_309, 2_137_797_441, 3, False, 784.400),
    ("rN.3", "cluster.nr90.members", 10_620_684, 8_278_537, 92_753_583, 1, False, 28.459),
    ("rN.3", "cluster.nr50.07.align", 10_620_684, 5_196_591, 167_430_531, 1, False, 60.396),
    ("rN.3", "cluster.nr50.11.prefilter", 10_620_684, 4_210_724, 1_266_193_392, 1, False, 34.138),
    ("rN.3", "cluster.nr50.12.align", 10_620_684, 4_210_724, 961_138_640, 1, False, 313.837),
    ("rN.3", "cluster.nr50.15.prefilter", 10_620_684, 3_141_911, 926_728_559, 1, False, 72.870),
    ("rN.3", "cluster.nr50.16.align", 10_620_684, 3_141_911, 869_122_135, 1, False, 210.883),
    ("rN.3", "cluster.nr50.19.prefilter", 10_620_684, 2_999_008, 893_814_095, 2, False, 166.864),
    ("rN.3", "cluster.nr50.20.align", 10_620_684, 2_999_008, 883_440_719, 1, False, 219.008),
    ("rN.3", "cluster.nr50.members", 10_620_684, 2_966_200, 87_441_246, 1, False, 67.594),
    # rM.1 (45,815,843 proteins; S1,P1..P6; 24 parts, shared index on nr50.15 and nr50.19)
    ("rM.1", "cluster.nr90.07.align", 45_815_843, 33_085_999, 1_356_134_142, 24, False, 326.44),
    ("rM.1", "cluster.nr90.11.prefilter", 45_815_843, 31_920_514, 10_194_781_853, 24, False, 3583.39),
    ("rM.1", "cluster.nr90.12.align", 45_815_843, 31_920_514, 8_616_200_048, 24, False, 3505.47),
    ("rM.1", "cluster.nr90.members", 45_815_843, 31_183_876, 432_415_353, 24, False, 148.77),
    ("rM.1", "cluster.nr50.07.align", 45_815_843, 17_348_804, 698_116_275, 24, False, 275.46),
    ("rM.1", "cluster.nr50.11.prefilter", 45_815_843, 13_520_716, 4_070_572_296, 24, False, 1287.42),
    ("rM.1", "cluster.nr50.12.align", 45_815_843, 13_520_716, 3_038_129_206, 24, False, 1253.90),
    ("rM.1", "cluster.nr50.15.prefilter", 45_815_843, 8_915_865, 2_601_223_781, 24, True, 402.44),
    ("rM.1", "cluster.nr50.16.align", 45_815_843, 8_915_865, 2_717_279_156, 24, False, 698.64),
    ("rM.1", "cluster.nr50.19.prefilter", 45_815_843, 8_312_546, 2_454_500_276, 24, True, 900.53),
    ("rM.1", "cluster.nr50.20.align", 45_815_843, 8_312_546, 2_611_996_035, 24, False, 679.65),
    ("rM.1", "cluster.nr50.members", 45_815_843, 8_192_065, 409_423_542, 24, False, 332.39),
    # rM.2 (45,815,843 proteins; S1,P1..P6; 24 parts, shared index on nr50.15 and nr50.19)
    ("rM.2", "cluster.nr90.07.align", 45_815_843, 33_085_999, 1_356_134_142, 24, False, 340.23),
    ("rM.2", "cluster.nr90.11.prefilter", 45_815_843, 31_920_514, 10_194_781_853, 24, False, 3598.11),
    ("rM.2", "cluster.nr90.12.align", 45_815_843, 31_920_514, 8_616_200_048, 24, False, 3500.65),
    ("rM.2", "cluster.nr90.members", 45_815_843, 31_183_876, 432_415_353, 24, False, 148.46),
    ("rM.2", "cluster.nr50.07.align", 45_815_843, 17_348_804, 698_116_275, 24, False, 279.40),
    ("rM.2", "cluster.nr50.11.prefilter", 45_815_843, 13_520_716, 4_070_572_296, 24, False, 1307.65),
    ("rM.2", "cluster.nr50.12.align", 45_815_843, 13_520_716, 3_038_129_206, 24, False, 1259.10),
    ("rM.2", "cluster.nr50.15.prefilter", 45_815_843, 8_915_865, 2_601_223_781, 24, True, 421.93),
    ("rM.2", "cluster.nr50.16.align", 45_815_843, 8_915_865, 2_717_279_156, 24, False, 707.66),
    ("rM.2", "cluster.nr50.19.prefilter", 45_815_843, 8_312_546, 2_454_500_276, 24, True, 957.96),
    ("rM.2", "cluster.nr50.20.align", 45_815_843, 8_312_546, 2_611_996_035, 24, False, 685.88),
    ("rM.2", "cluster.nr50.members", 45_815_843, 8_192_065, 409_423_542, 24, False, 332.32),
    # r26.11.1 (49,102,041 proteins; S1,P1..P6,H1,PL1; 24 parts, shared index on nr50.15 and nr50.19)
    ("r26.11.1", "cluster.nr90.07.align", 49_102_041, 35_153_162, 1_412_500_763, 24, False, 340.41),
    ("r26.11.1", "cluster.nr90.11.prefilter", 49_102_041, 33_882_143, 10_664_275_036, 24, False, 3832.62),
    ("r26.11.1", "cluster.nr90.12.align", 49_102_041, 33_882_143, 9_125_396_776, 24, False, 3656.60),
    ("r26.11.1", "cluster.nr90.members", 49_102_041, 33_117_773, 463_925_032, 24, False, 162.48),
    ("r26.11.1", "cluster.nr50.07.align", 49_102_041, 18_876_326, 758_111_957, 24, False, 289.75),
    ("r26.11.1", "cluster.nr50.11.prefilter", 49_102_041, 14_912_484, 4_385_758_277, 24, False, 1433.59),
    ("r26.11.1", "cluster.nr50.12.align", 49_102_041, 14_912_484, 3_333_344_123, 24, False, 1337.46),
    ("r26.11.1", "cluster.nr50.15.prefilter", 49_102_041, 10_064_063, 2_852_219_962, 24, True, 466.31),
    ("r26.11.1", "cluster.nr50.16.align", 49_102_041, 10_064_063, 3_073_734_693, 24, False, 755.21),
    ("r26.11.1", "cluster.nr50.19.prefilter", 49_102_041, 9_421_654, 2_698_588_280, 24, True, 1373.84),
    ("r26.11.1", "cluster.nr50.20.align", 49_102_041, 9_421_654, 2_964_472_779, 24, False, 743.50),
    ("r26.11.1", "cluster.nr50.members", 49_102_041, 9_291_312, 440_098_571, 24, False, 356.51),
]

# The DEC-55 trigger parts of r26.11.1's cluster.nr50.19.prefilter (part-00000, part-00001, and worst part-00010):
# (part_idx, part_units, measured_cpu_min) out of 2,698,588,280 total units and 9,421,654 surviving targets.
R26_NR50_19_PARTS = [
    (0, 112_441_227, 56.741),
    (1, 112_441_240, 55.825),
    (10, 112_441_069, 60.154),
]

# All 8 measured kmermatcher records across rN.3, rM.1, rM.2, and r26.11.1:
# (release, step, input_sequences, threads, measured_cpu_min)
KMERMATCHER_RECORDS = [
    ("rN.3", "cluster.nr90.00.kmermatcher", 10_620_684, 16, 6.3899),
    ("rN.3", "cluster.nr50.00.kmermatcher", 8_278_537, 16, 4.6200),
    ("rM.1", "cluster.nr90.00.kmermatcher", 45_815_843, 32, 29.1044),
    ("rM.1", "cluster.nr50.00.kmermatcher", 31_183_876, 32, 19.3616),
    ("rM.2", "cluster.nr90.00.kmermatcher", 45_815_843, 32, 28.7817),
    ("rM.2", "cluster.nr50.00.kmermatcher", 31_183_876, 48, 29.3050),
    ("r26.11.1", "cluster.nr90.00.kmermatcher", 49_102_041, 18, 30.5202),
    ("r26.11.1", "cluster.nr50.00.kmermatcher", 33_117_773, 18, 19.0934),
]

BY = {c.name: c for c in cluster.calls()}


@pytest.mark.parametrize(
    "release,call,proteins,targets,units,n_parts,shared,meas",
    CASCADE_CALLS,
    ids=[f"{r}:{c}" for r, c, *_ in CASCADE_CALLS],
)
def test_b54_each_cascade_call_is_within_20_percent_of_its_target_scaled_estimate(
    release, call, proteins, targets, units, n_parts, shared, meas
):
    c = BY[call]
    work = cluster.call_cpu_min(call, proteins, units, targets=targets)
    oh = cluster.part_overhead_cpu_min(c, targets, shared=shared) if c.module == "prefilter" else 0.0
    est = work + n_parts * oh
    assert abs(meas - est) <= 0.20 * est, (release, call, round(est, 2), meas, f"{(meas - est) / est:+.1%}")


@pytest.mark.parametrize("part_idx,part_units,meas", R26_NR50_19_PARTS, ids=["part-00000", "part-00001", "part-00010"])
def test_b54_r26_nr50_19_prefilter_dec55_parts_are_within_20_percent(part_idx, part_units, meas):
    call = "cluster.nr50.19.prefilter"
    proteins, targets, total_units = 49_102_041, 9_421_654, 2_698_588_280
    work = cluster.call_cpu_min(call, proteins, total_units, targets=targets)
    oh = cluster.part_overhead_cpu_min(BY[call], targets, shared=True)
    est = work * (part_units / total_units) + oh
    assert abs(meas - est) <= 0.20 * est, (part_idx, round(est, 3), meas, f"{(meas - est) / est:+.1%}")


@pytest.mark.parametrize(
    "release,step,seqs,threads,meas",
    KMERMATCHER_RECORDS,
    ids=[f"{r}:{s}:{t}t" for r, s, _, t, _ in KMERMATCHER_RECORDS],
)
def test_b54_kmermatcher_records_are_within_20_percent_including_48_vs_18_threads(release, step, seqs, threads, meas):
    est = cluster.kmermatcher_cpu_min(seqs, threads=threads)
    assert abs(meas - est) <= 0.20 * est, (release, step, threads, round(est, 3), meas, f"{(meas - est) / est:+.1%}")


def test_b54_cut_doc_scales_intermediate_prefilter_work_by_surviving_targets(tmp_path):
    """At fixed release `proteins` and fixed query `units`, doubling the surviving target count in `input_step2.index`
    must scale `_cut_doc`'s `cpu_min` by `2 ** g`, not leave it flat."""
    name = "cluster.nr50.19.prefilter"
    c = BY[name]
    units = 2_000_000_000
    cuts = []
    for label, n_targets in (("small", 1_000), ("double", 2_000)):
        d = tmp_path / label
        src = d / f"{c.positional[0]}.index"
        src.parent.mkdir(parents=True, exist_ok=True)
        src.write_text("".join(f"{k}\t{k * 10}\t{units // n_targets + 2}\n" for k in range(n_targets)))
        ctx = cluster._Ctx(root=d, out=None, release="rX", plan={"proteins": 49_102_041, "cut": {"parts": 24}},
                           plan_sha256="", w=d, threads=1)
        cuts.append(cluster._cut_doc(ctx, c))
    _, g = cluster.CALL_FIT[name]
    assert g > 0.5
    assert cuts[1]["cpu_min"] == pytest.approx(cuts[0]["cpu_min"] * (2.0 ** g), rel=1e-5)


def test_b54_call_cpu_min_without_targets_falls_back_to_release_proteins():
    """Callers that project before `.cut` runs (`c35_est.py`, `train_est.py`) pass three positional arguments."""
    name = "cluster.nr50.19.prefilter"
    rate, growth = cluster.CALL_FIT[name]
    assert cluster.call_cpu_min(name, cluster.CALL_N0, 10**9) == pytest.approx(rate)
    assert cluster.call_cpu_min(name, 2 * cluster.CALL_N0, 10**9) == pytest.approx(rate * (2.0 ** growth))


def test_b54_legacy_pre_b54_cut_records_remain_reusable_under_b49(tmp_path):
    """Published releases clustered before card B-54 (`rM.2`, `r26.11.1`) wrote `.cut.json` with the legacy
    `proteins`-only `cpu_min`. `_source_problem` must still accept them while rejecting an arbitrary wrong `cpu_min`."""
    root = tmp_path / "root"
    src, rel = "rM.2", "r26.11.2"
    proteins = 45_815_843
    tools = provenance.tool_versions("mmseqs2", "pyarrow")
    plan_doc = {
        "format": cluster.PLAN_FORMAT, "release": src, "stages": ["S1"], "gene_call_version": "gv",
        "mmseqs2": tools["mmseqs2"], "files": [], "cut": {"parts": 1}, "proteins": proteins,
        "release_cpu_min": 1.0, "outputs": {"db": "a", "db_h": "b", "prot.parquet": "c"},
    }
    plan_dir = root / cluster.plan_dir(src)
    plan_dir.mkdir(parents=True)
    raw_plan = cluster._json_bytes(plan_doc)
    (plan_dir / f"{cluster.PLAN_STEP}.json").write_bytes(raw_plan)
    sha1 = hashlib.sha256(raw_plan).hexdigest()

    for c in cluster.calls():
        if c.split:
            units = 2_454_500_276 if c.name == "cluster.nr50.19.prefilter" else 100_000_000
            legacy_est = round(cluster._legacy_call_cpu_min(c.name, proteins, units), 6)
            cut_doc = {
                "format": cluster.PLAN_FORMAT, "release": src, "step": f"{c.name}.cut", "plan_sha256": sha1,
                "module": c.module, "inputs": {}, "split": c.positional[cluster.PER_QUERY[c.module]],
                "units": "residues" if c.module == "prefilter" else "entry_bytes",
                "cpu_min": legacy_est, "overhead_cpu_min": 0.0,
                "parts": [{"first": 0, "last": 99, "keys": 100, "units": units, "cpu_min": legacy_est}],
            }
            (plan_dir / f"{c.name}.cut.json").write_bytes(cluster._json_bytes(cut_doc))
        call_doc = {
            "format": cluster.PLAN_FORMAT, "release": src, "step": c.name, "plan_sha256": sha1,
            "module": c.module, "argv": list(c.argv), "inputs": {}, "outputs": {},
        }
        (plan_dir / f"{c.name}.json").write_bytes(cluster._json_bytes(call_doc))

    # Minimal published tables + ledger + runs/stage=rM.2/cluster.json
    outputs = []
    for d in cluster._table_dirs(src)[:3]:
        p = root / d / "part-0.parquet"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"stub")
        outputs.append({"path": p.relative_to(root).as_posix(), "bytes": 4, "sha256": hashlib.sha256(b"stub").hexdigest()})
    ledger.write_release(root, src, "cluster", [
        ledger.row(b"nr50", "level", 10, 5, 5, tools=tools),
        ledger.row(b"nr90", "level", 20, 10, 10, tools=tools),
    ])
    led_p = ledger.release_path(root, src, "cluster")
    outputs.append({"path": led_p.relative_to(root).as_posix(), "bytes": led_p.stat().st_size,
                    "sha256": hashlib.sha256(led_p.read_bytes()).hexdigest()})
    run_p = root / "runs" / f"stage={src}" / "cluster.json"
    run_p.parent.mkdir(parents=True, exist_ok=True)
    run_p.write_text(json.dumps({"exit_code": 0, "stage": src, "step": "cluster", "outputs": outputs}))

    want = {k: plan_doc[k] for k in cluster._INPUT_KEYS if k != "release"}
    assert cluster._source_problem(root, rel, src, want) is None

    # Tampering with cpu_min to a value matching neither new nor legacy model must still be refused
    bad_path = plan_dir / "cluster.nr50.19.prefilter.cut.json"
    bad_doc = json.loads(bad_path.read_text())
    bad_doc["cpu_min"] = 12345.0
    bad_path.write_bytes(cluster._json_bytes(bad_doc))
    why = cluster._source_problem(root, rel, src, want)
    assert why is not None and "another model" in why
