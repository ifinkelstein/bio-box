"""Ring R1: the random-placement null of card E-16 (`scripts/stress/e16_crispr_near_cas_phages.py`) on constructed
genomes whose answer is counted by hand.

The null places a call of L bases uniformly over the starts at which it fits in its own genome (on a linear contig of
n bases, n - L + 1 starts; on a circular one, n starts, the call running across the origin when s + L > n) and keeps
it when the filtered producer's three rules hold at that place (card E-15, `biobox.l1.crispr` docstring): (a) it
shares no base with a default array, (b) its gap to a Cas gene is at most 1,000 bases, (c) it is not wholly inside a
called gene. The share of starts kept is written down here from the geometry (every interval is in the comments), so
the oracle is the hand count, not the script (AGENTS §5.2). The script's exact share, its Monte Carlo estimate (within
four standard errors) and its reading of E-13's window share are compared with it; E-13's own `placement_share` is
the second reading of the window share.
"""

from __future__ import annotations

import importlib.util
import itertools
import sys
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.covers("E-16")  # ring R1 by its directory (tests/conftest.py)

STRESS = Path(__file__).resolve().parents[2] / "scripts" / "stress"
L = 100
CAS4 = "PF01930"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, STRESS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod  # its dataclasses look their module up by name
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def e16():
    return _load("e16_crispr_near_cas_phages")


def _linear_genome(e16):
    """Contig A (100,000 bases, linear): a non-Cas gene [48000, 49500), a Cas4 gene [50000, 51000), a default array
    [51500, 51600). Contig B (50,000 bases, linear): nothing. Starts of a 100-base call:
    - domain: A 0..99900 (99,901 starts), B 0..49900 (49,901): 149,802;
    - (b) gap to the Cas gene at most 1,000: s + 100 >= 49000 and s <= 52000, so 48900..52000: 3,101 starts;
    - (c) inside the Cas gene 50000..50900 (901), inside the other gene 48000..49400, of which 48900..49400 (501) in
      the window;
    - (a) sharing a base with the default array: s + 100 > 51500 and s < 51600, so 51401..51599 (199);
    kept: 3,101 - 901 - 501 - 199 = 1,500 of 149,802."""
    a = e16.make_contig(100_000, False, genes=[(48_000, 49_500), (50_000, 51_000)], cas={1: {CAS4}},
                        defaults=[(51_500, 51_600)])
    b = e16.make_contig(50_000, False)
    return [a, b]


def _circular_genome(e16):
    """One circular contig of 10,000 bases with a Cas4 gene [100, 600). Starts of a 100-base call, 0..9999 (10,000):
    window of the gene as stored, -1000..1600, so 0..1600 (1,601); of its copy at +n, [10100, 10600):
    9000..11600, so 9000..9999 (1,000); 2,601 in all; inside the gene 100..500 (401); kept 2,200 of 10,000.
    E-13's window share reads every contig as linear: starts 0..9900 (9,901), window 0..1600 and 9000..9900 (901):
    2,502 of 9,901."""
    return [e16.make_contig(10_000, True, genes=[(100, 600)], cas={0: {CAS4}})]


def test_the_exact_share_kept_is_the_hand_count(e16):
    r = e16.exact_keep(_linear_genome(e16), L)
    assert (r["keep"], r["domain"]) == (1500, 149_802)
    assert r["family_keep"] == {CAS4: 1500}
    r = e16.exact_keep(_circular_genome(e16), L)
    assert (r["keep"], r["domain"]) == (2200, 10_000)


def test_e13s_window_share_is_the_hand_count_and_e13s_own(e16):
    e13 = _load("e13_crispr_min_repeats")
    assert e16.e13_window(_linear_genome(e16), L) == (3101, 149_802)
    assert e16.e13_window(_circular_genome(e16), L) == (2502, 9901)
    # E-13's own function on the same geometry
    starts, ends = np.array([48_000, 50_000]), np.array([49_500, 51_000])
    w10 = e13.Windows(starts, ends, 100_000, False, np.array([1]))
    w1 = e13.Windows(starts, ends, 100_000, False, np.array([1]), k=1)
    assert e13.placement_share(L, [(100_000, w10, w1), (50_000, None, None)])["1_kb"] == pytest.approx(3101 / 149_802)
    s, e = np.array([100]), np.array([600])
    w10 = e13.Windows(s, e, 10_000, True, np.array([0]))
    w1 = e13.Windows(s, e, 10_000, True, np.array([0]), k=1)
    assert e13.placement_share(L, [(10_000, w10, w1)])["1_kb"] == pytest.approx(2502 / 9901)


@pytest.mark.parametrize("where, verdict", [
    ((49_450, 49_550), "ok"),          # 450 bases before the Cas gene, past the other gene's end
    ((50_100, 50_200), "inside"),      # wholly inside the Cas gene
    ((51_550, 51_650), "overlap"),     # shares bases with the default array
    ((52_000, 52_100), "ok"),          # gap 1,000: within
    ((52_001, 52_101), "far"),         # gap 1,001
    ((48_900, 49_000), "inside"),      # gap 1,000, but wholly inside the other gene
    ((10_000, 10_100), "far"),         # contig A's far end
])
def test_the_rules_at_a_place_by_hand(e16, where, verdict):
    a = _linear_genome(e16)[0]
    ov, near, ins = e16.rules_at(a, np.array([where[0]]), np.array([where[1]]))
    got = "overlap" if ov[0] else "far" if not near[0] else "inside" if ins[0] else "ok"
    assert got == verdict


def test_the_rules_read_round_the_circle(e16):
    c = _circular_genome(e16)[0]
    ov, near, ins = e16.rules_at(c, np.array([9_800, 9_950, 5_000]), np.array([9_900, 10_050, 5_100]))
    assert near.tolist() == [True, True, False]  # 200 and 50 bases before the gene's copy at +n; far
    assert not ov.any() and not ins.any()


@pytest.mark.parametrize("genome, p", [("linear", 1500 / 149_802), ("circular", 2200 / 10_000)])
def test_the_monte_carlo_estimate_is_within_four_standard_errors_and_seeded(e16, genome, p):
    g = _linear_genome(e16) if genome == "linear" else _circular_genome(e16)
    draws = 200_000
    kept = e16.simulate(g, L, draws, np.random.default_rng(16))["kept"]
    se = (p * (1 - p) / draws) ** 0.5
    assert abs(kept.mean() - p) < 4 * se
    again = e16.simulate(g, L, draws, np.random.default_rng(16))["kept"]
    assert np.array_equal(kept, again)


def test_the_monte_carlo_window_share_reads_e13s_linear_placement(e16):
    draws = 200_000
    near = e16.simulate(_circular_genome(e16), L, draws, np.random.default_rng(3), linear=True)["near"]
    p = 2502 / 9901
    assert abs(near.mean() - p) < 4 * (p * (1 - p) / draws) ** 0.5


def test_the_poisson_binomial_tail_by_enumeration(e16):
    ps = [0.1, 0.5, 0.3, 0.05]
    for k in range(6):
        want = sum(np.prod([q if b else 1 - q for q, b in zip(ps, bits)])
                   for bits in itertools.product((0, 1), repeat=len(ps)) if sum(bits) >= k)
        assert e16.tail_at_least(ps, k) == pytest.approx(want, abs=1e-12)


# --- card E-17: a stage whose genome_input lists a contig that another stage stores (DEC-26) ----------------------

CAS_PROT, GENE_PROT = b"\x01" * 16, b"\x02" * 16
UID_A, UID_B, UID_C = b"\xa0" * 16, b"\xb0" * 16, b"\xc0" * 16  # C: listed by stage H, stored by stage S
GENOME, CALL_UID = b"\x77" * 16, b"\x55" * 16


def _put(path: Path, cols: dict, schema: dict) -> None:
    """One part file at `path` / part-0.parquet: `cols` {name: list} with the types of `schema` {name: pa type}."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    path.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.table({k: pa.array(cols.get(k, []), schema[k]) for k in schema}), path / "part-0.parquet")


def _stage_on_disk(e16, tmp_path: Path, listed_extra: bool, other_default: bool = True):
    """Stage H holds the genome of `_linear_genome` (contigs A and B) with one call of 100 bases on A at [60000, 60100),
    the default array [51500, 51600), the Cas4 gene [50000, 51000) and the other gene [48000, 49500). With
    `listed_extra`, genome_input lists a third contig C of 70,000 bases that stage S stores (`other_default`: as its
    default variant) and that H's contig_index does not hold."""
    import pyarrow as pa

    b16 = pa.binary(16)
    root, calls, kept = tmp_path / "root", tmp_path / "calls", tmp_path / "kept"
    listed = [(UID_A, False), (UID_B, False)] + ([(UID_C, True)] if listed_extra else [])
    _put(root / "l0seq/genome_input/stage=H",
         {"contig_uid": [u for u, _ in listed], "genome_uid": [GENOME] * len(listed), "is_circular": [c for _, c in listed]},
         {"contig_uid": b16, "genome_uid": b16, "is_circular": pa.bool_()})
    idx = {"contig_uid": b16, "contig_len": pa.int64(), "is_default_variant": pa.bool_()}
    _put(root / "l0seq/contig_index/stage=H", {"contig_uid": [UID_A, UID_B], "contig_len": [100_000, 50_000],
                                               "is_default_variant": [True, True]}, idx)
    _put(root / "l0seq/contig_index/stage=S", {"contig_uid": [UID_C], "contig_len": [70_000],
                                               "is_default_variant": [other_default]}, idx)
    feat = {"feature_uid": b16, "contig_uid": b16, "start": pa.int64(), "end": pa.int64()}
    one = {"feature_uid": [CALL_UID], "contig_uid": [UID_A], "start": [60_000], "end": [60_100]}
    arr = {"feature_uid": [b"\x66" * 16], "contig_uid": [UID_A], "start": [51_500], "end": [51_600]}
    _put(calls / f"l1/feature/producer={e16.NR2}/stage=H", one, feat)
    _put(calls / f"l1/feature/producer={e16.NR3}/stage=H", arr, feat)
    _put(kept / f"l1/feature/producer={e16.NEAR}/stage=H", {}, feat)
    _put(calls / f"l1/crispr/producer={e16.NR2}/stage=H", {"feature_uid": [CALL_UID], "contig_len": [100_000]},
         {"feature_uid": b16, "contig_len": pa.int64()})
    _put(root / f"l1/domains/library={e16.LIBRARY}/stage=H",
         {"protein_uid": [CAS_PROT], "model_acc": [f"{CAS4}.22"], "model_name": ["Cas4"]},
         {"protein_uid": b16, "model_acc": pa.string(), "model_name": pa.string()})
    _put(root / f"l1/feature_protein/gene_call_version={e16.GENE_CALL}/stage=H",
         {"contig_uid": [UID_A, UID_A], "start": [48_000, 50_000], "end": [49_500, 51_000], "strand": [1, 1],
          "protein_uid": [GENE_PROT, CAS_PROT]},
         {"contig_uid": b16, "start": pa.int64(), "end": pa.int64(), "strand": pa.int8(), "protein_uid": b16})
    return root, calls, kept


@pytest.mark.covers("E-17")
@pytest.mark.parametrize("listed_extra", [False, True])
def test_a_stage_that_lists_a_contig_another_stage_stores_keeps_the_hand_count(e16, tmp_path, listed_extra):
    """E-17. The null of stage H on the hand-counted genome: 1,500 of 149,802 starts kept, whether or not
    genome_input also lists contig C, which stage S stores (it is left out and counted, with its stage)."""
    root, calls, kept = _stage_on_disk(e16, tmp_path, listed_extra)
    fams = {CAS4: "cas4"}
    stg = e16.load_stage(root, calls, kept, "H", fams, seq_stages=["S"])
    assert set(stg.contigs) == {UID_A, UID_B}
    res, rows, _ = e16.measure_stage(stg, 20, 16, fams, _load("e13_crispr_min_repeats"))
    assert (rows[0]["keep_starts"], rows[0]["domain_starts"]) == (1500, 149_802)
    assert res["nulls"]["N2"]["expected_exact"] == round(1500 / 149_802, 4)
    assert res["contigs"] == 2 and res["bases"] == 150_000 and res["contigs_listed"] == 2 + listed_extra
    assert res["contigs_left_out"] == int(listed_extra)
    assert res["contigs_left_out_by_stage"] == ({"S": 1} if listed_extra else {})


@pytest.mark.covers("E-17")
def test_a_listed_contig_that_no_named_stage_stores_is_still_an_error(e16, tmp_path):
    root, calls, kept = _stage_on_disk(e16, tmp_path, True)
    with pytest.raises(SystemExit, match="stored by none of the stages"):
        e16.load_stage(root, calls, kept, "H", {CAS4: "cas4"})  # stage S not named
    with pytest.raises(SystemExit, match="stored by none of the stages"):
        e16.load_stage(root, calls, kept, "H", {CAS4: "cas4"}, seq_stages=["H", "T"])  # T has no contig_index


@pytest.mark.covers("E-17")
def test_a_non_default_variant_in_another_stage_does_not_store_the_contig(e16, tmp_path):
    root, calls, kept = _stage_on_disk(e16, tmp_path, True, other_default=False)
    with pytest.raises(SystemExit, match="stored by none of the stages"):
        e16.load_stage(root, calls, kept, "H", {CAS4: "cas4"}, seq_stages=["S"])
