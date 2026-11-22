"""Ring R1, card B-66: O(1)-stage SeqStore point lookups and per-stored-contig CRISPR CPU-minute term.

Tests:
1. `KeyLookup.find_first` probes files in order and stops at the first matching row without decoding later files'
   row groups.
2. On a 10-stage synthetic store (`S0`..`S9`), `SeqStore(root, stages).get(cu)` for a contig stored in `S0` reads
   `store.index.row_groups_read == 1` and `store.locator.row_groups_read == 1` (instead of all 10 stages), still
   resolves contigs stored in `S9`, raises the same `KeyError` on missing contigs/digests, and matches forward,
   reverse-complement (`strand=-1`), and `circular=True` origin-wrapping slices byte-for-byte.
3. `stage.crispr_cpu_min(bases, flank_bases=0, *, contigs=0, near_cas=False)` keeps `contigs=0` bit-identical to
   card B-53 (`CRISPR_CPU_S_PER_GBP = 168.1`, `CRISPR_NEAR_CAS_COST = 1.15` unchanged) and adds the per-stored-contig
   term when `contigs > 0`; `stage.stage_contigs(root, stage)` counts default-variant contigs from `ingest.json`'s
   explicit `contig_index` files one row group at a time; `stage.crispr_refusal` and `stage_est.py` (`before` and
   `rates`) pass stored contigs to `crispr_cpu_min`.
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import json
import random
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import fasta, ids, profile
from biobox import stage as bs
from biobox.l0seq import SeqStore, SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l0seq.lookup import KeyLookup
from biobox.l0seq.store import revcomp

pytestmark = pytest.mark.covers("B-66")

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"


def _seq(rng: random.Random, n: int) -> bytes:
    return rng.randbytes(n).translate(bytes(b"ACGT"[i % 4] for i in range(256)))


def _stage_est():
    spec = importlib.util.spec_from_file_location("stage_est", BSP / "stage_est.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["stage_est"] = mod
    spec.loader.exec_module(mod)
    return mod


def _ten_stage_store(tmp_path: Path) -> tuple[Path, list[str], dict[str, list[tuple[bytes, bytes]]]]:
    """Build a 10-stage store (`S0`..`S9`), each stage holding 4 distinct contigs so every stage's
    `contig_index` and `seq_locator` span hash space and have 1 row group each."""
    rng = random.Random(6601)
    root = tmp_path / "root"
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    stages = [f"S{i}" for i in range(10)]
    by_stage: dict[str, list[tuple[bytes, bytes]]] = {}
    for st in stages:
        seqs = [_seq(rng, 400 + 37 * k) for k in range(4)]
        fna = in_dir / f"{st}.fna"
        fna.write_bytes(b"\n".join(ln for k, s in enumerate(seqs) for ln in fasta.format_record(f"{st}_c{k}", s, 80)) + b"\n")
        ingest([SourceObject(fna, "b66", "1", "CC0", True)], root, st)
        by_stage[st] = [(ids.contig_uid(s), s) for s in seqs]
    return root, stages, by_stage


def test_b66_keylookup_find_first_stops_at_first_matching_file(tmp_path):
    """`KeyLookup.find_first` decodes only up to the first file with a matching row, leaving later files unread."""
    files = []
    target_key = b"\x80" * 16
    for i in range(5):
        p = tmp_path / f"part-{i}.parquet"
        keys = [b"\x00" * 16, target_key, b"\xff" * 16]
        flags = [False, i == 1, False]
        pq.write_table(
            pa.table({
                "k": pa.array(keys, pa.binary(16)),
                "flag": pa.array(flags, pa.bool_()),
                "src": pa.array([i, i, i], pa.int64()),
            }),
            p,
        )
        files.append(p)

    lk = KeyLookup(files, "k", ["flag", "src"])
    # Without predicate: matches in file 0, reads 1 row group
    row0 = lk.find_first(target_key)
    assert row0 is not None and row0["src"] == 0 and row0["flag"] is False
    assert lk.row_groups_read == 1

    # Fresh lookup with predicate flag=True: reads file 0 (no match) and file 1 (match), stops before files 2..4
    lk2 = KeyLookup(files, "k", ["flag", "src"])
    row1 = lk2.find_first(target_key, lambda r: r["flag"])
    assert row1 is not None and row1["src"] == 1 and row1["flag"] is True
    assert lk2.row_groups_read == 2

    # Missing key returns None
    assert lk2.find_first(b"\x77" * 16) is None


def test_b66_seqstore_get_is_o1_stage_and_single_locate(tmp_path):
    """On a 10-stage store (`S0`..`S9`), `SeqStore(root, stages).get(cu)` for a contig in `S0` reads exactly
    1 `contig_index` row group and 1 `seq_locator` row group, still resolves `S9`, raises the same `KeyError`
    on missing contigs/digests, and matches forward, reverse-complement, and circular origin-wrapping slices."""
    root, stages, by_stage = _ten_stage_store(tmp_path)
    cu0, seq0 = by_stage["S0"][0]
    cu9, seq9 = by_stage["S9"][0]

    store = SeqStore(root, stages)
    try:
        # 1. Contig in S0: O(1) stage lookup + single locate
        got0 = store.get(cu0)
        assert got0 == seq0
        assert store.index.row_groups_read == 1, f"expected 1 contig_index RG read, got {store.index.row_groups_read}"
        assert store.locator.row_groups_read == 1, f"expected 1 seq_locator RG read, got {store.locator.row_groups_read}"

        # 2. Forward slice, reverse complement, circular origin wrap
        assert store.get(cu0, 10, 80, strand=1) == seq0[10:80]
        assert store.get(cu0, 10, 80, strand=-1) == revcomp(seq0[10:80])
        n0 = len(seq0)
        wrap_start, wrap_end = n0 - 25, n0 + 35
        assert store.get(cu0, wrap_start, wrap_end, strand=1, circular=True) == seq0[wrap_start:] + seq0[:35]
        assert store.get(cu0, wrap_start, wrap_end, strand=-1, circular=True) == revcomp(seq0[wrap_start:] + seq0[:35])

        # Invalid slice / strand still raises ValueError
        with pytest.raises(ValueError):
            store.get(cu0, 0, n0 + 1, circular=False)
        with pytest.raises(ValueError):
            store.get(cu0, 0, 10, strand=0)
        with pytest.raises(ValueError):
            store.get(cu0, n0 - 10, 2 * n0, circular=True)

        # 3. Contig in last stage S9 still resolves
        assert store.get(cu9) == seq9

        # 4. Missing contig_uid and raw_digest raise identical KeyError
        missing_cu = b"\xab" * 16
        with pytest.raises(KeyError, match=f"contig_uid {missing_cu.hex()} is in none of the stages"):
            store.get(missing_cu)
        missing_rd = b"\xcd" * 32
        with pytest.raises(KeyError, match=f"raw_digest {missing_rd.hex()} is in none of the stages"):
            store.raw(missing_rd)
    finally:
        store.close()


def test_b66_crispr_cpu_min_keeps_b53_at_zero_contigs_and_adds_contig_term():
    """`contigs=0` is bit-identical to card B-53 (`CRISPR_CPU_S_PER_GBP = 168.1`, `CRISPR_NEAR_CAS_COST = 1.15`);
    `contigs > 0` adds a positive per-stored-contig term (`CRISPR_CPU_S_PER_K_CONTIGS`)."""
    assert bs.CRISPR_CPU_S_PER_GBP == 168.1
    assert bs.CRISPR_NEAR_CAS_COST == 1.15
    assert bs.CRISPR_CPU_S_PER_K_CONTIGS > 0

    bases, flank = 10_040_000_000, 123_456_789
    assert bs.crispr_cpu_min(bases, flank, contigs=0, near_cas=False) == bs.crispr_cpu_min(bases, flank, near_cas=False)
    assert bs.crispr_cpu_min(bases, flank, contigs=0, near_cas=True) == bs.crispr_cpu_min(bases, flank, near_cas=True)
    assert bs.crispr_cpu_min(bases, flank, contigs=0) == (bases + flank) * 168.1 / 1e9 / 60
    assert bs.crispr_cpu_min(bases, flank, contigs=0, near_cas=True) == (bases + flank) * (168.1 * 1.15) / 1e9 / 60

    with_c = bs.crispr_cpu_min(bases, flank, contigs=1_000_000, near_cas=False)
    without_c = bs.crispr_cpu_min(bases, flank, contigs=0, near_cas=False)
    assert with_c > without_c
    assert with_c - without_c == pytest.approx(1_000_000 * bs.CRISPR_CPU_S_PER_K_CONTIGS / 1000.0 / 60.0)


def _two_stage_crispr_root(tmp_path: Path) -> tuple[Path, Path, Path]:
    """Stage S0 owns 2 contigs; stage S1 lists those 2 (owned by S0) plus 5 new contigs stored by S1."""
    rng = random.Random(6602)
    root = tmp_path / "root"
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    s0_seqs = [_seq(rng, 5_000), _seq(rng, 6_000)]
    s1_new = [_seq(rng, 4_000 + 100 * i) for i in range(5)]

    def write_stage(st: str, seqs: list[bytes]) -> Path:
        fna = in_dir / f"{st}.fna"
        fna.write_bytes(b"\n".join(ln for k, s in enumerate(seqs) for ln in fasta.format_record(f"{st}_{k}", s, 80)) + b"\n")
        cus = [ids.contig_uid(s) for s in seqs]
        genomes = [GenomeInput((cu,), 11, kind="prokaryote", circular=()) for cu in cus]
        write_genome_input(root, st, genomes)
        ingest([SourceObject(fna, "b66", "1", "CC0", True)], root, st)
        runs = root / "runs" / f"stage={st}"
        runs.mkdir(parents=True, exist_ok=True)
        gi_files = sorted((root / "l0seq" / "genome_input" / f"stage={st}").glob("part-*.parquet"))
        (runs / "inputs.json").write_text(json.dumps({
            "outputs": [{"path": p.relative_to(root).as_posix()} for p in gi_files]
        }))
        ci_files = sorted((root / "l0seq" / "contig_index" / f"stage={st}").glob("part-*.parquet"))
        (runs / "ingest.json").write_text(json.dumps({
            "stage": st, "step": "ingest", "exit_code": 0, "cpu_min": 0.1,
            "outputs": [{"path": p.relative_to(root).as_posix()} for p in ci_files],
        }))
        (runs / "genecall.json").write_text(json.dumps({"stage": st, "step": "genecall", "exit_code": 0, "cpu_min": 0.1}))
        return fna

    write_stage("S0", s0_seqs)
    s1_fna = write_stage("S1", s0_seqs + s1_new)
    s1_in = tmp_path / "s1_inputs"
    s1_in.mkdir()
    s1_cus = [ids.contig_uid(s) for s in s0_seqs + s1_new]
    (s1_in / "genomes.json").write_text(json.dumps([
        {"contig_uids": [cu.hex()], "genetic_code": 11, "kind": "prokaryote", "circular": [], "assembly_quality": "contig"}
        for cu in s1_cus
    ]))
    (s1_in / "objects.json").write_text(json.dumps([
        {"path": str(s1_fna), "source_name": "b66", "source_release": "1", "license": "CC0", "redistributable": True}
    ]))
    return root, s1_in, s1_fna


def test_b66_stage_contigs_and_crispr_refusal_and_stage_est(tmp_path, capsys):
    """`stage_contigs` counts only default-variant contigs owned by the stage; `crispr_refusal` and
    `stage_est.py` (`before` and `rates`) pass stored contigs to `crispr_cpu_min`."""
    root, s1_in, _ = _two_stage_crispr_root(tmp_path)
    assert bs.stage_contigs(root, "S0") == 2
    assert bs.stage_contigs(root, "S1") == 5

    bases = bs.stage_bases(root, "S1")
    flank = bs.crispr_flank_bases(root, "S1")
    assert flank == 0
    est_no_c = bs.crispr_cpu_min(bases, flank, contigs=0, near_cas=False)
    est_with_c = bs.crispr_cpu_min(bases, flank, contigs=5, near_cas=False)
    assert est_with_c > est_no_c

    # Tight budget between est_no_c / 0.8 and est_with_c / 0.8: refused only when contigs are counted
    mid_budget = ((est_no_c + est_with_c) / 2.0) / bs.SHARD_MARGIN
    prof = profile.load(REPO / "profiles" / "wilkcomp02.yaml")
    tight = dataclasses.replace(prof, cpu_min_per_process=mid_budget)
    b = bs.Build(root=root, stage="S1", release="rB66.1", genomes=s1_in / "genomes.json",
                 objects=s1_in / "objects.json", stages=("S0", "S1"), crispr=True)
    step = bs.Step("S1", "crispr-minced", ("stage", "crispr"), check="crispr")
    why = bs.crispr_refusal(tight, b, root, step)
    assert why is not None and "5" in why and "no shards" in why, why

    # stage_est.py before & rates
    est = _stage_est()
    note = "planning/measurements/2026-11-20_stage_V1.md"
    conf_p = tmp_path / "S1.conf"
    conf_p.write_text(
        f"STAGE=S1\nNOTE={note}\nINPUTS={s1_in}\nSUM={'a' * 64} genomes.json\nSUM={'b' * 64} objects.json\n"
        "DEPTH=called\nCRISPR=yes\n"
    )
    prof_p = tmp_path / "prof.yaml"
    prof_p.write_text((REPO / "profiles" / "wilkcomp02.yaml").read_text())

    # before
    a = argparse.Namespace(conf=str(conf_p), stage="S1", root=str(root), profile=str(prof_p),
                           inputs=str(s1_in), tmp=str(tmp_path / "est_tmp"), out_plan=None)
    rc = est.before(a)
    assert rc == 0
    bdoc = json.loads(capsys.readouterr().out)
    assert bdoc["contigs_stored"] == 5
    assert bdoc["crispr_default_cpu_min"] == round(est_with_c, 3)

    # rates: record a crispr step and check its estimate matches crispr_cpu_min(..., contigs=5)
    from biobox.l1 import crispr as l1crispr
    runs = root / "runs" / "stage=S1"
    cname = f"crispr-minced-0.4.2-w8-nr3-r23x47-s26x50"
    (runs / f"{cname}.json").write_text(json.dumps({"stage": "S1", "step": cname, "exit_code": 0, "cpu_min": 0.5}))
    orig_prod = l1crispr.producer
    try:
        l1crispr.producer = lambda *a, **kw: "minced-0.4.2-w8-nr3-r23x47-s26x50"
        rdoc = est.rates_doc(est.read_conf(conf_p, "S1"), root, prof, s1_in / "objects.json", None, tmp_path / "rtmp")
    finally:
        l1crispr.producer = orig_prod
    crow = next(r for r in rdoc["rows"] if r["step"] == cname)
    assert crow["estimate"] == round(est_with_c, 3)


# Measured CRISPR records from /tmp/if2453/r10/root on wilkcomp02, 2026-11-22 (`b66_crispr_contig_arms.py summary`):
# stage -> (stored_bases, flank_bases, contigs_stored, default_cpu_min, near_cas_cpu_min)
COMPLETE_NON_C5_RECORDS = {
    "S1": (1_664_109_938, 0, 27_258, 5.085492, 5.486575),
    "P1": (1_016_625_803, 0, 539, 3.084777, 3.382703),
    "P2": (9_596_180_328, 0, 5_018, 26.48908, 31.794481),
    "P3": (9_949_077_522, 0, 9_107, 29.194508, 33.385812),
    "P4": (9_976_456_950, 0, 4_005, 27.585047, 32.606163),
    "P5": (9_967_711_632, 0, 4_111, 27.505385, 32.360891),
    "P6": (9_962_536_103, 0, 4_227, 27.937717, 32.796186),
    "H1": (1_819_831_338, 0, 6_714, 4.924614, 5.876108),
    "PL1": (5_924_552_875, 6_205_263_778, 61_647, 30.938485, 35.513282),
    "V1": (5_178_493_077, 1_818_968_830, 168_358, 20.560122, 21.893386),
}


@pytest.mark.parametrize("near_cas", [False, True], ids=["default", "near_cas"])
@pytest.mark.parametrize("stage", sorted(COMPLETE_NON_C5_RECORDS))
def test_b66_all_ten_non_c5_complete_stages_within_15_percent(stage, near_cas):
    """Every complete non-C5 stage (`S1..V1`, 20 records) stays within 15% of `crispr_cpu_min` with stored contigs."""
    bases, flank, contigs, def_min, nc_min = COMPLETE_NON_C5_RECORDS[stage]
    cpu_min = nc_min if near_cas else def_min
    est = bs.crispr_cpu_min(bases, flank, contigs=contigs, near_cas=near_cas)
    assert abs(cpu_min - est) <= 0.15 * est, (stage, near_cas, round(est, 3), cpu_min, round(cpu_min / est, 3))


def test_b66_post_fix_g1_and_spx1_samples_within_15_percent_and_spx1_under_240():
    """Post-fix pod scratch arms (`b66_crispr_contig_arms.py`, wilkcomp02, 2026-11-22) match `crispr_cpu_min`
    within 15%, and full-stage `SPX1` sits safely under the 240 CPU-minute margin."""
    # G1 1-in-20 sample (61,361 contigs, 510,039,671 bp, 12 stages opened)
    g1_def_est = bs.crispr_cpu_min(510_039_671, 0, contigs=61_361, near_cas=False)
    g1_nc_est = bs.crispr_cpu_min(510_039_671, 0, contigs=61_361, near_cas=True)
    assert abs(2.0179 - g1_def_est) <= 0.15 * g1_def_est, (round(g1_def_est, 4), 2.0179)
    assert abs(2.1339 - g1_nc_est) <= 0.15 * g1_nc_est, (round(g1_nc_est, 4), 2.1339)

    # SPX1 1-in-40 sample (49,153 contigs, 730,724,623 bp, 13 stages opened)
    spx_sample_est = bs.crispr_cpu_min(730_724_623, 0, contigs=49_153, near_cas=False)
    assert abs(2.3149 - spx_sample_est) <= 0.15 * spx_sample_est, (round(spx_sample_est, 4), 2.3149)

    # Full-stage SPX1 (2,829,368 contigs, 35,885,000,000 bp, linear)
    spx_full_est = bs.crispr_cpu_min(35_885_000_000, 0, contigs=2_829_368, near_cas=False)
    assert spx_full_est < 240.0, spx_full_est
