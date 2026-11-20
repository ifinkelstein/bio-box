"""Ring R2, card B-55: buffered `crispr._Scan` scratch Parquet row groups, 100x synthetic store bounded memory
(AGENTS §4.2), and byte-identity of published `feature`, `l1_crispr`, and `ledger` outputs for both `crispr.run` and
`crispr.run_near_cas` (including shared contigs and circular-origin arrays).
"""

from __future__ import annotations

import contextlib
import hashlib
import io as stdio
import json
import math
import runpy
import shutil
import subprocess
import sys
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from biobox import ids, io, runstep
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l1 import crispr

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "r1"))
import test_crispr_plant as plant  # noqa: E402
import test_e15_near_cas_producer as e15  # noqa: E402

MINCED = shutil.which("minced")
pytestmark = [pytest.mark.r2, pytest.mark.covers("B-55", "E-12", "E-15"),
              pytest.mark.skipif(MINCED is None, reason="minced (MinCED) not on PATH; it is in the "
                                                        "bio-box image, not in GitHub Actions")]

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts/stress/synth_scale.py"
SMALL, BIG = 10, 100
GROWTH_BOUND = 32 * 2**20
PEAK_BOUND = 768 * 2**20


def _digests(root: Path) -> dict[str, str]:
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(Path(root).rglob("*")) if p.is_file()}


def _assert_scratch_row_groups(paths: dict[str, Path]) -> dict[str, tuple[int, int]]:
    out: dict[str, tuple[int, int]] = {}
    for name, path in paths.items():
        pf = pq.ParquetFile(path)
        n_rows = pf.metadata.num_rows
        n_rg = pf.num_row_groups
        bound = max(1, math.ceil(n_rows / io.ROW_GROUP_ROWS))
        assert n_rg <= bound, f"{name}.parquet: {n_rg} row groups for {n_rows} rows (bound {bound})"
        out[name] = (n_rows, n_rg)
    return out


def test_real_minced_scan_scratch_row_groups(tmp_path, monkeypatch):
    """Run real MinCED through both `crispr.run` and `crispr.run_near_cas` with a small `batch_bases` (forcing
    multiple batches) and assert every scratch Parquet file obeys the `io.ROW_GROUP_ROWS` row-group bound."""
    layout = e15.cases()
    root = e15.make_stage(tmp_path / "stage", "Q", layout)
    observed: list[dict[str, tuple[int, int]]] = []

    real_publish = crispr._publish
    real_publish_nc = crispr._publish_near_cas

    def spy_publish(out, stage, prod, sc, scratch, tools, max_file_bytes):
        observed.append(_assert_scratch_row_groups(sc.paths))
        return real_publish(out, stage, prod, sc, scratch, tools, max_file_bytes)

    def spy_publish_nc(out, stage, prod, sc, scratch, tools, max_file_bytes, inputs, families):
        observed.append(_assert_scratch_row_groups(sc.paths))
        return real_publish_nc(out, stage, prod, sc, scratch, tools, max_file_bytes, inputs, families)

    monkeypatch.setattr(crispr, "_publish", spy_publish)
    monkeypatch.setattr(crispr, "_publish_near_cas", spy_publish_nc)

    out_run = tmp_path / "out_run"
    out_nc = tmp_path / "out_nc"
    scratch = tmp_path / "scratch"
    scratch.mkdir()

    # Run default producer into a fresh out dir (root already has default producer from make_stage, so use a new stage)
    root2 = tmp_path / "stage2"
    contigs = {k: v[0] for k, v in layout.items()}
    ingest(e15.write_objects(tmp_path / "obj2", contigs), root2, "R")
    write_genome_input(root2, "R", e15.genomes_of(contigs))
    crispr.run(root2, "R", 2, tmp=scratch, out=out_run, batch_bases=10_000)
    crispr.run_near_cas(root, "Q", 2, library_version=e15.LIB_VERSION, tmp=scratch, out=out_nc, batch_bases=10_000)

    assert len(observed) == 2
    assert set(observed[0]) == {"arrays", "elements", "contigs"}
    assert set(observed[1]) == {"arrays", "elements", "contigs", "topology"}
    assert observed[0]["contigs"][0] > len(layout)  # base rows + batch hit rows, still 1 row group
    assert observed[0]["contigs"][1] == 1
    assert observed[1]["topology"][1] == 1


def test_published_byte_identity_run_and_near_cas_with_shared_and_circular_contigs(tmp_path):
    """Published `feature`, `l1_crispr`, and `ledger` files of `crispr.run` and `crispr.run_near_cas` are
    byte-identical across thread counts and batch sizes, including circular-origin arrays and shared contigs
    (DEC-26 another-stage owner and DEC-27 shared across two genomes in the same stage)."""
    layout = e15.cases()
    # Add an earlier stage S0 that owns one contig (`near`), so stage Q skips it (DEC-26),
    # and put `across_origin` in two genomes of stage Q so it is shared within the stage (DEC-27).
    root = tmp_path / "root"
    contigs = {k: v[0] for k, v in layout.items()}
    ingest(e15.write_objects(tmp_path / "obj_s0", {"near": contigs["near"]}), root, "S0")
    write_genome_input(root, "S0", e15.genomes_of({"near": contigs["near"]}))

    ingest(e15.write_objects(tmp_path / "obj_q", contigs), root, "Q")
    genomes = e15.genomes_of(contigs)
    across_cu = ids.contig_uid(contigs["across_origin"].seq)
    plain_cu = ids.contig_uid(contigs["gap_1000_down"].seq)
    genomes.append(GenomeInput((across_cu, plain_cu), 11, kind="plasmid", circular=(across_cu,)))
    write_genome_input(root, "Q", genomes)

    fp, dom = e15.gene_rows({k: (v[0], v[1]) for k, v in layout.items()})
    e15.write_rows(root, "feature_protein", f"l1/feature_protein/gene_call_version={e15.V}/stage=Q", fp)
    e15.write_rows(root, "l1_domains", f"{e15.LIB_DIR}/stage=Q", dom)

    scratch = tmp_path / "scratch"
    scratch.mkdir()

    # 1. Default producer: cpus=1, default batch_bases vs cpus=2, tiny batch_bases=8_000
    out_a = tmp_path / "run_a"
    out_b = tmp_path / "run_b"
    rep_a = crispr.run(root, "Q", 1, tmp=scratch, out=out_a, seq_stages=["S0"])
    rep_b = crispr.run(root, "Q", 2, tmp=scratch, out=out_b, seq_stages=["S0"], batch_bases=8_000)
    assert rep_a["skipped_contigs"] == rep_b["skipped_contigs"]
    assert len(rep_a["skipped_contigs"]) >= 2  # at least S0-owned + DEC-27 shared
    assert _digests(out_a) == _digests(out_b)
    assert len(_digests(out_a)) == 3  # feature, l1_crispr, ledger

    # Publish default producer into root so run_near_cas can read it
    crispr.run(root, "Q", 1, tmp=scratch, seq_stages=["S0"])

    # 2. Filtered near_cas producer: cpus=1, default batch_bases vs cpus=2, tiny batch_bases=8_000
    nc_a = tmp_path / "nc_a"
    nc_b = tmp_path / "nc_b"
    nrep_a = crispr.run_near_cas(root, "Q", 1, library_version=e15.LIB_VERSION, tmp=scratch, out=nc_a,
                                 seq_stages=["S0"])
    nrep_b = crispr.run_near_cas(root, "Q", 2, library_version=e15.LIB_VERSION, tmp=scratch, out=nc_b,
                                 seq_stages=["S0"], batch_bases=8_000)
    assert nrep_a["kept"] == nrep_b["kept"] > 0
    assert nrep_a["skipped_contigs"] == nrep_b["skipped_contigs"]
    assert _digests(nc_a) == _digests(nc_b)
    assert len(_digests(nc_a)) == 3  # feature, l1_crispr, ledger


# sha256 of the three files the CRISPR step publishes on the `micro_s0` fixture (default producer, 2 MinCED
# processes, stage `S0`; fixture built through ingest, gene calls, IGRs and tapes as `test_crispr_producer.py`'s
# `fixture_stage` does), taken from the code BEFORE card B-55: commit 044586d97c194e96b1a73678d170d6db685df807
# (`git rev-parse 5647ce7^1`, the integration tree just before B-55 was merged), measured 2026-11-20 on wilkcomp02
# in the test image, by a throwaway pytest file built on this module's `_build_micro_s0_crispr` steps and run with
# `scripts/pod_pytest.sh D-07-old <file> -s` on a worktree of that commit; two builds of that tree gave the same
# three values. They record the old writer, not the function under test.
OLD_WRITER_PRODUCER = "minced-0.4.2-w8-nr3-r23x47-s26x50"
OLD_WRITER_SHA256 = {
    f"l1/feature/producer={OLD_WRITER_PRODUCER}/stage=S0/part-0.parquet":
        "9e5abc0fe516e24a1b496c5293a4f72d8c1dc82b64ce1799e3275fc97c183cd1",
    f"l1/crispr/producer={OLD_WRITER_PRODUCER}/stage=S0/part-0.parquet":
        "9e6dce0e270751638e53bed7748f49986fe26572bd7357a6a64af51dfb5dc706",
    f"ledger/stage=S0/step=crispr-{OLD_WRITER_PRODUCER}/part-0.parquet":
        "857cd28ffda2ee2c783c79880f790a293e1984d39a371c3b4a27308db8c26fa3",
}


def _build_micro_s0_crispr(base: Path) -> dict[str, str]:
    import test_builders_out as c08
    from biobox.l1 import igr
    from biobox.l1.genecall import run as call_genes
    from biobox.l3 import tapes

    root, scratch = base / "root", base / "scratch"
    scratch.mkdir(parents=True)
    contract = c08.CONTRACT
    ingest([SourceObject(c08.FIX / o["file"], contract["source_name"], contract["source_release"],
                         contract["license"], True) for o in contract["objects"]], root, "S0")
    genomes = c08.genome_inputs("S0") + c08.genome_inputs("S1")
    write_genome_input(root, "S0", genomes)
    call_genes(root, "S0", genomes, tmp=scratch)
    igr.run(root, "S0", tmp=scratch)
    tapes.run(root, "S0", tmp=scratch)
    before = set(_digests(root))
    crispr.run(root, "S0", 2, tmp=scratch)
    return {k: v for k, v in _digests(root).items() if k not in before}


def test_micro_s0_published_crispr_files_equal_the_pre_b55_writers_sha256(tmp_path):
    """Item (3) of card D-07: the `feature`, `l1_crispr` and `ledger` files of the fixture's CRISPR step are
    byte-identical to those the code before card B-55 published (pinned above)."""
    assert _build_micro_s0_crispr(tmp_path) == OLD_WRITER_SHA256


CHILD = r"""
import json, math, sys
from pathlib import Path
import pyarrow.parquet as pq
from biobox import io
from biobox.l1 import crispr

seen = {}
real_pub = crispr._publish
def spy_pub(out, stage, prod, sc, scratch, tools, max_file_bytes):
    for k, p in sc.paths.items():
        pf = pq.ParquetFile(p)
        seen[k] = [pf.metadata.num_rows, pf.num_row_groups]
    return real_pub(out, stage, prod, sc, scratch, tools, max_file_bytes)
crispr._publish = spy_pub

res = crispr.run(Path(sys.argv[1]), "S0", 1, tmp=Path(sys.argv[2]), out=Path(sys.argv[3]))
print(json.dumps({"res": {k: res[k] for k in ("contigs", "bases", "batches", "arrays")}, "scratch": seen}))
"""


def _store(out: Path, copies: int) -> Path:
    argv = sys.argv
    sys.argv = [str(SCRIPT), str(REPO), str(out), str(copies), "--shard-bases", str(2**26), "--seed", "1"]
    try:
        with contextlib.redirect_stdout(stdio.StringIO()):
            g = runpy.run_path(str(SCRIPT), run_name="__main__")
    finally:
        sys.argv = argv
    g["con"].close()
    return Path(g["root"])


def _measured(root: Path, tmp: Path) -> tuple[int, str, dict]:
    peak = tmp / f"peak-{root.parent.name}.json"
    scratch, out = tmp / f"scratch-{root.parent.name}", tmp / f"out-{root.parent.name}"
    scratch.mkdir()
    cmd = runstep.measured([sys.executable, "-c", CHILD, str(root), str(scratch), str(out)], peak, False)
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    assert p.returncode == 0, p.stderr[-4000:]
    rss, source = runstep.step_peak(peak)
    assert rss is not None, source
    return rss, source, json.loads(p.stdout.strip().splitlines()[-1])


def test_100x_synth_scale_memory_and_scratch_row_groups(tmp_path, record_property):
    """AGENTS §4.2: 100x synthetic store peak RSS stays bounded (growth <= 32 MiB over 10x, peak <= 768 MiB)
    and every scratch Parquet file has at most `max(1, ceil(num_rows / io.ROW_GROUP_ROWS))` row groups."""
    stores = {n: _store(tmp_path / f"b55x{n}", n) for n in (SMALL, BIG)}
    got = {n: _measured(root, tmp_path) for n, root in stores.items()}
    record_property("b55_peak_bytes", {n: g[0] for n, g in got.items()})
    record_property("b55_runs", {n: g[2] for n, g in got.items()})
    for n, (_, _, payload) in got.items():
        for name, (n_rows, n_rg) in payload["scratch"].items():
            bound = max(1, math.ceil(n_rows / io.ROW_GROUP_ROWS))
            assert n_rg <= bound, f"{n}x {name}.parquet: {n_rg} row groups for {n_rows} rows (bound {bound})"
    assert got[BIG][2]["res"]["bases"] >= 9 * got[SMALL][2]["res"]["bases"] > 0
    growth = got[BIG][0] - got[SMALL][0]
    assert growth <= GROWTH_BOUND, (f"peak RSS {got[SMALL][0] / 2**20:.1f} MiB at {SMALL}x, {got[BIG][0] / 2**20:.1f} "
                                    f"MiB at {BIG}x: grew {growth / 2**20:.1f} MiB, bound {GROWTH_BOUND / 2**20:.0f}")
    assert max(g[0] for g in got.values()) <= PEAK_BOUND, got
