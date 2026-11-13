"""Ring R2, card E-07: `build-stage --ncrna-cm` fills the `features` slot with the ncRNA producer.

On the C-04 fixture (lambda and M. genitalium, `test_builders_out.direct`, whose direct build is the build without
ncRNA, as `test_build_stage` shows), with the fixture library and a model set of two whole-genome and five
intergenic models (`test_ncrna_regions.MIXED`):

- a build with `--ncrna-cm` and one without give two tape versions, and the gene-call tables are byte-identical
  between them; the IGR, tape, token and posting steps are given both producers; its ncRNA files are the direct
  producer's on the same stage;
- with a budget below the ncRNA estimate the scan runs as plan, shard and merge steps, and the build's files are
  those of the build that scans in one process, byte for byte (shard records apart);
- the shard pool's memory term: shards at once are bounded by the memory cap at the counted bytes per cmsearch, and
  a cap below one cmsearch is refused before the scan.
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import pytest

from biobox import cli, runstep
from biobox import stage as bs
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l3 import tapes

sys.path.insert(0, str(Path(__file__).parent))
import test_build_stage as c04  # noqa: E402
import test_builders_out as c08  # noqa: E402
import test_ncrna_regions as rg  # noqa: E402

direct = c08.direct  # the module-scoped direct build of S0 and release c08.REL (a fixture)

pytestmark = [pytest.mark.r2, pytest.mark.covers("E-07", "D-01", "MR-20"), rg.needs_infernal,
              pytest.mark.skipif(c08.MMSEQS is None, reason="build-stage clusters with MMseqs2 (in the image)")]

GENE_CALL_DIRS = (f"l1/feature/producer={V}/", f"l0seq/protein_seq/gene_call_version={V}/",
                  f"l1/feature_protein/gene_call_version={V}/", f"l1/protein_occurrence/gene_call_version={V}/")
SHARD_RECORDS = ("l1/ncrna_scan/",)  # card E-08: the `ncrna_scan` artifact (card E-07 had them under domain_scan)


def set_file(direct) -> Path:
    p = direct["work"].base / "ncrna-mixed.tsv"
    if not p.exists():
        rg.write_set(p, rg.MIXED)
    return p


def with_ncrna(direct, root: Path) -> bs.Build:
    return dataclasses.replace(c04.build_of(direct, root), ncrna_cm=rg.LIB, ncrna_set=str(set_file(direct)))


def producer(direct) -> str:
    from biobox.l1 import ncrna

    return ncrna.producer(rg.LIB, set_file(direct))


@pytest.fixture(scope="module")
def with_cm(direct, tmp_path_factory):
    root = tmp_path_factory.mktemp("e07-build") / "root"
    argv = c04.argv_of(direct, root) + ["--ncrna-cm", str(rg.LIB), "--ncrna-set", str(set_file(direct))]
    assert cli.main(argv) == 0
    return root


def test_the_steps_with_ncrna_fill_the_features_slot_and_pass_both_producers(direct, tmp_path):
    prof = direct["prof"]
    p = producer(direct)
    got = bs.steps(prof, with_ncrna(direct, tmp_path / "root"))
    names = [s.name for s in got]
    assert names == [n if n != "features" else f"ncrna-{p}" for n in c04.NAMES]
    nc = next(s for s in got if s.name == f"ncrna-{p}")
    assert nc.check == "ncrna" and nc.argv[:2] == ("stage", "ncrna") and "--cm" in nc.argv
    for s in got:
        both = s.argv is not None and "--producers" in s.argv
        assert both == (s.name in ("igr", "tapes", "tokens", "postings")), s.name
        if both:
            assert s.argv[s.argv.index("--producers") + 1] == f"{V},{p}"
    assert got == bs.steps(prof, with_ncrna(direct, tmp_path / "root"))  # a pure function of its arguments


def test_with_and_without_ncrna_two_tape_versions_and_the_same_gene_call_tables(direct, with_cm, tmp_path):
    p = producer(direct)
    without, withc = c08.tree(direct["root"]), c08.tree(with_cm)
    gene = {k: v for k, v in without.items() if k.startswith(GENE_CALL_DIRS)}
    assert gene and {k: v for k, v in withc.items() if k.startswith(GENE_CALL_DIRS)} == gene
    tv_one, tv_two = tapes.tape_version([V]), tapes.tape_version([V, p])
    assert tv_one != tv_two
    assert {k.split("/")[2] for k in without if k.startswith("l3/tapes/")} == {f"tape_version={tv_one}"}
    assert {k.split("/")[2] for k in withc if k.startswith("l3/tapes/")} == {f"tape_version={tv_two}"}
    assert any(k.startswith(f"l1/ncrna/producer={p}/stage=S0/") for k in withc)
    assert (with_cm / "runs/stage=S0" / f"ncrna-{p}.json").is_file()
    # its ncRNA files are the direct producer's on the same stage
    from biobox.l1 import ncrna

    out = tmp_path / "out"
    ncrna.run(direct["root"], "S0", rg.LIB, 2, model_set=set_file(direct), tmp=tmp_path, out=out)
    want = {k: v for k, v in c08.tree(out).items()}
    assert want and {k: withc[k] for k in want} == want


def test_an_ncrna_scan_over_the_budget_runs_in_shards_whose_files_equal_the_one_process_build(
        direct, with_cm, tmp_path, monkeypatch):
    p = producer(direct)
    monkeypatch.setattr(bs, "NCRNA_CPU_H_PER_10GBP", 1e7)  # the fixture's bases far over 80% of the budget
    assert bs.ncrna_shard_bases(direct["prof"]) < 600_000
    root = tmp_path / "root"
    seen: list[str] = []
    assert bs.build_stage(direct["prof"], with_ncrna(direct, root), run=c04.runner, mem_cap=lambda: None,
                          _hook=seen.append) == 0
    names = c04.record_names(root)
    shards = sorted(n for n in names if n.startswith(f"stage=S0/ncrna-{p}.shard-"))
    assert len(shards) == 2 and f"stage=S0/ncrna-{p}.plan" in names and f"stage=S0/ncrna-{p}" in names
    assert f"after:ncrna-{p}.plan" in seen
    got = {k: v for k, v in c08.tree(root).items() if not k.startswith(SHARD_RECORDS)}
    want = c08.tree(with_cm)
    assert c08.data(got) == c08.data(want)


def test_the_memory_term_bounds_the_shards_at_once_and_refuses_a_cap_below_one(direct, tmp_path, capsys):
    gib = 2**30
    assert bs.ncrna_job_bytes(1) == bs.ncrna_job_bytes(2) == int(6.4 * gib)
    assert bs.ncrna_job_bytes(4) == int(2 * 6.4 * gib)
    assert bs.ncrna_jobs(10, 2, None) == 10
    assert bs.ncrna_jobs(10, 2, 16 * gib) == 2  # 16 GiB holds two cmsearch of 6.4 GiB
    assert bs.ncrna_jobs(10, 4, 16 * gib) == 1
    assert bs.ncrna_jobs(10, 2, 6 * gib) == 0
    # card E-08 (E-07 report (d)): shards are cut at the slower of the two measured shard rates, 304.0 CPU-hours
    # per 10 Gbp (derived in planning/measurements/2026-11-11_ncrna_producer.md), not the design note's 247.4
    # (card E-09: the rate is the set's; this is prok-1's, unchanged)
    assert bs.ncrna_shard_bases(direct["prof"], "prok-1") == int(0.8 * 300 * 1e10 / (304.0 * 60)) == 131_578_947
    root = tmp_path / "root"
    assert bs.build_stage(direct["prof"], with_ncrna(direct, root), run=c04.runner,
                          mem_cap=lambda: 4 * gib) == runstep.EXIT_REFUSED
    assert "above the memory cap" in capsys.readouterr().err
    assert not (root / "l1/ncrna").exists()
    assert not any(n.startswith("stage=S0/ncrna-") for n in c04.record_names(root))


@pytest.mark.covers("E-08")
def test_a_shard_planned_at_80_percent_stays_under_the_budget_at_the_slower_measured_rate(direct):
    """E-07 report (d): the two measured shards ran at 206.6 and 304.0 CPU-hours per 10 Gbp (derived, E-07's pod
    check). A shard cut at 80% of the 300 CPU-minute budget at the slower rate uses 240 CPU-minutes at that rate,
    and at the faster one less; at the design note's 247.4 it would have used 295 at the slower rate."""
    prof = direct["prof"]
    bases = bs.ncrna_shard_bases(prof, "prok-1")  # E-07's measurement was of prok-1 (card E-09: the rate is the set's)
    assert bs.NCRNA_CPU_H_PER_10GBP == 304.0
    for rate in (206.6, 254.5, 304.0):  # shard 0, both shards, shard 1 (E-07's measurement)
        assert bases * rate * 60 / 1e10 <= 0.8 * prof.cpu_min_per_process + 1e-6
    assert bs.ncrna_cpu_min(bases, "prok-1") == pytest.approx(240.0, abs=1e-4)


@pytest.mark.covers("C-08", "E-08")
def test_check_provenance_passes_on_a_release_with_ncrna_tapes(direct, with_cm):
    """E-07 report (a): a release whose tapes carry the ncRNA producer passes C-08. The fixture library is a subset
    of Rfam 15.1, so its version is `sha` + 12 hex digits; with that as the lock's `rfam` the release passes, and with
    the repository's lock (Rfam 15.1) the only problems are that version's."""
    from biobox import catalog, provenance
    from biobox.l1 import ncrna

    p = producer(direct)
    man = catalog.build_manifest(with_cm, c04.REL, ["S0"])
    defs = [f for f in man["files"] if f["table"] == "tape_version_def"]
    assert defs and any(f["partitions"]["tape_version"] == tapes.tape_version([V, p]) for f in defs)
    assert any(f["table"] == "ledger" and f["partitions"]["step"] == f"ncrna-{p}" for f in man["files"])
    lock = {**provenance.locked_versions(), "rfam": ncrna.library_version(rg.LIB)}
    assert provenance.check_provenance(man, lock) == []
    probs = provenance.check_provenance(man)
    assert probs and all("rfam" in x for x in probs), probs


def test_a_missing_set_is_refused_before_any_step(direct, tmp_path, capsys):
    root = tmp_path / "root"
    b = dataclasses.replace(with_ncrna(direct, root), ncrna_set="no-such-set")
    assert bs.build_stage(direct["prof"], b, run=c04.runner) == runstep.EXIT_REFUSED
    assert "no-such-set" in capsys.readouterr().err
    assert not (root / "runs").exists()
