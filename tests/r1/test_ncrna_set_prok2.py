"""Card E-09 (DEC-30): the frozen ncRNA model set `prok-2` is the default; `prok-1` stays readable.

The oracles are files this card does not write: the proposal `planning/design/0-08a_prok-2_proposed.tsv` (written by
measurement 0-08a's report script), the frozen `prok-1.tsv` (card E-07) and the design note's table of Rfam 15.1's
4,227 models (`planning/design/0-08_rfam15.1_models.tsv`). The set's relative cost is DEC-30's derived 299.7 against
254.5 CPU-hours per 10 Gbp (`planning/measurements/2026-11-12_rfam_sample.md`, item 6).
"""

from __future__ import annotations

import dataclasses
import hashlib
from pathlib import Path

import pytest

from biobox import cli, profile, provenance
from biobox import stage as bs
from biobox.l1 import ncrna

pytestmark = pytest.mark.covers("E-09")  # ring R1 (by directory)

REPO = Path(__file__).resolve().parents[2]
PROPOSED = REPO / "planning/design/0-08a_prok-2_proposed.tsv"
MODELS = REPO / "planning/design/0-08_rfam15.1_models.tsv"
LIB = REPO / "tests/fixtures/micro_s0/rfam15.1_micro_s0.cm.gz"
PROFILE = REPO / "profiles/wilkcomp02.yaml"
# item 6 of the measurement: the 11 models the sample rule scans on whole genomes; 4 of them are prok-1 models
SWITCHED = {"RF01670", "RF02111", "RF02820", "RF02908"}  # P17, IS009, V_AS9, rli30: intergenic in prok-1
NEW_WHOLE = {"RF01419", "RF01497", "RF02510", "RF02647", "RF02735", "RF02812", "RF03135"}


def data_lines(path: Path) -> list[str]:
    """The header and rows of a set file: every line that is not a comment."""
    return [ln for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]


def test_prok2_is_the_proposal_with_a_frozen_header():
    path = ncrna.set_path("prok-2")
    assert path == ncrna.SET_DIR / "prok-2.tsv"
    assert data_lines(path) == data_lines(PROPOSED)  # the proposal's rows, byte for byte, in its order
    head = "".join(ln for ln in path.read_text().splitlines(keepends=True) if ln.startswith("#"))
    assert "Frozen" in head and "DEC-30" in head and "PROPOSED" not in head
    assert "planning/design/0-08a_prok-2_proposed.tsv" in head
    got = ncrna.read_set("prok-2")
    assert len(got) == 770 and len({m.acc for m in got}) == 770
    assert sum(m.whole_genome for m in got) == 34


def test_every_prok2_model_is_an_rfam_15_1_model_with_its_name():
    table = {}
    for line in MODELS.read_text().splitlines()[1:]:
        f = line.split("\t")
        table[f[0]] = f[1]
    assert len(table) == 4227
    missing = [(m.acc, m.name) for m in ncrna.read_set("prok-2") if table.get(m.acc) != m.name]
    assert missing == []


def test_prok1_rows_are_in_prok2_with_their_flags_unchanged_but_the_four_the_sample_rule_moves():
    """Card E-09 as corrected 2026-11-12: every prok-1 row is in prok-2 with its name and flag, except four models
    the sample rule (more than half of their 0-08a hits outside every margin-200 region) moves to whole genome."""
    one = {m.acc: m for m in ncrna.read_set("prok-1")}
    two = {m.acc: m for m in ncrna.read_set("prok-2")}
    assert len(one) == 588
    want = {a: dataclasses.replace(m, whole_genome=True) if a in SWITCHED else m for a, m in one.items()}
    assert {a: two.get(a) for a in one} == want  # same name; the flag unchanged but for the four
    assert all(not one[a].whole_genome for a in SWITCHED)
    added = [m for a, m in two.items() if a not in one]
    assert len(added) == 182
    assert {m.acc for m in added if m.whole_genome} == NEW_WHOLE
    # DEC-30: 34 whole genome = prok-1's 23 + the 4 moved + 7 new
    assert sum(m.whole_genome for m in one.values()) == 23 and len(SWITCHED) == 4 and len(NEW_WHOLE) == 7


def test_prok2_is_frozen_by_its_sha256_and_is_the_default():
    path = ncrna.set_path("prok-2")
    assert hashlib.sha256(path.read_bytes()).hexdigest() == ncrna.FROZEN_SETS["prok-2"]
    assert ncrna.DEFAULT_SET == "prok-2"
    assert ncrna.set_version("prok-2") == ncrna.set_version(ncrna.DEFAULT_SET) == "prok2"


def test_prok1_still_loads_and_keeps_its_version():
    path = ncrna.set_path("prok-1")
    assert hashlib.sha256(path.read_bytes()).hexdigest() == ncrna.FROZEN_SETS["prok-1"] \
        == "092660962def1ec76466de3816558451f478637ffb43d141ef684273f68d6179"  # card E-07's file, unchanged
    assert ncrna.set_version("prok-1") == "prok1"
    assert len(ncrna.read_set("prok-1")) == 588


def test_a_changed_prok2_is_refused(monkeypatch):
    monkeypatch.setitem(ncrna.FROZEN_SETS, "prok-2", "0" * 64)
    with pytest.raises(ValueError, match="frozen model set"):
        ncrna.read_set("prok-2")


def test_the_default_producer_name_carries_prok2_and_prok1_keeps_its_old_name(monkeypatch):
    """The name is read from the binary, the library and the set: the binary's version is fixed here (the test image
    has it; `test_micro_s0_ncrna` checks it against `cmsearch -h`), the library is the micro fixture."""
    monkeypatch.setattr(ncrna, "infernal_version", lambda: "1.1.5")
    lib = "sha" + hashlib.sha256(LIB.read_bytes()).hexdigest()[:12]
    assert ncrna.producer(LIB) == f"infernal-1.1.5-rfam-{lib}-prok2-igr200nt-clan1"
    assert ncrna.producer(LIB, "prok-2") == ncrna.producer(LIB)
    # a stage scanned with prok-1 is still named, so still read, as card E-08 wrote it
    assert ncrna.producer(LIB, "prok-1") == f"infernal-1.1.5-rfam-{lib}-prok1-igr200nt-clan1"


def test_provenance_accepts_the_prok2_producer():
    lock = provenance.locked_versions()
    for s in ("prok2", "prok1"):
        p = f"infernal-{lock['infernal']}-rfam-{lock['rfam']}-{s}-igr200nt-clan1"
        assert provenance.NCRNA_RE.match(p)["set"] == s
        assert provenance._producer_problems("x", "p", p, lock) == []
    assert provenance._producer_problems("x", "p", "infernal-0.0-rfam-15.1-prok2-igr200nt-clan1", lock)


def test_the_cli_and_build_defaults_are_the_default_set(monkeypatch, tmp_path):
    assert bs.Build.__dataclass_fields__["ncrna_set"].default == ncrna.DEFAULT_SET
    seen = {}

    def fake_run(root, stage, cm, cpus, **kw):
        seen["stage"] = kw["model_set"]
        return {}

    monkeypatch.setattr(ncrna, "run", fake_run)
    assert cli.main(["stage", "ncrna", "--root", str(tmp_path), "--stage", "S", "--cm", str(LIB)]) == 0
    monkeypatch.setattr(profile, "load", lambda path: None)
    monkeypatch.setattr(bs, "build_stage", lambda prof, b: seen.__setitem__("build", b.ncrna_set) or 0)
    assert cli.main(["build-stage", "--root", str(tmp_path), "--profile", "p.yaml", "--stage", "S", "--release",
                     "r", "--genomes", "g.json", "--objects", "o.json"]) == 0
    assert seen == {"stage": ncrna.DEFAULT_SET, "build": ncrna.DEFAULT_SET}


# --- the shard planner's rate, scaled by the set ----------------------------------------------------------------


def test_the_planning_rate_is_prok1s_rate_times_the_sets_relative_cost():
    """rate(set) = NCRNA_CPU_H_PER_10GBP x cost(set) / cost(prok-1); DEC-30: prok-2 costs 299.7 against prok-1's 254.5
    CPU-hours per 10 Gbp (derived), so its rate is 304.0 x 299.7 / 254.5. prok-1's rate is unchanged."""
    assert bs.NCRNA_CPU_H_PER_10GBP == 304.0
    assert bs.ncrna_rate("prok-1") == 304.0
    assert bs.ncrna_rate("prok-2") == pytest.approx(304.0 * 299.7 / 254.5)
    assert bs.ncrna_rate(ncrna.DEFAULT_SET) == bs.ncrna_rate("prok-2")
    assert set(bs.NCRNA_SET_COST) == set(ncrna.FROZEN_SETS)  # every frozen set has its cost


def test_a_set_file_with_a_frozen_sets_bytes_has_that_sets_rate_and_another_file_the_dearest(tmp_path):
    copy = tmp_path / "copy.tsv"
    copy.write_bytes(ncrna.set_path("prok-2").read_bytes())
    assert bs.ncrna_rate(copy) == bs.ncrna_rate("prok-2")
    other = tmp_path / "mine.tsv"
    other.write_text("acc\tname\twhole_genome\nRF00005\ttRNA\t0\n")
    assert bs.ncrna_rate(other) == max(bs.ncrna_rate(s) for s in ncrna.FROZEN_SETS)


def test_shards_are_cut_at_the_sets_rate_and_stay_inside_80_percent_of_the_budget():
    prof = profile.load(PROFILE)
    assert prof.cpu_min_per_process == 300
    assert bs.ncrna_shard_bases(prof, "prok-1") == int(0.8 * 300 * 1e10 / (304.0 * 60)) == 131_578_947
    two = bs.ncrna_shard_bases(prof, "prok-2")
    assert two == int(0.8 * 300 * 1e10 / (304.0 * 299.7 / 254.5 * 60))
    assert bs.ncrna_shard_bases(prof) == two  # the default set
    # a prok-1-sized shard scanned with prok-2 would pass 80% of the budget; prok-2's own shard does not
    assert bs.ncrna_cpu_min(131_578_947, "prok-2") > 0.8 * 300
    assert bs.ncrna_cpu_min(two, "prok-2") <= 0.8 * 300 + 1e-6
    assert bs.ncrna_cpu_min(two, "prok-2") == pytest.approx(240.0, abs=1e-3)


def test_the_build_plans_ncrna_shards_at_its_sets_rate(tmp_path):
    prof = profile.load(PROFILE)
    b = bs.Build(root=tmp_path / "root", stage="S", release="r", genomes=tmp_path / "g.json",
                 objects=tmp_path / "o.json", ncrna_cm=LIB)
    step = bs.Step("S", "ncrna-x", ("stage", "ncrna"))
    for s in ("prok-1", "prok-2"):
        plan = bs.ncrna_scan_steps(prof, dataclasses.replace(b, ncrna_set=s), step, None, 2)[0]
        assert plan.argv[plan.argv.index("--scan-bases") + 1] == str(bs.ncrna_shard_bases(prof, s))
        assert plan.argv[plan.argv.index("--ncrna-set") + 1] == s
