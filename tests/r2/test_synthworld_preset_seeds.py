"""Ring R2: the `micro` preset holds at several seeds (synthworld_spec §2.1 S1-S2, §5.1, §4.1 N2, §6.1; DEC-36, DEC-37;
card F-02 sessions 10 and 11).

Session 10 drafted this test (planning/measurements/2026-11-14_synthworld_seeds/draft_test_synthworld_preset_seeds.py)
and it failed on the old preset: 3 of 3 failed (session 10: 648 s; session 11's rerun of the draft on the old
generator is in planning/measurements/2026-11-14_synthworld_redesign.md). Session 11 redesigned the preset as DEC-36
and DEC-37 decided and adapted the ranges to them. The calibration of the statistics needs the 20 N2 seeds and a nightly
rotating seed (§5.1), so the preset must give a world inside these ranges at any seed. This test builds the backbone
only (S1-S2, about 40 s a seed) at the default seed, session 9's second seed and the first four N2 seeds, each in a
child process, and counts everything from the emitted tables and ZOMBI2's events with loops written here; the species
tree is read from ZOMBI2's Newick with a parser of its own. The 21-seed scan stays a script
(scripts/stress/f0210_seed_scan.py).

Ranges, each with its source:

- chromosome mean within 1,500 +- 30% at every seed (DEC-37 (1); the band of `test_micro_backbone_parameters`); the
  1,000-4,000 spread is no longer asked per seed;
- planted share (planted, profile-only and decoy genes over all genes) at most PLANTED_SHARE_MAX, derived from §5.1's
  counts (session 10's derivation, for the PI to confirm);
- both plasmid types alive at the default seed, and the carriers of each type reported per seed (DEC-44 (3), which
  amends DEC-37 (3)'s survival at every seed; §5.1's per-genome range is asserted at the default seed by
  tests/r2/test_synthworld_backbone.py);
- the crown age, from the living lineages' most recent common ancestor to the present, equal to the preset's
  `species.crown_age` (DEC-37 (2), (3));
- in all but at most one cell of the strength grid per seed a module whose members are all present in at least
  MIN_GENOMES genomes (PLAN.md §4.3 min_detectable_count, 6-9 co-occurrences: the upper end);
- the grid's transfer axis realized as distinct levels (DEC-36): pooled over the seeds here, the median over a
  level's modules of realized transfer events per member family is at least LEVEL_RATIO times the level below's, and
  at every seed the highest level's median exceeds the lowest's (a level holds 8 or 9 modules, so the middle level's
  median at one seed is noisy: over the 21 seeds it lies strictly between the others at 20 of 22, measured);
- entry exclusion (DEC-36): no transfer of a module family into a lineage that held the module's first family;
- DEC-15: every adjacency-scored unit intact (the build itself refuses an unscrambled root).
"""

from __future__ import annotations

import collections
import csv
import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

from biobox.synthworld import backbone, nulls

pytestmark = [pytest.mark.r2, pytest.mark.covers("F-02", "B-SIM", "B-NULL")]

SEED2 = 20261113                       # session 9's second seed
#: DEC-37 (1): the chromosome mean within 30% of the preset's 1,500 genes, at every seed
CHROM_MEAN = (1050.0, 1950.0)
#: derived: §5.1 declares 25 modules of sizes 2-6 cycling (100 member families), 3 profile-only groups of 3 and
#: 5 decoys, 114 families; at one copy per genome against the §5.1 chromosome mean of 1,500 that is 7.6% of the
#: genes; twice that (room for paralogs and unequal genome sizes) is the bound
PLANTED_SHARE_MAX = 2 * 114 / 1500
#: PLAN.md §4.3 min_detectable_count: about 6-9 vehicle-independent co-occurrences (DS-08)
MIN_GENOMES = 9
#: DEC-36 "the three levels are distinct": each level's median realized transfers per member, pooled over the seeds,
#: at least this many times the level below's (asserted, session 11: the refitted rates give about 2x over the 21
#: seeds, planning/measurements/2026-11-14_synthworld_redesign.md; 1.5 leaves room for six seeds' sampling)
LEVEL_RATIO = 1.5
JOBS = max(1, min(3, len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else 1))

CHILD = r'''
import sys
from biobox.synthworld import backbone
p = backbone.load_preset("micro")
backbone.build(p, int(sys.argv[1]), sys.argv[2])
'''


def _seeds() -> list[int]:
    default = int(backbone.load_preset("micro")["default_seed"])
    return [default, SEED2] + nulls.n2_seeds(default, 20)[:4]


@pytest.fixture(scope="module")
def worlds(tmp_path_factory):
    if importlib.util.find_spec("zombi2") is None:
        pytest.skip("zombi2 0.54.0 is not importable here")
    base = tmp_path_factory.mktemp("seeds")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    pending = [(s, base / str(s)) for s in _seeds()]
    running, done = [], {}
    while pending or running:
        while pending and len(running) < JOBS:
            s, out = pending.pop(0)
            running.append((s, out, subprocess.Popen([sys.executable, "-c", CHILD, str(s), str(out)], env=env,
                                                     cwd=base, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                                     text=True)))
        s, out, p = running.pop(0)
        _, err = p.communicate()
        if p.returncode:
            raise AssertionError(f"seed {s} failed:\n{err[-4000:]}")
        done[s] = out
    return done


def _tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _tree(nwk: Path) -> tuple[dict[str, str | None], dict[str, float]]:
    """Parent and root-to-node depth of every labelled node of a Newick tree (root at 0), parsed here."""
    text = nwk.read_text(encoding="utf-8").strip()
    parent: dict[str, str | None] = {}
    length: dict[str, float] = {}
    opened: list[list[str]] = []
    closed = False
    for tok in re.findall(r"\(|\)|,|;|[^(),;]+", text):
        if tok == "(":
            opened.append([])
            closed = False
        elif tok == ")":
            closed = True
        elif tok in (",", ";"):
            closed = False
        else:
            label, _, ln = tok.partition(":")
            kids = opened.pop() if closed else []
            length[label] = float(ln or 0.0)
            parent[label] = None
            for k in kids:
                parent[k] = label
            if opened:
                opened[-1].append(label)
            closed = False
    depth: dict[str, float] = {}

    def d(x: str) -> float:
        if x not in depth:
            p = parent[x]
            depth[x] = 0.0 if p is None else d(p) + length[x]
        return depth[x]

    for x in parent:
        d(x)
    return parent, depth


def _crown_age(nwk: Path, tips: set[str]) -> float:
    """Present minus the depth of the deepest common ancestor of the living tips."""
    parent, depth = _tree(nwk)
    common = None
    for t in tips:
        anc, x = set(), t
        while x is not None:
            anc.add(x)
            x = parent[x]
        common = anc if common is None else common & anc
    return max(depth[t] for t in tips) - max(depth[x] for x in common)


def _exclusion_breaches(out: Path, first: dict[int, int]) -> tuple[int, int]:
    """Replay ZOMBI2's genome_events.tsv for the module families (all rows of one event checked before any is
    applied): (module transfer rows, rows into a lineage holding a live copy of the module's first family)."""
    live = collections.defaultdict(collections.Counter)
    n = bad = 0
    batch, last = [], None

    def apply():
        nonlocal n, bad
        for kind, fam, parents, kids in batch:
            if fam in first and kind.startswith("transfer"):
                n += 1
                bad += live[first[fam]][kids[1].split("_")[0]] > 0
        for kind, fam, parents, kids in batch:
            if fam in first:
                if parents:
                    live[fam][parents.split("_")[0]] -= 1
                for k in kids:
                    live[fam][k.split("_")[0]] += 1
        batch.clear()

    for r in _tsv(out / "zombi2" / "genome_events.tsv"):
        if r["event"] != last:
            apply()
            last = r["event"]
        batch.append((r["kind"], int(r["family"]), r["parents"], [c for c in r["children"].split(";") if c]))
    apply()
    return n, bad


def _world(out: Path, params: dict) -> dict:
    """Everything the tests assert, counted from one backbone's tables and ZOMBI2's files."""
    b = out / "backbone"
    genomes = _tsv(b / "genomes.tsv")
    gids = [g["genome_id"] for g in genomes]
    fams = _tsv(b / "families.tsv")
    kind = {r["family"]: r["kind"] for r in fams}
    name_id = {r["name"]: int(r["family"]) for r in fams if r["name"]}
    genes = _tsv(b / "genes.tsv")
    chrom = collections.Counter()
    planted = 0
    present = collections.defaultdict(set)
    for g in genes:
        if g["replicon"] == "0":
            chrom[g["genome_id"]] += 1
        if kind[g["family"]] in ("planted", "profile_only", "decoy"):
            planted += 1
        present[g["family"]].add(g["genome_id"])
    plasmids = collections.defaultdict(set)
    per_genome = collections.Counter()
    for g in genes:
        if g["replicon"] != "0":
            plasmids[g["group"]].add((g["genome_id"], g["unit"]))
    for name, units in plasmids.items():
        for gid, _u in units:
            per_genome[gid] += 1
    transfers = collections.Counter(e["family"] for e in _tsv(b / "events.tsv") if e["kind"].startswith("transfer"))
    m = params["modules"]
    cells = collections.defaultdict(list)
    level_tpm = collections.defaultdict(list)
    first = {}
    groups = _tsv(b / "groups.tsv")
    for g in groups:
        if g["kind"] in ("planted", "profile_only"):
            for f in g["members"].split(";"):
                first[int(f)] = name_id[f"{g['name']}.0"]
        if g["kind"] != "planted":
            continue
        i = int(g["name"][3:])
        mem = g["members"].split(";")
        co = set(gids)
        for f in mem:
            co &= present.get(f, set())
        lv = (i // len(m["together"])) % len(m["member_transfer"])
        cells[(m["together"][i % len(m["together"])], m["member_transfer"][lv])].append(len(co))
        level_tpm[lv].append(sum(transfers[f] for f in mem) / len(mem))
    scored = {g["name"] for g in groups if g["scored"] == "true"}
    units = [u for u in _tsv(b / "truth_units.tsv") if u["group"] in scored and int(u["n_members"]) >= 2]
    n_tr, bad = _exclusion_breaches(out, first)

    def median(v):
        v = sorted(v)
        return (v[(len(v) - 1) // 2] + v[len(v) // 2]) / 2

    return {"n_genomes": len(gids), "chrom_mean": sum(chrom[g] for g in gids) / len(gids),
            "planted_share": planted / len(genes), "plasmid_genomes": {k: len(v) for k, v in plasmids.items()},
            "max_plasmids": max(per_genome.values(), default=0),
            "crown_age": _crown_age(out / "zombi2" / "species_complete.nwk", {g["species_node"] for g in genomes}),
            "cells": cells, "level_median": [median(level_tpm[k]) for k in sorted(level_tpm)],
            "level_values": [list(level_tpm[k]) for k in sorted(level_tpm)],
            "module_transfers": n_tr, "exclusion_breaches": bad,
            "scored_units": len(units), "scored_intact": sum(u["intact"] == "true" for u in units)}


@pytest.fixture(scope="module")
def measured(worlds):
    params = backbone.load_preset("micro")
    out = {s: _world(w, params) for s, w in worlds.items()}
    for s, r in out.items():
        print(f"\nseed {s}: crown age {r['crown_age']:.4f}; chromosome mean {r['chrom_mean']:.1f}; planted share "
              f"{r['planted_share']:.4f}; plasmid genomes {dict(sorted(r['plasmid_genomes'].items()))}, at most "
              f"{r['max_plasmids']} per genome; realized transfers per member by level (median) "
              f"{[round(x, 2) for x in r['level_median']]}; module transfers {r['module_transfers']}, into a holder "
              f"{r['exclusion_breaches']}; module genomes by cell "
              f"{dict(sorted((k, sorted(v)) for k, v in r['cells'].items()))}")
    return params, out


@pytest.mark.covers("F-02", "B-SIM")
def test_genome_size_and_planted_share_hold_at_every_seed(measured):
    params, out = measured
    bad = {s: (round(r["chrom_mean"], 1), round(r["planted_share"], 4)) for s, r in out.items()
           if not (CHROM_MEAN[0] <= r["chrom_mean"] <= CHROM_MEAN[1] and r["planted_share"] <= PLANTED_SHARE_MAX)}
    assert all(r["n_genomes"] == int(params["species"]["n_extant"]) for r in out.values())
    assert not bad, f"(chromosome mean, planted share) outside {CHROM_MEAN} and <= {PLANTED_SHARE_MAX:.3f}: {bad}"


@pytest.mark.covers("F-02", "B-SIM")
def test_both_plasmid_types_live_at_the_default_seed_and_are_reported_per_seed(measured):
    """DEC-44 (3), amending DEC-37 (3): both plasmid types at the default seed, and reported per seed. At the rates
    that keep §5.1's per-genome range (0-2 plasmids, some genomes without one) a type is lost in the first few
    lineages at some seeds at every rate pair tried (session 11; planning/measurements/2026-11-14_synthworld_redesign.md
    §12 to §14), so no share over the seeds is asserted; the carriers of each type are printed per seed."""
    params, out = measured
    names = {f"pls{k}" for k in range(len(params["plasmids"]["sizes"]))}
    both = [s for s, r in out.items() if set(r["plasmid_genomes"]) == names]
    print(f"\nboth plasmid types at {len(both)} of {len(out)} seeds; carriers per type "
          f"{ {s: dict(sorted(r['plasmid_genomes'].items())) for s, r in out.items()} }; most plasmids in one genome "
          f"{ {s: r['max_plasmids'] for s, r in out.items()} }")
    default = int(params["default_seed"])
    assert default in out and set(out[default]["plasmid_genomes"]) == names, (default, out[default]["plasmid_genomes"])


@pytest.mark.covers("F-02", "B-SIM")
def test_every_seed_has_the_fixed_crown_age_from_the_living_mrca(measured):
    params, out = measured
    want = float(params["species"]["crown_age"])
    bad = {s: round(r["crown_age"], 6) for s, r in out.items() if abs(r["crown_age"] - want) > 1e-6 * want}
    assert not bad, (want, bad)


@pytest.mark.covers("F-02", "B-SIM")
def test_every_strength_cell_keeps_a_detectable_module_at_every_seed(measured):
    """The grid: together {0.3, 0.7, 1.0} (the PI's, 2026-11-12, kept by DEC-36) x three member transfer rates (refitted
    by session 11, DEC-36); each of the nine cells holds a module whose members are all present in at least
    MIN_GENOMES genomes, at all but at most one cell per seed (session 11 measured every cell powered at 18 of 22
    seeds and at most one cell without at the others, scan final: the lowest transfer level, 0.01, loses
    modules in some trees)."""
    params, out = measured
    m = params["modules"]
    assert m["together"] == [0.3, 0.7, 1.0] and len(m["member_transfer"]) == 3
    bad = {s: {c: v for c, v in r["cells"].items() if max(v) < MIN_GENOMES} for s, r in out.items()}
    assert all(len(r["cells"]) == 9 for r in out.values())
    print(f"\ncells without a powered module: { {s: v for s, v in bad.items() if v} }")
    assert all(len(v) <= 1 for v in bad.values()), {s: v for s, v in bad.items() if v}


@pytest.mark.covers("F-02", "B-SIM")
def test_the_transfer_axis_gives_distinct_realized_levels(measured):
    """DEC-36: the grid's transfer axis is the realized number of transfer events per member family (backbone
    events.tsv), and the three levels are distinct: pooled over the seeds, each level's median at least LEVEL_RATIO
    times the level below's; at every seed, the highest level's median above the lowest's."""
    _, out = measured

    def median(v):
        v = sorted(v)
        return (v[(len(v) - 1) // 2] + v[len(v) // 2]) / 2

    assert all(len(r["level_values"]) == 3 for r in out.values())
    pooled = [median([x for r in out.values() for x in r["level_values"][k]]) for k in range(3)]
    print(f"\nrealized transfers per member, pooled medians by level: {pooled}")
    assert all(b >= LEVEL_RATIO * a > 0 for a, b in zip(pooled, pooled[1:])), pooled
    bad = {s: [round(x, 2) for x in r["level_median"]] for s, r in out.items()
           if not r["level_median"][2] > r["level_median"][0]}
    assert not bad, bad


@pytest.mark.covers("F-02", "B-SIM", "DEC-36", "DEC-15")
def test_entry_exclusion_and_dec15_hold_at_every_seed(measured):
    _, out = measured
    assert all(r["module_transfers"] > 0 and r["exclusion_breaches"] == 0 for r in out.values()), \
        {s: (r["module_transfers"], r["exclusion_breaches"]) for s, r in out.items()}
    assert all(r["scored_units"] > 0 and r["scored_intact"] == r["scored_units"] for r in out.values()), \
        {s: (r["scored_units"], r["scored_intact"]) for s, r in out.items()}
