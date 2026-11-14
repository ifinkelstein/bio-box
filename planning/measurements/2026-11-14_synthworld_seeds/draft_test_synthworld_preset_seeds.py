# DRAFT (card F-02 session 10), kept beside planning/measurements/2026-11-14_synthworld_seeds.md and outside tests/:
# it fails on the current preset by design (3 of 3 failed, wilkcomp02, 2026-11-14, `scripts/pod_pytest.sh F-02-10
# tests/r2/test_synthworld_preset_seeds.py`, 648 s; N2 seed 6036017085589997374 alone takes about 10 minutes because
# its genomes reach 38,757 genes). The session that fixes the preset moves it to tests/r2/ and runs it first.
"""Ring R2: the `micro` preset holds the §5.1 genome size and its planted share at several seeds (synthworld_spec
§2.1 S1-S2, §5.1, §4.1 N2, §6.1; card F-02 session 10).

Session 9 found the preset full-size only at its default seed (701 genes per genome at seed 20261113, about 70% of
them planted families; planning/measurements/2026-11-13_synthworld_nulls.md). The calibration of the statistics needs
the 20 N2 seeds and a nightly rotating seed (§5.1), so the preset must give a world inside the spec's ranges at any
seed. This test builds the backbone only (S1-S2, about 30 s a seed) at the default seed, session 9's second seed and
the first four N2 seeds, each in a child process, and counts everything from the emitted tables with loops written
here; the species tree is read from ZOMBI2's Newick with a parser of its own.

Ranges (each stated where it is asserted): the chromosome mean of §5.1 (1,500 +- 30%, the band of
`test_micro_backbone_parameters`); the planted share at most PLANTED_SHARE_MAX, derived from §5.1's counts; both
plasmid types alive and at most two plasmids per genome (§2.1 S2, §5.1); in every cell of the strength grid a module
whose members are all present in at least MIN_GENOMES genomes (PLAN.md §4.3 min_detectable_count, 6-9
co-occurrences: the upper end). The crown age of each tree is printed.
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
#: §5.1 chromosome mean 1,500 genes, +- 30% (the band test_micro_backbone_parameters asserts at the default seed)
CHROM_MEAN = (1050.0, 1950.0)
#: derived: §5.1 declares 25 modules of sizes 2-6 cycling (100 member families), 3 profile-only groups of 3 and
#: 5 decoys, 114 families; at one copy per genome against the §5.1 chromosome mean of 1,500 that is 7.6% of the
#: genes; twice that (room for paralogs and unequal genome sizes) is the bound
PLANTED_SHARE_MAX = 2 * 114 / 1500
#: PLAN.md §4.3 min_detectable_count: about 6-9 vehicle-independent co-occurrences (DS-08)
MIN_GENOMES = 9
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


def _depths(nwk: Path) -> tuple[dict[str, float], float]:
    """Root-to-node depth of every labelled node of a Newick tree (root at 0) and the root's own branch (the stem),
    parsed here."""
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
    root = next(x for x, p in parent.items() if p is None)
    return depth, length[root]


def _world(out: Path, params: dict) -> dict:
    """Everything the tests assert, counted from one backbone's tables."""
    b = out / "backbone"
    genomes = _tsv(b / "genomes.tsv")
    gids = [g["genome_id"] for g in genomes]
    kind = {r["family"]: r["kind"] for r in _tsv(b / "families.tsv")}
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
    depth, _stem = _depths(out / "zombi2" / "species_complete.nwk")   # depths start at the first speciation
    present_time = max(depth[g["species_node"]] for g in genomes)
    m = params["modules"]
    cells = collections.defaultdict(list)
    for g in _tsv(b / "groups.tsv"):
        if g["kind"] != "planted":
            continue
        i = int(g["name"][3:])
        co = set(gids)
        for f in g["members"].split(";"):
            co &= present.get(f, set())
        cells[(m["together"][i % len(m["together"])],
               m["member_transfer"][(i // len(m["together"])) % len(m["member_transfer"])])].append(len(co))
    return {"n_genomes": len(gids), "chrom_mean": sum(chrom[g] for g in gids) / len(gids),
            "planted_share": planted / len(genes), "plasmid_genomes": {k: len(v) for k, v in plasmids.items()},
            "max_plasmids": max(per_genome.values(), default=0), "crown_age": present_time, "cells": cells}


@pytest.fixture(scope="module")
def measured(worlds):
    params = backbone.load_preset("micro")
    out = {s: _world(w, params) for s, w in worlds.items()}
    for s, r in out.items():
        print(f"\nseed {s}: crown age {r['crown_age']:.4f}; chromosome mean {r['chrom_mean']:.1f}; planted share "
              f"{r['planted_share']:.4f}; plasmid genomes {dict(sorted(r['plasmid_genomes'].items()))}; module "
              f"genomes by cell {dict(sorted((k, sorted(v)) for k, v in r['cells'].items()))}")
    return params, out


@pytest.mark.covers("F-02", "B-SIM")
def test_genome_size_and_planted_share_hold_at_every_seed(measured):
    params, out = measured
    bad = {s: (round(r["chrom_mean"], 1), round(r["planted_share"], 4)) for s, r in out.items()
           if not (CHROM_MEAN[0] <= r["chrom_mean"] <= CHROM_MEAN[1] and r["planted_share"] <= PLANTED_SHARE_MAX)}
    assert all(r["n_genomes"] == int(params["species"]["n_extant"]) for r in out.values())
    assert not bad, f"(chromosome mean, planted share) outside {CHROM_MEAN} and <= {PLANTED_SHARE_MAX:.3f}: {bad}"


@pytest.mark.covers("F-02", "B-SIM")
def test_both_plasmid_types_live_at_every_seed(measured):
    params, out = measured
    names = {f"pls{k}" for k in range(len(params["plasmids"]["sizes"]))}
    bad = {s: r["plasmid_genomes"] for s, r in out.items()
           if set(r["plasmid_genomes"]) != names or r["max_plasmids"] > 2}
    assert not bad, bad


@pytest.mark.covers("F-02", "B-SIM")
def test_every_strength_cell_keeps_a_detectable_module_at_every_seed(measured):
    """The PI's grid (2026-11-12): member transfer {0.05, 0.15, 0.30} x together {0.3, 0.7, 1.0}; each of the nine
    cells holds a module whose members are all present in at least MIN_GENOMES genomes."""
    params, out = measured
    m = params["modules"]
    assert (m["together"], m["member_transfer"]) == ([0.3, 0.7, 1.0], [0.05, 0.15, 0.30])
    bad = {s: {c: v for c, v in r["cells"].items() if max(v) < MIN_GENOMES} for s, r in out.items()}
    assert all(len(r["cells"]) == 9 for r in out.values())
    assert not any(bad.values()), {s: v for s, v in bad.items() if v}
