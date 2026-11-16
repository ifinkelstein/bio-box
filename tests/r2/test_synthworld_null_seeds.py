"""Ring R2: the N0 null worlds keep the root families' phyletic patterns in distribution over the 22 seeds
(synthworld_spec §4.1 N0, §5.1; DEC-48 (1)'s fallback; card F-02 session 13).

Until session 13 the two-sample Kolmogorov-Smirnov test of the root families' genomes-present counts (N0 against its
signal world, D below the alpha = 0.01 critical value) was asserted at the null module's two seeds
(tests/r2/test_synthworld_nulls.py, run at the default seed and at 20261113). DEC-48 (1): decoys reach 3 transfers
each at every seed only at rates that redraw the whole world, and no placement of their partners does it at a low
rate (planning/measurements/2026-11-14_synthworld_redesign.md §16), so the KS test is judged over the 22 seeds instead:
the default seed, 20261113 and the 20 N2 seeds. At alpha = 0.01 a test of a world drawn from the same process fails
with probability at most 0.01, so over 22 independent seeds the number of failing seeds X is at most
Binomial(22, 0.01) in distribution; the bound is the smallest k with P(X > k) <= 0.01, computed below (k = 2:
P(X >= 3) = 0.0013, P(X >= 2) = 0.0202). The bound applies to N0-full and to N0-strict separately.

Each seed's signal, N0-full and N0-strict backbones (S1-S2 only, about 40 s each) are built in a child process that
prints the three worlds' root presence counts and removes them; the KS statistic is computed here with a loop of its
own.
"""

from __future__ import annotations

import importlib.util
import json
import math
import os
import subprocess
import sys

import pytest

from biobox.synthworld import backbone, nulls

pytestmark = [pytest.mark.r2, pytest.mark.covers("F-02", "B-NULL", "MR-18", "DEC-48")]

SEED2 = 20261113                       # session 9's second seed
ALPHA = 0.01                           # the per-seed level of the KS bound (test_synthworld_nulls.py)
KS_C = 1.628                           # the asymptotic two-sample critical constant at alpha = 0.01
JOBS = max(1, min(3, len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else 1))

CHILD = r'''
import csv, json, shutil, sys
from pathlib import Path
from biobox.synthworld import backbone, nulls
seed, base = int(sys.argv[1]), Path(sys.argv[2])
p = backbone.load_preset("micro")
out = {}
for tag, params in (("signal", p), ("n0full", nulls.null_preset(p, "N0-full")),
                    ("n0strict", nulls.null_preset(p, "N0-strict"))):
    w = base / tag
    backbone.build(params, seed, w)
    with open(w / "backbone" / "families.tsv", newline="", encoding="utf-8") as fh:
        out[tag] = [int(r["n_genomes_present"]) for r in csv.DictReader(fh, delimiter="\t") if r["kind"] == "root"]
    shutil.rmtree(w)
print(json.dumps(out))
'''


def _seeds() -> list[int]:
    default = int(backbone.load_preset("micro")["default_seed"])
    return [default, SEED2] + nulls.n2_seeds(default, 20)


def allowed_failures(n: int, alpha: float) -> int:
    """The smallest k with P(Binomial(n, alpha) > k) <= alpha."""
    tail, k = 1.0, 0
    while True:
        tail -= math.comb(n, k) * alpha ** k * (1 - alpha) ** (n - k)
        if tail <= alpha + 1e-12:
            return k
        k += 1


def ks(a: list[int], b: list[int]) -> float:
    """Largest gap between the two empirical distribution functions, over every observed value."""
    best = 0.0
    for v in set(a) | set(b):
        best = max(best, abs(sum(x <= v for x in a) / len(a) - sum(x <= v for x in b) / len(b)))
    return best


@pytest.fixture(scope="module")
def presence(tmp_path_factory):
    if importlib.util.find_spec("zombi2") is None:
        pytest.skip("zombi2 0.54.0 is not importable here")
    base = tmp_path_factory.mktemp("null_seeds")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    pending = list(_seeds())
    running, done = [], {}
    while pending or running:
        while pending and len(running) < JOBS:
            s = pending.pop(0)
            running.append((s, subprocess.Popen([sys.executable, "-c", CHILD, str(s), str(base / str(s))], env=env,
                                                cwd=base, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)))
        s, p = running.pop(0)
        out, err = p.communicate()
        if p.returncode:
            raise AssertionError(f"seed {s} failed:\n{err[-4000:]}")
        done[s] = json.loads(out.strip().splitlines()[-1])
    return done


def test_allowed_failures_is_the_binomial_bound():
    assert allowed_failures(22, 0.01) == 2
    assert allowed_failures(1, 0.01) == 0 and allowed_failures(100, 0.05) == 9


@pytest.mark.covers("F-02", "B-NULL", "MR-18", "DEC-48")
def test_n0_root_presence_ks_fails_at_no_more_seeds_than_its_level_allows(presence, record_property):
    k = allowed_failures(len(presence), ALPHA)
    for null in ("n0full", "n0strict"):
        failing = {}
        for s, w in presence.items():
            a, b = w[null], w["signal"]
            d = ks(a, b)
            crit = KS_C * math.sqrt((len(a) + len(b)) / (len(a) * len(b)))
            print(f"\n{null} seed {s}: root presence KS D {d:.4f} (bound {crit:.4f})", end="")
            if d >= crit:
                failing[s] = round(d, 4)
        record_property(f"{null}_failing_seeds", len(failing))
        print(f"\n{null}: KS fails at {len(failing)} of {len(presence)} seeds {failing}; allowed {k}")
        assert len(failing) <= k, (null, failing, k)
