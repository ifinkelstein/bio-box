"""Card B-50: card B-48's fixture stage S1 gene-called by build-stage in 3 batches, and the sha256 of its root's files.

`build_in_batches` builds it with the code on the path (the test of card B-50 with this tree's code); `tree_digest` is
the summary the test compares: sha256 over the sorted lines `<path>\\t<sha256>` of every file of the root but
run-step's `runs/` and the build records `l1/genecall_plan/` (card B-50; none before it). Run as a script with the code
of another commit first on PYTHONPATH, it prints that commit's tree, so a stage called in batches before the card can be
compared with one called after it:

    PYTHONPATH=<parent tree>/src python tests/r2/b50_batched_tree.py <empty scratch directory>

It imports only the helpers of `test_b48_genecall_batches.py` (the fixture world) and `biobox` names that card B-48's
code already had. The gene-calling rates are raised (`rate_factor(LARGEST)`: S1's largest genome at 0.4 of the margin, so
the stage is over the margin in one process) and the plan is cut into 3 batches (`genecall.plan(..., batches=3)`), in
this process, where build-stage makes the plan; both are restored after.
"""

from __future__ import annotations

import functools
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import test_b48_genecall_batches as b48  # noqa: E402

from biobox import stage as bs  # noqa: E402
from biobox.l1 import genecall  # noqa: E402

LARGEST = 0.4  # S1's largest genome as a share of the margin: one process is over it, 3 batches each under it
N_BATCHES = 3
PLAN_RECORDS = "l1/genecall_plan/"  # card B-50's build records: a build's own, never in a manifest


def world_at(d: Path) -> dict:
    """The fixture world of card B-48's tests without its one-process call: the inputs of S0 and S1."""
    d.mkdir(parents=True, exist_ok=True)
    return {"d": d, "inputs": b48.write_inputs(d / "in")}


def build_in_batches(world: dict, base: Path, *, jobs: int = 1, hook=None, batch_run=None, n: int = N_BATCHES,
                     root: Path | None = None) -> tuple[Path, int | None]:
    """S0 ingested into `base/root` (or `root`, reused as it is when it has S0), then S1 built by build-stage with
    gene calling in `n` batches, `jobs` at once, stopped after the gene-call step (`b48.stop_after_genecall`, or `hook`,
    which must raise `b48.Killed` to stop). Returns the root and build-stage's code (None: stopped by the hook)."""
    root = base / "root" if root is None else root
    if not (root / "runs").exists() and not (root / "l0seq").exists():
        b48.s0_root(world, root)
    prof = b48.c08.profile_at(base / "prof")
    rates, plan = genecall.CALL_S_PER_MBP, genecall.plan
    k = b48.rate_factor(LARGEST)
    genecall.CALL_S_PER_MBP = {m: r * k for m, r in rates.items()}
    genecall.plan = functools.partial(plan, batches=n)
    try:
        b = b48.build_of(world, root, base / "scratch", scan_jobs=jobs, scan_cpus=1)
        kw = {} if batch_run is None else {"batch_run": batch_run}
        code = b48.build(prof, b, b48.stop_after_genecall() if hook is None else hook, **kw)
    finally:
        genecall.CALL_S_PER_MBP, genecall.plan = rates, plan
    return root, code


def tree(root: Path) -> dict[str, str]:
    """sha256 of every file under `root` but `runs/` and `l1/genecall_plan/`."""
    return {k: v for k, v in b48.tree(root).items() if not k.startswith(PLAN_RECORDS)}


def tree_digest(t: dict[str, str]) -> str:
    return hashlib.sha256("".join(f"{k}\t{v}\n" for k, v in sorted(t.items())).encode()).hexdigest()


if __name__ == "__main__":
    base = Path(sys.argv[1])
    w = world_at(base / "world")
    r, c = build_in_batches(w, base / "build")
    t = tree(r)
    print(json.dumps({"code": c, "digest": tree_digest(t), "files": t,
                      "batches": sorted(p.name for p in (r / "runs/stage=S1").iterdir())}, indent=1, sort_keys=True))
