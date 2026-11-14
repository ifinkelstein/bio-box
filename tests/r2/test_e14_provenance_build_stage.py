"""Ring R2, card E-14: the C-08 provenance check on a release whose stage ran the CRISPR array step for real.

The fixture is card E-12's: `build-stage --crispr` on the C-04 fixture (lambda and M. genitalium;
`test_crispr_build_stage.with_arrays`), with MinCED from the image. Its release manifest lists the producer's
`feature` and `l1_crispr` files and the ledger of step `crispr-<producer>`, whose rows carry MinCED's version, the
JVM's and MinCED's arguments as the producer recorded them. C-08 passes on it; a copy of the stage whose crispr
ledger names another MinCED, lacks `java` or carries other arguments fails, with a message naming the step; the
ledgers and files of every other step are checked as before.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pytest

from biobox import catalog, ledger, provenance

sys.path.insert(0, str(Path(__file__).parent))
import test_build_stage as c04  # noqa: E402
import test_builders_out as c08  # noqa: E402
from test_crispr_build_stage import direct, name, with_arrays  # noqa: E402,F401  (the E-12 fixtures)

pytestmark = [pytest.mark.r2, pytest.mark.covers("C-08", "E-14"),
              pytest.mark.skipif(shutil.which("minced") is None, reason="minced (MinCED) not on PATH; it is in the "
                                                                       "bio-box image, not in GitHub Actions"),
              pytest.mark.skipif(c08.MMSEQS is None, reason="build-stage clusters with MMseqs2 (in the image)")]


def test_c08_passes_on_a_release_whose_stage_ran_the_crispr_step(with_arrays):
    man = catalog.build_manifest(with_arrays, c04.REL, ["S0"])
    step = name()
    tables = {(f["table"], f["partitions"].get("step") or f["partitions"].get("producer")) for f in man["files"]}
    p = step.removeprefix("crispr-")
    assert {("ledger", step), ("feature", p), ("l1_crispr", p)} <= tables
    assert provenance.check_provenance(man) == []
    rows = ledger.read(with_arrays, "S0", step)
    tv = json.loads(rows[0]["tool_versions_json"])
    assert tv["minced_args"] == provenance.minced_args_of(p) and tv["java"]


@pytest.mark.parametrize("change, words", [({"minced": "0.4.1"}, "differ"), ({"java": None}, "missing ['java']"),
                                           ({"minced_args": "-minNR 2"}, "minced_args")])
def test_c08_fails_naming_the_step_when_the_crispr_ledger_is_wrong(with_arrays, tmp_path, change, words):
    root = tmp_path / "root"
    shutil.copytree(with_arrays, root)
    step = name()
    rows = ledger.read(root, "S0", step)
    for r in rows:
        tv = {**json.loads(r["tool_versions_json"]), **change}
        r["tool_versions_json"] = json.dumps({k: v for k, v in tv.items() if v is not None}, sort_keys=True)
    ledger.path(root, "S0", step).unlink()  # the writer never rewrites a ledger (card B-19)
    ledger.write(root, "S0", step, rows)
    man = catalog.build_manifest(root, c04.REL, ["S0"])
    probs = provenance.check_provenance(man)
    assert probs and all(step in p for p in probs), probs
    assert any(words in p for p in probs), probs
