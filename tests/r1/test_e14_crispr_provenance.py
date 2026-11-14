"""Ring R1, card E-14: the C-08 provenance check knows the CRISPR array producer (`biobox.l1.crispr`, card E-12).

The producer is `minced-<version>-w<searchWL>-nr<minNR>-r<minRL>x<maxRL>-s<minSL>x<maxSL>`; its step is
`crispr-<producer>` and its ledger rows record `minced`, `pyarrow`, `java` and `minced_args`. Here a release manifest
is made by hand around ledger files written with `biobox.ledger` (no MinCED needed): the lock is the repository's
(`container/environment.yml` pins `minced=0.4.2`). The step's ledger with the locked versions and the arguments its
producer names passes; a wrong MinCED version, a missing `java` or `minced_args`, or arguments that are not the ones
the producer's name encodes fail with a message that names the step. The producer name is checked wherever it
appears: the `producer=` partition of `feature` and `l1_crispr` files and a tape version's definition. The ledger of
a real run of the producer is checked in `tests/r2/test_e14_provenance_build_stage.py`.
"""

from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pytest

from biobox import io, ledger, provenance
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l3.tapes import tape_version
from biobox.registry.schemas import table

pytestmark = pytest.mark.covers("C-08", "E-14")  # ring R1 by its directory (tests/conftest.py)

P = "minced-0.4.2-w8-nr3-r23x47-s26x50"
ARGS = "-searchWL 8 -minNR 3 -minRL 23 -maxRL 47 -minSL 26 -maxSL 50"  # what crispr.minced_args() gives
STEP = f"crispr-{P}"


def tools(lock: dict, **change) -> dict:
    """The tool versions the producer writes into its ledger rows (`crispr._tools`), with changes."""
    out = {"minced": lock["minced"], "java": "25.0.2-internal", "pyarrow": lock["pyarrow"], "minced_args": ARGS}
    out.update(change)
    return {k: v for k, v in out.items() if v is not None}


def manifest(root: Path, tv: dict | None, step: str = STEP, producer: str = P) -> dict:
    """A manifest of one stage holding the producer's feature and l1_crispr partitions and its step's ledger."""
    files = [{"path": f"l1/feature/producer={producer}/stage=S0/part-0.parquet", "table": "feature",
              "partitions": {"producer": producer, "stage": "S0"}},
             {"path": f"l1/crispr/producer={producer}/stage=S0/part-0.parquet", "table": "l1_crispr",
              "partitions": {"producer": producer, "stage": "S0"}}]
    if tv is not None:
        rows = [ledger.row(bytes([k]) * 16, "contig", 1, 1, tools=tv) for k in range(1, 4)]
        ledger.write(root, "S0", step, rows)
        files.append({"path": ledger.path(root, "S0", step).relative_to(root).as_posix(), "table": "ledger",
                      "partitions": {"stage": "S0", "step": step}})
    return {"root": str(root), "release_id": "rX", "stages": ["S0"], "files": files}


@pytest.fixture
def lock() -> dict:
    lk = provenance.locked_versions()
    assert lk["minced"] == "0.4.2"  # the conda pin the image is built with
    return lk


def test_a_stage_that_ran_the_crispr_step_passes(lock, tmp_path):
    assert provenance.check_provenance(manifest(tmp_path, tools(lock)), lock) == []
    assert provenance._step_tools(STEP) == ("minced", "java", "minced_args", "pyarrow")


@pytest.mark.parametrize("change, words", [
    ({"minced": "0.4.1"}, "differ"),                 # a MinCED other than the pinned one
    ({"pyarrow": "1.0"}, "differ"),
    ({"java": None}, "missing ['java']"),            # a field the step must record
    ({"minced_args": None}, "missing ['minced_args']"),
    ({"java": ""}, "java"),                          # recorded, but empty
    ({"minced_args": "-searchWL 8 -minNR 2 -minRL 23 -maxRL 47 -minSL 26 -maxSL 50"}, "minced_args"),  # not nr3
    ({"minced_args": "-searchWL 8 -minNR 3"}, "minced_args"),           # an argument missing
    ({"minced_args": ARGS + " -spacers"}, "minced_args"),              # an argument the name does not encode
])
def test_a_wrong_or_missing_field_fails_naming_the_step(lock, tmp_path, change, words):
    probs = provenance.check_provenance(manifest(tmp_path, tools(lock, **change)), lock)
    assert probs and all(STEP in p for p in probs), probs
    assert any(words in p for p in probs), probs


def test_a_ledger_row_without_tool_versions_fails(lock, tmp_path):
    root = tmp_path
    ledger.write(root, "S0", STEP, [ledger.row(b"\x01" * 16, "contig", 1, 1, tools=None)])
    man = manifest(root, None)
    man["files"].append({"path": ledger.path(root, "S0", STEP).relative_to(root).as_posix(), "table": "ledger",
                         "partitions": {"stage": "S0", "step": STEP}})
    probs = provenance.check_provenance(man, lock)
    assert probs == [f"{man['files'][-1]['path']}: {STEP} ledger row without tool versions"]


def test_the_producer_name_is_checked_in_partitions_and_the_step(lock, tmp_path):
    # MinCED's version in the name differs from the lock: every place that names it says so
    bad = "minced-0.4.1-w8-nr3-r23x47-s26x50"
    probs = provenance.check_provenance(manifest(tmp_path, tools(lock), step=f"crispr-{bad}", producer=bad), lock)
    assert len(probs) == 3 and all("differs from the lock (minced 0.4.2)" in p for p in probs), probs
    assert any(f"crispr-{bad}" in p for p in probs) and any("l1_crispr" in p for p in probs)
    # a name that is not the producer's form
    for name in ("minced-0.4.2", "minced-0.4.2-nr3", "minced-0.4.2-w8-nr3-r23x47"):
        probs = provenance._producer_problems("x", "p", name, lock)
        assert probs and "does not name" in probs[0], (name, probs)
    assert provenance._producer_problems("x", "p", P, lock) == []
    assert provenance._producer_problems("x", "p", "minced-0.4.2-w8-nr2-r23x47-s26x50", lock) == []  # card E-13's nr2


def test_a_tape_version_with_the_crispr_producer_passes(lock, tmp_path):
    """Card E-12: `stage tapes --producers <gene caller>,<this producer>` builds a tape version of the two; C-08
    resolves it through `tape_version_def` and now accepts the producer."""
    tv = tape_version([V, P])
    rel = f"l3/tape_version_def/tape_version={tv}/part-0.parquet"
    io.write_table([pa.table({"producer": sorted([V, P])})], table("tape_version_def"), tmp_path / rel)
    man = {"root": str(tmp_path), "release_id": "rX", "stages": ["S0"],
           "files": [{"path": rel, "table": "tape_version_def", "partitions": {"tape_version": tv}}]}
    assert provenance.check_provenance(man, lock) == []
    probs = provenance.check_provenance(man, {**lock, "minced": "0.4.1"})
    assert probs and all(tv in p and "minced" in p for p in probs), probs


def test_existing_producers_are_checked_as_before(lock):
    """The gene caller's and the ncRNA producer's names keep their checks (no existing assertion changes)."""
    assert provenance._producer_problems("x", "p", V, lock) == []
    assert provenance._producer_problems("x", "p", "pyrodigal-0.0.1-gv-0.3.2-p2", lock)
    assert provenance._producer_problems("x", "p", f"infernal-{lock['infernal']}-rfam-15.1", lock) == []
    assert provenance._producer_problems("x", "p", "infernal-0.0-rfam-15.1", lock)
    assert provenance._producer_problems("x", "p", "test-ncrna", lock)  # neither form: not a gene-call version
