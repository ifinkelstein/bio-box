"""Ring R1, card E-15 (DEC-40): the producer of two-repeat CRISPR arrays beside Cas genes, without MinCED.

- The frozen list of Cas families (`src/biobox/l1/cas_families/cas-pfam-1.tsv`) is the 31 Pfam families of the CRISPR
  worked example's roles whose name starts with `cas`, read here from the example's own SQL
  (`examples/crispr_cas_neighbourhoods/sql/q01_loci_by_domains.sql`) by a parser of this file; a copy of the list
  with other bytes is refused under the frozen name.
- The producer's name: `minced-<v>-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-<library version>`; its tail names the
  filter, the family list and the domain library, and a library version that cannot be a path value is refused.
- The C-08 provenance check (`biobox.provenance`) accepts the name wherever it appears (step, `producer=` partitions,
  tape version) and checks the step's ledger: MinCED's fields as for the default producer (card E-14), and the
  fields of the filter (`gene_call_version`, `domain_library`, `cas_families`, `cas_filter`, `overlap_producer`),
  each recorded and equal to what the name encodes. The default producer is checked as before.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

import pyarrow as pa
import pytest

from biobox import io, ledger, provenance
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l3.tapes import tape_version
from biobox.registry.schemas import table

pytestmark = pytest.mark.covers("E-15", "C-08")  # ring R1 by its directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
Q01 = REPO / "examples/crispr_cas_neighbourhoods/sql/q01_loci_by_domains.sql"
DEFAULT = "minced-0.4.2-w8-nr3-r23x47-s26x50"
P = "minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2"
STEP = f"crispr-{P}"
ARGS2 = "-searchWL 8 -minNR 2 -minRL 23 -maxRL 47 -minSL 26 -maxSL 50"


def _crispr():
    from biobox.l1 import crispr  # imported here so that a missing name fails the tests, not the collection

    return crispr


def example_cas_roles() -> dict[str, str]:
    """accession -> role of every `('PFxxxxx', 'cas...')` pair of the example's q01 SQL (an independent reading)."""
    pairs = re.findall(r"\('(PF\d{5})',\s*'([A-Za-z0-9]+)'\)", Q01.read_text())
    return {pf: role for pf, role in pairs if role.startswith("cas")}


def test_the_frozen_family_list_is_the_examples_31_cas_families():
    cr = _crispr()
    want = example_cas_roles()
    assert len(want) == 31
    assert cr.cas_families() == want
    assert cr.CAS_SET == "cas-pfam-1" and set(cr.CAS_FAMILY_SETS) == {"cas-pfam-1"}
    assert cr.CAS_FAMILY_DIR / "cas-pfam-1.tsv" == REPO / "src/biobox/l1/cas_families/cas-pfam-1.tsv"


def test_a_frozen_list_with_other_bytes_is_refused(tmp_path, monkeypatch):
    cr = _crispr()
    d = tmp_path / "cas_families"
    shutil.copytree(cr.CAS_FAMILY_DIR, d)
    f = d / "cas-pfam-1.tsv"
    f.write_text(f.read_text() + "PF99999\tcas99\n")
    monkeypatch.setattr(cr, "CAS_FAMILY_DIR", d)
    with pytest.raises(ValueError, match="frozen"):
        cr.cas_families()
    with pytest.raises(FileNotFoundError):
        cr.cas_families("cas-pfam-9")


def test_the_name_tail_names_the_filter_the_family_list_and_the_library():
    cr = _crispr()
    assert cr.near_cas_tag("38.2") == "cas1kb-caspfam1-pfam-38.2"
    assert cr.FILTER_TAG == "cas1kb" and cr.CAS_BP == 1000 and cr.NEAR_CAS_MIN_REPEATS == 2
    assert cr.CAS_LIBRARY == "pfam"
    for bad in ("", "38-2", "38/2", None, 38.2):
        with pytest.raises(ValueError, match="library version"):
            cr.near_cas_tag(bad)
    assert provenance.minced_args_of(P) == ARGS2


# --- the provenance check -------------------------------------------------------------------------------------------


def tools(lock: dict, **change) -> dict:
    """The tool versions the filtered producer writes into its ledger rows, with changes (None removes a field)."""
    out = {"minced": lock["minced"], "java": "25.0.2-internal", "pyarrow": lock["pyarrow"], "minced_args": ARGS2,
           "gene_call_version": V, "domain_library": "pfam-38.2", "cas_families": "cas-pfam-1",
           "cas_filter": "cas1kb", "overlap_producer": DEFAULT}
    out.update(change)
    return {k: v for k, v in out.items() if v is not None}


def manifest(root: Path, tv: dict | None, producer: str = P) -> dict:
    step = f"crispr-{producer}"
    files = [{"path": f"l1/feature/producer={producer}/stage=S0/part-0.parquet", "table": "feature",
              "partitions": {"producer": producer, "stage": "S0"}},
             {"path": f"l1/crispr/producer={producer}/stage=S0/part-0.parquet", "table": "l1_crispr",
              "partitions": {"producer": producer, "stage": "S0"}}]
    if tv is not None:
        ledger.write(root, "S0", step, [ledger.row(bytes([k]) * 16, "contig", 1, 1, tools=tv) for k in range(1, 4)])
        files.append({"path": ledger.path(root, "S0", step).relative_to(root).as_posix(), "table": "ledger",
                      "partitions": {"stage": "S0", "step": step}})
    return {"root": str(root), "release_id": "rX", "stages": ["S0"], "files": files}


@pytest.fixture
def lock() -> dict:
    lk = provenance.locked_versions()
    assert lk["minced"] == "0.4.2" and lk["pfam"] == "38.2"  # the conda pin and refdb.yaml's Pfam
    return lk


def test_a_stage_that_ran_the_filtered_producer_passes(lock, tmp_path):
    assert provenance.check_provenance(manifest(tmp_path, tools(lock)), lock) == []
    assert provenance._producer_problems("x", "p", P, lock) == []


@pytest.mark.parametrize("change, words", [
    ({"minced_args": ARGS2.replace("-minNR 2", "-minNR 3")}, "minced_args"),  # not the nr2 its name encodes
    ({"gene_call_version": None}, "missing ['gene_call_version']"),
    ({"gene_call_version": "pyrodigal-0.0.1-gv-0.3.2-p2"}, "gene_call_version"),
    ({"domain_library": None}, "missing ['domain_library']"),
    ({"domain_library": "pfam-37.0"}, "domain_library"),
    ({"cas_families": None}, "missing ['cas_families']"),
    ({"cas_families": "cas-pfam-2"}, "cas_families"),
    ({"cas_filter": None}, "missing ['cas_filter']"),
    ({"cas_filter": "cas10kb"}, "cas_filter"),
    ({"overlap_producer": None}, "missing ['overlap_producer']"),
    ({"overlap_producer": "minced-0.4.2-w8-nr4-r23x47-s26x50"}, "overlap_producer"),
    ({"java": ""}, "java"),
    ({"minced": "0.4.1"}, "differ"),
])
def test_a_wrong_or_missing_field_of_the_filter_fails_naming_the_step(lock, tmp_path, change, words):
    probs = provenance.check_provenance(manifest(tmp_path, tools(lock, **change)), lock)
    assert probs and all(STEP in p for p in probs), probs
    assert any(words in p for p in probs), probs


def test_the_names_tail_is_checked_where_the_name_appears(lock, tmp_path):
    # a Pfam version other than refdb.yaml's: the step, the feature and the l1_crispr partitions say so
    bad = P.replace("pfam-38.2", "pfam-37.0")
    probs = provenance.check_provenance(manifest(tmp_path, tools(lock, domain_library="pfam-37.0"), bad), lock)
    assert len(probs) == 3 and all("differs from the lock (pfam 38.2)" in p for p in probs), probs
    # a family list that is not a frozen one of the repository
    probs = provenance._producer_problems("x", "p", P.replace("caspfam1", "caspfam7"), lock)
    assert probs and "family list" in probs[0], probs
    # a library the lock does not know
    probs = provenance._producer_problems("x", "p", P.replace("pfam-38.2", "tigrfam-15.0"), lock)
    assert probs and "tigrfam" in probs[0], probs
    # incomplete or malformed tails are not the producer's form
    for name in (DEFAULT + "-cas1kb", DEFAULT.replace("nr3", "nr2") + "-cas1kb-caspfam1",
                 DEFAULT.replace("nr3", "nr2") + "-cas1kb-caspfam1-pfam", P + "-x", P.replace("cas1kb", "near1kb")):
        probs = provenance._producer_problems("x", "p", name, lock)
        assert probs and "does not name" in probs[0], (name, probs)


def test_a_tape_version_with_both_producers_passes(lock, tmp_path):
    tv = tape_version([V, DEFAULT, P])
    rel = f"l3/tape_version_def/tape_version={tv}/part-0.parquet"
    io.write_table([pa.table({"producer": sorted([V, DEFAULT, P])})], table("tape_version_def"), tmp_path / rel)
    man = {"root": str(tmp_path), "release_id": "rX", "stages": ["S0"],
           "files": [{"path": rel, "table": "tape_version_def", "partitions": {"tape_version": tv}}]}
    assert provenance.check_provenance(man, lock) == []


def test_the_default_producer_is_checked_as_before(lock, tmp_path):
    """Card E-14's check of the default producer is unchanged: its step needs none of the filter's fields."""
    tv = {"minced": lock["minced"], "java": "25", "pyarrow": lock["pyarrow"],
          "minced_args": "-searchWL 8 -minNR 3 -minRL 23 -maxRL 47 -minSL 26 -maxSL 50"}
    assert provenance.check_provenance(manifest(tmp_path, tv, DEFAULT), lock) == []
    assert provenance._step_tools(f"crispr-{DEFAULT}") == ("minced", "java", "minced_args", "pyarrow")
