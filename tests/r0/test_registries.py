"""Ring R0: static checks on the schema and operation registries (TEST_ARCHITECTURE §8)."""

from __future__ import annotations

import re
from pathlib import Path

import pyarrow as pa
import pytest
import yaml

from biobox.registry.schemas import SHA256, TABLES, UID

pytestmark = pytest.mark.r0

ROOT = Path(__file__).resolve().parents[2]
REG = ROOT / "tests" / "registry"
PLAN_DOCS = [
    ROOT / "PLAN.md",
    ROOT / "STAGES.md",
    ROOT / "planning" / "data_layer_and_tests.md",
    ROOT / "planning" / "search" / "SEARCH_PLAN.md",
    ROOT / "planning" / "tests" / "TEST_ARCHITECTURE.md",
]
TEST_ID = re.compile(r"\b(?:[CLRADGKS]-\d{2}[a-z]?|LY-\d{2}|Q[1-8]|B[1-6]|B-NULL|B-SIM)\b")


def _load(name: str) -> dict:
    return yaml.safe_load((REG / name).read_text())


OPS = _load("operations.yaml")
PENDING = _load("pending_tables.yaml")["pending"]
REL = _load("relations.yaml")


def _documented_test_ids() -> set[str]:
    ids: set[str] = set()
    for doc in PLAN_DOCS:
        ids.update(TEST_ID.findall(doc.read_text()))
    return ids


# --- schema registry -------------------------------------------------------


@pytest.mark.parametrize("name", sorted(TABLES))
def test_keys_exist_and_are_not_nullable(name):
    spec = TABLES[name]
    cols = {f.name: f for f in spec.schema}
    for col in spec.primary_key:
        if col in spec.partition_keys:
            continue  # comes from the path; never null
        assert col in cols, f"{name}: primary key column {col} missing"
        assert not cols[col].nullable, f"{name}: primary key column {col} is nullable"
    for col in spec.sort_key:
        assert col in cols, f"{name}: sort key column {col} missing"


@pytest.mark.parametrize("name", sorted(TABLES))
def test_partition_keys_are_path_only(name):
    spec = TABLES[name]
    overlap = set(spec.partition_keys) & set(spec.schema.names)
    assert not overlap, f"{name}: partition keys duplicated as file columns: {overlap}"


@pytest.mark.parametrize("name", sorted(TABLES))
def test_uid_and_digest_widths(name):
    # LY-01 / F-02: UIDs are 16-byte BLAKE3, digests 32-byte sha256; never integers.
    for f in TABLES[name].schema:
        if f.name.endswith("_uid") or f.name in {"member_of", "pack_id", "run_id"}:
            if f.name == "source_uid":
                assert f.type == pa.uint64(), f"{name}.{f.name} must be uint64 (data_layer §2)"
            else:
                assert f.type == UID, f"{name}.{f.name} must be binary(16), got {f.type}"
        if f.name.endswith("_digest") and f.name != "runtime_image_digest" and f.name != "test_image_digest":
            assert f.type == SHA256, f"{name}.{f.name} must be binary(32), got {f.type}"


@pytest.mark.parametrize("name", sorted(TABLES))
def test_enums_refer_to_string_or_int_columns(name):
    spec = TABLES[name]
    for col, values in spec.enums.items():
        assert col in spec.schema.names, f"{name}: enum on missing column {col}"
        assert values and len(set(values)) == len(values), f"{name}.{col}: empty or duplicate enum values"


def test_small_ok_tables_match_ly04_list():
    # TEST_ARCHITECTURE LY-04 names these metadata tables; anything else small_ok must be a waived kind.
    named = {"l0_manifest", "seq_locator", "id_history"}
    waived_kinds = {"search", "tests"}
    for name, spec in TABLES.items():
        if spec.small_ok:
            assert name in named or spec.layer in waived_kinds, f"{name}: small_ok without an LY-04 waiver"


def test_specified_and_pending_are_disjoint():
    both = set(TABLES) & set(PENDING)
    assert not both, f"tables both specified and pending: {both}"


# --- operation registry ----------------------------------------------------


def test_every_operation_is_complete():
    for op, entry in OPS["operations"].items():
        for key in ("summary", "inputs", "outputs", "oracle", "mrs", "tests"):
            assert entry.get(key), f"{op}: missing {key}"


def test_operation_tables_are_known():
    known = set(TABLES) | set(PENDING) | set(OPS["external"])
    for op, entry in OPS["operations"].items():
        for t in entry["inputs"] + entry["outputs"]:
            assert t in known, f"{op}: unknown table or artifact {t}"


def test_operation_mrs_exist_and_every_mr_is_used():
    defined = set(REL["mrs"])
    used: set[str] = set()
    for op, entry in OPS["operations"].items():
        for mr in entry["mrs"]:
            assert mr in defined, f"{op}: unknown relation {mr}"
            used.add(mr)
    assert defined - used == set(), f"relations used by no operation: {sorted(defined - used)}"


def test_operation_tests_are_documented():
    documented = _documented_test_ids()
    for op, entry in OPS["operations"].items():
        for t in entry["tests"]:
            if t != "NEW":
                assert t in documented, f"{op}: test {t} is not defined in the plan docs"


def test_s1_failures_have_a_real_test():
    documented = _documented_test_ids() | set(REL["mrs"])
    for fid, f in REL["failures"].items():
        real = [t for t in f["caught_by"] if t != "NEW"]
        for t in real:
            assert t in documented, f"{fid}: {t} is not defined"
        if f["severity"] == "S1":
            assert real, f"{fid}: S1 failure with no defined test"


def test_operations_match_test_architecture():
    doc = (ROOT / "planning" / "tests" / "TEST_ARCHITECTURE.md").read_text()
    in_doc = set(re.findall(r"^\| (O-[A-Z]+) \|", doc, flags=re.M))
    assert in_doc == set(OPS["operations"]), (
        f"registry vs TEST_ARCHITECTURE §2: only in doc {in_doc - set(OPS['operations'])}, "
        f"only in registry {set(OPS['operations']) - in_doc}"
    )
