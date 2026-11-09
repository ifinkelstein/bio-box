"""Ring R1: contract checker behaves on tiny in-memory tables."""

from __future__ import annotations

import hashlib

import pyarrow as pa

import pyarrow.parquet as pq
import pytest
from blake3 import blake3
from hypothesis import assume, given, strategies as st

from biobox import catalog
from biobox.ids import KIND, feature_uid, region_uid
from biobox.io import write_table
from biobox.registry.schemas import TABLES, conformance_problems, table


@pytest.mark.parametrize("name", sorted(TABLES))
def test_registry_schema_survives_parquet_round_trip(name, tmp_path):
    # Parquet cannot hold some Arrow types (e.g. timestamp[s]); a schema that changes on
    # write would fail C-01 on every file we produce.
    spec = TABLES[name]
    path = tmp_path / f"{name}.parquet"
    pq.write_table(spec.schema.empty_table(), path)
    assert pq.read_table(path).schema.equals(spec.schema, check_metadata=False)


def _igr(strands):
    spec = table("igr")
    n = len(strands)
    return spec, pa.table(
        {
            "contig_uid": pa.array([bytes(16)] * n, pa.binary(16)),
            "start": pa.array(range(0, 10 * n, 10), pa.int64()),
            "end": pa.array(range(5, 10 * n + 5, 10), pa.int64()),
            "strand": pa.array(strands, pa.int8()),
            "left_slot": pa.array([None] * n, pa.uint32()),
            "right_slot": pa.array([0] * n, pa.uint32()),
        },
        schema=spec.schema,
    )


def test_conforming_table_has_no_problems():
    spec, t = _igr([1, -1, 1])
    assert conformance_problems(t, spec) == []


def test_enum_violation_is_reported():
    spec, t = _igr([1, 0])
    assert any("strand" in p for p in conformance_problems(t, spec))


def test_duplicate_key_is_reported():
    spec, t = _igr([1])
    assert any("duplicate" in p for p in conformance_problems(pa.concat_tables([t, t]), spec))


def test_wrong_uid_width_is_a_schema_problem():
    spec, t = _igr([1])
    bad = t.set_column(0, pa.field("contig_uid", pa.binary(8), nullable=False), pa.array([bytes(8)], pa.binary(8)))
    assert conformance_problems(bad, spec) == ["igr: schema differs from registry"]


_UID = st.binary(min_size=16, max_size=16)
_COORD = st.integers(min_value=0, max_value=2**40)
_FEATURE = st.tuples(
    _UID, _COORD, st.integers(min_value=1, max_value=2**40), st.sampled_from([1, -1]), st.sampled_from(sorted(KIND.values()))
).map(lambda t: (t[0], t[1], t[1] + t[2], t[3], t[4]))


@pytest.mark.covers("LY-01")
def test_kind_values():
    assert KIND == {"cds": 1, "ncrna": 2, "trna": 3, "crispr_array": 4, "other": 15}


@pytest.mark.covers("LY-01")
def test_feature_uid_known_answer():
    contig = bytes(range(16))
    layout = contig + (100).to_bytes(8, "big") + (400).to_bytes(8, "big") + b"\xff" + b"\x02"
    assert len(layout) == 34
    expected = blake3(layout).digest(length=16)
    assert feature_uid(contig, 100, 400, -1, KIND["ncrna"]) == expected
    assert len(feature_uid(contig, 100, 400, -1, 2)) == 16


@pytest.mark.covers("LY-01")
@given(_FEATURE, _FEATURE)
def test_feature_uid_distinct_inputs_give_distinct_uids(a, b):
    assume(a != b)
    assert feature_uid(*a) != feature_uid(*b)


@pytest.mark.covers("LY-01")
@given(_FEATURE)
def test_feature_uid_strand_kind_and_region_differ(f):
    contig, start, end, strand, kind = f
    assert feature_uid(contig, start, end, strand, kind) != feature_uid(contig, start, end, -strand, kind)
    other_kind = next(k for k in sorted(KIND.values()) if k != kind)
    assert feature_uid(contig, start, end, strand, kind) != feature_uid(contig, start, end, strand, other_kind)
    assert feature_uid(contig, start, end, strand, kind) != region_uid(contig, start, end, strand)


@pytest.mark.covers("LY-01")
def test_feature_uid_allows_origin_wrapping_end():
    assert len(feature_uid(bytes(16), 4_000_000, 4_700_000, 1, 1)) == 16


@pytest.mark.covers("LY-01")
@pytest.mark.parametrize(
    "args",
    [
        (bytes(15), 0, 5, 1, 1),
        (bytes(17), 0, 5, 1, 1),
        (bytes(16), -1, 5, 1, 1),
        (bytes(16), 5, 5, 1, 1),
        (bytes(16), 5, 4, 1, 1),
        (bytes(16), 0, 5, 0, 1),
        (bytes(16), 0, 5, 2, 1),
        (bytes(16), 0, 5, 1, 0),
        (bytes(16), 0, 5, 1, 5),
        (bytes(16), 0, 5, 1, 256),
    ],
)
def test_feature_uid_rejects_invalid_input(args):
    with pytest.raises(ValueError):
        feature_uid(*args)


# --- phase-B tables (card B-03; TASKS.md §2) ------------------------------------------------

PHASE_B = ("contig_dim", "feature", "feature_protein", "genome_input", "tape_version_def", "token_postings")
RELEASE = "rS0.test"
PRODUCERS = ("pyrodigal", "trnascan")
# TASKS.md §2: "t-" + first 12 hex digits of sha256 over the sorted producer names joined by \n
TAPE_VERSION = "t-" + hashlib.sha256("\n".join(sorted(PRODUCERS)).encode()).hexdigest()[:12]
TOKEN_BUCKET = "2a"  # both tokens below have low byte 0x2a


def _u(i: int) -> bytes:
    return bytes([i]) * 16


_C1, _C2, _G = _u(1), _u(2), _u(9)
_CDS = (_C1, 0, 300, 1, KIND["cds"])
_RNA = (_C1, 400, 500, -1, KIND["ncrna"])
_F_CDS, _F_RNA = feature_uid(*_CDS), feature_uid(*_RNA)


def _phase_b_tables() -> dict[str, tuple[str, pa.Table]]:
    """One small valid table per new spec, sorted by its sort key, and the root-relative path it goes to."""
    rows = {
        "genome_input": [
            {"contig_uid": _C1, "genome_uid": _G, "is_circular": True, "kind": "prokaryote", "genetic_code": 11,
             "assembly_quality": "complete"},
            {"contig_uid": _C2, "genome_uid": _G, "is_circular": False, "kind": "plasmid", "genetic_code": 11,
             "assembly_quality": "complete"},
        ],
        "feature": [
            {"feature_uid": _F_CDS, "contig_uid": _C1, "start": 0, "end": 300, "strand": 1, "kind": KIND["cds"]},
            {"feature_uid": _F_RNA, "contig_uid": _C1, "start": 400, "end": 500, "strand": -1, "kind": KIND["ncrna"]},
        ],
        "feature_protein": [
            {"contig_uid": _C1, "feature_uid": _F_CDS, "start": 0, "end": 300, "strand": 1, "protein_uid": _u(7),
             "genetic_code": 11, "edit_type": "none"},
        ],
        "contig_dim": [
            {"contig_uid": _C1, "genome_uid": _G, "contig_len": 1000, "is_circular": True, "tier": 1,
             "assembly_quality": "complete", "tax_order": 5, "is_species_rep": True},
            {"contig_uid": _C2, "genome_uid": _G, "contig_len": 2000, "is_circular": False, "tier": None,
             "assembly_quality": None, "tax_order": None, "is_species_rep": False},
        ],
        "tape_version_def": [{"tape_version": TAPE_VERSION, "producer": p} for p in PRODUCERS],
        "token_postings": [
            # a family token on two domains (seg_idx 0, 1) of the CDS slot, then a token on the ncRNA slot
            {"token_id": 0x2A, "contig_uid": _C1, "ordinal": 0, "seg_idx": 0, "feature_uid": _F_CDS, "start": 0,
             "end": 300, "strand": 1, "n_slots": 2, "is_circular": True, "contig_len": 1000},
            {"token_id": 0x2A, "contig_uid": _C1, "ordinal": 0, "seg_idx": 1, "feature_uid": _F_CDS, "start": 0,
             "end": 300, "strand": 1, "n_slots": 2, "is_circular": True, "contig_len": 1000},
            {"token_id": 0x12A, "contig_uid": _C1, "ordinal": 1, "seg_idx": 0, "feature_uid": _F_RNA, "start": 400,
             "end": 500, "strand": -1, "n_slots": 2, "is_circular": True, "contig_len": 1000},
        ],
    }
    parts = {
        "genome_input": "stage=S0",
        "feature": "producer=pyrodigal/stage=S0",
        "feature_protein": "gene_call_version=v1/stage=S0",
        "contig_dim": f"release={RELEASE}",
        "token_postings": f"release={RELEASE}/tape_version={TAPE_VERSION}/token_bucket={TOKEN_BUCKET}",
    }
    out = {}
    for name, r in rows.items():
        spec = table(name)
        rel = spec.location if spec.location.endswith(".parquet") else f"{spec.location}/{parts[name]}/part-0.parquet"
        out[name] = rel, pa.Table.from_pylist(r, schema=spec.schema)
    return out


def _write_phase_b(root) -> dict[str, tuple[str, pa.Table]]:
    tables = _phase_b_tables()
    for name, (rel, t) in tables.items():
        write_table(t, table(name), root / rel)
    return tables


@pytest.mark.covers("C-01", "LY-01")
@pytest.mark.parametrize("name", PHASE_B)
def test_phase_b_table_written_by_the_writer_conforms(name, tmp_path):
    rel, t = _phase_b_tables()[name]
    spec = table(name)
    assert conformance_problems(t, spec) == []
    write_table(t, spec, tmp_path / rel)
    back = pq.read_table(tmp_path / rel, partitioning=None)
    assert conformance_problems(back, spec) == []
    assert back.equals(t)


@pytest.mark.covers("LY-05", "LY-10")
def test_phase_b_tables_are_catalog_views_with_partition_columns(tmp_path):
    root = tmp_path / "root"
    tables = _write_phase_b(root)
    assert sorted(tables) == list(PHASE_B)
    man = catalog.build_manifest(root, RELEASE)
    assert sorted(f["table"] for f in man["files"]) == list(PHASE_B)
    cat = catalog.generate_catalog(man, tmp_path / "catalog.duckdb", root)
    con = catalog.open_catalog(cat, root=root)
    try:
        for name, (rel, t) in tables.items():
            spec = table(name)
            cols = [c[0] for c in con.execute(f"DESCRIBE {name}").fetchall()]
            assert cols == spec.schema.names + list(spec.partition_keys), name
            got = con.sql(f"SELECT * FROM {name}").to_arrow_table()
            assert got.num_rows == t.num_rows, name
            want = dict(p.split("=", 1) for p in rel.split("/")[:-1] if "=" in p)
            for k, v in want.items():
                assert set(got.column(k).to_pylist()) == {v}, f"{name}.{k}"
            assert sorted(got.select(spec.schema.names).to_pylist(), key=repr) == sorted(t.to_pylist(), key=repr), name
    finally:
        con.close()


@pytest.mark.covers("LY-05")
@pytest.mark.parametrize("bucket,ok", [("2a", True), ("2", True), ("all", True), ("ff", False), ("3", False)])
def test_token_postings_file_under_the_wrong_bucket_is_refused(tmp_path, bucket, ok):
    # token_bucket is the hex of the token's low byte (as slot_pairs.anchor_bucket); both tokens end in 0x2a
    rel, t = _phase_b_tables()["token_postings"]
    root = tmp_path / "root"
    write_table(t, table("token_postings"), root / rel.replace(f"token_bucket={TOKEN_BUCKET}", f"token_bucket={bucket}"))
    if ok:
        assert len(catalog.build_manifest(root, RELEASE)["files"]) == 1
    else:
        with pytest.raises(catalog.LayoutError, match="LY-05"):
            catalog.build_manifest(root, RELEASE)
