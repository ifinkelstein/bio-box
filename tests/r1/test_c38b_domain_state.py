"""Ring R1: card C-38b — domain state (`measured`, `inferred`, `none`, `not scanned`) unit tests.

Covers:
  - `catalog.annotated_stages_by_library(manifest)` derives stage coverage per library from the manifest's
    `l1_domains` file list (in release stage order), never from row counts, so an annotated stage with zero
    domain rows still covers its proteins as `none`.
  - `catalog.DomainStateLookup` resolves all four states (`measured`, `inferred`, `none`, `not scanned`,
    plus mixed measured+inferred -> `measured`, shared annotated+called -> `none`, unknown library ->
    `not scanned`) across multi-row-group Parquet files while reading O(1) row groups per lookup (AGENTS §4.4).
  - `nbr.summarize_rows` separates `not scanned` neighbours (`key="?"`, `unannotated=True`,
    `model_name="not scanned"`) at the `pfam` level when unscanned proteins exist, and emits no extra row when
    every neighbour protein was scanned.
  - `export.HITS_TSV_COLUMNS`, `export.LOCUS_TSV_COLUMNS`, `export.README_TEXT` and `GUIDE.md` document
    `domain_state` and `-` (`none`) vs `?` (`not scanned`).
"""

from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog
from biobox.registry.schemas import TABLES
from biobox.search import export, nbr

pytestmark = [pytest.mark.covers("C-38b", "LY-06", "LY-10", "W-05", "W-06")]

REPO = Path(__file__).resolve().parents[2]


def _uid(tag: int) -> bytes:
    return tag.to_bytes(16, "big")


@pytest.mark.covers("C-38b", "LY-10")
def test_c38b_annotated_stages_by_library_uses_manifest_files_not_rows():
    man = {
        "release_id": "r26.11.2",
        "stages": ["S0", "S1", "S2"],
        "files": [
            {"table": "protein_seq", "path": "l0seq/protein_seq/gene_call_version=v1/stage=S0/protein_prefix=all/part-0.parquet",
             "partitions": {"gene_call_version": "v1", "stage": "S0", "protein_prefix": "all"}, "rows": 100},
            {"table": "protein_seq", "path": "l0seq/protein_seq/gene_call_version=v1/stage=S1/protein_prefix=all/part-0.parquet",
             "partitions": {"gene_call_version": "v1", "stage": "S1", "protein_prefix": "all"}, "rows": 100},
            {"table": "protein_seq", "path": "l0seq/protein_seq/gene_call_version=v1/stage=S2/protein_prefix=all/part-0.parquet",
             "partitions": {"gene_call_version": "v1", "stage": "S2", "protein_prefix": "all"}, "rows": 50},
            {"table": "l1_domains", "path": "l1/domains/library=pfam-38.2/stage=S2/part-0.parquet",
             "partitions": {"library": "pfam-38.2", "stage": "S2"}, "rows": 0},  # 0 rows, still annotated!
            {"table": "l1_domains", "path": "l1/domains/library=pfam-38.2/stage=S0/part-0.parquet",
             "partitions": {"library": "pfam-38.2", "stage": "S0"}, "rows": 42},
            {"table": "l1_domains", "path": "l1/domains/library=tigr-15.0/stage=S0/part-0.parquet",
             "partitions": {"library": "tigr-15.0", "stage": "S0"}, "rows": 7},
        ],
    }
    assert catalog.annotated_stages_by_library(man) == {
        "pfam-38.2": ["S0", "S2"],
        "tigr-15.0": ["S0"],
    }
    assert catalog.annotated_stages_by_library({"release_id": "r0", "stages": ["S1"], "files": []}) == {}


def _dom_row(puid: bytes, acc: str, env_start: int, inferred_from: bytes | None = None) -> dict:
    return {
        "protein_uid": puid,
        "model_acc": acc,
        "model_name": acc.lower(),
        "model_len": 100,
        "seq_bitscore": 50.0,
        "seq_evalue_zref": 1e-12,
        "dom_bitscore": 45.0,
        "dom_ievalue_zref": 1e-12,
        "env_start": env_start,
        "env_end": env_start + 90,
        "ali_start": env_start,
        "ali_end": env_start + 90,
        "hmm_from": 0,
        "hmm_to": 90,
        "cutoff": "gathering",
        "z": 1e8,
        "dom_z": 1e8,
        "inferred_from": inferred_from,
    }


@pytest.mark.covers("C-38b", "LY-06")
def test_c38b_domain_state_lookup_all_four_states_and_prunes_row_groups(tmp_path):
    prot_schema = TABLES["protein_seq"].schema
    dom_schema = TABLES["l1_domains"].schema.append(pa.field("inferred_from", pa.binary(16), nullable=True))

    # S0 (annotated for libA): proteins 100..179 (80 rows -> 10 row groups at row_group_size=8)
    s0_uids = [_uid(i) for i in range(100, 180)]
    # S1 (called, not annotated): proteins 160..179 (shared with S0!) + 200..259 (60 S1-only proteins -> 10 row groups)
    s1_uids = [_uid(i) for i in range(160, 180)] + [_uid(i) for i in range(200, 260)]
    # S2 (annotated for libA with ZERO domain rows): proteins 300..315
    s2_uids = [_uid(i) for i in range(300, 316)]

    p_s0 = tmp_path / "p_s0.parquet"
    p_s1 = tmp_path / "p_s1.parquet"
    p_s2 = tmp_path / "p_s2.parquet"
    for path, uids in ((p_s0, s0_uids), (p_s1, s1_uids), (p_s2, s2_uids)):
        rows = [{"protein_uid": u, "protein_digest": b"\xaa" * 32, "length": 120, "seq": "M" * 120}
                for u in sorted(uids)]
        pq.write_table(pa.Table.from_pylist(rows, schema=prot_schema), path, row_group_size=8, write_statistics=True)

    # Domain rows in S0 for libA (64 rows -> 8 row groups):
    # - 100..139: measured (inferred_from=None), including 160 (shared with S1)
    # - 161 (shared with S1): BOTH measured and inferred -> state must be 'measured'
    # - 140..155 and 162..165 (shared with S1): inferred only (inferred_from=_uid(100)) -> 'inferred'
    # - 156..159 and 166..179 (shared with S1): no domain row, in S0 -> 'none'
    # - 300..315: in S2 (0-row domain file for libA) -> 'none'
    # - 200..259: in S1 only -> 'not scanned'
    dom_rows = []
    for i in [*range(100, 140), 160]:
        dom_rows.append(_dom_row(_uid(i), "PF00001", 0, None))
    dom_rows.append(_dom_row(_uid(161), "PF00001", 0, None))
    dom_rows.append(_dom_row(_uid(161), "PF00002", 100, _uid(100)))
    for i in [*range(140, 156), *range(162, 166)]:
        dom_rows.append(_dom_row(_uid(i), "PF00003", 0, _uid(100)))
    dom_rows.sort(key=lambda r: (r["protein_uid"], r["env_start"], r["model_acc"]))

    d_s0 = tmp_path / "d_s0.parquet"
    d_s2 = tmp_path / "d_s2.parquet"
    pq.write_table(pa.Table.from_pylist(dom_rows, schema=dom_schema), d_s0, row_group_size=8, write_statistics=True)
    pq.write_table(pa.Table.from_pylist([], schema=dom_schema), d_s2, row_group_size=8, write_statistics=True)

    assert pq.read_metadata(p_s0).num_row_groups == 10
    assert pq.read_metadata(d_s0).num_row_groups >= 7

    # Independent brute-force oracle dictionary
    want: dict[bytes, str] = {}
    for u in set(s0_uids) | set(s1_uids) | set(s2_uids) | {_uid(999)}:
        mine = [r for r in dom_rows if r["protein_uid"] == u]
        if any(r["inferred_from"] is None for r in mine):
            want[u] = "measured"
        elif mine:
            want[u] = "inferred"
        elif u in set(s0_uids) | set(s2_uids):
            want[u] = "none"
        else:
            want[u] = "not scanned"

    lookup = catalog.DomainStateLookup(
        domain_files=[(d_s0, "libA"), (d_s2, "libA")],
        annotated_protein_files={"libA": [p_s0, p_s2]},
    )
    assert lookup.row_groups_total >= 19

    # Single-key lookups touch O(1) row groups across every state (0 when footer min/max excludes every row group)
    for probe in (_uid(105), _uid(161), _uid(163), _uid(170), _uid(305), _uid(225), _uid(999)):
        mark = lookup.stamp
        got = lookup.states([(probe, "libA")])
        touched, _ = lookup.touched_since(mark)
        assert got[(probe, "libA")] == want[probe], probe
        lo = 1 if want[probe] != "not scanned" else 0
        assert lo <= len(touched) <= 3, (probe, len(touched), lookup.row_groups_total)

    # Batch check of all proteins + an unknown library
    pairs = [(u, "libA") for u in sorted(want)] + [(_uid(100), "unknown-lib")]
    got_all = lookup.states(pairs)
    for u, st in want.items():
        assert got_all[(u, "libA")] == st, (u, got_all[(u, "libA")], st)
    assert got_all[(_uid(100), "unknown-lib")] == "not scanned"


@pytest.mark.covers("C-38b", "W-06")
def test_c38b_nbr_summarize_rows_separates_not_scanned_only_when_present():
    anchor = _uid(1)
    p_fam = _uid(2)       # annotated with PF00001
    p_none = _uid(3)      # scanned in annotated stage, no domain
    p_unscanned = _uid(4) # in called stage only -> not scanned
    cu1, cu2 = _uid(10), _uid(11)
    fu_a1, fu_a2 = _uid(20), _uid(21)

    nbr_rows = [
        {"anchor_protein_uid": anchor, "contig_uid": cu1, "anchor_ordinal": 1, "ordinal": 0, "offset": -1,
         "offset_rel": -1, "protein_uid": p_fam, "start": 0, "end": 100, "strand": 1, "delta_bp": 10,
         "junction_topology": 1, "same_strand": True, "in_operon": None, "flags": 0, "policy_hash": "h",
         "anchor_feature_uid": fu_a1, "feature_uid": _uid(30), "tape_version": "tv1"},
        {"anchor_protein_uid": anchor, "contig_uid": cu1, "anchor_ordinal": 1, "ordinal": 1, "offset": 0,
         "offset_rel": 0, "protein_uid": anchor, "start": 110, "end": 210, "strand": 1, "delta_bp": None,
         "junction_topology": None, "same_strand": True, "in_operon": None, "flags": 0, "policy_hash": "h",
         "anchor_feature_uid": fu_a1, "feature_uid": fu_a1, "tape_version": "tv1"},
        {"anchor_protein_uid": anchor, "contig_uid": cu1, "anchor_ordinal": 1, "ordinal": 2, "offset": 1,
         "offset_rel": 1, "protein_uid": p_none, "start": 220, "end": 320, "strand": 1, "delta_bp": 10,
         "junction_topology": 1, "same_strand": True, "in_operon": None, "flags": 0, "policy_hash": "h",
         "anchor_feature_uid": fu_a1, "feature_uid": _uid(31), "tape_version": "tv1"},
        {"anchor_protein_uid": anchor, "contig_uid": cu2, "anchor_ordinal": 0, "ordinal": 0, "offset": 0,
         "offset_rel": 0, "protein_uid": anchor, "start": 0, "end": 100, "strand": 1, "delta_bp": None,
         "junction_topology": None, "same_strand": True, "in_operon": None, "flags": 0, "policy_hash": "h",
         "anchor_feature_uid": fu_a2, "feature_uid": fu_a2, "tape_version": "tv1"},
        {"anchor_protein_uid": anchor, "contig_uid": cu2, "anchor_ordinal": 0, "ordinal": 1, "offset": 1,
         "offset_rel": 1, "protein_uid": p_unscanned, "start": 110, "end": 210, "strand": -1, "delta_bp": 10,
         "junction_topology": 3, "same_strand": False, "in_operon": None, "flags": 0, "policy_hash": "h",
         "anchor_feature_uid": fu_a2, "feature_uid": _uid(32), "tape_version": "tv1"},
    ]
    doms = [{"protein_uid": p_fam, "library": "pfam-38.2", "model_acc": "PF00001", "model_name": "FamA"}]

    # 1. All-annotated (scanned_proteins=None or contains all proteins): NO '?' row emitted
    all_ann = nbr.summarize_rows(nbr_rows, by=("pfam",), domains=doms)
    assert [(r["key"], r["unannotated"], r["n_neighbors"]) for r in all_ann["rows"]] == [
        ("pfam-38.2:PF00001", False, 1),
        ("", True, 2),
    ]
    all_ann_explicit = nbr.summarize_rows(nbr_rows, by=("pfam",), domains=doms,
                                          scanned_proteins={anchor, p_fam, p_none, p_unscanned})
    assert all_ann_explicit == all_ann

    # 2. Mixed release (only anchor, p_fam, p_none scanned; p_unscanned not scanned): separates '' and '?'
    mixed = nbr.summarize_rows(nbr_rows, by=("pfam",), domains=doms, scanned_proteins={anchor, p_fam, p_none})
    assert [(r["key"], r["unannotated"], r["model_name"], r["n_loci"], r["n_neighbors"]) for r in mixed["rows"]] == [
        ("pfam-38.2:PF00001", False, "FamA", 1, 1),
        ("", True, None, 1, 1),
        ("?", True, "not scanned", 1, 1),
    ]
    tsv = "\n".join(export.nbr_summary_tsv({**mixed, "header": {"release": "r1", "tape_version": "tv1",
                                                                 "policy": {}, "policy_hash": "h",
                                                                 "n_anchor_proteins_given": 1,
                                                                 "anchors_dropped_by_edge": 0, **mixed["header"]}}))
    assert "\npfam\t?\ttrue\t\t\tnot scanned\t1\t0.5\t1\t" in "\n" + tsv


@pytest.mark.covers("C-38b", "W-05")
def test_c38b_readme_and_guide_explain_domain_state_and_question_mark():
    assert export.HITS_TSV_COLUMNS[export.HITS_TSV_COLUMNS.index("pfam_architecture") + 1] == "domain_state"
    assert export.LOCUS_TSV_COLUMNS[export.LOCUS_TSV_COLUMNS.index("pfam_architecture") + 1] == "domain_state"
    assert "domain_state" in export.README_TEXT
    assert "not scanned" in export.README_TEXT and "?" in export.README_TEXT
    guide = (REPO / "GUIDE.md").read_text()
    assert "domain_state" in guide and "not scanned" in guide and "`?`" in guide
