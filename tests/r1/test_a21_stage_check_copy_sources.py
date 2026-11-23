"""Ring R1: card A-21, `stage_check.py` `copy_check` accounts for shared-contig `owned` proteins (C-38h) and
`round` sources (C-38f), matching `biobox.l1.domains._measured_uids`.

All roots are constructed here (no data root, no MinCED, MMseqs2 or hmmsearch: the CRISPR producer's name comes
from a fixed MinCED version set with `monkeypatch`, as `tests/r1/test_stage_driver.py` does).
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("A-21", "A-12", "C-38f", "C-38h")

HERE = Path(__file__).resolve().parent


def _td():
    spec = importlib.util.spec_from_file_location("test_stage_driver", HERE / "test_stage_driver.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["test_stage_driver"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def no_minced(monkeypatch):
    """The producer names without the MinCED binary: a fixed version."""
    from biobox.l1 import crispr

    monkeypatch.setattr(crispr, "minced_version", lambda: "0.4.2")


def _owned_source_mut(td, *, alter_copied: bool = False):
    """Stage C5 (`"C5" < "P2"`) measured shared-contig protein `pr3` (on contig H owned by P2) via `s["owned"]`,
    while P2 also holds `pr2` (which ZQ stores on B). C5 appears before P2 in `copy["sources"]`."""
    from biobox.l1 import genecall

    gcv = genecall.GENE_CALL_VERSION

    def mut(S):
        H = S["u"]["H"]
        pr1 = S["u"]["pr1"]
        pr2 = td._u(22)
        pr3 = S["u"]["pr3"]
        pr_c5 = td._u(25)
        g_c5 = td._u(15)
        H2 = td._u(5)
        T = S["T"]

        # Replace OA with P2 as the owner of H (and also owner of H2 carrying pr2).
        for k in ("oaci", "oafp", "oaps", "oadom"):
            T.pop(k, None)
        S["records"].pop("OA", None)

        T["gcled"][1] = [
            dict(r, reason_code="owned_by:P2") if r["reason_code"] == "owned_by:OA" else r
            for r in T["gcled"][1]
        ]
        for ck in ("cled1", "cled2"):
            T[ck][1] = [
                dict(r, reason_code="owned_by:P2") if r["reason_code"] == "owned_by:OA" else r
                for r in T[ck][1]
            ]

        T["p2ci"] = ["l0seq/contig_index/stage=P2/part-0.parquet", [
            {"contig_uid": H, "contig_len": 3000, "is_default_variant": True},
            {"contig_uid": H2, "contig_len": 2500, "is_default_variant": True},
        ]]
        T["p2fp"] = [f"l1/feature_protein/gene_call_version={gcv}/stage=P2/part-0.parquet", [
            {"contig_uid": H, "protein_uid": pr3},
            {"contig_uid": H2, "protein_uid": pr2},
        ]]
        T["p2ps"] = [f"l0seq/protein_seq/gene_call_version={gcv}/stage=P2/part-0.parquet", [
            {"protein_uid": pr2, "length": 200},
            {"protein_uid": pr3, "length": 300},
        ]]
        T["p2dom"] = ["l1/domains/library=pfam-38.2/stage=P2/part-0.parquet", [
            td._dom(pr2, "PF2", 35.0),
            td._dom(pr3, "PF3", 40.0),
        ]]

        T["c5gi"] = ["l0seq/genome_input/stage=C5/part-0.parquet", [{
            "contig_uid": H, "genome_uid": g_c5, "kind": "virus", "genetic_code": 11,
            "is_circular": False, "assembly_quality": "complete",
        }]]
        T["c5ps"] = [f"l0seq/protein_seq/gene_call_version={gcv}/stage=C5/part-0.parquet", [
            {"protein_uid": pr_c5, "length": 150},
        ]]
        T["c5dom"] = ["l1/domains/library=pfam-38.2/stage=C5/part-0.parquet", [
            td._dom(pr3, "PF3", 40.0),
        ]]

        S["records"]["P2"] = {"ingest": ["p2ci"], "genecall": ["p2fp", "p2ps"], "domains-pfam-38.2": ["p2dom"]}
        S["records"]["C5"] = {"inputs": ["c5gi"], "genecall": ["c5ps"], "domains-pfam-38.2": ["c5dom"]}

        pr3_score = 41.0 if alter_copied else 40.0
        T["dom"][1] = [
            td._dom(pr1, "PF1", 30.0),
            td._dom(pr2, "PF2", 35.0),
            td._dom(pr3, "PF3", pr3_score),
        ]

        S["plan"] = {
            "stage": "ZQ",
            "proteins": 3,
            "residues": 600,
            "shards": [{"proteins": 1, "residues": 100}],
            "copy": {
                "proteins": 2,
                "residues": 500,
                "rows": 2,
                "proteins_with_hit": 2,
                "sources": [
                    {
                        "stage": "C5",
                        "protein_files": [{"path": T["c5ps"][0]}],
                        "owned": {
                            "genome_input": [{"path": T["c5gi"][0]}],
                            "stages": [{
                                "stage": "P2",
                                "feature_protein": [{"path": T["p2fp"][0]}],
                                "protein_files": [{"path": T["p2ps"][0]}],
                            }],
                        },
                        "files": [{"path": T["c5dom"][0]}],
                        "proteins": 1,
                        "residues": 300,
                        "rows": 1,
                        "proteins_with_hit": 1,
                    },
                    {
                        "stage": "P2",
                        "protein_files": [{"path": T["p2ps"][0]}],
                        "files": [{"path": T["p2dom"][0]}],
                        "proteins": 1,
                        "residues": 200,
                        "rows": 1,
                        "proteins_with_hit": 1,
                    },
                ],
            },
            "owned": {"stages": [{"stage": "P2", "feature_protein": [{"path": T["p2fp"][0]}]}]},
        }

    return mut


def _round_source_mut(td, *, alter_copied: bool = False):
    """Copy source is a round (`kind == "round"`, `round == "r-001"`, `coverage_files` instead of `protein_files`)."""
    def mut(S):
        pr1 = S["u"]["pr1"]
        pr3 = S["u"]["pr3"]
        T = S["T"]

        T["rcov"] = ["l1/round_coverage/library=pfam-38.2/round=r-001/part-0.parquet", [
            {"protein_uid": pr3},
        ]]
        T["rdom"] = ["l1/round_domains/library=pfam-38.2/round=r-001/part-0.parquet", [
            td._dom(pr3, "PF3", 40.0),
        ]]
        if alter_copied:
            T["dom"][1] = [td._dom(pr1, "PF1", 30.0), td._dom(pr3, "PF3", 41.0)]

        S["plan"]["copy"] = {
            "proteins": 1,
            "residues": 300,
            "rows": 1,
            "proteins_with_hit": 1,
            "sources": [{
                "kind": "round",
                "round": "r-001",
                "coverage_files": [{"path": T["rcov"][0]}],
                "files": [{"path": T["rdom"][0]}],
                "proteins": 1,
                "residues": 300,
                "rows": 1,
                "proteins_with_hit": 1,
            }],
        }

    return mut


def test_copy_check_accounts_for_shared_contig_owned_proteins(tmp_path, monkeypatch, no_minced):
    td = _td()
    r = td.build_root(tmp_path, _owned_source_mut(td))
    rc, s = td.run_check(r, monkeypatch)
    assert rc == 0, s["failed"]
    assert s["failed"] == []
    cc = s["pfam"]["copy_check"]
    assert cc["copied_proteins"] == 2
    by_stage = {row["stage"]: row for row in cc["sources"]}
    assert by_stage["C5"]["proteins"] == 1
    assert by_stage["C5"]["rows_in_stage"] == 1
    assert by_stage["C5"]["rows_in_source"] == 1
    assert by_stage["P2"]["proteins"] == 1
    assert by_stage["P2"]["rows_in_stage"] == 1
    assert by_stage["P2"]["rows_in_source"] == 1


def test_copy_check_fails_when_owned_copied_row_is_altered(tmp_path, monkeypatch, no_minced):
    td = _td()
    r = td.build_root(tmp_path, _owned_source_mut(td, alter_copied=True))
    rc, s = td.run_check(r, monkeypatch)
    assert rc == 1
    assert any("copied rows from C5 differ" in f for f in s["failed"]), s["failed"]


def test_copy_check_accounts_for_round_copy_source(tmp_path, monkeypatch, no_minced):
    td = _td()
    r = td.build_root(tmp_path, _round_source_mut(td))
    rc, s = td.run_check(r, monkeypatch)
    assert rc == 0, s["failed"]
    assert s["failed"] == []
    assert s["pfam"]["sources"][0]["stage"] == "r-001"
    cc = s["pfam"]["copy_check"]
    assert cc["copied_proteins"] == 1
    assert cc["sources"][0]["stage"] == "r-001"
    assert cc["sources"][0]["proteins"] == 1
    assert cc["sources"][0]["rows_in_stage"] == 1
    assert cc["sources"][0]["rows_in_source"] == 1


def test_copy_check_fails_when_round_copied_row_is_altered(tmp_path, monkeypatch, no_minced):
    td = _td()
    r = td.build_root(tmp_path, _round_source_mut(td, alter_copied=True))
    rc, s = td.run_check(r, monkeypatch)
    assert rc == 1
    assert any("copied rows from r-001 differ" in f for f in s["failed"]), s["failed"]
