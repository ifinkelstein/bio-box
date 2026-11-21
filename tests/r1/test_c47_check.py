"""Ring R1: card C-47's post-build check, its reading of the gene-call ledger's reasons.

The gene caller writes a skip's reason as `owned_by:` + `ContigOwners.owner` (src/biobox/l0seq/store.py): the owning
stage for a contig an earlier stage stores (DEC-26), or `<stage>/<genome_uid hex>` for a contig that another genome
of the caller's own stage owns (DEC-27). On stage C5 the first check counted 3 rows `owned_by:C5/<genome_uid>` as
earlier-stage skips and failed (2,553 against 2,550, 2026-11-20).
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-47")  # ring R1 by directory (tests/conftest.py)

CHECK = Path(__file__).resolve().parents[2] / "scripts" / "stress" / "build_stage_pod" / "c47_check.py"


def _load():
    spec = importlib.util.spec_from_file_location("c47_check_t", CHECK)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c47_check_t"] = mod
    spec.loader.exec_module(mod)
    return mod


G = "1cf299c3a49afa98e9d03d8e7af92deb"


@pytest.mark.parametrize("reason,want", [
    ("owned_by:PL1", "earlier"),
    ("owned_by:P2", "earlier"),
    ("owned_by:C5", "within"),
    (f"owned_by:C5/{G}", "within"),
    ("duplicate_contig_in_genome", "within"),
    ("ok", None),
    (f"owned_by:PL1/{G}", "earlier"),
    ("owned_by:C50", "earlier"),
])
def test_skip_kind(reason, want):
    assert _load().skip_kind(reason, "C5") == want


def test_classify_ledger_counts_within_stage_genome_owners_as_within():
    m = _load()
    keys = [b"\x01" * 16, b"\x02" * 16, b"\x03" * 16, b"\x04" * 16, b"\x05" * 16]
    reasons = ["ok", "owned_by:PL1", f"owned_by:C5/{G}", "owned_by:C5", "owned_by:S1"]
    earlier, within = m.classify_ledger(keys, reasons, "C5")
    assert earlier == {(b"\x02" * 16).hex(): "owned_by:PL1", (b"\x05" * 16).hex(): "owned_by:S1"}
    assert within == 2


def test_every_stored_contig_has_genes_or_a_ledger_row():
    m = _load()
    stored = {"aa", "bb", "cc", "dd"}
    with_genes = {"aa", "bb"}
    ledger = {"aa": "ok", "bb": "ok", "cc": "ok", "ee": "owned_by:PL1"}
    uncovered, reasons = m.uncovered_contigs(stored, with_genes, ledger)
    assert uncovered == ["dd"]  # stored, no gene row, no ledger row
    assert reasons == {"ok": 1}  # cc: called, no gene, listed
    assert m.uncovered_contigs({"aa"}, {"aa"}, {}) == ([], {})
