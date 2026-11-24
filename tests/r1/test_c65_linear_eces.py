"""Ring R1: Card C-65 — staging giant linear extrachromosomal elements (ECE1).

Tests:
- Exact TIR detection on linear replicons without EVER trimming a single base (DEC-65).
- Replicon class assignment (`borg`, `linear_megaplasmid`, `phage_plasmid`), host taxon extraction,
  and strict topology handling (`circular` only when explicit public metadata marks `circular`;
  Borgs and linear plasmids always staged `linear`, `kind="plasmid"`).
- Refusal to mix `redistributable: true` and `redistributable: false` sources in `ECE1` (DEC-51).
- Within-stage `contig_uid` deduplication (DEC-27) and cross-stage build-order ownership (DEC-26).
- Deterministic sorted `ece_metadata.tsv` sidecar, `inputs_SUMS`, and `trims.tsv` (0 trims).
"""
from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-65")

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"


def _mod():
    path = BSP / "c65_linear_eces.py"
    spec = importlib.util.spec_from_file_location("c65_linear_eces", path)
    assert spec is not None and spec.loader is not None, f"missing {path}"
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c65_linear_eces"] = mod
    spec.loader.exec_module(mod)
    return mod


def _revcomp(seq: bytes) -> bytes:
    return seq.translate(bytes.maketrans(b"ACGTacgt", b"TGCAtgca"))[::-1]


def test_detect_tir_finds_exact_terminal_inverted_repeat_without_trimming():
    m = _mod()
    tir = b"ACGTACGTACGTTTGCAACGTACGTACGT" * 2  # 58 bp
    body = b"AAACCCGGGTTTAGCA" * 20
    seq = tir + body + _revcomp(tir)
    assert m.detect_tir(seq) == len(tir)
    # Non-TIR sequence returns 0
    assert m.detect_tir(b"AAAAAAAAAA" + body + b"CCCCCCCCCC") == 0
    # Short (< min_len) inverted repeat returns 0
    short = b"ACGTACGTAC"  # 10 bp
    assert m.detect_tir(short + body + _revcomp(short), min_len=20) == 0


def test_classify_ece_record_enforces_linear_for_borgs_and_linear_plasmids():
    m = _mod()
    borg = m.classify_ece_record(
        "CP113311.1 MAG: Candidatus Methanoperedens sp. isolate ECE 1 plasmid unnamed",
        default_class="borg",
        topology="circular",  # even if mis-passed, Borgs are linear unless explicitly circular non-Borg
    )
    assert borg["ece_class"] == "borg"
    assert borg["circular"] is False
    assert "Methanoperedens" in borg["host_taxon"]

    lin_mp = m.classify_ece_record(
        "CP108999.1 Streptomyces sp. NBC_00056 plasmid p1, complete sequence",
        default_class="linear_megaplasmid",
        topology="linear",
    )
    assert lin_mp["ece_class"] == "linear_megaplasmid"
    assert lin_mp["circular"] is False
    assert "Streptomyces" in lin_mp["host_taxon"]

    circ_pp = m.classify_ece_record(
        "NZ_CP010398.1 Acinetobacter baumannii plasmid pAB1, complete sequence",
        default_class="phage_plasmid",
        topology="circular",
    )
    assert circ_pp["ece_class"] == "phage_plasmid"
    assert circ_pp["circular"] is True
    assert "Acinetobacter" in circ_pp["host_taxon"]


def test_stage_ece1_preserves_tirs_byte_identically_and_deduplicates(tmp_path):
    from biobox import ids

    m = _mod()
    tir = b"ACGTTGCAAACCGGTTACGTTGCAAACCGGTT"  # 32 bp
    borg_seq = tir + (b"ATGCGATCGATCGATC" * 100) + _revcomp(tir)
    mp_seq = (b"GCGCGCATATATATGC" * 80)
    pp_seq = (b"TTGGAACCGGTTAACC" * 60)

    fa_borg = tmp_path / "borg.fna"
    fa_borg.write_bytes(
        b">CP113311.1 MAG: Candidatus Methanoperedens sp. isolate ECE 1 plasmid unnamed\n"
        + borg_seq
        + b"\n"
    )
    fa_mp = tmp_path / "mp.fna"
    fa_mp.write_bytes(
        b">CP108999.1 Streptomyces sp. NBC_00056 plasmid p1, complete sequence\n"
        + mp_seq
        + b"\n"
        + b">CP108999_dup.1 Streptomyces sp. NBC_00056 plasmid p1 duplicate\n"
        + mp_seq
        + b"\n"
    )
    fa_pp = tmp_path / "pp.fna"
    fa_pp.write_bytes(
        b">NZ_CP010398.1 Acinetobacter baumannii plasmid pAB1, complete sequence\n"
        + pp_seq
        + b"\n"
    )

    # Simulate an earlier staged owner (PL2) holding pp_seq
    owner_pl2 = tmp_path / "PL2"
    owner_pl2.mkdir()
    pp_cu = ids.contig_uid(pp_seq).hex()
    (owner_pl2 / "genomes.json").write_text(json.dumps([{"contig_uids": [pp_cu]}]))

    sources = [
        {
            "source_name": "ncbi_eutils_c65_ece1",
            "source_release": "2026-11-24",
            "license": "INSDC (public)",
            "redistributable": True,
            "fasta": fa_borg,
            "output_name": "ece1_borgs.fna",
            "ece_class": "borg",
            "topology_map": {"CP113311.1": "linear"},
        },
        {
            "source_name": "ncbi_eutils_c65_ece1",
            "source_release": "2026-11-24",
            "license": "INSDC (public)",
            "redistributable": True,
            "fasta": fa_mp,
            "output_name": "ece1_actinomycetia.fna",
            "ece_class": "linear_megaplasmid",
            "topology_map": {"CP108999.1": "linear", "CP108999_dup.1": "linear"},
        },
        {
            "source_name": "pfeifer_phage_plasmids_github",
            "source_release": "2026-11-24",
            "license": "MIT",
            "redistributable": True,
            "fasta": fa_pp,
            "output_name": "ece1_phage_plasmids.fna",
            "ece_class": "phage_plasmid",
            "topology_map": {"NZ_CP010398.1": "circular"},
        },
    ]

    out_dir = tmp_path / "ECE1"
    rep = m.stage_ece1(
        sources,
        out_dir,
        root=None,
        root_stages=(),
        owner_inputs={"PL2": owner_pl2},
        retrieved="2026-11-24T00:00:00Z",
    )

    # Byte-identical TIR preservation (0 trims, exact contig_uid match)
    assert rep["trims"] == 0
    assert rep["records"] == 4
    assert rep["genomes"] == 3
    assert rep["repeated_records"] == 1
    assert rep["held_contigs"] == 1
    assert rep["new_contigs"] == 2
    assert rep["owned_by_stage"] == {"PL2": 1}

    genomes = json.loads((out_dir / "genomes.json").read_text())
    by_cu = {g["contig_uids"][0]: g for g in genomes}
    borg_cu = ids.contig_uid(borg_seq).hex()
    assert borg_cu in by_cu
    assert by_cu[borg_cu]["circular"] == []
    assert by_cu[borg_cu]["kind"] == "plasmid"
    assert by_cu[pp_cu]["circular"] == [pp_cu]
    assert by_cu[pp_cu]["kind"] == "plasmid"

    # Sidecar TSV checks: TIR length recorded for borg, owner_stage recorded for pp
    with open(out_dir / "ece_metadata.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 3
    by_acc = {r["accession"]: r for r in rows}
    assert int(by_acc["CP113311.1"]["tir_length"]) == len(tir)
    assert by_acc["CP113311.1"]["ece_class"] == "borg"
    assert by_acc["NZ_CP010398.1"]["owner_stage"] == "PL2"
    assert by_acc["NZ_CP010398.1"]["circular"] == "true"

    # inputs_SUMS exists and covers genomes.json, objects.json, and FASTAs
    sums_text = (out_dir / "inputs_SUMS").read_text()
    assert "genomes.json" in sums_text
    assert "objects.json" in sums_text
    assert "ece1_borgs.fna" in sums_text


def test_stage_ece1_refuses_non_redistributable_source(tmp_path):
    m = _mod()
    fa = tmp_path / "x.fna"
    fa.write_bytes(b">a\n" + b"ACGT" * 50 + b"\n")
    sources = [
        {
            "source_name": "restricted",
            "source_release": "1",
            "license": "restricted",
            "redistributable": False,
            "fasta": fa,
            "output_name": "x.fna",
            "ece_class": "linear_megaplasmid",
        }
    ]
    with pytest.raises(m.Refused):
        m.stage_ece1(
            sources,
            tmp_path / "ECE1",
            root=None,
            root_stages=(),
            owner_inputs={},
            retrieved="2026-11-24T00:00:00Z",
        )
