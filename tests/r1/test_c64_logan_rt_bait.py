"""Ring R1 tests for card C-64: Targeted Logan SRA contig bait-pull (LGRT1) via logan50
for ART-RT and novel prokaryotic RTs (DEC-69, DEC-71, DEC-65, DEC-50, DEC-26, DEC-27).
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import random
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-64")

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"


def _mod():
    spec = importlib.util.spec_from_file_location("c64_logan_rt_bait", BSP / "c64_logan_rt_bait.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c64_logan_rt_bait"] = mod
    spec.loader.exec_module(mod)
    return mod


def _make_seq(seed: int, length: int) -> bytes:
    rng = random.Random(seed)
    return bytes(rng.choice(b"ACGT") for _ in range(length))


def test_c64_parse_logan_header_standard_and_cl_tag():
    m = _mod()
    # Standard logan50 Prodigal header
    h1 = "ERR1124990_176_1 # 46 # 222 # -1 # ID=31_1;partial=00;start_type=GTG"
    p1 = m.parse_logan_protein_header(h1)
    assert p1["protein_id"] == "ERR1124990_176_1"
    assert p1["cluster_id"] == "ERR1124990_176_1"
    assert p1["run_acc"] == "ERR1124990"
    assert p1["contig_id"] == "ERR1124990_176"
    assert p1["gene_idx"] == 1

    # Explicit CL: tag in header
    h2 = "SRR22462154_15314_2 CL:DRR000020_1014_3 # 100 # 900 # 1"
    p2 = m.parse_logan_protein_header(h2)
    assert p2["protein_id"] == "SRR22462154_15314_2"
    assert p2["cluster_id"] == "DRR000020_1014_3"
    assert p2["run_acc"] == "SRR22462154"
    assert p2["contig_id"] == "SRR22462154_15314"
    assert p2["gene_idx"] == 2

    # Leading CL: token
    h3 = "CL:DRR000001_404_3 SRR5580900_41660_5"
    p3 = m.parse_logan_protein_header(h3)
    assert p3["cluster_id"] == "DRR000001_404_3"


def test_c64_locate_cluster_frames_and_resolve_members(tmp_path):
    m = _mod()
    idx_path = tmp_path / "cluster-to-protein.index.tsv"
    idx_path.write_text(
        "DRR000001_404_3\t0\n"
        "DRR000711_164182_1\t1000\n"
        "ERR100000_10_1\t2500\n"
        "SRR500000_99_2\t4000\n",
        encoding="utf-8",
    )
    index = m.load_cluster_index(idx_path)
    frames = m.locate_cluster_frames(
        index,
        ["DRR000020_1014_3", "DRR000001_404_3", "SRR500000_99_2", "SRR999999_1_1"],
        total_bytes=6000,
    )
    # First two map to frame [0, 999]; last two map to frame [4000, 5999]
    assert len(frames) == 2
    assert frames[0][0] == 0 and frames[0][1] == 999
    assert frames[0][2] == {"DRR000001_404_3", "DRR000020_1014_3"}
    assert frames[1][0] == 4000 and frames[1][1] == 5999
    assert frames[1][2] == {"SRR500000_99_2", "SRR999999_1_1"}

    lines = [
        "DRR000001_404_3 DRR000001_404_3\n",
        "DRR000001_404_3 DRR155047_455_3\n",
        "DRR000020_1014_3 DRR000020_1014_3\n",
        "DRR000020_1014_3 SRR22462154_15314_1\n",
        "DRR000020_1014_3 SRR22462154_15314_2\n",
        "OTHER_CLUSTER_1_1 OTHER_CLUSTER_1_1\n",
    ]
    resolved = m.resolve_cluster_members_from_lines(
        lines, {"DRR000001_404_3", "DRR000020_1014_3"}, max_members_per_cluster=2
    )
    assert resolved["DRR000001_404_3"] == ["DRR000001_404_3", "DRR155047_455_3"]
    assert len(resolved["DRR000020_1014_3"]) == 2
    assert "DRR000020_1014_3" in resolved["DRR000020_1014_3"]
    assert "OTHER_CLUSTER_1_1" not in resolved


def test_c64_extract_baited_contigs_ge_2kb_and_dec65_dtr():
    m = _mod()
    dtr = _make_seq(42, 31)  # minia3 k=31 terminal overlap
    c_circ_core = dtr + _make_seq(43, 2_500)
    c_circ_raw = c_circ_core + dtr
    c_lin_2500 = _make_seq(44, 2_500)
    c_short_1999 = _make_seq(45, 1_999)
    c_unwanted = _make_seq(46, 5_000)

    records = [
        ("ERR1124990_10 ka:f:25.0 circular=true", c_circ_raw),
        ("ERR1124990_11 ka:f:18.2 L:+:12:+", c_lin_2500),
        ("ERR1124990_12 ka:f:10.0", c_short_1999),
        ("ERR1124990_99 ka:f:30.0", c_unwanted),
    ]
    kept, stats = m.extract_baited_contigs_from_fasta(
        records,
        wanted_contigs={"ERR1124990_10", "ERR1124990_11", "ERR1124990_12"},
        min_len=2_000,
    )
    assert len(kept) == 2
    by_id = {r["contig_id"]: r for r in kept}
    assert by_id["ERR1124990_10"]["circular"] is True
    assert by_id["ERR1124990_10"]["trimmed_k"] == 31
    assert by_id["ERR1124990_10"]["seq"] == c_circ_core
    assert by_id["ERR1124990_11"]["circular"] is False
    assert by_id["ERR1124990_11"]["trimmed_k"] == 0
    assert by_id["ERR1124990_11"]["seq"] == c_lin_2500
    assert stats["short_dropped"] == 1
    assert stats["short_bases"] == 1_999
    assert stats["dtr_trimmed"] == 1


def test_c64_stage_lgrt1_owner_dedup_and_s2_check_inputs(tmp_path):
    m = _mod()
    from biobox import ids

    s_new1 = _make_seq(101, 2_400)
    s_new2 = _make_seq(102, 3_600)
    s_shared = _make_seq(103, 2_800)

    # Create a fake earlier stage owner_inputs dir owning s_shared
    earlier_dir = tmp_path / "earlier_SPX1"
    earlier_dir.mkdir()
    uid_shared = ids.contig_uid(s_shared).hex()
    (earlier_dir / "genomes.json").write_text(
        json.dumps([{"contig_uids": [uid_shared], "genetic_code": 11, "kind": "prokaryote", "circular": [], "assembly_quality": "contig"}]),
        encoding="utf-8",
    )

    baited = [
        {
            "run_acc": "ERR1124990",
            "contig_id": "ERR1124990_10",
            "seq": s_new1,
            "circular": False,
            "trimmed_k": 0,
            "raw_len": len(s_new1),
            "bait_models": ["NTERM", "clade11-A"],
            "cluster_ids": ["ERR1124990_10_1"],
            "biome": "soil metagenome",
        },
        {
            "run_acc": "SRR22462154",
            "contig_id": "SRR22462154_20",
            "seq": s_new2,
            "circular": True,
            "trimmed_k": 31,
            "raw_len": len(s_new2) + 31,
            "bait_models": ["DRT9__DRT9"],
            "cluster_ids": ["DRR000020_1014_3"],
            "biome": "aquatic metagenome",
        },
        # Duplicate sequence in another run -> deduplicated within LGRT1
        {
            "run_acc": "SRR99999999",
            "contig_id": "SRR99999999_5",
            "seq": s_new1,
            "circular": False,
            "trimmed_k": 0,
            "raw_len": len(s_new1),
            "bait_models": ["NTERM"],
            "cluster_ids": ["ERR1124990_10_1"],
            "biome": "soil metagenome",
        },
        # Contig owned by earlier stage SPX1 -> recorded in shared_records.tsv
        {
            "run_acc": "DRR000001",
            "contig_id": "DRR000001_404",
            "seq": s_shared,
            "circular": False,
            "trimmed_k": 0,
            "raw_len": len(s_shared),
            "bait_models": ["RVT-G2L4"],
            "cluster_ids": ["DRR000001_404_3"],
            "biome": "marine metagenome",
        },
    ]

    meta_dir = tmp_path / "meta_LGRT1"
    fasta_dir = tmp_path / "stor_LGRT1"
    rep = m.stage_lgrt1(
        baited_contigs=baited,
        out_meta_dir=meta_dir,
        out_fasta_dir=fasta_dir,
        final_fasta_dir=tmp_path / "final_LGRT1",
        fasta_name="lgrt1_baited.fna.gz",
        owner_inputs={"SPX1": earlier_dir},
        retrieved_at="2026-11-24T00:00:00Z",
    )

    assert rep["distinct_contigs"] == 3
    assert rep["repeated_contigs"] == 1
    assert rep["shared_contigs"] == 1
    assert rep["net_new_contigs"] == 2
    assert rep["shared_by_stage"] == {"SPX1": 1}
    assert rep["circular_contigs"] == 1

    # Verify s2_check_inputs.py passes with 0 problems
    proc = subprocess.run(
        [
            sys.executable,
            str(BSP / "s2_check_inputs.py"),
            "LGRT1",
            str(meta_dir / "genomes.json"),
            str(meta_dir / "objects_staging.json"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, f"s2_check_inputs failed: {proc.stdout}\n{proc.stderr}"
    chk = json.loads(proc.stdout)
    assert chk["problems"] == []
    assert chk["genomes"] == 3
    assert chk["circular_rows"] == 1
