"""Ring R1: Card C-70 — Tier-1 compact high-yield MGE (`ECE2`), divergent dsDNA virus (`MIR1`),
and post-GlobDB *Asgardarchaeota* (`ASG1`) pre-build staging, topology, dedup, and `.conf` validation.

Tests:
- License redistributability gate (`DEC-51`) across all three stages.
- `ECE2` topology rules (`Inocles` circular preserved; `Tycheposons` always linear untrimmed per `DEC-65`;
  `kind="plasmid"`, `0` trims, exact TIR detection, within-stage `DEC-27` and cross-stage `DEC-26` ownership).
- `MIR1` topology rules (*Mirusviricota* multi-contig bins + `PLVs` single-contig records, always linear
  untrimmed per `DEC-65`, `kind="virus"`, `0` trims, ownership tracking).
- `ASG1` staging (`kind="prokaryote"`, `genetic_code=11`, `skani >= 96%` ANI / `>= 50%` AF parser &
  reference-dedup filter against `G1`/GTDB archaea, `asg1_skani_dropped.tsv`).
- Determinism: identical inputs in reversed order produce byte-identical `genomes.json`, `objects.json`,
  metadata sidecars, and `.fna.gz` (`mtime=0`).
- `.conf` validation via `stage_est.read_conf(..., tree=REPO)` for `ECE2.conf`, `MIR1.conf`, `ASG1.conf`
  verifying `DEPTH=annotated`, `LIBRARY=pfam 38.2 gathering`, `CRISPR=yes`, and `EGV1 -> ECE2 -> MIR1 -> ASG1`.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-70")

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"


def _load(name: str, filename: str):
    path = BSP / filename
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None, f"missing {path}"
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _mod():
    return _load("c70_tier1_mge_dsdna_asgard", "c70_tier1_mge_dsdna_asgard.py")


def _stage_est():
    return _load("stage_est", "stage_est.py")


def _revcomp(seq: bytes) -> bytes:
    return seq.translate(bytes.maketrans(b"ACGTacgt", b"TGCAtgca"))[::-1]


def test_license_gate_enforces_dec51_redistributability():
    m = _mod()
    for ok in ("CC BY 4.0", "CC-BY-4.0", "CC BY-SA 4.0", "CC0 1.0", "MIT", "INSDC (public)"):
        assert m.is_redistributable_license(ok) is True, ok
    for bad in ("CC BY-NC 4.0", "CC-BY-NC-SA 4.0", "restricted", "All rights reserved", ""):
        assert m.is_redistributable_license(bad) is False, bad


def test_classify_ece2_and_mir1_topology_rules():
    m = _mod()
    # Inocles default circular unless explicitly linear/partial/fragment
    ino_circ = m.classify_ece2_record("INOCLE_001 Streptococcus salivarius closed circular ECE", "inocle")
    assert ino_circ["circular"] is True
    assert ino_circ["kind"] == "plasmid"

    ino_frag = m.classify_ece2_record("INOCLE_002 partial contig", "inocle")
    assert ino_frag["circular"] is False

    ino_lin_topo = m.classify_ece2_record("INOCLE_003 closed", "inocle", topology="linear")
    assert ino_lin_topo["circular"] is False

    # Tycheposons are ALWAYS linear even if topology="circular" is passed
    tych = m.classify_ece2_record("TYCHE_001 Prochlorococcus MIT0604", "tycheposon", topology="circular")
    assert tych["circular"] is False
    assert tych["kind"] == "plasmid"

    # MIR1 records are ALWAYS linear kind="virus"
    for vcls in ("mirusviricota", "mirus_goev", "plv", "plv_eve"):
        rec = m.classify_mir1_record(f"{vcls}_001 complete circular?", vcls)
        assert rec["circular"] is False
        assert rec["kind"] == "virus"


def test_stage_ece2_circular_inocles_linear_tycheposons_and_ownership(tmp_path):
    from biobox import ids

    m = _mod()
    ino_seq = b"ATGCGTACGTAGCTAGCTAG" * 80
    tir = b"ACGTTGCAAACCGGTTACGTTGCAAACCGGTT"  # 32 bp
    tych_seq = tir + (b"GCTAGCTAGCATCGATCGAT" * 60) + _revcomp(tir)

    fa_ino = tmp_path / "inocles.fna"
    fa_ino.write_bytes(
        b">INOCLE_01 Streptococcus oral giant circular ECE\n" + ino_seq + b"\n"
        + b">INOCLE_01_DUP duplicate circular ECE\n" + ino_seq + b"\n"
    )
    fa_tych = tmp_path / "tycheposons.fna"
    fa_tych.write_bytes(b">TYCHE_01 Prochlorococcus integrative element\n" + tych_seq + b"\n")

    owner_egv1 = tmp_path / "EGV1"
    owner_egv1.mkdir()
    tych_cu = ids.contig_uid(tych_seq).hex()
    (owner_egv1 / "genomes.json").write_text(json.dumps([{"contig_uids": [tych_cu]}]))

    sources = [
        {
            "source_name": "tycheposons_hackl_2023",
            "source_release": "github:thackl/pro-tycheposons",
            "license": "CC BY-SA 4.0",
            "redistributable": True,
            "fasta": fa_tych,
            "output_name": "ece2_tycheposons.fna.gz",
            "element_class": "tycheposon",
        },
        {
            "source_name": "inocles_kiguchi_2025",
            "source_release": "10.5281/zenodo.16656823",
            "license": "CC BY 4.0",
            "redistributable": True,
            "fasta": fa_ino,
            "output_name": "ece2_inocles.fna.gz",
            "element_class": "inocle",
        },
    ]

    out_dir = tmp_path / "ECE2"
    rep = m.stage_ece2(
        sources, out_dir, root=None, root_stages=(),
        owner_inputs={"EGV1": owner_egv1}, retrieved="2026-11-24T00:00:00Z",
    )

    assert rep["trims"] == 0
    assert rep["records"] == 3
    assert rep["genomes"] == 2
    assert rep["repeated_records"] == 1
    assert rep["circular_genomes"] == 1
    assert rep["linear_genomes"] == 1
    assert rep["genomes_with_tir"] == 1
    assert rep["held_contigs"] == 1
    assert rep["new_contigs"] == 1
    assert rep["owned_by_stage"] == {"EGV1": 1}

    genomes = json.loads((out_dir / "genomes.json").read_text())
    by_cu = {g["contig_uids"][0]: g for g in genomes}
    ino_cu = ids.contig_uid(ino_seq).hex()
    assert by_cu[ino_cu]["circular"] == [ino_cu]
    assert by_cu[ino_cu]["kind"] == "plasmid"
    assert by_cu[ino_cu]["assembly_quality"] == "complete"
    assert by_cu[tych_cu]["circular"] == []
    assert by_cu[tych_cu]["kind"] == "plasmid"

    with open(out_dir / "ece2_metadata.tsv", newline="") as fh:
        rows = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert rows["INOCLE_01"]["circular"] == "true"
    assert rows["TYCHE_01"]["circular"] == "false"
    assert int(rows["TYCHE_01"]["tir_length"]) == len(tir)
    assert rows["TYCHE_01"]["owner_stage"] == "EGV1"


def test_stage_mir1_multi_contig_mirus_and_single_contig_plvs(tmp_path):
    from biobox import ids

    m = _mod()
    c1 = b"AACCGGTTAACCGGTT" * 80
    c2 = b"TTGGCCAATTGGCCAA" * 60
    plv_seq = b"ACGTTGCAAACCGGTTACGTTGCAAA" + (b"CATGCATGCATGCATG" * 50) + _revcomp(b"ACGTTGCAAACCGGTTACGTTGCAAA")

    mirus_bin = tmp_path / "MIRUS_MAG_0001.fna"
    mirus_bin.write_bytes(b">mirus_c1\n" + c1 + b"\n>mirus_c2\n" + c2 + b"\n")
    plv_fa = tmp_path / "plvs.fna"
    plv_fa.write_bytes(b">PLV_GKS_001 aquatic PLV\n" + plv_seq + b"\n")

    sources = [
        {
            "source_name": "mirusviricota_medvedeva_2025",
            "source_release": "10.6084/m9.figshare.28955240",
            "license": "CC BY 4.0",
            "redistributable": True,
            "mode": "multi_contig_bins",
            "bin_fastas": [("MIRUS_MAG_0001", mirus_bin)],
            "output_name": "mir1_mirusviricota_2025.fna.gz",
            "virus_class": "mirusviricota",
        },
        {
            "source_name": "plv_bellas_2021",
            "source_release": "figshare.s.9a7a1d16d77ea9d658a1",
            "license": "CC BY 4.0",
            "redistributable": True,
            "mode": "single_contig",
            "fasta": plv_fa,
            "output_name": "mir1_plv_bellas_2021.fna.gz",
            "virus_class": "plv",
        },
    ]

    out_dir = tmp_path / "MIR1"
    rep = m.stage_mir1(
        sources, out_dir, root=None, root_stages=(),
        owner_inputs={}, retrieved="2026-11-24T00:00:00Z",
    )

    assert rep["trims"] == 0
    assert rep["genomes"] == 2
    assert rep["contigs"] == 3
    assert rep["circular_genomes"] == 0
    assert rep["linear_genomes"] == 2
    assert rep["genomes_by_class"] == {"mirusviricota": 1, "plv": 1}

    genomes = json.loads((out_dir / "genomes.json").read_text())
    assert all(g["kind"] == "virus" and g["circular"] == [] for g in genomes)
    lens = sorted(len(g["contig_uids"]) for g in genomes)
    assert lens == [1, 2]
    assert ids.contig_uid(plv_seq).hex() in {u for g in genomes for u in g["contig_uids"]}


def test_asg1_skani_parser_and_dedup_staging(tmp_path):
    m = _mod()
    mag_keep = tmp_path / "ASG_MAG_001.fna"
    mag_drop = tmp_path / "ASG_MAG_002.fna"
    mag_keep.write_bytes(b">asg1_c1\n" + (b"ATGCGTACGT" * 120) + b"\n>asg1_c2\n" + (b"GCTAGCATCG" * 100) + b"\n")
    mag_drop.write_bytes(b">asg2_c1\n" + (b"TTGGCCAATT" * 100) + b"\n")

    skani_tsv = tmp_path / "skani_vs_g1.tsv"
    skani_tsv.write_text(
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        f"/ref/G1_ARCH_001.fna\t{mag_drop}\t97.42\t62.1\t58.3\tG1_ARCH_001\tASG_MAG_002\n"
        f"/ref/G1_ARCH_002.fna\t{mag_keep}\t91.20\t70.0\t65.0\tG1_ARCH_002\tASG_MAG_001\n",
        encoding="utf-8",
    )

    matches = m.parse_skani_matches(skani_tsv, min_ani=96.0, min_af=50.0)
    assert list(matches) == ["ASG_MAG_002"]
    assert matches["ASG_MAG_002"]["ref"] == "G1_ARCH_001"

    entries = [
        {
            "genome_id": "ASG_MAG_002",
            "source_name": "appler_2026_asgard",
            "source_release": "10.6084/m9.figshare.26426707",
            "license": "CC BY 4.0",
            "redistributable": True,
            "assembly_quality": "mag_hq",
            "fasta": mag_drop,
        },
        {
            "genome_id": "ASG_MAG_001",
            "source_name": "appler_2026_asgard",
            "source_release": "10.6084/m9.figshare.26426707",
            "license": "CC BY 4.0",
            "redistributable": True,
            "assembly_quality": "mag_hq",
            "fasta": mag_keep,
        },
    ]

    out_dir = tmp_path / "ASG1"
    rep = m.stage_asg1(
        entries, out_dir, root=None, root_stages=(),
        owner_inputs={}, retrieved="2026-11-24T00:00:00Z",
        skani_ref_matches=matches, skani_min_ani=96.0,
    )

    assert rep["genomes"] == 1
    assert rep["contigs"] == 2
    assert rep["skani_dropped"] == 1
    assert rep["trims"] == 0
    genomes = json.loads((out_dir / "genomes.json").read_text())
    assert len(genomes) == 1
    assert genomes[0]["kind"] == "prokaryote"
    assert genomes[0]["genetic_code"] == 11
    assert genomes[0]["assembly_quality"] == "mag_hq"

    with open(out_dir / "asg1_skani_dropped.tsv", newline="") as fh:
        dropped = list(csv.DictReader(fh, delimiter="\t"))
    assert len(dropped) == 1
    assert dropped[0]["genome_id"] == "ASG_MAG_002"
    assert dropped[0]["ref_match"] == "G1_ARCH_001"


def test_determinism_independent_of_input_order(tmp_path):
    m = _mod()
    fa_a = tmp_path / "a.fna"
    fa_b = tmp_path / "b.fna"
    fa_a.write_bytes(b">INOCLE_A\n" + (b"ACGTACGTAC" * 50) + b"\n")
    fa_b.write_bytes(b">TYCHE_B\n" + (b"TGCATGCATG" * 50) + b"\n")
    src_a = {
        "source_name": "inocles_kiguchi_2025",
        "source_release": "10.5281/zenodo.16656823",
        "license": "CC BY 4.0",
        "redistributable": True,
        "fasta": fa_a,
        "output_name": "ece2_inocles.fna.gz",
        "element_class": "inocle",
    }
    src_b = {
        "source_name": "tycheposons_hackl_2023",
        "source_release": "github:thackl/pro-tycheposons",
        "license": "CC BY-SA 4.0",
        "redistributable": True,
        "fasta": fa_b,
        "output_name": "ece2_tycheposons.fna.gz",
        "element_class": "tycheposon",
    }

    d1 = tmp_path / "run1"
    d2 = tmp_path / "run2"
    final_p = Path("/tmp/if2453/s2/inputs/ECE2")
    m.stage_ece2([src_a, src_b], d1, root=None, root_stages=(), owner_inputs={}, retrieved="2026-11-24T00:00:00Z", final_dir=final_p)
    m.stage_ece2([src_b, src_a], d2, root=None, root_stages=(), owner_inputs={}, retrieved="2026-11-24T00:00:00Z", final_dir=final_p)

    for fname in ("genomes.json", "objects.json", "ece2_metadata.tsv", "ece2_inocles.fna.gz", "ece2_tycheposons.fna.gz", "inputs_SUMS"):
        assert (d1 / fname).read_bytes() == (d2 / fname).read_bytes(), f"non-deterministic {fname}"


def test_refuses_non_redistributable_sources_across_all_stages(tmp_path):
    m = _mod()
    fa = tmp_path / "x.fna"
    fa.write_bytes(b">seq1\n" + (b"ACGT" * 50) + b"\n")

    with pytest.raises(m.Refused):
        m.stage_ece2(
            [{"source_name": "bad", "source_release": "1", "license": "CC BY-NC 4.0", "redistributable": True,
              "fasta": fa, "output_name": "x.fna.gz", "element_class": "inocle"}],
            tmp_path / "E", root=None, root_stages=(), owner_inputs={}, retrieved="2026-11-24T00:00:00Z",
        )
    with pytest.raises(m.Refused):
        m.stage_mir1(
            [{"source_name": "bad", "source_release": "1", "license": "CC BY 4.0", "redistributable": False,
              "fasta": fa, "output_name": "x.fna.gz", "virus_class": "plv"}],
            tmp_path / "M", root=None, root_stages=(), owner_inputs={}, retrieved="2026-11-24T00:00:00Z",
        )
    with pytest.raises(m.Refused):
        m.stage_asg1(
            [{"genome_id": "g1", "source_name": "bad", "source_release": "1", "license": "restricted",
              "redistributable": True, "fasta": fa}],
            tmp_path / "A", root=None, root_stages=(), owner_inputs={}, retrieved="2026-11-24T00:00:00Z",
        )


@pytest.mark.parametrize(
    ("stage", "expected_after"),
    [("ECE2", ("EGV1",)), ("MIR1", ("ECE2",)), ("ASG1", ("MIR1",))],
)
def test_c70_stage_confs_validate_against_stage_est(stage: str, expected_after: tuple[str, ...]):
    se = _stage_est()
    conf_path = BSP / "stages" / "templates" / f"{stage}.conf"
    c = se.read_conf(conf_path, stage=stage, tree=REPO)
    assert c.stage == stage
    assert c.depth == "annotated"
    assert c.library == ("pfam", "38.2", "gathering")
    assert c.crispr is True
    assert c.after == expected_after
    assert c.inputs == Path(f"/tmp/if2453/s2/inputs/{stage}")
    assert "planning/measurements/2026-11-24_c70_tier1_prep.md" in c.notes
