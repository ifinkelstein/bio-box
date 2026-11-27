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


@pytest.fixture(params=["normal-path", "tools-hidden"], autouse=True)
def path_mode(request, monkeypatch):
    """Every test also runs with PATH reduced to /usr/bin:/bin (GitHub Actions has none of the pod's tools)."""
    if request.param == "tools-hidden":
        monkeypatch.setenv("PATH", "/usr/bin:/bin")
    return request.param


def _rand(n: int, seed: int) -> bytes:
    import random

    rng = random.Random(seed)
    return bytes(rng.choice(b"ACGT") for _ in range(n))


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


def test_dec80_an_unconfirmed_or_closed_licence_is_staged_with_redistributable_false(tmp_path):
    """Migrated from the prep's refusal test (DEC-80, 2026-11-26): a licence that no record confirms, or one that
    does not permit redistribution, no longer refuses the staging; the objects say `redistributable: false`."""
    m = _mod()
    fa = tmp_path / "x.fna"
    fa.write_bytes(b">seq1\n" + _rand(400, 1) + b"\n")

    rep = m.stage_ece2(
        [{"source_name": "bad", "source_release": "1", "license": "CC BY-NC 4.0", "redistributable": True,
          "fasta": fa, "output_name": "x.fna.gz", "element_class": "tycheposon"}],
        tmp_path / "E", root=None, root_stages=(), owner_inputs={}, retrieved="2026-11-24T00:00:00Z",
    )
    assert rep["redistributable"] is False
    assert json.loads((tmp_path / "E" / "objects.json").read_text())[0]["redistributable"] is False
    rep = m.stage_mir1(
        [{"source_name": "bad", "source_release": "1", "license": "CC BY 4.0", "redistributable": False,
          "fasta": fa, "output_name": "x.fna.gz", "virus_class": "plv"}],
        tmp_path / "M", root=None, root_stages=(), owner_inputs={}, retrieved="2026-11-24T00:00:00Z",
    )
    assert rep["redistributable"] is False
    assert json.loads((tmp_path / "M" / "objects.json").read_text())[0]["redistributable"] is False
    rep = m.stage_asg1(
        [{"genome_id": "g1", "source_name": "bad", "source_release": "1", "license": "",
          "redistributable": True, "fasta": fa}],
        tmp_path / "A", root=None, root_stages=(), owner_inputs={}, retrieved="2026-11-24T00:00:00Z",
    )
    assert rep["redistributable"] is False
    assert json.loads((tmp_path / "A" / "objects.json").read_text())[0]["redistributable"] is False


# --- card C-70 on the bytes of tower (2026-11-26) -------------------------------------------------------------------


def test_dec65_one_for_metagenome_assembled_circular_elements():
    """DEC-65 (1): an exact 50 bp terminal repeat loses its trailing copy and stays circular; no repeat stays circular
    untrimmed; a repeat longer than 1,000 bp makes the record linear and untrimmed."""
    m = _mod()
    body = _rand(5_000, 2)
    exact = body + body[:50]
    seq, circ, cls, k = m.apply_dec65_metagenomic(exact)
    assert (seq, circ, cls, k) == (body, True, "exact", 50)
    seq, circ, cls, k = m.apply_dec65_metagenomic(body)
    assert (seq, circ, cls, k) == (body, True, "none", 0)
    long_rep = body + body[:1_500]
    seq, circ, cls, k = m.apply_dec65_metagenomic(long_rep)
    assert (seq, circ, cls, k) == (long_rep, False, "long", 0)


def _ece2_fixture(tmp_path: Path):
    """Two Inocle files whose one record is named `sequence1` each (as on tower), the second with an exact 40 bp
    terminal repeat; a pro-tycheposons archive with the element FASTA, its table, host genomes and files that are not
    sequence."""
    import io
    import tarfile

    ino = tmp_path / "dl_ece2"
    ino.mkdir()
    i1 = _rand(3_000, 11)
    i2_body = _rand(3_200, 12)
    (ino / "Inocle_001.fasta").write_bytes(b">sequence1\n" + i1 + b"\n")
    (ino / "Inocle_002.fasta").write_bytes(b">sequence1\n" + i2_body + i2_body[:40] + b"\n")
    (ino / "SHA256SUMS").write_text("not a fasta\n")
    e1, e2 = _rand(1_200, 21), _rand(900, 22)
    host = _rand(4_000, 23)
    files = {
        "LICENSE": b"see LICENSE-MIT and LICENSE-CC-BY-SA\n",
        "elements-pro/pro-623-elements.fna": b">eAAA_00001\n" + e1 + b"\n>eBBB_00002\n" + e2 + b"\n",
        "elements-pro/pro-623-elements.tsv": b"contig_id\tstart\tend\telement_id\tset\n"
                                             b"c1\t0\t1200\teAAA_00001\tTycheposon\nc2\t0\t900\teBBB_00002\tOther\n",
        "data-pro/pro-623-allmaps.fna.xz": b"xz bytes of host genomes",
        "MIT0604-nanopore/MIT0604-2019-NH4.fa": b">MIT0604\n" + host + b"\n",
        "elements-gorg-single-cells/gorg-best-elements.tsv": b"contig_id\tstart\tend\n",
        "YR-phylogeny/yr-all-nr40.faa": b">p\nMKV\n",
    }
    tar = tmp_path / "pro-tycheposons-abc.tar.gz"
    with tarfile.open(tar, "w:gz") as tf:
        for name, data in files.items():
            ti = tarfile.TarInfo(f"pro-tycheposons-abc/{name}")
            ti.size = len(data)
            tf.addfile(ti, io.BytesIO(data))
    return ino, tar, {"i1": i1, "i2": i2_body, "e1": e1, "e2": e2, "host": host}


def test_ece2_from_the_bytes_names_inocles_trims_by_dec65_and_takes_only_the_element_fasta(tmp_path):
    from biobox import ids

    m = _mod()
    ino, tar, s = _ece2_fixture(tmp_path)
    holder = tmp_path / "EGV1"
    holder.mkdir()
    (holder / "genomes.json").write_text(json.dumps([{"contig_uids": [ids.contig_uid(s["e2"]).hex()]}]))
    later = tmp_path / "SAG1"
    later.mkdir()
    (later / "genomes.json").write_text(json.dumps([{"contig_uids": [ids.contig_uid(s["i1"]).hex()]}]))
    out = tmp_path / "meta" / "ECE2"
    rc = m.main([
        "stage-ece2", "--dl", str(ino), "--tycheposons-tar", str(tar), "--out", str(out),
        "--retrieved", "2026-11-25T21:56:05Z", "--owner-inputs", f"EGV1={holder}", "--later", f"SAG1={later}",
        "--final-dir", "/tmp/if2453/s2/inputs/ECE2", "--work", str(tmp_path / "work"),
        "--out-report", str(tmp_path / "rep.json"),
    ])
    assert rc == 0
    rep = json.loads((tmp_path / "rep.json").read_text())
    assert rep["s2_check"]["exit"] == 0 and rep["s2_check"]["problems"] == []
    assert rep["genomes"] == 4 and rep["contigs"] == 4
    assert rep["circular_genomes"] == 2 and rep["trims"] == 1
    assert rep["dec65_classes"] == {"exact": 1, "none": 1}
    assert rep["records_by_set"] == {"tycheposon:Other": 1, "tycheposon:Tycheposon": 1}
    assert rep["owned_by_stage"] == {"EGV1": 1}
    assert rep["later_held_by"] == {"SAG1": 1}
    assert rep["new_contigs_after_later"] == 2
    assert rep["redistributable"] is True

    ino_fa = (out / "ece2_inocles.fna").read_bytes()
    assert b">Inocle_001|sequence1\n" in ino_fa and b">Inocle_002|sequence1\n" in ino_fa
    genomes = json.loads((out / "genomes.json").read_text())
    uids = {u for g in genomes for u in g["contig_uids"]}
    assert ids.contig_uid(s["i2"]).hex() in uids  # the trimmed circle
    assert ids.contig_uid(s["host"]).hex() not in uids
    circ = {u for g in genomes for u in g["circular"]}
    assert circ == {ids.contig_uid(s["i1"]).hex(), ids.contig_uid(s["i2"]).hex()}
    assert all(g["kind"] == "plasmid" for g in genomes)

    with open(out / "skipped_members.tsv", newline="") as fh:
        skipped = {r["member"].split("/", 1)[1]: r["reason"] for r in csv.DictReader(fh, delimiter="\t")}
    assert skipped["data-pro/pro-623-allmaps.fna.xz"].startswith("host_genomes_not_the_subject")
    assert skipped["MIT0604-nanopore/MIT0604-2019-NH4.fa"].startswith("host_genome_not_the_subject")
    assert skipped["elements-gorg-single-cells/gorg-best-elements.tsv"].startswith("coordinates_only")
    assert "elements-pro/pro-623-elements.fna" not in skipped

    objs = json.loads((out / "objects.json").read_text())
    staging = json.loads((out / "objects_staging.json").read_text())
    assert sorted(Path(o["path"]).name for o in objs) == ["ece2_inocles.fna", "ece2_tycheposons.fna"]
    assert all(o["path"].startswith("/tmp/if2453/s2/inputs/ECE2/") for o in objs)
    assert all(Path(o["path"]).parent == out for o in staging)
    with open(out / "later_records.tsv", newline="") as fh:
        assert [r["holder"] for r in csv.DictReader(fh, delimiter="\t")] == ["SAG1"]
    assert "objects_staging.json" not in (out / "inputs_SUMS").read_text()


def test_ece2_from_the_bytes_is_byte_identical_twice_and_refuses_an_archive_without_the_elements(tmp_path):
    import io
    import tarfile

    m = _mod()
    ino, tar, _ = _ece2_fixture(tmp_path)
    outs = []
    for k in (1, 2):
        out = tmp_path / f"run{k}" / "ECE2"
        assert m.main(["stage-ece2", "--dl", str(ino), "--tycheposons-tar", str(tar), "--out", str(out),
                       "--retrieved", "2026-11-25T21:56:05Z", "--final-dir", "/tmp/if2453/s2/inputs/ECE2",
                       "--work", str(tmp_path / f"w{k}")]) == 0
        outs.append(out)
    for name in ("genomes.json", "objects.json", "records.tsv", "ece2_metadata.tsv", "ece2_inocles.fna",
                 "ece2_tycheposons.fna", "trims.tsv", "skipped_members.tsv", "inputs_SUMS"):
        assert (outs[0] / name).read_bytes() == (outs[1] / name).read_bytes(), name

    bad = tmp_path / "bad.tar.gz"
    with tarfile.open(bad, "w:gz") as tf:
        ti = tarfile.TarInfo("x/README.org")
        ti.size = 3
        tf.addfile(ti, io.BytesIO(b"abc"))
    with pytest.raises(m.Refused):
        m.main(["stage-ece2", "--dl", str(ino), "--tycheposons-tar", str(bad), "--out", str(tmp_path / "B" / "ECE2"),
                "--retrieved", "x", "--work", str(tmp_path / "wb")])
    assert not (tmp_path / "B" / "ECE2").exists()


def test_fasta_dir_moves_only_a_fasta_of_64_mib_and_more_and_objects_staging_names_where_each_is(tmp_path,
                                                                                                 monkeypatch):
    from biobox import ids

    m = _mod()
    monkeypatch.setattr(m, "STOR_MIN_BYTES", 1_000)  # the fixture's stand-in for 64 MiB
    big_fa, small_fa = tmp_path / "x.fna", tmp_path / "y.fna"
    big_fa.write_bytes(b">e1\n" + _rand(1_200, 31) + b"\n")
    small_fa.write_bytes(b">e2\n" + _rand(300, 32) + b"\n")
    out, stor = tmp_path / "meta" / "ECE2", tmp_path / "stor" / "ECE2"
    src = {"source_name": "s", "source_release": "1", "license": "CC BY 4.0", "redistributable": True,
           "element_class": "tycheposon"}
    rep = m.stage_ece2(
        [{**src, "fasta": big_fa, "output_name": "ece2_x.fna"}, {**src, "fasta": small_fa, "output_name": "ece2_y.fna"}],
        out, root=None, root_stages=(), owner_inputs={}, retrieved="t", final_dir=Path("/tmp/if2453/s2/inputs/ECE2"),
        fasta_dir=stor,
    )
    assert (stor / "ece2_x.fna").is_file() and not (out / "ece2_x.fna").exists()
    assert (out / "ece2_y.fna").is_file() and not (stor / "ece2_y.fna").exists()
    assert not list(stor.glob(".*partial"))
    staged = {Path(o["path"]).name: o["path"] for o in json.loads((out / "objects_staging.json").read_text())}
    assert staged == {"ece2_x.fna": str(stor / "ece2_x.fna"), "ece2_y.fna": str(out / "ece2_y.fna")}
    assert rep["sha256"]["ece2_x.fna"] == ids.sha256_file(stor / "ece2_x.fna").hex()


def test_goev_names_and_classes():
    m = _mod()
    assert m.goev_genome_name("GOEV_Genomic_DB/TARA_MED_NCLDV_00046-fixed.fa") == "TARA_MED_NCLDV_00046"
    assert m.goev_genome_name("GOEV_Genomic_DB/Schulz_GVMAG_M_3300009436_15-contigs.fa") == "Schulz_GVMAG_M_3300009436_15"
    assert m.goev_genome_name("GOEV_Genomic_DB/Moniruzzaman_MM05_ERX555935_132_dc.fa") == "Moniruzzaman_MM05_ERX555935_132"
    assert m.goev_class("Schulz_GVMAG_M_3300009436_15") == "environmental"
    assert m.goev_class("Moniruzzaman_MM05_ERX555935_132") == "environmental"
    assert m.goev_class("TARA_MED_NCLDV_00046") == "environmental"
    assert m.goev_class("Poxviridae_LT896721_Cowpox_virus_isolate_Ger_2010_MKY") == "Poxviridae"
    assert m.goev_class("Asfarviridae_AY261365_African_swine_fever_virus_isolate_Warmbaths") == "Asfarviridae"
    assert m.goev_class("Mimiviridae_KC008572_UNVERIFIED_Moumouvirus_goulette") == "Mimiviridae"


def _mir1_fixture(tmp_path: Path):
    """The ten files of tower_mir1 in small: two Mirusviricota genomes (Medvedeva 2025), two Tara Mirusviricota MAGs and
    a near-complete genome (Gaia 2023), a GOEV database whose Tara MAG repeats Gaia's names and bytes, a GVMAG byte-equal
    to a GV1 contig, a Moniruzzaman GVMAG, a Mimiviridae reference and the Poxviridae and Asfarviridae references, two
    EVE loci, a FLYE PLV file whose names hold a space, and the files that are not genomes."""
    import gzip as gz
    import io
    import tarfile
    import zipfile

    small, large = tmp_path / "dl_small", tmp_path / "dl_large"
    small.mkdir()
    large.mkdir()
    s = {k: _rand(n, seed) for seed, (k, n) in enumerate([
        ("m1a", 900), ("m1b", 700), ("m2", 800), ("t5a", 600), ("t5b", 650), ("t7", 500), ("nc", 1500),
        ("gv", 1000), ("mon", 700), ("mimi", 1100), ("pox", 900), ("asf", 950), ("eve1", 1300), ("eve2", 1200),
        ("fl1", 800), ("fl2", 700)], start=100)}

    def rec(name: str, key: str) -> bytes:
        return b">" + name.encode() + b"\n" + s[key] + b"\n"

    def zwrite(path: Path, members: dict[str, bytes]) -> None:
        with zipfile.ZipFile(path, "w") as zf:
            for name, data in members.items():
                zf.writestr(name, data)

    zwrite(large / "medvedeva2025_01_Genomes.zip", {
        "01_Genomes/MIRUS_G_0001.fa": rec("MIRUS_G_0001_c0001", "m1a") + rec("MIRUS_G_0001_c0002", "m1b"),
        "01_Genomes/MIRUS_G_0002.fa": rec("MIRUS_G_0002_c0001", "m2"),
        "__MACOSX/01_Genomes/._MIRUS_G_0001.fa": b"fork",
    })
    zwrite(small / "gaia2023_Mirusviricota_Genomes.zip", {
        "Mirus_111_MAGs.fa": rec("TARA_AON_NCLDV_00005_000000000001", "t5a")
                             + rec("TARA_AON_NCLDV_00005_000000000002", "t5b")
                             + rec("TARA_PSW_NCLDV_00007_000000000001", "t7"),
        "Mirus_near_complete_genome.fa": rec("TARA_SAMEA2619782_METAG_scaffold_2_length_432583", "nc"),
    })
    zwrite(large / "gaia2023_GOEV_Genomic_DB.zip", {
        "GOEV_Genomic_DB/TARA_AON_NCLDV_00005-fixed.fa": rec("TARA_AON_NCLDV_00005_000000000001", "t5a")
                                                          + rec("TARA_AON_NCLDV_00005_000000000002", "t5b"),
        "GOEV_Genomic_DB/Schulz_GVMAG_M_1_2-contigs.fa": rec("Schulz_GVMAG_M_1_2_000000000001", "gv"),
        "GOEV_Genomic_DB/Moniruzzaman_MM01_ERX1_3_dc.fa": rec("Moniruzzaman_MM01_ERX1_3_c1", "mon"),
        "GOEV_Genomic_DB/Mimiviridae_KC1_Mimivirus_x-contigs.fa": rec("KC1_1", "mimi"),
        "GOEV_Genomic_DB/Poxviridae_LT1_Cowpox_virus-contigs.fa": rec("LT1_1", "pox"),
        "GOEV_Genomic_DB/Asfarviridae_AY1_African_swine_fever_virus-contigs.fa": rec("AY1_1", "asf"),
        "__MACOSX/GOEV_Genomic_DB/._Poxviridae_LT1_Cowpox_virus-contigs.fa": b"fork",
    })
    (large / "bellas2023_mcp-loci-hq.fna.gz").write_bytes(gz.compress(rec("MJMQ01_1408", "eve1") + rec("BJOX01_1121", "eve2"), mtime=0))
    (small / "bellas2023_Synura_10000bp_FLYE_assembly_PLV.fasta").write_bytes(rec("Contig 317", "fl1") + rec("Contig 326", "fl2"))
    (small / "bellas2023_Readme.txt").write_text("Readme\n")
    zwrite(large / "medvedeva2025_Supplemental_Tables.zip", {"Table_S03_Mirusviricota_Genomes.xlsx": b"xlsx"})
    with tarfile.open(small / "bellas2023_MCPs_group_by_genome.fasta.tar.gz", "w:gz") as tf:
        ti = tarfile.TarInfo("MCP_grouped_by_genome/A_MCP_extract.fasta")
        ti.size = 8
        tf.addfile(ti, io.BytesIO(b">a\nACGT\n"))
    return small, large, s


def _mir1_args(small, large, out, work, *extra):
    return ["stage-mir1", "--dl", str(small), "--large-dir", str(large), "--out", str(out),
            "--retrieved", "2026-11-25T21:56:55Z", "--final-dir", "/tmp/if2453/s2/inputs/MIR1", "--work", str(work),
            *extra]


def test_mir1_from_the_bytes_splits_goev_leaves_out_pox_and_asfar_and_lets_gv1_own_its_contigs(tmp_path):
    from biobox import ids

    m = _mod()
    small, large, s = _mir1_fixture(tmp_path)
    gv1 = tmp_path / "GV1"
    gv1.mkdir()
    (gv1 / "genomes.json").write_text(json.dumps([{"contig_uids": [ids.contig_uid(s["gv"]).hex()]}]))
    out = tmp_path / "meta" / "MIR1"
    assert m.main([*_mir1_args(small, large, out, tmp_path / "work"), "--owner-inputs", f"GV1={gv1}",
                   "--out-report", str(tmp_path / "rep.json")]) == 0
    rep = json.loads((tmp_path / "rep.json").read_text())
    assert rep["s2_check"]["exit"] == 0 and rep["s2_check"]["problems"] == []
    # Medvedeva 2, GOEV 4 (Tara, Schulz, Moniruzzaman; Mimi), Gaia 2 of 3 (its Tara MAG repeats GOEV's contig set
    # exactly, so it is listed once, under GOEV, which comes first: DEC-27), 2 EVE loci, 2 FLYE contigs.
    assert rep["genomes_by_class"] == {"goev_environmental": 3, "goev_reference": 1, "mirus_gaia2023": 2,
                                       "mirusviricota": 2, "plv_eve": 2, "plv_flye": 2}
    assert rep["repeated_records"] == 1
    assert rep["owned_by_stage"] == {"GV1": 1}
    assert rep["held_genomes"] == {"tower_mir1:gaia2023_GOEV_Genomic_DB:environmental:wholly:GV1": 1}
    assert rep["left_out"]["isolate_reference_Poxviridae"] == 1
    assert rep["left_out"]["isolate_reference_Asfarviridae"] == 1
    assert rep["redistributable"] is False
    objs = {Path(o["path"]).name: o for o in json.loads((out / "objects.json").read_text())}
    assert objs["mir1_bellas2023_eve_loci.fna"]["redistributable"] is False
    assert objs["mir1_bellas2023_flye_plv.fna"]["redistributable"] is False
    assert objs["mir1_medvedeva2025.fna"]["redistributable"] is True
    uids = {u for g in json.loads((out / "genomes.json").read_text()) for u in g["contig_uids"]}
    assert ids.contig_uid(s["pox"]).hex() not in uids and ids.contig_uid(s["asf"]).hex() not in uids
    assert ids.contig_uid(s["gv"]).hex() in uids  # owned by GV1 by the chain, not removed
    gaia = (out / "mir1_gaia2023_mirus.fna").read_text()
    assert ">gaia2023_mirus|TARA_AON_NCLDV_00005_000000000001\n" in gaia
    flye = (out / "mir1_bellas2023_flye_plv.fna").read_text()
    assert ">bellas2023_Synura_10000bp_FLYE_assembly_PLV|Contig_317\n" in flye
    with open(out / "left_out.tsv", newline="") as fh:
        reasons = {(r["file"], r["member"]): r["reason"] for r in csv.DictReader(fh, delimiter="\t")}
    assert reasons[("bellas2023_MCPs_group_by_genome.fasta.tar.gz", "")].startswith("not_genome_sequence")
    assert reasons[("medvedeva2025_Supplemental_Tables.zip", "")].startswith("non_sequence_file")
    genomes = json.loads((out / "genomes.json").read_text())
    assert all(g["kind"] == "virus" and g["circular"] == [] for g in genomes)


def test_mir1_exclude_leaves_a_source_out_and_two_runs_are_byte_identical(tmp_path):
    m = _mod()
    small, large, _ = _mir1_fixture(tmp_path)
    outs = []
    for k in (1, 2):
        out = tmp_path / f"run{k}" / "MIR1"
        assert m.main([*_mir1_args(small, large, out, tmp_path / f"w{k}"), "--exclude", "bellas2023_eve_loci"]) == 0
        outs.append(out)
    names = {Path(o["path"]).name for o in json.loads((outs[0] / "objects.json").read_text())}
    assert "mir1_bellas2023_eve_loci.fna" not in names
    with open(outs[0] / "left_out.tsv", newline="") as fh:
        assert ("bellas2023_mcp-loci-hq.fna.gz", "excluded_by_option") in {
            (r["file"], r["reason"]) for r in csv.DictReader(fh, delimiter="\t")}
    for name in ("genomes.json", "objects.json", "records.tsv", "mir1_metadata.tsv", "mir1_medvedeva2025.fna",
                 "mir1_goev_environmental.fna", "left_out.tsv", "inputs_SUMS"):
        assert (outs[0] / name).read_bytes() == (outs[1] / name).read_bytes(), name


def test_mir1_refuses_a_record_name_twice_in_the_stage(tmp_path):
    m = _mod()
    fa1, fa2 = tmp_path / "a.fna", tmp_path / "b.fna"
    fa1.write_bytes(b">same\n" + _rand(300, 41) + b"\n")
    fa2.write_bytes(b">same\n" + _rand(300, 42) + b"\n")
    src = {"source_release": "1", "license": "CC BY 4.0", "redistributable": True, "virus_class": "plv"}
    with pytest.raises(m.Refused):
        m.stage_mir1([{**src, "source_name": "a", "fasta": fa1, "output_name": "a.fna"},
                      {**src, "source_name": "b", "fasta": fa2, "output_name": "b.fna"}],
                     tmp_path / "M", root=None, root_stages=(), owner_inputs={}, retrieved="t")
    assert not (tmp_path / "M").exists()


_ROOT13 = ("S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1")
_R26_11_4 = ("V2", "PL2", "PL3", "CG1", "MGE1", "SAT1", "ECE1", "O1", "O2", "LPL1", "SYM1", "LRV1", "LGRT1", "MGE2",
             "MV1", "V3", "LR1", "GART1", "GART2")


@pytest.mark.parametrize(("stage", "note", "holders"), [
    ("ECE2", "planning/measurements/2026-11-26_c70_tier1_inputs.md", (*_ROOT13, *_R26_11_4, "GV1", "EGV1")),
])
def test_c70_template_has_real_sums_found_in_the_note_and_expect_lines_over_its_holders(stage, note, holders):
    se = _stage_est()
    conf_path = BSP / "stages" / "templates" / f"{stage}.conf"
    c = se.read_conf(conf_path, stage=stage, tree=REPO)
    assert note in c.notes
    text = conf_path.read_text()
    sums = [ln.split()[0][4:] for ln in text.splitlines() if ln.startswith("SUM=")]
    assert sums and all(s != "0" * 64 for s in sums)
    note_text = (REPO / note).read_text()
    assert all(s in note_text for s in sums)
    expect = [ln.split()[1] for ln in text.splitlines() if ln.startswith("EXPECT=held_by")]
    assert tuple(expect) == holders


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
