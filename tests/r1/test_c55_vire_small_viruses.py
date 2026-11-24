"""Ring R1: card C-55 (`SAT1` + `V3` VIRE) parsing, INPHARED delta filtering, DEC-65 DTR trimming,
no-5-kb-floor small circular DNA viruses/satellites, redistributable check, contig_uid ownership
deduplication, and determinism.
"""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import io
import json
import sys
import tarfile
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-55")

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"


def _load_c55():
    spec = importlib.util.spec_from_file_location("c55_vire_small_viruses", BSP / "c55_vire_small_viruses.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c55_vire_small_viruses"] = mod
    spec.loader.exec_module(mod)
    return mod


def _make_seq(length: int, seed: int = 1) -> bytes:
    """Deterministic pseudo-random DNA sequence of given length (not low-complexity)."""
    bases = b"ACGT"
    x = seed & 0xFFFFFFFF
    out = bytearray(length)
    for i in range(length):
        x = (1664525 * x + 1013904223) & 0xFFFFFFFF
        out[i] = bases[(x >> 24) & 3]
    return bytes(out)


def test_inphared_delta_filters_against_s1_accessions(tmp_path):
    mod = _load_c55()
    s1_tsv = tmp_path / "14Apr2025_data.tsv.gz"
    new_tsv = tmp_path / "9Aug2025_data.tsv.gz"
    with gzip.open(s1_tsv, "wt") as fh:
        fh.write("Accession\tDescription\tGenome Length (bp)\n")
        fh.write("NC_001416\tEscherichia phage Lambda\t48502\n")
        fh.write("ON123456.1\tBacillus phage Old\t120000\n")
    with gzip.open(new_tsv, "wt") as fh:
        fh.write("Accession\tDescription\tGenome Length (bp)\n")
        fh.write("NC_001416.1\tEscherichia phage Lambda\t48502\n")
        fh.write("ON123456\tBacillus phage Old\t120000\n")
        fh.write("PX530545.1\tStaphylococcus phage NODE_1\t209876\n")
        fh.write("PQ999999\tPseudomonas phage New\t65000\n")

    delta = mod.filter_inphared_delta(new_tsv, s1_tsv)
    assert delta == {"PX530545", "PQ999999"}


def test_parse_genbank_and_classify_refseq_subsets_dna_only():
    mod = _load_c55()
    cress_seq = _make_seq(2100, seed=11).decode()
    arch_seq = _make_seq(15000, seed=12).decode()
    sat_seq = _make_seq(11624, seed=13).decode()
    rna_seq = _make_seq(1800, seed=14).decode()

    gb = (
        f"LOCUS       NC_099001               2100 bp    ss-DNA  circular VRL 10-NOV-2026\n"
        f"DEFINITION  Circovirus sp. isolate test, complete genome.\n"
        f"VERSION     NC_099001.1\n"
        f"SOURCE      Circovirus sp.\n"
        f"  ORGANISM  Circovirus sp.\n"
        f"            Viruses; Monodnaviria; Shotokuvirae; Cressdnaviricota; Arfiviricetes.\n"
        f"ORIGIN\n"
        f"        1 {cress_seq.lower()}\n"
        f"//\n"
        f"LOCUS       NC_099002              15000 bp    ds-DNA  circular VRL 10-NOV-2026\n"
        f"DEFINITION  Sulfolobus spindle-shaped virus 1, complete genome.\n"
        f"VERSION     NC_099002.1\n"
        f"SOURCE      Sulfolobus spindle-shaped virus 1\n"
        f"  ORGANISM  Sulfolobus spindle-shaped virus 1\n"
        f"            Viruses; Fuselloviridae; Alphafusellovirus.\n"
        f"FEATURES             Location/Qualifiers\n"
        f"     source          1..15000\n"
        f"                     /host=\"Sulfolobus shibatae\"\n"
        f"ORIGIN\n"
        f"        1 {arch_seq.lower()}\n"
        f"//\n"
        f"LOCUS       NC_099003              11624 bp    ds-DNA  linear   VRL 10-NOV-2026\n"
        f"DEFINITION  Enterobacteria phage P4, complete genome.\n"
        f"VERSION     NC_099003.1\n"
        f"SOURCE      Enterobacteria phage P4\n"
        f"  ORGANISM  Enterobacteria phage P4\n"
        f"            Viruses; Duplodnaviria; Heunggongvirae; Uroviricota.\n"
        f"ORIGIN\n"
        f"        1 {sat_seq.lower()}\n"
        f"//\n"
        f"LOCUS       NC_099004               1800 bp    ss-RNA  circular VRL 10-NOV-2026\n"
        f"DEFINITION  Hepatitis delta virus, complete genome.\n"
        f"VERSION     NC_099004.1\n"
        f"SOURCE      Hepatitis delta virus\n"
        f"  ORGANISM  Hepatitis delta virus\n"
        f"            Viruses; Ribozyviria.\n"
        f"ORIGIN\n"
        f"        1 {rna_seq.lower()}\n"
        f"//\n"
    )
    recs = list(mod.parse_genbank_dna_records(io.StringIO(gb)))
    assert [r["accession"] for r in recs] == ["NC_099001.1", "NC_099002.1", "NC_099003.1"]
    assert recs[0]["circular"] is True and recs[0]["seq"] == cress_seq.encode()
    assert recs[1]["circular"] is True
    assert recs[2]["circular"] is False
    assert mod.classify_refseq_viral_subset(recs[0]) == "small_circular_dna"
    assert mod.classify_refseq_viral_subset(recs[1]) == "archaeal_virus"
    assert mod.classify_refseq_viral_subset(recs[2]) == "phage_satellite"


def test_sat1_keeps_public_records_byte_identical_no_5kb_floor_and_enforces_redistributable(tmp_path):
    mod = _load_c55()
    # Public circular < 5 kb CRESS record WITH an exact 60 bp DTR at both ends:
    # under DEC-65 (2), public GenBank/RefSeq records MUST stay byte-identical (NEVER trimmed!)
    dtr = _make_seq(60, seed=21)
    core = _make_seq(2200, seed=22)
    cress_with_dtr = dtr + core + dtr  # 2320 bp (< 5 kb)
    sat_seq = _make_seq(12000, seed=23)

    fa1 = tmp_path / "cress.fna"
    fa1.write_bytes(b">NC_099001.1 Circovirus\n" + cress_with_dtr + b"\n")
    fa2 = tmp_path / "sat.fna"
    fa2.write_bytes(b">AF123456.1 SaPI1\n" + sat_seq + b"\n")

    # Non-redistributable source must be refused in SAT1 (DEC-51, DEC-63)
    bad_sources = [{
        "source_name": "bad_src",
        "source_release": "v1",
        "license": "restricted",
        "redistributable": False,
        "category": "phage_satellite",
        "fasta": fa2,
        "circular_accessions": set(),
    }]
    with pytest.raises(mod.Refused):
        mod.stage_sat1(
            bad_sources,
            tmp_path / "out_bad",
            root=None,
            root_stages=[],
            owner_inputs={},
            retrieved="2026-11-23T00:00:00Z",
        )

    # Also set up an earlier owner stage holding sat_seq to test DEC-26 ownership
    from biobox import ids
    owner_dir = tmp_path / "owner_S1"
    owner_dir.mkdir()
    sat_uid = ids.contig_uid(sat_seq).hex()
    (owner_dir / "genomes.json").write_text(json.dumps([{
        "contig_uids": [sat_uid], "genetic_code": 11, "kind": "virus",
        "assembly_quality": "complete", "circular": [],
    }]))

    good_sources = [
        {
            "source_name": "ncbi_refseq_viral_r238",
            "source_release": "RefSeq-238",
            "license": "INSDC",
            "redistributable": True,
            "category": "small_circular_dna",
            "fasta": fa1,
            "circular_accessions": {"NC_099001.1"},
        },
        {
            "source_name": "ncbi_eutils_c55_delta",
            "source_release": "2026-11-23",
            "license": "INSDC",
            "redistributable": True,
            "category": "phage_satellite",
            "fasta": fa2,
            "circular_accessions": set(),
        },
    ]

    out1 = tmp_path / "SAT1_a"
    out2 = tmp_path / "SAT1_b"
    rep1 = mod.stage_sat1(
        good_sources, out1, root=None, root_stages=[], owner_inputs={"S1": owner_dir},
        retrieved="2026-11-23T00:00:00Z",
    )
    rep2 = mod.stage_sat1(
        good_sources, out2, root=None, root_stages=[], owner_inputs={"S1": owner_dir},
        retrieved="2026-11-23T00:00:00Z",
    )
    assert rep1["genomes"] == 2
    assert rep1["contigs"] == 2
    assert rep1["lt_5kb_contigs"] == 1
    assert rep1["held_contigs"] == 1
    assert rep1["owner_contigs"] == {"S1": 1}
    assert rep1["new_contigs"] == 1
    # Byte-identical public record: contig_uid must equal untrimmed cress_with_dtr
    cress_uid = ids.contig_uid(cress_with_dtr).hex()
    genomes = json.loads((out1 / "genomes.json").read_text())
    assert genomes[0]["contig_uids"] == [cress_uid]
    assert genomes[0]["circular"] == [cress_uid]
    # Determinism across reruns
    for name in ("genomes.json", "records.tsv", "sat1_metadata.tsv", "shared_records.tsv"):
        assert (out1 / name).read_bytes() == (out2 / name).read_bytes()


def test_v3_vire_votu_selection_dec65_dtr_trimming_and_determinism(tmp_path):
    mod = _load_c55()
    # Construct synthetic VIRE metadata + tar with:
    # - species_cluster 1: two DNA members -> higher completeness rep wins, other goes to left_for_later;
    #   rep has DTR with exact 50 bp repeat -> trimmed and staged circular (DEC-65)
    # - species_cluster 2: DNA member flagged DTR but NO repeat -> staged circular untrimmed
    # - species_cluster 3: DNA member flagged DTR with partial/1-mismatch repeat -> staged linear untrimmed
    # - species_cluster 4: Riboviria member -> deferred to left_for_later (reason=riboviria_deferred_c50)
    dtr50 = _make_seq(50, seed=31)
    core1 = _make_seq(4000, seed=32)
    seq1 = dtr50 + core1 + dtr50  # exact 50 bp DTR
    seq1_other = _make_seq(3500, seed=33)
    seq2 = _make_seq(3200, seed=34)  # DTR flag, no repeat
    bad_dtr = bytearray(dtr50)
    bad_dtr[25] = ord("A") if bad_dtr[25] != ord("A") else ord("C")
    seq3 = dtr50 + _make_seq(3800, seed=35) + bytes(bad_dtr)  # partial DTR -> linear untrimmed
    seq4 = _make_seq(2500, seed=36)  # Riboviria

    meta = tmp_path / "vire_meta.txt"
    cols = [
        "genome_id", "biosample_name", "study_id", "contig_length",
        "checkv_completeness", "checkv_contamination", "checkv_checkv_quality",
        "genomad_topology", "species_cluster", "alignment_to_imgvr", "realm",
    ]
    rows = [
        ["vire_genome_1", "SAMN1", "1", str(len(seq1)), "100.0", "0.0", "Complete", "DTR", "1", "None", "Duplodnaviria"],
        ["vire_genome_1b", "SAMN1", "1", str(len(seq1_other)), "80.0", "0.0", "Medium-quality", "No terminal repeats", "1", "None", "Duplodnaviria"],
        ["vire_genome_2", "SAMN1", "1", str(len(seq2)), "95.0", "0.0", "High-quality", "DTR", "2", "True", "Monodnaviria"],
        ["vire_genome_3", "SAMN2", "2", str(len(seq3)), "92.0", "0.0", "High-quality", "DTR", "3", "None", "Duplodnaviria"],
        ["vire_genome_4", "SAMN2", "2", str(len(seq4)), "100.0", "0.0", "Complete", "DTR", "4", "None", "Riboviria"],
    ]
    with open(meta, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    tar_path = tmp_path / "vire.tar"
    with tarfile.open(tar_path, "w") as tf:
        for bname, recs in [
            ("genomes/vire_biosample_SAMN1.fna.gz", [
                (b"vire_genome_1", seq1), (b"vire_genome_1b", seq1_other), (b"vire_genome_2", seq2),
            ]),
            ("genomes/vire_biosample_SAMN2.fna.gz", [
                (b"vire_genome_3", seq3), (b"vire_genome_4", seq4),
            ]),
        ]:
            buf = io.BytesIO()
            with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:
                for gid, s in recs:
                    gz.write(b">" + gid + b"\n" + s + b"\n")
            data = buf.getvalue()
            ti = tarfile.TarInfo(name=bname)
            ti.size = len(data)
            tf.addfile(ti, io.BytesIO(data))

    out1_m, out1_f = tmp_path / "v3_m1", tmp_path / "v3_1.fna"
    out2_m, out2_f = tmp_path / "v3_m2", tmp_path / "v3_2.fna"
    rep1 = mod.stage_v3_from_tar(
        meta, tar_path, out1_m, out1_f, root=None, root_stages=[], owner_inputs={},
        retrieved="2026-11-23T00:00:00Z",
    )
    rep2 = mod.stage_v3_from_tar(
        meta, tar_path, out2_m, out2_f, root=None, root_stages=[], owner_inputs={},
        retrieved="2026-11-23T00:00:00Z",
    )
    assert rep1["genomes"] == 3
    assert rep1["trims"] == 1
    assert rep1["trimmed_bases"] == 50
    assert rep1["circular_genomes"] == 2  # vire_genome_1 (exact trimmed) and vire_genome_2 (none -> circular)
    assert rep1["left_for_later_rows"] == 2  # vire_genome_1b + vire_genome_4 (Riboviria)
    assert hashlib.sha256(out1_f.read_bytes()).hexdigest() == hashlib.sha256(out2_f.read_bytes()).hexdigest()
    assert (out1_m / "genomes.json").read_bytes() == (out2_m / "genomes.json").read_bytes()
    assert (out1_m / "trims.tsv").read_bytes() == (out2_m / "trims.tsv").read_bytes()
