"""Card C-56 (DEC-69, DEC-70, DEC-71): long-read (PacBio HiFi / ONT) metagenome assembly crosswalk,
biosample/run deduplication, ART-RT habitat stratification, DEC-57 >= 5,000 bp filtering,
DEC-65 circular/DTR trimming, and CATALOG_GOAL.md §4.1 remainder ordering.

Oracles are independent of the script: direct blake3 contig_uid hashing and brute-force Python
deduplication and habitat-priority ranking over shuffled input rows.
"""

from __future__ import annotations

import gzip
import importlib.util
import json
import random
import sys
from pathlib import Path

import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("C-56")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts" / "stress" / "build_stage_pod"


def _load():
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(
        "c56_longread_metagenomes", SCRIPTS / "c56_longread_metagenomes.py"
    )
    if spec is None or spec.loader is None:
        raise FileNotFoundError(SCRIPTS / "c56_longread_metagenomes.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c56_longread_metagenomes"] = mod
    spec.loader.exec_module(mod)
    return mod


def cuid(seq: str) -> str:
    return blake3(seq.upper().encode()).digest(length=16).hex()


def seq_of(seed: int, n: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_platform_classification():
    mod = _load()
    assert mod.classify_platform(["PACBIO_SMRT"], ["Revio"]) == ("pacbio", True)
    assert mod.classify_platform(["PACBIO_SMRT"], ["Sequel IIe"]) == ("pacbio", True)
    assert mod.classify_platform(["PACBIO_SMRT"], ["PacBio RS II"]) == ("pacbio", False)
    assert mod.classify_platform(["OXFORD_NANOPORE"], ["PromethION"]) == ("ont", False)
    assert mod.classify_platform(["PACBIO_SMRT", "ILLUMINA"], ["Sequel II", "NovaSeq 6000"]) == (
        "hybrid_pacbio",
        True,
    )
    assert mod.classify_platform(["OXFORD_NANOPORE", "ILLUMINA"], ["MinION", "NextSeq"]) == (
        "hybrid_ont",
        False,
    )
    assert mod.classify_platform(["PACBIO_SMRT", "OXFORD_NANOPORE"], ["Revio", "PromethION"]) == (
        "pacbio_ont",
        True,
    )
    assert mod.classify_platform(["ILLUMINA"], ["Illumina NovaSeq 6000"]) == ("short_only", False)
    # IMG/M comma-separated Sequencing Method strings
    assert mod.classify_imgm_method("PacBio Revio") == ("pacbio", True)
    assert mod.classify_imgm_method("PacBio Revio, Illumina NovaSeq S4") == ("hybrid_pacbio", True)
    assert mod.classify_imgm_method("Oxford Nanopore MinION") == ("ont", False)
    assert mod.classify_imgm_method("Illumina HiSeq 2500") == ("short_only", False)


def test_habitat_stratification_art_rt_priority():
    mod = _load()
    # DEC-71: soil/terrestrial, wastewater/sludge, freshwater/cyanobacterial lead
    s_soil, p_soil = mod.classify_habitat("soil metagenome")
    s_sludge, p_sludge = mod.classify_habitat("activated sludge metagenome")
    s_fresh, p_fresh = mod.classify_habitat("freshwater metagenome")
    s_cyano, p_cyano = mod.classify_habitat("microbial mat metagenome", "Environmental;Aquatic;Freshwater;Microbial mat")
    s_marine, p_marine = mod.classify_habitat("marine metagenome")
    s_gut, p_gut = mod.classify_habitat("human gut metagenome")
    s_other, p_other = mod.classify_habitat("metagenome")

    assert s_soil == "soil_terrestrial"
    assert s_sludge == "wastewater_sludge"
    assert s_fresh == "freshwater_aquatic"
    assert s_cyano == "freshwater_aquatic"
    assert s_marine == "marine_sediment"
    assert s_gut == "host_associated"
    assert s_other == "other_unlabelled"
    assert p_soil < p_marine < p_gut < p_other
    assert p_sludge < p_marine
    assert p_fresh < p_marine


def test_crosswalk_dedup_and_order_invariance():
    mod = _load()
    # Two assemblies of the same biosample SAMN001 (one HiFi ERZ, one smaller ONT/hybrid GCA):
    # winner must be the larger/HiFi primary assembly regardless of input row order.
    # Also a MAG row must be excluded, and a JGI IMG/M row on the same biosample must lose to open ENA/NCBI.
    candidates = [
        {
            "accession": "ERZ100001",
            "source": "ena_analysis",
            "access": "open",
            "assembly_type": "primary metagenome",
            "biosample": "SAMN001",
            "runs": ["SRR001"],
            "study": "PRJNA1",
            "platform_class": "pacbio",
            "is_hifi": True,
            "est_bp": 500_000_000,
            "ftp_url": "ftp.sra.ebi.ac.uk/vol1/ERZ100/ERZ100001/contigs.fa.gz",
            "scientific_name": "soil metagenome",
            "ecosystem": "",
        },
        {
            "accession": "GCA_000000015.1",
            "source": "ncbi_wgs",
            "access": "open",
            "assembly_type": "primary metagenome",
            "biosample": "SAMN001",
            "runs": ["SRR001"],
            "study": "PRJNA1",
            "platform_class": "hybrid_pacbio",
            "is_hifi": True,
            "est_bp": 400_000_000,
            "ftp_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/000/015/GCA_000000015.1_asm/GCA_000000015.1_asm_genomic.fna.gz",
            "scientific_name": "soil metagenome",
            "ecosystem": "",
        },
        {
            "accession": "3300099999",
            "source": "img_m",
            "access": "jgi_login",
            "assembly_type": "primary metagenome",
            "biosample": "SAMN001",
            "runs": ["SRR001"],
            "study": "Gs0001",
            "platform_class": "pacbio",
            "is_hifi": True,
            "est_bp": 600_000_000,
            "ftp_url": "",
            "scientific_name": "soil metagenome",
            "ecosystem": "Environmental;Terrestrial;Soil",
        },
        {
            "accession": "ERZ200002",
            "source": "ena_analysis",
            "access": "open",
            "assembly_type": "primary metagenome",
            "biosample": "SAMN002",
            "runs": ["SRR002"],
            "study": "PRJEB2",
            "platform_class": "ont",
            "is_hifi": False,
            "est_bp": 300_000_000,
            "ftp_url": "ftp.sra.ebi.ac.uk/vol1/ERZ200/ERZ200002/contigs.fa.gz",
            "scientific_name": "activated sludge metagenome",
            "ecosystem": "",
        },
        {
            "accession": "3300088888",
            "source": "img_m",
            "access": "jgi_login",
            "assembly_type": "primary metagenome",
            "biosample": "SAMN003",
            "runs": ["SRR003"],
            "study": "Gs0003",
            "platform_class": "pacbio",
            "is_hifi": True,
            "est_bp": 900_000_000,
            "ftp_url": "",
            "scientific_name": "freshwater metagenome",
            "ecosystem": "Environmental;Aquatic;Freshwater",
        },
    ]

    w1, d1 = mod.deduplicate_assemblies(candidates)
    shuffled = list(reversed(candidates))
    w2, d2 = mod.deduplicate_assemblies(shuffled)

    assert [x["accession"] for x in w1] == [x["accession"] for x in w2]
    assert [x["accession"] for x in d1] == [x["accession"] for x in d2]

    win_accs = {x["accession"] for x in w1}
    assert win_accs == {"ERZ100001", "ERZ200002", "3300088888"}
    dup_accs = {x["accession"]: x["winner_accession"] for x in d1}
    assert dup_accs == {"GCA_000000015.1": "ERZ100001", "3300099999": "ERZ100001"}


def test_tranche_selection_and_remainder_ordering():
    mod = _load()
    # 6 open assemblies across habitats: soil, sludge, freshwater must be prioritized before gut
    rows = [
        {"accession": "A_GUT1", "access": "open", "habitat": "host_associated", "habitat_priority": 5, "study": "S_G1", "est_bp": 100_000_000, "platform_class": "pacbio", "is_hifi": True},
        {"accession": "A_SOIL1", "access": "open", "habitat": "soil_terrestrial", "habitat_priority": 1, "study": "S_S1", "est_bp": 100_000_000, "platform_class": "pacbio", "is_hifi": True},
        {"accession": "A_SLUDGE1", "access": "open", "habitat": "wastewater_sludge", "habitat_priority": 2, "study": "S_W1", "est_bp": 100_000_000, "platform_class": "pacbio", "is_hifi": True},
        {"accession": "A_FRESH1", "access": "open", "habitat": "freshwater_aquatic", "habitat_priority": 3, "study": "S_F1", "est_bp": 100_000_000, "platform_class": "ont", "is_hifi": False},
        {"accession": "A_MARINE1", "access": "open", "habitat": "marine_sediment", "habitat_priority": 4, "study": "S_M1", "est_bp": 100_000_000, "platform_class": "pacbio", "is_hifi": True},
        {"accession": "A_SOIL2", "access": "open", "habitat": "soil_terrestrial", "habitat_priority": 1, "study": "S_S2", "est_bp": 100_000_000, "platform_class": "pacbio", "is_hifi": True},
    ]
    sel, rem = mod.select_tranche(rows, max_bases=350_000_000, seed=20261123)
    # Cap of 350 Mbp fits 3 assemblies; all 3 must come from ART-RT-priority habitats before gut
    assert len(sel) == 3
    assert len(rem) == 3
    assert "A_GUT1" not in {r["accession"] for r in sel}
    assert {r["accession"] for r in sel} | {r["accession"] for r in rem} == {r["accession"] for r in rows}
    # Remainder carries contiguous lr_rank starting after selection
    assert [r["lr_rank"] for r in sel] == [1, 2, 3]
    assert [r["lr_rank"] for r in rem] == [4, 5, 6]
    # Deterministic under row permutation
    sel2, rem2 = mod.select_tranche(list(reversed(rows)), max_bases=350_000_000, seed=20261123)
    assert [r["accession"] for r in sel] == [r["accession"] for r in sel2]
    assert [r["accession"] for r in rem] == [r["accession"] for r in rem2]

    # Submitter assemblies must be prioritized ahead of EMG re-assemblies
    rows_emg = rows + [
        {"accession": "A_EMG_SOIL", "access": "open", "center_name": "EMG", "habitat": "soil_terrestrial", "habitat_priority": 1, "study": "S_EMG", "est_bp": 50_000_000, "platform_class": "pacbio", "is_hifi": True},
    ]
    sel_e, rem_e = mod.select_tranche(rows_emg, max_bases=600_000_000, seed=20261123)
    assert "A_EMG_SOIL" not in {r["accession"] for r in sel_e}
    assert rem_e[-1]["accession"] == "A_EMG_SOIL"


def test_dec57_filter_and_dec65_circular_dtr_staging(tmp_path: Path):
    mod = _load()

    # Construct an assembly .fa.gz with 6 contigs:
    # 1. short (< 5,000 bp) -> dropped under DEC-57
    # 2. linear >= 5,000 bp -> kept linear untrimmed
    # 3. circular header + exact 60 bp DTR -> trimmed by 60 bp, staged circular, recorded in trims.tsv
    # 4. circular header + no DTR -> untrimmed, staged circular
    # 5. circular header + partial DTR (1 mismatch in 60 bp) -> untrimmed, staged linear
    # 6. circular header + low-complexity DTR (60x 'A') -> untrimmed, staged linear
    s_short = seq_of(1, 4_999)
    s_lin = seq_of(2, 6_000)
    dtr_exact = seq_of(3, 60)
    s_circ_exact_core = dtr_exact + seq_of(4, 5_500)
    s_circ_exact = s_circ_exact_core + dtr_exact
    s_circ_none = seq_of(5, 5_600)
    dtr_bad = dtr_exact[:30] + ("A" if dtr_exact[30] != "A" else "C") + dtr_exact[31:]
    s_circ_partial = dtr_exact + seq_of(6, 5_700) + dtr_bad
    s_circ_lowc = ("A" * 60) + seq_of(7, 5_800) + ("A" * 60)

    asm_dir = tmp_path / "dl"
    asm_dir.mkdir()
    fa_path = asm_dir / "ERZ999001.fa.gz"
    with gzip.open(fa_path, "wt") as fh:
        fh.write(f">ctg_short len=4999\n{s_short}\n")
        fh.write(f">ctg_lin len=6000\n{s_lin}\n")
        fh.write(f">ctg_circ_exact circular=true\n{s_circ_exact}\n")
        fh.write(f">ctg_circ_none topology=circular\n{s_circ_none}\n")
        fh.write(f">ctg_circ_partial circ=yes\n{s_circ_partial}\n")
        fh.write(f">ctg_circ_lowc circular=Y\n{s_circ_lowc}\n")

    # Also test owner overlap: mark s_lin as owned by SPX1
    owner_dir = tmp_path / "SPX1"
    owner_dir.mkdir()
    owner_uid = cuid(s_lin)
    (owner_dir / "genomes.json").write_text(
        json.dumps([{"contig_uids": [owner_uid], "genetic_code": 11, "kind": "prokaryote", "circular": [], "assembly_quality": "contig"}])
    )

    manifest = [
        {
            "accession": "ERZ999001",
            "path": str(fa_path),
            "source": "ena_analysis",
            "habitat": "soil_terrestrial",
            "platform_class": "pacbio",
            "retrieved_at": "2026-11-23T00:00:00Z",
        }
    ]
    out_fasta_dir = tmp_path / "staging_fasta"
    out_meta_dir = tmp_path / "staging_meta"
    rep = mod.stage_assemblies(
        manifest,
        out_fasta_dir=out_fasta_dir,
        out_meta_dir=out_meta_dir,
        final_fasta_dir=out_fasta_dir,
        min_contig_len=5_000,
        owner_inputs=[("SPX1", owner_dir)],
    )

    assert rep["records_read"] == 6
    assert rep["records_below_min_len"] == 1
    assert rep["bases_below_min_len"] == 4_999
    assert rep["distinct_contigs"] == 5
    assert rep["circular_contigs"] == 2
    assert rep["dtr_exact_trimmed"] == 1
    assert rep["dtr_trimmed_bases"] == 60
    assert rep["held_by_earlier"]["SPX1"] == 1
    assert rep["net_new_contigs"] == 4

    genomes = json.loads((out_meta_dir / "genomes.json").read_text())
    assert len(genomes) == 5
    by_uid = {g["contig_uids"][0]: g for g in genomes}
    assert cuid(s_circ_exact_core) in by_uid
    assert by_uid[cuid(s_circ_exact_core)]["circular"] == [cuid(s_circ_exact_core)]
    assert by_uid[cuid(s_circ_none)]["circular"] == [cuid(s_circ_none)]
    assert by_uid[cuid(s_circ_partial)]["circular"] == []
    assert by_uid[cuid(s_circ_lowc)]["circular"] == []
    assert by_uid[cuid(s_lin)]["circular"] == []
    assert all(g["kind"] == "prokaryote" and g["genetic_code"] == 11 and g["assembly_quality"] == "contig" for g in genomes)

    trims_lines = (out_meta_dir / "trims.tsv").read_text().strip().splitlines()
    assert len(trims_lines) == 2  # header + 1 exact trim
    assert "60" in trims_lines[1]
