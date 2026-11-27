"""Card C-85 (DEC-88 (1), (2), DEC-80, AGENTS §8): gentle next-wave `tower` fetch queues and GlobDB stream slicer.

Every test builds small constructed fixtures with hand-verified expected outputs (no snapshots). Each test runs under
both the normal PATH and a tools-hidden PATH (`/usr/bin:/bin`) so CI without pod binaries exercises the full logic.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import io
import json
import os
import subprocess
import sys
import tarfile
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-85")

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "orchestrator" / "c85_prepare_tower_queues.py"
CHAIN_SCRIPT = ROOT / "scripts" / "orchestrator" / "after_retry20_gentle_chain.sh"
FETCH_SCRIPT = ROOT / "scripts" / "orchestrator" / "fetch_to_tower.sh"


def _load():
    if not SCRIPT.is_file():
        raise FileNotFoundError(f"missing script under test: {SCRIPT}")
    spec = importlib.util.spec_from_file_location("c85_prepare_tower_queues", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c85_prepare_tower_queues"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(params=["normal-path", "tools-hidden"], autouse=True)
def path_mode(request, monkeypatch):
    if request.param == "tools-hidden":
        monkeypatch.setenv("PATH", "/usr/bin:/bin")
    return request.param


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _tar_gz_bytes(members: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:
        with tarfile.open(fileobj=gz, mode="w") as tf:
            for name, payload in members.items():
                info = tarfile.TarInfo(name=name)
                info.size = len(payload)
                info.mtime = 0
                tf.addfile(info, io.BytesIO(payload))
    return buf.getvalue()


def test_ncbi_assembly_url_builder_and_five_column_tsv_validation(tmp_path: Path):
    mod = _load()

    # From full NCBI ftp_path (both ftp:// and https:// normalized to https://)
    fname, url = mod.build_ncbi_fna_url(
        "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/053/992/195/GCA_053992195.1_ASM5399219v1/"
    )
    assert fname == "GCA_053992195.1_genomic.fna.gz"
    assert url == (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/053/992/195/"
        "GCA_053992195.1_ASM5399219v1/GCA_053992195.1_ASM5399219v1_genomic.fna.gz"
    )

    # From versioned GCF accession + asm_name with spaces
    fname2, url2 = mod.build_ncbi_fna_url("GCF_000005845.2", asm_name="ASM584 v2")
    assert fname2 == "GCF_000005845.2_genomic.fna.gz"
    assert url2 == (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/"
        "GCF_000005845.2_ASM584_v2/GCF_000005845.2_ASM584_v2_genomic.fna.gz"
    )

    # Invalid inputs refused
    for bad in ("na", "", "https://evil.example.org/GCA_000000015.1", "GCA_12345.1"):
        with pytest.raises(mod.Refused):
            mod.build_ncbi_fna_url(bad)

    # Five-column job validation & deterministic writer
    good_rows = [
        ("DSV70_NCBI/assemblies", fname2, url2, "0123456789abcdef0123456789abcdef", "INSDC"),
        ("DSV70_NCBI/assemblies", fname, url, "-", "INSDC"),
    ]
    out_tsv = tmp_path / "jobs.tsv"
    mod.write_jobs_tsv(out_tsv, ["# test header"], good_rows)
    lines = [ln for ln in out_tsv.read_text().splitlines() if ln and not ln.startswith("#")]
    assert len(lines) == 2
    assert [ln.split("\t")[1] for ln in lines] == sorted([fname, fname2])

    # Duplicate destination or shell metacharacters refused
    with pytest.raises(mod.Refused):
        mod.write_jobs_tsv(tmp_path / "dup.tsv", [], [good_rows[0], good_rows[0]])
    for bad_row in [
        ("../bad", fname, url, "-", "INSDC"),
        ("DSV70", "SHA256SUMS", url, "-", "INSDC"),
        ("DSV70", "x.partial", url, "-", "INSDC"),
        ("DSV70", fname, "http://ftp.ncbi.nlm.nih.gov/x.fna.gz", "-", "INSDC"),
        ("DSV70", fname, url, "NOT_MD5", "INSDC"),
        ("DSV70", fname, url, "-", "bad;lic"),
    ]:
        with pytest.raises(mod.Refused):
            mod.validate_job_row(bad_row)


def test_mgnify_biome_and_overlap_filter():
    mod = _load()

    head = (
        "Genome\tGenome_type\tLength\tN_contigs\tN50\tGC_content\tCompleteness\tContamination\t"
        "rRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tGenome_accession\tSpecies_rep\tLineage\t"
        "Sample_accession\tStudy_accession\tCountry\tContinent\tFTP_download\n"
    )
    rows_marine = (
        head
        # Non-rep row: must be skipped because Genome != Species_rep
        + "MGYG000296000\tMAG\t1700000\t30\t80000\t36.0\t90.0\t2.0\t0\t0\t0\t18\tCAJVWI01\tMGYG000444050\t"
        "d__Bacteria;s__A\tSAMEA1\tPRJEB1\tOcean\tNA\t"
        "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/marine/v3.0/all_genomes/MGYG0004440/MGYG000444050/genomes1/MGYG000296000.gff.gz\n"
        # Rep already in GlobDB by MGYG ID: must be excluded
        + "MGYG000444050\tMAG\t2100000\t20\t120000\t40.0\t95.0\t1.0\t0\t0\t0\t20\tCAJVWJ01\tMGYG000444050\t"
        "d__Bacteria;s__A\tSAMEA2\tPRJEB1\tOcean\tNA\t"
        "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/marine/v3.0/all_genomes/MGYG0004440/MGYG000444050/genomes1/MGYG000444050.gff.gz\n"
        # Rep already in GTDB R232 by WGS master prefix: must be excluded
        + "MGYG000444051\tMAG\t2200000\t15\t150000\t42.0\t96.0\t0.5\t0\t0\t0\t21\tCAJVWK01\tMGYG000444051\t"
        "d__Bacteria;s__B\tSAMEA3\tPRJEB1\tOcean\tNA\t"
        "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/marine/v3.0/all_genomes/MGYG0004440/MGYG000444051/genomes1/MGYG000444051.gff.gz\n"
        # Novel non-human marine rep: MUST BE KEPT
        + "MGYG000444052\tMAG\t2500000\t12\t200000\t44.0\t97.0\t0.8\t0\t0\t0\t22\tCAJVWL01\tMGYG000444052\t"
        "d__Archaea;s__Novel\tSAMEA4\tPRJEB1\tOcean\tNA\t"
        "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/marine/v3.0/all_genomes/MGYG0004440/MGYG000444052/genomes1/MGYG000444052.gff.gz\n"
    )
    rows_human = (
        head
        + "MGYG000000001\tMAG\t2400000\t10\t250000\t48.0\t99.0\t0.1\t0\t0\t0\t22\tGCA_000000001\tMGYG000000001\t"
        "d__Bacteria;s__Hum\tSAMEA9\tPRJEB9\tUSA\tNA\t"
        "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v3.0/all_genomes/MGYG0000000/MGYG000000001/genomes1/MGYG000000001.gff.gz\n"
    )

    kept, stats = mod.filter_mgnify_catalogs(
        catalog_texts={
            "marine_v3.0_genomes-all_metadata.tsv": rows_marine,
            "human-gut_v3.0_genomes-all_metadata.tsv": rows_human,
            "human-oral_v1.0_genomes-all_metadata.tsv": rows_human,
            "marine-eukaryotes_v1.0_genomes-all_metadata.tsv": rows_human,
        },
        known_mgyg_ids={"MGYG000444050"},
        known_insdc_accessions={"GCA_999999995"},
        known_wgs_prefixes={"CAJVWK01"},
        known_biosamples=set(),
    )
    assert [r["genome"] for r in kept] == ["MGYG000444052"]
    assert kept[0]["https_url"].startswith("https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/marine/v3.0/")
    assert stats["excluded_biomes"] == 3
    assert stats["kept_reps"] == 1


def test_gtdb_r232_unique_new_species_delta_filter():
    mod = _load()

    gtdb_r232_reps = [
        # Existing R226 species: excluded
        {"accession": "GCF_000005845", "accession_versioned": "RS_GCF_000005845.2", "genbank_accession": "GCA_000005845",
         "species": "Escherichia coli", "domain": "Bacteria", "is_representative": True},
        # New R232 species, but covered by GlobDB overlap set: excluded
        {"accession": "GCA_011111115", "accession_versioned": "GB_GCA_011111115.1", "genbank_accession": "GCA_011111115",
         "species": "Sp_globdb_covered", "domain": "Bacteria", "is_representative": True},
        # New R232 species, already in built stage SYM1: excluded
        {"accession": "GCA_022222225", "accession_versioned": "GB_GCA_022222225.1", "genbank_accession": "GCA_022222225",
         "species": "Sp_built_sym1", "domain": "Bacteria", "is_representative": True},
        # New R232 species, unique to GTDB R232: KEPT
        {"accession": "GCA_033333335", "accession_versioned": "GB_GCA_033333335.1", "genbank_accession": "GCA_033333335",
         "species": "Sp_gtdb_unique_new", "domain": "Archaea", "is_representative": True},
    ]
    post_r226_candidates = [
        # Absorbed into GTDB R232 non-rep or rep: excluded
        {"accession": "GCA_033333335.1", "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/033/333/335/GCA_033333335.1_ASM1"},
        {"accession": "GCA_044444445.1", "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/044/444/445/GCA_044444445.1_ASM2"},
        # True post-R232 species rep absent from GTDB R232 and GlobDB r232: KEPT
        {"accession": "GCF_055555555.1", "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/055/555/555/GCF_055555555.1_ASM3"},
    ]
    res = mod.select_gtdb_r232_delta(
        gtdb_r232_reps=gtdb_r232_reps,
        all_gtdb_r232_accessions={"GCF_000005845", "GCA_000005845", "GCA_011111115", "GCA_022222225", "GCA_033333335", "GCA_044444445"},
        r226_species={"Escherichia coli"},
        r226_rep_accessions={"GCF_000005845"},
        built_accessions={"GCF_000005845", "GCA_022222225"},
        globdb_covered_accessions={"GCA_011111115"},
        post_r226_candidates=post_r226_candidates,
    )
    assert [r["accession"] for r in res["gtdb_unique_new_reps"]] == ["GCA_033333335"]
    assert [r["accession"] for r in res["post_r232_ncbi_delta"]] == ["GCF_055555555"]


def test_extract_globdb_subset_streams_tar_gz_and_retains_only_targets(tmp_path: Path):
    mod = _load()

    g1 = b">c1\nACGTACGT\n"
    g2 = b">c2\nTTTTGGGG\n"
    g3 = b">c3\nCCCCAAAA\n"
    archive = tmp_path / "globdb_r232_genome_fasta.tar.gz"
    archive.write_bytes(
        _tar_gz_bytes(
            {
                "globdb_r232_genome_fasta/GCMETA_00001.fa": g1,
                "globdb_r232_genome_fasta/GCA_000005845.fa": g2,
                "globdb_r232_genome_fasta/NGDC_00002.fa.gz": gzip.compress(g3, mtime=0),
            }
        )
    )
    targets = tmp_path / "targets.txt"
    targets.write_text("GCMETA_00001\nNGDC_00002\nMISSING_99999\n")
    out_dir = tmp_path / "extracted"

    # Exercise CLI subcommand (works under both normal and tools-hidden PATH)
    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "extract-globdb-subset",
            "--tar-gz",
            str(archive),
            "--ids",
            str(targets),
            "--out-dir",
            str(out_dir),
        ],
        capture_output=True,
        text=True,
        env=dict(os.environ),
        timeout=60,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    summary = json.loads(proc.stdout)
    assert summary["matched"] == 2
    assert summary["missing"] == ["MISSING_99999"]
    assert (out_dir / "GCMETA_00001.fa").read_bytes() == g1
    assert not (out_dir / "GCA_000005845.fa").exists()
    assert gzip.decompress((out_dir / "NGDC_00002.fa.gz").read_bytes()) == g3

    manifest_lines = (out_dir / "extracted_manifest.tsv").read_text().splitlines()
    assert manifest_lines[0].split("\t") == ["genome_id", "filename", "bytes", "sha256"]
    assert len(manifest_lines) == 3
    assert _sha256(g1) in manifest_lines[1]


def test_after_retry20_gentle_chain_script_syntax_and_flags(tmp_path: Path):
    assert CHAIN_SCRIPT.is_file(), f"missing {CHAIN_SCRIPT}"

    # 1. bash -n syntax check (explicit /bin/bash so tools-hidden PATH works)
    r_syn = subprocess.run(["/bin/bash", "-n", str(CHAIN_SCRIPT)], capture_output=True, text=True, timeout=30)
    assert r_syn.returncode == 0, r_syn.stdout + r_syn.stderr

    # 2. Static policy inspection: ionice -c 3, nice -n 19, --jobs 2 --pause 1 and --jobs 1, STOP check
    text = CHAIN_SCRIPT.read_text()
    assert "ionice -c 3" in text
    assert "nice -n 19" in text
    assert "STOP" in text

    # 3. Behavioral dry-run / stub test verifying STOP halts before fetching and flags never exceed --jobs 2
    dest = tmp_path / "incoming"
    dest.mkdir()
    lists_dir = tmp_path / "lists"
    lists_dir.mkdir()
    for name in (
        "23_r26_11_6_asg1_fna.tsv",
        "24_r26_11_6_dsv70_prjna1244896.tsv",
        "25_r26_11_6_mgnify_nonhuman_reps.tsv",
        "26_r26_11_6_globdb_r232.tsv",
        "27_r26_11_6_gtdb_r232_delta.tsv",
    ):
        (lists_dir / name).write_text("# stub\n")
    (dest / "STOP").write_text("")
    dead = subprocess.Popen(["/bin/sh", "-c", "exit 0"])
    dead.wait()

    r_stop = subprocess.run(
        ["/bin/bash", str(CHAIN_SCRIPT), "--dest", str(dest), "--lists-dir", str(lists_dir), "--dry-run", str(dead.pid)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert r_stop.returncode == 0, r_stop.stdout + r_stop.stderr
    assert "STOP" in r_stop.stdout

    (dest / "STOP").unlink()
    r_dry = subprocess.run(
        ["/bin/bash", str(CHAIN_SCRIPT), "--dest", str(dest), "--lists-dir", str(lists_dir), "--dry-run", str(dead.pid)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert r_dry.returncode == 0, r_dry.stdout + r_dry.stderr
    dry_lines = [ln for ln in r_dry.stdout.splitlines() if "fetch_to_tower.sh" in ln]
    assert len(dry_lines) == 5, f"expected 5 dry-run commands in stdout, got:\n{r_dry.stdout}"
    for ln in dry_lines:
        assert "ionice -c 3" in ln and "nice -n 19" in ln, ln
        if "26_r26_11_6_globdb_r232.tsv" in ln:
            assert "--jobs 1" in ln, ln
        else:
            assert "--jobs 2" in ln and "--pause 1" in ln, ln


def test_asg1_study_bin_mapper_handles_clean_suffix_and_prefers_latest():
    mod = _load()

    study = ["HeGe__D4993_C5_H4_Bin_142", "Hela__Meg22_1012_Bin_214", "HeHo__Unmapped_Bin_99"]
    ncbi_rows = [
        {
            "acc_ver": "GCA_054387275.1",
            "isolate": "D4993_C5_H4_Bin_142_clean",
            "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/054/387/275/GCA_054387275.1_ASM5438727v1/",
            "seq_rel_date": "2026-01-06",
            "submitter": "University of Texas at Austin",
            "version_status": "latest",
        },
        {
            "acc_ver": "GCA_054384215.1",
            "isolate": "Meg22_1012_Bin_214",
            "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/054/384/215/GCA_054384215.1_ASM5438421v1/",
            "seq_rel_date": "2026-01-06",
            "submitter": "University of Texas at Austin",
            "version_status": "latest",
        },
        {
            "acc_ver": "GCA_057965895.1",
            "isolate": "Meg22_1012_Bin_214",
            "ftp_path": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/057/965/895/GCA_057965895.1_ASM5796589v1/",
            "seq_rel_date": "2026-06-04",
            "submitter": "University of Texas at Austin",
            "version_status": "latest",
        },
    ]
    chosen, unmapped = mod.map_asg1_study_bins(study, ncbi_rows)
    assert [r["acc_ver"] for r in chosen] == ["GCA_054387275.1", "GCA_057965895.1"]
    assert unmapped == ["HeHo__Unmapped_Bin_99"]


def test_committed_c85_job_lists_pass_fetch_to_tower_check():
    sources = ROOT.joinpath("sources.yaml")
    if not sources.is_file():
        pytest.skip("sources.yaml is not shipped inside the bare runtime image")
    lists = [
        ROOT / "scripts" / "orchestrator" / "tower_jobs" / name
        for name in (
            "23_r26_11_6_asg1_fna.tsv",
            "24_r26_11_6_dsv70_prjna1244896.tsv",
            "25_r26_11_6_mgnify_nonhuman_reps.tsv",
            "26_r26_11_6_globdb_r232.tsv",
            "27_r26_11_6_gtdb_r232_delta.tsv",
        )
    ]
    for p in lists:
        assert p.is_file(), f"missing committed job list {p}"
    r = subprocess.run(
        ["/bin/bash", str(FETCH_SCRIPT), "--check", "--sources", str(sources), *map(str, lists)],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert r.returncode == 0, r.stdout + r.stderr
    assert "check: 41134 URLs, all covered by sources.yaml; nothing fetched" in r.stdout

