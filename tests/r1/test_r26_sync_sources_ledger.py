"""Card R-26 (DEC-50, DEC-59, DEC-80, DEC-88, DEC-89, AGENTS §6, §8): ledger sync & verifier tests.

Every test builds synthetic fixtures in `tmp_path` with hand-verified expected outputs (no snapshots and no references
to top-level repo files outside `scripts/`). Each test runs under both the normal PATH and a tools-hidden PATH
(`/usr/bin:/bin`).
"""

from __future__ import annotations

import hashlib
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

pytestmark = pytest.mark.covers("R-26")

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "orchestrator" / "r26_sync_sources_ledger.py"
FETCH_SCRIPT = ROOT / "scripts" / "orchestrator" / "fetch_to_tower.sh"


def _load():
    if not SCRIPT.is_file():
        raise FileNotFoundError(f"missing script under test: {SCRIPT}")
    spec = importlib.util.spec_from_file_location("r26_sync_sources_ledger", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["r26_sync_sources_ledger"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(params=["normal-path", "tools-hidden"], autouse=True)
def path_mode(request, monkeypatch):
    if request.param == "tools-hidden":
        monkeypatch.setenv("PATH", "/usr/bin:/bin")
    return request.param


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_tower_dir(
    base: Path,
    rel_dir: str,
    rows: list[tuple[str, str, int, str, str, str]],
    *,
    partial_name: str | None = None,
) -> tuple[Path, str]:
    d = base / rel_dir
    d.mkdir(parents=True, exist_ok=True)
    urls_lines = ["filename\turl\tbytes\tsha256\tlicense\tdate_utc"]
    sums_lines = []
    for fname, url, nbytes, sha, lic, dt_utc in rows:
        urls_lines.append(f"{fname}\t{url}\t{nbytes}\t{sha}\t{lic}\t{dt_utc}")
        sums_lines.append(f"{sha}  {fname}")
    (d / "SOURCE_URLS.tsv").write_text("\n".join(urls_lines) + "\n", encoding="utf-8")
    sums_bytes = ("\n".join(sums_lines) + "\n").encode("utf-8")
    (d / "SHA256SUMS").write_bytes(sums_bytes)
    if partial_name:
        (d / partial_name).write_bytes(b"partial")
    return d, _sha256_bytes(sums_bytes)


def test_summarize_tower_dir_and_nested_subtrees(tmp_path: Path):
    mod = _load()

    inc = tmp_path / "sources_incoming"
    sha_a = "a" * 64
    sha_b = "b" * 64
    _, expected_sums_sha = _write_tower_dir(
        inc,
        "batch5_proks/HOT1",
        [
            ("hot1_a.fna.gz", "https://example.org/hot1_a.fna.gz", 1000, sha_a, "CC-BY-4.0", "2026-11-26T10:00:00Z"),
            # Duplicate retry row for hot1_a.fna.gz should deduplicate by filename
            ("hot1_a.fna.gz", "https://example.org/hot1_a.fna.gz", 1200, sha_a, "CC-BY-4.0", "2026-11-26T11:00:00Z"),
            ("hot1_b.fna.gz", "https://example.org/hot1_b.fna.gz", 2800, sha_b, "CC-BY-4.0", "2026-11-27T04:30:00Z"),
        ],
    )
    summary = mod.summarize_tower_dir(inc / "batch5_proks/HOT1", root_label="sources_incoming", rel_dir="batch5_proks/HOT1")
    assert summary["files"] == 2
    assert summary["bytes"] == 4000
    assert summary["date_min"] == "2026-11-26"
    assert summary["date_max"] == "2026-11-27"
    assert summary["date_range"] == "2026-11-26..2026-11-27"
    assert summary["sha256sums_sha256"] == expected_sums_sha
    assert summary["streaming"] is False
    assert summary["stored_path"] == "tower:/mnt/user/research-bulk/bio-box/sources_incoming/batch5_proks/HOT1/"

    # Parent + child directory tree (mirrors batch7_cngb/EEMC1 & GPB_PHAGE having top-level tables + child subdirs)
    _write_tower_dir(
        inc,
        "batch7_cngb/EEMC1",
        [("meta.tsv", "ftp://ftp.cngb.org/meta.tsv", 100, "e" * 64, "CNGB-CC-BY", "2026-11-27T00:30:00Z")],
    )
    _write_tower_dir(
        inc,
        "batch7_cngb/EEMC1/isolates",
        [("iso1.fa.gz", "ftp://ftp.cngb.org/iso1.fa.gz", 500, "c" * 64, "CNGB-CC-BY", "2026-11-27T01:00:00Z")],
    )
    _write_tower_dir(
        inc,
        "batch7_cngb/EEMC1/mags",
        [("mag1.fa.gz", "ftp://ftp.cngb.org/mag1.fa.gz", 750, "d" * 64, "CNGB-CC-BY", "2026-11-28T05:00:00Z")],
        partial_name=".mag2.fa.gz.partial",
    )
    nested = mod.summarize_tower_dir(inc / "batch7_cngb/EEMC1", root_label="sources_incoming", rel_dir="batch7_cngb/EEMC1")
    assert nested["files"] == 3
    assert nested["bytes"] == 1350
    assert nested["date_range"] == "2026-11-27..2026-11-28"
    assert nested["streaming"] is True
    assert len(nested["subdirs"]) == 2

    inv = mod.build_inventory_from_roots(sources_incoming=inc)
    rolled = inv["dirs"]["sources_incoming/batch7_cngb/EEMC1"]
    assert rolled["files"] == 3
    assert rolled["bytes"] == 1350
    assert rolled["streaming"] is True


def test_reconcile_sources_yaml_preserves_urls_and_approved_verbatim(tmp_path: Path):
    mod = _load()

    inc = tmp_path / "sources_incoming"
    inc_s2 = tmp_path / "sources_incoming_s2"
    _, sums_hot = _write_tower_dir(
        inc,
        "batch5_proks/HOT1",
        [("hot1.tar.gz", "https://zenodo.org/api/records/12345/files/hot1.tar.gz/content", 8192, "1" * 64, "CC-BY-4.0", "2026-11-26T14:22:00Z")],
    )
    _, sums_asg = _write_tower_dir(
        inc_s2,
        "ASG1_FNA/assemblies",
        [
            ("GCA_000002005.1_genomic.fna.gz", "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/002/005/GCA_000002005.1_ASM3/GCA_000002005.1_ASM3_genomic.fna.gz", 1500, "9" * 64, "INSDC-open", "2026-11-28T04:00:00Z"),
        ],
    )
    _, sums_delta = _write_tower_dir(
        inc_s2,
        "GTDB_R232_DELTA",
        [
            ("GCA_000001005.1_genomic.fna.gz", "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/005/GCA_000001005.1_ASM1/GCA_000001005.1_ASM1_genomic.fna.gz", 3000, "2" * 64, "INSDC-open", "2026-11-27T18:00:00Z"),
            ("GCA_000001015.1_genomic.fna.gz", "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/015/GCA_000001015.1_ASM2/GCA_000001015.1_ASM2_genomic.fna.gz", 5000, "3" * 64, "INSDC-open", "2026-11-28T02:00:00Z"),
        ],
    )
    _write_tower_dir(
        inc_s2,
        "MGNIFY_REPS",
        [("MGYG000000001.fa.gz", "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/fish-gut/v1.0/species_catalogue/MGYG00000000/MGYG000000001/genome/MGYG000000001.fa.gz", 4096, "4" * 64, "EMBL-EBI-Terms-of-Use", "2026-11-28T06:00:00Z")],
        partial_name=".MGYG000000002.fa.gz.partial",
    )
    _write_tower_dir(
        inc,
        "batch6b_ncbi/O4_ORG/O4_CRYP",
        [("cryp.fna.gz", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=1", 600, "6" * 64, "INSDC", "2026-11-26T19:00:00Z")],
    )
    _write_tower_dir(
        inc,
        "batch6b_ncbi/TICK_RNA/nuccore",
        [("tick.fna.gz", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=2", 400, "7" * 64, "INSDC", "2026-11-26T19:15:00Z")],
    )
    _write_tower_dir(
        inc,
        "GLOBDB_R232_TABLES",
        [("tax.tsv.gz", "https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/tax.tsv.gz", 900, "8" * 64, "CC-BY-SA-4.0", "2026-11-26T18:23:00Z")],
    )

    orig_yaml = (
        "# Synthetic sources.yaml header comment\n"
        "version: 1\n\n"
        "sources:\n"
        "  - name: tower_r26_11_6_ncbi_fna\n"
        "    url: https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{d1}/{d2}/{d3}/{asm_dir}/{file}\n"
        "    license: NCBI / INSDC public domain\n"
        "    redistributable: true\n"
        "    approved: \"PI 2026-11-27 DEC-88 (exact quote kept!)\"\n"
        "    fetched: nothing yet (added 2026-11-27; fetched onto tower)\n"
        "    stored: will be tower:/mnt/user/research-bulk/bio-box/sources_incoming/{ASG1_FNA,GTDB_R232_DELTA}/\n\n"
        "  - name: tower_mgnify_nonhuman_reps\n"
        "    url: https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/{biome}/{ver}/species_catalogue/{prefix}/{mgyg}/genome/{file}\n"
        "    license: EMBL-EBI Terms of Use\n"
        "    redistributable: true\n"
        "    approved: \"PI 2026-11-27 DEC-88 (MGnify quote kept!)\"\n"
        "    fetched: nothing yet (added 2026-11-27)\n"
        "    stored: will be tower:/mnt/user/research-bulk/bio-box/sources_incoming/MGNIFY_REPS/\n\n"
        "  - name: tower_hot1\n"
        "    url: https://zenodo.org/api/records/12345/files/hot1.tar.gz/content\n"
        "    url_record: https://zenodo.org/records/12345\n"
        "    license: CC-BY-4.0\n"
        "    redistributable: true\n"
        "    approved: \"PI 2026-11-26 (Batch 5 HOT1 quote kept!)\"\n"
        "    fetched: nothing yet (added 2026-11-26; fetched onto tower)\n"
        "    stored: will be tower:/mnt/user/research-bulk/bio-box/sources_incoming/batch5_proks/HOT1/\n\n"
        "  - name: tower_unlanded_sibling\n"
        "    url: https://zenodo.org/api/records/54321/files/unlanded.tar.gz/content\n"
        "    license: CC-BY-4.0\n"
        "    redistributable: true\n"
        "    approved: \"PI 2026-11-26\"\n"
        "    fetched: nothing yet (unlanded sibling under batch5_proks must NOT pick up HOT1!)\n"
        "    stored: will be tower:/mnt/user/research-bulk/bio-box/sources_incoming/batch5_proks/UNLANDED_SIBLING/\n\n"
        "  - name: tower_batch6b_ncbi_o4_org\n"
        "    url: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id={ids}\n"
        "    license: INSDC\n"
        "    redistributable: true\n"
        "    approved: \"PI 2026-11-25\"\n"
        "    fetched: nothing yet\n"
        "    stored: will be tower:/mnt/user/research-bulk/bio-box/sources_incoming/batch6b_ncbi/O4_ORG/ and .../TICK_RNA/nuccore/\n\n"
        "  - name: tower_globdb_r232_tables\n"
        "    url: https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/tax.tsv.gz\n"
        "    license: CC-BY-SA-4.0\n"
        "    redistributable: false\n"
        "    approved: \"PI 2026-11-26\"\n"
        "    fetched: nothing yet\n"
        "    stored: tower:/mnt/user/research-bulk/bio-box/sources_incoming/GLOBDB_R232_TABLES/ (SHA256SUMS there); copy on pod\n\n"
        "  - name: biosino_node_meer_oep004067\n"
        "    url: https://www.biosino.org/node/project/detail/OEP004067\n"
        "    url_download: https://download.node.ac.cn/{path}\n"
        "    url_ftps: ftps://fms.biosino.org:2122/{path}\n"
        "    license: Biosino NODE Terms of Use (redistributable false under DEC-80 (3) until open licence verified)\n"
        "    redistributable: false\n"
        "    approved: \"PI 2026-11-28 DEC-89 (2)\"\n"
        "    fetched: nothing yet (pre-fetch portal entry)\n"
        "    stored: will be tower:/mnt/user/research-bulk/bio-box/sources_incoming/BIOSINO_MEER/\n"
    )

    inventory = mod.build_inventory_from_roots(
        sources_incoming=inc,
        sources_incoming_s2=inc_s2,
    )
    updated_yaml, touched = mod.reconcile_sources_yaml(
        orig_yaml,
        inventory,
        as_of="2026-11-28",
        expected_counts={"tower_mgnify_nonhuman_reps": 10},
    )

    assert set(touched) == {
        "tower_r26_11_6_ncbi_fna",
        "tower_mgnify_nonhuman_reps",
        "tower_hot1",
        "tower_batch6b_ncbi_o4_org",
        "tower_globdb_r232_tables",
    }
    assert mod.extract_url_forms(updated_yaml) == mod.extract_url_forms(orig_yaml)
    assert mod.extract_all_url_lines(updated_yaml) == mod.extract_all_url_lines(orig_yaml)

    orig_parsed = {e["name"]: e for e in yaml.safe_load(orig_yaml)["sources"]}
    new_parsed = {e["name"]: e for e in yaml.safe_load(updated_yaml)["sources"]}
    for name, orig_entry in orig_parsed.items():
        assert new_parsed[name]["approved"] == orig_entry["approved"]

    # Completed entries (including `{ASG1_FNA,GTDB_R232_DELTA}/` brace expansion) have exact measured counts, bytes, dates, SHA256SUMS sha256, and actual stored root
    assert "3 files" in new_parsed["tower_r26_11_6_ncbi_fna"]["fetched"]
    assert "9,500 B" in new_parsed["tower_r26_11_6_ncbi_fna"]["fetched"]
    assert "2026-11-27..2026-11-28" in new_parsed["tower_r26_11_6_ncbi_fna"]["fetched"]
    assert sums_asg in new_parsed["tower_r26_11_6_ncbi_fna"]["fetched"]
    assert sums_delta in new_parsed["tower_r26_11_6_ncbi_fna"]["fetched"]
    assert new_parsed["tower_r26_11_6_ncbi_fna"]["stored"] == (
        "tower:/mnt/user/research-bulk/bio-box/sources_incoming_s2/ASG1_FNA/; "
        "tower:/mnt/user/research-bulk/bio-box/sources_incoming_s2/GTDB_R232_DELTA/"
    )

    assert "1 file" in new_parsed["tower_hot1"]["fetched"]
    assert "8,192 B" in new_parsed["tower_hot1"]["fetched"]
    assert sums_hot in new_parsed["tower_hot1"]["fetched"]
    assert new_parsed["tower_hot1"]["stored"] == (
        "tower:/mnt/user/research-bulk/bio-box/sources_incoming/batch5_proks/HOT1/"
    )

    # Shorthand `.../TICK_RNA/nuccore/` expands alongside `batch6b_ncbi/O4_ORG/`
    assert "2 files" in new_parsed["tower_batch6b_ncbi_o4_org"]["fetched"]
    assert "1,000 B" in new_parsed["tower_batch6b_ncbi_o4_org"]["fetched"]
    assert "TICK_RNA/nuccore" in new_parsed["tower_batch6b_ncbi_o4_org"]["stored"]

    # Existing non-`will be` `stored:` string is preserved verbatim
    assert new_parsed["tower_globdb_r232_tables"]["stored"] == orig_parsed["tower_globdb_r232_tables"]["stored"]
    assert "900 B" in new_parsed["tower_globdb_r232_tables"]["fetched"]

    # Unlanded sibling under batch5_proks stays untouched
    assert new_parsed["tower_unlanded_sibling"]["fetched"].startswith("nothing yet")
    assert new_parsed["biosino_node_meer_oep004067"]["fetched"].startswith("nothing yet")

    # Streaming entry records partial count and active streaming note
    assert "1 of 10" in new_parsed["tower_mgnify_nonhuman_reps"]["fetched"]
    assert "streaming" in new_parsed["tower_mgnify_nonhuman_reps"]["fetched"].lower()
    assert new_parsed["tower_mgnify_nonhuman_reps"]["stored"] == (
        "tower:/mnt/user/research-bulk/bio-box/sources_incoming_s2/MGNIFY_REPS/"
    )

    # Verify fetch_to_tower.sh --check accepts the reconciled YAML (including ftps:// on Biosino)
    yaml_file = tmp_path / "reconciled_sources.yaml"
    yaml_file.write_text(updated_yaml, encoding="utf-8")
    job_file = tmp_path / "jobs.tsv"
    job_file.write_text(
        "GTDB_R232_DELTA\tGCA_000001005.1_genomic.fna.gz\t"
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/005/GCA_000001005.1_ASM1/GCA_000001005.1_ASM1_genomic.fna.gz\t-\tINSDC-open\n"
        "batch5_proks/HOT1\thot1.tar.gz\thttps://zenodo.org/api/records/12345/files/hot1.tar.gz/content\t-\tCC-BY-4.0\n"
        "BIOSINO_MEER\tmeer.tar.gz\thttps://download.node.ac.cn/OEP004067/meer.tar.gz\t-\tBiosino-NODE\n",
        encoding="utf-8",
    )
    res = subprocess.run(
        ["/bin/bash", str(FETCH_SCRIPT), "--check", "--sources", str(yaml_file), str(job_file)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert res.returncode == 0, f"stderr={res.stderr} stdout={res.stdout}"
    assert "3 URLs, all covered" in res.stdout


def test_verify_sources_yaml_catches_regressions_and_stale_entries(tmp_path: Path):
    mod = _load()

    inc = tmp_path / "sources_incoming"
    _write_tower_dir(
        inc,
        "O3_FUNGAL",
        [("fungi.tar", "https://zenodo.org/api/records/99/files/fungi.tar/content", 2048, "5" * 64, "CC-BY-4.0", "2026-11-26T09:00:00Z")],
    )
    inv = mod.build_inventory_from_roots(sources_incoming=inc)

    baseline = tmp_path / "baseline.yaml"
    baseline.write_text(
        "version: 1\n"
        "sources:\n"
        "  - name: tower_o3_fungal\n"
        "    url: https://zenodo.org/api/records/99/files/fungi.tar/content\n"
        "    license: CC-BY-4.0\n"
        "    redistributable: true\n"
        "    approved: \"PI 2026-11-25\"\n"
        "    fetched: nothing yet\n"
        "    stored: will be tower:/mnt/user/research-bulk/bio-box/sources_incoming/O3_FUNGAL/\n",
        encoding="utf-8",
    )

    # Stale entry (single or `{...}` brace-expanded) fails when inventory shows landed files
    with pytest.raises(mod.LedgerError, match="O3_FUNGAL"):
        mod.verify_sources_yaml(baseline, baseline_sources_path=baseline, inventory=inv)
    brace_stale = tmp_path / "brace_stale.yaml"
    brace_stale.write_text(
        baseline.read_text(encoding="utf-8").replace("O3_FUNGAL/", "{UNLANDED_X,O3_FUNGAL}/"),
        encoding="utf-8",
    )
    with pytest.raises(mod.LedgerError, match="O3_FUNGAL"):
        mod.verify_sources_yaml(brace_stale, baseline_sources_path=baseline, inventory=inv)

    reconciled_text, _ = mod.reconcile_sources_yaml(baseline.read_text(encoding="utf-8"), inv, as_of="2026-11-28")
    good = tmp_path / "good.yaml"
    good.write_text(reconciled_text, encoding="utf-8")
    mod.verify_sources_yaml(good, baseline_sources_path=baseline, inventory=inv)

    # Mutated URL form fails against baseline
    bad_url = tmp_path / "bad_url.yaml"
    bad_url.write_text(
        reconciled_text.replace("records/99/files", "records/100/files"),
        encoding="utf-8",
    )
    with pytest.raises(mod.LedgerError, match="URL form"):
        mod.verify_sources_yaml(bad_url, baseline_sources_path=baseline, inventory=inv)

    # Bare `"unknown"` redistributable fails
    bad_redist = tmp_path / "bad_redist.yaml"
    bad_redist.write_text(
        reconciled_text.replace("redistributable: true", "redistributable: \"unknown\""),
        encoding="utf-8",
    )
    with pytest.raises(mod.LedgerError, match="redistributable"):
        mod.verify_sources_yaml(bad_redist)

    # Host placeholder in URL fails
    bad_host = tmp_path / "bad_host.yaml"
    bad_host.write_text(
        reconciled_text + (
            "\n  - name: bad_portal\n"
            "    url: https://{host}/data/file.tar.gz\n"
            "    license: test\n"
            "    redistributable: false\n"
            "    approved: \"PI\"\n"
            "    fetched: nothing yet\n"
            "    stored: will be tower:/mnt/user/research-bulk/bio-box/sources_incoming/BAD/\n"
        ),
        encoding="utf-8",
    )
    with pytest.raises(mod.LedgerError, match="host"):
        mod.verify_sources_yaml(bad_host)


def test_verify_catalog_goal_and_worker_rules_and_cli(tmp_path: Path):
    mod = _load()

    stages_dir = tmp_path / "stages"
    tpl_dir = stages_dir / "templates"
    tpl_dir.mkdir(parents=True)
    (stages_dir / "G1.conf").write_text("STAGE=G1\nDEPTH=annotated\n", encoding="utf-8")
    (stages_dir / "C1.conf").write_text("STAGE=C1\nDEPTH=annotated\n", encoding="utf-8")
    (stages_dir / "C2.conf").write_text("STAGE=C2\nDEPTH=annotated\n", encoding="utf-8")
    (stages_dir / "C3.conf").write_text("STAGE=C3\nDEPTH=annotated\n", encoding="utf-8")
    (stages_dir / "GRARE1.conf").write_text("STAGE=GRARE1\nDEPTH=called\n", encoding="utf-8")
    (tpl_dir / "O4.conf").write_text("STAGE=O4\nDEPTH=called\n", encoding="utf-8")
    (tpl_dir / "JGIUG1.conf").write_text("STAGE=JGIUG1\nDEPTH=called\n", encoding="utf-8")

    # Synthetic staging directory with non-empty left_out.tsv for O4 and JGIUG1
    staging_root = tmp_path / "staging"
    (staging_root / "O4").mkdir(parents=True)
    (staging_root / "O4" / "left_out.tsv").write_text("id\treason\nx1\tskip\n", encoding="utf-8")
    (staging_root / "JGIUG1").mkdir(parents=True)
    (staging_root / "JGIUG1" / "JGIUG_left_for_later.tsv").write_text("id\treason\nu1\tspillover\n", encoding="utf-8")

    sources_path = tmp_path / "sources.yaml"
    sources_path.write_text(
        "version: 1\n"
        "sources:\n"
        "  - name: globdb_r226\n"
        "    url: https://fileshare.lisc.univie.ac.at/globdb/globdb_r226/{file}\n"
        "    license: CC0\n"
        "    redistributable: true\n"
        "    fetched: 2026-11-12 (*measured*)\n"
        "    stored: wilkcomp02:/stor/work/Wilke/if2453/all-the-data/globdb_r226/\n",
        encoding="utf-8",
    )

    # Extra pre-merge rows (`PROMGE1`, `BBSAG1`) and `C1` to `C3` range in §3 must be accepted;
    # mentioning `UNLISTED1` only inside a non-Stage column must NOT satisfy §3!
    catalog_goal = tmp_path / "CATALOG_GOAL.md"
    catalog_goal.write_text(
        "# Bio-Box Catalog Goal\n\n"
        "## 2. Where we are\n\n"
        "Updated through `DEC-88` and `DEC-89`.\n\n"
        "## 3. What is in, and what is on disk and queued\n\n"
        "| Source | State | Stage |\n"
        "|---|---|---|\n"
        "| GlobDB r226 | built | `G1` |\n"
        "| GTDB complete genomes | deferred | `C1` to `C3` |\n"
        "| Rare biosphere slice 1 (after `UNLISTED1`) | staged | `GRARE1` |\n"
        "| Organelles and fungi | staged | `O4` |\n"
        "| JGI uncultivated HQ non-R232 | staged | `JGIUG1` |\n"
        "| Prokaryotic MGEs (C-89) | staged | `PROMGE1` |\n"
        "| Marine SAGs (C-89) | staged | `BBSAG1` |\n\n"
        "## 4. What is not yet on disk\n\n"
        "### 4.1 Left out of datasets we already pulled\n\n"
        "| Dataset | Left out | Pickup |\n"
        "|---|---|---|\n"
        "| **`O4`** | `431` records in `left_out.tsv` | later |\n"
        "| **`JGIUG1`** | `295,882` genomes in `JGIUG_left_for_later.tsv.gz` | later |\n\n"
        "### 4.3 New sources, by priority\n\n"
        "| Source | Status |\n"
        "|---|---|\n"
        "| `Biosino` | queued |\n",
        encoding="utf-8",
    )

    worker_rules = tmp_path / "STAGING_WORKER_RULES.md"
    worker_rules.write_text(
        "# Staging Worker Rules\n\n"
        "At every fetch completion and staging-card merge, synchronize `sources.yaml` and "
        "`planning/CATALOG_GOAL.md` and verify with `python scripts/orchestrator/r26_sync_sources_ledger.py verify`.\n",
        encoding="utf-8",
    )

    def _run_verify() -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "verify",
                "--sources",
                str(sources_path),
                "--catalog-goal",
                str(catalog_goal),
                "--stages-dir",
                str(stages_dir),
                "--worker-rules",
                str(worker_rules),
                "--staging-root",
                str(staging_root),
                "--require-remainder",
                "O4",
                "--require-remainder",
                "JGIUG1",
            ],
            capture_output=True,
            text=True,
            check=False,
        )

    ok = _run_verify()
    assert ok.returncode == 0, f"stderr={ok.stderr} stdout={ok.stdout}"

    # Adding UNLISTED1.conf (only mentioned in prose, NOT in Stage column) fails verify
    (tpl_dir / "UNLISTED1.conf").write_text("STAGE=UNLISTED1\nDEPTH=called\n", encoding="utf-8")
    bad_conf = _run_verify()
    assert bad_conf.returncode == 1
    assert "UNLISTED1" in bad_conf.stderr
    (tpl_dir / "UNLISTED1.conf").unlink()

    # Dropping O4 remainder from §4.1 fails verify
    orig_cg = catalog_goal.read_text(encoding="utf-8")
    catalog_goal.write_text(
        orig_cg.replace("| **`O4`** | `431` records in `left_out.tsv` | later |\n", ""),
        encoding="utf-8",
    )
    bad_rem = _run_verify()
    assert bad_rem.returncode == 1
    assert "O4" in bad_rem.stderr
    catalog_goal.write_text(orig_cg, encoding="utf-8")

    # Dropping sync command from STAGING_WORKER_RULES.md fails verify
    worker_rules.write_text("# Staging Worker Rules\n\nNo sync rule here.\n", encoding="utf-8")
    bad_rules = _run_verify()
    assert bad_rules.returncode == 1
    assert "r26_sync_sources_ledger.py" in bad_rules.stderr
