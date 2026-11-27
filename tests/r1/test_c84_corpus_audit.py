"""Card C-84: three-tier corpus inventory and r26.11.6+ roadmap audit tests.

Verifies `scripts/orchestrator/c84_corpus_audit.py` on synthetic fixtures (normal PATH and tools-hidden
`PATH=/usr/bin:/bin`) and on the committed inventory `planning/measurements/2026-11-27_corpus_inventory_and_roadmap.json`
and note `planning/measurements/2026-11-27_corpus_inventory_and_roadmap.md`.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-84")

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "orchestrator" / "c84_corpus_audit.py"
INVENTORY = REPO / "planning" / "measurements" / "2026-11-27_corpus_inventory_and_roadmap.json"
ROADMAP_MD = REPO / "planning" / "measurements" / "2026-11-27_corpus_inventory_and_roadmap.md"

EXPECTED_DISPOSITIONS = {
    "built",
    "building_or_queued_r26_11_4",
    "staged_r26_11_5",
    "staged_later_tranche",
    "landed_on_tower_unstaged",
    "actively_fetching_on_tower",
    "failed_or_corrupt_archive",
    "wrong_payload_non_dna",
    "non_sequence_aux",
    "superseded_duplicate",
}


def _load_module():
    assert SCRIPT.is_file(), f"missing script: {SCRIPT}"
    spec = importlib.util.spec_from_file_location("c84_corpus_audit", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c84_corpus_audit"] = mod
    spec.loader.exec_module(mod)
    return mod


def _num(val: int | float, label: str = "measured", source: str = "unit test fixture") -> dict:
    return {"value": val, "label": label, "source": source}


def _make_synthetic_fixture(tmp_path: Path) -> tuple[Path, Path, Path, Path, dict]:
    stages_dir = tmp_path / "stages"
    templates_dir = stages_dir / "templates"
    templates_dir.mkdir(parents=True)
    (stages_dir / "G1.conf").write_text("STAGE=G1\nDEPTH=annotated\n", encoding="utf-8")
    (stages_dir / "CG1.conf").write_text("STAGE=CG1\nDEPTH=called\n", encoding="utf-8")
    (stages_dir / "GRARE1.conf").write_text("STAGE=GRARE1\nDEPTH=called\n", encoding="utf-8")
    (stages_dir / "C1.conf").write_text("STAGE=C1\nDEPTH=annotated\n", encoding="utf-8")
    (templates_dir / "GRARE1.conf").write_text("STAGE=GRARE1\nDEPTH=called\n", encoding="utf-8")
    (templates_dir / "ET1.conf").write_text("STAGE=ET1\nDEPTH=called\n", encoding="utf-8")

    jobs_dir = tmp_path / "tower_jobs"
    jobs_dir.mkdir()
    (jobs_dir / "01_test.tsv").write_text(
        "# comment line\n"
        "O4_ORG\tchunk1.fna\thttps://example.org/chunk1.fna\t-\tINSDC\n"
        "O4_ORG\tchunk1.gbff\thttps://example.org/chunk1.gbff\t-\tINSDC\n",
        encoding="utf-8",
    )
    (jobs_dir / "22_eemc1.tsv").write_text(
        "batch7_cngb/EEMC1/isolates\tBGI1.fasta.gz\tftp://ftp.cngb.org/BGI1.fasta.gz\tabcd\tCNGBdb-CC-BY-4.0\n",
        encoding="utf-8",
    )

    sources_path = tmp_path / "sources.yaml"
    sources_path.write_text(
        "version: 1\n"
        "sources:\n"
        "  - name: globdb_r226\n"
        "    url: https://example.org/globdb/{file}\n"
        "  - name: tower_asg1\n"
        "    url: https://example.org/asg1/Data.zip\n"
        "  - name: tower_eemc1\n"
        "    url: https://db.cngb.org/data_resources/project/CNP0007106/\n",
        encoding="utf-8",
    )

    inv = {
        "schema_version": "c84.v1",
        "as_of": "2026-11-27",
        "observed": {
            "pod_root_stages": ["S1", "G1"],
            "pod_s2_inputs_dirs": ["CG1", "GV1"],
            "pod_biobox_staging_dirs": ["GRARE"],
            "tower_sources_incoming_dirs": [
                "EGV1",
                "GHM1_VHM1",
                "batch7_cngb/EEMC1",
                "DSV70",
                "ASG1",
                "CANDI1",
                "batch5_proks/EMERGE",
            ],
            "tower_sources_dirs": ["G1", "CG1", "GRARE1"],
            "stage_confs": ["C1", "CG1", "G1", "GRARE1"],
            "template_confs": ["ET1", "GRARE1"],
            "tower_job_lists": {
                "01_test.tsv": _num(2, "measured", "fixture 01_test.tsv line count"),
                "22_eemc1.tsv": _num(1, "measured", "fixture 22_eemc1.tsv line count"),
            },
            "sources_yaml_names": ["globdb_r226", "tower_asg1", "tower_eemc1"],
        },
        "summary_tiers": [
            {
                "tier_id": "state_1_built_pod",
                "title": "Built in pod root",
                "count": _num(2, "measured", "fixture root stages"),
                "disk_bytes": _num(1000, "derived", "sum of fixture stages"),
                "genomic_bp": _num(5000, "measured", "fixture contig_index"),
            }
        ],
        "pod_stages_and_assets": [
            {"name": "S1", "kind": "root_stage", "disposition": "built", "release_or_wave": "r26.11.3", "notes": "built phage stage", "bases_bp": _num(1000)},
            {"name": "G1", "kind": "conf_and_root_stage", "disposition": "built", "release_or_wave": "r26.11.3", "notes": "built archaeal stage", "bases_bp": _num(4000)},
            {"name": "CG1", "kind": "conf_and_s2_input", "disposition": "building_or_queued_r26_11_4", "release_or_wave": "r26.11.4", "notes": "wave-6 building", "bases_bp": _num(2000)},
            {"name": "GV1", "kind": "s2_input", "disposition": "staged_r26_11_5", "release_or_wave": "r26.11.5", "notes": "staged giant viruses", "bases_bp": _num(1500)},
            {"name": "GRARE1", "kind": "conf_and_template", "disposition": "staged_r26_11_5", "release_or_wave": "r26.11.5", "notes": "staged rare phyla", "bases_bp": _num(3000)},
            {"name": "GRARE", "kind": "biobox_staging_dir", "disposition": "staged_r26_11_5", "release_or_wave": "r26.11.5", "notes": "GRARE1..3 staging dir", "bytes": _num(900)},
            {"name": "C1", "kind": "conf_only", "disposition": "staged_later_tranche", "release_or_wave": "deferred_strains", "notes": "DEC-70 deferred strain stage", "bases_bp": _num(800, "estimate", "fixture estimate")},
            {"name": "ET1", "kind": "template_only", "disposition": "staged_later_tranche", "release_or_wave": "unassigned_template", "notes": "C-68 placeholder template", "bases_bp": _num(0, "derived", "placeholder")},
        ],
        "tower_sources_incoming": [
            {"dir": "EGV1", "disposition": "staged_r26_11_5", "files": _num(1), "bytes": _num(100), "summary": "staged in r26.11.5"},
            {"dir": "GHM1_VHM1", "disposition": "landed_on_tower_unstaged", "files": _num(3), "bytes": _num(300), "summary": "ready to stage for r26.11.6"},
            {"dir": "batch7_cngb/EEMC1", "disposition": "actively_fetching_on_tower", "files": _num(10), "bytes": _num(500), "summary": "list 22 running"},
            {"dir": "DSV70", "disposition": "failed_or_corrupt_archive", "files": _num(4), "bytes": _num(200), "summary": "DSV70_mag.tar.gz FAIL_GZIP"},
            {"dir": "ASG1", "disposition": "wrong_payload_non_dna", "files": _num(1), "bytes": _num(400), "summary": "Data.zip holds proteomes only"},
            {"dir": "CANDI1", "disposition": "non_sequence_aux", "files": _num(2), "bytes": _num(50), "summary": "workflow code zip only"},
            {"dir": "batch5_proks/EMERGE", "disposition": "superseded_duplicate", "files": _num(2), "bytes": _num(600), "summary": "duplicate of emerge_stordalen_v2_0_1"},
        ],
        "tower_sources_mirrors": [
            {"dir": "G1", "disposition": "built", "files": _num(3), "bytes": _num(400), "stage_or_role": "G1 mirror"},
            {"dir": "CG1", "disposition": "building_or_queued_r26_11_4", "files": _num(3), "bytes": _num(250), "stage_or_role": "CG1 mirror"},
            {"dir": "GRARE1", "disposition": "staged_r26_11_5", "files": _num(4), "bytes": _num(350), "stage_or_role": "GRARE1 mirror"},
        ],
        "tower_job_lists": [
            {"file": "01_test.tsv", "disposition": "landed_on_tower_unstaged", "jobs": _num(2), "status_detail": "landed on tower in 20-list run"},
            {"file": "22_eemc1.tsv", "disposition": "actively_fetching_on_tower", "jobs": _num(1), "status_detail": "actively running with --jobs 4"},
        ],
        "sources_yaml_reconciliation": [
            {"name": "globdb_r226", "disposition": "built", "reason_code": "built_or_staged_on_pod", "detail": "G1 built, GART1/2 and GRARE1..3 staged"},
            {"name": "tower_asg1", "disposition": "wrong_payload_non_dna", "reason_code": "wrong_payload_or_failed_archive", "detail": "Data.zip lacks nucleotide genomes; card C-82"},
            {"name": "tower_eemc1", "disposition": "actively_fetching_on_tower", "reason_code": "actively_fetching_or_queued_retry_on_tower", "detail": "CNGB ftp list 22 running"},
        ],
        "roadmap_priorities": [
            {
                "priority_group": "P1",
                "id": "GHM1",
                "title": "Global Hot Spring Microbiomes",
                "expected_gbp": {"value": 24.0, "low": 21.0, "high": 27.0, "label": "estimate", "source": "card C-71"},
                "how_to_grab": "already on tower under sources_incoming/GHM1_VHM1; relay_to_pod.sh",
            }
        ],
    }
    inv_path = tmp_path / "inventory.json"
    inv_path.write_text(json.dumps(inv, indent=2), encoding="utf-8")
    return inv_path, sources_path, jobs_dir, stages_dir, inv


def _run_cli(args: list[str], *, hide_pod_tools: bool = False) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ)
    if hide_pod_tools:
        env["PATH"] = "/usr/bin:/bin"
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )


def test_valid_synthetic_fixture_passes_normally_and_with_tools_hidden(tmp_path: Path):
    mod = _load_module()
    assert set(mod.DISPOSITIONS) == EXPECTED_DISPOSITIONS
    inv_path, sources_path, jobs_dir, stages_dir, _ = _make_synthetic_fixture(tmp_path)

    for hide in (False, True):
        res = _run_cli(
            [
                "--check",
                "--inventory",
                str(inv_path),
                "--sources",
                str(sources_path),
                "--jobs-dir",
                str(jobs_dir),
                "--stages-dir",
                str(stages_dir),
            ],
            hide_pod_tools=hide,
        )
        assert res.returncode == 0, f"hide={hide} stderr={res.stderr} stdout={res.stdout}"
        assert "unclassified=0" in res.stdout
        assert "number_problems=0" in res.stdout


def test_unclassified_or_missing_tower_directory_fails_exit_1(tmp_path: Path):
    inv_path, sources_path, jobs_dir, stages_dir, inv = _make_synthetic_fixture(tmp_path)

    # Case 1: observed nested directory under sources_incoming is missing from classification table
    bad1 = copy.deepcopy(inv)
    bad1["observed"]["tower_sources_incoming_dirs"].append("batch6_mges/PROMGE2")
    p1 = tmp_path / "bad_missing_incoming.json"
    p1.write_text(json.dumps(bad1), encoding="utf-8")
    r1 = _run_cli(["--check", "--inventory", str(p1), "--sources", str(sources_path), "--jobs-dir", str(jobs_dir), "--stages-dir", str(stages_dir)])
    assert r1.returncode == 1
    assert "batch6_mges/PROMGE2" in r1.stderr

    # Case 2: observed directory under tower sources/ mirror is missing
    bad2 = copy.deepcopy(inv)
    bad2["observed"]["tower_sources_dirs"].append("SYMX_UNMAPPED")
    p2 = tmp_path / "bad_missing_mirror.json"
    p2.write_text(json.dumps(bad2), encoding="utf-8")
    r2 = _run_cli(["--check", "--inventory", str(p2), "--sources", str(sources_path), "--jobs-dir", str(jobs_dir), "--stages-dir", str(stages_dir)])
    assert r2.returncode == 1
    assert "SYMX_UNMAPPED" in r2.stderr

    # Case 3: invalid disposition string
    bad3 = copy.deepcopy(inv)
    bad3["tower_sources_incoming"][0]["disposition"] = "mystery_bucket"
    p3 = tmp_path / "bad_disp.json"
    p3.write_text(json.dumps(bad3), encoding="utf-8")
    r3 = _run_cli(["--check", "--inventory", str(p3), "--sources", str(sources_path), "--jobs-dir", str(jobs_dir), "--stages-dir", str(stages_dir)])
    assert r3.returncode == 1
    assert "mystery_bucket" in r3.stderr


def test_unclassified_stage_conf_template_job_list_or_source_fails_exit_1(tmp_path: Path):
    inv_path, sources_path, jobs_dir, stages_dir, inv = _make_synthetic_fixture(tmp_path)

    # Unmapped stage .conf file added on disk
    (stages_dir / "UNMAPPED1.conf").write_text("STAGE=UNMAPPED1\n", encoding="utf-8")
    r_conf = _run_cli(["--check", "--inventory", str(inv_path), "--sources", str(sources_path), "--jobs-dir", str(jobs_dir), "--stages-dir", str(stages_dir)])
    assert r_conf.returncode == 1
    assert "UNMAPPED1" in r_conf.stderr
    (stages_dir / "UNMAPPED1.conf").unlink()

    # Unmapped template .conf file added on disk
    (stages_dir / "templates" / "UNMAPPED_TPL.conf").write_text("STAGE=UNMAPPED_TPL\n", encoding="utf-8")
    r_tpl = _run_cli(["--check", "--inventory", str(inv_path), "--sources", str(sources_path), "--jobs-dir", str(jobs_dir), "--stages-dir", str(stages_dir)])
    assert r_tpl.returncode == 1
    assert "UNMAPPED_TPL" in r_tpl.stderr
    (stages_dir / "templates" / "UNMAPPED_TPL.conf").unlink()

    # Unmapped job list .tsv added on disk
    extra_job = jobs_dir / "99_unmapped.tsv"
    extra_job.write_text("FOO\ta.fna\thttps://example.org/a.fna\t-\tCC0\n", encoding="utf-8")
    r_job = _run_cli(["--check", "--inventory", str(inv_path), "--sources", str(sources_path), "--jobs-dir", str(jobs_dir), "--stages-dir", str(stages_dir)])
    assert r_job.returncode == 1
    assert "99_unmapped.tsv" in r_job.stderr
    extra_job.unlink()

    # Wrong job count recorded for a job list
    bad_jobs = copy.deepcopy(inv)
    bad_jobs["tower_job_lists"][0]["jobs"]["value"] = 999
    p_jobs = tmp_path / "bad_job_count.json"
    p_jobs.write_text(json.dumps(bad_jobs), encoding="utf-8")
    r_jcnt = _run_cli(["--check", "--inventory", str(p_jobs), "--sources", str(sources_path), "--jobs-dir", str(jobs_dir), "--stages-dir", str(stages_dir)])
    assert r_jcnt.returncode == 1
    assert "01_test.tsv" in r_jcnt.stderr

    # Unmapped sources.yaml entry
    bad_src = copy.deepcopy(inv)
    bad_src["sources_yaml_reconciliation"].pop()
    p_src = tmp_path / "bad_src.json"
    p_src.write_text(json.dumps(bad_src), encoding="utf-8")
    r_src = _run_cli(["--check", "--inventory", str(p_src), "--sources", str(sources_path), "--jobs-dir", str(jobs_dir), "--stages-dir", str(stages_dir)])
    assert r_src.returncode == 1
    assert "tower_eemc1" in r_src.stderr


def test_unlabelled_or_unsourced_numbers_fail_exit_1(tmp_path: Path):
    _, sources_path, jobs_dir, stages_dir, inv = _make_synthetic_fixture(tmp_path)

    # 1. Bare unlabelled number
    bad_bare = copy.deepcopy(inv)
    bad_bare["raw_unlabelled_count"] = 42
    p_bare = tmp_path / "bad_bare.json"
    p_bare.write_text(json.dumps(bad_bare), encoding="utf-8")
    r_bare = _run_cli(["--check", "--inventory", str(p_bare), "--sources", str(sources_path), "--jobs-dir", str(jobs_dir), "--stages-dir", str(stages_dir)])
    assert r_bare.returncode == 1
    assert "raw_unlabelled_count" in r_bare.stderr

    # 2. Invalid label ("published" is not one of AGENTS.md §6 measured/derived/estimate)
    bad_lbl = copy.deepcopy(inv)
    bad_lbl["tower_sources_incoming"][0]["bytes"]["label"] = "published"
    p_lbl = tmp_path / "bad_lbl.json"
    p_lbl.write_text(json.dumps(bad_lbl), encoding="utf-8")
    r_lbl = _run_cli(["--check", "--inventory", str(p_lbl), "--sources", str(sources_path), "--jobs-dir", str(jobs_dir), "--stages-dir", str(stages_dir)])
    assert r_lbl.returncode == 1
    assert "published" in r_lbl.stderr

    # 3. Empty source string
    bad_src = copy.deepcopy(inv)
    bad_src["tower_sources_incoming"][0]["bytes"]["source"] = "   "
    p_src = tmp_path / "bad_empty_source.json"
    p_src.write_text(json.dumps(bad_src), encoding="utf-8")
    r_src = _run_cli(["--check", "--inventory", str(p_src), "--sources", str(sources_path), "--jobs-dir", str(jobs_dir), "--stages-dir", str(stages_dir)])
    assert r_src.returncode == 1
    assert "source" in r_src.stderr


def test_snapshot_only_ignores_new_disk_files_and_catches_observed_mismatches(tmp_path: Path):
    inv_path, sources_path, jobs_dir, stages_dir, inv = _make_synthetic_fixture(tmp_path)

    # Simulate future cards (e.g. A-32, C-82) adding new files/entries on disk without editing the snapshot JSON
    (stages_dir / "RMD1.conf").write_text("STAGE=RMD1\nDEPTH=called\n", encoding="utf-8")
    (stages_dir / "templates" / "GBRMGD.conf").write_text("STAGE=GBRMGD\nDEPTH=called\n", encoding="utf-8")
    (jobs_dir / "23_asg1_fna.tsv").write_text(
        "ASG1_FNA\ta.fna.gz\thttps://ftp.ncbi.nlm.nih.gov/a.fna.gz\t-\tINSDC\n",
        encoding="utf-8",
    )
    sources_path.write_text(
        sources_path.read_text(encoding="utf-8")
        + "  - name: ncbi_c82_asg1_fna\n"
        + "    url: https://ftp.ncbi.nlm.nih.gov/genomes/all/{file}\n",
        encoding="utf-8",
    )

    # Live --check fails because the new on-disk files are not in the historical snapshot
    r_live = _run_cli(
        [
            "--check",
            "--inventory",
            str(inv_path),
            "--sources",
            str(sources_path),
            "--jobs-dir",
            str(jobs_dir),
            "--stages-dir",
            str(stages_dir),
        ]
    )
    assert r_live.returncode == 1

    # --check --snapshot-only passes even with new on-disk files (and even when on-disk paths do not exist)
    for hide in (False, True):
        r_snap = _run_cli(
            [
                "--check",
                "--snapshot-only",
                "--inventory",
                str(inv_path),
                "--sources",
                str(tmp_path / "nonexistent_sources.yaml"),
                "--jobs-dir",
                str(tmp_path / "nonexistent_jobs"),
                "--stages-dir",
                str(tmp_path / "nonexistent_stages"),
            ],
            hide_pod_tools=hide,
        )
        assert r_snap.returncode == 0, f"hide={hide} stderr={r_snap.stderr} stdout={r_snap.stdout}"
        assert "unclassified=0" in r_snap.stdout
        assert "number_problems=0" in r_snap.stdout

    # --snapshot-only catches an unclassified stage conf in observed.stage_confs
    bad_sconf = copy.deepcopy(inv)
    bad_sconf["observed"]["stage_confs"].append("MISSING_STAGE_CONF")
    p_sconf = tmp_path / "bad_snap_sconf.json"
    p_sconf.write_text(json.dumps(bad_sconf), encoding="utf-8")
    r_sconf = _run_cli(["--check", "--snapshot-only", "--inventory", str(p_sconf)])
    assert r_sconf.returncode == 1
    assert "MISSING_STAGE_CONF" in r_sconf.stderr

    # --snapshot-only catches an unclassified template conf in observed.template_confs
    bad_tconf = copy.deepcopy(inv)
    bad_tconf["observed"]["template_confs"].append("MISSING_TPL_CONF")
    p_tconf = tmp_path / "bad_snap_tconf.json"
    p_tconf.write_text(json.dumps(bad_tconf), encoding="utf-8")
    r_tconf = _run_cli(["--check", "--snapshot-only", "--inventory", str(p_tconf)])
    assert r_tconf.returncode == 1
    assert "MISSING_TPL_CONF" in r_tconf.stderr

    # --snapshot-only catches job count mismatch against observed.tower_job_lists
    bad_jcnt = copy.deepcopy(inv)
    bad_jcnt["tower_job_lists"][0]["jobs"]["value"] = 77
    p_jcnt = tmp_path / "bad_snap_jcnt.json"
    p_jcnt.write_text(json.dumps(bad_jcnt), encoding="utf-8")
    r_jcnt = _run_cli(["--check", "--snapshot-only", "--inventory", str(p_jcnt)])
    assert r_jcnt.returncode == 1
    assert "01_test.tsv" in r_jcnt.stderr

    # --snapshot-only catches missing source reconciliation against observed.sources_yaml_names
    bad_src = copy.deepcopy(inv)
    bad_src["observed"]["sources_yaml_names"].append("unmapped_snapshot_source")
    p_src = tmp_path / "bad_snap_src.json"
    p_src.write_text(json.dumps(bad_src), encoding="utf-8")
    r_src = _run_cli(["--check", "--snapshot-only", "--inventory", str(p_src)])
    assert r_src.returncode == 1
    assert "unmapped_snapshot_source" in r_src.stderr


def test_committed_inventory_and_roadmap_pass_all_checks():
    assert INVENTORY.is_file(), f"missing committed inventory: {INVENTORY}"
    assert ROADMAP_MD.is_file(), f"missing committed roadmap note: {ROADMAP_MD}"

    for hide in (False, True):
        res = _run_cli(["--check", "--snapshot-only", "--inventory", str(INVENTORY)], hide_pod_tools=hide)
        assert res.returncode == 0, f"hide={hide} stderr={res.stderr} stdout={res.stdout}"
        assert "unclassified=0" in res.stdout
        assert "number_problems=0" in res.stdout

    raw_json = INVENTORY.read_text(encoding="utf-8")
    data = json.loads(raw_json)
    obs = data["observed"]
    assert isinstance(obs.get("stage_confs"), list) and len(obs["stage_confs"]) == 48
    assert isinstance(obs.get("template_confs"), list) and len(obs["template_confs"]) == 24
    assert isinstance(obs.get("tower_job_lists"), dict) and len(obs["tower_job_lists"]) == 23
    assert isinstance(obs.get("sources_yaml_names"), list) and len(obs["sources_yaml_names"]) == 227

    stage_names = {entry["name"] for entry in data["pod_stages_and_assets"]}
    r3_stages = {"S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1"}
    r4_stages = {
        "V2", "PL2", "PL3", "CG1", "MGE1", "SAT1", "ECE1", "O1", "O2", "LPL1",
        "SYM1", "LRV1", "LGRT1", "MGE2", "MV1", "V3", "LR1", "GART1", "GART2",
    }
    r5_stages = {
        "GRARE1", "GRARE2", "GRARE3", "LR2", "SAG1", "NX1", "GV1", "EGV1", "ECE2", "MIR1",
        "SAG2", "DPN1", "SYM2", "LRP1", "MIR2", "MV2", "MVEUKV", "KAVHIFI", "INPH2604",
        "GWV1", "SAG3", "WET1", "RMD1", "GBRMGD",
    }
    assert r3_stages <= stage_names
    assert r4_stages <= stage_names
    assert r5_stages <= stage_names

    md = ROADMAP_MD.read_text(encoding="utf-8")
    for heading in (
        "Executive Summary",
        "Section 1",
        "Section 2",
        "Section 3",
        "Section 4",
        "Priority 1",
        "Priority 2",
        "Priority 3",
        "Priority 4",
        "Priority 5",
        "Priority 6",
    ):
        assert heading in md, f"missing heading {heading!r} in {ROADMAP_MD.name}"

    # Verify ASG1 attribution (Appler/Baker, not Eme) and R-24 DSV70 diagnosis (PRJNA1244896 / S-BSST2227)
    assert "Eme " not in raw_json and "Eme " not in md
    assert "Appler" in raw_json and "Appler" in md
    assert "PRJNA1244896" in raw_json and "PRJNA1244896" in md
    assert "S-BSST2227" in raw_json and "S-BSST2227" in md
