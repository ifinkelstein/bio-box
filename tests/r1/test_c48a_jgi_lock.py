"""Card C-48a: unit tests for IMG/M jqdt JSON -> TSV.GZ conversion, metadata bundling,
and /stor/work archive + sources.yaml verification.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import sys
import tarfile
import zipfile
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("C-48a")  # ring R1 by its directory (tests/conftest.py)

ROOT = Path(__file__).resolve().parents[2]
IMG_M_SCRIPT = ROOT / "scripts" / "stress" / "c48a_img_m_tsv.py"
VERIFY_SCRIPT = ROOT / "scripts" / "stress" / "c48a_verify_stor.py"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None, f"missing script: {path}"
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_clean_and_parse_img_m_cells():
    m = _load("c48a_img_m_tsv", IMG_M_SCRIPT)
    assert m.clean_cell("*Microbiome#==#3300030970") == "*Microbiome"
    assert m.clean_cell("Ga0075381#==##==#analysis_projects?id=Ga0075381") == "Ga0075381"
    assert (
        m.clean_cell("<font color='red'> DOE Joint Genome Institute </font> (JGI)#==#DOE Joint Genome Institute (JGI)")
        == "DOE Joint Genome Institute (JGI)"
    )
    assert m.clean_cell('<a href="main.cgi?taxon_oid=12345">Sample &amp; 1</a>') == "Sample & 1"
    assert m.clean_cell(48.7184) == "48.7184"
    assert m.clean_cell(None) == ""
    assert m.clean_cell("line1\tline2\nline3") == "line1 line2 line3"

    text, extra = m.parse_cell('<a href="main.cgi?taxon_oid=12345">Sample 1</a>')
    assert text == "Sample 1"
    assert extra == "main.cgi?taxon_oid=12345"

    text2, extra2 = m.parse_cell("Ga0075381#==##==#analysis_projects?id=Ga0075381")
    assert text2 == "Ga0075381"
    assert extra2 == "analysis_projects?id=Ga0075381"


def test_pii_column_guard_blocks_email_and_contact():
    m = _load("c48a_img_m_tsv", IMG_M_SCRIPT)
    m.assert_no_pii_columns(["IMG Genome ID", "Study Name", "Geographic Location", "Host Name"])
    for bad in ["Contact Email", "PI_Email", "Submitter Name", "Principal Investigator", "Phone Number"]:
        with pytest.raises(ValueError, match="PII"):
            m.assert_no_pii_columns(["IMG Genome ID", bad])


def test_convert_img_m_json_deterministic_and_supports_data_and_aadata(tmp_path: Path):
    m = _load("c48a_img_m_tsv", IMG_M_SCRIPT)
    headers = ["", "NCBI Domain", "IMG Genome ID", "Sequencing Center", "GOLD Study ID", "Latitude"]
    row_a = [
        "*Microbiome#==#3300000002",
        "*Microbiome#==#3300000002",
        "3300000002#==#3300000002",
        "<font color='red'> DOE Joint Genome Institute </font> (JGI)#==#DOE Joint Genome Institute (JGI)",
        "Gs0000002#==##==#study?id=Gs0000002",
        12.5,
    ]
    row_b = [
        "*Microbiome#==#3300000001",
        "*Microbiome#==#3300000001",
        "3300000001#==#3300000001",
        "Broad\tInstitute\nUSA",
        "Gs0000001#==##==#study?id=Gs0000001",
        -3.25,
    ]

    p1 = tmp_path / "in1.json"
    p2 = tmp_path / "in2.json"
    p1.write_text(json.dumps({"data": [row_a, row_b]}), encoding="utf-8")
    p2.write_text(json.dumps({"aaData": [row_b, row_a]}), encoding="utf-8")

    out1 = tmp_path / "out1.tsv.gz"
    out2 = tmp_path / "out2.tsv.gz"
    s1 = m.convert_img_m_json(p1, out1, headers)
    s2 = m.convert_img_m_json(p2, out2, headers)

    assert s1["rows"] == 2
    assert s1["cols"] == 5
    assert s1["sha256"] == s2["sha256"]
    assert out1.read_bytes() == out2.read_bytes()

    with gzip.open(out1, "rt", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    assert lines[0] == "NCBI Domain\tIMG Genome ID\tSequencing Center\tGOLD Study ID\tLatitude"
    assert lines[1] == "*Microbiome\t3300000001\tBroad Institute USA\tGs0000001\t-3.25"
    assert lines[2] == "*Microbiome\t3300000002\tDOE Joint Genome Institute (JGI)\tGs0000002\t12.5"


def test_resolve_headers_handles_72_and_70_col_img_m_exports():
    m = _load("c48a_img_m_tsv", IMG_M_SCRIPT)
    cols_73 = ["", "NCBI Domain"] + [f"col_{i}" for i in range(63)] + [
        "Predicted Viruses",
        "Predicted Plasmids",
        "JGI Project ID / ITS SP ID",
        "ITS AP ID",
        "ITS SP ID",
        "SRA ID",
        "SRA Run",
        "JGI Data Utilization Status",
    ]
    assert len(cols_73) == 73
    h72 = m.resolve_headers(cols_73, raw_ncols=73)
    assert len(h72) == 73
    assert "Predicted Viruses" in h72 and "Predicted Plasmids" in h72

    h70 = m.resolve_headers(cols_73, raw_ncols=71)
    assert len(h70) == 71
    assert "Predicted Viruses" not in h70 and "Predicted Plasmids" not in h70
    assert h70[-1] == "JGI Data Utilization Status"


def test_cookie_exclusion_and_deterministic_metadata_tar(tmp_path: Path):
    m = _load("c48a_img_m_tsv", IMG_M_SCRIPT)
    assert m.is_excluded_bundle_path(Path(".cookie.cfg"))
    assert m.is_excluded_bundle_path(Path("sub/.cookie"))
    assert m.is_excluded_bundle_path(Path("foo.part"))
    assert m.is_excluded_bundle_path(Path(".tmp.bundle.tar"))
    assert not m.is_excluded_bundle_path(Path("gold/goldData.xlsx"))
    assert not m.is_excluded_bundle_path(Path("img_m/img_m_microbiomes.tsv.gz"))

    src = tmp_path / "tables"
    (src / "gold").mkdir(parents=True)
    (src / "img_m").mkdir(parents=True)
    (src / "gold" / "goldData.xlsx").write_bytes(b"PK\x03\x04fake-xlsx")
    (src / "gold" / ".cookie.cfg").write_text("cookie = secret", encoding="utf-8")
    (src / "img_m" / "microbiomes.tsv.gz").write_bytes(b"\x1f\x8b\x08fake-gz")
    (src / "img_m" / "temp.part").write_bytes(b"partial")

    t1 = tmp_path / "b1.tar"
    t2 = tmp_path / "b2.tar"
    m.build_metadata_tar([src / "gold", src / "img_m"], t1, min_bytes=0)
    m.build_metadata_tar([src / "img_m", src / "gold"], t2, min_bytes=0)

    assert t1.read_bytes() == t2.read_bytes()
    with tarfile.open(t1, "r:") as tf:
        names = [mem.name for mem in tf.getmembers()]
        assert names == ["gold/goldData.xlsx", "img_m/microbiomes.tsv.gz"]
        for mem in tf.getmembers():
            assert mem.uid == 0 and mem.gid == 0
            assert mem.mtime == m.DEFAULT_TAR_MTIME


def test_verify_stor_enforces_64mib_floor_sha256_and_archive_headers(tmp_path: Path):
    v = _load("c48a_verify_stor", VERIFY_SCRIPT)
    assert v.MIN_STOR_BYTES == 67_108_864
    assert v.validate_sha256("a" * 64)
    assert not v.validate_sha256("A" * 64)
    assert not v.validate_sha256("a" * 63)
    assert not v.validate_sha256("g" * 64)

    # Create a small valid zip and a small valid tar
    zpath = tmp_path / "sample.zip"
    with zipfile.ZipFile(zpath, "w", compression=zipfile.ZIP_STORED) as zf:
        zf.writestr("folder/data.tsv", b"col1\tcol2\n1\t2\n")

    tpath = tmp_path / "sample.tar"
    payload = tmp_path / "payload.bin"
    payload.write_bytes(b"0" * 1024)
    with tarfile.open(tpath, "w:") as tf:
        tf.add(payload, arcname="sub/payload.bin")
    payload.unlink()

    zmembers = v.inspect_archive_members(zpath)
    assert zmembers == [{"name": "folder/data.tsv", "size": 14}]
    tmembers = v.inspect_archive_members(tpath)
    assert tmembers == [{"name": "sub/payload.bin", "size": 1024}]

    # 64 MiB floor check rejects tiny files when default min_bytes is used
    issues = v.check_stor_dir(tmp_path)
    assert any("below 64 MiB" in msg for msg in issues)

    # With lowered min_bytes for unit test, clean directory passes
    assert v.check_stor_dir(tmp_path, min_bytes=100) == []

    # Cookie file inside directory fails check_stor_dir
    cookie = tmp_path / ".cookie.cfg"
    cookie.write_bytes(b"x" * 200)
    issues2 = v.check_stor_dir(tmp_path, min_bytes=100)
    assert any("cookie" in msg.lower() for msg in issues2)
    cookie.unlink()

    # Cookie member inside zip fails inspect_archive_members
    bad_zip = tmp_path / "bad.zip"
    with zipfile.ZipFile(bad_zip, "w") as zf:
        zf.writestr(".cookie.cfg", b"secret")
    with pytest.raises(ValueError, match="cookie"):
        v.inspect_archive_members(bad_zip)


def test_verify_sources_yaml_locks_archives(tmp_path: Path):
    v = _load("c48a_verify_stor", VERIFY_SCRIPT)
    zpath = tmp_path / "jgi_img_pr_20261118.zip"
    with zipfile.ZipFile(zpath, "w") as zf:
        zf.writestr("IMGPR_nucl.fna.gz", b"x" * 256)
    digest = hashlib.sha256(zpath.read_bytes()).hexdigest()
    size = zpath.stat().st_size

    sources_ok = {
        "sources": [
            {
                "name": "jgi_img_pr",
                "stored": [
                    {
                        "path": str(zpath),
                        "bytes": size,
                        "sha256": digest,
                    }
                ],
            }
        ]
    }
    assert v.verify_sources_yaml(sources_ok, required_names=["jgi_img_pr"], min_bytes=100) == []

    sources_bad = {
        "sources": [
            {
                "name": "jgi_img_pr",
                "stored": [
                    {
                        "path": str(zpath),
                        "bytes": size + 1,
                        "sha256": digest,
                    }
                ],
            }
        ]
    }
    errs = v.verify_sources_yaml(sources_bad, required_names=["jgi_img_pr"], min_bytes=100)
    assert any("size mismatch" in e for e in errs)

