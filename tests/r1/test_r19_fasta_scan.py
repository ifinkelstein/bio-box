"""Card R-19: a streaming scan of a staged FASTA (plain or gzip) for bytes that are not sequence.

`LR2`'s staged FASTA was found holding compressed data written into it (card C-81): headers that are not UTF-8 and
sequence bytes that are not IUPAC letters. The scan reads one file, one bounded line piece at a time, and reports the
records, the headers that are not UTF-8 or not printable ASCII, the records holding a byte outside the IUPAC nucleotide
letters, the empty records and the first offending record. Every fixture's truth is written out by hand in this file
(constructed truth, not a snapshot of the script's output). Each test runs twice, with the normal PATH and with PATH
reduced to `/usr/bin:/bin`: the scan needs nothing beyond the Python standard library, and GitHub Actions has none of
the pod's tools.
"""

from __future__ import annotations

import gzip
import importlib.util
import json
import os
import subprocess
import sys
import tracemalloc
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("R-19")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "stress" / "build_stage_pod" / "r19_fasta_scan.py"


def _load():
    spec = importlib.util.spec_from_file_location("r19_fasta_scan", SCRIPT)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["r19_fasta_scan"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(params=["normal-path", "tools-hidden"], autouse=True)
def path_mode(request, monkeypatch):
    if request.param == "tools-hidden":
        monkeypatch.setenv("PATH", "/usr/bin:/bin")
    return request.param


@pytest.fixture
def scan():
    return _load().scan_file


def _gz(data: bytes) -> bytes:
    return gzip.compress(data, mtime=0)


CLEAN = (
    b">c1 first contig\n"
    b"ACGTACGTAC\n"
    b"GTNNacgtRY\n"
    b">c2\n"
    b"SWKMBDHVUu\n"
)


def _counts(rep: dict) -> tuple:
    return (
        rep["records"],
        rep["bases"],
        rep["headers_not_utf8"],
        rep["headers_not_printable_ascii"],
        rep["records_with_non_iupac"],
        rep["empty_records"],
    )


def test_a_clean_plain_fasta_scans_clean(tmp_path, scan):
    p = tmp_path / "clean.fna"
    p.write_bytes(CLEAN)
    rep = scan(p)
    assert rep["status"] == "clean"
    assert rep["compression"] == "none"
    assert _counts(rep) == (2, 30, 0, 0, 0, 0)
    assert rep["first_offending"] is None
    assert rep["read_error"] is None
    assert rep["file_bytes"] == len(CLEAN)


def test_a_gzip_file_is_read_by_its_magic_and_gives_the_same_counts(tmp_path, scan):
    p = tmp_path / "clean.fna.gz"
    p.write_bytes(_gz(CLEAN))
    rep = scan(p)
    assert rep["status"] == "clean"
    assert rep["compression"] == "gzip"
    assert _counts(rep) == (2, 30, 0, 0, 0, 0)
    # A gzip file whose name does not say so is still read as gzip.
    q = tmp_path / "named_plain.fna"
    q.write_bytes(_gz(CLEAN))
    rep2 = scan(q)
    assert rep2["compression"] == "gzip"
    assert _counts(rep2) == (2, 30, 0, 0, 0, 0)


def test_crlf_and_blank_lines_are_read_as_biobox_fasta_reads_them(tmp_path, scan):
    p = tmp_path / "crlf.fna"
    p.write_bytes(b">c1\r\nACGT\r\n\r\nACGT\r\n\n>c2\r\nGG\r\n")
    rep = scan(p)
    assert rep["status"] == "clean"
    assert _counts(rep) == (2, 10, 0, 0, 0, 0)


def test_a_header_that_is_not_utf8_is_counted_and_is_the_first_offending_record(tmp_path, scan):
    p = tmp_path / "bad_header.fna.gz"
    p.write_bytes(_gz(b">ok\nACGT\n>bad\xff\xfeheader\nACGT\n>caf\xc3\xa9\nACGT\n"))
    rep = scan(p)
    assert rep["status"] == "damaged"
    # record 2: not UTF-8 (so not printable ASCII either); record 3: UTF-8 but not ASCII.
    assert _counts(rep) == (3, 12, 1, 2, 0, 0)
    first = rep["first_offending"]
    assert first["record"] == 2
    assert first["header"] == "bad\\xff\\xfeheader"
    assert first["reasons"] == ["header_not_utf8", "header_not_printable_ascii"]
    assert rep["first_record_of"]["header_not_utf8"] == 2
    assert rep["first_record_of"]["header_not_printable_ascii"] == 2


def test_a_tab_in_a_header_is_printable(tmp_path, scan):
    p = tmp_path / "tab.fna"
    p.write_bytes(b">c1\tlen=4\nACGT\n")
    rep = scan(p)
    assert rep["status"] == "clean"
    assert rep["headers_not_printable_ascii"] == 0


def test_a_byte_outside_the_iupac_letters_is_counted_with_its_value_and_offset(tmp_path, scan):
    p = tmp_path / "bad_seq.fna"
    # record 1 clean; record 2 holds 0xe5 at sequence offset 6; record 3 a gap '-' and a '*'; record 4 a digit.
    p.write_bytes(b">a\nACGT\n>b\nACGTAC\xe5GT\n>c\nAC-GT*\n>d\nAC1\n")
    rep = scan(p)
    assert rep["status"] == "damaged"
    assert _counts(rep) == (4, 4 + 9 + 6 + 3, 0, 0, 3, 0)
    first = rep["first_offending"]
    assert first["record"] == 2
    assert first["header"] == "b"
    assert first["reasons"] == ["non_iupac_byte"]
    assert first["first_bad_byte"] == "0xe5"
    assert first["first_bad_offset"] == 6
    assert rep["non_iupac_bytes"] == 4  # 0xe5, '-', '*', '1'
    assert rep["non_iupac_byte_values"] == {"0x2a": 1, "0x2d": 1, "0x31": 1, "0xe5": 1}
    assert rep["distinct_non_iupac_byte_values"] == 4
    q = tmp_path / "gaps.fna"
    q.write_bytes(b">a\nAC--GT\n>b\nA-C\n")
    assert scan(q)["non_iupac_byte_values"] == {"0x2d": 3}


def test_an_empty_record_is_counted(tmp_path, scan):
    p = tmp_path / "empty.fna"
    p.write_bytes(b">a\nACGT\n>b\n>c\nGG\n>d\n")
    rep = scan(p)
    assert rep["status"] == "damaged"
    assert _counts(rep) == (4, 6, 0, 0, 0, 2)
    assert rep["first_offending"]["record"] == 2
    assert rep["first_offending"]["reasons"] == ["empty_record"]
    assert rep["first_record_of"]["empty_record"] == 2


def test_sequence_before_the_first_header_is_reported(tmp_path, scan):
    p = tmp_path / "nohdr.fna"
    p.write_bytes(b"ACGT\n>a\nACGT\n")
    rep = scan(p)
    assert rep["status"] == "damaged"
    assert rep["lines_before_first_header"] == 1
    assert rep["records"] == 1
    assert rep["first_offending"]["record"] == 0
    assert rep["first_offending"]["reasons"] == ["sequence_before_first_header"]


def test_compressed_bytes_written_into_a_fasta_are_found_where_they_start(tmp_path, scan):
    # What card C-81 found in LR2: good records, then compressed data written as if it were sequence. The writer
    # wraps at 80, so a '>' that lands at the start of a line opens a record of its own.
    good = b"".join(b">%032x\n" % i + b"ACGT" * 30 + b"\n" for i in range(5))
    blob = _gz(bytes(range(256)) * 64 + b"ACGT" * 1000)
    wrapped = b"".join(blob[i : i + 80] + b"\n" for i in range(0, len(blob), 80))
    data = good + b">%032x\n" % 99 + wrapped + b">%032x\n" % 100 + b"ACGT\n"
    p = tmp_path / "lr2_like.fna.gz"
    p.write_bytes(_gz(data))
    rep = scan(p)
    assert rep["status"] == "damaged"
    assert rep["records"] >= 7
    assert rep["records_with_non_iupac"] >= 1
    first = rep["first_offending"]
    assert first["record"] == 6
    assert first["header"] == "%032x" % 99
    assert "non_iupac_byte" in first["reasons"]


def test_a_truncated_gzip_is_a_read_error(tmp_path, scan):
    blob = _gz(CLEAN * 200)
    p = tmp_path / "cut.fna.gz"
    p.write_bytes(blob[: len(blob) - 20])
    rep = scan(p)
    assert rep["status"] == "damaged"
    assert rep["read_error"] is not None
    assert "EOFError" in rep["read_error"] or "BadGzipFile" in rep["read_error"] or "error" in rep["read_error"]


def test_trailing_bytes_that_are_not_gzip_after_the_stream_are_a_read_error(tmp_path, scan):
    p = tmp_path / "trail.fna.gz"
    p.write_bytes(_gz(CLEAN) + b"this is not gzip at all")
    rep = scan(p)
    assert rep["status"] == "damaged"
    assert rep["read_error"] is not None


def test_concatenated_gzip_members_are_read_whole(tmp_path, scan):
    p = tmp_path / "multi.fna.gz"
    p.write_bytes(_gz(CLEAN) + _gz(b">c3\nACGT\n"))
    rep = scan(p)
    assert rep["status"] == "clean"
    assert _counts(rep) == (3, 34, 0, 0, 0, 0)


def test_a_long_unwrapped_record_is_scanned_in_bounded_memory(tmp_path, scan):
    seq = b"ACGT" * (4 * 1024 * 1024)  # one line of 16 MiB
    p = tmp_path / "long.fna"
    p.write_bytes(b">long\n" + seq + b"X\n>b\nAC\n")
    tracemalloc.start()
    rep = scan(p)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert rep["records"] == 2
    assert rep["bases"] == len(seq) + 1 + 2
    assert rep["records_with_non_iupac"] == 1
    assert rep["first_offending"]["first_bad_offset"] == len(seq)
    assert peak < 8 * 1024 * 1024, peak  # the 16 MiB line is never held whole


def test_the_sha256_of_the_file_is_reported(tmp_path, scan):
    import hashlib

    p = tmp_path / "clean.fna.gz"
    blob = _gz(CLEAN)
    p.write_bytes(blob)
    assert scan(p)["sha256"] == hashlib.sha256(blob).hexdigest()


def test_objects_json_paths_are_listed_absolute(tmp_path):
    mod = _load()
    d = tmp_path / "STAGE"
    d.mkdir()
    (d / "objects.json").write_text(
        json.dumps([{"path": str(tmp_path / "a.fna.gz")}, {"path": "b.fna"}]), encoding="utf-8"
    )
    assert mod.object_paths(d / "objects.json") == [tmp_path / "a.fna.gz", d / "b.fna"]


def test_the_command_line_writes_one_json_atomically_and_exits_by_status(tmp_path):
    good = tmp_path / "good.fna"
    good.write_bytes(CLEAN)
    bad = tmp_path / "bad.fna"
    bad.write_bytes(b">a\nAC\xe5\n")
    out_good = tmp_path / "out" / "good.json"
    out_bad = tmp_path / "out" / "bad.json"
    env = dict(os.environ)
    r1 = subprocess.run(
        [sys.executable, "-I", str(SCRIPT), "scan", str(good), str(out_good), "--stage", "T1"],
        capture_output=True, text=True, env=env,
    )
    assert r1.returncode == 0, r1.stderr
    r2 = subprocess.run(
        [sys.executable, "-I", str(SCRIPT), "scan", str(bad), str(out_bad), "--stage", "T1"],
        capture_output=True, text=True, env=env,
    )
    assert r2.returncode == 1, r2.stderr
    g = json.loads(out_good.read_text(encoding="utf-8"))
    b = json.loads(out_bad.read_text(encoding="utf-8"))
    assert g["status"] == "clean" and g["stage"] == "T1" and g["path"] == str(good)
    assert b["status"] == "damaged" and b["records_with_non_iupac"] == 1
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == ["bad.json", "good.json"]  # no temporary left


def test_scan_stage_scans_every_object_in_its_own_process_and_resumes(tmp_path):
    d = tmp_path / "S"
    d.mkdir()
    (d / "a.fna").write_bytes(CLEAN)
    (d / "b.fna.gz").write_bytes(_gz(b">x\nAC\xe5\n"))
    (d / "objects.json").write_text(
        json.dumps([{"path": str(d / "a.fna")}, {"path": str(d / "b.fna.gz")}]), encoding="utf-8"
    )
    out = tmp_path / "scan"
    cmd = [sys.executable, "-I", str(SCRIPT), "scan-stage", str(d / "objects.json"), str(out), "--stage", "S"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    summary = json.loads((out / "S.stage.json").read_text(encoding="utf-8"))
    assert summary["stage"] == "S"
    assert summary["files"] == 2
    assert summary["damaged"] == ["b.fna.gz"]
    assert [f["file"] for f in summary["per_file"]] == ["a.fna", "b.fna.gz"]
    assert summary["per_file"][0]["records"] == 2
    a_json = out / "S__a.fna.json"
    before = a_json.stat().st_mtime_ns
    r2 = subprocess.run(cmd, capture_output=True, text=True)
    assert r2.returncode == 0, r2.stderr
    assert a_json.stat().st_mtime_ns == before  # a finished file is not scanned again
    t = subprocess.run([sys.executable, "-I", str(SCRIPT), "table", str(out)], capture_output=True, text=True)
    assert t.returncode == 0, t.stderr
    rows = [ln for ln in t.stdout.splitlines() if ln.startswith("| `S`")]
    assert rows == [
        "| `S` | `a.fna` | clean | 2 | 30 | 0 | 0 | 0 | 0 | 0 | 0 | none |",
        "| `S` | `b.fna.gz` | damaged | 1 | 3 | 0 | 0 | 1 | 1 | 0 | 0 | none |",
    ]
