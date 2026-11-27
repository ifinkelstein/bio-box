"""Card R-21: `s2_check_inputs.py` reads the bytes of every object it checks.

`LR2` was staged with 691 gzip files written into its FASTA as sequence and passed the check twice, because the check
compared tables and checksums and never opened an object (card R-19, its note §5). The byte pass reads each object of
`objects.json` once, streaming, with card R-19's reader, and fails with a named problem for: a header that is not
UTF-8; a byte outside the IUPAC nucleotide letters (gap bytes `-` and `.` named apart from the others); an empty
record; a record whose `contig_uid` is in no genome of `genomes.json`; a contig of `genomes.json` that no object holds;
and, so that the check agrees with ingest, a line before the first header, a read error and a file whose compression
is not what its name says (ingest opens by the name). `--tables-only` skips the pass.

Truth is written by hand here: each fixture's contig uids are blake3-128 of the upper-cased sequence computed in this
file (the definition of `biobox.ids.contig_uid`, not a call of the check), and the agreement with ingest is checked
against `biobox.fasta.iter_records`, the reader ingest uses. Each test runs twice, with the normal PATH and with PATH
reduced to `/usr/bin:/bin`: GitHub Actions has none of the pod's tools.
"""

from __future__ import annotations

import gzip
import importlib.util
import json
import random
import subprocess
import sys
import tracemalloc
from pathlib import Path
from types import SimpleNamespace

import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("R-21")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
SCRIPT = BSP / "s2_check_inputs.py"


@pytest.fixture(params=["normal-path", "tools-hidden"], autouse=True)
def path_mode(request, monkeypatch):
    if request.param == "tools-hidden":
        monkeypatch.setenv("PATH", "/usr/bin:/bin")
    return request.param


def _uid(seq: bytes) -> str:
    return blake3(seq.upper()).digest(length=16).hex()


def _gz(data: bytes) -> bytes:
    return gzip.compress(data, mtime=0)


def _fasta(records: list[tuple[bytes, bytes]], width: int = 0, eol: bytes = b"\n") -> bytes:
    out = bytearray()
    for hdr, seq in records:
        out += b">" + hdr + eol
        if width:
            for i in range(0, len(seq), width):
                out += seq[i:i + width] + eol
        elif seq:
            out += seq + eol
    return bytes(out)


def _stage(tmp: Path, objects: list[tuple[str, bytes]], genomes: list[list[bytes]] | None = None,
           extra_uids: list[str] | None = None) -> tuple[Path, Path]:
    """Write `objects` (file name, file bytes) and a genomes.json; by default one genome per record, its uid made
    from the record's sequence as written, so the tables agree with the objects whatever the objects hold."""
    tmp.mkdir(parents=True, exist_ok=True)
    olist = []
    for name, data in objects:
        p = tmp / name
        p.write_bytes(data)
        olist.append({"path": str(p), "source_name": "r21_fixture", "source_release": "1",
                      "license": "test", "redistributable": False,
                      "retrieved_at": "2026-11-27T00:00:00+00:00"})
    if genomes is None:
        raise ValueError("give the genomes")
    glist = [{"contig_uids": [_uid(s) for s in g], "kind": "prokaryote", "genetic_code": 11,
              "assembly_quality": "mag_hq"} for g in genomes]
    if extra_uids:
        glist.append({"contig_uids": extra_uids, "kind": "prokaryote", "genetic_code": 11,
                      "assembly_quality": "mag_hq"})
    gp, op = tmp / "genomes.json", tmp / "objects.json"
    gp.write_text(json.dumps(glist), encoding="utf-8")
    op.write_text(json.dumps(olist), encoding="utf-8")
    return gp, op


def _check(gp: Path, op: Path, *opts: str, stage: str = "TST1") -> tuple[int, dict]:
    p = subprocess.run([sys.executable, str(SCRIPT), *opts, stage, str(gp), str(op)],
                       capture_output=True, text=True, check=False)
    assert p.stdout.strip().startswith("{"), p.stdout + p.stderr
    return p.returncode, json.loads(p.stdout)


def _named(doc: dict, name: str) -> list[str]:
    return [m for m in doc["problems"] if f": {name}:" in m]


A = b"ACGTACGTACGTNNNNacgtRYSWKMBDHVU"  # 31 IUPAC letters (12 + 4 + 4 + 11), both cases
B = b"GGGGCCCCAAAATTTTggggccccaaaatttt"  # 32
C = b"ACGTU" * 7  # 35


def test_a_clean_stage_passes_and_its_uids_agree_with_ingest(tmp_path):
    plain = _fasta([(b"c1 first", A), (b"c2", B)], width=7, eol=b"\r\n")
    gz = _gz(_fasta([(b"c3", C), (b"caf\xc3\xa9 en\xe2\x80\x93dash", A.lower())], width=10) + b"\n\n")
    gp, op = _stage(tmp_path, [("a.fna", plain), ("b.fna.gz", gz)], genomes=[[A, B], [C]])
    rc, doc = _check(gp, op)
    assert rc == 0, doc["problems"]
    assert doc["problems"] == []
    bp = doc["byte_pass"]
    # by hand: 2 + 2 records; 31 + 32 and 35 + 31 bases; A and A.lower() are one contig
    assert [o["records"] for o in bp["objects"]] == [2, 2]
    assert [o["bases"] for o in bp["objects"]] == [63, 66]
    assert [o["compression"] for o in bp["objects"]] == ["none", "gzip"]
    assert bp["records"] == 4 and bp["bases"] == 129
    assert bp["contigs"] == 3 and bp["contigs_in_no_object"] == 0
    # a UTF-8 header that is not ASCII is accepted by ingest: counted, not a problem
    assert bp["objects"][1]["headers_not_ascii"] == 1
    # the oracle: ingest's own reader and contig uid over the same files
    from biobox import fasta, ids

    got = {ids.contig_uid(r.seq).hex() for o in ("a.fna", "b.fna.gz") for r in fasta.iter_records(tmp_path / o)}
    assert got == {_uid(A), _uid(B), _uid(C)}


def test_the_lr2_damage_passed_the_tables_and_fails_the_byte_pass(tmp_path):
    # LR2's case (card R-19 §4.2): a gzip member read as text and written as "contigs", the tables made from what
    # was written, so genomes.json agrees with the FASTA. The records after a '>' in compressed bytes.
    blob = _gz(random.Random(21).randbytes(20000))
    junk = [seg for seg in blob.split(b">")[1:] if seg]
    assert junk, "the fixture needs a '>' in the compressed bytes"
    recs = [(b"ok1", A)]
    for k, seg in enumerate(junk):
        recs.append((b"g%d" % k, seg.replace(b"\n", b"").replace(b"\r", b"")))
    recs.append((b"bad\xff\xfe", B))
    data = _gz(_fasta([(h, s) for h, s in recs], width=80))
    gp, op = _stage(tmp_path, [("lr2.fna.gz", data)], genomes=[[s for _, s in recs if s]])
    rc_t, doc_t = _check(gp, op, "--tables-only")
    assert rc_t == 0 and doc_t["problems"] == [] and doc_t["byte_pass"] is None  # what the check said before R-21
    rc, doc = _check(gp, op)
    assert rc == 1
    assert len(_named(doc, "header_not_utf8")) == 1
    assert f"first record {len(recs)} " in _named(doc, "header_not_utf8")[0]
    assert _named(doc, "non_iupac_byte"), doc["problems"]
    o = doc["byte_pass"]["objects"][0]
    assert o["headers_not_utf8"] == 1
    assert o["records"] == len(recs)


def test_a_header_that_is_not_utf8_is_named_with_its_record(tmp_path):
    data = _fasta([(b"ok", A), (b"bad\xff\xfeheader", B), (b"ok3", C)])
    gp, op = _stage(tmp_path, [("x.fna", data)], genomes=[[A, B, C]])
    rc, doc = _check(gp, op)
    assert rc == 1
    msg = _named(doc, "header_not_utf8")
    assert len(msg) == 1 and "1 record" in msg[0] and "first record 2" in msg[0]
    assert doc["byte_pass"]["objects"][0]["headers_not_utf8"] == 1
    assert len(doc["problems"]) == 1


def test_gap_bytes_are_named_apart_from_other_bytes_outside_iupac(tmp_path):
    gapped = b"ACGT--AC.GT"      # 2 '-', 1 '.'
    gapped2 = b"AC-GT"           # 1 '-'
    other = b"ACGTXACGT*\xe5"   # 3 other bytes
    both = b"AC-GTZ"             # 1 '-', 1 other
    data = _fasta([(b"r1", A), (b"r2", gapped), (b"r3", gapped2), (b"r4", other), (b"r5", both)])
    gp, op = _stage(tmp_path, [("x.fna", data)], genomes=[[A, gapped, gapped2, other, both]])
    rc, doc = _check(gp, op)
    assert rc == 1
    gap = _named(doc, "gap_byte")
    oth = _named(doc, "non_iupac_byte")
    assert len(gap) == 1 and len(oth) == 1, doc["problems"]
    assert "3 record" in gap[0] and "4 '-'" in gap[0] and "1 '.'" in gap[0] and "first record 2" in gap[0]
    assert "2 record" in oth[0] and "4 byte" in oth[0] and "first record 4" in oth[0]
    o = doc["byte_pass"]["objects"][0]
    assert o["records_with_gap_bytes"] == 3 and o["gap_bytes"] == {"-": 4, ".": 1}
    assert o["records_with_other_non_iupac"] == 2 and o["other_non_iupac_bytes"] == 4
    assert len(doc["problems"]) == 2


def test_an_empty_record_is_named(tmp_path):
    data = b">r1\n" + A + b"\n>empty\n>r3\n" + B + b"\n"
    gp, op = _stage(tmp_path, [("x.fna", data)], genomes=[[A, B]], extra_uids=[_uid(b"")])
    rc, doc = _check(gp, op)
    assert rc == 1
    msg = _named(doc, "empty_record")
    assert len(msg) == 1 and "first record 2" in msg[0], doc["problems"]
    assert len(doc["problems"]) == 1


def test_a_record_in_no_genome_and_a_contig_in_no_object_are_named(tmp_path):
    data = _fasta([(b"r1", A), (b"stray", C), (b"r3", B)])
    gp, op = _stage(tmp_path, [("x.fna", data)], genomes=[[A, B]], extra_uids=[_uid(b"TTTTTTTTTT")])
    rc, doc = _check(gp, op)
    assert rc == 1
    stray = _named(doc, "record_not_in_genomes")
    assert len(stray) == 1 and "1 record" in stray[0] and "first record 2" in stray[0] and "stray" in stray[0]
    missing = [m for m in doc["problems"] if m.startswith("contigs_in_no_object:")]
    assert len(missing) == 1 and "1 contig" in missing[0] and _uid(b"TTTTTTTTTT") in missing[0]
    assert doc["byte_pass"]["contigs_in_no_object"] == 1
    assert doc["byte_pass"]["objects"][0]["records_not_in_genomes"] == 1
    assert len(doc["problems"]) == 2


def test_a_repeated_record_and_case_variants_are_not_problems(tmp_path):
    data = _fasta([(b"r1", A), (b"r1again", A), (b"r1lower", A.lower())])
    gp, op = _stage(tmp_path, [("x.fna", data)], genomes=[[A]])
    rc, doc = _check(gp, op)
    assert rc == 0 and doc["problems"] == []


def test_a_contig_is_found_in_any_object(tmp_path):
    gp, op = _stage(tmp_path, [("a.fna", _fasta([(b"a", A)])), ("b.fna", _fasta([(b"b", B)]))], genomes=[[A, B]])
    rc, doc = _check(gp, op)
    assert rc == 0 and doc["problems"] == []


@pytest.mark.parametrize("case", ["gzip_named_plain", "plain_named_gz"])
def test_a_file_whose_compression_is_not_its_name_is_named(tmp_path, case):
    plain = _fasta([(b"r1", A)])
    name, data = ("x.fna", _gz(plain)) if case == "gzip_named_plain" else ("x.fna.gz", plain)
    gp, op = _stage(tmp_path, [(name, data)], genomes=[[A]])
    rc, doc = _check(gp, op)
    assert rc == 1
    msg = _named(doc, "compression_not_as_named")
    assert len(msg) == 1 and "ingest" in msg[0], doc["problems"]


def test_a_line_before_the_first_header_and_a_truncated_gzip_are_named(tmp_path):
    gp, op = _stage(tmp_path / "before", [("x.fna", b"ACGT\n" + _fasta([(b"r1", A)]))], genomes=[[A]])
    rc, doc = _check(gp, op)
    assert rc == 1 and len(_named(doc, "sequence_before_first_header")) == 1, doc["problems"]
    whole = _gz(_fasta([(b"r1", A), (b"r2", B)] * 50))
    gp, op = _stage(tmp_path / "trunc", [("x.fna.gz", whole[: len(whole) // 2])], genomes=[[A, B]])
    rc, doc = _check(gp, op)
    assert rc == 1 and len(_named(doc, "read_error")) == 1, doc["problems"]
    assert doc["byte_pass"]["objects"][0]["read_error"]
    # an object with no record at all (an empty file), beside one that holds the contig
    gp, op = _stage(tmp_path / "none", [("x.fna", _fasta([(b"r1", A)])), ("empty.fna", b"")], genomes=[[A]])
    rc, doc = _check(gp, op)
    assert rc == 1 and len(_named(doc, "no_records")) == 1 and len(doc["problems"]) == 1, doc["problems"]


def test_tables_only_skips_the_byte_pass_and_keeps_every_table_check(tmp_path):
    data = _fasta([(b"r1", b"AC-GT")])
    gp, op = _stage(tmp_path, [("x.fna", data)], genomes=[[b"AC-GT"]])
    rc, doc = _check(gp, op, "--tables-only")
    assert rc == 0 and doc["problems"] == [] and doc["byte_pass"] is None
    rc2, doc2 = _check(gp, op)
    assert rc2 == 1 and _named(doc2, "gap_byte")
    # a table problem is still found with the option (an enum value outside the schema)
    g = json.loads(gp.read_text())
    g[0]["assembly_quality"] = "not_a_quality"
    gp.write_text(json.dumps(g))
    rc3, doc3 = _check(gp, op, "--tables-only")
    assert rc3 == 1 and any("assembly_quality" in m for m in doc3["problems"])


def test_usage_errors_exit_2(tmp_path):
    gp, op = _stage(tmp_path, [("x.fna", _fasta([(b"r1", A)]))], genomes=[[A]])
    for argv in ([str(SCRIPT), "TST1", str(gp)], [str(SCRIPT), "--bytes-too", "TST1", str(gp), str(op)]):
        p = subprocess.run([sys.executable, *argv], capture_output=True, text=True, check=False)
        assert p.returncode == 2, p.stdout + p.stderr


def _load():
    spec = importlib.util.spec_from_file_location("s2_check_inputs", SCRIPT)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["s2_check_inputs"] = mod
    spec.loader.exec_module(mod)
    return mod


def test_the_byte_pass_holds_no_whole_sequence(tmp_path):
    # one unwrapped record of 16 MiB in a gzip file: the pass may hold a line piece and a hash buffer, not the record
    big = (b"ACGTN" * (16 * 1024 * 1024 // 5 + 1))[: 16 * 1024 * 1024]
    p = tmp_path / "big.fna.gz"
    with gzip.open(p, "wb", compresslevel=1) as fh:
        fh.write(b">big\n" + big + b"\n>small\n" + A + b"\n")
    mod = _load()
    from biobox.l0seq.inputs import GenomeInput

    genomes = [GenomeInput(contig_uids=(bytes.fromhex(_uid(big)), bytes.fromhex(_uid(A))))]
    objs = [SimpleNamespace(path=p)]
    tracemalloc.start()
    try:
        rep, problems = mod.byte_pass(genomes, objs)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert problems == [], problems
    assert rep["records"] == 2 and rep["bases"] == len(big) + len(A)
    assert peak < 8 * 1024 * 1024, peak
