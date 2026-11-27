"""Card B-75 (DEC-84 (1)): nucleotide ingest refuses a byte that is not a base, and an empty record.

The alphabet is written out here by hand, not taken from `src/`: the IUPAC nucleotide letters `ACGTUNRYKMSWBDHV` of
either case (`U` for the RNA genomes planned). Every fixture's truth (the offending byte, its 1-based position in the
record's sequence, the record) is set by the test when it builds the file. Lower case (soft-masking) is ordinary input
and is stored as read, its uid that of the upper case (MR-14). The reader `biobox.fasta.iter_records` keeps no
alphabet: it also reads protein FASTA, so the check sits where nucleotide sequence enters the store, in pass 1 of
`biobox.l0seq.store` (both `ingest` and the batched `ingest_pass1`). Nothing here needs a tool beyond Python and
pyarrow (GitHub Actions has none of the pod's tools).
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest
from blake3 import blake3

pytestmark = pytest.mark.covers("B-75")  # ring R1 by directory (tests/conftest.py)

LETTERS = b"ACGTUNRYKMSWBDHV"
BASES = LETTERS + LETTERS.lower()


def dna(tag: str, n: int) -> bytes:
    out, h = bytearray(), hashlib.sha256(tag.encode()).digest()
    while len(out) < n:
        out += bytes(b"ACGT"[x & 3] for x in h)
        h = hashlib.sha256(h).digest()
    return bytes(out[:n])


def wrap(seq: bytes, w: int = 60, end: bytes = b"\n") -> bytes:
    return b"".join(seq[i:i + w] + end for i in range(0, len(seq), w))


def fasta_bytes(recs, end: bytes = b"\n") -> bytes:
    return b"".join(b">" + h.encode() + end + wrap(s, end=end) for h, s in recs)


def _obj(p: Path):
    from biobox.l0seq import SourceObject

    return SourceObject(p, "synth", "1", "public", True)


def _files(base: Path) -> list[Path]:
    return sorted(p for p in base.rglob("*") if p.is_file()) if base.exists() else []


def _all_bytes(base: Path) -> dict[str, bytes]:
    return {str(p.relative_to(base)): p.read_bytes() for p in _files(base)}


# --- what is refused -----------------------------------------------------------------------------------------------

# (case, record's sequence as a list of lines, the offending byte, its 1-based position in the joined sequence)
_A = dna("b75-a", 120)
BAD = {
    "gap at the end (GV1's Hoko_1)": ([_A[:60], _A[60:] + b"--"], 0x2D, 121),
    "gap inside (NX1's alignments)": ([_A[:60], _A[60:90] + b"-" + _A[90:]], 0x2D, 91),
    "first byte": ([b"*" + _A[:59], _A[59:]], 0x2A, 1),
    "stop of a protein": ([_A[:60], _A[60:] + b"*"], 0x2A, 121),
    "protein letter X": ([_A[:30] + b"X" + _A[30:60]], 0x58, 31),
    "protein letter E": ([_A[:60], b"E" + _A[60:]], 0x45, 61),
    "dot": ([_A[:10] + b"." + _A[10:]], 0x2E, 11),
    "digit": ([_A[:60], _A[60:70] + b"0" + _A[70:]], 0x30, 71),
    "space inside a line": ([_A[:20] + b" " + _A[20:]], 0x20, 21),
    "space at a line's end": ([_A[:60] + b" ", _A[60:]], 0x20, 61),
    "tab": ([_A[:60], _A[60:61] + b"\t" + _A[61:]], 0x09, 62),
    "stray CR inside a line": ([_A[:40] + b"\r" + _A[40:]], 0x0D, 41),
    "NUL": ([_A[:5] + b"\x00" + _A[5:]], 0x00, 6),
    "byte above ASCII (compressed data, LR2)": ([_A[:60], _A[60:100] + b"\xe5" + _A[100:]], 0xE5, 101),
    "IUPAC-like letter not IUPAC (I)": ([_A[:60], _A[60:] + b"I"], 0x49, 121),
    "lower case not IUPAC (x)": ([_A[:3] + b"x" + _A[3:]], 0x78, 4),
}


def _bad_file(tmp_path: Path, case: str) -> tuple[Path, str, int, int]:
    lines, byte, pos = BAD[case]
    hdr = "bad_rec some description"
    p = tmp_path / "src" / "bad.fna"
    p.parent.mkdir(parents=True)
    data = (b">clean_1\n" + wrap(dna("b75-c1", 150)) + b">" + hdr.encode() + b"\n" + b"".join(ln + b"\n" for ln in lines)
            + b">clean_2\n" + wrap(dna("b75-c2", 150)))
    p.write_bytes(data)
    return p, hdr, byte, pos


def _assert_message(e: Exception, p: Path, hdr: str, byte: int, pos: int) -> None:
    msg = str(e)
    assert str(p) in msg  # the object
    assert hdr in msg  # the record's header
    assert "record 2 " in msg  # its number in the object, 1-based
    assert f"position {pos} " in msg  # 1-based in the record's sequence
    assert f"0x{byte:02x}" in msg  # the byte value
    assert "DEC-84" in msg


@pytest.mark.parametrize("case", list(BAD))
@pytest.mark.parametrize("min_contig_len", [None, 1000], ids=["all-kept", "record-under-min-len"])
def test_b75_ingest_refuses_a_byte_that_is_not_a_base_and_stores_nothing(tmp_path, case, min_contig_len):
    """The second record holds the byte; ingest stops with an error naming the object, the header, the record, the
    1-based position and the byte value, and no file of the stage is in place: not in the root (no table, shard,
    ledger, manifest or store_params) and not in the scratch directory. A record that the least contig length would
    drop is refused the same way (the check reads every record)."""
    from biobox.l0seq import store

    p, hdr, byte, pos = _bad_file(tmp_path, case)
    root, scratch = tmp_path / "root", tmp_path / "scratch"
    scratch.mkdir()
    with pytest.raises(store.NotABase) as ei:
        store.ingest([_obj(p)], root, "S0", tmp=scratch, min_contig_len=min_contig_len)
    assert isinstance(ei.value, ValueError)  # `bio-box stage ingest` prints a ValueError after `refused:`, exit 2
    _assert_message(ei.value, p, hdr, byte, pos)
    assert _files(root) == []
    assert store.stage_files(root, "S0") == []
    assert _files(scratch) == []


@pytest.mark.parametrize("case", ["gap at the end (GV1's Hoko_1)", "byte above ASCII (compressed data, LR2)"])
def test_b75_batched_pass1_refuses_the_same_and_leaves_no_run(tmp_path, case):
    """`build-stage` ingests a large stage in batch processes (card B-30): pass 1 of a shard refuses the same record,
    and neither its run (scratch) nor its record (the root) is written."""
    from biobox.l0seq import store

    p, hdr, byte, pos = _bad_file(tmp_path, case)
    root, scratch, runs = tmp_path / "root", tmp_path / "scratch", tmp_path / "runs"
    scratch.mkdir()
    runs.mkdir()
    with pytest.raises(store.NotABase) as ei:
        store.ingest_pass1([_obj(p)], root, "S0", 0, 1 << 30, tmp=scratch, runs_tmp=runs)
    _assert_message(ei.value, p, hdr, byte, pos)
    assert _files(root) == []
    assert _files(runs) == []
    assert _files(scratch) == []


@pytest.mark.parametrize("where", ["between", "last", "only"])
def test_b75_ingest_refuses_an_empty_record(tmp_path, where):
    """A header with no sequence line after it is refused, naming the object, the header and the record."""
    from biobox.l0seq import store

    recs = {"between": [("a", dna("b75-e1", 90)), ("empty_one x", b""), ("c", dna("b75-e2", 90))],
            "last": [("a", dna("b75-e1", 90)), ("c", dna("b75-e2", 90)), ("empty_one x", b"")],
            "only": [("empty_one x", b"")]}[where]
    n = [h for h, _ in recs].index("empty_one x") + 1
    p = tmp_path / "e.fna"
    p.write_bytes(fasta_bytes(recs))
    root = tmp_path / "root"
    with pytest.raises(store.NotABase) as ei:
        store.ingest([_obj(p)], root, "S0")
    msg = str(ei.value)
    assert str(p) in msg and "empty_one x" in msg and f"record {n} " in msg and "no sequence" in msg
    assert _files(root) == []


def test_b75_the_cli_reports_the_refusal_and_exits_2(tmp_path, capsys, monkeypatch):
    """What `build-stage` and `run-step` show of it: the step's process prints `refused: <the error>` on stderr and
    exits 2 (run-step then records the failed attempt and publishes nothing)."""
    from biobox import cli

    monkeypatch.delenv("BIOBOX_STAGING", raising=False)
    p, hdr, byte, pos = _bad_file(tmp_path, "gap at the end (GV1's Hoko_1)")
    objs = tmp_path / "objects.json"
    objs.write_text(json.dumps([{"path": str(p), "source_name": "synth", "source_release": "1",
                                 "license": "public", "redistributable": True}]))
    root = tmp_path / "root"
    rc = cli.main(["stage", "ingest", "--root", str(root), "--stage", "S0", "--objects", str(objs),
                   "--tmp", str(tmp_path)])
    err = capsys.readouterr().err
    assert rc == 2
    assert err.startswith("refused: ")
    assert str(p) in err and hdr in err and f"position {pos} " in err and "0x2d" in err
    assert _files(root) == []


# --- what passes, unchanged -----------------------------------------------------------------------------------------


def test_b75_every_iupac_letter_of_either_case_passes_and_is_stored_as_read(tmp_path):
    """Every letter of the alphabet in both cases, `U` and `u` included, soft-masked runs and `N` runs, with LF and
    with CRLF line ends: ingest takes them, stores the bytes as read (lower case kept, MR-14), and the uid is that of
    the upper case."""
    from biobox import ids
    from biobox.l0seq import store

    every = (LETTERS + LETTERS.lower()) * 7
    soft = dna("b75-s1", 100) + dna("b75-s2", 80).lower() + b"N" * 50 + dna("b75-s3", 70)
    rna = dna("b75-r", 120).replace(b"T", b"U")
    rna_soft = rna[:60] + rna[60:].lower()
    recs = [("every_letter", every), ("soft_masked", soft), ("rna", rna), ("rna_soft", rna_soft)]
    for end, name in ((b"\n", "lf.fna"), (b"\r\n", "crlf.fna")):
        p = tmp_path / name
        p.write_bytes(fasta_bytes(recs, end=end))
        root = tmp_path / f"root_{name}"
        store.ingest([_obj(p)], root, "S0")
        s = store.SeqStore(root, ["S0"])
        try:
            for _, seq in recs:
                u = ids.contig_uid(seq)
                assert u == blake3(seq.upper()).digest(length=16)
                assert s.raw(ids.raw_digest(seq)) == seq  # the bytes as read, case kept
                assert s.length_of(s.default_digest(u)) == len(seq)
        finally:
            s.close()


def test_b75_line_ends_crlf_and_trailing_crs_are_stripped_as_before(tmp_path):
    """`iter_records` strips every trailing CR and LF of a line (`rstrip(b"\\r\\n")`), so CRLF, a doubled CR before
    the LF and a CR on a last line without LF never reach the sequence and are not refused; a CR inside a line would
    be stored, and is refused (`BAD`, "stray CR inside a line")."""
    from biobox import fasta, ids
    from biobox.l0seq import store

    a, b = dna("b75-l1", 100), dna("b75-l2", 70)
    p = tmp_path / "cr.fna"
    p.write_bytes(b">a\r\n" + a[:60] + b"\r\r\n" + a[60:] + b"\r\n\r\n>b\n" + b[:60] + b"\n" + b[60:] + b"\r")
    assert [r.seq for r in fasta.iter_records(p)] == [a, b]  # the reader as it was
    root = tmp_path / "root"
    store.ingest([_obj(p)], root, "S0")
    s = store.SeqStore(root, ["S0"])
    try:
        assert s.get(ids.contig_uid(a)) == a and s.get(ids.contig_uid(b)) == b
    finally:
        s.close()


def test_b75_the_reader_keeps_no_alphabet(tmp_path):
    """`biobox.fasta.iter_records` is unchanged: it reads a protein FASTA (stops, `X`) and gapped records as they are;
    the check belongs to nucleotide ingest only."""
    from biobox import fasta

    prot = b"MKVLAAGIVEEX*"
    p = tmp_path / "prot.faa"
    p.write_bytes(b">p1 protein\n" + prot + b"\n>p2\nAC--GT\n>p3\n")
    assert [(r.header, r.seq) for r in fasta.iter_records(p)] == [("p1 protein", prot), ("p2", b"AC--GT"), ("p3", b"")]


def test_b75_a_clean_stage_is_byte_identical_with_and_without_the_check(tmp_path, monkeypatch):
    """The check adds no output: the same clean objects ingested with it and with it replaced by the bare reader give
    the same files, byte for byte (two objects, soft-masked and CRLF records, a repeat, the least contig length)."""
    from biobox import fasta
    from biobox.l0seq import store

    src = tmp_path / "src"
    src.mkdir()
    one, two = src / "one.fna", src / "two.fna.gz"
    shared = dna("b75-shared", 400)
    one.write_bytes(fasta_bytes([("x1 a", dna("b75-x1", 900)), ("x2", dna("b75-x2", 500).lower()),
                                 ("x3", shared), ("tiny", dna("b75-t", 30))]))
    import gzip

    two.write_bytes(gzip.compress(fasta_bytes([("y1", dna("b75-y1", 700)), ("y2", shared)], end=b"\r\n"), mtime=0))
    objs = [_obj(one), _obj(two)]
    for mcl in (None, 100):
        with_check = tmp_path / f"with_{mcl}"
        store.ingest(objs, with_check, "S0", min_contig_len=mcl)
        with monkeypatch.context() as m:
            m.setattr(store, "_nucleotide_records", fasta.iter_records)
            without = tmp_path / f"without_{mcl}"
            store.ingest(objs, without, "S0", min_contig_len=mcl)
        got, want = _all_bytes(with_check), _all_bytes(without)
        assert got and got.keys() == want.keys()
        assert got == want
