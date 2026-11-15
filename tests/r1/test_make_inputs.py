"""Card C-21a (DEC-27): `scripts/stress/build_stage_pod/make_inputs.py` lists every record in its genome.

Before DEC-27 the script did not list a contig uid already listed by an earlier genome of the stage (6,818 of 34,076
INPHARED records, *measured*, `planning/measurements/2026-11-10_rebuild_10gbp.md`), which changed the later genome's
uid and its training set. Now a contig is listed in every genome that holds it; a genome whose contig set repeats an
earlier genome's is the same genome (`genome_uid` hashes the contig set) and is listed once; a contig repeated inside
one genome is listed once in it. `--drop-repeats` keeps the old rule byte for byte, because the inputs of stages S1,
P1 to P6 and H1 were made with it. The oracle for the old bytes is the rule written out here, not the script.
"""

from __future__ import annotations

import gzip
import json
import subprocess
import sys
from pathlib import Path

import pytest

from biobox import ids
from biobox.l0seq.inputs import GenomeInput, write_genome_input

pytestmark = pytest.mark.covers("D-01")  # ring R1 by its directory (tests/conftest.py)

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "stress" / "build_stage_pod" / "make_inputs.py"
X, Y, Z = b"ACGTACGTTTGACCA", b"GGGTTTAAACCCGTA", b"TTTTGGGGCCCCAAAT"
U = {k: ids.contig_uid(v) for k, v in {"X": X, "Y": Y, "Z": Z}.items()}


def _fasta(path: Path, recs: list[tuple[str, bytes]]) -> None:
    path.write_bytes(b"".join(b">" + h.encode() + b"\n" + s + b"\n" for h, s in recs))


def _prok(tmp_path: Path) -> tuple[Path, Path]:
    """ACC1 = (X, Y); ACC2 = (Y, Z, Z): Y is shared with ACC1, Z repeats inside ACC2; ACC3 = (X, Y), ACC1's contigs."""
    files = {"ACC1": [("a1 chr", X), ("a1 pl", Y)], "ACC2": [("a2 chr", Y), ("a2 pl", Z), ("a2 pl copy", Z)],
             "ACC3": [("a3 chr", X), ("a3 pl", Y)]}
    lock = {"genomes": {}}
    for acc, recs in files.items():
        p = tmp_path / f"{acc}.fna"
        _fasta(p, recs)
        lock["genomes"][acc] = {"file": str(p), "retrieved_utc": "2026-11-09T00:00:00+00:00"}
    (tmp_path / "lock.json").write_text(json.dumps(lock))
    (tmp_path / "sel.tsv").write_text("accession\ttranslation_table\nACC1\t11\nACC2\t\nACC3\tNone\n")
    return tmp_path / "lock.json", tmp_path / "sel.tsv"


def _run(*args) -> dict:
    done = subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], capture_output=True, text=True, check=True)
    return json.loads(done.stdout)


def _genome(*names: str) -> dict:
    return {"contig_uids": [U[n].hex() for n in names], "genetic_code": 11, "kind": "prokaryote", "circular": [],
            "assembly_quality": "complete"}


def test_prok_lists_a_shared_contig_in_every_genome_and_reports_it(tmp_path):
    lock, sel = _prok(tmp_path)
    rep = _run("prok", lock, sel, tmp_path / "new")
    got = json.loads((tmp_path / "new" / "genomes.json").read_text())
    assert got == [_genome("X", "Y"), _genome("Y", "Z")]
    assert rep["genomes"] == 2 and rep["contigs"] == 3 and rep["memberships"] == 4
    assert rep["contigs_in_several_genomes"] == 1
    assert rep["records_repeated_in_genome"] == 1 and rep["genomes_repeated"] == 1
    # the stage's inputs step takes them (DEC-27): Y is in both genomes
    gs = [GenomeInput(tuple(bytes.fromhex(u) for u in g["contig_uids"]), g["genetic_code"], kind=g["kind"])
          for g in got]
    assert write_genome_input(tmp_path / "root", "S9", gs) == {"genomes": 2, "contigs": 3, "rows": 4,
                                                                "shared_contigs": 1}


def test_drop_repeats_gives_the_old_bytes(tmp_path):
    """The old rule: a contig uid already listed by an earlier genome (accessions sorted) is not listed again, a genome
    left with none is not listed, and the report counts the uids not listed."""
    lock, sel = _prok(tmp_path)
    rep = _run("prok", lock, sel, tmp_path / "old", "--drop-repeats")
    want = [_genome("X", "Y"), _genome("Z")]
    assert (tmp_path / "old" / "genomes.json").read_text() == json.dumps(want)
    assert rep == {"genomes": 2, "contigs": 3, "contigs_not_listed_again": 4, "objects": 3, "codes": {"11": 2}}
    _run("prok", lock, sel, tmp_path / "new")
    assert (tmp_path / "old" / "objects.json").read_bytes() == (tmp_path / "new" / "objects.json").read_bytes()


def test_phage_lists_a_repeated_record_once_in_either_mode(tmp_path):
    """Every phage record is its own genome, so a record that repeats another's sequence is the same genome: listed
    once by both rules, and the files are the same."""
    src = tmp_path / "src.fa.gz"
    with gzip.open(src, "wb") as fh:
        fh.write(b">r1\n" + X + b"\n>r2\n" + Y + b"\n>r3 same as r1\n" + X + b"\n>r4\n" + Z + b"\n")
    new = _run("phage", src, 3, tmp_path / "new")
    old = _run("phage", src, 3, tmp_path / "old", "--drop-repeats")
    vir = [{**_genome(n), "kind": "virus"} for n in ("X", "Y")]
    for d in ("new", "old"):
        assert json.loads((tmp_path / d / "genomes.json").read_text()) == vir
    assert (tmp_path / "new" / "genomes.json").read_bytes() == (tmp_path / "old" / "genomes.json").read_bytes()
    assert new["genomes_repeated"] == 1 and old["contigs_not_listed_again"] == 1


def test_help_names_the_old_rule_and_the_stages_made_with_it():
    out = subprocess.run([sys.executable, str(SCRIPT), "--help"], capture_output=True, text=True, check=True).stdout
    assert "--drop-repeats" in out
    for s in ("S1", "P1", "P6", "H1"):
        assert s in out


# --- card C-26: the PLSDB plasmid stage PL1 ---------------------------------------------------------------------------
# `make_inputs.py plasmid FASTA META_CSV OUT_DIR [--root ROOT]`: one genome per record (kind plasmid, genetic code 11),
# topology from the metadata's NUCCORE_Topology (circular, linear; empty is linear and counted), assembly quality from
# NUCCORE_Completeness (complete -> complete, empty -> contig). No record is dropped for repeating a sequence (DEC-27):
# a record whose bytes repeat an earlier record's is the same genome (genome_uid hashes the contig set), listed once
# and written to repeated_records.tsv; a record whose bytes an earlier stage of ROOT stores is listed as a genome of
# its own and written to shared_records.tsv with the stage that owns the contig (DEC-26). The oracles are the tests'
# own loops over the records they construct.

W = b"CCCCAAAATTTTGGGGACGTAAC"  # a contig an earlier stage stores
V = b"GATTACAGATTACAGATTACAT"
META_HEAD = ("NUCCORE_UID,NUCCORE_ACC,NUCCORE_Description,NUCCORE_Completeness,NUCCORE_Genome,NUCCORE_Length,"
             "NUCCORE_Source,NUCCORE_Topology,STATUS\n")
# (accession, sequence, completeness, topology, status): r4 repeats r1's bytes with another topology, r6 repeats r2's
# with the same metadata, r5 is the earlier stage's contig W, r3 has an empty topology and completeness
PLASMIDS = [("NZ_AA000001.1", X, "complete", "circular", ""), ("NZ_AA000002.1", Y, "complete", "linear", ""),
            ("CP000003.1", Z, "", "", ""), ("CP000004.1", X, "complete", "linear", "suppressed"),
            ("NZ_AA000005.1", W, "complete", "circular", ""), ("NZ_AA000006.1", Y, "complete", "linear", ""),
            ("NZ_AA000007.1", V, "", "circular", "")]


def _plasmid_files(tmp_path: Path, recs=PLASMIDS) -> tuple[Path, Path]:
    fa, meta = tmp_path / "plsdb.fasta", tmp_path / "nuccore.csv"
    _fasta(fa, [(f"{a} Some bacterium plasmid p{i}, complete sequence", s) for i, (a, s, *_) in enumerate(recs)])
    rows = [f'{i},{a},"Some bacterium plasmid p{i}, complete sequence",{c},plasmid,{len(s)},refseq,{t},{st}\n'
            for i, (a, s, c, t, st) in enumerate(recs)]
    meta.write_text(META_HEAD + "".join(reversed(rows)))  # the metadata's order is not the FASTA's
    return fa, meta


def _earlier_stage(tmp_path: Path) -> Path:
    """A root whose stage S9 stores a chromosome and W, ingested by the store's own code."""
    from biobox.l0seq import SourceObject, ingest

    src = tmp_path / "s9.fna"
    _fasta(src, [("chr1 a chromosome", b"ACGT" * 40), ("pl1 its plasmid", W)])
    root = tmp_path / "root"
    ingest([SourceObject(src, "test", "1", "test", True)], root, "S9")
    return root


def _tsv(path: Path) -> list[dict]:
    head, *rows = path.read_text().splitlines()
    return [dict(zip(head.split("\t"), r.split("\t"))) for r in rows]


def _first(recs) -> dict[bytes, int]:
    first: dict[bytes, int] = {}
    for i, (_, s, *_) in enumerate(recs):
        first.setdefault(s, i)
    return first


def _oracle(recs=PLASMIDS, stored=()) -> dict:
    """The expected counts, by a loop over the constructed records: the first record of each sequence is its genome."""
    first = _first(recs)
    topo = {"circular": 0, "linear": 0, "empty": 0}
    qual = {"complete": 0, "contig": 0}
    for _, _, c, t, _ in recs:
        topo[t or "empty"] += 1
        qual["complete" if c == "complete" else "contig"] += 1
    reps = [i for i, (_, s, *_) in enumerate(recs) if first[s] != i]
    return {"records": len(recs), "genomes": len(first), "topology": topo, "quality": qual,
            "genomes_circular": sum(1 for i in first.values() if recs[i][3] == "circular"),
            "repeated_records": len(reps),
            "repeats_other_topology": sum(1 for i in reps if recs[i][3] != recs[first[recs[i][1]]][3]),
            "suppressed": sum(1 for r in recs if r[4] == "suppressed"),
            "record_bases": sum(len(s) for _, s, *_ in recs), "genome_bases": sum(len(s) for s in first),
            "shared_records": sum(1 for _, s, *_ in recs if s in stored),
            "new_contigs": sum(1 for s in first if s not in stored),
            "new_bases": sum(len(s) for s in first if s not in stored)}


@pytest.mark.covers("C-26")
def test_plasmid_one_genome_per_record_with_topology_and_quality_from_the_metadata(tmp_path):
    fa, meta = _plasmid_files(tmp_path)
    rep = _run("plasmid", fa, meta, tmp_path / "out")
    want = _oracle()
    for k in ("records", "genomes", "topology", "quality", "genomes_circular", "repeated_records",
              "repeats_other_topology", "suppressed", "record_bases", "genome_bases"):
        assert rep[k] == want[k], k
    got = json.loads((tmp_path / "out" / "genomes.json").read_text())
    exp = []
    for s, i in sorted(_first(PLASMIDS).items(), key=lambda t: t[1]):  # file order, the first record's attributes
        u = ids.contig_uid(s).hex()
        exp.append({"contig_uids": [u], "genetic_code": 11, "kind": "plasmid",
                    "circular": [u] if PLASMIDS[i][3] == "circular" else [],
                    "assembly_quality": "complete" if PLASMIDS[i][2] == "complete" else "contig"})
    assert got == exp
    objs = json.loads((tmp_path / "out" / "objects.json").read_text())
    assert [(o["path"], o["source_name"], o["source_release"]) for o in objs] == [(str(fa), "plsdb", "2024_05_31_v2")]
    reps = _tsv(tmp_path / "out" / "repeated_records.tsv")
    assert [(r["accession"], r["first_accession"], r["same_topology"]) for r in reps] == [
        ("CP000004.1", "NZ_AA000001.1", "false"), ("NZ_AA000006.1", "NZ_AA000002.1", "true")]
    # the stage's inputs step takes them: one row per genome, the circular ones circular
    gs = [GenomeInput(tuple(bytes.fromhex(u) for u in g["contig_uids"]), g["genetic_code"], kind=g["kind"],
                      circular=tuple(bytes.fromhex(u) for u in g["circular"]),
                      assembly_quality=g["assembly_quality"]) for g in got]
    assert write_genome_input(tmp_path / "r", "PL9", gs) == {"genomes": want["genomes"], "contigs": want["genomes"]}


@pytest.mark.covers("C-26")
def test_plasmid_record_stored_by_an_earlier_stage_is_a_genome_of_its_own_and_reported(tmp_path):
    root = _earlier_stage(tmp_path)
    fa, meta = _plasmid_files(tmp_path)
    rep = _run("plasmid", fa, meta, tmp_path / "out", "--root", root)
    want = _oracle(stored={W})
    assert rep["genomes"] == want["genomes"]  # W's record is listed, not dropped (DEC-27)
    assert rep["shared_records"] == want["shared_records"] and rep["owner_records"] == {"S9": 1}
    assert rep["owner_contigs"] == {"S9": 1}
    assert rep["new_contigs"] == want["new_contigs"] and rep["new_bases"] == want["new_bases"]
    got = json.loads((tmp_path / "out" / "genomes.json").read_text())
    assert [ids.contig_uid(W).hex()] in [g["contig_uids"] for g in got]
    shared = _tsv(tmp_path / "out" / "shared_records.tsv")
    assert [(r["accession"], r["contig_uid"], r["owner_stage"]) for r in shared] == [
        ("NZ_AA000005.1", ids.contig_uid(W).hex(), "S9")]
    # without --root nothing is said about earlier stages, and the genome list is the same
    _run("plasmid", fa, meta, tmp_path / "plain")
    assert (tmp_path / "plain" / "genomes.json").read_bytes() == (tmp_path / "out" / "genomes.json").read_bytes()
    assert not (tmp_path / "plain" / "shared_records.tsv").exists()


@pytest.mark.covers("C-26")
@pytest.mark.parametrize("bad", ["topology", "completeness", "missing"])
def test_plasmid_refuses_metadata_it_cannot_map(tmp_path, bad):
    recs = list(PLASMIDS)
    if bad == "topology":
        recs[1] = recs[1][:3] + ("unknown",) + recs[1][4:]
    elif bad == "completeness":
        recs[1] = recs[1][:2] + ("partial",) + recs[1][3:]
    fa, meta = _plasmid_files(tmp_path, recs)
    if bad == "missing":
        meta.write_text("".join(ln for ln in meta.read_text().splitlines(keepends=True) if "NZ_AA000002.1" not in ln))
    done = subprocess.run([sys.executable, str(SCRIPT), "plasmid", str(fa), str(meta), str(tmp_path / "out")],
                          capture_output=True, text=True)
    assert done.returncode == 2
    assert "NZ_AA000002.1" in done.stderr
    assert not (tmp_path / "out" / "genomes.json").exists()


@pytest.mark.covers("C-26")
def test_help_names_the_plasmid_mode():
    out = subprocess.run([sys.executable, str(SCRIPT), "--help"], capture_output=True, text=True, check=True).stdout
    for w in ("plasmid", "NUCCORE_Topology", "NUCCORE_Completeness", "PL1"):
        assert w in out
