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
