"""Card C-51: `scripts/stress/c51_deferred_list.py`, the per-accession list of what S2C phase 2's selection left out.

The rule (card C-27, DEC-34): in a GTDB species with more than `cap` complete genomes (`ncbi_assembly_level = 'Complete
Genome'`) the selection keeps `cap` (phase 1's genome and `cap - 1` picks, stages C1 to C6). The list is every other
complete genome of those species, one row per accession.

The oracle is a world small enough to work by hand (cap 2), written out below, with the answer written out too; the
script is never asked to reproduce its own earlier output.
"""

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("D-01")  # ring R1 by its directory (tests/conftest.py)

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "stress" / "c51_deferred_list.py"
HEADER = ["accession", "checkm2_completeness", "checkm2_contamination", "contig_count", "genome_size",
          "gtdb_representative", "gtdb_taxonomy", "ncbi_assembly_level"]
CAP = 2


def _g(acc, sp, level="Complete Genome"):
    return {"accession": acc, "checkm2_completeness": "99.0", "checkm2_contamination": "0.5", "contig_count": "1",
            "genome_size": "1000000", "gtdb_representative": "f",
            "gtdb_taxonomy": f"d__Bacteria;p__P;c__c;o__o;f__f;g__g;s__{sp}", "ncbi_assembly_level": level}


# Species A: 5 complete (a1..a5) and one scaffold (a6, never counted). phase 1 a1; picks a2 (C1). Left: a3 a4 a5.
# Species B: exactly CAP complete (b1 phase 1, b2 pick): nothing left, and not over the cap.
# Species C: 4 complete; phase 1 c1; pick c2; c3 was not served (a candidate): left, flagged; c4 left.
# Species D: 3 complete; phase 1 d0 is NOT among them (it is a contig-level assembly here); pick d1: left d2 d3, so
#   2 where complete - cap = 1: a difference to be named, with the reason.
# Species F: one complete genome, a pick that NCBI did not serve (so it is in no stage): not in the list.
# Species E: 3 complete; phase 1 e1; picks e2 (C2); left e3. In the other stage, so stage counts differ from A's.
ROWS = ([_g(f"RS_GCF_00000010{i}.1", "A") for i in range(1, 6)] + [_g("RS_GCF_000000169.1", "A", "Scaffold")]
        + [_g("RS_GCF_000000201.1", "B"), _g("RS_GCF_000000202.1", "B")]
        + [_g(f"RS_GCF_00000030{i}.1", "C") for i in range(1, 5)]
        + [_g("RS_GCF_000000400.1", "D", "Contig")] + [_g(f"RS_GCF_00000040{i}.1", "D") for i in range(1, 4)]
        + [_g(f"RS_GCF_00000050{i}.1", "E") for i in range(1, 4)]
        + [_g("RS_GCF_000000203.1", "F")])

PHASE1 = [("GCF_000000101.1", "A"), ("GCF_000000201.1", "B"), ("GCF_000000301.1", "C"), ("GCF_000000400.1", "D"),
          ("GCF_000000501.1", "E")]
PICKS = [("GCF_000000102.1", "C1", "A"), ("GCF_000000202.1", "C1", "B"), ("GCF_000000302.1", "C1", "C"),
         ("GCF_000000401.1", "C1", "D"), ("GCF_000000502.1", "C2", "E")]
# the genome-count table of card C-27 (complete genomes per species over the cap), as ani_species.tsv has it
ANI = [("A", 5), ("C", 4), ("D", 3), ("E", 3)]

EXPECTED_LEFT = {  # accession -> (species, not_served)
    "GCF_000000103.1": ("A", ""), "GCF_000000104.1": ("A", ""), "GCF_000000105.1": ("A", ""),
    "GCF_000000303.1": ("C", "y"), "GCF_000000304.1": ("C", ""),
    "GCF_000000402.1": ("D", ""), "GCF_000000403.1": ("D", ""),
    "GCF_000000503.1": ("E", ""),
}


def _tax(sp):
    return f"d__Bacteria;p__P;c__c;o__o;f__f;g__g;s__{sp}"


def _world(tmp: Path) -> dict[str, Path]:
    gtdb = tmp / "gtdb"
    gtdb.mkdir()
    with gzip.open(gtdb / "bac120_metadata_r226.tsv.gz", "wt", newline="") as fh:
        w = csv.DictWriter(fh, HEADER, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(ROWS)
    with gzip.open(gtdb / "ar53_metadata_r226.tsv.gz", "wt", newline="") as fh:
        fh.write("\t".join(HEADER) + "\n")
    p1 = tmp / "phase1.tsv"
    p1.write_text("accession\tgtdb_taxonomy\n" + "".join(f"{a}\t{_tax(s)}\n" for a, s in PHASE1))
    st = tmp / "stages.tsv"
    st.write_text("accession\tstage\tposition\tspecies\tmember\n"
                  + "".join(f"{a}\t{s}\t{i}\ts__{sp}\tx\n" for i, (a, s, sp) in enumerate(PICKS)))
    ani = tmp / "ani_species.tsv"
    ani.write_text("species\tcomplete_genomes\n" + "".join(f"{_tax(s)}\t{n}\n" for s, n in ANI))
    ns = tmp / "not_served_candidates.txt"
    ns.write_text("GCF_000000303.1\n")
    nsp = tmp / "not_served_picks.txt"
    nsp.write_text("GCF_000000203.1\n")
    return {"gtdb": gtdb, "p1": p1, "st": st, "ani": ani, "ns": ns, "nsp": nsp}


def _run(w, out: Path):
    return subprocess.run(
        [sys.executable, "-I", str(SCRIPT), str(w["gtdb"]), str(w["p1"]), str(w["st"]), str(out), "--cap", str(CAP),
         "--species-table", str(w["ani"]), "--not-served-candidates", str(w["ns"]),
         "--not-served-picks", str(w["nsp"])], capture_output=True, text=True)


def _read(p: Path):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_list_is_the_hand_worked_remainder(tmp_path):
    out = tmp_path / "out"
    r = _run(_world(tmp_path), out)
    assert r.returncode == 0, r.stderr
    rows = _read(out / "left_out.tsv")
    assert [x["accession"] for x in rows] == sorted(EXPECTED_LEFT)  # one row per accession, sorted
    for x in rows:
        sp, flag = EXPECTED_LEFT[x["accession"]]
        assert x["species"] == f"s__{sp}"
        assert x["assembly_level"] == "Complete Genome"
        assert x["not_served"] == flag
    assert list(rows[0]) == ["accession", "species", "assembly_level", "not_served"]


def test_counts_name_the_species_that_differ_and_why(tmp_path):
    out = tmp_path / "out"
    assert _run(_world(tmp_path), out).returncode == 0
    sp = {x["species"]: x for x in _read(out / "species_check.tsv")}
    assert set(sp) == {_tax(s) for s, _ in ANI}  # the species over the cap, and no other
    for s, n in ANI:
        x = sp[_tax(s)]
        assert int(x["complete"]) == n
        assert int(x["expected_left"]) == n - CAP
    assert [int(sp[_tax(s)]["left"]) for s in "ACDE"] == [3, 2, 2, 1]
    assert [sp[_tax(s)]["differs"] for s in "ACDE"] == ["", "", "y", ""]
    assert "phase-1 genome" in sp[_tax("D")]["why"]
    summary = json.loads((out / "summary.json").read_text())
    assert summary["left_out"] == 8 and summary["species_differ"] == 1
    assert summary["not_served_in_left_out"] == 1
    # the pick that failed to be served is not in the list; its species has no stage
    assert summary["not_served_picks_listed"] == 1
    assert summary["not_served_picks_by_stage"] == {"none": 1}
    assert summary["by_stage"] == {"C1": 7, "C2": 1}


def test_sorted_and_independent_of_input_order(tmp_path):
    a = tmp_path / "a"
    b = tmp_path / "b"
    a.mkdir()
    b.mkdir()
    wa = _world(a)
    wb = _world(b)
    # reverse the stage table's rows and the metadata's row order in the second world
    lines = wb["st"].read_text().splitlines()
    wb["st"].write_text("\n".join([lines[0]] + lines[:0:-1]) + "\n")
    with gzip.open(wb["gtdb"] / "bac120_metadata_r226.tsv.gz", "wt", newline="") as fh:
        w = csv.DictWriter(fh, HEADER, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(reversed(ROWS))
    assert _run(wa, a / "out").returncode == 0
    assert _run(wb, b / "out").returncode == 0
    for name in ("left_out.tsv", "species_check.tsv"):
        assert (a / "out" / name).read_bytes() == (b / "out" / name).read_bytes(), name
    # the summary differs only in the digests of the inputs, whose bytes were rewritten in another order
    sa, sb = (json.loads((x / "out" / "summary.json").read_text()) for x in (a, b))
    for s in (sa, sb):
        s.pop("gtdb_files_sha256")
        s.pop("inputs_sha256")
    assert sa == sb


def test_refuses_a_pick_that_is_in_the_list_of_the_left(tmp_path):
    w = _world(tmp_path)
    w["ns"].write_text("GCF_000000102.1\n")  # a pick named as a not-served candidate: the lists contradict
    r = _run(w, tmp_path / "out")
    assert r.returncode != 0
    assert "GCF_000000102.1" in r.stderr
