"""Card C-25 (1): `scripts/stress/build_stage_pod/c25_select.py`, the S2C selection by quality order.

The rule: complete genomes (`ncbi_assembly_level = 'Complete Genome'`), at most `cap` per `gtdb_taxonomy`, ranked
by `checkm2_completeness - 5 * checkm2_contamination` descending, a value that does not parse last, every tie by
GTDB accession. The oracle is a brute-force count written here: a genome's rank is one plus the number of complete
genomes of its species that beat it, by explicit pairwise comparison (no sort). The bounds on selected bp under any
other choice of the non-best genomes are checked against an enumeration of every allowed subset.
"""

from __future__ import annotations

import csv
import gzip
import itertools
import json
import random
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("D-01", "D-04")  # ring R1 by its directory (tests/conftest.py)

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "stress" / "build_stage_pod" / "c25_select.py"
HEADER = ["accession", "checkm2_completeness", "checkm2_contamination", "contig_count", "genome_size",
          "gtdb_genome_representative", "gtdb_representative", "gtdb_taxonomy", "ncbi_assembly_level",
          "ncbi_translation_table", "protein_count"]
PHASE1_COLS = ["accession", "genome_size", "contig_count", "protein_count", "translation_table", "gtdb_taxonomy",
               "gtdb_accession"]
CAP = 3


def _g(acc, sp, comp, cont, size, level="Complete Genome"):
    return {"accession": acc, "checkm2_completeness": comp, "checkm2_contamination": cont, "contig_count": "2",
            "genome_size": str(size), "gtdb_genome_representative": "RS_GCF_000000001.1", "gtdb_representative": "f",
            "gtdb_taxonomy": f"d__Bacteria;p__{sp[:1]}phylum;c__c;o__o;f__f;g__g;s__{sp}", "ncbi_assembly_level": level,
            "ncbi_translation_table": "11", "protein_count": str(size // 1000)}


def _hand_rows() -> list[dict]:
    # species A: 6 complete genomes; 100 - 5*0.2 == 99.5 - 5*0.1 == 99.0 - 5*0 == 99.0 exactly, so the second,
    # third and fourth tie and the GTDB accession decides (GB_ sorts before RS_): ranks 1, 4, 2, 3, 5; the sixth
    # has no CheckM2 value and sorts last (rank 6); the seventh is not complete.
    a = [_g("RS_GCF_000000010.1", "A", "100.0", "0.0", 5_000_000),
         _g("RS_GCF_000000020.1", "A", "100.0", "0.2", 5_100_000),
         _g("GB_GCA_000000030.1", "A", "99.5", "0.1", 4_900_000),
         _g("RS_GCF_000000005.1", "A", "99.0", "0.0", 5_300_000),
         _g("RS_GCF_000000040.1", "A", "90.0", "0.0", 4_000_000),
         _g("RS_GCF_000000001.1", "A", "", "", 6_000_000),
         _g("RS_GCF_000000002.1", "A", "100.0", "0.0", 5_000_000, level="Contig")]
    # species B: exactly CAP genomes; species C: one genome; species Dx: CAP + 1, two missing values.
    b = [_g(f"RS_GCF_00000010{i}.1", "B", "98.0", "1.0", 3_000_000 + i) for i in range(CAP)]
    c = [_g("GB_GCA_000000200.1", "C", "97.0", "0.5", 2_000_000)]
    d = [_g("RS_GCF_000000300.1", "Dx", "n/a", "0.0", 1_000_000), _g("RS_GCF_000000301.1", "Dx", "", "", 1_100_000),
         _g("RS_GCF_000000302.1", "Dx", "50.0", "0.0", 1_200_000), _g("RS_GCF_000000303.1", "Dx", "50", "0", 1_300_000)]
    return a + b + c + d


def _random_rows(seed: int = 20261112, n: int = 400) -> list[dict]:
    rng = random.Random(seed)
    rows = []
    for i in range(n):
        sp = f"R{rng.randrange(40)}"
        comp = rng.choice(["100.0", "99.5", "99.0", "", "95.25"])
        cont = rng.choice(["0.0", "0.1", "0.2", "1.5", ""])
        level = "Complete Genome" if rng.random() < 0.85 else "Scaffold"
        pre = rng.choice(["RS_GCF", "GB_GCA"])
        rows.append(_g(f"{pre}_{rng.randrange(10**9):09d}.{rng.randrange(1, 3)}", sp, comp, cont,
                       rng.randrange(500_000, 9_000_000), level))
    return rows


def _write_gz(path: Path, rows: list[dict]) -> None:
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.DictWriter(fh, HEADER, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def _score(r):
    try:
        return float(r["checkm2_completeness"]) - 5 * float(r["checkm2_contamination"])
    except ValueError:
        return None


def _beats(a, b) -> bool:
    """True if `a` is ranked before `b`: a parsed score beats a missing one; a higher score wins; else accession."""
    sa, sb = _score(a), _score(b)
    if sa is not None and sb is None:
        return True
    if sa is None and sb is not None:
        return False
    if sa is not None and sa != sb:
        return sa > sb
    return a["accession"] < b["accession"]


def _oracle(rows: list[dict], cap: int) -> dict[str, tuple[int, int]]:
    """accession (without prefix) -> (rank in species, complete genomes in species), for the selected genomes."""
    comp = [r for r in rows if r["ncbi_assembly_level"] == "Complete Genome"]
    out = {}
    for r in comp:
        same = [o for o in comp if o["gtdb_taxonomy"] == r["gtdb_taxonomy"]]
        rank = 1 + sum(_beats(o, r) for o in same if o is not r)
        if rank <= cap:
            out[r["accession"][3:]] = (rank, len(same))
    return out


def _run(tmp: Path, bac: list[dict], ar: list[dict], phase1: list[str], extra: list[str] = ()) -> Path:
    gt = tmp / "gtdb"
    gt.mkdir(parents=True, exist_ok=True)
    _write_gz(gt / "bac120_metadata_r226.tsv.gz", bac)
    _write_gz(gt / "ar53_metadata_r226.tsv.gz", ar)
    p1 = tmp / "phase1_all.tsv"
    p1.write_text("\t".join(PHASE1_COLS) + "\n" + "".join(line + "\n" for line in phase1))
    out = tmp / "out"
    subprocess.run([sys.executable, "-I", str(SCRIPT), str(gt), str(p1), str(out), "--cap", str(CAP), *extra],
                   check=True, capture_output=True, text=True)
    return out


def _phase1_lines(rows: list[dict]) -> list[str]:
    """The rank-1 genome of every species, in phase 1's columns, sorted by accession (the oracle's picks)."""
    best = {acc for acc, (rank, _) in _oracle(rows, 1).items()}
    lines = []
    for r in sorted((r for r in rows if r["accession"][3:] in best), key=lambda r: r["accession"][3:]):
        lines.append("\t".join([r["accession"][3:], r["genome_size"], r["contig_count"], r["protein_count"],
                                r["ncbi_translation_table"], r["gtdb_taxonomy"], r["accession"]]))
    return lines


def _read(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_cut_matches_brute_force_oracle(tmp_path):
    rows = _hand_rows() + _random_rows()
    bac = [r for r in rows if not r["gtdb_taxonomy"].endswith("s__C")]
    ar = [dict(r, gtdb_taxonomy=r["gtdb_taxonomy"].replace("d__Bacteria", "d__Archaea"))
          for r in rows if r["gtdb_taxonomy"].endswith("s__C")]
    out = _run(tmp_path, bac, ar, _phase1_lines(bac + ar))
    got = {r["accession"]: (int(r["species_rank"]), int(r["species_complete_genomes"]))
           for r in _read(out / "s2c_selected.tsv")}
    want = _oracle(bac + ar, CAP)
    assert got == want
    # the hand-made cases, spelled out
    assert {k: v[0] for k, v in got.items() if k[4:] < "000000400"} == {
        "GCF_000000010.1": 1, "GCA_000000030.1": 2, "GCF_000000005.1": 3,
        "GCF_000000100.1": 1, "GCF_000000101.1": 2, "GCF_000000102.1": 3, "GCA_000000200.1": 1,
        "GCF_000000302.1": 1, "GCF_000000303.1": 2, "GCF_000000300.1": 3}
    sel = _read(out / "s2c_selected.tsv")
    assert [r["accession"] for r in sel] == sorted(r["accession"] for r in sel)
    j = json.loads((out / "selection.json").read_text())
    assert j["phase1"]["rank1_equals_file"] is True
    assert j["phase1"]["not_in_selection"] == []
    rem = {r["accession"] for r in _read(out / "remainder.tsv")}
    assert rem == set(want) - {acc for acc, (rank, _) in want.items() if rank == 1}


def test_input_order_does_not_change_output(tmp_path):
    rows = _hand_rows() + _random_rows()
    shuffled = rows[:]
    random.Random(7).shuffle(shuffled)
    p1 = _phase1_lines(rows)
    a = _run(tmp_path / "a", rows, [], p1)
    b = _run(tmp_path / "b", shuffled, [], p1)
    for name in ("s2c_selected.tsv", "remainder.tsv", "by_phylum.tsv"):
        assert (a / name).read_bytes() == (b / name).read_bytes(), name
    ja, jb = (json.loads((d / "selection.json").read_text()) for d in (a, b))
    assert ja.pop("gtdb_files_sha256") != jb.pop("gtdb_files_sha256")  # the shuffled input file is other bytes
    assert ja == jb


def test_phase1_comparison_and_ani_bounds(tmp_path):
    rows = _hand_rows() + _random_rows(seed=3, n=120)
    p1 = _phase1_lines(rows)
    # a phase-1 genome outside the set (rank 5 of species A), one not in the metadata, one not complete
    p1 += ["GCF_000000040.1\t4000000\t2\t4000\t11\tx\tRS_GCF_000000040.1",
           "GCF_999999999.1\t1\t1\t1\t11\tx\tRS_GCF_999999999.1",
           "GCF_000000002.1\t5000000\t2\t5000\t11\tx\tRS_GCF_000000002.1"]
    p1.sort()
    out = _run(tmp_path, rows, [], p1, ["--phase1-not-served", "GCA_000000200.1,GCF_000000101.1"])
    j = json.loads((out / "selection.json").read_text())
    assert j["phase1"]["rank1_equals_file"] is False
    assert j["phase1"]["not_in_selection"] == [
        {"accession": "GCF_000000002.1", "reason": "not complete in the metadata (ncbi_assembly_level Contig)"},
        {"accession": "GCF_000000040.1", "reason": "rank 5 of 6 complete genomes of its species, beyond the cap"},
        {"accession": "GCF_999999999.1", "reason": "not in the metadata"}]
    ns = j["phase1"]["not_served"]
    assert ns == [{"accession": "GCA_000000200.1", "species_rank": 1, "species_complete_genomes": 1},
                  {"accession": "GCF_000000101.1", "species_rank": 2, "species_complete_genomes": CAP}]

    comp = [r for r in rows if r["ncbi_assembly_level"] == "Complete Genome"]
    by_sp: dict[str, list[dict]] = {}
    for r in comp:
        by_sp.setdefault(r["gtdb_taxonomy"], []).append(r)
    over = {sp: g for sp, g in by_sp.items() if len(g) > CAP}
    a = j["ani_step"]
    assert a["species_over_cap"] == len(over)
    assert a["genomes_in_species_over_cap"] == sum(len(g) for g in over.values())
    assert a["max_selected_differ_keep_rank1"] == sum(min(CAP - 1, len(g) - CAP) for g in over.values())
    assert a["max_selected_differ_any"] == sum(min(CAP, len(g) - CAP) for g in over.values())
    # bounds on selected bp over every allowed subset (the rank-1 genome kept), by enumeration
    lo = hi = 0
    for g in by_sp.values():
        best = next(r for r in g if not any(_beats(o, r) for o in g if o is not r))
        rest = [int(r["genome_size"]) for r in g if r is not best]
        sums = [sum(c) for c in itertools.combinations(rest, min(CAP - 1, len(rest)))]
        lo += int(best["genome_size"]) + min(sums)
        hi += int(best["genome_size"]) + max(sums)
    assert (a["bp_bounds_keep_rank1"]["min"], a["bp_bounds_keep_rank1"]["max"]) == (lo, hi)
    sel_bp = sum(int(r["genome_size"]) for r in _read(out / "s2c_selected.tsv"))
    assert a["bp_bounds_keep_rank1"]["quality_order"] == sel_bp == j["selected"]["bp"]
    assert sum(b["selected_genomes"] for b in j["by_species_size"]) == j["selected"]["genomes"]
    ph = _read(out / "by_phylum.tsv")
    assert sum(int(r["remainder_genomes"]) for r in ph) == j["remainder"]["genomes"]
