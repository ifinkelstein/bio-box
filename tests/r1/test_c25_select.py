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


NEEDS = SCRIPT.with_name("c25_ani_needs.py")


def test_ani_needs_counts_match_brute_force(tmp_path):
    """Card C-25 (1b): what a skani step would read and compute, against counts written here.

    For each species with more than CAP complete genomes: n genomes; all pairs, counted by enumerating
    itertools.combinations; and the pairs of a greedy farthest-first pick from rank 1, counted by simulating it
    (after each of the first CAP-1 picks, the distances from that pick to every genome not yet picked).
    """
    rows = _hand_rows() + _random_rows(seed=11, n=300)
    # species A's representative is one of its complete genomes; the other species' representative is not theirs
    rows = [dict(r, gtdb_genome_representative="RS_GCF_000000010.1") if r["gtdb_taxonomy"].endswith("s__A") else r
            for r in rows]
    shuffled = rows[:]
    random.Random(5).shuffle(shuffled)
    outs = []
    for name, data in (("a", rows), ("b", shuffled)):
        gt = tmp_path / name
        gt.mkdir()
        _write_gz(gt / "bac120_metadata_r226.tsv.gz", data)
        _write_gz(gt / "ar53_metadata_r226.tsv.gz", [])
        p = subprocess.run([sys.executable, "-I", str(NEEDS), str(gt), "--cap", str(CAP)],
                           check=True, capture_output=True, text=True)
        j = json.loads(p.stdout)
        j.pop("gtdb_files_sha256")
        outs.append(j)
    assert outs[0] == outs[1]
    j = outs[0]

    comp = [r for r in rows if r["ncbi_assembly_level"] == "Complete Genome"]
    by_sp: dict[str, list[dict]] = {}
    for r in comp:
        by_sp.setdefault(r["gtdb_taxonomy"], []).append(r)
    over = {sp: g for sp, g in by_sp.items() if len(g) > CAP}
    assert len(over) >= 3
    assert j["species_over_cap"] == len(over)
    assert j["genomes"] == sum(len(g) for g in over.values())
    assert j["bp"] == sum(int(r["genome_size"]) for g in over.values() for r in g)
    assert j["pairs_all"] == sum(1 for g in over.values() for _ in itertools.combinations(g, 2))
    ff = 0
    for g in over.values():
        picked, left = 1, len(g) - 1
        while picked < CAP:
            ff += left  # distances from the newest pick to every genome not yet picked
            picked, left = picked + 1, left - 1
    assert j["pairs_farthest_first"] == ff
    best = {sp: next(r for r in g if not any(_beats(o, r) for o in g if o is not r)) for sp, g in over.items()}
    assert j["rank1_held_by_phase1"] == {"genomes": len(over),
                                         "bp": sum(int(r["genome_size"]) for r in best.values())}
    assert j["to_download"] == {"genomes": j["genomes"] - len(over),
                                "bp": j["bp"] - j["rank1_held_by_phase1"]["bp"]}
    ranked = {sp: sorted(g, key=lambda r: sum(_beats(o, r) for o in g if o is not r)) for sp, g in over.items()}
    assert j["quality_ranks_2_to_cap"] == {"genomes": sum(len(g[1:CAP]) for g in ranked.values()),
                                           "bp": sum(int(r["genome_size"]) for g in ranked.values() for r in g[1:CAP])}
    assert j["representative_complete_species"] == sum(
        any(r["accession"] == r["gtdb_genome_representative"] for r in g) for g in over.values())
    assert j["representative_complete_species"] >= 1


# ---------------------------------------------------------------------------------------------------------------
# Card C-27: `c25_select.py --ani` (greedy farthest-first on skani ANI) and `--candidates`, the skani output
# normaliser and driver (`c27_skani.py`), and the fetch (`c27_fetch.py`).
# ---------------------------------------------------------------------------------------------------------------

import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import os  # noqa: E402
import zipfile  # noqa: E402
from fractions import Fraction  # noqa: E402

SKANI = SCRIPT.with_name("c27_skani.py")
FETCH = SCRIPT.with_name("c27_fetch.py")
SKANI_HEADER = ["Ref_file", "Query_file", "ANI", "Align_fraction_ref", "Align_fraction_query", "Ref_name",
                "Query_name"]

# Digests of the quality-only outputs of c25_select.py as it was before card C-27 (the file at commit 00e313e),
# on the fixture of test_quality_only_output_unchanged_without_ani. A no-change check, not an oracle: the
# correctness of those outputs is checked against the brute-force oracle by the tests above.
UNCHANGED = {
    "s2c_selected.tsv": "a60222bf19b0ce8943f27202d847eeb4de1aebf1dbe75eb27d442bdb4f303b58",
    "remainder.tsv": "f632660a867cea5a0854d9cdf1ad6ec6c0dd6dde0bbd09712c42806d8f2f06fc",
    "by_phylum.tsv": "a4687e468818fa6cb89aa95c7a89063e4e46c57a4002a6f38408f5ae2012d92f",
    "selection.json (sorted keys, without gtdb_files_sha256)":
        "54e55c87ec1a23a8b62f6caa4c5e7e7bcead91418aa66f83303bad952d21a2bc",
}


def _run_cap(tmp: Path, rows: list[dict], phase1: list[str], cap: int, extra: list[str] = (),
             check: bool = True) -> tuple[Path, subprocess.CompletedProcess]:
    gt = tmp / "gtdb"
    gt.mkdir(parents=True, exist_ok=True)
    _write_gz(gt / "bac120_metadata_r226.tsv.gz", rows)
    _write_gz(gt / "ar53_metadata_r226.tsv.gz", [])
    p1 = tmp / "phase1_all.tsv"
    p1.write_text("\t".join(PHASE1_COLS) + "\n" + "".join(line + "\n" for line in phase1))
    out = tmp / "out"
    p = subprocess.run([sys.executable, "-I", str(SCRIPT), str(gt), str(p1), str(out), "--cap", str(cap), *extra],
                       check=False, capture_output=True, text=True)
    if check and p.returncode != 0:
        raise AssertionError(f"exit {p.returncode}\n{p.stderr}")
    return out, p


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def test_quality_only_output_unchanged_without_ani(tmp_path):
    rows = _hand_rows() + _random_rows()
    out, _ = _run_cap(tmp_path, rows, _phase1_lines(rows), CAP, ["--phase1-not-served", "GCA_000000200.1"])
    j = json.loads((out / "selection.json").read_text())
    j.pop("gtdb_files_sha256")
    got = {name: _sha((out / name).read_bytes()) for name in ("s2c_selected.tsv", "remainder.tsv", "by_phylum.tsv")}
    got["selection.json (sorted keys, without gtdb_files_sha256)"] = _sha(json.dumps(j, sort_keys=True).encode())
    assert got == UNCHANGED
    assert sorted(p.name for p in out.iterdir()) == ["by_phylum.tsv", "remainder.tsv", "s2c_selected.tsv",
                                                     "selection.json"]


# --- an ANI world built here ---------------------------------------------------------------------------------

ANI_CAP = 5


def _ani_rows(seed: int = 27) -> list[dict]:
    """Species W (14 genomes, tied scores), X (9), Y (6, one over ANI_CAP), Z (ANI_CAP, not over), V (2)."""
    rng = random.Random(seed)
    rows = []
    for sp, n in (("W", 14), ("X", 9), ("Y", 6), ("Z", ANI_CAP), ("V", 2)):
        for _ in range(n):
            pre = rng.choice(["RS_GCF", "GB_GCA"])
            comp = rng.choice(["100.0", "99.0"]) if sp == "W" else f"{90 + rng.random() * 10:.2f}"
            rows.append(_g(f"{pre}_{rng.randrange(10**9):09d}.1", sp, comp, "0.0", rng.randrange(1_000_000, 6_000_000)))
    return rows


def _ani_world(rows: list[dict], seed: int = 99) -> dict[tuple[str, str], tuple]:
    """(a, b) with a < b -> (ANI string or None for 'skani reported nothing', af_a, af_b), within each species.

    ANI has two decimals as skani prints it and is drawn from a few values so that ties occur; some pairs are
    left unreported; a few have both alignment fractions under 15."""
    rng = random.Random(seed)
    comp = [r for r in rows if r["ncbi_assembly_level"] == "Complete Genome"]
    by_sp: dict[str, list[str]] = {}
    for r in comp:
        by_sp.setdefault(r["gtdb_taxonomy"], []).append(r["accession"][3:])
    world = {}
    for accs in by_sp.values():
        for a, b in itertools.combinations(sorted(accs), 2):
            u = rng.random()
            if u < 0.08:
                world[(a, b)] = (None, "", "")
            elif u < 0.12:
                world[(a, b)] = (rng.choice(["99.10", "98.75"]), f"{rng.uniform(1, 14.9):.2f}",
                                 f"{rng.uniform(1, 14.9):.2f}")
            else:
                world[(a, b)] = (rng.choice(["99.98", "99.50", "99.10", "98.75", "97.20", "96.01"]),
                                 f"{rng.uniform(15, 100):.2f}", f"{rng.uniform(15, 100):.2f}")
    return world


def _skani_text(world: dict, rng: random.Random, paths: dict[str, str]) -> str:
    """The world as skani's `dist` output: a header and one row per reported pair, ref and query in either order."""
    body = []
    for (a, b), (ani, af_a, af_b) in world.items():
        if ani is None:
            continue
        if rng.random() < 0.5:
            body.append("\t".join([paths[a], paths[b], ani, af_a, af_b, f"NZ_CP{a[4:10]}.1 {a} chromosome",
                                   f"NZ_CP{b[4:10]}.1 {b} chromosome, complete genome"]))
        else:
            body.append("\t".join([paths[b], paths[a], ani, af_b, af_a, f"NZ_CP{b[4:10]}.1 x", f"NZ_CP{a[4:10]}.1 y"]))
    rng.shuffle(body)
    return "\n".join(["\t".join(SKANI_HEADER), *body]) + "\n"


def _normalised_text(world: dict, rng: random.Random) -> str:
    """The world as the normalised table c25_select.py reads, rows in a random order and a/b in either order."""
    out = []
    for (a, b), (ani, af_a, af_b) in world.items():
        row = [a, b, ani or "", af_a, af_b] if rng.random() < 0.5 else [b, a, ani or "", af_b, af_a]
        out.append("\t".join(row))
    rng.shuffle(out)
    return "genome_a\tgenome_b\tani\taf_a\taf_b\n" + "".join(x + "\n" for x in out)


def _dist(world: dict, a: str, b: str, min_af: float) -> Fraction:
    """100 - ANI, exact; the largest distance (100) when skani reported nothing or both fractions are under min_af."""
    ani, af_a, af_b = world.get((min(a, b), max(a, b)), (None, "", ""))
    if ani is None or max(Fraction(af_a), Fraction(af_b)) < Fraction(min_af):
        return Fraction(100)
    return 100 - Fraction(ani)


def _ff_oracle(rows: list[dict], world: dict, cap: int, min_af: float = 15) -> dict[str, list[str]]:
    """species -> the picks, written without the script's code: quality order by the pairwise _beats count; at each
    step every unpicked genome's distance to the picks is recomputed from scratch over all picks, and a genome
    replaces the current best only if strictly farther, or equally far and earlier in quality order."""
    comp = [r for r in rows if r["ncbi_assembly_level"] == "Complete Genome"]
    by_sp: dict[str, list[dict]] = {}
    for r in comp:
        by_sp.setdefault(r["gtdb_taxonomy"], []).append(r)
    out = {}
    for sp, g in by_sp.items():
        rank = {r["accession"][3:]: 1 + sum(_beats(o, r) for o in g if o is not r) for r in g}
        if len(g) <= cap:
            out[sp] = sorted(rank, key=rank.get)
            continue
        picks = [next(a for a, k in rank.items() if k == 1)]
        while len(picks) < cap:
            best = best_d = None
            for c in rank:
                if c in picks:
                    continue
                d = min(_dist(world, p, c, min_af) for p in picks)
                if best is None or d > best_d or (d == best_d and rank[c] < rank[best]):
                    best, best_d = c, d
            picks.append(best)
        out[sp] = picks
    return out


def _picked(out: Path) -> dict[str, list[str]]:
    by_sp: dict[str, list[tuple]] = {}
    for r in _read(out / "s2c_selected.tsv"):
        by_sp.setdefault(r["gtdb_taxonomy"], []).append(
            (int(r["ani_pick"]) if r["ani_pick"] else int(r["species_rank"]), r["accession"]))
    return {sp: [a for _, a in sorted(v)] for sp, v in by_sp.items()}


def _species_of(rows: list[dict], acc: str) -> str:
    return next(r["gtdb_taxonomy"] for r in rows if r["accession"][3:] == acc)


def test_ani_farthest_first_matches_brute_force_oracle(tmp_path):
    rows = _ani_rows()
    world = _ani_world(rows)
    table = tmp_path / "pairs.tsv"
    table.write_text(_normalised_text(world, random.Random(1)))
    p1 = _phase1_lines(rows)
    out, _ = _run_cap(tmp_path, rows, p1, ANI_CAP, ["--ani", str(table)])
    want = _ff_oracle(rows, world, ANI_CAP)
    assert _picked(out) == want
    quality = _oracle(rows, ANI_CAP)
    q_by_sp = {sp: sorted((a for a in quality if _species_of(rows, a) == sp), key=lambda a: quality[a][0])
               for sp in want}
    n_sp = {sp: sum(r["gtdb_taxonomy"] == sp for r in rows) for sp in want}
    over = sorted(sp for sp, n in n_sp.items() if n > ANI_CAP)
    assert len(over) == 3
    # the rule differs from quality order somewhere in this world, or the test proves little
    assert any(set(want[sp]) != set(q_by_sp[sp]) for sp in over)
    j = json.loads((out / "selection.json").read_text())
    a = j["ani"]
    assert a["phase1_subset"] is True and j["phase1"]["not_in_selection"] == []
    assert a["min_af"] == 15 and a["species_over_cap"] == 3
    assert a["consulted_pairs"] == sum((ANI_CAP - 1) * n_sp[sp] - ANI_CAP * (ANI_CAP - 1) // 2 for sp in over)
    assert a["consulted_absent"] == 0
    asel = {r["accession"] for r in _read(out / "s2c_selected.tsv")}
    assert a["selected_differ_from_quality"] == len(asel - set(quality)) == len(set(quality) - asel)
    for sp in want:  # species at or under the cap keep every genome
        if sp not in over:
            assert set(want[sp]) == {acc for acc in asel if _species_of(rows, acc) == sp}
    p1_acc = {line.split("\t")[0] for line in p1}
    assert {r["accession"] for r in _read(out / "remainder.tsv")} == asel - p1_acc
    # per-species minimum pairwise ANI among the picks under each rule, from the world
    sp_rows = _read(out / "ani_species.tsv")
    assert sorted(r["species"] for r in sp_rows) == over
    for r in sp_rows:
        sp = r["species"]
        assert int(r["differ"]) == len(set(want[sp]) - set(q_by_sp[sp]))
        for rule, picks in (("ani", want[sp]), ("quality", q_by_sp[sp])):
            vals = [world[(min(x, y), max(x, y))] for x, y in itertools.combinations(picks, 2)]
            rep = [Fraction(v[0]) for v in vals if v[0] is not None and max(Fraction(v[1]), Fraction(v[2])) >= 15]
            assert int(r[f"{rule}_pairs_reported"]) == len(rep)
            assert int(r[f"{rule}_pairs_unreported"]) == len(vals) - len(rep)
            assert int(r[f"{rule}_pairs_absent"]) == 0
            got = r[f"{rule}_min_pairwise_ani"]
            assert (Fraction(got) if got else None) == (min(rep) if rep else None)


def test_ani_min_af_changes_the_rule_as_stated(tmp_path):
    rows = _ani_rows(seed=4)
    world = _ani_world(rows, seed=5)
    table = tmp_path / "pairs.tsv"
    table.write_text(_normalised_text(world, random.Random(2)))
    out, _ = _run_cap(tmp_path, rows, _phase1_lines(rows), ANI_CAP, ["--ani", str(table), "--min-af", "0"])
    assert _picked(out) == _ff_oracle(rows, world, ANI_CAP, min_af=0)
    assert json.loads((out / "selection.json").read_text())["ani"]["min_af"] == 0


def test_ani_input_order_does_not_change_output(tmp_path):
    rows = _ani_rows(seed=8)
    world = _ani_world(rows, seed=9)
    shuffled = rows[:]
    random.Random(3).shuffle(shuffled)
    p1 = _phase1_lines(rows)
    outs = []
    for name, data, trng in (("a", rows, 1), ("b", shuffled, 2)):
        table = tmp_path / f"pairs_{name}.tsv"
        table.write_text(_normalised_text(world, random.Random(trng)))
        out, _ = _run_cap(tmp_path / name, data, p1, ANI_CAP, ["--ani", str(table)])
        outs.append(out)
    for name in ("s2c_selected.tsv", "remainder.tsv", "by_phylum.tsv", "ani_species.tsv"):
        assert (outs[0] / name).read_bytes() == (outs[1] / name).read_bytes(), name
    ja, jb = (json.loads((d / "selection.json").read_text()) for d in outs)
    for j in (ja, jb):
        j.pop("gtdb_files_sha256")
        j["ani"].pop("table_sha256")
        j["ani"].pop("table")  # pairs_a.tsv, pairs_b.tsv
    assert ja == jb


def test_ani_absent_pairs_count_as_largest_or_refuse(tmp_path):
    rows = _ani_rows(seed=12)
    world = _ani_world(rows, seed=13)
    rng = random.Random(4)
    dropped = {k for k in world if rng.random() < 0.3}  # not even an unreported row: same picks as 'not reported'
    table = tmp_path / "pairs.tsv"
    table.write_text(_normalised_text({k: v for k, v in world.items() if k not in dropped}, random.Random(5)))
    out, _ = _run_cap(tmp_path / "a", rows, _phase1_lines(rows), ANI_CAP, ["--ani", str(table)])
    as_missing = {k: ((None, "", "") if k in dropped else v) for k, v in world.items()}
    assert _picked(out) == _ff_oracle(rows, as_missing, ANI_CAP)
    assert json.loads((out / "selection.json").read_text())["ani"]["consulted_absent"] > 0
    _, p = _run_cap(tmp_path / "b", rows, _phase1_lines(rows), ANI_CAP, ["--ani", str(table), "--require-pairs"],
                    check=False)
    assert p.returncode != 0 and "not in the ANI table" in p.stderr


def test_ani_table_with_a_stranger_is_refused(tmp_path):
    rows = _ani_rows(seed=14)
    world = _ani_world(rows, seed=15)
    by_sp: dict[str, list[str]] = {}
    for r in rows:
        by_sp.setdefault(r["gtdb_taxonomy"], []).append(r["accession"][3:])
    sw = next(v[0] for k, v in by_sp.items() if k.endswith("s__W"))
    sx = next(v[0] for k, v in by_sp.items() if k.endswith("s__X"))
    for extra, msg in (({(min(sw, sx), max(sw, sx)): ("99.00", "90", "90")}, "two species"),
                       ({("GCF_999999998.1", "GCF_999999999.1"): ("99.00", "90", "90")}, "not a complete genome")):
        table = tmp_path / "pairs.tsv"
        table.write_text(_normalised_text(world | extra, random.Random(6)))
        _, p = _run_cap(tmp_path / msg.replace(" ", "_"), rows, _phase1_lines(rows), ANI_CAP, ["--ani", str(table)],
                        check=False)
        assert p.returncode != 0 and msg in p.stderr, p.stderr


def test_candidates_and_fetch_list(tmp_path):
    rows = _hand_rows() + _random_rows(seed=21, n=200)
    p1 = _phase1_lines(rows)
    first = sorted(line.split("\t")[0] for line in p1)
    not_served = [first[0], first[1]]
    out, _ = _run_cap(tmp_path, rows, p1, CAP, ["--candidates", "--phase1-not-served", ",".join(not_served)])
    comp = [r for r in rows if r["ncbi_assembly_level"] == "Complete Genome"]
    size: dict[str, int] = {}
    for r in comp:
        size[r["gtdb_taxonomy"]] = size.get(r["gtdb_taxonomy"], 0) + 1
    cand = sorted(r["accession"][3:] for r in comp if size[r["gtdb_taxonomy"]] > CAP)
    got = _read(out / "ani_candidates.tsv")
    assert [r["accession"] for r in got] == cand
    assert all(int(r["species_complete_genomes"]) > CAP for r in got)
    p1_acc = {line.split("\t")[0] for line in p1}
    fetch = _read(out / "ani_fetch.tsv")
    want = sorted(set(cand) - p1_acc | set(not_served))
    assert [r["accession"] for r in fetch] == want
    assert {r["accession"] for r in fetch if r["fetch_reason"] == "phase1_not_served_retry"} == set(not_served)
    assert {r["fetch_reason"] for r in fetch if r["accession"] not in not_served} == {"ani_candidate"}
    gsize = {r["accession"][3:]: r["genome_size"] for r in comp}
    assert all(r["genome_size"] == gsize[r["accession"]] for r in fetch)
    j = json.loads((out / "selection.json").read_text())["candidates"]
    assert j["genomes"] == len(cand) and j["fetch_genomes"] == len(want)
    assert j["fetch_bp"] == sum(int(gsize[a]) for a in want)
    # the quality-only files do not change when --candidates is given
    out2, _ = _run_cap(tmp_path / "plain", rows, p1, CAP, ["--phase1-not-served", ",".join(not_served)])
    for name in ("s2c_selected.tsv", "remainder.tsv", "by_phylum.tsv"):
        assert (out / name).read_bytes() == (out2 / name).read_bytes()


# --- the skani output normaliser -------------------------------------------------------------------------------

def test_skani_normalise_hand_made_file(tmp_path):
    """A file in skani's documented `dist` format (README: Ref_file, Query_file, ANI, Align_fraction_ref,
    Align_fraction_query, Ref_name, Query_name; names hold spaces), with a reversed duplicate, a sketch path and
    an extra trailing column, normalised to genome_a < genome_b; with --expect-*, pairs asked for and not
    reported become rows with an empty ANI."""
    raw = tmp_path / "skani_dist.tsv"
    raw.write_text(
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\tExtra\n"
        "/g/B01/GCF_000005845.2.fna.gz\t/g/B02/GCA_000008865.2.fna.gz\t99.39\t93.95\t92.21\t"
        "NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome\tBA000007.3 Escherichia coli O157\tx\n"
        "/g/B02/GCA_000008865.2.fna.gz\t/g/B01/GCF_000005845.2.fna.gz\t99.41\t92.20\t93.90\tBA000007.3 E\tNC E\tx\n"
        "sk/GCF_000006925.2.fna.gz.sketch\t/g/B01/GCF_000005845.2.fna.gz\t97.05\t88.10\t80.00\tNC_004337.2 S\tNC E\tx\n")
    out = tmp_path / "norm.tsv"
    subprocess.run([sys.executable, "-I", str(SKANI), "normalise", str(raw), "-o", str(out)], check=True,
                   capture_output=True, text=True)
    assert out.read_text() == (
        "genome_a\tgenome_b\tani\taf_a\taf_b\n"
        "GCA_000008865.2\tGCF_000005845.2\t99.41\t92.20\t93.90\n"
        "GCF_000005845.2\tGCF_000006925.2\t97.05\t80.00\t88.10\n")
    refs = tmp_path / "refs.txt"
    refs.write_text("GCA_000008865.2\nGCF_000006925.2\nGCF_000010001.1\n")
    out2 = tmp_path / "norm2.tsv"
    subprocess.run([sys.executable, "-I", str(SKANI), "normalise", str(raw), "-o", str(out2),
                    "--expect-query", "GCF_000005845.2", "--expect-refs", str(refs)], check=True,
                   capture_output=True, text=True)
    assert out2.read_text().splitlines()[1:] == [
        "GCA_000008865.2\tGCF_000005845.2\t99.41\t92.20\t93.90",
        "GCF_000005845.2\tGCF_000006925.2\t97.05\t80.00\t88.10",
        "GCF_000005845.2\tGCF_000010001.1\t\t\t"]
    bad = tmp_path / "bad.tsv"
    bad.write_text("Ref_file\tQuery_file\tANI\n/x/GCF_000005845.2.fna\t/x/GCF_000006925.2.fna\t99\n")
    p = subprocess.run([sys.executable, "-I", str(SKANI), "normalise", str(bad), "-o", str(tmp_path / "n3.tsv")],
                       capture_output=True, text=True)
    assert p.returncode != 0 and "Align_fraction_ref" in p.stderr


def test_skani_normalised_world_selects_like_oracle(tmp_path):
    """skani-format text of a world -> normaliser -> c25_select.py --ani gives the oracle's picks."""
    rows = _ani_rows(seed=31)
    world = _ani_world(rows, seed=32)
    paths = {r["accession"][3:]: f"/tmp/g/B{i % 3:02d}/{r['accession'][3:]}.fna.gz" for i, r in enumerate(rows)}
    raw = tmp_path / "raw.tsv"
    raw.write_text(_skani_text(world, random.Random(7), paths))
    norm = tmp_path / "norm.tsv"
    subprocess.run([sys.executable, "-I", str(SKANI), "normalise", str(raw), "-o", str(norm)], check=True,
                   capture_output=True, text=True)
    out, _ = _run_cap(tmp_path, rows, _phase1_lines(rows), ANI_CAP, ["--ani", str(norm)])
    assert _picked(out) == _ff_oracle(rows, world, ANI_CAP)


# --- the skani driver: dry-run plan, and a run against a stand-in for the binary --------------------------------

def _driver_inputs(tmp: Path, rows: list[dict], cap: int) -> tuple[Path, Path]:
    """Candidates from c25_select.py --candidates; genome files (empty) listed in a lock.json; paths.tsv from it."""
    out, _ = _run_cap(tmp / "sel", rows, _phase1_lines(rows), cap, ["--candidates"])
    gdir = tmp / "genomes"
    gdir.mkdir()
    lock = {"genomes": {}}
    for r in _read(out / "ani_candidates.tsv"):
        f = gdir / f"{r['accession']}.fna.gz"
        f.write_bytes(b"")
        lock["genomes"][r["accession"]] = {"file": str(f), "file_name": f.name}
    (gdir / "lock.json").write_text(json.dumps(lock))
    paths = tmp / "paths.tsv"
    subprocess.run([sys.executable, "-I", str(SKANI), "paths", str(paths), str(gdir / "lock.json")], check=True,
                   capture_output=True, text=True)
    assert len(_read(paths)) == len(lock["genomes"])
    return out / "ani_candidates.tsv", paths


def test_skani_driver_dry_run_plan(tmp_path):
    rows = _ani_rows(seed=41) + _random_rows(seed=42, n=150)
    cand, paths = _driver_inputs(tmp_path, rows, ANI_CAP)
    work = tmp_path / "work"
    args = [str(cand), str(paths), str(work), "--cap", str(ANI_CAP), "--budget-cpu-min", "1.0", "--sketch-cpu-s", "2",
            "--pair-ms", "10", "--threads", "4", "--image", "/img/bio-box.sif"]
    p = subprocess.run([sys.executable, "-I", str(SKANI), "plan", *args, "--dry-run"], check=True, capture_output=True,
                       text=True)
    assert not work.exists()  # a dry run writes nothing
    plan = json.loads(p.stdout.split("\n#", 1)[0])
    comp = [r for r in rows if r["ncbi_assembly_level"] == "Complete Genome"]
    n_sp: dict[str, int] = {}
    for r in comp:
        n_sp[r["gtdb_taxonomy"]] = n_sp.get(r["gtdb_taxonomy"], 0) + 1
    over = {sp: n for sp, n in n_sp.items() if n > ANI_CAP}
    planned = [s for sh in plan["shards"] for s in sh["species"]]
    assert sorted(s["species"] for s in planned) == sorted(over) and plan["species"] == len(over)
    assert len({s["id"] for s in planned}) == len(planned)
    ff = 0
    for n in over.values():  # the pairs of the farthest-first rounds, by simulating them
        left = n - 1
        for _ in range(ANI_CAP - 1):
            ff, left = ff + left, left - 1
    assert plan["pairs_farthest_first"] == ff == sum(s["pairs_farthest_first"] for s in planned)
    assert plan["pairs_quality_top"] == len(over) * ANI_CAP * (ANI_CAP - 1) // 2
    assert plan["pairs_all"] == sum(n * (n - 1) // 2 for n in over.values())
    assert plan["genomes"] == sum(over.values()) and plan["missing_paths"] == 0
    assert len(plan["shards"]) >= 2 and any(len(sh["species"]) > 1 for sh in plan["shards"])
    for sh in plan["shards"]:
        est = sum(s["est_cpu_min"] for s in sh["species"])
        assert sh["est_cpu_min"] == pytest.approx(est)
        assert est <= 1.0 or len(sh["species"]) == 1  # a species over the budget alone is its own shard
        assert sh["unit"] == f"c27-skani-{sh['k']}"
        for s in sh["species"]:
            pairs = s["pairs_farthest_first"] + ANI_CAP * (ANI_CAP - 1) // 2
            assert s["est_cpu_min"] == pytest.approx((s["genomes"] * 2 + pairs * 0.010) / 60)
    text = p.stdout
    assert "systemd-run --user --unit c27-skani-0 " in text and "nice -n 19" in text
    assert "singularity exec" in text and "/img/bio-box.sif" in text and "run-shard" in text
    # the same plan written for real, then printed by `launch --dry-run`
    subprocess.run([sys.executable, "-I", str(SKANI), "plan", *args], check=True, capture_output=True, text=True)
    assert json.loads((work / "plan.json").read_text())["shards"] == plan["shards"]
    q = subprocess.run([sys.executable, "-I", str(SKANI), "launch", str(work), "--parallel", "2", "--dry-run"],
                       check=True, capture_output=True, text=True)
    assert q.stdout.count("systemd-run --user --unit c27-skani-") == 2


FAKE_SKANI = r'''
"""Stand-in for the skani binary in tests: sketch copies names; dist and triangle -E print the world's values."""
import json, os, re, sys
from pathlib import Path
W = {tuple(k.split("|")): v for k, v in json.loads(Path(os.environ["FAKE_SKANI_WORLD"]).read_text()).items()}
ACC = re.compile(r"(GC[AF]_\d{9}\.\d+)")
def acc(p): return ACC.search(Path(p).name).group(1)
def opt(name):
    return sys.argv[sys.argv.index(name) + 1]
cmd = sys.argv[1]
with open(os.environ["FAKE_SKANI_LOG"], "a") as fh:
    fh.write(" ".join(sys.argv[1:]) + "\n")
if cmd == "sketch":
    d = Path(opt("-o")); d.mkdir(parents=True)
    for ln in Path(opt("-l")).read_text().split():
        (d / (Path(ln).name + ".sketch")).write_text(ln)
    sys.exit(0)
def rows(q, refs):
    out = []
    for r in refs:
        a, b = sorted((acc(q), acc(r)))
        v = W.get((a, b))
        if v is None or v[0] is None:
            continue
        afq, afr = (v[1], v[2]) if acc(q) == a else (v[2], v[1])
        out.append("\t".join([r, q, v[0], afr, afq, "ref name", "query name"]))
    return out
head = "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name"
if cmd == "dist":
    lines = rows(opt("-q"), Path(opt("--rl")).read_text().split())
elif cmd == "triangle":
    assert "-E" in sys.argv
    lst = Path(opt("-l")).read_text().split()
    lines = [x for i, q in enumerate(lst) for x in rows(q, lst[:i])]
else:
    sys.exit(2)
Path(opt("-o")).write_text("\n".join([head] + lines) + "\n")
'''


def test_skani_driver_run_with_stand_in_binary(tmp_path):
    """run-shard and collect against a stand-in skani that prints a known world: the collected table, read by
    c25_select.py --ani --require-pairs, gives the oracle's picks, and the driver asked skani for the
    farthest-first rounds (cap - 1 `dist` calls per species) and one `triangle` of the quality top."""
    rows = _ani_rows(seed=51)
    world = _ani_world(rows, seed=52)
    cand, paths = _driver_inputs(tmp_path, rows, ANI_CAP)
    fake = tmp_path / "bin" / "skani"
    fake.parent.mkdir()
    fake.write_text(f"#!{sys.executable}" + FAKE_SKANI)
    fake.chmod(0o755)
    wfile = tmp_path / "world.json"
    wfile.write_text(json.dumps({f"{a}|{b}": list(v) for (a, b), v in world.items()}))
    env = dict(os.environ, FAKE_SKANI_WORLD=str(wfile), FAKE_SKANI_LOG=str(tmp_path / "calls.log"))
    work = tmp_path / "work"
    subprocess.run([sys.executable, "-I", str(SKANI), "plan", str(cand), str(paths), str(work), "--cap",
                    str(ANI_CAP), "--budget-cpu-min", "100", "--threads", "2", "--image", "/img/x.sif"],
                   check=True, capture_output=True, text=True)
    plan = json.loads((work / "plan.json").read_text())
    for sh in plan["shards"]:
        subprocess.run([sys.executable, "-I", str(SKANI), "run-shard", str(work), str(sh["k"]), "--skani", str(fake)],
                       check=True, capture_output=True, text=True, env=env)
        assert json.loads((work / "shards" / str(sh["k"]) / "DONE").read_text())["species"] == len(sh["species"])
    table = tmp_path / "ani_pairs.tsv"
    subprocess.run([sys.executable, "-I", str(SKANI), "collect", str(work), str(table)], check=True,
                   capture_output=True, text=True)
    out, _ = _run_cap(tmp_path / "s", rows, _phase1_lines(rows), ANI_CAP, ["--ani", str(table), "--require-pairs"])
    assert _picked(out) == _ff_oracle(rows, world, ANI_CAP)
    calls = (tmp_path / "calls.log").read_text().splitlines()
    n_over = plan["species"]
    assert sum(c.startswith("sketch ") for c in calls) == n_over
    assert sum(c.startswith("dist ") for c in calls) == n_over * (ANI_CAP - 1)
    assert sum(c.startswith("triangle ") for c in calls) == n_over
    t = _read(table)
    assert len(t) == len({(r["genome_a"], r["genome_b"]) for r in t})
    j = json.loads((out / "selection.json").read_text())["ani"]
    assert j["consulted_pairs"] == plan["pairs_farthest_first"] and j["consulted_absent"] == 0
    # a second run of a finished shard asks skani nothing (resumable by markers)
    subprocess.run([sys.executable, "-I", str(SKANI), "run-shard", str(work), "0", "--skani", str(fake)],
                   check=True, capture_output=True, text=True, env=env)
    assert len((tmp_path / "calls.log").read_text().splitlines()) == len(calls)


# --- the fetch: Datasets zip members kept compressed, verified against md5sum.txt --------------------------------

def _load(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(path.parent))
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.path.remove(str(path.parent))
    return mod


def _fasta(acc: str, n: int, seed: int) -> bytes:
    rng = random.Random(seed)
    out = []
    for k in range(n):
        seq = "".join(rng.choice("ACGT") for _ in range(rng.randrange(100, 3000)))
        out.append(f">NZ_CP{seed:06d}.{k} {acc} contig {k}\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)))
    return ("\n".join(out) + "\n").encode()


def _bp(data: bytes) -> int:
    return sum(len(x) for x in data.split(b"\n") if x and not x.startswith(b">"))


def _datasets_zip(path: Path, genomes: dict[str, bytes], bad=frozenset(), stored=frozenset()) -> None:
    md5 = []
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("README.md", "NCBI Datasets\n")
        for acc, data in genomes.items():
            name = f"ncbi_dataset/data/{acc}/{acc}_ASM{acc[4:8]}v1_genomic.fna"
            info = zipfile.ZipInfo(name, (2026, 11, 12, 0, 0, 0))
            info.compress_type = zipfile.ZIP_STORED if acc in stored else zipfile.ZIP_DEFLATED
            z.writestr(info, data)
            md5.append(f"{hashlib.md5(b'x' if acc in bad else data).hexdigest()}  {name}")
        z.writestr("md5sum.txt", "\n".join(md5) + "\n")


def test_fetch_keeps_members_compressed_and_verified(tmp_path):
    import gzip as gz
    f = _load(FETCH)
    genomes = {a: _fasta(a, 1 + i % 3, i) for i, a in
               enumerate(["GCF_000000010.1", "GCA_000000030.1", "GCF_000000040.1", "GCF_000000050.1"])}
    z = tmp_path / "batch.zip"
    _datasets_zip(z, genomes, bad={"GCF_000000040.1"}, stored={"GCA_000000030.1"})
    out = tmp_path / "out"
    out.mkdir()
    got = f.extract_zip_gz(z, set(genomes) - {"GCF_000000050.1"}, out)
    assert sorted(got) == ["GCA_000000030.1", "GCF_000000010.1", "GCF_000000040.1"]
    assert sorted(a for a, g in got.items() if g["ok"]) == ["GCA_000000030.1", "GCF_000000010.1"]
    assert not (out / "GCF_000000040.1.fna.gz").exists() and not (out / "GCF_000000050.1.fna.gz").exists()
    assert sorted(p.name for p in out.iterdir()) == ["GCA_000000030.1.fna.gz", "GCF_000000010.1.fna.gz"]
    for acc in ("GCA_000000030.1", "GCF_000000010.1"):
        p = out / f"{acc}.fna.gz"
        data = genomes[acc]
        assert gz.decompress(p.read_bytes()) == data
        sc = got[acc]["scan"]
        assert sc["sha256"] == _sha(data) and sc["md5"] == hashlib.md5(data).hexdigest()
        assert sc["bytes"] == len(data) and sc["gz_bytes"] == p.stat().st_size
        assert sc["gz_sha256"] == _sha(p.read_bytes())
        assert sc["records"] == data.count(b">") and sc["bp"] == _bp(data)
    out2 = tmp_path / "out2"  # deterministic: the same zip gives the same .fna.gz bytes
    out2.mkdir()
    f.extract_zip_gz(z, {"GCF_000000010.1"}, out2)
    assert (out2 / "GCF_000000010.1.fna.gz").read_bytes() == (out / "GCF_000000010.1.fna.gz").read_bytes()


def test_fetch_batch_retries_lists_failures_and_resumes(tmp_path):
    f = _load(FETCH)
    genomes = {f"GCF_0000001{i:02d}.1": _fasta(f"GCF_0000001{i:02d}.1", 1 + i % 3, 10 + i) for i in range(7)}
    absent = "GCF_000000103.1"   # never served
    flaky = "GCF_000000105.1"    # a bad md5 in its batch, good alone
    calls = []

    class FakeFetcher:
        def __init__(self):
            self.R = {"requests": 0, "retries": 0, "bytes_downloaded": 0, "refusal_codes": {}}

        def line(self, *a):
            pass

        def get(self, url, dest=None):
            accs = url.split("/accession/")[1].split("/download")[0].split(",")
            calls.append(accs)
            self.R["requests"] += 1
            _datasets_zip(dest, {a: genomes[a] for a in accs if a != absent}, bad={flaky} if len(accs) > 1 else set())
            self.R["bytes_downloaded"] += dest.stat().st_size

    tsv = tmp_path / "batch_B01.tsv"
    tsv.write_text("accession\tgenome_size\n" + "".join(f"{a}\t{_bp(d)}\n" for a, d in sorted(genomes.items())))
    out = tmp_path / "B01"
    res = f.run_batch(FakeFetcher(), "B01", tsv, out, batch=3, stop=tmp_path / "STOP")
    assert (res["genomes_wanted"], res["genomes_ok"], res["genomes_failed"]) == (7, 6, 1)
    r = json.loads((out / "fetch_result.json").read_text())
    assert list(r["failed"]) == [absent]
    assert r["md5_refetched"] == [flaky] and r["fallback_single"] == []
    lock = json.loads((out / "lock.json").read_text())["genomes"]
    assert sorted(lock) == sorted(set(genomes) - {absent})
    assert all(v["bp"] == v["gtdb_genome_size"] for v in lock.values())
    assert all(v["file_name"] == f"{a}.fna.gz" and v["md5_source"] == "datasets-zip md5sum.txt"
               for a, v in lock.items())
    assert [len(c) for c in calls] == [3, 3, 1, 1, 1]  # batches of 3, 3, 1; the absent and the flaky alone
    n = len(calls)
    f.run_batch(FakeFetcher(), "B01", tsv, out, batch=3, stop=tmp_path / "STOP")
    assert calls[n:] == [[absent], [absent]]  # resumed: only the failed genome is asked again (batch, then alone)
    (tmp_path / "STOP").write_text("")
    (out / "lock.json").unlink()
    with pytest.raises(f.Stopped):
        f.run_batch(FakeFetcher(), "B01", tsv, out, batch=3, stop=tmp_path / "STOP")


def test_fetch_split_into_batches(tmp_path):
    f = _load(FETCH)
    rng = random.Random(61)
    rows = [{"accession": f"GCF_{rng.randrange(10**9):09d}.1", "genome_size": str(rng.randrange(10**6, 9 * 10**6)),
             "fetch_reason": "ani_candidate"} for _ in range(300)]
    retry = [{"accession": "GCA_900637225.1", "genome_size": "2000000", "fetch_reason": "phase1_not_served_retry"}]
    allrows = rows + retry
    random.Random(1).shuffle(allrows)
    src = tmp_path / "ani_fetch.tsv"
    src.write_text("accession\tgenome_size\tfetch_reason\n" + "".join(
        f"{r['accession']}\t{r['genome_size']}\t{r['fetch_reason']}\n" for r in allrows))
    out = tmp_path / "batches"
    plan = f.split(src, out, batch_bp=10**8)
    assert [r["accession"] for r in _read(out / "batch_B00.tsv")] == ["GCA_900637225.1"]
    where = {r["accession"]: k for k in plan["batches"] for r in _read(out / f"batch_{k}.tsv")}
    later = sorted(rows, key=lambda r: r["accession"])
    total = sum(int(r["genome_size"]) for r in later)
    n = max(1, round(total / 10**8))
    assert len(plan["batches"]) == n + 1
    c = 0
    for r in later:  # the genome whose predecessors sum to c bp goes to batch 1 + floor(n c / total)
        assert where[r["accession"]] == f"B{1 + (n * c) // total:02d}"
        c += int(r["genome_size"])
    for k in plan["batches"]:
        accs = [r["accession"] for r in _read(out / f"batch_{k}.tsv")]
        assert accs == sorted(accs)
