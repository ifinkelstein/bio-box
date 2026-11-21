#!/usr/bin/env python3
"""The per-accession list of what S2C phase 2's selection left out (card C-51, DEC-59; host python3, stdlib only).

    python3 -I c51_deferred_list.py GTDB_DIR PHASE1_ALL_TSV STAGES_TSV OUT_DIR [--cap 20]
            [--species-table ANI_SPECIES_TSV] [--not-served-candidates FILE] [--not-served-picks FILE]

Card C-27 (DEC-31, DEC-34) kept at most CAP complete genomes in each GTDB species: phase 1's genome and CAP - 1 picks,
which stages C1 to C6 hold. This script writes the other complete genomes of every species with more than CAP.

GTDB_DIR holds bac120_metadata_r226.tsv.gz and ar53_metadata_r226.tsv.gz (the tables card C-27 read; sources.yaml
`gtdb_metadata`). A genome is complete when ncbi_assembly_level = 'Complete Genome', as in c25_select.py. PHASE1_ALL_TSV
is phase 1's all.tsv (columns accession, gtdb_taxonomy, ...); STAGES_TSV is card C-36's stages.tsv (accession, stage,
..., the picks of C1 to C6). --not-served-candidates: accessions NCBI did not serve that were candidates of the pick
(card C-27's 26); they stay in the list, marked not_served = y. --not-served-picks: picks NCBI did not serve (card
C-25's 8, in species at or under the cap); they are not in the list and are counted by the stage of their species.
--species-table: ani_species.tsv of card C-27 (columns species, complete_genomes); every species over the cap is
checked against it.

Writes, in OUT_DIR, with no timestamp, path or input-order dependence:
  left_out.tsv       accession, species, assembly_level, not_served; sorted by accession
  species_check.tsv  per species over the cap: species, table_complete, complete, phase1_in, picks_in, left,
                     expected_left (table_complete - cap), stage, differs (y), why; sorted by species
  summary.json       counts, by stage, the input digests
"""

import argparse
import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path

FILES = ["bac120_metadata_r226.tsv.gz", "ar53_metadata_r226.tsv.gz"]
LEFT_COLS = ["accession", "species", "assembly_level", "not_served"]
CHECK_COLS = ["species", "table_complete", "complete", "phase1_in", "picks_in", "left", "expected_left", "stage",
              "differs", "why"]


def die(msg: str):
    print(f"c51_deferred_list.py: error: {msg}", file=sys.stderr)
    raise SystemExit(2)


def bare(acc: str) -> str:
    return acc[3:] if acc[:3] in ("RS_", "GB_") else acc


def species_name(tax: str) -> str:
    return tax.rsplit(";", 1)[-1]


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_list(path: Path | None) -> set[str]:
    return set(path.read_text().split()) if path is not None else set()


def put(path: Path, text: str) -> str:
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(text)
    tmp.replace(path)
    return hashlib.sha256(text.encode()).hexdigest()


def tsv(rows: list[dict], cols: list[str]) -> str:
    return "\t".join(cols) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows)


def complete_genomes(gtdb: Path, sha: dict) -> dict[str, tuple[str, str]]:
    """bare accession -> (gtdb_taxonomy, assembly level), for every complete genome; each file's sha256 in `sha`."""
    out: dict[str, tuple[str, str]] = {}
    for name in FILES:
        p = gtdb / name
        sha[name] = hashlib.sha256(p.read_bytes()).hexdigest()
        with gzip.open(p, "rt", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["ncbi_assembly_level"] != "Complete Genome":
                    continue
                acc = bare(r["accession"])
                if acc in out:
                    die(f"{name}: {acc} appears twice among the complete genomes")
                out[acc] = (r["gtdb_taxonomy"], r["ncbi_assembly_level"])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("gtdb", type=Path)
    ap.add_argument("phase1", type=Path)
    ap.add_argument("stages", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--cap", type=int, default=20)
    ap.add_argument("--species-table", type=Path)
    ap.add_argument("--not-served-candidates", type=Path)
    ap.add_argument("--not-served-picks", type=Path)
    a = ap.parse_args()
    cap = a.cap

    sha: dict[str, str] = {}
    comp = complete_genomes(a.gtdb, sha)
    p1 = {r["accession"] for r in read_tsv(a.phase1)}
    pick_stage: dict[str, str] = {}
    for r in read_tsv(a.stages):
        if r["accession"] in pick_stage:
            die(f"{a.stages}: {r['accession']} is picked twice")
        pick_stage[r["accession"]] = r["stage"]
    ns_cand = read_list(a.not_served_candidates)
    ns_picks = read_list(a.not_served_picks)
    for acc in sorted(ns_cand):
        if acc in pick_stage or acc in p1:
            die(f"{acc} is listed as a not-served candidate and is also a pick or phase 1's genome")
    table: dict[str, int] = {}
    if a.species_table is not None:
        table = {r["species"]: int(r["complete_genomes"]) for r in read_tsv(a.species_table)}

    by_sp: dict[str, list[str]] = {}
    for acc, (tax, _) in comp.items():
        by_sp.setdefault(tax, []).append(acc)
    over = {tax: sorted(accs) for tax, accs in by_sp.items() if len(accs) > cap}
    for tax in sorted(set(table) - set(over)):
        die(f"{tax}: in the species table, but it has {len(by_sp.get(tax, []))} complete genomes in the metadata "
            f"(cap {cap})")
    for tax in sorted(set(over) - set(table)) if table else []:
        die(f"{tax}: {len(over[tax])} complete genomes in the metadata, over the cap, and not in the species table")

    stage_of_sp: dict[str, set[str]] = {}
    picks_not_complete = []
    for acc, st in pick_stage.items():
        if acc not in comp:
            picks_not_complete.append(acc)
            continue
        stage_of_sp.setdefault(comp[acc][0], set()).add(st)
    for tax, sts in stage_of_sp.items():
        if len(sts) > 1:
            die(f"{tax}: picks in more than one stage: {sorted(sts)}")
    stage = {tax: next(iter(sts)) for tax, sts in stage_of_sp.items()}

    left, checks = [], []
    for tax in sorted(over):
        accs = over[tax]
        p1_in = sum(x in p1 for x in accs)
        picks_in = sum(x in pick_stage for x in accs)
        rest = [x for x in accs if x not in p1 and x not in pick_stage]
        for x in rest:
            left.append({"accession": x, "species": species_name(tax), "assembly_level": comp[x][1],
                         "not_served": "y" if x in ns_cand else ""})
        tc = table.get(tax, len(accs))
        expected = tc - cap
        why = []
        if tc != len(accs):
            why.append(f"the metadata has {len(accs)} complete genomes, the table {tc}")
        if p1_in != 1:
            why.append("phase-1 genome not among the complete genomes" if p1_in == 0
                       else f"{p1_in} phase-1 genomes among the complete genomes")
        if picks_in != cap - 1 and p1_in == 1:
            why.append(f"{picks_in} picks in stages C1 to C6, not {cap - 1}")
        differs = len(rest) != expected
        if differs and not why:
            why.append("unexplained")
        checks.append({"species": tax, "table_complete": tc, "complete": len(accs), "phase1_in": p1_in,
                       "picks_in": picks_in, "left": len(rest), "expected_left": expected,
                       "stage": stage.get(tax, ""), "differs": "y" if differs else "",
                       "why": "; ".join(why) if differs else ""})
    left.sort(key=lambda r: r["accession"])

    by_stage: dict[str, int] = {}
    species_by_stage: dict[str, int] = {}
    for c in checks:
        k = c["stage"] or "none"
        by_stage[k] = by_stage.get(k, 0) + c["left"]
        species_by_stage[k] = species_by_stage.get(k, 0) + 1
    ns_in_left = {r["accession"] for r in left if r["not_served"]}
    ns_cand_by_stage: dict[str, int] = {}
    for acc in ns_in_left:
        k = stage.get(comp[acc][0], "none")
        ns_cand_by_stage[k] = ns_cand_by_stage.get(k, 0) + 1
    ns_pick_by_stage: dict[str, int] = {}
    for acc in sorted(ns_picks):
        k = stage.get(comp[acc][0], "none") if acc in comp else "none"
        ns_pick_by_stage[k] = ns_pick_by_stage.get(k, 0) + 1

    a.out.mkdir(parents=True, exist_ok=True)
    digests = {"left_out.tsv": put(a.out / "left_out.tsv", tsv(left, LEFT_COLS)),
               "species_check.tsv": put(a.out / "species_check.tsv", tsv(checks, CHECK_COLS))}
    summary = {
        "cap": cap, "gtdb_files_sha256": sha, "complete_genomes": len(comp), "species_with_complete_genome": len(by_sp),
        "species_over_cap": len(over), "complete_genomes_in_species_over_cap": sum(len(v) for v in over.values()),
        "phase1_genomes": len(p1), "picks": len(pick_stage), "picks_not_complete_in_metadata": sorted(picks_not_complete),
        "left_out": len(left), "species_differ": sum(c["differs"] == "y" for c in checks),
        "not_served_in_left_out": len(ns_in_left), "not_served_candidates_listed": len(ns_cand),
        "not_served_candidates_by_stage": dict(sorted(ns_cand_by_stage.items())),
        "not_served_picks_listed": len(ns_picks), "not_served_picks_by_stage": dict(sorted(ns_pick_by_stage.items())),
        "by_stage": dict(sorted(by_stage.items())), "species_over_cap_by_stage": dict(sorted(species_by_stage.items())),
        "inputs_sha256": {"phase1": hashlib.sha256(a.phase1.read_bytes()).hexdigest(),
                          "stages": hashlib.sha256(a.stages.read_bytes()).hexdigest()},
        "outputs_sha256": digests,
    }
    text = json.dumps(summary, indent=1, sort_keys=True) + "\n"
    put(a.out / "summary.json", text)
    print(text, end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
