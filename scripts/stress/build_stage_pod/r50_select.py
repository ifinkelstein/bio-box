#!/usr/bin/env python3
"""Genome selection for the 50 Gbp build (host python3, stdlib only).

    python3 r50_select.py GTDB_DIR LADDER_SELECTION_DIR LADDER_GENOMES_DIR OUT_DIR [--stage-bp 10000000000]

The import ladder's rule (planning/measurements/2026-11-10_import_ladder/scripts/select.py), re-implemented
without DuckDB: GTDB R226 bac120 + ar53 metadata; genomes with ncbi_assembly_level = 'Complete Genome';
one per species (gtdb_taxonomy): the highest checkm2_completeness - 5 * checkm2_contamination (a value that
does not parse sorts last), ties by GTDB accession (with its RS_/GB_ prefix). The ladder's P1 and P2 are
recomputed (species shuffled with random.Random(20261109) over the sorted species list; first 250, next
2,250) and must equal the committed stage_P1.tsv and stage_P2.tsv byte for byte, or the script exits 1.

New genomes: every selected accession that is not in the ladder's genomes/P1/lock.json or genomes/P2/lock.json
(the 2,497 already fetched and reused in place), sorted by accession. Stages (a pure function of that
sorted list and its GTDB genome_size): n = max(1, round(total_bp / stage_bp)); the genome whose preceding
genomes sum to c bp goes to stage P(3 + floor(n * c / total_bp)).

Writes OUT_DIR/all.tsv (every selected genome), OUT_DIR/new.tsv, OUT_DIR/stage_P<k>.tsv (the ladder's
columns: accession, genome_size, contig_count, protein_count, translation_table, gtdb_taxonomy,
gtdb_accession) and OUT_DIR/selection.json; prints selection.json.
"""

import argparse
import csv
import gzip
import hashlib
import json
import random
import sys
from pathlib import Path

COLS = ["accession", "genome_size", "contig_count", "protein_count", "translation_table", "gtdb_taxonomy",
        "gtdb_accession"]
FILES = ["bac120_metadata_r226.tsv.gz", "ar53_metadata_r226.tsv.gz"]


def num(s: str, kind=float):
    try:
        return kind(s)
    except (TypeError, ValueError):
        return None


def cell(v) -> str:
    return "" if v is None else str(v)


def write_tsv(path: Path, rows: list[dict]) -> None:
    tmp = path.with_name(f".{path.name}.tmp")
    with open(tmp, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(cell(r[c]) for c in COLS) + "\n")
    tmp.replace(path)


def summary(rows: list[dict]) -> dict:
    tt: dict[str, int] = {}
    for r in rows:
        tt[str(r["translation_table"])] = tt.get(str(r["translation_table"]), 0) + 1
    return {"genomes": len(rows), "bp": sum(r["genome_size"] for r in rows),
            "contigs": sum(r["contig_count"] for r in rows),
            "proteins": sum(r["protein_count"] or 0 for r in rows),
            "proteins_missing": sum(r["protein_count"] is None for r in rows),
            "archaea": sum(r["gtdb_taxonomy"].startswith("d__Archaea") for r in rows),
            "translation_tables": dict(sorted(tt.items())),
            "max_genome_bp": max(r["genome_size"] for r in rows) if rows else None,
            "min_genome_bp": min(r["genome_size"] for r in rows) if rows else None,
            "first": rows[0]["accession"] if rows else None, "last": rows[-1]["accession"] if rows else None}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("gtdb", type=Path)
    ap.add_argument("ladder_selection", type=Path)
    ap.add_argument("ladder_genomes", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--stage-bp", type=int, default=10_000_000_000)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    best: dict[str, tuple] = {}
    n_complete = n_complete_rep = 0
    sha = {}
    for name in FILES:
        p = a.gtdb / name
        sha[name] = hashlib.sha256(p.read_bytes()).hexdigest()
        with gzip.open(p, "rt", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["ncbi_assembly_level"] != "Complete Genome":
                    continue
                n_complete += 1
                n_complete_rep += r["gtdb_representative"] == "t"
                comp, cont = num(r["checkm2_completeness"]), num(r["checkm2_contamination"])
                score = None if comp is None or cont is None else comp - 5 * cont
                gacc = r["accession"]
                row = {"gtdb_accession": gacc, "accession": gacc[3:] if gacc[:3] in ("RS_", "GB_") else gacc,
                       "gtdb_taxonomy": r["gtdb_taxonomy"], "genome_size": int(r["genome_size"]),
                       "contig_count": int(r["contig_count"]), "protein_count": num(r["protein_count"], int),
                       "translation_table": num(r["ncbi_translation_table"], int),
                       "gtdb_representative": r["gtdb_representative"] == "t"}
                # sort key: score descending, missing last, then GTDB accession ascending
                key = (0 if score is not None else 1, -(score or 0.0), gacc)
                cur = best.get(row["gtdb_taxonomy"])
                if cur is None or key < cur[0]:
                    best[row["gtdb_taxonomy"]] = (key, row)
    picked = {sp: v[1] for sp, v in best.items()}

    # the ladder's P1 and P2, recomputed and compared byte for byte
    species = sorted(picked)
    order = species[:]
    random.Random(20261109).shuffle(order)
    ladder_ok = {}
    for stage, sp in (("P1", order[:250]), ("P2", order[250:2500])):
        mine = a.out / f".ladder_{stage}.tsv"
        write_tsv(mine, sorted((picked[s] for s in sp), key=lambda r: r["accession"]))
        theirs = a.ladder_selection / f"stage_{stage}.tsv"
        ladder_ok[stage] = mine.read_bytes() == theirs.read_bytes()
        mine.unlink()
    if not all(ladder_ok.values()):
        print(json.dumps({"error": "ladder selection not reproduced", "equal": ladder_ok}), file=sys.stderr)
        return 1

    present = set()
    for stage in ("P1", "P2"):
        present |= set(json.loads((a.ladder_genomes / stage / "lock.json").read_text())["genomes"])
    rows = sorted(picked.values(), key=lambda r: r["accession"])
    by_acc = {r["accession"]: r for r in rows}
    assert len(by_acc) == len(rows), "accession selected for two species"
    assert present <= set(by_acc), "a fetched ladder genome is not in the selection"
    new = [r for r in rows if r["accession"] not in present]
    ladder_unfetched = sorted({picked[s]["accession"] for s in order[:2500]} - present)

    total = sum(r["genome_size"] for r in new)
    n = max(1, round(total / a.stage_bp))
    stages: dict[str, list[dict]] = {}
    c = 0
    for r in new:
        stages.setdefault(f"P{3 + (n * c) // total}", []).append(r)
        c += r["genome_size"]

    write_tsv(a.out / "all.tsv", rows)
    write_tsv(a.out / "new.tsv", new)
    for s, rs in stages.items():
        write_tsv(a.out / f"stage_{s}.tsv", rs)
    out = {"gtdb_release": "R226", "gtdb_files_sha256": sha,
           "rule": "ncbi_assembly_level = 'Complete Genome'; one per gtdb_taxonomy by checkm2 completeness - 5 x "
                   "contamination desc (missing last), ties by GTDB accession",
           "complete_genomes": n_complete, "complete_genomes_gtdb_representative": n_complete_rep,
           "species_with_complete_genome": len(rows),
           "selected_that_are_gtdb_representative": sum(r["gtdb_representative"] for r in rows),
           "ladder_P1_P2_reproduced": ladder_ok, "ladder_present": len(present),
           "ladder_selected_not_fetched": ladder_unfetched,
           "all": summary(rows), "new": summary(new),
           "split": {"stage_bp": a.stage_bp, "n_stages": n, "total_bp": total,
                     "rule": "stage P(3 + floor(n * bp_before / total_bp)) over new genomes sorted by accession"},
           "stages": {s: summary(rs) for s, rs in stages.items()}}
    p = a.out / "selection.json"
    p.with_name(".selection.json.tmp").write_text(json.dumps(out, indent=1) + "\n")
    p.with_name(".selection.json.tmp").replace(p)
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
