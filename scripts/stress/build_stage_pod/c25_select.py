#!/usr/bin/env python3
"""S2C selection for card C-25 (host python3, stdlib only): at most CAP complete genomes per GTDB species.

    python3 -I c25_select.py GTDB_DIR PHASE1_ALL_TSV OUT_DIR [--cap 20] [--phase1-not-served ACC,ACC,...]

Rule, by quality order only (the plan's ANI-diversity step needs skani, which is not in the image):
GTDB R226 bac120 + ar53 metadata, read by r50_select.complete_genomes (phase 1's reader); a genome is complete when
ncbi_assembly_level = 'Complete Genome' (no other field is used); within each gtdb_taxonomy the genomes are ranked
by checkm2_completeness - 5 * checkm2_contamination descending, a value that does not parse last, every tie by
GTDB accession (with its RS_/GB_ prefix; this orders as the bare accession does, since GB_ < RS_ and GCA < GCF);
ranks 1 to CAP are selected. Rank 1 is phase 1's rule (r50_select.py), and the whole rule is the one that gave the
plan's 25,458 (scripts/stress/gtdb_counts.py, S2C_complete_le20_per_species).

PHASE1_ALL_TSV is phase 1's all.tsv (11,973 genomes); --phase1-not-served lists the phase-1 accessions NCBI did not
serve. Writes, sorted by accession, OUT_DIR/s2c_selected.tsv and OUT_DIR/remainder.tsv (selected minus phase 1),
with phase 1's columns then species_rank, species_complete_genomes, quality_score, checkm2_completeness,
checkm2_contamination, gtdb_representative; OUT_DIR/by_phylum.tsv; OUT_DIR/selection.json (counts, the phase-1
comparison, what an ANI step could change, file digests); prints selection.json.
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python3 -I leaves the script's directory off the path
from r50_select import COLS, FILES, cell, complete_genomes, summary  # noqa: E402

EXTRA = ["species_rank", "species_complete_genomes", "quality_score", "checkm2_completeness",
         "checkm2_contamination", "gtdb_representative"]
BINS = [(1, 1), (2, 2), (3, 5), (6, 10), (11, 20), (21, 50), (51, 100), (101, None)]
PHYLUM_COLS = ["phylum", "species", "selected_genomes", "selected_bp", "remainder_genomes", "remainder_bp",
               "remainder_proteins"]


def rank_species(items: list[tuple]) -> dict[str, list[dict]]:
    """(key, row) pairs -> species -> rows in quality order, each annotated with species_rank and the species size."""
    by_sp: dict[str, list[tuple]] = {}
    for key, row in items:
        by_sp.setdefault(row["gtdb_taxonomy"], []).append((key, row))
    out = {}
    for sp in sorted(by_sp):
        ranked = [row for _, row in sorted(by_sp[sp], key=lambda kr: kr[0])]
        for i, row in enumerate(ranked, 1):
            row["species_rank"], row["species_complete_genomes"] = i, len(ranked)
        out[sp] = ranked
    return out


def fmt(row: dict, col: str) -> str:
    v = row[col]
    if col == "gtdb_representative":
        return "t" if v else "f"
    return cell(v)


def tsv_text(rows: list[dict], cols: list[str]) -> str:
    return "\t".join(cols) + "\n" + "".join("\t".join(fmt(r, c) for c in cols) + "\n" for r in rows)


def put(path: Path, text: str) -> str:
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(text)
    tmp.replace(path)
    return hashlib.sha256(text.encode()).hexdigest()


def same_score(a: dict, b: dict) -> bool:
    return a["quality_score"] == b["quality_score"]  # None == None: both unparsed, ordered by accession alone


def bin_label(n: int) -> str:
    for lo, hi in BINS:
        if n >= lo and (hi is None or n <= hi):
            return str(lo) if lo == hi else (f"{lo}+" if hi is None else f"{lo}-{hi}")
    raise ValueError(n)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("gtdb", type=Path)
    ap.add_argument("phase1", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--cap", type=int, default=20)
    ap.add_argument("--phase1-not-served", default="")
    a = ap.parse_args()
    cap = a.cap
    a.out.mkdir(parents=True, exist_ok=True)

    p1_bytes = a.phase1.read_bytes()
    with open(a.phase1, newline="") as fh:
        p1 = [r["accession"] for r in csv.DictReader(fh, delimiter="\t")]
    p1_set = set(p1)
    not_served = sorted(x for x in a.phase1_not_served.split(",") if x)

    sha: dict[str, str] = {}
    levels: dict[str, str] = {}
    items = list(complete_genomes(a.gtdb, sha, watch=p1_set | set(not_served), levels=levels))
    header_ani = []
    for name in FILES:
        with gzip.open(a.gtdb / name, "rt") as fh:
            header_ani += [c for c in fh.readline().rstrip("\n").split("\t") if "ani" in c.lower().split("_")]
    unparsed = sum(r["quality_score"] is None for _, r in items)
    nonfinite = sum(r["quality_score"] is not None and not math.isfinite(r["quality_score"]) for _, r in items)

    ranked = rank_species(items)
    by_acc = {r["accession"]: r for rows in ranked.values() for r in rows}
    assert len(by_acc) == len(items), "an accession appears twice among the complete genomes"
    selected = sorted((r for rows in ranked.values() for r in rows[:cap]), key=lambda r: r["accession"])
    rank1 = sorted((rows[0] for rows in ranked.values()), key=lambda r: r["accession"])
    remainder = [r for r in selected if r["accession"] not in p1_set]

    # phase 1 against the set
    missing = []
    for acc in sorted(p1_set):
        r = by_acc.get(acc)
        if r is not None and r["species_rank"] <= cap:
            continue
        if acc not in levels:
            why = "not in the metadata"
        elif levels[acc] != "Complete Genome":
            why = f"not complete in the metadata (ncbi_assembly_level {levels[acc]})"
        else:
            why = f"rank {r['species_rank']} of {r['species_complete_genomes']} complete genomes of its species, " \
                  "beyond the cap"
        missing.append({"accession": acc, "reason": why})
    ns = [{"accession": acc, "species_rank": by_acc[acc]["species_rank"] if acc in by_acc else None,
           "species_complete_genomes": by_acc[acc]["species_complete_genomes"] if acc in by_acc else None}
          for acc in not_served]

    # what an ANI-diversity step could change: only species with more than cap complete genomes
    over = {sp: rows for sp, rows in ranked.items() if len(rows) > cap}
    lo = hi = 0
    tie_cut = cut_top = decided_by_acc = 0
    for sp, rows in ranked.items():
        best, rest = rows[0]["genome_size"], sorted(r["genome_size"] for r in rows[1:])
        k = min(cap - 1, len(rest))
        lo += best + sum(rest[:k])
        hi += best + sum(rest[len(rest) - k:])
        if sp in over:
            if same_score(rows[cap - 1], rows[cap]):
                tie_cut += 1
                decided_by_acc += sum(same_score(r, rows[cap - 1]) for r in rows[:cap])
            cut_top += same_score(rows[0], rows[cap])
    p1_bp = sum(r["genome_size"] for r in rank1)
    rem_over = [r for r in remainder if r["species_complete_genomes"] > cap]

    bins: dict[str, dict] = {}
    for lo_, hi_ in BINS:
        bins[bin_label(lo_)] = {"bin": bin_label(lo_), "species": 0, "selected_genomes": 0, "selected_bp": 0,
                                "remainder_genomes": 0, "remainder_bp": 0}
    for sp, rows in ranked.items():
        b = bins[bin_label(len(rows))]
        b["species"] += 1
        for r in rows[:cap]:
            b["selected_genomes"] += 1
            b["selected_bp"] += r["genome_size"]
            if r["accession"] not in p1_set:
                b["remainder_genomes"] += 1
                b["remainder_bp"] += r["genome_size"]

    phyla: dict[str, dict] = {}
    for sp, rows in ranked.items():
        ph = sp.split(";")[1] if ";" in sp else sp
        p = phyla.setdefault(ph, dict.fromkeys(PHYLUM_COLS, 0) | {"phylum": ph})
        p["species"] += 1
        for r in rows[:cap]:
            p["selected_genomes"] += 1
            p["selected_bp"] += r["genome_size"]
            if r["accession"] not in p1_set:
                p["remainder_genomes"] += 1
                p["remainder_bp"] += r["genome_size"]
                p["remainder_proteins"] += r["protein_count"] or 0
    phylum_rows = sorted(phyla.values(), key=lambda p: (-p["remainder_genomes"], -p["selected_genomes"], p["phylum"]))

    cols = COLS + EXTRA
    digests = {"s2c_selected.tsv": put(a.out / "s2c_selected.tsv", tsv_text(selected, cols)),
               "remainder.tsv": put(a.out / "remainder.tsv", tsv_text(remainder, cols)),
               "by_phylum.tsv": put(a.out / "by_phylum.tsv", tsv_text(phylum_rows, PHYLUM_COLS))}
    out = {
        "gtdb_release": "R226", "gtdb_files_sha256": sha, "cap": cap,
        "rule": "ncbi_assembly_level = 'Complete Genome'; per gtdb_taxonomy ranked by checkm2_completeness - 5 x "
                "checkm2_contamination desc (unparsed last), ties by GTDB accession; ranks 1..cap selected; no ANI step",
        "metadata_columns_named_ani": header_ani,
        "complete_genomes": len(items), "species_with_complete_genome": len(ranked),
        "complete_genomes_checkm2_unparsed": unparsed, "complete_genomes_checkm2_nonfinite": nonfinite,
        "selected": summary(selected),
        "selected_gtdb_representative": sum(r["gtdb_representative"] for r in selected),
        "phase1": {"file": a.phase1.name, "sha256": hashlib.sha256(p1_bytes).hexdigest(), "genomes": len(p1),
                   "rank1_equals_file": tsv_text(rank1, COLS).encode() == p1_bytes,
                   "rank1_bp": p1_bp, "not_in_selection": missing, "not_served": ns},
        "remainder": summary(remainder),
        "ani_step": {
            "species_over_cap": len(over),
            "genomes_in_species_over_cap": sum(len(v) for v in over.values()),
            "selected_in_species_over_cap": cap * len(over),
            "max_selected_differ_keep_rank1": sum(min(cap - 1, len(v) - cap) for v in over.values()),
            "max_selected_differ_any": sum(min(cap, len(v) - cap) for v in over.values()),
            "remainder_in_species_over_cap": {"genomes": len(rem_over), "bp": sum(r["genome_size"] for r in rem_over)},
            "bp_bounds_keep_rank1": {"min": lo, "max": hi, "quality_order": sum(r["genome_size"] for r in selected),
                                     "remainder_min": lo - p1_bp, "remainder_max": hi - p1_bp},
            "species_score_tie_at_cut": tie_cut,
            "selected_decided_by_accession_at_cut": decided_by_acc,
            "species_rank1_score_equals_rank_cap_plus_1": cut_top},
        "by_species_size": list(bins.values()),
        "outputs_sha256": digests,
    }
    text = json.dumps(out, indent=1) + "\n"
    put(a.out / "selection.json", text)
    print(text, end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
