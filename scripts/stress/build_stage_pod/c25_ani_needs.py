#!/usr/bin/env python3
"""What a skani step would need for card C-25 (1b) (host python3, stdlib only; counts only, nothing fetched).

    python3 -I c25_ani_needs.py GTDB_DIR [--cap 20]

Reads the GTDB R226 metadata with phase 1's reader (r50_select.complete_genomes, the same one c25_select.py
uses) and, for the species with more than CAP complete genomes, prints as JSON: the genomes and bp a
within-species ANI computation would have to read; how many of them phase 1 already holds (rank 1 of each
species); the pairs an all-pairs computation needs, sum n(n-1)/2; the pairs a greedy farthest-first pick of
CAP genomes from rank 1 needs, sum over the CAP-1 later picks of the distances from each pick to the genomes
not yet picked, sum (CAP-1) n - CAP(CAP-1)/2; whether the species representative is among the complete
genomes; and the ten largest species. Check: the counts equal step (1)'s selection.json
(genomes_in_species_over_cap, remainder_in_species_over_cap).
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python3 -I leaves the script's directory off the path
from r50_select import complete_genomes  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("gtdb", type=Path)
    ap.add_argument("--cap", type=int, default=20)
    a = ap.parse_args()
    cap = a.cap
    sha: dict[str, str] = {}
    by_sp: dict[str, list] = {}
    for key, row in complete_genomes(a.gtdb, sha):
        by_sp.setdefault(row["gtdb_taxonomy"], []).append((key, row))
    over = {}
    for sp in sorted(by_sp):
        if len(by_sp[sp]) > cap:
            over[sp] = [r for _, r in sorted(by_sp[sp], key=lambda kr: kr[0])]

    def sp_stats(sp: str, rows: list) -> dict:
        n = len(rows)
        rep = rows[0]["gtdb_genome_representative"]
        return {"species": sp.split(";")[-1], "genomes": n, "bp": sum(r["genome_size"] for r in rows),
                "pairs_all": n * (n - 1) // 2, "pairs_farthest_first": (cap - 1) * n - cap * (cap - 1) // 2,
                "representative_complete": any(r["gtdb_accession"] == rep for r in rows)}

    stats = [sp_stats(sp, rows) for sp, rows in over.items()]
    rank1 = [rows[0] for rows in over.values()]
    rest_after_cap = [r for rows in over.values() for r in rows[1:cap]]
    out = {
        "gtdb_files_sha256": sha, "cap": cap,
        "species_over_cap": len(over),
        "genomes": sum(s["genomes"] for s in stats),
        "bp": sum(s["bp"] for s in stats),
        "contigs": sum(r["contig_count"] for rows in over.values() for r in rows),
        "archaea_genomes": sum(len(rows) for sp, rows in over.items() if sp.startswith("d__Archaea")),
        "rank1_held_by_phase1": {"genomes": len(rank1), "bp": sum(r["genome_size"] for r in rank1)},
        "quality_ranks_2_to_cap": {"genomes": len(rest_after_cap), "bp": sum(r["genome_size"] for r in rest_after_cap)},
        "pairs_all": sum(s["pairs_all"] for s in stats),
        "pairs_farthest_first": sum(s["pairs_farthest_first"] for s in stats),
        "representative_complete_species": sum(s["representative_complete"] for s in stats),
        "largest_species_genomes": max(s["genomes"] for s in stats),
        "largest_species_pairs_all": max(s["pairs_all"] for s in stats),
        "top10": sorted(stats, key=lambda s: (-s["genomes"], s["species"]))[:10],
    }
    out["to_download"] = {"genomes": out["genomes"] - len(rank1), "bp": out["bp"] - out["rank1_held_by_phase1"]["bp"]}
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
