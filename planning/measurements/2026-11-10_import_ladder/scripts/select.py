"""Ladder genome selection from GTDB R226 metadata (bac120 + ar53).

    uv run --with duckdb python select.py GTDB_DIR OUT_DIR [--seed 20261109] [--b 250] [--c 2250]

Complete genomes (ncbi_assembly_level = 'Complete Genome'), one per species (gtdb_taxonomy):
best checkm2_completeness - 5 * checkm2_contamination, ties by accession. Species are shuffled
with a seeded RNG over the sorted species list; the first B are stage B (P1), the next C stage C (P2).
Writes OUT_DIR/stage_P1.tsv, stage_P2.tsv (accession without RS_/GB_, genome_size, contig_count,
protein_count, translation_table, gtdb_taxonomy, gtdb_accession) and selection.json (totals).
"""

import argparse
import json
import random
from pathlib import Path

import duckdb

ap = argparse.ArgumentParser()
ap.add_argument("gtdb", type=Path)
ap.add_argument("out", type=Path)
ap.add_argument("--seed", type=int, default=20261109)
ap.add_argument("--b", type=int, default=250)
ap.add_argument("--c", type=int, default=2250)
a = ap.parse_args()
a.out.mkdir(parents=True, exist_ok=True)

con = duckdb.connect()
files = [str(a.gtdb / "bac120_metadata_r226.tsv.gz"), str(a.gtdb / "ar53_metadata_r226.tsv.gz")]
con.execute(f"""
CREATE TABLE best AS SELECT * FROM (
  SELECT accession AS gtdb_accession, regexp_replace(accession, '^(RS_|GB_)', '') AS accession,
         gtdb_taxonomy, CAST(genome_size AS BIGINT) AS genome_size, CAST(contig_count AS BIGINT) AS contig_count,
         TRY_CAST(protein_count AS BIGINT) AS protein_count,
         TRY_CAST(ncbi_translation_table AS INTEGER) AS translation_table,
         TRY_CAST(checkm2_completeness AS DOUBLE) AS comp, TRY_CAST(checkm2_contamination AS DOUBLE) AS cont,
         row_number() OVER (PARTITION BY gtdb_taxonomy
                            ORDER BY TRY_CAST(checkm2_completeness AS DOUBLE) - 5 * TRY_CAST(checkm2_contamination AS DOUBLE)
                            DESC NULLS LAST, accession) AS rk
  FROM read_csv({files}, delim='\t', header=true, all_varchar=true, union_by_name=true)
  WHERE ncbi_assembly_level = 'Complete Genome') WHERE rk = 1
""")
species = [r[0] for r in con.execute("SELECT gtdb_taxonomy FROM best ORDER BY gtdb_taxonomy").fetchall()]
rng = random.Random(a.seed)
order = species[:]
rng.shuffle(order)
pick = {"P1": order[: a.b], "P2": order[a.b : a.b + a.c]}
assert not set(pick["P1"]) & set(pick["P2"])
summary = {"seed": a.seed, "species_with_complete_genome": len(species), "stages": {}}
cols = "accession, genome_size, contig_count, protein_count, translation_table, gtdb_taxonomy, gtdb_accession"
for stage, sp in pick.items():
    con.register("sp", __import__("pyarrow").table({"gtdb_taxonomy": sp}))
    rows = con.execute(f"SELECT {cols} FROM best JOIN sp USING (gtdb_taxonomy) ORDER BY accession").fetchall()
    con.unregister("sp")
    path = a.out / f"stage_{stage}.tsv"
    with open(path, "w") as fh:
        fh.write("\t".join(c.strip() for c in cols.split(",")) + "\n")
        for r in rows:
            fh.write("\t".join("" if v is None else str(v) for v in r) + "\n")
    tt = {}
    for r in rows:
        tt[str(r[4])] = tt.get(str(r[4]), 0) + 1
    summary["stages"][stage] = {
        "genomes": len(rows), "bp": sum(r[1] for r in rows), "contigs": sum(r[2] for r in rows),
        "proteins": sum(r[3] or 0 for r in rows), "proteins_missing": sum(r[3] is None for r in rows),
        "archaea": sum(r[5].startswith("d__Archaea") for r in rows), "translation_tables": tt,
        "max_genome_bp": max(r[1] for r in rows), "min_genome_bp": min(r[1] for r in rows), "file": str(path)}
(a.out / "selection.json").write_text(json.dumps(summary, indent=1))
print(json.dumps(summary, indent=1))
