"""Project every table to GTDB R226 scale from measured bytes per unit (card B-18; replaces gtdb_pairs.py).

    python gtdb_postings.py GTDB_DIR UNITS_JSON > projection.json

GTDB_DIR has `bac120_metadata_r226.tsv.gz` and `ar53_metadata_r226.tsv.gz` (as for gtdb_counts.py).
Unit counts per stage come from the metadata with the same selections as gtdb_counts.py: S2C (complete
genomes, at most 20 per species), S3 (species representatives under 2 Mbp), S5 (all species
representatives), S6 (all genomes). Units: `bp` (genome_size), `contig` (contig_count), `gene`
(protein_count; a tape slot is a gene while tapes hold only the gene-call producer), `genome`.

UNITS_JSON maps each table to {"unit": "bp" | "contig" | "gene" | "genome", "bytes_per_unit": B, "source": "..."};
an entry may give "rows_per_unit" and "bytes_per_row" instead of "bytes_per_unit" (token_postings:
rows per gene times bytes per row). The bytes per unit are measured elsewhere; this script only
multiplies: projected bytes = units x bytes per unit (*derived*). `token_postings` replaces the pair
projection of gtdb_pairs.py (AGENTS §3.5: pairs are not materialized), which is why that script was
deleted. No deduplication is assumed for `protein_seq` (one distinct protein per gene), an upper bound.
"""

import json
import sys

if len(sys.argv) < 3 or sys.argv[1] in ("-h", "--help"):
    print(__doc__)
    sys.exit(0 if sys.argv[1:2] in (["-h"], ["--help"]) else 2)

import duckdb  # noqa: E402

d, units = sys.argv[1], json.load(open(sys.argv[2]))
con = duckdb.connect()
con.execute(f"""
CREATE VIEW m AS SELECT accession, gtdb_taxonomy, gtdb_representative, ncbi_assembly_level,
       CAST(genome_size AS BIGINT) AS bp, CAST(contig_count AS BIGINT) AS contig,
       CAST(protein_count AS BIGINT) AS gene,
       TRY_CAST(checkm2_completeness AS DOUBLE) AS comp, TRY_CAST(checkm2_contamination AS DOUBLE) AS cont
FROM read_csv(['{d}/bac120_metadata_r226.tsv.gz', '{d}/ar53_metadata_r226.tsv.gz'], delim='\t', header=true,
              all_varchar=true, union_by_name=true)""")
agg = "count(*) AS genomes, sum(bp) AS bp, sum(contig) AS contig, sum(gene) AS gene"
stages = {
    "S2C": f"""SELECT {agg} FROM (SELECT *, row_number() OVER (PARTITION BY gtdb_taxonomy ORDER BY comp - 5 * cont DESC
               NULLS LAST, accession) AS rk FROM m WHERE ncbi_assembly_level = 'Complete Genome') WHERE rk <= 20""",
    "S3": f"SELECT {agg} FROM m WHERE gtdb_representative = 't' AND bp < 2000000",
    "S5": f"SELECT {agg} FROM m WHERE gtdb_representative = 't'",
    "S6": f"SELECT {agg} FROM m",
}
counts = {}
for k, sql in stages.items():
    row = con.execute(sql).fetchall()[0]
    counts[k] = dict(zip([c[0] for c in con.description], [int(x) for x in row]))
    counts[k]["genome"] = counts[k]["genomes"]
out = {"counts": counts, "tables": {}, "total_gb": {}}
for name, u in sorted(units.items()):
    if name.startswith("_"):
        continue
    bpu = u["bytes_per_unit"] if "bytes_per_unit" in u else u["rows_per_unit"] * u["bytes_per_row"]
    out["tables"][name] = {"unit": u["unit"], "bytes_per_unit": round(bpu, 4), "source": u.get("source", ""),
                           "gb": {k: round(c[u["unit"]] * bpu / 1e9, 2) for k, c in counts.items()}}
for k in counts:
    out["total_gb"][k] = round(sum(t["gb"][k] for t in out["tables"].values()), 1)
print(json.dumps(out, indent=1))
