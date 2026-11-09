"""Lower bound on slot_pairs rows per level from GTDB R226 metadata.

    python gtdb_pairs.py GTDB_DIR

Rows per slot on a linear contig of n genes within +-10: 20 - 110/n for n >= 11, else n - 1. Each
genome's genes are spread evenly over its contigs (a lower bound: real assemblies put most genes on
the long contigs). Upper bound is 20 per slot.
"""
import json
import sys

import duckdb

d = sys.argv[1]
con = duckdb.connect()
con.execute(f"""
CREATE VIEW m AS SELECT gtdb_taxonomy, gtdb_representative, ncbi_assembly_level,
       CAST(genome_size AS BIGINT) AS bp, CAST(contig_count AS DOUBLE) AS c, CAST(protein_count AS DOUBLE) AS p,
       TRY_CAST(checkm2_completeness AS DOUBLE) AS comp, TRY_CAST(checkm2_contamination AS DOUBLE) AS cont, accession
FROM read_csv(['{d}/bac120_metadata_r226.tsv.gz', '{d}/ar53_metadata_r226.tsv.gz'], delim='\t', header=true,
              all_varchar=true, union_by_name=true)""")
expr = "sum(p * CASE WHEN p / c >= 11 THEN 20 - 110 / (p / c) ELSE greatest(p / c - 1, 0) END)"
q = {
    "S2C": f"""SELECT sum(p), {expr} FROM (SELECT *, row_number() OVER (PARTITION BY gtdb_taxonomy ORDER BY comp - 5 * cont DESC NULLS LAST,
               accession) rk FROM m WHERE ncbi_assembly_level = 'Complete Genome') WHERE rk <= 20""",
    "S5": f"SELECT sum(p), {expr} FROM m WHERE gtdb_representative = 't'",
    "S6": f"SELECT sum(p), {expr} FROM m",
}
out = {}
for k, sql in q.items():
    genes, rows = con.execute(sql).fetchall()[0]
    out[k] = {"genes": genes, "pair_rows_per_level_lower": rows, "rows_per_gene_lower": round(rows / genes, 2), "rows_per_gene_upper": 20}
print(json.dumps(out, indent=1))
