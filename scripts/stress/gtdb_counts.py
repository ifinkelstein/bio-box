"""Stage unit counts from GTDB R226 metadata (bac120 + ar53).

    python gtdb_counts.py GTDB_DIR
"""

import json
import sys

import duckdb

if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
    print(__doc__)
    sys.exit(0 if sys.argv[1:2] in (["-h"], ["--help"]) else 2)

d = sys.argv[1]
con = duckdb.connect()
con.execute(f"""
CREATE VIEW m AS SELECT accession, gtdb_taxonomy, gtdb_representative, ncbi_assembly_level,
       CAST(genome_size AS BIGINT) AS genome_size, CAST(contig_count AS BIGINT) AS contig_count,
       CAST(protein_count AS BIGINT) AS protein_count,
       TRY_CAST(checkm2_completeness AS DOUBLE) AS comp, TRY_CAST(checkm2_contamination AS DOUBLE) AS cont
FROM read_csv(['{d}/bac120_metadata_r226.tsv.gz', '{d}/ar53_metadata_r226.tsv.gz'], delim='\t', header=true,
              all_varchar=true, union_by_name=true)
""")
agg = ("count(*) AS genomes, sum(genome_size) AS bp, sum(contig_count) AS contigs, sum(protein_count) AS proteins, "
       "median(contig_count) AS contigs_median, quantile_cont(contig_count, 0.99) AS contigs_p99, max(contig_count) AS contigs_max")
q = {
    "S6_all_genomes": f"SELECT {agg} FROM m",
    "S5_species_reps": f"SELECT {agg} FROM m WHERE gtdb_representative = 't'",
    "S3_reps_lt_2Mbp": f"SELECT {agg} FROM m WHERE gtdb_representative = 't' AND genome_size < 2000000",
    "complete_genomes_all": f"SELECT {agg} FROM m WHERE ncbi_assembly_level = 'Complete Genome'",
    "S2C_complete_le20_per_species": f"""SELECT {agg} FROM (SELECT *, row_number() OVER (PARTITION BY gtdb_taxonomy
         ORDER BY comp - 5 * cont DESC NULLS LAST, accession) AS rk FROM m WHERE ncbi_assembly_level = 'Complete Genome') WHERE rk <= 20""",
}
out = {}
for k, sql in q.items():
    row = con.execute(sql).fetchall()[0]
    out[k] = dict(zip([c[0] for c in con.description], [float(x) if x is not None else None for x in row]))
out["levels"] = dict(con.execute("SELECT ncbi_assembly_level, count(*) FROM m GROUP BY 1 ORDER BY 2 DESC").fetchall())
print(json.dumps(out, indent=1))
