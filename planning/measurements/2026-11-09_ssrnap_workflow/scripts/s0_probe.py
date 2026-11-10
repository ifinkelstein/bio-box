"""Probe the store: table sizes, and locate the well-known ssRNAP phages by header text."""
import time

from biobox.catalog import open_catalog

CAT = "/tmp/if2453/pilot_s1/root/catalog/release=rS1.pilot1.duckdb"
t0 = time.time()
con = open_catalog(CAT, local_dir="/tmp/if2453/ssrnap/cat", threads=2, memory_limit="8GB")
print("views:", sorted(r[0] for r in con.execute("SELECT table_name FROM information_schema.tables").fetchall()))
for t in ("protein_seq", "protein_occurrence", "contig_occurrence", "contig_tapes", "cluster_member", "contig_index"):
    print(t, con.execute(f"SELECT count(*) FROM {t}").fetchone()[0])
print(con.execute("SELECT level, count(*), count(DISTINCT cluster_uid) FROM cluster_member GROUP BY level").fetchall())
print(con.execute("SELECT * FROM contig_occurrence LIMIT 3").fetchall())
pats = {
    "T7": r"phage T7\b", "T3": r"phage T3\b", "SP6": r"phage SP6\b", "K11": r"phage K11\b",
    "phiKMV": r"phiKMV|phi[- ]?KMV", "N4": r"phage N4\b", "gh-1": r"gh-1", "phiYeO3-12": r"YeO3-12",
    "LKA1": r"LKA1", "Berlin": r"phage Berlin", "phiA1122": r"A1122", "K1F": r"K1F", "P60": r"P60\b",
    "Syn5": r"Syn5\b", "P-SSP7": r"P-SSP7",
}
for k, p in pats.items():
    rows = con.execute(
        "SELECT accession, original_header, contig_len_of(contig_uid) FROM contig_occurrence "
        "WHERE regexp_matches(original_header, ?) ORDER BY accession LIMIT 6", [p]).fetchall()
    print(f"== {k}: {len(rows)}")
    for r in rows:
        print("   ", r)
print(f"# {time.time() - t0:.1f} s")
