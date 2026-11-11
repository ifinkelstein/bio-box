"""Step 2 (direct path): every distinct protein of the release as FASTA, named by protein_uid hex."""
import time

from biobox.catalog import open_catalog

CAT = "/tmp/if2453/pilot_s1/root/catalog/release=rS1.pilot1.duckdb"
OUT = "/tmp/if2453/ssrnap/out/proteins.faa"
t0 = time.time()
con = open_catalog(CAT, local_dir="/tmp/if2453/ssrnap/cat", threads=2, memory_limit="8GB")
cur = con.execute("SELECT lower(hex(protein_uid)) AS n, replace(seq, '*', '') AS s FROM protein_seq ORDER BY protein_uid")
n = aa = 0
with open(OUT, "w") as fh:
    while True:
        rows = cur.fetchmany(200_000)
        if not rows:
            break
        for name, seq in rows:
            fh.write(f">{name}\n{seq}\n")
            n += 1
            aa += len(seq)
print(f"{n} proteins, {aa} residues, {time.time() - t0:.1f} s")
