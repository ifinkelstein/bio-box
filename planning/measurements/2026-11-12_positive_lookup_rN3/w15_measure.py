"""W-15 measurement: cost of resolving an accession (find --accession / search hmm --positive) on a real release.

Read-only on the root. Run inside the image with the W-15 tree first on PYTHONPATH.
usage: python w15_measure.py ROOT RELEASE OUTDIR
"""
import json
import subprocess
import sys
import time
from pathlib import Path

import pyarrow.parquet as pq

from biobox import catalog, cli
from biobox.search import export, nbr

root, release, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
THREADS = 2
files = export.release_files(root, release, export.release_manifest(root, release))
gcv = cli.GENE_CALL_VERSION
tv = nbr.resolve_tape_version(gcv, None)
rec: dict = {"release": release, "root": str(root), "threads": THREADS, "tape_version": tv}

# what the label relation reads on this release
for t in ("contig_label", "contig_occurrence"):
    ents = files.entries(t)
    rec[t] = [{"path": str(p), "parts": {k: str(v) for k, v in parts.items()},
               "rows": pq.read_metadata(p).num_rows, "row_groups": pq.read_metadata(p).num_row_groups}
              for p, parts in ents]

# pick accessions: per stage, the contig with the most genes (worst case for the tape read) and a one-record one
picks = []
with export.open_store(root, release=release, gene_call_version=gcv, manifest=files, threads=THREADS,
                       relations=("lab", "cidx", "dim")) as st:
    rows = st.con.execute("SELECT accession, contig_uid FROM lab WHERE accession IS NOT NULL "
                          "ORDER BY contig_uid LIMIT 3").fetchall()
    picks = [r[0] for r in rows]
    # one long prokaryote contig: the longest default-variant contig
    longest = st.con.execute("SELECT l.accession FROM lab l JOIN cidx c USING (contig_uid) WHERE c.is_default_variant "
                             "ORDER BY c.contig_len DESC, l.accession LIMIT 1").fetchone()[0]
    picks.append(longest)
    rec["picks"] = picks
    # profile the lookup itself: the SQL _contig_labels runs (warm-up first, then measured)
    prof = []
    for acc in picks:
        sql = f"SELECT DISTINCT contig_uid FROM lab WHERE accession = {export._q(acc)}"
        t0 = time.perf_counter()
        p = catalog.profile_query(st.con, sql, profile_dir=out)
        dt = time.perf_counter() - t0
        prof.append({"accession": acc, "seconds": round(dt, 4), "n_contigs": len(p["result"]),
                     "scans": p["scans"], "bytes_read": p.get("bytes_read")})
    # header substring: the slow path
    sql = "SELECT DISTINCT contig_uid FROM lab WHERE contains(original_header, 'Escherichia coli')"
    t0 = time.perf_counter()
    p = catalog.profile_query(st.con, sql, profile_dir=out)
    prof.append({"header_substring": "Escherichia coli", "seconds": round(time.perf_counter() - t0, 4),
                 "n_contigs": len(p["result"]), "scans": p["scans"], "bytes_read": p.get("bytes_read")})
    rec["profile"] = prof

# timings of the whole resolution as find and search hmm --positive call it (in-process)
tim = []
for acc in picks:
    t0 = time.perf_counter()
    labels = cli._contig_labels(root, files, gcv, THREADS, where=f"accession = {export._q(acc)}")
    t1 = time.perf_counter()
    genes, _ = cli.find_genes(root, files, accession=acc, tape_version=tv, gene_call_version=gcv, threads=THREADS)
    t2 = time.perf_counter()
    g = genes[len(genes) // 2] if genes else None
    spec = f"{acc}:{g['start']}-{g['end']}" if g else acc
    uid, lines = cli._positive_uid(root, files, spec, tv, gcv, THREADS)
    t3 = time.perf_counter()
    tim.append({"accession": acc, "spec": spec, "labels_s": round(t1 - t0, 4), "find_genes_s": round(t2 - t1, 4),
                "positive_uid_s": round(t3 - t2, 4), "n_genes": len(genes), "positive": uid, "refusal": lines[:2]})
rec["timings"] = tim

# end to end: the bio-box find command as a user runs it, a fresh process each time
e2e = []
for t in tim:
    cmd = [sys.executable, "-m", "biobox.cli", "find", "--root", str(root), "--release", release,
           "--accession", t["spec"], "--threads", str(THREADS)]
    t0 = time.perf_counter()
    r = subprocess.run(cmd, capture_output=True, text=True)
    e2e.append({"cmd": " ".join(cmd[1:]), "seconds": round(time.perf_counter() - t0, 3), "rc": r.returncode,
                "stdout_lines": len(r.stdout.splitlines()), "stderr": r.stderr.strip().splitlines()[:3]})
# python start-up and import alone, for comparison
t0 = time.perf_counter()
subprocess.run([sys.executable, "-c", "import biobox.cli"], check=True)
e2e.append({"cmd": "python -c 'import biobox.cli'", "seconds": round(time.perf_counter() - t0, 3)})
rec["end_to_end"] = e2e
(out / f"w15_{release}.json").write_text(json.dumps(rec, indent=1, default=str))
print(json.dumps(rec, indent=1, default=str))
