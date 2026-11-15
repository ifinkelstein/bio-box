"""Card B-48: what each gene-called stage of a root asked of gene calling, beside its run record (read-only).

    python b48_stage_inputs.py ROOT OUT_JSON [STAGE ...]

For every stage with `runs/stage=<S>/genecall.json` (or those named): its genomes from `genome_input`, each contig's
length from the `contig_index` of the stage that holds its default variant (any stage of the root), and the mode the
code calls each genome in (`genecall._finder`'s rule: virus gv-meta; plasmid, or under SINGLE_MODE_MIN_BP bases, meta;
else single). Per stage and mode, and by topology: genomes, contigs, bases, and the bases the caller reads (a circular
contig with its wrap flanks, `genecall._genes`: n + 2 min(n, WRAP_FLANK)); the same for the genomes that card B-48 (5)
would not call (every contig owned by another stage, DEC-26, or by another genome of the stage, DEC-27); the genes the
step wrote (rows of its `feature` files, from the footers); and the run record's cpu_min, wall_s, threads, workers,
started_at, image and peak. Also, per stage and mode, the called bases by genome length (bins), for a fit by length.
One DuckDB query per stage; nothing is written but OUT_JSON.
"""
import json
import sys
from pathlib import Path

import duckdb
import pyarrow.parquet as pq

from biobox import io
from biobox.l1 import genecall

root, out = Path(sys.argv[1]), Path(sys.argv[2])
names = sys.argv[3:] or sorted(p.name.removeprefix("stage=") for p in (root / "runs").iterdir()
                               if (p / "genecall.json").is_file())


def lits(files) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


all_index = {p.name.removeprefix("stage="): io.part_files(p) for p in sorted((root / "l0seq/contig_index").iterdir())}
con = duckdb.connect(config={"threads": 2, "memory_limit": "8GB"})
con.execute("CREATE TABLE ci AS " + " UNION ALL ".join(
    f"SELECT contig_uid, contig_len, '{s}' AS owner FROM read_parquet({lits(fs)}, hive_partitioning = false) "
    "WHERE is_default_variant" for s, fs in all_index.items() if fs))
res = {"root": str(root), "single_mode_min_bp": genecall.SINGLE_MODE_MIN_BP, "wrap_flank": genecall.WRAP_FLANK,
       "stages": {}}
BINS = [0, 10_000, 50_000, 100_000, 200_000, 500_000, 1_000_000, 3_000_000, 10**12]
for s in names:
    rec = json.loads((root / f"runs/stage={s}/genecall.json").read_text())
    gi = io.part_files(root / f"l0seq/genome_input/stage={s}")
    feats = [o["path"] for o in rec["outputs"] if o["path"].startswith("l1/feature/")]
    genes = sum(pq.read_metadata(root / p).num_rows for p in feats)
    argv = rec["argv"]
    q = f"""
    WITH gi AS (SELECT * FROM read_parquet({lits(gi)}, hive_partitioning = false)),
    least AS (SELECT contig_uid, min(genome_uid) AS g FROM gi GROUP BY contig_uid),
    c AS (SELECT gi.genome_uid, gi.contig_uid, gi.is_circular, gi.kind, ci.contig_len AS n, ci.owner,
                 (ci.owner = '{s}' AND least.g = gi.genome_uid) AS mine
          FROM gi JOIN ci USING (contig_uid) JOIN least USING (contig_uid)),
    g AS (SELECT genome_uid, any_value(kind) AS kind, count(*) AS contigs, sum(n) AS bases,
                 sum(CASE WHEN is_circular THEN n ELSE 0 END) AS circ_bases,
                 sum(n + CASE WHEN is_circular THEN 2 * least(n, {genecall.WRAP_FLANK}) ELSE 0 END) AS called,
                 bool_or(is_circular) AS any_circular, NOT bool_or(mine) AS all_owned_elsewhere
          FROM c GROUP BY genome_uid),
    m AS (SELECT *, CASE WHEN kind = 'virus' THEN 'gv-meta'
                         WHEN kind = 'plasmid' OR bases < {genecall.SINGLE_MODE_MIN_BP} THEN 'meta'
                         ELSE 'single' END AS mode FROM g)
    SELECT mode, any_circular, all_owned_elsewhere, count(*) AS genomes, sum(contigs) AS contigs, sum(bases) AS bases,
           sum(circ_bases) AS circ_bases, sum(called) AS called,
           {', '.join(f"sum(CASE WHEN bases >= {lo} AND bases < {hi} THEN called ELSE 0 END) AS called_bin{i}"
                      for i, (lo, hi) in enumerate(zip(BINS, BINS[1:])))},
           {', '.join(f"sum(CASE WHEN bases >= {lo} AND bases < {hi} THEN 1 ELSE 0 END) AS genomes_bin{i}"
                      for i, (lo, hi) in enumerate(zip(BINS, BINS[1:])))}
    FROM m GROUP BY ALL ORDER BY ALL"""
    rows = [dict(zip([d[0] for d in con.description], r)) for r in con.execute(q).fetchall()]
    listed = con.execute(f"SELECT count(DISTINCT contig_uid) FROM read_parquet({lits(gi)}, hive_partitioning = false)"
                         ).fetchone()[0]
    res["stages"][s] = {
        "record": {k: rec.get(k) for k in ("cpu_min", "wall_s", "threads", "started_at", "image_digest",
                                            "peak_rss_bytes", "over_cpu_budget")},
        "workers": int(argv[argv.index("--workers") + 1]) if "--workers" in argv else None,
        "seq_stages": argv[argv.index("--seq-stages") + 1] if "--seq-stages" in argv else None,
        "genes_written": genes, "contigs_listed": listed, "bins": BINS, "groups": rows}
    print(s, json.dumps(res["stages"][s]["record"]), genes, file=sys.stderr, flush=True)
out.write_text(json.dumps(res, indent=1, default=str))
print(json.dumps({s: {"cpu_min": v["record"]["cpu_min"], "genes": v["genes_written"]} for s, v in res["stages"].items()}))
