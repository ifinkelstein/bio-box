"""Memory and open-time of two in-memory structures, measured so they can be scaled.

    python builder_limits.py ROOT SCRATCH

1. Bytes of Python heap per slot_pairs row as `pairs.pairs_of_contig` builds them (E. coli, one level).
2. `SeqStore(root)` open time and heap for a synthetic store index of 1e6 contigs.
"""

import json
import os
import sys
import time
import tracemalloc
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from biobox.l0seq import SeqStore
from biobox.l3.pairs import pairs_of_contig
from biobox.registry.schemas import table

root, scratch = Path(sys.argv[1]), Path(sys.argv[2])
out = {}
tapes = pq.read_table(next((root / "l3/tapes").rglob("part-0.parquet"))).to_pylist()
t = max(tapes, key=lambda r: r["n_slots"])
n = t["n_slots"]
tracemalloc.start()
rows = pairs_of_contig(t, {"family": [None] * n, "nr50": [None] * n, "nr90": list(range(1, n + 1))})
cur, _ = tracemalloc.get_traced_memory()
tracemalloc.stop()
out["pairs_builder"] = {"contig_slots": n, "rows": len(rows), "heap_bytes_per_row": round(cur / len(rows))}
del rows

M = 1_000_000
rng = np.random.default_rng(0)
rd = rng.bytes(32 * M)
cu = rng.bytes(16 * M)
rds = sorted(rd[i * 32 : (i + 1) * 32] for i in range(M))
cus = sorted(cu[i * 16 : (i + 1) * 16] for i in range(M))
r2 = scratch / "locator_1e6"
(r2 / "l0seq/contig_index/stage=S0").mkdir(parents=True, exist_ok=True)
(r2 / "l0seq/seq_locator/stage=S0").mkdir(parents=True, exist_ok=True)
LOC = r2 / "l0seq/seq_locator/stage=S0/part-0.parquet"
pq.write_table(pa.table({"raw_digest": pa.array(rds, pa.binary(32)),
                         "shard": pa.array([f"l0seq/seq_shards/stage=S0/shard-{i % 400:04d}.fa.gz" for i in range(M)]),
                         "contig_len": pa.array(np.full(M, 20000, dtype=np.int64))}).cast(table("seq_locator").schema),
               LOC, compression="zstd", row_group_size=122_880)
pq.write_table(pa.table({"contig_uid": pa.array(cus, pa.binary(16)), "raw_digest": pa.array(rds, pa.binary(32)),
                         "contig_len": pa.array(np.full(M, 20000, dtype=np.int64)),
                         "is_default_variant": pa.array(np.ones(M, dtype=bool))}).cast(table("contig_index").schema),
               r2 / "l0seq/contig_index/stage=S0/part-0.parquet", compression="zstd", row_group_size=122_880)
del rds, cus, rd, cu
tracemalloc.start()
t0 = time.perf_counter()
s = SeqStore(r2)
dt = time.perf_counter() - t0
cur, peak = tracemalloc.get_traced_memory()
out["seqstore_open_1e6_contigs"] = {
    "seconds": round(dt, 2), "heap_mb_held": round(cur / 2**20), "heap_mb_peak": round(peak / 2**20),
    "locator_file_mb": round(os.path.getsize(LOC) / 2**20, 1),
    "contig_index_file_mb": round(os.path.getsize(r2 / "l0seq/contig_index/stage=S0/part-0.parquet") / 2**20, 1)}
print(json.dumps(out, indent=1))
