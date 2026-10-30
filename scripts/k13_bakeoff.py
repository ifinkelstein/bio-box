"""K-13: Parquet seq_chunks vs bgzip+faidx, scored by LY-06 (data_layer_and_tests.md §2A.1, §4.11).

    python scripts/k13_bakeoff.py GOLDEN_DIR OUT.json [--n-gene 500] [--n-nbr 200]

Every layout holds the whole golden set in one file, as a production shard would
(256 MB-1 GB holds many contigs), so per-file overheads do not dominate. For each
canonical query class (a 1 kb gene slice, a 20 kb neighbourhood slice, a whole contig)
it records storage units read (Parquet row groups from DuckDB's profiler; BGZF blocks
from the .gzi), compressed bytes read, an upper bound on bytes decoded, storage bytes/bp,
and wall time (recorded, not used to decide).
"""

from __future__ import annotations

import argparse
import json
import random
import statistics
import tempfile
import time
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import yaml

from biobox import fasta, ids
from biobox.catalog import profile_query
from biobox.l0seq.extract import BgzfFasta, write_bgzf_fasta
from biobox.l0seq.store import SEQ_CHUNK_STATS

REPO = Path(__file__).resolve().parents[1]
PARQUET_LAYOUTS = [  # (chunk bytes, chunks per row group)
    (65_536, 1_024),
    (65_536, 64),
    (65_536, 16),
    (16_384, 64),
    (16_384, 4),   # 64 KB of sequence per row group, the size of a BGZF block
    (4_096, 16),
]


def golden_contigs(golden_dir: Path) -> list[tuple[bytes, bytes]]:
    golden = yaml.safe_load((REPO / "golden.yaml").read_text())
    lock = yaml.safe_load((REPO / "golden.lock.yaml").read_text())
    seen, out = set(), []
    for e in golden["genomes"]:
        for r in fasta.read(golden_dir / lock["files"][e["id"]]["path"]):
            rd = ids.raw_digest(r.seq)
            if rd not in seen:
                seen.add(rd)
                out.append((rd, r.seq))
    return sorted(out)


def write_parquet(contigs, path: Path, chunk: int, rg_chunks: int) -> None:
    rows = {"raw_digest": [], "chunk_ix": [], "seq": [], "contig_len": []}
    for rd, seq in contigs:
        for ix in range(max(1, -(-len(seq) // chunk))):
            rows["raw_digest"].append(rd)
            rows["chunk_ix"].append(ix)
            rows["seq"].append(seq[ix * chunk : (ix + 1) * chunk].decode("ascii"))
            rows["contig_len"].append(len(seq))
    schema = pa.schema([("raw_digest", pa.binary(32)), ("chunk_ix", pa.uint32()), ("seq", pa.string()), ("contig_len", pa.int64())])
    pq.write_table(pa.table(rows, schema=schema), path, compression="zstd", row_group_size=rg_chunks,
                   write_statistics=SEQ_CHUNK_STATS)


def queries(contigs, rng: random.Random, n_gene: int, n_nbr: int) -> dict[str, list[tuple[int, int, int]]]:
    weights = [len(s) for _, s in contigs]
    out = {"gene_1kb": [], "nbr_20kb": [], "contig": []}
    for label, length, n in (("gene_1kb", 1_000, n_gene), ("nbr_20kb", 20_000, n_nbr)):
        while len(out[label]) < n:
            i = rng.choices(range(len(contigs)), weights)[0]
            L = len(contigs[i][1])
            if L <= length:
                continue
            a = rng.randrange(0, L - length)
            out[label].append((i, a, a + length))
    out["contig"] = [(i, 0, len(s)) for i, (_, s) in enumerate(contigs) if len(s) >= 1_000_000]
    return out


def _summ(xs):
    xs = sorted(xs)
    return {"mean": round(statistics.fmean(xs), 2), "p95": xs[int(0.95 * (len(xs) - 1))], "max": xs[-1]}


def run_parquet(path: Path, contigs, qs, chunk: int) -> dict:
    md = pq.ParquetFile(path).metadata
    seq_col = md.schema.to_arrow_schema().names.index("seq")
    rg_uncompressed = [md.row_group(i).column(seq_col).total_uncompressed_size for i in range(md.num_row_groups)]
    con = duckdb.connect(config={"threads": 1})
    con.execute(f"CREATE VIEW sc AS SELECT * FROM read_parquet('{path}')")
    sql = ("SELECT substr(string_agg(seq, '' ORDER BY chunk_ix), ? - min(chunk_ix) * ? + 1, ? - ?) FROM sc "
           "WHERE raw_digest = ? AND chunk_ix BETWEEN ? AND ?")
    res = {"bytes": path.stat().st_size, "row_groups": md.num_row_groups, "footer_bytes": md.serialized_size,
           "bytes_per_bp": round(path.stat().st_size / sum(len(s) for _, s in contigs), 4), "classes": {}}
    avg_rg = statistics.fmean(rg_uncompressed)
    for label, qq in qs.items():
        units, read, decoded, wall = [], [], [], []
        for i, a, b in qq:
            rd, seq = contigs[i]
            t0 = time.perf_counter()
            p = profile_query(con, sql, [a, chunk, b, a, rd, a // chunk, (b - 1) // chunk])
            wall.append(time.perf_counter() - t0)
            assert p["result"][0][0].encode() == seq[a:b]
            (scan,) = p["scans"]
            units.append(scan["row_groups_scanned"])
            read.append(p["bytes_read"])
            decoded.append(scan["row_groups_scanned"] * avg_rg)  # upper bound: whole seq column chunks
        res["classes"][label] = {"n": len(qq), "units": _summ(units), "bytes_read": _summ(read),
                                 "bytes_decoded_upper": _summ(decoded), "wall_ms": _summ([w * 1e3 for w in wall])}
    con.close()
    return res


def run_bgzf(path: Path, contigs, qs) -> dict:
    bg = BgzfFasta(path)
    size = sum(Path(str(path) + s).stat().st_size for s in ("", ".fai", ".gzi"))
    res = {"bytes": size, "bytes_per_bp": round(size / sum(len(s) for _, s in contigs), 4), "classes": {}}
    for label, qq in qs.items():
        units, read, decoded, wall = [], [], [], []
        for i, a, b in qq:
            rd, seq = contigs[i]
            b0, c0, u0 = bg.blocks_read, bg.compressed_bytes_read, bg.uncompressed_bytes_decoded
            t0 = time.perf_counter()
            got = bg.fetch(rd.hex(), a, b)
            wall.append(time.perf_counter() - t0)
            assert got == seq[a:b]
            units.append(bg.blocks_read - b0)
            read.append(bg.compressed_bytes_read - c0)
            decoded.append(bg.uncompressed_bytes_decoded - u0)
        res["classes"][label] = {"n": len(qq), "units": _summ(units), "bytes_read": _summ(read),
                                 "bytes_decoded_upper": _summ(decoded), "wall_ms": _summ([w * 1e3 for w in wall])}
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("golden_dir", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--n-gene", type=int, default=500)
    ap.add_argument("--n-nbr", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20261030)
    args = ap.parse_args()
    contigs = golden_contigs(args.golden_dir)
    qs = queries(contigs, random.Random(args.seed), args.n_gene, args.n_nbr)
    report = {"contigs": len(contigs), "bp": sum(len(s) for _, s in contigs), "duckdb": duckdb.__version__,
              "seed": args.seed, "layouts": {}}
    with tempfile.TemporaryDirectory() as d:
        for chunk, rg in PARQUET_LAYOUTS:
            p = Path(d) / f"p_{chunk}_{rg}.parquet"
            t0 = time.perf_counter()
            write_parquet(contigs, p, chunk, rg)
            r = run_parquet(p, contigs, qs, chunk)
            r["write_s"] = round(time.perf_counter() - t0, 1)
            report["layouts"][f"parquet_chunk{chunk // 1024}k_rg{rg}"] = r
            print(f"parquet {chunk}x{rg}", json.dumps(r["classes"]["gene_1kb"]["units"]), flush=True)
        bp = Path(d) / "golden.fa.gz"
        t0 = time.perf_counter()
        write_bgzf_fasta(((rd.hex(), s) for rd, s in contigs), bp)
        r = run_bgzf(bp, contigs, qs)
        r["write_s"] = round(time.perf_counter() - t0, 1)
        report["layouts"]["bgzf_64k_l6"] = r
    args.out.write_text(json.dumps(report, indent=1) + "\n")


if __name__ == "__main__":
    main()
