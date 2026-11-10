"""Step 3: scale the micro-S0 fixture COPIES-fold and count what canonical queries read.

    python synth_scale.py REPO OUT COPIES [--shard-bases N] [--pair-width {0,1,2}] [--seed S]

Sequences: every fixture contig copied COPIES times with 3% random substitutions (distinct digests),
ingested by the repo's own `ingest` (digest-sorted bgzip shards). Gene tables: the fixture's gene
calls replicated per copy (substitutions keep coordinates). Tokens are synthetic: nr90 = 10 copies
per cluster, nr50 = 100 copies per cluster, family = one family per fixture gene, 20% without a
family, 3% in one hub family. Tapes, tokens and slot_pairs are written in the repo's layout (same
schemas, sort keys, row-group sizes, bucketing rule) and accepted by `catalog.build_manifest`.
Queries run through the generated catalog and are scored with `catalog.profile_query`.
"""

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import yaml
from blake3 import blake3

ap = argparse.ArgumentParser()
ap.add_argument("repo", type=Path)
ap.add_argument("out", type=Path)
ap.add_argument("copies", type=int)
ap.add_argument("--shard-bases", type=int, default=2**30)
ap.add_argument("--pair-width", type=int, default=None)
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()

from biobox import catalog, fasta, ids, layout  # noqa: E402
from biobox.l0seq import SeqStore, SourceObject, ingest  # noqa: E402
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, run as call_genes  # noqa: E402
from biobox.l3 import tapes as tapes_mod  # noqa: E402

# card B-11 removed biobox.l3.pairs and the slot_pairs table; the enumeration is the test oracle now. The
# slot_pairs section below is stale (it writes a table the registry no longer has); card B-18 moves this
# script to token_postings
sys.path.insert(0, str(a.repo / "tests"))
from oracles import pairs_oracle as pairs_mod  # noqa: E402
from biobox.l3.tokens import token_id  # noqa: E402
from biobox.registry.schemas import table  # noqa: E402

N, REL, V = a.copies, "rS0.syn", GENE_CALL_VERSION
rep: dict = {"copies": N, "shard_bases": a.shard_bases}
T0 = time.perf_counter()


def rss_mb() -> int:
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 if sys.platform == "linux" else 1024 * 1024))


def tick(name: str, t: float) -> None:
    rep.setdefault("build_s", {})[name] = round(time.perf_counter() - t, 1)
    rep.setdefault("peak_rss_mb_after", {})[name] = rss_mb()


# --- base: the fixture through the repo's builders ---------------------------------------------
fx = a.repo / "tests/fixtures/micro_s0"
src = yaml.safe_load((fx / "sources.yaml").read_text())
base_root = a.out / "base"
ingest([SourceObject(fx / o["file"], "micro", "1") for o in src["objects"]], base_root, "S0")
base_recs, genomes = [], []
for o in src["objects"]:
    r = fasta.read(fx / o["file"])[0]
    cu = ids.contig_uid(r.seq)
    base_recs.append((o["accession"], r.seq.upper(), cu, o["topology"] == "circular"))
    genomes.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"], circular=(cu,) if o["topology"] == "circular" else ()))
call_genes(base_root, "S0", genomes)
tapes_mod.run(base_root, "S0", {g.circular[0] for g in genomes if g.circular})
base_tapes = {t["contig_uid"]: t for t in pq.read_table(next((base_root / "l3/tapes").rglob("part-0.parquet"))).to_pylist()}
base_occ = pq.read_table(next((base_root / "l1/protein_occurrence").rglob("part-0.parquet"))).to_pylist()
B = len(base_recs)
n_slots = [base_tapes[cu]["n_slots"] for _, _, cu, _ in base_recs]
g_off = np.concatenate([[0], np.cumsum(n_slots)])  # global fixture-gene index of each contig's slot 0
G = int(g_off[-1])
rep["fixture"] = {"contigs": B, "bp": sum(len(s) for _, s, _, _ in base_recs), "slots": G}

# --- sequences ---------------------------------------------------------------------------------
t = time.perf_counter()
root = a.out / f"x{N}" / "root"
srcdir = a.out / f"x{N}" / "src"
srcdir.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(a.seed)
ACGT = np.frombuffer(b"ACGT", dtype=np.uint8)
cuid = np.empty((N, B), dtype="S16")
objs, fh = [], None
for c in range(N):
    if c % 50 == 0:
        if fh:
            fh.close()
        p = srcdir / f"part-{c // 50:04d}.fna"
        fh = open(p, "wb")
        objs.append(SourceObject(p, "synth", "1"))
    for b, (acc, seq, _, _) in enumerate(base_recs):
        arr = np.frombuffer(seq, dtype=np.uint8).copy()
        pos = rng.choice(arr.size, size=int(arr.size * 0.03), replace=False)
        arr[pos] = ACGT[rng.integers(0, 4, size=pos.size)]
        s = arr.tobytes()
        cuid[c, b] = ids.contig_uid(s)
        fh.write(b">syn%d_%s\n" % (c, acc.encode()))
        fh.write(b"\n".join(s[i : i + 80] for i in range(0, len(s), 80)) + b"\n")
fh.close()
tick("write_fasta", t)
t = time.perf_counter()
rep["ingest"] = ingest(objs, root, "S0", shard_bases=a.shard_bases)
tick("ingest (repo builder)", t)

# --- gene tables -------------------------------------------------------------------------------
t = time.perf_counter()
uids = [[bytes(cuid[c, b]).ljust(16, b"\0") for b in range(B)] for c in range(N)]  # numpy S16 strips trailing NULs
order = sorted((uids[c][b], c, b) for c in range(N) for b in range(B))
rank = np.empty((N, B), dtype=np.int32)
for r, (_, c, b) in enumerate(order):
    rank[c, b] = r
sorted_uids = pa.array([u for u, _, _ in order], pa.binary(16))


def puid(base: bytes, c: int) -> bytes:
    return blake3(base + c.to_bytes(4, "big")).digest(length=16)


occ_by_contig = {}
for r in base_occ:
    occ_by_contig.setdefault(r["contig_uid"], []).append(r)
occ_rows, tape_rows = [], []
for c in range(N):
    for b, (_, _, bcu, _) in enumerate(base_recs):
        cu = uids[c][b]
        gu = ids.genome_uid([cu])
        for r in occ_by_contig[bcu]:
            occ_rows.append({**r, "protein_uid": puid(r["protein_uid"], c), "contig_uid": cu, "genome_uid": gu})
        bt = base_tapes[bcu]
        tape_rows.append({**bt, "contig_uid": cu, "genome_uid": gu, "protein_uids": [puid(p, c) for p in bt["protein_uids"]]})
occ_rows.sort(key=lambda r: (r["protein_uid"], r["contig_uid"], r["slot_ordinal"]))
p = root / f"l1/protein_occurrence/gene_call_version={V}/stage=S0/part-0.parquet"
p.parent.mkdir(parents=True)
pq.write_table(pa.Table.from_pylist(occ_rows, schema=table("protein_occurrence").schema), p, compression="zstd", row_group_size=122_880)
del occ_rows
tape_rows.sort(key=lambda r: r["contig_uid"])
p = root / f"l3/tapes/gene_call_version={V}/stage=S0/part-0.parquet"
p.parent.mkdir(parents=True)
pq.write_table(pa.Table.from_pylist(tape_rows, schema=table("contig_tapes").schema), p, compression="zstd",
               row_group_size=tapes_mod.ROW_GROUP_CONTIGS, write_statistics=["contig_uid", "tax_order_bucket", "n_slots", "contig_len"])
tick("occurrence+tapes", t)

# --- tokens ------------------------------------------------------------------------------------
t = time.perf_counter()
g = np.arange(G)
fam_kind = np.where(g % 33 == 0, 2, np.where(g % 5 == 4, 0, 1))  # 2 hub, 0 none, 1 own family
fam_keys = sorted({"HUB"} | {f"F{i:06d}" for i in g[fam_kind == 1]})
fam_local = {k: i for i, k in enumerate(fam_keys)}
fam_tok = np.array([0 if k == 0 else token_id("family", 1, fam_local["HUB" if k == 2 else f"F{i:06d}"]) for i, k in zip(g, fam_kind)],
                   dtype=np.uint64)
TOK = {"family": np.broadcast_to(fam_tok, (N, G))}
dict_rows = [{"token_id": token_id("family", 1, i), "level": "family", "source": "synfam-1", "key": k,
              "n_slots": int(N * ((fam_kind == 2).sum() if k == "HUB" else 1))} for k, i in fam_local.items()]
for lvl, per in (("nr50", 100), ("nr90", 10)):
    groups = (N + per - 1) // per
    keys = sorted(f"{i:06d}:{q:05d}" for i in range(G) for q in range(groups))
    local = {k: j for j, k in enumerate(keys)}
    lut = np.array([[token_id(lvl, 0, local[f"{i:06d}:{q:05d}"]) for i in range(G)] for q in range(groups)], dtype=np.uint64)
    TOK[lvl] = lut[np.arange(N) // per]
    sizes = np.bincount(np.arange(N) // per, minlength=groups)
    dict_rows += [{"token_id": int(lut[q, i]), "level": lvl, "source": "mmseqs2", "key": f"{i:06d}:{q:05d}", "n_slots": int(sizes[q])}
                  for q in range(groups) for i in range(G)]
dict_rows.sort(key=lambda r: r["token_id"])
p = root / f"l3/token_dict/release={REL}/part-0.parquet"
p.parent.mkdir(parents=True)
pq.write_table(pa.Table.from_pylist(dict_rows, schema=table("token_dict").schema), p, compression="zstd", row_group_size=122_880,
               write_statistics=["token_id", "n_slots"])
tok_rows = []
for _, c, b in order:
    sl = slice(int(g_off[b]), int(g_off[b + 1]))
    tok_rows.append({"contig_uid": uids[c][b], **{f"{lvl}_tokens": [int(x) or None for x in TOK[lvl][c, sl]] for lvl in TOK}})
p = root / f"l3/tape_tokens/gene_call_version={V}/release={REL}/part-0.parquet"
p.parent.mkdir(parents=True)
pq.write_table(pa.Table.from_pylist(tok_rows, schema=table("tape_tokens").schema), p, compression="zstd", row_group_size=8_192,
               write_statistics=["contig_uid"])
del tok_rows
tick("tokens", t)

# --- slot_pairs (vectorized replica of pairs.pairs_of_contig) ------------------------------------
t = time.perf_counter()
cols = {k: [] for k in ("a", "p", "rk", "i", "j", "d", "topo", "delta", "flags")}
for b, (_, _, bcu, _) in enumerate(base_recs):
    bt = base_tapes[bcu]
    n = bt["n_slots"]
    base_pairs = pairs_mod.pairs_of_contig(bt, {"family": [None] * n, "nr50": [None] * n, "nr90": list(range(1, n + 1))})
    I = np.array([r["anchor_ordinal"] for r in base_pairs], dtype=np.uint32)
    J = np.array([r["partner_ordinal"] for r in base_pairs], dtype=np.uint32)
    geo = {k: np.array([r[src_k] for r in base_pairs], dtype=dt) for k, src_k, dt in
           (("d", "orf_offset", np.int16), ("topo", "junction_topology", np.uint8), ("delta", "delta_bp", np.int32), ("flags", "flags", np.uint8))}
    for lvl in TOK:
        tk = TOK[lvl][:, int(g_off[b]) : int(g_off[b + 1])]
        A, P = tk[:, I], tk[:, J]
        keep = (A != 0) & (P != 0)
        cols["a"].append(A[keep])
        cols["p"].append(P[keep])
        cols["rk"].append(np.broadcast_to(rank[:, b : b + 1], A.shape)[keep])
        cols["i"].append(np.broadcast_to(I, A.shape)[keep])
        cols["j"].append(np.broadcast_to(J, A.shape)[keep])
        for k in geo:
            cols[k].append(np.broadcast_to(geo[k], A.shape)[keep])
cols = {k: np.concatenate(v) for k, v in cols.items()}
n_pairs = int(cols["a"].size)
perm = np.lexsort((cols["j"], cols["i"], cols["rk"], cols["p"], cols["a"]))
cols = {k: v[perm] for k, v in cols.items()}
del perm
width = layout.prefix_width(pairs_mod.BYTES_PER_ROW * n_pairs) if a.pair_width is None else a.pair_width
low = (cols["a"] & np.uint64(0xFF)).astype(np.uint8)
spec = table("slot_pairs")
max_rows = layout.MAX_FILE_BYTES // pairs_mod.BYTES_PER_ROW
n_files = 0
for bk in sorted({layout.bucket(bytes([x]), width) for x in range(256)}):
    sel = np.flatnonzero(np.array([layout.bucket(bytes([x]), width) == bk for x in range(256)])[low])
    if not sel.size:
        continue
    for part, lo in enumerate(range(0, sel.size, max_rows)):
        s = sel[lo : lo + max_rows]
        tbl = pa.table({
            "anchor_token": cols["a"][s], "partner_token": cols["p"][s], "contig_uid": sorted_uids.take(pa.array(cols["rk"][s])),
            "anchor_ordinal": cols["i"][s], "partner_ordinal": cols["j"][s], "orf_offset": cols["d"][s],
            "junction_topology": cols["topo"][s], "delta_bp": cols["delta"][s], "flags": cols["flags"][s],
            "tier": pa.array(np.ones(s.size, dtype=np.uint8)), "tax_order": pa.nulls(s.size, pa.uint32()),
        }).cast(spec.schema)
        d = root / f"l3/slot_pairs/release={REL}/anchor_bucket={bk}"
        d.mkdir(parents=True, exist_ok=True)
        pq.write_table(tbl, d / f"part-{part}.parquet", compression="zstd", row_group_size=122_880,
                       write_statistics=["anchor_token", "partner_token", "contig_uid", "anchor_ordinal", "partner_ordinal"])
        n_files += 1
rep["slot_pairs"] = {"rows": n_pairs, "rows_per_slot": round(n_pairs / (N * G), 2), "bucket_width": width, "files": n_files}
del cols, low
tick("slot_pairs", t)

# --- manifest and catalog ----------------------------------------------------------------------
t = time.perf_counter()
man = catalog.build_manifest(root, REL)
tick("build_manifest (validates layout)", t)
per = {}
for f in man["files"]:
    q = per.setdefault(f["table"], {"files": 0, "rows": 0, "bytes": 0, "row_groups": 0})
    for k in ("rows", "bytes", "row_groups"):
        q[k] += f[k]
    q["files"] += 1
shards = [x for x in man["artifacts"] if x["path"].endswith(".fa.gz")]
per["seq_shards"] = {"files": len(shards), "bytes": sum(x["bytes"] for x in shards)}
for k, v in per.items():
    if v.get("rows"):
        v["bytes_per_row"] = round(v["bytes"] / v["rows"], 2)
rep["tables"] = per
cat = catalog.generate_catalog(man, a.out / f"x{N}" / "catalog.duckdb")
r0 = rss_mb()
t = time.perf_counter()
con = catalog.open_catalog(cat)
rep["open_catalog"] = {"seconds": round(time.perf_counter() - t, 2), "contigs_in_locator": N * B}

# --- queries -----------------------------------------------------------------------------------


def run(name: str, sql: str, params=None) -> None:
    t = time.perf_counter()
    out = catalog.profile_query(con, sql, params)
    res = out["result"]
    scans = []
    for s in out["scans"]:
        if not s["files"]:
            continue
        rel = Path(s["files"][0]).relative_to(root).parts
        tname = rel[1] if rel[0] in ("l0seq", "l1", "l2", "l3") else rel[0]
        scans.append(f"{tname.removesuffix('.parquet')}: files {s['files_read']}, row groups {s['row_groups_scanned']}/"
                     f"{s['row_groups_total']}, rows {s['rows_scanned']}")
    rep.setdefault("queries", {})[name] = {
        "seconds": round(time.perf_counter() - t, 3), "result": str(res[0]) if len(res) == 1 else f"{len(res)} rows",
        "scans": scans, "parquet_bytes_read": out["bytes_read"], "bgzf": out["bgzf"]}


big = max(range(B), key=lambda b: n_slots[b])
gi = next(i for i in range(int(g_off[big]) + 50, int(g_off[big + 1])) if fam_kind[i] == 1 and fam_kind[i + 1] == 1)
fa, fb = int(fam_tok[gi]), int(fam_tok[gi + 1])
hub = token_id("family", 1, fam_local["HUB"])
na, nb = int(TOK["nr90"][N // 2, gi]), int(TOK["nr90"][N // 2, gi + 1])
run("Q1 family pair (each token on COPIES loci)", "SELECT q1_count(?, ?, 3, NULL)", [fa, fb])
run("Q1 nr90 pair (10 loci)", "SELECT q1_count(?, ?, 3, NULL)", [na, nb])
run("Q1 hub anchor x family partner", "SELECT q1_count(?, ?, 10, NULL)", [hub, fa])
run("Q4 family pair instances", "SELECT count(*), sum(delta_bp) FROM pairs_of(?, ?, 10)", [fa, fb])
run("Q4 nr90 pair instances", "SELECT count(*), sum(delta_bp) FROM pairs_of(?, ?, 10)", [na, nb])
run("Q3 partners of hub", "SELECT count(*), sum(n_pairs) FROM partners_of(?, 10)", [hub])
run("Q3 partners of one family", "SELECT count(*), sum(n_pairs) FROM partners_of(?, 10)", [fa])

qrng = np.random.default_rng(7)
anchors = []
for _ in range(1000):
    c, b = int(qrng.integers(N)), int(qrng.integers(B))
    bt = base_tapes[base_recs[b][2]]
    i = int(qrng.integers(bt["n_slots"]))
    lo, hi = max(i - 10, 0), min(i + 10, bt["n_slots"] - 1)
    anchors.append((uids[c][b], bt["start_bps"][lo], min(bt["end_bps"][hi], bt["contig_len"])))
con.register("anchors", pa.table({"contig_uid": pa.array([x[0] for x in anchors], pa.binary()),
                                  "s": pa.array([x[1] for x in anchors], pa.int64()), "e": pa.array([x[2] for x in anchors], pa.int64())}))
run("1,000 neighborhoods: sequence (+-10 genes)", "SELECT sum(length(seq_slice(contig_uid, s, e))) FROM anchors")
nb_bp = sum(x[2] - x[1] for x in anchors)
rep["queries"]["1,000 neighborhoods: sequence (+-10 genes)"]["bp_returned"] = nb_bp
run("1,000 neighborhoods: tapes and tokens", "SELECT sum(n_slots), sum(len(nr90_tokens)) FROM tapes_of(?)", [sorted({x[0] for x in anchors})])
run("1 neighborhood: tape and tokens", "SELECT sum(n_slots), sum(len(nr90_tokens)) FROM tapes_of(?)", [[anchors[0][0]]])
run("1 neighborhood: slots of a protein", "SELECT count(*), sum(start) FROM slots_of(?)", [puid(base_tapes[base_recs[big][2]]["protein_uids"][60], N // 2)])

store = SeqStore(root)
for k in (1, 50, 500):
    k = min(k, N * B)
    pick = [uids[int(qrng.integers(N))][int(qrng.integers(B))] for _ in range(k)]
    con.close()
    con = catalog.open_catalog(cat)  # fresh block cache
    con.register("g", pa.table({"contig_uid": pa.array(pick, pa.binary())}))
    name = f"export proxy: {k} whole contigs (a genome with {k} contigs)"
    run(name, "SELECT sum(length(seq_slice(contig_uid, 0, contig_len_of(contig_uid)))) FROM g")
    q = rep["queries"][name]
    q["shards_touched"] = len({store.shard_of(store.default_digest(u)) for u in pick})
    q["shards_total"] = len(shards)
    q["bp_returned"] = sum(store.length_of(store.default_digest(u)) for u in pick)
rep["peak_rss_mb_queries"] = rss_mb()
rep["rss_mb_before_open"] = r0
rep["total_s"] = round(time.perf_counter() - T0, 1)
print(json.dumps(rep, indent=1, default=str))
