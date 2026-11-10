"""Scale the micro-S0 fixture COPIES-fold and count what key lookups read (measurements file section 5).

    python synth_scale.py REPO OUT COPIES [--shard-bases N] [--postings-width {0,1,2}] [--seed S] > synth.json

Sequences: every fixture contig copied COPIES times with 3% random substitutions (distinct digests),
ingested by the repo's own `ingest`. Gene features: the fixture's gene calls replicated per copy
(substitutions keep coordinates; each copy's genes get their own `feature_uid` and a synthetic
`protein_uid`), written as `feature`, `feature_protein` and `protein_occurrence` through
`biobox.io.write_parts`. Clusters and families are synthetic: NR90 = one cluster per fixture protein per
10 copies, NR50 per 100 copies; families from an `l1_domains` library `synfam-1`, one family per fixture
protein, 20% of proteins without one and 3% in one hub family. Everything downstream is the repo's
own builders (card B-18): `tapes.run`, `tokens.run`, `postings.run`, `dim.run`, then
`catalog.build_manifest` (which validates the layout) and `generate_catalog`.

Queries run through the generated catalog and are scored with `catalog.profile_query`: for each key
lookup (sequence by contig, tape of a contig, postings of a token, neighbours of a token) with the key
as a literal and with the keys from a joined table (AGENTS §4.4), the row groups each Parquet scan read
of its total, and for sequence the index row groups the point lookups read (`profile_query`'s `index`).
`rg_per_key` is row groups read divided by the number of keys looked up. The script reports its own
peak RSS (`peak_rss_mb_after`); run it on the pod under a memory scope.
"""

import argparse
import json
import resource
import sys
import time
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("repo", type=Path)
ap.add_argument("out", type=Path)
ap.add_argument("copies", type=int)
ap.add_argument("--shard-bases", type=int, default=2**30)
ap.add_argument("--postings-width", type=int, choices=(0, 1, 2), default=None,
                help="token_bucket width for postings.run (default: its size rule)")
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()

import numpy as np  # noqa: E402
import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402
import yaml  # noqa: E402
from blake3 import blake3  # noqa: E402

from biobox import catalog, dim, fasta, ids  # noqa: E402
from biobox import io as bio  # noqa: E402
from biobox.l0seq import SeqStore, SourceObject, ingest  # noqa: E402
from biobox.l0seq.inputs import write_genome_input  # noqa: E402
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, run as call_genes  # noqa: E402
from biobox.l3 import postings, tokens  # noqa: E402
from biobox.l3 import tapes as tapes_mod  # noqa: E402
from biobox.registry.schemas import table  # noqa: E402

N, REL, V, LIB = a.copies, "rS0.syn", GENE_CALL_VERSION, "synfam-1"
rep: dict = {"copies": N, "shard_bases": a.shard_bases}
T0 = time.perf_counter()


def rss_mb() -> int:
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 if sys.platform == "linux" else 1024 * 1024))


def tick(name: str, t: float, result=None) -> None:
    rep.setdefault("build_s", {})[name] = round(time.perf_counter() - t, 1)
    rep.setdefault("peak_rss_mb_after", {})[name] = rss_mb()
    if result is not None:
        rep.setdefault("builders", {})[name] = result


def puid(base: bytes, c: int) -> bytes:
    return blake3(base + c.to_bytes(4, "big")).digest(length=16)


# --- base: the fixture through the repo's builders ---------------------------------------------
fx = a.repo / "tests/fixtures/micro_s0"
src = yaml.safe_load((fx / "sources.yaml").read_text())
base_root = a.out / "base"
ingest([SourceObject(fx / o["file"], "micro", "1") for o in src["objects"]], base_root, "S0")
base_recs, base_genomes = [], []
for o in src["objects"]:
    r = fasta.read(fx / o["file"])[0]
    cu = ids.contig_uid(r.seq)
    circ = o["topology"] == "circular"
    base_recs.append((o["accession"], r.seq.upper(), cu, circ, o["genetic_code"], o["kind"]))
    base_genomes.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"], circular=(cu,) if circ else ()))
write_genome_input(base_root, "S0", base_genomes)
call_genes(base_root, "S0")
fp = pa.concat_tables(pq.read_table(f, partitioning=None) for f in
                      bio.part_files(base_root / f"l1/feature_protein/gene_call_version={V}/stage=S0"))
occ = pa.concat_tables(pq.read_table(f, columns=["feature_uid", "dist_to_contig_end_bp"], partitioning=None) for f in
                       bio.part_files(base_root / f"l1/protein_occurrence/gene_call_version={V}/stage=S0"))
dist_of = dict(zip(occ.column("feature_uid").to_pylist(), occ.column("dist_to_contig_end_bp").to_pylist()))
genes = {}  # base contig_uid -> list of (start, end, strand, protein_uid, code, dist) in (start, end, strand) order
for r in fp.to_pylist():
    genes.setdefault(r["contig_uid"], []).append((r["start"], r["end"], r["strand"], r["protein_uid"], r["genetic_code"],
                                                  dist_of[r["feature_uid"]]))
for v in genes.values():
    v.sort(key=lambda g: g[:3])
B = len(base_recs)
n_slots = [len(genes.get(cu, [])) for _, _, cu, _, _, _ in base_recs]
G = sum(n_slots)
dprot = sorted({g[3] for v in genes.values() for g in v})  # distinct fixture proteins
pidx = {p: i for i, p in enumerate(dprot)}
P = len(dprot)
fam_kind = np.where(np.arange(P) % 33 == 0, 2, np.where(np.arange(P) % 5 == 4, 0, 1))  # 2 hub, 0 none, 1 own
rep["fixture"] = {"contigs": B, "bp": sum(len(x[1]) for x in base_recs), "slots": G, "distinct_proteins": P}

# --- sequences ---------------------------------------------------------------------------------
t = time.perf_counter()
root = a.out / f"x{N}" / "root"
srcdir = a.out / f"x{N}" / "src"
srcdir.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(a.seed)
ACGT = np.frombuffer(b"ACGT", dtype=np.uint8)
uids = [[b""] * B for _ in range(N)]
objs, fh = [], None
for c in range(N):
    if c % 50 == 0:
        if fh:
            fh.close()
        p = srcdir / f"part-{c // 50:04d}.fna"
        fh = open(p, "wb")
        objs.append(SourceObject(p, "synth", "1"))
    for b, (acc, seq, _, _, _, _) in enumerate(base_recs):
        arr = np.frombuffer(seq, dtype=np.uint8).copy()
        pos = rng.choice(arr.size, size=int(arr.size * 0.03), replace=False)
        arr[pos] = ACGT[rng.integers(0, 4, size=pos.size)]
        s = arr.tobytes()
        uids[c][b] = ids.contig_uid(s)
        fh.write(b">syn%d_%s\n" % (c, acc.encode()))
        fh.write(b"\n".join(s[i : i + 80] for i in range(0, len(s), 80)) + b"\n")
fh.close()
tick("write_fasta", t)
t = time.perf_counter()
tick("ingest", t, ingest(objs, root, "S0", shard_bases=a.shard_bases))
t = time.perf_counter()
gin = [GenomeInput((uids[c][b],), base_recs[b][4], kind=base_recs[b][5], circular=(uids[c][b],) if base_recs[b][3] else ())
       for c in range(N) for b in range(B)]
tick("genome_input", t, write_genome_input(root, "S0", gin))
del gin

# --- gene features (replicated), written in the registry layout ----------------------------------
t = time.perf_counter()
order = sorted((uids[c][b], c, b) for c in range(N) for b in range(B))  # contig_uid order


def contig_genes(u: bytes, c: int, b: int) -> dict:
    gs = genes.get(base_recs[b][2], [])
    return {"contig_uid": [u] * len(gs), "feature_uid": [ids.feature_uid(u, g[0], g[1], g[2], ids.KIND["cds"]) for g in gs],
            "start": [g[0] for g in gs], "end": [g[1] for g in gs], "strand": [g[2] for g in gs],
            "protein_uid": [puid(g[3], c) for g in gs], "genetic_code": [g[4] for g in gs],
            "dist_to_contig_end_bp": [g[5] for g in gs]}


def batches(cols: list[str], extra: dict | None = None):
    for u, c, b in order:
        d = contig_genes(u, c, b)
        if d["start"]:
            d.update({k: [v] * len(d["start"]) for k, v in (extra or {}).items()})
            yield pa.RecordBatch.from_pydict({k: d[k] for k in cols})


spec_f, spec_fp, spec_po = table("feature"), table("feature_protein"), table("protein_occurrence")
bio.write_parts(batches(spec_f.schema.names, {"kind": ids.KIND["cds"]}), spec_f, root / f"l1/feature/producer={V}/stage=S0", 24)
bio.write_parts(batches(spec_fp.schema.names, {"edit_type": "none"}), spec_fp,
                root / f"l1/feature_protein/gene_call_version={V}/stage=S0", 40)
po = pa.concat_tables(pa.Table.from_batches([x]) for x in batches(spec_po.schema.names))
po = po.sort_by([("protein_uid", "ascending"), ("contig_uid", "ascending"), ("start", "ascending"), ("end", "ascending"),
                 ("strand", "ascending")])
bio.write_parts(po, spec_po, root / f"l1/protein_occurrence/gene_call_version={V}/stage=S0", 48)
del po
tick("features (replicated)", t)
t = time.perf_counter()
tick("tapes.run", t, tapes_mod.run(root, "S0"))

# --- synthetic clusters and families -------------------------------------------------------------
t = time.perf_counter()
cl = {k: [] for k in ("member_uid", "level", "cluster_uid", "rep_uid")}
for level, per in (("nr50", 100), ("nr90", 10)):
    for p, base in enumerate(dprot):
        for c in range(N):
            q = c // per
            cl["member_uid"].append(puid(base, c))
            cl["level"].append(level)
            cl["cluster_uid"].append(blake3(level.encode() + p.to_bytes(4, "big") + q.to_bytes(4, "big")).digest(length=16))
            cl["rep_uid"].append(puid(base, q * per))
ct = pa.table({k: pa.array(v, pa.string() if k == "level" else pa.binary()) for k, v in cl.items()})
ct = ct.append_column("cluster_params_hash", pa.array(["synthetic"] * ct.num_rows))
ct = ct.sort_by([("level", "ascending"), ("cluster_uid", "ascending"), ("member_uid", "ascending")])
bio.write_parts(ct, table("cluster_member"), root / f"l2/cluster_member/release={REL}", 37)
del cl, ct
dom = []
for p, base in enumerate(dprot):
    if fam_kind[p] == 0:
        continue
    acc = "HUB" if fam_kind[p] == 2 else f"F{p:06d}"
    for c in range(N):
        dom.append({"protein_uid": puid(base, c), "model_acc": acc, "model_name": acc, "model_len": 100, "seq_bitscore": 50.0,
                    "seq_evalue_zref": 1e-10, "dom_bitscore": 50.0, "dom_ievalue_zref": 1e-10, "env_start": 0, "env_end": 100,
                    "ali_start": 0, "ali_end": 100, "hmm_from": 0, "hmm_to": 100, "cutoff": "gathering", "z": 1e8, "dom_z": 1e8})
dom.sort(key=lambda r: r["protein_uid"])
spec_d = table("l1_domains")
bio.write_parts([pa.RecordBatch.from_pylist(dom[i : i + 122_880], schema=spec_d.schema) for i in range(0, len(dom), 122_880)],
                spec_d, root / f"l1/domains/library={LIB}/stage=S0", 54)
del dom
tick("clusters and families (synthetic)", t)
t = time.perf_counter()
tick("tokens.run", t, tokens.run(root, REL, ["S0"], libraries=[LIB]))  # the release's stages (card B-25)
t = time.perf_counter()
tick("postings.run", t, postings.run(root, REL, ["S0"], prefix_width=a.postings_width))
t = time.perf_counter()
tick("dim.run", t, dim.run(root, REL, ["S0"]))

# --- manifest and catalog ----------------------------------------------------------------------
t = time.perf_counter()
man = catalog.build_manifest(root, REL, ["S0"])
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
rep["manifest_hash"] = catalog.manifest_hash(man)
cat = catalog.generate_catalog(man, a.out / f"x{N}" / "catalog.duckdb")
r0 = rss_mb()
t = time.perf_counter()
con = catalog.open_catalog(cat, threads=4)  # a guest on a shared host: DuckDB would take every core
rep["open_catalog"] = {"seconds": round(time.perf_counter() - t, 2), "contigs_in_locator": N * B}


# --- queries -----------------------------------------------------------------------------------
def fresh() -> None:
    """A new connection: empty block cache and empty index row-group cache, so counts are cold."""
    global con
    con.close()
    con = catalog.open_catalog(cat, threads=4)
    for name, tbl in registered.items():
        con.register(name, tbl)


def run(name: str, sql: str, params=None, keys: int = 1, cold: bool = False) -> dict:
    if cold:
        fresh()
    t = time.perf_counter()
    out = catalog.profile_query(con, sql, params)
    res = out["result"]
    scans = []
    for s in out["scans"]:
        if not s["files"]:
            continue
        rel = Path(s["files"][0]).relative_to(root.resolve()).parts
        tname = rel[1] if rel[0] in ("l0seq", "l1", "l2", "l3", "dim") else rel[0]
        scans.append({"table": tname.removesuffix(".parquet"), "files": s["files_read"],
                      "row_groups": s["row_groups_scanned"], "of": s["row_groups_total"], "rows": s["rows_scanned"],
                      "rg_per_key": round(s["row_groups_scanned"] / keys, 3)})
    idx = out["index"]
    q = {"keys": keys, "seconds": round(time.perf_counter() - t, 3),
         "result": str(res[0]) if len(res) == 1 else f"{len(res)} rows", "scans": scans,
         "parquet_bytes_read": out["bytes_read"], "bgzf": out["bgzf"], "index": idx,
         "index_rg_per_key": round((idx["locator_row_groups_read"] + idx["index_row_groups_read"]) / keys, 3)}
    q["max_rg_per_key"] = max([x["rg_per_key"] for x in scans] + [q["index_rg_per_key"]])
    rep.setdefault("queries", {})[name] = q
    return q


qrng = np.random.default_rng(7)
all_uids = [u for u, _, _ in order]
pick = {k: [all_uids[int(i)] for i in qrng.choice(len(all_uids), size=min(k, len(all_uids)), replace=False)] for k in (1, 10, 1000)}
registered = {f"pick{k}": pa.table({"contig_uid": pa.array(v, pa.binary())}) for k, v in pick.items()}
for name, tbl in registered.items():
    con.register(name, tbl)
one = pick[1][0]
lit = f"unhex('{one.hex()}')"

# 1. sequence by contig: the point lookups of bb_seq (seq_locator, contig_index) and the bgzf blocks
run("sequence by contig: literal key", f"SELECT length(seq_get({lit}, 0, 1000, 1))", cold=True)
for k in (1, 10, 1000):
    run(f"sequence by contig: {k} keys from a joined table", f"SELECT sum(length(seq_get(contig_uid, 0, 1000, 1))) FROM pick{k}",
        keys=k, cold=True)
run("sequence by contig: 1,000 keys from contig_dim joined to the pick table",
    "SELECT sum(length(seq_get(d.contig_uid, 0, 1000, 1))) FROM contig_dim_now() d WHERE d.contig_uid IN (SELECT contig_uid FROM pick1000)",
    keys=1000, cold=True)

# 2. tape of a contig
run("tape of a contig: literal key", f"SELECT sum(n_slots), sum(len(nr90_tokens)) FROM tapes_of([{lit}])")
for k in (1, 10, 1000):
    run(f"tape of a contig: {k} keys from a joined table (tapes_of)",
        f"SELECT sum(n_slots), sum(len(nr90_tokens)) FROM tapes_of((SELECT list(contig_uid) FROM pick{k}))", keys=k)
run("tape of a contig: 1 key, plain join to contig_tapes",
    "SELECT sum(t.n_slots) FROM pick1 p JOIN contig_tapes t ON t.contig_uid = p.contig_uid")

# tokens: a non-hub family pair on adjacent genes of the largest contig, their NR90 tokens, the hub
big = max(range(B), key=lambda b: n_slots[b])
gl = genes[base_recs[big][2]]
gi = next(i for i in range(50, len(gl) - 1) if fam_kind[pidx[gl[i][3]]] == 1 and fam_kind[pidx[gl[i + 1][3]]] == 1
          and pidx[gl[i][3]] != pidx[gl[i + 1][3]])
ka, kb = f"F{pidx[gl[gi][3]]:06d}", f"F{pidx[gl[gi + 1][3]]:06d}"
tok = {r[0]: r[1] for r in con.execute(
    "SELECT key, token_id FROM token_dict WHERE release = catalog_release() AND level = 'family' AND key IN (?, ?, 'HUB')",
    [ka, kb]).fetchall()}
fa, fb, hub = tok[ka], tok[kb], tok["HUB"]
mid = N // 2


def nr90_key(p: int, c: int) -> str:
    return blake3(b"nr90" + p.to_bytes(4, "big") + (c // 10).to_bytes(4, "big")).hexdigest(length=16)


k90a, k90b = nr90_key(pidx[gl[gi][3]], mid), nr90_key(pidx[gl[gi + 1][3]], mid)
na, nb = (con.execute("SELECT token_of('nr90', ?)", [k]).fetchall()[0][0] for k in (k90a, k90b))
rep["tokens"] = {"family_a": [ka, fa], "family_b": [kb, fb], "hub": hub, "nr90_a": [k90a, na], "nr90_b": [k90b, nb]}
tok100 = [r[0] for r in con.execute("SELECT token_id FROM token_dict WHERE release = catalog_release() AND level = 'nr90' "
                                    "ORDER BY hash(token_id), token_id LIMIT 100").fetchall()]
registered["picktok100"] = pa.table({"token_id": pa.array(tok100, pa.uint64())})
con.register("picktok100", registered["picktok100"])

# 3. postings of a token
run("postings of a token: literal key, family (COPIES loci)", "SELECT count(*), sum(ordinal) FROM token_rows(?)", [fa])
run("postings of a token: literal key, nr90 (10 loci)", "SELECT count(*), sum(ordinal) FROM token_rows(?)", [na])
run("postings of a token: key from token_dict (token_of), nr90",
    "SELECT count(*), sum(ordinal) FROM token_rows(token_of('nr90', ?))", [k90a])
run("postings of a token: key from token_dict, plain join, nr90",
    "SELECT count(*), sum(p.ordinal) FROM (SELECT token_id FROM token_dict WHERE release = catalog_release() "
    "AND level = 'nr90' AND key = ?) d JOIN token_postings p ON p.token_id = d.token_id "
    "WHERE p.release = catalog_release() AND p.tape_version = catalog_tape_version()", [k90a])
run("postings of a token: 100 nr90 keys from a joined table, plain join",
    "SELECT count(*), sum(p.ordinal) FROM picktok100 k JOIN token_postings p ON p.token_id = k.token_id "
    "WHERE p.release = catalog_release() AND p.tape_version = catalog_tape_version()", keys=100)

# 4. neighbours of a token
run("neighbours: Q1 family pair, literal keys", "SELECT q1_count(?, ?, 3, NULL)", [fa, fb])
run("neighbours: Q1 nr90 pair, literal keys", "SELECT q1_count(?, ?, 3, NULL)", [na, nb])
run("neighbours: Q1 nr90 pair, keys from token_dict (token_of)",
    "SELECT q1_count(token_of('nr90', ?), token_of('nr90', ?), 3, NULL)", [k90a, k90b])
run("neighbours: Q4 family pair instances, literal keys", "SELECT count(*), sum(delta_bp) FROM pairs_of(?, ?, 10)", [fa, fb])
run("neighbours: Q4 nr90 pair instances, keys from token_dict (token_of)",
    "SELECT count(*), sum(delta_bp) FROM pairs_of(token_of('nr90', ?), token_of('nr90', ?), 10)", [k90a, k90b])
run("neighbours: Q3 partners of one family, literal key", "SELECT count(*), sum(n_pairs) FROM partners_of(?, 10)", [fa])
run("neighbours: Q3 partners of one nr90, key from token_dict (token_of)",
    "SELECT count(*), sum(n_pairs) FROM partners_of(token_of('nr90', ?), 10)", [k90a])
run("neighbours: Q1 hub anchor x family partner (hub, not a key lookup)", "SELECT q1_count(?, ?, 10, NULL)", [hub, fa])
run("neighbours: Q3 partners of the hub (hub, not a key lookup)", "SELECT count(*), sum(n_pairs) FROM partners_of(?, 10)", [hub])

# 5. bulk reads, as in 2026-11-09: 1,000 neighbourhoods of +-10 genes and whole-contig exports
anchors = []
for _ in range(1000):
    c, b = int(qrng.integers(N)), int(qrng.integers(B))
    gs = genes[base_recs[b][2]]
    i = int(qrng.integers(len(gs)))
    lo, hi = max(i - 10, 0), min(i + 10, len(gs) - 1)
    anchors.append((uids[c][b], gs[lo][0], min(gs[hi][1], len(base_recs[b][1]))))
registered["anchors"] = pa.table({"contig_uid": pa.array([x[0] for x in anchors], pa.binary()),
                                  "s": pa.array([x[1] for x in anchors], pa.int64()), "e": pa.array([x[2] for x in anchors], pa.int64())})
run("1,000 neighbourhoods: sequence (+-10 genes)", "SELECT sum(length(seq_slice(contig_uid, s, e))) FROM anchors", keys=1000,
    cold=True)["bp_returned"] = sum(x[2] - x[1] for x in anchors)
run("1,000 neighbourhoods: tapes and tokens",
    "SELECT sum(n_slots), sum(len(nr90_tokens)) FROM tapes_of((SELECT list(DISTINCT contig_uid) FROM anchors))", keys=1000)
store = SeqStore(root)
for k in (1, 50, 500):
    k = min(k, N * B)
    pk = [uids[int(qrng.integers(N))][int(qrng.integers(B))] for _ in range(k)]
    registered["g"] = pa.table({"contig_uid": pa.array(pk, pa.binary())})
    q = run(f"export proxy: {k} whole contigs", "SELECT sum(length(seq_slice(contig_uid, 0, contig_len_of(contig_uid)))) FROM g",
            keys=k, cold=True)
    q["shards_touched"] = len({store.shard_of(store.default_digest(u)) for u in pk})
    q["shards_total"] = len(shards)
    q["bp_returned"] = sum(store.length_of(store.default_digest(u)) for u in pk)
store.close()

keyq = {k: v for k, v in rep["queries"].items() if not k.startswith(("1,000 neighbourhoods", "export proxy")) and "hub" not in k}
rep["key_lookups_max_rg_per_key"] = max(v["max_rg_per_key"] for v in keyq.values())
rep["key_lookups_over_2_rg_per_key"] = sorted(k for k, v in keyq.items() if v["max_rg_per_key"] > 2)
rep["peak_rss_mb_queries"] = rss_mb()
rep["rss_mb_before_open"] = r0
rep["total_s"] = round(time.perf_counter() - T0, 1)
print(json.dumps(rep, indent=1, default=str))
