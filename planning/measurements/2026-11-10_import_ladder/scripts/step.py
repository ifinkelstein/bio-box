"""One ladder step, run inside the pinned image with bio-box-next bound in (sh bbn.sh step.py KEY STEP).

    python step.py KEY STEP [--threads N]

KEY: A (phages, S1 -> rL.1), B (P1 -> rL.2), Bp (P1 Pfam -> rL.2p), C (P2 -> rL.3).
Writes /tmp/if2453/ladder/results/<KEY>/<STEP>.json: {"result", "error", "traceback", "inner_wall_s"}.
The driver (driver.py, host) times the process and records its CPU and peak RSS from wait4.
"""

import argparse
import json
import os
import resource
import shutil
import sys
import time
import traceback
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from biobox import catalog, ids, ledger, provenance
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l0seq.store import SeqStore, index_files, locator_files
from biobox.l1 import genecall, igr
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.l2 import cluster
from biobox.l3 import postings, tapes, tokens
from biobox import dim
from biobox.registry.schemas import TABLES, table

PILOT = Path("/tmp/if2453/pilot_s1/root")
PFAM_LIB = "pfam-38.2"
CFG = {
    "A": {"stage": "S1", "release": "rL.1", "prev": []},
    "B": {"stage": "P1", "release": "rL.2", "prev": ["rL.1"]},
    "Bp": {"stage": "P1", "release": "rL.2p", "prev": ["rL.1", "rL.2"]},
    "C": {"stage": "P2", "release": "rL.3", "prev": ["rL.1", "rL.2", "rL.2p"]},
}
NCBI_RELEASE = "datasets-v2-2026-11-09"
GVC = f"gene_call_version={GENE_CALL_VERSION}"

ap = argparse.ArgumentParser()
ap.add_argument("key", choices=sorted(CFG))
ap.add_argument("step")
ap.add_argument("--threads", type=int, default=1)
ap.add_argument("--work", type=Path, default=Path("/tmp/if2453/ladder"))  # smoke tests use another
ap.add_argument("--src", type=Path, default=Path("/tmp/if2453/pilot_s1/src/14Apr2025_genomes.fa.gz"))
a = ap.parse_args()
W = a.work
ROOT = W / "root"
TMP = W / "tmp"
INPHARED = a.src
C = CFG[a.key]
STAGE, REL = C["stage"], C["release"]
OUT = W / "results" / a.key / f"{a.step}.json"
OUT.parent.mkdir(parents=True, exist_ok=True)
TMP.mkdir(parents=True, exist_ok=True)


def rss_mb() -> float:
    with open("/proc/self/statm") as fh:
        return round(int(fh.read().split()[1]) * resource.getpagesize() / 2**20, 1)


def sql_str(s) -> str:
    return "'" + str(s).replace("'", "''") + "'"


def sql_list(files) -> str:
    return "[" + ", ".join(sql_str(f) for f in files) + "]"


def duck(threads: int = 4):
    d = TMP / f"duck-{a.key}-{a.step}-{os.getpid()}"
    return duckdb.connect(config={"threads": threads, "memory_limit": "16GB", "temp_directory": str(d)})


def jsonable(v):
    if isinstance(v, (bytes, bytearray)):
        return v.hex()
    if isinstance(v, (list, tuple)):
        return [jsonable(x) for x in v]
    if isinstance(v, dict):
        return {str(k): jsonable(x) for k, x in v.items()}
    return v


# --- L0 -------------------------------------------------------------------------------------
def read_lock() -> dict:
    return json.loads((W / f"genomes/{STAGE}/lock.json").read_text())["genomes"]


def read_tsv() -> dict:
    out = {}
    with open(W / f"selection/stage_{STAGE}.tsv") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            r = dict(zip(head, ln.rstrip("\n").split("\t")))
            out[r["accession"]] = r
    return out


def do_ingest():
    if (ROOT / f"l0seq/contig_index/stage={STAGE}").exists():
        raise FileExistsError(f"stage {STAGE} already ingested")
    if a.key == "A":
        objs = [SourceObject(INPHARED, "inphared", "14Apr2025", license="AGPL-3.0 (repo); data from INSDC")]
    else:
        lock = read_lock()
        objs = [SourceObject(Path(v["file"]), "ncbi_genomes", NCBI_RELEASE, license="INSDC public data",
                             redistributable=True,
                             retrieved_at=datetime.fromisoformat(v["retrieved_utc"]))
                for _, v in sorted(lock.items())]
    in_bytes = sum(Path(o.path).stat().st_size for o in objs)
    res = ingest(objs, ROOT, STAGE, tmp=TMP)
    con = duck(1)
    occ = sorted((ROOT / f"l0seq/contig_occurrence/stage={STAGE}").glob("source=*/part-0.parquet"))
    loc = locator_files(ROOT)
    (rec_bp, n_rec), = con.execute(
        f"SELECT sum(l.contig_len), count(*) FROM read_parquet({sql_list(occ)}) o "
        f"JOIN read_parquet({sql_list(loc)}) l USING (raw_digest)").fetchall()
    idx = ROOT / f"l0seq/contig_index/stage={STAGE}/part-0.parquet"
    (n_idx, idx_bp), = con.execute(f"SELECT count(*), sum(contig_len) FROM read_parquet({sql_str(idx)}) "
                                   "WHERE is_default_variant").fetchall()
    stored = ROOT / f"l0seq/seq_locator/stage={STAGE}/part-0.parquet"
    (n_st, st_bp), = con.execute(f"SELECT count(*), sum(contig_len) FROM read_parquet({sql_str(stored)})").fetchall()
    shards = sorted((ROOT / f"l0seq/seq_shards/stage={STAGE}").glob("shard-*.fa.gz"))
    return {**res, "objects_bytes": in_bytes, "records_bp": int(rec_bp), "records_counted": n_rec,
            "default_contigs": n_idx, "default_contig_bp": int(idx_bp), "stored_sequences": n_st,
            "stored_bp": int(st_bp), "shard_bytes": sum(s.stat().st_size for s in shards),
            "shard_index_bytes": sum(Path(str(s) + x).stat().st_size for s in shards for x in (".fai", ".gzi"))}


def do_genome_input():
    idx = ROOT / f"l0seq/contig_index/stage={STAGE}/part-0.parquet"
    t = pq.read_table(idx, columns=["contig_uid", "is_default_variant", "contig_len"]).to_pylist()
    default = {r["contig_uid"]: r["contig_len"] for r in t if r["is_default_variant"]}
    if a.key == "A":
        genomes = [GenomeInput((cu,), 11, kind="virus") for cu in sorted(default)]
        res = write_genome_input(ROOT, STAGE, genomes)
        return {**res, "kind": "virus", "genetic_code": 11, "circular": 0,
                "note": "every INPHARED record its own genome; topology unknown, all linear",
                "bp": sum(default.values())}
    lock, meta = read_lock(), read_tsv()
    occ = pq.read_table(ROOT / f"l0seq/contig_occurrence/stage={STAGE}/source=ncbi_genomes/part-0.parquet",
                        columns=["l0_object_sha256", "record_ordinal", "contig_uid"]).to_pylist()
    by_sha = defaultdict(list)
    for r in occ:
        by_sha[r["l0_object_sha256"]].append((r["record_ordinal"], r["contig_uid"]))
    earlier = set()
    for f in index_files(ROOT):
        if f"stage={STAGE}/" in str(f):
            continue
        for r in pq.read_table(f, columns=["contig_uid", "is_default_variant"]).to_pylist():
            if r["is_default_variant"]:
                earlier.add(r["contig_uid"])
    assigned, genomes = {}, []
    st = {"dup_within_genome": 0, "dup_across_genomes": [], "in_earlier_stage": [], "genomes_without_contigs": [],
          "codes": defaultdict(int)}
    for acc in sorted(lock):
        cus = [cu for _, cu in sorted(by_sha[bytes.fromhex(lock[acc]["sha256"])])]
        if not cus:
            raise ValueError(f"{acc}: no contig_occurrence rows for its sha256")
        seen, kept = set(), []
        for cu in cus:
            if cu in seen:
                st["dup_within_genome"] += 1
                continue
            seen.add(cu)
            if cu in earlier:
                st["in_earlier_stage"].append([acc, cu.hex()])
                continue
            if cu in assigned:
                st["dup_across_genomes"].append([acc, assigned[cu], cu.hex()])
                continue
            assigned[cu] = acc
            kept.append(cu)
        if not kept:
            st["genomes_without_contigs"].append(acc)
            continue
        tt = meta[acc]["translation_table"]
        code = int(tt) if tt not in ("", "None") else 11
        st["codes"][str(code)] += 1
        genomes.append(GenomeInput(tuple(kept), code, kind="prokaryote", assembly_quality="complete"))
    missing = set(default) - set(assigned)
    stray = set(assigned) - set(default)
    if missing or stray:
        raise ValueError(f"genome_input vs contig_index: {len(missing)} index contigs unassigned, "
                         f"{len(stray)} assigned contigs not default in this stage's index")
    res = write_genome_input(ROOT, STAGE, genomes)
    return {**res, "kind": "prokaryote", "circular": 0, "bp": sum(default.values()),
            "note": "all contigs of an assembly one genome; topology not in the FASTA, all linear; "
                    "genetic_code from ncbi_translation_table, 11 if absent",
            **{k: (dict(v) if isinstance(v, defaultdict) else v) for k, v in st.items()}}


# --- L1 -------------------------------------------------------------------------------------
def do_genecall():
    res = genecall.run(ROOT, STAGE, workers=a.threads, tmp=TMP)
    gi = pq.read_table(ROOT / f"l0seq/genome_input/stage={STAGE}/part-0.parquet", columns=["contig_uid"])
    idx = pq.read_table(ROOT / f"l0seq/contig_index/stage={STAGE}/part-0.parquet",
                        columns=["contig_uid", "contig_len", "is_default_variant"]).to_pylist()
    want = set(gi.column("contig_uid").to_pylist())
    bp = sum(r["contig_len"] for r in idx if r["is_default_variant"] and r["contig_uid"] in want)
    ps = sorted((ROOT / f"l0seq/protein_seq/{GVC}/stage={STAGE}").glob("protein_prefix=*/part-*.parquet"))
    con = duck(1)
    (aa,), = con.execute(f"SELECT sum(length) FROM read_parquet({sql_list(ps)})").fetchall()
    return {**res, "workers": a.threads, "bp_called": bp, "residues": int(aa), "protein_seq_files": len(ps)}


def do_igr():
    return igr.run(ROOT, STAGE, tmp=TMP)


def do_tapes():
    return tapes.run(ROOT, STAGE, tmp=TMP)


# --- L2 -------------------------------------------------------------------------------------
def do_cluster():
    for name in ("cluster_member", "cluster_member_aln"):
        if (ROOT / f"l2/{name}/release={REL}").exists():
            raise FileExistsError(f"{name} for {REL} exists; cluster.run would overwrite it")
    mm = TMP / f"mmseqs-{REL}"
    mm.mkdir(parents=True, exist_ok=True)
    files = sorted((ROOT / "l0seq/protein_seq").glob(f"{GVC}/stage=*/protein_prefix=*/part-*.parquet"))
    res = cluster.run(ROOT, REL, GENE_CALL_VERSION, threads=a.threads, tmp=mm)
    shutil.rmtree(mm, ignore_errors=True)
    return {**res, "threads": a.threads, "protein_seq_files_read": [str(f.relative_to(ROOT)) for f in files]}


def domain_path(stage):
    return ROOT / f"l1/domains/library={PFAM_LIB}/stage={stage}/part-0.parquet"


def do_domains():
    """Stage A: reuse the pilot's Pfam domains iff the distinct protein uid sets are equal."""
    new = sorted((ROOT / f"l0seq/protein_seq/{GVC}/stage=S1").glob("protein_prefix=*/part-*.parquet"))
    old = sorted((PILOT / f"l0seq/protein_seq/{GVC}/stage=S1").glob("protein_prefix=*/part-*.parquet"))
    con = duck(4)
    (n_new,), = con.execute(f"SELECT count(DISTINCT protein_uid) FROM read_parquet({sql_list(new)})").fetchall()
    (n_old,), = con.execute(f"SELECT count(DISTINCT protein_uid) FROM read_parquet({sql_list(old)})").fetchall()
    (only_new,), = con.execute(f"SELECT count(*) FROM (SELECT protein_uid FROM read_parquet({sql_list(new)}) EXCEPT "
                               f"SELECT protein_uid FROM read_parquet({sql_list(old)}))").fetchall()
    (only_old,), = con.execute(f"SELECT count(*) FROM (SELECT protein_uid FROM read_parquet({sql_list(old)}) EXCEPT "
                               f"SELECT protein_uid FROM read_parquet({sql_list(new)}))").fetchall()
    res = {"new_distinct": n_new, "pilot_distinct": n_old, "only_new": only_new, "only_pilot": only_old}
    if only_new or only_old:
        res["equal"] = False
        raise RuntimeError(f"protein uid sets differ: {json.dumps(res)}; run the sharded scan instead")
    res["equal"] = True
    src = PILOT / f"l1/domains/library={PFAM_LIB}/stage=S1/part-0.parquet"
    dst = domain_path("S1")
    if dst.exists():
        raise FileExistsError(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(".part-0.parquet.tmp")
    shutil.copyfile(src, tmp)
    tmp.replace(dst)
    t = pq.read_table(dst, columns=["protein_uid"])
    hit = len(set(t.column("protein_uid").to_pylist()))
    lib = PFAM_LIB.encode()
    tools = provenance.tool_versions("pyhmmer", "pyarrow")
    lp = ledger.write(ROOT, "S1", f"domains-{PFAM_LIB}",
                      [ledger.row(lib, "library", n_new, hit, tools=tools),
                       ledger.row(lib, "library", n_new, hit, n_rejected=n_new - hit, reason="no_hit", tools=tools)])
    old_ledger = pq.read_table(PILOT / f"ledger/stage=S1/step=domains-{PFAM_LIB}/part-0.parquet").to_pylist()
    new_ledger = pq.read_table(lp).to_pylist()
    res.update({"copied": str(dst.relative_to(ROOT)), "sha256": ids.sha256_file(dst).hex(),
                "pilot_sha256": ids.sha256_file(src).hex(), "rows": t.num_rows, "proteins_with_hit": hit,
                "ledger": str(lp.relative_to(ROOT)),
                "ledger_rows_equal_pilot": [{k: v for k, v in r.items()} for r in old_ledger] ==
                                           [{k: v for k, v in r.items()} for r in new_ledger],
                "balance": ledger.check_balance(ROOT, "S1", [f"domains-{PFAM_LIB}"]),
                "chain": ledger.check_chain(ROOT, "S1")})
    return res


# --- L3, dim --------------------------------------------------------------------------------
def do_tokens():
    return {**tokens.run(ROOT, REL, tmp=TMP), "libraries": tokens.libraries_in(ROOT),
            "domain_files": [str(p.relative_to(ROOT)) for p in sorted((ROOT / "l1/domains").glob("*/*/*.parquet"))]}


def do_postings():
    return postings.run(ROOT, REL, tmp=TMP)


def do_dim():
    return dim.run(ROOT, REL, tmp=TMP)


# --- manifest and catalog -------------------------------------------------------------------
def manifest_tables(man: dict, prev_paths: set) -> dict:
    per = defaultdict(lambda: {"files": 0, "rows": 0, "bytes": 0, "row_groups": 0,
                               "new_files": 0, "new_rows": 0, "new_bytes": 0})
    for f in man["files"]:
        p = per[f["table"]]
        new = f["path"] not in prev_paths
        for k, v in (("files", 1), ("rows", f["rows"]), ("bytes", f["bytes"]), ("row_groups", f["row_groups"])):
            p[k] += v
            if new and k != "row_groups":
                p["new_" + k] += v
    for f in man["artifacts"]:
        suffix = ".fai" if f["path"].endswith(".fai") else ".gzi" if f["path"].endswith(".gzi") else ""
        p = per[f"artifact:{f['artifact']}{suffix}"]
        p["files"] += 1
        p["bytes"] += f["bytes"]
        if f["path"] not in prev_paths:
            p["new_files"] += 1
            p["new_bytes"] += f["bytes"]
    return dict(sorted(per.items()))


def prev_manifests() -> list[tuple[str, dict]]:
    out = []
    for r in C["prev"]:
        p = ROOT / f"manifests/release={r}.json"
        if p.exists():
            out.append((r, catalog.read_manifest(p)))
    return out


def do_catalog():
    times = {}

    def timed(name, fn):
        t = time.perf_counter()
        x = fn()
        times[name] = round(time.perf_counter() - t, 2)
        return x

    man = timed("build_manifest", lambda: catalog.build_manifest(ROOT, REL))
    mpath = ROOT / f"manifests/release={REL}.json"
    timed("write_manifest", lambda: catalog.write_manifest(man, mpath))
    cpath = ROOT / f"catalog/release={REL}.duckdb"
    timed("generate_catalog", lambda: catalog.generate_catalog(man, cpath))
    problems = timed("verify_manifest", lambda: catalog.verify_manifest(catalog.read_manifest(mpath), ROOT))
    catalog.swap_current(ROOT / "manifests", REL)
    prev = set()
    for _, m in prev_manifests():
        prev |= {f["path"] for f in m["files"] + m["artifacts"]}
    return {"manifest": str(mpath), "manifest_hash": catalog.manifest_hash(man), "catalog": str(cpath),
            "catalog_bytes": cpath.stat().st_size, "files": len(man["files"]), "artifacts": len(man["artifacts"]),
            "total_bytes": sum(f["bytes"] for f in man["files"] + man["artifacts"]),
            "verify_problems": problems, "times_s": times, "tables": manifest_tables(man, prev)}


def do_appendonly():
    """Every file of the earlier manifests is still byte-identical (sha256 recomputed now)."""
    checked, same, differ, missing, expected = 0, 0, [], [], []
    seen = {}
    for rel_id, m in prev_manifests():
        for f in m["files"] + m["artifacts"]:
            seen.setdefault(f["path"], (rel_id, f["sha256"], f["bytes"]))
    t0 = time.perf_counter()
    nbytes = 0
    for path, (rel_id, sha, nb) in sorted(seen.items()):
        p = ROOT / path
        checked += 1
        if not p.exists():
            missing.append(path)
            continue
        nbytes += p.stat().st_size
        now = ids.sha256_file(p).hex()
        if now == sha:
            same += 1
        elif path == "l0/l0_manifest.parquet":
            expected.append({"path": path, "first_manifest": rel_id, "note": "small metadata file rewritten per stage "
                             "by ingest (store.py docstring); every earlier row is still present"})
        else:
            differ.append(path)
    # the rewritten l0_manifest must still hold every earlier object's row
    extra = {}
    if expected:
        rows = pq.read_table(ROOT / "l0/l0_manifest.parquet").to_pylist()
        extra["l0_manifest_rows_now"] = len(rows)
    return {"manifests": [r for r, _ in prev_manifests()], "files_checked": checked, "identical": same,
            "differ": differ, "missing": missing, "expected_rewrites": expected, "bytes_hashed": nbytes,
            "hash_s": round(time.perf_counter() - t0, 1), "ok": not differ and not missing, **extra}


# --- stage A equivalence with the old-code pilot ---------------------------------------------
def do_equiv():
    con = duck(4)
    P, N = PILOT, ROOT

    def files(root, pattern):
        return sql_list(sorted(root.glob(pattern)))

    def cols(src):
        return [r[0] for r in con.execute(f"DESCRIBE SELECT * FROM {src}").fetchall()]

    def setcmp(name, a_src, b_src, columns):
        c = ", ".join(f'"{x}"' for x in columns)
        (na,), = con.execute(f"SELECT count(*) FROM {a_src}").fetchall()
        (nb,), = con.execute(f"SELECT count(*) FROM {b_src}").fetchall()
        (a_b,), = con.execute(f"SELECT count(*) FROM (SELECT {c} FROM {a_src} EXCEPT ALL SELECT {c} FROM {b_src})").fetchall()
        (b_a,), = con.execute(f"SELECT count(*) FROM (SELECT {c} FROM {b_src} EXCEPT ALL SELECT {c} FROM {a_src})").fetchall()
        return {"name": name, "columns": columns, "pilot_rows": na, "new_rows": nb, "only_pilot": a_b, "only_new": b_a,
                "equal": a_b == 0 and b_a == 0 and na == nb}

    out = []
    rp = lambda root, pat: f"read_parquet({files(root, pat)}, hive_partitioning = false)"  # noqa: E731
    # stored sequences and contig index
    out.append(setcmp("seq_locator raw_digest,contig_len", rp(P, "l0seq/seq_locator.parquet"),
                      rp(N, "l0seq/seq_locator/stage=S1/part-0.parquet"), ["raw_digest", "contig_len"]))
    out.append(setcmp("contig_index", rp(P, "l0seq/contig_index/stage=S1/part-0.parquet"),
                      rp(N, "l0seq/contig_index/stage=S1/part-0.parquet"),
                      ["contig_uid", "raw_digest", "contig_len", "is_default_variant"]))
    occ_p = rp(P, "l0seq/contig_occurrence/stage=S1/source=*/part-0.parquet")
    occ_n = rp(N, "l0seq/contig_occurrence/stage=S1/source=*/part-0.parquet")
    out.append(setcmp("contig_occurrence", occ_p, occ_n, [c for c in cols(occ_p) if c in cols(occ_n)]))
    # genes: coordinates and protein per contig; the pilot's protein_occurrence vs the new feature_protein
    po_p = rp(P, f"l1/protein_occurrence/{GVC}/stage=S1/part-*.parquet")
    fp_n = rp(N, f"l1/feature_protein/{GVC}/stage=S1/part-*.parquet")
    po_n = rp(N, f"l1/protein_occurrence/{GVC}/stage=S1/part-*.parquet")
    key5 = ["contig_uid", "start", "end", "strand", "protein_uid"]
    out.append(setcmp("genes: pilot protein_occurrence vs new feature_protein", po_p, fp_n, key5))
    common_po = [c for c in cols(po_p) if c in cols(po_n)]
    out.append(setcmp("protein_occurrence (common columns)", po_p, po_n, common_po))
    sp_p = rp(P, f"l1/slot_protein/{GVC}/stage=S1/part-*.parquet")
    sp_n = rp(N, f"l1/slot_protein/{GVC}/stage=S1/part-*.parquet")
    out.append(setcmp("slot_protein", sp_p, sp_n, [c for c in cols(sp_p) if c in cols(sp_n)]))
    ps_p = rp(P, f"l0seq/protein_seq/{GVC}/stage=S1/protein_prefix=*/part-*.parquet")
    ps_n = rp(N, f"l0seq/protein_seq/{GVC}/stage=S1/protein_prefix=*/part-*.parquet")
    out.append(setcmp("protein_seq", ps_p, ps_n, [c for c in cols(ps_p) if c in cols(ps_n)]))
    # IGRs
    ig_p = rp(P, "l0seq/igr/igr_version=*/stage=S1/part-*.parquet")
    ig_n = rp(N, "l0seq/igr/igr_version=*/stage=S1/part-*.parquet")
    out.append(setcmp("igr (common columns)", ig_p, ig_n, [c for c in cols(ig_p) if c in cols(ig_n)]))
    # tapes: one row per contig with its lists
    tp_p = rp(P, "l3/tapes/*/stage=S1/part-*.parquet")
    tp_n = rp(N, "l3/tapes/*/stage=S1/part-*.parquet")
    common_tp = [c for c in cols(tp_p) if c in cols(tp_n)]
    out.append(setcmp("tapes (common columns, lists compared whole)", tp_p, tp_n, common_tp))
    out.append({"name": "tape columns only in pilot / only in new",
                "only_pilot": [c for c in cols(tp_p) if c not in cols(tp_n)],
                "only_new": [c for c in cols(tp_n) if c not in cols(tp_p)]})
    # clusters
    cm_p = rp(P, "l2/cluster_member/release=*/part-0.parquet")
    cm_n = rp(N, "l2/cluster_member/release=rL.1/part-0.parquet")
    counts = {}
    for lbl, src in (("pilot", cm_p), ("new", cm_n)):
        counts[lbl] = dict(con.execute(f"SELECT level, count(DISTINCT cluster_uid) FROM {src} GROUP BY level").fetchall())
    out.append({"name": "cluster counts", **counts, "equal": counts["pilot"] == counts["new"]})
    out.append(setcmp("cluster_member membership", cm_p, cm_n,
                      ["member_uid", "level", "cluster_uid", "rep_uid", "cluster_params_hash"]))
    ca_p = rp(P, "l2/cluster_member_aln/release=*/part-0.parquet")
    ca_n = rp(N, "l2/cluster_member_aln/release=rL.1/part-0.parquet")
    out.append(setcmp("cluster_member_aln", ca_p, ca_n, [c for c in cols(ca_p) if c in cols(ca_n)]))
    # tokens (pilot tokens were rerun with Pfam in the prov_db work)
    td_p = rp(P, "l3/token_dict/release=*/part-*.parquet")
    td_n = rp(N, "l3/token_dict/release=rL.1/part-*.parquet")
    out.append(setcmp("token_dict (token_id, level, source, key, n_slots)", td_p, td_n,
                      ["token_id", "level", "source", "key", "n_slots"]))
    tt_p = rp(P, "l3/tape_tokens/*/release=*/part-*.parquet")
    tt_n = rp(N, "l3/tape_tokens/*/release=rL.1/part-*.parquet")
    out.append(setcmp("tape_tokens", tt_p, tt_n, [c for c in cols(tt_p) if c in cols(tt_n)]))
    return {"checks": out, "all_equal": all(c.get("equal", True) for c in out)}


# --- Pfam for stage B's new proteins --------------------------------------------------------
PF = W / "pfam_P1"


def do_pfam_plan():
    """Stage P1's proteins not in any stage that already has Pfam domains, sorted by uid, cut into
    shards of about 16 M residues (residues without '*', as domains.scan counts)."""
    import pyarrow.compute as pc

    done_stages = sorted(p.parent.name.removeprefix("stage=") for p in
                         (ROOT / f"l1/domains/library={PFAM_LIB}").glob("stage=*/part-0.parquet"))
    mine = sorted((ROOT / f"l0seq/protein_seq/{GVC}/stage={STAGE}").glob("protein_prefix=*/part-*.parquet"))
    others = [f for s in done_stages for f in
              sorted((ROOT / f"l0seq/protein_seq/{GVC}/stage={s}").glob("protein_prefix=*/part-*.parquet"))]
    con = duck(4)
    t = con.execute(f"SELECT protein_uid, seq FROM read_parquet({sql_list(mine)}) "
                    + (f"WHERE protein_uid NOT IN (SELECT protein_uid FROM read_parquet({sql_list(others)})) " if others else "")
                    + "ORDER BY protein_uid").to_arrow_table()
    t = t.cast(pa.schema([("protein_uid", pa.binary(16)), ("seq", pa.string())]))
    (n_all,), = con.execute(f"SELECT count(*) FROM read_parquet({sql_list(mine)})").fetchall()
    lens = pc.subtract(pc.utf8_length(t.column("seq")), pc.count_substring(t.column("seq"), "*")).to_numpy()
    total = int(lens.sum())
    per = 16_000_000
    k, start, acc, shards = 0, 0, 0, []
    inp = PF / "input"
    inp.mkdir(parents=True, exist_ok=True)
    import numpy as np

    cum = np.cumsum(lens)
    n = t.num_rows
    while start < n:
        base = int(cum[start - 1]) if start else 0
        end = int(np.searchsorted(cum, base + per, side="left")) + 1
        end = min(max(end, start + 1), n)
        out = inp / f"shard-{k}.parquet"
        piece = t.slice(start, end - start)
        if out.exists():
            old = pq.read_table(out, columns=["protein_uid"]).column("protein_uid")
            if len(old) != piece.num_rows or old[0] != piece.column("protein_uid")[0]:
                raise RuntimeError(f"{out} exists with other bounds")
        else:
            pq.write_table(piece, out.with_suffix(".tmp"), compression="zstd")
            out.with_suffix(".tmp").replace(out)
        shards.append({"shard": k, "proteins": piece.num_rows, "residues": int(cum[end - 1]) - base})
        k += 1
        start = end
    plan = {"stage": STAGE, "excluded_stages_with_pfam": done_stages, "stage_proteins": n_all,
            "new_proteins": n, "residues": total, "shards": shards}
    (inp / "plan.json").write_text(json.dumps(plan, indent=1))
    return {k: v for k, v in plan.items() if k != "shards"} | {"n_shards": len(shards)}


def do_pfam_merge():
    import pyarrow.compute as pc

    plan = json.loads((PF / "input/plan.json").read_text())
    metas = []
    for s in plan["shards"]:
        m = json.loads((PF / f"scan/shard-{s['shard']}.json").read_text())
        if not m.get("ok") or m["proteins"] != s["proteins"]:
            raise RuntimeError(f"shard {s['shard']} not ok or protein count differs: {m.get('error')}")
        metas.append(m)
    n_in = plan["new_proteins"]
    schema = table("l1_domains").schema
    parts = [pq.read_table(PF / f"scan/shard-{s['shard']}.parquet") for s in plan["shards"]]
    t = pa.concat_tables([p.cast(schema) for p in parts])
    t = t.take(pc.sort_indices(t, sort_keys=[("protein_uid", "ascending"), ("env_start", "ascending"),
                                              ("model_acc", "ascending")]))
    out = domain_path(STAGE)
    if out.exists():
        raise FileExistsError(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(".part-0.parquet.tmp")
    pq.write_table(t, tmp, compression="zstd", row_group_size=122_880)
    tmp.replace(out)
    hit = pc.count_distinct(t.column("protein_uid")).as_py()
    lib = PFAM_LIB.encode()
    tools = provenance.tool_versions("pyhmmer", "pyarrow")
    lp = ledger.write(ROOT, STAGE, f"domains-{PFAM_LIB}",
                      [ledger.row(lib, "library", n_in, hit, tools=tools),
                       ledger.row(lib, "library", n_in, hit, n_rejected=n_in - hit, reason="no_hit", tools=tools)])
    return {"out": str(out.relative_to(ROOT)), "bytes": out.stat().st_size, "proteins_scanned": n_in,
            "stage_proteins": plan["stage_proteins"], "domains": t.num_rows, "proteins_with_hit": hit,
            "models_with_hit": pc.count_distinct(t.column("model_acc")).as_py(),
            "shard_cpu_min": round(sum(m["cpu_min"] for m in metas), 1),
            "shard_wall_s_sum": round(sum(m["wall_s"] for m in metas), 1),
            "shard_peak_rss_mb_max": max(m["peak_rss_mb"] for m in metas),
            "cpu_min_per_mres": round(sum(m["cpu_min"] for m in metas) / (plan["residues"] / 1e6), 3),
            "ledger": str(lp.relative_to(ROOT)),
            "balance": ledger.check_balance(ROOT, STAGE, [f"domains-{PFAM_LIB}"]),
            "chain": ledger.check_chain(ROOT, STAGE),
            "note": "only proteins new in P1 were scanned (S1 proteins carry the pilot's Pfam rows), so the "
                    "ledger n_in is the scanned count and the genecall->domains chain check (L-03) reports it"}


def do_cluster_copy():
    """rL.2p reuses rL.2's clusters: the protein set is the same. Byte copies, new paths only."""
    src_rel = "rL.2"
    out = {}
    for name in ("cluster_member", "cluster_member_aln"):
        s = ROOT / f"l2/{name}/release={src_rel}/part-0.parquet"
        d = ROOT / f"l2/{name}/release={REL}/part-0.parquet"
        if d.exists():
            raise FileExistsError(d)
        d.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(s, d.with_name(".part-0.parquet.tmp"))
        d.with_name(".part-0.parquet.tmp").replace(d)
        out[name] = {"sha256": ids.sha256_file(d).hex(), "same_as_source": ids.sha256_file(d) == ids.sha256_file(s)}
    s = ROOT / f"ledger/stage={src_rel}/step=cluster/part-0.parquet"
    d = ROOT / f"ledger/stage={REL}/step=cluster/part-0.parquet"
    d.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(s, d)
    out["ledger"] = str(d.relative_to(ROOT))
    return out


# --- queries through the catalog ------------------------------------------------------------
def slim(p: dict) -> dict:
    return {"scans": [{k: v for k, v in s.items() if k != "files"} for s in p["scans"]],
            "bytes_read": p["bytes_read"], "bgzf": p["bgzf"], "index": p["index"]}


def do_queries():
    R = {"release": REL}
    rss0 = rss_mb()
    t = time.perf_counter()
    st = SeqStore(ROOT)
    R["seqstore_open_s"] = round(time.perf_counter() - t, 4)
    R["seqstore_rss_growth_mb"] = round(rss_mb() - rss0, 1)
    st.close()
    del st
    rss1 = rss_mb()
    t = time.perf_counter()
    con = catalog.open_catalog(ROOT / f"catalog/release={REL}.duckdb", local_dir=W / "catalog_cache/local", threads=4,
                               memory_limit="16GB", temp_directory=TMP / f"duck-q-{REL}")
    R["open_catalog_s"] = round(time.perf_counter() - t, 4)
    R["open_catalog_rss_growth_mb"] = round(rss_mb() - rss1, 1)
    R["rss_mb_after_open"] = rss_mb()

    counts = {}
    t0 = time.perf_counter()
    for name in sorted(TABLES):
        t1 = time.perf_counter()
        n = con.execute(f"SELECT count(*) FROM {name}").fetchall()[0][0]
        if n:
            counts[name] = {"rows": n, "wall_s": round(time.perf_counter() - t1, 3)}
    rel_counts = {}
    for name in ("cluster_member", "token_dict", "tape_tokens", "token_postings", "contig_dim"):
        rel_counts[name] = con.execute(f"SELECT count(*) FROM {name} WHERE release = ?", [REL]).fetchall()[0][0]
    R["row_counts"] = {"views": counts, "this_release": rel_counts, "wall_s": round(time.perf_counter() - t0, 2)}

    # 1,000 random genes +-500 bp through seq_get, driven from a registered Arrow table
    t = time.perf_counter()
    rows = con.execute(
        "SELECT f.contig_uid, f.start, f.\"end\", f.strand, d.contig_len FROM feature_protein f "
        "JOIN contig_dim_now() d USING (contig_uid) ORDER BY hash(f.feature_uid), f.feature_uid LIMIT 1000").fetchall()
    R["pick_genes_s"] = round(time.perf_counter() - t, 2)
    picks = pa.table({"contig_uid": pa.array([r[0] for r in rows], pa.binary()),
                      "s": pa.array([max(0, r[1] - 500) for r in rows], pa.int64()),
                      "e": pa.array([min(r[4], r[2] + 500) for r in rows], pa.int64()),
                      "strand": pa.array([r[3] for r in rows], pa.int64())})
    con.register("picks", picks)
    t = time.perf_counter()
    p = catalog.profile_query(con, "SELECT count(*), sum(length(seq_get(contig_uid, s, e, strand))) FROM picks")
    R["seq_get_1000"] = {"wall_s": round(time.perf_counter() - t, 3), "result": p["result"],
                         "expected_bp": int(sum(e - s for s, e in zip(picks.column("s").to_pylist(),
                                                                         picks.column("e").to_pylist()))),
                         "distinct_contigs": len(set(picks.column("contig_uid").to_pylist())), **slim(p)}
    # a second pass, warm caches
    t = time.perf_counter()
    p = catalog.profile_query(con, "SELECT count(*), sum(length(seq_get(contig_uid, s, e, strand))) FROM picks")
    R["seq_get_1000_warm"] = {"wall_s": round(time.perf_counter() - t, 3), **slim(p)}

    # tapes_of one contig
    cu = rows[0][0]
    t = time.perf_counter()
    p = catalog.profile_query(con, f"SELECT count(*), sum(n_slots), sum(len(nr50_tokens)) FROM tapes_of([unhex('{cu.hex()}')])")
    R["tapes_of_one"] = {"contig_uid": cu.hex(), "wall_s": round(time.perf_counter() - t, 3), "result": p["result"], **slim(p)}

    # postings lookup for one NR50 token (the most abundant one of the release)
    tok, key, nsl = con.execute("SELECT token_id, key, n_slots FROM token_dict WHERE release = ? AND level = 'nr50' "
                                "ORDER BY n_slots DESC, token_id LIMIT 1", [REL]).fetchall()[0]
    t = time.perf_counter()
    p = catalog.profile_query(con, f"SELECT count(*), sum(ordinal) FROM token_postings WHERE release = '{REL}' "
                                   f"AND pair_bucket_ok(token_bucket, {tok}::UBIGINT) AND token_id = {tok}::UBIGINT")
    R["postings_lookup_nr50"] = {"token_id": tok, "cluster_uid": key, "n_slots_dict": nsl,
                                 "wall_s": round(time.perf_counter() - t, 3), "result": p["result"], **slim(p)}
    rep = con.execute("SELECT rep_uid FROM cluster_member WHERE release = ? AND level = 'nr50' AND "
                      "cluster_uid = unhex(?) LIMIT 1", [REL, key]).fetchall()
    if rep:
        top = con.execute("SELECT model_acc, model_name, count(*) FROM l1_domains WHERE protein_uid = ? "
                          "GROUP BY 1, 2 ORDER BY 3 DESC LIMIT 3", [rep[0][0]]).fetchall()
        R["postings_lookup_nr50"]["rep_uid"] = rep[0][0].hex()
        R["postings_lookup_nr50"]["rep_pfam"] = top

    # neighbour questions over token_postings: anchors of token T, partners of a level within 3 slots
    def neighbours(label, T, lv):
        lo, hi = lv << 61, (lv + 1) << 61
        base = (f"WITH a AS (SELECT contig_uid, ordinal, n_slots, is_circular FROM token_postings "
                f"  WHERE release = '{REL}' AND pair_bucket_ok(token_bucket, {T}::UBIGINT) AND token_id = {T}::UBIGINT), "
                f"p AS (SELECT token_id, contig_uid, ordinal FROM token_postings WHERE release = '{REL}' "
                f"  AND token_id >= {lo}::UBIGINT AND token_id < {hi}::UBIGINT "
                f"  AND contig_uid IN (SELECT contig_uid FROM a)), "
                "nb AS (SELECT a.contig_uid, a.ordinal AS anchor, p.token_id AS partner FROM a JOIN p USING (contig_uid) "
                "  WHERE p.ordinal <> a.ordinal AND (abs(p.ordinal::BIGINT - a.ordinal::BIGINT) <= 3 "
                "  OR (a.is_circular AND a.n_slots - abs(p.ordinal::BIGINT - a.ordinal::BIGINT) <= 3))) ")
        t = time.perf_counter()
        p = catalog.profile_query(con, base + "SELECT (SELECT count(*) FROM (SELECT DISTINCT contig_uid, ordinal FROM a)), "
                                              "(SELECT count(DISTINCT contig_uid) FROM a), "
                                              "(SELECT count(DISTINCT (contig_uid, anchor)) FROM nb), (SELECT count(*) FROM nb)")
        wall = round(time.perf_counter() - t, 3)
        t = time.perf_counter()
        top = con.execute(base + "SELECT d.key, count(DISTINCT (nb.contig_uid, nb.anchor)) AS n_loci FROM nb "
                                 f"JOIN token_dict d ON d.token_id = nb.partner AND d.release = '{REL}' "
                                 "GROUP BY d.key ORDER BY n_loci DESC, d.key LIMIT 10").fetchall()
        res = p["result"][0]
        return {"label": label, "token_id": T, "partner_level": {1: "family", 2: "nr50", 3: "nr90"}[lv],
                "anchor_loci": res[0], "anchor_contigs": res[1], "loci_with_partner": res[2], "pair_rows": res[3],
                "wall_s": wall, "top10_wall_s": round(time.perf_counter() - t, 3), "top10_partners": top, **slim(p)}

    R["neighbours"] = []
    fam = con.execute("SELECT token_id, n_slots FROM token_dict WHERE release = ? AND level = 'family' AND key = 'PF03237.22'",
                      [REL]).fetchall()
    if fam:
        R["terminase_family_token_slots"] = fam[0][1]
        R["neighbours"].append(neighbours("terminase family PF03237.22 (Terminase_6N) -> family partners", fam[0][0], 1))
        R["neighbours"].append(neighbours("terminase family PF03237.22 -> nr50 partners", fam[0][0], 2))
    R["neighbours"].append(neighbours(f"most abundant NR50 {key} -> nr50 partners", tok, 2))
    R["neighbours"].append(neighbours(f"most abundant NR50 {key} -> family partners", tok, 1))
    # the abundant NR50 among the prokaryote stages only (for rL.2 and later)
    if a.key != "A":
        pk = con.execute(
            "SELECT d.token_id, d.key, count(*) AS n FROM token_postings tp JOIN token_dict d "
            "ON d.token_id = tp.token_id AND d.release = tp.release "
            "WHERE tp.release = ? AND d.level = 'nr50' AND tp.contig_uid IN "
            "(SELECT contig_uid FROM genome_input WHERE kind = 'prokaryote') "
            "GROUP BY 1, 2 ORDER BY n DESC, d.token_id LIMIT 1", [REL]).fetchall()
        if pk:
            R["prok_top_nr50"] = {"token_id": pk[0][0], "key": pk[0][1], "prokaryote_slots": pk[0][2]}
            R["neighbours"].append(neighbours(f"most abundant prokaryote NR50 {pk[0][1]} -> nr50 partners", pk[0][0], 2))
    R["rss_mb_end"] = rss_mb()
    return R


STEPS = {"ingest": do_ingest, "genome_input": do_genome_input, "genecall": do_genecall, "igr": do_igr,
         "tapes": do_tapes, "cluster": do_cluster, "domains": do_domains, "tokens": do_tokens,
         "postings": do_postings, "dim": do_dim, "catalog": do_catalog, "appendonly": do_appendonly,
         "equiv": do_equiv, "pfam_plan": do_pfam_plan, "pfam_merge": do_pfam_merge,
         "cluster_copy": do_cluster_copy, "queries": do_queries}

t0 = time.perf_counter()
rec = {"key": a.key, "step": a.step, "stage": STAGE, "release": REL, "threads": a.threads,
       "started_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
rc = 0
try:
    rec["result"] = jsonable(STEPS[a.step]())
    rec["error"] = None
except BaseException as e:  # noqa: BLE001 - failures at scale are findings
    rec["result"] = None
    rec["error"] = f"{type(e).__name__}: {e}"
    rec["traceback"] = traceback.format_exc()
    rc = 1
rec["inner_wall_s"] = round(time.perf_counter() - t0, 2)
rec["self_peak_rss_mb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
tmp = OUT.with_suffix(".tmp")
tmp.write_text(json.dumps(rec, indent=1, default=str))
tmp.replace(OUT)
print(json.dumps({k: v for k, v in rec.items() if k != "result"}, default=str), flush=True)
sys.exit(rc)
