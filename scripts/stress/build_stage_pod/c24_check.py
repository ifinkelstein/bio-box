"""Card C-24: the checks after stage H1 was built into the data root. In the image (pyarrow, duckdb, numpy; the tree's
src on PYTHONPATH), read-only on the root:

    python c24_check.py --root ROOT --before DIR --after DIR --log UNIT.LOG --objects OBJECTS.JSON --out DIR
        --tmp DIR [--stage H1] [--release rM.3] [--earlier S1,P1,P2,P3,P4,P5,P6] [--view DIR]

Every table is read by the explicit file list of the run record that published it (AGENTS §3.9; no glob over a table).

1. Before against after (`listing.tsv`: path, bytes, mtime of every file of the root; `s1.sha256`: sha256 of every
   file under a `stage=S1` directory): no file that existed before is missing or has another size or mtime, no S1 file
   another sha256; every new file has a `stage=<stage>` path component, or is the release's definition
   (`release_def/release=<release>/`); no file of the release (`release=<release>`, `stage=<release>`) but that.
2. `root_added.tsv`: every new file with its bytes and sha256 (for the durable copy).
3. The stage's steps: each run record's CPU-minutes, wall seconds, peak RSS and threads, against the code's estimate
   where it makes one (ingest, the scan shards, the two CRISPR steps; `c24_est.py rates`), and the refusal line of
   the build's log.
4. The data: genomes and contigs of `genome_input`; the contigs another stage stores (against `contig_index` of the
   earlier stages) and the gene-call ledger's skipped contigs (`owned_by:<stage>`), which must be the same set;
   genes (`feature_protein` rows), distinct proteins (`protein_seq` rows) and those in no `protein_seq` of the earlier
   stages (DuckDB anti-join under a memory limit, spilling to --tmp); per gene the genetic code the caller chose;
   the Pfam scan (plan: proteins and residues scanned; `l1_domains`: rows, distinct proteins with a domain, models);
   the CRISPR arrays of each producer (feature rows, `l1_crispr` rows, the ledger's counts); the genetic-code flag:
   per genome the mean gene length, the stage's median gene length, and the share of genomes whose mean is under
   half of it (over the genomes whose contig the stage owns, and over all of them with the earlier stage's calls of
   the shared contigs). With --view (the card's E-16 view): the default producer's three files there against the
   root's, by sha256.

Writes into --out: report.txt, summary.json, root_added.tsv, steps.tsv, genomes.tsv. Exit 0 when every check of 1 and
4's equalities pass, 1 otherwise.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from biobox.l1 import crispr
from biobox.l1.genecall import GENE_CALL_VERSION as GCV

LIB = "pfam-38.2"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def listing(d: Path) -> dict[str, tuple[int, str]]:
    out = {}
    for line in (d / "listing.tsv").read_text().splitlines():
        path, size, mtime = line.split("\t")
        out[path] = (int(size), mtime)
    return out


def sums(d: Path) -> dict[str, str]:
    out = {}
    for line in (d / "s1.sha256").read_text().splitlines():
        h, path = line.split("  ", 1)
        out[path] = h
    return out


def record(root: Path, stage: str, step: str) -> dict | None:
    p = root / "runs" / f"stage={stage}" / f"{step}.json"
    return json.loads(p.read_text()) if p.exists() else None


def files(root: Path, stage: str, step: str, prefix: str) -> list[Path]:
    """The outputs of the step's run record under `prefix` (an explicit list)."""
    r = record(root, stage, step)
    if r is None:
        return []
    return [root / o["path"] for o in r["outputs"] if o["path"].startswith(prefix) and o["path"].endswith(".parquet")]


def read(paths: list[Path], cols: list[str]) -> pa.Table:
    parts = [pq.read_table(p, columns=cols, partitioning=None) for p in paths]
    return pa.concat_tables(parts) if parts else pa.table({c: [] for c in cols})


def hexs(v) -> str:
    return v.hex() if isinstance(v, (bytes, bytearray)) else str(v)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--before", required=True, type=Path)
    ap.add_argument("--after", required=True, type=Path)
    ap.add_argument("--log", required=True, type=Path)
    ap.add_argument("--objects", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--tmp", required=True, type=Path)
    ap.add_argument("--stage", default="H1")
    ap.add_argument("--release", default="rM.3")
    ap.add_argument("--earlier", default="S1,P1,P2,P3,P4,P5,P6")
    ap.add_argument("--view", type=Path, default=None)
    ap.add_argument("--est", type=Path, default=Path(__file__).with_name("c24_est.py"))
    a = ap.parse_args()
    root, out, st, rel = a.root, a.out, a.stage, a.release
    earlier = a.earlier.split(",")
    out.mkdir(parents=True, exist_ok=True)
    a.tmp.mkdir(parents=True, exist_ok=True)
    report: list[str] = []
    failed: list[str] = []
    summary: dict = {}

    def say(s: str) -> None:
        report.append(s)
        print(s, flush=True)

    # 1. before against after
    b, f = listing(a.before), listing(a.after)
    gone = sorted(p for p in b if p not in f)
    changed = sorted(p for p in b if p in f and f[p] != b[p])
    new = sorted(p for p in f if p not in b)
    mine = [p for p in new if f"stage={st}" in p.split("/")]
    defn = [p for p in new if p.startswith(f"release_def/release={rel}/")]
    unexpected = sorted(set(new) - set(mine) - set(defn))
    release_files = sorted(p for p in f if f"release={rel}" in p.split("/") or f"stage={rel}" in p.split("/")
                           or p in (f"manifests/release={rel}.json", f"catalog/release={rel}.duckdb"))
    release_steps = sorted(set(release_files) - set(defn))
    sb, sa = sums(a.before), sums(a.after)
    s1_changed = sorted(p for p in sb if sa.get(p) != sb[p])
    s1_new = sorted(p for p in sa if p not in sb)
    say(f"before: {len(b)} files; after: {len(f)}; new {len(new)} ({len(mine)} under stage={st}, {len(defn)} of "
        f"{rel}'s definition); gone {len(gone)}; size or mtime changed {len(changed)}; S1 files hashed before {len(sb)}, "
        f"after {len(sa)}, sha256 changed {len(s1_changed)}, new {len(s1_new)}")
    for what, lst in (("gone", gone), ("changed", changed), ("S1 sha256 changed", s1_changed), ("S1 new", s1_new),
                      ("unexpected new", unexpected), (f"files of a release step of {rel}", release_steps)):
        if lst:
            failed.append(f"{what}: {len(lst)}")
            for p in lst[:50]:
                say(f"  {what}: {p}")
    summary["files"] = {"before": len(b), "after": len(f), "new": len(new), "new_under_stage": len(mine),
                        "release_def": defn, "gone": gone, "changed": changed, "s1_sha256_changed": s1_changed,
                        "unexpected": unexpected, "release_step_files": release_steps}

    # 2. root_added.tsv
    added = [(p, (root / p).stat().st_size, sha256(root / p)) for p in new]
    with open(out / "root_added.tsv", "w") as fh:
        fh.write("path\tbytes\tsha256\n")
        for p, n, h in added:
            fh.write(f"{p}\t{n}\t{h}\n")
    summary["added"] = {"files": len(added), "bytes": sum(n for _, n, _ in added)}
    say(f"root_added.tsv: {len(added)} files, {sum(n for _, n, _ in added):,} bytes")
    for p, n, _ in added:
        if f[p][0] != n:
            failed.append(f"{p}: {n} bytes now, {f[p][0]} in the after listing")

    # 3. steps
    rows = []
    for p in sorted((root / "runs" / f"stage={st}").glob("*.json")):
        if p.name.endswith(".publishing.json"):
            failed.append(f"intent file {p.name}")
            continue
        r = json.loads(p.read_text())
        rows.append({"step": p.name[:-5], "cpu_min": round(r["cpu_min"], 3), "wall_s": round(r["wall_s"], 1),
                     "peak_rss_mib": round(r["peak_rss_bytes"] / 2**20), "threads": r["threads"],
                     "started_at": r["started_at"], "over_cpu_budget": r["over_cpu_budget"],
                     "image_digest": r["image_digest"], "code_id": r.get("code_id"),
                     "code_commit": r.get("code_commit"), "outputs": len(r["outputs"])})
    est = subprocess.run([sys.executable, str(a.est), "rates", "--root", str(root), "--objects", str(a.objects),
                          "--stage", st], capture_output=True, text=True)
    rates = json.loads(est.stdout) if est.stdout.strip() else {}
    by = {x["step"]: x for x in rates.get("rows", [])}
    for r in rows:
        if r["step"] in by:
            r["estimate_cpu_min"], r["ratio"] = by[r["step"]]["estimate"], by[r["step"]]["ratio"]
    if rates.get("outside_20pct"):
        failed.append(f"rates outside 20%: {rates['outside_20pct']}")
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(out / "steps.tsv", "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")
    shards = [r for r in rows if ".shard-" in r["step"]]
    summary["steps"] = {"records": len(rows), "names": [r["step"] for r in rows if ".shard-" not in r["step"]],
                        "shards": len(shards),
                        "shard_cpu_min": round(sum(r["cpu_min"] for r in shards), 2),
                        "shard_cpu_min_max": max((r["cpu_min"] for r in shards), default=None),
                        "shard_ratio_min": min((r.get("ratio", 0) for r in shards), default=None),
                        "shard_ratio_max": max((r.get("ratio", 0) for r in shards), default=None),
                        "cpu_min_total": round(sum(r["cpu_min"] for r in rows), 2),
                        "cpu_min_max": max((r["cpu_min"] for r in rows), default=None),
                        "over_budget": [r["step"] for r in rows if r["over_cpu_budget"]],
                        "images": sorted({r["image_digest"] for r in rows}),
                        "code_ids": sorted({str(r["code_id"]) for r in rows}),
                        "rates": rates}
    for r in rows:
        if ".shard-" not in r["step"]:
            say(f"step {r['step']}: {r['cpu_min']} CPU-min, {r['wall_s']} s, {r['peak_rss_mib']} MiB, "
                f"{r['threads']} threads" + (f", estimate {r['estimate_cpu_min']} (ratio {r['ratio']})"
                                             if "ratio" in r else ""))
    say(f"scan shards: {len(shards)}, {summary['steps']['shard_cpu_min']} CPU-min, largest "
        f"{summary['steps']['shard_cpu_min_max']}, ratio to estimate {summary['steps']['shard_ratio_min']} to "
        f"{summary['steps']['shard_ratio_max']}")
    log = a.log.read_text(errors="replace").splitlines() if a.log.exists() else []
    refusal = [x for x in log if "build-stage: refused" in x]
    gc_line = [x for x in log if '"step": "genecall"' in x]
    exit_line = [x for x in log if x.startswith("pod_run: ") and " exited " in x]
    usage = [x for x in log if x.startswith("pod_run: ") and " usage " in x]
    summary["log"] = {"refusal": refusal, "exit": exit_line, "usage": usage}
    say(f"log: {refusal[-1] if refusal else 'no refusal line'}; {exit_line[-1] if exit_line else 'no exit line'}")

    # 4. data
    gi = read([root / f"l0seq/genome_input/stage={st}/part-0.parquet"], ["contig_uid", "genome_uid", "kind",
                                                                          "genetic_code", "is_circular"])
    contigs = set(gi.column("contig_uid").to_pylist())
    genome_of = dict(zip(gi.column("contig_uid").to_pylist(), gi.column("genome_uid").to_pylist()))
    own_idx = files(root, st, "ingest", f"l0seq/contig_index/stage={st}/")
    ci = read(own_idx, ["contig_uid", "contig_len", "is_default_variant"])
    owned = set(ci.column("contig_uid").to_pylist())
    lens = dict(zip(ci.column("contig_uid").to_pylist(), ci.column("contig_len").to_pylist()))
    stored_by: dict[bytes, str] = {}
    for s in earlier:
        t = read(files(root, s, "ingest", f"l0seq/contig_index/stage={s}/"), ["contig_uid", "contig_len"])
        for cu, n in zip(t.column("contig_uid").to_pylist(), t.column("contig_len").to_pylist()):
            if cu in contigs and cu not in stored_by:
                stored_by[cu] = s
                lens.setdefault(cu, n)
    bases_owned = int(pc.sum(pc.filter(ci.column("contig_len"), ci.column("is_default_variant"))).as_py() or 0)
    led = read([root / f"ledger/stage={st}/step=genecall/part-0.parquet"], ["object_key", "reason_code", "n_in",
                                                                             "n_kept"])
    skipped = {hexs(k): why for k, why in zip(led.column("object_key").to_pylist(), led.column("reason_code").to_pylist())
               if why.startswith("owned_by:")}
    want = {cu.hex(): f"owned_by:{s}" for cu, s in stored_by.items()}
    same = skipped == want
    if not same:
        failed.append("the gene-call ledger's skipped contigs are not the contigs the earlier stages store")
    if owned & set(stored_by):
        failed.append("a contig of the stage's contig_index is also stored by an earlier stage")
    if contigs != owned | set(stored_by):
        failed.append("genome_input's contigs are not those of contig_index and the earlier stages")
    gc_res = {}
    if gc_line:
        try:
            gc_res = json.loads(gc_line[-1][gc_line[-1].index("{"):])
        except ValueError:
            gc_res = {}
    res_skipped = {e.get("contig_uid"): e.get("owner_stage") for e in
                   (gc_res.get("result", gc_res).get("skipped_contigs", []) if isinstance(gc_res, dict) else [])}
    summary["contigs"] = {"genome_input_rows": gi.num_rows, "genomes": len(set(genome_of.values())),
                          "contigs": len(contigs), "owned": len(owned), "bases_owned": bases_owned,
                          "stored_by_earlier": dict(Counter(stored_by.values())),
                          "ledger_skipped": len(skipped), "ledger_equals_stored_by": same,
                          "result_skipped": len(res_skipped),
                          "result_equals_stored_by": res_skipped == {k: v[9:] for k, v in want.items()},
                          "kinds": sorted(set(gi.column("kind").to_pylist())),
                          "codes": sorted(set(gi.column("genetic_code").to_pylist())),
                          "circular": int(pc.sum(gi.column("is_circular").cast(pa.int64())).as_py() or 0)}
    say(f"contigs: {summary['contigs']}")

    fp_paths = files(root, st, "genecall", f"l1/feature_protein/gene_call_version={GCV}/stage={st}/")
    ps_paths = files(root, st, "genecall", f"l0seq/protein_seq/gene_call_version={GCV}/stage={st}/")
    fp = read(fp_paths, ["contig_uid", "start", "end", "protein_uid", "genetic_code"])
    n_prot = sum(pq.read_metadata(p).num_rows for p in ps_paths)
    earlier_ps = [p for s in earlier for p in files(root, s, "genecall",
                                                     f"l0seq/protein_seq/gene_call_version={GCV}/stage={s}/")]
    con = duckdb.connect(config={"threads": 2, "memory_limit": "8GB", "temp_directory": str(a.tmp / "duckdb"),
                                 "preserve_insertion_order": False})
    lst = lambda ps: "[" + ", ".join("'" + str(p).replace("'", "''") + "'" for p in ps) + "]"  # noqa: E731
    new_prot, = con.execute(
        f"SELECT count(*) FROM (SELECT protein_uid FROM read_parquet({lst(ps_paths)}, hive_partitioning = false) "
        f"ANTI JOIN (SELECT protein_uid FROM read_parquet({lst(earlier_ps)}, hive_partitioning = false)) "
        f"USING (protein_uid))").fetchone()
    res_all, len_new = con.execute(
        f"SELECT sum(length), (SELECT sum(length) FROM read_parquet({lst(ps_paths)}, hive_partitioning = false) "
        f"ANTI JOIN (SELECT protein_uid FROM read_parquet({lst(earlier_ps)}, hive_partitioning = false)) "
        f"USING (protein_uid)) FROM read_parquet({lst(ps_paths)}, hive_partitioning = false)").fetchone()
    con.close()
    codes = Counter(fp.column("genetic_code").to_pylist())
    summary["genes"] = {"genes": fp.num_rows, "distinct_proteins": n_prot, "residues": int(res_all or 0),
                        "new_distinct_proteins": int(new_prot), "new_residues": int(len_new or 0),
                        "in_earlier_stages": n_prot - int(new_prot), "earlier_protein_seq_files": len(earlier_ps),
                        "genetic_code_per_gene": {str(k): v for k, v in sorted(codes.items())},
                        "result_line": gc_res}
    say(f"genes {fp.num_rows:,}; distinct proteins {n_prot:,}; new against protein_seq of {','.join(earlier)}: "
        f"{int(new_prot):,}; genetic code per gene {dict(codes)}")

    # the genetic-code flag
    glen = pc.subtract(fp.column("end"), fp.column("start")).to_numpy().astype(np.int64)
    gcon = fp.column("contig_uid").to_pylist()
    gcode = fp.column("genetic_code").to_pylist()
    shared_rows = []
    for s in sorted(set(stored_by.values())):
        t = read(files(root, s, "genecall", f"l1/feature_protein/gene_call_version={GCV}/stage={s}/"),
                 ["contig_uid", "start", "end", "genetic_code"])
        t = t.filter(pc.is_in(t.column("contig_uid"), value_set=pa.array(
            sorted(cu for cu, o in stored_by.items() if o == s), pa.binary(16))))
        shared_rows.append(t)
    sh = pa.concat_tables(shared_rows) if shared_rows else None

    def per_genome(cons, ls, cds) -> dict:
        acc: dict = {}
        for cu, n, c in zip(cons, ls, cds):
            g = acc.setdefault(genome_of[cu], [0, 0, 0])
            g[0] += 1
            g[1] += int(n)
            g[2] += int(c != 11)
        return acc

    own_g = per_genome(gcon, glen.tolist(), gcode)
    all_g = dict(own_g)
    if sh is not None:
        sl = pc.subtract(sh.column("end"), sh.column("start")).to_numpy().astype(np.int64)
        for k, v in per_genome(sh.column("contig_uid").to_pylist(), sl.tolist(),
                               sh.column("genetic_code").to_pylist()).items():
            all_g[k] = v
        all_len = np.concatenate([glen, sl])
    else:
        all_len = glen
    flag = {}
    for name, gd, ls in (("owned", own_g, glen), ("all", all_g, all_len)):
        med = float(np.median(ls))
        means = np.array([v[1] / v[0] for v in gd.values()])
        med_of_means = float(np.median(means))
        flag[name] = {"genomes": len(gd), "genes": int(len(ls)), "median_gene_bp": med,
                      "median_of_genome_means_bp": med_of_means,
                      "genomes_mean_under_half_median": int((means < med / 2).sum()),
                      "share_under_half_median": round(float((means < med / 2).mean()), 6),
                      "genomes_mean_under_half_median_of_means": int((means < med_of_means / 2).sum()),
                      "genomes_with_any_gene_not_code_11": sum(1 for v in gd.values() if v[2]),
                      "genomes_all_genes_not_code_11": sum(1 for v in gd.values() if v[2] == v[0])}
    no_genes = sorted(g.hex() for g in set(genome_of.values()) if g not in all_g)
    flag["genomes_without_genes"] = len(no_genes)
    summary["genetic_code_flag"] = flag
    say(f"genetic-code flag: {flag}")
    with open(out / "genomes.tsv", "w") as fh:
        fh.write("genome_uid\tcontig_len\tgenes\tmean_gene_bp\tgenes_not_code_11\tcalls_of\n")
        inv = {g: cu for cu, g in genome_of.items()}
        for g in sorted(all_g):
            n, tot, nc = all_g[g]
            cu = inv[g]
            fh.write(f"{g.hex()}\t{lens.get(cu)}\t{n}\t{tot / n:.1f}\t{nc}\t{stored_by.get(cu, st)}\n")

    # Pfam
    plan_p = root / f"l1/domain_scan/library={LIB}/stage={st}/plan.json"
    plan = json.loads(plan_p.read_text()) if plan_p.exists() else {}
    dom = read(files(root, st, f"domains-{LIB}", f"l1/domains/library={LIB}/stage={st}/"), ["protein_uid", "model_acc"])
    summary["pfam"] = {"proteins_scanned": plan.get("proteins"), "residues_scanned": plan.get("residues"),
                       "shards": len(plan.get("shards", [])), "shard_residues": plan.get("shard_residues"),
                       "domain_rows": dom.num_rows, "proteins_with_domain": pc.count_distinct(dom.column("protein_uid")).as_py(),
                       "models_with_hit": pc.count_distinct(dom.column("model_acc")).as_py()}
    if plan.get("proteins") is not None and plan["proteins"] != n_prot:
        failed.append(f"the scan plan holds {plan['proteins']} proteins, protein_seq {n_prot}")
    say(f"pfam: {summary['pfam']}")

    # CRISPR
    cr = {}
    for kind, prod in (("default", crispr.producer()), ("near_cas", crispr.near_cas_producer("38.2"))):
        step = f"crispr-{prod}"
        feats = files(root, st, step, f"l1/feature/producer={prod}/stage={st}/")
        elems = files(root, st, step, f"l1/crispr/producer={prod}/stage={st}/")
        ledp = files(root, st, step, f"ledger/stage={st}/step={step}/")
        lt = read(ledp, ["n_in", "n_kept", "n_rejected", "reason_code"]).to_pydict()
        c = {"producer": prod, "arrays": sum(pq.read_metadata(p).num_rows for p in feats),
             "elements": sum(pq.read_metadata(p).num_rows for p in elems), "contigs_searched": 0,
             "contigs_skipped_other_stage": 0, "calls_in": 0, "calls_kept": 0, "rejected": {}}
        for n_in, n_kept, n_rej, why in zip(lt["n_in"], lt["n_kept"], lt["n_rejected"], lt["reason_code"]):
            if why == "ok":
                c["contigs_searched"] += 1
                c["calls_in"] += n_in
                c["calls_kept"] += n_kept
            elif why.startswith("owned_by:"):
                c["contigs_skipped_other_stage"] += 1
            else:
                c["rejected"][why] = c["rejected"].get(why, 0) + n_rej
        if c["calls_kept"] != c["arrays"]:
            failed.append(f"{kind}: the ledger keeps {c['calls_kept']}, the feature files hold {c['arrays']}")
        c["arrays_per_gbp"] = round(c["arrays"] / (bases_owned / 1e9), 3) if bases_owned else None
        cr[kind] = c
        if a.view is not None and kind == "default":
            same_files = []
            for p in feats + elems + ledp:
                relp = p.relative_to(root)
                v = a.view / relp
                same_files.append((str(relp), v.exists() and sha256(v) == sha256(p)))
            c["view_byte_identical"] = all(x for _, x in same_files) and len(same_files) == 3
            if not c["view_byte_identical"]:
                failed.append(f"default producer: not byte-identical to the view {a.view}: {same_files}")
    summary["crispr"] = cr
    say(f"crispr: {cr}")

    summary["failed"] = failed
    (out / "summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True, default=str))
    say("RESULT: " + ("all checks pass" if not failed else "FAILED: " + "; ".join(failed)))
    (out / "report.txt").write_text("\n".join(report) + "\n")
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
