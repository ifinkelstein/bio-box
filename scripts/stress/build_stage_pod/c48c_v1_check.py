"""Card C-48c part 1: post-build checks for stage V1 (UHGV v1.0 vOTU representatives), in the form of card C-47's
c47_check.py, plus the counts card C-39's note fixes or estimates (`c48c_v1_est.expected_ranges`).

    python c48c_v1_check.py --root ROOT --before DIR --after DIR --log UNIT.LOG --objects OBJECTS.JSON
        --out DIR --tmp DIR [--stage V1] [--earlier S1,P1,P2,P3,P4,P5,P6,H1,PL1,C5] [--before-est FILE]
        [--no-expect]

Every table is read by the explicit file list of the run record that published it (AGENTS §3.9). --no-expect (the
rehearsal's sample) skips the comparison with C-39's counts. Writes into --out: report.txt, summary.json,
root_added.tsv, steps.tsv.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from biobox.l1 import crispr, genecall
from biobox.l1.genecall import GENE_CALL_VERSION as GCV

LIB = "pfam-38.2"
COLS = (
    "protein_uid", "model_acc", "model_name", "model_len", "seq_bitscore", "seq_evalue_zref",
    "dom_bitscore", "dom_ievalue_zref", "env_start", "env_end", "ali_start", "ali_end",
    "hmm_from", "hmm_to", "cutoff", "z", "dom_z",
)


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
    p = d / "earlier.sha256"
    for line in (p.read_text().splitlines() if p.exists() else []):
        h, path = line.split("  ", 1)
        out[path] = h
    return out


def record(root: Path, stage: str, step: str) -> dict | None:
    p = root / "runs" / f"stage={stage}" / f"{step}.json"
    return json.loads(p.read_text()) if p.exists() else None


def files(root: Path, stage: str, step: str, prefix: str) -> list[Path]:
    r = record(root, stage, step)
    if r is None:
        return []
    return [root / o["path"] for o in r["outputs"] if o["path"].startswith(prefix) and o["path"].endswith(".parquet")]


def read(paths: list[Path], cols: list[str]) -> pa.Table:
    parts = [pq.read_table(p, columns=cols, partitioning=None) for p in paths]
    return pa.concat_tables(parts) if parts else pa.table({c: [] for c in cols})


def hexs(v) -> str:
    return v.hex() if isinstance(v, (bytes, bytearray)) else str(v)


def lst(ps) -> str:
    return "[" + ", ".join("'" + str(p).replace("'", "''") + "'" for p in ps) + "]"


def rows_of(ps) -> int:
    return sum(pq.read_metadata(p).num_rows for p in ps)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--before", required=True, type=Path)
    ap.add_argument("--after", required=True, type=Path)
    ap.add_argument("--log", required=True, type=Path)
    ap.add_argument("--objects", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--tmp", required=True, type=Path)
    ap.add_argument("--stage", default="V1")
    ap.add_argument("--earlier", default="S1,P1,P2,P3,P4,P5,P6,H1,PL1,C5")
    ap.add_argument("--before-est", type=Path, default=None)
    ap.add_argument("--est", type=Path, default=Path(__file__).with_name("c48c_v1_est.py"))
    ap.add_argument("--no-expect", action="store_true")
    a = ap.parse_args()
    spec = importlib.util.spec_from_file_location("c48c_v1_est", a.est)
    est_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(est_mod)
    root, out, st = a.root, a.out, a.stage
    earlier = [s for s in a.earlier.split(",") if s]
    out.mkdir(parents=True, exist_ok=True)
    a.tmp.mkdir(parents=True, exist_ok=True)
    report: list[str] = []
    failed: list[str] = []
    summary: dict = {}

    def say(s: str) -> None:
        report.append(s)
        print(s, flush=True)

    # 1. before vs after
    b, f = listing(a.before), listing(a.after)
    gone = sorted(p for p in b if p not in f)
    changed = sorted(p for p in b if p in f and f[p] != b[p])
    new = sorted(p for p in f if p not in b)
    mine = [p for p in new if f"stage={st}" in p.split("/")]
    unexpected = sorted(set(new) - set(mine))
    sb, sa = sums(a.before), sums(a.after)
    e_changed = sorted(p for p in sb if sa.get(p) != sb[p])
    e_new = sorted(p for p in sa if p not in sb)
    stop_left = [p for p in f if p == "STOP"]
    say(
        f"before: {len(b)} files; after: {len(f)}; new {len(new)} ({len(mine)} under stage={st}); "
        f"gone {len(gone)}; size or mtime changed {len(changed)}; earlier hashed before {len(sb)}, "
        f"after {len(sa)}, sha256 changed {len(e_changed)}, new {len(e_new)}"
    )
    for what, items in (
        ("gone", gone),
        ("changed", changed),
        ("earlier sha256 changed", e_changed),
        ("earlier new", e_new),
        ("unexpected new", unexpected),
        ("STOP left in root", stop_left),
    ):
        if items:
            failed.append(f"{what}: {len(items)}")
            for p in items[:50]:
                say(f"  {what}: {p}")
    summary["files"] = {
        "before": len(b), "after": len(f), "new": len(new), "new_under_stage": len(mine),
        "gone": gone, "changed": changed, "earlier_hashed": len(sb),
        "earlier_sha256_changed": e_changed, "unexpected": unexpected,
    }

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
            failed.append(f"{p}: {n} bytes now, {f[p][0]} in after listing")

    # 3. steps
    rows = []
    for p in sorted((root / "runs" / f"stage={st}").glob("*.json")):
        if p.name.endswith(".publishing.json"):
            failed.append(f"intent file {p.name}")
            continue
        r = json.loads(p.read_text())
        rows.append({
            "step": p.name[:-5],
            "cpu_min": round(r["cpu_min"], 3),
            "wall_s": round(r["wall_s"], 1),
            "peak_rss_mib": round((r["peak_rss_bytes"] or 0) / 2**20),
            "threads": r["threads"],
            "started_at": r["started_at"],
            "over_cpu_budget": r["over_cpu_budget"],
            "image_digest": r["image_digest"],
            "code_id": r.get("code_id"),
            "code_commit": r.get("code_commit"),
            "outputs": len(r["outputs"]),
        })
    cmd = [sys.executable, str(a.est), "rates", "--root", str(root), "--objects", str(a.objects), "--stage", st]
    if a.before_est is not None:
        cmd.extend(["--before", str(a.before_est)])
    est = subprocess.run(cmd, capture_output=True, text=True)
    rates = json.loads(est.stdout) if est.stdout.strip() else {"error": est.stderr[-2000:]}
    by = {x["step"]: x for x in rates.get("rows", [])}
    for r in rows:
        if r["step"] in by:
            r["estimate_cpu_min"] = by[r["step"]]["estimate"]
            r["ratio"] = by[r["step"]]["ratio"]
            r["judged"] = by[r["step"]]["judged"]
    if rates.get("over_20pct"):
        failed.append(f"rates over +20% (DEC-53): {rates['over_20pct']}")
    if any(r["over_cpu_budget"] for r in rows):
        failed.append(f"over_cpu_budget: {[r['step'] for r in rows if r['over_cpu_budget']]}")
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(out / "steps.tsv", "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")

    def group(tag: str) -> dict:
        g = [r for r in rows if tag in r["step"]]
        return {
            "n": len(g),
            "cpu_min": round(sum(r["cpu_min"] for r in g), 2),
            "cpu_min_max": max((r["cpu_min"] for r in g), default=None),
            "ratio_min": min((r.get("ratio", 0) for r in g), default=None),
            "ratio_max": max((r.get("ratio", 0) for r in g), default=None),
            "wall_s_max": max((r["wall_s"] for r in g), default=None),
            "peak_rss_mib_max": max((r["peak_rss_mib"] for r in g), default=None),
        }

    summary["steps"] = {
        "records": len(rows),
        "names": [r["step"] for r in rows if ".shard-" not in r["step"]],
        "shards": group(".shard-"),
        "batches": group(".batch-"),
        "cpu_min_total": round(sum(r["cpu_min"] for r in rows), 2),
        "cpu_min_max": max((r["cpu_min"] for r in rows), default=None),
        "largest": max(rows, key=lambda r: r["cpu_min"])["step"] if rows else None,
        "over_budget": [r["step"] for r in rows if r["over_cpu_budget"]],
        "images": sorted({r["image_digest"] for r in rows}),
        "code_ids": sorted({str(r["code_id"]) for r in rows}),
        "rates": rates,
    }
    for r in rows:
        say(
            f"step {r['step']}: {r['cpu_min']} CPU-min, {r['wall_s']} s, {r['peak_rss_mib']} MiB, {r['threads']} threads"
            + (f", estimate {r['estimate_cpu_min']} (ratio {r['ratio']}"
               f"{'' if r['judged'] else ', not judged: DEC-43'})" if "ratio" in r else "")
        )
    say(
        f"batches {summary['steps']['batches']}; shards {summary['steps']['shards']}; total "
        f"{summary['steps']['cpu_min_total']} CPU-min, largest {summary['steps']['largest']} "
        f"{summary['steps']['cpu_min_max']}"
    )

    # 4. contigs & genes
    gi_rec = files(root, st, "inputs", f"l0seq/genome_input/stage={st}/")
    gi = read(gi_rec, ["contig_uid", "genome_uid", "kind", "genetic_code", "is_circular", "assembly_quality"])
    contigs = set(gi.column("contig_uid").to_pylist())
    own_idx = files(root, st, "ingest", f"l0seq/contig_index/stage={st}/")
    ci = read(own_idx, ["contig_uid", "contig_len", "is_default_variant"])
    owned = set(ci.column("contig_uid").to_pylist())
    stored_by: dict[bytes, str] = {}
    for s in earlier:
        t = read(files(root, s, "ingest", f"l0seq/contig_index/stage={s}/"), ["contig_uid", "is_default_variant"])
        for cu, dv in zip(t.column("contig_uid").to_pylist(), t.column("is_default_variant").to_pylist()):
            if dv and cu in contigs and cu not in stored_by:
                stored_by[cu] = s
    bases_owned = int(pc.sum(pc.filter(ci.column("contig_len"), ci.column("is_default_variant"))).as_py() or 0)
    led = read(files(root, st, "genecall", f"ledger/stage={st}/step=genecall/"), ["object_key", "reason_code"])
    skipped_earlier = {
        hexs(k): why for k, why in zip(led.column("object_key").to_pylist(), led.column("reason_code").to_pylist())
        if why.startswith("owned_by:") and why != f"owned_by:{st}"
    }
    skipped_within = sum(
        1 for why in led.column("reason_code").to_pylist()
        if why == f"owned_by:{st}" or why.startswith("duplicate_contig")
    )
    want = {cu.hex(): f"owned_by:{s}" for cu, s in stored_by.items()}
    if skipped_earlier != want:
        failed.append(f"genecall ledger earlier skips ({len(skipped_earlier)}) != stored_by ({len(want)})")
    if owned & set(stored_by):
        failed.append("contig in stage contig_index is also stored by earlier stage")
    if contigs != owned | set(stored_by):
        failed.append("genome_input contigs != contig_index | stored_by")
    occ = read(files(root, st, "ingest", f"l0seq/contig_occurrence/stage={st}/"), ["contig_uid"])
    occ_cu = occ.column("contig_uid").to_pylist()
    occ_shared = sum(1 for cu in occ_cu if cu in stored_by)
    lab = files(root, st, "ingest", f"l0seq/contig_label/stage={st}/")
    summary["contigs"] = {
        "genome_input_rows": gi.num_rows,
        "genomes": len(set(gi.column("genome_uid").to_pylist())),
        "contigs": len(contigs),
        "stored": len(owned),
        "bases_stored": bases_owned,
        "stored_by_earlier": len(stored_by),
        "stored_by_earlier_by_stage": dict(sorted(Counter(stored_by.values()).items())),
        "ledger_skipped_earlier": len(skipped_earlier),
        "ledger_skipped_within": skipped_within,
        "ledger_equals_stored_by": skipped_earlier == want,
        "contig_occurrence_rows": occ.num_rows,
        "contig_occurrence_contigs": len(set(occ_cu)),
        "contig_occurrence_rows_on_contigs_of_earlier_stages": occ_shared,
        "contig_label_rows": rows_of(lab),
        "kinds": sorted(set(gi.column("kind").to_pylist())),
        "codes": sorted(set(gi.column("genetic_code").to_pylist())),
        "quality": dict(Counter(gi.column("assembly_quality").to_pylist())),
        "circular": int(pc.sum(gi.column("is_circular").cast(pa.int64())).as_py() or 0),
    }
    say(f"contigs: {summary['contigs']}")

    fp_paths = files(root, st, "genecall", f"l1/feature_protein/gene_call_version={GCV}/stage={st}/")
    ps_paths = files(root, st, "genecall", f"l0seq/protein_seq/gene_call_version={GCV}/stage={st}/")
    n_genes, n_prot = rows_of(fp_paths), rows_of(ps_paths)
    earlier_ps = [
        p for s in earlier for p in files(root, s, "genecall", f"l0seq/protein_seq/gene_call_version={GCV}/stage={s}/")
    ]
    con = duckdb.connect(config={
        "threads": 2, "memory_limit": "8GB", "temp_directory": str(a.tmp / "duckdb"),
        "preserve_insertion_order": False,
    })
    rp = lambda ps: f"read_parquet({lst(ps)}, hive_partitioning = false)"  # noqa: E731
    new_prot, new_res = con.execute(
        f"SELECT count(*), sum(length) FROM (SELECT protein_uid, length FROM {rp(ps_paths)} "
        f"ANTI JOIN (SELECT protein_uid FROM {rp(earlier_ps)}) USING (protein_uid))"
    ).fetchone()
    res_all, = con.execute(f"SELECT sum(length) FROM {rp(ps_paths)}").fetchone()
    summary["genes"] = {
        "genes": n_genes,
        "distinct_proteins": n_prot,
        "residues": int(res_all or 0),
        "new_distinct_proteins": int(new_prot),
        "new_residues": int(new_res or 0),
        "in_earlier_stages": n_prot - int(new_prot),
        "earlier_protein_seq_files": len(earlier_ps),
    }
    say(
        f"genes {n_genes:,}; distinct proteins {n_prot:,} ({int(res_all or 0):,} residues); "
        f"new against {','.join(earlier)}: {int(new_prot):,} ({int(new_res or 0):,} residues)"
    )

    # 5. Pfam + independent B-47b copy verification (including C-38h / DEC-54(8) shared-contig proteins)
    plan_p = root / f"l1/domain_scan/library={LIB}/stage={st}/plan.json"
    plan = json.loads(plan_p.read_text()) if plan_p.exists() else {}
    dom_paths = files(root, st, f"domains-{LIB}", f"l1/domains/library={LIB}/stage={st}/")
    dom = read(dom_paths, ["protein_uid", "model_acc"])
    copy = plan.get("copy") or {}
    shards = plan.get("shards", [])
    owned_block = plan.get("owned") or {}
    owned_stage_fps = [
        root / x["path"] for o in owned_block.get("stages", []) for x in o.get("feature_protein", [])
    ]
    if owned_stage_fps:
        con.execute(
            f"CREATE TEMP TABLE mine AS "
            f"SELECT protein_uid FROM {rp(ps_paths)} "
            f"UNION "
            f"SELECT fp.protein_uid FROM {rp(owned_stage_fps)} fp "
            f"SEMI JOIN (SELECT DISTINCT contig_uid FROM {rp(gi_rec)}) g USING (contig_uid)"
        )
    else:
        con.execute(f"CREATE TEMP TABLE mine AS SELECT protein_uid FROM {rp(ps_paths)}")
    mine_prot, = con.execute("SELECT count(*) FROM mine").fetchone()
    summary["pfam"] = {
        "proteins": plan.get("proteins"),
        "stage_protein_seq": n_prot,
        "shared_contig_proteins_added": int(mine_prot) - n_prot,
        "residues": plan.get("residues"),
        "shards": len(shards),
        "shard_residues": plan.get("shard_residues"),
        "scanned_proteins": sum(s["proteins"] for s in shards),
        "scanned_residues": sum(s["residues"] for s in shards),
        "copied_proteins": copy.get("proteins"),
        "copied_residues": copy.get("residues"),
        "copied_rows": copy.get("rows"),
        "copied_proteins_with_hit": copy.get("proteins_with_hit"),
        "sources": [{k: s.get(k) for k in ("stage", "proteins", "residues", "rows", "proteins_with_hit")}
                    for s in copy.get("sources", [])],
        "passed_over": copy.get("passed_over"),
        "domain_rows": dom.num_rows,
        "proteins_with_domain": pc.count_distinct(dom.column("protein_uid")).as_py(),
        "models_with_hit": pc.count_distinct(dom.column("model_acc")).as_py(),
    }
    if plan.get("proteins") is not None and plan["proteins"] != int(mine_prot):
        failed.append(f"scan plan holds {plan['proteins']} proteins, expected {int(mine_prot)} (protein_seq {n_prot})")
    say(f"pfam: {summary['pfam']}")

    cc: dict = {"sources": []}
    if copy.get("sources"):
        con.execute("CREATE TEMP TABLE claimed (protein_uid BLOB, src INTEGER)")
        for k, s in enumerate(copy["sources"]):
            src_ps = [root / x["path"] for x in s["protein_files"]]
            con.execute(
                f"INSERT INTO claimed SELECT m.protein_uid, {k} FROM mine m SEMI JOIN "
                f"(SELECT protein_uid FROM {rp(src_ps)}) USING (protein_uid) "
                f"WHERE m.protein_uid NOT IN (SELECT protein_uid FROM claimed)"
            )
        total_copied, = con.execute("SELECT count(*) FROM claimed").fetchone()
        cc["copied_proteins"] = total_copied
        cc["equals_plan"] = total_copied == copy.get("proteins")
        if not cc["equals_plan"]:
            failed.append(f"copied proteins: {total_copied} independent, {copy.get('proteins')} plan")
        cols = ", ".join(COLS)
        all_rows = 0
        for k, s in enumerate(copy["sources"]):
            src_dom = [root / x["path"] for x in s["files"]]
            mine_rows = (
                f"SELECT {cols} FROM {rp(dom_paths)} d SEMI JOIN "
                f"(SELECT protein_uid FROM claimed WHERE src = {k}) USING (protein_uid)"
            )
            src_rows = (
                f"SELECT {cols} FROM {rp(src_dom)} d SEMI JOIN "
                f"(SELECT protein_uid FROM claimed WHERE src = {k}) USING (protein_uid)"
            )
            n_mine, = con.execute(f"SELECT count(*) FROM ({mine_rows})").fetchone()
            n_src, = con.execute(f"SELECT count(*) FROM ({src_rows})").fetchone()
            d1, = con.execute(f"SELECT count(*) FROM (({mine_rows}) EXCEPT ALL ({src_rows}))").fetchone()
            d2, = con.execute(f"SELECT count(*) FROM (({src_rows}) EXCEPT ALL ({mine_rows}))").fetchone()
            n_p, = con.execute(f"SELECT count(*) FROM claimed WHERE src = {k}").fetchone()
            all_rows += n_mine
            row = {
                "stage": s["stage"], "proteins": n_p, "plan_proteins": s.get("proteins"),
                "rows_in_stage": n_mine, "rows_in_source": n_src, "plan_rows": s.get("rows"),
                "stage_not_in_source": d1, "source_not_in_stage": d2,
            }
            cc["sources"].append(row)
            if d1 or d2 or n_mine != n_src or n_p != s.get("proteins") or n_mine != s.get("rows"):
                failed.append(f"copied rows from {s['stage']} differ: {row}")
        cc["rows"] = all_rows
        cc["columns_compared"] = list(COLS)
    con.close()
    summary["copy_check"] = cc
    say(f"copy check: {cc}")

    # 6. CRISPR
    cr = {}
    for kind, prod in (("default", crispr.producer()), ("near_cas", crispr.near_cas_producer("38.2"))):
        step = f"crispr-{prod}"
        feats = files(root, st, step, f"l1/feature/producer={prod}/stage={st}/")
        elems = files(root, st, step, f"l1/crispr/producer={prod}/stage={st}/")
        ledp = files(root, st, step, f"ledger/stage={st}/step={step}/")
        lt = read(ledp, ["n_in", "n_kept", "n_rejected", "reason_code"]).to_pydict()
        c = {
            "producer": prod, "recorded": record(root, st, step) is not None,
            "arrays": rows_of(feats), "elements": rows_of(elems),
            "contigs_searched": 0, "contigs_skipped_other_stage": 0,
            "calls_in": 0, "calls_kept": 0, "rejected": {},
        }
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
            failed.append(f"{kind}: ledger keeps {c['calls_kept']}, feature files hold {c['arrays']}")
        if not c["recorded"]:
            failed.append(f"{kind}: no run record {step}")
        c["arrays_per_gbp"] = round(c["arrays"] / (bases_owned / 1e9), 3) if bases_owned else None
        cr[kind] = c
    summary["crispr"] = cr
    say(f"crispr: {cr}")

    # 7. the stage's steps, and the counts card C-39 fixes or estimates
    want_steps = ["inputs", "ingest", "genecall", f"crispr-{crispr.producer()}", "igr", "tapes",
                  f"domains-{LIB}", f"domains-{LIB}.plan", f"crispr-{crispr.near_cas_producer('38.2')}"]
    missing = [s for s in want_steps if record(root, st, s) is None]
    if missing:
        failed.append(f"missing run records: {missing}")
    observed = {k: summary["contigs"][k] for k in (
        "genomes", "genome_input_rows", "contigs", "stored", "stored_by_earlier", "ledger_skipped_earlier",
        "bases_stored", "circular")}
    observed.update({"genes": n_genes, "distinct_proteins": n_prot, "new_distinct_proteins": int(new_prot),
                     "residues": int(res_all or 0), "added_bytes": summary["added"]["bytes"]})
    if a.no_expect:
        summary["expected"] = {"skipped": "--no-expect"}
    else:
        ranges = est_mod.expected_ranges()
        problems = est_mod.check_counts(observed, ranges)
        if summary["contigs"]["stored_by_earlier_by_stage"] != est_mod.C39_STORED_BY:
            problems.append(f"stored_by_earlier_by_stage: {summary['contigs']['stored_by_earlier_by_stage']} "
                            f"!= {est_mod.C39_STORED_BY}")
        summary["expected"] = {"observed": observed, "ranges": ranges, "problems": problems}
        if problems:
            failed.append(f"counts outside C-39's ranges (AGENTS §6): {problems}")
    say(f"expected: {summary['expected']}")

    summary["failed"] = failed
    (out / "summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True, default=str))
    say("RESULT: " + ("all checks pass" if not failed else "FAILED: " + "; ".join(failed)))
    (out / "report.txt").write_text("\n".join(report) + "\n")
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
