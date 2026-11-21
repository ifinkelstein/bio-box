"""Card A-12: the check after a stage's build, for every stage (made from card C-47's `c47_check.py` and card C-48c's
`c48c_v1_check.py`), reading the stage's choices from its conf.

    python stage_check.py --conf F --stage S --root ROOT --before DIR --after DIR --objects OBJECTS.JSON
        --out DIR --tmp DIR --profile P [--expect-copy]

What it asserts (each a line of `failed`, exit 1 when any):
1. Before against after (the driver's snapshots, never a fresh listing): no file gone, none changed in size or
   mtime, no file of another stage changed by sha256 and none added, every new file under a `stage=<S>` component,
   no STOP left. `before` lists the root without the stage's own files, so a stage resumed after a stop is counted
   whole.
2. `root_added.tsv`: every new file with its size and sha256 (read in 1 MiB blocks), for the durable copy.
3. The stage's run records: no intent file, none over its CPU budget, the rate guard (`stage_est.py rates`, the
   conf's accepted records included) with no step more than 20% over its estimate, and every step of the stage's
   depth recorded (`stage.stage_incomplete`, the release train's own check).
4. Contigs and genes, in DuckDB with a memory cap and spilling (a stage of millions of contigs is never held in
   Python): the contigs `genome_input` lists are exactly those the stage stores and those another stage holds
   (DEC-26), none both, none held by two stages; the gene-call ledger's skips for another stage are exactly those
   contigs with that stage's name, and its within-stage rows (`owned_by:<S>`, `owned_by:<S>/<genome_uid>`, card
   C-47's correction: the owner's stage is compared as a whole name) are counted apart; every stored contig has gene
   rows or a ledger row; no contig has gene rows that the stage does not store; with a least contig length, no
   stored contig under it. New distinct proteins against every other stage's `protein_seq`.
5. With a library: the scan plan's proteins are the stage's distinct proteins and the proteins of contigs other
   stages own (card C-38h); the rows copied from other stages (card B-47b) equal theirs row for row, column for
   column; with --expect-copy (the rehearsal) at least one protein is copied. Without a library: no domain record.
6. With CRISPR: each producer's ledger keeps as many arrays as its feature files hold; the near-Cas producer only
   at depth annotated with Pfam. Without: no CRISPR record.

Every table is read by the explicit file list of the run record that published it (AGENTS §3.9). Writes into --out:
report.txt, summary.json, root_added.tsv, steps.tsv.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pyarrow.parquet as pq

COLS = (
    "protein_uid", "model_acc", "model_name", "model_len", "seq_bitscore", "seq_evalue_zref",
    "dom_bitscore", "dom_ievalue_zref", "env_start", "env_end", "ali_start", "ali_end",
    "hmm_from", "hmm_to", "cutoff", "z", "dom_z",
)
OWNED_BY = "owned_by:"  # biobox.l0seq.store.OWNED_BY: a skip's ledger reason is this + ContigOwners.owner
OWNER_SEP = "/"  # biobox.l0seq.store.OWNER_SEP: `<stage>/<genome_uid hex>`, a contig another genome of the stage owns


def _est_module():
    spec = importlib.util.spec_from_file_location("stage_est", Path(__file__).with_name("stage_est.py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules["stage_est"] = mod  # a dataclass looks its module up by name
    spec.loader.exec_module(mod)
    return mod


def q(s: str) -> str:
    return "'" + str(s).replace("'", "''") + "'"


def skip_kind_sql(col: str, stage: str) -> str:
    """SQL for how a gene-call ledger reason counts for `stage`: 'earlier' for a contig another stage stores
    (`owned_by:<other stage>`, DEC-26), 'within' for a contig called once under another genome of the stage
    (`owned_by:<stage>` or `owned_by:<stage>/<genome_uid>`, DEC-27) or a duplicate record, NULL otherwise (`ok`).
    The owner is compared as a whole stage name: `owned_by:C50` is not stage C5's (card C-47)."""
    n = len(OWNED_BY) + 1
    owner = f"split_part(substr({col}, {n}), {q(OWNER_SEP)}, 1)"
    return (f"(CASE WHEN {col} LIKE 'duplicate_contig%' THEN 'within' "
            f"WHEN {col} LIKE {q(OWNED_BY + '%')} THEN (CASE WHEN {owner} = {q(stage)} THEN 'within' ELSE 'earlier' "
            "END) ELSE NULL END)")


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def listing(d: Path) -> dict[str, tuple[int, str]]:
    """The snapshot's listing of files (thousands of lines: the root's files, not a table's rows)."""
    out = {}
    with open(d / "listing.tsv") as fh:
        for line in fh:
            path, size, mtime = line.rstrip("\n").split("\t")
            out[path] = (int(size), mtime)
    return out


def sums(d: Path) -> dict[str, str]:
    out = {}
    p = d / "earlier.sha256"
    if p.exists():
        with open(p) as fh:
            for line in fh:
                h, path = line.rstrip("\n").split("  ", 1)
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


def lst(ps) -> str:
    return "[" + ", ".join(q(p) for p in ps) + "]"


def rows_of(ps) -> int:
    return sum(pq.read_metadata(p).num_rows for p in ps)


def main() -> int:
    ap = argparse.ArgumentParser()
    for k in ("conf", "stage", "root", "before", "after", "objects", "out", "tmp", "profile"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--expect-copy", action="store_true")
    ap.add_argument("--before-est", default=None, help="the driver's estimate/before.json (its fasta_bytes)")
    a = ap.parse_args()
    est = _est_module()
    conf = est.read_conf(a.conf, a.stage)
    import duckdb

    from biobox import stage as bs
    from biobox.l0seq.store import INDEX_DIR
    from biobox.l1 import crispr
    from biobox.l1.genecall import GENE_CALL_VERSION as GCV

    root, out, st, tmp = Path(a.root), Path(a.out), conf.stage, Path(a.tmp)
    out.mkdir(parents=True, exist_ok=True)
    tmp.mkdir(parents=True, exist_ok=True)
    report: list[str] = []
    failed: list[str] = []
    summary: dict = {"stage": st, "depth": conf.depth, "min_contig_len": conf.min_contig_len,
                     "library": None if conf.library is None else " ".join(conf.library), "crispr": conf.crispr}

    def say(s: str) -> None:
        report.append(s)
        print(s, flush=True)

    # 1. before against after
    b, f = listing(Path(a.before)), listing(Path(a.after))
    gone = sorted(p for p in b if p not in f)
    changed = sorted(p for p in b if p in f and f[p] != b[p])
    new = sorted(p for p in f if p not in b)
    mine = [p for p in new if f"stage={st}" in p.split("/")]
    unexpected = sorted(set(new) - set(mine))
    sb, sa = sums(Path(a.before)), sums(Path(a.after))
    e_changed = sorted(p for p in sb if sa.get(p) != sb[p])
    e_new = sorted(p for p in sa if p not in sb)
    stop_left = [p for p in f if p == "STOP"]
    say(f"before: {len(b)} files; after: {len(f)}; new {len(new)} ({len(mine)} under stage={st}); gone {len(gone)}; "
        f"size or mtime changed {len(changed)}; other stages' files hashed before {len(sb)}, after {len(sa)}, sha256 "
        f"changed {len(e_changed)}, new {len(e_new)}")
    for what, items in (("gone", gone), ("changed", changed), ("other stage's sha256 changed", e_changed),
                        ("other stage's file new", e_new), ("unexpected new", unexpected), ("STOP left", stop_left)):
        if items:
            failed.append(f"{what}: {len(items)}")
            for p in items[:50]:
                say(f"  {what}: {p}")
    if not mine:
        failed.append(f"no new file under stage={st}")
    summary["files"] = {"before": len(b), "after": len(f), "new": len(new), "new_under_stage": len(mine),
                        "gone": gone, "changed": changed, "other_hashed": len(sb),
                        "other_sha256_changed": e_changed, "unexpected": unexpected}

    # 2. root_added.tsv
    total = 0
    with open(out / "root_added.tsv", "w") as fh:
        fh.write("path\tbytes\tsha256\n")
        for p in new:
            n = (root / p).stat().st_size
            total += n
            fh.write(f"{p}\t{n}\t{sha256(root / p)}\n")
            if f[p][0] != n:
                failed.append(f"{p}: {n} bytes now, {f[p][0]} in the after listing")
    summary["added"] = {"files": len(new), "bytes": total}
    say(f"root_added.tsv: {len(new)} files, {total:,} bytes")

    # 3. the stage's run records
    rows = []
    for p in sorted((root / "runs" / f"stage={st}").glob("*.json")):
        if p.name.endswith(".publishing.json"):
            failed.append(f"intent file {p.name}")
            continue
        r = json.loads(p.read_text())
        rows.append({"step": p.name[:-5], "cpu_min": round(r["cpu_min"], 3), "wall_s": round(r["wall_s"], 1),
                     "peak_rss_mib": round((r.get("peak_rss_bytes") or 0) / 2**20), "threads": r.get("threads"),
                     "started_at": r.get("started_at"), "over_cpu_budget": r.get("over_cpu_budget"),
                     "image_digest": r.get("image_digest"), "code_id": r.get("code_id"),
                     "code_commit": r.get("code_commit"), "outputs": len(r.get("outputs", []))})
    cmd = [sys.executable, str(Path(__file__).with_name("stage_est.py")), "rates", "--conf", a.conf, "--stage", st,
           "--root", str(root), "--profile", a.profile, "--objects", a.objects, "--tmp", str(tmp / "rates")]
    if a.before_est:
        cmd += ["--before", a.before_est]
    ran = subprocess.run(cmd, capture_output=True, text=True)
    try:
        rates = json.loads(ran.stdout)
    except ValueError:
        rates = {"error": ran.stderr[-2000:]}
    if ran.returncode not in (0, 3) or "error" in rates:
        failed.append(f"stage_est.py rates exited {ran.returncode}: {ran.stderr[-400:]}")
    if rates.get("over_20pct"):
        failed.append(f"rates over +20% (DEC-53): {rates['over_20pct']}")
    by = {x["step"]: x for x in rates.get("rows", [])}
    for r in rows:
        if r["step"] in by:
            r.update(estimate_cpu_min=by[r["step"]]["estimate"], ratio=by[r["step"]]["ratio"],
                     judged=by[r["step"]]["judged"])
    over_budget = [r["step"] for r in rows if r["over_cpu_budget"]]
    if over_budget:
        failed.append(f"over_cpu_budget: {over_budget}")
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(out / "steps.tsv", "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")
    libs = [] if conf.library is None else [(conf.library[0], conf.library[1])]
    incomplete = bs.stage_incomplete(root, st, conf.depth, libraries=libs, crispr=conf.crispr)
    if incomplete:
        failed.append(f"the stage is not complete at depth {conf.depth}: {incomplete[:10]}")
    summary["steps"] = {"records": len(rows), "names": [r["step"] for r in rows if ".shard-" not in r["step"]],
                        "cpu_min_total": round(sum(r["cpu_min"] for r in rows), 2),
                        "cpu_min_max": max((r["cpu_min"] for r in rows), default=None),
                        "largest": max(rows, key=lambda r: r["cpu_min"])["step"] if rows else None,
                        "over_budget": over_budget, "images": sorted({str(r["image_digest"]) for r in rows}),
                        "code_ids": sorted({str(r["code_id"]) for r in rows}), "incomplete": incomplete,
                        "rates": rates}
    for r in rows:
        say(f"step {r['step']}: {r['cpu_min']} CPU-min, {r['wall_s']} s, {r['peak_rss_mib']} MiB, {r['threads']} "
            "threads" + (f", estimate {r['estimate_cpu_min']} (ratio {r['ratio']}"
                         f"{'' if r['judged'] else ', not judged: DEC-43'})" if "ratio" in r else ""))

    # 4. contigs and genes, in DuckDB
    con = duckdb.connect(config={"threads": 2, "memory_limit": "8GB", "temp_directory": str(tmp / "duckdb"),
                                 "preserve_insertion_order": False})
    rp = lambda ps: f"read_parquet({lst(ps)}, hive_partitioning = false)"  # noqa: E731
    empty = lambda cols: "SELECT " + ", ".join(f"NULL::{t} AS {c}" for c, t in cols) + " LIMIT 0"  # noqa: E731

    def src(ps, cols):
        return f"(SELECT {', '.join(c for c, _ in cols)} FROM {rp(ps)})" if ps else f"({empty(cols)})"

    gi_files = files(root, st, "inputs", f"l0seq/genome_input/stage={st}/")
    ci_files = files(root, st, "ingest", f"{INDEX_DIR}/stage={st}/")
    if not gi_files or not ci_files:
        failed.append("no genome_input or contig_index files in the inputs and ingest records")
    base = root / INDEX_DIR
    others = sorted(d.name.removeprefix("stage=") for d in base.glob("stage=*") if d.is_dir() and d.name != f"stage={st}")
    con.execute(f"CREATE TEMP TABLE gi AS SELECT contig_uid, genome_uid, kind, genetic_code, is_circular, "
                f"assembly_quality FROM {src(gi_files, [('contig_uid', 'BLOB'), ('genome_uid', 'BLOB'), ('kind', 'VARCHAR'), ('genetic_code', 'INTEGER'), ('is_circular', 'BOOLEAN'), ('assembly_quality', 'VARCHAR')])}")
    con.execute(f"CREATE TEMP TABLE ci AS SELECT contig_uid, contig_len FROM "
                f"{src(ci_files, [('contig_uid', 'BLOB'), ('contig_len', 'BIGINT'), ('is_default_variant', 'BOOLEAN')])} "
                "WHERE is_default_variant")
    parts = []
    for s in others:
        fs = files(root, s, "ingest", f"{INDEX_DIR}/stage={s}/")
        if fs:
            parts.append(f"SELECT contig_uid, {q(s)} AS st FROM {rp(fs)} WHERE is_default_variant")
    held_src = " UNION ALL ".join(parts) if parts else empty([("contig_uid", "BLOB"), ("st", "VARCHAR")])
    con.execute(f"CREATE TEMP TABLE held AS SELECT contig_uid, min(st) AS st, count(DISTINCT st) AS m FROM ({held_src}) "
                "SEMI JOIN (SELECT DISTINCT contig_uid FROM gi) USING (contig_uid) GROUP BY contig_uid")
    led_files = files(root, st, "genecall", f"ledger/stage={st}/step=genecall/")
    con.execute(f"CREATE TEMP TABLE led AS SELECT object_key AS contig_uid, reason_code, "
                f"{skip_kind_sql('reason_code', st)} AS kind FROM "
                f"{src(led_files, [('object_key', 'BLOB'), ('reason_code', 'VARCHAR')])}")
    one = lambda sql: con.execute(sql).fetchall()[0][0]  # noqa: E731  (one row, fetched whole)
    c = {
        "genome_input_rows": one("SELECT count(*) FROM gi"),
        "genomes": one("SELECT count(DISTINCT genome_uid) FROM gi"),
        "contigs": one("SELECT count(DISTINCT contig_uid) FROM gi"),
        "stored": one("SELECT count(*) FROM ci"),
        "bases_stored": int(one("SELECT coalesce(sum(contig_len), 0) FROM ci")),
        "held_by_other_stages": one("SELECT count(*) FROM held"),
        "held_by_stage": dict(con.execute("SELECT st, count(*) FROM held GROUP BY st ORDER BY st").fetchall()),
        "ledger_skipped_other": one("SELECT count(*) FROM led WHERE kind = 'earlier'"),
        "ledger_skipped_within": one("SELECT count(*) FROM led WHERE kind = 'within'"),
        "kinds": [r[0] for r in con.execute("SELECT DISTINCT kind FROM gi ORDER BY 1").fetchall()],
        "codes": [r[0] for r in con.execute("SELECT DISTINCT genetic_code FROM gi ORDER BY 1").fetchall()],
        "quality": dict(con.execute("SELECT assembly_quality, count(*) FROM gi GROUP BY 1 ORDER BY 1").fetchall()),
        "circular_rows": one("SELECT count(*) FROM gi WHERE is_circular"),
    }
    checks = {
        "contigs held by two other stages (DEC-26)": "SELECT count(*) FROM held WHERE m > 1",
        "contigs both stored by the stage and held by another": "SELECT count(*) FROM ci SEMI JOIN held USING (contig_uid)",
        "genome_input contigs neither stored nor held":
            "SELECT count(*) FROM (SELECT DISTINCT contig_uid FROM gi) ANTI JOIN ci USING (contig_uid) "
            "ANTI JOIN held USING (contig_uid)",
        "stored contigs not in genome_input":
            "SELECT count(*) FROM ci ANTI JOIN (SELECT DISTINCT contig_uid FROM gi) g USING (contig_uid)",
        "ledger skips for another stage that differ from the stage holding the contig":
            "SELECT count(*) FROM (SELECT contig_uid, split_part(substr(reason_code, "
            f"{len(OWNED_BY) + 1}), {q(OWNER_SEP)}, 1) AS st FROM led WHERE kind = 'earlier') l "
            "FULL OUTER JOIN held h USING (contig_uid) WHERE l.st IS DISTINCT FROM h.st",
    }
    if conf.min_contig_len is not None:
        checks[f"stored contigs under the least length {conf.min_contig_len}"] = \
            f"SELECT count(*) FROM ci WHERE contig_len < {int(conf.min_contig_len)}"
    for what, sql in checks.items():
        n = one(sql)
        c[what] = n
        if n:
            failed.append(f"{what}: {n}")
    ing_led = files(root, st, "ingest", f"ledger/stage={st}/step=ingest/")
    if ing_led:
        c["ingest_ledger"] = {r[0]: {"rows": r[1], "bases_rejected": int(r[2] or 0)} for r in con.execute(
            f"SELECT split_part(reason_code, ':', 1), count(*), sum(n_bases_rejected) FROM {rp(ing_led)} "
            "GROUP BY 1 ORDER BY 1").fetchall()} if "n_bases_rejected" in pq.read_schema(ing_led[0]).names else \
            dict(con.execute(f"SELECT split_part(reason_code, ':', 1), count(*) FROM {rp(ing_led)} GROUP BY 1 "
                             "ORDER BY 1").fetchall())
    occ = files(root, st, "ingest", f"l0seq/contig_occurrence/stage={st}/")
    c["contig_occurrence_rows"] = rows_of(occ)
    c["contig_label_rows"] = rows_of(files(root, st, "ingest", f"l0seq/contig_label/stage={st}/"))
    summary["contigs"] = c
    say(f"contigs: {c}")
    fp = files(root, st, "genecall", f"l1/feature_protein/gene_call_version={GCV}/stage={st}/")
    ps = files(root, st, "genecall", f"l0seq/protein_seq/gene_call_version={GCV}/stage={st}/")
    n_genes, n_prot = rows_of(fp), rows_of(ps)
    con.execute(f"CREATE TEMP TABLE with_genes AS SELECT DISTINCT contig_uid FROM "
                f"{src(fp, [('contig_uid', 'BLOB')])}")
    cov = {
        "stored_with_genes": one("SELECT count(*) FROM ci SEMI JOIN with_genes USING (contig_uid)"),
        "stored_without_genes_listed": dict(con.execute(
            "SELECT reason_code, count(*) FROM (SELECT * FROM ci ANTI JOIN with_genes USING (contig_uid)) x "
            "JOIN (SELECT contig_uid, min(reason_code) AS reason_code FROM led GROUP BY contig_uid) l "
            "USING (contig_uid) GROUP BY 1 ORDER BY 1").fetchall()),
        "stored_without_genes_or_ledger_row": one(
            "SELECT count(*) FROM ci ANTI JOIN with_genes USING (contig_uid) ANTI JOIN led USING (contig_uid)"),
        "contigs_with_genes_not_stored_here": one("SELECT count(*) FROM with_genes ANTI JOIN ci USING (contig_uid)"),
    }
    if cov["stored_without_genes_or_ledger_row"]:
        failed.append(f"{cov['stored_without_genes_or_ledger_row']} stored contigs have no gene row and no ledger row")
    if cov["contigs_with_genes_not_stored_here"]:
        failed.append(f"{cov['contigs_with_genes_not_stored_here']} contigs have gene rows in {st} but are not stored "
                      "by it")
    summary["coverage"] = cov
    say(f"coverage: {cov}")
    other_ps = [p for s in others for p in files(root, s, "genecall",
                                                 f"l0seq/protein_seq/gene_call_version={GCV}/stage={s}/")]
    pcols = [("protein_uid", "BLOB"), ("length", "BIGINT")]
    new_prot, new_res = con.execute(f"SELECT count(*), coalesce(sum(length), 0) FROM {src(ps, pcols)} x "
                                    f"ANTI JOIN {src(other_ps, [('protein_uid', 'BLOB')])} y USING (protein_uid)"
                                    ).fetchall()[0]
    res_all = one(f"SELECT coalesce(sum(length), 0) FROM {src(ps, pcols)} x")
    summary["genes"] = {"genes": n_genes, "distinct_proteins": n_prot, "residues": int(res_all),
                        "new_distinct_proteins": int(new_prot), "new_residues": int(new_res),
                        "in_other_stages": n_prot - int(new_prot), "other_protein_seq_files": len(other_ps)}
    say(f"genes {n_genes:,}; distinct proteins {n_prot:,} ({int(res_all):,} residues); new against "
        f"{','.join(others) or 'none'}: {int(new_prot):,} ({int(new_res):,} residues)")
    if n_genes == 0:
        failed.append("no gene rows")

    # 5. domains
    if conf.library is None:
        recs = [p.name for p in (root / "runs" / f"stage={st}").glob("domains-*.json")]
        if recs:
            failed.append(f"domain records without a library in the conf: {recs}")
        summary["pfam"] = None
        say("pfam: no library in the conf; no domain record")
    else:
        lib = f"{conf.library[0]}-{conf.library[1]}"
        plan_p = root / f"l1/domain_scan/library={lib}/stage={st}/plan.json"
        plan = json.loads(plan_p.read_text()) if plan_p.exists() else {}
        dom = files(root, st, f"domains-{lib}", f"l1/domains/library={lib}/stage={st}/")
        copy = plan.get("copy") or {}
        shards = plan.get("shards", [])
        owned_fp = [root / x["path"] for o in (plan.get("owned") or {}).get("stages", [])
                    for x in o.get("feature_protein", [])]
        if owned_fp:
            con.execute(f"CREATE TEMP TABLE mine AS SELECT protein_uid FROM {src(ps, [('protein_uid', 'BLOB')])} x "
                        f"UNION SELECT fp.protein_uid FROM {rp(owned_fp)} fp "
                        "SEMI JOIN (SELECT DISTINCT contig_uid FROM gi) g USING (contig_uid)")
        else:
            con.execute(f"CREATE TEMP TABLE mine AS SELECT DISTINCT protein_uid FROM "
                        f"{src(ps, [('protein_uid', 'BLOB')])} x")
        mine_prot = one("SELECT count(*) FROM mine")
        dcols = [("protein_uid", "BLOB"), ("model_acc", "VARCHAR")]
        pf = {"proteins": plan.get("proteins"), "stage_protein_seq": n_prot,
              "shared_contig_proteins_added": int(mine_prot) - n_prot, "residues": plan.get("residues"),
              "shards": len(shards), "scanned_proteins": sum(s["proteins"] for s in shards),
              "scanned_residues": sum(s["residues"] for s in shards), "copied_proteins": copy.get("proteins"),
              "copied_rows": copy.get("rows"),
              "sources": [{k: s.get(k) for k in ("stage", "proteins", "residues", "rows", "proteins_with_hit")}
                          for s in copy.get("sources", [])],
              "domain_rows": rows_of(dom),
              "proteins_with_domain": one(f"SELECT count(DISTINCT protein_uid) FROM {src(dom, dcols)} x"),
              "models_with_hit": one(f"SELECT count(DISTINCT model_acc) FROM {src(dom, dcols)} x")}
        if not dom:
            failed.append(f"no domain files in the record domains-{lib}")
        if plan.get("proteins") is not None and plan["proteins"] != int(mine_prot):
            failed.append(f"scan plan holds {plan['proteins']} proteins, expected {int(mine_prot)}")
        if a.expect_copy and not (copy.get("proteins") or 0) > 0:
            failed.append("the scan plan copied no protein (card B-47b); the rehearsal's fixture shares proteins")
        cc: dict = {"sources": []}
        if copy.get("sources") and dom:  # without domain files the failure is already listed
            con.execute("CREATE TEMP TABLE claimed (protein_uid BLOB, src INTEGER)")
            for k, s in enumerate(copy["sources"]):
                con.execute(f"INSERT INTO claimed SELECT m.protein_uid, {k} FROM mine m SEMI JOIN "
                            f"(SELECT protein_uid FROM {rp([root / x['path'] for x in s['protein_files']])}) "
                            "USING (protein_uid) WHERE m.protein_uid NOT IN (SELECT protein_uid FROM claimed)")
            cc["copied_proteins"] = one("SELECT count(*) FROM claimed")
            if cc["copied_proteins"] != copy.get("proteins"):
                failed.append(f"copied proteins: {cc['copied_proteins']} independent, {copy.get('proteins')} plan")
            cols = ", ".join(COLS)
            for k, s in enumerate(copy["sources"]):
                src_dom = [root / x["path"] for x in s["files"]]
                m_rows = (f"SELECT {cols} FROM {rp(dom)} d SEMI JOIN (SELECT protein_uid FROM claimed WHERE src = {k}) "
                          "USING (protein_uid)")
                s_rows = (f"SELECT {cols} FROM {rp(src_dom)} d SEMI JOIN (SELECT protein_uid FROM claimed WHERE "
                          f"src = {k}) USING (protein_uid)")
                row = {"stage": s["stage"], "proteins": one(f"SELECT count(*) FROM claimed WHERE src = {k}"),
                       "plan_proteins": s.get("proteins"), "rows_in_stage": one(f"SELECT count(*) FROM ({m_rows})"),
                       "rows_in_source": one(f"SELECT count(*) FROM ({s_rows})"), "plan_rows": s.get("rows"),
                       "stage_not_in_source": one(f"SELECT count(*) FROM (({m_rows}) EXCEPT ALL ({s_rows}))"),
                       "source_not_in_stage": one(f"SELECT count(*) FROM (({s_rows}) EXCEPT ALL ({m_rows}))")}
                cc["sources"].append(row)
                if (row["stage_not_in_source"] or row["source_not_in_stage"] or row["rows_in_stage"] !=
                        row["rows_in_source"] or row["proteins"] != s.get("proteins")
                        or row["rows_in_stage"] != s.get("rows")):
                    failed.append(f"copied rows from {s['stage']} differ: {row}")
            cc["columns_compared"] = list(COLS)
        pf["copy_check"] = cc
        summary["pfam"] = pf
        say(f"pfam: {pf}")

    # 6. CRISPR
    near_ok = conf.library is not None and conf.depth == "annotated" and conf.library[0] == crispr.CAS_LIBRARY
    want = []
    if conf.crispr:
        want.append(("default", crispr.producer()))
        if near_ok:
            want.append(("near_cas", crispr.near_cas_producer(conf.library[1])))
    stray = [p.name for p in (root / "runs" / f"stage={st}").glob("crispr-*.json")
             if p.name[len("crispr-"):-len(".json")] not in {prod for _, prod in want}]
    if stray:
        failed.append(f"CRISPR records the conf does not ask for: {stray}")
    cr = {}
    for kind, prod in want:
        step = f"crispr-{prod}"
        feats = files(root, st, step, f"l1/feature/producer={prod}/stage={st}/")
        ledp = files(root, st, step, f"ledger/stage={st}/step={step}/")
        lcols = [("n_in", "BIGINT"), ("n_kept", "BIGINT"), ("n_rejected", "BIGINT"), ("reason_code", "VARCHAR")]
        searched, calls_in, kept, skipped = con.execute(
            "SELECT count(*) FILTER (WHERE reason_code = 'ok'), coalesce(sum(n_in) FILTER (WHERE reason_code = 'ok'), 0), "
            "coalesce(sum(n_kept) FILTER (WHERE reason_code = 'ok'), 0), "
            f"count(*) FILTER (WHERE reason_code LIKE {q(OWNED_BY + '%')}) FROM {src(ledp, lcols)} x").fetchall()[0]
        rejected = dict(con.execute(f"SELECT reason_code, sum(n_rejected) FROM {src(ledp, lcols)} x WHERE reason_code "
                                    f"<> 'ok' AND reason_code NOT LIKE {q(OWNED_BY + '%')} GROUP BY 1 ORDER BY 1"
                                    ).fetchall())
        x = {"producer": prod, "recorded": record(root, st, step) is not None, "arrays": rows_of(feats),
             "elements": rows_of(files(root, st, step, f"l1/crispr/producer={prod}/stage={st}/")),
             "contigs_searched": searched, "contigs_skipped_other_stage": skipped, "calls_in": int(calls_in),
             "calls_kept": int(kept), "rejected": rejected}
        if x["calls_kept"] != x["arrays"]:
            failed.append(f"{kind}: ledger keeps {x['calls_kept']}, feature files hold {x['arrays']}")
        if not x["recorded"]:
            failed.append(f"{kind}: no run record {step}")
        x["arrays_per_gbp"] = round(x["arrays"] / (c["bases_stored"] / 1e9), 3) if c["bases_stored"] else None
        cr[kind] = x
    summary["crispr"] = cr or None
    say(f"crispr: {cr or 'not asked for'}")
    con.close()

    summary["failed"] = failed
    (out / "summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True, default=str))
    say("RESULT: " + ("all checks pass" if not failed else "FAILED: " + "; ".join(failed)))
    (out / "report.txt").write_text("\n".join(report) + "\n")
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
