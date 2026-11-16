"""Card C-35: the checks after release r26.11.1 was built into the data root (or the rehearsal's release into its scratch
root). In the runtime image (pyarrow, duckdb; the tree's src on PYTHONPATH), read-only on the root:

    python c35_check.py --root ROOT --before LISTING --after LISTING --release R --stages S1,...,PL1 --earlier rM.2
        --earlier-stages S1,...,P6 --est BEFORE.json --out DIR --tmp DIR [--dec27]

Every table is read by the explicit file list of the release's manifest (AGENTS §3.9; no glob over a table).

1. Before against after (each LISTING: path, bytes, mtime of every file of the root): no file that existed before is
   missing or has another size or mtime; every new file names the release in its path (`release=<R>` or `stage=<R>`
   as a component, or `<dir>/release=<R>.<ext>`): no stage file and no earlier release's file was written.
2. `root_added.tsv`: every new file with its bytes and sha256 (for the durable copy), and the totals.
3. The manifest: per stage, the tables it lists (with `contig_label` for every stage) and the CRISPR arrays per producer
   (the `feature` rows of each `minced-*` producer, and `l1_crispr` rows).
4. Clusters (DuckDB over `member_cluster` of both releases, under a memory limit, spilling to --tmp): members and
   clusters per level of each release; per stage of the release, its distinct proteins, those not in the earlier
   release, and those whose cluster at NR90 / NR50 holds no protein of the earlier release ("new at NR90 / NR50"),
   with the count of such clusters.
5. cluster_uid against the earlier release, per level: each earlier cluster by what became of its members: all in one
   cluster of the same members (that cluster must have the same cluster_uid: the uid is a hash of the member set and the
   parameters), in one cluster with new proteins only added, in one cluster with other earlier proteins, or split; and
   how many earlier uids appear in the release at all.
6. With --dec27 (DEC-27, the checks card C-26 moved here): genomes by an independent DuckDB query over the stages'
   `genome_input` files against the catalog's `genome_member` and `genome_dim`; contigs held by two or more genomes and
   those whose genomes span PL1 and an earlier stage; then on one such contig a pair of NR50 tokens on adjacent slots:
   the catalog's `q1_count` and `q1_genomes` against a count made here from `tape_tokens` and `genome_input` (loci
   counted once per contig and slot, genomes once each).

Writes into --out: report.txt, summary.json, root_added.tsv. Exit 0 when the checks of 1, 5 (same members, same uid)
and 6 (the q1 equalities) pass, 1 otherwise.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import duckdb


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def listing(p: Path) -> dict[str, tuple[int, str]]:
    out = {}
    for line in p.read_text().splitlines():
        path, size, mtime = line.split("\t")
        out[path] = (int(size), mtime)
    return out


def lst(ps) -> str:
    return "[" + ", ".join("'" + str(p).replace("'", "''") + "'" for p in ps) + "]"


def names_release(path: str, rel: str) -> bool:
    parts = path.split("/")
    return (f"release={rel}" in parts or f"stage={rel}" in parts
            or any(p.startswith(f"release={rel}.") for p in parts))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--before", required=True, type=Path)
    ap.add_argument("--after", required=True, type=Path)
    ap.add_argument("--release", required=True)
    ap.add_argument("--stages", required=True)
    ap.add_argument("--earlier", required=True)
    ap.add_argument("--earlier-stages", required=True)
    ap.add_argument("--est", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--tmp", required=True, type=Path)
    ap.add_argument("--dec27", action="store_true")
    ap.add_argument("--memory", default="10GB")
    a = ap.parse_args()
    root, rel, early = a.root, a.release, a.earlier
    stages, estages = a.stages.split(","), a.earlier_stages.split(",")
    a.out.mkdir(parents=True, exist_ok=True)
    a.tmp.mkdir(parents=True, exist_ok=True)
    report: list[str] = []
    failed: list[str] = []
    summary: dict = {"release": rel, "earlier": early}

    def say(s: str) -> None:
        report.append(s)
        print(s, flush=True)

    # 1. before against after
    b0, b1 = listing(a.before), listing(a.after)
    gone = sorted(p for p in b0 if p not in b1)
    changed = sorted(p for p in b0 if p in b1 and b0[p] != b1[p])
    new = sorted(p for p in b1 if p not in b0)
    stray = [p for p in new if not names_release(p, rel)]
    say(f"1. listing: {len(b0)} files before, {len(b1)} after; {len(new)} new, {len(gone)} gone, {len(changed)} "
        f"changed in size or mtime; {len(stray)} new files not naming {rel}")
    for p in (gone + changed + stray)[:20]:
        say(f"   PROBLEM: {p}")
    if gone or changed or stray:
        failed.append("listing")
    # 2. root_added.tsv
    tot = 0
    with open(a.out / "root_added.tsv", "w") as fh:
        fh.write("path\tbytes\tsha256\n")
        for p in new:
            n = (root / p).stat().st_size
            tot += n
            fh.write(f"{p}\t{n}\t{sha256(root / p)}\n")
    summary["added"] = {"files": len(new), "bytes": tot}
    say(f"2. root_added.tsv: {len(new)} files, {tot:,} bytes")
    by_top: dict[str, list[int]] = {}
    for p in new:
        k = "/".join(p.split("/")[:2])
        v = by_top.setdefault(k, [0, 0])
        v[0] += 1
        v[1] += b1[p][0]
    summary["added_by_dir"] = by_top
    for k, (n, s) in sorted(by_top.items()):
        say(f"   {k}: {n} files, {s:,} bytes")

    # 3. the manifest
    man = json.loads((root / "manifests" / f"release={rel}.json").read_text())
    files = man["files"]

    def of(table: str, **parts: str) -> list[Path]:
        return [root / f["path"] for f in files
                if f["table"] == table and all(f["partitions"].get(k) == v for k, v in parts.items())]

    def rows(table: str, **parts: str) -> int:
        return sum(f["rows"] for f in files
                   if f["table"] == table and all(f["partitions"].get(k) == v for k, v in parts.items()))

    summary["manifest"] = {"stages": man["stages"], "files": len(files), "artifacts": len(man["artifacts"])}
    say(f"3. manifest: stages {man['stages']}, {len(files)} table files, {len(man['artifacts'])} artifacts")
    if man["stages"] != stages:
        failed.append("manifest stages")
    labels, arrays = {}, {}
    for s in stages:
        labels[s] = rows("contig_label", stage=s)
        prods = sorted({f["partitions"].get("producer") for f in files if f["table"] == "feature"
                        and f["partitions"].get("stage") == s and str(f["partitions"].get("producer")).startswith("minced")})
        arrays[s] = {p: {"arrays": rows("feature", stage=s, producer=p), "l1_crispr_rows": rows("l1_crispr", stage=s, producer=p)}
                     for p in prods}
    summary["contig_label_rows"] = labels
    summary["arrays"] = arrays
    say(f"   contig_label rows per stage: {labels}" + ("" if all(labels.values()) else "  PROBLEM: a stage without labels"))
    if not all(labels.values()):
        failed.append("labels")
    totals: dict[str, int] = {}
    for s in stages:
        for p, v in arrays[s].items():
            totals[p] = totals.get(p, 0) + v["arrays"]
        say(f"   arrays {s}: " + "; ".join(f"{p}: {v['arrays']} arrays, {v['l1_crispr_rows']} l1_crispr rows"
                                         for p, v in arrays[s].items()))
    summary["arrays_total"] = totals
    say(f"   arrays in all: {totals}")

    # 4. clusters
    con = duckdb.connect(config={"threads": 2, "memory_limit": a.memory, "temp_directory": str(a.tmp),
                                 "preserve_insertion_order": False})
    emanf = json.loads((root / "manifests" / f"release={early}.json").read_text())["files"]
    mb = of("member_cluster")
    ma = [root / f["path"] for f in emanf if f["table"] == "member_cluster"]
    con.execute(f"CREATE VIEW b AS SELECT member_uid, level, cluster_uid FROM read_parquet({lst(mb)}, hive_partitioning = false)")
    con.execute(f"CREATE VIEW a AS SELECT member_uid, level, cluster_uid FROM read_parquet({lst(ma)}, hive_partitioning = false)")
    counts = {}
    for name in ("a", "b"):
        for lvl, m, c in con.execute(f"SELECT level, count(*), count(DISTINCT cluster_uid) FROM {name} GROUP BY level "
                                     "ORDER BY level").fetchall():
            counts.setdefault(early if name == "a" else rel, {})[lvl] = {"members": m, "clusters": c}
    summary["clusters"] = counts
    say(f"4. clusters: {json.dumps(counts)}")
    levels = sorted(counts.get(rel, {}))
    per_stage: dict = {}
    stab: dict = {}
    for lvl in levels:
        con.execute("DROP TABLE IF EXISTS bc")
        con.execute(f"""CREATE TEMP TABLE bc AS
            SELECT b.cluster_uid AS cb, count(*) AS nb, count(a.member_uid) AS nold
            FROM b LEFT JOIN a ON a.member_uid = b.member_uid AND a.level = '{lvl}'
            WHERE b.level = '{lvl}' GROUP BY b.cluster_uid""")
        nnew = con.execute("SELECT count(*), coalesce(sum(nb), 0) FROM bc WHERE nold = 0").fetchone()
        say(f"   {lvl}: {nnew[0]:,} clusters of the release hold no protein of {early} ({nnew[1]:,} members)")
        stab.setdefault(lvl, {})["new_clusters"] = nnew[0]
        stab[lvl]["new_cluster_members"] = nnew[1]
        for s in stages:
            ps = of("protein_seq", stage=s)
            r = con.execute(f"""
                SELECT count(*), count(*) FILTER (WHERE x.member_uid IS NULL),
                       count(*) FILTER (WHERE bc.nold = 0), count(DISTINCT CASE WHEN bc.nold = 0 THEN bc.cb END)
                FROM (SELECT DISTINCT protein_uid FROM read_parquet({lst(ps)}, hive_partitioning = false)) p
                JOIN b ON b.member_uid = p.protein_uid AND b.level = '{lvl}'
                JOIN bc ON bc.cb = b.cluster_uid
                LEFT JOIN (SELECT member_uid FROM a WHERE level = '{lvl}') x ON x.member_uid = p.protein_uid""").fetchone()
            d = per_stage.setdefault(s, {"proteins": r[0], "not_in_earlier": r[1]})
            d[f"new_at_{lvl}"] = r[2]
            d[f"new_{lvl}_clusters_with_a_member"] = r[3]
        # 5. cluster_uid against the earlier release
        r = con.execute(f"""
            WITH j AS (SELECT a.cluster_uid AS ca, b.cluster_uid AS cb FROM a LEFT JOIN b
                       ON b.member_uid = a.member_uid AND b.level = '{lvl}' WHERE a.level = '{lvl}'),
                 ab AS (SELECT ca, count(*) AS na, count(cb) AS present, count(DISTINCT cb) AS ncb, min(cb) AS cb1
                        FROM j GROUP BY ca)
            SELECT count(*), sum(na),
                   count(*) FILTER (WHERE present < na),
                   count(*) FILTER (WHERE ncb = 1 AND present = na AND bc.nb = na),
                   coalesce(sum(na) FILTER (WHERE ncb = 1 AND present = na AND bc.nb = na), 0),
                   count(*) FILTER (WHERE ncb = 1 AND present = na AND bc.nb = na AND cb1 = ca),
                   count(*) FILTER (WHERE ncb = 1 AND present = na AND bc.nb > na AND bc.nold = na),
                   count(*) FILTER (WHERE ncb = 1 AND present = na AND bc.nold > na),
                   count(*) FILTER (WHERE ncb > 1),
                   count(*) FILTER (WHERE cb1 = ca),
                   count(*) FILTER (WHERE cb1 = ca AND NOT (ncb = 1 AND present = na AND bc.nb = na))
            FROM ab LEFT JOIN bc ON bc.cb = ab.cb1""").fetchone()
        keys = ("earlier_clusters", "earlier_members", "lost_members", "same_members", "same_members_proteins",
                "same_members_same_uid", "grew_with_new_only", "merged_with_earlier", "split", "uid_kept",
                "uid_kept_other_members")
        st = dict(zip(keys, r))
        stab[lvl].update(st)
        say(f"5. {lvl}: {json.dumps(st)}")
        if st["same_members"] != st["same_members_same_uid"] or st["uid_kept_other_members"] or st["lost_members"]:
            failed.append(f"uid {lvl}")
    summary["per_stage"] = per_stage
    summary["uid_stability"] = stab
    for s in stages:
        say(f"   {s}: {json.dumps(per_stage.get(s))}")

    # 6. DEC-27
    if a.dec27:
        gi = {s: of("genome_input", stage=s) for s in stages}
        con.execute("CREATE VIEW gi AS " + " UNION ALL ".join(
            f"SELECT '{s}' AS stage, contig_uid, genome_uid FROM read_parquet({lst(ps)}, hive_partitioning = false)"
            for s, ps in gi.items() if ps))
        g_ind = con.execute("SELECT count(DISTINCT genome_uid) FROM gi").fetchone()[0]
        multi = con.execute("""SELECT count(*), count(*) FILTER (WHERE has_pl AND has_old) FROM (
            SELECT contig_uid, bool_or(stage = 'PL1') AS has_pl, bool_or(stage <> 'PL1') AS has_old
            FROM gi GROUP BY contig_uid HAVING count(DISTINCT genome_uid) > 1)""").fetchone()
        from biobox import catalog
        cat = catalog.open_catalog(root / "catalog" / f"release={rel}.duckdb", threads=2)
        g_cat = cat.execute("SELECT count(DISTINCT genome_uid) FROM genome_member").fetchone()[0]
        g_dim = rows("genome_dim")
        say(f"6. DEC-27: genomes {g_ind:,} (genome_input, independent), {g_cat:,} (catalog genome_member), "
            f"{g_dim:,} (genome_dim rows); contigs held by 2 or more genomes {multi[0]:,}, of them spanning PL1 and "
            f"an earlier stage {multi[1]:,}")
        d27 = {"genomes_independent": g_ind, "genomes_genome_member": g_cat, "genome_dim_rows": g_dim,
               "contigs_multi_genome": multi[0], "contigs_pl1_and_earlier": multi[1]}
        if g_ind != g_cat or g_ind != g_dim:
            failed.append("genome counts")
        tt = of("tape_tokens")
        pick = con.execute(f"""
            WITH s AS (SELECT contig_uid FROM gi GROUP BY contig_uid
                       HAVING bool_or(stage = 'PL1') AND bool_or(stage <> 'PL1') AND count(DISTINCT genome_uid) > 1)
            SELECT t.contig_uid, t.nr50_tokens[2], t.nr50_tokens[3] FROM read_parquet({lst(tt)}, hive_partitioning = false) t
            WHERE t.contig_uid IN (SELECT contig_uid FROM s) AND len(t.nr50_tokens) >= 4
              AND t.nr50_tokens[2] <> t.nr50_tokens[3]
            ORDER BY t.contig_uid LIMIT 1""").fetchone()
        if pick is None:
            say("   q1: no shared contig with 4 or more slots; not checked")
            failed.append("q1 pick")
        else:
            cu, ta, tb = pick
            # independent: every contig with token a at slot i and b at a slot within 1 of i (circular contigs wrap)
            cand = con.execute(f"""
                SELECT t.contig_uid, t.nr50_tokens, ci.is_circular FROM read_parquet({lst(tt)}, hive_partitioning = false) t
                LEFT JOIN (SELECT contig_uid, bool_or(is_circular) AS is_circular FROM read_parquet({lst([p for ps in gi.values() for p in ps])},
                           hive_partitioning = false) GROUP BY contig_uid) ci ON ci.contig_uid = t.contig_uid
                WHERE list_contains(t.nr50_tokens, {ta}::UBIGINT) AND list_contains(t.nr50_tokens, {tb}::UBIGINT)""").fetchall()
            loci, contigs = set(), set()
            for c, toks, circ in cand:
                n = len(toks)
                for i, x in enumerate(toks):
                    if x != ta:
                        continue
                    for j, y in enumerate(toks):
                        if y != tb:
                            continue
                        d = abs(i - j)
                        if circ:
                            d = min(d, n - d)
                        if d == 1:
                            loci.add((c, i))
                            contigs.add(c)
            con.execute("CREATE TEMP TABLE qc (contig_uid BLOB)")
            con.executemany("INSERT INTO qc VALUES (?)", [(c,) for c in sorted(contigs)])
            g_ind_q = con.execute("SELECT count(DISTINCT genome_uid) FROM gi WHERE contig_uid IN "
                                  "(SELECT contig_uid FROM qc)").fetchone()[0]
            q1c = cat.execute(f"SELECT q1_count({ta}::UBIGINT, {tb}::UBIGINT, 1, NULL)").fetchone()[0]
            q1g = cat.execute(f"SELECT q1_genomes({ta}::UBIGINT, {tb}::UBIGINT, 1)").fetchone()[0]
            say(f"   q1 on shared contig {cu.hex()}, NR50 tokens {ta} and {tb} on slots 1 and 2: q1_count {q1c}, "
                f"independent loci {len(loci)} on {len(contigs)} contigs; q1_genomes {q1g}, independent genomes {g_ind_q}")
            d27.update({"q1_contig": cu.hex(), "q1_tokens": [ta, tb], "q1_count": q1c, "loci_independent": len(loci),
                        "contigs": len(contigs), "q1_genomes": q1g, "genomes_independent_q1": g_ind_q})
            if q1c != len(loci) or q1g != g_ind_q or g_ind_q < len(contigs):
                failed.append("q1")
        summary["dec27"] = d27
    if a.est and a.est.exists():
        summary["estimate"] = json.loads(a.est.read_text()).get("projection", {}).get("cluster_cpu_min")
    summary["failed"] = failed
    (a.out / "summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True, default=str))
    say(f"checks: {'all pass' if not failed else 'FAILED: ' + ', '.join(failed)}")
    (a.out / "report.txt").write_text("\n".join(report) + "\n")
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
