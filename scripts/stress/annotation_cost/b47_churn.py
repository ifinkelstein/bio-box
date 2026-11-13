"""Card B-47 (design): what representatives-only Pfam does across releases and within one stage. Reads only;
writes one JSON file.

    python b47_churn.py <root> <out.json> [--memory 5GB] [--tmp DIR]

Releases read: rN.1a (S1), rN.1b (S1, P1), rN.3 (S1, P1, P2; its clusters equal rN.1's byte for byte). Every file
is one the release manifests list (never a glob), checked against the manifest's sha256: `protein_seq` (uid and
`length`; card 0-09b found `length` equal to the stripped sequence length on every row of these files),
`cluster_member` and `cluster_member_aln` at level nr90, and the three `l1_domains` files rN.3 lists (S1, P1, P2;
the S1 and P1 files are the ones 0-09b read).

The rule (DEC-32): a protein is scanned when it is its NR90 cluster's representative or its length differs from
its representative's by more than 20% (|L_m - L_r| > 0.2 L_r); every other member is inferred from its
representative. Definitions of right, false and missed pairs are 0-09b's (`representatives.py`): pairs are
(protein, model accession), distinct.

1. `decision`: per release and scope (each stage, all), the scanned residues and the pair counts, as 0-09b's
   decision table; rN.1b reproduces 0-09b's numbers (a check of this script). rN.1a's S1 scope is the stage-local
   case (S1 clustered alone). Also: inferred rows (each inferred member gets one row per domain row of its
   representative) and the same counts for members whose representative is in another stage.
2. `churn`: from release A to release B over A's proteins, NR90 transitions (rep/member, same or other
   representative), and the scan sets: proteins scanned under B that A did not scan (what B must scan if A's
   measured rows are kept), and the reverse; for B's proteins that A lacks, how many B infers from a protein A
   measured. Pairs: rN.1a -> rN.1b and rN.1b -> rN.3.
3. `coverage` (rN.3): each inferred row's representative envelope against the representative's aligned span in
   `cluster_member_aln` (t* are representative coordinates): the fraction of the envelope inside the span, binned,
   and per bin how many inferred rows are true pairs of the member's own scan.

DuckDB, one thread, `--memory` cap, spilling under `--tmp`. Deterministic: no clock or random input in the output.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
from pathlib import Path

import duckdb

RELEASES = ("rN.1a", "rN.1b", "rN.3")
STAGES = ("S1", "P1", "P2")
BIT = {s: 1 << i for i, s in enumerate(STAGES)}
LIBRARY = "pfam-38.2"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def q(s: str) -> str:
    return "'" + s.replace("'", "''") + "'"


def plist(paths: list[Path]) -> str:
    return "[" + ", ".join(q(str(p)) for p in paths) + "]"


def stage_of(path: str) -> str:
    return re.search(r"/stage=([^/]+)/", path).group(1)


def files_of(root: Path, manifest: dict, table: str, read: list) -> list[tuple[Path, dict]]:
    out = []
    for f in sorted((f for f in manifest["files"] if f["table"] == table), key=lambda f: f["path"]):
        p = root / f["path"]
        got = sha256(p)
        if got != f["sha256"]:
            raise SystemExit(f"{p}: sha256 {got} is not the manifest's {f['sha256']}")
        read.append({"path": f["path"], "bytes": p.stat().st_size, "sha256": got, "rows": f.get("rows")})
        out.append((p, f))
    return out


def one(con, sql: str):
    return con.execute(sql).fetchone()


def scopes(stages: list[str]) -> list[tuple[str, int]]:
    return [(s, BIT[s]) for s in stages] + [("all", sum(BIT[s] for s in stages))]


def decision(con, rel: str, stages: list[str]) -> dict:
    t = rel.replace(".", "_")
    out = {}
    for name, bits in scopes(stages):
        r = one(con, f"""
            SELECT count(*), sum(h.lm),
                   count(*) FILTER (WHERE h.scanned), coalesce(sum(h.lm) FILTER (WHERE h.scanned), 0),
                   count(*) FILTER (WHERE h.is_rep), count(*) FILTER (WHERE h.mismatch),
                   coalesce(sum(coalesce(dm.n, 0)), 0),
                   coalesce(sum(coalesce(dm.n, 0)) FILTER (WHERE h.scanned), 0),
                   coalesce(sum(coalesce(x.ni, 0)) FILTER (WHERE NOT h.scanned), 0),
                   coalesce(sum(coalesce(dr.n, 0) - coalesce(x.ni, 0)) FILTER (WHERE NOT h.scanned), 0),
                   coalesce(sum(coalesce(dm.n, 0) - coalesce(x.ni, 0)) FILTER (WHERE NOT h.scanned), 0),
                   coalesce(sum(coalesce(rr.n, 0)) FILTER (WHERE NOT h.scanned), 0),
                   count(*) FILTER (WHERE NOT h.scanned AND (h.rmask & {bits}) = 0),
                   coalesce(sum(coalesce(dr.n, 0) - coalesce(x.ni, 0))
                            FILTER (WHERE NOT h.scanned AND (h.rmask & {bits}) = 0), 0),
                   coalesce(sum(coalesce(x.ni, 0)) FILTER (WHERE NOT h.scanned AND (h.rmask & {bits}) = 0), 0),
                   coalesce(sum(coalesce(dm.n, 0) - coalesce(x.ni, 0))
                            FILTER (WHERE NOT h.scanned AND (h.rmask & {bits}) = 0), 0)
            FROM h_{t} h
            LEFT JOIN dn dm ON dm.p = h.m LEFT JOIN dn dr ON dr.p = h.r LEFT JOIN rn rr ON rr.p = h.r
            LEFT JOIN inter_{t} x ON x.m = h.m
            WHERE (h.mm & {bits}) <> 0""")
        (n, res, n_sc, res_sc, n_rep, n_mis, true_, right_sc, right_inf, false_, missed, inf_rows,
         n_other, false_other, right_other, missed_other) = r
        right = right_sc + right_inf
        out[name] = {
            "proteins": n, "residues": res, "representatives": n_rep, "members_scanned_for_length": n_mis,
            "scanned_proteins": n_sc, "scanned_residues": res_sc,
            "scanned_residue_fraction": round(res_sc / res, 6) if res else None,
            "inferred_members": n - n_sc, "inferred_rows": inf_rows,
            "true_pairs": true_, "right": right, "false": false_, "missed": missed,
            "right_of_true": round(right / true_, 6) if true_ else None,
            "false_of_reported": round(false_ / (right + false_), 6) if right + false_ else None,
            "missed_of_true": round(missed / true_, 6) if true_ else None,
            "inferred_with_rep_in_another_stage": {"members": n_other, "right": right_other, "false": false_other,
                                                   "missed": missed_other},
        }
    return out


def churn(con, a: str, b: str) -> dict:
    ta, tb = a.replace(".", "_"), b.replace(".", "_")
    rows = con.execute(f"""
        SELECT x.is_rep AS a_rep, y.is_rep AS b_rep, (x.r = y.r) AS same_rep, x.scanned AS a_scan, y.scanned AS b_scan,
               count(*) AS n, sum(y.lm) AS res
        FROM h_{ta} x JOIN h_{tb} y ON y.m = x.m
        GROUP BY ALL ORDER BY ALL""").fetchall()
    trans: dict[str, int] = {}
    scan: dict[str, list[int]] = {}
    for a_rep, b_rep, same, a_scan, b_scan, n, res in rows:
        k = ("rep" if a_rep else "member") + "->" + ("rep" if b_rep else "member")
        if not a_rep and not b_rep:
            k += ", same rep" if same else ", other rep"
        trans[k] = trans.get(k, 0) + n
        s = f"{'scanned' if a_scan else 'inferred'} in {a} -> {'scanned' if b_scan else 'inferred'} in {b}"
        scan.setdefault(s, [0, 0])
        scan[s][0] += n
        scan[s][1] += res
    new = one(con, f"""
        SELECT count(*), sum(y.lm), count(*) FILTER (WHERE y.scanned), coalesce(sum(y.lm) FILTER (WHERE y.scanned), 0),
               count(*) FILTER (WHERE NOT y.scanned AND xr.scanned),
               count(*) FILTER (WHERE NOT y.scanned AND xr.m IS NOT NULL AND NOT xr.scanned),
               count(*) FILTER (WHERE NOT y.scanned AND xr.m IS NULL)
        FROM h_{tb} y LEFT JOIN h_{ta} x ON x.m = y.m LEFT JOIN h_{ta} xr ON xr.m = y.r
        WHERE x.m IS NULL""")
    scratch_b = one(con, f"SELECT count(*) FILTER (WHERE scanned), sum(lm) FILTER (WHERE scanned), sum(lm) FROM h_{tb}")
    old_new_scans = sum(v[1] for k, v in scan.items() if k.startswith("inferred") and k.endswith(f"scanned in {b}"))
    old_new_n = sum(v[0] for k, v in scan.items() if k.startswith("inferred") and k.endswith(f"scanned in {b}"))
    return {
        "proteins_of_a": sum(trans.values()),
        "nr90_transitions": dict(sorted(trans.items())),
        "scan_sets_over_a_proteins": {k: {"proteins": v[0], "residues": v[1]} for k, v in sorted(scan.items())},
        "new_proteins": {"proteins": new[0], "residues": new[1], "scanned": new[2], "scanned_residues": new[3],
                         "inferred_from_a_protein_a_scanned": new[4],
                         "inferred_from_a_protein_a_inferred": new[5],
                         "inferred_from_a_new_protein": new[6]},
        "b_from_scratch": {"scanned_proteins": scratch_b[0], "scanned_residues": scratch_b[1],
                           "all_residues": scratch_b[2]},
        "b_keeping_a_measured": {"old_proteins_to_scan": old_new_n, "old_residues_to_scan": old_new_scans,
                                 "new_proteins_to_scan": new[2], "new_residues_to_scan": new[3],
                                 "residues_to_scan": old_new_scans + new[3]},
    }


def coverage(con, rel: str) -> dict:
    t = rel.replace(".", "_")
    rows = con.execute(f"""
        WITH inf AS (SELECT m, r FROM h_{t} WHERE NOT scanned),
        cov AS (SELECT i.m, d.acc, d.env_start, d.env_end, a.tstart, a.tend,
                       greatest(0, least(d.env_end, a.tend) - greatest(d.env_start, a.tstart))::DOUBLE
                         / (d.env_end - d.env_start) AS f
                FROM inf i JOIN aln_{t} a ON a.m = i.m JOIN dr d ON d.p = i.r)
        SELECT CASE WHEN f >= 1 THEN '1.0' WHEN f >= 0.8 THEN '[0.8,1)' WHEN f >= 0.5 THEN '[0.5,0.8)'
                    WHEN f > 0 THEN '(0,0.5)' ELSE '0' END AS bin,
               count(*) AS rows, count(z.p) AS true_pair
        FROM cov LEFT JOIN dpair z ON z.p = cov.m AND z.acc = cov.acc GROUP BY bin ORDER BY bin""").fetchall()
    return {b: {"rows": n, "true_pair": t_} for b, n, t_ in rows}


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("out")
    ap.add_argument("--memory", default="5GB")
    ap.add_argument("--tmp", default=None)
    a = ap.parse_args(argv)
    root = Path(a.root)
    t0 = time.monotonic()
    man = {r: json.loads((root / f"manifests/release={r}.json").read_bytes()) for r in RELEASES}
    read: list = []
    cfg = {"threads": 1, "memory_limit": a.memory, "preserve_insertion_order": False}
    db = ":memory:"
    if a.tmp:  # a database file there, so that tables as well as operators spill under the cap
        cfg["temp_directory"] = str(Path(a.tmp) / "spill")
        db = str(Path(a.tmp) / "b47.duckdb")
        Path(db).unlink(missing_ok=True)
    con = duckdb.connect(db, config=cfg)

    # proteins: rN.3's protein_seq files cover S1, P1, P2 (the files of rN.1a and rN.1b are among them)
    prot_files = files_of(root, man["rN.3"], "protein_seq", read)
    for r in ("rN.1a", "rN.1b"):
        theirs = {f["sha256"] for f in man[r]["files"] if f["table"] == "protein_seq"}
        if not theirs <= {f["sha256"] for _, f in prot_files}:
            raise SystemExit(f"{r}: protein_seq files not among rN.3's")
    parts = " UNION ALL ".join(
        f"SELECT protein_uid, length, {BIT[stage_of(f['path'])]} AS mask FROM read_parquet({q(str(p))})"
        for p, f in prot_files)
    con.execute(f"CREATE TABLE prot AS SELECT protein_uid AS uid, max(length)::BIGINT AS len, "
                f"count(DISTINCT length) AS nlen, bit_or(mask) AS mask FROM ({parts}) GROUP BY protein_uid")
    if one(con, "SELECT count(*) FROM prot WHERE nlen <> 1")[0]:
        raise SystemExit("a protein has two lengths")
    doc: dict = {"proteins": dict(zip(("distinct", "in_two_or_more_stages"),
                                      one(con, "SELECT count(*), count(*) FILTER (WHERE bit_count(mask) > 1) FROM prot"))),
                 "rule": "scan representatives and members with |L_m - L_r| > 0.2 L_r; infer the rest (DEC-32)"}

    # domains: rN.3's three l1_domains files (library pfam-38.2, stages S1, P1, P2)
    dom = [p for p, f in files_of(root, man["rN.3"], "l1_domains", read) if f"library={LIBRARY}/" in f["path"]]
    con.execute(f"CREATE TABLE dr AS SELECT DISTINCT protein_uid AS p, model_acc AS acc, env_start, env_end "
                f"FROM read_parquet({plist(dom)})")
    con.execute("CREATE TABLE dpair AS SELECT DISTINCT p, acc FROM dr")
    con.execute("CREATE TABLE dn AS SELECT p, count(*) AS n FROM dpair GROUP BY p")
    con.execute("CREATE TABLE rn AS SELECT p, count(*) AS n FROM dr GROUP BY p")
    doc["domains"] = dict(zip(("distinct_rows", "pairs", "proteins_with_a_domain"),
                              one(con, "SELECT (SELECT count(*) FROM dr), (SELECT count(*) FROM dpair), "
                                       "(SELECT count(*) FROM dn)")))
    doc["domains"]["rows_per_file"] = {stage_of(str(p)): one(con, f"SELECT count(*) FROM read_parquet({q(str(p))})")[0]
                                       for p in dom}

    for rel in RELEASES:
        t = rel.replace(".", "_")
        stages = man[rel]["stages"]
        (cm,) = [p for p, _ in files_of(root, man[rel], "cluster_member", read)]
        con.execute(f"CREATE TABLE cm_{t} AS SELECT member_uid AS m, rep_uid AS r FROM read_parquet({q(str(cm))}) "
                    "WHERE level = 'nr90'")
        mask = sum(BIT[s] for s in stages)
        n_cm, n_dist = one(con, f"SELECT count(*), count(DISTINCT m) FROM cm_{t}")
        n_prot = one(con, f"SELECT count(*) FROM prot WHERE (mask & {mask}) <> 0")[0]
        if not n_cm == n_dist == n_prot:
            raise SystemExit(f"{rel}: {n_cm} nr90 rows, {n_dist} members, {n_prot} proteins in its stages")
        con.execute(f"""CREATE TABLE h_{t} AS SELECT c.m, c.r, pm.len AS lm, pr.len AS lr, pm.mask & {mask} AS mm,
                        pr.mask & {mask} AS rmask, (c.m = c.r) AS is_rep,
                        (c.m <> c.r AND abs(pm.len - pr.len) > 0.2 * pr.len) AS mismatch,
                        (c.m = c.r OR abs(pm.len - pr.len) > 0.2 * pr.len) AS scanned
                        FROM cm_{t} c JOIN prot pm ON pm.uid = c.m JOIN prot pr ON pr.uid = c.r""")
        con.execute(f"DROP TABLE cm_{t}")
        con.execute(f"""CREATE TABLE inter_{t} AS SELECT h.m, count(*) AS ni FROM h_{t} h
                        JOIN dpair a ON a.p = h.m JOIN dpair b ON b.p = h.r AND b.acc = a.acc
                        WHERE NOT h.scanned GROUP BY h.m""")
        doc.setdefault("decision", {})[rel] = decision(con, rel, stages)
        print(f"b47: {rel} decision done at {time.monotonic() - t0:.0f} s", file=sys.stderr, flush=True)

    doc["churn"] = {"rN.1a->rN.1b": churn(con, "rN.1a", "rN.1b"), "rN.1b->rN.3": churn(con, "rN.1b", "rN.3")}
    print(f"b47: churn done at {time.monotonic() - t0:.0f} s", file=sys.stderr, flush=True)

    (aln,) = [p for p, _ in files_of(root, man["rN.3"], "cluster_member_aln", read)]
    con.execute(f"CREATE TABLE aln_rN_3 AS SELECT member_uid AS m, tstart::BIGINT AS tstart, tend::BIGINT AS tend "
                f"FROM read_parquet({q(str(aln))}) WHERE level = 'nr90' AND member_uid <> rep_uid")
    doc["coverage"] = {"rN.3": coverage(con, "rN.3")}
    doc["files_read"] = read
    con.close()
    if db != ":memory:":
        Path(db).unlink(missing_ok=True)
    out = Path(a.out)
    tmp = out.with_name(f".{out.name}.tmp")
    tmp.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
    tmp.replace(out)
    print(f"b47: done in {time.monotonic() - t0:.0f} s", file=sys.stderr, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
