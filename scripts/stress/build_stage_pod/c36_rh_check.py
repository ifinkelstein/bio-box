"""Card C-36 (6): the rehearsal's files against the sample's own inputs and the stored stages; read only.

    python c36_rh_check.py --root VIEW --data-root ROOT --stage S --inputs IN --sample SAMPLE --contigs CONTIGS_TSV --logs LOGS
                           --stored S1,P1,... --out OUT

The oracle is independent of the steps checked: the sample's records (CONTIGS_TSV, `c36_verify.py`: each record's
`ids.contig_uid` and length, filtered to the sample's accessions in SAMPLE/lock.json), its genomes (IN/genomes.json,
`make_inputs.py prok`) and the stored stages' `contig_index` in the data root, read by name. Expected:

- ingest stores exactly the sample's distinct contigs that no stored stage holds as a default variant (DEC-26);
- gene calling skips exactly (a) the sample's contigs a stored stage owns, with that stage, and (b) the contigs listed in
  several of the sample's genomes, in every genome but the one of least genome_uid (DEC-27), with that genome; and
  writes no feature_protein row on a contig of (a);
- counts by locus do not fan out: feature_protein holds one row per feature_uid.

Also measured, for the note's estimates: genes and distinct proteins (all, and those on the contigs of each species'
first genome in the stage's order), residues per distinct protein, the share of distinct proteins already in a
stored stage's `protein_seq`, and the run records' CPU, wall time and peak memory. Writes OUT/check.json and prints it;
exit 1 when a check fails.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import duckdb
import pyarrow.parquet as pq

from biobox import ids, io
from biobox.l0seq import store
from biobox.l1 import genecall


def tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def step_result(log: Path) -> dict:
    """The step's printed result: the line `{"result": ...}` its command wrote to the unit's log (pod_run.sh)."""
    for ln in log.read_text().splitlines():
        if ln.startswith('{"result"'):
            return json.loads(ln)["result"]
    return {}


def lit(files) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def main() -> int:
    ap = argparse.ArgumentParser()
    for x in ("--root", "--data-root", "--inputs", "--sample", "--contigs", "--logs", "--out"):
        ap.add_argument(x, required=True, type=Path)
    ap.add_argument("--stage", required=True)
    ap.add_argument("--stored", required=True)
    ap.add_argument("--pl1", type=Path, help="PL1 inputs: count the sample contigs its genomes list")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    stored = a.stored.split(",")
    S = a.stage
    problems: list[str] = []
    out: dict = {"stage": S}

    lock = json.loads((a.sample / "lock.json").read_text())["genomes"]
    sel = tsv(a.sample / "selection.tsv")
    recs = [r for r in tsv(a.contigs) if r["accession"] in lock]
    length = {bytes.fromhex(r["contig_uid"]): int(r["length"]) for r in recs}
    genomes = json.loads((a.inputs / "genomes.json").read_text())
    acc_contigs: dict[str, list[bytes]] = defaultdict(list)
    for r in recs:
        acc_contigs[r["accession"]].append(bytes.fromhex(r["contig_uid"]))

    owner: dict[bytes, str] = {}
    for s in stored:
        for f in io.part_files(a.data_root / store.INDEX_DIR / f"stage={s}"):
            t = pq.read_table(f, columns=["contig_uid", "is_default_variant"])
            for cu, d in zip(t.column("contig_uid").to_pylist(), t.column("is_default_variant").to_pylist()):
                if d and cu in length:
                    owner[cu] = s
    contigs = set(length)
    want_new = {c for c in contigs if c not in owner}
    members: dict[bytes, list[bytes]] = defaultdict(list)
    for g in genomes:
        cus = [bytes.fromhex(u) for u in g["contig_uids"]]
        guid = ids.genome_uid(cus)
        for c in cus:
            members[c].append(guid)
    within = {c: min(gs) for c, gs in members.items() if len(gs) > 1 and c not in owner}
    pl1 = set() if a.pl1 is None else {bytes.fromhex(u) for g in json.loads((a.pl1 / "genomes.json").read_text())
                                       for u in g["contig_uids"]}
    out["sample"] = {"genomes_in_lock": len(lock), "genomes_listed": len(genomes), "records": len(recs),
                     "distinct_contigs": len(contigs), "bases": sum(length.values()),
                     "owned_by_stored": len(owner), "owned_by_stored_bases": sum(length[c] for c in owner),
                     "owned_by_stored_by_stage": {s: sum(1 for v in owner.values() if v == s)
                                                  for s in stored if s in owner.values()},
                     "in_several_genomes_not_owned": len(within), "expected_new_contigs": len(want_new),
                     **({"in_pl1_inputs": len(contigs & pl1), "in_pl1_inputs_bases": sum(length[c] for c in contigs & pl1)}
                        if a.pl1 is not None else {}),
                     "expected_new_bases": sum(length[c] for c in want_new)}

    # ingest: the stage's contig_index holds exactly the new contigs
    ci = io.part_files(a.root / store.INDEX_DIR / f"stage={S}")
    t = pq.read_table(ci, columns=["contig_uid", "is_default_variant", "contig_len"]) if ci else None
    got = {cu for cu, d in zip(t.column("contig_uid").to_pylist(), t.column("is_default_variant").to_pylist()) if d}
    if got != want_new:
        problems.append(f"ingest stored {len(got)} contigs, expected {len(want_new)} "
                        f"({len(got - want_new)} unexpected, {len(want_new - got)} missing)")
    runs = {}
    for step in ("inputs", "ingest", "genecall"):
        rec = json.loads((a.root / "runs" / f"stage={S}" / f"{step}.json").read_text())
        runs[step] = {k: rec.get(k) for k in ("wall_s", "cpu_min", "peak_rss_bytes", "threads", "exit_code",
                                              "image_digest", "over_cpu_budget")}
        res = step_result(a.logs / f"c36-rh-{step}.log")
        runs[step]["result"] = {k: v for k, v in res.items() if k != "skipped_contigs"}
        if step == "genecall":
            skipped = res.get("skipped_contigs", [])
    out["runs"] = runs

    # gene calling: the skipped contigs, with their owners
    want_skip = {c.hex(): (s, None) for c, s in owner.items()}
    want_skip.update({c.hex(): (S, g.hex()) for c, g in within.items()})
    got_skip = {e["contig_uid"]: (e["owner_stage"], e.get("owner_genome")) for e in skipped}
    out["genecall_skipped"] = {"got": len(got_skip), "expected": len(want_skip),
                               "owned_by_stored": sum(1 for v in got_skip.values() if v[0] != S),
                               "owned_by_a_genome_of_the_stage": sum(1 for v in got_skip.values() if v[0] == S),
                               "equal": got_skip == want_skip}
    if got_skip != want_skip:
        diff = sorted(set(got_skip.items()) ^ set(want_skip.items()))[:5]
        problems.append(f"skipped contigs differ from the expected: {diff}")

    # rows: none on a stored stage's contig; one row per feature
    gcv = genecall.GENE_CALL_VERSION
    fp = io.part_files(a.root / "l1/feature_protein" / f"gene_call_version={gcv}" / f"stage={S}")
    ps = io.part_files(a.root / "l0seq/protein_seq" / f"gene_call_version={gcv}" / f"stage={S}", "**/part-*.parquet")
    con = duckdb.connect(config={"threads": 2, "memory_limit": "6GB", "temp_directory": str(a.out / "duckdb")})
    con.execute("CREATE TABLE owned (contig_uid BLOB)")
    if owner:  # executemany refuses an empty list
        con.executemany("INSERT INTO owned VALUES (?)", [(c,) for c in sorted(owner)])
    con.execute(f"CREATE VIEW fp AS SELECT * FROM read_parquet({lit(fp)}, hive_partitioning = false)")
    (n_rows, n_feat, n_on_owned), = con.execute(
        "SELECT count(*), count(DISTINCT feature_uid), count(*) FILTER (WHERE contig_uid IN (SELECT contig_uid FROM "
        "owned)) FROM fp").fetchall()
    if n_on_owned:
        problems.append(f"{n_on_owned} feature_protein rows on contigs a stored stage owns")
    if n_rows != n_feat:
        problems.append(f"feature_protein: {n_rows} rows, {n_feat} distinct features")
    if n_rows != runs["genecall"]["result"].get("genes"):
        problems.append(f"feature_protein rows {n_rows} != the run's genes {runs['genecall']['result'].get('genes')}")

    # the ratios for the estimates: genes and distinct proteins on the contigs of each species' first genome
    first_accs, seen_sp = [], set()
    for r in sel:
        if r["gtdb_taxonomy"] not in seen_sp:
            seen_sp.add(r["gtdb_taxonomy"])
            first_accs.append(r["accession"])
    first_contigs = {c for acc in first_accs for c in acc_contigs[acc]}
    con.execute("CREATE TABLE first (contig_uid BLOB)")
    con.executemany("INSERT INTO first VALUES (?)", [(c,) for c in sorted(first_contigs)])
    (g_all, d_all, g1, d1), = con.execute(
        "SELECT count(*), count(DISTINCT protein_uid), count(*) FILTER (WHERE contig_uid IN (SELECT contig_uid FROM "
        "first)), count(DISTINCT protein_uid) FILTER (WHERE contig_uid IN (SELECT contig_uid FROM first)) FROM fp"
    ).fetchall()
    (n_ps, res), = con.execute(f"SELECT count(*), sum(length) FROM read_parquet({lit(ps)}, hive_partitioning = false)"
                               ).fetchall()
    if n_ps != d_all:
        problems.append(f"protein_seq rows {n_ps} != distinct proteins of feature_protein {d_all}")
    # distinct proteins already in a stored stage's protein_seq (exact identity; what DEC-35 (2) would copy)
    con.execute(f"CREATE TABLE mine AS SELECT DISTINCT protein_uid FROM read_parquet({lit(ps)}, "
                "hive_partitioning = false)")
    old = [f for s in stored for f in io.part_files(
        a.data_root / "l0seq/protein_seq" / f"gene_call_version={gcv}" / f"stage={s}", "**/part-*.parquet")]
    (in_old,), = con.execute(f"SELECT count(*) FROM mine WHERE protein_uid IN (SELECT protein_uid FROM "
                             f"read_parquet({lit(old)}, hive_partitioning = false))").fetchall()
    bases_called = sum(length[c] for c in want_new)
    out["proteins"] = {"genes": g_all, "distinct": d_all, "genes_on_first_genomes": g1, "distinct_on_first_genomes": d1,
                       "first_genomes": len(first_accs), "residues": int(res or 0),
                       "distinct_in_a_stored_stage": in_old, "bases_with_rows": bases_called}
    out["ratios"] = {
        "distinct_per_gene_first": round(d1 / g1, 4) if g1 else None,
        "new_per_gene_later": round((d_all - d1) / (g_all - g1), 4) if g_all > g1 else None,
        "residues_per_distinct": round(int(res) / d_all, 2) if d_all else None,
        "genes_per_mbp": round(g_all / bases_called * 1e6, 1) if bases_called else None,
        "share_distinct_in_a_stored_stage": round(in_old / d_all, 4) if d_all else None}
    out["problems"] = problems
    (a.out / "check.json").write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    print(json.dumps(out, sort_keys=True))
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
