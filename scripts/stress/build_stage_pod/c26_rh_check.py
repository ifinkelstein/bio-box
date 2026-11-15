"""Card C-26 part 1: DEC-26 and DEC-27 on the rehearsal of the PLSDB stage (read only).

    python c26_rh_check.py --root RH_ROOT --data-root ROOT --inputs RH_INPUTS --stage PLR --release rC26.2
        --linked P1,P2 --log UNIT_LOG --out OUT

RH_ROOT is the rehearsal root (c26_view.py: P1 and P2 linked from ROOT); RH_INPUTS the sample's inputs
(`make_inputs.py plasmid ... --root RH_ROOT`: genomes.json, shared_records.tsv). The oracle of every check is the
sample's own inputs, read here, and the linked stages' files in ROOT, never the code under test:

1. DEC-26. The contigs of the sample that P1 or P2 store (shared_records.tsv) are not in the stage's contig_index as
   a default variant; the gene-call ledger lists exactly them with `owned_by:<that stage>`; no file the stage wrote
   (every Parquet file under a `stage=<STAGE>` directory of RH_ROOT with a `contig_uid` column) has a row for one; the
   linked stage's feature_protein holds their genes.
2. DEC-27. The stage's genome_input lists every genome of the sample, the shared ones too (their second genomes);
   contig_dim of the scratch release over the linked stages and the stage has one row per contig (no contig twice),
   `n_genomes` = the number of genomes the union of the stages' genome_input gives the contig, and its genome the
   least of them; genome_dim one row per genome of that union, the sample's genomes of shared contigs among them with
   kind plasmid and one contig.
3. Counts. By genome (distinct genome_uid over the union) the shared contigs count their second genomes; by locus
   (distinct feature_uid of the owning stage on those contigs) they do not fan out, while a join of the features to
   the membership rows counts each gene once per genome (the fan-out a count by locus must avoid, card C-21c).

Writes OUT/check.json and prints its last line; exit 1 when a check fails.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

import duckdb
import pyarrow.parquet as pq

from biobox import io
from biobox.l0seq import store


def sql_list(files) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def main() -> int:
    ap = argparse.ArgumentParser()
    for k in ("--root", "--data-root", "--inputs", "--out", "--log"):
        ap.add_argument(k, required=True, type=Path)
    ap.add_argument("--stage", required=True)
    ap.add_argument("--release", required=True)
    ap.add_argument("--linked", required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    linked = a.linked.split(",")
    stages = [*linked, a.stage]
    res: dict = {}
    bad: list[str] = []

    def need(ok: bool, what: str) -> None:
        if not ok:
            bad.append(what)

    genomes = json.loads((a.inputs / "genomes.json").read_text())
    sample = {bytes.fromhex(g["contig_uids"][0]) for g in genomes}
    with open(a.inputs / "shared_records.tsv", newline="") as fh:
        shared_rows = list(csv.DictReader(fh, delimiter="\t"))
    owned = {bytes.fromhex(r["contig_uid"]): r["owner_stage"] for r in shared_rows}
    res["sample"] = {"genomes": len(genomes), "contigs": len(sample), "shared_records": len(shared_rows),
                     "shared_contigs": len(owned), "shared_by_stage": {s: sum(1 for v in owned.values() if v == s)
                                                                      for s in linked},
                     "circular": sum(1 for g in genomes if g["circular"])}
    need(set(owned.values()) <= set(linked), "a shared contig is owned by a stage that is not linked")

    con = duckdb.connect(config={"threads": 2, "memory_limit": "4GB", "temp_directory": str(a.out / "duckdb")})

    def files(rel: str, s: str, pattern: str = "part-*.parquet") -> list[Path]:
        return io.part_files(a.root / rel / f"stage={s}", pattern)

    # 1. DEC-26 --------------------------------------------------------------------------------------------------
    ci = files("l0seq/contig_index", a.stage)
    own_stage = {r[0] for r in con.execute(f"SELECT contig_uid FROM read_parquet({sql_list(ci)}, hive_partitioning = "
                                           "false) WHERE is_default_variant").fetchall()}
    res["contig_index"] = {"default_variants": len(own_stage), "of_them_shared": len(own_stage & set(owned))}
    need(own_stage == sample - set(owned), "the stage's contig_index is not the sample less the shared contigs")
    led = pq.read_table(a.root / f"ledger/stage={a.stage}/step=genecall/part-0.parquet").to_pylist()
    skipped = {r["object_key"]: r["reason_code"].removeprefix(store.OWNED_BY) for r in led
               if r["reason_code"].startswith(store.OWNED_BY)}
    res["genecall_ledger_skipped"] = {s: sum(1 for v in skipped.values() if v == s) for s in sorted(set(skipped.values()))}
    need(skipped == owned, "the gene-call ledger's owned_by rows are not exactly the shared contigs and their stages")
    # the gene-call result line in the unit log (run-step prints the step's JSON result)
    for line in a.log.read_text().splitlines():
        if not (line.startswith("{") and '"skipped_contigs"' in line):
            continue
        try:
            r = json.loads(line)
        except ValueError:
            continue
        sk = {bytes.fromhex(x["contig_uid"]): x["owner_stage"] for x in r.get("skipped_contigs", [])}
        what = "genecall" if "genes" in r else "other"  # the CRISPR producers report their skips the same way
        res.setdefault("results_with_skips", []).append(
            {"producer": what, "skipped": len(sk), **{k: r[k] for k in r if k != "skipped_contigs"}})
        need(sk == owned, f"a {what} result's skipped_contigs are not the shared contigs and their stages")
    need(any(x["producer"] == "genecall" for x in res.get("results_with_skips", [])), "no gene-call result in the log")
    written: dict[str, dict] = {}
    for dp, dns, fns in os.walk(a.root):  # real directories only: os.walk does not follow the links to ROOT
        if f"stage={a.stage}" not in Path(dp).parts:
            continue
        for fn in fns:
            if not fn.endswith(".parquet"):
                continue
            p = Path(dp) / fn
            if "contig_uid" not in pq.read_schema(p).names:
                continue
            rel = str(p.relative_to(a.root))
            t = pq.read_table(p, columns=["contig_uid"]).column("contig_uid").to_pylist()
            written[rel] = {"rows": len(t), "rows_on_shared": sum(1 for c in t if c in owned)}
    res["stage_files_with_contig_uid"] = written
    need(written != {}, "no stage file with a contig_uid column was found")
    # ingest records every record (contig_occurrence, contig_label) and inputs every genome (genome_input); the
    # stage's contig_index is checked above for default variants
    allowed = ("contig_occurrence", "contig_label", "genome_input", "contig_index")
    need(all(v["rows_on_shared"] == 0 or any(t in k for t in allowed) for k, v in written.items()),
         "a feature, gene, IGR, tape or array file of the stage has a row on a shared contig")
    gcv = next(iter(sorted(os.listdir(a.data_root / "l1" / "feature_protein")))).removeprefix("gene_call_version=")
    con.execute("CREATE TABLE sh (contig_uid BLOB, owner VARCHAR)")
    con.executemany("INSERT INTO sh VALUES (?, ?)", sorted(owned.items()))
    genes_owner = {}
    for s in linked:
        fp = io.part_files(a.data_root / "l1/feature_protein" / f"gene_call_version={gcv}" / f"stage={s}")
        genes_owner[s] = con.execute(
            f"SELECT count(*), count(DISTINCT contig_uid) FROM read_parquet({sql_list(fp)}, hive_partitioning = false) "
            f"WHERE contig_uid IN (SELECT contig_uid FROM sh WHERE owner = '{s}')").fetchall()[0]
    res["owner_genes_on_shared"] = {s: {"genes": g, "contigs_with_genes": c} for s, (g, c) in genes_owner.items()}

    # 2. DEC-27 --------------------------------------------------------------------------------------------------
    gi = {s: files("l0seq/genome_input", s) for s in stages}
    union = " UNION ALL ".join(f"SELECT '{s}' AS stage, contig_uid, genome_uid, kind FROM read_parquet("
                               f"{sql_list(fs)}, hive_partitioning = false)" for s, fs in gi.items())
    con.execute(f"CREATE TABLE memb AS {union}")
    (n_stage_rows, n_stage_contigs), = con.execute(
        f"SELECT count(*), count(DISTINCT contig_uid) FROM memb WHERE stage = '{a.stage}'").fetchall()
    res["stage_genome_input"] = {"rows": n_stage_rows, "contigs": n_stage_contigs}
    need(n_stage_rows == len(genomes) and n_stage_contigs == len(sample),
         "the stage's genome_input does not list every genome of the sample")
    cd = io.part_files(a.root / "dim/contig_dim" / f"release={a.release}")
    gd = io.part_files(a.root / "dim/genome_dim" / f"release={a.release}")
    con.execute(f"CREATE TABLE cdim AS SELECT * FROM read_parquet({sql_list(cd)}, hive_partitioning = false)")
    con.execute(f"CREATE TABLE gdim AS SELECT * FROM read_parquet({sql_list(gd)}, hive_partitioning = false)")
    (rows, distinct), = con.execute("SELECT count(*), count(DISTINCT contig_uid) FROM cdim").fetchall()
    (want_contigs,), = con.execute("SELECT count(DISTINCT contig_uid) FROM memb").fetchall()
    res["contig_dim"] = {"rows": rows, "distinct_contigs": distinct, "membership_contigs": want_contigs}
    need(rows == distinct == want_contigs, "contig_dim is not one row per contig of the release's membership")
    (n_bad,), = con.execute(
        "SELECT count(*) FROM cdim c JOIN (SELECT contig_uid, count(DISTINCT genome_uid) AS n, min(genome_uid) AS g "
        "FROM memb GROUP BY contig_uid) m USING (contig_uid) WHERE c.n_genomes <> m.n OR c.genome_uid <> m.g").fetchall()
    res["contig_dim_disagreeing_with_membership"] = n_bad
    need(n_bad == 0, "contig_dim's n_genomes or genome_uid disagree with the membership")
    (sh_rows, sh_n2), = con.execute(
        "SELECT count(*), count(*) FILTER (WHERE n_genomes >= 2) FROM cdim WHERE contig_uid IN (SELECT contig_uid "
        "FROM sh)").fetchall()
    res["contig_dim_shared"] = {"rows": sh_rows, "with_two_or_more_genomes": sh_n2}
    need(sh_rows == len(owned) == sh_n2, "a shared contig does not have its second genome in contig_dim")
    (g_rows, g_distinct), = con.execute("SELECT count(*), count(DISTINCT genome_uid) FROM gdim").fetchall()
    (want_g,), = con.execute("SELECT count(DISTINCT genome_uid) FROM memb").fetchall()
    res["genome_dim"] = {"rows": g_rows, "distinct": g_distinct, "membership_genomes": want_g}
    need(g_rows == g_distinct == want_g, "genome_dim is not one row per genome of the membership")
    (second, second_ok), = con.execute(
        f"SELECT count(*), count(*) FILTER (WHERE g.kind = 'plasmid' AND g.n_contigs = 1) FROM gdim g WHERE "
        f"genome_uid IN (SELECT genome_uid FROM memb WHERE stage = '{a.stage}' AND contig_uid IN "
        "(SELECT contig_uid FROM sh))").fetchall()
    res["second_genomes_in_genome_dim"] = {"rows": second, "plasmid_with_one_contig": second_ok}
    need(second == second_ok == len(owned), "the sample's genomes of shared contigs are not all in genome_dim")

    # 3. counts by genome and by locus on the shared contigs --------------------------------------------------------
    fps = [f for s in linked for f in io.part_files(a.data_root / "l1/feature_protein" / f"gene_call_version={gcv}"
                                                     / f"stage={s}")]
    con.execute(f"CREATE TABLE gene AS SELECT contig_uid, feature_uid FROM read_parquet({sql_list(fps)}, "
                "hive_partitioning = false) WHERE contig_uid IN (SELECT contig_uid FROM sh)")
    (by_locus, joined), = con.execute(
        "SELECT (SELECT count(DISTINCT feature_uid) FROM gene), (SELECT count(*) FROM gene JOIN (SELECT DISTINCT "
        "contig_uid, genome_uid FROM memb) USING (contig_uid))").fetchall()
    (genomes_on_shared, genomes_owner_only), = con.execute(
        "SELECT count(DISTINCT genome_uid), count(DISTINCT genome_uid) FILTER (WHERE stage <> "
        f"'{a.stage}') FROM memb WHERE contig_uid IN (SELECT contig_uid FROM sh)").fetchall()
    (genes_rows,), = con.execute("SELECT count(*) FROM gene").fetchall()
    res["counts_on_shared"] = {"genes_rows": genes_rows, "by_locus_distinct_features": by_locus,
                               "features_joined_to_membership": joined, "genomes": genomes_on_shared,
                               "genomes_without_the_stage": genomes_owner_only}
    need(by_locus == genes_rows, "a feature of a shared contig is listed twice")
    need(genomes_on_shared == genomes_owner_only + len(owned), "counts by genome do not include the second genomes")
    res["problems"] = bad
    (a.out / "check.json").write_text(json.dumps(res, indent=1, sort_keys=True, default=str))
    print(json.dumps(res, sort_keys=True, default=str))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
