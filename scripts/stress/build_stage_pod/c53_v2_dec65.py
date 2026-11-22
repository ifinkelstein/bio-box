"""Card C-53 (DEC-65): exact build-order contig_uid ownership for re-staged V2 and GV1 after the 11 completed
root stages (S1, P1..P6, H1, PL1, C5, V1) plus G1 and SPX1 (read-only on the root and staged inputs).

    python c53_v2_dec65.py --root ROOT --root-stages S1,P1,P2,P3,P4,P5,P6,H1,PL1,C5,V1 \
        --g1 G1_DIR --spx1 SPX1_DIR --v2 V2_DIR --gv1 GV1_DIR --work WORK_DIR --out OUT_JSON
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb

from biobox import io as bbio
from biobox.l0seq import store


def _lit(files) -> str:
    return "[" + ", ".join("'" + str(f).replace("'", "''") + "'" for f in files) + "]"


def _genome_uids_sql(genomes_json: Path) -> str:
    p = str(genomes_json).replace("'", "''")
    return (f"SELECT DISTINCT from_hex(u) AS uid FROM (SELECT unnest(contig_uids) AS u FROM read_json('{p}', "
            "format = 'array', columns = {'contig_uids': 'VARCHAR[]'}))")


def compute_overlap(root: Path, root_stages: list[str], g1: Path, spx1: Path, v2: Path, gv1: Path,
                    work: Path) -> dict:
    work.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(config={"threads": 2, "memory_limit": "3GB", "temp_directory": str(work / "duck")})

    # Distinct contigs of V2 (as staged) and GV1 (as staged), plus lengths from records.tsv
    v2_rec = str(v2 / "records.tsv").replace("'", "''")
    gv1_rec = str(gv1 / "records.tsv").replace("'", "''")
    rd = "delim = '\t', header = true, quote = '', all_varchar = true"
    con.execute(f"CREATE TABLE v2_all AS SELECT DISTINCT from_hex(contig_uid) AS uid, length::BIGINT AS n "
                f"FROM read_csv('{v2_rec}', {rd})")
    con.execute(f"CREATE TABLE v2_src AS SELECT DISTINCT from_hex(source_contig_uid) AS uid, source_length::BIGINT AS n "
                f"FROM read_csv('{v2_rec}', {rd})")
    con.execute(f"CREATE TABLE v2_imgvr AS SELECT DISTINCT from_hex(contig_uid) AS uid, length::BIGINT AS n "
                f"FROM read_csv('{v2_rec}', {rd}) WHERE source = 'imgvr'")
    con.execute(f"CREATE TABLE v2_inovirus AS SELECT DISTINCT from_hex(contig_uid) AS uid, length::BIGINT AS n "
                f"FROM read_csv('{v2_rec}', {rd}) WHERE source = 'inovirus'")
    con.execute(f"CREATE TABLE gv1_all AS SELECT DISTINCT from_hex(contig_uid) AS uid, length::BIGINT AS n "
                f"FROM read_csv('{gv1_rec}', {rd})")

    # Verify genomes.json matches records.tsv distinct contig sets
    con.execute(f"CREATE TABLE v2_gj AS {_genome_uids_sql(v2 / 'genomes.json')}")
    con.execute(f"CREATE TABLE gv1_gj AS {_genome_uids_sql(gv1 / 'genomes.json')}")
    v2_diff = con.execute("SELECT (SELECT count(*) FROM v2_gj ANTI JOIN v2_all USING (uid)) + "
                          "(SELECT count(*) FROM v2_all ANTI JOIN v2_gj USING (uid))").fetchone()[0]
    gv1_diff = con.execute("SELECT (SELECT count(*) FROM gv1_gj ANTI JOIN gv1_all USING (uid)) + "
                           "(SELECT count(*) FROM gv1_all ANTI JOIN gv1_gj USING (uid))").fetchone()[0]
    if v2_diff or gv1_diff:
        raise RuntimeError(f"genomes.json vs records.tsv mismatch: v2={v2_diff}, gv1={gv1_diff}")

    con.execute("CREATE TABLE targets AS SELECT DISTINCT uid FROM "
                "(SELECT uid FROM v2_all UNION ALL SELECT uid FROM v2_src UNION ALL SELECT uid FROM gv1_all)")

    # Explicit 11 completed root stages only (ignores any in-flight stage=G1 directory in contig_index)
    con.execute("CREATE TABLE root_own (uid BLOB, st VARCHAR)")
    for s in root_stages:
        files = bbio.part_files(root / store.INDEX_DIR / f"stage={s}")
        if not files:
            raise RuntimeError(f"stage {s}: no contig_index part files under {root}")
        lit = "'" + s.replace("'", "''") + "'"
        con.execute(f"INSERT INTO root_own SELECT DISTINCT contig_uid, {lit} FROM read_parquet({_lit(files)}, "
                    "hive_partitioning = false) WHERE is_default_variant AND contig_uid IN (SELECT uid FROM targets)")
    twice = con.execute("SELECT count(*) FROM (SELECT uid FROM root_own GROUP BY uid HAVING count(DISTINCT st) > 1)"
                        ).fetchone()[0]
    if twice:
        raise RuntimeError(f"{twice} target contigs held by >1 root stage")

    con.execute(f"CREATE TABLE o_G1 AS SELECT uid FROM ({_genome_uids_sql(g1 / 'genomes.json')}) "
                "SEMI JOIN targets USING (uid)")
    con.execute(f"CREATE TABLE o_SPX1 AS SELECT uid FROM ({_genome_uids_sql(spx1 / 'genomes.json')}) "
                "SEMI JOIN targets USING (uid)")

    def build_order_stats(tbl: str, include_v2: bool) -> dict:
        contigs, bases = con.execute(f"SELECT count(*), coalesce(sum(n), 0)::BIGINT FROM {tbl}").fetchone()
        # 11 root stages only
        root_only = dict(con.execute(
            f"SELECT r.st, count(*) FROM {tbl} t JOIN root_own r USING (uid) GROUP BY r.st ORDER BY r.st").fetchall())
        root_only_bases = dict(con.execute(
            f"SELECT r.st, sum(t.n)::BIGINT FROM {tbl} t JOIN root_own r USING (uid) GROUP BY r.st ORDER BY r.st"
        ).fetchall())
        # Full build order: 11 root stages -> G1 -> SPX1 (-> V2 for GV1)
        whens = [
            "WHEN r.st IS NOT NULL THEN r.st",
            "WHEN t.uid IN (SELECT uid FROM o_G1) THEN 'G1'",
            "WHEN t.uid IN (SELECT uid FROM o_SPX1) THEN 'SPX1'",
        ]
        if include_v2:
            whens.append("WHEN t.uid IN (SELECT uid FROM v2_all) THEN 'V2'")
        case = "CASE " + " ".join(whens) + " END"
        held_by = dict(con.execute(
            f"SELECT h, count(*) FROM (SELECT {case} AS h FROM {tbl} t LEFT JOIN root_own r USING (uid)) "
            "WHERE h IS NOT NULL GROUP BY h ORDER BY h").fetchall())
        held_by_bases = dict(con.execute(
            f"SELECT h, sum(n)::BIGINT FROM (SELECT {case} AS h, t.n FROM {tbl} t LEFT JOIN root_own r USING (uid)) "
            "WHERE h IS NOT NULL GROUP BY h ORDER BY h").fetchall())
        held_c = sum(held_by.values())
        held_b = sum(held_by_bases.values())
        return {
            "contigs": int(contigs),
            "bases": int(bases),
            "root_11_held_by": {k: int(v) for k, v in root_only.items()},
            "root_11_held_by_bases": {k: int(v) for k, v in root_only_bases.items()},
            "root_11_held_contigs": int(sum(root_only.values())),
            "root_11_held_bases": int(sum(root_only_bases.values())),
            "build_order_held_by": {k: int(v) for k, v in held_by.items()},
            "build_order_held_by_bases": {k: int(v) for k, v in held_by_bases.items()},
            "build_order_held_contigs": int(held_c),
            "build_order_held_bases": int(held_b),
            "adds_contigs": int(contigs - held_c),
            "adds_bases": int(bases - held_b),
        }

    doc = {
        "root_stages": root_stages,
        "build_order": [*root_stages, "G1", "SPX1", "V2", "GV1"],
        "V2": build_order_stats("v2_all", include_v2=False),
        "V2_source_bytes": build_order_stats("v2_src", include_v2=False),
        "V2_imgvr": build_order_stats("v2_imgvr", include_v2=False),
        "V2_inovirus": build_order_stats("v2_inovirus", include_v2=False),
        "GV1": build_order_stats("gv1_all", include_v2=True),
    }
    con.close()
    return doc


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--root-stages", required=True)
    ap.add_argument("--g1", required=True, type=Path)
    ap.add_argument("--spx1", required=True, type=Path)
    ap.add_argument("--v2", required=True, type=Path)
    ap.add_argument("--gv1", required=True, type=Path)
    ap.add_argument("--work", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    doc = compute_overlap(a.root, a.root_stages.split(","), a.g1, a.spx1, a.v2, a.gv1, a.work)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
    print(json.dumps(doc, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
