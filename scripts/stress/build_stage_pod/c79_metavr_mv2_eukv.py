"""Card C-79 (DEC-79): stage MetaVR's `>= 30 kb` spillover as `MV2` and its eukaryotic-virus slice as `MVEUKV`.

    python c79_metavr_mv2_eukv.py profile PARQUET LEFT_GZ OUT_JSON
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import duckdb

PROBE_PATTERNS = ("Nucleocytoviricota", "Mirusviricota", "Preplasmiviricota", "Maveriviricetes", "Polintoviricetes",
                  "Lavidaviridae", "Tectiliviricetes", "virophage", "Polinton", "Megaviricetes", "Pokkesviricetes")


def _sql_lit(path: Path | str) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def profile(parquet: Path, left_gz: Path, *, threads: int = 2, memory_limit: str = "4GB",
            temp: Path | None = None) -> dict:
    """The table's own schema and its taxonomy column's values (read before the selection rule is fixed)."""
    cfg = {"threads": int(threads), "memory_limit": memory_limit}
    if temp is not None:
        Path(temp).mkdir(parents=True, exist_ok=True)
        cfg["temp_directory"] = str(temp)
    con = duckdb.connect(config=cfg)
    try:
        pq = _sql_lit(parquet)
        schema = [(r[0], r[1]) for r in con.execute(f"DESCRIBE SELECT * FROM read_parquet({pq})").fetchall()]
        out: dict = {"schema": schema}
        cols = {c for c, _ in schema}
        out["rows"] = int(con.execute(f"SELECT count(*) FROM read_parquet({pq})").fetchone()[0])
        if "ictv_taxonomy_method" in cols:
            out["ictv_taxonomy_method"] = con.execute(
                f"SELECT coalesce(ictv_taxonomy_method, '<null>'), count(*) FROM read_parquet({pq}) "
                "GROUP BY 1 ORDER BY 2 DESC, 1 LIMIT 50").fetchall()
        if "ictv_taxonomy" in cols:
            out["ictv_taxonomy_examples"] = [r[0] for r in con.execute(
                f"SELECT DISTINCT ictv_taxonomy FROM read_parquet({pq}) WHERE ictv_taxonomy IS NOT NULL "
                "ORDER BY 1 LIMIT 20").fetchall()]
            out["ictv_taxonomy_null"] = int(con.execute(
                f"SELECT count(*) FROM read_parquet({pq}) WHERE ictv_taxonomy IS NULL OR ictv_taxonomy = ''"
            ).fetchone()[0])
            ors = " OR ".join(f"ictv_taxonomy ILIKE '%{p}%'" for p in PROBE_PATTERNS)
            con.execute(f"""
                CREATE TEMP TABLE euk AS
                SELECT uvig, votu, quality, length::BIGINT AS length, coalesce(is_concatemer, false) AS is_concatemer,
                       checkv_contamination, ictv_taxonomy
                FROM read_parquet({pq}) WHERE {ors}
            """)
            con.execute(
                f"CREATE TEMP TABLE lft AS SELECT reason, mv2_rank::BIGINT AS mv2_rank, uvig FROM read_csv("
                f"{_sql_lit(left_gz)}, delim='\t', header=true, quote='', all_varchar=true)")
            out["left_rows"] = int(con.execute("SELECT count(*) FROM lft").fetchone()[0])
            out["patterns"] = {}
            for p in PROBE_PATTERNS:
                rows = con.execute(
                    f"SELECT coalesce(l.reason, 'MV1') AS r, count(*), count(DISTINCT e.votu), sum(e.length)::BIGINT "
                    f"FROM euk e LEFT JOIN lft l USING (uvig) WHERE e.ictv_taxonomy ILIKE '%{p}%' GROUP BY 1 ORDER BY 1"
                ).fetchall()
                out["patterns"][p] = {r[0]: {"uvigs": int(r[1]), "votus": int(r[2]), "bp": int(r[3])} for r in rows}
            out["lineages"] = [
                (r[0], int(r[1]), int(r[2]))
                for r in con.execute(
                    "SELECT ictv_taxonomy, count(*), sum(length)::BIGINT FROM euk GROUP BY 1 ORDER BY 2 DESC, 1 LIMIT 300"
                ).fetchall()
            ]
            out["quality_by_reason"] = [
                (r[0], r[1], int(r[2]), int(r[3]))
                for r in con.execute(
                    "SELECT coalesce(l.reason, 'MV1'), e.quality, count(*), sum(e.length)::BIGINT "
                    "FROM euk e LEFT JOIN lft l USING (uvig) GROUP BY 1, 2 ORDER BY 1, 2").fetchall()
            ]
            out["length_classes_lq_nd_lt30k"] = [
                (r[0], int(r[1]), int(r[2]))
                for r in con.execute(
                    "SELECT CASE WHEN e.length < 5000 THEN 'lt_5kb' WHEN e.length < 10000 THEN '5_to_10kb' "
                    "ELSE '10_to_30kb' END, count(*), sum(e.length)::BIGINT FROM euk e JOIN lft l USING (uvig) "
                    "WHERE l.reason = 'lq_nd_lt30k' GROUP BY 1 ORDER BY 1").fetchall()
            ]
        return out
    finally:
        con.close()


EUK_TAXA = (("Nucleocytoviricota", "p__Nucleocytoviricota"), ("Mirusviricota", "p__Mirusviricota"),
            ("Maveriviricetes", "c__Maveriviricetes"), ("Polintoviricetes", "c__Polintoviricetes"))


def slice_counts(parquet: Path, left_gz: Path, *, threads: int = 2, memory_limit: str = "4GB",
                 temp: Path | None = None) -> dict:
    """The card's MVEUKV rule counted (one representative per vOTU of EUK_TAXA, by C-62's order, among
    non-concatemer rows of checkv_contamination <= 10, in vOTUs none of whose UViGs MV1 or MV2 took), with
    breakdowns for the decision on its size."""
    cfg = {"threads": int(threads), "memory_limit": memory_limit}
    if temp is not None:
        Path(temp).mkdir(parents=True, exist_ok=True)
        cfg["temp_directory"] = str(temp)
    con = duckdb.connect(config=cfg)
    try:
        pq = _sql_lit(parquet)
        case = " ".join(f"WHEN list_contains(tok, '{t}') THEN '{n}'" for n, t in EUK_TAXA)
        con.execute(f"""
            CREATE TEMP TABLE raw AS
            SELECT uvig, votu, coalesce(quality, '') AS quality, coalesce(viral_confidence, '') AS viral_confidence,
                   coalesce(is_concatemer, false) AS is_concatemer, length::BIGINT AS length,
                   CASE quality WHEN 'Complete' THEN 0 WHEN 'High-quality' THEN 1 WHEN 'Medium-quality' THEN 2
                        WHEN 'Low-quality' THEN 3 WHEN 'Not-determined' THEN 4 ELSE 5 END AS qrank,
                   least(coalesce(completeness, -1.0), 100.0) AS comp_cap,
                   coalesce(checkv_contamination, 0.0) AS contam_eff,
                   coalesce(checkv_contamination, 999.0) AS contam_rep,
                   CASE {case} ELSE NULL END AS taxon
            FROM (SELECT *, string_split(coalesce(ictv_taxonomy, ''), ';') AS tok FROM read_parquet({pq}))
        """)
        con.execute(
            f"CREATE TEMP TABLE lft AS SELECT reason, mv2_rank::BIGINT AS mv2_rank, uvig FROM read_csv("
            f"{_sql_lit(left_gz)}, delim='\t', header=true, quote='', all_varchar=true)")
        out: dict = {}
        out["table_rows"], out["left_rows"], out["left_not_in_table"] = (int(x) for x in con.execute(
            "SELECT (SELECT count(*) FROM raw), (SELECT count(*) FROM lft), "
            "(SELECT count(*) FROM lft ANTI JOIN raw USING (uvig))").fetchone())
        out["mv2"] = dict(zip(("uvigs", "bp", "min_rank", "max_rank"), (int(x) for x in con.execute(
            "SELECT count(*), sum(r.length)::BIGINT, min(l.mv2_rank), max(l.mv2_rank) FROM lft l JOIN raw r "
            "USING (uvig) WHERE l.reason = 'mv2_spillover'").fetchone())))
        con.execute("""
            CREATE TEMP TABLE touched AS
            SELECT DISTINCT r.votu FROM raw r LEFT JOIN lft l USING (uvig)
            WHERE l.uvig IS NULL OR l.reason = 'mv2_spillover'
        """)
        con.execute("""
            CREATE TEMP TABLE pick AS
            SELECT * FROM (
                SELECT r.*, row_number() OVER (PARTITION BY r.votu
                    ORDER BY qrank, comp_cap DESC, contam_rep, length DESC, uvig) AS rn
                FROM raw r ANTI JOIN touched t USING (votu)
                WHERE r.taxon IS NOT NULL AND NOT r.is_concatemer AND r.contam_eff <= 10.0
            ) WHERE rn = 1
        """)
        q = lambda sql: [tuple(int(x) if isinstance(x, (int,)) else x for x in r) for r in con.execute(sql).fetchall()]
        out["by_taxon"] = q("SELECT taxon, count(*), sum(length)::BIGINT FROM pick GROUP BY 1 ORDER BY 1")
        out["total"] = q("SELECT count(*), sum(length)::BIGINT FROM pick")
        out["by_quality"] = q("SELECT taxon, quality, count(*), sum(length)::BIGINT FROM pick GROUP BY 1, 2 ORDER BY 1, 2")
        out["by_confidence"] = q(
            "SELECT taxon, viral_confidence, count(*), sum(length)::BIGINT FROM pick GROUP BY 1, 2 ORDER BY 1, 2")
        out["by_length"] = q(
            "SELECT taxon, CASE WHEN length < 5000 THEN 'a_lt5kb' WHEN length < 10000 THEN 'b_5_10kb' "
            "WHEN length < 20000 THEN 'c_10_20kb' WHEN length < 30000 THEN 'd_20_30kb' ELSE 'e_ge30kb' END, "
            "count(*), sum(length)::BIGINT FROM pick GROUP BY 1, 2 ORDER BY 1, 2")
        out["floors"] = {
            str(f): q(f"SELECT count(*), coalesce(sum(length), 0)::BIGINT FROM pick WHERE length >= {f}")[0]
            for f in (5000, 10000, 15000, 20000)}
        out["high_confidence"] = q("SELECT count(*), coalesce(sum(length), 0)::BIGINT FROM pick "
                                   "WHERE viral_confidence = 'High'")[0]
        out["high_confidence_ge10kb"] = q("SELECT count(*), coalesce(sum(length), 0)::BIGINT FROM pick "
                                          "WHERE viral_confidence = 'High' AND length >= 10000")[0]
        out["euk_rows_any"] = q("SELECT taxon, count(*), count(DISTINCT votu), sum(length)::BIGINT FROM raw "
                                "WHERE taxon IS NOT NULL GROUP BY 1 ORDER BY 1")
        out["euk_votus_untouched"] = q(
            "SELECT taxon, count(DISTINCT votu) FROM raw ANTI JOIN touched USING (votu) "
            "WHERE taxon IS NOT NULL GROUP BY 1 ORDER BY 1")
        out["mirus_rows"] = q("SELECT count(*) FROM read_parquet(" + pq + ") WHERE ictv_taxonomy ILIKE '%mirus%' "
                              "OR ictv_taxonomy_method ILIKE '%mirus%'")[0]
        return out
    finally:
        con.close()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Stage MetaVR MV2 and MVEUKV inputs (card C-79, DEC-79)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("profile")
    p.add_argument("parquet", type=Path)
    p.add_argument("left_gz", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--memory-limit", default="4GB")
    p.add_argument("--temp", type=Path, default=None)
    s = sub.add_parser("slice")
    s.add_argument("parquet", type=Path)
    s.add_argument("left_gz", type=Path)
    s.add_argument("out", type=Path)
    s.add_argument("--threads", type=int, default=2)
    s.add_argument("--memory-limit", default="4GB")
    s.add_argument("--temp", type=Path, default=None)
    a = ap.parse_args(argv)
    if a.cmd in ("profile", "slice"):
        fn = profile if a.cmd == "profile" else slice_counts
        rep = fn(a.parquet, a.left_gz, threads=a.threads, memory_limit=a.memory_limit, temp=a.temp)
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(json.dumps(rep, indent=1, sort_keys=True) + "\n")
        print(json.dumps({"rows": rep.get("rows"), "out": str(a.out)}))
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
