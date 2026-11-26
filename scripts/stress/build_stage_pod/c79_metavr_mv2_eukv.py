"""Card C-79 (DEC-79): stage MetaVR's `>= 30 kb` spillover as `MV2` and its eukaryotic-virus slice as `MVEUKV`.

    python c79_metavr_mv2_eukv.py profile PARQUET LEFT_GZ OUT_JSON      (the table's schema and taxonomy values)
    python c79_metavr_mv2_eukv.py slice   PARQUET LEFT_GZ OUT_JSON      (the slice's size before its rule was fixed)
    python c79_metavr_mv2_eukv.py stage   PARQUET LEFT_GZ FNA_BGZ META_OUT FASTA_OUT [--root ROOT]
        [--chain-home /tmp/if2453/s2/inputs] [--later GV1] [--inputs-home /tmp/if2453/s2/inputs] [--work DIR]
        [--bgzip] [--report-out F]
    python c79_metavr_mv2_eukv.py ani     META_DIR FASTA GV1_FASTA OUT [--threads 3] [--every 25]

Selection (deterministic, independent of row order; LEFT_GZ is card C-62's `MV1_left_for_later.tsv.gz`, every
row of the table that MV1 did not take):
- `MV2`: the rows of LEFT_GZ with reason `mv2_spillover`, in `mv2_rank` order (refused unless the ranks are 1..N).
- `MVEUKV` (rule decided 2026-11-26 by the orchestrating session after the slice count, card C-79): the vOTUs of
  `ictv_taxonomy` tokens `p__Nucleocytoviricota`, `c__Maveriviricetes` (virophages) and `c__Polintoviricetes`
  (polintons and PLVs) none of whose UViGs MV1 or MV2 took; among their rows of those taxa that are not concatemers,
  have `checkv_contamination <= 10` (NULL counts as 0, as C-62), `viral_confidence = 'High'` and `length >= 5000`,
  one per vOTU by C-62's representative order (quality tier, capped completeness, contamination, length, uvig).
  Ranked by (taxon in that order, then the same key). `c__Tectiliviricetes` (tectiviruses, adenoviruses) is not
  taken. The table's `ictv_taxonomy` names no Mirusviricota (0 rows, profile 2026-11-26).
- Every other vOTU of those taxa that MV1 and MV2 left goes to `MVEUKV_left_for_later.tsv.gz` with its best eligible
  row and reason `lt5kb` or `low_confidence_ge5kb`, or `no_eligible_member` (only concatemers or contamination
  over 10) with its best row (DEC-59).

Staging: one stream of FNA_BGZ writes both FASTA files (dot-partial names, renamed at the end) under
FASTA_OUT/<STAGE>/ and the metadata under META_OUT/<STAGE>/ (written as <STAGE>.tmp, renamed at the end); DEC-65
exact 20..1,000 bp DTR trimming as `c62_metavr_inputs.py`; DEC-27 repeated contig sets in `repeated_records.tsv`;
`shared_records.tsv` the contigs the root's stored stages or the `r26.11.4` stages (in their chain order) hold;
`later_records.tsv` those that only a stage of the `r26.11.5` chain before this one holds (GV1; and MV2 for MVEUKV).
`objects.json` names `<inputs-home>/<STAGE>/<file>` (card A-29 copies there), `objects_staging.json` the file on
/stor. Both stages stay `redistributable: false`.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import duckdb

sys.path.insert(0, str(Path(__file__).resolve().parent))

import a26_r26_11_4_chain as a26  # noqa: E402
import c39_ani  # noqa: E402
import c62_metavr_inputs as c62  # noqa: E402
import make_inputs  # noqa: E402
from biobox import fasta, ids  # noqa: E402

Refused = make_inputs.Refused

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


# the slice as counted on 2026-11-26 (the card's four taxa; Mirusviricota matched no row)
SLICE_TAXA = (("Nucleocytoviricota", "p__Nucleocytoviricota"), ("Mirusviricota", "p__Mirusviricota"),
              ("Maveriviricetes", "c__Maveriviricetes"), ("Polintoviricetes", "c__Polintoviricetes"))
EUK_TAXA = (("Nucleocytoviricota", "p__Nucleocytoviricota"), ("Maveriviricetes", "c__Maveriviricetes"),
            ("Polintoviricetes", "c__Polintoviricetes"))
MIN_LEN = 5000  # DEC-60's floor for metagenome contigs
CONFIDENCE = "High"
MAX_CONTAMINATION = c62.MAX_CONTAMINATION
ROOT_STAGES = ("S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1", "C5", "V1", "G1", "SPX1")
FASTA_NAME = {"MV2": "mv2.fna", "MVEUKV": "mveukv.fna"}
LEFT_EUKV = "MVEUKV_left_for_later.tsv.gz"
RULE = {
    "MV2": "the 98,920 mv2_spillover rows of MV1_left_for_later.tsv.gz (card C-62) in mv2_rank order",
    "MVEUKV": "one representative per MetaVR vOTU of p__Nucleocytoviricota, c__Maveriviricetes, c__Polintoviricetes "
              "untouched by MV1 and MV2, viral_confidence High, length >= 5000, not concatemer, contamination <= 10",
}


def _source(stage: str) -> dict:
    return {
        "source_name": "metavr",
        "source_release": "MetaVR (2026-08-26 METAVR.fna.bgz / 2026-06-24 METAVR_main_table.parquet); "
                          f"{stage} selection by c79_metavr_mv2_eukv.py (DEC-65, DEC-79)",
        "license": c62.METAVR_LICENSE,
        "redistributable": False,
    }


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
        case = " ".join(f"WHEN list_contains(tok, '{t}') THEN '{n}'" for n, t in SLICE_TAXA)
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


# --- selection ----------------------------------------------------------------------------------------------------

REP_ORDER = "qrank ASC, comp_cap DESC, contam_rep ASC, length DESC, uvig ASC"  # C-62's representative order
SEL_COLS = ("uvig, votu, quality, viral_confidence, uvig_topology, completeness, checkv_contamination, length, "
            "genetic_code, genome_type, source, ecosystem, host_taxonomy, taxon")


def _connect(work: Path | None, threads: int, memory_limit: str, name: str) -> duckdb.DuckDBPyConnection:
    cfg = {"threads": int(threads), "memory_limit": memory_limit}
    if work is not None:
        d = Path(work) / name
        d.mkdir(parents=True, exist_ok=True)
        cfg["temp_directory"] = str(d)
    return duckdb.connect(config=cfg)


def _build_tables(con: duckdb.DuckDBPyConnection, parquet: Path, left_gz: Path) -> dict:
    """Tables `mv2`, `sel` (MVEUKV) and `rest` (its remainder); refuses a left list that does not fit the table."""
    case = " ".join(f"WHEN list_contains(tok, '{t}') THEN '{n}'" for n, t in EUK_TAXA)
    ordc = " ".join(f"WHEN '{n}' THEN {i}" for i, (n, _) in enumerate(EUK_TAXA))
    con.execute(f"""
        CREATE TEMP TABLE raw AS SELECT *, CASE taxon {ordc} ELSE NULL END AS taxon_ord FROM (
        SELECT uvig, votu, coalesce(quality, '') AS quality, coalesce(viral_confidence, '') AS viral_confidence,
               coalesce(uvig_topology, '') AS uvig_topology, coalesce(is_concatemer, false) AS is_concatemer,
               completeness, checkv_contamination, length::BIGINT AS length,
               coalesce(CAST(genetic_code AS VARCHAR), '') AS genetic_code, coalesce(genome_type, '') AS genome_type,
               coalesce(source, '') AS source, coalesce(Ecosystem, '') AS ecosystem,
               coalesce(host_taxonomy, '') AS host_taxonomy,
               CASE quality WHEN 'Complete' THEN 0 WHEN 'High-quality' THEN 1 WHEN 'Medium-quality' THEN 2
                    WHEN 'Low-quality' THEN 3 WHEN 'Not-determined' THEN 4 ELSE 5 END AS qrank,
               least(coalesce(completeness, -1.0), 100.0) AS comp_cap,
               coalesce(checkv_contamination, 0.0) AS contam_eff,
               coalesce(checkv_contamination, 999.0) AS contam_rep,
               CASE {case} ELSE NULL END AS taxon
        FROM (SELECT *, string_split(coalesce(ictv_taxonomy, ''), ';') AS tok FROM read_parquet({_sql_lit(parquet)})))
    """)
    con.execute(
        f"CREATE TEMP TABLE lft AS SELECT reason, mv2_rank::BIGINT AS mv2_rank, uvig FROM read_csv("
        f"{_sql_lit(left_gz)}, delim='\t', header=true, quote='', all_varchar=true)")
    n_raw, d_raw = con.execute("SELECT count(*), count(DISTINCT uvig) FROM raw").fetchone()
    if n_raw != d_raw:
        raise Refused(f"the table lists {n_raw - d_raw} uvig more than once")
    n_l, d_l = con.execute("SELECT count(*), count(DISTINCT uvig) FROM lft").fetchone()
    if n_l != d_l:
        raise Refused(f"the left list lists {n_l - d_l} uvig more than once")
    stray = con.execute("SELECT uvig FROM lft ANTI JOIN raw USING (uvig) ORDER BY uvig").fetchall()
    if stray:
        raise Refused(f"{len(stray)} UViGs of the left list are not in the table; first {stray[0][0]}")
    n, lo, hi, dn = con.execute(
        "SELECT count(*), min(mv2_rank), max(mv2_rank), count(DISTINCT mv2_rank) FROM lft "
        "WHERE reason = 'mv2_spillover'").fetchone()
    if n and (lo != 1 or hi != n or dn != n):
        raise Refused(f"mv2_rank of the {n} mv2_spillover rows is not 1..{n} (min {lo}, max {hi}, distinct {dn})")
    con.execute(f"""
        CREATE TEMP TABLE mv2 AS
        SELECT r.*, l.mv2_rank AS rank FROM lft l JOIN raw r USING (uvig) WHERE l.reason = 'mv2_spillover'
    """)
    con.execute("""
        CREATE TEMP TABLE touched AS
        SELECT DISTINCT r.votu FROM raw r LEFT JOIN lft l USING (uvig)
        WHERE l.uvig IS NULL OR l.reason = 'mv2_spillover'
    """)
    con.execute("CREATE TEMP TABLE eukrows AS SELECT r.* FROM raw r ANTI JOIN touched t USING (votu) "
                "WHERE r.taxon IS NOT NULL")
    con.execute(f"CREATE TEMP TABLE elig AS SELECT * FROM eukrows WHERE NOT is_concatemer "
                f"AND contam_eff <= {MAX_CONTAMINATION}")
    con.execute(f"""
        CREATE TEMP TABLE sel0 AS
        SELECT * EXCLUDE (rn) FROM (
            SELECT *, row_number() OVER (PARTITION BY votu ORDER BY {REP_ORDER}) AS rn FROM elig
            WHERE viral_confidence = '{CONFIDENCE}' AND length >= {MIN_LEN}
        ) WHERE rn = 1
    """)
    con.execute(f"CREATE TEMP TABLE sel AS SELECT *, row_number() OVER (ORDER BY taxon_ord, {REP_ORDER})::BIGINT "
                "AS rank FROM sel0")
    con.execute(f"""
        CREATE TEMP TABLE rest AS
        SELECT *, CASE WHEN eligible = 0 THEN 'no_eligible_member' WHEN length < {MIN_LEN} THEN 'lt5kb'
                       ELSE 'low_confidence_ge5kb' END AS reason
        FROM (
            SELECT e.*, row_number() OVER (PARTITION BY e.votu ORDER BY eligible DESC, {REP_ORDER}) AS rn
            FROM (SELECT *, CASE WHEN NOT is_concatemer AND contam_eff <= {MAX_CONTAMINATION} THEN 1 ELSE 0 END
                  AS eligible FROM eukrows) e
            ANTI JOIN sel s USING (votu)
        ) WHERE rn = 1
    """)
    q = con.execute
    stats = {
        "table_rows": int(n_raw), "left_rows": int(n_l), "mv1_rows": int(n_raw - n_l),
        "mv2_uvigs": int(q("SELECT count(*) FROM mv2").fetchone()[0]),
        "mv2_bases": int(q("SELECT coalesce(sum(length), 0) FROM mv2").fetchone()[0]),
        "mveukv_uvigs": int(q("SELECT count(*) FROM sel").fetchone()[0]),
        "mveukv_bases": int(q("SELECT coalesce(sum(length), 0) FROM sel").fetchone()[0]),
        "mveukv_by_taxon": {t: {"votus": int(c), "bp": int(b)} for t, c, b in q(
            "SELECT taxon, count(*), sum(length) FROM sel GROUP BY 1 ORDER BY 1").fetchall()},
        "mveukv_by_quality": {t: {"votus": int(c), "bp": int(b)} for t, c, b in q(
            "SELECT quality, count(*), sum(length) FROM sel GROUP BY 1 ORDER BY 1").fetchall()},
        "euk_votus_untouched": int(q("SELECT count(DISTINCT votu) FROM eukrows").fetchone()[0]),
        "left_by_reason": {r: {"votus": int(c), "bp": int(b)} for r, c, b in q(
            "SELECT reason, count(*), sum(length) FROM rest GROUP BY 1 ORDER BY 1").fetchall()},
        "left_by_reason_taxon": [(r, t, int(c), int(b)) for r, t, c, b in q(
            "SELECT reason, taxon, count(*), sum(length) FROM rest GROUP BY 1, 2 ORDER BY 1, 2").fetchall()],
    }
    if stats["mveukv_uvigs"] + int(q("SELECT count(*) FROM rest").fetchone()[0]) != stats["euk_votus_untouched"]:
        raise Refused("MVEUKV accounting: selected + left != untouched eukaryotic-virus vOTUs")
    return stats


def _rows(con, sql: str) -> list[dict]:
    cur = con.execute(sql)
    cols = [d[0] for d in cur.description]
    return [dict(zip(cols, r)) for r in cur.fetchall()]


LEFT_EUKV_HEAD = ["reason", "taxon", "votu", "uvig", "quality", "viral_confidence", "length", "completeness",
                  "checkv_contamination", "is_concatemer"]


def _write_left(con, path: Path) -> int:
    cur = con.execute("SELECT reason, taxon, votu, uvig, quality, viral_confidence, length, completeness, "
                      "checkv_contamination, is_concatemer FROM rest ORDER BY reason, taxon_ord, votu")
    n = 0
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        gz.write(("\t".join(LEFT_EUKV_HEAD) + "\n").encode())
        while batch := cur.fetchmany(100_000):
            gz.write("".join(
                f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]}\t{c62._fmt_num(r[7])}\t"
                f"{c62._fmt_num(r[8])}\t{str(bool(r[9])).lower()}\n" for r in batch).encode())
            n += len(batch)
    return n


def select_lists(parquet: Path, left_gz: Path, *, work: Path | None = None, threads: int = 2,
                 memory_limit: str = "4GB") -> tuple[list[dict], list[dict], list[dict], dict]:
    """(MV2 rows in mv2_rank order, MVEUKV rows in rank order, MVEUKV's remainder by vOTU, stats): for tests and
    small tables (the pod's run streams them)."""
    con = _connect(work, threads, memory_limit, "duck_sel")
    try:
        stats = _build_tables(con, Path(parquet), Path(left_gz))
        mv2 = _rows(con, f"SELECT {SEL_COLS}, rank AS mv2_rank FROM mv2 ORDER BY rank")
        sel = _rows(con, f"SELECT {SEL_COLS}, rank FROM sel ORDER BY rank")
        rest = _rows(con, "SELECT votu, reason, uvig, taxon, length FROM rest ORDER BY reason, taxon_ord, votu")
        return mv2, sel, rest, stats
    finally:
        con.close()


# --- staging ------------------------------------------------------------------------------------------------------

RECORDS_HEAD = [
    "ordinal", "source", "genome", "accession", "votu", "rank", "taxon", "contig_uid", "length", "quality",
    "viral_confidence", "completeness", "contamination", "topology", "source_genetic_code", "genetic_code",
    "circular", "genome_type", "origin_source", "ecosystem", "host_taxonomy", "dtr_class", "source_contig_uid",
    "source_length", "trimmed_bases",
]
META_FILES = ("genomes.json", "objects.json", "objects_staging.json", "records.tsv", "repeated_records.tsv",
              "shared_records.tsv", "later_records.tsv", "trims.tsv")


def stage(parquet: Path, left_gz: Path, fna: Path, meta_out: Path, fasta_out: Path, *, root: Path | None,
          root_stages, owner_inputs: dict[str, Path], later_inputs: dict[str, Path], inputs_home: Path,
          retrieved_at: str = "2026-10-23T15:22:00Z", work: Path | None = None, threads: int = 2,
          memory_limit: str = "4GB", use_bgzip: bool = False) -> dict:
    """Select MV2 and MVEUKV, stream FNA once and write both stages; returns the report (module docstring)."""
    meta_out, fasta_out = Path(meta_out), Path(fasta_out)
    stages = tuple(FASTA_NAME)
    for st in stages:
        for p in (meta_out / st, fasta_out / st / FASTA_NAME[st]):
            if p.exists():
                raise Refused(f"{p} exists")
    tmp = {st: meta_out / f"{st}.tmp" for st in stages}
    partial = {st: fasta_out / st / f".{FASTA_NAME[st]}.partial" for st in stages}
    work_dir = Path(work) if work is not None else meta_out / "_work"
    for st in stages:
        shutil.rmtree(tmp[st], ignore_errors=True)
        tmp[st].mkdir(parents=True)
        partial[st].parent.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    try:
        rep = _stage(Path(parquet), Path(left_gz), Path(fna), tmp, partial, work_dir, root, list(root_stages),
                     owner_inputs, later_inputs, Path(inputs_home), fasta_out, retrieved_at, threads, memory_limit,
                     use_bgzip)
    except BaseException:
        for st in stages:
            shutil.rmtree(tmp[st], ignore_errors=True)
            partial[st].unlink(missing_ok=True)
        raise
    if work is None:
        shutil.rmtree(work_dir, ignore_errors=True)
    for st in stages:
        partial[st].rename(fasta_out / st / FASTA_NAME[st])
        tmp[st].rename(meta_out / st)
    for st in stages:
        files = {n: meta_out / st / n for n in META_FILES}
        if st == "MVEUKV":
            files[LEFT_EUKV] = meta_out / st / LEFT_EUKV
        files[FASTA_NAME[st]] = fasta_out / st / FASTA_NAME[st]
        rep["stages"][st]["sha256"] = {n: ids.sha256_file(p).hex() for n, p in files.items()}
        rep["stages"][st]["bytes"] = {n: p.stat().st_size for n, p in files.items()}
    return rep


def _stage(parquet, left_gz, fna, tmp, partial, work_dir, root, root_stages, owner_inputs, later_inputs,
           inputs_home, fasta_out, retrieved_at, threads, memory_limit, use_bgzip) -> dict:
    con = _connect(work_dir, threads, memory_limit, "duck_sel")
    chosen: dict[str, tuple] = {}  # uvig -> (stage, rank, votu, quality, conf, topology, comp, contam, length,
    #                                         genetic_code, genome_type, source, ecosystem, host, taxon)
    try:
        stats = _build_tables(con, parquet, left_gz)
        left_rows = _write_left(con, tmp["MVEUKV"] / LEFT_EUKV)
        for st, table in (("MV2", "mv2"), ("MVEUKV", "sel")):
            cur = con.execute(f"SELECT rank, votu, quality, viral_confidence, uvig_topology, completeness, "
                              f"checkv_contamination, length, genetic_code, genome_type, source, ecosystem, "
                              f"host_taxonomy, coalesce(taxon, ''), uvig FROM {table}")
            while batch := cur.fetchmany(100_000):
                for r in batch:
                    if r[14] in chosen:
                        raise Refused(f"{r[14]} selected twice")
                    chosen[r[14]] = (st, *r[:14])
    finally:
        con.close()

    parts: dict[str, list[tuple]] = {}  # uvig -> [(cu, n, src_n, rec_ord, dtr_cls, src_cu, k, accession)]
    circ: dict[str, bool] = {}
    trims: dict[str, list[tuple]] = {st: [] for st in FASTA_NAME}
    dtr_classes: dict[str, Counter] = {st: Counter() for st in FASTA_NAME}
    n_scanned = n_written = 0
    stream, close_stream = c62._open_fasta_stream(fna, use_bgzip=use_bgzip)
    writers = {st: open(partial[st], "wb", buffering=1 << 22) for st in FASTA_NAME}
    try:
        for header, r in c62.v2mod.iter_fasta_taken(stream, lambda hd: c62._uvig_from_header(hd) in chosen):
            n_scanned += 1
            if r is None:
                continue
            uvig = c62._uvig_from_header(header)
            meta = chosen[uvig]
            st, topo, expected_len = meta[0], meta[5], int(meta[8])
            src_len = len(r.seq)
            src_cu = ids.contig_uid(r.seq)
            if uvig not in parts and src_len == expected_len and topo == "Direct terminal repeat":
                seq, is_c, dtr_cls, k = c62.apply_dec65(topo, r.seq)
                dtr_classes[st][dtr_cls] += 1
                circ[uvig] = is_c
            else:
                seq, dtr_cls, k = r.seq, "", 0
                circ[uvig] = False
            writers[st].write(b"\n".join(fasta.format_record(r.header, seq, r.line_width)) + b"\n")
            cu = ids.contig_uid(seq) if k else src_cu
            if k:
                trims[st].append(("metavr", uvig, r.accession, src_cu.hex(), cu.hex(), src_len, k, k, len(seq)))
            parts.setdefault(uvig, []).append((cu, len(seq), src_len, n_written, dtr_cls, src_cu, k, r.accession))
            n_written += 1
    finally:
        for w in writers.values():
            w.close()
        close_stream()

    missing = sorted(set(chosen) - set(parts))
    if missing:
        raise Refused(f"{len(missing)} selected UViGs have no FASTA record; first {missing[0]}")
    multi_dtr = sorted(u for u, p in parts.items() if chosen[u][5] == "Direct terminal repeat" and len(p) > 1)
    if multi_dtr:
        raise Refused(f"{len(multi_dtr)} Direct terminal repeat UViGs have >1 FASTA record; first {multi_dtr[0]}")

    order = {st: sorted((u for u in parts if chosen[u][0] == st), key=lambda u: chosen[u][1]) for st in FASTA_NAME}
    all_uids = work_dir / "all_uids.tsv"
    with open(all_uids, "w") as fh:
        fh.write("contig_uid\n")
        for cu in sorted({p[0] for ps in parts.values() for p in ps}):
            fh.write(cu.hex() + "\n")
    owners = c62._owners_from_records(all_uids, root, root_stages, owner_inputs, work=work_dir, threads=threads,
                                      memory_limit=memory_limit)
    later = c62._owners_from_records(all_uids, None, [], later_inputs, work=work_dir, threads=threads,
                                     memory_limit=memory_limit) if later_inputs else {}
    mv2_uids = {p[0] for u in order["MV2"] for p in parts[u]}

    report: dict = {"selection": stats, "streams": 1, "fasta_records_scanned": n_scanned,
                    "fasta_records_written": n_written, "mveukv_left_rows": left_rows, "stages": {}}
    for st in FASTA_NAME:
        report["stages"][st] = _write_stage(
            st, order[st], parts, circ, chosen, trims[st], dtr_classes[st], owners, later,
            mv2_uids if st == "MVEUKV" else set(), tmp[st], inputs_home, fasta_out, retrieved_at)
    return report


def _write_stage(st, order, parts, circ, chosen, trims, dtr_classes, owners, later, earlier_uids, out,
                 inputs_home, fasta_out, retrieved_at) -> dict:
    with open(out / "records.tsv", "w", encoding="utf-8", buffering=1 << 20) as rf:
        rf.write("\t".join(RECORDS_HEAD) + "\n")
        ordinal = 0
        for uvig in order:
            (_, rank, votu, qual, conf, topo, comp, contam, _len, src_gc, gtype, orig, eco, host, taxon) = chosen[uvig]
            for cu, n, src_n, _, dtr_cls, src_cu, k, acc in parts[uvig]:
                rf.write(f"{ordinal}\tmetavr\t{uvig}\t{acc}\t{votu}\t{rank}\t{taxon}\t{cu.hex()}\t{n}\t{qual}\t{conf}\t"
                         f"{c62._fmt_num(comp)}\t{c62._fmt_num(contam)}\t{topo}\t{src_gc}\t11\t"
                         f"{str(circ[uvig]).lower()}\t{gtype}\t{orig}\t{eco}\t{host}\t{dtr_cls}\t{src_cu.hex()}\t"
                         f"{src_n}\t{k}\n")
                ordinal += 1
    listing = make_inputs.Listing()
    genomes, repeated, shared, held_later = [], [], [], []
    first_of: dict[bytes, str] = {}
    contig_len: dict[bytes, int] = {}
    by_taxon: dict[str, list[int]] = {}
    ordinal = 0
    for uvig in order:
        uids = [p[0] for p in parts[uvig]]
        t = by_taxon.setdefault(chosen[uvig][14] or "-", [0, 0])
        t[0] += 1
        t[1] += sum(p[1] for p in parts[uvig])
        for cu, n, *_ in parts[uvig]:
            contig_len.setdefault(cu, n)
            if cu in owners:
                shared.append((ordinal, uvig, cu.hex(), owners[cu]))
            elif cu in later:
                held_later.append((ordinal, uvig, cu.hex(), later[cu]))
            elif cu in earlier_uids:
                held_later.append((ordinal, uvig, cu.hex(), "MV2"))
            ordinal += 1
        g = listing.add(uids, 11, "virus", "vOTU")
        guid = ids.genome_uid(list(dict.fromkeys(uids)))
        if g is None:
            repeated.append((ordinal - len(uids), uvig, guid.hex(), first_of[guid]))
            continue
        g["circular"] = [uids[0].hex()] if circ[uvig] else []
        genomes.append(g)
        first_of[guid] = uvig
    make_inputs._write_tsv(out / "repeated_records.tsv", ["ordinal", "genome", "genome_uid", "first_genome"], repeated)
    make_inputs._write_tsv(out / "shared_records.tsv", ["ordinal", "genome", "contig_uid", "owner_stage"], shared)
    make_inputs._write_tsv(out / "later_records.tsv", ["ordinal", "genome", "contig_uid", "holder"], held_later)
    make_inputs._write_tsv(out / "trims.tsv", c62.TRIM_HEAD, sorted(trims, key=lambda t: (t[0], t[1], t[2])))
    obj = {**_source(st), "retrieved_at": retrieved_at}
    (out / "objects.json").write_text(json.dumps([{"path": str(inputs_home / st / FASTA_NAME[st]), **obj}], indent=1))
    (out / "objects_staging.json").write_text(
        json.dumps([{"path": str(fasta_out / st / FASTA_NAME[st]), **obj}], indent=1))
    (out / "genomes.json").write_text(json.dumps(genomes))

    def tally(holders: dict[bytes, str]):
        c, b = Counter(), Counter()
        for cu, n in contig_len.items():
            if cu in holders:
                c[holders[cu]] += 1
                b[holders[cu]] += n
        return dict(sorted(c.items())), dict(sorted(b.items()))

    own_c, own_b = tally(owners)
    lat = {cu: later.get(cu) or "MV2" for cu in contig_len
           if cu not in owners and (cu in later or cu in earlier_uids)}
    lat_c, lat_b = tally(lat)
    mism = [u for u in order if sum(p[2] for p in parts[u]) != int(chosen[u][8])]
    return {
        "rule": RULE[st],
        "uvigs": len(order),
        "records": sum(len(parts[u]) for u in order),
        "by_taxon": {k: {"uvigs": v[0], "bp": v[1]} for k, v in sorted(by_taxon.items())},
        "source_bases": sum(p[2] for u in order for p in parts[u]),
        "staged_bases": sum(p[1] for u in order for p in parts[u]),
        "trims": len(trims), "trimmed_bases": sum(t[7] for t in trims),
        "dtr_classes": dict(sorted(dtr_classes.items())),
        "length_mismatches_count": len(mism), "length_mismatches": mism[:50],
        "genomes": len(genomes), "contigs": len(contig_len), "contig_bases": sum(contig_len.values()),
        "circular_genomes": sum(1 for g in genomes if g["circular"]),
        "multi_record_genomes": sum(1 for u in order if len(parts[u]) > 1),
        "repeated_genomes": len(repeated), "shared_records": len(shared), "later_records": len(held_later),
        "owner_contigs": own_c, "owner_bases": own_b, "later_contigs": lat_c, "later_bases": lat_b,
        "new_contigs": sum(1 for cu in contig_len if cu not in owners),
        "new_bases": sum(n for cu, n in contig_len.items() if cu not in owners),
        "new_after_later_contigs": sum(1 for cu in contig_len if cu not in owners and cu not in lat),
        "new_after_later_bases": sum(n for cu, n in contig_len.items() if cu not in owners and cu not in lat),
        **listing.report(),
    }


# --- the ANI sample against GV1 (card C-48c's method: c39_ani's skani settings and rules) --------------------------

def ani_summary(table: Path, sample: dict[str, int], taxon: dict[str, str], refs: dict[str, int]) -> dict:
    """Per rule of `c39_ani.hits`, how many sampled queries meet the vOTU threshold with a reference, by taxon and
    length class."""
    hh = c39_ani.hits(table, {**refs, **sample})
    ref_names = set(refs)
    out: dict = {"queries": len(sample)}
    for rule, h in hh.items():
        met = {q for q in sample if (h.get(q, set()) - {q}) & ref_names}
        by_t, by_l = {}, {}
        for q, n in sample.items():
            for d, k in ((by_t, taxon.get(q, "-")), (by_l, c39_ani.length_class(n))):
                e = d.setdefault(k, {"queries": 0, "meet": 0})
                e["queries"] += 1
                e["meet"] += q in met
        out[rule] = {"meet": len(met), "by_taxon": dict(sorted(by_t.items())), "by_length": dict(sorted(by_l.items()))}
    return out


def ani_sample(meta: Path, fasta_path: Path, gv1_fasta: Path, out: Path, *, threads: int = 3,
               every: int = 25) -> dict:
    """Every EVERY-th single-record, non-repeated MVEUKV genome (FASTA order) against all of GV1's records."""
    meta, out = Path(meta), Path(out)
    work = out / "work"
    work.mkdir(parents=True, exist_ok=True)
    with open(meta / "records.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    per_genome = Counter(r["genome"] for r in rows)
    tax_of = {r["genome"]: r.get("taxon", "") for r in rows}
    with open(meta / "repeated_records.tsv", newline="") as fh:
        repeated = {r["genome"] for r in csv.DictReader(fh, delimiter="\t")}
    sample, taxon, k = {}, {}, 0
    sample_fa = work / "mveukv_sample.fna"
    with open(sample_fa, "w") as w:
        for r in fasta.iter_records(Path(fasta_path)):
            g = r.accession.split("|", 1)[0]
            if per_genome.get(g) != 1 or g in repeated:
                continue
            k += 1
            if k % every:
                continue
            w.write(f">{r.accession}\n{r.seq.decode()}\n")
            sample[r.accession] = len(r.seq)
            taxon[r.accession] = tax_of[g]
    refs = c39_ani.all_lengths(Path(gv1_fasta))
    if set(refs) & set(sample):
        raise Refused("record names are shared between the sample and GV1; names would be ambiguous")
    cpu = c39_ani.skani(threads, [sample_fa], [Path(gv1_fasta)], out / "a.tsv")
    res = {"every": every, "threads": threads, "samples": len(sample), "sample_bases": sum(sample.values()),
           "eligible_genomes": k, "reference_records": len(refs), "reference_bases": sum(refs.values()),
           "threshold": {"ani_min": c39_ani.ANI_MIN, "af_shorter_min": c39_ani.AF_MIN}, "cpu_s": round(cpu, 1),
           "vs_gv1": ani_summary(out / "a.tsv", sample, taxon, refs)}
    (out / "summary.json").write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    return res


# --- command line -------------------------------------------------------------------------------------------------

def _check(stage_name: str, meta: Path) -> dict:
    p = subprocess.run([sys.executable, str(Path(__file__).resolve().parent / "s2_check_inputs.py"), stage_name,
                        str(meta / "genomes.json"), str(meta / "objects_staging.json")],
                       capture_output=True, text=True, check=False)
    doc = json.loads(p.stdout) if p.stdout.strip().startswith("{") else {"stdout": p.stdout}
    return {"exit": p.returncode, "problems": doc.get("problems"), "gene_call_mode": doc.get("gene_call_mode"),
            "genomes": doc.get("genomes"), "genome_input_rows": doc.get("genome_input_rows"),
            "stderr": p.stderr[-2000:]}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Stage MetaVR MV2 and MVEUKV inputs (card C-79, DEC-79)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    st = sub.add_parser("stage")
    for n in ("parquet", "left_gz", "fna", "meta_out", "fasta_out"):
        st.add_argument(n, type=Path)
    st.add_argument("--root", type=Path, default=None)
    st.add_argument("--root-stages", default=",".join(ROOT_STAGES))
    st.add_argument("--chain-home", type=Path, default=None,
                    help="the r26.11.4 stages' inputs, read in a26's chain order as owners")
    st.add_argument("--later", action="append", default=[], help="an r26.11.5 stage before MV2 (under chain home)")
    st.add_argument("--inputs-home", type=Path, default=Path("/tmp/if2453/s2/inputs"))
    st.add_argument("--retrieved-at", default="2026-10-23T15:22:00Z")
    st.add_argument("--work", type=Path, default=None)
    st.add_argument("--threads", type=int, default=2)
    st.add_argument("--memory-limit", default="4GB")
    st.add_argument("--bgzip", action="store_true")
    st.add_argument("--report-out", type=Path, default=None)
    an = sub.add_parser("ani")
    for n in ("meta", "fasta", "gv1_fasta", "out"):
        an.add_argument(n, type=Path)
    an.add_argument("--threads", type=int, default=3)
    an.add_argument("--every", type=int, default=25)
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
    try:
        if a.cmd == "stage":
            owners = {s: a.chain_home / s for s, _ in a26.R26_11_4_CHAIN} if a.chain_home is not None else {}
            later = {s: a.chain_home / s for s in a.later} if a.chain_home is not None else {}
            rep = stage(a.parquet, a.left_gz, a.fna, a.meta_out, a.fasta_out, root=a.root,
                        root_stages=[s for s in a.root_stages.split(",") if s], owner_inputs=owners,
                        later_inputs=later, inputs_home=a.inputs_home, retrieved_at=a.retrieved_at, work=a.work,
                        threads=a.threads, memory_limit=a.memory_limit, use_bgzip=a.bgzip)
            rep["owner_stages"] = [*([s for s in a.root_stages.split(",") if s] if a.root else []), *owners]
            rep["later_stages"] = [*later]
            rep["s2_check"] = {s: _check(s, a.meta_out / s) for s in FASTA_NAME}
            txt = json.dumps(rep, indent=1, sort_keys=True)
            if a.report_out is not None:
                a.report_out.parent.mkdir(parents=True, exist_ok=True)
                a.report_out.write_text(txt + "\n")
            print(txt)
            return 0 if all(c["exit"] == 0 and c["problems"] == [] for c in rep["s2_check"].values()) else 1
        if a.cmd == "ani":
            print(json.dumps(ani_sample(a.meta, a.fasta, a.gv1_fasta, a.out, threads=a.threads, every=a.every),
                             sort_keys=True))
            return 0
    except Refused as exc:
        print(f"refused: {exc}", file=sys.stderr)
        return 2
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
