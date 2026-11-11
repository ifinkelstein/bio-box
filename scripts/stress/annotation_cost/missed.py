"""Which truth pairs an MMseqs2 run misses (card 0-09a diagnostics; reads only).

    python missed.py <export_dir> <models.tsv> <hits.tsv> [--margin-min 50]

Prints, for truth pairs with margin (seq_bits - GA_seq) >= margin-min: recall by model length bin, and
the models with the most missed pairs.
"""

from __future__ import annotations

import argparse

import duckdb


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("export")
    ap.add_argument("models")
    ap.add_argument("hits")
    ap.add_argument("--margin-min", type=float, default=50)
    a = ap.parse_args()
    con = duckdb.connect()
    con.execute("SET threads=4")
    con.execute(f"""CREATE TABLE models AS SELECT acc, name, leng::INT leng, ga_seq::DOUBLE ga_seq
        FROM read_csv('{a.models}', delim='\t', header=true, all_varchar=true)""")
    con.execute(f"""CREATE TABLE tp AS SELECT uid, acc, max(seq_bits::DOUBLE) seq_bits FROM
        read_csv('{a.export}/truth.tsv', delim='\t', header=true, all_varchar=true) GROUP BY uid, acc""")
    con.execute(f"""CREATE TABLE mp AS SELECT DISTINCT column0 uid, column1 acc
        FROM read_csv('{a.hits}', delim='\t', header=false, all_varchar=true)""")
    con.execute(f"""CREATE TABLE t AS SELECT tp.*, m.name, m.leng, tp.seq_bits - m.ga_seq margin,
        (mp.uid IS NOT NULL) found FROM tp JOIN models m USING (acc) LEFT JOIN mp USING (uid, acc)
        WHERE tp.seq_bits - m.ga_seq >= {a.margin_min}""")
    print("model length bin, truth pairs, recall")
    for row in con.execute("""SELECT CASE WHEN leng < 50 THEN '<50' WHEN leng < 100 THEN '50-99'
            WHEN leng < 200 THEN '100-199' WHEN leng < 400 THEN '200-399' ELSE '>=400' END b, count(*),
            avg(found::INT) FROM t GROUP BY b ORDER BY min(leng)""").fetchall():
        print(row)
    print("models with most missed pairs: acc, name, leng, truth pairs, missed, median margin of missed")
    for row in con.execute("""SELECT acc, name, any_value(leng), count(*), sum((NOT found)::INT) miss,
            median(margin) FILTER (WHERE NOT found) FROM t GROUP BY acc, name ORDER BY miss DESC LIMIT 20""").fetchall():
        print(row)
    print("missed pairs:", con.execute("SELECT sum((NOT found)::INT), count(*) FROM t").fetchone(),
          "models with any miss:", con.execute("SELECT count(DISTINCT acc) FROM t WHERE NOT found").fetchone()[0],
          "models with every pair missed:", con.execute(
              "SELECT count(*) FROM (SELECT acc FROM t GROUP BY acc HAVING bool_and(NOT found))").fetchone()[0],
          "models:", con.execute("SELECT count(DISTINCT acc) FROM t").fetchone()[0])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
