"""What the genes that p3 changes look like, against the stored calls (card C-21f; after compare.py).

    python detail.py --feature-protein FP_DIR --domains DOM_DIR --p2 RUN_DIR --p3 RUN_DIR --out OUT.json

FP_DIR is the stage's stored `l1/feature_protein/gene_call_version=<v>/stage=<S>` directory, DOM_DIR its
`l1/domains/library=<lib>/stage=<S>` directory (both read-only). RUN_DIRs are recall.py outputs (one shard each).
Classes and the stop rule are compare.py's.

No truth set is available (card C-21f), so a stored Pfam domain is used only as a hint that a gene codes for a
real protein: a lost gene with a domain is more likely a real gene lost; a start change that moves the start past
the first residue of a stored domain cuts that domain. Proteins shared with an earlier stage have no domain rows
in this stage's partition, so the share is a lower bound.

Reports, per class:
- lost, gained and all stored genes: count, median length in nucleotides, share with a stored domain (lost, all);
- start-only changes: how many make the gene shorter or longer under p3, and how many shorter ones cut a domain;
- per contig, the share of stored genes identical under p3, by contig length bin: minimum, median, and contigs under
  99%.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import duckdb

from compare import CLASS_SQL, _lits  # same directory


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--feature-protein", type=Path, required=True)
    ap.add_argument("--domains", type=Path, required=True)
    ap.add_argument("--p2", type=Path, required=True)
    ap.add_argument("--p3", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    fp = sorted(a.feature_protein.glob("part-*.parquet"))
    dom = sorted(a.domains.glob("part-*.parquet"))
    con = duckdb.connect(config={"threads": 1})
    con.execute(f"CREATE VIEW contigs AS SELECT contig_uid, contig_len, {CLASS_SQL} AS cls "
                f"FROM read_parquet({_lits([a.p2 / 'contigs.parquet'])})")
    stop = ('CASE WHEN g.strand = 1 THEN (CASE WHEN g."end" > c.contig_len THEN g."end" - c.contig_len '
            'ELSE g."end" END) ELSE g.start END')
    con.execute(f"CREATE TABLE dp AS SELECT protein_uid, min(env_start) AS dom_first FROM read_parquet({_lits(dom)}) "
                "GROUP BY protein_uid")
    con.execute(f"""CREATE TABLE s AS SELECT g.contig_uid, g.start, g."end", g.strand, g.protein_uid, c.cls,
        c.contig_len, {stop} AS stop, dp.dom_first IS NOT NULL AS has_dom, dp.dom_first
        FROM read_parquet({_lits(fp)}) g JOIN contigs c USING (contig_uid) LEFT JOIN dp USING (protein_uid)""")
    con.execute(f"""CREATE TABLE n AS SELECT g.contig_uid, g.start, g."end", g.strand, c.cls, {stop} AS stop
        FROM read_parquet({_lits([a.p3 / 'genes.parquet'])}) g JOIN contigs c USING (contig_uid)""")
    res: dict = {"feature_protein_files": [str(p) for p in fp], "domain_files": [str(p) for p in dom]}

    def q(sql: str) -> list[tuple]:
        return con.execute(sql).fetchall()

    res["stored"] = {cls: {"n": n, "median_nt": med, "with_domain": d} for cls, n, med, d in q(
        'SELECT coalesce(cls, \'all\'), count(*), median("end" - start), avg(has_dom::INT) FROM s '
        "GROUP BY ROLLUP (cls) ORDER BY 1")}
    res["lost"] = {cls: {"n": n, "median_nt": med, "with_domain": d, "with_domain_n": dn} for cls, n, med, d, dn in q(
        'SELECT coalesce(cls, \'all\'), count(*), median("end" - start), avg(has_dom::INT), sum(has_dom::INT) '
        "FROM s ANTI JOIN n USING (contig_uid, strand, stop) GROUP BY ROLLUP (cls) ORDER BY 1")}
    res["gained"] = {cls: {"n": n, "median_nt": med} for cls, n, med in q(
        'SELECT coalesce(cls, \'all\'), count(*), median("end" - start) '
        "FROM n ANTI JOIN s USING (contig_uid, strand, stop) GROUP BY ROLLUP (cls) ORDER BY 1")}
    # start-only: same stop, different start; shift > 0 when p3's gene is shorter (its start lies downstream)
    res["start_only"] = {cls: {"n": k, "shorter": sh, "longer": lo, "median_abs_shift_nt": med,
                               "shorter_cutting_a_domain": cut}
                         for cls, k, sh, lo, med, cut in q("""
        WITH j AS (SELECT s.cls, s.dom_first,
                     CASE WHEN s.strand = 1 THEN n.start - s.start ELSE s."end" - n."end" END AS shift
                   FROM s JOIN n USING (contig_uid, strand, stop) WHERE s.start <> n.start OR s."end" <> n."end")
        SELECT coalesce(cls, 'all'), count(*), count(*) FILTER (WHERE shift > 0), count(*) FILTER (WHERE shift < 0),
               median(abs(shift)), count(*) FILTER (WHERE shift > 0 AND dom_first IS NOT NULL
                                                    AND dom_first < shift / 3)
        FROM j GROUP BY ROLLUP (cls) ORDER BY 1""")}
    res["per_contig"] = [dict(zip(("cls", "bin", "contigs", "min_identical", "median_identical", "under_99"), r))
                         for r in q("""
        WITH per AS (
          SELECT c.cls, c.contig_len, count(*) AS k, count(n.stop) AS ident
          FROM contigs c JOIN s USING (contig_uid) LEFT JOIN n USING (contig_uid, start, "end", strand)
          GROUP BY c.cls, c.contig_uid, c.contig_len)
        SELECT cls, CASE WHEN contig_len < 10000 THEN 'a <10k' WHEN contig_len < 100000 THEN 'b 10k-100k'
                         WHEN contig_len < 200000 THEN 'c 100k-200k' WHEN contig_len < 500000 THEN 'd 200k-500k'
                         WHEN contig_len < 1000000 THEN 'e 500k-1M' WHEN contig_len < 3000000 THEN 'f 1M-3M'
                         ELSE 'g >=3M' END AS bin,
               count(*), min(ident / k), median(ident / k), count(*) FILTER (WHERE ident / k < 0.99)
        FROM per GROUP BY ALL ORDER BY cls, bin""")]
    a.out.write_text(json.dumps(res, indent=1, sort_keys=True, default=str) + "\n")
    print(json.dumps({k: v for k, v in res.items() if not k.endswith("_files")}, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
