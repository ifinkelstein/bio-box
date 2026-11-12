"""Compare re-called genes with a stage's stored gene calls, per contig class (card C-21f).

    python compare.py --stored FEATURE_DIR --p2 RUN_DIR [RUN_DIR ...] --p3 RUN_DIR [RUN_DIR ...] --out OUT.json

FEATURE_DIR is the stage's `l1/feature/producer=<v>/stage=<S>` directory (read-only); each RUN_DIR is one shard
written by `recall.py`. Prints a Markdown summary and writes every number to OUT.json.

Genes are compared on the same contig. A stored gene and a re-called gene are:

- identical: same (start, end, strand);
- start-only different: same strand and same stop, different start. The stop is `end` on the + strand and `start`
  on the - strand; a gene that wraps a circular contig's origin (`end > contig_len`) has its + strand stop taken
  modulo the contig length, so a stop is one position on the circle;
- lost: a stored gene whose stop no re-called gene has (same contig and strand);
- gained: a re-called gene whose stop no stored gene has.

Prodigal reports at most one gene per stop, so on each side (contig, strand, stop) is unique; the script checks this
and stops if it is not. Then identical + start-only + lost = stored and identical + start-only + gained = re-called.

Contig classes, from the genome's kind, its total bases and the contig's length:

- `virus`: a contig of a viral genome (gv-meta under both policies);
- `small_genome`: a contig of a non-viral genome under 100 kb in all (meta under both);
- `plasmid_genome`: a contig of a genome of kind plasmid (meta under p2; single, self-trained, under p3 if 100 kb or
  more);
- `chromosome`: the longest contig of a prokaryote genome of 100 kb or more, itself 100 kb or more (single, trained
  on the genome, under p2; single, trained on itself, under p3);
- `other_ge100k`: any other contig of 100 kb or more of such a genome (as `chromosome`);
- `lt100k_circular`, `lt100k_linear`: a contig under 100 kb of such a genome (genome-trained single under p2; meta
  under p3), split by topology.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import duckdb

CLASS_SQL = """
CASE
  WHEN kind = 'virus' THEN 'virus'
  WHEN genome_bases < 100000 THEN 'small_genome'
  WHEN kind = 'plasmid' THEN 'plasmid_genome'
  WHEN contig_len >= 100000 AND is_largest THEN 'chromosome'
  WHEN contig_len >= 100000 THEN 'other_ge100k'
  WHEN is_circular THEN 'lt100k_circular'
  ELSE 'lt100k_linear'
END"""
ORDER = ["chromosome", "other_ge100k", "lt100k_circular", "lt100k_linear", "plasmid_genome", "small_genome",
         "virus", "all"]


def _lits(paths: list[Path]) -> str:
    return "[" + ", ".join("'" + str(p).replace("'", "''") + "'" for p in paths) + "]"


def _files(dirs: list[Path], name: str) -> list[Path]:
    out = []
    for d in dirs:
        if not (d / "done.json").is_file():
            raise SystemExit(f"{d}: no done.json; the run did not finish")
        out.append(d / name)
    return out


def _stops(con, rel: str) -> None:
    """A view `<rel>_s` of `rel` with the stop of each gene, and a check that (contig, strand, stop) is unique."""
    con.execute(f"""CREATE VIEW {rel}_s AS SELECT g.contig_uid, g.start, g."end", g.strand,
        CASE WHEN g.strand = 1 THEN (CASE WHEN g."end" > c.contig_len THEN g."end" - c.contig_len ELSE g."end" END)
             ELSE g.start END AS stop, c.cls
        FROM {rel} g JOIN contigs c USING (contig_uid)""")
    n, k = con.execute(f"SELECT count(*), count(DISTINCT (contig_uid, strand, stop)) FROM {rel}_s").fetchone()
    if n != k:
        raise SystemExit(f"{rel}: {n - k} genes share a (contig, strand, stop) with another; the comparison assumes "
                         "one gene per stop")


def compare(con, new: str) -> dict:
    """Per class: counts of stored and `new` genes, identical, start-only, lost, gained; contigs with any change."""
    rows = con.execute(f"""
        WITH s AS (SELECT contig_uid, start, "end", strand, stop FROM stored_s),
        n AS (SELECT contig_uid, start, "end", strand, stop FROM {new}_s),
        ns AS (SELECT contig_uid, count(*) AS k FROM s GROUP BY contig_uid),
        nn AS (SELECT contig_uid, count(*) AS k FROM n GROUP BY contig_uid),
        ni AS (SELECT contig_uid, count(*) AS k FROM s JOIN n USING (contig_uid, start, "end", strand)
               GROUP BY contig_uid),
        np AS (SELECT contig_uid, count(*) AS k FROM s JOIN n USING (contig_uid, strand, stop) GROUP BY contig_uid),
        per AS (
          SELECT c.cls, c.contig_uid, coalesce(ns.k, 0) AS n_stored, coalesce(nn.k, 0) AS n_new,
                 coalesce(ni.k, 0) AS n_ident, coalesce(np.k, 0) AS n_stop
          FROM contigs c LEFT JOIN ns USING (contig_uid) LEFT JOIN nn USING (contig_uid)
          LEFT JOIN ni USING (contig_uid) LEFT JOIN np USING (contig_uid))
        SELECT cls, count(*) AS contigs, sum(n_stored), sum(n_new), sum(n_ident), sum(n_stop),
               count(*) FILTER (WHERE n_ident <> n_stored OR n_ident <> n_new) AS contigs_changed
        FROM per GROUP BY ROLLUP (cls) ORDER BY cls""").fetchall()
    out = {}
    for cls, contigs, n_s, n_n, ident, stop, changed in rows:
        n_s, n_n, ident, stop = int(n_s or 0), int(n_n or 0), int(ident or 0), int(stop or 0)
        d = {"contigs": contigs, "contigs_changed": changed, "stored": n_s, "recalled": n_n, "identical": ident,
             "start_only": stop - ident, "lost": n_s - stop, "gained": n_n - stop}
        for key in ("identical", "start_only", "lost", "gained"):
            d[f"{key}_frac"] = d[key] / n_s if n_s else None
        out["all" if cls is None else cls] = d
    return out


def costs(con) -> dict:
    """Per class and policy: bases, CPU seconds of training and of finding genes (thread time, recall.py)."""
    out: dict = {}
    for pol in ("p2", "p3"):
        rows = con.execute(f"""SELECT cls, sum(contig_len), sum(train_cpu_s), sum(find_cpu_s),
                               string_agg(DISTINCT mode, ',' ORDER BY mode)
                               FROM {pol}_contigs GROUP BY ROLLUP (cls) ORDER BY cls""").fetchall()
        for cls, bases, train, find, modes in rows:
            out.setdefault("all" if cls is None else cls, {})[pol] = {
                "bases": int(bases), "train_cpu_s": train, "find_cpu_s": find, "cpu_s": train + find, "modes": modes}
    return out


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--stored", type=Path, required=True)
    ap.add_argument("--p2", type=Path, nargs="+", required=True)
    ap.add_argument("--p3", type=Path, nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    stored = sorted(a.stored.glob("part-*.parquet"))  # one partition directory, listed explicitly in the output
    con = duckdb.connect(config={"threads": 1})
    con.execute(f"CREATE VIEW stored AS SELECT contig_uid, start, \"end\", strand FROM read_parquet({_lits(stored)})")
    for pol, dirs in (("p2", a.p2), ("p3", a.p3)):
        con.execute(f"CREATE VIEW {pol} AS SELECT * FROM read_parquet({_lits(_files(dirs, 'genes.parquet'))})")
        con.execute(f"CREATE VIEW {pol}_contigs AS SELECT *, {CLASS_SQL} AS cls "
                    f"FROM read_parquet({_lits(_files(dirs, 'contigs.parquet'))})")
    # contig attributes are the same in both runs; the comparison needs every contig of the stage in each
    for q in ("SELECT count(*) FROM p2_contigs", "SELECT count(*) FROM p3_contigs",
              "SELECT count(DISTINCT contig_uid) FROM p2_contigs",
              "SELECT count(*) FROM p2_contigs a JOIN p3_contigs b USING (contig_uid) WHERE a.contig_len = b.contig_len "
              "AND a.cls = b.cls"):
        print(f"check: {q}: {con.execute(q).fetchone()[0]}")
    con.execute("CREATE VIEW contigs AS SELECT contig_uid, contig_len, cls FROM p2_contigs")
    n_unknown = con.execute("SELECT count(*) FROM stored WHERE contig_uid NOT IN (SELECT contig_uid FROM contigs)"
                            ).fetchone()[0]
    if n_unknown:
        raise SystemExit(f"{n_unknown} stored genes are on contigs no run called")
    for rel in ("stored", "p2", "p3"):
        _stops(con, rel)
    res = {"stored_files": [str(p) for p in stored], "p2_runs": [str(p) for p in a.p2],
           "p3_runs": [str(p) for p in a.p3], "p2_vs_stored": compare(con, "p2"), "p3_vs_stored": compare(con, "p3"),
           "cost": costs(con)}
    for d in (*a.p2, *a.p3):
        res.setdefault("done", {})[str(d)] = json.loads((d / "done.json").read_text())
    # the share of bases under 100 kb in prokaryote genomes of 100 kb or more: the meta-mode CPU p3 adds
    res["bases_by_class"] = dict(con.execute("SELECT cls, sum(contig_len) FROM p2_contigs GROUP BY cls ORDER BY cls"
                                             ).fetchall())
    tmp = a.out.with_name(a.out.name + ".tmp")
    tmp.write_text(json.dumps(res, indent=1, sort_keys=True, default=str) + "\n")
    tmp.replace(a.out)

    def table(title: str, d: dict) -> None:
        print(f"\n### {title}\n")
        print("| class | contigs | changed contigs | stored genes | re-called | identical | start only | lost | gained |")
        print("|---|---|---|---|---|---|---|---|---|")
        for cls in ORDER:
            if cls not in d:
                continue
            r = d[cls]
            pct = {k: "n/a" if r[f"{k}_frac"] is None else f"{100 * r[f'{k}_frac']:.2f}%"
                   for k in ("identical", "start_only", "lost", "gained")}
            print(f"| {cls} | {r['contigs']:,} | {r['contigs_changed']:,} | {r['stored']:,} | {r['recalled']:,} | "
                  f"{r['identical']:,} ({pct['identical']}) | {r['start_only']:,} ({pct['start_only']}) | "
                  f"{r['lost']:,} ({pct['lost']}) | {r['gained']:,} ({pct['gained']}) |")

    table("p2 re-call against the stored calls (control)", res["p2_vs_stored"])
    table("p3 against the stored calls", res["p3_vs_stored"])
    print("\n### CPU by class (thread CPU seconds: training + finding)\n")
    print("| class | bases | p2 modes | p2 train s | p2 find s | p3 modes | p3 train s | p3 find s | p3 / p2 |")
    print("|---|---|---|---|---|---|---|---|---|")
    for cls in ORDER:
        if cls not in res["cost"]:
            continue
        p2, p3 = res["cost"][cls]["p2"], res["cost"][cls]["p3"]
        ratio = f"{p3['cpu_s'] / p2['cpu_s']:.2f}" if p2["cpu_s"] else "n/a"
        print(f"| {cls} | {p2['bases']:,} | {p2['modes']} | {p2['train_cpu_s']:.1f} | {p2['find_cpu_s']:.1f} | "
              f"{p3['modes']} | {p3['train_cpu_s']:.1f} | {p3['find_cpu_s']:.1f} | {ratio} |")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
