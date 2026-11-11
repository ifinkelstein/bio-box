"""Compare MMseqs2 profile-search hits with the producer's hmmsearch domains (card 0-09a).

    python compare.py <export_dir> <models.tsv> <hits.tsv> <out.json> [--label L]

<export_dir> is export_stage.py's output (proteins.tsv, truth.tsv); <models.tsv> is pfam_models.py's;
<hits.tsv> is mm_search.sh's convertalis output (query uid, target accession, qstart, qend, tstart, tend,
evalue, bits; 1-based inclusive query coordinates). Everything is computed from these files in DuckDB;
nothing here calls MMseqs2 or HMMER.

Definitions (the measurement file repeats them):
- truth pair: a (protein, model) with at least one row in the producer's l1_domains (GA cutoffs, Z=1e8).
  Its margin is seq_bits - GA_seq (the sequence-level cutoff that admits the pair).
- MMseqs2 pair: a (protein, model) with at least one MMseqs2 alignment; its score is the best bits, its
  E-value the smallest.
- At threshold t, a pair is called when its bits >= t. recall = called truth pairs / truth pairs,
  precision = called truth pairs / called pairs. The curve runs over every distinct bits value.
- Boundaries: for a called truth pair, MMseqs2's best alignment (0-based half-open [qstart-1, qend)) against
  the pair's truth domain that overlaps it most; absolute start and end differences, against env and ali.
- Domain recall: a truth domain row is found when the pair's best MMseqs2 alignment overlaps it.
- Prefilter (derived, not run): proteins with any MMseqs2 hit at E <= e go to hmmsearch against every
  model, work = their residues / all residues, ceiling = truth pairs of those proteins / truth pairs; or
  only the called (protein, model) pairs go, work = sum of protein length x model length over them /
  (all residues x all model lengths), ceiling = called truth pairs / truth pairs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb
import numpy as np

MARGIN_BINS = [0, 1, 2, 5, 10, 20, 50, 1e9]
EVALUES = [1e-10, 1e-5, 1e-3, 1e-2, 0.1, 1, 10, 100]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("export", type=Path)
    ap.add_argument("models", type=Path)
    ap.add_argument("hits", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--label", default="")
    a = ap.parse_args()
    con = duckdb.connect()
    con.execute("SET threads=4")
    con.execute("SET memory_limit='20GB'")
    con.execute(f"""CREATE TABLE models AS SELECT acc, name, leng::INT leng, ga_seq::DOUBLE ga_seq,
        ga_dom::DOUBLE ga_dom FROM read_csv('{a.models}', delim='\t', header=true, all_varchar=true)""")
    con.execute(f"""CREATE TABLE prot AS SELECT uid, length::BIGINT length
        FROM read_csv('{a.export}/proteins.tsv', delim='\t', header=true, all_varchar=true)""")
    con.execute(f"""CREATE TABLE truth AS SELECT uid, acc, env_start::INT env_start, env_end::INT env_end,
        ali_start::INT ali_start, ali_end::INT ali_end, seq_bits::DOUBLE seq_bits, dom_bits::DOUBLE dom_bits
        FROM read_csv('{a.export}/truth.tsv', delim='\t', header=true, all_varchar=true)""")
    con.execute(f"""CREATE TABLE mm AS SELECT column0 uid, column1 acc, column2::INT - 1 qs, column3::INT qe,
        column6::DOUBLE evalue, column7::DOUBLE bits
        FROM read_csv('{a.hits}', delim='\t', header=false, all_varchar=true)""")
    unknown = con.execute("SELECT count(*) FROM (SELECT DISTINCT acc FROM mm) m ANTI JOIN models USING (acc)").fetchone()[0]
    if unknown:
        raise SystemExit(f"{unknown} MMseqs2 target accessions are not in {a.models}")
    con.execute("""CREATE TABLE tp AS SELECT t.uid, t.acc, max(seq_bits) seq_bits, count(*) n_dom,
        max(seq_bits) - any_value(m.ga_seq) margin FROM truth t JOIN models m USING (acc) GROUP BY t.uid, t.acc""")
    con.execute("""CREATE TABLE mp AS SELECT uid, acc, max(bits) bits, min(evalue) evalue,
        arg_max(qs, bits) qs, arg_max(qe, bits) qe FROM mm GROUP BY uid, acc""")
    con.execute("""CREATE TABLE j AS SELECT mp.*, (tp.uid IS NOT NULL) is_true, tp.margin
        FROM mp LEFT JOIN tp USING (uid, acc)""")
    n_truth = con.execute("SELECT count(*) FROM tp").fetchone()[0]
    n_rows = con.execute("SELECT count(*) FROM truth").fetchone()[0]
    n_prot, residues = con.execute("SELECT count(*), sum(length) FROM prot").fetchone()
    n_models, model_len = con.execute("SELECT count(*), sum(leng) FROM models").fetchone()
    n_mm_rows = con.execute("SELECT count(*) FROM mm").fetchone()[0]
    n_mm_pairs = con.execute("SELECT count(*) FROM mp").fetchone()[0]

    # the curve: called pairs and called truth pairs at each distinct bits value, descending
    cur = con.execute("""SELECT bits, sum(is_true::INT) tp, count(*) n FROM j GROUP BY bits ORDER BY bits DESC""").fetchnumpy()
    bits = cur["bits"]
    ctp = np.cumsum(cur["tp"])
    cn = np.cumsum(cur["n"])
    recall = ctp / n_truth
    precision = ctp / cn
    f1 = np.where(ctp > 0, 2 * precision * recall / (precision + recall), 0)
    best = int(np.argmax(f1))
    hi = np.where(precision >= 0.99)[0]
    r99 = (float(recall[hi].max()), float(bits[hi[np.argmax(recall[hi])]])) if hi.size else (0.0, None)
    hi95 = np.where(precision >= 0.95)[0]
    r95 = (float(recall[hi95].max()), float(bits[hi95[np.argmax(recall[hi95])]])) if hi95.size else (0.0, None)
    curve_pts = []
    for t in [15, 20, 25, 30, 35, 40, 50, 60, 80, 100]:
        k = np.searchsorted(-bits, -t, side="right") - 1
        if k >= 0:
            curve_pts.append({"bits_ge": t, "recall": float(recall[k]), "precision": float(precision[k]),
                              "called": int(cn[k])})
    ev_pts = []
    for e in EVALUES:
        c, t = con.execute(f"SELECT count(*), coalesce(sum(is_true::INT), 0) FROM j WHERE evalue <= {e}").fetchone()
        ev_pts.append({"evalue_le": e, "recall": t / n_truth, "precision": (t / c) if c else None, "called": c})
    t_best = float(bits[best])

    def margins(where: str) -> list[dict]:
        out = []
        for lo, hi_ in zip(MARGIN_BINS, MARGIN_BINS[1:]):
            n, f = con.execute(f"""SELECT count(*), coalesce(sum((j.uid IS NOT NULL)::INT), 0) FROM tp
                LEFT JOIN (SELECT uid, acc FROM j WHERE {where}) j USING (uid, acc)
                WHERE tp.margin >= {lo} AND tp.margin < {hi_}""").fetchone()
            out.append({"margin_lo": lo, "margin_hi": hi_ if hi_ < 1e9 else None, "truth_pairs": n,
                        "recall": (f / n) if n else None})
        return out

    margin_tables = {"any_hit": margins("true"), "evalue_le_1e-3": margins("evalue <= 1e-3"),
                     f"bits_ge_best_f1_{t_best}": margins(f"bits >= {t_best}")}
    neg = con.execute("SELECT count(*) FROM tp WHERE margin < 0").fetchone()[0]

    # boundaries of called truth pairs (any hit), against the overlapping truth domain
    bnd = con.execute("""WITH c AS (SELECT j.uid, j.acc, j.qs, j.qe, t.env_start, t.env_end, t.ali_start, t.ali_end,
            least(j.qe, t.env_end) - greatest(j.qs, t.env_start) ov FROM j JOIN truth t USING (uid, acc) WHERE j.is_true
            QUALIFY row_number() OVER (PARTITION BY j.uid, j.acc ORDER BY ov DESC, t.env_start) = 1)
        SELECT count(*) n, median(abs(qs - env_start)), median(abs(qe - env_end)), quantile_cont(abs(qs - env_start), 0.9),
            quantile_cont(abs(qe - env_end), 0.9), median(abs(qs - ali_start)), median(abs(qe - ali_end)),
            avg((ov > 0)::INT), median(qs - env_start), median(qe - env_end),
            median((qe - qs)::DOUBLE / (env_end - env_start))
        FROM c""").fetchone()
    boundaries = dict(zip(["pairs", "median_abs_start_env", "median_abs_end_env", "p90_abs_start_env",
                           "p90_abs_end_env", "median_abs_start_ali", "median_abs_end_ali", "fraction_overlapping",
                           "median_signed_start_env", "median_signed_end_env", "median_len_ratio_mm_over_env"], bnd))
    dom_found = con.execute("""SELECT avg((j.uid IS NOT NULL AND j.qs < t.env_end AND t.env_start < j.qe)::INT)
        FROM truth t LEFT JOIN j USING (uid, acc)""").fetchone()[0]
    multi = con.execute("SELECT count(*), sum(n_dom) FROM tp WHERE n_dom > 1").fetchone()

    # called pairs that are not truth pairs: do they sit on a truth domain of another model?
    fp = con.execute(f"""WITH f AS (SELECT * FROM j WHERE NOT is_true AND bits >= {t_best})
        SELECT count(*), coalesce(sum(o::INT), 0), coalesce(sum(np_::INT), 0) FROM (
          SELECT f.uid, f.acc, bool_or(t.uid IS NOT NULL AND f.qs < t.env_end AND t.env_start < f.qe) o,
                 bool_and(t.uid IS NULL) np_
          FROM f LEFT JOIN truth t ON t.uid = f.uid GROUP BY f.uid, f.acc)""").fetchone()

    # family-agnostic regions (a stand-in for clan collapsing, which needs Pfam-A.clans.tsv, not on the pod):
    # a truth domain is covered when an MMseqs2 alignment of ANY model on its protein overlaps it by at least
    # half of the shorter of the two; an MMseqs2 alignment is on a truth region when it so overlaps any truth
    # domain of its protein
    def regions(where: str) -> tuple[float, float]:
        rec = con.execute(f"""SELECT avg(c::INT) FROM (SELECT t.uid, t.env_start, t.acc,
            bool_or(m.uid IS NOT NULL AND least(m.qe, t.env_end) - greatest(m.qs, t.env_start)
                    >= 0.5 * least(m.qe - m.qs, t.env_end - t.env_start)) c
            FROM truth t LEFT JOIN (SELECT * FROM mm WHERE {where}) m ON m.uid = t.uid
            GROUP BY t.uid, t.env_start, t.acc)""").fetchone()[0]
        prec = con.execute(f"""SELECT avg(c::INT) FROM (SELECT m.uid, m.acc, m.qs,
            bool_or(t.uid IS NOT NULL AND least(m.qe, t.env_end) - greatest(m.qs, t.env_start)
                    >= 0.5 * least(m.qe - m.qs, t.env_end - t.env_start)) c
            FROM (SELECT * FROM mm WHERE {where}) m LEFT JOIN truth t ON m.uid = t.uid
            GROUP BY m.uid, m.acc, m.qs)""").fetchone()[0]
        return rec, prec
    region = {}
    for name, where in [("any_hit", "true"), (f"bits_ge_{t_best}", f"bits >= {t_best}"),
                        ("evalue_le_1e-3", "evalue <= 1e-3")]:
        r, p = regions(where)
        region[name] = {"truth_domains_covered": r, "alignments_on_truth_region": p}

    # prefilter (derived)
    pref = []
    for e in [100, 10, 1, 0.1, 1e-3]:
        sel_res, sel_n = con.execute(f"""SELECT coalesce(sum(length), 0), count(*) FROM prot SEMI JOIN
            (SELECT uid FROM mp WHERE evalue <= {e}) s USING (uid)""").fetchone()
        ceil_p = con.execute(f"""SELECT count(*) FROM tp SEMI JOIN (SELECT uid FROM mp WHERE evalue <= {e}) s
            USING (uid)""").fetchone()[0]
        pair_work, pair_n, ceil_pair = con.execute(f"""SELECT coalesce(sum(p.length * m.leng), 0), count(*),
            coalesce(sum(j.is_true::INT), 0) FROM j JOIN prot p USING (uid) JOIN models m USING (acc)
            WHERE j.evalue <= {e}""").fetchone()
        pref.append({"evalue_le": e, "proteins_selected": sel_n, "protein_work_fraction": sel_res / residues,
                     "protein_recall_ceiling": ceil_p / n_truth, "pairs_selected": pair_n,
                     "pair_work_fraction": pair_work / (residues * model_len),
                     "pair_recall_ceiling": ceil_pair / n_truth})
    prot_with_truth, res_with_truth = con.execute("""SELECT count(*), sum(length) FROM prot SEMI JOIN tp USING (uid)""").fetchone()

    doc = {
        "label": a.label, "hits": str(a.hits), "export": str(a.export),
        "proteins": n_prot, "residues": int(residues), "models": n_models, "model_match_states": int(model_len),
        "truth_rows": n_rows, "truth_pairs": n_truth, "truth_pairs_negative_margin": neg,
        "truth_pairs_multi_domain": int(multi[0]), "truth_proteins": prot_with_truth,
        "truth_protein_residue_fraction": res_with_truth / residues,
        "mm_rows": n_mm_rows, "mm_pairs": n_mm_pairs,
        "best_f1": {"f1": float(f1[best]), "bits_ge": t_best, "recall": float(recall[best]),
                    "precision": float(precision[best]), "called": int(cn[best])},
        "recall_at_precision_0.99": {"recall": r99[0], "bits_ge": r99[1]},
        "recall_at_precision_0.95": {"recall": r95[0], "bits_ge": r95[1]},
        "recall_any_hit": float(recall[-1]), "precision_any_hit": float(precision[-1]),
        "curve_bits": curve_pts, "curve_evalue": ev_pts, "recall_by_margin": margin_tables,
        "boundaries": boundaries, "domain_rows_found_any_hit": dom_found,
        "false_pairs_at_best_f1": {"pairs": fp[0], "overlapping_a_truth_domain_of_another_model": fp[1],
                                   "on_proteins_without_truth": fp[2]},
        "prefilter": pref, "regions_any_model": region,
    }
    a.out.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: doc[k] for k in ("label", "truth_pairs", "mm_pairs", "best_f1", "recall_at_precision_0.99",
                                           "recall_any_hit")}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
