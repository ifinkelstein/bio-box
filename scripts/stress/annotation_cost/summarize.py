"""Markdown tables of card 0-09a from compare.py's JSON files (host python3, reads only).

    python3 summarize.py <compare.json>...
"""

from __future__ import annotations

import json
import sys


def pct(x) -> str:
    return "-" if x is None else f"{100 * x:.1f}%"


def main() -> int:
    docs = [json.load(open(p)) for p in sys.argv[1:]]
    print("| run | truth pairs | MMseqs2 pairs | recall, any hit | best F1 (bits >=, P, R) | R at P>=0.99 (bits >=) "
          "| R at P>=0.95 | E<=1e-3: R, P | E<=1: R, P |")
    print("|---|---|---|---|---|---|---|---|---|")
    for d in docs:
        b, r99, r95 = d["best_f1"], d["recall_at_precision_0.99"], d["recall_at_precision_0.95"]
        ev = {e["evalue_le"]: e for e in d["curve_evalue"]}
        print(f"| {d['label']} | {d['truth_pairs']:,} | {d['mm_pairs']:,} | {pct(d['recall_any_hit'])} | "
              f"{b['f1']:.3f} ({b['bits_ge']:g}, {pct(b['precision'])}, {pct(b['recall'])}) | "
              f"{pct(r99['recall'])} ({r99['bits_ge']}) | {pct(r95['recall'])} | "
              f"{pct(ev[0.001]['recall'])}, {pct(ev[0.001]['precision'])} | {pct(ev[1]['recall'])}, {pct(ev[1]['precision'])} |")
    for d in docs:
        print(f"\n**{d['label']}**: recall by truth margin (seq_bits - GA_seq, bits)\n")
        keys = list(d["recall_by_margin"])
        print("| margin | truth pairs | " + " | ".join(keys) + " |")
        print("|---|---|" + "---|" * len(keys))
        rows = d["recall_by_margin"][keys[0]]
        for i, row in enumerate(rows):
            hi = row["margin_hi"]
            label = f"[{row['margin_lo']:g}, {hi:g})" if hi is not None else f">= {row['margin_lo']:g}"
            print(f"| {label} | {row['truth_pairs']:,} | " +
                  " | ".join(pct(d["recall_by_margin"][k][i]["recall"]) for k in keys) + " |")
        print("\ncurve (bits >= t): " + "; ".join(f"{c['bits_ge']}: R {pct(c['recall'])} P {pct(c['precision'])}"
                                                   for c in d["curve_bits"]))
        print("\ncurve (E <= e): " + "; ".join(f"{c['evalue_le']:g}: R {pct(c['recall'])} P {pct(c['precision'])}"
                                                 for c in d["curve_evalue"]))
        bd = d["boundaries"]
        print(f"\nboundaries ({bd['pairs']:,} called truth pairs): median |start-env_start| {bd['median_abs_start_env']}, "
              f"median |end-env_end| {bd['median_abs_end_env']}, p90 {bd['p90_abs_start_env']} / {bd['p90_abs_end_env']}; "
              f"vs ali: {bd['median_abs_start_ali']} / {bd['median_abs_end_ali']}; signed median start {bd['median_signed_start_env']}, "
              f"end {bd['median_signed_end_env']}; median length ratio {bd['median_len_ratio_mm_over_env']:.3f}; "
              f"overlapping {pct(bd['fraction_overlapping'])}; truth domain rows overlapped by the pair's best alignment "
              f"{pct(d['domain_rows_found_any_hit'])} ({d['truth_pairs_multi_domain']:,} truth pairs have >1 domain)")
        f = d["false_pairs_at_best_f1"]
        print(f"\ncalled non-truth pairs at best F1: {f['pairs']:,}; on a truth domain of another model "
              f"{f['overlapping_a_truth_domain_of_another_model']:,}; on proteins with no truth domain {f['on_proteins_without_truth']:,}")
        if "regions_any_model" in d:
            print("\nregions, any model (truth domains covered / MMseqs2 alignments on a truth region): " + "; ".join(
                f"{k}: {pct(v['truth_domains_covered'])} / {pct(v['alignments_on_truth_region'])}"
                for k, v in d["regions_any_model"].items()))
        print("\n| prefilter E <= | proteins kept | hmmsearch work (proteins) | ceiling | pairs kept | work (pairs) | ceiling |")
        print("|---|---|---|---|---|---|---|")
        for p in d["prefilter"]:
            print(f"| {p['evalue_le']:g} | {p['proteins_selected']:,} | {pct(p['protein_work_fraction'])} | "
                  f"{pct(p['protein_recall_ceiling'])} | {p['pairs_selected']:,} | {100 * p['pair_work_fraction']:.4f}% | "
                  f"{pct(p['pair_recall_ceiling'])} |")
        print(f"\nproteins {d['proteins']:,}, residues {d['residues']:,}, with a truth domain {d['truth_proteins']:,} "
              f"({pct(d['truth_protein_residue_fraction'])} of residues); models {d['models']:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
