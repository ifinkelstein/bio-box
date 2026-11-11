"""Print the markdown tables of card 0-09b from the JSON files of `representatives.py`.

    python rep_tables.py <infer.json> <saving.json>

<infer.json> is the run with `--domains` (release rN.1b); <saving.json> the run without (release rN.1).
Every number printed is read from the two files or derived from them by the formula in the table's caption.
"""

from __future__ import annotations

import json
import sys

LEVELS = ("nr90", "nr50")
ORDER = {
    "by_pident": ["none", "[0.00,0.50)", "[0.50,0.60)", "[0.60,0.70)", "[0.70,0.80)", "[0.80,0.90)", "[0.90,0.95)",
                  "[0.95,0.99)", "[0.99,1.00)", "1.00"],
    "by_cluster_size": ["2-2", "3-4", "5-9", "10-99", ">=100"],
    "by_length": ["shorter_gt20", "within20", "longer_gt20"],
}


def pct(x, d=1):
    return "-" if x is None else f"{100 * x:.{d}f}%"


def n(x):
    return f"{x:,}"


def main(argv: list[str]) -> int:
    inf = json.load(open(argv[0]))
    sav = json.load(open(argv[1]))
    scopes = inf["stages"] + ["all"]

    print("### Decision table (rN.1b)\n")
    print("| level | scope | mode | residues scanned | pairs right (of true) | pairs wrong (of reported) | "
          "pairs missed (of true) | M residues not scanned | CPU-min saved at 7.44 / 7.96 |")
    print("|---|---|---|---|---|---|---|---|---|")
    for lv in LEVELS:
        for sc in scopes:
            for mode in ("reps_only", "hybrid"):
                d = inf["inference"][lv]["scopes"][sc]["decision"][mode]
                print(f"| {lv} | {sc} | {mode.replace('_', ' ')} | {pct(d['residue_fraction_scanned'])} | "
                      f"{pct(d['right_fraction_of_true'], 2)} | {pct(d['false_fraction_of_reported'], 2)} | "
                      f"{pct(d['missed_fraction_of_true'], 2)} | {d['residues_not_scanned_millions']:.1f} | "
                      f"{d['cpu_min_saved_at_7.44']:,.0f} / {d['cpu_min_saved_at_7.96']:,.0f} |")
    print()
    print("### Hybrid: what it adds (rN.1b)\n")
    print("| level | scope | members scanned too | residues added (of scope) | true pairs |")
    print("|---|---|---|---|---|")
    for lv in LEVELS:
        for sc in scopes:
            dd = inf["inference"][lv]["scopes"][sc]["decision"]
            add = dd["hybrid"]["residue_fraction_scanned"] - dd["reps_only"]["residue_fraction_scanned"]
            print(f"| {lv} | {sc} | {n(dd['hybrid']['members_scanned_too'])} | {pct(add, 2)} | {n(dd['hybrid']['true_pairs'])} |")
    print()

    for name, doc in (("rN.1b", inf), ("rN.1", sav)):
        print(f"### Saving ({name})\n")
        print("| level | scope | proteins | residues | representatives | clusters touched | not representatives: proteins | residues | "
              "members whose rep is in another stage | hybrid: residues scanned |")
        print("|---|---|---|---|---|---|---|---|---|---|")
        for lv in LEVELS:
            for sc in doc["stages"] + ["all"]:
                s = doc["saving"][f"{lv}/{sc}"]
                print(f"| {lv} | {sc} | {n(s['proteins'])} | {n(s['residues'])} | {n(s['reps'])} | {n(s['clusters_touched'])} | "
                      f"{pct(s['nonrep_protein_fraction'])} | {pct(s['nonrep_residue_fraction'])} | "
                      f"{n(s['members_whose_rep_is_outside_scope'])} | {pct(s['hybrid_residue_fraction_scanned'])} |")
        print()
    print("### rN.1: CPU-minutes not spent (*derived*: (residues - residues scanned) / 1e6 x rate)\n")
    print("| level | scope | mode | M residues not scanned | CPU-min at 7.44 / 7.96 |")
    print("|---|---|---|---|---|")
    for lv in LEVELS:
        for sc in sav["stages"] + ["all"]:
            s = sav["saving"][f"{lv}/{sc}"]
            for mode, scanned in (("reps only", s["rep_residues"]),
                                  ("hybrid", s["rep_residues"] + s["members_length_differs_gt20_residues"])):
                m = (s["residues"] - scanned) / 1e6
                print(f"| {lv} | {sc} | {mode} | {m:,.1f} | {m * 7.44:,.0f} / {m * 7.96:,.0f} |")
    print()

    print("### Members against their representative (rN.1b)\n")
    print("| level | scope | members | clusters | equal | subset | (of which member has none) | superset | (of which rep has none) | disjoint | overlap | "
          "equal, members with any domain | architecture equal | false of inferred pairs | missed of true pairs | "
          "per cluster: equal, false, missed |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for lv in LEVELS:
        for sc in scopes:
            s = inf["inference"][lv]["scopes"][sc]["all"]
            c = s["cluster_weighted"]
            print(f"| {lv} | {sc} | {n(s['members'])} | {n(s['clusters'])} | {pct(s['equal'])} | {pct(s['subset'], 2)} | "
                  f"{pct(s['member_empty_rep_not'], 2)} | {pct(s['superset'], 2)} | {pct(s['rep_empty_member_not'], 2)} | "
                  f"{pct(s['disjoint'], 2)} | {pct(s['overlap'], 2)} | {pct(s['equal_among_members_with_any_domain'])} | "
                  f"{pct(s['architecture_equal'])} | {pct(s['false_fraction_of_inferred'], 2)} | "
                  f"{pct(s['missed_fraction_of_true'], 2)} | {pct(c['equal'])}, {pct(c['false_fraction_of_inferred'], 2)}, "
                  f"{pct(c['missed_fraction_of_true'], 2)} |")
    print()

    for key, title in (("by_pident", "identity to the representative (`cluster_member_aln.pident`)"),
                       ("by_cluster_size", "cluster size (all members, representative included)"),
                       ("by_length", "member length against the representative's")):
        for sc in scopes:
            print(f"### By {title}, scope {sc} (rN.1b)\n")
            print("| level | bin | members | equal | subset | superset | disjoint or overlap | false of inferred | missed of true |")
            print("|---|---|---|---|---|---|---|---|---|")
            for lv in LEVELS:
                bins = inf["inference"][lv]["scopes"][sc][key]
                for b in ORDER[key]:
                    s = bins.get(b, {"members": 0})
                    if not s["members"]:
                        continue
                    print(f"| {lv} | {b} | {n(s['members'])} | {pct(s['equal'])} | {pct(s['subset'], 2)} | "
                          f"{pct(s['superset'], 2)} | {pct(s['disjoint'] + s['overlap'], 2)} | "
                          f"{pct(s['false_fraction_of_inferred'], 2)} | {pct(s['missed_fraction_of_true'], 2)} |")
            print()

    print("### Coordinates where the sets agree (rN.1b)\n")
    print("Residues: |start| and |end| differences, member's own envelope against the representative's mapped "
          "through the alignment, and copied without mapping. Fractions: |s_m/L_m - s_r/L_r| and the same for ends.\n")
    print("| level | scope | subset | domains | mapped start: median, p95 | mapped end: median, p95 | both ends within 5 | "
          "copied start: median, p95 | copied end: median, p95 | fraction start: median, p95 | fraction end: median, p95 |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for lv in LEVELS:
        for sc in scopes:
            for sub, s in inf["inference"][lv]["scopes"][sc]["coordinates"].items():
                if not s["ds"]["n"]:
                    continue
                q = lambda k, f="{:.0f}": f"{f.format(s[k]['median'])}, {f.format(s[k]['p95'])}"
                print(f"| {lv} | {sc} | {sub.replace('_', ' ')} | {n(s['ds']['n'])} | {q('ds')} | {q('de')} | "
                      f"{pct(s['both_ends_within_5_mapped'])} | {q('cs_')} | {q('ce')} | {q('fs', '{:.4f}')} | {q('fe', '{:.4f}')} |")
    print()
    print("### Checks\n")
    print(f"- proteins rN.1b: {inf['proteins']}")
    print(f"- proteins rN.1: {sav['proteins']}")
    print(f"- domains: {inf['domains']}")
    for lv in LEVELS:
        L = inf["inference"][lv]
        print(f"- {lv}: no alignment {L['no_alignment']}, CIGAR span mismatch {L['cigar_span_mismatch']}, "
              f"members longer than their rep {L['members_longer_than_rep']}; alignment check {inf['alignment_check'][lv]}")
    print(f"- elapsed: {inf['elapsed_s']} s and {sav['elapsed_s']} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
