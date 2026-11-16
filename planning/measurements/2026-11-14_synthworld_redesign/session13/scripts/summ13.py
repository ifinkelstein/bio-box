"""summ13.py DIR LABEL...: one Markdown row per scan (card F-02 session 13, scratch; session 12's summ.py columns
plus the DEC-48 (1) decoy columns read from the decoy test's replay in the scan: seeds whose every decoy has >= 3
transfers into a partner carrier, the fewest over the seeds and at the default seed; transfers off the partner; the
partner's mean share of the other lineages (the test asks <= 0.9) over the seeds and at the default; decoy genes per
genome and copies per genome where present; and, with --n0, the root-presence KS of both N0 worlds)."""
import csv
import sys
from pathlib import Path


def med(v):
    v = sorted(v)
    return (v[(len(v) - 1) // 2] + v[len(v) // 2]) / 2 if v else float("nan")


cols = ["scan", "chromosome mean (in band)", "planted share %", "copies per present member", "exclusion breaches",
        "both plasmid types (default)", "some genome without a plasmid", "a genome > 2 plasmids", "all cells powered",
        "levels pooled (increasing)", "core >= 40 genomes: default, range (>= 100)", "N1 unmatched at default",
        "decoys >= 3 transfers: seeds, fewest, default", "off partner", "partner share: max, default (<= 0.9)",
        "decoy genes per genome: range, default", "decoy copies per present genome: range, max",
        "N0 ratio (in [0.9, 1.1])", "root KS full / strict below bound; default D full, strict; 20261113",
        "per-seed test's six seeds: pooled levels (steps >= 1.5), most unpowered cells"]
print("| " + " | ".join(cols) + " |")
print("|" + "---|" * len(cols))
for label in sys.argv[2:]:
    d = Path(sys.argv[1])
    rows = list(csv.DictReader(open(d / f"{label}.tsv"), delimiter="\t"))
    grid = list(csv.DictReader(open(d / f"{label}_grid.tsv"), delimiter="\t"))
    df = next(r for r in rows if r["seed"] == "20261112")
    cm = [float(r["chrom_mean"]) for r in rows]
    sh = [100 * float(r["share_planted_all"]) for r in rows]
    cp = [float(r["module_copies_mean"]) for r in rows]
    br = max(int(r["module_transfers_into_first_holder"]) for r in rows)
    both = sum(int(r["pls0_genomes"]) > 0 and int(r["pls1_genomes"]) > 0 for r in rows)
    bdf = int(df["pls0_genomes"]) > 0 and int(df["pls1_genomes"]) > 0
    over = without = 0
    for r in rows:
        z, o, t = (int(x) for x in r["genomes_0_1_2_plasmids"].split("/"))
        over += int(r["genomes"]) - (z + o + t) > 0
        without += z > 0
    pw = sum(int(r["cells_level0_unpowered"]) + int(r["cells_level1_unpowered"]) + int(r["cells_level2_unpowered"])
             == 0 for r in rows)
    lv = sorted({float(x["member_transfer"]) for x in grid})
    pooled = [med([float(x["transfers_per_member"]) for x in grid if float(x["member_transfer"]) == v]) for v in lv]
    inc = sum(float(r["tpm_level0_median"]) < float(r["tpm_level1_median"]) < float(r["tpm_level2_median"])
              for r in rows)
    core = [int(r["root_families_ge40_genomes"]) for r in rows]
    per = [[int(x) for x in r["decoy_transfers_per_decoy"].split("/")] for r in rows]
    fewest = [min(p) for p in per]
    fd = min(int(x) for x in df["decoy_transfers_per_decoy"].split("/"))
    off = sum(int(r["decoy_transfers_off_partner"]) for r in rows)
    shares = [float(r["decoy_partner_share_mean"]) for r in rows if r["decoy_partner_share_mean"] != ""]
    dg = [float(r["per_genome_decoy"]) for r in rows]
    dc = [float(r["decoy_copies_per_present_genome"]) for r in rows]
    dmax = max(int(r["decoy_copies_max"]) for r in rows)
    n0 = ks = "-"
    if df.get("n0_bg_ratio"):
        q = [float(r["n0_bg_ratio"]) for r in rows]
        n0 = f"{min(q):.3f}-{max(q):.3f} ({sum(0.9 <= x <= 1.1 for x in q)}/{len(q)})"
        okf = sum(float(r["n0full_root_ks_d"]) < float(r["n0full_root_ks_crit"]) for r in rows)
        oks = sum(float(r["n0strict_root_ks_d"]) < float(r["n0strict_root_ks_crit"]) for r in rows)
        s2 = next(r for r in rows if r["seed"] == "20261113")
        ks = (f"{okf}/{len(rows)} / {oks}/{len(rows)}; {df['n0full_root_ks_d']}, {df['n0strict_root_ks_d']} "
              f"({df['n0full_root_ks_crit']}); {s2['n0full_root_ks_d']}, {s2['n0strict_root_ks_d']} "
              f"({s2['n0full_root_ks_crit']})")
    six = [r["seed"] for r in rows[:6]]          # default, 20261113, the first four N2 seeds (the scan's order)
    p6 = [med([float(x["transfers_per_member"]) for x in grid if x["seed"] in six and float(x["member_transfer"]) == v])
          for v in lv]
    unp6 = max(int(r["cells_level0_unpowered"]) + int(r["cells_level1_unpowered"]) + int(r["cells_level2_unpowered"])
               for r in rows[:6])
    sixcol = f"{', '.join(f'{m:.2f}' for m in p6)} ({', '.join(f'{b / a:.2f}x' for a, b in zip(p6, p6[1:]))}), {unp6}"
    print(f"| {label} | {min(cm):.0f}-{max(cm):.0f} ({sum(1050 <= x <= 1950 for x in cm)}/{len(cm)}) | "
          f"{min(sh):.1f}-{max(sh):.1f} | {min(cp):.2f}-{max(cp):.2f} | {br} | "
          f"{both}/{len(rows)} ({'yes' if bdf else 'no'}) | {without}/{len(rows)} | {over} | {pw}/{len(rows)} | "
          f"{', '.join(f'{m:.1f}' for m in pooled)} ({inc}/{len(rows)}) | {int(df['root_families_ge40_genomes'])}, "
          f"{min(core)}-{max(core)} ({sum(c >= 100 for c in core)}) | {df.get('n1_planted_without_match', '-')} | "
          f"{sum(x >= 3 for x in fewest)}/{len(rows)}, {min(fewest)}, {fd} | {off} | "
          f"{max(shares):.3f}, {df['decoy_partner_share_mean']} | {min(dg):.1f}-{max(dg):.1f}, {df['per_genome_decoy']} | "
          f"{min(dc):.2f}-{max(dc):.2f}, {dmax} | {n0} | {ks} | {sixcol} |")
    bad = {r["seed"]: r["decoy_transfers_per_decoy"] for r in rows if min(int(x) for x in r["decoy_transfers_per_decoy"].split("/")) < 3}
    if bad:
        print(f"<!-- {label}: seeds with a decoy below 3 transfers: {bad} -->")
