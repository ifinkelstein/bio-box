"""summ.py DIR LABEL...: one Markdown row per scan (card F-02 session 12, scratch; session 11's vsumm.py columns plus
the DEC-44 columns: seeds whose fewest decoy transfers are >= 3, the fewest over the seeds and at the default seed;
seeds with some genome without a plasmid; both types at the default seed)."""
import csv
import sys
from pathlib import Path


def med(v):
    v = sorted(v)
    return (v[(len(v) - 1) // 2] + v[len(v) // 2]) / 2 if v else float("nan")


print("| scan | chromosome mean (in band) | default mean | planted share % | copies per present member | exclusion "
      "breaches | both plasmid types (default) | some genome without a plasmid | a genome > 2 plasmids | all cells powered | "
      "levels pooled (increasing) | core >= 40 genomes: default, range (>= 100) | decoys >= 3 transfers: seeds, min, "
      "default | N0 ratio range (in [0.9, 1.1]) | N1 unmatched planted at default |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
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
    dec = [int(r["decoy_transfers_min"]) for r in rows]
    n0 = "-"
    if df.get("n0_bg_ratio"):
        q = [float(r["n0_bg_ratio"]) for r in rows]
        n0 = f"{min(q):.3f}-{max(q):.3f} ({sum(0.9 <= x <= 1.1 for x in q)}/{len(q)})"
    print(f"| {label} | {min(cm):.0f}-{max(cm):.0f} ({sum(1050 <= x <= 1950 for x in cm)}/{len(cm)}) | "
          f"{float(df['chrom_mean']):.0f} | {min(sh):.1f}-{max(sh):.1f} | {min(cp):.2f}-{max(cp):.2f} | {br} | "
          f"{both}/{len(rows)} ({'yes' if bdf else 'no'}) | {without}/{len(rows)} | {over} | {pw}/{len(rows)} | "
          f"{', '.join(f'{m:.1f}' for m in pooled)} ({inc}/{len(rows)}) | {int(df['root_families_ge40_genomes'])}, "
          f"{min(core)}-{max(core)} ({sum(c >= 100 for c in core)}) | {sum(x >= 3 for x in dec)}/{len(rows)}, "
          f"{min(dec)}, {df['decoy_transfers_min']} | {n0} | {df.get('n1_planted_without_match', '-')} |")
