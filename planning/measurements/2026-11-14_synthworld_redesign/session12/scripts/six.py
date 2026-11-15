"""six.py DIR LABEL: the per-seed test's view of a scan (card F-02 session 12, scratch): pooled medians of realized
transfers per member by level over the test's six seeds (default, 20261113, the first four N2 seeds) and their ratios,
seeds whose top level median does not exceed the bottom's, and cells without a module whole in >= 9 genomes."""
import csv
import sys
from pathlib import Path

d, label = Path(sys.argv[1]), sys.argv[2]
rows = list(csv.DictReader(open(d / f"{label}.tsv"), delimiter="\t"))
six = [r["seed"] for r in rows[:6]]
grid = [g for g in csv.DictReader(open(d / f"{label}_grid.tsv"), delimiter="\t") if g["seed"] in six]


def med(v):
    v = sorted(v)
    return (v[(len(v) - 1) // 2] + v[len(v) // 2]) / 2


lv = sorted({float(g["member_transfer"]) for g in grid})
pooled = [med([float(g["transfers_per_member"]) for g in grid if float(g["member_transfer"]) == x]) for x in lv]
print("seeds", six)
print("pooled", pooled, "ratios", [round(b / a, 3) for a, b in zip(pooled, pooled[1:])])
print("top <= bottom at", [r["seed"] for r in rows[:6] if not float(r["tpm_level2_median"]) > float(r["tpm_level0_median"])])
print("unpowered cells", {r["seed"]: int(r["cells_level0_unpowered"]) + int(r["cells_level1_unpowered"]) +
                          int(r["cells_level2_unpowered"]) for r in rows[:6]})
