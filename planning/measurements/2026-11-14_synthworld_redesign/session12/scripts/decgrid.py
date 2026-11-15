"""decgrid.py DIR LABEL: per seed, the fewest decoy transfers and the partner modules' genomes with every member and
realized transfers per member (from the scan's grid). Card F-02 session 12, scratch."""
import csv
import sys
from pathlib import Path

d, label = Path(sys.argv[1]), sys.argv[2]
partners = ["mod01", "mod09", "mod11", "mod19", "mod22"]
rows = list(csv.DictReader(open(d / f"{label}.tsv"), delimiter="\t"))
grid = list(csv.DictReader(open(d / f"{label}_grid.tsv"), delimiter="\t"))
for r in rows:
    g = {x["module"]: x for x in grid if x["seed"] == r["seed"]}
    print(r["seed"], r["decoy_transfers_min"],
          " ".join(f"{m}:{g[m]['genomes_all_members']}g/{g[m]['transfers_per_member']}t" for m in partners))
