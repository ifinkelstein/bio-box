"""decoys.py WORLD...: per decoy of a kept micro backbone, replayed from ZOMBI2's genome_events.tsv (as
tests/r2/test_synthworld_backbone.py replays it): transfer rows, loss rows, the time its last copy died (or 'alive'),
its partner's (module first family) live-lineage count at the first speciation after the rate onset, and the time from
which at most one lineage holds the partner. Card F-02 session 12, scratch, measurement only."""
import collections
import csv
import sys
from pathlib import Path


def tsv(p):
    with open(p, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


for w in map(Path, sys.argv[1:]):
    fams = tsv(w / "backbone" / "families.tsv")
    name = {r["name"]: int(r["family"]) for r in fams if r["name"]}
    dec = {name[f"dec{i}"]: i for i in range(5)}
    partner = {int(r["family"]): int(r["partner"]) for r in fams if r["kind"] == "decoy"}
    watch = set(dec) | set(partner.values())
    live = {f: collections.Counter() for f in watch}
    stats = {f: {"tr": 0, "loss": 0, "dead": None} for f in dec}
    holders_le1 = {f: None for f in partner.values()}
    for r in tsv(w / "zombi2" / "genome_events.tsv"):
        f = int(r["family"])
        if f not in watch:
            continue
        t = float(r["time"])
        kids = [c for c in r["children"].split(";") if c]
        if f in dec:
            if r["kind"].startswith("transfer"):
                stats[f]["tr"] += 1
            elif r["kind"] == "loss":
                stats[f]["loss"] += 1
        if r["parents"]:
            live[f][r["parents"].split("_")[0]] -= 1
        for k in kids:
            live[f][k.split("_")[0]] += 1
        n = sum(1 for v in live[f].values() if v > 0)
        if f in dec and n == 0 and stats[f]["dead"] is None:
            stats[f]["dead"] = round(t, 3)
        if f in holders_le1:
            if n <= 1 and holders_le1[f] is None:
                holders_le1[f] = round(t, 3)
            elif n > 1:
                holders_le1[f] = None
    for f, i in sorted(dec.items(), key=lambda x: x[1]):
        p = partner[f]
        print(w.name, f"dec{i}", f"transfers {stats[f]['tr']}", f"losses {stats[f]['loss']}",
              f"extinct at {stats[f]['dead'] if stats[f]['dead'] is not None else 'alive'}",
              f"partner {p} held by <= 1 lineage from {holders_le1[p] if holders_le1[p] is not None else 'never'}")
