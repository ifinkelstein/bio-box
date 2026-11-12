"""W-25 real run: summary of `nbr flanks` outputs (groups, context types of the largest groups, phage/prokaryote sharing).

python analyse.py OUT_DIR [OUT_DIR ...]
"""
import csv
import sys
from collections import Counter, defaultdict


def rows(path):
    with open(path, newline="") as fh:
        lines = [x for x in fh if not x.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t", quoting=csv.QUOTE_NONE))


for out in sys.argv[1:]:
    groups = rows(f"{out}/groups.tsv")
    ctx = rows(f"{out}/context_types.tsv")
    print(f"== {out}")
    print("groups", len(groups), "loci", len(ctx), "in_tree", sum(r["in_tree"] == "true" for r in ctx))
    sizes = Counter(int(g["n_loci"]) for g in groups)
    print("groups with 1 locus", sizes[1], "with >= 10 loci", sum(v for k, v in sizes.items() if k >= 10))
    print("context types total", sum(int(g["n_context_types"]) for g in groups))
    for g in groups[:10]:
        kinds = Counter(r["kind"] for r in ctx if r["group"] == g["group"])
        print("  group", g["group"][:8], "loci", g["n_loci"], "leaves", g["n_leaves"], "types", g["n_context_types"],
              "kinds", dict(sorted(kinds.items())))
    by_type = defaultdict(Counter)
    for r in ctx:
        if r["in_tree"] == "true":
            by_type[(r["group"], r["context_type"])][r["kind"]] += 1
    shared = {k: v for k, v in by_type.items() if v["virus"] and v["prokaryote"]}
    print("context types with virus and prokaryote loci", len(shared), "of", len(by_type))
    for k, v in sorted(shared.items(), key=lambda kv: -sum(kv[1].values()))[:10]:
        print("  shared", k[0][:8], "type", k[1], dict(sorted(v.items())))
    pro = [k for k, v in by_type.items() if v["prokaryote"]]
    print("context types with prokaryote loci", len(pro), "prokaryote loci", sum(by_type[k]["prokaryote"] for k in pro))
    top = Counter({k: sum(v.values()) for k, v in by_type.items()}).most_common(5)
    print("largest context types", [(k[0][:8], k[1], n) for k, n in top])
