"""Card W-31: context types and figure colours of `nbr flanks` under the card's alternatives, on one run.

Run in the image (PYTHONPATH = the tree's src), from `run.sh`:

    python measure.py OUT NAME ANCHOR -- <nbr flanks arguments, without --out>

ANCHOR names the anchor's family (TniQ, Cas7f or other). The script runs the unflagged `nbr flanks` once through
`biobox.cli.main` into OUT/NAME/cli_plain (the loci are captured from `flanks.load_loci` on the way), then recomputes
in-process, with the functions `nbr flanks` itself uses (`leaves_of`, `matrices`, `upgma`, `cut_types`, `auto_cut`,
`group_shares`, `_fallback_archs`), the context types under each cut rule and the colour of every flanking gene under
each colouring. It checks that its fixed-cut types equal the CLI's (groups.tsv) and writes OUT/NAME/types.tsv,
colours.tsv and summary.json.

Classes of loci: phage (genome kind virus), cast (a Tn7-like CRISPR-associated locus: the anchor or a flanking gene
carries a TniQ domain and a flanking gene a Cas domain, any Pfam name starting with "Cas"), prok (the others).
"""

from __future__ import annotations

import json
import statistics
import sys
from collections import Counter
from pathlib import Path

from biobox import cli
from biobox.search import flanks

RULES = [("fixed 0.5", 0.5), ("fixed 0.7", 0.7), ("fixed 0.8", 0.8), ("gap", "gap"), ("N=10", 10), ("N=20", 20)]
COLOURINGS = [(w, f) for f in ("none", "pfam") for w in flanks.COLOUR_WEIGHTS]


def domains(arch: str | None) -> list[str]:
    return arch.split("+") if arch else []


def is_cas(g) -> bool:
    return any(d.startswith("Cas") for d in domains(g.arch))


def locus_class(loc, anchor: str) -> str:
    if loc.kind == "virus":
        return "phage"
    genes = loc.up + loc.down
    tniq = anchor == "TniQ" or any("TniQ" in domains(g.arch) for g in genes)
    if tniq and any(is_cas(g) for g in genes):
        return "cast"
    return "prok"


def arg(args: list[str], flag: str, default: str) -> str:
    return args[args.index(flag) + 1] if flag in args else default


def main() -> int:
    out, name, anchor = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
    args = sys.argv[sys.argv.index("--") + 1:]
    d = out / name
    d.mkdir(parents=True, exist_ok=True)
    cap = {}
    real = flanks.load_loci

    def spy(*a, **k):
        cap["loci"], cap["facts"] = real(*a, **k)
        return cap["loci"], cap["facts"]

    flanks.load_loci = spy
    rc = cli.main(args + ["--out", str(d / "cli_plain")])
    assert rc == 0, rc
    loci = cap["loci"]
    window = flanks.parse_window(arg(args, "--window", flanks.DEFAULT_WINDOW))
    group_by = arg(args, "--group-by", flanks.DEFAULT_GROUP_BY)
    level = arg(args, "--level", flanks.DEFAULT_LEVEL)
    cls = {loc.feature_uid: locus_class(loc, anchor) for loc in loci}
    groups: dict[str, list] = {}
    for loc in loci:
        groups.setdefault(flanks.group_of(loc, group_by), []).append(loc)
    order = sorted(groups, key=lambda g: (-len(groups[g]), g))

    # --- context types ---------------------------------------------------------------------------------------------
    types: dict[str, list[list]] = {label: [] for label, _ in RULES}  # rule -> list of types (lists of loci)
    cuts: dict[str, list[float]] = {label: [] for label, _ in RULES}
    for gid in order:
        leaves = flanks.leaves_of(groups[gid], window)
        kept = leaves[:flanks.DEFAULT_MAX_LEAVES]
        assert len(kept) == len(leaves), (gid, len(leaves))  # no group of these runs is sampled
        mj, _, _ = flanks.matrices(kept)
        merges = flanks.upgma(mj, [len(leaf.loci) for leaf in kept])
        del mj
        for label, rule in RULES:
            cut = rule if isinstance(rule, float) else flanks.auto_cut(merges, rule)
            if len(kept) > 1:
                cuts[label].append(cut)
            members: dict[int, list] = {}
            for i, s in enumerate(flanks.cut_types(len(kept), merges, cut)):
                members.setdefault(s, []).extend(kept[i].loci)
            types[label].extend(members.values())
    rows = [r.split("\t") for r in (d / "cli_plain" / "groups.tsv").read_text().splitlines() if not r.startswith("#")]
    at = rows[0].index("n_context_types")
    cli_types = sum(int(r[at]) for r in rows[1:])
    assert cli_types == len(types["fixed 0.5"]), (cli_types, len(types["fixed 0.5"]))

    type_rows = []
    n_by_class = Counter(cls.values())
    for label, _ in RULES:
        sizes = sorted((len(t) for t in types[label]), reverse=True)
        row = {"run": name, "rule": label, "groups": len(order), "loci": len(loci), "types": len(sizes),
               "types_of_1": sum(s == 1 for s in sizes), "types_2_4": sum(2 <= s <= 4 for s in sizes),
               "types_5_9": sum(5 <= s <= 9 for s in sizes), "types_10_49": sum(10 <= s <= 49 for s in sizes),
               "types_50_up": sum(s >= 50 for s in sizes), "largest": ",".join(map(str, sizes[:5])),
               "cut_min": round(min(cuts[label]), 4) if cuts[label] else "",
               "cut_median": round(statistics.median(cuts[label]), 4) if cuts[label] else "",
               "cut_max": round(max(cuts[label]), 4) if cuts[label] else ""}
        for c in ("cast", "phage", "prok"):
            held = [Counter(cls[x.feature_uid] for x in t)[c] for t in types[label]]
            mixed = sum(1 for t in types[label] if Counter(cls[x.feature_uid] for x in t)[c] not in (0, len(t)))
            row[f"{c}_loci"] = n_by_class[c]
            row[f"{c}_types"] = sum(h > 0 for h in held)
            row[f"{c}_largest_type"] = max(held, default=0)
            row[f"{c}_mixed_types"] = mixed
        type_rows.append(row)

    # --- colours ---------------------------------------------------------------------------------------------------
    colour_rows = []
    for weight, fallback in COLOURINGS:
        if fallback == "pfam" and level == "pfam":
            continue
        w = None if weight == "locus" else weight
        fb = None if fallback == "none" else fallback
        palette = flanks.run_colours(groups, w, fb)
        n = Counter()
        families = set()
        for gid in order:
            gl = groups[gid]
            shares = flanks.group_shares(gl, w)
            coloured = {k for k, s in shares.items() if s >= flanks.FIGURE_MIN_SHARE}
            stand_in = flanks._fallback_archs(gl, coloured, w) if fb else set()
            for loc in gl:
                c = cls[loc.feature_uid]
                for g in loc.up + loc.down:
                    if g.key is None:
                        kind = "none"
                    elif g.key in coloured:
                        kind = "cluster"
                        families.add(g.key)
                    elif g.arch in stand_in:
                        kind = "pfam"
                        families.add(flanks.FALLBACK_PREFIX + g.arch)
                    else:
                        kind = "grey"
                    n[("all", kind)] += 1
                    n[(c, kind)] += 1
                    if is_cas(g):
                        n[(c + "_cas", kind)] += 1
        row = {"run": name, "weight": weight, "fallback": fallback, "colours_in_run": len(palette),
               "families_coloured": len(families)}
        for scope in ("all", "cast", "cast_cas", "phage", "phage_cas", "prok", "prok_cas"):
            tot = sum(n[(scope, k)] for k in ("cluster", "pfam", "grey", "none"))
            row[f"{scope}_genes"] = tot
            row[f"{scope}_coloured"] = n[(scope, "cluster")] + n[(scope, "pfam")]
            row[f"{scope}_by_pfam"] = n[(scope, "pfam")]
        colour_rows.append(row)

    for fname, rows in (("types.tsv", type_rows), ("colours.tsv", colour_rows)):
        cols = list(rows[0])
        lines = ["\t".join(cols)] + ["\t".join(str(r[c]) for c in cols) for r in rows]
        (d / fname).write_text("\n".join(lines) + "\n")
    (d / "summary.json").write_text(json.dumps({"run": name, "args": args, "loci": len(loci), "groups": len(order),
                                                "classes": dict(sorted(n_by_class.items())),
                                                "cli_types_fixed": cli_types}, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
