"""Card W-32: what one figure of `nbr flanks` draws, under the defaults and under the old colours.

Run in the runtime image with the card's src on PYTHONPATH (not `python -I`):

    python figure_stats.py OUT GROUP -- <nbr flanks arguments without --out>

It runs `nbr flanks ... --figure` twice in-process, with no colour options and with --colour-weight locus
--colour-fallback none, each into OUT/<default|old>, captures the plan of GROUP that `render_figure` receives (the
figures are still written) and prints one JSON object: per setting, the rows drawn, the flanking genes drawn (grey,
white, coloured by a cluster's colour, coloured by a Pfam colour), the colours used and the legend.
"""

import json
import sys
from collections import Counter
from pathlib import Path

from biobox import cli
from biobox.search import flanks


def stats(plan) -> dict:
    pfam = {c for c, label in plan.legend if label.startswith("Pfam ")}
    n = Counter()
    colours = set()
    for row in plan.rows:
        for g in row.genes:
            if g.anchor:
                continue
            if g.colour == flanks.FIGURE_GREY:
                n["grey"] += 1
            elif g.colour == flanks.FIGURE_NONE:
                n["white"] += 1
            else:
                n["pfam_colour" if g.colour in pfam else "cluster_colour"] += 1
                colours.add(g.colour)
    return {"rows": len(plan.rows), "genes": dict(sorted(n.items())), "colours": len(colours),
            "legend": [label for _, label in plan.legend], "title": plan.title}


def main() -> None:
    out, group = Path(sys.argv[1]), sys.argv[2]
    args = sys.argv[sys.argv.index("--") + 1:]
    got = {}
    real = flanks.render_figure

    def spy(plan, gdir):
        if plan.group == group:
            got["plan"] = plan
        return real(plan, gdir)

    flanks.render_figure = spy
    result = {}
    for name, extra in (("default", []), ("old", ["--colour-weight", "locus", "--colour-fallback", "none"])):
        got.clear()
        rc = cli.main(["nbr", "flanks", *args, "--figure", *extra, "--out", str(out / name)])
        result[name] = {"rc": rc, **(stats(got["plan"]) if "plan" in got else {})}
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
