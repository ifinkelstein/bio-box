"""Card W-31: `nbr flanks` on the constructed releases of tests/r2/test_w25_flanks.py, without the card's flags.

    python fixture_identity.py TREE OUT

Builds the W-25 release and the W-31 one (MOBILE) from TREE's test module under OUT, then runs `nbr flanks` (the
`biobox` on PYTHONPATH) with several unflagged settings into OUT/runs/<name>. Run once with the code before the card
and once with the card's code, each with its own tree first on PYTHONPATH; `diff -r` of the two OUT/runs must be
empty.
"""

import importlib.util
import sys
from pathlib import Path

from biobox import cli

tree, out = Path(sys.argv[1]), Path(sys.argv[2])
spec = importlib.util.spec_from_file_location("w25_world", tree / "tests" / "r2" / "test_w25_flanks.py")
w25 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w25)
root, mroot = out / "root", out / "mobile_root"
loci = w25.build(root)
mloci = w25.build(mroot, w25.MOBILE, w25.MOBILE_GENOMES)
RUNS = {
    "default": (root, w25.anchors(loci), []),
    "figure_gff": (root, w25.anchors(loci), ["--figure", "--export", "gff"]),
    "pfam_figure": (root, w25.anchors(loci), ["--level", "pfam", "--figure"]),
    "none_combined_figure": (root, w25.anchors(loci), ["--group-by", "none", "--distance", "combined", "--figure"]),
    "nr90_a_cut03": (root, w25.anchors(loci), ["--level", "nr90", "--distance", "a", "--cut", "0.3"]),
    "mobile_colours_figure": (mroot, w25.anchors(mloci, w25.COLOUR_LOCI), ["--group-by", "none", "--figure"]),
    "mobile_cut_figure": (mroot, w25.anchors(mloci, w25.CUT_LOCI), ["--figure", "--export", "gff"]),
}
for name, (r, anchors, extra) in RUNS.items():
    rc = cli.main(["nbr", "flanks", anchors, "--root", str(r), "--release", w25.REL, "--out", str(out / "runs" / name),
                   *extra])
    assert rc == 0, (name, rc)
print(f"fixture_identity: {len(RUNS)} runs under {out / 'runs'}")
