"""Ring R2: the synteny figure of `bio-box nbr flanks --figure` (card W-27).

The constructed release is card W-25's (tests/r2/test_w25_flanks.py, loaded here as a separate module copy) with two
more loci, r3_a and r3_b, whose anchors are in NR90 cluster R2 (the group of U1 and U2) and whose flanks carry T1's
upstream families A1..A5 and T2's downstream families D1..D5 (r3_b with D5 unannotated and its anchor on the minus
strand), so that the same families appear in two groups. Genes are 300 bp with 60 bp gaps, so every drawn position is
known exactly.

The oracle is the planted lists: the families, strands and order of each locus's flanks, the share of a group's loci
that hold each family (counted by a loop over the planted lists), the contig start that cuts `t1_cens`, and the
coordinates the builder laid out. The figure's content is checked on the plan `render_figure` receives (captured by
wrapping it), its files on disk (SVG parsed as XML, PDF header and trailer, no date), and their bytes across processes.
"""

from __future__ import annotations

import importlib.util
import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from biobox import cli
from biobox.search import flanks

pytestmark = [pytest.mark.r2, pytest.mark.covers("W-27")]

HAS_MPL = importlib.util.find_spec("matplotlib") is not None
needs_mpl = pytest.mark.skipif(not HAS_MPL, reason="matplotlib is not importable here (the image built from "
                                                   "`main` has it; where it is absent, for example on a laptop, the test skips); "
                                                   "GitHub Actions installs it from uv.lock and runs this test")

_spec = importlib.util.spec_from_file_location("w25_world_for_w27", Path(__file__).with_name("test_w25_flanks.py"))
w25 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(w25)
w25.PLANTED = {**w25.PLANTED,
               "r3_a": ("P1", "prokaryote", "R2", 1, w25.A + w25.EXTRA, w25.D + w25.EXTRA, False),
               # a leaf of its own (D5 replaced by an unannotated gene) with the anchor on the minus strand
               "r3_b": ("S1", "virus", "R2", -1, w25.A + w25.EXTRA, w25.D[:4] + [(None, 1)] + w25.EXTRA, False)}
FAMILIES = sorted({f for v in w25.PLANTED.values() for side in (v[4], v[5]) for f, _ in side if f})
HALF = w25.GENE / 2


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    root = tmp_path_factory.mktemp("w27") / "root"
    return root, w25.build(root)


def run(world, out: Path, *extra: str) -> int:
    root, loci = world
    return cli.main(["nbr", "flanks", w25.anchors(loci), "--root", str(root), "--release", w25.REL, "--out", str(out),
                     *extra])


@pytest.fixture
def plans(monkeypatch):
    """The plans `render_figure` draws, by group, captured while the real function still writes the files."""
    got: dict[str, flanks.FigurePlan] = {}
    real = flanks.render_figure

    def spy(plan, gdir):
        got[plan.group] = plan
        return real(plan, gdir)

    monkeypatch.setattr(flanks, "render_figure", spy)
    return got


def family_of_key(level: str) -> dict[str, str]:
    """The key a family's genes carry at `level` -> the family's planted name."""
    if level == "pfam":
        return {f: f for f in FAMILIES}
    return {w25._nr50(f).hex(): f for f in FAMILIES}


def oracle_shares(loci: dict, names: list[str]) -> dict[str, float]:
    """Per planted family: the share of the loci `names` whose window (W genes on each side) holds it."""
    out = {}
    for fam in FAMILIES:
        holders = 0
        for n in names:
            found = False
            for side in ("up", "down"):
                for f, _ in loci[n][side]:
                    if f == fam:
                        found = True
            holders += found
        if holders:
            out[fam] = holders / len(names)
    return out


def groups_of(loci: dict) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for n, f in loci.items():
        out.setdefault(w25._nr90(f["anchor_family"]).hex(), []).append(n)
    return out


# --- tests ----------------------------------------------------------------------------------------------------------


@needs_mpl
def test_w27_figure_writes_an_svg_and_a_pdf_per_group_without_a_date(world, tmp_path):
    """--figure: every group directory holds figure.svg (an SVG document) and figure.pdf (a PDF), with no creation
    date in either; without --figure no figure is written."""
    out = tmp_path / "out"
    assert run(world, out, "--figure") == 0
    gdirs = sorted(p for p in (out / "groups").iterdir() if p.is_dir())
    assert len(gdirs) == 2
    for d in gdirs:
        svg, pdf = (d / "figure.svg").read_bytes(), (d / "figure.pdf").read_bytes()
        assert ET.fromstring(svg).tag == "{http://www.w3.org/2000/svg}svg"
        assert b"<dc:date>" not in svg
        assert pdf.startswith(b"%PDF-") and pdf.rstrip().endswith(b"%%EOF")
        assert b"/CreationDate" not in pdf and b"/ModDate" not in pdf
    plain = tmp_path / "plain"
    assert run(world, plain) == 0
    assert not list(plain.rglob("figure.*"))
    assert "figure.svg" in (out / "README.txt").read_text()


@needs_mpl
def test_w27_rows_are_the_leaves_in_tree_order_beside_the_tree(world, tmp_path, plans):
    """One row per leaf in the tree, in the order of leaves.tsv's tree_order; the tree drawn beside them has one
    junction per merge (rows - 1) and its root at twice the Newick root height (heights are merge distance / 2)."""
    root, loci = world
    out = tmp_path / "out"
    assert run(world, out, "--figure") == 0
    for gid, plan in plans.items():
        leaves = [r for r in w25.read_tsv(out / "groups" / gid / "leaves.tsv") if r["in_tree"] == "true"]
        want = [r["leaf"] for r in sorted(leaves, key=lambda r: int(r["tree_order"]))]
        assert [row.leaf for row in plan.rows] == want, gid
        assert [row.n_loci for row in plan.rows] == [int(r["n_loci"]) for r in sorted(
            leaves, key=lambda r: int(r["tree_order"]))]
        junctions = [s for s in plan.segments if s[0] == s[2] and s[1] != s[3]]  # one per merge, joining two rows
        assert len(junctions) == len(plan.rows) - 1, gid
        if len(plan.rows) > 1:
            text = (out / "groups" / gid / "tree.nwk").read_text().strip()
            first = re.search(r"\(([0-9a-f]{32}):([0-9.]+)", text)
            depth, height = 0, float(first.group(2))
            for ch, m in ((text[i], i) for i in range(first.end(), len(text))):
                if ch == "(":
                    depth += 1
                elif ch == ")":
                    if depth == 0:
                        bl = re.match(r"\):([0-9.]+)", text[m:])
                        height += float(bl.group(1)) if bl else 0.0
                    else:
                        depth -= 1
            assert max(s[0] for s in junctions) == pytest.approx(2 * height, abs=1e-5), gid


@needs_mpl
@pytest.mark.parametrize("level", ["nr50", "pfam"])
def test_w27_genes_are_arrows_along_the_anchor_strand_with_the_anchor_centred(world, tmp_path, plans, level):
    """Each row draws its leaf's first locus: the anchor centred on 0 and marked, its W upstream genes to the left and
    downstream genes to the right at the planted positions (300 bp genes, 60 bp gaps), whatever the anchor's strand;
    each arrow points right when the gene is on the anchor's strand; each gene's family is the planted one (none for
    an unannotated gene). A contig start cutting t1_cens upstream is a break mark at its distance; no other row has
    one."""
    root, loci = world
    out = tmp_path / "out"
    assert run(world, out, "--figure", "--level", level) == 0
    names = w25.by_name(loci)
    fam = family_of_key(level)
    strands = set()
    for plan in plans.values():
        for row in plan.rows:
            n = names[row.drawn]
            assert row.drawn == row.leaf  # the leaf's name is its least feature_uid, the locus drawn
            strands.add(w25.PLANTED[n][3])
            anchor = [g for g in row.genes if g.anchor]
            assert len(anchor) == 1 and (anchor[0].x0, anchor[0].x1) == (-HALF, HALF) and anchor[0].forward
            assert anchor[0].colour == flanks.FIGURE_ANCHOR and anchor[0].feature_uid == row.drawn
            for side, sign in (("up", -1), ("down", 1)):
                drawn = sorted((g for g in row.genes if not g.anchor and (g.x0 + g.x1) * sign > 0),
                               key=lambda g: abs(g.x0 + g.x1))
                planted = loci[n][side]
                assert len(drawn) == len(planted), (n, side)
                for k, (g, (f, strand)) in enumerate(zip(drawn, planted)):
                    inner = HALF + w25.GAP + k * (w25.GENE + w25.GAP)
                    assert sorted((g.x0, g.x1)) == sorted((sign * inner, sign * (inner + w25.GENE))), (n, side, k)
                    assert g.x0 < g.x1
                    assert g.forward == (strand == 1), (n, side, k)
                    assert (fam.get(g.key) if g.key is not None else None) == f, (n, side, k)
                    if f is None:
                        assert g.colour == flanks.FIGURE_NONE
            if n == "t1_cens":
                assert row.breaks == [-(HALF + 2 * (w25.GENE + w25.GAP))]
            else:
                assert row.breaks == []
    assert any(names[r.drawn] == "t1_cens" for p in plans.values() for r in p.rows)
    assert strands == {1, -1}  # r3_b, alone in its leaf, has its anchor on the minus strand


@needs_mpl
def test_w27_a_family_has_one_colour_in_every_group_and_rare_families_are_grey(world, tmp_path, plans):
    """Grouped by the anchor's NR90 cluster, A1..A5 and D1..D5 lie in two groups (R and R2) and carry the same colour
    in both, a different colour per family. In one group of every locus (--group-by none), a family held by fewer
    than FIGURE_MIN_SHARE of the loci (X, one locus of 14) is grey and every other family is coloured, as the
    planted shares say; the legend names the coloured families (the release's Pfam names at nr50) with their shares."""
    root, loci = world
    out = tmp_path / "out"
    assert run(world, out, "--figure") == 0
    fam = family_of_key("nr50")
    colour: dict[str, set[str]] = {}
    for plan in plans.values():
        for row in plan.rows:
            for g in row.genes:
                if not g.anchor and g.key is not None:
                    colour.setdefault(fam[g.key], set()).add(g.colour)
    for f in [f"A{k}" for k in range(1, 6)] + [f"D{k}" for k in range(1, 6)]:
        assert len(colour[f]) == 1, f
    assert len({c for f in colour for c in colour[f]}) == len(colour)  # one colour per family, none shared
    groups = groups_of(loci)
    r2 = plans[w25._nr90("R2").hex()]
    shares = oracle_shares(loci, groups[w25._nr90("R2").hex()])
    labels = [label for _, label in r2.legend]
    for f, s in shares.items():
        assert any(label.startswith(w25._nr50(f).hex()[:8]) and f in label and f"{round(100 * s)}%" in label
                   for label in labels), f

    plans.clear()
    out_all = tmp_path / "all"
    assert run(world, out_all, "--figure", "--group-by", "none") == 0
    (plan,) = plans.values()
    shares = oracle_shares(loci, list(loci))
    assert shares["X"] < flanks.FIGURE_MIN_SHARE
    for row in plan.rows:
        for g in row.genes:
            if g.anchor or g.key is None:
                continue
            f = fam[g.key]
            assert (g.colour == flanks.FIGURE_GREY) == (shares[f] < flanks.FIGURE_MIN_SHARE), f
    assert any(c == flanks.FIGURE_GREY for c, _ in plan.legend)


@needs_mpl
def test_w27_above_the_row_bound_the_figure_shows_the_stated_leaves_and_says_so(world, tmp_path, plans):
    """--figure-rows 2: the R group's six leaves form two context types; the figure draws the leaf with the most loci
    of type 1 (T1's four identical loci) and of type 2 (T2's two), in tree order, says '2 of 6 leaves' on the figure
    (in the SVG's text) and the README states the rule."""
    root, loci = world
    out = tmp_path / "out"
    assert run(world, out, "--figure", "--figure-rows", "2") == 0
    gid = w25._nr90("R").hex()
    plan = plans[gid]
    names = w25.by_name(loci)
    first = lambda group: min(group, key=lambda n: loci[n]["feature_uid"])  # noqa: E731
    assert {names[r.drawn] for r in plan.rows} == {first(("t1_v1", "t1_v2", "t1_p1", "t1_p2")),
                                                   first(("t2_p1", "t2_p2"))}
    assert sorted(r.n_loci for r in plan.rows) == [2, 4]
    assert sorted(r.context_type for r in plan.rows) == [1, 2]
    leaves = {r["leaf"]: int(r["tree_order"]) for r in w25.read_tsv(out / "groups" / gid / "leaves.tsv")}
    assert [r.leaf for r in plan.rows] == sorted((r.leaf for r in plan.rows), key=leaves.get)
    assert len([s for s in plan.segments if s[0] == s[2] and s[1] != s[3]]) == 1
    svg = (out / "groups" / gid / "figure.svg").read_text()
    assert "2 of 6 leaves" in svg
    assert "2 of 6 leaves" in " ".join(plan.title)
    readme = (out / "README.txt").read_text()
    assert "--figure-rows" in readme and "most loci of each context type" in readme


@needs_mpl
def test_w27_columns_name_the_context_type_kind_stage_taxonomy_and_collapsed_loci(world, tmp_path, plans):
    """The columns at the right of each row: the leaf's loci (collapsed count), its context type, its loci's genome
    kinds and stages (each value with its count when they differ), taxonomy (none in this release) and the anchor's
    NR90 cluster; equal to the per-locus rows of context_types.tsv."""
    root, loci = world
    out = tmp_path / "out"
    assert run(world, out, "--figure") == 0
    rows = w25.read_tsv(out / "context_types.tsv")
    for gid, plan in plans.items():
        assert plan.columns[:2] == ["accession", "loci"]
        for row in plan.rows:
            mine = [r for r in rows if r["leaf"] == row.leaf]
            col = dict(zip(plan.columns, row.columns))
            assert col["loci"] == str(len(mine)) == str(row.n_loci)
            assert col["type"] == mine[0]["context_type"] == str(row.context_type)
            for name in ("kind", "stage"):
                values = sorted({r[name] for r in mine})
                if len(values) == 1:
                    assert col[name] == values[0]
                else:
                    for v in values:
                        assert f"{v}:{sum(r[name] == v for r in mine)}" in col[name]
            assert col["taxon"] == "-"
            assert col["anchor NR90"] == mine[0]["anchor_nr90"][:8]
            assert col["accession"] == min(mine, key=lambda r: r["feature_uid"])["accession"]


@needs_mpl
def test_w27_figure_bytes_are_identical_across_hash_seeds_and_threads(world, tmp_path):
    """Two processes, one with 1 DuckDB thread and PYTHONHASHSEED=0, the other with 3 threads and PYTHONHASHSEED=4242,
    each with its own empty matplotlib config directory: every figure.svg and figure.pdf is byte-identical."""
    root, loci = world
    outs = []
    for k, (threads, seed) in enumerate([("1", "0"), ("3", "4242")]):
        out = tmp_path / f"out{k}"
        mpl_dir = tmp_path / f"mpl{k}"
        mpl_dir.mkdir()
        env = {**os.environ, "PYTHONHASHSEED": seed, "MPLCONFIGDIR": str(mpl_dir)}
        got = subprocess.run([sys.executable, "-m", "biobox.cli", "nbr", "flanks", w25.anchors(loci), "--root",
                              str(root), "--release", w25.REL, "--out", str(out), "--threads", threads, "--figure"],
                             env=env, capture_output=True, text=True)
        assert got.returncode == 0, got.stderr
        outs.append(out)
    figs = sorted(p.relative_to(outs[0]) for p in outs[0].rglob("figure.*"))
    assert len(figs) == 4  # two groups
    for p in figs:
        assert (outs[0] / p).read_bytes() == (outs[1] / p).read_bytes(), p


def test_w27_figure_is_refused_without_matplotlib_and_nothing_is_written(world, tmp_path, capsys, monkeypatch):
    """Where matplotlib cannot be imported, --figure is refused (exit 2) before any output is written."""
    for name in [m for m in sys.modules if m == "matplotlib" or m.startswith("matplotlib.")]:
        monkeypatch.delitem(sys.modules, name)
    monkeypatch.setitem(sys.modules, "matplotlib", None)
    out = tmp_path / "out"
    assert run(world, out, "--figure") == 2
    assert "matplotlib" in capsys.readouterr().err
    assert not out.exists()


@pytest.mark.parametrize("extra, needle", [
    (["--figure", "--figure-rows", "0"], "--figure-rows"),
    (["--figure-rows", "5"], "--figure"),
])
def test_w27_figure_rows_refusals(world, tmp_path, capsys, extra, needle):
    """--figure-rows below 1, or without --figure, is refused (exit 2) and nothing is written."""
    out = tmp_path / "out"
    assert run(world, out, *extra) == 2
    assert needle in capsys.readouterr().err
    assert not out.exists()


# --- card W-31: colours that count near-identical loci once, and the Pfam fallback ----------------------------------
#
# The release of w25.MOBILE: 30 near-identical phage loci (one anchor NR90 cluster, one genome, one leaf) and 12
# transposon-like loci of 12 species (each its own anchor NR90 cluster, genome and leaf), in one group (--group-by
# none). By construction every weighting (leaf, genome, anchor NR90 cluster) has the same 13 units: the phage loci one,
# each tn locus one. The oracle counts the planted lists by unit.

UNIT = {n: ("phage" if n.startswith("ph_") else n) for n in w25.COLOUR_LOCI}


@pytest.fixture(scope="module")
def mobile(tmp_path_factory):
    root = tmp_path_factory.mktemp("w31fig") / "root"
    return root, w25.build(root, w25.MOBILE, w25.MOBILE_GENOMES)


def run_colours(mobile, out: Path, *extra: str) -> int:
    root, loci = mobile
    return cli.main(["nbr", "flanks", w25.anchors(loci, w25.COLOUR_LOCI), "--root", str(root), "--release", w25.REL,
                     "--out", str(out), "--group-by", "none", "--figure", *extra])


def oracle_weighted_shares(loci: dict, names: list[str], unit: dict | None, of=lambda f: f) -> dict[str, float]:
    """Per family (or `of(family)`, e.g. its domain): the mean over units of the share of the unit's loci whose window
    holds it; each locus its own unit when `unit` is None."""
    units: dict[str, list[str]] = {}
    for n in names:
        units.setdefault(unit[n] if unit else n, []).append(n)
    out: dict[str, float] = {}
    for members in units.values():
        for n in members:
            held = {of(f) for side in ("up", "down") for f, _ in loci[n][side] if f}
            for f in held:
                out[f] = out.get(f, 0.0) + 1 / len(members) / len(units)
    return out


def drawn_colours(plan, fam: dict[str, str]) -> dict[str, set[str]]:
    """Planted family -> the colours its drawn genes carry."""
    got: dict[str, set[str]] = {}
    for row in plan.rows:
        for g in row.genes:
            if not g.anchor and g.key is not None:
                got.setdefault(fam[g.key], set()).add(g.colour)
    return got


def mobile_families(loci: dict) -> dict[str, str]:
    names = {f for n in w25.COLOUR_LOCI for side in ("up", "down") for f, _ in loci[n][side] if f}
    return {w25._nr50(f).hex(): f for f in names}


@needs_mpl
@pytest.mark.covers("W-31")
@pytest.mark.parametrize("weight", ["leaf", "genome", "anchor-nr90"])
def test_w31_colour_weight_counts_near_identical_loci_once(mobile, tmp_path, plans, weight):
    """Unweighted, the 30 phage loci colour P1..P10 (30 of 42 loci) while Cas6, in 3 of the 12 species, is grey (3 of
    42 < 10%). --colour-weight leaf, genome or anchor-nr90 counts the phage loci once: 13 units, so P1..P10 are grey
    (1/13) and Cas6 coloured (3/13), as the oracle's weighted shares say for every drawn gene; the legend gives the
    weighted share."""
    root, loci = mobile
    fam = mobile_families(loci)
    assert run_colours(mobile, tmp_path / "plain") == 0
    (plain,) = plans.values()
    plain_colours = drawn_colours(plain, fam)
    assert plain_colours["P1"] != {flanks.FIGURE_GREY} and plain_colours["Cas6"] == {flanks.FIGURE_GREY}
    plans.clear()
    assert run_colours(mobile, tmp_path / "w", "--colour-weight", weight) == 0
    (plan,) = plans.values()
    shares = oracle_weighted_shares(loci, w25.COLOUR_LOCI, UNIT)
    assert shares["P1"] == pytest.approx(1 / 13) and shares["Cas6"] == pytest.approx(3 / 13)
    got = drawn_colours(plan, fam)
    for f, colours in got.items():
        assert (colours == {flanks.FIGURE_GREY}) == (shares[f] < flanks.FIGURE_MIN_SHARE), f
    assert len(got["Cas6"]) == 1 and len(got["TnsB"]) == 1 and got["TnsB"] != {flanks.FIGURE_GREY}
    labels = [label for _, label in plan.legend]
    assert any(f"{round(100 * 3 / 13)}%" in label and "Cas6" in label for label in labels)
    assert any(weight in line for line in plan.title)


@needs_mpl
@pytest.mark.covers("W-31")
def test_w31_colour_fallback_colours_a_rare_cluster_by_its_common_pfam_family(mobile, tmp_path, plans):
    """Each species' Cas7 is its own NR50 cluster (1 of 42 loci, grey) but the 12 share the Pfam-like domain Cas7 (12 of
    42 loci). --colour-fallback pfam colours the 12 Cas7 genes with one colour of their family, a colour no cluster
    has; Cas6 stays grey (its family too is in 3 of 42); clusters coloured without the fallback keep their colours;
    the legend names the family."""
    root, loci = mobile
    fam = mobile_families(loci)
    assert run_colours(mobile, tmp_path / "plain") == 0
    (plain,) = plans.values()
    before = drawn_colours(plain, fam)
    for k in range(1, 13):
        assert before[f"Cas7@sp{k}"] == {flanks.FIGURE_GREY}
    plans.clear()
    assert run_colours(mobile, tmp_path / "fb", "--colour-fallback", "pfam") == 0
    (plan,) = plans.values()
    after = drawn_colours(plan, fam)
    domain = oracle_weighted_shares(loci, w25.COLOUR_LOCI, None, of=lambda f: f.split("@")[0])
    assert domain["Cas7"] == pytest.approx(12 / 42) and domain["Cas6"] == pytest.approx(3 / 42)
    cas7 = set().union(*(after[f"Cas7@sp{k}"] for k in range(1, 13)))
    assert len(cas7) == 1 and cas7 != {flanks.FIGURE_GREY}
    others = {c for f, cs in after.items() if not f.startswith("Cas7@") for c in cs}
    assert not cas7 & others
    assert after["Cas6"] == {flanks.FIGURE_GREY}
    for f in before:
        if not f.startswith("Cas7@"):
            assert after[f] == before[f], f
    assert any("Cas7" in label and colour in cas7 for colour, label in plan.legend)


@needs_mpl
@pytest.mark.covers("W-31")
def test_w31_the_default_values_of_the_new_flags_give_the_unflagged_bytes(world, tmp_path):
    """--colour-weight locus and --colour-fallback none are today's colours: with --figure and --export gff every file
    is byte-identical to the unflagged call's, built here from it."""
    plain, explicit = tmp_path / "plain", tmp_path / "explicit"
    assert run(world, plain, "--figure", "--export", "gff") == 0
    assert run(world, explicit, "--figure", "--export", "gff", "--colour-weight", "locus",
               "--colour-fallback", "none") == 0
    a = sorted(p.relative_to(plain) for p in plain.rglob("*") if p.is_file())
    b = sorted(p.relative_to(explicit) for p in explicit.rglob("*") if p.is_file())
    assert a == b and len(a) > 10
    for p in a:
        assert (plain / p).read_bytes() == (explicit / p).read_bytes(), p


@needs_mpl
@pytest.mark.covers("W-31")
def test_w31_the_figure_draws_the_cut_auto_cut_chose(mobile, tmp_path, plans):
    """With --auto-cut gap the figure's dashed line is at the group's chosen cut (groups.tsv) and its title says how
    the cut was chosen."""
    root, loci = mobile
    out = tmp_path / "out"
    assert cli.main(["nbr", "flanks", w25.anchors(loci, w25.CUT_LOCI), "--root", str(root), "--release", w25.REL,
                     "--out", str(out), "--figure", "--auto-cut", "gap"]) == 0
    (plan,) = plans.values()
    g = {r["group"]: r for r in w25.read_tsv(out / "groups.tsv")}[plan.group]
    assert plan.cut == pytest.approx(float(g["cut"]), abs=1e-6)
    assert any("largest gap" in line for line in plan.title)
