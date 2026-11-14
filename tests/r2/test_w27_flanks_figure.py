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
    planted shares say; the legend names the coloured families (the release's Pfam names at nr50) with their shares.
    Card W-32: the release has a domain library, so the defaults weight each locus by its anchor's NR90 cluster
    (DEC-41): grouped by that cluster the weighting changes nothing (one unit per group); in one group the shares are
    the oracle's weighted ones (two units, R and R2). The unweighted assertion is kept under --colour-weight locus
    --colour-fallback none, the colours before the card."""
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

    for (colour, label) in r2.legend:  # the grouped legend's shares are the unweighted ones (one unit per group)
        for f, s in shares.items():
            if label.startswith(w25._nr50(f).hex()[:8]):
                assert label.endswith(f"{round(100 * s)}%"), (f, label)

    # before card W-32 (now behind the explicit flags): unweighted shares over the 14 loci
    plans.clear()
    out_all = tmp_path / "all"
    assert run(world, out_all, "--figure", "--group-by", "none", "--colour-weight", "locus",
               "--colour-fallback", "none") == 0
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

    # the defaults since card W-32: each locus weighs 1 / (the loci of its anchor's NR90 cluster), from
    # context_types.tsv's anchor_nr90 column
    plans.clear()
    out_def = tmp_path / "all_default"
    assert run(world, out_def, "--figure", "--group-by", "none") == 0
    (plan,) = plans.values()
    unit = units_from_tables(out_def, loci)
    assert len(set(unit.values())) == 2
    weighted = oracle_weighted_shares(loci, list(loci), unit)
    assert weighted["X"] < flanks.FIGURE_MIN_SHARE
    for row in plan.rows:
        for g in row.genes:
            if g.anchor or g.key is None:
                continue
            f = fam[g.key]
            assert (g.colour == flanks.FIGURE_GREY) == (weighted[f] < flanks.FIGURE_MIN_SHARE), f
    assert any(c == flanks.FIGURE_GREY for c, _ in plan.legend)
    assert any("--colour-weight anchor-nr90" in line for line in plan.title)


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


OLD_COLOURS = ("--colour-weight", "locus", "--colour-fallback", "none")  # the defaults before card W-32


@needs_mpl
@pytest.mark.covers("W-31")
@pytest.mark.parametrize("weight", ["leaf", "genome", "anchor-nr90"])
def test_w31_colour_weight_counts_near_identical_loci_once(mobile, tmp_path, plans, weight):
    """Unweighted, the 30 phage loci colour P1..P10 (30 of 42 loci) while Cas6, in 3 of the 12 species, is grey (3 of
    42 < 10%). --colour-weight leaf, genome or anchor-nr90 counts the phage loci once: 13 units, so P1..P10 are grey
    (1/13) and Cas6 coloured (3/13), as the oracle's weighted shares say for every drawn gene; the legend gives the
    weighted share. Card W-32: the unweighted run is now the explicit --colour-weight locus --colour-fallback none
    (it was the unflagged call), and the weighted run gives --colour-fallback none (the Pfam fallback is a default
    since DEC-41 and would colour the Cas7 clusters, which this test does not ask about)."""
    root, loci = mobile
    fam = mobile_families(loci)
    assert run_colours(mobile, tmp_path / "plain", *OLD_COLOURS) == 0
    (plain,) = plans.values()
    plain_colours = drawn_colours(plain, fam)
    assert plain_colours["P1"] != {flanks.FIGURE_GREY} and plain_colours["Cas6"] == {flanks.FIGURE_GREY}
    plans.clear()
    assert run_colours(mobile, tmp_path / "w", "--colour-weight", weight, "--colour-fallback", "none") == 0
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
    the legend names the family. Card W-32: the run without the fallback is now the explicit --colour-weight locus
    --colour-fallback none (it was the unflagged call), and the run with it gives --colour-weight locus (anchor-nr90 is
    a default since DEC-41 and would colour Cas6, 3 of 13 units)."""
    root, loci = mobile
    fam = mobile_families(loci)
    assert run_colours(mobile, tmp_path / "plain", *OLD_COLOURS) == 0
    (plain,) = plans.values()
    before = drawn_colours(plain, fam)
    for k in range(1, 13):
        assert before[f"Cas7@sp{k}"] == {flanks.FIGURE_GREY}
    plans.clear()
    assert run_colours(mobile, tmp_path / "fb", "--colour-fallback", "pfam", "--colour-weight", "locus") == 0
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
@pytest.mark.covers("W-31", "W-32")
def test_w31_the_default_values_of_the_new_flags_give_the_unflagged_bytes(world, tmp_path, plans):
    """--colour-weight locus and --colour-fallback none are the colours before card W-32: with --figure and --export
    gff every file but the two figures has the bytes the code before the card wrote for the unflagged call (sha256,
    BASE_TABLES, README.txt included), and the figure is drawn as before (no colour line in its title). Before W-32
    this test compared the explicit flags with the unflagged call, which they then equalled; since DEC-41 the unflagged
    call at nr50 on a release with a domain library is weighted and falls back to Pfam: its tables are the same bytes
    (test_w32_tables_are_byte_identical_to_the_code_before_the_card) and its figures are not. The figures' bytes under
    the old flags against the old code are compared on the pod
    (planning/measurements/2026-11-14_flanks_figure_defaults.md), not here: they depend on the matplotlib build."""
    explicit = tmp_path / "explicit"
    assert run(world, explicit, "--figure", "--export", "gff", *OLD_COLOURS) == 0
    n, digest, readme = tables_digest(explicit, readme=False)  # README.txt by its own sha256
    want = BASE_TABLES["world_nr50_gff_figure"]
    assert (n, digest, readme) == (want[0], want[1], want[2])
    for plan in plans.values():
        assert not any("--colour-" in line for line in plan.title)
    assert sorted(p.name for p in explicit.rglob("figure.*")) == ["figure.pdf", "figure.pdf", "figure.svg",
                                                                   "figure.svg"]


@needs_mpl
@pytest.mark.covers("W-31", "W-32")
def test_w31_the_figure_draws_the_cut_auto_cut_chose(mobile, tmp_path, plans):
    """With --auto-cut 2 the figure's dashed line is at the group's chosen cut (groups.tsv) and its title says how
    the cut was chosen. (Before card W-32 this used --auto-cut gap, removed by DEC-41, and found "largest gap" in the
    title.)"""
    root, loci = mobile
    out = tmp_path / "out"
    assert cli.main(["nbr", "flanks", w25.anchors(loci, w25.CUT_LOCI), "--root", str(root), "--release", w25.REL,
                     "--out", str(out), "--figure", "--auto-cut", "2"]) == 0
    (plan,) = plans.values()
    g = {r["group"]: r for r in w25.read_tsv(out / "groups.tsv")}[plan.group]
    assert plan.cut == pytest.approx(float(g["cut"]), abs=1e-6)
    assert any("at most 2 types" in line and "--auto-cut 2" in line for line in plan.title)


# --- card W-32 (DEC-41): the figure's defaults -----------------------------------------------------------------------
#
# At --level nr50 and nr90, when the release's manifest lists a domain library, --colour-weight anchor-nr90 and
# --colour-fallback pfam are the defaults, and a family that would show in two colours takes its Pfam colour
# throughout. The oracle (`oracle_classes`) computes from the planted lists, and the units from context_types.tsv,
# the class every drawn gene must have: its cluster's colour, its Pfam domain's colour, or grey; genes of one class
# share one colour in the run, classes differ in colour. The tables' bytes are compared with the code before the
# card by sha256 (BASE_TABLES, written by that code on the pod:
# planning/measurements/2026-11-14_flanks_figure_defaults/).


def units_from_tables(out: Path, loci: dict) -> dict[str, str]:
    """Locus name -> its anchor's NR90 cluster, as context_types.tsv gives it."""
    names = w25.by_name(loci)
    return {names[r["feature_uid"]]: r["anchor_nr90"] for r in w25.read_tsv(out / "context_types.tsv")}


def oracle_classes(loci: dict, groups: dict[str, list[str]], unit: dict | None, fallback: bool) -> dict:
    """Per group, per planted flanking family (a cluster, `F@x` or `F`): ("cluster", family), ("pfam", domain) or
    "grey". Shares are exact fractions: per unit (each locus its own without `unit`), the share of its loci holding
    the cluster (or its domain, the name before "@"), averaged over the units. A cluster held by at least
    FIGURE_MIN_SHARE is coloured; a grey cluster whose domain is held by at least FIGURE_MIN_SHARE takes the domain's
    colour (fallback); a domain that does so in any group of the run takes over every coloured cluster of that domain,
    in every group (the one-colour rule)."""
    from fractions import Fraction

    dom = lambda f: f.split("@")[0]  # noqa: E731
    least = Fraction(flanks.FIGURE_MIN_SHARE).limit_denominator(1000)

    def shares(names: list[str], of) -> dict[str, Fraction]:
        units: dict[str, list[str]] = {}
        for n in names:
            units.setdefault(unit[n] if unit else n, []).append(n)
        out: dict[str, Fraction] = {}
        for members in units.values():
            for n in members:
                held = set()
                for side in ("up", "down"):
                    for f, _ in loci[n][side]:
                        if f:
                            held.add(of(f))
                for x in held:
                    out[x] = out.get(x, Fraction(0)) + Fraction(1, len(members) * len(units))
        return out

    per = {}
    stand_in = set()
    for g, names in groups.items():
        cs, ds = shares(names, lambda f: f), shares(names, dom)
        coloured = {f for f, s in cs.items() if s >= least}
        per[g] = (cs, ds, coloured)
        if fallback:
            for f in cs:
                if f not in coloured and ds[dom(f)] >= least:
                    stand_in.add(dom(f))
    out = {}
    for g, (cs, ds, coloured) in per.items():
        cls = {}
        for f in cs:
            if f in coloured:
                cls[f] = ("pfam", dom(f)) if dom(f) in stand_in else ("cluster", f)
            elif fallback and ds[dom(f)] >= least:
                cls[f] = ("pfam", dom(f))
            else:
                cls[f] = "grey"
        out[g] = cls
    return out


def check_classes(plans: dict, classes: dict, fam: dict[str, str], uid=None) -> dict:
    """Every drawn flanking gene has its oracle class's colour: grey for "grey"; one colour per class across the run;
    no two classes share a colour; each group's legend has one row per class drawn there (a cluster by its uid's
    first 8 hex digits, a domain as "Pfam <domain>"), none for a cluster drawn in its domain's colour. Returns the
    colour of each class."""
    colour: dict = {}
    for gid, plan in plans.items():
        drawn = set()
        for row in plan.rows:
            for g in row.genes:
                if g.anchor or g.key is None:
                    continue
                c = classes[gid][fam[g.key]]
                if c == "grey":
                    assert g.colour == flanks.FIGURE_GREY, (gid, fam[g.key])
                    continue
                assert g.colour != flanks.FIGURE_GREY, (gid, fam[g.key])
                colour.setdefault(c, set()).add(g.colour)
                drawn.add(c)
        labels = [label for c, label in plan.legend if c not in (flanks.FIGURE_ANCHOR, flanks.FIGURE_GREY,
                                                                    flanks.FIGURE_NONE)]
        want = sorted(("Pfam " + c[1]) if c[0] == "pfam" else (uid or w25._nr50)(c[1]).hex()[:8] for c in drawn)
        assert sorted(label.split(" ")[0] if not label.startswith("Pfam ") else " ".join(label.split(" ")[:2])
                      for label in labels) == want, gid
    for c, cs in colour.items():
        assert len(cs) == 1, c
    assert len({next(iter(cs)) for cs in colour.values()}) == len(colour)
    return {c: next(iter(cs)) for c, cs in colour.items()}


def run_one(world_, out: Path, *extra: str, names=None) -> int:
    root, loci = world_
    return cli.main(["nbr", "flanks", w25.anchors(loci, names), "--root", str(root), "--release", w25.REL, "--out",
                     str(out), *extra])


# The calls whose tables are compared with the code before the card. name: (fixture, anchors, arguments, README.txt
# compared). README.txt gains a paragraph on the colours when a figure is drawn with them (card W-32), so it is
# compared only where the colours do not change: without --figure and at --level pfam.
CALLS = {
    "world_nr50_gff_figure": ("world", None, ["--figure", "--export", "gff"], False),
    "world_nr90_none_combined_figure": ("world", None, ["--figure", "--level", "nr90", "--group-by", "none",
                                                        "--distance", "combined"], False),
    "world_pfam_figure": ("world", None, ["--figure", "--level", "pfam"], True),
    "world_nr50_gff": ("world", None, ["--export", "gff"], True),
    "mobile_colours_none_figure": ("mobile", "colour", ["--figure", "--group-by", "none"], False),
    "one_colour_figure": ("one", None, ["--figure"], False),
}


def build_fixture(kind: str, root: Path) -> dict:
    if kind == "world":
        return w25.build(root)
    if kind == "mobile":
        return w25.build(root, w25.MOBILE, w25.MOBILE_GENOMES)
    return w25.build(root, w25.ONE_COLOUR)


def tables_digest(out: Path, readme: bool) -> tuple[int, str, str | None]:
    """(files, sha256 of the sorted `path<TAB>sha256` lines of every file but figure.svg and figure.pdf, and README.txt
    unless `readme`; sha256 of README.txt)."""
    import hashlib

    lines = []
    for p in sorted(out.rglob("*")):
        rel = p.relative_to(out).as_posix()
        if not p.is_file() or p.name in ("figure.svg", "figure.pdf") or (rel == "README.txt" and not readme):
            continue
        lines.append(f"{rel}\t{hashlib.sha256(p.read_bytes()).hexdigest()}\n")
    text = (out / "README.txt").read_bytes()
    return len(lines), hashlib.sha256("".join(lines).encode()).hexdigest(), hashlib.sha256(text).hexdigest()


def base_digests(scratch: Path) -> dict:
    """Every call of CALLS on fresh fixtures under `scratch`: name -> [files, digest, README sha256]. Run with the code
    before the card (src of 4469fab) to write BASE_TABLES."""
    fixtures: dict[str, tuple] = {}
    out = {}
    for name, (kind, which, args, readme) in CALLS.items():
        if kind not in fixtures:
            root = scratch / kind / "root"
            fixtures[kind] = (root, build_fixture(kind, root))
        names = w25.COLOUR_LOCI if which == "colour" else None
        d = scratch / "out" / name
        assert run_one(fixtures[kind], d, *args, names=names) == 0, name
        out[name] = list(tables_digest(d, readme))
    return out


# written by base_digests with the code of 4469fab (before card W-32) on wilkcomp02, 2026-11-14, in the test image
# (planning/measurements/2026-11-14_flanks_figure_defaults/base_digests.py)
_README_BEFORE = "6ac3c22c979b02fec1ef8425fe9332e3127efc7cabe2742ae36291e9f79cd13e"  # README.txt of every call
BASE_TABLES = {
    "mobile_colours_none_figure": [8, "7d07d85dd3dd8de2250ee3b3132875ba9a6abbb9a4d22295d6dd4c5275d0818b",
                                   _README_BEFORE],
    "one_colour_figure": [13, "3e27c32b63bd4bc242aaae42d6bb0e969fa96af6b53d6af5af8edd523ed20e0f", _README_BEFORE],
    "world_nr50_gff": [29, "355555ac6402ccb3370d75087b0ef0c1f400e615ca885fed131965e6641485db", _README_BEFORE],
    "world_nr50_gff_figure": [28, "f4d5c8d3677ea301228d4a2c31aa66c3f2a263e7c7a37795eac4b1334c498869", _README_BEFORE],
    "world_nr90_none_combined_figure": [8, "8feac075325824bce4e72f6b36b1c58e1f477ca26285820298d068cc095bbb0f",
                                        _README_BEFORE],
    "world_pfam_figure": [14, "3e9497be3a8325c619f57e377c0b84b6e1ddf06969fbceddc9f4fb3c0b8d23d6", _README_BEFORE],
}


@needs_mpl
@pytest.mark.covers("W-32")
@pytest.mark.parametrize("name", list(CALLS))
def test_w32_tables_are_byte_identical_to_the_code_before_the_card(tmp_path, name):
    """The new defaults change the figures only: every other file a call writes (README.txt where the colours do not
    change) has the bytes the code before card W-32 wrote for the same call on the same constructed release."""
    kind, which, args, readme = CALLS[name]
    root = tmp_path / "root"
    world_ = (root, build_fixture(kind, root))
    out = tmp_path / "out"
    assert run_one(world_, out, *args, names=w25.COLOUR_LOCI if which == "colour" else None) == 0
    n, digest, readme_sha = tables_digest(out, readme)
    want = BASE_TABLES[name]
    assert (n, digest) == (want[0], want[1])
    if readme:
        assert readme_sha == want[2]
    else:  # a figure drawn with the new colours: README.txt says how they were chosen
        text = (out / "README.txt").read_text()
        assert readme_sha != want[2] and "--colour-weight" in text and "anchor-nr90" in text


@needs_mpl
@pytest.mark.covers("W-32")
@pytest.mark.parametrize("level", ["nr50", "nr90"])
def test_w32_the_defaults_weight_by_the_anchor_nr90_cluster_and_fall_back_to_pfam(mobile, tmp_path, plans, level):
    """At nr50 and nr90 on a release whose manifest lists a domain library, the unflagged figure of the 42 loci in one
    group counts the 30 near-identical phage loci once (13 units by anchor NR90 cluster, from context_types.tsv):
    P1..P10 are grey (1/13), Cas6 has its cluster's colour (3/13) and the twelve rare Cas7 clusters one colour, that of the
    domain Cas7 (12 of 13 units), as the oracle's classes say; the title names both settings and README.txt says
    which were used. The explicit flags give the same figure plan."""
    root, loci = mobile
    fam = mobile_families(loci) if level == "nr50" else {w25._nr90(f).hex(): f for f in mobile_families(loci).values()}
    out = tmp_path / "out"
    assert run_colours(mobile, out, "--level", level) == 0
    (plan,) = plans.values()
    unit = units_from_tables(out, loci)
    assert len({unit[n] for n in w25.COLOUR_LOCI}) == 13
    classes = oracle_classes(loci, {"all": w25.COLOUR_LOCI}, unit, fallback=True)
    assert classes["all"]["P1"] == "grey" and classes["all"]["Cas6"] == ("cluster", "Cas6")
    assert classes["all"]["Cas7@sp1"] == ("pfam", "Cas7") and classes["all"]["TnsB"] == ("cluster", "TnsB")
    check_classes(plans, classes, fam, uid=w25._nr50 if level == "nr50" else w25._nr90)  # the legend's uid prefix
    assert any("--colour-weight anchor-nr90" in line for line in plan.title)
    assert any("--colour-fallback pfam" in line for line in plan.title)
    readme = (out / "README.txt").read_text()
    assert "--colour-weight anchor-nr90" in readme and "--colour-fallback pfam" in readme
    first = plan
    plans.clear()
    assert run_colours(mobile, tmp_path / "explicit", "--level", level, "--colour-weight", "anchor-nr90",
                       "--colour-fallback", "pfam") == 0
    (again,) = plans.values()
    assert [(r.leaf, [g.colour for g in r.genes]) for r in again.rows] == \
           [(r.leaf, [g.colour for g in r.genes]) for r in first.rows]
    assert again.legend == first.legend


@pytest.fixture(scope="module")
def one_colour(tmp_path_factory):
    root = tmp_path_factory.mktemp("w32one") / "root"
    return root, w25.build(root, w25.ONE_COLOUR)


@needs_mpl
@pytest.mark.covers("W-32")
@pytest.mark.parametrize("extra", [[], ["--colour-fallback", "pfam"]])
def test_w32_a_family_that_would_show_in_two_colours_takes_its_pfam_colour_throughout(one_colour, tmp_path, plans,
                                                                                       extra):
    """In group 2 the rare Cas1 clusters (1 of 12 each) take the colour of the domain Cas1 (12 of 12) while the common
    Cas1 cluster (6 of 12) is coloured as a cluster; the rule draws the common cluster in the domain's colour too, in
    group 2 and in group 1 (where every Cas1 gene is of the common cluster): one colour for Cas1 in the run, a "Pfam
    Cas1" row in both legends and none for the cluster; TnsB keeps its cluster's colour, the same in both groups. As
    the oracle's classes say, with the default flags and with --colour-fallback pfam given. Under --colour-fallback
    none the common cluster has its own colour and the rare ones are grey, as before."""
    root, loci = one_colour
    fam = {w25._nr50(f).hex(): f for n in loci for side in ("up", "down") for f, _ in loci[n][side] if f}
    out = tmp_path / "out"
    assert run_one(one_colour, out, "--figure", *extra) == 0
    assert len(plans) == 2
    groups: dict[str, list[str]] = {}
    for n, f in loci.items():
        groups.setdefault(w25._nr90(f["anchor_family"]).hex(), []).append(n)
    classes = oracle_classes(loci, groups, units_from_tables(out, loci), fallback=True)
    g1, g2 = w25._nr90("O@g1").hex(), w25._nr90("O@g2").hex()
    assert classes[g1]["Cas1@common"] == classes[g2]["Cas1@common"] == classes[g2]["Cas1@sp7"] == ("pfam", "Cas1")
    assert classes[g1]["TnsB"] == ("cluster", "TnsB")
    colour = check_classes(plans, classes, fam)
    cas1 = {g.colour for p in plans.values() for r in p.rows for g in r.genes
            if g.key is not None and not g.anchor and fam[g.key].startswith("Cas1@")}
    assert cas1 == {colour[("pfam", "Cas1")]} and colour[("cluster", "TnsB")] not in cas1
    common = w25._nr50("Cas1@common").hex()[:8]
    for p in plans.values():
        labels = [label for _, label in p.legend]
        assert any(label.startswith("Pfam Cas1 ") for label in labels)
        assert not any(label.startswith(common) for label in labels)

    plans.clear()
    assert run_one(one_colour, tmp_path / "none", "--figure", "--colour-fallback", "none") == 0
    old = oracle_classes(loci, groups, units_from_tables(out, loci), fallback=False)
    assert old[g2]["Cas1@common"] == ("cluster", "Cas1@common") and old[g2]["Cas1@sp7"] == "grey"
    check_classes(plans, old, fam)


@needs_mpl
@pytest.mark.covers("W-32")
@pytest.mark.parametrize("how", ["no manifest", "no domain library"])
def test_w32_without_a_domain_library_in_the_manifest_the_defaults_are_the_old_ones_and_stderr_says_so(
        tmp_path_factory, tmp_path, plans, capsys, how):
    """A release read from its release_def only (no manifest), or whose manifest lists no l1_domains file: the figure
    is drawn with --colour-weight locus --colour-fallback none (P1..P10 coloured by the 30 phage loci, Cas6 and every
    Cas7 cluster grey, as the unweighted oracle says) and one stderr line says why; without --figure no line."""
    root = tmp_path_factory.mktemp("w32nolib") / "root"
    loci = w25.build(root, w25.MOBILE, w25.MOBILE_GENOMES, manifest=how != "no manifest",
                     domains=how != "no domain library")
    mobile_ = (root, loci)
    capsys.readouterr()
    assert run_colours(mobile_, tmp_path / "out") == 0
    err = capsys.readouterr().err
    lines = [x for x in err.splitlines() if "domain library" in x]
    assert len(lines) == 1, err
    assert "--colour-weight locus --colour-fallback none" in lines[0]
    assert ("manifest" in lines[0])
    (plan,) = plans.values()
    classes = oracle_classes(loci, {"all": w25.COLOUR_LOCI}, None, fallback=False)
    assert classes["all"]["P1"] == ("cluster", "P1") and classes["all"]["Cas6"] == "grey"
    check_classes(plans, classes, mobile_families(loci))
    assert not any("--colour-" in line for line in plan.title)
    assert "--colour-weight" not in (tmp_path / "out" / "README.txt").read_text()
    assert run_one(mobile_, tmp_path / "plain", "--group-by", "none", names=w25.COLOUR_LOCI) == 0
    assert "domain library" not in capsys.readouterr().err


@needs_mpl
@pytest.mark.covers("W-32")
def test_w32_at_level_pfam_the_defaults_do_not_change(mobile, tmp_path, plans, capsys):
    """At --level pfam the families are Pfam architectures already: the defaults stay locus and none, the title
    carries no colour line and stderr no note; an explicit --colour-weight anchor-nr90 is still honoured there."""
    root, loci = mobile
    capsys.readouterr()
    assert run_colours(mobile, tmp_path / "out", "--level", "pfam") == 0
    assert "domain library" not in capsys.readouterr().err
    (plan,) = plans.values()
    assert not any("--colour-" in line for line in plan.title)
    assert "--colour-weight" not in (tmp_path / "out" / "README.txt").read_text()
    names = {f.split("@")[0] for n in w25.COLOUR_LOCI for side in ("up", "down") for f, _ in loci[n][side] if f}
    shares = oracle_weighted_shares(loci, w25.COLOUR_LOCI, None, of=lambda f: f.split("@")[0])
    got = drawn_colours(plan, {d: d for d in names})
    assert got["P1"] != {flanks.FIGURE_GREY} and shares["P1"] >= flanks.FIGURE_MIN_SHARE
    assert got["Cas6"] == {flanks.FIGURE_GREY} and shares["Cas6"] < flanks.FIGURE_MIN_SHARE
    plans.clear()
    assert run_colours(mobile, tmp_path / "w", "--level", "pfam", "--colour-weight", "anchor-nr90") == 0
    (weighted,) = plans.values()
    got = drawn_colours(weighted, {d: d for d in names})
    assert got["P1"] == {flanks.FIGURE_GREY} and got["Cas6"] != {flanks.FIGURE_GREY}
