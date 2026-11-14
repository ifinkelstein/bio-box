"""Ring R2, card E-15: `bio-box nbr` and `nbr summarize` show the arrays of both CRISPR producers, the default
(`minced-...-nr3-...`) and the filtered one of two-repeat arrays beside Cas genes (`minced-...-nr2-...-cas1kb-...`),
each row with its producer, and the summary counts them per producer.

The world and the oracle are card E-14's (`tests/r2/test_e14_arrays_nbr.py`): the micro-S0 search world with arrays
planted as the producers write them, and the window and array rows computed by loops over the tapes and the planted
lists. Here the planted arrays are split between the two producers' partitions (the filtered producer's are the
two-repeat ones), and nothing of `biobox.search.arrays` is told about the filtered producer by name.
"""

from __future__ import annotations

import random
import shutil
import sys
from pathlib import Path

import pytest

from biobox.search import arrays, nbr

sys.path.insert(0, str(Path(__file__).parent))
import test_e14_arrays_nbr as e14  # noqa: E402
import test_micro_s0_search as ms  # noqa: E402
from test_micro_s0_search import world  # noqa: E402,F401  (the module-scoped search world)

pytestmark = [pytest.mark.r2, pytest.mark.covers("E-15", "E-14", "S-09")]

P = e14.P
NEAR = "minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2"


def plan(tps: dict[bytes, dict]) -> list[dict]:
    """Arrays of the default producer (3 to 6 repeats) after every seventh gene, and two-repeat arrays of the
    filtered producer 150 bases after every eleventh gene, one across the circular origin."""
    rng = random.Random(1514)
    out = []
    for cu, t in sorted(tps.items()):
        n, circ = t["contig_len"], bool(t["is_circular"])
        sl = list(zip(t["start_bps"], t["end_bps"]))
        for k in range(2, len(sl) - 2, 7):
            out.append(e14._array(rng, cu, n, circ, sl[k][1] + rng.randint(5, 60), rng.randint(3, 6), P))
        for k in range(5, len(sl) - 2, 11):
            out.append(e14._array(rng, cu, n, circ, sl[k][1] + 150, 2, NEAR))
        if circ:
            out.append(e14._array(rng, cu, n, circ, n - 60, 2, NEAR))  # across the origin: end = n + 35
    return out


@pytest.fixture(scope="module")
def both(world, tmp_path_factory):
    root = tmp_path_factory.mktemp("e15-nbr") / "a"
    shutil.copytree(world[0], root)
    shutil.copytree(ms._packs(world), ms._packs_of(root))
    tps = ms._tape_rows(root, e14.TV)
    planted = plan(tps)
    e14.plant(root, planted)
    return {"root": root, "tapes": tps, "planted": planted}


def test_locus_tables_show_the_arrays_of_both_producers(world, both, capsys):
    uids = e14.anchors(both["tapes"], both["root"])
    policy = nbr.Policy.parse(window="3g")
    base = e14._nbr_tsv(world[0], uids, "3g", capsys).rstrip("\n").split("\n\n")
    new = e14._nbr_tsv(both["root"], uids, "3g", capsys).rstrip("\n").split("\n\n")
    assert len(base) == len(new) > 20
    seen = []
    for b, g in zip(base, new):
        want, rows = e14.oracle_table(both, b.split("\n"), policy)
        assert g.split("\n") == want
        seen += rows
    pcol, ccol = arrays.COLUMNS.index("producer"), arrays.COLUMNS.index("crosses_origin")
    assert {r[pcol] for r in seen} == {P, NEAR}
    assert any(r[pcol] == NEAR and r[ccol] == "true" for r in seen)
    assert arrays.producers(_files(both["root"])) == sorted([P, NEAR])


def _files(root: Path):
    from biobox import catalog

    return catalog.ReleaseFiles.of(root, e14.REL)


def test_the_summary_counts_the_arrays_of_each_producer(world, both, capsys):
    uids = e14.anchors(both["tapes"], both["root"])
    policy = nbr.Policy.parse(window="3g")
    base = e14._summary(world[0], uids, "3g", capsys, ["--format", "tsv"])
    new = e14._summary(both["root"], uids, "3g", capsys, ["--format", "tsv"])
    loci = [e14.oracle_table(both, b.split("\n"), policy)[1]
            for b in e14._nbr_tsv(world[0], uids, "3g", capsys).rstrip("\n").split("\n\n")]
    lines = ["", "\t".join(arrays.SUMMARY_COLUMNS)]
    pcol, partcol = arrays.COLUMNS.index("producer"), arrays.COLUMNS.index("partial")
    for p in sorted([P, NEAR]):
        with_ = [rs for rs in loci if any(r[pcol] == p for r in rs)]
        partial_only = [rs for rs in with_ if all(r[partcol] == "true" for r in rs if r[pcol] == p)]
        n_arrays = len({r[1] for rs in loci for r in rs if r[pcol] == p})
        lines.append("\t".join(["crispr_array", p, str(len(loci)), str(len(with_)), str(round(len(with_) / len(loci), 4)),
                                str(len(partial_only)), str(n_arrays)]))
        assert len(with_) > 0, p
    assert new == base + "\n".join(lines) + "\n"
