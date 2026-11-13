"""Card W-18 (5): `bio-box find`, and `search hmm --positive` resolving through the same function.

Oracle: a loop over the `protein_occurrence` files and the contig labels (`_w07_genes` of the search tests), written
independently of `cli.find_genes`.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from test_micro_s0_search import REL, V, _w07_genes, world  # noqa: E402,F401  (the module-scoped search world)

from biobox import cli  # noqa: E402
from biobox.l3 import tapes  # noqa: E402

pytestmark = [pytest.mark.r2, pytest.mark.covers("W-18")]

# card W-23 appended the gene's span in both forms: 0-based half-open and 1-based inclusive (GenBank)
FIND_COLUMNS = ["accession", "contig_uid", "start", "end", "strand", "feature_uid", "protein_uid", "span",
                "genbank_span"]


def _oracle(root: Path, keep) -> list[list[str]]:
    """The rows `find` must print, in coordinate order within each contig."""
    rows = {(g["contig_uid"], g["feature_uid"]): [g["accession"], g["contig_uid"].hex(), str(g["start"]), str(g["end"]),
                                                  f"{g['strand']:+d}", g["feature_uid"].hex(), g["protein_uid"].hex(),
                                                  f"{g['accession']}:{g['start']}-{g['end']}",
                                                  f"{g['accession']}:{g['start'] + 1}..{g['end']}"]
            for g in _w07_genes(root) if keep(g)}
    return sorted(rows.values(), key=lambda r: (r[0], r[1], int(r[2]), int(r[3]), r[4]))


def _find(root: Path, capsys, *args) -> tuple[int, list[list[str]], str]:
    code = cli.main(["find", "--root", str(root), "--release", REL, *args])
    got = capsys.readouterr()
    lines = [x.split("\t") for x in got.out.splitlines()]
    assert not lines or lines[0] == FIND_COLUMNS, lines[:1]
    return code, lines[1:], got.err


def _first_lambda_gene(root: Path) -> dict:
    return min((g for g in _w07_genes(root) if "phage lambda" in g["header"]), key=lambda g: g["start"])


def test_w18_find_prints_a_gene_by_accession_span_and_by_protein_against_the_occurrence_table(world, capsys):
    root = world[0]
    g = _first_lambda_gene(root)
    acc = g["accession"]
    # one gene: contig, 0-based half-open coordinates, strand, feature_uid and protein uid
    code, rows, err = _find(root, capsys, "--accession", f"{acc}:{g['start']}-{g['end']}")
    assert code == 0, err
    assert rows == _oracle(root, lambda x: x["accession"] == acc and (x["start"], x["end"]) == (g["start"], g["end"]))
    assert len(rows) == 1 and rows[0][5] == g["feature_uid"].hex() and rows[0][6] == g["protein_uid"].hex()
    # the whole contig, in coordinate order
    code, rows, err = _find(root, capsys, "--accession", acc)
    assert code == 0, err
    assert rows == _oracle(root, lambda x: x["accession"] == acc) and len(rows) > 20
    # a protein: every gene that carries it
    code, rows, err = _find(root, capsys, "--protein", g["protein_uid"].hex())
    assert code == 0, err
    assert rows == _oracle(root, lambda x: x["protein_uid"] == g["protein_uid"]) and rows
    # no ordinal is printed (AGENTS 3.3): exactly the columns of FIND_COLUMNS
    assert all(len(r) == len(FIND_COLUMNS) for r in rows)


def test_w18_find_refuses_what_it_cannot_resolve_with_exit_2(world, capsys):
    root = world[0]
    g = _first_lambda_gene(root)
    acc = g["accession"]
    for args, text in [((), "required"), (("--accession", acc, "--protein", "ab" * 16), "not allowed")]:
        with pytest.raises(SystemExit) as e:  # argparse: exactly one of the two
            cli.main(["find", "--root", str(root), "--release", REL, *args])
        assert e.value.code == 2
        assert text in capsys.readouterr().err
    for args, text in [(("--protein", "not-a-uid"), "not-a-uid"), (("--protein", "00" * 16), "no gene"),
                       (("--accession", "NC_999999.9"), "no contig"),
                       (("--accession", f"{acc}:{g['start'] + 1}-{g['end']}"), "no gene")]:
        code, rows, err = _find(root, capsys, *args)
        assert code == 2 and not rows and text in err, (args, code, err)


def test_w18_search_hmm_positive_resolves_through_the_function_find_uses(world, monkeypatch):
    """`search hmm --positive` and `find` are one lookup, not two copies."""
    root = world[0]
    g = _first_lambda_gene(root)
    files, _ = cli._release_files_or_refusal(argparse.Namespace(release=REL), root)
    calls = []
    real = cli.find_genes
    monkeypatch.setattr(cli, "find_genes", lambda *a, **k: calls.append(k) or real(*a, **k))
    uid, lines = cli._positive_uid(root, files, f"{g['accession']}:{g['start']}-{g['end']}", tapes.tape_version([V]), V)
    assert (uid, lines) == (g["protein_uid"].hex(), []) and len(calls) == 1
