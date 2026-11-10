"""Ring R2: the layout-at-scale measurement as a test (card D-05; AGENTS §4.4 and §4.6).

`scripts/stress/synth_scale.py` is run unchanged (by path, as `__main__`, with its command line) on
the micro-S0 fixture at COPIES = 100 in `tmp_path`: sequences through `ingest`, gene features
replicated, synthetic clusters and families, then the repo's own tape, token, postings and
`contig_dim` builders, the manifest, the catalog and the key lookups scored with
`catalog.profile_query` (2026-11-10 phase B exit, section 2).

Scale of the row groups. At the real limits (122,880 rows, 131,072 slots) a 100x store has one or two
row groups per table, and 8 would take about 1,100x for the tapes and about 10^6 contigs for
`contig_index` and `seq_locator` (data_layer_and_tests.md and the exit measurement: 9 of 9 at 10^6
contigs), which no test can build. So the build runs with the row-group caps of `io.write_table`
scaled down, test side only (nothing in `src/` changes): row tables 1/16 of ROW_GROUP_ROWS (7,680
rows), list tables 1/16 of MAX_SLOTS (8,192 slots), and the three contig-level tables 32 rows. The
pruning behaviour (zone maps on the sort key) is what is under test; the assertions state the caps
they ran with and also the real AGENTS bounds. Every probed table must have at least MIN_RG row
groups, or the test fails before it measures anything.
"""

from __future__ import annotations

import contextlib
import io as stdio
import json
import runpy
import sys
import time
from pathlib import Path

import pyarrow.compute as pc
import pyarrow.parquet as pq
import pytest

from biobox import catalog
from biobox import io as bio
from biobox.l3 import tapes

pytestmark = [pytest.mark.r2, pytest.mark.covers("D-05", "LY-06")]

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts/stress/synth_scale.py"
COPIES = 100
SEED = 1
MIN_RG = 8  # row groups a probed table must have for pruning to mean something
MAX_RG_PER_KEY = 2  # AGENTS §4.4, as the exit measurement states it
ROW_CAP = bio.ROW_GROUP_ROWS // 16  # 7,680 rows
SLOT_CAP = tapes.MAX_SLOTS // 16  # 8,192 slots
CONTIG_CAP = 32  # contig_index, seq_locator, contig_dim: 300 contigs at 100x
CONTIG_TABLES = ("contig_index", "seq_locator", "contig_dim")
LIST_TABLES = {"contig_tapes": "n_slots", "tape_tokens": "nr90_tokens"}  # table -> column that counts slots
PROBED = ("seq_locator", "contig_index", "contig_dim", "contig_tapes", "tape_tokens", "token_postings")


def _scaled_write_table(real):
    def write_table(batches, spec, path, *, row_group_rows=bio.ROW_GROUP_ROWS, max_slots=None, **kw):
        cap = CONTIG_CAP if spec.name in CONTIG_TABLES else ROW_CAP
        if max_slots is not None:
            max_slots = min(max_slots, SLOT_CAP)
        return real(batches, spec, path, row_group_rows=min(row_group_rows, cap), max_slots=max_slots, **kw)
    return write_table


@pytest.fixture(scope="module")
def store(tmp_path_factory):
    out = tmp_path_factory.mktemp("scale")
    mp = pytest.MonkeyPatch()
    mp.setattr(bio, "write_table", _scaled_write_table(bio.write_table))
    t0 = time.perf_counter()
    argv = sys.argv
    sys.argv = [str(SCRIPT), str(REPO), str(out), str(COPIES), "--shard-bases", str(2**26), "--seed", str(SEED)]
    try:
        with contextlib.redirect_stdout(stdio.StringIO()) as buf:
            g = runpy.run_path(str(SCRIPT), run_name="__main__")
    finally:
        sys.argv = argv
        mp.undo()
    rep = json.loads(buf.getvalue())
    yield {"g": g, "rep": rep, "con": g["con"], "man": g["man"], "root": g["root"], "build_s": time.perf_counter() - t0}
    g["con"].close()


def _rows(rep: dict, name: str) -> dict:
    return rep["queries"][name]


def _scans(q: dict) -> list[dict]:
    return q["scans"]


def _show(q: dict) -> str:
    return "; ".join(f"{s['table']} {s['row_groups']} of {s['of']}" for s in q["scans"]) + \
        f"; index {q['index']['locator_row_groups_read']}+{q['index']['index_row_groups_read']}"


def _assert_prunes(q: dict) -> None:
    got = q["max_rg_per_key"]
    assert got <= MAX_RG_PER_KEY, f"{got} row groups per key (> {MAX_RG_PER_KEY}); read: {_show(q)}"


# --- the store has enough row groups for pruning to mean something and respects AGENTS 4.6 --------


def _files(store, table: str) -> list[Path]:
    root = store["root"]
    return [root / f["path"] for f in store["man"]["files"] if f["table"] == table]


def test_probed_tables_have_enough_row_groups(store):
    n = {t: sum(f["row_groups"] for f in store["man"]["files"] if f["table"] == t) for t in PROBED}
    short = {t: k for t, k in n.items() if k < MIN_RG}
    assert not short, f"fewer than {MIN_RG} row groups at {COPIES}x (caps {ROW_CAP}/{SLOT_CAP}/{CONTIG_CAP}): {short}; all: {n}"


def test_row_group_sizes_agents_4_6(store):
    assert (bio.ROW_GROUP_ROWS, tapes.MAX_SLOTS) == (122_880, 131_072)
    tables = {f["table"] for f in store["man"]["files"]}
    checked = {}
    for t in sorted(tables):
        for p in _files(store, t):
            pf = pq.ParquetFile(p)
            for i in range(pf.metadata.num_row_groups):
                rows = pf.metadata.row_group(i).num_rows
                if t in LIST_TABLES:
                    col = pf.read_row_group(i, columns=[LIST_TABLES[t]]).column(0)
                    slots = pc.sum(col).as_py() if t == "contig_tapes" else pc.sum(pc.list_value_length(col)).as_py()
                    assert slots <= tapes.MAX_SLOTS, f"{t} {p.name} row group {i}: {slots} slots > 131,072"
                    assert slots <= SLOT_CAP or rows == 1, f"{t} {p.name} row group {i}: {slots} slots > cap {SLOT_CAP}"
                    checked[t] = max(checked.get(t, 0), slots)
                else:
                    assert rows <= bio.ROW_GROUP_ROWS, f"{t} {p.name} row group {i}: {rows} rows > 122,880"
                    cap = CONTIG_CAP if t in CONTIG_TABLES else ROW_CAP
                    assert rows <= cap, f"{t} {p.name} row group {i}: {rows} rows > cap {cap}"
                    checked[t] = max(checked.get(t, 0), rows)
    assert set(LIST_TABLES) <= set(checked)


# --- AGENTS 4.4: lookups the exit measurement found pruning --------------------------------------

PRUNING = [
    ("sequence by contig: literal key", 1),
    ("sequence by contig: 1 keys from a joined table", 1),
    ("sequence by contig: 10 keys from a joined table", 10),
    ("sequence by contig: 1000 keys from a joined table", 1000),
    ("sequence by contig: 1,000 keys from contig_dim joined to the pick table", 1000),
    ("tape of a contig: literal key", 1),
    ("tape of a contig: 1 keys from a joined table (tapes_of)", 1),
    ("tape of a contig: 10 keys from a joined table (tapes_of)", 10),
    ("tape of a contig: 1000 keys from a joined table (tapes_of)", 1000),
    ("tape of a contig: 1 key, plain join to contig_tapes", 1),
    ("postings of a token: literal key, family (COPIES loci)", 1),
    ("postings of a token: literal key, nr90 (10 loci)", 1),
    ("postings of a token: key from token_dict, plain join, nr90", 1),
    ("postings of a token: 100 nr90 keys from a joined table, plain join", 100),
]


@pytest.mark.parametrize("name,keys", PRUNING, ids=[p[0] for p in PRUNING])
def test_lookup_prunes(store, name, keys):
    q = _rows(store["rep"], name)
    assert q["keys"] == keys
    _assert_prunes(q)


# --- lookups the exit measurement found NOT pruning: card B-31 ------------------------------------

B31 = "card B-31"


def _scans_of(q: dict, table: str) -> list[dict]:
    return [s for s in q["scans"] if s["table"] == table]


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=B31)
def test_token_of_key_prunes_postings(store):
    q = _rows(store["rep"], "postings of a token: key from token_dict (token_of), nr90")
    assert _scans_of(q, "token_postings")
    _assert_prunes(q)


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=B31)
def test_q1_token_of_keys_prune(store):
    _assert_prunes(_rows(store["rep"], "neighbours: Q1 nr90 pair, keys from token_dict (token_of)"))


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=B31)
def test_q4_token_of_keys_prune(store):
    _assert_prunes(_rows(store["rep"], "neighbours: Q4 nr90 pair instances, keys from token_dict (token_of)"))


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=B31)
def test_partners_of_near_side_prunes_literal_key(store):
    q = _rows(store["rep"], "neighbours: Q3 partners of one family, literal key")
    assert len(_scans_of(q, "token_postings")) >= 2, "the anchor and the near side are separate scans"
    _assert_prunes(q)


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=B31)
def test_partners_of_token_of_key_prunes(store):
    _assert_prunes(_rows(store["rep"], "neighbours: Q3 partners of one nr90, key from token_dict (token_of)"))


def _macro(store, sql: str) -> dict:
    return catalog.profile_query(store["con"], sql)


def _macro_prunes(store, sql: str) -> None:
    out = _macro(store, sql)
    assert out["scans"], "the query must scan Parquet for the row group count to mean something"
    worst = max(s["row_groups_scanned"] for s in out["scans"])
    assert worst <= MAX_RG_PER_KEY, "; ".join(f"{s['row_groups_scanned']} of {s['row_groups_total']}" for s in out["scans"])


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=B31)
def test_default_digest_joined_key_prunes(store):
    _macro_prunes(store, "SELECT count(default_digest(contig_uid)) FROM pick1")


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=B31)
def test_digest_len_joined_key_prunes(store):
    _macro_prunes(store, "SELECT sum(digest_len(default_digest(contig_uid))) FROM pick1")


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=B31)
def test_contig_len_of_joined_key_prunes(store):
    _macro_prunes(store, "SELECT sum(contig_len_of(contig_uid)) FROM pick1")


def test_build_time(store, record_property):
    record_property("build_s", round(store["build_s"], 1))
    assert store["build_s"] < 300, f"{store['build_s']:.0f} s for the {COPIES}x build and queries (card D-05: under 5 minutes)"
