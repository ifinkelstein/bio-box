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


def _build(out: Path, copies: int) -> dict:
    """`synth_scale.py` at `copies` into `out`, with the scaled row-group caps."""
    mp = pytest.MonkeyPatch()
    mp.setattr(bio, "write_table", _scaled_write_table(bio.write_table))
    t0 = time.perf_counter()
    argv = sys.argv
    sys.argv = [str(SCRIPT), str(REPO), str(out), str(copies), "--shard-bases", str(2**26), "--seed", str(SEED)]
    try:
        with contextlib.redirect_stdout(stdio.StringIO()) as buf:
            g = runpy.run_path(str(SCRIPT), run_name="__main__")
    finally:
        sys.argv = argv
        mp.undo()
    rep = json.loads(buf.getvalue())
    return {"g": g, "rep": rep, "con": g["con"], "man": g["man"], "root": g["root"], "build_s": time.perf_counter() - t0}


@pytest.fixture(scope="module")
def store(tmp_path_factory):
    s = _build(tmp_path_factory.mktemp("scale"), COPIES)
    yield s
    s["con"].close()


@pytest.fixture(scope="module")
def store2x(tmp_path_factory):
    """The same build at twice the copies: copies 0 to COPIES - 1 are the same contigs (same seed, same order)."""
    s = _build(tmp_path_factory.mktemp("scale2x"), 2 * COPIES)
    yield s
    s["con"].close()


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


def _fresh(store):
    """A fresh connection in place of the store's: empty lookup caches, so counts are cold (synth_scale's
    `fresh()`). One open catalog per database file in a process (its functions are the instance's), so the
    old connection is closed first; the script's tables are registered again."""
    g = store["g"]
    store["con"].close()
    store["con"] = con = catalog.open_catalog(g["cat"], threads=4)
    for name, tbl in g["registered"].items():
        con.register(name, tbl)
    return con


def _cold(store, sql: str, params=None, keys: int = 1) -> dict:
    """synth_scale's `run` scoring of `sql` on a fresh connection. Since card B-31 `partners_of` reads tape
    tokens by point lookup with a row-group cache on the connection, so the script's reports of partner
    queries after the first are warm (the family query's contigs hold the NR90 query's): these are re-run cold."""
    out = catalog.profile_query(_fresh(store), sql, params)
    scans = [{"table": _table_of(s["files"][0]), "row_groups": s["row_groups_scanned"], "of": s["row_groups_total"],
              "rg_per_key": s["row_groups_scanned"] / keys} for s in out["scans"] if s["files"]]
    idx = out["index"]
    return {"keys": keys, "scans": scans, "index": idx, "result": out["result"],
            "max_rg_per_key": max([s["rg_per_key"] for s in scans] + [sum(idx.values()) / keys])}


def test_token_of_key_prunes_postings(store):
    q = _rows(store["rep"], "postings of a token: key from token_dict (token_of), nr90")
    assert _scans_of(q, "token_postings")
    _assert_prunes(q)


def test_q1_token_of_keys_prune(store):
    _assert_prunes(_rows(store["rep"], "neighbours: Q1 nr90 pair, keys from token_dict (token_of)"))


def test_q4_token_of_keys_prune(store):
    _assert_prunes(_rows(store["rep"], "neighbours: Q4 nr90 pair instances, keys from token_dict (token_of)"))


# The two partners_of lookups below keep their markers after card B-31. Their key is one token, but the partner
# side reads the tape tokens of every anchor contig (TASKS §2), one row group per contig at most: 100 contigs of
# 300 for the family token, 10 for the NR90 token, spread by contig_uid over the 11 row groups of tape_tokens, so
# no layout can serve them in 2 row groups per token. What B-31 asks of them (reads in proportion to the anchors'
# contigs, not to the release) is tested by test_partners_of_reads_grow_with_the_anchors_not_the_release.
@pytest.mark.xfail(strict=True, raises=AssertionError, reason=B31)
def test_partners_of_near_side_prunes_literal_key(store):
    q = _cold(store, "SELECT count(*), sum(n_pairs) FROM partners_of(?, 10)", [store["rep"]["tokens"]["family_a"][1]])
    assert len(_scans_of(q, "token_postings")) >= 2, "the anchor and the near side are separate scans"
    _assert_prunes(q)


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=B31)
def test_partners_of_token_of_key_prunes(store):
    _assert_prunes(_cold(store, "SELECT count(*), sum(n_pairs) FROM partners_of(token_of('nr90', ?), 10)",
                         [store["rep"]["tokens"]["nr90_a"][0]]))


def _macro(store, sql: str) -> dict:
    return catalog.profile_query(_fresh(store), sql)  # cold: the point lookups cache row groups per connection


def _macro_prunes(store, sql: str) -> None:
    """Row groups read per key, by Parquet scans and by the point lookups of the catalog's Python functions
    (`profile_query`'s `index`: `seq_locator` and `contig_index`, card B-31 as B-14)."""
    out = _macro(store, sql)
    reads = [s["row_groups_scanned"] for s in out["scans"]] + list(out["index"].values())
    assert sum(reads) > 0, "the query must read row groups for the count to mean something"
    worst = max(reads)
    assert worst <= MAX_RG_PER_KEY, "; ".join(f"{s['row_groups_scanned']} of {s['row_groups_total']}" for s in out["scans"]) + \
        f"; index {out['index']}"


def test_default_digest_joined_key_prunes(store):
    _macro_prunes(store, "SELECT count(default_digest(contig_uid)) FROM pick1")


def test_digest_len_joined_key_prunes(store):
    _macro_prunes(store, "SELECT sum(digest_len(default_digest(contig_uid))) FROM pick1")


def test_contig_len_of_joined_key_prunes(store):
    _macro_prunes(store, "SELECT sum(contig_len_of(contig_uid)) FROM pick1")


def _table_of(path: str) -> str:
    parts = Path(path).parts
    return parts[parts.index("l3") + 1] if "l3" in parts else Path(path).parent.name


def _footer_row_groups(store, table: str, keys: list[bytes]) -> int:
    """Oracle: the row groups of `table` whose contig_uid min/max holds one of `keys`, from the Parquet footers."""
    n = 0
    for p in _files(store, table):
        md = pq.ParquetFile(p).metadata
        ci = next(i for i in range(md.num_columns) if md.schema.column(i).path == "contig_uid")
        for g in range(md.num_row_groups):
            st = md.row_group(g).column(ci).statistics
            n += any(st.min <= k <= st.max for k in keys)
    return n


def _partner_reads(s: dict) -> dict:
    """partners_of(<NR90 token of copies 0-9>, 10) on a fresh connection (cold lookup caches): the row groups
    each table was read for, the anchors' contigs, and the result with partner tokens resolved to stable keys."""
    g = s["g"]
    key = g["nr90_key"](g["pidx"][g["gl"][g["gi"]][3]], 0)  # copies 0-9: the same contigs at every scale
    con = _fresh(s)
    tok =con.execute("SELECT token_of('nr90', ?)", [key]).fetchall()[0][0]  # resolved first, passed as a literal
    prof = catalog.profile_query(con, "SELECT partner_token, n_loci, n_pairs FROM partners_of(?, 10)", [tok])
    anchors = sorted(r[0] for r in con.execute(
        "SELECT DISTINCT contig_uid FROM token_postings WHERE release = catalog_release() AND token_id = ?", [tok]).fetchall())
    keyed = sorted(con.execute(
        "SELECT d.level, d.key, p.n_loci, p.n_pairs FROM partners_of(?, 10) p JOIN token_dict d "
        "ON d.token_id = p.partner_token AND d.release = catalog_release()", [tok]).fetchall())
    reads: dict[str, list[tuple[int, int]]] = {}
    for sc in prof["scans"]:
        if sc["files"]:
            reads.setdefault(_table_of(sc["files"][0]), []).append((sc["row_groups_scanned"], sc["row_groups_total"]))
    return {"reads": reads, "anchors": anchors, "keyed": keyed, "n_rows": len(prof["result"]),
            "total": {t: sum(f["row_groups"] for f in s["man"]["files"] if f["table"] == t) for t in ("tape_tokens", "token_postings")}}


@pytest.mark.covers("Q3", "LY-06")
def test_partners_of_reads_grow_with_the_anchors_not_the_release(store, store2x, record_property):
    """Card B-31: the partner side of partners_of reads the anchors' contigs' row groups of `tape_tokens` (exactly
    the row groups whose footer range holds an anchor contig, at most one per contig), whatever the size of the
    release; the anchors' postings stay at most 2 row groups; doubling the store changes neither the reads' bound
    nor the answer (partner tokens resolved to their stable keys through token_dict)."""
    one, two = _partner_reads(store), _partner_reads(store2x)
    assert one["anchors"] == two["anchors"] and len(one["anchors"]) == 10
    assert one["keyed"] == two["keyed"] and one["n_rows"] == len(one["keyed"]) > 0
    assert all(two["total"][t] >= 1.8 * one["total"][t] for t in one["total"]), (one["total"], two["total"])
    for s, r in ((store, one), (store2x, two)):
        assert set(r["reads"]) <= {"token_postings", "tape_tokens"}, r["reads"]
        assert all(n <= MAX_RG_PER_KEY for n, _ in r["reads"]["token_postings"]), r["reads"]
        near = sum(n for n, _ in r["reads"].get("tape_tokens", []))
        assert near == _footer_row_groups(s, "tape_tokens", r["anchors"]) <= len(r["anchors"]), r["reads"]
    record_property("partner_reads", json.dumps({"1x": one["reads"], "2x": two["reads"], "totals": [one["total"], two["total"]]}))


def test_build_time(store, record_property):
    record_property("build_s", round(store["build_s"], 1))
    assert store["build_s"] < 300, f"{store['build_s']:.0f} s for the {COPIES}x build and queries (card D-05: under 5 minutes)"
