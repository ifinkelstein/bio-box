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
rows), list tables 1/16 of MAX_SLOTS (8,192 slots), the three contig-level tables 32 rows, and
`token_dict` 1/128 of ROW_GROUP_ROWS (960 rows; about 10,000 tokens at 100x, card B-34). The
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
TOKEN_CAP = bio.ROW_GROUP_ROWS // 128  # token_dict: 960 rows (card B-34)
LIST_TABLES = {"contig_tapes": "n_slots", "tape_tokens": "nr90_tokens"}  # table -> column that counts slots
PROBED = ("seq_locator", "contig_index", "contig_dim", "contig_tapes", "tape_tokens", "token_postings", "token_dict")


def _cap(table: str) -> int:
    return CONTIG_CAP if table in CONTIG_TABLES else TOKEN_CAP if table == "token_dict" else ROW_CAP


def _scaled_write_table(real):
    def write_table(batches, spec, path, *, row_group_rows=bio.ROW_GROUP_ROWS, max_slots=None, **kw):
        cap = _cap(spec.name)
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
    assert not short, (f"fewer than {MIN_RG} row groups at {COPIES}x (caps {ROW_CAP}/{SLOT_CAP}/{CONTIG_CAP}/{TOKEN_CAP}): "
                       f"{short}; all: {n}")


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
                    cap = _cap(t)
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


# --- lookups the exit measurement found NOT pruning, fixed by card B-31 ----------------------------

def _scans_of(q: dict, table: str) -> list[dict]:
    return [s for s in q["scans"] if s["table"] == table]


def _fresh(store, tables: dict | None = None):
    """A fresh connection in place of the store's: empty lookup caches, so counts are cold (synth_scale's
    `fresh()`). The old connection is closed first (since card B-34 two could be open at once); the script's
    tables and `tables` are registered on the new one."""
    g = store["g"]
    store["con"].close()
    store["con"] = con = catalog.open_catalog(g["cat"], threads=4)
    for name, tbl in {**g["registered"], **(tables or {})}.items():
        con.register(name, tbl)
    return con


def _cold(store, sql: str, params=None, keys: int = 1, tables: dict | None = None) -> dict:
    """synth_scale's `run` scoring of `sql` on a fresh connection. Since card B-31 `partners_of` reads tape
    tokens by point lookup with a row-group cache on the connection, so the script's reports of partner
    queries after the first are warm (the family query's contigs hold the NR90 query's): these are re-run cold."""
    out = catalog.profile_query(_fresh(store, tables), sql, params)
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


def _assert_partner_bound(store, sql: str, params: list, tok: int) -> dict:
    """The bound of a neighbour query (TASKS §2, card B-31), run cold: the anchor side reads the postings of
    one token (at most 2 row groups per scan; two scans, `token_slots` and the hub check, and no third one
    for the partner side); the partner side reads `tape_tokens` only, exactly the row groups whose footer
    range holds an anchor contig, so at most one per anchor contig."""
    q = _cold(store, sql, params)
    con = store["con"]
    anchors = sorted(r[0] for r in con.execute(
        "SELECT DISTINCT contig_uid FROM token_postings WHERE release = catalog_release() AND token_id = ?", [tok]).fetchall())
    assert {s["table"] for s in q["scans"]} <= {"token_postings", "tape_tokens", "token_dict"}, q["scans"]
    post = _scans_of(q, "token_postings")
    assert 1 <= len(post) <= 2 and all(s["row_groups"] <= MAX_RG_PER_KEY for s in post), _show(q)
    near = sum(s["row_groups"] for s in _scans_of(q, "tape_tokens"))
    assert near == _footer_row_groups(store, "tape_tokens", anchors) <= len(anchors), (_show(q), len(anchors))
    assert q["result"][0][0] > 0
    return q


# card D-05 first scored these two lookups at 2 row groups per token, which a neighbour query cannot meet: its
# partner side reads the tape tokens of every anchor contig (TASKS §2), one row group per contig at most, and
# the anchors of one token can lie on any number of contigs (100 of 300 for the family token here, 10 for the
# NR90 token, spread by contig_uid over the 11 row groups). So the bound is per anchor contig on the partner
# side and per token on the anchor side; that it does not grow with the release is
# test_partners_of_reads_grow_with_the_anchors_not_the_release.
def test_partners_of_literal_key_reads_per_anchor_contig(store):
    """partners_of with a literal family token (COPIES loci): the bound of `_assert_partner_bound`."""
    tok = store["rep"]["tokens"]["family_a"][1]
    _assert_partner_bound(store, "SELECT count(*), sum(n_pairs) FROM partners_of(?, 10)", [tok], tok)


def test_partners_of_token_of_key_reads_per_anchor_contig(store):
    """partners_of with an NR90 key through token_of (10 loci): the bound of `_assert_partner_bound`."""
    key = store["rep"]["tokens"]["nr90_a"][0]
    tok = store["con"].execute("SELECT token_of('nr90', ?)", [key]).fetchall()[0][0]
    _assert_partner_bound(store, "SELECT count(*), sum(n_pairs) FROM partners_of(token_of('nr90', ?), 10)", [key], tok)


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


# --- token_of resolves a stable key by point lookup (card B-34) ------------------------------------


def _token_oracle(store) -> dict:
    """Oracle: every `token_dict` row of the release read whole with pyarrow, with the row group that holds it
    (numbered in token_id order): `by` maps (level, key) to [(token_id, row group)]; `n` row groups."""
    groups = []
    for p in _files(store, "token_dict"):
        pf = pq.ParquetFile(p)
        groups += [pf.read_row_group(g, columns=["token_id", "level", "key"]).to_pylist()
                   for g in range(pf.metadata.num_row_groups)]
    groups.sort(key=lambda rows: rows[0]["token_id"])
    by: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for i, rows in enumerate(groups):
        for r in rows:
            by.setdefault((r["level"], r["key"]), []).append((r["token_id"], i))
    return {"by": by, "n": len(groups)}


def _token_reads(q: dict) -> int:
    """Row groups of `token_dict` a scored query read: DuckDB's scans and token_of's point lookups."""
    return sum(s["row_groups"] for s in q["scans"] if s["table"] == "token_dict")


@pytest.mark.covers("LY-06")
@pytest.mark.parametrize("level", ["family", "nr50", "nr90"])
def test_token_of_literal_key_reads_at_most_2_row_groups(store, level, record_property):
    """Card B-34: token_of(level, key) with a literal key reads at most 2 of token_dict's row groups (at least
    MIN_RG), and gives the token the oracle reads for the key."""
    o = _token_oracle(store)
    assert o["n"] >= MIN_RG, f"token_dict has {o['n']} row groups at cap {TOKEN_CAP}"
    keys = sorted(k for lv, k in o["by"] if lv == level)
    key = keys[len(keys) // 2]
    (tok, _), = o["by"][(level, key)]
    q = _cold(store, "SELECT token_of(?, ?)", [level, key])
    assert q["result"] == [(tok,)]
    record_property("token_dict_reads", json.dumps({"read": _token_reads(q), "of": o["n"]}))
    assert 1 <= _token_reads(q) <= MAX_RG_PER_KEY, _show(q)


def _keys_in_two_row_groups(o: dict, n: int) -> tuple[list[str], set[int]]:
    """`n` nr90 keys, half from the first row group that holds at least n/2 nr90 keys and half from the last
    such one, and those two row groups. A lookup that prunes reads the two; a scan, or a join filter on the
    keys' min and max, reads every row group between them as well."""
    per: dict[int, list[str]] = {}
    for (lv, k), hits in o["by"].items():
        if lv == "nr90":
            per.setdefault(hits[0][1], []).append(k)
    full = sorted(rg for rg, ks in per.items() if len(ks) >= n // 2)
    a, b = full[0], full[-1]
    assert b - a >= 4, f"the nr90 keys lie in row groups {sorted(per)}: too few for the test to mean something"
    return sorted(per[a])[: n // 2] + sorted(per[b])[-(n - n // 2):], {a, b}


@pytest.mark.covers("LY-06")
def test_token_of_100_joined_keys_read_only_the_row_groups_that_hold_them(store, record_property):
    """Card B-34: token_of for 100 keys from a joined table reads at most 2 row groups per key, and exactly the
    row groups of token_dict that hold the keys (the oracle), which here are 2 of at least MIN_RG."""
    import pyarrow as pa

    o = _token_oracle(store)
    keys, holds = _keys_in_two_row_groups(o, 100)
    want = sum(o["by"][("nr90", k)][0][0] for k in keys)
    tables = {"k100": pa.table({"key": pa.array(keys, pa.string())})}
    q = _cold(store, "SELECT sum(token_of('nr90', key)), count(token_of('nr90', key)) FROM k100", keys=100, tables=tables)
    assert q["result"] == [(want, 100)]
    n = _token_reads(q)
    record_property("token_dict_reads", json.dumps({"read": n, "hold": sorted(holds), "of": o["n"]}))
    assert n <= MAX_RG_PER_KEY * 100 and n == len(holds), f"{n} row groups read; they hold the keys: {sorted(holds)}; {_show(q)}"


@pytest.mark.covers("LY-06")
def test_token_of_1_joined_key_reads_at_most_2_row_groups(store):
    import pyarrow as pa

    o = _token_oracle(store)
    keys, _ = _keys_in_two_row_groups(o, 2)
    q = _cold(store, "SELECT sum(token_of('nr90', key)) FROM k1", tables={"k1": pa.table({"key": pa.array(keys[-1:], pa.string())})})
    assert q["result"] == [(o["by"][("nr90", keys[-1])][0][0],)]
    assert 1 <= _token_reads(q) <= MAX_RG_PER_KEY, _show(q)


def test_token_of_answers_every_key_as_token_dict_holds_it(store):
    """Card B-34: token_of on open_catalog's connection gives, for every (level, key) of token_dict, the token
    the oracle reads, and NULL for an unknown key or level and for a NULL argument."""
    import pyarrow as pa

    o = _token_oracle(store)
    k0 = next(iter(o["by"]))[1]
    rows = [(lv, k, hits[0][0]) for (lv, k), hits in sorted(o["by"].items())]
    rows += [("nr90", "no such key", None), ("nr91", k0, None), (None, k0, None), ("nr90", None, None)]
    assert all(len(h) == 1 for h in o["by"].values())
    tbl = pa.table({"level": pa.array([r[0] for r in rows], pa.string()), "key": pa.array([r[1] for r in rows], pa.string()),
                    "want": pa.array([r[2] for r in rows], pa.uint64())})
    con = _fresh(store, {"allk": tbl})
    assert con.execute("SELECT count(*), count(*) FILTER (WHERE token_of(level, key) IS DISTINCT FROM want) "
                       "FROM allk").fetchall() == [(len(rows), 0)]


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
