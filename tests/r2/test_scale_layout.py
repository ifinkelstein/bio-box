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
rows), list tables 1/16 of MAX_SLOTS (8,192 slots), the contig-level tables 32 rows (`contig_index`,
`seq_locator`, `contig_dim`, and since card B-27 `contig_occurrence` and `contig_label`), and
`token_dict` 1/128 of ROW_GROUP_ROWS (960 rows; about 10,000 tokens at 100x, card B-34). The
pruning behaviour (zone maps on the sort key) is what is under test; the assertions state the caps
they ran with and also the real AGENTS bounds. Every probed table must have at least MIN_RG row
groups, or the test fails before it measures anything.

Card B-27: the label and membership lookups. `synth_scale.py` writes `cluster_member` by hand (no
MMseqs2), so the `lookups` fixture writes the release's `member_cluster` from it with the builder's
own function (`cluster.write_member_cluster`, which `cluster.run` calls) and builds a second manifest
and catalog; `contig_label` is written by ingest inside the script.
"""

from __future__ import annotations

import contextlib
import io as stdio
import json
import runpy
import sys
import time
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pytest

from biobox import catalog
from biobox import io as bio
from biobox.l3 import tapes
from biobox.registry.schemas import conformance_problems, table

pytestmark = [pytest.mark.r2, pytest.mark.covers("D-05", "LY-06")]

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts/stress/synth_scale.py"
COPIES = 100
SEED = 1
MIN_RG = 8  # row groups a probed table must have for pruning to mean something
MAX_RG_PER_KEY = 2  # AGENTS §4.4, as the exit measurement states it
ROW_CAP = bio.ROW_GROUP_ROWS // 16  # 7,680 rows
SLOT_CAP = tapes.MAX_SLOTS // 16  # 8,192 slots
CONTIG_CAP = 32  # contig_index, seq_locator, contig_dim, contig_occurrence, contig_label: 300 contigs at 100x
CONTIG_TABLES = ("contig_index", "seq_locator", "contig_dim", "contig_occurrence", "contig_label")
TOKEN_CAP = bio.ROW_GROUP_ROWS // 128  # token_dict: 960 rows (card B-34)
LIST_TABLES = {"contig_tapes": "n_slots", "tape_tokens": "nr90_tokens"}  # table -> column that counts slots
PROBED = ("seq_locator", "contig_index", "contig_dim", "contig_tapes", "tape_tokens", "token_postings", "token_dict",
          "contig_occurrence", "contig_label")


def _cap(table: str) -> int:
    return CONTIG_CAP if table in CONTIG_TABLES else TOKEN_CAP if table == "token_dict" else ROW_CAP


def _scaled_write_table(real, slot_caps: dict | None = None):
    """`real` with the scaled caps; `slot_caps` overrides SLOT_CAP for the list tables it names (card B-38)."""
    def write_table(batches, spec, path, *, row_group_rows=bio.ROW_GROUP_ROWS, max_slots=None, **kw):
        cap = _cap(spec.name)
        if max_slots is not None:
            max_slots = min(max_slots, (slot_caps or {}).get(spec.name, SLOT_CAP))
        return real(batches, spec, path, row_group_rows=min(row_group_rows, cap), max_slots=max_slots, **kw)
    return write_table


def _build(out: Path, copies: int, slot_caps: dict | None = None) -> dict:
    """`synth_scale.py` at `copies` into `out`, with the scaled row-group caps."""
    mp = pytest.MonkeyPatch()
    mp.setattr(bio, "write_table", _scaled_write_table(bio.write_table, slot_caps))
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


# --- card B-35: what the partner side reads, with few anchor contigs and with many -------------------

PARTNERS_SQL = "SELECT partner_token, n_loci, n_pairs FROM partners_of(?, 10)"


def _lookup_scans(prof: dict, table: str) -> list[dict]:
    """The profiled scans of `table` that a point lookup of the catalog made (`"lookup"` set)."""
    return [s for s in prof["scans"] if s.get("lookup") and s["files"] and _table_of(s["files"][0]) == table]


def _anchor_contigs(con, tok: int) -> list[bytes]:
    return sorted(r[0] for r in con.execute(
        "SELECT DISTINCT contig_uid FROM token_postings WHERE release = catalog_release() AND token_id = ?", [tok]).fetchall())


@pytest.mark.covers("Q3", "LY-06")
def test_b35_partner_reads_served_by_the_cache_are_counted(store, record_property):
    """Card B-35: `profile_query` counts the `tape_tokens` row groups a `partners_of` query touches, read from the
    file or served by the lookup's row-group cache. A query with 10 anchor contigs touches exactly the row groups
    whose footer range holds one of them (the oracle), at most 10, cold on a fresh connection and warm after a query
    that cached row groups on the same connection; warm, the cached ones are counted as such, and the answer is the
    same. The same holds for `token_of`'s lookup of `token_dict`."""
    fam, tok = store["rep"]["tokens"]["family_a"][1], store["rep"]["tokens"]["nr90_a"][1]
    key = store["rep"]["tokens"]["nr90_a"][0]
    con = _fresh(store)
    anchors = _anchor_contigs(con, tok)
    assert len(anchors) == 10
    hold = _footer_row_groups(store, "tape_tokens", anchors)
    cold = catalog.profile_query(con, PARTNERS_SQL, [tok])
    many = catalog.profile_query(con, PARTNERS_SQL, [fam])  # caches the row groups of the family's 100 contigs
    warm = catalog.profile_query(con, PARTNERS_SQL, [tok])
    assert warm["result"] == cold["result"] and len(cold["result"]) > 0
    got = {name: _lookup_scans(p, "tape_tokens") for name, p in (("cold", cold), ("many", many), ("warm", warm))}
    record_property("b35_tape_tokens", json.dumps({n: [{k: s.get(k) for k in ("row_groups_scanned", "row_groups_read",
                                                    "row_groups_cached", "row_groups_total")} for s in v]
                                                   for n, v in got.items()} | {"hold": hold}))
    for name in ("cold", "warm"):
        assert len(got[name]) == 1 and got[name][0]["row_groups_scanned"] == hold <= 10, (name, got[name], hold)
    for name in ("cold", "warm"):
        (s,) = got[name]
        assert s["row_groups_read"] + s["row_groups_cached"] == hold, (name, s, hold)
    assert got["cold"][0]["row_groups_cached"] == 0
    assert got["warm"][0]["row_groups_cached"] >= 1, got["warm"]
    # token_of: its row group is cached after the first lookup on the connection; a repeat still reads it
    t_cold = _lookup_scans(catalog.profile_query(con, "SELECT token_of('nr90', ?)", [key]), "token_dict")
    t_warm = _lookup_scans(catalog.profile_query(con, "SELECT token_of('nr90', ?)", [key]), "token_dict")
    assert len(t_warm) == 1 and 1 <= t_warm[0]["row_groups_scanned"] <= MAX_RG_PER_KEY, (t_cold, t_warm)
    assert t_warm[0]["row_groups_cached"] == t_warm[0]["row_groups_scanned"] and t_warm[0]["row_groups_read"] == 0


def _partners_oracle(store, tok: int, k: int) -> list[tuple[int, int, int]]:
    """Oracle for partners_of(tok, k), by brute force over the release's files read whole with pyarrow: the anchor
    slots are the (contig, ordinal) of every `token_postings` row of `tok`; a partner is every non-null token at
    every level of `tape_tokens` on a slot within k of an anchor (`pairs_oracle.window`, slot count and topology
    from `contig_tapes`); (token, loci, pairs) ranked as the macro ranks them."""
    from collections import Counter, defaultdict

    from oracles.pairs_oracle import window

    anchors: dict[bytes, set[int]] = defaultdict(set)
    for p in _files(store, "token_postings"):
        t = pq.read_table(p, columns=["token_id", "contig_uid", "ordinal"], partitioning=None)
        t = t.filter(pc.equal(t.column("token_id"), pa.scalar(tok, pa.uint64())))
        for cu, o in zip(t.column("contig_uid").to_pylist(), t.column("ordinal").to_pylist()):
            anchors[cu].add(o)
    shape = {}
    for p in _files(store, "contig_tapes"):
        for r in pq.read_table(p, columns=["contig_uid", "n_slots", "is_circular"], partitioning=None).to_pylist():
            shape[r["contig_uid"]] = (r["n_slots"], r["is_circular"])
    loci: dict[int, set] = defaultdict(set)
    pairs: Counter = Counter()
    for p in _files(store, "tape_tokens"):
        for r in pq.read_table(p, partitioning=None).to_pylist():
            cu = r["contig_uid"]
            if cu not in anchors:
                continue
            n, circ = shape[cu]
            cols = [r[f"{lv}_tokens"] for lv in ("family", "nr50", "nr90")]
            for i in sorted(anchors[cu]):
                for j, _ in window(i, n, circ, k):
                    for c in cols:
                        if c[j] is not None:
                            loci[c[j]].add((cu, i))
                            pairs[c[j]] += 1
    return sorted(((b, len(s), pairs[b]) for b, s in loci.items()), key=lambda g: (-g[1], g[0]))


@pytest.mark.covers("Q3", "LY-06")
def test_b35_partner_reads_with_many_anchor_contigs_are_every_row_group(store, record_property):
    """Card B-35, the bound when the anchor contigs are many: the partner side touches exactly the `tape_tokens` row
    groups whose footer range holds an anchor contig (the oracle), never more than there are anchor contigs or row
    groups. Contig uids are random, so C anchor contigs over R row groups leave about R * exp(-C / R) row groups
    untouched (derived), fewer than one once C > R ln R: then every row group is read, as the 10 Gbp rebuild
    measured (19 of 19, 26 of 26). The family token here lies on 100 contigs. The answer equals a brute-force count."""
    import math

    fam = store["rep"]["tokens"]["family_a"][1]
    con = _fresh(store)
    anchors = _anchor_contigs(con, fam)
    prof = catalog.profile_query(con, PARTNERS_SQL, [fam])
    (s,) = _lookup_scans(prof, "tape_tokens")
    hold, total = _footer_row_groups(store, "tape_tokens", anchors), s["row_groups_total"]
    record_property("b35_many", json.dumps({"anchor_contigs": len(anchors), "read": s["row_groups_scanned"], "hold": hold,
                                            "of": total, "expected_untouched": total * math.exp(-len(anchors) / total)}))
    assert total >= MIN_RG and len(anchors) > total * math.log(total), (len(anchors), total)
    assert s["row_groups_scanned"] == hold == total, (s, hold)
    assert prof["result"] == _partners_oracle(store, fam, 10)


# --- card B-38: a many-anchor partners_of reads each tape_tokens row group once; its time against the join ------

B38_SLOT_CAP = 1024  # tape_tokens slots per row group in the B-38 store, so the table has more row groups than the
# partner lookup's cache holds (the fixture's three contigs carry about 850 slots per copy, so about 85 row groups
# at 100x, *derived*; the test checks the count before it measures)
B38_RUNS = 5  # timed runs per query and token, each on a fresh connection; the median is compared
B38_TIME_MARGIN = 1.25  # partners_of's median may be this many times the join's (a timing tolerance, not a measurement)
B38_K = 10
# The same answer as partners_of(?, 10) by plain SQL joins, for the timing (card B-38). TAPE_JOIN: the anchors'
# windows joined to tape_tokens on contig_uid (DuckDB's own scan, which prunes by the join filter). POSTINGS_JOIN:
# the ladder's form, the postings self-join on contig_uid (r10_queries.py), at every level; it equals the macro
# here because every synthetic protein has at most one domain (card B-39 is about the difference).
TAPE_JOIN_SQL = f"""
    WITH _anchor AS (SELECT * FROM token_slots(?)),
         _off AS (SELECT contig_uid, ordinal::BIGINT AS _i, n_slots::BIGINT AS _n, is_circular AS _c,
                         unnest(range(-least({B38_K}, n_slots::BIGINT), least({B38_K}, n_slots::BIGINT) + 1)) AS _d
                  FROM _anchor),
         _win AS (SELECT DISTINCT contig_uid, _i, _n, CASE WHEN _c THEN ((_i + _d) % _n + _n) % _n ELSE _i + _d END AS _j
                  FROM _off WHERE _d <> 0),
         _near AS (SELECT w.contig_uid, w._i,
                          unnest([t.family_tokens[w._j + 1], t.nr50_tokens[w._j + 1], t.nr90_tokens[w._j + 1]]) AS _tok
                   FROM _win w JOIN tape_tokens t ON t.contig_uid = w.contig_uid
                   WHERE w._j >= 0 AND w._j < w._n AND w._j <> w._i
                     AND t.release = catalog_release() AND t.tape_version = catalog_tape_version())
    SELECT _tok AS partner_token, count(DISTINCT (contig_uid, _i)) AS n_loci, count(*) AS n_pairs
    FROM _near WHERE _tok IS NOT NULL GROUP BY _tok ORDER BY n_loci DESC, partner_token"""
POSTINGS_JOIN_SQL = f"""
    WITH a AS (SELECT DISTINCT contig_uid, ordinal, n_slots, is_circular FROM token_rows(?)),
         p AS (SELECT DISTINCT token_id, contig_uid, ordinal FROM token_postings
               WHERE release = catalog_release() AND tape_version = catalog_tape_version()
                 AND contig_uid IN (SELECT contig_uid FROM a)),
         nb AS (SELECT a.contig_uid, a.ordinal AS _i, p.token_id AS _tok FROM a JOIN p USING (contig_uid)
                WHERE slot_dist(a.ordinal, p.ordinal, a.n_slots, a.is_circular) BETWEEN 1 AND {B38_K})
    SELECT _tok AS partner_token, count(DISTINCT (contig_uid, _i)) AS n_loci, count(*) AS n_pairs
    FROM nb GROUP BY _tok ORDER BY n_loci DESC, partner_token"""


@pytest.fixture(scope="module")
def store_b38(tmp_path_factory):
    """The 100x store with `tape_tokens` cut at B38_SLOT_CAP slots per row group (every other cap as `store`)."""
    s = _build(tmp_path_factory.mktemp("scale_b38"), COPIES, {"tape_tokens": B38_SLOT_CAP})
    yield s
    s["con"].close()


def _tape_lookup(con):
    """The connection's `tape_tokens` point lookup (built on first use)."""
    return dict(catalog._LOOKUPS[id(con)])["bb_partner_counts"].get()


def _b38_tokens(s: dict) -> dict[str, int]:
    t = s["rep"]["tokens"]
    return {"nr90 (10 loci)": t["nr90_a"][1], "family (COPIES loci)": t["family_a"][1], "hub": t["hub"]}


@pytest.mark.covers("Q3", "LY-06")
def test_b38_many_anchor_partners_of_reads_each_row_group_once(store_b38, record_property):
    """Card B-38: `tape_tokens` has more row groups than the partner lookup's cache holds, and the anchors of the
    family token and of the hub touch more of them than that. Cold, a `partners_of` query decodes each row group it
    touches exactly once (`row_groups_read` counts every decode, card B-35); warm, on the same connection, each is
    decoded at most once or served by the cache. The touched row groups are exactly those whose footer range holds
    an anchor contig, and the answer equals the B-35 oracle."""
    rec = {}
    for name in ("family (COPIES loci)", "hub"):
        tok = _b38_tokens(store_b38)[name]
        con = _fresh(store_b38)
        anchors = _anchor_contigs(con, tok)
        hold = _footer_row_groups(store_b38, "tape_tokens", anchors)
        cold = catalog.profile_query(con, PARTNERS_SQL, [tok])
        warm = catalog.profile_query(con, PARTNERS_SQL, [tok])
        cache = _tape_lookup(con).cache_row_groups
        (c,), (w,) = _lookup_scans(cold, "tape_tokens"), _lookup_scans(warm, "tape_tokens")
        rec[name] = {"anchor_contigs": len(anchors), "hold": hold, "cache": cache,
                     **{f"{n}_{k}": s[k] for n, s in (("cold", c), ("warm", w))
                        for k in ("row_groups_scanned", "row_groups_read", "row_groups_cached", "row_groups_total")}}
        assert c["row_groups_total"] > cache and hold > cache, rec[name]  # the case under test
        assert c["row_groups_scanned"] == hold and c["row_groups_read"] == hold, rec[name]
        assert w["row_groups_scanned"] == hold and w["row_groups_read"] + w["row_groups_cached"] == hold, rec[name]
        assert cold["result"] == warm["result"] == _partners_oracle(store_b38, tok, B38_K)
    record_property("b38_reads", json.dumps(rec))
    print("b38_reads", json.dumps(rec))  # shown with pytest -s


@pytest.mark.covers("Q3", "LY-06")
def test_b38_many_anchor_partners_of_is_no_slower_than_the_join(store_b38, record_property):
    """Card B-38: wall time of `partners_of` against the same answer by plain SQL joins (TAPE_JOIN_SQL,
    POSTINGS_JOIN_SQL), each run B38_RUNS times on a fresh connection (cold lookup caches), alternating, for anchor
    sets of 10, COPIES and the hub's loci. All three give the oracle's answer. With the most anchors (the hub) the
    macro's median is at most B38_TIME_MARGIN times the median of TAPE_JOIN_SQL, the join that answers from the
    same table; the postings self-join is recorded for information only (since card B-21 it answers from other
    rows, card B-39). Wall time on a shared host varies by more than the 1% the two were apart in one run, hence
    the margin; the deterministic acceptance of card B-38 is the row-group counts of
    test_b38_many_anchor_partners_of_reads_each_row_group_once (each touched row group read once), which stay
    exact. Every median is recorded (`b38_times`)."""
    import statistics

    out = {}
    for name, tok in _b38_tokens(store_b38).items():
        times: dict[str, list[float]] = {"partners_of": [], "tape_join": [], "postings_join": []}
        got = {}
        for _ in range(B38_RUNS):
            for q, sql in (("partners_of", PARTNERS_SQL), ("tape_join", TAPE_JOIN_SQL), ("postings_join", POSTINGS_JOIN_SQL)):
                con = _fresh(store_b38)
                t = time.perf_counter()
                got[q] = con.execute(sql, [tok]).fetchall()
                times[q].append(time.perf_counter() - t)
        oracle = _partners_oracle(store_b38, tok, B38_K)
        assert {q: r == oracle for q, r in got.items()} == dict.fromkeys(times, True), name
        con = _fresh(store_b38)
        (loci, contigs), = con.execute("SELECT count(*), count(DISTINCT contig_uid) FROM token_slots(?)", [tok]).fetchall()
        out[name] = {"anchors": loci, "anchor_contigs": contigs,
                     **{f"{q}_median_s": round(statistics.median(ts), 4) for q, ts in times.items()},
                     **{f"{q}_s": [round(x, 4) for x in ts] for q, ts in times.items()}}
    record_property("b38_times", json.dumps(out))
    print("b38_times", json.dumps(out))  # shown with pytest -s
    hub = out["hub"]
    assert hub["partners_of_median_s"] <= B38_TIME_MARGIN * hub["tape_join_median_s"], json.dumps(out)


def test_build_time(store, record_property):
    record_property("build_s", round(store["build_s"], 1))
    assert store["build_s"] < 300, f"{store['build_s']:.0f} s for the {COPIES}x build and queries (card D-05: under 5 minutes)"


# --- card B-27: contig labels and cluster membership by their lookup keys ---------------------------

LABEL_SQL = "SELECT accession, original_header, n_occurrences FROM contig_label WHERE contig_uid = ?"
MEMBER_SQL = "SELECT level, cluster_uid, rep_uid FROM member_cluster WHERE member_uid = ? ORDER BY level"


@pytest.fixture(scope="module")
def lookups(store, tmp_path_factory):
    """The store with its release's `member_cluster`, written from the synthetic `cluster_member` by the
    builder's own function at the scaled caps, and a manifest and catalog that list it."""
    from biobox.l2 import cluster

    g, root = store["g"], store["root"]
    rel = g["REL"]
    mp = pytest.MonkeyPatch()
    mp.setattr(bio, "write_table", _scaled_write_table(bio.write_table))
    try:
        cluster.write_member_cluster(bio.part_files(root / f"l2/cluster_member/release={rel}"),
                                     root / f"l2/member_cluster/release={rel}", tmp=tmp_path_factory.mktemp("mc"))
    finally:
        mp.undo()
    man = catalog.build_manifest(root, rel, ["S0"])
    cat = catalog.generate_catalog(man, tmp_path_factory.mktemp("b27") / "catalog.duckdb")
    return {"root": root, "rel": rel, "V": g["V"], "man": man, "cat": cat}


def _table_name(root: Path, path: str) -> str:
    """The registry table of a file the profiler names (an absolute path under the store's root)."""
    return catalog.classify(Path(path).resolve().relative_to(Path(root).resolve()).as_posix())[0].name


def _reads_of(root: Path, scans: list[dict]) -> dict[str, int]:
    """Row groups read per registry table, summed over the profiled scans."""
    out: dict[str, int] = {}
    for s in scans:
        if s["files"]:
            t = _table_name(root, s["files"][0])
            out[t] = out.get(t, 0) + s["row_groups_scanned"]
    return out


def _query(root: Path, cat: Path, sql: str, params=None, tables: dict | None = None) -> tuple[list, dict[str, int], int]:
    """`sql` on a fresh `open_catalog` connection (cold): its rows, the row groups read per table, and the
    connection's `dynamic_or_filter_threshold` (`catalog.set_join_key_filter`)."""
    con = catalog.open_catalog(cat, threads=4)
    try:
        for name, tbl in (tables or {}).items():
            con.register(name, tbl)
        out = catalog.profile_query(con, sql, params)
        thr = int(con.execute("SELECT current_setting('dynamic_or_filter_threshold')").fetchall()[0][0])
    finally:
        con.close()
    return out["result"], _reads_of(root, out["scans"]), thr


def _manifest_files(root: Path, man: dict, table_name: str) -> list[Path]:
    return [root / f["path"] for f in man["files"] if f["table"] == table_name]


def _footer_ranges(files: list[Path], column: str) -> list[tuple]:
    """Oracle: (min, max) of `column` of every row group of `files`, in file then row-group order, from the
    Parquet footers."""
    out = []
    for p in files:
        md = pq.ParquetFile(p).metadata
        ci = next(i for i in range(md.num_columns) if md.schema.column(i).path == column)
        out += [(md.row_group(g).column(ci).statistics.min, md.row_group(g).column(ci).statistics.max)
                for g in range(md.num_row_groups)]
    return out


def _holding(files: list[Path], column: str, keys) -> int:
    """Oracle: the row groups whose footer range of `column` holds one of `keys`."""
    return sum(any(lo <= k <= hi for k in keys) for lo, hi in _footer_ranges(files, column))


def _spread(files: list[Path], column: str, keys) -> int:
    """Oracle: the row groups whose footer range of `column` meets the keys' range, min(keys) to max(keys)."""
    lo_k, hi_k = min(keys), max(keys)
    return sum(lo <= hi_k and lo_k <= hi for lo, hi in _footer_ranges(files, column))


def _assert_key_bound(read: int, keys: list, hold: int, spread: int, total: int, thr: int, levels: int = 1) -> None:
    """The bound of card B-27 for keys from a joined table: at most 2 row groups per key and level looked up (AGENTS
    §4.4; `member_cluster` is sorted by level, then member, so a member at both levels is two lookups) and no more
    than the keys' spread; with no more keys than the connection's join filter lists (`thr`), exactly the row
    groups that hold a key."""
    assert total >= MIN_RG, f"{total} row groups in the table"
    assert 1 <= read <= MAX_RG_PER_KEY * levels * len(keys) and read <= spread, (read, hold, spread, total)
    if len(keys) <= thr:
        assert read == hold, f"{read} row groups read; {hold} hold the {len(keys)} keys (spread {spread} of {total})"


def _occurrence_oracle(root: Path, man: dict) -> dict[bytes, tuple]:
    """Brute force over every `contig_occurrence` row: per contig, the accession, header, l0_object_sha256 and
    record_ordinal of its first record by (l0_object_sha256, record_ordinal), and its number of records."""
    first: dict[bytes, list] = {}
    cols = ["contig_uid", "l0_object_sha256", "record_ordinal", "accession", "original_header"]
    for f in _manifest_files(root, man, "contig_occurrence"):
        for r in pq.read_table(f, columns=cols, partitioning=None).to_pylist():
            k = (r["l0_object_sha256"], r["record_ordinal"])
            cur = first.get(r["contig_uid"])
            if cur is None:
                first[r["contig_uid"]] = [k, r["accession"], r["original_header"], 1]
                continue
            cur[3] += 1
            if k < cur[0]:
                cur[0], cur[1], cur[2] = k, r["accession"], r["original_header"]
    return {cu: (acc, hdr, n, k[0], k[1]) for cu, (k, acc, hdr, n) in first.items()}


def _membership_oracle(root: Path, man: dict) -> dict[tuple[bytes, str], list[tuple[bytes, bytes]]]:
    """Brute force over every `cluster_member` row of the release: (member, level) -> [(cluster, representative)]."""
    out: dict[tuple[bytes, str], list[tuple[bytes, bytes]]] = {}
    for f in _manifest_files(root, man, "cluster_member"):
        t = pq.read_table(f, columns=["member_uid", "level", "cluster_uid", "rep_uid"], partitioning=None)
        for m, lv, c, r in zip(*(t.column(i).to_pylist() for i in range(4))):
            out.setdefault((m, lv), []).append((c, r))
    return out


def _ends(ordered: list, n: int) -> list:
    """`n` keys from both ends of the key order, half the first of `ordered` and half the last, so the keys'
    range spans the table while the keys lie in few row groups."""
    return ordered[: n // 2] + ordered[-(n - n // 2):]


LABEL_PICKS = [("40 at both ends", 40), ("100 at both ends", 100), ("100 random", 100)]


@pytest.mark.covers("C-01")
def test_b27_contig_label_conforms_and_is_each_contigs_first_record(store):
    """Card B-27: `contig_label` (written by ingest) conforms to its spec, holds one row per contig sorted by
    contig_uid, and each row is the contig's first `contig_occurrence` record by (l0_object_sha256,
    record_ordinal) with its number of records (the brute-force oracle)."""
    files = _files(store, "contig_label")
    assert files, "ingest wrote no contig_label files"
    t = pa.concat_tables(pq.read_table(f, partitioning=None) for f in files)
    assert conformance_problems(t, table("contig_label")) == []
    uids = t.column("contig_uid").to_pylist()
    assert uids == sorted(uids) and len(set(uids)) == len(uids)
    got = {r["contig_uid"]: (r["accession"], r["original_header"], r["n_records"], r["l0_object_sha256"],
                             r["record_ordinal"]) for r in t.to_pylist()}
    assert got == _occurrence_oracle(store["root"], store["man"])


@pytest.mark.covers("C-01")
def test_b27_member_cluster_conforms_and_equals_cluster_member(lookups):
    """Card B-27: `member_cluster` conforms to its spec, is sorted by (level, member_uid), and every (member,
    level) of `cluster_member` appears in it exactly once, with the same cluster and representative."""
    root, man = lookups["root"], lookups["man"]
    files = _manifest_files(root, man, "member_cluster")
    assert files
    t = pa.concat_tables(pq.read_table(f, partitioning=None) for f in files)
    assert conformance_problems(t, table("member_cluster")) == []
    keys = list(zip(t.column("level").to_pylist(), t.column("member_uid").to_pylist()))
    assert keys == sorted(keys) and len(set(keys)) == len(keys)
    got: dict[tuple[bytes, str], list[tuple[bytes, bytes]]] = {}
    for r in t.to_pylist():
        got.setdefault((r["member_uid"], r["level"]), []).append((r["cluster_uid"], r["rep_uid"]))
    want = _membership_oracle(root, man)
    assert all(len(v) == 1 for v in want.values())
    assert got == want


@pytest.mark.covers("LY-06")
def test_b27_label_literal_key_reads_at_most_2_row_groups(store, record_property):
    """Card B-27: the catalog's `contig_label` by a literal contig_uid reads at most 2 row groups of
    `contig_label` (of at least MIN_RG) and nothing of `contig_occurrence`, and gives the oracle's label."""
    o = _occurrence_oracle(store["root"], store["man"])
    key = sorted(o)[len(o) // 2]
    rows, reads, _ = _query(store["root"], store["g"]["cat"], LABEL_SQL, [key])
    assert rows == [o[key][:3]]
    n = sum(f["row_groups"] for f in store["man"]["files"] if f["table"] == "contig_label")
    record_property("label_reads", json.dumps({"read": reads, "of": n}))
    assert set(reads) == {"contig_label"} and 1 <= reads["contig_label"] <= MAX_RG_PER_KEY, reads
    assert n >= MIN_RG, f"contig_label has {n} row groups at cap {CONTIG_CAP}"


@pytest.mark.covers("LY-06")
@pytest.mark.parametrize("pick,n", LABEL_PICKS, ids=[p[0] for p in LABEL_PICKS])
def test_b27_label_joined_keys_read_the_row_groups_that_hold_them(store, pick, n, record_property):
    """Card B-27: contig_uids from a joined table read nothing of `contig_occurrence`, and of `contig_label`
    the bound of `_assert_key_bound`. "at both ends" takes half the keys from the first row groups and half
    from the last, so their range spans the table: 40 such keys are fewer than DuckDB's join filter lists
    and read only the row groups that hold them, not the spread."""
    import random

    o = _occurrence_oracle(store["root"], store["man"])
    ordered = sorted(o)
    keys = _ends(ordered, n) if "ends" in pick else sorted(random.Random(SEED).sample(ordered, n))
    tbl = pa.table({"contig_uid": pa.array(keys, pa.binary(16))})
    rows, reads, thr = _query(store["root"], store["g"]["cat"],
                              "SELECT l.contig_uid, l.accession, l.original_header, l.n_occurrences FROM k k "
                              "JOIN contig_label l USING (contig_uid) ORDER BY l.contig_uid", tables={"k": tbl})
    assert rows == [(k, *o[k][:3]) for k in keys]
    files = _files(store, "contig_label")
    hold, spread, total = (_holding(files, "contig_uid", keys), _spread(files, "contig_uid", keys),
                           len(_footer_ranges(files, "contig_uid")))
    record_property("label_reads", json.dumps({"read": reads, "hold": hold, "spread": spread, "of": total, "thr": thr}))
    assert set(reads) == {"contig_label"}, reads
    _assert_key_bound(reads["contig_label"], keys, hold, spread, total, thr)
    if pick == "40 at both ends":
        assert hold < spread == total and len(keys) <= thr, (hold, spread, total, thr)


def _member_keys(files: list[Path], n: int) -> tuple[list[bytes], int]:
    """`n` member_uids, half from the second row group of `member_cluster` and half from the second to last (the
    row groups numbered over the files in order; a row group's first and last member, which the next row group
    may hold too, are left out), and how many row groups the table has."""
    groups = []
    for p in files:
        pf = pq.ParquetFile(p)
        groups += [pf.read_row_group(g, columns=["member_uid"]).column(0).to_pylist() for g in range(pf.metadata.num_row_groups)]
    assert len(groups) >= MIN_RG, f"member_cluster has {len(groups)} row groups at cap {ROW_CAP}"
    a, b = sorted(set(groups[1]))[1:-1], sorted(set(groups[-2]))[1:-1]
    return a[: n // 2] + b[-(n - n // 2):], len(groups)


@pytest.mark.covers("LY-06")
def test_b27_membership_literal_key_reads_at_most_2_row_groups_per_level(lookups, record_property):
    """Card B-27: the catalog's `member_cluster` (sorted by level, then member) by a literal member_uid reads at
    most 2 row groups at one level and at most 2 per level at both (of at least MIN_RG), nothing of
    `cluster_member`, and gives the member's NR50 and NR90 clusters as the oracle."""
    root, man = lookups["root"], lookups["man"]
    want = _membership_oracle(root, man)
    members = sorted({m for m, _ in want})
    key = members[len(members) // 3]
    rows, reads, _ = _query(root, lookups["cat"], MEMBER_SQL, [key])
    assert rows == [(lv, *want[(key, lv)][0]) for lv in ("nr50", "nr90")]
    one, reads1, _ = _query(root, lookups["cat"], "SELECT cluster_uid, rep_uid FROM member_cluster "
                                                  "WHERE level = 'nr90' AND member_uid = ?", [key])
    assert one == [want[(key, "nr90")][0]]
    _, n = _member_keys(_manifest_files(root, man, "member_cluster"), 2)
    record_property("member_reads", json.dumps({"both levels": reads, "nr90": reads1, "of": n}))
    assert set(reads) == {"member_cluster"} and 1 <= reads["member_cluster"] <= 2 * MAX_RG_PER_KEY, reads
    assert set(reads1) == {"member_cluster"} and 1 <= reads1["member_cluster"] <= MAX_RG_PER_KEY, reads1


MEMBER_PICKS = [("40 from two row groups", 40), ("100 from two row groups", 100), ("100 random", 100)]


@pytest.mark.covers("LY-06")
@pytest.mark.parametrize("pick,n", MEMBER_PICKS, ids=[p[0] for p in MEMBER_PICKS])
def test_b27_membership_joined_keys_read_the_row_groups_that_hold_them(lookups, pick, n, record_property):
    """Card B-27: member_uids from a joined table, looked up at both levels, read nothing of `cluster_member`,
    and of `member_cluster` the bound of `_assert_key_bound` (2 row groups per key and level). "from two row
    groups" takes them from the second and the second to last row group, so their range spans the table while
    few row groups hold them (their own two, those that hold the same members at the other level, and the one
    where the levels meet)."""
    import random

    root, man = lookups["root"], lookups["man"]
    files = _manifest_files(root, man, "member_cluster")
    want = _membership_oracle(root, man)
    keys = (_member_keys(files, n)[0] if "two" in pick
            else sorted(random.Random(SEED).sample(sorted({m for m, _ in want}), n)))
    rows, reads, thr = _query(root, lookups["cat"],
                              "SELECT m.member_uid, m.level, m.cluster_uid, m.rep_uid FROM k JOIN member_cluster m "
                              "USING (member_uid) ORDER BY m.member_uid, m.level",
                              tables={"k": pa.table({"member_uid": pa.array(keys, pa.binary(16))})})
    assert rows == [(m, lv, *want[(m, lv)][0]) for m in keys for lv in ("nr50", "nr90")]
    hold, spread, total = (_holding(files, "member_uid", keys), _spread(files, "member_uid", keys),
                           len(_footer_ranges(files, "member_uid")))
    record_property("member_reads", json.dumps({"read": reads, "hold": hold, "spread": spread, "of": total, "thr": thr}))
    assert set(reads) == {"member_cluster"}, reads
    _assert_key_bound(reads["member_cluster"], keys, hold, spread, total, thr, levels=2)
    if pick == "40 from two row groups":
        assert hold < spread and len(keys) <= thr, (hold, spread, total, thr)


class _Profiled:
    """A connection whose `execute` runs each statement through `catalog.profile_query` and adds up the row
    groups it read per registry table: the export functions take it in place of their connection."""

    def __init__(self, con, root: Path):
        self.con, self.root, self.reads = con, root, {}

    def execute(self, sql: str, params=None):
        for t, n in _reads_of(self.root, catalog.profile_query(self.con, sql, params)["scans"]).items():
            self.reads[t] = self.reads.get(t, 0) + n


@pytest.mark.covers("LY-06")
@pytest.mark.parametrize("n", [40, 100])
def test_b27_export_lookups_read_the_row_groups_that_hold_the_keys(lookups, n, record_property):
    """Card B-27: the hit and locus exports' label and cluster lookups (`export._annotate_contigs` and
    `_annotate_proteins` on `open_store`'s connection) for `n` contigs and `n` proteins whose ranges span the
    tables never read `contig_occurrence` or `cluster_member`, read `contig_label` and `member_cluster` within
    the bound of `_assert_key_bound`, and give the oracles' labels and clusters."""
    from biobox.l3 import tapes as tapes_mod
    from biobox.search import export

    root, man, rel, v = lookups["root"], lookups["man"], lookups["rel"], lookups["V"]
    o = _occurrence_oracle(root, man)
    keys = {"contig_label": ("contig_uid", _ends(sorted(o), n)),
            "member_cluster": ("member_uid", _member_keys(_manifest_files(root, man, "member_cluster"), n)[0])}
    want = _membership_oracle(root, man)
    with export.open_store(root, release=rel, gene_call_version=v, tape_version=tapes_mod.tape_version([v]), manifest=man,
                           threads=4, relations=("prot", "dom", "clus", "lab", "cidx", "dim")) as st:
        con = st.con
        thr = int(con.execute("SELECT current_setting('dynamic_or_filter_threshold')").fetchall()[0][0])
        con.register("_cin", pa.table({"contig_uid": pa.array(keys["contig_label"][1], pa.binary(16))}))
        con.execute("CREATE TEMP TABLE _c AS SELECT contig_uid FROM _cin")
        con.register("_pin", pa.table({"protein_uid": pa.array(keys["member_cluster"][1], pa.binary(16))}))
        con.execute("CREATE TEMP TABLE _p AS SELECT protein_uid FROM _pin")
        px = _Profiled(con, root)
        export._annotate_contigs(px, False)
        export._annotate_proteins(px, None)
        lab = con.execute("SELECT contig_uid, accession, header, n_records FROM _lab ORDER BY contig_uid").fetchall()
        cl = con.execute("SELECT protein_uid, nr90, nr50, nr50_rep FROM _cl ORDER BY protein_uid").fetchall()
    assert lab == [(k, *o[k][:3]) for k in keys["contig_label"][1]]
    assert cl == [(m, want[(m, "nr90")][0][0], want[(m, "nr50")][0][0], want[(m, "nr50")][0][1])
                  for m in keys["member_cluster"][1]]
    record_property("export_reads", json.dumps({"read": px.reads, "thr": thr}))
    assert not {"contig_occurrence", "cluster_member"} & set(px.reads), px.reads
    for t, (col, ks) in keys.items():
        files = _manifest_files(root, man, t)
        _assert_key_bound(px.reads.get(t, 0), ks, _holding(files, col, ks), _spread(files, col, ks),
                          len(_footer_ranges(files, col)), thr, levels=2 if t == "member_cluster" else 1)


def test_b27_contig_label_view_merges_the_release_stages(tmp_path):
    """Card B-27: `contig_label` is written per stage; the catalog's view and the exports' label of a release over
    two stages give, per contig, the first record by (l0_object_sha256, record_ordinal) over both stages and the
    records of both (the brute-force oracle over `contig_occurrence`). Stage S1 holds a second copy of S0's
    lambda genome under another name and source, and T4."""
    import yaml

    from biobox import fasta
    from biobox.l0seq import SourceObject, ingest
    from biobox.search import export

    fix = REPO / "tests/fixtures/micro_s0"
    objs = yaml.safe_load((fix / "sources.yaml").read_text())["objects"]
    lam_path, t4_path = fix / objs[0]["file"], fix / objs[1]["file"]
    dup = tmp_path / "lambda_again.fna"
    dup.write_bytes(b"\n".join(fasta.format_record("lambda_again", fasta.read(lam_path)[0].seq, 70)) + b"\n")
    root = tmp_path / "root"
    ingest([SourceObject(lam_path, "micro", "1")], root, "S0")
    ingest([SourceObject(dup, "other", "1"), SourceObject(t4_path, "micro", "1")], root, "S1")
    for s in ("S0", "S1"):
        assert bio.part_files(root / f"l0seq/contig_label/stage={s}"), f"no contig_label files in stage {s}"
    rel = "rS0.b27"
    catalog.define_release(root, rel, ["S0", "S1"])
    man = catalog.build_manifest(root, rel, ["S0", "S1"])
    o = _occurrence_oracle(root, man)
    assert sorted(v[2] for v in o.values()) == [1, 2]
    con = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / "c.duckdb"), threads=2)
    try:
        view = con.execute("SELECT contig_uid, accession, original_header, n_occurrences FROM contig_label "
                           "ORDER BY contig_uid").fetchall()
    finally:
        con.close()
    assert view == [(k, *o[k][:3]) for k in sorted(o)]
    with export.open_store(root, release=rel, gene_call_version="none", manifest=man, threads=2,
                           relations=("lab", "cidx", "dim")) as st:
        st.con.register("_cin", pa.table({"contig_uid": pa.array(sorted(o), pa.binary(16))}))
        st.con.execute("CREATE TEMP TABLE _c AS SELECT contig_uid FROM _cin")
        export._annotate_contigs(st.con, False)
        lab = st.con.execute("SELECT contig_uid, accession, header, n_records FROM _lab ORDER BY contig_uid").fetchall()
    assert lab == view
