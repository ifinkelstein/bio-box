"""Ring R2: `token_postings` on micro-S0 (card B-10; planning/tasks/TASKS.md §2).

Micro-S0 is built as in test_micro_s0_l3 (a 6-model domain library from lambda proteins for the
family tokens, NR50/NR90 clusters from MMseqs2 or the deterministic stand-in). The oracle is a
plain Python loop over the tape rows and the tape-token rows read with pyarrow: one posting per
(slot, level) with a token, every one at `seg_idx` 0 (per-domain family postings with
`seg_idx` > 0 are a later card). Neighbour questions over the catalog's `token_postings` view are
checked against `pairs_of_contig` (tests/oracles/pairs_oracle.py), which survives as the pair oracle
(AGENTS §3.5); the pair macros' pruning and the hub refusal are tested here on postings written for the
purpose (card B-11).
"""

from __future__ import annotations

import hashlib
import random
import shutil
import sys
from collections import Counter
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog, contracts, ids, io, provenance
from biobox.l3 import postings, tapes, tokens
from biobox.registry.schemas import conformance_problems, table
from oracles.pairs_oracle import pairs_of_contig

sys.path.insert(0, str(Path(__file__).parent))
import test_micro_s0_l3 as micro  # noqa: E402

pytestmark = pytest.mark.r2

REL, TV, V = micro.REL, micro.TV, micro.V
LEVELS = ("family", "nr50", "nr90")
LEVEL_BITS = {"family": 1, "nr50": 2, "nr90": 3}  # PLAN.md §4.1: token = level (3) | source (5) | local (56)
COLS = ("token_id", "contig_uid", "ordinal", "seg_idx", "feature_uid", "start", "end", "strand", "n_slots",
        "is_circular", "contig_len")
SORT = ("token_id", "contig_uid", "ordinal", "seg_idx")  # the registry's sort key, written out here
EST_ROW = 20  # the builder's per-row estimate (bytes), which cuts parts and picks the bucket width
POST = "l3/token_postings"
NR90_BASE = LEVEL_BITS["nr90"] << 61


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    d = tmp_path_factory.mktemp("postings")
    root = d / "a"
    circ = micro.build(root, with_postings=False)
    res = postings.run(root, REL, ["S0"])
    man = catalog.build_manifest(root, REL, ["S0"])
    con = catalog.open_catalog(catalog.generate_catalog(man, d / "c.duckdb"))
    yield root, man, con, circ, res
    con.close()


# --- oracle --------------------------------------------------------------------------------------


def _rows(root: Path, glob: str) -> list[dict]:
    return [r for f in sorted(root.glob(glob)) for r in pq.read_table(f, partitioning=None).to_pylist()]


def _oracle(root: Path, release: str = REL, tv: str = TV) -> list[dict]:
    """Every (slot, level) with a token, from the tapes and the tape tokens read straight from the files."""
    tps = {r["contig_uid"]: r for r in _rows(root, f"l3/tapes/tape_version={tv}/stage=*/part-*.parquet")}
    out = []
    for tok in _rows(root, f"l3/tape_tokens/tape_version={tv}/release={release}/part-*.parquet"):
        t = tps[tok["contig_uid"]]
        for i in range(t["n_slots"]):
            for lvl in LEVELS:
                tid = tok[f"{lvl}_tokens"][i]
                if tid is None:
                    continue
                out.append({"token_id": tid, "contig_uid": t["contig_uid"], "ordinal": i, "seg_idx": 0,
                            "feature_uid": t["feature_uids"][i], "start": t["start_bps"][i], "end": t["end_bps"][i],
                            "strand": t["strands"][i], "n_slots": t["n_slots"], "is_circular": t["is_circular"],
                            "contig_len": t["contig_len"]})
    return out


def _key(r: dict) -> tuple:
    return tuple(r[c] for c in COLS)


def _part_no(p: Path) -> int:
    return int(p.stem.removeprefix("part-"))


def _files(root: Path, release: str = REL, tv: str = TV) -> dict[str, list[Path]]:
    """Bucket -> its part files in part order."""
    base = root / f"{POST}/release={release}/tape_version={tv}"
    out = {}
    for d in sorted(base.iterdir()):
        assert d.name.startswith("token_bucket="), d
        out[d.name.removeprefix("token_bucket=")] = sorted(d.glob("part-*.parquet"), key=_part_no)
    return out


def _written(root: Path, release: str = REL, tv: str = TV) -> list[dict]:
    return [r for fs in _files(root, release, tv).values() for f in fs
            for r in pq.read_table(f, partitioning=None).to_pylist()]


def _level(tid: int) -> str:
    return {v: k for k, v in LEVEL_BITS.items()}[tid >> 61]


def _fresh(root: Path, tmp_path: Path, name: str, drop=(POST,)) -> Path:
    fresh = tmp_path / name
    shutil.copytree(root, fresh)
    for d in drop:
        shutil.rmtree(fresh / d, ignore_errors=True)
    return fresh


def _digests(root: Path) -> dict[str, str]:
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted((root / POST).rglob("*.parquet"))}


# --- content -------------------------------------------------------------------------------------


@pytest.mark.covers("LY-02")
def test_b10_postings_equal_the_oracle(built):
    root, _, _, _, res = built
    want, got = _oracle(root), _written(root)
    assert sorted(map(_key, got)) == sorted(map(_key, want))
    assert len(got) == len(set(map(_key, got)))
    # the row count is the number of non-null tokens over every slot and level of the tape tokens
    toks = _rows(root, f"l3/tape_tokens/tape_version={TV}/release={REL}/part-*.parquet")
    by_level = {lvl: sum(t is not None for r in toks for t in r[f"{lvl}_tokens"]) for lvl in LEVELS}
    assert res["rows"] == len(got) == sum(by_level.values())
    assert res["rows_by_level"] == by_level == dict(Counter(_level(r["token_id"]) for r in got))
    assert by_level["family"] > 0 and by_level["nr50"] > 0 and by_level["nr90"] > 0  # the synthetic library maps
    assert res["tape_version"] == TV
    files = _files(root)
    assert res["buckets"] == len(files) and res["files"] == sum(len(fs) for fs in files.values())
    # every file and every bucket in part order is sorted by the spec's sort key
    for fs in files.values():
        rows = [r for f in fs for r in pq.read_table(f, partitioning=None).to_pylist()]
        keys = [tuple(r[c] for c in SORT) for r in rows]
        assert keys == sorted(keys)
    assert all(r["seg_idx"] == 0 for r in got)


def test_b10_ordinals_name_the_features_of_the_tape(built):
    """The ordinal of a posting is its feature's index in the tape, and its coordinates are the feature's."""
    root, _, _, _, _ = built
    feats = {r["feature_uid"]: r for r in _rows(root, f"l1/feature/producer={V}/stage=S0/part-*.parquet")}
    tps = {r["contig_uid"]: r for r in _rows(root, f"l3/tapes/tape_version={TV}/stage=S0/part-*.parquet")}
    for r in _written(root):
        t, f = tps[r["contig_uid"]], feats[r["feature_uid"]]
        assert t["feature_uids"][r["ordinal"]] == r["feature_uid"]
        assert (f["contig_uid"], f["start"], f["end"], f["strand"]) == (r["contig_uid"], r["start"], r["end"], r["strand"])
        assert (r["n_slots"], r["is_circular"], r["contig_len"]) == (t["n_slots"], t["is_circular"], t["contig_len"])


def test_b10_without_a_domain_library_the_family_level_adds_nothing(built, tmp_path):
    root, _, _, _, _ = built
    fresh = _fresh(root, tmp_path, "s", (POST, "l1/domains", "l3/token_dict", "l3/tape_tokens"))
    tokens.run(fresh, REL, ["S0"])
    res = postings.run(fresh, REL, ["S0"])
    got = _written(fresh)
    assert sorted(map(_key, got)) == sorted(map(_key, _oracle(fresh)))
    assert res["rows_by_level"]["family"] == 0 and res["rows"] == len(got) > 0
    assert not [r for r in got if r["token_id"] >> 61 == LEVEL_BITS["family"]]
    # cluster tokens do not depend on the families: the same NR50/NR90 postings as with the library
    clusters = [r for r in _written(root) if r["token_id"] >> 61 != LEVEL_BITS["family"]]
    assert sorted(map(_key, got)) == sorted(map(_key, clusters))


def test_b10_a_non_cds_slot_has_no_postings_and_shifts_later_ordinals(built, tmp_path):
    """One ncRNA (second producer `test-ncrna`) on the circular contig: it has no token, so no posting.
    The ordinals of the slots after it grow by one; this is expected, because ordinals are valid in one
    tape version only, and `feature_uid` is the stable key (AGENTS §3.3)."""
    root, _, _, circ, _ = built
    fresh = _fresh(root, tmp_path, "s")
    mg = next(iter(circ))
    fu = micro._add_ncrna(fresh, mg)
    producers = ["test-ncrna", V]
    tv2 = tapes.run(fresh, "S0", circ, producers=producers)["tape_version"]
    rel = "rS0.b10nc"
    shutil.copytree(fresh / f"l2/cluster_member/release={REL}", fresh / f"l2/cluster_member/release={rel}")
    tokens.run(fresh, rel, ["S0"], producers=producers)
    res = postings.run(fresh, rel, ["S0"], producers=producers)
    assert res["tape_version"] == tv2 != TV
    got = _written(fresh, rel, tv2)
    assert sorted(map(_key, got)) == sorted(map(_key, _oracle(fresh, rel, tv2)))
    assert not [r for r in got if r["feature_uid"] == fu]
    i = _tapes_of(fresh, tv2)[mg]["feature_uids"].index(fu)
    # the token dictionary is the same as the first release's (one more slot without a protein)
    old = _written(root)
    assert len(got) == len(old)

    def shifted(r):
        if r["contig_uid"] != mg:
            return r
        return {**r, "ordinal": r["ordinal"] + (r["ordinal"] >= i), "n_slots": r["n_slots"] + 1}

    assert sorted(map(_key, got)) == sorted(_key(shifted(r)) for r in old)
    assert any(r["contig_uid"] == mg and r["ordinal"] > i for r in got)
    # keyed by feature_uid, the postings are unchanged apart from the ordinals and the slot count
    strip = lambda r: (r["token_id"], r["contig_uid"], r["feature_uid"], r["seg_idx"], r["start"], r["end"])  # noqa: E731
    assert sorted(map(strip, got)) == sorted(map(strip, old))


def _tapes_of(root: Path, tv: str) -> dict[bytes, dict]:
    return {r["contig_uid"]: r for r in _rows(root, f"l3/tapes/tape_version={tv}/stage=*/part-*.parquet")}


# --- layout --------------------------------------------------------------------------------------


@pytest.mark.covers("C-01", "LY-02", "LY-03", "LY-05", "C-03", "C-08")
def test_b10_files_conform_and_the_manifest_accepts_them(built):
    root, man, con, _, _ = built
    files = [f for f in man["files"] if f["table"] == "token_postings"]
    assert files and all(f["partitions"]["release"] == REL and f["partitions"]["tape_version"] == TV for f in files)
    assert [f["path"] for f in files] == [f"{POST}/release={REL}/tape_version={TV}/token_bucket=all/part-0.parquet"]
    for f in files:
        assert conformance_problems(pq.read_table(root / f["path"], partitioning=None), table("token_postings")) == []
    only = {**man, "files": files}
    assert contracts.check_sort_order(only) == [] and contracts.check_statistics(only) == []
    assert contracts.check_sort_order(man) == [] and contracts.check_statistics(man) == []
    assert contracts.check_foreign_keys(con) == []
    assert provenance.check_provenance(man) == []
    assert ("token_postings", ("contig_uid",), "contig_tapes", ("contig_uid",)) in contracts.FOREIGN_KEYS
    assert ("token_postings", ("token_id", "release"), "token_dict", ("token_id", "release")) in contracts.FOREIGN_KEYS
    # without its release partition C-08 objects
    bad = {**files[0], "partitions": {k: v for k, v in files[0]["partitions"].items() if k != "release"}}
    assert any("token_postings without its release partition" in p
               for p in provenance.check_provenance({**man, "files": [bad]}))
    n = con.execute("SELECT count(*) FROM token_postings WHERE release = ?", [REL]).fetchall()[0][0]
    assert n == len(_written(root))


@pytest.mark.covers("C-03")
def test_b10_foreign_keys_catch_an_unknown_token(built, tmp_path):
    root, _, _, _, _ = built
    fresh = _fresh(root, tmp_path, "s")
    spec = table("token_postings")
    row = {**_written(root)[0], "token_id": 7 << 61}  # level 7 is no token level: not in token_dict
    io.write_table([pa.Table.from_pylist([row], schema=spec.schema)], spec,
                   fresh / f"{POST}/release={REL}/tape_version={TV}/token_bucket=all/part-0.parquet")
    man = catalog.build_manifest(fresh, REL, ["S0"])
    con = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / "c.duckdb"))
    try:
        assert any(p.startswith("token_postings('token_id', 'release')") for p in contracts.check_foreign_keys(con))
    finally:
        con.close()


@pytest.mark.covers("LY-05", "LY-06")
@pytest.mark.parametrize("width", [1, 2])
def test_b10_rows_sit_in_the_bucket_of_their_low_byte(built, tmp_path, width):
    root, _, _, _, _ = built
    fresh = _fresh(root, tmp_path, "s")
    res = postings.run(fresh, REL, ["S0"], prefix_width=width)
    files = _files(fresh)
    assert len(files) > 1 and res["buckets"] == len(files)
    for b, fs in files.items():
        assert len(b) == width
        for f in fs:
            for tid in pq.read_table(f, columns=["token_id"]).column(0).to_pylist():
                assert f"{tid & 0xFF:02x}"[:width] == b, (f, tid)
    assert sorted(map(_key, _written(fresh))) == sorted(map(_key, _written(root)))
    man = catalog.build_manifest(fresh, REL, ["S0"])  # LY-05: the derived token_bucket partition is checked
    only = {**man, "files": [f for f in man["files"] if f["table"] == "token_postings"]}
    assert contracts.check_sort_order(only) == [] and contracts.check_statistics(only) == []
    # a lookup by token opens the one bucket that can hold it
    con = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / "c.duckdb"))
    try:
        tid = _written(fresh)[0]["token_id"]
        b = f"{tid & 0xFF:02x}"[:width]
        prof = catalog.profile_query(con, "SELECT sum(ordinal) FROM token_postings WHERE token_bucket = ? "
                                          "AND token_id = ?", [b, tid])
        (scan,) = [s for s in prof["scans"] if any(f"/{POST}/" in x for x in s["files"])]
        assert scan["files_read"] == 1
    finally:
        con.close()


@pytest.mark.covers("LY-04")
@pytest.mark.parametrize("width", [0, 1])
def test_b10_parts_are_cut_at_max_file_bytes(built, tmp_path, width):
    root, _, _, _, _ = built
    [single] = _files(root)["all"]
    want = pq.read_table(single, partitioning=None).to_pylist()
    per_part = 7
    fresh = _fresh(root, tmp_path, "s")
    res = postings.run(fresh, REL, ["S0"], prefix_width=width, max_file_bytes=EST_ROW * per_part + EST_ROW // 2)
    files = _files(fresh)
    assert res["files"] == sum(len(fs) for fs in files.values())
    assert any(len(fs) > 1 for fs in files.values())
    for b, fs in files.items():
        assert [p.name for p in fs] == [f"part-{k}.parquet" for k in range(len(fs))]
        parts = [pq.read_table(f, partitioning=None).to_pylist() for f in fs]
        assert all(len(p) == per_part for p in parts[:-1]) and 1 <= len(parts[-1]) <= per_part
        mine = [r for r in want if b == "all" or f"{r['token_id'] & 0xFF:02x}"[:width] == b]
        assert [r for p in parts for r in p] == mine, b
    man = catalog.build_manifest(fresh, REL, ["S0"])
    only = {**man, "files": [f for f in man["files"] if f["table"] == "token_postings"]}
    assert contracts.check_sort_order(only) == [] and contracts.check_statistics(only) == []


# --- determinism and refusals --------------------------------------------------------------------


@pytest.mark.covers("D-01", "D-04", "LY-08")
def test_b10_bytes_do_not_depend_on_the_run(built, tmp_path):
    root, _, _, _, _ = built
    want = _digests(root)
    assert want
    for name in ("again", "other"):
        fresh = _fresh(root, tmp_path, name)
        postings.run(fresh, REL, ["S0"], tmp=tmp_path)
        assert _digests(fresh) == want, name
    # a second run into the same root refuses to overwrite and leaves the files as they were
    with pytest.raises(FileExistsError, match="overwrite"):
        postings.run(tmp_path / "again", REL, ["S0"])
    assert _digests(tmp_path / "again") == want


def test_b10_refusals_and_scratch(built, tmp_path):
    root, _, _, _, _ = built
    fresh = _fresh(root, tmp_path, "s")
    with pytest.raises(ValueError, match="scratch"):
        postings.run(fresh, REL, ["S0"], tmp=fresh / "tmp")
    with pytest.raises(FileNotFoundError, match="tape_tokens"):
        postings.run(fresh, "no-such-release", ["S0"])
    assert not (fresh / POST).exists()
    parent = tmp_path / "scratch-parent"
    parent.mkdir()
    postings.run(fresh, REL, ["S0"], tmp=parent)
    assert list(parent.iterdir()) == []  # the private scratch directory is removed


def test_b10_token_lists_of_another_length_than_the_tape_are_an_error(built, tmp_path):
    root, _, _, _, _ = built
    fresh = _fresh(root, tmp_path, "s", (POST, "l3/tape_tokens"))
    toks = _rows(root, f"l3/tape_tokens/tape_version={TV}/release={REL}/part-*.parquet")
    victim = toks[1]
    victim["nr50_tokens"] = victim["nr50_tokens"][:-1]
    spec = table("tape_tokens")
    io.write_table([pa.Table.from_pylist(toks, schema=spec.schema)], spec,
                   fresh / f"l3/tape_tokens/tape_version={TV}/release={REL}/part-0.parquet")
    with pytest.raises(ValueError, match=victim["contig_uid"].hex()):
        postings.run(fresh, REL, ["S0"])
    assert not (fresh / POST).exists()


# --- sufficiency: a neighbour question over the catalog --------------------------------------------

NEAR = """
    SELECT count(DISTINCT (a.contig_uid, a.ordinal)) FROM token_postings a
    JOIN token_postings b ON b.contig_uid = a.contig_uid AND b.release = a.release AND b.tape_version = a.tape_version
    WHERE a.release = $rel AND a.tape_version = $tv AND a.token_id = $a AND b.token_id = $b
      AND (CASE WHEN a.is_circular
                THEN least(abs(a.ordinal::BIGINT - b.ordinal), a.n_slots - abs(a.ordinal::BIGINT - b.ordinal))
                ELSE abs(a.ordinal::BIGINT - b.ordinal) END) BETWEEN 1 AND $k
"""


def _pair_oracle(root: Path, a: int, b: int, k: int) -> int:
    """Loci (contig, anchor slot) of A with B within k slots, from `pairs_of_contig` on every tape."""
    tps = _tapes_of(root, TV)
    loci = set()
    for tok in _rows(root, f"l3/tape_tokens/tape_version={TV}/release={REL}/part-*.parquet"):
        for p in pairs_of_contig(tps[tok["contig_uid"]], {lvl: tok[f"{lvl}_tokens"] for lvl in LEVELS}, k):
            if p["anchor_token"] == a and p["partner_token"] == b:
                loci.add((p["contig_uid"], p["anchor_ordinal"]))
    return len(loci)


@pytest.mark.covers("Q1")
def test_b10_a_self_join_answers_loci_of_a_near_b(built):
    root, _, con, circ, _ = built
    tps = _tapes_of(root, TV)
    toks = {r["contig_uid"]: r for r in _rows(root, f"l3/tape_tokens/tape_version={TV}/release={REL}/part-*.parquet")}
    mg = next(iter(circ))
    lin = next(cu for cu in sorted(tps) if not tps[cu]["is_circular"])
    cases = []
    # circular: the first slot and the last, one slot apart across the origin and n - 1 apart on a line
    col = toks[mg]["nr90_tokens"]
    n = tps[mg]["n_slots"]
    assert n > 20 and col[0] is not None and col[n - 1] is not None
    cases.append(("circular", col[0], col[n - 1]))
    # linear: two slots three apart in the middle of the contig
    col = toks[lin]["nr50_tokens"]
    i = next(i for i in range(len(col) // 2, len(col) - 3) if col[i] is not None and col[i + 3] is not None)
    cases.append(("linear", col[i], col[i + 3]))
    checked = 0
    for what, a, b in cases:
        for k in (1, 3, 10):
            want = _pair_oracle(root, a, b, k)
            got = con.execute(NEAR, {"rel": REL, "tv": TV, "a": a, "b": b, "k": k}).fetchall()[0][0]
            assert got == want, (what, k)
            checked += want > 0
        assert _pair_oracle(root, a, b, 3) >= 1, what
    assert checked >= 4
    # the circular case needs the distance round the origin: on a line, slots 0 and n - 1 are not neighbours
    _, a, b = cases[0]
    linear = NEAR.replace("CASE WHEN a.is_circular", "CASE WHEN false")
    assert con.execute(linear, {"rel": REL, "tv": TV, "a": a, "b": b, "k": 1}).fetchall()[0][0] < _pair_oracle(root, a, b, 1)


# --- pruning at several row groups -----------------------------------------------------------------

N_CONTIGS, N_SLOTS, N_TOKENS = 260, 1000, 997  # 260,000 postings: three row groups of 122,880 rows


def _synthetic(root: Path) -> tuple[str, list[dict], list[dict]]:
    """Tapes and NR90-only tape tokens written directly (no gene calls), enough for several row groups."""
    rng = random.Random(10)
    tv = tapes.tape_version(["synthetic-b10"])
    tape_rows, tok_rows = [], []
    for c in range(N_CONTIGS):
        cu = ids.contig_uid(b"ACGT" * 25 + b"%d" % c)  # a uid only: no sequence is stored
        slots = [(100 * j, 100 * j + 90, 1 if rng.random() < 0.5 else -1) for j in range(N_SLOTS)]
        circular = c % 7 == 0
        tape_rows.append({"contig_uid": cu, "contig_len": 100 * N_SLOTS, "is_circular": circular, "n_slots": N_SLOTS,
                          **tapes.build_tape(slots, 100 * N_SLOTS, circular),
                          "protein_uids": [None] * N_SLOTS,
                          "feature_uids": [ids.feature_uid(cu, s, e, st, ids.KIND["cds"]) for s, e, st in slots]})
        nr90 = [None if rng.random() < 0.05 else 3 << 61 | rng.randrange(N_TOKENS) for _ in range(N_SLOTS)]
        tok_rows.append({"contig_uid": cu, "family_tokens": [None] * N_SLOTS, "nr50_tokens": [None] * N_SLOTS,
                         "nr90_tokens": nr90})
    tape_rows.sort(key=lambda r: r["contig_uid"])
    tok_rows.sort(key=lambda r: r["contig_uid"])
    spec = table("contig_tapes")
    io.write_table([pa.Table.from_pylist(tape_rows, schema=spec.schema)], spec,
                   root / f"l3/tapes/tape_version={tv}/stage=S0/part-0.parquet", max_slots=131_072,
                   slots_column="n_slots")
    spec = table("tape_tokens")
    io.write_table([pa.Table.from_pylist(tok_rows, schema=spec.schema)], spec,
                   root / f"l3/tape_tokens/tape_version={tv}/release={REL}/part-0.parquet", max_slots=131_072,
                   slots_column="nr90_tokens")
    return tv, tape_rows, tok_rows


@pytest.mark.covers("LY-06", "LY-02", "LY-03")
def test_b10_a_token_lookup_scans_at_most_two_row_groups(tmp_path):
    root = tmp_path / "r"
    tv, tape_rows, tok_rows = _synthetic(root)
    res = postings.run(root, REL, ["S0"], producers=["synthetic-b10"])
    n_tok = Counter(t for r in tok_rows for t in r["nr90_tokens"] if t is not None)
    assert res["rows"] == sum(n_tok.values()) > 2 * io.ROW_GROUP_ROWS
    [path] = _files(root, REL, tv)["all"]
    md = pq.ParquetFile(path).metadata
    assert md.num_row_groups >= 3
    assert [md.row_group(i).num_rows for i in range(md.num_row_groups - 1)] == [io.ROW_GROUP_ROWS] * (md.num_row_groups - 1)
    # the oracle at this size: posting rows equal the loop over tapes and tokens
    want = sorted(_key(r) for r in _oracle(root, REL, tv))
    assert sorted(map(_key, _written(root, REL, tv))) == want
    man = catalog.build_manifest(root, REL, ["S0"])
    assert contracts.check_sort_order(man) == [] and contracts.check_statistics(man) == []
    con = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / "c.duckdb"))
    try:
        # the tokens at the row-group boundaries, the most frequent one, and a sample
        first = [pq.ParquetFile(path).read_row_group(i, columns=["token_id"]).column(0) for i in range(md.num_row_groups)]
        probe = {g[0].as_py() for g in first} | {g[-1].as_py() for g in first} | {n_tok.most_common(1)[0][0]}
        probe |= set(random.Random(4).sample(sorted(n_tok), 10))
        for tid in sorted(probe):
            prof = catalog.profile_query(con, "SELECT count(*), sum(ordinal) FROM token_postings WHERE token_id = ?", [tid])
            (scan,) = [s for s in prof["scans"] if any(f"/{POST}/" in x for x in s["files"])]
            assert prof["result"][0][0] == n_tok[tid]
            assert scan["row_groups_scanned"] <= 2 and scan["row_groups_total"] == md.num_row_groups, tid
    finally:
        con.close()


# --- card B-11: the pair macros prune; partners_of refuses a hub ------------------------------------


def _near(tok_rows: list[dict], tape_rows: list[dict], a: int, b: int, k: int) -> tuple[int, int]:
    """(loci of A with B within k slots, (A slot, B slot) pairs), by a loop over the token lists of each contig."""
    circ = {t["contig_uid"]: t["is_circular"] for t in tape_rows}
    loci = pairs = 0
    for r in tok_rows:
        col, n = r["nr90_tokens"], len(r["nr90_tokens"])
        for i in (i for i, t in enumerate(col) if t == a):
            near = [j for j, t in enumerate(col) if t == b and j != i and
                    (min(abs(i - j), n - abs(i - j)) if circ[r["contig_uid"]] else abs(i - j)) <= k]
            loci += bool(near)
            pairs += len(near)
    return loci, pairs


@pytest.mark.covers("LY-06", "Q1", "Q4")
def test_b11_pair_macros_scan_at_most_two_row_groups_per_side(tmp_path):
    """q1_count and pairs_of for non-hub token pairs on postings of several row groups: each side of the self-join
    reads at most 2 row groups (AGENTS §4.4), including tokens whose rows straddle a row-group boundary."""
    root = tmp_path / "r"
    tv, tape_rows, tok_rows = _synthetic(root)
    postings.run(root, REL, ["S0"], producers=["synthetic-b10"])
    [path] = _files(root, REL, tv)["all"]
    md = pq.ParquetFile(path).metadata
    assert md.num_row_groups >= 3
    con = catalog.open_catalog(catalog.generate_catalog(catalog.build_manifest(root, REL, ["S0"]), tmp_path / "c.duckdb"))
    try:
        first = [pq.ParquetFile(path).read_row_group(i, columns=["token_id"]).column(0) for i in range(md.num_row_groups)]
        anchors = sorted({g[0].as_py() for g in first[1:]} | {g[-1].as_py() for g in first[:-1]}
                         | set(random.Random(5).sample(sorted({t for r in tok_rows for t in r["nr90_tokens"] if t}), 6)))
        checked = 0
        for a in anchors:
            # the partner: the token most often within 3 slots of A
            near = Counter(r["nr90_tokens"][j] for r in tok_rows for i, t in enumerate(r["nr90_tokens"]) if t == a
                           for j in range(max(0, i - 3), min(N_SLOTS, i + 4)) if j != i and r["nr90_tokens"][j])
            b = near.most_common(1)[0][0]
            loci, n_pairs = _near(tok_rows, tape_rows, a, b, 10)
            for sql, want in (("SELECT q1_count(?, ?, 10, NULL)", (loci,)),
                              ("SELECT count(*), sum(delta_bp) IS NOT NULL FROM pairs_of(?, ?, 10)", (n_pairs, True))):
                prof = catalog.profile_query(con, sql, [a, b])
                assert prof["result"][0] == want, (sql, a, b)
                scans = [s for s in prof["scans"] if any(f"/{POST}/" in x for x in s["files"])]
                assert len(scans) == 2, (sql, scans)  # one per side
                assert all(s["row_groups_scanned"] <= 2 and s["row_groups_total"] == md.num_row_groups
                           for s in scans), (sql, a, b, scans)
                checked += 1
        assert checked == 2 * len(anchors)
        straddle = {first[i][-1].as_py() for i in range(md.num_row_groups - 1)} & {g[0].as_py() for g in first[1:]}
        assert straddle and straddle <= set(anchors)  # tokens whose rows span two row groups were looked up
    finally:
        con.close()


@pytest.mark.covers("Q3")
def test_b11_partners_of_refuses_a_hub_token(tmp_path):
    """More than 100,000 posting rows: partners_of refuses and names pair_stats; at exactly 100,000 it answers."""
    root = tmp_path / "r"
    tv = tapes.tape_version(["synthetic-hub"])
    hub, big, partner, n = NR90_BASE | 1, NR90_BASE | 2, NR90_BASE | 3, 100_000
    cu = lambda tag, i: tag + i.to_bytes(15, "big")  # noqa: E731
    rows = [(hub, cu(b"h", i), 0) for i in range(n + 1)] + [(big, cu(b"d", i), 0) for i in range(n)] + \
        [(partner, cu(b"d", i), 1) for i in range(n)]
    tok, cus, ords = zip(*rows)
    m = len(rows)
    spec = table("token_postings")
    tbl = pa.table({"token_id": pa.array(tok, pa.uint64()), "contig_uid": pa.array(cus, pa.binary(16)),
                    "ordinal": pa.array(ords, pa.uint32()), "seg_idx": pa.array([0] * m, pa.uint16()),
                    "feature_uid": pa.array([c[:15] + bytes([o]) for c, o in zip(cus, ords)], pa.binary(16)),
                    "start": pa.array([100 * o for o in ords], pa.int64()),
                    "end": pa.array([100 * o + 90 for o in ords], pa.int64()), "strand": pa.array([1] * m, pa.int8()),
                    "n_slots": pa.array([2] * m, pa.uint32()), "is_circular": pa.array([False] * m),
                    "contig_len": pa.array([200] * m, pa.int64())}).cast(spec.schema)
    io.write_table([tbl], spec, root / f"{POST}/release={REL}/tape_version={tv}/token_bucket=all/part-0.parquet")
    con = catalog.open_catalog(catalog.generate_catalog(catalog.build_manifest(root, REL, []), tmp_path / "c.duckdb"))
    try:
        with pytest.raises(duckdb.Error, match=r"partners_of: token has more than 100000 loci; use pair_stats \(hub routing\)"):
            con.execute("SELECT * FROM partners_of(?, 10)", [hub]).fetchall()
        assert con.execute("SELECT * FROM partners_of(?, 1)", [big]).fetchall() == [(partner, n, n)]
        assert con.execute("SELECT q1_count(?, ?, 1, NULL)", [big, partner]).fetchall() == [(n,)]
    finally:
        con.close()
