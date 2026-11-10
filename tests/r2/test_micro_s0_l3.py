"""Ring R2: L3 tapes, tokens and postings on micro-S0, and the Q1-Q4, Q7, Q8 oracles (data_layer_and_tests.md §4.3).

The build: ingest, gene calls, a 6-model domain library from lambda proteins (family
tokens) and NR50/NR90 clusters from MMseqs2 where it is on PATH, else a deterministic
stand-in grouping (the oracles read the same tokens, so what they test does not depend on
how the clusters were made), tapes, tape tokens and token postings. Each production path
(catalog macros) is compared with an oracle: the pair macros over `token_postings` (card
B-11) with `pairs_of_contig` (tests/oracles/pairs_oracle.py) run on the tape rows read
straight from the Parquet shards, the partner counts with a brute-force loop over the
oracle's windows, and the tapes with a rebuild of the tape lists from protein_occurrence.
Covers Q1, Q2, Q3 (partner counts; z_score ranking needs pair_stats), Q4, Q7, Q8, C-04,
C-07, K-07, LY-06 for pair lookups, and D-01 for the L3 tables.
"""

from __future__ import annotations

import hashlib
import inspect
import random
import shutil
from collections import Counter, defaultdict
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pyhmmer
import pytest
import yaml

from hypothesis import HealthCheck, given, settings, strategies as st

from biobox import catalog, contracts, dim, fasta, ids, io, layout
from biobox.l0seq import SourceObject, ingest, revcomp
from biobox.l0seq.inputs import write_genome_input
from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION as V, GenomeInput, run as call_genes
from biobox.l3 import geometry, postings, tapes, tokens
from biobox.registry.schemas import conformance_problems, table
from oracles.pairs_oracle import FLAG_CENSORED, FLAG_CROSSES_ORIGIN, K_MAX, pairs_of_contig, window

pytestmark = pytest.mark.r2

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
MMSEQS = shutil.which("mmseqs")
SAMTOOLS = shutil.which("samtools")
REL = "rS0.l3"
LEVELS = ("family", "nr50", "nr90")
# the tape version of the single gene-call producer, computed here from TASKS.md §2 (not by biobox)
TV = "t-" + hashlib.sha256(V.encode("utf-8")).hexdigest()[:12]


def _library(root: Path, path: Path) -> None:
    occ = pq.read_table(root / f"l1/protein_occurrence/gene_call_version={V}/stage=S0/part-0.parquet")
    lam = ids.contig_uid(fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq)
    prots = domains.load_proteins(root, V, "S0")
    uids = sorted({r["protein_uid"] for r in occ.to_pylist() if r["contig_uid"] == lam}, key=lambda u: (-len(prots[u]), u))[:6]
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder, bg = pyhmmer.plan7.Builder(alphabet, seed=42), pyhmmer.plan7.Background(alphabet)
    with open(path, "wb") as fh:
        for i, u in enumerate(uids):
            seq = pyhmmer.easel.TextSequence(name=b"s", sequence=prots[u].replace("*", "")).digitize(alphabet)
            hmm, _, _ = builder.build(seq, bg)
            hmm.name, hmm.accession = f"lam{i}".encode(), f"TEST{i:04d}".encode()
            hmm.cutoffs.gathering = (25.0, 25.0)
            hmm.write(fh)


def _stand_in_clusters(root: Path, release: str) -> None:
    """NR90/NR50 stand-ins when MMseqs2 is absent: group by the first UID byte (97 and 23 groups)."""
    rows = []
    for p in sorted(domains.load_proteins(root, V, "S0")):
        for lvl, g in (("nr90", 97), ("nr50", 23)):
            rows.append({"member_uid": p, "level": lvl, "cluster_uid": bytes([p[0] % g, len(lvl)]) * 8,
                         "rep_uid": p, "cluster_params_hash": "stand-in"})
    rows.sort(key=lambda r: (r["level"], r["cluster_uid"], r["member_uid"]))
    d = root / f"l2/cluster_member/release={release}"
    d.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=table("cluster_member").schema), d / "part-0.parquet")


def build(root: Path, order: int = 1, with_postings: bool = True) -> set[bytes]:
    objs = CONTRACT["objects"][::order]
    ingest([SourceObject(FIX / o["file"], CONTRACT["source_name"], CONTRACT["source_release"], CONTRACT["license"], True)
            for o in objs], root, "S0")
    genomes, circ = [], set()
    for o in objs:
        cu = ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)
        c = (cu,) if o["topology"] == "circular" else ()
        circ |= set(c)
        genomes.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"], circular=c))
    call_genes(root, "S0", genomes)
    lib = root.parent / f"{root.name}.hmm"
    _library(root, lib)
    domains.run(root, "S0", V, "test", "1", lib, "gathering")
    if MMSEQS:
        from biobox.l2 import cluster
        tmp = root.parent / f"{root.name}.mmseqs"
        tmp.mkdir(parents=True, exist_ok=True)
        cluster.run(root, REL, V, ["S0"], tmp=tmp)
    else:
        _stand_in_clusters(root, REL)
    tapes.run(root, "S0", circ)
    tokens.run(root, REL, ["S0"])
    if with_postings:
        postings.run(root, REL, ["S0"])
    return circ


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    d = tmp_path_factory.mktemp("l3")
    circ = build(d / "a")
    man = catalog.build_manifest(d / "a", REL, ["S0"])
    con = catalog.open_catalog(catalog.generate_catalog(man, d / "c.duckdb"))
    yield d / "a", man, con, circ
    con.close()


# --- oracles -------------------------------------------------------------------------------


def _read(root: Path, glob: str) -> list[dict]:
    return [r for f in sorted(root.glob(glob)) for r in pq.read_table(f).to_pylist()]


def _oracle_rows(root: Path, release: str = REL, tv: str = TV) -> list[tuple[dict, dict]]:
    """(tape, tokens) per contig, straight from the shards."""
    t = {r["contig_uid"]: r for r in _read(root, f"l3/tapes/tape_version={tv}/stage=*/part-*.parquet")}
    k = {r["contig_uid"]: r for r in _read(root, f"l3/tape_tokens/tape_version={tv}/release={release}/part-*.parquet")}
    return [(t[cu], k[cu]) for cu in sorted(t)]


PAIR_COLS = ("contig_uid", "anchor_feature_uid", "partner_feature_uid", "anchor_ordinal", "partner_ordinal",
             "orf_offset", "junction_topology", "delta_bp", "flags")
PAIRS_SQL = f"SELECT {', '.join(PAIR_COLS)} FROM pairs_of(?, ?, ?)"


def _oracle_pairs(rows, k: int = K_MAX) -> dict[tuple[int, int], list[tuple]]:
    """`pairs_of_contig` on every tape, as rows of PAIR_COLS keyed by (anchor token, partner token)."""
    out: dict[tuple[int, int], list[tuple]] = defaultdict(list)
    for tape, tok in rows:
        fu = tape["feature_uids"]
        for p in pairs_of_contig(tape, {lvl: tok[f"{lvl}_tokens"] for lvl in LEVELS}, k):
            i, j = p["anchor_ordinal"], p["partner_ordinal"]
            out[(p["anchor_token"], p["partner_token"])].append(
                (p["contig_uid"], fu[i], fu[j], i, j, p["orf_offset"], p["junction_topology"], p["delta_bp"], p["flags"]))
    return out


def _loci(rows) -> int:
    """Distinct (contig, anchor slot) among oracle rows."""
    return len({(r[0], r[3]) for r in rows})


def _sample_tokens(con, n: int, seed: int) -> list[int]:
    """Stratified by document frequency decile, so the most frequent tokens are always sampled."""
    rows = con.execute("SELECT token_id, n_slots FROM token_dict WHERE release = ? ORDER BY n_slots, token_id",
                       [REL]).fetchall()
    rng = random.Random(seed)
    deciles = [rows[i * len(rows) // 10:(i + 1) * len(rows) // 10] for i in range(10)]
    picked = [rng.choice(d)[0] for d in deciles if d] + [rows[-1][0]]
    return sorted(set(picked))[:n]


# --- Q1-Q4 ---------------------------------------------------------------------------------


_PAIR_LIST = "(SELECT unnest($a)::UBIGINT AS a, unnest($b)::UBIGINT AS b) l"
# the macros for many token pairs in one statement: the same macro text, its arguments taken from a joined list
PAIRS_BULK = f"SELECT l.a, l.b, {', '.join('p.' + c for c in PAIR_COLS)} FROM {_PAIR_LIST}, " \
             "LATERAL (SELECT * FROM pairs_of(l.a, l.b, $k)) p"
Q1_BULK = f"SELECT l.a, l.b, q1_count(l.a, l.b, 1, NULL), q1_count(l.a, l.b, 3, NULL), " \
          f"q1_count(l.a, l.b, {K_MAX}, NULL) FROM {_PAIR_LIST}"


def _check_pairs_against_the_oracle(con, rows, pairs, literal: int = 40) -> dict[str, int]:
    """pairs_of (k = 10 and k = 3) and q1_count (k = 1, 3, 10) equal the oracle for every ordered token pair in
    `pairs`, column for column; `literal` of them also one statement at a time with the tokens as parameters,
    which also checks the order of pairs_of. Returns what the oracle rows exercised."""
    by_k = {k: _oracle_pairs(rows, k) for k in (1, 3, K_MAX)}
    full = by_k[K_MAX]
    lists = {"a": [a for a, _ in pairs], "b": [b for _, b in pairs]}
    for k in (K_MAX, 3):
        got = defaultdict(list)
        for a, b, *r in con.execute(PAIRS_BULK, {**lists, "k": k}).fetchall():
            got[(a, b)].append(tuple(r))
        assert set(got) <= set(pairs)
        for p in pairs:
            # with k = 3 the rows within 3 slots are unchanged: the censoring flag is about the +-K_MAX window
            assert sorted(got.get(p, [])) == sorted(r for r in full.get(p, []) if abs(r[5]) <= k), (p, k)
    for a, b, *counts in con.execute(Q1_BULK, lists).fetchall():
        assert counts == [_loci(by_k[k].get((a, b), [])) for k in (1, 3, K_MAX)], (a, b)
    for a, b in random.Random(15).sample(pairs, min(literal, len(pairs))):
        got = con.execute(PAIRS_SQL, [a, b, K_MAX]).fetchall()
        assert got == sorted(got, key=lambda r: (r[0], r[3], r[4])), (a, b)  # ordered by contig, anchor, partner
        assert sorted(got) == sorted(full.get((a, b), [])), (a, b)
        assert con.execute("SELECT q1_count(?, ?, ?, NULL)", [a, b, K_MAX]).fetchall()[0][0] == \
            _loci(full.get((a, b), [])), (a, b)
    seen = Counter()
    for p in pairs:
        for r in full.get(p, []):
            seen["rows"] += 1
            seen["crosses_origin"] += bool(r[8] & FLAG_CROSSES_ORIGIN)
            seen["censored"] += bool(r[8] & FLAG_CENSORED)
            seen[f"topology_{r[6]}"] += 1
    return seen


@pytest.mark.covers("Q1", "Q4")
def test_b11_pairs_of_and_q1_count_equal_the_pair_oracle(built):
    """Card B-11: every ordered token pair that co-occurs within 10 slots at a level of the fixture (two linear
    contigs and the odd circle of M. genitalium), column for column against `pairs_of_contig`, and q1_count for
    k in 1, 3, 10."""
    root, _, con, circ = built
    rows = _oracle_rows(root)
    assert sum(t["is_circular"] for t, _ in rows) == 1 and len(rows) == 3
    full = _oracle_pairs(rows)
    per_level = Counter(a >> 61 for a, _ in full)
    assert all(per_level[bits] > 0 for bits in (1, 2, 3)), per_level  # family, nr50, nr90
    seen = _check_pairs_against_the_oracle(con, rows, sorted(full), literal=60)
    assert seen["rows"] == sum(len(v) for v in full.values())
    assert seen["crosses_origin"] > 0 and seen["censored"] > 0, seen
    assert seen["topology_1"] and seen["topology_2"] and seen["topology_3"], seen
    # token pairs of one level that never meet within 10 slots: no rows, zero loci
    by_level = defaultdict(set)
    for _, k in rows:
        for lvl in LEVELS:
            by_level[lvl] |= {t for t in k[f"{lvl}_tokens"] if t is not None}
    rng = random.Random(12)
    apart = sorted({p for lvl in LEVELS for p in (tuple(rng.sample(sorted(by_level[lvl]), 2)) for _ in range(200))
                    if p not in full})
    assert len(apart) > 50
    assert _check_pairs_against_the_oracle(con, rows, apart)["rows"] == 0
    assert con.execute(PAIRS_SQL, [sorted(by_level["nr90"])[0], 7 << 61, K_MAX]).fetchall() == []  # an unknown token


@pytest.mark.covers("Q3")
def test_q3_partners_of_equals_a_brute_force_count(built):
    """partners_of: every token at every level within k slots of the loci of A, counted by a loop over the oracle's
    windows; ranked by loci."""
    root, _, con, _ = built
    rows = _oracle_rows(root)
    other_level = 0
    for a in _sample_tokens(con, 10, 3):
        for k in (1, K_MAX):
            loci: dict[int, set] = defaultdict(set)
            n_pairs: Counter = Counter()
            for tape, tok in rows:
                n, circ = tape["n_slots"], tape["is_circular"]
                cols = [tok[f"{lvl}_tokens"] for lvl in LEVELS]
                for i in range(n):
                    if a not in [c[i] for c in cols]:
                        continue
                    for j, _ in window(i, n, circ, k):
                        for c in cols:
                            if c[j] is not None:
                                loci[c[j]].add((tape["contig_uid"], i))
                                n_pairs[c[j]] += 1
            got = con.execute("SELECT partner_token, n_loci, n_pairs FROM partners_of(?, ?)", [a, k]).fetchall()
            assert sorted(got) == sorted((b, len(s), n_pairs[b]) for b, s in loci.items()), (a, k)
            assert got == sorted(got, key=lambda g: (-g[1], g[0]))  # ranked by loci, then token
            other_level += sum(b >> 61 != a >> 61 for b, _, _ in got)
    assert other_level > 0  # partners at the other token levels are counted too


@pytest.mark.covers("Q2", "C-05")
def test_q2_fetch_equals_shards_and_l1_rebuild(built):
    root, _, con, circ = built
    rows = {t["contig_uid"]: (t, k) for t, k in _oracle_rows(root)}
    occ = defaultdict(list)
    for r in _read(root, f"l1/protein_occurrence/gene_call_version={V}/stage=S0/part-0.parquet"):
        occ[r["contig_uid"]].append(r)
    cuids = sorted(rows)
    got = con.sql("SELECT * FROM tapes_of(?)", params=[cuids]).to_arrow_table().to_pylist()
    assert [g["contig_uid"] for g in got] == cuids
    for g in got:
        tape, tok = rows[g["contig_uid"]]
        for col, v in tape.items():
            assert g[col] == v, col
        for lvl in LEVELS:
            assert g[f"{lvl}_tokens"] == tok[f"{lvl}_tokens"]
        # independent rebuild from L1: slots in coordinate order (start, end, strand), feature and protein
        # UIDs, strands, coordinates
        slots = sorted(occ[g["contig_uid"]], key=lambda r: (r["start"], r["end"], r["strand"]))
        assert g["n_slots"] == len(slots)
        assert g["feature_uids"] == [r["feature_uid"] for r in slots]
        assert g["protein_uids"] == [r["protein_uid"] for r in slots]
        assert g["start_bps"] == [r["start"] for r in slots] and g["end_bps"] == [r["end"] for r in slots]
        assert g["strands"] == [r["strand"] for r in slots]
        assert g["is_circular"] == (g["contig_uid"] in circ)
    assert con.sql("SELECT count(*) FROM tapes_of(?)", params=[[b"\0" * 16]]).fetchall()[0][0] == 0


@pytest.mark.covers("Q7")
def test_q7_slice_equals_the_decompressed_source(built, tmp_path):
    root, _, con, _ = built
    rng = random.Random(7)
    for o in CONTRACT["objects"]:
        seq = fasta.read(FIX / o["file"])[0].seq
        cu = ids.contig_uid(seq)
        plain = tmp_path / f"{o['accession']}.fa"
        plain.write_bytes(b">" + o["accession"].encode() + b"\n" + seq + b"\n")
        for _ in range(30):
            a = rng.randrange(0, len(seq) - 1)
            b = min(len(seq), a + rng.randrange(1, 100_000))
            for strand in (1, -1):
                got = con.execute("SELECT seq_get(?, ?, ?, ?)", [cu, a, b, strand]).fetchall()[0][0].encode()
                assert got == (seq[a:b] if strand == 1 else revcomp(seq[a:b]))
                if SAMTOOLS and strand == 1 and _ < 5:
                    import subprocess
                    res = subprocess.run([SAMTOOLS, "faidx", str(plain), f"{o['accession']}:{a + 1}-{b}"],
                                         check=True, capture_output=True)
                    assert b"".join(res.stdout.splitlines()[1:]) == got


Q8_QUERIES = [
    "SELECT level, count(*), sum(n_slots) FROM token_dict GROUP BY level ORDER BY level",
    "SELECT strand, count(*), sum(\"end\" - start) FROM token_postings GROUP BY 1 ORDER BY 1",
    "SELECT is_circular, count(*), sum(n_slots), max(contig_len) FROM contig_tapes GROUP BY 1 ORDER BY 1",
    "SELECT is_circular, count(DISTINCT contig_uid), count(*) FROM token_postings WHERE ordinal <= 2 GROUP BY 1 ORDER BY 1",
    "SELECT strand, count(*), sum(\"end\" - start) FROM protein_occurrence GROUP BY 1 ORDER BY 1",
    "SELECT model_acc, count(DISTINCT protein_uid) FROM l1_domains GROUP BY 1 ORDER BY 1",
    "SELECT d.level, count(*) FROM token_postings p JOIN token_dict d ON d.token_id = p.token_id "
    "GROUP BY 1 ORDER BY 1",
]


@pytest.mark.covers("Q8")
@pytest.mark.parametrize("sql", Q8_QUERIES)
def test_q8_ad_hoc_sql_equals_an_unindexed_copy(built, sql):
    """The catalog's answer equals the same SQL over plain read_parquet views of a copy of the files."""
    root, man, con, _ = built
    plain = duckdb.connect()
    by_table = defaultdict(list)
    for f in man["files"]:
        by_table[f["table"]].append(str(root / f["path"]))
    for name, files in by_table.items():
        plain.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet({files!r}, hive_partitioning = true, "
                      "union_by_name = true)")
    assert con.execute(sql).fetchall() == plain.execute(sql).fetchall()
    plain.close()


# --- contracts and layout ------------------------------------------------------------------


@pytest.mark.covers("C-04", "C-07", "C-03", "C-05")
def test_l3_contracts_hold(built):
    _, man, con, _ = built
    assert contracts.check_list_lengths(con) == []
    assert contracts.check_saturation(con) == []
    assert contracts.check_foreign_keys(con) == []
    assert contracts.check_coordinates(con) == []
    assert contracts.check_sort_order(man) == [] and contracts.check_statistics(man) == []
    assert contracts.check_file_sizes(man) == []


@pytest.mark.covers("C-04", "C-07")
def test_l3_faults_are_caught(built):
    _, _, con, _ = built
    t = con.sql("SELECT * EXCLUDE (tape_version, stage) FROM contig_tapes").to_arrow_table().to_pylist()
    short = [dict(r) for r in t]
    short[0]["strands"] = short[0]["strands"][:-1]
    assert any("parallel lists" in p for p in contracts.check_list_lengths(_with_tv(con, short)))
    flagged = [dict(r) for r in t]
    lin = next(r for r in flagged if r["orfs_to_contig_starts"] and r["n_slots"])
    lin["slot_flags"] = [1] * lin["n_slots"]  # flagged as clipped but not saturated
    assert any("not saturated" in p for p in contracts.check_saturation(_with_tv(con, flagged)))
    # tokens whose lists differ in length from their tape (joined on contig_uid and tape_version)
    shorter = [dict(r) for r in t]
    longer = next(r for r in shorter if r["n_slots"])
    longer["n_slots"] += 1
    for col in table("contig_tapes").list_columns:
        if longer[col] is not None:
            longer[col] = longer[col] + longer[col][:1]
    assert any("differ in length from their tape" in p for p in contracts.check_list_lengths(_with_tv(con, shorter)))


def _with_tv(con, tape_rows):
    c = duckdb.connect()
    tbl = pa.Table.from_pylist(tape_rows)
    c.register("_t", tbl)
    c.execute(f"CREATE VIEW contig_tapes AS SELECT *, '{TV}' AS tape_version FROM _t")
    c.register("tape_tokens", con.sql("SELECT * FROM tape_tokens").to_arrow_table())
    return c


@pytest.mark.covers("LY-06", "LY-05")
@pytest.mark.parametrize("width", [1, 2])
def test_ly06_pair_lookup_opens_one_bucket(built, tmp_path, width):
    """Q1 and Q4 over postings bucketed by the token's low byte: each side of the self-join opens one file."""
    root, _, _, _ = built
    copy = tmp_path / "s"
    shutil.copytree(root, copy, ignore=shutil.ignore_patterns("token_postings"))
    postings.run(copy, REL, ["S0"], prefix_width=width)
    man = catalog.build_manifest(copy, REL, ["S0"])  # LY-05: every row's token matches its bucket
    files = [f for f in man["files"] if f["table"] == "token_postings"]
    assert len(files) > 1
    full = _oracle_pairs(_oracle_rows(root))
    bucket = lambda t: f"{t & 0xFF:02x}"[:width]  # noqa: E731
    a, b = next(p for p in sorted(full) if bucket(p[0]) != bucket(p[1]))
    c = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / "c.duckdb"))
    try:
        for sql, want in (("SELECT q1_count(?, ?, 10, NULL)", _loci(full[(a, b)])),
                          ("SELECT count(*) FROM pairs_of(?, ?, 10)", len(full[(a, b)]))):
            prof = catalog.profile_query(c, sql, [a, b])
            scans = [s for s in prof["scans"] if any("/token_postings/" in f for f in s["files"])]
            assert prof["result"][0][0] == want >= 1, sql
            assert len(scans) == 2, sql  # one per side
            assert sorted(s["files"][0].split("token_bucket=")[1].split("/")[0] for s in scans) == \
                sorted([bucket(a), bucket(b)]), sql
            assert all(s["files_read"] == 1 and s["row_groups_scanned"] == 1 for s in scans), (sql, scans)
    finally:
        c.close()


@pytest.mark.covers("K-07")
def test_k07_list_projection_cost_is_recorded(built, record_property):
    """Bytes decoded per Q2 fetch grow with the list columns projected; the curve is recorded (decision at S2).

    Decoded bytes are the projected columns' uncompressed chunk sizes in the row groups the
    profiler says were scanned. At micro scale DuckDB fetches a small file whole, so bytes
    read do not separate the columns; they are recorded alongside.
    """
    root, man, con, _ = built
    files = [root / f["path"] for f in man["files"] if f["table"] == "contig_tapes"]
    cuids = [r[0] for r in con.execute("SELECT contig_uid FROM contig_tapes").fetchall()]
    lists = list(table("contig_tapes").list_columns)
    curve = {}
    for n in (1, 3, len(lists)):
        cols = ["contig_uid"] + lists[:n]
        prof = catalog.profile_query(con, f"SELECT {', '.join(cols)} FROM contig_tapes WHERE list_contains(?, contig_uid)",
                                     [cuids])
        (scan,) = [s for s in prof["scans"] if any("/l3/tapes/" in f for f in s["files"])]
        decoded = 0
        for f in files:
            md = pq.ParquetFile(f).metadata
            for rg in range(md.num_row_groups):
                for i in range(md.num_columns):
                    if md.schema.column(i).path.split(".")[0] in cols:
                        decoded += md.row_group(rg).column(i).total_uncompressed_size
        curve[n] = {"bytes_decoded": decoded, "bytes_read": prof["bytes_read"], "row_groups": scan["row_groups_scanned"]}
    record_property("k07_curve", curve)
    assert curve[1]["bytes_decoded"] < curve[3]["bytes_decoded"] < curve[len(lists)]["bytes_decoded"]
    assert len({c["row_groups"] for c in curve.values()}) == 1  # projection does not change the row groups read


@pytest.mark.covers("D-01", "D-04", "LY-08")
def test_l3_is_deterministic_and_order_free(built, tmp_path):
    root, man, _, _ = built
    other = tmp_path / "b"
    build(other, order=-1)
    h1 = contracts.content_hashes(man)
    h2 = contracts.content_hashes(catalog.build_manifest(other, REL, ["S0"]))
    for name in ("contig_tapes", "tape_tokens", "token_dict", "token_postings", "tape_version_def"):
        assert h1[name] == h2[name], name
    for rel in sorted(p.relative_to(root).as_posix() for p in (root / "l3").rglob("*.parquet")):
        assert (root / rel).read_bytes() == (other / rel).read_bytes(), rel


# --- card B-07: tapes by tape version --------------------------------------------------------------

TAPE_DIRS = ("l3/tapes", "l3/tape_version_def")


def _tapes_of_tv(root: Path, tv: str = TV) -> dict[bytes, dict]:
    return {r["contig_uid"]: r for r in _read(root, f"l3/tapes/tape_version={tv}/stage=S0/part-*.parquet")}


def _drop_tape_ledgers(root: Path) -> None:
    """The tapes' ledgers go with the tapes: the ledger writer never rewrites a step either (card B-19)."""
    for d in sorted(root.glob("ledger/stage=*/step=tape-*")):
        shutil.rmtree(d)


def _fresh(root: Path, tmp_path: Path, name: str, drop=TAPE_DIRS) -> Path:
    """A copy of the built root without the tables the test rebuilds (the writer never overwrites)."""
    fresh = tmp_path / name
    shutil.copytree(root, fresh)
    for d in drop:
        shutil.rmtree(fresh / d)
    if "l3/tapes" in drop:
        _drop_tape_ledgers(fresh)
    return fresh


def _digests(root: Path, dirs) -> dict[str, str]:
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for d in dirs for p in sorted((root / d).rglob("*.parquet"))}


@pytest.mark.covers("C-04")
def test_b07_feature_uids_follow_the_genes_in_coordinate_order(built):
    """With the gene-call producer alone, slot i of a tape is the i-th gene of its contig in (start, end, strand)
    order (the rule the legacy slot ordinal followed until card B-17 removed it)."""
    root, _, _, _ = built
    occ = defaultdict(list)
    for r in _read(root, f"l1/protein_occurrence/gene_call_version={V}/stage=S0/part-0.parquet"):
        occ[r["contig_uid"]].append(r)
    feats = {r["feature_uid"]: r for r in _read(root, f"l1/feature/producer={V}/stage=S0/part-*.parquet")}
    tps = _tapes_of_tv(root)
    assert len(tps) == len(CONTRACT["objects"]) and all(t["n_slots"] for t in tps.values())
    for cu, t in tps.items():
        genes = sorted(occ[cu], key=lambda r: (r["start"], r["end"], r["strand"]))
        assert len(genes) == t["n_slots"]
        assert t["feature_uids"] == [r["feature_uid"] for r in genes]
        assert t["protein_uids"] == [r["protein_uid"] for r in genes]
        assert t["feature_kinds"] == [ids.KIND["cds"]] * t["n_slots"]
        for fu, s, e, strand in zip(t["feature_uids"], t["start_bps"], t["end_bps"], t["strands"]):
            f = feats[fu]
            assert (f["contig_uid"], f["start"], f["end"], f["strand"], f["kind"]) == (cu, s, e, strand, ids.KIND["cds"])
            assert fu == ids.feature_uid(cu, s, e, strand, ids.KIND["cds"])


@pytest.mark.covers("C-01", "LY-05", "LY-10")
def test_b07_layout_conforms_and_the_manifest_accepts_it(built):
    root, man, con, _ = built
    assert tapes.tape_version([V]) == TV
    by = defaultdict(list)
    for f in man["files"]:
        by[f["table"]].append(f)
    assert [f["path"] for f in by["tape_version_def"]] == [f"l3/tape_version_def/tape_version={TV}/part-0.parquet"]
    assert by["contig_tapes"] and all(f["partitions"] == {"tape_version": TV, "stage": "S0"} for f in by["contig_tapes"])
    assert by["tape_tokens"] and all(f["partitions"] == {"tape_version": TV, "release": REL} for f in by["tape_tokens"])
    for name in ("contig_tapes", "tape_tokens", "tape_version_def"):
        for f in by[name]:
            assert conformance_problems(pq.read_table(root / f["path"], partitioning=None), table(name)) == [], f["path"]
    assert pq.read_table(root / by["tape_version_def"][0]["path"], partitioning=None).to_pylist() == [{"producer": V}]
    assert pq.read_schema(root / by["contig_tapes"][0]["path"]).names[-1] == "feature_uids"
    assert con.execute("SELECT tape_version, producer FROM tape_version_def").fetchall() == [(TV, V)]
    assert con.execute("SELECT DISTINCT tape_version FROM contig_tapes").fetchall() == [(TV,)]
    assert con.execute("SELECT DISTINCT tape_version FROM tape_tokens").fetchall() == [(TV,)]
    assert contracts.check_append_only(man, man) == []


@pytest.mark.covers("LY-04", "LY-06", "Q2")
def test_b07_slot_bounded_row_groups_and_a_one_contig_lookup_reads_one(built, tmp_path, record_property):
    root, _, _, circ = built
    old = _tapes_of_tv(root)
    n_by = {cu: t["n_slots"] for cu, t in old.items()}
    cap = min(n_by.values())  # no two contigs fit together, so each is its own row group; the others exceed the cap
    fresh = _fresh(root, tmp_path, "s")
    res = tapes.run(fresh, "S0", circ, max_slots=cap)
    assert res["tape_version"] == TV and res["files"] == 1 and res["contigs"] == len(n_by)
    [path] = sorted((fresh / f"l3/tapes/tape_version={TV}/stage=S0").glob("*.parquet"))
    pf = pq.ParquetFile(path)
    groups = [pf.read_row_group(i, columns=["n_slots"]).column(0).to_pylist() for i in range(pf.metadata.num_row_groups)]
    assert len(groups) == len(n_by) >= 3
    assert all(sum(g) <= cap or len(g) == 1 for g in groups)
    assert any(sum(g) > cap for g in groups)  # an oversized contig is a row group of its own
    assert pq.read_table(path).to_pylist() == [old[cu] for cu in sorted(old)]
    # the default build: every row group within AGENTS §4.6
    for f in sorted((root / "l3/tapes").rglob("*.parquet")):
        p = pq.ParquetFile(f)
        assert all(sum(p.read_row_group(i, columns=["n_slots"]).column(0).to_pylist()) <= 131_072
                   for i in range(p.metadata.num_row_groups))
    # tape tokens with one contig per row group, so that their pruning is visible too
    [tok] = sorted((fresh / f"l3/tape_tokens/tape_version={TV}/release={REL}").glob("*.parquet"))
    pq.write_table(pq.read_table(tok, partitioning=None), tok, row_group_size=1, write_statistics=["contig_uid"])
    man = catalog.build_manifest(fresh, REL, ["S0"])
    assert contracts.check_sort_order(man) == [] and contracts.check_statistics(man) == []
    con = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / "c.duckdb"))

    def scanned(prof, d):
        (s,) = [s for s in prof["scans"] if any(f"/l3/{d}/" in x for x in s["files"])]
        return s["row_groups_scanned"], s["row_groups_total"]

    try:
        sql = "SELECT sum(n_slots), sum(len(nr90_tokens)) FROM tapes_of({})"
        for cu in sorted(n_by):
            for arg, params in (("?", [[cu]]),  # a literal key, then a key that comes from a joined table
                                ("(SELECT list(contig_uid) FROM contig_index WHERE contig_uid = ?)", [cu])):
                prof = catalog.profile_query(con, sql.format(arg), params)
                assert prof["result"][0] == (n_by[cu], n_by[cu]), arg
                assert scanned(prof, "tapes") == (1, len(groups)), arg
                assert scanned(prof, "tape_tokens") == (1, len(n_by)), arg
        prof = catalog.profile_query(con, sql.format("?"), [sorted(n_by)])
        record_property("b07_three_contig_lookup", {"tapes": scanned(prof, "tapes"),
                                                    "tape_tokens": scanned(prof, "tape_tokens")})
        assert scanned(prof, "tapes")[0] <= len(n_by)
    finally:
        con.close()


@pytest.mark.covers("LY-04")
def test_b07_parts_are_cut_at_max_file_bytes(built, tmp_path):
    root, _, _, circ = built
    old = _tapes_of_tv(root)
    order = sorted(old)
    est = [200 + 30 * old[cu]["n_slots"] for cu in order]  # the builder's estimate per contig row
    fresh = _fresh(root, tmp_path, "s")
    res = tapes.run(fresh, "S0", circ, max_file_bytes=est[0] + est[1])  # the first two rows fill part 0
    assert res["files"] == 2
    d = fresh / f"l3/tapes/tape_version={TV}/stage=S0"
    assert sorted(p.name for p in d.iterdir()) == ["part-0.parquet", "part-1.parquet"]
    assert [r["contig_uid"] for r in pq.read_table(d / "part-0.parquet").to_pylist()] == order[:2]
    assert [r["contig_uid"] for r in pq.read_table(d / "part-1.parquet").to_pylist()] == order[2:]
    assert _tapes_of_tv(fresh) == old


@pytest.mark.covers("R-09", "C-04", "C-05")
def test_b07_second_producer_is_a_new_tape_version_beside_the_first(built, tmp_path):
    root, _, _, circ = built
    fresh = tmp_path / "s"
    shutil.copytree(root, fresh)
    first = _digests(fresh, ("l3",))
    assert any("/tape_version_def/" in k for k in first) and any("/tapes/" in k for k in first)
    mg = next(iter(circ))  # the only circular contig is M. genitalium (NC_000908.2)
    genes = sorted((r for r in _read(fresh, f"l1/feature/producer={V}/stage=S0/part-*.parquet") if r["contig_uid"] == mg),
                   key=lambda r: (r["start"], r["end"]))
    reach, gap = 0, None
    for g in genes:
        if g["start"] - reach >= 100 and reach > 0:
            gap = reach
            break
        reach = max(reach, g["end"])
    assert gap is not None
    s, e, kind = gap + 30, gap + 70, ids.KIND["ncrna"]
    fu = ids.feature_uid(mg, s, e, 1, kind)
    spec = table("feature")
    io.write_table([pa.RecordBatch.from_pylist(
        [{"feature_uid": fu, "contig_uid": mg, "start": s, "end": e, "strand": 1, "kind": kind}], schema=spec.schema)],
        spec, fresh / "l1/feature/producer=test-ncrna/stage=S0/part-0.parquet")
    res = tapes.run(fresh, "S0", circ, producers=["test-ncrna", V])
    tv2 = "t-" + hashlib.sha256("\n".join(sorted([V, "test-ncrna"])).encode("utf-8")).hexdigest()[:12]
    assert res["tape_version"] == tv2 != TV
    assert res["slots"] == sum(t["n_slots"] for t in _tapes_of_tv(root).values()) + 1
    # the first tape version's files (and its tokens) are untouched
    assert {k: v for k, v in _digests(fresh, ("l3",)).items() if k in first} == first
    old, new = _tapes_of_tv(fresh), _tapes_of_tv(fresh, tv2)
    assert set(new) == set(old)
    for cu in set(old) - {mg}:
        assert new[cu] == old[cu]
    o, n = old[mg], new[mg]
    assert n["n_slots"] == o["n_slots"] + 1
    i = n["feature_uids"].index(fu)
    assert (n["start_bps"][i], n["end_bps"][i], n["strands"][i]) == (s, e, 1)
    assert n["protein_uids"][i] is None and n["feature_kinds"][i] == kind
    drop = lambda xs: xs[:i] + xs[i + 1:]  # noqa: E731
    for col in ("feature_uids", "protein_uids", "feature_kinds", "start_bps", "end_bps", "strands"):
        assert drop(n[col]) == o[col], col
    keys = list(zip(n["start_bps"], n["end_bps"], n["strands"], n["feature_kinds"]))
    assert keys == sorted(keys)
    assert pq.read_table(fresh / f"l3/tape_version_def/tape_version={tv2}/part-0.parquet").column(
        "producer").to_pylist() == sorted([V, "test-ncrna"])
    man = catalog.build_manifest(fresh, REL, ["S0"])
    assert contracts.check_sort_order(man) == [] and contracts.check_statistics(man) == []
    con = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / "c.duckdb"))
    try:
        assert contracts.check_list_lengths(con) == [] and contracts.check_coordinates(con) == []
        assert contracts.check_foreign_keys(con) == []
        got = con.execute("SELECT tape_version, n_slots, nr90_tokens IS NULL FROM tapes_of(?)", [[mg]]).fetchall()
        assert got == sorted([(TV, o["n_slots"], False), (tv2, n["n_slots"], True)])  # tokens exist for TV only
    finally:
        con.close()


@pytest.mark.covers("C-04")
def test_b07_a_contig_without_features_has_an_empty_tape(built, tmp_path):
    root, _, _, circ = built
    fresh = tmp_path / "s"
    shutil.copytree(root, fresh)
    io.write_table([], table("feature"), fresh / "l1/feature/producer=test-none/stage=S0/part-0.parquet")
    res = tapes.run(fresh, "S0", circ, producers=["test-none"])
    got, old = _tapes_of_tv(fresh, res["tape_version"]), _tapes_of_tv(root)
    assert set(got) == set(old) and res["slots"] == 0 and res["contigs"] == len(old)
    for cu, t in got.items():
        assert (t["n_slots"], t["feature_uids"], t["protein_uids"], t["feature_kinds"]) == (0, [], [], [])
        assert (t["contig_len"], t["is_circular"]) == (old[cu]["contig_len"], cu in circ)
        assert set(t) == set(table("contig_tapes").schema.names)  # card B-17: no dimension column on a tape


def test_b07_tape_version_def_refuses_a_mismatching_producer_list(built, tmp_path):
    root, _, _, circ = built
    fresh = _fresh(root, tmp_path, "s")
    spec = table("tape_version_def")
    io.write_table([pa.table({"producer": ["somebody-else"]})], spec,
                   fresh / f"l3/tape_version_def/tape_version={TV}/part-0.parquet")
    with pytest.raises(ValueError, match="tape_version_def"):
        tapes.run(fresh, "S0", circ)
    assert not (fresh / "l3/tapes").exists()
    # the matching definition is accepted and left as it is
    shutil.rmtree(fresh / "l3/tape_version_def")
    io.write_table([pa.table({"producer": [V]})], spec, fresh / f"l3/tape_version_def/tape_version={TV}/part-0.parquet")
    kept = _digests(fresh, ("l3/tape_version_def",))
    tapes.run(fresh, "S0", circ)
    assert _digests(fresh, ("l3/tape_version_def",)) == kept == _digests(root, ("l3/tape_version_def",))


@pytest.mark.covers("D-01", "D-04", "LY-08")
def test_b07_tape_bytes_do_not_depend_on_input_row_order_or_the_run(built, tmp_path):
    root, _, _, circ = built
    want = _digests(root, TAPE_DIRS)
    rng = random.Random(3)
    for name, shuffle in (("again", False), ("shuffled", True)):
        fresh = _fresh(root, tmp_path, name)
        if shuffle:
            for rel in (f"l1/feature/producer={V}/stage=S0/part-0.parquet",
                        f"l1/feature_protein/gene_call_version={V}/stage=S0/part-0.parquet",
                        f"l1/protein_occurrence/gene_call_version={V}/stage=S0/part-0.parquet",
                        "l0seq/contig_index/stage=S0/part-0.parquet"):
                t = pq.read_table(fresh / rel, partitioning=None)
                idx = list(range(t.num_rows))
                rng.shuffle(idx)
                pq.write_table(t.take(idx), fresh / rel, row_group_size=max(1, t.num_rows // 7))
        tapes.run(fresh, "S0", circ)
        assert _digests(fresh, TAPE_DIRS) == want, name


def test_b07_default_topology_comes_from_genome_input(built, tmp_path):
    root, _, _, circ = built
    fresh = _fresh(root, tmp_path, "s")
    write_genome_input(fresh, "S0", [GenomeInput((cu,), circular=(cu,) if cu in circ else ()) for cu in _tapes_of_tv(root)])
    tapes.run(fresh, "S0")  # circular=None
    assert _digests(fresh, TAPE_DIRS) == _digests(root, TAPE_DIRS)


# --- card B-09: tokens joined, ranked and aggregated in DuckDB -------------------------------------

TOKEN_DIRS = ("l3/token_dict", "l3/tape_tokens")
LEVEL_BITS = {"family": 1, "nr50": 2, "nr90": 3}  # PLAN.md §4.1: token = level (3) | source (5) | local (56)


def _rows(root: Path, glob: str, columns=None) -> list[dict]:
    """The rows of the files of `glob` in path order, without partition columns."""
    return [r for f in sorted(root.glob(glob)) for r in pq.read_table(f, columns=columns, partitioning=None).to_pylist()]


def _oracle_tokens(root: Path, release: str, libraries: list[str], tv: str = TV) -> tuple[list[dict], list[dict]]:
    """The in-memory algorithm of tokens.run before card B-09, kept here as a brute-force reference:
    every tape, domain and cluster row in Python dicts. Returns (token_dict rows, tape_tokens rows)."""
    tps = _rows(root, f"l3/tapes/tape_version={tv}/stage=*/part-*.parquet", ["contig_uid", "protein_uids"])
    best: dict[bytes, tuple] = {}
    for rank, lib in enumerate(libraries):
        for r in _rows(root, f"l1/domains/library={lib}/stage=*/part-*.parquet"):
            key = (rank, -r["dom_bitscore"], r["model_acc"], r["env_start"])
            if r["protein_uid"] not in best or key < best[r["protein_uid"]][0]:
                best[r["protein_uid"]] = (key, lib, r["model_acc"])
    clusters: dict[str, dict[bytes, bytes]] = {"nr50": {}, "nr90": {}}
    for r in _rows(root, f"l2/cluster_member/release={release}/part-*.parquet"):
        if r["level"] in clusters:
            clusters[r["level"]][r["member_uid"]] = r["cluster_uid"]

    def keys(pu):
        fam = best.get(pu) if pu is not None else None
        out = {"family": (libraries.index(fam[1]) + 1, fam[2]) if fam else None}
        for lvl in ("nr50", "nr90"):
            c = clusters[lvl].get(pu) if pu is not None else None
            out[lvl] = (0, c.hex()) if c else None
        return out

    per_contig = [(t["contig_uid"], [keys(pu) for pu in t["protein_uids"]]) for t in tps]
    uses = Counter((lvl, *k) for _, slots in per_contig for s in slots for lvl, k in s.items() if k is not None)
    tid: dict[tuple, int] = {}
    for lvl, bits in LEVEL_BITS.items():
        for src in sorted({s for (lv, s, _) in uses if lv == lvl}):
            for local, key in enumerate(sorted(k for (lv, s, k) in uses if lv == lvl and s == src)):
                tid[(lvl, src, key)] = bits << 61 | src << 56 | local
    names = {0: "mmseqs2", **{i + 1: lib for i, lib in enumerate(libraries)}}
    dict_rows = sorted(({"token_id": t, "level": lv, "source": names[s], "key": k, "n_slots": uses[(lv, s, k)]}
                        for (lv, s, k), t in tid.items()), key=lambda r: r["token_id"])
    tok_rows = [{"contig_uid": cu, **{f"{lvl}_tokens": [None if s[lvl] is None else tid[(lvl, *s[lvl])] for s in slots]
                                      for lvl in LEVEL_BITS}}
                for cu, slots in sorted(per_contig, key=lambda x: x[0])]
    return dict_rows, tok_rows


def _written_tokens(root: Path, release: str, tv: str = TV) -> tuple[list[dict], list[dict]]:
    return (_rows(root, f"l3/token_dict/release={release}/part-*.parquet"),
            _rows(root, f"l3/tape_tokens/tape_version={tv}/release={release}/part-*.parquet"))


def _row_group_slots(path: Path) -> list[list[int]]:
    pf = pq.ParquetFile(path)
    return [[len(x) for x in pf.read_row_group(i, columns=["nr90_tokens"]).column(0).to_pylist()]
            for i in range(pf.metadata.num_row_groups)]


@pytest.mark.covers("C-04")
def test_b09_tokens_equal_the_in_memory_oracle(built):
    root, _, _, _ = built
    got = _written_tokens(root, REL)
    assert got == _oracle_tokens(root, REL, ["test-1"])
    dict_rows, tok_rows = got
    assert {r["level"] for r in dict_rows} == set(LEVELS) and len(tok_rows) == len(CONTRACT["objects"])
    fams = [x for r in tok_rows for x in r["family_tokens"]]
    assert None in fams and any(x is not None for x in fams)  # slots without a family keep a null in the list
    n_slots = {cu: t["n_slots"] for cu, t in _tapes_of_tv(root).items()}
    for r in tok_rows:
        assert all(len(r[f"{lvl}_tokens"]) == n_slots[r["contig_uid"]] for lvl in LEVELS)


def _dom(pu: bytes, acc: str, score: float, env: int) -> dict:
    return {"protein_uid": pu, "model_acc": acc, "model_name": acc.lower(), "model_len": 100, "seq_bitscore": score,
            "seq_evalue_zref": 1e-9, "dom_bitscore": score, "dom_ievalue_zref": 1e-9, "env_start": env,
            "env_end": env + 20, "ali_start": env, "ali_end": env + 20, "hmm_from": 0, "hmm_to": 20,
            "cutoff": "gathering", "z": 1e8, "dom_z": 1e8}


def _write_sorted(spec_name: str, rows: list[dict], path: Path) -> None:
    spec = table(spec_name)
    rows = sorted(rows, key=lambda r: tuple(r[c] for c in spec.sort_key if c in spec.schema.names))
    io.write_table([pa.Table.from_pylist(rows, schema=spec.schema)], spec, path)


@pytest.mark.covers("C-04")
def test_b09_canonical_families_and_cluster_keys_on_synthetic_mappings(built, tmp_path):
    """Library precedence, bitscore and accession ties, key order and unused keys, on real protein uids."""
    root, _, _, _ = built
    fresh = tmp_path / "s"
    shutil.copytree(root, fresh)
    rel = "rS0.b09"
    tps = _tapes_of_tv(root)
    p = sorted({pu for t in tps.values() for pu in t["protein_uids"]})[:8]
    libs = {
        "zzz-1": [_dom(p[0], "M2", 50.0, 10), _dom(p[0], "M1", 50.0, 40),  # equal bitscores: the smaller accession
                  _dom(p[1], "Z9", 40.0, 5), _dom(p[1], "A1", 60.0, 90),  # the higher bitscore before the accession
                  _dom(p[2], "Q1", -5.0, 3), _dom(p[2], "Q2", 3.0, 7),  # negative bitscores
                  _dom(p[3], "b10", 20.0, 1), _dom(p[3], "b10", 20.0, 0)],  # one model twice
        "mmm-1": [_dom(p[0], "TOP", 999.0, 1), _dom(p[4], "b9", 10.0, 2), _dom(p[5], "B2", 10.0, 2)],
        "aaa-1": [_dom(p[6], "a", 1.0, 0), _dom(p[4], "zzz", 500.0, 0)],
    }
    for lib, rows in libs.items():
        _write_sorted("l1_domains", rows, fresh / f"l1/domains/library={lib}/stage=S0/part-0.parquet")
    order = ["zzz-1", "mmm-1", "aaa-1"]  # precedence is this order, not the names' order
    want_family = {p[0]: ("zzz-1", "M1"), p[1]: ("zzz-1", "A1"), p[2]: ("zzz-1", "Q2"), p[3]: ("zzz-1", "b10"),
                   p[4]: ("mmm-1", "b9"), p[5]: ("mmm-1", "B2"), p[6]: ("aaa-1", "a")}
    hi, mid, lo = b"\xff" * 16, b"\x80" + b"\0" * 15, b"\x01" * 16  # unsigned order: lo < mid < hi
    cl = [("nr90", p[0], hi), ("nr90", p[1], mid), ("nr90", p[2], lo), ("nr90", p[3], hi),
          ("nr50", p[0], b"\x7f" * 16), ("nr50", p[4], b"\0" * 16), ("subfam", p[5], b"\x05" * 16)]
    _write_sorted("cluster_member", [{"member_uid": m, "level": lv, "cluster_uid": c, "rep_uid": m,
                                      "cluster_params_hash": "synthetic"} for lv, m, c in cl],
                  fresh / f"l2/cluster_member/release={rel}/part-0.parquet")
    res = tokens.run(fresh, rel, ["S0"], libraries=order)
    dict_rows, tok_rows = got = _written_tokens(fresh, rel)
    assert got == _oracle_tokens(fresh, rel, order)
    by_key = {(r["level"], r["source"], r["key"]): r for r in dict_rows}
    uses = Counter(pu for t in tps.values() for pu in t["protein_uids"])
    fam = [k for k in by_key if k[0] == "family"]
    assert sorted(fam, key=lambda k: by_key[k]["token_id"]) == [  # by source, then key in code point order
        ("family", "zzz-1", "A1"), ("family", "zzz-1", "M1"), ("family", "zzz-1", "Q2"), ("family", "zzz-1", "b10"),
        ("family", "mmm-1", "B2"), ("family", "mmm-1", "b9"), ("family", "aaa-1", "a")]  # "TOP" is on no slot
    for (lv, src, key), r in by_key.items():
        tid = r["token_id"]
        if lv == "family":
            assert tid >> 61 == 1 and (tid >> 56) & 31 == order.index(src) + 1
            assert r["n_slots"] == sum(uses[pu] for pu, f in want_family.items() if f == (src, key))
    assert [(by_key[("nr90", "mmseqs2", c.hex())]["token_id"]) for c in (lo, mid, hi)] == [3 << 61 | i for i in range(3)]
    assert by_key[("nr90", "mmseqs2", hi.hex())]["n_slots"] == uses[p[0]] + uses[p[3]]
    assert res["tokens"] == {"family": 7, "nr50": 2, "nr90": 3} and res["contigs"] == len(tps)
    tid_of = {(r["level"], r["token_id"]): (r["source"], r["key"]) for r in dict_rows}
    for r in tok_rows:
        for i, pu in enumerate(tps[r["contig_uid"]]["protein_uids"]):
            f = r["family_tokens"][i]
            assert (tid_of[("family", f)] if f is not None else None) == want_family.get(pu)
            n50 = r["nr50_tokens"][i]
            assert (n50 is None) == (pu not in (p[0], p[4]))  # null where the protein has no cluster at that level


def _add_ncrna(root: Path, mg: bytes) -> bytes:
    """One kind-2 feature in an intergenic gap of contig `mg` as producer `test-ncrna` (as in the B-07 test)."""
    genes = sorted((r for r in _read(root, f"l1/feature/producer={V}/stage=S0/part-*.parquet") if r["contig_uid"] == mg),
                   key=lambda r: (r["start"], r["end"]))
    reach, gap = 0, None
    for g in genes:
        if g["start"] - reach >= 100 and reach > 0:
            gap = reach
            break
        reach = max(reach, g["end"])
    assert gap is not None
    s, e, kind = gap + 30, gap + 70, ids.KIND["ncrna"]
    fu = ids.feature_uid(mg, s, e, 1, kind)
    spec = table("feature")
    io.write_table([pa.RecordBatch.from_pylist(
        [{"feature_uid": fu, "contig_uid": mg, "start": s, "end": e, "strand": 1, "kind": kind}], schema=spec.schema)],
        spec, root / "l1/feature/producer=test-ncrna/stage=S0/part-0.parquet")
    return fu


@pytest.mark.covers("C-04")
def test_b09_a_non_cds_slot_has_null_tokens_at_every_level(built, tmp_path):
    root, _, _, circ = built
    fresh = tmp_path / "s"
    shutil.copytree(root, fresh)
    mg = next(iter(circ))  # the only circular contig is M. genitalium (NC_000908.2)
    fu = _add_ncrna(fresh, mg)
    producers = ["test-ncrna", V]
    tv2 = tapes.run(fresh, "S0", circ, producers=producers)["tape_version"]
    rel = "rS0.b09nc"
    shutil.copytree(fresh / f"l2/cluster_member/release={REL}", fresh / f"l2/cluster_member/release={rel}")
    res = tokens.run(fresh, rel, ["S0"], producers=producers)
    got = _written_tokens(fresh, rel, tv2)
    assert got == _oracle_tokens(fresh, rel, ["test-1"], tv2)
    old_dict, old_tok = _written_tokens(root, REL)
    assert got[0] == old_dict  # one more slot without a protein: the same tokens and counts
    old, new = ({r["contig_uid"]: r for r in rows} for rows in (old_tok, got[1]))
    i = _tapes_of_tv(fresh, tv2)[mg]["feature_uids"].index(fu)
    for lvl in LEVELS:
        col = f"{lvl}_tokens"
        assert new[mg][col][i] is None
        assert new[mg][col][:i] + new[mg][col][i + 1:] == old[mg][col]
        assert all(new[cu][col] == old[cu][col] for cu in set(old) - {mg})
    assert res["contigs"] == len(old)


@pytest.mark.covers("C-04")
def test_b09_tapes_without_slots_get_three_empty_lists(built, tmp_path):
    root, _, _, circ = built
    fresh = tmp_path / "s"
    shutil.copytree(root, fresh)
    io.write_table([], table("feature"), fresh / "l1/feature/producer=test-none/stage=S0/part-0.parquet")
    tv0 = tapes.run(fresh, "S0", circ, producers=["test-none"])["tape_version"]
    res = tokens.run(fresh, REL + "-empty", ["S0"], producers=["test-none"])
    dict_rows, tok_rows = got = _written_tokens(fresh, REL + "-empty", tv0)
    assert got == _oracle_tokens(fresh, REL + "-empty", ["test-1"], tv0)
    assert dict_rows == []
    assert tok_rows == [{"contig_uid": cu, "family_tokens": [], "nr50_tokens": [], "nr90_tokens": []}
                        for cu in sorted(_tapes_of_tv(root))]
    assert res == {"tokens": {"family": 0, "nr50": 0, "nr90": 0}, "contigs": len(tok_rows), "files": 2}


@pytest.mark.covers("LY-04")
def test_b09_tape_token_row_groups_hold_whole_contigs_within_max_slots(built, tmp_path):
    root, _, _, _ = built
    n_by = {cu: t["n_slots"] for cu, t in _tapes_of_tv(root).items()}
    cap = min(n_by.values())  # no two contigs fit together; the larger ones exceed the cap alone
    fresh = _fresh(root, tmp_path, "s", TOKEN_DIRS)
    tokens.run(fresh, REL, ["S0"], max_slots=cap)
    [path] = sorted((fresh / f"l3/tape_tokens/tape_version={TV}/release={REL}").glob("*.parquet"))
    groups = _row_group_slots(path)
    assert groups == [[n_by[cu]] for cu in sorted(n_by)]
    assert all(sum(g) <= cap or len(g) == 1 for g in groups) and any(sum(g) > cap for g in groups)
    assert _written_tokens(fresh, REL) == _written_tokens(root, REL)
    for f in sorted((root / "l3/tape_tokens").rglob("*.parquet")):  # the default build, AGENTS §4.6
        assert all(sum(g) <= 131_072 for g in _row_group_slots(f))


@pytest.mark.covers("LY-04")
def test_b09_parts_are_cut_at_max_file_bytes(built, tmp_path):
    root, _, _, _ = built
    want_dict, want_tok = _written_tokens(root, REL)
    n_slots = [t["n_slots"] for t in _tapes_of_tv(root).values()]
    est_tok = 100 + 30 * max(1, sum(n_slots) // len(n_slots))  # the estimates per row, as before card B-09
    fresh = _fresh(root, tmp_path, "s", TOKEN_DIRS)
    cap = 3_000  # 50 dictionary rows per part; every tape-token row is a part of its own
    res = tokens.run(fresh, REL, ["S0"], max_file_bytes=cap)
    for d, rows, est in ((f"l3/token_dict/release={REL}", want_dict, 60),
                         (f"l3/tape_tokens/tape_version={TV}/release={REL}", want_tok, est_tok)):
        starts = layout.cut([est] * len(rows), cap)
        assert len(starts) > 1
        names = sorted((fresh / d).iterdir(), key=lambda p: int(p.stem.removeprefix("part-")))
        assert [p.name for p in names] == [f"part-{k}.parquet" for k in range(len(starts))]
        for p, a, b in zip(names, starts, starts[1:] + [len(rows)]):
            assert pq.read_table(p, partitioning=None).to_pylist() == rows[a:b], p.name
    assert res["files"] == len(layout.cut([60] * len(want_dict), cap)) + len(want_tok)
    man = catalog.build_manifest(fresh, REL, ["S0"])
    assert contracts.check_sort_order(man) == [] and contracts.check_statistics(man) == []


@pytest.mark.covers("D-01", "D-04", "LY-08")
def test_b09_token_bytes_do_not_depend_on_the_run_or_on_input_layout(built, tmp_path):
    root, _, _, circ = built
    want = _digests(root, TOKEN_DIRS)
    rng = random.Random(9)
    for name in ("again", "shuffled", "split"):
        fresh = _fresh(root, tmp_path, name, TOKEN_DIRS)
        if name == "shuffled":  # domain and cluster rows in another order, in several row groups
            for f in sorted(fresh.glob("l1/domains/**/*.parquet")) + sorted(fresh.glob("l2/cluster_member/**/*.parquet")):
                t = pq.read_table(f, partitioning=None)
                idx = list(range(t.num_rows))
                rng.shuffle(idx)
                pq.write_table(t.take(idx), f, row_group_size=max(1, t.num_rows // 5))
        if name == "split":  # the tapes in three parts
            for d in TAPE_DIRS:
                shutil.rmtree(fresh / d)
            _drop_tape_ledgers(fresh)
            n = sorted(t["n_slots"] for t in _tapes_of_tv(root).values())
            assert tapes.run(fresh, "S0", circ, max_file_bytes=200 + 30 * n[0])["files"] == 3
        tokens.run(fresh, REL, ["S0"])
        assert _digests(fresh, TOKEN_DIRS) == want, name


@pytest.mark.covers("C-01", "LY-05", "LY-02", "LY-03")
def test_b09_token_files_conform_and_the_manifest_accepts_them(built):
    root, man, con, _ = built
    by = defaultdict(list)
    for f in man["files"]:
        by[f["table"]].append(f)
    assert [f["partitions"] for f in by["token_dict"]] == [{"release": REL}]
    for name in ("token_dict", "tape_tokens"):
        for f in by[name]:
            assert conformance_problems(pq.read_table(root / f["path"], partitioning=None), table(name)) == [], f["path"]
    assert contracts.check_sort_order(man) == [] and contracts.check_statistics(man) == []
    assert contracts.check_list_lengths(con) == [] and contracts.check_foreign_keys(con) == []
    counted = con.execute("SELECT level, sum(n_slots) FROM token_dict GROUP BY 1 ORDER BY 1").fetchall()
    on_slots = [(lvl, con.execute(f"SELECT count(t) FROM (SELECT unnest({lvl}_tokens) AS t FROM tape_tokens)").fetchall()[0][0])
                for lvl in sorted(LEVELS)]
    assert counted == on_slots


def test_b09_tokens_never_load_every_tape(built, tmp_path, monkeypatch):
    root, _, _, _ = built
    fresh = _fresh(root, tmp_path, "s", TOKEN_DIRS)

    def refuse(*a, **k):
        raise AssertionError("tokens.run read every tape into Python")

    monkeypatch.setattr(tokens, "read_tapes", refuse)
    parent = tmp_path / "scratch-parent"
    parent.mkdir()
    tokens.run(fresh, REL, ["S0"], tmp=parent)
    assert _digests(fresh, TOKEN_DIRS) == _digests(root, TOKEN_DIRS)
    assert list(parent.iterdir()) == []  # the private scratch directory is removed


def test_b09_refusals(built, tmp_path):
    root, _, _, _ = built
    fresh = _fresh(root, tmp_path, "s", TOKEN_DIRS)
    with pytest.raises(ValueError, match="31 libraries"):
        tokens.run(fresh, REL, ["S0"], libraries=[f"lib{i}-1" for i in range(32)])
    with pytest.raises(ValueError, match="scratch"):
        tokens.run(fresh, REL, ["S0"], tmp=fresh / "tmp")
    assert not (fresh / "l3/token_dict").exists() and not (fresh / "l3/tape_tokens").exists()


def test_b09_iter_tapes_yields_rows_in_contig_order_across_files_and_row_groups(built, tmp_path):
    root, _, _, _ = built
    rows = sorted(_tapes_of_tv(root).values(), key=lambda r: r["contig_uid"])
    rows = [{k: v for k, v in r.items() if k in table("contig_tapes").schema.names} for r in rows]
    fake = tmp_path / "r"
    spec = table("contig_tapes")
    for stage, part in (("S0", rows[0::2]), ("S1", rows[1::2])):  # interleaved in contig order across stages
        io.write_table([pa.Table.from_pylist(part, schema=spec.schema)], spec,
                       fake / f"l3/tapes/tape_version={TV}/stage={stage}/part-0.parquet", row_group_rows=1)
    it = tokens.iter_tapes(fake, TV, ["S0", "S1"])
    assert inspect.isgenerator(it)
    assert list(it) == rows
    assert list(tokens.iter_tapes(fake, TV, ["S0", "S1"], ["n_slots"])) == [{"n_slots": r["n_slots"]} for r in rows]
    # card B-25: only the stages given are read, and a stage given without tapes of the version is an error
    assert list(tokens.iter_tapes(fake, TV, ["S1"])) == rows[1::2]
    with pytest.raises(FileNotFoundError, match="S0"):
        list(tokens.iter_tapes(fake, "t-000000000000", ["S0"]))
    assert [r["contig_uid"] for r in tokens.read_tapes(root, V, ["contig_uid"], stages=["S0"])] == [
        r["contig_uid"] for r in rows]


# --- card B-11: pair macros over token postings ------------------------------------------------------

SYN_PRODUCER = "synthetic-b11"
SYN_REL = "rS0.b11syn"
FAM, NR50, NR90 = 1 << 61 | 1 << 56, 2 << 61, 3 << 61  # token id bases (PLAN.md §4.1): level | source


def _syn_slots(rng: random.Random, n: int, length: int, circular: bool, wrap: bool) -> list[tuple[int, int, int]]:
    """n slots in coordinate order with overlaps and slots nested in the previous one on the other strand; with
    `wrap`, the last slot of a circular contig crosses the origin."""
    step, slots = length // n, []
    for idx in range(n):
        prev = slots[-1] if slots else None
        if prev and prev[1] - prev[0] > 60 and rng.random() < 0.25:
            s, e, strand = prev[0] + 10, prev[1] - 10, -prev[2]
        else:
            s = step * idx + rng.randrange(0, step // 3)
            e, strand = s + rng.randrange(20, 2 * step), rng.choice((1, -1))
        slots.append((s, min(e, length), strand))
    slots.sort()
    if wrap:
        s, _, strand = slots[-1]
        slots[-1] = (s, s + 3 * step, strand)
        assert slots[-1][1] > length
    return slots


def _write_syn_tapes(root: Path, contigs: list[tuple[str, list[tuple[int, int, int]], int, bool]], seed: int) -> str:
    """contig_tapes and tape_tokens for the synthetic tape version, written directly; returns the tape version."""
    rng = random.Random(seed)
    tv = tapes.tape_version([SYN_PRODUCER])
    tape_rows, tok_rows = [], []
    for name, slots, length, circular in contigs:
        cu, n = ids.contig_uid(name.encode()), len(slots)
        tape_rows.append({"contig_uid": cu, "contig_len": length, "is_circular": circular, "n_slots": n,
                          **tapes.build_tape(slots, length, circular), "protein_uids": [None] * n,
                          "feature_uids": [ids.feature_uid(cu, s, e, st, ids.KIND["cds"]) for s, e, st in slots]})
        pick = lambda base, k, p_none: [None if rng.random() < p_none else base | rng.randrange(k) for _ in slots]  # noqa: E731
        tok_rows.append({"contig_uid": cu, "family_tokens": pick(FAM, 3, 0.3), "nr50_tokens": pick(NR50, 3, 0.2),
                         "nr90_tokens": pick(NR90, 4, 0.1)})
    for name, rows, path in (("contig_tapes", tape_rows, f"l3/tapes/tape_version={tv}/stage=S0/part-0.parquet"),
                             ("tape_tokens", tok_rows, f"l3/tape_tokens/tape_version={tv}/release={SYN_REL}/part-0.parquet")):
        spec = table(name)
        io.write_table([pa.Table.from_pylist(sorted(rows, key=lambda r: r["contig_uid"]), schema=spec.schema)], spec,
                       root / path)
    return tv


@pytest.mark.covers("Q1", "Q4")
def test_b11_pairs_of_equals_the_oracle_on_even_circles_and_short_lines(tmp_path):
    """Synthetic tapes through `postings.run`: circles with an even slot count of at most 20 (a partner exactly half
    way round is within 10 slots), a 2-slot circle whose two ways round have equal gaps (forward on a tie), slots
    across the origin, nested antisense slots, and linear contigs shorter than 21 slots (every anchor censored).
    Every ordered token pair against `pairs_of_contig`."""
    rng = random.Random(1111)
    contigs = [("even-8", _syn_slots(rng, 8, 4_000, True, False), 4_000, True),
               ("even-20-wrap", _syn_slots(rng, 20, 9_000, True, True), 9_000, True),
               ("even-12-wrap", _syn_slots(rng, 12, 6_000, True, True), 6_000, True),
               ("even-2", _syn_slots(rng, 2, 1_500, True, False), 1_500, True),
               ("tie-2", [(0, 100, 1), (500, 600, -1)], 1_000, True),
               ("odd-9-wrap", _syn_slots(rng, 9, 5_000, True, True), 5_000, True),
               ("linear-15", _syn_slots(rng, 15, 7_000, False, False), 7_000, False),
               ("linear-1", [(10, 300, -1)], 400, False),
               ("circle-1", [(10, 300, 1)], 400, True)]
    root = tmp_path / "r"
    tv = _write_syn_tapes(root, contigs, 1112)
    postings.run(root, SYN_REL, ["S0"], producers=[SYN_PRODUCER])
    man = catalog.build_manifest(root, SYN_REL, ["S0"])
    con = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / "c.duckdb"))
    try:
        assert con.execute("SELECT catalog_tape_version()").fetchall() == [(tv,)]
        rows = _oracle_rows(root, SYN_REL, tv)
        toks = sorted({t for _, k in rows for lvl in LEVELS for t in k[f"{lvl}_tokens"] if t is not None})
        assert len(toks) == 10
        same_level = [(a, b) for a in toks for b in toks if a >> 61 == b >> 61]  # the oracle pairs within a level
        seen = _check_pairs_against_the_oracle(con, rows, same_level)
    finally:
        con.close()
    assert seen["crosses_origin"] and seen["censored"] and seen["topology_4"], seen
    n_of = {t["contig_uid"]: (t["n_slots"], t["is_circular"]) for t, _ in rows}
    half = [r for rs in _oracle_pairs(rows).values() for r in rs if n_of[r[0]][1] and 2 * abs(r[5]) == n_of[r[0]][0]]
    assert any(r[5] < 0 for r in half) and any(r[5] > 0 for r in half)  # the half-way rule takes both ways
    tie = ids.contig_uid(b"tie-2")
    assert {r[5] for r in half if r[0] == tie} == {1}  # equal gaps both ways: forward
    lin = ids.contig_uid(b"linear-15")
    assert all(r[8] & FLAG_CENSORED for rs in _oracle_pairs(rows).values() for r in rs if r[0] == lin)


@pytest.mark.covers("Q1")
def test_b11_q1_count_filters_tiers_through_contig_dim(built, tmp_path):
    """The tier filter joins contig_dim of the catalog's release: the builder's null tiers match no tier; tiers set
    by hand in a second release select their contigs, and a contig without a dimension row matches none."""
    root, _, _, circ = built
    fresh = tmp_path / "s"
    shutil.copytree(root, fresh)
    cus = sorted(_tapes_of_tv(root))
    write_genome_input(fresh, "S0", [GenomeInput((cu,), circular=(cu,) if cu in circ else ()) for cu in cus])
    dim.run(fresh, REL, ["S0"])
    dim_rows = _rows(fresh, f"dim/contig_dim/release={REL}/part-*.parquet")
    assert [r["contig_uid"] for r in dim_rows] == cus and {r["tier"] for r in dim_rows} == {None}
    rows = _oracle_rows(root)
    full = _oracle_pairs(rows)
    pairs = sorted(full)
    lists = {"a": [a for a, _ in pairs], "b": [b for _, b in pairs]}
    by_tier = f"SELECT l.a, l.b, q1_count(l.a, l.b, 10, NULL), q1_count(l.a, l.b, 10, 1), q1_count(l.a, l.b, 10, 2), " \
              f"q1_count(l.a, l.b, 10, 3) FROM {_PAIR_LIST}"
    con = catalog.open_catalog(catalog.generate_catalog(catalog.build_manifest(fresh, REL, ["S0"]), tmp_path / "a.duckdb"))
    try:
        got = con.execute(by_tier, lists).fetchall()
        assert sorted(got) == [(a, b, _loci(full[(a, b)]), 0, 0, 0) for a, b in pairs]  # null tiers match no tier
        a, b = pairs[0]
        assert con.execute("SELECT q1_count(?, ?, 10, 1)", [a, b]).fetchall() == [(0,)]
    finally:
        con.close()
    rel2 = "rS0.b11tier"
    tiers = {cus[0]: 1, cus[1]: 2}  # cus[2] has no dimension row in rel2
    _write_sorted("contig_dim", [{**{k: v for k, v in r.items() if k != "release"}, "tier": tiers[r["contig_uid"]]}
                                 for r in dim_rows if r["contig_uid"] in tiers],
                  fresh / f"dim/contig_dim/release={rel2}/part-0.parquet")
    shutil.copytree(fresh / f"l3/token_postings/release={REL}", fresh / f"l3/token_postings/release={rel2}")
    con = catalog.open_catalog(catalog.generate_catalog(catalog.build_manifest(fresh, rel2, ["S0"]), tmp_path / "b.duckdb"))
    try:
        want = [(a, b, _loci(full[(a, b)]), *(_loci([r for r in full[(a, b)] if tiers.get(r[0]) == t]) for t in (1, 2, 3)))
                for a, b in pairs]
        assert sorted(con.execute(by_tier, lists).fetchall()) == want
        for a, b, _, t1, t2, _ in want[:20]:
            assert con.execute("SELECT q1_count(?, ?, 10, 1), q1_count(?, ?, 10, 2)", [a, b, a, b]).fetchall() == [(t1, t2)]
    finally:
        con.close()
    assert any(w[3] for w in want) and any(w[4] for w in want)
    assert any(r[0] == cus[2] for rs in full.values() for r in rs)  # the contig without a row had loci to drop


def _write_postings(root: Path, release: str, tv: str, rows: list[dict]) -> None:
    """A hand-written token_postings file (bucket `all`) with the registry's defaults for the columns not given."""
    full = [{"seg_idx": 0, "start": 100 * r["ordinal"], "end": 100 * r["ordinal"] + 90, "strand": 1, "is_circular": False,
             "feature_uid": ids.feature_uid(r["contig_uid"], 100 * r["ordinal"], 100 * r["ordinal"] + 90, 1, 1),
             "contig_len": 100 * r["n_slots"], **r} for r in rows]
    _write_sorted("token_postings", full, root / f"l3/token_postings/release={release}/tape_version={tv}/"
                                              "token_bucket=all/part-0.parquet")


@pytest.mark.covers("Q4")
def test_b11_fusions_of_finds_two_tokens_on_one_slot_with_different_segments(built, tmp_path):
    root, _, con, _ = built
    # the fixture's postings carry seg_idx 0 only, so no fusion exists yet (card B-21 adds per-domain postings)
    full = _oracle_pairs(_oracle_rows(root))
    for a, b in random.Random(14).sample(sorted(full), 20):
        assert con.execute("SELECT * FROM fusions_of(?, ?)", [a, b]).fetchall() == []
    tv = tapes.tape_version([SYN_PRODUCER])
    x, y = ids.contig_uid(b"fusion-x"), ids.contig_uid(b"fusion-y")
    f1, f2 = FAM | 1, FAM | 2
    rows = [{"token_id": f1, "contig_uid": x, "ordinal": 2, "seg_idx": 0, "n_slots": 10},
            {"token_id": f2, "contig_uid": x, "ordinal": 2, "seg_idx": 1, "n_slots": 10},  # a fusion at slot 2
            {"token_id": f1, "contig_uid": x, "ordinal": 5, "seg_idx": 0, "n_slots": 10},
            {"token_id": f2, "contig_uid": x, "ordinal": 6, "seg_idx": 1, "n_slots": 10},  # neighbours, no fusion
            {"token_id": f1, "contig_uid": y, "ordinal": 7, "seg_idx": 0, "n_slots": 10},  # f1 twice on slot 7
            {"token_id": f1, "contig_uid": y, "ordinal": 7, "seg_idx": 2, "n_slots": 10},
            {"token_id": f2, "contig_uid": y, "ordinal": 7, "seg_idx": 1, "n_slots": 10}]
    r = tmp_path / "r"
    _write_postings(r, SYN_REL, tv, rows)
    c = catalog.open_catalog(catalog.generate_catalog(catalog.build_manifest(r, SYN_REL, []), tmp_path / "c.duckdb"))
    fu = lambda cu, i: ids.feature_uid(cu, 100 * i, 100 * i + 90, 1, 1)  # noqa: E731
    try:
        got = c.execute("SELECT contig_uid, feature_uid, ordinal, anchor_seg_idx, partner_seg_idx "
                        "FROM fusions_of(?, ?)", [f1, f2]).fetchall()
        want = sorted([(x, fu(x, 2), 2, 0, 1), (y, fu(y, 7), 7, 0, 1), (y, fu(y, 7), 7, 2, 1)])
        assert got == want
        assert c.execute("SELECT count(*) FROM fusions_of(?, ?)", [f1, f1]).fetchall()[0][0] == 2  # seg 0 and seg 2
        # a slot with a token on two segments is one slot for the pair macros
        assert c.execute("SELECT anchor_ordinal, partner_ordinal FROM pairs_of(?, ?, 10)", [f1, f2]).fetchall() == \
            [(2, 6), (5, 2), (5, 6)]
        assert c.execute("SELECT q1_count(?, ?, 10, NULL)", [f1, f2]).fetchall()[0][0] == 2
    finally:
        c.close()


def test_b11_the_catalog_serves_one_tape_version(tmp_path):
    x = ids.contig_uid(b"tv-x")
    tv1, tv2 = tapes.tape_version(["p-one"]), tapes.tape_version(["p-two"])
    root = tmp_path / "r"
    for tv, (i, j) in ((tv1, (0, 1)), (tv2, (0, 3))):
        _write_postings(root, SYN_REL, tv, [{"token_id": NR90 | 1, "contig_uid": x, "ordinal": i, "n_slots": 9},
                                            {"token_id": NR90 | 2, "contig_uid": x, "ordinal": j, "n_slots": 9}])
    man = catalog.build_manifest(root, SYN_REL, [])  # postings only: no stage file
    with pytest.raises(ValueError, match="tape version"):
        catalog.generate_catalog(man, tmp_path / "a.duckdb")
    for tv, d in ((tv1, 1), (tv2, 3)):
        c = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / f"{tv}.duckdb", tape_version=tv))
        try:
            assert c.execute("SELECT catalog_tape_version()").fetchall() == [(tv,)]
            assert c.execute("SELECT orf_offset FROM pairs_of(?, ?, 10)", [NR90 | 1, NR90 | 2]).fetchall() == [(d,)]
        finally:
            c.close()
    # postings of another release only: no tape version is served and the pair macros return nothing
    c = catalog.open_catalog(catalog.generate_catalog(catalog.build_manifest(root, "rOther", []), tmp_path / "o.duckdb",
                                                      tape_version=None))
    try:
        assert c.execute("SELECT catalog_tape_version()").fetchall() == [("",)]
        assert c.execute("SELECT q1_count(?, ?, 10, NULL)", [NR90 | 1, NR90 | 2]).fetchall() == [(0,)]
        assert c.execute("SELECT count(*) FROM pairs_of(?, ?, 10)", [NR90 | 1, NR90 | 2]).fetchall() == [(0,)]
        assert c.execute("SELECT count(*) FROM partners_of(?, 10)", [NR90 | 1]).fetchall() == [(0,)]
    finally:
        c.close()


@st.composite
def _contig_slots(draw):
    circ = draw(st.booleans())
    length = draw(st.integers(min_value=200, max_value=5_000))
    n = draw(st.integers(min_value=2, max_value=24))
    starts = sorted(draw(st.lists(st.integers(min_value=0, max_value=length - 2), min_size=n, max_size=n, unique=True)))
    slots = []
    for idx, s in enumerate(starts):
        e = s + draw(st.integers(min_value=1, max_value=600))
        # only the last slot of a circle may cross the origin (tapes order slots by start)
        slots.append((s, min(e, s + length) if circ and idx == n - 1 else min(e, length), draw(st.sampled_from([1, -1]))))
    return slots, length, circ


HELPERS_SQL = """
    SELECT slot_dist(i, j, n, circ), pair_offset(i, a_s, a_e, j, b_s, b_e, n, circ, len),
           pair_delta_bp(d, i, a_s, a_e, j, b_s, b_e, circ, len),
           pair_junction(d, i, a_s, a_e, a_str, j, b_s, b_e, b_str, circ, len), pair_flags(d, i, j, n, circ)
    FROM (SELECT unnest($i)::UINTEGER AS i, unnest($j)::UINTEGER AS j, unnest($d)::BIGINT AS d,
                 unnest($s1)::BIGINT AS a_s, unnest($e1)::BIGINT AS a_e, unnest($t1)::TINYINT AS a_str,
                 unnest($s2)::BIGINT AS b_s, unnest($e2)::BIGINT AS b_e, unnest($t2)::TINYINT AS b_str,
                 $n::UINTEGER AS n, $circ::BOOLEAN AS circ, $clen::BIGINT AS len, unnest(range(len($i))) AS rn)
    ORDER BY rn
"""


@pytest.mark.covers("Q4")
@settings(max_examples=150, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(_contig_slots())
def test_b11_pair_helper_macros_equal_geometry_and_the_oracle(built, c):
    """slot_dist, pair_offset, pair_delta_bp, pair_junction and pair_flags on every ordered slot pair of a contig,
    against biobox.l3.geometry for both signs of the offset, and against `pairs_of_contig` for the rows it gives."""
    _, _, con, _ = built
    slots, length, circ = c
    n = len(slots)
    cases = []
    for i in range(n):
        for j in range(n):
            if i != j:
                off = geometry.cyclic_offset(i, j, n) if circ else j - i
                cases += [(i, j, off), (i, j, -off)]
    params = {"i": [x[0] for x in cases], "j": [x[1] for x in cases], "d": [x[2] for x in cases],
              "s1": [slots[x[0]][0] for x in cases], "e1": [slots[x[0]][1] for x in cases],
              "t1": [slots[x[0]][2] for x in cases], "s2": [slots[x[1]][0] for x in cases],
              "e2": [slots[x[1]][1] for x in cases], "t2": [slots[x[1]][2] for x in cases],
              "n": n, "circ": circ, "clen": length}
    got = con.execute(HELPERS_SQL, params).fetchall()
    assert len(got) == len(cases)
    chosen = {}
    for (i, j, d), (dist, offset, delta, topo, flags) in zip(cases, got):
        assert dist == (min(abs(i - j), n - abs(i - j)) if circ else abs(i - j))
        want_d = geometry.cyclic_offset(i, j, n) if circ else j - i
        if circ and 2 * abs(want_d) == n and \
                geometry.pair_geometry(slots, i, j, -want_d, length, circ)[1] < geometry.pair_geometry(slots, i, j, want_d, length, circ)[1]:
            want_d = -want_d
        assert offset == want_d, (i, j)
        g_topo, g_delta, g_cross = geometry.pair_geometry(slots, i, j, d, length, circ)
        assert (topo, delta) == (g_topo, g_delta), (i, j, d)
        censored = not circ and (i < K_MAX or n - 1 - i < K_MAX)
        assert flags == (FLAG_CENSORED if censored else 0) | (FLAG_CROSSES_ORIGIN if g_cross else 0), (i, j, d)
        if d == offset:
            chosen[(i, j)] = (dist, offset, topo, delta, flags)
    tape = {"contig_uid": b"h" * 16, "is_circular": circ, "contig_len": length, **tapes.build_tape(slots, length, circ)}
    oracle = pairs_of_contig(tape, {"family": [None] * n, "nr50": [None] * n, "nr90": [1] * n}, K_MAX)
    assert {(r["anchor_ordinal"], r["partner_ordinal"]) for r in oracle} == \
        {p for p, v in chosen.items() if 1 <= v[0] <= K_MAX}
    for r in oracle:
        _, offset, topo, delta, flags = chosen[(r["anchor_ordinal"], r["partner_ordinal"])]
        assert (offset, topo, delta, flags) == (r["orf_offset"], r["junction_topology"], r["delta_bp"], r["flags"])
