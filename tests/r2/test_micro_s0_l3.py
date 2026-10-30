"""Ring R2: L3 tapes, tokens and pairs on micro-S0, and the Q1-Q4, Q7, Q8 oracles (data_layer_and_tests.md §4.3).

The build: ingest, gene calls, a 6-model domain library from lambda proteins (family
tokens) and NR50/NR90 clusters from MMseqs2 where it is on PATH, else a deterministic
stand-in grouping (the oracles read the same tokens, so what they test does not depend on
how the clusters were made). Each production path (catalog macros) is compared with a
brute-force oracle written here from the definitions: a nested loop over tape rows read
straight from the Parquet shards, and a rebuild of the tape lists from protein_occurrence.
Covers Q1, Q2, Q3 (partner counts; z_score ranking needs pair_stats), Q4, Q7, Q8, C-04,
C-07, K-07, LY-06 for pair lookups, and D-01 for the L3 tables.
"""

from __future__ import annotations

import random
import shutil
from collections import defaultdict
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pyhmmer
import pytest
import yaml

from biobox import catalog, contracts, fasta, ids
from biobox.l0seq import SourceObject, ingest, revcomp
from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION as V, GenomeInput, run as call_genes
from biobox.l3 import pairs, tapes, tokens
from biobox.registry.schemas import table

pytestmark = pytest.mark.r2

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
MMSEQS = shutil.which("mmseqs")
SAMTOOLS = shutil.which("samtools")
REL = "rS0.l3"
LEVELS = ("family", "nr50", "nr90")


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


def build(root: Path, order: int = 1) -> set[bytes]:
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
        cluster.run(root, REL, V, tmp=tmp)
    else:
        _stand_in_clusters(root, REL)
    tapes.run(root, "S0", circ)
    tokens.run(root, REL)
    pairs.run(root, REL)
    return circ


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    d = tmp_path_factory.mktemp("l3")
    circ = build(d / "a")
    man = catalog.build_manifest(d / "a", REL)
    con = catalog.open_catalog(catalog.generate_catalog(man, d / "c.duckdb"))
    yield d / "a", man, con, circ
    con.close()


# --- oracles -------------------------------------------------------------------------------


def _read(root: Path, glob: str) -> list[dict]:
    return [r for f in sorted(root.glob(glob)) for r in pq.read_table(f).to_pylist()]


def _oracle_rows(root: Path) -> list[tuple[dict, dict]]:
    """(tape, tokens) per contig, straight from the shards."""
    t = {r["contig_uid"]: r for r in _read(root, f"l3/tapes/gene_call_version={V}/stage=*/part-*.parquet")}
    k = {r["contig_uid"]: r for r in _read(root, f"l3/tape_tokens/gene_call_version={V}/release={REL}/part-*.parquet")}
    return [(t[cu], k[cu]) for cu in sorted(t)]


def _within(i: int, j: int, n: int, circular: bool) -> int:
    if not circular:
        return abs(j - i)
    return min((j - i) % n, (i - j) % n)


def _oracle_instances(rows, a: int, b: int, k: int, tier=None) -> list[tuple]:
    """Every (contig, anchor slot, partner slot) with token a at the anchor and b within k slots."""
    out = []
    for tape, tok in rows:
        if tier is not None and tape["tier"] != tier:
            continue
        n = tape["n_slots"]
        for lvl in LEVELS:
            col = tok[f"{lvl}_tokens"]
            for i in range(n):
                if col[i] != a:
                    continue
                for j in range(n):
                    if j != i and col[j] == b and _within(i, j, n, tape["is_circular"]) <= k:
                        out.append((tape["contig_uid"], i, j))
    return sorted(out)


def _sample_tokens(con, n: int, seed: int) -> list[int]:
    """Stratified by document frequency decile, so the most frequent tokens are always sampled."""
    rows = con.execute("SELECT token_id, n_slots FROM token_dict WHERE release = ? ORDER BY n_slots, token_id",
                       [REL]).fetchall()
    rng = random.Random(seed)
    deciles = [rows[i * len(rows) // 10:(i + 1) * len(rows) // 10] for i in range(10)]
    picked = [rng.choice(d)[0] for d in deciles if d] + [rows[-1][0]]
    return sorted(set(picked))[:n]


# --- Q1-Q4 ---------------------------------------------------------------------------------


@pytest.mark.covers("Q1")
@pytest.mark.parametrize("k", [1, 3, 10])
def test_q1_count_equals_nested_loop(built, k):
    root, _, con, _ = built
    rows = _oracle_rows(root)
    anchors = _sample_tokens(con, 12, k)
    checked = nonzero = 0
    for a in anchors:
        partners = {b for _, b in con.execute("SELECT DISTINCT 0, partner_token FROM slot_pairs WHERE anchor_token = ?",
                                              [a]).fetchall()}
        for b in sorted(partners)[:6] + [anchors[0]]:
            for tier in (None, 1):
                want = len({(c, i) for c, i, _ in _oracle_instances(rows, a, b, k, tier)})
                got = con.execute("SELECT q1_count(?, ?, ?, ?)", [a, b, k, tier]).fetchall()[0][0]
                assert got == want, (a, b, k, tier)
                checked += 1
                nonzero += want > 0
    assert checked > 50 and nonzero > 10


@pytest.mark.covers("Q3")
def test_q3_partner_profile_equals_nested_loop(built):
    root, _, con, _ = built
    rows = _oracle_rows(root)
    for a in _sample_tokens(con, 10, 3):
        want: dict[int, set] = defaultdict(set)
        n_pairs: dict[int, int] = defaultdict(int)
        for tape, tok in rows:
            n = tape["n_slots"]
            for lvl in LEVELS:
                col = tok[f"{lvl}_tokens"]
                for i in range(n):
                    if col[i] != a:
                        continue
                    for j in range(n):
                        if j != i and col[j] is not None and _within(i, j, n, tape["is_circular"]) <= 10:
                            want[col[j]].add((tape["contig_uid"], i))
                            n_pairs[col[j]] += 1
        got = con.execute("SELECT partner_token, n_loci, n_pairs FROM partners_of(?, 10)", [a]).fetchall()
        assert sorted(got) == sorted((b, len(s), n_pairs[b]) for b, s in want.items()), a
        assert [g[1] for g in got] == sorted((g[1] for g in got), reverse=True)  # ranked


@pytest.mark.covers("Q4")
def test_q4_drill_down_equals_nested_loop(built):
    root, _, con, _ = built
    rows = _oracle_rows(root)
    by_contig = {t["contig_uid"]: t for t, _ in rows}
    for a in _sample_tokens(con, 8, 4):
        for (b,) in con.execute("SELECT DISTINCT partner_token FROM slot_pairs WHERE anchor_token = ? ORDER BY 1 LIMIT 4",
                                [a]).fetchall():
            got = con.execute("SELECT contig_uid, anchor_ordinal, partner_ordinal, orf_offset, junction_topology, "
                              "delta_bp FROM pairs_of(?, ?, 10)", [a, b]).fetchall()
            assert sorted((c, i, j) for c, i, j, *_ in got) == _oracle_instances(rows, a, b, 10)
            for c, i, j, d, topo, delta in got:
                t = by_contig[c]
                s = list(zip(t["start_bps"], t["end_bps"], t["strands"]))
                left, right = (s[i], s[j]) if d > 0 else (s[j], s[i])
                if t["is_circular"] and ((d > 0) != (j > i)):
                    right = (right[0] + t["contig_len"], right[1] + t["contig_len"], right[2])
                assert delta == max(left[0], right[0]) - min(left[1], right[1])
                assert (topo == 1) == (left[2] == right[2])


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
        # independent rebuild from L1: slots in ordinal order, protein UIDs, strands, coordinates
        slots = sorted(occ[g["contig_uid"]], key=lambda r: r["slot_ordinal"])
        assert g["n_slots"] == len(slots)
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
    "SELECT junction_topology, count(*), sum(delta_bp) FROM slot_pairs GROUP BY 1 ORDER BY 1",
    "SELECT is_circular, count(*), sum(n_slots), max(contig_len) FROM contig_tapes GROUP BY 1 ORDER BY 1",
    "SELECT tier, count(DISTINCT contig_uid), count(*) FROM slot_pairs WHERE abs(orf_offset) <= 2 GROUP BY 1 ORDER BY 1",
    "SELECT strand, count(*), sum(\"end\" - start) FROM protein_occurrence GROUP BY 1 ORDER BY 1",
    "SELECT model_acc, count(DISTINCT protein_uid) FROM l1_domains GROUP BY 1 ORDER BY 1",
    "SELECT d.level, count(*) FROM slot_pairs p JOIN token_dict d ON d.token_id = p.anchor_token "
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
    t = con.sql("SELECT * EXCLUDE (gene_call_version, stage) FROM contig_tapes").to_arrow_table().to_pylist()
    short = [dict(r) for r in t]
    short[0]["strands"] = short[0]["strands"][:-1]
    assert any("parallel lists" in p for p in contracts.check_list_lengths(_with_gcv(con, short)))
    flagged = [dict(r) for r in t]
    lin = next(r for r in flagged if r["orfs_to_contig_starts"] and r["n_slots"])
    lin["slot_flags"] = [1] * lin["n_slots"]  # flagged as clipped but not saturated
    assert any("not saturated" in p for p in contracts.check_saturation(_with_gcv(con, flagged)))


def _with_gcv(con, tape_rows):
    c = duckdb.connect()
    tbl = pa.Table.from_pylist(tape_rows)
    c.register("_t", tbl)
    c.execute(f"CREATE VIEW contig_tapes AS SELECT *, '{V}' AS gene_call_version FROM _t")
    c.register("tape_tokens", con.sql("SELECT * FROM tape_tokens").to_arrow_table())
    return c


@pytest.mark.covers("LY-06", "LY-05")
@pytest.mark.parametrize("width", [1, 2])
def test_ly06_pair_lookup_opens_one_bucket(built, tmp_path, width):
    root, _, _, _ = built
    copy = tmp_path / "s"
    shutil.copytree(root, copy, ignore=shutil.ignore_patterns("slot_pairs"))
    pairs.run(copy, REL, prefix_width=width)
    man = catalog.build_manifest(copy, REL)  # LY-05: every row's token matches its bucket
    files = [f for f in man["files"] if f["table"] == "slot_pairs"]
    assert len(files) > 1
    c = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / "c.duckdb"))
    try:
        a, b = c.execute("SELECT anchor_token, partner_token FROM slot_pairs ORDER BY anchor_token LIMIT 1").fetchall()[0]
        prof = catalog.profile_query(c, "SELECT q1_count(?, ?, 10, NULL)", [a, b])
        (scan,) = [s for s in prof["scans"] if any("/slot_pairs/" in f for f in s["files"])]
        assert prof["result"][0][0] >= 1
        assert scan["files_read"] == 1 and scan["row_groups_scanned"] == 1
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
    h2 = contracts.content_hashes(catalog.build_manifest(other, REL))
    for name in ("contig_tapes", "tape_tokens", "token_dict", "slot_pairs"):
        assert h1[name] == h2[name], name
    for rel in sorted(p.relative_to(root).as_posix() for p in (root / "l3").rglob("*.parquet")):
        assert (root / rel).read_bytes() == (other / rel).read_bytes(), rel
