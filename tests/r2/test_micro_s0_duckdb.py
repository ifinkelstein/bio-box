"""Ring R2: DuckDB behavior (K tests, data_layer_and_tests.md §4.7) and LY-06 on micro-S0.

These are the micro-scale forms. The pod forms (K-02 at 150 GB inside the 200 GB scope,
K-04 CPU tracking at 8/20/40 threads, K-10 NFS vs local) run in ring R4 on wilkcomp02.

LY-06 counts storage units, not time: row groups scanned and files read per canonical
query, from DuckDB's profiler (`catalog.profile_query`), against the count expected from
the Parquet min/max statistics.
"""

from __future__ import annotations

import hashlib
import subprocess
import sys
import threading
import time
import tomllib
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path

import duckdb
import pyarrow.parquet as pq
import pytest
import yaml

from biobox import catalog, fasta, ids
from biobox.l0seq import SourceObject, ingest
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, run as call_genes

pytestmark = pytest.mark.r2

REPO = Path(__file__).resolve().parents[2]
FIX = REPO / "tests" / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
CHUNK, RG_CHUNKS = 4096, 16


def _obj(o):
    return SourceObject(FIX / o["file"], CONTRACT["source_name"], CONTRACT["source_release"], CONTRACT["license"], True)


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """lambda and M. genitalium in S0, T4 in S1, so stage pruning has something to prune."""
    root = tmp_path_factory.mktemp("duck")
    objs = CONTRACT["objects"]
    s0 = [o for o in objs if o["accession"] != "NC_000866.4"]
    s1 = [o for o in objs if o["accession"] == "NC_000866.4"]
    ingest([_obj(o) for o in s0], root, "S0", chunk_size=CHUNK, row_group_chunks=RG_CHUNKS)
    ingest([_obj(o) for o in s1], root, "S1", chunk_size=CHUNK, row_group_chunks=RG_CHUNKS)
    for stage, group in (("S0", s0), ("S1", s1)):
        genomes = []
        for o in group:
            cu = ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)
            genomes.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"],
                                       circular=(cu,) if o["topology"] == "circular" else ()))
        call_genes(root, stage, genomes)
    man = catalog.build_manifest(root, "rS1.test")
    cat = catalog.generate_catalog(man, root / "catalog/release=rS1.test.duckdb")
    return root, man, cat


@pytest.fixture()
def con(built):
    c = catalog.open_catalog(built[2])
    yield c
    c.close()


def _seq_scan(prof: dict, table_dir: str) -> dict:
    scans = [s for s in prof["scans"] if any(table_dir in f for f in s["files"])]
    assert len(scans) == 1, scans
    return scans[0]


def _expected_row_groups(root: Path, man: dict, rd: bytes, first_ix: int, last_ix: int) -> int:
    """Row groups whose (raw_digest, chunk_ix) min/max range overlaps the wanted chunks."""
    n = 0
    for f in man["files"]:
        if f["table"] != "seq_chunks" or f["partitions"]["digest_prefix"] != rd[:1].hex():
            continue
        md = pq.ParquetFile(root / f["path"]).metadata
        for rg in range(md.num_row_groups):
            d, c = md.row_group(rg).column(0).statistics, md.row_group(rg).column(1).statistics
            lo, hi = (bytes(d.min), c.min), (bytes(d.max), c.max)
            if lo <= (rd, last_ix) and (rd, first_ix) <= hi:
                n += 1
    return n


# --- LY-06 / K-06: canonical point queries --------------------------------------------


@pytest.mark.parametrize("start,end", [
    (300_000, 301_000),                                   # inside one chunk
    (CHUNK * RG_CHUNKS * 3 - 10, CHUNK * RG_CHUNKS * 3 + 10),  # across a row-group boundary
    (0, 1),
    (580_000, 580_076),                                   # last chunk
])
def test_ly06_slice_reads_one_file_and_only_the_row_groups_it_needs(built, con, start, end):
    root, man, _ = built
    seq = fasta.read(FIX / "NC_000908.2.fna.gz")[0].seq
    cu, rd = ids.contig_uid(seq), ids.raw_digest(seq)
    prof = catalog.profile_query(con, "SELECT seq_get(?, ?, ?, 1)", [cu, start, end])
    assert prof["result"][0][0].encode() == seq[start:end]
    scan = _seq_scan(prof, "/seq_chunks/")
    assert scan["files_read"] == 1  # digest_prefix pruning (K-05)
    assert scan["row_groups_scanned"] == _expected_row_groups(root, man, rd, start // CHUNK, (end - 1) // CHUNK)


def test_ly06_protein_lookup_reads_one_protein_file(built, con):
    root, man, _ = built
    f = next(x for x in man["files"] if x["table"] == "protein_seq")
    t = pq.read_table(root / f["path"])
    puid, want = t.column("protein_uid")[0].as_py(), t.column("seq")[0].as_py()
    prof = catalog.profile_query(con, "SELECT protein_of(?)", [puid])
    assert prof["result"][0][0] == want
    scan = _seq_scan(prof, "/protein_seq/")
    n_files = sum(1 for x in man["files"] if x["table"] == "protein_seq")
    assert n_files > 1 and scan["files_read"] == 1 and scan["row_groups_scanned"] == 1


def test_ly06_reverse_lookup_reads_one_row_group_per_file(built, con):
    root, man, _ = built
    occ_files = [x for x in man["files"] if x["table"] == "protein_occurrence"]
    puid = pq.read_table(root / occ_files[0]["path"]).column("protein_uid")[5].as_py()
    prof = catalog.profile_query(con, "SELECT count(*) FROM slots_of(?)", [puid])
    assert prof["result"][0][0] >= 1
    scan = _seq_scan(prof, "/protein_occurrence/")
    assert scan["row_groups_scanned"] <= len(occ_files)


# --- K-05: hive partition pruning ------------------------------------------------------


def test_k05_stage_filter_reads_only_that_stage(built, con):
    _, man, _ = built
    s1 = [f for f in man["files"] if f["table"] == "seq_chunks" and f["partitions"]["stage"] == "S1"]
    allf = [f for f in man["files"] if f["table"] == "seq_chunks"]
    assert 0 < len(s1) < len(allf)
    prof = catalog.profile_query(con, "SELECT count(chunk_ix) FROM seq_chunks WHERE stage = 'S1'")
    assert prof["result"][0][0] == sum(f["rows"] for f in s1)
    assert _seq_scan(prof, "/seq_chunks/")["files_read"] == len(s1)


def test_k05_gene_call_version_and_stage_filter(built, con):
    root, man, _ = built
    want = [f for f in man["files"] if f["table"] == "protein_occurrence" and f["partitions"]["stage"] == "S0"]
    # a sum, not count(*): DuckDB answers count(*) on hive keys without a profiled scan
    prof = catalog.profile_query(con, "SELECT sum(slot_ordinal) FROM protein_occurrence WHERE gene_call_version = ? AND stage = 'S0'",
                                 [GENE_CALL_VERSION])
    assert prof["result"][0][0] == sum(sum(pq.read_table(root / f["path"]).column("slot_ordinal").to_pylist()) for f in want)
    assert _seq_scan(prof, "/protein_occurrence/")["files_read"] == len(want)


# --- K-01: version pin --------------------------------------------------------------------


def test_k01_duckdb_version_matches_the_lock(con):
    lock = tomllib.loads((REPO / "uv.lock").read_text())
    locked = next(p["version"] for p in lock["package"] if p["name"] == "duckdb")
    assert duckdb.__version__ == locked
    assert con.execute("SELECT version()").fetchone()[0].lstrip("v") == locked
    assert dict(con.execute("SELECT * FROM catalog_info").fetchall())["duckdb_version"] == locked


# --- K-02 / K-03: memory limit spills to a local temp directory ----------------------------


def test_k02_k03_oversized_group_by_spills_to_temp_directory(tmp_path):
    spill = tmp_path / "spill"
    seen = {"max": 0}
    stop = threading.Event()

    def watch():
        while not stop.is_set():
            if spill.exists():
                seen["max"] = max(seen["max"], sum(p.stat().st_size for p in spill.rglob("*") if p.is_file()))
            time.sleep(0.01)

    t = threading.Thread(target=watch)
    t.start()
    try:
        c = duckdb.connect(config={"memory_limit": "100MB", "temp_directory": str(spill), "threads": 2})
        n, total = c.execute(
            "SELECT count(*), sum(c) FROM (SELECT i % 5000000 AS k, count(*) AS c, any_value(md5(i::VARCHAR)) "
            "FROM range(15000000) t(i) GROUP BY k)").fetchone()
        c.close()
    finally:
        stop.set()
        t.join()
    assert (n, total) == (5_000_000, 15_000_000)
    assert seen["max"] > 100 * 2**20, "the query fit in memory, so nothing was tested"
    assert not list(tmp_path.glob("*.tmp")) and not (REPO / ".tmp").exists()


def test_k04_thread_setting_is_honored(built):
    c = catalog.open_catalog(built[2], threads=3, memory_limit="1GB")
    assert c.execute("SELECT current_setting('threads')").fetchone()[0] == 3
    c.close()


# --- K-08: read-only concurrency ---------------------------------------------------------


def _client(args):
    cat, cu, n = args
    c = catalog.open_catalog(cat)
    times = []
    for i in range(n):
        t0 = time.perf_counter()
        c.execute("SELECT seq_get(?, ?, ?, 1)", [cu, i * 1000, i * 1000 + 500]).fetchone()
        times.append(time.perf_counter() - t0)
    c.close()
    return sorted(times)[int(0.95 * (n - 1))]


@pytest.mark.parametrize("clients", [1, 4, 16])
def test_k08_concurrent_read_only_clients(built, clients):
    cat = built[2]
    cu = ids.contig_uid(fasta.read(FIX / "NC_000908.2.fna.gz")[0].seq)
    with ProcessPoolExecutor(max_workers=clients, mp_context=get_context("spawn")) as ex:
        p95 = list(ex.map(_client, [(cat, cu, 20)] * clients))
    assert len(p95) == clients  # no lock errors; p95 per client is recorded by R4 on the pod


# --- K-09: writer exclusivity ------------------------------------------------------------


def test_k09_second_writer_fails_fast_and_nothing_is_corrupted(built, tmp_path):
    build = tmp_path / "build.duckdb"
    build.write_bytes(Path(built[2]).read_bytes())
    writer = duckdb.connect(str(build))
    try:
        code = "import duckdb, sys\ntry:\n    duckdb.connect(sys.argv[1])\nexcept duckdb.IOException as e:\n    print('locked', e)\n"
        t0 = time.perf_counter()
        out = subprocess.run([sys.executable, "-c", code, str(build)], capture_output=True, text=True, timeout=60)
        elapsed = time.perf_counter() - t0
    finally:
        writer.close()
    assert out.stdout.startswith("locked") and "lock" in out.stdout.lower(), out.stdout + out.stderr
    assert elapsed < 30
    c = duckdb.connect(str(build), read_only=True)
    assert dict(c.execute("SELECT * FROM catalog_info").fetchall())["release_id"] == "rS1.test"
    c.close()


# --- K-11: sorted COPY is deterministic across thread counts --------------------------------


def test_k11_sorted_copy_is_identical_across_threads(tmp_path):
    hashes, row_hashes = set(), set()
    for threads in (1, 4):
        c = duckdb.connect(config={"threads": threads})
        out = tmp_path / f"t{threads}.parquet"
        c.execute(f"COPY (SELECT (i * 7919) % 100003 AS k, md5(i::VARCHAR) AS v FROM range(300000) t(i) ORDER BY k, v) "
                  f"TO '{out}' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 50000)")
        rows = c.execute(f"SELECT k, v FROM read_parquet('{out}')").fetchall()
        c.close()
        hashes.add(hashlib.sha256(out.read_bytes()).hexdigest())
        row_hashes.add(hashlib.sha256(repr(rows).encode()).hexdigest())
    assert len(row_hashes) == 1  # K-11 / D-02: sorted-row hashes
    assert len(hashes) == 1      # LY-08 holds for the DuckDB writer too at this size


# --- K-12: IN-list and semi-join agree for growing k ---------------------------------------


def test_k12_in_list_and_semi_join_agree(built, con):
    root, man, _ = built
    occ = [r for f in man["files"] if f["table"] == "protein_occurrence" for r in pq.read_table(root / f["path"]).to_pylist()]
    uids = sorted({r["protein_uid"] for r in occ})
    import pyarrow as pa
    for k in (10, 100, len(uids)):
        want = sorted((r["contig_uid"], r["slot_ordinal"]) for r in occ if r["protein_uid"] in set(uids[:k]))
        got_in = sorted(con.execute(
            "SELECT contig_uid, slot_ordinal FROM protein_occurrence WHERE list_contains(?, protein_uid)", [uids[:k]]).fetchall())
        q = pa.table({"protein_uid": pa.array(uids[:k], pa.binary(16))})
        con.register("q", q)
        got_semi = sorted(con.execute(
            "SELECT contig_uid, slot_ordinal FROM protein_occurrence SEMI JOIN q USING (protein_uid)").fetchall())
        con.unregister("q")
        assert got_in == got_semi == want, k
