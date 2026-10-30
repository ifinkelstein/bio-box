"""Ring R2: contracts, layout and determinism over the whole micro-S0 build.

The build is O-INGEST, O-CALL, O-SCAN (a 6-model library made from lambda proteins) and,
where MMseqs2 is on PATH, O-CLUSTER. Covers LY-02, LY-03, LY-04, C-03, C-05, R-06, R-09,
D-01 (same inputs twice), D-02 (thread counts) and D-04 (input order), each with an
injected fault that the check must catch.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pyhmmer
import pytest
import yaml

from biobox import catalog, contracts, fasta, ids
from biobox.l0seq import SourceObject, ingest
from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, run as call_genes

pytestmark = pytest.mark.r2

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
MMSEQS = shutil.which("mmseqs")


def _library(root: Path, path: Path) -> None:
    occ = pq.read_table(root / f"l1/protein_occurrence/gene_call_version={GENE_CALL_VERSION}/stage=S0/part-0.parquet")
    lam = ids.contig_uid(fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq)
    prots = domains.load_proteins(root, GENE_CALL_VERSION, "S0")
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


def build(root: Path, order: int = 1, threads: int = 1) -> dict:
    objs = CONTRACT["objects"][::order]
    ingest([SourceObject(FIX / o["file"], CONTRACT["source_name"], CONTRACT["source_release"], CONTRACT["license"], True)
            for o in objs], root, "S0", chunk_size=4096, row_group_chunks=16)
    genomes = []
    for o in objs:
        cu = ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)
        genomes.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"],
                                   circular=(cu,) if o["topology"] == "circular" else ()))
    call_genes(root, "S0", genomes)
    lib = root.parent / f"{root.name}.hmm"
    _library(root, lib)
    domains.run(root, "S0", GENE_CALL_VERSION, "test", "1", lib, "gathering", cpus=threads)
    if MMSEQS:
        from biobox.l2 import cluster
        tmp = root.parent / f"{root.name}.mmseqs"
        cluster.run(root, "rS0.test", GENE_CALL_VERSION, threads=threads, tmp=tmp)
        shutil.rmtree(tmp, ignore_errors=True)
    return catalog.build_manifest(root, "rS0.test")


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    d = tmp_path_factory.mktemp("contracts")
    man = build(d / "a")
    cat = catalog.generate_catalog(man, d / "c.duckdb")
    con = catalog.open_catalog(cat)
    yield d, man, con
    con.close()


# --- layout ----------------------------------------------------------------------------


def test_ly02_every_file_is_in_declared_sort_order(built):
    _, man, _ = built
    assert contracts.check_sort_order(man) == []


def test_ly02_catches_an_unsorted_file(built, tmp_path):
    _, man, _ = built
    f = next(x for x in man["files"] if x["table"] == "protein_occurrence")
    t = pq.read_table(Path(man["root"]) / f["path"])
    bad = tmp_path / "x.parquet"
    pq.write_table(t.take(list(range(t.num_rows - 1, -1, -1))), bad)
    assert contracts.check_sort_order({"root": str(tmp_path), "files": [{**f, "path": "x.parquet"}]})


def test_ly03_sort_keys_have_statistics(built):
    _, man, _ = built
    assert contracts.check_statistics(man) == []


def test_ly03_catches_missing_statistics(built, tmp_path):
    _, man, _ = built
    f = next(x for x in man["files"] if x["table"] == "seq_chunks")
    pq.write_table(pq.read_table(Path(man["root"]) / f["path"]), tmp_path / "x.parquet", write_statistics=False)
    assert contracts.check_statistics({"root": str(tmp_path), "files": [{**f, "path": "x.parquet"}]})


@pytest.mark.xfail(strict=True, reason="LY-04 finding (2026-10-30): 256 fixed digest/protein prefix buckets split "
                   "a small stage into many tiny files; the bucket count has to scale with the stage")
def test_ly04_no_small_files_outside_waivers(built):
    _, man, _ = built
    assert contracts.check_file_sizes(man) == []


def test_size_check_allows_one_small_file_per_stage():
    files = [{"table": "slot_protein", "partitions": {"stage": "S0"}, "bytes": 10}]
    assert contracts.check_file_sizes({"files": files}) == []
    assert contracts.check_file_sizes({"files": files * 2})


# --- contracts ---------------------------------------------------------------------------


def test_c03_foreign_keys_resolve(built):
    _, _, con = built
    assert contracts.check_foreign_keys(con) == []


def test_c05_coordinates_in_bounds(built):
    _, _, con = built
    assert contracts.check_coordinates(con) == []  # the wrap rule is exercised on the golden set (R3)


def test_r06_lineage_is_complete(built):
    _, _, con = built
    assert contracts.check_lineage(con) == []


def _fault_con(con, **replace) -> duckdb.DuckDBPyConnection:
    """An in-memory copy of the catalog's tables with some replaced by Arrow tables (fault injection)."""
    c = duckdb.connect()
    for name in ("protein_occurrence", "contig_index", "slot_protein", "protein_seq", "l1_domains", "contig_occurrence",
                 "seq_chunks", "l0_manifest", "seq_locator", "cluster_member", "cluster_member_aln"):
        tbl = replace.get(name)
        if tbl is None:
            tbl = con.sql(f"SELECT * FROM {name}").to_arrow_table()
        c.register(name, tbl)
    return c


def test_injected_faults_are_caught(built):
    _, _, con = built
    occ = con.sql("SELECT * FROM protein_occurrence").to_arrow_table()
    # a slot ending past the end of a linear contig (C-05)
    lin = con.execute("SELECT contig_uid FROM contig_index WHERE contig_len < 100000").fetchall()[0][0]
    rows = occ.to_pylist()
    k = next(n for n, r in enumerate(rows) if r["contig_uid"] == lin)
    rows[k]["end"] = rows[k]["start"] + 10**6
    bad = pa.Table.from_pylist(rows, schema=occ.schema)
    assert contracts.check_coordinates(_fault_con(con, protein_occurrence=bad))
    # a slot whose protein is missing (C-03)
    ps = con.sql("SELECT * FROM protein_seq").to_arrow_table()
    assert contracts.check_foreign_keys(_fault_con(con, protein_seq=ps.slice(1)))
    # an occurrence row lost, so a chunk is orphaned (R-06)
    co = con.sql("SELECT * FROM contig_occurrence").to_arrow_table()
    assert contracts.check_lineage(_fault_con(con, contig_occurrence=co.slice(1)))


# --- append-only and determinism -----------------------------------------------------------


def test_r09_a_later_stage_keeps_earlier_files(built, tmp_path):
    _, man, _ = built
    root = tmp_path / "s"
    shutil.copytree(man["root"], root)
    before = catalog.build_manifest(root, "rS0.test")
    t4 = next(o for o in CONTRACT["objects"] if o["accession"] == "NC_000866.4")
    dup = tmp_path / "t4_again.fna.gz"
    shutil.copy(FIX / t4["file"], dup)
    lam_masked = tmp_path / "new.fna"
    lam = fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq
    lam_masked.write_bytes(b"\n".join(fasta.format_record("lam_variant", lam[:2000].lower() + lam[2000:], 70)) + b"\n")
    ingest([SourceObject(dup, "s1_source", "1"), SourceObject(lam_masked, "s1_source", "1")], root, "S1",
           chunk_size=4096, row_group_chunks=16)
    after = catalog.build_manifest(root, "rS1.test")
    assert contracts.check_append_only(before, after) == []
    assert len(after["files"]) > len(before["files"])
    # and the check notices a rewrite
    f = next(x for x in after["files"] if x["table"] == "seq_chunks")
    assert contracts.check_append_only(before, {"files": [g for g in after["files"] if g["path"] != f["path"]]})


def test_d01_same_inputs_same_content(built, tmp_path):
    _, man, _ = built
    assert contracts.content_hashes(build(tmp_path / "b")) == contracts.content_hashes(man)


def test_d02_thread_count_does_not_change_content(built, tmp_path):
    _, man, _ = built
    assert contracts.content_hashes(build(tmp_path / "t4", threads=4)) == contracts.content_hashes(man)


def test_d04_input_order_does_not_change_content(built, tmp_path):
    _, man, _ = built
    assert contracts.content_hashes(build(tmp_path / "rev", order=-1)) == contracts.content_hashes(man)


def test_content_hash_sees_a_changed_row(built, tmp_path):
    _, man, _ = built
    root = tmp_path / "s"
    shutil.copytree(man["root"], root)
    f = next(x for x in man["files"] if x["table"] == "l1_domains")
    t = pq.read_table(root / f["path"]).to_pylist()
    t[0]["dom_bitscore"] += 1.0
    pq.write_table(pa.Table.from_pylist(t, schema=pq.read_schema(root / f["path"])), root / f["path"])
    changed = contracts.content_hashes({**man, "root": str(root)})
    assert changed["l1_domains"] != contracts.content_hashes(man)["l1_domains"]
    assert {k: v for k, v in changed.items() if k != "l1_domains"} == \
        {k: v for k, v in contracts.content_hashes(man).items() if k != "l1_domains"}
