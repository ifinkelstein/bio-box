"""Ring R2: contracts, layout and determinism over the whole micro-S0 build.

The build is O-INGEST, O-CALL, O-SCAN (a 6-model library made from lambda proteins) and,
where MMseqs2 is on PATH, O-CLUSTER. Covers LY-02, LY-03, LY-04, C-03 (tables and shards), C-05, R-06, R-09,
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

from biobox import catalog, contracts, fasta, ids, ledger
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
            for o in objs], root, "S0", shard_bases=10**6)
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
        tmp.mkdir(parents=True, exist_ok=True)
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
    f = next(x for x in man["files"] if x["table"] == "contig_occurrence")
    pq.write_table(pq.read_table(Path(man["root"]) / f["path"]), tmp_path / "x.parquet", write_statistics=False)
    assert contracts.check_statistics({"root": str(tmp_path), "files": [{**f, "path": "x.parquet"}]})


@pytest.mark.xfail(strict=True, reason="LY-04 finding (2026-10-30): 256 fixed protein prefix buckets split "
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


@pytest.mark.covers("C-03", "LY-02")
def test_c03_shards_agree_with_the_locator(built):
    _, man, _ = built
    assert contracts.check_shards(man) == []


@pytest.mark.covers("C-03", "LY-02")
@pytest.mark.parametrize("fault", ["dropped fai line", "wrong length", "unsorted", "bad gzi", "locator row lost"])
def test_c03_shard_faults_are_caught(built, tmp_path, fault):
    import struct

    from biobox.l0seq.bgzf import read_fai, read_gzi

    _, man, _ = built
    root = tmp_path / "s"
    shutil.copytree(man["root"], root)
    shard = next(a["path"] for a in man["artifacts"] if a["path"].endswith(".fa.gz") and
                 len(read_fai(root / (a["path"] + ".fai"))) > 1)
    fai = root / (shard + ".fai")
    lines = fai.read_text().splitlines(keepends=True)
    if fault == "dropped fai line":
        fai.write_text("".join(lines[1:]))
    elif fault == "wrong length":
        f = lines[0].split("\t")
        fai.write_text("".join(["\t".join([f[0], str(int(f[1]) + 1)] + f[2:])] + lines[1:]))
    elif fault == "unsorted":
        fai.write_text("".join(lines[::-1]))
    elif fault == "bad gzi":
        gzi = root / (shard + ".gzi")
        pairs = read_gzi(gzi)[1:]  # a block start missing
        gzi.write_bytes(struct.pack("<Q", len(pairs)) + b"".join(struct.pack("<QQ", c, u) for c, u in pairs))
    else:
        loc = next(f for f in man["files"] if f["table"] == "seq_locator")
        pq.write_table(pq.read_table(root / loc["path"]).slice(1), root / loc["path"])
    assert contracts.check_shards({**man, "root": str(root)})


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
                 "l0_manifest", "seq_locator", "cluster_member", "cluster_member_aln", "igr"):
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
    # an occurrence row lost, so a stored sequence is orphaned (R-06)
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
           shard_bases=10**6)
    after = catalog.build_manifest(root, "rS1.test")
    assert contracts.check_append_only(before, after) == []
    assert len(after["files"]) > len(before["files"])
    # and the check notices a rewrite
    assert len(after["artifacts"]) > len(before["artifacts"])  # the lambda variant is a new S1 shard
    f = next(x for x in after["files"] if x["table"] == "contig_occurrence")
    assert contracts.check_append_only(before, {**after, "files": [g for g in after["files"] if g["path"] != f["path"]]})
    shard = next(a for a in after["artifacts"] if "stage=S0" in a["path"])
    rewritten = [{**a, "sha256": "0" * 64} if a is shard else a for a in after["artifacts"]]
    assert contracts.check_append_only(before, {**after, "artifacts": rewritten}) == [f"{shard['path']}: rewritten"]


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
    assert "seq_shards" in changed


# --- ledger (L-01..L-03) --------------------------------------------------------------------


def _steps(root: Path) -> list[str]:
    return sorted(p.name.removeprefix("step=") for p in (root / "ledger/stage=S0").iterdir())


def test_l01_every_step_balances(built):
    _, man, _ = built
    root = Path(man["root"])
    steps = _steps(root)
    assert {"ingest", "genecall", "domains-test-1"} <= set(steps)
    assert ledger.check_balance(root, "S0", steps) == []


def test_l02_ingest_ledger_lists_the_input_objects(built):
    _, man, _ = built
    expected = {ids.sha256_file(FIX / o["file"]) for o in CONTRACT["objects"]}
    assert ledger.check_inputs(Path(man["root"]), "S0", expected) == []
    assert ledger.check_inputs(Path(man["root"]), "S0", expected | {b"x" * 32})


def test_l03_counts_carry_along_the_dag(built):
    _, man, _ = built
    root = Path(man["root"])
    assert ledger.check_chain(root, "S0") == []
    gc = ledger.read(root, "S0", "genecall")
    distinct = len(domains.load_proteins(root, GENE_CALL_VERSION, "S0"))
    assert sum(r["n_kept"] for r in gc) == distinct
    if MMSEQS:
        assert {r["object_key"] for r in ledger.read(root, "S0", "cluster")} == {b"nr90", b"nr50"}


def test_ledger_faults_are_caught(built, tmp_path):
    _, man, _ = built
    root = tmp_path / "s"
    shutil.copytree(man["root"], root)
    gc = ledger.read(root, "S0", "genecall")
    gc[0]["n_kept"] += 1  # unbalanced (L-01) and breaks the chain (L-03)
    ledger.write(root, "S0", "genecall", gc)
    assert ledger.check_balance(root, "S0", ["genecall"])
    assert ledger.check_chain(root, "S0")
