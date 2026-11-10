"""Ring R2: contracts, layout and determinism over the whole micro-S0 build.

The build is O-INGEST, O-CALL, O-SCAN (a 6-model library made from lambda proteins) and,
where MMseqs2 is on PATH, O-CLUSTER, then O-TAPE and O-PAIRS. Covers LY-02, LY-03, LY-04, C-03 (tables and shards), C-05, R-06, R-09,
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

from biobox import catalog, contracts, fasta, ids, io, ledger, provenance
from biobox.l0seq import SourceObject, ingest
from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, run as call_genes
from biobox.l3 import postings, tapes, tokens
from biobox.registry.schemas import table

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
    tapes.run(root, "S0", {u for g in genomes for u in g.circular})
    tokens.run(root, "rS0.test")
    postings.run(root, "rS0.test")
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


@pytest.mark.covers("LY-04")
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


@pytest.mark.covers("C-03", "LY-02")
def test_c03_shards_are_not_ordered_across_each_other(tmp_path):
    """card B-15: one batch writes one shard, sorted inside; there is no order across shards."""
    from biobox.l0seq.bgzf import read_fai

    root = tmp_path / "s"
    ingest([SourceObject(FIX / o["file"], "s", "1") for o in CONTRACT["objects"]], root, "S0", shard_bases=200_000)
    man = catalog.build_manifest(root, "x")
    shards = sorted(a["path"] for a in man["artifacts"] if a["path"].endswith(".fa.gz"))
    names = [[e.name for e in read_fai(root / (s + ".fai"))] for s in shards]
    assert len(shards) == 3
    assert any(names[i][-1] >= names[i + 1][0] for i in range(len(shards) - 1))  # precondition
    assert contracts.check_shards(man) == []
    # a locator row naming a shard that does not hold the record is caught, on both shards
    loc = root / next(f["path"] for f in man["files"] if f["table"] == "seq_locator")
    rows = pq.read_table(loc).to_pylist()
    rows[0]["shard"], rows[1]["shard"] = rows[1]["shard"], rows[0]["shard"]
    pq.write_table(pa.Table.from_pylist(rows, schema=table("seq_locator").schema), loc)
    probs = contracts.check_shards(man)
    assert len(probs) == 2 and all(".fai vs seq_locator" in p for p in probs), probs


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
                 "l0_manifest", "seq_locator", "cluster_member", "cluster_member_aln", "igr",
                 "contig_tapes", "tape_tokens", "token_dict", "token_postings"):
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
    # a posting whose ordinal lies past the end of its tape (C-05; the pair macros read slots from postings)
    post = con.sql("SELECT * FROM token_postings").to_arrow_table()
    rows = post.to_pylist()
    rows[0]["ordinal"] = rows[0]["n_slots"]
    probs = contracts.check_coordinates(_fault_con(con, token_postings=pa.Table.from_pylist(rows, schema=post.schema)))
    assert probs == ["token_postings: 1 rows with ordinals outside their tape"]
    # a slot whose protein is missing (C-03)
    ps = con.sql("SELECT * FROM protein_seq").to_arrow_table()
    assert contracts.check_foreign_keys(_fault_con(con, protein_seq=ps.slice(1)))
    # an occurrence row lost, so a stored sequence is orphaned (R-06)
    co = con.sql("SELECT * FROM contig_occurrence").to_arrow_table()
    assert contracts.check_lineage(_fault_con(con, contig_occurrence=co.slice(1)))


@pytest.mark.covers("C-08")
def test_c08_provenance_matches_the_container_lock(built):
    _, man, _ = built
    assert provenance.check_provenance(man) == []
    steps = {f["partitions"]["step"] for f in man["files"] if f["table"] == "ledger"}
    assert {"ingest", "genecall", "domains-test-1"} <= steps


@pytest.mark.covers("C-08")
def test_c08_faults_are_caught(built, tmp_path):
    _, man, _ = built
    lock = provenance.locked_versions()
    # a pyrodigal upgrade in the image that the gene-call version string does not name
    assert any("differs from the lock" in p for p in provenance.check_provenance(man, {**lock, "pyrodigal": "0.0.1"}))
    # a step that ran with a tool other than the locked one
    root = tmp_path / "s"
    shutil.copytree(man["root"], root)
    rows = ledger.read(root, "S0", "domains-test-1")
    for r in rows:
        r["tool_versions_json"] = '{"pyarrow": "1.0", "pyhmmer": "0.0.1"}'
    ledger.write(root, "S0", "domains-test-1", rows)
    probs = provenance.check_provenance({**man, "root": str(root)})
    assert any("domains-test-1 tool versions" in p for p in probs)
    # a step that recorded no versions
    rows = ledger.read(root, "S0", "genecall")
    for r in rows:
        r["tool_versions_json"] = None
    ledger.write(root, "S0", "genecall", rows)
    assert any("genecall ledger row without tool versions" in p for p in provenance.check_provenance({**man, "root": str(root)}))
    # a versioned table outside its version partition
    f = next(x for x in man["files"] if x["table"] == "protein_seq")
    bad = {**f, "partitions": {k: v for k, v in f["partitions"].items() if k != "gene_call_version"}}
    assert provenance.check_provenance({**man, "files": [bad]})


@pytest.mark.covers("C-08")
def test_c08_a_tape_version_is_checked_through_its_producers(built, tmp_path):
    """card B-07: a tape_version names no tool versions; C-08 resolves its producers in tape_version_def."""
    import hashlib

    _, man, _ = built
    lock = provenance.locked_versions()
    tv = "t-" + hashlib.sha256(GENE_CALL_VERSION.encode()).hexdigest()[:12]
    tape_files = [f for f in man["files"] if f["table"] in ("contig_tapes", "tape_tokens")]
    defs = [f for f in man["files"] if f["table"] == "tape_version_def"]
    assert tape_files and [f["partitions"] for f in defs] == [{"tape_version": tv}]
    only = {**man, "files": tape_files + defs}
    assert provenance.check_provenance(only) == []
    # the tapes alone hold no gene_call_version: the stale pyrodigal is found through the definition
    probs = provenance.check_provenance(only, {**lock, "pyrodigal": "0.0.1"})
    assert probs and all("differs from the lock" in p and tv in p for p in probs)
    # a tape version without its definition
    probs = provenance.check_provenance({**man, "files": tape_files})
    assert probs and all("no tape_version_def" in p for p in probs)
    # definitions that do not hash to their tape_version, or name a producer that names no tool versions
    root = tmp_path / "s"
    spec = table("tape_version_def")
    cases = {"x": ([GENE_CALL_VERSION, "test-ncrna"], tv, "does not hash"),
             "y": (["test-ncrna"], "t-" + hashlib.sha256(b"test-ncrna").hexdigest()[:12], "does not name its tool versions")}
    for name, (producers, version, msg) in cases.items():
        rel = f"{name}/l3/tape_version_def/tape_version={version}/part-0.parquet"
        io.write_table([pa.table({"producer": sorted(producers)})], spec, root / rel)
        entry = {"path": rel, "table": "tape_version_def", "partitions": {"tape_version": version}}
        probs = provenance.check_provenance({**man, "root": str(root), "files": [entry]})
        assert probs and all(msg in p for p in probs), probs


@pytest.mark.covers("C-08")
def test_c08_an_igr_version_is_a_tape_version_plus_igr1(built):
    """card B-20: `igr_version` = `<tape_version>-igr1`; its tape version must be defined in the manifest."""
    import hashlib

    _, man, _ = built
    tv = "t-" + hashlib.sha256(GENE_CALL_VERSION.encode()).hexdigest()[:12]
    defs = [f for f in man["files"] if f["table"] == "tape_version_def"]
    assert [f["partitions"] for f in defs] == [{"tape_version": tv}]

    def igr_file(version):
        return {"path": f"l0seq/igr/igr_version={version}/stage=S0/part-0.parquet", "table": "igr",
                "partitions": {"igr_version": version, "stage": "S0"}}

    ok = {**man, "files": defs + [igr_file(f"{tv}-igr1")]}
    assert provenance.check_provenance(ok) == []
    # a stale producer is found through the IGR's tape version
    lock = provenance.locked_versions()
    probs = provenance.check_provenance(ok, {**lock, "pyrodigal": "0.0.1"})
    assert probs and all("differs from the lock" in p for p in probs)
    # a tape version without a definition
    probs = provenance.check_provenance({**man, "files": [igr_file(f"{tv}-igr1")]})
    assert probs and all("no tape_version_def" in p and "igr_version" in p for p in probs), probs
    # a malformed suffix, and the old form (a bare gene-call version)
    for bad in (f"{tv}-igr2", tv, GENE_CALL_VERSION + "-igr1", "-igr1"):
        probs = provenance.check_provenance({**man, "files": defs + [igr_file(bad)]})
        assert any("igr_version" in p and "<tape_version>-igr1" in p for p in probs), (bad, probs)


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
