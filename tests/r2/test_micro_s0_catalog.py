"""Ring R2: O-CATALOG (manifest, generated read-only catalog, sequence macros) on micro-S0.

Covers LY-05 (partition paths and derived partition values), LY-07 (DuckDB, pyarrow and
polars read identical rows from every table), LY-10 (explicit file lists, sha256 per file,
strays refused, shards with their indexes), C-01 (schema drift refused at manifest time),
R-03 (the DuckDB macros return the same bytes as the Python path, both strands, BGZF block
and shard edges, origin wrap), MR-01
(manifest hash independent of input order and of the data root) and the atomic swap of
`current.json` (data_layer_and_tests.md §3.3).
"""

from __future__ import annotations

import hashlib
import random
import shutil
from datetime import datetime
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

from biobox import catalog, fasta, ids
from biobox.l0seq import SeqStore, SourceObject, ingest, revcomp
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, run as call_genes
from biobox.l3 import tapes
from biobox.registry.schemas import TABLES

pytestmark = pytest.mark.r2

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
SHARD_BASES = 200_000  # the micro set spans several shards


def _build(root: Path, order=1) -> Path:
    objs = CONTRACT["objects"][::order]
    ingest([SourceObject(FIX / o["file"], CONTRACT["source_name"], CONTRACT["source_release"], CONTRACT["license"], True)
            for o in objs], root, "S0", shard_bases=SHARD_BASES)
    genomes = []
    for o in objs:
        cu = ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)
        genomes.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"],
                                   circular=(cu,) if o["topology"] == "circular" else ()))
    call_genes(root, "S0", genomes)
    tapes.run(root, "S0", {u for g in genomes for u in g.circular})  # LIST columns for LY-07
    return root


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = _build(tmp_path_factory.mktemp("catalog"))
    man = catalog.build_manifest(root, "rS0.test")
    catalog.write_manifest(man, root / "manifests/release=rS0.test.json")
    cat = catalog.generate_catalog(man, root / "catalog/release=rS0.test.duckdb")
    return root, man, cat


@pytest.fixture()
def con(built, tmp_path):
    _, _, cat = built
    c = catalog.open_catalog(cat, local_dir=tmp_path / "local")
    yield c
    c.close()


def _contigs():
    out = []
    for o in CONTRACT["objects"]:
        seq = fasta.read(FIX / o["file"])[0].seq
        out.append((ids.contig_uid(seq), seq, o["topology"] == "circular"))
    return out


# --- manifest (LY-10, LY-05, C-01) ------------------------------------------------


@pytest.mark.covers("LY-10")
def test_manifest_lists_every_parquet_file_and_verifies(built):
    root, man, _ = built
    on_disk = sorted(p.relative_to(root).as_posix() for p in root.rglob("*.parquet"))
    assert [f["path"] for f in man["files"]] == on_disk
    assert catalog.verify_manifest(man) == []
    for f in man["files"]:
        assert f["rows"] == pq.ParquetFile(root / f["path"]).metadata.num_rows
    assert {f["table"] for f in man["files"]} >= {"seq_locator", "contig_index", "protein_seq", "feature", "feature_protein",
                                                  "protein_occurrence"}
    shard_files = sorted(p.relative_to(root).as_posix() for p in (root / "l0seq/seq_shards").rglob("*") if p.is_file())
    assert [a["path"] for a in man["artifacts"]] == shard_files
    assert len(shard_files) == 3 * len([a for a in shard_files if a.endswith(".fa.gz")]) > 3
    from biobox import contracts

    assert contracts.check_shards(man) == []  # order across several shards (LY-02)


@pytest.mark.covers("LY-10")
def test_verify_detects_tampering_and_strays(built, tmp_path):
    root, man, _ = built
    copy = tmp_path / "copy"
    shutil.copytree(root, copy)
    f = next(x for x in man["files"] if x["table"] == "seq_locator")
    data = bytearray((copy / f["path"]).read_bytes())
    data[-20] ^= 0xFF
    (copy / f["path"]).write_bytes(bytes(data))
    shutil.copy(copy / f["path"], copy / "l0seq/stray.parquet")
    shard = next(a for a in man["artifacts"] if a["path"].endswith(".fa.gz"))
    data = bytearray((copy / shard["path"]).read_bytes())
    data[100] ^= 0xFF
    (copy / shard["path"]).write_bytes(bytes(data))
    extra = shard["path"].replace("shard-0000", "shard-0099")
    shutil.copy(copy / shard["path"], copy / extra)
    problems = catalog.verify_manifest(man, copy)
    assert f"sha256 differs: {f['path']}" in problems
    assert "not in manifest: l0seq/stray.parquet" in problems
    assert f"sha256 differs: {shard['path']}" in problems
    assert f"not in manifest: {extra}" in problems


@pytest.mark.covers("C-01", "LY-05", "LY-10")
def test_stray_wrong_partition_and_schema_drift_are_refused(built, tmp_path):
    root, man, _ = built
    occ = next(x for x in man["files"] if x["table"] == "contig_occurrence")
    cases = {
        "stray": "l9/unknown/part-0.parquet",
        "missing partition": "l0seq/contig_occurrence/stage=S0/part-0.parquet",
        "partition order": f"l0seq/contig_occurrence/source={occ['partitions']['source']}/stage=S0/part-0.parquet",
    }
    for why, rel in cases.items():
        d = tmp_path / why.replace(" ", "_")
        (d / rel).parent.mkdir(parents=True)
        shutil.copy(root / occ["path"], d / rel)
        with pytest.raises(catalog.LayoutError):
            catalog.build_manifest(d, "x")
    # a protein_seq file under the wrong protein_prefix (LY-05)
    prot = next(x for x in man["files"] if x["table"] == "protein_seq")
    d = tmp_path / "prefix"
    have = prot["partitions"]["protein_prefix"]
    rel = prot["path"].replace(f"protein_prefix={have}", f"protein_prefix={'ff' if have != 'ff' else '00'}")
    (d / rel).parent.mkdir(parents=True)
    shutil.copy(root / prot["path"], d / rel)
    with pytest.raises(catalog.LayoutError, match="LY-05"):
        catalog.build_manifest(d, "x")
    # a file whose schema is not the registered one (C-01)
    d = tmp_path / "schema"
    occ = next(x for x in man["files"] if x["table"] == "contig_index")
    (d / occ["path"]).parent.mkdir(parents=True)
    t = pq.read_table(root / occ["path"])
    pq.write_table(t.drop_columns(["is_default_variant"]), d / occ["path"])
    with pytest.raises(catalog.LayoutError, match="C-01"):
        catalog.build_manifest(d, "x")


def test_mr01_manifest_hash_ignores_input_order_and_root(built, tmp_path):
    _, man, _ = built
    other = _build(tmp_path / "reversed", order=-1)
    m2 = catalog.build_manifest(other, "rS0.test")
    assert m2["root"] != man["root"]
    assert catalog.manifest_hash(m2) == catalog.manifest_hash(man)


@pytest.mark.covers("LY-05", "LY-10")
@pytest.mark.parametrize("case", ["no gzi", "orphan fai", "misnamed", "no stage", "parquet in shards"])
def test_bad_shard_layouts_are_refused(built, tmp_path, case):
    root, man, _ = built
    shard = next(a["path"] for a in man["artifacts"] if a["path"].endswith(".fa.gz"))
    d = tmp_path / "s"
    shutil.copytree(root / "l0seq/seq_shards", d / "l0seq/seq_shards")
    if case == "no gzi":
        (d / (shard + ".gzi")).unlink()
    elif case == "orphan fai":
        (d / shard).unlink()
    elif case == "misnamed":
        shutil.copy(d / shard, d / shard.replace("shard-0000", "contigs"))
    elif case == "no stage":
        shutil.copy(d / shard, d / "l0seq/seq_shards/shard-0000.fa.gz")
    else:
        shutil.copy(root / man["files"][0]["path"], d / "l0seq/seq_shards/stage=S0/part-0.parquet")
    with pytest.raises(catalog.LayoutError):
        catalog.build_manifest(d, "x")


# --- catalog views (LY-10, LY-07) ---------------------------------------------------


@pytest.mark.covers("LY-10")
def test_views_use_explicit_file_lists(built, con):
    _, man, _ = built
    views = dict(con.execute("SELECT view_name, sql FROM duckdb_views() WHERE NOT internal").fetchall())
    assert set(TABLES) <= set(views)
    for name, sql in views.items():
        assert "*" not in sql.split("FROM", 1)[-1].split("read_parquet", 1)[-1].split(",")[0], f"{name}: glob in view"
    listed = {f["path"] for f in man["files"]}
    for f in man["files"]:
        assert f"/{f['path']}'" in views[f["table"]]
    assert len(listed) == con.execute("SELECT count(*) FROM catalog_files").fetchone()[0]
    assert len(man["artifacts"]) == con.execute("SELECT count(*) FROM catalog_artifacts").fetchall()[0][0]
    info = dict(con.execute("SELECT * FROM catalog_info").fetchall())
    assert info["manifest_hash"] == catalog.manifest_hash(man)
    assert info["duckdb_version"] == duckdb.__version__


def test_view_row_counts_equal_manifest(built, con):
    _, man, _ = built
    for name in TABLES:
        want = sum(f["rows"] for f in man["files"] if f["table"] == name)
        assert con.execute(f"SELECT count(*) FROM {name}").fetchone()[0] == want, name


def test_empty_tables_are_typed_views(con):
    cols = con.execute("DESCRIBE cluster_member").fetchall()
    assert [c[0] for c in cols] == [f.name for f in TABLES["cluster_member"].schema] + ["release"]
    assert con.execute("SELECT count(*) FROM search_hit").fetchone()[0] == 0


def _canon(v):
    # timestamps compare as epoch ms: readers differ only in the time zone they attach
    return round(v.timestamp() * 1000) if isinstance(v, datetime) else v


def _digest(rows) -> str:
    return hashlib.sha256(repr(sorted(tuple(_canon(v) for v in r) for r in rows)).encode()).hexdigest()


def test_ly07_duckdb_pyarrow_polars_read_identical_rows(built, con):
    root, man, _ = built
    polars = pytest.importorskip("polars")
    for name, spec in TABLES.items():
        files = [root / f["path"] for f in man["files"] if f["table"] == name]
        if not files:
            continue
        cols = ", ".join(f'"{c}"' for c in spec.schema.names)
        # through Arrow: DuckDB's Python row fetch of TIMESTAMPTZ needs pytz
        duck = [tuple(r.values()) for r in con.sql(f"SELECT {cols} FROM {name}").to_arrow_table().to_pylist()]
        arrow = [tuple(r.values()) for f in files for r in pq.read_table(f).to_pylist()]
        pl = [r for f in files for r in polars.read_parquet(f).rows()]
        assert len(duck) == len(arrow) == len(pl), name
        assert _digest(duck) == _digest(arrow) == _digest(pl), name


def test_catalog_is_read_only(con):
    with pytest.raises(duckdb.Error):
        con.execute("CREATE TABLE x (a INTEGER)")


@pytest.mark.covers("LY-05")
def test_partition_columns_come_from_paths(con):
    stages = {r[0] for r in con.execute("SELECT DISTINCT stage FROM contig_index").fetchall()}
    assert stages == {"S0"}
    bad = con.execute("SELECT count(*) FROM protein_seq WHERE protein_prefix <> 'all' "
                      "AND NOT starts_with(lower(hex(protein_uid[1:1])), protein_prefix)").fetchall()[0][0]
    assert bad == 0


# --- macros vs the Python path (R-03) ------------------------------------------------


def _block_edges(store: SeqStore, rd: bytes) -> list[int]:
    """Sequence positions where one BGZF block of the record's shard ends and the next begins."""
    bg = store._shard(store.shard_of(rd))
    e = bg.entry(rd.hex())
    out = []
    for u in bg._ustart[1:]:
        rel = u - e.offset
        if 0 < rel < e.byte_of(e.length - 1) - e.offset:
            out.append(rel - rel // e.line_width)  # newlines before the edge are not bases
    return out


def test_r03_macros_equal_python_path(built, con):
    root, _, _ = built
    store = SeqStore(root)
    rng = random.Random(20261030)
    n_edges = 0
    for cu, seq, _ in _contigs():
        n = len(seq)
        edges = _block_edges(store, ids.raw_digest(seq))
        n_edges += len(edges)
        cases = [(0, n), (0, 0), (n, n), (0, 1), (n - 1, n)]
        cases += [(max(0, p - k), min(n, p + k)) for p in edges for k in (1, 81)]
        cases += [(p, p) for p in edges]
        cases += [tuple(sorted(rng.randrange(0, n + 1) for _ in range(2))) for _ in range(60)]
        for a, b in cases:
            for strand in (1, -1):
                got = con.execute("SELECT seq_get(?, ?, ?, ?)", [cu, a, b, strand]).fetchone()[0]
                assert got.encode() == store.get(cu, a, b, strand), (a, b, strand)
                want = seq[a:b] if strand == 1 else revcomp(seq[a:b])
                assert got.encode() == want
    assert n_edges > 5


def test_r03_origin_wrap_on_circular_contig(built, con):
    root, _, _ = built
    store = SeqStore(root)
    cu, seq, circ = next(c for c in _contigs() if c[2])
    n = len(seq)
    for a, b in [(n - 100, n + 50), (n - 1, n + 1), (5, n + 4)]:
        for strand in (1, -1):
            got = con.execute("SELECT seq_get_circular(?, ?, ?, ?)", [cu, a, b, strand]).fetchone()[0]
            assert got.encode() == store.get(cu, a, b, strand, circular=True)
    lin = con.execute("SELECT seq_get_circular(?, 10, 20, 1)", [cu]).fetchone()[0]
    assert lin.encode() == seq[10:20]


def test_mr14_macro_reads_only_the_default_variant(tmp_path):
    lam_path = FIX / CONTRACT["objects"][0]["file"]
    lam = fasta.read(lam_path)[0]
    masked = lam.seq[:5000].lower() + lam.seq[5000:]
    src = tmp_path / "src"
    src.mkdir()
    mpath = src / "lambda_masked.fna"
    mpath.write_bytes(b"\n".join(fasta.format_record("lambda_masked", masked, 70)) + b"\n")
    root = tmp_path / "s"
    ingest([SourceObject(lam_path, CONTRACT["source_name"], CONTRACT["source_release"], CONTRACT["license"], True),
            SourceObject(mpath, "other_source", "1", "INSDC (public)", True, precedence=200)],
           root, "S0", shard_bases=SHARD_BASES)
    cat = catalog.generate_catalog(catalog.build_manifest(root, "x"), tmp_path / "c.duckdb")
    c = catalog.open_catalog(cat)
    cu = ids.contig_uid(lam.seq)
    assert c.execute("SELECT seq_get(?, 0, 6000, 1)", [cu]).fetchone()[0].encode() == lam.seq[:6000]
    assert c.execute("SELECT seq_raw(?, 0, 6000)", [ids.raw_digest(masked)]).fetchone()[0].encode() == masked[:6000]
    c.close()


def test_w08_contig_label_has_one_row_per_contig(tmp_path):
    lam_path = FIX / CONTRACT["objects"][0]["file"]
    lam = fasta.read(lam_path)[0]
    src = tmp_path / "src"
    src.mkdir()
    dup = src / "lambda_again.fna"
    dup.write_bytes(b"\n".join(fasta.format_record("lambda_again", lam.seq, 70)) + b"\n")
    root = tmp_path / "s"
    ingest([SourceObject(lam_path, CONTRACT["source_name"], CONTRACT["source_release"], CONTRACT["license"], True),
            SourceObject(dup, "other_source", "1", "INSDC (public)", True, precedence=200)],
           root, "S0", shard_bases=SHARD_BASES)
    man = catalog.build_manifest(root, "x")
    c = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / "c.duckdb"))
    occ = c.execute("SELECT l0_object_sha256, record_ordinal, accession, original_header, contig_uid "
                    "FROM contig_occurrence ORDER BY l0_object_sha256, record_ordinal").fetchall()
    assert len(occ) == 2 and len({r[4] for r in occ}) == 1
    rows = c.execute("SELECT contig_uid, accession, original_header, n_occurrences FROM contig_label").fetchall()
    assert rows == [(occ[0][4], occ[0][2], occ[0][3], 2)]
    assert c.execute("SELECT count(*) FROM contig_label").fetchone()[0] == c.execute(
        "SELECT count(DISTINCT contig_uid) FROM contig_occurrence").fetchone()[0]
    c.close()


@pytest.mark.parametrize("args", [(-1, 5, 1), (5, 4, 1), (0, 10**9, 1), (0, 5, 0)])
def test_macro_rejects_bad_slices(con, args):
    cu = _contigs()[0][0]
    with pytest.raises(duckdb.Error):
        con.execute("SELECT seq_get(?, ?, ?, ?)", [cu, *args]).fetchone()


def test_sequence_macros_need_open_catalog(built):
    # the message the module docstring documents: the Python function is not registered
    _, _, cat = built
    plain = duckdb.connect(str(cat), read_only=True)
    cu, seq, _ = _contigs()[0]
    for sql, args, fn in [("SELECT seq_get(?, 0, 10, 1)", [cu], "bb_seq"),
                          ("SELECT seq_slice(?, 0, 10)", [cu], "bb_seq"),
                          ("SELECT seq_get_circular(?, 0, 10, 1)", [cu], "bb_seq_circular"),
                          ("SELECT seq_raw(?, 0, 10)", [ids.raw_digest(seq)], "bb_seq_raw")]:
        with pytest.raises(duckdb.Error, match=f"Scalar Function with name {fn} does not exist"):
            plain.execute(sql, args).fetchall()
    assert plain.execute("SELECT contig_len_of(?)", [cu]).fetchall()[0][0] == len(_contigs()[0][1])
    plain.close()


def test_open_catalog_on_another_root(built, tmp_path):
    # the views name the build root's Parquet files; `root` moves only the shard reader
    root, man, cat = built
    moved = tmp_path / "moved"
    shutil.copytree(root / "l0seq", moved / "l0seq")
    cu, seq, _ = _contigs()[0]
    c = catalog.open_catalog(cat, root=moved)
    assert catalog.seq_reader(c).root == moved
    assert c.execute("SELECT seq_get(?, 5, 50, 1)", [cu]).fetchall()[0][0].encode() == seq[5:50]
    c.close()
    shutil.rmtree(moved / "l0seq/seq_shards")
    c = catalog.open_catalog(cat, root=moved)
    with pytest.raises(duckdb.Error, match="No such file"):
        c.execute("SELECT seq_get(?, 0, 10, 1)", [cu]).fetchall()
    c.close()


def test_unknown_contig_raises(con):
    with pytest.raises(duckdb.Error):
        con.execute("SELECT seq_get(?, 0, 1, 1)", [b"\x00" * 16]).fetchone()


# --- B-14: sequence macros over bb_seq, no index scans (AGENTS §4.4) -------------------------


def _oracle(seq: bytes, a: int, b: int, strand: int) -> bytes:
    """The FASTA bytes themselves, wrapping the origin when b > len(seq)."""
    x = seq[a:b] if b <= len(seq) else seq[a:] + seq[: b - len(seq)]
    return x if strand == 1 else revcomp(x)


def _genes(root: Path) -> list[dict]:
    f = next(root.glob("l1/protein_occurrence/gene_call_version=*/stage=S0/part-0.parquet"))
    return pq.read_table(f, columns=["contig_uid", "start", "end", "strand"]).to_pylist()


def _slices_table(rows) -> pa.Table:
    """rows: (contig_uid, start, end, strand); `i` keeps the input order."""
    return pa.table({"i": pa.array(range(len(rows)), pa.int64()),
                     "contig_uid": pa.array([r[0] for r in rows], pa.binary(16)),
                     "s": pa.array([r[1] for r in rows], pa.int64()),
                     "e": pa.array([r[2] for r in rows], pa.int64()),
                     "strand": pa.array([r[3] for r in rows], pa.int8())})


def _random_slices(n: int, seed: int, only=None, max_len: int = 3000) -> list[tuple]:
    rng = random.Random(seed)
    contigs = [c for c in _contigs() if only is None or c[0] == only]
    out = []
    for _ in range(n):
        cu, seq, _ = rng.choice(contigs)
        a = rng.randrange(0, len(seq))
        out.append((cu, a, min(len(seq), a + rng.randrange(0, max_len)), rng.choice((1, -1))))
    return out


def _index_scans(prof: dict) -> list[dict]:
    return [s for s in prof["scans"] if any("/contig_index/" in f or "/seq_locator/" in f for f in s["files"])]


def _check_genes(con, store, seqs, macro, group):
    """Every gene of `group` through `macro` in one table-driven query, against FASTA and SeqStore."""
    con.register("g", _slices_table([(g["contig_uid"], g["start"], g["end"], g["strand"]) for g in group]))
    got = con.execute(f"SELECT i, {macro}(contig_uid, s, e, strand) FROM g ORDER BY i").fetchall()
    con.unregister("g")
    assert [r[0] for r in got] == list(range(len(group)))
    for g, (_, x) in zip(group, got):
        seq, is_circ = seqs[g["contig_uid"]]
        assert x.encode() == _oracle(seq, g["start"], g["end"], g["strand"]), (macro, g)
        assert x.encode() == store.get(g["contig_uid"], g["start"], g["end"], g["strand"], circular=is_circ)


@pytest.mark.covers("R-03")
def test_r03_every_gene_through_the_macros_equals_store_and_fasta(built, con):
    root, _, _ = built
    store = SeqStore(root)
    seqs = {cu: (seq, circ) for cu, seq, circ in _contigs()}
    genes = _genes(root)
    assert len(genes) > 500 and {g["strand"] for g in genes} == {1, -1}
    _check_genes(con, store, seqs, "seq_get", genes)
    _check_genes(con, store, seqs, "seq_get_circular", [g for g in genes if seqs[g["contig_uid"]][1]])
    # forward slices of the same genes through seq_slice
    con.register("g", _slices_table([(g["contig_uid"], g["start"], g["end"], 1) for g in genes]))
    got = con.execute("SELECT seq_slice(contig_uid, s, e) FROM g ORDER BY i").fetchall()
    con.unregister("g")
    assert [r[0].encode() for r in got] == [seqs[g["contig_uid"]][0][g["start"]:g["end"]] for g in genes]


@pytest.mark.covers("R-03", "C-05")
def test_r03_origin_wrapping_genes_through_seq_get_circular(built, tmp_path):
    # micro-S0 M. genitalium has no gene across its origin, so rotate it into one (as the C-05 gene-call test does)
    root, _, _ = built
    o = next(o for o in CONTRACT["objects"] if o["topology"] == "circular")
    seq = fasta.read(FIX / o["file"])[0].seq
    on_mg = [g for g in _genes(root) if g["contig_uid"] == ids.contig_uid(seq)]
    mid = on_mg[len(on_mg) // 2]
    r = (mid["start"] + mid["end"]) // 2
    rot = seq[r:] + seq[:r]
    src = tmp_path / "src"
    src.mkdir()
    (src / "rot.fna").write_bytes(b"\n".join(fasta.format_record("rot", rot, 70)) + b"\n")
    s = tmp_path / "s"
    ingest([SourceObject(src / "rot.fna", "rot", "1")], s, "S0", shard_bases=SHARD_BASES)
    cu = ids.contig_uid(rot)
    call_genes(s, "S0", [GenomeInput((cu,), o["genetic_code"], kind=o["kind"], circular=(cu,))])
    genes = _genes(s)
    wrap = [g for g in genes if g["end"] > len(rot)]
    assert len(wrap) >= 1
    c = catalog.open_catalog(catalog.generate_catalog(catalog.build_manifest(s, "x"), tmp_path / "c.duckdb"))
    try:
        _check_genes(c, SeqStore(s), {cu: (rot, True)}, "seq_get_circular", genes)
        _check_genes(c, SeqStore(s), {cu: (rot, True)}, "seq_get_circular", wrap)
        for g in wrap:  # strand -1 of a wrapping gene too, whatever the gene's own strand
            got = c.execute("SELECT seq_get_circular(?, ?, ?, ?)", [cu, g["start"], g["end"], -g["strand"]]).fetchone()[0]
            assert got.encode() == _oracle(rot, g["start"], g["end"], -g["strand"])
    finally:
        c.close()


@pytest.mark.covers("R-03")
@pytest.mark.parametrize("order", ["random", "descending", "one contig repeated"])
def test_r03_table_driven_fetch_equals_per_row_results_in_order(built, con, order):
    root, _, _ = built
    store = SeqStore(root)
    mg = next(c for c in _contigs() if c[2])[0]
    rows = _random_slices(1000, 20261109, only=mg if order == "one contig repeated" else None)
    if order == "descending":
        rows.sort(key=lambda r: (r[0], r[1], r[2]), reverse=True)
    assert len({r[0] for r in rows}) < len(rows)  # contigs repeat
    seqs = {cu: seq for cu, seq, _ in _contigs()}
    con.register("t", _slices_table(rows))
    got = con.execute("SELECT i, seq_get(contig_uid, s, e, strand) FROM t ORDER BY i").fetchall()
    con.unregister("t")
    assert [r[0] for r in got] == list(range(len(rows)))
    for k, (r, (_, x)) in enumerate(zip(rows, got)):
        assert x.encode() == _oracle(seqs[r[0]], r[1], r[2], r[3]) == store.get(*r), k
        if k % 10 == 0:  # the same slice asked alone
            assert con.execute("SELECT seq_get(?, ?, ?, ?)", list(r)).fetchone()[0] == x, k


@pytest.mark.covers("LY-06")
@pytest.mark.parametrize("n", [1, 1000])
def test_ly06_joined_fetch_scans_no_index_table(built, tmp_path, n):
    root, man, cat = built
    c = catalog.open_catalog(cat, local_dir=tmp_path / "local")
    try:
        # the macros hold no subquery on the index tables at all
        defs = dict(c.execute("SELECT function_name, macro_definition FROM duckdb_functions() "
                              "WHERE function_name IN ('seq_get', 'seq_slice', 'seq_get_circular', 'seq_raw')").fetchall())
        assert len(defs) == 4
        for name, sql in defs.items():
            assert not any(t in sql for t in ("contig_index", "seq_locator", "default_digest", "digest_len",
                                              "contig_len_of", "SELECT")), (name, sql)
        rows = _random_slices(n, 7)
        c.register("t", _slices_table(rows))
        sql = "SELECT sum(length(seq_get(contig_uid, s, e, strand))) FROM t"
        prof = catalog.profile_query(c, sql)
        assert prof["result"][0][0] == sum(r[2] - r[1] for r in rows)
        assert _index_scans(prof) == [], prof["scans"]
        touched = len({r[0] for r in rows})
        rg_total = {t: sum(f["row_groups"] for f in man["files"] if f["table"] == t) for t in ("seq_locator", "contig_index")}
        assert set(prof["index"]) == {"locator_row_groups_read", "index_row_groups_read"}
        assert 1 <= prof["index"]["locator_row_groups_read"] <= min(touched, rg_total["seq_locator"])
        assert 1 <= prof["index"]["index_row_groups_read"] <= min(touched, rg_total["contig_index"])
        again = catalog.profile_query(c, sql)  # every row group is in the lookup cache now
        assert again["result"] == prof["result"] and _index_scans(again) == []
        assert again["index"] == {"locator_row_groups_read": 0, "index_row_groups_read": 0}
        # the same through seq_raw (keys from a joined table of digests)
        c.register("r", pa.table({"rd": pa.array([ids.raw_digest(s) for _, s, _ in _contigs()], pa.binary(32))}))
        raw = catalog.profile_query(c, "SELECT sum(length(seq_raw(rd, 0, 100))) FROM r")
        assert raw["result"][0][0] == 300 and _index_scans(raw) == []
    finally:
        c.close()


@pytest.mark.covers("LY-06")
def test_ly06_a_batch_reads_blocks_in_file_order_whatever_the_row_order(built, record_property):
    """With a one-block cache, slices in random order read exactly as many blocks as the same
    slices given in (shard, raw_digest, start) order, which is file order since a shard's records
    are sorted by raw_digest; a direct loop in the random order reads more."""
    root, _, cat = built
    rows = [(cu, a, b, 1) for cu, a, b, _ in _random_slices(1000, 31, max_len=2000)]
    probe = SeqStore(root)
    key = {cu: (probe.shard_of(ids.raw_digest(seq)), ids.raw_digest(seq)) for cu, seq, _ in _contigs()}
    in_order = sorted(rows, key=lambda r: (*key[r[0]], r[1], r[2]))
    blocks = {}
    for name, rs in (("random", rows), ("sorted", in_order)):
        c = catalog.open_catalog(cat, cache_blocks=1)
        c.register("t", _slices_table(rs))
        prof = catalog.profile_query(c, "SELECT sum(length(seq_get(contig_uid, s, e, strand))) FROM t")
        c.close()
        assert prof["result"][0][0] == sum(r[2] - r[1] for r in rs)
        blocks[name] = prof["bgzf"]["blocks_read"]
    for name, rs in (("direct random", rows), ("direct sorted", in_order)):
        s = SeqStore(root, cache_blocks=1)
        for r in rs:
            s.get(*r)
        blocks[name] = s.counters()["blocks_read"]
    for name, v in blocks.items():
        record_property(f"blocks_read {name}", v)
    assert blocks["random"] == blocks["sorted"] == blocks["direct sorted"], blocks
    assert blocks["direct random"] > 2 * blocks["sorted"], blocks


def test_sequence_macros_return_null_for_a_null_argument(built, con):
    cu, seq, _ = next(c for c in _contigs() if c[2])
    rd = ids.raw_digest(seq)
    calls = {"seq_get": [cu, 0, 10, 1], "seq_get_circular": [cu, 0, 10, 1], "seq_slice": [cu, 0, 10], "seq_raw": [rd, 0, 10]}
    for macro, args in calls.items():
        marks = ", ".join("?" * len(args))
        for k in range(len(args)):
            a = list(args)
            a[k] = None
            assert con.execute(f"SELECT {macro}({marks})", a).fetchall() == [(None,)], (macro, k)
    # a NULL row among real rows: the others are still fetched, in order
    con.register("t", pa.table({"i": [0, 1, 2], "contig_uid": pa.array([cu, None, cu], pa.binary(16)),
                                "s": [0, 0, None], "e": [5, 5, 5], "strand": pa.array([1, 1, -1], pa.int8())}))
    got = con.execute("SELECT seq_get(contig_uid, s, e, strand) FROM t ORDER BY i").fetchall()
    con.unregister("t")
    assert got == [(seq[:5].decode(),), (None,), (None,)]


@pytest.mark.parametrize("sql, args, message", [
    ("SELECT seq_get(?, 0, 1, 1)", [b"\x00" * 16], "seq_get: unknown contig_uid"),
    ("SELECT seq_slice(?, 0, 1)", [b"\x00" * 16], "seq_get: unknown contig_uid"),
    ("SELECT seq_get_circular(?, 0, 1, 1)", [b"\x00" * 16], "seq_get_circular: unknown contig_uid"),
    ("SELECT seq_get(?, -1, 5, 1)", None, "seq_get: slice outside 0:contig_len"),
    ("SELECT seq_get(?, 5, 4, -1)", None, "seq_get: slice outside 0:contig_len"),
    ("SELECT seq_get(?, 0, 10000000, 1)", None, "seq_get: slice outside 0:contig_len"),
    ("SELECT seq_slice(?, 0, 10000000)", None, "seq_get: slice outside 0:contig_len"),
    ("SELECT seq_get(?, 0, 5, 0)", None, "seq_get: strand must be 1 or -1"),
    ("SELECT seq_get(?, 0, 5, 2)", None, "seq_get: strand must be 1 or -1"),
    ("SELECT seq_get_circular(?, 0, 5, 0)", None, "seq_get_circular: strand must be 1 or -1"),
    ("SELECT seq_get_circular(?, -1, 5, 1)", None, "seq_get_circular: slice outside 0:contig_len"),
    ("SELECT seq_get_circular(?, 10, 5, 1)", None, "seq_get_circular: slice outside 0:contig_len"),
    ("SELECT seq_get_circular(?, 580000, 2000000, 1)", None, "seq_get_circular: invalid wrapping slice"),
    ("SELECT seq_raw(?, 0, 1)", [b"\x00" * 32], "seq_raw: unknown raw_digest"),
    ("SELECT seq_raw(?, 0, 10000000)", "rd", "seq_raw: slice outside 0:contig_len"),
    ("SELECT seq_raw(?, 3, 2)", "rd", "seq_raw: slice outside 0:contig_len"),
])
def test_sequence_macro_errors_name_the_macro_and_the_fault(con, sql, args, message):
    _, seq, _ = next(c for c in _contigs() if c[2])  # M. genitalium, 580,076 bp
    if args is None:
        args = [ids.contig_uid(seq)]
    elif args == "rd":
        args = [ids.raw_digest(seq)]
    with pytest.raises(duckdb.Error, match=message):
        con.execute(sql, args).fetchall()


def test_a_bad_row_fails_the_whole_table_query(con):
    cu, seq, _ = _contigs()[0]
    con.register("t", _slices_table([(cu, 0, 5, 1), (cu, 0, len(seq) + 1, 1), (cu, 5, 9, -1)]))
    with pytest.raises(duckdb.Error, match="seq_get: slice outside 0:contig_len"):
        con.execute("SELECT seq_get(contig_uid, s, e, strand) FROM t").fetchall()
    con.unregister("t")


def test_protein_and_slot_macros(built, con):
    root, _, _ = built
    occ = pq.read_table(root / f"l1/protein_occurrence/gene_call_version={GENE_CALL_VERSION}/stage=S0/part-0.parquet").to_pylist()
    puid = occ[len(occ) // 2]["protein_uid"]
    want = sorted((r["contig_uid"], r["feature_uid"]) for r in occ if r["protein_uid"] == puid)
    got = sorted(con.execute("SELECT contig_uid, feature_uid FROM slots_of(?)", [puid]).fetchall())
    assert got == want
    prot = {}
    for f in root.glob("l0seq/protein_seq/**/part-0.parquet"):
        prot.update({r["protein_uid"]: r["seq"] for r in pq.read_table(f).to_pylist()})
    assert con.execute("SELECT protein_of(?)", [puid]).fetchone()[0] == prot[puid]


# --- atomic swap --------------------------------------------------------------------


def test_swap_current_is_a_rename(built, tmp_path):
    _, man, _ = built
    d = tmp_path / "manifests"
    for r in ("rS0.1", "rS0.2"):
        catalog.write_manifest({**man, "release_id": r}, d / f"release={r}.json")
    catalog.swap_current(d, "rS0.1")
    assert catalog.read_manifest(d / "current.json")["release_id"] == "rS0.1"
    catalog.swap_current(d, "rS0.2")
    assert (d / "current.json").is_symlink()
    assert catalog.read_manifest(d / "current.json")["release_id"] == "rS0.2"
    assert sorted(p.name for p in d.iterdir()) == ["current.json", "release=rS0.1.json", "release=rS0.2.json"]
    with pytest.raises(FileNotFoundError):
        catalog.swap_current(d, "rS0.9")


def test_regenerating_the_catalog_replaces_it_whole(built, tmp_path):
    _, man, cat = built
    out = tmp_path / "c.duckdb"
    catalog.generate_catalog(man, out)
    first = catalog.open_catalog(out)
    n = first.execute("SELECT count(*) FROM contig_index").fetchall()[0][0]
    catalog.generate_catalog({**man, "release_id": "rS0.next"}, out)
    # an open reader keeps its snapshot; a new connection sees the new catalog
    assert first.execute("SELECT count(*) FROM contig_index").fetchall()[0][0] == n
    first.close()
    second = catalog.open_catalog(out)
    assert dict(second.execute("SELECT * FROM catalog_info").fetchall())["release_id"] == "rS0.next"
    second.close()
    assert sorted(p.name for p in tmp_path.iterdir() if p.name.startswith(".")) == []


# --- LY-09: a nullable column added later --------------------------------------------------


def test_ly09_catalog_reads_files_from_before_a_nullable_column(built, tmp_path, monkeypatch):
    import dataclasses

    from biobox.registry import schemas

    root, man, _ = built
    spec = schemas.TABLES["contig_index"]
    wider = dataclasses.replace(spec, schema=spec.schema.append(pa.field("topology", pa.string(), nullable=True)))
    monkeypatch.setitem(schemas.TABLES, "contig_index", wider)
    store = tmp_path / "s"
    shutil.copytree(root, store, ignore=shutil.ignore_patterns("catalog", "manifests"))
    old = next(store.glob("l0seq/contig_index/stage=S0/part-0.parquet"))
    # a later stage writes the wider schema
    t = pq.read_table(old).append_column("topology", pa.array(["linear"] * pq.read_table(old).num_rows, pa.string()))
    new = store / "l0seq/contig_index/stage=S1/part-0.parquet"
    new.parent.mkdir(parents=True)
    pq.write_table(t.cast(wider.schema), new)
    m = catalog.build_manifest(store, "x")
    c = catalog.open_catalog(catalog.generate_catalog(m, tmp_path / "c.duckdb"))
    rows = c.execute("SELECT stage, topology FROM contig_index ORDER BY stage").fetchall()
    assert {r for r in rows} == {("S0", None), ("S1", "linear")}
    c.close()
    # with only old files, the new column is a typed NULL
    new.unlink()
    m = catalog.build_manifest(store, "x")
    c = catalog.open_catalog(catalog.generate_catalog(m, tmp_path / "c2.duckdb"))
    assert c.execute("SELECT count(*) FROM contig_index WHERE topology IS NULL").fetchall()[0][0] > 0
    c.close()
    # a retyped column is still refused
    bad = pq.read_table(old)
    bad = bad.set_column(bad.column_names.index("contig_len"), "contig_len", bad.column("contig_len").cast(pa.int32()))
    pq.write_table(bad, old)
    with pytest.raises(catalog.LayoutError, match="C-01"):
        catalog.build_manifest(store, "x")
