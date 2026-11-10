"""Ring R2: micro-S0 fixture through O-INGEST and O-SLICE.

Covers R-01, R-02, R-03 (data_layer_and_tests.md §2A.5), C-01/C-02/C-06 via the
registry contract checker, LY-07 and LY-08 (TEST_ARCHITECTURE §7), and
MR-01, MR-02, MR-03, MR-14 (TEST_ARCHITECTURE §4).
"""

from __future__ import annotations

import gzip
import hashlib
import random
import shutil
from pathlib import Path

import duckdb
import pyarrow.parquet as pq
import pytest
import yaml

from biobox import fasta, ids
from biobox.l0seq import SeqStore, SourceObject, ingest, revcomp
from biobox.registry.schemas import conformance_problems, table

pytestmark = pytest.mark.r2

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())


def _objects(paths):
    return [
        SourceObject(Path(p), CONTRACT["source_name"], CONTRACT["source_release"], CONTRACT["license"], True)
        for p in paths
    ]


def _fixture_paths():
    return [FIX / o["file"] for o in CONTRACT["objects"]]


def _files(root: Path) -> dict[str, bytes]:
    """Every data file: Parquet tables and the bgzip shards with their indexes."""
    return {str(p.relative_to(root)): p.read_bytes() for p in sorted(root.rglob("*"))
            if p.is_file() and (p.suffix == ".parquet" or "/seq_shards/" in p.as_posix())}


@pytest.fixture(scope="module")
def store_root(tmp_path_factory):
    root = tmp_path_factory.mktemp("micro_s0")
    # small shards so the three genomes span several shards and slices cross BGZF blocks
    ingest(_objects(_fixture_paths()), root, stage="S0", shard_bases=200_000)
    return root


# --- contract -----------------------------------------------------------------


def test_fixture_matches_contract():
    for o in CONTRACT["objects"]:
        p = FIX / o["file"]
        assert ids.sha256_file(p).hex() == o["sha256_gz"]
        with gzip.open(p, "rb") as fh:
            assert hashlib.sha256(fh.read()).hexdigest() == o["sha256_fna"]
        recs = fasta.read(p)
        assert [r.accession for r in recs] == [o["accession"]]
        assert len(recs[0].seq) == o["length"]


@pytest.mark.covers("C-01", "C-02", "C-06", "LY-01")
@pytest.mark.parametrize(
    "name,glob",
    [
        ("contig_occurrence", "l0seq/contig_occurrence/stage=*/source=*/part-0.parquet"),
        ("contig_index", "l0seq/contig_index/stage=*/part-0.parquet"),
        ("seq_locator", "l0seq/seq_locator/stage=*/part-0.parquet"),
        ("l0_manifest", "l0/l0_manifest.parquet"),
    ],
)
def test_written_tables_conform(store_root, name, glob):
    files = sorted(store_root.glob(glob))
    assert files, f"no files for {name}"
    for f in files:
        assert conformance_problems(pq.read_table(f), table(name)) == [], f


# --- R-01, R-02, R-03 -----------------------------------------------------------


def test_r01_record_digests(store_root):
    store = SeqStore(store_root)
    for p in _fixture_paths():
        for rec in fasta.read(p):
            rd = ids.raw_digest(rec.seq)
            assert ids.raw_digest(store.raw(rd)) == rd
            assert store.get(ids.contig_uid(rec.seq)) == rec.seq


def test_r02_object_round_trip(store_root):
    store = SeqStore(store_root)
    for p in _fixture_paths():
        expected = hashlib.sha256(b"\n".join(fasta.normalized_lines(p)) + b"\n").digest()
        got = store.reconstruct_object(ids.sha256_file(p))
        assert hashlib.sha256(b"\n".join(got) + b"\n").digest() == expected, p.name


def test_r03_random_slices_both_strands(store_root):
    store = SeqStore(store_root)
    rng = random.Random(20261029)
    for p in _fixture_paths():
        seq = fasta.read(p)[0].seq
        uid = ids.contig_uid(seq)
        for _ in range(200):
            a, b = sorted(rng.randrange(0, len(seq) + 1) for _ in range(2))
            assert store.get(uid, a, b, 1) == seq[a:b]
            assert store.get(uid, a, b, -1) == revcomp(seq[a:b])
        assert store.get(uid, len(seq) - 1, len(seq)) == seq[-1:]


def test_slice_out_of_bounds_raises(store_root):
    store = SeqStore(store_root)
    seq = fasta.read(_fixture_paths()[0])[0].seq
    with pytest.raises(ValueError):
        store.get(ids.contig_uid(seq), 0, len(seq) + 1)


# --- layout tests ----------------------------------------------------------------


def test_ly07_duckdb_reads_the_same_rows(store_root):
    files = sorted(str(f) for f in store_root.glob("l0seq/contig_occurrence/stage=*/source=*/part-0.parquet"))
    con = duckdb.connect()
    n_duck, bytes_duck = con.execute(
        "SELECT count(*), sum(length(original_header)) FROM read_parquet(?)", [files]
    ).fetchall()[0]
    arrow = [pq.read_table(f) for f in files]
    assert n_duck == sum(t.num_rows for t in arrow)
    assert bytes_duck == sum(sum(len(s) for s in t.column("original_header").to_pylist()) for t in arrow)
    typ = con.execute("SELECT typeof(raw_digest) FROM read_parquet(?) LIMIT 1", [files]).fetchall()[0][0]
    assert typ == "BLOB"


def test_shards_hold_each_sequence_once_sorted_by_digest(store_root):
    """Each shard is sorted by digest; there is no order across shards (card B-15, TASKS §2)."""
    from biobox.l0seq.bgzf import read_fai

    shards = sorted(store_root.glob("l0seq/seq_shards/stage=S0/shard-*.fa.gz"))
    assert len(shards) > 1  # shard_bases=200_000 splits the micro set
    per_shard = [[e.name for e in read_fai(str(s) + ".fai")] for s in shards]
    for names in per_shard:
        assert names == sorted(names)
    names = sorted(n for ns in per_shard for n in ns)
    want = sorted(ids.raw_digest(r.seq).hex() for p in _fixture_paths() for r in fasta.read(p))
    assert names == want


def test_ly08_rewrite_is_byte_identical(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    for r in (a, b):
        ingest(_objects(_fixture_paths()), r, stage="S0", shard_bases=200_000)
    assert _files(a) == _files(b)


# --- metamorphic relations -------------------------------------------------------


def test_mr01_input_order_does_not_matter(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    ingest(_objects(_fixture_paths()), a, stage="S0", shard_bases=200_000)
    ingest(_objects(list(reversed(_fixture_paths()))), b, stage="S0", shard_bases=200_000)
    assert _files(a) == _files(b)


def _write_fasta(path: Path, header: str, seq: bytes, width: int = 70) -> Path:
    path.write_bytes(b"\n".join(fasta.format_record(header, seq, width)) + b"\n")
    return path


def test_mr02_renamed_header_keeps_uids_and_chunks(tmp_path):
    lam = fasta.read(_fixture_paths()[0])[0]
    src = tmp_path / "src"
    src.mkdir()
    renamed = _write_fasta(src / "renamed.fna", "XX_000001.1 renamed lambda", lam.seq)
    a, b = tmp_path / "a", tmp_path / "b"
    ingest(_objects([_fixture_paths()[0]]), a, stage="S0", shard_bases=200_000)
    ingest(_objects([renamed]), b, stage="S0", shard_bases=200_000)
    fa, fb = _files(a), _files(b)
    for k in fa:
        if "seq_shards" in k or "contig_index" in k or "seq_locator" in k:
            assert fa[k] == fb[k], k
    occ_a = pq.read_table(next(a.glob("l0seq/contig_occurrence/**/*.parquet"))).to_pylist()[0]
    occ_b = pq.read_table(next(b.glob("l0seq/contig_occurrence/**/*.parquet"))).to_pylist()[0]
    assert occ_a["contig_uid"] == occ_b["contig_uid"]
    assert occ_a["accession"] != occ_b["accession"]


def test_mr03_reverse_complement_slices_mirror(tmp_path):
    lam = fasta.read(_fixture_paths()[0])[0]
    src = tmp_path / "src"
    src.mkdir()
    rc_path = _write_fasta(src / "lambda_rc.fna", "lambda_rc", revcomp(lam.seq))
    root = tmp_path / "s"
    ingest(_objects([_fixture_paths()[0], rc_path]), root, stage="S0", shard_bases=200_000)
    store = SeqStore(root)
    fwd, rc = ids.contig_uid(lam.seq), ids.contig_uid(revcomp(lam.seq))
    n = len(lam.seq)
    rng = random.Random(3)
    for _ in range(200):
        a, b = sorted(rng.randrange(0, n + 1) for _ in range(2))
        assert store.get(rc, n - b, n - a, -1) == store.get(fwd, a, b, 1)


def test_mr14_soft_mask_keeps_contig_uid(tmp_path):
    lam = fasta.read(_fixture_paths()[0])[0]
    masked = lam.seq[:1000].lower() + lam.seq[1000:]
    src = tmp_path / "src"
    src.mkdir()
    mpath = _write_fasta(src / "lambda_masked.fna", "lambda_masked", masked)
    root = tmp_path / "s"
    objs = _objects([_fixture_paths()[0]]) + [
        SourceObject(mpath, "other_source", "1", "INSDC (public)", True, precedence=200)
    ]
    ingest(objs, root, stage="S0", shard_bases=200_000)
    store = SeqStore(root)
    assert ids.contig_uid(masked) == ids.contig_uid(lam.seq)
    assert ids.raw_digest(masked) != ids.raw_digest(lam.seq)
    assert store.raw(ids.raw_digest(masked)) == masked
    assert store.get(ids.contig_uid(lam.seq)) == lam.seq  # default variant: higher-precedence source


def test_identical_sequence_stored_once_across_sources_and_stages(tmp_path):
    root = tmp_path / "s"
    ingest(_objects([_fixture_paths()[0]]), root, stage="S0", shard_bases=200_000)
    before = _files(root)
    copy_dir = tmp_path / "copy"
    copy_dir.mkdir()
    dup = copy_dir / "lambda_again.fna.gz"
    shutil.copy(_fixture_paths()[0], dup)
    report = ingest(
        [SourceObject(dup, "second_source", "1", "INSDC (public)", True)], root, stage="S1", shard_bases=200_000
    )
    assert report["new_sequences"] == 0 and report["index_rows"] == 0
    after = _files(root)
    for k, v in before.items():
        if "seq_shards" in k or "contig_index" in k or "contig_occurrence" in k:
            assert after[k] == v, f"S1 rewrote an S0 file: {k}"
    assert not list(root.glob("l0seq/seq_shards/stage=S1/*"))
    assert list(root.glob("l0seq/contig_occurrence/stage=S1/source=second_source/part-0.parquet"))


# --- the index is looked up, not loaded (card B-13) --------------------------------


@pytest.mark.covers("LY-06", "R-03")
def test_store_open_reads_no_row_group_and_one_get_reads_one_of_each_index(store_root):
    store = SeqStore(store_root)
    assert store.locator.row_groups_read == 0 and store.index.row_groups_read == 0
    assert store.locator.footers_read == 1 and store.index.footers_read == 1  # one stage, one file each
    seq = fasta.read(_fixture_paths()[1])[0].seq
    cu, rd = ids.contig_uid(seq), ids.raw_digest(seq)
    assert store.get(cu, 1000, 3000, -1) == revcomp(seq[1000:3000])
    assert store.locator.row_groups_read == 1 and store.index.row_groups_read == 1
    assert store.default_digest(cu) == rd and store.length_of(rd) == len(seq)
    assert store.shard_of(rd).startswith("l0seq/seq_shards/stage=S0/shard-")
    assert store.locator.row_groups_read == 1 and store.index.row_groups_read == 1  # served by the cache


def test_unknown_keys_raise_key_error(store_root):
    store = SeqStore(store_root)
    with pytest.raises(KeyError):
        store.get(b"\x00" * 16)
    with pytest.raises(KeyError):
        store.default_digest(b"\xff" * 16)
    with pytest.raises(KeyError):
        store.raw(b"\x00" * 32)
    with pytest.raises(KeyError):
        store.length_of(b"\xff" * 32)


@pytest.mark.covers("R-09")
def test_second_stage_adds_locator_files_and_keeps_the_first_byte_identical(tmp_path):
    lam = fasta.read(_fixture_paths()[0])[0]
    other = fasta.read(_fixture_paths()[1])[0]
    root = tmp_path / "s"
    ingest(_objects([_fixture_paths()[0]]), root, stage="S0", shard_bases=200_000)
    s0 = root / "l0seq/seq_locator/stage=S0/part-0.parquet"
    before = s0.read_bytes()
    assert pq.read_table(s0).column("raw_digest").to_pylist() == [ids.raw_digest(lam.seq)]
    src = tmp_path / "src"
    src.mkdir()
    rc_path = _write_fasta(src / "lambda_rc.fna", "lambda_rc", revcomp(lam.seq))
    ingest(_objects([rc_path, _fixture_paths()[1], _fixture_paths()[0]]), root, stage="S1", shard_bases=200_000)
    assert s0.read_bytes() == before
    assert sorted(p.relative_to(root).as_posix() for p in (root / "l0seq").rglob("seq_locator*")) == [
        "l0seq/seq_locator"]
    s1 = root / "l0seq/seq_locator/stage=S1/part-0.parquet"
    assert sorted(p.relative_to(root / "l0seq/seq_locator").as_posix() for p in s0.parents[1].rglob("*.parquet")) == [
        "stage=S0/part-0.parquet", "stage=S1/part-0.parquet"]
    t1 = pq.read_table(s1)
    # only the sequences first stored in S1, sorted by raw_digest
    assert t1.column("raw_digest").to_pylist() == sorted([ids.raw_digest(revcomp(lam.seq)), ids.raw_digest(other.seq)])
    assert conformance_problems(t1, table("seq_locator")) == []
    md = pq.read_metadata(s1)
    assert all(md.row_group(g).column(0).statistics.has_min_max for g in range(md.num_row_groups))
    store = SeqStore(root)
    assert store.locator.footers_read == 2
    for s in (lam.seq, revcomp(lam.seq), other.seq):
        assert store.get(ids.contig_uid(s), 5, 905) == s[5:905]
    assert store.shard_of(ids.raw_digest(lam.seq)).startswith("l0seq/seq_shards/stage=S0/")
    assert store.shard_of(ids.raw_digest(other.seq)).startswith("l0seq/seq_shards/stage=S1/")


def _many_records(n: int, seed: int) -> list[tuple[str, bytes]]:
    rng = random.Random(seed)
    return sorted((hashlib.sha256(str(i).encode()).hexdigest(),
                   bytes(rng.choice(b"ACGT") for _ in range(rng.randrange(1, 300)))) for i in range(n))


@pytest.mark.covers("LY-06")
def test_fai_lookup_touches_about_log2_n_lines(tmp_path):
    import math

    from biobox.l0seq.bgzf import BgzfFasta, read_fai, write_bgzf_fasta

    n = 4096
    recs = _many_records(n, 9)
    out = tmp_path / "s.fa.gz"
    write_bgzf_fasta(recs, out)
    bg = BgzfFasta(out)
    assert bg.fai_lines_read == 0  # opening reads no line of the index
    bound = math.ceil(math.log2(n)) + 3
    rng = random.Random(10)
    worst = 0
    for name, seq in rng.sample(recs, 300):
        before = bg.fai_lines_read
        assert bg.fetch(name, 0, len(seq)) == seq
        worst = max(worst, bg.fai_lines_read - before)
    assert 0 < worst <= bound, worst
    before = bg.fai_lines_read
    assert bg.fetch(name, 0, 1) == seq[:1]  # a repeat is served by the entry cache
    assert bg.fai_lines_read == before
    full = read_fai(str(out) + ".fai")  # oracle: the whole index, read line by line
    for e in [full[0], full[-1], *rng.sample(full, 50)]:
        assert bg.entry(e.name) == e
    for absent in ("0" * 64, "f" * 65, "", recs[7][0][:-1], recs[7][0] + "0"):
        with pytest.raises(KeyError):
            bg.entry(absent)
    assert bg.check_sorted()


@pytest.mark.covers("LY-02")
def test_unsorted_fai_is_reported(tmp_path):
    from biobox.l0seq.bgzf import BgzfFasta, write_bgzf_fasta

    recs = _many_records(200, 11)
    out = tmp_path / "s.fa.gz"
    write_bgzf_fasta(recs, out)
    fai = Path(str(out) + ".fai")
    fai.write_text("".join(fai.read_text().splitlines(keepends=True)[::-1]))
    bg = BgzfFasta(out)
    assert not bg.check_sorted()
    raised = 0
    for name, _ in recs[:40]:
        try:
            e = bg.entry(name)
        except ValueError as err:
            assert "not sorted by name" in str(err)
            raised += 1
        except KeyError:
            pass
        else:
            assert e.name == name
    assert raised > 0


# --- streaming FASTA reading (card B-15) ------------------------------------------------


def _oracle_records(path: Path) -> list[fasta.Record]:
    """An independent whole-file parse: CRLF folded, blank lines dropped, width 0 if irregular."""
    data = Path(path).read_bytes()
    if str(path).endswith(".gz"):
        data = gzip.decompress(data)
    lines = [ln.rstrip(b"\r") for ln in data.split(b"\n")]
    parsed: list[list] = []
    for ln in lines:
        if not ln.strip():
            continue
        if ln.startswith(b">"):
            parsed.append([ln[1:].decode(), []])
        else:
            parsed[-1][1].append(ln)
    out = []
    for header, ls in parsed:
        w = len(ls[0]) if ls else 0
        regular = bool(ls) and all(len(x) == w for x in ls[:-1]) and len(ls[-1]) <= w
        out.append(fasta.Record(header, b"".join(ls), w if regular else 0))
    return out


_FASTA_CASES = {
    "crlf_and_blank_lines.fna": (
        b">a desc\r\nACGT\r\nAC\r\n\r\n>b\r\n\r\nGGGG\r\nGG\r\n",
        [("a desc", b"ACGTAC", 4), ("b", b"GGGGGG", 4)]),
    "irregular.fna": (
        b">x\nACG\nACGTA\nA\n>y\nAC\nACG\n>z\nACGT\nAC\n>s\nACGTACGT\n",
        [("x", b"ACGACGTAA", 0), ("y", b"ACACG", 0), ("z", b"ACGTAC", 4), ("s", b"ACGTACGT", 8)]),
    "empty_record_no_final_newline.fna": (
        b">e\n>f\nA\n>m\nacgtNRY\n>n\nAC\nG",
        [("e", b"", 0), ("f", b"A", 1), ("m", b"acgtNRY", 7), ("n", b"ACG", 2)]),
}


@pytest.mark.covers("R-01", "R-02")
@pytest.mark.parametrize("gz", [False, True])
@pytest.mark.parametrize("name", sorted(_FASTA_CASES))
def test_iter_records_keeps_exact_bytes_and_line_width(tmp_path, name, gz):
    data, want = _FASTA_CASES[name]
    path = tmp_path / (name + (".gz" if gz else ""))
    path.write_bytes(gzip.compress(data, mtime=0) if gz else data)
    it = fasta.iter_records(path)
    assert not isinstance(it, list)
    got = list(it)
    assert [(r.header, r.seq, r.line_width) for r in got] == want  # constructed truth
    assert got == _oracle_records(path) == fasta.read(path)


@pytest.mark.covers("R-01")
def test_iter_records_on_the_fixtures_equals_a_whole_file_parse():
    for o in CONTRACT["objects"]:
        p = FIX / o["file"]
        got = list(fasta.iter_records(p))
        assert got == _oracle_records(p) == fasta.read(p)
        assert [len(r.seq) for r in got] == [o["length"]] and got[0].line_width > 0


def test_iter_records_streams_and_refuses_sequence_before_a_header(tmp_path):
    lazy = tmp_path / "lazy.fna"
    lazy.write_bytes(b">ok\nAC\n>\xff\nGG\n")  # the second header is not UTF-8
    assert next(fasta.iter_records(lazy)) == fasta.Record("ok", b"AC", 2)  # yielded before the bad line is read
    with pytest.raises(UnicodeDecodeError):
        fasta.read(lazy)
    bad = tmp_path / "bad.fna"
    bad.write_bytes(b"ACGT\n>a\nAC\n")
    with pytest.raises(ValueError, match="sequence before the first header"):
        next(fasta.iter_records(bad))


# --- batch ingest (card B-15) -----------------------------------------------------------


def _all_files(root: Path) -> dict[str, bytes]:
    return {p.relative_to(root).as_posix(): p.read_bytes() for p in sorted(root.rglob("*")) if p.is_file()}


def _rand(rng: random.Random, n: int) -> bytes:
    return bytes(b"ACGT"[x & 3] for x in rng.randbytes(n))


def _write_records(path: Path, recs: list[tuple[str, bytes]], width: int = 70) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"".join(b"\n".join(fasta.format_record(h, s, width)) + b"\n" for h, s in recs))
    return path


def _world(d: Path) -> list[SourceObject]:
    """The fixture plus five multi-record objects with duplicates within and across objects, a copy of
    lambda, and a soft-masked variant from a lower-precedence source."""
    rng = random.Random(20261109)
    a, b, c, dd, e, f, g = (_rand(rng, n) for n in (30_000, 45_000, 25_000, 60_000, 20_000, 35_000, 50_000))
    lam = fasta.read(_fixture_paths()[0])[0].seq
    synth = [
        _write_records(d / "s0.fna", [("a1", a), ("b1", b)]),
        _write_records(d / "s1.fna", [("c1", c), ("a2 again", a)]),
        _write_records(d / "s2.fna", [("d1", dd), ("d2 same", dd), ("lam copy", lam)], width=60),
        _write_records(d / "s3.fna", [("e1", e), ("f1", f), ("g1", g)]),
    ]
    masked = _write_records(d / "s4.fna", [("b masked", b[:500].lower() + b[500:])])
    return (_objects(_fixture_paths()) + _objects(synth)
            + [SourceObject(masked, "other_source", "1", "INSDC (public)", True, precedence=200)])


def _sha(o: SourceObject) -> bytes:
    return hashlib.sha256(Path(o.path).read_bytes()).digest()


def _expected_batches(objs: list[SourceObject], shard_bases: int, known: set[bytes] = frozenset()
                      ) -> list[dict[bytes, bytes]]:
    """The batch rule, computed independently: objects by sha256; a digest is stored by its first
    occurrence unless an earlier stage holds it; a batch is a maximal run of objects whose new bases
    total at most shard_bases (one larger object is a batch of its own)."""
    seen = set(known)
    batches: list[dict[bytes, bytes]] = []
    held = 0
    open_batch = False
    for o in sorted(objs, key=_sha):
        new = {}
        for r in _oracle_records(o.path):
            rd = hashlib.sha256(r.seq).digest()
            if rd not in seen:
                seen.add(rd)
                new[rd] = r.seq
        if not new:
            continue
        n = sum(len(s) for s in new.values())
        if not open_batch or held + n > shard_bases:
            batches.append({})
            held = 0
            open_batch = True
        batches[-1].update(new)
        held += n
    return batches


def _expected_index(objs: list[SourceObject]) -> list[tuple]:
    """contig_index rows of a first stage, by brute force: one row per (contig_uid, raw_digest);
    the default is the variant with the lowest (precedence, object sha256, record ordinal)."""
    best: dict[bytes, tuple] = {}
    pairs: dict[tuple, int] = {}
    for o in objs:
        for k, r in enumerate(_oracle_records(o.path)):
            cu, rd = ids.contig_uid(r.seq), hashlib.sha256(r.seq).digest()
            pairs[(cu, rd)] = len(r.seq)
            key = (o.precedence, _sha(o), k)
            if cu not in best or key < best[cu][0]:
                best[cu] = (key, rd)
    return [(cu, rd, n, best[cu][1] == rd) for (cu, rd), n in sorted(pairs.items())]


SMALL_SHARDS = 100_000


@pytest.mark.covers("MR-01", "LY-08")
def test_mr01_shuffled_objects_give_byte_identical_files(tmp_path):
    objs = _world(tmp_path / "src")
    later = _write_records(tmp_path / "src2" / "later.fna", [("later", _rand(random.Random(1), 9_000))])
    s1 = _objects([later]) + objs[:2] + objs[-2:]  # one new object, four whose sequences S0 holds
    results = []
    for seed in (None, 1, 2, 3):
        order, order1 = list(objs), list(s1)
        if seed is not None:
            random.Random(seed).shuffle(order)
            random.Random(seed).shuffle(order1)
        root = tmp_path / f"r{seed}"
        ingest(order, root, stage="S0", shard_bases=SMALL_SHARDS)
        ingest(order1, root, stage="S1", shard_bases=SMALL_SHARDS)
        results.append(_all_files(root))
    assert len(results[0]) > 20
    for r in results[1:]:
        assert r.keys() == results[0].keys()
        for k in r:
            assert r[k] == results[0][k], k


@pytest.mark.covers("R-02", "R-03", "LY-02")
def test_shards_follow_the_batch_rule(tmp_path):
    from biobox.l0seq.bgzf import read_fai

    objs = _world(tmp_path / "src")
    root = tmp_path / "s"
    report = ingest(objs, root, stage="S0", shard_bases=SMALL_SHARDS)
    batches = _expected_batches(objs, SMALL_SHARDS)
    assert len(batches) >= 5 and any(len(b) > 2 for b in batches)  # several shards, some of several objects
    shards = sorted(p.relative_to(root).as_posix() for p in root.glob("l0seq/seq_shards/stage=S0/shard-*.fa.gz"))
    assert shards == [f"l0seq/seq_shards/stage=S0/shard-{i:04d}.fa.gz" for i in range(len(batches))]
    for shard, batch in zip(shards, batches):
        fai = read_fai(root / (shard + ".fai"))
        assert [e.name for e in fai] == sorted(rd.hex() for rd in batch)  # sorted by digest inside the shard
        assert [e.length for e in fai] == [len(batch[rd]) for rd in sorted(batch)]
    loc = pq.read_table(root / "l0seq/seq_locator/stage=S0/part-0.parquet").to_pylist()
    want_loc = sorted((rd, shard, len(s)) for shard, b in zip(shards, batches) for rd, s in b.items())
    assert [(r["raw_digest"], r["shard"], r["contig_len"]) for r in loc] == want_loc
    assert report["new_sequences"] == len(want_loc)
    store = SeqStore(root)
    for b in batches:
        for rd, s in b.items():
            assert store.raw(rd) == s
    for o in objs:  # R-02: every source object is rebuilt from occurrence rows and shards
        expected = hashlib.sha256(b"\n".join(fasta.normalized_lines(o.path)) + b"\n").digest()
        got = store.reconstruct_object(_sha(o))
        assert hashlib.sha256(b"\n".join(got) + b"\n").digest() == expected, o.path
    idx = pq.read_table(root / "l0seq/contig_index/stage=S0/part-0.parquet").to_pylist()
    assert [(r["contig_uid"], r["raw_digest"], r["contig_len"], r["is_default_variant"]) for r in idx] == \
        _expected_index(objs)
    occ = [r for f in sorted(root.glob("l0seq/contig_occurrence/stage=S0/source=*/part-0.parquet"))
           for r in pq.read_table(f).to_pylist()]
    assert len(occ) == sum(len(_oracle_records(o.path)) for o in objs)


@pytest.mark.covers("R-09")
def test_duplicates_are_stored_once_in_the_smaller_sha_batch_and_not_again_later(tmp_path):
    from biobox import ledger
    from biobox.l0seq.bgzf import read_fai

    rng = random.Random(5)
    s, u, v, t = (_rand(rng, n) for n in (40_000, 70_000, 70_000, 30_000))
    a = _write_records(tmp_path / "src/a.fna", [("s in a", s), ("u", u)])
    b = _write_records(tmp_path / "src/b.fna", [("v", v), ("s in b", s)])
    objs = _objects([a, b])
    lo, hi = sorted(objs, key=_sha)
    root = tmp_path / "s"
    ingest(objs, root, stage="S0", shard_bases=SMALL_SHARDS)
    store = SeqStore(root)
    rd = ids.raw_digest
    first_lo = rd(u) if lo.path == a else rd(v)  # lo's own sequence
    assert store.shard_of(rd(s)) == store.shard_of(first_lo) == "l0seq/seq_shards/stage=S0/shard-0000.fa.gz"
    names = [e.name for sh in sorted(root.glob("l0seq/seq_shards/stage=S0/*.fa.gz")) for e in read_fai(str(sh) + ".fai")]
    assert sorted(names) == sorted({rd(x).hex() for x in (s, u, v)})  # s once
    rows = {r["object_key"]: (r["n_in"], r["n_kept"], r["n_merged"]) for r in ledger.read(root, "S0", "ingest")}
    assert rows == {_sha(lo): (2, 2, 0), _sha(hi): (2, 1, 1)}

    before = _all_files(root)
    c = _write_records(tmp_path / "src/c.fna", [("s again", s), ("t", t), ("u again", u)])
    report = ingest(_objects([c]), root, stage="S1", shard_bases=SMALL_SHARDS)
    after = _all_files(root)
    for k, val in before.items():
        if k != "l0/l0_manifest.parquet":  # rewritten per stage (a finding of card B-15, not changed here)
            assert after[k] == val, f"S1 rewrote {k}"
    assert report["new_sequences"] == 1
    s1 = sorted(root.glob("l0seq/seq_shards/stage=S1/*.fa.gz"))
    assert [sh.name for sh in s1] == ["shard-0000.fa.gz"]
    assert [e.name for e in read_fai(str(s1[0]) + ".fai")] == [rd(t).hex()]
    loc1 = pq.read_table(root / "l0seq/seq_locator/stage=S1/part-0.parquet").to_pylist()
    assert [(r["raw_digest"], r["shard"]) for r in loc1] == [(rd(t), "l0seq/seq_shards/stage=S1/shard-0000.fa.gz")]
    idx1 = pq.read_table(root / "l0seq/contig_index/stage=S1/part-0.parquet").to_pylist()
    assert [(r["contig_uid"], r["is_default_variant"]) for r in idx1] == [(ids.contig_uid(t), True)]
    assert {r["object_key"]: (r["n_in"], r["n_kept"], r["n_merged"]) for r in ledger.read(root, "S1", "ingest")} == {
        _sha(_objects([c])[0]): (3, 1, 2)}
    store = SeqStore(root)
    for x in (s, u, v, t):
        assert store.get(ids.contig_uid(x)) == x


def test_memory_holds_one_batch_not_the_stage(tmp_path, monkeypatch):
    """Sequence bytes alive at once, counted on the records the FASTA reader hands out, stay under
    shard_bases plus the largest object; so does every batch handed to the shard writer."""
    import biobox.l0seq.store as store_mod

    objs = _world(tmp_path / "src")
    live = {"now": 0, "max": 0}

    class Tracked(bytes):
        def __del__(self):
            live["now"] -= len(self)

    def track(r: fasta.Record) -> fasta.Record:
        if isinstance(r.seq, Tracked):
            return r
        seq = Tracked(r.seq)
        live["now"] += len(seq)
        live["max"] = max(live["max"], live["now"])
        return fasta.Record(r.header, seq, r.line_width)

    orig_read, orig_iter = fasta.read, getattr(fasta, "iter_records", None)
    monkeypatch.setattr(fasta, "read", lambda p: [track(r) for r in orig_read(p)])
    if orig_iter is not None:
        monkeypatch.setattr(fasta, "iter_records", lambda p: (track(r) for r in orig_iter(p)))
    handed = []
    orig_write = store_mod.write_bgzf_fasta

    def write(records, out, *args, **kw):
        records = list(records)
        handed.append((sum(len(s) for _, s in records), live["now"]))
        return orig_write(records, out, *args, **kw)

    monkeypatch.setattr(store_mod, "write_bgzf_fasta", write)
    ingest(objs, tmp_path / "s", stage="S0", shard_bases=SMALL_SHARDS)
    largest = max(sum(len(r.seq) for r in _oracle_records(o.path)) for o in objs)
    stage = sum(len(s) for b in _expected_batches(objs, SMALL_SHARDS) for s in b.values())
    bound = SMALL_SHARDS + largest
    assert stage > bound  # the stage would not fit under the bound
    assert len(handed) == len(_expected_batches(objs, SMALL_SHARDS))
    assert all(n <= bound and alive <= bound for n, alive in handed), handed
    assert live["max"] <= bound, (live["max"], bound)
