"""Card W-36: a release whose stages are only partly labelled (AGENTS §3.9, card B-27, card W-18).

`rM.2` lists seven stages; S1, P1 and P2 were ingested before card B-27 and have no `contig_label` file, P3 to P6
have one. Card W-18's fallback to `contig_occurrence` applied only when a release had no label file at all, so the
contigs of the unlabelled stages had no accession: `bio-box find --accession` refused them, and the release check
of `rM.2` stopped at its first step.

The world: three roots built from the same inputs, each with stage S0 (the micro-S0 lambda and T4 genomes) and
stage S1 (a lambda genome with one base changed, under a name of its own, given by two source objects so that its
contig has two records), gene-called and taped, and release REL over both with a manifest and a catalog:

- `full`: both stages labelled, as ingest writes them since card B-27;
- `part`: S0 labelled, S1 built without a `contig_label` file (ingest run with the label writer left out, as the
  code before card B-27 ingested; never by deleting a published file);
- `none`: neither stage labelled (as `rN.3`).

S0's label file is written with one row per row group, so a lookup can show that it prunes. The oracle of a
contig's label is a loop over the stage's own `contig_occurrence` files (the first record by (l0_object_sha256,
record_ordinal) and the number of records), and of its bases the FASTA written here; neither goes through a
reader under test. Every reader must give the same bytes on the three roots.
"""

from __future__ import annotations

import json
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

from biobox import catalog, cli, fasta, ids, io
from biobox.l0seq import SourceObject, ingest
from biobox.l1.genecall import GENE_CALL_VERSION as V, GenomeInput, run as call_genes
from biobox.l3 import tapes
from biobox.search import export

pytestmark = [pytest.mark.r2, pytest.mark.covers("W-36")]

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
REL = "rW36.a"
STAGES = ["S0", "S1"]
KINDS = ("full", "part", "none")
LABELS = {"full": {"S0", "S1"}, "part": {"S0"}, "none": set()}
TV = tapes.tape_version([V])
S1_NAME, S1_AGAIN = "lambda_w36", "lambda_w36_again"


def _s1_seq() -> bytes:
    lam = fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq
    pos = 20_000
    return lam[:pos] + (b"C" if lam[pos:pos + 1].upper() != b"C" else b"G") + lam[pos + 1:]


def _ingest(root: Path, stage: str, objects, *, labels: bool, rows_per_group: int | None) -> None:
    """Ingest `objects` into `stage`; without `labels` the `contig_label` writer is left out (as before card B-27);
    with `rows_per_group` the label file is cut into row groups of that many rows."""
    mp = pytest.MonkeyPatch()
    real_parts, real_table = io.write_parts, io.write_table

    def write_parts(batches, spec, base_dir, *a, **kw):
        if spec.name == "contig_label" and not labels:
            return []
        return real_parts(batches, spec, base_dir, *a, **kw)

    def write_table(batches, spec, path, **kw):
        if spec.name == "contig_label" and rows_per_group is not None:
            kw["row_group_rows"] = rows_per_group
        return real_table(batches, spec, path, **kw)

    mp.setattr(io, "write_parts", write_parts)
    mp.setattr(io, "write_table", write_table)
    try:
        ingest(objects, root, stage)
    finally:
        mp.undo()


def _build(d: Path, kind: str) -> Path:
    root = d / kind
    lab = LABELS[kind]
    s0 = [o for o in CONTRACT["objects"] if o["kind"] == "virus"]  # lambda and T4
    _ingest(root, "S0", [SourceObject(FIX / o["file"], CONTRACT["source_name"], CONTRACT["source_release"],
                                      CONTRACT["license"], True) for o in s0], labels="S0" in lab, rows_per_group=1)
    call_genes(root, "S0", [GenomeInput((ids.contig_uid(fasta.read(FIX / o["file"])[0].seq),), o["genetic_code"],
                                        kind=o["kind"]) for o in s0])
    tapes.run(root, "S0", set())
    seq = _s1_seq()
    src = d / "src"
    src.mkdir(exist_ok=True)
    a, b = src / "s1.fna", src / "s1_again.fna"
    a.write_bytes(b"\n".join(fasta.format_record(f"{S1_NAME} planted lambda", seq, 70)) + b"\n")
    b.write_bytes(b"\n".join(fasta.format_record(f"{S1_AGAIN} the same bases", seq, 70)) + b"\n")
    _ingest(root, "S1", [SourceObject(a, "w36_one", "1", "test", True),
                         SourceObject(b, "w36_two", "1", "test", True, precedence=200)],
            labels="S1" in lab, rows_per_group=None)
    call_genes(root, "S1", [GenomeInput((ids.contig_uid(seq),), 11, kind="virus")])
    tapes.run(root, "S1", set())
    catalog.define_release(root, REL, STAGES)
    man = catalog.build_manifest(root, REL, STAGES)
    catalog.write_manifest(man)
    catalog.generate_catalog(man, root=root)
    return root


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    d = tmp_path_factory.mktemp("w36")
    return {k: _build(d, k) for k in KINDS}


def _occ_rows(root: Path, stage: str) -> list[dict]:
    out = []
    for f in sorted((root / f"l0seq/contig_occurrence/stage={stage}").rglob("*.parquet")):
        out += pq.read_table(f, partitioning=None).to_pylist()
    return out


def _label_oracle(root: Path) -> dict[bytes, tuple]:
    """{contig_uid: (accession, header, n_records, l0_object_sha256, record_ordinal)} over the release's stages:
    the first record by (l0_object_sha256, record_ordinal) and the count, by a loop over `contig_occurrence`."""
    best: dict[bytes, dict] = {}
    n: dict[bytes, int] = {}
    for s in STAGES:
        for r in _occ_rows(root, s):
            cu = r["contig_uid"]
            n[cu] = n.get(cu, 0) + 1
            if cu not in best or (r["l0_object_sha256"], r["record_ordinal"]) < (best[cu]["l0_object_sha256"],
                                                                               best[cu]["record_ordinal"]):
                best[cu] = r
    return {cu: (r["accession"], r["original_header"], n[cu], r["l0_object_sha256"], r["record_ordinal"])
            for cu, r in best.items()}


def _s1_contig() -> bytes:
    return ids.contig_uid(_s1_seq())


def _s1_accession(root: Path) -> str:
    return _label_oracle(root)[_s1_contig()][0]


def _occurrences(root: Path, stage: str) -> list[dict]:
    d = root / f"l1/protein_occurrence/gene_call_version={V}/stage={stage}"
    return [r for f in sorted(d.glob("*.parquet")) for r in pq.read_table(f, partitioning=None).to_pylist()]


def test_w36_the_world_is_built_with_and_without_label_files(world):
    """The construction: S1 of `part` and both stages of `none` have no `contig_label` file, and the manifests
    list exactly the label files of the labelled stages; the S1 contig has two records under two accessions."""
    for kind, root in world.items():
        for s in STAGES:
            assert bool(io.part_files(root / f"l0seq/contig_label/stage={s}")) == (s in LABELS[kind]), (kind, s)
        man = catalog.read_manifest(catalog.manifest_path(root, REL))
        assert {f["partitions"]["stage"] for f in man["files"] if f["table"] == "contig_label"} == LABELS[kind]
    rows = [r for r in _occ_rows(world["part"], "S1") if r["contig_uid"] == _s1_contig()]
    assert sorted(r["accession"] for r in rows) == sorted([S1_NAME, S1_AGAIN])
    md = pq.read_metadata(io.part_files(world["part"] / "l0seq/contig_label/stage=S0")[0])
    assert md.num_row_groups == 2  # lambda and T4, one row group each


@pytest.mark.parametrize("which", ["S1", "S0"])
def test_w36_find_by_accession_finds_the_genes_of_every_stage(world, which, capsys):
    """`bio-box find --accession` of a contig of the unlabelled stage (and, the control, of the labelled one)
    returns its genes, the accession and contig the oracle gives, the same bytes on the three roots."""
    outs = {}
    for kind, root in world.items():
        o = _label_oracle(root)
        cu = _s1_contig() if which == "S1" else ids.contig_uid(fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq)
        acc = o[cu][0]
        capsys.readouterr()
        code = cli.main(["find", "--root", str(root), "--release", REL, "--accession", acc])
        got = capsys.readouterr()
        assert code == 0, (kind, got.err)
        lines = got.out.splitlines()
        assert lines[0].split("\t")[:2] == ["accession", "contig_uid"]
        genes = {(r["start"], r["end"]) for r in _occurrences(root, which) if r["contig_uid"] == cu}
        rows = [x.split("\t") for x in lines[1:]]
        assert rows and {(x[0], x[1]) for x in rows} == {(acc, cu.hex())}, (kind, rows[:3])
        assert {(int(x[2]), int(x[3])) for x in rows} == genes
        outs[kind] = got.out
    assert outs["part"] == outs["full"] == outs["none"]


def test_w36_seq_by_accession_reads_a_contig_of_the_unlabelled_stage(world, capsys):
    want = _s1_seq()[100:160].decode()
    outs = {}
    for kind, root in world.items():
        acc = _s1_accession(root)
        capsys.readouterr()
        code = cli.main(["seq", "--root", str(root), "--release", REL, "--contig", f"{acc}:100-160"])
        got = capsys.readouterr()
        assert code == 0, (kind, got.err)
        head, body = got.out.splitlines()
        assert f"accession={acc} " in head and head.startswith(f">{_s1_contig().hex()} ")
        assert body == want
        outs[kind] = got.out
    assert outs["part"] == outs["full"] == outs["none"]


def _hits(root: Path) -> list[dict]:
    """One hit for each protein of the S1 contig and for the first protein of S0's lambda."""
    lam = ids.contig_uid(fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq)
    s1 = sorted({r["protein_uid"] for r in _occurrences(root, "S1") if r["contig_uid"] == _s1_contig()})
    s0 = sorted({r["protein_uid"] for r in _occurrences(root, "S0") if r["contig_uid"] == lam})[:1]
    return [{"target_uid": u, "target_kind": "protein", "seq_bitscore": 50.0, "seq_evalue_zref": 1e-10}
            for u in s1 + s0]


def test_w36_hit_exports_name_the_accession_of_every_stage(world):
    """`search hits --format tsv` and `--format faa` (`export.hits_tsv`, `hits_faa`): every gene of a hit protein
    has its contig's accession, header and record count from the oracle, the same bytes on the three roots."""
    tsv, faa = {}, {}
    for kind, root in world.items():
        hits = _hits(root)
        o = _label_oracle(root)
        lines = list(export.hits_tsv(root, hits, release=REL, tape_version=TV, gene_call_version=V))
        cols = lines[0].split("\t")
        rows = [dict(zip(cols, x.split("\t"))) for x in lines[1:]]
        assert rows
        for r in rows:
            acc, header, n = o[bytes.fromhex(r["contig_uid"])][:3]
            assert (r["accession"], r["n_records"]) == (acc, str(n)), (kind, r["contig_uid"], r["accession"])
            assert r["header"] == header
        assert _s1_accession(root) in {r["accession"] for r in rows}, kind
        tsv[kind] = "\n".join(lines)
        faa[kind] = "\n".join(export.hits_faa(root, hits, release=REL, tape_version=TV, gene_call_version=V))
        assert f"example={_s1_accession(root)}:" in faa[kind], kind
    assert tsv["part"] == tsv["full"] == tsv["none"]
    assert faa["part"] == faa["full"] == faa["none"]


LABEL_COLS = "contig_uid, accession, original_header, n_records, l0_object_sha256, record_ordinal, n_occurrences"


def _view_rows(cat: Path) -> list[tuple]:
    con = catalog.open_catalog(cat, threads=1)
    try:
        return con.execute(f"SELECT {LABEL_COLS} FROM contig_label ORDER BY contig_uid").fetchall()
    finally:
        con.close()


def test_w36_the_catalog_label_view_holds_every_contig_of_the_release(world):
    """The catalog's `contig_label` view: one row per contig of both stages, as the oracle, on the three roots."""
    got = {kind: _view_rows(catalog.catalog_path(root, REL)) for kind, root in world.items()}
    for kind, root in world.items():
        want = [(cu, *v[:2], v[2], v[3], v[4], v[2]) for cu, v in sorted(_label_oracle(root).items())]
        assert got[kind] == want, kind
    assert got["part"] == got["full"] == got["none"]


def test_w36_a_catalog_built_before_the_card_is_read_with_every_stage_labelled(world, tmp_path, monkeypatch):
    """A catalog of the partly labelled release made by the code before this card (its stored `contig_label` view
    reads the label files only): `open_catalog` gives the view of every stage all the same."""
    root = world["part"]
    real = catalog.label_view_sql
    monkeypatch.setattr(catalog, "label_view_sql", lambda files, n_stages, *a, **kw: real(files, n_stages))
    man = catalog.read_manifest(catalog.manifest_path(root, REL))
    old = catalog.generate_catalog(man, tmp_path / "old.duckdb")
    monkeypatch.undo()
    raw = duckdb.connect(str(old), read_only=True)
    try:  # the stored view is the old one: S1's contig is missing from it
        stored = {r[0] for r in raw.execute("SELECT contig_uid FROM contig_label").fetchall()}
    finally:
        raw.close()
    assert _s1_contig() not in stored
    assert _view_rows(old) == _view_rows(catalog.catalog_path(world["full"], REL))


def _reads(root: Path, scans: list[dict]) -> dict[tuple[str, str], int]:
    """Row groups read per (table, stage) of the profiled scans."""
    out: dict[tuple[str, str], int] = {}
    for s in scans:
        for f in s["files"]:
            spec, parts = catalog.classify(Path(f).resolve().relative_to(Path(root).resolve()).as_posix())
            key = (spec.name, parts.get("stage"))
            out[key] = out.get(key, 0) + s["row_groups_scanned"]
            break  # a scan's row groups are counted once, under its first file's table and stage
    return out


def test_w36_a_label_lookup_by_contig_still_prunes_on_the_labelled_stage(world, tmp_path):
    """AGENTS §4.4: in the partly labelled release a lookup of one S0 contig by contig_uid, from a joined table (the
    exports) and as a literal (the catalog's view), reads one of the two row groups of S0's `contig_label` and no
    `contig_occurrence` row group of S0; the unlabelled S1 is read from its `contig_occurrence` (card W-18's cost:
    every row group of that stage's file, one here)."""
    root = world["part"]
    lam = ids.contig_uid(fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq)
    with export.open_store(root, release=REL, gene_call_version=V, threads=1, relations=("lab",)) as st:
        st.con.register("_cin", pa.table({"contig_uid": pa.array([lam], pa.binary(16))}))
        st.con.execute("CREATE TEMP TABLE _c AS SELECT contig_uid FROM _cin")
        q = catalog.profile_query(st.con, "SELECT contig_uid, accession FROM lab WHERE contig_uid IN "
                                          "(SELECT contig_uid FROM _c)", profile_dir=tmp_path)
    assert [r[1] for r in q["result"]] == [_label_oracle(root)[lam][0]]
    reads = _reads(root, q["scans"])
    assert reads.get(("contig_label", "S0")) == 1, reads
    assert ("contig_occurrence", "S0") not in reads, reads
    con = catalog.open_catalog(catalog.catalog_path(root, REL), threads=1)
    try:
        q = catalog.profile_query(con, "SELECT accession FROM contig_label WHERE contig_uid = ?", [lam],
                                  profile_dir=tmp_path)
    finally:
        con.close()
    assert q["result"] == [(_label_oracle(root)[lam][0],)]
    reads = _reads(root, q["scans"])
    assert reads.get(("contig_label", "S0")) == 1, reads
    assert ("contig_occurrence", "S0") not in reads, reads


def test_w36_release_verify_names_the_stages_without_label_files(world, capsys):
    """`release verify` of the partly labelled release passes (exit 0) and warns, naming S1; of the fully labelled
    one it says nothing of labels; of the one without labels it keeps card W-17's warning."""
    errs = {}
    for kind, root in world.items():
        capsys.readouterr()
        assert cli.main(["release", "verify", "--root", str(root), "--release", REL]) == 0, kind
        got = capsys.readouterr()
        assert f"release {REL}: verified" in got.out
        errs[kind] = [x for x in got.err.splitlines() if x.startswith("warning") and "contig_label" in x]
    assert errs["full"] == []
    assert len(errs["part"]) == 1 and "stages S1 " in errs["part"][0] and "S0" not in errs["part"][0], errs["part"]
    assert len(errs["none"]) == 1 and "built before card B-27" in errs["none"][0], errs["none"]


def test_w36_releases_shows_a_partly_labelled_release_as_partial(world, capsys, monkeypatch):
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)
    shown = {}
    for kind, root in world.items():
        capsys.readouterr()
        assert cli.main(["releases", "--root", str(root), "--json"]) == 0
        rec = json.loads(capsys.readouterr().out.splitlines()[0])
        assert rec["lookups"]["contig_label"] is (kind != "none")
        assert rec.get("unlabelled_stages") == (["S1"] if kind == "part" else None), (kind, rec)
        assert cli.main(["releases", "--root", str(root)]) == 0
        lines = [x.split("\t") for x in capsys.readouterr().out.splitlines()]
        shown[kind] = dict(zip(lines[0], lines[1]))["contig_label"]
    assert shown == {"full": "yes", "part": "partial(not S1)", "none": "no"}
