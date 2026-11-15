"""Card C-34: `bio-box stage labels --stage S`, the label files of a stage ingested before card B-27.

Stages S1, P1 and P2 of the data root have no `contig_label` file: they were ingested before card B-27, which made
ingest write that table. The step writes `l0seq/contig_label/stage=<S>/` as new files from the stage's own
`contig_occurrence`, with the code ingest uses (`biobox.l0seq.store._label_batches`, `io.write_parts`), and a ledger
`ledger/stage=<S>/step=labels`; under `run-step` it has a run record. It refuses a stage that has label files.

The world is card W-36's (`tests/r2/test_w36_partial_labels.py`, `_build`): roots built from the same inputs with
stage S0 (micro-S0 lambda and T4) and S1 (a lambda genome with one base changed, two source objects), gene-called and
taped, and release `rW36.a` over both with its manifest and catalog: `part` (S1 ingested with the label writer left
out, as before card B-27), `none` (neither stage labelled, as `rN.3`) and `full` (both labelled by ingest). Nothing is
deleted from a published release to make a stage unlabelled.

Oracles, none of them the step under test:
- a re-ingest of the same inputs with the label writer in place (ingest of card B-27): its label files must have the
  same sha256 as the step's, on W-36's stages and on a constructed stage with a contig another stage owns (DEC-26),
  one source object listed twice, a contig of three records and parts cut small;
- a loop over the stage's `contig_occurrence` files for each label row and for the ledger's counts.

What a published release says after its stages gain label files (card C-34's rule, `catalog.STAGE_ADDED_TABLES`; the
rule is stated and tested in `tests/r1/test_c34_verify_rule.py`): `release verify` of `rW36.a` on `part` and `none`
verifies as before, and every reader of it (`find`, `seq`, `search hits` TSV and FASTA, the catalog's `contig_label`
view, `release verify`, `releases`) gives the same bytes before and after: readers take the manifest's file lists
(AGENTS §3.9), which do not name the new files. A release defined afterwards over the same stages lists them and is
fully labelled (`releases`: `contig_label yes`).
"""

from __future__ import annotations

import contextlib
import hashlib
import io as stdio
import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from biobox import catalog, cli, fasta, ids, io, ledger, runstep
from biobox.l0seq import SourceObject, ingest

sys.path.insert(0, str(Path(__file__).parent))
import test_builders_out as c08  # noqa: E402
import test_w36_partial_labels as w36  # noqa: E402

pytestmark = [pytest.mark.r2, pytest.mark.covers("C-34")]

REL, STAGES = w36.REL, w36.STAGES
LABEL = "l0seq/contig_label"
NEW_REL = "rC34.n"
NOON = datetime(2026, 11, 15, 12, 0)  # a Sunday: the night cap applies all day (only the threads depend on it)


def _cli(*argv: str) -> tuple[int, str, str]:
    out, err = stdio.StringIO(), stdio.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            code = cli.main(list(argv))
        except SystemExit as e:  # argparse: an unknown step
            code = e.code if isinstance(e.code, int) else 2
    return code, out.getvalue(), err.getvalue()


@contextlib.contextmanager
def _no_env(*names: str):
    old = {n: os.environ.pop(n) for n in names if n in os.environ}
    try:
        yield
    finally:
        os.environ.update(old)


def _sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _label_files(root: Path, stage: str) -> list[Path]:
    return io.part_files(Path(root) / f"{LABEL}/stage={stage}")


def _listing(root: Path) -> dict[str, tuple[int, int]]:
    return {p.relative_to(root).as_posix(): (p.stat().st_size, p.stat().st_mtime_ns)
            for p in sorted(Path(root).rglob("*")) if p.is_file()}


def _reads(root: Path) -> dict:
    """What every reader of the published release `rW36.a` gives (stdout and stderr, or rows)."""
    lam = ids.contig_uid(fasta.read(w36.FIX / "NC_001416.1.fna.gz")[0].seq)
    o = w36._label_oracle(root)
    got = {}
    for name, cu in (("S1", w36._s1_contig()), ("S0", lam)):
        got[f"find {name}"] = _cli("find", "--root", str(root), "--release", REL, "--accession", o[cu][0])
    got["seq"] = _cli("seq", "--root", str(root), "--release", REL, "--contig", f"{w36._s1_accession(root)}:100-160")
    hits = w36._hits(root)
    err = stdio.StringIO()
    with contextlib.redirect_stderr(err):
        got["hits tsv"] = list(w36.export.hits_tsv(root, hits, release=REL, tape_version=w36.TV,
                                                    gene_call_version=w36.V))
        got["hits faa"] = list(w36.export.hits_faa(root, hits, release=REL, tape_version=w36.TV,
                                                    gene_call_version=w36.V))
    got["hits stderr"] = err.getvalue()
    got["view"] = w36._view_rows(catalog.catalog_path(root, REL))
    got["verify"] = _cli("release", "verify", "--root", str(root), "--release", REL)
    got["releases"] = _cli("releases", "--root", str(root))
    got["releases json"] = _cli("releases", "--root", str(root), "--json")
    return got


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    d = tmp_path_factory.mktemp("c34")
    return {k: w36._build(d, k) for k in ("part", "none", "full")}


@pytest.fixture(scope="module")
def labelled(world):
    """On `part` and `none`: every reader before, `stage labels` on each unlabelled stage, every reader after."""
    out = {}
    with _no_env("BIOBOX_STAGING", "BIOBOX_PROFILE"):
        for kind in ("part", "none"):
            root = world[kind]
            stages = [s for s in STAGES if not _label_files(root, s)]
            before, listing = _reads(root), _listing(root)
            runs = {s: _cli("stage", "labels", "--root", str(root), "--stage", s) for s in stages}
            out[kind] = {"stages": stages, "before": before, "listing": listing, "runs": runs,
                         "after": _reads(root), "after_listing": _listing(root)}
    return out


def _stage_oracle(root: Path, stage: str) -> list[tuple]:
    """The label rows of one stage by a loop over its `contig_occurrence` files: per contig, in contig_uid order,
    the first record by (l0_object_sha256, record_ordinal) and the stage's number of records of the contig."""
    best: dict[bytes, dict] = {}
    n: dict[bytes, int] = {}
    for f in sorted((Path(root) / f"l0seq/contig_occurrence/stage={stage}").rglob("*.parquet")):
        for r in pq.read_table(f, partitioning=None).to_pylist():
            cu = r["contig_uid"]
            n[cu] = n.get(cu, 0) + 1
            if cu not in best or (r["l0_object_sha256"], r["record_ordinal"]) < (best[cu]["l0_object_sha256"],
                                                                               best[cu]["record_ordinal"]):
                best[cu] = r
    return [(cu, best[cu]["accession"], best[cu]["original_header"], n[cu], best[cu]["l0_object_sha256"],
             best[cu]["record_ordinal"]) for cu in sorted(best)]


def _label_rows(root: Path, stage: str) -> list[tuple]:
    cols = ("contig_uid", "accession", "original_header", "n_records", "l0_object_sha256", "record_ordinal")
    rows = [r for f in _label_files(root, stage) for r in pq.read_table(f, partitioning=None).to_pylist()]
    return [tuple(r[c] for c in cols) for r in rows]


def test_c34_the_step_writes_the_label_files_of_each_unlabelled_stage(labelled, world):
    """Each unlabelled stage gets new label files, one row per contig of its `contig_occurrence` (the oracle loop),
    and its ledger; no file that existed before changed."""
    assert labelled["part"]["stages"] == ["S1"] and labelled["none"]["stages"] == ["S0", "S1"]
    for kind, got in labelled.items():
        for s, (code, out, err) in got["runs"].items():
            assert code == 0, (kind, s, err)
            res = json.loads(out.splitlines()[-1])
            assert (res["stage"], res["step"]) == (s, "labels")
            root = world[kind]
            assert _label_rows(root, s) == _stage_oracle(root, s), (kind, s)
            assert res["result"]["contigs"] == len(_stage_oracle(root, s))
            assert ledger.path(root, s, "labels").is_file()
        before, after = got["listing"], got["after_listing"]
        assert {k: v for k, v in after.items() if k in before} == before, kind
        new = sorted(set(after) - set(before))
        want = sorted([f"{LABEL}/stage={s}/part-0.parquet" for s in got["stages"]]
                      + [f"ledger/stage={s}/step=labels/part-0.parquet" for s in got["stages"]])
        assert new == want, (kind, new)


def _oracle_root(d: Path, stage: str, objects, **kw) -> Path:
    """A re-ingest of `objects` as `stage` into a new root, with the label writer of card B-27."""
    root = d / f"oracle-{stage}"
    ingest(objects, root, stage, **kw)
    return root


def _w36_objects(stage: str, src: Path) -> list[SourceObject]:
    if stage == "S0":
        return [SourceObject(w36.FIX / o["file"], w36.CONTRACT["source_name"], w36.CONTRACT["source_release"],
                             w36.CONTRACT["license"], True) for o in w36.CONTRACT["objects"] if o["kind"] == "virus"]
    return [SourceObject(src / "s1.fna", "w36_one", "1", "test", True),
            SourceObject(src / "s1_again.fna", "w36_two", "1", "test", True, precedence=200)]


def test_c34_labels_equal_those_of_a_reingest_of_the_same_inputs(labelled, world, tmp_path):
    """The oracle: for every stage the step labelled, a re-ingest of the stage's inputs with ingest's label writer
    gives label files of the same sha256 (and, for S1 of `part`, the label file ingest wrote in `full`)."""
    src = world["part"].parent / "src"
    for kind in ("part", "none"):
        for s in labelled[kind]["stages"]:
            want = _label_files(_oracle_root(tmp_path / kind, s, _w36_objects(s, src)), s)
            got = _label_files(world[kind], s)
            assert [p.name for p in got] == [p.name for p in want] and want, (kind, s)
            assert [_sha(p) for p in got] == [_sha(p) for p in want], (kind, s)
    assert [_sha(p) for p in _label_files(world["part"], "S1")] == [_sha(p) for p in _label_files(world["full"], "S1")]


def test_c34_a_constructed_stage_with_shared_contigs_and_small_parts(tmp_path):
    """DEC-26 and ties: stage S2 holds T4 (a contig S0 stored first and owns), lambda's planted copy in two source
    objects of different names (one contig, three records: one object is the same bytes listed under two sources), and
    its parts are cut small. The step's label files equal a re-ingest's, part for part, by sha256; the T4 contig has its
    label row in S2 as ingest gives it; the ledger balances (L-01) and counts the records and contigs."""
    src = tmp_path / "src"
    src.mkdir()
    seq = w36._s1_seq()
    a, b = src / "s1.fna", src / "s1_again.fna"
    a.write_bytes(b"\n".join(fasta.format_record("c34_planted planted lambda", seq, 70)) + b"\n")
    b.write_bytes(b"\n".join(fasta.format_record("c34_planted_again the same bases", seq, 60)) + b"\n")
    t4 = next(o for o in w36.CONTRACT["objects"] if o["file"].startswith("NC_000866"))
    objs = [SourceObject(w36.FIX / t4["file"], "c34_t4", "1", "test", True),
            SourceObject(a, "c34_one", "1", "test", True), SourceObject(a, "c34_one_too", "1", "test", True),
            SourceObject(b, "c34_two", "1", "test", True, precedence=200)]
    s0 = _w36_objects("S0", src)
    mfb = 50  # one label row per part (48 estimated bytes a row, store.EST_LABEL_ROW)
    roots = {}
    for kind, labels in (("step", False), ("oracle", True)):
        r = tmp_path / kind
        w36._ingest(r, "S0", s0, labels=True, rows_per_group=None)
        mp = pytest.MonkeyPatch()
        real = io.write_parts
        mp.setattr(io, "write_parts", lambda batches, spec, *x, **kw: [] if spec.name == "contig_label" and not labels
                   else real(batches, spec, *x, **kw))
        try:
            ingest(objs, r, "S2", max_file_bytes=mfb)
        finally:
            mp.undo()
        roots[kind] = r
    assert not _label_files(roots["step"], "S2")
    with _no_env("BIOBOX_STAGING"):
        code, out, err = _cli("stage", "labels", "--root", str(roots["step"]), "--stage", "S2",
                              "--max-file-bytes", str(mfb))
    assert code == 0, err
    got, want = _label_files(roots["step"], "S2"), _label_files(roots["oracle"], "S2")
    assert len(want) >= 2 and [p.name for p in got] == [p.name for p in want]
    assert [_sha(p) for p in got] == [_sha(p) for p in want]
    rows = _label_rows(roots["step"], "S2")
    assert rows == _stage_oracle(roots["step"], "S2")
    t4_uid = ids.contig_uid(fasta.read(w36.FIX / t4["file"])[0].seq)
    assert t4_uid in {r[0] for r in rows}  # a contig S0 owns: labelled in S2 too, as ingest labels it
    assert {r[0]: r[3] for r in rows}[ids.contig_uid(seq)] == 3
    lrows = ledger.read(roots["step"], "S2", "labels")
    assert ledger.check_balance(roots["step"], "S2", ["labels"]) == []
    assert sum(r["n_in"] for r in lrows) == sum(r[3] for r in rows) == 4
    assert sum(r["n_kept"] for r in lrows) == len(rows) == 2
    assert {r["object_kind"] for r in lrows} == {"l0_object"}
    tools = json.loads(lrows[0]["tool_versions_json"])
    assert {"pyarrow", "duckdb"} <= set(tools)


CHILD = r"""
import sys
from pathlib import Path
from biobox.l0seq import store
store.write_labels(Path(sys.argv[1]), sys.argv[2], out=Path(sys.argv[3]), threads=int(sys.argv[4]))
"""


def test_c34_the_files_are_the_same_at_any_thread_count_input_order_and_process(tmp_path, monkeypatch):
    """AGENTS §4.3: the label and ledger files written at 1 and at 4 DuckDB threads, with the stage's
    `contig_occurrence` files listed in reverse, and in a child process, have the same sha256."""
    from biobox.l0seq import store

    src = tmp_path / "src"
    src.mkdir()
    seq = w36._s1_seq()
    for name, head in (("s1.fna", "w36"), ("s1_again.fna", "w36 again")):
        (src / name).write_bytes(b"\n".join(fasta.format_record(head, seq, 70)) + b"\n")
    root = tmp_path / "root"
    w36._ingest(root, "S1", [*_w36_objects("S0", src), *_w36_objects("S1", src)], labels=False, rows_per_group=None)
    assert len(store.occurrence_files(root, "S1")) >= 2  # several sources: an order to reverse
    outs = {}
    for name, threads in (("t1", 1), ("t4", 4)):
        outs[name] = tmp_path / name
        store.write_labels(root, "S1", out=outs[name], threads=threads)
    real = store.occurrence_files
    monkeypatch.setattr(store, "occurrence_files", lambda r, s: list(reversed(real(r, s))))
    outs["reversed"] = tmp_path / "reversed"
    store.write_labels(root, "S1", out=outs["reversed"])
    monkeypatch.undo()
    outs["child"] = tmp_path / "child"
    p = subprocess.run([sys.executable, "-c", CHILD, str(root), "S1", str(outs["child"]), "2"], capture_output=True,
                       text=True, timeout=600)
    assert p.returncode == 0, p.stderr[-3000:]
    shas = {k: {f.relative_to(d).as_posix(): _sha(f) for f in sorted(d.rglob("*.parquet"))} for k, d in outs.items()}
    assert len(shas["t1"]) == 2, shas["t1"]  # the label file and the ledger
    assert shas["t1"] == shas["t4"] == shas["reversed"] == shas["child"], shas


def test_c34_refuses_a_stage_that_has_label_files(labelled, world):
    """A stage labelled by ingest (`full`), or by the step already (`part`'s S1), is refused (exit 2) and nothing is
    written; so is a stage with no `contig_occurrence`."""
    with _no_env("BIOBOX_STAGING"):
        for root, s in ((world["full"], "S0"), (world["full"], "S1"), (world["part"], "S1"), (world["part"], "S0")):
            before = _listing(root)
            code, out, err = _cli("stage", "labels", "--root", str(root), "--stage", s)
            assert code == 2 and err.startswith("refused: ") and "contig_label" in err, (root.name, s, err)
            assert _listing(root) == before
        code, _, err = _cli("stage", "labels", "--root", str(world["part"]), "--stage", "S9")
        assert code == 2 and "S9" in err, err


def test_c34_a_published_release_reads_as_before(labelled):
    """No reader of `rW36.a` changes its output when its stages gain label files: `find`, `seq`, the hit TSV and
    FASTA, the catalog's `contig_label` view, `release verify` (exit 0, the same lines and warnings) and `releases`."""
    for kind, got in labelled.items():
        before, after = got["before"], got["after"]
        assert before["verify"][0] == 0, before["verify"]
        for key in before:
            assert after[key] == before[key], (kind, key)
    part = labelled["part"]["after"]
    assert "stages S1 have no contig_label files" in part["verify"][2]  # W-36's warning: the manifest is as it was
    assert json.loads(part["releases json"][1].splitlines()[0])["unlabelled_stages"] == ["S1"]


def test_c34_a_release_defined_afterwards_is_fully_labelled(labelled, world):
    """A release defined over the same stages after the step lists every stage's label files: it verifies with no
    label warning, `releases` shows `contig_label yes` for it (the published one still `partial(not S1)` and `no`),
    and its catalog's `contig_label` view has the same rows as the published release's."""
    for kind, old in (("part", "partial(not S1)"), ("none", "no")):
        root = world[kind]
        catalog.define_release(root, NEW_REL, STAGES)
        man = catalog.build_manifest(root, NEW_REL, STAGES)
        catalog.write_manifest(man)
        catalog.generate_catalog(man, root=root)
        assert {f["partitions"]["stage"] for f in man["files"] if f["table"] == "contig_label"} == set(STAGES)
        assert catalog.unlabelled_stages(man) == []
        code, out, err = _cli("release", "verify", "--root", str(root), "--release", NEW_REL)
        assert code == 0 and f"release {NEW_REL}: verified" in out, (kind, out, err)
        assert "contig_label" not in err, (kind, err)
        with _no_env("BIOBOX_PROFILE"):
            code, out, _ = _cli("releases", "--root", str(root))
        assert code == 0
        lines = [x.split("\t") for x in out.splitlines()]
        shown = {r[0]: dict(zip(lines[0], r))["contig_label"] for r in lines[1:]}
        assert shown == {REL: old, NEW_REL: "yes"}, (kind, shown)
        assert (w36._view_rows(catalog.catalog_path(root, NEW_REL))
                == w36._view_rows(catalog.catalog_path(root, REL))), kind
        assert catalog.verify_manifest(catalog.read_manifest(catalog.manifest_path(root, REL)), root) == []


def test_c34_run_step_publishes_the_files_with_a_run_record_and_ledger(tmp_path):
    """Under `run-step` the step publishes the label file and its ledger and writes its run record; a second call is
    a resume (exit 0, nothing run); the step called again directly is refused."""
    root = tmp_path / "root"
    w36._ingest(root, "S0", _w36_objects("S0", tmp_path), labels=False, rows_per_group=None)
    prof = c08.profile_at(tmp_path / "prof")
    cmd = [sys.executable, "-m", "biobox.cli", "stage", "labels", "--root", str(root), "--stage", "S0"]
    with _no_env("BIOBOX_STAGING"):
        assert runstep.run_step(prof, root, "S0", "labels", cmd, clock=lambda: NOON, load=lambda: 0.0) == 0
        rec_path = root / "runs/stage=S0/labels.json"
        rec = json.loads(rec_path.read_text())
        outputs = sorted(o["path"] for o in rec["outputs"])
        assert outputs == [f"{LABEL}/stage=S0/part-0.parquet", "ledger/stage=S0/step=labels/part-0.parquet"], outputs
        assert rec["step"] == "labels" and rec["stage"] == "S0"
        before = _listing(root)
        assert runstep.run_step(prof, root, "S0", "labels", cmd, clock=lambda: NOON, load=lambda: 0.0) == 0
        assert _listing(root) == before
        code, _, err = _cli(*cmd[3:])
        assert code == 2 and "refused" in err
    assert _label_rows(root, "S0") == _stage_oracle(root, "S0")
    assert ledger.check_balance(root, "S0", ["labels"]) == []
