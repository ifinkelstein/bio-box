"""Ring R2: a release names its stages (card B-25; import ladder 2026-11-10, finding 1).

A release is defined once, by a `release_def` file listing its stages in order. Every release-level
builder (clustering, tokens, postings, `contig_dim`) and the manifest read only those stages: a stage
directory on disk that the release does not list is ignored, and a listed stage that is missing is
an error.

The acceptance build: release R1 over stage S0 (the two phages of micro-S0), then a second stage S1
(M. genitalium) is added to the store: ingested, gene-called, scanned with the S0 library and with a
second library that only S1 has, and taped, so a builder that reads every `stage=*` or `library=*`
directory on disk would take S1 in. R1 is then rebuilt into a fresh copy of the store (its stage
files and its `release_def`, not its release files) from the stages that `release_def` names. The
oracle is the first build, made before S1 existed: every R1 file and the manifest must be
byte-identical to it.

Covers D-01 (the same inputs give the same bytes) and LY-10 (the manifest is an explicit file list).
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pyarrow.parquet as pq
import pyhmmer
import pytest
import yaml

from biobox import catalog, cli, dim, fasta, ids
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l1.genecall import run as call_genes
from biobox.l3 import postings, tapes, tokens
from biobox.registry.schemas import conformance_problems, table

pytestmark = pytest.mark.r2

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
MMSEQS = shutil.which("mmseqs")
REL = "rS0.b25"
# S0: the two phages; S1: the prokaryote. S1 alone also has library a-1, which sorts before test-1, so
# a builder that took every library on disk would renumber the family tokens of S0.
STAGE_OBJECTS = {"S0": ("NC_001416.1", "NC_000866.4"), "S1": ("NC_000908.2",)}
LIBRARIES = {"S0": ("test",), "S1": ("test", "a")}


def _objects(stage: str) -> list[dict]:
    return [o for o in CONTRACT["objects"] if o["accession"] in STAGE_OBJECTS[stage]]


def _library(root: Path, path: Path) -> None:
    """Six HMMs built from the longest lambda proteins of S0 (as test_micro_s0_l3 does)."""
    lam = ids.contig_uid(fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq)
    prots = domains.load_proteins(root, V, "S0")
    occ = pq.read_table(root / f"l1/protein_occurrence/gene_call_version={V}/stage=S0/part-0.parquet").to_pylist()
    uids = sorted({r["protein_uid"] for r in occ if r["contig_uid"] == lam}, key=lambda u: (-len(prots[u]), u))[:6]
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder, bg = pyhmmer.plan7.Builder(alphabet, seed=42), pyhmmer.plan7.Background(alphabet)
    with open(path, "wb") as fh:
        for i, u in enumerate(uids):
            seq = pyhmmer.easel.TextSequence(name=b"s", sequence=prots[u].replace("*", "")).digitize(alphabet)
            hmm, _, _ = builder.build(seq, bg)
            hmm.name, hmm.accession = f"lam{i}".encode(), f"TEST{i:04d}".encode()
            hmm.cutoffs.gathering = (25.0, 25.0)
            hmm.write(fh)


def add_stage(root: Path, stage: str, lib: Path) -> None:
    """Every per-stage step of one stage: ingest, genome inputs, gene calls, domain scans, tapes."""
    objs = _objects(stage)
    ingest([SourceObject(FIX / o["file"], CONTRACT["source_name"], CONTRACT["source_release"], CONTRACT["license"], True)
            for o in objs], root, stage)
    genomes = []
    for o in objs:
        cu = ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)
        genomes.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"],
                                   circular=(cu,) if o["topology"] == "circular" else ()))
    write_genome_input(root, stage, genomes)
    call_genes(root, stage)
    if not lib.exists():
        _library(root, lib)
    for name in LIBRARIES[stage]:
        domains.run(root, stage, V, name, "1", lib, "gathering")
    tapes.run(root, stage)


def build_release(root: Path, release: str, stages: list[str], work: Path) -> dict:
    """Every release-level builder, each given the release's stages."""
    out = {}
    if MMSEQS:
        from biobox.l2 import cluster

        out["cluster"] = cluster.run(root, release, V, stages, tmp=work)
    out["tokens"] = tokens.run(root, release, stages)
    out["postings"] = postings.run(root, release, stages)
    out["dim"] = dim.run(root, release, stages, tmp=work)
    return out


def _stage_files_only(src: Path):
    """copytree filter: leave out the files R1's builders wrote (its `release=` directories other than
    `release_def`, and the clustering ledger), so the copy holds the stages and the release's definition."""
    defs = src / table("release_def").location

    def ignore(d, names):
        return {n for n in names if (n == f"release={REL}" and Path(d) != defs) or n == "step=cluster"}
    return ignore


def _portable(man: dict) -> dict:
    return {k: v for k, v in man.items() if k != "root"}


@pytest.fixture(scope="module")
def ladder(tmp_path_factory):
    d = tmp_path_factory.mktemp("b25")
    a, b, work, lib = d / "a", d / "b", d / "work", d / "lib.hmm"
    work.mkdir()
    add_stage(a, "S0", lib)
    assert catalog.define_release(a, REL, ["S0"]) == ["S0"]
    first = build_release(a, REL, ["S0"], work)
    m1 = catalog.build_manifest(a, REL, ["S0"])
    # the oracle: R1's files as the first build left them, before S1 existed
    snap = {f["path"]: (a / f["path"]).read_bytes() for f in m1["files"] + m1["artifacts"]}
    add_stage(a, "S1", lib)
    shutil.copytree(a, b, ignore=_stage_files_only(a))
    stages = catalog.release_stages(b, REL)
    second = build_release(b, REL, stages, work)
    m2 = catalog.build_manifest(b, REL, stages)
    return {"a": a, "b": b, "work": work, "m1": m1, "m2": m2, "snap": snap, "first": first, "second": second,
            "stages": stages}


@pytest.mark.covers("D-01")
def test_b25_rebuild_after_a_later_stage_gives_byte_identical_release_files(ladder):
    a, b, snap = ladder["a"], ladder["b"], ladder["snap"]
    assert ladder["stages"] == ["S0"]
    assert ladder["second"] == ladder["first"]
    # S1 is really in both stores, at every place a glob would have found it
    for rel in ("l0seq/contig_index/stage=S1", f"l1/feature/producer={V}/stage=S1", "l1/domains/library=a-1/stage=S1",
                f"l0seq/protein_seq/gene_call_version={V}/stage=S1", f"l3/tapes/tape_version={tapes.tape_version([V])}/stage=S1"):
        assert (a / rel).is_dir() and (b / rel).is_dir(), rel
    release_files = sorted(p for p in snap if f"release={REL}/" in p or "/step=cluster/" in p)
    want = {"release_def", "token_dict", "tape_tokens", "token_postings", "contig_dim"} | (
        {"cluster_member", "cluster_member_aln", "ledger"} if MMSEQS else set())
    assert {f["table"] for f in ladder["m1"]["files"] if f["path"] in release_files} == want
    for p in release_files:
        assert (b / p).read_bytes() == snap[p], p
    # nothing else of R1 was written in the rebuild
    rebuilt = sorted(p.relative_to(b).as_posix() for p in b.rglob("*")
                     if p.is_file() and (f"release={REL}/" in p.relative_to(b).as_posix()))
    assert rebuilt == sorted(p for p in release_files if f"release={REL}/" in p)


@pytest.mark.covers("D-01", "LY-10")
def test_b25_manifest_of_the_rebuild_equals_the_first(ladder):
    m1, m2 = ladder["m1"], ladder["m2"]
    assert m1["stages"] == m2["stages"] == ["S0"]
    assert not [f["path"] for f in m2["files"] + m2["artifacts"] if f["partitions"].get("stage") == "S1"]
    assert not [f["path"] for f in m2["files"] if f["partitions"].get("library") == "a-1"]

    def entries(m):
        return {f["path"]: f for f in m["files"] + m["artifacts"]}

    e1, e2 = entries(m1), entries(m2)
    assert sorted(e1) == sorted(e2)
    differ = sorted(p for p in e1 if e1[p] != e2[p])
    assert differ == []
    same = {k: v for k, v in _portable(m1).items() if k not in ("files", "artifacts")}
    assert same == {k: v for k, v in _portable(m2).items() if k not in ("files", "artifacts")}
    assert catalog.manifest_bytes(_portable(m1)) == catalog.manifest_bytes(_portable(m2))


@pytest.mark.covers("LY-10")
def test_b25_verify_ignores_a_stage_the_release_does_not_list(ladder, tmp_path):
    a, m1 = ladder["a"], ladder["m1"]
    assert catalog.verify_manifest(m1, a) == []
    # a stray file in a listed stage is still reported
    c = tmp_path / "c"
    shutil.copytree(a, c)
    stray = c / "l0seq/contig_index/stage=S0/part-9.parquet"
    shutil.copy(c / "l0seq/contig_index/stage=S0/part-0.parquet", stray)
    assert "not in manifest: l0seq/contig_index/stage=S0/part-9.parquet" in catalog.verify_manifest(m1, c)


def test_b25_a_listed_stage_that_is_missing_raises(ladder, tmp_path):
    b, work = ladder["b"], ladder["work"]
    c = tmp_path / "c"
    shutil.copytree(b, c)
    rel, stages = "rS0.missing", ["S0", "S9"]
    if MMSEQS:
        from biobox.l2 import cluster

        with pytest.raises(FileNotFoundError, match="S9"):
            cluster.run(c, rel, V, stages, tmp=work)
    with pytest.raises(FileNotFoundError, match="S9"):
        tokens.run(c, rel, stages)
    with pytest.raises(FileNotFoundError, match="S9"):
        postings.run(c, rel, stages)
    with pytest.raises(FileNotFoundError, match="S9"):
        dim.run(c, rel, stages, tmp=work)
    with pytest.raises(FileNotFoundError, match="S9"):
        catalog.build_manifest(c, rel, stages)
    # nothing was written for the release, not even its definition: a typo is not made permanent
    assert not [p for p in c.rglob("*") if f"release={rel}" in p.as_posix()]


def test_b25_release_def_is_written_once_and_verified(ladder, tmp_path):
    c = tmp_path / "c"
    shutil.copytree(ladder["b"], c)
    path = catalog.release_def_path(c, "rS0.order")
    assert catalog.define_release(c, "rS0.order", ["S1", "S0"]) == ["S1", "S0"]
    t = pq.read_table(path, partitioning=None)
    assert t.to_pylist() == [{"position": 0, "stage": "S1"}, {"position": 1, "stage": "S0"}]
    assert conformance_problems(t, table("release_def")) == []
    before = (path.read_bytes(), path.stat().st_mtime_ns)
    assert catalog.define_release(c, "rS0.order", ["S1", "S0"]) == ["S1", "S0"]  # the same stages: no-op
    assert (path.read_bytes(), path.stat().st_mtime_ns) == before
    assert catalog.release_stages(c, "rS0.order") == ["S1", "S0"]  # in the order defined, not sorted
    for other in (["S0", "S1"], ["S0"], ["S1", "S0", "S2"]):
        with pytest.raises(ValueError, match="release_def"):
            catalog.define_release(c, "rS0.order", other)
    assert path.read_bytes() == before[0]
    with pytest.raises(FileNotFoundError, match="rS0.undefined"):
        catalog.release_stages(c, "rS0.undefined")
    with pytest.raises(TypeError):
        catalog.define_release(c, "rS0.bad", "S0")
    for bad in ([], ["S0", "S0"], [""]):
        with pytest.raises(ValueError):
            catalog.define_release(c, "rS0.bad", bad)
    # a builder given other stages than the release's definition refuses before writing
    with pytest.raises(ValueError, match="release_def"):
        tokens.run(c, REL, ["S0", "S1"])
    with pytest.raises(ValueError, match="release_def"):
        dim.run(c, "rS0.order", ["S0", "S1"], tmp=ladder["work"])
    with pytest.raises(ValueError, match="release_def"):
        catalog.build_manifest(c, REL, ["S0", "S1"])
    assert not (c / f"dim/contig_dim/release=rS0.order").exists()


def test_b25_cli_release_steps_define_once_and_read_the_definition_after(ladder, tmp_path, capsys):
    c = tmp_path / "c"
    shutil.copytree(ladder["b"], c)
    # each release has its own clustering ledger, ledger/release=<R>/step=cluster (card B-28)
    rel = "rCLI.1"
    base = ["--root", str(c), "--release", rel]
    assert cli.main(["release", "tokens", *base]) == 2  # not defined, and no --stages
    assert "release_def" in capsys.readouterr().err
    assert cli.main(["release", "define", *base, "--stages", "S0,S9"]) == 2  # S9 was never ingested
    assert "S9" in capsys.readouterr().err
    assert not catalog.release_def_path(c, rel).exists()
    assert cli.main(["release", "define", *base, "--stages", "S0"]) == 0
    assert cli.main(["release", "define", *base, "--stages", "S0"]) == 0  # the same stages: no-op
    assert cli.main(["release", "dim", *base, "--stages", "S0,S1"]) == 2  # not the defined stages
    assert "release_def" in capsys.readouterr().err
    for step in (("cluster",) if MMSEQS else ()) + ("tokens", "postings", "dim", "manifest"):
        assert cli.main(["release", step, *base, "--tmp", str(ladder["work"])]) == 0, step
    man = json.loads((c / f"manifests/release={rel}.json").read_text())
    assert man["release_id"] == rel and man["stages"] == ["S0"]
    assert not [f for f in man["files"] if f["partitions"].get("stage") == "S1"]
    assert {"token_dict", "tape_tokens", "token_postings", "contig_dim", "release_def"} <= {
        f["table"] for f in man["files"] if f["partitions"].get("release") == rel}
    # the CLI's tokens equal the library call's for the same stages (the R1 build)
    for t in ("token_dict", "tape_tokens"):
        mine = sorted(f["sha256"] for f in man["files"] if f["table"] == t and f["partitions"].get("release") == rel)
        r1 = sorted(f["sha256"] for f in ladder["m2"]["files"] if f["table"] == t and f["partitions"].get("release") == REL)
        assert mine == r1, t
