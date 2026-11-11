"""C-08: every builder writes to a separate output root, so each runs as a `bio-box run-step` step.

A step command reads the data root (`BIOBOX_ROOT`) and writes its outputs under `BIOBOX_STAGING`;
run-step publishes them (its module docstring is the contract). Every builder takes `out=` for
that, and the CLI's `bio-box stage ...` and `bio-box release ...` subcommands default `--out` to
`$BIOBOX_STAGING`. Each step here is a real subprocess, `python -m biobox.cli ...`, started by
`runstep.run_step`, as `tests/r2/test_ingest_restart.py` does for ingest.

The oracle is the direct build: every builder called as the existing code calls it (no `out=`),
one after the other, in one root, with a snapshot of the root after each step. Then:

- each builder alone, as a run-step on a copy of the snapshot before it, gives the snapshot after
  it, file for file by sha256;
- the whole chain inputs -> ingest -> genecall -> igr -> tapes -> domains -> cluster -> tokens ->
  postings -> dim -> manifest -> catalog, each its own run-step in one root, gives the direct root's
  files, the same manifest hash and the same catalog file list;
- a second stage and a second release built the same way change no file of the first.

The manifest and the catalog hold the data root's absolute path (the manifest's `root`; the
catalog's `catalog_info.root` and every file path in its views). Built in staging, they must still
name the final root, never the staging directory: that is checked on the published files, and the
published catalog is opened and read from.

A builder refuses, before any work, when its partitions already hold a file in the root or in `out`;
its scratch directory may lie under neither; a definition row (`tape_version_def`, `release_def`)
found in the root is verified and not written again into `out`.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import sys
from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import duckdb
import pyarrow.parquet as pq
import pyhmmer
import pytest
import yaml

from biobox import catalog, cli, dim, fasta, ids, profile, runstep
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l1 import domains, igr
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l1.genecall import run as call_genes
from biobox.l3 import postings, tapes, tokens
from biobox.registry.schemas import table

pytestmark = [pytest.mark.r2, pytest.mark.covers("D-01", "LY-10")]

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
MMSEQS = shutil.which("mmseqs")
NOON_MONDAY = datetime(2026, 11, 9, 12, 0)
REL, REL2 = "rS0.c08", "rS1.c08"
# S0: lambda and M. genitalium (circular, code 4, single-mode training); S1: T4
STAGE_OBJECTS = {"S0": ("NC_001416.1", "NC_000908.2"), "S1": ("NC_000866.4",)}
STAGE_STEPS = ["inputs", "ingest", "genecall", "igr", "tapes", "domains"]
RELEASE_STEPS = (["cluster"] if MMSEQS else []) + ["tokens", "postings", "dim", "manifest", "catalog"]
NAMES = STAGE_STEPS + RELEASE_STEPS
RECORDS = ("manifests/", "catalog/")  # files that hold the root's absolute path
REFUSED = r"writes its partitions once|clean stage directory|written once"
DEFS = ("l3/tape_version_def/", "release_def/")  # verified when present, never refused


def _objects(stage: str) -> list[dict]:
    return [o for o in CONTRACT["objects"] if o["accession"] in STAGE_OBJECTS[stage]]


def source_objects(stage: str) -> list[SourceObject]:
    return [SourceObject(FIX / o["file"], CONTRACT["source_name"], CONTRACT["source_release"], CONTRACT["license"], True)
            for o in _objects(stage)]


def genome_inputs(stage: str) -> list[GenomeInput]:
    out = []
    for o in _objects(stage):
        cu = ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)
        out.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"],
                               circular=(cu,) if o["topology"] == "circular" else ()))
    return out


def write_library(root: Path, path: Path) -> None:
    """Six HMMs built from the longest lambda proteins (as test_micro_s0_release_def does)."""
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


@dataclass
class Work:
    """What every step reads from outside the root: scratch, the HMM library, the CLI's input files."""
    base: Path

    @property
    def scratch(self) -> Path:
        return self.base / "scratch"

    @property
    def lib(self) -> Path:
        return self.base / "lib.hmm"

    def genomes(self, stage: str) -> Path:
        return self.base / f"genomes-{stage}.json"

    def objects(self, stage: str) -> Path:
        return self.base / f"objects-{stage}.json"

    def prepare(self) -> None:
        self.scratch.mkdir(parents=True, exist_ok=True)
        for stage in STAGE_OBJECTS:
            self.genomes(stage).write_text(json.dumps([
                {"contig_uids": [u.hex() for u in g.contig_uids], "genetic_code": g.genetic_code, "kind": g.kind,
                 "circular": [u.hex() for u in g.circular], "assembly_quality": g.assembly_quality}
                for g in genome_inputs(stage)]))
            self.objects(stage).write_text(json.dumps([
                {"path": str(o.path), "source_name": o.source_name, "source_release": o.source_release,
                 "license": o.license, "redistributable": o.redistributable} for o in source_objects(stage)]))


def _kw(out) -> dict:
    """`out=` only when given: the oracle calls each builder exactly as the existing code does."""
    return {} if out is None else {"out": out}


@dataclass
class Step:
    run_stage: str  # run-step's --stage (the release, for a release step)
    name: str
    argv: Callable[[Path], list[str]]  # the CLI arguments, given the root
    direct: Callable[..., object]  # direct(root, out=None)


def stage_steps(stage: str, w: Work) -> list[Step]:
    def a(sub: str, *extra: str) -> Callable[[Path], list[str]]:
        return lambda r: ["stage", sub, "--root", str(r), "--stage", stage, *extra]

    s = str(w.scratch)

    def scan(r, out=None):
        if not w.lib.exists():
            write_library(r, w.lib)
        return domains.run(r, stage, V, "test", "1", w.lib, "gathering", **_kw(out))

    return [
        Step(stage, "inputs", a("inputs", "--genomes", str(w.genomes(stage))),
             lambda r, out=None: write_genome_input(r, stage, genome_inputs(stage), **_kw(out))),
        Step(stage, "ingest", a("ingest", "--objects", str(w.objects(stage)), "--tmp", s),
             lambda r, out=None: ingest(source_objects(stage), r, stage, tmp=w.scratch, **_kw(out))),
        Step(stage, "genecall", a("genecall", "--tmp", s),
             lambda r, out=None: call_genes(r, stage, tmp=w.scratch, **_kw(out))),
        Step(stage, "igr", a("igr", "--tmp", s), lambda r, out=None: igr.run(r, stage, tmp=w.scratch, **_kw(out))),
        Step(stage, "tapes", a("tapes", "--tmp", s), lambda r, out=None: tapes.run(r, stage, tmp=w.scratch, **_kw(out))),
        Step(stage, "domains", a("domains", "--library", "test", "--library-version", "1", "--hmm", str(w.lib),
                                 "--cutoff", "gathering"), scan),
    ]


def release_steps(release: str, stages: list[str], w: Work) -> list[Step]:
    csv, s = ",".join(stages), str(w.scratch)

    def a(sub: str, *extra: str) -> Callable[[Path], list[str]]:
        return lambda r: ["release", sub, "--root", str(r), "--release", release, *extra]

    def manifest(r, out=None):
        man = catalog.build_manifest(r, release, stages)
        if out is None:
            return catalog.write_manifest(man, catalog.manifest_path(r, release))
        return catalog.write_manifest(man, out=out)

    def cat(r, out=None):
        man = catalog.read_manifest(catalog.manifest_path(r, release))
        if out is None:
            return catalog.generate_catalog(man, Path(r) / f"catalog/release={release}.duckdb")
        return catalog.generate_catalog(man, root=r, out=out)

    steps = []
    if MMSEQS:
        from biobox.l2 import cluster

        steps.append(Step(release, "cluster", a("cluster", "--stages", csv, "--tmp", s),
                          lambda r, out=None: cluster.run(r, release, V, stages, tmp=w.scratch, **_kw(out))))
    return steps + [
        Step(release, "tokens", a("tokens", "--stages", csv, "--tmp", s),
             lambda r, out=None: tokens.run(r, release, stages, tmp=w.scratch, **_kw(out))),
        Step(release, "postings", a("postings", "--stages", csv, "--tmp", s),
             lambda r, out=None: postings.run(r, release, stages, tmp=w.scratch, **_kw(out))),
        Step(release, "dim", a("dim", "--stages", csv, "--tmp", s),
             lambda r, out=None: dim.run(r, release, stages, tmp=w.scratch, **_kw(out))),
        Step(release, "manifest", a("manifest"), manifest),
        Step(release, "catalog", a("catalog"), cat),
    ]


def tree(root: Path) -> dict[str, str]:
    """sha256 of every file under the root except run-step's `runs/`."""
    root = Path(root)
    if not root.exists():
        return {}
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file() and p.relative_to(root).parts[0] != "runs"}


def data(t: dict[str, str]) -> dict[str, str]:
    return {k: v for k, v in t.items() if not k.startswith(RECORDS)}


def profile_at(base: Path) -> profile.Profile:
    data_ = {"max_threads_day": 4, "max_threads_night": 6, "mem_scope": "4G",
             "tmp": str(base / "tmp"), "durable": str(base / "durable"), "bulk": str(base / "bulk"),
             "nice_monolithic": 19, "nice_fanout": 15, "day_start_hour": 8, "day_end_hour": 20,
             "cpu_min_per_process": 300, "small_file_floor_bytes": 1000}
    base.mkdir(parents=True, exist_ok=True)
    (base / "profile.yaml").write_text(yaml.safe_dump(data_))
    return profile.load(base / "profile.yaml")


def run_cli(prof: profile.Profile, root: Path, step: Step) -> int:
    cmd = [sys.executable, "-m", "biobox.cli", *step.argv(root)]
    return runstep.run_step(prof, root, step.run_stage, step.name, cmd, clock=lambda: NOON_MONDAY, load=lambda: 0.0)


@pytest.fixture(scope="module")
def direct(tmp_path_factory):
    """The direct build of S0 and release REL, with the root as it stood after each step."""
    d = tmp_path_factory.mktemp("c08")
    w = Work(d / "work")
    w.prepare()
    root = d / "direct"
    steps = stage_steps("S0", w) + release_steps(REL, ["S0"], w)
    assert [s.name for s in steps] == NAMES
    snaps = []
    for k, step in enumerate(steps):
        step.direct(root)
        snaps.append(d / f"snap-{k:02d}-{step.name}")
        shutil.copytree(root, snaps[-1], symlinks=True)
    return {"root": root, "work": w, "steps": steps, "snaps": snaps,
            "manifest": catalog.read_manifest(catalog.manifest_path(root, REL)), "prof": profile_at(d / "prof")}


def _snap_before(direct, k: int) -> Path | None:
    return direct["snaps"][k - 1] if k else None


def check_manifest(root: Path, release: str, want: dict) -> dict:
    """The published manifest names the final root and is otherwise the direct build's."""
    man = catalog.read_manifest(catalog.manifest_path(root, release))
    assert man["root"] == str(Path(root).resolve())
    assert catalog.manifest_hash(man) == catalog.manifest_hash(want)
    assert {k: v for k, v in man.items() if k != "root"} == {k: v for k, v in want.items() if k != "root"}
    return man


def check_catalog(root: Path, release: str, man: dict, staging: Path) -> None:
    """The published catalog embeds the final root, not the staging directory, and reads from it."""
    root = Path(root).resolve()
    path = root / f"catalog/release={release}.duckdb"
    con = duckdb.connect(str(path), read_only=True)
    try:
        info = dict(con.execute("SELECT key, value FROM catalog_info").fetchall())
        sql = " ".join(s for (s,) in con.execute("SELECT sql FROM duckdb_views() WHERE NOT internal").fetchall())
        files = con.execute('SELECT path, "table", sha256 FROM catalog_files ORDER BY path').fetchall()
    finally:
        con.close()
    assert info["root"] == str(root)
    assert info["manifest_hash"] == catalog.manifest_hash(man) and info["release_id"] == release
    paths = re.findall(r"'(/[^']*\.parquet)'", sql)
    assert paths and all(p.startswith(str(root) + "/") for p in paths)
    assert {p.removeprefix(str(root) + "/") for p in paths} == {f["path"] for f in man["files"]}
    assert str(Path(staging).resolve()) not in sql + json.dumps(info)
    assert files == sorted((f["path"], f["table"], f["sha256"]) for f in man["files"])
    # opened as a reader opens it, the catalog reads the published root's files and shards
    o = _objects(man["stages"][0])[0]
    seq = fasta.read(FIX / o["file"])[0].seq
    c = catalog.open_catalog(path)
    try:
        assert catalog.seq_reader(c).root == root
        got = c.execute("SELECT seq_get(?, 0, 60, 1)", [ids.contig_uid(seq)]).fetchone()[0]
        assert got.encode() == seq[:60]
        assert c.execute("SELECT count(*) FROM contig_index").fetchone()[0] == sum(
            f["rows"] for f in man["files"] if f["table"] == "contig_index")
    finally:
        c.close()


def published(root: Path, step: Step) -> list[str]:
    rec = json.loads((Path(root) / f"runs/stage={step.run_stage}/{step.name}.json").read_text())
    return [o["path"] for o in rec["outputs"]]


# --- each builder alone ------------------------------------------------------------------------


@pytest.mark.parametrize("name", NAMES)
def test_each_builder_as_a_run_step_equals_the_direct_call(direct, tmp_path, name):
    k = NAMES.index(name)
    step, before, after = direct["steps"][k], _snap_before(direct, k), direct["snaps"][k]
    root = tmp_path / "root"
    if before is not None:
        shutil.copytree(before, root, symlinks=True)
    assert run_cli(direct["prof"], root, step) == 0
    got, want = tree(root), tree(after)
    assert data(got) == data(want)
    assert sorted(got) == sorted(want)
    new = sorted(set(want) - set(tree(before) if before is not None else {}))
    assert published(root, step) == new  # only the step's own files, nothing rewritten
    if name == "manifest":
        check_manifest(root, REL, direct["manifest"])
    if name == "catalog":  # from the manifest copied in with the snapshot; the root is this one
        check_catalog(root, REL, direct["manifest"], Path(direct["prof"].tmp) / "staging")


# --- the whole chain, step by step -------------------------------------------------------------


@pytest.fixture(scope="module")
def chain(direct, tmp_path_factory):
    root = tmp_path_factory.mktemp("c08-chain") / "root"
    for step in direct["steps"]:
        assert run_cli(direct["prof"], root, step) == 0, step.name
    return root


def test_the_whole_chain_under_run_step_equals_the_direct_build(direct, chain):
    got, want = tree(chain), tree(direct["root"])
    assert sorted(got) == sorted(want)
    assert data(got) == data(want)
    man = check_manifest(chain, REL, direct["manifest"])
    assert catalog.verify_manifest(man) == []
    check_catalog(chain, REL, man, Path(direct["prof"].tmp) / "staging")
    # every step published its own files and nothing else: the run records partition the root
    outs = [p for step in direct["steps"] for p in published(chain, step)]
    assert sorted(outs) == sorted(got) and len(outs) == len(set(outs))


@pytest.mark.covers("R-09")
def test_a_second_stage_and_release_change_no_file_of_the_first(direct, chain, tmp_path):
    w, prof = direct["work"], direct["prof"]
    later = stage_steps("S1", w) + release_steps(REL2, ["S0", "S1"], w)
    root = tmp_path / "root"
    shutil.copytree(chain, root, symlinks=True)
    first = tree(root)
    for step in later:
        assert run_cli(prof, root, step) == 0, step.name
    after = tree(root)
    assert {k: after[k] for k in first} == first  # every first-stage and first-release file unchanged
    assert all("stage=S0" not in p for step in later for p in published(root, step))
    # the same steps called directly on the direct store give the same files and the same manifest
    oracle = tmp_path / "direct"
    shutil.copytree(direct["root"], oracle, symlinks=True)
    for step in later:
        step.direct(oracle)
    assert data(after) == data(tree(oracle))
    man = check_manifest(root, REL2, catalog.read_manifest(catalog.manifest_path(oracle, REL2)))
    check_catalog(root, REL2, man, Path(prof.tmp) / "staging")


# --- refusals ----------------------------------------------------------------------------------


def _targets(direct, k: int) -> list[str]:
    """The files step k writes into partitions it owns (definition rows are verified, never refused)."""
    before = _snap_before(direct, k)
    new = set(tree(direct["snaps"][k])) - set(tree(before) if before is not None else {})
    return sorted(p for p in new if not p.startswith(DEFS))


@pytest.mark.parametrize("name", NAMES)
def test_a_builder_refuses_before_any_work_when_its_partition_holds_files(direct, tmp_path, name):
    k = NAMES.index(name)
    step, before, after = direct["steps"][k], _snap_before(direct, k), direct["snaps"][k]
    targets = _targets(direct, k)
    assert targets
    # the step's files are in the root: an empty output directory stays empty
    root, out = tmp_path / "root", tmp_path / "out"
    shutil.copytree(after, root, symlinks=True)
    with pytest.raises(FileExistsError, match=REFUSED):
        step.direct(root, out=out)
    assert tree(out) == {} and tree(root) == tree(after)
    # one of the step's files is in the output directory (a killed attempt's): nothing else is written
    shutil.rmtree(root)
    if before is not None:
        shutil.copytree(before, root, symlinks=True)
    planted = targets[-1]
    (out / planted).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(after / planted, out / planted)
    with pytest.raises(FileExistsError, match=REFUSED):
        step.direct(root, out=out)
    assert sorted(tree(out)) == [planted]
    assert tree(root) == (tree(before) if before is not None else {})


@pytest.mark.parametrize("name", ["ingest", "genecall", "igr", "tapes", *(["cluster"] if MMSEQS else []),
                                  "tokens", "postings", "dim"])
def test_scratch_may_lie_under_neither_the_root_nor_the_output(direct, tmp_path, name):
    k = NAMES.index(name)
    step, before = direct["steps"][k], _snap_before(direct, k)
    root, out = tmp_path / "root", tmp_path / "out"
    if before is not None:
        shutil.copytree(before, root, symlinks=True)
    for bad in (out / "scratch", root / "scratch"):
        bad.mkdir(parents=True, exist_ok=True)
        argv = step.argv(root)
        argv[argv.index("--tmp") + 1] = str(bad)
        assert cli.main([*argv, "--out", str(out)]) == 2, (name, bad)
        assert tree(out) == {} and tree(bad) == {}
        shutil.rmtree(bad)
    assert tree(root) == (tree(before) if before is not None else {})


def test_definition_rows_in_the_root_are_verified_not_written_into_out(direct, tmp_path):
    # tapes after igr: the tape version's definition is the root's; tokens after cluster (or the first
    # release step): the release's definition is the root's. Neither is copied into `out`.
    tv_def = f"l3/tape_version_def/tape_version={tapes.tape_version([V])}/part-0.parquet"
    rel_def = f"release_def/release={REL}/part-0.parquet"
    k_tapes, k_post = NAMES.index("tapes"), NAMES.index("postings")
    root, out = tmp_path / "root", tmp_path / "out"
    shutil.copytree(direct["snaps"][k_tapes - 1], root, symlinks=True)
    assert (root / tv_def).exists()
    direct["steps"][k_tapes].direct(root, out=out)
    assert not (out / tv_def).exists() and (out / f"ledger/stage=S0/step=tape-{tapes.tape_version([V])}").is_dir()
    shutil.rmtree(root)
    shutil.rmtree(out)
    shutil.copytree(direct["snaps"][k_post - 1], root, symlinks=True)
    assert (root / rel_def).exists()
    direct["steps"][k_post].direct(root, out=out)
    assert not (out / "release_def").exists() and (out / f"l3/token_postings/release={REL}").is_dir()
    # a definition that disagrees is refused before any work, with `out` or without
    with pytest.raises(ValueError, match="release_def"):
        catalog.define_release(root, REL, ["S0", "S1"], out=tmp_path / "o2")
    assert not (tmp_path / "o2").exists()
    # with no definition in the root, the release's first step writes it into `out`, once
    fresh = tmp_path / "fresh"
    assert catalog.define_release(root, "rS0.new", ["S0"], out=fresh) == ["S0"]
    assert (fresh / "release_def/release=rS0.new/part-0.parquet").exists()
    assert not catalog.release_def_path(root, "rS0.new").exists()
    assert catalog.define_release(root, "rS0.new", ["S0"], out=fresh) == ["S0"]  # verified in `out`: no-op
    with pytest.raises(ValueError, match="release_def"):
        catalog.define_release(root, "rS0.new", ["S0", "S1"], out=fresh)


def test_out_defaults_to_staging_when_set_and_to_the_root_otherwise(direct, tmp_path, monkeypatch, capsys):
    w = direct["work"]
    root, staging = tmp_path / "root", tmp_path / "staging"
    gi = f"{table('genome_input').location}/stage=S9/part-0.parquet"
    args = ["stage", "inputs", "--root", str(root), "--stage", "S9", "--genomes", str(w.genomes("S0"))]
    monkeypatch.setenv("BIOBOX_STAGING", str(staging))
    assert cli.main(args) == 0
    assert sorted(tree(staging)) == [gi] and tree(root) == {}
    monkeypatch.delenv("BIOBOX_STAGING")
    assert cli.main(args) == 0  # no staging: into the root
    assert sorted(tree(root)) == [gi]
    assert cli.main([*args, "--out", str(tmp_path / "elsewhere")]) == 2  # the root holds the partition
    assert not (tmp_path / "elsewhere").exists()
    assert "refusing to overwrite" in capsys.readouterr().err
