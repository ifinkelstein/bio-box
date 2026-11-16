"""Card B-49 (DEC-45 (2)): a release takes an earlier release's clustering when its proteins are the same.

Release R1 over stage S0 is built by `build-stage` with its clustering in parts (`--cluster-parts 1`: a clustering in
one process writes no plan records, so it cannot be a source, `biobox.l2.cluster` "Reuse"). From a copy of that root:

- R2 over the same stage, by `build-stage`: its cluster step runs no clustering call and copies R1's files;
- R2 over the same stage, by `build-stage --no-reuse`: clustered. This is the oracle: every file of the release (the
  cluster tables, tokens, postings, dimension tables, `release_def`, R2's own cluster plan records) is the same by
  sha256 in both roots, the manifests are the same but for the fields that name the run (the root, and the clustering
  ledger, which names R1 in the reused build), and `release verify` passes on both;
- R3 over S0 and S1 (a stage added), and R4 over a stage S0 whose one protein's sequence was changed in place: nothing
  is reused and the clustering runs. A brute-force comparison of the releases' (protein_uid, seq) sets decides "the
  same proteins" independently of the plan records;
- a clustering stopped after its plan step is finished, not reused over; R1's files are never changed; a source file
  whose bytes are not its run record's is refused, and nothing is written for the release.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog, cli, ids, io, ledger, provenance, runstep
from biobox import stage as bs
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l2 import cluster as l2cluster
from biobox.registry.schemas import table

sys.path.insert(0, str(Path(__file__).parent))
import test_builders_out as c08  # noqa: E402

direct = c08.direct  # the module-scoped direct build of S0 (its profile, inputs and HMM library are used here)

pytestmark = [
    pytest.mark.r2,
    pytest.mark.covers("B-49"),
    pytest.mark.skipif(c08.MMSEQS is None, reason="build-stage clusters with MMseqs2 (in the image)"),
]

R1, R2, R3, R4 = "rB49.1", "rB49.2", "rB49.3", "rB49.4"
PARTS = 1  # every per-query call in one part: the plan and its records, at the fewest steps
PLAN = "l2/cluster_plan/"
TABLE_DIRS = ("l2/cluster_member/", "l2/cluster_member_aln/", "l2/member_cluster/")
RELEASE_STEPS = ["catalog", "cluster", "dim", "manifest", "postings", "tokens"]


class Killed(Exception):
    """Raised by the injected hook: the build process dies between two steps."""


def build(direct, root: Path, release: str, stage: str = "S0", stages=("S0",), **kw) -> bs.Build:
    w = direct["work"]
    return bs.Build(root=Path(root), stage=stage, release=release, stages=tuple(stages), genomes=w.genomes(stage),
                    objects=w.objects(stage), libraries=(bs.Library("test", "1", w.lib, "gathering"),), threads=1,
                    cluster_parts=PARTS, **kw)


def runner(prof, root, stage, step, cmd, **kw) -> int:
    """run-step at a fixed hour and load, as `test_build_stage.runner`."""
    return runstep.run_step(prof, root, stage, step, cmd, clock=lambda: c08.NOON_MONDAY, load=lambda: 0.0, **kw)


def copy(src: Path, dst: Path) -> Path:
    shutil.copytree(src, dst, symlinks=True)
    return dst


def run_names(root: Path, release: str) -> list[str]:
    d = Path(root) / "runs" / f"stage={release}"
    return sorted(p.name.removesuffix(".json") for p in d.glob("*.json")) if d.is_dir() else []


def clustering_calls(root: Path, release: str) -> list[str]:
    """The run records of the release's clustering calls (each module call, cut and part is a run-step)."""
    return [n for n in run_names(root, release) if n.startswith("cluster.") and n != l2cluster.PLAN_STEP]


def protein_set(root: Path, stages: list[str]) -> set[tuple[bytes, str]]:
    """Brute force, independent of the plan records: every (protein_uid, seq) row of the stages' `protein_seq`."""
    out = set()
    for s in stages:
        for f in sorted((Path(root) / f"l0seq/protein_seq/gene_call_version={V}/stage={s}").rglob("*.parquet")):
            out |= {(r["protein_uid"], r["seq"]) for r in pq.read_table(f, columns=["protein_uid", "seq"]).to_pylist()}
    return out


def reused_dirs(release: str) -> tuple[str, ...]:
    return (f"{PLAN}release={release}/{l2cluster.REUSE_RECORD}.json", f"ledger/release={release}/step=cluster/",
            "manifests/", "catalog/")


def release_tree(root: Path, release: str) -> dict[str, str]:
    """sha256 of every file of the root but run records, the manifest and catalog (which name the root) and the two
    files that name the source of a reused clustering: its reuse record and its clustering ledger."""
    return {k: v for k, v in c08.tree(root).items() if not k.startswith(reused_dirs(release))}


def ledger_rows(root: Path, release: str) -> list[dict]:
    return pq.ParquetFile(ledger.release_path(root, release, "cluster")).read().to_pylist()


def manifest_but_the_run(root: Path, release: str) -> dict:
    man = catalog.read_manifest(catalog.manifest_path(root, release))
    led = ledger.release_path(Path(""), release, "cluster").as_posix()
    return {**{k: v for k, v in man.items() if k != "root"}, "files": [f for f in man["files"] if f["path"] != led]}


def cut_of(direct, root: Path, release: str) -> dict:
    return bs.cluster_cut(direct["prof"], build(direct, root, release))


@pytest.fixture(scope="module")
def first(direct, tmp_path_factory):
    """R1 over S0, clustered in parts by build-stage."""
    root = tmp_path_factory.mktemp("b49-first") / "root"
    assert bs.build_stage(direct["prof"], build(direct, root, R1), run=runner) == 0
    return root, c08.tree(root)


@pytest.fixture(scope="module")
def reused(direct, first, tmp_path_factory):
    """R2 over S0 by build-stage, from a copy of R1's root: the earlier clustering is taken."""
    root = copy(first[0], tmp_path_factory.mktemp("b49-reused") / "root")
    assert bs.build_stage(direct["prof"], build(direct, root, R2), run=runner) == 0
    return root


@pytest.fixture(scope="module")
def clustered(direct, first, tmp_path_factory):
    """R2 over S0 by build-stage --no-reuse, from a copy of R1's root: the oracle."""
    root = copy(first[0], tmp_path_factory.mktemp("b49-clustered") / "root")
    assert bs.build_stage(direct["prof"], build(direct, root, R2, reuse=False), run=runner) == 0
    return root


def test_b49_the_same_proteins_take_the_earlier_clustering_and_no_clustering_call_runs(direct, first, reused):
    root = reused
    assert protein_set(root, ["S0"]) == protein_set(first[0], ["S0"])  # the brute-force oracle: the same proteins
    assert clustering_calls(root, R2) == [] and run_names(root, R2) == RELEASE_STEPS
    rec = json.loads((root / f"runs/stage={R2}/cluster.json").read_text())
    i = rec["argv"].index("reuse")
    assert rec["argv"][i - 2:i] == ["release", "cluster"] and rec["argv"][rec["argv"].index("--from") + 1] == R1
    assert rec["cpu_min"] < 1.0  # seconds, not a clustering (the fixture's clustering in one process is one CPU-minute)
    # the reuse record names the source and both plans: the plan records' own hashes, which name the release, differ;
    # their input hashes, which do not, are equal
    doc = l2cluster.read_record(root, R2, l2cluster.REUSE_RECORD)
    raw1 = (root / PLAN / f"release={R1}/{l2cluster.PLAN_STEP}.json").read_bytes()
    raw2 = (root / PLAN / f"release={R2}/{l2cluster.PLAN_STEP}.json").read_bytes()
    assert doc["from"] == R1 and doc["release"] == R2
    assert doc["from_plan_sha256"] == hashlib.sha256(raw1).hexdigest() != doc["plan_sha256"] == \
        hashlib.sha256(raw2).hexdigest()
    h = l2cluster.inputs_sha256(json.loads(raw1))
    assert doc["from_inputs_sha256"] == doc["inputs_sha256"] == h == l2cluster.inputs_sha256(json.loads(raw2))
    # the run record lists the reuse record among its outputs
    assert f"{PLAN}release={R2}/{l2cluster.REUSE_RECORD}.json" in {o["path"] for o in rec["outputs"]}
    # the ledger: the clustering's rows, then one per level that names the source and the plans' input hash
    rows = ledger_rows(root, R2)
    ok, more = [r for r in rows if r["reason_code"] == "ok"], [r for r in rows if r["reason_code"] != "ok"]
    assert [r["object_key"] for r in ok] == [b"nr50", b"nr90"] == [r["object_key"] for r in more]
    assert {r["reason_code"] for r in more} == {f"reused:{R1}:{h}"}
    assert [{k: v for k, v in r.items() if k != "reason_code"} for r in more] == \
        [{k: v for k, v in r.items() if k != "reason_code"} for r in ok]
    assert ledger.check_balance(root, "S0", ["cluster"], release=R2) == []


def test_b49_reuse_equals_a_clustering_of_the_release_file_for_file(direct, first, reused, clustered):
    """Every file of R2 by reuse is the file of R2 by clustering, by sha256: the cluster tables, tokens, postings,
    dimension tables, release_def and R2's cluster plan records (rewritten for R2 as its clustering wrote them)."""
    got, want = release_tree(reused, R2), release_tree(clustered, R2)
    assert sorted(got) == sorted(want)
    assert {k: v for k, v in got.items() if got[k] != want[k]} == {}
    mine = [k for k in got if f"release={R2}" in k]
    for d in (*TABLE_DIRS, "l3/token_dict/", "l3/tape_tokens/", "l3/token_postings/", "dim/contig_dim/", f"{PLAN}"):
        assert any(k.startswith(d) for k in mine), d
    plan_records = [k for k in mine if k.startswith(f"{PLAN}release={R2}/")]
    assert len(plan_records) == len(run_names(clustered, R2)) - len(RELEASE_STEPS)  # each step's record but `cluster`
    # the manifests are the same but for the fields that name the run; both releases verify
    assert manifest_but_the_run(reused, R2) == manifest_but_the_run(clustered, R2)
    for root in (reused, clustered):
        assert cli.main(["release", "verify", "--root", str(root), "--release", R2]) == 0
        assert provenance.check_provenance(catalog.read_manifest(catalog.manifest_path(root, R2))) == []
    # the clustering ledger: the clustered build's rows are the reused build's rows of the clustering
    assert ledger_rows(clustered, R2) == [r for r in ledger_rows(reused, R2) if r["reason_code"] == "ok"]
    assert (clustered / PLAN / f"release={R2}/{l2cluster.REUSE_RECORD}.json").exists() is False


def test_b49_the_copies_are_files_of_their_own_and_r1_is_not_touched(first, reused):
    root, before = reused, first[1]
    after = c08.tree(root)
    assert {k: after[k] for k in before} == before  # every file of R1's root, R1's own among them, unchanged
    copies = [k for k in after if k.startswith(TABLE_DIRS) and f"release={R2}/" in k]
    assert len(copies) >= 3
    for k in copies:
        p = root / k
        src = k.replace(f"release={R2}/", f"release={R1}/")
        assert after[k] == before[src]  # byte for byte the source's
        assert not p.is_symlink() and p.stat().st_nlink == 1 and p.stat().st_ino != (root / src).stat().st_ino


def test_b49_the_choice_is_deterministic_and_no_reuse_forces_a_clustering(direct, first, reused, clustered, tmp_path):
    a = first[0]
    cut = cut_of(direct, a, R2)
    assert l2cluster.reuse_source(a, R2, V, ["S0"], cut, tmp=tmp_path) == R1
    # --no-reuse: the clustering ran, every call a run-step of its own, and nothing names a source
    assert l2cluster.PLAN_STEP in run_names(clustered, R2) and len(clustering_calls(clustered, R2)) > 40
    assert {r["reason_code"] for r in ledger_rows(clustered, R2)} == {"ok"}
    # two finished clusterings of the same proteins: both qualify, the first by name is taken
    inputs = l2cluster.plan_inputs(reused, R3, V, ["S0"], cut)
    assert l2cluster.reuse_sources(reused, R3, inputs) == [R1, R2]
    assert l2cluster.reuse_source(reused, R3, V, ["S0"], cut, tmp=tmp_path) == R1
    # another cut is another plan: a clustering of R3 at two parts would write other records, so nothing is taken
    two = l2cluster.cut_of(parts=2)
    assert l2cluster.reuse_sources(reused, R3, l2cluster.plan_inputs(reused, R3, V, ["S0"], two)) == []
    # the cut build-stage gives the plan step is the one it compares
    assert l2cluster.read_record(clustered, R2, l2cluster.PLAN_STEP)["cut"] == cut


def test_b49_a_release_over_more_stages_reuses_nothing_and_is_clustered(direct, first, tmp_path):
    root = copy(first[0], tmp_path / "root")
    assert bs.build_stage(direct["prof"], build(direct, root, R3, stage="S1", stages=("S0", "S1")), run=runner) == 0
    assert protein_set(root, ["S0"]) < protein_set(root, ["S0", "S1"])  # the brute-force oracle: more proteins
    assert l2cluster.PLAN_STEP in run_names(root, R3) and len(clustering_calls(root, R3)) > 40
    assert not (root / PLAN / f"release={R3}/{l2cluster.REUSE_RECORD}.json").exists()
    assert {r["reason_code"] for r in ledger_rows(root, R3)} == {"ok"}
    assert cli.main(["release", "verify", "--root", str(root), "--release", R3]) == 0


def _change_one_sequence(root: Path) -> None:
    """Stage S0's `protein_seq` with one protein's sequence changed (a residue appended; its uid, digest and length
    made again; the rows sorted by uid again): the constructed case of a changed protein."""
    [f] = sorted((root / f"l0seq/protein_seq/gene_call_version={V}/stage=S0").rglob("part-*.parquet"))
    rows = pq.read_table(f).to_pylist()
    r = rows[len(rows) // 2]
    s = r["seq"] + "W"
    rows[len(rows) // 2] = {**r, "protein_uid": ids.protein_uid(s), "protein_digest": hashlib.sha256(s.encode()).digest(),
                            "length": len(s), "seq": s}
    rows.sort(key=lambda x: x["protein_uid"])
    f.unlink()
    io.write_table(pa.Table.from_pylist(rows, schema=table("protein_seq").schema), table("protein_seq"), f)


def test_b49_a_changed_sequence_reuses_nothing_and_is_clustered(direct, first, tmp_path):
    root = copy(first[0], tmp_path / "root")
    before = protein_set(root, ["S0"])
    _change_one_sequence(root)
    assert protein_set(root, ["S0"]) != before and len(protein_set(root, ["S0"])) == len(before)
    cut = cut_of(direct, root, R4)
    assert l2cluster.reuse_source(root, R4, V, ["S0"], cut, tmp=tmp_path) is None
    with pytest.raises(ValueError, match=R1):
        l2cluster.reuse(root, R4, V, ["S0"], tmp_path, source=R1, parts=PARTS)
    assert not [p for p in root.rglob("*") if f"release={R4}" in p.as_posix()]  # nothing written for R4
    res = l2cluster.run_split(root, R4, V, ["S0"], tmp=tmp_path / "w", parts=PARTS, threads=1)
    assert res["proteins"] == len(before)
    calls = [p.name for p in (root / PLAN / f"release={R4}").glob("*.json")]
    assert len(calls) > 40 and f"{l2cluster.REUSE_RECORD}.json" not in calls  # the clustering ran
    new = {k for k in c08.tree(root) if k.startswith(TABLE_DIRS) and f"release={R4}/" in k}
    old = {k.replace(f"release={R1}/", f"release={R4}/"): v for k, v in first[1].items() if k.startswith(TABLE_DIRS)}
    assert new and any(c08.tree(root)[k] != old.get(k) for k in new)  # other clusters than R1's


def test_b49_a_clustering_stopped_after_its_plan_is_finished_not_reused(direct, first, clustered, tmp_path):
    root = copy(first[0], tmp_path / "root")

    def kill(point: str) -> None:
        if point == f"after:{l2cluster.PLAN_STEP}":
            raise Killed(point)

    with pytest.raises(Killed):
        bs.build_stage(direct["prof"], build(direct, root, R2, reuse=False), run=runner, _hook=kill)
    assert run_names(root, R2) == [l2cluster.PLAN_STEP]
    assert bs.build_stage(direct["prof"], build(direct, root, R2), run=runner) == 0  # reuse allowed: not taken
    assert run_names(root, R2) == run_names(clustered, R2)
    assert release_tree(root, R2) == release_tree(clustered, R2)
    assert not (root / PLAN / f"release={R2}/{l2cluster.REUSE_RECORD}.json").exists()


def test_b49_a_source_file_that_is_not_its_run_records_is_refused_and_nothing_is_written(direct, first, tmp_path):
    root = copy(first[0], tmp_path / "root")
    rec = json.loads((root / f"runs/stage={R1}/cluster.json").read_text())
    victim = next(o for o in rec["outputs"] if o["path"].startswith(f"l2/member_cluster/release={R1}/"))
    p = root / victim["path"]
    data = bytearray(p.read_bytes())
    data[len(data) // 2] ^= 0xFF  # the same size, another sha256
    p.write_bytes(bytes(data))
    with pytest.raises(ValueError, match="sha256"):
        l2cluster.reuse(root, R2, V, ["S0"], tmp_path / "w", source=R1, parts=PARTS)
    assert not [q for q in root.rglob("*") if f"release={R2}" in q.as_posix()]


def test_b49_the_cli_refuses_a_source_that_does_not_qualify(direct, first, tmp_path, capsys):
    root = copy(first[0], tmp_path / "root")
    args = ["release", "cluster", "reuse", "--root", str(root), "--release", R2, "--stages", "S0",
            "--tmp", str(tmp_path / "w")]
    assert cli.main([*args, "--from", R1, "--parts", "2"]) == 2  # another cut
    assert "cut" in capsys.readouterr().err
    assert cli.main([*args, "--from", "rB49.none", "--parts", str(PARTS)]) == 2  # no such clustering
    assert "rB49.none" in capsys.readouterr().err
    assert not [q for q in root.rglob("*") if f"release={R2}" in q.as_posix()]
    out = tmp_path / "out"
    assert cli.main([*args, "--from", R1, "--parts", str(PARTS), "--out", str(out)]) == 0
    res = json.loads(capsys.readouterr().out)["result"]
    assert res["reused_from"] == R1 and res["files"] >= 3
