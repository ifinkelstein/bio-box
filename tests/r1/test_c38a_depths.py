"""Ring R1: card C-38a — stage depth chosen at build, derived for a release (DEC-54).

Covers:
  (1) `stage.steps` at each depth (`stored`, `called`, `annotated`), default depth when `--depth` is
      omitted (`annotated` with `--library`, else `called`), optional `--release` on stage-only builds,
      and exit-2 refusals before any step runs (`--depth annotated` without `--library`, `--depth called`
      with `--library` or `--ncrna-cm`, `--depth stored` with `--library`, `--ncrna-cm`, `--crispr` or any
      release option, and `--stages` / `--cluster-parts` without `--release`).
  (3) `catalog.stage_depths(manifest)` and `catalog.describe_depths` on constructed manifests (stored,
      called, single/multi-library annotated, ncRNA, mixed release) and on exact path copies of `rM.2`
      (all seven stages `annotated pfam-38.2`) and `rN.1` (all three stages `called`).
  (4) A release over a stage at depth `stored` (ingested, no `contig_tapes`) is refused (exit 2), naming
      the stage and its depth `stored`, with nothing written.
  (5) `stage.stage_steps` and `stage.stage_complete` take the depth (`stored`, `called`, `annotated`),
      derive the stage step list from `stage.steps` (never a copy), and check run records, `.publishing.json`
      intent files, and recorded outputs.
"""

from __future__ import annotations

import dataclasses
import json
import shutil
from datetime import datetime
from pathlib import Path

import pytest

from biobox import catalog, cli, fasta, ids, runstep
from biobox import stage as bs
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l1 import crispr, ncrna
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.profile import Profile

pytestmark = [pytest.mark.covers("C-38a", "D-01", "LY-10")]

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
NOON = datetime(2026, 11, 9, 12, 0)

# Exact path list of `manifests/release=rM.2.json` on wilkcomp02 (/tmp/if2453/r10/root, 7 stages, 158 files,
# 154 artifacts; measured 2026-11-18). `classify` / `classify_artifact` reconstruct the exact table/artifact
# and partitions of every entry.
RM2_STAGES = ["S1", "P1", "P2", "P3", "P4", "P5", "P6"]
RM2_SHARDS = {"S1": 1, "P1": 1, "P2": 9, "P3": 10, "P4": 10, "P5": 10, "P6": 10}
RM2_PROTEIN_PARTS = {"S1": 1, "P1": 1, "P2": 2, "P3": 2, "P4": 2, "P5": 2, "P6": 2}
RM2_LABEL_STAGES = ("P3", "P4", "P5", "P6")

RN1_STAGES = ["S1", "P1", "P2"]
RN1_SHARDS = {"S1": 1, "P1": 1, "P2": 9}
RN1_PROTEIN_PARTS = {"S1": 1, "P1": 1, "P2": 2}


def _reconstruct_manifest(release_id: str, stages: list[str], file_paths: list[str],
                          artifact_paths: list[str]) -> dict:
    files = []
    for rel in sorted(file_paths):
        hit = catalog.classify(rel)
        assert hit is not None, rel
        spec, parts = hit
        files.append({"path": rel, "table": spec.name, "partitions": parts})
    artifacts = []
    for rel in sorted(artifact_paths):
        hit = catalog.classify_artifact(rel)
        assert hit is not None, rel
        spec, parts, _ = hit
        artifacts.append({"path": rel, "artifact": spec.name, "partitions": parts})
    return {"manifest_version": 3, "release_id": release_id, "stages": list(stages),
            "files": files, "artifacts": artifacts}


def rm2_manifest_copy() -> dict:
    tv = "t-ca596dd85040"
    art = ["l0seq/store_params.json"]
    for s, n in RM2_SHARDS.items():
        for k in range(n):
            for suf in ("", ".fai", ".gzi"):
                art.append(f"l0seq/seq_shards/stage={s}/shard-{k:04d}.fa.gz{suf}")
    paths = [
        "dim/contig_dim/release=rM.2/part-0.parquet",
        "dim/genome_dim/release=rM.2/part-0.parquet",
        "l3/tape_version_def/tape_version=t-ca596dd85040/part-0.parquet",
        "ledger/release=rM.2/step=cluster/part-0.parquet",
        "release_def/release=rM.2/part-0.parquet",
        *(f"l2/cluster_member/release=rM.2/part-{i}.parquet" for i in range(4)),
        *(f"l2/cluster_member_aln/release=rM.2/part-{i}.parquet" for i in range(3)),
        *(f"l2/member_cluster/release=rM.2/part-{i}.parquet" for i in range(4)),
        *(f"l3/tape_tokens/tape_version={tv}/release=rM.2/part-{i}.parquet" for i in range(2)),
        *(f"l3/token_dict/release=rM.2/part-{i}.parquet" for i in range(3)),
        *(f"l3/token_postings/release=rM.2/tape_version={tv}/token_bucket={b:x}/part-0.parquet" for b in range(16)),
        *(f"l0seq/contig_label/stage={s}/part-0.parquet" for s in RM2_LABEL_STAGES),
    ]
    for s in RM2_STAGES:
        src = "inphared" if s == "S1" else "ncbi_genomes"
        paths += [
            f"l0/l0_manifest/stage={s}/part-0.parquet",
            f"l0seq/contig_index/stage={s}/part-0.parquet",
            f"l0seq/contig_occurrence/stage={s}/source={src}/part-0.parquet",
            f"l0seq/genome_input/stage={s}/part-0.parquet",
            f"l0seq/igr/igr_version={tv}-igr1/stage={s}/part-0.parquet",
            *(f"l0seq/protein_seq/gene_call_version={V}/stage={s}/protein_prefix=all/part-{i}.parquet"
              for i in range(RM2_PROTEIN_PARTS[s])),
            f"l0seq/seq_locator/stage={s}/part-0.parquet",
            f"l1/domains/library=pfam-38.2/stage={s}/part-0.parquet",
            f"l1/feature/producer={V}/stage={s}/part-0.parquet",
            f"l1/feature_protein/gene_call_version={V}/stage={s}/part-0.parquet",
            f"l1/protein_occurrence/gene_call_version={V}/stage={s}/part-0.parquet",
            f"l3/tapes/tape_version={tv}/stage={s}/part-0.parquet",
            f"ledger/stage={s}/step=domains-pfam-38.2/part-0.parquet",
            f"ledger/stage={s}/step=genecall/part-0.parquet",
            f"ledger/stage={s}/step=ingest/part-0.parquet",
            f"ledger/stage={s}/step=tape-{tv}/part-0.parquet",
        ]
    man = _reconstruct_manifest("rM.2", RM2_STAGES, paths, art)
    assert len(man["files"]) == 158 and len(man["artifacts"]) == 154
    return man


def rn1_manifest_copy() -> dict:
    tv = "t-ca596dd85040"
    art = ["l0seq/store_params.json"]
    for s, n in RN1_SHARDS.items():
        for k in range(n):
            for suf in ("", ".fai", ".gzi"):
                art.append(f"l0seq/seq_shards/stage={s}/shard-{k:04d}.fa.gz{suf}")
    paths = [
        "dim/contig_dim/release=rN.1/part-0.parquet",
        "l2/cluster_member/release=rN.1/part-0.parquet",
        "l2/cluster_member_aln/release=rN.1/part-0.parquet",
        f"l3/tape_tokens/tape_version={tv}/release=rN.1/part-0.parquet",
        f"l3/tape_version_def/tape_version={tv}/part-0.parquet",
        "l3/token_dict/release=rN.1/part-0.parquet",
        f"l3/token_postings/release=rN.1/tape_version={tv}/token_bucket=all/part-0.parquet",
        "ledger/release=rN.1/step=cluster/part-0.parquet",
        "release_def/release=rN.1/part-0.parquet",
    ]
    for s in RN1_STAGES:
        src = "inphared" if s == "S1" else "ncbi_genomes"
        paths += [
            f"l0/l0_manifest/stage={s}/part-0.parquet",
            f"l0seq/contig_index/stage={s}/part-0.parquet",
            f"l0seq/contig_occurrence/stage={s}/source={src}/part-0.parquet",
            f"l0seq/genome_input/stage={s}/part-0.parquet",
            f"l0seq/igr/igr_version={tv}-igr1/stage={s}/part-0.parquet",
            *(f"l0seq/protein_seq/gene_call_version={V}/stage={s}/protein_prefix=all/part-{i}.parquet"
              for i in range(RN1_PROTEIN_PARTS[s])),
            f"l0seq/seq_locator/stage={s}/part-0.parquet",
            f"l1/feature/producer={V}/stage={s}/part-0.parquet",
            f"l1/feature_protein/gene_call_version={V}/stage={s}/part-0.parquet",
            f"l1/protein_occurrence/gene_call_version={V}/stage={s}/part-0.parquet",
            f"l3/tapes/tape_version={tv}/stage={s}/part-0.parquet",
            f"ledger/stage={s}/step=genecall/part-0.parquet",
            f"ledger/stage={s}/step=ingest/part-0.parquet",
            f"ledger/stage={s}/step=tape-{tv}/part-0.parquet",
        ]
    return _reconstruct_manifest("rN.1", RN1_STAGES, paths, art)


def _profile(tmp_path: Path) -> Profile:
    return Profile(
        max_threads_day=4, max_threads_night=6, day_start_hour=8, day_end_hour=20,
        mem_scope="4G", require_mem_cap=False, nice_monolithic=0, nice_fanout=0,
        cpu_min_per_process=300, tmp=str(tmp_path / "tmp"),
        durable=str(tmp_path / "durable"), bulk=str(tmp_path / "bulk"),
        small_file_floor_bytes=1000,
    )


def dummy_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    g = tmp_path / "genomes.json"
    o = tmp_path / "objects.json"
    hmm = tmp_path / "lib.hmm"
    fna = FIX / "NC_001416.1.fna.gz"
    cu = ids.contig_uid(fasta.read(fna)[0].seq)
    g.write_text(json.dumps([{"contig_uids": [cu.hex()], "genetic_code": 11, "kind": "virus"}]))
    o.write_text(json.dumps([{"path": str(fna), "source_name": "refseq", "source_release": "1",
                              "license": "PD", "redistributable": True}]))
    hmm.write_text("HMMER3/f [3.4 | Aug 2023]\nNAME x\nLENG 10\n//\n")
    return g, o, hmm


# --- (1) step lists and flag refusals per depth ----------------------------------------------------------


@pytest.mark.covers("C-38a", "D-01")
def test_c38a_step_list_at_each_depth_and_default_depth(tmp_path, monkeypatch):
    if shutil.which("minced") is None:  # GitHub Actions R1 has no MinCED binary (card E-12)
        monkeypatch.setattr(crispr, "minced_version", lambda: "0.4.2")
    prof = _profile(tmp_path)
    g, o, hmm = dummy_inputs(tmp_path)
    root = tmp_path / "root"
    pfam = bs.Library("pfam", "38.2", hmm, "gathering")

    # Default when --depth is omitted: called without a library, annotated with one (preserves pre-C-38a behavior)
    b_default_called = bs.Build(root=root, stage="S0", release="rT.1", stages=("S0",), genomes=g, objects=o)
    assert bs.depth_of(b_default_called) == "called"
    assert [s.name for s in bs.steps(prof, b_default_called)] == [
        "inputs", "ingest", "genecall", "features", "igr", "tapes",
        "cluster", "tokens", "postings", "dim", "manifest", "catalog",
    ]

    b_default_ann = bs.Build(root=root, stage="S0", release="rT.1", stages=("S0",), genomes=g, objects=o,
                             libraries=(pfam,))
    assert bs.depth_of(b_default_ann) == "annotated"
    assert [s.name for s in bs.steps(prof, b_default_ann)] == [
        "inputs", "ingest", "genecall", "features", "igr", "tapes", "domains-pfam-38.2",
        "cluster", "tokens", "postings", "dim", "manifest", "catalog",
    ]

    # Explicit --depth stored (stage-only: no release)
    b_stored = bs.Build(root=root, stage="S0", release=None, genomes=g, objects=o, depth="stored")
    assert bs.depth_of(b_stored) == "stored"
    got_stored = bs.steps(prof, b_stored)
    assert [(s.run_stage, s.name) for s in got_stored] == [("S0", "inputs"), ("S0", "ingest")]

    # Explicit --depth called without --release, with and without --crispr
    b_called = bs.Build(root=root, stage="S0", release=None, genomes=g, objects=o, depth="called")
    assert [s.name for s in bs.steps(prof, b_called)] == ["inputs", "ingest", "genecall", "features", "igr", "tapes"]

    b_called_crispr = dataclasses.replace(b_called, crispr=True)
    assert [s.name for s in bs.steps(prof, b_called_crispr)] == [
        "inputs", "ingest", "genecall", "features", bs.crispr_name(), "igr", "tapes",
    ]

    # Explicit --depth annotated without --release, with --crispr and pfam (runs near-Cas producer after Pfam scan)
    b_ann_stage = bs.Build(root=root, stage="S0", release=None, genomes=g, objects=o, depth="annotated",
                           libraries=(pfam,), crispr=True)
    assert [s.name for s in bs.steps(prof, b_ann_stage)] == [
        "inputs", "ingest", "genecall", "features", bs.crispr_name(), "igr", "tapes",
        "domains-pfam-38.2", bs.crispr_near_name(pfam),
    ]


@pytest.mark.covers("C-38a", "D-01")
@pytest.mark.parametrize("bad_kwargs, match", [
    ({"depth": "bogus", "release": None}, "--depth"),
    ({"depth": "annotated", "release": None, "libraries": ()}, "annotated"),
    ({"depth": "called", "release": None, "libraries": "pfam"}, "--depth called"),
    ({"depth": "called", "release": None, "ncrna": True}, "--depth called"),
    ({"depth": "stored", "release": None, "libraries": "pfam"}, "--depth stored"),
    ({"depth": "stored", "release": None, "ncrna": True}, "--depth stored"),
    ({"depth": "stored", "release": None, "crispr": True}, "--depth stored"),
    ({"depth": "stored", "release": "rT.1", "stages": ("S0",)}, "--depth stored"),
    ({"depth": "stored", "release": None, "stages": ("S0",)}, "stored"),
    ({"depth": "stored", "release": None, "cluster_parts": 2}, "stored"),
    ({"depth": "called", "release": None, "stages": ("S0",)}, "--release"),
    ({"depth": "called", "release": None, "cluster_parts": 2}, "--release"),
])
def test_c38a_incompatible_depth_options_refused_before_any_step(tmp_path, capsys, bad_kwargs, match):
    prof = _profile(tmp_path)
    g, o, hmm = dummy_inputs(tmp_path)
    cm = tmp_path / "rfam.cm"
    cm.write_text("INFERNAL1/a\n")
    root = tmp_path / "root"
    kw = dict(bad_kwargs)
    libs = (bs.Library("pfam", "38.2", hmm, "gathering"),) if kw.pop("libraries", ()) == "pfam" else ()
    ncrna_cm = cm if kw.pop("ncrna", False) else None
    b = bs.Build(root=root, stage="S0", genomes=g, objects=o, libraries=libs, ncrna_cm=ncrna_cm, **kw)

    with pytest.raises(ValueError, match=match):
        bs.steps(prof, b)

    called = []
    code = bs.build_stage(prof, b, run=lambda *a, **k: called.append(a) or 0, clock=lambda: NOON, load=lambda: 0.0)
    err = capsys.readouterr().err
    assert code == runstep.EXIT_REFUSED
    assert match in err
    assert called == [] and not root.exists()


@pytest.mark.covers("C-38a")
def test_c38a_cli_build_stage_parses_depth_and_refuses_bad_combinations(tmp_path, capsys):
    import yaml

    prof = _profile(tmp_path)
    pf = tmp_path / "profile.yaml"
    pf.write_text(yaml.safe_dump(dataclasses.asdict(prof)))
    g, o, hmm = dummy_inputs(tmp_path)
    root = tmp_path / "root"
    base = ["build-stage", "--profile", str(pf), "--root", str(root), "--stage", "S0",
            "--genomes", str(g), "--objects", str(o)]

    # --depth annotated without --library -> exit 2, nothing created
    assert cli.main([*base, "--depth", "annotated"]) == 2
    assert "annotated" in capsys.readouterr().err and not root.exists()

    # --depth called with --library -> exit 2, nothing created
    assert cli.main([*base, "--depth", "called", "--library", "pfam", "38.2", str(hmm), "gathering"]) == 2
    assert "--depth called" in capsys.readouterr().err and not root.exists()

    # --depth stored with --release -> exit 2, nothing created
    assert cli.main([*base, "--depth", "stored", "--release", "rT.1", "--stages", "S0"]) == 2
    assert "--depth stored" in capsys.readouterr().err and not root.exists()


# --- (3) catalog.stage_depths and describe_depths --------------------------------------------------------


@pytest.mark.covers("C-38a", "LY-10")
def test_c38a_stage_depths_on_rm2_and_rn1_and_mixed_manifests():
    rm2 = rm2_manifest_copy()
    got_rm2 = catalog.stage_depths(rm2)
    assert list(got_rm2) == RM2_STAGES
    assert got_rm2 == {s: "annotated pfam-38.2" for s in RM2_STAGES}
    assert catalog.describe_depths(got_rm2) == "S1,P1,P2,P3,P4,P5,P6 annotated pfam-38.2"

    rn1 = rn1_manifest_copy()
    got_rn1 = catalog.stage_depths(rn1)
    assert list(got_rn1) == RN1_STAGES
    assert got_rn1 == {s: "called" for s in RN1_STAGES}
    assert catalog.describe_depths(got_rn1) == "S1,P1,P2 called"

    # Mixed manifest: S1,P1 annotated (with two libraries and ncRNA on S1), SP1 called, MG0 stored
    tv = "t-ca596dd85040"
    mixed = _reconstruct_manifest(
        "r26.11.2",
        ["S1", "P1", "SP1", "MG0"],
        [
            f"l0seq/contig_index/stage=S1/part-0.parquet",
            f"l0seq/protein_seq/gene_call_version={V}/stage=S1/protein_prefix=all/part-0.parquet",
            f"l3/tapes/tape_version={tv}/stage=S1/part-0.parquet",
            "l1/domains/library=pfam-38.2/stage=S1/part-0.parquet",
            "l1/domains/library=defense-1.0/stage=S1/part-0.parquet",
            "l1/ncrna/producer=infernal-1.1.5-rfam15-prok2/stage=S1/part-0.parquet",
            f"l0seq/contig_index/stage=P1/part-0.parquet",
            f"l0seq/protein_seq/gene_call_version={V}/stage=P1/protein_prefix=all/part-0.parquet",
            f"l3/tapes/tape_version={tv}/stage=P1/part-0.parquet",
            "l1/domains/library=pfam-38.2/stage=P1/part-0.parquet",
            f"l0seq/contig_index/stage=SP1/part-0.parquet",
            f"l0seq/protein_seq/gene_call_version={V}/stage=SP1/protein_prefix=all/part-0.parquet",
            f"l3/tapes/tape_version={tv}/stage=SP1/part-0.parquet",
            f"l0seq/contig_index/stage=MG0/part-0.parquet",
        ],
        [],
    )
    depths = catalog.stage_depths(mixed)
    assert depths == {
        "S1": "annotated defense-1.0,pfam-38.2,ncrna:infernal-1.1.5-rfam15-prok2",
        "P1": "annotated pfam-38.2",
        "SP1": "called",
        "MG0": "stored",
    }
    assert catalog.describe_depths(depths) == (
        "S1 annotated defense-1.0,pfam-38.2,ncrna:infernal-1.1.5-rfam15-prok2; "
        "P1 annotated pfam-38.2; SP1 called; MG0 stored"
    )

    # Empty synthetic manifest (as in W-37 unit tests) degrades cleanly
    empty = {"release_id": "rEmpty.1", "stages": ["S0"], "files": [], "artifacts": []}
    assert catalog.describe_depths(catalog.stage_depths(empty)) == "-"


# --- (4) release over a stored stage refused -------------------------------------------------------------


@pytest.mark.covers("C-38a")
def test_c38a_release_over_stored_stage_is_refused_naming_stage_and_depth(tmp_path, capsys):
    root = tmp_path / "root"
    fna = FIX / "NC_001416.1.fna.gz"
    seq = fasta.read(fna)[0].seq
    cu = ids.contig_uid(seq)
    write_genome_input(root, "MG0", [GenomeInput((cu,), 11, kind="virus")])
    ingest([SourceObject(fna, "refseq", "1", "PD", True)], root, "MG0")

    assert catalog.stage_depth_on_disk(root, "MG0") == "stored"
    assert catalog.stage_depth_on_disk(root, "MISSING") is None

    with pytest.raises(ValueError, match=r"MG0.*stored"):
        catalog.require_release_tapes(root, "r26.11.9", ["MG0"])

    with pytest.raises(FileNotFoundError, match="not ingested"):
        catalog.require_release_tapes(root, "r26.11.9", ["MISSING"])

    for step in ("define", "cluster", "tokens", "postings", "dim"):
        extra = ["--without-clusters"] if step == "tokens" else []
        code = cli.main(["release", step, "--root", str(root), "--release", "r26.11.9", "--stages", "MG0", *extra])
        err = capsys.readouterr().err
        assert code == 2, (step, err)
        assert "MG0" in err and "stored" in err, (step, err)
    assert not catalog.release_def_path(root, "r26.11.9").exists()


# --- (5) stage completeness helper for the release train (card C-46) -------------------------------------


def _write_run_record(root: Path, stage: str, name: str, outputs: list[str] = (), exit_code: int = 0) -> None:
    d = root / "runs" / f"stage={stage}"
    d.mkdir(parents=True, exist_ok=True)
    recs = []
    for rel in outputs:
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        if not p.exists():
            p.write_bytes(b"ok")
        recs.append({"path": rel, "bytes": p.stat().st_size})
    (d / f"{name}.json").write_text(json.dumps({"stage": stage, "step": name, "exit_code": exit_code, "outputs": recs}))


@pytest.mark.covers("C-38a")
def test_c38a_stage_complete_judges_by_build_stage_step_list_for_the_given_depth(tmp_path, monkeypatch):
    if shutil.which("minced") is None:  # GitHub Actions R1 has no MinCED binary (card E-12)
        monkeypatch.setattr(crispr, "minced_version", lambda: "0.4.2")
    prof = _profile(tmp_path)
    root = tmp_path / "root"

    # Derived from stage.steps (never a copy): stored -> [inputs, ingest]
    assert [s.name for s in bs.stage_steps("S0", depth="stored", prof=prof)] == ["inputs", "ingest"]
    assert [s.name for s in bs.stage_steps("S0", depth="called", prof=prof, crispr=True)] == [
        "inputs", "ingest", "genecall", bs.crispr_name(), "igr", "tapes",
    ]
    assert [s.name for s in bs.stage_steps("S0", depth="annotated", prof=prof,
                                           libraries=[("pfam", "38.2")], crispr=True)] == [
        "inputs", "ingest", "genecall", bs.crispr_name(), "igr", "tapes",
        "domains-pfam-38.2", f"crispr-{crispr.near_cas_producer('38.2')}",
    ]

    # Empty stage is incomplete at every depth
    assert not bs.stage_complete(root, "S0", "stored")
    assert bs.stage_incomplete(root, "S0", "stored")

    # Populate stored records + outputs
    _write_run_record(root, "S0", "inputs", ["l0seq/genome_input/stage=S0/part-0.parquet"])
    assert not bs.stage_complete(root, "S0", "stored")
    _write_run_record(root, "S0", "ingest", ["l0seq/contig_index/stage=S0/part-0.parquet"])
    assert bs.stage_complete(root, "S0", "stored")
    assert bs.stage_incomplete(root, "S0", "stored") == []

    # A stored stage is NOT complete at depth `called` or `annotated`
    assert not bs.stage_complete(root, "S0", "called")
    missing_called = bs.stage_incomplete(root, "S0", "called")
    assert any("genecall" in m for m in missing_called) and any("tapes" in m for m in missing_called)

    # An intent file (.publishing.json) on any step or batch makes the stage incomplete
    intent = root / "runs" / "stage=S0" / "ingest.publishing.json"
    intent.write_text("{}")
    assert not bs.stage_complete(root, "S0", "stored")
    assert any("publishing" in m or "intent" in m for m in bs.stage_incomplete(root, "S0", "stored"))
    intent.unlink()
    assert bs.stage_complete(root, "S0", "stored")

    sub_intent = root / "runs" / "stage=S0" / "ingest.batch-00000.publishing.json"
    sub_intent.write_text("{}")
    assert not bs.stage_complete(root, "S0", "stored")
    sub_intent.unlink()

    # Populate remaining `called` steps (with crispr=True)
    for step in ("genecall", bs.crispr_name(), "igr", "tapes"):
        _write_run_record(root, "S0", step, [f"out/{step}.parquet"])
    assert bs.stage_complete(root, "S0", "called", crispr=True)
    assert not bs.stage_complete(root, "S0", "annotated", libraries=[("pfam", "38.2")], crispr=True)

    # Populate `annotated` steps (domains + near-Cas CRISPR)
    for step in ("domains-pfam-38.2", f"crispr-{crispr.near_cas_producer('38.2')}"):
        _write_run_record(root, "S0", step, [f"out/{step}.parquet"])
    assert bs.stage_complete(root, "S0", "annotated", libraries=[("pfam", "38.2")], crispr=True)

    # Missing output file invalidates completeness
    (root / "out/tapes.parquet").unlink()
    assert not bs.stage_complete(root, "S0", "annotated", libraries=[("pfam", "38.2")], crispr=True)

    # Unknown depth raises ValueError
    with pytest.raises(ValueError, match="depth"):
        bs.stage_complete(root, "S0", "bogus")
