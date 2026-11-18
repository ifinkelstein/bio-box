"""Ring R2: card C-38a — stage depth at build (`stored`, `called`, `annotated`) and mixed-depth releases (DEC-54).

Covers:
  (2) A stage built `--depth called` (no `--release`) and later scanned with the library — both by
      resuming `build-stage --depth annotated --library ... --release <R> --stages S0` and by running
      `stage domains` + `release` steps — produces files byte-identical to building `annotated` at once
      (`test_builders_out.direct`), and the cross-host micro-S0 fixture (`scripts/stress/build_fixture.py`)
      keeps its manifest hash (`28ffa9808ecd4829e52f6892f94f7b85`) with `stage_depths == {"S0": "called"}`.
  (3) `catalog.stage_depths(manifest)` on test manifests; on a mixed release over `S0` (`annotated test-1`)
      and `S1` (`called`, T4 unscanned even though the HMM file holds T4 models), `release verify` prints
      `depths: S0 annotated test-1; S1 called` on stdout (with no C-19 warning on stderr), `bio-box releases`
      prints the identical `depths` column, and `token_dict`, `tape_tokens` and `token_postings` family
      tokens equal an independent brute-force count written separately over `S0`'s `l1_domains` Parquet file
      and the stages' `contig_tapes` Parquet files (AGENTS §5.2).
  (4) A stage built with `build-stage --depth stored` completes at `stored` (`stage_complete`), and both
      `bio-box release` and `bio-box build-stage` refuse a release over it (exit 2), naming the stage and
      its depth `stored`, with no release file published.
"""

from __future__ import annotations

import dataclasses
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pyarrow.parquet as pq
import pyhmmer
import pytest
import yaml

from biobox import catalog, cli, fasta, ids, runstep
from biobox import stage as bs
from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION as V

sys.path.insert(0, str(Path(__file__).parent))
import test_builders_out as c08  # noqa: E402

direct = c08.direct

pytestmark = [
    pytest.mark.r2,
    pytest.mark.covers("C-38a", "D-01", "D-03", "LY-10"),
    pytest.mark.skipif(c08.MMSEQS is None, reason="build-stage clusters with MMseqs2 (in the image)"),
]

REPO = Path(__file__).resolve().parents[2]
FIXTURE_HASH = "28ffa9808ecd4829e52f6892f94f7b85"


def runner(prof, root, stage, step, cmd, **kw) -> int:
    return runstep.run_step(prof, root, stage, step, cmd, clock=lambda: c08.NOON_MONDAY, load=lambda: 0.0, **kw)


def profile_path(prof, path: Path) -> Path:
    path.write_text(yaml.safe_dump(dataclasses.asdict(prof)))
    return path


def record_names(root: Path) -> set[str]:
    runs = Path(root) / "runs"
    return {p.relative_to(runs).as_posix().removesuffix(".json") for p in sorted(runs.rglob("*.json"))} \
        if runs.exists() else set()


# --- (2) called then scanned == annotated at once; cross-host fixture hash unchanged ----------------------


@pytest.mark.covers("C-38a", "D-01", "D-03")
def test_c38a_stage_built_called_then_scanned_equals_annotated_at_once(direct, tmp_path):
    prof, w = direct["prof"], direct["work"]
    pf = profile_path(prof, tmp_path / "profile.yaml")

    # Path A: CLI `build-stage --depth called` (no --release), then `build-stage --depth annotated --library ... --release`
    root_a = tmp_path / "called_then_build_stage" / "root"
    assert cli.main([
        "build-stage", "--profile", str(pf), "--root", str(root_a), "--stage", "S0",
        "--depth", "called", "--genomes", str(w.genomes("S0")), "--objects", str(w.objects("S0")),
        "--threads", "1",
    ]) == 0
    assert record_names(root_a) == {f"stage=S0/{n}" for n in ("inputs", "ingest", "genecall", "igr", "tapes")}
    assert bs.stage_complete(root_a, "S0", "called")
    assert not bs.stage_complete(root_a, "S0", "annotated", libraries=[("test", "1")])
    assert not (root_a / "l1/domains").exists() and not (root_a / "manifests").exists()

    called_files = c08.tree(root_a)
    assert cli.main([
        "build-stage", "--profile", str(pf), "--root", str(root_a), "--stage", "S0",
        "--depth", "annotated", "--release", c08.REL, "--stages", "S0",
        "--genomes", str(w.genomes("S0")), "--objects", str(w.objects("S0")),
        "--library", "test", "1", str(w.lib), "gathering", "--threads", "1",
    ]) == 0
    assert bs.stage_complete(root_a, "S0", "annotated", libraries=[("test", "1")])
    # Earlier stage files were not rewritten (AGENTS §3.10)
    after_a = c08.tree(root_a)
    assert {k: after_a[k] for k in called_files} == called_files
    # Byte-identical to the direct annotated-at-once build, including manifest hash
    want = c08.tree(direct["root"])
    assert sorted(after_a) == sorted(want)
    assert c08.data(after_a) == c08.data(want)
    man_a = c08.check_manifest(root_a, c08.REL, direct["manifest"])
    assert catalog.stage_depths(man_a) == {"S0": "annotated test-1"}

    # Path B: `build-stage --depth called`, then standalone `stage domains` + `release` steps under run-step
    root_b = tmp_path / "called_then_direct_scan" / "root"
    b_called = bs.Build(root=root_b, stage="S0", release=None, genomes=w.genomes("S0"),
                        objects=w.objects("S0"), depth="called", threads=1)
    assert bs.build_stage(prof, b_called, run=runner, clock=lambda: c08.NOON_MONDAY, load=lambda: 0.0) == 0
    ann_plan = bs.steps(prof, bs.Build(root=root_b, stage="S0", release=c08.REL, stages=("S0",),
                                       genomes=w.genomes("S0"), objects=w.objects("S0"), depth="annotated",
                                       libraries=(bs.Library("test", "1", w.lib, "gathering"),), threads=1))
    for step in ann_plan:
        if step.argv is not None and not (root_b / "runs" / f"stage={step.run_stage}" / f"{step.name}.json").exists():
            if step.scratch is not None:
                Path(step.scratch).mkdir(parents=True, exist_ok=True)
            assert runner(prof, root_b, step.run_stage, step.name,
                          [sys.executable, "-m", step.module, *step.argv], threads=1) == 0
    after_b = c08.tree(root_b)
    assert sorted(after_b) == sorted(want)
    assert c08.data(after_b) == c08.data(want)
    c08.check_manifest(root_b, c08.REL, direct["manifest"])


@pytest.mark.covers("C-38a", "D-05")
def test_c38a_micro_s0_fixture_manifest_hash_unchanged_and_depth_is_called(tmp_path):
    fx_root = tmp_path / "fixture_root"
    out = subprocess.run(
        [sys.executable, str(REPO / "scripts/stress/build_fixture.py"), str(REPO), str(fx_root), "rX.fixture"],
        capture_output=True, text=True, check=True,
    ).stdout
    report = json.loads(out)
    assert report["manifest_hash"] == FIXTURE_HASH
    man = catalog.build_manifest(fx_root, "rX.fixture", ["S0"])
    assert catalog.manifest_hash(man) == FIXTURE_HASH
    assert catalog.stage_depths(man) == {"S0": "called"}
    assert catalog.describe_depths(catalog.stage_depths(man)) == "S0 called"


# --- (3) mixed release: one annotated stage + one called stage, brute-force family token oracle -----------


def _write_lambda_and_t4_library(root_s0: Path, root_s1: Path, path: Path) -> tuple[set[str], set[str]]:
    """Build 3 HMMs from S0's longest lambda proteins and 3 HMMs from S1's longest T4 proteins.
    Returns ({lambda model_accs}, {t4 model_accs})."""
    lam = ids.contig_uid(fasta.read(c08.FIX / "NC_001416.1.fna.gz")[0].seq)
    t4 = ids.contig_uid(fasta.read(c08.FIX / "NC_000866.4.fna.gz")[0].seq)
    prots_s0 = domains.load_proteins(root_s0, V, "S0")
    occ_s0 = pq.read_table(root_s0 / f"l1/protein_occurrence/gene_call_version={V}/stage=S0/part-0.parquet").to_pylist()
    lam_uids = sorted({r["protein_uid"] for r in occ_s0 if r["contig_uid"] == lam},
                      key=lambda u: (-len(prots_s0[u]), u))[:3]
    prots_s1 = domains.load_proteins(root_s1, V, "S1")
    occ_s1 = pq.read_table(root_s1 / f"l1/protein_occurrence/gene_call_version={V}/stage=S1/part-0.parquet").to_pylist()
    t4_uids = sorted({r["protein_uid"] for r in occ_s1 if r["contig_uid"] == t4},
                     key=lambda u: (-len(prots_s1[u]), u))[:3]
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder, bg = pyhmmer.plan7.Builder(alphabet, seed=42), pyhmmer.plan7.Background(alphabet)
    lam_accs, t4_accs = set(), set()
    with open(path, "wb") as fh:
        for i, u in enumerate(lam_uids):
            seq = pyhmmer.easel.TextSequence(name=b"s", sequence=prots_s0[u].replace("*", "")).digitize(alphabet)
            hmm, _, _ = builder.build(seq, bg)
            acc = f"LAM{i:04d}"
            hmm.name, hmm.accession = f"lam{i}".encode(), acc.encode()
            hmm.cutoffs.gathering = (25.0, 25.0)
            hmm.write(fh)
            lam_accs.add(acc)
        for i, u in enumerate(t4_uids):
            seq = pyhmmer.easel.TextSequence(name=b"s", sequence=prots_s1[u].replace("*", "")).digitize(alphabet)
            hmm, _, _ = builder.build(seq, bg)
            acc = f"T4M{i:04d}"
            hmm.name, hmm.accession = f"t4m{i}".encode(), acc.encode()
            hmm.cutoffs.gathering = (200.0, 200.0)  # self-hit >1000 bits on T4 (S1), above weak cross-hits in S0
            hmm.write(fh)
            t4_accs.add(acc)
    return lam_accs, t4_accs


def _brute_force_family_oracle(root: Path, annotated_stage: str, library: str,
                               stages: tuple[str, ...]) -> tuple[dict[str, int], dict[bytes, dict], set[tuple[str, bytes, bytes]]]:
    """Independent brute-force oracle over the annotated stage's `l1_domains` Parquet files and all stages'
    `contig_tapes` Parquet files (AGENTS §5.2: never calling `biobox.l3.tokens` or `biobox.l3.postings`).

    Returns:
      - slot_counts: {model_acc: number of tape slots whose protein_uid has >=1 domain row of model_acc}
      - per_contig: {contig_uid: {"canonical": [model_acc | None, ...], "domains": [[model_acc, ...], ...]}}
      - postings_set: {(model_acc, contig_uid, feature_uid), ...} across all slots
    """
    dom_files = sorted((root / f"l1/domains/library={library}/stage={annotated_stage}").glob("*.parquet"))
    assert dom_files
    # protein_uid -> list of (dom_bitscore, model_acc, env_start)
    by_protein: dict[bytes, list[tuple[float, str, int]]] = {}
    for f in dom_files:
        for r in pq.read_table(f, columns=["protein_uid", "model_acc", "dom_bitscore", "env_start"]).to_pylist():
            by_protein.setdefault(r["protein_uid"], []).append(
                (float(r["dom_bitscore"]), str(r["model_acc"]), int(r["env_start"]))
            )

    tape_rows = []
    for s in stages:
        tfiles = sorted((root / "l3/tapes").glob(f"tape_version=*/stage={s}/*.parquet"))
        assert tfiles
        for f in tfiles:
            tape_rows.extend(pq.read_table(f, columns=["contig_uid", "feature_uids", "protein_uids"]).to_pylist())

    slot_counts: dict[str, int] = {}
    per_contig: dict[bytes, dict] = {}
    postings_set: set[tuple[str, bytes, bytes]] = set()
    for row in tape_rows:
        cu = row["contig_uid"]
        can_list: list[str | None] = []
        dom_list: list[list[str]] = []
        for fu, pu in zip(row["feature_uids"], row["protein_uids"]):
            hits = by_protein.get(pu, []) if pu is not None else []
            if not hits:
                can_list.append(None)
                dom_list.append([])
                continue
            # Canonical family: smallest (-dom_bitscore, model_acc, env_start)
            best = min(hits, key=lambda t: (-t[0], t[1], t[2]))[1]
            can_list.append(best)
            # Distinct (env_start, model_acc) ordered by (env_start, model_acc)
            distinct_env = sorted({(env_s, acc) for _, acc, env_s in hits})
            dom_list.append([acc for _, acc in distinct_env])
            for acc in {acc for _, acc, _ in hits}:
                slot_counts[acc] = slot_counts.get(acc, 0) + 1
                postings_set.add((acc, cu, fu))
        per_contig[cu] = {"canonical": can_list, "domains": dom_list}
    return slot_counts, per_contig, postings_set


@pytest.mark.covers("C-38a", "LY-10")
def test_c38a_mixed_release_over_called_and_annotated_stage_matches_brute_force_oracle(tmp_path, capsys):
    w = c08.Work(tmp_path / "work")
    w.prepare()
    prof = c08.profile_at(tmp_path / "prof")
    pf = profile_path(prof, tmp_path / "profile.yaml")
    root = tmp_path / "root"
    rel = "r26.11.1"

    # 1. Build S0 and S1 to `called` (stage-only, no --release) so we can make HMMs from both lambda and T4
    for s in ("S0", "S1"):
        assert cli.main([
            "build-stage", "--profile", str(pf), "--root", str(root), "--stage", s,
            "--depth", "called", "--genomes", str(w.genomes(s)), "--objects", str(w.objects(s)),
            "--threads", "1",
        ]) == 0
        assert bs.stage_complete(root, s, "called")

    hmm = tmp_path / "mixed.hmm"
    lam_accs, t4_accs = _write_lambda_and_t4_library(root, root, hmm)
    assert len(lam_accs) == 3 and len(t4_accs) == 3

    # 2. Annotate ONLY S0 and cut the mixed release over (S0, S1); S1 stays at depth `called` (unscanned)
    assert cli.main([
        "build-stage", "--profile", str(pf), "--root", str(root), "--stage", "S0",
        "--depth", "annotated", "--release", rel, "--stages", "S0,S1",
        "--genomes", str(w.genomes("S0")), "--objects", str(w.objects("S0")),
        "--library", "test", "1", str(hmm), "gathering", "--threads", "1",
    ]) == 0
    assert bs.stage_complete(root, "S0", "annotated", libraries=[("test", "1")])
    assert bs.stage_complete(root, "S1", "called")
    assert not bs.stage_complete(root, "S1", "annotated", libraries=[("test", "1")])
    assert not (root / "l1/domains/library=test-1/stage=S1").exists()
    capsys.readouterr()

    # 3. stage_depths, `release verify` (`depths:` line on stdout, no C-19 warning on stderr), and `bio-box releases`
    man = catalog.read_manifest(catalog.manifest_path(root, rel))
    assert catalog.stage_depths(man) == {"S0": "annotated test-1", "S1": "called"}
    want_depths = "S0 annotated test-1; S1 called"
    assert catalog.describe_depths(catalog.stage_depths(man)) == want_depths

    assert cli.main(["release", "verify", "--root", str(root), "--release", rel]) == 0
    v_io = capsys.readouterr()
    assert f"depths: {want_depths}" in v_io.out.splitlines()
    assert v_io.err == ""  # no C-19 partial-library warning on stderr

    assert cli.main(["releases", "--root", str(root)]) == 0
    r_lines = [ln.split("\t") for ln in capsys.readouterr().out.splitlines()]
    assert "depths" in r_lines[0]
    row = dict(zip(r_lines[0], r_lines[1]))
    assert row["release"] == rel and row["depths"] == want_depths

    assert cli.main(["releases", "--root", str(root), "--json"]) == 0
    jrec = json.loads(capsys.readouterr().out.strip())
    assert jrec["depths"] == want_depths
    assert jrec["stage_depths"] == {"S0": "annotated test-1", "S1": "called"}

    # 4. Brute-force oracle over S0's domain files and S0+S1's tapes
    want_counts, want_contig, want_postings = _brute_force_family_oracle(root, "S0", "test-1", ("S0", "S1"))
    assert set(want_counts) >= lam_accs  # all 3 lambda families present from S0
    assert not (set(want_counts) & t4_accs)  # T4 was not scanned (S1 is `called`), so T4 families are absent

    tdict = [r for f in sorted((root / f"l3/token_dict/release={rel}").glob("*.parquet"))
             for r in pq.read_table(f).to_pylist()]
    fam_rows = [r for r in tdict if r["level"] == "family"]
    assert {r["source"] for r in fam_rows} == {"test-1"}
    assert {r["key"]: r["n_slots"] for r in fam_rows} == want_counts
    id_to_acc = {r["token_id"]: r["key"] for r in fam_rows}

    ttok_rows = [r for f in sorted((root / "l3/tape_tokens").glob(f"tape_version=*/release={rel}/*.parquet"))
                 for r in pq.read_table(f).to_pylist()]
    assert len(ttok_rows) == len(want_contig)
    t4_cu = ids.contig_uid(fasta.read(c08.FIX / "NC_000866.4.fna.gz")[0].seq)
    for r in ttok_rows:
        cu = r["contig_uid"]
        assert [None if t is None else id_to_acc[t] for t in r["family_tokens"]] == want_contig[cu]["canonical"]
        # Map expected domain accessions to token_ids and sort by (env_start) which matches token_id tiebreak
        got_doms = [[id_to_acc[t] for t in lst] for lst in r["domain_family_tokens"]]
        assert [sorted(x) for x in got_doms] == [sorted(x) for x in want_contig[cu]["domains"]]
        if cu == t4_cu:
            # S1 (called): zero family tokens, but every protein slot still gets NR50 and NR90 cluster tokens
            assert all(t is None for t in r["family_tokens"])
            assert all(lst == [] for lst in r["domain_family_tokens"])
            assert len(r["nr90_tokens"]) > 0 and all(t is not None for t in r["nr90_tokens"])
            assert len(r["nr50_tokens"]) > 0 and all(t is not None for t in r["nr50_tokens"])

    post_files = sorted((root / f"l3/token_postings/release={rel}").rglob("*.parquet"))
    got_postings = {
        (id_to_acc[r["token_id"]], r["contig_uid"], r["feature_uid"])
        for f in post_files
        for r in pq.read_table(f, columns=["token_id", "contig_uid", "feature_uid"]).to_pylist()
        if r["token_id"] in id_to_acc
    }
    assert got_postings == want_postings


# --- (4) end-to-end refusal of a release over a stage built at `--depth stored` ---------------------------


@pytest.mark.covers("C-38a")
def test_c38a_stage_built_stored_completes_ingest_and_cannot_enter_a_release(tmp_path, capsys):
    w = c08.Work(tmp_path / "work")
    w.prepare()
    prof = c08.profile_at(tmp_path / "prof")
    pf = profile_path(prof, tmp_path / "profile.yaml")
    root = tmp_path / "root"

    # Build S0 at `--depth stored` (inputs + ingest only)
    assert cli.main([
        "build-stage", "--profile", str(pf), "--root", str(root), "--stage", "S0",
        "--depth", "stored", "--genomes", str(w.genomes("S0")), "--objects", str(w.objects("S0")),
        "--threads", "1",
    ]) == 0
    assert record_names(root) == {"stage=S0/inputs", "stage=S0/ingest"}
    assert bs.stage_complete(root, "S0", "stored")
    assert not bs.stage_complete(root, "S0", "called")
    assert catalog.stage_depth_on_disk(root, "S0") == "stored"
    capsys.readouterr()

    # CLI `bio-box release define` refuses S0, naming S0 and `stored`
    assert cli.main(["release", "define", "--root", str(root), "--release", "rBad.1", "--stages", "S0"]) == 2
    err = capsys.readouterr().err
    assert "S0" in err and "stored" in err
    assert not catalog.release_def_path(root, "rBad.1").exists()

    # `build-stage` building S1 at `called` with `--release rBad.2 --stages S0,S1` runs S1's stage steps
    # and refuses before clustering, naming S0 and `stored`, publishing no release file
    assert cli.main([
        "build-stage", "--profile", str(pf), "--root", str(root), "--stage", "S1",
        "--depth", "called", "--release", "rBad.2", "--stages", "S0,S1",
        "--genomes", str(w.genomes("S1")), "--objects", str(w.objects("S1")), "--threads", "1",
    ]) == 2
    err2 = capsys.readouterr().err
    assert "S0" in err2 and "stored" in err2
    assert bs.stage_complete(root, "S1", "called")
    assert not any("release=rBad.2" in p.as_posix() for p in root.rglob("*") if p.is_file())
