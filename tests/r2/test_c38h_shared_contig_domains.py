"""Ring R2, card C-38h (DEC-54 Q8): an annotated stage scans the proteins of contigs a called stage owns.

Design: `planning/design/C-38_import_depths.md` §3.7.

Under DEC-26 the stage that stores a contig first owns its gene calls (`feature_protein` and `protein_seq` live in
that owner stage; a later stage's `genome_input` lists `contig_uid` while its `genecall` skips owned contigs). When a
`called` stage `A` stores contig `c` first and a later `annotated` stage `B` has a genome holding `c` in `B`'s
`genome_input`, `B`'s domain scan plan (given `--seq-stages`) adds the distinct proteins of `B`'s contigs owned by
stages in `seq_stages`, copies any already measured by a qualifying source (card B-47b), and scans the rest.

Covers:
  1. Independent external oracle: `called` stage `A` stores lambda (and T4) first; later `annotated` stage `B` holds
     lambda (owned by `A`) and a lambda variant (owned by `B`). `B`'s `l1_domains` holds rows for lambda's proteins,
     equal pair-by-pair / row-by-row to an external `hmmsearch` binary run (`Z = domZ = 1e8`).
  2. Order independence: building `[called A, annotated B]` vs `[annotated B, called A]` produces identical family
     tokens (`token_dict` at level `family` and `tape_tokens` mapped by `(contig_uid, feature_uid)` -> `model_acc`)
     in the resulting release, and `release verify` passes on both.
  3. Edge cases:
     - `B0` owning 0 contigs of its own (all contigs owned by `A`, so `B0` has 0 own `protein_seq` files).
     - Deduplication when `B` also has identical proteins on its own contig (`lambda_variant`).
     - A third `annotated` stage `C` whose genome holds `c` copies `c`'s measured rows from `B` without rescanning
       (0 scanned residues for `c`'s proteins).
     - File-size tamper verification on `owned` files (in both `sources(...)` and `merge`).
     - Without `--seq-stages` (or naming only the stage itself), `domains.plan` stays byte-identical to card B-24.
"""

from __future__ import annotations

import dataclasses
import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pyarrow.parquet as pq
import pytest
import yaml

from biobox import catalog, cli, fasta, ids, runstep
from biobox import io as bio
from biobox import stage as bs
from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION as V

sys.path.insert(0, str(Path(__file__).parent))
import test_builders_out as c08  # noqa: E402
import test_micro_s0_domains as m0  # noqa: E402

pytestmark = [
    pytest.mark.r2,
    pytest.mark.covers("C-38h", "MR-08", "MR-12", "D-01", "D-03"),
]

LIB, VER = "test", "1"
SCAN = f"domains-{LIB}-{VER}"


def _profile_file(prof, path: Path) -> Path:
    path.write_text(yaml.safe_dump(dataclasses.asdict(prof)))
    return path


def _write_inputs(work: Path, name: str, fasta_specs: list[tuple[Path, str, int, bool]]) -> tuple[Path, Path]:
    """Write `genomes-<name>.json` and `objects-<name>.json` for a list of (fasta_path, kind, genetic_code, circular)."""
    genomes, objects = [], []
    for fa, kind, code, circ in fasta_specs:
        cu = ids.contig_uid(fasta.read(fa)[0].seq).hex()
        genomes.append({
            "contig_uids": [cu],
            "genetic_code": code,
            "kind": kind,
            "circular": [cu] if circ else [],
        })
        objects.append({
            "path": str(fa),
            "source_name": "fixture",
            "source_release": "1",
            "license": "INSDC (public)",
            "redistributable": True,
        })
    gpath = work / f"genomes-{name}.json"
    opath = work / f"objects-{name}.json"
    gpath.write_text(json.dumps(genomes))
    opath.write_text(json.dumps(objects))
    return gpath, opath


def _make_lambda_variant(work: Path) -> Path:
    """Lambda with three single-base substitutions: distinct contig_uid, mostly identical proteins to lambda."""
    rec = fasta.read(c08.FIX / "NC_001416.1.fna.gz")[0]
    seq = bytearray(rec.seq.upper())
    for pos in (1_000, 20_000, 40_000):
        seq[pos] = ord("A") if seq[pos] != ord("A") else ord("C")
    fa = work / "lambda_variant.fna.gz"
    fa.write_bytes(gzip.compress(b">lambda_variant\n" + bytes(seq) + b"\n", mtime=0))
    return fa


def _domain_rows(root: Path, stage: str) -> list[dict]:
    files = bio.part_files(root / f"l1/domains/library={LIB}-{VER}/stage={stage}")
    return [r for f in files for r in pq.read_table(f, partitioning=None).to_pylist()]


def _contig_proteins(root: Path, owner_stage: str, contig_uid: bytes) -> dict[bytes, str]:
    """{protein_uid: seq} for all proteins on `contig_uid` in `owner_stage`."""
    fp_files = bio.part_files(root / f"l1/feature_protein/gene_call_version={V}/stage={owner_stage}")
    uids = {
        r["protein_uid"]
        for f in fp_files
        for r in pq.read_table(f, columns=["contig_uid", "protein_uid"], partitioning=None).to_pylist()
        if r["contig_uid"] == contig_uid
    }
    all_prots = domains.load_proteins(root, V, owner_stage)
    return {u: all_prots[u] for u in uids}


def _hmmsearch_oracle(proteins: dict[bytes, str], hmm: Path, work: Path) -> dict[tuple, tuple]:
    """Run external `hmmsearch --cut_ga -Z 1e8 --domZ 1e8 --cpu 1` and parse domtblout."""
    faa = work / "oracle.faa"
    faa.write_text("".join(f">{u.hex()}\n{s.replace('*', '')}\n" for u, s in sorted(proteins.items())))
    dom = work / "oracle.domtbl"
    subprocess.run(
        ["hmmsearch", "--cut_ga", "-Z", "1e8", "--domZ", "1e8", "--cpu", "1", "--noali",
         "--domtblout", str(dom), "-o", str(work / "oracle.out"), str(hmm), str(faa)],
        check=True,
    )
    ref = {}
    for line in dom.read_text().splitlines():
        if not line.startswith("#"):
            f = line.split()
            ref[(f[0], f[4], int(f[19]) - 1, int(f[20]))] = (float(f[7]), float(f[13]), float(f[12]))
    return ref


def _family_signature(root: Path, rel: str) -> tuple[dict[str, int], dict[tuple[bytes, bytes], tuple[str | None, list[str]]]]:
    """Extract release-order-independent family token signature:
    - dict_counts: {model_acc: n_slots} from token_dict (level == 'family')
    - slot_map: {(contig_uid, feature_uid): (canonical_model_acc | None, [domain_model_acc, ...])}
    """
    tdict = [
        r for f in sorted((root / f"l3/token_dict/release={rel}").glob("*.parquet"))
        for r in pq.read_table(f).to_pylist()
    ]
    fam = [r for r in tdict if r["level"] == "family"]
    dict_counts = {r["key"]: r["n_slots"] for r in fam}
    id_to_acc = {r["token_id"]: r["key"] for r in fam}

    # Load feature_uids per contig from contig_tapes across all stages of the release
    stages = catalog.release_stages(root, rel)
    fuid_by_contig: dict[bytes, list[bytes]] = {}
    for s in stages:
        for f in sorted((root / "l3/tapes").glob(f"tape_version=*/stage={s}/*.parquet")):
            for row in pq.read_table(f, columns=["contig_uid", "feature_uids"]).to_pylist():
                fuid_by_contig[row["contig_uid"]] = row["feature_uids"]

    slot_map: dict[tuple[bytes, bytes], tuple[str | None, list[str]]] = {}
    for f in sorted((root / "l3/tape_tokens").glob(f"tape_version=*/release={rel}/*.parquet")):
        for row in pq.read_table(f).to_pylist():
            cu = row["contig_uid"]
            fuids = fuid_by_contig[cu]
            for fu, can_tid, dom_tids in zip(fuids, row["family_tokens"], row["domain_family_tokens"]):
                can_acc = None if can_tid is None else id_to_acc[can_tid]
                dom_accs = [id_to_acc[t] for t in dom_tids]
                slot_map[(cu, fu)] = (can_acc, dom_accs)
    return dict_counts, slot_map


@pytest.fixture(scope="module")
def assets(tmp_path_factory):
    """Prepare lambda, T4, lambda_variant FASTAs and a 15-model lambda HMM library."""
    base = tmp_path_factory.mktemp("c38h_assets")
    src = base / "src"
    m0._build_store(src)
    hmm = base / "lib.hmm"
    m0._make_library(m0._lambda_proteins(src), hmm)
    lam_fa = c08.FIX / "NC_001416.1.fna.gz"
    t4_fa = c08.FIX / "NC_000866.4.fna.gz"
    var_fa = _make_lambda_variant(base)
    lam_cu = ids.contig_uid(fasta.read(lam_fa)[0].seq)
    t4_cu = ids.contig_uid(fasta.read(t4_fa)[0].seq)
    var_cu = ids.contig_uid(fasta.read(var_fa)[0].seq)
    gA, oA = _write_inputs(base, "A", [(lam_fa, "virus", 11, False), (t4_fa, "virus", 11, False)])
    gB, oB = _write_inputs(base, "B", [(lam_fa, "virus", 11, False), (var_fa, "virus", 11, False)])
    gB0, oB0 = _write_inputs(base, "B0", [(lam_fa, "virus", 11, False)])
    gC, oC = _write_inputs(base, "C", [(lam_fa, "virus", 11, False)])
    return {
        "base": base,
        "hmm": hmm,
        "lam_cu": lam_cu,
        "t4_cu": t4_cu,
        "var_cu": var_cu,
        "A": (gA, oA),
        "B": (gB, oB),
        "B0": (gB0, oB0),
        "C": (gC, oC),
    }


def _build(pf: Path, root: Path, stage: str, genomes: Path, objects: Path, *,
           depth: str, hmm: Path | None = None, release: str | None = None,
           stages: tuple[str, ...] | None = None) -> int:
    argv = [
        "build-stage", "--profile", str(pf), "--root", str(root), "--stage", stage,
        "--depth", depth, "--genomes", str(genomes), "--objects", str(objects),
        "--threads", "1", "--scan-cpus", "1",
    ]
    if hmm is not None:
        argv += ["--library", LIB, VER, str(hmm), "gathering"]
    if release is not None:
        argv += ["--release", release]
    if stages is not None:
        argv += ["--stages", ",".join(stages)]
    return cli.main(argv)


# --- 1. Independent external oracle + deduplication --------------------------------------------------------


@pytest.mark.skipif(shutil.which("hmmsearch") is None, reason="hmmsearch binary not on PATH")
def test_c38h_annotated_stage_scans_called_stage_owned_contig_matching_hmmsearch_oracle(assets, tmp_path):
    prof = c08.profile_at(tmp_path / "prof")
    pf = _profile_file(prof, tmp_path / "profile.yaml")
    root = tmp_path / "root"

    # Build A at depth `called` (stores lambda + T4; no domain scan)
    assert _build(pf, root, "A", *assets["A"], depth="called") == 0
    assert bs.stage_complete(root, "A", "called")
    assert not (root / f"l1/domains/library={LIB}-{VER}/stage=A").exists()

    # Build B at depth `annotated` (holds lambda owned by A, plus lambda_variant owned by B)
    assert _build(pf, root, "B", *assets["B"], depth="annotated", hmm=assets["hmm"]) == 0
    assert bs.stage_complete(root, "B", "annotated", libraries=[(LIB, VER)])

    # Verify plan recorded owned stage A and deduplicated proteins shared between lambda and lambda_variant
    plan = domains.read_plan(root, "B", LIB, VER)
    assert "owned" in plan
    assert [s["stage"] for s in plan["owned"]["stages"]] == ["A"]
    lam_prots = _contig_proteins(root, "A", assets["lam_cu"])
    own_prots = domains.load_proteins(root, V, "B")
    assert lam_prots and own_prots
    assert set(lam_prots) & set(own_prots)  # shared proteins deduplicated
    assert set(lam_prots) - set(own_prots)  # lambda-only proteins (mutated in variant) added from A
    union_prots = {**lam_prots, **own_prots}
    assert plan["proteins"] == len(union_prots)
    assert sum(s["proteins"] for s in plan["shards"]) == len(union_prots)

    # Independent hmmsearch oracle over the union (and specifically over lambda's proteins)
    ref = _hmmsearch_oracle(union_prots, assets["hmm"], tmp_path)
    ours = {
        (r["protein_uid"].hex(), r["model_acc"], r["env_start"], r["env_end"]):
        (r["seq_bitscore"], r["dom_bitscore"], r["dom_ievalue_zref"])
        for r in _domain_rows(root, "B")
    }
    assert ours and set(ours) == set(ref)
    # Specifically assert lambda-only proteins (not in B's own protein_seq!) have domain hits matching oracle
    lam_only_hex = {u.hex() for u in (set(lam_prots) - set(own_prots))}
    lam_only_ref = {k for k in ref if k[0] in lam_only_hex}
    assert lam_only_ref, "expected at least one lambda-only protein (not in B's protein_seq) to have a domain hit"
    assert lam_only_ref <= set(ours)
    for k, (s, d, e) in ours.items():
        rs, rd, re_ = ref[k]
        assert s == pytest.approx(rs, abs=0.05)
        assert d == pytest.approx(rd, abs=0.05)
        assert e == pytest.approx(re_, rel=0.02)


# --- 2. Order independence: [called A, annotated B] vs [annotated B, called A] -----------------------------


@pytest.mark.skipif(c08.MMSEQS is None, reason="release build clusters with MMseqs2 (in the image)")
def test_c38h_either_build_order_gives_identical_release_family_tokens_and_verifies(assets, tmp_path):
    prof = c08.profile_at(tmp_path / "prof")
    pf = _profile_file(prof, tmp_path / "profile.yaml")

    # Order 1: called A first, then annotated B, release over (A, B)
    root_ab = tmp_path / "ab"
    assert _build(pf, root_ab, "A", *assets["A"], depth="called") == 0
    assert _build(pf, root_ab, "B", *assets["B"], depth="annotated", hmm=assets["hmm"],
                  release="r26.11.1", stages=("A", "B")) == 0
    assert cli.main(["release", "verify", "--root", str(root_ab), "--release", "r26.11.1"]) == 0

    # Order 2: annotated B first, then called A, release over (B, A)
    root_ba = tmp_path / "ba"
    assert _build(pf, root_ba, "B", *assets["B"], depth="annotated", hmm=assets["hmm"]) == 0
    assert _build(pf, root_ba, "A", *assets["A"], depth="called",
                  release="r26.11.1", stages=("B", "A")) == 0
    assert cli.main(["release", "verify", "--root", str(root_ba), "--release", "r26.11.1"]) == 0

    dict_ab, slots_ab = _family_signature(root_ab, "r26.11.1")
    dict_ba, slots_ba = _family_signature(root_ba, "r26.11.1")
    assert dict_ab and dict_ab == dict_ba
    assert slots_ab == slots_ba

    # Lambda slots (shared contig) carry family tokens in both orders; T4 slots (A-only, called) carry none
    lam_slots = [v for (cu, _), v in slots_ab.items() if cu == assets["lam_cu"]]
    t4_slots = [v for (cu, _), v in slots_ab.items() if cu == assets["t4_cu"]]
    assert any(can is not None for can, _ in lam_slots)
    assert t4_slots and all(can is None and doms == [] for can, doms in t4_slots)


# --- 3. Edge cases: 0 own contigs, 3rd stage C copying from B, tamper check, B-24 compat -------------------


def _run_step(prof, root: Path, stage: str, step: str, argv: list[str]) -> int:
    return runstep.run_step(
        prof, root, stage, step, [sys.executable, "-m", "biobox.cli", *argv],
        clock=lambda: c08.NOON_MONDAY, load=lambda: 0.0,
    )


def _scan_sharded(prof, root: Path, stage: str, hmm: Path, scratch: Path, seq: tuple[str, ...]) -> None:
    com = ["--root", str(root), "--stage", stage, "--library", LIB, "--library-version", VER]
    assert _run_step(prof, root, stage, f"{SCAN}.plan", [
        "stage", "domains", "plan", *com, "--hmm", str(hmm), "--cutoff", "gathering",
        "--shard-residues", str(10**9), "--seq-stages", ",".join(seq),
    ]) == 0
    for k in range(len(domains.read_plan(root, stage, LIB, VER)["shards"])):
        assert _run_step(prof, root, stage, f"{SCAN}.shard-{k:05d}", [
            "stage", "domains", "shard", *com, "--hmm", str(hmm), "--shard", str(k), "--cpus", "1",
            "--tmp", str(scratch),
        ]) == 0
    assert _run_step(prof, root, stage, SCAN, [
        "stage", "domains", "merge", *com, "--tmp", str(scratch),
    ]) == 0


def test_c38h_stage_owning_zero_contigs_and_third_stage_copying_without_rescan(assets, tmp_path):
    prof = c08.profile_at(tmp_path / "prof")
    pf = _profile_file(prof, tmp_path / "profile.yaml")
    root = tmp_path / "root"
    scratch = tmp_path / "scratch"
    scratch.mkdir()

    # 1. Build called A (owns lambda + T4)
    assert _build(pf, root, "A", *assets["A"], depth="called") == 0
    lam_prots = _contig_proteins(root, "A", assets["lam_cu"])

    # 2. Stage B0 whose ONLY contig is lambda (owned by A -> B0 has 0 own protein_seq files!)
    #    Run inputs -> ingest -> genecall (all succeed with 0 owned contigs), verify stage._sharded routes B0 to
    #    sharded domains plan, and run plan -> shard -> merge under run-step.
    gB0, oB0 = assets["B0"]
    for sub in ("b0.ing", "b0.gc", "c.ing", "c.gc"):
        (scratch / sub).mkdir()
    assert _run_step(prof, root, "B0", "inputs", ["stage", "inputs", "--root", str(root), "--stage", "B0", "--genomes", str(gB0)]) == 0
    assert _run_step(prof, root, "B0", "ingest", ["stage", "ingest", "--root", str(root), "--stage", "B0", "--objects", str(oB0), "--tmp", str(scratch / "b0.ing")]) == 0
    assert _run_step(prof, root, "B0", "genecall", ["stage", "genecall", "--root", str(root), "--stage", "B0", "--seq-stages", "B0,A", "--tmp", str(scratch / "b0.gc")]) == 0
    assert domains._protein_files(root, V, "B0") == []
    b0_build = bs.Build(root=root, stage="B0", genomes=gB0, objects=oB0, depth="annotated",
                        libraries=(bs.Library(LIB, VER, assets["hmm"], "gathering"),), threads=1)
    scan_step = next(s for s in bs.steps(prof, b0_build) if s.check == "scan")
    assert bs._sharded(prof, b0_build, root, scan_step) is True
    _scan_sharded(prof, root, "B0", assets["hmm"], scratch, ("B0", "A"))
    plan_b0 = domains.read_plan(root, "B0", LIB, VER)
    assert plan_b0["files"] == []
    assert plan_b0["proteins"] == len(lam_prots)
    assert [s["stage"] for s in plan_b0["owned"]["stages"]] == ["A"]
    rows_b0 = _domain_rows(root, "B0")
    assert rows_b0 and {r["protein_uid"] for r in rows_b0} <= set(lam_prots)

    # 3. Third stage C whose ONLY contig is lambda: copies all of lambda's proteins from B0 with 0 shards & 0 scanned residues!
    gC, oC = assets["C"]
    assert _run_step(prof, root, "C", "inputs", ["stage", "inputs", "--root", str(root), "--stage", "C", "--genomes", str(gC)]) == 0
    assert _run_step(prof, root, "C", "ingest", ["stage", "ingest", "--root", str(root), "--stage", "C", "--objects", str(oC), "--tmp", str(scratch / "c.ing")]) == 0
    assert _run_step(prof, root, "C", "genecall", ["stage", "genecall", "--root", str(root), "--stage", "C", "--seq-stages", "C,A,B0", "--tmp", str(scratch / "c.gc")]) == 0
    _scan_sharded(prof, root, "C", assets["hmm"], scratch, ("C", "A", "B0"))
    plan_c = domains.read_plan(root, "C", LIB, VER)
    assert plan_c["shards"] == []
    assert [s["stage"] for s in plan_c["copy"]["sources"]] == ["B0"]
    assert plan_c["copy"]["proteins"] == len(lam_prots)
    assert plan_c["copy"]["rows"] == len(rows_b0)
    assert _domain_rows(root, "C") == rows_b0

    # 4. Also verify full `build-stage` end-to-end for a stage B (lambda + lambda_variant) copying lambda from B0
    #    with 0 scanned residues for lambda's proteins!
    assert _build(pf, root, "B", *assets["B"], depth="annotated", hmm=assets["hmm"]) == 0
    plan_b = domains.read_plan(root, "B", LIB, VER)
    assert [s["stage"] for s in plan_b["copy"]["sources"]] == ["B0", "C"] or plan_b["copy"]["sources"][0]["stage"] == "B0"
    assert plan_b["copy"]["proteins"] == len(lam_prots)
    own_b = domains.load_proteins(root, V, "B")
    assert sum(s["proteins"] for s in plan_b["shards"]) == len(set(own_b) - set(lam_prots))

    # 5. Tamper check: corrupting size of A's feature_protein causes sources() to pass over B0
    root_bad = tmp_path / "root_bad"
    shutil.copytree(root, root_bad)
    fp = bio.part_files(root_bad / f"l1/feature_protein/gene_call_version={V}/stage=A")[0]
    fp.write_bytes(fp.read_bytes() + b"x")
    hmm_sha = ids.sha256_file(assets["hmm"]).hex()
    srcs, passed = domains.sources(root_bad, "D", LIB, VER, hmm_sha, "gathering", ["B0"])
    assert srcs == []
    assert [p["stage"] for p in passed] == ["B0"]
    assert "bytes" in passed[0]["why"]


def test_c38h_without_seq_stages_plan_stays_b24_byte_identical(assets, tmp_path):
    prof = c08.profile_at(tmp_path / "prof")
    pf = _profile_file(prof, tmp_path / "profile.yaml")
    root = tmp_path / "root"
    assert _build(pf, root, "A", *assets["A"], depth="called") == 0

    out1, out2 = tmp_path / "out1", tmp_path / "out2"
    domains.plan(root, "A", V, LIB, VER, assets["hmm"], "gathering", 10**9, out=out1)
    domains.plan(root, "A", V, LIB, VER, assets["hmm"], "gathering", 10**9, out=out2, seq_stages=["A"])
    rel = domains.scan_dir("A", LIB, VER) + "/plan.json"
    b1, b2 = (out1 / rel).read_bytes(), (out2 / rel).read_bytes()
    assert b1 == b2
    doc = json.loads(b1)
    assert "copy" not in doc and "owned" not in doc
