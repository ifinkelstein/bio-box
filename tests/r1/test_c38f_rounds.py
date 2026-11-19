"""Ring R1, card C-38f: annotation rounds (planning/design/C-38_import_depths.md §3.3, §3.6, §3.8).

Covers:
  - §3.8 test (1): independent `hmmsearch` binary oracle on a planted mixed (`called` + `annotated`) release
    with a hit, neighbours with and without domains, a neighbour in an annotated stage, and a protein shared
    between both stages.
  - §3.8 test (2): content-derived round id `r-<16 hex>`, idempotence (asking twice is a no-op), two overlapping
    rounds built in either order give byte-identical `round_domains` and `round_coverage` while the second scans
    only proteins the first did not request and copies the rest (verified from plan records).
  - §3.6: a round copies proteins measured by a stage built after R; `domains.sources` includes completed rounds
    so a later annotated stage copies round-measured proteins without rescanning.
  - §3.3: round ledger under `ledger/round=<id>/step=domains-<L>/` with `n_in`, `n_kept`, `no_hit`,
    `settled_elsewhere`, L-01 balanced, `provenance.check_provenance` and `catalog.verify_manifest` clean;
    `bio-box annotate round` with `--profile` runs shards via `run-step`; `bio-box search hmm --annotate` runs
    the round when root is writable and prints the exact `bio-box annotate round` command when root is read-only.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pyhmmer
import pytest

from biobox import catalog, cli, dim, ids, io, ledger, profile, provenance, runstep
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l1 import domains, genecall, igr
from biobox.l1.genecall import GENE_CALL_VERSION as GCV
from biobox.l2 import cluster
from biobox.l3 import postings, tapes, tokens
from biobox.registry.schemas import TABLES
from biobox.search import pack

pytestmark = [pytest.mark.covers("C-38f")]

REPO = Path(__file__).resolve().parents[2]
MICRO = REPO / "tests/fixtures/micro_s0"
LIB, VER = "testfam", "1"
LIB_FULL = f"{LIB}-{VER}"
STEP = f"domains-{LIB_FULL}"


def _build_test_hmm_ordered(picked: list[tuple[bytes, str]], out: Path) -> Path:
    """Build a small HMMER3 library from planted proteins in the given order with GA cutoffs."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet)
    bg = pyhmmer.plan7.Background(alphabet)
    hmms = []
    for idx, (uid, seq) in enumerate(picked):
        clean = seq.replace("*", "")
        msa = pyhmmer.easel.TextMSA(
            name=f"fam_{idx}".encode(),
            sequences=[pyhmmer.easel.TextSequence(name=uid.hex().encode(), sequence=clean)],
        ).digitize(alphabet)
        hmm, _, _ = builder.build_msa(msa, bg)
        hmm.accession = f"TF{idx:05d}.1".encode()
        hmm.name = f"fam_{idx}".encode()
        hmm.cutoffs.gathering = (25.0, 25.0)
        hmms.append(hmm)
    with open(out, "wb") as fh:
        for hmm in hmms:
            hmm.write(fh)
    return out


def _hmmsearch_oracle(hmm_path: Path, proteins: dict[bytes, str], tmp: Path) -> set[tuple[bytes, str, int, int]]:
    """Independent oracle: run the system `hmmsearch` binary directly with `-Z 1e8 --domZ 1e8 --cut_ga`."""
    exe = shutil.which("hmmsearch")
    if exe is None:
        pytest.skip("hmmsearch binary not on PATH")
    faa = tmp / "oracle_in.faa"
    domtbl = tmp / "oracle.domtbl"
    with open(faa, "w", encoding="utf-8") as fh:
        for uid in sorted(proteins):
            fh.write(f">{uid.hex()}\n{proteins[uid].replace('*', '')}\n")
    subprocess.run(
        [exe, "--cpu", "1", "-Z", "1e8", "--domZ", "1e8", "--cut_ga", "--domtblout", str(domtbl),
         str(hmm_path), str(faa)],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    hits: set[tuple[bytes, str, int, int]] = set()
    for line in domtbl.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split()
        puid = bytes.fromhex(cols[0])
        acc = cols[4] if cols[4] != "-" else cols[3]
        env_from = int(cols[19]) - 1
        env_to = int(cols[20])
        hits.add((puid, acc, env_from, env_to))
    return hits


def _make_profile(tmp: Path) -> Path:
    p = tmp / "profile.yaml"
    p.write_text(
        f"max_threads_day: 40\nmax_threads_night: 56\nmem_scope: 4G\n"
        f"tmp: {tmp / 'prof_tmp'}\ndurable: {tmp / 'durable'}\nbulk: {tmp / 'bulk'}\n"
        f"nice_monolithic: 0\nnice_fanout: 0\nday_start_hour: 8\nday_end_hour: 20\n"
        f"cpu_min_per_process: 300\nsmall_file_floor_bytes: 1000\nrequire_mem_cap: false\n",
        encoding="utf-8",
    )
    return p


def _make_mixed_fixture(tmp_path: Path) -> dict:
    """Create a mixed release `r1` with:
      - S0 (annotated with `testfam-1`): NC_001416.1 (lambda)
      - S1 (called only, NOT annotated): NC_000866.4 (T4) + NC_001416.1 (shared contig/proteins with S0!)
    """
    root = tmp_path / "root"
    lam = MICRO / "NC_001416.1.fna.gz"
    t4 = MICRO / "NC_000866.4.fna.gz"

    # Stage S0: lambda
    s0_objs = [SourceObject(lam, "micro", "1", "public", True, 100)]
    ingest(s0_objs, root, "S0", shard_bases=200_000)
    ci_s0 = pq.read_table(root / "l0seq/contig_index/stage=S0/part-0.parquet").to_pylist()
    s0_cuids = [r["contig_uid"] for r in ci_s0]
    write_genome_input(root, "S0", [
        GenomeInput(contig_uids=tuple(s0_cuids), kind="virus"),
    ])
    genecall.run(root, "S0", version=GCV, prefix_width=0)
    igr.run(root, "S0", gene_call_version=GCV)
    tapes.run(root, "S0", gene_call_version=GCV)

    # Stage S1: T4 + lambda (lambda already owned by S0 -> shared contig; T4 owned by S1)
    s1_objs = [
        SourceObject(t4, "micro", "1", "public", True, 100),
        SourceObject(lam, "micro", "1", "public", True, 100),
    ]
    ingest(s1_objs, root, "S1", shard_bases=200_000)
    ci_s1 = pq.read_table(root / "l0seq/contig_index/stage=S1/part-0.parquet").to_pylist()
    s1_own_cuids = [r["contig_uid"] for r in ci_s1]
    write_genome_input(root, "S1", [
        GenomeInput(contig_uids=tuple(s1_own_cuids + s0_cuids), kind="virus"),
    ])
    genecall.run(root, "S1", version=GCV, prefix_width=0, seq_stages=["S1", "S0"])
    igr.run(root, "S1", gene_call_version=GCV)
    tapes.run(root, "S1", gene_call_version=GCV)

    # Build an HMM from 1 protein of S0 (`fam_0`) and 1 protein of S1 (`fam_1`) in fixed order
    s0_prots = domains.load_proteins(root, GCV, "S0")
    s1_prots = domains.load_proteins(root, GCV, "S1")
    s0_seed = sorted(s0_prots)[0]
    s1_seed = sorted(set(s1_prots) - set(s0_prots))[0]
    hmm = _build_test_hmm_ordered([(s0_seed, s0_prots[s0_seed]), (s1_seed, s1_prots[s1_seed])], tmp_path / "testfam.hmm")

    # Annotate S0 only (via sharded plan+shard+merge + run record so S0 is also a valid copy source)
    prof_path = _make_profile(tmp_path)
    prof = profile.load(prof_path)
    scan_tmp = tmp_path / "scan_tmp"
    scan_tmp.mkdir()
    assert runstep.run_step(
        prof, root, "S0", f"{STEP}.plan",
        ["python", "-m", "biobox.cli", "stage", "domains", "plan",
         "--root", str(root), "--stage", "S0", "--gene-call-version", GCV,
         "--library", LIB, "--library-version", VER, "--hmm", str(hmm),
         "--cutoff", "gathering", "--shard-residues", "10000000"],
        threads=1,
    ) == 0
    assert runstep.run_step(
        prof, root, "S0", f"{STEP}.shard-00000",
        ["python", "-m", "biobox.cli", "stage", "domains", "shard",
         "--root", str(root), "--stage", "S0", "--gene-call-version", GCV,
         "--library", LIB, "--library-version", VER, "--hmm", str(hmm),
         "--shard", "0", "--cpus", "1", "--tmp", str(scan_tmp)],
        threads=1,
    ) == 0
    assert runstep.run_step(
        prof, root, "S0", STEP,
        ["python", "-m", "biobox.cli", "stage", "domains", "merge",
         "--root", str(root), "--stage", "S0", "--gene-call-version", GCV,
         "--library", LIB, "--library-version", VER, "--tmp", str(scan_tmp)],
        threads=1,
    ) == 0

    # Cut mixed release r1 over S0 (annotated) and S1 (called) without requiring mmseqs in the base fixture
    catalog.define_release(root, "r1", ["S0", "S1"])
    dim.run(root, "r1", stages=["S0", "S1"])
    man = catalog.build_manifest(root, "r1", ["S0", "S1"])
    catalog.write_manifest(man)
    catalog.generate_catalog(man)

    return {
        "root": root,
        "hmm": hmm,
        "prof_path": prof_path,
        "s0_prots": s0_prots,
        "s1_prots": s1_prots,
        "s1_seed": s1_seed,
        "tmp": tmp_path,
    }


@pytest.fixture(scope="module")
def mixed_base(tmp_path_factory):
    """Build the mixed S0 (annotated) + S1 (called) fixture once per module; each test copies it."""
    return _make_mixed_fixture(tmp_path_factory.mktemp("c38f_base"))


def _clone_fixture(base: dict, dest: Path) -> dict:
    root = dest / "root"
    shutil.copytree(base["root"], root)
    prof_path = _make_profile(dest)
    return {**base, "root": root, "prof_path": prof_path, "tmp": dest}


def test_c38f_round_id_content_derived_and_order_independent():
    u1 = ids.protein_uid("MAAAAAA")
    u2 = ids.protein_uid("MCCCCCC")
    u3 = ids.protein_uid("MGGGGGG")
    hmm_sha = hashlib.sha256(b"dummy-hmm").hexdigest()

    rid1 = domains.round_id(LIB, VER, hmm_sha, "gathering", [u2, u1, u3])
    rid2 = domains.round_id(LIB, VER, hmm_sha, "gathering", [u3, u1, u2, u1])
    assert rid1 == rid2
    assert rid1.startswith("r-") and len(rid1) == 18
    int(rid1[2:], 16)  # valid 16 hex digits

    # Different parameter or uid set -> different id
    assert domains.round_id(LIB, VER, hmm_sha, "trusted", [u1, u2, u3]) != rid1
    assert domains.round_id(LIB, VER, hmm_sha, "gathering", [u1, u2]) != rid1
    assert domains.round_id(LIB, "2", hmm_sha, "gathering", [u1, u2, u3]) != rid1


def test_c38f_oracle_mixed_release_matches_hmmsearch_binary(mixed_base, tmp_path):
    """§3.8 test (1): filter to `not scanned` in R, match independent `hmmsearch` binary oracle pair by pair,
    verify round_coverage (`n_rows = 0` for non-hits), ledger L-01 balance, provenance, and release verify."""
    fx = _clone_fixture(mixed_base, tmp_path)
    root, hmm, s0_prots, s1_prots = fx["root"], fx["hmm"], fx["s0_prots"], fx["s1_prots"]

    # Request ALL proteins from both S0 (annotated) and S1 (called)
    all_uids = sorted(set(s0_prots) | set(s1_prots))
    req_file = tmp_path / "proteins.txt"
    req_file.write_text("\n".join(u.hex() for u in all_uids) + "\n", encoding="utf-8")

    # Expected unscanned proteins: those in S1 that are NOT in S0
    unscanned = {u: s1_prots[u] for u in sorted(set(s1_prots) - set(s0_prots))}
    assert len(unscanned) > 0

    code = cli.main([
        "annotate", "round",
        "--root", str(root), "--release", "r1",
        "--library", LIB, VER, str(hmm), "gathering",
        "--proteins", str(req_file),
        "--profile", str(fx["prof_path"]),
        "--scan-jobs", "1", "--scan-cpus", "1",
    ])
    assert code == 0

    hmm_sha = ids.sha256_file(hmm).hex()
    rid = domains.round_id(LIB, VER, hmm_sha, "gathering", sorted(unscanned))
    assert rid in catalog.complete_rounds(root)

    # Check round_domains against independent hmmsearch binary oracle (when hmmsearch is on PATH; AGENTS §5.5)
    dom_files = io.part_files(root / f"l1/round_domains/round={rid}/library={LIB_FULL}")
    assert dom_files
    dom_rows = [r for f in dom_files for r in pq.read_table(f).to_pylist()]
    got_pairs = {(r["protein_uid"], r["model_acc"], r["env_start"], r["env_end"]) for r in dom_rows}
    if shutil.which("hmmsearch") is not None:
        oracle_pairs = _hmmsearch_oracle(hmm, unscanned, tmp_path)
        assert got_pairs == oracle_pairs and len(got_pairs) > 0
    else:
        assert len(got_pairs) > 0

    # Check round_coverage: 1 row per unscanned protein, n_rows == 0 for non-hits, > 0 for hits
    cov_files = io.part_files(root / f"l1/round_coverage/round={rid}/library={LIB_FULL}")
    assert cov_files
    cov_rows = [r for f in cov_files for r in pq.read_table(f).to_pylist()]
    cov_map = {r["protein_uid"]: r["n_rows"] for r in cov_rows}
    assert set(cov_map) == set(unscanned)
    hit_counts: dict[bytes, int] = {}
    for r in dom_rows:
        hit_counts[r["protein_uid"]] = hit_counts.get(r["protein_uid"], 0) + 1
    assert all(cov_map[u] == hit_counts.get(u, 0) for u in unscanned)
    assert any(v == 0 for v in cov_map.values())
    assert any(v > 0 for v in cov_map.values())

    # Check round ledger: L-01 balanced and valid provenance
    led_path = root / f"ledger/round={rid}/step={STEP}/part-0.parquet"
    assert led_path.is_file()
    led_rows = pq.read_table(led_path).to_pylist()
    by_reason = {r["reason_code"]: r for r in led_rows}
    assert by_reason["ok"]["n_in"] == len(unscanned)
    assert by_reason["ok"]["n_kept"] == len(hit_counts)
    assert by_reason["no_hit"]["n_rejected"] == len(unscanned) - len(hit_counts)
    assert sum(r["n_rejected"] for r in led_rows) + by_reason["ok"]["n_kept"] == by_reason["ok"]["n_in"]

    # Earlier release r1 still verifies unchanged (ignores rounds it does not list)
    m1 = catalog.read_manifest(catalog.manifest_path(root, "r1"))
    assert catalog.verify_manifest(m1, root) == []

    # A new release r2 lists the round and passes verify_manifest + check_provenance
    catalog.define_release(root, "r2", ["S0", "S1"])
    assert catalog.release_rounds(root, "r2") == [rid]
    dim.run(root, "r2", stages=["S0", "S1"])
    m2 = catalog.build_manifest(root, "r2", ["S0", "S1"])
    catalog.write_manifest(m2)
    assert catalog.verify_manifest(m2, root) == []
    assert provenance.check_provenance(m2) == []


def test_c38f_idempotence_overlap_order_and_copy_parity(mixed_base, tmp_path):
    """§3.8 test (2): asking for the same round twice is a no-op; two overlapping rounds built in either order
    produce byte-identical `round_domains` and `round_coverage`, while the second scans only the new proteins
    and copies the overlap from the first."""
    fx = mixed_base
    unscanned = sorted(set(fx["s1_prots"]) - set(fx["s0_prots"]))
    assert len(unscanned) >= 6
    # Two overlapping sets A and B
    mid = len(unscanned) // 2
    set_a = unscanned[: mid + 2]
    set_b = unscanned[mid - 2 :]
    overlap = set(set_a) & set(set_b)
    assert len(overlap) >= 4

    hmm = fx["hmm"]
    hmm_sha = ids.sha256_file(hmm).hex()
    rid_a = domains.round_id(LIB, VER, hmm_sha, "gathering", set_a)
    rid_b = domains.round_id(LIB, VER, hmm_sha, "gathering", set_b)

    def run_order(tag: str, first_uids: list[bytes], second_uids: list[bytes]) -> dict:
        work = tmp_path / tag
        root = work / "root"
        shutil.copytree(fx["root"], root)
        prof_path = _make_profile(work)
        f1 = work / "first.txt"
        f2 = work / "second.txt"
        f1.write_text("\n".join(u.hex() for u in first_uids) + "\n", encoding="utf-8")
        f2.write_text("\n".join(u.hex() for u in second_uids) + "\n", encoding="utf-8")

        for f in (f1, f1, f2):  # f1 asked twice -> second is a no-op
            assert cli.main([
                "annotate", "round",
                "--root", str(root), "--release", "r1",
                "--library", LIB, VER, str(hmm), "gathering",
                "--proteins", str(f),
                "--profile", str(prof_path),
                "--scan-jobs", "1", "--scan-cpus", "1",
            ]) == 0

        out = {}
        for rid in (rid_a, rid_b):
            out[rid] = {
                "domains": [p.read_bytes() for p in io.part_files(root / f"l1/round_domains/round={rid}/library={LIB_FULL}")],
                "coverage": [p.read_bytes() for p in io.part_files(root / f"l1/round_coverage/round={rid}/library={LIB_FULL}")],
                "plan": domains.read_round_plan(root, rid, LIB, VER),
                "ledger": pq.read_table(root / f"ledger/round={rid}/step={STEP}/part-0.parquet").to_pylist(),
            }
        return out

    ab = run_order("order_ab", set_a, set_b)
    ba = run_order("order_ba", set_b, set_a)

    # Byte-identical round_domains and round_coverage in either order
    for rid in (rid_a, rid_b):
        assert ab[rid]["domains"] == ba[rid]["domains"]
        assert ab[rid]["coverage"] == ba[rid]["coverage"]

    # Second round scans only new proteins and copies the overlap
    assert ab[rid_b]["plan"]["copy"]["proteins"] == len(overlap)
    assert sum(s["proteins"] for s in ab[rid_b]["plan"]["shards"]) == len(set_b) - len(overlap)
    assert ba[rid_a]["plan"]["copy"]["proteins"] == len(overlap)
    assert sum(s["proteins"] for s in ba[rid_a]["plan"]["shards"]) == len(set_a) - len(overlap)

    # Ledger records settled_elsewhere for copied proteins and balances L-01
    led_b = {r["reason_code"]: r for r in ab[rid_b]["ledger"]}
    assert led_b["settled_elsewhere"]["n_rejected"] == len(overlap)
    assert led_b["ok"]["n_kept"] + sum(r["n_rejected"] for r in led_b.values()) == led_b["ok"]["n_in"]


def test_c38f_round_copies_later_stage_and_stage_copies_completed_round(mixed_base, tmp_path):
    """§3.6: (1) a round copies proteins measured by a stage built after R; (2) `domains.sources` discovers
    completed rounds so a later annotated stage copies round-measured proteins without rescanning."""
    fx = _clone_fixture(mixed_base, tmp_path)
    root, hmm = fx["root"], fx["hmm"]
    hmm_sha = ids.sha256_file(hmm).hex()
    unscanned = sorted(set(fx["s1_prots"]) - set(fx["s0_prots"]))

    # Scan S1 in a separate copy WITHOUT any round (so S1 scans its 265 proteins + copies S0's shared-contig
    # proteins via C-38h) as our byte-identity reference
    ref_root = tmp_path / "ref_s1" / "root"
    shutil.copytree(root, ref_root)
    ref_tmp = tmp_path / "ref_s1" / "tmp"
    ref_tmp.mkdir()
    ref_plan = domains.plan(ref_root, "S1", GCV, LIB, VER, hmm, "gathering", 10_000_000, seq_stages=["S1", "S0"])
    assert ref_plan["shards"] == 1
    domains.scan_shard(ref_root, "S1", LIB, VER, hmm, 0, cpus=1, tmp=ref_tmp)
    domains.merge(ref_root, "S1", LIB, VER, tmp=ref_tmp)

    # Run a round over all unscanned proteins of S1
    req_file = tmp_path / "req.txt"
    req_file.write_text("\n".join(u.hex() for u in unscanned) + "\n", encoding="utf-8")
    assert cli.main([
        "annotate", "round",
        "--root", str(root), "--release", "r1",
        "--library", LIB, VER, str(hmm), "gathering",
        "--proteins", str(req_file),
        "--profile", str(fx["prof_path"]),
        "--scan-jobs", "1", "--scan-cpus", "1",
    ]) == 0
    rid = domains.round_id(LIB, VER, hmm_sha, "gathering", unscanned)

    # `domains.sources` discovers the completed round as a copy source for S1
    srcs, _ = domains.sources(root, "S1", LIB, VER, hmm_sha, "gathering", ["S0"])
    src_ids = [s.get("round") or s.get("stage") for s in srcs]
    assert "S0" in src_ids and rid in src_ids

    # Planning S1's domain scan copies 100% of S1's proteins (from S0 + round) and scans 0 shards!
    plan_res = domains.plan(root, "S1", GCV, LIB, VER, hmm, "gathering", 10_000_000, seq_stages=["S1", "S0"])
    assert plan_res["shards"] == 0
    assert plan_res["copied_proteins"] == plan_res["proteins"] > 0
    merge_tmp = tmp_path / "merge_tmp"
    merge_tmp.mkdir()
    domains.merge(root, "S1", LIB, VER, tmp=merge_tmp)

    # Compare S1's copied l1_domains (0 shards scanned) against the reference scan (1 shard scanned)
    got = [p.read_bytes() for p in io.part_files(root / f"l1/domains/library={LIB_FULL}/stage=S1")]
    want = [p.read_bytes() for p in io.part_files(ref_root / f"l1/domains/library={LIB_FULL}/stage=S1")]
    assert got == want


def test_c38f_search_hmm_annotate_writable_and_readonly(mixed_base, tmp_path):
    """§3.3: `search hmm --annotate` runs the round over unscanned window proteins when root is writable,
    and prints the exact `bio-box annotate round` command without touching root when root is read-only."""
    if shutil.which("mmseqs") is None or shutil.which("hmmsearch") is None:
        pytest.skip("mmseqs and hmmsearch binaries required for end-to-end search hmm --annotate")
    fx = _clone_fixture(mixed_base, tmp_path)
    root, hmm = fx["root"], fx["hmm"]
    # Build clustered release r1s (search hmm needs cluster/tokens/postings + pack)
    catalog.define_release(root, "r1s", ["S0", "S1"])
    cluster.run(root, "r1s", stages=["S0", "S1"], gene_call_version=GCV, threads=1)
    tokens.run(root, "r1s", stages=["S0", "S1"], gene_call_version=GCV)
    postings.run(root, "r1s", stages=["S0", "S1"], gene_call_version=GCV)
    dim.run(root, "r1s", stages=["S0", "S1"])
    man = catalog.build_manifest(root, "r1s", ["S0", "S1"])
    catalog.write_manifest(man)
    catalog.generate_catalog(man)

    packs = tmp_path / "packs"
    runs = tmp_path / "runs"
    work = tmp_path / "work"
    pack.build(root, "r1s", GCV, {}, packs=packs)

    # Pick the planted positive control gene from S1 (called stage) matching fam_1
    s1_target_uid = fx["s1_seed"]
    occ = pq.read_table(root / f"l1/protein_occurrence/gene_call_version={GCV}/stage=S1/part-0.parquet").to_pylist()
    hit_occ = next(r for r in occ if r["protein_uid"] == s1_target_uid)
    lab = pq.read_table(root / "l0seq/contig_label/stage=S1/part-0.parquet").to_pylist()
    acc = next(r["accession"] for r in lab if r["contig_uid"] == hit_occ["contig_uid"])
    pos = f"{acc}:{hit_occ['start']}-{hit_occ['end']}"

    # 1. Read-only root: writes nothing to root, prints exact `bio-box annotate round` command
    ro_root = tmp_path / "ro_root"
    shutil.copytree(root, ro_root)
    ro_out = tmp_path / "search_ro"
    os.chmod(ro_root, 0o555)
    try:
        res = subprocess.run(
            ["python", "-m", "biobox.cli", "search", "hmm",
             "--root", str(ro_root), "--release", "r1s",
             "--hmm", str(hmm), "--model", "fam_1", "--bits", "20",
             "--neighbors", "3", "--positive", pos,
             "--pack-dir", str(packs), "--runs", str(runs), "--work", str(work),
             "--out", str(ro_out), "--cpus", "1",
             "--annotate", LIB, VER, str(hmm), "gathering"],
            capture_output=True, text=True, check=False,
        )
        assert res.returncode == 0, res.stderr
        assert "bio-box annotate round" in res.stdout or "bio-box annotate round" in res.stderr
        assert (ro_out / "round_proteins.txt").is_file()
        assert catalog.complete_rounds(ro_root) == []
    finally:
        os.chmod(ro_root, 0o755)

    # 2. Writable root: runs the round over unscanned window proteins
    rw_out = tmp_path / "search_rw"
    res2 = subprocess.run(
        ["python", "-m", "biobox.cli", "search", "hmm",
         "--root", str(root), "--release", "r1s",
         "--hmm", str(hmm), "--model", "fam_1", "--bits", "20",
         "--neighbors", "3", "--positive", pos,
         "--pack-dir", str(packs), "--runs", str(runs), "--work", str(work),
         "--out", str(rw_out), "--cpus", "1",
         "--annotate", LIB, VER, str(hmm), "gathering"],
        capture_output=True, text=True, check=False,
    )
    assert res2.returncode == 0, res2.stderr
    assert len(catalog.complete_rounds(root)) == 1


def test_c38f_round_copies_stage_built_after_release_and_multishard_parity(mixed_base, tmp_path):
    """§3.6 + multi-shard parity: (1) a round against release r1 copies proteins measured by a stage S1
    annotated AFTER r1 was cut (0 shards scanned); (2) multi-shard `annotate round --profile` produces
    byte-identical `round_domains` and `round_coverage` to single-shard."""
    fx = mixed_base
    unscanned = sorted(set(fx["s1_prots"]) - set(fx["s0_prots"]))
    hmm = fx["hmm"]
    hmm_sha = ids.sha256_file(hmm).hex()
    rid = domains.round_id(LIB, VER, hmm_sha, "gathering", unscanned)

    # Part A: multi-shard vs single-shard parity
    root_1 = tmp_path / "one_shard" / "root"
    root_n = tmp_path / "multi_shard" / "root"
    shutil.copytree(fx["root"], root_1)
    shutil.copytree(fx["root"], root_n)
    req = tmp_path / "req.txt"
    req.write_text("\n".join(u.hex() for u in unscanned) + "\n", encoding="utf-8")

    for rdir, limit in ((root_1, 10_000_000), (root_n, 15_000)):
        prof_path = _make_profile(rdir.parent)
        assert cli.main([
            "annotate", "round",
            "--root", str(rdir), "--release", "r1",
            "--library", LIB, VER, str(hmm), "gathering",
            "--proteins", str(req),
            "--profile", str(prof_path),
            "--shard-residues", str(limit),
            "--scan-jobs", "1", "--scan-cpus", "1",
        ]) == 0

    p_multi = domains.read_round_plan(root_n, rid, LIB, VER)
    assert len(p_multi["shards"]) >= 2
    assert ([p.read_bytes() for p in io.part_files(root_1 / f"l1/round_domains/round={rid}/library={LIB_FULL}")]
            == [p.read_bytes() for p in io.part_files(root_n / f"l1/round_domains/round={rid}/library={LIB_FULL}")])
    assert ([p.read_bytes() for p in io.part_files(root_1 / f"l1/round_coverage/round={rid}/library={LIB_FULL}")]
            == [p.read_bytes() for p in io.part_files(root_n / f"l1/round_coverage/round={rid}/library={LIB_FULL}")])

    # Part B: annotate S1 AFTER r1 was cut, then run a round against r1 -> copies 100% from S1, scans 0 shards
    root_after = tmp_path / "after_r1" / "root"
    shutil.copytree(fx["root"], root_after)
    prof_after = profile.load(_make_profile(root_after.parent))
    scan_tmp = root_after.parent / "scan_tmp"
    scan_tmp.mkdir()
    assert runstep.run_step(
        prof_after, root_after, "S1", f"{STEP}.plan",
        ["python", "-m", "biobox.cli", "stage", "domains", "plan",
         "--root", str(root_after), "--stage", "S1", "--gene-call-version", GCV,
         "--library", LIB, "--library-version", VER, "--hmm", str(hmm),
         "--cutoff", "gathering", "--shard-residues", "10000000"],
        threads=1,
    ) == 0
    assert runstep.run_step(
        prof_after, root_after, "S1", f"{STEP}.shard-00000",
        ["python", "-m", "biobox.cli", "stage", "domains", "shard",
         "--root", str(root_after), "--stage", "S1", "--gene-call-version", GCV,
         "--library", LIB, "--library-version", VER, "--hmm", str(hmm),
         "--shard", "0", "--cpus", "1", "--tmp", str(scan_tmp)],
        threads=1,
    ) == 0
    assert runstep.run_step(
        prof_after, root_after, "S1", STEP,
        ["python", "-m", "biobox.cli", "stage", "domains", "merge",
         "--root", str(root_after), "--stage", "S1", "--gene-call-version", GCV,
         "--library", LIB, "--library-version", VER, "--tmp", str(scan_tmp)],
        threads=1,
    ) == 0

    assert cli.main([
        "annotate", "round",
        "--root", str(root_after), "--release", "r1",
        "--library", LIB, VER, str(hmm), "gathering",
        "--proteins", str(req),
        "--profile", str(root_after.parent / "profile.yaml"),
        "--scan-jobs", "1", "--scan-cpus", "1",
    ]) == 0
    p_after = domains.read_round_plan(root_after, rid, LIB, VER)
    assert len(p_after["shards"]) == 0
    assert p_after["copy"]["proteins"] == len(unscanned)
    assert ([p.read_bytes() for p in io.part_files(root_1 / f"l1/round_domains/round={rid}/library={LIB_FULL}")]
            == [p.read_bytes() for p in io.part_files(root_after / f"l1/round_domains/round={rid}/library={LIB_FULL}")])
    assert ([p.read_bytes() for p in io.part_files(root_1 / f"l1/round_coverage/round={rid}/library={LIB_FULL}")]
            == [p.read_bytes() for p in io.part_files(root_after / f"l1/round_coverage/round={rid}/library={LIB_FULL}")])
