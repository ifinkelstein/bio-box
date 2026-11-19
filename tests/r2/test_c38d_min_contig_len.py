"""Ring R2: card C-38d — `--min-contig-len` at ingest (`planning/design/C-38_import_depths.md` §5.3 (1)–(5)).

Covers:
  (1) Multi-object FASTAs with records `< L`, `== L` (boundary), and `> L` in two input orders,
      via both single-process `ingest` and batched `pass1`/`plan`/`batch`/`merge` (+ `rescan`):
      byte-identical files across orders and between single-process and batched ingest, and ledger
      `below:<L>` and `ok` rows equal to an independent brute-force oracle loop (AGENTS §5.2).
  (2) A contig `< L` held by another stage is counted under `below:<L>` (not as `n_merged`),
      independent of whether the unfiltered stage or the `--min-contig-len L` stage is built first.
  (3) `genome_input` that disagrees with `--min-contig-len L` (a contig `< L` in `genome_input`,
      a `< L` contig held by an earlier stage in `genome_input`, or a stored `>= L` contig missing
      from `genome_input`) is refused before any ingest file is written (both single-process and
      batched `build-stage`); matching `genomes.json` succeeds through `build-stage --depth called`.
  (4) Without `--min-contig-len`, `n_bases_in` and `n_bases_rejected` are null on every ledger row
      and all pre-existing columns and L0/L0-seq tables match the unfiltered store contract.
  (5) Bounded peak RSS (`<= 1.25 GiB`, *derived* from `DUCK_MEMORY_LIMIT` 1 GiB + 256 MiB batch/Python
      headroom) when ingesting the 100x synthetic store (`scripts/stress/synth_scale.py` scale) with
      short `< L` contigs interleaved in every source object, measured via `runstep.measured` /
      `runstep.step_peak` (AGENTS §4.1, §4.2).
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import random
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
import pytest
import yaml

from biobox import cli, fasta, ids, io, ledger, profile, runstep
from biobox import stage as bs
from biobox.l0seq import SeqStore, SourceObject, ingest, store
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.registry.schemas import conformance_problems, table

pytestmark = [
    pytest.mark.r2,
    pytest.mark.covers("C-38d", "D-01", "D-03", "L-01", "L-02", "LY-09", "MR-01", "R-09"),
]

REPO = Path(__file__).resolve().parents[2]
FIX = REPO / "tests" / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
NOON_MONDAY = datetime(2026, 11, 9, 12, 0)
PEAK_BOUND = (1 << 30) + (256 << 20)  # 1.25 GiB derived: 1 GiB DuckDB cap + 256 MiB shard/interpreter


def _rand_seq(rng: random.Random, n: int) -> bytes:
    return bytes(b"ACGT"[b & 3] for b in rng.randbytes(n))


def _write_fasta(path: Path, recs: list[tuple[str, bytes]], width: int = 70) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"".join(b"\n".join(fasta.format_record(h, s, width)) + b"\n" for h, s in recs))
    return path


def _tree_sha(root: Path, *, exclude_plan: bool = False) -> dict[str, str]:
    root = Path(root)
    if not root.exists():
        return {}
    out: dict[str, str] = {}
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(root).as_posix()
        if rel.startswith("runs/"):
            continue
        if exclude_plan and rel.startswith("l0seq/ingest_plan/"):
            continue
        out[rel] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out


def _profile_and_path(tmp_path: Path) -> tuple[profile.Profile, Path]:
    base = tmp_path / "prof"
    for sub in ("tmp", "durable", "bulk"):
        (base / sub).mkdir(parents=True, exist_ok=True)
    prof = profile.Profile(
        max_threads_day=4,
        max_threads_night=8,
        mem_scope="8G",
        tmp=str(base / "tmp"),
        durable=str(base / "durable"),
        bulk=str(base / "bulk"),
        nice_monolithic=0,
        nice_fanout=0,
        day_start_hour=8,
        day_end_hour=20,
        cpu_min_per_process=300,
        small_file_floor_bytes=0,
    )
    pf = tmp_path / "profile.yaml"
    pf.write_text(yaml.safe_dump(dataclasses.asdict(prof)))
    return prof, pf


def _brute_force_ingest_oracle(
    objs: list[SourceObject],
    min_contig_len: int | None,
    earlier_digests: set[bytes] = frozenset(),
) -> tuple[list[tuple], set[bytes], set[bytes]]:
    """Independent brute-force oracle for ingest ledger rows and stored sequences/contigs (AGENTS §5.2).

    Returns:
      - expected_ledger: sorted list of tuples
        (object_key, reason_code, n_in, n_kept, n_merged, n_rejected, n_bases_in, n_bases_rejected)
      - kept_digests: set of raw_digest values newly stored in this stage
      - stage_contigs: set of contig_uid values appearing in this stage's kept (>= L) records
    """
    seen_digests = set(earlier_digests)
    ordered = sorted(
        objs,
        key=lambda o: (
            hashlib.sha256(Path(o.path).read_bytes()).digest(),
            o.source_name,
            o.source_release,
            o.precedence,
            Path(o.path).name,
        ),
    )
    expected_ledger: list[tuple] = []
    new_digests: set[bytes] = set()
    stage_contigs: set[bytes] = set()
    for o in ordered:
        raw_file = Path(o.path).read_bytes()
        sha = hashlib.sha256(raw_file).digest()
        recs = fasta.read(o.path)
        n_in = len(recs)
        bases_in = sum(len(r.seq) for r in recs)
        n_rej = 0
        bases_rej = 0
        kept = 0
        merged = 0
        for r in recs:
            if min_contig_len is not None and len(r.seq) < min_contig_len:
                n_rej += 1
                bases_rej += len(r.seq)
                continue
            rd = hashlib.sha256(r.seq).digest()
            stage_contigs.add(ids.contig_uid(r.seq))
            if rd not in seen_digests:
                seen_digests.add(rd)
                new_digests.add(rd)
                kept += 1
            else:
                merged += 1
        if min_contig_len is None:
            expected_ledger.append((sha, "ok", n_in, kept, merged, 0, None, None))
        else:
            expected_ledger.append((sha, f"below:{min_contig_len}", n_in, kept, merged, n_rej, bases_in, bases_rej))
            expected_ledger.append((sha, "ok", n_in, kept, merged, 0, bases_in, 0))
    expected_ledger.sort(key=lambda r: (r[0], r[1]))
    return expected_ledger, new_digests, stage_contigs


# --- §5.3 (1): records on both sides of L, two input orders, single-process & batched + rescan -------------


@pytest.mark.covers("C-38d", "D-01", "D-03", "L-01", "L-02", "MR-01")
def test_c38d_1_both_sides_of_min_len_two_orders_single_and_batched_match_oracle(tmp_path):
    rng = random.Random(20261118)
    min_len = 1_500
    s_1 = _rand_seq(rng, 1)              # 1 bp (< L)
    s_1499 = _rand_seq(rng, 1_499)       # L - 1 bp (< L, boundary!)
    s_1500 = _rand_seq(rng, 1_500)       # L bp (== L, boundary kept!)
    s_2400 = _rand_seq(rng, 2_400)       # > L
    s_3200 = _rand_seq(rng, 3_200)       # > L
    s_900 = _rand_seq(rng, 900)          # < L, repeated across objects

    f0 = _write_fasta(tmp_path / "src/obj0.fna", [
        ("a_short", s_1),
        ("a_edge_below", s_1499),
        ("a_edge_at", s_1500),
        ("a_long", s_2400),
    ])
    f1 = _write_fasta(tmp_path / "src/obj1.fna", [
        ("b_dup_long", s_2400),          # duplicate of obj0's >= L sequence
        ("b_short_dup", s_900),
        ("b_other_long", s_3200),
        ("b_short_dup2", s_900),
    ])
    f2 = _write_fasta(tmp_path / "src/obj2_all_short.fna", [
        ("c_only_short1", s_900),
        ("c_only_short2", s_1499),
    ])
    objs_a = [SourceObject(p, "metag", "1", "CC0", True) for p in (f0, f1, f2)]
    objs_b = list(reversed(objs_a))

    root_a = tmp_path / "order_a"
    root_b = tmp_path / "order_b"
    ingest(objs_a, root_a, "S0", shard_bases=4_000, min_contig_len=min_len)
    ingest(objs_b, root_b, "S0", shard_bases=4_000, min_contig_len=min_len)

    # Order independence (MR-01)
    assert _tree_sha(root_a) == _tree_sha(root_b)

    # Batched ingest (pass1 per object -> plan -> lost scratch rescan -> batch -> merge) equals single-process byte-for-byte
    prof, _ = _profile_and_path(tmp_path)
    root_batch = tmp_path / "batched" / "root"
    scratch = tmp_path / "batched" / "scratch"
    objs_json = tmp_path / "batched" / "objects.json"
    objs_json.parent.mkdir(parents=True, exist_ok=True)
    objs_json.write_text(json.dumps([
        {"path": str(o.path), "source_name": o.source_name, "source_release": o.source_release,
         "license": o.license, "redistributable": o.redistributable}
        for o in objs_b
    ]))
    n_p1 = len(store.pass1_shards(objs_b, 1))
    assert n_p1 == 3
    for k in range(n_p1):
        assert runstep.run_step(
            prof, root_batch, "S0", f"ingest.pass1-{k:05d}",
            [sys.executable, "-m", "biobox.cli", "stage", "ingest", "pass1",
             "--root", str(root_batch), "--stage", "S0", "--objects", str(objs_json),
             "--tmp", str(scratch), "--shard", str(k), "--pass1-bytes", "1",
             "--shard-bases", "4000", "--min-contig-len", str(min_len)],
            clock=lambda: NOON_MONDAY, load=lambda: 0.0,
        ) == 0
    assert runstep.run_step(
        prof, root_batch, "S0", "ingest.plan",
        [sys.executable, "-m", "biobox.cli", "stage", "ingest", "plan",
         "--root", str(root_batch), "--stage", "S0", "--objects", str(objs_json),
         "--tmp", str(scratch), "--shard-bases", "4000", "--min-contig-len", str(min_len)],
        clock=lambda: NOON_MONDAY, load=lambda: 0.0,
    ) == 0

    # Wipe scratch runs and verify rescan reproduces exact run_sha256 and new_seq_sha256 with min_contig_len
    shutil.rmtree(store.ingest_runs_dir(scratch, root_batch, "S0"))
    for k in range(n_p1):
        assert subprocess.run(
            [sys.executable, "-m", "biobox.l0seq.store", "rescan",
             "--root", str(root_batch), "--stage", "S0", "--objects", str(objs_json),
             "--tmp", str(scratch), "--pass1", str(k)],
            check=False,
        ).returncode == 0
    assert subprocess.run(
        [sys.executable, "-m", "biobox.l0seq.store", "rescan",
         "--root", str(root_batch), "--stage", "S0", "--objects", str(objs_json),
         "--tmp", str(scratch), "--plan"],
        check=False,
    ).returncode == 0

    plan_doc = store.read_ingest_plan(root_batch, "S0")
    for j in range(len(plan_doc["batches"])):
        assert runstep.run_step(
            prof, root_batch, "S0", f"ingest.batch-{j:04d}",
            [sys.executable, "-m", "biobox.cli", "stage", "ingest", "batch",
             "--root", str(root_batch), "--stage", "S0", "--objects", str(objs_json),
             "--tmp", str(scratch), "--batch", str(j)],
            clock=lambda: NOON_MONDAY, load=lambda: 0.0,
        ) == 0
    assert runstep.run_step(
        prof, root_batch, "S0", "ingest",
        [sys.executable, "-m", "biobox.cli", "stage", "ingest", "merge",
         "--root", str(root_batch), "--stage", "S0", "--objects", str(objs_json),
         "--tmp", str(scratch)],
        clock=lambda: NOON_MONDAY, load=lambda: 0.0,
    ) == 0
    assert _tree_sha(root_batch, exclude_plan=True) == _tree_sha(root_a)

    # Brute-force oracle comparison
    want_ledger, want_digests, want_contigs = _brute_force_ingest_oracle(objs_a, min_len)
    got_rows = ledger.read(root_a, "S0", "ingest")
    assert conformance_problems(pq.read_table(ledger.path(root_a, "S0", "ingest")), table("ledger")) == []
    got_tuples = [
        (r["object_key"], r["reason_code"], r["n_in"], r["n_kept"], r["n_merged"],
         r["n_rejected"], r["n_bases_in"], r["n_bases_rejected"])
        for r in got_rows
    ]
    assert got_tuples == want_ledger
    assert ledger.check_balance(root_a, "S0", ["ingest"]) == []
    assert ledger.check_inputs(root_a, "S0", {ids.sha256_file(o.path) for o in objs_a}) == []

    # Stored sequences and occurrence/label/index tables contain ONLY >= L contigs
    st = SeqStore(root_a, ["S0"])
    try:
        for s_kept in (s_1500, s_2400, s_3200):
            assert st.get(ids.contig_uid(s_kept)) == s_kept
        for s_dropped in (s_1, s_900, s_1499):
            with pytest.raises(KeyError):
                st.get(ids.contig_uid(s_dropped))
    finally:
        st.close()

    occ_files = sorted(root_a.glob("l0seq/contig_occurrence/stage=S0/source=*/part-*.parquet"))
    occ_rows = [r for f in occ_files for r in pq.read_table(f).to_pylist()]
    assert {r["contig_uid"] for r in occ_rows} == want_contigs
    assert len(occ_rows) == 4  # a_edge_at, a_long, b_dup_long, b_other_long
    # True 0-based FASTA record ordinals preserved on kept records (in obj1: indices 0 and 2)
    f1_sha = ids.sha256_file(f1)
    assert [r["record_ordinal"] for r in occ_rows if r["l0_object_sha256"] == f1_sha] == [0, 2]

    loc_rows = pq.read_table(root_a / "l0seq/seq_locator/stage=S0/part-0.parquet").to_pylist()
    assert {r["raw_digest"] for r in loc_rows} == want_digests


# --- §5.3 (2): a contig < L held by another stage is counted under below:<L>, in both stage orders ---------


@pytest.mark.covers("C-38d", "D-01", "L-01", "R-09")
def test_c38d_2_contig_below_min_len_held_by_another_stage_counted_below_len_in_both_orders(tmp_path):
    rng = random.Random(20261119)
    min_len = 1_500
    s_short_shared = _rand_seq(rng, 800)      # < L, present in BOTH S_unfiltered and S_filtered
    s_long_shared = _rand_seq(rng, 2_200)     # >= L, present in BOTH S_unfiltered and S_filtered
    s_long_unf = _rand_seq(rng, 1_800)        # >= L, only in S_unfiltered
    s_long_flt = _rand_seq(rng, 2_600)        # >= L, only in S_filtered

    f_unf = _write_fasta(tmp_path / "src/unfiltered.fna", [
        ("u_short", s_short_shared),
        ("u_long_shared", s_long_shared),
        ("u_long_own", s_long_unf),
    ])
    f_flt = _write_fasta(tmp_path / "src/filtered.fna", [
        ("f_short_shared", s_short_shared),
        ("f_long_shared", s_long_shared),
        ("f_long_own", s_long_flt),
    ])
    obj_unf = [SourceObject(f_unf, "plsdb", "1", "CC0", True)]
    obj_flt = [SourceObject(f_flt, "spire", "1", "CC0", True)]

    # Order 1: unfiltered stage PL1 first, then filtered stage SP1 (--min-contig-len 1500)
    r1 = tmp_path / "unfiltered_first"
    ingest(obj_unf, r1, "PL1")
    ingest(obj_flt, r1, "SP1", min_contig_len=min_len)

    # Order 2: filtered stage SP1 (--min-contig-len 1500) first, then unfiltered stage PL1
    r2 = tmp_path / "filtered_first"
    ingest(obj_flt, r2, "SP1", min_contig_len=min_len)
    ingest(obj_unf, r2, "PL1")

    # Check SP1's ledger in BOTH orders:
    # - s_short_shared (< 1500) is ALWAYS counted in `below:1500` (n_rejected=1, n_bases_rejected=800), NEVER in n_merged!
    # - s_long_shared (>= 1500) follows DEC-26: merged (1) when PL1 went first, kept (1) when SP1 went first.
    total_flt_bases = len(s_short_shared) + len(s_long_shared) + len(s_long_flt)
    for root_x, exp_kept, exp_merged in ((r1, 1, 1), (r2, 2, 0)):
        rows = {(r["reason_code"]): r for r in ledger.read(root_x, "SP1", "ingest")}
        assert set(rows) == {f"below:{min_len}", "ok"}
        below = rows[f"below:{min_len}"]
        ok = rows["ok"]
        assert (below["n_in"], below["n_kept"], below["n_merged"], below["n_rejected"]) == (3, exp_kept, exp_merged, 1)
        assert (below["n_bases_in"], below["n_bases_rejected"]) == (total_flt_bases, len(s_short_shared))
        assert (ok["n_in"], ok["n_kept"], ok["n_merged"], ok["n_rejected"]) == (3, exp_kept, exp_merged, 0)
        assert (ok["n_bases_in"], ok["n_bases_rejected"]) == (total_flt_bases, 0)
        assert ledger.check_balance(root_x, "SP1", ["ingest"]) == []
        assert ledger.check_balance(root_x, "PL1", ["ingest"]) == []

        # SP1 never records s_short_shared in contig_occurrence or contig_index
        sp1_occ = pq.read_table(root_x / "l0seq/contig_occurrence/stage=SP1/source=spire/part-0.parquet").to_pylist()
        assert ids.contig_uid(s_short_shared) not in {r["contig_uid"] for r in sp1_occ}

        # PL1 stores s_short_shared in BOTH orders (either as first stage or because SP1 dropped it!)
        pl1_loc = pq.read_table(root_x / "l0seq/seq_locator/stage=PL1/part-0.parquet").to_pylist()
        assert ids.raw_digest(s_short_shared) in {r["raw_digest"] for r in pl1_loc}

        st = SeqStore(root_x, ["PL1", "SP1"])
        try:
            for seq in (s_short_shared, s_long_shared, s_long_unf, s_long_flt):
                assert st.get(ids.contig_uid(seq)) == seq
        finally:
            st.close()


# --- §5.3 (3): disagreeing genome_input refused before anything is written; agreeing build-stage succeeds ---


@pytest.mark.covers("C-38d", "D-01", "D-03")
def test_c38d_3_disagreeing_inputs_refused_before_any_ingest_file_is_written_and_agreeing_succeeds(tmp_path):
    lam = fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq[:8_000]   # 8,000 bp (>= L)
    t4 = fasta.read(FIX / "NC_000866.4.fna.gz")[0].seq[:10_000]  # 10,000 bp (>= L)
    s_short = b"ACGT" * 200                                       # 800 bp (< L)
    s_earlier_short = b"TGCA" * 250                               # 1,000 bp (< L, held by S0)
    min_len = 1_500

    f_s0 = _write_fasta(tmp_path / "src/s0.fna", [("s0_short", s_earlier_short), ("s0_lam", lam)])
    f_s1 = _write_fasta(tmp_path / "src/s1.fna", [("s1_short", s_short), ("s1_t4", t4), ("s1_lam", lam)])
    objs_s1 = [SourceObject(f_s1, "metag", "1", "CC0", True)]

    cu_lam = ids.contig_uid(lam)
    cu_t4 = ids.contig_uid(t4)
    cu_short = ids.contig_uid(s_short)
    cu_earlier_short = ids.contig_uid(s_earlier_short)
    cu_bogus = ids.contig_uid(b"CCCC" * 600)

    def _base_root(name: str) -> Path:
        r = tmp_path / name / "root"
        write_genome_input(r, "S0", [GenomeInput((cu_earlier_short, cu_lam), 11, kind="virus")])
        ingest([SourceObject(f_s0, "src0", "1", "CC0", True)], r, "S0")
        return r

    # Case (a): genome_input lists s_short (800 bp < 1500 bp) from S1's FASTA
    ra, out_a = _base_root("case_a"), tmp_path / "case_a" / "out"
    write_genome_input(ra, "S1", [GenomeInput((cu_short, cu_t4, cu_lam), 11, kind="virus")])
    with pytest.raises(ValueError, match="min-contig-len|shorter|<"):
        ingest(objs_s1, ra, "S1", min_contig_len=min_len, out=out_a)
    assert store.stage_files(ra, "S1") == [] and store.stage_files(out_a, "S1") == []

    # Case (b): genome_input lists s_earlier_short (1000 bp < 1500 bp) held by earlier unfiltered stage S0
    rb, out_b = _base_root("case_b"), tmp_path / "case_b" / "out"
    write_genome_input(rb, "S1", [GenomeInput((cu_earlier_short, cu_t4, cu_lam), 11, kind="virus")])
    with pytest.raises(ValueError, match="min-contig-len|shorter|<"):
        ingest(objs_s1, rb, "S1", min_contig_len=min_len, out=out_b)
    assert store.stage_files(rb, "S1") == [] and store.stage_files(out_b, "S1") == []

    # Case (c): genome_input omits cu_t4 (a >= 1500 bp contig kept by S1)
    rc, out_c = _base_root("case_c"), tmp_path / "case_c" / "out"
    write_genome_input(rc, "S1", [GenomeInput((cu_lam,), 11, kind="virus")])
    with pytest.raises(ValueError, match="genome_input|missing"):
        ingest(objs_s1, rc, "S1", min_contig_len=min_len, out=out_c)
    assert store.stage_files(rc, "S1") == [] and store.stage_files(out_c, "S1") == []

    # Case (d): genome_input lists a contig absent from both S1's >= L records and earlier stages
    rd, out_d = _base_root("case_d"), tmp_path / "case_d" / "out"
    write_genome_input(rd, "S1", [GenomeInput((cu_t4, cu_lam, cu_bogus), 11, kind="virus")])
    with pytest.raises(ValueError, match="genome_input"):
        ingest(objs_s1, rd, "S1", min_contig_len=min_len, out=out_d)
    assert store.stage_files(rd, "S1") == [] and store.stage_files(out_d, "S1") == []

    # Case (e): CLI `build-stage --min-contig-len` refuses disagreeing genomes.json in BOTH single-process and
    # batched ingest (forced by tiny cpu_min_per_process) without publishing any ingest/plan file
    prof, pf = _profile_and_path(tmp_path)
    bad_gen = tmp_path / "bad_genomes.json"
    bad_gen.write_text(json.dumps([{"contig_uids": [cu_short.hex(), cu_t4.hex(), cu_lam.hex()],
                                    "genetic_code": 11, "kind": "virus"}]))
    objs_p = tmp_path / "s1_objects.json"
    objs_p.write_text(json.dumps([{"path": str(f_s1), "source_name": "metag", "source_release": "1",
                                   "license": "CC0", "redistributable": True}]))
    r_cli = _base_root("cli_refuse")
    assert cli.main([
        "build-stage", "--profile", str(pf), "--root", str(r_cli), "--stage", "S1",
        "--depth", "stored", "--genomes", str(bad_gen), "--objects", str(objs_p),
        "--min-contig-len", str(min_len), "--threads", "1",
    ]) != 0
    assert store.stage_files(r_cli, "S1") == []
    assert not (r_cli / "l0seq/ingest_plan/stage=S1").exists()

    # Batched CLI refusal: tiny cpu_min_per_process forces `ingest_steps` (`pass1`, `plan`, ...)
    tiny_prof = dataclasses.replace(prof, cpu_min_per_process=1e-9)
    tiny_pf = tmp_path / "tiny_profile.yaml"
    tiny_pf.write_text(yaml.safe_dump({k: str(v) if isinstance(v, Path) else v
                                       for k, v in dataclasses.asdict(tiny_prof).items()}))
    r_bcli = _base_root("batched_cli_refuse")
    assert bs.build_stage(
        tiny_prof,
        bs.Build(root=r_bcli, stage="S1", release=None, genomes=bad_gen, objects=objs_p,
                 depth="stored", min_contig_len=min_len, threads=1),
        clock=lambda: NOON_MONDAY, load=lambda: 0.0,
    ) != 0
    assert store.stage_files(r_bcli, "S1") == []
    assert not (r_bcli / "l0seq/ingest_plan/stage=S1").exists()

    # Case (f): Matching `genomes.json` succeeds end-to-end through `build-stage --depth called --min-contig-len`
    good_gen = tmp_path / "good_genomes.json"
    good_gen.write_text(json.dumps([
        {"contig_uids": [cu_t4.hex()], "genetic_code": 11, "kind": "virus"},
        {"contig_uids": [cu_lam.hex()], "genetic_code": 11, "kind": "virus"},
    ]))
    r_ok = _base_root("cli_ok")
    assert cli.main([
        "build-stage", "--profile", str(pf), "--root", str(r_ok), "--stage", "S1",
        "--depth", "called", "--genomes", str(good_gen), "--objects", str(objs_p),
        "--min-contig-len", str(min_len), "--threads", "1",
    ]) == 0
    assert bs.stage_complete(r_ok, "S1", "called")
    assert ledger.check_balance(r_ok, "S1", ["ingest", "genecall"]) == []


# --- §5.3 (4): stage without --min-contig-len keeps null base columns and exact unfiltered behavior ---------


@pytest.mark.covers("C-38d", "L-01", "LY-09")
def test_c38d_4_stage_without_filter_writes_null_base_columns_and_matches_unfiltered_oracle(tmp_path):
    paths = [FIX / o["file"] for o in CONTRACT["objects"]]
    objs = [SourceObject(p, CONTRACT["source_name"], CONTRACT["source_release"], CONTRACT["license"], True)
            for p in paths]
    root = tmp_path / "unfiltered"
    ingest(objs, root, "S0", shard_bases=200_000)

    want_ledger, want_digests, want_contigs = _brute_force_ingest_oracle(objs, None)
    got_rows = ledger.read(root, "S0", "ingest")
    assert len(got_rows) == len(objs)
    for r in got_rows:
        assert r["reason_code"] == "ok"
        assert r["n_bases_in"] is None
        assert r["n_bases_rejected"] is None
    got_tuples = [
        (r["object_key"], r["reason_code"], r["n_in"], r["n_kept"], r["n_merged"],
         r["n_rejected"], r["n_bases_in"], r["n_bases_rejected"])
        for r in got_rows
    ]
    assert got_tuples == want_ledger
    loc = pq.read_table(root / "l0seq/seq_locator/stage=S0/part-0.parquet").to_pylist()
    assert {r["raw_digest"] for r in loc} == want_digests
    idx = pq.read_table(root / "l0seq/contig_index/stage=S0/part-0.parquet").to_pylist()
    assert {r["contig_uid"] for r in idx} == want_contigs


# --- §5.3 (5): bounded memory on the 100x synthetic store with short contigs interleaved --------------------


CHILD_INGEST = r"""
import json, sys
from pathlib import Path
from biobox.l0seq import SourceObject, ingest
objs_doc = json.loads(Path(sys.argv[1]).read_text())
objs = [SourceObject(Path(o["path"]), o["source_name"], o["source_release"]) for o in objs_doc]
res = ingest(objs, Path(sys.argv[2]), "S0", shard_bases=int(sys.argv[3]),
             min_contig_len=int(sys.argv[4]), tmp=Path(sys.argv[5]))
print(json.dumps(res, sort_keys=True))
"""


@pytest.mark.covers("C-38d", "LY-06")
def test_c38d_5_100x_synthetic_store_ingest_with_min_contig_len_in_bounded_memory(tmp_path, record_property):
    """Scale micro_s0 100x as `scripts/stress/synth_scale.py` does (~79.7 Mbp in 300 >= L contigs across 2
    source objects), plus 600 short `< L` contigs (200..1,499 bp) interleaved in those objects, and ingest in a
    child process measured by `runstep.measured` / `runstep.step_peak` (AGENTS §4.1, §4.2)."""
    copies = 100
    min_len = 1_500
    shard_bases = 2**26  # 64 MiB per shard, matching test_c34_labels_memory.py
    base_seqs = [fasta.read(FIX / o["file"])[0].seq.upper() for o in CONTRACT["objects"]]

    srcdir = tmp_path / "src"
    srcdir.mkdir()
    rng = np.random.default_rng(1)
    acgt = np.frombuffer(b"ACGT", dtype=np.uint8)

    objs: list[SourceObject] = []
    fh = None
    want_by_obj: list[tuple[int, int, int, int, int]] = []  # (n_in, kept, n_rej, bases_in, bases_rej)
    cur_in = cur_kept = cur_rej = cur_bases_in = cur_bases_rej = 0
    seen_rd: set[bytes] = set()

    for c in range(copies):
        if c % 50 == 0:
            if fh is not None:
                fh.close()
                want_by_obj.append((cur_in, cur_kept, cur_rej, cur_bases_in, cur_bases_rej))
                cur_in = cur_kept = cur_rej = cur_bases_in = cur_bases_rej = 0
            p = srcdir / f"part-{c // 50:04d}.fna"
            fh = open(p, "wb")
            objs.append(SourceObject(p, "synth", "1"))
        # Two short `< 1500` records per copy (including exact `min_len - 1` boundary on every 10th copy)
        for j, slen in enumerate((250 + (c % 500), 1_499 if c % 10 == 0 else 900)):
            s_short = acgt[rng.integers(0, 4, size=slen)].tobytes()
            fh.write(b">short_%d_%d\n%s\n" % (c, j, s_short))
            cur_in += 1
            cur_rej += 1
            cur_bases_in += slen
            cur_bases_rej += slen
        # Three full-length substituted micro_s0 contigs (`>= 48,502 bp >> 1,500 bp`)
        for b_idx, seq in enumerate(base_seqs):
            arr = np.frombuffer(seq, dtype=np.uint8).copy()
            pos = rng.choice(arr.size, size=int(arr.size * 0.03), replace=False)
            arr[pos] = acgt[rng.integers(0, 4, size=pos.size)]
            s = arr.tobytes()
            rd = hashlib.sha256(s).digest()
            fh.write(b">syn_%d_%d\n" % (c, b_idx))
            fh.write(b"\n".join(s[i : i + 80] for i in range(0, len(s), 80)) + b"\n")
            cur_in += 1
            cur_bases_in += len(s)
            if rd not in seen_rd:
                seen_rd.add(rd)
                cur_kept += 1
    assert fh is not None
    fh.close()
    want_by_obj.append((cur_in, cur_kept, cur_rej, cur_bases_in, cur_bases_rej))

    objs_json = tmp_path / "objs.json"
    objs_json.write_text(json.dumps([
        {"path": str(o.path), "source_name": o.source_name, "source_release": o.source_release}
        for o in objs
    ]))
    root = tmp_path / "x100_root"
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    peak_file = tmp_path / "peak-x100.json"

    cmd = runstep.measured(
        [sys.executable, "-c", CHILD_INGEST, str(objs_json), str(root),
         str(shard_bases), str(min_len), str(scratch)],
        peak_file,
        False,
    )
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    assert proc.returncode == 0, proc.stderr[-4000:]
    rss, source = runstep.step_peak(peak_file)
    assert rss is not None, source
    record_property("c38d_x100_peak_bytes", rss)
    print(f"C-38d 100x ingest (--min-contig-len {min_len}): peak {rss / 2**20:.1f} MiB ({source})")
    assert rss <= PEAK_BOUND, (rss, PEAK_BOUND, source)

    # Verify ledger counts against our independent oracle
    expected_by_sha = {}
    for o, counts in zip(objs, want_by_obj):
        sha = ids.sha256_file(o.path)
        n_in, kept, n_rej, b_in, b_rej = counts
        merged = n_in - kept - n_rej
        expected_by_sha[(sha, f"below:{min_len}")] = (n_in, kept, merged, n_rej, b_in, b_rej)
        expected_by_sha[(sha, "ok")] = (n_in, kept, merged, 0, b_in, 0)

    got_ledger = {
        (r["object_key"], r["reason_code"]): (
            r["n_in"], r["n_kept"], r["n_merged"], r["n_rejected"], r["n_bases_in"], r["n_bases_rejected"]
        )
        for r in ledger.read(root, "S0", "ingest")
    }
    assert got_ledger == expected_by_sha
    assert ledger.check_balance(root, "S0", ["ingest"]) == []


# --- write_labels regression: duplicate SourceObject bytes under --min-contig-len ----------------------------


@pytest.mark.covers("C-38d", "C-34", "L-01")
def test_c38d_write_labels_on_filtered_stage_with_duplicate_source_object_bytes(tmp_path):
    """When two SourceObjects share identical file bytes (`object_key`) under distinct `source_name`s on a
    `--min-contig-len`-filtered stage, `ingest` writes sorted `(object_key, reason_code)` ledger rows and
    `store.write_labels` reconciles the stage's `contig_occurrence` record count without collapsing the two
    objects and reproduces ingest's own `contig_label` file byte-for-byte."""
    rng = random.Random(20261120)
    min_len = 1_500
    s_short = _rand_seq(rng, 700)    # < L, dropped in each object
    s_long1 = _rand_seq(rng, 1_800)  # >= L, kept
    s_long2 = _rand_seq(rng, 2_100)  # >= L, kept
    s_other = _rand_seq(rng, 1_600)  # >= L, in a distinct object

    fa_dup = _write_fasta(tmp_path / "src/shared_bytes.fna", [
        ("dup_short", s_short),
        ("dup_long1", s_long1),
        ("dup_long2", s_long2),
    ])
    fa_other = _write_fasta(tmp_path / "src/other.fna", [
        ("other_short", s_short),
        ("other_long", s_other),
    ])
    objs = [
        SourceObject(fa_dup, "src_alpha", "1", "CC0", True),
        SourceObject(fa_dup, "src_beta", "1", "CC0", True),
        SourceObject(fa_other, "src_gamma", "1", "CC0", True),
    ]

    # Oracle root: ingested normally with contig_label written by ingest
    oracle_root = tmp_path / "oracle"
    ingest(objs, oracle_root, "S0", min_contig_len=min_len)
    want_labels = io.part_files(oracle_root / "l0seq/contig_label/stage=S0")
    assert want_labels

    # Unlabelled root: ingested with contig_label suppressed, as before card B-27
    step_root = tmp_path / "unlabelled"
    mp = pytest.MonkeyPatch()
    real_write_parts = io.write_parts
    mp.setattr(
        io,
        "write_parts",
        lambda batches, spec, *x, **kw: [] if spec.name == "contig_label" else real_write_parts(batches, spec, *x, **kw),
    )
    try:
        ingest(objs, step_root, "S0", min_contig_len=min_len)
    finally:
        mp.undo()
    assert io.part_files(step_root / "l0seq/contig_label/stage=S0") == []

    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir()
    res = store.write_labels(step_root, "S0", tmp=scratch_dir)
    # 2 kept records in src_alpha + 2 in src_beta + 1 in src_gamma = 5 occurrence rows across 3 distinct >= L contigs
    assert res["records"] == 5 and res["contigs"] == 3
    got_labels = io.part_files(step_root / "l0seq/contig_label/stage=S0")
    assert [p.read_bytes() for p in got_labels] == [p.read_bytes() for p in want_labels]
    assert ledger.check_balance(step_root, "S0", ["labels"]) == []

