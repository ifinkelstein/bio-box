"""Ring R2: bounded memory in `_family_run` and bounded open files in `postings.run` (card B-62).

Why existing `test_b29_long_contigs_stay_inside_the_memory_limit` did not catch `_family_run` OOM:
that fixture gives only 3 CDSs in 1,000 a family domain (~1,967 domain rows in 115 MB), so the
single-statement `_family_run` hash join stayed tiny. Real annotated releases (`r26.11.2`) have
tens of millions of domain rows.

Both faults are reproduced small here:
1. A store with 400,000 domain rows over 200,000 slots at `memory_limit="25MB"` (before B-62,
   `_family_run` raises `duckdb.OutOfMemoryException`).
2. A store whose `tape_tokens` has 180 row groups, run in a child process with soft
   `RLIMIT_NOFILE=128` (before B-62, `_explode_runs` writes 180 `run-*.parquet` and 180
   `canon-*.parquet` files and `_check_canonical` raises `duckdb.IOException: Too many open files`).
"""

from __future__ import annotations

import hashlib
import os
import resource
import subprocess
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import io
from biobox.l3 import postings, tapes
from biobox.registry.schemas import table

pytestmark = pytest.mark.r2

REL = "rS0.b62"
PROD = "synthetic-b62"
TV = tapes.tape_version([PROD])
LIB = "pfam-38.2"
LEVEL_BITS = {"family": 1, "nr50": 2, "nr90": 3}
POST = "l3/token_postings"


def _digests(root: Path) -> dict[str, str]:
    return {
        p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted((root / POST).rglob("*.parquet"))
    }


def _build_dense_family_store(
    root: Path,
    n_contigs: int,
    slots_per_contig: int,
    n_models: int = 200,
    domains_per_slot: int = 2,
    tape_max_slots: int = 131_072,
    extra_unknown_domain: bool = False,
) -> tuple[list[str], int]:
    """Construct tapes, l1_domains, token_dict, and tape_tokens using vectorized Arrow arrays.

    Every slot is a CDS with a distinct protein_uid, `domains_per_slot` distinct Pfam domains
    (canonical = first domain, higher bitscore), plus NR50 and NR90 tokens.
    """
    stages = ["S0", "S1"] if n_contigs >= 2 else ["S0"]
    n_slots_total = n_contigs * slots_per_contig
    # deterministic 16-byte contig_uids sorted lexicographically (big-endian counter)
    cu_bytes = [b"B62C" + c.to_bytes(12, "big") for c in range(n_contigs)]
    # deterministic 16-byte protein_uids and feature_uids
    slot_ids = np.arange(n_slots_total, dtype=np.uint64)
    # 16-byte fixed binary via structured array: 8-byte prefix + 8-byte big-endian id
    pu_buf = np.empty(n_slots_total, dtype=[("tag", "S8"), ("id", ">u8")])
    pu_buf["tag"] = b"B62PROT!"
    pu_buf["id"] = slot_ids
    pu_flat = pa.FixedSizeBinaryArray.from_buffers(
        pa.binary(16), n_slots_total, [None, pa.py_buffer(pu_buf.tobytes())]
    )
    fu_buf = np.empty(n_slots_total, dtype=[("tag", "S8"), ("id", ">u8")])
    fu_buf["tag"] = b"B62FEAT!"
    fu_buf["id"] = slot_ids
    fu_flat = pa.FixedSizeBinaryArray.from_buffers(
        pa.binary(16), n_slots_total, [None, pa.py_buffer(fu_buf.tobytes())]
    )

    # token_dict: family (source=1), nr50 (source=0), nr90 (source=0)
    fam_base = (LEVEL_BITS["family"] << 61) | (1 << 56)
    nr50_base = LEVEL_BITS["nr50"] << 61
    nr90_base = LEVEL_BITS["nr90"] << 61
    models = [f"PF{m:05d}" for m in range(n_models)]

    # Assign 2 domains per protein: m0 = slot_id % n_models (env_start=10, bitscore=80, canonical),
    # m1 = (slot_id + 1) % n_models (env_start=200, bitscore=40)
    m0 = (slot_ids % n_models).astype(np.uint64)
    fam_tok_flat = pa.array(fam_base | m0, type=pa.uint64())
    nr50_tok_flat = pa.array(nr50_base | (slot_ids % 500), type=pa.uint64())
    nr90_tok_flat = pa.array(nr90_base | (slot_ids % 2000), type=pa.uint64())

    # Write token_dict
    td_spec = table("token_dict")
    fam_counts = np.bincount(m0, minlength=n_models)
    for d_i in range(1, domains_per_slot):
        mi = ((slot_ids + d_i) % n_models).astype(np.uint64)
        fam_counts = fam_counts + np.bincount(mi, minlength=n_models)
    nr50_counts = np.bincount(slot_ids % 500, minlength=500)
    nr90_counts = np.bincount(slot_ids % 2000, minlength=2000)

    td_tids, td_lvls, td_srcs, td_keys, td_ns = [], [], [], [], []
    for m in range(n_models):
        if fam_counts[m]:
            td_tids.append(int(fam_base | m))
            td_lvls.append("family")
            td_srcs.append(LIB)
            td_keys.append(models[m])
            td_ns.append(int(fam_counts[m]))
    for k in range(500):
        if nr50_counts[k]:
            td_tids.append(int(nr50_base | k))
            td_lvls.append("nr50")
            td_srcs.append("mmseqs2")
            td_keys.append(f"{k:032x}")
            td_ns.append(int(nr50_counts[k]))
    for k in range(2000):
        if nr90_counts[k]:
            td_tids.append(int(nr90_base | k))
            td_lvls.append("nr90")
            td_srcs.append("mmseqs2")
            td_keys.append(f"{k:032x}")
            td_ns.append(int(nr90_counts[k]))
    td_tbl = pa.table(
        {
            "token_id": pa.array(td_tids, pa.uint64()),
            "level": pa.array(td_lvls, pa.string()),
            "source": pa.array(td_srcs, pa.string()),
            "key": pa.array(td_keys, pa.string()),
            "n_slots": pa.array(td_ns, pa.uint64()),
        },
        schema=td_spec.schema,
    )
    io.write_table([td_tbl], td_spec, root / f"l3/token_dict/release={REL}/part-0.parquet")

    # Write tapes (split across stages) and l1_domains (split across stages)
    tape_spec = table("contig_tapes")
    dom_spec = table("l1_domains")
    starts_one = np.arange(slots_per_contig, dtype=np.uint32) * 1000
    ends_one = starts_one + 900
    strands_one = np.ones(slots_per_contig, dtype=np.int8)
    delta_one = np.full(slots_per_contig, 100, dtype=np.int32)
    delta_one[0] = 0
    junc_one = np.ones(slots_per_contig, dtype=np.uint8)
    orfs_start_one = np.arange(slots_per_contig, dtype=np.uint16)
    orfs_end_one = np.arange(slots_per_contig - 1, -1, -1, dtype=np.uint16)
    flags_one = np.zeros(slots_per_contig, dtype=np.uint8)
    kinds_one = np.ones(slots_per_contig, dtype=np.uint8)

    mid = n_contigs // 2 if len(stages) == 2 else n_contigs
    stage_ranges = [(stages[0], 0, mid)] + ([(stages[1], mid, n_contigs)] if len(stages) == 2 else [])
    for st, c0, c1 in stage_ranges:
        nc = c1 - c0
        ns = nc * slots_per_contig
        s0, s1 = c0 * slots_per_contig, c1 * slots_per_contig
        offsets = pa.array(np.arange(0, ns + 1, slots_per_contig, dtype=np.int32))
        tape_tbl = pa.table(
            {
                "contig_uid": pa.array(cu_bytes[c0:c1], pa.binary(16)),
                "contig_len": pa.array(np.full(nc, slots_per_contig * 1000, dtype=np.int64)),
                "is_circular": pa.array(np.zeros(nc, dtype=bool)),
                "n_slots": pa.array(np.full(nc, slots_per_contig, dtype=np.uint32)),
                "feature_kinds": pa.ListArray.from_arrays(offsets, pa.array(np.tile(kinds_one, nc), pa.uint8())),
                "strands": pa.ListArray.from_arrays(offsets, pa.array(np.tile(strands_one, nc), pa.int8())),
                "start_bps": pa.ListArray.from_arrays(offsets, pa.array(np.tile(starts_one, nc), pa.uint32())),
                "end_bps": pa.ListArray.from_arrays(offsets, pa.array(np.tile(ends_one, nc), pa.uint32())),
                "delta_bp_prev_slots": pa.ListArray.from_arrays(offsets, pa.array(np.tile(delta_one, nc), pa.int32())),
                "junction_topology_prevs": pa.ListArray.from_arrays(offsets, pa.array(np.tile(junc_one, nc), pa.uint8())),
                "orfs_to_contig_starts": pa.ListArray.from_arrays(offsets, pa.array(np.tile(orfs_start_one, nc), pa.uint16())),
                "orfs_to_contig_ends": pa.ListArray.from_arrays(offsets, pa.array(np.tile(orfs_end_one, nc), pa.uint16())),
                "slot_flags": pa.ListArray.from_arrays(offsets, pa.array(np.tile(flags_one, nc), pa.uint8())),
                "protein_uids": pa.ListArray.from_arrays(offsets, pu_flat.slice(s0, ns)),
                "feature_uids": pa.ListArray.from_arrays(offsets, fu_flat.slice(s0, ns)),
            },
            schema=tape_spec.schema,
        )
        io.write_table(
            [tape_tbl],
            tape_spec,
            root / f"l3/tapes/tape_version={TV}/stage={st}/part-0.parquet",
            max_slots=tape_max_slots,
            slots_column="n_slots",
        )

        # Build l1_domains for this stage, already sorted by (protein_uid, env_start, model_acc)
        dom_pus, dom_accs, dom_envs, dom_scores = [], [], [], []
        stage_ids = slot_ids[s0:s1]
        pu_stage = pu_flat.slice(s0, ns)
        # Interleave domains per protein so sort order (protein_uid, env_start, model_acc) holds directly
        rep_idx = np.repeat(np.arange(ns, dtype=np.int64), domains_per_slot)
        dom_pu = pu_stage.take(pa.array(rep_idx))
        dom_idx = np.tile(np.arange(domains_per_slot, dtype=np.uint64), ns)
        dom_model_num = (np.repeat(stage_ids, domains_per_slot) + dom_idx) % n_models
        model_arr = pa.array(models, pa.string())
        dom_acc = model_arr.take(pa.array(dom_model_num, pa.int64()))
        dom_env = pa.array((10 + dom_idx * 190).astype(np.uint32))
        dom_end = pa.array((110 + dom_idx * 190).astype(np.uint32))
        dom_score = pa.array(np.where(dom_idx == 0, 80.0, 40.0).astype(np.float32))
        nd = ns * domains_per_slot
        dom_tbl = pa.table(
            {
                "protein_uid": dom_pu,
                "model_acc": dom_acc,
                "model_name": dom_acc,
                "model_len": pa.array(np.full(nd, 100, dtype=np.uint32)),
                "seq_bitscore": dom_score,
                "seq_evalue_zref": pa.array(np.full(nd, 1e-10, dtype=np.float64)),
                "dom_bitscore": dom_score,
                "dom_ievalue_zref": pa.array(np.full(nd, 1e-10, dtype=np.float64)),
                "env_start": dom_env,
                "env_end": dom_end,
                "ali_start": dom_env,
                "ali_end": dom_end,
                "hmm_from": pa.array(np.zeros(nd, dtype=np.uint32)),
                "hmm_to": pa.array(np.full(nd, 100, dtype=np.uint32)),
                "cutoff": pa.array(["gathering"] * nd, pa.string()),
                "z": pa.array(np.full(nd, 1e8, dtype=np.float64)),
                "dom_z": pa.array(np.full(nd, 1e8, dtype=np.float64)),
            },
            schema=dom_spec.schema,
        )
        if extra_unknown_domain and st == stages[-1]:
            extra = pa.table(
                {
                    "protein_uid": pu_stage.slice(ns - 1, 1),
                    "model_acc": pa.array(["PF99999"], pa.string()),
                    "model_name": pa.array(["PF99999"], pa.string()),
                    "model_len": pa.array([100], pa.uint32()),
                    "seq_bitscore": pa.array([30.0], pa.float32()),
                    "seq_evalue_zref": pa.array([1e-5], pa.float64()),
                    "dom_bitscore": pa.array([30.0], pa.float32()),
                    "dom_ievalue_zref": pa.array([1e-5], pa.float64()),
                    "env_start": pa.array([500], pa.uint32()),
                    "env_end": pa.array([600], pa.uint32()),
                    "ali_start": pa.array([500], pa.uint32()),
                    "ali_end": pa.array([600], pa.uint32()),
                    "hmm_from": pa.array([0], pa.uint32()),
                    "hmm_to": pa.array([100], pa.uint32()),
                    "cutoff": pa.array(["gathering"], pa.string()),
                    "z": pa.array([1e8], pa.float64()),
                    "dom_z": pa.array([1e8], pa.float64()),
                },
                schema=dom_spec.schema,
            )
            dom_tbl = pa.concat_tables([dom_tbl, extra])
        io.write_table([dom_tbl], dom_spec, root / f"l1/domains/library={LIB}/stage={st}/part-0.parquet")

    # Write tape_tokens
    tok_spec = table("tape_tokens")
    all_offsets = pa.array(np.arange(0, n_slots_total + 1, slots_per_contig, dtype=np.int32))
    # domain_family_tokens: list<list<uint64>>
    inner_vals = np.empty(n_slots_total * domains_per_slot, dtype=np.uint64)
    for d_i in range(domains_per_slot):
        inner_vals[d_i::domains_per_slot] = fam_base | ((slot_ids + d_i) % n_models).astype(np.uint64)
    inner_off = pa.array(np.arange(0, n_slots_total * domains_per_slot + 1, domains_per_slot, dtype=np.int32))
    dom_fam_flat = pa.ListArray.from_arrays(inner_off, pa.array(inner_vals, pa.uint64()))
    tok_tbl = pa.table(
        {
            "contig_uid": pa.array(cu_bytes, pa.binary(16)),
            "family_tokens": pa.ListArray.from_arrays(all_offsets, fam_tok_flat),
            "nr50_tokens": pa.ListArray.from_arrays(all_offsets, nr50_tok_flat),
            "nr90_tokens": pa.ListArray.from_arrays(all_offsets, nr90_tok_flat),
            "domain_family_tokens": pa.ListArray.from_arrays(all_offsets, dom_fam_flat),
        },
        schema=tok_spec.schema,
    )
    io.write_table(
        [tok_tbl],
        tok_spec,
        root / f"l3/tape_tokens/tape_version={TV}/release={REL}/part-0.parquet",
        max_slots=tape_max_slots,
        slots_column="nr90_tokens",
    )
    return stages, n_slots_total


@pytest.mark.covers("K-02", "LY-02", "D-01")
def test_b62_family_run_stays_inside_memory_limit(tmp_path):
    """Fault 1 (card B-62): 1,600,000 domain rows over 400,000 slots at memory_limit='55MB'.

    Before B-62, `_family_run` ran a single DuckDB statement over all domain rows and raised
    `duckdb.OutOfMemoryException` at 55 MB (while `_check` passed). After B-62, `_family_run`
    partitions by protein range and writes postings byte-identical to a 2 GB run.
    """
    root = tmp_path / "root"
    stages, n_slots = _build_dense_family_store(root, n_contigs=400, slots_per_contig=1000, domains_per_slot=4)
    out_low = tmp_path / "out_low"
    out_hi = tmp_path / "out_hi"
    res_low = postings.run(root, REL, stages, producers=[PROD], memory_limit="55MB", tmp=tmp_path, out=out_low)
    assert res_low["rows_by_level"] == {"family": 4 * n_slots, "nr50": n_slots, "nr90": n_slots}
    assert res_low["rows"] == 6 * n_slots
    res_hi = postings.run(root, REL, stages, producers=[PROD], memory_limit="2GB", tmp=tmp_path, out=out_hi)
    assert res_hi == res_low
    assert _digests(out_low) == _digests(out_hi)


@pytest.mark.covers("K-02", "LY-02", "D-01")
def test_b62_many_row_groups_stay_under_nofile_limit(tmp_path):
    """Fault 2 (card B-62): 180 tape_tokens row groups under soft RLIMIT_NOFILE=128.

    Before B-62, `_explode_runs` wrote one `run-*.parquet` and one `canon-*.parquet` per
    tape_tokens row group (180 each), and `_check_canonical` / final sort opened all of them at
    once, raising `duckdb.IOException: Too many open files`.
    """
    root = tmp_path / "root"
    # 180 contigs x 8 slots each, tape_max_slots=8 -> 180 row groups in tape_tokens
    stages, n_slots = _build_dense_family_store(
        root, n_contigs=180, slots_per_contig=8, n_models=20, domains_per_slot=2, tape_max_slots=8
    )
    tok_pf = pq.ParquetFile(root / f"l3/tape_tokens/tape_version={TV}/release={REL}/part-0.parquet")
    assert tok_pf.metadata.num_row_groups == 180

    out_ref = tmp_path / "out_ref"
    res_ref = postings.run(root, REL, stages, producers=[PROD], tmp=tmp_path, out=out_ref)
    assert res_ref["rows"] == 4 * n_slots

    out_sub = tmp_path / "out_sub"
    code = (
        "import resource, sys\n"
        "from pathlib import Path\n"
        "from biobox.l3 import postings\n"
        "soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)\n"
        "resource.setrlimit(resource.RLIMIT_NOFILE, (min(128, hard), hard))\n"
        "res = postings.run(Path(sys.argv[1]), sys.argv[2], sys.argv[3].split(','),\n"
        "                   producers=[sys.argv[4]], tmp=Path(sys.argv[5]), out=Path(sys.argv[6]))\n"
        "assert res['rows'] == int(sys.argv[7]), res\n"
    )
    env = {**os.environ}
    proc = subprocess.run(
        [sys.executable, "-c", code, str(root), REL, ",".join(stages), PROD, str(tmp_path), str(out_sub), str(4 * n_slots)],
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )
    assert proc.returncode == 0, f"stdout={proc.stdout}\nstderr={proc.stderr}"
    assert _digests(out_sub) == _digests(out_ref)


@pytest.mark.covers("C-03", "LY-02")
def test_b62_missing_family_token_still_raises_across_ranges(tmp_path):
    """A domain whose model has no token_dict entry raises ValueError even when `_family_run` splits into k > 1 ranges."""
    root = tmp_path / "root"
    stages, _ = _build_dense_family_store(
        root, n_contigs=60, slots_per_contig=1000, n_models=50, domains_per_slot=1, extra_unknown_domain=True
    )
    out = tmp_path / "out"
    with pytest.raises(ValueError, match=r"PF99999.*token_dict"):
        postings.run(root, REL, stages, producers=[PROD], memory_limit="50MB", tmp=tmp_path, out=out)
    assert not (out / POST).exists()


@pytest.mark.covers("K-02", "LY-02", "D-01")
def test_b62_coalesce_bounds_open_files_even_when_many_runs_are_produced(tmp_path):
    """Even when `_RUN_MAX_ROWS` forces 180 run/canon files, `_coalesce_parquet` keeps open descriptors
    bounded under `RLIMIT_NOFILE=128` and produces byte-identical postings."""
    root = tmp_path / "root"
    stages, n_slots = _build_dense_family_store(
        root, n_contigs=180, slots_per_contig=8, n_models=20, domains_per_slot=2, tape_max_slots=8
    )
    out_ref = tmp_path / "out_ref"
    postings.run(root, REL, stages, producers=[PROD], tmp=tmp_path, out=out_ref)

    out_sub = tmp_path / "out_sub"
    code = (
        "import resource, sys\n"
        "from pathlib import Path\n"
        "from biobox.l3 import postings\n"
        "postings._RUN_MAX_ROWS = 1\n"
        "postings._MAX_READ_FILES = 16\n"
        "soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)\n"
        "resource.setrlimit(resource.RLIMIT_NOFILE, (min(128, hard), hard))\n"
        "res = postings.run(Path(sys.argv[1]), sys.argv[2], sys.argv[3].split(','),\n"
        "                   producers=[sys.argv[4]], tmp=Path(sys.argv[5]), out=Path(sys.argv[6]))\n"
        "assert res['rows'] == int(sys.argv[7]), res\n"
    )
    proc = subprocess.run(
        [sys.executable, "-c", code, str(root), REL, ",".join(stages), PROD, str(tmp_path), str(out_sub), str(4 * n_slots)],
        capture_output=True,
        text=True,
        env={**os.environ},
        check=False,
    )
    assert proc.returncode == 0, f"stdout={proc.stdout}\nstderr={proc.stderr}"
    assert _digests(out_sub) == _digests(out_ref)
