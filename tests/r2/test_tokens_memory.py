"""Ring R2: `tokens.run` inside a DuckDB memory limit scaled from the 50 Gbp release (card B-46).

Release `rM.1` on the pod (54,265 contigs, 48,438,455 slots, 45,815,843 proteins, 31,183,876 NR90 and
8,192,065 NR50 clusters, 16,766,950 Pfam domain rows on 8,233,292 proteins; measured, wilkcomp02,
2026-11-12) ran the tokens step out of its 2 GB DuckDB limit in the slot join: one LEFT JOIN of the
unnested slots to the 45.8 M cluster rows of a level failed alone, with its inputs materialised as
Parquet (card B-46's report). This store has that shape at about 1/40:
long and short contigs over two stages, a protein per slot with a few shared, two cluster rows per
protein at the release's NR90 and NR50 ratios, domains on about a fifth of the proteins from two
libraries, non-CDS slots and tapes without slots; the DuckDB limit is scaled from 2 GB by the store's
share of rM.1's slots. The oracle is the in-memory loop of test_micro_s0_l3 (`_oracle_tokens`).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pytest

from biobox import io
from biobox.l3 import tapes, tokens
from biobox.registry.schemas import table

sys.path.insert(0, str(Path(__file__).parent))
import test_micro_s0_l3 as micro  # noqa: E402

pytestmark = pytest.mark.r2

REL = "rS0.b46"
PRODUCER = "synthetic-b46"
GROUPS = (("S0", 150, 4_000), ("S1", 300, 2_000))  # (stage, contigs, slots each): 1,200,000 slots
RM1_SLOTS = 48_438_455  # rM.1 (measured, wilkcomp02, 2026-11-12)
LIMIT = "50MB"  # 2 GB x 1,200,000 / 48,438,455 slots = 49.5 MB, rounded up (derived)
LIBS = ("aaa-1", "synfam-1")


def _uids(rng: np.random.Generator, n: int) -> pa.Array:
    raw = rng.bytes(16 * n)
    return pa.array([raw[i:i + 16] for i in range(0, 16 * n, 16)], pa.binary())


def _list(offsets: np.ndarray, values: pa.Array) -> pa.Array:
    return pa.ListArray.from_arrays(pa.array(offsets, pa.int32()), values)


def _write(spec_name: str, t: pa.Table, path: Path, sort: list[str], **kw) -> None:
    spec = table(spec_name)
    io.write_table(t.sort_by([(c, "ascending") for c in sort]).cast(spec.schema), spec, path, **kw)


def _synthetic(root: Path, seed: int = 46) -> tuple[str, list[str]]:
    """Tapes, cluster members and domains of `GROUPS`, written directly. Returns (tape version, stages)."""
    rng = np.random.default_rng(seed)
    tv = tapes.tape_version([PRODUCER])
    n_all = sum(c * s for _, c, s in GROUPS)
    # a protein per slot, 5% of slots reusing an earlier protein (rM.1: 1.057 slots per protein), 1% non-CDS
    pool = _uids(rng, n_all)
    reuse = rng.random(n_all) < 0.05
    pick = np.where(reuse, rng.integers(0, n_all, n_all), np.arange(n_all))
    prot = pool.take(pa.array(pick))
    prot = pa.array([None if x else p for p, x in zip(prot.to_pylist(), rng.random(n_all) < 0.01)], pa.binary())
    at, by_stage = 0, {}
    for stage, n_contigs, n_slots in GROUPS:
        n = n_contigs + 3  # three tapes without slots per stage
        lens = np.r_[np.full(n_contigs, n_slots), np.zeros(3, int)]
        off = np.r_[0, np.cumsum(lens)]
        z = lambda typ: _list(off, pa.array(np.zeros(off[-1], np.int64)).cast(typ))  # noqa: E731
        by_stage[stage] = (stage, pa.table({
            "contig_uid": _uids(rng, n), "contig_len": pa.array(lens * 1000, pa.int64()),
            "is_circular": pa.array(np.arange(n) % 3 == 0), "n_slots": pa.array(lens, pa.uint32()),
            "feature_kinds": z(pa.uint8()), "strands": _list(off, pa.array(np.ones(off[-1], np.int8))),
            "start_bps": z(pa.uint32()), "end_bps": z(pa.uint32()), "delta_bp_prev_slots": z(pa.int32()),
            "junction_topology_prevs": z(pa.uint8()), "orfs_to_contig_starts": z(pa.uint16()),
            "orfs_to_contig_ends": z(pa.uint16()), "slot_flags": z(pa.uint8()),
            "protein_uids": _list(off, prot.slice(at, off[-1])), "feature_uids": _list(off, _uids(rng, off[-1]))}))
        at += off[-1]
    for stage, t in by_stage.values():
        _write("contig_tapes", t, root / f"l3/tapes/tape_version={tv}/stage={stage}/part-0.parquet", ["contig_uid"],
               max_slots=131_072, slots_column="n_slots")  # AGENTS §4.6, as tapes.run writes them
    # clusters: every protein of the pool at NR90 and NR50 (rM.1: 0.68 and 0.18 clusters per protein), and 1% more
    # members on no slot
    members = pa.concat_arrays([pa.array(sorted(set(prot.drop_null().to_pylist())), pa.binary()), _uids(rng, n_all // 100)])
    n_m = len(members)
    rows = []
    for level, ratio in (("nr90", 0.68), ("nr50", 0.18)):
        cl = _uids(rng, int(n_m * ratio))
        rows.append(pa.table({"member_uid": members, "level": pa.array([level] * n_m),
                              "cluster_uid": cl.take(pa.array(rng.integers(0, len(cl), n_m))),
                              "rep_uid": members, "cluster_params_hash": pa.array(["synthetic"] * n_m)}))
    _write("cluster_member", pa.concat_tables(rows), root / f"l2/cluster_member/release={REL}/part-0.parquet",
           ["level", "cluster_uid", "member_uid"])
    # domains on 20% of the proteins, one to three each, from 3,000 models of one library and 40 of a second;
    # some domains with equal scores, and a model twice on one protein
    with_dom = members.filter(pa.array(rng.random(n_m) < 0.2))
    per = rng.integers(1, 4, len(with_dom))
    pu = with_dom.take(pa.array(np.repeat(np.arange(len(with_dom)), per)))
    n_d = len(pu)
    lib = rng.random(n_d) < 0.05
    acc = np.where(lib, rng.integers(0, 40, n_d), rng.integers(0, 3_000, n_d))
    score = rng.integers(10, 30, n_d).astype(np.float32)
    env = rng.integers(0, 400, n_d).astype(np.uint32)
    for name, mask in ((LIBS[0], lib), (LIBS[1], ~lib)):
        for k, (stage, _) in enumerate(by_stage.values()):
            m = mask & (np.arange(n_d) % 2 == k)
            n = int(m.sum())
            accs = pa.array([f"{'A' if name == LIBS[0] else 'M'}{a:05d}" for a in acc[m]])
            full = lambda v, typ: pa.array(np.full(n, v)).cast(typ)  # noqa: E731
            t = pa.table({
                "protein_uid": pu.filter(pa.array(m)), "model_acc": accs, "model_name": accs,
                "model_len": full(100, pa.uint32()), "seq_bitscore": pa.array(score[m]),
                "seq_evalue_zref": full(1e-9, pa.float64()), "dom_bitscore": pa.array(score[m]),
                "dom_ievalue_zref": full(1e-9, pa.float64()), "env_start": pa.array(env[m]),
                "env_end": pa.array(env[m] + 20), "ali_start": pa.array(env[m]), "ali_end": pa.array(env[m] + 20),
                "hmm_from": full(0, pa.uint32()), "hmm_to": full(20, pa.uint32()),
                "cutoff": pa.array(["gathering"] * n), "z": full(1e8, pa.float64()), "dom_z": full(1e8, pa.float64())})
            _write("l1_domains", t, root / f"l1/domains/library={name}/stage={stage}/part-0.parquet",
                   ["protein_uid", "env_start", "model_acc"])
    return tv, list(by_stage)


@pytest.mark.covers("K-02", "C-04")
def test_b46_tokens_stay_inside_a_scaled_memory_limit(tmp_path):
    """1.2 M slots, 2.3 M cluster rows: the tokens are built inside the scaled DuckDB limit (the builder before
    card B-46 raised duckdb.OutOfMemoryException in its slot join here) and equal the oracle."""
    root = tmp_path / "r"
    tv, stages = _synthetic(root)
    res = tokens.run(root, REL, stages, producers=[PRODUCER], memory_limit=LIMIT, tmp=tmp_path)
    got = micro._written_tokens(root, REL, tv)
    want = micro._oracle_tokens(root, REL, list(LIBS), tv)
    assert len(got[1]) == len(want[1]) == sum(c + 3 for _, c, _ in GROUPS)
    assert got == want
    assert res["tokens"]["nr90"] > res["tokens"]["nr50"] > res["tokens"]["family"] == 3_040
