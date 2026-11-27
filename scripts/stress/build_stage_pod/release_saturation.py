#!/usr/bin/env python3
"""Release-level protein-family and syntenic-pair saturation tracker (card A-33, DEC-86 (1)).

Three subcommands:

    release_saturation.py measure --root ROOT --release REL [--earlier PREV] \
        --out DIR --tmp DIR [--threads 4] [--buckets B] [--stats-out FILE]
        Read-only on ROOT, reading ONLY explicit file lists from manifests/release=<REL>.json
        (AGENTS.md §3.9; no globs). Computes deterministic stage-prefix collector curves and
        marginal yields for:
          (a) NR50 clusters over all ORFs and complete ORFs (`dist_to_contig_end_bp > 0`),
              deduplicated to distinct genome incidences (`contig_dim.genome_uid`), with
              Good-Turing sample coverage `c_hat` and bias-corrected Chao1 richness;
          (b) marginal actionable multi-genome complete-ORF NR50 family yield (`k >= 3` and
              `k >= 5` distinct genomes) decomposed into `brand_new_ge_k` + `promoted_to_ge_k`;
          (c) strand-canonical adjacent CDS NR50 gene pairs from `contig_tapes` + `tape_tokens` +
              `contig_dim` (`k >= 1`, `f1`, `f2`, `k >= 3`, `k >= 5`, `c_hat_pairs`, `chao1_pairs`).
        Writes deterministic `saturation.json` and `saturation.txt` under `--out`.

    release_saturation.py history [--data JSON] [--add SAT.json --host H --date YYYY-MM-DD \
        --command CMD [--stats STATS.json] [--set-collector]] [--diversity JSON] [--check]
        Updates or verifies `planning/measurements/release_saturation_history.json` with
        AGENTS.md §6 labels (`measured`, `derived`, `estimate`) and sources.

    release_saturation.py render --out PREFIX [--data JSON]
        Renders `PREFIX.png` (300 DPI) and deterministic `PREFIX.svg` (Wilke / Okabe-Ito
        4-panel figure, 1 thread in the pod image).
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import contextlib
from dataclasses import dataclass
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time
from typing import Iterator

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_HISTORY = REPO_ROOT / "planning" / "measurements" / "release_saturation_history.json"
DEFAULT_DIVERSITY = REPO_ROOT / "planning" / "measurements" / "release_diversity.json"

LABELS = ("measured", "derived", "estimate")
UID_DTYPE = np.dtype([("hi", "<u8"), ("lo", "<u8")])

CLASS_COLORS = {
    "isolate": "#8C96A0",
    "uncultured": "#0072B2",
    "virus": "#E69F00",
    "plasmid": "#009E73",
    "metagenome": "#CC79A7",
    "organelle": "#F0E442",
    "other": "#555555",
}

SLICE_MC_ROWS = 1_500_000
TARGET_BUCKET_ROWS = 2_500_000
STREAM_BATCH_ROWS = 65_536
TAPE_BATCH_ROWS = 16_384


# ---------------------------------------------------------------------------------------------------------------
# Resource measurement (/proc/self/status VmHWM) and atomic file writing
# ---------------------------------------------------------------------------------------------------------------


def read_vmhwm_kib() -> int:
    """Peak resident set size (`VmHWM`) in KiB from `/proc/self/status` on Linux, else `ru_maxrss` fallback."""
    status = Path("/proc/self/status")
    if status.is_file():
        for line in status.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("VmHWM:"):
                parts = line.split()
                if len(parts) >= 2 and parts[1].isdigit():
                    return int(parts[1])
    import resource

    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(rss // 1024) if sys.platform == "darwin" else int(rss)


def atomic_write_text(path: Path, text: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


# ---------------------------------------------------------------------------------------------------------------
# 128-bit BLAKE3 UID conversion and primitive-<u8 accelerated exact lookup
# ---------------------------------------------------------------------------------------------------------------


def _bytes_to_uid_scalar(raw16: bytes) -> np.void:
    arr = np.empty(1, dtype=UID_DTYPE)
    arr["hi"][0] = int.from_bytes(raw16[:8], "big")
    arr["lo"][0] = int.from_bytes(raw16[8:], "big")
    return arr[0]


def _uid_to_pair(u: np.void) -> tuple[int, int]:
    return int(u["hi"]), int(u["lo"])


def _pair_to_uid(pair: tuple[int, int] | None) -> np.void | None:
    if pair is None:
        return None
    arr = np.empty(1, dtype=UID_DTYPE)
    arr["hi"][0] = pair[0]
    arr["lo"][0] = pair[1]
    return arr[0]


def _uid_array(col: pa.Array | pa.ChunkedArray) -> np.ndarray:
    if isinstance(col, pa.ChunkedArray):
        col = col.combine_chunks()
    n = len(col)
    if n == 0:
        return np.empty(0, dtype=UID_DTYPE)
    if col.null_count:
        raise ValueError("unexpected null in 16-byte UID column")
    buf = col.buffers()[1]
    be = np.frombuffer(buf, dtype=">u8", offset=col.offset * 16, count=n * 2).reshape(n, 2)
    out = np.empty(n, dtype=UID_DTYPE)
    out["hi"] = be[:, 0]
    out["lo"] = be[:, 1]
    return out


def _uid_eq(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return (a["hi"] == b["hi"]) & (a["lo"] == b["lo"])


def _uid_scalar_gt(a_elem: np.void, b_elem: np.void) -> bool:
    ah, bh = int(a_elem["hi"]), int(b_elem["hi"])
    if ah != bh:
        return ah > bh
    return int(a_elem["lo"]) > int(b_elem["lo"])


def _uid_scalar_le(a_elem: np.void, b_elem: np.void) -> bool:
    return not _uid_scalar_gt(a_elem, b_elem)


def _lookup_sorted_uids(
    sorted_keys: np.ndarray,
    queries: np.ndarray,
    label: str,
    *,
    sorted_hi: np.ndarray | None = None,
    sorted_lo: np.ndarray | None = None,
) -> np.ndarray:
    """Exact binary search of `queries` in `sorted_keys` (`UID_DTYPE`) via primitive `<u8` fast path."""
    if len(queries) == 0:
        return np.empty(0, dtype=np.int64)
    n_keys = len(sorted_keys)
    if n_keys == 0:
        raise ValueError(f"{label}: lookup against empty key table ({len(queries)} queries)")
    sk_hi = sorted_keys["hi"] if sorted_hi is None else sorted_hi
    sk_lo = sorted_keys["lo"] if sorted_lo is None else sorted_lo
    q_hi = queries["hi"]
    q_lo = queries["lo"]

    pos = np.searchsorted(sk_hi, q_hi, side="left")
    clipped = np.minimum(pos, n_keys - 1)
    exact = (pos < n_keys) & (sk_hi[clipped] == q_hi) & (sk_lo[clipped] == q_lo)
    if not bool(np.all(exact)):
        rem = np.flatnonzero(~exact)
        pos_rem = np.searchsorted(sorted_keys, queries[rem], side="left")
        clip_rem = np.minimum(pos_rem, n_keys - 1)
        ok_rem = (pos_rem < n_keys) & _uid_eq(sorted_keys[clip_rem], queries[rem])
        if not bool(np.all(ok_rem)):
            bad_n = int(len(ok_rem) - np.count_nonzero(ok_rem))
            raise ValueError(f"{label}: {bad_n} of {len(queries)} UIDs missing from sorted target")
        pos[rem] = pos_rem
    return pos


# ---------------------------------------------------------------------------------------------------------------
# Sequential sorted Parquet cursor with footer min/max row-group pruning
# ---------------------------------------------------------------------------------------------------------------


class _SortedParquetStream:
    """Sequential cursor over pre-sorted Parquet files in `(lo_uid_exclusive, hi_uid_inclusive]`."""

    def __init__(
        self,
        files: list[Path],
        columns: list[str],
        uid_col: str,
        *,
        lo_uid_exclusive: np.void | None = None,
        hi_uid_inclusive: np.void | None = None,
        batch_size: int = STREAM_BATCH_ROWS,
    ) -> None:
        self._files = list(files)
        self._columns = list(columns)
        self._uid_col = uid_col
        self._lo = lo_uid_exclusive
        self._hi = hi_uid_inclusive
        self._batch_size = batch_size
        self._iter = self._iter_batches()
        self._batch: pa.RecordBatch | None = None
        self._uids: np.ndarray = np.empty(0, dtype=UID_DTYPE)
        self._pos = 0
        self._advance()

    def _iter_batches(self) -> Iterator[pa.RecordBatch]:
        for f in self._files:
            pf = pq.ParquetFile(f)
            meta = pf.metadata
            uid_idx = pf.schema_arrow.names.index(self._uid_col)
            rgs: list[int] = []
            for rg_i in range(meta.num_row_groups):
                st = meta.row_group(rg_i).column(uid_idx).statistics
                if st is not None and st.has_min_max and isinstance(st.min, bytes) and len(st.min) == 16:
                    rg_min = _bytes_to_uid_scalar(st.min)
                    rg_max = _bytes_to_uid_scalar(st.max)
                    if self._lo is not None and _uid_scalar_le(rg_max, self._lo):
                        continue
                    if self._hi is not None and _uid_scalar_gt(rg_min, self._hi):
                        continue
                rgs.append(rg_i)
            if not rgs:
                continue
            for b in pf.iter_batches(batch_size=self._batch_size, row_groups=rgs, columns=self._columns):
                if b.num_rows:
                    yield b

    def _advance(self) -> None:
        while True:
            try:
                b = next(self._iter)
            except StopIteration:
                self._batch = None
                self._uids = np.empty(0, dtype=UID_DTYPE)
                self._pos = 0
                return
            uids = _uid_array(b.column(self._uid_col))
            start = 0
            end = len(uids)
            if self._lo is not None:
                if _uid_scalar_le(uids[-1], self._lo):
                    continue
                if _uid_scalar_le(uids[0], self._lo):
                    start = int(np.searchsorted(uids, self._lo, side="right"))
            if self._hi is not None:
                if _uid_scalar_gt(uids[start], self._hi):
                    self._batch = None
                    self._uids = np.empty(0, dtype=UID_DTYPE)
                    self._pos = 0
                    return
                if _uid_scalar_gt(uids[-1], self._hi):
                    end = int(np.searchsorted(uids, self._hi, side="right"))
            if start >= end:
                continue
            if end < len(uids):
                b = b.slice(0, end)
                uids = uids[:end]
            self._batch = b
            self._uids = uids
            self._pos = start
            return

    def has_more(self) -> bool:
        return self._batch is not None and self._pos < self._batch.num_rows

    def take_upto(self, hi_uid: np.void) -> Iterator[tuple[pa.RecordBatch, np.ndarray]]:
        """Yield `(batch_slice, uid_slice)` for all remaining rows in window with `uid <= hi_uid`."""
        while self._batch is not None:
            n = self._batch.num_rows
            if self._pos >= n:
                self._advance()
                continue
            if _uid_scalar_gt(self._uids[self._pos], hi_uid):
                break
            if _uid_scalar_le(self._uids[n - 1], hi_uid):
                p = self._pos
                b_slice = self._batch if p == 0 else self._batch.slice(p, n - p)
                u_slice = self._uids if p == 0 else self._uids[p:n]
                self._advance()
                yield b_slice, u_slice
            else:
                end = int(np.searchsorted(self._uids, hi_uid, side="right"))
                if end > self._pos:
                    p = self._pos
                    b_slice = self._batch.slice(p, end - p)
                    u_slice = self._uids[p:end]
                    self._pos = end
                    yield b_slice, u_slice
                break


def _list_nr50_row_groups(mc_files: list[Path]) -> list[tuple[Path, int]]:
    out: list[tuple[Path, int]] = []
    for f in mc_files:
        pf = pq.ParquetFile(f)
        meta = pf.metadata
        lvl_idx = pf.schema_arrow.names.index("level")
        for rg_i in range(meta.num_row_groups):
            st = meta.row_group(rg_i).column(lvl_idx).statistics
            if st is not None and st.has_min_max:
                mn = st.min.decode("utf-8") if isinstance(st.min, bytes) else str(st.min)
                mx = st.max.decode("utf-8") if isinstance(st.max, bytes) else str(st.max)
                if mn > "nr50" or mx < "nr50":
                    continue
            out.append((f, rg_i))
    return out


def _read_nr50_rg(f: Path, rg_i: int) -> tuple[np.ndarray, np.ndarray]:
    pf = pq.ParquetFile(f)
    meta = pf.metadata
    lvl_idx = pf.schema_arrow.names.index("level")
    st = meta.row_group(rg_i).column(lvl_idx).statistics
    all_nr50 = (
        st is not None
        and st.has_min_max
        and (st.min == b"nr50" or st.min == "nr50")
        and (st.max == b"nr50" or st.max == "nr50")
    )
    if all_nr50:
        tbl = pf.read_row_group(rg_i, columns=["member_uid", "cluster_uid"])
    else:
        tbl = pf.read_row_group(rg_i, columns=["level", "member_uid", "cluster_uid"])
        if tbl.num_rows:
            tbl = tbl.filter(pc.equal(tbl.column("level"), "nr50"))
    if tbl.num_rows == 0:
        return np.empty(0, dtype=UID_DTYPE), np.empty(0, dtype=UID_DTYPE)
    return _uid_array(tbl.column("member_uid")), _uid_array(tbl.column("cluster_uid"))


def _iter_nr50_chunks_from_rgs(
    rgs: list[tuple[Path, int]],
    target_rows: int = SLICE_MC_ROWS,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    buf_m: list[np.ndarray] = []
    buf_c: list[np.ndarray] = []
    buffered = 0
    for f, rg_i in rgs:
        m, c = _read_nr50_rg(f, rg_i)
        if len(m) == 0:
            continue
        buf_m.append(m)
        buf_c.append(c)
        buffered += len(m)
        if buffered >= target_rows:
            yield np.concatenate(buf_m), np.concatenate(buf_c)
            buf_m.clear()
            buf_c.clear()
            buffered = 0
    if buffered > 0:
        yield np.concatenate(buf_m), np.concatenate(buf_c)


# ---------------------------------------------------------------------------------------------------------------
# Collector curve & marginal yield accumulator (exact integer counts across buckets)
# ---------------------------------------------------------------------------------------------------------------


@dataclass
class _CollectorAcc:
    n_stages: int
    earlier_prefix_len: int | None
    inc_by_stage: np.ndarray
    ge1_by_stage: np.ndarray
    ge2_by_stage: np.ndarray
    ge3_by_stage: np.ndarray
    ge5_by_stage: np.ndarray
    brand_new_ge3: np.ndarray
    promoted_ge3: np.ndarray
    brand_new_ge5: np.ndarray
    promoted_ge5: np.ndarray
    earlier_brand_new_ge3: int = 0
    earlier_promoted_ge3: int = 0
    earlier_brand_new_ge5: int = 0
    earlier_promoted_ge5: int = 0

    @classmethod
    def zeros(cls, n_stages: int, earlier_prefix_len: int | None = None) -> "_CollectorAcc":
        z = lambda: np.zeros(n_stages, dtype=np.int64)
        return cls(
            n_stages=n_stages,
            earlier_prefix_len=earlier_prefix_len,
            inc_by_stage=z(),
            ge1_by_stage=z(),
            ge2_by_stage=z(),
            ge3_by_stage=z(),
            ge5_by_stage=z(),
            brand_new_ge3=z(),
            promoted_ge3=z(),
            brand_new_ge5=z(),
            promoted_ge5=z(),
        )

    def iadd(self, other: "_CollectorAcc") -> None:
        self.inc_by_stage += other.inc_by_stage
        self.ge1_by_stage += other.ge1_by_stage
        self.ge2_by_stage += other.ge2_by_stage
        self.ge3_by_stage += other.ge3_by_stage
        self.ge5_by_stage += other.ge5_by_stage
        self.brand_new_ge3 += other.brand_new_ge3
        self.promoted_ge3 += other.promoted_ge3
        self.brand_new_ge5 += other.brand_new_ge5
        self.promoted_ge5 += other.promoted_ge5
        self.earlier_brand_new_ge3 += other.earlier_brand_new_ge3
        self.earlier_promoted_ge3 += other.earlier_promoted_ge3
        self.earlier_brand_new_ge5 += other.earlier_brand_new_ge5
        self.earlier_promoted_ge5 += other.earlier_promoted_ge5


def _accumulate_item_stages(
    item_new: np.ndarray,
    first_stages: np.ndarray,
    n_stages: int,
    earlier_prefix_len: int | None,
) -> _CollectorAcc:
    acc = _CollectorAcc.zeros(n_stages, earlier_prefix_len)
    n_inc = len(first_stages)
    if n_inc == 0:
        return acc

    item_id = np.cumsum(item_new, dtype=np.uint64) - np.uint64(1)
    packed = (item_id << np.uint64(16)) | first_stages.astype(np.uint64)
    packed.sort()
    sorted_stages = (packed & np.uint64(0xFFFF)).astype(np.int64)
    sorted_items = packed >> np.uint64(16)

    new_mask = np.empty(n_inc, dtype=bool)
    new_mask[0] = True
    new_mask[1:] = sorted_items[1:] != sorted_items[:-1]
    starts = np.flatnonzero(new_mask)
    counts = np.diff(np.r_[starts, n_inc])

    acc.inc_by_stage += np.bincount(sorted_stages, minlength=n_stages)
    j1 = sorted_stages[starts]
    acc.ge1_by_stage += np.bincount(j1, minlength=n_stages)

    m2 = counts >= 2
    if np.any(m2):
        acc.ge2_by_stage += np.bincount(sorted_stages[starts[m2] + 1], minlength=n_stages)

    for k, ge_arr, bn_arr, pr_arr, attr_bn, attr_pr in (
        (3, acc.ge3_by_stage, acc.brand_new_ge3, acc.promoted_ge3, "earlier_brand_new_ge3", "earlier_promoted_ge3"),
        (5, acc.ge5_by_stage, acc.brand_new_ge5, acc.promoted_ge5, "earlier_brand_new_ge5", "earlier_promoted_ge5"),
    ):
        mk = counts >= k
        if not np.any(mk):
            continue
        j1_k = j1[mk]
        jk = sorted_stages[starts[mk] + (k - 1)]
        ge_arr += np.bincount(jk, minlength=n_stages)
        bn_mask = j1_k == jk
        if np.any(bn_mask):
            bn_arr += np.bincount(jk[bn_mask], minlength=n_stages)
        if np.any(~bn_mask):
            pr_arr += np.bincount(jk[~bn_mask], minlength=n_stages)
        if earlier_prefix_len is not None:
            m_early = earlier_prefix_len - 1
            crossed = jk > m_early
            if np.any(crossed):
                bn_early = crossed & (j1_k > m_early)
                pr_early = crossed & (j1_k <= m_early)
                setattr(acc, attr_bn, int(np.count_nonzero(bn_early)))
                setattr(acc, attr_pr, int(np.count_nonzero(pr_early)))

    return acc


def _coverage_and_chao1(s_obs: int, f1: int, f2: int, ge3: int, ge5: int, n_inc: int) -> dict:
    c_hat = 1.0 if n_inc == 0 else 1.0 - (f1 / n_inc)
    chao1 = float(s_obs) + (f1 * (f1 - 1)) / (2.0 * (f2 + 1))
    return {
        "s_obs": int(s_obs),
        "f1": int(f1),
        "f2": int(f2),
        "ge3": int(ge3),
        "ge5": int(ge5),
        "n_incidences": int(n_inc),
        "c_hat": c_hat,
        "chao1": chao1,
    }


def _pair_dict_from_stats(st: dict) -> dict:
    return {
        "pairs_ge1": st["s_obs"],
        "pairs_f1": st["f1"],
        "pairs_f2": st["f2"],
        "pairs_ge3": st["ge3"],
        "pairs_ge5": st["ge5"],
        "n_incidences": st["n_incidences"],
        "c_hat_pairs": st["c_hat"],
        "chao1_pairs": st["chao1"],
    }


def _cumulative_series(acc: _CollectorAcc) -> list[dict]:
    cum_inc = np.cumsum(acc.inc_by_stage)
    cum_ge1 = np.cumsum(acc.ge1_by_stage)
    cum_ge2 = np.cumsum(acc.ge2_by_stage)
    cum_ge3 = np.cumsum(acc.ge3_by_stage)
    cum_ge5 = np.cumsum(acc.ge5_by_stage)
    out = []
    for m in range(acc.n_stages):
        s_obs = int(cum_ge1[m])
        ge2 = int(cum_ge2[m])
        ge3 = int(cum_ge3[m])
        ge5 = int(cum_ge5[m])
        f1 = s_obs - ge2
        f2 = ge2 - ge3
        out.append(_coverage_and_chao1(s_obs, f1, f2, ge3, ge5, int(cum_inc[m])))
    return out


# ---------------------------------------------------------------------------------------------------------------
# Bucket scatter & multi-process worker functions
# ---------------------------------------------------------------------------------------------------------------


def _choose_buckets(total_rows: int, explicit: int | None) -> int:
    if explicit is not None and explicit > 0:
        b = 1
        while b < explicit:
            b <<= 1
        return min(b, 256)
    if total_rows <= 500_000:
        return 1
    b = 8
    while b < 256 and b * TARGET_BUCKET_ROWS < total_rows:
        b <<= 1
    return b


def _scatter_rows(handles: list, bucket_idx: np.ndarray, words3: np.ndarray, n_buckets: int) -> None:
    if len(bucket_idx) == 0:
        return
    if n_buckets == 1:
        handles[0].write(np.ascontiguousarray(words3).tobytes())
        return
    order = np.argsort(bucket_idx, kind="quicksort")
    sorted_b = bucket_idx[order]
    sorted_w = np.ascontiguousarray(words3[order])
    cuts = np.searchsorted(sorted_b, np.arange(n_buckets + 1, dtype=sorted_b.dtype))
    for b in np.flatnonzero(np.diff(cuts)):
        lo, hi = int(cuts[b]), int(cuts[b + 1])
        handles[int(b)].write(sorted_w[lo:hi].tobytes())


def _read_concat_and_unlink(paths: list[Path]) -> np.ndarray:
    chunks: list[np.ndarray] = []
    try:
        for p in paths:
            if p.is_file() and p.stat().st_size > 0:
                chunks.append(np.fromfile(p, dtype="<u8").reshape(-1, 3))
    finally:
        for p in paths:
            p.unlink(missing_ok=True)
    if not chunks:
        return np.empty((0, 3), dtype="<u8")
    if len(chunks) == 1:
        return chunks[0]
    out = np.concatenate(chunks, axis=0)
    chunks.clear()
    return out


def _sort_3u64_fast(raw: np.ndarray) -> np.ndarray:
    """Sort `N x 3` `<u8` array ascending by `(col0, col1, col2)` via single-pass primary sort + tie fixup."""
    n = len(raw)
    if n <= 1:
        return raw
    # Fast single-pass primitive uint64 introsort on col0
    o0 = np.argsort(raw[:, 0], kind="quicksort")
    s = raw[o0]
    del o0
    same0 = s[1:, 0] == s[:-1, 0]
    if not np.any(same0):
        return s
    # Check whether col1 ever differs inside equal-col0 runs
    has_col1_tie = bool(np.any(same0 & (s[1:, 1] != s[:-1, 1])))
    if not has_col1_tie:
        # Common ORF-cluster path: col0 uniquely identifies each cluster in the bucket!
        # Assign local cluster rank (<= N <= ~3M < 2^23) and pack `(rank << 40) | col2` when col2 < 2^40
        max_w2 = int(s[:, 2].max())
        new_c = np.empty(n, dtype=bool)
        new_c[0] = True
        new_c[1:] = ~same0
        n_items = int(np.count_nonzero(new_c))
        if max_w2 < (1 << 40) and n_items < (1 << 24):
            rank = np.cumsum(new_c, dtype=np.uint64) - np.uint64(1)
            key = (rank << np.uint64(40)) | s[:, 2]
            o2 = np.argsort(key, kind="quicksort")
            return s[o2]
    # General fallback (including pairs where one t1 pairs with multiple t2):
    # Check if `(t1, t2)` share upper 34 bits so `(local1 << 32) | (local2 << 2) | scode` fits in 64 bits
    c0_xor = int(s[:, 0].min()) ^ int(s[:, 0].max())
    c1_xor = int(s[:, 1].min()) ^ int(s[:, 1].max())
    if (c0_xor >> 30) == 0 and (c1_xor >> 30) == 0:
        k64 = ((s[:, 0] & np.uint64(0x3FFFFFFF)) << np.uint64(34)) | (
            (s[:, 1] & np.uint64(0x3FFFFFFF)) << np.uint64(4)
        ) | (s[:, 2] >> np.uint64(48))
        o_pair = np.argsort(k64, kind="quicksort")
        s = s[o_pair]
        del k64, o_pair
        # Now `(col0, col1, scode)` are sorted; sort `(gid, stage)` within equal pairs in one 64-bit pass!
        scode = s[:, 2] >> np.uint64(48)
        new_p = np.empty(n, dtype=bool)
        new_p[0] = True
        new_p[1:] = (s[1:, 0] != s[:-1, 0]) | (s[1:, 1] != s[:-1, 1]) | (scode[1:] != scode[:-1])
        n_p = int(np.count_nonzero(new_p))
        low48 = s[:, 2] & np.uint64(0xFFFFFFFFFFFF)
        if int(low48.max()) < (1 << 40) and n_p < (1 << 24):
            prank = np.cumsum(new_p, dtype=np.uint64) - np.uint64(1)
            k2 = (prank << np.uint64(40)) | low48
            return s[np.argsort(k2, kind="quicksort")]
    order = np.lexsort((s[:, 2], s[:, 1], s[:, 0]))
    return s[order]


def _reduce_orf_bucket_job(args: tuple[list[Path], int, int | None]) -> tuple[_CollectorAcc, _CollectorAcc, int]:
    pa.set_cpu_count(1)
    pa.set_io_thread_count(1)
    paths, n_stages, earlier_prefix_len = args
    raw = _read_concat_and_unlink(paths)
    if len(raw) == 0:
        return (
            _CollectorAcc.zeros(n_stages, earlier_prefix_len),
            _CollectorAcc.zeros(n_stages, earlier_prefix_len),
            read_vmhwm_kib(),
        )
    s = _sort_3u64_fast(raw)
    del raw
    n = len(s)
    w0, w1, w2 = s[:, 0], s[:, 1], s[:, 2]
    gid = w2 >> np.uint64(17)

    new_c = np.empty(n, dtype=bool)
    new_c[0] = True
    new_c[1:] = (w0[1:] != w0[:-1]) | (w1[1:] != w1[:-1])

    new_cg = np.empty(n, dtype=bool)
    new_cg[0] = True
    new_cg[1:] = new_c[1:] | (gid[1:] != gid[:-1])
    cg_starts = np.flatnonzero(new_cg)

    stage_rows = (w2 & np.uint64(0xFFFF)).astype(np.int64)
    first_stage_all = np.minimum.reduceat(stage_rows, cg_starts)
    c_new_all = new_c[cg_starts]
    acc_all = _accumulate_item_stages(c_new_all, first_stage_all, n_stages, earlier_prefix_len)

    has_comp = ((w2[cg_starts] >> np.uint64(16)) & np.uint64(1)) == 0
    if not np.any(has_comp):
        return acc_all, _CollectorAcc.zeros(n_stages, earlier_prefix_len), read_vmhwm_kib()

    comp_idx = cg_starts[has_comp]
    cw0, cw1 = w0[comp_idx], w1[comp_idx]
    c_new_comp = np.empty(len(comp_idx), dtype=bool)
    c_new_comp[0] = True
    c_new_comp[1:] = (cw0[1:] != cw0[:-1]) | (cw1[1:] != cw1[:-1])
    first_stage_comp = stage_rows[comp_idx]
    acc_comp = _accumulate_item_stages(c_new_comp, first_stage_comp, n_stages, earlier_prefix_len)
    return acc_all, acc_comp, read_vmhwm_kib()


def _reduce_pair_bucket_job(args: tuple[list[Path], int, int | None]) -> tuple[_CollectorAcc, int]:
    pa.set_cpu_count(1)
    pa.set_io_thread_count(1)
    paths, n_stages, earlier_prefix_len = args
    raw = _read_concat_and_unlink(paths)
    if len(raw) == 0:
        return _CollectorAcc.zeros(n_stages, earlier_prefix_len), read_vmhwm_kib()

    s = _sort_3u64_fast(raw)
    del raw
    n = len(s)
    w0, w1, w2 = s[:, 0], s[:, 1], s[:, 2]
    scode = w2 >> np.uint64(48)
    gid = (w2 >> np.uint64(16)) & np.uint64(0xFFFFFFFF)

    new_p = np.empty(n, dtype=bool)
    new_p[0] = True
    new_p[1:] = (w0[1:] != w0[:-1]) | (w1[1:] != w1[:-1]) | (scode[1:] != scode[:-1])

    new_pg = np.empty(n, dtype=bool)
    new_pg[0] = True
    new_pg[1:] = new_p[1:] | (gid[1:] != gid[:-1])
    pg_starts = np.flatnonzero(new_pg)

    first_stage = (w2[pg_starts] & np.uint64(0xFFFF)).astype(np.int64)
    p_new_at_pg = new_p[pg_starts]
    return _accumulate_item_stages(p_new_at_pg, first_stage, n_stages, earlier_prefix_len), read_vmhwm_kib()


def _scatter_orf_worker(task: dict) -> dict:
    pa.set_cpu_count(1)
    pa.set_io_thread_count(1)
    w_id: int = task["worker_id"]
    stages: list[str] = task["stages"]
    n_stages = len(stages)
    mc_rgs: list[tuple[Path, int]] = task["mc_rgs"]
    po_by_stage: dict[str, list[Path]] = task["po_by_stage"]
    lo_pu = _pair_to_uid(task["lo_pu"])
    hi_pu = _pair_to_uid(task["hi_pu"])
    is_last: bool = task["is_last"]
    n_buckets: int = task["n_buckets"]
    work_dir: Path = task["work_dir"]
    n_cdim: int = task["n_cdim"]

    if n_cdim > 0:
        cdim_cu = np.memmap(work_dir / "cdim_cu.bin", dtype=UID_DTYPE, mode="r", shape=(n_cdim,))
        cdim_hi = np.memmap(work_dir / "cdim_hi.bin", dtype="<u8", mode="r", shape=(n_cdim,))
        cdim_lo = np.memmap(work_dir / "cdim_lo.bin", dtype="<u8", mode="r", shape=(n_cdim,))
        cdim_gid = np.memmap(work_dir / "cdim_gid.bin", dtype="<u4", mode="r", shape=(n_cdim,))
    else:
        cdim_cu = np.empty(0, dtype=UID_DTYPE)
        cdim_hi = cdim_lo = np.empty(0, dtype="<u8")
        cdim_gid = np.empty(0, dtype="<u4")

    stage_occ_all = np.zeros(n_stages, dtype=np.int64)
    stage_occ_comp = np.zeros(n_stages, dtype=np.int64)
    nr50_members = 0

    po_streams = [
        _SortedParquetStream(
            po_by_stage[s],
            ["protein_uid", "contig_uid", "dist_to_contig_end_bp"],
            "protein_uid",
            lo_uid_exclusive=lo_pu,
            hi_uid_inclusive=None if is_last else hi_pu,
        )
        for s in stages
    ]

    mask_u64 = np.uint64(n_buckets - 1)
    c1_u64 = np.uint64(0x9E3779B97F4A7C15)
    c2_u64 = np.uint64(0xBF58476D1CE4E5B9)
    out_paths = [work_dir / f"orf-w{w_id:02d}-b{b:03d}.bin" for b in range(n_buckets)]

    with contextlib.ExitStack() as stack:
        handles = [stack.enter_context(open(p, "wb")) for p in out_paths]
        for mc_pu, mc_cu in _iter_nr50_chunks_from_rgs(mc_rgs):
            nr50_members += len(mc_pu)
            mc_hi = mc_pu["hi"]
            mc_lo = mc_pu["lo"]
            chunk_hi = mc_pu[-1]
            for s_idx, stream in enumerate(po_streams):
                for b_slice, pu_slice in stream.take_upto(chunk_hi):
                    n_rows = b_slice.num_rows
                    stage_occ_all[s_idx] += n_rows
                    dist_np = b_slice.column("dist_to_contig_end_bp").to_numpy(zero_copy_only=False)
                    comp_mask = dist_np > 0
                    stage_occ_comp[s_idx] += int(np.count_nonzero(comp_mask))

                    mc_pos = _lookup_sorted_uids(
                        mc_pu, pu_slice, f"protein_occurrence[{stages[s_idx]}]->member_cluster",
                        sorted_hi=mc_hi, sorted_lo=mc_lo,
                    )
                    cu_arr = _uid_array(b_slice.column("contig_uid"))
                    cd_pos = _lookup_sorted_uids(
                        cdim_cu, cu_arr, f"protein_occurrence[{stages[s_idx]}]->contig_dim",
                        sorted_hi=cdim_hi, sorted_lo=cdim_lo,
                    )

                    cl_hi = mc_cu["hi"][mc_pos]
                    cl_lo = mc_cu["lo"][mc_pos]
                    gid = cdim_gid[cd_pos].astype(np.uint64)
                    not_comp = (~comp_mask).astype(np.uint64)
                    w2 = (gid << np.uint64(17)) | (not_comp << np.uint64(16)) | np.uint64(s_idx)

                    words = np.empty((n_rows, 3), dtype="<u8")
                    words[:, 0] = cl_hi
                    words[:, 1] = cl_lo
                    words[:, 2] = w2

                    with np.errstate(over="ignore"):
                        mix = (cl_hi ^ (cl_lo * c1_u64)) * c2_u64
                        mix ^= mix >> np.uint64(30)
                    _scatter_rows(handles, mix & mask_u64, words, n_buckets)

        for s_idx, stream in enumerate(po_streams):
            if stream.has_more():
                raise ValueError(f"protein_occurrence of stage {stages[s_idx]} has protein_uid missing from member_cluster")

    return {
        "nr50_members": nr50_members,
        "stage_occ_all": stage_occ_all,
        "stage_occ_comp": stage_occ_comp,
        "vmhwm_kib": read_vmhwm_kib(),
    }


def _scatter_pair_worker(task: dict) -> int:
    pa.set_cpu_count(1)
    pa.set_io_thread_count(1)
    w_id: int = task["worker_id"]
    stages: list[str] = task["stages"]
    tt_rgs: list[tuple[Path, int]] = task["tt_rgs"]
    ct_by_stage: dict[str, list[Path]] = task["ct_by_stage"]
    lo_cu = _pair_to_uid(task["lo_cu"])
    hi_cu = _pair_to_uid(task["hi_cu"])
    is_last: bool = task["is_last"]
    n_buckets: int = task["n_buckets"]
    work_dir: Path = task["work_dir"]
    n_cdim: int = task["n_cdim"]

    if n_cdim > 0:
        cdim_cu = np.memmap(work_dir / "cdim_cu.bin", dtype=UID_DTYPE, mode="r", shape=(n_cdim,))
        cdim_hi = np.memmap(work_dir / "cdim_hi.bin", dtype="<u8", mode="r", shape=(n_cdim,))
        cdim_lo = np.memmap(work_dir / "cdim_lo.bin", dtype="<u8", mode="r", shape=(n_cdim,))
        cdim_gid = np.memmap(work_dir / "cdim_gid.bin", dtype="<u4", mode="r", shape=(n_cdim,))
    else:
        cdim_cu = np.empty(0, dtype=UID_DTYPE)
        cdim_hi = cdim_lo = np.empty(0, dtype="<u8")
        cdim_gid = np.empty(0, dtype="<u4")

    ct_streams = [
        _SortedParquetStream(
            ct_by_stage[s],
            ["contig_uid", "is_circular", "n_slots", "strands"],
            "contig_uid",
            lo_uid_exclusive=lo_cu,
            hi_uid_inclusive=None if is_last else hi_cu,
            batch_size=TAPE_BATCH_ROWS,
        )
        for s in stages
    ]

    mask_u64 = np.uint64(n_buckets - 1)
    c1_u64 = np.uint64(0x9E3779B97F4A7C15)
    c2_u64 = np.uint64(0xBF58476D1CE4E5B9)
    out_paths = [work_dir / f"pair-w{w_id:02d}-b{b:03d}.bin" for b in range(n_buckets)]

    with contextlib.ExitStack() as stack:
        handles = [stack.enter_context(open(p, "wb")) for p in out_paths]
        for tt_file, rg_i in tt_rgs:
            pf_tt = pq.ParquetFile(tt_file)
            for tt_batch in pf_tt.iter_batches(
                batch_size=TAPE_BATCH_ROWS, row_groups=[rg_i], columns=["contig_uid", "nr50_tokens"]
            ):
                if not tt_batch.num_rows:
                    continue
                tt_cu = _uid_array(tt_batch.column("contig_uid"))
                tt_hi = tt_cu["hi"]
                tt_lo = tt_cu["lo"]
                chunk_hi = tt_cu[-1]
                nr50_list = tt_batch.column("nr50_tokens")
                tt_off = nr50_list.offsets.to_numpy(zero_copy_only=False).astype(np.int64)
                tok_vals = nr50_list.values
                if tok_vals.null_count > 0:
                    tok_valid = tok_vals.is_valid().to_numpy(zero_copy_only=False)
                    tok_np = pc.fill_null(tok_vals, 0).to_numpy(zero_copy_only=False).astype(np.uint64)
                else:
                    tok_valid = None
                    tok_np = tok_vals.to_numpy(zero_copy_only=False).astype(np.uint64)

                for s_idx, stream in enumerate(ct_streams):
                    for ct_slice, ct_cu in stream.take_upto(chunk_hi):
                        pos_tt = _lookup_sorted_uids(
                            tt_cu, ct_cu, f"contig_tapes[{stages[s_idx]}]->tape_tokens",
                            sorted_hi=tt_hi, sorted_lo=tt_lo,
                        )
                        pos_cd = _lookup_sorted_uids(
                            cdim_cu, ct_cu, f"contig_tapes[{stages[s_idx]}]->contig_dim",
                            sorted_hi=cdim_hi, sorted_lo=cdim_lo,
                        )
                        gid_c = cdim_gid[pos_cd].astype(np.uint64)

                        strands_list = ct_slice.column("strands")
                        ct_off = strands_list.offsets.to_numpy(zero_copy_only=False).astype(np.int64)
                        n_slots = ct_off[1:] - ct_off[:-1]
                        tt_n_slots = tt_off[pos_tt + 1] - tt_off[pos_tt]
                        if not np.array_equal(n_slots, tt_n_slots):
                            raise ValueError(f"slot count mismatch between contig_tapes and tape_tokens in {stages[s_idx]}")

                        n_adj = np.maximum(0, n_slots - 1)
                        m_lin = int(n_adj.sum())
                        circ_np = ct_slice.column("is_circular").to_numpy(zero_copy_only=False)
                        c_wrap = np.flatnonzero(circ_np & (n_slots >= 2))
                        m_wrap = len(c_wrap)
                        if m_lin + m_wrap == 0:
                            continue

                        parts_ct_a: list[np.ndarray] = []
                        parts_ct_b: list[np.ndarray] = []
                        parts_tt_a: list[np.ndarray] = []
                        parts_tt_b: list[np.ndarray] = []
                        parts_gid: list[np.ndarray] = []

                        if m_lin > 0:
                            c_idx = np.repeat(np.arange(len(n_slots), dtype=np.int64), n_adj)
                            p_starts = np.zeros(len(n_adj) + 1, dtype=np.int64)
                            np.cumsum(n_adj, out=p_starts[1:])
                            r = np.arange(m_lin, dtype=np.int64) - np.repeat(p_starts[:-1], n_adj)
                            ca = ct_off[c_idx] + r
                            ta = tt_off[pos_tt[c_idx]] + r
                            parts_ct_a.append(ca)
                            parts_ct_b.append(ca + 1)
                            parts_tt_a.append(ta)
                            parts_tt_b.append(ta + 1)
                            parts_gid.append(gid_c[c_idx])

                        if m_wrap > 0:
                            parts_ct_a.append(ct_off[c_wrap + 1] - 1)
                            parts_ct_b.append(ct_off[c_wrap])
                            parts_tt_a.append(tt_off[pos_tt[c_wrap] + 1] - 1)
                            parts_tt_b.append(tt_off[pos_tt[c_wrap]])
                            parts_gid.append(gid_c[c_wrap])

                        idx_ct_a = parts_ct_a[0] if len(parts_ct_a) == 1 else np.concatenate(parts_ct_a)
                        idx_ct_b = parts_ct_b[0] if len(parts_ct_b) == 1 else np.concatenate(parts_ct_b)
                        idx_tt_a = parts_tt_a[0] if len(parts_tt_a) == 1 else np.concatenate(parts_tt_a)
                        idx_tt_b = parts_tt_b[0] if len(parts_tt_b) == 1 else np.concatenate(parts_tt_b)
                        p_gid = parts_gid[0] if len(parts_gid) == 1 else np.concatenate(parts_gid)

                        if tok_valid is not None:
                            keep = tok_valid[idx_tt_a] & tok_valid[idx_tt_b]
                            if not np.any(keep):
                                continue
                            idx_ct_a = idx_ct_a[keep]
                            idx_ct_b = idx_ct_b[keep]
                            idx_tt_a = idx_tt_a[keep]
                            idx_tt_b = idx_tt_b[keep]
                            p_gid = p_gid[keep]

                        strands_vals = strands_list.values.to_numpy(zero_copy_only=False).astype(np.int8)
                        sa = strands_vals[idx_ct_a]
                        sb = strands_vals[idx_ct_b]
                        ta_v = tok_np[idx_tt_a]
                        tb_v = tok_np[idx_tt_b]

                        nsa = -sa
                        nsb = -sb
                        swap = (tb_v < ta_v) | ((tb_v == ta_v) & (nsb < sa))
                        t1 = np.where(swap, tb_v, ta_v)
                        t2 = np.where(swap, ta_v, tb_v)
                        s1 = np.where(swap, nsb, sa)
                        s2 = np.where(swap, nsa, sb)
                        scode = ((s1 > 0).astype(np.uint64) << np.uint64(1)) | (s2 > 0).astype(np.uint64)
                        w2_p = (scode << np.uint64(48)) | (p_gid << np.uint64(16)) | np.uint64(s_idx)

                        n_emit = len(t1)
                        words_p = np.empty((n_emit, 3), dtype="<u8")
                        words_p[:, 0] = t1
                        words_p[:, 1] = t2
                        words_p[:, 2] = w2_p

                        with np.errstate(over="ignore"):
                            pmix = (t1 * c1_u64) ^ (t2 * c2_u64) ^ scode
                            pmix ^= pmix >> np.uint64(30)
                        _scatter_rows(handles, pmix & mask_u64, words_p, n_buckets)

        for s_idx, stream in enumerate(ct_streams):
            if stream.has_more():
                raise ValueError(f"contig_tapes of stage {stages[s_idx]} has contig_uid missing from tape_tokens")

    return read_vmhwm_kib()


# ---------------------------------------------------------------------------------------------------------------
# Core `measure` implementation
# ---------------------------------------------------------------------------------------------------------------


def _check_not_inside(child: Path, parent: Path, label: str) -> None:
    c = child.resolve()
    p = parent.resolve()
    if c == p or p in c.parents:
        raise ValueError(f"{label} ({child}) must not lie inside the data root ({parent})")


def _manifest_files(root: Path, files: list[dict], tbl_name: str, **parts: str) -> list[Path]:
    matched = [
        root / f["path"]
        for f in files
        if f.get("table") == tbl_name and all(f.get("partitions", {}).get(k) == v for k, v in parts.items())
    ]
    for p in matched:
        if not p.is_file():
            raise FileNotFoundError(f"manifested {tbl_name} file missing on disk: {p}")
    return matched


def _resolve_tape_and_gcv(root: Path, rel: str, files: list[dict]) -> tuple[str, str]:
    tvs = sorted({
        f["partitions"]["tape_version"]
        for f in files
        if f.get("table") == "tape_tokens" and f.get("partitions", {}).get("release") == rel
    })
    if len(tvs) != 1:
        raise ValueError(f"expected exactly 1 tape_version in tape_tokens for release {rel}, found {tvs}")
    tv = tvs[0]

    gcvs = sorted({
        f["partitions"]["gene_call_version"]
        for f in files
        if f.get("table") == "protein_occurrence" and "gene_call_version" in f.get("partitions", {})
    })
    if len(gcvs) == 1:
        return tv, gcvs[0]
    if not gcvs:
        return tv, ""
    tv_defs = _manifest_files(root, files, "tape_version_def", tape_version=tv)
    if tv_defs:
        prods = set(pq.read_table(tv_defs[0], columns=["producer"]).column("producer").to_pylist())
        inter = [g for g in gcvs if g in prods]
        if len(inter) == 1:
            return tv, inter[0]
    raise ValueError(f"cannot resolve unique gene_call_version among {gcvs} for tape_version {tv}")


def _write_contig_dim_memmap(cdim_files: list[Path], work_dir: Path) -> tuple[int, int]:
    """Build read-only memmappable `cdim_cu`, `cdim_hi`, `cdim_lo`, `cdim_gid` files; return `(n_cdim, n_genomes)`."""
    cu_parts: list[np.ndarray] = []
    gu_parts: list[np.ndarray] = []
    for f in cdim_files:
        pf = pq.ParquetFile(f)
        for b in pf.iter_batches(batch_size=STREAM_BATCH_ROWS, columns=["contig_uid", "genome_uid"]):
            if b.num_rows:
                cu_parts.append(_uid_array(b.column("contig_uid")))
                gu_parts.append(_uid_array(b.column("genome_uid")))
    if not cu_parts:
        return 0, 0
    cu = np.concatenate(cu_parts)
    gu = np.concatenate(gu_parts)
    del cu_parts, gu_parts
    n_cdim = len(cu)

    g_order = np.lexsort((gu["lo"], gu["hi"]))
    s_gu = gu[g_order]
    new_g = np.empty(len(s_gu), dtype=bool)
    new_g[0] = True
    new_g[1:] = ~_uid_eq(s_gu[1:], s_gu[:-1])
    n_genomes = int(np.count_nonzero(new_g))
    gid_sorted = np.cumsum(new_g, dtype=np.uint32) - np.uint32(1)
    gid = np.empty(n_cdim, dtype=np.uint32)
    gid[g_order] = gid_sorted
    del gu, g_order, s_gu, new_g, gid_sorted

    cu.tofile(work_dir / "cdim_cu.bin")
    np.ascontiguousarray(cu["hi"]).tofile(work_dir / "cdim_hi.bin")
    np.ascontiguousarray(cu["lo"]).tofile(work_dir / "cdim_lo.bin")
    gid.tofile(work_dir / "cdim_gid.bin")
    del cu, gid
    return n_cdim, n_genomes


def _stage_stored_bases(ci_files: list[Path]) -> tuple[int, int]:
    contigs = 0
    bases = 0
    for f in ci_files:
        pf = pq.ParquetFile(f)
        for b in pf.iter_batches(batch_size=STREAM_BATCH_ROWS, columns=["contig_len", "is_default_variant"]):
            if not b.num_rows:
                continue
            sub = b.filter(b.column("is_default_variant"))
            if sub.num_rows:
                contigs += sub.num_rows
                bases += int(pc.sum(sub.column("contig_len")).as_py() or 0)
    return contigs, bases


def _partition_rgs(rgs: list[tuple[Path, int]], n_parts: int) -> list[list[tuple[Path, int]]]:
    if not rgs:
        return []
    k = max(1, min(n_parts, len(rgs)))
    out: list[list[tuple[Path, int]]] = []
    for i in range(k):
        lo = (i * len(rgs)) // k
        hi = ((i + 1) * len(rgs)) // k
        if hi > lo:
            out.append(rgs[lo:hi])
    return out


def _last_nr50_member_uid(rgs: list[tuple[Path, int]]) -> tuple[int, int] | None:
    for f, rg_i in reversed(rgs):
        m, _ = _read_nr50_rg(f, rg_i)
        if len(m) > 0:
            return _uid_to_pair(m[-1])
    return None


def _last_tt_contig_uid(rgs: list[tuple[Path, int]]) -> tuple[int, int] | None:
    for f, rg_i in reversed(rgs):
        pf = pq.ParquetFile(f)
        meta = pf.metadata
        cu_idx = pf.schema_arrow.names.index("contig_uid")
        st = meta.row_group(rg_i).column(cu_idx).statistics
        if st is not None and st.has_min_max and isinstance(st.max, bytes) and len(st.max) == 16:
            return _uid_to_pair(_bytes_to_uid_scalar(st.max))
        tbl = pf.read_row_group(rg_i, columns=["contig_uid"])
        if tbl.num_rows > 0:
            return _uid_to_pair(_uid_array(tbl.column("contig_uid"))[-1])
    return None


def format_saturation_txt(doc: dict) -> str:
    tot = doc["totals"]
    comp = tot["complete_orf"]
    all_o = tot["all_orf"]
    pairs = tot["syntenic_pairs"]
    lines = [
        f"Release saturation summary: {doc['release']} (earlier: {doc['earlier'] or 'none'})",
        f"Stages ({len(doc['stages'])}): {','.join(doc['stages'])}",
        f"Stored sequence: {tot['stored_bases']:,} bp ({tot['stored_gbp']:.3f} Gbp) across {tot['stored_contigs']:,} contigs "
        f"({tot['distinct_genomes_contig_dim']:,} contig_dim genomes)",
        f"ORF occurrences: {tot['protein_occurrences_all']:,} total, {tot['protein_occurrences_complete']:,} complete, "
        f"{tot['protein_occurrences_edge_truncated']:,} linear-edge-truncated",
        f"NR50 all ORFs:      S_obs={all_o['s_obs']:,}  f1={all_o['f1']:,}  f2={all_o['f2']:,}  "
        f"k>=3={all_o['ge3']:,}  k>=5={all_o['ge5']:,}  C_hat={all_o['c_hat']:.6f}  Chao1={all_o['chao1']:,.1f}",
        f"NR50 complete ORFs: S_obs={comp['s_obs']:,}  f1={comp['f1']:,}  f2={comp['f2']:,}  "
        f"k>=3={comp['ge3']:,}  k>=5={comp['ge5']:,}  C_hat={comp['c_hat']:.6f}  Chao1={comp['chao1']:,.1f}",
        f"Syntenic NR50 pairs: S_obs={pairs['pairs_ge1']:,}  f1={pairs['pairs_f1']:,}  f2={pairs['pairs_f2']:,}  "
        f"k>=3={pairs['pairs_ge3']:,}  k>=5={pairs['pairs_ge5']:,}  C_hat={pairs['c_hat_pairs']:.6f}  "
        f"Chao1={pairs['chao1_pairs']:,.1f}",
    ]
    se = doc.get("since_earlier")
    if se is not None:
        rate3 = f"{se['net_new_ge3_per_gbp']:,.1f}/Gbp" if se["net_new_ge3_per_gbp"] is not None else "n/a"
        prate3 = f"{se['net_new_pairs_ge3_per_gbp']:,.1f}/Gbp" if se["net_new_pairs_ge3_per_gbp"] is not None else "n/a"
        lines.append(
            f"Since {se['earlier']} (new stages {','.join(se['new_stages']) or 'none'}, +{se['new_stored_gbp']:.3f} Gbp): "
            f"+{se['net_new_ge3']:,} k>=3 complete-ORF NR50 ({rate3}: {se['brand_new_ge3']:,} brand-new + "
            f"{se['promoted_to_ge3']:,} promoted); +{se['net_new_pairs_ge3']:,} k>=3 syntenic pairs ({prate3})"
        )
    lines.append("")
    lines.append(
        "stage\tstored_gbp\tcum_gbp\tall_nr50\tcomp_nr50\tcomp_ge3\tcomp_ge5\tc_hat_comp\t"
        "net_ge3_per_gbp\tbrand_new_ge3\tpromoted_ge3\tpairs_ge3\tc_hat_pairs\tnet_pairs_ge3_per_gbp"
    )
    for row in doc["stage_curve"]:
        co = row["complete_orf"]
        ao = row["all_orf"]
        sp = row["syntenic_pairs"]
        mg = row["marginal"]
        r3 = "-" if mg["net_new_ge3_per_gbp"] is None else f"{mg['net_new_ge3_per_gbp']:.1f}"
        pr3 = "-" if mg["net_new_pairs_ge3_per_gbp"] is None else f"{mg['net_new_pairs_ge3_per_gbp']:.1f}"
        lines.append(
            f"{row['stage']}\t{row['stored_gbp']:.3f}\t{row['cum_stored_gbp']:.3f}\t"
            f"{ao['s_obs']}\t{co['s_obs']}\t{co['ge3']}\t{co['ge5']}\t{co['c_hat']:.6f}\t"
            f"{r3}\t{mg['brand_new_ge3']}\t{mg['promoted_to_ge3']}\t"
            f"{sp['pairs_ge3']}\t{sp['c_hat_pairs']:.6f}\t{pr3}"
        )
    return "\n".join(lines) + "\n"


def measure_release(
    *,
    root: Path,
    release: str,
    out: Path,
    tmp: Path,
    earlier: str | None = None,
    threads: int = 4,
    buckets: int | None = None,
    stats_out: Path | None = None,
) -> dict:
    t0 = time.monotonic()
    cpu0 = time.process_time()
    import resource

    rc0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    child_cpu0 = float(rc0.ru_utime + rc0.ru_stime)
    root = Path(root)
    out = Path(out)
    tmp = Path(tmp)
    _check_not_inside(tmp, root, "--tmp")
    _check_not_inside(out, root, "--out")

    threads = max(1, min(4, int(threads)))
    pa.set_cpu_count(1)
    pa.set_io_thread_count(1)

    arch_path = root / "archived" / f"release={release}.json"
    if arch_path.is_file():
        raise FileNotFoundError(f"release {release} is archived ({arch_path}); restore its release= tables before measuring")
    man_path = root / "manifests" / f"release={release}.json"
    if not man_path.is_file():
        raise FileNotFoundError(f"missing manifest: {man_path}")

    man = json.loads(man_path.read_text(encoding="utf-8"))
    man_rel = man.get("release_id") or man.get("release")
    if man_rel and man_rel != release:
        raise ValueError(f"{man_path}: manifest names release {man_rel!r}, expected {release!r}")
    stages = list(man["stages"])
    n_stages = len(stages)
    if n_stages == 0:
        raise ValueError(f"release {release} has 0 stages in {man_path}")
    files = man["files"]

    earlier_stages: list[str] | None = None
    earlier_prefix_len: int | None = None
    if earlier:
        eman_path = root / "manifests" / f"release={earlier}.json"
        if not eman_path.is_file():
            raise FileNotFoundError(f"missing earlier manifest: {eman_path}")
        eman = json.loads(eman_path.read_text(encoding="utf-8"))
        estages = list(eman["stages"])
        if stages[: len(estages)] == estages:
            earlier_stages = estages
            earlier_prefix_len = len(estages)

    tv, gcv = _resolve_tape_and_gcv(root, release, files)
    mc_files = _manifest_files(root, files, "member_cluster", release=release)
    cdim_files = _manifest_files(root, files, "contig_dim", release=release)
    tt_files = _manifest_files(root, files, "tape_tokens", release=release, tape_version=tv)
    if not mc_files or not cdim_files or not tt_files:
        raise FileNotFoundError(
            f"release {release} manifest lacks required release tables "
            f"(member_cluster={len(mc_files)}, contig_dim={len(cdim_files)}, tape_tokens={len(tt_files)})"
        )

    ci_by_stage = {s: _manifest_files(root, files, "contig_index", stage=s) for s in stages}
    po_by_stage = {
        s: (_manifest_files(root, files, "protein_occurrence", stage=s, gene_call_version=gcv) if gcv else [])
        for s in stages
    }
    ct_by_stage = {s: _manifest_files(root, files, "contig_tapes", stage=s, tape_version=tv) for s in stages}

    total_occ_rows = sum(
        int(f.get("rows", 0))
        for f in files
        if f.get("table") == "protein_occurrence"
        and (not gcv or f.get("partitions", {}).get("gene_call_version") == gcv)
    )
    n_buckets = _choose_buckets(total_occ_rows, buckets)

    tmp.mkdir(parents=True, exist_ok=True)
    work_dir = Path(tempfile.mkdtemp(prefix=f"sat-{release}-", dir=tmp))
    worker_peaks_kib: list[int] = []
    use_pool = False
    added_syspath = False
    script_dir = str(Path(__file__).resolve().parent)
    try:
        stage_contigs = np.zeros(n_stages, dtype=np.int64)
        stage_bases = np.zeros(n_stages, dtype=np.int64)
        for idx, s in enumerate(stages):
            c_cnt, b_cnt = _stage_stored_bases(ci_by_stage[s])
            stage_contigs[idx] = c_cnt
            stage_bases[idx] = b_cnt

        n_cdim, n_genomes = _write_contig_dim_memmap(cdim_files, work_dir)

        mc_rgs = _list_nr50_row_groups(mc_files)
        tt_rgs = [(f, rg_i) for f in tt_files for rg_i in range(pq.ParquetFile(f).metadata.num_row_groups)]

        use_pool = threads > 1 and (len(mc_rgs) >= 2 or len(tt_rgs) >= 2 or n_buckets >= 2)
        pool_ctx = mp.get_context("spawn") if use_pool else None
        if use_pool and script_dir not in sys.path:
            sys.path.append(script_dir)
            added_syspath = True

        with contextlib.ExitStack() as pstack:
            ex = pstack.enter_context(ProcessPoolExecutor(max_workers=threads, mp_context=pool_ctx)) if use_pool else None

            # --- Phase 2A: ORF scatter across worker slices ---
            mc_slices_raw = _partition_rgs(mc_rgs, threads if ex is not None else 1)
            mc_slices: list[tuple[list[tuple[Path, int]], tuple[int, int] | None]] = []
            for sl in mc_slices_raw:
                last_u = _last_nr50_member_uid(sl)
                if last_u is not None:
                    mc_slices.append((sl, last_u))

            stage_occ_all = np.zeros(n_stages, dtype=np.int64)
            stage_occ_comp = np.zeros(n_stages, dtype=np.int64)
            nr50_members = 0
            n_orf_workers = len(mc_slices)

            if n_orf_workers == 0:
                # Zero nr50 proteins: still verify protein_occurrence is empty
                res0 = _scatter_orf_worker({
                    "worker_id": 0, "stages": stages, "mc_rgs": [], "po_by_stage": po_by_stage,
                    "lo_pu": None, "hi_pu": None, "is_last": True,
                    "n_buckets": n_buckets, "work_dir": work_dir, "n_cdim": n_cdim,
                })
                worker_peaks_kib.append(res0["vmhwm_kib"])
                n_orf_workers = 1
            else:
                orf_tasks = []
                for w_i, (sl, hi_pair) in enumerate(mc_slices):
                    lo_pair = None if w_i == 0 else mc_slices[w_i - 1][1]
                    orf_tasks.append({
                        "worker_id": w_i,
                        "stages": stages,
                        "mc_rgs": sl,
                        "po_by_stage": po_by_stage,
                        "lo_pu": lo_pair,
                        "hi_pu": hi_pair,
                        "is_last": w_i == n_orf_workers - 1,
                        "n_buckets": n_buckets,
                        "work_dir": work_dir,
                        "n_cdim": n_cdim,
                    })
                results = list(ex.map(_scatter_orf_worker, orf_tasks)) if (ex is not None and n_orf_workers > 1) else [
                    _scatter_orf_worker(t) for t in orf_tasks
                ]
                for r in results:
                    nr50_members += int(r["nr50_members"])
                    stage_occ_all += r["stage_occ_all"]
                    stage_occ_comp += r["stage_occ_comp"]
                    worker_peaks_kib.append(int(r["vmhwm_kib"]))

            # --- Phase 2B: ORF bucket reduction ---
            orf_reduce_jobs = [
                ([work_dir / f"orf-w{w:02d}-b{b:03d}.bin" for w in range(n_orf_workers)], n_stages, earlier_prefix_len)
                for b in range(n_buckets)
            ]
            orf_acc_all = _CollectorAcc.zeros(n_stages, earlier_prefix_len)
            orf_acc_comp = _CollectorAcc.zeros(n_stages, earlier_prefix_len)
            orf_red_res = (
                list(ex.map(_reduce_orf_bucket_job, orf_reduce_jobs))
                if (ex is not None and n_buckets > 1)
                else [_reduce_orf_bucket_job(j) for j in orf_reduce_jobs]
            )
            for a_all, a_comp, pk in orf_red_res:
                orf_acc_all.iadd(a_all)
                orf_acc_comp.iadd(a_comp)
                worker_peaks_kib.append(pk)

            # --- Phase 3A: Syntenic pair scatter across worker slices ---
            tt_slices_raw = _partition_rgs(tt_rgs, threads if ex is not None else 1)
            tt_slices: list[tuple[list[tuple[Path, int]], tuple[int, int] | None]] = []
            for sl in tt_slices_raw:
                last_c = _last_tt_contig_uid(sl)
                if last_c is not None:
                    tt_slices.append((sl, last_c))

            n_pair_workers = len(tt_slices)
            if n_pair_workers == 0:
                pk0 = _scatter_pair_worker({
                    "worker_id": 0, "stages": stages, "tt_rgs": [], "ct_by_stage": ct_by_stage,
                    "lo_cu": None, "hi_cu": None, "is_last": True,
                    "n_buckets": n_buckets, "work_dir": work_dir, "n_cdim": n_cdim,
                })
                worker_peaks_kib.append(pk0)
                n_pair_workers = 1
            else:
                pair_tasks = []
                for w_i, (sl, hi_pair) in enumerate(tt_slices):
                    lo_pair = None if w_i == 0 else tt_slices[w_i - 1][1]
                    pair_tasks.append({
                        "worker_id": w_i,
                        "stages": stages,
                        "tt_rgs": sl,
                        "ct_by_stage": ct_by_stage,
                        "lo_cu": lo_pair,
                        "hi_cu": hi_pair,
                        "is_last": w_i == n_pair_workers - 1,
                        "n_buckets": n_buckets,
                        "work_dir": work_dir,
                        "n_cdim": n_cdim,
                    })
                p_peaks = list(ex.map(_scatter_pair_worker, pair_tasks)) if (ex is not None and n_pair_workers > 1) else [
                    _scatter_pair_worker(t) for t in pair_tasks
                ]
                worker_peaks_kib.extend(int(pk) for pk in p_peaks)

            # --- Phase 3B: Syntenic pair bucket reduction ---
            pair_reduce_jobs = [
                ([work_dir / f"pair-w{w:02d}-b{b:03d}.bin" for w in range(n_pair_workers)], n_stages, earlier_prefix_len)
                for b in range(n_buckets)
            ]
            pair_acc = _CollectorAcc.zeros(n_stages, earlier_prefix_len)
            pair_red_res = (
                list(ex.map(_reduce_pair_bucket_job, pair_reduce_jobs))
                if (ex is not None and n_buckets > 1)
                else [_reduce_pair_bucket_job(j) for j in pair_reduce_jobs]
            )
            for sub_acc, pk in pair_red_res:
                pair_acc.iadd(sub_acc)
                worker_peaks_kib.append(pk)
    finally:
        if added_syspath:
            with contextlib.suppress(ValueError):
                sys.path.remove(script_dir)
        shutil.rmtree(work_dir, ignore_errors=True)

    # Assemble deterministic stage_curve, totals, and since_earlier
    all_series = _cumulative_series(orf_acc_all)
    comp_series = _cumulative_series(orf_acc_comp)
    pair_series_raw = _cumulative_series(pair_acc)

    cum_contigs = np.cumsum(stage_contigs)
    cum_bases = np.cumsum(stage_bases)

    def _rate(val: int, bases: int) -> float | None:
        return None if bases == 0 else float(val) / (float(bases) / 1e9)

    stage_curve: list[dict] = []
    for m in range(n_stages):
        sb = int(stage_bases[m])
        cb = int(cum_bases[m])
        net3 = int(orf_acc_comp.ge3_by_stage[m])
        b3 = int(orf_acc_comp.brand_new_ge3[m])
        p3 = int(orf_acc_comp.promoted_ge3[m])
        net5 = int(orf_acc_comp.ge5_by_stage[m])
        b5 = int(orf_acc_comp.brand_new_ge5[m])
        p5 = int(orf_acc_comp.promoted_ge5[m])

        pnet3 = int(pair_acc.ge3_by_stage[m])
        pb3 = int(pair_acc.brand_new_ge3[m])
        pp3 = int(pair_acc.promoted_ge3[m])
        pnet5 = int(pair_acc.ge5_by_stage[m])
        pb5 = int(pair_acc.brand_new_ge5[m])
        pp5 = int(pair_acc.promoted_ge5[m])

        occ_a = int(stage_occ_all[m])
        occ_c = int(stage_occ_comp[m])
        stage_curve.append({
            "stage": stages[m],
            "stage_index": m,
            "stored_contigs": int(stage_contigs[m]),
            "stored_bases": sb,
            "stored_gbp": sb / 1e9,
            "cum_stored_contigs": int(cum_contigs[m]),
            "cum_stored_bases": cb,
            "cum_stored_gbp": cb / 1e9,
            "protein_occurrences_all": occ_a,
            "protein_occurrences_complete": occ_c,
            "protein_occurrences_edge_truncated": occ_a - occ_c,
            "all_orf": all_series[m],
            "complete_orf": comp_series[m],
            "syntenic_pairs": _pair_dict_from_stats(pair_series_raw[m]),
            "marginal": {
                "new_all_orf_ge1": int(orf_acc_all.ge1_by_stage[m]),
                "new_complete_orf_ge1": int(orf_acc_comp.ge1_by_stage[m]),
                "net_new_ge3": net3,
                "brand_new_ge3": b3,
                "promoted_to_ge3": p3,
                "net_new_ge3_per_gbp": _rate(net3, sb),
                "brand_new_ge3_per_gbp": _rate(b3, sb),
                "promoted_to_ge3_per_gbp": _rate(p3, sb),
                "net_new_ge5": net5,
                "brand_new_ge5": b5,
                "promoted_to_ge5": p5,
                "net_new_ge5_per_gbp": _rate(net5, sb),
                "brand_new_ge5_per_gbp": _rate(b5, sb),
                "promoted_to_ge5_per_gbp": _rate(p5, sb),
                "new_pairs_ge1": int(pair_acc.ge1_by_stage[m]),
                "net_new_pairs_ge3": pnet3,
                "brand_new_pairs_ge3": pb3,
                "promoted_to_pairs_ge3": pp3,
                "net_new_pairs_ge3_per_gbp": _rate(pnet3, sb),
                "brand_new_pairs_ge3_per_gbp": _rate(pb3, sb),
                "promoted_to_pairs_ge3_per_gbp": _rate(pp3, sb),
                "net_new_pairs_ge5": pnet5,
                "brand_new_pairs_ge5": pb5,
                "promoted_to_pairs_ge5": pp5,
                "net_new_pairs_ge5_per_gbp": _rate(pnet5, sb),
            },
        })

    since_earlier = None
    if earlier_stages is not None and earlier_prefix_len is not None:
        prev_bases = 0 if earlier_prefix_len == 0 else int(cum_bases[earlier_prefix_len - 1])
        new_bases = int(cum_bases[-1]) - prev_bases
        prev_comp_s = 0 if earlier_prefix_len == 0 else comp_series[earlier_prefix_len - 1]["s_obs"]
        prev_pair_s = 0 if earlier_prefix_len == 0 else pair_series_raw[earlier_prefix_len - 1]["s_obs"]
        prev_c_hat = 1.0 if earlier_prefix_len == 0 else comp_series[earlier_prefix_len - 1]["c_hat"]
        prev_p_chat = 1.0 if earlier_prefix_len == 0 else pair_series_raw[earlier_prefix_len - 1]["c_hat"]

        e_b3 = int(orf_acc_comp.earlier_brand_new_ge3)
        e_p3 = int(orf_acc_comp.earlier_promoted_ge3)
        e_net3 = e_b3 + e_p3
        e_b5 = int(orf_acc_comp.earlier_brand_new_ge5)
        e_p5 = int(orf_acc_comp.earlier_promoted_ge5)
        e_net5 = e_b5 + e_p5
        ep_b3 = int(pair_acc.earlier_brand_new_ge3)
        ep_p3 = int(pair_acc.earlier_promoted_ge3)
        ep_net3 = ep_b3 + ep_p3
        ep_b5 = int(pair_acc.earlier_brand_new_ge5)
        ep_p5 = int(pair_acc.earlier_promoted_ge5)
        ep_net5 = ep_b5 + ep_p5

        since_earlier = {
            "earlier": earlier,
            "earlier_stages": earlier_stages,
            "new_stages": stages[earlier_prefix_len:],
            "new_stored_bases": new_bases,
            "new_stored_gbp": new_bases / 1e9,
            "net_new_complete_ge1": comp_series[-1]["s_obs"] - prev_comp_s,
            "net_new_ge3": e_net3,
            "brand_new_ge3": e_b3,
            "promoted_to_ge3": e_p3,
            "net_new_ge3_per_gbp": _rate(e_net3, new_bases),
            "brand_new_ge3_per_gbp": _rate(e_b3, new_bases),
            "promoted_to_ge3_per_gbp": _rate(e_p3, new_bases),
            "net_new_ge5": e_net5,
            "brand_new_ge5": e_b5,
            "promoted_to_ge5": e_p5,
            "net_new_ge5_per_gbp": _rate(e_net5, new_bases),
            "net_new_pairs_ge1": pair_series_raw[-1]["s_obs"] - prev_pair_s,
            "net_new_pairs_ge3": ep_net3,
            "brand_new_pairs_ge3": ep_b3,
            "promoted_to_pairs_ge3": ep_p3,
            "net_new_pairs_ge3_per_gbp": _rate(ep_net3, new_bases),
            "net_new_pairs_ge5": ep_net5,
            "brand_new_pairs_ge5": ep_b5,
            "promoted_to_pairs_ge5": ep_p5,
            "c_hat_prefix_before": prev_c_hat,
            "c_hat_after": comp_series[-1]["c_hat"],
            "c_hat_pairs_prefix_before": prev_p_chat,
            "c_hat_pairs_after": pair_series_raw[-1]["c_hat"],
        }

    tot_occ_all = int(stage_occ_all.sum())
    tot_occ_comp = int(stage_occ_comp.sum())
    tot_bases = int(cum_bases[-1])
    doc = {
        "schema_version": 1,
        "release": release,
        "earlier": earlier,
        "tape_version": tv,
        "gene_call_version": gcv,
        "stages": stages,
        "totals": {
            "stages_count": n_stages,
            "stored_contigs": int(cum_contigs[-1]),
            "stored_bases": tot_bases,
            "stored_gbp": tot_bases / 1e9,
            "distinct_genomes_contig_dim": n_genomes,
            "nr50_members": nr50_members,
            "protein_occurrences_all": tot_occ_all,
            "protein_occurrences_complete": tot_occ_comp,
            "protein_occurrences_edge_truncated": tot_occ_all - tot_occ_comp,
            "all_orf": all_series[-1],
            "complete_orf": comp_series[-1],
            "syntenic_pairs": _pair_dict_from_stats(pair_series_raw[-1]),
        },
        "since_earlier": since_earlier,
        "stage_curve": stage_curve,
    }

    out.mkdir(parents=True, exist_ok=True)
    atomic_write_text(out / "saturation.json", json.dumps(doc, indent=2, sort_keys=True) + "\n")
    atomic_write_text(out / "saturation.txt", format_saturation_txt(doc))

    wall_s = round(time.monotonic() - t0, 3)
    r_child = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu_total_s = round(
        (time.process_time() - cpu0) + max(0.0, float(r_child.ru_utime + r_child.ru_stime) - child_cpu0),
        3,
    )
    self_vmhwm_kib = read_vmhwm_kib()
    max_worker_kib = max(worker_peaks_kib, default=0) if use_pool else 0
    # Conservative simultaneous resident upper bound: parent VmHWM + threads * max worker VmHWM when pool used
    est_peak_kib = self_vmhwm_kib + (threads * max_worker_kib if use_pool else 0)
    vmhwm_report_kib = max(self_vmhwm_kib, max_worker_kib)
    stats = {
        "release": release,
        "threads": threads,
        "buckets": n_buckets,
        "wall_s": wall_s,
        "cpu_s": cpu_total_s,
        "cpu_min": round(cpu_total_s / 60.0, 4),
        "vmhwm_self_kib": self_vmhwm_kib,
        "vmhwm_max_worker_kib": max_worker_kib,
        "vmhwm_upper_bound_kib": est_peak_kib,
        "vmhwm_kib": vmhwm_report_kib,
        "vmhwm_mib": round(vmhwm_report_kib / 1024.0, 2),
        "vmhwm_upper_bound_mib": round(est_peak_kib / 1024.0, 2),
        "vmhwm_bytes": vmhwm_report_kib * 1024,
    }
    if stats_out is not None:
        atomic_write_text(stats_out, json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(
        f"release_saturation: {release} done in {wall_s:.2f}s wall ({cpu_total_s:.2f}s CPU), "
        f"max_proc_VmHWM={stats['vmhwm_mib']:.1f} MiB (upper_bound={stats['vmhwm_upper_bound_mib']:.1f} MiB), "
        f"threads={threads}, buckets={n_buckets}",
        file=sys.stderr,
    )
    return doc


# ---------------------------------------------------------------------------------------------------------------
# `history`: build and validate `planning/measurements/release_saturation_history.json` (AGENTS.md §6)
# ---------------------------------------------------------------------------------------------------------------


def _meas(val: int | float, host: str, date_str: str, cmd: str, source: str) -> dict:
    return {"value": val, "label": "measured", "host": host, "date": date_str, "command": cmd, "source": source}


def _deriv(val: int | float | None, formula: str, source: str) -> dict | None:
    if val is None:
        return None
    return {"value": val, "label": "derived", "formula": formula, "source": source}


def _wrap_stats_block(st: dict, host: str, date_str: str, cmd: str, src_prefix: str, *, pair_keys: bool = False) -> dict:
    if pair_keys:
        return {
            "pairs_ge1": _meas(st["pairs_ge1"], host, date_str, cmd, f"{src_prefix}: distinct canonical adjacent NR50 pairs (k>=1 genomes)"),
            "pairs_f1": _meas(st["pairs_f1"], host, date_str, cmd, f"{src_prefix}: singleton canonical pairs (1 genome)"),
            "pairs_f2": _meas(st["pairs_f2"], host, date_str, cmd, f"{src_prefix}: doubleton canonical pairs (2 genomes)"),
            "pairs_ge3": _meas(st["pairs_ge3"], host, date_str, cmd, f"{src_prefix}: canonical pairs in >=3 distinct genomes"),
            "pairs_ge5": _meas(st["pairs_ge5"], host, date_str, cmd, f"{src_prefix}: canonical pairs in >=5 distinct genomes"),
            "n_incidences": _meas(st["n_incidences"], host, date_str, cmd, f"{src_prefix}: distinct (canonical_pair, genome_uid) incidences"),
            "c_hat_pairs": _deriv(st["c_hat_pairs"], "1 - pairs_f1 / n_incidences (1.0 when n_incidences == 0)", f"Good-Turing pair coverage ({src_prefix})"),
            "chao1_pairs": _deriv(st["chao1_pairs"], "pairs_ge1 + pairs_f1 * (pairs_f1 - 1) / (2 * (pairs_f2 + 1))", f"bias-corrected Chao1 pair richness ({src_prefix})"),
        }
    return {
        "s_obs": _meas(st["s_obs"], host, date_str, cmd, f"{src_prefix}: distinct NR50 clusters observed (k>=1 genomes)"),
        "f1": _meas(st["f1"], host, date_str, cmd, f"{src_prefix}: singleton NR50 clusters (1 genome)"),
        "f2": _meas(st["f2"], host, date_str, cmd, f"{src_prefix}: doubleton NR50 clusters (2 genomes)"),
        "ge3": _meas(st["ge3"], host, date_str, cmd, f"{src_prefix}: NR50 clusters in >=3 distinct genomes"),
        "ge5": _meas(st["ge5"], host, date_str, cmd, f"{src_prefix}: NR50 clusters in >=5 distinct genomes"),
        "n_incidences": _meas(st["n_incidences"], host, date_str, cmd, f"{src_prefix}: distinct (cluster_uid, genome_uid) incidences"),
        "c_hat": _deriv(st["c_hat"], "1 - f1 / n_incidences (1.0 when n_incidences == 0)", f"Good-Turing sample coverage ({src_prefix})"),
        "chao1": _deriv(st["chao1"], "s_obs + f1 * (f1 - 1) / (2 * (f2 + 1))", f"bias-corrected Chao1 richness ({src_prefix})"),
    }


def _stage_class_map(diversity_path: Path | None) -> tuple[dict[str, str], dict[str, str]]:
    if diversity_path is not None and Path(diversity_path).is_file():
        div = json.loads(Path(diversity_path).read_text(encoding="utf-8"))
        st_cls = {s["stage"]: s["class"] for s in div.get("stages", [])}
        cls_color = {c["id"]: c["color"] for c in div.get("classes", [])}
        return st_cls, {**CLASS_COLORS, **cls_color}
    return {}, dict(CLASS_COLORS)


def add_to_history(
    history_path: Path,
    sat_doc: dict,
    *,
    host: str,
    date_str: str,
    command: str,
    diversity_path: Path | None = DEFAULT_DIVERSITY,
    set_collector: bool = False,
    stats: dict | None = None,
) -> dict:
    st_cls, cls_colors = _stage_class_map(diversity_path)
    if history_path.is_file():
        hist = json.loads(history_path.read_text(encoding="utf-8"))
    else:
        hist = {
            "schema_version": 1,
            "description": (
                "Release-level protein-family and syntenic-pair saturation history (card A-33, DEC-86 (1)). "
                "Every numeric entry carries value, label (measured, derived, or estimate; AGENTS.md section 6), and source."
            ),
            "figure_title": "bio-box NR50 family and syntenic-pair saturation across stages and releases",
            "footnote": (
                "Complete ORFs exclude linear-contig edge fragments (dist_to_contig_end_bp == 0) and retain all circular-contig "
                "CDS; incidences count distinct contig_dim.genome_uid per NR50 cluster or strand-canonical adjacent NR50 pair."
            ),
            "class_colors": cls_colors,
            "collector_release": None,
            "releases": [],
        }
    hist["class_colors"] = cls_colors

    rel = sat_doc["release"]
    tot = sat_doc["totals"]
    src_base = f"release_saturation.py measure on {rel} ({host}, {date_str})"

    stage_entries = []
    for row in sat_doc["stage_curve"]:
        st = row["stage"]
        mg = row["marginal"]
        s_src = f"{src_base}, stage {st} (prefix 0..{row['stage_index']})"
        stage_entries.append({
            "stage": st,
            "stage_index": row["stage_index"],
            "class": st_cls.get(st, "other"),
            "stored_contigs": _meas(row["stored_contigs"], host, date_str, command, f"{s_src}: default-variant contigs in contig_index"),
            "stored_bases": _meas(row["stored_bases"], host, date_str, command, f"{s_src}: default-variant bases in contig_index"),
            "stored_gbp": _deriv(row["stored_gbp"], "stored_bases / 1e9", f"{s_src}: stage stored Gbp"),
            "cum_stored_contigs": _deriv(row["cum_stored_contigs"], "sum(stored_contigs[0..m])", f"{s_src}: cumulative stored contigs"),
            "cum_stored_bases": _deriv(row["cum_stored_bases"], "sum(stored_bases[0..m])", f"{s_src}: cumulative stored bases"),
            "cum_stored_gbp": _deriv(row["cum_stored_gbp"], "cum_stored_bases / 1e9", f"{s_src}: cumulative stored Gbp"),
            "protein_occurrences_all": _meas(row["protein_occurrences_all"], host, date_str, command, f"{s_src}: protein_occurrence rows"),
            "protein_occurrences_complete": _meas(row["protein_occurrences_complete"], host, date_str, command, f"{s_src}: protein_occurrence rows with dist_to_contig_end_bp > 0"),
            "protein_occurrences_edge_truncated": _deriv(row["protein_occurrences_edge_truncated"], "protein_occurrences_all - protein_occurrences_complete", f"{s_src}: linear-contig edge-truncated CDS occurrences"),
            "all_orf": _wrap_stats_block(row["all_orf"], host, date_str, command, f"{s_src} all_orf"),
            "complete_orf": _wrap_stats_block(row["complete_orf"], host, date_str, command, f"{s_src} complete_orf"),
            "syntenic_pairs": _wrap_stats_block(row["syntenic_pairs"], host, date_str, command, f"{s_src} syntenic_pairs", pair_keys=True),
            "marginal": {
                "new_all_orf_ge1": _meas(mg["new_all_orf_ge1"], host, date_str, command, f"{s_src}: marginal new NR50 clusters (all ORFs)"),
                "new_complete_orf_ge1": _meas(mg["new_complete_orf_ge1"], host, date_str, command, f"{s_src}: marginal new complete-ORF NR50 clusters"),
                "brand_new_ge3": _meas(mg["brand_new_ge3"], host, date_str, command, f"{s_src}: complete-ORF NR50 clusters jumping 0 -> >=3 genomes"),
                "promoted_to_ge3": _meas(mg["promoted_to_ge3"], host, date_str, command, f"{s_src}: complete-ORF NR50 clusters promoted 1..2 -> >=3 genomes"),
                "net_new_ge3": _deriv(mg["net_new_ge3"], "brand_new_ge3 + promoted_to_ge3", f"{s_src}: marginal k>=3 complete-ORF NR50 clusters"),
                "net_new_ge3_per_gbp": _deriv(mg["net_new_ge3_per_gbp"], "net_new_ge3 / (stored_bases / 1e9)", f"{s_src}: marginal k>=3 complete-ORF NR50 per Gbp"),
                "brand_new_ge3_per_gbp": _deriv(mg["brand_new_ge3_per_gbp"], "brand_new_ge3 / (stored_bases / 1e9)", f"{s_src}: brand-new k>=3 complete-ORF NR50 per Gbp"),
                "promoted_to_ge3_per_gbp": _deriv(mg["promoted_to_ge3_per_gbp"], "promoted_to_ge3 / (stored_bases / 1e9)", f"{s_src}: promoted k>=3 complete-ORF NR50 per Gbp"),
                "brand_new_ge5": _meas(mg["brand_new_ge5"], host, date_str, command, f"{s_src}: complete-ORF NR50 clusters jumping 0 -> >=5 genomes"),
                "promoted_to_ge5": _meas(mg["promoted_to_ge5"], host, date_str, command, f"{s_src}: complete-ORF NR50 clusters promoted 1..4 -> >=5 genomes"),
                "net_new_ge5": _deriv(mg["net_new_ge5"], "brand_new_ge5 + promoted_to_ge5", f"{s_src}: marginal k>=5 complete-ORF NR50 clusters"),
                "net_new_ge5_per_gbp": _deriv(mg["net_new_ge5_per_gbp"], "net_new_ge5 / (stored_bases / 1e9)", f"{s_src}: marginal k>=5 complete-ORF NR50 per Gbp"),
                "new_pairs_ge1": _meas(mg["new_pairs_ge1"], host, date_str, command, f"{s_src}: marginal new canonical adjacent NR50 pairs"),
                "brand_new_pairs_ge3": _meas(mg["brand_new_pairs_ge3"], host, date_str, command, f"{s_src}: canonical pairs jumping 0 -> >=3 genomes"),
                "promoted_to_pairs_ge3": _meas(mg["promoted_to_pairs_ge3"], host, date_str, command, f"{s_src}: canonical pairs promoted 1..2 -> >=3 genomes"),
                "net_new_pairs_ge3": _deriv(mg["net_new_pairs_ge3"], "brand_new_pairs_ge3 + promoted_to_pairs_ge3", f"{s_src}: marginal k>=3 canonical pairs"),
                "net_new_pairs_ge3_per_gbp": _deriv(mg["net_new_pairs_ge3_per_gbp"], "net_new_pairs_ge3 / (stored_bases / 1e9)", f"{s_src}: marginal k>=3 canonical pairs per Gbp"),
                "brand_new_pairs_ge3_per_gbp": _deriv(mg["brand_new_pairs_ge3_per_gbp"], "brand_new_pairs_ge3 / (stored_bases / 1e9)", f"{s_src}: brand-new k>=3 canonical pairs per Gbp"),
                "promoted_to_pairs_ge3_per_gbp": _deriv(mg["promoted_to_pairs_ge3_per_gbp"], "promoted_to_pairs_ge3 / (stored_bases / 1e9)", f"{s_src}: promoted k>=3 canonical pairs per Gbp"),
                "brand_new_pairs_ge5": _meas(mg["brand_new_pairs_ge5"], host, date_str, command, f"{s_src}: canonical pairs jumping 0 -> >=5 genomes"),
                "promoted_to_pairs_ge5": _meas(mg["promoted_to_pairs_ge5"], host, date_str, command, f"{s_src}: canonical pairs promoted 1..4 -> >=5 genomes"),
                "net_new_pairs_ge5": _deriv(mg["net_new_pairs_ge5"], "brand_new_pairs_ge5 + promoted_to_pairs_ge5", f"{s_src}: marginal k>=5 canonical pairs"),
            },
        })

    se = sat_doc.get("since_earlier")
    se_wrapped = None
    if se is not None:
        se_src = f"{src_base} vs earlier prefix {se['earlier']}"
        se_wrapped = {
            "earlier": se["earlier"],
            "earlier_stages": se["earlier_stages"],
            "new_stages": se["new_stages"],
            "new_stored_bases": _deriv(se["new_stored_bases"], "tot_stored_bases - earlier_prefix_bases", f"{se_src}: new stored bases"),
            "new_stored_gbp": _deriv(se["new_stored_gbp"], "new_stored_bases / 1e9", f"{se_src}: new stored Gbp"),
            "brand_new_ge3": _meas(se["brand_new_ge3"], host, date_str, command, f"{se_src}: complete-ORF NR50 jumping 0 -> >=3 genomes"),
            "promoted_to_ge3": _meas(se["promoted_to_ge3"], host, date_str, command, f"{se_src}: complete-ORF NR50 promoted 1..2 -> >=3 genomes"),
            "net_new_ge3": _deriv(se["net_new_ge3"], "brand_new_ge3 + promoted_to_ge3", f"{se_src}: net new k>=3 complete-ORF NR50"),
            "net_new_ge3_per_gbp": _deriv(se["net_new_ge3_per_gbp"], "net_new_ge3 / (new_stored_bases / 1e9)", f"{se_src}: net new k>=3 NR50 per Gbp"),
            "brand_new_ge5": _meas(se["brand_new_ge5"], host, date_str, command, f"{se_src}: complete-ORF NR50 jumping 0 -> >=5 genomes"),
            "promoted_to_ge5": _meas(se["promoted_to_ge5"], host, date_str, command, f"{se_src}: complete-ORF NR50 promoted 1..4 -> >=5 genomes"),
            "net_new_ge5": _deriv(se["net_new_ge5"], "brand_new_ge5 + promoted_to_ge5", f"{se_src}: net new k>=5 complete-ORF NR50"),
            "net_new_ge5_per_gbp": _deriv(se["net_new_ge5_per_gbp"], "net_new_ge5 / (new_stored_bases / 1e9)", f"{se_src}: net new k>=5 NR50 per Gbp"),
            "brand_new_pairs_ge3": _meas(se["brand_new_pairs_ge3"], host, date_str, command, f"{se_src}: canonical pairs jumping 0 -> >=3 genomes"),
            "promoted_to_pairs_ge3": _meas(se["promoted_to_pairs_ge3"], host, date_str, command, f"{se_src}: canonical pairs promoted 1..2 -> >=3 genomes"),
            "net_new_pairs_ge3": _deriv(se["net_new_pairs_ge3"], "brand_new_pairs_ge3 + promoted_to_pairs_ge3", f"{se_src}: net new k>=3 canonical pairs"),
            "net_new_pairs_ge3_per_gbp": _deriv(se["net_new_pairs_ge3_per_gbp"], "net_new_pairs_ge3 / (new_stored_bases / 1e9)", f"{se_src}: net new k>=3 canonical pairs per Gbp"),
        }

    run_stats_wrapped = None
    if stats is not None:
        run_stats_wrapped = {
            "threads": _meas(int(stats["threads"]), host, date_str, command, f"{src_base}: worker thread cap"),
            "buckets": _meas(int(stats["buckets"]), host, date_str, command, f"{src_base}: hash partition count"),
            "wall_s": _meas(float(stats["wall_s"]), host, date_str, command, f"{src_base}: wall-clock seconds"),
            "cpu_s": _meas(float(stats["cpu_s"]), host, date_str, command, f"{src_base}: total CPU seconds"),
            "cpu_min": _deriv(float(stats["cpu_s"]) / 60.0, "cpu_s / 60", f"{src_base}: total CPU minutes"),
            "vmhwm_kib": _meas(int(stats["vmhwm_kib"]), host, date_str, command, f"{src_base}: max per-process /proc/self/status VmHWM KiB"),
            "vmhwm_mib": _deriv(float(stats["vmhwm_kib"]) / 1024.0, "vmhwm_kib / 1024", f"{src_base}: max per-process VmHWM MiB"),
        }
        if "cgroup_peak_bytes" in stats:
            run_stats_wrapped["cgroup_peak_bytes"] = _meas(
                int(stats["cgroup_peak_bytes"]), host, date_str, command, f"{src_base}: systemd unit MemoryPeak bytes"
            )
            run_stats_wrapped["cgroup_peak_gib"] = _deriv(
                float(stats["cgroup_peak_bytes"]) / (1024.0 ** 3), "cgroup_peak_bytes / 2^30", f"{src_base}: systemd unit MemoryPeak GiB"
            )

    entry = {
        "release": rel,
        "earlier": sat_doc["earlier"],
        "tape_version": sat_doc["tape_version"],
        "gene_call_version": sat_doc["gene_call_version"],
        "stages": sat_doc["stages"],
        "totals": {
            "stages_count": _meas(tot["stages_count"], host, date_str, command, f"{src_base}: number of stages in manifest"),
            "stored_contigs": _meas(tot["stored_contigs"], host, date_str, command, f"{src_base}: default-variant contigs"),
            "stored_bases": _meas(tot["stored_bases"], host, date_str, command, f"{src_base}: default-variant bases"),
            "stored_gbp": _deriv(tot["stored_gbp"], "stored_bases / 1e9", f"{src_base}: stored sequence in Gbp"),
            "distinct_genomes_contig_dim": _meas(tot["distinct_genomes_contig_dim"], host, date_str, command, f"{src_base}: distinct genome_uid in contig_dim"),
            "nr50_members": _meas(tot["nr50_members"], host, date_str, command, f"{src_base}: nr50 rows in member_cluster"),
            "protein_occurrences_all": _meas(tot["protein_occurrences_all"], host, date_str, command, f"{src_base}: total protein_occurrence rows"),
            "protein_occurrences_complete": _meas(tot["protein_occurrences_complete"], host, date_str, command, f"{src_base}: complete-ORF protein_occurrence rows"),
            "protein_occurrences_edge_truncated": _deriv(tot["protein_occurrences_edge_truncated"], "protein_occurrences_all - protein_occurrences_complete", f"{src_base}: linear-edge-truncated protein_occurrence rows"),
            "all_orf": _wrap_stats_block(tot["all_orf"], host, date_str, command, f"{src_base} all_orf"),
            "complete_orf": _wrap_stats_block(tot["complete_orf"], host, date_str, command, f"{src_base} complete_orf"),
            "syntenic_pairs": _wrap_stats_block(tot["syntenic_pairs"], host, date_str, command, f"{src_base} syntenic_pairs", pair_keys=True),
        },
        "since_earlier": se_wrapped,
        "run_stats": run_stats_wrapped,
        "stage_curve": stage_entries,
    }

    existing = [r for r in hist.get("releases", []) if r.get("release") != rel]
    existing.append(entry)
    hist["releases"] = existing
    if set_collector or not hist.get("collector_release"):
        hist["collector_release"] = rel

    atomic_write_text(history_path, json.dumps(hist, indent=2) + "\n")
    return hist


def walk_history_numbers(obj: object, where: str = "") -> list[tuple[str, dict]]:
    found: list[tuple[str, dict]] = []
    if isinstance(obj, dict):
        if "value" in obj and "label" in obj:
            found.append((where, obj))
            return found
        for k, v in obj.items():
            found.extend(walk_history_numbers(v, f"{where}.{k}" if where else str(k)))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            tag = (v.get("release") or v.get("stage")) if isinstance(v, dict) else None
            found.extend(walk_history_numbers(v, f"{where}[{tag or i}]"))
    return found


def _check_stats_consistency(blk: dict, where: str, *, pair_keys: bool = False) -> list[str]:
    problems: list[str] = []
    s_obs = blk["pairs_ge1" if pair_keys else "s_obs"]["value"]
    f1 = blk["pairs_f1" if pair_keys else "f1"]["value"]
    f2 = blk["pairs_f2" if pair_keys else "f2"]["value"]
    n_inc = blk["n_incidences"]["value"]
    c_hat = blk["c_hat_pairs" if pair_keys else "c_hat"]["value"]
    chao1 = blk["chao1_pairs" if pair_keys else "chao1"]["value"]
    exp_chat = 1.0 if n_inc == 0 else 1.0 - (f1 / n_inc)
    exp_chao1 = float(s_obs) + (f1 * (f1 - 1)) / (2.0 * (f2 + 1))
    if not math.isclose(c_hat, exp_chat, rel_tol=1e-9, abs_tol=1e-12):
        problems.append(f"{where}: c_hat {c_hat} != expected {exp_chat}")
    if not math.isclose(chao1, exp_chao1, rel_tol=1e-9, abs_tol=1e-6):
        problems.append(f"{where}: chao1 {chao1} != expected {exp_chao1}")
    return problems


def check_history(hist: dict) -> list[str]:
    problems: list[str] = []
    nums = walk_history_numbers(hist)
    if not nums:
        return ["history JSON contains zero labelled number records"]

    for where, num in nums:
        lbl = num.get("label")
        if lbl not in LABELS:
            problems.append(f"{where}: invalid label {lbl!r}")
        if not str(num.get("source", "")).strip():
            problems.append(f"{where}: missing source")
        val = num.get("value")
        if not isinstance(val, (int, float)) or isinstance(val, bool):
            problems.append(f"{where}: non-numeric value {val!r}")
        if lbl == "measured" and not (num.get("host") and num.get("date") and num.get("command")):
            problems.append(f"{where}: measured number missing host, date, or command")
        if lbl == "derived" and not str(num.get("formula", "")).strip():
            problems.append(f"{where}: derived number missing formula")

    for tkey in ("figure_title", "footnote"):
        txt = str(hist.get(tkey, ""))
        if "—" in txt or "–" in txt:
            problems.append(f"{tkey}: contains prose dash")

    if problems:
        return problems

    rels = {r["release"]: r for r in hist.get("releases", [])}
    coll = hist.get("collector_release")
    if coll not in rels:
        problems.append(f"collector_release {coll!r} not in releases {sorted(rels)}")

    for r in hist.get("releases", []):
        rid = r["release"]
        tot = r["totals"]
        problems.extend(_check_stats_consistency(tot["all_orf"], f"{rid}.totals.all_orf"))
        problems.extend(_check_stats_consistency(tot["complete_orf"], f"{rid}.totals.complete_orf"))
        problems.extend(_check_stats_consistency(tot["syntenic_pairs"], f"{rid}.totals.syntenic_pairs", pair_keys=True))

        cum_b = 0
        for row in r["stage_curve"]:
            st = row["stage"]
            w = f"{rid}.stage_curve[{st}]"
            sb = row["stored_bases"]["value"]
            cum_b += sb
            if row["cum_stored_bases"]["value"] != cum_b:
                problems.append(f"{w}: cum_stored_bases {row['cum_stored_bases']['value']} != {cum_b}")
            problems.extend(_check_stats_consistency(row["all_orf"], f"{w}.all_orf"))
            problems.extend(_check_stats_consistency(row["complete_orf"], f"{w}.complete_orf"))
            problems.extend(_check_stats_consistency(row["syntenic_pairs"], f"{w}.syntenic_pairs", pair_keys=True))
            mg = row["marginal"]
            if mg["brand_new_ge3"]["value"] + mg["promoted_to_ge3"]["value"] != mg["net_new_ge3"]["value"]:
                problems.append(f"{w}.marginal: brand_new_ge3 + promoted_to_ge3 != net_new_ge3")
            if mg["brand_new_pairs_ge3"]["value"] + mg["promoted_to_pairs_ge3"]["value"] != mg["net_new_pairs_ge3"]["value"]:
                problems.append(f"{w}.marginal: brand_new_pairs_ge3 + promoted_to_pairs_ge3 != net_new_pairs_ge3")
        if cum_b != tot["stored_bases"]["value"]:
            problems.append(f"{rid}: stage stored_bases sum {cum_b} != totals.stored_bases {tot['stored_bases']['value']}")
    return problems


# ---------------------------------------------------------------------------------------------------------------
# `render`: 4-panel Wilke / Okabe-Ito publication figure (.png + deterministic .svg)
# ---------------------------------------------------------------------------------------------------------------


def _style() -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Source Sans 3", "Source Sans Pro", "DejaVu Sans", "Arial"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#333333",
        "axes.linewidth": 0.8,
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "svg.fonttype": "none",
        "svg.hashsalt": "bio-box-a33",
        "savefig.dpi": 300,
    })


def _smart_stage_ticks(ax, x_gbp: np.ndarray, stages: list[str], ymax: float) -> None:
    last_x = -1e9
    span = max(float(x_gbp[-1]), 1e-6)
    min_dx = 0.065 * span
    for i, (xv, st) in enumerate(zip(x_gbp, stages)):
        is_end = i in (0, len(stages) - 1)
        if is_end or (xv - last_x >= min_dx and (x_gbp[-1] - xv >= min_dx * 0.7)):
            ax.axvline(xv, color="#EAEAEA", lw=0.6, zorder=1)
            ax.text(xv, ymax * 0.985, st, ha="center", va="top", fontsize=7.1, color="#555555")
            last_x = xv


def build_figure(hist: dict):
    import matplotlib.pyplot as plt

    rels = {r["release"]: r for r in hist["releases"]}
    coll_id = hist["collector_release"]
    coll = rels[coll_id]
    curve = coll["stage_curve"]
    n_st = len(curve)

    stages = [r["stage"] for r in curve]
    cum_gbp = np.array([float(r["cum_stored_gbp"]["value"]) for r in curve])
    scale_fam = 1e6 if float(coll["totals"]["all_orf"]["s_obs"]["value"]) >= 50_000 else 1.0
    unit_fam = "millions" if scale_fam == 1e6 else "count"

    all_ge1 = np.array([float(r["all_orf"]["s_obs"]["value"]) / scale_fam for r in curve])
    comp_ge1 = np.array([float(r["complete_orf"]["s_obs"]["value"]) / scale_fam for r in curve])
    comp_ge3 = np.array([float(r["complete_orf"]["ge3"]["value"]) / scale_fam for r in curve])
    comp_ge5 = np.array([float(r["complete_orf"]["ge5"]["value"]) / scale_fam for r in curve])

    chat_comp = np.array([100.0 * float(r["complete_orf"]["c_hat"]["value"]) for r in curve])
    chat_all = np.array([100.0 * float(r["all_orf"]["c_hat"]["value"]) for r in curve])
    chat_pair = np.array([100.0 * float(r["syntenic_pairs"]["c_hat_pairs"]["value"]) for r in curve])

    fig = plt.figure(figsize=(16.0, 11.6))
    gs = fig.add_gridspec(
        2, 2, height_ratios=[1.0, 1.12], width_ratios=[1.0, 1.0],
        hspace=0.38, wspace=0.25, left=0.085, right=0.965, top=0.895, bottom=0.095,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    c_gray = "#7F8C8D"
    c_blue = "#0072B2"
    c_verm = "#D55E00"
    c_green = "#009E73"
    c_purple = "#CC79A7"
    c_orange = "#E69F00"

    # Panel a: collector curves of NR50 families across cumulative stored Gbp
    x0 = np.r_[0.0, cum_gbp]
    ymax_a = max(float(all_ge1.max()), 1e-6) * 1.22
    _smart_stage_ticks(ax_a, cum_gbp, stages, ymax_a)
    ax_a.plot(x0, np.r_[0.0, all_ge1], color=c_gray, lw=1.7, ls="--", marker="o", ms=3.8, label="All ORFs (incl. edge fragments)", zorder=3)
    ax_a.plot(x0, np.r_[0.0, comp_ge1], color=c_blue, lw=2.2, marker="o", ms=4.5, label="Complete ORFs (k >= 1 genome)", zorder=4)
    ax_a.plot(x0, np.r_[0.0, comp_ge3], color=c_verm, lw=2.2, marker="s", ms=4.5, label="Complete ORFs, k >= 3 genomes", zorder=5)
    ax_a.plot(x0, np.r_[0.0, comp_ge5], color=c_green, lw=2.0, marker="^", ms=4.6, label="Complete ORFs, k >= 5 genomes", zorder=5)

    other_ids = sorted(k for k in rels if k != coll_id)
    for rid in other_ids:
        rdoc = rels[rid]
        xg = float(rdoc["totals"]["stored_gbp"]["value"])
        yc1 = float(rdoc["totals"]["complete_orf"]["s_obs"]["value"]) / scale_fam
        yc3 = float(rdoc["totals"]["complete_orf"]["ge3"]["value"]) / scale_fam
        ax_a.scatter([xg], [yc1], s=46, facecolors="white", edgecolors=c_blue, linewidths=1.6, zorder=6)
        ax_a.scatter([xg], [yc3], s=46, facecolors="white", edgecolors=c_verm, linewidths=1.6, zorder=6,
                     label=f"{rid} own clustering (open)" if rid == other_ids[0] else None)

    chao_comp = float(coll["totals"]["complete_orf"]["chao1"]["value"]) / scale_fam
    fmt_v = (lambda v: f"{v:.2f}M") if scale_fam == 1e6 else (lambda v: f"{v:.0f}")
    ax_a.text(
        0.02, 0.56,
        f"At {coll_id} ({cum_gbp[-1]:.1f} Gbp):\n"
        f"All ORFs: {fmt_v(all_ge1[-1])}\n"
        f"Complete ORFs: {fmt_v(comp_ge1[-1])} (Chao1 {fmt_v(chao_comp)})\n"
        f"k >= 3 genomes: {fmt_v(comp_ge3[-1])}\n"
        f"k >= 5 genomes: {fmt_v(comp_ge5[-1])}",
        transform=ax_a.transAxes, ha="left", va="top", fontsize=7.4, color="#222222", linespacing=1.26,
        bbox=dict(boxstyle="round,pad=0.30", facecolor="#F8F9FA", edgecolor="#CCCCCC", lw=0.7), zorder=7,
    )
    ax_a.set_xlim(0, max(float(cum_gbp[-1]) * 1.06, 1e-3))
    ax_a.set_ylim(0, ymax_a)
    ax_a.set_xlabel("Cumulative stored sequence (Gbp)", fontsize=9.2)
    ax_a.set_ylabel(f"Cumulative NR50 families ({unit_fam})", fontsize=9.2)
    ax_a.set_title(f"a   NR50 family collector curves across {coll_id} stages", loc="left", fontsize=10.2, fontweight="bold", pad=10)
    ax_a.legend(loc="upper left", bbox_to_anchor=(0.01, 0.94), frameon=True, facecolor="white", framealpha=0.90, edgecolor="none", fontsize=7.3)
    ax_a.grid(axis="y", color="#E5E5E5", lw=0.6)
    ax_a.set_axisbelow(True)

    # Panel b: Good-Turing sample coverage C_hat (%)
    _smart_stage_ticks(ax_b, cum_gbp, stages, 108.0)
    ax_b.plot(cum_gbp, chat_comp, color=c_blue, lw=2.3, marker="o", ms=4.8, label="Complete-ORF NR50 families", zorder=5)
    ax_b.plot(cum_gbp, chat_all, color=c_verm, lw=1.5, ls="--", marker="s", ms=3.4, label="All-ORF NR50 families (with edge fragments)", zorder=6)
    ax_b.plot(cum_gbp, chat_pair, color=c_purple, lw=2.2, marker="D", ms=4.4, label="Strand-canonical adjacent NR50 pairs", zorder=5)

    for rid in other_ids:
        rdoc = rels[rid]
        xg = float(rdoc["totals"]["stored_gbp"]["value"])
        cc = 100.0 * float(rdoc["totals"]["complete_orf"]["c_hat"]["value"])
        cp = 100.0 * float(rdoc["totals"]["syntenic_pairs"]["c_hat_pairs"]["value"])
        ax_b.scatter([xg], [cc], s=46, facecolors="white", edgecolors=c_blue, linewidths=1.6, zorder=6)
        ax_b.scatter([xg], [cp], s=46, facecolors="white", edgecolors=c_purple, linewidths=1.6, zorder=6)

    ax_b.text(
        0.02, 0.05,
        f"{coll_id} Good-Turing coverage C_hat:\n"
        f"Complete-ORF NR50: {chat_comp[-1]:.1f}%\n"
        f"All-ORF NR50: {chat_all[-1]:.1f}%\n"
        f"Syntenic NR50 pairs: {chat_pair[-1]:.1f}%",
        transform=ax_b.transAxes, ha="left", va="bottom", fontsize=7.4, color="#222222", linespacing=1.26,
        bbox=dict(boxstyle="round,pad=0.30", facecolor="#F8F9FA", edgecolor="#CCCCCC", lw=0.7), zorder=7,
    )
    ax_b.set_xlim(0, max(float(cum_gbp[-1]) * 1.06, 1e-3))
    ax_b.set_ylim(0, 108.0)
    ax_b.set_yticks([0, 20, 40, 60, 80, 100])
    ax_b.set_yticklabels(["0%", "20%", "40%", "60%", "80%", "100%"])
    ax_b.set_xlabel("Cumulative stored sequence (Gbp)", fontsize=9.2)
    ax_b.set_ylabel("Good-Turing genome-incidence coverage C_hat", fontsize=9.2)
    ax_b.set_title("b   Sample coverage of families and syntenic gene pairs", loc="left", fontsize=10.2, fontweight="bold", pad=10)
    ax_b.legend(loc="lower right", frameon=True, facecolor="white", framealpha=0.90, edgecolor="none", fontsize=7.3)
    ax_b.grid(axis="y", color="#E5E5E5", lw=0.6)
    ax_b.set_axisbelow(True)

    # Panel c: Marginal k >= 3 complete-ORF NR50 families per stored Gbp (brand-new vs promoted)
    yp = np.arange(n_st)[::-1]
    raw_bn3 = np.array([float((r["marginal"]["brand_new_ge3_per_gbp"] or {"value": 0.0})["value"]) for r in curve])
    raw_pr3 = np.array([float((r["marginal"]["promoted_to_ge3_per_gbp"] or {"value": 0.0})["value"]) for r in curve])
    rate_scale_c = 1e3 if max(float((raw_bn3 + raw_pr3).max()), 0.0) >= 1_000 else 1.0
    unit_c = "thousands / Gbp" if rate_scale_c == 1e3 else "count / Gbp"
    bn3 = raw_bn3 / rate_scale_c
    pr3 = raw_pr3 / rate_scale_c
    tot3 = bn3 + pr3

    stage_labels = [f"{r['stage']} ({r['class']})" if r.get("class") and r["class"] != "other" else r["stage"] for r in curve]
    ax_c.barh(yp, bn3, height=0.66, color=c_blue, edgecolor="white", lw=0.5, label="Brand-new (0 -> >=3 genomes)")
    ax_c.barh(yp, pr3, left=bn3, height=0.66, color=c_orange, edgecolor="white", lw=0.5, label="Promoted (1..2 -> >=3 genomes)")
    vmax_c = max(float(tot3.max()), 1.0)
    for y, v, p_v in zip(yp, tot3, pr3):
        lbl = f"{v:,.1f}k" if rate_scale_c == 1e3 else f"{v:,.0f}"
        if v > 0 and p_v / v >= 0.15:
            lbl += f" ({100.0 * p_v / v:.0f}% prom.)"
        ax_c.text(v + 0.016 * vmax_c, y, lbl, va="center", ha="left", fontsize=7.3, color="#222222")
    ax_c.set_yticks(yp)
    ax_c.set_yticklabels(stage_labels, fontsize=7.6)
    ax_c.set_xlim(0, vmax_c * 1.52)
    ax_c.set_xlabel(f"Marginal k >= 3 complete-ORF NR50 families ({unit_c})", fontsize=9.2)
    ax_c.set_title("c   Marginal multi-genome NR50 yield per Gbp, by stage", loc="left", fontsize=10.2, fontweight="bold", pad=10)
    ax_c.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.92, edgecolor="none", fontsize=7.3)
    ax_c.grid(axis="x", color="#E5E5E5", lw=0.6)
    ax_c.set_axisbelow(True)

    # Panel d: Marginal k >= 3 strand-canonical syntenic NR50 pairs per stored Gbp
    raw_pbn3 = np.array([float((r["marginal"]["brand_new_pairs_ge3_per_gbp"] or {"value": 0.0})["value"]) for r in curve])
    raw_ppr3 = np.array([float((r["marginal"]["promoted_to_pairs_ge3_per_gbp"] or {"value": 0.0})["value"]) for r in curve])
    rate_scale_d = 1e3 if max(float((raw_pbn3 + raw_ppr3).max()), 0.0) >= 1_000 else 1.0
    unit_d = "thousands / Gbp" if rate_scale_d == 1e3 else "count / Gbp"
    pbn3 = raw_pbn3 / rate_scale_d
    ppr3 = raw_ppr3 / rate_scale_d
    ptot3 = pbn3 + ppr3

    ax_d.barh(yp, pbn3, height=0.66, color=c_purple, edgecolor="white", lw=0.5, label="Brand-new pairs (0 -> >=3 genomes)")
    ax_d.barh(yp, ppr3, left=pbn3, height=0.66, color=c_green, edgecolor="white", lw=0.5, label="Promoted pairs (1..2 -> >=3 genomes)")
    vmax_d = max(float(ptot3.max()), 1.0)
    for y, v, p_v in zip(yp, ptot3, ppr3):
        lbl = f"{v:,.1f}k" if rate_scale_d == 1e3 else f"{v:,.0f}"
        if v > 0 and p_v / v >= 0.15:
            lbl += f" ({100.0 * p_v / v:.0f}% prom.)"
        ax_d.text(v + 0.016 * vmax_d, y, lbl, va="center", ha="left", fontsize=7.3, color="#222222")
    ax_d.set_yticks(yp)
    ax_d.set_yticklabels([r["stage"] for r in curve], fontsize=7.6)
    ax_d.set_xlim(0, vmax_d * 1.62)
    ax_d.set_xlabel(f"Marginal k >= 3 syntenic NR50 pairs ({unit_d})", fontsize=9.2)
    ax_d.set_title("d   Marginal multi-genome syntenic pair yield per Gbp, by stage", loc="left", fontsize=10.2, fontweight="bold", pad=10)
    ax_d.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.92, edgecolor="none", fontsize=7.3)
    ax_d.grid(axis="x", color="#E5E5E5", lw=0.6)
    ax_d.set_axisbelow(True)

    fig.suptitle(hist["figure_title"], x=0.085, y=0.968, ha="left", fontsize=12.6, fontweight="bold", color="#111111")
    fig.text(0.085, 0.016, hist["footnote"], fontsize=7.4, color="#555555", ha="left", va="bottom", linespacing=1.25)
    return fig


def render_figure(hist: dict, out_prefix: Path) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _style()
    fig = build_figure(hist)
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    png = Path(f"{out_prefix}.png")
    svg = Path(f"{out_prefix}.svg")
    fig.savefig(png, dpi=300, metadata={"Software": None})
    fig.savefig(svg, metadata={"Date": None, "Creator": None})
    plt.close(fig)
    return png, svg


# ---------------------------------------------------------------------------------------------------------------
# CLI entrypoint
# ---------------------------------------------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Card A-33: release-level saturation curve tracker")
    sub = ap.add_subparsers(dest="cmd", required=True)

    pm = sub.add_parser("measure", help="Measure deterministic saturation summary from a release manifest")
    pm.add_argument("--root", required=True, type=Path)
    pm.add_argument("--release", required=True)
    pm.add_argument("--earlier", default=None)
    pm.add_argument("--out", required=True, type=Path)
    pm.add_argument("--tmp", required=True, type=Path)
    pm.add_argument("--threads", type=int, default=4)
    pm.add_argument("--buckets", type=int, default=None)
    pm.add_argument("--stats-out", type=Path, default=None)

    ph = sub.add_parser("history", help="Update or validate release_saturation_history.json")
    ph.add_argument("--data", type=Path, default=DEFAULT_HISTORY)
    ph.add_argument("--diversity", type=Path, default=DEFAULT_DIVERSITY)
    ph.add_argument("--add", type=Path, default=None, help="saturation.json from `measure` to insert/update")
    ph.add_argument("--stats", type=Path, default=None, help="optional stats.json from `measure`")
    ph.add_argument("--host", default="wilkcomp02")
    ph.add_argument("--date", default="2026-11-27")
    ph.add_argument("--command", default="")
    ph.add_argument("--set-collector", action="store_true")
    ph.add_argument("--check", action="store_true")

    pr = sub.add_parser("render", help="Render saturation figure (.png + .svg)")
    pr.add_argument("--data", type=Path, default=DEFAULT_HISTORY)
    pr.add_argument("--out", required=True, type=Path, help="output path prefix (writes .png and .svg)")

    args = ap.parse_args(argv)

    if args.cmd == "measure":
        measure_release(
            root=args.root,
            release=args.release,
            earlier=args.earlier,
            out=args.out,
            tmp=args.tmp,
            threads=args.threads,
            buckets=args.buckets,
            stats_out=args.stats_out,
        )
        return 0

    if args.cmd == "history":
        if args.add is not None:
            if not args.command.strip():
                raise ValueError("--command is required when --add is given (AGENTS.md §6)")
            sat_doc = json.loads(args.add.read_text(encoding="utf-8"))
            st_doc = json.loads(args.stats.read_text(encoding="utf-8")) if args.stats and args.stats.is_file() else None
            add_to_history(
                args.data,
                sat_doc,
                host=args.host,
                date_str=args.date,
                command=args.command,
                diversity_path=args.diversity,
                set_collector=args.set_collector,
                stats=st_doc,
            )
        if args.check or args.add is None:
            if not args.data.is_file():
                print(f"problem: missing history file {args.data}")
                return 1
            hist = json.loads(args.data.read_text(encoding="utf-8"))
            probs = check_history(hist)
            for p in probs:
                print(f"problem: {p}")
            print(f"checks: {'all pass' if not probs else f'{len(probs)} problems'}")
            return 0 if not probs else 1
        return 0

    if args.cmd == "render":
        hist = json.loads(args.data.read_text(encoding="utf-8"))
        probs = check_history(hist)
        if probs:
            for p in probs:
                print(f"problem: {p}", file=sys.stderr)
            return 1
        png, svg = render_figure(hist, args.out)
        print(f"wrote {png} and {svg}")
        return 0

    return 2


if __name__ == "__main__":
    sys.exit(main())
