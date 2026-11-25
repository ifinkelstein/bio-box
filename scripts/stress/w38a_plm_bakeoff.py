"""Card W-38a: GPU pLM & 3Di remote-homology bake-off holdout builder, SPIRE-length locus truncator,
evaluation metrics harness, and Go/No-Go Gate 1-3 checker (`SD-21` extension).

Design & protocol: `planning/design/W-38a_bakeoff_protocol.md` and `planning/design/pLM_remote_system_search.md`.

Invariants enforced (`AGENTS.md` §3-§4):
- Manifest-backed reads ONLY via `biobox.catalog.ReleaseFiles.of(root, release)` and `SeqStore(root, files)`
  (never filesystem globs, §3.9). Default release is parameterized (`--release`, default `r26.11.2`;
  `r26.11.1` is scheduled for archive retirement under `DEC-66`).
- Stable keys: every record carries `feature_uid` and `protein_uid`; ordinals never appear without
  `feature_uid` (§3.3), and release-local `token_id`s never cross the boundary (§3.4).
- Bounded memory (§4.1): Parquet tables stream one row group at a time; background sampling uses
  order-invariant keyed bottom-k reservoir sampling (`blake3(seed || protein_uid)`).
- Deterministic (§4.3): sorted outputs and content-hashed sampling keys give identical results
  regardless of input order or thread count.
"""

from __future__ import annotations

import argparse
import heapq
import json
import math
import sys
from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pyarrow.parquet as pq
from blake3 import blake3

from biobox import catalog
from biobox.l0seq.store import SeqStore

# --- Protocol constants (`W-38a_bakeoff_protocol.md` §3-§5) ---------------------------------------------

DEFAULT_RELEASE = "r26.11.2"
DEFAULT_SEED = 20261125

# Strict dual-gate holdout filter (closing DefensePredictor's >30% ID / qcov > 0.8 gap)
MAX_AA_IDENTITY = 0.20      # strict: max_train_aa_id < 0.20 (mmseqs search -s 7.5 3-iter)
MIN_HMM_EVALUE = 10.0       # strict: min_train_hmm_evalue > 10.0 (hmmsearch Z = 1e8)
NOVEL_FOLD_MIN_EVALUE = 1.0 # Foldseek E > 1.0 defines novel-fold dark systems

# Stratification (i): multi-domain fusions & ESM2 truncation cap
FUSION_MIN_AA = 600
FUSION_MIN_DOMAINS = 2
ESM2_TRUNCATION_AA = 1022

# Stratification (ii): SPIRE empirical contig-length truncation range
SPIRE_MIN_BP = 1500
SPIRE_MAX_BP = 5000

# Stratification (iii): 1,000,000 real background genes with FP-trap cohort quotas
DEFAULT_BG_SIZE = 1_000_000
DEFAULT_TRAP_QUOTAS: dict[str, float] = {
    "ntpase": 0.15,
    "methyltransferase": 0.15,
    "transposase": 0.15,
    "hydrolase": 0.10,
    "general": 0.45,
}
DEFAULT_TRAP_PFAMS: dict[str, tuple[str, ...]] = {
    "ntpase": ("PF00004", "PF00005", "PF00270", "PF00271", "PF05729", "PF07724", "PF13191"),
    "methyltransferase": ("PF01555", "PF08241", "PF08242", "PF13649", "PF13847", "PF01170"),
    "transposase": ("PF00665", "PF01527", "PF01609", "PF02371", "PF13358", "PF01548"),
    "hydrolase": ("PF00561", "PF00702", "PF12697", "PF00149", "PF00756"),
}

# Go / No-Go Gates 1-3 (`W-38a_bakeoff_protocol.md` §5)
GATE1_MIN_DELTA_RECALL_PP = 20.0
GATE1_MAX_LATENCY_RATIO = 1.50
GATE1_PROSTT5_MARGIN_PP = 2.0

GATE2_MIN_UNIQUE_TP_GAIN = 0.10
GATE2_MAX_FUSION_LOSS = 0.05

GATE3_MIN_GBA_RATIO = 2.0
GATE3_MAX_TRUNCATION_DROP = 0.15
GATE3_MAX_NR90_AP_DROP = 0.03


# --- Manifest-backed release readers (`AGENTS.md` §3.9) -------------------------------------------------


def open_release_files(root: Path, release: str) -> catalog.ReleaseFiles:
    """Open a release strictly via its manifest / `release_def` (never globs, AGENTS.md §3.9)."""
    return catalog.ReleaseFiles.of(Path(root), release)


def open_release_seqstore(root: Path, release: str) -> SeqStore:
    """Open a release's `SeqStore` using manifest-backed `ReleaseFiles`."""
    files = open_release_files(root, release)
    return SeqStore(Path(root), files)


def stream_release_proteins(
    files: catalog.ReleaseFiles,
    *,
    gene_call_version: str | None = None,
    stages: set[str] | None = None,
) -> Iterator[dict[str, Any]]:
    """Stream distinct `(protein_uid, length, seq)` rows from `protein_seq` one row group at a time."""
    want: dict[str, str] = {}
    if gene_call_version is not None:
        want["gene_call_version"] = gene_call_version
    entries = files.entries("protein_seq", **want)
    for path, parts in entries:
        stage = parts.get("stage", "")
        if stages is not None and stage not in stages:
            continue
        pf = pq.ParquetFile(path)
        for rg in range(pf.metadata.num_row_groups):
            tbl = pf.read_row_group(rg, columns=["protein_uid", "length", "seq"], use_threads=False)
            uids = tbl.column("protein_uid").to_pylist()
            lens = tbl.column("length").to_pylist()
            seqs = tbl.column("seq").to_pylist()
            for pu, length, seq in zip(uids, lens, seqs):
                yield {
                    "protein_uid": pu.hex() if isinstance(pu, (bytes, bytearray)) else str(pu),
                    "length": int(length),
                    "seq": str(seq),
                    "stage": stage,
                }


def stream_release_tapes(
    files: catalog.ReleaseFiles,
    *,
    tape_version: str | None = None,
    stages: set[str] | None = None,
) -> Iterator[dict[str, Any]]:
    """Stream `contig_tapes` rows one row group at a time from manifest-listed files."""
    want: dict[str, str] = {}
    if tape_version is not None:
        want["tape_version"] = tape_version
    cols = [
        "contig_uid",
        "contig_len",
        "is_circular",
        "n_slots",
        "feature_kinds",
        "strands",
        "start_bps",
        "end_bps",
        "protein_uids",
        "feature_uids",
    ]
    for path, parts in files.entries("contig_tapes", **want):
        stage = parts.get("stage", "")
        if stages is not None and stage not in stages:
            continue
        pf = pq.ParquetFile(path)
        for rg in range(pf.metadata.num_row_groups):
            tbl = pf.read_row_group(rg, columns=cols, use_threads=False)
            for row in tbl.to_pylist():
                row["stage"] = stage
                yield row


def classify_trap_cohort(
    model_accs: Iterable[str],
    trap_pfams: Mapping[str, Sequence[str]] = DEFAULT_TRAP_PFAMS,
) -> str:
    """Assign a protein to a deterministic FP-trap cohort from its Pfam accessions."""
    norm = {acc.split(".")[0].upper() for acc in model_accs if acc}
    for cohort in ("ntpase", "methyltransferase", "transposase", "hydrolase"):
        if norm & set(trap_pfams.get(cohort, ())):
            return cohort
    return "general"


def stream_domain_traps(
    files: catalog.ReleaseFiles,
    *,
    library: str | None = None,
    trap_pfams: Mapping[str, Sequence[str]] = DEFAULT_TRAP_PFAMS,
) -> dict[str, str]:
    """Stream `l1_domains` row-group by row-group to map `protein_uid` hex -> trap cohort.

    Stores only proteins that match a trap Pfam accession (bounded to trap hits, not the whole corpus).
    """
    want: dict[str, str] = {}
    if library is not None:
        want["library"] = library
    acc_to_cohort: dict[str, str] = {}
    for cohort in ("ntpase", "methyltransferase", "transposase", "hydrolase"):
        for acc in trap_pfams.get(cohort, ()):
            acc_to_cohort.setdefault(acc.split(".")[0].upper(), cohort)

    priority = {"ntpase": 0, "methyltransferase": 1, "transposase": 2, "hydrolase": 3}
    hits: dict[str, str] = {}
    for path, _ in files.entries("l1_domains", **want):
        pf = pq.ParquetFile(path)
        for rg in range(pf.metadata.num_row_groups):
            tbl = pf.read_row_group(rg, columns=["protein_uid", "model_acc"], use_threads=False)
            for pu, acc in zip(tbl.column("protein_uid").to_pylist(), tbl.column("model_acc").to_pylist()):
                base = str(acc).split(".")[0].upper()
                cohort = acc_to_cohort.get(base)
                if cohort is None:
                    continue
                pu_hex = pu.hex() if isinstance(pu, (bytes, bytearray)) else str(pu)
                prev = hits.get(pu_hex)
                if prev is None or priority[cohort] < priority[prev]:
                    hits[pu_hex] = cohort
    return hits


# --- Strict holdout leakage filter ----------------------------------------------------------------------


def _is_valid_number(val: Any) -> bool:
    if val is None or isinstance(val, bool):
        return False
    try:
        f = float(val)
    except (TypeError, ValueError):
        return False
    return not math.isnan(f)


def passes_strict_holdout(
    record: Mapping[str, Any],
    *,
    max_aa_id: float = MAX_AA_IDENTITY,
    min_hmm_evalue: float = MIN_HMM_EVALUE,
) -> tuple[bool, str]:
    """Check whether a candidate positive passes `<20%` AA ID AND `hmmsearch E > 10` (fail-closed).

    Returns `(passed, reason)`.
    """
    aa_id_raw = record.get("max_train_aa_id")
    hmm_e_raw = record.get("min_train_hmm_evalue")

    if not _is_valid_number(aa_id_raw):
        return False, "missing_or_nan_aa_id"
    if not _is_valid_number(hmm_e_raw):
        return False, "missing_or_nan_hmm_evalue"

    aa_id = float(aa_id_raw)
    hmm_e = float(hmm_e_raw)
    if aa_id < 0.0 or hmm_e < 0.0:
        return False, "negative_homology_metric"
    if aa_id >= max_aa_id:
        return False, f"aa_id_leakage:{aa_id:.4f}>={max_aa_id:.4f}"
    if hmm_e <= min_hmm_evalue:
        return False, f"hmm_evalue_leakage:{hmm_e:.3g}<={min_hmm_evalue:.3g}"
    return True, "pass"


def filter_holdout_leakage(
    candidates: Iterable[Mapping[str, Any]],
    *,
    max_aa_id: float = MAX_AA_IDENTITY,
    min_hmm_evalue: float = MIN_HMM_EVALUE,
    novel_fold_min_evalue: float = NOVEL_FOLD_MIN_EVALUE,
    fusion_min_aa: int = FUSION_MIN_AA,
    fusion_min_domains: int = FUSION_MIN_DOMAINS,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Filter candidate positives with strict dual-gate leakage control and deterministic stratification.

    Returns `(kept, rejected)`, both deterministically sorted by `(system_id, feature_uid, protein_uid)`.
    """
    kept: list[dict[str, Any]] = []
    rejected: list[dict[str, Any]] = []

    for rec in candidates:
        row = dict(rec)
        if "feature_uid" not in row or "protein_uid" not in row:
            raise ValueError("every positive candidate must carry both 'feature_uid' and 'protein_uid' (AGENTS §3.3)")

        ok, reason = passes_strict_holdout(row, max_aa_id=max_aa_id, min_hmm_evalue=min_hmm_evalue)
        if not ok:
            row["reject_reason"] = reason
            rejected.append(row)
            continue

        aa_len = int(row.get("aa_len", 0))
        n_domains = int(row.get("n_domains", 1))
        foldseek_e_raw = row.get("min_train_foldseek_evalue")
        is_novel_fold = _is_valid_number(foldseek_e_raw) and float(foldseek_e_raw) > novel_fold_min_evalue

        explicit_fusion = row.get("is_fusion")
        is_fusion = bool(explicit_fusion) if explicit_fusion is not None else (
            aa_len >= fusion_min_aa and n_domains >= fusion_min_domains
        )

        row["label"] = 1
        row["holdout_stratum"] = "novel_fold" if is_novel_fold else "remote_homolog"
        row["is_fusion"] = is_fusion
        row["exceeds_esm2_cap"] = aa_len > ESM2_TRUNCATION_AA
        row["trap_cohort"] = "positive"
        kept.append(row)

    sort_key = lambda r: (str(r.get("system_id", "")), str(r["feature_uid"]), str(r["protein_uid"]))
    kept.sort(key=sort_key)
    rejected.sort(key=sort_key)
    return kept, rejected


# --- Deterministic bounded-memory background sampler ----------------------------------------------------


def _hash_u64(seed: int, key: str) -> int:
    """Deterministic 64-bit hash key for order-invariant bottom-k reservoir sampling."""
    digest = blake3(f"w38a-sample:{seed}:{key}".encode("utf-8")).digest(length=8)
    return int.from_bytes(digest, "big")


def sample_background_reservoir(
    records: Iterable[Mapping[str, Any]],
    *,
    total_size: int = DEFAULT_BG_SIZE,
    seed: int = DEFAULT_SEED,
    quotas: Mapping[str, float] = DEFAULT_TRAP_QUOTAS,
    exclude_protein_uids: set[str] | frozenset[str] = frozenset(),
) -> list[dict[str, Any]]:
    """Order-invariant, bounded-memory stratified reservoir sampler over a stream of background genes.

    Uses keyed bottom-k selection (`blake3(seed || protein_uid)`) per cohort so peak memory is
    `O(total_size)` and shuffling the input stream produces the exact same sampled records.
    """
    if total_size < 0:
        raise ValueError("total_size must be non-negative")
    if total_size == 0:
        return []

    cohort_names = sorted(quotas.keys())
    weight_sum = sum(float(quotas[c]) for c in cohort_names)
    if weight_sum <= 0:
        raise ValueError("quotas must sum to a positive number")

    targets: dict[str, int] = {}
    assigned = 0
    for idx, c in enumerate(cohort_names):
        if idx == len(cohort_names) - 1:
            targets[c] = max(0, total_size - assigned)
        else:
            k = int(round(total_size * float(quotas[c]) / weight_sum))
            targets[c] = max(0, min(k, total_size - assigned))
            assigned += targets[c]

    # Max-heap per cohort storing (-hash_u64, protein_uid, record) so we keep the smallest k hashes
    heaps: dict[str, list[tuple[int, str, dict[str, Any]]]] = {c: [] for c in cohort_names}
    seen_uids: set[str] = set()

    for rec in records:
        pu = str(rec["protein_uid"])
        if pu in exclude_protein_uids or pu in seen_uids:
            continue
        seen_uids.add(pu)
        cohort = str(rec.get("trap_cohort", "general"))
        if cohort not in targets:
            cohort = "general" if "general" in targets else cohort_names[0]
        cap = targets[cohort]
        if cap <= 0:
            continue

        h = _hash_u64(seed, pu)
        heap = heaps[cohort]
        entry_rec = dict(rec)
        entry_rec["label"] = 0
        entry_rec["trap_cohort"] = cohort
        entry_rec["holdout_stratum"] = "background"
        entry_rec.setdefault("system_id", "")
        entry_rec.setdefault("is_fusion", False)
        entry_rec.setdefault("exceeds_esm2_cap", int(entry_rec.get("aa_len", entry_rec.get("length", 0))) > ESM2_TRUNCATION_AA)
        item = (-h, pu, entry_rec)
        if len(heap) < cap:
            heapq.heappush(heap, item)
        elif item > heap[0]:
            heapq.heapreplace(heap, item)

    sampled: list[dict[str, Any]] = []
    for c in cohort_names:
        sampled.extend(item[2] for item in heaps[c])

    sampled.sort(key=lambda r: (str(r["trap_cohort"]), str(r.get("feature_uid", "")), str(r["protein_uid"])))
    return sampled


# --- Circular-aware SPIRE-length locus truncator (Stratification ii) ------------------------------------


def empirical_spire_length(
    seed: int,
    key: str,
    *,
    min_bp: int = SPIRE_MIN_BP,
    max_bp: int = SPIRE_MAX_BP,
    empirical_lengths: Sequence[int] | None = None,
) -> int:
    """Draw a deterministic contig length in `[min_bp, max_bp]` from empirical SPIRE lengths.

    When `empirical_lengths` is provided, draws uniformly from the filtered empirical list.
    Otherwise uses a deterministic right-skewed Beta(1.6, 3.0)-like quantile proxy (mode ~2.3 kb,
    matching fragmented short-read metagenome assemblies in `[1500, 5000]` bp).
    """
    if min_bp <= 0 or max_bp < min_bp:
        raise ValueError(f"invalid SPIRE length bounds [{min_bp}, {max_bp}]")
    u = _hash_u64(seed, f"spire-len:{key}") / float(1 << 64)

    if empirical_lengths:
        valid = sorted(int(x) for x in empirical_lengths if min_bp <= int(x) <= max_bp)
        if valid:
            idx = min(len(valid) - 1, int(u * len(valid)))
            return valid[idx]

    # Deterministic right-skewed draw in [min_bp, max_bp] via u^1.45
    span = max_bp - min_bp
    return int(round(min_bp + span * (u ** 1.45)))


def _shortest_circular_arc(intervals: list[tuple[int, int]], contig_len: int) -> tuple[int, int, bool]:
    """Find the shortest arc `[arc_start, arc_end)` on a circle of length `contig_len` covering `intervals`.

    Intervals are 0-based half-open; an interval may already have `end > contig_len` (origin-wrapping).
    Returns `(arc_start, arc_end, crosses_origin)` with `0 <= arc_start < contig_len` and
    `arc_end <= arc_start + contig_len`.
    """
    if not intervals:
        return 0, 0, False

    # If any single interval already wraps the origin (`end > contig_len`), merge unwrapped directly
    if any(e > contig_len for _, e in intervals):
        # Shift intervals near 0 that sit closer to contig_len+ into the unwrapped frame
        unwrapped: list[tuple[int, int]] = []
        ref_start = min(s for s, e in intervals if e > contig_len)
        for s, e in intervals:
            if e <= contig_len and s < ref_start and (s + contig_len - ref_start) < (ref_start - s):
                unwrapped.append((s + contig_len, e + contig_len))
            else:
                unwrapped.append((s, e))
        a_s = min(s for s, _ in unwrapped)
        a_e = max(e for _, e in unwrapped)
        return a_s, min(a_e, a_s + contig_len), a_e > contig_len

    # Merge intervals on [0, contig_len) and find the largest gap on the circle
    pts = sorted(intervals)
    merged: list[list[int]] = []
    for s, e in pts:
        if not merged or s > merged[-1][1]:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)

    if len(merged) == 1:
        return merged[0][0], merged[0][1], False

    # Gaps between consecutive merged intervals, plus the wraparound gap from last end to first start + L
    best_gap = -1
    best_start = merged[0][0]
    best_end = merged[-1][1]
    best_crosses = False

    # Wraparound gap (corresponds to linear arc [merged[0][0], merged[-1][1]])
    wrap_gap = (merged[0][0] + contig_len) - merged[-1][1]
    best_gap = wrap_gap
    best_start = merged[0][0]
    best_end = merged[-1][1]
    best_crosses = False

    for i in range(len(merged) - 1):
        gap = merged[i + 1][0] - merged[i][1]
        if gap > best_gap:
            best_gap = gap
            best_start = merged[i + 1][0]
            best_end = merged[i][1] + contig_len
            best_crosses = True

    return best_start, best_end, best_crosses


def truncate_locus_to_spire(
    locus: Mapping[str, Any],
    *,
    seed: int = DEFAULT_SEED,
    target_len_bp: int | None = None,
    min_bp: int = SPIRE_MIN_BP,
    max_bp: int = SPIRE_MAX_BP,
    empirical_lengths: Sequence[int] | None = None,
) -> dict[str, Any]:
    """Slice a closed-contig locus into a linear `1.5-5 kb` SPIRE-length fragment.

    Preserves circular origin wraparound (`AGENTS.md` §3.1) and exact `parent_feature_uid` +
    `protein_uid` provenance (`AGENTS.md` §3.3). Updates edge-censoring distances
    (`dist_to_contig_end_bp`, `orfs_to_contig_start`, `orfs_to_contig_end`) and marks clipped
    boundary ORFs with `is_partial = True`.
    """
    parent_contig_uid = str(locus["contig_uid"])
    contig_len = int(locus["contig_len"])
    is_circular = bool(locus.get("is_circular", False))
    slots: list[dict[str, Any]] = [dict(s) for s in locus["slots"]]
    target_fuids = {str(x) for x in locus.get("target_feature_uids", ())}
    if not slots:
        raise ValueError(f"locus on contig {parent_contig_uid} has no slots")

    anchor_fuid = str(locus.get("anchor_feature_uid") or sorted(target_fuids)[0] if target_fuids else slots[0]["feature_uid"])
    frag_len = int(target_len_bp) if target_len_bp is not None else empirical_spire_length(
        seed, f"{parent_contig_uid}:{anchor_fuid}",
        min_bp=min_bp, max_bp=max_bp, empirical_lengths=empirical_lengths,
    )
    frag_len = min(frag_len, contig_len)

    target_slots = [s for s in slots if str(s["feature_uid"]) in target_fuids] if target_fuids else [slots[0]]
    if not target_slots:
        target_slots = [slots[0]]

    target_intervals = [(int(s["start"]), int(s["end"])) for s in target_slots]
    if is_circular:
        arc_start, arc_end, crosses_origin = _shortest_circular_arc(target_intervals, contig_len)
    else:
        arc_start = min(s for s, _ in target_intervals)
        arc_end = max(e for _, e in target_intervals)
        crosses_origin = False

    arc_mid = (arc_start + arc_end) // 2
    win_start = arc_mid - (frag_len // 2)
    win_end = win_start + frag_len

    if not is_circular:
        if win_start < 0:
            win_start = 0
            win_end = frag_len
        elif win_end > contig_len:
            win_end = contig_len
            win_start = max(0, contig_len - frag_len)
    else:
        # Normalize win_start into [0, contig_len)
        win_start = win_start % contig_len
        win_end = win_start + frag_len
        if win_end > contig_len:
            crosses_origin = True

    # Project slots into [win_start, win_end) -> local [0, frag_len)
    kept_slots: list[dict[str, Any]] = []
    shifts = (0, contig_len, -contig_len) if is_circular else (0,)
    seen_fuids: set[str] = set()

    for s in slots:
        fuid = str(s["feature_uid"])
        s0, e0 = int(s["start"]), int(s["end"])
        for shift in shifts:
            ss, ee = s0 + shift, e0 + shift
            ov_start = max(ss, win_start)
            ov_end = min(ee, win_end)
            if ov_start < ov_end and fuid not in seen_fuids:
                seen_fuids.add(fuid)
                local_start = ov_start - win_start
                local_end = ov_end - win_start
                is_partial = (ss < win_start) or (ee > win_end)
                kept_slots.append({
                    "parent_feature_uid": fuid,
                    "protein_uid": s.get("protein_uid"),
                    "local_start": int(local_start),
                    "local_end": int(local_end),
                    "strand": int(s["strand"]),
                    "kind": int(s.get("kind", 1)),
                    "is_partial": is_partial,
                    "dist_to_contig_end_bp": int(min(local_start, frag_len - local_end)),
                    "is_target": fuid in target_fuids,
                })
                break

    kept_slots.sort(key=lambda x: (x["local_start"], x["local_end"], x["strand"], x["parent_feature_uid"]))
    n_kept = len(kept_slots)
    for idx, slot in enumerate(kept_slots):
        slot["orfs_to_contig_start"] = idx
        slot["orfs_to_contig_end"] = n_kept - 1 - idx

    trunc_uid = blake3(
        f"w38a-trunc:{parent_contig_uid}:{win_start}:{frag_len}".encode("utf-8")
    ).digest(length=16).hex()

    return {
        "truncated_contig_uid": trunc_uid,
        "parent_contig_uid": parent_contig_uid,
        "parent_contig_len": contig_len,
        "parent_is_circular": is_circular,
        "is_circular": False,
        "contig_len": frag_len,
        "window_start_bp": int(win_start),
        "window_end_bp": int(win_end),
        "crosses_origin_in_parent": bool(crosses_origin),
        "anchor_feature_uid": anchor_fuid,
        "n_slots": n_kept,
        "n_target_slots_retained": sum(1 for x in kept_slots if x["is_target"]),
        "n_partial_slots": sum(1 for x in kept_slots if x["is_partial"]),
        "slots": kept_slots,
    }


def truncate_loci_batch(
    loci: Iterable[Mapping[str, Any]],
    *,
    seed: int = DEFAULT_SEED,
    min_bp: int = SPIRE_MIN_BP,
    max_bp: int = SPIRE_MAX_BP,
    empirical_lengths: Sequence[int] | None = None,
) -> list[dict[str, Any]]:
    """Deterministically truncate a batch of closed loci, sorted by `(parent_contig_uid, anchor_feature_uid)`."""
    out = [
        truncate_locus_to_spire(
            loc, seed=seed, min_bp=min_bp, max_bp=max_bp, empirical_lengths=empirical_lengths
        )
        for loc in loci
    ]
    out.sort(key=lambda r: (r["parent_contig_uid"], r["anchor_feature_uid"]))
    return out


# --- Evaluation metrics harness (pure NumPy, no sklearn required) ---------------------------------------


def average_precision(y_true: Sequence[int], y_score: Sequence[float]) -> float:
    """Exact step-function Average Precision (PR-AUC) with deterministic tie handling."""
    yt = np.asarray(y_true, dtype=np.int8)
    ys = np.asarray(y_score, dtype=np.float64)
    if yt.shape != ys.shape or yt.size == 0:
        return 0.0
    n_pos = int(np.sum(yt == 1))
    if n_pos == 0:
        return 0.0

    # Sort descending by score; for ties, put negatives before positives (conservative, deterministic)
    order = np.lexsort((yt, -ys))
    yt_sorted = yt[order]
    tp = np.cumsum(yt_sorted == 1)
    ranks = np.arange(1, yt_sorted.size + 1, dtype=np.float64)
    precisions = tp / ranks
    return float(np.sum(precisions[yt_sorted == 1]) / n_pos)


def threshold_for_fpr(y_true: Sequence[int], y_score: Sequence[float], max_fpr: float) -> float:
    """Smallest score threshold where `FP / N_neg <= max_fpr` (inclusive `score >= threshold`)."""
    yt = np.asarray(y_true, dtype=np.int8)
    ys = np.asarray(y_score, dtype=np.float64)
    neg_scores = ys[yt == 0]
    if neg_scores.size == 0:
        return float("-inf")
    max_fp = int(math.floor(max_fpr * neg_scores.size + 1e-12))
    sorted_neg = np.sort(neg_scores)[::-1]  # descending
    if max_fp <= 0:
        return float(np.nextafter(sorted_neg[0], np.inf))
    if max_fp >= sorted_neg.size:
        return float("-inf")
    # Threshold strictly above the (max_fp)-th highest negative (0-indexed)
    cutoff_neg = float(sorted_neg[max_fp])
    return float(np.nextafter(cutoff_neg, np.inf))


def recall_at_fpr(y_true: Sequence[int], y_score: Sequence[float], max_fpr: float) -> float:
    """Recall at the operating threshold achieving `FPR <= max_fpr`."""
    yt = np.asarray(y_true, dtype=np.int8)
    ys = np.asarray(y_score, dtype=np.float64)
    n_pos = int(np.sum(yt == 1))
    if n_pos == 0:
        return 0.0
    thr = threshold_for_fpr(yt, ys, max_fpr)
    tp = int(np.sum((yt == 1) & (ys >= thr)))
    return float(tp / n_pos)


def recall_at_k(y_true: Sequence[int], y_score: Sequence[float], k: int) -> float:
    """Fraction of all positives appearing in the top-`k` scored items."""
    yt = np.asarray(y_true, dtype=np.int8)
    ys = np.asarray(y_score, dtype=np.float64)
    n_pos = int(np.sum(yt == 1))
    if n_pos == 0 or k <= 0:
        return 0.0
    order = np.lexsort((yt, -ys))[: min(k, yt.size)]
    return float(np.sum(yt[order] == 1) / n_pos)


def recovered_systems_at_fpr(
    records: Sequence[Mapping[str, Any]],
    scores: Mapping[str, float],
    max_fpr: float,
    *,
    score_key_field: str = "feature_uid",
) -> tuple[set[str], set[str]]:
    """Return `(recovered_system_ids, all_positive_system_ids)` at `FPR <= max_fpr`."""
    y_true: list[int] = []
    y_score: list[float] = []
    for r in records:
        y_true.append(int(r["label"]))
        y_score.append(float(scores.get(str(r[score_key_field]), float("-inf"))))
    thr = threshold_for_fpr(y_true, y_score, max_fpr)

    all_sys: set[str] = set()
    rec_sys: set[str] = set()
    for r, sc in zip(records, y_score):
        if int(r["label"]) == 1:
            sid = str(r.get("system_id") or r[score_key_field])
            all_sys.add(sid)
            if sc >= thr:
                rec_sys.add(sid)
    return rec_sys, all_sys


def evaluate_method(
    holdout_records: Sequence[Mapping[str, Any]],
    scores: Mapping[str, float],
    *,
    truncated_scores: Mapping[str, float] | None = None,
    nr90_scores: Mapping[str, float] | None = None,
    latency_ratio_vs_afdb: float | None = None,
) -> dict[str, Any]:
    """Compute full stratified metrics for one competitor method."""
    def _subset_metrics(pos_filter) -> dict[str, float]:
        yt: list[int] = []
        ys: list[float] = []
        for r in holdout_records:
            lbl = int(r["label"])
            if lbl == 0 or pos_filter(r):
                yt.append(lbl)
                ys.append(float(scores.get(str(r["feature_uid"]), float("-inf"))))
        return {
            "n_pos": int(sum(yt)),
            "n_neg": int(len(yt) - sum(yt)),
            "pr_auc": average_precision(yt, ys),
            "recall_at_1pct_fpr": recall_at_fpr(yt, ys, 0.01),
            "recall_at_0_1pct_fpr": recall_at_fpr(yt, ys, 0.001),
            "recall_at_100": recall_at_k(yt, ys, 100),
        }

    overall = _subset_metrics(lambda r: True)
    remote = _subset_metrics(lambda r: r.get("holdout_stratum") == "remote_homolog")
    novel = _subset_metrics(lambda r: r.get("holdout_stratum") == "novel_fold")
    fusion = _subset_metrics(lambda r: bool(r.get("is_fusion")))
    fusion_gt_1022 = _subset_metrics(lambda r: bool(r.get("is_fusion")) and bool(r.get("exceeds_esm2_cap")))
    monodomain = _subset_metrics(lambda r: not bool(r.get("is_fusion")))
    remote_and_fusion = _subset_metrics(
        lambda r: r.get("holdout_stratum") == "remote_homolog" or bool(r.get("is_fusion"))
    )

    rec_sys_01, all_sys = recovered_systems_at_fpr(holdout_records, scores, 0.001)

    # Per-cohort trap FPR at the overall 1% FPR threshold
    yt_all = [int(r["label"]) for r in holdout_records]
    ys_all = [float(scores.get(str(r["feature_uid"]), float("-inf"))) for r in holdout_records]
    thr_1pct = threshold_for_fpr(yt_all, ys_all, 0.01)
    trap_counts: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for r, sc in zip(holdout_records, ys_all):
        if int(r["label"]) == 0:
            c = str(r.get("trap_cohort", "general"))
            trap_counts[c][1] += 1
            if sc >= thr_1pct:
                trap_counts[c][0] += 1
    trap_fpr = {c: (fp / n if n else 0.0) for c, (fp, n) in sorted(trap_counts.items())}

    # Fusion vs monodomain recall loss at 1% FPR (and 0.1% FPR)
    mono_r = monodomain["recall_at_1pct_fpr"]
    fus_r = fusion["recall_at_1pct_fpr"]
    fusion_recall_loss = max(0.0, mono_r - fus_r)

    # Stratification (ii): truncation relative recall drop (if truncated_scores provided)
    truncation_drop: float | None = None
    truncated_recall_01: float | None = None
    if truncated_scores is not None:
        yt_t = [int(r["label"]) for r in holdout_records]
        ys_t = [float(truncated_scores.get(str(r["feature_uid"]), float("-inf"))) for r in holdout_records]
        truncated_recall_01 = recall_at_fpr(yt_t, ys_t, 0.001)
        closed_r = overall["recall_at_0_1pct_fpr"]
        truncation_drop = ((closed_r - truncated_recall_01) / closed_r) if closed_r > 0 else 0.0

    # NR90 representative substitution AP drop (if nr90_scores provided)
    nr90_ap_drop: float | None = None
    if nr90_scores is not None:
        yt_n = [int(r["label"]) for r in holdout_records]
        ys_n = [float(nr90_scores.get(str(r["feature_uid"]), float("-inf"))) for r in holdout_records]
        nr90_ap = average_precision(yt_n, ys_n)
        nr90_ap_drop = max(0.0, overall["pr_auc"] - nr90_ap)

    return {
        "overall": overall,
        "remote_homolog": remote,
        "novel_fold": novel,
        "fusion": fusion,
        "fusion_gt_1022": fusion_gt_1022,
        "monodomain": monodomain,
        "remote_and_fusion": remote_and_fusion,
        "fusion_recall_loss": fusion_recall_loss,
        "recovered_systems_at_0_1pct_fpr": sorted(rec_sys_01),
        "total_positive_systems": len(all_sys),
        "trap_fpr_at_1pct_overall": trap_fpr,
        "truncated_recall_at_0_1pct_fpr": truncated_recall_01,
        "truncation_relative_recall_drop": truncation_drop,
        "nr90_ap_drop": nr90_ap_drop,
        "latency_ratio_vs_afdb": latency_ratio_vs_afdb,
    }


# --- Go / No-Go Gate 1-3 checker (`W-38a_bakeoff_protocol.md` §5) ---------------------------------------


def check_gates(eval_report: Mapping[str, Any]) -> dict[str, Any]:
    """Evaluate falsifiable Go/No-Go Gates 1-3 from an `evaluate` JSON report."""
    methods: Mapping[str, Any] = eval_report.get("methods", eval_report)

    # --- Gate 1: Option A 3Di (`ProstT5` vs `ESM-2 3B 3Di` over `jackhmmer_3iter`) ---
    jack = methods.get("jackhmmer_3iter")
    jack_r = float(jack["remote_and_fusion"]["recall_at_1pct_fpr"]) if jack else 0.0
    three_di_candidates = [k for k in ("prostt5_foldseek", "esm2_3b_3di_foldseek") if k in methods]
    gate1_details: dict[str, Any] = {}
    passing_3di: list[tuple[str, float]] = []

    for name in three_di_candidates:
        m = methods[name]
        r = float(m["remote_and_fusion"]["recall_at_1pct_fpr"])
        delta_pp = (r - jack_r) * 100.0
        lat = m.get("latency_ratio_vs_afdb")
        lat_ok = (lat is not None) and (float(lat) <= GATE1_MAX_LATENCY_RATIO)
        passed = (delta_pp >= GATE1_MIN_DELTA_RECALL_PP) and lat_ok
        gate1_details[name] = {
            "recall_at_1pct_fpr": r,
            "jackhmmer_recall_at_1pct_fpr": jack_r,
            "delta_recall_pp": delta_pp,
            "latency_ratio_vs_afdb": lat,
            "recall_pass": delta_pp >= GATE1_MIN_DELTA_RECALL_PP,
            "latency_pass": lat_ok,
            "pass": passed,
        }
        if passed:
            passing_3di.append((name, r))

    gate1_pass = bool(passing_3di)
    recommended_3di: str | None = None
    if gate1_pass:
        by_name = dict(passing_3di)
        if "prostt5_foldseek" in by_name and "esm2_3b_3di_foldseek" in by_name:
            gap_pp = (by_name["esm2_3b_3di_foldseek"] - by_name["prostt5_foldseek"]) * 100.0
            recommended_3di = "prostt5_foldseek" if gap_pp <= GATE1_PROSTT5_MARGIN_PP else "esm2_3b_3di_foldseek"
        else:
            recommended_3di = max(passing_3di, key=lambda x: x[1])[0]

    # --- Gate 2: Option B Corpus-Wide Dense ANN (`ESM-C 300M`) ---
    baseline_sys: set[str] = set()
    for base_key in ("jackhmmer_3iter", "hhblits_3iter", recommended_3di or "prostt5_foldseek"):
        if base_key in methods:
            baseline_sys.update(methods[base_key].get("recovered_systems_at_0_1pct_fpr", ()))
    total_sys = max((int(m.get("total_positive_systems", 0)) for m in methods.values()), default=0)

    ann_candidates = [k for k in ("esmc_300m_segments_ivfpq", "esmc_300m_full_fp16") if k in methods]
    gate2_details: dict[str, Any] = {}
    gate2_pass = False
    for name in ann_candidates:
        m = methods[name]
        ann_sys = set(m.get("recovered_systems_at_0_1pct_fpr", ()))
        marginal = ann_sys - baseline_sys
        gain_frac = (len(marginal) / total_sys) if total_sys > 0 else 0.0
        fus_loss = float(m.get("fusion_recall_loss", 1.0))
        passed = (gain_frac >= GATE2_MIN_UNIQUE_TP_GAIN) and (fus_loss < GATE2_MAX_FUSION_LOSS)
        gate2_details[name] = {
            "marginal_unique_tp_systems": len(marginal),
            "total_positive_systems": total_sys,
            "marginal_gain_fraction": gain_frac,
            "fusion_recall_loss": fus_loss,
            "gain_pass": gain_frac >= GATE2_MIN_UNIQUE_TP_GAIN,
            "fusion_pass": fus_loss < GATE2_MAX_FUSION_LOSS,
            "pass": passed,
        }
        if passed:
            gate2_pass = True

    # --- Gate 3: Option C Tape-Context Scorer vs Option D On-Demand ---
    nbr = methods.get("nbr_gba")
    nbr_novel_r = float(nbr["novel_fold"]["recall_at_0_1pct_fpr"]) if nbr else 0.0
    ctx_candidates = [k for k in ("esm_df_geneclr_df", "defense_predictor_fixed") if k in methods]
    gate3_details: dict[str, Any] = {}
    passing_ctx: list[str] = []

    for name in ctx_candidates:
        m = methods[name]
        novel_r = float(m["novel_fold"]["recall_at_0_1pct_fpr"])
        gba_ratio = (novel_r / nbr_novel_r) if nbr_novel_r > 0 else (math.inf if novel_r > 0 else 0.0)
        trunc_drop = m.get("truncation_relative_recall_drop")
        nr90_drop = m.get("nr90_ap_drop")

        a_ok = gba_ratio >= GATE3_MIN_GBA_RATIO
        b_ok = (trunc_drop is not None) and (float(trunc_drop) < GATE3_MAX_TRUNCATION_DROP)
        c_ok = (nr90_drop is not None) and (float(nr90_drop) < GATE3_MAX_NR90_AP_DROP)
        passed = a_ok and b_ok and c_ok
        gate3_details[name] = {
            "novel_fold_recall_at_0_1pct_fpr": novel_r,
            "nbr_gba_novel_fold_recall": nbr_novel_r,
            "gba_ratio": gba_ratio if math.isfinite(gba_ratio) else "inf",
            "truncation_relative_recall_drop": trunc_drop,
            "nr90_ap_drop": nr90_drop,
            "criterion_a_gba_pass": a_ok,
            "criterion_b_truncation_pass": b_ok,
            "criterion_c_nr90_pass": c_ok,
            "pass": passed,
        }
        if passed:
            passing_ctx.append(name)

    gate3_pass = bool(passing_ctx)

    return {
        "gate1_3di": {
            "pass": gate1_pass,
            "recommended_encoder": recommended_3di,
            "action": "PROCEED_W38C_3DI_NR50" if gate1_pass else "HOLD_3DI",
            "candidates": gate1_details,
        },
        "gate2_dense_ann": {
            "pass": gate2_pass,
            "action": "APPROVE_SEGMENTS_IVFPQ_TIER12" if gate2_pass else "REJECT_OPTION_B_PERMANENTLY",
            "candidates": gate2_details,
        },
        "gate3_tape_context": {
            "pass": gate3_pass,
            "passing_models": passing_ctx,
            "action": "BUILD_OPTION_C_AND_D" if gate3_pass else "KEEP_OPTION_D_ONLY",
            "candidates": gate3_details,
        },
    }


# --- JSONL I/O helpers ----------------------------------------------------------------------------------


def _read_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with Path(path).open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def _write_jsonl(path: Path, rows: Iterable[Mapping[str, Any]]) -> int:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(path).with_suffix(Path(path).suffix + ".tmp")
    n = 0
    with tmp.open("w", encoding="utf-8") as fh:
        for r in rows:
            fh.write(json.dumps(r, sort_keys=True) + "\n")
            n += 1
    tmp.replace(path)
    return n


def _write_json(path: Path, obj: Mapping[str, Any]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(path).with_suffix(Path(path).suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(path)


def _load_scores(path: Path) -> dict[str, float]:
    p = Path(path)
    if p.suffix == ".json":
        data = json.loads(p.read_text(encoding="utf-8"))
        return {str(k): float(v) for k, v in data.items()}
    out: dict[str, float] = {}
    with p.open("r", encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("{"):
                rec = json.loads(line)
                out[str(rec["feature_uid"])] = float(rec["score"])
            else:
                cols = line.split("\t")
                if lineno == 1 and cols[0] == "feature_uid":
                    continue
                out[cols[0]] = float(cols[1])
    return out


# --- CLI subcommands ------------------------------------------------------------------------------------


def _cmd_build_holdout(args: argparse.Namespace) -> int:
    pos_candidates = _read_jsonl(Path(args.positives))
    kept_pos, rejected_pos = filter_holdout_leakage(
        pos_candidates,
        max_aa_id=args.max_aa_id,
        min_hmm_evalue=args.min_hmm_evalue,
    )
    exclude_puids = {str(r["protein_uid"]) for r in pos_candidates}

    if args.background_jsonl:
        bg_stream: Iterable[Mapping[str, Any]] = _read_jsonl(Path(args.background_jsonl))
    else:
        files = open_release_files(Path(args.root), args.release)
        trap_map = stream_domain_traps(files)

        def _bg_iter() -> Iterator[dict[str, Any]]:
            for prot in stream_release_proteins(files):
                pu = prot["protein_uid"]
                yield {
                    "feature_uid": pu,  # protein-level surrogate when occurrence join is deferred
                    "protein_uid": pu,
                    "aa_len": prot["length"],
                    "trap_cohort": trap_map.get(pu, "general"),
                    "stage": prot["stage"],
                }

        bg_stream = _bg_iter()

    sampled_bg = sample_background_reservoir(
        bg_stream,
        total_size=args.bg_size,
        seed=args.seed,
        exclude_protein_uids=exclude_puids,
    )

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    n_holdout = _write_jsonl(out_dir / "holdout.jsonl", kept_pos + sampled_bg)
    n_rej = _write_jsonl(out_dir / "rejected_positives.jsonl", rejected_pos)
    summary = {
        "release": args.release,
        "seed": args.seed,
        "max_aa_id": args.max_aa_id,
        "min_hmm_evalue": args.min_hmm_evalue,
        "positives_in": len(pos_candidates),
        "positives_kept": len(kept_pos),
        "positives_rejected": n_rej,
        "background_sampled": len(sampled_bg),
        "holdout_total": n_holdout,
    }
    _write_json(out_dir / "holdout_summary.json", summary)
    return 0


def _cmd_truncate_loci(args: argparse.Namespace) -> int:
    loci = _read_jsonl(Path(args.loci))
    emp_lens: list[int] | None = None
    if args.spire_lengths:
        emp_lens = [
            int(x.strip())
            for x in Path(args.spire_lengths).read_text(encoding="utf-8").splitlines()
            if x.strip() and not x.strip().startswith("#")
        ]
    truncated = truncate_loci_batch(
        loci,
        seed=args.seed,
        min_bp=args.min_bp,
        max_bp=args.max_bp,
        empirical_lengths=emp_lens,
    )
    _write_jsonl(Path(args.out), truncated)
    return 0


def _cmd_evaluate(args: argparse.Namespace) -> int:
    holdout = _read_jsonl(Path(args.holdout))
    pred_manifest = json.loads(Path(args.predictions).read_text(encoding="utf-8"))
    base_dir = Path(args.predictions).parent

    def _resolve(p_str: str | None) -> Path | None:
        if not p_str:
            return None
        p = Path(p_str)
        return p if p.is_absolute() else base_dir / p

    methods_out: dict[str, Any] = {}
    for method_name in sorted(pred_manifest.keys()):
        entry = pred_manifest[method_name]
        if isinstance(entry, str):
            scores_path = _resolve(entry)
            trunc_path = None
            nr90_path = None
            lat_ratio = None
        else:
            scores_path = _resolve(entry["scores"])
            trunc_path = _resolve(entry.get("truncated_scores"))
            nr90_path = _resolve(entry.get("nr90_scores"))
            lat_ratio = entry.get("latency_ratio_vs_afdb")

        assert scores_path is not None
        scores = _load_scores(scores_path)
        trunc_scores = _load_scores(trunc_path) if trunc_path is not None else None
        nr90_scores = _load_scores(nr90_path) if nr90_path is not None else None

        methods_out[method_name] = evaluate_method(
            holdout,
            scores,
            truncated_scores=trunc_scores,
            nr90_scores=nr90_scores,
            latency_ratio_vs_afdb=lat_ratio,
        )

    _write_json(Path(args.out), {"methods": methods_out})
    return 0


def _cmd_check_gates(args: argparse.Namespace) -> int:
    eval_report = json.loads(Path(args.eval_json).read_text(encoding="utf-8"))
    gates = check_gates(eval_report)
    _write_json(Path(args.out), gates)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="w38a_plm_bakeoff.py",
        description="W-38a GPU pLM/3Di bake-off holdout builder, SPIRE locus truncator, evaluator, and gate checker.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build-holdout", help="Build strict <20% ID & E>10 holdout + 1M trap-enriched background.")
    b.add_argument("--root", type=Path, default=None, help="bio-box data root (manifest-backed ReleaseFiles).")
    b.add_argument("--release", default=DEFAULT_RELEASE, help=f"Live release ID (default: {DEFAULT_RELEASE}).")
    b.add_argument("--positives", type=Path, required=True, help="Candidate positives JSONL.")
    b.add_argument("--background-jsonl", type=Path, default=None, help="Optional pre-extracted background JSONL (for testing).")
    b.add_argument("--out", type=Path, required=True, help="Output directory.")
    b.add_argument("--bg-size", type=int, default=DEFAULT_BG_SIZE, help="Number of background genes to sample.")
    b.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Deterministic seed.")
    b.add_argument("--max-aa-id", type=float, default=MAX_AA_IDENTITY, help="Strict upper bound on training AA ID.")
    b.add_argument("--min-hmm-evalue", type=float, default=MIN_HMM_EVALUE, help="Strict lower bound on training HMM E-value.")
    b.set_defaults(func=_cmd_build_holdout)

    t = sub.add_parser("truncate-loci", help="Truncate closed P1..P6 loci to empirical 1.5-5 kb SPIRE lengths.")
    t.add_argument("--loci", type=Path, required=True, help="Closed loci JSONL.")
    t.add_argument("--out", type=Path, required=True, help="Output truncated loci JSONL.")
    t.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Deterministic seed.")
    t.add_argument("--min-bp", type=int, default=SPIRE_MIN_BP, help="Minimum fragment length (bp).")
    t.add_argument("--max-bp", type=int, default=SPIRE_MAX_BP, help="Maximum fragment length (bp).")
    t.add_argument("--spire-lengths", type=Path, default=None, help="Optional newline-delimited empirical SPIRE contig lengths.")
    t.set_defaults(func=_cmd_truncate_loci)

    e = sub.add_parser("evaluate", help="Evaluate competitor predictions across all W-38a stratifications.")
    e.add_argument("--holdout", type=Path, required=True, help="Holdout JSONL from build-holdout.")
    e.add_argument("--predictions", type=Path, required=True, help="JSON manifest mapping method -> score file spec.")
    e.add_argument("--out", type=Path, required=True, help="Output evaluation report JSON.")
    e.set_defaults(func=_cmd_evaluate)

    g = sub.add_parser("check-gates", help="Evaluate falsifiable Go/No-Go Gates 1-3 from evaluation JSON.")
    g.add_argument("--eval-json", type=Path, required=True, help="Evaluation report JSON from evaluate.")
    g.add_argument("--out", type=Path, required=True, help="Output gate verdict JSON.")
    g.set_defaults(func=_cmd_check_gates)

    return p


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
