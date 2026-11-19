"""Release-level merge of annotation rounds into O(1)-pruning lookup copies (card C-38g, DEC-54).

A release lists the complete rounds it takes in `release_round` (card C-38e). Per-round files live under:

    l1/round_domains/round=<rid>/library=<L>/part-<k>.parquet
    l1/round_coverage/round=<rid>/library=<L>/part-<k>.parquet

The release step `rounds` (after `cluster`, before `tokens`) performs a bounded-memory streaming k-way
merge by `protein_uid` over the rounds listed in `release_round` for each library, deduplicating proteins
covered by overlapping rounds deterministically (lowest round index in `release_round` order; by MR-12
their domain rows are identical), and writes the merged sorted copies under the alternative partition keys:

    l1/round_domains/release=<R>/library=<L>/part-<k>.parquet
    l1/round_coverage/release=<R>/library=<L>/part-<k>.parquet

When `release_round` lists 0 rounds, `merge_release` is a clean no-op and writes nothing.

CLI entrypoint (invoked by `stage.build_stage` via `run-step` with `$BIOBOX_STAGING`):

    python -m biobox.l1.rounds merge-release --root R --release REL [--out OUT] [--tmp TMP]
"""

from __future__ import annotations

import argparse
import heapq
import json
import os
import sys
from collections.abc import Iterator
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from biobox import catalog, layout, ledger
from biobox import io as bio
from biobox.l1.domains import EST_COV_ROW_BYTES, EST_ROW_BYTES
from biobox.registry.schemas import table

_STREAM_BATCH = 4_096  # rows per Arrow batch read/flushed during k-way merge (bounds memory per round)
_FLUSH_PIECES = 128


def release_round_domains_dir(release: str, library: str) -> str:
    return f"{table('round_domains').location}/release={release}/library={library}"


def release_round_coverage_dir(release: str, library: str) -> str:
    return f"{table('round_coverage').location}/release={release}/library={library}"


def _round_libraries(root: Path, rid: str) -> dict[str, tuple[list[Path], list[Path]]]:
    """{library: (round_domains_files, round_coverage_files)} for one round `rid` under `root`."""
    dom_base = root / table("round_domains").location / f"round={rid}"
    cov_base = root / table("round_coverage").location / f"round={rid}"
    dom_libs: dict[str, list[Path]] = {}
    if dom_base.is_dir():
        for d in sorted(dom_base.iterdir()):
            if d.is_dir() and d.name.startswith("library="):
                files = bio.part_files(d)
                if files:
                    dom_libs[d.name.removeprefix("library=")] = files
    cov_libs: dict[str, list[Path]] = {}
    if cov_base.is_dir():
        for d in sorted(cov_base.iterdir()):
            if d.is_dir() and d.name.startswith("library="):
                files = bio.part_files(d)
                if files:
                    cov_libs[d.name.removeprefix("library=")] = files
    if not dom_libs or set(dom_libs) != set(cov_libs):
        raise FileNotFoundError(
            f"round {rid} under {root} is incomplete: round_domains libraries {sorted(dom_libs)} "
            f"vs round_coverage libraries {sorted(cov_libs)}"
        )
    return {lib: (dom_libs[lib], cov_libs[lib]) for lib in sorted(dom_libs)}


def _iter_domains_by_protein(files: list[Path], schema: pa.Schema) -> Iterator[tuple[bytes, list[pa.RecordBatch]]]:
    """Stream `(protein_uid, record_batch_slices)` in strictly increasing `protein_uid` order, one small batch at a time."""
    held_uid: bytes | None = None
    held_parts: list[pa.RecordBatch] = []
    last_emitted: bytes | None = None

    for path in files:
        pf = pq.ParquetFile(path)
        for batch in pf.iter_batches(batch_size=_STREAM_BATCH, columns=schema.names, use_threads=False):
            if batch.num_rows == 0:
                continue
            uids = batch.column("protein_uid").to_pylist()
            n = len(uids)
            start = 0
            for i in range(1, n):
                if uids[i] != uids[start]:
                    u = uids[start]
                    sl = batch.slice(start, i - start)
                    if held_uid == u:
                        held_parts.append(sl)
                    else:
                        if held_uid is not None:
                            if last_emitted is not None and held_uid <= last_emitted:
                                raise ValueError(f"{path}: round_domains not sorted by protein_uid")
                            last_emitted = held_uid
                            yield held_uid, held_parts
                        held_uid = u
                        held_parts = [sl]
                    start = i
            u = uids[start]
            sl = batch.slice(start, n - start)
            if held_uid == u:
                held_parts.append(sl)
            else:
                if held_uid is not None:
                    if last_emitted is not None and held_uid <= last_emitted:
                        raise ValueError(f"{path}: round_domains not sorted by protein_uid")
                    last_emitted = held_uid
                    yield held_uid, held_parts
                held_uid = u
                held_parts = [sl]

    if held_uid is not None:
        if last_emitted is not None and held_uid <= last_emitted:
            raise ValueError("round_domains not sorted by protein_uid")
        yield held_uid, held_parts


def _iter_coverage(files: list[Path]) -> Iterator[tuple[bytes, int]]:
    """Stream `(protein_uid, n_rows)` in strictly increasing `protein_uid` order, one small batch at a time."""
    last: bytes | None = None
    for path in files:
        pf = pq.ParquetFile(path)
        for batch in pf.iter_batches(batch_size=_STREAM_BATCH, columns=["protein_uid", "n_rows"], use_threads=False):
            if batch.num_rows == 0:
                continue
            for u, nr in zip(batch.column("protein_uid").to_pylist(), batch.column("n_rows").to_pylist()):
                if last is not None and u <= last:
                    raise ValueError(f"{path}: round_coverage not strictly increasing by protein_uid")
                last = u
                yield u, int(nr)


def _merge_domains_stream(
    round_dom_files: list[list[Path]],
    spec,
    row_group_rows: int,
    stats: dict[str, int],
) -> Iterator[pa.Table]:
    """Streaming k-way merge of per-round `round_domains` by `protein_uid`; lowest round index wins on overlap."""
    iters = [_iter_domains_by_protein(files, spec.schema) for files in round_dom_files]
    heap: list[tuple[bytes, int, list[pa.RecordBatch]]] = []
    for idx, it in enumerate(iters):
        first = next(it, None)
        if first is not None:
            heapq.heappush(heap, (first[0], idx, first[1]))

    pend: list[pa.RecordBatch] = []
    held = 0
    last_uid: bytes | None = None
    flush_rows = min(row_group_rows, _STREAM_BATCH)

    while heap:
        u, idx, batches = heapq.heappop(heap)
        nxt = next(iters[idx], None)
        if nxt is not None:
            heapq.heappush(heap, (nxt[0], idx, nxt[1]))
        if u == last_uid:
            continue
        last_uid = u
        pend.extend(batches)
        n_u = sum(b.num_rows for b in batches)
        held += n_u
        stats["domains"] += n_u
        stats["proteins_with_hit"] += 1
        if held >= flush_rows or len(pend) >= _FLUSH_PIECES:
            yield pa.Table.from_batches(pend, schema=spec.schema).combine_chunks()
            pend, held = [], 0

    if pend:
        yield pa.Table.from_batches(pend, schema=spec.schema).combine_chunks()


def _merge_coverage_stream(
    round_cov_files: list[list[Path]],
    spec,
    row_group_rows: int,
    stats: dict[str, int],
) -> Iterator[pa.RecordBatch]:
    """Streaming k-way merge of per-round `round_coverage` by `protein_uid`; lowest round index wins on overlap."""
    iters = [_iter_coverage(files) for files in round_cov_files]
    heap: list[tuple[bytes, int, int]] = []
    for idx, it in enumerate(iters):
        first = next(it, None)
        if first is not None:
            heapq.heappush(heap, (first[0], idx, first[1]))

    u_buf: list[bytes] = []
    n_buf: list[int] = []
    last_uid: bytes | None = None
    flush_rows = min(row_group_rows, _STREAM_BATCH)

    while heap:
        u, idx, nr = heapq.heappop(heap)
        nxt = next(iters[idx], None)
        if nxt is not None:
            heapq.heappush(heap, (nxt[0], idx, nxt[1]))
        if u == last_uid:
            continue
        last_uid = u
        u_buf.append(u)
        n_buf.append(nr)
        stats["proteins"] += 1
        if len(u_buf) >= flush_rows:
            yield pa.RecordBatch.from_arrays(
                [pa.array(u_buf, pa.binary(16)), pa.array(n_buf, pa.uint32())],
                schema=spec.schema,
            )
            u_buf, n_buf = [], []

    if u_buf:
        yield pa.RecordBatch.from_arrays(
            [pa.array(u_buf, pa.binary(16)), pa.array(n_buf, pa.uint32())],
            schema=spec.schema,
        )


def merge_release(
    root: Path,
    release: str,
    *,
    out: Path | None = None,
    max_file_bytes: int = layout.MAX_FILE_BYTES,
    row_group_rows: int = bio.ROW_GROUP_ROWS,
    tmp: Path | None = None,
) -> dict:
    """Merge the rounds listed in `release`'s `release_round` into sorted `release=` lookup files.

    Reads per-round `round_domains` and `round_coverage` under `root`, and writes the merged copies under
    `out` (default `root`, or `$BIOBOX_STAGING` under `run-step`). Clean no-op when `release` lists 0 rounds.
    Bounded memory: holds O(num_rounds) row groups at a time (AGENTS §4.1).
    """
    root = Path(root)
    out_root = root if out is None else Path(out)
    ledger.check_scratch(tmp, root, out_root)
    catalog.release_stages(root, release)
    rounds = catalog.release_rounds(root, release)
    if not rounds:
        return {"release": release, "rounds": 0, "libraries": {}, "skipped": True}

    by_lib: dict[str, tuple[list[list[Path]], list[list[Path]]]] = {}
    for rid in rounds:
        libs = _round_libraries(root, rid)
        for lib, (dom_files, cov_files) in libs.items():
            dom_list, cov_list = by_lib.setdefault(lib, ([], []))
            dom_list.append(dom_files)
            cov_list.append(cov_files)

    targets = [
        d for lib in sorted(by_lib)
        for d in (release_round_domains_dir(release, lib), release_round_coverage_dir(release, lib))
    ]
    ledger.refuse_written(root, out_root, targets, f"merged rounds of release {release}")

    dom_spec = table("round_domains")
    cov_spec = table("round_coverage")
    lib_stats: dict[str, dict[str, int]] = {}
    for lib in sorted(by_lib):
        dom_lists, cov_lists = by_lib[lib]
        st = {"rounds": len(dom_lists), "proteins": 0, "proteins_with_hit": 0, "domains": 0}
        dom_parts = bio.write_parts(
            _merge_domains_stream(dom_lists, dom_spec, row_group_rows, st),
            dom_spec,
            out_root / release_round_domains_dir(release, lib),
            EST_ROW_BYTES,
            max_file_bytes,
            row_group_rows=row_group_rows,
        )
        pa.default_memory_pool().release_unused()
        cov_parts = bio.write_parts(
            _merge_coverage_stream(cov_lists, cov_spec, row_group_rows, st),
            cov_spec,
            out_root / release_round_coverage_dir(release, lib),
            EST_COV_ROW_BYTES,
            max_file_bytes,
            row_group_rows=row_group_rows,
        )
        pa.default_memory_pool().release_unused()
        st["domain_files"] = len(dom_parts)
        st["coverage_files"] = len(cov_parts)
        lib_stats[lib] = st

    return {"release": release, "rounds": len(rounds), "libraries": lib_stats, "skipped": False}


def main(argv: list[str] | None = None) -> int:
    """`python -m biobox.l1.rounds merge-release --root R --release REL [--out OUT] [--tmp TMP]`."""
    ap = argparse.ArgumentParser(prog="python -m biobox.l1.rounds")
    ap.add_argument("action", choices=["merge-release"])
    ap.add_argument("--root", required=True)
    ap.add_argument("--release", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--tmp", default=None)
    ap.add_argument("--max-file-bytes", type=int, default=layout.MAX_FILE_BYTES)
    a = ap.parse_args(argv)
    out = Path(a.out) if a.out else (Path(os.environ["BIOBOX_STAGING"]) if os.environ.get("BIOBOX_STAGING") else None)
    try:
        res = merge_release(
            Path(a.root),
            a.release,
            out=out,
            max_file_bytes=a.max_file_bytes,
            tmp=Path(a.tmp) if a.tmp else None,
        )
    except (ValueError, FileNotFoundError, FileExistsError) as e:
        print(f"rounds {a.action}: {e}", file=sys.stderr)
        return 2
    print(json.dumps({"action": a.action, "result": res}, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
