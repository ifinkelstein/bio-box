"""Contig tapes, coordinates only (operation O-TAPE; PLAN.md §4.2; planning/tasks/TASKS.md §2; card B-07).

A tape version is a set of feature producers: `tape_version(producers)` is "t-" + the first 12
hex digits of sha256 over the sorted producer names joined by "\\n". A tape is one row per contig
of a stage (every default variant in `contig_index`, also a contig with no features) and tape
version, with parallel lists holding one entry per slot. The slots are the union of the `feature`
rows of the producers, ordered by (start, end, strand, kind); the list index is the slot's ordinal,
valid inside one tape version only, and `feature_uids` names each slot (AGENTS §3.3). With the
gene-call producer alone this is genecall's `slot_ordinal` order. Tapes are append-only per stage
and carry stable keys only (`feature_uid`, `protein_uid`); tokens are a release-time join
(`biobox.l3.tokens`).

S0 fills the columns L1 can supply: kinds, strands, coordinates, the junction to the previous
slot, censoring distances and the protein of each CDS slot (null for any other kind). Columns that
need later sidecars (segments, directons, regimes, taxonomy, structure) are added as nullable
columns when their producers exist (LY-09). `tax_order_bucket`, `genome_uid`, `tier` and
`assembly_quality` are legacy (card B-17 removes them): taken from the contig's first gene in
`protein_occurrence`, null for a contig without genes.

Junctions and distances on a linear contig: the first slot has no previous junction (null),
and `orfs_to_contig_starts/ends` count the slots between a slot and the contig ends,
saturating at 65,535 with `slot_flags` bit 0 set (C-07). On a circular contig the first
slot's previous slot is the last one across the origin, and the distance lists are null
(no edge).

Slots are in start order, so where features nest or share a start, reverse complementing a
contig reorders them and their slot offsets change; MR-03 holds exactly for slots whose start
and end orders agree (tests/r1/test_l3.py).

Writes, all through `biobox.io.write_table`:

    l3/tapes/tape_version=<tv>/stage=<S>/part-<k>.parquet   sorted by contig_uid, cut at
                                                            `max_file_bytes` (estimated, `io.write_parts`)
    l3/tape_version_def/tape_version=<tv>/part-0.parquet    the producers, once per tape version
    ledger/stage=<S>/step=tape-<tv>/part-0.parquet

Bounded memory (AGENTS §4.1): DuckDB sorts the features (spilling into a private temporary
directory) and an Arrow reader hands them over in `contig_uid` order; one contig's features are
held at a time, and a row group holds whole contigs and at most `max_slots` slots (AGENTS §4.6).
The ledger (one row per contig) streams to its writer from the tape parts just written (card B-19).
"""

from __future__ import annotations

import hashlib
import shutil
import tempfile
from collections.abc import Iterable, Iterator
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from biobox import io, layout, ledger, provenance
from biobox.ids import KIND
from biobox.l1.genecall import DUCK_MEMORY_LIMIT, GENE_CALL_VERSION
from biobox.l3.geometry import pair_geometry
from biobox.registry.schemas import table

KIND_CDS = KIND["cds"]
U16_MAX = 65_535
FLAG_DIST_CLIPPED = 1
MAX_SLOTS = 131_072  # slots per row group (AGENTS §4.6)
EST_ROW_BYTES, EST_SLOT_BYTES = 200, 30  # zstd estimate per contig row and per slot; cuts parts only
BATCH_ROWS = 4_096  # contig rows handed to the writer at a time


def tape_version(producers: Iterable[str]) -> str:
    """"t-" + the first 12 hex digits of sha256 over the sorted producer names joined by "\\n" (UTF-8)."""
    names = list(producers)
    if not names:
        raise ValueError("a tape version needs at least one producer")
    if len(set(names)) != len(names):
        raise ValueError(f"a producer is listed twice: {sorted(names)}")
    return "t-" + hashlib.sha256("\n".join(sorted(names)).encode("utf-8")).hexdigest()[:12]


def build_tape(slots: list[tuple[int, int, int]], contig_len: int, circular: bool,
               kinds: list[int] | None = None) -> dict:
    """List columns for one contig from its slots (start, end, strand) in ordinal order; `kinds` are
    the slots' feature kinds (default: every slot a CDS)."""
    n = len(slots)
    if kinds is None:
        kinds = [KIND_CDS] * n
    elif len(kinds) != n:
        raise ValueError(f"{len(kinds)} kinds for {n} slots")
    delta: list[int | None] = []
    junc: list[int | None] = []
    for i in range(n):
        if i == 0 and (not circular or n == 1):
            delta.append(None)
            junc.append(None)
            continue
        t, d, _ = pair_geometry(slots, i, (i - 1) % n, -1, contig_len, circular)
        junc.append(t)
        delta.append(d)
    flags = [0] * n
    to_start = to_end = None
    if not circular:
        to_start, to_end = [], []
        for i in range(n):
            a, b = i, n - 1 - i
            if a > U16_MAX or b > U16_MAX:
                flags[i] |= FLAG_DIST_CLIPPED
            to_start.append(min(a, U16_MAX))
            to_end.append(min(b, U16_MAX))
    return {
        "feature_kinds": list(kinds),
        "strands": [s[2] for s in slots],
        "start_bps": [s[0] for s in slots],
        "end_bps": [s[1] for s in slots],
        "delta_bp_prev_slots": delta,
        "junction_topology_prevs": junc,
        "orfs_to_contig_starts": to_start,
        "orfs_to_contig_ends": to_end,
        "slot_flags": flags,
    }


def _sql_path(p: Path) -> str:
    return "'" + str(p).replace("'", "''") + "'"


def _source(files: list[Path], cols: str, empty: str) -> str:
    """A subquery over `files`, or a typed empty one when there are none (read_parquet needs a file)."""
    if not files:
        return f"(SELECT {empty} WHERE false)"
    return f"(SELECT {cols} FROM read_parquet([{', '.join(_sql_path(f) for f in files)}]))"


def def_path(root: Path, tv: str) -> Path:
    return Path(root) / table("tape_version_def").location / f"tape_version={tv}" / "part-0.parquet"


def _stored_producers(root: Path, tv: str, producers: list[str]) -> bool:
    """True if the tape version's definition exists; raises if it lists other producers."""
    path = def_path(root, tv)
    if not path.exists():
        return False
    have = sorted(pq.read_table(path, columns=["producer"], partitioning=None).column("producer").to_pylist())
    if have != sorted(producers):
        raise ValueError(f"{path}: tape_version_def of {tv} lists {have}, not {sorted(producers)}")
    return True


def ensure_def(root: Path, producers: list[str]) -> str:
    """Make sure the tape version of `producers` has its `tape_version_def` file: write it if absent,
    verify it lists the same producers if present (ValueError otherwise). Returns the tape version.
    Shared by `tapes.run` and `biobox.l1.igr.run` (card B-20)."""
    tv = tape_version(producers)
    if not _stored_producers(root, tv, producers):
        io.write_table([pa.table({"producer": sorted(producers)})], table("tape_version_def"), def_path(root, tv))
    return tv


def run(root: Path, stage: str, circular: set[bytes] | frozenset | None = None,
        gene_call_version: str = GENE_CALL_VERSION, producers: list[str] | None = None,
        max_file_bytes: int = layout.MAX_FILE_BYTES, max_slots: int = MAX_SLOTS, tmp=None) -> dict:
    """Write the stage's tapes for the tape version of `producers` (default `[gene_call_version]`).

    `circular`: the circular contig_uids; None reads the stage's `genome_input` table. Proteins come
    from `feature_protein` of `gene_call_version` (CDS slots only), the legacy columns from its
    `protein_occurrence`. `tmp` holds DuckDB's private spill directory (system temp by default).
    """
    root = Path(root)
    producers = [gene_call_version] if producers is None else list(producers)
    tv = tape_version(producers)
    if max_slots < 1:
        raise ValueError("max_slots must be positive")
    if tmp is not None and Path(tmp).resolve().is_relative_to(root.resolve()):
        raise ValueError(f"{tmp}: the scratch directory may not lie under the data root")
    _stored_producers(root, tv, producers)  # fail before any work if the stored definition disagrees
    feats = []
    for p in sorted(producers):
        found = io.part_files(root / f"l1/feature/producer={p}/stage={stage}")
        if not found:
            raise FileNotFoundError(f"no feature files for producer {p} at stage {stage}")
        feats += found
    gcv = f"gene_call_version={gene_call_version}/stage={stage}"
    fp = io.part_files(root / f"l1/feature_protein/{gcv}")
    occ = io.part_files(root / f"l1/protein_occurrence/{gcv}")
    idx = io.part_files(root / f"l0seq/contig_index/stage={stage}")
    if not idx:
        raise FileNotFoundError(f"no contig_index files at stage {stage}")
    gi = []
    if circular is None:
        gi = io.part_files(root / table("genome_input").location / f"stage={stage}")
        if not gi:
            raise FileNotFoundError(f"no genome_input for stage {stage}: pass `circular` or write it first")
    sql = (
        f"WITH c AS (SELECT contig_uid, contig_len FROM read_parquet([{', '.join(_sql_path(f) for f in idx)}]) "
        "WHERE is_default_variant), "
        f"f AS (SELECT DISTINCT * FROM {_source(feats, 'feature_uid, contig_uid, start, \"end\", strand, kind', '')}), "
        f"fp AS {_source(fp, 'feature_uid, protein_uid', 'NULL::BLOB AS feature_uid, NULL::BLOB AS protein_uid')}, "
        "o AS " + _source(occ, "feature_uid, slot_ordinal, genome_uid, tier, assembly_quality",
                          "NULL::BLOB AS feature_uid, NULL::UINTEGER AS slot_ordinal, NULL::BLOB AS genome_uid, "
                          "NULL::UTINYINT AS tier, NULL::VARCHAR AS assembly_quality") + ", "
        "g AS " + _source(gi, "contig_uid, is_circular", "NULL::BLOB AS contig_uid, NULL::BOOLEAN AS is_circular") + " "
        "SELECT coalesce(c.contig_uid, f.contig_uid) AS contig_uid, c.contig_len, g.is_circular, f.feature_uid, "
        f'f.start, f."end", f.strand, f.kind, CASE WHEN f.kind = {KIND_CDS} THEN fp.protein_uid END AS protein_uid, '
        "o.slot_ordinal, o.genome_uid, o.tier, o.assembly_quality "
        "FROM c FULL OUTER JOIN f ON f.contig_uid = c.contig_uid "
        "LEFT JOIN fp ON fp.feature_uid = f.feature_uid "
        "LEFT JOIN o ON o.feature_uid = f.feature_uid "
        "LEFT JOIN g ON g.contig_uid = coalesce(c.contig_uid, f.contig_uid) "
        'ORDER BY 1, f.start, f."end", f.strand, f.kind'
    )
    spec = table("contig_tapes")
    tools = provenance.tool_versions("pyarrow")
    counts = {"contigs": 0, "slots": 0}
    scratch = Path(tempfile.mkdtemp(prefix="biobox-tapes-", dir=tmp))
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb"), "preserve_insertion_order": True})

    def contig_row(cu: bytes, length: int | None, circ: bool, fs: list[tuple]) -> dict:
        if length is None:
            raise KeyError(f"features on contig {cu.hex()} that contig_index does not hold")
        genes = [x for x in fs if x[6] is not None]
        first = min(genes, key=lambda x: x[6]) if genes else (None,) * 10
        n = len(fs)
        counts["contigs"] += 1
        counts["slots"] += n
        return {
            "contig_uid": cu,
            "tax_order_bucket": None,
            "genome_uid": first[7],
            "contig_len": length,
            "is_circular": circ,
            "n_slots": n,
            "tier": first[8],
            "assembly_quality": first[9],
            **build_tape([(x[1], x[2], x[3]) for x in fs], length, circ, [x[4] for x in fs]),
            "protein_uids": [x[5] for x in fs],
            "feature_uids": [x[0] for x in fs],
        }

    def rows() -> Iterator[dict]:
        cur, length, circ, fs = None, None, False, []
        reader = con.sql(sql).to_arrow_reader(io.ROW_GROUP_ROWS)
        try:
            for b in reader:
                cols = [b.column(i).to_pylist() for i in range(b.num_columns)]
                for cu, n, gcirc, fu, s, e, strand, kind, pu, slot, guid, tier, aq in zip(*cols):
                    if cu != cur:
                        if cur is not None:
                            yield contig_row(cur, length, circ, fs)
                        cur, length, fs = cu, n, []
                        circ = bool(gcirc) if circular is None else cu in circular
                    if fu is not None:
                        fs.append((fu, s, e, strand, kind, pu, slot, guid, tier, aq))
            if cur is not None:
                yield contig_row(cur, length, circ, fs)
        finally:
            reader.close()

    def batches() -> Iterator[pa.RecordBatch]:
        """The tape rows, at most `BATCH_ROWS` rows or `max_slots` slots held at a time."""
        buf, buf_slots = [], 0
        for r in rows():
            buf.append(r)
            buf_slots += r["n_slots"]
            if len(buf) >= BATCH_ROWS or buf_slots >= max_slots:
                yield pa.RecordBatch.from_pylist(buf, schema=spec.schema)
                buf, buf_slots = [], 0
        if buf:
            yield pa.RecordBatch.from_pylist(buf, schema=spec.schema)

    def est(b: pa.RecordBatch) -> np.ndarray:
        return EST_ROW_BYTES + EST_SLOT_BYTES * b.column("n_slots").to_numpy(zero_copy_only=False).astype(np.float64)

    base = root / f"l3/tapes/tape_version={tv}/stage={stage}"
    try:
        paths = io.write_parts(batches(), spec, base, est, max_file_bytes, max_slots=max_slots,
                               slots_column="n_slots", statistics=["n_slots", "contig_len"])
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)
    ensure_def(root, producers)

    def ledger_rows() -> Iterator[dict]:
        """Per contig: slots in, all kept; read back from the parts just written, in contig_uid order."""
        for f in paths:
            for b in pq.ParquetFile(f).iter_batches(batch_size=io.ROW_GROUP_ROWS, columns=["contig_uid", "n_slots"]):
                for cu, n in zip(b.column("contig_uid").to_pylist(), b.column("n_slots").to_pylist()):
                    yield ledger.row(cu, "contig", n, n, tools=tools)

    ledger.write(root, stage, f"tape-{tv}", ledger_rows())  # per tape version: a second version never rewrites it
    return {**counts, "files": len(paths), "tape_version": tv}
