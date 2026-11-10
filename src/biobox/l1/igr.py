"""Intergenic regions as coordinates (data_layer_and_tests.md §2A.1; R-05).

An IGR is a maximal run of bases covered by no feature of the tape's producers (a feature of
any kind covers its bases), so overlapping and nested features produce none between them.
Every run of >= 1 nt is stored; the >= 30 nt filter (Weinberg convention, SEARCH_PLAN SD-07) is applied when the search
pack is built. Coordinates are 0-based half-open on the + strand (`strand` = 1);
orientation is derived from the flanking features (`left_feature_uid`, `right_feature_uid`, of
any kind). On linear contigs the end runs have a null flank; on circular contigs the run across the origin has `end > contig_len`
(C-05), and a contig with no genes is one IGR with both flanks null.

Writes `<root>/l0seq/igr/igr_version=<v>/stage=<S>/part-<k>.parquet` (cut at `max_file_bytes` of
estimated size, `io.write_parts`, card B-19), where `<v>` is `<tape_version>-igr1` (card B-20) and the
tape version's `tape_version_def` file is made sure to exist.
"""

from __future__ import annotations

import shutil
import tempfile
from collections.abc import Iterable
from pathlib import Path

import duckdb
import pyarrow as pa

from biobox import io, layout
from biobox.l0seq.inputs import iter_genome_input
from biobox.l1.genecall import DUCK_MEMORY_LIMIT, GENE_CALL_VERSION
from biobox.l3.tapes import ensure_def, tape_version  # tapes.py does not import this module: no cycle
from biobox.registry.schemas import table

IGR_RULE = "igr1"
# estimate: zstd bytes per IGR row; cuts parts only. Measured 12.7 B/row on the golden set and 9.7 on the
# S1 pilot (planning/measurements/2026-11-09_design_stress_test.md §3, §6a) before the two flanking
# feature_uid columns (16 random bytes each) were added: 13 + 32.
EST_ROW_BYTES = 45


def igr_version(producers: str | Iterable[str] = GENE_CALL_VERSION) -> str:
    """`<tape_version>-igr1`: the IGR version names the producer set through its tape version, so IGRs of
    a second producer set land beside the first. `producers` is a list of producer names, or (the old
    single-argument meaning) a plain string, which is the single producer `[gene_call_version]`."""
    names = [producers] if isinstance(producers, str) else list(producers)
    return f"{tape_version(names)}-{IGR_RULE}"


def igrs_of_contig(slots: list[tuple[int, int, int]], length: int, circular: bool) -> list[tuple[int, int, int | None, int | None]]:
    """(start, end, left, right) for one contig; slots are (name, start, end), and `left`/`right` are the
    names of the flanking slots (None at a linear contig end). The caller picks the names: `_contig_igrs`
    uses each feature's position in its sorted list and maps it back to the feature_uid."""
    if not slots:
        return [(0, length, None, None)] if length else []
    genes = sorted(slots, key=lambda s: (s[1], s[2], s[0]))
    out = []
    wraps = [(e - length, o) for o, s, e in genes if e > length] if circular else []
    if wraps:
        # an origin-wrapping gene also covers [0, end - length); scan from its tail, and the
        # run before the first gene closes the circle, so no closing run is added
        reach, reach_slot = max(wraps)
        rest = genes
    else:
        reach, reach_slot = genes[0][2], genes[0][0]
        rest = genes[1:]
        if not circular and genes[0][1] > 0:
            out.append((0, genes[0][1], None, genes[0][0]))
    for ordinal, s, e in rest:
        if s > reach:
            out.append((reach, s, reach_slot, ordinal))
        if e > reach:
            reach, reach_slot = e, ordinal
    if circular and not wraps:
        first_ordinal, first_start = genes[0][0], genes[0][1]
        if reach < first_start + length:
            start, end = reach, first_start + length if first_start else length
            if start >= length:  # the last gene ends exactly at the origin
                start, end = start - length, end - length
            out.append((start, end, reach_slot, first_ordinal))
    elif not circular and reach < length:
        out.append((reach, length, reach_slot, None))
    return sorted(out)


def _sql_path(p: Path) -> str:
    return "'" + str(p).replace("'", "''") + "'"


def _contig_igrs(cu: bytes, length: int, feats: list[tuple], circular: bool) -> list[dict]:
    """IGR rows of one contig; `feats` are (feature_uid, start, end) in sort order. The position in
    `feats` names a flank for `igrs_of_contig` only; the row stores its feature_uid."""
    runs = igrs_of_contig([(i, f[1], f[2]) for i, f in enumerate(feats)], length, circular)
    return [{"contig_uid": cu, "start": s, "end": e, "strand": 1,
             "left_feature_uid": None if left is None else feats[left][0],
             "right_feature_uid": None if right is None else feats[right][0]} for s, e, left, right in runs]


def run(root: Path, stage: str, circular: set[bytes] | frozenset | None = None,
        gene_call_version: str = GENE_CALL_VERSION, producers: list[str] | None = None, tmp=None,
        max_file_bytes: int = layout.MAX_FILE_BYTES) -> dict:
    """Build IGRs for every contig of a stage from the `feature` rows of `producers` and contig_index.

    `circular`: the circular contig_uids; None reads the stage's `genome_input` table. `producers`:
    default `[gene_call_version]`; `gene_call_version` names the default producer only.
    Contigs are streamed in `contig_uid` order; one contig's features are held at a time.
    """
    root = Path(root)
    producers = [gene_call_version] if producers is None else list(producers)
    v = igr_version(producers)
    files = []
    for p in producers:
        found = io.part_files(root / f"l1/feature/producer={p}/stage={stage}")
        if not found:
            raise FileNotFoundError(f"no feature files for producer {p} at stage {stage}")
        files += found
    ensure_def(root, producers)  # write the tape version's definition, or verify the stored one
    idx = io.part_files(root / f"l0seq/contig_index/stage={stage}")
    if not idx:
        raise FileNotFoundError(f"no contig_index files at stage {stage}")
    scratch = Path(tempfile.mkdtemp(prefix="biobox-igr-", dir=tmp))
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb")})
    try:
        if circular is None:
            circular = frozenset(cu for g in iter_genome_input(root, stage, tmp=scratch) for cu in g.circular)
        flist = "[" + ", ".join(_sql_path(f) for f in files) + "]"
        ilist = "[" + ", ".join(_sql_path(f) for f in idx) + "]"
        sql = (
            "SELECT coalesce(c.contig_uid, f.contig_uid) AS contig_uid, c.contig_len, f.feature_uid, f.start, "
            'f."end", f.strand, f.kind '
            f"FROM (SELECT contig_uid, contig_len FROM read_parquet({ilist}) WHERE is_default_variant) c "
            f"FULL OUTER JOIN read_parquet({flist}) f ON f.contig_uid = c.contig_uid "
            'ORDER BY 1, f.start, f."end", f.strand, f.kind'
        )
        spec = table("igr")
        stats = {"igrs": 0, "igrs_ge_30": 0, "bp": 0}

        def batches():
            out: list[dict] = []
            cur, length, feats = None, 0, []

            def close():
                if cur is None:
                    return
                if length is None:
                    raise KeyError(f"features on contig {cur.hex()} that contig_index does not hold")
                for r in _contig_igrs(cur, length, feats, cur in circular):
                    out.append(r)
                    stats["igrs"] += 1
                    stats["igrs_ge_30"] += r["end"] - r["start"] >= 30
                    stats["bp"] += r["end"] - r["start"]

            reader = con.sql(sql).to_arrow_reader(io.ROW_GROUP_ROWS)
            for b in reader:
                for cu, n, fu, s, e in zip(*(b.column(i).to_pylist() for i in range(5))):
                    if cu != cur:
                        close()
                        cur, length, feats = cu, n, []
                    if fu is not None:
                        feats.append((fu, s, e))
                if len(out) >= io.ROW_GROUP_ROWS:
                    yield pa.RecordBatch.from_pylist(out, schema=spec.schema)
                    out = []
            close()
            if out:
                yield pa.RecordBatch.from_pylist(out, schema=spec.schema)

        io.write_parts(batches(), spec, root / f"l0seq/igr/igr_version={v}/stage={stage}", EST_ROW_BYTES, max_file_bytes)
        return stats
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)
