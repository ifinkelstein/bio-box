"""Intergenic regions as coordinates (data_layer_and_tests.md §2A.1; R-05).

An IGR is a maximal run of bases covered by no feature of the tape's producers (a gene of
any kind of feature covers its bases), so overlapping and nested features produce none
between them. Every run of >= 1 nt is stored;
the >= 30 nt filter (Weinberg convention, SEARCH_PLAN SD-07) is applied when the search
pack is built. Coordinates are 0-based half-open on the + strand (`strand` = 1);
orientation is derived from the flanking features (`left_feature_uid`, `right_feature_uid`;
the legacy `left_slot`, `right_slot` stay until card B-17). On linear contigs the end runs have a
null flank; on circular contigs the run across the origin has `end > contig_len`
(C-05), and a contig with no genes is one IGR with both flanks null.

Writes `<root>/l0seq/igr/igr_version=<v>/stage=<S>/part-0.parquet`.
"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import duckdb
import pyarrow as pa

from biobox import io
from biobox.l0seq.inputs import iter_genome_input
from biobox.l1.genecall import DUCK_MEMORY_LIMIT, GENE_CALL_VERSION
from biobox.registry.schemas import table

IGR_RULE = "igr1"


def igr_version(gene_call_version: str = GENE_CALL_VERSION) -> str:
    return f"{gene_call_version}-{IGR_RULE}"


def igrs_of_contig(slots: list[tuple[int, int, int]], length: int, circular: bool) -> list[tuple[int, int, int | None, int | None]]:
    """(start, end, left_slot, right_slot) for one contig; slots are (slot_ordinal, start, end)."""
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
    """IGR rows of one contig; `feats` are (feature_uid, start, end, slot_ordinal) in sort order. The
    position in `feats` stands in for the ordinal that `igrs_of_contig` uses to name the flanks."""
    runs = igrs_of_contig([(i, f[1], f[2]) for i, f in enumerate(feats)], length, circular)
    return [{"contig_uid": cu, "start": s, "end": e, "strand": 1,
             "left_slot": None if left is None else feats[left][3],
             "right_slot": None if right is None else feats[right][3],
             "left_feature_uid": None if left is None else feats[left][0],
             "right_feature_uid": None if right is None else feats[right][0]} for s, e, left, right in runs]


def run(root: Path, stage: str, circular: set[bytes] | frozenset | None = None,
        gene_call_version: str = GENE_CALL_VERSION, producers: list[str] | None = None, tmp=None) -> dict:
    """Build IGRs for every contig of a stage from the `feature` rows of `producers` and contig_index.

    `circular`: the circular contig_uids; None reads the stage's `genome_input` table. `producers`:
    default `[gene_call_version]`. `left_slot`/`right_slot` (legacy, B-17 removes them) are the
    `slot_ordinal` of a gene flank taken from protein_occurrence, null for any other flank.
    Contigs are streamed in `contig_uid` order; one contig's features are held at a time.
    """
    root = Path(root)
    producers = [gene_call_version] if producers is None else list(producers)
    files = []
    for p in producers:
        found = sorted((root / f"l1/feature/producer={p}/stage={stage}").glob("part-*.parquet"))
        if not found:
            raise FileNotFoundError(f"no feature files for producer {p} at stage {stage}")
        files += found
    occ = sorted((root / f"l1/protein_occurrence/gene_call_version={gene_call_version}/stage={stage}").glob("part-*.parquet"))
    idx = root / f"l0seq/contig_index/stage={stage}/part-0.parquet"
    scratch = Path(tempfile.mkdtemp(prefix="biobox-igr-", dir=tmp))
    con = duckdb.connect(config={"threads": 1, "memory_limit": DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb")})
    try:
        if circular is None:
            circular = frozenset(cu for g in iter_genome_input(root, stage, tmp=scratch) for cu in g.circular)
        flist = "[" + ", ".join(_sql_path(f) for f in files) + "]"
        olist = "[" + ", ".join(_sql_path(f) for f in occ) + "]"
        sql = (
            "SELECT coalesce(c.contig_uid, f.contig_uid) AS contig_uid, c.contig_len, f.feature_uid, f.start, "
            'f."end", f.strand, f.kind, o.slot_ordinal '
            f"FROM (SELECT contig_uid, contig_len FROM read_parquet({_sql_path(idx)}) WHERE is_default_variant) c "
            f"FULL OUTER JOIN read_parquet({flist}) f ON f.contig_uid = c.contig_uid "
            f"LEFT JOIN read_parquet({olist}) o ON o.feature_uid = f.feature_uid "
            'ORDER BY 1, f.start, f."end", f.strand, f.kind, o.slot_ordinal'
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
                for cu, n, fu, s, e, _, _, slot in zip(*(b.column(i).to_pylist() for i in range(8))):
                    if cu != cur:
                        close()
                        cur, length, feats = cu, n, []
                    if fu is not None:
                        feats.append((fu, s, e, slot))
                if len(out) >= io.ROW_GROUP_ROWS:
                    yield pa.RecordBatch.from_pylist(out, schema=spec.schema)
                    out = []
            close()
            if out:
                yield pa.RecordBatch.from_pylist(out, schema=spec.schema)

        v = igr_version(gene_call_version)
        io.write_table(batches(), spec, root / f"l0seq/igr/igr_version={v}/stage={stage}/part-0.parquet")
        return stats
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)
