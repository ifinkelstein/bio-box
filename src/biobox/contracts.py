"""Layout and contract checks over a release (data_layer_and_tests.md §4.1, §2A.5, §4.5, §4.11).

Each check returns a list of problems; an empty list passes.

    LY-02  declared sort order holds within each file
    LY-03  min/max statistics on every sort-key column in every row group
    LY-04  no file under 64 MB outside the waivers, unless it is its table's only file in the stage
    C-03   foreign keys resolve (SQL anti-joins over the catalog)
    C-05   coordinate bounds, including origin wrap on circular contigs
    R-06   lineage: every chunk, slot and occurrence reaches an L0 object; no orphan chunks
    R-09   append-only: a later release keeps every earlier stage file byte-identical
    D-01/02/04  `content_hashes`: an order-independent hash of every table's rows
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import duckdb
import pyarrow.parquet as pq

from biobox.registry.schemas import TABLES

NFS_MIN_BYTES = 64 * 2**20


def _sort_ok(rows: list[tuple]) -> int | None:
    """Index of the first row that breaks non-decreasing order, or None."""
    for i in range(1, len(rows)):
        if rows[i] < rows[i - 1]:
            return i
    return None


def check_sort_order(manifest: dict) -> list[str]:
    """LY-02. Sort keys are compared as Python values (bytes compare as unsigned, like Parquet)."""
    root = Path(manifest["root"])
    problems = []
    for f in manifest["files"]:
        spec = TABLES[f["table"]]
        cols = list(spec.sort_key)
        t = pq.read_table(root / f["path"], columns=cols)
        rows = list(zip(*(t.column(c).to_pylist() for c in cols)))
        i = _sort_ok(rows)
        if i is not None:
            problems.append(f"{f['path']}: row {i} breaks sort order {cols}")
    return problems


def check_statistics(manifest: dict) -> list[str]:
    """LY-03."""
    root = Path(manifest["root"])
    problems = []
    for f in manifest["files"]:
        spec = TABLES[f["table"]]
        md = pq.ParquetFile(root / f["path"]).metadata
        names = md.schema.to_arrow_schema().names
        for rg in range(md.num_row_groups):
            if md.row_group(rg).num_rows == 0:
                continue
            for col in spec.sort_key:
                st = md.row_group(rg).column(names.index(col)).statistics
                if st is None or not st.has_min_max:
                    problems.append(f"{f['path']}: row group {rg} has no min/max on {col}")
    return problems


def check_file_sizes(manifest: dict, min_bytes: int = NFS_MIN_BYTES) -> list[str]:
    """LY-04. Waived: registry `small_ok` tables. A table whose stage fits in one small file may be one file."""
    problems = []
    groups: dict[tuple, list[dict]] = {}
    for f in manifest["files"]:
        spec = TABLES[f["table"]]
        if spec.small_ok:
            continue
        stage = f["partitions"].get("stage") or f["partitions"].get("release") or ""
        groups.setdefault((f["table"], stage), []).append(f)
    for (name, stage), files in sorted(groups.items()):
        small = [f for f in files if f["bytes"] < min_bytes]
        if small and len(files) > 1:
            total = sum(f["bytes"] for f in files)
            problems.append(f"{name} {stage}: {len(small)} of {len(files)} files under {min_bytes >> 20} MB "
                            f"({total / 2**20:.1f} MB in all); split only above the band")
    return problems


# Foreign keys that exist at S0 (C-03): (child, child cols, parent, parent cols)
FOREIGN_KEYS = [
    ("contig_occurrence", ("raw_digest",), "seq_chunks", ("raw_digest",)),
    ("contig_occurrence", ("l0_object_sha256",), "l0_manifest", ("object_sha256",)),
    ("contig_index", ("raw_digest",), "seq_chunks", ("raw_digest",)),
    ("seq_locator", ("raw_digest",), "seq_chunks", ("raw_digest",)),
    ("slot_protein", ("contig_uid",), "contig_index", ("contig_uid",)),
    ("slot_protein", ("protein_uid", "gene_call_version"), "protein_seq", ("protein_uid", "gene_call_version")),
    ("protein_occurrence", ("contig_uid", "slot_ordinal", "gene_call_version"),
     "slot_protein", ("contig_uid", "slot_ordinal", "gene_call_version")),
    ("l1_domains", ("protein_uid",), "protein_seq", ("protein_uid",)),
    ("cluster_member", ("member_uid",), "protein_seq", ("protein_uid",)),
    ("cluster_member", ("rep_uid",), "protein_seq", ("protein_uid",)),
    ("cluster_member_aln", ("member_uid", "level", "release"), "cluster_member", ("member_uid", "level", "release")),
]


def _q(cols) -> list[str]:
    return [f'"{c}"' for c in cols]


def check_foreign_keys(con: duckdb.DuckDBPyConnection) -> list[str]:
    """C-03."""
    problems = []
    for child, ccols, parent, pcols in FOREIGN_KEYS:
        on = " AND ".join(f"c.{a} = p.{b}" for a, b in zip(_q(ccols), _q(pcols)))
        n = con.execute(f"SELECT count(*) FROM {child} c ANTI JOIN {parent} p ON {on}").fetchall()[0][0]
        if n:
            problems.append(f"{child}{ccols} -> {parent}{pcols}: {n} unresolved rows")
    return problems


def check_coordinates(con: duckdb.DuckDBPyConnection) -> list[str]:
    """C-05. A slot may end past the contig end only on a circular contig, which genecall marks with
    dist_to_contig_end_bp = contig_len (no edge), and only within one turn."""
    problems = []
    bad = con.execute("""
        SELECT count(*) FROM protein_occurrence o JOIN contig_index ci
          ON o.contig_uid = ci.contig_uid AND ci.is_default_variant
        WHERE NOT (0 <= o.start AND o.start < o."end"
                   AND (o."end" <= ci.contig_len
                        OR (o.start < ci.contig_len AND o."end" <= o.start + ci.contig_len
                            AND o.dist_to_contig_end_bp = ci.contig_len)))
    """).fetchall()[0][0]
    if bad:
        problems.append(f"protein_occurrence: {bad} slots outside contig bounds")
    bad = con.execute("""
        SELECT count(*) FROM l1_domains d JOIN protein_seq p USING (protein_uid)
        WHERE NOT (d.env_start < d.env_end AND d.env_end <= length(replace(p.seq, '*', ''))
                   AND d.ali_start >= d.env_start AND d.ali_end <= d.env_end
                   AND d.hmm_from < d.hmm_to AND d.hmm_to <= d.model_len)
    """).fetchall()[0][0]
    if bad:
        problems.append(f"l1_domains: {bad} rows outside protein or model bounds")
    bad = con.execute("""
        SELECT count(*) FROM cluster_member_aln a
        JOIN protein_seq m ON m.protein_uid = a.member_uid
        JOIN protein_seq r ON r.protein_uid = a.rep_uid
        WHERE NOT (a.qstart < a.qend AND a.qend <= m.length AND a.tstart < a.tend AND a.tend <= r.length)
    """).fetchall()[0][0]
    if bad:
        problems.append(f"cluster_member_aln: {bad} rows outside sequence bounds")
    return problems


def check_lineage(con: duckdb.DuckDBPyConnection) -> list[str]:
    """R-06: no orphan chunks; every slot's contig has an occurrence row; every occurrence an L0 object."""
    problems = []
    checks = {
        "seq_chunks without an occurrence":
            "SELECT count(DISTINCT raw_digest) FROM seq_chunks s ANTI JOIN contig_occurrence o USING (raw_digest)",
        "slots whose contig has no occurrence":
            "SELECT count(*) FROM slot_protein s ANTI JOIN contig_occurrence o USING (contig_uid)",
        "occurrences without an L0 object":
            "SELECT count(*) FROM contig_occurrence o ANTI JOIN l0_manifest m ON o.l0_object_sha256 = m.object_sha256",
    }
    for what, sql in checks.items():
        n = con.execute(sql).fetchall()[0][0]
        if n:
            problems.append(f"{what}: {n}")
    return problems


def check_append_only(before: dict, after: dict) -> list[str]:
    """R-09: every non-metadata file of the earlier manifest is in the later one with the same sha256."""
    later = {f["path"]: f["sha256"] for f in after["files"]}
    problems = []
    for f in before["files"]:
        if TABLES[f["table"]].small_ok:
            continue  # small metadata files (seq_locator, l0_manifest) are rewritten per stage by design
        if later.get(f["path"]) != f["sha256"]:
            problems.append(f"{f['path']}: {'removed' if f['path'] not in later else 'rewritten'}")
    return problems


def content_hashes(manifest: dict) -> dict[str, str]:
    """D-01/D-02/D-04: per table, sha256 over the sorted per-row hashes (partition values included).

    Independent of file layout, row order and thread count; equal hashes mean equal row multisets.
    """
    root = Path(manifest["root"])
    per: dict[str, list[bytes]] = {}
    for f in manifest["files"]:
        t = pq.read_table(root / f["path"])
        part = repr(sorted(f["partitions"].items())).encode()
        cols = [t.column(c).to_pylist() for c in t.column_names]
        rows = per.setdefault(f["table"], [])
        for r in zip(*cols):
            rows.append(hashlib.sha256(part + repr(r).encode()).digest())
    return {name: hashlib.sha256(b"".join(sorted(rows))).hexdigest() for name, rows in sorted(per.items())}
