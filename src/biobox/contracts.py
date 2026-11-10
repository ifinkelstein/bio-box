"""Layout and contract checks over a release (data_layer_and_tests.md §4.1, §2A.5, §4.5, §4.11).

Each check returns a list of problems; an empty list passes.

    LY-02  declared sort order holds within each file (each shard too; none across shards)
    LY-03  min/max statistics on every sort-key column in every row group
    LY-04  no file under 64 MB outside the waivers, unless it is its table's only file in the stage
    C-03   foreign keys resolve (SQL anti-joins over the catalog); `check_shards`: seq_locator and
           the shards' `.fai` agree record for record, and each `.gzi` lists the shard's blocks
    C-04   parallel LIST columns of a row have equal lengths (tapes: = n_slots; tokens: = the tape's)
    C-05   coordinate bounds, including origin wrap on circular contigs
    C-07   saturating fields carry their overflow flag; clipped share < 0.1%
    R-06   lineage: every stored sequence, slot and occurrence reaches an L0 object
    R-09   append-only: a later release keeps every earlier stage file byte-identical
    D-01/02/04  `content_hashes`: an order-independent hash of every table's rows
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import duckdb
import pyarrow.parquet as pq

from biobox.l0seq.bgzf import BgzfFasta, block_starts, read_fai, read_gzi
from biobox.registry.schemas import ARTIFACTS, TABLES

NFS_MIN_BYTES = 64 * 2**20


def _sort_ok(rows: list[tuple]) -> int | None:
    """Index of the first row that breaks non-decreasing order (nulls first), or None."""
    keyed = [tuple((0,) if v is None else (1, v) for v in r) for r in rows]
    for i in range(1, len(keyed)):
        if keyed[i] < keyed[i - 1]:
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
        leaf = {md.schema.column(i).path: i for i in range(md.num_columns)}  # leaf columns, lists flattened
        for rg in range(md.num_row_groups):
            if md.row_group(rg).num_rows == 0:
                continue
            for col in spec.sort_key:
                st = md.row_group(rg).column(leaf[col]).statistics
                if st is not None and st.has_null_count and st.null_count == md.row_group(rg).num_rows:
                    continue  # an all-null column has no min/max; its null count prunes it
                if st is None or not st.has_min_max:
                    problems.append(f"{f['path']}: row group {rg} has no min/max on {col}")
    return problems


def check_file_sizes(manifest: dict, min_bytes: int = NFS_MIN_BYTES) -> list[str]:
    """LY-04. Waived: registry `small_ok` tables and shard indexes. A table whose stage fits in one
    small file may be one file."""
    problems = []
    groups: dict[tuple, list[dict]] = {}
    for f in manifest["files"]:
        spec = TABLES[f["table"]]
        if spec.small_ok:
            continue
        stage = f["partitions"].get("stage") or f["partitions"].get("release") or ""
        groups.setdefault((f["table"], stage), []).append(f)
    for a in manifest.get("artifacts", []):
        if any(a["path"].endswith(c) for c in ARTIFACTS[a["artifact"]].companions):
            continue
        groups.setdefault((a["artifact"], a["partitions"].get("stage", "")), []).append(a)
    for (name, stage), files in sorted(groups.items()):
        small = [f for f in files if f["bytes"] < min_bytes]
        if small and len(files) > 1:
            total = sum(f["bytes"] for f in files)
            problems.append(f"{name} {stage}: {len(small)} of {len(files)} files under {min_bytes >> 20} MB "
                            f"({total / 2**20:.1f} MB in all); split only above the band")
    return problems


# Foreign keys that exist at S0 (C-03): (child, child cols, parent, parent cols). Slots are named by
# feature_uid, never by an ordinal (AGENTS §3.3; card B-17); a gene-call version is the producer of its
# CDS features.
FOREIGN_KEYS = [
    ("contig_occurrence", ("raw_digest",), "seq_locator", ("raw_digest",)),
    ("contig_occurrence", ("l0_object_sha256",), "l0_manifest", ("object_sha256",)),
    ("contig_index", ("raw_digest",), "seq_locator", ("raw_digest",)),
    ("contig_dim", ("contig_uid",), "contig_index", ("contig_uid",)),
    ("feature", ("contig_uid",), "contig_index", ("contig_uid",)),
    ("feature_protein", ("feature_uid", "gene_call_version"), "feature", ("feature_uid", "producer")),
    ("feature_protein", ("protein_uid", "gene_call_version"), "protein_seq", ("protein_uid", "gene_call_version")),
    ("protein_occurrence", ("feature_uid", "gene_call_version"),
     "feature_protein", ("feature_uid", "gene_call_version")),
    ("l1_domains", ("protein_uid",), "protein_seq", ("protein_uid",)),
    ("cluster_member", ("member_uid",), "protein_seq", ("protein_uid",)),
    ("cluster_member", ("rep_uid",), "protein_seq", ("protein_uid",)),
    ("cluster_member_aln", ("member_uid", "level", "release"), "cluster_member", ("member_uid", "level", "release")),
    ("igr", ("contig_uid",), "contig_index", ("contig_uid",)),
    ("contig_tapes", ("contig_uid",), "contig_index", ("contig_uid",)),
    ("tape_tokens", ("contig_uid", "tape_version"), "contig_tapes", ("contig_uid", "tape_version")),
    ("token_postings", ("token_id", "release"), "token_dict", ("token_id", "release")),
    ("token_postings", ("contig_uid",), "contig_tapes", ("contig_uid",)),
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


# Search outputs (C-03, S-05): run files resolve to their run, hits to their targets. Over the
# catalog's search views (catalog.attach_search), where `run` is the run directory's id.
SEARCH_LINKS = [
    ("search_hit -> search_run",
     "SELECT count(*) FROM search_hit h ANTI JOIN search_run r ON h.run_id = r.run_id"),
    ("search_profile -> search_run",
     "SELECT count(*) FROM search_profile p ANTI JOIN search_run r ON p.run_id = r.run_id"),
    ("neighborhoods -> search_run",
     "SELECT count(*) FROM neighborhoods n ANTI JOIN search_run r ON n.run = r.run"),
    ("run_id equals its directory",
     "SELECT (SELECT count(*) FROM search_run WHERE hex(run_id) <> upper(run)) "
     "+ (SELECT count(*) FROM search_hit WHERE hex(run_id) <> upper(run)) "
     "+ (SELECT count(*) FROM search_profile WHERE hex(run_id) <> upper(run))"),
    ("protein hits -> protein_seq",
     "SELECT count(*) FROM search_hit h WHERE h.target_kind = 'protein' "
     "AND NOT EXISTS (SELECT 1 FROM protein_seq p WHERE p.protein_uid = h.target_uid)"),
    ("region hits -> search_region",
     "SELECT count(*) FROM search_hit h WHERE h.target_kind IN ('dark_orf', 'igr') "
     "AND NOT EXISTS (SELECT 1 FROM search_region s WHERE s.region_uid = h.target_uid)"),
    ("neighborhood anchors -> final hits",
     "SELECT count(*) FROM neighborhoods n JOIN search_run r ON n.run = r.run WHERE NOT EXISTS ("
     "SELECT 1 FROM search_hit h WHERE h.run = n.run AND h.round = r.final_round "
     "AND h.curation_label = 'keep' AND h.target_uid = n.anchor_protein_uid)"),
]


def check_search_links(con: duckdb.DuckDBPyConnection) -> list[str]:
    """C-03 for search runs: every row resolves (catalog opened with its data root)."""
    problems = []
    for name, sql in SEARCH_LINKS:
        n = con.execute(sql).fetchall()[0][0]
        if n:
            problems.append(f"{name}: {n} unresolved rows")
    return problems


def check_shards(manifest: dict) -> list[str]:
    """C-03 for the shards, and LY-02 inside each.

    Every seq_locator row names a shard of the manifest whose `.fai` has a record
    `raw_digest.hex()` of the same length, and every `.fai` record is a locator row for that
    shard. Records are sorted by name within a shard; there is no order across shards (one
    ingest batch writes one shard, card B-15). Each `.gzi` lists exactly the block starts found
    by walking the block headers.
    """
    root = Path(manifest["root"])
    problems = []
    shards = sorted(a["path"] for a in manifest.get("artifacts", [])
                    if a["artifact"] == "seq_shards" and a["path"].endswith(".fa.gz"))
    loc_files = [f for f in manifest["files"] if f["table"] == "seq_locator"]
    locator: dict[str, int] = {}
    by_shard: dict[str, dict[str, int]] = {}
    for f in loc_files:
        t = pq.read_table(root / f["path"])
        for rd, shard, n in zip(*(t.column(c).to_pylist() for c in ("raw_digest", "shard", "contig_len"))):
            by_shard.setdefault(shard, {})[rd.hex()] = n
            locator[rd.hex()] = n
    for shard in sorted(set(by_shard) - set(shards)):
        problems.append(f"seq_locator names {shard}, which is not in the manifest ({len(by_shard[shard])} rows)")
    for shard in shards:
        fai = read_fai(root / (shard + ".fai"))
        bg = BgzfFasta(root / shard)
        if not bg.check_sorted():  # the readers bisect the .fai, so its names must increase
            problems.append(f"{shard}: records not sorted by raw_digest (LY-02)")
        bg.close()
        want = by_shard.get(shard, {})
        got = {e.name: e.length for e in fai}
        if got != want:
            extra, missing = set(got) - set(want), set(want) - set(got)
            wrong = {k for k in set(got) & set(want) if got[k] != want[k]}
            problems.append(f"{shard}: .fai vs seq_locator: {len(extra)} unlisted, {len(missing)} missing, "
                            f"{len(wrong)} lengths differ")
        if read_gzi(root / (shard + ".gzi")) != block_starts(root / shard)[1:]:
            problems.append(f"{shard}: .gzi does not list the shard's blocks")
    return problems


def check_list_lengths(con: duckdb.DuckDBPyConnection) -> list[str]:
    """C-04."""
    problems = []
    for name, spec in sorted(TABLES.items()):
        if not spec.list_columns:
            continue
        first = spec.list_columns[0]
        conds = [f'("{c}" IS NULL OR len("{c}") = len("{first}"))' for c in spec.list_columns[1:]]
        if "n_slots" in spec.schema.names:
            conds.append(f'len("{first}") = n_slots')
        n = con.execute(f"SELECT count(*) FROM {name} WHERE NOT ({' AND '.join(conds)})").fetchall()[0][0]
        if n:
            problems.append(f"{name}: {n} rows with parallel lists of different lengths")
    n = con.execute("""
        SELECT count(*) FROM tape_tokens k JOIN contig_tapes t USING (contig_uid, tape_version)
        WHERE len(k.family_tokens) <> t.n_slots OR len(k.nr50_tokens) <> t.n_slots OR len(k.nr90_tokens) <> t.n_slots
    """).fetchall()[0][0]
    if n:
        problems.append(f"tape_tokens: {n} rows whose lists differ in length from their tape")
    return problems


SATURATION_MAX_SHARE = 0.001


def check_saturation(con: duckdb.DuckDBPyConnection) -> list[str]:
    """C-07 for the tapes' uint16 censoring distances (slot_flags bit 0 marks a clipped value)."""
    problems = []
    bad, clipped, total = con.execute("""
        SELECT count(*) FILTER (WHERE f & 1 = 1 AND NOT (s = 65535 OR e = 65535)),
               count(*) FILTER (WHERE f & 1 = 1), count(*)
        FROM (SELECT unnest(slot_flags) f, unnest(orfs_to_contig_starts) s, unnest(orfs_to_contig_ends) e
              FROM contig_tapes WHERE orfs_to_contig_starts IS NOT NULL)
    """).fetchall()[0]
    if bad:
        problems.append(f"contig_tapes: {bad} slots flagged clipped whose distances are not saturated")
    if total and clipped / total > SATURATION_MAX_SHARE:
        problems.append(f"contig_tapes: {clipped} of {total} censoring distances clipped (> 0.1%)")
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
        SELECT count(*) FROM igr i JOIN contig_index ci
          ON i.contig_uid = ci.contig_uid AND ci.is_default_variant
        WHERE NOT (0 <= i.start AND i.start < i."end" AND i.start < ci.contig_len
                   AND i."end" <= i.start + ci.contig_len AND i.strand = 1)
    """).fetchall()[0][0]
    if bad:
        problems.append(f"igr: {bad} regions outside contig bounds")
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
    bad = con.execute("""
        SELECT count(*) FROM (
          SELECT contig_len, is_circular, unnest(start_bps) s, unnest(end_bps) e, unnest(strands) st
          FROM contig_tapes)
        WHERE NOT (s < e AND s < contig_len AND st IN (1, -1)
                   AND (e <= contig_len OR (is_circular AND e <= s + contig_len)))
    """).fetchall()[0][0]
    if bad:
        problems.append(f"contig_tapes: {bad} slots outside contig bounds")
    # the pair macros read slot positions from the postings (card B-11): each ordinal lies in its tape
    bad = con.execute("""
        SELECT count(*) FROM token_postings p JOIN contig_tapes t USING (contig_uid, tape_version)
        WHERE NOT (p.ordinal < t.n_slots AND p.n_slots = t.n_slots)
    """).fetchall()[0][0]
    if bad:
        problems.append(f"token_postings: {bad} rows with ordinals outside their tape")
    return problems


def check_lineage(con: duckdb.DuckDBPyConnection) -> list[str]:
    """R-06: no orphan sequence; every feature's contig has an occurrence row; every occurrence an L0 object."""
    problems = []
    checks = {
        "stored sequences without an occurrence":
            "SELECT count(*) FROM seq_locator s ANTI JOIN contig_occurrence o USING (raw_digest)",
        "features whose contig has no occurrence":
            "SELECT count(*) FROM feature f ANTI JOIN contig_occurrence o USING (contig_uid)",
        "occurrences without an L0 object":
            "SELECT count(*) FROM contig_occurrence o ANTI JOIN l0_manifest m ON o.l0_object_sha256 = m.object_sha256",
    }
    for what, sql in checks.items():
        n = con.execute(sql).fetchall()[0][0]
        if n:
            problems.append(f"{what}: {n}")
    return problems


def check_append_only(before: dict, after: dict) -> list[str]:
    """R-09: every non-metadata file and artifact of the earlier manifest is in the later one with the same sha256."""
    later = {f["path"]: f["sha256"] for f in after["files"] + after.get("artifacts", [])}
    problems = []
    for f in before["files"] + before.get("artifacts", []):
        if "table" in f and TABLES[f["table"]].small_ok:
            continue  # small metadata files (l0_manifest) are rewritten per stage by design
        if later.get(f["path"]) != f["sha256"]:
            problems.append(f"{f['path']}: {'removed' if f['path'] not in later else 'rewritten'}")
    return problems


def content_hashes(manifest: dict) -> dict[str, str]:
    """D-01/D-02/D-04: per table, sha256 over the sorted per-row hashes (partition values included).

    Independent of file layout, row order and thread count; equal hashes mean equal row multisets.
    Artifacts are hashed by (path, sha256): shard bytes depend only on the set of sequences.
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
    for a in manifest.get("artifacts", []):
        per.setdefault(a["artifact"], []).append(hashlib.sha256(f"{a['path']}\t{a['sha256']}".encode()).digest())
    return {name: hashlib.sha256(b"".join(sorted(rows))).hexdigest() for name, rows in sorted(per.items())}
