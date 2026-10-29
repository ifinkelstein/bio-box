"""Schema registry, v0.

One `TableSpec` per planned Parquet table: Arrow schema of the file columns,
keys, sort order, hive partition keys, enum domains and whether the table is a
small metadata file exempt from the NFS >=64 MB rule (LY-04).

Partition keys live only in the path (`gene_call_version=v1/...`), never as file
columns, so a reader never sees two copies that could disagree. A primary key
may include a partition key; readers get it from the path.

Sources: planning/data_layer_and_tests.md §1, §2, §2A.1, §2B.1;
planning/tests/TEST_ARCHITECTURE.md §7, §8. Tables planned but not yet
specified here are listed in tests/registry/pending_tables.yaml.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pyarrow as pa

UID = pa.binary(16)  # BLAKE3-128 content UIDs
SHA256 = pa.binary(32)  # raw_digest and object checksums
STRAND = pa.int8()  # +1 / -1
COORD = pa.int64()  # 0-based half-open contig coordinates

STRAND_VALUES = (1, -1)


@dataclass(frozen=True)
class TableSpec:
    name: str
    layer: str
    schema: pa.Schema
    primary_key: tuple[str, ...]
    sort_key: tuple[str, ...]
    partition_keys: tuple[str, ...] = ()
    enums: dict[str, tuple] = field(default_factory=dict)
    small_ok: bool = False
    doc: str = ""


def _f(name: str, typ: pa.DataType, nullable: bool = False) -> pa.Field:
    return pa.field(name, typ, nullable=nullable)


_SPECS = [
    TableSpec(
        name="l0_manifest",
        layer="L0",
        schema=pa.schema([
            _f("source_uid", pa.uint64()),
            _f("source_name", pa.string()),
            _f("source_release", pa.string()),
            _f("object_path", pa.string()),
            _f("object_sha256", SHA256),
            _f("bytes", pa.int64()),
            _f("license", pa.string()),
            _f("redistributable", pa.bool_()),
            _f("retrieved_at", pa.timestamp("s", tz="UTC")),
        ]),
        primary_key=("source_uid", "object_path"),
        sort_key=("source_uid", "object_path"),
        small_ok=True,
        doc="data_layer_and_tests.md §1 (L0)",
    ),
    TableSpec(
        name="seq_chunks",
        layer="L0-seq",
        schema=pa.schema([
            _f("raw_digest", SHA256),
            _f("chunk_ix", pa.uint32()),
            _f("seq", pa.string()),
            _f("contig_len", COORD),
        ]),
        primary_key=("raw_digest", "chunk_ix"),
        sort_key=("raw_digest", "chunk_ix"),
        partition_keys=("stage", "digest_prefix"),
        doc="data_layer_and_tests.md §2A.1",
    ),
    TableSpec(
        name="contig_occurrence",
        layer="L0-seq",
        schema=pa.schema([
            _f("source_uid", pa.uint64()),
            _f("source_release", pa.string()),
            _f("l0_object_sha256", SHA256),
            _f("record_ordinal", pa.uint64()),
            _f("accession", pa.string()),
            _f("original_header", pa.string()),
            _f("raw_digest", SHA256),
            _f("contig_uid", UID),
            _f("line_width", pa.uint32()),
        ]),
        primary_key=("l0_object_sha256", "record_ordinal"),
        sort_key=("l0_object_sha256", "record_ordinal"),
        partition_keys=("stage", "source"),
        doc="data_layer_and_tests.md §2A.1",
    ),
    TableSpec(
        name="contig_index",
        layer="L0-seq",
        schema=pa.schema([
            _f("contig_uid", UID),
            _f("raw_digest", SHA256),
            _f("contig_len", COORD),
            _f("is_default_variant", pa.bool_()),
        ]),
        primary_key=("contig_uid", "raw_digest"),
        sort_key=("contig_uid", "raw_digest"),
        doc="data_layer_and_tests.md §2A.3 (contig_uid -> raw_digest; default = earliest occurrence by precedence)",
    ),
    TableSpec(
        name="seq_locator",
        layer="L0-seq",
        schema=pa.schema([
            _f("raw_digest", SHA256),
            _f("file", pa.string()),
            _f("row_group", pa.uint32()),
            _f("first_chunk_row", pa.uint64()),
        ]),
        primary_key=("raw_digest",),
        sort_key=("raw_digest",),
        small_ok=True,
        doc="data_layer_and_tests.md §2A.1",
    ),
    TableSpec(
        name="protein_seq",
        layer="L0-seq",
        schema=pa.schema([
            _f("protein_uid", UID),
            _f("protein_digest", SHA256),
            _f("length", pa.uint32()),
            _f("seq", pa.string()),
        ]),
        primary_key=("protein_uid", "gene_call_version"),
        sort_key=("protein_uid",),
        partition_keys=("gene_call_version", "stage", "protein_prefix"),
        doc="data_layer_and_tests.md §2A.1",
    ),
    TableSpec(
        name="slot_protein",
        layer="L0-seq",
        schema=pa.schema([
            _f("contig_uid", UID),
            _f("slot_ordinal", pa.uint32()),
            _f("gene_call_version", pa.string()),
            _f("protein_uid", UID),
            _f("genetic_code", pa.uint8()),
            _f("edit_type", pa.string()),
        ]),
        primary_key=("contig_uid", "slot_ordinal", "gene_call_version"),
        sort_key=("contig_uid", "slot_ordinal"),
        enums={"edit_type": ("none", "frameshift", "selenocysteine", "pyrrolysine", "intron")},
        doc="data_layer_and_tests.md §2A.1; slot_uid is the composite (contig_uid, slot_ordinal)",
    ),
    TableSpec(
        name="igr",
        layer="L0-seq",
        schema=pa.schema([
            _f("contig_uid", UID),
            _f("start", COORD),
            _f("end", COORD),
            _f("strand", STRAND),
            _f("left_slot", pa.uint32(), nullable=True),
            _f("right_slot", pa.uint32(), nullable=True),
        ]),
        primary_key=("contig_uid", "start", "end", "strand", "igr_version"),
        sort_key=("contig_uid", "start"),
        partition_keys=("igr_version", "stage"),
        enums={"strand": STRAND_VALUES},
        doc="data_layer_and_tests.md §2A.1 (coordinates only; null slot at a contig end)",
    ),
    TableSpec(
        name="genome_identity",
        layer="L0.5",
        schema=pa.schema([
            _f("assembly_accession", pa.string()),
            _f("genome_uid", UID),
            _f("paired_accession", pa.string(), nullable=True),
            _f("biosample", pa.string(), nullable=True),
            _f("source_uid", pa.uint64()),
            _f("precedence_rank", pa.uint16()),
            _f("member_of", UID, nullable=True),
            _f("containment_to_rep", pa.float32(), nullable=True),
            _f("fcs_gx_status", pa.string()),
            _f("identity_policy_hash", pa.string()),
            _f("release_id", pa.string()),
        ]),
        primary_key=("assembly_accession", "release_id"),
        sort_key=("genome_uid", "assembly_accession"),
        enums={"fcs_gx_status": ("ncbi_screened", "screened_here", "not_screened")},
        doc="data_layer_and_tests.md §1 (L0.5); fcs_gx_status from preflight_gaps.md #4 (2026-10-29)",
    ),
    TableSpec(
        name="contig_member",
        layer="L0.5",
        schema=pa.schema([
            _f("member_contig_uid", UID),
            _f("rep_contig_uid", UID),
            _f("method", pa.string()),
            _f("ani", pa.float32(), nullable=True),
            _f("containment", pa.float32(), nullable=True),
            _f("release_id", pa.string()),
        ]),
        primary_key=("member_contig_uid", "release_id"),
        sort_key=("rep_contig_uid", "member_contig_uid"),
        enums={"method": ("self", "exact", "containment", "ani")},
        doc="data_layer_and_tests.md §2A.1; representatives carry a self row",
    ),
    TableSpec(
        name="protein_occurrence",
        layer="L1-derived",
        schema=pa.schema([
            _f("protein_uid", UID),
            _f("contig_uid", UID),
            _f("slot_ordinal", pa.uint32()),
            _f("genome_uid", UID),
            _f("start", COORD),
            _f("end", COORD),
            _f("strand", STRAND),
            _f("tax_order", pa.uint32(), nullable=True),
            _f("is_species_rep", pa.bool_()),
            _f("tier", pa.uint8()),
            _f("assembly_quality", pa.string()),
            _f("dist_to_contig_end_bp", pa.uint32()),
        ]),
        primary_key=("contig_uid", "slot_ordinal", "gene_call_version"),
        sort_key=("protein_uid",),
        partition_keys=("gene_call_version", "stage"),
        enums={
            "strand": STRAND_VALUES,
            "assembly_quality": ("complete", "chromosome", "scaffold", "contig", "mag_hq", "mag_mq", "vOTU", "plasmid"),
        },
        doc="data_layer_and_tests.md §2B.1",
    ),
    TableSpec(
        name="cluster_member",
        layer="L2",
        schema=pa.schema([
            _f("member_uid", UID),
            _f("level", pa.string()),
            _f("cluster_uid", UID),
            _f("rep_uid", UID),
            _f("release_id", pa.string()),
        ]),
        primary_key=("member_uid", "level", "release_id"),
        sort_key=("level", "cluster_uid", "member_uid"),
        enums={"level": ("nr90", "nr50", "subfam")},
        doc="data_layer_and_tests.md §1 (L2), §2 (cluster_uid)",
    ),
    TableSpec(
        name="cluster_member_aln",
        layer="L2",
        schema=pa.schema([
            _f("member_uid", UID),
            _f("level", pa.string()),
            _f("rep_uid", UID),
            _f("qstart", pa.uint32()),
            _f("qend", pa.uint32()),
            _f("tstart", pa.uint32()),
            _f("tend", pa.uint32()),
            _f("pident", pa.float32()),
            _f("alnlen", pa.uint32()),
            _f("cigar", pa.binary()),
            _f("release_id", pa.string()),
        ]),
        primary_key=("member_uid", "level", "release_id"),
        sort_key=("level", "rep_uid", "member_uid"),
        enums={"level": ("nr90", "nr50", "subfam")},
        doc="data_layer_and_tests.md §2B.1; representatives carry a self row",
    ),
    TableSpec(
        name="id_history",
        layer="L2",
        schema=pa.schema([
            _f("entity_type", pa.string()),
            _f("old_uid", UID),
            _f("new_uid", UID),
            _f("relation", pa.string()),
            _f("from_release", pa.string()),
            _f("to_release", pa.string()),
            _f("evidence", pa.string(), nullable=True),
        ]),
        primary_key=("entity_type", "old_uid", "new_uid", "from_release", "to_release"),
        sort_key=("entity_type", "old_uid"),
        enums={
            "entity_type": ("genome", "contig", "protein", "segment", "cluster", "token"),
            "relation": ("same", "split", "merged", "retired", "renamed"),
        },
        small_ok=True,
        doc="data_layer_and_tests.md §2; for `retired`, new_uid equals old_uid",
    ),
    TableSpec(
        name="search_run",
        layer="search",
        schema=pa.schema([
            _f("run_id", UID),
            _f("recipe_hash", pa.string()),
            _f("recipe_yaml", pa.string()),
            _f("argv", pa.list_(pa.string())),
            _f("tool_versions_json", pa.string()),
            _f("runtime_image_digest", pa.string()),
            _f("written_against", pa.string()),
            _f("ran_against", pa.string()),
            _f("pack_id", UID),
            _f("shard_sha256", pa.list_(SHA256)),
            _f("z_prot", pa.float64()),
            _f("z_nt_mb", pa.float64()),
            _f("query_hash", pa.string()),
            _f("submitted_at", pa.timestamp("s", tz="UTC")),
            _f("wall_s", pa.float64(), nullable=True),
            _f("peak_rss_bytes", pa.int64(), nullable=True),
            _f("cpu_min_max_pid", pa.float64(), nullable=True),
            _f("host", pa.string()),
            _f("exit_code", pa.int32(), nullable=True),
            _f("status", pa.string()),
            _f("partial_policy", pa.string()),
            _f("max_output_gb", pa.float32()),
        ]),
        primary_key=("run_id",),
        sort_key=("run_id",),
        enums={
            "status": ("ok", "failed", "timeout", "oom", "stopped", "output_cap"),
            "partial_policy": ("exclude", "include", "flag"),
        },
        small_ok=True,
        doc="data_layer_and_tests.md §2B.1; SEARCH_PLAN SD-20",
    ),
    TableSpec(
        name="search_hit",
        layer="search",
        schema=pa.schema([
            _f("run_id", UID),
            _f("round", pa.uint16()),
            _f("target_kind", pa.string()),
            _f("target_uid", UID),
            _f("bitscore", pa.float32()),
            _f("evalue_zref", pa.float64()),
            _f("target_size", pa.uint64()),
            _f("env_start", pa.uint32()),
            _f("env_end", pa.uint32()),
            _f("ali_start", pa.uint32()),
            _f("ali_end", pa.uint32()),
            _f("hmm_from", pa.uint32(), nullable=True),
            _f("hmm_to", pa.uint32(), nullable=True),
            _f("model_len", pa.uint32(), nullable=True),
            _f("is_partial", pa.bool_()),
            _f("first_inclusion_round", pa.uint16()),
            _f("score_kind", pa.string()),
            _f("curation_label", pa.string(), nullable=True),
            _f("curation_reason", pa.string(), nullable=True),
        ]),
        primary_key=("run_id", "round", "target_kind", "target_uid", "env_start"),
        sort_key=("round", "target_kind", "target_uid"),
        enums={
            "target_kind": ("protein", "segment", "dark_orf", "igr", "contig_region", "spacer"),
            "score_kind": ("exact", "projected"),
        },
        small_ok=True,
        doc="data_layer_and_tests.md §2B.1; profile coordinates per SEARCH_PLAN SD-17",
    ),
    TableSpec(
        name="test_run",
        layer="tests",
        schema=pa.schema([
            _f("test_id", pa.string()),
            _f("ring", pa.string()),
            _f("operation", pa.string(), nullable=True),
            _f("workflow", pa.string(), nullable=True),
            _f("params_hash", pa.string()),
            _f("outcome", pa.string()),
            _f("metrics_json", pa.string(), nullable=True),
            _f("duration_s", pa.float64()),
            _f("started_at", pa.timestamp("ms", tz="UTC")),
            _f("code_commit", pa.string()),
            _f("runtime_image_digest", pa.string(), nullable=True),
            _f("test_image_digest", pa.string(), nullable=True),
            _f("release_id", pa.string(), nullable=True),
            _f("host", pa.string()),
            _f("seed", pa.int64(), nullable=True),
        ]),
        primary_key=("test_id", "params_hash", "started_at", "host"),
        sort_key=("started_at", "test_id"),
        enums={
            "ring": ("R0", "R1", "R2", "R3", "R4", "R5"),
            "outcome": ("pass", "fail", "error", "skip"),
        },
        small_ok=True,
        doc="TEST_ARCHITECTURE.md §8",
    ),
]

TABLES: dict[str, TableSpec] = {s.name: s for s in _SPECS}


def table(name: str) -> TableSpec:
    return TABLES[name]


def conformance_problems(data: pa.Table, spec: TableSpec) -> list[str]:
    """Contract check for one Arrow table against its spec (C-01, C-02, C-06, LY-01).

    Returns human-readable problems; an empty list means the table conforms.
    """
    problems: list[str] = []
    if not data.schema.equals(spec.schema, check_metadata=False):
        problems.append(f"{spec.name}: schema differs from registry")
        return problems
    for col, allowed in spec.enums.items():
        values = set(data.column(col).to_pylist())
        values.discard(None)
        bad = values - set(allowed)
        if bad:
            problems.append(f"{spec.name}.{col}: values outside enum: {sorted(map(str, bad))}")
    file_key = [c for c in spec.primary_key if c not in spec.partition_keys]
    if data.num_rows:
        keys = data.select(file_key)
        if keys.group_by(file_key).aggregate([]).num_rows != data.num_rows:
            problems.append(f"{spec.name}: duplicate primary key")
    return problems
