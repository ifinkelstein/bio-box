"""Schema registry, v0.

One `TableSpec` per planned Parquet table: Arrow schema of the file columns,
keys, sort order, hive partition keys, enum domains and whether the table is a
small metadata file exempt from the NFS >=64 MB rule (LY-04), and where its
files live under the data root (`location`, read by the manifest builder).

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
    # Parallel LIST columns: one entry per slot, so equal lengths in every row (C-04).
    list_columns: tuple[str, ...] = ()
    # Path under the data root: a directory holding `<k>=<v>/.../*.parquet` in
    # partition_keys order, or a single `.parquet` file. Empty: not yet placed
    # (the catalog exposes an empty view).
    location: str = ""


def _f(name: str, typ: pa.DataType, nullable: bool = False) -> pa.Field:
    return pa.field(name, typ, nullable=nullable)


_SPECS = [
    TableSpec(
        name="l0_manifest",
        location="l0/l0_manifest.parquet",
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
            _f("retrieved_at", pa.timestamp("ms", tz="UTC")),
        ]),
        primary_key=("source_uid", "object_path"),
        sort_key=("source_uid", "object_path"),
        small_ok=True,
        doc="data_layer_and_tests.md §1 (L0)",
    ),
    TableSpec(
        name="contig_occurrence",
        location="l0seq/contig_occurrence",
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
        location="l0seq/contig_index",
        layer="L0-seq",
        schema=pa.schema([
            _f("contig_uid", UID),
            _f("raw_digest", SHA256),
            _f("contig_len", COORD),
            _f("is_default_variant", pa.bool_()),
        ]),
        primary_key=("contig_uid", "raw_digest"),
        sort_key=("contig_uid", "raw_digest"),
        partition_keys=("stage",),
        doc=(
            "data_layer_and_tests.md §2A.3 (contig_uid -> raw_digest). Default variant: lowest "
            "(source precedence, l0_object_sha256, record_ordinal), so it does not depend on input order (MR-01)"
        ),
    ),
    TableSpec(
        name="seq_locator",
        location="l0seq/seq_locator.parquet",
        layer="L0-seq",
        schema=pa.schema([
            _f("raw_digest", SHA256),
            _f("shard", pa.string()),  # root-relative path of the .fa.gz; the record name is raw_digest hex
            _f("contig_len", COORD),
        ]),
        primary_key=("raw_digest",),
        sort_key=("raw_digest",),
        small_ok=True,
        doc="data_layer_and_tests.md §2A.1; DEC-17",
    ),
    TableSpec(
        name="protein_seq",
        location="l0seq/protein_seq",
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
        location="l1/slot_protein",
        layer="L0-seq",
        schema=pa.schema([
            _f("contig_uid", UID),
            _f("slot_ordinal", pa.uint32()),
            _f("protein_uid", UID),
            _f("genetic_code", pa.uint8()),
            _f("edit_type", pa.string()),
        ]),
        primary_key=("contig_uid", "slot_ordinal", "gene_call_version"),
        sort_key=("contig_uid", "slot_ordinal"),
        partition_keys=("gene_call_version", "stage"),
        enums={"edit_type": ("none", "frameshift", "selenocysteine", "pyrrolysine", "intron")},
        doc="data_layer_and_tests.md §2A.1; slot_uid is the composite (contig_uid, slot_ordinal)",
    ),
    TableSpec(
        name="igr",
        location="l0seq/igr",
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
        location="l1/protein_occurrence",
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
        name="l1_domains",
        location="l1/domains",
        layer="L1",
        schema=pa.schema([
            _f("protein_uid", UID),
            _f("model_acc", pa.string()),
            _f("model_name", pa.string()),
            _f("model_len", pa.uint32()),
            _f("seq_bitscore", pa.float32()),
            _f("seq_evalue_zref", pa.float64()),
            _f("dom_bitscore", pa.float32()),
            _f("dom_ievalue_zref", pa.float64()),
            _f("env_start", pa.uint32()),
            _f("env_end", pa.uint32()),
            _f("ali_start", pa.uint32()),
            _f("ali_end", pa.uint32()),
            _f("hmm_from", pa.uint32()),
            _f("hmm_to", pa.uint32()),
            _f("cutoff", pa.string()),
            _f("z", pa.float64()),
            _f("dom_z", pa.float64()),
        ]),
        primary_key=("protein_uid", "model_acc", "env_start", "library", "stage"),
        sort_key=("protein_uid", "env_start", "model_acc"),
        partition_keys=("library", "stage"),
        enums={"cutoff": ("gathering", "trusted", "noise", "none")},
        doc=(
            "data_layer_and_tests.md §1 (L1 kind=domains), SEARCH_PLAN SD-04/SD-05. One row per included "
            "domain; `library` path value is <name>-<version> (e.g. pfam-38.2). Coordinates 0-based half-open "
            "on the protein (env/ali) and on the model (hmm). E-values at the fixed Z_prot = domZ = 1e8"
        ),
    ),
    TableSpec(
        name="cluster_member",
        location="l2/cluster_member",
        layer="L2",
        schema=pa.schema([
            _f("member_uid", UID),
            _f("level", pa.string()),
            _f("cluster_uid", UID),
            _f("rep_uid", UID),
            _f("cluster_params_hash", pa.string()),
        ]),
        primary_key=("member_uid", "level", "release"),
        sort_key=("level", "cluster_uid", "member_uid"),
        partition_keys=("release",),
        enums={"level": ("nr90", "nr50", "subfam")},
        doc="data_layer_and_tests.md §1 (L2), §2 (cluster_uid); one row per protein per level, reps included",
    ),
    TableSpec(
        name="cluster_member_aln",
        location="l2/cluster_member_aln",
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
        ]),
        primary_key=("member_uid", "level", "release"),
        sort_key=("level", "rep_uid", "member_uid"),
        partition_keys=("release",),
        enums={"level": ("nr90", "nr50", "subfam")},
        doc=(
            "data_layer_and_tests.md §2B.1; representatives carry a self row. q* are member coordinates, "
            "t* representative coordinates, 0-based half-open; pident is the identity fraction (0-1) over "
            "the alignment; cigar is the MMseqs2 CIGAR as ASCII bytes"
        ),
    ),
    TableSpec(
        name="id_history",
        location="l2/id_history.parquet",
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
            _f("submitted_at", pa.timestamp("ms", tz="UTC")),
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
        name="ledger",
        layer="ledger",
        location="ledger",
        schema=pa.schema([
            _f("object_key", pa.binary()),
            _f("object_kind", pa.string()),
            _f("n_in", pa.uint64()),
            _f("n_kept", pa.uint64()),
            _f("n_merged", pa.uint64()),
            _f("n_rejected", pa.uint64()),
            _f("reason_code", pa.string()),
            _f("tool_versions_json", pa.string(), nullable=True),  # C-08; appended 2026-10-30 (LY-09)
        ]),
        primary_key=("object_key", "reason_code", "stage", "step"),
        sort_key=("object_key", "reason_code"),
        partition_keys=("stage", "step"),
        enums={"object_kind": ("l0_object", "contig", "library", "level")},
        small_ok=True,
        doc=(
            "data_layer_and_tests.md §4.2. One row per (object, reason); n_rejected counts that reason only, so "
            "L-01 sums rows per object. object_key is the L0 object sha256, a contig_uid, or the UTF-8 name of a "
            "library or cluster level"
        ),
    ),
    TableSpec(
        name="contig_tapes",
        location="l3/tapes",
        layer="L3",
        schema=pa.schema([
            _f("contig_uid", UID),
            _f("tax_order_bucket", pa.uint32(), nullable=True),
            _f("genome_uid", UID, nullable=True),
            _f("contig_len", COORD),
            _f("is_circular", pa.bool_()),
            _f("n_slots", pa.uint32()),
            _f("tier", pa.uint8(), nullable=True),
            _f("assembly_quality", pa.string(), nullable=True),
            _f("feature_kinds", pa.list_(pa.uint8())),
            _f("strands", pa.list_(STRAND)),
            _f("start_bps", pa.list_(pa.uint32())),
            _f("end_bps", pa.list_(pa.uint32())),
            _f("delta_bp_prev_slots", pa.list_(pa.int32())),
            _f("junction_topology_prevs", pa.list_(pa.uint8())),
            _f("orfs_to_contig_starts", pa.list_(pa.uint16()), nullable=True),
            _f("orfs_to_contig_ends", pa.list_(pa.uint16()), nullable=True),
            _f("slot_flags", pa.list_(pa.uint8())),
            _f("protein_uids", pa.list_(UID)),
        ]),
        primary_key=("contig_uid", "gene_call_version"),
        sort_key=("tax_order_bucket", "contig_uid"),
        partition_keys=("gene_call_version", "stage"),
        list_columns=("feature_kinds", "strands", "start_bps", "end_bps", "delta_bp_prev_slots",
                      "junction_topology_prevs", "orfs_to_contig_starts", "orfs_to_contig_ends", "slot_flags",
                      "protein_uids"),
        doc=(
            "PLAN.md §4.2; biobox/l3/tapes.py. One row per contig, list index = slot_ordinal; coordinates "
            "only (append-only per stage). feature_kinds 1 = CDS; junction codes 1 co-directional, 2 divergent, "
            "3 convergent, 4 nested antisense, null for the first slot of a linear contig"
        ),
    ),
    TableSpec(
        name="token_dict",
        location="l3/token_dict",
        layer="L3",
        schema=pa.schema([
            _f("token_id", pa.uint64()),
            _f("level", pa.string()),
            _f("source", pa.string()),
            _f("key", pa.string()),
            _f("n_slots", pa.uint64()),
        ]),
        primary_key=("token_id", "release"),
        sort_key=("token_id",),
        partition_keys=("release",),
        enums={"level": ("family", "nr50", "nr90", "subfamily")},
        doc="biobox/l3/tokens.py: uint64 token = level(3) | source(5) | local(56), release-local",
    ),
    TableSpec(
        name="tape_tokens",
        location="l3/tape_tokens",
        layer="L3",
        schema=pa.schema([
            _f("contig_uid", UID),
            _f("family_tokens", pa.list_(pa.uint64())),
            _f("nr50_tokens", pa.list_(pa.uint64())),
            _f("nr90_tokens", pa.list_(pa.uint64())),
        ]),
        primary_key=("contig_uid", "gene_call_version", "release"),
        sort_key=("contig_uid",),
        partition_keys=("gene_call_version", "release"),
        list_columns=("family_tokens", "nr50_tokens", "nr90_tokens"),
        doc="data_layer_and_tests.md §1 L3 (`tape_tokens_<release>`): token columns, recomputed per release",
    ),
    TableSpec(
        name="slot_pairs",
        location="l3/slot_pairs",
        layer="L3",
        schema=pa.schema([
            _f("anchor_token", pa.uint64()),
            _f("partner_token", pa.uint64()),
            _f("contig_uid", UID),
            _f("anchor_ordinal", pa.uint32()),
            _f("partner_ordinal", pa.uint32()),
            _f("orf_offset", pa.int16()),
            _f("junction_topology", pa.uint8()),
            _f("delta_bp", pa.int32()),
            _f("flags", pa.uint8()),
            _f("tier", pa.uint8(), nullable=True),
            _f("tax_order", pa.uint32(), nullable=True),
        ]),
        primary_key=("anchor_token", "partner_token", "contig_uid", "anchor_ordinal", "partner_ordinal", "release"),
        sort_key=("anchor_token", "partner_token", "contig_uid", "anchor_ordinal", "partner_ordinal"),
        partition_keys=("release", "anchor_bucket"),
        enums={"junction_topology": (0, 1, 2, 3, 4, 5)},
        doc=(
            "PLAN.md §4.2; biobox/l3/pairs.py. Directed pairs within +-10 slots at one token level; flags bit 0 "
            "anchor window cut by a linear end, bit 1 pair crosses the origin; anchor_bucket = hex of the "
            "anchor token's low byte at the release's width"
        ),
    ),
    TableSpec(
        name="test_run",
        location="tests/test_run",
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


@dataclass(frozen=True)
class ArtifactSpec:
    """Non-Parquet files of a release, listed in the manifest like tables (LY-10).

    Files live at `<location>/<k>=<v>/.../<name>` in partition_keys order; `<name>` matches
    `pattern`, and each such file has one companion per suffix in `companions`.
    """

    name: str
    layer: str
    location: str
    partition_keys: tuple[str, ...]
    pattern: str
    companions: tuple[str, ...] = ()
    doc: str = ""


ARTIFACTS: dict[str, ArtifactSpec] = {
    "seq_shards": ArtifactSpec(
        name="seq_shards",
        layer="L0-seq",
        location="l0seq/seq_shards",
        partition_keys=("stage",),
        pattern=r"shard-\d{4}\.fa\.gz",
        companions=(".fai", ".gzi"),
        doc="bgzip FASTA of every stored byte sequence, records named by raw_digest hex and sorted "
            "by it (data_layer_and_tests.md §2A.1; DEC-17)",
    ),
}


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


def schema_lock() -> dict[str, list[list]]:
    """The registry as plain data, for tests/registry/schema_lock.json (LY-09)."""
    return {name: [[f.name, str(f.type), f.nullable] for f in spec.schema] for name, spec in sorted(TABLES.items())}


def evolution_problems(locked: dict[str, list[list]], current: dict[str, list[list]] | None = None) -> list[str]:
    """LY-09: a locked column may never be removed, renamed, retyped or reordered; new columns are
    appended and nullable, so files and catalogs from before the change stay readable."""
    current = current if current is not None else schema_lock()
    problems = []
    for name, fields in locked.items():
        if name not in current:
            problems.append(f"{name}: table removed")
            continue
        now = current[name]
        if now[: len(fields)] != fields:
            problems.append(f"{name}: a locked column was removed, renamed, retyped or reordered")
        for f in now[len(fields):]:
            if not f[2]:
                problems.append(f"{name}.{f[0]}: an added column must be nullable")
    return problems
