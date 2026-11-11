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

import re
from dataclasses import dataclass, field

import pyarrow as pa

from biobox.ids import KIND

UID = pa.binary(16)  # BLAKE3-128 content UIDs
SHA256 = pa.binary(32)  # raw_digest and object checksums
STRAND = pa.int8()  # +1 / -1
COORD = pa.int64()  # 0-based half-open contig coordinates

STRAND_VALUES = (1, -1)
# the domains of assembly quality and of a CDS's edit type, as the pre-phase-B tables had them
# (R0 checks the tables that carry them keep these values)
_ASSEMBLY_QUALITY = ("complete", "chromosome", "scaffold", "contig", "mag_hq", "mag_mq", "vOTU", "plasmid")
_EDIT_TYPES = ("none", "frameshift", "selenocysteine", "pyrrolysine", "intron")


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
    # Other partition key tuples a path under `location` may use instead of `partition_keys`; the first
    # key of each differs, so a path still maps to one key set (the ledger, card B-28).
    alt_partition_keys: tuple[tuple[str, ...], ...] = ()
    # Partition keys that only split one stage's or release's rows by key prefix (a width a builder
    # picks); they never separate LY-04 size groups (card C-09).
    bucket_keys: tuple[str, ...] = ()


def _f(name: str, typ: pa.DataType, nullable: bool = False) -> pa.Field:
    return pa.field(name, typ, nullable=nullable)


_SPECS = [
    TableSpec(
        name="l0_manifest",
        location="l0/l0_manifest",
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
        primary_key=("source_uid", "object_path", "stage"),
        sort_key=("source_uid", "object_path"),
        # the objects one stage ingested, written once by that stage's ingest (AGENTS §3.10, card B-22)
        partition_keys=("stage",),
        small_ok=True,
        doc="data_layer_and_tests.md §1 (L0); one directory per stage (card B-22)",
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
        location="l0seq/seq_locator",
        layer="L0-seq",
        schema=pa.schema([
            _f("raw_digest", SHA256),
            _f("shard", pa.string()),  # root-relative path of the .fa.gz; the record name is raw_digest hex
            _f("contig_len", COORD),
        ]),
        primary_key=("raw_digest",),
        sort_key=("raw_digest",),
        # one file per stage with the sequences first stored in it, never rewritten (AGENTS §3.10, card B-13)
        partition_keys=("stage",),
        doc="data_layer_and_tests.md §2A.1; DEC-17; planning/tasks/TASKS.md §2 (card B-13)",
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
        bucket_keys=("protein_prefix",),
        doc="data_layer_and_tests.md §2A.1",
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
            _f("left_feature_uid", UID, nullable=True),
            _f("right_feature_uid", UID, nullable=True),
        ]),
        primary_key=("contig_uid", "start", "end", "strand", "igr_version"),
        sort_key=("contig_uid", "start"),
        partition_keys=("igr_version", "stage"),
        enums={"strand": STRAND_VALUES},
        doc=(
            "data_layer_and_tests.md §2A.1; planning/tasks/TASKS.md §2. Coordinates only; the flanks are the "
            "feature uids of the tape version's producers, null at a linear contig end"
        ),
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
            _f("start", COORD),
            _f("end", COORD),
            _f("strand", STRAND),
            _f("dist_to_contig_end_bp", pa.uint32()),
            _f("feature_uid", UID),  # card B-05
        ]),
        primary_key=("feature_uid", "gene_call_version"),
        sort_key=("protein_uid",),
        partition_keys=("gene_call_version", "stage"),
        enums={"strand": STRAND_VALUES},
        doc=(
            "data_layer_and_tests.md §2B.1; planning/tasks/TASKS.md §2. Rows sorted by (protein_uid, contig_uid, "
            "start, end, strand). Genome, tier, taxonomy and assembly quality live in contig_dim (AGENTS §3.6; "
            "card B-17 removed them here with the ordinal column)"
        ),
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
    # card E-01 (biobox.l1.ncrna): one row per cmsearch hit at or above its model's GA; the hit's locus is a
    # `feature` row of kind ncrna under the same producer
    TableSpec(
        name="l1_ncrna",
        location="l1/ncrna",
        layer="L1",
        schema=pa.schema([
            _f("feature_uid", UID),
            _f("contig_uid", UID),
            _f("start", COORD),
            _f("end", COORD),
            _f("strand", STRAND),
            _f("model_acc", pa.string()),
            _f("model_name", pa.string()),
            _f("model_len", pa.uint32()),
            _f("mdl_start", pa.uint32()),
            _f("mdl_end", pa.uint32()),
            _f("mdl_kind", pa.string()),
            _f("bitscore", pa.float32()),
            _f("evalue_zref", pa.float64()),
            _f("ga_bits", pa.float32()),
            _f("z_mb", pa.float64()),
            _f("trunc_5p", pa.bool_()),
            _f("trunc_3p", pa.bool_()),
            _f("gc", pa.float32()),
            _f("contig_len", COORD),
        ]),
        primary_key=("feature_uid", "model_acc", "producer"),
        sort_key=("contig_uid", "start", "end", "strand", "model_acc"),
        partition_keys=("producer", "stage"),
        enums={"strand": STRAND_VALUES, "mdl_kind": ("cm", "hmm")},
        doc=(
            "card E-01; data_layer_and_tests.md §1 (L1 kind=ncrna); DEC-04; 00-spec.md §3.15. One row per Infernal "
            "cmsearch hit at or above its model's GA (`--cut_ga`), E-values at the fixed Z_nt = 1e5 Mb. "
            "feature_uid = ids.feature_uid(contig_uid, start, end, strand, KIND['ncrna']), the `feature` row of "
            "the same producer (`infernal-<v>-rfam-<v>`); two models that hit the same coordinates share it. "
            "start/end 0-based half-open on the contig, end > contig_len for a hit across the origin of a circular "
            "contig; mdl_start/mdl_end 0-based half-open on the model (of model_len consensus columns, CLEN); "
            "mdl_kind 'hmm' for a model searched by its profile HMM alone; trunc_5p/trunc_3p for a hit cut by a "
            "linear contig's end; gc the hit's GC fraction; contig_len the real target size (bases searched)"
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
    # --- card B-27: lookup tables sorted by the key they are looked up by, written beside their sources ---
    TableSpec(
        name="contig_label",
        location="l0seq/contig_label",
        layer="L0-seq",
        schema=pa.schema([
            _f("contig_uid", UID),
            _f("accession", pa.string()),
            _f("original_header", pa.string()),
            _f("n_records", pa.uint32()),
            _f("l0_object_sha256", SHA256),
            _f("record_ordinal", pa.uint64()),
        ]),
        primary_key=("contig_uid", "stage"),
        sort_key=("contig_uid",),
        partition_keys=("stage",),
        doc=(
            "card B-27: one row per contig of the stage's contig_occurrence (every source), written by ingest "
            "with it: the accession and header of the contig's first record by (l0_object_sha256, "
            "record_ordinal), that record's object and ordinal (so a release over several stages takes the "
            "first over all of them), and its number of records in the stage"
        ),
    ),
    TableSpec(
        name="member_cluster",
        location="l2/member_cluster",
        layer="L2",
        schema=pa.schema([
            _f("member_uid", UID),
            _f("level", pa.string()),
            _f("cluster_uid", UID),
            _f("rep_uid", UID),
        ]),
        primary_key=("member_uid", "level", "release"),
        sort_key=("level", "member_uid"),
        partition_keys=("release",),
        enums={"level": ("nr90", "nr50", "subfam")},
        doc=(
            "card B-27: cluster_member re-sorted by level and member, written by the clustering step with it: a "
            "protein's cluster and representative at each level, looked up by member_uid (one row group per "
            "level; a lookup at one level reads that level's rows only)"
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
            # appended 2026-10-30 with the runner (LY-09)
            _f("target", pa.string(), nullable=True),
            _f("target_size", pa.uint64(), nullable=True),
            _f("rounds", pa.uint16(), nullable=True),
            _f("final_round", pa.uint16(), nullable=True),
            _f("rounds_json", pa.string(), nullable=True),
            _f("stop_reason", pa.string(), nullable=True),
            _f("controls_json", pa.string(), nullable=True),
            _f("warnings_json", pa.string(), nullable=True),
            _f("limits_json", pa.string(), nullable=True),
            _f("error", pa.string(), nullable=True),
            # appended 2026-11-09 with W-02 (LY-09): the seed as resolved at submit; a file seed carries
            # its path as given, sha256 and what it holds (model name, accession, length; or seed names)
            _f("seed_json", pa.string(), nullable=True),
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
            # appended 2026-10-30 with the runner (LY-09): the sequence-level score that decides inclusion
            _f("seq_bitscore", pa.float32(), nullable=True),
            _f("seq_evalue_zref", pa.float64(), nullable=True),
        ]),
        primary_key=("run_id", "round", "target_kind", "target_uid", "env_start"),
        sort_key=("round", "target_kind", "target_uid"),
        enums={
            "target_kind": ("protein", "segment", "dark_orf", "igr", "contig_region", "spacer"),
            "score_kind": ("exact", "projected"),
        },
        small_ok=True,
        doc=(
            "data_layer_and_tests.md §2B.1; profile coordinates per SEARCH_PLAN SD-17. One row per included "
            "domain (env/ali on the target, hmm on the profile, 0-based half-open); round counts from 1"
        ),
    ),
    TableSpec(
        name="search_profile",
        layer="search",
        schema=pa.schema([
            _f("run_id", UID),
            _f("round", pa.uint16()),
            _f("kind", pa.string()),
            _f("profile_zst", pa.binary()),
            _f("msa_zst", pa.binary(), nullable=True),
            _f("n_seq", pa.uint32()),
            _f("n_eff", pa.float32(), nullable=True),
        ]),
        primary_key=("run_id", "round"),
        sort_key=("round",),
        enums={"kind": ("hmm", "cm", "pssm")},
        small_ok=True,
        doc="data_layer_and_tests.md §2B.1: the profile each round searched with (HMMER3 text, zstd) and its MSA "
            "(Stockholm, zstd)",
    ),
    TableSpec(
        name="neighborhoods",
        layer="search",
        schema=pa.schema([
            _f("anchor_protein_uid", UID),
            _f("contig_uid", UID),
            _f("anchor_ordinal", pa.uint32()),
            _f("ordinal", pa.uint32()),
            _f("offset", pa.int16()),
            _f("offset_rel", pa.int16()),
            _f("protein_uid", UID),
            _f("start", COORD),
            _f("end", COORD),
            _f("strand", STRAND),
            _f("delta_bp", pa.int32(), nullable=True),
            _f("junction_topology", pa.uint8(), nullable=True),
            _f("same_strand", pa.bool_()),
            _f("in_operon", pa.bool_(), nullable=True),
            _f("flags", pa.uint8()),
            _f("policy_hash", pa.string()),
            # card B-16: the slots by feature_uid (primary key, so non-null: appended before the first
            # release, AGENTS §3.8) and the tape version the ordinals belong to; the ordinals are plain
            # columns, never a key (TASKS.md §2, card B-17)
            _f("anchor_feature_uid", UID),
            _f("feature_uid", UID),
            _f("tape_version", pa.string(), nullable=True),
        ]),
        primary_key=("contig_uid", "anchor_feature_uid", "feature_uid", "policy_hash"),
        sort_key=("anchor_protein_uid", "contig_uid", "anchor_ordinal", "offset"),
        enums={"strand": STRAND_VALUES},
        small_ok=True,
        doc=(
            "SEARCH_PLAN SD-08; biobox/search/nbr.py. One row per slot in an anchor's window, the anchor "
            "itself at offset 0. offset is along the contig, offset_rel along the anchor's strand (positive = "
            "downstream); flags bit 0 window cut by a linear contig end, bit 1 crosses the origin. "
            "anchor_ordinal and ordinal are valid inside tape_version only (AGENTS §3.3)"
        ),
    ),
    TableSpec(
        name="search_region",
        layer="search",
        schema=pa.schema([
            _f("region_uid", UID),
            _f("target_kind", pa.string()),
            _f("contig_uid", UID),
            _f("start", COORD),
            _f("end", COORD),
            _f("strand", STRAND),
            _f("n_residues", pa.uint32()),
            _f("orientation", pa.uint8(), nullable=True),
            _f("gc", pa.float32(), nullable=True),
            _f("dist_to_contig_end_bp", pa.uint32(), nullable=True),
            _f("partial", pa.bool_()),
            _f("at_contig_end", pa.bool_()),
            _f("gene_overlap_frac", pa.float32(), nullable=True),
            _f("genetic_code", pa.uint8(), nullable=True),
            # card B-16: the flanks as features of the pack's tape version (pack_manifest.json
            # `tape_version`); they replace the ordinal flanks that card B-17 removed
            _f("left_feature_uid", UID, nullable=True),
            _f("right_feature_uid", UID, nullable=True),
        ]),
        primary_key=("region_uid",),
        sort_key=("region_uid",),
        enums={"target_kind": ("dark_orf", "igr", "contig_region", "spacer"), "strand": STRAND_VALUES},
        small_ok=True,
        doc=(
            "data_layer_and_tests.md §2B.2: `igr_meta.parquet` and `dark_orf_meta.parquet` of a search pack; "
            "region_uid = blake3(contig_uid, start, end, strand). n_residues counts nt for IGRs, aa for dark "
            "ORFs; orientation is the flanks' junction code (1 co-directional, 2 divergent, 3 convergent)"
        ),
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
        # a release's own steps (clustering) are not a stage's: ledger/release=<R>/step=<name>/ (card B-28)
        alt_partition_keys=(("release", "step"),),
        enums={"object_kind": ("l0_object", "contig", "library", "level")},
        small_ok=True,
        doc=(
            "data_layer_and_tests.md §4.2. One row per (object, reason); n_rejected counts that reason only, so "
            "L-01 sums rows per object. object_key is the L0 object sha256, a contig_uid, or the UTF-8 name of a "
            "library or cluster level. A stage's steps are partitioned by (stage, step); a release's own steps "
            "(clustering) by (release, step), so two releases over one stage never share a ledger path"
        ),
    ),
    TableSpec(
        name="contig_tapes",
        location="l3/tapes",
        layer="L3",
        schema=pa.schema([
            _f("contig_uid", UID),
            _f("contig_len", COORD),
            _f("is_circular", pa.bool_()),
            _f("n_slots", pa.uint32()),
            _f("feature_kinds", pa.list_(pa.uint8())),
            _f("strands", pa.list_(STRAND)),
            _f("start_bps", pa.list_(pa.uint32())),
            _f("end_bps", pa.list_(pa.uint32())),
            _f("delta_bp_prev_slots", pa.list_(pa.int32())),
            _f("junction_topology_prevs", pa.list_(pa.uint8())),
            _f("orfs_to_contig_starts", pa.list_(pa.uint16()), nullable=True),
            _f("orfs_to_contig_ends", pa.list_(pa.uint16()), nullable=True),
            _f("slot_flags", pa.list_(pa.uint8())),
            _f("protein_uids", pa.list_(UID)),  # entries nullable: a non-CDS slot has no protein
            _f("feature_uids", pa.list_(UID)),  # card B-07: one entry per slot, never null
        ]),
        primary_key=("contig_uid", "tape_version"),
        sort_key=("contig_uid",),
        partition_keys=("tape_version", "stage"),
        list_columns=("feature_kinds", "strands", "start_bps", "end_bps", "delta_bp_prev_slots",
                      "junction_topology_prevs", "orfs_to_contig_starts", "orfs_to_contig_ends", "slot_flags",
                      "protein_uids", "feature_uids"),
        doc=(
            "PLAN.md §4.2; planning/tasks/TASKS.md §2; biobox/l3/tapes.py. One row per contig and tape version; "
            "the slots are the features of the version's producers ordered by (start, end, strand, kind), and "
            "the list index (ordinal) is valid inside one tape_version only (AGENTS §3.3). feature_kinds from "
            "ids.KIND; protein_uids null for a non-CDS slot; junction codes 1 co-directional, 2 divergent, "
            "3 convergent, 4 nested antisense, null for the first slot of a linear contig. Genome, tier, "
            "taxonomy and assembly quality live in contig_dim (AGENTS §3.6; card B-17 removed them here)"
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
        primary_key=("contig_uid", "tape_version", "release"),
        sort_key=("contig_uid",),
        partition_keys=("tape_version", "release"),
        list_columns=("family_tokens", "nr50_tokens", "nr90_tokens"),
        doc="data_layer_and_tests.md §1 L3 (`tape_tokens_<release>`): token columns, recomputed per release",
    ),
    # --- phase B (card B-03): added beside the legacy tables they replaced; B-17 removed those ---
    TableSpec(
        name="genome_input",
        location="l0seq/genome_input",
        layer="L0-seq",
        schema=pa.schema([
            _f("contig_uid", UID),
            _f("genome_uid", UID),
            _f("is_circular", pa.bool_()),
            _f("kind", pa.string()),
            _f("genetic_code", pa.uint8()),
            _f("assembly_quality", pa.string()),
        ]),
        primary_key=("contig_uid", "stage"),
        sort_key=("contig_uid",),
        partition_keys=("stage",),
        enums={"kind": ("prokaryote", "virus", "plasmid"), "assembly_quality": _ASSEMBLY_QUALITY},
        doc="planning/tasks/TASKS.md §2; DEC-18. The gene caller's per-contig inputs, persisted per stage",
    ),
    TableSpec(
        name="feature",
        location="l1/feature",
        layer="L1",
        schema=pa.schema([
            _f("feature_uid", UID),
            _f("contig_uid", UID),
            _f("start", COORD),
            _f("end", COORD),
            _f("strand", STRAND),
            _f("kind", pa.uint8()),
        ]),
        primary_key=("feature_uid", "producer"),
        sort_key=("contig_uid", "start", "end", "strand", "kind"),
        partition_keys=("producer", "stage"),
        enums={"strand": STRAND_VALUES, "kind": tuple(sorted(KIND.values()))},
        doc=(
            "planning/tasks/TASKS.md §2; DEC-18. feature_uid = ids.feature_uid(contig_uid, start, end, strand, "
            "kind), kind from ids.KIND; 0-based half-open, end > contig_len for an origin-crossing feature"
        ),
    ),
    TableSpec(
        name="feature_protein",
        location="l1/feature_protein",
        layer="L1",
        schema=pa.schema([
            _f("contig_uid", UID),
            _f("feature_uid", UID),
            _f("start", COORD),
            _f("end", COORD),
            _f("strand", STRAND),
            _f("protein_uid", UID),
            _f("genetic_code", pa.uint8()),
            _f("edit_type", pa.string()),
        ]),
        primary_key=("feature_uid", "gene_call_version"),
        sort_key=("contig_uid", "start", "end"),
        partition_keys=("gene_call_version", "stage"),
        enums={"strand": STRAND_VALUES, "edit_type": _EDIT_TYPES},
        doc="planning/tasks/TASKS.md §2; DEC-18. The protein of each CDS feature (replaced the ordinal-keyed table)",
    ),
    TableSpec(
        name="contig_dim",
        location="dim/contig_dim",
        layer="dim",
        schema=pa.schema([
            _f("contig_uid", UID),
            _f("genome_uid", UID),
            _f("contig_len", COORD),
            _f("is_circular", pa.bool_()),
            _f("tier", pa.uint8(), nullable=True),
            _f("assembly_quality", pa.string(), nullable=True),
            _f("tax_order", pa.uint32(), nullable=True),
            _f("is_species_rep", pa.bool_()),
        ]),
        primary_key=("contig_uid", "release"),
        sort_key=("contig_uid",),
        partition_keys=("release",),
        enums={"assembly_quality": _ASSEMBLY_QUALITY},
        doc=(
            "planning/tasks/TASKS.md §2; DEC-20. Contig attributes that change without the sequence changing; "
            "the only place they live (AGENTS §3.6)"
        ),
    ),
    TableSpec(
        name="tape_version_def",
        location="l3/tape_version_def",
        layer="L3",
        schema=pa.schema([
            _f("producer", pa.string()),
        ]),
        primary_key=("tape_version", "producer"),
        sort_key=("producer",),
        partition_keys=("tape_version",),
        small_ok=True,
        doc=(
            "planning/tasks/TASKS.md §2; DEC-18. tape_version = 't-' + first 12 hex digits of sha256 over the "
            "sorted producer names joined by \\n; one row per (tape_version, producer), one directory per "
            "tape version so that a new version never rewrites a file (card B-07, AGENTS §3.10)"
        ),
    ),
    TableSpec(
        name="release_def",
        location="release_def",
        layer="release",
        schema=pa.schema([
            _f("position", pa.uint16()),
            _f("stage", pa.string()),
        ]),
        primary_key=("release", "stage"),
        sort_key=("position",),
        partition_keys=("release",),
        small_ok=True,
        doc=(
            "card B-25: the stages a release is built over, one row per stage, `position` its rank in the "
            "release's order (from 0). Written once per release (biobox.catalog.define_release); every "
            "release-level builder and the manifest read only these stages (AGENTS §3.9)"
        ),
    ),
    TableSpec(
        name="token_postings",
        location="l3/token_postings",
        layer="L3",
        schema=pa.schema([
            _f("token_id", pa.uint64()),
            _f("contig_uid", UID),
            _f("ordinal", pa.uint32()),
            _f("seg_idx", pa.uint16()),
            _f("feature_uid", UID),
            _f("start", COORD),
            _f("end", COORD),
            _f("strand", STRAND),
            _f("n_slots", pa.uint32()),
            _f("is_circular", pa.bool_()),
            _f("contig_len", COORD),
        ]),
        primary_key=("token_id", "contig_uid", "feature_uid", "seg_idx", "release", "tape_version"),
        sort_key=("token_id", "contig_uid", "ordinal", "seg_idx"),
        partition_keys=("release", "tape_version", "token_bucket"),
        bucket_keys=("token_bucket",),
        enums={"strand": STRAND_VALUES},
        doc=(
            "planning/tasks/TASKS.md §2; DEC-19. Replaces slot_pairs (removed by card B-11; pairs are computed "
            "at query time by the catalog's pair macros): one row per slot and token (seg_idx 0 for "
            "NR90/NR50, the domain's rank by env_start for family tokens). ordinal orders the rows within one "
            "tape_version but is not part of the key (AGENTS §3.3); token_bucket = hex of the token's low byte "
            "at the release's width"
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
    `pattern`, and each such file has one companion per suffix in `companions`. An artifact with
    no partition keys whose `location` ends in a name `pattern` matches is that one file
    (`single_file`), as a table whose location ends in `.parquet` is.

    `in_manifest` False marks the files that describe a release rather than belong to it: its
    manifest and its catalog (card C-08). They are claimed here so that `bio-box run-step` can
    publish them (LY-10), and are never listed in a release manifest nor reported as its strays.
    """

    name: str
    layer: str
    location: str
    partition_keys: tuple[str, ...]
    pattern: str
    companions: tuple[str, ...] = ()
    doc: str = ""
    in_manifest: bool = True
    # Small by nature: exempt from the LY-04 floor, as `TableSpec.small_ok` is (card C-09).
    small_ok: bool = False

    @property
    def single_file(self) -> bool:
        return not self.partition_keys and re.fullmatch(self.pattern, self.location.rsplit("/", 1)[-1]) is not None


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
    "store_params": ArtifactSpec(
        name="store_params",
        layer="L0-seq",
        location="l0seq/store_params.json",
        partition_keys=(),
        pattern=r"store_params\.json",
        doc="the sequence store's parameters (format, shard_bases, line_width, bgzf_level) as sorted-key JSON; "
            "written once, by the first stage's ingest, and checked, never rewritten, by every later one "
            "(card B-22)",
        small_ok=True,
    ),
    # A stage's ingest by batch processes (card B-30, `biobox.l0seq.store`): one record per pass-1 shard (its
    # objects and its scratch run's sha256) and the plan (the batches and the decision table's sha256), written
    # once each by the pass-1 and plan steps so that `bio-box run-step` can publish them. A build's own records,
    # like `domain_scan`: never in a release manifest, so a release is the same whether ingest was batched or not.
    "ingest_plan": ArtifactSpec(
        name="ingest_plan",
        layer="L0-seq",
        location="l0seq/ingest_plan",
        partition_keys=("stage",),
        pattern=r"(plan|pass1-\d{5})\.json",
        doc="a stage's ingest in batch processes: each pass-1 shard's record (pass1-<k>.json: its source objects "
            "by sha256, size and record count, and its scratch run's sha256) and the plan (plan.json: the "
            "objects in sha256 order, the batches of at most shard_bases new bases, one shard each), as "
            "sorted-key JSON",
        in_manifest=False,
        small_ok=True,
    ),
    # A stage's domain scan in shards (card B-24, `biobox.l1.domains`): the plan and one record per shard
    # (its scratch run's sha256 and counts), written once each by the plan and shard steps so that
    # `bio-box run-step` can publish them. A build's own records, like `runs/`: never in a release manifest,
    # so a release is the same whether its scans were sharded or not.
    "domain_scan": ArtifactSpec(
        name="domain_scan",
        layer="L1",
        location="l1/domain_scan",
        partition_keys=("library", "stage"),
        pattern=r"(plan|shard-\d{5})\.json",
        doc="a stage's shard plan for one domain library (plan.json: proteins by protein_uid cut by cumulative "
            "residues) and each scanned shard's record (shard-<k>.json), as sorted-key JSON",
        in_manifest=False,
        small_ok=True,
    ),
    # A release's own records (card C-08): one file per release, the release in the file name
    # (`release=<R>.json`, `release=<R>.duckdb`), so the artifact has no partition directories. Both
    # hold the data root's absolute path (the manifest's `root`; the catalog's `catalog_info.root` and
    # its views' file lists) and are built for the final root even when written into staging.
    "release_manifest": ArtifactSpec(
        name="release_manifest",
        layer="release",
        location="manifests",
        partition_keys=(),
        pattern=r"release=[^/=]+\.json",
        doc="a release's manifest, the system of record (biobox.catalog.build_manifest, write_manifest); "
            "written once per release",
        in_manifest=False,
        small_ok=True,
    ),
    "catalog": ArtifactSpec(
        name="catalog",
        layer="release",
        location="catalog",
        partition_keys=(),
        pattern=r"release=[^/=]+\.duckdb",
        doc="a release's generated read-only DuckDB catalog (biobox.catalog.generate_catalog): views over "
            "the manifest's explicit file lists",
        in_manifest=False,
        small_ok=True,
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
