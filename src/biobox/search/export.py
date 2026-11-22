"""Exports (operation O-EXPORT; SEARCH_PLAN SD-18; S-13; card W-05).

`--json`: NDJSON with versioned record schemas. Every record carries `"schema": "<kind>/<version>"`.
The rule, checked against `tests/registry/json_schema_lock.json`: keys may be added at the same
version; removing or renaming a key needs a version bump. Bytes are lowercase hex, timestamps ISO
8601 UTC. Lists are capped at query time (`top_n`, default 50), so one stored result serves any cap.

`search_hit/2` (card W-05): the nested occurrence objects are `contig_uid, ordinal, feature_uid, start,
end, strand, genome_uid, tier`. Card B-17 renamed `slot_ordinal` to `ordinal` and dropped `genome_uid`
and `tier` from the stored rows; version 2 records that change, and the two come back joined from the
run's release's `contig_dim` at export time (null when the release has no row for the contig).
Card C-21c (DEC-27) appends `n_genomes` and `genome_uids` to every occurrence (keys added at the same
version): the genomes that hold the contig, its `contig_dim` genome and every genome the release's
`genome_input` lists it in, in uid order (null when the release records none); `genome_uid` stays the
least of them. hits.tsv and each locus table's facts line end with the same two fields, on every
release; on a release whose contigs each lie in one genome they are 1 and that genome.

An ordinal is never exported alone (AGENTS §3.3, card B-16): a record or table row with an ordinal
carries the `feature_uid` beside it and the `tape_version` it belongs to. Cluster ids are cluster uids,
never token ids (AGENTS §3.4); genome and tier come from `contig_dim` (AGENTS §3.6).

Tables a person opens (card W-05): `hits_tsv` (one row per gene of a hit protein), `hits_faa` (one
record per hit protein), `locus_tables` (one table per anchor locus) and `protein_fasta` (every distinct
protein of a release). Scores are rounded to one decimal: bit scores as `%.1f`, E-values as `%.1e`.
They are DuckDB queries over explicit file lists (`catalog.ReleaseFiles`: the release manifest when
there is one, else the stages of the release's `release_def`; card W-12; AGENTS §3.9), each table filtered by key with `IN (SELECT ...)` before a join, and
the results stream out in batches.
"""

from __future__ import annotations

import bisect
import json
import os
import re
import shutil
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Iterable, Iterator

import pyarrow as pa

SCHEMAS: dict[str, dict] = {
    "search_hit": {"version": 2, "keys": [
        "schema", "run_id", "round", "target_kind", "target_uid", "bitscore", "evalue_zref", "target_size",
        "seq_bitscore", "seq_evalue_zref", "env_start", "env_end", "ali_start", "ali_end", "hmm_from", "hmm_to",
        "model_len", "is_partial", "first_inclusion_round", "score_kind", "curation_label", "curation_reason",
        "n_occurrences", "occurrences", "tape_version"]},
    "search_run": {"version": 1, "keys": [
        "schema", "run_id", "status", "error", "recipe_hash", "written_against", "ran_against", "pack_id", "target",
        "target_size", "rounds", "final_round", "stop_reason", "controls", "warnings", "submitted_at", "wall_s",
        "cpu_min_max_pid", "peak_rss_bytes", "z_prot", "z_nt_mb", "partial_policy", "max_output_gb"]},
    "neighborhood": {"version": 1, "keys": [
        "schema", "anchor_protein_uid", "contig_uid", "anchor_ordinal", "ordinal", "offset", "offset_rel",
        "protein_uid", "start", "end", "strand", "delta_bp", "junction_topology", "same_strand", "in_operon",
        "censored", "crosses_origin", "policy_hash", "anchor_feature_uid", "feature_uid", "tape_version"]},
    "pack_target": {"version": 1, "keys": ["schema", "release_id", "pack_id", "target", "kind", "n_seq",
                                           "n_residues", "shards"]},
    # card W-06: `bio-box nbr summarize`, one header record, then one record per key and level
    "nbr_summary_header": {"version": 1, "keys": [
        "schema", "release", "tape_version", "policy", "policy_hash", "levels", "top_n", "n_anchor_proteins_given",
        "n_anchor_proteins", "anchors_dropped_by_edge", "n_anchor_loci", "n_contigs", "n_anchor_genomes",
        "genome_basis", "n_censored_loci", "n_neighbors", "n_keys"]},
    "nbr_summary": {"version": 1, "keys": [
        "schema", "release", "tape_version", "policy_hash", "level", "key", "unannotated", "library", "model_acc",
        "model_name", "cluster_uid", "rep_uid", "rep_length", "n_loci", "share_of_loci", "n_neighbors", "n_genomes",
        "genome_basis", "offsets", "modal_offset", "same_strand_share"]},
}
TOP_N = 50
OCCURRENCES_SHOWN = 5


def _plain(v):
    if isinstance(v, bytes):
        return v.hex()
    if isinstance(v, datetime):
        return v.isoformat()
    if isinstance(v, list):
        return [_plain(x) for x in v]
    if isinstance(v, dict):
        return {k: _plain(x) for k, x in v.items()}
    return v


def record(kind: str, data: dict) -> dict:
    spec = SCHEMAS[kind]
    out = {"schema": f"{kind}/{spec['version']}"}
    for k in spec["keys"][1:]:
        out[k] = _plain(data.get(k))
    return out


def hit_records(hits: list[dict], occurrences: dict[bytes, list[dict]] | None = None,
                top_n: int | None = TOP_N, members: bool = False, tape_version: str | None = None,
                dims: dict[bytes, dict] | None = None) -> Iterator[dict]:
    """Best domain first per target, targets by descending sequence score. Each occurrence's `ordinal`
    is valid inside `tape_version` (nbr.occurrences looks it up there) and travels with its `feature_uid`.
    `genome_uid` and `tier` are not stored on an occurrence (AGENTS §3.6): they come from `dims`
    (`contig_dims`, the release's `contig_dim` rows by contig), null for a contig it lacks; `n_genomes` and
    `genome_uids` (card C-21c) every genome that holds the contig, null when the release records none."""
    rows = sorted(hits, key=lambda r: (-(r["seq_bitscore"] or 0), r["target_uid"], -r["bitscore"]))
    dims = dims or {}
    for r in rows[:top_n] if top_n else rows:
        occ = (occurrences or {}).get(r["target_uid"], [])
        shown = occ if members else occ[:OCCURRENCES_SHOWN]
        yield record("search_hit", {**r, "n_occurrences": len(occ), "tape_version": tape_version, "occurrences": [
            {"contig_uid": o["contig_uid"], "ordinal": o["ordinal"], "feature_uid": o["feature_uid"],
             "start": o["start"], "end": o["end"], "strand": o["strand"],
             "genome_uid": dims.get(o["contig_uid"], {}).get("genome_uid"),
             "tier": dims.get(o["contig_uid"], {}).get("tier"),
             "n_genomes": len(dims.get(o["contig_uid"], {}).get("genome_uids") or []) or None,
             "genome_uids": dims.get(o["contig_uid"], {}).get("genome_uids") or None} for o in shown]})


def run_record(row: dict) -> dict:
    return record("search_run", {**row, "controls": json.loads(row["controls_json"]) if row.get("controls_json") else None,
                                 "warnings": json.loads(row["warnings_json"]) if row.get("warnings_json") else []})


def nbr_records(rows: list[dict], top_n: int | None = None) -> Iterator[dict]:
    for r in rows[:top_n] if top_n else rows:
        yield record("neighborhood", {**r, "censored": bool(r["flags"] & 1), "crosses_origin": bool(r["flags"] & 2)})


def nbr_summary_records(summary: dict) -> Iterator[dict]:
    """`nbr.summarize` as records: the header (run-level facts), then each level's rows in order."""
    h = summary["header"]
    yield record("nbr_summary_header", h)
    for r in summary["rows"]:
        yield record("nbr_summary", {"release": h.get("release"), "tape_version": h.get("tape_version"),
                                     "policy_hash": h.get("policy_hash"), **r})


# At the cluster levels `key` is the cluster uid: the TSV does not repeat it as `cluster_uid` (card W-24); the
# `nbr_summary/1` records keep that key (a key leaves a record only with a version bump, S-13)
SUMMARY_TSV_COLUMNS = {
    "pfam": ["level", "key", "unannotated", "library", "model_acc", "model_name"],
    "nr": ["level", "key", "unannotated", "rep_uid", "rep_length", "model_name"],
}
SUMMARY_TSV_STATS = ["n_loci", "share_of_loci", "n_neighbors", "n_genomes", "genome_basis", "offsets", "modal_offset",
                     "same_strand_share"]


def _tsv_cell(col: str, v) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if col == "offsets":
        return ",".join(f"{o}:{n}" for o, n in v.items())
    return str(v)


def nbr_summary_tsv(summary: dict) -> Iterator[str]:
    """`--format tsv`: the header's facts as `# key<TAB>json` lines, then one table per level (a column line and
    its rows; offsets as `offset:count,...`), tables separated by a blank line. On a release with CRISPR arrays
    (`summary["arrays"]`, card E-14) one more table after them: `arrays.SUMMARY_COLUMNS`, one row per producer."""
    for k, v in record("nbr_summary_header", summary["header"]).items():
        yield f"# {k}\t{json.dumps(v, separators=(',', ':'))}"
    for level in summary["header"]["levels"]:
        level_rows = [r for r in summary["rows"] if r["level"] == level]
        base_cols = list(SUMMARY_TSV_COLUMNS["pfam" if level == "pfam" else "nr"])
        if level == "pfam" and any("source" in r for r in level_rows):
            base_cols.append("source")
        cols = base_cols + SUMMARY_TSV_STATS
        yield ""
        yield "\t".join(cols)
        for r in level_rows:
            yield "\t".join(_tsv_cell(c, r.get(c)) for c in cols)
    if summary.get("arrays") is not None:  # card E-14: only on a release that holds CRISPR arrays
        from biobox.search import arrays

        yield ""
        yield "\t".join(arrays.SUMMARY_COLUMNS)
        for r in summary["arrays"]:
            yield "\t".join(_tsv_cell(c, r.get(c)) for c in arrays.SUMMARY_COLUMNS)


README_TEXT = """\
Files written by `bio-box search hmm` (README.txt; card W-20)

recipe.yaml     the recipe the run ran.
hits.tsv        one row per gene of a hit protein.
hits.faa        one record per hit protein: >protein_uid, then length, scores, n_genes and one example gene.
neighbors.tsv   which protein families and clusters sit next to the hits (one table per level).
loci/           one table per anchor locus (a hit gene and the genes around it), locus_<accession>_<feature_uid>.tsv.
RUN             the run's id and its release, `key<TAB>value` per line (`bio-box results archive --under` finds this
                directory by it).

Conventions, the same in every file
  Coordinates (start, end) are 0-based and half-open: the first base of a contig is 0 and `end` is the first base
  after the gene. A GenBank range 3171..5822 is start 3170, end 5822.
  Strand is +1 or -1 (never 1, (+1) or "plus"). A gene on a circular contig that crosses the origin has end greater
  than the contig length.
  Uids (protein_uid, contig_uid, feature_uid, genome_uid, nr90, nr50, nr50_rep, rep_uid, and key at the cluster
  levels of neighbors.tsv) are 32 hex digits. nr90 and nr50 are cluster uids of the release at 90% and 50% identity, not numbers to compare across
  releases.
  An ordinal (position of a gene in its contig's tape) means something only inside its tape_version and is always
  written beside the feature_uid of the gene; use feature_uid to name a gene.
  Scores: seq_bitscore and hit_seq_bitscore in bits, one decimal; seq_evalue_zref at the reference database size,
  scientific notation with two significant digits. An empty cell means no value.

hits.tsv columns
  accession, header   the first source record of the contig (by source file hash, then record number).
  n_records           how many source records the contig has: records of the input FASTA that hold this identical
                      sequence and were stored as one contig. It is not a count of hits.
  contig_len          length of the contig's default variant, in bases.
  start, end, strand  the gene (see conventions).
  length              the protein's length in amino acids.
  seq_bitscore        the protein's best sequence score in the run; seq_evalue_zref its E-value.
  pfam_architecture   names of the domain library's hits on the protein, in order along it, a run of one name
                      written once, joined with +; - when scanned with no hit (domain_state is none), ? when not
                      scanned (including a release with no domain library).
  domain_state        how the protein's domain state is known in the release: measured, inferred, none (in an
                      annotated stage, no hit) or not scanned.
  protein_uid         the protein. nr90, nr50: its clusters; nr50_rep: the representative protein of its nr50 cluster.
  genome_uid          the genome of the contig, from the release's contig_dim: when the contig lies in several
                      genomes (an identical plasmid in two isolates, say), the least of their uids.
  tier                present only when the release's contig_dim has tiers (the tiering step does not exist yet, so
                      it is left out of every table built so far).
  contig_uid, feature_uid, ordinal, tape_version   the gene's identity; see conventions.
  dist_to_contig_end_bp   from the gene caller: 0 for a gene it reports as partial at a linear contig end (it runs
                      off the end or stops within a codon of it, whatever the distance in bases); the contig's
                      length for any gene of a circular contig, which has no end; otherwise the distance in bases
                      from the gene to the nearer end.
  stage               the build stage whose gene calls hold the gene (a release is built over one or more stages,
                      for example a phage stage and a prokaryote stage), so hits split by source with no SQL.
  kind                the kind of the genome the contig belongs to in that stage, as its inputs declared it:
                      virus, prokaryote or plasmid; empty when the stage recorded no genome for the contig.
  n_genomes           how many genomes of the release hold this contig (1 for most; more when the identical
                      sequence is a contig of several genomes); empty when the release records no genome for it.
  genome_uids         those genomes, comma-separated, in uid order. A gene is still one row, however many genomes
                      hold its contig: count genomes from this column, not from rows.

loci/*.tsv: the first line is the contig's header, the second `key=value` facts of the locus (anchor_feature_uid,
anchor_ordinal, anchor_protein_uid, censored = a linear contig end cut the window, ..., n_genomes and genome_uids as
in hits.tsv), the third the column names.
  offset              the slot's place relative to the anchor in coordinate order (lower coordinates negative).
  offset_rel          the same relative to the anchor's strand: positive is downstream of the anchor along the
                      anchor's strand, so on a -1 anchor the signs are the opposite of offset.
  ordinal, feature_uid, start, end, strand, protein_uid   the slot's gene.
  length_aa, pfam_architecture, domain_state, nr90, nr50   as in hits.tsv, empty for a slot that is not a protein.
  hit_seq_bitscore    the slot protein's score in the run, empty when it is not a hit.
  delta_bp            the gap in bases to the anchor; empty at the anchor.
  crosses_origin      true when the slot lies across the origin of a circular contig.

neighbors.tsv: a block of header lines (# key<TAB>json), then one table per level (pfam, nr50 or nr90), blank-line
separated. The offsets are relative to the anchor's strand (as offset_rel above), not coordinate order.
  level, key          the level, and what is counted: at the pfam level key is the family (library:model_acc); at
                      the cluster levels key is the cluster uid (written once: there is no separate cluster_uid
                      column; the --json records of `bio-box nbr summarize` carry it again as cluster_uid).
  unannotated         true for the rows of neighbours that have no family or cluster at that level (at pfam, key is
                      empty for scanned neighbours with no hit or non-protein slots, and ? with model_name "not
                      scanned" when the window holds proteins not scanned in the release).
  library, model_acc, model_name   the family (pfam level). At the cluster levels model_name is the domain name most
                      members of the cluster carry, empty when none has one, and so empty on a release without a
                      domain library; the cluster is then named only by its uid and its representative.
  rep_uid, rep_length the cluster's representative protein and that protein's length in amino acids (cluster levels).
  n_loci              anchor loci with such a neighbour; share_of_loci is n_loci over all anchor loci.
  n_neighbors         neighbour genes counted (a gene in the windows of two anchors counts twice).
  n_genomes           genomes with such a neighbour; genome_basis says what was counted: genome_uid from contig_dim,
                      else contig_uid. A neighbour on a contig of several genomes counts in each of them here, and
                      once in n_loci and n_neighbors.
  offsets             offset:count pairs, comma-separated, offset relative to the anchor's strand.
  modal_offset        the commonest of those offsets.
  same_strand_share   the share of those neighbours on the anchor's strand.
"""

# card E-14: appended to README.txt only when the release holds CRISPR arrays, so a README of a release without them
# is the same bytes as before
README_ARRAYS = """
CRISPR arrays (this release holds them; card E-14)
  Arrays are not genes and not slots: they are found by coordinates in each locus's window in bases. The window
  of a gene window runs from the least start to the greatest end of the locus's genes (a gene across a circular
  origin moved to the anchor's side), that of a window in bases from the anchor's start minus the window to its end
  plus the window. An array is listed when it shares at least one base with that span.
  loci/*.tsv end with `# crispr_arrays<TAB>window_start=..<TAB>window_end=..` (the span, in the anchor's frame), a
  column line and one row per array within the window (none: the two lines alone):
  kind                crispr_array.
  feature_uid, producer   the array (its stable key) and the MinCED setting that called it.
  start, end, strand  the array on the contig, as stored (strand +1: MinCED gives no orientation).
  n_repeats, n_spacers, repeat   its repeats and spacers, and its commonest repeat sequence.
  delta_bp            the gap in bases to the anchor; negative when the array overlaps it.
  side_rel            +1 downstream of the anchor along its strand, -1 upstream, 0 overlapping it.
  partial             true when part of the array lies outside the window.
  crosses_origin      true when the array is reached across the origin of a circular contig.
  neighbors.tsv ends with one more table, one row per array producer: kind, producer, n_loci, n_loci_with_array,
  share_of_loci, n_loci_partial_only (loci whose only arrays in the window are partial), n_arrays.
"""


def ndjson(records: Iterable[dict]) -> Iterator[str]:
    for r in records:
        yield json.dumps(r, sort_keys=False, separators=(",", ":"))


def lock_path(repo: Path) -> Path:
    return repo / "tests/registry/json_schema_lock.json"


def lock_problems(locked: dict, current: dict | None = None) -> list[str]:
    """S-13: a key may be removed or renamed only with a version bump."""
    current = SCHEMAS if current is None else current
    problems = []
    for kind, spec in locked.items():
        if kind not in current:
            problems.append(f"{kind}: schema removed")
            continue
        now = current[kind]
        if now["version"] < spec["version"]:
            problems.append(f"{kind}: version went down")
        elif now["version"] == spec["version"]:
            gone = [k for k in spec["keys"] if k not in now["keys"]]
            if gone:
                problems.append(f"{kind}/{spec['version']}: keys {gone} removed or renamed without a version bump")
    return problems


# --- tables a person opens (card W-05) -------------------------------------------------------------------

HITS_TSV_COLUMNS = ["accession", "header", "contig_len", "start", "end", "strand", "length", "seq_bitscore",
                    "seq_evalue_zref", "pfam_architecture", "domain_state", "protein_uid", "nr90", "nr50", "nr50_rep",
                    "genome_uid", "tier", "contig_uid", "feature_uid", "ordinal", "tape_version",
                    "dist_to_contig_end_bp", "n_records", "stage", "kind",  # card W-23 appended stage and kind
                    "n_genomes", "genome_uids"]  # card C-21c: every genome that holds the contig
LOCUS_TSV_COLUMNS = ["offset", "offset_rel", "ordinal", "feature_uid", "start", "end", "strand", "protein_uid",
                     "length_aa", "pfam_architecture", "domain_state", "nr90", "nr50", "hit_seq_bitscore", "delta_bp",
                     "crosses_origin"]
_BITS = {"seq_bitscore", "hit_seq_bitscore"}
_EVALUE = {"seq_evalue_zref"}
_STRAND = {"strand"}  # card W-20: `+1` or `-1` in every table a person opens (AGENTS §3.1), as `bio-box find` writes it
EXPORT_THREADS = 4  # DuckDB threads of an export without --threads (DuckDB's own default is every core of the host)
EXPORT_MEMORY_LIMIT = "2GB"  # larger sorts and aggregations spill into a scratch directory outside the data root
_BATCH_ROWS = 8192
_SAFE = re.compile(r"[^A-Za-z0-9._-]")

# The store relations the exports read: registry table, the partitions a file must match (values from
# open_store's arguments; a null value matches no file), the columns, and partition values kept as columns.
_RELATIONS = {
    "occ": ("protein_occurrence", ("gene_call_version",),
            ("protein_uid", "contig_uid", "feature_uid", "start", "end", "strand", "dist_to_contig_end_bp"), ("stage",)),
    # card W-23: the kind of the genome a contig belongs to in a stage (virus, prokaryote, plasmid), sorted by contig
    "gin": ("genome_input", (), ("contig_uid", "genome_uid", "kind"), ("stage",)),  # card C-21c: and its genomes
    "prot": ("protein_seq", ("gene_call_version",), ("protein_uid", "length", "seq"), ("stage",)),
    "dom": ("l1_domains", (), ("protein_uid", "model_acc", "model_name", "env_start"), ("library",)),
    # by member and by contig: the lookup tables of card B-27, sorted by those keys (cluster_member and
    # contig_occurrence are sorted otherwise, so a lookup by member or contig read every row group of them)
    "clus": ("member_cluster", ("release",), ("member_uid", "level", "cluster_uid", "rep_uid"), ()),
    "dim": ("contig_dim", ("release",), ("contig_uid", "genome_uid", "tier"), ()),
    "lab": ("contig_label", (), ("contig_uid", "accession", "original_header", "n_records", "l0_object_sha256",
                                 "record_ordinal"), ()),
    "cidx": ("contig_index", (), ("contig_uid", "contig_len", "is_default_variant"), ()),
    "tape": ("contig_tapes", ("tape_version",), ("contig_uid", "feature_uids"), ()),
}
# Where a release built before card B-27 (no `member_cluster`, no `contig_label`) is read instead: the table the
# lookup is sorted from and its columns. Slower (it cannot prune by member or contig), never empty (card W-17). A
# release whose stages are only partly labelled (card W-36) reads the label files of the stages that have them and
# `contig_occurrence` for the others (`_unlabelled`), as one relation.
_PRE_B27 = {
    "clus": ("cluster_member", ("member_uid", "level", "cluster_uid", "rep_uid")),
    "lab": ("contig_occurrence", ("contig_uid", "accession", "original_header", "l0_object_sha256", "record_ordinal")),
}


def _unlabelled(root: Path, files, labels: list[tuple[Path, dict]]) -> tuple[list[str], list[tuple[Path, dict]]] | None:
    """Card W-36: for a release some of whose stages have `contig_label` files (`labels`) and others none (ingested
    before card B-27, as S1, P1 and P2 of `rM.2`), (those other stages in the release's order, their
    `contig_occurrence` files), from which their labels are read as card W-18 reads a release with no label file at
    all; None when every stage has label files or none has (`labels` empty: the W-18 path)."""
    from biobox.search.nbr import _store_files

    if not labels:
        return None
    have = {str(p.get("stage")) for _, p in labels}
    missing = [s for s in files.stages if s not in have]
    if not missing:
        return None
    occ = [(p, parts) for p, parts in _store_files(root, "contig_occurrence", {}, files)
           if str(parts.get("stage")) in missing]
    return missing, occ


def _partial_label_note(release: str | None, stages: list[str]) -> None:
    """One stderr line (card W-36, as `nbr.slow_path_note`, once per command inside `nbr.noted_once`): the stages
    without label files are read from their `contig_occurrence`, which cannot prune a lookup by contig."""
    from biobox.search import nbr

    key = (release, "contig_label", ",".join(stages))
    if nbr._NOTED is not None:
        if key in nbr._NOTED:
            return
        nbr._NOTED.add(key)
    print(f"note: release {release}: stages {','.join(stages)} have no contig_label files (ingested before card "
          f"B-27); their contigs' labels are read from their contig_occurrence instead (the slow path for those "
          f"stages: each of their row groups is scanned)", file=sys.stderr)


def cell(col: str, v) -> str:
    """One TSV cell: empty for null; bit scores `%.1f`, E-values `%.1e`; bytes as lowercase hex; booleans
    `true`/`false`; tabs and line breaks inside text become spaces."""
    if v is None:
        return ""
    if col in _BITS:
        return f"{v:.1f}"
    if col in _EVALUE:
        return f"{v:.1e}"
    if col in _STRAND:
        return f"{v:+d}"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, bytes):
        return v.hex()
    return str(v).replace("\t", " ").replace("\n", " ").replace("\r", " ")


def _tsv_line(cols: list[str], row: dict) -> str:
    return "\t".join(cell(c, row[c]) for c in cols)


def _q(s) -> str:
    return "'" + str(s).replace("'", "''") + "'"


def release_manifest(root: Path, release: str | None) -> dict | None:
    """`manifests/release=<release>.json` under root, or None when there is none."""
    path = Path(root) / f"manifests/release={release}.json"
    if release is None or not path.exists():
        return None
    from biobox import catalog

    return catalog.read_manifest(path)


def release_files(root: Path, release: str | None, manifest=None):
    """The release's files for a reader (card W-12): a `catalog.ReleaseFiles` as given, else
    `ReleaseFiles.of(root, release, manifest)` (the manifest, else the `release_def` stages). ValueError, naming
    the release and `bio-box release define`, when the release has no manifest and no `release_def`."""
    from biobox import catalog

    if isinstance(manifest, catalog.ReleaseFiles):
        return manifest
    if release is None:
        raise ValueError("no release: a hit or locus export reads one release's files; name it (bio-box release define "
                         "writes its stages)")
    try:
        return catalog.ReleaseFiles.of(Path(root), release, manifest)
    except FileNotFoundError:
        raise ValueError(f"release {release} is not defined under {root}: it has no manifest and no release_def; "
                         f"define it with: bio-box release define --root {root} --release {release} "
                         f"--stages S0,...") from None


class _Store:
    def __init__(self, con, libraries: list[str], annotated_stages: dict[str, tuple[str, ...]] | None = None):
        self.con, self.libraries = con, libraries
        self.annotated_stages = annotated_stages or {}

    def library(self, wanted: str | None) -> str | None:
        """The domain library of the architecture column: `wanted`, else the only one, else the only one whose
        name starts with `pfam`; None when there is none (every architecture is then `?`, `not scanned`)."""
        libs = self.libraries
        if wanted is not None:
            if wanted not in libs:
                raise ValueError(f"domain library {wanted!r} is not in the store (found: {', '.join(libs) or 'none'})")
            return wanted
        if len(libs) <= 1:
            return libs[0] if libs else None
        pfam = [x for x in libs if x.lower().startswith("pfam")]
        if len(pfam) == 1:
            return pfam[0]
        raise ValueError(f"several domain libraries ({', '.join(libs)}): name one with --library")

    def stages_of(self, library: str | None) -> tuple[str, ...]:
        """The release's stages that hold an `l1_domains` file for `library` (card C-38b)."""
        return () if library is None else self.annotated_stages.get(library, ())


def _dom_view(con, files: list[tuple[Path, dict]], rdom_files: list[tuple[Path, dict]] = ()) -> None:
    """Temp view `dom` with `(protein_uid, model_acc, model_name, env_start, inferred_from, library, source)`
    over explicit `l1_domains` (`source = 'stage'`) and merged `round_domains` (`source = 'round'`), filling
    `inferred_from` as NULL for pre-B-47c and round files (LY-09, card C-38g)."""
    import pyarrow.parquet as pq

    if not files and not rdom_files:
        con.execute('CREATE TEMP VIEW dom AS SELECT CAST(NULL AS BLOB) AS "protein_uid", '
                    'CAST(NULL AS VARCHAR) AS "model_acc", CAST(NULL AS VARCHAR) AS "model_name", '
                    'CAST(NULL AS UINTEGER) AS "env_start", CAST(NULL AS BLOB) AS "inferred_from", '
                    'CAST(NULL AS VARCHAR) AS "library", CAST(NULL AS VARCHAR) AS "source" WHERE false')
        return
    parts_sql = []
    for path, parts in sorted(files, key=lambda x: (str(x[1].get("library", "")), str(x[0]))):
        has_inf = "inferred_from" in pq.read_schema(path).names
        inf_col = '"inferred_from"' if has_inf else 'CAST(NULL AS BLOB) AS "inferred_from"'
        parts_sql.append(
            f'SELECT "protein_uid", "model_acc", "model_name", "env_start", {inf_col}, '
            f'{_q(str(parts.get("library", "")))} AS "library", \'stage\' AS "source" '
            f'FROM read_parquet([{_q(str(path))}], hive_partitioning = false)'
        )
    for path, parts in sorted(rdom_files, key=lambda x: (str(x[1].get("library", "")), str(x[0]))):
        parts_sql.append(
            f'SELECT "protein_uid", "model_acc", "model_name", "env_start", CAST(NULL AS BLOB) AS "inferred_from", '
            f'{_q(str(parts.get("library", "")))} AS "library", \'round\' AS "source" '
            f'FROM read_parquet([{_q(str(path))}], hive_partitioning = false)'
        )
    con.execute(f"CREATE TEMP VIEW dom AS {' UNION ALL '.join(parts_sql)}")


def _view(con, name: str, spec_name: str, files: list[tuple[Path, dict]], cols: tuple[str, ...],
          keep: tuple[str, ...]) -> None:
    """A temp view over explicit files (one read_parquet per set of kept partition values, those values as
    literal columns); typed and empty without files."""
    from biobox.catalog import duck_type
    from biobox.registry.schemas import table

    types = {f.name: f.type for f in table(spec_name).schema}
    if not files:
        sel = [f'CAST(NULL AS {duck_type(types[c])}) AS "{c}"' for c in cols]
        sel += [f'CAST(NULL AS VARCHAR) AS "{k}"' for k in keep]
        con.execute(f"CREATE TEMP VIEW {name} AS SELECT {', '.join(sel)} WHERE false")
        return
    groups: dict[tuple, list[str]] = {}
    for path, parts in files:
        groups.setdefault(tuple(str(parts.get(k, "")) for k in keep), []).append(str(path))
    names = ", ".join(f'"{c}"' for c in cols)
    parts_sql = []
    for vals, paths in sorted(groups.items()):
        lits = "".join(f', {_q(v)} AS "{k}"' for k, v in zip(keep, vals))
        listing = ", ".join(_q(p) for p in sorted(paths))
        parts_sql.append(f"SELECT {names}{lits} FROM read_parquet([{listing}], hive_partitioning = false)")
    con.execute(f"CREATE TEMP VIEW {name} AS {' UNION ALL '.join(parts_sql)}")


class open_store:
    """A DuckDB connection with temp views over the store's relations (`_RELATIONS`) for one release, gene-call
    version and tape version; files from `release_files` (`manifest`, else the release's stages). Threads and memory
    are capped, and spills go to a private scratch directory under `tmp` (system temp by default). Keys taken from a
    temp table (`IN (SELECT ...)`) prune each relation to the row groups that hold them (card B-27): the
    connection's join filter lists as many keys as the largest relation has row groups
    (`catalog.set_join_key_filter`)."""

    def __init__(self, root: Path, *, release: str | None, gene_call_version: str, tape_version: str | None = None,
                 manifest: dict | None = None, threads: int | None = None, tmp=None,
                 relations: tuple[str, ...] = tuple(_RELATIONS)):
        self.args = (Path(root), release, gene_call_version, tape_version, manifest, threads, tmp, relations)

    def __enter__(self) -> _Store:
        import duckdb

        from biobox.catalog import row_groups_of, set_join_key_filter
        from biobox.search.nbr import _merged_round_files, _store_files, slow_path_note

        root, release, gcv, tv, manifest, threads, tmp, relations = self.args
        manifest = release_files(root, release, manifest)  # refuses an undefined release before anything is opened
        self.scratch = Path(tempfile.mkdtemp(prefix="biobox-export-", dir=tmp))
        self.con = duckdb.connect(config={"threads": threads or EXPORT_THREADS, "memory_limit": EXPORT_MEMORY_LIMIT,
                                          "preserve_insertion_order": False,
                                          "temp_directory": str(self.scratch / "duckdb")})
        values = {"gene_call_version": gcv, "release": release, "tape_version": tv}
        libraries: set[str] = set()
        annotated: dict[str, set[str]] = {}
        try:
            largest = 0
            for name in relations:
                spec_name, filt, cols, keep = _RELATIONS[name]
                want = {k: values[k] for k in filt}
                files = [] if None in want.values() else _store_files(root, spec_name, want, manifest)
                if name in _PRE_B27 and not files and None not in want.values():
                    # a release built before card B-27 has no lookup table: read its source, say so (card W-17)
                    spec_name, cols = _PRE_B27[name]
                    files = _store_files(root, spec_name, want if name == "clus" else {}, manifest)
                    if files:
                        slow_path_note(release, _RELATIONS[name][0], spec_name)
                if name == "dom":
                    rdom_files = _merged_round_files(root, "round_domains", release, manifest) if release is not None else []
                    rcov_files = _merged_round_files(root, "round_coverage", release, manifest) if release is not None else []
                    libraries = {str(p["library"]) for _, p in (*files, *rdom_files, *rcov_files) if "library" in p}
                    for _, p in files:
                        if "stage" in p:
                            annotated.setdefault(str(p["library"]), set()).add(str(p["stage"]))
                gaps = _unlabelled(root, manifest, files) if name == "lab" and spec_name == "contig_label" else None
                if name == "dom":
                    _dom_view(self.con, files, rdom_files)
                    _view(self.con, "rcov", "round_coverage", rcov_files, ("protein_uid", "n_rows"), ("library",))
                    files = files + rdom_files + rcov_files
                elif spec_name == "contig_occurrence":  # one row per source record: the label's n_records is 1
                    _view(self.con, f"{name}_raw", spec_name, files, cols, keep)
                    self.con.execute(f"CREATE TEMP VIEW {name} AS SELECT *, 1::UINTEGER AS n_records FROM {name}_raw")
                elif gaps:  # card W-36: some stages have label files, others not
                    stages, occ = gaps
                    _partial_label_note(release, stages)
                    _view(self.con, f"{name}_lab", spec_name, files, cols, keep)
                    occ_cols = _PRE_B27[name][1]
                    _view(self.con, f"{name}_occ", "contig_occurrence", occ, occ_cols, keep)
                    names = ", ".join(f'"{c}"' for c in cols)
                    self.con.execute(f"CREATE TEMP VIEW {name} AS SELECT {names} FROM {name}_lab UNION ALL "
                                     f"SELECT {names} FROM (SELECT *, 1::UINTEGER AS n_records FROM {name}_occ)")
                    files = files + occ
                else:
                    _view(self.con, name, spec_name, files, cols, keep)
                largest = max(largest, row_groups_of(p for p, _ in files))
            set_join_key_filter(self.con, largest)
        except BaseException:
            self.__exit__(None, None, None)
            raise
        return _Store(self.con, sorted(libraries), {lib: tuple(sorted(ss)) for lib, ss in sorted(annotated.items())})

    def __exit__(self, *exc) -> None:
        self.con.close()
        shutil.rmtree(self.scratch, ignore_errors=True)


def _stream(con, sql: str) -> Iterator[dict]:
    """The rows of one query, a batch at a time through an Arrow reader (as `biobox.io.duck_batches`)."""
    reader = con.sql(sql).to_arrow_reader(_BATCH_ROWS)
    try:
        for batch in reader:
            yield from batch.to_pylist()
    finally:
        reader.close()


def _annotate_proteins(con, library: str | None, annotated_stages: tuple[str, ...] = ()) -> None:
    """Temp tables `_plen`, `_cl` (NR90 and NR50 cluster uids, the NR50 representative), `_arch` (domain
    names of `library` by env_start, then model accession, a run of one name written once, joined by `+`) and
    `_dstate` (card C-38b: `measured`, `inferred`, `none` or `not scanned`) for the protein uids of temp table
    `_p`. Each base relation is filtered by key before any join."""
    con.execute("CREATE OR REPLACE TEMP TABLE _plen AS SELECT protein_uid, max(length) AS length FROM prot "
                "WHERE protein_uid IN (SELECT protein_uid FROM _p) GROUP BY protein_uid")
    con.execute("""CREATE OR REPLACE TEMP TABLE _cl AS
                   SELECT member_uid AS protein_uid, min(cluster_uid) FILTER (WHERE level = 'nr90') AS nr90,
                          min(cluster_uid) FILTER (WHERE level = 'nr50') AS nr50,
                          min(rep_uid) FILTER (WHERE level = 'nr50') AS nr50_rep
                   FROM clus WHERE member_uid IN (SELECT protein_uid FROM _p) AND level IN ('nr90', 'nr50')
                   GROUP BY member_uid""")
    lib = "NULL" if library is None else _q(library)
    con.execute(f"""CREATE OR REPLACE TEMP TABLE _arch AS
                    WITH st AS (SELECT DISTINCT protein_uid FROM dom
                                WHERE library = {lib} AND source = 'stage'
                                  AND protein_uid IN (SELECT protein_uid FROM _p)),
                         d AS (SELECT DISTINCT protein_uid, env_start, model_acc, model_name FROM dom
                               WHERE library = {lib} AND protein_uid IN (SELECT protein_uid FROM _p)
                                 AND (source = 'stage' OR protein_uid NOT IN (SELECT protein_uid FROM st))),
                         l AS (SELECT *, lag(model_name) OVER (PARTITION BY protein_uid
                                                               ORDER BY env_start, model_acc) AS prev FROM d)
                    SELECT protein_uid, string_agg(model_name, '+' ORDER BY env_start, model_acc) AS arch
                    FROM l WHERE prev IS DISTINCT FROM model_name GROUP BY protein_uid""")
    if library is None:
        con.execute("CREATE OR REPLACE TEMP TABLE _dstate AS "
                    "SELECT DISTINCT protein_uid, 'not scanned' AS domain_state FROM _p")
    else:
        cov = (f"CASE WHEN stage IN ({', '.join(_q(s) for s in annotated_stages)}) THEN 1 ELSE 0 END"
               if annotated_stages else "0")
        con.execute(f"""CREATE OR REPLACE TEMP TABLE _dstate AS
                        SELECT p.protein_uid,
                               CASE coalesce(max(r.rank), 0)
                                 WHEN 3 THEN 'measured'
                                 WHEN 2 THEN 'inferred'
                                 WHEN 1 THEN 'none'
                                 ELSE 'not scanned'
                               END AS domain_state
                        FROM (SELECT DISTINCT protein_uid FROM _p) p
                        LEFT JOIN (
                            SELECT protein_uid, CASE WHEN inferred_from IS NULL THEN 3 ELSE 2 END AS rank
                            FROM dom WHERE library = {lib} AND protein_uid IN (SELECT protein_uid FROM _p)
                            UNION ALL
                            SELECT protein_uid, CASE WHEN n_rows > 0 THEN 3 ELSE 1 END AS rank
                            FROM rcov WHERE library = {lib} AND protein_uid IN (SELECT protein_uid FROM _p)
                            UNION ALL
                            SELECT protein_uid, {cov} AS rank
                            FROM prot WHERE protein_uid IN (SELECT protein_uid FROM _p)
                        ) r ON r.protein_uid = p.protein_uid
                        GROUP BY p.protein_uid""")


def _annotate_genomes(con) -> None:
    """Temp table `_gens` for the contigs of `_c` (card C-21c)."""
    con.execute("""CREATE OR REPLACE TEMP TABLE _gens AS
                   SELECT contig_uid, count(*) AS n_genomes,
                          string_agg(lower(hex(genome_uid)), ',' ORDER BY genome_uid) AS genome_uids
                   FROM (SELECT DISTINCT contig_uid, genome_uid FROM (
                           SELECT contig_uid, genome_uid FROM dim WHERE contig_uid IN (SELECT contig_uid FROM _c)
                           UNION ALL
                           SELECT contig_uid, genome_uid FROM gin WHERE contig_uid IN (SELECT contig_uid FROM _c)))
                   GROUP BY contig_uid""")


def _annotate_contigs(con, with_ordinals: bool, genomes: bool = False) -> None:
    """Temp tables for the contig uids of temp table `_c`: `_lab` (accession and header of the contig's first
    source record by (l0_object_sha256, record_ordinal), as the catalog's `contig_label`, and the number of
    records), `_clen` (default variant's length), `_dim` (genome_uid, tier of the release's `contig_dim`), with
    `genomes` (the relations `dim` and `gin` open) `_gens` (card C-21c: `n_genomes` and `genome_uids`, the contig's
    `contig_dim` genome and every genome of the release's `genome_input` that holds it, comma-separated hex in uid
    order; one row per contig, so no join fans out) and,
    `with_ordinals`, `_ord` (each slot's feature_uid and ordinal in the tape version's tapes).

    `_lab` reads `contig_label` (card B-27): one row per contig and stage, with that stage's first record and
    its records. The rows of the contigs are taken first, which prunes the table to the row groups that hold
    them (`catalog.set_join_key_filter`), and only then merged over the release's stages."""
    con.execute("""CREATE OR REPLACE TEMP TABLE _lab AS
                   SELECT contig_uid, arg_min(accession, (l0_object_sha256, record_ordinal)) AS accession,
                          arg_min(original_header, (l0_object_sha256, record_ordinal)) AS header,
                          sum(n_records)::BIGINT AS n_records
                   FROM lab WHERE contig_uid IN (SELECT contig_uid FROM _c) GROUP BY contig_uid""")
    con.execute("CREATE OR REPLACE TEMP TABLE _clen AS SELECT contig_uid, max(contig_len) AS contig_len FROM cidx "
                "WHERE is_default_variant AND contig_uid IN (SELECT contig_uid FROM _c) GROUP BY contig_uid")
    con.execute("CREATE OR REPLACE TEMP TABLE _dim AS SELECT contig_uid, min(genome_uid) AS genome_uid, "
                "min(tier) AS tier FROM dim WHERE contig_uid IN (SELECT contig_uid FROM _c) GROUP BY contig_uid")
    if genomes:
        _annotate_genomes(con)
    if with_ordinals:
        con.execute("""CREATE OR REPLACE TEMP TABLE _ord AS
                       SELECT DISTINCT * FROM (
                           SELECT contig_uid, unnest(feature_uids) AS feature_uid,
                                  unnest(range(len(feature_uids))) AS ordinal
                           FROM tape WHERE contig_uid IN (SELECT contig_uid FROM _c))""")


def _load_hits(con, hits: list[dict]) -> None:
    """Temp table `_h`: one row per hit protein with its sequence score and E-value at the reference Z."""
    if any(h["target_kind"] != "protein" for h in hits):
        kinds = sorted({h["target_kind"] for h in hits} - {"protein"})
        raise ValueError(f"the run has {', '.join(kinds)} hits; only protein hits have genes (use --json)")
    t = pa.table({"protein_uid": pa.array([h["target_uid"] for h in hits], pa.binary(16)),
                  "seq_bitscore": pa.array([h.get("seq_bitscore") for h in hits], pa.float64()),
                  "seq_evalue_zref": pa.array([h.get("seq_evalue_zref") for h in hits], pa.float64())})
    con.register("_hits_in", t)
    con.execute("CREATE OR REPLACE TEMP TABLE _h AS SELECT protein_uid, max(seq_bitscore) AS seq_bitscore, "
                "min(seq_evalue_zref) AS seq_evalue_zref FROM _hits_in GROUP BY protein_uid")
    con.unregister("_hits_in")


_GENE_ORDER = 'l.accession NULLS LAST, o.contig_uid, o.start, o."end", o.strand'


def _hit_genes(con, hits: list[dict], tape_version: str, library: str | None,
               prefer: dict[bytes, tuple[str | None, str | None]] | None = None,
               annotated_stages: tuple[str, ...] = ()) -> None:
    """Temp table `_genes`: one row per protein_occurrence of a hit protein with the columns of HITS_TSV_COLUMNS,
    plus `gene_rank` (1 for the protein's first gene in the table's order) and `n_genes`.

    `prefer` (card W-20): {protein_uid: (accession, header text)}, the contig a positive control names. A gene on a
    contig with that accession (or whose header contains that text) ranks before the protein's others, so the
    FASTA's `example=` can be found by the control's own accession."""
    _load_hits(con, hits)
    con.register("_pref_in", pa.table({
        "protein_uid": pa.array(sorted(prefer or {}), pa.binary(16)),
        "accession": pa.array([(prefer or {})[u][0] for u in sorted(prefer or {})], pa.string()),
        "text": pa.array([(prefer or {})[u][1] for u in sorted(prefer or {})], pa.string())}))
    con.execute("CREATE OR REPLACE TEMP TABLE _pref AS SELECT * FROM _pref_in")
    con.unregister("_pref_in")
    con.execute("CREATE OR REPLACE TEMP TABLE _p AS SELECT protein_uid FROM _h")
    # a gene's stage is the stage of its protein_occurrence file (card W-23); a contig's genes are written by one
    # stage only (DEC-26), and should two stages hold the same gene the first stage by name is kept
    con.execute("""CREATE OR REPLACE TEMP TABLE _occ AS
                   SELECT protein_uid, contig_uid, feature_uid, start, "end", strand, dist_to_contig_end_bp,
                          min(stage) AS stage
                   FROM occ WHERE protein_uid IN (SELECT protein_uid FROM _p)
                   GROUP BY protein_uid, contig_uid, feature_uid, start, "end", strand, dist_to_contig_end_bp""")
    con.execute("CREATE OR REPLACE TEMP TABLE _c AS SELECT DISTINCT contig_uid FROM _occ")
    _annotate_contigs(con, with_ordinals=True, genomes=True)
    _annotate_proteins(con, library, annotated_stages)
    # card W-23: the genome kind of each contig in each stage (genome_input; empty where a stage has none)
    con.execute("CREATE OR REPLACE TEMP TABLE _kind AS SELECT contig_uid, stage, min(kind) AS kind FROM gin "
                "WHERE contig_uid IN (SELECT contig_uid FROM _c) GROUP BY contig_uid, stage")
    lost = con.execute("SELECT o.feature_uid, o.contig_uid FROM _occ o ANTI JOIN _ord r "
                       "ON r.contig_uid = o.contig_uid AND r.feature_uid = o.feature_uid "
                       "ORDER BY o.contig_uid, o.feature_uid").fetchall()
    if lost:
        raise ValueError(f"{len(lost)} genes of hit proteins are not slots of tape_version {tape_version}; first: "
                         f"feature {lost[0][0].hex()} of contig {lost[0][1].hex()}")
    con.execute(f"""CREATE OR REPLACE TEMP TABLE _genes AS
                    SELECT l.accession, l.header, c.contig_len, o.start, o."end", o.strand, p.length,
                           h.seq_bitscore, h.seq_evalue_zref,
                           coalesce(a.arch, CASE WHEN coalesce(ds.domain_state, 'not scanned') = 'not scanned'
                                                 THEN '?' ELSE '-' END) AS pfam_architecture,
                           coalesce(ds.domain_state, 'not scanned') AS domain_state,
                           o.protein_uid, cl.nr90, cl.nr50, cl.nr50_rep, d.genome_uid, d.tier, o.contig_uid,
                           o.feature_uid, r.ordinal, {_q(tape_version)} AS tape_version, o.dist_to_contig_end_bp,
                           l.n_records, o.stage, k.kind, gs.n_genomes, gs.genome_uids,
                           row_number() OVER (PARTITION BY o.protein_uid ORDER BY
                               CASE WHEN pf.protein_uid IS NOT NULL
                                         AND (l.accession = pf.accession
                                              OR (pf.text IS NOT NULL AND contains(coalesce(l.header, ''), pf.text)))
                                    THEN 0 ELSE 1 END, {_GENE_ORDER}) AS gene_rank,
                           count(*) OVER (PARTITION BY o.protein_uid) AS n_genes
                    FROM _occ o JOIN _h h ON h.protein_uid = o.protein_uid
                    LEFT JOIN _lab l ON l.contig_uid = o.contig_uid
                    LEFT JOIN _clen c ON c.contig_uid = o.contig_uid
                    LEFT JOIN _dim d ON d.contig_uid = o.contig_uid
                    LEFT JOIN _ord r ON r.contig_uid = o.contig_uid AND r.feature_uid = o.feature_uid
                    LEFT JOIN _plen p ON p.protein_uid = o.protein_uid
                    LEFT JOIN _cl cl ON cl.protein_uid = o.protein_uid
                    LEFT JOIN _arch a ON a.protein_uid = o.protein_uid
                    LEFT JOIN _dstate ds ON ds.protein_uid = o.protein_uid
                    LEFT JOIN _pref pf ON pf.protein_uid = o.protein_uid
                    LEFT JOIN _kind k ON k.contig_uid = o.contig_uid AND k.stage = o.stage
                    LEFT JOIN _gens gs ON gs.contig_uid = o.contig_uid""")


def hits_tsv(root: Path, hits: list[dict], *, release: str | None, tape_version: str, gene_call_version: str,
             manifest: dict | None = None, library: str | None = None, threads: int | None = None,
             tmp=None) -> Iterator[str]:
    """`search hits --format tsv`: a column line, then one row per gene of a hit protein (`protein_occurrence`),
    rows by accession (genes of contigs without one last), contig_uid, start, end, strand, protein_uid.

    accession, header, n_records: the contig's first source record and how many it has (`contig_label`, card
    B-27: `contig_occurrence` by contig);
    contig_len: `contig_index`; length: the protein's (`protein_seq`); seq_bitscore, seq_evalue_zref: the
    protein's sequence score in the run and its E-value at the reference Z; pfam_architecture: domain names of
    the library (`_Store.library`) in order, `-` for none; nr90, nr50, nr50_rep: cluster uids and the NR50
    representative (`member_cluster` of the release: `cluster_member` by member); genome_uid, tier: `contig_dim` of
    the release; ordinal:
    the gene's slot in `tape_version`, beside its feature_uid (AGENTS §3.3); stage, kind (card W-23): the stage
    of the gene's `protein_occurrence` file and the `genome_input` kind of its contig in that stage; n_genomes,
    genome_uids (card C-21c, last): every genome that holds the contig (`_annotate_contigs`). One row per gene: a
    contig in several genomes lists them in genome_uids and is not repeated."""
    with open_store(root, release=release, gene_call_version=gene_call_version, tape_version=tape_version,
                    manifest=manifest, threads=threads, tmp=tmp) as st:
        lib = st.library(library)
        _hit_genes(st.con, hits, tape_version, lib, annotated_stages=st.stages_of(lib))
        # card W-20: `tier` is a column only when the release's contig_dim has one for some gene (the tiering step
        # does not exist yet, so contig_dim.tier is NULL in every release built so far: an always-empty column)
        cols = HITS_TSV_COLUMNS if st.con.execute("SELECT count(tier) FROM _genes").fetchall()[0][0] else \
            [c for c in HITS_TSV_COLUMNS if c != "tier"]
        yield "\t".join(cols)
        for r in _stream(st.con, 'SELECT * FROM _genes ORDER BY accession NULLS LAST, contig_uid, start, "end", '
                                 "strand, protein_uid"):
            yield _tsv_line(cols, r)


def _format_hits_faa_chunk(chunk: list[dict], seqs: dict[bytes, tuple[int, str]]) -> Iterator[str]:
    for r in chunk:
        got = seqs.get(r["protein_uid"])
        if got is None:
            raise ValueError(f"hit protein {r['protein_uid'].hex()} has no protein_seq row")
        desc = (f"length={cell('length', r['length'])} seq_bitscore={cell('seq_bitscore', r['seq_bitscore'])} "
                f"seq_evalue_zref={cell('seq_evalue_zref', r['seq_evalue_zref'])} n_genes={r['n_genes']}")
        if r["contig_uid"] is not None:
            desc += (f" example={r['accession'] or r['contig_uid'].hex()}:{r['start']}-{r['end']}"
                     f" strand={cell('strand', r['strand'])}")
            if r["header"]:
                desc += " " + cell("header", r["header"])
        yield f">{r['protein_uid'].hex()} {desc}"
        yield got[1]


def hits_faa(root: Path, hits: list[dict], *, release: str | None, tape_version: str, gene_call_version: str,
             manifest: dict | None = None, library: str | None = None, threads: int | None = None,
             tmp=None, prefer: dict[bytes, tuple[str | None, str | None]] | None = None) -> Iterator[str]:
    """`search hits --format faa`: one record per hit protein by protein_uid, named by its uid, the sequence as
    `protein_seq` stores it on one line. The description: length, seq_bitscore, seq_evalue_zref, the number of
    genes and one gene as `example=ACCESSION:start-end strand=+1 header` (start-end 0-based, half-open): the first
    gene of `hits_tsv`, or, for a protein in `prefer` (`_hit_genes`: the contig a positive control names), the
    first gene on that contig."""
    files = release_files(root, release, manifest)
    with open_store(root, release=release, gene_call_version=gene_call_version, tape_version=tape_version,
                    manifest=files, threads=threads, tmp=tmp) as st:
        con = st.con
        lib = st.library(library)
        _hit_genes(con, hits, tape_version, lib, prefer, st.stages_of(lib))
        sql = """SELECT h.protein_uid, p.length, h.seq_bitscore, h.seq_evalue_zref, coalesce(g.n_genes, 0) AS n_genes,
                        g.accession, g.contig_uid, g.start, g."end", g.strand, g.header
                 FROM _h h LEFT JOIN _genes g ON g.protein_uid = h.protein_uid AND g.gene_rank = 1
                 LEFT JOIN _plen p ON p.protein_uid = h.protein_uid
                 ORDER BY h.protein_uid"""
        chunk: list[dict] = []
        for r in _stream(con, sql):
            chunk.append(r)
            if len(chunk) >= _BATCH_ROWS:
                seqs, _ = protein_seqs(root, [x["protein_uid"] for x in chunk], release=release,
                                       gene_call_version=gene_call_version, manifest=files, threads=threads, tmp=tmp)
                yield from _format_hits_faa_chunk(chunk, seqs)
                chunk.clear()
        if chunk:
            seqs, _ = protein_seqs(root, [x["protein_uid"] for x in chunk], release=release,
                                   gene_call_version=gene_call_version, manifest=files, threads=threads, tmp=tmp)
            yield from _format_hits_faa_chunk(chunk, seqs)


def _locus_head(r: dict) -> tuple[str, list[str]]:
    title = cell("header", r["header"]) if r["header"] else r["contig_uid"].hex()
    facts = [("contig_uid", r["contig_uid"]), ("anchor_feature_uid", r["anchor_feature_uid"]),
             ("anchor_ordinal", r["anchor_ordinal"]), ("anchor_protein_uid", r["anchor_protein_uid"]),
             ("tape_version", r["tape_version"]), ("policy_hash", r["policy_hash"]),
             ("censored", bool(r["flags"] & 1)),
             ("n_genomes", r["n_genomes"]), ("genome_uids", r["genome_uids"])]  # card C-21c
    name = f"locus_{_SAFE.sub('_', r['accession'] or r['contig_uid'].hex())}_{r['anchor_feature_uid'].hex()}.tsv"
    return name, [f"# {title}", "# " + "\t".join(f"{k}={cell(k, v)}" for k, v in facts), "\t".join(LOCUS_TSV_COLUMNS)]


def _array_lines(found: dict) -> list[str]:
    """A locus table's CRISPR array rows (card E-14, `biobox.search.arrays`): a marker line with the window's base
    span, the column line, one row per array within the window (none: the two lines alone)."""
    from biobox.search import arrays

    lo, hi = found["window"]
    return ([f"{arrays.MARKER}\twindow_start={lo}\twindow_end={hi}", "\t".join(arrays.COLUMNS)]
            + [_tsv_line(arrays.COLUMNS, a) for a in found["arrays"]])


def locus_table_rows(lines: list[str]) -> int:
    """The data rows of one locus table (`locus_tables`'s lines): its gene rows, and its array rows when it has an
    array section (card E-14), without the header, facts, column and marker lines."""
    from biobox.search import arrays

    k = next((i for i, line in enumerate(lines) if line.startswith(arrays.MARKER)), None)
    return len(lines) - 3 if k is None else len(lines) - 5


def locus_tables(root: Path, rows: list[dict], *, release: str | None, gene_call_version: str,
                 hit_bits: dict[bytes, float] | None = None, manifest: dict | None = None,
                 library: str | None = None, threads: int | None = None,
                 tmp=None, arrays: dict | None = None) -> Iterator[tuple[str, list[str]]]:
    """`nbr --format tsv`: (file name, lines) per anchor locus of `rows` (`nbr.fetch`'s), loci by accession
    (contigs without one last), contig_uid and anchor ordinal; rows by offset.

    Lines: `# <header of the contig's first source record>` (the contig uid without one), `# key=value` facts
    of the locus (contig_uid, anchor_feature_uid, anchor_ordinal, anchor_protein_uid, tape_version,
    policy_hash, censored: a linear contig end cut the window; n_genomes, genome_uids: every genome that holds
    the contig, card C-21c), the column line and one row per slot.
    hit_seq_bitscore: the slot protein's score in the run (`hit_bits`); length_aa, pfam_architecture, nr90,
    nr50 as in `hits_tsv` (empty for a non-CDS slot); delta_bp: the gap to the anchor (empty at the anchor);
    crosses_origin: flag bit 1. Name: `locus_<accession>_<anchor_feature_uid>.tsv`, characters other than
    letters, digits, `.`, `_`, `-` as `_`.

    `arrays` (card E-14): `biobox.search.arrays.near`'s result for `rows`, or None (a release without CRISPR arrays:
    the tables are exactly as before). Given, each table ends with its array section after the gene rows: the line
    `# crispr_arrays<TAB>window_start=LO<TAB>window_end=HI` (the window's base span in the anchor's frame), the
    column line `arrays.COLUMNS` and one row per array within the window, no blank line (stdout separates tables
    with one)."""
    from biobox.registry.schemas import table

    with open_store(root, release=release, gene_call_version=gene_call_version, manifest=manifest, threads=threads,
                    tmp=tmp, relations=("prot", "dom", "clus", "lab", "cidx", "dim", "gin")) as st:
        con = st.con
        library = st.library(library)
        con.register("_nb_in", pa.Table.from_pylist(rows, schema=table("neighborhoods").schema))
        con.execute("CREATE TEMP TABLE _nb AS SELECT * FROM _nb_in")
        con.unregister("_nb_in")
        bits = sorted((hit_bits or {}).items())
        con.register("_hb_in", pa.table({"protein_uid": pa.array([u for u, _ in bits], pa.binary(16)),
                                         "seq_bitscore": pa.array([b for _, b in bits], pa.float64())}))
        con.execute("CREATE TEMP TABLE _hb AS SELECT * FROM _hb_in")
        con.unregister("_hb_in")
        con.execute("CREATE TEMP TABLE _p AS SELECT DISTINCT protein_uid FROM _nb WHERE protein_uid IS NOT NULL")
        con.execute("CREATE TEMP TABLE _c AS SELECT DISTINCT contig_uid FROM _nb")
        _annotate_contigs(con, with_ordinals=False, genomes=True)
        _annotate_proteins(con, library, st.stages_of(library))
        sql = """SELECT n.contig_uid, n.anchor_feature_uid, n.anchor_ordinal, n.anchor_protein_uid, n.tape_version,
                        n.policy_hash, n.flags, l.accession, l.header, n."offset", n.offset_rel, n.ordinal,
                        n.feature_uid, n.start, n."end", n.strand, n.protein_uid, p.length AS length_aa,
                        CASE WHEN n.protein_uid IS NULL THEN NULL
                             ELSE coalesce(a.arch, CASE WHEN coalesce(ds.domain_state, 'not scanned') = 'not scanned'
                                                        THEN '?' ELSE '-' END) END AS pfam_architecture,
                        CASE WHEN n.protein_uid IS NULL THEN NULL
                             ELSE coalesce(ds.domain_state, 'not scanned') END AS domain_state,
                        cl.nr90, cl.nr50, b.seq_bitscore AS hit_seq_bitscore, n.delta_bp,
                        (n.flags::INTEGER & 2) = 2 AS crosses_origin, gs.n_genomes, gs.genome_uids
                 FROM _nb n LEFT JOIN _lab l ON l.contig_uid = n.contig_uid
                 LEFT JOIN _gens gs ON gs.contig_uid = n.contig_uid
                 LEFT JOIN _plen p ON p.protein_uid = n.protein_uid LEFT JOIN _arch a ON a.protein_uid = n.protein_uid
                 LEFT JOIN _dstate ds ON ds.protein_uid = n.protein_uid
                 LEFT JOIN _cl cl ON cl.protein_uid = n.protein_uid LEFT JOIN _hb b ON b.protein_uid = n.protein_uid
                 ORDER BY l.accession NULLS LAST, n.contig_uid, n.anchor_ordinal, n."offset" """
        cur, name, lines = None, "", []
        for r in _stream(con, sql):
            key = (r["contig_uid"], r["anchor_feature_uid"])
            if key != cur:
                if cur is not None:
                    yield name, lines + ([] if arrays is None else _array_lines(arrays[cur]))
                cur = key
                name, lines = _locus_head(r)
            lines.append(_tsv_line(LOCUS_TSV_COLUMNS, r))
        if cur is not None:
            yield name, lines + ([] if arrays is None else _array_lines(arrays[cur]))


def protein_fasta(root: Path, *, release: str, gene_call_version: str, manifest: dict | None = None,
                  threads: int | None = None, tmp=None) -> Iterator[str]:
    """`export proteins`: every distinct protein of `protein_seq` (the release manifest's files), by
    protein_uid, named by its uid, the sequence as stored on one line. Without a manifest it raises
    ValueError at the call: readers use the manifest's file list, never a glob (AGENTS 3.9)."""
    if manifest is None:
        raise ValueError(f"export proteins needs manifests/release={release}.json (write it with "
                         f"`bio-box release manifest --release {release}`); readers never glob")
    return _protein_fasta(root, release, gene_call_version, manifest, threads, tmp)


def _protein_fasta(root, release, gene_call_version, manifest, threads, tmp) -> Iterator[str]:
    with open_store(root, release=release, gene_call_version=gene_call_version, manifest=manifest, threads=threads,
                    tmp=tmp, relations=("prot",)) as st:
        for r in _stream(st.con, "SELECT protein_uid, min(seq) AS seq FROM prot GROUP BY protein_uid "
                                 "ORDER BY protein_uid"):
            yield ">" + r["protein_uid"].hex()
            yield r["seq"]


def protein_seqs(root: Path, uids, *, release: str, gene_call_version: str, manifest=None, threads: int | None = None,
                 tmp=None) -> tuple[dict[bytes, tuple[int, str]], list[dict]]:
    """({protein_uid: (length, sequence)} for those of `uids` the release's `protein_seq` holds, and what the lookup
    read), for `bio-box seq --protein` (card W-23) and `hits_faa` (card B-65).

    Reads Parquet footers first, enforces `LY-02`/`LY-03` sort and min/max statistics on `protein_uid`, bisects the
    sorted unique keys against each row group's `[min, max]`, and decodes at most one file and one matching row group
    at a time so memory stays bounded regardless of release stage count or key count (AGENTS §4.1, §4.3, §4.4)."""
    import pyarrow.compute as pc
    import pyarrow.parquet as pq

    from biobox.search.nbr import _store_files

    files = release_files(root, release, manifest)
    wanted = sorted(set(uids))
    pfiles = sorted(_store_files(Path(root), "protein_seq", {"gene_call_version": gene_call_version}, files),
                    key=lambda x: str(x[0]))
    if not pfiles:
        return {}, []

    out: dict[bytes, tuple[int, str]] = {}
    files_read = 0
    rgs_scanned = 0
    rgs_total = 0
    rows_scanned = 0

    for path, _ in pfiles:
        md = pq.read_metadata(path)
        leaf = [i for i in range(md.num_columns) if md.schema.column(i).path == "protein_uid"]
        if not leaf:
            raise ValueError(f"{path}: no column 'protein_uid'")
        col_idx = leaf[0]
        rgs_total += md.num_row_groups
        prev_max = None
        matching: list[tuple[int, int, int]] = []
        for g in range(md.num_row_groups):
            rg = md.row_group(g)
            if rg.num_rows == 0:
                continue
            st = rg.column(col_idx).statistics
            if st is None or not st.has_min_max:
                raise ValueError(f"{path}: row group {g} has no min/max statistics on protein_uid; "
                                 "a key lookup cannot prune it (LY-03)")
            lo, hi = st.min, st.max
            if lo > hi or (prev_max is not None and lo < prev_max):
                raise ValueError(f"{path}: row group {g} is not sorted by protein_uid (LY-02)")
            prev_max = hi
            if wanted:
                a = bisect.bisect_left(wanted, lo)
                b = bisect.bisect_right(wanted, hi)
                if a < b:
                    matching.append((g, a, b))
        if not matching:
            continue
        files_read += 1
        pf = pq.ParquetFile(path, metadata=md, pre_buffer=False)
        try:
            for g, a, b in matching:
                rgs_scanned += 1
                rows_scanned += md.row_group(g).num_rows
                sub = pa.array(wanted[a:b], type=pa.binary(16))
                u_tbl = pf.read_row_group(g, columns=["protein_uid"], use_threads=False)
                has_hit = pc.any(pc.is_in(u_tbl.column("protein_uid"), value_set=sub)).as_py()
                del u_tbl
                if not has_hit:
                    del sub
                    continue
                for batch in pf.iter_batches(batch_size=_BATCH_ROWS, row_groups=[g],
                                             columns=["protein_uid", "length", "seq"], use_threads=False):
                    mask = pc.is_in(batch.column("protein_uid"), value_set=sub)
                    hit_batch = batch.filter(mask)
                    del mask
                    if hit_batch.num_rows == 0:
                        del hit_batch
                        continue
                    u_col = hit_batch.column("protein_uid").to_pylist()
                    n_col = hit_batch.column("length").to_pylist()
                    s_col = hit_batch.column("seq").to_pylist()
                    del hit_batch
                    for u, n, s in zip(u_col, n_col, s_col):
                        n_int, s_str = int(n), str(s)
                        prev = out.get(u)
                        if prev is None:
                            out[u] = (n_int, s_str)
                        else:
                            out[u] = (max(prev[0], n_int), min(prev[1], s_str))
                del sub
        finally:
            del pf
            pa.default_memory_pool().release_unused()

    scans = [{"files": [str(p) for p, _ in pfiles], "files_read": files_read, "row_groups_scanned": rgs_scanned,
              "row_groups_total": rgs_total, "rows_scanned": rows_scanned}]
    return out, scans


def contig_dims(root: Path, release: str | None, contigs: set[bytes], manifest: dict | None = None) -> dict[bytes, dict]:
    """{contig_uid: {genome_uid, tier, is_circular, genome_uids}} for the contigs (files from `release_files`):
    `genome_uid`, `tier`, `is_circular` from the release's `contig_dim` (null for a contig without a row), and
    `genome_uids` (card C-21c) every genome that holds the contig, in uid order: its `contig_dim` genome and every
    genome the release's `genome_input` lists it in. A contig in neither table is absent."""
    import pyarrow.parquet as pq

    from biobox.search.nbr import _store_files

    if release is None or not contigs:
        return {}
    files = release_files(root, release, manifest)
    out: dict[bytes, dict] = {}
    genomes: dict[bytes, set[bytes]] = {}
    for path, _ in _store_files(Path(root), "contig_dim", {"release": release}, files):
        t = pq.read_table(path, columns=["contig_uid", "genome_uid", "tier", "is_circular"], partitioning=None,
                          filters=[("contig_uid", "in", sorted(contigs))])
        for r in t.to_pylist():
            out.setdefault(r["contig_uid"], {"genome_uid": r["genome_uid"], "tier": r["tier"],
                                             "is_circular": r["is_circular"]})
            genomes.setdefault(r["contig_uid"], set()).add(r["genome_uid"])
    for path, _ in _store_files(Path(root), "genome_input", {}, files):
        t = pq.read_table(path, columns=["contig_uid", "genome_uid"], partitioning=None,
                          filters=[("contig_uid", "in", sorted(contigs))])
        for r in t.to_pylist():
            genomes.setdefault(r["contig_uid"], set()).add(r["genome_uid"])
    for cu, gs in genomes.items():
        out.setdefault(cu, {"genome_uid": None, "tier": None, "is_circular": None})["genome_uids"] = sorted(gs)
    return out


RUN_FILE = "RUN"


def run_file_text(run_id: str, release: str) -> str:
    """The text of a `RUN` file (card W-29): the run id and the release, one `key<TAB>value` per line, nothing else
    (no path, host or date), so that the file is the same wherever and whenever the export is written."""
    return f"run\t{run_id}\nrelease\t{release}\n"


def parse_run_file(text: str) -> tuple[str, str]:
    """(run id, release) of a `RUN` file's text; ValueError when it is not exactly what `run_file_text` writes."""
    rows = text.split("\n")
    if len(rows) == 3 and rows[2] == "":
        k = [r.split("\t") for r in rows[:2]]
        if [len(x) for x in k] == [2, 2] and k[0][0] == "run" and k[1][0] == "release" and k[0][1] and k[1][1]:
            return k[0][1], k[1][1]
    raise ValueError("not a RUN file: expected the two lines `run<TAB>RUN_ID` and `release<TAB>RELEASE`")


def write_run_file(directory: Path, run_id: str, release: str) -> Path:
    """Write `RUN` into `directory` (card W-29). A `RUN` already there that names this run and release is left as it
    is; one that names another is ValueError, before anything is written."""
    path = Path(directory) / RUN_FILE
    want = run_file_text(run_id, release)
    if path.is_file():
        try:
            have = parse_run_file(path.read_text(encoding="utf-8"))
        except (ValueError, UnicodeDecodeError) as e:
            raise ValueError(f"{path} is not a RUN file ({e})") from None
        if have != (run_id, release):
            raise ValueError(f"{directory} holds a RUN file that names run {have[0]} (release {have[1]}), not run "
                             f"{run_id} (release {release}): it is another run's output directory")
        return path
    return write_lines(want.splitlines(), path)


def write_lines(lines: Iterable[str], path: Path) -> Path:
    """Write the lines to `path` under a temporary name and rename it into place; nothing is left on error."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
            for line in lines:
                fh.write(line + "\n")
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    return path
