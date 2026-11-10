"""`--json` exports: NDJSON with versioned record schemas (operation O-EXPORT; SEARCH_PLAN SD-18; S-13).

Every record carries `"schema": "<kind>/<version>"`. The rule, checked against
`tests/registry/json_schema_lock.json`: keys may be added at the same version; removing or
renaming a key needs a version bump. Bytes are lowercase hex, timestamps ISO 8601 UTC.
Lists are capped at query time (`top_n`, default 50), so one stored result serves any cap.

An ordinal is never exported alone (AGENTS §3.3, card B-16): a record with an ordinal, or with
occurrences that have one, carries the `feature_uid` beside it and the `tape_version` it belongs to.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Iterable, Iterator

SCHEMAS: dict[str, dict] = {
    "search_hit": {"version": 1, "keys": [
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
                top_n: int | None = TOP_N, members: bool = False, tape_version: str | None = None) -> Iterator[dict]:
    """Best domain first per target, targets by descending sequence score. Each occurrence's
    `slot_ordinal` is valid inside `tape_version` (nbr.occurrences looks it up there)."""
    rows = sorted(hits, key=lambda r: (-(r["seq_bitscore"] or 0), r["target_uid"], -r["bitscore"]))
    for r in rows[:top_n] if top_n else rows:
        occ = (occurrences or {}).get(r["target_uid"], [])
        shown = occ if members else occ[:OCCURRENCES_SHOWN]
        yield record("search_hit", {**r, "n_occurrences": len(occ), "tape_version": tape_version, "occurrences": [
            {"contig_uid": o["contig_uid"], "slot_ordinal": o["slot_ordinal"], "feature_uid": o["feature_uid"],
             "genome_uid": o["genome_uid"], "start": o["start"], "end": o["end"], "strand": o["strand"],
             "tier": o["tier"]} for o in shown]})


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


SUMMARY_TSV_COLUMNS = {
    "pfam": ["level", "key", "unannotated", "library", "model_acc", "model_name"],
    "nr": ["level", "key", "unannotated", "cluster_uid", "rep_uid", "rep_length", "model_name"],
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
    its rows; offsets as `offset:count,...`), tables separated by a blank line."""
    for k, v in record("nbr_summary_header", summary["header"]).items():
        yield f"# {k}\t{json.dumps(v, separators=(',', ':'))}"
    for level in summary["header"]["levels"]:
        cols = SUMMARY_TSV_COLUMNS["pfam" if level == "pfam" else "nr"] + SUMMARY_TSV_STATS
        yield ""
        yield "\t".join(cols)
        for r in summary["rows"]:
            if r["level"] == level:
                yield "\t".join(_tsv_cell(c, r.get(c)) for c in cols)


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
