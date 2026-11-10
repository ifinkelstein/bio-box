"""Exports (operation O-EXPORT; SEARCH_PLAN SD-18; S-13; card W-05).

`--json`: NDJSON with versioned record schemas. Every record carries `"schema": "<kind>/<version>"`.
The rule, checked against `tests/registry/json_schema_lock.json`: keys may be added at the same
version; removing or renaming a key needs a version bump. Bytes are lowercase hex, timestamps ISO
8601 UTC. Lists are capped at query time (`top_n`, default 50), so one stored result serves any cap.

`search_hit/2` (card W-05): the nested occurrence objects are `contig_uid, ordinal, feature_uid, start,
end, strand, genome_uid, tier`. Card B-17 renamed `slot_ordinal` to `ordinal` and dropped `genome_uid`
and `tier` from the stored rows; version 2 records that change, and the two come back joined from the
run's release's `contig_dim` at export time (null when the release has no row for the contig).

An ordinal is never exported alone (AGENTS §3.3, card B-16): a record or table row with an ordinal
carries the `feature_uid` beside it and the `tape_version` it belongs to. Cluster ids are cluster uids,
never token ids (AGENTS §3.4); genome and tier come from `contig_dim` (AGENTS §3.6).

Tables a person opens (card W-05): `hits_tsv` (one row per gene of a hit protein), `hits_faa` (one
record per hit protein), `locus_tables` (one table per anchor locus) and `protein_fasta` (every distinct
protein of a release). Scores are rounded to one decimal: bit scores as `%.1f`, E-values as `%.1e`.
They are DuckDB queries over explicit file lists (the release manifest when there is one, else the
registry layout; AGENTS §3.9), each table filtered by key with `IN (SELECT ...)` before a join, and
the results stream out in batches.
"""

from __future__ import annotations

import json
import os
import re
import shutil
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
    (`contig_dims`, the release's `contig_dim` rows by contig), null for a contig it lacks."""
    rows = sorted(hits, key=lambda r: (-(r["seq_bitscore"] or 0), r["target_uid"], -r["bitscore"]))
    dims = dims or {}
    for r in rows[:top_n] if top_n else rows:
        occ = (occurrences or {}).get(r["target_uid"], [])
        shown = occ if members else occ[:OCCURRENCES_SHOWN]
        yield record("search_hit", {**r, "n_occurrences": len(occ), "tape_version": tape_version, "occurrences": [
            {"contig_uid": o["contig_uid"], "ordinal": o["ordinal"], "feature_uid": o["feature_uid"],
             "start": o["start"], "end": o["end"], "strand": o["strand"],
             "genome_uid": dims.get(o["contig_uid"], {}).get("genome_uid"),
             "tier": dims.get(o["contig_uid"], {}).get("tier")} for o in shown]})


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


# --- tables a person opens (card W-05) -------------------------------------------------------------------

HITS_TSV_COLUMNS = ["accession", "header", "contig_len", "start", "end", "strand", "length", "seq_bitscore",
                    "seq_evalue_zref", "pfam_architecture", "protein_uid", "nr90", "nr50", "nr50_rep", "genome_uid",
                    "tier", "contig_uid", "feature_uid", "ordinal", "tape_version", "dist_to_contig_end_bp",
                    "n_records"]
LOCUS_TSV_COLUMNS = ["offset", "offset_rel", "ordinal", "feature_uid", "start", "end", "strand", "protein_uid",
                     "length_aa", "pfam_architecture", "nr90", "nr50", "hit_seq_bitscore", "delta_bp", "crosses_origin"]
_BITS = {"seq_bitscore", "hit_seq_bitscore"}
_EVALUE = {"seq_evalue_zref"}
EXPORT_THREADS = 4  # DuckDB threads of an export without --threads (DuckDB's own default is every core of the host)
EXPORT_MEMORY_LIMIT = "2GB"  # larger sorts and aggregations spill into a scratch directory outside the data root
_BATCH_ROWS = 8192
_SAFE = re.compile(r"[^A-Za-z0-9._-]")

# The store relations the exports read: registry table, the partitions a file must match (values from
# open_store's arguments; a null value matches no file), the columns, and partition values kept as columns.
_RELATIONS = {
    "occ": ("protein_occurrence", ("gene_call_version",),
            ("protein_uid", "contig_uid", "feature_uid", "start", "end", "strand", "dist_to_contig_end_bp"), ()),
    "prot": ("protein_seq", ("gene_call_version",), ("protein_uid", "length", "seq"), ()),
    "dom": ("l1_domains", (), ("protein_uid", "model_acc", "model_name", "env_start"), ("library",)),
    "clus": ("cluster_member", ("release",), ("member_uid", "level", "cluster_uid", "rep_uid"), ()),
    "dim": ("contig_dim", ("release",), ("contig_uid", "genome_uid", "tier"), ()),
    "lab": ("contig_occurrence", (), ("contig_uid", "accession", "original_header", "l0_object_sha256",
                                      "record_ordinal"), ()),
    "cidx": ("contig_index", (), ("contig_uid", "contig_len", "is_default_variant"), ()),
    "tape": ("contig_tapes", ("tape_version",), ("contig_uid", "feature_uids"), ()),
}


def cell(col: str, v) -> str:
    """One TSV cell: empty for null; bit scores `%.1f`, E-values `%.1e`; bytes as lowercase hex; booleans
    `true`/`false`; tabs and line breaks inside text become spaces."""
    if v is None:
        return ""
    if col in _BITS:
        return f"{v:.1f}"
    if col in _EVALUE:
        return f"{v:.1e}"
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


class _Store:
    def __init__(self, con, libraries: list[str]):
        self.con, self.libraries = con, libraries

    def library(self, wanted: str | None) -> str | None:
        """The domain library of the architecture column: `wanted`, else the only one, else the only one whose
        name starts with `pfam`; None when there is none (every architecture is then `-`)."""
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
    version and tape version; files from `manifest` when given, else the registry layout. Threads and memory
    are capped, and spills go to a private scratch directory under `tmp` (system temp by default)."""

    def __init__(self, root: Path, *, release: str | None, gene_call_version: str, tape_version: str | None = None,
                 manifest: dict | None = None, threads: int | None = None, tmp=None,
                 relations: tuple[str, ...] = tuple(_RELATIONS)):
        self.args = (Path(root), release, gene_call_version, tape_version, manifest, threads, tmp, relations)

    def __enter__(self) -> _Store:
        import duckdb

        from biobox.search.nbr import _store_files

        root, release, gcv, tv, manifest, threads, tmp, relations = self.args
        self.scratch = Path(tempfile.mkdtemp(prefix="biobox-export-", dir=tmp))
        self.con = duckdb.connect(config={"threads": threads or EXPORT_THREADS, "memory_limit": EXPORT_MEMORY_LIMIT,
                                          "temp_directory": str(self.scratch / "duckdb")})
        values = {"gene_call_version": gcv, "release": release, "tape_version": tv}
        libraries: set[str] = set()
        try:
            for name in relations:
                spec_name, filt, cols, keep = _RELATIONS[name]
                want = {k: values[k] for k in filt}
                files = [] if None in want.values() else _store_files(root, spec_name, want, manifest)
                if name == "dom":
                    libraries = {str(p["library"]) for _, p in files}
                _view(self.con, name, spec_name, files, cols, keep)
        except BaseException:
            self.__exit__(None, None, None)
            raise
        return _Store(self.con, sorted(libraries))

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


def _annotate_proteins(con, library: str | None) -> None:
    """Temp tables `_plen`, `_cl` (NR90 and NR50 cluster uids, the NR50 representative) and `_arch` (domain
    names of `library` by env_start, then model accession, a run of one name written once, joined by `+`) for
    the protein uids of temp table `_p`. Each base relation is filtered by key before any join."""
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
                    WITH d AS (SELECT DISTINCT protein_uid, env_start, model_acc, model_name FROM dom
                               WHERE library = {lib} AND protein_uid IN (SELECT protein_uid FROM _p)),
                         l AS (SELECT *, lag(model_name) OVER (PARTITION BY protein_uid
                                                               ORDER BY env_start, model_acc) AS prev FROM d)
                    SELECT protein_uid, string_agg(model_name, '+' ORDER BY env_start, model_acc) AS arch
                    FROM l WHERE prev IS DISTINCT FROM model_name GROUP BY protein_uid""")


def _annotate_contigs(con, with_ordinals: bool) -> None:
    """Temp tables for the contig uids of temp table `_c`: `_lab` (accession and header of the contig's first
    source record by (l0_object_sha256, record_ordinal), as the catalog's `contig_label`, and the number of
    records), `_clen` (default variant's length), `_dim` (genome_uid, tier of the release's `contig_dim`) and,
    `with_ordinals`, `_ord` (each slot's feature_uid and ordinal in the tape version's tapes)."""
    con.execute("""CREATE OR REPLACE TEMP TABLE _lab AS
                   SELECT contig_uid, arg_min(accession, (l0_object_sha256, record_ordinal)) AS accession,
                          arg_min(original_header, (l0_object_sha256, record_ordinal)) AS header,
                          count(*) AS n_records
                   FROM lab WHERE contig_uid IN (SELECT contig_uid FROM _c) GROUP BY contig_uid""")
    con.execute("CREATE OR REPLACE TEMP TABLE _clen AS SELECT contig_uid, max(contig_len) AS contig_len FROM cidx "
                "WHERE is_default_variant AND contig_uid IN (SELECT contig_uid FROM _c) GROUP BY contig_uid")
    con.execute("CREATE OR REPLACE TEMP TABLE _dim AS SELECT contig_uid, min(genome_uid) AS genome_uid, "
                "min(tier) AS tier FROM dim WHERE contig_uid IN (SELECT contig_uid FROM _c) GROUP BY contig_uid")
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


def _hit_genes(con, hits: list[dict], tape_version: str, library: str | None) -> None:
    """Temp table `_genes`: one row per protein_occurrence of a hit protein with the columns of HITS_TSV_COLUMNS,
    plus `gene_rank` (1 for the protein's first gene in the table's order) and `n_genes`."""
    _load_hits(con, hits)
    con.execute("CREATE OR REPLACE TEMP TABLE _p AS SELECT protein_uid FROM _h")
    con.execute("""CREATE OR REPLACE TEMP TABLE _occ AS
                   SELECT DISTINCT protein_uid, contig_uid, feature_uid, start, "end", strand, dist_to_contig_end_bp
                   FROM occ WHERE protein_uid IN (SELECT protein_uid FROM _p)""")
    con.execute("CREATE OR REPLACE TEMP TABLE _c AS SELECT DISTINCT contig_uid FROM _occ")
    _annotate_contigs(con, with_ordinals=True)
    _annotate_proteins(con, library)
    lost = con.execute("SELECT o.feature_uid, o.contig_uid FROM _occ o ANTI JOIN _ord r "
                       "ON r.contig_uid = o.contig_uid AND r.feature_uid = o.feature_uid "
                       "ORDER BY o.contig_uid, o.feature_uid").fetchall()
    if lost:
        raise ValueError(f"{len(lost)} genes of hit proteins are not slots of tape_version {tape_version}; first: "
                         f"feature {lost[0][0].hex()} of contig {lost[0][1].hex()}")
    con.execute(f"""CREATE OR REPLACE TEMP TABLE _genes AS
                    SELECT l.accession, l.header, c.contig_len, o.start, o."end", o.strand, p.length,
                           h.seq_bitscore, h.seq_evalue_zref, coalesce(a.arch, '-') AS pfam_architecture,
                           o.protein_uid, cl.nr90, cl.nr50, cl.nr50_rep, d.genome_uid, d.tier, o.contig_uid,
                           o.feature_uid, r.ordinal, {_q(tape_version)} AS tape_version, o.dist_to_contig_end_bp,
                           l.n_records,
                           row_number() OVER (PARTITION BY o.protein_uid ORDER BY {_GENE_ORDER}) AS gene_rank,
                           count(*) OVER (PARTITION BY o.protein_uid) AS n_genes
                    FROM _occ o JOIN _h h ON h.protein_uid = o.protein_uid
                    LEFT JOIN _lab l ON l.contig_uid = o.contig_uid
                    LEFT JOIN _clen c ON c.contig_uid = o.contig_uid
                    LEFT JOIN _dim d ON d.contig_uid = o.contig_uid
                    LEFT JOIN _ord r ON r.contig_uid = o.contig_uid AND r.feature_uid = o.feature_uid
                    LEFT JOIN _plen p ON p.protein_uid = o.protein_uid
                    LEFT JOIN _cl cl ON cl.protein_uid = o.protein_uid
                    LEFT JOIN _arch a ON a.protein_uid = o.protein_uid""")


def hits_tsv(root: Path, hits: list[dict], *, release: str | None, tape_version: str, gene_call_version: str,
             manifest: dict | None = None, library: str | None = None, threads: int | None = None,
             tmp=None) -> Iterator[str]:
    """`search hits --format tsv`: a column line, then one row per gene of a hit protein (`protein_occurrence`),
    rows by accession (genes of contigs without one last), contig_uid, start, end, strand, protein_uid.

    accession, header, n_records: the contig's first source record and how many it has (`contig_occurrence`);
    contig_len: `contig_index`; length: the protein's (`protein_seq`); seq_bitscore, seq_evalue_zref: the
    protein's sequence score in the run and its E-value at the reference Z; pfam_architecture: domain names of
    the library (`_Store.library`) in order, `-` for none; nr90, nr50, nr50_rep: cluster uids and the NR50
    representative (`cluster_member` of the release); genome_uid, tier: `contig_dim` of the release; ordinal:
    the gene's slot in `tape_version`, beside its feature_uid (AGENTS §3.3)."""
    with open_store(root, release=release, gene_call_version=gene_call_version, tape_version=tape_version,
                    manifest=manifest, threads=threads, tmp=tmp) as st:
        _hit_genes(st.con, hits, tape_version, st.library(library))
        yield "\t".join(HITS_TSV_COLUMNS)
        for r in _stream(st.con, 'SELECT * FROM _genes ORDER BY accession NULLS LAST, contig_uid, start, "end", '
                                 "strand, protein_uid"):
            yield _tsv_line(HITS_TSV_COLUMNS, r)


def hits_faa(root: Path, hits: list[dict], *, release: str | None, tape_version: str, gene_call_version: str,
             manifest: dict | None = None, library: str | None = None, threads: int | None = None,
             tmp=None) -> Iterator[str]:
    """`search hits --format faa`: one record per hit protein by protein_uid, named by its uid, the sequence as
    `protein_seq` stores it on one line. The description: length, seq_bitscore, seq_evalue_zref, the number of
    genes and the first gene of `hits_tsv` as `example=ACCESSION:start-end(strand) header`."""
    with open_store(root, release=release, gene_call_version=gene_call_version, tape_version=tape_version,
                    manifest=manifest, threads=threads, tmp=tmp) as st:
        con = st.con
        _hit_genes(con, hits, tape_version, st.library(library))
        con.execute("CREATE OR REPLACE TEMP TABLE _seq AS SELECT protein_uid, min(seq) AS seq FROM prot "
                    "WHERE protein_uid IN (SELECT protein_uid FROM _p) GROUP BY protein_uid")
        sql = """SELECT h.protein_uid, p.length, h.seq_bitscore, h.seq_evalue_zref, coalesce(g.n_genes, 0) AS n_genes,
                        g.accession, g.contig_uid, g.start, g."end", g.strand, g.header, s.seq
                 FROM _h h LEFT JOIN _genes g ON g.protein_uid = h.protein_uid AND g.gene_rank = 1
                 LEFT JOIN _plen p ON p.protein_uid = h.protein_uid LEFT JOIN _seq s ON s.protein_uid = h.protein_uid
                 ORDER BY h.protein_uid"""
        for r in _stream(con, sql):
            if r["seq"] is None:
                raise ValueError(f"hit protein {r['protein_uid'].hex()} has no protein_seq row")
            desc = (f"length={cell('length', r['length'])} seq_bitscore={cell('seq_bitscore', r['seq_bitscore'])} "
                    f"seq_evalue_zref={cell('seq_evalue_zref', r['seq_evalue_zref'])} n_genes={r['n_genes']}")
            if r["contig_uid"] is not None:
                desc += (f" example={r['accession'] or r['contig_uid'].hex()}:{r['start']}-{r['end']}"
                         f"({r['strand']:+d})")
                if r["header"]:
                    desc += " " + cell("header", r["header"])
            yield f">{r['protein_uid'].hex()} {desc}"
            yield r["seq"]


def _locus_head(r: dict) -> tuple[str, list[str]]:
    title = cell("header", r["header"]) if r["header"] else r["contig_uid"].hex()
    facts = [("contig_uid", r["contig_uid"]), ("anchor_feature_uid", r["anchor_feature_uid"]),
             ("anchor_ordinal", r["anchor_ordinal"]), ("anchor_protein_uid", r["anchor_protein_uid"]),
             ("tape_version", r["tape_version"]), ("policy_hash", r["policy_hash"]),
             ("censored", bool(r["flags"] & 1))]
    name = f"locus_{_SAFE.sub('_', r['accession'] or r['contig_uid'].hex())}_{r['anchor_feature_uid'].hex()}.tsv"
    return name, [f"# {title}", "# " + "\t".join(f"{k}={cell(k, v)}" for k, v in facts), "\t".join(LOCUS_TSV_COLUMNS)]


def locus_tables(root: Path, rows: list[dict], *, release: str | None, gene_call_version: str,
                 hit_bits: dict[bytes, float] | None = None, manifest: dict | None = None,
                 library: str | None = None, threads: int | None = None,
                 tmp=None) -> Iterator[tuple[str, list[str]]]:
    """`nbr --format tsv`: (file name, lines) per anchor locus of `rows` (`nbr.fetch`'s), loci by accession
    (contigs without one last), contig_uid and anchor ordinal; rows by offset.

    Lines: `# <header of the contig's first source record>` (the contig uid without one), `# key=value` facts
    of the locus (contig_uid, anchor_feature_uid, anchor_ordinal, anchor_protein_uid, tape_version,
    policy_hash, censored: a linear contig end cut the window), the column line and one row per slot.
    hit_seq_bitscore: the slot protein's score in the run (`hit_bits`); length_aa, pfam_architecture, nr90,
    nr50 as in `hits_tsv` (empty for a non-CDS slot); delta_bp: the gap to the anchor (empty at the anchor);
    crosses_origin: flag bit 1. Name: `locus_<accession>_<anchor_feature_uid>.tsv`, characters other than
    letters, digits, `.`, `_`, `-` as `_`."""
    from biobox.registry.schemas import table

    with open_store(root, release=release, gene_call_version=gene_call_version, manifest=manifest, threads=threads,
                    tmp=tmp, relations=("prot", "dom", "clus", "lab", "cidx", "dim")) as st:
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
        _annotate_contigs(con, with_ordinals=False)
        _annotate_proteins(con, library)
        sql = """SELECT n.contig_uid, n.anchor_feature_uid, n.anchor_ordinal, n.anchor_protein_uid, n.tape_version,
                        n.policy_hash, n.flags, l.accession, l.header, n."offset", n.offset_rel, n.ordinal,
                        n.feature_uid, n.start, n."end", n.strand, n.protein_uid, p.length AS length_aa,
                        CASE WHEN n.protein_uid IS NULL THEN NULL ELSE coalesce(a.arch, '-') END AS pfam_architecture,
                        cl.nr90, cl.nr50, b.seq_bitscore AS hit_seq_bitscore, n.delta_bp,
                        (n.flags::INTEGER & 2) = 2 AS crosses_origin
                 FROM _nb n LEFT JOIN _lab l ON l.contig_uid = n.contig_uid
                 LEFT JOIN _plen p ON p.protein_uid = n.protein_uid LEFT JOIN _arch a ON a.protein_uid = n.protein_uid
                 LEFT JOIN _cl cl ON cl.protein_uid = n.protein_uid LEFT JOIN _hb b ON b.protein_uid = n.protein_uid
                 ORDER BY l.accession NULLS LAST, n.contig_uid, n.anchor_ordinal, n."offset" """
        cur, name, lines = None, "", []
        for r in _stream(con, sql):
            key = (r["contig_uid"], r["anchor_feature_uid"])
            if key != cur:
                if cur is not None:
                    yield name, lines
                cur = key
                name, lines = _locus_head(r)
            lines.append(_tsv_line(LOCUS_TSV_COLUMNS, r))
        if cur is not None:
            yield name, lines


def protein_fasta(root: Path, *, release: str, gene_call_version: str, manifest: dict | None = None,
                  threads: int | None = None, tmp=None) -> Iterator[str]:
    """`export proteins`: every distinct protein of `protein_seq` (the release manifest's files when given, else
    every file of the gene-call version), by protein_uid, named by its uid, the sequence as stored on one line."""
    with open_store(root, release=release, gene_call_version=gene_call_version, manifest=manifest, threads=threads,
                    tmp=tmp, relations=("prot",)) as st:
        for r in _stream(st.con, "SELECT protein_uid, min(seq) AS seq FROM prot GROUP BY protein_uid "
                                 "ORDER BY protein_uid"):
            yield ">" + r["protein_uid"].hex()
            yield r["seq"]


def contig_dims(root: Path, release: str | None, contigs: set[bytes], manifest: dict | None = None) -> dict[bytes, dict]:
    """{contig_uid: {genome_uid, tier}} from the release's `contig_dim` for the contigs (files from `manifest`
    when given, else the registry layout); a contig without a row is absent."""
    import pyarrow.parquet as pq

    from biobox.search.nbr import _store_files

    if release is None or not contigs:
        return {}
    out: dict[bytes, dict] = {}
    for path, _ in _store_files(Path(root), "contig_dim", {"release": release}, manifest):
        t = pq.read_table(path, columns=["contig_uid", "genome_uid", "tier"], partitioning=None,
                          filters=[("contig_uid", "in", sorted(contigs))])
        for r in t.to_pylist():
            out.setdefault(r["contig_uid"], {"genome_uid": r["genome_uid"], "tier": r["tier"]})
    return out


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
