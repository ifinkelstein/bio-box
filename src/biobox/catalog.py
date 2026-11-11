"""Release manifest and generated read-only DuckDB catalog (operation O-CATALOG).

The manifest is the system of record (data_layer_and_tests.md §3.3): every
Parquet file of a release with its table, partitions, sha256, bytes, rows,
row groups and schema hash, and every artifact file (the bgzip sequence shards
and their indexes) with its sha256 and bytes. The catalog is disposable: a
DuckDB file holding only views over explicit file lists (no globs, LY-10) and
the sequence macros of §2A.3. It is built under a temporary name and moved into
place with rename(2); `current.json` is a symlink swapped the same way.

DuckDB cannot read BGZF, so the sequence macros are thin wrappers over Python functions
that `open_catalog` registers on each connection (DEC-17): `seq_get` and `seq_slice` call
`bb_seq`, `seq_get_circular` calls `bb_seq_circular`, `seq_raw` calls `bb_seq_raw`. They
resolve keys through `SeqStore` point lookups, never through SQL on `contig_index` or
`seq_locator`, so a fetch reads O(1) index row groups also when the keys come from a joined
table (AGENTS §4.4, card B-14). `partners_of` calls `bb_partner_counts`, a point lookup of
`tape_tokens` (card B-31), once per query with every partner slot, so each row group is read once
(card B-38).

The stored catalog holds `bb_seq`, `bb_seq_circular`, `bb_seq_raw` and `bb_partner_counts` as stand-in
macros that raise an error naming `biobox.catalog.open_catalog` (card B-34), so a catalog opened with
plain `duckdb` has every view and macro, and the sequence macros and `partners_of` fail there with that
message. `open_catalog` shadows each stand-in on its connection with a temporary macro of the same name
over the real Python function (`<name>_udf`). The lookup macros `default_digest`, `digest_len`,
`contig_len_of` and `token_of` work on a plain connection as SQL that scans the index views or
`token_dict`; `open_catalog` shadows them too, with point lookups over footer statistics (the
`SeqStore` for the first three, card B-31; `TokenKeyLookup` over `token_dict` for `token_of`, card
B-34), so a joined key reads O(1) row groups, not all.

Lookups by a key that a table is not sorted by cannot prune (card B-27): `contig_label` and `member_cluster`
are lookup tables sorted by `contig_uid` and by (`level`, `member_uid`), written beside `contig_occurrence` and
`cluster_member`. The catalog's `contig_label` view reads them (`label_view_sql`), and `open_catalog` raises
the connection's join filter so that keys from a joined table prune to the row groups that hold them
(`set_join_key_filter`).

`open_catalog` attaches the catalog file read-only to a private in-memory DuckDB instance (card B-34).
Python functions are registered in an instance's catalog, and `duckdb.connect` on one file path hands
every connection of the process the same instance, so a second connection could not register them
again ("Scalar Function with name bb_seq already exists"). With an instance per connection, any number
of connections to one catalog, or to several, live side by side in one process, each with its own
functions, settings and caches, and closing one leaves the others working.

A release names its stages (card B-25): `define_release` writes its `release_def` once, and
`build_manifest` lists only the files of those stages, so a stage written later, or being written
while the release is built, never enters it. A file under `release=<R>` is R's only (card B-33): a
release's manifest lists its stages' files, its own `release=` files, the `tape_version_def` of each tape
version those use, and the unpartitioned files, so a later release's files neither enter an earlier
release's manifest nor stop it verifying.

    stages = define_release(root, "rS0.1", ["S0"])   # later: release_stages(root, "rS0.1")
    manifest = build_manifest(root, "rS0.1", stages)
    write_manifest(manifest, root / "manifests/release=rS0.1.json")
    generate_catalog(manifest, root / "catalog/release=rS0.1.duckdb")
    swap_current(root / "manifests", "rS0.1")
    con = open_catalog(root / "catalog/release=rS0.1.duckdb", local_dir="/tmp/if2453")

Query-time readers of a release (the search pack, `nbr`, the search worker, `SeqStore`) open only the
release's files too (card B-28): `ReleaseFiles.of(root, release)` lists them from the release's manifest
when it has one, else from the stages its `release_def` names, one explicit stage directory at a time.

    files = ReleaseFiles.of(root, "rS0.1")
    files.table("protein_occurrence", gene_call_version="pyrodigal-3.6.3")   # [Path, ...]
"""

from __future__ import annotations

import bisect
import hashlib
import json
import os
import re
import shutil
import tempfile
import threading
from collections import OrderedDict
from pathlib import Path

import blake3
import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from biobox import ids, io, layout
from biobox.l0seq.lookup import KeyLookup
from biobox.registry.schemas import ARTIFACTS, TABLES, ArtifactSpec, TableSpec

MANIFEST_VERSION = 3  # 3: the manifest records the release's stages (card B-25)

# Partitions whose value is derived from a file column (biobox.layout buckets): the column and
# how a value maps to its bucket key. A prefix of the sort key is checked from row-group
# min/max; the low byte of a token is not monotone in the sort order, so it is checked by
# reading the column's distinct values.
DERIVED_PARTITIONS = {
    "protein_prefix": ("protein_uid", "prefix"),
    "token_bucket": ("token_id", "low_byte"),
}

# The pair macros' censoring window (card B-11): an anchor on a linear contig is flagged when its +-10-slot
# window is cut by a contig end, whatever k the query asks for (K_MAX of the pair oracle, tests/oracles).
PAIR_K_MAX = 10
PAIR_FLAG_CENSORED, PAIR_FLAG_CROSSES_ORIGIN = 1, 2
# partners_of reads the tape tokens of every anchor's contig; above this many anchor postings it refuses
PARTNERS_MAX_ROWS = 100_000

# The Python functions behind the catalog's macros (card B-34): name -> (parameters, SQL return type, the
# macros that need it). The stored catalog holds each as a stand-in macro that raises; `open_catalog` shadows
# it with a temporary macro over the function `<name>_udf`.
PY_FUNCTIONS = {
    "bb_seq": ("_cuid, _s, _e, _strand", "VARCHAR", "seq_get and seq_slice"),
    "bb_seq_circular": ("_cuid, _s, _e, _strand", "VARCHAR", "seq_get_circular"),
    "bb_seq_raw": ("_rd, _s, _e", "VARCHAR", "seq_raw"),
    # a whole query per call (card B-38; card B-31's bb_slot_tokens was called per contig)
    "bb_partner_counts": ("_w", "STRUCT(t UBIGINT, l BIGINT, p BIGINT)[]", "partners_of"),
}


# Trees under the data root that are not release files: run-step records (runs/stage=<S>/<step>.json,
# card C-02), which hold timestamps and host names; and two trees of roots written before card W-03,
# whose search path wrote the search pack (search_pack/release=<R>/) and search runs
# (search/run=<run_id>/ and search/run_index.json) into the root. Since W-03 nothing writes either
# tree: the pack and the runs live outside the data root (`biobox.search.runner`, "Where things
# live"). They stay listed so that an old root still builds and verifies its manifest.
EXCLUDED_TREES = ("search_pack", "search", "runs")


class LayoutError(ValueError):
    """A file under the data root that does not fit the registry (stray, bad path, bad schema)."""


def schema_hash(schema: pa.Schema) -> str:
    lines = [f"{f.name}\t{f.type}\t{'null' if f.nullable else 'notnull'}" for f in schema]
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


def compatible_columns(file_schema: pa.Schema, registry: pa.Schema) -> int | None:
    """Number of registry columns a file holds, or None if incompatible (LY-09).

    A file written before a nullable column was appended holds a prefix of the registry's
    fields; every field it lacks must be nullable. Anything else (a removed, renamed,
    retyped or reordered column) is incompatible.
    """
    k = len(file_schema)
    if k > len(registry) or not file_schema.equals(pa.schema(list(registry)[:k]), check_metadata=False):
        return None
    if not all(f.nullable for f in list(registry)[k:]):
        return None
    return k


def classify(rel: str) -> tuple[TableSpec, dict[str, str]] | None:
    """Map a root-relative path to (table, partitions). None: no table claims the path."""
    for spec in TABLES.values():
        loc = spec.location
        if not loc:
            continue
        if loc.endswith(".parquet"):
            if rel == loc:
                return spec, {}
            continue
        if not rel.startswith(loc + "/"):
            continue
        parts = rel[len(loc) + 1 :].split("/")
        dirs, fname = parts[:-1], parts[-1]
        if not fname.endswith(".parquet"):
            raise LayoutError(f"{rel}: not a .parquet file")
        kv = [d.split("=", 1) for d in dirs]
        if any(len(x) != 2 for x in kv) or tuple(k for k, _ in kv) not in key_sets(spec):
            raise LayoutError(f"{rel}: partition path {dirs} does not match {spec.name} keys "
                              f"{' or '.join(map(str, key_sets(spec)))} (LY-05)")
        return spec, dict(kv)
    return None


def key_sets(spec: TableSpec) -> tuple[tuple[str, ...], ...]:
    """The partition key tuples a path of `spec` may use: its `partition_keys`, then its
    `alt_partition_keys` (the ledger: per stage, or per release for a release's own steps, card B-28)."""
    return (spec.partition_keys, *spec.alt_partition_keys)


def all_keys(spec: TableSpec) -> tuple[str, ...]:
    """Every partition key of `spec` over all its key sets, in first-seen order (its view's key columns)."""
    return tuple(dict.fromkeys(k for ks in key_sets(spec) for k in ks))


def classify_artifact(rel: str) -> tuple[ArtifactSpec, dict[str, str], str] | None:
    """Map a path under an artifact location to (artifact, partitions, companion suffix or "")."""
    for spec in ARTIFACTS.values():
        if spec.single_file:
            if rel == spec.location:
                return spec, {}, ""
            continue
        if not rel.startswith(spec.location + "/"):
            continue
        parts = rel[len(spec.location) + 1 :].split("/")
        dirs, fname = parts[:-1], parts[-1]
        kv = [d.split("=", 1) for d in dirs]
        if any(len(x) != 2 for x in kv) or tuple(k for k, _ in kv) != spec.partition_keys:
            raise LayoutError(f"{rel}: partition path {dirs} does not match {spec.name} keys {spec.partition_keys} (LY-05)")
        for suffix in ("",) + spec.companions:
            if re.fullmatch(spec.pattern + re.escape(suffix), fname):
                return spec, dict(kv), suffix
        raise LayoutError(f"{rel}: not a {spec.name} file (LY-10)")
    return None


def _path_parts(rel: str) -> dict[str, str]:
    """The `<k>=<v>` directories of a root-relative path, read from the path alone (no registry lookup)."""
    return dict(d.split("=", 1) for d in rel.split("/")[:-1] if "=" in d)


def _artifact_files(root: Path, spec: ArtifactSpec) -> list[Path]:
    """The files under the artifact's location (the one file itself for a `single_file` artifact); none
    for an artifact that is not a release file (`in_manifest` False: the manifests and catalogs)."""
    if not spec.in_manifest:
        return []
    base = Path(root) / spec.location
    if spec.single_file:
        return [base] if base.is_file() else []
    return sorted(p for p in base.rglob("*") if p.is_file()) if base.exists() else []


def _artifact_entries(root: Path, in_release) -> list[dict]:
    """The artifact files `in_release(name, partitions)` keeps; a file of another stage is skipped before
    it is classified, so a stage still being written cannot fail the manifest of a release without it."""
    out = []
    for spec in ARTIFACTS.values():
        found = {}
        for path in _artifact_files(root, spec):
            rel = path.relative_to(root).as_posix()
            if not in_release(spec.name, _path_parts(rel)):
                continue
            _, parts, suffix = classify_artifact(rel)
            found[rel] = suffix
            out.append({"path": rel, "artifact": spec.name, "partitions": parts,
                        "sha256": ids.sha256_file(path).hex(), "bytes": path.stat().st_size})
        for rel, suffix in found.items():
            if suffix == "":
                missing = [s for s in spec.companions if rel + s not in found]
                if missing:
                    raise LayoutError(f"{rel}: {spec.name} file without {', '.join(missing)}")
            elif rel.removesuffix(suffix) not in found:
                raise LayoutError(f"{rel}: index without its {spec.name} file")
    return out


def _check_derived_partitions(rel: str, parts: dict[str, str], md: pq.FileMetaData, path: Path) -> None:
    """LY-05: every row's key matches the derived partition value.

    For a key prefix, rows are sorted by key, so each row group's min and max check every row.
    """
    leaf = {md.schema.column(i).path: i for i in range(md.num_columns)}
    for key, value in parts.items():
        if key not in DERIVED_PARTITIONS:
            continue
        col, how = DERIVED_PARTITIONS[key]
        if how == "low_byte":
            import pyarrow.compute as pc

            for v in pc.unique(pq.read_table(path, columns=[col]).column(col)).to_pylist():
                if not layout.bucket_matches(value, bytes([v & 0xFF])):
                    raise LayoutError(f"{rel}: {key}={value} but {col} {v} has low byte {v & 0xFF:02x} (LY-05)")
            continue
        ci = leaf[col]
        for rg in range(md.num_row_groups):
            st = md.row_group(rg).column(ci).statistics
            if st is None or not st.has_min_max:
                raise LayoutError(f"{rel}: row group {rg} has no min/max on {col} (LY-03)")
            for v in (st.min, st.max):
                if not layout.bucket_matches(value, bytes(v)):
                    raise LayoutError(f"{rel}: {key}={value} but {col} starts with {bytes(v)[:1].hex()} (LY-05)")


# --- release definition (card B-25) --------------------------------------------


def release_def_path(root: Path, release: str) -> Path:
    return Path(root) / TABLES["release_def"].location / f"release={release}" / "part-0.parquet"


def check_stages(stages) -> list[str]:
    """`stages` as a list in the order given: stage names, none empty and none twice. A bare string
    is refused (it would read as its characters)."""
    if isinstance(stages, (str, bytes)) or not isinstance(stages, (list, tuple)):
        raise TypeError(f"stages must be a list of stage names, got {stages!r}")
    out = list(stages)
    if not all(isinstance(s, str) and s and "/" not in s and "=" not in s for s in out):
        raise ValueError(f"stages {out}: each must be a non-empty name without '/' or '='")
    if len(set(out)) != len(out):
        raise ValueError(f"stages {out}: a stage is listed twice")
    return out


def read_release_def(root: Path, release: str) -> list[str] | None:
    """The stages of `release` in its order, or None when it has no `release_def` file."""
    path = release_def_path(root, release)
    if not path.exists():
        return None
    t = pq.read_table(path, columns=["position", "stage"], partitioning=None)
    pos = t.column("position").to_pylist()
    if pos != list(range(len(pos))):
        raise ValueError(f"{path}: release_def positions {pos} are not 0..{len(pos) - 1}")
    return t.column("stage").to_pylist()


def release_stages(root: Path, release: str) -> list[str]:
    """The stages `release_def` names for `release`; FileNotFoundError when the release is not defined."""
    stages = read_release_def(root, release)
    if stages is None:
        raise FileNotFoundError(f"{release_def_path(root, release)}: release {release} has no release_def; "
                                "define it once with its stages (define_release, or --stages)")
    return stages


def verify_release_def(root: Path, release: str, stages) -> list[str]:
    """`stages` checked (`check_stages`) and, when `release` has a `release_def`, equal to it in order
    (ValueError otherwise). Writes nothing."""
    stages = check_stages(stages)
    have = read_release_def(root, release)
    if have is not None and have != stages:
        raise ValueError(f"{release_def_path(root, release)}: release_def of {release} lists stages {have}, "
                         f"not {stages}")
    return stages


def define_release(root: Path, release: str, stages, out: Path | None = None) -> list[str]:
    """Make sure `release` is defined over `stages`: write its `release_def` if it has none, verify it
    lists the same stages in the same order if it has one (ValueError otherwise). Returns the stages.

    Written once, as a tape version's `tape_version_def` is: defining a release again with the same
    stages is a no-op, so every release-level builder calls this after checking its inputs.

    `out` (card C-08) is the output root of the step that calls this (default `root`). A definition in
    `root` is verified and never written again; else one already in `out` (written earlier in the same
    step) is verified; else it is written under `out`, and run-step publishes it with the step."""
    stages = verify_release_def(root, release, stages)
    if out is not None and Path(out).resolve() != Path(root).resolve() and read_release_def(root, release) is None:
        stages = verify_release_def(out, release, stages)
        root = out  # the definition, if absent from the root, is the output root's to hold
    if not stages:
        raise ValueError(f"release {release}: a release needs at least one stage")
    path = release_def_path(root, release)
    if not path.exists():
        spec = TABLES["release_def"]
        io.write_table(pa.table({"position": pa.array(range(len(stages)), pa.uint16()),
                                 "stage": pa.array(stages, pa.string())}, schema=spec.schema), spec, path)
    return stages


# A release's scope (card B-33). A file belongs to a release by its partitions:
#   (a) `stage=<S>`: to every release whose `release_def` lists S;
#   (b) `release=<R>`: to R only (its tokens, postings, clusters, `contig_dim`, `release_def`, and the
#       ledger of its own steps, `ledger/release=<R>/`, card B-28);
#   (c) a version alone, neither stage nor release (`tape_version_def`): to the releases whose (a) and (b)
#       files use that version, through a `tape_version` partition or the tape version an `igr_version`
#       names (`<tape_version>-igr1`, card B-20), so that C-08 finds every definition it needs;
#   (d) no partition (`store_params`, `id_history`, `test_run`): to every release.
# The table of each partitioned file in (c), and the version key that places it.
VERSION_ONLY = {"tape_version_def": "tape_version"}


def _scope(name: str, parts: dict[str, str]) -> tuple[str, str] | None:
    """Where a file of table or artifact `name` with partitions `parts` falls in a release's scope:
    ("stage", S), ("release", R), (version key, value) for (c), or None for an unpartitioned file (d).
    A partitioned file in none of these raises LayoutError: no rule says which release it belongs to."""
    if "stage" in parts:
        return "stage", parts["stage"]
    if "release" in parts:
        return "release", parts["release"]
    if not parts:
        return None
    key = VERSION_ONLY.get(name)
    if key is None or set(parts) != {key}:
        raise LayoutError(f"{name} partitions {parts}: neither a stage, a release nor a version alone "
                          f"({', '.join(VERSION_ONLY)}): no release scope claims it (card B-33)")
    return key, parts[key]


def _used_tape_versions(partitions) -> set[str]:
    """The tape versions an iterable of files' partitions name: `tape_version`, and the tape version of an
    `igr_version` (a malformed one names none; C-08 reports it)."""
    from biobox.provenance import _igr_tape_version  # provenance parses igr versions (card B-20)

    out = set()
    for parts in partitions:
        if "tape_version" in parts:
            out.add(parts["tape_version"])
        tv = _igr_tape_version(parts["igr_version"]) if "igr_version" in parts else None
        if tv is not None:
            out.add(tv)
    return out


def _in_release(stages: set[str], release_id: str, tape_versions: set[str] = frozenset()):
    """Whether a file of table or artifact `name` with partitions `parts` is in the release's scope (above):
    a listed stage, the release's own `release=` partition, a tape version of `tape_versions`, or none."""

    def keep(name: str, parts: dict[str, str]) -> bool:
        scope = _scope(name, parts)
        if scope is None:
            return True
        kind, value = scope
        if kind == "stage":
            return value in stages
        if kind == "release":
            return value == release_id
        return value in tape_versions
    return keep


def manifest_path(root: Path, release: str) -> Path:
    """Where a release's manifest is written: `manifests/release=<R>.json` under the data root."""
    return Path(root) / f"manifests/release={release}.json"


class ReleaseFiles:
    """The files of one release that a query-time reader may open (card B-28; AGENTS §3.9).

    With a manifest (`of` reads `manifests/release=<R>.json` when it exists, or takes the one given),
    exactly the manifest's files and artifacts: a file written after it is not the release's. Without
    one, the files of the stages the release's `release_def` names: each listed stage's own directory
    of a table is listed (a glob over buckets or libraries inside it, never over stage directories), so
    a stage on disk the release does not list is never opened. A table partitioned by release is listed
    in the release's own `release=<R>` directory only (card B-33: tokens, postings, clusters, `contig_dim`,
    `ledger/release=<R>/`); other files outside any stage partition (`tape_version_def`, unpartitioned
    tables) are listed by their location, filtered by the partitions a caller asks for. Files come in the
    release's stage order, then part order.

        files = ReleaseFiles.of(root, "rS0.1")
        files.table("contig_tapes", tape_version=tv)      # [Path]
        files.entries("l1_domains")                       # [(Path, {"library": ..., "stage": ...})]
        files.artifact("seq_shards")                      # the shard files, without their indexes
    """

    def __init__(self, root: Path, release: str, stages, manifest: dict | None = None):
        self.root = Path(root)
        self.release = release
        self.stages = check_stages(stages)
        self.manifest = manifest
        self._rank = {s: i for i, s in enumerate(self.stages)}

    @classmethod
    def of(cls, root: Path, release: str, manifest: dict | None = None) -> "ReleaseFiles":
        """The release's files: from `manifest`, else from `manifests/release=<R>.json` if present, else from
        `release_stages` (FileNotFoundError when the release has neither a manifest nor a `release_def`)."""
        root = Path(root)
        if manifest is None and manifest_path(root, release).exists():
            manifest = read_manifest(manifest_path(root, release))
        if manifest is not None:
            if manifest["release_id"] != release:
                raise ValueError(f"manifest of release {manifest['release_id']} given for release {release}")
            return cls(root, release, manifest["stages"], manifest)
        return cls(root, release, release_stages(root, release))

    @property
    def source(self) -> str:
        """`manifest` or `release_def`: where the file lists come from."""
        return "manifest" if self.manifest is not None else "release_def"

    def _order(self, path: Path, parts: dict) -> tuple:
        return (self._rank.get(parts.get("stage"), -1), io.part_key(path))

    def entries(self, name: str, **want) -> list[tuple[Path, dict]]:
        """(path, partitions) of table `name`'s files of the release whose partitions include `want`."""
        want = {k: str(v) for k, v in want.items()}
        if self.manifest is not None:
            found = [(self.root / f["path"], f["partitions"]) for f in self.manifest["files"] if f["table"] == name]
        else:
            found = self._layout(name, want)
        keep = [(p, parts) for p, parts in found if all(parts.get(k) == v for k, v in want.items())]
        return sorted(keep, key=lambda x: self._order(*x))

    def table(self, name: str, **want) -> list[Path]:
        """The paths of `entries(name, **want)`."""
        return [p for p, _ in self.entries(name, **want)]

    def _layout(self, name: str, want: dict[str, str]) -> list[tuple[Path, dict]]:
        """A table's files from the registry layout: one listing per key set and per listed stage."""
        spec = TABLES[name]
        base = self.root / spec.location
        if spec.location.endswith(".parquet"):
            return [(base, {})] if base.is_file() else []
        out = []
        fixed = {"release": self.release}  # a `release=` partition holds that release's files only (card B-33)
        # limitation: a version-only table (tape_version_def) is listed by `want`, not by the versions the release uses
        for keys in key_sets(spec):
            for stage in (self.stages if "stage" in keys else [None]):
                dirs = "/".join(f"{k}={stage if k == 'stage' else fixed.get(k, want.get(k, '*'))}" for k in keys)
                for p in io.part_files(base, f"{dirs}/*.parquet" if dirs else "*.parquet"):
                    hit = classify(p.relative_to(self.root).as_posix())
                    if hit is not None and hit[0].name == name:
                        out.append((p, hit[1]))
        return out

    def artifact(self, name: str) -> list[Path]:
        """The files of artifact `name` in the release (its main files; companions such as a shard's `.fai`
        and `.gzi` are left out), in stage order, then by path."""
        spec = ARTIFACTS[name]
        if self.manifest is not None:
            found = [(self.root / a["path"], a["partitions"]) for a in self.manifest.get("artifacts", [])
                     if a["artifact"] == name]
        elif spec.single_file:
            base = self.root / spec.location
            found = [(base, {})] if base.is_file() else []
        else:
            found = []
            for stage in (self.stages if "stage" in spec.partition_keys else [None]):
                dirs = "/".join(f"{k}={stage if k == 'stage' else '*'}" for k in spec.partition_keys)
                base = self.root / spec.location
                for p in sorted(base.glob(f"{dirs}/*" if dirs else "*")) if base.exists() else []:
                    if p.is_file() and re.fullmatch(spec.pattern, p.name):
                        found.append((p, classify_artifact(p.relative_to(self.root).as_posix())[1]))
        main = [(p, parts) for p, parts in found if re.fullmatch(spec.pattern, p.name)]
        return [p for p, parts in sorted(main, key=lambda x: (self._rank.get(x[1].get("stage"), -1), str(x[0])))]


def build_manifest(root: Path, release_id: str, stages) -> dict:
    """The release's manifest: every .parquet and artifact file under `root` in the release's scope (card
    B-33, `VERSION_ONLY` above): the files of `stages` (the release's, card B-25), its own `release=<R>`
    files, the `tape_version_def` of each tape version those files use, and the unpartitioned files. A
    stage directory the release does not list is ignored without being read, as are other releases' files
    and definitions of tape versions the release does not use; a listed stage with no file raises
    FileNotFoundError; `stages` must equal the release's `release_def` when it has one. Refuses strays
    and bad paths anywhere under `root`, and schema drift in the release's files."""
    root = Path(root).resolve()
    stages = verify_release_def(root, release_id, stages)
    hits = []
    for path in sorted(root.rglob("*.parquet")):
        rel = path.relative_to(root).as_posix()
        if rel.split("/", 1)[0] in EXCLUDED_TREES:
            continue
        hit = classify(rel)
        if hit is None:
            raise LayoutError(f"{rel}: no registered table claims this path (LY-10)")
        hits.append((path, rel, *hit))
    # (a), (b) and (d) first (no tape version given: no (c) file is kept): the tape versions their files
    # use decide which definitions (c) are the release's
    direct = _in_release(set(stages), release_id)
    used = _used_tape_versions(parts for _, _, spec, parts in hits if direct(spec.name, parts))
    in_release = _in_release(set(stages), release_id, used)
    picked = [h for h in hits if in_release(h[2].name, h[3])]
    present = {parts.get("stage") for _, _, _, parts in picked} | {
        s for spec in ARTIFACTS.values() if spec.partition_keys[:1] == ("stage",)
        for s in stages if (root / spec.location / f"stage={s}").is_dir()}
    missing = [s for s in stages if s not in present]
    if missing:
        raise FileNotFoundError(f"release {release_id}: listed stages {missing} have no files under {root}")
    files = []
    for path, rel, spec, parts in picked:
        pf = pq.ParquetFile(path)
        md = pf.metadata
        n_cols = compatible_columns(pf.schema_arrow, spec.schema)
        if n_cols is None:
            raise LayoutError(f"{rel}: schema differs from registry table {spec.name} (C-01)")
        _check_derived_partitions(rel, parts, md, path)
        files.append({
            "path": rel,
            "table": spec.name,
            "partitions": parts,
            "sha256": ids.sha256_file(path).hex(),
            "bytes": path.stat().st_size,
            "rows": md.num_rows,
            "row_groups": md.num_row_groups,
            "schema_hash": schema_hash(pf.schema_arrow),
            "n_columns": n_cols,
        })
    params_path = root / ARTIFACTS["store_params"].location
    params = json.loads(params_path.read_text()) if params_path.exists() else {}
    return {
        "manifest_version": MANIFEST_VERSION,
        "release_id": release_id,
        "stages": stages,
        "root": str(root),
        "params": params,
        "tables": {name: schema_hash(spec.schema) for name, spec in sorted(TABLES.items())},
        "files": files,
        "artifacts": _artifact_entries(root, in_release),
    }


def manifest_bytes(manifest: dict) -> bytes:
    return json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode()


def manifest_hash(manifest: dict) -> str:
    """BLAKE3-128 of the canonical manifest, excluding `root` so it is the same on every host."""
    portable = {k: v for k, v in manifest.items() if k != "root"}
    return blake3.blake3(manifest_bytes(portable)).hexdigest(length=16)


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def _refuse_record(root: Path, out: Path, rel: str, what: str) -> Path:
    """`out / rel`, after refusing (FileExistsError) when `root` or `out` already holds `rel`."""
    for base in dict.fromkeys([Path(root).resolve(), Path(out).resolve()]):
        if (base / rel).exists():
            raise FileExistsError(f"{base / rel}: {what} is written once; refusing to overwrite")
    return Path(out) / rel


def write_manifest(manifest: dict, path: Path | None = None, out: Path | None = None) -> Path:
    """Write `manifest` as sorted-key JSON under a temporary name renamed into place; return its path.

    With `path`: at that path. Without (card C-08): as the `release_manifest` artifact,
    `manifests/release=<R>.json` under `out` (default: the manifest's own `root`; under run-step
    `$BIOBOX_STAGING`, from which run-step publishes it), written once: FileExistsError, before
    anything is written, when the root or `out` already holds the release's manifest. The file holds
    the manifest's `root`, the data root it lists (`build_manifest`), wherever it is written."""
    if path is None:
        root = Path(manifest["root"])
        rel = manifest_path("", manifest["release_id"]).as_posix()
        path = _refuse_record(root, root if out is None else out, rel, f"the manifest of release {manifest['release_id']}")
    _atomic_write(Path(path), json.dumps(manifest, sort_keys=True, indent=1).encode() + b"\n")
    return Path(path)


def read_manifest(path: Path) -> dict:
    return json.loads(Path(path).read_text())


def _claimed_by(rel: str) -> tuple[str, dict[str, str]] | None:
    """The table or artifact that claims a root-relative path, with its partitions; None for a stray."""
    try:
        hit = classify(rel) or classify_artifact(rel)
    except LayoutError:
        return None
    return (hit[0].name, hit[1]) if hit else None


def verify_manifest(manifest: dict, root: Path | None = None) -> list[str]:
    """LY-10: every listed file exists with its sha256 and size, and no unlisted .parquet or artifact file
    exists in the release's scope (card B-33, as `build_manifest` draws it). A file of a stage the release
    does not list (card B-25), of another release's `release=` partition, or defining a tape version none
    of the listed files uses is not the release's and is ignored; any other stray is reported, also one no
    table claims, except in the trees `build_manifest` leaves out (`EXCLUDED_TREES`: the search pack and
    runs, written after the release, card B-28)."""
    root = Path(root or manifest["root"])
    used = _used_tape_versions(f["partitions"] for f in manifest["files"])
    in_release = _in_release(set(manifest["stages"]), manifest["release_id"], used)
    problems = []
    listed = set()
    for f in manifest["files"] + manifest.get("artifacts", []):
        listed.add(f["path"])
        p = root / f["path"]
        if not p.exists():
            problems.append(f"missing: {f['path']}")
            continue
        if p.stat().st_size != f["bytes"]:
            problems.append(f"size differs: {f['path']}")
        elif ids.sha256_file(p).hex() != f["sha256"]:
            problems.append(f"sha256 differs: {f['path']}")
    present = list(root.rglob("*.parquet"))
    for spec in ARTIFACTS.values():
        present += _artifact_files(root, spec)
    for p in sorted(set(present)):
        rel = p.relative_to(root).as_posix()
        if rel in listed or rel.split("/", 1)[0] in EXCLUDED_TREES:  # as build_manifest: runs carry their own records
            continue
        hit = _claimed_by(rel)
        if hit is not None and not in_release(*hit):
            continue
        problems.append(f"not in manifest: {rel}")
    return problems


def swap_current(manifests_dir: Path, release_id: str) -> Path:
    """Point `current.json` at `release=<R>.json` with an atomic rename of a fresh symlink."""
    d = Path(manifests_dir)
    target = f"release={release_id}.json"
    if not (d / target).exists():
        raise FileNotFoundError(d / target)
    tmp = d / f".current.json.tmp-{os.getpid()}"
    if tmp.is_symlink() or tmp.exists():
        tmp.unlink()
    os.symlink(target, tmp)
    os.replace(tmp, d / "current.json")
    return d / "current.json"


# --- catalog -----------------------------------------------------------------

_DUCK = {
    pa.string(): "VARCHAR", pa.bool_(): "BOOLEAN", pa.binary(): "BLOB",
    pa.int8(): "TINYINT", pa.int16(): "SMALLINT", pa.int32(): "INTEGER", pa.int64(): "BIGINT",
    pa.uint8(): "UTINYINT", pa.uint16(): "USMALLINT", pa.uint32(): "UINTEGER", pa.uint64(): "UBIGINT",
    pa.float32(): "FLOAT", pa.float64(): "DOUBLE",
}


def duck_type(t: pa.DataType) -> str:
    if pa.types.is_fixed_size_binary(t):
        return "BLOB"  # DuckDB has no fixed-size binary (LY-07)
    if pa.types.is_timestamp(t):
        return "TIMESTAMPTZ" if t.tz else "TIMESTAMP"
    if pa.types.is_list(t):
        return duck_type(t.value_type) + "[]"
    return _DUCK[t]


def _sql_str(s: str) -> str:
    return "'" + s.replace("'", "''") + "'"


def _view_sql(spec: TableSpec, files: list[tuple[str, tuple[str, ...]]], n_columns: int | None = None) -> str:
    """The view over `files`, (path, the partition keys of its path). One `read_parquet` per key set: a
    table with several (the ledger, card B-28) is their UNION ALL, each key a set lacks a typed NULL."""
    keys_all = all_keys(spec)
    if not files:
        cols = [f'CAST(NULL AS {duck_type(f.type)}) AS "{f.name}"' for f in spec.schema]
        cols += [f'CAST(NULL AS VARCHAR) AS "{k}"' for k in keys_all]
        return f"CREATE VIEW {spec.name} AS SELECT {', '.join(cols)} WHERE false"
    # union_by_name reads files written before a nullable column was added (LY-09); a column no
    # file has yet is a typed NULL
    n = len(spec.schema) if n_columns is None else n_columns
    data = [f'"{f.name}"' if i < n else f'CAST(NULL AS {duck_type(f.type)}) AS "{f.name}"'
            for i, f in enumerate(spec.schema)]
    groups: dict[tuple[str, ...], list[str]] = {}
    for path, keys in files:
        groups.setdefault(tuple(keys), []).append(path)
    selects = []
    for keys in sorted(groups, key=lambda ks: key_sets(spec).index(ks) if ks in key_sets(spec) else len(key_sets(spec))):
        listing = "[" + ", ".join(_sql_str(p) for p in groups[keys]) + "]"
        if keys:
            types = "{" + ", ".join(f"{_sql_str(k)}: 'VARCHAR'" for k in keys) + "}"
            opts = f"hive_partitioning = true, hive_types = {types}"
        else:
            opts = "hive_partitioning = false"
        cols = data + [f'"{k}"' if k in keys else f'CAST(NULL AS VARCHAR) AS "{k}"' for k in keys_all]
        selects.append(f"SELECT {', '.join(cols)} FROM read_parquet({listing}, {opts}, union_by_name = true)")
    return f"CREATE VIEW {spec.name} AS " + " UNION ALL ".join(selects)


_FIRST = "(l0_object_sha256, record_ordinal)"


def label_view_sql(files: list[str], n_stages: int) -> str:
    """The catalog's `contig_label` (card B-27; card W-08 made it a view over `contig_occurrence`): one row
    per contig of the release, the accession and header of its first source record by (l0_object_sha256,
    record_ordinal) over the release's stages, and its number of records (`n_occurrences`, as W-08 named it;
    also `n_records`). The columns are the registry table's, then `n_occurrences`; no `stage` column.

    The rows come from the registry table `contig_label` (one row per contig and stage, sorted by
    contig_uid), never from `contig_occurrence`, so a lookup by contig_uid prunes. With the files of one stage
    the view is a projection of them: a literal key and keys from a joined table both read only the row
    groups that hold them. With several stages it merges each contig's rows with a GROUP BY: a literal key
    still prunes, but DuckDB 2.0 pushes no join filter through an aggregate (16 of 16 row groups for one
    joined key, measured on wilkcomp02, 2026-11-10, planning/measurements/2026-11-10_label_and_membership_lookups.md),
    so there joined keys scan the label table; the exports filter the rows by key before they merge them
    (`biobox.search.export`) and prune in both cases."""
    cols = ("contig_uid, accession, original_header, n_records, l0_object_sha256, record_ordinal, "
            "n_records::BIGINT AS n_occurrences")
    if not files:
        # a release built before card B-27 has no contig_label files (card W-18): the same rows from
        # `contig_occurrence` (a view created before this one; empty and typed when the release has none),
        # one per contig, as the registry table would hold them. Not pruned by contig_uid: the slow path.
        return ("CREATE VIEW contig_label AS SELECT contig_uid, "
                f"arg_min(accession, {_FIRST}) AS accession, arg_min(original_header, {_FIRST}) AS original_header, "
                f"count(*)::UINTEGER AS n_records, arg_min(l0_object_sha256, {_FIRST}) AS l0_object_sha256, "
                f"arg_min(record_ordinal, {_FIRST}) AS record_ordinal, count(*)::BIGINT AS n_occurrences "
                "FROM contig_occurrence GROUP BY contig_uid")
    rows = f"read_parquet([{', '.join(_sql_str(p) for p in files)}], hive_partitioning = false)"
    if n_stages <= 1:
        return f"CREATE VIEW contig_label AS SELECT {cols} FROM {rows}"
    return ("CREATE VIEW contig_label AS SELECT contig_uid, "
            f"arg_min(accession, {_FIRST}) AS accession, arg_min(original_header, {_FIRST}) AS original_header, "
            f"sum(n_records)::UINTEGER AS n_records, arg_min(l0_object_sha256, {_FIRST}) AS l0_object_sha256, "
            f"arg_min(record_ordinal, {_FIRST}) AS record_ordinal, sum(n_records)::BIGINT AS n_occurrences "
            f"FROM {rows} GROUP BY contig_uid")


# DuckDB 2.0 turns the keys of a hash join's build side into a filter on the probe side's scan. Up to
# `dynamic_or_filter_threshold` keys (default 50) the filter is the list of keys, and a scan of a table sorted
# by the key reads only the row groups whose min/max statistics hold one of them; above it only the keys' min
# and max are pushed, and random keys then read every row group between them (59 of 1,024 row groups for 60
# random keys with the threshold raised, 1,015 at the default; 631 and 1,022 for 1,000 keys; measured on
# wilkcomp02, 2026-11-10, planning/measurements/2026-11-10_label_and_membership_lookups.md). Card B-27 raises it
# to the row groups a connection's largest table has: below that many keys the list prunes to the row groups
# that hold a key (at most two per key, AGENTS §4.4); above it the table has fewer row groups than there are
# keys, so a scan reads less than one row group per key, and the cheaper min/max filter is used. The list
# costs time on every row the scan reads (0.52 s against 0.16 s for 10,000 keys over 2^20 rows, measured
# as above), so it is capped at JOIN_KEY_FILTER_MAX keys. The cap holds the lookup tables at GTDB scale
# (derived: `member_cluster` 190 M rows / 122,880 = 1,547 row groups; `contig_label` 111.5 M contigs =
# 908); a table with more row groups (the postings) gets the list for up to that many keys and the min/max
# filter above it.
JOIN_KEY_FILTER_MIN = 50  # DuckDB's own default
JOIN_KEY_FILTER_MAX = 4096


def set_join_key_filter(con: duckdb.DuckDBPyConnection, row_groups: int) -> int:
    """Set the connection's `dynamic_or_filter_threshold` to `row_groups`, the row groups of the largest table
    its queries read (see above), at least DuckDB's default and at most JOIN_KEY_FILTER_MAX. Returns it."""
    n = max(JOIN_KEY_FILTER_MIN, min(JOIN_KEY_FILTER_MAX, int(row_groups)))
    con.execute(f"SET dynamic_or_filter_threshold = {n}")
    return n


def row_groups_of(files) -> int:
    """The row groups of Parquet `files`, from their footers."""
    return sum(pq.read_metadata(p).num_row_groups for p in files)


COMPLEMENT =("ACGTRYKMBVDHNacgtrykmbvdhn", "TGCAYRMKVBHDNtgcayrmkvbhdn")


def pair_macro_sql() -> list[str]:
    """Neighbour macros over `token_postings` (card B-11; planning/tasks/TASKS.md §2; AGENTS §3.5).

    Pairs are never stored: Q1 and Q4 self-join the postings of this catalog's release and tape
    version (`catalog_tape_version()`) on `contig_uid`. Each side of a lookup by token reads one
    `token_bucket` and the row groups whose `token_id` range holds the token, also when the token comes
    from `token_of` or a joined table (card B-31: `token_rows` joins the token in). Q3 (`partners_of`)
    takes its anchors from the postings and the tokens near them from `tape_tokens` of the anchors'
    contigs (`bb_partner_counts`), so its partner side reads row groups in proportion to the anchors'
    contigs, not to the release: exactly the `tape_tokens` row groups whose footer `contig_uid` range holds
    an anchor contig, at most min(C, R) for C anchor contigs and R row groups (10 anchor contigs: at most 10).
    That is the bound when the anchors are few. When they are many it is every row group: contig uids are
    random, so C contigs leave about R * exp(-C / R) of the R row groups untouched (*derived*, uniform keys),
    fewer than one once C > R ln R. The 10 Gbp rebuild read 26 of 26 for 345 anchors on 243 contigs
    (`rN.1b`) and 19 of 19 for an NR50 cluster of 2,139 slots (`rN.1a`) (*measured*, wilkcomp02,
    2026-11-10, `scripts/stress/build_stage_pod/r10_queries.py`); only a `tape_tokens` sorted otherwise
    could do better (card B-35). Each of those row groups is read once per query (card B-38). The geometry helpers
    reproduce `pairs_of_contig` (tests/oracles/pairs_oracle.py) and `biobox.l3.geometry`:

        pair_offset    signed shortest offset from slot _i to slot _j (j - i on a linear contig); on a
                       circle a partner exactly half way round takes the way with the smaller bp gap,
                       forward on a tie
        pair_crosses   the path of offset _d crosses the origin (circular only)
        pair_delta_bp  gap between the two slots ordered left to right, or minus their overlap, the right
                       slot unwrapped by one contig length when the path crosses the origin
        pair_junction  1 co-directional, 4 nested antisense, 2 divergent, 3 convergent
        pair_flags     1: the anchor's +-10-slot window is cut by a linear contig end (K_MAX, not the
                       query's k); 2: the path crosses the origin

    Slot arguments come as (ordinal, start, end[, strand]) of the anchor, then of the partner.
    """
    k_max, cens, cross = PAIR_K_MAX, PAIR_FLAG_CENSORED, PAIR_FLAG_CROSSES_ORIGIN
    left = "CASE WHEN _d > 0 THEN {a} ELSE {b} END"  # the slot to the left along the contig
    right = "(CASE WHEN _d > 0 THEN {b} ELSE {a} END + pair_shift(_d, _i, _j, _circ, _len))"
    ls, le, lstr = left.format(a="_as", b="_bs"), left.format(a="_ae", b="_be"), left.format(a="_astr", b="_bstr")
    rs, re_ = right.format(a="_as", b="_bs"), right.format(a="_ae", b="_be")
    fwd = "ring_fwd(_i, _j, _n)"
    delta = "pair_delta_bp({d}, _i, _as, _ae, _j, _bs, _be, _circ, _len)"
    return [
        "CREATE MACRO catalog_tape_version() AS (SELECT value FROM catalog_info WHERE key = 'tape_version')",
        # a token's bucket at every width a release may use (biobox.layout, biobox.l3.postings): 'all', or the
        # token's low byte as hex cut to the bucket's width. _tok appears once, so a key given as a subquery
        # (token_of) is evaluated once here
        "CREATE MACRO token_bucket_ok(_bucket, _tok) AS "
        "(_bucket = 'all' OR _bucket = left(printf('%02x', (_tok % 256)::INTEGER), length(_bucket)))",
        # every posting of a token in this catalog's release and tape version, and its slots (one row per
        # slot whatever the number of segments the token has there). The token is joined in (card B-31): a
        # literal, a scalar subquery (token_of) or a column of an outer query is one build-side row whose
        # value DuckDB pushes into the postings scan as a join filter, so the scan reads the row groups whose
        # token_id range holds it. `token_id = _tok` with a subquery read every row group (20 of 20 at 1,000x,
        # 2026-11-10 phase B exit); with a literal the bucket filter also prunes files, as before
        """CREATE MACRO token_rows(_tok) AS TABLE
              SELECT p.* FROM (SELECT (_tok)::UBIGINT AS _t) _k
              JOIN token_postings p ON p.token_id = _k._t
              WHERE p.release = catalog_release() AND p.tape_version = catalog_tape_version()
                AND token_bucket_ok(p.token_bucket, _tok)""",
        """CREATE MACRO token_slots(_tok) AS TABLE
              SELECT DISTINCT contig_uid, ordinal, feature_uid, start, "end", strand, n_slots, is_circular, contig_len
              FROM token_rows(_tok)""",
        "CREATE MACRO slot_dist(_i, _j, _n, _circ) AS CASE WHEN _circ "
        "THEN least(abs(_i::BIGINT - _j::BIGINT), _n::BIGINT - abs(_i::BIGINT - _j::BIGINT)) "
        "ELSE abs(_i::BIGINT - _j::BIGINT) END",
        "CREATE MACRO ring_fwd(_i, _j, _n) AS ((_j::BIGINT - _i::BIGINT) % _n::BIGINT + _n::BIGINT) % _n::BIGINT",
        "CREATE MACRO pair_crosses(_d, _i, _j, _circ) AS (_circ AND ((_d > 0 AND _j < _i) OR (_d < 0 AND _i < _j)))",
        "CREATE MACRO pair_shift(_d, _i, _j, _circ, _len) AS "
        "CASE WHEN pair_crosses(_d, _i, _j, _circ) THEN _len::BIGINT ELSE 0 END",
        f"CREATE MACRO pair_delta_bp(_d, _i, _as, _ae, _j, _bs, _be, _circ, _len) AS "
        f"(greatest({ls}, {rs}) - least({le}, {re_}))::BIGINT",
        f"""CREATE MACRO pair_junction(_d, _i, _as, _ae, _astr, _j, _bs, _be, _bstr, _circ, _len) AS CASE
              WHEN _astr = _bstr THEN 1
              WHEN ({ls} <= {rs} AND {re_} <= {le}) OR ({rs} <= {ls} AND {le} <= {re_}) THEN 4
              WHEN {lstr} = -1 THEN 2 ELSE 3 END""",
        f"""CREATE MACRO pair_offset(_i, _as, _ae, _j, _bs, _be, _n, _circ, _len) AS CASE
              WHEN NOT _circ THEN _j::BIGINT - _i::BIGINT
              WHEN 2 * {fwd} < _n THEN {fwd}
              WHEN 2 * {fwd} > _n THEN {fwd} - _n
              WHEN {delta.format(d=fwd + " - _n")} < {delta.format(d=fwd)} THEN {fwd} - _n
              ELSE {fwd} END""",
        f"""CREATE MACRO pair_flags(_d, _i, _j, _n, _circ) AS
              (CASE WHEN NOT _circ AND (_i::BIGINT < {k_max} OR _n::BIGINT - 1 - _i::BIGINT < {k_max}) THEN {cens} ELSE 0 END)
              + (CASE WHEN pair_crosses(_d, _i, _j, _circ) THEN {cross} ELSE 0 END)""",
        # Q1: loci (anchor slots) of A with B on another slot within k; _tier NULL for every tier, else
        # the contigs of that tier in this release's contig_dim (a contig without a row matches no tier)
        """CREATE MACRO q1_count(_a, _b, _k, _tier) AS (
              SELECT count(DISTINCT (a.contig_uid, a.ordinal))
              FROM token_slots(_a) a JOIN token_slots(_b) b ON b.contig_uid = a.contig_uid
              WHERE slot_dist(a.ordinal, b.ordinal, a.n_slots, a.is_circular) BETWEEN 1 AND _k
                AND (_tier IS NULL OR a.contig_uid IN (SELECT contig_uid FROM contig_dim_now() WHERE tier = _tier)))""",
        # Q4: every (A slot, B slot) within k, with the geometry of pairs_of_contig
        """CREATE MACRO pairs_of(_a, _b, _k) AS TABLE
              SELECT contig_uid, anchor_feature_uid, partner_feature_uid, anchor_ordinal, partner_ordinal,
                     orf_offset::SMALLINT AS orf_offset,
                     pair_junction(orf_offset, anchor_ordinal, a_s, a_e, a_str, partner_ordinal, b_s, b_e, b_str,
                                   c_circ, c_len)::UTINYINT AS junction_topology,
                     pair_delta_bp(orf_offset, anchor_ordinal, a_s, a_e, partner_ordinal, b_s, b_e, c_circ,
                                   c_len)::INTEGER AS delta_bp,
                     pair_flags(orf_offset, anchor_ordinal, partner_ordinal, c_n, c_circ)::UTINYINT AS flags
              FROM (SELECT a.contig_uid, a.feature_uid AS anchor_feature_uid, b.feature_uid AS partner_feature_uid,
                           a.ordinal AS anchor_ordinal, b.ordinal AS partner_ordinal,
                           a.start AS a_s, a."end" AS a_e, a.strand AS a_str,
                           b.start AS b_s, b."end" AS b_e, b.strand AS b_str,
                           a.n_slots AS c_n, a.is_circular AS c_circ, a.contig_len AS c_len,
                           pair_offset(a.ordinal, a.start, a."end", b.ordinal, b.start, b."end", a.n_slots,
                                       a.is_circular, a.contig_len) AS orf_offset
                    FROM token_slots(_a) a JOIN token_slots(_b) b ON b.contig_uid = a.contig_uid
                    WHERE slot_dist(a.ordinal, b.ordinal, a.n_slots, a.is_circular) BETWEEN 1 AND _k)
              ORDER BY contig_uid, anchor_ordinal, partner_ordinal""",
        # both tokens on one slot with different segments (domain postings, seg_idx > 0)
        """CREATE MACRO fusions_of(_a, _b) AS TABLE
              SELECT a.contig_uid, a.feature_uid, a.ordinal, a.seg_idx AS anchor_seg_idx, b.seg_idx AS partner_seg_idx
              FROM token_rows(_a) a JOIN token_rows(_b) b
                ON b.contig_uid = a.contig_uid AND b.ordinal = a.ordinal AND b.seg_idx <> a.seg_idx
              ORDER BY a.contig_uid, a.ordinal, anchor_seg_idx, partner_seg_idx""",
        partners_of_sql(),
    ]


def partners_of_sql(create: str = "CREATE MACRO") -> str:
    """Q3 `partners_of(_a, _k)` (counts; z_score ranking arrives with pair_stats): every token at every level
    within k slots of the loci of A. `open_catalog` creates it again as a temporary macro (`create`), so a
    catalog built before card B-38, whose stored macro calls card B-31's `bb_slot_tokens` per contig, runs
    this one.

    Anchors come from the postings of A (by token_id); the partner slots are the slots within k of each
    anchor on its contig (slot_dist 1..k, each slot once per anchor), and their tokens come from `tape_tokens`
    of the anchors' contigs through a point lookup by contig_uid (card B-31) that reads only the row groups
    holding those contigs: at most one per anchor contig, and every row group once the anchor contigs
    outnumber the row groups by about ln(row groups) (see pair_macro_sql; card B-35). Every partner slot of
    the query goes to `bb_partner_counts` in one call, as one list of (contig, anchor slot, partner slot)
    structs (`_all`), so each of those row groups is read once (card B-38), and the function counts per token
    the loci (distinct (contig, anchor slot)) and the pairs (partner slots that hold it). Postings are sorted
    by token_id, so the earlier near side (postings filtered by contig_uid) read the whole token level of the
    release. A hub is still refused here (its anchors are many).

    No switch to a plain SQL join (card B-38): neither join prunes, and the lookup wins or ties at every anchor
    count measured on `rN.1` (10 Gbp, 93 `tape_tokens` row groups, k = 3; seconds, median of 3 runs on a fresh
    connection, 4 threads; *measured*, wilkcomp02, 2026-11-10, card B-38 probe `b38_rn1.py`, units b38-base and
    b38-new; tape join = the anchor windows joined to `tape_tokens`, the same answer from the same table, 93 of
    93 row groups every time; postings self-join = the ladder's form, 197 of 197 `token_postings` row groups):

        anchors (contigs)   before B-38            after          tape join   postings self-join
        10 (8)              0.33 (8 reads)         0.14 (8)       0.47        0.39
        100 (29)            0.84 (27)              0.23 (27)      0.49        0.44
        300 (291)           4.61 (191 of 92)       0.47 (92)      0.45        0.40
        989 (645)           7.57 (316 of 93)       0.47 (93)      0.44        0.44
        2,950 (2,331)       13.24 (588 of 93)      0.54 (93)      0.76        2.38
        3,509 (2,399)       12.35 (543 of 93)      0.52 (93)      0.78        2.40

    Every answer equal (sha256 of the rows). The tape join was 3 and 6% faster at 291 and 645 anchor contigs,
    where both read every row group: a tie within the timing tolerance of the B-38 test (1.25x), and no count
    above which a join stays faster."""
    return f"""{create} partners_of(_a, _k) AS TABLE
              WITH _anchor AS (
                     SELECT * FROM token_slots(_a)
                     WHERE (SELECT CASE WHEN count(*) > {PARTNERS_MAX_ROWS} THEN error('partners_of: token has more '
                                   || 'than {PARTNERS_MAX_ROWS} loci; use pair_stats (hub routing)') ELSE true END
                            FROM token_rows(_a))),
                   _off AS (
                     SELECT contig_uid, ordinal::BIGINT AS _i, n_slots::BIGINT AS _n, is_circular AS _c,
                            unnest(range(-least(_k::BIGINT, n_slots::BIGINT), least(_k::BIGINT, n_slots::BIGINT) + 1)) AS _d
                     FROM _anchor),
                   _win AS (
                     SELECT DISTINCT contig_uid, _i, _n,
                            CASE WHEN _c THEN ((_i + _d) % _n + _n) % _n ELSE _i + _d END AS _j
                     FROM _off WHERE _d <> 0),
                   _all AS (
                     SELECT list({{'c': contig_uid, 'i': _i, 'j': _j}}) AS _w
                     FROM _win WHERE _j >= 0 AND _j < _n AND _j <> _i),
                   _n AS (SELECT unnest(bb_partner_counts(_w)) AS _r FROM _all)
              SELECT _r.t AS partner_token, _r.l AS n_loci, _r.p AS n_pairs
              FROM _n ORDER BY n_loci DESC, partner_token"""


def macro_sql() -> list[str]:
    """Sequence macros (data_layer_and_tests.md §2A.3). 0-based half-open; strand +1/-1.

    Parameter names start with an underscore so they never collide with column names. The
    sequence macros only call `bb_seq`, `bb_seq_circular` and `bb_seq_raw` (registered by
    `open_catalog`), which check bounds and strand, resolve keys and read the shards; the lookup
    macros `default_digest`, `digest_len` and `contig_len_of` are for users. As stored here they scan
    the index views (for a plain `duckdb` connection); `open_catalog` replaces them on its connection
    with point lookups (`_register_lookup_udf`, card B-31). The sequence macros do not call them.

    The first macros are the stand-ins of `PY_FUNCTIONS` (card B-34): on a plain connection each raises
    an error that names `biobox.catalog.open_catalog`. The error is cast to the function's return type:
    an untyped `error()` passed to another function is folded to NULL by DuckDB's binder (2.0.0), and the
    caller would get NULL instead of the error.
    """
    stand_ins = [
        f"CREATE MACRO {name}({params}) AS CAST(error('{users}: open this catalog with "
        f"biobox.catalog.open_catalog(path); {name} is a Python function of biobox that a plain duckdb "
        f"connection does not have') AS {rtype})"
        for name, (params, rtype, users) in PY_FUNCTIONS.items()]
    return stand_ins + [
        f"CREATE MACRO revcomp(_x) AS reverse(translate(_x, '{COMPLEMENT[0]}', '{COMPLEMENT[1]}'))",
        "CREATE MACRO default_digest(_cuid) AS "
        "(SELECT raw_digest FROM contig_index WHERE contig_uid = _cuid AND is_default_variant)",
        "CREATE MACRO digest_len(_rd) AS (SELECT contig_len FROM seq_locator WHERE raw_digest = _rd)",
        "CREATE MACRO contig_len_of(_cuid) AS "
        "(SELECT contig_len FROM contig_index WHERE contig_uid = _cuid AND is_default_variant)",
        # sequence: thin wrappers over the Python functions (card B-14; see _register_seq_udf).
        # forward bytes of one stored variant, by exact raw_digest
        "CREATE MACRO seq_raw(_rd, _s, _e) AS bb_seq_raw(_rd, _s, _e)",
        # the default variant of a contig, forward strand
        "CREATE MACRO seq_slice(_cuid, _s, _e) AS bb_seq(_cuid, _s, _e, 1)",
        "CREATE MACRO seq_get(_cuid, _s, _e, _strand) AS bb_seq(_cuid, _s, _e, _strand)",
        # circular contigs: an origin-wrapping slice has contig_len < _e <= _s + contig_len (C-05, DEC-13)
        "CREATE MACRO seq_get_circular(_cuid, _s, _e, _strand) AS bb_seq_circular(_cuid, _s, _e, _strand)",
        "CREATE MACRO protein_of(_puid) AS (SELECT any_value(seq) FROM protein_seq "
        # the key's bucket at every width a stage may use (biobox.layout); prunes files
        "WHERE protein_prefix IN ('all', left(lower(hex((_puid)[1:1])), 1), lower(hex((_puid)[1:1]))) "
        "AND protein_uid = _puid)",
        "CREATE MACRO slots_of(_puid) AS TABLE SELECT * FROM protein_occurrence WHERE protein_uid = _puid",
        "CREATE MACRO domains_of(_puid) AS TABLE SELECT * FROM l1_domains WHERE protein_uid = _puid ORDER BY env_start, model_acc",
        # --- L3 (Q1-Q4); release-scoped tables are read for this catalog's release only ---
        "CREATE MACRO catalog_release() AS (SELECT value FROM catalog_info WHERE key = 'release_id')",
        "CREATE MACRO token_of(_level, _key) AS "
        "(SELECT token_id FROM token_dict WHERE release = catalog_release() AND level = _level AND key = _key)",
        # dimension rows of this catalog's release only (card B-12)
        "CREATE MACRO contig_dim_now() AS TABLE SELECT * FROM contig_dim WHERE release = catalog_release()",
        *pair_macro_sql(),
        # Q2: tapes of a list of contigs (every tape version) with this release's tokens. Each side is
        # filtered by `IN` over the unnested list before the join, which prunes row groups by contig_uid
        # in both tables, also when the list comes from a joined table (AGENTS §4.4, card B-07). In DuckDB
        # 2.0.0 `list_contains` prunes only a literal list, and the same `IN` written after the join
        # scanned every tapes row group on the micro fixture
        """CREATE MACRO tapes_of(_cuids) AS TABLE
              WITH _q AS (SELECT DISTINCT unnest(_cuids) AS _cu)
              SELECT t.*, k.family_tokens, k.nr50_tokens, k.nr90_tokens
              FROM (SELECT * FROM contig_tapes WHERE contig_uid IN (SELECT _cu FROM _q)) t
              LEFT JOIN (SELECT * FROM tape_tokens
                         WHERE release = catalog_release() AND contig_uid IN (SELECT _cu FROM _q)) k
                ON k.contig_uid = t.contig_uid AND k.tape_version = t.tape_version
              ORDER BY t.contig_uid, t.tape_version""",
    ]


def served_tape_version(manifest: dict) -> str:
    """The one tape version of the release's `token_postings` in the manifest; "" when it has none.

    Raises ValueError when the release holds postings of several tape versions: the caller must name one.
    """
    rel = manifest["release_id"]
    tvs = sorted({f["partitions"]["tape_version"] for f in manifest["files"]
                  if f["table"] == "token_postings" and f["partitions"].get("release") == rel})
    if len(tvs) > 1:
        raise ValueError(f"token_postings of release {rel} hold {len(tvs)} tape versions {tvs}; "
                         "pass tape_version to generate_catalog")
    return tvs[0] if tvs else ""


def catalog_path(root: Path, release: str) -> Path:
    """Where a release's catalog is published: `catalog/release=<R>.duckdb` under the data root."""
    return Path(root) / f"catalog/release={release}.duckdb"


def generate_catalog(manifest: dict, out_path: Path | None = None, root: Path | None = None,
                     tape_version: str | None = None, out: Path | None = None) -> Path:
    """Write the catalog under a temporary name and rename it into place (atomic on one filesystem).

    `tape_version`: the tape version the pair macros serve (`catalog_tape_version()`); by default the
    single one of the release's `token_postings` (`served_tape_version`).

    The catalog holds absolute paths: `catalog_info.root` and every file in its views are under `root`
    (default: the manifest's `root`), the data root its readers open. Where the file goes is separate
    (card C-08): `out_path` when given; else the `catalog` artifact, `catalog/release=<R>.duckdb` under
    `out` (default `root`; under run-step `$BIOBOX_STAGING`, so a catalog built in staging still names
    the final root), written once: FileExistsError, before any work, when the root or `out` holds it.
    """
    root = Path(root or manifest["root"]).resolve()
    if out_path is None:
        rel = catalog_path("", manifest["release_id"]).as_posix()
        out_path = _refuse_record(root, root if out is None else out, rel, f"the catalog of release {manifest['release_id']}")
    if tape_version is None:
        tape_version = served_tape_version(manifest)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(f".{out_path.name}.tmp-{os.getpid()}")
    for p in (tmp, Path(str(tmp) + ".wal")):
        if p.exists():
            p.unlink()
    by_table: dict[str, list[tuple[str, tuple[str, ...]]]] = {name: [] for name in TABLES}
    for f in manifest["files"]:
        by_table[f["table"]].append((str(root / f["path"]), tuple(f["partitions"])))
    con = duckdb.connect(str(tmp))
    try:
        widest = {}
        for f in manifest["files"]:
            widest[f["table"]] = max(widest.get(f["table"], 0), f.get("n_columns", len(TABLES[f["table"]].schema)))
        for name in sorted(TABLES):
            if name != "contig_label":  # the release's labels, merged over its stages (label_view_sql)
                con.execute(_view_sql(TABLES[name], sorted(by_table[name]), widest.get(name)))
        label_stages = {f["partitions"].get("stage") for f in manifest["files"] if f["table"] == "contig_label"}
        con.execute(label_view_sql([p for p, _ in sorted(by_table["contig_label"])], len(label_stages)))
        con.execute("CREATE TABLE catalog_info (key VARCHAR, value VARCHAR)")
        con.executemany("INSERT INTO catalog_info VALUES (?, ?)", [
            ("release_id", manifest["release_id"]),
            ("manifest_hash", manifest_hash(manifest)),
            ("root", str(root)),
            ("duckdb_version", duckdb.__version__),
            ("tape_version", tape_version),
            ("stages", json.dumps(manifest.get("stages", []))),  # the release's stages, in order (card B-28)
        ])
        con.execute("CREATE TABLE catalog_files AS SELECT * FROM (VALUES (NULL::VARCHAR, NULL::VARCHAR, NULL::VARCHAR, NULL::BIGINT, NULL::BIGINT)) "
                    "t(path, \"table\", sha256, \"rows\", row_groups) WHERE false")
        rows = [(f["path"], f["table"], f["sha256"], f["rows"], f["row_groups"]) for f in manifest["files"]]
        if rows:  # a release with no Parquet file of its own (card B-33): executemany refuses an empty list
            con.executemany("INSERT INTO catalog_files VALUES (?, ?, ?, ?, ?)", rows)
        con.execute("CREATE TABLE catalog_artifacts (path VARCHAR, artifact VARCHAR, sha256 VARCHAR, bytes BIGINT)")
        arts = manifest.get("artifacts", [])
        if arts:
            con.executemany("INSERT INTO catalog_artifacts VALUES (?, ?, ?, ?)",
                            [(a["path"], a["artifact"], a["sha256"], a["bytes"]) for a in arts])
        # DuckDB binds a macro's functions when it is created: the stand-in macros of PY_FUNCTIONS come
        # first in macro_sql and are stored (card B-34; open_catalog shadows them with the real readers);
        # macros come last because they read catalog_info
        for stmt in macro_sql():
            con.execute(stmt)
        con.execute("CHECKPOINT")
    finally:
        con.close()
    os.replace(tmp, out_path)
    return out_path


_READERS: dict[int, object] = {}
# id(connection) -> its point lookups outside DuckDB's scans, as (function name, _Lazy): the tape tokens of
# bb_partner_counts (cards B-31, B-38) and the token keys of token_of (card B-34); `profile_query` reports their reads
_LOOKUPS: dict[int, list] = {}


def _shadow(con: duckdb.DuckDBPyConnection, name: str) -> None:
    """Shadow the stored stand-in macro `name` of PY_FUNCTIONS, on this connection only, with a temporary
    macro of the same name over the registered Python function `<name>_udf` (card B-34). Temporary objects
    come first in DuckDB's search path, also inside the stored macros that call `name`."""
    params = PY_FUNCTIONS[name][0]
    con.execute(f"CREATE OR REPLACE TEMP MACRO {name}({params}) AS {name}_udf({params})")


def _read_in_file_order(store, n_rows: int, jobs: list) -> list:
    """Read one vector's slices and return them as strings in row order (None for NULL rows).

    `jobs` holds (row, strand, pieces), pieces being (shard, raw_digest, start, end) already
    bounds-checked; an origin-wrapping slice is two pieces. The pieces are read sorted by
    (shard, raw_digest, start, end), which is (shard, file offset) order because a shard's records
    are sorted by raw_digest (the bisected `.fai` requires it): a batch of fetches decompresses
    each BGZF block once and front to back, whatever order the rows came in. The results are then
    put back in the rows' original order.
    """
    from biobox.l0seq.store import revcomp

    order = sorted(((p, row, k) for row, _, pieces in jobs for k, p in enumerate(pieces)), key=lambda x: x[0])
    got = {}
    for (shard, rd, a, b), row, k in order:
        got[(row, k)] = store.read(shard, rd, a, b) if a < b else b""
    out: list = [None] * n_rows
    for row, strand, pieces in jobs:
        x = b"".join(got[(row, k)] for k in range(len(pieces)))
        out[row] = (x if strand == 1 else revcomp(x)).decode("ascii")
    return out


def _register_seq_udf(con: duckdb.DuckDBPyConnection, store) -> None:
    """Register `bb_seq`, `bb_seq_circular` and `bb_seq_raw` over `store` (card B-14).

    Each is vectorized (Arrow in, Arrow out) and resolves keys with SeqStore point lookups
    (`default_digest`, `locate`: one row group per index per key, cached), never with SQL on the
    index views. A NULL in any argument gives NULL. With `store` None (a release without
    `seq_locator` files) they raise. Registered as `<name>_udf` and shadowing the stored stand-ins
    (`_shadow`, card B-34).
    """

    def need() -> None:
        if store is None:
            raise RuntimeError("bb_seq: this catalog's release lists no seq_locator files, so it has no sequence to read")

    def contig_jobs(who: str, circular: bool, cuids: pa.Array, starts: pa.Array, ends: pa.Array,
                    strands: pa.Array) -> pa.Array:
        need()
        rows = list(zip(cuids.to_pylist(), starts.to_pylist(), ends.to_pylist(), strands.to_pylist()))
        resolved: dict[bytes, tuple[bytes, str, int]] = {}  # contig_uid -> (raw_digest, shard, contig_len)
        jobs = []
        for i, (cu, s, e, t) in enumerate(rows):
            if cu is None or s is None or e is None or t is None:
                continue
            if t not in (1, -1):
                raise ValueError(f"{who}: strand must be 1 or -1 (got {t})")
            if cu not in resolved:
                try:
                    rd = store.default_digest(cu)
                except KeyError:
                    raise ValueError(f"{who}: unknown contig_uid {cu.hex()} (in none of the stages "
                                     f"{store.stages} this catalog's sequence store reads)") from None
                resolved[cu] = (rd, *store.locate(rd))
            rd, shard, n = resolved[cu]
            if circular and e > n:
                if not 0 <= s < n < e <= s + n:
                    raise ValueError(f"{who}: invalid wrapping slice {s}:{e} (contig_len {n})")
                pieces = [(shard, rd, s, n), (shard, rd, 0, e - n)]
            elif not 0 <= s <= e <= n:
                raise ValueError(f"{who}: slice outside 0:contig_len ({s}:{e}, contig_len {n})")
            else:
                pieces = [(shard, rd, s, e)]
            jobs.append((i, t, pieces))
        return pa.array(_read_in_file_order(store, len(rows), jobs), pa.string())

    def bb_seq(cuids, starts, ends, strands):
        return contig_jobs("seq_get", False, cuids, starts, ends, strands)

    def bb_seq_circular(cuids, starts, ends, strands):
        return contig_jobs("seq_get_circular", True, cuids, starts, ends, strands)

    def bb_seq_raw(rds: pa.Array, starts: pa.Array, ends: pa.Array) -> pa.Array:
        need()
        rows = list(zip(rds.to_pylist(), starts.to_pylist(), ends.to_pylist()))
        located: dict[bytes, tuple[str, int]] = {}
        jobs = []
        for i, (rd, s, e) in enumerate(rows):
            if rd is None or s is None or e is None:
                continue
            if rd not in located:
                try:
                    located[rd] = store.locate(rd)
                except KeyError:
                    raise ValueError(f"seq_raw: unknown raw_digest {rd.hex()} (in none of the stages "
                                     f"{store.stages} this catalog's sequence store reads)") from None
            shard, n = located[rd]
            if not 0 <= s <= e <= n:
                raise ValueError(f"seq_raw: slice outside 0:contig_len ({s}:{e}, contig_len {n})")
            jobs.append((i, 1, [(shard, rd, s, e)]))
        return pa.array(_read_in_file_order(store, len(rows), jobs), pa.string())

    opts = {"type": "arrow", "null_handling": "special", "side_effects": False}
    con.create_function("bb_seq_udf", bb_seq, ["BLOB", "BIGINT", "BIGINT", "BIGINT"], "VARCHAR", **opts)
    con.create_function("bb_seq_circular_udf", bb_seq_circular, ["BLOB", "BIGINT", "BIGINT", "BIGINT"], "VARCHAR", **opts)
    con.create_function("bb_seq_raw_udf", bb_seq_raw, ["BLOB", "BIGINT", "BIGINT"], "VARCHAR", **opts)
    for name in ("bb_seq", "bb_seq_circular", "bb_seq_raw"):
        _shadow(con, name)


def _register_lookup_udf(con: duckdb.DuckDBPyConnection, store) -> None:
    """`default_digest`, `digest_len` and `contig_len_of` as point lookups over `store` (card B-31).

    The macros the catalog stores scan `contig_index` or `seq_locator`: a key from a joined table is a
    correlated subquery, which DuckDB 2.0 answers with a left join that pushes no filter into the index
    scan (10 of 10 row groups for one joined key at 100x; every rewrite as a subquery read the same,
    measured on wilkcomp02, 2026-11-10). `open_catalog` shadows them, on its connection only, with temporary
    macros over `bb_default_digest`, `bb_digest_len` and `bb_contig_len_of`, which resolve each key with
    SeqStore's `KeyLookup`s (one row group per index per key, counted in `profile_query`'s `index`), as the
    sequence functions do (card B-14). Unknown keys and NULLs give NULL, as the SQL macros do. A catalog
    opened with plain `duckdb` keeps the SQL macros: the same answers, without the pruning.
    """

    def per_key(fn, rtype):
        def udf(keys: pa.Array) -> pa.Array:
            got: dict[bytes, object] = {}
            out = []
            for k in keys.to_pylist():
                if k is None:
                    out.append(None)
                    continue
                if k not in got:
                    try:
                        got[k] = fn(k)
                    except KeyError:
                        got[k] = None
                out.append(got[k])
            return pa.array(out, rtype)
        return udf

    opts = {"type": "arrow", "null_handling": "special", "side_effects": False}
    con.create_function("bb_default_digest", per_key(store.default_digest, pa.binary()), ["BLOB"], "BLOB", **opts)
    con.create_function("bb_digest_len", per_key(store.length_of, pa.int64()), ["BLOB"], "BIGINT", **opts)
    con.create_function("bb_contig_len_of", per_key(lambda cu: store.length_of(store.default_digest(cu)), pa.int64()),
                        ["BLOB"], "BIGINT", **opts)
    for stmt in ("CREATE OR REPLACE TEMP MACRO default_digest(_cuid) AS bb_default_digest(_cuid)",
                 "CREATE OR REPLACE TEMP MACRO digest_len(_rd) AS bb_digest_len(_rd)",
                 "CREATE OR REPLACE TEMP MACRO contig_len_of(_cuid) AS bb_contig_len_of(_cuid)"):
        con.execute(stmt)


TAPE_TOKEN_LISTS = ("family_tokens", "nr50_tokens", "nr90_tokens")  # a slot's levels, in this order
# threads that read tape_tokens row groups for one partners_of query (card B-38): the connection's DuckDB threads,
# at most this many; a cap that bounds the row groups in flight (TapeTokenLookup), not a measured optimum
TAPE_READ_THREADS_MAX = 8
# a batch of tape_tokens row groups read and decoded together (card B-38): caps that bound memory, not measured optima
TAPE_BATCH_ROW_GROUPS = 16
TAPE_BATCH_BYTES = 4 << 20


class _Touches:
    """The row groups a point lookup touched, for `profile_query` (card B-35).

    Every access to a row group, decoded from its file or served by the lookup's row-group cache, stamps it
    with a running count (`stamp`); `touched_since(mark)` gives the row groups accessed since `mark` (an earlier
    value of `stamp`) with their files, and those of them decoded in that span. A cache hit is a read the
    query made: counting decodes only, a query whose row groups an earlier query on the same connection had
    cached reported no read at all (10 Gbp rebuild, 2026-11-10). Holds one entry per row group touched, so
    at most the table's row groups."""

    def _touches_init(self) -> None:
        self.stamp = 0
        self._last_touch: dict = {}  # row group -> (stamp of its last access, file index)
        self._last_decode: dict = {}  # row group -> stamp of its last decode

    def _touch(self, rg, fi: int, decoded: bool) -> None:
        self._last_touch[rg] = (self.stamp, fi)
        if decoded:
            self._last_decode[rg] = self.stamp
        self.stamp += 1

    def touched_since(self, mark: int) -> tuple[dict, set]:
        """({row group: file index} accessed at or after `mark`, the set of those decoded at or after it)."""
        with self._lock:
            touched = {rg: fi for rg, (s, fi) in self._last_touch.items() if s >= mark}
            decoded = {rg for rg, s in self._last_decode.items() if s >= mark}
        return touched, decoded


class TapeTokenLookup(_Touches, KeyLookup):
    """`tape_tokens` rows by `contig_uid` from the release's files of one tape version (card B-31).

    The footer statistics of `KeyLookup` (a contig is looked for only in the row groups whose `contig_uid`
    range holds it), with its own decoding: a row group is held as numpy arrays per token list (offsets,
    values with nulls as 0, validity) and its keys as one Arrow array, so finding slots' tokens is array
    indexing, never a Python object per slot or per contig. `slot_tokens` answers a whole query in one call
    (card B-38): it visits the row groups that may hold the query's contigs once each, in file and row-group
    order, whatever order the contigs came in, so no row group is decoded twice in one call. Before card B-38
    the lookup was called per contig in DuckDB's order and kept `cache_row_groups` (32) row groups: with more
    row groups than that, a query decoded some again, 556 decodes of 93 row groups for 3,509 anchors on `rN.1`
    in 9.5 s (*measured*, wilkcomp02, 2026-11-10, card B-35). `get` of the base class is not available.

    Row groups the cache does not hold are read in batches of consecutive row groups of one file (one Parquet
    read and one numpy decoding per batch, at most TAPE_BATCH_ROW_GROUPS row groups and, past the first, at most
    TAPE_BATCH_BYTES decoded bytes as the footer counts them), by `threads` threads (`open_catalog` passes the
    connection's DuckDB thread count, at most TAPE_READ_THREADS_MAX), at most `threads` batches ahead of the one
    in use. A batch is cached as one entry and counts as its row groups. On small row groups the time was the
    fixed cost per row group in Python, about 1.9 ms each (0.23 s for 119 row groups of at most 1,024 slots,
    *measured*, wilkcomp02, 2026-11-10, cProfile of the B-38 test store), which batching pays once per batch.

    Memory (AGENTS §4.1): the LRU of at most `cache_row_groups` decoded row groups, which lives as long as the
    connection and serves repeated queries, plus at most `threads` + 1 batches being read or used. A decoded row
    group holds 9 bytes per slot and level (value, validity) and 28 bytes per contig (key, three int32 offsets):
    at most 131,072 x (3 x 9 + 28) bytes, about 7.2 MB, for a row group at the slot cap of AGENTS §4.6 with one
    slot per contig, and about 3.5 MB with long contigs (*derived*). A batch holds at most max(TAPE_BATCH_BYTES,
    one row group), so at most 32 x 7.2 + 9 x 7.2 MB, about 300 MB, with the default cache and the most threads
    (*derived*). The query's own arrays are 8 bytes per slot asked (in), 27 bytes per slot (out) and 16 per
    contig, which the macro bounds by its hub refusal (PARTNERS_MAX_ROWS anchors, each with at most 2k slots).

    `row_groups_read` and `files_hit` (per file) count decodes (a cache hit counts nothing); `touched_since`
    (`_Touches`) also counts the row groups the cache served.
    """

    def __init__(self, files, cache_row_groups: int = 32, threads: int = 1):
        super().__init__(files, "contig_uid", TAPE_TOKEN_LISTS, cache_row_groups)
        self._touches_init()
        self.threads = max(1, min(int(threads), TAPE_READ_THREADS_MAX))
        self.files_hit = [0] * len(self._files)
        self.row_groups_total = sum(len(f.row_groups) for f in self._files)
        self._units: OrderedDict = OrderedDict()  # cached batches: tuple of (file, row group) -> decoded
        self._unit_of: dict = {}  # (file, row group) -> the cached batch that holds it
        self._held = 0  # row groups the cached batches hold
        self._bytes: dict = {}  # (file, row group) -> its decoded bytes, from the footer (see the class docstring)
        for fi, f in enumerate(self._files):
            md = f.metadata
            lists = [i for i in range(md.num_columns) if md.schema.column(i).path.split(".")[0] in TAPE_TOKEN_LISTS]
            for g in f.row_groups:
                rg = md.row_group(g)
                self._bytes[(fi, g)] = 9 * sum(rg.column(i).num_values for i in lists) + 28 * rg.num_rows

    def _batches(self, locs: list) -> list[tuple]:
        """`locs` (sorted) cut into batches: one file, at most TAPE_BATCH_ROW_GROUPS row groups and, past the first,
        at most TAPE_BATCH_BYTES."""
        out, cur, size = [], [], 0
        for loc in locs:
            b = self._bytes[loc]
            if cur and (loc[0] != cur[0][0] or len(cur) == TAPE_BATCH_ROW_GROUPS or size + b > TAPE_BATCH_BYTES):
                out.append(tuple(cur))
                cur, size = [], 0
            cur.append(loc)
            size += b
        if cur:
            out.append(tuple(cur))
        return out

    def _read(self, batch: tuple):
        """The row groups of `batch` (one file) decoded together: ([(offsets, values, valid or None) per level],
        keys). Touches no state of the lookup, so several threads may run it."""
        fi = batch[0][0]
        f = self._files[fi]
        t = pq.ParquetFile(f.path, metadata=f.metadata).read_row_groups([g for _, g in batch], columns=self.columns,
                                                                         use_threads=False)
        levels = []
        for c in TAPE_TOKEN_LISTS:
            a = t.column(c).combine_chunks()
            v = a.values  # row groups read whole and combined: offsets start at 0 and index `values`
            vals = (np.frombuffer(v.buffers()[1], np.uint64, len(v), v.offset * 8) if len(v)
                    else np.zeros(0, np.uint64))  # a null token's value is never used (valid is False)
            valid = None if v.null_count == 0 else v.is_valid().to_numpy(zero_copy_only=False)
            levels.append((a.offsets.to_numpy(zero_copy_only=False), vals, valid))
        return levels, t.column(self.key_column).combine_chunks().cast(pa.binary())

    def _decoded_from(self, batch: tuple, out):
        """Count the decode of each row group of `batch` and keep the batch in the cache."""
        for loc in batch:
            self.row_groups_read += 1
            self.files_hit[loc[0]] += 1
            self._touch(loc, loc[0], True)
        if 0 < len(batch) <= self.cache_row_groups:
            self._units[batch] = out
            self._unit_of.update(dict.fromkeys(batch, batch))
            self._held += len(batch)
            while self._held > self.cache_row_groups:
                old, _ = self._units.popitem(last=False)
                for loc in old:
                    del self._unit_of[loc]
                self._held -= len(old)
        return out

    def _visit(self, locs: list, use) -> None:
        """Call `use(decoded, row_groups)` so that every row group of `locs` is in exactly one call: first the
        cached batches that hold some of them, then the others in batches (`_batches`), read by `threads`
        threads at most `threads` batches ahead (the memory bound above)."""
        hits: dict[tuple, list] = {}
        misses = []
        for loc in locs:
            unit = self._unit_of.get(loc)
            if unit is None:
                misses.append(loc)
            else:
                hits.setdefault(unit, []).append(loc)
        for unit, ls in hits.items():
            self._units.move_to_end(unit)
            for loc in ls:
                self._touch(loc, loc[0], False)
            use(self._units[unit], ls)
        batches = self._batches(misses)
        if self.threads == 1 or len(batches) < 2:
            for b in batches:
                use(self._decoded_from(b, self._read(b)), b)
            return
        from collections import deque
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(self.threads) as ex:
            todo, ahead = iter(batches), deque()
            for b in todo:
                ahead.append((b, ex.submit(self._read, b)))
                if len(ahead) == self.threads:
                    break
            while ahead:
                b, fut = ahead.popleft()
                got = fut.result()
                nxt = next(todo, None)
                if nxt is not None:
                    ahead.append((nxt, ex.submit(self._read, nxt)))
                use(self._decoded_from(b, got), b)

    def get(self, key):  # noqa: D102
        raise NotImplementedError("TapeTokenLookup: use slot_tokens()")

    def slot_tokens(self, keys: list, lens: np.ndarray, js: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """The tokens of slots of contigs: contig c (`keys[c]`, bytes or None) owns the next `lens[c]` entries of
        `js`, 0-based slots of its tape. Returns (tok, ok), each of shape (len(js), 3), levels in TAPE_TOKEN_LISTS
        order: tok[s, L] is slot s's level-L token where ok[s, L]; a slot outside its tape, a null token and a
        contig no row group holds give none. Each row group whose footer range holds one of the keys is visited
        once (`_visit`): served by the cache, or decoded."""
        lens = np.asarray(lens, np.int64)
        js = np.asarray(js, np.int64)
        starts = np.concatenate([[0], np.cumsum(lens)]).astype(np.int64)
        tok = np.zeros((len(js), len(TAPE_TOKEN_LISTS)), np.uint64)
        ok = np.zeros(tok.shape, bool)
        need: dict[tuple[int, int], list[int]] = {}  # (file, row group) -> the contigs it may hold
        for c, key in enumerate(keys):
            if key is None or lens[c] == 0:
                continue
            for fi, f in enumerate(self._files):
                for j in range(bisect.bisect_left(f.maxs, key), bisect.bisect_right(f.mins, key)):
                    need.setdefault((fi, f.row_groups[j]), []).append(c)

        def use(decoded, locs) -> None:
            levels, rg_keys = decoded
            cs = np.unique(np.concatenate([np.asarray(need[loc], np.int64) for loc in locs]))
            idx = pc.index_in(pa.array([keys[c] for c in cs], pa.binary()), value_set=rg_keys)
            found = idx.is_valid().to_numpy(zero_copy_only=False)
            cs, rows = cs[found], idx.drop_null().to_numpy(zero_copy_only=False).astype(np.int64)
            if not len(cs):
                return
            n = lens[cs]
            # the slots of the found contigs, as positions in js, and the tape row of each
            sid = np.repeat(starts[cs] - np.concatenate([[0], np.cumsum(n)[:-1]]), n) + np.arange(int(n.sum()))
            row, j = np.repeat(rows, n), js[sid]
            for L, (offs, vals, valid) in enumerate(levels):
                o = offs[row]
                inside = (j >= 0) & (j < offs[row + 1] - o)
                if len(vals):
                    pos = np.minimum(o + np.where(inside, j, 0), len(vals) - 1)
                    tok[sid, L] = vals[pos]
                    ok[sid, L] = inside if valid is None else inside & valid[pos]

        with self._lock:
            self._visit(sorted(need), use)
        return tok, ok


class _Lazy:
    """A connection's point lookup over `files` (`make(files)`: `TapeTokenLookup`, `TokenKeyLookup`), built
    (footers read) on first use, so a catalog whose queries never need it never opens the files."""

    def __init__(self, make, files: list[Path]):
        self.make = make
        self.files = files
        self._lookup = None
        self._lock = threading.Lock()

    def get(self):
        with self._lock:
            if self._lookup is None:
                self._lookup = self.make(self.files)
            return self._lookup

    def peek(self):
        return self._lookup


def _offsets(lengths) -> pa.Array:
    """List offsets (int32) from list lengths; refuses more than int32 can index."""
    out = np.concatenate([[0], np.cumsum(np.asarray(lengths, np.int64))])
    if out[-1] > np.iinfo(np.int32).max:
        raise ValueError(f"bb_partner_counts: {out[-1]} values are more than one list can hold")
    return pa.array(out.astype(np.int32), pa.int32())


def _register_tape_udf(con: duckdb.DuckDBPyConnection, tapes: _Lazy) -> None:
    """Register `bb_partner_counts(slots)` over `tapes` (card B-38; card B-31 had `bb_slot_tokens`, per contig).

    Per row: a list of partner slots as structs (c: contig_uid, i: the anchor's slot, j: the partner's slot,
    0-based); returns one struct (t: token, l: loci, p: pairs) per token found on those slots in `tape_tokens`
    (family, NR50 and NR90; null tokens and slots outside the tape give none): `l` the distinct (c, i) whose
    slots hold `t`, `p` the slots that hold it. An empty list when no slot has a token, NULL for a NULL row.
    `partners_of` passes every partner slot of a query in one row, so `TapeTokenLookup.slot_tokens` sees the
    whole query at once and reads each row group once. Counting here (`counts`) replaced an unnest of every
    (slot, token) and two group-bys that took most of the macro's DuckDB time (card B-38 probe of the hub on
    the B-38 test store). Registered as `bb_partner_counts_udf`, shadowing the stored stand-in (`_shadow`,
    card B-34)."""

    out_type = pa.struct([("t", pa.uint64()), ("l", pa.int64()), ("p", pa.int64())])

    def whole(a):  # DuckDB may pass, and pyarrow return, a ChunkedArray
        return a.combine_chunks() if isinstance(a, pa.ChunkedArray) else a

    def counts(flat) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(tokens, loci, pairs) of one row's slots, in no order. The slots are sorted by contig for the lookup;
        the counts are one Arrow hash aggregate (count, count_distinct of the locus) per token: sorting the
        tokens instead took 18 ms of 69 for the hub, and numpy's np.unique 44 ms of 110 (card B-38 probe)."""
        enc = whole(pc.dictionary_encode(whole(pc.struct_field(flat, "c"))))  # the distinct contigs, each slot's
        codes = enc.indices.fill_null(-1).to_numpy(zero_copy_only=False).astype(np.int64)
        js = whole(pc.struct_field(flat, "j")).fill_null(-1).to_numpy(zero_copy_only=False).astype(np.int64)
        iis = whole(pc.struct_field(flat, "i")).fill_null(-1).to_numpy(zero_copy_only=False).astype(np.int64)
        known = np.flatnonzero(codes >= 0)
        if not len(known):
            return np.zeros(0, np.uint64), np.zeros(0, np.int64), np.zeros(0, np.int64)
        order = known[np.argsort(codes[known], kind="stable")]  # the slots grouped by contig
        lens = np.bincount(codes[known], minlength=len(enc.dictionary))
        tok, ok = tapes.get().slot_tokens(enc.dictionary.to_pylist(), lens, js[order])
        t = tok[ok]  # row-major: each slot's tokens in turn
        if not len(t):
            return np.zeros(0, np.uint64), np.zeros(0, np.int64), np.zeros(0, np.int64)
        locus = np.repeat(codes[order] * (int(iis.max()) + 2) + (iis[order] + 1), ok.sum(axis=1))
        g = pa.table({"t": t, "l": locus}).group_by("t", use_threads=False).aggregate(
            [("l", "count"), ("l", "count_distinct")])
        return (g.column("t").to_numpy(), g.column("l_count_distinct").to_numpy().astype(np.int64),
                g.column("l_count").to_numpy().astype(np.int64))

    def bb_partner_counts(wss: pa.Array) -> pa.Array:
        wss = whole(wss)
        valid = wss.is_valid().to_numpy(zero_copy_only=False)
        parts, sizes = [], []
        for r in range(len(wss)):  # one row per partners_of query
            if not valid[r]:
                sizes.append(0)
                continue
            t, loci, pairs = counts(whole(wss[r].values))
            parts.append(pa.StructArray.from_arrays([pa.array(t, pa.uint64()), pa.array(loci, pa.int64()),
                                                     pa.array(pairs, pa.int64())], fields=list(out_type)))
            sizes.append(len(t))
        values = pa.concat_arrays(parts) if parts else pa.array([], out_type)
        return pa.ListArray.from_arrays(_offsets(sizes), values, mask=pa.array(~valid))

    con.create_function("bb_partner_counts_udf", bb_partner_counts, ["STRUCT(c BLOB, i BIGINT, j BIGINT)[]"],
                        "STRUCT(t UBIGINT, l BIGINT, p BIGINT)[]", type="arrow", null_handling="special",
                        side_effects=False)
    _shadow(con, "bb_partner_counts")


def _tape_token_files(con: duckdb.DuckDBPyConnection, root: Path) -> list[Path]:
    """The catalog's `tape_tokens` files of its release and served tape version (`catalog_files`), under `root`."""
    info = dict(con.execute("SELECT key, value FROM catalog_info").fetchall())
    rows = con.execute("SELECT path FROM catalog_files WHERE \"table\" = 'tape_tokens' ORDER BY path").fetchall()
    out = []
    for (p,) in rows:
        parts = classify(p)[1]
        if parts.get("release") == info["release_id"] and parts.get("tape_version") == info.get("tape_version"):
            out.append(Path(root) / p)
    return out


def _token_dict_files(con: duckdb.DuckDBPyConnection, root: Path) -> list[Path]:
    """The catalog's `token_dict` files of its release (`catalog_files`), under `root`."""
    rel = con.execute("SELECT value FROM catalog_info WHERE key = 'release_id'").fetchall()[0][0]
    rows = con.execute("SELECT path FROM catalog_files WHERE \"table\" = 'token_dict' ORDER BY path").fetchall()
    return [Path(root) / p for (p,) in rows if classify(p)[1].get("release") == rel]


class TokenKeyLookup(_Touches):
    """`token_dict` rows by (level, key) from footer statistics, without loading the table (card B-34).

    `token_dict` is sorted by `token_id` = level | source | local, and `local` is the key's rank within
    (level, source), so the rows of each (level, source) run are sorted by key. A row group's `token_id`
    min and max are its first and last rows, which give the runs it starts and ends in (`token_id >> 56`);
    its `key` min and max (statistics the builder writes since card B-34) bound every key in it. For a key
    at a level, a row group is read when its `token_id` range meets the level, its `key` range holds the
    key, and it is either inside one run, or holds part of a run whose other rows do not exclude the key:
    the tail of a run is skipped when the key is not above the run's last row group inside the run, its
    head when the key is not below the run's first such row group. A run that starts and ends inside one
    row group cannot be excluded that way. So a key reads one row group of a run that spans whole row
    groups, and at most two when its run is cut by the row-group boundaries or shares a row group with a
    run of unknown keys (`tests/r2/test_scale_layout.py`). A row group without `key` statistics (a
    `token_dict` written before card B-34) is read whenever it meets the level: right answers, no pruning.

    Counters as `TapeTokenLookup`: `row_groups_read`, `files_hit` per file, `row_groups_total`, and
    `touched_since` (card B-35). At most `cache_row_groups` decoded row groups are held.
    """

    COLUMNS = ("token_id", "level", "key")

    def __init__(self, files, cache_row_groups: int = 16):
        from biobox.l3.tokens import LEVELS  # tokens imports this module

        self.levels = dict(LEVELS)
        self.files = [Path(p) for p in files]
        self.cache_row_groups = cache_row_groups
        self.row_groups_read = 0
        self.files_hit = [0] * len(self.files)
        self._cache: OrderedDict = OrderedDict()
        self._lock = threading.Lock()
        self._touches_init()
        rgs = []  # (token_id min, token_id max, key min, key max, file index, row group)
        self._metadata = []
        for fi, path in enumerate(self.files):
            md = pq.read_metadata(path)
            self._metadata.append(md)
            leaf = {md.schema.column(i).path: i for i in range(md.num_columns)}
            for g in range(md.num_row_groups):
                rg = md.row_group(g)
                if rg.num_rows == 0:
                    continue
                t, k = rg.column(leaf["token_id"]).statistics, rg.column(leaf["key"]).statistics
                if t is None or not t.has_min_max:
                    raise ValueError(f"{path}: row group {g} has no min/max statistics on token_id (LY-03)")
                # a token_dict written before card B-34 has no key statistics: its row groups are read
                # whenever they meet the level, as the stored SQL macro reads them
                k = (k.min, k.max) if k is not None and k.has_min_max else (None, None)
                rgs.append((t.min, t.max, *k, fi, g))
        rgs.sort(key=lambda r: r[0])
        for prev, nxt in zip(rgs, rgs[1:]):
            if nxt[0] <= prev[1]:
                raise ValueError(f"token_dict row groups overlap in token_id ({prev[0]}..{prev[1]} and "
                                 f"{nxt[0]}..{nxt[1]}): the files are not one sorted table (LY-02)")
        self.row_groups_total = len(rgs)
        self._tmin = [r[0] for r in rgs]
        self._tmax = [r[1] for r in rgs]
        self._kmin = [r[2] for r in rgs]
        self._kmax = [r[3] for r in rgs]
        self._where = [(r[4], r[5]) for r in rgs]
        self._first = [r[0] >> 56 for r in rgs]  # the run (level << 5 | source) of the row group's first row
        self._last = [r[1] >> 56 for r in rgs]
        self._first_pure: dict[int, int] = {}  # run -> its first and last row group that holds that run only
        self._last_pure: dict[int, int] = {}
        for i, (a, b) in enumerate(zip(self._first, self._last)):
            if a == b:
                self._first_pure.setdefault(a, i)
                self._last_pure[a] = i

    def candidates(self, level: str, key: str) -> list[int]:
        """The row groups (indices in token_id order) that may hold `key` at `level`; see the class docstring."""
        lv = self.levels.get(level)
        if lv is None:
            return []
        lo, hi = lv << 5, (lv << 5) | 31  # the runs of the level
        first = bisect.bisect_left(self._last, lo)  # row groups ending before the level do not meet it
        out = []
        for i in range(first, len(self._first)):
            a, b = self._first[i], self._last[i]
            if a > hi:
                break
            if self._kmin[i] is not None and not self._kmin[i] <= key <= self._kmax[i]:
                continue
            if a == b or self._mixed_may_hold(i, key, max(a, lo), min(b, hi)):
                out.append(i)
        return out

    def _mixed_may_hold(self, i: int, key: str, lo: int, hi: int) -> bool:
        """Whether row group `i`, which holds rows of several runs, may hold `key` in one of the runs `lo`..`hi`."""
        for r in range(lo, hi + 1):
            if r == self._first[i] and self._last_pure.get(r) == i - 1 and self._kmax[i - 1] is not None:
                if key > self._kmax[i - 1]:  # the run's tail: keys above its row groups before this one
                    return True
            elif r == self._last[i] and self._first_pure.get(r) == i + 1 and self._kmin[i + 1] is not None:
                if key < self._kmin[i + 1]:  # the run's head: keys below its row groups after this one
                    return True
            else:
                return True  # a run (or a run number without rows) whose keys here nothing bounds
        return False

    def _decoded(self, i: int) -> pa.Table:
        fi, g = self._where[i]
        hit = self._cache.pop(i, None)
        if hit is not None:
            self._cache[i] = hit
            self._touch(i, fi, False)
            return hit
        t = pq.ParquetFile(self.files[fi], metadata=self._metadata[fi]).read_row_group(g, columns=list(self.COLUMNS))
        self.row_groups_read += 1
        self.files_hit[fi] += 1
        self._touch(i, fi, True)
        if self.cache_row_groups > 0:
            self._cache[i] = t
            if len(self._cache) > self.cache_row_groups:
                self._cache.popitem(last=False)
        return t

    def tokens(self, pairs) -> dict[tuple[str, str], list[int]]:
        """The token ids of each distinct (level, key) of `pairs` (an empty list for an unknown one); each row
        group is read once for all the keys that need it."""
        need: dict[int, dict[str, list[str]]] = {}  # row group -> level -> keys
        out: dict[tuple[str, str], list[int]] = {}
        for level, key in dict.fromkeys(pairs):
            out[(level, key)] = []
            for i in self.candidates(level, key):
                need.setdefault(i, {}).setdefault(level, []).append(key)
        with self._lock:
            for i in sorted(need):
                t = self._decoded(i)
                for level, keys in need[i].items():
                    m = pc.and_(pc.equal(t.column("level"), level), pc.is_in(t.column("key"), value_set=pa.array(keys)))
                    hit = t.filter(m)
                    for k, tid in zip(hit.column("key").to_pylist(), hit.column("token_id").to_pylist()):
                        out[(level, k)].append(tid)
        return out


def _register_token_udf(con: duckdb.DuckDBPyConnection, tokens: _Lazy) -> None:
    """`token_of(level, key)` as a point lookup over the release's `token_dict` footers (card B-34).

    The stored `token_of` is a scalar subquery on `token_dict`, which is sorted by `token_id`: by key it can
    prune only with a literal and only on the `key` statistics, and for keys from a joined table it reads every
    row group. `open_catalog` shadows it, on its connection only, with a temporary macro over
    `bb_token_of_udf`, which resolves each distinct key of a vector with `TokenKeyLookup`. The macro calls the
    function directly: wrapped in a scalar subquery, 100 joined keys took 0.059 s in 8 calls of 20 to 30 rows,
    against 0.016 s in 2 calls of 100 (measured, wilkcomp02, 2026-11-10, 1,000x synthetic store), and for a
    literal or a parameter DuckDB calls it as often either way (`token_rows`, `q1_count`). An unknown
    level or key, or a NULL, gives NULL; a key with tokens under two sources of one level raises (there is no
    one token to give)."""

    def bb_token_of(levels: pa.Array, keys: pa.Array) -> pa.Array:
        rows = list(zip(levels.to_pylist(), keys.to_pylist()))
        got = tokens.get().tokens([r for r in rows if r[0] is not None and r[1] is not None])
        out = []
        for level, key in rows:
            tids = got.get((level, key), [])
            if len(tids) > 1:
                raise ValueError(f"token_of: key {key!r} names {len(tids)} tokens at level {level} "
                                 f"(one per source): {sorted(tids)}")
            out.append(tids[0] if tids else None)
        return pa.array(out, pa.uint64())

    con.create_function("bb_token_of_udf", bb_token_of, ["VARCHAR", "VARCHAR"], "UBIGINT",
                        type="arrow", null_handling="special", side_effects=False)
    con.execute("CREATE OR REPLACE TEMP MACRO token_of(_level, _key) AS bb_token_of_udf(_level, _key)")


SEQ_TABLES = ("seq_locator", "contig_index", "contig_occurrence")  # what SeqStore reads


def catalog_release_files(con: duckdb.DuckDBPyConnection, root: Path) -> ReleaseFiles:
    """The release's sequence-store files as the catalog lists them (`catalog_files`, from its manifest),
    under `root`: the sequence functions read only these (card B-28)."""
    info = dict(con.execute("SELECT key, value FROM catalog_info").fetchall())
    rows = con.execute('SELECT path, "table" FROM catalog_files WHERE "table" IN (SELECT unnest(?::VARCHAR[])) '
                       "ORDER BY path", [list(SEQ_TABLES)]).fetchall()
    files = [{"path": p, "table": t, "partitions": classify(p)[1]} for p, t in rows]
    if "stages" in info:
        stages = json.loads(info["stages"])
    else:  # a catalog from before card B-28: the stages its files name
        stages = sorted({f["partitions"]["stage"] for f in files if "stage" in f["partitions"]})
    rel = info["release_id"]
    return ReleaseFiles(root, rel, stages, {"release_id": rel, "stages": stages, "files": files, "artifacts": []})


def seq_reader(con: duckdb.DuckDBPyConnection):
    """The SeqStore behind a connection's sequence functions, for its LY-06 counters."""
    return _READERS.get(id(con))


def open_catalog(path: Path, local_dir: Path | None = None, threads: int | None = None,
                 memory_limit: str | None = None, temp_directory: Path | None = None,
                 root: Path | None = None, cache_blocks: int = 64, runs=(),
                 pack: Path | None = None) -> duckdb.DuckDBPyConnection:
    """Open read-only, optionally from a copy on local disk (never open a DuckDB file over NFS).

    The catalog file is attached read-only, under its file name without `.duckdb` (the name a plain
    `duckdb.connect` gives it), to a DuckDB instance of this connection's own, and made the default
    database (card B-34; module docstring). `threads`, `memory_limit` and `temp_directory` configure that
    instance. DuckDB spills into `temp_directory`, by default a directory of this process under the system
    temporary directory, never a path under the data root or beside the catalog file (card W-18, W-03's rule:
    a reader may not be able to write the root, and one who can must not write into it).

    Registers the Python functions of `PY_FUNCTIONS` over the store under `root` (default: the root the
    catalog was built from, in `catalog_info`), with an LRU of `cache_blocks` decoded blocks per shard, and
    shadows the stored stand-ins and the lookup macros with them. Sets the join filter to the row groups of
    the release's largest table (`set_join_key_filter`, card B-27).

    `runs` (search run directories) and `pack` (the release's search pack directory), both outside the data
    root (card W-03), are given to `attach_search`.
    """
    path = Path(path)
    if local_dir is not None:
        local = Path(local_dir) / path.name
        local.parent.mkdir(parents=True, exist_ok=True)
        if not local.exists() or local.stat().st_size != path.stat().st_size or local.stat().st_mtime < path.stat().st_mtime:
            tmp = local.with_name(f".{local.name}.tmp-{os.getpid()}")
            shutil.copyfile(path, tmp)
            os.replace(tmp, local)
        path = local
    spill = temp_directory if temp_directory is not None else Path(tempfile.gettempdir()) / f"biobox-duckdb-{os.getpid()}"
    config = {"temp_directory": str(spill)}
    if threads is not None:
        config["threads"] = threads
    if memory_limit is not None:
        config["memory_limit"] = memory_limit
    con = duckdb.connect(":memory:", config=config)  # a new instance: its functions are this connection's
    name = '"' + path.name.removesuffix(".duckdb").replace('"', '""') + '"'
    try:
        con.execute(f"ATTACH {_sql_str(str(path))} AS {name} (READ_ONLY)")
        con.execute(f"USE {name}")
    except BaseException:
        con.close()
        raise
    if root is None:
        root = Path(con.execute("SELECT value FROM catalog_info WHERE key = 'root'").fetchall()[0][0])
    # keys from a joined table prune to the row groups that hold them (card B-27; set_join_key_filter)
    (largest,), = con.execute('SELECT coalesce(max(n), 0)::BIGINT FROM (SELECT sum(row_groups) AS n FROM catalog_files '
                              'GROUP BY "table")').fetchall()
    set_join_key_filter(con, largest)
    from biobox.l0seq.store import SeqStore

    files = catalog_release_files(con, Path(root))
    store = SeqStore(root, files, cache_blocks) if files.table("seq_locator") else None
    _register_seq_udf(con, store)
    _READERS[id(con)] = store
    if store is not None:
        _register_lookup_udf(con, store)
    # the tape tokens and token keys the views read: under the build root, as every view's files (`root`
    # moves the shards only)
    built = Path(con.execute("SELECT value FROM catalog_info WHERE key = 'root'").fetchall()[0][0])
    (n_threads,), = con.execute("SELECT current_setting('threads')::BIGINT").fetchall()
    tapes = _Lazy(lambda files: TapeTokenLookup(files, threads=n_threads), _tape_token_files(con, built))
    _register_tape_udf(con, tapes)
    # a catalog built before card B-38 stores a partners_of that calls bb_slot_tokens per contig
    con.execute(partners_of_sql("CREATE OR REPLACE TEMP MACRO"))
    keys = _Lazy(TokenKeyLookup, _token_dict_files(con, built))
    _register_token_udf(con, keys)
    _LOOKUPS[id(con)] = [("bb_partner_counts", tapes), ("bb_token_of", keys)]
    attach_search(con, runs, pack)
    return con


SEARCH_RUN_FILES = {"search_run": "run.parquet", "search_hit": "hits.parquet", "search_profile": "profiles.parquet",
                    "neighborhoods": "neighborhoods.parquet"}


def _run_label(d: Path) -> str:
    """A run directory's `run` in the search views: the id of its name `run=<run_id>`, else (a directory
    named freely with `bio-box search run --out`) the run id its record holds."""
    if d.name.startswith("run="):
        return d.name.removeprefix("run=")
    return pq.read_table(d / "run.parquet", columns=["run_id"]).column("run_id")[0].as_py().hex()


def attach_search(con: duckdb.DuckDBPyConnection, runs=(), pack: Path | None = None) -> dict[str, int]:
    """Temporary views over closed search runs and the release's pack regions, shadowing the catalog's
    empty views (data_layer_and_tests.md §2B.3). Runs and packs live outside the data root (card W-03) and
    close after the catalog is built, so they are given when the catalog is opened: `runs` is a list of run
    directories (each holds its `run.parquet`; FileNotFoundError otherwise), `pack` the directory of the
    release's search pack (ValueError for a pack of another release). Each run file view carries its run's
    `run` (`_run_label`)."""
    dirs = sorted({Path(d) for d in runs})
    for d in dirs:
        if not (d / "run.parquet").is_file():
            raise FileNotFoundError(f"{d}: not a closed search run (no run.parquet)")
    labels = {d: _run_label(d) for d in dirs}
    counts = {}
    for name, fname in SEARCH_RUN_FILES.items():
        spec = TABLES[name]
        found = [(d / fname, labels[d]) for d in dirs if (d / fname).is_file()]
        counts[name] = len(found)
        if found:
            cols = [f'"{f.name}"' for f in spec.schema] + ['"run"']
            parts = " UNION ALL BY NAME ".join(
                f"SELECT *, {_sql_str(label)} AS \"run\" FROM read_parquet({_sql_str(str(p))}, hive_partitioning = false)"
                for p, label in found)
            con.execute(f"CREATE OR REPLACE TEMP VIEW {name} AS SELECT {', '.join(cols)} FROM ({parts})")
        else:  # same columns with no runs, so queries over `run` bind either way
            cols = [f'CAST(NULL AS {duck_type(f.type)}) AS "{f.name}"' for f in spec.schema]
            con.execute(f"CREATE OR REPLACE TEMP VIEW {name} AS SELECT {', '.join(cols)}, "
                        "CAST(NULL AS VARCHAR) AS \"run\" WHERE false")
    if pack is not None:
        pack = Path(pack)
        man = json.loads((pack / "pack_manifest.json").read_text())
        rel = con.execute("SELECT value FROM catalog_info WHERE key = 'release_id'").fetchall()
        if rel and man["release_id"] != rel[0][0]:
            raise ValueError(f"{pack}: a search pack of release {man['release_id']}, not of the catalog's {rel[0][0]}")
        metas = [pack / m["path"] for m in sorted(man["meta"], key=lambda m: m["path"])]  # the pack's own list
        counts["search_region"] = len(metas)
        if metas:
            files = "[" + ", ".join(_sql_str(str(x)) for x in metas) + "]"
            con.execute(f"CREATE OR REPLACE TEMP VIEW search_region AS SELECT * FROM read_parquet({files}, "
                        "hive_partitioning = false)")
    return counts


# --- LY-06: storage units read per query -------------------------------------


def profile_query(con: duckdb.DuckDBPyConnection, sql: str, params: list | None = None) -> dict:
    """Run `sql` with profiling and return what it read, from DuckDB's own profiler.

    `scans` has one entry per Parquet scan: files read, row groups scanned and the total the
    scan could have read. `bytes_read` is DuckDB's I/O counter for the query. `bgzf` is what
    the sequence functions read from the shards (blocks decompressed, compressed bytes read,
    bytes decoded), and `index` the row groups their point lookups read from `seq_locator` and
    `contig_index` (`KeyLookup.row_groups_read`; a cache hit counts nothing); both are zero when
    the query reads no sequence. The `tape_tokens` row groups that `bb_partner_counts` reads by point lookup
    (`partners_of`, card B-31) are one more entry of `scans`, with `"lookup": "bb_partner_counts"` and
    `rows_scanned` None, and so are the `token_dict` row groups of `token_of` (`"lookup": "bb_token_of"`,
    card B-34). For these two, `row_groups_scanned` is the distinct row groups the query touched, whether
    the lookup decoded them from the file (`row_groups_read`, decodes) or its row-group cache, which lives
    as long as the connection, served them (`row_groups_cached`), so a query reports the same rows read
    cold or after another query on the connection (card B-35). The lookups `default_digest`, `digest_len`
    and `contig_len_of` of `open_catalog` count in `index` (decodes only). Counts, not time, decide
    layouts (LY-06).
    """
    store = seq_reader(con)
    before = store.counters() if store is not None else None

    def index_counts() -> dict[str, int]:
        if store is None:
            return {"locator_row_groups_read": 0, "index_row_groups_read": 0}
        return {"locator_row_groups_read": store.locator.row_groups_read,
                "index_row_groups_read": store.index.row_groups_read}

    index_before = index_counts()
    lookups = _LOOKUPS.get(id(con), [])

    def lookup_mark(lazy) -> tuple[int, int]:  # (stamp, decodes) before the query; a lookup built by it starts at 0
        t = lazy.peek()
        return (0, 0) if t is None else (t.stamp, t.row_groups_read)

    lookups_before = [lookup_mark(lazy) for _, lazy in lookups]
    with tempfile.TemporaryDirectory() as d:
        out = Path(d) / "profile.json"
        con.execute("PRAGMA enable_profiling = 'json'")
        con.execute(f"PRAGMA profiling_output = {_sql_str(str(out))}")  # keeps the profile off stdout
        try:
            result = con.execute(sql, params or []).fetchall()
        finally:
            con.execute("PRAGMA disable_profiling")
        # DuckDB 2.0 writes no profile for some plans (count(*) filtered only on hive keys);
        # callers that count units must use a query that scans (e.g. sum of a column)
        info = json.loads(out.read_text()) if out.exists() else {}
    scans = []

    def walk(node):
        if node.get("type") == "TABLE_SCAN" and "row_groups_scanned" in node:
            extra = node.get("extra_info", {})
            names = extra.get("Filename(s)", "")
            scans.append({
                "files": [n for n in names.split("\n") if n] if isinstance(names, str) else list(names),
                "files_read": int(extra.get("Total Files Read", 0) or 0),
                "row_groups_scanned": node["row_groups_scanned"],
                "row_groups_total": node.get("total_row_groups_to_scan"),
                "rows_scanned": node.get("rows_scanned"),
            })
        for c in node.get("children", []):
            walk(c)

    for op in info.get("operator", []) if isinstance(info.get("operator"), list) else [info]:
        walk(op)
    for c in info.get("children", []):
        walk(c)
    # the point lookups of bb_partner_counts (partners_of, cards B-31, B-38) and of token_of (card B-34) read Parquet
    # outside DuckDB's scans: each is listed as one more scan of its table (`tape_tokens`, `token_dict`), the
    # distinct row groups the query touched (decoded or served by the cache, card B-35) of the files' total
    for (fn, lazy), (mark, n_before) in zip(lookups, lookups_before):
        t = lazy.peek()
        if t is None:
            continue
        touched, decoded = t.touched_since(mark)
        if touched:
            scans.append({
                "files": [str(p) for p in t.files],
                "files_read": len(set(touched.values())),
                "row_groups_scanned": len(touched),
                "row_groups_total": t.row_groups_total,
                "rows_scanned": None,
                "lookup": fn,
                "row_groups_read": t.row_groups_read - n_before,
                "row_groups_cached": len(touched.keys() - decoded),
            })
    bgzf =({k: v - before[k] for k, v in store.counters().items()} if store is not None
            else {"blocks_read": 0, "compressed_bytes_read": 0, "uncompressed_bytes_decoded": 0})
    index = {k: v - index_before[k] for k, v in index_counts().items()}
    return {"result": result, "scans": scans, "bytes_read": info.get("io", {}).get("total_bytes_read"), "bgzf": bgzf,
            "index": index}
