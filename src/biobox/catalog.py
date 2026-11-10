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
table (AGENTS §4.4, card B-14). A catalog opened with plain `duckdb` has every view and
macro, but the sequence macros fail with "Catalog Error: Scalar Function with name bb_seq
does not exist!" (`bb_seq_circular`, `bb_seq_raw` for the other two), and so does `partners_of`
(`bb_slot_tokens`, its point lookup of `tape_tokens`, card B-31). The lookup macros `default_digest`,
`digest_len` and `contig_len_of` work there as SQL that scans the index views; `open_catalog` shadows
them on its connection with point lookups over the same `SeqStore` (card B-31), so a joined key reads
one index row group, not all.

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


# Trees under the data root that carry their own manifests and are not release files: the
# search pack (search_pack/release=<R>/pack_manifest.json, verified by C-12) and search runs
# (search/run_index.json), which are written after the release is signed; and run-step records
# (runs/stage=<S>/<step>.json, card C-02), which hold timestamps and host names.
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


COMPLEMENT = ("ACGTRYKMBVDHNacgtrykmbvdhn", "TGCAYRMKVBHDNtgcayrmkvbhdn")


def pair_macro_sql() -> list[str]:
    """Neighbour macros over `token_postings` (card B-11; planning/tasks/TASKS.md §2; AGENTS §3.5).

    Pairs are never stored: Q1 and Q4 self-join the postings of this catalog's release and tape
    version (`catalog_tape_version()`) on `contig_uid`. Each side of a lookup by token reads one
    `token_bucket` and the row groups whose `token_id` range holds the token, also when the token comes
    from `token_of` or a joined table (card B-31: `token_rows` joins the token in). Q3 (`partners_of`)
    takes its anchors from the postings and the tokens near them from `tape_tokens` of the anchors'
    contigs (`bb_slot_tokens`), so its partner side reads row groups in proportion to the anchors'
    contigs, not to the release. The geometry helpers
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
        # Q3 (counts; z_score ranking arrives with pair_stats): every token at every level within k slots
        # of the loci of A. Anchors come from the postings of A (by token_id); the partner slots are the
        # slots within k of each anchor on its contig (slot_dist 1..k, each slot once per anchor), and their
        # tokens come from `tape_tokens` of the anchors' contigs through `bb_slot_tokens` (card B-31), a
        # point lookup by contig_uid that reads only the row groups holding those contigs. Postings are
        # sorted by token_id, so the earlier near side (postings filtered by contig_uid) read the whole
        # token level of the release. A hub is still refused here (its anchors are many)
        f"""CREATE MACRO partners_of(_a, _k) AS TABLE
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
                   _per AS (
                     SELECT contig_uid, list(_i ORDER BY _i, _j) AS _is, list(_j ORDER BY _i, _j) AS _js
                     FROM _win WHERE _j >= 0 AND _j < _n AND _j <> _i GROUP BY contig_uid),
                   _slots AS (
                     SELECT contig_uid, unnest(_is) AS _i, unnest(bb_slot_tokens(contig_uid, _js)) AS _toks FROM _per),
                   _near AS (SELECT contig_uid, _i, unnest(_toks) AS _tok FROM _slots)
              SELECT _tok AS partner_token, count(DISTINCT (contig_uid, _i)) AS n_loci, count(*) AS n_pairs
              FROM _near GROUP BY _tok ORDER BY n_loci DESC, partner_token""",
    ]


def macro_sql() -> list[str]:
    """Sequence macros (data_layer_and_tests.md §2A.3). 0-based half-open; strand +1/-1.

    Parameter names start with an underscore so they never collide with column names. The
    sequence macros only call `bb_seq`, `bb_seq_circular` and `bb_seq_raw` (registered by
    `open_catalog`), which check bounds and strand, resolve keys and read the shards; the lookup
    macros `default_digest`, `digest_len` and `contig_len_of` are for users. As stored here they scan
    the index views (for a plain `duckdb` connection); `open_catalog` replaces them on its connection
    with point lookups (`_register_lookup_udf`, card B-31). The sequence macros do not call them.
    """
    return [
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
            con.execute(_view_sql(TABLES[name], sorted(by_table[name]), widest.get(name)))
        # derived view, not a registry table: one row per contig_uid (contig_occurrence has one per source record)
        con.execute(
            "CREATE VIEW contig_label AS SELECT contig_uid, "
            "arg_min(accession, (l0_object_sha256, record_ordinal)) AS accession, "
            "arg_min(original_header, (l0_object_sha256, record_ordinal)) AS original_header, "
            "count(*) AS n_occurrences FROM contig_occurrence GROUP BY contig_uid")
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
        # DuckDB binds a macro's functions when it is created: stand-ins for bb_seq,
        # bb_seq_circular and bb_seq_raw let the macros be stored without a store; open_catalog
        # registers the real readers; macros come last because they read catalog_info
        _register_seq_udf(con, None)
        _register_tape_udf(con, None)  # partners_of calls bb_slot_tokens (card B-31)
        for stmt in macro_sql():
            con.execute(stmt)
        con.execute("CHECKPOINT")
    finally:
        con.close()
    os.replace(tmp, out_path)
    return out_path


_READERS: dict[int, object] = {}
_TAPES: dict[int, object] = {}  # id(connection) -> its _LazyTapes (card B-31)


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
    index views. A NULL in any argument gives NULL. With `store` None (the stand-ins
    `generate_catalog` stores the macros against, or a root without sequence) they raise.
    """

    def need() -> None:
        if store is None:
            raise RuntimeError("bb_seq: open the catalog with biobox.catalog.open_catalog to read sequence")

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
    con.create_function("bb_seq", bb_seq, ["BLOB", "BIGINT", "BIGINT", "BIGINT"], "VARCHAR", **opts)
    con.create_function("bb_seq_circular", bb_seq_circular, ["BLOB", "BIGINT", "BIGINT", "BIGINT"], "VARCHAR", **opts)
    con.create_function("bb_seq_raw", bb_seq_raw, ["BLOB", "BIGINT", "BIGINT"], "VARCHAR", **opts)


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


TAPE_TOKEN_LISTS = ("family_tokens", "nr50_tokens", "nr90_tokens")  # bb_slot_tokens returns a slot's tokens in this order


class TapeTokenLookup(KeyLookup):
    """`tape_tokens` rows by `contig_uid` from the release's files of one tape version (card B-31).

    The footer statistics and binary search of `KeyLookup` (one row group read per contig, at most
    `cache_row_groups` decoded row groups held), with its own decoding: a row group is held as numpy arrays
    per token list (offsets, values with nulls as 0, validity), so finding a slot's tokens is array
    indexing, never a pyarrow call or a Python object per slot. `row(key)` gives (levels, row index);
    `get` of the base class is not available. `files_hit` counts row groups read per file (a cache hit
    counts nothing).
    """

    def __init__(self, files, cache_row_groups: int = 32):
        super().__init__(files, "contig_uid", TAPE_TOKEN_LISTS, cache_row_groups)
        self.files_hit = [0] * len(self._files)
        self.row_groups_total = sum(len(f.row_groups) for f in self._files)

    def _decoded(self, fi: int, g: int):
        hit = self._cache.pop((fi, g), None)
        if hit is not None:
            self._cache[(fi, g)] = hit
            return hit
        f = self._files[fi]
        t = pq.ParquetFile(f.path, metadata=f.metadata).read_row_group(g, columns=self.columns)
        self.row_groups_read += 1
        self.files_hit[fi] += 1
        levels = []
        for c in TAPE_TOKEN_LISTS:
            a = t.column(c).combine_chunks()
            v = a.values  # a row group read whole: offsets start at 0 and index `values`
            levels.append((a.offsets.to_numpy(zero_copy_only=False).astype(np.int64),
                           v.fill_null(0).to_numpy(zero_copy_only=False).astype(np.uint64),
                           v.is_valid().to_numpy(zero_copy_only=False)))
        out = (levels, t.column(self.key_column).to_pylist())
        if self.cache_row_groups > 0:
            self._cache[(fi, g)] = out
            if len(self._cache) > self.cache_row_groups:
                self._cache.popitem(last=False)
        return out

    def get(self, key):  # noqa: D102
        raise NotImplementedError("TapeTokenLookup: use row()")

    def row(self, key: bytes) -> tuple[list, int] | None:
        """(the row group's token lists as (offsets, values, valid) per level, the contig's row in it), or None
        when no file holds the contig."""
        with self._lock:
            for fi, f in enumerate(self._files):
                first, last = bisect.bisect_left(f.maxs, key), bisect.bisect_right(f.mins, key)
                for j in range(first, last):
                    levels, keys = self._decoded(fi, f.row_groups[j])
                    a = bisect.bisect_left(keys, key)
                    if a < len(keys) and keys[a] == key:
                        return levels, a
        return None


class _LazyTapes:
    """The connection's `TapeTokenLookup`, built (footers read) on first use, so a catalog whose queries never
    ask for tape tokens never opens them."""

    def __init__(self, files: list[Path]):
        self.files = files
        self._lookup: TapeTokenLookup | None = None
        self._lock = threading.Lock()

    def get(self) -> TapeTokenLookup:
        with self._lock:
            if self._lookup is None:
                self._lookup = TapeTokenLookup(self.files)
            return self._lookup

    def peek(self) -> TapeTokenLookup | None:
        return self._lookup


def _register_tape_udf(con: duckdb.DuckDBPyConnection, tapes: _LazyTapes | None) -> None:
    """Register `bb_slot_tokens(contig_uid, js)` over `tapes` (card B-31): for each 0-based slot of the list
    `js`, the non-null tokens of that slot of the contig's tape tokens, family then NR50 then NR90; an empty
    list for a slot outside the tape or a contig without tape tokens, NULL for a NULL argument. One call per
    contig (the macro groups its slots), so a contig is looked up once and its slots are read with numpy.
    With `tapes` None (the stand-in `generate_catalog` stores the macros against) it raises."""

    def bb_slot_tokens(cuids: pa.Array, jss: pa.Array) -> pa.Array:
        if tapes is None:
            raise RuntimeError("bb_slot_tokens: open the catalog with biobox.catalog.open_catalog to read tape tokens")
        lookup = tapes.get()
        cus = cuids.to_pylist()
        lens = pc.list_value_length(jss).fill_null(0).to_numpy(zero_copy_only=False).astype(np.int64)
        flat = pc.list_flatten(jss).fill_null(-1).to_numpy(zero_copy_only=False).astype(np.int64)
        starts = np.concatenate([[0], np.cumsum(lens)])
        null = np.array([cu is None or not jss[r].is_valid for r, cu in enumerate(cus)], dtype=bool)
        toks, masks = [], []
        for r, cu in enumerate(cus):
            if null[r]:
                continue
            js = flat[starts[r]:starts[r + 1]]
            hit = lookup.row(cu)
            tok = np.zeros((len(js), len(TAPE_TOKEN_LISTS)), np.uint64)
            ok = np.zeros(tok.shape, bool)
            if hit is not None:
                levels, i = hit
                for L, (offs, vals, valid) in enumerate(levels):
                    o, n = offs[i], offs[i + 1] - offs[i]
                    inside = (js >= 0) & (js < n)
                    pos = o + np.where(inside, js, 0)
                    if len(vals):
                        tok[:, L] = vals[np.minimum(pos, len(vals) - 1)]
                        ok[:, L] = inside & valid[np.minimum(pos, len(vals) - 1)]
            toks.append(tok)
            masks.append(ok)
        tok = np.concatenate(toks) if toks else np.zeros((0, len(TAPE_TOKEN_LISTS)), np.uint64)
        ok = np.concatenate(masks) if masks else np.zeros(tok.shape, bool)
        inner_offsets = np.concatenate([[0], np.cumsum(ok.sum(axis=1))]).astype(np.int32)
        inner = pa.ListArray.from_arrays(pa.array(inner_offsets, pa.int32()), pa.array(tok[ok], pa.uint64()))
        outer_lens = np.where(null, 0, lens)
        outer_offsets = np.concatenate([[0], np.cumsum(outer_lens)]).astype(np.int32)
        return pa.ListArray.from_arrays(pa.array(outer_offsets, pa.int32()), inner, mask=pa.array(null))

    con.create_function("bb_slot_tokens", bb_slot_tokens, ["BLOB", "BIGINT[]"], "UBIGINT[][]",
                        type="arrow", null_handling="special", side_effects=False)


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
                 root: Path | None = None, cache_blocks: int = 64) -> duckdb.DuckDBPyConnection:
    """Open read-only, optionally from a copy on local disk (never open a DuckDB file over NFS).

    Registers `bb_seq`, `bb_seq_circular` and `bb_seq_raw` over the store under `root`
    (default: the root the catalog was built from, in `catalog_info`), with an LRU of
    `cache_blocks` decoded blocks per shard.
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
    config = {"access_mode": "READ_ONLY"}
    if threads is not None:
        config["threads"] = threads
    if memory_limit is not None:
        config["memory_limit"] = memory_limit
    if temp_directory is not None:
        config["temp_directory"] = str(temp_directory)
    con = duckdb.connect(str(path), read_only=True, config=config)
    if root is None:
        root = Path(con.execute("SELECT value FROM catalog_info WHERE key = 'root'").fetchall()[0][0])
    from biobox.l0seq.store import SeqStore

    files = catalog_release_files(con, Path(root))
    store = SeqStore(root, files, cache_blocks) if files.table("seq_locator") else None
    _register_seq_udf(con, store)
    _READERS[id(con)] = store
    if store is not None:
        _register_lookup_udf(con, store)
    # the tape tokens the views read: under the build root, as every view's files (`root` moves the shards only)
    built = Path(con.execute("SELECT value FROM catalog_info WHERE key = 'root'").fetchall()[0][0])
    tapes = _LazyTapes(_tape_token_files(con, built))
    _register_tape_udf(con, tapes)
    _TAPES[id(con)] = tapes
    attach_search(con, Path(root))
    return con


SEARCH_RUN_FILES = {"search_run": "run.parquet", "search_hit": "hits.parquet", "search_profile": "profiles.parquet",
                    "neighborhoods": "neighborhoods.parquet"}


def attach_search(con: duckdb.DuckDBPyConnection, root: Path) -> dict[str, int]:
    """Temporary views over closed search runs (`search/run_index.json`) and the release's pack
    regions, shadowing the catalog's empty views (data_layer_and_tests.md §2B.3). Runs close
    after the catalog is built, so the views are made when the catalog is opened; each run
    file view carries the run id from its path as `run`."""
    idx_path = root / "search/run_index.json"
    runs = json.loads(idx_path.read_text())["runs"] if idx_path.exists() else []
    counts = {}
    for name, fname in SEARCH_RUN_FILES.items():
        spec = TABLES[name]
        paths = [str(root / f"search/run={r['run_id']}/{fname}") for r in runs if fname in r["files"]]
        counts[name] = len(paths)
        if paths:
            cols = [f'"{f.name}"' for f in spec.schema] + ['"run"']
            files = "[" + ", ".join(_sql_str(x) for x in sorted(paths)) + "]"
            con.execute(f"CREATE OR REPLACE TEMP VIEW {name} AS SELECT {', '.join(cols)} FROM read_parquet("
                        f"{files}, hive_partitioning = true, hive_types = {{'run': 'VARCHAR'}}, union_by_name = true)")
        else:  # same columns with no runs, so queries over `run` bind either way
            cols = [f'CAST(NULL AS {duck_type(f.type)}) AS "{f.name}"' for f in spec.schema]
            con.execute(f"CREATE OR REPLACE TEMP VIEW {name} AS SELECT {', '.join(cols)}, "
                        "CAST(NULL AS VARCHAR) AS \"run\" WHERE false")
    rel = con.execute("SELECT value FROM catalog_info WHERE key = 'release_id'").fetchall()
    if rel:
        metas = sorted((root / f"search_pack/release={rel[0][0]}").glob("*/*_meta.parquet"))
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
    the query reads no sequence. The `tape_tokens` row groups that `bb_slot_tokens` reads by point lookup
    (`partners_of`, card B-31) are one more entry of `scans`, with `"lookup": "bb_slot_tokens"` and
    `rows_scanned` None. The lookups `default_digest`, `digest_len` and `contig_len_of` of `open_catalog`
    count in `index`. Counts, not time, decide layouts (LY-06).
    """
    store = seq_reader(con)
    before = store.counters() if store is not None else None

    def index_counts() -> dict[str, int]:
        if store is None:
            return {"locator_row_groups_read": 0, "index_row_groups_read": 0}
        return {"locator_row_groups_read": store.locator.row_groups_read,
                "index_row_groups_read": store.index.row_groups_read}

    index_before = index_counts()
    lazy = _TAPES.get(id(con))

    def tape_counts() -> tuple[int, list[int]]:
        t = lazy.peek() if lazy is not None else None
        return (0, []) if t is None else (t.row_groups_read, list(t.files_hit))

    tapes_before = tape_counts()
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
    # the tape-token point lookups of bb_slot_tokens (partners_of, card B-31) read Parquet outside DuckDB's
    # scans: they are listed as one more scan of `tape_tokens`, row groups read of the files' total
    tapes_after = tape_counts()
    if tapes_after[0] > tapes_before[0]:
        t = lazy.peek()
        hit_before = tapes_before[1] or [0] * len(tapes_after[1])
        scans.append({
            "files": [str(p) for p in t.files],
            "files_read": sum(a > b for a, b in zip(tapes_after[1], hit_before)),
            "row_groups_scanned": tapes_after[0] - tapes_before[0],
            "row_groups_total": t.row_groups_total,
            "rows_scanned": None,
            "lookup": "bb_slot_tokens",
        })
    bgzf =({k: v - before[k] for k, v in store.counters().items()} if store is not None
            else {"blocks_read": 0, "compressed_bytes_read": 0, "uncompressed_bytes_decoded": 0})
    index = {k: v - index_before[k] for k, v in index_counts().items()}
    return {"result": result, "scans": scans, "bytes_read": info.get("io", {}).get("total_bytes_read"), "bgzf": bgzf,
            "index": index}
