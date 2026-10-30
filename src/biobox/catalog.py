"""Release manifest and generated read-only DuckDB catalog (operation O-CATALOG).

The manifest is the system of record (data_layer_and_tests.md §3.3): every
Parquet file of a release with its table, partitions, sha256, bytes, rows,
row groups and schema hash, and every artifact file (the bgzip sequence shards
and their indexes) with its sha256 and bytes. The catalog is disposable: a
DuckDB file holding only views over explicit file lists (no globs, LY-10) and
the sequence macros of §2A.3. It is built under a temporary name and moved into
place with rename(2); `current.json` is a symlink swapped the same way.

DuckDB cannot read BGZF, so `seq_raw` calls `bb_seq_raw`, a Python function that
`open_catalog` registers on each connection (DEC-17). A catalog opened with plain
`duckdb` has every view and macro, but the sequence macros fail with "bb_seq_raw
does not exist".

    manifest = build_manifest(root, "rS0.1")
    write_manifest(manifest, root / "manifests/release=rS0.1.json")
    generate_catalog(manifest, root / "catalog/release=rS0.1.duckdb")
    swap_current(root / "manifests", "rS0.1")
    con = open_catalog(root / "catalog/release=rS0.1.duckdb", local_dir="/tmp/if2453")
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import tempfile
from pathlib import Path

import blake3
import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from biobox import ids, layout
from biobox.registry.schemas import ARTIFACTS, TABLES, ArtifactSpec, TableSpec

MANIFEST_VERSION = 2

# Partitions whose value is derived from a file column (biobox.layout buckets): the column and
# how a value maps to its bucket key. A prefix of the sort key is checked from row-group
# min/max; the low byte of a token is not monotone in the sort order, so it is checked by
# reading the column's distinct values.
DERIVED_PARTITIONS = {
    "protein_prefix": ("protein_uid", "prefix"),
    "anchor_bucket": ("anchor_token", "low_byte"),
}


# Trees under the data root that carry their own manifests and are not release files: the
# search pack (search_pack/release=<R>/pack_manifest.json, verified by C-12) and search runs
# (search/run_index.json), which are written after the release is signed.
EXCLUDED_TREES = ("search_pack", "search")


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
        if any(len(x) != 2 for x in kv) or tuple(k for k, _ in kv) != spec.partition_keys:
            raise LayoutError(f"{rel}: partition path {dirs} does not match {spec.name} keys {spec.partition_keys} (LY-05)")
        return spec, dict(kv)
    return None


def classify_artifact(rel: str) -> tuple[ArtifactSpec, dict[str, str], str] | None:
    """Map a path under an artifact location to (artifact, partitions, companion suffix or "")."""
    for spec in ARTIFACTS.values():
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


def _artifact_entries(root: Path) -> list[dict]:
    out = []
    for spec in ARTIFACTS.values():
        base = root / spec.location
        if not base.exists():
            continue
        found = {}
        for path in sorted(p for p in base.rglob("*") if p.is_file()):
            rel = path.relative_to(root).as_posix()
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


def build_manifest(root: Path, release_id: str) -> dict:
    """Scan every .parquet file and artifact file under `root`; refuse strays, bad paths and schema drift."""
    root = Path(root).resolve()
    files = []
    for path in sorted(root.rglob("*.parquet")):
        rel = path.relative_to(root).as_posix()
        if rel.split("/", 1)[0] in EXCLUDED_TREES:
            continue
        hit = classify(rel)
        if hit is None:
            raise LayoutError(f"{rel}: no registered table claims this path (LY-10)")
        spec, parts = hit
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
    params_path = root / "l0seq/store_params.json"
    params = json.loads(params_path.read_text()) if params_path.exists() else {}
    return {
        "manifest_version": MANIFEST_VERSION,
        "release_id": release_id,
        "root": str(root),
        "params": params,
        "tables": {name: schema_hash(spec.schema) for name, spec in sorted(TABLES.items())},
        "files": files,
        "artifacts": _artifact_entries(root),
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


def write_manifest(manifest: dict, path: Path) -> None:
    _atomic_write(Path(path), json.dumps(manifest, sort_keys=True, indent=1).encode() + b"\n")


def read_manifest(path: Path) -> dict:
    return json.loads(Path(path).read_text())


def verify_manifest(manifest: dict, root: Path | None = None) -> list[str]:
    """LY-10: every listed file exists with its sha256 and size, and no unlisted .parquet or artifact file exists."""
    root = Path(root or manifest["root"])
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
        if (root / spec.location).exists():
            present += [p for p in (root / spec.location).rglob("*") if p.is_file()]
    for p in sorted(set(present)):
        rel = p.relative_to(root).as_posix()
        if rel not in listed:
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


def _view_sql(spec: TableSpec, paths: list[str], n_columns: int | None = None) -> str:
    if not paths:
        cols = [f'CAST(NULL AS {duck_type(f.type)}) AS "{f.name}"' for f in spec.schema]
        cols += [f'CAST(NULL AS VARCHAR) AS "{k}"' for k in spec.partition_keys]
        return f"CREATE VIEW {spec.name} AS SELECT {', '.join(cols)} WHERE false"
    files = "[" + ", ".join(_sql_str(p) for p in paths) + "]"
    if spec.partition_keys:
        types = "{" + ", ".join(f"{_sql_str(k)}: 'VARCHAR'" for k in spec.partition_keys) + "}"
        opts = f"hive_partitioning = true, hive_types = {types}"
    else:
        opts = "hive_partitioning = false"
    # union_by_name reads files written before a nullable column was added (LY-09); a column no
    # file has yet is a typed NULL
    n = len(spec.schema) if n_columns is None else n_columns
    cols = [f'"{f.name}"' if i < n else f'CAST(NULL AS {duck_type(f.type)}) AS "{f.name}"'
            for i, f in enumerate(spec.schema)] + [f'"{k}"' for k in spec.partition_keys]
    return (f"CREATE VIEW {spec.name} AS SELECT {', '.join(cols)} "
            f"FROM read_parquet({files}, {opts}, union_by_name = true)")


COMPLEMENT = ("ACGTRYKMBVDHNacgtrykmbvdhn", "TGCAYRMKVBHDNtgcayrmkvbhdn")


def macro_sql() -> list[str]:
    """Sequence macros (data_layer_and_tests.md §2A.3). 0-based half-open; strand +1/-1.

    Parameter names start with an underscore so they never collide with column names. Bounds
    are checked in SQL; `bb_seq_raw` (registered by `open_catalog`) only reads the shard.
    """
    return [
        f"CREATE MACRO revcomp(_x) AS reverse(translate(_x, '{COMPLEMENT[0]}', '{COMPLEMENT[1]}'))",
        "CREATE MACRO default_digest(_cuid) AS "
        "(SELECT raw_digest FROM contig_index WHERE contig_uid = _cuid AND is_default_variant)",
        "CREATE MACRO digest_len(_rd) AS (SELECT contig_len FROM seq_locator WHERE raw_digest = _rd)",
        "CREATE MACRO contig_len_of(_cuid) AS "
        "(SELECT contig_len FROM contig_index WHERE contig_uid = _cuid AND is_default_variant)",
        # forward bytes of one stored variant
        """CREATE MACRO seq_raw(_rd, _s, _e) AS CASE
              WHEN digest_len(_rd) IS NULL THEN error('seq_raw: unknown raw_digest')
              WHEN _s < 0 OR _e < _s OR _e > digest_len(_rd) THEN error('seq_raw: slice outside 0:contig_len')
              WHEN _s = _e THEN ''
              ELSE bb_seq_raw(_rd, _s, _e)
            END""",
        "CREATE MACRO seq_slice(_cuid, _s, _e) AS seq_raw(default_digest(_cuid), _s, _e)",
        # each macro below evaluates its slice once: a macro argument is substituted as an
        # expression, so passing seq_slice(...) to two CASE branches would read the shard twice
        """CREATE MACRO seq_get(_cuid, _s, _e, _strand) AS (
              SELECT CASE WHEN _strand = 1 THEN _x WHEN _strand = -1 THEN revcomp(_x)
                          ELSE error('seq_get: strand must be 1 or -1') END
              FROM (SELECT seq_slice(_cuid, _s, _e) AS _x))""",
        # circular contigs: an origin-wrapping slice has contig_len < _e <= _s + contig_len (C-05, DEC-13)
        """CREATE MACRO seq_get_circular(_cuid, _s, _e, _strand) AS (
              SELECT CASE
                WHEN _strand NOT IN (1, -1) THEN error('seq_get_circular: strand must be 1 or -1')
                WHEN _e > _n AND NOT (_s < _n AND _e <= _s + _n) THEN error('seq_get_circular: invalid wrapping slice')
                WHEN _strand = 1 THEN _x ELSE revcomp(_x) END
              FROM (SELECT _n, seq_slice(_cuid, _s, least(_e, _n)) || seq_slice(_cuid, 0, greatest(_e - _n, 0)) AS _x
                    FROM (SELECT contig_len_of(_cuid) AS _n)))""",
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
        # the anchor's bucket at every width a release may use (biobox.layout, biobox.l3.pairs)
        "CREATE MACRO pair_bucket_ok(_bucket, _tok) AS _bucket IN "
        "('all', left(printf('%02x', (_tok % 256)::INTEGER), 1), printf('%02x', (_tok % 256)::INTEGER))",
        # Q1: loci (anchor slots) of A with B within +-k slots; _tier NULL for every tier
        """CREATE MACRO q1_count(_a, _b, _k, _tier) AS (
              SELECT count(DISTINCT (contig_uid, anchor_ordinal)) FROM slot_pairs
              WHERE release = catalog_release() AND pair_bucket_ok(anchor_bucket, _a)
                AND anchor_token = _a AND partner_token = _b AND abs(orf_offset) <= _k
                AND (_tier IS NULL OR tier = _tier))""",
        # Q3 (counts; z_score ranking arrives with pair_stats): partners of A within +-k slots
        """CREATE MACRO partners_of(_a, _k) AS TABLE
              SELECT partner_token, count(DISTINCT (contig_uid, anchor_ordinal)) AS n_loci, count(*) AS n_pairs
              FROM slot_pairs
              WHERE release = catalog_release() AND pair_bucket_ok(anchor_bucket, _a)
                AND anchor_token = _a AND abs(orf_offset) <= _k
              GROUP BY partner_token ORDER BY n_loci DESC, partner_token""",
        # Q4: every instance of (A, B) within +-k slots, with geometry
        """CREATE MACRO pairs_of(_a, _b, _k) AS TABLE
              SELECT contig_uid, anchor_ordinal, partner_ordinal, orf_offset, junction_topology, delta_bp, flags, tier
              FROM slot_pairs
              WHERE release = catalog_release() AND pair_bucket_ok(anchor_bucket, _a)
                AND anchor_token = _a AND partner_token = _b AND abs(orf_offset) <= _k
              ORDER BY contig_uid, anchor_ordinal, partner_ordinal""",
        # Q2: tapes of a list of contigs with this release's tokens
        """CREATE MACRO tapes_of(_cuids) AS TABLE
              SELECT t.*, k.family_tokens, k.nr50_tokens, k.nr90_tokens
              FROM contig_tapes t LEFT JOIN tape_tokens k
                ON k.contig_uid = t.contig_uid AND k.gene_call_version = t.gene_call_version
                   AND k.release = catalog_release()
              WHERE list_contains(_cuids, t.contig_uid)
              ORDER BY t.contig_uid, t.gene_call_version""",
    ]


def generate_catalog(manifest: dict, out_path: Path, root: Path | None = None) -> Path:
    """Write the catalog under a temporary name and rename it into place (atomic on one filesystem)."""
    root = Path(root or manifest["root"]).resolve()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(f".{out_path.name}.tmp-{os.getpid()}")
    for p in (tmp, Path(str(tmp) + ".wal")):
        if p.exists():
            p.unlink()
    by_table: dict[str, list[str]] = {name: [] for name in TABLES}
    for f in manifest["files"]:
        by_table[f["table"]].append(str(root / f["path"]))
    con = duckdb.connect(str(tmp))
    try:
        widest = {}
        for f in manifest["files"]:
            widest[f["table"]] = max(widest.get(f["table"], 0), f.get("n_columns", len(TABLES[f["table"]].schema)))
        for name in sorted(TABLES):
            con.execute(_view_sql(TABLES[name], sorted(by_table[name]), widest.get(name)))
        con.execute("CREATE TABLE catalog_info (key VARCHAR, value VARCHAR)")
        con.executemany("INSERT INTO catalog_info VALUES (?, ?)", [
            ("release_id", manifest["release_id"]),
            ("manifest_hash", manifest_hash(manifest)),
            ("root", str(root)),
            ("duckdb_version", duckdb.__version__),
        ])
        con.execute("CREATE TABLE catalog_files AS SELECT * FROM (VALUES (NULL::VARCHAR, NULL::VARCHAR, NULL::VARCHAR, NULL::BIGINT, NULL::BIGINT)) "
                    "t(path, \"table\", sha256, \"rows\", row_groups) WHERE false")
        con.executemany("INSERT INTO catalog_files VALUES (?, ?, ?, ?, ?)",
                        [(f["path"], f["table"], f["sha256"], f["rows"], f["row_groups"]) for f in manifest["files"]])
        con.execute("CREATE TABLE catalog_artifacts (path VARCHAR, artifact VARCHAR, sha256 VARCHAR, bytes BIGINT)")
        arts = manifest.get("artifacts", [])
        if arts:
            con.executemany("INSERT INTO catalog_artifacts VALUES (?, ?, ?, ?)",
                            [(a["path"], a["artifact"], a["sha256"], a["bytes"]) for a in arts])
        # DuckDB binds a macro's functions when it is created: a stand-in lets the macros be
        # stored; open_catalog registers the real reader;
        # macros come last because they read catalog_info
        _register_seq_udf(con, None)
        for stmt in macro_sql():
            con.execute(stmt)
        con.execute("CHECKPOINT")
    finally:
        con.close()
    os.replace(tmp, out_path)
    return out_path


_READERS: dict[int, object] = {}


def _register_seq_udf(con: duckdb.DuckDBPyConnection, store) -> None:
    def bb_seq_raw(rd: pa.Array, s: pa.Array, e: pa.Array) -> pa.Array:
        if store is None:
            raise RuntimeError("bb_seq_raw: open the catalog with biobox.catalog.open_catalog to read sequence")
        out = [None if None in (r, a, b) else store.raw(r, a, b).decode("ascii")
               for r, a, b in zip(rd.to_pylist(), s.to_pylist(), e.to_pylist())]
        return pa.array(out, pa.string())

    con.create_function("bb_seq_raw", bb_seq_raw, ["BLOB", "BIGINT", "BIGINT"], "VARCHAR",
                        type="arrow", null_handling="special", side_effects=False)


def seq_reader(con: duckdb.DuckDBPyConnection):
    """The SeqStore behind a connection's `bb_seq_raw`, for its LY-06 counters."""
    return _READERS.get(id(con))


def open_catalog(path: Path, local_dir: Path | None = None, threads: int | None = None,
                 memory_limit: str | None = None, temp_directory: Path | None = None,
                 root: Path | None = None, cache_blocks: int = 64) -> duckdb.DuckDBPyConnection:
    """Open read-only, optionally from a copy on local disk (never open a DuckDB file over NFS).

    Registers `bb_seq_raw` over the shards under `root` (default: the root the catalog was
    built from, in `catalog_info`), with an LRU of `cache_blocks` decoded blocks per shard.
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

    store = SeqStore(root, cache_blocks) if (Path(root) / "l0seq/seq_locator.parquet").exists() else None
    _register_seq_udf(con, store)
    _READERS[id(con)] = store
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
    counts = {}
    if idx_path.exists():
        runs = json.loads(idx_path.read_text())["runs"]
        for name, fname in SEARCH_RUN_FILES.items():
            paths = [str(root / f"search/run={r['run_id']}/{fname}") for r in runs if fname in r["files"]]
            counts[name] = len(paths)
            if paths:
                spec = TABLES[name]
                cols = [f'"{f.name}"' for f in spec.schema] + ['"run"']
                files = "[" + ", ".join(_sql_str(x) for x in sorted(paths)) + "]"
                con.execute(f"CREATE OR REPLACE TEMP VIEW {name} AS SELECT {', '.join(cols)} FROM read_parquet("
                            f"{files}, hive_partitioning = true, hive_types = {{'run': 'VARCHAR'}}, union_by_name = true)")
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
    `bb_seq_raw` read from the shards (blocks decompressed, compressed bytes read, bytes
    decoded), zero when the query reads no sequence. Counts, not time, decide layouts (LY-06).
    """
    store = seq_reader(con)
    before = store.counters() if store is not None else None
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
    bgzf = ({k: v - before[k] for k, v in store.counters().items()} if store is not None
            else {"blocks_read": 0, "compressed_bytes_read": 0, "uncompressed_bytes_decoded": 0})
    return {"result": result, "scans": scans, "bytes_read": info.get("io", {}).get("total_bytes_read"), "bgzf": bgzf}
