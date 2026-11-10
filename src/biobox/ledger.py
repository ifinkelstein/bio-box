"""Step ledgers and their reconciliation (data_layer_and_tests.md §4.2: L-01, L-02, L-03).

Each step writes `<root>/ledger/stage=<S>/step=<name>/part-0.parquet`: one row per
(object, reason) with n_in, n_kept, n_merged and n_rejected for that reason. A release's own steps
(`RELEASE_STEPS`: clustering) write `<root>/ledger/release=<R>/step=<name>/part-0.parquet` instead
(`write_release`, card B-28), so two releases over one stage (rS0.1, rS0.2) each have their own.

S0 steps and what their counts mean:

    ingest      per L0 object: records in; kept = byte sequences stored new; merged = records
                whose bytes were already stored (an occurrence row only)
    genecall    per contig: genes in; kept = genes whose protein is the first of its sequence in
                the stage; merged = genes repeating an earlier protein
    domains-<lib>  per library: distinct proteins in; kept = with >= 1 included domain;
                rejected (no_hit) = the rest
    cluster     per level: proteins in; kept = representatives; merged = other members

The DAG edges L-03 checks are `CHAIN`: a downstream n_in equals the upstream n_kept.

Bounded memory (card B-19; AGENTS §4.1): `write` streams any iterable of rows, which must arrive
sorted by (object_key, reason_code), through `biobox.io.write_table` in batches of `BATCH_ROWS`;
a step's ledger is written once (append-only). The checks are DuckDB aggregates over the ledger
files, so no ledger is loaded; only `read` returns a whole ledger as a list.
"""

from __future__ import annotations

import contextlib
import json
import shutil
import tempfile
from collections.abc import Iterable, Iterator
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from biobox import io
from biobox.registry.schemas import table

# (upstream step, upstream object filter, downstream step, downstream object filter)
CHAIN = [
    ("genecall", None, "domains-*", None),      # distinct proteins scanned, one step per library
    ("genecall", None, "cluster", b"nr90"),     # distinct proteins clustered at NR90
    ("cluster", b"nr90", "cluster", b"nr50"),   # NR50 runs on the NR90 representatives
]
BATCH_ROWS = 65_536  # ledger rows handed to the writer at a time
_DUCK_MEMORY_LIMIT = "1GB"  # the checks' DuckDB cap; larger aggregates spill to a private temporary directory


def row(object_key: bytes, kind: str, n_in: int, n_kept: int, n_merged: int = 0, n_rejected: int = 0,
        reason: str = "ok", tools: dict[str, str] | None = None) -> dict:
    """`tools`: versions of the tools that made this step's output (C-08, `provenance.tool_versions`)."""
    return {"object_key": object_key, "object_kind": kind, "n_in": n_in, "n_kept": n_kept,
            "n_merged": n_merged, "n_rejected": n_rejected, "reason_code": reason,
            "tool_versions_json": json.dumps(tools, sort_keys=True) if tools is not None else None}


RELEASE_STEPS = frozenset({"cluster"})  # steps of a release, not of a stage (card B-28)


def path(root: Path, stage: str, step: str) -> Path:
    return Path(root) / f"ledger/stage={stage}/step={step}/part-0.parquet"


def release_path(root: Path, release: str, step: str) -> Path:
    """The ledger of a release's own step (`RELEASE_STEPS`)."""
    return Path(root) / f"ledger/release={release}/step={step}/part-0.parquet"


def _step_path(root: Path, stage: str, step: str, release: str | None) -> Path | None:
    """Where `step`'s ledger is for `stage`: a release step's is its release's (None without a release)."""
    if step in RELEASE_STEPS:
        return None if release is None else release_path(root, release, step)
    return path(root, stage, step)


def write(root: Path, stage: str, step: str, rows: Iterable[dict]) -> Path:
    """Write the step's ledger from `rows`, any iterable consumed once and never materialized.

    Rows must come sorted by (object_key, reason_code) (the writer refuses them otherwise; callers
    `ORDER BY` in SQL or sort a handful of rows). A ledger that exists is never rewritten. A release
    step (`RELEASE_STEPS`) is refused here: it goes to `write_release`.
    """
    if step in RELEASE_STEPS:
        raise ValueError(f"{step} is a release step: its ledger is the release's (ledger.write_release)")
    return _write(path(root, stage, step), rows)


def write_release(root: Path, release: str, step: str, rows: Iterable[dict]) -> Path:
    """As `write`, for a release's own step (`RELEASE_STEPS`), under `ledger/release=<R>/`."""
    if step not in RELEASE_STEPS:
        raise ValueError(f"{step} is not a release step {sorted(RELEASE_STEPS)}: its ledger is a stage's (ledger.write)")
    return _write(release_path(root, release, step), rows)


def _write(out: Path, rows: Iterable[dict]) -> Path:
    spec = table("ledger")

    def batches() -> Iterator[pa.RecordBatch]:
        buf: list[dict] = []
        for r in rows:
            buf.append(r)
            if len(buf) >= BATCH_ROWS:
                yield pa.RecordBatch.from_pylist(buf, schema=spec.schema)
                buf = []
        if buf:
            yield pa.RecordBatch.from_pylist(buf, schema=spec.schema)

    io.write_table(batches(), spec, out)
    return out


def read(root: Path, stage: str, step: str) -> list[dict]:
    """The step's ledger rows as a list (tests and small ledgers; the checks do not use it)."""
    p = path(root, stage, step)
    return pq.read_table(p).to_pylist() if p.exists() else []


def read_release(root: Path, release: str, step: str) -> list[dict]:
    """As `read`, for a release's own step."""
    p = release_path(root, release, step)
    return pq.read_table(p).to_pylist() if p.exists() else []


def _one_stage_release(root: Path, stage: str, release: str | None) -> None:
    """A release step's counts follow from one stage's only when the release is that stage alone."""
    if release is None:
        return
    from biobox import catalog

    stages = catalog.read_release_def(root, release)
    if stages is not None and stages != [stage]:
        raise ValueError(f"release {release} has stages {stages}: its steps are checked against stage {stage} "
                         "only when it is that stage alone")


@contextlib.contextmanager
def _duck() -> Iterator[duckdb.DuckDBPyConnection]:
    scratch = Path(tempfile.mkdtemp(prefix="biobox-ledger-"))
    con = duckdb.connect(config={"threads": 1, "memory_limit": _DUCK_MEMORY_LIMIT,
                                 "temp_directory": str(scratch / "duckdb")})
    try:
        yield con
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)


def _src(files: list[Path]) -> str:
    """The ledger rows of `files` (each with `si`, its position in the list) as a subquery."""
    return "(" + " UNION ALL ".join(
        f"SELECT {i} AS si, * FROM read_parquet('{str(f).replace(chr(39), chr(39) * 2)}', hive_partitioning = false)"
        for i, f in enumerate(files)) + ")"


def check_balance(root: Path, stage: str, steps: list[str], release: str | None = None) -> list[str]:
    """L-01: per object, n_in = n_kept + n_merged + sum of n_rejected over reasons. A release step
    (`RELEASE_STEPS`) is read from `release`'s ledger, and skipped without one."""
    problems = []
    with _duck() as con:
        for step in steps:
            p = _step_path(root, stage, step, release)
            if p is None or not p.exists():
                continue
            for uid, same, n_in, kept, merged, rejected in con.execute(
                    "SELECT object_key, count(DISTINCT n_in) = 1 AND count(DISTINCT n_kept) = 1 "
                    "AND count(DISTINCT n_merged) = 1, min(n_in), min(n_kept), min(n_merged), sum(n_rejected) "
                    f"FROM {_src([p])} GROUP BY object_key "
                    "HAVING NOT (count(DISTINCT n_in) = 1 AND count(DISTINCT n_kept) = 1 AND count(DISTINCT n_merged) = 1) "
                    "OR min(n_in)::HUGEINT <> min(n_kept)::HUGEINT + min(n_merged)::HUGEINT + sum(n_rejected) "
                    "ORDER BY object_key").fetchall():
                uid = bytes(uid)
                if not same:
                    problems.append(f"{step} {uid!r}: n_in/n_kept/n_merged differ between reason rows")
                else:
                    problems.append(f"{step} {uid!r}: {n_in} in != {kept} kept + {merged} merged + {rejected} rejected")
    return problems


def check_inputs(root: Path, stage: str, expected_objects: set[bytes]) -> list[str]:
    """L-02: the first ledger lists exactly the stage's L0 objects."""
    p = path(root, stage, "ingest")
    expected = pa.table({"k": pa.array(sorted(expected_objects), pa.binary())})
    with _duck() as con:
        con.register("expected", expected)
        got = f"(SELECT DISTINCT object_key AS k FROM {_src([p])})" if p.exists() else "(SELECT NULL::BLOB AS k WHERE false)"
        (missing,), = con.execute(f"SELECT count(*) FROM expected e ANTI JOIN {got} g ON g.k = e.k").fetchall()
        (extra,), = con.execute(f"SELECT count(*) FROM {got} g ANTI JOIN expected e ON g.k = e.k").fetchall()
    problems = []
    if missing:
        problems.append(f"ingest ledger misses {missing} L0 objects")
    if extra:
        problems.append(f"ingest ledger has {extra} objects not in the input")
    return problems


def check_chain(root: Path, stage: str, release: str | None = None) -> list[str]:
    """L-03: along each DAG edge, downstream n_in equals upstream n_kept (summed over the filter).

    An object's count is taken from its rows (equal across its reason rows when L-01 holds; the
    smallest otherwise). Problems come in the order of the downstream objects' first appearance:
    steps in name order, objects in key order within a step. The release steps (`RELEASE_STEPS`) are
    read from `release`'s ledger (card B-28), a release whose `release_def` is `stage` alone
    (ValueError otherwise); without `release` their edges are not checked.
    """
    _one_stage_release(root, stage, release)
    problems = []
    with _duck() as con:
        for up, up_obj, down, down_obj in CHAIN:
            up_path = _step_path(root, stage, up, release)
            if up_path is None or not up_path.exists():
                continue
            (rows_up,), = con.execute(f"SELECT count(*) FROM {_src([up_path])}").fetchall()
            if not rows_up:
                continue
            filt = "" if up_obj is None else "WHERE object_key = $obj"
            (kept,), = con.execute(
                f"SELECT coalesce(sum(k), 0) FROM (SELECT object_key, min(n_kept)::HUGEINT AS k "
                f"FROM {_src([up_path])} GROUP BY object_key) {filt}",
                {} if up_obj is None else {"obj": up_obj}).fetchall()
            downs = ([p.name.removeprefix("step=") for p in sorted((Path(root) / f"ledger/stage={stage}").glob(f"step={down}"))]
                     if down.endswith("*") else [down])
            files = [f for f in (_step_path(root, stage, d, release) for d in downs) if f is not None and f.exists()]
            if not files:
                continue
            dfilt = "" if down_obj is None else "HAVING object_key = $dobj"
            params = {"kept": kept, **({} if down_obj is None else {"dobj": down_obj})}
            for u, n_in in con.execute(
                    f"SELECT object_key, n_in FROM (SELECT object_key, min(n_in)::HUGEINT AS n_in, min(si) AS first "
                    f"FROM {_src(files)} GROUP BY object_key {dfilt}) WHERE n_in <> $kept::HUGEINT "
                    "ORDER BY first, object_key", params).fetchall():
                problems.append(f"{up}->{down} {bytes(u)!r}: n_in {n_in} != upstream n_kept {kept}")
    return problems
