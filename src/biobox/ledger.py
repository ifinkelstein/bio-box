"""Step ledgers and their reconciliation (data_layer_and_tests.md §4.2: L-01, L-02, L-03).

Each step writes `<root>/ledger/stage=<S>/step=<name>/part-0.parquet`: one row per
(object, reason) with n_in, n_kept, n_merged and n_rejected for that reason.

S0 steps and what their counts mean:

    ingest      per L0 object: records in; kept = byte sequences stored new; merged = records
                whose bytes were already stored (an occurrence row only)
    genecall    per contig: genes in; kept = genes whose protein is the first of its sequence in
                the stage; merged = genes repeating an earlier protein
    domains-<lib>  per library: distinct proteins in; kept = with >= 1 included domain;
                rejected (no_hit) = the rest
    cluster     per level: proteins in; kept = representatives; merged = other members

The DAG edges L-03 checks are `CHAIN`: a downstream n_in equals the upstream n_kept.
"""

from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from biobox.registry.schemas import table

# (upstream step, upstream object filter, downstream step, downstream object filter)
CHAIN = [
    ("genecall", None, "domains-*", None),      # distinct proteins scanned, one step per library
    ("genecall", None, "cluster", b"nr90"),     # distinct proteins clustered at NR90
    ("cluster", b"nr90", "cluster", b"nr50"),   # NR50 runs on the NR90 representatives
]


def row(object_key: bytes, kind: str, n_in: int, n_kept: int, n_merged: int = 0, n_rejected: int = 0,
        reason: str = "ok") -> dict:
    return {"object_key": object_key, "object_kind": kind, "n_in": n_in, "n_kept": n_kept,
            "n_merged": n_merged, "n_rejected": n_rejected, "reason_code": reason}


def write(root: Path, stage: str, step: str, rows: list[dict]) -> Path:
    path = Path(root) / f"ledger/stage={stage}/step={step}/part-0.parquet"
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(rows, key=lambda r: (r["object_key"], r["reason_code"]))
    pq.write_table(pa.Table.from_pylist(rows, schema=table("ledger").schema), path, compression="zstd")
    return path


def read(root: Path, stage: str, step: str) -> list[dict]:
    path = Path(root) / f"ledger/stage={stage}/step={step}/part-0.parquet"
    return pq.read_table(path).to_pylist() if path.exists() else []


def _per_object(rows: list[dict]) -> dict[bytes, dict]:
    out: dict[bytes, dict] = {}
    for r in rows:
        o = out.setdefault(r["object_key"], {"n_in": set(), "n_kept": set(), "n_merged": set(), "rejected": 0})
        for k in ("n_in", "n_kept", "n_merged"):
            o[k].add(r[k])
        o["rejected"] += r["n_rejected"]
    return out


def check_balance(root: Path, stage: str, steps: list[str]) -> list[str]:
    """L-01: per object, n_in = n_kept + n_merged + sum of n_rejected over reasons."""
    problems = []
    for step in steps:
        for uid, o in _per_object(read(root, stage, step)).items():
            if len(o["n_in"]) != 1 or len(o["n_kept"]) != 1 or len(o["n_merged"]) != 1:
                problems.append(f"{step} {uid!r}: n_in/n_kept/n_merged differ between reason rows")
                continue
            (n_in,), (kept,), (merged,) = o["n_in"], o["n_kept"], o["n_merged"]
            if n_in != kept + merged + o["rejected"]:
                problems.append(f"{step} {uid!r}: {n_in} in != {kept} kept + {merged} merged + {o['rejected']} rejected")
    return problems


def check_inputs(root: Path, stage: str, expected_objects: set[bytes]) -> list[str]:
    """L-02: the first ledger lists exactly the stage's L0 objects."""
    got = {r["object_key"] for r in read(root, stage, "ingest")}
    problems = []
    if expected_objects - got:
        problems.append(f"ingest ledger misses {len(expected_objects - got)} L0 objects")
    if got - expected_objects:
        problems.append(f"ingest ledger has {len(got - expected_objects)} objects not in the input")
    return problems


def check_chain(root: Path, stage: str) -> list[str]:
    """L-03: along each DAG edge, downstream n_in equals upstream n_kept (summed over the filter)."""
    problems = []
    for up, up_obj, down, down_obj in CHAIN:
        urows = read(root, stage, up)
        if not urows:
            continue
        kept = sum(next(iter(o["n_kept"])) for u, o in _per_object(urows).items() if up_obj is None or u == up_obj)
        downs = ([p.name.removeprefix("step=") for p in sorted((Path(root) / f"ledger/stage={stage}").glob(f"step={down}"))]
                 if down.endswith("*") else [down])
        drows = [r for d in downs for r in read(root, stage, d)]
        for u, o in _per_object(drows).items():
            if down_obj is not None and u != down_obj:
                continue
            n_in = next(iter(o["n_in"]))
            if n_in != kept:
                problems.append(f"{up}->{down} {u!r}: n_in {n_in} != upstream n_kept {kept}")
    return problems
