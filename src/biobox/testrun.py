"""`test_run` rows (data_layer_and_tests.md §4.12; TEST_ARCHITECTURE §8).

One row per (test execution, plan test ID). A pytest test names the IDs it covers in its
function name (`test_r03_...`, `test_k02_k03_...`, `test_ly06_...`, `test_mr14_...`) or
with `@pytest.mark.covers("C-01", "LY-01")`; a test that covers no plan ID gets its node
ID as `test_id`, so every execution is still recorded.

Rows go to `<dir>/test_run/<yyyy-mm-dd>-<k>.parquet`, one file per pytest session, and
`<dir>/test_run_index.json` is rewritten by atomic rename to list the files.
"""

from __future__ import annotations

import json
import os
import re
import socket
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import blake3
import pyarrow as pa
import pyarrow.parquet as pq

from biobox.registry.schemas import table

_PREFIX = {"r": "R", "ly": "LY", "k": "K", "mr": "MR", "c": "C", "l": "L", "d": "D", "g": "G", "a": "A", "s": "S"}
_TOKEN = re.compile(r"(ly|mr|[rkcldgas])(\d{2})([a-z]?)$")


def ids_from_name(name: str) -> list[str]:
    """Plan test IDs at the start of a test function name: test_k02_k03_x -> [K-02, K-03]."""
    out = []
    for tok in name.removeprefix("test_").split("_"):
        m = _TOKEN.match(tok)
        if not m:
            break
        out.append(f"{_PREFIX[m.group(1)]}-{m.group(2)}{m.group(3)}")
    return out


def code_commit(repo: Path | None = None) -> str:
    if os.environ.get("BIOBOX_CODE_COMMIT"):
        return os.environ["BIOBOX_CODE_COMMIT"]
    try:
        out = subprocess.run(["git", "-C", str(repo or Path.cwd()), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=10)
        return out.stdout.strip() or "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def params_hash(nodeid: str) -> str:
    return blake3.blake3(nodeid.encode()).hexdigest(length=16)


def make_row(test_id: str, ring: str, outcome: str, duration_s: float, started_at: datetime, nodeid: str,
             metrics: dict | None = None, operation: str | None = None, workflow: str | None = None,
             commit: str = "unknown", seed: int | None = None) -> dict:
    return {
        "test_id": test_id,
        "ring": ring,
        "operation": operation,
        "workflow": workflow,
        "params_hash": params_hash(nodeid),
        "outcome": outcome,
        "metrics_json": json.dumps({"nodeid": nodeid, **(metrics or {})}, sort_keys=True, default=str),
        "duration_s": float(duration_s),
        "started_at": started_at,
        "code_commit": commit,
        "runtime_image_digest": os.environ.get("BIOBOX_IMAGE_DIGEST") or None,
        "test_image_digest": os.environ.get("BIOBOX_TEST_IMAGE_DIGEST") or None,
        "release_id": os.environ.get("BIOBOX_RELEASE_ID") or None,
        "host": socket.gethostname(),
        "seed": seed,
    }


def write_session(rows: list[dict], out_dir: Path, day: datetime | None = None) -> Path:
    """Write one session's rows as the next `<date>-<k>.parquet` and refresh the index."""
    out_dir = Path(out_dir)
    d = out_dir / "test_run"
    d.mkdir(parents=True, exist_ok=True)
    day = (day or datetime.now(timezone.utc)).strftime("%Y-%m-%d")
    k = 0
    while (d / f"{day}-{k}.parquet").exists():
        k += 1
    path = d / f"{day}-{k}.parquet"
    spec = table("test_run")
    t = pa.Table.from_pylist(rows, schema=spec.schema).sort_by([("started_at", "ascending"), ("test_id", "ascending")])
    tmp = d / f".{path.name}.tmp-{os.getpid()}"
    pq.write_table(t, tmp, compression="zstd")
    os.replace(tmp, path)
    index = sorted(p.relative_to(out_dir).as_posix() for p in d.glob("*.parquet"))
    tmp_idx = out_dir / f".test_run_index.json.tmp-{os.getpid()}"
    tmp_idx.write_text(json.dumps({"files": index}, indent=1) + "\n")
    os.replace(tmp_idx, out_dir / "test_run_index.json")
    return path
