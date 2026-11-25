"""Ring R1, card W-42: DuckDB memory_limit scales with threads in `export.open_store` and `arrays._connect`
without exceeding half of an active cgroup v2 `memory.max` cap.
"""

from __future__ import annotations

import re
import shutil

import pytest

from biobox.search import arrays, export

pytestmark = pytest.mark.covers("W-42", "A-25", "W-05")

GIB = 2**30
_SIZE_RE = re.compile(r"^\s*([0-9]+(?:\.[0-9]+)?)\s*([A-Za-z]+)\s*$")
_UNITS = {
    "B": 1,
    "KB": 10**3, "KIB": 2**10,
    "MB": 10**6, "MIB": 2**20,
    "GB": 10**9, "GIB": 2**30,
    "TB": 10**12, "TIB": 2**40,
}


def _duck_bytes(con) -> float:
    raw = str(con.execute("SELECT current_setting('memory_limit')").fetchone()[0])
    m = _SIZE_RE.match(raw)
    assert m is not None, f"unrecognised DuckDB memory_limit format: {raw!r}"
    return float(m.group(1)) * _UNITS[m.group(2).upper()]


@pytest.mark.parametrize("threads, expected", [
    (0, "2GB"),
    (1, "2GB"),
    (2, "2GB"),
    (4, "2GB"),
    (8, "4GB"),
    (16, "8GB"),
    (24, "12GB"),
    (32, "16GB"),
    (40, "16GB"),
])
def test_duck_memory_limit_scales_with_threads_without_cgroup(threads: int, expected: str):
    assert export.duck_memory_limit(threads, cgroup_bytes=None) == expected
    assert export.duck_memory_limit(threads, cgroup_bytes=64 * GIB) == expected


@pytest.mark.parametrize("threads", [1, 2, 4, 8, 16, 24, 40])
def test_duck_memory_limit_clamps_to_half_of_tight_cgroup(threads: int):
    # 4 GiB cgroup -> 50% cap is 2 GB: matches default at <=4 threads and never exceeds 2GB
    assert export.duck_memory_limit(threads, cgroup_bytes=4 * GIB) == "2GB"
    # 2 GiB cgroup -> 50% cap is 1 GB: prevents OOM in tight test cgroups
    assert export.duck_memory_limit(threads, cgroup_bytes=2 * GIB) == "1GB"


def test_open_store_scales_memory_limit_with_threads_and_supports_override(tmp_path):
    manifest = {"release_id": "r-w42", "stages": [], "files": []}
    with export.open_store(tmp_path, release="r-w42", gene_call_version="gcv-1", manifest=manifest,
                           threads=4, tmp=tmp_path, relations=()) as st:
        mem4 = _duck_bytes(st.con)
    assert 1.5 * GIB <= mem4 <= 2.2 * GIB

    with export.open_store(tmp_path, release="r-w42", gene_call_version="gcv-1", manifest=manifest,
                           threads=24, tmp=tmp_path, relations=()) as st:
        mem24 = _duck_bytes(st.con)
    assert mem24 >= 10 * GIB

    with export.open_store(tmp_path, release="r-w42", gene_call_version="gcv-1", manifest=manifest,
                           threads=24, memory_limit="3GB", tmp=tmp_path, relations=()) as st:
        mem_override = _duck_bytes(st.con)
    assert 2.5 * GIB <= mem_override <= 3.2 * GIB


def test_arrays_connect_scales_memory_limit_with_threads(tmp_path):
    con24, scratch24 = arrays._connect(tmp_path, threads=24)
    try:
        mem24 = _duck_bytes(con24)
    finally:
        con24.close()
        shutil.rmtree(scratch24, ignore_errors=True)
    assert mem24 >= 10 * GIB

    con_def, scratch_def = arrays._connect(tmp_path, threads=None)
    try:
        mem_def = _duck_bytes(con_def)
    finally:
        con_def.close()
        shutil.rmtree(scratch_def, ignore_errors=True)
    assert 1.5 * GIB <= mem_def <= 2.2 * GIB
