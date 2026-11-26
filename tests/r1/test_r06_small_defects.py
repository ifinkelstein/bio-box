"""Ring R1, card R-06: four small defects of the audit of 2026-11-25.

(a) `export.duck_memory_limit` reads a cgroup v1 limit too; (b) `search_mem_gb` refuses a `$MEM` it cannot parse;
(c) `panel._resolve_library_path` orders library versions numerically; (d) the comment on `MIN_SHARD_BASES`.
No test reads the host's cgroup, memory or thread count.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

from biobox.l0seq import store
from biobox.search import export, panel

pytestmark = pytest.mark.covers("W-42", "A-25", "W-39", "B-73")

GIB = 2**30
REPO = Path(__file__).resolve().parents[2]
LAUNCH = REPO / "examples" / "_shared" / "launch.sh"


# (a) cgroup v1 ------------------------------------------------------------------------------------------------------

def _v1_world(tmp_path: Path, value: str, rel: str = "/job") -> tuple[Path, Path]:
    proc = tmp_path / "proc"
    (proc / "self").mkdir(parents=True)
    (proc / "self" / "cgroup").write_text(f"4:memory:{rel}\n3:cpu,cpuacct:{rel}\n")
    cg = tmp_path / "cg"
    d = cg / "memory" / rel.strip("/")
    d.mkdir(parents=True)
    (d / "memory.limit_in_bytes").write_text(value + "\n")
    (cg / "memory" / "memory.limit_in_bytes").write_text("9223372036854771712\n")
    return proc, cg


def test_cgroup_v1_limit_is_read(tmp_path):
    proc, cg = _v1_world(tmp_path, str(4 * GIB))
    assert export._cgroup_memory_bytes(proc_root=proc, cgroup_root=cg) == 4 * GIB


def test_cgroup_v1_limit_clamps_duck_memory(tmp_path):
    proc, cg = _v1_world(tmp_path, str(2 * GIB))
    assert export.duck_memory_limit(40, proc_root=proc, cgroup_root=cg) == "1GB"
    proc, cg = _v1_world(tmp_path / "b", str(4 * GIB))
    assert export.duck_memory_limit(40, proc_root=proc, cgroup_root=cg) == "2GB"


def test_cgroup_v1_unlimited_is_no_limit(tmp_path):
    proc, cg = _v1_world(tmp_path, "9223372036854771712")
    assert export._cgroup_memory_bytes(proc_root=proc, cgroup_root=cg) is None
    assert export.duck_memory_limit(40, proc_root=proc, cgroup_root=cg) == "16GB"


def test_no_cgroup_files_is_no_limit(tmp_path):
    assert export._cgroup_memory_bytes(proc_root=tmp_path / "none", cgroup_root=tmp_path / "none") is None


# (b) search_mem_gb --------------------------------------------------------------------------------------------------

def _mem_gb(env_extra: dict[str, str]) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if k not in ("MEM", "SEARCH_MEM_GB")}
    env.update(env_extra)
    return subprocess.run(["bash", "-c", f'. "{LAUNCH}" && search_mem_gb'], env=env, capture_output=True,
                          text=True, timeout=30)


@pytest.mark.parametrize("mem", ["1.5T", "lots", "64", "G", "-4G", "64K", "12 G", "1.5G"])
def test_search_mem_gb_refuses_unparseable_mem(mem):
    got = _mem_gb({"MEM": mem})
    assert got.returncode != 0, got.stdout
    assert got.stdout.strip() == ""
    assert repr(mem) in got.stderr or mem in got.stderr, got.stderr


@pytest.mark.parametrize("mem,expected", [("64G", "62"), ("750G", "748"), ("32000M", "29"), ("16G", "14"),
                                          ("1T", "1022"), ("64GiB", "62")])
def test_search_mem_gb_accepted_values_unchanged(mem, expected):
    got = _mem_gb({"MEM": mem})
    assert got.returncode == 0, got.stderr
    assert got.stdout.strip() == expected


def test_search_mem_gb_override_wins_even_over_bad_mem():
    got = _mem_gb({"MEM": "1.5T", "SEARCH_MEM_GB": "48"})
    assert got.returncode == 0, got.stderr
    assert got.stdout.strip() == "48"


def test_search_mem_gb_failure_stops_a_set_e_caller():
    script = f'set -euo pipefail\n. "{LAUNCH}"\nX=$(search_mem_gb)\necho started $X\n'
    env = {k: v for k, v in os.environ.items() if k != "SEARCH_MEM_GB"}
    env["MEM"] = "1.5T"
    got = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True, timeout=30)
    assert got.returncode != 0
    assert "started" not in got.stdout
    assert "1.5T" in got.stderr


# (c) library versions -----------------------------------------------------------------------------------------------

def _pack(tmp_path: Path, names: list[str]) -> Path:
    m = tmp_path / "pack" / "models"
    m.mkdir(parents=True)
    for n in names:
        (m / n).write_text("x")
    return tmp_path / "pack"


def test_resolve_library_orders_versions_numerically(tmp_path):
    pack = _pack(tmp_path, ["pfam-38.2.hmm.gz", "pfam-38.10.hmm.gz", "pfam-37.4.hmm.gz"])
    got = panel._resolve_library_path("pfam", pack, None, [])
    assert got.name == "pfam-38.10.hmm.gz"


def test_resolve_library_unparseable_version_does_not_win(tmp_path):
    pack = _pack(tmp_path, ["pfam-38.2.hmm.gz", "pfam-zzz.hmm.gz", "pfam-old.hmm.gz"])
    assert panel._resolve_library_path("pfam", pack, None, []).name == "pfam-38.2.hmm.gz"


def test_resolve_library_refdb_orders_versions_numerically(tmp_path):
    base = tmp_path / "refdb"
    d = base / "pfam"
    d.mkdir(parents=True)
    for n in ("pfam-38.2.hmm", "pfam-38.10.hmm"):
        (d / n).write_text("x")
    assert panel._resolve_library_path("pfam", None, None, [base]).name == "pfam-38.10.hmm"


def test_resolve_library_single_unversioned_and_missing(tmp_path):
    pack = _pack(tmp_path, ["ncbifam.hmm.gz"])
    assert panel._resolve_library_path("ncbifam", pack, None, []).name == "ncbifam.hmm.gz"
    with pytest.raises(panel.PanelError, match="not found"):
        panel._resolve_library_path("nothing", pack, None, [])


# (d) the comment on MIN_SHARD_BASES ---------------------------------------------------------------------------------

def test_min_shard_bases_comment_states_the_exact_bound():
    assert store.MIN_SHARD_BASES == 320 * 2**20
    assert store.MIN_SHARD_BASES * 20 == 64 * 2**20 * 100  # exactly 64 MiB at 0.20 B/bp, not above it
    line = next(ln for ln in (REPO / "src" / "biobox" / "l0seq" / "store.py").read_text().splitlines()
                if ln.startswith("MIN_SHARD_BASES ="))
    assert "> 64 MiB" not in line
    assert "64 MiB" in line and "0.20" in line
