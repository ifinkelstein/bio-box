"""C-23: a release records the code that built it (ring R2 by directory).

The code identity (`biobox.codeid`) is a sha256 over the package's own files, so it is the same for the same source
at any path or host; run records carry it; `catalog_info.built_by` reports it; and the MANIFEST does not.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import duckdb
import pytest

import biobox
from biobox import catalog, codeid

pytestmark = pytest.mark.r2

PKG = Path(biobox.__file__).resolve().parent

# The manifest hash and the sha256 over every manifest-listed file of the micro-S0 fixture, measured on the pod on
# the tree before card C-23 (the same fixture build as tests/r2/test_micro_s0_catalog.py): this card must not move them.
BASELINE_MANIFEST_HASH = "e42fea7f024ff237a6cfa1102ff54493"
BASELINE_FILES_SHA256 = "8ff62c46a63ee6246bec2c0ca2379feb14d3cd26a4d29b87210dc8abc22d6e13"

BUILD = r"""
import hashlib, json, sys
from pathlib import Path
sys.path.insert(0, sys.argv[2])
import biobox
assert Path(biobox.__file__).resolve().is_relative_to(Path(sys.argv[3]).resolve()), biobox.__file__
from test_micro_s0_catalog import _build
from biobox import catalog
root = _build(Path(sys.argv[1]))
man = catalog.build_manifest(root, "rS0.test", ["S0"])
h = hashlib.sha256()
for f in man["files"] + man["artifacts"]:
    h.update(f["path"].encode() + (root / f["path"]).read_bytes())
print("RESULT", catalog.manifest_hash(man), h.hexdigest())
"""


def _copy_pkg(dst: Path) -> Path:
    shutil.copytree(PKG, dst / "biobox", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    return dst / "biobox"


# --- the identity ---------------------------------------------------------------------


def test_code_id_is_a_function_of_the_source_not_of_where_it_is(tmp_path):
    a = _copy_pkg(tmp_path / "one")
    b = _copy_pkg(tmp_path / "two" / "deeper")
    ia, ib = codeid.code_identity(a), codeid.code_identity(b)
    assert ia["code_id"] == ib["code_id"] == codeid.code_identity(PKG)["code_id"]
    assert len(ia["code_id"]) == 64 and int(ia["code_id"], 16) >= 0
    assert ia["code_commit"] is None  # a copy under tmp is no work tree


def test_code_id_changes_with_any_source_byte_or_file(tmp_path):
    a = _copy_pkg(tmp_path / "one")
    base = codeid.code_identity(a)["code_id"]
    target = a / "ids.py"
    original = target.read_bytes()
    target.write_bytes(original + b"\n")
    changed = codeid.code_identity(a)["code_id"]
    assert changed != base
    target.write_bytes(original)
    assert codeid.code_identity(a)["code_id"] == base
    (a / "extra_data.txt").write_bytes(b"x")  # a data file counts, as does a new name
    with_data = codeid.code_identity(a)["code_id"]
    assert with_data != base
    (a / "extra_data.txt").rename(a / "other_name.txt")
    assert codeid.code_identity(a)["code_id"] not in (base, with_data)


def test_code_id_ignores_bytecode_caches(tmp_path):
    a = _copy_pkg(tmp_path / "one")
    base = codeid.code_identity(a)["code_id"]
    (a / "__pycache__").mkdir(exist_ok=True)
    (a / "__pycache__" / "ids.cpython-312.pyc").write_bytes(b"junk")
    assert codeid.code_identity(a)["code_id"] == base


def test_commit_from_a_git_that_answers(tmp_path, monkeypatch):
    """A stand-in `git` (the image has none): inside a work tree it prints `true` and a commit."""
    a = _copy_pkg(tmp_path / "one")
    bin_ = tmp_path / "bin"
    bin_.mkdir()
    (bin_ / "git").write_text("#!/bin/sh\nprintf 'true\\n0123456789abcdef0123456789abcdef01234567\\n'\n")
    (bin_ / "git").chmod(0o755)
    monkeypatch.setenv("PATH", str(bin_))
    assert codeid.code_identity(a)["code_commit"] == "0123456789abcdef0123456789abcdef01234567"
    (bin_ / "git").write_text("#!/bin/sh\necho 'fatal: not a git repository' >&2\nexit 128\n")
    assert codeid.code_identity(a)["code_commit"] is None


@pytest.mark.skipif(shutil.which("git") is None, reason="needs a real git")
def test_code_version_and_commit_from_a_work_tree(tmp_path):
    work = tmp_path / "work"
    work.mkdir()
    a = _copy_pkg(work)
    git = ["git", "-C", str(work), "-c", "user.name=t", "-c", "user.email=t@example.org"]
    subprocess.run([*git, "init", "-q"], check=True)
    subprocess.run([*git, "add", "."], check=True)
    subprocess.run([*git, "commit", "-q", "-m", "x"], check=True)
    want = subprocess.run([*git, "rev-parse", "HEAD"], check=True, capture_output=True, text=True).stdout.strip()
    got = codeid.code_identity(a)
    assert got["code_commit"] == want
    assert got["code_version"] == codeid.code_identity(PKG)["code_version"]


def test_the_default_identity_is_computed_once():
    first = codeid.code_identity()
    assert codeid.code_identity() is first
    assert first["code_id"] == codeid.code_identity(PKG)["code_id"]
    assert set(first) == {"code_id", "code_version", "code_commit"}


def test_a_missing_git_is_not_an_error(tmp_path, monkeypatch):
    a = _copy_pkg(tmp_path / "one")
    monkeypatch.setenv("PATH", str(tmp_path / "nowhere"))
    assert codeid.code_identity(a)["code_commit"] is None


# --- built_by -------------------------------------------------------------------------


def _record(root: Path, stage: str, step: str, **meta) -> None:
    p = root / "runs" / f"stage={stage}" / f"{step}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"stage": stage, "step": step, "outputs": [], **meta}))


def test_built_by_lists_distinct_identities_and_counts_unknown(tmp_path):
    ident = {"code_id": "b" * 64, "code_version": "0.1.0", "code_commit": "abc"}
    _record(tmp_path, "S0", "ingest", **ident)
    _record(tmp_path, "S0", "genecall", **ident)
    _record(tmp_path, "S0", "tapes", code_id="a" * 64, code_version="0.1.0", code_commit=None)
    _record(tmp_path, "S9", "ingest", **ident)  # not in the release
    _record(tmp_path, "rX", "cluster", **ident)
    _record(tmp_path, "rX", "tokens")  # a record from before card C-23
    (tmp_path / "runs" / "stage=S0" / "tapes.publishing.json").write_text("{}")  # an intent file, not a record
    got = json.loads(catalog.built_by(tmp_path, "rX", ["S0"]))
    assert got == {
        "code": [
            {"code_id": "a" * 64, "code_version": "0.1.0", "code_commit": None, "steps": 1},
            {"code_id": "b" * 64, "code_version": "0.1.0", "code_commit": "abc", "steps": 3},
        ],
        "unknown_steps": 1,
    }


def test_built_by_of_a_release_with_no_records_is_unknown_not_an_error(tmp_path):
    assert json.loads(catalog.built_by(tmp_path, "rX", ["S0"])) == {"code": [], "unknown_steps": 0}


def test_catalog_info_carries_built_by(tmp_path):
    _record(tmp_path, "S0", "ingest", code_id="c" * 64, code_version="0.1.0", code_commit=None)
    _record(tmp_path, "rX", "cluster")
    man = {"release_id": "rX", "root": str(tmp_path), "files": [], "artifacts": [], "stages": ["S0"]}
    cat = catalog.generate_catalog(man, tmp_path / "cat.duckdb", tape_version="t")
    con = duckdb.connect(str(cat), read_only=True)
    try:
        info = dict(con.execute("SELECT key, value FROM catalog_info").fetchall())
    finally:
        con.close()
    assert json.loads(info["built_by"]) == json.loads(catalog.built_by(tmp_path, "rX", ["S0"]))
    assert json.loads(info["built_by"])["unknown_steps"] == 1


def test_describe_built_by_lines():
    assert "unknown: built before card C-23" in catalog.describe_built_by(json.dumps({"code": [], "unknown_steps": 0}))
    line = catalog.describe_built_by(json.dumps({"code": [
        {"code_id": "ab" * 32, "code_version": "0.1.0", "code_commit": "740d049", "steps": 2}], "unknown_steps": 3}))
    assert "\n" not in line and "740d049" in line and "abababababab" in line
    assert "3 steps" in line and "unknown" in line


# --- the manifest is not touched -----------------------------------------------------------


def _build_with(src: Path, out: Path) -> tuple[str, str]:
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(src.parent)] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
    tests = Path(__file__).resolve().parent
    r = subprocess.run([sys.executable, "-c", BUILD, str(out), str(tests), str(src)], env=env, capture_output=True,
                       text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    line = next(x for x in r.stdout.splitlines() if x.startswith("RESULT "))
    return tuple(line.split()[1:])


def test_manifest_is_a_function_of_the_data_not_of_the_code(tmp_path):
    one = _build_with(_copy_pkg(tmp_path / "src1"), tmp_path / "r1")
    two = _build_with(_copy_pkg(tmp_path / "other" / "src2"), tmp_path / "r2")
    assert one == two
    assert one == (BASELINE_MANIFEST_HASH, BASELINE_FILES_SHA256)
