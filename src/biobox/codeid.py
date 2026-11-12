"""The identity of the code that runs (card C-23).

`code_identity()` gives `code_id` (sha256 over the `biobox` package's files as imported: the relative path and the
bytes of every file under the package directory except bytecode caches, in sorted path order), `code_version` (the
package version from its metadata, else the package's `__version__`) and `code_commit` (the git commit when the
package directory is inside a git work tree and `git` answers; None otherwise, never an error). `code_id` depends
on the source bytes only: the same on any host and at any install path. Run records carry all three
(`runstep.run_step`); release data never does, so the manifest stays a function of the data.
"""

from __future__ import annotations

import hashlib
import subprocess
from functools import lru_cache
from importlib import metadata
from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent
_SKIP_DIRS = {"__pycache__"}
_SKIP_SUFFIXES = {".pyc", ".pyo"}


def _files(pkg: Path) -> list[Path]:
    out = [p for p in pkg.rglob("*")
           if p.is_file() and p.suffix not in _SKIP_SUFFIXES
           and not _SKIP_DIRS.intersection(p.relative_to(pkg).parts)]
    return sorted(out, key=lambda p: p.relative_to(pkg).as_posix())


def tree_hash(pkg: Path) -> str:
    """sha256 over the relative path and the bytes of every file under `pkg` (sorted by path)."""
    h = hashlib.sha256()
    for p in _files(pkg):
        rel = p.relative_to(pkg).as_posix().encode()
        data = p.read_bytes()
        h.update(len(rel).to_bytes(4, "big") + rel + len(data).to_bytes(8, "big") + data)
    return h.hexdigest()


def _version() -> str | None:
    try:
        return metadata.version("biobox")
    except metadata.PackageNotFoundError:
        pass
    try:
        from biobox import __version__
        return __version__
    except Exception:
        return None


def _commit(pkg: Path) -> str | None:
    """The git HEAD commit when `pkg` is inside a work tree and git answers; else None."""
    try:
        r = subprocess.run(["git", "-C", str(pkg), "rev-parse", "--is-inside-work-tree", "HEAD"],
                           capture_output=True, text=True, timeout=5)
    except (OSError, subprocess.SubprocessError):
        return None
    lines = r.stdout.split()
    if r.returncode != 0 or len(lines) != 2 or lines[0] != "true":
        return None
    return lines[1]


def code_identity(pkg: Path | str | None = None) -> dict:
    """`{"code_id", "code_version", "code_commit"}` of the package at `pkg` (default: this one, computed once)."""
    if pkg is None:
        return _own()
    return _identity(Path(pkg))


def _identity(pkg: Path) -> dict:
    return {"code_id": tree_hash(pkg), "code_version": _version(), "code_commit": _commit(pkg)}


@lru_cache(maxsize=1)
def _own() -> dict:
    return _identity(PACKAGE_DIR)
