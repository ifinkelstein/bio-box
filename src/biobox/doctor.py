"""`python -m biobox.doctor`: record the host and image the pipeline actually runs on.

Prints one JSON object: CPU model and SIMD flags, cores, memory, the cgroup
memory limit, ulimits, filesystems of the working paths, and versions of the
native tools and Python packages. The S0 gate stores this next to `test_run`.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import resource
import shutil
import subprocess
import sys
import tempfile
from importlib import metadata
from pathlib import Path

from biobox.profile import load as load_profile

TOOLS = {
    "mmseqs": ["mmseqs", "version"],
    "hmmsearch": ["hmmsearch", "-h"],
    "cmsearch": ["cmsearch", "-h"],
    "minced": ["minced", "--version"],  # card E-12; prints `minced <version>` first
    "singularity": ["singularity", "--version"],
}
PACKAGES = ["bio-box", "pyarrow", "duckdb", "pyrodigal", "pyrodigal-gv", "blake3", "pyyaml"]


def _cpu() -> dict:
    info = {"machine": platform.machine(), "logical_cpus": os.cpu_count()}
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        text = cpuinfo.read_text()
        for line in text.splitlines():
            if line.startswith("model name"):
                info["model"] = line.split(":", 1)[1].strip()
                break
        flags = next((ln.split(":", 1)[1].split() for ln in text.splitlines() if ln.startswith("flags")), [])
        info["simd"] = sorted(f for f in flags if f in {"sse4_1", "sse4_2", "avx", "avx2", "avx512f", "avx512bw"})
    return info


def _memory() -> dict:
    out: dict = {}
    meminfo = Path("/proc/meminfo")
    if meminfo.exists():
        for line in meminfo.read_text().splitlines():
            if line.startswith("MemTotal"):
                out["total_gib"] = round(int(line.split()[1]) / 2**20, 1)
    cg = cgroup_memory()
    if cg["available"]:
        out["cgroup_memory_max"] = cg["memory_max"]
    return out


def _read(path: Path) -> str | None:
    try:
        return path.read_text().strip()
    except OSError:
        return None


def _limit(raw: str | None) -> int | None:
    """A cgroup limit file's value: bytes, or None for `max`, a missing file or unparsable text."""
    if raw is None or not raw.isdigit():
        return None
    return int(raw)


def _gib(n: int) -> str:
    return f"{n / 2**30:.1f} GiB".replace(".0 GiB", " GiB")


def cgroup_memory(proc_root: str | Path = "/proc", cgroup_root: str | Path = "/sys/fs/cgroup") -> dict:
    """Memory limit of the cgroup this process runs in (cgroup v2 only).

    Walks from the process's cgroup up to the root and reports the smallest numeric `memory.max`
    (`max` means unlimited at that level), plus `memory.high` and `memory.current` of the
    process's own cgroup. Never raises: cgroup v1 or a missing tree gives `available: False`
    with a reason.
    """
    proc_root, cgroup_root = Path(proc_root), Path(cgroup_root)
    text = _read(proc_root / "self" / "cgroup")
    if text is None:
        return {"available": False, "reason": f"{proc_root}/self/cgroup is not readable"}
    rel = None
    for line in text.splitlines():
        if line.startswith("0::"):
            rel = line[3:].strip()
            break
    if rel is None:
        return {"available": False, "reason": "cgroup v1 (no unified '0::' hierarchy in /proc/self/cgroup)"}
    if not cgroup_root.is_dir():
        return {"available": False, "reason": f"{cgroup_root} does not exist"}
    parts = [p for p in rel.split("/") if p]
    best: int | None = None
    best_at: str | None = None
    found = False
    for depth in range(len(parts), -1, -1):
        raw = _read(cgroup_root.joinpath(*parts[:depth]) / "memory.max")
        if raw is None:
            continue
        found = True
        value = _limit(raw)
        if value is not None and (best is None or value < best):
            best, best_at = value, "/" + "/".join(parts[:depth])
    if not found:
        return {"available": False, "reason": f"no memory.max under {cgroup_root}/{'/'.join(parts)} or its parents"}
    own = cgroup_root.joinpath(*parts)
    return {
        "available": True,
        "version": 2,
        "path": "/" + "/".join(parts),
        "memory_max_bytes": best,
        "memory_max": "max (unlimited)" if best is None else _gib(best),
        "limited_at": best_at,
        "memory_high_bytes": _limit(_read(own / "memory.high")),
        "memory_current_bytes": _limit(_read(own / "memory.current")),
    }


def _limits() -> dict:
    names = {"nofile": resource.RLIMIT_NOFILE, "as": resource.RLIMIT_AS}
    return {k: resource.getrlimit(v) for k, v in names.items()}


def _tool(cmd: list[str]) -> str | None:
    if shutil.which(cmd[0]) is None:
        return None
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=30, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return f"error: {exc}"
    text = (res.stdout or res.stderr).strip().splitlines()
    for line in text:  # HMMER/Infernal print the version on a '# ... 3.4 (...)' line
        if line.startswith("# HMMER") or line.startswith("# INFERNAL"):
            return line.lstrip("# ").strip()
    return text[0].strip() if text else ""


def _packages() -> dict:
    out = {}
    for p in PACKAGES:
        try:
            out[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            out[p] = None
    return out


def _filesystems(paths: list[str]) -> dict:
    out = {}
    for p in paths:
        if Path(p).exists():
            usage = shutil.disk_usage(p)
            out[p] = {"free_gib": round(usage.free / 2**30), "total_gib": round(usage.total / 2**30)}
    return out


def _mount_table(mounts_file: str | Path) -> list[tuple[str, str]]:
    """(mount point, fs type) pairs from a /proc/mounts-style file; [] if unreadable."""
    rows = []
    for line in (_read(Path(mounts_file)) or "").splitlines():
        f = line.split()
        if len(f) >= 3:
            rows.append((re.sub(r"\\([0-7]{3})", lambda m: chr(int(m.group(1), 8)), f[1]), f[2]))
    return rows


def fs_type(path: str, mounts_file: str | Path = "/proc/mounts") -> str | None:
    """Filesystem type of the longest mount point containing `path`, or None if unknown."""
    target = os.path.abspath(path)
    best: tuple[int, str] | None = None
    for mount, kind in _mount_table(mounts_file):
        inside = mount == "/" or target == mount or target.startswith(mount.rstrip("/") + "/")
        if inside and (best is None or len(mount) >= best[0]):
            best = (len(mount), kind)
    return best[1] if best else None


def _can_write(directory: str) -> bool:
    """True if a uniquely named temporary file can be created inside `directory`; it is removed."""
    try:
        fd, name = tempfile.mkstemp(prefix=".biobox-doctor-", dir=directory)
    except OSError:
        return False
    os.close(fd)
    try:
        os.unlink(name)
    except OSError:
        pass
    return True


def probe_path(path: str, mounts_file: str | Path = "/proc/mounts") -> dict:
    """Existence, kind, access, filesystem type and free space of one path. Never creates it."""
    out: dict = {"exists": os.path.exists(path), "is_dir": False, "readable": False, "writable": False}
    if not out["exists"]:
        out.update(status="missing", fs_type=None, free_gib=None)
        return out
    out["is_dir"] = os.path.isdir(path)
    out["readable"] = os.access(path, os.R_OK)
    out["writable"] = _can_write(path) if out["is_dir"] else os.access(path, os.W_OK)
    out["status"] = "ok"
    out["fs_type"] = fs_type(path, mounts_file)
    try:
        out["free_gib"] = round(shutil.disk_usage(path).free / 2**30, 1)
    except OSError:
        out["free_gib"] = None
    return out


def profile_paths(profile: str | Path) -> list[str]:
    """The `tmp`, `durable` and `bulk` paths of a host profile."""
    p = load_profile(profile)
    return [p.tmp, p.durable, p.bulk]


def report(paths: list[str] | None = None, profile: str | Path | None = None) -> dict:
    """`paths` and the profile's paths are probed under `paths`; `filesystems` keeps its old meaning."""
    probe = list(paths or [])
    if profile is not None:
        probe += profile_paths(profile)
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "cpu": _cpu(),
        "memory": _memory(),
        "cgroup": cgroup_memory(),
        "ulimits": _limits(),
        "filesystems": _filesystems(paths or ["/tmp", os.getcwd()]),
        "paths": {p: probe_path(p) for p in probe},
        "tools": {k: _tool(v) for k, v in TOOLS.items()},
        "packages": _packages(),
        "image": os.environ.get("BIOBOX_IMAGE_DIGEST"),
    }


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="biobox.doctor", description=__doc__.splitlines()[0])
    ap.add_argument("--paths", nargs="*", default=[], metavar="P", help="paths to probe (exists, writable, fs, free)")
    ap.add_argument("--profile", metavar="FILE", help="host profile whose tmp, durable and bulk paths are probed")
    args = ap.parse_args(argv)
    print(json.dumps(report(args.paths or None, args.profile), indent=2, default=str))


if __name__ == "__main__":
    main(sys.argv[1:])
