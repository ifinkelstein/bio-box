"""`python -m biobox.doctor`: record the host and image the pipeline actually runs on.

Prints one JSON object: CPU model and SIMD flags, cores, memory, the cgroup
memory limit, ulimits, filesystems of the working paths, and versions of the
native tools and Python packages. The S0 gate stores this next to `test_run`.
"""

from __future__ import annotations

import json
import os
import platform
import resource
import shutil
import subprocess
import sys
from importlib import metadata
from pathlib import Path

TOOLS = {
    "mmseqs": ["mmseqs", "version"],
    "hmmsearch": ["hmmsearch", "-h"],
    "cmsearch": ["cmsearch", "-h"],
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
    cg = Path("/proc/self/cgroup")
    if cg.exists():
        rel = cg.read_text().strip().splitlines()[-1].split(":", 2)[-1]
        limit = Path("/sys/fs/cgroup") / rel.lstrip("/") / "memory.max"
        if limit.exists():
            raw = limit.read_text().strip()
            out["cgroup_memory_max"] = raw if raw == "max" else f"{int(raw) / 2**30:.0f} GiB"
    return out


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


def report(paths: list[str] | None = None) -> dict:
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "cpu": _cpu(),
        "memory": _memory(),
        "ulimits": _limits(),
        "filesystems": _filesystems(paths or ["/tmp", os.getcwd()]),
        "tools": {k: _tool(v) for k, v in TOOLS.items()},
        "packages": _packages(),
        "image": os.environ.get("BIOBOX_IMAGE_DIGEST"),
    }


def main(argv: list[str] | None = None) -> None:
    print(json.dumps(report(argv if argv else None), indent=2, default=str))


if __name__ == "__main__":
    main(sys.argv[1:])
