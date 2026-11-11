"""Is a path read-only here? r50_probe.py PATH... (50 Gbp build, planning/measurements/2026-11-11_build_r50.md).

Run inside the image with the real root bound read-only (`--bind ROOT:ROOT:ro`). For each path prints whether
its filesystem is mounted read-only (statvfs ST_RDONLY), whether it is writable by access(W_OK), and the
/proc/self/mounts line of the longest mount point above it. Writes nothing.
"""
import os
import sys


def mount_of(path: str) -> str:
    best = ""
    line_best = ""
    with open("/proc/self/mounts") as fh:
        for line in fh:
            mp = line.split()[1]
            if (path == mp or path.startswith(mp.rstrip("/") + "/")) and len(mp) > len(best):
                best, line_best = mp, line.strip()
    return line_best


for p in sys.argv[1:]:
    rp = os.path.realpath(p)
    ro = bool(os.statvfs(rp).f_flag & os.ST_RDONLY)
    print(f"{p} -> {rp}: ST_RDONLY={ro} W_OK={os.access(rp, os.W_OK)} mount: {mount_of(rp)}")
