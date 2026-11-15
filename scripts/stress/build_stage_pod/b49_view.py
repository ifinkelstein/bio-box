"""Card B-49: a scratch root of file links over the data root, read-only on it, for a release that reuses a clustering.

    python3 -I b49_view.py SRC DST S1,P1,P2,P3,P4,P5,P6 rM.1,rM.2

DST mirrors SRC with real directories and one symlink per file, for the files of the stages named and of the releases
named (their `release=<R>` partitions and `runs/stage=<R>`), and the unpartitioned files (`l0seq/store_params.json`,
...). Files of other stages and releases, and every top-level entry whose name starts with a dot (`.incoming`) or is
`STOP`, are left out. Python's rglob does not descend into a linked directory (card C-34), so links are made per file:
the manifest and `release verify` of the new release then see the stages' files. The new release's files, records,
ledger and manifest land in new real directories of DST (`release=<new>`, `runs/stage=<new>`), so nothing is written
through a link; the steps also run with SRC bound read-only. Prints the files linked and their bytes. Refuses when DST
exists.
"""
import os
import sys

src, dst = sys.argv[1], sys.argv[2]
stages, releases = set(sys.argv[3].split(",")), set(sys.argv[4].split(","))
os.makedirs(dst)  # refuses an existing DST
n = size = 0


def keep(rel_parts: list[str]) -> bool:
    for p in rel_parts:
        if p.startswith("stage=") and p.removeprefix("stage=") not in stages | releases:
            return False
        if p.startswith("release=") and p.removeprefix("release=") not in releases:
            return False
    return True


for top in sorted(os.listdir(src)):
    if top.startswith(".") or top == "STOP":
        continue
    for d, subdirs, files in os.walk(os.path.join(src, top)):
        rel = os.path.relpath(d, src)
        parts = rel.split(os.sep)
        if not keep(parts):
            subdirs[:] = []
            continue
        subdirs.sort()
        os.makedirs(os.path.join(dst, rel), exist_ok=True)
        for f in sorted(files):
            if f.startswith(".") or not keep(parts + [f]):
                continue
            os.symlink(os.path.join(d, f), os.path.join(dst, rel, f))
            n += 1
            size += os.path.getsize(os.path.join(d, f))
print(f"{n} files linked into {dst} ({size:,} bytes in SRC) for stages {sorted(stages)} and releases {sorted(releases)}")
