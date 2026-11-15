"""Card C-34: the before-and-after check of the data root and the list of the files the step added.

    python -I c34_check.py --root ROOT --before DIR --after DIR --rehearsal REH --out DIR [--stage S ...]

Reads the two snapshots the driver wrote (`listing.tsv`: path, bytes, mtime of every file of ROOT; `stages.sha256`:
sha256 of every file under a stage=<S> directory of the stages named) and checks that every file of the first is in
the second with the same bytes, mtime and (where hashed) sha256. The files only in the second must be exactly, for
each stage, its label part files, its `labels` ledger and its `labels` run record; each is hashed, and the label and
ledger files must have the sha256 of the same files in the rehearsal root REH (the same step on the same inputs).
Writes OUT/root_added.tsv (path, bytes, sha256, sorted by path) and prints one line per finding, then `check: ok` or
`check: FAILED (<n> problems)` last; exit 0 only when ok. Reads only; writes only OUT.
"""
import argparse
import hashlib
import os
import sys

ap = argparse.ArgumentParser()
ap.add_argument("--root", required=True)
ap.add_argument("--before", required=True)
ap.add_argument("--after", required=True)
ap.add_argument("--rehearsal", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--stage", action="append", default=[])
a = ap.parse_args()


def listing(d: str) -> dict[str, tuple[str, str]]:
    out = {}
    with open(os.path.join(d, "listing.tsv")) as fh:
        for line in fh:
            path, size, mtime = line.rstrip("\n").split("\t")
            out[path] = (size, mtime)
    return out


def shas(d: str) -> dict[str, str]:
    out = {}
    with open(os.path.join(d, "stages.sha256")) as fh:
        for line in fh:
            sha, path = line.rstrip("\n").split("  ", 1)
            out[path.removeprefix("./")] = sha
    return out


def sha256(p: str) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


problems = []
b, af = listing(a.before), listing(a.after)
for path, v in sorted(b.items()):
    if path not in af:
        problems.append(f"gone: {path}")
    elif af[path] != v:
        problems.append(f"changed (bytes, mtime): {path}: {v} -> {af[path]}")
sb, sa = shas(a.before), shas(a.after)
for path, sha in sorted(sb.items()):
    if sa.get(path) != sha:
        problems.append(f"sha256 changed: {path}")
print(f"before: {len(b)} files; after: {len(af)}; {len(sb)} stage files hashed before, {len(sa)} after")
new = sorted(set(af) - set(b))
want_dirs = {f"l0seq/contig_label/stage={s}/" for s in a.stage} | {f"ledger/stage={s}/step=labels/" for s in a.stage}
records = {f"runs/stage={s}/labels.json" for s in a.stage}
for path in new:
    if path in records:
        continue
    if not any(path.startswith(d) and path.endswith(".parquet") and "/" not in path[len(d):] for d in want_dirs):
        problems.append(f"unexpected new file: {path}")
for r in sorted(records):
    if r not in new:
        problems.append(f"no run record: {r}")
for d in sorted(want_dirs):
    if not any(p.startswith(d) for p in new):
        problems.append(f"no new file under {d}")
os.makedirs(a.out, exist_ok=True)
rows = []
for path in new:
    sha = sha256(os.path.join(a.root, path))
    rows.append((path, af[path][0], sha))
    if path.endswith(".parquet"):
        twin = os.path.join(a.rehearsal, path)
        if not os.path.isfile(twin):
            problems.append(f"not in the rehearsal: {path}")
        elif sha256(twin) != sha:
            problems.append(f"differs from the rehearsal: {path}")
        else:
            print(f"same sha256 as the rehearsal: {path} {sha}")
with open(os.path.join(a.out, "root_added.tsv"), "w") as fh:
    fh.write("path\tbytes\tsha256\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")
print(f"added: {len(rows)} files, {sum(int(r[1]) for r in rows)} bytes; list in {os.path.join(a.out, 'root_added.tsv')}")
for p in problems:
    print(p)
print("check: ok" if not problems else f"check: FAILED ({len(problems)} problems)")
sys.exit(0 if not problems else 1)
