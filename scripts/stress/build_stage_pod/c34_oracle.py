"""Card C-34: the oracle on the pod. The label files `stage labels` wrote for a stage, against those a re-ingest of the
stage's inputs with the current code wrote into a scratch root.

    python -I c34_oracle.py --stage S1 --step REH --ingest ORACLE --root ROOT

Prints, for the stage, the sha256 of each label part of REH (the step's) and of ORACLE (the re-ingest's): they must be
equal. Then the `contig_occurrence` of ORACLE against ROOT's (that the re-ingest read the same inputs as the stage's
ingest did): their rows must be equal; their sha256 are printed, and may differ only if the Parquet writer differs
(ROOT's stage was written in an older image). On a difference, the first differing rows. Last line `oracle: same` or
`oracle: DIFFERENT`; exit 0 only for `same`. Reads only.
"""
import argparse
import hashlib
import sys
from pathlib import Path

import pyarrow.parquet as pq

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True)
ap.add_argument("--step", required=True, type=Path)
ap.add_argument("--ingest", required=True, type=Path)
ap.add_argument("--root", required=True, type=Path)
a = ap.parse_args()


def files(base: Path, rel: str) -> dict[str, str]:
    d = base / rel
    return {p.relative_to(d).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(d.rglob("*.parquet"))} if d.is_dir() else {}


def rows(base: Path, rel: str) -> list[dict]:
    return [r for f in sorted((base / rel).rglob("*.parquet")) for r in pq.read_table(f, partitioning=None).to_pylist()]


def first_differences(x: list[dict], y: list[dict]) -> None:
    print(f"rows: {len(x)} and {len(y)}")
    shown = 0
    for i, (r, s) in enumerate(zip(x, y)):
        if r != s:
            print(f"row {i}: {r}\n    != {s}")
            shown += 1
            if shown == 5:
                break


ok = True
for what, rel, left, right, need in (
        ("contig_label", f"l0seq/contig_label/stage={a.stage}", a.step, a.ingest, "sha256"),
        ("contig_occurrence", f"l0seq/contig_occurrence/stage={a.stage}", a.ingest, a.root, "rows")):
    x, y = files(left, rel), files(right, rel)
    for name in sorted(set(x) | set(y)):
        print(f"{what} {name}: {left} {x.get(name)}; {right} {y.get(name)}")
    same_sha = bool(x) and x == y
    rx, ry = rows(left, rel), rows(right, rel)
    same_rows = bool(rx) and rx == ry
    print(f"{what} of stage {a.stage}: {len(x)} and {len(y)} files, {len(rx)} and {len(ry)} rows; sha256 "
          f"{'same' if same_sha else 'different'}; rows {'same' if same_rows else 'different'}")
    if not (same_sha if need == "sha256" else same_rows):
        ok = False
        first_differences(rx, ry)
print("oracle: same" if ok else "oracle: DIFFERENT")
sys.exit(0 if ok else 1)
