"""Stage gate as a query over `test_run` (data_layer_and_tests.md §4.12).

    python -m biobox.gate S0 --dir $ATD_ROOT/tests [--image DIGEST] [--release rS0.1]

A stage passes when every blocking test has at least one row and all its rows pass
(for the given runtime image and release; R0-R2 rows count without a release), and
every required ring is green in its latest session. Prints a JSON report.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import duckdb
import yaml

from biobox.catalog import _sql_str

REPO = Path(__file__).resolve().parents[2]
GATES = REPO / "tests" / "registry" / "gates.yaml"
_RANGE = re.compile(r"^([A-Z]+-?)(\d+)\.\.\1(\d+)$")

GATE_SQL = """
SELECT g.test_id,
       count(t.test_id) AS n_rows,
       coalesce(bool_and(t.outcome = 'pass'), false) AS all_pass
FROM gate_spec g
LEFT JOIN test_run t
  ON t.test_id = g.test_id
 AND t.outcome <> 'skip'  -- a skipped test was not run: it neither passes nor fails the gate
 AND ($image::VARCHAR IS NULL OR t.runtime_image_digest = $image::VARCHAR)
 AND ($release::VARCHAR IS NULL OR t.release_id = $release::VARCHAR OR t.ring IN ('R0', 'R1', 'R2'))
WHERE g.stage = $stage AND g.gate = $gate
GROUP BY g.test_id
ORDER BY g.test_id
"""

RING_SQL = """
WITH latest AS (
  SELECT ring, arg_max(filename, started_at) AS f FROM test_run GROUP BY ring
)
SELECT t.ring, count(*) AS n,
       count(*) FILTER (WHERE outcome = 'pass') AS passed,
       count(*) FILTER (WHERE outcome = 'skip') AS skipped,
       count(*) FILTER (WHERE outcome IN ('fail', 'error')) AS failed
FROM test_run t JOIN latest l ON t.ring = l.ring AND t.filename = l.f
GROUP BY t.ring ORDER BY t.ring
"""


def expand(item: str) -> list[str]:
    m = _RANGE.match(item)
    if not m:
        return [item]
    prefix, lo, hi = m.group(1), m.group(2), m.group(3)
    return [f"{prefix}{i:0{len(lo)}d}" for i in range(int(lo), int(hi) + 1)]


def load_gates(path: Path = GATES) -> dict:
    raw = yaml.safe_load(Path(path).read_text())["stages"]
    out = {}
    for stage, g in raw.items():
        out[stage] = {
            "B": [t for item in g.get("B", []) for t in expand(item)],
            "W": [t for item in g.get("W", []) for t in expand(item)],
            "rings": g.get("rings", []),
            "milestones": g.get("milestones", []),
        }
    return out


def _files(test_dir: Path) -> list[str]:
    idx = Path(test_dir) / "test_run_index.json"
    if idx.exists():
        return [str(Path(test_dir) / f) for f in json.loads(idx.read_text())["files"]]
    return sorted(str(p) for p in (Path(test_dir) / "test_run").glob("*.parquet"))


def gate_report(test_dir: Path, stage: str, image: str | None = None, release: str | None = None) -> dict:
    gates = load_gates()[stage]
    con = duckdb.connect()
    files = _files(test_dir)
    if not files:
        raise FileNotFoundError(f"no test_run files under {test_dir}")
    listed = "[" + ", ".join(_sql_str(f) for f in files) + "]"
    con.execute(f"CREATE VIEW test_run AS SELECT * FROM read_parquet({listed}, filename = true)")
    con.execute("CREATE TABLE gate_spec (stage VARCHAR, test_id VARCHAR, gate VARCHAR)")
    con.executemany("INSERT INTO gate_spec VALUES (?, ?, ?)",
                    [(stage, t, g) for g in ("B", "W") for t in gates[g]])
    report = {"stage": stage, "image": image, "release": release, "files": len(files)}
    for g in ("B", "W"):
        rows = con.execute(GATE_SQL, {"image": image, "release": release, "stage": stage, "gate": g}).fetchall()
        report[g] = {
            "passing": [t for t, n, ok in rows if n and ok],
            "failing": [t for t, n, ok in rows if n and not ok],
            "missing": [t for t, n, ok in rows if not n],
        }
    rings = {r: {"n": n, "passed": p, "skipped": s, "failed": f}
             for r, n, p, s, f in con.execute(RING_SQL).fetchall()}
    report["rings"] = {r: {**rings.get(r, {"n": 0}),
                           "green": r in rings and rings[r]["failed"] == 0 and rings[r]["passed"] > 0}
                       for r in gates["rings"]}
    report["milestones"] = gates["milestones"]
    report["passes"] = (not report["B"]["failing"] and not report["B"]["missing"]
                        and all(v["green"] for v in report["rings"].values()))
    return report


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage")
    ap.add_argument("--dir", type=Path, required=True, help="the tests/ directory holding test_run/")
    ap.add_argument("--image")
    ap.add_argument("--release")
    a = ap.parse_args(argv)
    rep = gate_report(a.dir, a.stage, a.image, a.release)
    print(json.dumps(rep, indent=1))
    return 0 if rep["passes"] else 1


if __name__ == "__main__":
    sys.exit(main())
