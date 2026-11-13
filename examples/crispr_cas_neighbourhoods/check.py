"""Compare a run of the CRISPR-Cas example with what it must show (standard library only).

    python -I check.py --expected expected.tsv [--recorded recorded/<release>.tsv] [--write-recorded FILE]
                       OBSERVED.tsv [OBSERVED.tsv ...]

OBSERVED files are TSV with a column line `check_id<TAB>value` (as `bio-box sql --out` and `oracle.py` write
them); a check id may appear in one file only. `expected.tsv` has one row per check: `check_id`, `asserts` (what
it asserts, in words), `kind` and `expectation`. Two kinds:

- `invariant`: must hold on every release that holds the control genomes. The expectation is `== V` (equal as
  text, or as numbers when both parse), `>= N`, `<= N` or `in A,B,...`. A failed or missing invariant fails.
- `recorded`: a value of one release only (a count). Its expectation is `recorded`; the value is compared with
  the release's file `recorded/<release>.tsv` (`check_id<TAB>value`) when there is one. A difference is reported,
  not failed: a new release is expected to differ, and adds its own file (`--write-recorded`) instead of editing
  an old one.

Prints one line per check, `<status>  <check_id>  observed=<v>  expected=<e>`, status one of ok, FAIL (an
invariant), differs (a recorded value), new (no recorded value for this release), missing (expected but not
observed) or unexpected (observed, in no row of expected.tsv); then a summary line. Exits 1 when an invariant
fails or is missing, 2 on a malformed input, else 0.
"""

from __future__ import annotations

import argparse
import csv
import sys

KINDS = ("invariant", "recorded")


class Malformed(Exception):
    pass


def read_tsv(path: str, columns: tuple[str, ...]) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))
    rows = [r for r in rows if r and any(c.strip() for c in r)]
    if not rows or tuple(c.lstrip("#") for c in rows[0][:len(columns)]) != columns:
        raise Malformed(f"{path}: the first line must be the columns {', '.join(columns)}")
    out = []
    for k, r in enumerate(rows[1:], 2):
        if len(r) < len(columns):
            raise Malformed(f"{path}: line {k} has {len(r)} fields, not {len(columns)}")
        out.append(dict(zip(columns, r)))
    return out


def values(paths: list[str]) -> dict[str, str]:
    """{check_id: value} of the observed files; a check id twice is malformed."""
    got: dict[str, str] = {}
    for p in paths:
        for r in read_tsv(p, ("check_id", "value")):
            if r["check_id"] in got:
                raise Malformed(f"{p}: check {r['check_id']} observed twice")
            got[r["check_id"]] = r["value"]
    return got


def _num(s: str) -> float | None:
    try:
        return float(s)
    except ValueError:
        return None


def holds(value: str, expectation: str) -> bool:
    """Whether an observed value meets an invariant's expectation (`== V`, `>= N`, `<= N`, `in A,B`)."""
    op, _, rhs = expectation.strip().partition(" ")
    rhs = rhs.strip()
    if op == "==":
        a, b = _num(value), _num(rhs)
        return a == b if a is not None and b is not None else value == rhs
    if op in (">=", "<="):
        a, b = _num(value), _num(rhs)
        if b is None:
            raise Malformed(f"expectation {expectation!r}: {rhs!r} is not a number")
        if a is None:
            return False
        return a >= b if op == ">=" else a <= b
    if op == "in":
        return value in [x.strip() for x in rhs.split(",")]
    raise Malformed(f"expectation {expectation!r}: the operator must be ==, >=, <= or in")


def same(a: str, b: str) -> bool:
    x, y = _num(a), _num(b)
    return x == y if x is not None and y is not None else a == b


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--expected", required=True)
    ap.add_argument("--recorded", help="recorded/<release>.tsv; absent or missing file: recorded checks are new")
    ap.add_argument("--write-recorded", help="write this run's recorded values to this new file")
    ap.add_argument("observed", nargs="+")
    a = ap.parse_args(argv)
    try:
        expected = read_tsv(a.expected, ("check_id", "asserts", "kind", "expectation"))
        seen: set[str] = set()
        for r in expected:
            if r["kind"] not in KINDS:
                raise Malformed(f"{a.expected}: check {r['check_id']} has kind {r['kind']!r}, not one of {KINDS}")
            if r["check_id"] in seen:
                raise Malformed(f"{a.expected}: check {r['check_id']} listed twice")
            if r["kind"] == "recorded" and r["expectation"] != "recorded":
                raise Malformed(f"{a.expected}: recorded check {r['check_id']} must have the expectation 'recorded'")
            if r["kind"] == "invariant":
                holds("0", r["expectation"])  # malformed expectations fail here, before any line is printed
            seen.add(r["check_id"])
        got = values(a.observed)
        rec: dict[str, str] = {}
        if a.recorded:
            try:
                rec = {r["check_id"]: r["value"] for r in read_tsv(a.recorded, ("check_id", "value"))}
            except FileNotFoundError:
                rec = {}
    except (Malformed, OSError) as e:
        print(f"check: {e}", file=sys.stderr)
        return 2
    counts = {"ok": 0, "FAIL": 0, "differs": 0, "new": 0, "missing": 0, "unexpected": 0}
    failed = 0
    lines = []
    for r in expected:
        cid, kind, exp = r["check_id"], r["kind"], r["expectation"]
        if cid not in got:
            status = "missing"
            failed += kind == "invariant"
            shown = "-"
        else:
            shown = got[cid]
            if kind == "invariant":
                status = "ok" if holds(shown, exp) else "FAIL"
                failed += status == "FAIL"
            else:
                exp = rec.get(cid, "-")
                status = "new" if cid not in rec else ("ok" if same(shown, rec[cid]) else "differs")
        counts[status] += 1
        lines.append(f"{status}\t{cid}\tobserved={shown}\texpected={exp}")
    for cid in sorted(got.keys() - seen):
        counts["unexpected"] += 1
        lines.append(f"unexpected\t{cid}\tobserved={got[cid]}\texpected=-")
    for line in lines:
        print(line)
    print("check: " + ", ".join(f"{n} {k}" for k, n in counts.items()) + (f"; {failed} invariants failed" if failed
                                                                            else "; invariants hold"))
    if a.write_recorded:
        kinds = {r["check_id"]: r["kind"] for r in expected}
        try:
            with open(a.write_recorded, "x", newline="") as fh:
                fh.write("check_id\tvalue\n")
                for cid in sorted(got):
                    if kinds.get(cid) == "recorded":
                        fh.write(f"{cid}\t{got[cid]}\n")
        except OSError as e:
            print(f"check: --write-recorded: {e}", file=sys.stderr)
            return 2
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
