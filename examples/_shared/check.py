"""Compare a run of a worked example with what it must show (standard library only; shared by every example under
`examples/`, card E-11; written for the CRISPR-Cas example).

    python -I check.py --expected expected.tsv [--recorded recorded/<release>.tsv] [--write-recorded FILE]
                       [--table FILE] OBSERVED.tsv [OBSERVED.tsv ...]

OBSERVED files are TSV with a column line `check_id<TAB>value` (as `bio-box sql --out` and `oracle.py` write
them); a check id may appear in one file only. `expected.tsv` has one row per check: `check_id`, `asserts` (what
it asserts, in words), `kind` and `expectation`. Two kinds:

- `invariant`: must hold on every release that holds the control genomes. The expectation is `== V` (equal as
  text, or as numbers when both parse), `>= N`, `<= N` or `in A,B,...`. A failed or missing invariant fails.
- `recorded`: a value of one release only (a count). Its expectation is `recorded`; the value is compared with
  the release's file `recorded/<release>.tsv` (`check_id<TAB>value`) when there is one. A difference is reported,
  not failed: a new release is expected to differ, and adds its own file (`--write-recorded`) instead of editing
  an old one.

A check that needs a table not every release holds (card E-14: the CRISPR array checks need `l1_crispr`) starts its
expectation with `[requires <table>] `, for example `[requires l1_crispr] == true` or `[requires l1_crispr]
recorded`. On a release without that table its SQL observes the value `not applicable`; such a check, or one not
observed at all (a run made before the check existed), is reported `n/a`: neither ok nor failed, left out of the
`--table` file (which lists the checks that apply; `scripts/release_examples.sh` counts every invariant row there
that is not ok as failed), and never written by `--write-recorded`, so a release without the table gets no recorded
value for it. Any other value is checked as usual. Only a check that says `[requires ...]` can be `n/a`: an
invariant without it that observes `not applicable` fails, as before.

Prints one line per check, `<status>  <check_id>  observed=<v>  expected=<e>`, status one of ok, FAIL (an
invariant), differs (a recorded value), new (no recorded value for this release), missing (expected but not
observed), n/a (above) or unexpected (observed, in no row of expected.tsv); then a summary line (it names n/a only
when some check is). `--table FILE` also writes the
same lines as a new TSV file `check, kind, expected, observed, verdict` (`-` for none), which
`scripts/release_examples.sh` gathers into one summary. Exits 1 when an invariant fails or is missing, 2 on a
malformed input, else 0.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys

KINDS = ("invariant", "recorded")
NOT_APPLICABLE = "not applicable"  # the value a check's SQL observes on a release without the table it requires
_REQUIRES = re.compile(r"^\[requires ([A-Za-z0-9_]+)\]\s+(\S.*)$")


def requires(expectation: str) -> tuple[str | None, str]:
    """(the table a check requires, or None; the expectation without the `[requires <table>] ` prefix)."""
    m = _REQUIRES.match(expectation.strip())
    if m is None:
        if expectation.strip().startswith("[requires"):
            raise Malformed(f"expectation {expectation!r}: write `[requires <table>] <expectation>`")
        return None, expectation
    return m[1], m[2]


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
    ap.add_argument("--table", help="write every check's line as TSV (check, kind, expected, observed, verdict) to "
                                    "this new file")
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
            _, bare = requires(r["expectation"])
            if r["kind"] == "recorded" and bare != "recorded":
                raise Malformed(f"{a.expected}: recorded check {r['check_id']} must have the expectation 'recorded'")
            if r["kind"] == "invariant":
                holds("0", bare)  # malformed expectations fail here, before any line is printed
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
    counts = {"ok": 0, "n/a": 0, "FAIL": 0, "differs": 0, "new": 0, "missing": 0, "unexpected": 0}
    failed = 0
    lines = []
    table = []  # (check, kind, expected, observed, verdict)
    for r in expected:
        cid, kind = r["check_id"], r["kind"]
        need, exp = requires(r["expectation"])
        if need is not None and got.get(cid, NOT_APPLICABLE) == NOT_APPLICABLE:
            # card E-14: the release lacks the table the check requires (or the run predates the check)
            status, shown, exp = "n/a", got.get(cid, "-"), r["expectation"]
        elif cid not in got:
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
        if status != "n/a":  # the table holds the checks that apply (scripts/release_examples.sh counts it)
            table.append((cid, kind, exp, shown, status))
    for cid in sorted(got.keys() - seen):
        counts["unexpected"] += 1
        lines.append(f"unexpected\t{cid}\tobserved={got[cid]}\texpected=-")
        table.append((cid, "-", "-", got[cid], "unexpected"))
    if a.table:
        try:
            with open(a.table, "x", newline="") as fh:
                fh.write("check\tkind\texpected\tobserved\tverdict\n")
                fh.writelines("\t".join(row) + "\n" for row in table)
        except OSError as e:
            print(f"check: --table: {e}", file=sys.stderr)
            return 2
    for line in lines:
        print(line)
    shown_counts = {k: n for k, n in counts.items() if k != "n/a" or n}  # n/a only when some check is (card E-14)
    print("check: " + ", ".join(f"{n} {k}" for k, n in shown_counts.items())
          + (f"; {failed} invariants failed" if failed else "; invariants hold"))
    if a.write_recorded:
        kinds = {r["check_id"]: r["kind"] for r in expected}
        conditional = {r["check_id"] for r in expected if requires(r["expectation"])[0] is not None}
        try:
            with open(a.write_recorded, "x", newline="") as fh:
                fh.write("check_id\tvalue\n")
                for cid in sorted(got):
                    if kinds.get(cid) == "recorded" and not (cid in conditional and got[cid] == NOT_APPLICABLE):
                        fh.write(f"{cid}\t{got[cid]}\n")
        except OSError as e:
            print(f"check: --write-recorded: {e}", file=sys.stderr)
            return 2
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
