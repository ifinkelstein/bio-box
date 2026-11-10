"""Second pass of the mutation pipeline (card D-01).

`mutmut run` attributes a test to a function only when the test calls it in-process during mutmut's
stats run. Tests that reach the code through a subprocess, or through a module-scoped fixture that
only the first test is charged for, are not attributed, so a mutant that mutmut leaves alive (or
marks "no tests") may still be killed by the suite. The verdict for those mutants comes from here:
each one is activated (MUTANT_UNDER_TEST, the mechanism mutmut's own trampoline reads, which child
processes inherit) in the mutated copy under ./mutants and the FULL selection runs against it.

    uv run mutmut run                                   # pass 1 (fast, per-mutant test selection)
    uv run python scripts/mutation_recheck.py recheck   # pass 2, resumable, writes the TSV
    uv run python scripts/mutation_recheck.py summary   # counts and the one-line diffs of survivors

TSV columns: mutant id, module, function, killed-by-full-run (yes/no), first failing test.
Run from the repository root with BIOBOX_TEST_RUN_DIR unset (this script unsets it for children).
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MUTANTS = REPO / "mutants"
# the same selection as [tool.mutmut] pytest_add_cli_args_test_selection in pyproject.toml
SELECTION = ["tests/r1", "tests/r2"]
# mutmut exit codes of the fast pass (mutmut.stats.status_by_exit_code): 0 alive, 5 and 33 no tests,
# 35 suspicious. Timeouts and segfaults count as detected and are not rechecked.
RECHECK_CODES = {0, 5, 33, 35}
STATUS = {0: "survived", 5: "no tests", 33: "no tests", 35: "suspicious"}
FAILED = re.compile(r"^(?:FAILED|ERROR) (\S+)", re.M)
TIMEOUT_S = 900


def fast_pass_alive(patterns: list[str]) -> dict[str, int]:
    """mutant id -> exit code, for every mutant the fast pass did not detect."""
    out: dict[str, int] = {}
    for meta in sorted(MUTANTS.rglob("*.meta")):
        for name, code in json.loads(meta.read_text())["exit_code_by_key"].items():
            if code in RECHECK_CODES and (not patterns or any(fnmatch.fnmatch(name, p) for p in patterns)):
                out[name] = code
    return out


def split_id(name: str) -> tuple[str, str]:
    module, _, rest = name.rpartition(".")
    return module, rest.rsplit("__mutmut_", 1)[0].removeprefix("x_")


def run_full(name: str) -> tuple[str, str]:
    env = dict(os.environ, MUTANT_UNDER_TEST=name, PY_IGNORE_IMPORTMISMATCH="1",
               PYTHONPATH=os.pathsep.join([str(MUTANTS / "src"), str(MUTANTS)]))
    env.pop("BIOBOX_TEST_RUN_DIR", None)
    cmd = [sys.executable, "-m", "pytest", "--rootdir=.", "--tb=native", "-x", "-q", "-rfE",
           "-p", "no:cacheprovider", "-p", "no:randomly", *SELECTION]
    try:
        p = subprocess.run(cmd, cwd=MUTANTS, env=env, capture_output=True, text=True, timeout=TIMEOUT_S)
    except subprocess.TimeoutExpired:
        return "yes", "TIMEOUT"
    if p.returncode == 0:
        return "no", ""
    if p.returncode in (1, 2, 3):
        m = FAILED.search(p.stdout)
        return "yes", m.group(1) if m else f"exit {p.returncode}"
    return "error", f"pytest exit {p.returncode}: {p.stdout[-200:].strip()!r}"


def recheck(args: argparse.Namespace) -> None:
    if not MUTANTS.is_dir():
        sys.exit("no ./mutants directory: run `uv run mutmut run` first")
    tsv = Path(args.tsv)
    done = set()
    if tsv.exists():
        done = {line.split("\t", 1)[0] for line in tsv.read_text().splitlines()[1:] if line}
    else:
        tsv.write_text("mutant\tmodule\tfunction\tkilled_by_full_run\tfirst_failing_test\n")
    todo = [n for n in fast_pass_alive(args.mutants) if n not in done]
    print(f"{len(done)} already in {tsv}, {len(todo)} to run", flush=True)

    def one(name: str) -> str:
        killed, test = run_full(name)
        module, func = split_id(name)
        return f"{name}\t{module}\t{func}\t{killed}\t{test}\n"

    with ThreadPoolExecutor(args.jobs) as ex, tsv.open("a") as f:
        for line in ex.map(one, todo):
            f.write(line)
            f.flush()
            print(line, end="", flush=True)


def one_line_diff(name: str) -> str:
    out = subprocess.run([sys.executable, "-m", "mutmut", "show", name], cwd=REPO, capture_output=True, text=True).stdout
    minus = [x[1:].strip() for x in out.splitlines() if x.startswith("-") and not x.startswith("---")]
    plus = [x[1:].strip() for x in out.splitlines() if x.startswith("+") and not x.startswith("+++")]
    return f"- {' / '.join(minus)}   + {' / '.join(plus)}"


def summary(args: argparse.Namespace) -> None:
    status: dict[str, int] = {}
    total = 0
    for meta in MUTANTS.rglob("*.meta"):
        for name, code in json.loads(meta.read_text())["exit_code_by_key"].items():
            if args.mutants and not any(fnmatch.fnmatch(name, p) for p in args.mutants):
                continue
            total += 1
            s = STATUS.get(code) or {1: "killed", 3: "killed"}.get(code) or ("timeout/segfault" if code in (36, 24, -24, 152, 255, -11, -9) else "not checked" if code is None else f"other({code})")
            status[s] = status.get(s, 0) + 1
    print(f"mutants {total}; fast pass: " + ", ".join(f"{k} {v}" for k, v in sorted(status.items())))
    tsv = Path(args.tsv)
    rows = [l.split("\t") for l in tsv.read_text().splitlines()[1:]] if tsv.exists() else []
    rows = [r for r in rows if not args.mutants or any(fnmatch.fnmatch(r[0], p) for p in args.mutants)]
    yes = [r for r in rows if r[3] == "yes"]
    no = [r for r in rows if r[3] == "no"]
    err = [r for r in rows if r[3] == "error"]
    print(f"recheck: {len(rows)} rows; killed by the full run {len(yes)}; final survivors {len(no)}; errors {len(err)}")
    for r in no + err:
        print(f"{r[0]}  [{r[3]}]  {one_line_diff(r[0]) if r[3] == 'no' else r[4]}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", default="mutation_recheck.tsv")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("recheck", help="run the full selection against each mutant the fast pass left alive")
    r.add_argument("mutants", nargs="*", help="optional fnmatch patterns on mutant ids, e.g. 'biobox.layout.*'")
    r.add_argument("--jobs", type=int, default=1)
    r.set_defaults(fn=recheck)
    s = sub.add_parser("summary", help="counts, and the diff of every final survivor")
    s.add_argument("mutants", nargs="*", help="optional fnmatch patterns on mutant ids")
    s.set_defaults(fn=summary)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
