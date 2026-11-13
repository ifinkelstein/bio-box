"""The release check of the CRISPR-Cas example (`examples/crispr_cas_neighbourhoods/check.py` and `oracle.py`).

`check.py` compares a run's observed values with `expected.tsv` (invariants) and the release's `recorded/` file.
The files here are written by hand; the expected outcome of each case is stated beside it, not computed by the
code under test. `oracle.py` is checked on a hand-written `hmmsearch --tblout` and hits table.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("E-11")  # ring R1 by its directory (tests/conftest.py)

EX = Path(__file__).resolve().parents[2] / "examples" / "crispr_cas_neighbourhoods"
CHECK = EX / "check.py"
ORACLE = EX / "oracle.py"

EXPECTED = (
    "check_id\tasserts\tkind\texpectation\n"
    "zero\tnone of them\tinvariant\t== 0\n"
    "floor\tat least 27 bits\tinvariant\t>= 27\n"
    "ceiling\tat most 0.1 bit\tinvariant\t<= 0.1\n"
    "text\tthe order\tinvariant\t== tnsB ^ cas5f\n"
    "flag\ttrue or yes\tinvariant\tin true,yes\n"
    "count\ta count of one release\trecorded\trecorded\n"
    "other\tanother count\trecorded\trecorded\n"
)


def _write(path: Path, text: str) -> str:
    path.write_text(text)
    return str(path)


def _observed(tmp_path: Path, **values: str) -> str:
    base = {"zero": "0", "floor": "27.0", "ceiling": "0.10", "text": "tnsB ^ cas5f", "flag": "true",
            "count": "12", "other": "3"}
    base.update(values)
    rows = "".join(f"{k}\t{v}\n" for k, v in base.items() if v is not None)
    return _write(tmp_path / "observed.tsv", "check_id\tvalue\n" + rows)


def _run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-I", str(CHECK), *args], capture_output=True, text=True)


def _status(out: str) -> dict[str, str]:
    return {line.split("\t")[1]: line.split("\t")[0] for line in out.splitlines() if "\t" in line}


def test_everything_holds_and_matches_the_recorded_file(tmp_path):
    exp = _write(tmp_path / "expected.tsv", EXPECTED)
    rec = _write(tmp_path / "rX.tsv", "check_id\tvalue\ncount\t12\nother\t3.0\n")
    p = _run("--expected", exp, "--recorded", rec, _observed(tmp_path))
    assert p.returncode == 0, p.stdout + p.stderr
    # 27.0 >= 27, "0.10" <= 0.1 and 3 == 3.0 compare as numbers; the text compares as text
    assert _status(p.stdout) == {k: "ok" for k in ("zero", "floor", "ceiling", "text", "flag", "count", "other")}
    assert p.stdout.splitlines()[-1] == "check: 7 ok, 0 FAIL, 0 differs, 0 new, 0 missing, 0 unexpected; invariants hold"


@pytest.mark.parametrize("cid,value", [("zero", "1"), ("floor", "26.9"), ("ceiling", "0.2"),
                                       ("text", "tnsB ^ cas7f"), ("flag", "false"), ("floor", "n/a")])
def test_a_failed_invariant_fails_the_check(tmp_path, cid, value):
    exp = _write(tmp_path / "expected.tsv", EXPECTED)
    p = _run("--expected", exp, _observed(tmp_path, **{cid: value}))
    assert p.returncode == 1
    assert _status(p.stdout)[cid] == "FAIL"
    assert "1 invariants failed" in p.stdout.splitlines()[-1]


def test_a_missing_invariant_fails_and_a_missing_recorded_value_does_not(tmp_path):
    exp = _write(tmp_path / "expected.tsv", EXPECTED)
    p = _run("--expected", exp, _observed(tmp_path, zero=None))
    assert p.returncode == 1 and _status(p.stdout)["zero"] == "missing"
    p = _run("--expected", exp, _observed(tmp_path, count=None))
    assert p.returncode == 0 and _status(p.stdout)["count"] == "missing"


def test_a_recorded_value_that_differs_is_reported_not_failed(tmp_path):
    exp = _write(tmp_path / "expected.tsv", EXPECTED)
    rec = _write(tmp_path / "rX.tsv", "check_id\tvalue\ncount\t11\nother\t3\n")
    p = _run("--expected", exp, "--recorded", rec, _observed(tmp_path))
    assert p.returncode == 0
    assert _status(p.stdout)["count"] == "differs"
    assert "count\tobserved=12\texpected=11" in p.stdout


def test_a_release_without_a_recorded_file_reports_new_and_writes_one(tmp_path):
    exp = _write(tmp_path / "expected.tsv", EXPECTED)
    obs = _observed(tmp_path, extra=None)
    new = tmp_path / "recorded_new.tsv"
    p = _run("--expected", exp, "--recorded", str(tmp_path / "absent.tsv"), "--write-recorded", str(new), obs)
    assert p.returncode == 0
    assert _status(p.stdout)["count"] == "new" and _status(p.stdout)["other"] == "new"
    assert new.read_text() == "check_id\tvalue\ncount\t12\nother\t3\n"  # recorded checks only, by id
    p = _run("--expected", exp, "--write-recorded", str(new), obs)  # never overwrites
    assert p.returncode == 2


def test_an_observed_value_no_row_expects_is_listed_and_does_not_fail(tmp_path):
    exp = _write(tmp_path / "expected.tsv", EXPECTED)
    p = _run("--expected", exp, _observed(tmp_path, surprise="5"))
    assert p.returncode == 0
    assert _status(p.stdout)["surprise"] == "unexpected"


@pytest.mark.parametrize("expected", [
    EXPECTED.replace("== 0", "= 0"),                      # an unknown operator
    EXPECTED.replace(">= 27", ">= many"),                 # a bound that is not a number
    EXPECTED.replace("\trecorded\trecorded\n", "\trecorded\t== 3\n", 1),  # a recorded check with an expectation
    EXPECTED.replace("\tinvariant\t== 0", "\tsometimes\t== 0"),          # an unknown kind
    EXPECTED + "zero\tagain\tinvariant\t== 0\n",          # a check listed twice
    EXPECTED.replace("check_id\tasserts", "id\tasserts"),  # a wrong column line
])
def test_a_malformed_expected_file_is_refused(tmp_path, expected):
    exp = _write(tmp_path / "expected.tsv", expected)
    p = _run("--expected", exp, _observed(tmp_path))
    assert p.returncode == 2
    assert p.stdout == ""
    assert p.stderr.startswith("check: ")  # refused by check.py, not a crash


def test_a_check_observed_twice_is_refused(tmp_path):
    exp = _write(tmp_path / "expected.tsv", EXPECTED)
    a = _observed(tmp_path)
    b = _write(tmp_path / "more.tsv", "check_id\tvalue\nzero\t0\n")
    p = _run("--expected", exp, a, b)
    assert p.returncode == 2
    assert "observed twice" in p.stderr


def test_the_example_expected_file_is_well_formed_and_its_invariants_are_named():
    rows = [line.split("\t") for line in (EX / "expected.tsv").read_text().splitlines()]
    assert rows[0] == ["check_id", "asserts", "kind", "expectation"]
    kinds = {r[0]: r[2] for r in rows[1:]}
    assert len(kinds) == len(rows) - 1
    assert set(kinds.values()) == {"invariant", "recorded"}
    for cid in ("tniq_control_bits", "tniq_control_order", "tniq_cas7f_loci_with_cas1_cas2_cas3",
                "icp1_cas7f_window_has_cas1_and_cas3", "oracle_tniq_only_hmmsearch"):
        assert kinds[cid] == "invariant"
    assert all(r[3] != "== TBD" for r in rows[1:])


def test_the_oracle_compares_hits_with_hmmsearch_at_the_gathering_cutoff(tmp_path):
    hmm = _write(tmp_path / "a.hmm", "HMMER3/f [3.4 | Aug 2023]\nNAME  TniQ\nACC   PF06527.17\nGA    27.00 27.00;\n//\n")
    # target, accession, query name, accession, E-value, full-sequence score, ...
    tbl = _write(tmp_path / "a.tbl", "# comment\n"
                 "p1 - TniQ PF06527.17 1e-30 120.4 0.1 1e-29 119.0 0.1 1.0 1 1 0 1 1 1 1 -\n"
                 "p2 - TniQ PF06527.17 1e-5 27.0 0.1 1e-5 26.0 0.1 1.0 1 1 0 1 1 1 1 -\n"
                 "p3 - TniQ PF06527.17 1e-3 26.9 0.1 1e-3 26.0 0.1 1.0 1 1 0 1 1 1 1 -\n"
                 "p4 - TniQ PF06527.17 1e-6 30.0 0.1 1e-6 29.0 0.1 1.0 1 1 0 1 1 1 1 -\n")
    # p1 at 120.6 (0.2 bit off), p2 twice (two genes), p3 below GA in hmmsearch, p4 not a hit
    hits = _write(tmp_path / "hits.tsv", "accession\tprotein_uid\tseq_bitscore\n"
                  "A\tp1\t120.6\nA\tp2\t27.0\nB\tp2\t27.0\nC\tp3\t26.9\n")
    p = subprocess.run([sys.executable, "-I", str(ORACLE), hmm, tbl, f"tniq={hits}"], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    assert p.stdout.splitlines() == ["check_id\tvalue", "oracle_tniq_only_search\t1", "oracle_tniq_only_hmmsearch\t1",
                                     "oracle_tniq_max_bits_diff\t0.20"]
