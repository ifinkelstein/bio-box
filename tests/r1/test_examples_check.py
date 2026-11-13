"""The release check of the worked examples: the shared checker (`examples/_shared/check.py`, card E-11; it was
the CRISPR-Cas example's own `check.py` until then) and each example's oracle and expectations.

`check.py` compares a run's observed values with `expected.tsv` (invariants) and the release's `recorded/` file.
The files here are written by hand; the expected outcome of each case is stated beside it, not computed by the
code under test. The oracles (`examples/crispr_cas_neighbourhoods/oracle.py`, `examples/_shared/oracle.py`) are
checked on a hand-written `hmmsearch --tblout` and hits table.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("E-11")  # ring R1 by its directory (tests/conftest.py)

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
EX = EXAMPLES / "crispr_cas_neighbourhoods"
SHARED = EXAMPLES / "_shared"
CHECK = SHARED / "check.py"
ORACLE = EX / "oracle.py"
SHARED_ORACLE = SHARED / "oracle.py"
RNAP = EXAMPLES / "rnap"

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


# --- card E-11: the checker is shared, and writes a table for the release summary --------------------------------


def test_the_checker_is_shared_and_the_crispr_example_has_no_copy_of_its_own():
    assert CHECK.is_file()
    assert not (EX / "check.py").exists()


def test_the_table_has_one_row_per_check_with_its_kind_expectation_observation_and_verdict(tmp_path):
    exp = _write(tmp_path / "expected.tsv", EXPECTED)
    rec = _write(tmp_path / "rX.tsv", "check_id\tvalue\ncount\t11\n")
    table = tmp_path / "table.tsv"
    p = _run("--expected", exp, "--recorded", rec, "--table", str(table),
             _observed(tmp_path, floor="26.9", zero=None, surprise="5"))
    assert p.returncode == 1
    rows = [line.split("\t") for line in table.read_text().splitlines()]
    assert rows[0] == ["check", "kind", "expected", "observed", "verdict"]
    assert rows[1:] == [
        ["zero", "invariant", "== 0", "-", "missing"],
        ["floor", "invariant", ">= 27", "26.9", "FAIL"],
        ["ceiling", "invariant", "<= 0.1", "0.10", "ok"],
        ["text", "invariant", "== tnsB ^ cas5f", "tnsB ^ cas5f", "ok"],
        ["flag", "invariant", "in true,yes", "true", "ok"],
        ["count", "recorded", "11", "12", "differs"],
        ["other", "recorded", "-", "3", "new"],
        ["surprise", "-", "-", "5", "unexpected"],
    ]
    p = _run("--expected", exp, "--table", str(table), _observed(tmp_path))  # never overwrites
    assert p.returncode == 2 and p.stderr.startswith("check: ")


# --- the shared oracle: a run's hits against hmmsearch at a stated cutoff ------------------------------------------

HMM_TEXT = "HMMER3/f [3.4 | Aug 2023]\nNAME  RNA_pol\nACC   PF00940.25\nGA    22.00 22.00;\n//\n"
TBL_TEXT = ("# comment\n"
            "p1 - RNA_pol PF00940.25 1e-90 367.7 0.1 1e-89 360.0 0.1 1.0 1 1 0 1 1 1 1 -\n"
            "p2 - RNA_pol PF00940.25 1e-5 30.0 0.1 1e-5 29.0 0.1 1.0 1 1 0 1 1 1 1 -\n"
            "p3 - RNA_pol PF00940.25 1e-3 25.0 0.1 1e-3 24.0 0.1 1.0 1 1 0 1 1 1 1 -\n"
            "p4 - RNA_pol PF00940.25 1e-6 31.0 0.1 1e-6 30.0 0.1 1.0 1 1 0 1 1 1 1 -\n"
            "q1 - Other PF00001.1 1e-6 99.0 0.1 1e-6 99.0 0.1 1.0 1 1 0 1 1 1 1 -\n")
# p1 at 367.9 in the run (0.2 bit off), p2 on two genes, p3 below 30 in hmmsearch but above GA, p4 not a hit, q1
# another model's row
HITS_TEXT = "accession\tprotein_uid\tseq_bitscore\nA\tp1\t367.9\nA\tp2\t30.0\nB\tp2\t30.0\nC\tp3\t25.0\n"


def _oracle(tmp_path: Path, bits: str, model: str = "RNA_pol") -> subprocess.CompletedProcess:
    hmm = _write(tmp_path / "a.hmm", HMM_TEXT)
    tbl = _write(tmp_path / "a.tbl", TBL_TEXT)
    hits = _write(tmp_path / "hits.tsv", HITS_TEXT)
    return subprocess.run([sys.executable, "-I", str(SHARED_ORACLE), "--hmm", hmm, "--tblout", tbl, "--prefix",
                           "rnap", "--model", model, "--bits", bits, hits], capture_output=True, text=True)


def test_the_shared_oracle_compares_hits_with_hmmsearch_at_a_cutoff_in_bits(tmp_path):
    p = _oracle(tmp_path, "30")
    assert p.returncode == 0, p.stderr
    # at 30 bits: p3 is a hit hmmsearch scores below the cutoff; p4 scores above it and is not a hit
    assert p.stdout.splitlines() == ["check_id\tvalue", "oracle_rnap_only_search\t1", "oracle_rnap_only_hmmsearch\t1",
                                     "oracle_rnap_max_bits_diff\t0.20"]


def test_the_shared_oracle_takes_the_gathering_cutoff_from_the_model(tmp_path):
    p = _oracle(tmp_path, "ga")
    assert p.returncode == 0, p.stderr
    # at GA 22: p3 is above it in hmmsearch too; p4 is still not a hit
    assert p.stdout.splitlines() == ["check_id\tvalue", "oracle_rnap_only_search\t0", "oracle_rnap_only_hmmsearch\t1",
                                     "oracle_rnap_max_bits_diff\t0.20"]


@pytest.mark.parametrize("model,bits", [("TniQ", "ga"), ("RNA_pol", "many")])
def test_the_shared_oracle_refuses_a_model_it_lacks_or_a_cutoff_that_is_not_one(tmp_path, model, bits):
    p = _oracle(tmp_path, bits, model)
    assert p.returncode == 2 and p.stdout == ""
    assert p.stderr.startswith("oracle: ") or "usage" in p.stderr


# --- the RNA polymerase example -------------------------------------------------------------------------------------


def test_the_rnap_expected_file_is_well_formed_and_names_its_invariants():
    rows = [line.split("\t") for line in (RNAP / "expected.tsv").read_text().splitlines()]
    assert rows[0] == ["check_id", "asserts", "kind", "expectation"]
    assert all(len(r) == 4 for r in rows)
    kinds = {r[0]: r[2] for r in rows[1:]}
    assert len(kinds) == len(rows) - 1
    assert set(kinds.values()) == {"invariant", "recorded"}
    for cid in ("t7_control_bits", "t7_control_architecture", "t7_control_order", "t7_window_is_the_contig_order",
                "coord_t7_kinase_offset", "coord_t7_ligase_offset", "rnap_negative_control_hits",
                "rnap_loci_virus", "rnap_loci_not_virus", "oracle_rnap_only_search", "oracle_rnap_only_hmmsearch",
                "oracle_rnap_max_bits_diff", "rnap_flank_loci_missing"):
        assert kinds[cid] == "invariant", cid
    for cid in ("rnap_hit_proteins", "rnap_hit_genes", "rnap_flank_groups"):
        assert kinds[cid] == "recorded", cid
    assert all(r[3] != "== TBD" for r in rows[1:])


def test_the_rnap_recorded_file_holds_exactly_the_recorded_checks():
    kinds = {r.split("\t")[0]: r.split("\t")[2] for r in (RNAP / "expected.tsv").read_text().splitlines()[1:]}
    rows = [line.split("\t") for line in (RNAP / "recorded" / "rN.3.tsv").read_text().splitlines()]
    assert rows[0] == ["check_id", "value"]
    assert {r[0] for r in rows[1:]} == {c for c, k in kinds.items() if k == "recorded"}
    assert [r[0] for r in rows[1:]] == sorted(r[0] for r in rows[1:])  # as check.py --write-recorded writes it


@pytest.mark.parametrize("example", ["crispr_cas_neighbourhoods", "rnap"])
def test_each_example_passes_its_check_on_its_own_recorded_run(example):
    """The observed values of the recorded run (`results/observed.tsv`) meet every invariant and equal the recorded
    file: the committed files agree with each other."""
    ex = EXAMPLES / example
    p = _run("--expected", str(ex / "expected.tsv"), "--recorded", str(ex / "recorded" / "rN.3.tsv"),
             str(ex / "results" / "observed.tsv"))
    assert p.returncode == 0, p.stdout + p.stderr
    assert p.stdout.splitlines()[-1].endswith("0 FAIL, 0 differs, 0 new, 0 missing, 0 unexpected; invariants hold")


NEIGHBORS_TEXT = (
    '# schema\t"nbr_summary_header/1"\n# n_anchor_loci\t4\n# n_contigs\t3\n# n_censored_loci\t1\n'
    '# n_neighbors\t37\n# n_keys\t{"pfam":3,"nr50":5}\n'
    "\n"
    "level\tkey\tunannotated\tlibrary\tmodel_acc\tmodel_name\tn_loci\tshare_of_loci\tn_neighbors\tn_genomes\t"
    "genome_basis\toffsets\tmodal_offset\tsame_strand_share\n"
    "pfam\tPF01068.27\tfalse\tpfam-38.2\tPF01068.27\tDNA_ligase_A_M\t2\t0.5\t2\t2\tgenome_uid\t3:2\t3\t1.0\n"
    "pfam\t\ttrue\tpfam-38.2\t\t\t4\t1.0\t20\t3\tgenome_uid\t1:4\t1\t1.0\n"
    "pfam\tPF08200.18\tfalse\tpfam-38.2\tPF08200.18\tPhage_T7_1_1\t3\t0.75\t3\t3\tgenome_uid\t1:3\t1\t1.0\n"
    "\n"
    "level\tkey\tunannotated\trep_uid\trep_length\tmodel_name\tn_loci\tshare_of_loci\tn_neighbors\tn_genomes\t"
    "genome_basis\toffsets\tmodal_offset\tsame_strand_share\n"
    "nr50\taa\tfalse\tbb\t90\tPhage_T7_1_1\t9\t1.0\t9\t3\tgenome_uid\t1:9\t1\t1.0\n"
)
GROUPS_TEXT = ('# release\t"rX"\n# n_groups\t2\n'
               "group\tn_loci\tn_leaves\tn_leaves_in_tree\tsampled\tn_context_types\n"
               "g1\t3\t3\t3\tfalse\t2\ng2\t1\t1\t1\tfalse\t1\n")
CONTEXT_TEXT = ("feature_uid\tgroup\tcontext_type\tleaf\tin_tree\tkind\n"
                "f1\tg1\t1\tf1\ttrue\tvirus\nf2\tg1\t1\tf1\ttrue\tprokaryote\nf5\tg1\t2\tf5\ttrue\tvirus\n"
                "f3\tg2\t1\tf3\ttrue\tvirus\n")
HITS_GENES_TEXT = "accession\tfeature_uid\tprotein_uid\nA\tf1\tp1\nA\tf2\tp1\nB\tf3\tp2\nC\tf4\tp3\nC\tf5\tp3\n"


def test_the_rnap_observer_reads_the_summary_the_flanks_and_the_model(tmp_path):
    """`examples/rnap/observe.py`: the values the release check takes from the files of the run (not from SQL).
    Expected values are counted by hand from the files above."""
    flanks = tmp_path / "flanks"
    flanks.mkdir()
    _write(flanks / "groups.tsv", GROUPS_TEXT)
    _write(flanks / "context_types.tsv", CONTEXT_TEXT)
    hmm = _write(tmp_path / "m.hmm", "abc")
    p = subprocess.run([sys.executable, "-I", str(RNAP / "observe.py"), "--hmm", hmm,
                        "--hits", _write(tmp_path / "hits.tsv", HITS_GENES_TEXT),
                        "--neighbors", _write(tmp_path / "neighbors.tsv", NEIGHBORS_TEXT),
                        "--flanks", str(flanks), "--control", "f1"], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    lines = p.stdout.splitlines()
    assert lines[0] == "check_id\tvalue"
    got = dict(line.split("\t") for line in lines[1:])
    assert len(got) == len(lines) - 1
    assert got == {
        "rnap_hmm_sha256": "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",  # sha256("abc")
        "rnap_summary_loci": "4", "rnap_summary_contigs": "3", "rnap_summary_censored_loci": "1",
        "rnap_summary_neighbours": "37",
        "rnap_top_pfam_neighbour": "Phage_T7_1_1", "rnap_top_pfam_neighbour_loci": "3",  # the unannotated row is not one
        "rnap_flank_groups": "2", "rnap_flank_leaves": "4", "rnap_flank_context_types": "3",
        "rnap_flank_types_mixed_kind": "1",  # g1 type 1 holds a virus and a prokaryote locus
        "rnap_flank_loci_missing": "1",  # f4 is a hit gene with no row in context_types.tsv
        "t7_flank_group_loci": "3", "t7_flank_group_types": "2", "t7_flank_type_kinds": "prokaryote,virus",
    }


def test_the_rnap_observer_leaves_out_the_control_rows_when_the_flanks_do_not_hold_it(tmp_path):
    flanks = tmp_path / "flanks"
    flanks.mkdir()
    _write(flanks / "groups.tsv", GROUPS_TEXT)
    _write(flanks / "context_types.tsv", CONTEXT_TEXT)
    p = subprocess.run([sys.executable, "-I", str(RNAP / "observe.py"), "--hmm", _write(tmp_path / "m.hmm", "abc"),
                        "--hits", _write(tmp_path / "hits.tsv", HITS_GENES_TEXT),
                        "--neighbors", _write(tmp_path / "neighbors.tsv", NEIGHBORS_TEXT),
                        "--flanks", str(flanks), "--control", "f9"], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr  # the control's rows are left out: check.py then reports them missing
    got = dict(line.split("\t") for line in p.stdout.splitlines()[1:])
    assert "t7_flank_group_loci" not in got and got["rnap_flank_groups"] == "2"
