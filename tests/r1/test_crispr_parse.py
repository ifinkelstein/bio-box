"""Ring R1: the CRISPR producer's parsing, coordinate conversion and origin rule, its table and its commands (card
E-12), without the MinCED binary (GitHub Actions has none; these run everywhere).

The GFF lines below are MinCED 0.4.2's own `-gffFull` output for chromosome 1 of *Photobacterium damselae* 04Ya311
(NZ_AP026780.1; run on wilkcomp02, 2026-11-13, `minced -gffFull GCF_026001825.1.fna pdam.txt pdam.gff` in the test
image), copied verbatim, plus a hand-made contig whose repeat and spacer positions are written down 1-based, as GFF3
counts, and checked 0-based half-open after the conversion.
"""

from __future__ import annotations

import pytest

from biobox import cli, ids
from biobox.l1 import crispr
from biobox.registry.schemas import TABLES, table

pytestmark = pytest.mark.covers("E-12")

PDAM = """##gff-version 3
NZ_AP026780.1\tminced:0.4.2\trepeat_region\t2224314\t2224582\t5\t.\t.\tID=CRISPR1;rpt_type=direct;rpt_family=CRISPR;rpt_unit_seq=ATTTTTAGCTACCTATGCGGCAGGACAC
NZ_AP026780.1\tminced:0.4.2\trepeat_unit\t2224314\t2224341\t1\t.\t.\tParent=CRISPR1;ID=DR.CRISPR1.1
NZ_AP026780.1\tminced:0.4.2\trepeat_unit\t2224374\t2224401\t1\t.\t.\tParent=CRISPR1;ID=DR.CRISPR1.2
NZ_AP026780.1\tminced:0.4.2\trepeat_unit\t2224434\t2224461\t1\t.\t.\tParent=CRISPR1;ID=DR.CRISPR1.3
NZ_AP026780.1\tminced:0.4.2\trepeat_unit\t2224495\t2224522\t1\t.\t.\tParent=CRISPR1;ID=DR.CRISPR1.4
NZ_AP026780.1\tminced:0.4.2\trepeat_unit\t2224555\t2224582\t1\t.\t.\tParent=CRISPR1;ID=DR.CRISPR1.5
##gff-version 3
"""

R = b"GTTCACTGCCGTACAGGCAGCTTAGAAA"  # 28
S1 = b"ACCGTTAGCATGCAAGTCCTGATCGGATACGTT"  # 33
S2 = b"CTTGAGGCATTACGCAACGTGTCAGTCCATAGG"  # 33
HAND = b"ACGTACGTACGT" + R + S1 + R + S2 + R + b"TTTTGGGG"
# 1-based inclusive, as a GFF3 reader counts: the first repeat begins at the 13th base
HAND_GFF = """##gff-version 3
h\tminced:0.4.2\trepeat_region\t13\t162\t3\t.\t.\tID=CRISPR1;rpt_type=direct
h\tminced:0.4.2\trepeat_unit\t13\t40\t1\t.\t.\tParent=CRISPR1;ID=DR.CRISPR1.1
h\tminced:0.4.2\trepeat_unit\t74\t101\t1\t.\t.\tParent=CRISPR1;ID=DR.CRISPR1.2
h\tminced:0.4.2\trepeat_unit\t135\t162\t1\t.\t.\tParent=CRISPR1;ID=DR.CRISPR1.3
"""


def test_the_parser_reads_minceds_own_gff_0_based_half_open():
    (a,) = list(crispr.parse_gff(PDAM.splitlines()))
    assert (a.record, a.start, a.end) == ("NZ_AP026780.1", 2224313, 2224582)
    assert a.units[0] == (2224313, 2224341) and a.units[-1] == (2224554, 2224582)
    assert len(a.units) == 5 and all(e - s == 28 for s, e in a.units)
    els = crispr.elements(a.units)
    assert [e[0] for e in els] == ["repeat", "spacer"] * 4 + ["repeat"]
    assert [e[1] for e in els] == list(range(9))
    # MinCED's text report gives these spacers 32, 32, 33 and 32 bases (the same run)
    assert [e[3] - e[2] for e in els if e[0] == "spacer"] == [32, 32, 33, 32]


@pytest.mark.covers("C-05")
def test_a_hand_made_contig_converts_to_the_planted_bases():
    assert len(HAND) == 170
    (a,) = list(crispr.parse_gff(HAND_GFF.splitlines()))
    assert (a.start, a.end) == (12, len(HAND) - 8)
    got = [(k, HAND[s:e]) for k, _, s, e in crispr.elements(a.units)]
    assert got == [("repeat", R), ("spacer", S1), ("repeat", R), ("spacer", S2), ("repeat", R)]


@pytest.mark.parametrize("bad, why", [
    ("\tParent=CRISPR1;", "unit of another parent"),
    ("units 3", "score"),
    ("overlap", "overlap"),
    ("range", "range"),
])
def test_the_parser_refuses_malformed_arrays(bad, why):
    lines = HAND_GFF.splitlines()
    if bad.startswith("\tParent"):
        lines[2] = lines[2].replace("Parent=CRISPR1", "Parent=CRISPR9")
    elif bad == "units 3":
        lines[1] = lines[1].replace("\t3\t.", "\t4\t.")
    elif bad == "overlap":
        lines[3] = lines[3].replace("\t74\t101", "\t30\t57")
    else:
        lines[1] = lines[1].replace("\t13\t", "\t12\t")
    with pytest.raises(ValueError):
        list(crispr.parse_gff(lines))


def _arr(start: int, end: int, n_units: int = 3) -> crispr.Array:
    step = (end - start) // n_units
    units = tuple((start + i * step, start + i * step + 28) for i in range(n_units - 1)) + ((end - 28, end),)
    return crispr.Array("r", start, end, units)


@pytest.mark.covers("C-05", "MR-04")
def test_the_leftmost_base_rule_finds_an_array_once_with_end_past_the_origin():
    n, f = 10_000, 10_000  # a contig shorter than CIRC_FLANK: F = n
    rec = crispr.Rec(b"\0" * 16, n, f, True, n + 2 * f)
    # an array across the origin: in the record once with its leftmost base on the contig, once in the left flank
    assert crispr.place(_arr(f + n - 50, f + n + 150), rec) == (n - 50, n + 150)
    assert crispr.place(_arr(f - 50, f + 150), rec) is None
    assert crispr.place(_arr(f + n + 10, f + n + 200), rec) is None  # its copy in the right flank
    assert crispr.place(_arr(f + 10, f + 300), rec) == (10, 300)
    # longer than the circle: impossible on the molecule
    assert crispr.place(_arr(f + 5, f + n + 50, 60), rec) == crispr.REJECT_CIRCLE
    # a linear record: as reported
    lin = crispr.Rec(b"\0" * 16, n, 0, False, n)
    assert crispr.place(_arr(0, 200), lin) == (0, 200)


def test_an_array_a_long_circular_records_end_may_have_cut_stops_the_run():
    n, f = 1_000_000, crispr.CIRC_FLANK
    rec = crispr.Rec(b"\0" * 16, n, f, True, n + 2 * f)
    with pytest.raises(RuntimeError, match="CIRC_FLANK"):
        crispr.place(_arr(f + n - 10, n + 2 * f - 5, 400), rec)


def test_extended_records_hold_the_flanks_around_the_contig():
    seq = b"ABCDEFGHIJ"
    assert crispr.extended(seq, 0) == seq
    assert crispr.extended(seq, 3) == b"HIJ" + seq + b"ABC"
    assert crispr.extended(seq, 10) == seq * 3


def test_parameters_are_minceds_defaults_and_named_in_the_producer_tag():
    assert dict(crispr.MINCED_ARGS) == {"-searchWL": 8, "-minNR": 3, "-minRL": 23, "-maxRL": 47, "-minSL": 26,
                                         "-maxSL": 50}
    assert crispr.PARAM_TAG == "w8-nr3-r23x47-s26x50"
    assert crispr.MAX_UNIT == 47 + 50
    assert crispr.KIND_CRISPR == ids.KIND["crispr_array"] == 4 and crispr.STRAND == 1


@pytest.mark.covers("C-01", "LY-01", "LY-09")
def test_the_l1_crispr_table_is_registered_with_a_stable_key_and_no_ordinal_in_it():
    t = table("l1_crispr")
    assert (t.location, t.layer, t.partition_keys) == ("l1/crispr", "L1", ("producer", "stage"))
    assert t.schema.names == ["feature_uid", "contig_uid", "region_uid", "element", "ordinal", "start", "end",
                              "strand", "seq", "contig_len"]
    assert t.primary_key == ("feature_uid", "region_uid", "producer") and "ordinal" not in t.primary_key
    assert t.sort_key == ("contig_uid", "feature_uid", "ordinal")
    assert t.enums == {"strand": (1, -1), "element": ("repeat", "spacer")}
    assert not any(f.nullable for f in t.schema)
    assert "l1_crispr" in TABLES and not t.small_ok


def test_the_commands_parse():
    import argparse

    seen = {}
    real = argparse.ArgumentParser.parse_args

    def capture(self, args=None, namespace=None):
        ns = real(self, args, namespace)
        seen["ns"] = ns
        raise SystemExit(0)

    mp = pytest.MonkeyPatch()
    mp.setattr(argparse.ArgumentParser, "parse_args", capture)
    try:
        with pytest.raises(SystemExit):
            cli.main(["stage", "crispr", "--root", "/r", "--stage", "S0", "--cpus", "4", "--tmp", "/t"])
        assert (seen["ns"].step, seen["ns"].cpus) == ("crispr", 4)
        with pytest.raises(SystemExit):
            cli.main(["build-stage", "--profile", "p", "--root", "/r", "--stage", "S0", "--release", "R",
                      "--genomes", "g", "--objects", "o", "--crispr"])
        assert seen["ns"].crispr is True
    finally:
        mp.undo()


def test_the_build_stage_estimate_and_refusal_threshold():
    from biobox import stage as bs

    assert bs.crispr_cpu_min(int(1e9)) == pytest.approx(bs.CRISPR_CPU_S_PER_GBP / 60)
    assert bs.CRISPR_CPU_S_PER_GBP > 0  # measured on stage P1 (planning/measurements/2026-11-13_crispr_arrays.md)
    assert bs.Build.__dataclass_fields__["crispr"].default is False  # off by default, as the ncRNA step is
