"""Card C-52 (DEC-63): `scripts/stress/build_stage_pod/c52_pl23.py` makes stages PL2 and PL3 from card C-40's staging.

Every rule is tried on constructed inputs. The oracle for a terminal repeat is the sequence the test built (a core
`C` and a repeat `R` written at both ends), never the script's own finding; the oracle for the split is the objects'
`redistributable` flag the test wrote; for the deferred list, the rows of the table the test wrote.
"""

from __future__ import annotations

import importlib.util
import json
import random
import subprocess
import sys
from pathlib import Path

import pytest

from biobox import ids

pytestmark = pytest.mark.covers("C-52")  # ring R1 by its directory (tests/conftest.py)

BSP = Path(__file__).resolve().parents[2] / "scripts" / "stress" / "build_stage_pod"
SCRIPT = BSP / "c52_pl23.py"


def _module():
    spec = importlib.util.spec_from_file_location("c52_pl23", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c52_pl23"] = mod  # a dataclass looks its module up by name
    spec.loader.exec_module(mod)
    return mod


def _seq(n: int, seed: int) -> bytes:
    """Random ACGT of length n that neither starts nor ends with A (so an A run at the ends stays its own length)."""
    rng = random.Random(seed)
    s = bytearray(rng.choice(b"ACGT") for _ in range(n))
    s[0], s[-1] = ord("C"), ord("G")
    return bytes(s)


C = _seq(3000, 1)
R = _seq(57, 2)


def _rot(s: bytes) -> set[bytes]:
    return {s[i:] + s[:i] for i in range(len(s))}


# --------------------------------------------------------------------------------------------------------------
# the two ends of a record


def test_an_exact_terminal_repeat_is_found_and_its_trim_leaves_the_circle_once():
    m = _module()
    seq = R + C + R
    e = m.ends(seq)
    assert e.exact == len(R)
    assert m.classify(e) == "exact"
    trimmed = m.trim(seq, e.exact)
    assert trimmed == R + C
    # the circle the source meant (R once, then C) is the trimmed record read round from any start
    assert C + R in _rot(trimmed)


def test_case_does_not_hide_a_repeat():
    m = _module()
    e = m.ends(R + C + R.lower())
    assert e.exact == len(R) and m.classify(e) == "exact"


def test_a_record_whose_ends_do_not_repeat_is_none():
    m = _module()
    e = m.ends(C)
    assert (e.exact, e.near, e.far) == (0, 0, 0)
    assert m.classify(e) == "none"


def test_a_record_shorter_than_two_seeds_is_none():
    m = _module()
    e = m.ends(R[:30])
    assert e.exact == 0 and m.classify(e) == "none"


def test_ends_that_agree_only_in_part_are_partial_with_the_span_and_mismatches():
    m = _module()
    r2 = bytearray(R)
    r2[40] = ord("A") if r2[40] != ord("A") else ord("C")  # one mismatch in the second copy, after the seed
    e = m.ends(R + C + bytes(r2))
    assert e.exact == 0
    assert e.near >= 1 and (e.near_len, e.near_mm, e.near_run) == (len(R), 1, 40)
    assert m.classify(e) == "partial"


def test_a_second_copy_with_a_changed_start_is_partial_from_the_other_end():
    m = _module()
    r2 = bytearray(R)
    r2[3] = ord("A") if r2[3] != ord("A") else ord("C")  # the start seed of the second copy differs
    e = m.ends(R + C + bytes(r2))
    assert e.exact == 0 and e.near >= 1 and (e.near_len, e.near_mm, e.near_run) == (len(R), 1, len(R) - 4)
    assert m.classify(e) == "partial"


def test_a_second_copy_broken_by_an_indel_is_partial():
    m = _module()
    t = _seq(300, 6)
    t2 = t[:200] + t[201:]  # the second copy lacks one base 200 bases in: an unbroken run of 200 over a span of 299
    e = m.ends(t + C + t2)
    assert e.exact == 0 and e.near >= 1 and e.near_run >= 200
    assert m.classify(e) == "partial"


def test_a_longer_near_copy_beside_a_shorter_exact_repeat_is_partial():
    m = _module()
    d = _seq(30, 7)
    d2 = bytearray(d)
    d2[10] = ord("A") if d2[10] != ord("A") else ord("C")
    # the ends are R D R ... R D' R: R is an exact terminal repeat; R D R against R D' R is one with a substitution
    e = m.ends(R + d + R + C + R + bytes(d2) + R)
    assert e.exact == len(R)
    assert e.near >= 1 and (e.near_len, e.near_mm) == (2 * len(R) + len(d), 1)
    assert m.classify(e) == "partial"


def test_a_shared_stretch_that_then_differs_is_far_and_leaves_the_class():
    m = _module()
    # a copy of R's first 30 bases 500 bases before the end, then other bases: an internal repeat, not a terminal one
    tail = R[:30] + _seq(470, 3)
    e = m.ends(R + C + tail + R)
    assert e.exact == len(R) and (e.near, e.far) == (0, 1)
    assert m.classify(e) == "exact"
    e = m.ends(R + C + tail)
    assert (e.exact, e.near, e.far) == (0, 0, 1)
    assert m.classify(e) == "none"


def test_a_long_repeat_is_its_own_class():
    m = _module()
    big = _seq(2500, 4)
    e = m.ends(big + C + big)
    assert e.exact == 2500
    assert m.classify(e, long_bound=1000) == "long"
    assert m.classify(e, long_bound=5000) == "exact"


@pytest.mark.parametrize("run", [b"A" * 30, b"N" * 25, b"AT" * 15])
def test_a_low_complexity_repeat_is_its_own_class(run):
    m = _module()
    e = m.ends(run + C + run)
    assert e.exact >= len(run)
    assert m.classify(e) == "low_complexity"


def test_a_repeat_longer_than_half_the_record_is_not_a_terminal_repeat():
    m = _module()
    # (XY)(XY)(X): the longest border is longer than half the record; the rule looks only up to half
    x = _seq(40, 5)
    e = m.ends(x + x + x)
    assert e.exact <= len(x + x + x) // 2


# --------------------------------------------------------------------------------------------------------------
# the sources' flags


@pytest.mark.parametrize("source,value,circ", [
    ("gut_plasmids_yu2024", "True", True),
    ("gut_plasmids_yu2024", "False", False),
    ("plasmidscope_IMG-PR", "direct terminal repeat", True),
    ("plasmidscope_IMG-PR", "linear", False),
    ("plasmidscope_IMG-PR", "concatemer", False),
    ("plasmidscope_IMG-PR", "inverted terminal repeat", False),
    ("plasmidscope_mMGE", "circular", True),
    ("plasmidscope_RefSeq", "circular", True),
    ("plasmidscope_GenBank", "linear", False),
])
def test_the_source_flag_of_each_set(source, value, circ):
    assert _module().source_flag(source, value) is circ


def test_a_source_or_value_the_rule_does_not_know_is_refused():
    m = _module()
    with pytest.raises(m.Refused):
        m.source_flag("someone_else", "circular")
    with pytest.raises(m.Refused):
        m.source_flag("plasmidscope_RefSeq", "round")


# --------------------------------------------------------------------------------------------------------------
# a constructed staging in card C-40's form

REC_HEAD = ["ordinal", "source", "accession", "contig_uid", "length", "topology_value", "completeness_value",
            "circular", "assembly_quality", "fragment", "first_source", "first_accession", "owner_stage"]


def _staging(tmp: Path, sources: list[tuple[str, bool, list[tuple[str, bytes, str, str]]]], owners=None) -> Path:
    """sources: (name, redistributable, [(accession, seq, topology_value, completeness_value)]), in staging order.
    Writes the FASTA objects, objects.json, records.tsv and skipped_records.tsv as make_inputs.py plasmid-sources
    does (first record of a sequence, owner stage from `owners`)."""
    d = tmp / "PL2"
    d.mkdir()
    objs, rows, first = [], [], {}
    owners = owners or {}
    k = 0
    for name, red, recs in sources:
        with open(d / f"{name}.fasta", "wb") as w:
            for acc, seq, tv, cv in recs:
                w.write(b">" + acc.encode() + b" some description\n" + seq + b"\n")
                cu = ids.contig_uid(seq)
                frag = cv == "incomplete"
                circ = (tv in ("circular", "True")) and not frag
                q = "complete" if cv in ("complete", "True") else "contig"
                f0 = first.get(cu)
                rows.append([k, name, acc, cu.hex(), len(seq), tv, cv, str(circ).lower(), q, str(frag).lower(),
                             "" if f0 is None else f0[0], "" if f0 is None else f0[1], owners.get(acc, "")])
                first.setdefault(cu, (name, acc))
                k += 1
        objs.append({"path": str(d / f"{name}.fasta"), "source_name": name.split("_")[0], "source_release": "r1",
                     "license": "L", "redistributable": red, "retrieved_at": "2026-11-17T00:00:00Z"})
    (d / "objects.json").write_text(json.dumps(objs, indent=1))
    (d / "records.tsv").write_text("".join("\t".join(map(str, r)) + "\n" for r in [REC_HEAD, *rows]))
    (d / "skipped_records.tsv").write_text("source\taccession\tcontig_uid\tlength\treason\n")
    return d


def _run(*args, check=True):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], capture_output=True, text=True, check=check)


def _tsv(path: Path) -> list[dict]:
    lines = path.read_text().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]


def _world(tmp: Path) -> Path:
    big = _seq(1500, 9)
    return _staging(tmp, [
        ("gut_plasmids_yu2024", True, [
            ("g1", R + C + R, "True", "True"),           # circular, exact repeat
            ("g2", C[::-1], "True", "True"),             # circular, ends do not repeat
            ("g3", R + C[:900] + R, "False", "False"),  # linear: not looked at
        ]),
        ("plasmidscope_IMG-PR", False, [
            ("i1", R + C[:1200] + R, "direct terminal repeat", "complete"),
            ("i2", big + C[:800] + big, "direct terminal repeat", "complete"),   # a long repeat
            ("i3", C[:1500], "linear", "incomplete"),
        ]),
        ("plasmidscope_RefSeq", True, [
            ("r1", C[100:2100], "circular", "complete"),
            ("r2", R + C[:700] + R, "circular", "incomplete"),  # flagged circular and incomplete
        ]),
    ])


def test_survey_reads_every_flagged_record_and_counts_its_class(tmp_path):
    staging = _world(tmp_path)
    out = tmp_path / "survey"
    _run("survey", "--inputs", staging, "--out", out, "--long-bound", "1000")
    rows = {r["accession"]: r for r in _tsv(out / "survey.tsv")}
    assert set(rows) == {"g1", "g2", "i1", "i2", "r1", "r2"}  # the flagged ones only
    assert rows["g1"]["class"] == "exact" and int(rows["g1"]["exact"]) == len(R)
    assert rows["g2"]["class"] == "none"
    assert rows["i1"]["class"] == "exact"
    assert rows["i2"]["class"] == "long" and int(rows["i2"]["exact"]) == 1500
    assert rows["r1"]["class"] == "none"
    assert rows["r2"]["fragment"] == "true"
    s = json.loads((out / "summary.json").read_text())
    assert s["sources"]["gut_plasmids_yu2024"]["records"] == 3
    assert s["sources"]["gut_plasmids_yu2024"]["flagged"] == 2
    assert s["sources"]["plasmidscope_IMG-PR"]["classes"] == {"exact": 1, "long": 1}
    assert s["sources"]["plasmidscope_RefSeq"]["flagged_and_fragment"] == 1
    assert s["sources"]["plasmidscope_IMG-PR"]["examples"]["long"] == ["i2"]


def test_survey_refuses_records_that_do_not_line_up_with_the_objects(tmp_path):
    staging = _world(tmp_path)
    lines = (staging / "records.tsv").read_text().splitlines(keepends=True)
    lines[1], lines[2] = lines[2], lines[1]
    (staging / "records.tsv").write_text("".join(lines))
    got = _run("survey", "--inputs", staging, "--out", tmp_path / "survey", check=False)
    assert got.returncode == 2 and "line up" in got.stderr
    assert not (tmp_path / "survey" / "summary.json").exists()


# --------------------------------------------------------------------------------------------------------------
# the two stagings

CG, CI, CM, CR = _seq(2000, 21), _seq(1800, 22), _seq(1600, 23), _seq(2200, 24)


def _w2(tmp: Path, extra_refseq=(), owners=None, extra_img=()) -> Path:
    return _staging(tmp, [
        ("gut_plasmids_yu2024", True, [
            ("g1", R + CG + R, "True", "True"),
            ("g2", CG[::-1], "True", "True"),
            ("g3", R + CG[:900] + R, "False", "False"),
        ]),
        ("plasmidscope_IMG-PR", False, [
            ("i1", R + CI + R, "direct terminal repeat", "complete"),
            ("i1b", R + CI, "direct terminal repeat", "complete"),     # i1 as its source might have trimmed it
            ("i2", R + CI + R, "direct terminal repeat", "complete"),  # i1's bytes again
            ("i3", CI[:1500], "linear", "incomplete"),
            *extra_img,
        ]),
        ("plasmidscope_mMGE", False, [("m1", CM, "circular", "complete")]),
        ("plasmidscope_RefSeq", True, [("r1", CR, "circular", "complete"), ("r2", CR[:1000], "linear", "-"),
                                       *extra_refseq]),
    ], owners=owners)


def _root(tmp: Path, held: dict[str, list[bytes]]) -> Path:
    import pyarrow as pa
    import pyarrow.parquet as pq

    root = tmp / "root"
    for st, uids in held.items():
        d = root / "l0seq" / "contig_index" / f"stage={st}"
        d.mkdir(parents=True)
        pq.write_table(pa.table({"contig_uid": pa.array(uids, pa.binary(16)),
                                 "is_default_variant": pa.array([True] * len(uids))}), d / "part-0.parquet")
    root.mkdir(exist_ok=True)
    return root


def _stage(tmp: Path, staging: Path, root: Path, *extra, check=True):
    _run("survey", "--inputs", staging, "--out", tmp / "survey")
    return _run("stage", "--inputs", staging, "--survey", tmp / "survey", "--root", root,
                "--out-pl2", tmp / "out" / "PL2.new", "--final-pl2", tmp / "final" / "PL2",
                "--out-pl3", tmp / "out" / "PL3", "--work", tmp / "work", *extra, check=check)


def _genome(seq: bytes, circ: bool, q: str) -> dict:
    u = ids.contig_uid(seq).hex()
    return {"contig_uids": [u], "genetic_code": 11, "kind": "plasmid", "circular": [u] if circ else [],
            "assembly_quality": q}


def _hashes(d: Path) -> dict:
    import hashlib

    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(d.iterdir())}


def test_the_split_follows_redistributable_and_keeps_every_record(tmp_path):
    staging = _w2(tmp_path)
    before = _hashes(staging)
    root = _root(tmp_path, {})
    _stage(tmp_path, staging, root)
    pl2, pl3 = tmp_path / "out" / "PL2.new", tmp_path / "out" / "PL3"
    orig = json.loads((staging / "objects.json").read_text())
    o2 = json.loads((pl2 / "objects.json").read_text())
    o3 = json.loads((pl3 / "objects.json").read_text())
    assert [Path(o["path"]).name for o in o2] == ["gut_plasmids_yu2024.fasta", "plasmidscope_RefSeq.fasta"]
    assert [Path(o["path"]).name for o in o3] == ["plasmidscope_IMG-PR.fasta", "plasmidscope_mMGE.fasta"]
    # the paths name the final directory (the orchestrator renames PL2.new into place); the rest is C-40's entry
    assert {Path(o["path"]).parent for o in o2} == {tmp_path / "final" / "PL2"}
    assert {Path(o["path"]).parent for o in o3} == {pl3}
    for o in o2 + o3:
        src = next(x for x in orig if Path(x["path"]).name == Path(o["path"]).name)
        assert {k: v for k, v in o.items() if k != "path"} == {k: v for k, v in src.items() if k != "path"}
        assert list(o) == list(src)
    # every record is in exactly one stage, and the bases add up once the trims are counted
    recs = _tsv(pl2 / "records.tsv") + _tsv(pl3 / "records.tsv")
    want = [(r["source"], r["accession"], int(r["length"])) for r in _tsv(staging / "records.tsv")]
    assert sorted((r["source"], r["accession"]) for r in recs) == sorted((s, a) for s, a, _ in want)
    assert sum(int(r["length"]) + int(r["trimmed_bases"]) for r in recs) == sum(n for *_, n in want)
    assert all(int(r["source_length"]) == int(r["length"]) + int(r["trimmed_bases"]) for r in recs)
    assert _hashes(staging) == before  # card C-40's staging is only read


def test_a_flagged_record_with_an_exact_repeat_is_trimmed_once_and_recorded(tmp_path):
    staging = _w2(tmp_path)
    _stage(tmp_path, staging, _root(tmp_path, {}))
    pl2, pl3 = tmp_path / "out" / "PL2.new", tmp_path / "out" / "PL3"
    # the rewritten objects: every record byte for byte, but the trimmed ones
    assert (pl2 / "gut_plasmids_yu2024.fasta").read_bytes() == (
        b">g1 some description\n" + R + CG + b"\n>g2 some description\n" + CG[::-1] + b"\n"
        b">g3 some description\n" + R + CG[:900] + R + b"\n")
    assert (pl3 / "plasmidscope_IMG-PR.fasta").read_bytes() == b"".join(
        b">" + a + b" some description\n" + s + b"\n"
        for a, s in [(b"i1", R + CI), (b"i1b", R + CI), (b"i2", R + CI), (b"i3", CI[:1500])])
    trims = {s: _tsv(d / "trims.tsv") for s, d in (("PL2", pl2), ("PL3", pl3))}
    assert [(t["accession"], t["repeat_length"], t["bases_removed"], t["length"]) for t in trims["PL2"]] == [
        ("g1", "57", "57", str(len(R + CG)))]
    assert [t["accession"] for t in trims["PL3"]] == ["i1", "i2"]
    assert trims["PL3"][0]["source_contig_uid"] == ids.contig_uid(R + CI + R).hex()
    assert trims["PL3"][0]["contig_uid"] == ids.contig_uid(R + CI).hex()


def test_unchanged_objects_are_hard_links(tmp_path):
    staging = _w2(tmp_path)
    _stage(tmp_path, staging, _root(tmp_path, {}))
    for stage, name in (("PL2.new", "plasmidscope_RefSeq.fasta"), ("PL3", "plasmidscope_mMGE.fasta")):
        assert (tmp_path / "out" / stage / name).stat().st_ino == (staging / name).stat().st_ino
    assert (tmp_path / "out" / "PL3" / "plasmidscope_IMG-PR.fasta").stat().st_ino != \
        (staging / "plasmidscope_IMG-PR.fasta").stat().st_ino


def test_genomes_are_one_per_staged_sequence_with_the_sources_flag(tmp_path):
    staging = _w2(tmp_path)
    _stage(tmp_path, staging, _root(tmp_path, {}))
    pl2, pl3 = tmp_path / "out" / "PL2.new", tmp_path / "out" / "PL3"
    assert json.loads((pl2 / "genomes.json").read_text()) == [
        _genome(R + CG, True, "complete"), _genome(CG[::-1], True, "complete"),
        _genome(R + CG[:900] + R, False, "contig"),  # linear by its flag: its repeat stays
        _genome(CR, True, "complete"), _genome(CR[:1000], False, "contig")]
    # i1 trimmed is i1b's bytes, and i2 is i1's: one genome, from the first record
    assert json.loads((pl3 / "genomes.json").read_text()) == [
        _genome(R + CI, True, "complete"), _genome(CI[:1500], False, "contig"), _genome(CM, True, "complete")]
    r3 = {r["accession"]: r for r in _tsv(pl3 / "records.tsv")}
    assert (r3["i1b"]["first_source"], r3["i1b"]["first_accession"]) == ("plasmidscope_IMG-PR", "i1")
    assert r3["i2"]["first_accession"] == "i1" and r3["i1"]["first_accession"] == ""
    assert r3["i1"]["circular"] == "true" and r3["i3"]["circular"] == "false"
    # the file is written as json.dumps of the list would write it (card C-40's staging)
    text = (pl3 / "genomes.json").read_text()
    assert text == json.dumps(json.loads(text))


def test_a_class_a_trim_cannot_settle_stops_the_staging_until_a_rule_is_given(tmp_path):
    r2 = bytearray(R)
    r2[40] = ord("A") if r2[40] != ord("A") else ord("C")
    p = R + CR[:1500] + bytes(r2)
    staging = _w2(tmp_path, extra_refseq=[("r3", p, "circular", "complete")])
    got = _stage(tmp_path, staging, _root(tmp_path, {}), check=False)
    assert got.returncode == 2
    assert "plasmidscope_RefSeq partial: 1 (e.g. r3)" in got.stderr and "--ambiguous" in got.stderr
    assert not (tmp_path / "out" / "PL2.new").exists() and not (tmp_path / "out" / "PL3").exists()
    for rule, circ in (("linear", False), ("circular", True), ("asis", True)):  # C-40 staged RefSeq's flag
        d = tmp_path / rule
        d.mkdir()
        _stage(d, staging, _root(d, {}), "--ambiguous", rule)
        g = json.loads((d / "out" / "PL2.new" / "genomes.json").read_text())
        assert g[-1] == _genome(p, circ, "complete")  # untrimmed in every case
        assert json.loads((d / "out" / "PL2.new" / "stage_report.json").read_text())["sources"][
            "plasmidscope_RefSeq"]["ambiguous"] == {"partial": 1}


def test_asis_keeps_the_topology_card_c40_gave_an_ambiguous_record(tmp_path):
    r2 = bytearray(R)
    r2[40] = ord("A") if r2[40] != ord("A") else ord("C")
    q = R + CI[:1300] + bytes(r2)
    staging = _w2(tmp_path, extra_img=[("i9", q, "direct terminal repeat", "complete")])
    _stage(tmp_path, staging, _root(tmp_path, {}), "--ambiguous", "asis")
    g = json.loads((tmp_path / "out" / "PL3" / "genomes.json").read_text())
    assert _genome(q, False, "complete") in g  # C-40 staged every IMG-PR record linear


def test_a_record_flagged_circular_and_incomplete_needs_a_rule(tmp_path):
    staging = _w2(tmp_path, extra_refseq=[("r4", R + CR[:1400] + R, "circular", "incomplete")])
    got = _stage(tmp_path, staging, _root(tmp_path, {}), check=False)
    assert got.returncode == 2 and "incomplete" in got.stderr and "r4" in got.stderr
    for rule, seq, circ in (("linear", R + CR[:1400] + R, False), ("circular", R + CR[:1400], True)):
        d = tmp_path / rule
        d.mkdir()
        _stage(d, staging, _root(d, {}), "--fragments", rule)
        g = json.loads((d / "out" / "PL2.new" / "genomes.json").read_text())
        assert g[-1] == _genome(seq, circ, "contig")


def test_the_owner_of_a_contig_is_the_root_stage_that_stores_it(tmp_path):
    staging = _w2(tmp_path)
    root = _root(tmp_path, {"PL1": [ids.contig_uid(CR)], "V1": [ids.contig_uid(R + CI + R)]})
    _stage(tmp_path, staging, root)
    r2 = {r["accession"]: r for r in _tsv(tmp_path / "out" / "PL2.new" / "records.tsv")}
    r3 = {r["accession"]: r for r in _tsv(tmp_path / "out" / "PL3" / "records.tsv")}
    assert r2["r1"]["owner_stage"] == "PL1"
    assert r3["i1"]["owner_stage"] == ""  # V1 holds i1's source bytes, not the trimmed record
    rep = json.loads((tmp_path / "out" / "PL2.new" / "stage_report.json").read_text())
    assert rep["owned_by_stage"] == {"PL1": 1}
    assert rep["new_contigs"] == 4 and rep["contigs"] == 5


def test_overlap_counts_the_staged_sequences_and_the_source_bytes(tmp_path):
    staging = _w2(tmp_path)
    root = _root(tmp_path, {"PL1": [ids.contig_uid(CR)], "V1": [ids.contig_uid(R + CI + R)]})
    _stage(tmp_path, staging, root)
    other = tmp_path / "V2in"
    other.mkdir()
    (other / "genomes.json").write_text(json.dumps([_genome(CM, False, "contig"), _genome(R + CG + R, True, "c")]))
    out = tmp_path / "overlap.json"
    _run("overlap", "--staging", f"PL2={tmp_path / 'out' / 'PL2.new'}", "--staging", f"PL3={tmp_path / 'out' / 'PL3'}",
         "--root", root, "--other", f"V2={other}", "--work", tmp_path / "w2", "--out", out)
    doc = json.loads(out.read_text())
    assert doc["stages"]["PL3"]["as_staged"]["held_by"] == {}
    assert doc["stages"]["PL3"]["source_bytes"]["held_by"] == {"V1": 1}  # the trim stops i1 equalling V1's contig
    assert doc["stages"]["PL3"]["as_staged"]["in_staged"] == {"V2": 1}  # m1
    assert doc["stages"]["PL2"]["as_staged"]["in_staged"] == {"V2": 0}
    assert doc["stages"]["PL2"]["source_bytes"]["in_staged"] == {"V2": 1}  # g1 before its trim
    assert doc["stages"]["PL2"]["as_staged"]["held_by"] == {"PL1": 1}
    assert doc["stages"]["PL2"]["as_staged"]["contigs"] == 5 and doc["stages"]["PL3"]["as_staged"]["contigs"] == 3
    assert doc["shared_PL2_PL3"] == 0


def test_overlap_says_who_holds_each_contig_in_the_build_order(tmp_path):
    # i8 repeats PL2's g2; r1 is PL1's; m1 is listed by V2, built before PL2 and PL3; CR is also in V2 but PL1 holds it
    staging = _w2(tmp_path, extra_img=[("i8", CG[::-1], "linear", "complete")])
    root = _root(tmp_path, {"PL1": [ids.contig_uid(CR)]})
    _stage(tmp_path, staging, root)
    other = tmp_path / "V2in"
    other.mkdir()
    (other / "genomes.json").write_text(json.dumps([_genome(CM, False, "contig"), _genome(CR, True, "complete")]))
    out = tmp_path / "overlap.json"
    _run("overlap", "--staging", f"PL2={tmp_path / 'out' / 'PL2.new'}", "--staging", f"PL3={tmp_path / 'out' / 'PL3'}",
         "--root", root, "--other", f"V2={other}", "--work", tmp_path / "w2", "--out", out)
    doc = json.loads(out.read_text())
    assert doc["build_order"] == ["V2", "PL2", "PL3"]
    assert doc["stages"]["PL2"]["held_in_build_order"] == {"PL1": 1}
    assert doc["stages"]["PL3"]["held_in_build_order"] == {"PL2": 1, "V2": 1}
    assert doc["shared_PL2_PL3"] == 1


# --------------------------------------------------------------------------------------------------------------
# the deferred gut plasmids (DEC-59)

GUT_HEAD = ["contig", "contig_length", "cluster_ID", "cluster_type", "in_dereplicated_plasmids",
            "cluster_dereplication_representative", "circular"]


def test_the_deferred_list_is_the_circular_non_representatives_sorted_with_cluster_and_uid(tmp_path):
    rows = [  # contig, length, cluster, type, representative?, cluster representative, circular
        ("Z_2", 40, "7", "maximal", "False", "Z_9", "True"),
        ("A_1", 40, "3", "backbone", "False", "A_5", "True"),
        ("B_1", 40, "3", "backbone", "True", "B_1", "True"),   # a representative: taken, not deferred
        ("C_1", 40, "4", "fragment", "False", "C_2", "False"),  # linear: not in this list
    ]
    table = tmp_path / "summary.txt"
    table.write_text("# a comment\n# another\n" + "\t".join(GUT_HEAD) + "\n" +
                     "".join("\t".join(map(str, r)) + "\n" for r in rows))
    seqs = {"Z_2": _seq(40, 11), "A_1": _seq(40, 12), "C_1": _seq(40, 13)}
    skipped = tmp_path / "skipped_records.tsv"
    skipped.write_text("source\taccession\tcontig_uid\tlength\treason\n" + "".join(
        f"gut_plasmids_yu2024\t{a}\t{ids.contig_uid(s).hex()}\t{len(s)}\tnot_selected\n" for a, s in seqs.items()) +
        "plasmidscope_GenBank\tX1\t00\t5\tno_row\n")
    out = tmp_path / "deferred.tsv"
    _run("deferred", "--gut-table", table, "--skipped", skipped, "--out", out)
    got = _tsv(out)
    assert [r["contig"] for r in got] == ["A_1", "Z_2"]
    assert got[0] == {"contig": "A_1", "cluster_ID": "3", "cluster_type": "backbone",
                      "cluster_dereplication_representative": "A_5", "contig_length": "40",
                      "contig_uid": ids.contig_uid(seqs["A_1"]).hex()}


def test_the_deferred_list_refuses_a_contig_with_no_record_in_the_fasta(tmp_path):
    table = tmp_path / "summary.txt"
    table.write_text("\t".join(GUT_HEAD) + "\nQ_1\t40\t1\tmaximal\tFalse\tQ_2\tTrue\n")
    skipped = tmp_path / "skipped_records.tsv"
    skipped.write_text("source\taccession\tcontig_uid\tlength\treason\n")
    got = _run("deferred", "--gut-table", table, "--skipped", skipped, "--out", tmp_path / "d.tsv", check=False)
    assert got.returncode == 2 and "Q_1" in got.stderr
    assert not (tmp_path / "d.tsv").exists()
