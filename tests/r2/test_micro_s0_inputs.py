"""Ring R2: persisting genome inputs (card B-04) on the micro-S0 fixture genomes.

Covers C-01 (registry contract) and D-01 (same inputs give the same bytes, in any input order).
"""

from __future__ import annotations

import random
from pathlib import Path

import pyarrow.parquet as pq
import pytest
import yaml

from biobox import catalog, fasta, ids
from biobox.l0seq.inputs import read_genome_input, write_genome_input
from biobox.l1.genecall import GenomeInput
from biobox.registry.schemas import conformance_problems, table

pytestmark = pytest.mark.r2

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
OBJS = yaml.safe_load((FIX / "sources.yaml").read_text())["objects"]


def _fixture_genomes() -> list[GenomeInput]:
    out = []
    for o in OBJS:
        cu = ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)
        circ = (cu,) if o["topology"] == "circular" else ()
        out.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"], circular=circ))
    return out


def _synthetic() -> GenomeInput:
    a, b = ids.contig_uid(b"ACGTACGTAAAA"), ids.contig_uid(b"TTGCATTGCAGG")
    return GenomeInput((b, a), 11, kind="plasmid", circular=(a,), assembly_quality="mag_hq")


def _genomes() -> list[GenomeInput]:
    return [*_fixture_genomes(), _synthetic()]


def _key(g: GenomeInput):
    return (sorted(g.contig_uids), g.genetic_code, g.kind, sorted(g.circular), g.assembly_quality)


def _path(root: Path, stage: str = "S0") -> Path:
    return root / "l0seq" / "genome_input" / f"stage={stage}" / "part-0.parquet"


@pytest.mark.covers("C-01")
def test_round_trip(tmp_path):
    genomes = _genomes()
    stats = write_genome_input(tmp_path, "S0", genomes)
    assert stats == {"genomes": len(genomes), "contigs": sum(len(g.contig_uids) for g in genomes)}
    back = read_genome_input(tmp_path, "S0")
    assert {str(_key(g)) for g in back} == {str(_key(g)) for g in genomes}
    assert len(back) == len(genomes)
    assert [ids.genome_uid(list(g.contig_uids)) for g in back] == sorted(
        ids.genome_uid(list(g.contig_uids)) for g in genomes)
    for g in back:
        assert list(g.contig_uids) == sorted(g.contig_uids)
        assert list(g.circular) == sorted(g.circular)


@pytest.mark.covers("C-01")
def test_file_conforms_to_registry(tmp_path):
    write_genome_input(tmp_path, "S0", _genomes())
    t = pq.read_table(_path(tmp_path))
    assert conformance_problems(t, table("genome_input")) == []
    uids = t.column("contig_uid").to_pylist()
    assert uids == sorted(uids)


@pytest.mark.covers("D-01")
def test_bytes_do_not_depend_on_input_order(tmp_path):
    shuffled = _genomes()
    random.Random(3).shuffle(shuffled)
    shuffled = [GenomeInput(tuple(reversed(g.contig_uids)), g.genetic_code, g.kind, g.circular,
                            g.assembly_quality) for g in shuffled]
    write_genome_input(tmp_path / "a", "S0", _genomes())
    write_genome_input(tmp_path / "b", "S0", shuffled)
    assert _path(tmp_path / "a").read_bytes() == _path(tmp_path / "b").read_bytes()


def _sharing() -> tuple[GenomeInput, GenomeInput, bytes]:
    """`_synthetic` and a second genome that holds its circular contig `a` beside a new one, with the same call
    inputs for `a` (plasmid, code 11, circular) and another assembly quality, which is the genome's own."""
    g = _synthetic()
    a = ids.contig_uid(b"ACGTACGTAAAA")
    other = GenomeInput((ids.contig_uid(b"GGGGCCCCAAAA"), a), 11, kind="plasmid", circular=(a,),
                        assembly_quality="complete")
    return g, other, a


@pytest.mark.covers("C-01")
def test_a_contig_in_two_genomes_of_a_stage_is_listed_in_each(tmp_path):
    """DEC-27 (card C-21a) replaces B-04's refusal of a contig in two genomes of one stage (this test was
    `test_duplicate_contig_is_refused`, whose two genomes shared a contig and disagreed on its call inputs; that
    refusal is kept below, and the agreeing case is now accepted). `genome_input` lists the contig once per genome:
    key (contig_uid, genome_uid, stage), rows sorted by (contig_uid, genome_uid)."""
    g, other, a = _sharing()
    stats = write_genome_input(tmp_path, "S0", [g, other])
    assert stats == {"genomes": 2, "contigs": 3, "rows": 4, "shared_contigs": 1}
    t = pq.read_table(_path(tmp_path))
    assert conformance_problems(t, table("genome_input")) == []
    keys = list(zip(t.column("contig_uid").to_pylist(), t.column("genome_uid").to_pylist()))
    assert keys == sorted(keys) and len(set(keys)) == 4
    guids = {ids.genome_uid(list(x.contig_uids)) for x in (g, other)}
    assert {gu for cu, gu in keys if cu == a} == guids
    back = read_genome_input(tmp_path, "S0")
    assert {str(_key(x)) for x in back} == {str(_key(x)) for x in (g, other)}


@pytest.mark.covers("D-01")
def test_bytes_with_a_shared_contig_do_not_depend_on_input_order(tmp_path):
    g, other, _ = _sharing()
    write_genome_input(tmp_path / "a", "S0", [*_fixture_genomes(), g, other])
    flipped = [GenomeInput(tuple(reversed(x.contig_uids)), x.genetic_code, x.kind, x.circular, x.assembly_quality)
               for x in (other, g, *_fixture_genomes())]
    write_genome_input(tmp_path / "b", "S0", flipped)
    assert _path(tmp_path / "a").read_bytes() == _path(tmp_path / "b").read_bytes()


@pytest.mark.parametrize("change, field", [
    ({"genetic_code": 4}, "genetic_code"),
    ({"kind": "prokaryote"}, "kind"),
    ({"circular": ()}, "is_circular"),
])
def test_two_genomes_that_disagree_on_a_shared_contigs_call_inputs_are_refused_naming_both(tmp_path, change, field):
    """The design note's rule for C-21a (planning/design/C-21_shared_contigs.md §3, §6): a contig has one set of call
    inputs in a stage (kind, genetic code, topology); genomes that disagree on them are refused, naming both."""
    g, other, a = _sharing()
    kw = {"genetic_code": other.genetic_code, "kind": other.kind, "circular": other.circular,
          "assembly_quality": other.assembly_quality, **change}
    bad = GenomeInput(other.contig_uids, **kw)
    with pytest.raises(ValueError, match=field) as e:
        write_genome_input(tmp_path, "S0", [g, bad])
    for x in (g, bad):
        assert ids.genome_uid(list(x.contig_uids)).hex() in str(e.value)
    assert a.hex() in str(e.value)
    assert not _path(tmp_path).exists()


def test_the_same_contig_twice_in_one_genome_is_refused(tmp_path):
    a = ids.contig_uid(b"ACGTACGTAAAA")
    with pytest.raises(ValueError, match="twice"):
        write_genome_input(tmp_path, "S0", [GenomeInput((a, a), 11)])
    assert not _path(tmp_path).exists()


def test_the_same_genome_listed_twice_is_refused(tmp_path):
    """Two entries with the same contigs are one genome (`genome_uid` hashes the contig set): its rows would repeat
    the key."""
    g = _synthetic()
    with pytest.raises(ValueError, match="listed twice"):
        write_genome_input(tmp_path, "S0", [g, GenomeInput(tuple(reversed(g.contig_uids)), 11, kind="plasmid",
                                                           circular=g.circular, assembly_quality="mag_hq")])
    assert not _path(tmp_path).exists()


def test_circular_contig_outside_genome_is_refused(tmp_path):
    g = _synthetic()
    bad = GenomeInput(g.contig_uids, 11, circular=(ids.contig_uid(b"CCCCGGGGTTTT"),))
    with pytest.raises(ValueError, match="circular"):
        write_genome_input(tmp_path, "S0", [bad])
    assert not _path(tmp_path).exists()


def test_second_write_for_a_stage_is_refused(tmp_path):
    write_genome_input(tmp_path, "S0", _genomes())
    before = _path(tmp_path).read_bytes()
    with pytest.raises(FileExistsError):
        write_genome_input(tmp_path, "S0", _genomes())
    assert _path(tmp_path).read_bytes() == before


@pytest.mark.covers("C-01")
def test_iter_streams_the_same_genomes_in_genome_uid_order(tmp_path):
    from biobox.l0seq.inputs import iter_genome_input

    write_genome_input(tmp_path, "S0", _genomes())
    it = iter_genome_input(tmp_path, "S0")
    assert next(it) == read_genome_input(tmp_path, "S0")[0]  # a lazy iterator, not a list
    streamed = list(iter_genome_input(tmp_path, "S0"))
    assert streamed == read_genome_input(tmp_path, "S0")
    assert {str(_key(g)) for g in streamed} == {str(_key(g)) for g in _genomes()}
    guids = [ids.genome_uid(list(g.contig_uids)) for g in streamed]
    assert guids == sorted(guids)


def test_iter_refuses_a_genome_whose_contigs_disagree(tmp_path):
    import pyarrow as pa

    from biobox.l0seq.inputs import iter_genome_input

    write_genome_input(tmp_path / "ok", "S0", [_synthetic()])
    t = pq.read_table(_path(tmp_path / "ok"))
    codes = t.column("genetic_code").to_pylist()
    bad = t.set_column(t.schema.get_field_index("genetic_code"), "genetic_code",
                       pa.array([codes[0], codes[0] + 1], type=pa.uint8()))
    _path(tmp_path / "bad").parent.mkdir(parents=True)
    pq.write_table(bad, _path(tmp_path / "bad"))
    with pytest.raises(ValueError, match="disagree on genetic_code"):
        list(iter_genome_input(tmp_path / "bad", "S0"))


def test_genome_input_lives_below_l1():
    # B-05: l0seq.inputs owns GenomeInput and imports nothing from l1; genecall re-exports it
    import subprocess
    import sys

    code = ("import sys, biobox.l0seq.inputs as i; "
            "bad = sorted(m for m in sys.modules if m.startswith('biobox.l1')); assert not bad, bad; "
            "import biobox.l1.genecall as g; assert g.GenomeInput is i.GenomeInput")
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.mark.covers("C-01")
def test_manifest_accepts_a_root_holding_only_this_file(tmp_path):
    write_genome_input(tmp_path, "S0", _genomes())
    man = catalog.build_manifest(tmp_path, "rS0.test", ["S0"])
    assert [f["table"] for f in man["files"]] == ["genome_input"]
