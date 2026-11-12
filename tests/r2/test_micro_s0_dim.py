"""Ring R2: `contig_dim` on the micro-S0 fixture (card B-12).

Covers C-01 (registry contract) and D-01 (same inputs give the same bytes).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

from biobox import catalog, dim, fasta, ids
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.registry.schemas import conformance_problems, table

pytestmark = pytest.mark.r2

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
REL = "rS0.dim"


def _genomes(quality: str = "complete") -> list[GenomeInput]:
    out = []
    for o in CONTRACT["objects"]:
        cu = ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)
        out.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"],
                               circular=(cu,) if o["topology"] == "circular" else (), assembly_quality=quality))
    return out


def _store(root: Path, genomes: list[GenomeInput] | None = None) -> Path:
    ingest([SourceObject(FIX / o["file"], "s", "1") for o in CONTRACT["objects"]], root, "S0", shard_bases=200_000)
    write_genome_input(root, "S0", _genomes() if genomes is None else genomes)
    return root


def _out(root: Path, rel: str = REL) -> Path:
    return root / "dim" / "contig_dim" / f"release={rel}" / "part-0.parquet"


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = _store(tmp_path_factory.mktemp("dim"))
    stats = dim.run(root, REL, ["S0"])
    return root, stats


def _expected():
    rows = {}
    for g, o in zip(_genomes(), CONTRACT["objects"]):
        (cu,) = g.contig_uids
        rows[cu] = (ids.genome_uid([cu]), len(fasta.read(FIX / o["file"])[0].seq), o["topology"] == "circular")
    return rows


@pytest.mark.covers("C-01")
def test_one_row_per_default_contig(built):
    root, stats = built
    t = pq.read_table(_out(root)).to_pylist()
    want = _expected()
    assert [r["contig_uid"] for r in t] == sorted(want)
    for r in t:
        guid, n, circ = want[r["contig_uid"]]
        assert (r["genome_uid"], r["contig_len"], r["is_circular"]) == (guid, n, circ)
        assert r["assembly_quality"] == "complete"
        assert (r["tier"], r["tax_order"], r["is_species_rep"]) == (None, None, False)
    assert sum(r["is_circular"] for r in t) == 1 and len(t) == 3
    assert stats == {"contigs": 3, "genomes": 3, "circular": 1}


@pytest.mark.covers("C-01")
def test_file_conforms_and_enters_manifest(built):
    root, _ = built
    assert conformance_problems(pq.read_table(_out(root)), table("contig_dim")) == []
    man = catalog.build_manifest(root, REL, ["S0"])
    assert "contig_dim" in {f["table"] for f in man["files"]}


@pytest.mark.covers("C-01")
def test_contig_dim_now_reads_only_the_catalog_release(built, tmp_path):
    root, _ = built
    want = pq.read_table(_out(root)).to_pylist()
    con = catalog.open_catalog(catalog.generate_catalog(catalog.build_manifest(root, REL, ["S0"]), tmp_path / "a.duckdb"))
    cols = ", ".join(table("contig_dim").schema.names)
    got = con.execute(f"SELECT {cols} FROM contig_dim_now() ORDER BY contig_uid").fetchall()
    assert got == [tuple(r.values()) for r in want]
    con.close()
    # card B-33: REL's contig_dim file is REL's only, so another release's manifest does not list it
    mine = catalog.build_manifest(root, REL, ["S0"])
    theirs = catalog.build_manifest(root, "rOther", ["S0"])
    assert not [f["path"] for f in theirs["files"] if f["table"] == "contig_dim"]
    other = catalog.open_catalog(catalog.generate_catalog(theirs, tmp_path / "b.duckdb"))
    assert other.execute("SELECT count(*) FROM contig_dim_now()").fetchone()[0] == 0
    assert other.execute("SELECT count(*) FROM contig_dim").fetchone()[0] == 0
    other.close()
    # and contig_dim_now() filters by release when the view does hold another release's rows (a hand-made list)
    mixed = {**theirs, "files": theirs["files"] + [f for f in mine["files"] if f["table"] == "contig_dim"]}
    other = catalog.open_catalog(catalog.generate_catalog(mixed, tmp_path / "m.duckdb"))
    assert other.execute("SELECT count(*) FROM contig_dim_now()").fetchone()[0] == 0
    assert other.execute("SELECT count(*) FROM contig_dim").fetchone()[0] == 3
    other.close()


def test_contig_dim_now_is_empty_without_any_dim_file(tmp_path):
    root = _store(tmp_path / "s")
    con = catalog.open_catalog(catalog.generate_catalog(catalog.build_manifest(root, REL, ["S0"]), tmp_path / "c.duckdb"))
    assert con.execute("SELECT count(*) FROM contig_dim_now()").fetchone()[0] == 0
    con.close()


def test_missing_genome_input_is_refused(tmp_path):
    root = _store(tmp_path, _genomes()[:1])
    with pytest.raises(ValueError, match=r"2 default-variant contigs have no genome_input row; first: [0-9a-f]{32}"):
        dim.run(root, REL, ["S0"])
    assert not _out(root).exists()


def test_stray_genome_input_contig_is_refused(tmp_path):
    stray = GenomeInput((ids.contig_uid(b"ACGTACGTAAAA"),), 11, kind="plasmid")
    root = _store(tmp_path, [*_genomes(), stray])
    with pytest.raises(ValueError, match="genome_input rows name a contig that contig_index does not hold"):
        dim.run(root, REL, ["S0"])
    assert not _out(root).exists()


def test_conflicting_duplicate_across_stages_is_refused_and_identical_collapses(tmp_path):
    root = _store(tmp_path)
    shutil.copytree(root / "l0seq/contig_index/stage=S0", root / "l0seq/contig_index/stage=S1")
    write_genome_input(root, "S1", _genomes())  # identical rows: one row per contig
    assert dim.run(root, REL, ["S0", "S1"], tmp=tmp_path) == {"contigs": 3, "genomes": 3, "circular": 1}
    assert pq.read_table(_out(root)).num_rows == 3
    shutil.rmtree(root / "l0seq/genome_input/stage=S1")
    write_genome_input(root, "S1", _genomes("mag_hq"))
    with pytest.raises(ValueError, match="conflicting genome_input rows"):
        dim.run(root, "rConflict", ["S0", "S1"])
    assert not _out(root, "rConflict").exists()


def test_second_build_for_a_release_refuses_to_overwrite(built):
    root, _ = built
    before = _out(root).read_bytes()
    with pytest.raises(FileExistsError):
        dim.run(root, REL, ["S0"])
    assert _out(root).read_bytes() == before


def _shared_store(root: Path, tmp_path: Path) -> tuple[Path, bytes, bytes]:
    """S0 as `_store`; S1 ingests T4 again beside a new contig N and lists both in one genome (card C-20). S0 stored
    T4 first, so S0 owns it (DEC-26); S1's contig_index holds N only."""
    _store(root)
    t4 = fasta.read(FIX / "NC_000866.4.fna.gz")[0].seq
    new = t4[:5000][::-1]
    src = tmp_path / "s1.fna"
    src.write_bytes(b">t4 again\n" + t4 + b"\n>new\n" + new + b"\n")
    ingest([SourceObject(src, "s", "1")], root, "S1", shard_bases=200_000)  # the store's parameters, as `_store`
    cu_t4, cu_new = ids.contig_uid(t4), ids.contig_uid(new)
    write_genome_input(root, "S1", [GenomeInput((cu_t4, cu_new), 11, kind="prokaryote")])
    return root, cu_t4, cu_new


@pytest.mark.covers("C-01")
def test_a_contig_shared_across_stages_keeps_its_owning_stages_genome_and_the_later_row_is_listed_skipped(tmp_path):
    root, cu_t4, cu_new = _shared_store(tmp_path / "s", tmp_path)
    res = dim.run(root, "rShared", ["S0", "S1"], tmp=tmp_path)
    later = ids.genome_uid([cu_t4, cu_new])
    assert res["contigs"] == 4 and res["genomes"] == 4
    assert res["skipped"] == [{"contig_uid": cu_t4.hex(), "stage": "S1", "genome_uid": later.hex(),
                               "owner_stage": "S0", "kept_genome_uid": ids.genome_uid([cu_t4]).hex()}]
    rows = {r["contig_uid"]: r for r in pq.read_table(_out(root, "rShared")).to_pylist()}
    assert rows[cu_t4]["genome_uid"] == ids.genome_uid([cu_t4])  # S0's genome, not S1's
    assert rows[cu_new]["genome_uid"] == later
    assert len(rows) == 4


def test_a_release_with_the_later_stage_but_not_the_owning_stage_is_refused_naming_it(tmp_path):
    root, cu_t4, _ = _shared_store(tmp_path / "s", tmp_path)
    with pytest.raises(ValueError, match=rf"genome_input rows name a contig that contig_index does not hold.*"
                                         rf"{cu_t4.hex()}.*stage S0 .*DEC-26"):
        dim.run(root, "rLate", ["S1"], tmp=tmp_path)
    assert not _out(root, "rLate").exists()


@pytest.mark.covers("C-01")
def test_a_contig_in_two_genomes_of_one_stage_keeps_the_least_genome_and_lists_the_other_skipped(tmp_path):
    """DEC-27 (card C-21a) migrates C-20's `test_a_contig_in_two_genomes_of_one_stage_is_still_refused`: a contig may
    be in several genomes of one stage. Until card C-21b adds `n_genomes` and `genome_dim`, `contig_dim` keeps one row
    per contig, its genome the one with the least `genome_uid` (the genome that writes the contig's gene calls), and
    the other genomes' rows are listed under `skipped`, as C-20 lists a later stage's. The refusal of a real conflict
    within one stage is the next test."""
    root = tmp_path / "s"
    ingest([SourceObject(FIX / o["file"], "s", "1") for o in CONTRACT["objects"]], root, "S0", shard_bases=200_000)
    single = _genomes()
    lam, t4 = single[0].contig_uids[0], single[1].contig_uids[0]
    both = GenomeInput((lam, t4), 11, kind="virus", assembly_quality="mag_hq")  # lambda and T4 in one more genome
    write_genome_input(root, "S0", [*single, both])
    res = dim.run(root, "rTwo", ["S0"], tmp=tmp_path)
    g_both = ids.genome_uid([lam, t4])
    kept = {cu: min(ids.genome_uid([cu]), g_both) for cu in (lam, t4)}
    other = {cu: max(ids.genome_uid([cu]), g_both) for cu in (lam, t4)}
    rows = {r["contig_uid"]: r for r in pq.read_table(_out(root, "rTwo")).to_pylist()}
    assert len(rows) == 3
    for cu in (lam, t4):
        assert rows[cu]["genome_uid"] == kept[cu]
        assert rows[cu]["assembly_quality"] == ("mag_hq" if kept[cu] == g_both else "complete")
    m = single[2].contig_uids[0]
    assert rows[m]["genome_uid"] == ids.genome_uid([m])
    assert res["skipped"] == sorted(
        ({"contig_uid": cu.hex(), "stage": "S0", "genome_uid": other[cu].hex(), "owner_stage": "S0",
          "kept_genome_uid": kept[cu].hex()} for cu in (lam, t4)), key=lambda r: r["contig_uid"])
    assert res["contigs"] == 3 and res["genomes"] == len({r["genome_uid"] for r in rows.values()})


def test_two_genomes_of_one_stage_that_disagree_on_a_contigs_topology_are_refused(tmp_path):
    """A real conflict within one stage is still refused: two rows that disagree on the contig's own fact
    (`is_circular`). `write_genome_input` refuses to write them (C-21a), so the file is written past it."""
    root = tmp_path / "s"
    ingest([SourceObject(FIX / o["file"], "s", "1") for o in CONTRACT["objects"]], root, "S0", shard_bases=200_000)
    spec = table("genome_input")
    rows = [{"contig_uid": g.contig_uids[0], "genome_uid": ids.genome_uid(list(g.contig_uids)),
             "is_circular": bool(g.circular), "kind": g.kind, "genetic_code": g.genetic_code,
             "assembly_quality": g.assembly_quality} for g in _genomes()]
    rows.append({**rows[0], "genome_uid": b"\x01" * 16, "is_circular": not rows[0]["is_circular"]})
    d = root / spec.location / "stage=S0"
    d.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist(sorted(rows, key=lambda r: (r["contig_uid"], r["genome_uid"])),
                                        schema=spec.schema), d / "part-0.parquet")
    with pytest.raises(ValueError, match="conflicting genome_input rows.*is_circular"):
        dim.run(root, "rTwo", ["S0"])
    assert not _out(root, "rTwo").exists()


@pytest.mark.covers("D-01")
def test_output_bytes_are_identical_across_fresh_builds(tmp_path):
    a, b = _store(tmp_path / "a"), _store(tmp_path / "b")
    dim.run(a, REL, ["S0"])
    dim.run(b, REL, ["S0"])
    assert _out(a).read_bytes() == _out(b).read_bytes()
