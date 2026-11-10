"""Ring R2: IGRs on micro-S0 (R-05; C-03 and C-05 for the igr table; LY-08)."""

from __future__ import annotations

import shutil
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

from biobox import catalog, contracts, fasta, ids, io
from biobox.l0seq import SeqStore, SourceObject, ingest
from biobox.l0seq.inputs import write_genome_input
from biobox.l1 import igr
from biobox.l1.igr import igrs_of_contig
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, run as call_genes
from biobox.l3 import tapes
from biobox.registry.schemas import conformance_problems, table

pytestmark = pytest.mark.r2

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())


def _build(root: Path) -> set[bytes]:
    ingest([SourceObject(FIX / o["file"], "s", "1") for o in CONTRACT["objects"]], root, "S0", shard_bases=200_000)
    genomes, circular = [], set()
    for o in CONTRACT["objects"]:
        cu = ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)
        if o["topology"] == "circular":
            circular.add(cu)
        genomes.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"], circular=(cu,) if cu in circular else ()))
    call_genes(root, "S0", genomes)
    igr.run(root, "S0", circular)
    return circular


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = tmp_path_factory.mktemp("igr")
    circular = _build(root)
    rows = pq.read_table(root / f"l0seq/igr/igr_version={igr.igr_version()}/stage=S0/part-0.parquet").to_pylist()
    return root, circular, rows


def _slots(root):
    occ = pq.read_table(root / f"l1/protein_occurrence/gene_call_version={GENE_CALL_VERSION}/stage=S0/part-0.parquet")
    return {(r["contig_uid"], r["slot_ordinal"]): r for r in occ.to_pylist()}


def test_r05_igrs_lie_between_their_flanking_slots(built):
    root, circular, rows = built
    slots = _slots(root)
    store = SeqStore(root)
    assert rows
    for r in rows:
        n = store.length_of(store._default[r["contig_uid"]])
        if r["left_slot"] is not None:
            assert slots[(r["contig_uid"], r["left_slot"])]["end"] % n == r["start"] % n
        if r["right_slot"] is not None:
            assert slots[(r["contig_uid"], r["right_slot"])]["start"] == r["end"] % n or \
                slots[(r["contig_uid"], r["right_slot"])]["start"] == r["end"]
    # no slot overlaps any IGR interior
    by_contig = {}
    for (cu, _), s in slots.items():
        by_contig.setdefault(cu, []).append(s)
    for r in rows:
        n = store.length_of(store._default[r["contig_uid"]])
        span = {p % n for p in range(r["start"], r["end"])}
        for s in by_contig.get(r["contig_uid"], []):
            assert not span & {p % n for p in range(s["start"], s["end"])}


def test_r05_igr_slices_equal_the_source(built, tmp_path):
    root, circular, rows = built
    store = SeqStore(root)
    con = catalog.open_catalog(catalog.generate_catalog(catalog.build_manifest(root, "x"), tmp_path / "c.duckdb"))
    src = {ids.contig_uid(fasta.read(FIX / o["file"])[0].seq): fasta.read(FIX / o["file"])[0].seq for o in CONTRACT["objects"]}
    for r in rows[::7]:
        seq = src[r["contig_uid"]]
        want = (seq + seq)[r["start"]:r["end"]]
        assert store.get(r["contig_uid"], r["start"], r["end"], 1, circular=r["contig_uid"] in circular) == want
        got = con.execute("SELECT seq_get_circular(?, ?, ?, 1)", [r["contig_uid"], r["start"], r["end"]]).fetchall()[0][0]
        assert got.encode() == want
    con.close()


def test_igr_table_passes_contracts(built, tmp_path):
    root, _, _ = built
    man = catalog.build_manifest(root, "x")
    con = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / "c.duckdb"))
    assert contracts.check_foreign_keys(con) == []
    assert contracts.check_coordinates(con) == []
    con.close()
    assert contracts.check_sort_order(man) == []


def test_ly08_igr_rewrite_is_byte_identical(built, tmp_path):
    root, circular, _ = built
    out = root / f"l0seq/igr/igr_version={igr.igr_version()}/stage=S0/part-0.parquet"
    before = out.read_bytes()
    fresh = tmp_path / "root"
    shutil.copytree(root, fresh)  # the writer refuses to overwrite, so rebuild in a copy without the table
    shutil.rmtree(fresh / "l0seq/igr")
    igr.run(fresh, "S0", circular)
    assert (fresh / out.relative_to(root)).read_bytes() == before


# --- card B-06: flanks reference features ------------------------------------------------------

def _igr_path(root, producers=(GENE_CALL_VERSION,)):
    return root / f"l0seq/igr/igr_version={igr.igr_version(list(producers))}/stage=S0/part-0.parquet"


@pytest.mark.covers("R-05")
def test_b20_igr_version_names_the_tape_version_and_the_definition_is_written(built):
    root, _, _ = built
    tv = tapes.tape_version([GENE_CALL_VERSION])
    assert igr.igr_version([GENE_CALL_VERSION]) == f"{tv}-igr1"
    assert igr.igr_version(GENE_CALL_VERSION) == igr.igr_version() == f"{tv}-igr1"  # a string is one producer
    assert igr.igr_version(["b", "a"]) == igr.igr_version(["a", "b"]) == f"{tapes.tape_version(['a', 'b'])}-igr1"
    assert _igr_path(root).exists()
    d = pq.read_table(tapes.def_path(root, tv), columns=["producer"], partitioning=None)
    assert d.column("producer").to_pylist() == [GENE_CALL_VERSION]


def _contig_lengths(root):
    idx = pq.read_table(root / "l0seq/contig_index/stage=S0/part-0.parquet").to_pylist()
    return {r["contig_uid"]: r["contig_len"] for r in idx if r["is_default_variant"]}


def _reference(root, circular):
    """(contig_uid, start, end, left_slot, right_slot) the old way: protein_occurrence + igrs_of_contig."""
    by = {}
    for r in _slots(root).values():
        by.setdefault(r["contig_uid"], []).append((r["slot_ordinal"], r["start"], r["end"]))
    out = []
    for cu, n in _contig_lengths(root).items():
        out += [(cu, s, e, left, right) for s, e, left, right in igrs_of_contig(by.get(cu, []), n, cu in circular)]
    return sorted(out, key=lambda t: (t[0], t[1]))


@pytest.mark.covers("R-05")
def test_b06_igr_conforms_and_matches_the_old_computation(built):
    root, circular, rows = built
    spec = table("igr")
    t = pq.read_table(_igr_path(root))
    assert t.schema.equals(spec.schema)
    assert t.schema.names[-2:] == ["left_feature_uid", "right_feature_uid"]
    assert conformance_problems(t, spec) == []
    got = [(r["contig_uid"], r["start"], r["end"], r["left_slot"], r["right_slot"]) for r in rows]
    assert got == _reference(root, circular)
    assert {cu for cu, *_ in got} == set(_contig_lengths(root))


@pytest.mark.covers("R-05")
def test_b06_flank_feature_uids_exist_and_touch_the_igr(built):
    """Asserted relation, with n the contig length: the left feature's `end % n` equals the IGR's `start % n`
    and the right feature's `start` equals the IGR's `end % n` (the origin-crossing run has end > n)."""
    root, circular, rows = built
    feats = {r["feature_uid"]: r for r in pq.read_table(
        root / f"l1/feature/producer={GENE_CALL_VERSION}/stage=S0/part-0.parquet").to_pylist()}
    occ = {r["feature_uid"]: r["slot_ordinal"] for r in _slots(root).values()}
    n_of = _contig_lengths(root)
    wrapped = 0
    for r in rows:
        n = n_of[r["contig_uid"]]
        wrapped += r["end"] > n
        for side, slot in (("left", "left_slot"), ("right", "right_slot")):
            fu = r[f"{side}_feature_uid"]
            assert (fu is None) == (r[slot] is None)  # nulls exactly where the legacy slot is null
            if fu is None:
                continue
            f = feats[fu]
            assert f["contig_uid"] == r["contig_uid"]
            assert occ[fu] == r[slot]
            if side == "left":
                assert f["end"] % n == r["start"] % n
            else:
                assert f["start"] == r["end"] % n
    assert wrapped >= 1  # the circular contig has an origin-crossing IGR
    # a gene-free or linear-end flank is null: linear contigs have one IGR with a null left at 0
    assert any(r["left_feature_uid"] is None for r in rows)


@pytest.mark.covers("R-05")
def test_b06_second_producer_splits_one_igr(built, tmp_path):
    root, circular, rows = built
    fresh = tmp_path / "root"
    shutil.copytree(root, fresh)  # B-20: the second producer set writes beside the first, no removal
    n_of = _contig_lengths(fresh)
    mg = next(cu for cu in circular)  # the only circular contig is M. genitalium (NC_000908.2)
    target = next(r for r in rows if r["contig_uid"] == mg and r["end"] - r["start"] >= 100 and r["end"] <= n_of[mg]
                  and r["left_feature_uid"] is not None and r["right_feature_uid"] is not None)
    s, e = target["start"] + 30, target["start"] + 70
    fu = ids.feature_uid(mg, s, e, 1, ids.KIND["ncrna"])
    spec = table("feature")
    batch = pa.RecordBatch.from_pylist(
        [{"feature_uid": fu, "contig_uid": mg, "start": s, "end": e, "strand": 1, "kind": ids.KIND["ncrna"]}],
        schema=spec.schema)
    io.write_table([batch], spec, fresh / "l1/feature/producer=test-ncrna/stage=S0/part-0.parquet")
    both = (GENE_CALL_VERSION, "test-ncrna")
    igr.run(fresh, "S0", circular, producers=list(both))
    assert _igr_path(fresh).read_bytes() == _igr_path(root).read_bytes()  # the first file is untouched
    new = pq.read_table(_igr_path(fresh, both)).to_pylist()
    old_keys = {(r["contig_uid"], r["start"]): r for r in rows}
    new_keys = {(r["contig_uid"], r["start"]): r for r in new}
    assert len(new) == len(rows) + 1
    gone = set(old_keys) - set(new_keys)
    added = set(new_keys) - set(old_keys)
    assert gone == set()  # the first half keeps the replaced IGR's key (same start)
    assert added == {(mg, e)}
    for k in set(old_keys) - {(mg, target["start"])}:
        assert new_keys[k] == old_keys[k]
    first, second = new_keys[(mg, target["start"])], new_keys[(mg, e)]
    assert (first["end"], second["end"]) == (s, target["end"])
    assert first["left_feature_uid"] == target["left_feature_uid"] and first["left_slot"] == target["left_slot"]
    assert first["right_feature_uid"] == fu and first["right_slot"] is None
    assert second["left_feature_uid"] == fu and second["left_slot"] is None
    assert second["right_feature_uid"] == target["right_feature_uid"] and second["right_slot"] == target["right_slot"]


def test_b06_default_topology_comes_from_genome_input(built, tmp_path):
    root, circular, _ = built
    fresh = tmp_path / "root"
    shutil.copytree(root, fresh)
    shutil.rmtree(fresh / "l0seq/igr")
    write_genome_input(fresh, "S0", [GenomeInput((cu,), circular=(cu,) if cu in circular else ()) for cu in _contig_lengths(fresh)])
    igr.run(fresh, "S0")  # circular=None
    assert _igr_path(fresh).read_bytes() == _igr_path(root).read_bytes()
