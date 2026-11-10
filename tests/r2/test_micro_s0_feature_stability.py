"""Ring R2: feature addition stability on micro-S0 (MR-20, card B-08; AGENTS §3.3 "an ordinal is never a key").

Build the store once with the gene-call producer alone (ingest, genome inputs, gene calls, IGRs, tapes), then
add one 40 nt kind-2 feature inside an IGR of the M. genitalium contig as a second producer `test-ncrna` and
rebuild IGRs and tapes for both producers. Everything keyed by `feature_uid` must be unchanged; only the
ordinals inside the new tape version move (by one, after the new feature).

Card B-20: `igr.igr_version` names the producer set (`<tape_version>-igr1`), so both producer sets are built in
ONE store. Every file that existed before the second build is byte-identical afterwards (the first IGR file
included); the second set's IGR file sits beside it under its own `igr_version`.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

from biobox import fasta, ids, io
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import write_genome_input
from biobox.l1 import igr
from biobox.l1.genecall import GENE_CALL_VERSION as V, GenomeInput, run as call_genes
from biobox.l3 import tapes
from biobox.registry.schemas import table

pytestmark = pytest.mark.r2

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
NEW = "test-ncrna"
LISTS = ("feature_uids", "protein_uids", "strands", "start_bps", "end_bps", "feature_kinds")


def _rows(root: Path, glob: str) -> list[dict]:
    return [r for f in sorted(root.glob(glob)) for r in pq.read_table(f, partitioning=None).to_pylist()]


def _snapshot(root: Path, producers: list[str]) -> dict:
    return {
        "feature": _rows(root, "l1/feature/producer=*/stage=S0/part-*.parquet"),
        "feature_protein": _rows(root, f"l1/feature_protein/gene_call_version={V}/stage=S0/part-*.parquet"),
        "protein_occurrence": _rows(root, f"l1/protein_occurrence/gene_call_version={V}/stage=S0/part-*.parquet"),
        "igr": _rows(root, f"l0seq/igr/igr_version={igr.igr_version(producers)}/stage=S0/part-*.parquet"),
    }


def _digests(root: Path) -> dict[str, str]:
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file()}


def _tapes(root: Path, producers: list[str]) -> dict[bytes, dict]:
    tv = tapes.tape_version(producers)
    return {r["contig_uid"]: r for r in _rows(root, f"l3/tapes/tape_version={tv}/stage=S0/part-*.parquet")}


def _by_feature_uid(tape: dict) -> dict[bytes, tuple]:
    """slot rows of one tape keyed by feature_uid (never by position)."""
    cols = zip(tape["feature_uids"], tape["start_bps"], tape["end_bps"], tape["strands"], tape["feature_kinds"],
               tape["protein_uids"])
    return {fu: rest for fu, *rest in map(tuple, cols)}


@pytest.fixture(scope="module")
def stability(tmp_path_factory):
    base = tmp_path_factory.mktemp("stability")
    root = base / "first"
    ingest([SourceObject(FIX / o["file"], CONTRACT["source_name"], CONTRACT["source_release"], CONTRACT["license"], True)
            for o in CONTRACT["objects"]], root, "S0", shard_bases=200_000)
    genomes, circ = [], set()
    for o in CONTRACT["objects"]:
        cu = ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)
        c = (cu,) if o["topology"] == "circular" else ()
        circ |= set(c)
        genomes.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"], circular=c))
    write_genome_input(root, "S0", genomes)
    call_genes(root, "S0", genomes)
    igr.run(root, "S0", circ)
    tapes.run(root, "S0", circ)
    before = _snapshot(root, [V])
    before_tapes = _tapes(root, [V])
    before_files = _digests(root)

    mg = next(iter(circ))  # the only circular contig is M. genitalium (NC_000908.2)
    n_mg = before_tapes[mg]["contig_len"]
    target = next(r for r in sorted(before["igr"], key=lambda r: r["start"])
                  if r["contig_uid"] == mg and r["end"] - r["start"] >= 200 and r["end"] <= n_mg)
    s, e, kind = target["start"] + 80, target["start"] + 120, ids.KIND["ncrna"]
    fu = ids.feature_uid(mg, s, e, 1, kind)

    second = root  # one store: the second producer set is built beside the first
    spec = table("feature")
    io.write_table([pa.RecordBatch.from_pylist(
        [{"feature_uid": fu, "contig_uid": mg, "start": s, "end": e, "strand": 1, "kind": kind}], schema=spec.schema)],
        spec, second / f"l1/feature/producer={NEW}/stage=S0/part-0.parquet")
    producers = [V, NEW]
    igr.run(second, "S0", circ, producers=producers)
    tapes.run(second, "S0", circ, producers=producers)
    return {"mg": mg, "target": target, "s": s, "e": e, "kind": kind, "fu": fu, "circ": circ,
            "before": before, "before_tapes": before_tapes, "before_files": before_files,
            "after": _snapshot(second, producers), "after_tapes": _tapes(second, producers), "after_files": _digests(second),
            "second": second, "first": root}


@pytest.mark.covers("MR-20")
def test_mr20_feature_addition_changes_nothing_it_does_not_touch(stability):
    d = stability
    mg, fu, s, e = d["mg"], d["fu"], d["s"], d["e"]
    assert tapes.tape_version([V]) != tapes.tape_version([V, NEW])

    # 1. every file that existed before the second build is byte-identical afterwards, with no exception
    assert len([k for k in d["before_files"] if k.startswith("l0seq/igr/")]) == 1
    for k, h in d["before_files"].items():
        assert d["after_files"][k] == h, k
    added_files = set(d["after_files"]) - set(d["before_files"])
    assert added_files and all(k.startswith(("l3/", "l1/feature/producer=test-ncrna/", "ledger/", "l0seq/igr/"))
                               for k in added_files), added_files
    assert not any("tape_version=" + tapes.tape_version([V]) in k for k in added_files)

    # 2. feature: every earlier row unchanged, exactly one new row across producers
    old_f, new_f = d["before"]["feature"], d["after"]["feature"]
    assert len(new_f) == len(old_f) + 1
    assert all(r in new_f for r in old_f)
    [added] = [r for r in new_f if r not in old_f]
    assert (added["feature_uid"], added["contig_uid"], added["start"], added["end"], added["strand"], added["kind"]) \
        == (fu, mg, s, e, 1, ids.KIND["ncrna"])

    # 3. feature_protein and protein_occurrence: identical row sets
    for name in ("feature_protein", "protein_occurrence"):
        assert d["after"][name] == d["before"][name], name
        assert d["before"][name]

    # 4. igr: one IGR (the one containing the new feature) is replaced by two; the others are identical
    target = d["target"]
    old_i, new_i = d["before"]["igr"], d["after"]["igr"]
    assert target["start"] < s and e < target["end"]
    assert len(new_i) == len(old_i) + 1
    others_old = [r for r in old_i if r != target]
    split = [r for r in new_i if r["contig_uid"] == mg and target["start"] <= r["start"] < target["end"]]
    assert [(r["start"], r["end"]) for r in split] == [(target["start"], s), (e, target["end"])]
    others_new = [r for r in new_i if r not in split]
    assert len(others_old) == len(others_new) == len(old_i) - 1
    # `left_slot` and `right_slot` are legacy ordinal-like values (B-17 removes them). Today they are the gene-call
    # slot_ordinals from protein_occurrence, which the second producer does not touch, so unrelated IGRs keep them;
    # that is all that is asserted about them. Everything else, including both flank feature uids, is compared.
    assert others_new == others_old
    stable = [c for c in old_i[0] if c not in ("left_slot", "right_slot")]
    assert [{c: r[c] for c in stable} for r in others_new] == [{c: r[c] for c in stable} for r in others_old]
    first, second = split
    assert first["left_feature_uid"] == target["left_feature_uid"] and first["right_feature_uid"] == fu
    assert second["left_feature_uid"] == fu and second["right_feature_uid"] == target["right_feature_uid"]
    assert first["left_slot"] == target["left_slot"] and second["right_slot"] == target["right_slot"]
    assert first["right_slot"] is None and second["left_slot"] is None  # the ncRNA is not a gene


@pytest.mark.covers("MR-20")
def test_b20_the_two_igr_files_sit_under_different_igr_versions(stability):
    d = stability
    v1, v2 = igr.igr_version([V]), igr.igr_version([V, NEW])
    assert v1 != v2 and v1 == f"{tapes.tape_version([V])}-igr1" and v2 == f"{tapes.tape_version([V, NEW])}-igr1"
    p1 = f"l0seq/igr/igr_version={v1}/stage=S0/part-0.parquet"
    p2 = f"l0seq/igr/igr_version={v2}/stage=S0/part-0.parquet"
    assert p1 in d["before_files"] and p2 not in d["before_files"]
    assert p2 in d["after_files"]
    assert d["after_files"][p1] == d["before_files"][p1]  # the first file's sha256 is unchanged
    assert {k for k in d["after_files"] if k.startswith("l0seq/igr/")} == {p1, p2}


@pytest.mark.covers("MR-20")
def test_mr20_ordinals_move_inside_the_new_tape_but_feature_uids_do_not(stability):
    d = stability
    mg, fu = d["mg"], d["fu"]
    old, new = d["before_tapes"], d["after_tapes"]
    assert set(new) == set(old) and len(old) == len(CONTRACT["objects"])
    # the other contigs: equal in every list column
    for cu in set(old) - {mg}:
        for col in table("contig_tapes").list_columns:
            assert new[cu][col] == old[cu][col], col
        assert new[cu]["n_slots"] == old[cu]["n_slots"]
    o, n = old[mg], new[mg]
    assert n["n_slots"] == o["n_slots"] + 1
    i = n["feature_uids"].index(fu)
    assert (n["start_bps"][i], n["end_bps"][i], n["strands"][i], n["feature_kinds"][i]) == (d["s"], d["e"], 1, d["kind"])
    assert n["protein_uids"][i] is None
    for col in LISTS:
        assert n[col][:i] + n[col][i + 1:] == o[col], col
    # the point of the test: ordinals after the new feature rise by exactly one, feature uids stay
    after = range(i + 1, n["n_slots"])
    assert len(after) > 0 and i > 0
    old_ord = {u: k for k, u in enumerate(o["feature_uids"])}
    new_ord = {u: k for k, u in enumerate(n["feature_uids"])}
    assert set(old_ord) <= set(new_ord) and len(new_ord) == len(old_ord) + 1
    for k in range(i):
        assert new_ord[o["feature_uids"][k]] == k
    for k in after:
        u = n["feature_uids"][k]
        assert u in old_ord  # the same feature_uid is still on the tape ...
        assert new_ord[u] == old_ord[u] + 1  # ... at an ordinal larger by exactly one
    # and every slot row keyed by feature_uid is unchanged
    old_rows, new_rows = _by_feature_uid(o), _by_feature_uid(n)
    assert all(new_rows[u] == r for u, r in old_rows.items())
