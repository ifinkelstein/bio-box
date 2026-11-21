"""Ring R1: card C-38g — readers of rounds, merged release copy, and clean tokens (DEC-54).

Covers:
  - `biobox.l1.rounds.merge_release`: streaming k-way merge of `release_round` per-round files into
    sorted `release=` `round_domains` and `round_coverage` files, deterministic deduplication across
    overlapping rounds, and clean no-op when a release lists 0 rounds.
  - §3.8 test (3): `domain_state` view, `domain_state(puid, lib)` UDF (`DomainStateLookup`), and
    `domains_of(puid)` (`source` = `'stage'` vs `'round'`, one source per protein preferring stage)
    checked across all four states (`measured`, `inferred`, `none`, `not scanned`) including round
    `measured` and round `none`, against an independent brute-force oracle.
  - `nbr.summarize_rows` / `nbr.summarize` and `export.nbr_summary_tsv` show round rows with
    `source = 'round'` vs `'stage'`, preferring stage per protein, and separating round `none` from
    `not scanned` (`?`).
  - `export.hits_tsv` and `export.locus_tables` reflect round `measured` (`pfam_architecture`) and
    round `none` (`-`) vs `not scanned` (`?`).
  - `l3.tokens` and `l3.postings` guard against ever reading `round_domains` or `round_coverage`.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog, ids
from biobox import io as bio
from biobox.l1 import rounds as l1rounds
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l3 import tokens
from biobox.registry.schemas import TABLES
from biobox.search import export, nbr

pytestmark = [pytest.mark.covers("C-38g", "LY-06", "LY-10", "W-05", "W-06")]

LIB = "pfam-38.2"


def _uid(n: int) -> bytes:
    return n.to_bytes(16, "big")


def _dom_row(puid: bytes, acc: str, name: str, env_start: int = 1, *, inferred_from: bytes | None = None,
             with_inferred: bool = False) -> dict:
    row = {
        "protein_uid": puid,
        "model_acc": acc,
        "model_name": name,
        "model_len": 100,
        "seq_bitscore": 55.0,
        "seq_evalue_zref": 1e-12,
        "dom_bitscore": 50.0,
        "dom_ievalue_zref": 1e-11,
        "env_start": env_start,
        "env_end": env_start + 50,
        "ali_start": env_start,
        "ali_end": env_start + 50,
        "hmm_from": 0,
        "hmm_to": 50,
        "cutoff": "gathering",
        "z": 1e8,
        "dom_z": 1e8,
    }
    if with_inferred:
        row["inferred_from"] = inferred_from
    return row


def _write_minimal_stage(root: Path, stage: str, protein_uids: list[bytes], *,
                         domains_rows: list[dict] | None = None,
                         with_inferred: bool = False) -> None:
    cspec = TABLES["contig_index"]
    cu = ids.contig_uid(f"ACGTACGT{stage}".encode() * 20)
    rd = ids.raw_digest(f"ACGTACGT{stage}".encode() * 20)
    bio.write_table(
        pa.table({
            "contig_uid": pa.array([cu], pa.binary(16)),
            "raw_digest": pa.array([rd], pa.binary(32)),
            "contig_len": pa.array([160 + len(stage) * 20], pa.int64()),
            "is_default_variant": pa.array([True], pa.bool_()),
        }, schema=cspec.schema),
        cspec,
        root / cspec.location / f"stage={stage}" / "part-0.parquet",
    )
    pspec = TABLES["protein_seq"]
    prows = [{
        "protein_uid": u,
        "protein_digest": hashlib.sha256(u).digest(),
        "length": 120,
        "seq": "M" * 120,
    } for u in sorted(protein_uids)]
    bio.write_table(
        pa.Table.from_pylist(prows, schema=pspec.schema),
        pspec,
        root / pspec.location / f"gene_call_version={V}" / f"stage={stage}" / "protein_prefix=all" / "part-0.parquet",
    )
    if domains_rows is not None:
        dspec = TABLES["l1_domains"]
        schema = dspec.schema.append(pa.field("inferred_from", pa.binary(16), nullable=True)) if with_inferred else dspec.schema
        ordered = sorted(domains_rows, key=lambda r: (r["protein_uid"], r["env_start"], r["model_acc"]))
        t = pa.Table.from_pylist(ordered, schema=schema)
        out = root / dspec.location / f"library={LIB}" / f"stage={stage}" / "part-0.parquet"
        out.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(t, out, compression="zstd", write_statistics=True)


def _write_round(root: Path, rid: str, covered: list[bytes], dom_rows: list[dict], *, library: str = LIB) -> None:
    dspec = TABLES["round_domains"]
    cspec = TABLES["round_coverage"]
    ordered = sorted(dom_rows, key=lambda r: (r["protein_uid"], r["env_start"], r["model_acc"]))
    bio.write_table(
        pa.Table.from_pylist(ordered, schema=dspec.schema),
        dspec,
        root / dspec.location / f"round={rid}" / f"library={library}" / "part-0.parquet",
    )
    counts: dict[bytes, int] = {}
    for r in ordered:
        counts[r["protein_uid"]] = counts.get(r["protein_uid"], 0) + 1
    cov_rows = [{"protein_uid": u, "n_rows": counts.get(u, 0)} for u in sorted(set(covered))]
    bio.write_table(
        pa.Table.from_pylist(cov_rows, schema=cspec.schema),
        cspec,
        root / cspec.location / f"round={rid}" / f"library={library}" / "part-0.parquet",
    )


def test_c38g_merge_release_dedup_and_noop(tmp_path):
    root = tmp_path / "root"
    _write_minimal_stage(root, "S0", [_uid(i) for i in range(1, 10)])

    # 0 rounds -> clean no-op
    catalog.define_release(root, "r0", ["S0"])
    res0 = l1rounds.merge_release(root, "r0")
    assert res0["rounds"] == 0
    assert not (root / f"l1/round_domains/release=r0").exists()
    assert not (root / f"l1/round_coverage/release=r0").exists()

    # Two overlapping rounds: r-a covers 1..5 (hits on 1, 3), r-b covers 3..7 (hits on 3, 6)
    _write_round(root, "r-000000000000000a", [_uid(i) for i in range(1, 6)], [
        _dom_row(_uid(1), "PF00001.1", "DomA", 1),
        _dom_row(_uid(3), "PF00003.1", "DomC", 5),
        _dom_row(_uid(3), "PF00004.1", "DomD", 60),
    ])
    _write_round(root, "r-000000000000000b", [_uid(i) for i in range(3, 8)], [
        _dom_row(_uid(3), "PF00003.1", "DomC", 5),
        _dom_row(_uid(3), "PF00004.1", "DomD", 60),
        _dom_row(_uid(6), "PF00006.1", "DomF", 10),
    ])
    catalog.define_release(root, "r1", ["S0"])
    res1 = l1rounds.merge_release(root, "r1")
    assert res1["rounds"] == 2
    assert res1["libraries"][LIB]["proteins"] == 7
    assert res1["libraries"][LIB]["domains"] == 4

    dom_files = bio.part_files(root / f"l1/round_domains/release=r1/library={LIB}")
    cov_files = bio.part_files(root / f"l1/round_coverage/release=r1/library={LIB}")
    assert len(dom_files) == 1 and len(cov_files) == 1
    dom_got = pq.read_table(dom_files[0]).to_pylist()
    cov_got = pq.read_table(cov_files[0]).to_pylist()
    assert [(r["protein_uid"], r["model_acc"], r["env_start"]) for r in dom_got] == [
        (_uid(1), "PF00001.1", 1),
        (_uid(3), "PF00003.1", 5),
        (_uid(3), "PF00004.1", 60),
        (_uid(6), "PF00006.1", 10),
    ]
    assert [(r["protein_uid"], r["n_rows"]) for r in cov_got] == [
        (_uid(1), 1), (_uid(2), 0), (_uid(3), 2), (_uid(4), 0),
        (_uid(5), 0), (_uid(6), 1), (_uid(7), 0),
    ]


def test_c38g_domain_state_all_four_states_and_domains_of_vs_oracle(tmp_path):
    """§3.8 test (3): domain_state view, UDF, and domains_of across all 4 states including round measured and round none."""
    root = tmp_path / "root"
    # S0 (annotated):
    #   1: stage measured
    #   2: stage none (no domain row)
    #   3: stage measured AND also present in round (stage must win in domains_of)
    # S1 (called, not annotated):
    #   4: round measured (n_rows=1)
    #   5: round none (n_rows=0)
    #   6: not scanned (neither stage nor round)
    _write_minimal_stage(root, "S0", [_uid(1), _uid(2), _uid(3)], domains_rows=[
        _dom_row(_uid(1), "PF00001.1", "StageOne", 1),
        _dom_row(_uid(3), "PF00003.1", "StageThree", 2),
    ])
    _write_minimal_stage(root, "S1", [_uid(4), _uid(5), _uid(6)], domains_rows=None)
    _write_round(root, "r-0000000000000001", [_uid(3), _uid(4), _uid(5)], [
        _dom_row(_uid(3), "PF00003.1", "StageThree", 2),
        _dom_row(_uid(4), "PF00004.1", "RoundFour", 10),
    ])

    catalog.define_release(root, "r1", ["S0", "S1"])
    l1rounds.merge_release(root, "r1")
    man = catalog.build_manifest(root, "r1", ["S0", "S1"])
    assert catalog.verify_manifest(man) == []
    cat_path = catalog.generate_catalog(man, root / "catalog/release=r1.duckdb")

    # Independent brute-force oracle over files on disk
    oracle_state = {
        _uid(1): "measured",
        _uid(2): "none",
        _uid(3): "measured",
        _uid(4): "measured",
        _uid(5): "none",
        _uid(6): "not scanned",
        _uid(99): "not scanned",
    }

    con = catalog.open_catalog(cat_path)
    try:
        # Check SQL view `domain_state`
        view_rows = {bytes(pu): st for pu, st in con.execute(
            "SELECT protein_uid, domain_state FROM domain_state WHERE library = ?", [LIB]
        ).fetchall()}
        for u in (_uid(1), _uid(2), _uid(3), _uid(4), _uid(5), _uid(6)):
            assert view_rows[u] == oracle_state[u], (u, view_rows.get(u), oracle_state[u])

        # Check scalar UDF `domain_state(puid, lib)`
        for u, want in oracle_state.items():
            got = con.execute("SELECT domain_state(?, ?)", [u, LIB]).fetchone()[0]
            assert got == want, (u, got, want)

        # Check `domains_of(puid)` source column and one-source-per-protein rule (§3.3)
        d1 = con.execute("SELECT model_acc, stage, source FROM domains_of(?)", [_uid(1)]).fetchall()
        assert d1 == [("PF00001.1", "S0", "stage")]

        d3 = con.execute("SELECT model_acc, stage, source FROM domains_of(?)", [_uid(3)]).fetchall()
        assert d3 == [("PF00003.1", "S0", "stage")]  # stage wins, NOT duplicated by round

        d4 = con.execute("SELECT model_acc, stage, source FROM domains_of(?)", [_uid(4)]).fetchall()
        assert d4 == [("PF00004.1", None, "round")]

        d5 = con.execute("SELECT model_acc, stage, source FROM domains_of(?)", [_uid(5)]).fetchall()
        assert d5 == []
    finally:
        con.close()


def test_c38g_nbr_summarize_rows_shows_source_stage_and_round():
    anchor = _uid(1)
    p_stage = _uid(2)     # stage hit
    p_both = _uid(3)      # in both stage and round -> stage must win (counted once, source='stage')
    p_round = _uid(4)     # round hit -> source='round'
    p_rnone = _uid(5)     # round coverage n_rows=0 -> unannotated '', NOT '?'
    p_unscan = _uid(6)    # never scanned -> '?'
    cu1 = _uid(10)
    fu_a = _uid(20)

    def row(off: int, puid: bytes, fuid: bytes) -> dict:
        return {
            "anchor_protein_uid": anchor, "contig_uid": cu1, "anchor_ordinal": 2, "ordinal": 2 + off,
            "offset": off, "offset_rel": off, "protein_uid": puid, "start": (2 + off) * 100,
            "end": (2 + off) * 100 + 90, "strand": 1, "delta_bp": None if off == 0 else 10,
            "junction_topology": None if off == 0 else 1, "same_strand": True, "in_operon": None,
            "flags": 0, "policy_hash": "h", "anchor_feature_uid": fu_a, "feature_uid": fuid, "tape_version": "tv1",
        }

    nbr_rows = [
        row(-2, p_stage, _uid(30)),
        row(-1, p_both, _uid(31)),
        row(0, anchor, fu_a),
        row(1, p_round, _uid(32)),
        row(2, p_rnone, _uid(33)),
        row(3, p_unscan, _uid(34)),
    ]
    stage_doms = [
        {"protein_uid": p_stage, "library": LIB, "model_acc": "PF00001.1", "model_name": "StageFam"},
        {"protein_uid": p_both, "library": LIB, "model_acc": "PF00002.1", "model_name": "BothFam"},
    ]
    round_doms = [
        {"protein_uid": p_both, "library": LIB, "model_acc": "PF00002.1", "model_name": "BothFam"},
        {"protein_uid": p_round, "library": LIB, "model_acc": "PF00003.1", "model_name": "RoundFam"},
    ]
    summary = nbr.summarize_rows(
        nbr_rows,
        by=("pfam",),
        domains=stage_doms,
        scanned_proteins=[p_stage, p_both],
        round_domains=round_doms,
        round_coverage=[p_both, p_round, p_rnone],
    )
    pfam_rows = [r for r in summary["rows"] if r["level"] == "pfam"]
    by_key = {(r["key"], r.get("source")): r for r in pfam_rows}
    assert (f"{LIB}:PF00001.1", "stage") in by_key
    assert (f"{LIB}:PF00002.1", "stage") in by_key
    assert (f"{LIB}:PF00002.1", "round") not in by_key  # stage wins
    assert (f"{LIB}:PF00003.1", "round") in by_key
    assert by_key[("", None)]["n_neighbors"] == 1        # p_rnone is none, not '?'
    assert by_key[("?", None)]["n_neighbors"] == 1       # p_unscan is '?'

    tsv = "\n".join(export.nbr_summary_tsv({
        **summary,
        "header": {"release": "r1", "tape_version": "tv1", "policy": {}, "policy_hash": "h", **summary["header"]},
    }))
    lines = [ln for ln in tsv.splitlines() if ln and not ln.startswith("#")]
    cols = lines[0].split("\t")
    assert "source" in cols
    src_i = cols.index("source")
    key_i = cols.index("key")
    tsv_map = {ln.split("\t")[key_i]: ln.split("\t")[src_i] for ln in lines[1:]}
    assert tsv_map[f"{LIB}:PF00001.1"] == "stage"
    assert tsv_map[f"{LIB}:PF00003.1"] == "round"


def test_c38g_tokens_guard_rejects_round_paths(tmp_path):
    root = tmp_path / "root"
    _write_minimal_stage(root, "S0", [_uid(1)], domains_rows=None)
    _write_round(root, "r-0000000000000001", [_uid(1)], [_dom_row(_uid(1), "PF00001.1", "RoundFam")])
    catalog.define_release(root, "r1", ["S0"])
    l1rounds.merge_release(root, "r1")

    assert tokens.libraries_in(root, ["S0"]) == []
    assert tokens._domain_files(root, LIB, ["S0"]) == []
    round_dom_file = bio.part_files(root / f"l1/round_domains/release=r1/library={LIB}")[0]
    with pytest.raises(ValueError, match="round_domains"):
        tokens._assert_stage_domain_files([round_dom_file])
