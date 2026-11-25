"""Ring R1: card B-70 — normalize MMseqs2 `0M` empty member alignments in `cluster_member_aln`
and validate canonical empty rows in `contracts.check_coordinates`.
"""

from __future__ import annotations

from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import contracts, ids, layout
from biobox.l2 import cluster
from biobox.registry.schemas import conformance_problems, table


class _StubMMseqs:
    """Stub for `cluster.MMseqs` that writes prebuilt `clu*.tsv` and `aln_*.tsv` files."""

    def __init__(self, workdir: Path, clu_rows: dict[str, list[tuple[str, str]]],
                 aln_rows: dict[str, list[tuple]]):
        self.workdir = workdir
        self.clu_rows = clu_rows
        self.aln_rows = aln_rows
        self.log = open(workdir / "mmseqs.log", "w")

    def __call__(self, *args: str, threaded: bool = True) -> None:
        cmd = args[0]
        if cmd == "createtsv":
            # createtsv db db <clu> <out_tsv>
            clu_name, out_rel = args[3], args[4]
            lvl = "nr90" if clu_name == "clu90" else "nr50"
            (self.workdir / out_rel).write_text(
                "".join(f"{rep}\t{mem}\n" for rep, mem in self.clu_rows[lvl])
            )
        elif cmd == "convertalis":
            # convertalis db db aln_<lvl> <out_tsv> --format-output ...
            aln_name, out_rel = args[3], args[4]
            lvl = aln_name.removeprefix("aln_")
            (self.workdir / out_rel).write_text(
                "".join("\t".join(str(x) for x in row) + "\n" for row in self.aln_rows[lvl])
            )
        else:
            raise AssertionError(f"unexpected mmseqs subcommand in stub: {args}")


def _write_prot_parquet(path: Path, uids: list[bytes]) -> None:
    hex_uids = sorted(u.hex() for u in uids)
    pq.write_table(
        pa.table({"uid": pa.array(hex_uids, type=pa.string())},
                 schema=pa.schema([pa.field("uid", pa.string(), nullable=False)])),
        path,
    )


def _coords_con(prots: list[tuple[bytes, str]], aln_table: pa.Table) -> duckdb.DuckDBPyConnection:
    """In-memory DuckDB connection with the tables `contracts.check_coordinates` queries."""
    con = duckdb.connect()
    empty_tables = ("protein_occurrence", "contig_index", "igr", "l1_domains", "contig_tapes", "token_postings")
    for name in empty_tables:
        spec = table(name)
        schema = spec.schema
        for pk in spec.partition_keys:
            if pk not in schema.names:
                schema = schema.append(pa.field(pk, pa.string()))
        con.register(name, pa.Table.from_pylist([], schema=schema))
    p_rows = [
        {"protein_uid": u, "seq": s, "length": len(s), "md5": b"\x00" * 16}
        for u, s in prots
    ]
    con.register("protein_seq", pa.Table.from_pylist(p_rows, schema=table("protein_seq").schema))
    con.register("cluster_member_aln", aln_table)
    return con


@pytest.mark.covers("C-05", "C-11")
def test_b70_empty_0m_alignments_normalized_and_conformant(tmp_path, monkeypatch):
    """Reproduces `r26.11.3` `mmseqs convertalis` `0M` rows (`qstart=0, qend=1, tstart=0, tend=0,
    alnlen=0, cigar="0M"`, with `fident=0.5` and `fident=0.0`) alongside normal `alnlen > 0` rows."""
    monkeypatch.setattr(cluster.provenance, "tool_versions", lambda *names: {n: "test" for n in names})

    rep_seq = "M" + "A" * 103  # 104 aa representative
    norm_seq = "M" + "A" * 97  # 98 aa normal member
    z1_seq = "M" + "C" * 57    # 58 aa member (like r26.11.3 protein 80107049)
    z2_seq = "M" + "D" * 37    # 38 aa member (like r26.11.3 protein 38179080)
    rep_u, norm_u, z1_u, z2_u = (ids.protein_uid(s) for s in (rep_seq, norm_seq, z1_seq, z2_seq))
    rep_h, norm_h, z1_h, z2_h = (u.hex() for u in (rep_u, norm_u, z1_u, z2_u))

    w = tmp_path / "work"
    out = tmp_path / "out"
    w.mkdir()
    out.mkdir()
    _write_prot_parquet(w / "prot.parquet", [rep_u, norm_u, z1_u, z2_u])

    clu_rows = {
        lvl: [(rep_h, rep_h), (rep_h, norm_h), (rep_h, z1_h), (rep_h, z2_h)]
        for lvl in ("nr90", "nr50")
    }
    # convertalis columns: query (rep), target (member), qstart, qend, tstart, tend, fident, alnlen, cigar
    aln_rows = {
        lvl: [
            (rep_h, rep_h, 1, 104, 1, 104, 1.0, 104, "104M"),
            (rep_h, norm_h, 2, 101, 1, 98, 0.92, 100, "50M2I48M"),
            (rep_h, z1_h, 0, 1, 0, 0, 0.5, 0, "0M"),
            (rep_h, z2_h, 0, 1, 0, 0, 0.0, 0, "0M"),
        ]
        for lvl in ("nr90", "nr50")
    }

    mm = _StubMMseqs(w, clu_rows, aln_rows)
    res = cluster._write_tables(mm, w, w, out, "rB70.test", 4, "256MB", layout.MAX_FILE_BYTES, member_align=False)
    assert res["proteins"] == 4
    assert res["aligned"] == {"nr90": 4, "nr50": 4}

    aln = pq.read_table(out / "l2/cluster_member_aln/release=rB70.test/part-0.parquet")
    assert conformance_problems(aln, table("cluster_member_aln")) == []

    by_key = {(r["level"], r["member_uid"]): r for r in aln.to_pylist()}
    for lvl in ("nr90", "nr50"):
        for z_u in (z1_u, z2_u):
            row = by_key[(lvl, z_u)]
            assert row["rep_uid"] == rep_u
            assert row["qstart"] == 0
            assert row["qend"] == 0
            assert row["tstart"] == 0
            assert row["tend"] == 0
            assert row["pident"] == 0.0
            assert row["alnlen"] == 0
            assert row["cigar"] == b""

        norm_row = by_key[(lvl, norm_u)]
        assert (norm_row["qstart"], norm_row["qend"], norm_row["tstart"], norm_row["tend"]) == (0, 98, 1, 101)
        assert norm_row["pident"] == pytest.approx(0.92, abs=1e-6)
        assert norm_row["alnlen"] == 100
        assert norm_row["cigar"] == b"50M2I48M"

        self_row = by_key[(lvl, rep_u)]
        assert (self_row["qstart"], self_row["qend"], self_row["tstart"], self_row["tend"]) == (0, 104, 0, 104)
        assert self_row["pident"] == pytest.approx(1.0)
        assert self_row["alnlen"] == 104
        assert self_row["cigar"] == b"104M"

    con = _coords_con([(rep_u, rep_seq), (norm_u, norm_seq), (z1_u, z1_seq), (z2_u, z2_seq)], aln)
    try:
        assert contracts.check_coordinates(con) == []
    finally:
        con.close()


@pytest.mark.covers("C-11")
def test_b70_nonempty_alignment_output_and_plan_format_unchanged(tmp_path, monkeypatch):
    """`PLAN_FORMAT` is 2 since card B-71 (DEC-77, with format-1 read compatibility),
    and non-empty `alnlen > 0` rows produce identical bytes."""
    assert cluster.PLAN_FORMAT == 2
    monkeypatch.setattr(cluster.provenance, "tool_versions", lambda *names: {n: "test" for n in names})

    rep_seq = "M" + "K" * 79
    mem_seq = "M" + "K" * 74
    rep_u, mem_u = ids.protein_uid(rep_seq), ids.protein_uid(mem_seq)
    rep_h, mem_h = rep_u.hex(), mem_u.hex()

    w, out = tmp_path / "work", tmp_path / "out"
    w.mkdir()
    out.mkdir()
    _write_prot_parquet(w / "prot.parquet", [rep_u, mem_u])

    clu_rows = {lvl: [(rep_h, rep_h), (rep_h, mem_h)] for lvl in ("nr90", "nr50")}
    aln_rows = {
        lvl: [
            (rep_h, rep_h, 1, 80, 1, 80, 1.0, 80, "80M"),
            (rep_h, mem_h, 3, 77, 1, 75, 0.96, 75, "75M"),
        ]
        for lvl in ("nr90", "nr50")
    }
    mm = _StubMMseqs(w, clu_rows, aln_rows)
    cluster._write_tables(mm, w, w, out, "rB70.norm", 2, "256MB", layout.MAX_FILE_BYTES, member_align=False)
    aln = pq.read_table(out / "l2/cluster_member_aln/release=rB70.norm/part-0.parquet")
    assert conformance_problems(aln, table("cluster_member_aln")) == []
    by_key = {(r["level"], r["member_uid"]): r for r in aln.to_pylist()}
    for lvl in ("nr90", "nr50"):
        r = by_key[(lvl, mem_u)]
        assert (r["qstart"], r["qend"], r["tstart"], r["tend"], r["alnlen"], r["cigar"]) == (0, 75, 2, 77, 75, b"75M")
        assert r["pident"] == pytest.approx(0.96, abs=1e-6)


@pytest.mark.covers("C-05")
@pytest.mark.parametrize(
    "bad_fields",
    [
        {"alnlen": 0, "qstart": 0, "qend": 0, "tstart": 0, "tend": 0, "pident": 0.0, "cigar": b"0M"},
        {"alnlen": 0, "qstart": 0, "qend": 1, "tstart": 0, "tend": 0, "pident": 0.0, "cigar": b""},
        {"alnlen": 0, "qstart": 1, "qend": 0, "tstart": 0, "tend": 0, "pident": 0.0, "cigar": b""},
        {"alnlen": 0, "qstart": 0, "qend": 0, "tstart": 1, "tend": 0, "pident": 0.0, "cigar": b""},
        {"alnlen": 0, "qstart": 0, "qend": 0, "tstart": 0, "tend": 1, "pident": 0.0, "cigar": b""},
        {"alnlen": 0, "qstart": 0, "qend": 0, "tstart": 0, "tend": 0, "pident": 0.5, "cigar": b""},
        {"alnlen": 10, "qstart": 0, "qend": 10, "tstart": 0, "tend": 10, "pident": 0.9, "cigar": b""},
        {"alnlen": 10, "qstart": 0, "qend": 60, "tstart": 0, "tend": 10, "pident": 0.9, "cigar": b"10M"},
        {"alnlen": 10, "qstart": 0, "qend": 10, "tstart": 0, "tend": 120, "pident": 0.9, "cigar": b"10M"},
    ],
)
def test_b70_check_coordinates_rejects_noncanonical_empty_and_invalid_nonempty(bad_fields):
    rep_seq = "M" + "A" * 99   # 100 aa
    mem_seq = "M" + "C" * 49   # 50 aa
    rep_u, mem_u = ids.protein_uid(rep_seq), ids.protein_uid(mem_seq)
    spec = table("cluster_member_aln")

    # Canonical empty alignment is accepted
    ok_tbl = pa.Table.from_pylist([{
        "member_uid": mem_u, "level": "nr50", "rep_uid": rep_u,
        "qstart": 0, "qend": 0, "tstart": 0, "tend": 0, "pident": 0.0, "alnlen": 0, "cigar": b"",
    }], schema=spec.schema)
    con_ok = _coords_con([(rep_u, rep_seq), (mem_u, mem_seq)], ok_tbl)
    try:
        assert contracts.check_coordinates(con_ok) == []
    finally:
        con_ok.close()

    # Non-canonical empty or invalid non-empty row is rejected
    bad_tbl = pa.Table.from_pylist([{
        "member_uid": mem_u, "level": "nr50", "rep_uid": rep_u, **bad_fields,
    }], schema=spec.schema)
    con_bad = _coords_con([(rep_u, rep_seq), (mem_u, mem_seq)], bad_tbl)
    try:
        assert contracts.check_coordinates(con_bad) == ["cluster_member_aln: 1 rows outside sequence bounds"]
    finally:
        con_bad.close()
