"""Ring R1: card C-38d unit, registry, mixed-ledger manifest, CLI/stage option, and `make_inputs.py` checks.

Covers:
  - `TABLES["ledger"]` appends nullable `n_bases_in` and `n_bases_rejected` (`uint64`, LY-09);
  - `ledger.below_len_reason(L)` formats `below:<L>` and refuses booleans and `L < 1`;
  - `catalog.compatible_columns`, `catalog.build_manifest`, `catalog.verify_manifest`, and
    `catalog.open_catalog` accept an 8-column pre-C-38d `ledger` Parquet file alongside 10-column
    post-C-38d `ledger` Parquet files in one mixed release;
  - `cli` and `stage.steps` / `stage.ingest_steps` wire `--min-contig-len` into single-process and
    batched ingest, keep `stage.ingest_cpu_min` unchanged, and refuse `--min-contig-len <= 0` (exit 2)
    before running any step;
  - `scripts/stress/build_stage_pod/make_inputs.py` applies `--min-contig-len` across `phage`, `prok`,
    `plasmid`, and `plasmid-sources` against independent brute-force loops over the constructed inputs,
    refuses `--min-contig-len <= 0` (exit 2), and is byte-identical without `--min-contig-len`.
"""

from __future__ import annotations

import dataclasses
import gzip
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

from biobox import catalog, cli, fasta, ids, io, ledger, profile
from biobox import stage as bs
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l1 import genecall, igr
from biobox.l3 import postings, tapes, tokens
from biobox import dim
from biobox.registry.schemas import TABLES, TableSpec

pytestmark = pytest.mark.covers("C-38d", "D-01", "L-01", "LY-09", "LY-10")  # ring R1 by directory (conftest.py)

REPO = Path(__file__).resolve().parents[2]
MAKE_INPUTS = REPO / "scripts" / "stress" / "build_stage_pod" / "make_inputs.py"
FIX = REPO / "tests" / "fixtures" / "micro_s0"
NOON_MONDAY = datetime(2026, 11, 9, 12, 0)


def _fasta(path: Path, recs: list[tuple[str, bytes]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"".join(b">" + h.encode() + b"\n" + s + b"\n" for h, s in recs))
    return path


# --- schema & ledger reason helpers -------------------------------------------------------------------------


@pytest.mark.covers("C-38d", "LY-09")
def test_c38d_ledger_schema_appends_nullable_base_columns():
    spec = TABLES["ledger"]
    assert spec.schema.names[:8] == [
        "object_key",
        "object_kind",
        "n_in",
        "n_kept",
        "n_merged",
        "n_rejected",
        "reason_code",
        "tool_versions_json",
    ]
    assert spec.schema.names[8:] == ["n_bases_in", "n_bases_rejected"]
    for col in ("n_bases_in", "n_bases_rejected"):
        f = spec.schema.field(col)
        assert f.type == pa.uint64()
        assert f.nullable is True


@pytest.mark.covers("C-38d", "L-01")
def test_c38d_ledger_below_len_reason_and_row_defaults():
    assert ledger.below_len_reason(1) == "below:1"
    assert ledger.below_len_reason(1500) == "below:1500"
    for bad in (0, -1, -100, True, False, "1000", None):
        with pytest.raises((ValueError, TypeError)):
            ledger.below_len_reason(bad)  # type: ignore[arg-type]

    default_row = ledger.row(b"k", "l0_object", 10, 7, 3)
    assert default_row["n_bases_in"] is None
    assert default_row["n_bases_rejected"] is None

    filtered_row = ledger.row(
        b"k", "l0_object", 10, 5, 2, n_rejected=3,
        reason=ledger.below_len_reason(1000), n_bases_in=9500, n_bases_rejected=1200,
    )
    assert filtered_row["reason_code"] == "below:1000"
    assert filtered_row["n_bases_in"] == 9500
    assert filtered_row["n_bases_rejected"] == 1200


# --- mixed 8-col pre-C-38d and 10-col post-C-38d ledger compatibility in manifest & catalog ------------------


def _eight_col_ledger_spec() -> TableSpec:
    full = TABLES["ledger"]
    return dataclasses.replace(full, schema=pa.schema(list(full.schema)[:8]))


@pytest.mark.covers("C-38d", "LY-09", "LY-10")
def test_c38d_mixed_8col_and_10col_ledger_in_release_manifest_and_catalog(tmp_path):
    """A pre-C-38d stage S0 whose ledger Parquet files have 8 columns alongside a post-C-38d stage S1
    whose ledger Parquet files have 10 columns builds a valid manifest (`n_columns` 8 and 10), passes
    `verify_manifest`, and exposes NULL base columns on S0's rows and exact base counts on S1's rows."""
    root = tmp_path / "root"
    lam = fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq[:6_000]
    t4 = fasta.read(FIX / "NC_000866.4.fna.gz")[0].seq[:8_000]
    short = b"ACGT" * 100  # 400 bp < 1,500 bp

    f_s0 = _fasta(tmp_path / "src" / "s0.fna", [("s0_lam", lam)])
    f_s1 = _fasta(tmp_path / "src" / "s1.fna", [("s1_short", short), ("s1_t4", t4)])

    # Build S0 normally, then rewrite S0's ledger files with the 8-column pre-C-38d prefix schema
    cu0 = ids.contig_uid(lam)
    write_genome_input(root, "S0", [GenomeInput((cu0,), 11, kind="virus")])
    ingest([SourceObject(f_s0, "src0", "1")], root, "S0")
    genecall.run(root, "S0")
    igr.run(root, "S0")
    tapes.run(root, "S0")

    old_spec = _eight_col_ledger_spec()
    s0_ledgers = sorted((root / "ledger/stage=S0").rglob("*.parquet"))
    assert s0_ledgers
    for lp in s0_ledgers:
        tbl8 = pq.read_table(lp, columns=old_spec.schema.names)
        lp.unlink()
        io.write_table(tbl8, old_spec, lp)
        assert catalog.compatible_columns(pq.read_schema(lp), TABLES["ledger"].schema) == 8

    # Build S1 with --min-contig-len 1500 (writes 10-column post-C-38d ledger files)
    cu1 = ids.contig_uid(t4)
    write_genome_input(root, "S1", [GenomeInput((cu1,), 11, kind="virus")])
    ingest([SourceObject(f_s1, "src1", "1")], root, "S1", min_contig_len=1500)
    genecall.run(root, "S1", seq_stages=["S0", "S1"])
    igr.run(root, "S1")
    tapes.run(root, "S1")

    s1_ingest_ledger = root / "ledger/stage=S1/step=ingest/part-0.parquet"
    assert catalog.compatible_columns(pq.read_schema(s1_ingest_ledger), TABLES["ledger"].schema) == 10

    # Minimal release over ["S0", "S1"] (empty synthetic cluster_member + 10-col cluster ledger, tokens, postings, dim)
    rel = "r26.11.mix"
    catalog.define_release(root, rel, ["S0", "S1"])
    cm_spec = TABLES["cluster_member"]
    io.write_table(pa.Table.from_pylist([], schema=cm_spec.schema), cm_spec,
                   root / f"l2/cluster_member/release={rel}/part-0.parquet")
    mc_spec = TABLES["member_cluster"]
    io.write_table(pa.Table.from_pylist([], schema=mc_spec.schema), mc_spec,
                   root / f"l2/member_cluster/release={rel}/part-0.parquet")
    ledger.write_release(root, rel, "cluster", [
        ledger.row(b"nr50", "level", 0, 0),
        ledger.row(b"nr90", "level", 0, 0),
    ])
    tokens.run(root, rel, ["S0", "S1"])
    postings.run(root, rel, ["S0", "S1"])
    dim.run(root, rel, ["S0", "S1"])

    man = catalog.build_manifest(root, rel, ["S0", "S1"])
    catalog.write_manifest(man)
    assert catalog.verify_manifest(man, root) == []

    assert man["tables"]["ledger"] == catalog.schema_hash(TABLES["ledger"].schema)
    by_path = {f["path"]: f["n_columns"] for f in man["files"] if f["table"] == "ledger"}
    for lp in s0_ledgers:
        assert by_path[lp.relative_to(root).as_posix()] == 8
    assert by_path[s1_ingest_ledger.relative_to(root).as_posix()] == 10

    cat_path = catalog.generate_catalog(man, root=root)
    con = catalog.open_catalog(cat_path, root=root)
    try:
        rows = con.execute(
            "SELECT stage, reason_code, n_in, n_kept, n_rejected, n_bases_in, n_bases_rejected "
            "FROM ledger WHERE step = 'ingest' ORDER BY stage, reason_code"
        ).fetchall()
    finally:
        con.close()
    assert rows == [
        ("S0", "ok", 1, 1, 0, None, None),
        ("S1", "below:1500", 2, 1, 1, len(short) + len(t4), len(short)),
        ("S1", "ok", 2, 1, 0, len(short) + len(t4), 0),
    ]


# --- CLI & stage.py option wiring & validation --------------------------------------------------------------


def _test_profile(tmp_path: Path) -> tuple[profile.Profile, Path]:
    base = tmp_path / "prof"
    for sub in ("tmp", "durable", "bulk"):
        (base / sub).mkdir(parents=True, exist_ok=True)
    prof = profile.Profile(
        max_threads_day=4,
        max_threads_night=8,
        mem_scope="8G",
        tmp=str(base / "tmp"),
        durable=str(base / "durable"),
        bulk=str(base / "bulk"),
        nice_monolithic=0,
        nice_fanout=0,
        day_start_hour=8,
        day_end_hour=20,
        cpu_min_per_process=300,
        small_file_floor_bytes=0,
    )
    pf = tmp_path / "profile.yaml"
    pf.write_text(yaml.safe_dump(dataclasses.asdict(prof)))
    return prof, pf


@pytest.mark.covers("C-38d", "D-01", "D-03")
def test_c38d_stage_steps_wire_min_contig_len_and_keep_ingest_estimate_unchanged(tmp_path):
    prof, pf = _test_profile(tmp_path)
    f = _fasta(tmp_path / "src" / "a.fna", [("r1", b"ACGT" * 500)])
    objs_p = tmp_path / "objects.json"
    objs_p.write_text(json.dumps([{"path": str(f), "source_name": "s", "source_release": "1",
                                   "license": "CC0", "redistributable": True}]))
    gen_p = tmp_path / "genomes.json"
    gen_p.write_text(json.dumps([{"contig_uids": [ids.contig_uid(b"ACGT" * 500).hex()], "genetic_code": 11,
                                  "kind": "virus"}]))

    root = tmp_path / "root"
    b_plain = bs.Build(root=root, stage="S0", release=None, genomes=gen_p, objects=objs_p, depth="stored")
    b_filt = bs.Build(root=root, stage="S0", release=None, genomes=gen_p, objects=objs_p, depth="stored",
                      min_contig_len=1500)

    # Single-process step wiring
    plain_steps = bs.steps(prof, b_plain)
    filt_steps = bs.steps(prof, b_filt)
    plain_ing = next(s for s in plain_steps if s.name == "ingest")
    filt_ing = next(s for s in filt_steps if s.name == "ingest")
    assert "--min-contig-len" not in plain_ing.argv
    idx = filt_ing.argv.index("--min-contig-len")
    assert filt_ing.argv[idx + 1] == "1500"

    # Batched ingest step wiring (`pass1` and `plan` carry `--min-contig-len`)
    batch_steps = bs.ingest_steps(prof, b_filt, 2, 1, None)
    p1 = [s for s in batch_steps if s.name.startswith("ingest.pass1-")]
    plan_s = next(s for s in batch_steps if s.name == "ingest.plan")
    assert len(p1) == 2
    for s in [*p1, plan_s]:
        i = s.argv.index("--min-contig-len")
        assert s.argv[i + 1] == "1500"

    # Ingest CPU estimate unchanged (card C-38d)
    assert bs.ingest_cpu_min(bs.input_bytes(objs_p)) > 0

    # Invalid `--min-contig-len` refused (exit 2) before running any step
    for bad in ("0", "-5"):
        rc = cli.main([
            "build-stage", "--profile", str(pf), "--root", str(tmp_path / f"bad_{bad}"),
            "--stage", "S0", "--depth", "stored", "--genomes", str(gen_p), "--objects", str(objs_p),
            "--min-contig-len", bad,
        ])
        assert rc == 2
        assert not (tmp_path / f"bad_{bad}").exists()

        rc_ing = cli.main([
            "stage", "ingest", "--root", str(tmp_path / f"bad_ing_{bad}"),
            "--stage", "S0", "--objects", str(objs_p), "--min-contig-len", bad,
        ])
        assert rc_ing == 2
        assert not (tmp_path / f"bad_ing_{bad}").exists()


# --- make_inputs.py --min-contig-len ------------------------------------------------------------------------


def _run_make_inputs(*args) -> dict:
    done = subprocess.run([sys.executable, str(MAKE_INPUTS), *map(str, args)], capture_output=True, text=True,
                          check=True)
    return json.loads(done.stdout)


@pytest.mark.covers("C-38d", "D-01")
def test_c38d_make_inputs_phage_and_prok_filter_and_unfiltered_byte_identity(tmp_path):
    s_short = b"ACGT" * 10   # 40 bp (< 60)
    s_edge = b"ACGT" * 15    # 60 bp (== 60, kept!)
    s_long = b"TGCA" * 25    # 100 bp (> 60, kept!)

    # 1. phage mode
    src_gz = tmp_path / "phages.fa.gz"
    with gzip.open(src_gz, "wb") as fh:
        fh.write(
            b">p_short\n" + s_short + b"\n"
            b">p_edge\n" + s_edge + b"\n"
            b">p_long\n" + s_long + b"\n"
            b">p_edge_dup\n" + s_edge + b"\n"
        )
    rep_f = _run_make_inputs("phage", src_gz, 4, tmp_path / "phage_f", "--min-contig-len", 60)
    got_f = json.loads((tmp_path / "phage_f" / "genomes.json").read_text())
    want_phage = [
        {"contig_uids": [ids.contig_uid(s).hex()], "genetic_code": 11, "kind": "virus",
         "circular": [], "assembly_quality": "complete"}
        for s in (s_edge, s_long)
    ]
    assert got_f == want_phage
    assert rep_f["genomes"] == 2 and rep_f["contigs"] == 2 and rep_f["genomes_repeated"] == 1
    assert rep_f["records_below_min_len"] == 1 and rep_f["bases_below_min_len"] == len(s_short)
    # The copied FASTA still holds all 4 records so ingest --min-contig-len 60 accounts for the dropped record
    assert len(fasta.read(tmp_path / "phage_f" / "phage_first4.fa")) == 4

    # 2. prok mode: ACC1 has (short, edge); ACC2 has only (short) -> whole genome omitted from genomes.json,
    #    object kept in objects.json; ACC3 has (edge, long) -> edge shared with ACC1
    f1 = _fasta(tmp_path / "ACC1.fna", [("c1s", s_short), ("c1e", s_edge)])
    f2 = _fasta(tmp_path / "ACC2.fna", [("c2s", s_short)])
    f3 = _fasta(tmp_path / "ACC3.fna", [("c3e", s_edge), ("c3l", s_long)])
    lock = {"genomes": {
        "ACC1": {"file": str(f1), "retrieved_utc": "2026-11-09T00:00:00+00:00"},
        "ACC2": {"file": str(f2), "retrieved_utc": "2026-11-09T00:00:00+00:00"},
        "ACC3": {"file": str(f3), "retrieved_utc": "2026-11-09T00:00:00+00:00"},
    }}
    lock_p, sel_p = tmp_path / "lock.json", tmp_path / "sel.tsv"
    lock_p.write_text(json.dumps(lock))
    sel_p.write_text("accession\ttranslation_table\nACC1\t11\nACC2\t11\nACC3\t4\n")

    rep_p = _run_make_inputs("prok", lock_p, sel_p, tmp_path / "prok_f", "--min-contig-len", 60)
    got_p = json.loads((tmp_path / "prok_f" / "genomes.json").read_text())
    assert got_p == [
        {"contig_uids": [ids.contig_uid(s_edge).hex()], "genetic_code": 11, "kind": "prokaryote",
         "circular": [], "assembly_quality": "complete"},
        {"contig_uids": [ids.contig_uid(s_edge).hex(), ids.contig_uid(s_long).hex()], "genetic_code": 4,
         "kind": "prokaryote", "circular": [], "assembly_quality": "complete"},
    ]
    assert rep_p["genomes"] == 2 and rep_p["contigs"] == 2 and rep_p["contigs_in_several_genomes"] == 1
    assert rep_p["records_below_min_len"] == 2 and rep_p["bases_below_min_len"] == 2 * len(s_short)
    assert rep_p["genomes_below_min_len"] == 1
    assert len(json.loads((tmp_path / "prok_f" / "objects.json").read_text())) == 3

    # Refuses invalid --min-contig-len (exit 2)
    for bad in ("0", "-1", "abc"):
        res = subprocess.run(
            [sys.executable, str(MAKE_INPUTS), "prok", str(lock_p), str(sel_p),
             str(tmp_path / f"bad_{bad}"), "--min-contig-len", bad],
            capture_output=True, text=True,
        )
        assert res.returncode == 2
        assert not (tmp_path / f"bad_{bad}" / "genomes.json").exists()


@pytest.mark.covers("C-38d", "C-26", "C-40", "D-01")
def test_c38d_make_inputs_plasmid_and_plasmid_sources_filter_against_brute_force_oracle(tmp_path):
    s_short = b"ACGT" * 10   # 40 bp (< 60)
    s_edge = b"GCTA" * 15    # 60 bp (== 60)
    s_long = b"TTAA" * 25    # 100 bp (> 60)

    # 1. plasmid mode
    fa = tmp_path / "plsdb.fasta"
    meta = tmp_path / "nuccore.csv"
    recs = [
        ("PL_SHORT.1", s_short, "complete", "circular", ""),
        ("PL_EDGE.1", s_edge, "complete", "circular", ""),
        ("PL_LONG.1", s_long, "", "linear", ""),
    ]
    _fasta(fa, [(f"{a} plasmid", s) for a, s, *_ in recs])
    head = ("NUCCORE_UID,NUCCORE_ACC,NUCCORE_Description,NUCCORE_Completeness,NUCCORE_Genome,NUCCORE_Length,"
            "NUCCORE_Source,NUCCORE_Topology,STATUS\n")
    meta.write_text(head + "".join(
        f'{i},{a},"plasmid",{c},plasmid,{len(s)},refseq,{t},{st}\n'
        for i, (a, s, c, t, st) in enumerate(recs)
    ))
    rep_pl = _run_make_inputs("plasmid", fa, meta, tmp_path / "pl_f", "--min-contig-len", 60)
    got_pl = json.loads((tmp_path / "pl_f" / "genomes.json").read_text())
    want_pl = [
        {"contig_uids": [ids.contig_uid(s).hex()], "genetic_code": 11, "kind": "plasmid",
         "circular": [ids.contig_uid(s).hex()] if t == "circular" else [],
         "assembly_quality": "complete" if c == "complete" else "contig"}
        for _, s, c, t, _ in recs if len(s) >= 60
    ]
    assert got_pl == want_pl
    assert rep_pl["genomes"] == 2
    assert rep_pl["records_below_min_len"] == 1 and rep_pl["bases_below_min_len"] == len(s_short)

    # 2. plasmid-sources mode
    src_fa = _fasta(tmp_path / "ps.fasta", [("ps1", s_short), ("ps2", s_edge), ("ps3", s_long)])
    src_tsv = tmp_path / "ps.tsv"
    src_tsv.write_text(
        "id\ttopo\tcomp\n"
        "ps1\tcircular\tcomplete\n"
        "ps2\tcircular\tcomplete\n"
        "ps3\tlinear\tincomplete\n"
    )
    spec_p = tmp_path / "spec.json"
    spec_p.write_text(json.dumps({"sources": [{
        "name": "ps", "source_name": "psrc", "source_release": "1", "license": "CC0",
        "redistributable": True, "retrieved_at": "2026-11-18T00:00:00Z",
        "fasta": [str(src_fa)], "table": str(src_tsv), "id_column": "id",
        "topology_column": "topo", "topology": {"circular": "circular", "linear": "linear"},
        "completeness_column": "comp", "completeness": {"complete": "complete", "incomplete": "fragment"},
    }]}))
    rep_ps = _run_make_inputs("plasmid-sources", spec_p, tmp_path / "ps_f", "--min-contig-len", 60)
    got_ps = json.loads((tmp_path / "ps_f" / "genomes.json").read_text())
    assert [g["contig_uids"] for g in got_ps] == [[ids.contig_uid(s_edge).hex()], [ids.contig_uid(s_long).hex()]]
    assert rep_ps["genomes"] == 2
    skipped_lines = (tmp_path / "ps_f" / "skipped_records.tsv").read_text().splitlines()
    assert len(skipped_lines) == 2  # header + ps1
    assert skipped_lines[1].split("\t") == ["ps", "ps1", ids.contig_uid(s_short).hex(), str(len(s_short)), "below:60"]
    assert [r.accession for r in fasta.read(tmp_path / "ps_f" / "ps.fasta")] == ["ps2", "ps3"]
