"""Ring R1: card R-17, `stage_check.py` on a stage that holds a non-default variant of a contig another stage stores.

On 2026-11-26 `V2`'s check failed with "contigs both stored by the stage and held by another: 43": for 43 contigs that
`S1` stores, `V2`'s `contig_index` holds one row each with `is_default_variant` false and a `raw_digest` no earlier
stage has. That is the data model: `contig_index` maps `contig_uid` to the raw variants of a contig, "a later stage
adds only new (contig_uid, raw_digest) pairs and never changes an earlier default" (planning/data_layer_and_tests.md
§2A.3; `store._write_tables`), and `contig_uid` is the digest of the upper-case sequence (`ids.contig_uid`), so a
soft-masked copy is the same contig with other bytes (`tests/r2/test_micro_s0_l0seq.py::test_mr14_*`, and
`tests/r2/test_micro_s0_contracts.py::test_r09_a_later_stage_keeps_earlier_files`, whose later stage ingests a
soft-masked lambda that the first stage stores).

The fixture root is built with the library: `build-stage --depth called` (no library, no CRISPR, no release: no tool
beyond Python and DuckDB), every step a real run-step subprocess. Stage `OA` stores phage lambda (`L`); stage `ZQ`
holds a genome of `L` written with its first 2,000 bases in lower case (the same `contig_uid`, another `raw_digest`)
and a genome of the first 20,000 bases of phage T4 (`Y`), which it stores. Each test copies the root, may change one
table as a defect would (the run records' sizes and sha256 follow, so that only the defect is new), takes the two
snapshots as the driver's `before` and `after` would hold them, and runs the check's own `main`.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import sys
from datetime import datetime
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

from biobox import fasta, ids, profile, runstep
from biobox import stage as bs

pytestmark = pytest.mark.covers("R-17", "A-12", "R-09")  # ring R1 by directory (tests/conftest.py)

HERE = Path(__file__).resolve().parent
FIX = HERE.parent / "fixtures" / "micro_s0"
NOON_MONDAY = datetime(2026, 11, 9, 12, 0)

L = fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq  # phage lambda, as deposited (upper case)
L_MASKED = L[:2000].lower() + L[2000:]
Y = fasta.read(FIX / "NC_000866.4.fna.gz")[0].seq[:20_000]
CU_L, CU_Y = ids.contig_uid(L), ids.contig_uid(Y)
INDEX = "l0seq/contig_index/stage={}/part-0.parquet"

# today's message, which V2 got (/tmp/if2453/stages/v2/check/1/check.out)
BOTH = "contigs both stored by the stage and held by another"
ONE_DEFAULT = "contigs of the stage's index without exactly one default variant among all stages (DEC-26)"
REPEAT = "rows of the stage's index whose (contig_uid, raw_digest) another stage's index holds"
LEDGER = "ledger skips for another stage that differ from the stage holding the contig"


def _td():
    spec = importlib.util.spec_from_file_location("test_stage_driver", HERE / "test_stage_driver.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["test_stage_driver"] = mod
    spec.loader.exec_module(mod)
    return mod


def _write_stage(d: Path, stage: str, records: list[tuple[str, bytes]], genomes: list[bytes]) -> tuple[Path, Path]:
    """One FASTA object and the build-stage `--genomes` and `--objects` files of `stage` under `d`."""
    d.mkdir(parents=True, exist_ok=True)
    obj = d / f"{stage.lower()}.fna"
    obj.write_bytes(b"\n".join(ln for h, s in records for ln in fasta.format_record(h, s, 80)) + b"\n")
    g = d / f"genomes-{stage}.json"
    g.write_text(json.dumps([{"contig_uids": [cu.hex()], "genetic_code": 11, "kind": "virus", "circular": [],
                              "assembly_quality": "complete"} for cu in genomes]))
    o = d / f"objects-{stage}.json"
    o.write_text(json.dumps([{"path": str(obj), "source_name": "r17_fixture", "source_release": "1",
                              "license": "INSDC (public)", "redistributable": True}]))
    return g, o


def _runner(prof, root, stage, step, cmd, **kw) -> int:
    """run-step at a fixed hour and load, as tests/r2/test_shared_contigs.py runs build-stage."""
    return runstep.run_step(prof, root, stage, step, cmd, clock=lambda: NOON_MONDAY, load=lambda: 0.0, **kw)


@pytest.fixture(scope="module")
def built(tmp_path_factory) -> Path:
    d = tmp_path_factory.mktemp("r17")
    pdir = d / "prof"
    pdir.mkdir()
    (pdir / "profile.yaml").write_text(yaml.safe_dump({
        "max_threads_day": 4, "max_threads_night": 6, "mem_scope": "4G", "tmp": str(d / "tmp"),
        "durable": str(d / "durable"), "bulk": str(d / "bulk"), "nice_monolithic": 19, "nice_fanout": 15,
        "day_start_hour": 8, "day_end_hour": 20, "cpu_min_per_process": 300, "small_file_floor_bytes": 1000}))
    prof = profile.load(pdir / "profile.yaml")
    root = d / "root"
    for stage, records, genomes in (("OA", [("phage lambda", L)], [CU_L]),
                                    ("ZQ", [("lambda soft-masked", L_MASKED), ("T4 first 20 kb", Y)], [CU_L, CU_Y])):
        g, o = _write_stage(d / "inputs", stage, records, genomes)
        b = bs.Build(root=root, stage=stage, genomes=g, objects=o, depth="called", threads=1, tmp=d / "scratch")
        assert bs.build_stage(prof, b, run=_runner) == 0, stage
    return root


def _rows(root: Path, rel: str) -> list[dict]:
    return pq.read_table(root / rel).to_pylist()


def _rewrite(root: Path, rel: str, change) -> None:
    """Rewrite one published Parquet file with `change(rows)`, in its own schema, and make every run record that
    lists it name its new size and sha256 (a defect of the table, not a size that differs from its record)."""
    p = root / rel
    t = pq.read_table(p)
    pq.write_table(pa.Table.from_pylist(change(t.to_pylist()), schema=t.schema), p)
    n, h = p.stat().st_size, hashlib.sha256(p.read_bytes()).hexdigest()
    hits = 0
    for rec_p in sorted((root / "runs").glob("stage=*/*.json")):
        rec = json.loads(rec_p.read_text())
        for o in rec.get("outputs", []):
            if o["path"] == rel:
                o["bytes"] = n
                if "sha256" in o:
                    o["sha256"] = h
                hits += 1
        rec_p.write_text(json.dumps(rec))
    assert hits == 1, f"{rel}: listed by {hits} run records"


def _check(built: Path, tmp: Path, monkeypatch, change=None) -> tuple[int, dict]:
    """The check of stage ZQ on a copy of the built root, after `change(root)`, with the snapshots taken after it
    (`before` without ZQ's own files, as the driver's `before` lists the root)."""
    td = _td()
    root = tmp / "root"
    shutil.copytree(built, root, symlinks=True)
    if change is not None:
        change(root)
    before, after = tmp / "before", tmp / "after"
    before.mkdir()
    after.mkdir()
    lines = sorted(f"{p.relative_to(root)}\t{p.stat().st_size}\t{p.stat().st_mtime}" for p in root.rglob("*")
                   if p.is_file())
    (after / "listing.tsv").write_text("".join(ln + "\n" for ln in lines))
    (before / "listing.tsv").write_text("".join(ln + "\n" for ln in lines if "/stage=ZQ/" not in "/" + ln.split("\t")[0]))
    conf = tmp / "ZQ.conf"
    conf.write_text(td.conf_text("ZQ", tmp / "inputs", td.GOOD, depth="called", library=None, crispr="no"))
    prof = tmp / "profile.yaml"
    prof.write_text(td.PROFILE.replace("tmp: /x", f"tmp: {tmp}"))
    bj = tmp / "before.json"
    bj.write_text(json.dumps({"fasta_bytes": len(L) + len(Y), "held_bases": len(L)}))
    r = {"root": root, "conf": conf, "before": before, "after": after, "profile": prof, "before_json": bj,
         "out": tmp / "out"}
    return td.run_check(r, monkeypatch)


def test_the_fixture_is_a_later_stage_holding_a_new_raw_variant_of_a_contig_the_earlier_stage_stores(built):
    """The oracle for the fixture, independent of the check: the tables as ingest wrote them, read with pyarrow, and
    the digests computed here with hashlib (`ids.raw_digest` is the sha256 of the bytes, AGENTS §3.2)."""
    oa, zq = _rows(built, INDEX.format("OA")), _rows(built, INDEX.format("ZQ"))
    assert [(r["contig_uid"], r["raw_digest"], r["is_default_variant"]) for r in oa] == \
        [(CU_L, hashlib.sha256(L).digest(), True)]
    assert sorted((r["contig_uid"], r["raw_digest"], r["contig_len"], r["is_default_variant"]) for r in zq) == sorted([
        (CU_L, hashlib.sha256(L_MASKED).digest(), len(L), False), (CU_Y, hashlib.sha256(Y).digest(), len(Y), True)])
    assert CU_L == ids.contig_uid(L_MASKED) and hashlib.sha256(L_MASKED).digest() != hashlib.sha256(L).digest()
    led = [r for p in sorted(built.glob("ledger/stage=ZQ/step=genecall/*.parquet")) for r in pq.read_table(p).to_pylist()]
    assert {(r["object_key"], r["reason_code"]) for r in led if r["object_key"] == CU_L} == {(CU_L, "owned_by:OA")}


def test_a_non_default_variant_of_a_contig_another_stage_stores_passes_and_is_counted_apart(built, tmp_path,
                                                                                              monkeypatch):
    rc, s = _check(built, tmp_path, monkeypatch)
    assert rc == 0, s["failed"]
    assert s["failed"] == []
    c = s["contigs"]
    assert (c["stored"], c["index_contigs"], c["index_rows"]) == (1, 2, 2)
    assert c["bases_stored"] == len(Y)
    assert c["variants_of_contigs_held_by_other_stages"] == {"contigs": 1, "rows": 1, "bases": len(L),
                                                             "by_stage": {"OA": 1}}
    assert c["non_default_rows_of_stored_contigs"] == 0
    assert c["held_by_other_stages"] == 1 and c["held_by_stage"] == {"OA": 1}
    assert c["ledger_skipped_other"] == 1
    assert c[BOTH] == 0 and c[ONE_DEFAULT] == 0 and c[REPEAT] == 0
    assert s["coverage"]["stored_with_genes"] == 1 and s["coverage"]["contigs_with_genes_not_stored_here"] == 0


def _set_default(stage: str, cu: bytes, value: bool):
    def change(root: Path) -> None:
        _rewrite(root, INDEX.format(stage), lambda rows: [dict(r, is_default_variant=value) if r["contig_uid"] == cu
                                                          else r for r in rows])
    return change


def test_a_default_variant_row_of_a_contig_another_stage_stores_still_fails_with_the_same_message(built, tmp_path,
                                                                                                   monkeypatch):
    rc, s = _check(built, tmp_path, monkeypatch, _set_default("ZQ", CU_L, True))
    assert rc == 1
    assert f"{BOTH}: 1" in s["failed"], s["failed"]
    assert f"{ONE_DEFAULT}: 1" in s["failed"], s["failed"]


def test_a_variant_of_a_contig_that_no_stage_holds_as_default_fails(built, tmp_path, monkeypatch):
    """A non-default row of a contig no stage stores. (Built by removing `OA`'s default of `L` instead, the rate
    helper's `genecall.plan` raises on the contig of `genome_input` that no stage stores before the check's own
    assertions run: a crash, exit 1 all the same.)"""
    z = fasta.read(FIX / "NC_000866.4.fna.gz")[0].seq[20_000:30_000]

    def stray(root: Path) -> None:
        _rewrite(root, INDEX.format("ZQ"), lambda rows: sorted(
            rows + [{"contig_uid": ids.contig_uid(z), "raw_digest": hashlib.sha256(z).digest(), "contig_len": len(z),
                     "is_default_variant": False}], key=lambda r: (r["contig_uid"], r["raw_digest"])))

    rc, s = _check(built, tmp_path, monkeypatch, stray)
    assert rc == 1
    assert f"{ONE_DEFAULT}: 1" in s["failed"], s["failed"]


def test_two_default_variants_of_one_contig_in_the_stage_fail(built, tmp_path, monkeypatch):
    def second_default(root: Path) -> None:
        _rewrite(root, INDEX.format("ZQ"), lambda rows: sorted(
            rows + [{"contig_uid": CU_Y, "raw_digest": hashlib.sha256(Y.lower()).digest(), "contig_len": len(Y),
                     "is_default_variant": True}], key=lambda r: (r["contig_uid"], r["raw_digest"])))

    rc, s = _check(built, tmp_path, monkeypatch, second_default)
    assert rc == 1
    assert f"{ONE_DEFAULT}: 1" in s["failed"], s["failed"]


def test_a_variant_row_that_repeats_the_raw_digest_the_earlier_stage_has_for_the_contig_fails(built, tmp_path,
                                                                                                monkeypatch):
    def repeat(root: Path) -> None:
        _rewrite(root, INDEX.format("ZQ"), lambda rows: [dict(r, raw_digest=hashlib.sha256(L).digest())
                                                         if r["contig_uid"] == CU_L else r for r in rows])

    rc, s = _check(built, tmp_path, monkeypatch, repeat)
    assert rc == 1
    assert f"{REPEAT}: 1" in s["failed"], s["failed"]


def test_the_ledger_must_still_skip_a_contig_whose_variant_the_stage_holds(built, tmp_path, monkeypatch):
    def drop(root: Path) -> None:
        (p,) = sorted(root.glob("ledger/stage=ZQ/step=genecall/*.parquet"))
        _rewrite(root, p.relative_to(root).as_posix(), lambda rows: [r for r in rows if r["object_key"] != CU_L])

    rc, s = _check(built, tmp_path, monkeypatch, drop)
    assert rc == 1
    assert f"{LEDGER}: 1" in s["failed"], s["failed"]


def test_a_non_default_row_without_a_raw_digest_to_compare_fails(tmp_path, monkeypatch):
    """The hand-built root of tests/r1/test_stage_driver.py writes `contig_index` without `raw_digest`; a stage there
    with a non-default row cannot be checked against the earlier stage's digests, so it fails."""
    td = _td()

    def variant_of_h(S):
        S["T"]["ci"][1].append({"contig_uid": S["u"]["H"], "contig_len": 3000, "is_default_variant": False})

    from biobox.l1 import crispr

    monkeypatch.setattr(crispr, "minced_version", lambda: "0.4.2")
    rc, s = td.run_check(td.build_root(tmp_path, variant_of_h), monkeypatch)
    assert rc == 1
    assert any("non-default rows in the stage's index, but the index files of" in f and "raw_digest" in f
               for f in s["failed"]), s["failed"]
