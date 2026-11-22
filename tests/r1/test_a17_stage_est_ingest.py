"""Ring R1: the stage driver judges ingest by build-stage's own estimate (card A-17, 2026-11-21).

`stage_est.py before` must record the ingest estimate `build-stage` plans by (`stage._batched`: `stage.ingest_cpu_min`
over the bytes of FASTA, the bases another stage holds, and card B-61's two terms, the bytes out of `.gz` objects and
the bytes in records under the least contig length, counted in the one read `stage.object_sizes(objects,
min_contig_len)` makes of the objects), whether ingest runs in batches by the same comparison, and the first-pass cut
`build-stage` would make; `stage_est.py rates` judges the ingest group against that estimate, read from
`estimate/before.json`.

The constructed stage: two `.gz` objects and one plain one, records over and under a least length of 5,000 bp (one at
4,999 and one at exactly 5,000), sequences wrapped at 60 so that bytes of FASTA and bases differ, and one contig that
an earlier stage `OA` holds. The oracle for the terms is what this file wrote (the bytes of each record, the
uncompressed bytes of each object), and for the figure and the batch decision the code's own functions, called with
the arguments `stage._batched` passes them on a root where the stage's `genome_input` is published, as it is when
`build-stage` decides. Nothing of the pod is read; no tool beyond Python and DuckDB is needed.
"""
from __future__ import annotations

import argparse
import dataclasses
import gzip
import importlib.util
import json
import random
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import ids
from biobox import profile as bprofile
from biobox import stage as bs
from biobox.l0seq import store
from biobox.registry.schemas import table

pytestmark = pytest.mark.covers("A-17")  # ring R1 by directory (tests/conftest.py)

REPO = Path(__file__).resolve().parents[2]
BSP = REPO / "scripts" / "stress" / "build_stage_pod"
MIN = 5000
NOTE = "planning/measurements/2026-11-20_stage_V1.md"  # a note of the tree, as test_stage_driver.py's confs name
PROFILE = """max_threads_day: 40
max_threads_night: 56
mem_scope: 200G
tmp: {tmp}
durable: /x
bulk: /x
nice_monolithic: 19
nice_fanout: 15
day_start_hour: 8
day_end_hour: 20
cpu_min_per_process: {budget}
cluster_cpu_min_max: 20000
small_file_floor_bytes: 67108864
require_mem_cap: true
"""


def _est():
    spec = importlib.util.spec_from_file_location("stage_est", BSP / "stage_est.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["stage_est"] = mod  # a dataclass looks its module up by name
    spec.loader.exec_module(mod)
    return mod


def _seq(rng: random.Random, n: int) -> bytes:
    return rng.randbytes(n).translate(bytes(b"ACGT"[i % 4] for i in range(256)))


def _record(name: str, s: bytes) -> bytes:
    return b">" + name.encode() + b" a header\n" + b"".join(s[i:i + 60] + b"\n" for i in range(0, len(s), 60))


@dataclasses.dataclass
class Stage:
    inputs: Path
    objects: Path
    root: Path
    conf: Path
    fasta: int  # bytes of FASTA in every object (uncompressed), counted as written
    gz: int  # of them, in the .gz objects
    below: int  # of them, in records under MIN bases
    held: int  # bases of the contig OA holds
    kept: dict  # name: sequence of every record at or over MIN


def build_stage(tmp: Path, *, minlen: int | None = MIN, gz: bool = True, crispr: str = "no") -> Stage:
    """The constructed stage `ZZ1` (module docstring) and a root with stage `OA`, which holds contig `held` and is
    gene-called, and with `ZZ1`'s `genome_input` published (as after the inputs step; `stage_est.py` does not read it)."""
    rng = random.Random(1717)
    held = _seq(rng, 5_500)
    shorts = [_seq(rng, n) for n in [4_999, 300] + [4_000] * 499]  # about 2 MB under MIN: a term over 0.001 CPU-min
    objs = {  # file: records
        "a.fna.gz" if gz else "a.fna": [("l1", _seq(rng, 7_000)), ("s0", shorts[0]), ("h", held),
                                        *((f"s{i}", s) for i, s in enumerate(shorts[2:252], 2))],
        "b.fna.gz" if gz else "b.fna": [("l3", _seq(rng, 5_000)), ("s1", shorts[1]), ("l2", _seq(rng, 6_000)),
                                        *((f"s{i}", s) for i, s in enumerate(shorts[252:], 252))],
        "c.fna": [("l4", _seq(rng, 8_000)), ("x", _seq(rng, 30_000))],  # a plain object: not in the gzip term
    }
    inp = tmp / "inputs"
    inp.mkdir(parents=True)
    fasta = gzb = below = 0
    kept = {}
    for name, recs in objs.items():
        data = b"".join(_record(h, s) for h, s in recs)
        fasta += len(data)
        if name.endswith(".gz"):
            gzb += len(data)
            with gzip.open(inp / name, "wb") as fh:
                fh.write(data)
        else:
            (inp / name).write_bytes(data)
        for h, s in recs:
            if minlen is not None and len(s) < minlen:
                below += len(_record(h, s))
            else:
                kept[h] = s
    genomes = [{"contig_uids": [ids.contig_uid(s).hex()], "genetic_code": 11, "kind": "virus", "circular": [],
                "assembly_quality": "contig"} for _, s in sorted(kept.items())]
    (inp / "genomes.json").write_text(json.dumps(genomes))
    objects = inp / "objects.json"
    objects.write_text(json.dumps([{"path": str(inp / n), "source_name": "t", "source_release": "1", "license": "t",
                                    "redistributable": True} for n in objs]))
    root = tmp / "root"
    oa = root / "l0seq" / "contig_index" / "stage=OA"
    oa.mkdir(parents=True)
    pq.write_table(pa.table({"contig_uid": pa.array([ids.contig_uid(held)], pa.binary(16)),
                             "contig_len": pa.array([len(held)], pa.int64()),
                             "is_default_variant": pa.array([True])}), oa / "part-0.parquet")
    (root / "runs" / "stage=OA").mkdir(parents=True)
    (root / "runs" / "stage=OA" / "genecall.json").write_text("{}")  # gene-called: not a parked stage
    gi = root / table("genome_input").location / "stage=ZZ1"
    gi.mkdir(parents=True)
    pq.write_table(pa.table({"contig_uid": pa.array(sorted(ids.contig_uid(s) for s in kept.values()), pa.binary(16))}),
                   gi / "part-0.parquet")
    lines = ["STAGE=ZZ1", f"NOTE={NOTE}", f"INPUTS={inp}", f"SUM={'a' * 64} genomes.json", f"SUM={'b' * 64} objects.json",
             "DEPTH=called", f"CRISPR={crispr}"] + ([f"MIN_CONTIG_LEN={minlen}"] if minlen is not None else [])
    conf = tmp / "ZZ1.conf"
    conf.write_text("\n".join(lines) + "\n")
    return Stage(inp, objects, root, conf, fasta, gzb, below, len(held), kept)


def profile(tmp: Path, budget: int) -> Path:
    p = tmp / f"profile-{budget}.yaml"
    p.write_text(PROFILE.format(tmp=tmp, budget=budget))
    return p


def run_before(est, s: Stage, prof: Path, tmp: Path, capsys) -> tuple[int, dict]:
    a = argparse.Namespace(conf=str(s.conf), stage="ZZ1", root=str(s.root), profile=str(prof), inputs=str(s.inputs),
                           tmp=str(tmp / "est"), out_plan=None)
    rc = est.before(a)
    return rc, json.loads(capsys.readouterr().out)


def build(s: Stage, tmp: Path, minlen: int | None = MIN) -> bs.Build:
    return bs.Build(root=s.root, stage="ZZ1", genomes=s.inputs / "genomes.json", objects=s.objects, tmp=tmp / "bs",
                    min_contig_len=minlen, depth="called")


def code_estimate(s: Stage, tmp: Path, minlen: int | None = MIN) -> float:
    """build-stage's one-process figure: what `stage._batched` computes, with its own calls and arguments."""
    sizes = bs.object_sizes(s.objects, minlen)
    return bs.ingest_cpu_min(sum(f for _, _, f in sizes), bs.held_bases(s.root, "ZZ1", tmp / "held"),
                             gz_bytes=bs.gz_bytes(sizes), below_bytes=bs.below_bytes(s.objects, minlen))


# --- before: the code's own figure -------------------------------------------------------------------------------

def test_before_records_build_stages_own_ingest_estimate_for_gzip_objects_a_least_length_and_a_held_contig(
        tmp_path, capsys):
    est = _est()
    s = build_stage(tmp_path)
    rc, doc = run_before(est, s, profile(tmp_path, 300), tmp_path, capsys)
    assert rc == 0
    # the terms against what this file wrote
    assert doc["fasta_bytes"] == s.fasta and doc["held_bases"] == s.held
    assert doc["gz_bytes"] == s.gz and doc["below_bytes"] == s.below
    # ... and against the code's own counts (the arguments `stage._batched` passes)
    sizes = bs.object_sizes(s.objects, MIN)
    assert (bs.gz_bytes(sizes), bs.below_bytes(s.objects, MIN)) == (s.gz, s.below)
    assert bs.held_bases(s.root, "ZZ1", tmp_path / "held") == s.held
    code = code_estimate(s, tmp_path)
    assert doc["ingest_cpu_min"] == round(code, 3), "stage_est.py's ingest figure is build-stage's own"
    # the terms move the figure here: card B-52's figure would not pass
    assert round(bs.ingest_cpu_min(s.fasta, s.held), 3) != doc["ingest_cpu_min"]
    f = doc["ingest_formula"]
    assert f"{s.below} x {bs.INGEST_BELOW_CPU_MIN}" in f and f"{s.gz} x {bs.INGEST_GZ_CPU_MIN}" in f
    assert f"max(0, {s.fasta} - {s.held} - {s.below})" in f


@pytest.mark.parametrize("budget", [0, 1])
def test_before_decides_batches_as_build_stage_does_and_records_its_first_pass_cut(tmp_path, capsys, budget):
    """`stage._batched` itself on the same stage and profile: at a budget of 0 every estimate is over the margin (and
    build-stage would refuse the stage: exit 7, the document still written), at 1 this stage is under it."""
    est = _est()
    s = build_stage(tmp_path)
    prof = profile(tmp_path, budget)
    rc, doc = run_before(est, s, prof, tmp_path, capsys)
    p = bprofile.load(prof)
    batched = bs._batched(p, build(s, tmp_path), s.root)
    assert batched is (budget == 0)
    assert doc["ingest_in_batches"] is batched
    assert rc == (7 if budget == 0 else 0)
    if batched:
        per = bs.pass1_bytes(p, build(s, tmp_path), s.root)
        from biobox.cli import _source_objects

        assert doc["ingest_pass1_bytes"] == per
        assert doc["ingest_pass1_shards"] == len(store.pass1_shards(_source_objects(str(s.objects.resolve())), per))
    else:
        assert "ingest_pass1_bytes" not in doc and "ingest_pass1_shards" not in doc


def test_before_reads_each_object_once_for_its_sizes_and_terms(tmp_path, capsys, monkeypatch):
    """No second read of the objects for the terms (on SPX1 one read is 257 GB of FASTA, about 28 CPU-minutes): the
    bytes under the least length come out of the read `object_sizes` makes (`fasta_counts`), not a read of their own."""
    est = _est()
    s = build_stage(tmp_path)
    monkeypatch.setattr(bs, "_SIZES", {})
    monkeypatch.setattr(bs, "_BELOW", {})
    seen = {"counts": [], "bytes": []}
    real_counts, real_bytes = bs.fasta_counts, bs.fasta_bytes
    monkeypatch.setattr(bs, "fasta_counts", lambda p, n: seen["counts"].append(Path(p).name) or real_counts(p, n))
    monkeypatch.setattr(bs, "fasta_bytes", lambda p: seen["bytes"].append(Path(p).name) or real_bytes(p))
    rc, doc = run_before(est, s, profile(tmp_path, 300), tmp_path, capsys)
    assert rc == 0 and doc["below_bytes"] == s.below
    assert sorted(seen["counts"]) == ["a.fna.gz", "b.fna.gz", "c.fna"] and seen["bytes"] == []


def test_a_stage_of_plain_fasta_without_a_least_length_is_estimated_and_said_as_before(tmp_path, capsys):
    """Card B-52's figure and formula, and no new key in the document (V2's and GV1's documents equal before and
    after, shown on their staged inputs on the pod)."""
    est = _est()
    s = build_stage(tmp_path, minlen=None, gz=False)
    assert s.gz == 0 and s.below == 0
    rc, doc = run_before(est, s, profile(tmp_path, 300), tmp_path, capsys)
    assert rc == 0
    assert doc["ingest_cpu_min"] == round(bs.ingest_cpu_min(s.fasta, s.held), 3) == round(
        code_estimate(s, tmp_path, None), 3)
    assert doc["ingest_formula"] == (f"{bs.INGEST_FIXED_CPU_MIN} + max(0, {s.fasta} - {s.held}) x {bs.INGEST_CPU_MIN}"
                                     f" / {bs.INGEST_BP}")
    assert not {"gz_bytes", "below_bytes", "ingest_pass1_bytes", "ingest_pass1_shards"} & set(doc)


# --- the confs of G1 and SPX1, from their notes ------------------------------------------------------------------

def test_the_confs_of_g1_and_spx1_hold_the_values_of_their_notes():
    est = _est()
    g1 = est.read_conf(BSP / "stages" / "G1.conf", "G1", tree=REPO)
    assert (g1.depth, g1.library, g1.crispr, g1.min_contig_len) == ("annotated", ("pfam", "38.2", "gathering"), True,
                                                                    None)
    assert g1.after == () and g1.judge == {}
    assert g1.accepted == {("crispr-minced-0.4.2-w8-nr3-r23x47-s26x50", 66.117, 38.363): "DEC-68"}
    assert {s: n for s, (n, _) in g1.expect.items()} == {"P2": 1, "P4": 6, "H1": 1}
    spx = est.read_conf(BSP / "stages" / "SPX1.conf", "SPX1", tree=REPO)
    assert (spx.depth, spx.library, spx.crispr, spx.min_contig_len) == ("called", None, True, 5000)
    assert spx.after == ("G1",) and spx.judge == {} and not spx.accepted
    assert {s: n for s, (n, _) in spx.expect.items()} == {"G1": 12052}
    assert {n for _, n in spx.sums} == {"genomes.json", "objects.json", "samples.tsv", "input_report.json",
                                        "assemblies.sha256"}


# --- rates: the ingest group against that estimate ---------------------------------------------------------------

# SPX1's counts by the code (card B-61, planning/measurements/2026-11-21_ingest_terms.md §4, *measured* there):
# 256,948,810,775 bytes of FASTA, every one out of a .gz object; 220,918,826,783 of them in records under 5,000 bp;
# none held. The one-process figure is 245.9 CPU-minutes (card B-52's 1,056.1).
SPX = {"fasta_bytes": 256_948_810_775, "held_bases": 0, "gz_bytes": 256_948_810_775, "below_bytes": 220_918_826_783,
       "ingest_in_batches": True}
SPX_EST = 245.9


def _rates_root(tmp: Path, cpu: dict) -> Path:
    root = tmp / "root"
    d = root / "runs" / "stage=ZZ1"
    d.mkdir(parents=True)
    for step, c in cpu.items():
        (d / f"{step}.json").write_text(json.dumps({"stage": "ZZ1", "step": step, "exit_code": 0, "cpu_min": c}))
    return root


def _conf(tmp: Path) -> Path:
    c = tmp / "ZZ1.conf"
    c.write_text(f"STAGE=ZZ1\nNOTE={NOTE}\nINPUTS={tmp / 'inputs'}\nSUM={'a' * 64} genomes.json\n"
                 f"SUM={'b' * 64} objects.json\nDEPTH=called\nCRISPR=no\nMIN_CONTIG_LEN={MIN}\n")
    return c


# the pieces of a batched ingest: two pass-1 shards, the plan, two batches, the merge (`stage.ingest_steps`)
PIECES = ("ingest.pass1-00000", "ingest.pass1-00001", "ingest.plan", "ingest.batch-0000", "ingest.batch-0001", "ingest")


@pytest.mark.parametrize("cpus,ratio,stops", [
    ((150.0, 5.0, 1.0, 60.0, 60.0, 4.0), 1.139, False),  # 280 against 245.9: inside 20%
    ((150.0, 5.0, 1.0, 100.0, 120.0, 4.0), 1.545, True),  # 380: over DEC-64's 1.50, which 1,056.1 would never see
])
def test_rates_judges_every_ingest_piece_as_one_group_against_the_estimate_in_before_json(tmp_path, cpus, ratio, stops):
    est = _est()
    root = _rates_root(tmp_path, dict(zip(PIECES, cpus)))
    bj = tmp_path / "before.json"
    bj.write_text(json.dumps(SPX))
    doc = est.rates_doc(est.read_conf(_conf(tmp_path), "ZZ1"), root, bprofile.load(profile(tmp_path, 300)),
                        tmp_path / "no-objects.json", str(bj), tmp_path / "t")  # the objects are never read
    ing = next(x for x in doc["rows"] if x["step"] == "ingest")
    assert ing["records"] == 6 and ing["cpu_min"] == pytest.approx(sum(cpus))
    code = bs.ingest_cpu_min(SPX["fasta_bytes"], 0, gz_bytes=SPX["gz_bytes"], below_bytes=SPX["below_bytes"])
    assert ing["estimate"] == round(code, 3) and round(code, 1) == SPX_EST
    assert ing["ratio"] == ratio and doc["not_judged"] == []
    assert f"gz_bytes={SPX['gz_bytes']}" in ing["basis"] and f"below_bytes={SPX['below_bytes']}" in ing["basis"]
    assert "in batches" in ing["basis"]
    assert est.rates_exit(doc) == (3 if stops else 0)
    assert (doc["over_20pct"] == ["ingest"]) is stops


def test_rates_without_before_json_measures_the_terms_with_the_codes_own_read(tmp_path):
    est = _est()
    s = build_stage(tmp_path)
    d = s.root / "runs" / "stage=ZZ1"
    d.mkdir(parents=True)
    (d / "ingest.json").write_text(json.dumps({"stage": "ZZ1", "step": "ingest", "exit_code": 0, "cpu_min": 0.1}))
    doc = est.rates_doc(est.read_conf(s.conf, "ZZ1"), s.root, bprofile.load(profile(tmp_path, 300)), s.objects,
                        None, tmp_path / "t")
    ing = next(x for x in doc["rows"] if x["step"] == "ingest")
    assert ing["estimate"] == round(code_estimate(s, tmp_path), 3)
    assert f"gz_bytes={s.gz}, below_bytes={s.below})" in ing["basis"]
