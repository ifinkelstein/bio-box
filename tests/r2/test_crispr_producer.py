"""Ring R2: the CRISPR array producer (card E-12, operation O-ARRAY; `biobox.l1.crispr`) against two oracles.

Oracle 1 (AGENTS §5.2): MinCED's own default text report for the same contigs, run here through the conda wrapper
`minced` with no option (so it also checks that the producer's explicit parameters are MinCED's defaults) and parsed
by a parser of this file, not the producer's. Every array the binary reports is one array of the producer with the
same repeats and spacers (positions converted from 1-based, sequences equal), and the producer has no other array.
On real genomes too, when the golden set is fetched (`BIOBOX_GOLDEN_DIR`; `scripts/pod_pytest.sh --golden`).

Oracle 2: contigs with planted arrays whose repeats and spacers are known by construction
(`tests/r1/test_crispr_plant.py`): a plain array, one inserted reverse-complemented, two arrays on one contig, an
array across the origin of a circular contig (the same array as on a linear twin, rotated), and a circular contig
with none. Also: a contig reverse-complemented gives the same array mirrored (MR-03).

Then: bytes independent of the number of MinCED processes, the batches and the genome order, and of a second run
(MR-01, MR-16); a run under `bio-box run-step` equals the direct call; refusals before work; the ledger's tool
versions and balance; a contig another stage owns, or another genome of the stage, is skipped and ledgered (DEC-26,
DEC-27); on the micro-S0 fixture the producer changes no other file (MR-20); its arrays go on a tape built with its
producer (kind crispr_array, no protein); a lookup by contig, and by (contig, feature_uid), prunes (AGENTS §4.4).
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import duckdb
import pyarrow.parquet as pq
import pytest
import yaml

from biobox import catalog, fasta, ids, io, ledger, provenance
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l1 import igr
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l1.genecall import run as call_genes
from biobox.l3 import tapes
from biobox.registry.schemas import conformance_problems, table

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "r1"))
import test_builders_out as c08  # noqa: E402
import test_crispr_plant as plant  # noqa: E402

MINCED = shutil.which("minced")
pytestmark = [pytest.mark.r2, pytest.mark.covers("E-12"),
              pytest.mark.skipif(MINCED is None, reason="minced (MinCED, a Java tool) not on PATH; it is in the "
                                                        "bio-box image, not in the uv environment of GitHub Actions")]

CRISPR = ids.KIND["crispr_array"]
FIX = c08.FIX
CONTRACT = c08.CONTRACT


def _crispr():
    from biobox.l1 import crispr  # imported here so that a missing module fails the tests, not the collection

    return crispr


# --- an independent reading of the binary's text report ------------------------------------------------------------


_SEQ = re.compile(r"^Sequence '(\S+)' \((\d+) bp\)")
_RANGE = re.compile(r"^CRISPR (\d+)\s+Range: (\d+) - (\d+)")


def oracle_text(path: Path) -> dict[str, list[dict]]:
    """MinCED's default report: per record name, its arrays, each with its 1-based range and its rows (position,
    repeat, spacer or None for the last repeat)."""
    out: dict[str, list[dict]] = {}
    name = arr = None
    for line in path.read_text().splitlines():
        if m := _SEQ.match(line):
            name = m[1]
            out.setdefault(name, [])
        elif m := _RANGE.match(line):
            arr = {"range": (int(m[2]), int(m[3])), "rows": []}
            out[name].append(arr)
        elif arr is not None and line[:1].isdigit():
            f = line.split()
            arr["rows"].append((int(f[0]), f[1], f[2] if len(f) > 2 and f[2] != "[" else None))
    return out


def oracle_elements(arr: dict) -> list[tuple[str, int, int, str]]:
    """(element, start, end, seq), 0-based half-open on the record, of one reported array."""
    out = []
    for pos, rep, sp in arr["rows"]:
        out.append(("repeat", pos - 1, pos - 1 + len(rep), rep))
        if sp is not None:
            out.append(("spacer", pos - 1 + len(rep), pos - 1 + len(rep) + len(sp), sp))
    return out


def run_minced(records: dict[str, bytes], d: Path) -> dict[str, list[dict]]:
    """The binary through its wrapper, no option, on the records (name -> sequence)."""
    d.mkdir(parents=True, exist_ok=True)
    fna = d / "in.fna"
    fna.write_bytes(b"".join(b">" + k.encode() + b"\n" + v + b"\n" for k, v in records.items()))
    subprocess.run([MINCED, str(fna), str(d / "out.txt")], check=True, capture_output=True, cwd=d)
    return oracle_text(d / "out.txt")


# --- the stores -----------------------------------------------------------------------------------------------------


def _rows(root: Path, rel: str) -> list[dict]:
    return [r for f in io.part_files(root / rel) for r in pq.read_table(f, partitioning=None).to_pylist()]


def digests(root: Path) -> dict[str, str]:
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(Path(root).rglob("*")) if p.is_file()}


def write_objects(d: Path, contigs: dict[str, plant.Contig]) -> list[SourceObject]:
    d.mkdir(parents=True, exist_ok=True)
    objs = []
    for name, c in contigs.items():
        path = d / f"{name}.fna"
        path.write_bytes(b"\n".join(fasta.format_record(f"planted_{name}", c.seq, 80)) + b"\n")
        objs.append(SourceObject(path, "planted", "1", "INSDC (public)", True))
    return objs


def genomes_of(contigs: dict[str, plant.Contig]) -> list[GenomeInput]:
    out = []
    for c in contigs.values():
        cu = ids.contig_uid(c.seq)
        out.append(GenomeInput((cu,), 11, kind="plasmid", circular=(cu,) if c.circular else ()))
    return out


def arrays_of(rows_f: list[dict], rows_e: list[dict]) -> dict[bytes, list[tuple]]:
    """Per contig, its arrays as (start, end, [(element, ordinal, start, end, seq)]) sorted by start."""
    els: dict[bytes, list[tuple]] = {}
    for r in rows_e:
        els.setdefault(r["feature_uid"], []).append((r["element"], r["ordinal"], r["start"], r["end"], r["seq"]))
    out: dict[bytes, list[tuple]] = {}
    for f in rows_f:
        out.setdefault(f["contig_uid"], []).append((f["start"], f["end"], sorted(els[f["feature_uid"]],
                                                                                key=lambda t: t[1])))
    return {k: sorted(v) for k, v in out.items()}


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """The planted contigs in one stage `P`, gene-called (for the tape test), then the producer with one process."""
    base = tmp_path_factory.mktemp("crispr")
    root, scratch = base / "root", base / "scratch"
    scratch.mkdir()
    cs = plant.contigs()
    ingest(write_objects(base / "objects", cs), root, "P")
    genomes = genomes_of(cs)
    write_genome_input(root, "P", genomes)
    call_genes(root, "P", genomes, tmp=scratch)
    before = base / "before"
    shutil.copytree(root, before)
    report = _crispr().run(root, "P", 1, tmp=scratch)
    p = report["producer"]
    new = {k: h for k, h in digests(root).items() if k not in digests(before)}
    return {"root": root, "before": before, "scratch": scratch, "base": base, "contigs": cs, "genomes": genomes,
            "report": report, "producer": p, "files": new,
            "features": _rows(root, f"l1/feature/producer={p}/stage=P"),
            "elements": _rows(root, f"l1/crispr/producer={p}/stage=P")}


# --- oracle 2: the planted truth ------------------------------------------------------------------------------------


@pytest.mark.covers("O-ARRAY", "C-05", "MR-04")
def test_the_stored_arrays_are_the_planted_ones(world):
    got = arrays_of(world["features"], world["elements"])
    want = {}
    for c in world["contigs"].values():
        cu = ids.contig_uid(c.seq)
        if c.arrays:
            want[cu] = sorted((a.start, a.end, [(e.kind, e.ordinal, e.start, e.end, e.seq.decode()) for e in a.elements])
                              for a in c.arrays)
    assert got == want
    across = world["contigs"]["across"]
    n = len(across.seq)
    (a,) = got[ids.contig_uid(across.seq)]
    assert a[0] < n < a[1]  # found once, with end past the contig (AGENTS §3.1)
    assert ids.contig_uid(world["contigs"]["none"].seq) not in got
    assert world["report"]["arrays"] == 6 and world["report"]["rejected"] == 0  # 1 + 1 + 2 + 1 + 1, none on `none`


@pytest.mark.covers("C-01", "C-02", "LY-01")
def test_rows_are_keyed_conform_and_every_seq_is_the_contigs_bases(world):
    root, p = world["root"], world["producer"]
    for name, rel in (("l1_crispr", f"l1/crispr/producer={p}/stage=P"), ("feature", f"l1/feature/producer={p}/stage=P")):
        (f,) = io.part_files(root / rel)
        t = pq.read_table(f, partitioning=None)
        assert t.num_rows > 0 and conformance_problems(t, table(name)) == []
        md = pq.ParquetFile(f).metadata
        for col in table(name).sort_key:
            i = next(k for k in range(md.num_columns) if md.schema.column(k).path == col)
            assert md.row_group(0).column(i).statistics.has_min_max, (name, col)
    by_cu = {ids.contig_uid(c.seq): c for c in world["contigs"].values()}
    feats = {f["feature_uid"]: f for f in world["features"]}
    for f in world["features"]:
        assert f["kind"] == CRISPR and f["strand"] == 1
        assert f["feature_uid"] == ids.feature_uid(f["contig_uid"], f["start"], f["end"], 1, CRISPR)
    for r in world["elements"]:
        c = by_cu[r["contig_uid"]]
        assert r["seq"].encode() == plant.bases(c, r["start"], r["end"])
        assert r["region_uid"] == ids.region_uid(r["contig_uid"], r["start"], r["end"], 1)
        assert r["strand"] == 1 and r["contig_len"] == len(c.seq) and 0 <= r["start"] < len(c.seq)
        assert r["feature_uid"] in feats


@pytest.mark.covers("MR-03")
def test_a_reverse_complemented_contig_gives_the_same_array_mirrored(world, tmp_path):
    c = world["contigs"]["plain"]
    rc = plant.Contig(plant.revcomp(c.seq), [])
    root = tmp_path / "root"
    ingest(write_objects(tmp_path / "o", {"rc": rc}), root, "R")
    write_genome_input(root, "R", genomes_of({"rc": rc}))
    rep = _crispr().run(root, "R", 1, tmp=world["scratch"])
    p = rep["producer"]
    (arr,) = arrays_of(_rows(root, f"l1/feature/producer={p}/stage=R"), _rows(root, f"l1/crispr/producer={p}/stage=R")
                       ).values()
    n = len(c.seq)
    (a,) = c.arrays
    assert [(s, e) for s, e, _ in arr] == [(n - a.end, n - a.start)]
    assert [x[4] for x in arr[0][2]] == [plant.revcomp(e.seq).decode() for e in reversed(a.elements)]


# --- oracle 1: the binary's own report ------------------------------------------------------------------------------


@pytest.mark.covers("O-ARRAY")
def test_every_array_equals_the_binarys_text_report(world):
    cs = world["contigs"]
    linear = {ids.contig_uid(c.seq).hex(): c.seq for name, c in cs.items() if name != "across"}
    rep = run_minced(linear, world["base"] / "oracle")
    got = arrays_of(world["features"], world["elements"])
    for name, seq in linear.items():
        want = sorted((a["range"][0] - 1, a["range"][1], oracle_elements(a)) for a in rep.get(name, []))
        have = [(s, e, [(k, a, b, q) for k, _, a, b, q in els]) for s, e, els in got.get(bytes.fromhex(name), [])]
        assert have == want, name
    assert sum(len(v) for v in rep.values()) == 5  # plain, reverse, two (2) and the linear twin
    # the circular contig: the binary on its linear twin, its coordinates rotated as the contig was
    twin, across = cs["linear_twin"], cs["across"]
    n = len(twin.seq)
    k = (twin.arrays[0].start - across.arrays[0].start) % n
    assert across.seq == twin.seq[k:] + twin.seq[:k]
    (a,) = rep[ids.contig_uid(twin.seq).hex()]
    s0 = (a["range"][0] - 1 - k) % n
    want = (s0, s0 + a["range"][1] - a["range"][0] + 1,
            [(el, (s - k) % n, (s - k) % n + e - s, q) for el, s, e, q in oracle_elements(a)])
    (have,) = got[ids.contig_uid(across.seq)]
    assert (have[0], have[1], [(x, s, e, q) for x, _, s, e, q in have[2]]) == want


GOLDEN = os.environ.get("BIOBOX_GOLDEN_DIR")
# golden genomes with CRISPR arrays (E. coli K-12: I-E; S. pyogenes SF370: II-A; M. tuberculosis: III-A;
# S. acidocaldarius: several long arrays)
GOLDEN_IDS = ("ecoli_mg1655", "spyogenes_sf370", "mtb_h37rv", "sacidocaldarius_dsm639")


@pytest.mark.covers("O-ARRAY", "G-04")
@pytest.mark.skipif(not GOLDEN, reason="golden set not fetched (BIOBOX_GOLDEN_DIR; scripts/pod_pytest.sh --golden)")
def test_on_real_genomes_every_array_equals_the_binarys_text_report(tmp_path, record_property):
    lock = yaml.safe_load((c08.FIX.parents[2] / "golden.lock.yaml").read_text())
    gold = {g["id"]: g for g in yaml.safe_load((c08.FIX.parents[2] / "golden.yaml").read_text())["genomes"]}
    objs, genomes, records = [], [], {}
    for gid in GOLDEN_IDS:
        path = Path(GOLDEN) / lock["files"][gid]["path"]
        recs = fasta.read(path)
        uids = tuple(ids.contig_uid(r.seq) for r in recs)
        circ = tuple(u for u, spec in zip(uids, gold[gid]["records"]) if spec["topology"] == "circular")
        objs.append(SourceObject(path, "golden", "1", "INSDC (public)", True))
        genomes.append(GenomeInput(uids, gold[gid]["genetic_code"], kind=gold[gid]["kind"], circular=circ))
        records.update({u.hex(): r.seq.upper() for u, r in zip(uids, recs)})
    root = tmp_path / "root"
    ingest(objs, root, "G")
    write_genome_input(root, "G", genomes)
    rep = _crispr().run(root, "G", 2, tmp=tmp_path / "scratch")
    p = rep["producer"]
    got = arrays_of(_rows(root, f"l1/feature/producer={p}/stage=G"), _rows(root, f"l1/crispr/producer={p}/stage=G"))
    oracle = run_minced(records, tmp_path / "oracle")
    n_arrays = 0
    for name, arrs in sorted(oracle.items()):
        want = sorted((a["range"][0] - 1, a["range"][1], oracle_elements(a)) for a in arrs)
        have = [(s, e, [(k, a, b, q) for k, _, a, b, q in els]) for s, e, els in got.get(bytes.fromhex(name), [])]
        assert have == want, name
        n_arrays += len(want)
    assert set(got) <= {bytes.fromhex(k) for k in oracle}  # no array the binary did not report
    record_property("golden_arrays", n_arrays)
    print("golden arrays", n_arrays, {k: len(v) for k, v in oracle.items()})  # shown with pytest -rA
    assert n_arrays >= len(GOLDEN_IDS) and rep["arrays"] == n_arrays


# --- determinism ---------------------------------------------------------------------------------------------------


@pytest.mark.covers("D-01", "D-02", "MR-01", "MR-16")
@pytest.mark.parametrize("variant", ["again", "4 processes, one contig per batch, reverse order"])
def test_bytes_do_not_depend_on_processes_batches_order_or_the_run(world, tmp_path, variant):
    kw = {} if variant == "again" else {"genomes": world["genomes"][::-1], "batch_bases": 1}
    cpus = 1 if variant == "again" else 4
    out = tmp_path / "out"
    rep = _crispr().run(world["before"], "P", cpus, tmp=world["scratch"], out=out, **kw)
    if "batch_bases" in kw:
        assert rep["batches"] == len(world["genomes"])
    assert digests(out) == world["files"]


@pytest.mark.covers("D-01")
def test_a_contigs_arrays_do_not_depend_on_the_other_records_of_its_batch(world, tmp_path):
    """MinCED alone on each planted contig gives the GFF lines it gives for that contig in one batch of all."""
    cr = _crispr()
    java, jar = cr._tool()
    recs = {ids.contig_uid(c.seq).hex(): c.seq for c in world["contigs"].values() if not c.circular}

    def gff(names: list[str], d: Path) -> list[str]:
        d.mkdir(parents=True, exist_ok=True)
        (d / "in.fna").write_bytes(b"".join(b">" + n.encode() + b"\n" + recs[n] + b"\n" for n in names))
        subprocess.run(cr._argv(java, jar, d / "in.fna", d / "o.txt", d / "o.gff", d), check=True, capture_output=True)
        return [ln.split("\t", 8)[:8] for ln in (d / "o.gff").read_text().splitlines() if not ln.startswith("#")]

    together = gff(sorted(recs), tmp_path / "all")
    alone = [x for k, n in enumerate(sorted(recs)) for x in gff([n], tmp_path / f"one{k}")]
    assert together == alone and together


# --- run-step, refusals, ledger -------------------------------------------------------------------------------------


@pytest.mark.covers("LY-10", "D-01")
def test_under_run_step_equals_the_direct_call(world, tmp_path):
    root = tmp_path / "root"
    shutil.copytree(world["before"], root)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    step = c08.Step("P", "crispr", lambda r: ["stage", "crispr", "--root", str(r), "--stage", "P", "--cpus", "2",
                                              "--tmp", str(scratch)], None)
    prof = c08.profile_at(tmp_path / "prof")
    assert c08.run_cli(prof, root, step) == 0
    got, before = c08.tree(root), c08.tree(world["before"])
    assert {k: v for k, v in got.items() if k not in before} == world["files"]
    assert {k: got[k] for k in before} == before
    assert c08.published(root, step) == sorted(world["files"])


@pytest.mark.covers("D-03")
def test_refuses_before_any_work(world, tmp_path):
    cr = _crispr()
    out = tmp_path / "out"
    with pytest.raises(FileExistsError, match="writes its partitions once"):
        cr.run(world["root"], "P", 1, tmp=world["scratch"], out=out)
    assert not out.exists() or c08.tree(out) == {}
    for bad in (0, None, 1.5, True):
        with pytest.raises(ValueError, match="cpus"):
            cr.run(world["before"], "P", bad, tmp=world["scratch"], out=tmp_path / "o2")
    with pytest.raises(ValueError, match="scratch"):
        cr.run(world["before"], "P", 1, tmp=world["before"] / "tmp", out=tmp_path / "o3")
    assert not (tmp_path / "o2").exists() and not (tmp_path / "o3").exists()


@pytest.mark.covers("C-08", "L-01")
def test_the_ledger_names_minceds_version_and_arguments_and_balances(world):
    root, p = world["root"], world["producer"]
    step = f"crispr-{p}"
    rows = ledger.read(root, "P", step)
    assert {r["object_key"] for r in rows} == {cu for g in world["genomes"] for cu in g.contig_uids}
    lock = provenance.locked_versions()
    version = re.search(r"^minced (\S+)", subprocess.run([MINCED, "--version"], capture_output=True, text=True).stdout,
                        re.M)[1]
    assert version == lock["minced"] == "0.4.2"  # container/environment.yml's pin, read from the binary
    assert p == f"minced-{version}-w8-nr3-r23x47-s26x50"
    for r in rows:
        tv = json.loads(r["tool_versions_json"])
        assert tv["minced"] == lock["minced"] and tv["pyarrow"] == lock["pyarrow"] and tv["java"]
        assert tv["minced_args"] == "-searchWL 8 -minNR 3 -minRL 23 -maxRL 47 -minSL 26 -maxSL 50"
    assert sum(r["n_kept"] for r in rows) == len(world["features"])
    assert ledger.check_balance(root, "P", [step]) == []


# --- shared contigs (DEC-26, DEC-27) ----------------------------------------------------------------------------------


@pytest.mark.covers("C-20", "C-21a")
def test_a_contig_another_stage_or_genome_owns_is_skipped_and_ledgered(world, tmp_path):
    cr = _crispr()
    cs = world["contigs"]
    root, scratch = tmp_path / "root", world["scratch"]
    a = {"plain": cs["plain"]}
    ingest(write_objects(tmp_path / "oA", a), root, "A")
    write_genome_input(root, "A", genomes_of(a))
    b = {"plain": cs["plain"], "two": cs["two"]}
    ingest(write_objects(tmp_path / "oB", b), root, "B")
    cu = {k: ids.contig_uid(c.seq) for k, c in cs.items()}
    write_genome_input(root, "B", [GenomeInput((cu["plain"], cu["two"]), 11, kind="prokaryote")])
    res = cr.run(root, "B", 2, tmp=scratch, seq_stages=["B", "A"])
    p = res["producer"]
    assert res["skipped_contigs"] == [{"contig_uid": cu["plain"].hex(), "owner_stage": "A"}]
    feats = _rows(root, f"l1/feature/producer={p}/stage=B")
    assert {f["contig_uid"] for f in feats} == {cu["two"]}
    rows = ledger.read(root, "B", f"crispr-{p}")
    assert {r["object_key"]: r["reason_code"] for r in rows} == {cu["plain"]: "owned_by:A", cu["two"]: "ok"}
    assert ledger.check_balance(root, "B", [f"crispr-{p}"]) == []
    # within one stage: the contig of two genomes is searched once, for the genome with the least genome_uid
    c = {"reverse": cs["reverse"], "none": cs["none"], "linear_twin": cs["linear_twin"]}
    ingest(write_objects(tmp_path / "oC", c), root, "C")
    g1, g2 = (cu["reverse"], cu["none"]), (cu["reverse"], cu["linear_twin"])
    write_genome_input(root, "C", [GenomeInput(g1, 11, kind="prokaryote", circular=(cu["none"],)),
                                   GenomeInput(g2, 11, kind="prokaryote")])
    res = cr.run(root, "C", 2, tmp=scratch)
    owner = min(ids.genome_uid(list(g1)), ids.genome_uid(list(g2)))
    assert res["skipped_contigs"] == [{"contig_uid": cu["reverse"].hex(), "owner_stage": "C", "owner_genome": owner.hex()}]
    rows = ledger.read(root, "C", f"crispr-{p}")
    mine = sorted(r["reason_code"] for r in rows if r["object_key"] == cu["reverse"])
    assert mine == ["ok", f"owned_by:C/{owner.hex()}"]
    feats = _rows(root, f"l1/feature/producer={p}/stage=C")
    assert sorted(f["contig_uid"] for f in feats) == sorted([cu["reverse"], cu["linear_twin"]])
    assert ledger.check_balance(root, "C", [f"crispr-{p}"]) == []


# --- the fixture: nothing else changes (MR-20) ----------------------------------------------------------------------


@pytest.fixture(scope="module")
def fixture_stage(tmp_path_factory):
    """The micro-S0 fixture through inputs, ingest, genes, IGRs and tapes; then the producer."""
    base = tmp_path_factory.mktemp("crispr-s0")
    root, scratch = base / "root", base / "scratch"
    scratch.mkdir()
    ingest([SourceObject(FIX / o["file"], CONTRACT["source_name"], CONTRACT["source_release"], CONTRACT["license"],
                         True) for o in CONTRACT["objects"]], root, "S0")
    genomes = c08.genome_inputs("S0") + c08.genome_inputs("S1")
    write_genome_input(root, "S0", genomes)
    call_genes(root, "S0", genomes, tmp=scratch)
    igr.run(root, "S0", tmp=scratch)
    tapes.run(root, "S0", tmp=scratch)
    before = digests(root)
    rep = _crispr().run(root, "S0", 2, tmp=scratch)
    return {"root": root, "before": before, "after": digests(root), "report": rep}


@pytest.mark.covers("MR-20", "D-01")
def test_on_the_fixture_every_other_file_is_byte_identical(fixture_stage, record_property):
    before, after, p = fixture_stage["before"], fixture_stage["after"], fixture_stage["report"]["producer"]
    assert {k: after[k] for k in before} == before
    new = sorted(k for k in after if k not in before)
    assert new == sorted([f"l1/feature/producer={p}/stage=S0/part-0.parquet",
                          f"l1/crispr/producer={p}/stage=S0/part-0.parquet",
                          f"ledger/stage=S0/step=crispr-{p}/part-0.parquet"])
    record_property("fixture_arrays", fixture_stage["report"]["arrays"])
    print("fixture", fixture_stage["report"])  # shown with pytest -rA


# --- tapes ---------------------------------------------------------------------------------------------------------


@pytest.mark.covers("C-04", "MR-20")
def test_a_tape_built_with_the_producer_carries_the_arrays_and_keeps_every_gene_slot(world, tmp_path):
    root = tmp_path / "root"
    shutil.copytree(world["root"], root)
    p = world["producer"]
    tapes.run(root, "P", tmp=world["scratch"])
    tapes.run(root, "P", producers=[V, p], tmp=world["scratch"])
    one = {r["contig_uid"]: r for r in _rows(root, f"l3/tapes/tape_version={tapes.tape_version([V])}/stage=P")}
    both = {r["contig_uid"]: r for r in _rows(root, f"l3/tapes/tape_version={tapes.tape_version([V, p])}/stage=P")}
    feats = {f["feature_uid"]: f for f in world["features"]}
    seen = set()
    for cu, t in both.items():
        slots = dict(zip(t["feature_uids"], zip(t["feature_kinds"], t["protein_uids"], t["start_bps"], t["end_bps"])))
        for fu, (kind, pu, s, e) in slots.items():
            if kind == CRISPR:
                assert pu is None and (feats[fu]["start"], feats[fu]["end"]) == (s, e)
                seen.add(fu)
        old = dict(zip(one[cu]["feature_uids"], one[cu]["protein_uids"]))
        assert {fu: slots[fu][1] for fu in old} == old  # every gene slot, keyed by feature_uid, is unchanged
    assert seen == set(feats)


# --- pruning (AGENTS §4.4) -----------------------------------------------------------------------------------------

ROW_CAP = 64  # rows per row group in the pruning store (test side only), so the table has many row groups


@pytest.fixture(scope="module")
def many(tmp_path_factory):
    """48 planted contigs of two arrays each, the producer's files written with ROW_CAP rows per row group."""
    base = tmp_path_factory.mktemp("crispr-prune")
    root, scratch = base / "root", base / "scratch"
    scratch.mkdir()
    cs = {f"c{i:02d}": plant.plant(100 + i, 3000, [(800, plant.REPEAT_IF, 5, False), (2100, plant.REPEAT_IIA, 4, i % 2 == 1)])
          for i in range(48)}
    ingest(write_objects(base / "o", cs), root, "M")
    write_genome_input(root, "M", genomes_of(cs))
    real = io.write_table

    def capped(batches, spec, path, *, row_group_rows=io.ROW_GROUP_ROWS, **kw):
        return real(batches, spec, path, row_group_rows=min(row_group_rows, ROW_CAP), **kw)

    mp = pytest.MonkeyPatch()
    mp.setattr(io, "write_table", capped)
    try:
        rep = _crispr().run(root, "M", 2, tmp=scratch)
    finally:
        mp.undo()
    p = rep["producer"]
    files = io.part_files(root / f"l1/crispr/producer={p}/stage=M")
    con = duckdb.connect()
    con.execute(f"CREATE VIEW l1_crispr AS SELECT * FROM read_parquet([{', '.join(repr(str(f)) for f in files)}])")
    return {"con": con, "files": files, "rows": _rows(root, f"l1/crispr/producer={p}/stage=M"), "base": base}


def _scanned(prof: dict) -> list[tuple[int, int]]:
    return [(s["row_groups_scanned"], s["row_groups_total"]) for s in prof["scans"] if s["files"]]


@pytest.mark.covers("LY-06")
def test_a_lookup_by_contig_and_by_contig_and_feature_uid_reads_at_most_two_row_groups(many, record_property):
    rows, con = many["rows"], many["con"]
    total = sum(pq.ParquetFile(f).metadata.num_row_groups for f in many["files"])
    assert total >= 12, total
    cu = sorted({r["contig_uid"] for r in rows})[17]
    fu = next(r["feature_uid"] for r in rows if r["contig_uid"] == cu)
    by_contig = catalog.profile_query(con, "SELECT count(*), sum(ordinal) FROM l1_crispr WHERE contig_uid = ?", [cu])
    by_feature = catalog.profile_query(
        con, "SELECT count(*), sum(ordinal) FROM l1_crispr WHERE contig_uid = ? AND feature_uid = ?", [cu, fu])
    bare = catalog.profile_query(con, "SELECT count(*), sum(ordinal) FROM l1_crispr WHERE feature_uid = ?", [fu])
    seen = json.dumps({"total": total, "contig": _scanned(by_contig), "contig_feature": _scanned(by_feature),
                       "feature_only": _scanned(bare)})
    record_property("prune", seen)
    print("prune", seen)  # shown with pytest -rA
    assert by_contig["result"][0][0] == sum(r["contig_uid"] == cu for r in rows)
    assert by_feature["result"][0][0] == sum(r["feature_uid"] == fu for r in rows) == bare["result"][0][0]
    for prof in (by_contig, by_feature):
        (s,) = _scanned(prof)
        assert 1 <= s[0] <= 2 and s[1] == total, _scanned(prof)
