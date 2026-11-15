"""Ring R2, card E-15 (DEC-40): the producer of two-repeat CRISPR arrays beside Cas genes
(`biobox.l1.crispr.run_near_cas`) on constructed contigs.

The truth is the constructor's (AGENTS §5.2): arrays planted by `tests/r1/test_crispr_plant.py` (repeats and spacers
known as the sequence is made), gene calls written by hand as the stage's `feature_protein` rows, and the Pfam scan
written by hand as its `l1_domains` rows (library `pfam-38.2`), so every gap and containment is known by construction.
The default producer (MinCED at 3 repeats, card E-12) runs for real first; then the filtered one. A planted array is
expected among the kept calls exactly when the test's own reading of the rule (`expected`, written here: loops over the
planted genes and the circle's shifts) keeps it:

- a two-repeat array 128 bases from a Cas gene (kept); gaps of exactly 1,000 bases up- and downstream (kept) and of
  1,001 (dropped); 5 kb from the nearest Cas gene (dropped);
- wholly inside a gene, and exactly the span of a gene (dropped); overlapping part of a Cas gene (kept);
- a three-repeat array the default producer calls, beside a Cas gene (dropped: it overlaps the default's array);
- across the origin of a circular contig beside a Cas gene (kept, found once, end past the contig); just before the
  origin with its Cas gene just after it (kept: the gap is measured round the circle), and the same layout on a linear
  contig (dropped);
- beside a gene whose protein has no domain row, and beside one with a domain of a family not on the list (dropped);
- a stage whose domain scan has no row at all (nothing kept, no crash).

Then: the default producer's files are never rewritten (the root's bytes before and after); the bytes do not depend on
the number of MinCED processes, the batches, the genome order or the run; under `bio-box run-step` (`stage crispr
--near-cas`) the files equal the direct call's; refusals before any work; the ledger balances, names every reason and
passes the C-08 check; the lookups of genes, default arrays and domains by contig prune (AGENTS §4.4).
"""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import fasta, ids, io, ledger, provenance
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.registry.schemas import conformance_problems, table

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "r1"))
import test_builders_out as c08  # noqa: E402
import test_crispr_plant as plant  # noqa: E402

MINCED = shutil.which("minced")
pytestmark = [pytest.mark.r2, pytest.mark.covers("E-15", "O-ARRAY"),
              pytest.mark.skipif(MINCED is None, reason="minced (MinCED, a Java tool) not on PATH; it is in the "
                                                        "bio-box image, not in the uv environment of GitHub Actions")]

CDS = ids.KIND["cds"]
LIB_VERSION = "38.2"
LIB_DIR = "l1/domains/library=pfam-38.2"
CAS6F, CAS7F, CAS1, TNIQ = "PF09618.14", "PF09615.12", "PF01867.19", "PF06527.15"  # model_acc with versions
BP = 1000


def _crispr():
    from biobox.l1 import crispr

    return crispr


# --- the constructed stage ------------------------------------------------------------------------------------------


def two_repeat(seed: int, length: int, pos: int, repeat: bytes = plant.REPEAT_IF, end_gap: int | None = None
               ) -> plant.Contig:
    """A contig with one planted array of two repeats; the first seed from `seed` whose bases beside the two copies
    differ (before both, and after both), so CRT cannot extend the repeat (card E-13's construction). `end_gap`: put
    the array that many bases before the contig's end instead of at `pos`."""
    for s in range(seed, seed + 200):
        at = length - end_gap if end_gap is not None else pos
        c = plant.plant(s, length, [(at, repeat, 2, False)])
        r1, _, r2 = c.arrays[0].elements
        after = c.seq[r1.end:r1.end + 1] != c.seq[r2.end:r2.end + 1] if r2.end < len(c.seq) else True
        if c.seq[r1.start - 1] != c.seq[r2.start - 1] and after:
            return c
    raise AssertionError(f"no seed in {seed}..{seed + 199} gives distinct flanking bases")


class Gene:
    def __init__(self, start: int, end: int, strand: int = 1, domain: str | None = None):
        self.start, self.end, self.strand, self.domain = start, end, strand, domain


def cases() -> dict[str, tuple[plant.Contig, list[Gene], bool]]:
    """name -> (contig, its genes, whether the rule keeps its planted array). The genes are placed from the planted
    array's coordinates; `False` cases say why in their name."""
    out: dict[str, tuple[plant.Contig, list[Gene], bool]] = {}

    def add(name: str, c: plant.Contig, genes, keep: bool) -> None:
        out[name] = (c, genes, keep)

    c = two_repeat(101, 8000, 3000)
    a = c.arrays[0]
    add("near", c, [Gene(100, 1000), Gene(a.end + 128, a.end + 728, -1, CAS6F)], True)
    c = two_repeat(201, 8000, 3000)
    a = c.arrays[0]
    add("gap_1000_down", c, [Gene(a.end + BP, a.end + BP + 600, 1, CAS7F)], True)
    c = two_repeat(301, 8000, 3000)
    a = c.arrays[0]
    add("gap_1000_up", c, [Gene(a.start - BP - 600, a.start - BP, -1, CAS1)], True)
    c = two_repeat(401, 8000, 3000)
    a = c.arrays[0]
    add("far_gap_1001", c, [Gene(a.start - BP - 601, a.start - BP - 1, 1, CAS6F),
                            Gene(a.end + BP + 1, a.end + BP + 601, 1, CAS7F)], False)
    c = two_repeat(501, 12000, 3000)
    a = c.arrays[0]
    add("far_5kb", c, [Gene(a.end + 5000, a.end + 5600, 1, CAS6F), Gene(a.end + 100, a.end + 400)], False)
    c = two_repeat(601, 8000, 3000)
    a = c.arrays[0]
    add("inside_gene", c, [Gene(a.start - 40, a.end + 40), Gene(a.end + 100, a.end + 700, 1, CAS6F)], False)
    c = two_repeat(701, 8000, 3000)
    a = c.arrays[0]
    add("inside_exactly", c, [Gene(a.start, a.end, -1), Gene(a.end + 100, a.end + 700, 1, CAS6F)], False)
    c = two_repeat(801, 8000, 3000)
    a = c.arrays[0]
    add("overlaps_cas_gene", c, [Gene(a.start + 10, a.end + 500, 1, CAS6F)], True)
    c = plant.plant(901, 8000, [(3000, plant.REPEAT_IIA, 3, False)])
    a = c.arrays[0]
    add("three_repeats_default", c, [Gene(a.end + 128, a.end + 728, 1, CAS6F)], False)
    # across the origin: a linear twin rotated so the origin falls inside the spacer; its Cas gene after the array
    twin = two_repeat(1001, 7000, 3000)
    a = twin.arrays[0]
    k = (a.elements[1].start + a.elements[1].end) // 2
    across = plant.rotate(twin, k)
    n = len(across.seq)
    add("across_origin", across, [Gene(a.end + 128 - k, a.end + 728 - k, 1, CAS6F), Gene(2000, 2900)], True)
    # just before the origin, the Cas gene just after it: 80 bases round the circle
    c = two_repeat(1101, 6000, 0, plant.REPEAT_IE, end_gap=50)
    c.circular = True
    add("gap_round_origin", c, [Gene(30, 630, 1, CAS6F), Gene(2000, 2600)], True)
    c = two_repeat(1201, 6000, 0, plant.REPEAT_IE, end_gap=50)
    add("gap_round_origin_but_linear", c, [Gene(30, 630, 1, CAS6F), Gene(2000, 2600)], False)
    c = two_repeat(1301, 8000, 3000)
    a = c.arrays[0]
    add("no_domain_row", c, [Gene(a.end + 128, a.end + 728)], False)
    c = two_repeat(1401, 8000, 3000)
    a = c.arrays[0]
    add("not_a_cas_family", c, [Gene(a.end + 128, a.end + 728, 1, TNIQ)], False)
    assert across.arrays[0].start < n < across.arrays[0].end
    return out


def gap(a: tuple[int, int], g: tuple[int, int]) -> int:
    return max(0, g[0] - a[1], a[0] - g[1])


def expected(c: plant.Contig, genes: list[Gene], default: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """The rule read independently of the producer: per planted two-repeat array, every gene and default array at
    every shift of the circle (0, and -n, +n on a circular contig)."""
    n = len(c.seq)
    shifts = (-n, 0, n) if c.circular else (0,)
    keep = []
    for arr in c.arrays:
        s, e = arr.start, arr.end
        if any(ds + k < e and s < de + k for ds, de in default for k in shifts):
            continue
        if not any(g.domain in (CAS6F, CAS7F, CAS1) and gap((s, e), (g.start + k, g.end + k)) <= BP
                   for g in genes for k in shifts):
            continue
        if any(g.start + k <= s and g.end + k >= e for g in genes for k in shifts):
            continue
        keep.append((s, e))
    return keep


def write_objects(d: Path, contigs: dict[str, plant.Contig]) -> list[SourceObject]:
    d.mkdir(parents=True, exist_ok=True)
    objs = []
    for name, c in contigs.items():
        path = d / f"{name}.fna"
        path.write_bytes(b"\n".join(fasta.format_record(f"e15_{name}", c.seq, 80)) + b"\n")
        objs.append(SourceObject(path, "planted", "1", "INSDC (public)", True))
    return objs


def genomes_of(contigs: dict[str, plant.Contig]) -> list[GenomeInput]:
    out = []
    for c in contigs.values():
        cu = ids.contig_uid(c.seq)
        out.append(GenomeInput((cu,), 11, kind="plasmid", circular=(cu,) if c.circular else ()))
    return out


def gene_rows(contigs: dict[str, tuple[plant.Contig, list[Gene]]]) -> tuple[list[dict], list[dict]]:
    """(`feature_protein` rows, `l1_domains` rows) of the constructed genes."""
    fp, dom = [], []
    for name, (c, genes) in sorted(contigs.items()):
        cu = ids.contig_uid(c.seq)
        for k, g in enumerate(genes):
            pu = ids.protein_uid(f"M{name}{k}")
            fp.append({"contig_uid": cu, "feature_uid": ids.feature_uid(cu, g.start, g.end, g.strand, CDS),
                       "start": g.start, "end": g.end, "strand": g.strand, "protein_uid": pu, "genetic_code": 11,
                       "edit_type": "none"})
            if g.domain is not None:
                dom.append({"protein_uid": pu, "model_acc": g.domain, "model_name": "m", "model_len": 100,
                            "seq_bitscore": 50.0, "seq_evalue_zref": 1e-10, "dom_bitscore": 50.0,
                            "dom_ievalue_zref": 1e-10, "env_start": 0, "env_end": 90, "ali_start": 1,
                            "ali_end": 89, "hmm_from": 0, "hmm_to": 99, "cutoff": "gathering", "z": 1e8,
                            "dom_z": 1e8})
    fp.sort(key=lambda r: (r["contig_uid"], r["start"], r["end"]))
    dom.sort(key=lambda r: (r["protein_uid"], r["env_start"], r["model_acc"]))
    return fp, dom


def write_rows(root: Path, name: str, rel: str, rows: list[dict], row_group_rows: int = io.ROW_GROUP_ROWS) -> None:
    t = pa.Table.from_pylist(rows, schema=table(name).schema)
    assert conformance_problems(t, table(name)) == []
    io.write_table([t], table(name), root / rel / "part-0.parquet", row_group_rows=row_group_rows)


def make_stage(base: Path, stage: str, layout: dict[str, tuple[plant.Contig, list[Gene], bool]],
               domains: bool = True) -> Path:
    """A root with the stage's contigs ingested, its genome_input, the hand-made gene calls and Pfam scan (with no
    rows when `domains` is False) and the default producer's files."""
    root, scratch = base / "root", base / "scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    contigs = {k: v[0] for k, v in layout.items()}
    ingest(write_objects(base / "objects", contigs), root, stage)
    write_genome_input(root, stage, genomes_of(contigs))
    fp, dom = gene_rows({k: (v[0], v[1]) for k, v in layout.items()})
    write_rows(root, "feature_protein", f"l1/feature_protein/gene_call_version={V}/stage={stage}", fp)
    write_rows(root, "l1_domains", f"{LIB_DIR}/stage={stage}", dom if domains else [])
    _crispr().run(root, stage, 1, tmp=scratch)  # the default producer, for real
    return root


def digests(root: Path) -> dict[str, str]:
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(Path(root).rglob("*")) if p.is_file()}


def _rows(root: Path, rel: str) -> list[dict]:
    return [r for f in io.part_files(root / rel) for r in pq.read_table(f, partitioning=None).to_pylist()]


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    base = tmp_path_factory.mktemp("e15")
    layout = cases()
    root = make_stage(base, "Q", layout)
    before = digests(root)
    out = base / "out"
    cr = _crispr()
    report: dict = {}
    res = cr.run_near_cas(root, "Q", 1, library_version=LIB_VERSION, tmp=base / "scratch", out=out, report=report)
    p = res["producer"]
    return {"base": base, "root": root, "before": before, "out": out, "res": res, "producer": p, "layout": layout,
            "report": report, "files": digests(out),
            "features": _rows(out, f"l1/feature/producer={p}/stage=Q"),
            "elements": _rows(out, f"l1/crispr/producer={p}/stage=Q")}


def default_arrays(root: Path, stage: str) -> dict[bytes, list[tuple[int, int]]]:
    out: dict[bytes, list[tuple[int, int]]] = {}
    for r in _rows(root, f"l1/feature/producer={_crispr().producer()}/stage={stage}"):
        out.setdefault(r["contig_uid"], []).append((r["start"], r["end"]))
    return out


# --- the planted truth ----------------------------------------------------------------------------------------------


def test_the_kept_calls_are_exactly_the_planted_arrays_the_rule_keeps(world):
    d3 = default_arrays(world["root"], "Q")
    want: dict[bytes, list[tuple]] = {}
    keep_names = set()
    for name, (c, genes, keep) in world["layout"].items():
        cu = ids.contig_uid(c.seq)
        mine = expected(c, genes, d3.get(cu, []))
        assert bool(mine) == keep, name  # the test's own reading agrees with the case's label
        if keep:
            keep_names.add(name)
            want[cu] = [(a.start, a.end, [(e.kind, e.ordinal, e.start, e.end, e.seq.decode()) for e in a.elements])
                        for a in c.arrays]
    els: dict[bytes, list[tuple]] = {}
    for r in world["elements"]:
        els.setdefault(r["feature_uid"], []).append((r["element"], r["ordinal"], r["start"], r["end"], r["seq"]))
    got: dict[bytes, list[tuple]] = {}
    for f in world["features"]:
        got.setdefault(f["contig_uid"], []).append((f["start"], f["end"], sorted(els[f["feature_uid"]],
                                                                                key=lambda t: t[1])))
    assert got == want
    assert len(keep_names) == 6
    # the array across the origin is found once, its end past the contig (AGENTS §3.1)
    across = world["layout"]["across_origin"][0]
    (arr,) = got[ids.contig_uid(across.seq)]
    assert arr[0] < len(across.seq) < arr[1]
    # the default producer called the three-repeat array, so that case is the overlap rule's
    three = world["layout"]["three_repeats_default"][0]
    assert d3[ids.contig_uid(three.seq)] == [(three.arrays[0].start, three.arrays[0].end)]


@pytest.mark.covers("C-01", "C-02", "LY-01")
def test_rows_are_keyed_conform_and_name_the_producer(world):
    out, p = world["out"], world["producer"]
    cr = _crispr()
    assert p == cr.near_cas_producer(LIB_VERSION) == \
        f"minced-{cr.minced_version()}-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2"
    for name, rel in (("l1_crispr", f"l1/crispr/producer={p}/stage=Q"), ("feature", f"l1/feature/producer={p}/stage=Q")):
        (f,) = io.part_files(out / rel)
        t = pq.read_table(f, partitioning=None)
        assert t.num_rows > 0 and conformance_problems(t, table(name)) == []
    for f in world["features"]:
        assert f["kind"] == ids.KIND["crispr_array"] and f["strand"] == 1
        assert f["feature_uid"] == ids.feature_uid(f["contig_uid"], f["start"], f["end"], 1, f["kind"])
    assert sorted(world["files"]) == sorted([f"l1/feature/producer={p}/stage=Q/part-0.parquet",
                                             f"l1/crispr/producer={p}/stage=Q/part-0.parquet",
                                             f"ledger/stage=Q/step=crispr-{p}/part-0.parquet"])


def test_a_stage_whose_domain_scan_has_no_row_keeps_nothing(tmp_path):
    layout = {k: v for k, v in cases().items() if k in ("near", "across_origin", "gap_1000_up")}
    root = make_stage(tmp_path, "E", layout, domains=False)
    out = tmp_path / "out"
    res = _crispr().run_near_cas(root, "E", 1, library_version=LIB_VERSION, tmp=tmp_path / "scratch", out=out)
    p = res["producer"]
    assert res["kept"] == 0 and res["arrays"] == 3
    assert _rows(out, f"l1/feature/producer={p}/stage=E") == [] and _rows(out, f"l1/crispr/producer={p}/stage=E") == []
    rows = ledger.read(out, "E", f"crispr-{p}")
    assert {r["reason_code"] for r in rows} == {"ok", "far_from_cas"}
    assert ledger.check_balance(out, "E", [f"crispr-{p}"]) == []


# --- the default producer's bytes, determinism, run-step -----------------------------------------------------------


@pytest.mark.covers("MR-20", "D-01")
def test_the_default_producers_files_are_never_rewritten(world):
    assert digests(world["root"]) == world["before"]  # the root, its default files included, as before
    assert not set(world["files"]) & set(world["before"])


@pytest.mark.covers("D-01", "D-02", "MR-01", "MR-16")
@pytest.mark.parametrize("variant", ["again", "4 processes, one contig per batch, reverse order"])
def test_bytes_do_not_depend_on_processes_batches_order_or_the_run(world, tmp_path, variant):
    kw = {} if variant == "again" else {"genomes": genomes_of({k: v[0] for k, v in world["layout"].items()})[::-1],
                                        "batch_bases": 1}
    cpus = 1 if variant == "again" else 4
    out = tmp_path / "out"
    _crispr().run_near_cas(world["root"], "Q", cpus, library_version=LIB_VERSION, tmp=world["base"] / "scratch",
                           out=out, **kw)
    assert digests(out) == world["files"]


@pytest.mark.covers("LY-10", "D-01")
def test_under_run_step_equals_the_direct_call(world, tmp_path):
    root = tmp_path / "root"
    shutil.copytree(world["root"], root)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    step = c08.Step("Q", "crispr-near-cas", lambda r: ["stage", "crispr", "--root", str(r), "--stage", "Q", "--near-cas",
                                                       "--library", "pfam", "--library-version", LIB_VERSION,
                                                       "--cpus", "2", "--tmp", str(scratch)], None)
    prof = c08.profile_at(tmp_path / "prof")
    assert c08.run_cli(prof, root, step) == 0
    got, before = c08.tree(root), c08.tree(world["root"])
    assert {k: got[k] for k in before} == before
    assert {k: v for k, v in got.items() if k not in before} == world["files"]


@pytest.mark.covers("D-03")
def test_refuses_before_any_work(world, tmp_path):
    cr = _crispr()
    scratch = world["base"] / "scratch"
    with pytest.raises(FileExistsError, match="writes its partitions once"):
        cr.run_near_cas(world["root"], "Q", 1, library_version=LIB_VERSION, tmp=scratch, out=world["out"])
    for bad in ("", "38-2", None):
        with pytest.raises(ValueError, match="library version"):
            cr.run_near_cas(world["root"], "Q", 1, library_version=bad, tmp=scratch, out=tmp_path / "o1")
    with pytest.raises(ValueError, match="cpus"):
        cr.run_near_cas(world["root"], "Q", 0, library_version=LIB_VERSION, tmp=scratch, out=tmp_path / "o1")
    # each input the filter reads must be there: the default producer's arrays, the gene calls, the domain scan
    for missing, words in ((f"l1/feature/producer={cr.producer()}/stage=Q", "default producer"),
                           (f"l1/feature_protein/gene_call_version={V}/stage=Q", "gene calls"),
                           (f"{LIB_DIR}/stage=Q", "domain scan")):
        root = tmp_path / f"r-{len(missing)}"
        shutil.copytree(world["root"], root)
        shutil.rmtree(root / missing)
        with pytest.raises(FileNotFoundError, match=words):
            cr.run_near_cas(root, "Q", 1, library_version=LIB_VERSION, tmp=scratch, out=tmp_path / "o2")
    assert not (tmp_path / "o1").exists() and not (tmp_path / "o2").exists()


# --- the ledger -----------------------------------------------------------------------------------------------------


@pytest.mark.covers("C-08", "L-01")
def test_the_ledger_names_every_reason_its_inputs_and_balances(world):
    out, p = world["out"], world["producer"]
    cr = _crispr()
    step = f"crispr-{p}"
    rows = ledger.read(out, "Q", step)
    by: dict[bytes, dict[str, dict]] = {}
    for r in rows:
        by.setdefault(r["object_key"], {})[r["reason_code"]] = r
    want = {"three_repeats_default": "overlaps_default", "far_gap_1001": "far_from_cas", "far_5kb": "far_from_cas",
            "gap_round_origin_but_linear": "far_from_cas", "no_domain_row": "far_from_cas",
            "not_a_cas_family": "far_from_cas", "inside_gene": "inside_gene", "inside_exactly": "inside_gene"}
    for name, (c, _, keep) in world["layout"].items():
        mine = by[ids.contig_uid(c.seq)]
        assert mine["ok"]["n_in"] == 1 and mine["ok"]["n_kept"] == (1 if keep else 0), name
        assert set(mine) == ({"ok"} if keep else {"ok", want[name]}), (name, set(mine))
    assert ledger.check_balance(out, "Q", [step]) == []
    lock = provenance.locked_versions()
    for r in rows:
        tv = json.loads(r["tool_versions_json"])
        assert tv["minced"] == lock["minced"] and tv["pyarrow"] == lock["pyarrow"] and tv["java"]
        assert tv["minced_args"] == "-searchWL 8 -minNR 2 -minRL 23 -maxRL 47 -minSL 26 -maxSL 50"
        assert tv["gene_call_version"] == V and tv["domain_library"] == "pfam-38.2"
        assert tv["cas_families"] == "cas-pfam-1" and tv["cas_filter"] == "cas1kb"
        assert tv["overlap_producer"] == cr.producer()
    files = [{"path": f"l1/feature/producer={p}/stage=Q/part-0.parquet", "table": "feature",
              "partitions": {"producer": p, "stage": "Q"}},
             {"path": f"l1/crispr/producer={p}/stage=Q/part-0.parquet", "table": "l1_crispr",
              "partitions": {"producer": p, "stage": "Q"}},
             {"path": f"ledger/stage=Q/step={step}/part-0.parquet", "table": "ledger",
              "partitions": {"stage": "Q", "step": step}}]
    assert provenance.check_provenance({"root": str(out), "release_id": "rX", "stages": ["Q"], "files": files}) == []
    res = world["res"]
    assert res["kept"] == 6 and res["arrays"] == len(world["layout"])
    assert res["rejected_by_reason"] == {"exceeds_circle": 0, "far_from_cas": 5, "inside_gene": 2,
                                         "overlaps_default": 1}


# --- pruning (AGENTS §4.4) -----------------------------------------------------------------------------------------

ROW_CAP = 8  # rows per row group of the pruning store's input tables


@pytest.fixture(scope="module")
def many(tmp_path_factory):
    """48 contigs, ten genes each (one a Cas gene), all in the stage's gene calls and Pfam scan with ROW_CAP rows per
    row group; four carry a two-repeat array beside their Cas gene, four a three-repeat array the default producer
    calls (its file written with 1 row per row group; MinCED at 2 calls it too), the other 40 none."""
    base = tmp_path_factory.mktemp("e15-prune")
    root, scratch = base / "root", base / "scratch"
    scratch.mkdir()
    layout = {}
    for i in range(48):
        if i % 12 == 0:
            c = two_repeat(5000 + 100 * i, 9000, 3000)
        elif i % 12 == 6:
            c = plant.plant(9000 + i, 9000, [(3000, plant.REPEAT_IIA, 3, False)])
        else:
            c = plant.plant(9000 + i, 9000, [])
        a_end = c.arrays[0].end if c.arrays else 3100
        genes = [Gene(a_end + 128, a_end + 728, 1, CAS6F)] + [Gene(5000 + 300 * k, 5200 + 300 * k) for k in range(9)]
        layout[f"m{i:02d}"] = (c, genes, False)
    contigs = {k: v[0] for k, v in layout.items()}
    ingest(write_objects(base / "objects", contigs), root, "M")
    write_genome_input(root, "M", genomes_of(contigs))
    fp, dom = gene_rows({k: (v[0], v[1]) for k, v in layout.items()})
    write_rows(root, "feature_protein", f"l1/feature_protein/gene_call_version={V}/stage=M", fp, ROW_CAP)
    # every gene has a domain row here (the Cas genes PF09618, the others a family not on the list), so the domain
    # table has as many row groups as the gene table
    dom = sorted(dom + [{**dom[0], "protein_uid": r["protein_uid"], "model_acc": TNIQ} for r in fp
                        if r["protein_uid"] not in {d["protein_uid"] for d in dom}],
                 key=lambda r: (r["protein_uid"], r["env_start"], r["model_acc"]))
    write_rows(root, "l1_domains", f"{LIB_DIR}/stage=M", dom, ROW_CAP)
    real = io.write_table

    def capped(batches, spec, path, *, row_group_rows=io.ROW_GROUP_ROWS, **kw):
        return real(batches, spec, path, row_group_rows=1, **kw)

    mp = pytest.MonkeyPatch()
    mp.setattr(io, "write_table", capped)
    try:
        _crispr().run(root, "M", 1, tmp=scratch)
    finally:
        mp.undo()
    report: dict = {}
    res = _crispr().run_near_cas(root, "M", 1, library_version=LIB_VERSION, tmp=scratch, out=base / "out",
                                 report=report)
    return {"root": root, "layout": layout, "res": res, "report": report}


def _holding(files: list[Path], col: str, keys: set[bytes]) -> tuple[int, int]:
    """(row groups whose footer range of `col` holds one of `keys`, row groups in all)."""
    hold = total = 0
    for f in files:
        md = pq.ParquetFile(f).metadata
        i = next(k for k in range(md.num_columns) if md.schema.column(k).path == col)
        for g in range(md.num_row_groups):
            st = md.row_group(g).column(i).statistics
            total += 1
            hold += any(bytes(st.min) <= k <= bytes(st.max) for k in keys)
    return hold, total


@pytest.mark.covers("LY-06")
def test_the_filters_lookups_by_contig_read_the_row_groups_of_those_contigs(many, record_property):
    root, cr = many["root"], _crispr()
    # the contigs with a call at 2 repeats: the planted two- and three-repeat arrays
    cand = {ids.contig_uid(c.seq) for name, (c, _, _) in many["layout"].items() if int(name[1:]) % 6 == 0}
    near_proteins = {ids.protein_uid(f"M{name}0") for name, (c, _, _) in many["layout"].items()
                     if ids.contig_uid(c.seq) in cand}  # the Cas gene of each, gene 0: the only gene within 1 kb
    assert many["res"]["kept"] == 4 and many["res"]["arrays"] == 8 and len(cand) == 8
    scans = {s["table"]: s for s in many["report"]["scans"]}
    assert set(scans) == {"feature_protein", "feature", "l1_domains"}
    want = {"feature_protein": _holding(io.part_files(root / f"l1/feature_protein/gene_call_version={V}/stage=M"),
                                        "contig_uid", cand),
            "feature": _holding(io.part_files(root / f"l1/feature/producer={cr.producer()}/stage=M"), "contig_uid",
                                cand),
            "l1_domains": _holding(io.part_files(root / f"{LIB_DIR}/stage=M"), "protein_uid", near_proteins)}
    seen = {t: (s["row_groups_scanned"], s["row_groups_total"]) for t, s in scans.items()}
    record_property("e15_prune", json.dumps({"seen": seen, "holding": want}))
    print("e15 prune", seen, want)  # shown with pytest -rA
    for t, (hold, total) in want.items():
        assert seen[t] == (hold, total), (t, seen[t], hold)
    # the gene calls and the domain scan have 60 row groups each; the lookups read a fraction of them (a contig's ten
    # gene rows lie in at most 3 row groups of 8; each looked-up protein in one). Every default array is also a call
    # at 2, so its lookup holds every row group of its file (4 here)
    assert want["feature_protein"][1] == 60 and want["feature_protein"][0] <= 3 * len(cand) < 60
    assert want["l1_domains"][1] == 60 and want["l1_domains"][0] <= len(near_proteins) < 60
    assert want["feature"] == (4, 4)
