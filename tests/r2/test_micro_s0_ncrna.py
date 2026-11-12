"""Ring R2: the ncRNA producer (card E-01, operation O-SCAN; `biobox.l1.ncrna`) on the micro-S0 fixture.

The CM library is `tests/fixtures/micro_s0/rfam15.1_micro_s0.cm.gz`: nine Rfam 15.1 families that occur in
the fixture genomes (RF00001 5S_rRNA, RF00005 tRNA, RF00011 RNaseP_bact_b, RF00023 tmRNA, RF00028
Intron_gpI, RF00169 Bacteria_small_SRP, RF01804 Lambda_thermo, RF01852 tRNA-Sec, RF02111 IS009), fetched
verbatim with `cmfetch -f Rfam.cm keys.txt` (Infernal 1.1.5, in the test image) from Rfam 15.1's `Rfam.cm.gz`
(sha256 e2e636e9..., as refdb.yaml lists it) and compressed with `gzip -n -9`; the uncompressed file's
sha256 is d9b4b15c753c12cdecfbbb40bae6527ab9c7d67bb5f155c9605eb2ac1e808f83. They were chosen from a
full-Rfam `cmsearch --cut_ga -Z 100000` over the three contigs (planning/measurements/2026-11-10_ncrna_rate.md),
leaving out the rRNA subunit models (thousands of columns, slow). tRNA and tRNA-Sec hit overlapping loci of
M. genitalium, so overlapping hits of two models are exercised.

Card E-07 migration: the producer now scans a model set, its intergenic models on intergenic regions and its
whole-genome models on whole contigs (DEC-28). These tests give it the set `ncrna_set_whole.tsv`, the nine models
above all scanned on whole contigs, so it scans as card E-01 did and the oracle below is unchanged; the producer's
name now carries the set and the region rule. Regions, the set and shards are tested in `test_ncrna_regions.py`,
`test_ncrna_shards.py` and `test_ncrna_build_stage.py`.

The oracle (AGENTS §5.2) is the `cmsearch` binary itself, run here on the fixture contigs read from the
source FASTA, its `--tblout` parsed by a parser of this file. The builder's rows, after its coordinate
conversion, must equal it. A tRNA of M. genitalium is planted across the origin of two circular contigs (one
longer, one shorter than the library's longest window), which must each report it once with
`end > contig_len`. Output bytes must not depend on the thread count, the genome order or the batches; a run
under `bio-box run-step` must equal the direct call (helpers of `test_builders_out` by import); IGRs and tapes
with both producers must build with the ncRNA slots on the tape (kind ncrna, null protein); and adding the
producer must change no row keyed by `feature_uid` in the gene-call tables (MR-20 with a real producer).
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pyarrow.parquet as pq
import pytest
import yaml

from biobox import fasta, ids, io, ledger, provenance
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, read_genome_input, write_genome_input
from biobox.l1 import igr
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l1.genecall import run as call_genes
from biobox.l3 import tapes
from biobox.registry.schemas import conformance_problems, table

sys.path.insert(0, str(Path(__file__).parent))
import test_builders_out as c08  # noqa: E402

CMSEARCH = shutil.which("cmsearch")
pytestmark = [pytest.mark.r2, pytest.mark.covers("E-01"),
              pytest.mark.skipif(CMSEARCH is None, reason="cmsearch (Infernal) not on PATH; it is in the test image")]

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
LIB = FIX / "rfam15.1_micro_s0.cm.gz"
SET = FIX / "ncrna_set_whole.tsv"  # card E-07: every model of LIB on whole contigs, as card E-01 scanned
NCRNA = ids.KIND["ncrna"]
Z = "100000"  # DEC-04: Z_nt = 1e5 Mb
FAMILIES = {"RF00001", "RF00005", "RF00011", "RF00023", "RF00028", "RF00169", "RF01804", "RF01852", "RF02111"}


def _ncrna():
    from biobox.l1 import ncrna  # imported here so that a missing module fails the tests, not the collection

    return ncrna


# --- an independent reading of the binary's output ----------------------------------------------------


def oracle_tblout(path: Path) -> list[dict]:
    """`cmsearch --tblout` (format 1) rows: the 17 whitespace-separated fields before the description."""
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split()
        out.append({"target": f[0], "name": f[2], "acc": f[3], "mdl": f[4], "mdl_from": int(f[5]),
                    "mdl_to": int(f[6]), "from": int(f[7]), "to": int(f[8]), "strand": f[9], "trunc": f[10],
                    "gc": float(f[12]), "score": float(f[14]), "evalue": float(f[15]), "inc": f[16]})
    return out


def oracle_ga(cm_text: str) -> dict[str, float]:
    """GA of every model: the `GA` line of each CM header (between `ACC` and the `CM` line)."""
    out, acc = {}, None
    for line in cm_text.splitlines():
        if line.startswith("ACC "):
            acc = line.split()[1]
        elif line.startswith("GA ") and acc is not None:
            out.setdefault(acc, float(line.split()[1]))
        elif line.startswith("HMMER3/"):
            acc = None
    return out


def run_cmsearch(cm: Path, fna: Path, tbl: Path, cpus: int = 2) -> list[dict]:
    subprocess.run([CMSEARCH, "--cpu", str(cpus), "-Z", Z, "--cut_ga", "--noali", "--tblout", str(tbl),
                    "-o", str(tbl.with_suffix(".out")), str(cm), str(fna)], check=True)
    return oracle_tblout(tbl)


def converted(h: dict, cu: bytes) -> tuple:
    """The builder's row key and values for one oracle hit: 1-based inclusive (from > to on the minus strand)
    to 0-based half-open."""
    lo, hi = sorted((h["from"], h["to"]))
    key = (cu, lo - 1, hi, 1 if h["strand"] == "+" else -1, h["acc"])
    return key, (h["name"], h["mdl_from"] - 1, h["mdl_to"], h["mdl"], round(h["score"], 1), h["evalue"],
                 "5'" in h["trunc"], "3'" in h["trunc"], round(h["gc"], 2))


def row_key(r: dict) -> tuple:
    return (r["contig_uid"], r["start"], r["end"], r["strand"], r["model_acc"]), (
        r["model_name"], r["mdl_start"], r["mdl_end"], r["mdl_kind"], round(r["bitscore"], 1), r["evalue_zref"],
        r["trunc_5p"], r["trunc_3p"], round(r["gc"], 2))


# --- the store -----------------------------------------------------------------------------------------


def _genomes() -> list[GenomeInput]:
    out = []
    for o in CONTRACT["objects"]:
        cu = ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)
        out.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"],
                               circular=(cu,) if o["topology"] == "circular" else ()))
    return out


def _rows(root: Path, rel: str) -> list[dict]:
    return [r for f in io.part_files(root / rel) for r in pq.read_table(f, partitioning=None).to_pylist()]


def _digests(root: Path) -> dict[str, str]:
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file()}


def _tapes(root: Path, producers: list[str]) -> dict[bytes, dict]:
    tv = tapes.tape_version(producers)
    return {r["contig_uid"]: r for r in _rows(root, f"l3/tapes/tape_version={tv}/stage=S0")}


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """The fixture stage built with the gene caller alone (inputs, ingest, genes, IGRs, tapes), then the ncRNA
    producer, then IGRs and tapes of both producers, in one root; with copies of the root before the ncRNA
    step and the files that step wrote."""
    base = tmp_path_factory.mktemp("ncrna")
    root, scratch = base / "root", base / "scratch"
    scratch.mkdir()
    ingest([SourceObject(FIX / o["file"], CONTRACT["source_name"], CONTRACT["source_release"], CONTRACT["license"],
                         True) for o in CONTRACT["objects"]], root, "S0")
    genomes = _genomes()
    write_genome_input(root, "S0", genomes)
    call_genes(root, "S0", genomes, tmp=scratch)
    igr.run(root, "S0", tmp=scratch)
    tapes.run(root, "S0", tmp=scratch)
    before = base / "before"  # the root as the ncRNA step finds it
    shutil.copytree(root, before)
    report = _ncrna().run(root, "S0", LIB, cpus=1, tmp=scratch, model_set=SET)
    after_ncrna = _digests(root)
    p = report["producer"]
    igr.run(root, "S0", producers=[V, p], tmp=scratch)
    tapes.run(root, "S0", producers=[V, p], tmp=scratch)
    plain = base / "lib.cm"
    plain.write_bytes(gzip.decompress(LIB.read_bytes()))
    new = {k: h for k, h in after_ncrna.items() if k not in _digests(before)}
    return {"root": root, "before": before, "report": report, "producer": p, "plain": plain, "base": base,
            "scratch": scratch, "ncrna_files": new, "genomes": genomes,
            "hits": _rows(root, f"l1/ncrna/producer={p}/stage=S0"),
            "features": _rows(root, f"l1/feature/producer={p}/stage=S0")}


@pytest.fixture(scope="module")
def oracle(world):
    """The binary's hits on the three fixture contigs read from the source FASTA (named by contig_uid)."""
    d = world["base"] / "oracle"
    d.mkdir()
    fna = d / "contigs.fna"
    with open(fna, "wb") as fh:
        for o in CONTRACT["objects"]:
            seq = fasta.read(FIX / o["file"])[0].seq
            fh.write(b">" + ids.contig_uid(seq).hex().encode() + b"\n" + seq + b"\n")
    return run_cmsearch(world["plain"], fna, d / "hits.tbl")


# --- the oracle ----------------------------------------------------------------------------------------


@pytest.mark.covers("C-01", "C-06", "LY-01")
def test_tables_conform(world):
    root, p = world["root"], world["producer"]
    for name, rel in (("l1_ncrna", f"l1/ncrna/producer={p}/stage=S0"), ("feature", f"l1/feature/producer={p}/stage=S0")):
        t = pq.read_table(io.part_files(root / rel)[0], partitioning=None)
        assert t.num_rows > 0
        assert conformance_problems(t, table(name)) == []


@pytest.mark.covers("O-SCAN", "MR-12")
def test_rows_equal_the_cmsearch_binary_after_the_coordinate_conversion(world, oracle):
    """Every hit the binary reports at or above GA (`--cut_ga`), on the plain contigs, is one builder row with the
    same coordinates (0-based half-open), strand, model, model coordinates, score, E-value at Z = 1e5 Mb,
    truncation and GC; and the builder has no other row."""
    ga = oracle_ga(world["plain"].read_text())
    assert set(ga) == FAMILIES
    assert oracle and {h["inc"] for h in oracle} == {"!"}
    assert all(h["score"] >= ga[h["acc"]] - 0.05 for h in oracle)  # printed with one decimal, GA with two
    want = dict(converted(h, bytes.fromhex(h["target"])) for h in oracle)
    assert len(want) == len(oracle)
    got = dict(row_key(r) for r in world["hits"])
    assert got == want
    # every family of the library is found, and the fixture holds hits on both strands, on all three contigs
    assert {k[4] for k in got} == FAMILIES
    assert {k[3] for k in got} == {1, -1} and len({k[0] for k in got}) == 3
    for r in world["hits"]:
        assert r["z_mb"] == 1e5 and r["ga_bits"] == pytest.approx(ga[r["model_acc"]])
        assert r["feature_uid"] == ids.feature_uid(r["contig_uid"], r["start"], r["end"], r["strand"], NCRNA)
        assert 0 <= r["start"] < r["end"] <= r["contig_len"]  # no fixture hit crosses an origin (C-05)


@pytest.mark.covers("C-02", "C-05")
def test_features_are_the_distinct_loci_of_the_hits(world):
    feats = world["features"]
    want = sorted({(r["feature_uid"], r["contig_uid"], r["start"], r["end"], r["strand"]) for r in world["hits"]},
                  key=lambda t: (t[1], t[2], t[3], t[4]))
    assert [(f["feature_uid"], f["contig_uid"], f["start"], f["end"], f["strand"]) for f in feats] == want
    assert {f["kind"] for f in feats} == {NCRNA}
    assert len({f["feature_uid"] for f in feats}) == len(feats)
    # tRNA and tRNA-Sec on overlapping loci of M. genitalium: both kept (no clan rule in the spec or decisions)
    sec = [r for r in world["hits"] if r["model_acc"] == "RF01852"]
    assert sec and all(any(t["model_acc"] == "RF00005" and t["contig_uid"] == s["contig_uid"]
                           and t["start"] < s["end"] and s["start"] < t["end"] for t in world["hits"]) for s in sec)


def test_producer_names_infernal_and_the_library(world, monkeypatch):
    help_text = subprocess.run([CMSEARCH, "-h"], capture_output=True, text=True, check=True).stdout
    infernal = re.search(r"^# INFERNAL (\S+)", help_text, re.M)[1]
    assert infernal == provenance.locked_versions()["infernal"]
    sha = hashlib.sha256(LIB.read_bytes()).hexdigest()
    set_sha = hashlib.sha256(SET.read_bytes()).hexdigest()  # card E-07: the set and the region rule
    assert world["producer"] == f"infernal-{infernal}-rfam-sha{sha[:12]}-sha{set_sha[:12]}-igr200nt"
    nc = _ncrna()
    # a file whose sha256 refdb.yaml lists is named by its version there (Rfam.cm.gz of 15.1 gives 15.1)
    rfam = next(db for db in provenance.refdb_libraries() if db["name"] == "rfam")
    assert rfam["version"] == "15.1"
    monkeypatch.setattr(provenance, "refdb_libraries", lambda repo=None: [{**rfam, "sha256": sha}])
    assert nc.library_version(LIB) == "15.1"
    assert nc.producer(LIB, SET) == f"infernal-{infernal}-rfam-15.1-sha{set_sha[:12]}-igr200nt"
    assert nc.producer(LIB) == f"infernal-{infernal}-rfam-15.1-prok1-igr200nt"  # the frozen set by default


@pytest.mark.covers("C-08", "L-01")
def test_ledger_records_tool_and_library_versions_and_balances(world):
    root, p = world["root"], world["producer"]
    step = f"ncrna-{p}"
    rows = ledger.read(root, "S0", step)
    assert {r["object_key"] for r in rows} == {cu for g in world["genomes"] for cu in g.contig_uids}
    lock = provenance.locked_versions()
    for r in rows:
        tv = json.loads(r["tool_versions_json"])
        assert tv["infernal"] == lock["infernal"] and tv["pyarrow"] == lock["pyarrow"]
        assert tv["rfam"] == _ncrna().library_version(LIB)  # card E-07: no longer the producer's last field
    by = {r["object_key"]: r["n_in"] for r in rows}
    assert sum(by.values()) == len(world["hits"])
    assert ledger.check_balance(root, "S0", [step]) == []
    # C-08 through a tape version: the producer of the refdb library passes, the fixture's subset is reported
    assert provenance._producer_problems("x", "p", f"infernal-{lock['infernal']}-rfam-15.1", lock) == []
    assert provenance._producer_problems("x", "p", p, lock)
    assert provenance._producer_problems("x", "p", "infernal-0.0-rfam-15.1", lock)


# --- circular contigs ----------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def planted(world, oracle):
    """A tRNA of M. genitalium (the best plus-strand RF00005 hit) planted across the origin of two circular
    contigs: a region around it, rotated so that the origin falls in its middle; one region longer than the
    library's longest window (flank) and one shorter, so the flanks repeat the contig."""
    mg = next(g.contig_uids[0] for g in world["genomes"] if g.circular)
    seq = next(fasta.read(FIX / o["file"])[0].seq for o in CONTRACT["objects"] if o["topology"] == "circular")
    trna = max((h for h in oracle if h["acc"] == "RF00005" and h["strand"] == "+" and h["target"] == mg.hex()),
               key=lambda h: h["score"])
    a, b = trna["from"] - 1, trna["to"]
    flank = world["report"]["flank"]
    lengths = (3000, (b - a) + 20)
    assert lengths[0] > flank > lengths[1]  # both ways of building a circular record are exercised
    d = world["base"] / "planted"
    d.mkdir()
    objs, genomes, want = [], [], {}
    for n in lengths:
        r0 = a - (n - (b - a)) // 2
        region = seq[r0:r0 + n]
        k = (a - r0) + (b - a) // 2
        rotated = region[k:] + region[:k]
        path = d / f"planted_{n}.fna"
        path.write_bytes(b">planted_" + str(n).encode() + b"\n" + rotated + b"\n")
        objs.append(SourceObject(path, "planted", "1", "INSDC (public)", True))
        cu = ids.contig_uid(rotated)
        genomes.append(GenomeInput((cu,), 4, kind="plasmid", circular=(cu,)))
        start = (a - r0 - k) % n
        want[cu] = (n, start, start + (b - a))
    root = d / "root"
    ingest(objs, root, "P")
    write_genome_input(root, "P", genomes)
    report = _ncrna().run(root, "P", LIB, cpus=1, tmp=world["scratch"], model_set=SET)
    hits = _rows(root, f"l1/ncrna/producer={report['producer']}/stage=P")
    return {"want": want, "hits": hits, "trna": trna, "lengths": lengths}


@pytest.mark.covers("C-05", "MR-04")
def test_a_trna_planted_across_the_origin_is_found_once_with_end_past_the_contig(planted):
    for cu, (n, start, end) in planted["want"].items():
        got = [r for r in planted["hits"] if r["contig_uid"] == cu and r["model_acc"] == "RF00005"]
        # the region holds other tRNAs of the M. genitalium cluster; the planted one is the only one across
        # the origin, and no other row covers any of its bases (it is not found a second time in a flank)
        cross = [r for r in got if r["end"] > n]
        assert len(cross) == 1, (n, got)
        r = cross[0]
        covered = {x % n for x in range(r["start"], r["end"])}
        assert not [o for o in got if o is not r and covered & {x % n for x in range(o["start"], o["end"])}]
        assert len({(o["start"], o["end"], o["strand"]) for o in got}) == len(got)
        assert r["contig_len"] == n and r["strand"] == 1 and not r["trunc_5p"] and not r["trunc_3p"]
        assert 0 <= r["start"] < n < r["end"] <= r["start"] + n  # C-05 on a circular contig
        if n == planted["lengths"][0]:  # with the whole context: the same hit as on the genome
            assert (r["start"], r["end"]) == (start, end)
            assert r["bitscore"] == pytest.approx(planted["trna"]["score"], abs=0.05)
        assert r["feature_uid"] == ids.feature_uid(cu, r["start"], r["end"], 1, NCRNA)


# --- determinism -----------------------------------------------------------------------------------------


def _outputs(base: Path, files: dict[str, str]) -> dict[str, str]:
    return {k: hashlib.sha256((base / k).read_bytes()).hexdigest() for k in files if (base / k).exists()}


@pytest.mark.covers("D-01", "D-02", "D-04", "MR-01", "MR-16")
@pytest.mark.parametrize("variant", ["4 threads", "reverse order, one contig per batch"])
def test_bytes_do_not_depend_on_threads_order_or_batches(world, tmp_path, variant):
    genomes = list(world["genomes"])
    kw = {"cpus": 4} if variant == "4 threads" else {"cpus": 1, "genomes": genomes[::-1], "batch_bases": 1}
    out = tmp_path / "out"
    rep = _ncrna().run(world["before"], "S0", LIB, tmp=world["scratch"], out=out, model_set=SET, **kw)
    if "batch_bases" in kw:
        assert rep["batches"] == len(genomes)
    assert sorted(_digests(out)) == sorted(world["ncrna_files"])
    assert _digests(out) == world["ncrna_files"]


# --- run-step ----------------------------------------------------------------------------------------------


@pytest.mark.covers("LY-10", "D-01")
def test_under_run_step_equals_the_direct_call(world, tmp_path):
    root = tmp_path / "root"
    shutil.copytree(world["before"], root)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    step = c08.Step("S0", "ncrna", lambda r: ["stage", "ncrna", "--root", str(r), "--stage", "S0", "--cm", str(LIB),
                                              "--ncrna-set", str(SET), "--cpus", "2", "--tmp", str(scratch)], None)
    prof = c08.profile_at(tmp_path / "prof")
    assert c08.run_cli(prof, root, step) == 0
    got, before = c08.tree(root), c08.tree(world["before"])
    assert {k: v for k, v in got.items() if k not in before} == world["ncrna_files"]
    assert {k: got[k] for k in before} == before
    assert c08.published(root, step) == sorted(world["ncrna_files"])


@pytest.mark.covers("D-03")
def test_refuses_before_any_work(world, tmp_path):
    nc = _ncrna()
    out = tmp_path / "out"
    with pytest.raises(FileExistsError, match="writes its partitions once"):
        nc.run(world["root"], "S0", LIB, cpus=1, tmp=world["scratch"], out=out, model_set=SET)
    assert not out.exists() or c08.tree(out) == {}
    planted = sorted(world["ncrna_files"])[-1]  # a killed attempt's file in the output directory
    (out / planted).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(world["root"] / planted, out / planted)
    with pytest.raises(FileExistsError, match="writes its partitions once"):
        nc.run(world["before"], "S0", LIB, cpus=1, tmp=world["scratch"], out=out, model_set=SET)
    assert sorted(c08.tree(out)) == [planted]
    for bad in (0, None, 1.5, True):
        with pytest.raises(ValueError, match="cpus"):
            nc.run(world["before"], "S0", LIB, cpus=bad, tmp=world["scratch"], out=tmp_path / "o2", model_set=SET)
    with pytest.raises(ValueError, match="scratch"):
        nc.run(world["before"], "S0", LIB, cpus=1, tmp=world["before"] / "tmp", out=tmp_path / "o3", model_set=SET)
    assert not (tmp_path / "o2").exists() and not (tmp_path / "o3").exists()


# --- IGRs, tapes and stability with both producers ------------------------------------------------------------


@pytest.mark.covers("C-04", "R-05")
def test_igr_and_tapes_build_with_both_producers_and_carry_the_ncrna_slots(world):
    root, p = world["root"], world["producer"]
    assert io.part_files(root / f"l0seq/igr/igr_version={igr.igr_version([V, p])}/stage=S0")
    tp = _tapes(root, [V, p])
    feats = {f["feature_uid"]: f for f in world["features"]}
    seen = set()
    for t in tp.values():
        for fu, kind, pu, s, e, strand in zip(t["feature_uids"], t["feature_kinds"], t["protein_uids"], t["start_bps"],
                                               t["end_bps"], t["strands"]):
            if kind == NCRNA:
                assert pu is None
                f = feats[fu]
                assert (f["contig_uid"], f["start"], f["end"], f["strand"]) == (t["contig_uid"], s, e, strand)
                seen.add(fu)
            else:
                assert pu is not None
    assert seen == set(feats)
    # an ncRNA inside an IGR of the gene caller's version splits it: both flanks name it
    flanks = {u for r in _rows(root, f"l0seq/igr/igr_version={igr.igr_version([V, p])}/stage=S0")
              for u in (r["left_feature_uid"], r["right_feature_uid"]) if u is not None}
    assert flanks & set(feats)


@pytest.mark.covers("MR-20")
def test_adding_the_ncrna_producer_changes_no_gene_call_row(world):
    before, root, p = world["before"], world["root"], world["producer"]
    # every file of the gene-call build is byte-identical afterwards
    old = _digests(before)
    new = _digests(root)
    assert {k: new[k] for k in old} == old
    # rows keyed by feature_uid in the gene-call tables: none changed, none added, none removed
    changed = 0
    for rel in (f"l1/feature/producer={V}/stage=S0", f"l1/feature_protein/gene_call_version={V}/stage=S0",
                f"l1/protein_occurrence/gene_call_version={V}/stage=S0"):
        a = {r["feature_uid"]: r for r in _rows(before, rel)}
        b = {r["feature_uid"]: r for r in _rows(root, rel)}
        assert a
        changed += len(set(a) ^ set(b)) + sum(a[u] != b[u] for u in set(a) & set(b))
    assert changed == 0
    # on the tape of both producers every CDS slot, keyed by feature_uid, is the slot of the gene-call tape
    one, both = _tapes(root, [V]), _tapes(root, [V, p])
    assert set(one) == set(both)
    for cu, t in one.items():
        u = both[cu]
        old_slots = {fu: (s, e, st, k, pu) for fu, s, e, st, k, pu in zip(
            t["feature_uids"], t["start_bps"], t["end_bps"], t["strands"], t["feature_kinds"], t["protein_uids"])}
        new_slots = {fu: (s, e, st, k, pu) for fu, s, e, st, k, pu in zip(
            u["feature_uids"], u["start_bps"], u["end_bps"], u["strands"], u["feature_kinds"], u["protein_uids"])}
        assert {fu: new_slots[fu] for fu in old_slots} == old_slots
        assert all(v[3] == NCRNA for fu, v in new_slots.items() if fu not in old_slots)
    assert any(both[cu]["n_slots"] > one[cu]["n_slots"] for cu in one)


def test_the_genome_inputs_are_read_when_none_are_given(world):
    """`genomes=None` reads `genome_input` (sorted by genome_uid): the run above passed none."""
    assert {g.contig_uids for g in read_genome_input(world["root"], "S0")} == {g.contig_uids for g in world["genomes"]}
    assert world["report"]["contigs"] == 3 and world["report"]["batches"] == 1
