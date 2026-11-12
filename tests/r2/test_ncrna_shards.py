"""Ring R2, card E-07: the ncRNA producer's model subset, and its scan in shards.

The stage: lambda (linear), M. genitalium (circular) and T4 cut into two linear contigs of one genome, gene-called;
the model set is a proper subset of the fixture library (`test_ncrna_regions.MIXED`: two whole-genome models, five
intergenic ones).

- The subset's hits equal the full fixture library's hits restricted to the subset. Oracle: the `cmsearch` binary
  with all nine models, on the whole contigs (for the whole-genome models) and on the per-base oracle's regions cut
  from the source FASTA with `--notrunc` (for the intergenic ones), its hits mapped back here.
- Sharded output is byte-identical to the unsharded run's at 1, 2 and 4 shards; the shards are consecutive
  `contig_uid` ranges cut at the limit (a pure function of the sorted contigs and the limit).
- A shard whose scratch runs are lost after its record was published is scanned again to the same sha256, and the
  merge then gives the same files.
"""

from __future__ import annotations

import gzip
import hashlib
import sys
from pathlib import Path

import pytest

from biobox import fasta, ids
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l1.genecall import run as call_genes

sys.path.insert(0, str(Path(__file__).parent))
import test_ncrna_regions as rg  # noqa: E402

pytestmark = [pytest.mark.r2, pytest.mark.covers("E-07", "D-01", "MR-12"), rg.needs_infernal]


def _ncrna():
    from biobox.l1 import ncrna

    return ncrna


@pytest.fixture(scope="module")
def stage(tmp_path_factory):
    base = tmp_path_factory.mktemp("e07-shards")
    root, scratch, src = base / "root", base / "scratch", base / "src"
    scratch.mkdir()
    src.mkdir()
    objs, genomes, seqs = [], [], {}
    for o in rg.CONTRACT["objects"]:
        seq = fasta.read(rg.FIX / o["file"])[0].seq.upper()
        if o["topology"] == "circular" or len(seq) < 100_000:
            parts = [seq]
            objs.append(SourceObject(rg.FIX / o["file"], rg.CONTRACT["source_name"], rg.CONTRACT["source_release"],
                                     rg.CONTRACT["license"], True))
        else:  # T4: two linear contigs of one genome
            parts = [seq[:len(seq) // 2], seq[len(seq) // 2:]]
            path = src / f"{o['accession']}.split.fna"
            path.write_bytes(b"".join(b">%s_%d\n%s\n" % (o["accession"].encode(), i, p) for i, p in enumerate(parts)))
            objs.append(SourceObject(path, "split", "1", rg.CONTRACT["license"], True))
        cus = tuple(ids.contig_uid(p) for p in parts)
        seqs.update(zip(cus, parts))
        genomes.append(GenomeInput(cus, o["genetic_code"], kind=o["kind"],
                                   circular=cus if o["topology"] == "circular" else ()))
    ingest(objs, root, "S")
    write_genome_input(root, "S", genomes)
    call_genes(root, "S", genomes, tmp=scratch)
    model_set = rg.write_set(base / "mixed.tsv", rg.MIXED)
    out = base / "unsharded"
    report = _ncrna().run(root, "S", rg.LIB, cpus=2, tmp=scratch, out=out, model_set=model_set)
    plain = base / "lib.cm"
    plain.write_bytes(gzip.decompress(rg.LIB.read_bytes()))
    circ = {cu for g in genomes for cu in g.circular}
    return {"root": root, "base": base, "scratch": scratch, "set": model_set, "out": out, "report": report,
            "seqs": seqs, "circular": circ, "plain": plain, "producer": report["producer"]}


def files(base: Path, prod: str) -> dict[str, str]:
    """sha256 of the producer's files (features, hits, ledger) under `base`."""
    out = {}
    for rel in (f"l1/feature/producer={prod}", f"l1/ncrna/producer={prod}", f"ledger/stage=S/step=ncrna-{prod}"):
        for p in sorted((base / rel).rglob("*")):
            if p.is_file():
                out[p.relative_to(base).as_posix()] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out


# --- the subset against the full library -------------------------------------------------------------------------


def test_the_subsets_hits_are_the_full_librarys_hits_restricted_to_it(stage, tmp_path):
    seqs, circ = stage["seqs"], stage["circular"]
    whole = {a for a, (_, w) in rg.MIXED.items() if w}
    inter = set(rg.MIXED) - whole
    # whole contigs, all nine models
    fna = tmp_path / "whole.fna"
    fna.write_bytes(b"".join(b">" + cu.hex().encode() + b"\n" + s + b"\n" for cu, s in sorted(seqs.items())))
    want = {(bytes.fromhex(h["target"]), h["start"], h["end"], h["strand"], h["acc"], h["score"])
            for h in rg.binary_hits(stage["plain"], fna, tmp_path / "whole.tbl") if h["acc"] in whole}
    # the oracle's regions, all nine models, --notrunc
    genes = rg.gene_calls(stage["root"], "S")
    recs, where = [], {}
    for cu, s in sorted(seqs.items()):
        ring = s + s if cu in circ else s
        for a, b in rg.oracle_regions(genes.get(cu, []), len(s), cu in circ, rg.MARGIN):
            name = f"r{len(recs)}"
            recs.append(b">" + name.encode() + b"\n" + ring[a:b] + b"\n")
            where[name] = (cu, a, len(s))
    fna = tmp_path / "regions.fna"
    fna.write_bytes(b"".join(recs))
    for h in rg.binary_hits(stage["plain"], fna, tmp_path / "regions.tbl", "--notrunc"):
        if h["acc"] not in inter:
            continue
        cu, a, n = where[h["target"]]
        s, e = a + h["start"], a + h["end"]
        if s >= n:
            s, e = s - n, e - n
        want.add((cu, s, e, h["strand"], h["acc"], h["score"]))
    got = {(r["contig_uid"], r["start"], r["end"], r["strand"], r["model_acc"], round(r["bitscore"], 1))
           for r in rg.rows(stage["out"], f"l1/ncrna/producer={stage['producer']}/stage=S")}
    assert want and {k[4] for k in want} >= {"RF00005", "RF00023", "RF01804"}
    assert got == want


# --- shards -------------------------------------------------------------------------------------------------------


def greedy(lengths: list[int], limit: int) -> list[list[int]]:
    """Consecutive groups of at most `limit` bases; a contig above it is a group of its own."""
    out: list[list[int]] = []
    for n in lengths:
        if not out or sum(out[-1]) + n > limit:
            out.append([])
        out[-1].append(n)
    return out


def limit_for(lengths: list[int], k: int) -> int:
    for limit in sorted({sum(lengths[i:j]) for i in range(len(lengths)) for j in range(i + 1, len(lengths) + 1)}):
        if len(greedy(lengths, limit)) == k:
            return limit
    if len(greedy(lengths, 1)) == k:
        return 1
    raise AssertionError(f"no limit cuts {lengths} into {k} shards")


@pytest.mark.parametrize("k", [1, 2, 4])
def test_sharded_is_byte_identical_to_unsharded(stage, tmp_path, k):
    nc = _ncrna()
    root, ms = stage["root"], stage["set"]
    order = sorted(stage["seqs"])
    lengths = [len(stage["seqs"][cu]) for cu in order]
    assert len(order) == 4
    limit = limit_for(lengths, k)
    out, scratch = tmp_path / "out", tmp_path / "scratch"
    scratch.mkdir()
    res = nc.plan(root, "S", rg.LIB, limit, model_set=ms, tmp=scratch, out=out)
    assert res["shards"] == k and res["contigs"] == 4 and res["bases"] == sum(lengths)
    plan = nc.read_plan(root, "S", stage["producer"], out)
    groups = greedy(lengths, limit)
    i = 0
    for s, g in zip(plan["shards"], groups):
        assert (s["first"], s["last"], s["contigs"], s["bases"]) == (order[i].hex(), order[i + len(g) - 1].hex(),
                                                                       len(g), sum(g))
        i += len(g)
    for j in reversed(range(k)):  # any order
        nc.scan_shard(root, "S", rg.LIB, j, 1, model_set=ms, tmp=scratch, out=out)
    nc.merge(root, "S", rg.LIB, model_set=ms, tmp=scratch, out=out)
    want = files(stage["out"], stage["producer"])
    assert want and files(out, stage["producer"]) == want


def test_a_shard_whose_runs_were_lost_is_scanned_again_to_its_record(stage, tmp_path):
    nc = _ncrna()
    root, ms, prod = stage["root"], stage["set"], stage["producer"]
    lengths = [len(stage["seqs"][cu]) for cu in sorted(stage["seqs"])]
    # the plan and records are published in a copy of the root, as build-stage publishes them
    import shutil

    copy = tmp_path / "root"
    shutil.copytree(root, copy, symlinks=True)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    nc.plan(copy, "S", rg.LIB, limit_for(lengths, 2), model_set=ms, tmp=scratch)
    for j in range(2):
        nc.scan_shard(copy, "S", rg.LIB, j, 1, model_set=ms, tmp=scratch)
    hits, contigs = nc.run_paths(scratch, copy, "S", prod, 1)
    assert nc.run_problem(copy, "S", prod, 1, scratch) is None
    hits.unlink()
    assert nc.run_problem(copy, "S", prod, 1, scratch) == "is missing"
    with pytest.raises(FileNotFoundError, match="not in the scratch directory"):
        nc.merge(copy, "S", rg.LIB, model_set=ms, tmp=scratch, out=tmp_path / "o")
    res = nc.rescan_shard(copy, "S", rg.LIB, 1, 2, model_set=ms, tmp=scratch)
    assert res["hits_sha256"] == nc.shard_record(copy, "S", prod, 1)["hits_sha256"]
    assert nc.run_problem(copy, "S", prod, 1, scratch) is None
    contigs.write_bytes(contigs.read_bytes() + b"x")  # damaged
    assert nc.run_problem(copy, "S", prod, 1, scratch) == "differs from the sha256 its record holds"
    contigs.unlink()
    nc.rescan_shard(copy, "S", rg.LIB, 1, 1, model_set=ms, tmp=scratch)
    nc.merge(copy, "S", rg.LIB, model_set=ms, tmp=scratch)
    assert files(copy, prod) == files(stage["out"], prod)
    with pytest.raises(FileExistsError):
        nc.scan_shard(copy, "S", rg.LIB, 0, 1, model_set=ms, tmp=scratch)


@pytest.mark.covers("E-08", "LY-10")
def test_records_are_the_ncrna_scan_artifact_and_e07_records_under_domain_scan_stay_readable(stage, tmp_path):
    """E-07 report (b): the plan and shard records are the registry's `ncrna_scan` artifact
    (`l1/ncrna_scan/producer=<p>/stage=<S>/`). A stage planned by card E-07's code has them under the `domain_scan`
    artifact (`l1/domain_scan/library=<p>/stage=<S>/`); they are read from there, a plan there counts as written, and
    the merge gives the same files."""
    import shutil

    from biobox.registry.schemas import ARTIFACTS

    nc = _ncrna()
    spec = ARTIFACTS["ncrna_scan"]
    assert (spec.location, spec.partition_keys, spec.in_manifest, spec.small_ok) == (
        "l1/ncrna_scan", ("producer", "stage"), False, True)
    root, ms, prod = stage["root"], stage["set"], stage["producer"]
    lengths = [len(stage["seqs"][cu]) for cu in sorted(stage["seqs"])]
    copy = tmp_path / "root"
    shutil.copytree(root, copy, symlinks=True)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    nc.plan(copy, "S", rg.LIB, limit_for(lengths, 2), model_set=ms, tmp=scratch)
    for j in range(2):
        nc.scan_shard(copy, "S", rg.LIB, j, 1, model_set=ms, tmp=scratch)
    new = copy / f"l1/ncrna_scan/producer={prod}/stage=S"
    assert sorted(p.name for p in new.iterdir()) == ["plan.json", "shard-00000.json", "shard-00001.json"]
    assert not (copy / "l1/domain_scan").exists()
    old = copy / f"l1/domain_scan/library={prod}/stage=S"  # where card E-07 wrote them
    old.parent.mkdir(parents=True)
    new.rename(old)
    assert nc.read_plan(copy, "S", prod)["producer"] == prod
    assert nc.find_plan(copy, "S", prod) == old / "plan.json"
    assert nc.shard_record(copy, "S", prod, 1)["shard"] == 1
    assert nc.run_problem(copy, "S", prod, 1, scratch) is None
    with pytest.raises(FileExistsError):
        nc.plan(copy, "S", rg.LIB, 1, model_set=ms, tmp=scratch)
    with pytest.raises(FileExistsError):
        nc.scan_shard(copy, "S", rg.LIB, 0, 1, model_set=ms, tmp=scratch)
    res = nc.merge(copy, "S", rg.LIB, model_set=ms, tmp=scratch)
    assert files(copy, prod) == files(stage["out"], prod)
    assert res["gene_call_version"] == V and not new.exists()


@pytest.mark.covers("E-08", "C-08")
def test_the_ledger_and_the_result_name_the_gene_call_version_scanned(stage):
    """E-07 report (c): the gene-call version whose intergenic regions were scanned is in the result and in every
    ledger row's tool versions (not in the producer name)."""
    import json

    from biobox import ledger

    assert stage["report"]["gene_call_version"] == V
    assert V not in stage["producer"]
    rows = ledger.read(stage["out"], "S", f"ncrna-{stage['producer']}")
    assert rows and all(json.loads(r["tool_versions_json"])["gene_call_version"] == V for r in rows)
    assert all(json.loads(r["tool_versions_json"])["rfam_clanin"] == "15.1" for r in rows)


def test_a_shard_refuses_another_library_or_set(stage, tmp_path):
    nc = _ncrna()
    root, ms = stage["root"], stage["set"]
    lengths = [len(stage["seqs"][cu]) for cu in sorted(stage["seqs"])]
    out, scratch = tmp_path / "out", tmp_path / "scratch"
    scratch.mkdir()
    nc.plan(root, "S", rg.LIB, sum(lengths), model_set=ms, tmp=scratch, out=out)
    with pytest.raises(ValueError, match="shard 3"):
        nc.scan_shard(root, "S", rg.LIB, 3, 1, model_set=ms, tmp=scratch, out=out)
    other = rg.write_set(tmp_path / "other.tsv", {**rg.MIXED, "RF00001": ("5S_rRNA", 0)})
    with pytest.raises(FileNotFoundError, match="no shard plan"):  # another set is another producer: no plan
        nc.scan_shard(root, "S", rg.LIB, 0, 1, model_set=other, tmp=scratch, out=out)
    with pytest.raises(FileExistsError):
        nc.plan(root, "S", rg.LIB, 1, model_set=ms, tmp=scratch, out=out)
    with pytest.raises(ValueError, match="not in the library"):
        nc.run(root, "S", rg.LIB, 1, model_set=rg.write_set(tmp_path / "x.tsv", {"RF99999": ("none", 0)}),
               tmp=scratch, out=tmp_path / "o2")
