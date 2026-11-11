"""Ring R2, card E-07: the ncRNA producer's intergenic regions, what it sends to cmsearch, and how hits map back.

Oracles (AGENTS §5.2), each written here and independent of `biobox.l1.ncrna`:

- regions: a per-base loop. A base is in a region when it lies within the margin of a base that no gene covers
  (on a circular contig, distance around the circle; on a linear one, inside the contig); regions are the maximal
  runs of such bases, a run across the origin of a circular contig written with end > n, the whole circle as
  (0, n). It is compared with `ncrna.regions` on seeded random layouts and edge cases, and with the records the
  producer actually hands to cmsearch for the micro-S0 fixture's gene calls (read from the stage's feature table).
- what is searched: a stand-in `cmsearch` first on PATH records each call's arguments, CM file and FASTA records,
  then runs the binary. Every intergenic model is searched only on the oracle's regions, cut from the source FASTA,
  with `--notrunc`; every whole-genome model only on whole contigs; no model in both; no region or contig twice.
- mapping: the `cmsearch` binary's own coordinates on the whole contig. On the fixture, each intergenic hit of the
  producer is a hit the binary reports on the whole contig (same model, strand, start and end). Across the origin:
  a tRNA of M. genitalium planted across the origin of a circular contig whose only intergenic region crosses it,
  against the binary's hit on the same sequence unrotated, moved by the known rotation.
"""

from __future__ import annotations

import gzip
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

from biobox import fasta, ids, io
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l1.genecall import run as call_genes
from biobox.registry.schemas import table

CMSEARCH = shutil.which("cmsearch")
pytestmark = [pytest.mark.r2, pytest.mark.covers("E-07")]
needs_infernal = pytest.mark.skipif(CMSEARCH is None or shutil.which("cmfetch") is None,
                                    reason="Infernal (cmsearch, cmfetch) not on PATH; it is in the test image")

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
LIB = FIX / "rfam15.1_micro_s0.cm.gz"
MARGIN = 200  # DEC-28
Z = "100000"  # DEC-04
# a proper subset of the fixture library's nine models: two whole-genome, five intergenic (RF00001 and RF00028 left out)
MIXED = {"RF01804": ("Lambda_thermo", 1), "RF02111": ("IS009", 1), "RF00005": ("tRNA", 0),
         "RF01852": ("tRNA-Sec", 0), "RF00023": ("tmRNA", 0), "RF00169": ("Bacteria_small_SRP", 0),
         "RF00011": ("RNaseP_bact_b", 0)}


def _ncrna():
    from biobox.l1 import ncrna  # imported here so that a missing function fails the tests, not the collection

    return ncrna


# --- the per-base oracle ---------------------------------------------------------------------------------------


def oracle_regions(genes, n: int, circular: bool, margin: int) -> list[tuple[int, int]]:
    """Maximal runs of bases within `margin` of a base no gene covers, by a loop over bases."""
    covered = [False] * n
    for s, e in genes:
        for x in range(s, e):
            covered[x % n] = True
    near = [False] * n
    for y in range(n):
        if covered[y]:
            continue
        for x in range(y - margin, y + margin + 1):
            if circular:
                near[x % n] = True
            elif 0 <= x < n:
                near[x] = True
    if not any(near):
        return []
    if all(near):
        return [(0, n)]
    out = []
    for x in range(n):
        if near[x] and (not near[(x - 1) % n] if circular else (x == 0 or not near[x - 1])):
            e = x
            while (near[e % n] if circular else e < n and near[e]):
                e += 1
            out.append((x, e))
    return sorted(out)


def random_layout(rng: random.Random, n: int, circular: bool) -> list[tuple[int, int]]:
    genes = []
    for _ in range(rng.randint(0, 12)):
        s = rng.randrange(n)
        length = rng.randint(1, max(1, n // 3))
        e = s + length if circular else min(n, s + length)
        genes.append((s, e))
    return genes


EDGES = [  # (genes, n, circular)
    ([], 500, False), ([], 500, True), ([(0, 500)], 500, False), ([(0, 500)], 500, True),
    ([(100, 200), (200, 300)], 1000, False), ([(0, 100)], 1000, False), ([(900, 1000)], 1000, False),
    ([(900, 1100)], 1000, True), ([(950, 1000)], 1000, True), ([(0, 50), (950, 1000)], 1000, True),
    ([(10, 990)], 1000, True), ([(10, 990)], 1000, False), ([(300, 1290)], 1000, True),
    ([(500, 600), (100, 450), (650, 900)], 1000, True), ([(1, 999)], 1000, True), ([(0, 999)], 1000, True),
]


@pytest.mark.parametrize("margin", [0, 1, 7, 200])
def test_regions_equal_the_per_base_loop(margin):
    ncrna = _ncrna()
    rng = random.Random(20261111 + margin)
    cases = list(EDGES)
    for _ in range(60):
        circular = rng.random() < 0.5
        n = rng.choice([1, 2, 37, 400, 1000, 2500])
        cases.append((random_layout(rng, n, circular), n, circular))
    for genes, n, circular in cases:
        want = oracle_regions(genes, n, circular, margin)
        assert ncrna.regions(genes, n, circular, margin) == want, (genes, n, circular, margin)


# --- a stand-in cmsearch that records what it is given --------------------------------------------------------


def shim(bin_dir: Path, log: Path) -> None:
    """`bin_dir/cmsearch`: records argv, the CM file and the FASTA of each search into its own directory under
    `log`, then runs the binary."""
    bin_dir.mkdir(parents=True)
    log.mkdir(parents=True)
    (bin_dir / "cmsearch").write_text(f"""#!{sys.executable}
import os, shutil, sys, uuid
args = sys.argv[1:]
if "--tblout" in args:
    d = os.path.join({str(log)!r}, uuid.uuid4().hex)
    os.makedirs(d)
    shutil.copy(args[-2], os.path.join(d, "cm"))
    shutil.copy(args[-1], os.path.join(d, "fna"))
    with open(os.path.join(d, "argv"), "w") as fh:
        fh.write("\\n".join(args))
os.execv({CMSEARCH!r}, [{CMSEARCH!r}] + args)
""")
    (bin_dir / "cmsearch").chmod(0o755)


def calls(log: Path) -> list[dict]:
    """Each recorded search: whether it ran with --notrunc, its models' accessions, its records (name -> bases)."""
    out = []
    for d in sorted(log.iterdir()):
        accs, inside = [], False  # each model's ACC in its CM header (its filter HMM repeats it)
        for ln in (d / "cm").read_text().splitlines():
            if ln.startswith("INFERNAL1/"):
                inside = True
            elif ln.startswith("HMMER3/"):
                inside = False
            elif inside and ln.startswith("ACC "):
                accs.append(ln.split()[1])
        recs = {r.accession: r.seq for r in fasta.read(d / "fna")}
        out.append({"notrunc": "--notrunc" in (d / "argv").read_text().split("\n"), "models": accs, "records": recs})
    return out


def write_set(path: Path, models: dict[str, tuple[str, int]]) -> Path:
    path.write_text("acc\tname\twhole_genome\n" + "".join(f"{a}\t{n}\t{w}\n" for a, (n, w) in models.items()))
    return path


# --- the cmsearch binary, run here ----------------------------------------------------------------------------


def binary_hits(cm: Path, fna: Path, tbl: Path, *flags: str) -> list[dict]:
    subprocess.run([CMSEARCH, "--cpu", "2", "-Z", Z, "--cut_ga", *flags, "--noali", "--tblout", str(tbl),
                    "-o", os.devnull, str(cm), str(fna)], check=True)
    out = []
    for line in tbl.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split()
        lo, hi = sorted((int(f[7]), int(f[8])))
        out.append({"target": f[0], "acc": f[3], "start": lo - 1, "end": hi, "strand": 1 if f[9] == "+" else -1,
                    "score": float(f[14]), "evalue": float(f[15])})
    return out


# --- the fixture -------------------------------------------------------------------------------------------------


def fixture_genomes() -> list[tuple[GenomeInput, bytes]]:
    out = []
    for o in CONTRACT["objects"]:
        seq = fasta.read(FIX / o["file"])[0].seq
        cu = ids.contig_uid(seq)
        out.append((GenomeInput((cu,), o["genetic_code"], kind=o["kind"],
                                circular=(cu,) if o["topology"] == "circular" else ()), seq.upper()))
    return out


def gene_calls(root: Path, stage: str) -> dict[bytes, list[tuple[int, int]]]:
    out: dict[bytes, list] = {}
    for f in io.part_files(root / f"l1/feature/producer={V}/stage={stage}"):
        for r in pq.read_table(f, partitioning=None).to_pylist():
            out.setdefault(r["contig_uid"], []).append((r["start"], r["end"]))
    return out


def write_genes(root: Path, stage: str, cu: bytes, genes: list[tuple[int, int, int]]) -> None:
    """Constructed gene calls (start, end, strand) of one contig as the stage's gene-call feature rows; none
    writes the empty partition of a stage whose gene caller found no gene."""
    spec = table("feature")
    rs = [{"feature_uid": ids.feature_uid(cu, s, e, st, ids.KIND["cds"]), "contig_uid": cu, "start": s, "end": e,
           "strand": st, "kind": ids.KIND["cds"]} for s, e, st in sorted(genes)]
    io.write_parts([pa.RecordBatch.from_pylist(rs, schema=spec.schema)], spec,
                   root / f"l1/feature/producer={V}/stage={stage}", 24)


def rows(root: Path, rel: str) -> list[dict]:
    return [r for f in io.part_files(root / rel) for r in pq.read_table(f, partitioning=None).to_pylist()]


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """The micro-S0 stage with its gene calls, scanned by the producer with the MIXED set under the recording shim."""
    if CMSEARCH is None or shutil.which("cmfetch") is None:
        pytest.skip("Infernal not on PATH; it is in the test image")
    base = tmp_path_factory.mktemp("e07-regions")
    root, scratch = base / "root", base / "scratch"
    scratch.mkdir()
    ingest([SourceObject(FIX / o["file"], CONTRACT["source_name"], CONTRACT["source_release"], CONTRACT["license"],
                         True) for o in CONTRACT["objects"]], root, "S0")
    gs = fixture_genomes()
    write_genome_input(root, "S0", [g for g, _ in gs])
    call_genes(root, "S0", [g for g, _ in gs], tmp=scratch)
    model_set = write_set(base / "mixed.tsv", MIXED)
    shim(base / "bin", base / "log")
    old = os.environ["PATH"]
    os.environ["PATH"] = f"{base / 'bin'}{os.pathsep}{old}"
    try:
        report = _ncrna().run(root, "S0", LIB, cpus=2, tmp=scratch, model_set=model_set)
    finally:
        os.environ["PATH"] = old
    plain = base / "lib.cm"
    plain.write_bytes(gzip.decompress(LIB.read_bytes()))
    return {"root": root, "base": base, "report": report, "set": model_set, "genomes": gs, "plain": plain,
            "calls": calls(base / "log"), "genes": gene_calls(root, "S0"),
            "hits": rows(root, f"l1/ncrna/producer={report['producer']}/stage=S0")}


@needs_infernal
def test_intergenic_models_search_exactly_the_oracle_regions_with_notrunc(world):
    want_models = sorted(a for a, (_, w) in MIXED.items() if not w)
    igr = [c for c in world["calls"] if c["notrunc"]]
    assert igr and all(sorted(c["models"]) == want_models for c in igr)
    seen: dict[tuple[bytes, int, int], int] = {}
    for c in igr:
        for name, rec in c["records"].items():
            hexu, s = name.rsplit("_", 1)
            key = (bytes.fromhex(hexu), int(s), int(s) + len(rec))
            seen[key] = seen.get(key, 0) + 1
            g, seq = next((g, q) for g, q in world["genomes"] if g.contig_uids[0] == key[0])
            ring = seq + seq if g.circular else seq
            assert rec == ring[key[1]:key[2]], name  # the record is the region's bases, cut from the source FASTA
    want = {}
    for g, seq in world["genomes"]:
        cu = g.contig_uids[0]
        for s, e in oracle_regions(world["genes"].get(cu, []), len(seq), bool(g.circular), MARGIN):
            want[(cu, s, e)] = 1
    assert seen == want  # each oracle region once, nothing else
    assert world["report"]["regions"] == len(want)


@needs_infernal
def test_whole_genome_models_search_each_whole_contig_once_and_no_model_searches_twice(world):
    want_models = sorted(a for a, (_, w) in MIXED.items() if w)
    whole = [c for c in world["calls"] if not c["notrunc"]]
    assert whole and all(sorted(c["models"]) == want_models for c in whole)
    names = [n for c in whole for n in c["records"]]
    assert sorted(names) == sorted(g.contig_uids[0].hex() for g, _ in world["genomes"])  # each contig once
    for c in whole:
        for name, rec in c["records"].items():
            g, seq = next((g, q) for g, q in world["genomes"] if g.contig_uids[0].hex() == name)
            if not g.circular:
                assert rec == seq
            else:  # the contig with the same number of its bases on each side
                f = (len(rec) - len(seq)) // 2
                assert f > 0 and rec[f:f + len(seq)] == seq and rec[:f] == (seq * 2)[-f:] and rec[-f:] == (seq * 2)[:f]
    igr_models = {a for c in world["calls"] if c["notrunc"] for a in c["models"]}
    whole_models = {a for c in whole for a in c["models"]}
    assert not igr_models & whole_models and igr_models | whole_models == set(MIXED)


@needs_infernal
def test_intergenic_hits_are_the_binarys_hits_on_the_whole_contig(world, tmp_path):
    fna = tmp_path / "contigs.fna"
    fna.write_bytes(b"".join(b">" + g.contig_uids[0].hex().encode() + b"\n" + s + b"\n" for g, s in world["genomes"]))
    truth = {(bytes.fromhex(h["target"]), h["acc"], h["strand"], h["start"], h["end"])
             for h in binary_hits(world["plain"], fna, tmp_path / "wg.tbl")}
    got = [h for h in world["hits"] if not MIXED[h["model_acc"]][1]]
    assert got and {h["model_acc"] for h in got} >= {"RF00005", "RF00023"}
    missing = [h for h in got if (h["contig_uid"], h["model_acc"], h["strand"], h["start"], h["end"]) not in truth]
    assert not missing, missing
    assert all(not h["trunc_5p"] and not h["trunc_3p"] for h in got)  # --notrunc


# --- across the origin ---------------------------------------------------------------------------------------


@needs_infernal
def test_a_hit_in_a_region_across_the_origin_maps_to_the_binarys_coordinates(world, tmp_path):
    """A 3,000-base piece of M. genitalium around its best plus-strand tRNA, rotated so that the origin falls in the
    tRNA's middle; two touching genes cover [500, 2400), so the one intergenic region, [2200, 3700), crosses the
    origin and holds the tRNA."""
    g, mg = next((g, s) for g, s in world["genomes"] if g.circular)
    fna = tmp_path / "mg.fna"
    fna.write_bytes(b">mg\n" + mg + b"\n")
    trna = max((h for h in binary_hits(world["plain"], fna, tmp_path / "mg.tbl")
                if h["acc"] == "RF00005" and h["strand"] == 1), key=lambda h: h["score"])
    a, b, n = trna["start"], trna["end"], 3000
    r0 = a - (n - (b - a)) // 2
    region = mg[r0:r0 + n]
    k = (a - r0) + (b - a) // 2
    rotated = region[k:] + region[:k]
    start = (a - r0 - k) % n  # where the tRNA starts on the rotated contig
    assert start > 2400 and start + (b - a) > n  # it crosses the origin, outside the genes
    # the binary on the unrotated piece (a linear record holding the whole tRNA)
    piece = tmp_path / "piece.fna"
    piece.write_bytes(b">piece\n" + region + b"\n")
    one = [h for h in binary_hits(world["plain"], piece, tmp_path / "piece.tbl")
           if h["acc"] == "RF00005" and h["strand"] == 1 and h["start"] == a - r0]
    assert len(one) == 1
    # the stage: the rotated contig, circular, with two constructed gene calls
    root, d = tmp_path / "root", tmp_path / "in"
    d.mkdir()
    (d / "rot.fna").write_bytes(b">rot\n" + rotated + b"\n")
    ingest([SourceObject(d / "rot.fna", "planted", "1", "INSDC (public)", True)], root, "P")
    cu = ids.contig_uid(rotated)
    write_genome_input(root, "P", [GenomeInput((cu,), 4, kind="plasmid", circular=(cu,))])
    write_genes(root, "P", cu, [(500, 1400, 1), (1400, 2400, -1)])
    genes = [(500, 1400, 1), (1400, 2400, -1)]
    assert oracle_regions([(s, e) for s, e, _ in genes], n, True, MARGIN) == [(2200, 3700)]
    rep = _ncrna().run(root, "P", LIB, cpus=1, tmp=tmp_path, model_set=world["set"])
    got = [r for r in rows(root, f"l1/ncrna/producer={rep['producer']}/stage=P")
           if r["model_acc"] == "RF00005" and r["end"] > n]
    assert len(got) == 1
    r = got[0]
    assert (r["start"], r["end"], r["strand"]) == (start, start + (b - a), 1)
    assert r["bitscore"] == pytest.approx(one[0]["score"], abs=0.05) and r["contig_len"] == n
    assert r["feature_uid"] == ids.feature_uid(cu, start, start + (b - a), 1, ids.KIND["ncrna"])
