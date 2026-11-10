"""Ring R3: the golden store checked by tools the project did not write (card D-03).

Ports the hand-run check of 2026-11-09 (scripts/stress/step6_*.py; measurements file §2). Runs only when
BIOBOX_GOLDEN_DIR is set, like test_golden.py. Each external binary that is absent skips its tests with a reason;
the pure-Python comparison logic has unit tests that need no binary and no golden set (class "Unit" below).

Measured on 2026-11-09 (laptop-built store, tools on the pod: samtools 1.19.2, Prodigal 2.6.3):
slices 33,713 of 33,713 identical; .fai and .gzi byte-identical; Prodigal 32,195 of 32,198 genes exact.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import pyarrow.parquet as pq
import pytest
import yaml

from biobox import fasta
from biobox.l0seq import SeqStore
from biobox.l1.genecall import GENE_CALL_VERSION

pytestmark = pytest.mark.r3

ROOT = Path(__file__).resolve().parents[2]
GOLDEN = yaml.safe_load((ROOT / "golden.yaml").read_text())
LOCK = yaml.safe_load((ROOT / "golden.lock.yaml").read_text()) if (ROOT / "golden.lock.yaml").exists() else None
_DIR = os.environ.get("BIOBOX_GOLDEN_DIR")

REFSEQ_GFF = "truth/GCF_000005845.2_ASM584v2_genomic.gff.gz"
# sha256 of the RefSeq GFF for E. coli K-12 MG1655 (GCF_000005845.2_ASM584v2), measured 2026-11-09 on the copy
# used for the stress test. Interim guard: the proper home is an L0 source with a contract in sources.yaml (E-03).
REFSEQ_GFF_SHA256 = "afdf03dc1d06e423d874ee29d9e0df14f5d32baf893f4da9f93aa721eb5f495e"
EDGE_BP = 50_000  # distance from a circular contig end that counts as "at the end" for a Prodigal difference
MAX_DIFF_PER_CIRCULAR_CONTIG = 6  # measured 2026-11-09: at most 4 (B. theta chromosome: 2 bio-box-only + 2 Prodigal-only)


# ---------------------------------------------------------------- pure logic (unit-tested below)


def region(rd_hex: str, start: int, end: int) -> str:
    """samtools region string for a 0-based half-open slice: `name:start+1-end` (1-based inclusive)."""
    return f"{rd_hex}:{start + 1}-{end}"


def parse_fasta_bytes(data: bytes) -> list[tuple[str, bytes]]:
    """(name, sequence) per record of a FASTA stream; name is the first word of the header."""
    out: list[tuple[str, bytes]] = []
    for chunk in data.split(b"\n>"):
        chunk = chunk.lstrip(b">")
        if not chunk.strip():
            continue
        head, _, body = chunk.partition(b"\n")
        out.append((head.split()[0].decode(), body.replace(b"\n", b"")))
    return out


def parse_gff_cds(lines, drop_flagged: bool) -> dict[str, set[tuple[int, int, int]]]:
    """CDS features as {seqid: {(start, end, strand)}}, 0-based half-open, strand +1/-1. `drop_flagged` drops
    pseudo and partial features (RefSeq annotates them in column 9; Prodigal's `partial=01` fields are kept)."""
    genes: dict[str, set[tuple[int, int, int]]] = defaultdict(set)
    for ln in lines:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        if drop_flagged and ("pseudo=true" in f[8] or "partial=true" in f[8]):
            continue
        genes[f[0]].add((int(f[3]) - 1, int(f[4]), 1 if f[6] == "+" else -1))
    return genes


def stop_key(gene: tuple[int, int, int], n: int) -> tuple[int, int]:
    """(3' end coordinate, strand): the end for + strand genes, the start for - strand genes. A + strand gene
    that wraps the origin (end > n) has its stop at end - n."""
    s, e, st = gene
    return (e % n if e > n else e, st) if st == 1 else (s, st)


def compare_by_stop(a: set, b: set, n: int) -> dict[str, int]:
    """a = bio-box genes, b = reference genes on a contig of length n."""
    sa, sb = {stop_key(g, n) for g in a}, {stop_key(g, n) for g in b}
    return {"biobox": len(a), "reference": len(b), "exact": len(a & b), "same_stop": len(sa & sb),
            "only_biobox_by_stop": len(sa - sb), "only_reference_by_stop": len(sb - sa)}


def near_end(gene: tuple[int, int, int], n: int, margin: int = EDGE_BP) -> bool:
    """True if the gene wraps the origin or lies within `margin` bp of either end of a contig of length n."""
    s, e, _ = gene
    return e > n or s < margin or e > n - margin


def prodigal_problems(biobox: set, prodigal: set, n: int, circular: bool) -> list[str]:
    """Rule for one contig. Linear: the gene sets are equal. Circular: every differing gene is near an end or
    wraps the origin, and at most MAX_DIFF_PER_CIRCULAR_CONTIG genes differ in total."""
    only_bb, only_pr = sorted(biobox - prodigal), sorted(prodigal - biobox)
    if not circular:
        return [f"linear contig differs: only bio-box {only_bb[:5]}, only Prodigal {only_pr[:5]}"] if only_bb or only_pr else []
    problems = [f"difference away from the ends: {g}" for g in only_bb + only_pr if not near_end(g, n)]
    if len(only_bb) + len(only_pr) > MAX_DIFF_PER_CIRCULAR_CONTIG:
        problems.append(f"{len(only_bb) + len(only_pr)} genes differ (limit {MAX_DIFF_PER_CIRCULAR_CONTIG}): "
                        f"only bio-box {only_bb[:5]}, only Prodigal {only_pr[:5]}")
    return problems


# ---------------------------------------------------------------- shared fixtures


@pytest.fixture(scope="module")
def genes(store_root):
    """Every gene of the golden store: dicts with accession, contig_uid, start, end, strand, rd (default raw_digest),
    length (contig length). Read from protein_occurrence, as step6_export.py did."""
    store = SeqStore(store_root)
    occ = pq.read_table(
        store_root / f"l1/protein_occurrence/gene_call_version={GENE_CALL_VERSION}/stage=S0/part-0.parquet").to_pylist()
    co = pq.read_table(next((store_root / "l0seq/contig_occurrence").rglob("part-0.parquet"))).to_pylist()
    acc = {r["contig_uid"]: r["accession"] for r in co}
    occ.sort(key=lambda r: (acc[r["contig_uid"]], r["start"], r["end"], r["strand"]))  # coordinate order (B-17)
    out = []
    for r in occ:
        rd = store.default_digest(r["contig_uid"])
        out.append({"accession": acc[r["contig_uid"]], "contig_uid": r["contig_uid"], "start": r["start"],
                    "end": r["end"], "strand": r["strand"], "rd": rd, "length": store.length_of(rd)})
    return out


def _shards(store_root: Path) -> list[Path]:
    return sorted(store_root.glob("l0seq/seq_shards/stage=*/shard-*.fa.gz"))


def _run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, check=True, capture_output=True, **kw)


# ---------------------------------------------------------------- golden-set tests


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not on PATH (it is in the image)")
@pytest.mark.covers("R-03")
def test_every_gene_slice_equals_samtools_faidx(store_root, genes, tmp_path):
    """Every non-origin-wrapping gene: bio-box's bytes for (contig, start, end, strand) equal `samtools faidx`
    on the shard with the repo-written .fai/.gzi (`-i` for the minus strand). Measured 2026-11-09: 33,713 of 33,713."""
    store = SeqStore(store_root)
    groups: dict[tuple[str, int], list[dict]] = defaultdict(list)
    for g in genes:
        if g["end"] <= g["length"]:
            groups[(store.shard_of(g["rd"]), g["strand"])].append(g)
    checked = 0
    for (shard, strand), gs in sorted(groups.items()):
        regions = tmp_path / "regions.txt"
        regions.write_text("".join(region(g["rd"].hex(), g["start"], g["end"]) + "\n" for g in gs))
        cmd = ["samtools", "faidx", str(store_root / shard), "-r", str(regions)]
        if strand == -1:
            cmd += ["-i", "--mark-strand", "no"]
        got = parse_fasta_bytes(_run(cmd).stdout)
        assert len(got) == len(gs), (shard, strand)
        for g, (name, seq) in zip(gs, got):
            assert name == region(g["rd"].hex(), g["start"], g["end"]), (name, g)
            want = store.get(g["contig_uid"], g["start"], g["end"], g["strand"])
            assert seq == want, (shard, strand, g["accession"], g["start"], g["end"])
            checked += 1
    assert checked >= 1
    print(f"\nsamtools slices: {checked} of {checked} identical")


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not on PATH (it is in the image)")
@pytest.mark.covers("R-07")
def test_samtools_builds_identical_indexes(store_root, tmp_path):
    """samtools faidx on a copy of each shard without its indexes writes .fai and .gzi byte-identical to the
    repo's (measured 2026-11-09, samtools 1.19.2: identical)."""
    shards = _shards(store_root)
    assert shards
    for shard in shards:
        copy = tmp_path / shard.parent.name / shard.name
        copy.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(shard, copy)
        _run(["samtools", "faidx", str(copy)])
        for suffix in (".fai", ".gzi"):
            assert Path(str(copy) + suffix).read_bytes() == Path(str(shard) + suffix).read_bytes(), (shard.name, suffix)


@pytest.mark.skipif(shutil.which("bgzip") is None, reason="bgzip not on PATH (it ships with samtools/htslib in the image)")
@pytest.mark.covers("R-07")
def test_bgzip_accepts_every_shard(store_root):
    """`bgzip -t` passes on every shard (measured 2026-11-09: passes)."""
    shards = _shards(store_root)
    assert shards
    for shard in shards:
        _run(["bgzip", "-t", str(shard)])


@pytest.mark.skipif(
    shutil.which("prodigal") is None,
    reason="prodigal not on PATH: the Prodigal binary is not in container/environment.yml yet, so this skips in CI "
           "until the PI adds prodigal=2.6.3 to the image")
@pytest.mark.covers("G-01")
def test_gene_calls_equal_prodigal_binary(genes, tmp_path):
    """bio-box vs the Prodigal 2.6.3 binary (`-g CODE -f gff -q`, single mode) on each golden prokaryote.
    Linear contigs: gene sets (start, end, strand) are equal. Circular contigs: bio-box calls genes across the
    origin and Prodigal calls fragments at the contig ends, so every difference must be an origin-wrapping gene or
    lie within 50 kb of an end, and at most 6 genes may differ per circular contig (bio-box-only plus Prodigal-only; the largest measured was 4). Measured 2026-11-09: 32,195 of
    32,198 bio-box genes exact across ten genomes; the 3 differences were at the ends of three circular contigs."""
    by_acc: dict[str, set] = defaultdict(set)
    length: dict[str, int] = {}
    for g in genes:
        by_acc[g["accession"]].add((g["start"], g["end"], g["strand"]))
        length[g["accession"]] = g["length"]
    problems: list[str] = []
    n_contigs = 0
    for e in GOLDEN["genomes"]:
        if e["kind"] != "prokaryote":
            continue
        fna = tmp_path / f"{e['id']}.fna"
        with gzip.open(Path(_DIR) / LOCK["files"][e["id"]]["path"], "rb") as fi, open(fna, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        gff = tmp_path / f"{e['id']}.gff"
        _run(["prodigal", "-i", str(fna), "-f", "gff", "-o", str(gff), "-g", str(e["genetic_code"]), "-q"])
        ref = parse_gff_cds(gff.read_text().splitlines(), drop_flagged=False)
        circular = {r["accession"] for r in e["records"] if r["topology"] == "circular"}
        for r in e["records"]:
            acc = r["accession"]
            n_contigs += 1
            problems += [f"{e['id']} {acc}: {p}"
                         for p in prodigal_problems(by_acc[acc], ref.get(acc, set()), length[acc], acc in circular)]
    assert n_contigs >= 10
    assert problems == []


def _refseq_gff() -> Path:
    if not _DIR:
        pytest.skip("golden set not fetched (set BIOBOX_GOLDEN_DIR)")
    path = Path(_DIR) / REFSEQ_GFF
    if not path.exists():
        pytest.skip(f"RefSeq truth file {path} absent: place GCF_000005845.2_ASM584v2_genomic.gff.gz under "
                    "$BIOBOX_GOLDEN_DIR/truth/ (no download in tests; its home will be an L0 source, card E-03)")
    return path


@pytest.mark.covers("G-01")
def test_refseq_cds_for_ecoli_k12(genes):
    """RefSeq CDS (non-pseudo, non-partial) of NC_000913.3 in GCF_000005845.2: a RefSeq CDS is found when a bio-box
    gene has the same stop coordinate and strand. Measured 2026-11-09: 4,119 of 4,303 found by stop (0.957) and
    3,889 of those with the same start (0.944); thresholds 0.95 and 0.93. The file is pinned by sha256."""
    path = _refseq_gff()
    assert hashlib.sha256(path.read_bytes()).hexdigest() == REFSEQ_GFF_SHA256
    with gzip.open(path, "rt") as fh:
        ref = parse_gff_cds(fh, drop_flagged=True)["NC_000913.3"]
    bb = {(g["start"], g["end"], g["strand"]) for g in genes if g["accession"] == "NC_000913.3"}
    n = next(g["length"] for g in genes if g["accession"] == "NC_000913.3")
    c = compare_by_stop(bb, ref, n)
    sensitivity, start_given_stop = c["same_stop"] / c["reference"], c["exact"] / c["same_stop"]
    print(f"\nRefSeq K-12: {c}, sensitivity by stop {sensitivity:.4f}, same start given stop {start_given_stop:.4f}")
    assert c["reference"] > 4000
    assert sensitivity >= 0.95
    assert start_given_stop >= 0.93


# ---------------------------------------------------------------- unit tests: no binary, no golden set


class TestUnit:
    def test_region_is_one_based_inclusive(self):
        assert region("ab", 0, 3) == "ab:1-3"
        assert region("ab", 99, 100) == "ab:100-100"

    def test_parse_fasta_bytes_joins_wrapped_lines(self):
        data = b">r1:1-6 extra\nACG\nTTT\n>r2:5-6\nGA\n"
        assert parse_fasta_bytes(data) == [("r1:1-6", b"ACGTTT"), ("r2:5-6", b"GA")]
        assert parse_fasta_bytes(b"") == []

    def test_parse_gff_cds(self):
        lines = [
            "##gff-version 3",
            "c1\tRefSeq\tgene\t1\t9\t.\t+\t.\tID=g1",
            "c1\tRefSeq\tCDS\t1\t9\t.\t+\t0\tID=a",
            "c1\tRefSeq\tCDS\t20\t40\t.\t-\t0\tID=b;pseudo=true",
            "c1\tRefSeq\tCDS\t50\t70\t.\t-\t0\tID=c;partial=true",
            "c2\tProdigal\tCDS\t5\t10\t.\t-\t0\tID=d;partial=01",
        ]
        assert parse_gff_cds(lines, drop_flagged=True) == {"c1": {(0, 9, 1)}, "c2": {(4, 10, -1)}}
        assert parse_gff_cds(lines, drop_flagged=False)["c1"] == {(0, 9, 1), (19, 40, -1), (49, 70, -1)}

    def test_stop_key(self):
        assert stop_key((10, 40, 1), 100) == (40, 1)
        assert stop_key((10, 40, -1), 100) == (10, -1)
        assert stop_key((90, 130, 1), 100) == (30, 1)  # wraps the origin: stop at end - n

    def test_compare_by_stop(self):
        bb = {(0, 30, 1), (40, 70, -1), (80, 99, 1)}
        ref = {(3, 30, 1), (40, 70, -1), (10, 20, 1)}
        c = compare_by_stop(bb, ref, 100)
        assert c == {"biobox": 3, "reference": 3, "exact": 1, "same_stop": 2,
                     "only_biobox_by_stop": 1, "only_reference_by_stop": 1}

    def test_near_end(self):
        n = 1_000_000
        assert near_end((900, 1200, 1), n)
        assert near_end((990_000, 1_000_050, 1), n)  # wraps
        assert near_end((960_000, 970_000, -1), n)  # within 50 kb of the end
        assert not near_end((500_000, 501_000, 1), n)
        assert not near_end((50_000, 51_000, 1), n)  # start exactly at the margin is not near

    def test_linear_contigs_must_be_equal(self):
        a = {(0, 30, 1), (40, 70, -1)}
        assert prodigal_problems(a, set(a), 1000, circular=False) == []
        assert prodigal_problems(a, a - {(0, 30, 1)}, 1000, circular=False)
        assert prodigal_problems(a, a | {(80, 90, 1)}, 1000, circular=False)

    def test_circular_tolerates_wrapping_gene_and_edge_fragment(self):
        n = 1_000_000
        common = {(500_000, 501_000, 1), (600_000, 601_000, -1)}
        wrap = (999_000, 1_000_300, 1)
        frag = (100, 700, 1)
        assert prodigal_problems(common | {wrap}, common | {frag}, n, circular=True) == []
        assert prodigal_problems(common | {wrap}, common, n, circular=True) == []

    def test_circular_rejects_interior_difference(self):
        n = 1_000_000
        common = {(500_000, 501_000, 1)}
        interior = (300_000, 301_000, -1)
        problems = prodigal_problems(common | {interior}, common, n, circular=True)
        assert any("away from the ends" in p for p in problems)

    def test_circular_rejects_too_many_differences(self):
        n = 1_000_000
        # all near the end, but one more than the limit
        extra = {(100 + 400 * i, 400 + 400 * i, 1) for i in range(MAX_DIFF_PER_CIRCULAR_CONTIG + 1)}
        problems = prodigal_problems(extra, set(), n, circular=True)
        assert len(problems) == 1 and f"limit {MAX_DIFF_PER_CIRCULAR_CONTIG}" in problems[0]
        assert prodigal_problems({(100, 400, 1), (500, 800, 1)}, set(), n, circular=True) == []
