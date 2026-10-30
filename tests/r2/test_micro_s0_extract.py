"""Ring R2: O-EXTRACT (bgzip FASTA + .fai + .gzi) on micro-S0.

Covers R-07 (extract round trip: re-ingested digests equal, samtools accepts the indexes),
R-03 against `samtools faidx` (with `-i` for minus strand) where samtools is on PATH,
LY-08 (byte-identical rewrite), and the block counting used by K-13.
"""

from __future__ import annotations

import random
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

from biobox import fasta, ids
from biobox.l0seq import SeqStore, SourceObject, ingest, revcomp
from biobox.l0seq.extract import BgzfFasta, extract

pytestmark = pytest.mark.r2

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
SAMTOOLS = shutil.which("samtools")


def _seqs():
    return [fasta.read(FIX / o["file"])[0].seq for o in CONTRACT["objects"]]


@pytest.fixture(scope="module")
def exported(tmp_path_factory):
    d = tmp_path_factory.mktemp("extract")
    root = d / "store"
    ingest([SourceObject(FIX / o["file"], "s", "1") for o in CONTRACT["objects"]], root, "S0", chunk_size=4096, row_group_chunks=16)
    store = SeqStore(root)
    contigs = [ids.contig_uid(s) for s in _seqs()]
    out = d / "out.fa.gz"
    report = extract(store, contigs, out)
    return d, store, contigs, out, report


def test_r07_reingest_gives_identical_digests(exported, tmp_path):
    _, _, contigs, out, report = exported
    recs = fasta.read(out)  # BGZF is multi-member gzip
    assert [r.accession for r in recs] == [c.hex() for c in contigs]
    assert [ids.raw_digest(r.seq) for r in recs] == [ids.raw_digest(s) for s in _seqs()]
    assert report["records"] == 3 and report["blocks"] > 3


@pytest.mark.covers("LY-08")
def test_rewrite_is_byte_identical(exported, tmp_path):
    _, store, contigs, out, _ = exported
    again = tmp_path / "again.fa.gz"
    extract(store, contigs, again)
    for suffix in ("", ".fai", ".gzi"):
        assert Path(str(again) + suffix).read_bytes() == Path(str(out) + suffix).read_bytes(), suffix


def test_index_slices_equal_store_and_count_blocks(exported):
    _, store, contigs, out, _ = exported
    bg = BgzfFasta(out)
    rng = random.Random(7)
    for cu in contigs:
        n = len(store.get(cu))
        for _ in range(100):
            a, b = sorted(rng.randrange(0, n + 1) for _ in range(2))
            before = bg.blocks_read
            assert bg.fetch(cu.hex(), a, b) == store.get(cu, a, b)
            assert bg.blocks_read - before == bg.units_for(cu.hex(), a, b)
    with pytest.raises(ValueError):
        bg.fetch(contigs[0].hex(), 0, 10**9)


@pytest.mark.skipif(SAMTOOLS is None, reason="samtools not on PATH")
def test_r07_samtools_builds_the_same_indexes(exported, tmp_path):
    _, _, _, out, _ = exported
    copy = tmp_path / "copy.fa.gz"
    shutil.copy(out, copy)
    subprocess.run([SAMTOOLS, "faidx", str(copy)], check=True)
    assert Path(str(copy) + ".fai").read_bytes() == Path(str(out) + ".fai").read_bytes()
    assert Path(str(copy) + ".gzi").read_bytes() == Path(str(out) + ".gzi").read_bytes()


@pytest.mark.skipif(SAMTOOLS is None, reason="samtools not on PATH")
def test_r03_slices_equal_samtools_faidx_both_strands(exported):
    _, store, contigs, out, _ = exported
    rng = random.Random(11)
    for cu in contigs:
        n = len(store.get(cu))
        regions = []
        for _ in range(40):
            a, b = sorted(rng.randrange(0, n) for _ in range(2))
            regions.append((a, b + 1))
        regions += [(0, 1), (n - 1, n), (0, n)]
        for a, b in regions:
            reg = f"{cu.hex()}:{a + 1}-{b}"  # samtools regions are 1-based inclusive
            for strand, flag in ((1, []), (-1, ["-i"])):
                res = subprocess.run([SAMTOOLS, "faidx", *flag, str(out), reg], check=True, capture_output=True)
                got = b"".join(res.stdout.splitlines()[1:])
                want = store.get(cu, a, b, strand)
                assert got == want, (reg, strand)
                if strand == -1:
                    assert got == revcomp(store.get(cu, a, b))
