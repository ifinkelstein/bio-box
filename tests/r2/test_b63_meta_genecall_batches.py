"""Ring R2: byte-identity of gene calling with interleaved short and long plasmid/viral/prokaryotic
genomes across worker counts and batch sizes (card B-63; D-01).

Runs with no external tools (only pyrodigal / pyrodigal-gv), so it passes when MinCED, MMseqs2,
hmmsearch, and samtools are absent (GitHub Actions compatibility).
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from biobox import fasta, ids
from biobox.l0seq import SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l1 import genecall

pytestmark = [pytest.mark.r2, pytest.mark.covers("B-63", "D-01")]

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"


def _sha_tree(root: Path) -> dict[str, str]:
    return {
        p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def test_b63_interleaved_short_and_long_plasmids_are_byte_identical_across_workers(tmp_path):
    """Interleaved <10 kb and >=10 kb linear & circular plasmid contigs (plus a viral and a single-mode
    prokaryote genome) produce byte-identical Parquet tables and ledger across workers=1 vs workers=4
    and different `batch_bases`."""
    m = fasta.read(FIX / "NC_000908.2.fna.gz")[0].seq
    t4 = fasta.read(FIX / "NC_000866.4.fna.gz")[0].seq
    lam = fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq

    # Build 20 short (3 kb) and 10 long (15 kb) distinct plasmid fragments, interleaved
    seqs: dict[str, bytes] = {"M": m, "L": lam}
    genomes: list[GenomeInput] = []
    for i in range(10):
        s_short_a = t4[i * 3_200 : i * 3_200 + 3_000]
        s_short_b = m[i * 3_500 : i * 3_500 + 3_200]
        s_long = t4[40_000 + i * 12_000 : 40_000 + i * 12_000 + 15_000]
        for tag, s, circ in (
            (f"sa{i}", s_short_a, i % 3 == 0),
            (f"lg{i}", s_long, i % 4 == 0),
            (f"sb{i}", s_short_b, False),
        ):
            seqs[tag] = s
            cu = ids.contig_uid(s)
            genomes.append(GenomeInput((cu,), 11, kind="plasmid", circular=(cu,) if circ else ()))
    # Also include a viral (gv-meta) and a >=100 kb prokaryote (single mode) genome
    genomes.append(GenomeInput((ids.contig_uid(lam),), 11, kind="virus"))
    genomes.append(GenomeInput((ids.contig_uid(m),), 4, kind="prokaryote", circular=(ids.contig_uid(m),)))

    fna = tmp_path / "interleaved.fna"
    fna.write_bytes(b"".join(b"\n".join(fasta.format_record(k, v, 80)) + b"\n" for k, v in seqs.items()))

    root = tmp_path / "root"
    write_genome_input(root, "PL3T", genomes)
    ingest([SourceObject(fna, "b63_test", "1", "INSDC (public)", True)], root, "PL3T")

    out1 = tmp_path / "out_w1"
    res1 = genecall.run(root, "PL3T", out=out1, tmp=tmp_path, workers=1, batch_bases=20_000)
    hashes1 = _sha_tree(out1)
    assert len(hashes1) == 5  # feature, feature_protein, protein_occurrence, protein_seq, ledger
    assert set(res1["modes"]) == {"gv-meta", "meta", "single"}

    for tag, kw in (
        ("w4_b20k", {"workers": 4, "batch_bases": 20_000}),
        ("w4_b1m", {"workers": 4, "batch_bases": 1_000_000}),
        ("w2_b10k", {"workers": 2, "batch_bases": 10_000}),
    ):
        out_k = tmp_path / f"out_{tag}"
        res_k = genecall.run(root, "PL3T", out=out_k, tmp=tmp_path, **kw)
        assert res_k == res1, tag
        assert _sha_tree(out_k) == hashes1, tag
