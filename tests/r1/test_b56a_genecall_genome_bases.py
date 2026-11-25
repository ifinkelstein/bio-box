"""Ring R1, card B-56a: `KeyLookup.find_first` in `genecall._genome_bases` (O(1) stage lookup).

Tests:
1. On a 10-stage synthetic store (`S0`..`S9`, opened with `stages[0]` owning the queried contig),
   `genecall._genome_bases(store, lengths, genome)` returns the exact base count and reads
   `lengths.row_groups_read == 1` (instead of scanning every earlier stage's `contig_index`).
2. `genecall._genome_bases` still resolves contigs owned by the last stage `stages[-1]` (and
   multi-contig genomes spanning `stages[0]` and `stages[-1]`), and raises the same `KeyError`
   when a contig is absent from `store.stages` (both nonexistent and owned by an unopened stage).
3. `genecall.run` on a 10-stage store (`seq_stages` spanning all earlier stages) produces
   100% byte-identical output files and result dict compared to the single-stage run.
"""

from __future__ import annotations

import hashlib
import random
from pathlib import Path

import pytest

from biobox import fasta, ids
from biobox.l0seq import SeqStore, SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput
from biobox.l0seq.lookup import KeyLookup
from biobox.l1 import genecall

pytestmark = pytest.mark.covers("B-56a", "B-66")

REPO = Path(__file__).resolve().parents[2]


def _seq(rng: random.Random, n: int) -> bytes:
    return rng.randbytes(n).translate(bytes(b"ACGT"[i % 4] for i in range(256)))


def _ten_stage_store(tmp_path: Path) -> tuple[Path, list[str], dict[str, list[tuple[bytes, bytes]]]]:
    """Build a 10-stage store (`S0`..`S9`), each stage holding 8 distinct contigs so every stage's
    `contig_index` spans hash space with 1 row group each. Stage `S9` includes a real coding slice
    from micro_s0 T4 so `genecall.run` produces non-empty gene calls."""
    rng = random.Random(5601)
    root = tmp_path / "root"
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    t4_path = REPO / "tests" / "fixtures" / "micro_s0" / "NC_000866.4.fna.gz"
    t4_slice = fasta.read(t4_path)[0].seq[:25_000]

    stage_names = [f"S{i}" for i in range(10)]
    by_stage: dict[str, list[tuple[bytes, bytes]]] = {}
    for st in stage_names:
        seqs = [_seq(rng, 500 + 41 * k) for k in range(8)]
        if st == "S9":
            seqs[0] = t4_slice
        fna = in_dir / f"{st}.fna"
        fna.write_bytes(
            b"\n".join(ln for k, s in enumerate(seqs) for ln in fasta.format_record(f"{st}_c{k}", s, 80)) + b"\n"
        )
        ingest([SourceObject(fna, "b56a", "1", "CC0", True)], root, st)
        by_stage[st] = [(ids.contig_uid(s), s) for s in seqs]

    # Current stage first (`stages[0] == "S9"`), followed by earlier stages (`stages[-1] == "S8"`),
    # matching `genecall.run`'s `[stage, *seq_stages]` ordering.
    stages = ["S9", "S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"]
    return root, stages, by_stage


def _tree_hashes(root: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for p in sorted(root.rglob("*")):
        if p.is_file():
            out[p.relative_to(root).as_posix()] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out


def test_b56a_genome_bases_reads_single_row_group_and_handles_earlier_and_missing(tmp_path):
    """`_genome_bases` reads 1 `contig_index` row group for a contig owned by `stages[0]`,
    resolves contigs owned by `stages[-1]`, and raises the exact `KeyError` when absent."""
    root, stages, by_stage = _ten_stage_store(tmp_path)
    cu_first, seq_first = by_stage[stages[0]][0]
    cu_last, seq_last = by_stage[stages[-1]][0]

    store = SeqStore(root, stages)
    try:
        lengths = KeyLookup(store.index_files, "contig_uid", ["contig_len", "is_default_variant"])
        g_first = GenomeInput((cu_first,), 11, kind="prokaryote", circular=())
        assert genecall._genome_bases(store, lengths, g_first) == len(seq_first)
        assert lengths.row_groups_read == 1, (
            f"expected 1 contig_index RG read for contig in stages[0], got {lengths.row_groups_read}"
        )

        # Contig owned by the last stage `stages[-1]` resolves with exact base count
        g_last = GenomeInput((cu_last,), 11, kind="prokaryote", circular=())
        assert genecall._genome_bases(store, lengths, g_last) == len(seq_last)

        # Multi-contig genome spanning stages[0] and stages[-1] sums both lengths
        g_both = GenomeInput((cu_first, cu_last), 11, kind="prokaryote", circular=())
        assert genecall._genome_bases(store, lengths, g_both) == len(seq_first) + len(seq_last)

        # Missing contig not on disk anywhere raises KeyError
        missing_cu = b"\xab" * 16
        g_missing = GenomeInput((missing_cu,), 11, kind="prokaryote", circular=())
        with pytest.raises(KeyError, match=f"contig_uid {missing_cu.hex()} is in none of the stages"):
            genecall._genome_bases(store, lengths, g_missing)
    finally:
        store.close()

    # Contig present on disk in an unopened stage (`stages[-1]`) names that stage (DEC-26)
    store_subset = SeqStore(root, stages[:-1])
    try:
        lengths_subset = KeyLookup(store_subset.index_files, "contig_uid", ["contig_len", "is_default_variant"])
        with pytest.raises(KeyError, match=stages[-1]):
            genecall._genome_bases(store_subset, lengths_subset, GenomeInput((cu_last,), 11, kind="prokaryote", circular=()))
    finally:
        store_subset.close()


def test_b56a_genecall_run_byte_identical_across_multistage_store(tmp_path):
    """`genecall.run` with 9 earlier `seq_stages` produces byte-identical outputs and identical
    result dict compared to running with only the owning stage."""
    root, stages, by_stage = _ten_stage_store(tmp_path)
    cu0, _ = by_stage["S9"][0]
    cu1, _ = by_stage["S9"][1]
    genomes = [
        GenomeInput((cu0, cu1), 11, kind="prokaryote", circular=()),
    ]

    out_multi = tmp_path / "out_multi"
    out_single = tmp_path / "out_single"
    tmp_multi = tmp_path / "tmp_multi"
    tmp_single = tmp_path / "tmp_single"
    tmp_multi.mkdir()
    tmp_single.mkdir()

    res_multi = genecall.run(
        root, "S9", genomes=genomes, seq_stages=stages[1:], out=out_multi, tmp=tmp_multi, workers=1
    )
    res_single = genecall.run(
        root, "S9", genomes=genomes, seq_stages=None, out=out_single, tmp=tmp_single, workers=1
    )

    assert res_multi == res_single
    assert res_multi["genes"] > 0
    hashes_multi = _tree_hashes(out_multi)
    hashes_single = _tree_hashes(out_single)
    assert hashes_multi, "expected published genecall output files"
    assert hashes_multi == hashes_single
