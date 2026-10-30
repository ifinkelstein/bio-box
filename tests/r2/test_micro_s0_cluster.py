"""Ring R2: O-CLUSTER (MMseqs2 NR90/NR50) on the micro-S0 fixture plus a mutated lambda.

Covers C-01 (registry contract), C-11 (one alignment per member per level),
an independent check of every stored alignment against the sequences, cluster
nesting, cluster_uid recomputation, and MR-16 (thread count does not change output).
Needs `mmseqs` on PATH (the bio-box image); skipped otherwise.
"""

from __future__ import annotations

import random
import re
import shutil
from pathlib import Path

import pyarrow.parquet as pq
import pytest
import yaml

from biobox import fasta, ids
from biobox.l0seq import SourceObject, ingest
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, run as call_genes
from biobox.l2.cluster import LEVELS, cluster_uid, params_hash, run as cluster
from biobox.registry.schemas import conformance_problems, table

pytestmark = [pytest.mark.r2, pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs not on PATH")]

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
RELEASE = "rS0.0-dev"


def _mutate(seq: bytes, rate: float, seed: int) -> bytes:
    rng = random.Random(seed)
    s = bytearray(seq)
    for i in range(len(s)):
        if rng.random() < rate:
            s[i] = rng.choice([b for b in b"ACGT" if b != s[i]])
    return bytes(s)


def _build(root: Path, threads: int) -> dict:
    lam = fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq
    src = root / "src"
    src.mkdir(parents=True, exist_ok=True)
    mut = src / "lambda_mut.fna"
    mut.write_bytes(b"\n".join(fasta.format_record("lambda_mut", _mutate(lam, 0.02, 11), 70)) + b"\n")
    objs = [SourceObject(FIX / o["file"], "fixture", "1") for o in CONTRACT["objects"]] + [SourceObject(mut, "mut", "1")]
    ingest(objs, root, "S0")
    genomes = []
    for o in CONTRACT["objects"]:
        cu = ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)
        genomes.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"],
                                   circular=(cu,) if o["topology"] == "circular" else ()))
    genomes.append(GenomeInput((ids.contig_uid(fasta.read(mut)[0].seq),), 11, kind="virus"))
    call_genes(root, "S0", genomes)
    return cluster(root, RELEASE, GENE_CALL_VERSION, threads=threads)


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = tmp_path_factory.mktemp("cluster")
    report = _build(root, threads=1)
    mem = pq.read_table(root / f"l2/cluster_member/release={RELEASE}/part-0.parquet")
    aln = pq.read_table(root / f"l2/cluster_member_aln/release={RELEASE}/part-0.parquet")
    prot = {}
    for f in root.glob(f"l0seq/protein_seq/gene_call_version={GENE_CALL_VERSION}/stage=*/protein_prefix=*/part-0.parquet"):
        t = pq.read_table(f)
        prot.update(zip(t.column("protein_uid").to_pylist(), t.column("seq").to_pylist()))
    return root, report, mem, aln, prot


def test_tables_conform(built):
    _, _, mem, aln, _ = built
    assert conformance_problems(mem, table("cluster_member")) == []
    assert conformance_problems(aln, table("cluster_member_aln")) == []


def test_every_protein_once_per_level_and_reps_are_members(built):
    _, _, mem, _, prot = built
    for level in LEVELS:
        rows = [r for r in mem.to_pylist() if r["level"] == level]
        assert sorted(r["member_uid"] for r in rows) == sorted(prot)
        reps = {r["rep_uid"] for r in rows}
        assert all(any(x["member_uid"] == rep and x["rep_uid"] == rep for x in rows) for rep in reps)


def test_nr90_nests_inside_nr50(built):
    _, _, mem, _, _ = built
    nr90 = {r["member_uid"]: r["cluster_uid"] for r in mem.to_pylist() if r["level"] == "nr90"}
    nr50 = {r["member_uid"]: r["cluster_uid"] for r in mem.to_pylist() if r["level"] == "nr50"}
    parent = {}
    for m, c90 in nr90.items():
        assert parent.setdefault(c90, nr50[m]) == nr50[m]


def test_cluster_uid_recomputes(built):
    _, _, mem, _, _ = built
    for level in LEVELS:
        rows = [r for r in mem.to_pylist() if r["level"] == level]
        by_c = {}
        for r in rows:
            by_c.setdefault(r["cluster_uid"], []).append(r["member_uid"])
        for c, ms in by_c.items():
            assert cluster_uid(ms, params_hash(level)) == c


def test_c11_one_alignment_per_member_per_level(built):
    _, _, mem, aln, _ = built
    for level in LEVELS:
        m = {(r["member_uid"], r["rep_uid"]) for r in mem.to_pylist() if r["level"] == level}
        a = [(r["member_uid"], r["rep_uid"]) for r in aln.to_pylist() if r["level"] == level]
        assert len(a) == len(set(a)) and set(a) == m
    for r in aln.to_pylist():
        if r["member_uid"] == r["rep_uid"]:
            assert r["pident"] == pytest.approx(1.0)


def _walk(cigar: str, rep: str, member: str, rs: int, ms: int):
    """Count identities by walking the CIGAR; M consumes both, I the representative, D the member."""
    ident = cols = 0
    i, j = rs, ms
    for n, op in re.findall(r"(\d+)([MID])", cigar):
        n = int(n)
        if op == "M":
            ident += sum(rep[i + k] == member[j + k] for k in range(n))
            i, j = i + n, j + n
        elif op == "I":
            i += n
        else:
            j += n
        cols += n
    return ident, cols, i, j


def test_stored_alignments_check_out_against_sequences(built):
    _, _, _, aln, prot = built
    for r in aln.to_pylist():
        rep, member = prot[r["rep_uid"]], prot[r["member_uid"]]
        ident, cols, i, j = _walk(r["cigar"].decode(), rep, member, r["tstart"], r["qstart"])
        assert (i, j) == (r["tend"], r["qend"]), "CIGAR spans must end at the stored coordinates"
        assert cols == r["alnlen"]
        assert ident / cols == pytest.approx(r["pident"], abs=0.006)


def test_members_meet_level_thresholds(built):
    _, _, _, aln, prot = built
    for level, p in LEVELS.items():
        rows = [r for r in aln.to_pylist() if r["level"] == level]
        ok = sum(
            r["pident"] >= p["min_seq_id"] - 0.005
            and (r["qend"] - r["qstart"]) / len(prot[r["member_uid"]]) >= p["coverage"] - 0.005
            for r in rows
        )
        assert ok / len(rows) >= 0.99, (level, ok, len(rows))


def test_mutated_lambda_proteins_join_lambda_clusters(built):
    root, _, mem, _, _ = built
    occ = pq.read_table(root / f"l1/protein_occurrence/gene_call_version={GENE_CALL_VERSION}/stage=S0/part-0.parquet").to_pylist()
    lam = ids.contig_uid(fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq)
    mut = ids.contig_uid(fasta.read(root / "src" / "lambda_mut.fna")[0].seq)
    nr90 = {r["member_uid"]: r["cluster_uid"] for r in mem.to_pylist() if r["level"] == "nr90"}
    lam_clusters = {nr90[o["protein_uid"]] for o in occ if o["contig_uid"] == lam}
    mut_prots = {o["protein_uid"] for o in occ if o["contig_uid"] == mut}
    shared = sum(nr90[p] in lam_clusters for p in mut_prots)
    assert shared / len(mut_prots) >= 0.8, (shared, len(mut_prots))


def test_mr16_thread_count_does_not_change_output(tmp_path, built):
    root1 = built[0]
    root4 = tmp_path / "t4"
    _build(root4, threads=4)
    for name in ("cluster_member", "cluster_member_aln"):
        rel = f"l2/{name}/release={RELEASE}/part-0.parquet"
        assert (root1 / rel).read_bytes() == (root4 / rel).read_bytes(), name
