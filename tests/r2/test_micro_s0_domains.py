"""Ring R2: O-SCAN (pyhmmer, fixed Z) on the micro-S0 fixture.

The HMM library is built at test time from single lambda proteins (pyhmmer Builder, fixed
seed), with gathering cutoffs set, so no external database is needed. Covers C-01, a
positive control (each model finds its source protein), the `hmmsearch` binary as an
independent implementation, MR-12 (searching a subset gives identical scores at fixed Z),
and MR-16 (cpus 1 vs 4 give identical bytes).
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pyarrow.parquet as pq
import pyhmmer
import pytest
import yaml

from biobox import fasta, ids
from biobox.l0seq import SourceObject, ingest
from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, run as call_genes
from biobox.registry.schemas import conformance_problems, table

pytestmark = pytest.mark.r2

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
N_MODELS = 15


def _build_store(root: Path) -> None:
    ingest([SourceObject(FIX / o["file"], "fixture", "1") for o in CONTRACT["objects"]], root, "S0")
    genomes = []
    for o in CONTRACT["objects"]:
        cu = ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)
        genomes.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"],
                                   circular=(cu,) if o["topology"] == "circular" else ()))
    call_genes(root, "S0", genomes)


def _lambda_proteins(root: Path) -> dict[bytes, str]:
    occ = pq.read_table(root / f"l1/protein_occurrence/gene_call_version={GENE_CALL_VERSION}/stage=S0/part-0.parquet").to_pylist()
    lam = ids.contig_uid(fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq)
    allp = domains.load_proteins(root, GENE_CALL_VERSION, "S0")
    return {o["protein_uid"]: allp[o["protein_uid"]] for o in occ if o["contig_uid"] == lam}


def _make_library(prots: dict[bytes, str], path: Path) -> dict[str, bytes]:
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet, seed=42)
    background = pyhmmer.plan7.Background(alphabet)
    chosen = sorted(prots, key=lambda u: (-len(prots[u]), u))[:N_MODELS]
    source = {}
    with open(path, "wb") as fh:
        for i, uid in enumerate(chosen):
            seq = pyhmmer.easel.TextSequence(name=f"src{i}".encode(), sequence=prots[uid]).digitize(alphabet)
            hmm, _, _ = builder.build(seq, background)
            hmm.name = f"lam{i:02d}".encode()
            hmm.accession = f"TEST{i:04d}".encode()
            hmm.cutoffs.gathering = (25.0, 25.0)
            hmm.write(fh)
            source[f"TEST{i:04d}"] = uid
    return source


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = tmp_path_factory.mktemp("domains")
    _build_store(root)
    lib = root / "lib.hmm"
    source = _make_library(_lambda_proteins(root), lib)
    report = domains.run(root, "S0", GENE_CALL_VERSION, "test", "1", lib, "gathering", cpus=1)
    rows = pq.read_table(root / "l1/domains/library=test-1/stage=S0/part-0.parquet")
    return root, lib, source, report, rows


@pytest.mark.covers("C-01", "C-06", "LY-01")
def test_table_conforms(built):
    _, _, _, _, rows = built
    assert rows.num_rows > 0
    assert conformance_problems(rows, table("l1_domains")) == []


def test_each_model_finds_its_source_protein(built):
    root, _, source, _, rows = built
    prots = domains.load_proteins(root, GENE_CALL_VERSION, "S0")
    by = {}
    for r in rows.to_pylist():
        by.setdefault((r["model_acc"], r["protein_uid"]), []).append(r)
    for acc, uid in source.items():
        hits = by.get((acc, uid), [])
        assert hits, acc
        covered = max(h["env_end"] - h["env_start"] for h in hits)
        assert covered >= 0.9 * len(prots[uid]), (acc, covered, len(prots[uid]))


def test_fixed_z_is_recorded(built):
    _, _, _, _, rows = built
    assert set(rows.column("z").to_pylist()) == {1e8}
    assert set(rows.column("dom_z").to_pylist()) == {1e8}


def test_mr12_subset_search_gives_identical_scores(built):
    root, lib, _, _, rows = built
    lam = _lambda_proteins(root)
    sub = domains.scan(lam, lib, "gathering", cpus=1)
    full = [r for r in rows.to_pylist() if r["protein_uid"] in lam]
    assert sub == full


def test_mr16_cpus_do_not_change_output(built, tmp_path):
    root, lib, _, _, _ = built
    out1 = root / "l1/domains/library=test-1/stage=S0/part-0.parquet"
    a = out1.read_bytes()
    # into a fresh root holding the same proteins: the writers never rewrite a file (card B-19)
    fresh = tmp_path / "s"
    shutil.copytree(root / "l0seq/protein_seq", fresh / "l0seq/protein_seq")
    domains.run(fresh, "S0", GENE_CALL_VERSION, "test", "1", lib, "gathering", cpus=4)
    assert (fresh / "l1/domains/library=test-1/stage=S0/part-0.parquet").read_bytes() == a


@pytest.mark.skipif(shutil.which("hmmsearch") is None, reason="hmmsearch binary not on PATH")
def test_matches_hmmsearch_binary(built, tmp_path):
    root, lib, _, _, rows = built
    prots = domains.load_proteins(root, GENE_CALL_VERSION, "S0")
    faa = tmp_path / "p.faa"
    faa.write_text("".join(f">{u.hex()}\n{s.replace('*', '')}\n" for u, s in sorted(prots.items())))
    dom = tmp_path / "dom.tsv"
    subprocess.run(["hmmsearch", "--cut_ga", "-Z", "1e8", "--domZ", "1e8", "--cpu", "1", "--noali",
                    "--domtblout", str(dom), "-o", str(tmp_path / "out.txt"), str(lib), str(faa)], check=True)
    ref = {}
    for line in dom.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split()
        # target, tacc, tlen, query, qacc, qlen, E, score, bias, #, of, c-E, i-E, dscore, dbias, hmmfrom, hmmto, alifrom, alito, envfrom, envto
        ref[(f[0], f[4], int(f[19]) - 1, int(f[20]))] = (float(f[7]), float(f[13]), float(f[12]))
    ours = {(r["protein_uid"].hex(), r["model_acc"], r["env_start"], r["env_end"]):
            (r["seq_bitscore"], r["dom_bitscore"], r["dom_ievalue_zref"]) for r in rows.to_pylist()}
    assert set(ours) == set(ref)
    for k, (s, d, e) in ours.items():
        rs, rd, re_ = ref[k]
        assert s == pytest.approx(rs, abs=0.05) and d == pytest.approx(rd, abs=0.05)
        assert e == pytest.approx(re_, rel=0.02)


# --- card B-24: the scan in shards ----------------------------------------------------------------------
#
# `domains.plan` cuts the stage's proteins, sorted by protein_uid, into shards by cumulative residues;
# `domains.scan_shard` scans one shard in its own call into a sorted scratch run; `domains.merge` writes
# the runs as the stage's domain files and ledger. The oracle is the unsharded `domains.run` of the
# `built` fixture: every file must be byte-identical, for any shard size.

OUTS = ("l1/domains/library=test-1/stage=S0", "ledger/stage=S0/step=domains-test-1")


def _outputs(root: Path) -> dict[str, bytes]:
    return {p.relative_to(root).as_posix(): p.read_bytes() for d in OUTS for p in sorted((root / d).rglob("*"))
            if p.is_file()}


def _fresh(root: Path, at: Path) -> Path:
    """A new root holding only the fixture's proteins (the writers never rewrite a file, card B-19)."""
    shutil.copytree(root / "l0seq/protein_seq", at / "l0seq/protein_seq")
    return at


def _lengths(root: Path) -> list[tuple[bytes, int]]:
    """(protein_uid, residues) of every distinct protein, in uid order, from the sequences themselves."""
    return sorted((u, len(s)) for u, s in domains.load_proteins(root, GENE_CALL_VERSION, "S0").items())


def _sharded(root: Path, lib: Path, size: int, tmp: Path, cpus: int = 1) -> dict:
    plan = domains.plan(root, "S0", GENE_CALL_VERSION, "test", "1", lib, "gathering", size)
    for k in range(plan["shards"]):
        domains.scan_shard(root, "S0", "test", "1", lib, k, cpus=cpus, tmp=tmp)
    domains.merge(root, "S0", "test", "1", tmp=tmp)
    return plan


@pytest.mark.covers("D-01", "MR-12")
@pytest.mark.parametrize("want", ["one", "two", "many"])
def test_b24_sharded_scan_equals_the_unsharded_run(built, tmp_path, want):
    root, lib = built[0], built[1]
    lengths = _lengths(root)
    total, longest = sum(n for _, n in lengths), max(n for _, n in lengths)
    size = {"one": total, "two": -(-total // 2) + longest, "many": total // 8}[want]
    fresh = _fresh(root, tmp_path / "s")
    _sharded(fresh, lib, size, tmp_path / "scratch", cpus=2 if want == "many" else 1)
    assert _outputs(fresh) == _outputs(root)
    # the plan: shards in uid order cover every protein once, each within `size` unless it is one protein,
    # and each closed only when its next protein would not fit (cut by cumulative residues)
    shards = domains.read_plan(fresh, "S0", "test", "1")["shards"]
    if want == "one":
        assert len(shards) == 1
    elif want == "two":
        assert len(shards) == 2
    else:
        assert len(shards) >= 5
    at = 0
    for s in shards:
        assert s["first"] == lengths[at][0].hex()
        took = lengths[at:at + s["proteins"]]
        assert s["last"] == took[-1][0].hex() and s["residues"] == sum(n for _, n in took)
        assert s["residues"] <= size or s["proteins"] == 1
        at += s["proteins"]
        if at < len(lengths):
            assert s["residues"] + lengths[at][1] > size
    assert at == len(lengths)


def test_b24_a_shard_refuses_another_library_file(built, tmp_path):
    root, lib = built[0], built[1]
    fresh = _fresh(root, tmp_path / "s")
    domains.plan(fresh, "S0", GENE_CALL_VERSION, "test", "1", lib, "gathering", 10**9)
    other = tmp_path / "other.hmm"
    other.write_bytes(lib.read_bytes().split(b"//\n", 1)[1])  # every model but the first
    with pytest.raises(ValueError, match="sha256"):
        domains.scan_shard(fresh, "S0", "test", "1", other, 0, cpus=1, tmp=tmp_path / "scratch")
    with pytest.raises(FileExistsError):  # a plan is written once
        domains.plan(fresh, "S0", GENE_CALL_VERSION, "test", "1", lib, "gathering", 10**9)


def test_b24_merge_refuses_a_missing_shard_or_a_changed_run(built, tmp_path):
    root, lib = built[0], built[1]
    fresh = _fresh(root, tmp_path / "s")
    scratch = tmp_path / "scratch"
    plan = domains.plan(fresh, "S0", GENE_CALL_VERSION, "test", "1", lib, "gathering", sum(n for _, n in _lengths(root)) // 3)
    assert plan["shards"] >= 3
    for k in range(plan["shards"] - 1):
        domains.scan_shard(fresh, "S0", "test", "1", lib, k, cpus=1, tmp=scratch)
    with pytest.raises(FileNotFoundError, match=f"shard {plan['shards'] - 1}"):
        domains.merge(fresh, "S0", "test", "1", tmp=scratch)
    domains.scan_shard(fresh, "S0", "test", "1", lib, plan["shards"] - 1, cpus=1, tmp=scratch)
    with pytest.raises(FileExistsError):  # a shard is scanned once
        domains.scan_shard(fresh, "S0", "test", "1", lib, 0, cpus=1, tmp=scratch)
    run0 = domains.run_path(scratch, fresh, "S0", "test", "1", 0)
    run0.write_bytes(run0.read_bytes() + b"\0")
    with pytest.raises(ValueError, match="sha256"):
        domains.merge(fresh, "S0", "test", "1", tmp=scratch)
    assert _outputs(fresh) == {}  # nothing was written
