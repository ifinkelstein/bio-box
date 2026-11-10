"""Ring R2: query-time readers take the release's stages (card B-28; AGENTS §3.9).

The store: micro-S0 as stage S0 with release R1 over it (test_micro_s0_l3's build, plus IGRs), then a
second stage S1 that plants what a search of R1 would find if it read S1: two copies of the lambda
genome, one with a codon of protein P changed (P', a protein only S1 has and that a search seeded with
P hits) and one with a codon of another lambda protein changed (so P itself sits on an S1 contig, where
a neighborhood fetch of P would find it). S1 is gene-called, scanned, given IGRs and taped, so every
reader that lists `stage=*` directories finds it. Release R2 is defined over S0 and S1: it is the
control that shows the planted genes are reachable.

Readers of R1 (the search pack, a search run's hits and neighborhoods, `nbr`, `SeqStore` and the
catalog's sequence functions) must return nothing of S1. Separately, two releases over one stage
(rS0.b28a, rS0.b28b) each cluster, each with its own ledger under `ledger/release=<R>/`, and both
manifests verify.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from biobox import catalog, fasta, ids, ledger
from biobox.l0seq import SeqStore, SourceObject, ingest, revcomp
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l1 import domains, igr
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l1.genecall import run as call_genes
from biobox.l3 import tapes
from biobox.search import fastz, nbr, pack, runner

sys.path.insert(0, str(Path(__file__).parent))
import test_micro_s0_l3 as micro  # noqa: E402

pytestmark = pytest.mark.r2

FIX = micro.FIX
R1 = micro.REL  # over S0 only (micro.build defines it)
R2 = "rS1.b28"  # over S0 and S1: the control
MMSEQS = shutil.which("mmseqs")
LIMITS = runner.Limits(cpus=2, wall_s=120)
ALA = {b"GCT", b"GCC", b"GCA", b"GCG"}


def _mutate(seq: bytes, start: int, end: int, strand: int) -> bytes:
    """`seq` with the middle codon of the gene (start, end, strand) changed to another amino acid
    (alanine, or tryptophan where it was alanine): one substitution, never a stop."""
    k = (end - start) // 3 // 2
    pos = start + 3 * k if strand == 1 else end - 3 * (k + 1)
    codon = seq[pos:pos + 3] if strand == 1 else revcomp(seq[pos:pos + 3])
    new = b"TGG" if codon.upper() in ALA else b"GCT"
    return seq[:pos] + (new if strand == 1 else revcomp(new)) + seq[pos + 3:]


def _occurrences(root: Path, stage: str) -> list[dict]:
    d = root / f"l1/protein_occurrence/gene_call_version={V}/stage={stage}"
    return [r for f in sorted(d.glob("part-*.parquet")) for r in pq.read_table(f, partitioning=None).to_pylist()]


@pytest.fixture(scope="module")
def store(tmp_path_factory):
    d = tmp_path_factory.mktemp("b28")
    root = d / "a"
    circ = micro.build(root)
    igr.run(root, "S0", circ)
    lam = fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq
    lam_cu = ids.contig_uid(lam)
    prots = domains.load_proteins(root, V, "S0")
    genes = sorted((r for r in _occurrences(root, "S0") if r["contig_uid"] == lam_cu), key=lambda r: r["start"])
    p = max(genes, key=lambda r: (len(prots[r["protein_uid"]]), r["protein_uid"]))
    q = max(genes, key=lambda r: (abs(r["start"] - p["start"]), r["protein_uid"]))  # far from P
    seq_a = _mutate(lam, p["start"], p["end"], p["strand"])  # P -> P'
    seq_b = _mutate(lam, q["start"], q["end"], q["strand"])  # Q -> Q'; P unchanged
    src = d / "s1.fna"
    src.write_bytes(b"\n".join(fasta.format_record("lambda_p_mut", seq_a, 70) + fasta.format_record("lambda_q_mut", seq_b, 70))
                    + b"\n")
    ingest([SourceObject(src, "b28_planted", "1", "test", True)], root, "S1")
    cu_a, cu_b = ids.contig_uid(seq_a), ids.contig_uid(seq_b)
    write_genome_input(root, "S1", [GenomeInput((cu,), 11, kind="virus") for cu in (cu_a, cu_b)])
    call_genes(root, "S1")
    domains.run(root, "S1", V, "test", "1", d / "a.hmm", "gathering")
    igr.run(root, "S1", set())
    tapes.run(root, "S1", set())
    catalog.define_release(root, R2, ["S0", "S1"])
    s1 = _occurrences(root, "S1")
    p_mut = {r["protein_uid"] for r in s1 if r["contig_uid"] == cu_a and (r["start"], r["end"]) == (p["start"], p["end"])}
    assert len(p_mut) == 1 and p["protein_uid"] not in p_mut, "the planted P' was not called as one new protein"
    assert any(r["protein_uid"] == p["protein_uid"] and r["contig_uid"] == cu_b for r in s1), "P is not on the S1 contig"
    models = {"test-1": d / "a.hmm"}
    return {"root": root, "d": d, "models": models, "lam": lam_cu, "s0_contigs": {r["contig_uid"] for r in _occurrences(root, "S0")},
            "s1_contigs": {cu_a, cu_b}, "cu_a": cu_a, "cu_b": cu_b, "p": p["protein_uid"], "p_mut": next(iter(p_mut)),
            "s0_proteins": set(prots), "s1_only": {r["protein_uid"] for r in s1} - set(prots),
            "packs": {rel: pack.build(root, rel, models=models) for rel in (R1, R2)}}


def _pack_ids(root: Path, release: str, target: str) -> set[str]:
    man, _ = pack.read_manifest(root, release)
    return {n for s in man["shards"] if s["target"] == target
            for n, _ in fastz.read_records(pack.pack_dir(root, release) / s["path"])}


def _meta_contigs(root: Path, release: str, target: str) -> set[bytes]:
    t = pq.read_table(pack.pack_dir(root, release) / f"{target}/{target}_meta.parquet", columns=["contig_uid"])
    return set(t.column("contig_uid").to_pylist())


# --- the search pack ------------------------------------------------------------------------------


@pytest.mark.covers("C-12")
def test_b28_pack_of_r1_holds_only_r1_stages(store):
    root = store["root"]
    assert store["s1_only"] and store["p_mut"] in store["s1_only"]
    # R1: exactly the S0 proteins, and regions of S0 contigs only
    assert _pack_ids(root, R1, "unique") == {u.hex() for u in store["s0_proteins"]}
    for target in ("igr", "dark_orf"):
        got = _meta_contigs(root, R1, target)
        assert got and got <= store["s0_contigs"], target
    assert pack.verify(root, R1) == []
    # the control: R2's pack has the planted proteins and S1 regions
    assert _pack_ids(root, R2, "unique") == {u.hex() for u in store["s0_proteins"] | store["s1_only"]}
    assert _meta_contigs(root, R2, "igr") & store["s1_contigs"]
    assert pack.verify(root, R2) == []
    # the pack's inputs are the release's files: no S1 path is hashed into R1's sources
    m1, m2 = store["packs"][R1], store["packs"][R2]
    assert m1["targets"]["unique"]["sources"] != m2["targets"]["unique"]["sources"]


@pytest.mark.covers("C-12", "LY-10")
def test_b28_pack_from_the_manifest_equals_the_pack_from_release_def(store, tmp_path):
    root = tmp_path / "a"
    shutil.copytree(store["root"], root)
    man = catalog.build_manifest(root, R1, catalog.release_stages(root, R1))
    catalog.write_manifest(man, catalog.manifest_path(root, R1))
    files = catalog.ReleaseFiles.of(root, R1)
    assert files.source == "manifest" and files.stages == ["S0"]
    by_def = catalog.ReleaseFiles(root, R1, ["S0"])
    for name in ("protein_seq", "protein_occurrence", "feature_protein", "igr", "contig_tapes", "l1_domains",
                 "seq_locator", "contig_index", "contig_occurrence"):
        assert files.table(name) == by_def.table(name), name
        assert files.table(name) and not [p for p in files.table(name) if "stage=S1" in p.as_posix()], name
    assert files.artifact("seq_shards") == by_def.artifact("seq_shards") != []
    again = pack.build(root, R1, models=store["models"])
    assert again["pack_id"] == store["packs"][R1]["pack_id"]


# --- SeqStore and the catalog's sequence functions ----------------------------------------------------


@pytest.mark.covers("LY-10")
def test_b28_seqstore_for_r1_cannot_fetch_an_s1_contig_and_says_so(store):
    root = store["root"]
    for stages in (catalog.ReleaseFiles.of(root, R1), ["S0"]):
        s = SeqStore(root, stages)
        try:
            assert s.locator.footers_read == 1 and s.index.footers_read == 1  # S0's files only, footers only
            assert s.locator.row_groups_read == 0 and s.index.row_groups_read == 0
            assert len(s.get(store["lam"], 0, 10)) == 10
            with pytest.raises(KeyError, match=r"none of the stages \['S0'\]"):
                s.get(store["cu_a"], 0, 10)
        finally:
            s.close()
    both = SeqStore(root, catalog.ReleaseFiles.of(root, R2))
    try:
        assert len(both.get(store["cu_a"], 0, 10)) == 10  # the control: R2 reads S1
    finally:
        both.close()
    with pytest.raises(TypeError, match="stages"):
        SeqStore(root)  # no default: a store never lists every stage on disk
    with pytest.raises(TypeError, match="stages"):
        SeqStore(root, 64)


@pytest.mark.covers("LY-10")
def test_b28_catalog_of_r1_reads_no_s1_sequence(store, tmp_path):
    root = store["root"]
    con = catalog.open_catalog(catalog.generate_catalog(catalog.build_manifest(root, R1, ["S0"]), tmp_path / "c.duckdb"))
    try:
        assert catalog.seq_reader(con).stages == ["S0"]
        assert len(con.execute("SELECT seq_get(?, 0, 12, 1)", [store["lam"]]).fetchone()[0]) == 12
        with pytest.raises(Exception, match=r"unknown contig_uid .*\['S0'\]"):
            con.execute("SELECT seq_get(?, 0, 12, 1)", [store["cu_a"]]).fetchall()
    finally:
        con.close()


# --- nbr and a search run -----------------------------------------------------------------------------


@pytest.mark.covers("S-09")
def test_b28_nbr_of_r1_returns_nothing_from_s1(store):
    root, policy = store["root"], nbr.Policy.parse(window="3g")
    rows, counts = nbr.fetch(root, {store["p"]}, policy, files=catalog.ReleaseFiles.of(root, R1))
    assert rows and {r["contig_uid"] for r in rows} == {store["lam"]} and counts["anchors"] == 1
    occ = nbr.occurrences(root, {store["p"]}, files=catalog.ReleaseFiles.of(root, R1))
    assert {o["contig_uid"] for o in occ} == {store["lam"]}
    s1 = nbr.summarize(root, {store["p"]}, policy, release=R1, top_n=0)
    assert s1["header"]["n_contigs"] == 1 and s1["header"]["n_anchor_loci"] == 1
    # the control: R2 reads S1, where P sits on a second contig
    rows2, _ = nbr.fetch(root, {store["p"]}, policy, files=catalog.ReleaseFiles.of(root, R2))
    assert {r["contig_uid"] for r in rows2} == {store["lam"], store["cu_b"]}
    assert nbr.summarize(root, {store["p"]}, policy, release=R2, top_n=0)["header"]["n_contigs"] == 2
    with pytest.raises(TypeError, match="ReleaseFiles"):
        nbr.fetch(root, {store["p"]}, policy, files=None)


def _recipe(release: str, seed: bytes) -> str:
    return (f"recipe: b28_phmmer\nwritten_against: {release}\nseed: [{{protein_uid: {seed.hex()}}}]\n"
            "target: unique\ntool: phmmer\niterations: {max: 1, inclusion_bits: 25}\n"
            f"controls: {{positive: [{{protein_uid: {seed.hex()}}}]}}\ncontext: {{window: 3g, edge: truncate}}\n")


def _search(store, release: str) -> tuple[set[bytes], list[dict]]:
    root = store["root"]
    row = runner.run(root, _recipe(release, store["p"]), release, work_base=store["d"] / "work", limits=LIMITS)
    assert row["status"] == "ok", row.get("error")
    hits = {h["target_uid"] for h in runner.read_hits(root, row["run_id"].hex())}
    nb = pq.read_table(runner.run_dir(root, row["run_id"]) / "neighborhoods.parquet").to_pylist()
    return hits, nb


@pytest.mark.covers("S-09")
def test_b28_search_run_against_r1_returns_nothing_from_s1(store):
    hits, nb = _search(store, R1)
    assert store["p"] in hits and store["p_mut"] not in hits
    assert hits <= store["s0_proteins"]
    assert nb and {r["contig_uid"] for r in nb} <= store["s0_contigs"]
    assert not {r["contig_uid"] for r in nb} & store["s1_contigs"]
    # the control: the same search against R2 finds P' and the S1 copy of P
    hits2, nb2 = _search(store, R2)
    assert store["p_mut"] in hits2
    assert {store["cu_a"], store["cu_b"]} <= {r["contig_uid"] for r in nb2}


# --- two releases over one stage both cluster ---------------------------------------------------------


@pytest.mark.skipif(MMSEQS is None, reason="needs MMseqs2 (runs in the image)")
@pytest.mark.covers("LY-10", "L-03")
def test_b28_two_releases_over_one_stage_both_cluster_and_verify(store, tmp_path):
    from biobox.l2 import cluster

    root = tmp_path / "a"
    shutil.copytree(store["root"], root)
    rels = ("rS0.b28a", "rS0.b28b")  # one stage prefix: before B-28 both ledgers were ledger/stage=S0/step=cluster
    work = tmp_path / "mm"
    work.mkdir()
    got = [cluster.run(root, rel, V, ["S0"], tmp=work) for rel in rels]
    assert got[0] == got[1]
    for rel in rels:
        assert ledger.release_path(root, rel, "cluster").is_file()
        assert {r["object_key"] for r in ledger.read_release(root, rel, "cluster")} == {b"nr90", b"nr50"}
        assert ledger.check_chain(root, "S0", release=rel) == []
        assert ledger.check_balance(root, "S0", ["genecall", "cluster"], release=rel) == []
    assert ledger.read_release(root, rels[0], "cluster") == ledger.read_release(root, rels[1], "cluster")
    assert not (root / "ledger/stage=S0/step=cluster").exists()
    mans = {rel: catalog.build_manifest(root, rel, ["S0"]) for rel in rels}
    for rel, man in mans.items():
        own = [f for f in man["files"] if f["table"] == "ledger" and f["partitions"].get("release") == rel]
        assert [f["path"] for f in own] == [f"ledger/release={rel}/step=cluster/part-0.parquet"]
        assert own[0]["partitions"] == {"release": rel, "step": "cluster"}
        assert not [f for f in man["files"] if f["partitions"].get("stage") == "S1"]
        assert catalog.verify_manifest(man, root) == []
    # the catalog's ledger view holds stage and release ledgers side by side
    con = catalog.open_catalog(catalog.generate_catalog(mans[rels[0]], tmp_path / "c.duckdb"))
    try:
        steps = set(con.execute("SELECT stage, release, step FROM ledger GROUP BY ALL").fetchall())
        assert ("S0", None, "genecall") in steps and (None, rels[0], "cluster") in steps
    finally:
        con.close()
