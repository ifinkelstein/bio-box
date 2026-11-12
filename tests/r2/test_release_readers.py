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

Card B-33: a file under `release=<R>` is R's only. R1's manifest is built, then a later release writes its
files (rS0.b33 over S0 alone, and R2 over S0 and S1): R1's manifest still verifies and is rebuilt
byte-identical, the later manifest lists no R1 file, a stray under `release=R1` is reported by R1's
manifest and not by the later one's, and a `tape_version_def` is in a manifest only when the release's
files use that tape version.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from biobox import catalog, contracts, dim, fasta, ids, ledger
from biobox.l0seq import SeqStore, SourceObject, ingest, revcomp
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l1 import domains, igr
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l1.genecall import run as call_genes
from biobox.l3 import postings, tapes, tokens
from biobox.registry.schemas import TABLES
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
    # card W-03: the packs and the runs live beside the data root, never under it
    pack_dir, runs = d / "a.packs", d / "a.runs"
    return {"root": root, "d": d, "models": models, "lam": lam_cu, "s0_contigs": {r["contig_uid"] for r in _occurrences(root, "S0")},
            "s1_contigs": {cu_a, cu_b}, "cu_a": cu_a, "cu_b": cu_b, "p": p["protein_uid"], "p_mut": next(iter(p_mut)),
            "s0_proteins": set(prots), "s1_only": {r["protein_uid"] for r in s1} - set(prots),
            "packs": {rel: pack.build(root, rel, models=models, packs=pack_dir) for rel in (R1, R2)},
            "pack_dir": pack_dir, "runs": runs}


def at(store) -> dict:
    """The runner's keyword arguments that name the store's pack and runs directories (card W-03)."""
    return {"packs": store["pack_dir"], "runs": store["runs"]}


def _pack_ids(store, release: str, target: str) -> set[str]:
    pdir = pack.find(store["pack_dir"], release)
    man, _ = pack.read_manifest(pdir)
    return {n for s in man["shards"] if s["target"] == target for n, _ in fastz.read_records(pdir / s["path"])}


def _meta_contigs(store, release: str, target: str) -> set[bytes]:
    t = pq.read_table(pack.find(store["pack_dir"], release) / f"{target}/{target}_meta.parquet", columns=["contig_uid"])
    return set(t.column("contig_uid").to_pylist())


# --- the search pack ------------------------------------------------------------------------------


@pytest.mark.covers("C-12")
def test_b28_pack_of_r1_holds_only_r1_stages(store):
    root = store["root"]
    assert store["s1_only"] and store["p_mut"] in store["s1_only"]
    # R1: exactly the S0 proteins, and regions of S0 contigs only
    assert _pack_ids(store, R1, "unique") == {u.hex() for u in store["s0_proteins"]}
    for target in ("igr", "dark_orf"):
        got = _meta_contigs(store, R1, target)
        assert got and got <= store["s0_contigs"], target
    assert pack.verify(root, R1, packs=store["pack_dir"]) == []
    # the control: R2's pack has the planted proteins and S1 regions
    assert _pack_ids(store, R2, "unique") == {u.hex() for u in store["s0_proteins"] | store["s1_only"]}
    assert _meta_contigs(store, R2, "igr") & store["s1_contigs"]
    assert pack.verify(root, R2, packs=store["pack_dir"]) == []
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
    again = pack.build(root, R1, models=store["models"], packs=tmp_path / "packs")  # card W-03: its own pack dir
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
    row = runner.run(root, _recipe(release, store["p"]), release, work_base=store["d"] / "work", limits=LIMITS,
                     **at(store))
    assert row["status"] == "ok", row.get("error")
    d = runner.run_dir(store["runs"], row["run_id"])
    hits = {h["target_uid"] for h in runner.read_hits(d)}
    nb = pq.read_table(d / "neighborhoods.parquet").to_pylist()
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


# --- each release's manifest covers that release's files only (card B-33) -------------------------------

R2_SAME = "rS0.b33"  # a second release over S0 alone
# member_cluster: card B-27; genome_dim: card C-21b
OWN_TABLES = {"release_def", "token_dict", "tape_tokens", "token_postings", "contig_dim", "genome_dim",
              "cluster_member", "member_cluster"} | ({"cluster_member_aln", "ledger"} if MMSEQS else set())


def _release_files(root: Path, release: str, stages: list[str], work: Path) -> None:
    """Every release-level builder of `release` over `stages` (as test_micro_s0_release_def builds one)."""
    if MMSEQS:
        from biobox.l2 import cluster

        cluster.run(root, release, V, stages, tmp=work)
    else:
        micro._stand_in_clusters(root, release)  # S0's proteins only: the S0+S1 case needs MMseqs2
    tokens.run(root, release, stages)
    postings.run(root, release, stages)
    dim.run(root, release, stages, tmp=work)


def _s0_genomes() -> list[GenomeInput]:
    out = []
    for o in micro.CONTRACT["objects"]:
        cu = ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)
        out.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"],
                               circular=(cu,) if o["topology"] == "circular" else ()))
    return out


@pytest.fixture(scope="module", params=[(R2_SAME, ("S0",)), (R2, ("S0", "S1"))], ids=["same-stage", "more-stages"])
def later(store, tmp_path_factory, request):
    """R1's manifest, built before a later release over the same stage (or over S0 and S1) writes its files."""
    rel2, stages2 = request.param[0], list(request.param[1])
    if "S1" in stages2 and not MMSEQS:
        pytest.skip("the stand-in clusters hold S0's proteins only; S0+S1 needs MMseqs2 (runs in the image)")
    d = tmp_path_factory.mktemp("b33")
    root, work = d / "a", d / "work"
    shutil.copytree(store["root"], root)
    work.mkdir()
    write_genome_input(root, "S0", _s0_genomes())  # micro.build passes them to the gene caller directly
    dim.run(root, R1, ["S0"], tmp=work)  # R1's contig_dim beside its clusters, tokens and postings
    m1 = catalog.build_manifest(root, R1, ["S0"])
    catalog.define_release(root, rel2, stages2)  # R2 is already defined in the store: a no-op there
    _release_files(root, rel2, stages2, work)
    m2 = catalog.build_manifest(root, rel2, stages2)
    return {"root": root, "rel2": rel2, "stages2": stages2, "m1": m1, "m2": m2}


def _copy(later, tmp_path: Path) -> Path:
    c = tmp_path / "c"
    shutil.copytree(later["root"], c)
    return c


@pytest.mark.covers("LY-10", "D-01")
def test_b33_an_earlier_manifest_still_verifies_and_rebuilds_identically_after_a_later_release(later):
    root, m1, m2, rel2 = later["root"], later["m1"], later["m2"], later["rel2"]
    for m, rel in ((m1, R1), (m2, rel2)):  # both releases wrote every release-level table
        assert {f["table"] for f in m["files"] if f["partitions"].get("release") == rel} == OWN_TABLES, rel
    assert catalog.verify_manifest(m1, root) == []
    assert catalog.manifest_bytes(catalog.build_manifest(root, R1, ["S0"])) == catalog.manifest_bytes(m1)
    assert catalog.verify_manifest(m2, root) == []
    assert catalog.manifest_bytes(catalog.build_manifest(root, rel2, later["stages2"])) == catalog.manifest_bytes(m2)


@pytest.mark.covers("LY-10")
def test_b33_a_later_manifest_lists_no_file_of_an_earlier_release(later):
    m1, m2, rel2 = later["m1"], later["m2"], later["rel2"]
    assert not [f["path"] for f in m2["files"] + m2["artifacts"] if f"release={R1}/" in f["path"]]
    assert {f["partitions"]["release"] for f in m2["files"] if "release" in f["partitions"]} == {rel2}
    assert {f["partitions"]["release"] for f in m1["files"] if "release" in f["partitions"]} == {R1}
    # the stage files: the later release holds every file of R1's stage with the same bytes
    def s0(m):
        return {f["path"]: f["sha256"] for f in m["files"] + m["artifacts"] if f["partitions"].get("stage") == "S0"}
    assert s0(m1) and s0(m1) == s0(m2)
    assert bool([f for f in m2["files"] if f["partitions"].get("stage") == "S1"]) == ("S1" in later["stages2"])
    assert not [f for f in m1["files"] if f["partitions"].get("stage") == "S1"]
    # the unpartitioned artifact and the one tape version both use are in both
    tv_defs = [f"l3/tape_version_def/tape_version={tapes.tape_version([V])}/part-0.parquet"]
    for m in (m1, m2):
        assert [f["path"] for f in m["files"] if f["table"] == "tape_version_def"] == tv_defs
        assert [a["path"] for a in m["artifacts"] if a["artifact"] == "store_params"] == ["l0seq/store_params.json"]


@pytest.mark.covers("LY-10")
def test_b33_a_stray_of_a_release_is_reported_by_that_release_only(later, tmp_path):
    c = _copy(later, tmp_path)
    m1, m2, rel2 = later["m1"], later["m2"], later["rel2"]
    for rel, table_dir in ((R1, "l3/token_dict"), (rel2, "dim/contig_dim")):
        src = c / f"{table_dir}/release={rel}/part-0.parquet"
        shutil.copy(src, src.with_name("part-9.parquet"))
    stray1, stray2 = f"l3/token_dict/release={R1}/part-9.parquet", f"dim/contig_dim/release={rel2}/part-9.parquet"
    assert catalog.verify_manifest(m1, c) == [f"not in manifest: {stray1}"]
    assert catalog.verify_manifest(m2, c) == [f"not in manifest: {stray2}"]
    # a stray in a stage both releases list is reported by both
    src = c / "l0seq/contig_index/stage=S0/part-0.parquet"
    shutil.copy(src, src.with_name("part-9.parquet"))
    stray0 = "l0seq/contig_index/stage=S0/part-9.parquet"
    assert catalog.verify_manifest(m1, c) == [f"not in manifest: {p}" for p in sorted([stray0, stray1])]
    assert catalog.verify_manifest(m2, c) == [f"not in manifest: {p}" for p in sorted([stray0, stray2])]


@pytest.mark.covers("LY-10")
def test_b33_a_tape_version_def_is_in_the_manifests_of_releases_that_use_it(later, tmp_path):
    c = _copy(later, tmp_path)
    tv = tapes.tape_version([V])
    other = tapes.ensure_def(c, ["b33-unused-producer"])  # a tape version no file of either release names
    assert other != tv and tapes.def_path(c, other).is_file()
    for rel, stages, m in ((R1, ["S0"], later["m1"]), (later["rel2"], later["stages2"], later["m2"])):
        assert catalog.verify_manifest(m, c) == []
        again = catalog.build_manifest(c, rel, stages)
        assert catalog.manifest_hash(again) == catalog.manifest_hash(m), rel
    # a stray beside the definition of the version the releases use is reported
    stray = tapes.def_path(c, tv).with_name("part-1.parquet")
    shutil.copy(tapes.def_path(c, tv), stray)
    for m in (later["m1"], later["m2"]):
        assert catalog.verify_manifest(m, c) == [f"not in manifest: {stray.relative_to(c).as_posix()}"]


@pytest.mark.covers("LY-10")
def test_b33_release_files_from_release_def_list_only_the_release_own_files(later):
    root = later["root"]
    for rel, stages in ((R1, ["S0"]), (later["rel2"], later["stages2"])):
        by_def = catalog.ReleaseFiles(root, rel, stages)
        by_man = catalog.ReleaseFiles.of(root, rel, catalog.build_manifest(root, rel, stages))
        for name in sorted(OWN_TABLES - {"ledger"}):
            got = by_def.table(name)
            assert got and all(f"release={rel}/" in p.as_posix() for p in got), (rel, name)
            assert got == by_man.table(name), (rel, name)
        assert by_def.table("ledger", release=rel) == by_man.table("ledger", release=rel)


@pytest.mark.covers("LY-10")
def test_b33_a_manifest_with_no_parquet_file_gives_a_catalog_with_empty_views(tmp_path):
    root = tmp_path / "empty"
    root.mkdir()
    man = catalog.build_manifest(root, "rS0.empty", [])
    assert man["files"] == [] and man["artifacts"] == []
    con = catalog.open_catalog(catalog.generate_catalog(man, tmp_path / "e.duckdb"))
    try:
        assert con.execute("SELECT count(*) FROM catalog_files").fetchone()[0] == 0
        views = {r[0] for r in con.execute("SELECT view_name FROM duckdb_views() WHERE NOT internal").fetchall()}
        placed = sorted(name for name, spec in TABLES.items() if spec.location)
        assert set(placed) <= views
        for name in placed:
            assert con.execute(f"SELECT count(*) FROM {name}").fetchone()[0] == 0, name
    finally:
        con.close()


@pytest.mark.covers("R-09", "D-06")
def test_b33_append_only_across_releases_reports_a_rewritten_stage_file_only(later, tmp_path):
    m1, m2, rel2 = later["m1"], later["m2"], later["rel2"]
    assert contracts.check_append_only(m1, m2) == []
    # a rewritten S0 file is still found: the same rows in other bytes
    c = _copy(later, tmp_path)
    rel = next(f["path"] for f in m1["files"] if f["table"] == "contig_index")
    p = c / rel
    old = p.read_bytes()
    codec = pq.ParquetFile(p).metadata.row_group(0).column(0).compression
    pq.write_table(pq.read_table(p, partitioning=None), p, compression="gzip" if codec == "ZSTD" else "zstd")
    assert p.read_bytes() != old
    assert contracts.check_append_only(m1, catalog.build_manifest(c, rel2, later["stages2"])) == [f"{rel}: rewritten"]
    # within one release, its own files are still compared
    gone = [f for f in m1["files"] if f["table"] == "token_dict"]
    assert gone
    assert contracts.check_append_only(m1, {**m1, "files": [f for f in m1["files"] if f not in gone]}) == [
        f"{f['path']}: removed" for f in gone]
