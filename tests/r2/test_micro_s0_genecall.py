"""Ring R2: O-CALL on the micro-S0 fixture and the protein tables it writes.

Covers C-01/C-06 (registry contract), C-10 (protein_occurrence completeness),
R-04 (protein consistency against an independent translation), LY-08, MR-03
(reverse complement, with a measured tolerance), and a G-02-style genetic-code
check on M. genitalium (code 4).

Card B-05: the `feature` and `feature_protein` tables against `protein_occurrence` and a direct
`call_genome` (C-03), `feature_uid` recomputed from each row, reading `genome_input` when no genome
list is passed, and output bytes that depend neither on the genome order (MR-01) nor on `workers`
(D-01). Card B-17 removed the slot-keyed table these were first compared with.
"""

from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

from biobox import catalog, fasta, ids
from biobox.l0seq import SeqStore, SourceObject, ingest, revcomp
from biobox.l0seq.inputs import write_genome_input
from biobox.l1.genecall import GENE_CALL_VERSION, GenomeInput, call_genome, run
from biobox.registry.schemas import conformance_problems, table

pytestmark = pytest.mark.r2

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
OBJS = CONTRACT["objects"]

# NCBI codon tables in TCAG order (oracle independent of pyrodigal)
_BASES = "TCAG"
_AA = {
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    15: "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}


def _translate(nt: bytes, code: int) -> str:
    aa = _AA[code]
    out = []
    for i in range(0, len(nt) - 2, 3):
        c = nt[i : i + 3].decode().upper()
        if any(b not in _BASES for b in c):
            out.append("X")
        else:
            out.append(aa[16 * _BASES.index(c[0]) + 4 * _BASES.index(c[1]) + _BASES.index(c[2])])
    return "".join(out)


def _genome(o, seq=None):
    seq = seq if seq is not None else fasta.read(FIX / o["file"])[0].seq
    cu = ids.contig_uid(seq)
    circ = (cu,) if o["topology"] == "circular" else ()
    return GenomeInput((cu,), o["genetic_code"], kind=o["kind"], circular=circ)


def _genomes():
    return [_genome(o) for o in OBJS]


CIRCULAR = {ids.contig_uid(fasta.read(FIX / o["file"])[0].seq) for o in OBJS if o["topology"] == "circular"}


def _build(root: Path):
    objs = [SourceObject(FIX / o["file"], CONTRACT["source_name"], CONTRACT["source_release"]) for o in OBJS]
    ingest(objs, root, "S0")
    return run(root, "S0", _genomes())


def _read(root: Path, rel_glob: str) -> pa.Table:
    return pa.concat_tables(pq.read_table(f) for f in sorted(root.glob(rel_glob)))


@pytest.fixture(scope="module")
def root(tmp_path_factory):
    r = tmp_path_factory.mktemp("genecall")
    _build(r)
    return r


V = GENE_CALL_VERSION
GLOBS = {
    "protein_occurrence": f"l1/protein_occurrence/gene_call_version={V}/stage=S0/part-0.parquet",
    "protein_seq": f"l0seq/protein_seq/gene_call_version={V}/stage=S0/protein_prefix=*/part-*.parquet",
    "feature": f"l1/feature/producer={V}/stage=S0/part-0.parquet",
    "feature_protein": f"l1/feature_protein/gene_call_version={V}/stage=S0/part-0.parquet",
}
# every file gene calling writes for a stage (B-05): the tables above and its ledger
OUTPUTS = ("l1/**/*.parquet", "l0seq/protein_seq/**/*.parquet", "ledger/stage=S0/step=genecall/*.parquet")


def _outputs(root: Path) -> dict[str, bytes]:
    return {str(p.relative_to(root)): p.read_bytes() for g in OUTPUTS for p in sorted(root.glob(g))}


@pytest.mark.covers("C-01", "C-02", "C-06", "LY-01")
@pytest.mark.parametrize("name", sorted(GLOBS))
def test_tables_conform(root, name):
    files = sorted(root.glob(GLOBS[name]))
    assert files
    for f in files:
        assert conformance_problems(pq.read_table(f), table(name)) == [], f


def test_c10_occurrence_matches_feature_protein(root):
    cds = {r["feature_uid"]: r["protein_uid"] for r in _read(root, GLOBS["feature_protein"]).to_pylist()}
    occ = _read(root, GLOBS["protein_occurrence"]).to_pylist()
    occ_keys = [r["feature_uid"] for r in occ]
    assert len(occ_keys) == len(set(occ_keys)) == len(cds)
    assert {k: r["protein_uid"] for k, r in zip(occ_keys, occ)} == cds


def _direct_calls(root: Path) -> dict[tuple, dict]:
    """(contig_uid, start, end, strand) -> the gene as `call_genome` returns it, without the tables."""
    store = SeqStore(root, ["S0"])
    try:
        return {(c["contig_uid"], c["start"], c["end"], c["strand"]): c for g in _genomes() for c in call_genome(store, g)}
    finally:
        store.close()


def _legacy_genes(root: Path) -> dict[tuple, bytes]:
    """(contig_uid, start, end, strand) -> protein_uid from protein_occurrence (the slot-keyed table it was once
    joined with is gone, card B-17), checked against the genes `call_genome` returns directly."""
    out = {}
    for r in _read(root, GLOBS["protein_occurrence"]).to_pylist():
        out[(r["contig_uid"], r["start"], r["end"], r["strand"])] = r["protein_uid"]
    direct = _direct_calls(root)
    assert set(out) == set(direct) and len(out) == len(direct)
    assert out == {k: ids.protein_uid(c["aa"]) for k, c in direct.items()}
    return out


@pytest.mark.covers("C-01", "LY-10")
def test_new_tables_enter_the_manifest(root):
    man = catalog.build_manifest(root, "rS0.test", ["S0"])
    by_table = {f["table"]: f for f in man["files"]}
    for name in ("feature", "feature_protein"):
        assert by_table[name]["path"] == GLOBS[name]
        assert by_table[name]["rows"] == len(_legacy_genes(root))
    assert by_table["feature"]["partitions"] == {"producer": V, "stage": "S0"}


@pytest.mark.covers("C-03")
def test_features_equal_the_legacy_gene_set(root):
    legacy = _legacy_genes(root)
    feats = _read(root, GLOBS["feature"]).to_pylist()
    assert {(r["contig_uid"], r["start"], r["end"], r["strand"]) for r in feats} == set(legacy)
    assert len(feats) == len(legacy)
    for r in feats:
        assert r["kind"] == ids.KIND["cds"]
        assert r["feature_uid"] == ids.feature_uid(r["contig_uid"], r["start"], r["end"], r["strand"], r["kind"])


@pytest.mark.covers("C-03")
def test_feature_protein_equals_the_legacy_mapping(root):
    legacy = _legacy_genes(root)
    fp = _read(root, GLOBS["feature_protein"]).to_pylist()
    assert {(r["contig_uid"], r["start"], r["end"], r["strand"]): r["protein_uid"] for r in fp} == legacy
    feature_uids = {r["feature_uid"] for r in _read(root, GLOBS["feature"]).to_pylist()}
    assert {r["feature_uid"] for r in fp} == feature_uids
    direct = _direct_calls(root)
    assert {(r["contig_uid"], r["start"], r["end"], r["strand"]): r["genetic_code"] for r in fp} == \
        {k: c["genetic_code"] for k, c in direct.items()}
    assert {r["edit_type"] for r in fp} == {"none"}


@pytest.mark.covers("C-03")
def test_occurrence_feature_uid_recomputes_from_its_row(root):
    occ = _read(root, GLOBS["protein_occurrence"]).to_pylist()
    feature_uids = {r["feature_uid"] for r in _read(root, GLOBS["feature"]).to_pylist()}
    for r in occ:
        assert r["feature_uid"] == ids.feature_uid(r["contig_uid"], r["start"], r["end"], r["strand"], ids.KIND["cds"])
    assert {r["feature_uid"] for r in occ} == feature_uids


@pytest.mark.covers("LY-08")
def test_run_reads_genome_input_when_no_list_is_given(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    want = _build(a)
    objs = [SourceObject(FIX / o["file"], CONTRACT["source_name"], CONTRACT["source_release"]) for o in OBJS]
    ingest(objs, b, "S0")
    write_genome_input(b, "S0", _genomes())
    assert run(b, "S0") == want
    assert _outputs(a) == _outputs(b)
    assert len(_outputs(a)) == 5  # feature, feature_protein, protein_occurrence, protein_seq, ledger (B-17)


@pytest.mark.covers("MR-01", "D-01")
def test_mr01_bytes_do_not_depend_on_genome_order(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    want = _build(a)
    objs = [SourceObject(FIX / o["file"], CONTRACT["source_name"], CONTRACT["source_release"]) for o in OBJS]
    ingest(objs, b, "S0")
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    assert run(b, "S0", iter(_genomes()[::-1]), tmp=scratch) == want  # any iterable, consumed once
    assert _outputs(a) == _outputs(b)
    assert list(scratch.iterdir()) == []  # the scratch directory is removed


@pytest.mark.covers("D-01", "LY-08")
def test_bytes_do_not_depend_on_workers(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    want = _build(a)
    objs = [SourceObject(FIX / o["file"], CONTRACT["source_name"], CONTRACT["source_release"]) for o in OBJS]
    ingest(objs, b, "S0")
    assert run(b, "S0", _genomes(), workers=4) == want
    assert _outputs(a) == _outputs(b)


@pytest.mark.covers("D-01")
def test_worker_exception_propagates_and_publishes_nothing(tmp_path):
    from dataclasses import replace

    root, scratch = tmp_path / "s", tmp_path / "scratch"
    objs = [SourceObject(FIX / o["file"], CONTRACT["source_name"], CONTRACT["source_release"]) for o in OBJS]
    ingest(objs, root, "S0")
    genomes = _genomes()
    genomes[1] = replace(genomes[1], kind="bogus")  # the gene finder of this genome raises in its worker
    scratch.mkdir()
    with pytest.raises(ValueError, match="unknown genome kind 'bogus'"):
        run(root, "S0", genomes, workers=4, tmp=scratch)
    assert _outputs(root) == {}
    assert list(scratch.iterdir()) == []
    with pytest.raises(ValueError, match="under the data root"):
        run(root, "S0", _genomes(), tmp=root / "scratch")


@pytest.mark.covers("C-05")
def test_origin_wrapping_gene_has_end_past_contig_len_in_feature(tmp_path):
    # micro-S0 M. genitalium has no gene across its origin, so rotate it into one (as MR-04 does)
    o = next(o for o in OBJS if o["accession"] == "NC_000908.2")
    seq = fasta.read(FIX / o["file"])[0].seq
    n = len(seq)
    root = tmp_path / "s"
    ingest([SourceObject(FIX / o["file"], "a", "1")], root, "S0")
    base = call_genome(SeqStore(root, ["S0"]), _genome(o, seq))
    target = base[len(base) // 2]
    r = (target["start"] + target["end"]) // 2
    rot_file = tmp_path / "rot.fna"
    rot_file.write_bytes(b"\n".join(fasta.format_record("rot", seq[r:] + seq[:r], 70)) + b"\n")
    ingest([SourceObject(rot_file, "b", "1")], root, "S1")
    rotated = _genome(o, seq[r:] + seq[:r])
    report = run(root, "S1", [rotated])
    assert report["wrapping_genes"] >= 1
    cu = rotated.contig_uids[0]
    feats = _read(root, f"l1/feature/producer={V}/stage=S1/part-0.parquet").to_pylist()
    wraps = [f for f in feats if f["end"] > n]
    assert len(wraps) == report["wrapping_genes"]
    for f in wraps:
        assert f["contig_uid"] == cu and f["start"] < n and f["end"] - f["start"] < n
        assert f["feature_uid"] == ids.feature_uid(cu, f["start"], f["end"], f["strand"], ids.KIND["cds"])
    assert conformance_problems(_read(root, f"l1/feature/producer={V}/stage=S1/part-0.parquet"), table("feature")) == []


def test_protein_seq_is_content_addressed(root):
    prots = _read(root, GLOBS["protein_seq"]).to_pylist()
    for p in prots:
        assert p["protein_uid"] == ids.protein_uid(p["seq"])
        assert p["protein_digest"] == ids.raw_digest(p["seq"].encode())
        assert p["length"] == len(p["seq"])
    used = {r["protein_uid"] for r in _read(root, GLOBS["feature_protein"]).to_pylist()}
    assert used == {p["protein_uid"] for p in prots}


def test_r04_proteins_match_independent_translation(root):
    store = SeqStore(root, ["S0"])
    prot = {p["protein_uid"]: p["seq"] for p in _read(root, GLOBS["protein_seq"]).to_pylist()}
    codes = {r["feature_uid"]: r["genetic_code"] for r in _read(root, GLOBS["feature_protein"]).to_pylist()}
    checked = 0
    for r in _read(root, GLOBS["protein_occurrence"]).to_pylist():
        if r["dist_to_contig_end_bp"] == 0:
            continue  # partial gene at a contig edge
        nt = store.get(r["contig_uid"], r["start"], r["end"], r["strand"], circular=r["contig_uid"] in CIRCULAR)
        assert len(nt) % 3 == 0
        aa = _translate(nt, codes[r["feature_uid"]])
        assert aa.endswith("*"), "complete gene must end in a stop codon"
        expected = "M" + aa[1:-1]  # Prodigal reports alternative starts as M
        assert prot[r["protein_uid"]] == expected
        checked += 1
    assert checked > 700


def test_gene_counts_are_plausible(root):
    # Sanity bands against NCBI annotation (~73, ~278, ~480 protein genes), not truth.
    occ = _read(root, GLOBS["protein_occurrence"]).to_pylist()
    by_contig = {}
    for r in occ:
        by_contig[r["contig_uid"]] = by_contig.get(r["contig_uid"], 0) + 1
    bands = {"NC_001416.1": (50, 85), "NC_000866.4": (230, 310), "NC_000908.2": (450, 570)}
    for o in OBJS:
        n = by_contig[ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)]
        lo, hi = bands[o["accession"]]
        assert lo <= n <= hi, (o["accession"], n)


def test_mycoplasma_needs_code_4(root):
    store = SeqStore(root, ["S0"])
    mg = ids.contig_uid(fasta.read(FIX / "NC_000908.2.fna.gz")[0].seq)
    c4 = call_genome(store, GenomeInput((mg,), 4, circular=(mg,)))
    c11 = call_genome(store, GenomeInput((mg,), 11, circular=(mg,)))
    mean = lambda cs: sum(len(c["aa"]) for c in cs) / len(cs)  # noqa: E731
    assert mean(c4) > 1.8 * mean(c11)
    assert len(c11) > 1.5 * len(c4)


def test_ly08_rerun_is_byte_identical(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    _build(a)
    _build(b)
    fa = {str(p.relative_to(a)): p.read_bytes() for p in sorted(a.rglob("*.parquet"))}
    fb = {str(p.relative_to(b)): p.read_bytes() for p in sorted(b.rglob("*.parquet"))}
    assert fa == fb


@pytest.mark.parametrize("acc", ["NC_000866.4", "NC_000908.2"])
def test_mr03_reverse_complement_keeps_proteins(tmp_path, acc):
    # Prodigal is not exactly strand-symmetric (tie-breaking in its dynamic programming),
    # so the relation is checked with a measured tolerance: >= 98% of proteins shared.
    o = next(o for o in OBJS if o["accession"] == acc)
    seq = fasta.read(FIX / o["file"])[0].seq
    src = tmp_path / "src"
    src.mkdir()
    rc_file = src / "rc.fna"
    rc_file.write_bytes(b"\n".join(fasta.format_record("rc", revcomp(seq), 70)) + b"\n")
    root = tmp_path / "s"
    ingest([SourceObject(FIX / o["file"], "a", "1"), SourceObject(rc_file, "b", "1")], root, "S0")
    store = SeqStore(root, ["S0"])
    fwd = call_genome(store, _genome(o, seq))
    rev = call_genome(store, _genome(o, revcomp(seq)))
    pf, pr = {c["aa"] for c in fwd}, {c["aa"] for c in rev}
    assert len(pf & pr) / max(len(pf), len(pr)) >= 0.98
    n = len(seq)

    def canon(s, e, strand):  # circular coordinates mod n
        return (s % n, e - s, strand) if o["topology"] == "circular" else (s, e - s, strand)

    mirrored = {canon(n - c["end"], n - c["start"], -c["strand"]) for c in rev}
    fwd_coords = {canon(c["start"], c["end"], c["strand"]) for c in fwd}
    assert len(fwd_coords & mirrored) / len(fwd_coords) >= 0.98


def test_viruses_route_to_gv_and_mycoplasma_to_single(root):
    store = SeqStore(root, ["S0"])
    modes = {}
    for o in OBJS:
        calls = call_genome(store, _genome(o))
        modes[o["accession"]] = {c["mode"] for c in calls}
        codes = {c["genetic_code"] for c in calls}
        if o["kind"] == "prokaryote":
            assert codes == {o["genetic_code"]}
    assert modes == {"NC_001416.1": {"gv-meta"}, "NC_000866.4": {"gv-meta"}, "NC_000908.2": {"single"}}


def test_mr04_rotating_a_circular_genome_keeps_genes(tmp_path):
    # Rotate M. genitalium so that the new origin falls inside a gene; that gene must be called
    # once, wrapping the origin, and the gene set must match the unrotated genome.
    o = next(o for o in OBJS if o["accession"] == "NC_000908.2")
    seq = fasta.read(FIX / o["file"])[0].seq
    n = len(seq)
    src = tmp_path / "src"
    src.mkdir()
    root = tmp_path / "s"
    ingest([SourceObject(FIX / o["file"], "a", "1")], root, "S0")
    base = call_genome(SeqStore(root, ["S0"]), _genome(o, seq))
    target = base[len(base) // 2]
    r = (target["start"] + target["end"]) // 2  # inside a mid-genome gene
    rotated = seq[r:] + seq[:r]
    rot_file = src / "rot.fna"
    rot_file.write_bytes(b"\n".join(fasta.format_record("rot", rotated, 70)) + b"\n")
    ingest([SourceObject(rot_file, "b", "1")], root, "S1")
    rot = call_genome(SeqStore(root, ["S0", "S1"]), _genome(o, rotated))
    wrap = [c for c in rot if c["wraps_origin"]]
    assert any(c["aa"] == target["aa"] for c in wrap), "the split gene must be called whole across the origin"
    key = lambda cs, shift: {((c["start"] + shift) % n, c["end"] - c["start"], c["strand"], c["aa"]) for c in cs}  # noqa: E731
    a, b = key(base, 0), key(rot, r)
    assert len(a & b) / len(a) >= 0.99
    assert len([c for c in rot if c["aa"] == target["aa"]]) == 1


def test_wrapping_slice_round_trips():
    import random

    root_seq = fasta.read(FIX / "NC_000908.2.fna.gz")[0].seq
    n = len(root_seq)
    rng = random.Random(7)
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        ingest([SourceObject(FIX / "NC_000908.2.fna.gz", "a", "1")], Path(d), "S0", shard_bases=200_000)
        store = SeqStore(Path(d), ["S0"])
        cu = ids.contig_uid(root_seq)
        doubled = root_seq + root_seq
        for _ in range(100):
            s = rng.randrange(0, n)
            e = rng.randrange(n + 1, s + n + 1) if s + n > n else n
            assert store.get(cu, s, e, 1, circular=True) == doubled[s:e]
            assert store.get(cu, s, e, -1, circular=True) == revcomp(doubled[s:e])
