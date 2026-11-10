"""Ring R2: tables cut into several parts, and every reader takes all of them (card B-19).

The micro-S0 fixture is built twice. Store A uses the defaults, so every table is one
`part-0.parquet`. Store B passes a small `max_file_bytes` to every builder, so each table of the card
is cut into several `part-<k>.parquet` files. The oracle for a split table is store A's single file:
the parts of B, read in part order, must equal it row for row. The downstream builders (IGRs, tapes,
tokens, postings, the dimension table) and the checks then run on B, and their outputs must equal
A's, which proves that each reader took every part.

The ledger checks are compared with an independent pure-Python implementation of L-01..L-03 (the
row-by-row code they replaced) on the fixture and on injected faults.

Covers LY-02, LY-03, LY-04 (parts), C-03, L-01, L-02, L-03, D-01 (content does not depend on layout).
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

import pyarrow.parquet as pq
import pyhmmer
import pytest
import yaml

from biobox import catalog, contracts, dim, fasta, ids, ledger, provenance
from biobox.l0seq import SeqStore, SourceObject, ingest
from biobox.l0seq.inputs import GenomeInput, iter_genome_input, write_genome_input
from biobox.l1 import domains, igr
from biobox.l1.genecall import GENE_CALL_VERSION as V
from biobox.l1.genecall import run as call_genes
from biobox.l3 import postings, tapes, tokens
from biobox.registry.schemas import conformance_problems, table

pytestmark = pytest.mark.r2

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "micro_s0"
CONTRACT = yaml.safe_load((FIX / "sources.yaml").read_text())
MMSEQS = shutil.which("mmseqs")
REL = "rS0.parts"
TV = tapes.tape_version([V])
IGRV = igr.igr_version([V])
SRC = CONTRACT["source_name"]

# max_file_bytes per builder for store B: small enough that each table it writes has several parts
SMALL = {"ingest": 1, "genome_input": 1, "genecall": 8_192, "domains": 60, "cluster": 8_192, "igr": 4_096,
         "tapes": 1, "dim": 1}

# table -> its directory for the fixture's stage, release and versions
SPLIT = {
    "contig_occurrence": f"l0seq/contig_occurrence/stage=S0/source={SRC}",
    "contig_index": "l0seq/contig_index/stage=S0",
    "seq_locator": "l0seq/seq_locator/stage=S0",
    "genome_input": "l0seq/genome_input/stage=S0",
    "feature": f"l1/feature/producer={V}/stage=S0",
    "feature_protein": f"l1/feature_protein/gene_call_version={V}/stage=S0",
    "protein_occurrence": f"l1/protein_occurrence/gene_call_version={V}/stage=S0",
    "slot_protein": f"l1/slot_protein/gene_call_version={V}/stage=S0",
    "igr": f"l0seq/igr/igr_version={IGRV}/stage=S0",
    "l1_domains": "l1/domains/library=test-1/stage=S0",
    "cluster_member": f"l2/cluster_member/release={REL}",
    "cluster_member_aln": f"l2/cluster_member_aln/release={REL}",
    "contig_dim": f"dim/contig_dim/release={REL}",
    "contig_tapes": f"l3/tapes/tape_version={TV}/stage=S0",
}
NEEDS_MMSEQS = {"cluster_member", "cluster_member_aln"}


def _genomes() -> list[GenomeInput]:
    out = []
    for o in CONTRACT["objects"]:
        cu = ids.contig_uid(fasta.read(FIX / o["file"])[0].seq)
        out.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"],
                               circular=(cu,) if o["topology"] == "circular" else ()))
    return out


def _part_order(p: Path) -> tuple[str, int]:
    """Directory, then k of `part-<k>.parquet` as a number (written here, not taken from biobox.io)."""
    m = re.fullmatch(r"part-(\d+)\.parquet", p.name)
    assert m, p
    return p.parent.as_posix(), int(m[1])


def _parts(base: Path, pattern: str = "part-*.parquet") -> list[Path]:
    return sorted(base.glob(pattern), key=_part_order)


def _rows(paths: list[Path]) -> list[dict]:
    return [r for p in paths for r in pq.read_table(p).to_pylist()]


def _library(root: Path, path: Path) -> None:
    """Six HMMs built from the longest lambda proteins (as test_micro_s0_contracts does)."""
    lam = ids.contig_uid(fasta.read(FIX / "NC_001416.1.fna.gz")[0].seq)
    prots = domains.load_proteins(root, V, "S0")
    occ = _rows(_parts(root / SPLIT["protein_occurrence"]))
    uids = sorted({r["protein_uid"] for r in occ if r["contig_uid"] == lam}, key=lambda u: (-len(prots[u]), u))[:6]
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder, bg = pyhmmer.plan7.Builder(alphabet, seed=42), pyhmmer.plan7.Background(alphabet)
    with open(path, "wb") as fh:
        for i, u in enumerate(uids):
            seq = pyhmmer.easel.TextSequence(name=b"s", sequence=prots[u].replace("*", "")).digitize(alphabet)
            hmm, _, _ = builder.build(seq, bg)
            hmm.name, hmm.accession = f"lam{i}".encode(), f"TEST{i:04d}".encode()
            hmm.cutoffs.gathering = (25.0, 25.0)
            hmm.write(fh)


def build(root: Path, lib: Path | None, small: dict[str, int] | None) -> dict:
    """The whole fixture build; `small` gives each builder's max_file_bytes (None: the defaults)."""
    kw = {k: {"max_file_bytes": v} for k, v in (small or {}).items()}
    out: dict = {}
    ingest([SourceObject(FIX / o["file"], SRC, CONTRACT["source_release"], CONTRACT["license"], True)
            for o in CONTRACT["objects"]], root, "S0", shard_bases=10**6, **kw.get("ingest", {}))
    write_genome_input(root, "S0", _genomes(), **kw.get("genome_input", {}))
    out["genecall"] = call_genes(root, "S0", **kw.get("genecall", {}))
    if lib is None:
        lib = root.parent / f"{root.name}.hmm"
        _library(root, lib)
    out["domains"] = domains.run(root, "S0", V, "test", "1", lib, "gathering", **kw.get("domains", {}))
    if MMSEQS:
        from biobox.l2 import cluster

        work = root.parent / f"{root.name}.mmseqs"
        work.mkdir(parents=True, exist_ok=True)
        out["cluster"] = cluster.run(root, REL, V, tmp=work, **kw.get("cluster", {}))
        shutil.rmtree(work, ignore_errors=True)
    out["igr"] = igr.run(root, "S0", **kw.get("igr", {}))
    out["tapes"] = tapes.run(root, "S0", **kw.get("tapes", {}))
    out["tokens"] = tokens.run(root, REL)
    out["postings"] = postings.run(root, REL)
    out["dim"] = dim.run(root, REL, **kw.get("dim", {}))
    out["lib"] = lib
    return out


@pytest.fixture(scope="module")
def stores(tmp_path_factory):
    d = tmp_path_factory.mktemp("parts")
    a = build(d / "a", None, None)
    b = build(d / "b", a["lib"], SMALL)
    man_a, man_b = catalog.build_manifest(d / "a", REL), catalog.build_manifest(d / "b", REL)
    cons = [catalog.open_catalog(catalog.generate_catalog(m, d / f"{n}.duckdb")) for n, m in (("a", man_a), ("b", man_b))]
    yield {"a": d / "a", "b": d / "b", "out_a": a, "out_b": b, "man_a": man_a, "man_b": man_b,
           "con_a": cons[0], "con_b": cons[1]}
    for c in cons:
        c.close()


# --- every table of the card splits, and its parts concatenate to the single file -----------------


@pytest.mark.covers("LY-02", "LY-03", "LY-04")
@pytest.mark.parametrize("name", sorted(SPLIT))
def test_each_table_splits_into_sorted_conforming_parts_equal_to_the_single_file(stores, name):
    if name in NEEDS_MMSEQS and not MMSEQS:
        pytest.skip("MMseqs2 not on PATH (runs in the image)")
    single = _parts(stores["a"] / SPLIT[name])
    parts = _parts(stores["b"] / SPLIT[name])
    assert [p.name for p in single] == ["part-0.parquet"]  # the default build: one file
    assert len(parts) > 1, f"{name}: {len(parts)} part"
    assert [p.name for p in parts] == [f"part-{k}.parquet" for k in range(len(parts))]
    spec = table(name)
    for p in parts:
        t = pq.read_table(p)
        assert t.num_rows > 0
        assert conformance_problems(t, spec) == [], p
    files = [{"table": name, "path": p.relative_to(stores["b"]).as_posix()} for p in parts]
    m = {"root": str(stores["b"]), "files": files}
    assert contracts.check_sort_order(m) == [] and contracts.check_statistics(m) == []
    assert _rows(parts) == _rows(single)  # in part order, row for row


@pytest.mark.covers("LY-04")
def test_protein_seq_buckets_split_too(stores):
    """`protein_seq` already cut parts (genecall._write_protein_seq); the small max_file_bytes reaches it."""
    base = f"l0seq/protein_seq/gene_call_version={V}/stage=S0"
    a = _parts(stores["a"] / base, "protein_prefix=*/part-*.parquet")
    b = _parts(stores["b"] / base, "protein_prefix=*/part-*.parquet")
    assert len(a) == 1 and len(b) > 1
    assert _rows(b) == _rows(a)


# --- downstream builders and readers on the multi-part store ---------------------------------------


@pytest.mark.covers("D-01")
def test_downstream_builders_give_the_same_results(stores):
    a, b = stores["out_a"], stores["out_b"]
    for step in ("genecall", "domains", "igr", "tokens", "postings", "dim") + (("cluster",) if MMSEQS else ()):
        assert a[step] == b[step], step
    assert {k: v for k, v in a["tapes"].items() if k != "files"} == {k: v for k, v in b["tapes"].items() if k != "files"}
    assert a["tapes"]["files"] == 1 and b["tapes"]["files"] == 3


@pytest.mark.covers("D-01")
def test_every_table_has_the_same_rows(stores):
    """Layout-independent content of every table and artifact (`contracts.content_hashes`), then the
    release-level tables row for row in file order."""
    ha, hb = contracts.content_hashes(stores["man_a"]), contracts.content_hashes(stores["man_b"])
    assert ha == hb
    for name in ("token_dict", "tape_tokens", "token_postings", "contig_dim", "tape_version_def") + (
            ("cluster_member", "cluster_member_aln") if MMSEQS else ()):
        fa = [stores["a"] / f["path"] for f in stores["man_a"]["files"] if f["table"] == name]
        fb = [stores["b"] / f["path"] for f in stores["man_b"]["files"] if f["table"] == name]
        assert fa and fb, name
        assert _rows(sorted(fa, key=_part_order)) == _rows(sorted(fb, key=_part_order)), name


@pytest.mark.covers("C-03", "LY-02", "LY-03")
def test_checks_pass_on_the_multi_part_store(stores):
    man, con = stores["man_b"], stores["con_b"]
    assert len([f for f in man["files"] if f["table"] == "seq_locator"]) == 3
    assert contracts.check_foreign_keys(con) == []
    assert contracts.check_sort_order(man) == []
    assert contracts.check_statistics(man) == []
    assert contracts.check_shards(man) == []
    assert contracts.check_list_lengths(con) == []
    assert contracts.check_coordinates(con) == []
    assert contracts.check_lineage(con) == []
    assert provenance.check_provenance(man) == []
    assert catalog.verify_manifest(man) == []


def test_catalog_views_hold_the_same_rows(stores):
    ca, cb = stores["con_a"], stores["con_b"]
    names = sorted({f["table"] for f in stores["man_a"]["files"]})
    assert len(names) > 15
    for name in names:
        ra = sorted(map(repr, ca.sql(f"SELECT * FROM {name}").to_arrow_table().to_pylist()))
        rb = sorted(map(repr, cb.sql(f"SELECT * FROM {name}").to_arrow_table().to_pylist()))
        assert ra and ra == rb, name
    # the sequence macros resolve keys through SeqStore over every part of both indexes
    sql = ("SELECT contig_uid, seq_get(contig_uid, 100, 400, -1) FROM contig_index WHERE is_default_variant "
           "ORDER BY contig_uid")
    assert ca.execute(sql).fetchall() == cb.execute(sql).fetchall()


@pytest.mark.covers("L-01", "L-02", "L-03")
def test_ledger_checks_on_the_multi_part_store(stores):
    expected = {ids.sha256_file(FIX / o["file"]) for o in CONTRACT["objects"]}
    for root in (stores["a"], stores["b"]):
        steps = sorted(p.name.removeprefix("step=") for p in (root / "ledger/stage=S0").iterdir())
        assert {"ingest", "genecall", "domains-test-1", f"tape-{TV}"} <= set(steps)
        assert ledger.check_balance(root, "S0", steps) == []
        assert ledger.check_inputs(root, "S0", expected) == []
        assert ledger.check_chain(root, "S0") == []
        for step in steps:
            assert ledger.read(root, "S0", step) == ledger.read(stores["a"], "S0", step), step
    tape_ledger = ledger.read(stores["b"], "S0", f"tape-{TV}")
    assert [r["object_key"] for r in tape_ledger] == sorted(g.contig_uids[0] for g in _genomes())


def test_seqstore_and_genome_input_read_every_part(stores):
    sa, sb = SeqStore(stores["a"]), SeqStore(stores["b"])
    try:
        assert len(sb.locator_files) == len(sb.index_files) == 3 and len(sa.locator_files) == 1
        for g in _genomes():
            (cu,) = g.contig_uids
            assert sb.get(cu) == sa.get(cu)
            assert sb.get(cu, 1_000, 2_500, -1) == sa.get(cu, 1_000, 2_500, -1)
        for o in CONTRACT["objects"]:
            sha = ids.sha256_file(FIX / o["file"])
            assert sb.reconstruct_object(sha) == sa.reconstruct_object(sha)
    finally:
        sa.close()
        sb.close()
    assert list(iter_genome_input(stores["b"], "S0")) == list(iter_genome_input(stores["a"], "S0"))
    assert len(_parts(stores["b"] / SPLIT["genome_input"])) == 3


# --- the ledger checks against an independent implementation ---------------------------------------


def _oracle_per_object(rows):
    out = {}
    for r in rows:
        o = out.setdefault(r["object_key"], {"n_in": set(), "n_kept": set(), "n_merged": set(), "rejected": 0})
        for k in ("n_in", "n_kept", "n_merged"):
            o[k].add(r[k])
        o["rejected"] += r["n_rejected"]
    return out


def _oracle_balance(root, stage, steps):
    problems = []
    for step in steps:
        for uid, o in _oracle_per_object(ledger.read(root, stage, step)).items():
            if len(o["n_in"]) != 1 or len(o["n_kept"]) != 1 or len(o["n_merged"]) != 1:
                problems.append(f"{step} {uid!r}: n_in/n_kept/n_merged differ between reason rows")
                continue
            (n_in,), (kept,), (merged,) = o["n_in"], o["n_kept"], o["n_merged"]
            if n_in != kept + merged + o["rejected"]:
                problems.append(f"{step} {uid!r}: {n_in} in != {kept} kept + {merged} merged + {o['rejected']} rejected")
    return problems


def _oracle_inputs(root, stage, expected):
    got = {r["object_key"] for r in ledger.read(root, stage, "ingest")}
    problems = []
    if expected - got:
        problems.append(f"ingest ledger misses {len(expected - got)} L0 objects")
    if got - expected:
        problems.append(f"ingest ledger has {len(got - expected)} objects not in the input")
    return problems


def _oracle_chain(root, stage):
    problems = []
    for up, up_obj, down, down_obj in ledger.CHAIN:
        urows = ledger.read(root, stage, up)
        if not urows:
            continue
        kept = sum(next(iter(o["n_kept"])) for u, o in _oracle_per_object(urows).items() if up_obj is None or u == up_obj)
        downs = ([p.name.removeprefix("step=") for p in sorted((Path(root) / f"ledger/stage={stage}").glob(f"step={down}"))]
                 if down.endswith("*") else [down])
        for u, o in _oracle_per_object([r for d in downs for r in ledger.read(root, stage, d)]).items():
            if down_obj is not None and u != down_obj:
                continue
            if next(iter(o["n_in"])) != kept:
                problems.append(f"{up}->{down} {u!r}: n_in {next(iter(o['n_in']))} != upstream n_kept {kept}")
    return problems


def _rewrite(root: Path, step: str, rows: list[dict]) -> None:
    ledger.path(root, "S0", step).unlink()
    ledger.write(root, "S0", step, sorted(rows, key=lambda r: (r["object_key"], r["reason_code"])))


@pytest.mark.covers("L-01", "L-02", "L-03")
@pytest.mark.parametrize("fault", ["none", "unbalanced contig", "reason rows differ", "object missing",
                                   "object extra", "domains chain", "cluster chain", "second library"])
def test_ledger_checks_agree_with_the_row_by_row_oracle(stores, tmp_path, fault):
    root = tmp_path / "s"
    shutil.copytree(stores["b"] / "ledger", root / "ledger")
    expected = {ids.sha256_file(FIX / o["file"]) for o in CONTRACT["objects"]}
    gc, dom, ing = (ledger.read(root, "S0", s) for s in ("genecall", "domains-test-1", "ingest"))
    if fault == "unbalanced contig":
        gc[1]["n_kept"] += 1
        _rewrite(root, "genecall", gc)
    elif fault == "reason rows differ":
        dom[0]["n_kept"] += 2  # two reason rows of one object disagree; not on the chain's counts
        _rewrite(root, "domains-test-1", dom)
    elif fault == "object missing":
        _rewrite(root, "ingest", ing[1:])
    elif fault == "object extra":
        _rewrite(root, "ingest", ing + [{**ing[0], "object_key": b"\xff" * 32}])
    elif fault == "domains chain":
        _rewrite(root, "domains-test-1", [{**r, "n_in": r["n_in"] + 5} for r in dom])
    elif fault == "cluster chain":
        kept = sum(r["n_kept"] for r in gc)
        shutil.rmtree(root / "ledger/stage=S0/step=cluster", ignore_errors=True)  # present when MMseqs2 ran
        ledger.write(root, "S0", "cluster", [ledger.row(b"nr50", "level", 7, 3, 4), ledger.row(b"nr90", "level", kept, 9, kept - 9)])
    elif fault == "second library":
        ledger.write(root, "S0", "domains-a-1", [ledger.row(b"a-1", "library", 3, 1, n_rejected=2)])
    steps = sorted(p.name.removeprefix("step=") for p in (root / "ledger/stage=S0").iterdir())
    got = (ledger.check_balance(root, "S0", steps), ledger.check_inputs(root, "S0", expected), ledger.check_chain(root, "S0"))
    want = (_oracle_balance(root, "S0", steps), _oracle_inputs(root, "S0", expected), _oracle_chain(root, "S0"))
    assert got == want
    assert any(got) == (fault != "none")
