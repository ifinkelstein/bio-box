"""Card W-23: commands a user needs and does not have.

(1) `bio-box releases`, (2) `bio-box seq`, (3) `bio-box sql`, (4) `stage` and `kind` beside each gene of
`hits.tsv`, (5) 1-based inclusive GenBank coordinates (`ACC:START..END`) in `find` and `search hmm --positive`,
(6) one stderr line per output file with its count.

Oracles are loops over the Parquet files and the fixture's FASTA files, written here, independently of the code
under test. The releases are the micro-S0 search world (`test_micro_s0_search.world`), a copy of it with a
`contig_dim`, a manifest and a catalog (`rel`), and the two-stage store of card B-28 (`test_release_readers.store`).
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

import duckdb
import pyarrow.parquet as pq
import pytest

sys.path.insert(0, str(Path(__file__).parent))
import test_micro_s0_l3 as micro  # noqa: E402
import test_release_readers as rr  # noqa: E402
from test_micro_s0_search import (REL, V, _copy, _packs, _packs_of, _recipe, _run, _runs, _runs_of,  # noqa: E402
                                  _seed, _w03_profile, _w06_anchors, _w07_genes, _w07_hmm, _w07_labels, _w07_spec,
                                  world)  # noqa: F401  (the module-scoped search world)
from test_release_readers import store  # noqa: E402,F401  (the two-stage store: R1 over S0, R2 over S0 and S1)

from biobox import catalog, cli, dim, fasta, ids  # noqa: E402
from biobox.l0seq.inputs import GenomeInput, write_genome_input  # noqa: E402
from biobox.l3 import tapes  # noqa: E402
from biobox.search import export, runner  # noqa: E402

pytestmark = [pytest.mark.r2, pytest.mark.covers("W-23")]

NODEF = "rX.w23"  # a release with a release_def and no manifest
TV = tapes.tape_version([V])
WROTE = re.compile(r"^wrote (?P<path>.+): (?P<n>\d+) (?P<unit>rows|sequences|records|tables, (?P<rows>\d+) rows)$")


def _s0_genome_input(root: Path) -> None:
    """Stage S0's `genome_input`, as `bio-box stage inputs` writes it (micro.build hands the genomes to the gene
    caller without persisting them): one genome per fixture record, with the kind and topology of sources.yaml."""
    gs = []
    for o in micro.CONTRACT["objects"]:
        cu = ids.contig_uid(fasta.read(micro.FIX / o["file"])[0].seq)
        gs.append(GenomeInput((cu,), o["genetic_code"], kind=o["kind"],
                              circular=(cu,) if o["topology"] == "circular" else ()))
    write_genome_input(root, "S0", gs)


@pytest.fixture(scope="module")
def rel(world, tmp_path_factory):
    """A copy of the world with S0's `genome_input`, REL's `contig_dim` (the dim step), manifest and catalog, and a
    second release NODEF over the same stage with only a `release_def`."""
    d = tmp_path_factory.mktemp("w23")
    root = _copy(world, d)
    _s0_genome_input(root)
    (d / "dimtmp").mkdir()
    dim.run(root, REL, ["S0"], tmp=d / "dimtmp")
    man = catalog.build_manifest(root, REL, ["S0"])
    catalog.write_manifest(man)
    catalog.generate_catalog(man, root=root)
    catalog.define_release(root, NODEF, ["S0"])
    return root


def _tree(root: Path) -> dict[str, str]:
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else "dir"
            for p in sorted(root.rglob("*"))}


def _rows(root: Path, base: str, columns=None) -> list[dict]:
    """Rows of every Parquet file under `root/base`, each with the `stage` of its path when it has one."""
    out = []
    for f in sorted((root / base).rglob("*.parquet")):
        stage = next((p.split("=", 1)[1] for p in f.relative_to(root).parts if p.startswith("stage=")), None)
        out += [{**r, "_stage": stage} for r in pq.read_table(f, columns=columns, partitioning=None).to_pylist()]
    return out


def _proteins(root: Path) -> dict[bytes, dict]:
    return {r["protein_uid"]: r for r in _rows(root, "l0seq/protein_seq")}


def _contig_seqs(root: Path) -> dict[str, tuple[bytes, bytes]]:
    """{accession: (contig_uid, sequence)} from the fixture's own FASTA files, matched to the labels by header."""
    seqs = {r.header: r.seq for f in sorted(micro.FIX.glob("*.fna.gz")) for r in fasta.read(f)}
    return {acc: (cu, seqs[hdr]) for cu, (acc, hdr) in _w07_labels(root).items()}


def _revcomp(s: bytes) -> bytes:
    return s[::-1].translate(bytes.maketrans(b"ACGTacgtNn", b"TGCAtgcaNn"))


def _wrote(err: str) -> dict[str, re.Match]:
    return {m["path"]: m for m in (WROTE.match(x) for x in err.splitlines()) if m}


# --- (1) bio-box releases -----------------------------------------------------------------------------------


def _release_oracle(root: Path) -> dict:
    dims = _rows(root, f"dim/contig_dim/release={REL}")
    return {"genomes": len({r["genome_uid"] for r in dims}), "contigs": len(dims),
            "proteins": len({r["protein_uid"] for r in _rows(root, "l0seq/protein_seq", ["protein_uid"])})}


def _releases(capsys, *args) -> dict[str, dict]:
    assert cli.main(["releases", *args, "--json"]) == 0
    return {r["release"]: r for r in map(json.loads, capsys.readouterr().out.splitlines())}


def test_w23_releases_lists_each_release_with_its_stages_counts_libraries_and_lookups(rel, capsys, monkeypatch):
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)
    want = _release_oracle(rel)
    assert want["genomes"] == want["contigs"] == 3 and want["proteins"] > 100
    recs = _releases(capsys, "--root", str(rel))
    assert sorted(recs) == sorted([REL, NODEF])
    r = recs[REL]
    assert r["manifest"] is True and r["stages"] == ["S0"]
    assert {k: r[k] for k in want} == want
    assert r["libraries"] == ["test-1"] and r["ncrna"] == []
    assert r["lookups"] == {"member_cluster": True, "contig_label": True, "domain_family_tokens": True}
    assert "pack" not in r  # neither a profile nor a pack directory: no pack is looked for
    n = recs[NODEF]  # listed, not skipped: no manifest, its stages from its release_def
    assert n["manifest"] is False and n["stages"] == ["S0"]
    assert (n["contigs"], n["proteins"], n["lookups"]["member_cluster"]) == (None, None, False)
    # the table for a person: a column line, then one line per release
    assert cli.main(["releases", "--root", str(rel)]) == 0
    lines = [x.split("\t") for x in capsys.readouterr().out.splitlines()]
    assert lines[0][:3] == ["release", "manifest", "stages"] and len(lines) == 3
    row = dict(zip(lines[0], next(x for x in lines[1:] if x[0] == REL)))
    assert (row["manifest"], row["stages"], row["proteins"], row["libraries"]) == ("yes", "S0", str(want["proteins"]),
                                                                                    "test-1")
    assert dict(zip(lines[0], next(x for x in lines[1:] if x[0] == NODEF)))["manifest"] == "no"


def test_w23_releases_counts_from_footers_and_says_which_lookups_an_older_release_lacks(rel, tmp_path, capsys):
    """A release built before cards B-27 and B-39 (no `member_cluster`, no `domain_family_tokens` column): the
    lookups say so, and the protein count comes from `cluster_member` instead."""
    root = tmp_path / "old"
    shutil.copytree(rel, root)
    shutil.rmtree(root / "l2/member_cluster")
    for f in sorted((root / "l3/tape_tokens").rglob("*.parquet")):
        t = pq.read_table(f, partitioning=None)
        pq.write_table(t.drop_columns([catalog.TAPE_DOMAIN_FAMILIES]), f)
    catalog.manifest_path(root, REL).unlink()
    catalog.write_manifest(catalog.build_manifest(root, REL, ["S0"]))
    r = _releases(capsys, "--root", str(root))[REL]
    assert r["lookups"] == {"member_cluster": False, "contig_label": True, "domain_family_tokens": False}
    assert r["proteins"] == _release_oracle(root)["proteins"]


def test_w23_releases_names_the_pack_the_profiles_pack_directory_holds(rel, world, tmp_path, capsys, monkeypatch):
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)
    prof = _w03_profile(tmp_path)
    assert _releases(capsys, "--root", str(rel), "--profile", str(prof))[REL]["pack"] is None
    shutil.copytree(_packs_of(world[0]), tmp_path / "ptmp" / "packs")
    recs = _releases(capsys, "--root", str(rel), "--profile", str(prof))
    assert recs[REL]["pack"] == world[2]["pack_id"] and recs[NODEF]["pack"] is None
    assert _releases(capsys, "--root", str(rel), "--pack-dir", str(_packs_of(world[0])))[REL]["pack"] \
        == world[2]["pack_id"]


# --- (2) bio-box seq ----------------------------------------------------------------------------------------


def _protein_fasta(prots: dict, uids) -> str:
    return "".join(f">{u.hex()} release={REL} length={prots[u]['length']}\n{prots[u]['seq']}\n" for u in uids)


def test_w23_seq_prints_proteins_by_uid_as_fasta(rel, tmp_path, capsys):
    prots = _proteins(rel)
    uids = sorted(prots)[::len(prots) // 3][:3]
    base = ["seq", "--root", str(rel), "--release", REL]
    assert cli.main([*base, "--protein", *[u.hex() for u in uids]]) == 0
    assert capsys.readouterr().out == _protein_fasta(prots, uids)
    assert cli.main([*base, "--protein", uids[0].hex(), "--protein", uids[1].hex()]) == 0  # repeated, in that order
    assert capsys.readouterr().out == _protein_fasta(prots, uids[:2])
    f = tmp_path / "uids.txt"
    f.write_text("# proteins\n" + "\n".join(u.hex() for u in reversed(uids)) + "\n\n")
    out = tmp_path / "p.faa"
    assert cli.main([*base, "--proteins-file", str(f), "--out", str(out)]) == 0
    assert out.read_text() == _protein_fasta(prots, list(reversed(uids)))
    # an unknown uid: exit 2, the uid named, nothing written
    capsys.readouterr()
    assert cli.main([*base, "--protein", uids[0].hex(), "00" * 16, "--out", str(tmp_path / "q.faa")]) == 2
    err = capsys.readouterr().err
    assert "00" * 16 in err and not (tmp_path / "q.faa").exists()
    assert cli.main([*base, "--protein", "not-a-uid"]) == 2


def test_w23_seq_reads_only_the_protein_seq_row_groups_that_hold_the_keys(rel, tmp_path):
    """AGENTS §4.4: a lookup by key reads O(1) row groups. The copy's `protein_seq` is rewritten with 4 rows per row
    group, so that reading every row group and reading one per key differ."""
    root = tmp_path / "rg"
    shutil.copytree(rel, root)
    for f in sorted((root / "l0seq/protein_seq").rglob("*.parquet")):
        pq.write_table(pq.read_table(f, partitioning=None), f, row_group_size=4)
    catalog.manifest_path(root, REL).unlink()
    catalog.write_manifest(catalog.build_manifest(root, REL, ["S0"]))
    files = catalog.ReleaseFiles.of(root, REL)
    total = sum(pq.ParquetFile(p).metadata.num_row_groups for p in files.table("protein_seq", gene_call_version=V))
    prots = _proteins(root)
    uids = sorted(prots)[::len(prots) // 3][:3]
    got, scans = export.protein_seqs(root, uids, release=REL, gene_call_version=V, manifest=files)
    assert got == {u: (prots[u]["length"], prots[u]["seq"]) for u in uids}
    read = sum(s["row_groups_scanned"] for s in scans)
    assert total >= 20 and 1 <= read <= 2 * len(uids), (read, total, scans)


def test_w23_seq_slices_a_contig_by_accession_with_coordinates_and_strand(rel, capsys):
    seqs = _contig_seqs(rel)
    acc = next(a for a in seqs if a.startswith("NC_001416"))  # lambda: linear
    cu, s = seqs[acc]
    base = ["seq", "--root", str(rel), "--release", REL]

    def got(*args) -> tuple[str, str]:
        assert cli.main([*base, *args]) == 0
        head, body = capsys.readouterr().out.splitlines()
        return head, body

    head, body = got("--contig", acc)
    assert body == s.decode() and head == (f">{cu.hex()} release={REL} accession={acc} start=0 end={len(s)} "
                                           f"strand=+1 contig_len={len(s)}")
    head, body = got("--contig", f"{acc}:100-250")
    assert body == s[100:250].decode()
    assert head == f">{cu.hex()} release={REL} accession={acc} start=100 end=250 strand=+1 contig_len={len(s)}"
    assert got("--contig", f"{acc}:101..250") == (head, body)  # the GenBank form of the same bases
    assert got("--contig", f"{cu.hex()}:100-250") == (head, body)  # a contig uid in place of the accession
    head, body = got("--contig", f"{acc}:100-250", "--strand", "-1")
    assert body == _revcomp(s[100:250]).decode() and " strand=-1 " in head
    # across the origin of a circular contig (end > contig_len), and refused on a linear one
    circ = next(a for a in seqs if a.startswith("NC_000908"))
    ccu, cs = seqs[circ]
    n = len(cs)
    _, body = got("--contig", f"{circ}:{n - 50}-{n + 50}")
    assert body == (cs[n - 50:] + cs[:50]).decode()
    for bad in (f"{acc}:0-{len(s) + 1}", "NC_999999.9", f"{acc}:250-100", f"{acc}:1.-5"):
        assert cli.main([*base, "--contig", bad]) == 2, bad
        assert capsys.readouterr().out == "", bad


def test_w23_seq_of_a_contig_slice_reads_one_slice_of_one_shard(rel):
    seqs = _contig_seqs(rel)
    acc = next(a for a in seqs if a.startswith("NC_000908"))  # 580 kb: many BGZF blocks of at most 64 KiB
    cu, s = seqs[acc]
    assert len(s) > 4 * 65536
    got, n, reads = cli.contig_slice(rel, catalog.ReleaseFiles.of(rel, REL), cu, 300_000, 300_150, 1, False)
    assert got == s[300_000:300_150] and n == len(s)
    assert reads["blocks_read"] <= 2 and reads["index_row_groups_read"] <= 2 and reads["locator_row_groups_read"] <= 2, \
        reads


# --- (3) bio-box sql ----------------------------------------------------------------------------------------


def _sql(rel_root: Path, *args) -> list[str]:
    return ["sql", "--root", str(rel_root), "--release", REL, *args]


def test_w23_sql_runs_selects_through_the_catalog_and_prints_the_row_groups_it_read(rel, tmp_path, capsys):
    before = _tree(rel)
    prots = _proteins(rel)
    u = sorted(prots)[0]
    q = f"SELECT protein_uid, length FROM protein_seq WHERE protein_uid = unhex('{u.hex()}')"
    assert cli.main(_sql(rel, "-e", q)) == 0
    got = capsys.readouterr()
    assert got.out.splitlines() == ["protein_uid\tlength", f"{u.hex()}\t{prots[u]['length']}"]
    assert re.search(r"protein_seq \d+ of \d+ row groups", got.err), got.err
    # a file of two statements: two tables, a blank line between
    f = tmp_path / "q.sql"
    f.write_text(q + ";\n-- a comment\nSELECT key, value FROM catalog_info WHERE key = 'release_id';\n")
    assert cli.main(_sql(rel, str(f))) == 0
    blocks = capsys.readouterr().out.rstrip("\n").split("\n\n")
    assert len(blocks) == 2 and blocks[1].splitlines() == ["key\tvalue", f"release_id\t{REL}"]
    # JSON, one object per row, uids as hex
    assert cli.main(_sql(rel, "-e", q, "--format", "json")) == 0
    assert [json.loads(x) for x in capsys.readouterr().out.splitlines()] == [{"protein_uid": u.hex(),
                                                                              "length": prots[u]["length"]}]
    # the catalog's functions are there (open_catalog registers the sequence readers)
    acc, (cu, s) = next(iter(sorted(_contig_seqs(rel).items())))
    assert cli.main(_sql(rel, "-e", f"SELECT seq_get(unhex('{cu.hex()}'), 10, 40, -1) AS s")) == 0
    assert capsys.readouterr().out.splitlines() == ["s", _revcomp(s[10:40]).decode()]
    assert cli.main(_sql(rel, "-e", q, "--out", str(tmp_path / "o.tsv"))) == 0
    assert (tmp_path / "o.tsv").read_text().splitlines()[1] == f"{u.hex()}\t{prots[u]['length']}"
    assert _tree(rel) == before  # nothing written under the root, no spill or copy left there


@pytest.mark.parametrize("stmt", [
    "COPY (SELECT 1 AS x) TO '{root}/w23.csv'",
    "SELECT 1; COPY (SELECT 1 AS x) TO '{root}/w23.csv'",
    "CREATE TABLE w23 AS SELECT 1",
    "CREATE TEMP TABLE w23 AS SELECT 1",
    "INSERT INTO protein_seq SELECT * FROM protein_seq",
    "ATTACH '{root}/w23.duckdb' AS w23",
    "INSTALL httpfs",
    "LOAD httpfs",
    "PRAGMA enable_profiling",
    "SET temp_directory = '{root}'",
    "EXPORT DATABASE '{root}/w23'",
    "CALL pragma_version()",
])
def test_w23_sql_refuses_a_statement_that_writes_or_changes_a_setting(rel, stmt, capsys):
    before = _tree(rel)
    assert cli.main(_sql(rel, "-e", stmt.format(root=rel))) == 2
    got = capsys.readouterr()
    assert got.out == "" and "refused" in got.err and "SELECT" in got.err, got.err
    assert _tree(rel) == before


def test_w23_sql_connection_cannot_write_into_the_root_past_the_statement_check(rel, tmp_path):
    """The statement check is not the only protection: the connection itself is read-only. Each statement here is run
    on it directly, as if the check had let it through."""
    before = _tree(rel)
    some = next(iter(sorted((rel / "l0seq/protein_seq").rglob("*.parquet"))))
    con, scratch = cli.sql_connection(rel, REL, runs=[], tmp=tmp_path, threads=1, memory_limit="1GB")
    try:
        assert con.execute("SELECT count(*) FROM protein_seq").fetchall() == [(len(_proteins(rel)),)]
        for stmt in (f"COPY (SELECT 1 AS x) TO '{rel / 'w23.csv'}'",
                     f"COPY (SELECT * FROM protein_seq) TO '{some}' (FORMAT parquet)",  # a file the views read
                     f"COPY (SELECT 1 AS x) TO '{tmp_path / 'elsewhere.csv'}'",
                     "CREATE TABLE w23 AS SELECT 1",
                     f"ATTACH '{tmp_path / 'n.duckdb'}' AS n",
                     "SET enable_external_access = true",
                     "SET allowed_directories = ['/']",
                     "INSTALL httpfs",
                     f"SELECT * FROM read_csv('{tmp_path / 'any.csv'}')"):
            with pytest.raises(duckdb.Error):
                con.execute(stmt)
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)
    assert _tree(rel) == before and not (tmp_path / "elsewhere.csv").exists()
    assert not scratch.is_relative_to(rel)


def test_w23_sql_refuses_an_output_under_the_root_and_a_release_without_a_catalog(rel, capsys):
    assert cli.main(_sql(rel, "-e", "SELECT 1 AS x", "--out", str(rel / "x.tsv"))) == 2
    assert "data root" in capsys.readouterr().err and not (rel / "x.tsv").exists()
    assert cli.main(["sql", "--root", str(rel), "--release", NODEF, "-e", "SELECT 1 AS x"]) == 2
    assert "bio-box release catalog" in capsys.readouterr().err


def test_w23_sql_attaches_the_runs_of_the_release(rel, world, tmp_path, capsys, monkeypatch):
    hit, _ = _seed(world[0])
    row = _run(world, _recipe(hit, hit, tool="jackhmmer", mx=2))
    assert row["status"] == "ok", row["error"]
    want = sorted(d.name.removeprefix("run=") for d in _runs(world).glob("run=*")
                  if (d / "run.parquet").exists() and runner.read_run(d)["ran_against"] == REL)
    q = "SELECT DISTINCT run FROM search_run ORDER BY run"
    assert cli.main(_sql(rel, "--runs", str(_runs(world)), "-e", q)) == 0
    assert capsys.readouterr().out.splitlines() == ["run", *want] and want
    one = runner.run_dir(_runs(world), row["run_id"])
    assert cli.main(_sql(rel, "--runs", str(one), "-e", "SELECT count(*) AS n FROM search_hit WHERE run = "
                                                    f"'{row['run_id'].hex()}'")) == 0
    n = len(pq.read_table(one / "hits.parquet"))
    assert capsys.readouterr().out.splitlines() == ["n", str(n)]
    # with a profile and no --runs: the runs directory the profile gives this user
    monkeypatch.delenv("BIOBOX_PROFILE", raising=False)
    monkeypatch.setenv("LOGNAME", "w23")
    prof = _w03_profile(tmp_path)
    shutil.copytree(_runs(world), tmp_path / "ptmp/search/w23/runs")
    assert cli.main(_sql(rel, "--profile", str(prof), "-e", q)) == 0
    assert capsys.readouterr().out.splitlines() == ["run", *want]
    assert not list((tmp_path / "ptmp").glob("biobox-sql-*"))  # its scratch under the profile's tmp is removed


# --- (4) stage and kind in hits.tsv ---------------------------------------------------------------------------


def test_w23_hits_tsv_names_each_genes_stage_and_genome_kind(store, tmp_path):
    """R2 is over S0 and S1; protein P sits on an S0 contig and on an S1 contig. Each gene's stage is the stage its
    `protein_occurrence` row is in, and its kind is the `genome_input` kind of its contig in that stage (a copy of the
    store, given S0's `genome_input`; S1's is the store's own, all virus)."""
    root = tmp_path / "a"
    shutil.copytree(store["root"], root)
    _s0_genome_input(root)
    occ = _rows(root, f"l1/protein_occurrence/gene_call_version={V}")
    kind = {(r["contig_uid"], r["_stage"]): r["kind"] for r in _rows(root, "l0seq/genome_input")}
    prok = next(r["protein_uid"] for r in occ if kind[(r["contig_uid"], r["_stage"])] == "prokaryote")
    picked = {store["p"], prok}
    hits = [{"target_kind": "protein", "target_uid": u, "seq_bitscore": 40.0 + k, "seq_evalue_zref": 1e-9}
            for k, u in enumerate(sorted(picked))]
    lines = list(export.hits_tsv(root, hits, release=rr.R2, tape_version=TV, gene_call_version=V))
    cols = lines[0].split("\t")
    # card C-21c appended n_genomes and genome_uids after them
    assert cols[-4:-2] == ["stage", "kind"] and cols == [c for c in export.HITS_TSV_COLUMNS if c != "tier"]
    rows = [dict(zip(cols, x.split("\t"))) for x in lines[1:]]
    want = {(r["feature_uid"].hex(), r["_stage"], kind[(r["contig_uid"], r["_stage"])]) for r in occ
            if r["protein_uid"] in picked}
    assert {(r["feature_uid"], r["stage"], r["kind"]) for r in rows} == want
    assert {s for _, s, _ in want} == {"S0", "S1"} and {k for _, _, k in want} == {"virus", "prokaryote"}
    assert "stage" in export.README_TEXT and "kind" in export.README_TEXT


# --- (5) GenBank coordinates --------------------------------------------------------------------------------


def _lambda_gene(root: Path) -> dict:
    return min((g for g in _w07_genes(root) if "phage lambda" in g["header"]), key=lambda g: g["start"])


def _find(root: Path, capsys, spec: str) -> tuple[int, list[dict], str]:
    code = cli.main(["find", "--root", str(root), "--release", REL, "--accession", spec])
    got = capsys.readouterr()
    lines = [x.split("\t") for x in got.out.splitlines()]
    return code, [dict(zip(lines[0], x)) for x in lines[1:]], got.err


def test_w23_find_takes_genbank_coordinates_and_prints_both_forms(world, capsys):
    root = world[0]
    g = _lambda_gene(root)
    acc, s, e = g["accession"], g["start"], g["end"]
    code, zero, err = _find(root, capsys, f"{acc}:{s}-{e}")
    assert code == 0, err
    code, one, err = _find(root, capsys, f"{acc}:{s + 1}..{e}")
    assert code == 0, err
    assert zero == one and len(one) == 1
    assert one[0]["feature_uid"] == g["feature_uid"].hex() and (one[0]["start"], one[0]["end"]) == (str(s), str(e))
    assert one[0]["span"] == f"{acc}:{s}-{e}" and one[0]["genbank_span"] == f"{acc}:{s + 1}..{e}"
    code, rows, _ = _find(root, capsys, acc)  # every gene of the contig, each in both forms
    assert code == 0 and len(rows) > 20
    assert all(r["span"] == f"{acc}:{r['start']}-{r['end']}" and r["genbank_span"] == f"{acc}:{int(r['start']) + 1}.."
               f"{r['end']}" for r in rows)


def test_w23_search_hmm_positive_takes_genbank_coordinates(world):
    root = world[0]
    g = _lambda_gene(root)
    files = catalog.ReleaseFiles.of(root, REL)
    want = (g["protein_uid"].hex(), [])
    assert cli._positive_uid(root, files, f"{g['accession']}:{g['start'] + 1}..{g['end']}", TV, V) == want
    assert cli._positive_uid(root, files, f"{g['accession']}:{g['start']}-{g['end']}", TV, V) == want
    # the label of a GenBank-form control still names its accession (the FASTA's example=, card W-20)
    ctl = json.dumps({"positive": [{"uid": "ab" * 16, "label": "NC_001604:3171..5822"}]})
    assert cli._positive_prefs(ctl) == {bytes.fromhex("ab" * 16): ("NC_001604", None)}


BAD_SPANS = ["{acc}:{s1}.-{e}", "{acc}:{s1}...{e}", "{acc}:{s1}.{e}", "{acc}:{s1}-{e}..", "{acc}:{s1}..{e}-{e}",
             "{acc}:0..{e}", "{acc}:{e}..{s1}", "{acc}:{s1}", "{acc}:{s1}--{e}"]


@pytest.mark.parametrize("bad", BAD_SPANS)
def test_w23_a_mixed_or_ambiguous_span_is_refused_with_both_spellings(world, tmp_path, capsys, bad):
    root = world[0]
    g = _lambda_gene(root)
    spec = bad.format(acc=g["accession"], s1=g["start"] + 1, e=g["end"])
    a = re.escape(g["accession"])
    for argv in (["find", "--root", str(root), "--release", REL, "--accession", spec],
                 ["search", "hmm", "--root", str(root), "--release", REL, "--hmm", str(_w07_hmm(world, tmp_path)[0]),
                  "--bits", "25", "--positive", spec, "--out", str(tmp_path / "out"), "--cpus", "1", "--pack-dir",
                  str(_packs(world)), "--runs", str(tmp_path / "runs"), "--work", str(tmp_path / "work")]):
        assert cli.main(argv) == 2, (argv[0], spec)
        err = capsys.readouterr().err
        assert re.search(rf"{a}:\d+-\d+", err) and re.search(rf"{a}:\d+\.\.\d+", err), (spec, err)
        assert "GenBank" in err and "0-based" in err, err
    assert not (tmp_path / "out").exists() and not (tmp_path / "runs").exists()


# --- (6) one stderr line per output file, with its count ------------------------------------------------------


def _tsv_rows(text: str) -> int:
    """Data rows of a TSV of blank-line-separated tables, each a column line and its rows, with `#` lines."""
    n, fresh = 0, True
    for line in text.splitlines():
        if not line:
            fresh = True
        elif line.startswith("#"):
            continue
        elif fresh:
            fresh = False
        else:
            n += 1
    return n


def test_w23_search_hmm_prints_one_line_per_output_file_with_its_count(world, tmp_path, capsys):
    root = world[0]
    hmm, pos = _w07_hmm(world, tmp_path)
    out = tmp_path / "out"
    capsys.readouterr()
    code = cli.main(["search", "hmm", "--root", str(root), "--release", REL, "--hmm", str(hmm), "--bits", "25",
                     "--positive", _w07_spec(root, pos), "--out", str(out), "--cpus", "1", "--pack-dir",
                     str(_packs(world)), "--runs", str(tmp_path / "runs"), "--work", str(tmp_path / "work")])
    got = capsys.readouterr()
    assert code == 0, got.err
    wrote = _wrote(got.err)
    loci = sorted((out / "loci").iterdir())
    want = {str(out / "hits.tsv"): ("rows", len((out / "hits.tsv").read_text().splitlines()) - 1),
            str(out / "hits.faa"): ("sequences", (out / "hits.faa").read_text().count(">")),
            str(out / "neighbors.tsv"): ("rows", _tsv_rows((out / "neighbors.tsv").read_text())),
            str(out / "loci"): ("tables", len(loci))}
    assert set(wrote) == set(want), got.err  # the final paths, never the temporary directory's
    for path, (unit, n) in want.items():
        assert wrote[path]["unit"].startswith(unit) and int(wrote[path]["n"]) == n, (path, wrote[path].group(0))
    assert int(wrote[str(out / "loci")]["rows"]) == sum(len(f.read_text().splitlines()) - 3 for f in loci)
    assert all(int(m["n"]) > 0 for m in wrote.values())


def test_w23_each_export_given_out_prints_its_count_and_none_without(rel, world, tmp_path, capsys):
    root = world[0]
    hit, _ = _seed(root)
    row = _run(world, _recipe(hit, hit, tool="jackhmmer", mx=2))
    assert row["status"] == "ok", row["error"]
    rid = row["run_id"].hex()
    where = ["--root", str(root), "--runs", str(_runs(world)), "--pack-dir", str(_packs(world))]
    uids = ",".join(sorted(u.hex() for u in _w06_anchors(root)))
    o = tmp_path
    commands = {
        o / "h.tsv": (["search", "hits", rid, *where, "--format", "tsv"], "rows"),
        o / "h.faa": (["search", "hits", rid, *where, "--format", "faa"], "sequences"),
        o / "n.tsv": (["nbr", "summarize", rid, *where, "--window", "5g", "--format", "tsv"], "rows"),
        o / "n.json": (["nbr", "summarize", uids, *where, "--release", REL, "--window", "5g", "--json"], "records"),
        o / "loci": (["nbr", rid, *where, "--window", "3g", "--format", "tsv"], "tables"),
        o / "p.faa": (["export", "proteins", "--root", str(rel), "--release", REL], "sequences"),
        o / "s.faa": (["seq", "--root", str(rel), "--release", REL, "--protein", sorted(_proteins(rel))[0].hex()],
                      "sequences"),
        o / "q.tsv": (_sql(rel, "-e", "SELECT key, value FROM catalog_info"), "rows"),
        o / "q.json": (_sql(rel, "-e", "SELECT key, value FROM catalog_info", "--format", "json"), "records"),
    }
    for path, (argv, unit) in commands.items():
        assert cli.main(argv) == 0, argv
        assert not _wrote(capsys.readouterr().err), argv  # to stdout: no count line
        assert cli.main([*argv, "--out", str(path)]) == 0, argv
        wrote = _wrote(capsys.readouterr().err)
        assert list(wrote) == [str(path)], (argv, wrote)
        m = wrote[str(path)]
        if unit == "tables":
            files = sorted(p for p in path.iterdir() if p.name != "RUN")  # card W-29: RUN is no table
            assert int(m["n"]) == len(files) > 0
            assert int(m["rows"]) == sum(len(f.read_text().splitlines()) - 3 for f in files)
            continue
        text = path.read_text()
        n = {"rows": _tsv_rows(text), "sequences": text.count("\n>") + text.startswith(">"),
             "records": len(text.splitlines())}[unit]
        assert m["unit"] == unit and int(m["n"]) == n > 0, (argv, m.group(0))
