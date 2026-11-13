"""Card W-24: silent fallbacks in the search readers, each a wrong or incomplete answer with exit 0.

(1) `search hits` and `nbr` on a run whose pack is not in the pack directory read the gene-call tape version
instead of the pack's: refused, naming the pack. (2) `nbr summarize RUN --release OTHER` read one release's run
against another release's files: refused unless the releases are equal. (3) On a release built before card B-21
the family postings hold each slot's canonical family only, so `fusions_of` finds nothing and family lookups miss
a protein's other domains: one note on stderr per connection, as `partners_of` gives for card B-39's column.
(4) `neighbors.tsv` repeated the cluster uid in `key` and `cluster_uid` at the cluster levels: the TSV drops the
copy, and the README says so.

The releases are the micro-S0 search world (`test_micro_s0_search.world`) and the B-21 postings fixture with a
two-domain protein (`test_micro_s0_postings.pfam`); a release as card B-21 left it is made here from the latter, its
family postings rewritten as a loop over the files says the B-10 builder wrote them (one per slot with a canonical
family, that family, segment 0).
"""

from __future__ import annotations

import json
import re
import shutil
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

sys.path.insert(0, str(Path(__file__).parent))
import test_micro_s0_postings as post  # noqa: E402
from test_micro_s0_postings import built, pfam  # noqa: E402,F401  (module-scoped postings fixtures)
from test_micro_s0_search import (LIMITS, REL, _copy, _dirs, _packs_of, _recipe, _runs_of, _seed,  # noqa: E402
                                  _w06_anchors, world)  # noqa: F401  (the module-scoped search world)

from biobox import catalog, cli, io  # noqa: E402
from biobox.registry.schemas import table  # noqa: E402
from biobox.search import export, runner  # noqa: E402

pytestmark = [pytest.mark.r2, pytest.mark.covers("W-24")]

OTHER = "rX.w24"  # a second release over the world's stage


def _search(root: Path, work: Path) -> dict:
    """One closed phmmer run on `root`'s release REL, its pack and runs beside the root (card W-03)."""
    hit, _ = _seed(root)
    row = runner.run(root, _recipe(hit, hit), REL, work_base=work, limits=LIMITS, packs=_packs_of(root),
                     runs=_runs_of(root))
    assert row["status"] == "ok", row["error"]
    return row


# --- (1) a run whose pack is not found ----------------------------------------------------------------------


@pytest.mark.covers("W-24", "S-09", "W-03")
def test_w24_hits_and_nbr_of_a_run_whose_pack_is_not_found_are_refused_naming_the_pack(world, tmp_path, capsys):
    """The tape version that names a run's genes is its pack's (card B-16), found by the run's pack id in the pack
    directory. Before W-24 a pack directory without that pack gave the gene-call producer's tape version, exit 0. Now
    each command that reads a run's genes exits 2, names the run's release, its pack id and the directory looked in,
    and writes nothing on stdout; with the pack directory that holds the pack, each answers."""
    root = world[0]
    row = _search(root, tmp_path / "work")
    rid, pid = row["run_id"].hex(), row["pack_id"].hex()
    empty = tmp_path / "no-packs"
    empty.mkdir()
    runs = ["--root", str(root), "--runs", str(_runs_of(root))]
    commands = {"hits-table": ["search", "hits", rid],
                "hits-json": ["search", "hits", rid, "--json"],
                "hits-tsv": ["search", "hits", rid, "--format", "tsv"],
                "hits-faa": ["search", "hits", rid, "--format", "faa"],
                "nbr": ["nbr", rid, "--window", "3g"],
                "nbr-tsv": ["nbr", rid, "--window", "3g", "--format", "tsv"],
                "nbr-anchors": ["nbr", "--anchors", rid, "--window", "3g", "--json"],
                "summarize": ["nbr", "summarize", rid, "--window", "5g", "--by", "nr50"]}
    for name, argv in commands.items():
        capsys.readouterr()
        code = cli.main([*argv, *runs, "--pack-dir", str(empty)])
        got = capsys.readouterr()
        assert code == 2, (name, got.err)
        assert got.out == "", name
        assert "refused" in got.err and pid in got.err and REL in got.err and str(empty) in got.err, (name, got.err)
    for name, argv in commands.items():
        assert cli.main([*argv, *runs, "--pack-dir", str(_packs_of(root))]) == 0, name
        assert capsys.readouterr().out.strip(), name


# --- (2) a run read against another release -----------------------------------------------------------------


@pytest.mark.covers("W-24", "W-12")
def test_w24_nbr_of_a_run_against_another_release_is_refused(world, tmp_path, capsys):
    """A run's anchors, tape version and hits are its own release's; `--release` naming another release made `nbr`
    read that release's tapes, clusters and contig_dim with them, exit 0. Refused now (exit 2, naming both releases
    and the uid-list form that does summarize given proteins on another release); `--release` equal to the run's is
    accepted, and protein uids with `--release OTHER` still answer."""
    root = _copy(world, tmp_path)
    catalog.define_release(root, OTHER, ["S0"])  # the same stage: every file of REL is OTHER's too
    row = _search(root, tmp_path / "work")
    rid = row["run_id"].hex()
    where = ["--root", str(root), *_dirs(root)]
    commands = {"summarize": ["nbr", "summarize", rid, "--window", "5g", "--by", "pfam,nr50", "--format", "tsv"],
                "nbr-tsv": ["nbr", rid, "--window", "3g", "--format", "tsv"],
                "nbr": ["nbr", rid, "--window", "3g"]}
    for name, argv in commands.items():
        capsys.readouterr()
        code = cli.main([*argv, *where, "--release", OTHER])
        got = capsys.readouterr()
        assert code == 2, (name, got.err)
        assert got.out == "" and "refused" in got.err and OTHER in got.err and REL in got.err, (name, got.err)
        assert "--anchors-file" in got.err, got.err
        assert cli.main([*argv, *where, "--release", REL]) == 0, name  # the run's own release, named
        named = capsys.readouterr().out
        assert cli.main([*argv, *where]) == 0, name
        assert capsys.readouterr().out == named, name
    uids = ",".join(sorted(h["target_uid"].hex() for h in runner.read_hits(runner.run_dir(_runs_of(root), rid))
                           if h["target_kind"] == "protein"))
    assert cli.main(["nbr", "summarize", uids, *where, "--release", OTHER, "--window", "5g", "--by", "nr50"]) == 0


# --- (3) family postings of a release built before card B-21 ------------------------------------------------


def _pre_b21(root: Path, dst: Path) -> Path:
    """A copy of `root` as card B-21 found a release: its tape_tokens without `domain_family_tokens` (card B-39's
    column, `post._without_domain_families`) and one family posting per slot with a canonical family: that family
    (`tape_tokens.family_tokens`), segment 0. The cluster postings are kept; the rows are re-sorted by the registry's
    sort key into one bucket, `all`, which the readers accept at every width."""
    old = post._without_domain_families(root, dst)
    canon = {(r["contig_uid"], i): t for r in post._rows(old, f"l3/tape_tokens/tape_version={post.TV}/release={REL}/*.parquet")
             for i, t in enumerate(r["family_tokens"])}
    rows = []
    for r in post._written(old):
        if post._level(r["token_id"]) != "family":
            rows.append(r)
        elif r["seg_idx"] == 0:
            rows.append({**r, "token_id": canon[(r["contig_uid"], r["ordinal"])]})
    assert all(t is not None for t in (r["token_id"] for r in rows))
    rows.sort(key=lambda r: tuple(r[c] for c in post.SORT))
    base = old / f"{post.POST}/release={REL}/tape_version={post.TV}"
    shutil.rmtree(base)
    spec = table("token_postings")
    io.write_table(pa.Table.from_pylist(rows, schema=spec.schema), spec, base / "token_bucket=all/part-0.parquet")
    return old


def _catalog(root: Path, path: Path):
    return catalog.open_catalog(catalog.generate_catalog(catalog.build_manifest(root, REL, ["S0"]), path))


def _notes(err: str) -> list[str]:
    return [x for x in err.splitlines() if "B-21" in x]


@pytest.mark.covers("W-24", "Q3", "Q4", "LY-09")
def test_w24_family_lookups_on_a_release_built_before_b21_say_so_once(pfam, built, tmp_path, capfd):
    """Card B-21 put one family posting per included domain on a slot; before it a slot had its canonical family
    only, and `fusions_of` on such a release returns no row for two families that are fused, exit 0. The release
    holds domains, and the domain table has a protein with two, but no family posting has a segment above 0: every
    family lookup (fusions_of, token_rows, q1_count, partners_of's anchors) says so on stderr, once per connection.
    Cluster lookups say nothing, nor does a release whose postings are complete: built today (card B-39's column), built
    between B-21 and B-39 (no column, segments above 0), or one without a protein of two domains."""
    root, p, _ = pfam
    tid = {k: r["token_id"] for k, r in post._family_dict(root).items()}
    fused, partner = tid["PF00664"], tid["PF00005"]
    toks = post._rows(root, f"l3/tape_tokens/tape_version={post.TV}/release={REL}/*.parquet")
    nr90 = next(t for r in toks for t in r["nr90_tokens"] if t is not None)
    nr50 = next(t for r in toks for t in r["nr50_tokens"] if t is not None)
    old = _pre_b21(root, tmp_path / "old")
    between = post._without_domain_families(root, tmp_path / "between")
    plain = post._without_domain_families(built[0], tmp_path / "plain")  # no protein has two domains
    family_queries = [("SELECT * FROM fusions_of(?, ?)", [fused, partner]),
                      ("SELECT * FROM fusions_of(?, ?)", [partner, fused]),
                      ("SELECT count(*) FROM token_rows(?)", [partner]),
                      ("SELECT q1_count(?, ?, 3, NULL)", [partner, nr90]),
                      ("SELECT * FROM partners_of(?, 3)", [partner])]
    cluster_queries = [("SELECT count(*) FROM token_rows(?)", [nr90]),
                       ("SELECT q1_count(?, ?, 3, NULL)", [nr50, nr90]),
                       ("SELECT * FROM pairs_of(?, ?, 3)", [nr50, nr90])]

    def run(r: Path, name: str, queries) -> tuple[list, str]:
        con = _catalog(r, tmp_path / f"{name}.duckdb")
        try:
            capfd.readouterr()
            got = [con.execute(sql, args).fetchall() for sql, args in queries]
            return got, capfd.readouterr().err
        finally:
            con.close()

    new, err = run(root, "new", family_queries)
    assert new[0], "the two-domain protein is fused today"
    assert not _notes(err), err
    got, err = run(old, "old", family_queries)
    assert got[0] == [] and got[1] == []  # the incomplete answer the note is for
    notes = _notes(err)
    assert len(notes) == 1, err
    assert REL in notes[0] and "fusions_of" in notes[0] and "canonical" in notes[0], notes
    _, err = run(old, "old-clusters", cluster_queries)
    assert not _notes(err), err
    got, err = run(between, "between", family_queries)
    assert got[0] == new[0] and not _notes(err), err
    plain_family = sorted(r["token_id"] for r in post._family_dict(built[0]).values())
    assert plain_family
    got, err = run(plain, "plain", [("SELECT count(*) FROM token_rows(?)", [t]) for t in plain_family])
    assert all(n for (n,), in got) and not _notes(err), err


# --- (4) neighbors.tsv at the cluster levels ----------------------------------------------------------------


def _tables(text: str) -> dict[str, tuple[list[str], list[list[str]]]]:
    """`nbr summarize --format tsv` read back: level -> (column line, rows), header lines skipped."""
    out, cols = {}, None
    for line in text.splitlines():
        if not line:
            cols = None
        elif line.startswith("#"):
            continue
        elif cols is None:
            cols = line.split("\t")
        else:
            row = line.split("\t")
            out.setdefault(row[cols.index("level")], (cols, []))[1].append(row)
    return out


@pytest.mark.covers("W-24", "W-06", "S-13")
def test_w24_neighbors_tsv_names_a_cluster_once_and_the_readme_says_where(world, capsys):
    """At the NR50 and NR90 levels `key` is the cluster uid; the TSV no longer repeats it as `cluster_uid` (no column
    twice, no other column equal to `key` in every row). The JSON records keep `cluster_uid` (`nbr_summary/1`, locked), equal to
    `key`; the README says the cluster level's key is the cluster uid and that model_name is empty without a domain
    library."""
    root = world[0]
    uids = ",".join(sorted(u.hex() for u in _w06_anchors(root)))
    args = ["nbr", "summarize", uids, "--root", str(root), *_dirs(root), "--window", "5g", "--release", REL,
            "--by", "pfam,nr50,nr90"]
    assert cli.main([*args, "--format", "tsv"]) == 0
    tables = _tables(capsys.readouterr().out)
    assert {"nr50", "nr90"} <= set(tables)
    for level in ("nr50", "nr90"):
        cols, rows = tables[level]
        assert "cluster_uid" not in cols and len(cols) == len(set(cols)), cols
        keyed = [r for r in rows if r[cols.index("unannotated")] == "false"]
        assert keyed and all(re.fullmatch(r"[0-9a-f]{32}", r[cols.index("key")]) for r in keyed)
        k = cols.index("key")
        for j, c in enumerate(cols):  # no other column repeats the cluster uid in every row
            assert j == k or any(r[j] != r[k] for r in keyed), (level, c)
    assert cli.main([*args, "--json"]) == 0
    recs = [json.loads(x) for x in capsys.readouterr().out.splitlines()[1:]]
    nr = [r for r in recs if r["level"] != "pfam" and not r["unannotated"]]
    assert nr and all(r["cluster_uid"] == r["key"] for r in nr)
    readme = " ".join(export.README_TEXT.split())
    assert "at the cluster levels key is the cluster uid" in readme, readme
    assert "empty on a release without a domain library" in readme, readme
