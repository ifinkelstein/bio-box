"""Ring R2: counts by genome fan out; counts by contig, locus and neighbour do not (card C-21c, DEC-27).

A metamorphic relation on the micro-S0 search world (test_micro_s0_search): two copies of its root that differ only in
`genome_input`. In `base` each of the three contigs is a genome of its own (lambda, T4, M. genitalium). In `shared` one
more genome repeats two existing contigs: GX = (lambda, T4), a virus genome of code 11 like theirs, so lambda and T4
each lie in two genomes and M. genitalium in one. Tapes, postings, tokens, clusters, the pack and the search run are
the same files in both: a contig has one tape and one call set whatever the genomes that hold it (DEC-26, DEC-27).

Each root's `contig_dim` and `genome_dim` are written by hand from that construction (DEC-27's rules: a contig's
genome is the least of its genomes, `n_genomes` their number, a genome attribute null where they disagree) with
tiers by genome, LAM and T4 1, MG 2, GX 3, so tier is constructed truth; `dim.run` is checked to write the same rows
but the tiers (no tiering step exists yet). Then, against oracles that read the fixture's genome lists and the
Parquet shards (the pair oracle of tests/oracles, the S-09 window oracle of tests/r1), never the code under test:

    unchanged       postings, q1_count, partners_of (n_loci, n_pairs), nbr summarize's n_loci, n_neighbors,
                    offsets and header counts, the hits.tsv rows, the JSON records and the locus tables but for
                    what follows
    by genome       q1_genomes (new), every n_genomes of nbr summarize, `bio-box releases`' genomes: the distinct
                    genomes of the contigs counted, so a count rises by one where GX holds a contig counted
    semi-join       q1_count's tier filter: a locus on lambda counts in tier 1 (LAM) and tier 3 (GX), once in each
    listed          a hit or locus lists every genome that carries its contig (hits.tsv `n_genomes`, `genome_uids`;
                    the JSON occurrence; the locus table's facts line), on every release: one schema everywhere

On a release without a shared contig (`base`) those appended columns, keys and facts are 1 and the contig's one
genome, and every other output is what it was before this card (checked against the code before it when the card
was made; a test cannot hold that code).
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog, cli, dim, fasta, ids
from biobox import io as bio
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.l3 import tapes
from biobox.registry.schemas import table
from biobox.search import export, nbr

sys.path.insert(0, str(Path(__file__).parent))
import test_micro_s0_search as s  # noqa: E402
from test_micro_s0_search import world  # noqa: E402,F401  (the micro-S0 search world, built once for this module)

pytestmark = [pytest.mark.r2, pytest.mark.covers("C-01")]

REL = s.REL
micro = s.micro
TV = tapes.tape_version([s.V])
K = 10
WINDOW = "5g"
LEVELS = ("pfam", "nr50", "nr90")
GI = "/l0seq/genome_input/"


def _cu(acc: str) -> bytes:
    return ids.contig_uid(fasta.read(micro.FIX / f"{acc}.fna.gz")[0].seq)


LAM, T4, MG = _cu("NC_001416.1"), _cu("NC_000866.4"), _cu("NC_000908.2")
OBJ = {_cu(o["accession"]): o for o in micro.CONTRACT["objects"]}
BASE = [GenomeInput((cu,), o["genetic_code"], kind=o["kind"], circular=(cu,) if o["topology"] == "circular" else ())
        for cu, o in sorted(OBJ.items())]
GX = GenomeInput(tuple(sorted((LAM, T4))), 11, kind="virus")  # repeats two existing contigs
GUID = {g.contig_uids: ids.genome_uid(list(g.contig_uids)) for g in [*BASE, GX]}
TIER = {GUID[(LAM,)]: 1, GUID[(T4,)]: 1, GUID[(MG,)]: 2, GUID[GX.contig_uids]: 3}
WORLDS = {"base": BASE, "shared": [*BASE, GX]}


def genomes_of(name: str) -> dict[bytes, set[bytes]]:
    """{contig: its genomes} in world `name`, from the fixture's genome lists."""
    out: dict[bytes, set[bytes]] = {}
    for g in WORLDS[name]:
        for cu in g.contig_uids:
            out.setdefault(cu, set()).add(GUID[g.contig_uids])
    return out


def _dim_rows(name: str, lens: dict[bytes, int]) -> tuple[list[dict], list[dict]]:
    """`contig_dim` and `genome_dim` rows of world `name` by DEC-27's rules, with the tiers of TIER."""
    of = genomes_of(name)
    cdim = []
    for cu in sorted(of):
        tiers = {TIER[g] for g in of[cu]}
        cdim.append({"contig_uid": cu, "genome_uid": min(of[cu]), "contig_len": lens[cu],
                     "is_circular": OBJ[cu]["topology"] == "circular", "tier": tiers.pop() if len(tiers) == 1 else None,
                     "assembly_quality": "complete", "tax_order": None, "is_species_rep": False,
                     "n_genomes": len(of[cu])})
    gdim = [{"genome_uid": GUID[g.contig_uids], "n_contigs": len(g.contig_uids),
             "genome_len": sum(lens[cu] for cu in g.contig_uids), "kind": g.kind, "genetic_code": g.genetic_code,
             "assembly_quality": "complete", "tier": TIER[GUID[g.contig_uids]], "tax_order": None,
             "is_species_rep": False} for g in WORLDS[name]]
    return cdim, sorted(gdim, key=lambda r: r["genome_uid"])


def _write(name: str, rows: list[dict], path: Path) -> None:
    spec = table(name)
    bio.write_table([pa.Table.from_pylist(rows, schema=spec.schema)], spec, path)


def _parts(root: Path, rel: str) -> list[dict]:
    return [r for p in bio.part_files(root / rel) for r in pq.read_table(p, partitioning=None).to_pylist()]


@pytest.fixture(scope="module")
def c21c(world, tmp_path_factory):  # noqa: F811
    """The two roots, each with the world's search run and pack beside it, its `genome_input`, hand-written dims and a
    catalog; and what `dim.run` writes for each, in a scratch output root."""
    hit, _ = s._seed(world[0])
    row = s._run(world, s._recipe(hit, hit))
    assert row["status"] == "ok", row["error"]
    rid = row["run_id"]
    d = tmp_path_factory.mktemp("c21c")
    lens = {cu: t["contig_len"] for cu, t in s._tape_rows(world[0], TV).items()}
    out = {"rid": rid.hex(), "d": d}
    for name, genomes in WORLDS.items():
        root = d / name / "a"
        shutil.copytree(world[0], root)
        shutil.copytree(s._rd(world, rid), s._rdir(root, rid))
        shutil.copytree(s._packs(world), s._packs_of(root))
        write_genome_input(root, "S0", genomes)
        made = d / name / "dim_run"
        dim.run(root, REL, ["S0"], tmp=d, out=made)
        cdim, gdim = _dim_rows(name, lens)
        _write("contig_dim", cdim, root / f"dim/contig_dim/release={REL}/part-0.parquet")
        _write("genome_dim", gdim, root / f"dim/genome_dim/release={REL}/part-0.parquet")
        cat = catalog.generate_catalog(catalog.build_manifest(root, REL, ["S0"]), d / name / "c.duckdb")
        out[name] = {"root": root, "made": made, "cdim": cdim, "gdim": gdim, "catalog": cat}
    return out


def _con(c21c, name: str):
    return catalog.open_catalog(c21c[name]["catalog"])


def _tier(c21c, name: str, cu: bytes) -> int | None:
    """The contig's `contig_dim` tier in world `name`: null in the shared world where its genomes' tiers differ."""
    return next(r["tier"] for r in c21c[name]["cdim"] if r["contig_uid"] == cu)


def test_c21c_the_two_worlds_differ_in_membership_only(c21c):
    """The fixture: `dim.run` writes the hand-written dimension rows but the tiers; the shared world's membership holds
    lambda and T4 in two genomes; every file but `genome_input` and the dims is the same in both roots."""
    for name in WORLDS:
        w = c21c[name]
        for tbl, want in (("contig_dim", w["cdim"]), ("genome_dim", w["gdim"])):
            got = _parts(w["made"], f"dim/{tbl}/release={REL}")
            assert got == [{**r, "tier": None} for r in want], (name, tbl)
    con = _con(c21c, "shared")
    try:
        got = con.execute("SELECT contig_uid, count(DISTINCT genome_uid) FROM genome_member GROUP BY ALL").fetchall()
    finally:
        con.close()
    assert dict(got) == {LAM: 2, T4: 2, MG: 1}
    a, b = c21c["base"]["root"], c21c["shared"]["root"]
    for p in sorted(x for x in a.rglob("*.parquet") if not x.relative_to(a).as_posix().startswith(
            ("l0seq/genome_input/", "dim/"))):
        assert p.read_bytes() == (b / p.relative_to(a)).read_bytes(), p.relative_to(a)


# --- the catalog's neighbour macros (Q1, Q3) -------------------------------------------------------------------------


def _pairs(root: Path) -> dict[tuple[int, int], list[tuple]]:
    """Every (A token, B token) within K slots, as rows of the pair oracle (tests/oracles/pairs_oracle.py)."""
    return micro._oracle_pairs(micro._oracle_rows(root), K)


_PAIR_LIST = micro._PAIR_LIST


@pytest.mark.covers("Q1", "Q3")
def test_c21c_postings_q1_count_and_partners_of_do_not_change_when_a_genome_repeats_a_contig(c21c):
    """Counts by locus: the same in both worlds and equal to the pair oracle; partners_of's n_loci and n_pairs too."""
    a, b = c21c["base"]["root"], c21c["shared"]["root"]
    for p in sorted((a / "l3/token_postings").rglob("*.parquet")):
        assert p.read_bytes() == (b / p.relative_to(a)).read_bytes()
    full = _pairs(a)
    pairs = sorted(full)
    lists = {"a": [x for x, _ in pairs], "b": [y for _, y in pairs]}
    tokens = sorted({x for x, _ in pairs})
    got = {}
    for name in WORLDS:
        con = _con(c21c, name)
        try:
            got[name] = (sorted(con.execute(f"SELECT l.a, l.b, q1_count(l.a, l.b, {K}, NULL) FROM {_PAIR_LIST}",
                                            lists).fetchall()),
                         {t: con.execute("SELECT partner_token, n_loci, n_pairs FROM partners_of(?, ?)",
                                         [t, K]).fetchall() for t in tokens})
        finally:
            con.close()
    assert got["base"][0] == [(x, y, micro._loci(full[(x, y)])) for x, y in pairs]
    assert got["shared"] == got["base"]
    assert any(rows for rows in got["base"][1].values())


@pytest.mark.covers("Q1")
def test_c21c_q1_counts_tier_filter_is_a_semi_join_through_the_genomes_tiers(c21c):
    """A locus counts in tier t when one of its contig's genomes has tier t (`genome_dim`) or its `contig_dim` row
    has (a release whose tiers are set by contig only, test_micro_s0_l3): a locus on lambda counts once in tier 1
    and once in tier 3, never twice in one. `contig_dim`'s tier of lambda is null in the shared world (its genomes
    disagree), so a filter on `contig_dim` alone counts none of its loci."""
    full = _pairs(c21c["base"]["root"])
    pairs = sorted(full)
    lists = {"a": [x for x, _ in pairs], "b": [y for _, y in pairs]}
    sql = (f"SELECT l.a, l.b, q1_count(l.a, l.b, {K}, NULL), q1_count(l.a, l.b, {K}, 1), q1_count(l.a, l.b, {K}, 2), "
           f"q1_count(l.a, l.b, {K}, 3) FROM {_PAIR_LIST}")
    for name in WORLDS:
        of = genomes_of(name)
        tiers = {cu: {TIER[g] for g in gs} for cu, gs in of.items()}
        want = [(x, y, micro._loci(full[(x, y)]),
                 *(micro._loci([r for r in full[(x, y)] if t in tiers[r[0]]]) for t in (1, 2, 3))) for x, y in pairs]
        con = _con(c21c, name)
        try:
            assert sorted(con.execute(sql, lists).fetchall()) == want, name
            x, y = next(p for p in pairs if any(r[0] == LAM for r in full[p]))  # one statement, the tier a parameter
            got = con.execute(f"SELECT q1_count(?, ?, {K}, ?), q1_count(?, ?, {K}, ?)", [x, y, 1, x, y, 3]).fetchone()
        finally:
            con.close()
        assert list(got) == [want[pairs.index((x, y))][3], want[pairs.index((x, y))][5]], name
    # in the shared world a locus counts in two tiers (lambda's and T4's: tier 1 and tier 3), each once
    assert any(w[5] for w in want) and any(w[4] for w in want) and any(w[3] + w[4] + w[5] > w[2] for w in want)


@pytest.mark.covers("Q1")
def test_c21c_q1_genomes_counts_the_distinct_genomes_of_the_loci(c21c, tmp_path):
    """q1_genomes(A, B, k): the genomes that hold a contig with a locus of A with B within k, each once."""
    full = _pairs(c21c["base"]["root"])
    pairs = sorted(full)
    lists = {"a": [x for x, _ in pairs], "b": [y for _, y in pairs]}
    sql = f"SELECT l.a, l.b, q1_genomes(l.a, l.b, {K}) FROM {_PAIR_LIST}"
    seen = {}
    for name in WORLDS:
        of = genomes_of(name)
        want = sorted((x, y, len(set().union(*(of[r[0]] for r in full[(x, y)])))) for x, y in pairs)
        con = _con(c21c, name)
        try:
            assert sorted(con.execute(sql, lists).fetchall()) == want, name
            x, y = pairs[0]
            assert con.execute(f"SELECT q1_genomes(?, ?, {K})", [x, y]).fetchone()[0] == want[0][2]
            assert con.execute(f"SELECT q1_genomes(?, 0, {K})", [x]).fetchone()[0] == 0  # token 0: no locus
        finally:
            con.close()
        seen[name] = want
    rises = [(b_[2], s_[2]) for b_, s_ in zip(seen["base"], seen["shared"]) if b_ != s_]
    assert rises and all(sv == bv + 1 for bv, sv in rises)  # GX is one more genome wherever lambda or T4 is counted
    # a catalog built before this card holds no q1_genomes: open_catalog gives it one
    old = tmp_path / "old.duckdb"
    shutil.copyfile(c21c["shared"]["catalog"], old)
    raw = duckdb.connect(str(old))
    try:
        raw.execute("DROP MACRO q1_genomes")
    finally:
        raw.close()
    con = catalog.open_catalog(old)
    try:
        assert sorted(con.execute(sql, lists).fetchall()) == seen["shared"]
    finally:
        con.close()


ROW_CAP = 16  # rows per row group of the padded genome_input (test_genome_member_scale's method)


@pytest.mark.covers("LY-06")
def test_c21c_q1_genomes_reads_only_the_genome_input_row_groups_of_its_loci_contigs(c21c, tmp_path, record_property):
    """AGENTS §4.4: the loci's contigs reach the `genome_input` scan as a join filter. On a copy of the shared world
    whose `genome_input` also lists 3,000 made-up genomes in row groups of ROW_CAP rows, q1_genomes reads at most
    two row groups per loci contig (its rows, or two where they straddle a cut), not the file."""
    root = tmp_path / "a"
    shutil.copytree(c21c["shared"]["root"], root)
    shutil.rmtree(root / "l0seq/genome_input/stage=S0")
    rng = np.random.default_rng(2121)
    pad = [GenomeInput((rng.bytes(16),), 11, kind="virus") for _ in range(3_000)]
    real = bio.write_table

    def capped(batches, spec, path, *, row_group_rows=bio.ROW_GROUP_ROWS, **kw):
        return real(batches, spec, path, row_group_rows=min(row_group_rows, ROW_CAP), **kw)

    mp = pytest.MonkeyPatch()
    mp.setattr(bio, "write_table", capped)
    try:
        write_genome_input(root, "S0", [*WORLDS["shared"], *pad])
    finally:
        mp.undo()
    total = sum(pq.read_metadata(p).num_row_groups for p in bio.part_files(root / "l0seq/genome_input/stage=S0"))
    assert total >= 100
    full = _pairs(root)
    of = genomes_of("shared")
    con = catalog.open_catalog(catalog.generate_catalog(catalog.build_manifest(root, REL, ["S0"]), tmp_path / "c.duckdb"))
    try:
        checked = 0
        for x, y in sorted(full)[:: max(1, len(full) // 20)]:
            contigs = {r[0] for r in full[(x, y)]}
            prof = catalog.profile_query(con, f"SELECT q1_genomes(?, ?, {K})", [x, y])
            assert prof["result"] == [(len(set().union(*(of[c] for c in contigs))),)]
            scans = [sc for sc in prof["scans"] if any(GI in f for f in sc["files"])]
            assert scans, prof["scans"]
            read = sum(sc["row_groups_scanned"] for sc in scans)
            record_property(f"rg_{x}_{y}", f"{read} of {total} for {len(contigs)} contigs")
            assert read <= 2 * len(contigs) < total, (x, y, read, total)
            checked += 1
    finally:
        con.close()
    assert checked >= 5


# --- nbr summarize (W-06) --------------------------------------------------------------------------------------------


def _key_contigs(root: Path, anchors: set[bytes], window: str, levels) -> tuple[dict, set]:
    """{(level, key): the contigs of its neighbour rows} ('' the unannotated key) and the contigs of the anchor loci, by
    the S-09 window oracle over every tape slot and the annotations read from the shards (as `_oracle_summary`)."""
    r1 = s._r1()
    policy = nbr.Policy.parse(window=window)
    st = s._w06_store(root)
    out: dict[tuple[str, str], set] = {}
    loci: set = set()
    for cu, t in sorted(s._tape_rows(root, TV).items()):
        slots = list(zip(t["start_bps"], t["end_bps"], t["strands"]))
        for i, pu in enumerate(t["protein_uids"]):
            if pu not in anchors:
                continue
            got = r1._oracle_nbr(slots, t["contig_len"], t["is_circular"], i, policy)
            if got is None:
                continue
            loci.add(cu)
            for j, dd, *_ in got:
                if dd == 0:
                    continue
                u = t["protein_uids"][j]
                for level in levels:
                    if level == "pfam":
                        keys = {f"{lib}:{acc}" for lib, acc in st["dom"].get(u, {})}
                    else:
                        keys = {st["clus"][level][u][0].hex()} if u in st["clus"][level] else set()
                    for k in keys or {""}:
                        out.setdefault((level, k), set()).add(cu)
    return out, loci


@pytest.mark.covers("W-06")
def test_c21c_nbr_summarize_counts_every_genome_of_a_shared_contig_and_no_neighbour_twice(c21c):
    """Per key, `n_genomes` is the distinct genomes of the contigs its neighbours lie on, and the header's
    `n_anchor_genomes` those of the anchor loci; every other number is the base world's, which is the W-06 oracle's.
    (`genome_of` pointed at the membership without care would count each neighbour once per genome of its contig.)"""
    base_root = c21c["base"]["root"]
    anchors = s._w06_anchors(base_root)
    policy = nbr.Policy.parse(window=WINDOW)
    keys, loci = _key_contigs(base_root, anchors, WINDOW, LEVELS)
    got = {name: nbr.summarize(c21c[name]["root"], anchors, policy, by=LEVELS, release=REL, top_n=0)
           for name in WORLDS}
    s._compare(got["base"], s._oracle_summary(base_root, anchors, WINDOW, levels=LEVELS))
    rises = 0
    for name in WORLDS:
        of = genomes_of(name)
        h = got[name]["header"]
        assert h["genome_basis"] == "genome_uid"
        assert h["n_anchor_genomes"] == len(set().union(*(of[c] for c in loci))), name
        for r in got[name]["rows"]:
            cs = keys.get((r["level"], r["key"]), set())
            assert r["n_genomes"] == len(set().union(set(), *(of[c] for c in cs))), (name, r["level"], r["key"])
    strip_h = ("n_anchor_genomes",)
    assert {k: v for k, v in got["shared"]["header"].items() if k not in strip_h} == \
        {k: v for k, v in got["base"]["header"].items() if k not in strip_h}
    assert len(got["shared"]["rows"]) == len(got["base"]["rows"])
    for a, b in zip(got["base"]["rows"], got["shared"]["rows"]):
        assert {k: v for k, v in a.items() if k != "n_genomes"} == {k: v for k, v in b.items() if k != "n_genomes"}
        rises += b["n_genomes"] == a["n_genomes"] + 1
        assert b["n_genomes"] in (a["n_genomes"], a["n_genomes"] + 1)
    assert rises and got["shared"]["header"]["n_anchor_genomes"] == got["base"]["header"]["n_anchor_genomes"] + 1


# --- exports (W-05): hits.tsv, the JSON hits, the locus tables ------------------------------------------------------


def _hits_tsv(c21c, name: str) -> list[list[str]]:
    root = c21c[name]["root"]
    lines = export.hits_tsv(root, s._w05_hits(root), release=REL, tape_version=TV, gene_call_version=s.V)
    return [x.split("\t") for x in lines]


def _gl(of: dict[bytes, set[bytes]], cu: bytes) -> list[str]:
    """The contig's genomes as hex, in uid order."""
    return sorted(g.hex() for g in of[cu])


@pytest.mark.covers("W-05")
def test_c21c_hits_tsv_lists_every_genome_of_a_shared_contig_and_keeps_one_row_per_gene(c21c):
    """Both tables end with `n_genomes` and `genome_uids` (comma-separated, in uid order), on every release: in the
    base world 1 and the contig's one genome. The shared world's table holds the base world's rows in the base world's
    order; `genome_uid` and `tier` are the contig's `contig_dim` row (its least genome; the tier its genomes agree on,
    else empty), and on lambda and T4 `n_genomes` is 2 and `genome_uids` lists LAM or T4 and GX."""
    base, shared = _hits_tsv(c21c, "base"), _hits_tsv(c21c, "shared")
    assert base[0] == shared[0] == [*export.HITS_TSV_COLUMNS]  # tier included: the hand-written dims hold tiers
    assert base[0][-2:] == ["n_genomes", "genome_uids"]
    assert len(shared) == len(base) > 1
    cols = base[0]
    seen = set()
    for name, rows in (("base", base), ("shared", shared)):
        of = genomes_of(name)
        for line in rows[1:]:
            r = dict(zip(cols, line))
            cu = bytes.fromhex(r["contig_uid"])
            seen.add(cu)
            t = _tier(c21c, name, cu)
            assert (r["genome_uid"], r["tier"]) == (min(of[cu]).hex(), "" if t is None else str(t)), name
            assert (r["n_genomes"], r["genome_uids"]) == (str(len(of[cu])), ",".join(_gl(of, cu))), name
    keep = [k for k, c in enumerate(cols) if c not in ("genome_uid", "tier", "n_genomes", "genome_uids")]
    assert [[x[k] for k in keep] for x in base] == [[x[k] for k in keep] for x in shared]
    assert seen == {LAM, T4, MG}
    assert {x[cols.index("n_genomes")] for x in shared[1:]} == {"1", "2"}


def _json_hits(c21c, name: str, capsys) -> list[dict]:
    root = c21c[name]["root"]
    capsys.readouterr()
    assert cli.main(["search", "hits", c21c["rid"], "--root", str(root), *s._dirs(root), "--json", "--members",
                     "--top_n", "0"]) == 0
    return [json.loads(x) for x in capsys.readouterr().out.splitlines()]


@pytest.mark.covers("W-05", "S-13")
def test_c21c_the_json_occurrence_lists_every_genome_of_its_contig(c21c, capsys):
    """`search_hit/2`: every occurrence ends with `n_genomes` and `genome_uids` (keys added at the same version, S-13;
    in the base world 1 and the one genome); `genome_uid` and `tier` are its `contig_dim` row's (the least genome; the
    tier its genomes agree on, else null). Every other key and record is the same in both worlds."""
    got = {name: _json_hits(c21c, name, capsys) for name in WORLDS}
    assert len(got["base"]) == len(got["shared"]) > 0
    keys = ["contig_uid", "ordinal", "feature_uid", "start", "end", "strand", "genome_uid", "tier", "n_genomes",
            "genome_uids"]
    n_shared = 0
    for name, recs in got.items():
        of = genomes_of(name)
        for o in (o for h in recs for o in h["occurrences"]):
            cu = bytes.fromhex(o["contig_uid"])
            assert list(o) == keys
            assert (o["genome_uid"], o["tier"], o["n_genomes"], o["genome_uids"]) == (
                min(of[cu]).hex(), _tier(c21c, name, cu), len(of[cu]), _gl(of, cu)), name
            n_shared += name == "shared" and len(of[cu]) > 1
    assert n_shared

    def strip(recs):
        return [{**h, "occurrences": [{k: v for k, v in o.items() if k not in ("genome_uid", "tier", "n_genomes",
                                                                              "genome_uids")} for o in h["occurrences"]]}
                for h in recs]

    assert strip(got["base"]) == strip(got["shared"])
    locked = json.loads(export.lock_path(s.REPO).read_text())
    assert export.lock_problems(locked) == []


def _loci_files(c21c, name: str, tmp_path: Path, capsys) -> dict[str, list[str]]:
    root = c21c[name]["root"]
    out = tmp_path / name
    uids = ",".join(sorted(u.hex() for u in s._w06_anchors(c21c["base"]["root"])))
    assert cli.main(["nbr", uids, "--release", REL, "--root", str(root), *s._dirs(root), "--window", WINDOW,
                     "--format", "tsv", "--out", str(out)]) == 0
    capsys.readouterr()
    return {p.name: p.read_text().splitlines() for p in sorted(out.iterdir())}


@pytest.mark.covers("W-05")
def test_c21c_a_locus_table_lists_every_genome_of_its_contig_in_its_facts(c21c, tmp_path, capsys):
    """The facts line of every locus ends with `n_genomes=` and `genome_uids=` (comma-separated, in uid order): 1 and
    the one genome in the base world, 2 and both genomes on lambda and T4 in the shared world; every other line is the
    same in both."""
    got = {name: _loci_files(c21c, name, tmp_path, capsys) for name in WORLDS}
    assert sorted(got["base"]) == sorted(got["shared"]) and got["base"]
    on = set()
    for name, files in got.items():
        of = genomes_of(name)
        for f, lines in files.items():
            facts = [x.split("=", 1) for x in lines[1].removeprefix("# ").split("\t")]
            cu = bytes.fromhex(dict(facts)["contig_uid"])
            on.add((name, len(of[cu])))
            assert facts[-2:] == [["n_genomes", str(len(of[cu]))], ["genome_uids", ",".join(_gl(of, cu))]], (name, f)
    for f, lines in got["base"].items():
        other = got["shared"][f]
        assert lines[0] == other[0] and lines[2:] == other[2:], f
        assert lines[1].rsplit("\tn_genomes=", 1)[0] == other[1].rsplit("\tn_genomes=", 1)[0], f
    assert on == {("base", 1), ("shared", 1), ("shared", 2)}


# --- `bio-box releases` (W-23) ---------------------------------------------------------------------------------------


@pytest.mark.covers("W-23")
def test_c21c_releases_counts_every_genome_of_the_release(c21c, tmp_path):
    """A genome whose every contig another genome also holds (GX) is a genome of the release: 3 in base, 4 in shared;
    from `genome_dim`, and from the membership (`genome_input`) for a release that has no `genome_dim`."""
    assert cli.release_summary(c21c["base"]["root"], REL)["genomes"] == 3
    assert cli.release_summary(c21c["shared"]["root"], REL)["genomes"] == 4
    root = tmp_path / "a"
    shutil.copytree(c21c["shared"]["root"], root)
    shutil.rmtree(root / f"dim/genome_dim/release={REL}")
    assert cli.release_summary(root, REL)["genomes"] == 4
